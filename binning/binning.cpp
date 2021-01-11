#include <omp.h>
#include <algorithm>
#include <tuple>
#include <cstring>

#include "CSC.hpp"
#include "CSR.hpp"
#include "util.hpp"
#include "radix_sort.hpp"

#define MATRIX_TYPE double
#define SIZE 16
static uint32_t nrows_per_blocker;

template <class T>
void symbolic(const CSC<T> &A, const CSR<T> &B, int startIdx, int endIdx,
              uint32_t nrows_per_blocker, uint16_t num_blockers,
              int *flops_by_row_blockers, int &total_flops)
{
    #pragma omp parallel for reduction(+ : flops_by_row_blockers[:num_blockers])
    for (int i = startIdx; i < endIdx; ++i)
    {
        int rownnz = B.rowptr[i + 1] - B.rowptr[i];
        for (int j = A.colptr[i]; j < A.colptr[i + 1]; ++j)
        {
            uint16_t row_blocker_id = A.rowids[j] / nrows_per_blocker;
            flops_by_row_blockers[row_blocker_id] += rownnz;
        }
    }
    for (int i = 0; i < num_blockers; ++i)
    {
        total_flops += flops_by_row_blockers[i];
    }
}

struct ExtractKey
{
    inline int64_t operator()(tuple<int, int, MATRIX_TYPE> tup)
    {
        int64_t res = std::get<0>(tup);
        res = (res << 32);
        res = res | (int64_t)(uint32_t) std::get<1>(tup);
        return res;
    }
};

// Key for Radix Sort
struct ExtractKey2
{
    inline uint32_t operator()(tuple<int, int, MATRIX_TYPE> tup)
    {
		return ((std::get<0>(tup) % nrows_per_blocker) << 20 | (uint32_t) std::get<1>(tup));
    }
};


// Prefix sum (Sequential)
template <typename T> void seq_scan(T *in, T *out, T N) {
  out[0] = 0;
  for (T i = 0; i < N - 1; ++i) {
    out[i + 1] = out[i] + in[i];
  }
}


// Prefix sum (Thread parallel)
template <typename T> void scan(T *in, T *out, T N) {
  // if the array is comparatively small, use sequential scan instead
  if (N < (1 << 17)) {
    seq_scan(in, out, N);
  } else {
    int tnum = 1;
    #pragma omp parallel
    { tnum = omp_get_num_threads(); }
    T each_n = N / tnum;
    T *partial_sum = _malloc<T>(tnum);
    #pragma omp parallel
    {
      // thead level prefix summing
      int tid = omp_get_thread_num();
      T start = each_n * tid;
      T end = (tid < tnum - 1) ? start + each_n : N;
      out[start] = 0;
      for (T i = start; i < end - 1; ++i) {
        out[i + 1] = out[i] + in[i];
      }
      // calculate offset in every thread
      partial_sum[tid] = out[end - 1] + in[end - 1];
    #pragma omp barrier

      T offset = 0;
      for (int i = 0; i < tid; ++i) {
        offset += partial_sum[i];
      }
      for (T i = start; i < end; ++i) {
        out[i] += offset;
      }
    }
    _free<T>(partial_sum);
  }
}


template <typename T>
inline bool isTupleEqual (tuple<int, int, T> t1, tuple<int, int, T> t2)
{
    if (std::get<1>(t1) != std::get<1>(t2))
        return false;
    if (std::get<0>(t1) != std::get<0>(t2))
        return false;
    return true;
}


template <typename T>
inline int doMerge(tuple<int, int, T>* vec, int length)
{
    if (length == 0) return 0;
    ExtractKey op = ExtractKey();
    int i = 0;
    int j = 1;

    while (i < length && j < length)
    {
        if (j < length && isTupleEqual (vec[i], vec[j]))
            std::get<2>(vec[i]) += std::get<2>(vec[j]);
        else
        {
            ++i;
            std::get<0>(vec[i]) = std::get<0>(vec[j]);
            std::get<1>(vec[i]) = std::get<1>(vec[j]);
            std::get<2>(vec[i]) = std::get<2>(vec[j]);
        }
        ++j;
    }
    return i + 1;
}


template <typename T>
inline void doRadixSort(tuple<int, int, NT>* begin, tuple<int, int, T>* end, tuple<int, int T>* buffer)
{
    radix_sort(begin, end, buffer, ExtractKey2());
}


template <class T>
void BinSpGEMM(const CSC<T> &A, const CSR<T> &B, CSR<T> &C)
{
    typedef tuple<int, int, T> TripleNode;
    const uint16_t nthreads = omp_get_max_threads();
    uint16_t num_blockers = 256;
    const uint16_t block_width = 128;
    int nrows_of_A = A.rows;
    nrows_per_blocker = A.rows <= num_blockers * 64 ? 64 : (A.rows + num_blockers - 1) / num_blockers;

    int total_flop = 0;

    int *flops_by_row_blockers = _malloc<int>(num_blockers);
    int *flops_by_rows = _malloc<int>(A.rows);
    int *nnz_by_row = _malloc<int>(A.rows);

    int *global_blocker_counters = _malloc<int>(num_blockers);
    TripleNode **global_blockers = _malloc<TripleNode *>(num_blockers);
    int **local_blocker_counters = _malloc<int *>(nthreads);
    TripleNode **local_blockers = _malloc<TripleNode *>(nthreads);
    TripleNode **sorting_buffer = _malloc<TripleNode *>(nthreads);
    int *nnz_per_row_blocker = _malloc<int>(num_blockers);

    // Symbolic Phase
    symbolic(A, B, 0, A.cols, nrows_per_blocker, num_blockers, flops_by_row_blockers, total_flop);
    for (uint16_t blocker_id = 0; blocker_id < num_blockers; ++blocker_id){
        global_blockers[blocker_id] = static_cast<TripleNode *>(::operator new(SIZE * flops_by_row_blockers[blocker_id]));
    }
    int max_flops_in_row_blockers = *std::max_element(flops_by_row_blockers, flops_by_row_blockers + num_blockers);

    #pragma omp parallel
    {
        uint16_t thread_id = omp_get_thread_num();
        TripleNode *begin_local_blockers, *cur_local_blockers, *end_local_blockers, *cur_global_blockers;
        local_blockers[thread_id] = static_cast<TripleNode *>(::operator new(SIZE *num_blockers *block_width));
        local_blocker_counters[thread_id] = _malloc<IT>(num_blockers);
        sorting_buffer[thread_id] = static_cast<TripleNode *>(::operator new(SIZE *max_flops_in_row_blockers));

        // computing phase
        #pragma omp for nowait
        for (int idx = 0; idx < A.cols; ++idx)
        {
            for (int j = A.colptr[idx]; j < A.colptr[idx + 1]; ++j) // ncols(A) * 4
            {
                int rowid = A.rowids[j]; // nnz(A) * 4
                uint16_t row_blocker_id = rowid / nrows_per_blocker;
                begin_local_blockers = local_blockers[thread_id] + row_blocker_id * block_width;
                cur_local_blockers = begin_local_blockers + local_blocker_counters[thread_id][row_blocker_id];
                end_local_blockers = begin_local_blockers + block_width;
                for (int k = B.rowptr[idx]; k < B.rowptr[idx + 1]; ++k) // nrows(B) * 4
                {

                    std::get<0>(*cur_local_blockers) = A.rowids[j];
                    std::get<1>(*cur_local_blockers) = B.colids[k];
                    std::get<2>(*cur_local_blockers) = A.values[j] * B.values[k];
                    cur_local_blockers++;
                    if (cur_local_blockers == end_local_blockers) // flop * 16
                    {

                        std::memcpy(
                            global_blockers[row_blocker_id] + __sync_fetch_and_add(&global_blocker_counters[row_blocker_id], block_width),
                            begin_local_blockers,
                            block_width * SIZE);
                        cur_local_blockers = begin_local_blockers;
                    }
                }
                local_blocker_counters[thread_id][row_blocker_id] = cur_local_blockers - begin_local_blockers;
            }
        }
        for (uint16_t row_blocker_id = 0; row_blocker_id < num_blockers; row_blocker_id++)
        {

            std::memcpy(
                global_blockers[row_blocker_id] + __sync_fetch_and_add(&global_blocker_counters[row_blocker_id], local_blocker_counters[thread_id][row_blocker_id]),
                local_blockers[thread_id] + row_blocker_id * block_width,
                local_blocker_counters[thread_id][row_blocker_id] * SIZE);
        }
    }

    #pragma omp parallel
    {
        uint16_t thread_id = omp_get_thread_num();
        #pragma omp for reduction(+ : nnz_per_row_blocker[:num_blockers])
        for (uint16_t row_blocker_id = 0; row_blocker_id < num_blockers; ++row_blocker_id)
        {
            doRadixSort(global_blockers[row_blocker_id],
                        global_blockers[row_blocker_id] + global_blocker_counters[row_blocker_id],
                        sorting_buffer[thread_id]);
            int after = doMerge(global_blockers[row_blocker_id], global_blocker_counters[row_blocker_id]);
            nnz_per_row_blocker[row_blocker_id] += after;
        }
    }

    int *cumulative_row_indices = _malloc<int>(num_blockers + 1);
    scan(nnz_per_row_blocker, cumulative_row_indices, (int)(num_blockers) + 1);
    int total_nnz = cumulative_row_indices[num_blockers];

    if (C.isEmpty())
    {
        C.clear();
    }
    C.rows = A.rows;
    C.cols = B.cols;

    C.colids = static_cast<int *>(::operator new(sizeof(int[total_nnz])));
    C.rowptr = static_cast<int *>(::operator new(sizeof(int[C.rows + 1])));
    C.values = static_cast<T *>(::operator new(sizeof(T[total_nnz])));

    C.rowptr[0] = 0;

    #pragma omp parallel for
    for (uint16_t row_blocker_id = 0; row_blocker_id < num_blockers; ++row_blocker_id)
    {
        int base = cumulative_row_indices[row_blocker_id];
        TripleNode *this_blocker = global_blockers[row_blocker_id];
        for (IT i = 0; i < nnz_per_row_blocker[row_blocker_id]; ++i)
        {
            ++nnz_by_row[std::get<0>(this_blocker[i])];
            C.colids[base + i] = std::get<1>(this_blocker[i]);
            C.values[base + i] = std::get<2>(this_blocker[i]);
        }
    }

    scan(nnz_by_row, C.rowptr, C.rows + 1);
    C.nnz = total_nnz;

    _free<int>(flops_by_row_blockers);
    _free<int>(nnz_by_row);
    _free<int>(nnz_per_row_blocker);
    _free<int>(cumulative_row_indices);

    for (uint16_t row_blocker_id = 0; row_blocker_id < num_blockers; ++row_blocker_id)
    {
        _free<TripleNode>(global_blockers[row_blocker_id]);
    }
    _free<TripleNode *>(global_blockers);
    for (uint16_t thread_id = 0; thread_id < nthreads; ++thread_id)
    {
        _free<TripleNode>(local_blockers[thread_id]);
        _free<IT>(local_blocker_counters[thread_id]);
    }
    _free<TripleNode *>(local_blockers);
    _free<IT *>(local_blocker_counters);
}


int main(int argc, char *argv[])
{

    CSC<MATRIX_TYPE> A_csc, B_csc;
    inputToCSC(A_csc, argv[1]);
    inputToCSC(B_csc, argv[2]);

    CSR<MATRIX_TYPE> B_csr(B_csc);
    CSR<MATRIX_TYPE> C_csr;

    A_csc.sort();
    B_csr.sort();

    auto nflop = get_flop(A_csc, B_csr);

    BinSpGEMM(A_csc, B_csr, C_csr);

    C_csr.clear();
    A_csc.clear();
    B_csc.clear();
    B_csr.clear();

    return 0;
}