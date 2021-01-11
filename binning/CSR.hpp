#include <cstdint>
#include <vector>
#include <utility>
#include <algorithm>
#include <functional>

#include "CSC.hpp"
#include "util.hpp"

template <class T> 
class CSR {
public:
    CSR() : nnz(0), rows(0), cols(0) {}
    CSR(const CSC<T> &csc); // CSC -> CSR conversion
    CSR<T> &operator=(const CSR<T> &rhs); // assignment operator
    
    void clear();
    ~CSR() {clear();}
    bool isEmpty() { return (nnz == 0); }
    void sort();

    int rows;
    int cols;
    int nnz; // number of nonzeros

    int *rowptr;
    int *colids;
    T *values;
};

template <class T>
CSR<T>::CSR(const CSC<T> &csc): nnz(csc.nnz), rows(csc.rows), cols(csc.cols) {
    rowptr = _malloc<int>(rows + 1);
    colids = _malloc<int>(nnz);
    values = _malloc<T>(nnz);

    int *work = _malloc<int>(rows);
    std::fill(work, work + rows, 0); // initilized to zero
    for (int k = 0; k < nnz; ++k) {
        int tmp = csc.rowids[k];
        work[tmp]++; // row counts (i.e, w holds the "row difference array")
    }

    int last;

    if (nnz > 0) {
        rowptr[rows] = SumPrev(work, rows); // cumulative sum of work
        std::copy(work, work + rows, rowptr);

        for (int i = 0; i < cols; ++i) {
        for (int j = csc.colptr[i]; j < csc.colptr[i + 1]; ++j) {
            colids[last = work[csc.rowids[j]]++] = i;
            values[last] = csc.values[j];
        }
        }
    }
    _free<int>(work);
}

template <class T>
void CSR<T>::clear(){
    if (nnz > 0) {
        my_free<int>(colids);
        my_free<T>(values);
        nnz = 0;
    }
    if (rows > 0) {
        my_free<int>(rowptr);
        rows = 0;
    }
    cols = 0;
}


template <class T> 
void CSR<T>::sort() {
    bool sorted = true;
    for (int i = 0; i < cols; ++i) {
        sorted &= _is_sorted(colids + rowptr[i], colids + rowptr[i + 1], std::less<int>());
    }
}