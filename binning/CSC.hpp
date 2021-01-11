#include <cstdint>
#include <vector>
#include <utility>
#include <algorithm>
#include <functional>

#include "util.hpp"
#include "CSR.hpp"

template <class T> 
struct Triple {
  Triple() : row(0), col(0), val(0){};
  Triple(int myrow, int mycol, T myval) : row(myrow), col(mycol), val(myval) {};

  int row; 
  int col;
  T val; 
};

template <class T>
class CSC {
    public:
        // Constructors
        CSC(): nnz(0), rows(0), cols(0) {}
        CSC(Triple<T> *triples, int _nnz, int m, int n);    

        // Methods
        void sort();
        void clear();

        CSC<T> &operator = (const CSC<T> &rhs); // assignment operator
        ~CSC() {clear();} // Destructor

    // Attributes
    int rows;
    int cols;
    int nnz;

    int *colptr;
    int *rowids;
    T *values;
};


template <class T>
CSC<T>::CSC(Triple<T> *triples, int _nnz, int m, int n) : nnz(_nnz), rows(m), cols(n) {
    colptr = _malloc<int>(cols + 1);
    rowids = _malloc<int>(nnz);
    values = _malloc<T>(nnz);

    vector<pair<int, T>> tosort(nnz);

    int *nnz_per_col = _malloc<int>(cols);
    std::fill(nnz_per_col, nnz_per_col + cols, 0);

    // Count the number of entries per column.
    for (int k = 0; k < nnz; ++k) {
        int tmp = triples[k].col;
        nnz_per_col[tmp]++;
    }


    // Init values
    if (nnz > 0) {
        colptr[cols] = SumPrev(nnz_per_col, cols); 
        std::copy(nnz_per_col, nnz_per_col + cols, colptr);

        int last;
        for (int k = 0; k < nnz; ++k) {
            tosort[nnz_per_col[triples[k].col]++] = make_pair(triples[k].row, triples[k].val);
        }

        #pragma omp parallel for
        for (int i = 0; i < cols; ++i) {
            std::sort(tosort.begin() + colptr[i], tosort.begin() + colptr[i + 1]);

            typename vector<pair<int, T>>::iterator itr;
            int ind;
            for (itr = tosort.begin() + colptr[i], ind = colptr[i];
                itr != tosort.begin() + colptr[i + 1]; ++itr, ++ind) {
                rowids[ind] = itr->first;
                values[ind] = itr->second;
            }
        }
    }

    _free<int>(nnz_per_col);
}


template <class T> 
void CSC<T>::sort() {
    bool sorted = true;
    for (int i = 0; i < cols; ++i) {
        sorted &= _is_sorted(rowids + colptr[i], rowids + colptr[i + 1], std::less<T>());
    }
}


template <class T>
void CSC<T>::clear(){
    if (nnz > 0){
        _free<int>(rowids);
        _free<T>(values);
        nnz = 0;
    }
    if (cols > 0){
        _free<int>(colptr);
        cols = 0;
    }
    rows = 0;
}


template <class T>
CSC<T> &CSC<T>::operator=(const CSC<T> &rhs){
  if (this != &rhs) {
    if (nnz > 0) // if the existing object is not empty
    {
      _free<int>(rowids);
      _free<T>(values);
    }
    if (cols > 0) {
      _free<int>(colptr);
    }

    nnz = rhs.nnz;
    rows = rhs.rows;
    cols = rhs.cols;
    if (rhs.nnz > 0) // if the copied object is not empty
    {
      values = _malloc<T>(nnz);
      rowids = _malloc<int>(nnz);
      std::copy(rhs.values, rhs.values + nnz, values);
      std::copy(rhs.rowids, rhs.rowids + nnz, rowids);
    }
    if (rhs.cols > 0) {
      colptr = _malloc<int>(cols + 1);
      std::copy(rhs.colptr, rhs.colptr + cols + 1, colptr);
    }
  }
  return *this;
}

template <class T>
void inputToCSC(CSC<T> &A, char *input){
    A.clear();
    ifstream infile((*input).c_str());

    char line[256];
    char c = infile.get();

    // Skip comments
    while(c == '%')
    {
        infile.getline(line,256);
        c = infile.get();
    }
    infile.unget();

    // Scan first line.
    infile.getline(line,256);
    int m,n,nnz;
    sscanf(line, "%d %d %d", &m, &n, &nnz);
    
    // Store values.
    Triple<T> *triples = new Triple<T>[nnz];
    if (infile.is_open())
    {
        int cnz = 0; // current number of nonzeros
        while (!infile.eof() && cnz < nnz)
        {
            infile.getline(line,256);
            char *ch = line;
            triples[cnz].row = atoi(ch);
            ch = strchr(ch, ' ');
            ch++;
            triples[cnz].col = atoi(ch);
            ch = strchr(ch, ' ');
            if (ch != NULL) {
                ch++;
                /* Read third word (value data)*/
                triples[cnz].val = (T) (atoi(ch));
                ch = strchr(ch, ' ');
            }
            else {
                triples[cnz].val = 1.0;
            }

            // Account for 1-indexing.
            triples[cnz].row--;
            triples[cnz].col--;
            ++cnz;
        }
    }
    
    A = *(new CSC<T>(triples, nnz, m, n));
    delete [] triples;
    return 1;
}


template <class T>
long long int get_flop(const CSC<T> &A, const CSR<T> &B)
{
    long long int flops = 0;

#pragma omp parallel for reduction(+: flops)
    for (int i = 0; i < A.cols; ++i){
        int colnnz = A.colptr[i + 1] - A.colptr[i];
        int rownnz = B.rowptr[i + 1] - B.rowptr[i];
        flops += (colnnz * rownnz);
    }
    return (flops * 2);
}