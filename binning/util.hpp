#include <string>
#include <fstream>

using namespace std;


template <typename T> 
inline T *_malloc(int array_size) {

    T *a = new T[array_size];

    #pragma omp parallel for
    for(int i=0; i<array_size; i++){
        a[i] = T();
    }
    return a;
}

template <typename T> 
inline void _free(T *a) {
    delete [] a;
}


// Add the sum of previous entries to the current entry.
// Ex: SumPrev([1, 2, 3]) = [1, 3, 6]
int SumPrev(int *arr, int size) {
  int prev;
  int tempnz = 0;
  for (int i = 0; i < size; ++i) {
    prev = arr[i];
    arr[i] = tempnz;
    tempnz += prev;
  }
  return tempnz; // return sum
}


template <typename _ForwardIterator, typename _StrictWeakOrdering>
bool _is_sorted(_ForwardIterator first, _ForwardIterator last, _StrictWeakOrdering comp) {
    if (first == last)
        return true;

    _ForwardIterator next = first;
    for (++next; next != last; first = next, ++next)
        if (comp(*next, *first)) return false;
    return true;
};
