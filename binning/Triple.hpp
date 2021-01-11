template <class T> 
struct Triple {
  Triple() : row(0), col(0), val(0){};
  Triple(int myrow, int mycol, T myval) : row(myrow), col(mycol), val(myval) {};

  int row; 
  int col;
  T val; 
};