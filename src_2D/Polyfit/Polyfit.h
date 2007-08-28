/* C interface for Fortran DPOLFT subroutine */

#ifndef _POLYFIT_H
#define _POLYFIT_H

#include <vector>
#include <iostream>

using namespace std;


extern "C" {
  int dpolft_(int* n, double* x, double* y, double* w, 
	      int* maxdeg, int* ndeg, double* eps, double* r,
	      int* ierr, double* a);

  int dpcoef_(int* l, double* c, double* tc, double* a);
}



double polyfit(int n, double* x, double* y, double* w,
	       int maxdeg, int& ndeg, double eps, double* r);



template<class D, class R>
  R poly6(D x, R* c) {
  return ((((((c[6]*x + c[5])*x + c[4])*x + c[3])*x + 
	    c[2])*x + c[1])*x + c[0]);
}

template<class D, class R>
  R poly8(D x, R* c) {
  return ((((((((c[8]*x + c[7])*x + c[6])*x + c[5])*x + c[4])*x + c[3])*x + 
	    c[2])*x + c[1])*x + c[0]);
}

template<class D, class R>
  R poly10(D x, R* c) {
  return ((((((((((c[10]*x + c[9])*x + c[8])*x + c[7])*x 
		+ c[6])*x + c[5])*x + c[4])*x + c[3])*x 
	    + c[2])*x + c[1])*x + c[0]);
}
    
template<class D, class R>
  R poly5(D x, R* c) {
  return (((((c[5]*x + c[4])*x + c[3])*x + 
	    c[2])*x + c[1])*x + c[0]);
}
    
template<class D, class R>
  R poly4(D x, R* c) {
  return ((((c[4]*x + c[3])*x + 
	    c[2])*x + c[1])*x + c[0]);
}
    
template<class D, class R>
  R poly3(D x, R* c) {
  return (((c[3]*x + c[2])*x + c[1])*x + c[0]);
}
    
#endif

 
