/* C interface for Fortran DPOLFT subroutine */

#ifndef _POLYFIT_H
#define _POLYFIT_H

#include <vector>
#include <cmath>
#include <iostream>

using namespace std;


#include "../../../bpkit/src/arch.h"

extern "C" {
  int F77NAME(dpolft)(int* n, double* x, double* y, double* w, 
	      int* maxdeg, int* ndeg, double* eps, double* r,
	      int* ierr, double* a);

  int F77NAME(dpcoef)(int* l, double* c, double* tc, double* a);
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

template<class D, class R>
  R polyN(D x, R* c, int N) {
    
    double polynomial = c[0];
      for(int i=1; i<=N; i++){
          double temp = 1;
          for (int j=0; j<i; j++)
              temp *= x;
          polynomial += c[i]*temp;
      }
    return polynomial;
}


/*!
 * \brief polyfit_smoothing : Smooths the given data by using polynomials.
 *
 * For each data point, a polynomial will be fitted using the points around
 * that data point lying in a window around the point. The smoothed data point
 * is then the point of the polynomial with the same x-coordinate.
 *
 * \param [in] N        The number of data points in the data
 * \param [in] x        The x-coordinates of the unsmooth data
 * \param [in] y        The y-coordinates of the unsmooth data
 * \param [in] order    The degree of the polynomials used
 * \param [in] window   Defines the amount of points used around the center of
 *                      each polynomial
 * \param [in] window_flag   true: window is defined in physical space \n
 *                           false: window is defined in indexes
 * \param [out] y_smooth The smoothed data
 * \return the maximum eps as defined in polyfit
 */
double polyfit_smoothing(int N, double *x, double *y, int order, double window, bool window_flag, double *y_smooth) ;

/*!
 * \brief polyfit_smoothing : Smooths the given data by using polynomials.
 *
 * For each data point, a polynomial will be fitted using the points around
 * that data point lying in a window around the point. The smoothed data point
 * is then the point of the polynomial with the same x-coordinate.
 *
 * \param [in] N        The number of data points in the data
 * \param [in] x        The x-coordinates of the unsmooth data
 * \param [in] y        The y-coordinates of the unsmooth data
 * \param [in] order    The degree of the polynomials used
 * \param [in] window   Defines the amount of points used around the center of
 *                      each polynomial. In case of a physical window, it must 
 *                      be given logscale 10:   1.0 for 1 decade window
 * \param [in] window_flag   true: window is defined in physical space \n
 *                           false: window is defined in indexes
 * \param [out] y_smooth The smoothed data
 * \return the maximum eps as defined in polyfit
 */
double polyfit_smoothing_logscale(int N, double *x, double *y, int order, double window, bool window_flag, double *y_smooth) ;
#endif

 
