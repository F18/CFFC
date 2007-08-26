/********************************************************************
 * Program:  Spline1D.h                                             *
 *                                                                  *
 * Description:  This program generates a cubic spline interpolant. *
 *               See:                                               *   
 *                 Numerical Recipes in C, Second Edition (1992).   *
 *                 Available at http://www.nr.com/                  *
 *                                                                  *
 ********************************************************************/
#ifndef _SPLINE1D_INCLUDED
#define _SPLINE1D_INCLUDED

//--------------------------------------------------//
// SAMPLE USAGE
//--------------------------------------------------//
// double f(const double x) {
//   double y;
//   if (x<1.0) y = -1.0+4.0*x-pow(x,2);
//   else if (x>=1.0 && x<2.0) y = 2.0*x;
//   else y = 2.0-x+pow(x,2.0);
//   return y;
// }

// double f1(const double x) {
//   return sin(x);
// }


// int main() {

//   const int n = 11;
//   double a = 0.0, b = 10.0;
//   double x[n], y[n], y2[n], dx;
//   double yp0 = 1.0e30, ypn = 1.0e30;
//   double xi, yi, dxi;

  
//   dx = (b-a)/double(n-1);
//   for (int i=0; i<n; i++){
//     x[i] = a + double(i)*dx;
//     y[i] = f1(x[i]);
//   }

//   spline( x, y, n, yp0, ypn, y2 );


//   dxi = (b-a)/double(100-1); 
//   for (int i=0; i<100; i++){
//     xi = a + dxi*double(i);
//     splint( x, y, y2, n, xi, yi );
//     cout << setw(18) << xi << setw(18) << yi << setw(18) << f1(xi) << endl;
//   }

//   return 0;
// }
//--------------------------------------------------//
// END SAMPLE USAGE
//--------------------------------------------------//


/* Include required C++ libraries. */

#include <math.h>
#include <iostream>
#include <iomanip>
using namespace std;


/********************************************************************
 * spline
 *
 * Given arrays x[1..n] and y[1..n], containing a tabulated function,
 * i.e., y(i)=f(x(i)), with x1<x2<..<xn, and given values of yp1 
 * and ypn for the first derivative of the interpolating function 
 * at points 1 and n, respectively, this routine returns an array
 * y2[1..n] that contains the second derivatives of the 
 * interpolating function at the tabulated points x(i).  If yp1 
 * and/or ypn are equal to 1E+30 or larger, the routine is signaled 
 * to set the corresponding boundary condition for a natural spline,
 * with zero second derivative on that boundary.
 ********************************************************************/
inline void spline( const double *x, const double *y, const int n, 
		    const double yp0, const double ypn, double *y2)
{
  // declares
  int i, k;
  double p, qn, sig, un, *u;

  // allocate temporary storage
  u = new double[n];


  // the lower boundary condition is set either to be "natural 
  if (yp0>=1.0e30) {
    y2[0] = u[0] = 0.0;

  // or else to have a specified first derivative
  } else {
    y2[0] = -0.5;
    u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp0);
  } /* endif */


  // this is the decomposition loop of the tridiagonal algorithm.
  // y2 and u are used for temporary storage of the decomposed 
  // factors.
  for(i=1; i<=n-2; i++) {
    sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p = sig*y2[i-1]+2.0;
    y2[i] = (sig-1.0)/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  } /* endfor */


  // the upper boundary condition is set either to be "natural"
  if (ypn>=1.0e30) {
    qn = un = 0.0;

  // or else to have a specified first derivative
  } else {
    qn = 0.5;
    un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  } /* endif */


  // this is the backsubstitution loop of the tridiagonal algorithm
  y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k=n-2; k>=0; k--) y2[k] = y2[k]*y2[k+1]+u[k];


  // clean up memory
  delete[] u;  u = NULL;

}


/********************************************************************
 * splint
 *
 * Given the arrays xa[1..n] and ya[1..n], which tabulate a 
 * function (with the xa(i)'s in order), and given the array 
 * y2a[1..n], which is the output from 'spline' above, and given a 
 * value of x, this routine returns a cubic-spline interpolated 
 * value.
 ********************************************************************/
inline void splint( const double *xa, const double *ya, const double *y2a, 
		    const int n, const double x, double &y)
{
  // declares
  int klo, khi, k;
  double h, b, a;

  // We will find the right place in the table by means of bisection.
  // This is optimal if sequential calls to this routine are at random 
  // values of x.  If sequential calls are in order, and closely spaced,
  // one would do better to store previous values of klo and khi and
  // test if they remain appropriate on the next call.
  klo = 0;
  khi = n-1;
  while (khi-klo>1) {
    k = (khi+klo) >> 1;
    if (xa[k]>x) khi = k;
    else klo = k;
  } /* endwhile */
  // klo and khi now bracket the input value of x


  // the xa's must be distinct
  h = xa[khi]-xa[klo];
  if (h==0.0) {
    cerr << "Bad xa input to routine splint\n";
    exit(-1);
  }

  // evaluate cubic spline polynomial is now evaluated
  a = (xa[khi]-x)/h;
  b = (x-xa[klo])/h;
  y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}



/********************************************************************
 * splint
 *
 * Given the arrays xa[1..n] and ya[1..n], which tabulate a 
 * function (with the xa(i)'s uniformly spaced by dxa in order), 
 * and given the array  y2a[1..n], which is the output from 'spline' 
 * above, and given a value of x, this routine returns a cubic-spline 
 * interpolated value.
 ********************************************************************/
inline void splint( const double *xa, const double *ya, const double *y2a, 
		    const double dxa, const int n, const double x, double &y)
{
  // declares
  int klo, khi, k;
  double h, b, a;

  // find the right place in the table
  // This is simple for uniformly spaced values.
  if (x<=xa[0]) {
    klo = 0;
  } else if (x>=xa[n-1]){
    klo = n-2;
  } else {
    klo = int( (x-xa[0])/dxa );
  }
  khi = klo+1;
  // klo and khi now bracket the input value of x


  // the xa's must be distinct
  h = xa[khi]-xa[klo];
  if (h==0.0) {
    cerr << "Bad xa input to routine splint\n";
    exit(-1);
  }

  // evaluate cubic spline polynomial is now evaluated
  a = (xa[khi]-x)/h;
  b = (x-xa[klo])/h;
  y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}



#endif //_SPLINE1D_INCLUDED
