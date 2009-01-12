/* FunctionTest.cc: Source file defining
                    the forms of the test functions for reconstruction 2D*/

/* Include header files. */

#ifndef _TESTFUNCTIONS_3D_INCLUDED
#include "TestFunctions_3D.h"
#endif // _TESTFUNCTIONS_3D_INCLUDED

#include <cmath>
#include <cassert>

/********************************************************
 * Function Test_Default1D :                            *
 *          f(x) = z^6 - 4z^4 + 2z                      *
 *******************************************************/

double Test_Default3D (double x, double y, double z){

  double f;
  f = x; f = y; 		// that's for the compiler!
  f = pow(z,6)-4*pow(z,4)+2*z;
  return f;
}

//this expression's only for cartesian mesh !!!!!!!!!!!!!!!!!
double Test_Default3D_Integral (double x1, double x2, double y1, double y2, double z1, double z2) {
  assert((x1<x2)&&(y1<y2)&&(z1<z2));
  double f;
  f = 1.0/7.0*(pow(z2,7) - pow(z1,7)) - 4.0/5.0*(pow(z2,5)-pow(z1,5)) + 
    (z2*z2 - z1*z1);
  f *= (x2 - x1)*(y2 - y1);
  return f;
}

/********************************************************
 * Function Test_Example1:                              *
 *          f(x,y,z) = z
 *******************************************************/

double Test_Example1 (double x, double y, double z) {
  double f;
  f = x; f = y;  	// that's for the compiler !
  f=z;
  return f;
}

double Test_Example1_Integral (double x1, double x2, double y1, double y2, double z1, double z2) {
  assert((x1<x2)&&(y1<y2)&&(z1<z2));
  double f;
  f = 0.5*(z2*z2 - z1*z1);
  return f;
}
/**********************************************************************
 * Function Test_Example2:
 *          f(x,y,z) = y
 *********************************************************************/

double Test_Example2 (double x, double y, double z) {

  double f;
  f = x; f = z;  	// that's for the compiler !
  f=y;
  return f;
}

double Test_Example2_Integral (double x1, double x2, double y1, double y2, double z1, double z2) {
  assert((x1<x2)&&(y1<y2)&&(z1<z2));
  double f;
  f = 0.5*(y2*y2 - y1*y1);
  return f;
}

/*********************************************************
 * Function Test_Example3:
 *          f(x,y,z) = x
 ********************************************************/

double Test_Example3 (double x, double y, double z) {

  double f;
  f = y; f = z;  	// that's for the compiler !
  f = x;
return f;
}

double Test_Example3_Integral (double x1, double x2, double y1, double y2, double z1, double z2) {
  assert((x1<x2)&&(y1<y2)&&(z1<z2));
  double f;
  f = 0.5*(x2*x2 - x1*x1);
  return f;
}

/********************************************************
 * Function Test_Example4:                              *
 *          f(x,y,z) = sin(z)
 ********************************************************/

double Test_Example4 (double x, double y, double z) {

  double f;
  f = x; f = y; 	// that's for the compiler !
  f = sin(z);
  return f;
}

double Test_Example4_Integral (double x1, double x2, double y1, double y2, double z1, double z2) {
  assert((x1<x2)&&(y1<y2)&&(z1<z2));
  double f;
  f = cos(z2)-cos(z1);
  return f;
}

/********************************************************
 * Function Test_Example5:                              *
 *          f(x,y,z) = sin(y)
 *******************************************************/

double Test_Example5 (double x, double y, double z) {
  double f;
  f = x; f = z;  	// that's for the compiler !
  f = sin(y);
  return f;
}

double Test_Example5_Integral (double x1, double x2, double y1, double y2, double z1, double z2) {
  assert((x1<x2)&&(y1<y2)&&(z1<z2));
  double f;
  f = cos(y2)-cos(y1);
  return f;
}

/********************************************************
 * Function Test_Example6:                              *
 *          f(x,y,z) = sin(x)
 *******************************************************/

double Test_Example6 (double x, double y, double z) {

  double f;
  f = y; f=z;  	// that's for the compiler !
  f = sin(x);
  return f;
}

double Test_Example6_Integral (double x1, double x2, double y1, double y2, double z1, double z2) {
  assert((x1<x2)&&(y1<y2)&&(z1<z2));
  double f;
  f = cos(x2)-cos(x1);
}

/********************************************************
 * Function Test_Example7:                              *

 *******************************************************/

double Test_Example7 (double x, double y, double z) {

  double f;
  f = cos(x)*cos(x) + 3*sin(y) + 0.5*cos(z);
  return f;
}

double Test_Example7_Integral (double x1, double x2, double y1, double y2, double z1, double z2) {
  assert((x1<x2)&&(y1<y2)&&(z1<z2));
  double f;
  f = -sin(2*x2) + sin(2*x1) + 3*cos(y2) - 3*cos(y1) - 0.5*sin(z2) + 0.5*sin(z1);
  return f;
}

/********************************************************
 * Function Test_Example8:                              *
 *          f(x) = step function at y = 0.5             *
 *******************************************************/

double Test_Example8 (double x, double y, double z)
{
  double f;
  double a = 0.5;
  //f = x;
  if (y <= a)
    f = 1.0e5;
  else 
    f = 1.0;
  return f;
}

double Test_Example8_Integral (double x1, double x2, double y1, double y2, double z1, double z2) {
  assert((x1<x2)&&(y1<y2)&&(z1<z2));
  double tol = 0.1e-6;
  double f;

  double a = 0.5;

  if (y2 <= (a-tol))
    f = (y2-y1)*Test_Example8(x2,y2,z2);
  else if ((y1<=a)&&(y2>a))
    f = (a-y1)*Test_Example8(x1,y1,z1) + (y2-a)*Test_Example8(x2,y2,z2);
  else
    f = (y2 - y1)*Test_Example8(x1,y1,z1);
  f *= (x2 - x1);
  return f;

}
