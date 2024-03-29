/*!\file TestFunctions_3D.cc
  \brief Source file to implement the test functions prototyped in TestFunctions_3D.h. */

/* Include header files. */
#include <cmath>
#include <cassert>

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestFunctions_3D.h"

/********************************************************
 * Function Test_Default1D :                            *
 *          f(x) = z^6 - 4z^4 + 2z                      *
 *******************************************************/

double Test_Default3D (double x, double y, double z){

  double f;
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
 *          f(x,y,z) = 
 *******************************************************/

double Test_Example1 (double x, double y, double z) {
  double f(0.0);

  return f;
}

/**********************************************************************
 * Function Test_Example2:
 *          f(x,y,z) = 
 *********************************************************************/

double Test_Example2 (double x, double y, double z) {

  double f(0.0);

  return f;
}

/*********************************************************
 * Function Test_Example3:
 *          f(x,y,z) = 
 ********************************************************/

double Test_Example3 (double x, double y, double z) {

  double f(0.0);

  return f;
}

/********************************************************
 * Function Test_Example4:                              *

 *******************************************************/

double Test_Example4 (double x, double y, double z) {

  double f(0.0);
  return f;
}

/********************************************************
 * Function Test_Example5:                              *

 *******************************************************/

double Test_Example5 (double x, double y, double z) {
  double f(0.0);
  return f;
}

/********************************************************
 * Function Test_Example6:                              *

 *******************************************************/

double Test_Example6 (double x, double y, double z) {
  double f(0.0);
  return f;
}

/********************************************************
 * Function Test_Example7:                              *

 *******************************************************/

double Test_Example7 (double x, double y, double z) {
  double f(0.0);
  return f;
}

/********************************************************
 * Function Test_Example8:                              *
 *          f(x) =                    *
 *******************************************************/

double Test_Example8 (double x, double y, double z)
{
  double f(0.0);
  return f;
}
