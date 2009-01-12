/* FunctionTest.h: Header file giving the definitions
   of the test case functions for 3D reconstruction */

#ifndef _TESTFUNCTIONS_3D_INCLUDED
#define _TESTFUNCTIONS_3D_INCLUDED

/********************************************************
 * Function Test_Default3D :                            *
 *          f(x,y,z) = z^6 - 4z^4 + 2z                  *
 *******************************************************/

double Test_Default3D (double x, double y, double z);
double Test_Default3D_Integral (double x1, double x2, double y1, double y2, double z1, double z2);

/********************************************************
 * Function Test_Example1:                              *
 *          f(x,y,z) = x
 *******************************************************/

double Test_Example1 (double x, double y, double z); 
double Test_Example1_Integral (double x1, double x2, double y1, double y2, double z1, double z2);

/**********************************************************************
 * Function Test_Example2:                                            *
 *          f(x,y,z) = y
 ********************************************************************+
*/

double Test_Example2 (double x, double y, double z);
double Test_Example2_Integral (double x1, double x2, double y1, double y2, double z1, double z2);

/********************************************************
 * Function Test_Example3:                              *
 *          f(x,y,z) = z
 *******************************************************/

double Test_Example3 (double x, double y, double z);
double Test_Example3_Integral (double x1, double x2, double y1, double y2, double z1, double z2);

/********************************************************
 * Function Test_Example4:                              *
 *          f(x,y,z) = sin(10x)
 *******************************************************/

double Test_Example4 (double x, double y, double z);
double Test_Example4_Integral (double x1, double x2, double y1, double y2, double z1, double z2);

/********************************************************
 * Function Test_Example5:                              *
 *          f(x,y,z) = sin(10y)
 *******************************************************/

double Test_Example5 (double x, double y, double z);
double Test_Example5_Integral (double x1, double x2, double y1, double y2, double z1, double z2);

/********************************************************
 * Function Test_Example6:                              *
 *          f(x,y,z) = sin(10z)
 *******************************************************/

double Test_Example6 (double x, double y, double z);
double Test_Example6_Integral (double x1, double x2, double y1, double y2, double z1, double z2);

/********************************************************
 * Function Test_Example7:                              *
 *          f(x,y,z) = step at y = 0.5
 *******************************************************/

double Test_Example7 (double x, double y, double z);
double Test_Example7_Integral (double x1, double x2, double y1, double y2, double z1, double z2);

/********************************************************
 * Function Test_Example8:                              *
 *          f(x,y,z) = cos(10*x)*cos(10*x) + 3*sin(10*y)
 *******************************************************/

double Test_Example8 (double x, double y, double z);
double Test_Example8_Integral (double x1, double x2, double y1, double y2, double z1, double z2);

/********************************************************
 * Function Test_Example9:                              *
 *          f(x,y,z) = x + y + z
 *******************************************************/

double Test_Example9 (double x, double y, double z);
double Test_Example9_Integral (double x1, double x2, double y1, double y2, double z1, double z2);

/*************************************************************************
 * Function Test_Example10:                                              *
 *          f(x,y,z) = cos(10*x)*cos(10*x) + 3*sin(10*y) + 0.5*cos(10*z)
 ************************************************************************/

double Test_Example10 (double x, double y, double z);
double Test_Example10_Integral (double x1, double x2, double y1, double y2, double z1, double z2);

/*************************************************************************
 * Function Test_Example11:                                              *
 *          f(x,y,z) = cos(10*x)*cos(10*x) + 3*sin(10*y) + 0.5*cos(10*z)
 ************************************************************************/

double Test_Example11 (double x, double y, double z);
double Test_Example11_Integral (double x1, double x2, double y1, double y2, double z1, double z2);

/*************************************************************************
 * Function Test_Example12:                                              *
 *          f(x,y,z) = cos(10*x)*cos(10*x) + 3*sin(10*y) + 0.5*cos(10*z)
 ************************************************************************/

double Test_Example12 (double x, double y, double z);
double Test_Example12_Integral (double x1, double x2, double y1, double y2, double z1, double z2);

#endif // _TESTFUNCTIONS_2D_INCLUDED

  
