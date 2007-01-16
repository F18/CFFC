/* FunctionTest.h: Header file defining
                   the forms of the test functions for reconstruction 2D*/

#ifndef _TESTFUNCTIONS_2D_INCLUDED
#define _TESTFUNCTIONS_2D_INCLUDED

#include "Grid/Grid2D/Cell2D.h"
#include "Math/NumericalLibrary.h"

using namespace std;

/********************************************************
 * Function Test_Default2D :                            *
 *          f(x) = x^6 - 4x^4 + 2x                      *
 *******************************************************/

double Test_Default2D (double x, double y);

double Test_Default2D_Integral (double x1, double x2, double y1, double y2);

/********************************************************
 * Function Test_Example1:                              *
 *          f(x) = x*x + y*y + x*y                      *
 *******************************************************/

double Test_Example1 (double x, double y); 

double Test_Example1_Integral (double x1, double x2, double y1, double y2);

/**********************************************************************
 * Function Test_Example2:                                            *
 *          f(x) = (x-0.25)*(x-0.5)*(x+0.25)*(y-0.24)*(y-1.5)*(y+0.25)*
 ********************************************************************+
*/

double Test_Example2 (double x, double y);

double Test_Example2_g (double x, double y);

double Test_Example2_Integral (double x1, double x2, double y1, double y2);

/********************************************************
 * Function Test_Example3:                              *
 *          f(x) = (x-0.5)*(y-0.25)+(y+2)^3             *
 *******************************************************/

double Test_Example3 (double x, double y);

double Test_Example3_Integral (double x1, double x2, double y1, double y2);

/********************************************************
 * Function Test_Example4:                              *
 *          f(x) =  x*(x-0.25)*(x+0.25)*pow(y,3)        *
 *******************************************************/

double Test_Example4 (double x, double y);

double Test_Example4_Integral (double x1, double x2, double y1, double y2);

/********************************************************
 * Function Test_Example5:                              *
 *          f(x) =                    *
 *******************************************************/

double Test_Example5 (double x, double y);

double Test_Example5_Integral (double x1, double x2, double y1, double y2);

/********************************************************
 * Function Test_Example6:                              *
 *          f(x) =                    *
 *******************************************************/

double Test_Example6 (double x, double y);

double Test_Example6_Integral (double x1, double x2, double y1, double y2);
// Node2D SW, Node2D SE, Node2D NE, Node2D NW)

/********************************************************
 * Function Test_Example7:                              *
 *          f(x) =                    *
 *******************************************************/

double Test_Example7 (double x, double y);

double Test_Example7_Integral (double x1, double x2, double y1, double y2);

/********************************************************
 * Function Test_Example8:                              *
 *          f(x) =                    *
 *******************************************************/

double Test_Example8 (double x, double y);

double Test_Example8_f (double r);

double Test_Example8_Integral (double x1, double x2, double y1, double y2);

double Test_Example9 (double x, double y);

/********************************************************
 * Function Test_Example10:                             *
 *          f(x,y) = exp((x+y)/(x-y))                   *
 *******************************************************/

double Test_Example10(double x, double y);

/********************************************************
 * Function Test_Example11:                             *
 *          f(x,y) = 1/50*exp(-x/10)*exp(-y/5)          *
 *******************************************************/

double Test_Example11(double x, double y);

/********************************************************
 * Function Test_Example12:                             *
 *          f(x,y) = 1.0                                *
 *******************************************************/

double Test_Example12(double x, double y);


/********************************************************
 * Function Test_Example12:                             *
 *          f(x,y) = x^3*y^5                            *
 *******************************************************/

double Test_Example13(double x, double y);


#endif // _TESTFUNCTIONS_2D_INCLUDED

  
