/* FunctionTest.h: Header file defining
                   the forms of the test functions for reconstruction 1D and 2D*/

#ifndef _TESTFUNCTIONS_INCLUDED
#define _TESTFUNCTIONS_INCLUDED

#include "Math/Math.h"
#include <iostream>
#include <cmath>

/********************************************************
 * Function Test_Default1D :                            *
 *          f(x) = x^6 - 4x^4 + 2x                      *
 *******************************************************/

double Test_Default1D (double x);

double Test_Default1D_Integral (double x1, double x2);

double Test_Default1D_FirstDeriv (double x);

double Test_Default1D_SecondDeriv (double x);

double Test_Default1D_ThirdDeriv (double x);

double Test_Default1D_FourthDeriv (double x);

/********************************************************
 * Function Test_Example1:                              *
 *          f(x) = sin(2x)+2*cos(x)                     *
 *******************************************************/

double Test_Example1 (double x); 

double Test_Example1_Integral (double x1, double x2);

double Test_Example1_FirstDeriv (double x);

double Test_Example1_SecondDeriv (double x);

double Test_Example1_ThirdDeriv (double x);

double Test_Example1_FourthDeriv (double x);


/********************************************************
 * Function Test_Example2:                              *
 *          f(x) = exp(-4*x)*sin(5*x)                   *
 *******************************************************/

double Test_Example2 (double x);

double Test_Example2_g (double x);

double Test_Example2_Integral (double x1, double x2);

double Test_Example2_FirstDeriv (double x);

double Test_Example2_SecondDeriv (double x);

double Test_Example2_ThirdDeriv (double x);

double Test_Example2_FourthDeriv (double x);

/********************************************************
 * Function Test_Example3:                              *
 *          f(x) = step function                        *
 *******************************************************/

double Test_Example3 (double x);

double Test_Example3_Integral (double x1, double x2);

/********************************************************
 * Function Test_Example4:                              *
 *          f(x) = 2.0 + sin((x+1) * PI)                *
 *******************************************************/

double Test_Example4 (double x);

double Test_Example4_Integral (double x1, double x2);

/**********************************************************
 * Function Test_Example5:                                *
 *          f(x) -- test from Ollivier-Gooch paper        *
 *  Quasi-ENO Schemes for Unstructured Meshes Based on    *
 *  Unlimited Data-Dependent Least-Squares Reconstruction *
 *  Journal of Comp.Phy. 133 (1997)                       *
 *  In text is referred as Abgrall's function (pg. 13)    *
 *********************************************************/

double Test_Example5 (double x);

double Test_Example5_Integral (double x1, double x2);

/********************************************************
 * Function Test_Example6:                              *
 *          f(x) =  0.345*x + 0.5                       *
 *******************************************************/

double Test_Example6 (double x);

double Test_Example6_Integral (double x1, double x2);

/********************************************************
 * Function Test_Example7:                              *
 *          f(x) = sin(PI*x)                            *
 *******************************************************/

double Test_Example7 (double x);

double Test_Example7_Integral (double x1, double x2);

/********************************************************
 * Function Test_Example8:                              *
 *          f(x) =                    *
 *******************************************************/

double Test_Example8 (double x);

double Test_Example8_Integral (double x1, double x2);

/********************************************************
 * Function Test_Example9:                              *
 *          f(x) =                    *
 *******************************************************/

double Test_Example9 (double x);

double Test_Example9_Integral (double x1, double x2);


/********************************************************
 * Function Test_Example10:                              *
 *          f(x) =                    *
 *******************************************************/

double Test_Example10 (double x);

double Test_Example10_Integral (double x1, double x2);

double JIANG_Function (double x);

#endif // _TESTFUNCTIONS_INCLUDED

  
