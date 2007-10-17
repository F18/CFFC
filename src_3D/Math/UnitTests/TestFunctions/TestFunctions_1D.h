/*!\file TestFunctions_1D.h
  \brief Header file defining the prototype of 1D test functions */

#ifndef _TESTFUNCTIONS_1D_INCLUDED
#define _TESTFUNCTIONS_1D_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
// None


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
 *          f(x) = 2 steps of different intensities     *
 *******************************************************/

double Test_Example8 (double x);

double Test_Example8_Integral (double x1, double x2);

/********************************************************
 * Function Test_Example9:                              *
 *          f(x) = sin(x*log(x))                        *
 *******************************************************/

double Test_Example9 (double x);

/********************************************************
 * Function Test_Example10:                             *
 *          f(x) = exp(x)                               *
 *******************************************************/

double Test_Example10 (double x);

double Test_Example10_Integral (double x1, double x2);

/********************************************************
 * Function Test_Example11:                             *
 *          f(x) =  sqrt(x)                             *
 *******************************************************/

double Test_Example11(double x);

double Test_Example11_Integral (double x1, double x2);

/********************************************************
 * Function Test_Example12:                             *
 *          f(x) =  1.0/(pow(x,4) + x*x + 0.9)          *
 *******************************************************/

double Test_Example12 (double x);

double Test_Example12_Integral (double x1, double x2);

/********************************************************
 * Function Test_Example13:                             *
 *          f(x) =  (23.0/25.0)*cosh(x) - cos(x)        *
 *******************************************************/

double Test_Example13 (double x);

double Test_Example13_Integral (double x1, double x2);

/********************************************************
 * Function Test_Example14:                             *
 *          f(x) =  sqrt(x*x*x)                         *
 *******************************************************/

double Test_Example14 (double x);

double Test_Example14_Integral (double x1, double x2);

/********************************************************
 * Function Test_Example15:                             *
 *          f(x) =  1.0/(1.0 + x)                       *
 *******************************************************/

double Test_Example15 (double x);

double Test_Example15_Integral (double x1, double x2);

/********************************************************
 * Function Test_Example16:                             *
 *          f(x) =  1.0/(1.0 + exp(x))                  *
 *******************************************************/

double Test_Example16 (double x);

double Test_Example16_Integral (double x1, double x2);

/********************************************************
 * Function Test_Example17:                             *
 *          f(x) =  1.0/(1.0 + pow(x,4))                *
 *******************************************************/

double Test_Example17 (double x);

double Test_Example17_Integral (double x1, double x2);

double Test_Example18 (double x);

double Test_Example18_Integral (double x1, double x2);

#endif // _TESTFUNCTIONS_INCLUDED

  
