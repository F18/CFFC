/*!file ExactSolutions.h
  \brief Subroutines defining exact solutions for 1D Euler equations, including density convection.*/

#ifndef _EXACTSOLUTIONS_EULER1D_INCLUDED
#define _EXACTSOLUTIONS_EULER1D_INCLUDED

/***********************************************************
 * Routine: SIN_WAVE_Solution                              *
 *                                                         *
 * Define an initial density variation having a sinusoidal *
 * shape.                                                  *
 *                                                         *
 **********************************************************/
double SIN_WAVE_Solution (double x);

/********************************************************
 * Routine: Integral_SIN_WAVE                           *
 *                                                      *
 * Determines the value of the integral of:             *
 *                                                      *
 * rho(x,0) = 2 + sin(2*Pi*(x-xmin)/(xmax-xmin))        *
 *                                                      *
 * in point "x"                                         *
 *                                                      *
 ********************************************************/
double Integral_SIN_WAVE(double xmin, double xmax, double x);


/********************************************************
 * Routine: JIANG_IVP_Solution                          *
 *                                                      *
 * Define an initial density variation given by the     *
 * function proposed by Jiang.                          *
 *                                                      *
 ********************************************************/
double JIANG_IVP_Solution (double x);

double JIANG_Function(double x);


/****************************************************************************
 * Routine:  ConvectionShapes                                               *
 *                                                                          *
 * Define an initial density variation which formed by                      *
 * different shapes (e.g triangle, rectangle etc.)                          *
 *                                                                          *
 * see 'Efficient Implementation of Weighted ENO Schemes', JCP 126, 202-228 *
 *                                                                          *
 ***************************************************************************/
double ConvectionShapes(double x);

double ConvectionShapes_HelperFunction_G (double x, double beta, double z);

double ConvectionShapes_HelperFunction_F (double x, double alpha, double a);

#endif
