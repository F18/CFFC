/*!\file ExactSolutions.cc
  \brief Implementation of subroutines prototyped in ExactSolutions.h */

#include "../../Math/NumericalLibrary.h"
#include "ExactSolutions.h"

/***********************************************************
 * Routine: SIN_WAVE_Solution                              *
 *                                                         *
 * Define an initial density variation having a sinusoidal *
 * shape.                                                  *
 *                                                         *
 **********************************************************/
double SIN_WAVE_Solution (double x){
  return (2.0 + sin((x+1)*PI));
}

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
double Integral_SIN_WAVE(double xmin, double xmax, double x){

  double integral;

  integral = cos(-2*PI*x/(xmin-xmax)+2*PI*xmin/(xmin-xmax));
  integral *= 0.5*(xmin-xmax)/PI;
  integral += 2*x;

  return integral; 
}

/********************************************************
 * Routine: JIANG_IVP_Solution                          *
 *                                                      *
 * Define an initial density variation given by the     *
 * function proposed by Jiang.                          *
 *                                                      *
 ********************************************************/
double JIANG_IVP_Solution (double x){

  /*Obs: The function is shifted on "X" and "Y"
    X-shift -> in order to position the singular points inside the interval
    Y-shift -> in order to use this data distribution as initial condition for density in Euler1D
  */

  return 6.0 + JIANG_Function(x-0.5);
}

double JIANG_Function (double x){

  double point1, point2, point3, point4, point5;
  point1 = -1.0;
  point2 = - 1.0/3.0;
  point3 =  0.0;
  point4 = 1.0/3.0;
  point5 = 1.0;

  if ( x >= -1.5 && x < point1){
    x += 2.0; // periodic function
  }

  double solution = -(sqrt(3.0)/2.0 + 4.5 + 2.0*PI/3.0)*(x+1);

  if( point1 <= x && x < point2 ) {
    solution += 2.0*cos(3.0*PI*x*x/2.0) - sqrt(3.0);
    
  } else if( point2 <= x && x < point3 ) {
    solution += 3.0/2.0 + 3.0*cos(2.0*PI*x);

  } else if( point3 <= x && x < point4 ){
    solution += 15.0/2.0 - 3.0*cos(2.0*PI*x);

  } else if( point4 <= x && x <= point5 ) {
    solution += (28.0 + 4*PI + cos(3.0*PI*x))/3.0 + 6*PI*x*(x-1);

  } else {
    solution = 0.0;
  } //endif

  return solution;
}

/********************************************************
 * Routine:  ConvectionShapes                           *
 *                                                      *
 * Define an initial density variation which formed by  *
 * different shapes (e.g triangle, rectangle etc.)      *
 *                                                      *
 ********************************************************/
double ConvectionShapes(double x)
{

  double a = 0.5;
  double z = -0.7;
  double delta = 0.005;
  double alpha = 10;
  double beta = log10(2.0)/(36*delta*delta);

  if ((x>=-0.8)&&(x<=-0.6)){
    return 1.0 + (1.0/6.0)*(ConvectionShapes_HelperFunction_G(x,beta,z-delta)+
			 ConvectionShapes_HelperFunction_G(x,beta,z+delta)+
			 4*ConvectionShapes_HelperFunction_G(x,beta,z));
  }

  if ((x>=-0.4)&&(x<=-0.2)){
    return 1.0 + 1.0;
  }

  if ((x>= 0.0)&&(x<= 0.2)){
    return 1.0 + 1.0-fabs(10*(x-0.1));
  }

  if ((x>= 0.4)&&(x<=0.6)){
    return 1.0 + (1.0/6.0)*(ConvectionShapes_HelperFunction_F(x,alpha,a-delta)+
			    ConvectionShapes_HelperFunction_F(x,alpha,a+delta)+
			    4*ConvectionShapes_HelperFunction_F(x,alpha,a));
  }

  return 1.0;			// 1.0 is introduced in order to avoid density to be zero
}


double ConvectionShapes_HelperFunction_G (double x, double beta, double z)
{
  // see 'Efficient Implementation of Weighted ENO Schemes', JCP 126, 202-228
  return exp(-beta*(x-z)*(x-z));
}

double ConvectionShapes_HelperFunction_F (double x, double alpha, double a)
{
  // see 'Efficient Implementation of Weighted ENO Schemes', JCP 126, 202-228
  return sqrt(max(1-alpha*alpha*(x-a)*(x-a),0.0));
}
