/*!\file ExactSolutions.h
  \brief Header file defining 2D advection-diffusion exact solutions. */

#ifndef _ADVECTDIFFUSE2D_EXACTSOLUTION_INCLUDED
#define _ADVECTDIFFUSE2D_EXACTSOLUTION_INCLUDED

/* Include required C++ libraries. */
#include <cmath>
#include <iostream>

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Math/Math.h"
#include "../Math/Vector2D.h"
#include "../Math/NumericalLibrary.h"
#include "../Utilities/TypeDefinition.h"
#include "../Utilities/Utilities.h"

/* Laplace_Solutions: this class defines solutions to the Laplace eqn: "(nabla)^2 w(x,y) = 0.0" */
class Laplace_Solutions{

 public:

  enum Accuracy_Measurement_Type {Soln=1, Grad=2 };
  static Accuracy_Measurement_Type Accuracy_Parameter;

  static double A,B,C,D, miu;



  // ============ INITIAL CONDITIONS ===============
  static double IC_1(double x, double y);
  static Vector2D Grad_IC_1(double x, double y);

  static double IC_2(double x, double y);
  static Vector2D Grad_IC_2(double x, double y);

  static double IC_3(double x, double y);
  static Vector2D Grad_IC_3(double x, double y);

  static double IC_4(double x, double y);
  static Vector2D Grad_IC_4(double x, double y);

  static double IC_5(double x, double y);
  static Vector2D Grad_IC_5(double x, double y);

  static double Result(double x,double y);

 private:
  Laplace_Solutions(){ };

};

// Definitions of member functions of Laplace_Solutions
inline double Laplace_Solutions::IC_1(double x, double y){
  return A*x + B*y + C;
}

inline Vector2D Laplace_Solutions::Grad_IC_1(double x, double y){
  return Vector2D(A,B);
}

inline double Laplace_Solutions::IC_2(double x, double y){
  return A*(x*x - y*y) + B*x*y;
}

inline Vector2D Laplace_Solutions::Grad_IC_2(double x, double y){
  return Vector2D(2.0*A*x + B*y, -2.0*A*y + B*x);
}

inline double Laplace_Solutions::IC_3(double x, double y){
  return A*(x*x*x - 3.0*x*y*y) + B*(3*x*x*y - y*y*y);
}

inline Vector2D Laplace_Solutions::Grad_IC_3(double x, double y){
  return Vector2D(3.0*A*(x*x - y*y) + 6.0*B*x*y , -6.0*A*x*y + 3.0*B*(x*x - y*y) ) ;
}

inline double Laplace_Solutions::IC_4(double x, double y){
  return exp(miu*x)*(A*cos(miu*y) + B*sin(miu*y));
}

inline Vector2D Laplace_Solutions::Grad_IC_4(double x, double y){
  return (miu*exp(miu*x))*Vector2D( A*cos(miu*y) + B*sin(miu*y), -A*sin(miu*y) + B*cos(miu*y));
}

inline double Laplace_Solutions::IC_5(double x, double y){
  return ( A*sinh(miu*x) + B*cosh(miu*x)) * (C*cos(miu*y) + D*sin(miu*y));
}

inline Vector2D Laplace_Solutions::Grad_IC_5(double x, double y){
  return miu*Vector2D( (A*cosh(miu*x) + B*sinh(miu*x)) * (C*cos(miu*y) + D*sin(miu*y)),
		       (A*sinh(miu*x) + B*cosh(miu*x)) *(-C*sin(miu*y) + D*cos(miu*y)) );
}

inline double Laplace_Solutions::Result(double x,double y){
  return 0.0;
}

/* Poisson_NonlinearSource_Solutions: this class defines solutions to the Poisson eqn with non-linear source term:
   "(nabla)^2 w(x,y) = f(w)", "f(w) = a*exp(beta*w)" */
class Poisson_NonlinearSource_Solutions{

 public:

  enum Accuracy_Measurement_Type {Soln=1, Grad=2 };
  static Accuracy_Measurement_Type Accuracy_Parameter;

  static double a, beta;
  static double A,B,C;

  /* pointer to the analytic solution: used for estimating the RHS of the eqn. */
  static FunctionType2D AnalyticSoln;


  // ============ INITIAL CONDITIONS ===============
  static double IC_1(double x, double y);
  static double IC_2(double x, double y);
  static double IC_3(double x, double y);
  static double IC_4(double x, double y);
  static double IC_5(double x, double y);
  static double GradX_IC_5(double x, double y);
  static double GradY_IC_5(double x, double y);
  static double Result(double x,double y);

 private:
  Poisson_NonlinearSource_Solutions(){ };

};

// Definitions of member functions of Poisson_NonlinearSource_Solutions
inline double Poisson_NonlinearSource_Solutions::IC_1(double x, double y){
  return log( (2.0*(A*A + B*B) )/( a*beta*(A*x + B*y + C)*(A*x + B*y + C) ) )/beta;
}

inline double Poisson_NonlinearSource_Solutions::IC_2(double x, double y){
  return log( (2.0*(A*A + B*B) )/( a*beta* sinh(A*x + B*y + C)* sinh(A*x + B*y + C) ) )/beta;
}

inline double Poisson_NonlinearSource_Solutions::IC_3(double x, double y){
  return log( (-2.0*(A*A + B*B) )/( a*beta* cosh(A*x + B*y + C)* cosh(A*x + B*y + C) ) )/beta;
}

inline double Poisson_NonlinearSource_Solutions::IC_4(double x, double y){
  return log( (2.0*(A*A + B*B) )/( a*beta* cos(A*x + B*y + C)* cos(A*x + B*y + C) ) )/beta;
}

inline double Poisson_NonlinearSource_Solutions::IC_5(double x, double y){
  return (log(8.0*C/a/beta) - 2.0*log(fabs( (x+A)*(x+A) + (y+B)*(y+B) - C )) )/beta;
}

inline double Poisson_NonlinearSource_Solutions::Result(double x,double y){
  return a*exp(beta* AnalyticSoln(x,y) );
}


inline double Poisson_NonlinearSource_Solutions::GradX_IC_5(double x, double y){

  double expression ((x+A)*(x+A) + (y+B)*(y+B) - C);

  return -2.0*(sign(expression) * 2.0 * (x + A) )/ beta / fabs(expression);

}
 
inline double Poisson_NonlinearSource_Solutions::GradY_IC_5(double x, double y){

  double expression ((x+A)*(x+A) + (y+B)*(y+B) - C);

  return -2.0*(sign(expression) * 2.0 * (y + B) )/ beta / fabs(expression);
}


/* Solutions to Diffusion Eqn.: this class defines solutions to the diffusion eqn with non-linear source term:
   "(nabla)^2 u(x,y) = f(u)", "f(u) = -lambda*u" */
class DiffusionEqn_NonlinearSource_Solutions{

 public:

  static double lambda, lambda_sqrt; 	/* lambda = -1/(kappa*tau) */
  static double CoordA, CoordB, SolnA, SolnB;	/* parameters for imposing BCs */
  static double C1, C2; 	/* particular solution constants, defined by the boundary conditions */

  // ============ SOLUTION FOR DIFFUSION IN 1D  ===============
  /* Assume Dirichlet BCs, having the form:
     u(CoordA) = SolnA  , the coordinate in the other direction is free
     u(CoordB) = SolnB  , 
  */

  /* X-DIRECTION */
  static double X_Diffusion(double x, double y);
  static Vector2D Grad_X_Diffusion(double x, double y);

  /* Y-DIRECTION */
  static double Y_Diffusion(double x, double y);
  static Vector2D Grad_Y_Diffusion(double x, double y);

  /* 1D-Solution (based on what existed in the old code) */
  static double Diffusion_Old(double Var1D);

  /* 1D-Solution derived by L. Ivan */
  static double Diffusion(double Var1D);
  
  /* Set particular solution constants */
  static void Set_ParticularSolution_Parameters(void);


 private:
  DiffusionEqn_NonlinearSource_Solutions(){ };

};

inline double DiffusionEqn_NonlinearSource_Solutions::X_Diffusion(double x, double y){
  return Diffusion(x);
}

inline double DiffusionEqn_NonlinearSource_Solutions::Y_Diffusion(double x, double y){
  return Diffusion(y);
}

inline double DiffusionEqn_NonlinearSource_Solutions::Diffusion_Old(double Var1D){
  double C1, C2;		/* integration constants */
  double sqld, detA;

  if (lambda > ZERO) {
    sqld = sqrt(lambda);
    detA = cos(sqld*CoordA)*sin(sqld*CoordB)-cos(sqld*CoordB)*sin(sqld*CoordA);
    C1 = ( SolnA*sin(sqld*CoordB)-SolnB*sin(sqld*CoordA))/detA;
    C2 = (-SolnA*cos(sqld*CoordB)+SolnB*cos(sqld*CoordA))/detA;
    return C1*cos(sqld * Var1D) + C2*sin(sqld * Var1D);

  } else if (lambda < ZERO) {
    sqld = sqrt(-lambda);
    detA = cosh(sqld*CoordA)*sinh(sqld*CoordB)-cosh(sqld*CoordB)*sinh(sqld*CoordA);
    C1 = ( SolnA*sinh(sqld*CoordB)-SolnB*sinh(sqld*CoordA))/detA;
    C2 = (-SolnA*cosh(sqld*CoordB)+SolnB*cosh(sqld*CoordA))/detA;
    return C1*cosh(sqld * Var1D)+C2*sinh(sqld * Var1D);

  } else { /* (lambda == ZERO) */ 
    detA = CoordB-CoordA;
    C1 = (SolnA*CoordB-SolnB*CoordA)/detA;
    C2 = (-SolnA+SolnB)/detA;
    return C1 + C2 * Var1D;
  } /* endif */
}


#endif
