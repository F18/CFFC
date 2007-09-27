#include "AdvectDiffuse2DExactSolutions.h"
#include "../Utilities/Utilities.h"

// Initialization of static member variables in Laplace_Solutions class
double Laplace_Solutions::A = 2.0;
double Laplace_Solutions::B = 3.0;
double Laplace_Solutions::C = 4.0;
double Laplace_Solutions::D = 5.0;
double Laplace_Solutions::miu = 0.25;
Laplace_Solutions::Accuracy_Measurement_Type Laplace_Solutions::Accuracy_Parameter = Soln;


// Initialization of static member variables in Poison_NonlinearSource_Solutions class
double Poisson_NonlinearSource_Solutions::A = 2.0;
double Poisson_NonlinearSource_Solutions::B = 3.0;
double Poisson_NonlinearSource_Solutions::C = 4.0;
double Poisson_NonlinearSource_Solutions::a = 0.15;
double Poisson_NonlinearSource_Solutions::beta = 0.25;
FunctionType2D Poisson_NonlinearSource_Solutions::AnalyticSoln = Poisson_NonlinearSource_Solutions::IC_1;
Poisson_NonlinearSource_Solutions::Accuracy_Measurement_Type Poisson_NonlinearSource_Solutions::Accuracy_Parameter = Soln;


// Initialization of static member variables in DiffusionEqn_NonlinearSource_Solutions class
double DiffusionEqn_NonlinearSource_Solutions::lambda = 1.0;
double DiffusionEqn_NonlinearSource_Solutions::lambda_sqrt = 1.0;
double DiffusionEqn_NonlinearSource_Solutions::C1 = 1.0;
double DiffusionEqn_NonlinearSource_Solutions::C2 = 1.0;
double DiffusionEqn_NonlinearSource_Solutions::CoordA = 1.0;
double DiffusionEqn_NonlinearSource_Solutions::CoordB = 1.0;
double DiffusionEqn_NonlinearSource_Solutions::SolnA = 1.0;
double DiffusionEqn_NonlinearSource_Solutions::SolnB = 1.0;

double DiffusionEqn_NonlinearSource_Solutions::Diffusion(double Var1D){
  /* assume lambda > ZERO */
  return C1*exp(lambda_sqrt * Var1D) + C2*exp(-lambda_sqrt * Var1D);
}

Vector2D DiffusionEqn_NonlinearSource_Solutions::Grad_X_Diffusion(double x, double y){
  return Vector2D(lambda_sqrt*(C1*exp(lambda_sqrt*x) - C2*exp(-lambda_sqrt*x)),
		  0.0);
}

Vector2D DiffusionEqn_NonlinearSource_Solutions::Grad_Y_Diffusion(double x, double y){
  return Vector2D(0.0,
		  lambda_sqrt*(C1*exp(lambda_sqrt*y) - C2*exp(-lambda_sqrt*y)) );
}

void DiffusionEqn_NonlinearSource_Solutions::Set_ParticularSolution_Parameters(void){
  double V1, V2, V3, V4, detV;
  
  /* assume lambda has already been set and is greater than ZERO */
  lambda = -lambda;

  /* set lambda_sqrt */
  lambda_sqrt = sqrt(lambda);

  /* compute the LHS coefficient matrix */
  V1 = exp( lambda_sqrt * CoordA);
  V2 = exp(-lambda_sqrt * CoordA);
  V3 = exp( lambda_sqrt * CoordB);
  V4 = exp(-lambda_sqrt * CoordB);

  /* compute the determinent of LHS */
  detV = V1*V4 - V3*V2;

  /* compute the constants */
  C1 = (V4*SolnA - V2*SolnB)/detV;
  C2 = (V1*SolnB - V3*SolnA)/detV;
}

