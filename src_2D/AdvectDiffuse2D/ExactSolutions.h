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
#include "../CFD/CFD.h"
#include "AdvectDiffuse2DInflowField.h"


/*******************************************************************
 *                                                                 *
 *      Advection-diffusion Exact Solution Type Indexes            *
 *                                                                 *
 *******************************************************************/
#define AD2D_NO_EXACT_SOLUTION              -1
// ========== PARTICULAR SOLUTIONS FOR LAPLACE EQUATION ===========
#define AD2D_EXACT_SOLUTION_LAPLACE_I        0
#define AD2D_EXACT_SOLUTION_LAPLACE_II       1
#define AD2D_EXACT_SOLUTION_LAPLACE_III      2
#define AD2D_EXACT_SOLUTION_LAPLACE_IV       3
#define AD2D_EXACT_SOLUTION_LAPLACE_V        4
// ========== PARTICULAR SOLUTIONS FOR POISSON EQUATION ===========
#define AD2D_EXACT_SOLUTION_POISSON_I        10
#define AD2D_EXACT_SOLUTION_POISSON_II       11
#define AD2D_EXACT_SOLUTION_POISSON_III      12
#define AD2D_EXACT_SOLUTION_POISSON_IV       13
#define AD2D_EXACT_SOLUTION_POISSON_V        14
// ========== PARTICULAR SOLUTIONS FOR PURE ADVECTION ===========
#define AD2D_EXACT_SOLUTION_PURE_CIRCULAR_ADVECTION_AT_CONSTANT_SPIN  20
#define AD2D_EXACT_SOLUTION_SINUSOIDAL_WAVE                           21
#define AD2D_EXACT_SOLUTION_ABGRAL_FUNCTION                           22
#define AD2D_EXACT_SOLUTION_SINUSOIDAL_EXPONENTIAL_WAVE               23
#define AD2D_EXACT_SOLUTION_SINUSOIDAL_EXPONENTIAL_ROTATED_WAVE       24
#define AD2D_EXACT_SOLUTION_COSINE_HILL                               25
#define AD2D_EXACT_SOLUTION_HYPER_TANGENT                             26
#define AD2D_EXACT_SOLUTION_MODULATED_SINUSOIDAL_WAVE                 27
// ========== PARTICULAR SOLUTIONS FOR ADVECTION-DIFFUSION ===========
#define AD2D_EXACT_SOLUTION_ADVECTION_DIFFUSION_IN_ANNULUS  30
#define AD2D_EXACT_SOLUTION_ADVECTION_DIFFUSION_IN_RECTANGULAR_CHANNEL  31
// ========== PARTICULAR SOLUTIONS FOR STATIONARY HEAT TRANSFER EQUATION WITH SOURCE ===========
#define AD2D_EXACT_SOLUTION_STATIONARY_HEAT_TRANSFER_WITH_LINEAR_SOURCE   40




// Declare the input parameters class
class AdvectDiffuse2D_Input_Parameters;

/*! 
 * \class ExactSolutionBasicType
 *
 * \brief Basic data type for any exact solution type.
 *
 * This is an abstract data type (ADT).
 */
class ExactSolutionBasicType{
public:

  //! @name Public types
  //@{
  enum Accuracy_Measurement_Type {Soln=1, Grad=2 }; // Solution or gradient
  //@}

  //! Default ctor
  ExactSolutionBasicType(void);

  //! @name Virtual member functions
  //@{
  //! Declare a pure virtual destructor
  virtual ~ExactSolutionBasicType(void) = 0;

  /*! Calculate the exact solution based on the location of interest */
  virtual double EvaluateSolutionAt(const double &x, const double &y) const = 0;

  /*! Calculate the exact solution gradient based on the location of interest */
  virtual Vector2D EvaluateGradientAt(const double &x, const double &y) const {
    throw runtime_error("EvaluateGradientAt() ERROR! The current exact solution doesn't have the gradient defined");
  }

  /*! Calculate the right hand side (RHS) of the partial differential equation (PDE) at the location of interest */
  virtual double PDE_RightHandSide(const double &x, const double &y) const  = 0;

  //! Update internal variables
  virtual void Set_ParticularSolution_Parameters(void){ };

  //! @name Functions for input-output and broadcast
  //@{
  virtual void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command) = 0;
  virtual void Print_Info(std::ostream & out_file) = 0;
  virtual void Broadcast(void) = 0;
  //@}
  //@}

  const string & whatExactSolution(void) const { return ExactSolutionName; }   //!< Get the name of the exact solution type

protected:
  string ExactSolutionName;		//!< storage for the name of the source field type
  Accuracy_Measurement_Type Accuracy_Parameter;	//!< indicator for the parameter used to estimate the accuracy
                                                // (e.g solution or gradient)

};


/*! 
 * \class Laplace_I_ExactSolution
 * 
 * \brief Implements a particular exact solution to the Laplace equation: 
 *        \f$ \frac{\partial^2 w}{\partial x^2} + \frac{\partial^2 w}{\partial y^2} = 0 \f$
 *
 * The solution has the expression: \f$ w(x,y) = A*x + B*y + C \f$
 */
class Laplace_I_ExactSolution: public ExactSolutionBasicType{
public:

  //! Basic Constructor
  Laplace_I_ExactSolution(void): A(0.0), B(0.0), C(0.0){ 
    ExactSolutionName = "Laplace 1, w(x,y) = A*x + B*y + C";	// Name the exact solution
  }

  //! Return exact solution
  double EvaluateSolutionAt(const double &x, const double &y) const {return A*x + B*y + C;}

  //! Return the exact solution gradient
  Vector2D EvaluateGradientAt(const double &x, const double &y) const {return Vector2D(A,B); }

  //! Calculate the PDE RHS
  double PDE_RightHandSide(const double &x, const double &y) const {return 0.0; }

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:
  double A, B, C;		//!< coefficients of the exact solution
};


/*! 
 * \class Laplace_II_ExactSolution
 * 
 * \brief Implements a particular exact solution to the Laplace equation: 
 *        \f$ \frac{\partial^2 w}{\partial x^2} + \frac{\partial^2 w}{\partial y^2} = 0 \f$
 *
 * The solution has the expression: \f$ w(x,y) = A*(x*x - y*y) + B*x*y \f$
 */
class Laplace_II_ExactSolution: public ExactSolutionBasicType{
public:

  //! Basic Constructor
  Laplace_II_ExactSolution(void): A(0.0), B(0.0) {
    ExactSolutionName = "Laplace 2, w(x,y) = A*(x*x - y*y) + B*x*y";	// Name the exact solution
  }

  //! Return exact solution 
  double EvaluateSolutionAt(const double &x, const double &y) const {return A*(x*x - y*y) + B*x*y;}

  //! Return the exact solution gradient
  Vector2D EvaluateGradientAt(const double &x, const double &y) const {return Vector2D(2.0*A*x + B*y, -2.0*A*y + B*x); }

  //! Calculate the PDE RHS
  double PDE_RightHandSide(const double &x, const double &y) const {return 0.0; }

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:
  double A, B;		//!< coefficients of the exact solution
};


/*! 
 * \class Laplace_III_ExactSolution
 * 
 * \brief Implements a particular exact solution to the Laplace equation: 
 *        \f$ \frac{\partial^2 w}{\partial x^2} + \frac{\partial^2 w}{\partial y^2} = 0 \f$
 *
 * The solution has the expression: \f$ w(x,y) = A*(x*x*x - 3.0*x*y*y) + B*(3*x*x*y - y*y*y) \f$
 */
class Laplace_III_ExactSolution: public ExactSolutionBasicType{
public:

  //! Basic Constructor
  Laplace_III_ExactSolution(void): A(0.0), B(0.0) {
    ExactSolutionName = "Laplace 3, w(x,y) = A*(x*x*x - 3.0*x*y*y) + B*(3*x*x*y - y*y*y)";	// Name the exact solution
  }

  //! Return exact solution 
  double EvaluateSolutionAt(const double &x, const double &y) const {return A*(x*x*x - 3.0*x*y*y) + B*(3*x*x*y - y*y*y);}

  //! Return the exact solution gradient
  Vector2D EvaluateGradientAt(const double &x, const double &y) const {
    return Vector2D(3.0*A*(x*x - y*y) + 6.0*B*x*y , -6.0*A*x*y + 3.0*B*(x*x - y*y) ) ;
 }

  //! Calculate the PDE RHS
  double PDE_RightHandSide(const double &x, const double &y) const {return 0.0; }

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:
  double A, B;		//!< coefficients of the exact solution
};


/*! 
 * \class Laplace_IV_ExactSolution
 * 
 * \brief Implements a particular exact solution to the Laplace equation: 
 *        \f$ \frac{\partial^2 w}{\partial x^2} + \frac{\partial^2 w}{\partial y^2} = 0 \f$
 *
 * The solution has the expression: \f$ w(x,y) = e^{(\mu \, x)}\, (A \, \cos(\mu \, y) + B\, \sin(\mu \, y)) \f$
 */
class Laplace_IV_ExactSolution: public ExactSolutionBasicType{
public:

  //! Basic Constructor
  Laplace_IV_ExactSolution(void): A(0.0), B(0.0), mu(0.0) {
    ExactSolutionName = "Laplace 4, w(x,y) = exp(mu*x)*(A*cos(mu*y) + B*sin(mu*y))";	// Name the exact solution
  }

  //! Return exact solution 
  double EvaluateSolutionAt(const double &x, const double &y) const {return exp(mu*x)*(A*cos(mu*y) + B*sin(mu*y));}

  //! Return the exact solution gradient
  Vector2D EvaluateGradientAt(const double &x, const double &y) const {
    return (mu*exp(mu*x))*Vector2D( A*cos(mu*y) + B*sin(mu*y), -A*sin(mu*y) + B*cos(mu*y));
  }

  //! Calculate the PDE RHS
  double PDE_RightHandSide(const double &x, const double &y) const {return 0.0; }

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:
  double A, B, mu;		//!< coefficients of the exact solution
};


/*! 
 * \class Laplace_V_ExactSolution
 * 
 * \brief Implements a particular exact solution to the Laplace equation: 
 *        \f$ \frac{\partial^2 w}{\partial x^2} + \frac{\partial^2 w}{\partial y^2} = 0 \f$
 *
 * The solution has the expression: 
 *  \f$ w(x,y) = A \, \sinh(\mu \, x) + B \, \cosh(\mu \, x))  \,  (C \, \cos(\mu \, y) + D \, \sin(\mu \, y) \f$
 */
class Laplace_V_ExactSolution: public ExactSolutionBasicType{
public:

  //! Basic Constructor
  Laplace_V_ExactSolution(void): A(0.0), B(0.0), C(0.0), D(0.0), mu(0.0) {
    // Name the exact solution
    ExactSolutionName = "Laplace 5, w(x,y) = A*sinh(mu*x) + B*cosh(mu*x)) * (C*cos(mu*y) + D*sin(mu*y)";
  }

  //! Return exact solution 
  double EvaluateSolutionAt(const double &x, const double &y) const {
    return ( A*sinh(mu*x) + B*cosh(mu*x)) * (C*cos(mu*y) + D*sin(mu*y));
  }

  //! Return the exact solution gradient
  Vector2D EvaluateGradientAt(const double &x, const double &y) const {
    return mu*Vector2D( (A*cosh(mu*x) + B*sinh(mu*x)) * (C*cos(mu*y) + D*sin(mu*y)),
			 (A*sinh(mu*x) + B*cosh(mu*x)) *(-C*sin(mu*y) + D*cos(mu*y)) );
  }

  //! Calculate the PDE RHS
  double PDE_RightHandSide(const double &x, const double &y) const {return 0.0; }

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:
  double A, B, C, D, mu;		//!< coefficients of the exact solution
};


/*! 
 * \class Poisson_I_ExactSolution
 * 
 * \brief Implements a particular exact solution to the non-linear Poisson equation: 
 *        \f$ \frac{\partial^2 w}{\partial x^2} + \frac{\partial^2 w}{\partial y^2} = a\,e^{\beta \, w} \f$
 *
 * The solution has the expression: 
 * \f$ w(x,y) = \frac{1}{\beta} \, \ln \left[ \frac{2(A^2 + B^2)}{ a \beta (Ax + By + C)^2 } \right], \,\, a\beta>0  \f$
 */
class Poisson_I_ExactSolution: public ExactSolutionBasicType{
public:

  //! Basic Constructor
  Poisson_I_ExactSolution(void): A(0.0), B(0.0), C(0.0), a(0.0), beta(0.0){ 
    // Name the exact solution
    ExactSolutionName = "Poisson 1, w(x,y) = log( (2.0*(A*A + B*B) )/( a*beta*(A*x + B*y + C)*(A*x + B*y + C) ) )/beta";
  }

  //! Return exact solution
  double EvaluateSolutionAt(const double &x, const double &y) const {
    double Temp(A*x + B*y + C);
    return log( (2.0*(A*A + B*B) )/( a*beta*Temp*Temp ) )/beta;
  }

  //! Calculate the PDE RHS
  double PDE_RightHandSide(const double &x, const double &y) const {return a*exp( beta * EvaluateSolutionAt(x,y) ); }

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:
  double A, B, C, a, beta;		//!< coefficients of the exact solution
};


/*! 
 * \class Poisson_II_ExactSolution
 * 
 * \brief Implements a particular exact solution to the non-linear Poisson equation: 
 *        \f$ \frac{\partial^2 w}{\partial x^2} + \frac{\partial^2 w}{\partial y^2} = a\,e^{\beta \, w} \f$
 *
 * The solution has the expression: 
 *  \f$ w(x,y) = \frac{1}{\beta} \, \ln \left[ \frac{2(A^2 + B^2)}{ a \beta \sinh^2(Ax + By + C) } \right], \,\, a\beta>0  \f$
 */
class Poisson_II_ExactSolution: public ExactSolutionBasicType{
public:

  //! Basic Constructor
  Poisson_II_ExactSolution(void): A(0.0), B(0.0), C(0.0), a(0.0), beta(0.0) {
    // Name the exact solution
    ExactSolutionName = "Poisson 2, w(x,y) = log( (2.0*(A*A + B*B) )/( a*beta* sinh^2(A*x + B*y + C) ) )/beta";
  }

  //! Return exact solution 
  double EvaluateSolutionAt(const double &x, const double &y) const {
    double Temp(sinh(A*x + B*y + C));
    return log( (2.0*(A*A + B*B) )/( a*beta* Temp* Temp ) )/beta;
  }

  //! Calculate the PDE RHS
  double PDE_RightHandSide(const double &x, const double &y) const {return a*exp( beta * EvaluateSolutionAt(x,y) ); }

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:
  double A, B, C, a, beta;		//!< coefficients of the exact solution
};


/*! 
 * \class Poisson_III_ExactSolution
 * 
 * \brief Implements a particular exact solution to the non-linear Poisson equation: 
 *        \f$ \frac{\partial^2 w}{\partial x^2} + \frac{\partial^2 w}{\partial y^2} = a\,e^{\beta \, w} \f$
 *
 * The solution has the expression:
 * \f$ w(x,y) = \frac{1}{\beta} \, \ln \left[ \frac{-2(A^2 + B^2)}{ a \beta \cosh^2(Ax + By + C) } \right], \,\, a\beta<0  \f$
 */
class Poisson_III_ExactSolution: public ExactSolutionBasicType{
public:

  //! Basic Constructor
  Poisson_III_ExactSolution(void): A(0.0), B(0.0), C(0.0), a(0.0), beta(0.0) {
    // Name the exact solution
    ExactSolutionName = "Poisson 3, w(x,y) = log( (-2.0*(A*A + B*B) )/( a*beta* cosh^2(A*x + B*y + C) ) )/beta";
  }

  //! Return exact solution 
  double EvaluateSolutionAt(const double &x, const double &y) const {
    double Temp(cosh(A*x + B*y + C));
    return log( (-2.0*(A*A + B*B) )/( a*beta* Temp* Temp ) )/beta;
  }

  //! Calculate the PDE RHS
  double PDE_RightHandSide(const double &x, const double &y) const {return a*exp( beta * EvaluateSolutionAt(x,y) ); }

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:
  double A, B, C, a, beta;		//!< coefficients of the exact solution
};


/*! 
 * \class Poisson_IV_ExactSolution
 * 
 * \brief Implements a particular exact solution to the non-linear Poisson equation: 
 *        \f$ \frac{\partial^2 w}{\partial x^2} + \frac{\partial^2 w}{\partial y^2} = a\,e^{\beta \, w} \f$
 *
 * The solution has the expression:
 *   \f$ w(x,y) = \frac{1}{\beta} \, \ln \left[ \frac{2(A^2 + B^2)}{ a \beta \cos^2(Ax + By + C) } \right], \,\, a\beta>0  \f$
 */
class Poisson_IV_ExactSolution: public ExactSolutionBasicType{
public:

  //! Basic Constructor
  Poisson_IV_ExactSolution(void): A(0.0), B(0.0), C(0.0), a(0.0), beta(0.0) {
    // Name the exact solution
    ExactSolutionName = "Poisson 4, w(x,y) = log( (2.0*(A*A + B*B) )/( a*beta* cos^2(A*x + B*y + C) ) )/beta";
  }

  //! Return exact solution 
  double EvaluateSolutionAt(const double &x, const double &y) const {
    double Temp(cos(A*x + B*y + C));
    return log( (2.0*(A*A + B*B) )/( a*beta* Temp* Temp ) )/beta;
  }

  //! Calculate the PDE RHS
  double PDE_RightHandSide(const double &x, const double &y) const {return a*exp( beta * EvaluateSolutionAt(x,y) ); }

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:
  double A, B, C, a, beta;		//!< coefficients of the exact solution
};


/*! 
 * \class Poisson_V_ExactSolution
 * 
 * \brief Implements a particular exact solution to the non-linear Poisson equation: 
 *        \f$ \frac{\partial^2 w}{\partial x^2} + \frac{\partial^2 w}{\partial y^2} = a\,e^{\beta \, w} \f$
 *
 * The solution has the expression:
 *  \f$ w(x,y) = \frac{1}{\beta}\left[\ln \left(\frac{8C}{a\beta} \right) - 2\ln | (x+A)^2 + (y+B)^2 -C | \right] \f$
 */
class Poisson_V_ExactSolution: public ExactSolutionBasicType{
public:

  //! Basic Constructor
  Poisson_V_ExactSolution(void): A(0.0), B(0.0), C(0.0), a(0.0), beta(0.0) {
    // Name the exact solution
    ExactSolutionName = "Poisson 5, w(x,y) = (log(8.0*C/a/beta) - 2.0*log(fabs( (x+A)*(x+A) + (y+B)*(y+B) - C )) )/beta";
  };

  //! Return exact solution 
  double EvaluateSolutionAt(const double &x, const double &y) const {
    return (log(8.0*C/a/beta) - 2.0*log(fabs( (x+A)*(x+A) + (y+B)*(y+B) - C )) )/beta;
  }

  //! Return the exact solution gradient
  Vector2D EvaluateGradientAt(const double &x, const double &y) const {
    double Temp((x+A)*(x+A) + (y+B)*(y+B) - C);
    double SignTemp(sign(Temp));
    return Vector2D( -4.0*(SignTemp*(x+A)) /beta /fabs(Temp),
		     -4.0*(SignTemp*(y+B)) /beta /fabs(Temp) );
  }

  //! Calculate the PDE RHS
  double PDE_RightHandSide(const double &x, const double &y) const {return a*exp( beta * EvaluateSolutionAt(x,y) ); }

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:
  double A, B, C, a, beta;		//!< coefficients of the exact solution
};


/*! 
 * \class StationaryHeatEqnWithLinearSource_ExactSolution
 * 
 * \brief Implements a particular exact solution to the stationary heat equation with linear source: 
 *   \f$ \frac{\partial^2 w}{\partial x^2} + \frac{\partial^2 w}{\partial y^2} = \lambda w \f$
 *
 * The solution is written for a rectangular domain on which the following boundary conditions have been imposed: \n
 *        -> Dirichlet BC on oposite edges (Edge I: w(CoordA) = SolnA; Edge II:  w(CoordB) = SolnB) \n
 *        -> Neumann BC on oposite faces
 */
class StationaryHeatEqnWithLinearSource_ExactSolution: public ExactSolutionBasicType{
public:

  //! Basic Constructor
  StationaryHeatEqnWithLinearSource_ExactSolution(void): lambda(1.0), lambda_sqrt(1.0),
							 CoordA(1.0), SolnA(1.0),
							 CoordB(1.0), SolnB(1.0),
							 C1(1.0), C2(1.0), detA(1.0),
							 Direction_Of_Variation(X_DIRECTION) {
    // Name the exact solution
    ExactSolutionName = "Stationary heat equation with linear source";
  };

  //! Return exact solution 
  double EvaluateSolutionAt(const double &x, const double &y) const;

  //! Calculate the one directional variation
  double Helper_EvaluateSolutionAt(const double &Var1D) const;

  //! Return the exact solution gradient
  Vector2D EvaluateGradientAt(const double &x, const double &y) const ;

  //! Calculate the PDE RHS
  double PDE_RightHandSide(const double &x, const double &y) const {return lambda*EvaluateSolutionAt(x,y); }

  //! Update internal variables
  void Set_ParticularSolution_Parameters(void);

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:
  short Direction_Of_Variation;	        //!< parameter that defines the direction on which Dirichlet BC are imposed (x-axis = 0)
  double lambda, lambda_sqrt;		//!< the coefficient of the source term
  double CoordA, CoordB, SolnA, SolnB;	//!< parameters related to the imposed BCs
  double C1, C2; 	                //!< particular integration constants, defined by the boundary conditions
  double detA;			        //!< internal variable 
};

//! Update internal variables
inline void StationaryHeatEqnWithLinearSource_ExactSolution::
Set_ParticularSolution_Parameters(void){

  if (lambda > ZERO){
    lambda_sqrt = sqrt(lambda);
    detA = cos(lambda_sqrt*CoordA)*sin(lambda_sqrt*CoordB)-cos(lambda_sqrt*CoordB)*sin(lambda_sqrt*CoordA);
    C1 = ( SolnA*sin(lambda_sqrt*CoordB)-SolnB*sin(lambda_sqrt*CoordA))/detA;
    C2 = (-SolnA*cos(lambda_sqrt*CoordB)+SolnB*cos(lambda_sqrt*CoordA))/detA;
  } else if (lambda < ZERO){
    lambda_sqrt = sqrt(-lambda);
    detA = cosh(lambda_sqrt*CoordA)*sinh(lambda_sqrt*CoordB)-cosh(lambda_sqrt*CoordB)*sinh(lambda_sqrt*CoordA);
    C1 = ( SolnA*sinh(lambda_sqrt*CoordB)-SolnB*sinh(lambda_sqrt*CoordA))/detA;
    C2 = (-SolnA*cosh(lambda_sqrt*CoordB)+SolnB*cosh(lambda_sqrt*CoordA))/detA;
  } else {
    detA = CoordB-CoordA;
    C1 = (SolnA*CoordB-SolnB*CoordA)/detA;
    C2 = (-SolnA+SolnB)/detA;
  }

}

//! Return exact solution 
inline double StationaryHeatEqnWithLinearSource_ExactSolution::
EvaluateSolutionAt(const double &x, const double &y) const {

  if (Direction_Of_Variation == X_DIRECTION){
    // the solution is build using x-axis
    return Helper_EvaluateSolutionAt(x);
  } else {
    // the solution is build using y-axis
    return Helper_EvaluateSolutionAt(y);
  }
}

//! Calculate the one directional variation
inline double StationaryHeatEqnWithLinearSource_ExactSolution::
Helper_EvaluateSolutionAt(const double &Var1D) const{

  if (lambda > ZERO) {
    return C1*cos(lambda_sqrt * Var1D) + C2*sin(lambda_sqrt * Var1D);

  } else if (lambda < ZERO) {
    return C1*cosh(lambda_sqrt * Var1D)+C2*sinh(lambda_sqrt * Var1D);

  } else { /* (lambda == ZERO) */ 
    return C1 + C2 * Var1D;
  } /* endif */
}

//! Return the exact solution gradient
inline Vector2D StationaryHeatEqnWithLinearSource_ExactSolution::
EvaluateGradientAt(const double &x, const double &y) const {

  if (Direction_Of_Variation == X_DIRECTION){
    // return the gradient based on the x coordinate
    return Vector2D(lambda_sqrt*(C1*exp(lambda_sqrt*x) - C2*exp(-lambda_sqrt*x)),
		    0.0);
  } else {
    // return the gradient based on the y coordinate
    return Vector2D(0.0,
		    lambda_sqrt*(C1*exp(lambda_sqrt*y) - C2*exp(-lambda_sqrt*y)) ); 
  }
}

/*! 
 * \class PureCircularAdvectionAtConstantSpin_ExactSolution
 * 
 * \brief Implements the exact solution for a circular advection problem:
 *   \f$ \frac{\partial u}{\partial t} + \nabla \cdot (\vec{V} \, u)  = 0 \f$
 * 
 * The solution accepts different inflow conditions and considers them 
 * implemented along the x-axis.
 */
class PureCircularAdvectionAtConstantSpin_ExactSolution: public ExactSolutionBasicType{
public:

  //! Basic Constructor
  PureCircularAdvectionAtConstantSpin_ExactSolution(void): CenterOfRotation(0.0) {
    // Name the exact solution
    ExactSolutionName = "Pure circular advection at constant spin";

    // Get access to the AdvectDiffuse2D_InflowField object
    Inflow = &AdvectDiffuse2D_InflowField::getInstance();
  };

  //! Return exact solution 
  double EvaluateSolutionAt(const double &x, const double &y) const;

  //! Calculate the PDE RHS
  double PDE_RightHandSide(const double &x, const double &y) const {return 0.0; }

  //! Update internal variables
  void Set_ParticularSolution_Parameters(void){ };

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:
  Vector2D CenterOfRotation;	//!< location of the center of rotation
  AdvectDiffuse2D_InflowField *Inflow;		//!< pointer to the inflow field
};

//! Return exact solution 
inline double PureCircularAdvectionAtConstantSpin_ExactSolution::
EvaluateSolutionAt(const double &x, const double &y) const {
  double R;

  // Calculate the distance R relative to the center of rotation
  R = abs(Vector2D(x,y) - CenterOfRotation);

  if (Inflow->IsInflowFieldSet()){
    return Inflow->Solution(CenterOfRotation.x + R, CenterOfRotation.y);
  } else {
    throw runtime_error("PureCircularAdvectionAtConstantSpin_ExactSolution::EvaluateSolutionAt() ERROR! Inflow field not set.");
  }
}

/*! 
 * \class AdvectionDiffusionInRectangularChannel_ExactSolution
 * 
 * \brief Implements the exact solution for an advection-diffusion problem of constant velocity V(uo,0):
 *   \f$ \frac{\partial u}{\partial t} + \nabla \cdot (\vec{V} u) = \nabla \cdot (k \, \nabla u) \f$
 * 
 * The inflow condition is considered to be \f$ \sin(\pi \, y) \f$ and diffusion coefficient constant.
 * The velocity is along the x-axis.
 * The solution is written for a rectangular domain on which the following boundary conditions have been imposed: \n
 * \f$ T(x,0) = T(x,1) = 0.0 \f$ ; \f$ T(0,y) = sin(\pi \, y) \f$ ; \f$ \frac{\partial T(L,y)}{\partial x} = 0 \f$
 *
 * This test case is presented by C. Ollivier-Gooch and M.Van Altena in JCP 181, 729-752 (2002).
 */
class AdvectionDiffusionInRectangularChannel_ExactSolution: public ExactSolutionBasicType{
public:

  //! Basic Constructor
  AdvectionDiffusionInRectangularChannel_ExactSolution(void): XVelocity(1.0),
							      k(0.01),
							      L(3.0){
    // Name the exact solution
    ExactSolutionName = "Advection diffusion in a rectangular channel along x-axis";
  };

  //! Return exact solution 
  double EvaluateSolutionAt(const double &x, const double &y) const;

  //! Calculate the PDE RHS
  double PDE_RightHandSide(const double &x, const double &y) const {return 0.0; }

  //! Update internal variables
  void Set_ParticularSolution_Parameters(void);

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:
  double XVelocity;	//!< velocity in the x-direction
  double k;		//!< diffusion coefficient
  double L;             //!< channel length

  double r1,r2,R,Term;	//!< internal variables
};

//! Return exact solution 
inline double AdvectionDiffusionInRectangularChannel_ExactSolution::
EvaluateSolutionAt(const double &x, const double &y) const {

  if (x <= L){
    /* The expression below is the exact implementation of the analytic solution.
       However, it poses some numerical problems for high ratios between velocity and diffusion coefficient. */
    //  return sin(PI*y)*( (r2*exp(r1*x + r2*L) - r1*exp(r1*L+r2*x)) / (r2*exp(r2*L) - r1*exp(r1*L)) );
    
    // Rearrange the terms to avoid numerical problems
    return sin(PI*y)*( (R*pow(Term*exp(x),r1) - exp(r2*x)) / (R*pow(Term,r1) - 1.0) );
  } else {
    return 0.0;
  }

}

//! Set particular parameters/variables
inline void AdvectionDiffusionInRectangularChannel_ExactSolution::
Set_ParticularSolution_Parameters(void){
  double Temp(HALF*XVelocity/k);

  r1 = Temp + sqrt(sqr(Temp) + sqr(PI));
  r2 = Temp - sqrt(sqr(Temp) + sqr(PI));

  R = r2/r1;

  Term = exp(L*(R-1));
}

/*! 
 * \class SinusoidalVariation_ExactSolution
 * 
 * \brief Implements a sinusoidal wave exact solution. 
 * 
 * This variation can represent an exact solution to an unsteady periodic problem, \n
 * which has as initial condition the variation:
 *   \f$ u(x,y) = \sin(\pi \, \phi(Var)) \f$, where \f$ \phi(Var) \f$ represents the mapping of Var into the domain [-1:1] \n
 * Var is either x or y, depending on the chosen direction. \n
 */
class SinusoidalVariation_ExactSolution: public ExactSolutionBasicType{
public:

  //! Basic Constructor
  SinusoidalVariation_ExactSolution(void): MinCoord(-1.0),
					   MaxCoord( 1.0),
					   Direction(X_DIRECTION){
    // Name the exact solution
    ExactSolutionName = "Sinusoidal wave, w(x,y) = sin(PI*Var), Var in [Min Var-Coordinate,Max Var-Coordinate]";
  };

  //! Return exact solution 
  double EvaluateSolutionAt(const double &x, const double &y) const;

  //! Return the one dimensional solution
  double EvaluateSolutionAt_OneVariable(const double &Var) const;

  //! Calculate the PDE RHS
  double PDE_RightHandSide(const double &x, const double &y) const {return 0.0; }

  //! Update internal variables
  void Set_ParticularSolution_Parameters(void){ };

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:
  double MinCoord;      //!< minimum coordinate
  double MaxCoord;	//!< maximum coordinate
  short Direction;	//!< the direction of interest ('x' or 'y')
};

//! Return exact solution 
inline double SinusoidalVariation_ExactSolution::
EvaluateSolutionAt(const double &x, const double &y) const {
  if (Direction == X_DIRECTION) {
    return EvaluateSolutionAt_OneVariable(x);
  } else {
    return EvaluateSolutionAt_OneVariable(y);
  }
}

//! Return the one dimensional solution
inline double SinusoidalVariation_ExactSolution::
EvaluateSolutionAt_OneVariable(const double &Var) const {
  if (Var<MinCoord || Var>MaxCoord){
    return 0.0;
  } else {
    return sin((ConvertDomainToMinusOneOne(MinCoord,MaxCoord,Var))*PI);
  }
}


/*! 
 * \class ModulatedSinusoidalVariation_ExactSolution
 * 
 * \brief Implements a modulated sinusoidal wave exact solution. 
 * 
 * This variation can represent an exact solution to an unsteady periodic problem, \n
 * which has as initial condition the variation:
 *   \f$ u(x,y) = A \, \cos(B\,\phi(Var)) \sin(C \,\pi \, \phi(Var)) \f$,
 * where \f$ \phi(Var) \f$ represents the mapping of Var into the domain [-1:1] \n
 * Var is either x or y, depending on the chosen direction. \n
 */
class ModulatedSinusoidalVariation_ExactSolution: public ExactSolutionBasicType{
public:

  //! Basic Constructor
  ModulatedSinusoidalVariation_ExactSolution(void): MinCoord(-1.0),
						    MaxCoord( 1.0),
						    Direction(X_DIRECTION),
						    A(1.0), B(3.0), C(5.0){
    // Name the exact solution
    ExactSolutionName = "Modulated sinusoidal wave, w(x,y) = A*cos(B*Var)*sin(C*PI*Var), Var in [Min Var,Max Var]";
  };

  //! Return exact solution 
  double EvaluateSolutionAt(const double &x, const double &y) const;

  //! Return the one dimensional solution
  double EvaluateSolutionAt_OneVariable(const double &Var) const;

  //! Calculate the PDE RHS
  double PDE_RightHandSide(const double &x, const double &y) const {return 0.0; }

  //! Update internal variables
  void Set_ParticularSolution_Parameters(void){ };

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:
  double MinCoord;      //!< minimum coordinate
  double MaxCoord;	//!< maximum coordinate
  short Direction;	//!< the direction of interest ('x' or 'y')
  double A, B, C;       //!< variation parameters
};

//! Return exact solution 
inline double ModulatedSinusoidalVariation_ExactSolution::
EvaluateSolutionAt(const double &x, const double &y) const {
  if (Direction == X_DIRECTION) {
    return EvaluateSolutionAt_OneVariable(x);
  } else {
    return EvaluateSolutionAt_OneVariable(y);
  }
}

//! Return the one dimensional solution
inline double ModulatedSinusoidalVariation_ExactSolution::
EvaluateSolutionAt_OneVariable(const double &Var) const {
  if (Var<MinCoord || Var>MaxCoord){
    return 0.0;
  } else {
    return ( A * cos(B * (ConvertDomainToMinusOneOne(MinCoord,MaxCoord,Var))) *
	     sin(C * (ConvertDomainToMinusOneOne(MinCoord,MaxCoord,Var))*PI) );
  }
}

#endif
