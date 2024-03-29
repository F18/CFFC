/*!\file ExactSolutions.h
  \brief Header file defining 2D Euler exact solutions. */

#ifndef _EULER2D_EXACTSOLUTION_INCLUDED
#define _EULER2D_EXACTSOLUTION_INCLUDED

/* Include required C++ libraries. */
#include <cmath>
#include <iostream>

/* Using std namespace functions */
using std::sin;
using std::fabs;
using std::cos;
using std::tan;
using std::sqrt;
using std::tanh;

/* Include CFFC header files */
#include "../Math/Math.h"
#include "../Math/Vector2D.h"
#include "../Math/NumericalLibrary.h"
#include "../Utilities/TypeDefinition.h"
#include "../Utilities/Utilities.h"
#include "../CFD/CFD.h"
#include "Euler2DState.h"


/*******************************************************
 *                                                     *
 *      Euler2D Exact Solution Type Indexes            *
 *                                                     *
 ******************************************************/
#define EULER2D_NO_EXACT_SOLUTION              -1
// ========== PARTICULAR SOLUTIONS FOR EULER EQUATIONS OR RECONSTRUCTIONS PERFORMED WITH THE EULER SOLVER ===========
#define EULER2D_EXACT_SOLUTION_RINGLEB_FLOW     0
#define EULER2D_EXACT_SOLUTION_ABGRALL_FUNCTION 1
#define EULER2D_EXACT_SOLUTION_SINUSOIDAL_FUNCTION 2
#define EULER2D_EXACT_SOLUTION_COSSIN_FUNCTION 3
#define EULER2D_EXACT_SOLUTION_HYPERTANGENT_FUNCTION 4
#define EULER2D_EXACT_SOLUTION_UNITTEST_FUNCTION 5


// Declare the input parameters class
class Euler2D_Input_Parameters;

/*! 
 * \class ExactSolutionBasicType_Euler2D
 *
 * \brief Basic data type for any exact solution type.
 *
 * This is an abstract data type (ADT).
 */
class ExactSolutionBasicType_Euler2D{
public:

  //! @name Public types
  //@{
  enum Accuracy_Measurement_Type {Soln=1, Grad=2 }; // Solution or gradient
  //@}

  //! Default ctor
  ExactSolutionBasicType_Euler2D(void);

  //! @name Virtual member functions
  //@{
  //! Declare a pure virtual destructor
  virtual ~ExactSolutionBasicType_Euler2D(void) = 0;

  /*! Calculate the exact solution based on the location of interest */
  virtual Euler2D_pState EvaluateSolutionAt(const double &x, const double &y) = 0;

  /*! Calculate the integrated dependency of the solution with respect to x at the location of interest.
    This function can be used to integrate the RHS over domain with curved boundaries,
    by using a Gauss quadrature along the domain contour.
  */
  virtual Euler2D_pState XDependencyIntegrated_Solution(const double &x, const double &y) const {
    throw runtime_error("XDependencyIntegrated_Solution() ERROR! The current exact solution doesn't have the x-dependency integrated function of the solution defined");
  }

  //! Update internal variables
  virtual void Set_ParticularSolution_Parameters(const Euler2D_Input_Parameters & IP){ };

  //! @name Functions for input-output and broadcast
  //@{
  virtual void Parse_Next_Input_Control_Parameter(Euler2D_Input_Parameters & IP, int & i_command) = 0;
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
 * \class Ringleb_Flow_ExactSolution_Euler
 * 
 * \brief Implements a particular exact solution to the Euler equations: 
 *
 */
class Ringleb_Flow_ExactSolution_Euler: public ExactSolutionBasicType_Euler2D{
public:

  //! Basic Constructor
  Ringleb_Flow_ExactSolution_Euler(void);

  //! Return exact solution
  Euler2D_pState EvaluateSolutionAt(const double &x, const double &y) ;

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(Euler2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

  /*! @brief Function for determining c inside Ringleb's flow at location xLoc, yLoc. */
  double RinglebResidual(const double &c) const ;

private:
  double c_a, c_b;
  double g, po, rhoo;
  int Precision, MaxIter;
  double TwoOver_gm1;
  double xLoc, yLoc;

  //! Determine the Cartesian coordinates (X_Coord,Y_Coord) corresponding to the pair of 
  //  non-dimensional velocity q and the streamline parameter k.
  void RinglebFlowCoordinates(const double &q , const double &k , double &X_Coord, double &Y_Coord);

  double J_Func(const double &c) const;
  double rho_Func(const double &c) const;
  double q_Func(const double &c) const;
  
};

/*
 * Default constructor
 */
inline Ringleb_Flow_ExactSolution_Euler::Ringleb_Flow_ExactSolution_Euler(void):
  g(GAMMA_AIR), po(PRESSURE_STDATM), rhoo(DENSITY_STDATM),
  c_a(0.7), c_b(0.999999),
  Precision(14), MaxIter(200),
  TwoOver_gm1(TWO/(g-ONE)),
  xLoc(0.0), yLoc(0.0)
{ 
  ExactSolutionName = "Ringleb Flow";	// Name the exact solution
}


inline double Ringleb_Flow_ExactSolution_Euler::J_Func(const double &c) const{
  return ONE/c + ONE/(THREE*c*c*c) + ONE/(FIVE*c*c*c*c*c) - HALF*log((ONE+c)/(ONE-c));
}

inline double Ringleb_Flow_ExactSolution_Euler::rho_Func(const double &c) const{
  return pow(c,TwoOver_gm1);
}

inline double Ringleb_Flow_ExactSolution_Euler::q_Func(const double &c) const{
  return sqrt(TwoOver_gm1*(ONE - c*c));
}


/*
 * This subroutine determines X_Coord and Y_Coord based on
 * the non-dimensional velocity q and the streamline parameter k.
 */
inline void Ringleb_Flow_ExactSolution_Euler::RinglebFlowCoordinates(const double &q , const double &k,
							       double &X_Coord, double &Y_Coord) {
 
  double c, rho, J;
  double g = 1.40;

  c = sqrt(ONE - ((g-ONE)/TWO)*q*q);
  rho = rho_Func(c);
  J = J_Func(c);

  X_Coord = (HALF/rho)*(TWO/(k*k) - ONE/(q*q)) - HALF*J;
  Y_Coord = (ONE/(k*rho*q))*sqrt(ONE - (q*q)/(k*k));  

}

/*! Function for determining c inside Ringleb's flow at location xLoc, yLoc.
 *  Returns residual of (x+J/2)^2 + y^2 - 1/(4*rho^2*q^4) based on a guess of c.
 *  Parameters xLoc and yLoc must be preset in the class.
 */
inline double Ringleb_Flow_ExactSolution_Euler::RinglebResidual(const double &c) const {

  const double J = J_Func(c);
  const double rho = rho_Func(c);
  const double q = q_Func(c);

  return (xLoc + HALF*J)*(xLoc + 0.5*J) + yLoc*yLoc - 0.25/(rho*rho *q*q*q*q);
}


/*! 
 * \class Abgrall_Function_ExactSolution_Euler
 * 
 * \brief Implements a 2D non-smooth function used for testing the robustness of the reconstruction procedure;
 *
 */
class Abgrall_Function_ExactSolution_Euler: public ExactSolutionBasicType_Euler2D{
public:

  //! Basic Constructor
  Abgrall_Function_ExactSolution_Euler(void){ 
    ExactSolutionName = "Abgrall Function";	// Name the exact solution
  }

  //! Return exact solution
  Euler2D_pState EvaluateSolutionAt(const double &x, const double &y) ;

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(Euler2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:

  double Abgrall_Function(const double & x, const double & y) const;
  double Basic_1D_Variation(const double & r) const;
  
};

//! Return exact solution 
inline Euler2D_pState Abgrall_Function_ExactSolution_Euler::EvaluateSolutionAt(const double &x, const double &y) {
  return Euler2D_pState(Abgrall_Function(x,y), ZERO, ZERO, PRESSURE_STDATM);
}

/*!
 * Calculate the value of the function for a given location.
 */
inline double Abgrall_Function_ExactSolution_Euler::Abgrall_Function(const double & x, const double & y) const {

  if ( x <= 0.5*cos(PI * y) ){
    return (2.0 + Basic_1D_Variation(x - (1.0/tan(sqrt(0.5*PI))) * y) );

  } else if ( x > 0.5*cos(PI * y) ) {
    return (2.0 + Basic_1D_Variation(x + (1.0/tan(sqrt(0.5*PI))) * y) + cos(2.0 * PI * y));

  } else {
    return 0.0e0;
  }

}

/*!
 * Implement the basic 1D variation used by Abgrall's function.
 */
inline double Abgrall_Function_ExactSolution_Euler::Basic_1D_Variation(const double & r) const {

  if (r <= -1.0/3.0){
    return (-r * sin(1.5 * PI * r * r));
    
  } else if (fabs(r) < 1.0/3.0){
    return (fabs(sin(2.0 * PI * r)));

  } else if (r >= 1.0/3.0){
    return (2.0*r - 1.0 + sin(3.0 * PI * r)/6.0);

  } else {
    return 0.0;
  }

}


/*! 
 * \class Sinusoidal_Function_ExactSolution_Euler
 * 
 * \brief Implements a 2D sinusoidal function, \f$ \sin[(x+1)\pi] \f$, in either x- or y-direction;
 *
 */
class Sinusoidal_Function_ExactSolution_Euler: public ExactSolutionBasicType_Euler2D{
public:

  //! Basic Constructor
  Sinusoidal_Function_ExactSolution_Euler(void);

  //! Return exact solution
  Euler2D_pState EvaluateSolutionAt(const double &x, const double &y) ;

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(Euler2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:

  double DomainMinLimit,        //!< minimum limit of the coordinate
    DomainMaxLimit;		//!< maximum limit of the coordinate
  double ReferenceValue;	//!< magnitude of the function
  double Velocity;		//!< convection velocity
  short Direction;		//!< the direction of interest ('x' or 'y')
  
  double BaseSinusoidalVariation(const double & Coord) const;
};

// Basic Constructor
inline Sinusoidal_Function_ExactSolution_Euler::Sinusoidal_Function_ExactSolution_Euler(void):
  DomainMinLimit(-100.0), DomainMaxLimit(100.0),
  ReferenceValue(2.0), Velocity(100.0), Direction(X_DIRECTION)
{
  // Name the exact solution
  ExactSolutionName = "Sinusoidal Function, rho(x,y) = RefVal + sin[(Var+1)PI] for Var in [Min Var, Max Var], otherwise 0";
}

//! Return exact solution 
inline Euler2D_pState Sinusoidal_Function_ExactSolution_Euler::EvaluateSolutionAt(const double &x, const double &y) {
  if (Direction == X_DIRECTION){ // X-direction variation
    return Euler2D_pState(BaseSinusoidalVariation(x), Velocity, ZERO, PRESSURE_STDATM);
  } else {			// Y-direction variation
    return Euler2D_pState(BaseSinusoidalVariation(y), ZERO, Velocity, PRESSURE_STDATM);
  }
}

//! Return the basic sinusoidal variation
inline double Sinusoidal_Function_ExactSolution_Euler::BaseSinusoidalVariation(const double & Coord) const{
  if (Coord < DomainMinLimit || Coord > DomainMaxLimit) {
    return ReferenceValue;
  } else {
    return ReferenceValue + sin((ConvertDomainToMinusOneOne(DomainMinLimit,DomainMaxLimit,Coord)+1.0)*PI);
  }
}


/*! 
 * \class CosSin_Function_ExactSolution_Euler
 * 
 * \brief Implements a 2D trigonometric function, 
 *        \f$ 1.0 + 0.5 \, \cos(\pi x) \sin(5 \pi x) \f$,
 *        in either x- or y-direction;
 *
 */
class CosSin_Function_ExactSolution_Euler: public ExactSolutionBasicType_Euler2D{
public:

  //! Basic Constructor
  CosSin_Function_ExactSolution_Euler(void);

  //! Return exact solution
  Euler2D_pState EvaluateSolutionAt(const double &x, const double &y) ;

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(Euler2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:

  double DomainMinLimit,        //!< Sets the minimum limit of the variable domain
    DomainMaxLimit;		//!< Sets the maximum limit of the variable domain
  double ReferenceValue;	//!< Sets the magnitude of the function
  double Amplitude;		//!< Sets the amplitude of the function
  double Velocity;		//!< Sets the velocity in the direction of variation
  int Direction;		//!< Indicates the direction of variation (0 for X-dir, 1 for Y-dir)
  
  double BaseCosSinVariation(const double & Coord) const;
};

// Basic Constructor
inline CosSin_Function_ExactSolution_Euler::CosSin_Function_ExactSolution_Euler(void):
  DomainMinLimit(-100.0), DomainMaxLimit(100.0),
  ReferenceValue(1.0), Amplitude(0.5), Velocity(100.0),
  Direction(X_DIRECTION)
{
  // Name the exact solution
  ExactSolutionName = "CosSin Function, rho(x,y) = RefVal + A * cos(PI*Var)*sin(5*PI*Var) for Var in [Min Var, Max Var], otherwise 0";	
}

//! Return exact solution 
inline Euler2D_pState CosSin_Function_ExactSolution_Euler::EvaluateSolutionAt(const double &x, const double &y) {
  if (Direction == X_DIRECTION){		// X-direction variation
    return Euler2D_pState(BaseCosSinVariation(x), Velocity, ZERO, PRESSURE_STDATM);
  } else {			// Y-direction variation
    return Euler2D_pState(BaseCosSinVariation(y), ZERO, Velocity, PRESSURE_STDATM);
  }
}

//! Return the basic sinusoidal variation
inline double CosSin_Function_ExactSolution_Euler::BaseCosSinVariation(const double & Coord) const {
  if (Coord < DomainMinLimit || Coord > DomainMaxLimit) {
    return ReferenceValue;
  } else {
    return ReferenceValue + ( Amplitude * cos(PI * ConvertDomainToMinusOneOne(DomainMinLimit,DomainMaxLimit,Coord))*
			      sin(5.0 * PI * ConvertDomainToMinusOneOne(DomainMinLimit,DomainMaxLimit,Coord)) );
  }
}


/*! 
 * \class HyperTangent_Function_ExactSolution_Euler
 * 
 * \brief Implements the 2D trigonometric function 
 *        \f$ rho(Val) = RefVal + A (1.0 - \tanh^2 (S Val)) \f$,
 *        where Val is the distance between the calculation point
 *        and the centroid of the function.
 */
class HyperTangent_Function_ExactSolution_Euler: public ExactSolutionBasicType_Euler2D{
public:

  //! Basic Constructor
  HyperTangent_Function_ExactSolution_Euler(void);

  //! Return exact solution
  Euler2D_pState EvaluateSolutionAt(const double &x, const double &y) ;

  //! Update the translation point location based on the flow input parameters
  void Set_ParticularSolution_Parameters(const Euler2D_Input_Parameters & IP);

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(Euler2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:

  double ReferenceValue;	//!< Sets the magnitude of the function
  double Amplitude,		//!< Sets the amplitude of the function
    Steepness;			//!< Sets the steepness of the function
  Vector2D Centroid;		//!< Sets the reference point for distance calculation
  bool SetCentroid;		//!< Indicates if the centroid position is the initial (false) or the final one (true)
  
  double BaseHyperTangentVariation(const double & DistanceToCentroid) const;
};

// Basic Constructor
inline HyperTangent_Function_ExactSolution_Euler::HyperTangent_Function_ExactSolution_Euler(void):
  ReferenceValue(1.0), Amplitude(1.0), Steepness(1.0),
  Centroid(0.0,0.0), SetCentroid(false)
{
  // Name the exact solution
  ExactSolutionName = "HyperTangent Function";	
}

//! Return exact solution 
inline Euler2D_pState HyperTangent_Function_ExactSolution_Euler::EvaluateSolutionAt(const double &x, const double &y) {
  // Compute distance between centroid and the current position
  double Distance(abs(Centroid - Vector2D(x,y)));

  return Euler2D_pState(BaseHyperTangentVariation(Distance), ZERO, ZERO, PRESSURE_STDATM);
}

//! Return the basic sinusoidal variation
inline double HyperTangent_Function_ExactSolution_Euler::BaseHyperTangentVariation(const double & DistanceToCentroid) const {

  // Function: "1 - tanh(x)^2" modulated for amplitude and steepness
  double Func;

  Func = tanh(Steepness * DistanceToCentroid);
  Func *= Func; 		// square hyperbolic tangent
  Func = 1.0 - Func;

  return ReferenceValue + Amplitude*Func;
}


/*! 
 * \class UnitTest_Function_ExactSolution_Euler
 * 
 * \brief Implements a 2D function used for testing purposes.
 *        
 */
class UnitTest_Function_ExactSolution_Euler: public ExactSolutionBasicType_Euler2D{
public:

  //! Basic Constructor
  UnitTest_Function_ExactSolution_Euler(void);

  //! Return exact solution
  Euler2D_pState EvaluateSolutionAt(const double &x, const double &y);

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(Euler2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:

};

// Basic Constructor
inline UnitTest_Function_ExactSolution_Euler::UnitTest_Function_ExactSolution_Euler(void) {
  // Name the exact solution
  ExactSolutionName = "UnitTest Function";	
}

//! Return exact solution 
inline Euler2D_pState UnitTest_Function_ExactSolution_Euler::EvaluateSolutionAt(const double &x, const double &y) {
  return Euler2D_pState(Euler2D_W_STDATM.d   * (1.1 + sin(x)*cos(y)),
			Euler2D_W_STDATM.a() * sin(x),
			Euler2D_W_STDATM.a() * cos(y),
			Euler2D_W_STDATM.p   * exp(x-y)*(5.1 + x*sin(y)) );
}

#endif
