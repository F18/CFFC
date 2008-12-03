/*!\file ExactSolutions.h
  \brief Header file defining 2D NavierStokes exact solutions. */

#ifndef _NAVIERSTOKES2D_EXACTSOLUTION_INCLUDED
#define _NAVIERSTOKES2D_EXACTSOLUTION_INCLUDED

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
#include "NavierStokes2DState.h"


/*******************************************************
 *                                                     *
 *    NavierStokes2D Exact Solution Type Indexes       *
 *                                                     *
 ******************************************************/
#define NAVIERSTOKES2D_NO_EXACT_SOLUTION              -1
// === PARTICULAR SOLUTIONS FOR NAVIER-STOKES EQUATIONS OR RECONSTRUCTIONS PERFORMED WITH THE NAVIERSTOKES SOLVER ===
#define NAVIERSTOKES2D_EXACT_SOLUTION_ABGRALL_FUNCTION 0
#define NAVIERSTOKES2D_EXACT_SOLUTION_SINUSOIDAL_FUNCTION 1
#define NAVIERSTOKES2D_EXACT_SOLUTION_COSSIN_FUNCTION 2
#define NAVIERSTOKES2D_EXACT_SOLUTION_UNITTEST_FUNCTION 3


// Declare the input parameters class
class NavierStokes2D_Input_Parameters;

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
  virtual NavierStokes2D_pState EvaluateSolutionAt(const double &x, const double &y) = 0;

  /*! Calculate the integrated dependency of the solution with respect to x at the location of interest.
    This function can be used to integrate the RHS over domain with curved boundaries,
    by using a Gauss quadrature along the domain contour.
  */
  virtual NavierStokes2D_pState XDependencyIntegrated_Solution(const double &x, const double &y) const {
    throw runtime_error("XDependencyIntegrated_Solution() ERROR! The current exact solution doesn't have the x-dependency integrated function of the solution defined");
  }

  //! Update internal variables
  virtual void Set_ParticularSolution_Parameters(const NavierStokes2D_Input_Parameters & IP){ };

  //! @name Functions for input-output and broadcast
  //@{
  virtual void Parse_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters & IP, int & i_command) = 0;
  virtual void Print_Info(std::ostream & out_file) = 0;
  virtual void Broadcast(void) = 0;
  //@}

  //! @name Tecplot output functions
  //@{
  /*! @brief Output the variable names that defined the exact solution in a format 
   * suitable for Tecplot to the provided output stream.
   * Overwrite the standard format if customized output is desired. */
  virtual void Output_Tecplot_Title(std::ostream & out_file) const;
  /*! @brief Write the Tecplot double precision format for the associated exact solution. */
  virtual void Output_Tecplot_Double_Precision(std::ostream & out_file) const;
  /*! @brief Write the exact solution at the given location in a format suitable for plotting with Tecplot. */
  virtual void Output_Tecplot_Solution(std::ostream & out_file,
				       const double &x, const double &y) const {
    // Overwrite me if you want special format and set IsSpecialTecplotFormatDefined() to true!
  };
  /*!@brief Indicate if the current exact solution is output in a special format.
   * For special formats, overwrite this function and return TRUE.
   */
  virtual bool IsSpecialTecplotFormatDefined(void) const { return false; }
  //@}

  //@}

  const string & whatExactSolution(void) const { return ExactSolutionName; }   //!< Get the name of the exact solution type

protected:
  string ExactSolutionName;		//!< storage for the name of the source field type
  Accuracy_Measurement_Type Accuracy_Parameter;	//!< indicator for the parameter used to estimate the accuracy
                                                // (e.g solution or gradient)

};


/*! 
 * \class Abgrall_Function_ExactSolution
 * 
 * \brief Implements a 2D non-smooth function used for testing the robustness of the reconstruction procedure;
 *
 */
class Abgrall_Function_ExactSolution: public ExactSolutionBasicType{
public:

  //! Basic Constructor
  Abgrall_Function_ExactSolution(void){ 
    ExactSolutionName = "Abgrall Function";	// Name the exact solution
  }

  //! Return exact solution
  NavierStokes2D_pState EvaluateSolutionAt(const double &x, const double &y) ;

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:

  double Abgrall_Function(const double & x, const double & y) const;
  double Basic_1D_Variation(const double & r) const;
  
};

//! Return exact solution 
inline NavierStokes2D_pState Abgrall_Function_ExactSolution::EvaluateSolutionAt(const double &x, const double &y) {
  return NavierStokes2D_pState(Abgrall_Function(x,y), ZERO, ZERO, PRESSURE_STDATM);
}

/*!
 * Calculate the value of the function for a given location.
 */
inline double Abgrall_Function_ExactSolution::Abgrall_Function(const double & x, const double & y) const {

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
inline double Abgrall_Function_ExactSolution::Basic_1D_Variation(const double & r) const {

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
 * \class Sinusoidal_Function_ExactSolution
 * 
 * \brief Implements a 2D sinusoidal function, \f$ \sin[(x+1)\pi] \f$, in either x- or y-direction;
 *
 */
class Sinusoidal_Function_ExactSolution: public ExactSolutionBasicType{
public:

  //! Basic Constructor
  Sinusoidal_Function_ExactSolution(void);

  //! Return exact solution
  NavierStokes2D_pState EvaluateSolutionAt(const double &x, const double &y) ;

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters & IP, int & i_command);

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
inline Sinusoidal_Function_ExactSolution::Sinusoidal_Function_ExactSolution(void):
  DomainMinLimit(-100.0), DomainMaxLimit(100.0),
  ReferenceValue(2.0), Velocity(100.0), Direction(X_DIRECTION)
{
  // Name the exact solution
  ExactSolutionName = "Sinusoidal Function, rho(x,y) = RefVal + sin[(Var+1)PI] for Var in [Min Var, Max Var], otherwise 0";
}

//! Return exact solution 
inline NavierStokes2D_pState Sinusoidal_Function_ExactSolution::EvaluateSolutionAt(const double &x, const double &y) {
  if (Direction == X_DIRECTION){ // X-direction variation
    return NavierStokes2D_pState(BaseSinusoidalVariation(x), Velocity, ZERO, PRESSURE_STDATM);
  } else {			// Y-direction variation
    return NavierStokes2D_pState(BaseSinusoidalVariation(y), ZERO, Velocity, PRESSURE_STDATM);
  }
}

//! Return the basic sinusoidal variation
inline double Sinusoidal_Function_ExactSolution::BaseSinusoidalVariation(const double & Coord) const{
  if (Coord < DomainMinLimit || Coord > DomainMaxLimit) {
    return ReferenceValue;
  } else {
    return ReferenceValue + sin((ConvertDomainToMinusOneOne(DomainMinLimit,DomainMaxLimit,Coord)+1.0)*PI);
  }
}


/*! 
 * \class CosSin_Function_ExactSolution
 * 
 * \brief Implements a 2D trigonometric function, 
 *        \f$ 1.0 + 0.5 \, \cos(\pi x) \sin(5 \pi x) \f$,
 *        in either x- or y-direction;
 *
 */
class CosSin_Function_ExactSolution: public ExactSolutionBasicType{
public:

  //! Basic Constructor
  CosSin_Function_ExactSolution(void);

  //! Return exact solution
  NavierStokes2D_pState EvaluateSolutionAt(const double &x, const double &y) ;

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters & IP, int & i_command);

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
inline CosSin_Function_ExactSolution::CosSin_Function_ExactSolution(void):
  DomainMinLimit(-100.0), DomainMaxLimit(100.0),
  ReferenceValue(1.0), Amplitude(0.5), Velocity(100.0),
  Direction(X_DIRECTION)
{
  // Name the exact solution
  ExactSolutionName = "CosSin Function, rho(x,y) = RefVal + A * cos(PI*Var)*sin(5*PI*Var) for Var in [Min Var, Max Var], otherwise 0";	
}

//! Return exact solution 
inline NavierStokes2D_pState CosSin_Function_ExactSolution::EvaluateSolutionAt(const double &x, const double &y) {
  if (Direction == X_DIRECTION){ // X-direction variation
    return NavierStokes2D_pState(BaseCosSinVariation(x), Velocity, ZERO, PRESSURE_STDATM);
  } else {			// Y-direction variation
    return NavierStokes2D_pState(BaseCosSinVariation(y), ZERO, Velocity, PRESSURE_STDATM);
  }
}

//! Return the basic sinusoidal variation
inline double CosSin_Function_ExactSolution::BaseCosSinVariation(const double & Coord) const {
  if (Coord < DomainMinLimit || Coord > DomainMaxLimit) {
    return ReferenceValue;
  } else {
    return ReferenceValue + ( Amplitude * cos(PI * ConvertDomainToMinusOneOne(DomainMinLimit,DomainMaxLimit,Coord))*
			      sin(5.0 * PI * ConvertDomainToMinusOneOne(DomainMinLimit,DomainMaxLimit,Coord)) );
  }
}


/*! 
 * \class UnitTest_Function_ExactSolution
 * 
 * \brief Implements a 2D function used for testing purposes.
 *        
 */
class UnitTest_Function_ExactSolution: public ExactSolutionBasicType{
public:

  //! Basic Constructor
  UnitTest_Function_ExactSolution(void);

  //! Return exact solution
  NavierStokes2D_pState EvaluateSolutionAt(const double &x, const double &y);

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:

};

// Basic Constructor
inline UnitTest_Function_ExactSolution::UnitTest_Function_ExactSolution(void) {
  // Name the exact solution
  ExactSolutionName = "UnitTest Function";	
}

//! Return exact solution 
inline NavierStokes2D_pState UnitTest_Function_ExactSolution::EvaluateSolutionAt(const double &x, const double &y) {
  NavierStokes2D_pState NavierStokes2D_W_STDATM;
  NavierStokes2D_W_STDATM.Standard_Atmosphere();
  return NavierStokes2D_pState(NavierStokes2D_W_STDATM.rho * (1.1 + sin(x)*cos(y)),
			       NavierStokes2D_W_STDATM.a() * sin(x),
			       NavierStokes2D_W_STDATM.a() * cos(y),
			       NavierStokes2D_W_STDATM.p   * exp(x-y)*(5.1 + x*sin(y)) );
}

#endif
