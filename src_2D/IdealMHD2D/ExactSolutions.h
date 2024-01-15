/*!\file ExactSolutions.h
  \brief Header file defining 2D MHD exact/reference solutions. */

#ifndef _MHD2D_EXACTSOLUTION_INCLUDED
#define _MHD2D_EXACTSOLUTION_INCLUDED

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
#include "MHD3DState.h"


/*******************************************************
 *                                                     *
 *        MHD2D Exact Solution Type Indexes            *
 *                                                     *
 ******************************************************/
#define MHD2D_NO_EXACT_SOLUTION              -1
// ========== PARTICULAR SOLUTIONS FOR MHD EQUATIONS OR RECONSTRUCTIONS PERFORMED WITH THE MHD SOLVER ===========
#define MHD2D_EXACT_SOLUTION_UNITTEST_FUNCTION 0


// Declare the input parameters class
class MHD2D_Input_Parameters;

/*! 
 * \class ExactSolutionBasicType_MHD2D
 *
 * \brief Basic data type for any exact solution type.
 *
 * This is an abstract data type (ADT).
 */
class ExactSolutionBasicType_MHD2D{
public:

  //! @name Public types
  //@{
  enum Accuracy_Measurement_Type {Soln=1, Grad=2 }; // Solution or gradient
  //@}

  //! Default ctor
  ExactSolutionBasicType_MHD2D(void);

  //! @name Virtual member functions
  //@{
  //! Declare a pure virtual destructor
  virtual ~ExactSolutionBasicType_MHD2D(void) = 0;

  /*! Calculate the exact solution based on the location of interest */
  virtual MHD3D_pState EvaluateSolutionAt(const double &x, const double &y) = 0;

  /*! Calculate the integrated dependency of the solution with respect to x at the location of interest.
    This function can be used to integrate the RHS over domain with curved boundaries,
    by using a Gauss quadrature along the domain contour.
  */
  virtual MHD3D_pState XDependencyIntegrated_Solution(const double &x, const double &y) const {
    throw runtime_error("XDependencyIntegrated_Solution() ERROR! The current exact solution doesn't have the x-dependency integrated function of the solution defined");
  }

  //! Update internal variables
  virtual void Set_ParticularSolution_Parameters(const MHD2D_Input_Parameters & IP){ };

  //! @name Functions for input-output and broadcast
  //@{
  virtual void Parse_Next_Input_Control_Parameter(MHD2D_Input_Parameters & IP, int & i_command) = 0;
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
 * \class UnitTest_Function_ExactSolution_MHD
 * 
 * \brief Implements a 2D function used for testing purposes.
 *        
 */
class UnitTest_Function_ExactSolution_MHD: public ExactSolutionBasicType_MHD2D{
public:

  //! Basic Constructor
  UnitTest_Function_ExactSolution_MHD(void);

  //! Return exact solution
  MHD3D_pState EvaluateSolutionAt(const double &x, const double &y);

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(MHD2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:

};

// Basic Constructor
inline UnitTest_Function_ExactSolution_MHD::UnitTest_Function_ExactSolution_MHD(void) {
  // Name the exact solution
  ExactSolutionName = "UnitTest Function";	
}

//! Return exact solution 
inline MHD3D_pState UnitTest_Function_ExactSolution_MHD::EvaluateSolutionAt(const double &x, const double &y) {
  return MHD3D_pState(MHD3D_pState::MHD3D_W_REF.d() * (1.1 + sin(x)*cos(y)),
		      MHD3D_pState::MHD3D_W_REF.a() * sin(x),
		      MHD3D_pState::MHD3D_W_REF.a() * cos(y),
		      ZERO,
		      sin(x)*cos(y),
		      -cos(x)*sin(y),
		      ZERO,
		      ZERO, ZERO, ZERO,
		      MHD3D_pState::MHD3D_W_REF.p() * exp(x-y)*(5.1 + x*sin(y)) );
}

#endif
