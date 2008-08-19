/*!\file ExactSolutions.h
  \brief Header file defining 2D Euler exact solutions. */

#ifndef _EULER2D_EXACTSOLUTION_INCLUDED
#define _EULER2D_EXACTSOLUTION_INCLUDED

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
#include "Euler2DState.h"


/*******************************************************
 *                                                     *
 *      Euler2D Exact Solution Type Indexes            *
 *                                                     *
 ******************************************************/
#define EULER2D_NO_EXACT_SOLUTION              -1
// ========== PARTICULAR SOLUTIONS FOR ... EQUATION ===========
#define EULER2D_EXACT_SOLUTION_RINGLEB_FLOW     0



// Declare the input parameters class
class Euler2D_Input_Parameters;

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
  virtual Euler2D_pState EvaluateSolutionAt(const double &x, const double &y) const = 0;

  /*! Calculate the integrated dependency of the solution with respect to x at the location of interest.
    This function can be used to integrate the RHS over domain with curved boundaries,
    by using a Gauss quadrature along the domain contour.
  */
  virtual Euler2D_pState XDependencyIntegrated_Solution(const double &x, const double &y) const {
    throw runtime_error("XDependencyIntegrated_Solution() ERROR! The current exact solution doesn't have the x-dependency integrated function of the solution defined");
  }

  //! Update internal variables
  virtual void Set_ParticularSolution_Parameters(void){ };

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
 * \class Ringleb_Flow_ExactSolution
 * 
 * \brief Implements a particular exact solution to the Euler equations: 
 *
 */
class Ringleb_Flow_ExactSolution: public ExactSolutionBasicType{
public:

  //! Basic Constructor
  Ringleb_Flow_ExactSolution(void){ 
    ExactSolutionName = "Ringleb Flow";	// Name the exact solution
  }

  //! Return exact solution
  Euler2D_pState EvaluateSolutionAt(const double &x, const double &y) const {return Euler2D_pState(0);}

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(Euler2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);

  //! Broadcast relevant parameters
  void Broadcast(void);

private:

};

#endif
