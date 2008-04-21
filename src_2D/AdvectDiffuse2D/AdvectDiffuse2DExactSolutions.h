/*!\file AdvectDiffuse2DExactSolutions.h
  \brief Header file defining the singleton advection-diffusion exact solution class. */

#ifndef _ADVECTDIFFUSE2D_EXACTSOLUTION_SINGLETON_INCLUDED
#define _ADVECTDIFFUSE2D_EXACTSOLUTION_SINGLETON_INCLUDED

/* Include required C++ libraries. */
#include <cmath>
#include <iostream>

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Utilities/Utilities.h"
#include "ExactSolutions.h"

/*! 
 * \class AdvectDiffuse2D_ExactSolutions
 * 
 * \brief Collection of exact solutions for Advection-Diffusion 2D problems
 *
 * The class generates only one object (Meyers singleton pattern)
 * which can be access through getInstance() function.
 *
 * \nosubgrouping
 */
class AdvectDiffuse2D_ExactSolutions{
public:
 
  //! @name Access, deallocate and setting functions
  //@{ 
  static AdvectDiffuse2D_ExactSolutions& getInstance(void);
  void SetExactSolution(const short &SolutionIndex);
  static void DestroyExactSolutionObject(void);
  //@}

  //! @name Evaluation of exact solution
  //@{ 
  double operator()(const double &x, const double &y) const;
  double Solution(const double &x, const double &y) const;
  Vector2D Gradient(const double &x, const double &y) const;
  double PDE_RightHandSide(const double &x, const double &y) const;
  double XDependencyIntegrated_PDE_RightHandSide(const double &x, const double &y) const;
  //@}

  //! Update the internal variables of the exact solution
  void Set_ParticularSolution_Parameters(void) const { 
    if ( IsExactSolutionSet() ){ ExactSoln->Set_ParticularSolution_Parameters();}
  }

  //! @name Functions for input-output and broadcast
  //@{
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);
  void Print_Info(std::ostream & out_file);
  void Broadcast(void);
  //@}

  //! Indicate whether the exact solution has been set or not
  bool IsExactSolutionSet(void) const { return (ExactSoln != NULL)? true: false; }

protected:
  AdvectDiffuse2D_ExactSolutions();		//!< Private constructor
  AdvectDiffuse2D_ExactSolutions(const AdvectDiffuse2D_ExactSolutions&); //!< Private copy constructor
  AdvectDiffuse2D_ExactSolutions& operator=(const AdvectDiffuse2D_ExactSolutions&); //!< Private assignment operator
  
private:
  static ExactSolutionBasicType *ExactSoln; //!< pointer to the exact solution
  static short i_Exact_Solution_Type;       //!< type indicator for the exact solution
};

/*! 
 * Evaluate the exact solution 
 * \param x,y Cartesian coordinates
 */
inline double AdvectDiffuse2D_ExactSolutions::operator()(const double &x, const double &y) const{
  return ExactSoln->EvaluateSolutionAt(x,y);
}

/*! 
 * Evaluate the exact solution 
 * \param x,y Cartesian coordinates
 */
inline double AdvectDiffuse2D_ExactSolutions::Solution(const double &x, const double &y) const{
  return ExactSoln->EvaluateSolutionAt(x,y);
}

/*! 
 * Evaluate the exact solution gradient
 * \param x,y Cartesian coordinates
 */
inline Vector2D AdvectDiffuse2D_ExactSolutions::Gradient(const double &x, const double &y) const{
  return ExactSoln->EvaluateGradientAt(x,y);
}

/*! 
 * Evaluate the right-hand-side term of
 * the equation to which the provided
 * function is a solution.
 * \param x,y Cartesian coordinates
 */
inline double AdvectDiffuse2D_ExactSolutions::PDE_RightHandSide(const double &x, const double &y) const{
  return ExactSoln->PDE_RightHandSide(x,y);
}

/*! 
 * Evaluate the integral of right-hand-side term 
 * with respect to x-coordinate of the equation 
 * to which the provided function is a solution.
 * This function can be use to calculate the 
 * integral of the right-hand-side term over 
 * domains with curved boundaries.
 * \param x,y Cartesian coordinates
 */
inline double AdvectDiffuse2D_ExactSolutions::XDependencyIntegrated_PDE_RightHandSide(const double &x,
										      const double &y) const{
  return ExactSoln->XDependencyIntegrated_PDE_RightHandSide(x,y);
}
#endif
