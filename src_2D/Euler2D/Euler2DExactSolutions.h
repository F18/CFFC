/*!\file Euler2DExactSolutions.h
  \brief Header file defining the singleton Euler2D exact solution class. */

#ifndef _EULER2D_EXACTSOLUTION_SINGLETON_INCLUDED
#define _EULER2D_EXACTSOLUTION_SINGLETON_INCLUDED

/* Include required C++ libraries. */
#include <cmath>
#include <iostream>

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Utilities/Utilities.h"
#include "ExactSolutions.h"

/*! 
 * \class Euler2D_ExactSolutions
 * 
 * \brief Collection of exact solutions for Euler 2D problems
 *
 * The class generates only one object (Meyers singleton pattern)
 * which can be access through getInstance() function.
 *
 * \nosubgrouping
 */
class Euler2D_ExactSolutions{
public:
 
  //! @name Access, deallocate and setting functions
  //@{ 
  static Euler2D_ExactSolutions& getInstance(void);
  void SetExactSolution(const short &SolutionIndex);
  static void DestroyExactSolutionObject(void);
  //@}

  //! @name Evaluation of exact solution
  //@{ 
  Euler2D_pState operator()(const double &x, const double &y) const;
  Euler2D_pState Solution(const double &x, const double &y) const;
  Euler2D_pState XDependencyIntegrated_Solution(const double &x, const double &y) const;
  //@}

  //! Update the internal variables of the exact solution
  void Set_ParticularSolution_Parameters(void) const { 
    if ( IsExactSolutionSet() ){ ExactSoln->Set_ParticularSolution_Parameters();}
  }

  //! @name Functions for input-output and broadcast
  //@{
  void Parse_Next_Input_Control_Parameter(Euler2D_Input_Parameters & IP, int & i_command);
  void Print_Info(std::ostream & out_file);
  void Broadcast(void);
  //@}

  //! Indicate whether the exact solution has been set or not
  bool IsExactSolutionSet(void) const { return (ExactSoln != NULL)? true: false; }

protected:
  Euler2D_ExactSolutions();		//!< Private constructor
  Euler2D_ExactSolutions(const Euler2D_ExactSolutions&); //!< Private copy constructor
  Euler2D_ExactSolutions& operator=(const Euler2D_ExactSolutions&); //!< Private assignment operator
  
private:
  static ExactSolutionBasicType *ExactSoln; //!< pointer to the exact solution
  static short i_Exact_Solution_Type;       //!< type indicator for the exact solution
};

/*! 
 * Evaluate the exact solution 
 * \param x,y Cartesian coordinates
 */
inline Euler2D_pState Euler2D_ExactSolutions::operator()(const double &x, const double &y) const{
  return ExactSoln->EvaluateSolutionAt(x,y);
}

/*! 
 * Evaluate the exact solution 
 * \param x,y Cartesian coordinates
 */
inline Euler2D_pState Euler2D_ExactSolutions::Solution(const double &x, const double &y) const{
  return ExactSoln->EvaluateSolutionAt(x,y);
}

/*! 
 * Evaluate the integral of the solution 
 * with respect to x-coordinate.
 * This function can be use to calculate the 
 * integral of the solution over domains 
 * with curved boundaries.
 * \param x,y Cartesian coordinates
 */
inline Euler2D_pState Euler2D_ExactSolutions::XDependencyIntegrated_Solution(const double &x, const double &y) const{
  return ExactSoln->XDependencyIntegrated_Solution(x,y);
}

#endif