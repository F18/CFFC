/*!\file MHD2DExactSolutions.h
  \brief Header file defining the singleton MHD2D exact solution class. */

#ifndef _MHD2D_EXACTSOLUTION_SINGLETON_INCLUDED
#define _MHD2D_EXACTSOLUTION_SINGLETON_INCLUDED

/* Include required C++ libraries. */
#include <cmath>
#include <iostream>

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Utilities/Utilities.h"
#include "ExactSolutions.h"

/*! 
 * \class MHD2D_ExactSolutions
 * 
 * \brief Collection of exact/reference solutions for MHD 2D problems
 *
 * The class generates only one object (Meyers singleton pattern)
 * which can be access through getInstance() function.
 *
 * \nosubgrouping
 */
class MHD2D_ExactSolutions{
public:
 
  //! @name Access, deallocate and setting functions
  //@{ 
  static MHD2D_ExactSolutions& getInstance(void);
  void SetExactSolution(const short &SolutionIndex);
  static void DestroyExactSolutionObject(void);
  //@}

  //! @name Evaluation of exact solution
  //@{ 
  MHD3D_pState operator()(const double &x, const double &y) const;
  MHD3D_pState Solution(const double &x, const double &y) const;
  double SolutionForParameter(const double &x, const double &y, const unsigned &parameter) const;
  MHD3D_pState XDependencyIntegrated_Solution(const double &x, const double &y) const;
  double XDependencyIntegrated_SolutionForParameter(const double &x, const double &y, const unsigned &parameter) const;
  //@}

  //! Update the internal variables of the exact solution
  void Set_ParticularSolution_Parameters(MHD2D_Input_Parameters & IP) const { 
    if ( IsExactSolutionSet() ){ ExactSoln->Set_ParticularSolution_Parameters(IP);}
  }

  //! @name Functions for input-output and broadcast
  //@{
  void Parse_Next_Input_Control_Parameter(MHD2D_Input_Parameters & IP, int & i_command);
  void Print_Info(std::ostream & out_file);
  void Broadcast(void);
  //@}

  //! Indicate whether the exact solution has been set or not
  bool IsExactSolutionSet(void) const { return (ExactSoln != NULL)? true: false; }

protected:
  MHD2D_ExactSolutions();		//!< Private constructor
  MHD2D_ExactSolutions(const MHD2D_ExactSolutions&); //!< Private copy constructor
  MHD2D_ExactSolutions& operator=(const MHD2D_ExactSolutions&); //!< Private assignment operator
  
private:
  static ExactSolutionBasicType_MHD2D *ExactSoln; //!< pointer to the exact solution
  static short i_Exact_Solution_Type;       //!< type indicator for the exact solution
};

/*! 
 * Evaluate the exact solution 
 * \param x,y Cartesian coordinates
 */
inline MHD3D_pState MHD2D_ExactSolutions::operator()(const double &x, const double &y) const{
  return ExactSoln->EvaluateSolutionAt(x,y);
}

/*! 
 * Evaluate the exact solution 
 * \param x,y Cartesian coordinates
 */
inline MHD3D_pState MHD2D_ExactSolutions::Solution(const double &x, const double &y) const{
  return ExactSoln->EvaluateSolutionAt(x,y);
}

/*! 
 * Evaluate the exact solution for a specified parameter.
 * \param x,y Cartesian coordinates
 * \param parameter the state variable for which the exact solution is evaluated
 */
inline double MHD2D_ExactSolutions::SolutionForParameter(const double &x, const double &y,
							   const unsigned &parameter) const{
  return ExactSoln->EvaluateSolutionAt(x,y)[parameter];
}

/*! 
 * Evaluate the integral of the solution 
 * with respect to x-coordinate.
 * This function can be use to calculate the 
 * integral of the solution over domains 
 * with curved boundaries.
 * \param x,y Cartesian coordinates
 */
inline MHD3D_pState MHD2D_ExactSolutions::XDependencyIntegrated_Solution(const double &x, const double &y) const{
  return ExactSoln->XDependencyIntegrated_Solution(x,y);
}

/*! 
 * Evaluate the integrated X dependency of the exact solution for a specified parameter.
 * \param x,y Cartesian coordinates
 * \param parameter the state variable for which the integrated x dependency is evaluated
 */
inline double MHD2D_ExactSolutions::XDependencyIntegrated_SolutionForParameter(const double &x, const double &y,
										 const unsigned &parameter) const {
  return ExactSoln->XDependencyIntegrated_Solution(x,y)[parameter];
}
#endif
