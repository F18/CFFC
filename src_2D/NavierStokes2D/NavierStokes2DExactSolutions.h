/*!\file NavierStokes2DExactSolutions.h
  \brief Header file defining the singleton NavierStokes2D exact solution class. */

#ifndef _NAVIERSTOKES2D_EXACTSOLUTION_SINGLETON_INCLUDED
#define _NAVIERSTOKES2D_EXACTSOLUTION_SINGLETON_INCLUDED

/* Include required C++ libraries. */
#include <cmath>
#include <iostream>

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Utilities/Utilities.h"
#include "ExactSolutions.h"

/*! 
 * \class NavierStokes2D_ExactSolutions
 * 
 * \brief Collection of exact solutions for Navier-Stokes 2D problems
 *
 * The class generates only one object (Meyers singleton pattern)
 * which can be access through getInstance() function.
 *
 * \nosubgrouping
 */
class NavierStokes2D_ExactSolutions{
public:
 
  //! @name Access, deallocate and setting functions
  //@{ 
  static NavierStokes2D_ExactSolutions& getInstance(void);
  void SetExactSolution(const short &SolutionIndex);
  static void DestroyExactSolutionObject(void);
  //@}

  //! @name Evaluation of exact solution
  //@{ 
  NavierStokes2D_pState operator()(const double &x, const double &y) const;
  NavierStokes2D_pState Solution(const double &x, const double &y) const;
  double SolutionForParameter(const double &x, const double &y, const unsigned &parameter) const;
  NavierStokes2D_pState XDependencyIntegrated_Solution(const double &x, const double &y) const;
  double XDependencyIntegrated_SolutionForParameter(const double &x, const double &y, const unsigned &parameter) const;
  //@}

  //! Update the internal variables of the exact solution
  void Set_ParticularSolution_Parameters(NavierStokes2D_Input_Parameters & IP) const { 
    if ( IsExactSolutionSet() ){ ExactSoln->Set_ParticularSolution_Parameters(IP);}
  }

  //! @name Functions for input-output and broadcast
  //@{
  void Parse_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters & IP, int & i_command);
  void Print_Info(std::ostream & out_file);
  void Broadcast(void);
  //@}

  //! Indicate whether the exact solution has been set or not
  bool IsExactSolutionSet(void) const { return (ExactSoln != NULL)? true: false; }

  //! @name Tecplot output functions
  //@{
  void Output_Tecplot_Title(std::ostream & out_file) const;
  void Output_Tecplot_Double_Precision(std::ostream & out_file) const;
  void Output_Tecplot_Solution(std::ostream & out_file,
			       const double &x, const double &y) const;
  //@}

protected:
  NavierStokes2D_ExactSolutions();		//!< Private constructor
  NavierStokes2D_ExactSolutions(const NavierStokes2D_ExactSolutions&); //!< Private copy constructor
  NavierStokes2D_ExactSolutions& operator=(const NavierStokes2D_ExactSolutions&); //!< Private assignment operator
  
private:
  static ExactSolutionBasicType_NavierStokes2D *ExactSoln; //!< pointer to the exact solution
  static short i_Exact_Solution_Type;       //!< type indicator for the exact solution
};

/*! 
 * Evaluate the exact solution 
 * \param x,y Cartesian coordinates
 */
inline NavierStokes2D_pState NavierStokes2D_ExactSolutions::operator()(const double &x, const double &y) const{
  return ExactSoln->EvaluateSolutionAt(x,y);
}

/*! 
 * Evaluate the exact solution 
 * \param x,y Cartesian coordinates
 */
inline NavierStokes2D_pState NavierStokes2D_ExactSolutions::Solution(const double &x, const double &y) const{
  return ExactSoln->EvaluateSolutionAt(x,y);
}

/*! 
 * Evaluate the exact solution for a specified parameter.
 * \param x,y Cartesian coordinates
 * \param parameter the state variable for which the exact solution is evaluated
 */
inline double NavierStokes2D_ExactSolutions::SolutionForParameter(const double &x, const double &y,
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
inline NavierStokes2D_pState NavierStokes2D_ExactSolutions::XDependencyIntegrated_Solution(const double &x,
											   const double &y) const{
  return ExactSoln->XDependencyIntegrated_Solution(x,y);
}

/*! 
 * Evaluate the integrated X dependency of the exact solution for a specified parameter.
 * \param x,y Cartesian coordinates
 * \param parameter the state variable for which the integrated x dependency is evaluated
 */
inline double NavierStokes2D_ExactSolutions::XDependencyIntegrated_SolutionForParameter(const double &x, const double &y,
										 const unsigned &parameter) const {
  return ExactSoln->XDependencyIntegrated_Solution(x,y)[parameter];
}

/*! 
 * Output the title of the associated solution to a Tecplot file.
 */
inline void NavierStokes2D_ExactSolutions::Output_Tecplot_Title(std::ostream & out_file) const {
  return ExactSoln->Output_Tecplot_Title(out_file);
}

/*! 
 * Output the double precision format for the associated solution to a Tecplot file.
 */
inline void NavierStokes2D_ExactSolutions::Output_Tecplot_Double_Precision(std::ostream & out_file) const {
  return ExactSoln->Output_Tecplot_Double_Precision(out_file);
}

/*! Write the exact solution in a format suitable for plotting with Tecplot.
 * \param out_file output stream
 * \param x the x-coordinate of the location of interest
 * \param y the y-coordinate of the location of interest
 */
inline void NavierStokes2D_ExactSolutions::Output_Tecplot_Solution(std::ostream & out_file,
								   const double &x, const double &y) const{

  NavierStokes2D_pState Soln;

  if (ExactSoln->IsSpecialTecplotFormatDefined()){
    // redirect the output to the particular solution
    ExactSoln->Output_Tecplot_Solution(out_file,
				       x, y);
  } else {
    // use a generic output
    Soln = Solution(x,y);

    // Output the four basic primitive variables (i.e. rho, u, v, p)
    out_file << " " << Soln.rho 
	     << " " << Soln.v
	     << " " << Soln.p;
  } // endif
}

#endif
