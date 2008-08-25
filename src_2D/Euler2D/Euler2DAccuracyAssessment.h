/*!\file Euler2DAccuracyAssessment.h
  \brief Header file defining specializations of 2D accuracy assessment class for Euler solution. */

#ifndef _EULER2D_ACCURACY_ASSESSMENT_2D_INCLUDED
#define _EULER2D_ACCURACY_ASSESSMENT_2D_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../HighOrderReconstruction/AccuracyAssessment2D.h" /* Include 2D accuracy assessment header file. */

/*******************************************
 *          SPECIALIZATIONS                *
 ******************************************/


/*!
 * Compute solution errors relative to the exact solution
 * set in Euler2D_Quad_Block.
 *
 * \param parameter the state variable which is used for computing the errors
 * \accuracy_digits the number of exact digits with which the error is sought
 * \param AccessToHighOrderVars Euler2D_Quad_Block member function which returns the 
 *                              array of high-order objects in the block.
 * \param Pos  The specific high-order object which is used to compute the error.
 * \throw ArgumentNullException when the exact solution pointer is NULL
 */
template<> inline
void AccuracyAssessment2D<Euler2D_Quad_Block>::
ComputeSolutionErrorsHighOrder(const unsigned int &parameter,
			       const unsigned int &accuracy_digits,
			       const unsigned short int &Pos)
  throw(ArgumentNullException)
{

  // Set the type of the returned value
  double _dummy_param(0.0);
  
  if (SolnBlk->ExactSolution()->IsExactSolutionSet()){  // exact solution is set

    // Compute the errors
    SolnBlk->HighOrderVariable(Pos).
      ComputeSolutionErrors(wrapped_member_function_one_parameter(SolnBlk->ExactSolution(),
								  &Euler2D_Quad_Block::Exact_Solution_Type::SolutionForParameter,
								  parameter,
								  _dummy_param),
			    parameter,
			    accuracy_digits);
    
    // Write the final information in the designated variables
    L1() = SolnBlk->HighOrderVariable(Pos).L1();
    L2() = SolnBlk->HighOrderVariable(Pos).L2();
    LMax() = SolnBlk->HighOrderVariable(Pos).LMax();
    TotalBlockArea = SolnBlk->HighOrderVariable(Pos).BlockArea();
    CellsUsed = SolnBlk->HighOrderVariable(Pos).UsedCells();
    
  } else {
    // exact solution is not set
    throw ArgumentNullException("AccuracyAssessment2D::ComputeSolutionErrors() ERROR! There is no exact solution set!");
  }  
}

/*!
 * Compute solution errors relative to the exact solution
 * set in Euler2D_Quad_Block.
 *
 * \param parameter the state variable which is used for computing the errors
 * \param ComputeSolutionAt Euler2D_Quad_Block member function which returns the solution 
 *                          at a given location (PointOfInterest) using the reconstruction 
 *                          of cell (i,j) for a specified parameter.
 * \throw ArgumentNullException when the exact solution pointer is NULL
 *
 */
template<> inline
void AccuracyAssessment2D<Euler2D_Quad_Block>::
ComputeSolutionErrors(const unsigned int &parameter,
		      const unsigned int &accuracy_digits,
		      double (Euler2D_Quad_Block::*ComputeSolutionAt)(const int &, const int &,
								      const Vector2D &,
								      const unsigned int &) const )
  throw(ArgumentNullException)
{

  // Set the type of the returned value
  double _dummy_param(0.0);
  
  if (SolnBlk->ExactSolution()->IsExactSolutionSet()){  // exact solution is set
    ComputeSolutionErrors(wrapped_member_function_one_parameter(SolnBlk->ExactSolution(),
								&Euler2D_Quad_Block::Exact_Solution_Type::SolutionForParameter,
								parameter,
								_dummy_param),
			  parameter,
			  accuracy_digits,
			  ComputeSolutionAt);
  } else {
    // exact solution is not set
    throw ArgumentNullException("AccuracyAssessment2D::ComputeSolutionErrors() ERROR! There is no exact solution set!");
  }
  
}

#endif
