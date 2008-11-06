/*!\file AdvectDiffuse2DAccuracyAssessment.h
  \brief Header file defining specializations of 2D accuracy assessment class for advection-diffusion solution. */

#ifndef _ADVECTDIFFUSE2D_ACCURACY_ASSESSMENT_2D_INCLUDED
#define _ADVECTDIFFUSE2D_ACCURACY_ASSESSMENT_2D_INCLUDED

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
 * Compute solution entropy errors.
 * This type of error assessment is not valid for advection-diffusion.
 * Just throw an error!
 */
template<>
template<> inline
double AccuracyAssessment2D<AdvectDiffuse2D_Quad_Block>::
ComputeSolutionEntropyErrorsHighOrder(const unsigned int &parameter,
				      const unsigned int &accuracy_digits,
				      const AdvectDiffuse2D_Input_Parameters & IP,
				      const unsigned short int &Pos)
  throw(ArgumentNullException)
{

  // not applicable to advection-diffusion
  throw ArgumentNullException("AccuracyAssessment2D::ComputeSolutionEntropyErrorsHighOrder() ERROR! This type of error assessment doesn't apply to advection-diffusion!");
}

/*!
 * Compute solution entropy errors.
 * This type of error assessment is not valid for advection-diffusion.
 * Just throw an error!
 */
template<>
template<> inline
double AccuracyAssessment2D<AdvectDiffuse2D_Quad_Block>::
ComputeSolutionEntropyErrors(const unsigned int &parameter,
			     const unsigned int &accuracy_digits,
			     const AdvectDiffuse2D_Input_Parameters & IP,
			     double (AdvectDiffuse2D_Quad_Block::*ComputeSolutionEntropyAt)(const int &, const int &,
											    const Vector2D &,
											    const unsigned int &) const)
  throw(ArgumentNullException)
{

  // not applicable to advection-diffusion
  throw ArgumentNullException("AccuracyAssessment2D::ComputeSolutionEntropyErrors() ERROR! This type of error assessment doesn't apply to advection-diffusion!");
}


#endif
