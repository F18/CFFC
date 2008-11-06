/*!\file AdvectDiffuse2DAccuracyAssessmentMultiBlock.h
  \brief Header file defining specializations of 2D accuracy assessment class for advection-diffusion solution. */

#ifndef _ADVECTDIFFUSE2D_ACCURACY_ASSESSMENT_MULTIBLOCK_2D_INCLUDED
#define _ADVECTDIFFUSE2D_ACCURACY_ASSESSMENT_MULTIBLOCK_2D_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../HighOrderReconstruction/AccuracyAssessment2DMultiBlock.h" /* Include 2D accuracy assessment header file. */

/*******************************************
 *          SPECIALIZATIONS                *
 ******************************************/

/*!
 * Assess the solution accuracy of the problem based on entropy variation
 * relative to a reference state.
 * This type of error assessment criterion is not valid for advection-diffusion.
 * Just return an error!
 */
template<> inline
int AccuracyAssessment2D_MultiBlock::
AssessSolutionAccuracyBasedOnEntropyVariation(AdvectDiffuse2D_Quad_Block *SolnBlk,
					      const AdaptiveBlock2D_List &Soln_Block_List,
					      const AdvectDiffuse2D_Input_Parameters &IP){

  std::cerr << endl
	    << " ================================================================== " 
	    << endl
	    << " ERROR: The accuracy of the solution couldn't be determined!" << endl
	    << " The entropy variation criterion cannot be used for advection-diffusion! " << endl
	    << " Don't require accuracy assessment with this criterion!" << endl
	    << " ================================================================== " 
	    << endl;  
  return 1;
}

#endif
