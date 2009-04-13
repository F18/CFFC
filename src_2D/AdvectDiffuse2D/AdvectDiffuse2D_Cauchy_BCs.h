/*!\file AdvectDiffuse2D_Cauchy_BCs.h
  \brief Header file defining specializations of Cauchy_BCs template class for AdvectDiffuse2D solution states. */

#ifndef _EULER2D_CAUCHY_BCS_INCLUDED
#define _EULER2D_CAUCHY_BCS_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "AdvectDiffuse2DState.h"				  /* Include AdvectDiffuse2D solution states header file */
#include "../HighOrderReconstruction/Cauchy_BoundaryConditions.h" /* Include Cauchy_BCs header file */

///////////////////////////////////////////////////////////////
// Specialized routines for AdvectDiffuse2D solution states  //
///////////////////////////////////////////////////////////////

/*!
 * Specialization of SetCharacteristicConstraintsBasedOnBCtype()
 * from Cauchy_BCs class.
 */
template<> inline
void Cauchy_BCs<AdvectDiffuse2D_State>::SetCharacteristicConstraintsBasedOnBCtype(const int & BCtype){
  
  switch(BCtype){
    // No special case yet
  default:
    // Impose No relational constraints
    RelationalConstraints[0] = 0;
    RelationalConstraintsFlag = false;
    // Impose Individual constraints
    IndividualConstraints[0] = 1;
    IndividualConstraintsFlag = true;
  } // endswitch
}

#endif
