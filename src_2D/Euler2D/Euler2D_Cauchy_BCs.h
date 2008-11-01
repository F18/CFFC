/*!\file Euler2D_Cauchy_BCs.h
  \brief Header file defining specializations of Cauchy_BCs template class for Euler2D solution states. */

#ifndef _EULER2D_CAUCHY_BCS_INCLUDED
#define _EULER2D_CAUCHY_BCS_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "Euler2DState.h"					  /* Include Euler2D solution states header file */
#include "../HighOrderReconstruction/Cauchy_BoundaryConditions.h" /* Include Cauchy_BCs header file */

//////////////////////////////////////////////////////////
// Specialized routines for Euler2D solution states     //
//////////////////////////////////////////////////////////

/*!
 * Specialization of SetCharacteristicConstraintsBasedOnBCtype()
 * from Cauchy_BCs class.
 */
template<> inline
void Cauchy_BCs<Euler2D_pState>::SetCharacteristicConstraintsBasedOnBCtype(const int & BCtype){
  
  switch(BCtype){

  case BC_REFLECTION:		// treat it in the same way as BC_WALL_INVISCID
  case BC_WALL_INVISCID:
    // Impose relational constraint between 'u' and 'v'
    RelationalConstraints[0] = 0; // rho
    RelationalConstraints[1] = 1; // u
    RelationalConstraints[2] = 1; // v
    RelationalConstraints[3] = 0; // p
    RelationalConstraintsFlag = true;
    // Impose No individual constraints
    IndividualConstraints[0] = 0;
    IndividualConstraints[1] = 0;
    IndividualConstraints[2] = 0;
    IndividualConstraints[3] = 0;
    IndividualConstraintsFlag = false;
    break;

  case BC_EXACT_SOLUTION:
    // Impose NO relational constraints
    RelationalConstraints[0] = 0;
    RelationalConstraints[1] = 0;
    RelationalConstraints[2] = 0;
    RelationalConstraints[3] = 0;
    RelationalConstraintsFlag = false;
    // Assume individual constraints for each parameter.
    // (the final implementation might be different depending on the flow conditions)
    IndividualConstraints[0] = 1;
    IndividualConstraints[1] = 1;
    IndividualConstraints[2] = 1;
    IndividualConstraints[3] = 1;
    IndividualConstraintsFlag = true;
    break;

  default:
    // Impose NO relational constraints
    RelationalConstraints[0] = 0;
    RelationalConstraints[1] = 0;
    RelationalConstraints[2] = 0;
    RelationalConstraints[3] = 0;
    RelationalConstraintsFlag = false;
    // Impose No individual constraints
    IndividualConstraints[0] = 0;
    IndividualConstraints[1] = 0;
    IndividualConstraints[2] = 0;
    IndividualConstraints[3] = 0;
    IndividualConstraintsFlag = false;
  } // endswitch
}

#endif
