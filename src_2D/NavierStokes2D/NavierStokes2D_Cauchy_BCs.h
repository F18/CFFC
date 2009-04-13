/*!\file NavierStokes2D_Cauchy_BCs.h
  \brief Header file defining specializations of Cauchy_BCs template class for NavierStokes2D solution states. */

#ifndef _NAVIERSTOKES2D_CAUCHY_BCS_INCLUDED
#define _NAVIERSTOKES2D_CAUCHY_BCS_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "NavierStokes2DState.h"				  /* Include NavierStokes2D solution states header file */
#include "../HighOrderReconstruction/Cauchy_BoundaryConditions.h" /* Include Cauchy_BCs header file */

/////////////////////////////////////////////////////////////////
// Specialized routines for NavierStokes2D solution states     //
////////////////////////////////////////////////////////////////

/*!
 * Specialization of SetCharacteristicConstraintsBasedOnBCtype()
 * from Cauchy_BCs class.
 */
template<> inline
void Cauchy_BCs<NavierStokes2D_pState>::SetCharacteristicConstraintsBasedOnBCtype(const int & BCtype){

  // Set default constraints
  // Impose NO relational constraints
  RelationalConstraints[0] = 0; // rho
  RelationalConstraints[1] = 0;	// u 
  RelationalConstraints[2] = 0;	// v 
  RelationalConstraints[3] = 0;	// p 
  RelationalConstraints[4] = 0;	// k 
  RelationalConstraints[5] = 0;	// omega
  RelationalConstraints[6] = 0;	// ke
  RelationalConstraints[7] = 0;	// ee
  RelationalConstraintsFlag = false;
  // Impose No individual constraints
  IndividualConstraints[0] = 0; // rho
  IndividualConstraints[1] = 0;	// u 
  IndividualConstraints[2] = 0;	// v 
  IndividualConstraints[3] = 0;	// p 
  IndividualConstraints[4] = 0;	// k 
  IndividualConstraints[5] = 0;	// omega
  IndividualConstraints[6] = 0;	// ke
  IndividualConstraints[7] = 0;	// ee
  IndividualConstraintsFlag = false;


  // Set particular constraints based on the boundary condition type
  switch(BCtype){

  case BC_REFLECTION:		// treat it in the same way as BC_WALL_INVISCID
  case BC_WALL_INVISCID:
    // Impose relational constraint between 'u' and 'v'
    RelationalConstraints[1] = 1; // u
    RelationalConstraints[2] = 1; // v
    RelationalConstraintsFlag = true;
    break;

  case BC_EXACT_SOLUTION:
    // Assume individual constraints for each parameter.
    // (the final implementation might be different depending on the flow conditions)
    IndividualConstraints[0] = 1;
    IndividualConstraints[1] = 1;
    IndividualConstraints[2] = 1;
    IndividualConstraints[3] = 1;
    IndividualConstraints[4] = 1;
    IndividualConstraints[5] = 1;
    IndividualConstraints[6] = 1;
    IndividualConstraints[7] = 1;
    IndividualConstraintsFlag = true;
    break;

  case BC_FARFIELD:
    // Assume individual constraints for velocity and pressure.
    // (the final implementation might be different depending on the flow conditions)
    IndividualConstraints[0] = 0;
    IndividualConstraints[1] = 1;
    IndividualConstraints[2] = 1;
    IndividualConstraints[3] = 1;
    IndividualConstraintsFlag = true;
    break;

  case BC_CONSTANT_EXTRAPOLATION:
    // Assume zero gradient for all variables.
    // (the final implementation might be different depending on the flow conditions)
    IndividualConstraints[0] = 1;
    IndividualConstraints[1] = 1;
    IndividualConstraints[2] = 1;
    IndividualConstraints[3] = 1;
    IndividualConstraints[4] = 1;
    IndividualConstraints[5] = 1;
    IndividualConstraints[6] = 1;
    IndividualConstraints[7] = 1;
    IndividualConstraintsFlag = true;
    break;

  case BC_DIRICHLET:		// same as BC_FIXED
    // Assume individual constraints for each parameter.
    // (the final implementation might be different depending on the flow conditions)
    IndividualConstraints[0] = 1;
    IndividualConstraints[1] = 1;
    IndividualConstraints[2] = 1;
    IndividualConstraints[3] = 1;
    IndividualConstraints[4] = 1;
    IndividualConstraints[5] = 1;
    IndividualConstraints[6] = 1;
    IndividualConstraints[7] = 1;
    IndividualConstraintsFlag = true;
    break;

  case BC_WALL_VISCOUS_HEATFLUX:
    // Individual constraints for both velocities (i.e. Dirichlet bc) as well as density and pressure (i.e. Neumann bc)
    // The current setup is for laminar flows!
    IndividualConstraints[0] = 1;
    IndividualConstraints[1] = 1;
    IndividualConstraints[2] = 1;
    IndividualConstraints[3] = 1;
    IndividualConstraintsFlag = true;
    break;

  case BC_OUTFLOW_SUBSONIC :
    // Impose individual constraints (i.e. Dirichlet and Neumann) to all variables.
    // The current setup is for laminar flows!
    IndividualConstraints[0] = 1;
    IndividualConstraints[1] = 1;
    IndividualConstraints[2] = 1;
    IndividualConstraints[3] = 1;
    IndividualConstraintsFlag = true;

  case BC_INFLOW_SUBSONIC :
    // Impose individual constraints (i.e. Dirichlet and Neumann) to all variables.
    // The current setup is for laminar flows!
    IndividualConstraints[0] = 1;
    IndividualConstraints[1] = 1;
    IndividualConstraints[2] = 1;
    IndividualConstraints[3] = 1;
    IndividualConstraintsFlag = true;

  default:
    // Leave the default values unchanged
    break;

  } // endswitch

}

#endif
