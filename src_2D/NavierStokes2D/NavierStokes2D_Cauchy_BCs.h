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

}

#endif
