/* Cell3D.cc:  Subroutines and data for 3D cells. */

/* Include 3D cell header file. */

#ifndef _CELL3D_INCLUDED
#include "Cell3D.h"
#endif // _CELL3D_INCLUDED

/********************************************************************
 * Cell3D_Cartesian -- Create storage and set cell lengths.         *
 ********************************************************************/
Vector3D Cell3D_Hexa::dx = Vector3D_ZERO;
Vector3D Cell3D_Hexa::dy = Vector3D_ZERO;
Vector3D Cell3D_Hexa::dz = Vector3D_ZERO;

/*************************************************************************
 * Cell3D_Quad -- Static constant data memebers.                         *
 *************************************************************************/
const double Cell3D_Hexa::Gauss2QuadPoints_CoeffPoint1 = 0.2113248654051871177454256; /* 0.5*(1 - 1/sqrt(3)) */
const double Cell3D_Hexa::Gauss2QuadPoints_CoeffPoint2 = 0.7886751345948128822545744; /* 0.5*(1 + 1/sqrt(3)) */
