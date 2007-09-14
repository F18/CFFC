/* Cell2D.cc:  Subroutines and data for 2D cells. */

/* Include 2D cell header file. */

#ifndef _CELL2D_INCLUDED
#include "Cell2D.h"
#endif // _CELL2D_INCLUDED

/********************************************************************
 * Cell2D_Cartesian -- Create storage and set cell lengths.         *
 ********************************************************************/
Vector2D Cell2D_Cartesian::dx = Vector2D_ZERO;

/*************************************************************************
 * Cell2D_Quad -- Static constant data memebers.                         *
 *************************************************************************/
const double Cell2D_Quad::Gauss2QuadPoints_CoeffPoint1 = 0.2113248654051871177454256; /* 0.5*(1 - 1/sqrt(3)) */
const double Cell2D_Quad::Gauss2QuadPoints_CoeffPoint2 = 0.7886751345948128822545744; /* 0.5*(1 + 1/sqrt(3)) */
