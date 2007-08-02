/**********************************************************************
 * Dusty2DQuad.cc: Subroutines for the 2D dusty quadrilateral         *
 *                 solution block class.                              *
 **********************************************************************/

// Include 2D Dusty quadrilateral solution block class.

#ifndef _DUSTY2D_QUAD_INCLUDED
#include "Dusty2DQuad.h"
#endif // _DUSTY2D_QUAD_INCLUDED

/**********************************************************************
 * Dusty2D_Quad_Block -- Create storage for the static variables.     *
 **********************************************************************/
int Dusty2D_Quad_Block::NUM_VAR_DUSTY2D = NUM_VAR_BASE;
int Dusty2D_Quad_Block::residual_variable = 1;
