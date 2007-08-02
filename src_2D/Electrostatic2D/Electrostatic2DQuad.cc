/**********************************************************************
 * Electrostatic2DQuad.cc: Subroutines for the 2D electrostatic       *
 *                         quadrilateral solution block class.        *
 **********************************************************************/

// Include 2D Electrostatic quadrilateral solution block class.

#ifndef _ELECTROSTATIC2D_QUAD_INCLUDED
#include "Electrostatic2DQuad.h"
#endif // _ELECTROSTATIC2D_QUAD_INCLUDED

/**********************************************************************
 * Electrostatic2D_Quad_Block -- Create storage for the static        *
 *                               variables.                           *
 **********************************************************************/
int Electrostatic2D_Quad_Block::residual_variable = 1;
int Electrostatic2D_Quad_Block::Flow_Type = FLOWTYPE_LAMINAR;
