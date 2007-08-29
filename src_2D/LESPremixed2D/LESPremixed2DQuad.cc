/**********************************************************************
 * LESPremixed2DQuad.cc: Subroutines for the 2D LESPremixed           *
 *                       quadrilateral solution  block class.         *
 **********************************************************************/

// Include 2D LESPremixed quadrilateral solution block class.

#ifndef _LESPREMIXED2D_QUAD_INCLUDED
#include "LESPremixed2DQuad.h"
#endif // _LESPREMIXED2D_QUAD_INCLUDED

/**********************************************************************
 * LESPremixed2D_Quad_Block -- Create storage for the total number of *
 *                             variables.                             *
 **********************************************************************/
int LESPremixed2D_Quad_Block::residual_variable = 1;
int LESPremixed2D_Quad_Block::Number_of_Residual_Norms = 4;
//int LESPremixed_Quad_Block::Flow_Type = FLOWTYPE_INVISCID;
