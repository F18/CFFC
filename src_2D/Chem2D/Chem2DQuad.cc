/**********************************************************************
 * Chem2DQuad.cc: Subroutines for the 2D chem quadrilateral solution  *
 *                block class.                                        *
 **********************************************************************/

// Include 2D Chem quadrilateral solution block class.

#ifndef _CHEM2D_QUAD_INCLUDED
#include "Chem2DQuad.h"
#endif // _CHEM2D_QUAD_INCLUDED

/**********************************************************************
 * Chem2D_Quad_Block -- Create storage for the total number of        *
 *                      variables.                                    *
 **********************************************************************/
int Chem2D_Quad_Block::residual_variable = 1;
int Chem2D_Quad_Block::Number_of_Residual_Norms = 4;
//int Chem2D_Quad_Block::Flow_Type = FLOWTYPE_INVISCID;
