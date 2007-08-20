/**********************************************************************
 *      HighTemp2DQuad.cc: Subroutines for the 2D N-S HighTemp        *
 *                        quadrilateral solution block class.         *
 **********************************************************************/

#include "HighTemp2DQuad.h"

/**********************************************************************
 * HighTemp2D_Quad_Block -- Create storage for the static         *
 *                              variables.                            *
 **********************************************************************/
int HighTemp2D_Quad_Block::residual_variable = 1;
int HighTemp2D_Quad_Block::Number_of_Residual_Norms = 4;

