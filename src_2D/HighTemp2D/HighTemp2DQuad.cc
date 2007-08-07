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

enum current_solver_types HighTemp2D_Quad_Block::current_solver_type = CST_EXPLICIT;

