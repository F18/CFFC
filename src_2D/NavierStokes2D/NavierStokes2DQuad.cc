/**********************************************************************
 * NavierStokes2DQuad.cc: Subroutines for the 2D Navier-Stokes        *
 *                        quadrilateral solution block class.         *
 **********************************************************************/

// Include 2D NavierStokes quadrilateral solution block class.

#ifndef _NAVIERSTOKES2D_QUAD_INCLUDED
#include "NavierStokes2DQuad.h"
#endif // _NAVIERSTOKES2D_QUAD_INCLUDED

/**********************************************************************
 * NavierStokes2D_Quad_Block -- Create storage for the static         *
 *                              variables.                            *
 **********************************************************************/
int NavierStokes2D_Quad_Block::residual_variable = 1;
#ifdef _NS_PARALLEL_DEBUG_
ofstream NavierStokes2D_Quad_Block::dout;
#endif
int NavierStokes2D_Quad_Block::Number_of_Residual_Norms = 4;
