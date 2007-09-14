/* Grid2DQuad.cc:  Member functions implementation for 
                   2D quadrilateral block grid class. */

/* Include 2D quadrilateral block grid type header file. */

#ifndef _GRID2D_QUAD_BLOCK_INCLUDED
#include "Grid2DQuad.h"
#endif // _GRID2D_QUAD_BLOCK_INCLUDED

/*************************************************************************
 * Grid2D_Quad_Block -- Static constant data memebers.                   *
 *************************************************************************/
const double Grid2D_Quad_Block::Gauss2QuadPoints_CoeffPoint1 = 0.2113248654051871177454256; /* 0.5*(1 - 1/sqrt(3)) */
const double Grid2D_Quad_Block::Gauss2QuadPoints_CoeffPoint2 = 0.7886751345948128822545744; /* 0.5*(1 + 1/sqrt(3)) */
