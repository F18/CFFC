/*! \file NavierStokes2DQuad.cc
  @brief Subroutines for 2D Navier-Stokes quadrilateral solution block class. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "NavierStokes2DQuad.h"  // Include 2D NavierStokes quadrilateral solution block class.


/**********************************************************************
 * NavierStokes2D_Quad_Block -- Create storage for the static         *
 *                              variables.                            *
 **********************************************************************/
int NavierStokes2D_Quad_Block::residual_variable = 1;
#ifdef _NS_PARALLEL_DEBUG_
ofstream NavierStokes2D_Quad_Block::dout;
#endif
int NavierStokes2D_Quad_Block::Number_of_Residual_Norms = 4;
// Initialize RefW
NavierStokes2D_pState NavierStokes2D_Quad_Block::RefW(1.0);
// Initialize ExactSoln
NavierStokes2D_ExactSolutions *NavierStokes2D_Quad_Block::ExactSoln = NULL;

/***********************************************************************
 * NavierStokes2D_Quad_Block -- Single Block Member Functions.         *
 **********************************************************************/

/*****************************************************//**
 * Copy the solution information of quadrilateral solution 
 * block SolnBlk to the current solution block.
 ********************************************************/
NavierStokes2D_Quad_Block & NavierStokes2D_Quad_Block::operator =(const NavierStokes2D_Quad_Block &Soln){

  // Handle self-assignment:
  if (this == & Soln) return *this;

  return *this;
}
