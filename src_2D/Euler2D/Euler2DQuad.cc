/*! \file Euler2DQuad.cc
  @brief Subroutines for 2D Euler Quadrilateral Mesh Solution Class. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "Euler2DQuad.h"                // Euler2D_Quad_Block class

/**********************************************************************
 * Euler2D_Quad_Block -- Create storage for the static variables.     *
 **********************************************************************/
// Initialize residual_variable
int Euler2D_Quad_Block::residual_variable = 1;
// Initialize Number_of_Residual_Norms
int Euler2D_Quad_Block::Number_of_Residual_Norms = 1;
// Initialize Flow_Type
int Euler2D_Quad_Block::Flow_Type = FLOWTYPE_INVISCID;
// Initialize RefU
Euler2D_pState Euler2D_Quad_Block::RefU(1.0);
// Initialize ExactSoln
Euler2D_ExactSolutions *Euler2D_Quad_Block::ExactSoln = NULL;



