/*! \file AdvectDiffuse2DState.cc
  @brief Subroutines for 2D Advection Diffusion Equation Quadrilateral Mesh Solution Classes. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "New_AdvectDiffuse2DQuad.h"        // AdvectDiffuse2D_Quad_Block class

/*******************************************************************************
 * Static Variable Initialization for AdvectDiffuse2D_Quad_Block_New
 ******************************************************************************/

/********************************************************************//**
 * AdvectDiffuse2D_Quad_Block_New -- Create storage for the static
 *                                   variables.                   
 **********************************************************************/
int AdvectDiffuse2D_Quad_Block_New::residual_variable = 1;


// Initialize ExactGrad
AdvectDiffuse2D_Quad_Block_New::Exact_Gradient_Function AdvectDiffuse2D_Quad_Block_New::ExactGrad = NULL;

// Initialize ExactSoln
FunctionType2D AdvectDiffuse2D_Quad_Block_New::ExactSoln = NULL;


/*******************************************************************************
 * AdvectDiffuse2D_Quad_Block -- Single Block Member Functions.
 ******************************************************************************/
