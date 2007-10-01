/*******************************************************************************
 *
 * Long member functions for class AdvectDiffuse2D_Quad_Block
 *
 ******************************************************************************/

//=====Included Files=====//

//-----Standard Library-----//

//-----External-----//

//-----Internal-----//

#include "AdvectDiffuse2DQuad.h"        // AdvectDiffuse2D_Quad_Block class

//=====End of Includes=====//


/*******************************************************************************
 * Static Variable Initialization for AdvectDiffuse2D_Quad_Block
 ******************************************************************************/

/**********************************************************************
 * AdvectDiffuse2D_Quad_Block -- Create storage for the static        *
 *                               variables.                           *
 **********************************************************************/
int AdvectDiffuse2D_Quad_Block::residual_variable = 1;


// Initialize ExactGrad
AdvectDiffuse2D_Quad_Block::Exact_Gradient_Function AdvectDiffuse2D_Quad_Block::ExactGrad = NULL;

// Initialize ExactSoln
FunctionType2D AdvectDiffuse2D_Quad_Block::ExactSoln = NULL;


/*******************************************************************************
 * AdvectDiffuse2D_Quad_Block -- Single Block Member Functions.
 ******************************************************************************/
