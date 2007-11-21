/**********************************************************************
 * Flame2DQuad.cc: Subroutines for the 2D chem quadrilateral solution  *
 *                block class.                                        *
 **********************************************************************/

// Include 2D Chem quadrilateral solution block class.

#ifndef _FLAME2D_QUAD_INCLUDED
#include "Flame2DQuad.h"
#endif // _FLAME2D_QUAD_INCLUDED

/**********************************************************************
 * Flame2D_Quad_Block -- Create storage for the total number of        *
 *                      variables.                                    *
 **********************************************************************/
int Flame2D_Quad_Block::residual_variable = 1;
int Flame2D_Quad_Block::Number_of_Residual_Norms = 4;
//int Flame2D_Quad_Block::Flow_Type = FLOWTYPE_INVISCID;

//SNBCK data object
PlanckMean* Flame2D_Quad_Block::PlanckMean_data=NULL;
