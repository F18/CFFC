/**********************************************************************
 * Dusty2DQuadTurbulenceMultiBlock.cc:                                *
 *                                                                    *
 * Multi-block versions of turbulence subroutines for 2D Navier-      *
 * Stokes multi-block quadrilateral mesh solution classes.            *
 *                                                                    *
 **********************************************************************/

// Include 2D Dusty quadrilateral mesh solution header file.

#ifndef _DUSTY2D_QUAD_INCLUDED
#include "Dusty2DQuad.h"
#endif // _DUSTY2D_QUAD_INCLUDED

/**********************************************************************
 * Dusty2D_Quad_Block -- Turbulence Multiple Block External           *
 *                       Subroutines.                                 *
 **********************************************************************/

/**********************************************************************
 * Routine: Turbulence_BCs                                            *
 *                                                                    *
 * Apply turbulence boundary conditions (ie wall functions) at the    *
 * boundaries of a 1D array of 2D quadrilateral multi-block solution  *
 * blocks.                                                            *
 *                                                                    *
 **********************************************************************/
int Turbulence_BCs(Dusty2D_Quad_Block *Soln_ptr,
		   AdaptiveBlock2D_List &Soln_Block_List,
		   Dusty2D_Input_Parameters &Input_Parameters) {

  // Exit immediately if not a turbulent flow.
  if (Input_Parameters.FlowType != FLOWTYPE_TURBULENT_RANS_K_OMEGA) return 0;

  int error_flag;

  // Set turbulent boundary data for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Turbulence_BCs(Soln_ptr[nb],
				  Input_Parameters);
      if (error_flag) return error_flag;
    }
  }

  // Turbulent boundary conditions successfully applied.
  return 0;

}
