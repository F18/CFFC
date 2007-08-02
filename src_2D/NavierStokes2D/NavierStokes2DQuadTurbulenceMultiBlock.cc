/**********************************************************************
 * NavierStokes2DQuadTurbulenceMultiBlock.cc:                         *
 *                                                                    *
 * Multi-block versions of turbulence subroutines for 2D Navier-      *
 * Stokes multi-block quadrilateral mesh solution classes.            *
 *                                                                    *
 **********************************************************************/

// Include 2D NavierStokes quadrilateral mesh solution header file.

#ifndef _NAVIERSTOKES2D_QUAD_INCLUDED
#include "NavierStokes2DQuad.h"
#endif // _NAVIERSTOKES2D_QUAD_INCLUDED

/**********************************************************************
 * NavierStokes2D_Quad_Block -- Turbulence Multiple Block External    *
 *                              Subroutines.                          *
 **********************************************************************/

/**********************************************************************
 * Routine: Turbulence_BCs                                            *
 *                                                                    *
 * Apply turbulence boundary conditions (ie wall functions) at the    *
 * boundaries of a 1D array of 2D quadrilateral multi-block solution  *
 * blocks.                                                            *
 *                                                                    *
 **********************************************************************/
int Turbulence_BCs(NavierStokes2D_Quad_Block *Soln_ptr,
		   AdaptiveBlock2D_List &Soln_Block_List,
		   NavierStokes2D_Input_Parameters &Input_Parameters) {

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
