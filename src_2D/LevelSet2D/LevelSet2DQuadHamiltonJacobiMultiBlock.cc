/**********************************************************************
 * LevelSet2DQuadHamiltonJacobiMultiBlock.cc                          *
 *                                                                    *
 * Multi-block versions of subroutines for the solution of the 2D     *
 * Hamilton-Jacobi-type equations for the 2D Level Set multi-block    *
 * quadrilateral mesh solution classes.                               *
 *                                                                    *
 **********************************************************************/

// Include 2D LevelSet quadrilateral mesh solution header file.

#ifndef _LEVELSET2D_QUAD_INCLUDED
#include "LevelSet2DQuad.h"
#endif // _LEVELSET2D_QUAD_INCLUDED

/**********************************************************************
 * LevelSet2D_Quad_Block -- Hamilton-Jacobi Multiple Block External   *
 *                          Subroutines.                              *
 **********************************************************************/

/**********************************************************************
 * Routine: CFL_Hamilton_Jacobi                                       *
 *                                                                    *
 * Determines the allowable global and local time steps (for explicit *
 * LevelSet time stepping scheme) for a 1D array of 2D quadrilateral  *
 * multi-block solution blocks according to the                       *
 * Courant-Friedrichs-Lewy condition.                                 *
 *                                                                    *
 **********************************************************************/
double CFL_Hamilton_Jacobi(LevelSet2D_Quad_Block *Soln_ptr,
			   AdaptiveBlock2D_List &Soln_Block_List,
			   LevelSet2D_Input_Parameters &Input_Parameters) {

  double dtMin = MILLION;
  
  // Determine the allowable time step for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      dtMin = min(dtMin,CFL_Hamilton_Jacobi(Soln_ptr[nb],Input_Parameters));
    }
  }

  // Return the global time step.
  return dtMin;
  
}

/**********************************************************************
 * Routine: dUdt_Multistage_Hamilton_Jacobi                           *
 *                                                                    *
 * This routine evaluates the stage solution residual for the level   *
 * set equation for a 1D array of 2D quadrilateral multi-block        *
 * solution blocks.  A variety of multistage explicit time            *
 * integration and upwind finite-volume spatial discretization        *
 * procedures can be used depending on the specified input data.      *
 *                                                                    *
 **********************************************************************/
int dUdt_Multistage_Hamilton_Jacobi(LevelSet2D_Quad_Block *Soln_ptr,
				    AdaptiveBlock2D_List &Soln_Block_List,
				    LevelSet2D_Input_Parameters &Input_Parameters,
				    const int I_Stage) {

  int error_flag;

  // Update the solution for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = dUdt_Multistage_Hamilton_Jacobi(Soln_ptr[nb],
						   I_Stage,
						   Input_Parameters);
      if (error_flag) return error_flag;
    }
  }
  
  // Residuals for each quadrilateral multi-block solution block 
  // successfully calculated.
  return 0;

}

/**********************************************************************
 * Routine: Update_Solution_Multistage_Hamilton_Jacobi                *
 *                                                                    *
 * This routine updates the solution for a 1D array of 2D             *
 * quadrilateral multi-block Level Set solution blocks.  Second-order *
 * multistage explicit time integration and a finite-volume spatial   *
 * discretization procedure is used.                                  *
 *                                                                    *
 **********************************************************************/
int Update_Solution_Multistage_Hamilton_Jacobi(LevelSet2D_Quad_Block *Soln_ptr,
					       AdaptiveBlock2D_List &Soln_Block_List,
					       LevelSet2D_Input_Parameters &Input_Parameters,
					       const int I_Stage) {

  int error_flag;

  // Update the solution for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = Update_Solution_Multistage_Hamilton_Jacobi(Soln_ptr[nb],
							      I_Stage,
							      Input_Parameters);
      if (error_flag) return error_flag;
    }
  }
  
  // Quadrilateral multi-block solution blocks successfully updated.
  return 0;
 
}
