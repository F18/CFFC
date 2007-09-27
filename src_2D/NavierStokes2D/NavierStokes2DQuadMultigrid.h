/**********************************************************************
 * NavierStokes2DQuadMultigrid.h: NavierStokes2D versions of some     *
 *                                multigrid functions that could not  *
 *                                use templated version.              *
 **********************************************************************/

#ifndef _NAVIERSTOKES2D_QUAD_MULTIGRID_INCLUDED
#define _NAVIERSTOKES2D_QUAD_MULTIGRID_INCLUDED

// Include 2D Navier-Stokes quadrilateral mesh solution header file.

#ifndef _NAVIERSTOKES2D_QUAD_INCLUDED
#include "NavierStokes2DQuad.h"
#endif // _NAVIERSTOKES2D_QUAD_INCLUDED

// Include 2D Navier-Stokes quadrilateral mesh solution header file.

#ifndef _FASMULTIGRID2D_INCLUDED
#include "../FASMultigrid2D/FASMultigrid2D.h"
#endif // _FASMULTIGRID2D_INCLUDED

/**********************************************************************
 * NavierStokes2D_Quad_Block -- NavierStokes2D Multigrid Subroutines. *
 **********************************************************************/

/**********************************************************************
 * FAS_Multigrid2D_Solver::Determine_Wall_Distance_on_Coarse_Grids -- *
 *                                                                    *
 * Perform distance to wall calculation on coarse grids as required.  *
 *                                                                    *
 **********************************************************************/
template <> int FAS_Multigrid2D_Solver<NavierStokes2D_cState,
				       NavierStokes2D_Quad_Block,
				       NavierStokes2D_Input_Parameters>::
Determine_Wall_Distance_on_Coarse_Grids(void) {

  // Exit immediately if not a turbulent flow.
  if (IP->FlowType != FLOWTYPE_TURBULENT_RANS_K_OMEGA) return 0;

  int error_flag;

  // Determine the wall distance and wall distance location if required.
  for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {
    error_flag = Determine_Wall_Distance(Local_SolnBlk[level],
					 *QuadTree,
					 List_of_Local_Solution_Blocks[level],
					 *IP);
    if (error_flag) return error_flag;
  }

  // All required functions computed successfully.
  return 0;

}

/**********************************************************************
 * FAS_Multigrid2D_Solver::Zero_Residuals_on_Coarse_Grid --           *
 *                                                                    *
 * Zero the required residuals on coarse grids.                       *
 *                                                                    *
 **********************************************************************/
template <> int FAS_Multigrid2D_Solver<NavierStokes2D_cState,
				       NavierStokes2D_Quad_Block,
				       NavierStokes2D_Input_Parameters>::
Zero_Residuals_on_Coarse_Grid(const int &Level,
			      const int &i_stage) {

  // Exit immediately if not a turbulent flow.
  if (IP->FlowType != FLOWTYPE_TURBULENT_RANS_K_OMEGA) return 0;

  int k_residual;

  switch(IP->i_Time_Integration) {
  case TIME_STEPPING_EXPLICIT_EULER :
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
    k_residual = 0;
    if (IP->N_Stage == 4) {
      if (i_stage == 4) k_residual = 0;
      else k_residual = i_stage - 1;
    }
    return 7777777;
    break;
  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
    k_residual = 0;
    break;
  default:
    k_residual = 0;
    break;
  };

  // Zero the required residuals on each solution block in use.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = Local_SolnBlk[Level][nb].JCl; j <= Local_SolnBlk[Level][nb].JCu; j++) {
	for (int i = Local_SolnBlk[Level][nb].ICl; i <= Local_SolnBlk[Level][nb].ICu; i++) {
	  Local_SolnBlk[Level][nb].dUdt[i][j][k_residual].dk = ZERO;
	  Local_SolnBlk[Level][nb].dUdt[i][j][k_residual].domega = ZERO;
	  Local_SolnBlk[Level][nb].dUdt[i][j][k_residual].dke = ZERO;
	  Local_SolnBlk[Level][nb].dUdt[i][j][k_residual].dee = ZERO;
	}
      }
    }
  }

  // Residuals for k and omega and ke and ee zeroed.
  return 0;

}

/**********************************************************************
 * Routine: Apply_Turbulence_Boundary_Conditions                      *
 *                                                                    *
 * Apply turbulence boundary conditions on finest mesh.  Explicit     *
 * specialization of this routine is required for the solution of     *
 * systems of equations that include turbulent flow.                  *
 *                                                                    *
 **********************************************************************/
template <> int FAS_Multigrid2D_Solver<NavierStokes2D_cState,
				       NavierStokes2D_Quad_Block,
				       NavierStokes2D_Input_Parameters>::
Apply_Turbulence_Boundary_Conditions(const int &Level) {

  if (Level != FINEST_LEVEL) return 0;

  return Turbulence_BCs(Local_SolnBlk[Level],
			List_of_Local_Solution_Blocks[Level],
			*IP);

}

#endif // _NAVIERSTOKES2D_QUAD_MULTIGRID_INCLUDED
