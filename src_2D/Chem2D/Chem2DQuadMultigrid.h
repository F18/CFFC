/**********************************************************************
 * Chem2DQuadMultigrid.h: Chem2D versions of some multigrid functions *
 *                        that could not use the templated version.   *
 **********************************************************************/

#ifndef _CHEM_QUAD_MULTIGRID_INCLUDED
#define _CHEM_QUAD_MULTIGRID_INCLUDED

// Include 2D Chemistry quadrilateral mesh solution header file.

#ifndef _CHEM_QUAD_INCLUDED
#include "Chem2DQuad.h"
#endif // _CHEM_QUAD_INCLUDED

// Include 2D Chemistry quadrilateral mesh solution header file.

#ifndef _FASMULTIGRID2D_INCLUDED
#include "../FASMultigrid2D/FASMultigrid2D.h"
#endif // _FASMULTIGRID2D_INCLUDED

/**********************************************************************
 * Chem2D_Quad_Block -- Chem2D Multigrid Subroutines.                 *
 **********************************************************************/

/**********************************************************************
 * FAS_Multigrid2D_Solver::Determine_Wall_Distance_on_Coarse_Grids -- *
 *                                                                    *
 * Perform distance to wall calculation on coarse grids as required.  *
 *                                                                    *
 **********************************************************************/
template <> int FAS_Multigrid2D_Solver<Chem2D_cState,
				       Chem2D_Quad_Block,
				       Chem2D_Input_Parameters>::
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
template <> int FAS_Multigrid2D_Solver<Chem2D_cState,
				       Chem2D_Quad_Block,
				       Chem2D_Input_Parameters>::
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
	  Local_SolnBlk[Level][nb].dUdt[i][j][k_residual].rhok = ZERO;
	  Local_SolnBlk[Level][nb].dUdt[i][j][k_residual].rhoomega = ZERO;
	}
      }
    }
  }

  // Residuals for k and omega zeroed.
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
template <> int FAS_Multigrid2D_Solver<Chem2D_cState,
				       Chem2D_Quad_Block,
				       Chem2D_Input_Parameters>::
Apply_Turbulence_Boundary_Conditions(const int &Level) {

  int error_flag;

  if (Level == FINEST_LEVEL) {
    cout << endl << " Need to apply turbulence boundary conditions.";
//     error_flag = Turbulence_BCs(Local_SolnBlk[Level],
// 			       List_of_Local_Solution_Blocks[Level],
// 			       *IP);
//     if (error_flag) return error_flag;
  }

  return 0;

}

template <> inline void FAS_Multigrid2D_Solver<Chem2D_cState,
					       Chem2D_Quad_Block,				
					       Chem2D_Input_Parameters>::
CFL_Multigrid(const int &Level_Coarse) {

  int i_fine,j_fine,Level_Fine;
  double dt_NE,dt_SE,dt_NW,dt_SW,A_coarse;
  Level_Fine = Level_Coarse - 1;

  // loop through used blocks only
  for (int b = 0 ; b <= List_of_Local_Solution_Blocks[Level_Fine].Nblk-1 ; ++b ) {
    if (List_of_Local_Solution_Blocks[Level_Fine].Block[b].used == ADAPTIVEBLOCK2D_USED) {

      /* Get number of ghost cells */
      int Nghost = Local_SolnBlk[Level_Fine][b].Nghost;

      /* Loop through coarse cells */
      for (int i_coarse = Local_SolnBlk[Level_Coarse][b].ICl; i_coarse <= Local_SolnBlk[Level_Coarse][b].ICu; i_coarse++) {
	for (int j_coarse = Local_SolnBlk[Level_Coarse][b].JCl; j_coarse <= Local_SolnBlk[Level_Coarse][b].JCu; j_coarse++) {
	  
	  i_fine = 2*(i_coarse-Nghost)+Nghost; // SW Corner i index (fine)
	  j_fine = 2*(j_coarse-Nghost)+Nghost; // SW Corner j index (fine)

	  // ORIGINAL 
// 	  A_coarse = Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].A;
// 	  dt_NE = sqrt(A_coarse/Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][j_fine+1].A)*
//                   Local_SolnBlk[Level_Fine][b].dt[i_fine+1][j_fine+1];
// 	  dt_SE = sqrt(A_coarse/Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][j_fine].A)*
//                   Local_SolnBlk[Level_Fine][b].dt[i_fine+1][j_fine];
// 	  dt_NW = sqrt(A_coarse/Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine+1].A)*
//                   Local_SolnBlk[Level_Fine][b].dt[i_fine][j_fine+1];
// 	  dt_SW = sqrt(A_coarse/Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].A)*
//                   Local_SolnBlk[Level_Fine][b].dt[i_fine][j_fine];
// 	 //  Local_SolnBlk[Level_Coarse][b].dt[i_coarse][j_coarse] = 
// //                   min(Local_SolnBlk[Level_Coarse][b].dt[i_coarse][j_coarse],
// //                       min(dt_NE,min(dt_SE,min(dt_NW,dt_SW))));

	  //CHEM2D
	  Local_SolnBlk[Level_Coarse][b].dt[i_coarse][j_coarse] = Local_SolnBlk[Level_Fine][b].dt[i_fine][j_fine]; 

	} /* endfor */
      } /* endfor */
    } /* end if */     
  } /* endfor */
}

template <> inline void FAS_Multigrid2D_Solver<Chem2D_cState,
					       Chem2D_Quad_Block,
					       Chem2D_Input_Parameters>::
Update_Primitive_Variables(const int &Level) {

  /* Loop through each block */
  for (int b = 0 ; b <= List_of_Local_Solution_Blocks[Level].Nblk-1 ; ++b ) {
    if (List_of_Local_Solution_Blocks[Level].Block[b].used == ADAPTIVEBLOCK2D_USED) {
      /* Loop through all interior cells - ghost cells are taken care of by BCs()
	 and message-passing functions */
      for (int i = Local_SolnBlk[Level][b].ICl; i <= Local_SolnBlk[Level][b].ICu; i++) {
	for (int j = Local_SolnBlk[Level][b].JCl; j <= Local_SolnBlk[Level][b].JCu; j++) {
	  
	  // ORIGINAL 
	  //Local_SolnBlk[Level][b].W[i][j] = W(Local_SolnBlk[Level][b].U[i][j]);
	  
	  // ADDED FOR CHEM2D to check for negative species concentrations
	  if(Local_SolnBlk[Level][b].U[i][j].negative_speccheck(0)){
	    Local_SolnBlk[Level][b].W[i][j] = W(Local_SolnBlk[Level][b].U[i][j]);
	  } else {
	    cerr<<"\n Failed negative_speccheck() in Update_Primitive_Variables() FASMultigrid";
	    exit(1);
	  }    	  
	}
      }
    } 
  }  
}

#endif // _CHEM_QUAD_MULTIGRID_INCLUDED
