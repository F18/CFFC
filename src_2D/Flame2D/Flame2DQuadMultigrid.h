/**********************************************************************
 * Flame2DQuadMultigrid.h: Flame2D versions of some multigrid functions *
 *                        that could not use the templated version.   *
 **********************************************************************/

#ifndef _CHEM_QUAD_MULTIGRID_INCLUDED
#define _CHEM_QUAD_MULTIGRID_INCLUDED

// Include 2D Chemistry quadrilateral mesh solution header file.
#include "Flame2DQuad.h"

// Include 2D Chemistry quadrilateral mesh solution header file.
#include "../FASMultigrid2D/FASMultigrid2D.h"

/**********************************************************************
 * Flame2D_Quad_Block -- Flame2D Multigrid Subroutines.                 *
 **********************************************************************/



/**********************************************************************
 * Routine: CFL_Multigrid                                             *
 *                                                                    *
 * This routine sets the time step for each cell on the current       *
 * coarse grid level such that it is the minimum of the computed time *
 * step on the coarse grid and the time steps of the associated finer *
 * grid cells.                                                        *
 *                                                                    *
 **********************************************************************/
template <> inline void FAS_Multigrid2D_Solver<Flame2D_cState,
					       Flame2D_Quad_Block,				
					       Flame2D_Input_Parameters>::
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

	  //FLAME2D
	  Local_SolnBlk[Level_Coarse][b].dt[i_coarse][j_coarse] = Local_SolnBlk[Level_Fine][b].dt[i_fine][j_fine]; 

	} /* endfor */
      } /* endfor */
    } /* end if */     
  } /* endfor */
}

/**********************************************************************
 * Routine: Update_Primitive_Variables                                *
 *                                                                    *
 * This routine updates all primitive variables, W, from the          *
 * conserved variables, U, for all solution blocks on the given grid  *
 * level.                                                             *
 *                                                                    *
 **********************************************************************/
template <> inline void FAS_Multigrid2D_Solver<Flame2D_cState,
					       Flame2D_Quad_Block,
					       Flame2D_Input_Parameters>::
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
	  
	  // ADDED FOR FLAME2D to check for negative species concentrations
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
