/**********************************************************************
 **********************************************************************
 **                                                                  **
 ** File: Rte2DQuadMultigrid.cc                                      **
 **                                                                  **
 ** Description: FAS Multigrid 2D Explicit specialization for Rte2D. **
 **              This module defines the overloaded functions.       **
 **                                                                  **
 ** Author: Marc "T-Bone" Charest                                    **
 **                                                                  **
 ** Revision:  Date        Initials   Change                         **
 **            04/24/2007  MRC        Original creation              **
 **                                                                  **
 **********************************************************************
 **********************************************************************/
#ifndef _RTE2D_MULTIGRID_INCLUDED 
#define _RTE2D_MULTIGRID_INCLUDED 

/* Include 2D Rte equation quadrilateral mesh solution header file. */
#include "Rte2DQuad.h"

/* Include 2D Rte equation input header file. */
#include "Rte2DInput.h"

/* Include FAS Multigrid header file */
#include "../FASMultigrid2D/FASMultigrid2D.h"


/**********************************************************************
 * Specialization: Additional_Solution_Block_Setup                    *
 *                                                                    *
 * Perform some additional setup on the local solution block before   *
 * beginning Multi Grid computions.                                   *
 *                                                                    *
 **********************************************************************/
template <> void FAS_Multigrid2D_Solver<Rte2D_State,
                                        Rte2D_Quad_Block,
                                        Rte2D_Input_Parameters>::
Additional_Solution_Block_Setup(Rte2D_Quad_Block &SolnBlk) 
{
  SolnBlk.ScaleGridTo3D(IP->Axisymmetric);
}



/**********************************************************************
 * Routine: Restrict_Solution_Blocks (for Multigrid)                  *
 *                                                                    *
 * Restrict solution from Level_Fine to Level_Coarse for all blocks   *
 * on the local solution block list.  Note that the solution at the   *
 * coarse grid level is overwritten.                                  *
 *                                                                    *
 **********************************************************************/
template <> void FAS_Multigrid2D_Solver<Rte2D_State,
                                        Rte2D_Quad_Block,
                                        Rte2D_Input_Parameters>::
Restrict_Solution_Blocks(const int &Level_Fine) {

  int i_fine, j_fine, Nghost, nghost;
  int Level_Coarse = Level_Fine + 1;

  // Determine if the restriction includes the ghost cells.
  if (IP->Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions) nghost = 0;
  else nghost = 1;

  // Loop through each solution block.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level_Fine].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level_Fine].Block[nb].used == ADAPTIVEBLOCK2D_USED) { 

      // Get the number of ghost cells on the fine grid.
      Nghost = Local_SolnBlk[Level_Fine][nb].Nghost;

      // Loop through the coarse grid cells.
      for (int i_coarse = Local_SolnBlk[Level_Coarse][nb].ICl-nghost; i_coarse <= Local_SolnBlk[Level_Coarse][nb].ICu+nghost; i_coarse++) {
	for (int j_coarse = Local_SolnBlk[Level_Coarse][nb].JCl-nghost; j_coarse <= Local_SolnBlk[Level_Coarse][nb].JCu+nghost; j_coarse++) {
	  // Determine the (i,j) index of the SW corner fine cell.
	  i_fine = 2*(i_coarse-Nghost)+Nghost;
	  j_fine = 2*(j_coarse-Nghost)+Nghost;
	  // Determine the solution state of the coarse grid cell by
	  // a area-weighted average of the associated fine grid cells.

	  /***********************************************************************
	   *************************** RTE SPECIFIC ******************************
	  // ORIGINAL : 
	  Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse] = 
	    (Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine]*
	     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].A +
	     Local_SolnBlk[Level_Fine][nb].U[i_fine+1][j_fine]*
	     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine].A +
	     Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine+1]*
	     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine+1].A +
	     Local_SolnBlk[Level_Fine][nb].U[i_fine+1][j_fine+1]*
	     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine+1].A)/
 	    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].A;
	  *************************************************************************/
	  // RTE2D: Added Quasi-3D scaling param 
	  Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse] = 
	    (Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] *
	     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].A * 
	     Local_SolnBlk[Level_Fine][nb].Sp[i_fine][j_fine] +
	     Local_SolnBlk[Level_Fine][nb].U[i_fine+1][j_fine] *
	     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine].A * 
	     Local_SolnBlk[Level_Fine][nb].Sp[i_fine+1][j_fine] +
	     Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine+1] *
	     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine+1].A * 
	     Local_SolnBlk[Level_Fine][nb].Sp[i_fine][j_fine+1] +
	     Local_SolnBlk[Level_Fine][nb].U[i_fine+1][j_fine+1] *
	     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine+1].A * 
	     Local_SolnBlk[Level_Fine][nb].Sp[i_fine+1][j_fine+1]) /
	    (Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].A * 
	     Local_SolnBlk[Level_Coarse][nb].Sp[i_coarse][j_coarse]);	
	  /*************************************************************************
	   *************************************************************************/

	}
      }

    }
  }

}



/**********************************************************************
 * Routine: Restrict_Residuals (for Multigrid)                        *
 *                                                                    *
 * This routine restricts the residual dUdt from Level_Fine to        *
 * Level_Coarse (stored in P) for all blocks on the local solution    *
 * block list.  Note that the residual values on the coarse level are *
 * overwritten by this routine.                                       *
 *                                                                    *
 **********************************************************************/
template <> void FAS_Multigrid2D_Solver<Rte2D_State,
		  		        Rte2D_Quad_Block,
					Rte2D_Input_Parameters>::
Restrict_Residuals(const int &Level_Fine) {

  int i_fine, j_fine, Nghost, nghost;
  int Level_Coarse = Level_Fine + 1;

  // Determine if the restriction includes the ghost cells.
  if (IP->Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions) nghost = 0;
  else nghost = 1;

  // Loop through each local solution block.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level_Fine].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level_Fine].Block[nb].used == ADAPTIVEBLOCK2D_USED) { 

      // Get the number of ghost cells.
      Nghost = Local_SolnBlk[Level_Fine][nb].Nghost;

      // Loop through the coarse cells.
      for (int i_coarse = Local_SolnBlk[Level_Coarse][nb].ICl-nghost; i_coarse <= Local_SolnBlk[Level_Coarse][nb].ICu+nghost; i_coarse++) {
	for (int j_coarse = Local_SolnBlk[Level_Coarse][nb].JCl-nghost; j_coarse <= Local_SolnBlk[Level_Coarse][nb].JCu+nghost; j_coarse++) {
	  // Determine the (i,j) index of the corresponding SW corner
	  // fine cell.
	  i_fine = 2*(i_coarse-Nghost)+Nghost;
	  j_fine = 2*(j_coarse-Nghost)+Nghost;
	  // Restrict the residuals.
	  /***********************************************************************
	   *************************** RTE SPECIFIC ******************************
	  // ORIGINAL : 
	  MG_SolnBlk[Level_Coarse][nb].P[i_coarse][j_coarse] = 
	    (Local_SolnBlk[Level_Fine][nb].dUdt[i_fine][j_fine][0]*
	     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].A +
	     Local_SolnBlk[Level_Fine][nb].dUdt[i_fine+1][j_fine][0]*
	     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine].A +
	     Local_SolnBlk[Level_Fine][nb].dUdt[i_fine][j_fine+1][0]*
	     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine+1].A +
	     Local_SolnBlk[Level_Fine][nb].dUdt[i_fine+1][j_fine+1][0]*
	     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine+1].A)/
 	    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].A;
	  *************************************************************************/
	  // RTE2D: Added Quasi-3D scaling param 
	  MG_SolnBlk[Level_Coarse][nb].P[i_coarse][j_coarse] = 
	    (Local_SolnBlk[Level_Fine][nb].dUdt[i_fine][j_fine][0] *
	     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].A *
	     Local_SolnBlk[Level_Fine][nb].Sp[i_fine][j_fine] +
	     Local_SolnBlk[Level_Fine][nb].dUdt[i_fine+1][j_fine][0] *
	     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine].A *
	     Local_SolnBlk[Level_Fine][nb].Sp[i_fine+1][j_fine] +
	     Local_SolnBlk[Level_Fine][nb].dUdt[i_fine][j_fine+1][0] *
	     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine+1].A *
	     Local_SolnBlk[Level_Fine][nb].Sp[i_fine][j_fine+1] +
	     Local_SolnBlk[Level_Fine][nb].dUdt[i_fine+1][j_fine+1][0] *
	     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine+1].A *
	     Local_SolnBlk[Level_Fine][nb].Sp[i_fine+1][j_fine+1] ) /
	    (Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].A * 
	     Local_SolnBlk[Level_Coarse][nb].Sp[i_coarse][j_coarse] );
	/*************************************************************************
	 *************************************************************************/
	}
      }

    }
  }

}




/**********************************************************************
 * Routine: Restrict_Boundary_Ref_States (for Multigrid)              *
 *                                                                    *
 * This routine restricts the boundary refeference states (UoN/S/E/W) *
 * from Level_Fine to Level_Coarse for all solution blocks.  The      *
 * values on the coarse grid level are overwritten by this routine.   *
 * The restriction operator used is area weighted average.            *
 *                                                                    *
 **********************************************************************/
template <> void FAS_Multigrid2D_Solver<Rte2D_State,
		  		        Rte2D_Quad_Block,
				        Rte2D_Input_Parameters>::
Restrict_Boundary_Ref_States(const int &Level_Fine) {

  int i_fine, j_fine, ICl, ICu, JCl, JCu, Nghost;
  int Level_Coarse = Level_Fine + 1;
  
  // Loop through each local solution block.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level_Fine].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level_Fine].Block[nb].used == ADAPTIVEBLOCK2D_USED) { 

      // Get the number of ghost cells.
      Nghost = Local_SolnBlk[Level_Fine][nb].Nghost;

      // Store the lower and upper j-indices on the fine grid.
      JCl = Local_SolnBlk[Level_Fine][nb].JCl;
      JCu = Local_SolnBlk[Level_Fine][nb].JCu;

      // Loop through the i-direction cells of the coarse grid.
      for (int i_coarse = Local_SolnBlk[Level_Coarse][nb].ICl; i_coarse <= Local_SolnBlk[Level_Coarse][nb].ICu; i_coarse++) {
	// Determine the i-index of the corresponding SW corner
	// fine cell.
	i_fine = 2*(i_coarse-Nghost)+Nghost;

	/***********************************************************************
	 *************************** RTE SPECIFIC ******************************
	// ORIGINAL : 
	// Restrict the north boundary reference state.
	Local_SolnBlk[Level_Coarse][nb].WoN[i_coarse] = 
	  (Local_SolnBlk[Level_Fine][nb].WoN[i_fine]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][JCu].A +
	   Local_SolnBlk[Level_Fine][nb].WoN[i_fine+1]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][JCu].A)/
	  (Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][JCu].A +
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][JCu].A);
	// Restrict the south boundary reference state.
	Local_SolnBlk[Level_Coarse][nb].WoS[i_coarse] = 
	  (Local_SolnBlk[Level_Fine][nb].WoS[i_fine]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][JCl].A +
	   Local_SolnBlk[Level_Fine][nb].WoS[i_fine+1]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][JCl].A)/
	  (Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][JCl].A +
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][JCl].A);
	*************************************************************************/
	// RTE2D - added grid thickness scaling parameter
	//       - apply to conserved state 
	// Restrict the north boundary reference state.
	Local_SolnBlk[Level_Coarse][nb].UoN[i_coarse] = 
	  (Local_SolnBlk[Level_Fine][nb].UoN[i_fine] *
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][JCu].A * 
	   Local_SolnBlk[Level_Fine][nb].Sp[i_fine][JCu] +
	   Local_SolnBlk[Level_Fine][nb].UoN[i_fine+1] *
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][JCu].A * 
	   Local_SolnBlk[Level_Fine][nb].Sp[i_fine+1][JCu]) /
	  (Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][JCu].A * 
	   Local_SolnBlk[Level_Fine][nb].Sp[i_fine][JCu] +
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][JCu].A * 
	   Local_SolnBlk[Level_Fine][nb].Sp[i_fine+1][JCu]);

	// Restrict the south boundary reference state.
	Local_SolnBlk[Level_Coarse][nb].UoS[i_coarse] = 
	  (Local_SolnBlk[Level_Fine][nb].UoS[i_fine] *
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][JCl].A * 
	   Local_SolnBlk[Level_Fine][nb].Sp[i_fine][JCl] +
	   Local_SolnBlk[Level_Fine][nb].UoS[i_fine+1] *
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][JCl].A * 
	   Local_SolnBlk[Level_Fine][nb].Sp[i_fine+1][JCl]) /
	  (Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][JCl].A * 
	   Local_SolnBlk[Level_Fine][nb].Sp[i_fine][JCl] +
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][JCl].A * 
	   Local_SolnBlk[Level_Fine][nb].Sp[i_fine+1][JCl]);
	/*************************************************************************
	 *************************************************************************/
      }

      // Store the lower and upper i-indices on the fine grid.
      ICl = Local_SolnBlk[Level_Fine][nb].ICl;
      ICu = Local_SolnBlk[Level_Fine][nb].ICu;

      // Loop through the j-direction cells of the coarse grid.
      for (int j_coarse = Local_SolnBlk[Level_Coarse][nb].JCl; j_coarse <= Local_SolnBlk[Level_Coarse][nb].JCu; j_coarse++) {
	// Determine the j-index of the corresponding SW corner
	// fine cell.
	j_fine = 2*(j_coarse-Nghost)+Nghost;

	/***********************************************************************
	 *************************** RTE SPECIFIC ******************************
	// ORIGINAL : 
	// Restrict the west boundary reference state.
	Local_SolnBlk[Level_Coarse][nb].WoW[j_coarse] = 
	  (Local_SolnBlk[Level_Fine][nb].WoW[j_fine]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICl][j_fine].A +
	   Local_SolnBlk[Level_Fine][nb].WoW[j_fine+1]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICl][j_fine+1].A)/
	  (Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICl][j_fine].A +
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICl][j_fine+1].A);
	// Restrict the east boundary reference state.
	Local_SolnBlk[Level_Coarse][nb].WoE[j_coarse] = 
	  (Local_SolnBlk[Level_Fine][nb].WoE[j_fine]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICu][j_fine].A +
	   Local_SolnBlk[Level_Fine][nb].WoE[j_fine+1]*
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICu][j_fine+1].A)/
	  (Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICu][j_fine].A +
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICu][j_fine+1].A);
	*************************************************************************/
	// RTE2D - added grid thickness scaling parameter : 
	//       - apply to conserved state 
	// Restrict the west boundary reference state.
	Local_SolnBlk[Level_Coarse][nb].UoW[j_coarse] = 
	  (Local_SolnBlk[Level_Fine][nb].UoW[j_fine] *
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICl][j_fine].A * 
	   Local_SolnBlk[Level_Fine][nb].Sp[ICl][j_fine] +
	   Local_SolnBlk[Level_Fine][nb].UoW[j_fine+1] *
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICl][j_fine+1].A * 
	   Local_SolnBlk[Level_Fine][nb].Sp[ICl][j_fine+1]) /
	  (Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICl][j_fine].A * 
	   Local_SolnBlk[Level_Fine][nb].Sp[ICl][j_fine] +
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICl][j_fine+1].A * 
	   Local_SolnBlk[Level_Fine][nb].Sp[ICl][j_fine+1]);

	// Restrict the east boundary reference state.
	Local_SolnBlk[Level_Coarse][nb].UoE[j_coarse] = 
	  (Local_SolnBlk[Level_Fine][nb].UoE[j_fine] *
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICu][j_fine].A * 
	   Local_SolnBlk[Level_Fine][nb].Sp[ICu][j_fine] +
	   Local_SolnBlk[Level_Fine][nb].UoE[j_fine+1] *
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICu][j_fine+1].A * 
	   Local_SolnBlk[Level_Fine][nb].Sp[ICu][j_fine+1]) /
	  (Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICu][j_fine].A * 
	   Local_SolnBlk[Level_Fine][nb].Sp[ICu][j_fine] +
	   Local_SolnBlk[Level_Fine][nb].Grid.Cell[ICu][j_fine+1].A * 
	   Local_SolnBlk[Level_Fine][nb].Sp[ICu][j_fine+1]);
	/*************************************************************************
	 *************************************************************************/
      }

    }
  }

}





/**********************************************************************
 * Routine: CFL_Multigrid                                             *
 *                                                                    *
 * This routine sets the time step for each cell on the current       *
 * coarse grid level such that it is the minimum of the computed time *
 * step on the coarse grid and the time steps of the associated finer *
 * grid cells.                                                        *
 *                                                                    *
 **********************************************************************/
template <> void FAS_Multigrid2D_Solver<Rte2D_State,
		  		        Rte2D_Quad_Block,
				        Rte2D_Input_Parameters>::
CFL_Multigrid(const int &Level_Coarse) {

  int i_fine, j_fine, Level_Fine, Nghost;
  double dt_NE, dt_SE, dt_NW, dt_SW, A_coarse;
  Level_Fine = Level_Coarse - 1;

  // Loop through each local solution block.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level_Fine].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level_Fine].Block[nb].used == ADAPTIVEBLOCK2D_USED) {

      // Get the number of ghost cells.
      Nghost = Local_SolnBlk[Level_Fine][nb].Nghost;

      // Loop through the coarse cells.
      for (int i_coarse = Local_SolnBlk[Level_Coarse][nb].ICl; i_coarse <= Local_SolnBlk[Level_Coarse][nb].ICu; i_coarse++) {
	for (int j_coarse = Local_SolnBlk[Level_Coarse][nb].JCl; j_coarse <= Local_SolnBlk[Level_Coarse][nb].JCu; j_coarse++) {
	  // Determine the (i,j) index of the corresponding SW corner
	  // fine cell.
	  i_fine = 2*(i_coarse-Nghost)+Nghost;
	  j_fine = 2*(j_coarse-Nghost)+Nghost;
	  /***********************************************************************
	   *************************** RTE SPECIFIC ******************************
	  // ORIGINAL : 
	  // Determine the area of the coarse grid cell.
	  A_coarse = Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].A;
	  // Determine the time-steps of each of the associated fine
	  // grid cells.
	  dt_NE = sqrt(A_coarse/Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine+1].A)*
                  Local_SolnBlk[Level_Fine][nb].dt[i_fine+1][j_fine+1];
	  dt_SE = sqrt(A_coarse/Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine].A)*
                  Local_SolnBlk[Level_Fine][nb].dt[i_fine+1][j_fine];
	  dt_NW = sqrt(A_coarse/Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine+1].A)*
                  Local_SolnBlk[Level_Fine][nb].dt[i_fine][j_fine+1];
	  dt_SW = sqrt(A_coarse/Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].A)*
                  Local_SolnBlk[Level_Fine][nb].dt[i_fine][j_fine];
	  *************************************************************************/
	  // RTE2D: Added Quasi-3D scaling param 
	  // Determine the area of the coarse grid cell.
	  A_coarse = Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].A*
	             Local_SolnBlk[Level_Coarse][nb].Sp[i_coarse][j_coarse];

	  // Determine the time-steps of each of the associated fine
	  // grid cells.
	  dt_NE = sqrt(A_coarse/(Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine+1].A*
				 Local_SolnBlk[Level_Fine][nb].Sp[i_fine+1][j_fine+1])) *
	          Local_SolnBlk[Level_Fine][nb].dt[i_fine+1][j_fine+1];
	  dt_SE = sqrt(A_coarse/(Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine].A*
				 Local_SolnBlk[Level_Fine][nb].Sp[i_fine+1][j_fine])) *
	          Local_SolnBlk[Level_Fine][nb].dt[i_fine+1][j_fine];
	  dt_NW = sqrt(A_coarse/(Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine+1].A*
				 Local_SolnBlk[Level_Fine][nb].Sp[i_fine][j_fine+1])) *
	          Local_SolnBlk[Level_Fine][nb].dt[i_fine][j_fine+1];
	  dt_SW = sqrt(A_coarse/(Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].A*
				 Local_SolnBlk[Level_Fine][nb].Sp[i_fine][j_fine])) *
   	          Local_SolnBlk[Level_Fine][nb].dt[i_fine][j_fine];
	  /*************************************************************************
	   *************************************************************************/
	  // Determine the coarse grid cell time-step by choosing the 
	  // minimum of the fine grid cell time-steps.
	  Local_SolnBlk[Level_Coarse][nb].dt[i_coarse][j_coarse] = 
	    min(Local_SolnBlk[Level_Coarse][nb].dt[i_coarse][j_coarse],
		min(dt_NE,min(dt_SE,min(dt_NW,dt_SW))) );
	}
      }

    }
  }

}




/**********************************************************************
 * Routine: Update_Primitive_Variables                                *
 *                                                                    *
 * This routine updates all primitive variables, W, from the          *
 * conserved variables, U, for all solution blocks on the given grid  *
 * level.                                                             *
 *                                                                    *
 **********************************************************************/
template <> void FAS_Multigrid2D_Solver<Rte2D_State,
			 	        Rte2D_Quad_Block,
				        Rte2D_Input_Parameters>::
Update_Primitive_Variables(const int &Level) {  /* DO NOTHING */ }



/**********************************************************************
 * Routine: Restrict_NonSolution_Blocks                               *
 *                                                                    *
 * Restrict solution from Level_Fine to Level_Coarse for all blocks   *
 * on the local solution block list.  Note that the solution at the   *
 * coarse grid level is overwritten.  This is for the medium state    *
 * which is not part of the overall solution state.  This function    *
 * only called to initialize the coarse grid solutions.               *
 *                                                                    *
 **********************************************************************/
void Restrict_NonSolution_Blocks( Rte2D_Quad_Block **Local_SolnBlk,
				  AdaptiveBlock2D_List *List_of_Local_Solution_Blocks,
				  const int &Level_Fine )
{
				  
  //
  // Declares
  //
  int i_fine, j_fine, Nghost;

  // the coarse level
  int Level_Coarse = Level_Fine + 1;

  // Determine if the restriction includes the ghost cells.
  // For the medium state, we always include the ghost cells
  //if (IP->Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions) nghost = 0;
  //else nghost = 1;
  int nghost(1);

  //
  // Loop through each solution block.
  //
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level_Fine].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level_Fine].Block[nb].used == ADAPTIVEBLOCK2D_USED) { 

      // Get the number of ghost cells on the fine grid.
      Nghost = Local_SolnBlk[Level_Fine][nb].Nghost;

      //
      // Loop through the coarse grid cells.
      //
      for (int i_coarse = Local_SolnBlk[Level_Coarse][nb].ICl-nghost; 
	   i_coarse <= Local_SolnBlk[Level_Coarse][nb].ICu+nghost; 
	   i_coarse++) 
	for (int j_coarse = Local_SolnBlk[Level_Coarse][nb].JCl-nghost; 
	     j_coarse <= Local_SolnBlk[Level_Coarse][nb].JCu+nghost; 
	     j_coarse++) {

	  // Determine the (i,j) index of the SW corner fine cell.
	  i_fine = 2*(i_coarse-Nghost)+Nghost;
	  j_fine = 2*(j_coarse-Nghost)+Nghost;

	  // Determine the solution state of the coarse grid cell by
	  // a area-weighted average of the associated fine grid cells.
	  Local_SolnBlk[Level_Coarse][nb].M[i_coarse][j_coarse] = 
	    (Local_SolnBlk[Level_Fine][nb].M[i_fine][j_fine] *
	     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].A * 
	     Local_SolnBlk[Level_Fine][nb].Sp[i_fine][j_fine] +
	     Local_SolnBlk[Level_Fine][nb].M[i_fine+1][j_fine] *
	     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine].A * 
	     Local_SolnBlk[Level_Fine][nb].Sp[i_fine+1][j_fine] +
	     Local_SolnBlk[Level_Fine][nb].M[i_fine][j_fine+1] *
	     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine+1].A * 
	     Local_SolnBlk[Level_Fine][nb].Sp[i_fine][j_fine+1] +
	     Local_SolnBlk[Level_Fine][nb].M[i_fine+1][j_fine+1] *
	     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine+1].A * 
	     Local_SolnBlk[Level_Fine][nb].Sp[i_fine+1][j_fine+1]) /
	    (Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].A * 
	     Local_SolnBlk[Level_Coarse][nb].Sp[i_coarse][j_coarse]);	

	} // endfor - cells

    } // endif - used
  } // endfor - blocks

}



/**********************************************************************
 * Routine: Apply_ICs                                                 *
 *                                                                    *
 * This routing applies the initial conditions on all coarse grid     *
 * levels.                                                            *
 *                                                                    *
 **********************************************************************/
template <> void FAS_Multigrid2D_Solver<Rte2D_State,
                                        Rte2D_Quad_Block,
                                        Rte2D_Input_Parameters>::
Apply_ICs(const int &level) 
{ 
  //
  // conserved solution state
  //
  // specified ics
  ICs(Local_SolnBlk[level],
      List_of_Local_Solution_Blocks[level],
      *IP);

  // if restart, restrict fine solution
  if (IP->i_ICs == IC_RESTART) {
    Restrict_Solution_Blocks(level-1);
    // Update_Primitive_Variables(level); // <- don't need it
  }

  // If this is a discretely specified medium field,
  // then we will have to restrict it as well.
  if ( Local_SolnBlk[level]->Medium_Field_Type == MEDIUM2D_FIELD_DISCRETE )
    Restrict_NonSolution_Blocks(Local_SolnBlk, 
				List_of_Local_Solution_Blocks, 
				level-1);
  
  // apply coarse mesh BCs
  if (IP->Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions)
    Restrict_Boundary_Ref_States(level-1);

}



#endif // _RTE2D_MULTIGRID_INCLUDED 
