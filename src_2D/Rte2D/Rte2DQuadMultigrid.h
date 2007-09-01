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
 * Routine:allocate                                                   *
 *                                                                    *
 * This routine performs the memory allocation and initialization for *
 * all grid levels of the FAS multigrid solution class.               *
 *                                                                    *
 **********************************************************************/
template <> int FAS_Multigrid2D_Solver<Rte2D_State,
                                       Rte2D_Quad_Block,
                                       Rte2D_Input_Parameters>::
allocate(Rte2D_Quad_Block *FinestBlks,
	 QuadTreeBlock_DataStructure *FinestQuadTree,
	 AdaptiveBlockResourceList *FinestGlobalList,
	 AdaptiveBlock2D_List *FinestLocalList,
	 Rte2D_Input_Parameters *ip) {

  int error_flag;

  // Point the input parameters to the given input parameters.
  IP = ip;

  // Point the quadtree to the input quadtree.
  QuadTree = FinestQuadTree;

  // Point the global solution block list to the input list.
  List_of_Global_Solution_Blocks = FinestGlobalList;

  // Create list of local solution blocks for each multigrid level.
  List_of_Local_Solution_Blocks = new AdaptiveBlock2D_List[IP->Multigrid_IP.Levels];

  // Point the list of local solution blocks for the finest level to
  // the input list of local solution blocks.
  List_of_Local_Solution_Blocks[FINEST_LEVEL] = *FinestLocalList;
  // Allocate memory for the coarse grid list of local solution blocks
  // and set the CPU number.
  for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {
    List_of_Local_Solution_Blocks[level].allocate(IP->Number_of_Blocks_Per_Processor);
    List_of_Local_Solution_Blocks[level].ThisCPU = List_of_Local_Solution_Blocks[FINEST_LEVEL].ThisCPU;
  }

  // Allocate memory for the local solution blocks on each level.
  Local_SolnBlk = new Rte2D_Quad_Block*[IP->Multigrid_IP.Levels];
  // Point the local solution block for the finest level to the input
  // local solution block.
  Local_SolnBlk[FINEST_LEVEL] = FinestBlks;
  // Allocate memory for the coarse grid local solution blocks.
  for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {
    Local_SolnBlk[level] = new Rte2D_Quad_Block[IP->Number_of_Blocks_Per_Processor];
  }

  // Allocate memory for the FAS multigrid solution blocks on each level.
  MG_SolnBlk = new FAS_Multigrid_Quad_Block<Rte2D_State>*[IP->Multigrid_IP.Levels];
  // Allocate memory for the coarse grid FAS multigrid solution blocks.
  for (int level = 0; level < IP->Multigrid_IP.Levels; level++) {
    MG_SolnBlk[level] = new FAS_Multigrid_Quad_Block<Rte2D_State>[IP->Number_of_Blocks_Per_Processor];
  }

  // Allocate memory for the DTS multigrid solution blocks on each level.
  if (IP->i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
    DTS_SolnBlk = new DTS_Multigrid_Quad_Block<Rte2D_State>*[IP->Multigrid_IP.Levels];
    // Allocate memory for the coarse grid DTS multigrid solution blocks.
    for (int level = 0; level < IP->Multigrid_IP.Levels; level++) {
      DTS_SolnBlk[level] = new DTS_Multigrid_Quad_Block<Rte2D_State>[IP->Number_of_Blocks_Per_Processor];
    }
  }

  // Allocate memory and set data for all coarse mesh variables on all
  // blocks on this processor.
  for (int nb = 0; nb < IP->Number_of_Blocks_Per_Processor; nb++) {
    if (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].used == ADAPTIVEBLOCK2D_USED) {

      // Ensure that the number of cells in each direction is even.
      if (Local_SolnBlk[FINEST_LEVEL][nb].NCi % 2 != 0) {
	cout << "\nFASMultigrid2D Error: block #" << nb << " on processor " 
             << List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.cpu 
             << ", on level 0 has " << Local_SolnBlk[FINEST_LEVEL][nb].NCi 
             << "cells in the x-direction, which is odd; cannot coarsen any further." << endl;
	return 1101;
      }
      if (Local_SolnBlk[FINEST_LEVEL][nb].NCj % 2 != 0) {
	cout << "\nFASMultigrid2D Error: block #" << nb << " on processor " 
             << List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.cpu 
             << ", on level 0 has " << Local_SolnBlk[FINEST_LEVEL][nb].NCj 
             << " cells in the y-direction, which is odd; cannot coarsen any further." << endl;
	return 1102;
      }

      // Allocate memory for the forcing term and the uo storage on the
      // finest level.
      MG_SolnBlk[FINEST_LEVEL][nb].allocate(Local_SolnBlk[FINEST_LEVEL][nb].NCi,
					    Local_SolnBlk[FINEST_LEVEL][nb].NCj);

      // Allocate memory for the DTS solution block finest level.
      if (IP->i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
	DTS_SolnBlk[FINEST_LEVEL][nb].allocate(Local_SolnBlk[FINEST_LEVEL][nb].NCi,
					       Local_SolnBlk[FINEST_LEVEL][nb].NCj);
      }

      // Allocate memory and set data for the coarse mesh levels.
      for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {

	// Copy the list of local solution block neighbour information
	// and calculate the coarse grid resolution.
	List_of_Local_Solution_Blocks[level].Block[nb] =
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb];
	List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.i = 
	  int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.dimen.i/pow(2.0,double(level)));
	List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.j = 
	  int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.dimen.j/pow(2.0,double(level)));
	
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nS; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoS[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoS[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoS[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoS[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nN; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoN[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoN[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoN[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoN[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nE; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoE[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoE[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoE[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoE[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nW; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoW[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoW[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoW[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoW[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nSE; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoSE[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoSE[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoSE[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoSE[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nSW; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoSW[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoSW[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoSW[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoSW[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nNE; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoNE[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoNE[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoNE[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoNE[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[nb].nNW; n++) {
	  List_of_Local_Solution_Blocks[level].Block[nb].infoNW[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoNW[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[nb].infoNW[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].infoNW[n].dimen.j/pow(2.0,double(level)));
	}

	// Set-up the local block list for each level.
	List_of_Local_Solution_Blocks[level].Block[nb].used = ADAPTIVEBLOCK2D_USED;
	List_of_Local_Solution_Blocks[level].Block[nb].gblknum = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].gblknum;
	List_of_Local_Solution_Blocks[level].Block[nb].info.cpu = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.cpu;
	List_of_Local_Solution_Blocks[level].Block[nb].info.blknum = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.blknum;
	List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.i = 
	  int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.dimen.i/pow(2.0,double(level)));
	List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.j = 
	  int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.dimen.j/pow(2.0,double(level)));

	// If the current level is not the coarsest level then ensure
	// that the number of cells is even.
	if (level != IP->Multigrid_IP.Levels-1) {
	  if (List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.i % 2 != 0) {
	    cout << "\nFASMultigrid2D Error: block #" << nb << " on processor " 
                 << List_of_Local_Solution_Blocks[level].Block[nb].info.cpu 
                 << ", on level " << level << " has " 
                 << List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.i 
                 << " cells in the x-direction, which is odd; cannot coarsen any further." << endl;
	    return 1103;
	  }
	}
	if (level != IP->Multigrid_IP.Levels-1) {
	  if (List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.j % 2 != 0) {
	    cout << "\nFASMultigrid2D Error: block #" << nb << " on processor " 
                 << List_of_Local_Solution_Blocks[level].Block[nb].info.cpu << ", on level " 
                 << level << " has " << List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.j 
                 << " cells in the y-direction, which is odd; cannot coarsen any further." << endl;
	    return 1104;
	  }
	}
	
	List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.ghost = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.dimen.ghost;
	List_of_Local_Solution_Blocks[level].Block[nb].info.sector = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.sector;
	List_of_Local_Solution_Blocks[level].Block[nb].info.level = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].info.level;
	
	// Allocate the coarse grid local solution block.
	Local_SolnBlk[level][nb].allocate(List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.i,
					  List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.j,
					  List_of_Local_Solution_Blocks[level].Block[nb].info.dimen.ghost);

	// Create the coarse grid mesh.
	Half_Mesh_Resolution(Local_SolnBlk[level][nb].Grid,
			     Local_SolnBlk[level-1][nb].Grid);

	/***********************************************************************
	 *************************** RTE SPECIFIC ******************************/
	// added call to compute 2D to Quasi-3D scaling param for each grid
	Local_SolnBlk[level][nb].ScaleGridTo3D(ip->Axisymmetric);
	/*************************************************************************
	 *************************************************************************/

	// Allocate the coarse grid FAS multigrid solution block.
	MG_SolnBlk[level][nb].allocate(Local_SolnBlk[level][nb].NCi,
				       Local_SolnBlk[level][nb].NCj);

	// Allocate the coarse grid DTS multigrid solution block.
	if (IP->i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
	  DTS_SolnBlk[level][nb].allocate(Local_SolnBlk[level][nb].NCi,
					  Local_SolnBlk[level][nb].NCj);
	}

	// Allocate memory for the message passing buffers used to send solution
	// information between neighbouring blocks for the coarse grid.
	Allocate_Message_Buffers(List_of_Local_Solution_Blocks[level],
				 Local_SolnBlk[level][nb].NumVar()+NUM_COMP_VECTOR2D);

      }

    }
  }

  // Perform distance to wall calcuation on coarse grids as required.
  error_flag = Determine_Wall_Distance_on_Coarse_Grids();
  if (error_flag) return error_flag;

  // Solution block allocation and assignment successful.
  return 0;

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





#endif // _RTE2D_MULTIGRID_INCLUDED 