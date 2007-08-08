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
#ifndef _RTE2D_MULTIGRID_SPECIALIZATION_INCLUDED 
#define _RTE2D_MULTIGRID_INCLUDED 

/* Include 2D Rte equation quadrilateral mesh solution header file. */

#ifndef _RTE2D_QUAD_INCLUDED
#include "Rte2DQuad.h"
#endif // _RTE2D_QUAD_INCLUDED

/* Include 2D Rte equation input header file. */

#ifndef _RTE2D_INPUT_INCLUDED
#include "Rte2D/Rte2DInput.h"
#endif // _RTE2D_INPUT_INCLUDED

/* Include FAS Multigrid header file */

#ifndef _FASMULTIGRID2D_INCLUDED
#include "../FASMultigrid2D.h"
#endif // _FASMULTIGRID2D_INCLUDED



/******************************************************************
 * Routine:  FAS_Multigrid2D_Solver class memory allocation       *
 *                                                                *
 * constructs class using data inside argument Input_Parameters   *
 *                                                                *
 ******************************************************************/
template <> void FAS_Multigrid2D_Solver<Rte2D_State,
                                        Rte2D_Quad_Block,
                                        Rte2D_Input_Parameters>::
allocate(int &setbatchflag,
	 Rte2D_Quad_Block *FinestBlks,
	 QuadTreeBlock_DataStructure *FinestQuadTree,
	 AdaptiveBlock2D_List *FinestLocalList,
	 Rte2D_Input_Parameters *IP,
	 int num_other_vars) {

  int level;
  
  Input_Parameters = IP;
  error_flag = 0;
  batch_flag = setbatchflag;
  residual_l1_norm = 0.00;
  residual_l2_norm = 0.00;
  residual_max_norm = 0.00;
  Num_Other_Vars = num_other_vars;
  flag_for_First_Time_Restriction = true;

  /* Pointer to QuadTree */
  
  QuadTree = FinestQuadTree;
  
  List_of_Local_Solution_Blocks = new AdaptiveBlock2D_List[Input_Parameters->Multigrid_IP.Levels];
    
  List_of_Local_Solution_Blocks[FINEST_LEVEL] = *FinestLocalList;  // Make a copy of the finest local list
  for (level = 1; level < Input_Parameters->Multigrid_IP.Levels; level++) {
    List_of_Local_Solution_Blocks[level].allocate(Input_Parameters->Number_of_Blocks_Per_Processor);
    List_of_Local_Solution_Blocks[level].ThisCPU = List_of_Local_Solution_Blocks[FINEST_LEVEL].ThisCPU;
  } /* end for */

  /* Allocate storage for solution blocks on each level */
  Local_SolnBlk = new Rte2D_Quad_Block*[Input_Parameters->Multigrid_IP.Levels];
  MG_Blk = new Additional_Storage_for_Multigrid<Rte2D_State>*[Input_Parameters->Multigrid_IP.Levels];
  
  Local_SolnBlk[FINEST_LEVEL] = FinestBlks;  // Point finest level pointer to FinestBlks
  /* Create coarser level blocks */
  for (level = 1; level < Input_Parameters->Multigrid_IP.Levels; level++) {
    Local_SolnBlk[level] = new Rte2D_Quad_Block[Input_Parameters->Number_of_Blocks_Per_Processor];
  }
  
  for (level = 0; level < Input_Parameters->Multigrid_IP.Levels; level++) {
    MG_Blk[level] = new Additional_Storage_for_Multigrid<Rte2D_State>[Input_Parameters->Number_of_Blocks_Per_Processor];
  } /* end for */

  /* Loop over all blocks on this processor */
  for (int b = 0; b < Input_Parameters->Number_of_Blocks_Per_Processor; ++b) {
    if (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].used == ADAPTIVEBLOCK2D_USED) {
      /* Check if the number of cells in each direction is even */
      if (Local_SolnBlk[FINEST_LEVEL][b].NCi % 2 != 0) {
	error_flag = 1;
	cout << "\nFASMultigrid2D Error: block #" << b << " on processor " 
             << List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].info.cpu 
             << ", on level 0 has " << Local_SolnBlk[FINEST_LEVEL][b].NCi 
             << "cells in the x-direction, which is odd; cannot coarsen any further." << endl;
	//return (error_flag);
      }
      if (Local_SolnBlk[FINEST_LEVEL][b].NCj % 2 != 0) {
	error_flag = 1;
	cout << "\nFASMultigrid2D Error: block #" << b << " on processor " 
             << List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].info.cpu 
             << ", on level 0 has " << Local_SolnBlk[FINEST_LEVEL][b].NCj 
             << " cells in the y-direction, which is odd; cannot coarsen any further." << endl;
	//return (error_flag);
      }

      /* allocate memory for forcing term and uo storage, finest level*/
      MG_Blk[FINEST_LEVEL][b].allocate(Local_SolnBlk[FINEST_LEVEL][b].NCi,
				       Local_SolnBlk[FINEST_LEVEL][b].NCj);

      /* allocate memory for coarser levels */
      for (level = 1; level < Input_Parameters->Multigrid_IP.Levels; level++) {

	/* Copy Local_List neighbour information */
	/* and calculate coarser grid resolution */
	
	List_of_Local_Solution_Blocks[level].Block[b] =
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b];
	List_of_Local_Solution_Blocks[level].Block[b].info.dimen.i = 
	  int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].info.dimen.i/pow(2.0,double(level)));
	List_of_Local_Solution_Blocks[level].Block[b].info.dimen.j = 
	  int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].info.dimen.j/pow(2.0,double(level)));
	
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[b].nS; n++) {
	  List_of_Local_Solution_Blocks[level].Block[b].infoS[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].infoS[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[b].infoS[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].infoS[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[b].nN; n++) {
	  List_of_Local_Solution_Blocks[level].Block[b].infoN[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].infoN[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[b].infoN[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].infoN[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[b].nE; n++) {
	  List_of_Local_Solution_Blocks[level].Block[b].infoE[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].infoE[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[b].infoE[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].infoE[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[b].nW; n++) {
	  List_of_Local_Solution_Blocks[level].Block[b].infoW[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].infoW[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[b].infoW[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].infoW[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[b].nSE; n++) {
	  List_of_Local_Solution_Blocks[level].Block[b].infoSE[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].infoSE[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[b].infoSE[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].infoSE[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[b].nSW; n++) {
	  List_of_Local_Solution_Blocks[level].Block[b].infoSW[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].infoSW[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[b].infoSW[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].infoSW[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[b].nNE; n++) {
	  List_of_Local_Solution_Blocks[level].Block[b].infoNE[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].infoNE[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[b].infoNE[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].infoNE[n].dimen.j/pow(2.0,double(level)));
	}
	for (int n = 0; n < List_of_Local_Solution_Blocks[level].Block[b].nNW; n++) {
	  List_of_Local_Solution_Blocks[level].Block[b].infoNW[n].dimen.i = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].infoNW[n].dimen.i/pow(2.0,double(level)));
	  List_of_Local_Solution_Blocks[level].Block[b].infoNW[n].dimen.j = 
	    int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].infoNW[n].dimen.j/pow(2.0,double(level)));
	}

	/* First setup local block list for each level */
	List_of_Local_Solution_Blocks[level].Block[b].used = ADAPTIVEBLOCK2D_USED;
	List_of_Local_Solution_Blocks[level].Block[b].gblknum = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].gblknum;
	List_of_Local_Solution_Blocks[level].Block[b].info.cpu = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].info.cpu;
	List_of_Local_Solution_Blocks[level].Block[b].info.blknum = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].info.blknum;
	List_of_Local_Solution_Blocks[level].Block[b].info.dimen.i = 
	  int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].info.dimen.i/pow(2.0,double(level)));

	/* If this is not the coarsest level, check if the # of cells is even */
	if (level != Input_Parameters->Multigrid_IP.Levels-1) {
	  if (List_of_Local_Solution_Blocks[level].Block[b].info.dimen.i % 2 != 0) {
	    error_flag = 1;
	    cout << "\nFASMultigrid2D Error: block #" << b << " on processor " 
                 << List_of_Local_Solution_Blocks[level].Block[b].info.cpu 
                 << ", on level " << level << " has " 
                 << List_of_Local_Solution_Blocks[level].Block[b].info.dimen.i 
                 << " cells in the x-direction, which is odd; cannot coarsen any further." << endl;
	    //return (error_flag);
	  }
	}
	
	List_of_Local_Solution_Blocks[level].Block[b].info.dimen.j = 
	  int (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].info.dimen.j/pow(2.0,double(level)));
	
	/* If this is not the coarsest level, check if the # of cells is even */
	if (level != Input_Parameters->Multigrid_IP.Levels-1) {
	  if (List_of_Local_Solution_Blocks[level].Block[b].info.dimen.j % 2 != 0) {
	    error_flag = 1;
	    cout << "\nFASMultigrid2D Error: block #" << b << " on processor " 
                 << List_of_Local_Solution_Blocks[level].Block[b].info.cpu << ", on level " 
                 << level << " has " << List_of_Local_Solution_Blocks[level].Block[b].info.dimen.j 
                 << " cells in the y-direction, which is odd; cannot coarsen any further." << endl;
	    //return (error_flag);
	  }
	}
	
	List_of_Local_Solution_Blocks[level].Block[b].info.dimen.ghost = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].info.dimen.ghost;
	List_of_Local_Solution_Blocks[level].Block[b].info.sector = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].info.sector;
	List_of_Local_Solution_Blocks[level].Block[b].info.level = 
	  List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[b].info.level;
	
	/* allocate coarser soln blocks and MG blocks */
	Local_SolnBlk[level][b].allocate(
	  List_of_Local_Solution_Blocks[level].Block[b].info.dimen.i, 
	  List_of_Local_Solution_Blocks[level].Block[b].info.dimen.j,
	  List_of_Local_Solution_Blocks[level].Block[b].info.dimen.ghost);
	MG_Blk[level][b].allocate(Local_SolnBlk[level][b].NCi,
				  Local_SolnBlk[level][b].NCj);
	/* Coarsen grid and copy to coarse blocks */
	Half_Mesh_Resolution(Local_SolnBlk[level][b].Grid,
			     Local_SolnBlk[level-1][b].Grid);

	/* Rte2D Specific - Need to compute 2D to Quasi-3D scaling param */
	Local_SolnBlk[level][b].ScaleGridTo3D(Input_Parameters->Axisymmetric);
	
	/* Allocates memory for all message passing buffers used  
	   to send solution information between neighbouring solnblks. */
	
	Allocate_Message_Buffers(List_of_Local_Solution_Blocks[level],
				 Local_SolnBlk[level][b].NumVar()+Num_Other_Vars);
      } /* endfor */
    }
  }
  /* Solution block allocation and assignment complete. */  
}




/********************************************************
 * Routine: Restrict_Solution_Blocks (for Multigrid)    *
 *                                                      *
 * Restrict solution from Level_Fine to                 *
 * Level_Coarse for all blocks on local block list      *
 * **Soln_ptr (overwrites any solution                  *
 * present on Level_Coarse)                             *
 *                                                      *
 ********************************************************/
template <> inline void FAS_Multigrid2D_Solver<Rte2D_State,
			  		       Rte2D_Quad_Block,
					       Rte2D_Input_Parameters>::
Restrict_Solution_Blocks(const int Level_Fine) {

  int i_fine,j_fine;
  int Level_Coarse = Level_Fine + 1;
  
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
	  Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse] = 
	    (Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine] *
	     Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].A*Local_SolnBlk[Level_Fine][b].Sp[i_fine][j_fine] +
	     Local_SolnBlk[Level_Fine][b].U[i_fine+1][j_fine] *
	     Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][j_fine].A*Local_SolnBlk[Level_Fine][b].Sp[i_fine][j_fine] +
	     Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine+1] *
	     Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine+1].A*Local_SolnBlk[Level_Fine][b].Sp[i_fine][j_fine] +
	     Local_SolnBlk[Level_Fine][b].U[i_fine+1][j_fine+1] *
	     Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][j_fine+1].A*Local_SolnBlk[Level_Fine][b].Sp[i_fine][j_fine]) /
	    (Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].A*Local_SolnBlk[Level_Coarse][b].Sp[i_coarse][j_coarse]);	
	       
	  //------------------------ Rte2D Specific -------------------------//

	  // For the RTE, we also need to restrict and prolong the
	  // coefficients and the blackbody intensity
	  
// 	  Restrict_NonSol( Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 			   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].A, 
// 			   Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine],
// 			   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].A,
// 			   Local_SolnBlk[Level_Fine][b].U[i_fine+1][j_fine],
// 			   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][j_fine].A,
// 			   Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine+1],
// 			   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine+1].A,
// 			   Local_SolnBlk[Level_Fine][b].U[i_fine+1][j_fine+1],
// 			   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][j_fine+1].A );
	  
	  //---------------------- End Rte2D Specific -----------------------//

	} /* endfor */
      } /* endfor */
    } /* end if */     
  } /* endfor */
}

/********************************************************
 * Routine: CFL_Multigrid                               *
 *                                                      *
 * Calculate Time Steps for Level_Coarse for all blocks *
 *                                                      *
 ********************************************************/
template <> inline void FAS_Multigrid2D_Solver<Rte2D_State,
			  		       Rte2D_Quad_Block,
					       Rte2D_Input_Parameters>::
CFL_Multigrid(const int Level_Coarse) {

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
	  
	  A_coarse = Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].A*
	    Local_SolnBlk[Level_Coarse][b].Sp[i_coarse][j_coarse];

	  dt_NE = sqrt(A_coarse/(Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][j_fine+1].A*
				 Local_SolnBlk[Level_Fine][b].Sp[i_fine+1][j_fine+1]));
	  dt_NE *= Local_SolnBlk[Level_Fine][b].dt[i_fine+1][j_fine+1];
	  dt_SE = sqrt(A_coarse/(Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][j_fine].A*
				 Local_SolnBlk[Level_Fine][b].Sp[i_fine+1][j_fine]));
	  dt_SE *= Local_SolnBlk[Level_Fine][b].dt[i_fine+1][j_fine];
	  dt_NW = sqrt(A_coarse/(Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine+1].A*
				 Local_SolnBlk[Level_Fine][b].Sp[i_fine][j_fine+1]));
	  dt_NW *= Local_SolnBlk[Level_Fine][b].dt[i_fine][j_fine+1];
	  dt_SW = sqrt(A_coarse/(Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].A*
				 Local_SolnBlk[Level_Fine][b].Sp[i_fine][j_fine]));
	  dt_SW *= Local_SolnBlk[Level_Fine][b].dt[i_fine][j_fine];

	  Local_SolnBlk[Level_Coarse][b].dt[i_coarse][j_coarse] = 
                  min(Local_SolnBlk[Level_Coarse][b].dt[i_coarse][j_coarse],
                      min(dt_NE,min(dt_SE,min(dt_NW,dt_SW))));

	} /* endfor */
      } /* endfor */
    } /* end if */     
  } /* endfor */
}


/********************************************************
 * Routine: Restrict_Residuals (for Multigrid)          *
 *                                                      *
 * Restrict residual dUdt from Level_Fine to            *
 * Level_Coarse (stored in P) for all blocks on local   *
 * block list **Soln_ptr (overwrites any value stored   *
 * in P on Level_Coarse)                                *
 *                                                      *
 ********************************************************/
template <> inline void FAS_Multigrid2D_Solver<Rte2D_State,
			  		       Rte2D_Quad_Block,
					       Rte2D_Input_Parameters>::
Restrict_Residuals(const int Level_Fine) {

  int i_fine,j_fine,Level_Coarse;
  Level_Coarse = Level_Fine + 1;
  
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

	  /* Restrict Residuals */
	  
	  MG_Blk[Level_Coarse][b].P[i_coarse][j_coarse] = 
	    (Local_SolnBlk[Level_Fine][b].dUdt[i_fine][j_fine][0] *
	     Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].A *
	     Local_SolnBlk[Level_Fine][b].Sp[i_fine][j_fine] +
	     Local_SolnBlk[Level_Fine][b].dUdt[i_fine+1][j_fine][0] *
	     Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][j_fine].A *
	     Local_SolnBlk[Level_Fine][b].Sp[i_fine+1][j_fine] +
	     Local_SolnBlk[Level_Fine][b].dUdt[i_fine][j_fine+1][0] *
	     Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine+1].A *
	     Local_SolnBlk[Level_Fine][b].Sp[i_fine][j_fine+1] +
	     Local_SolnBlk[Level_Fine][b].dUdt[i_fine+1][j_fine+1][0] *
	     Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][j_fine+1].A *
	     Local_SolnBlk[Level_Fine][b].Sp[i_fine+1][j_fine+1] ) /
	    (Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].A * 
	     Local_SolnBlk[Level_Coarse][b].Sp[i_coarse][j_coarse] );
			       
	} /* endfor */
      } /* endfor */
    } /* end if */     
  } /* endfor */
}



/*********************************************************
 * Routine: Restrict_Boundary_Ref_States (for Multigrid) *
 *                                                       *
 * Restrict boundary_ref_states UoN/S/E/W from Level_Fine*
 * to Level_Coarse for all blocks (overwrites any value  *
 * stored in UoN/S/E/W on Level_Coarse). The restriction *
 * operator used is area weighted average                *
 *                                                       *
 *********************************************************/
template <> inline void FAS_Multigrid2D_Solver<Rte2D_State,
			  		       Rte2D_Quad_Block,
					       Rte2D_Input_Parameters>::
Restrict_Boundary_Ref_States(const int Level_Fine) {

  int i_fine,j_fine;
  int Level_Coarse = Level_Fine + 1;
  
  // loop through used blocks only
  for (int b = 0 ; b <= List_of_Local_Solution_Blocks[Level_Fine].Nblk-1 ; ++b ) {
    if (List_of_Local_Solution_Blocks[Level_Fine].Block[b].used == ADAPTIVEBLOCK2D_USED) { 

      /* Get number of ghost cells */
      int Nghost = Local_SolnBlk[Level_Fine][b].Nghost;

      /* Get lower and upper indices */
      int JCl = Local_SolnBlk[Level_Fine][b].JCl;
      int JCu = Local_SolnBlk[Level_Fine][b].JCu;

      /* Loop through i-direction cells */
      for (int i_coarse = Local_SolnBlk[Level_Coarse][b].ICl; i_coarse <= Local_SolnBlk[Level_Coarse][b].ICu; i_coarse++) {
	i_fine = 2*(i_coarse-Nghost)+Nghost;

	Local_SolnBlk[Level_Coarse][b].UoN[i_coarse] = 
	  (Local_SolnBlk[Level_Fine][b].UoN[i_fine] *
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][JCu].A*Local_SolnBlk[Level_Fine][b].Sp[i_fine][JCu] +
	   Local_SolnBlk[Level_Fine][b].UoN[i_fine+1] *
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][JCu].A*Local_SolnBlk[Level_Fine][b].Sp[i_fine+1][JCu]) /
	  (Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][JCu].A*Local_SolnBlk[Level_Fine][b].Sp[i_fine][JCu] +
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][JCu].A*Local_SolnBlk[Level_Fine][b].Sp[i_fine+1][JCu]);

	Local_SolnBlk[Level_Coarse][b].UoS[i_coarse] = 
	  (Local_SolnBlk[Level_Fine][b].UoS[i_fine] *
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][JCl].A*Local_SolnBlk[Level_Fine][b].Sp[i_fine][JCl] +
	   Local_SolnBlk[Level_Fine][b].UoS[i_fine+1] *
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][JCl].A*Local_SolnBlk[Level_Fine][b].Sp[i_fine+1][JCl]) /
	  (Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][JCl].A*Local_SolnBlk[Level_Fine][b].Sp[i_fine][JCl] +
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][JCl].A*Local_SolnBlk[Level_Fine][b].Sp[i_fine+1][JCl]);

	  //------------------------ Rte2D Specific -------------------------//

	  // For the RTE, we also need to restrict and prolong the
	  // coefficients and the blackbody intensity

// 	Restrict_NonSol_Boundary_Ref_States( 
// 		    Local_SolnBlk[Level_Coarse][b].UoN[i_coarse],
// 		    Local_SolnBlk[Level_Fine][b].UoN[i_fine],
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][JCu].A,
// 		    Local_SolnBlk[Level_Fine][b].UoN[i_fine+1],
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][JCu].A );
   
// 	Restrict_NonSol_Boundary_Ref_States( 
// 		    Local_SolnBlk[Level_Coarse][b].UoS[i_coarse],
// 		    Local_SolnBlk[Level_Fine][b].UoS[i_fine],
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][JCl].A,
// 		    Local_SolnBlk[Level_Fine][b].UoS[i_fine+1],
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][JCl].A );
	  
	  //---------------------- End Rte2D Specific -----------------------//


      } /* endfor */

      /* Get lower and upper indices */
      int ICl = Local_SolnBlk[Level_Fine][b].ICl;
      int ICu = Local_SolnBlk[Level_Fine][b].ICu;

      /* Loop through j-direction cells */
      for (int j_coarse = Local_SolnBlk[Level_Coarse][b].JCl; j_coarse <= Local_SolnBlk[Level_Coarse][b].JCu; j_coarse++) {
	j_fine = 2*(j_coarse-Nghost)+Nghost;

	Local_SolnBlk[Level_Coarse][b].UoW[j_coarse] = 
	  (Local_SolnBlk[Level_Fine][b].UoW[j_fine] *
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[ICl][j_fine].A*Local_SolnBlk[Level_Fine][b].Sp[ICl][j_fine] +
	   Local_SolnBlk[Level_Fine][b].UoW[j_fine+1] *
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[ICl][j_fine+1].A*Local_SolnBlk[Level_Fine][b].Sp[ICl][j_fine+1]) /
	  (Local_SolnBlk[Level_Fine][b].Grid.Cell[ICl][j_fine].A*Local_SolnBlk[Level_Fine][b].Sp[ICl][j_fine] +
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[ICl][j_fine+1].A*Local_SolnBlk[Level_Fine][b].Sp[ICl][j_fine+1]);

	Local_SolnBlk[Level_Coarse][b].UoE[j_coarse] = 
	  (Local_SolnBlk[Level_Fine][b].UoE[j_fine] *
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[ICu][j_fine].A*Local_SolnBlk[Level_Fine][b].Sp[ICu][j_fine] +
	   Local_SolnBlk[Level_Fine][b].UoE[j_fine+1] *
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[ICu][j_fine+1].A*Local_SolnBlk[Level_Fine][b].Sp[ICu][j_fine+1]) /
	  (Local_SolnBlk[Level_Fine][b].Grid.Cell[ICu][j_fine].A*Local_SolnBlk[Level_Fine][b].Sp[ICu][j_fine] +
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[ICu][j_fine+1].A*Local_SolnBlk[Level_Fine][b].Sp[ICu][j_fine+1]);

	  //------------------------ Rte2D Specific -------------------------//

	  // For the RTE, we also need to restrict and prolong the
	  // coefficients and the blackbody intensity

// 	Restrict_NonSol_Boundary_Ref_States( 
// 		    Local_SolnBlk[Level_Coarse][b].UoW[j_coarse],
// 		    Local_SolnBlk[Level_Fine][b].UoW[j_fine],
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[ICl][j_fine].A,
// 		    Local_SolnBlk[Level_Fine][b].UoW[j_fine+1],
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[ICl][j_fine+1].A );
   
// 	Restrict_NonSol_Boundary_Ref_States( 
// 		    Local_SolnBlk[Level_Coarse][b].UoE[j_coarse],
// 		    Local_SolnBlk[Level_Fine][b].UoE[j_fine],
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[ICu][j_fine].A,
// 		    Local_SolnBlk[Level_Fine][b].UoE[j_fine+1],
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[ICu][j_fine+1].A );
	  
	  //---------------------- End Rte2D Specific -----------------------//



      } /* endfor */
      
    } /* end if */     
  } /* endfor */
}


/********************************************************************
 * Routine: Prolong_Solution_Blocks (for Multigrid)                 *
 *                                                                  *
 * Prolong solution states stored in U[][] from Level_Coarse        *
 * to Level_Fine for all blocks on local block list **Soln_ptr      *
 *                                                                  *
 ********************************************************************/
template <> inline int FAS_Multigrid2D_Solver<Rte2D_State,
				              Rte2D_Quad_Block,
				              Rte2D_Input_Parameters>::
Prolong_Solution_Blocks(const int Level_Coarse) {

  int i_fine,j_fine;
  Rte2D_State solution;
  int Level_Fine = Level_Coarse - 1;
  bool injection_for_this_coarse_cell;
  injection_for_this_coarse_cell = false;

  // loop through used blocks only
  for (int b = 0 ; b <= List_of_Local_Solution_Blocks[Level_Coarse].Nblk-1 ; ++b ) {
    if (List_of_Local_Solution_Blocks[Level_Coarse].Block[b].used == ADAPTIVEBLOCK2D_USED) { 

      /* Get number of ghost cells */
      int Nghost = Local_SolnBlk[Level_Coarse][b].Nghost;

      /* Loop through coarse cells */
      for (int i_coarse = Local_SolnBlk[Level_Coarse][b].ICl; i_coarse <= Local_SolnBlk[Level_Coarse][b].ICu; i_coarse++) {
	for (int j_coarse = Local_SolnBlk[Level_Coarse][b].JCl; j_coarse <= Local_SolnBlk[Level_Coarse][b].JCu; j_coarse++) {

	  if (Input_Parameters->Multigrid_IP.PROLONGATE_USING_INJECTION) injection_for_this_coarse_cell = true;
	  
	  restart: ;
	
	  /*******************************************************************/
	  /*** Calculate cell-centred solution at SW Corner of coarse cell ***/
	  /*******************************************************************/
	  i_fine = 2*(i_coarse-Nghost)+Nghost; // SW Corner i index (fine)
	  j_fine = 2*(j_coarse-Nghost)+Nghost; // SW Corner j index (fine)
	  
	  if (!injection_for_this_coarse_cell) {
	    error_flag = Bilinear_Interpolation(
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1],
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
	       Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc,
	       solution);
// 	    if (error_flag == 0) {
// 	      error_flag = Prolong_NonSol(
// 		  Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1],
// 		  Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
// 		  Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
// 		  Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 		  Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		  Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		  Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
// 		  Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 		  Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc,
// 		  solution );
// 	    }
	    
	    /* If fine grid cell-centre falls outside of quad in SW dir
	     formed by coarse grid cell-centres, try NW one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
// 	      if (error_flag == 0) {
// 		error_flag = Prolong_NonSol(
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		    solution);
// 	      }
	    }
	    /* If fine grid cell-centre falls outside of quad in NW dir
	       formed by coarse grid cell-centres, try SE one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
// 	      if (error_flag == 0) {
// 		error_flag = Prolong_NonSol(
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		    solution);
// 	      }
	    }
	    /* If fine grid cell-centre falls outside of quad in SE dir
	       formed by coarse grid cell-centres, try NE one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
// 	      if (error_flag == 0) {
// 		error_flag = Prolong_NonSol(
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		    solution);
// 	      }
	    }
	    if (error_flag != 0) {
	      injection_for_this_coarse_cell = true;
	      goto restart;
	    }
	  } else {
	    // use simple injection instead of bilinear interpolation
	    solution = Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse];
	    error_flag = 0;
	  }
	
	  Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine] = solution;
	
	  /**********************************************/
	  /*** Finished with SW fine-grid cell centre ***/
	  /**********************************************/
	  
	  /*******************************************************************/
	  /*** Calculate cell-centred solution at SE Corner of coarse cell ***/
	  /*******************************************************************/
	  
	  i_fine = 2*(i_coarse-Nghost)+Nghost+1; // SE Corner i index (fine)
	  j_fine = 2*(j_coarse-Nghost)+Nghost;   // SE Corner j index (fine)
	  
	  if (!injection_for_this_coarse_cell) {
	    error_flag = Bilinear_Interpolation(
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1],
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
	       Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
	       solution);
	    if (error_flag == 0) {
	      error_flag = Prolong_NonSol(
		  Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
		  Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		  Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		  Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		  Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
		  Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		  Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1],
		  Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
		  Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		  solution);
	    }
	    
	    /* If fine grid cell-centre falls outside of quad in SE dir
	       formed by coarse grid cell-centres, try NE one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
// 	      if (error_flag == 0) {
// 		error_flag = Prolong_NonSol(
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		    solution);
// 	      }
	    }
	    /* If fine grid cell-centre falls outside of quad in NE dir
	       formed by coarse grid cell-centres, try SW one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
// 	      if (error_flag == 0) {
// 		error_flag = Prolong_NonSol(
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		    solution);
// 	      }
	    }
	    /* If fine grid cell-centre falls outside of quad in SW dir
	       formed by coarse grid cell-centres, try NW one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
// 	      if (error_flag == 0) {
// 		error_flag = Prolong_NonSol(
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		    solution);
// 	      }
	    }

	    if (error_flag != 0) {
	      injection_for_this_coarse_cell = true;
	      goto restart;
	    }
	  } else {
	    // use simple injection instead of bilinear interpolation
	    solution = Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse];
	    error_flag = 0;
	  }
	  
	  Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine] = solution;
	  
	  /**********************************************/
	  /*** Finished with SE fine-grid cell centre ***/
	  /**********************************************/
	  
	  /*******************************************************************/
	  /*** Calculate cell-centred solution at NE Corner of coarse cell ***/
	  /*******************************************************************/
	  
	  i_fine = 2*(i_coarse-Nghost)+Nghost+1; // NE Corner i index (fine)
	  j_fine = 2*(j_coarse-Nghost)+Nghost+1; // NE Corner j index (fine)
	  
	  if (!injection_for_this_coarse_cell) {
	    error_flag = Bilinear_Interpolation(
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1],
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
	       Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
	       solution);
// 	    if (error_flag == 0) {
// 	      error_flag = Prolong_NonSol(
// 		  Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		  Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		  Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
// 		  Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 		  Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1],
// 		  Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
// 		  Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
// 		  Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 		  Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		  solution);
// 	    }
					  
	    /* If fine grid cell-centre falls outside of quad in NE dir
	       formed by coarse grid cell-centres, try SE one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
// 	      if (error_flag == 0) {
// 		error_flag = Prolong_NonSol(
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		    solution);
// 	      }
	    }
	    /* If fine grid cell-centre falls outside of quad in SE dir
	       formed by coarse grid cell-centres, try NW one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
// 	      if (error_flag == 0) {
// 		error_flag = Prolong_NonSol(
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		    solution);
// 	      }
	    }
	    /* If fine grid cell-centre falls outside of quad in NW dir
	       formed by coarse grid cell-centres, try SW one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
// 	      if (error_flag == 0) {
// 		error_flag = Prolong_NonSol(
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		    solution);
// 	      }
	    }
	    if (error_flag != 0) {
	      injection_for_this_coarse_cell = true;
	      goto restart;
	    }
	  } else {
	    // use simple injection instead of bilinear interpolation
	    solution = Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse];
	    error_flag = 0;
	  }
	
	  Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine] = solution;
	  
	  /**********************************************/
	  /*** Finished with NE fine-grid cell centre ***/
	  /**********************************************/
	  
	  /*******************************************************************/
	  /*** Calculate cell-centred solution at NW Corner of coarse cell ***/
	  /*******************************************************************/
	
	  i_fine = 2*(i_coarse-Nghost)+Nghost;   // NW Corner i index (fine)
	  j_fine = 2*(j_coarse-Nghost)+Nghost+1; // NW Corner j index (fine)
	  
	  if (!injection_for_this_coarse_cell) {
	    error_flag = Bilinear_Interpolation(
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1],
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
	       Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
	       solution);
// 	    if (error_flag == 0) {
// 	      error_flag = Prolong_NonSol(
// 		  Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
// 		  Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 		  Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1],
// 		  Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
// 		  Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
// 		  Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 		  Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		  Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		  Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		  solution);
// 	    }
	    
	    /* If fine grid cell-centre falls outside of quad in NW dir
	       formed by coarse grid cell-centres, try NE one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
// 	      if (error_flag == 0) {
// 		error_flag = Prolong_NonSol(
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		    solution);
// 	      }
	    }
	    /* If fine grid cell-centre falls outside of quad in NE dir
	       formed by coarse grid cell-centres, try SW one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
// 	      if (error_flag == 0) {
// 		error_flag = Prolong_NonSol(
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		    solution);
// 	      }
	    }
	    /* If fine grid cell-centre falls outside of quad in SW dir
	       formed by coarse grid cell-centres, try SE one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
// 	      if (error_flag == 0) {
// 		error_flag = Prolong_NonSol(
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		    solution);
// 	      }
	    }
	    if (error_flag != 0) {
	      injection_for_this_coarse_cell = true;
	      goto restart;
	    }
	  } else {
	    // use simple injection instead of bilinear interpolation
	    solution = Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse];
	    error_flag = 0;
	  }
	  
	  Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine] = solution;

	  /**********************************************/
	  /*** Finished with NW fine-grid cell centre ***/
	  /**********************************************/
	  
	  // reset injection flag
	  injection_for_this_coarse_cell = false;
	  
	} /* end if */
      } /* end for */
    } /* endfor */
  } /* endif */
  
  /* No Errors... prolongation completed */
  return(error_flag);
}


/********************************************************************
 * Routine: Prolong_and_Update_Solution_Blocks (for Multigrid)      *
 *                                                                  *
 * Prolong solution changes stored in U[][].u from Level_Coarse     *
 * to Level_Fine for all blocks on local block list **Soln_ptr and  *
 * add the changes to the solution U[][].u on Level_Fine blocks     *
 * (modifies any value present in U[][] on Level_Fine)              *
 *                                                                  *
 ********************************************************************/
template <> inline int FAS_Multigrid2D_Solver<Rte2D_State,
				              Rte2D_Quad_Block,
				              Rte2D_Input_Parameters>::
Prolong_and_Update_Solution_Blocks(const int Level_Coarse) {

  int i_fine,j_fine;
  Rte2D_State sol_change;
  int Level_Fine = Level_Coarse - 1;
  bool injection_for_this_coarse_cell;
  injection_for_this_coarse_cell = false;

  // loop through used blocks only
  for (int b = 0 ; b <= List_of_Local_Solution_Blocks[Level_Coarse].Nblk-1 ; ++b ) {
    if (List_of_Local_Solution_Blocks[Level_Coarse].Block[b].used == ADAPTIVEBLOCK2D_USED) { 

      /* Get number of ghost cells */
      int Nghost = Local_SolnBlk[Level_Coarse][b].Nghost;

      /* Enforce Neumann boundary conditions as appropriate */
      if (!Input_Parameters->Multigrid_IP.PROLONGATE_USING_INJECTION) {
	for (int j_coarse = Local_SolnBlk[Level_Coarse][b].JCl; j_coarse <= Local_SolnBlk[Level_Coarse][b].JCu; j_coarse++) {
	  
	  // West Boundary
	  if (Local_SolnBlk[Level_Coarse][b].Grid.BCtypeW[j_coarse] == BC_CONSTANT_EXTRAPOLATION) {
	    
	    Local_SolnBlk[Level_Coarse][b].U[Local_SolnBlk[Level_Coarse][b].ICl-1][j_coarse] =
	      Local_SolnBlk[Level_Coarse][b].U[Local_SolnBlk[Level_Coarse][b].ICl][j_coarse];
	  } /* end if */
	  // East Boundary
	  if (Local_SolnBlk[Level_Coarse][b].Grid.BCtypeE[j_coarse] == BC_CONSTANT_EXTRAPOLATION) {
	    Local_SolnBlk[Level_Coarse][b].U[Local_SolnBlk[Level_Coarse][b].ICu+1][j_coarse] =
	      Local_SolnBlk[Level_Coarse][b].U[Local_SolnBlk[Level_Coarse][b].ICu][j_coarse];
	  } /* end if */
	} /* end for */
	for (int i_coarse = Local_SolnBlk[Level_Coarse][b].ICl; i_coarse <= Local_SolnBlk[Level_Coarse][b].ICu; i_coarse++) {
	  // South Boundary
	  if (Local_SolnBlk[Level_Coarse][b].Grid.BCtypeS[i_coarse] == BC_CONSTANT_EXTRAPOLATION) {
	    
	    Local_SolnBlk[Level_Coarse][b].U[i_coarse][Local_SolnBlk[Level_Coarse][b].JCl-1] =
	      Local_SolnBlk[Level_Coarse][b].U[i_coarse][Local_SolnBlk[Level_Coarse][b].JCl];
	  } /* end if */
	  // North Boundary
	  if (Local_SolnBlk[Level_Coarse][b].Grid.BCtypeN[i_coarse] == BC_CONSTANT_EXTRAPOLATION) {
	    Local_SolnBlk[Level_Coarse][b].U[i_coarse][Local_SolnBlk[Level_Coarse][b].JCu+1] =
	      Local_SolnBlk[Level_Coarse][b].U[i_coarse][Local_SolnBlk[Level_Coarse][b].JCu];
	  } /* end if */
	} /* end for */
      } /* end if */

      /* Loop through coarse cells */
      for (int i_coarse = Local_SolnBlk[Level_Coarse][b].ICl; i_coarse <= Local_SolnBlk[Level_Coarse][b].ICu; i_coarse++) {
	for (int j_coarse = Local_SolnBlk[Level_Coarse][b].JCl; j_coarse <= Local_SolnBlk[Level_Coarse][b].JCu; j_coarse++) {

	  /* Check if coarse grid cell is on a Dirichlet type Boundary.  If so, use
	     simple injection instead of bilinear interpolation for all
	     4 related fine cells. */
	  if ((i_coarse == Local_SolnBlk[Level_Coarse][b].ICl &&
	       Local_SolnBlk[Level_Coarse][b].Grid.BCtypeW[j_coarse] == BC_FIXED) ||
	      (i_coarse == Local_SolnBlk[Level_Coarse][b].ICu &&
	       Local_SolnBlk[Level_Coarse][b].Grid.BCtypeE[j_coarse] == BC_FIXED) ||
	      (j_coarse == Local_SolnBlk[Level_Coarse][b].JCl &&
	       Local_SolnBlk[Level_Coarse][b].Grid.BCtypeS[i_coarse] == BC_FIXED) ||
	      (j_coarse == Local_SolnBlk[Level_Coarse][b].JCu &&
	       Local_SolnBlk[Level_Coarse][b].Grid.BCtypeN[i_coarse] == BC_FIXED) ||
	      (i_coarse == Local_SolnBlk[Level_Coarse][b].ICl &&
	       Local_SolnBlk[Level_Coarse][b].Grid.BCtypeW[j_coarse] == BC_FIXED_TEMP_WALL) ||
	      (i_coarse == Local_SolnBlk[Level_Coarse][b].ICu &&
	       Local_SolnBlk[Level_Coarse][b].Grid.BCtypeE[j_coarse] == BC_FIXED_TEMP_WALL) ||
	      (j_coarse == Local_SolnBlk[Level_Coarse][b].JCl &&
	       Local_SolnBlk[Level_Coarse][b].Grid.BCtypeS[i_coarse] == BC_FIXED_TEMP_WALL) ||
	      (j_coarse == Local_SolnBlk[Level_Coarse][b].JCu &&
	       Local_SolnBlk[Level_Coarse][b].Grid.BCtypeN[i_coarse] == BC_FIXED_TEMP_WALL) ||
	      (i_coarse == Local_SolnBlk[Level_Coarse][b].ICl &&
	       Local_SolnBlk[Level_Coarse][b].Grid.BCtypeW[j_coarse] == BC_ADIABATIC_WALL) ||
	      (i_coarse == Local_SolnBlk[Level_Coarse][b].ICu &&
	       Local_SolnBlk[Level_Coarse][b].Grid.BCtypeE[j_coarse] == BC_ADIABATIC_WALL) ||
	      (j_coarse == Local_SolnBlk[Level_Coarse][b].JCl &&
	       Local_SolnBlk[Level_Coarse][b].Grid.BCtypeS[i_coarse] == BC_ADIABATIC_WALL) ||
	      (j_coarse == Local_SolnBlk[Level_Coarse][b].JCu &&
	       Local_SolnBlk[Level_Coarse][b].Grid.BCtypeN[i_coarse] == BC_ADIABATIC_WALL)) {
	    /* Find i,j indices of the SW Corner of coarse cell */
	    i_fine = 2*(i_coarse-Nghost)+Nghost; // SW Corner i index (fine)
	    j_fine = 2*(j_coarse-Nghost)+Nghost; // SW Corner j index (fine)
	    
	    sol_change = Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse];
	    Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine] += sol_change;
	    Local_SolnBlk[Level_Fine][b].U[i_fine+1][j_fine] += sol_change;
	    Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine+1] += sol_change;
	    Local_SolnBlk[Level_Fine][b].U[i_fine+1][j_fine+1] += sol_change;
	    error_flag = 0;	    	    

	  } else {

	    if (Input_Parameters->Multigrid_IP.PROLONGATE_USING_INJECTION) injection_for_this_coarse_cell = true;
	    
	    restart: ;

	    /*******************************************************************/
	    /*** Calculate cell-centred solution at SW Corner of coarse cell ***/
	    /*******************************************************************/
	    i_fine = 2*(i_coarse-Nghost)+Nghost; // SW Corner i index (fine)
	    j_fine = 2*(j_coarse-Nghost)+Nghost; // SW Corner j index (fine)
	    
	    if (!injection_for_this_coarse_cell) {
	      error_flag = Bilinear_Interpolation(
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc,
		 sol_change);
// 	      if (error_flag == 0) {
// 		error_flag = Prolong_NonSol(
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc,
// 		    sol_change);
// 	      }
	      
	      /* If fine grid cell-centre falls outside of quad in SW dir
		 formed by coarse grid cell-centres, try NW one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
// 		if (error_flag == 0) {
// 		  error_flag = Prolong_NonSol(
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		      sol_change);
// 		}
	      }
	      /* If fine grid cell-centre falls outside of quad in NW dir
		 formed by coarse grid cell-centres, try SE one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
	           Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
// 		if (error_flag == 0) {
// 		  error_flag = Prolong_NonSol(
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
// 		      Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		      sol_change);
// 		}
	      }
	      /* If fine grid cell-centre falls outside of quad in SE dir
		 formed by coarse grid cell-centres, try NE one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
// 		if (error_flag == 0) {
// 		  error_flag = Prolong_NonSol(
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		      sol_change);
// 		}
	      }

	      if (error_flag != 0) {
		injection_for_this_coarse_cell = true;
		goto restart;
	      }
	    } else {
	      // use simple injection instead of bilinear interpolation
	      sol_change = Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse];
	      error_flag = 0;
	    }

	    //------------------------ Rte2D Specific -------------------------//

	    // Check if update pushes a cell into negative density/pressure
	    // If so, keep on halfing update until after several times, set update
	    // to zero and continue.

            //sol_change = PROLONGATION_REDUCTION_FACTOR*sol_change;
	    if (Input_Parameters->Multigrid_IP.UPDATE_STABILITY_SWITCH) {
	      Rte2D_State new_sol;
	      new_sol = Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine] + sol_change;
	      if (new_sol.NegIntensity()) {
		int num_it = 1;
		while(num_it <= Input_Parameters->Multigrid_IP.NUMBER_OF_UPDATE_REDUCTIONS) {
		  sol_change = HALF*sol_change;
		  new_sol = Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine] + sol_change;
		  if (!new_sol.NegIntensity()) break;
		  num_it++;
		} /* end while */
		if (num_it > Input_Parameters->Multigrid_IP.NUMBER_OF_UPDATE_REDUCTIONS) sol_change = sol_change*ZERO;
		cout << num_it;
	      } /* endif */
	    } /* endif */
	    //---------------------- End Rte2D Specific -----------------------//

	    Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine] += sol_change;
	  
	    /**********************************************/
	    /*** Finished with SW fine-grid cell centre ***/
	    /**********************************************/
	    
	    /*******************************************************************/
	    /*** Calculate cell-centred solution at SE Corner of coarse cell ***/
	    /*******************************************************************/
	    
	    i_fine = 2*(i_coarse-Nghost)+Nghost+1; // SE Corner i index (fine)
	    j_fine = 2*(j_coarse-Nghost)+Nghost;   // SE Corner j index (fine)
	    
	    if (!injection_for_this_coarse_cell) {
	      error_flag = Bilinear_Interpolation(
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 sol_change);
// 	      if (error_flag == 0) {
// 		error_flag = Prolong_NonSol(
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		    sol_change);
// 	      }	     
 
	      /* If fine grid cell-centre falls outside of quad in SE dir
		 formed by coarse grid cell-centres, try NE one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
// 		if (error_flag == 0) {
// 		  error_flag = Prolong_NonSol(
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		      sol_change);
// 		}
	      }
	      /* If fine grid cell-centre falls outside of quad in NE dir
		 formed by coarse grid cell-centres, try SW one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
// 		if (error_flag == 0) {
// 		  error_flag = Prolong_NonSol(
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 		      Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		      sol_change);
// 		}
	      }
	      /* If fine grid cell-centre falls outside of quad in SW dir
		 formed by coarse grid cell-centres, try NW one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
// 		if (error_flag == 0) {
// 		  error_flag = Prolong_NonSol(
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		      sol_change);
// 		}
	      }

	      if (error_flag != 0) {
		injection_for_this_coarse_cell = true;
		goto restart;
	      }
	    } else {
	      // use simple injection instead of bilinear interpolation
	      sol_change = Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse];
	      error_flag = 0;
	    }

	    //------------------------ Rte2D Specific -------------------------//

	    // Check if update pushes a cell into negative density/pressure
	    // If so, keep on halfing update until after several times, set update
	    // to zero and continue.

            //sol_change = PROLONGATION_REDUCTION_FACTOR*sol_change;
	    if (Input_Parameters->Multigrid_IP.UPDATE_STABILITY_SWITCH) {
	      Rte2D_State new_sol;
	      new_sol = Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine] + sol_change;
	      if (new_sol.NegIntensity()) {
		int num_it = 1;
		while(num_it <= Input_Parameters->Multigrid_IP.NUMBER_OF_UPDATE_REDUCTIONS) {
		  sol_change = HALF*sol_change;
		  new_sol = Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine] + sol_change;
		  if (!new_sol.NegIntensity()) break;
		  num_it++;
		} /* end while */
		if (num_it > Input_Parameters->Multigrid_IP.NUMBER_OF_UPDATE_REDUCTIONS) sol_change = sol_change*ZERO;
		cout << num_it;
	      } /* endif */
	    } /* endif */
	    //---------------------- End Rte2D Specific -----------------------//

	    Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine] += sol_change;
	    
	    /**********************************************/
	    /*** Finished with SE fine-grid cell centre ***/
	    /**********************************************/
	    
	    /*******************************************************************/
	    /*** Calculate cell-centred solution at NE Corner of coarse cell ***/
	    /*******************************************************************/
	    
	    i_fine = 2*(i_coarse-Nghost)+Nghost+1; // NE Corner i index (fine)
	    j_fine = 2*(j_coarse-Nghost)+Nghost+1; // NE Corner j index (fine)
	    
	    if (!injection_for_this_coarse_cell) {
	      error_flag = Bilinear_Interpolation(
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 sol_change);
// 		if (error_flag == 0) {
// 		  error_flag = Prolong_NonSol(
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		      sol_change);
// 		}
	      
	      /* If fine grid cell-centre falls outside of quad in NE dir
		 formed by coarse grid cell-centres, try SE one */
	      if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 sol_change);
// 		if (error_flag == 0) {
// 		  error_flag = Prolong_NonSol(
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
// 		      Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		      sol_change);
// 		}
	      }
	      /* If fine grid cell-centre falls outside of quad in SE dir
		 formed by coarse grid cell-centres, try NW one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
// 		if (error_flag == 0) {
// 		  error_flag = Prolong_NonSol(
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		      sol_change);
// 		}
	      }
	    /* If fine grid cell-centre falls outside of quad in NW dir
	       formed by coarse grid cell-centres, try SW one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
// 		if (error_flag == 0) {
// 		  error_flag = Prolong_NonSol(
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 		      Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		      sol_change);
// 		}
	      }

	      if (error_flag != 0) {
		injection_for_this_coarse_cell = true;
		goto restart;
	      }
	    } else {
	      // use simple injection instead of bilinear interpolation
	      sol_change = Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse];
	      error_flag = 0;
	    }

	    //------------------------ Rte2D Specific -------------------------//

	    // Check if update pushes a cell into negative density/pressure
	    // If so, keep on halfing update until after several times, set update
	    // to zero and continue.

            //sol_change = PROLONGATION_REDUCTION_FACTOR*sol_change;
	    if (Input_Parameters->Multigrid_IP.UPDATE_STABILITY_SWITCH) {
	      Rte2D_State new_sol;
	      new_sol = Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine] + sol_change;
	      if (new_sol.NegIntensity()) {
		int num_it = 1;
		while(num_it <= Input_Parameters->Multigrid_IP.NUMBER_OF_UPDATE_REDUCTIONS) {
		  sol_change = HALF*sol_change;
		  new_sol = Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine] + sol_change;
		  if (!new_sol.NegIntensity()) break;
		  num_it++;
		} /* end while */
		if (num_it > Input_Parameters->Multigrid_IP.NUMBER_OF_UPDATE_REDUCTIONS) sol_change = sol_change*ZERO;
		cout << num_it;
	      } /* endif */
	    } /* endif */
	    //---------------------- End Rte2D Specific -----------------------//

	    Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine] += sol_change;
	    
	    /**********************************************/
	    /*** Finished with NE fine-grid cell centre ***/
	    /**********************************************/
	    
	    /*******************************************************************/
	    /*** Calculate cell-centred solution at NW Corner of coarse cell ***/
	    /*******************************************************************/
	    
	    i_fine = 2*(i_coarse-Nghost)+Nghost;   // NW Corner i index (fine)
	    j_fine = 2*(j_coarse-Nghost)+Nghost+1; // NW Corner j index (fine)
	    
	    if (!injection_for_this_coarse_cell) {
	      error_flag = Bilinear_Interpolation(
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 sol_change);
// 	      if (error_flag == 0) {
// 		error_flag = Prolong_NonSol(
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 		    Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		    Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		    sol_change);
// 	      }
	      
	      /* If fine grid cell-centre falls outside of quad in NW dir
		 formed by coarse grid cell-centres, try NE one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
// 		if (error_flag == 0) {
// 		  error_flag = Prolong_NonSol(
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		      sol_change);
// 		}
	      }
	      /* If fine grid cell-centre falls outside of quad in NE dir
		 formed by coarse grid cell-centres, try SW one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
// 		if (error_flag == 0) {
// 		  error_flag = Prolong_NonSol(
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 		      Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		      sol_change);
// 		}
	      }
	      /* If fine grid cell-centre falls outside of quad in SW dir
		 formed by coarse grid cell-centres, try SE one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1],
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
// 		if (error_flag == 0) {
// 		  error_flag = Prolong_NonSol(
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 		      Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1],
// 		      Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
// 		      Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
// 		      sol_change);
// 		}
	      }

	      if (error_flag != 0) {
		injection_for_this_coarse_cell = true;
		goto restart;
	      }
	    } else {

	      // use simple injection instead of bilinear interpolation
	      sol_change = Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse];
	      error_flag = 0;
	    }

	    //------------------------ Rte2D Specific -------------------------//

	    // Check if update pushes a cell into negative intensity
	    // If so, keep on halfing update until after several times, set update
	    // to zero and continue.

            //sol_change = PROLONGATION_REDUCTION_FACTOR*sol_change;
	    if (Input_Parameters->Multigrid_IP.UPDATE_STABILITY_SWITCH) {
	      Rte2D_State new_sol;
	      new_sol = Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine] + sol_change;

	      // if there is a negative value
	      if (new_sol.NegIntensity()) {
		int num_it = 1;
		while(num_it <= Input_Parameters->Multigrid_IP.NUMBER_OF_UPDATE_REDUCTIONS) {
		  sol_change = HALF*sol_change;
		  new_sol = Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine] + sol_change;
		  if (!new_sol.NegIntensity()) break;
		  num_it++;
		} /* end while */
		if (num_it > Input_Parameters->Multigrid_IP.NUMBER_OF_UPDATE_REDUCTIONS) sol_change = sol_change*ZERO;
		cout << num_it;
	      } /* end if */
	    

	    } /* end if */
	    //---------------------- End Rte2D Specific -----------------------//


	    Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine] += sol_change;
	  
	    /**********************************************/
	    /*** Finished with NW fine-grid cell centre ***/
	    /**********************************************/
	    
	    // reset injection flag
	    injection_for_this_coarse_cell = false;
	    
	  } /* end if */
	} /* end for */
      } /* endfor */
    } /* endif */
  } /* endfor */

  /* No Errors... prolongation completed */
  return(error_flag);
}




/**********************************************************
 * Routine: Update_Primitive_Variables                    *
 *                                                        *
 * This routine updates all primitive variables W using   *
 * the conserved variables U, for all solution blocks on  *
 * a level.                                               *
 **********************************************************/
template <> inline void FAS_Multigrid2D_Solver<Rte2D_State,
			 		       Rte2D_Quad_Block,
					       Rte2D_Input_Parameters>::
Update_Primitive_Variables(const int Level) {  /* DO NOTHING */ }


/********************************************************
 * Routine: dUdt_Multistage_Explicit_for_Multigrid      *
 *                                                      *
 * This routine updates the solution for a 1D array of  *
 * 2D quadrilateral multi-block solution blocks.  A     *
 * variety of multistage explicit time integration      *
 * and a 2nd-ororder limited upwind finite-volume       *
 * spatial discretization scheme for the convective     *
 * flux coupled with a centrally-weighted finite-volume *
 * discretization for the diffused flux can be used     *
 * depending on the specified input values.             *
 *                                                      *
 ********************************************************/
template <> inline int FAS_Multigrid2D_Solver<Rte2D_State,
			                      Rte2D_Quad_Block,
			                      Rte2D_Input_Parameters>::
dUdt_Multistage_Explicit_for_Multigrid(const int Level,
				       const int i_stage) {

  int flux_function_type, limiter_type, k_residual;

  /* Force use of first-order spatial reconstruction and HLLE
     flux function on coarse levels. */

//   flux_function_type = Input_Parameters->i_Flux_Function;
  limiter_type = Input_Parameters->i_Limiter;
  if (Level > Top_Level &&
      Input_Parameters->Multigrid_IP.FORCE_FIRST_ORDER_SPATIAL_RECONSTRUCTION == true) {
//     if (Input_Parameters->Local_Time_Stepping != LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER ||
//         ((Input_Parameters->Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) &&
//          (Input_Parameters->i_Flux_Function != FLUX_FUNCTION_ROE_PRECON_WS)) &&
//         ((Input_Parameters->Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) &&
//          (Input_Parameters->i_Flux_Function != FLUX_FUNCTION_HLLE_PRECON_WS)))
//       Input_Parameters->i_Flux_Function = FLUX_FUNCTION_HLLE;
    Input_Parameters->i_Limiter = LIMITER_ZERO;
  }

  /* Temporarily overwrite Input_Parameters->i_Time_Integration. */
  Input_Parameters->i_Time_Integration = Input_Parameters->Multigrid_IP.i_Smoothing;

  /* Evaluate the residual for each solution block. */

  error_flag = dUdt_Multistage_Explicit(Local_SolnBlk[Level],
					List_of_Local_Solution_Blocks[Level],
					*Input_Parameters,
					i_stage);
  if (error_flag) return (error_flag);

  /* Add the defect correction forcing term as necessary. */

  if (Level > Top_Level) {
     switch(Input_Parameters->i_Time_Integration) {
     case TIME_STEPPING_EXPLICIT_EULER :
       k_residual = 0;
       break;
     case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
       k_residual = 0;
       break;
     case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
       k_residual = 0;
       if (Input_Parameters->N_Stage == 4) {
         if (i_stage == 4) {
	   k_residual = 0;
         } else {
	   k_residual = i_stage - 1;
         } /* endif */
       } /* endif */
       break;
     case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
       k_residual = 0;
       break;
     default:
       k_residual = 0;
       break;
     } /* endswitch */

     for (int b = 0 ; b <= List_of_Local_Solution_Blocks[Level].Nblk-1 ; ++b ) {
       if (List_of_Local_Solution_Blocks[Level].Block[b].used == ADAPTIVEBLOCK2D_USED) {      
         /* Add Multigrid Forcing Term */
         for (int j  = Local_SolnBlk[Level][b].JCl ; j <= Local_SolnBlk[Level][b].JCu ; ++j ) {
	   for (int i = Local_SolnBlk[Level][b].ICl ; i <= Local_SolnBlk[Level][b].ICu ; ++i ) {
	     Local_SolnBlk[Level][b].dUdt[i][j][k_residual] +=
	       (Input_Parameters->CFL_Number*Local_SolnBlk[Level][b].dt[i][j])*
	       MG_Blk[Level][b].P[i][j];
   	   } /* endfor */
         } /* endfor */
       } /* endif */
     }  /* endfor */
  } /* endif */

  /* Reset flux function, limiter, and time integration parameters. */

  if (Level > Top_Level &&
      Input_Parameters->Multigrid_IP.FORCE_FIRST_ORDER_SPATIAL_RECONSTRUCTION == true) {
//     Input_Parameters->i_Flux_Function = flux_function_type;
    Input_Parameters->i_Limiter = limiter_type;
  }
  Input_Parameters->i_Time_Integration = TIME_STEPPING_MULTIGRID;

  /* Quadrilateral multi-block solution residuals
     successfully calculated.  Return. */
  return (error_flag = 0);
}



/******************************************************************************
 * Routine: Apply_Boundary_Flux_Corrections_Multistage_Explicit_for_Multigrid *
 *                                                                            *
 * Apply flux corrections at boundaries of a 1D array                         *
 * of 2D quadrilateral multi-block solution blocks to                         *
 * ensure that the scheme is conservative at boundaries                       *
 * with resolution mesh changes.                                              *
 *                                                                            *
 ******************************************************************************/
template <> inline void FAS_Multigrid2D_Solver<Rte2D_State,
					       Rte2D_Quad_Block,
					       Rte2D_Input_Parameters>::
Apply_Boundary_Flux_Corrections_Multistage_Explicit_for_Multigrid(const int Level,
								  const int i_stage) {

  int flux_function_type, limiter_type;
  
  /* Force use of first-order spatial reconstruction and HLLE
     flux function on coarse levels. */

//   flux_function_type = Input_Parameters->i_Flux_Function;
  limiter_type = Input_Parameters->i_Limiter;
  if (Level > Top_Level &&
      Input_Parameters->Multigrid_IP.FORCE_FIRST_ORDER_SPATIAL_RECONSTRUCTION == true) {
//     if (Input_Parameters->Local_Time_Stepping != LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER ||
//         ((Input_Parameters->Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) &&
//          (Input_Parameters->i_Flux_Function != FLUX_FUNCTION_ROE_PRECON_WS)) &&
//         ((Input_Parameters->Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) &&
//          (Input_Parameters->i_Flux_Function != FLUX_FUNCTION_HLLE_PRECON_WS)))
//        Input_Parameters->i_Flux_Function = FLUX_FUNCTION_HLLE;
    Input_Parameters->i_Limiter = LIMITER_ZERO;
  }

  /* Temporarily overwrite Input_Parameters->i_Time_Integration. */
  Input_Parameters->i_Time_Integration = Input_Parameters->Multigrid_IP.i_Smoothing;
        
  /* Prescribe boundary data for each solution block. */

  Apply_Boundary_Flux_Corrections_Multistage_Explicit(Local_SolnBlk[Level],
						      List_of_Local_Solution_Blocks[Level],
						      *Input_Parameters,
						      i_stage);

  /* Reset flux function, limiter, and time integration parameters. */

  if (Level > Top_Level &&
      Input_Parameters->Multigrid_IP.FORCE_FIRST_ORDER_SPATIAL_RECONSTRUCTION == true) {
//     Input_Parameters->i_Flux_Function = flux_function_type;
    Input_Parameters->i_Limiter = limiter_type;
  }
  Input_Parameters->i_Time_Integration = TIME_STEPPING_MULTIGRID;

}



/**************************************************************
 * Routine: Update_Solution_Multistage_Explicit_for Multigrid *
 *                                                            *
 * This routine updates the solution for a 1D array of        *
 * 2D quadrilateral multi-block solution blocks.  A           *
 * variety of multistage explicit time integration            *
 * and upwind finite-volume spatial discretization            *
 * procedures can be used depending on the specified          *
 * input values.                                              *
 *                                                            *
 **************************************************************/
template <> inline int FAS_Multigrid2D_Solver<Rte2D_State,
				              Rte2D_Quad_Block,
				              Rte2D_Input_Parameters>::
Update_Solution_Multistage_Explicit_for_Multigrid(const int Level,
						  const int i_stage) {
  
  int b, error_flag, flux_function_type, limiter_type;
  
  error_flag = 0;
  
  /* Force use of first-order spatial reconstruction and HLLE
     flux function on coarse levels. */

//   flux_function_type = Input_Parameters->i_Flux_Function;
  limiter_type = Input_Parameters->i_Limiter;
  if (Level > Top_Level &&
      Input_Parameters->Multigrid_IP.FORCE_FIRST_ORDER_SPATIAL_RECONSTRUCTION == true) {
//     if (Input_Parameters->Local_Time_Stepping != LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER ||
//         ((Input_Parameters->Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) &&
//          (Input_Parameters->i_Flux_Function != FLUX_FUNCTION_ROE_PRECON_WS)) &&
//         ((Input_Parameters->Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) &&
//          (Input_Parameters->i_Flux_Function != FLUX_FUNCTION_HLLE_PRECON_WS)))
//        Input_Parameters->i_Flux_Function = FLUX_FUNCTION_HLLE;
    Input_Parameters->i_Limiter = LIMITER_ZERO;
  }

  /* Temporarily overwrite Input_Parameters->i_Time_Integration. */
  Input_Parameters->i_Time_Integration = Input_Parameters->Multigrid_IP.i_Smoothing;
  
  /* Update the solution for each solution block. */

  error_flag = Update_Solution_Multistage_Explicit(Local_SolnBlk[Level],
						   List_of_Local_Solution_Blocks[Level],
						   *Input_Parameters,
						   i_stage);
  if (error_flag) return (error_flag);

  /* Reset flux function, limiter, and time integration parameters. */

  if (Level > Top_Level &&
      Input_Parameters->Multigrid_IP.FORCE_FIRST_ORDER_SPATIAL_RECONSTRUCTION == true) {
//     Input_Parameters->i_Flux_Function = flux_function_type;
    Input_Parameters->i_Limiter = limiter_type;
  }
  Input_Parameters->i_Time_Integration = TIME_STEPPING_MULTIGRID;

  /* Quadrilateral multi-block solution blocks
     successfully updated.  Return. */
  
  return(error_flag);

}


/********************************************************
 * Routine: Residual_Evaluation_for_Multigrid           *
 *                                                      *
 * This routine evaluates the residual for a 1D array   *
 * of solution blocks given the                         *
 * solution U using a higher-order limited upwind       *
 * finite-volume spatial discretization scheme for the  *
 * convective flux coupled with a centrally-weighted    *
 * finite-volume discretization for the diffused flux.  *
 *                                                      *
 ********************************************************/
template <> inline int FAS_Multigrid2D_Solver<Rte2D_State,
				              Rte2D_Quad_Block,
				              Rte2D_Input_Parameters>::
Residual_Evaluation_for_Multigrid(const int Level,
				  const int flag_for_P_term) {

  int flux_function_type, limiter_type;

  /* Force use of first-order spatial reconstruction and HLLE
     flux function on coarse levels. */

//   flux_function_type = Input_Parameters->i_Flux_Function;
  limiter_type = Input_Parameters->i_Limiter;
  if (Level > Top_Level &&
      Input_Parameters->Multigrid_IP.FORCE_FIRST_ORDER_SPATIAL_RECONSTRUCTION == true) {
//     if (Input_Parameters->Local_Time_Stepping != LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER ||
//         ((Input_Parameters->Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) &&
//          (Input_Parameters->i_Flux_Function != FLUX_FUNCTION_ROE_PRECON_WS)) &&
//         ((Input_Parameters->Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) &&
//          (Input_Parameters->i_Flux_Function != FLUX_FUNCTION_HLLE_PRECON_WS)))
//        Input_Parameters->i_Flux_Function = FLUX_FUNCTION_HLLE;
    Input_Parameters->i_Limiter = LIMITER_ZERO;
  }

  /* Evaluate the residual for each solution block. */

  for (int b = 0 ; b <= List_of_Local_Solution_Blocks[Level].Nblk-1 ; ++b ) {
    if (List_of_Local_Solution_Blocks[Level].Block[b].used == ADAPTIVEBLOCK2D_USED) {
      // Evaluate Residual
      dUdt_Residual_Evaluation(Local_SolnBlk[Level][b],
			       *(Input_Parameters));

      // Add in forcing term if needed
      if (flag_for_P_term > 0) {
	for (int j  = Local_SolnBlk[Level][b].JCl ; j <= Local_SolnBlk[Level][b].JCu ; ++j ) {
	  for (int i = Local_SolnBlk[Level][b].ICl ; i <= Local_SolnBlk[Level][b].ICu ; ++i ) {
	    Local_SolnBlk[Level][b].dUdt[i][j][0] += MG_Blk[Level][b].P[i][j];
	  } /* endfor */
	} /* endfor */
      } else {
      } /* endif */

    } /* endif */
  }  /* endfor */


  /* Reset flux function and limiter. */

  if (Level > Top_Level &&
      Input_Parameters->Multigrid_IP.FORCE_FIRST_ORDER_SPATIAL_RECONSTRUCTION == true) {
//     Input_Parameters->i_Flux_Function = flux_function_type;
    Input_Parameters->i_Limiter = limiter_type;
  }

  // Send boundary flux corrections at block interfaces with resolution changes.
  error_flag = Send_Conservative_Flux_Corrections(Local_SolnBlk[Level],
						  List_of_Local_Solution_Blocks[Level],
						  Local_SolnBlk[Level][0].NumVar());
  if (error_flag) {
    cout << "\n FASMultigrid2D ERROR: flux correction message passing error on processor "
	 << List_of_Local_Solution_Blocks[Level].ThisCPU
	 << ".\n";
    cout.flush();
  } /* endif */
  error_flag = CFDkit_OR_MPI(error_flag);
  if (error_flag) return (error_flag);
  
  // Apply boundary flux corrections to ensure that method is conservative.
  Apply_Boundary_Flux_Corrections(Local_SolnBlk[Level],
				  List_of_Local_Solution_Blocks[Level]);
  
  /* Residuals successfully calculated for quadrilateral multi-block solution blocks
     Return. */

  return (error_flag);

}


/********************************************************
 * Routine: Coarse_Grid_Correction                      *
 *                                                      *
 * Implements the Full Approximation Storage Multigrid  *
 * Algorithm                                            *
 *                                                      *
 ********************************************************/
template <> int FAS_Multigrid2D_Solver<Rte2D_State,
				       Rte2D_Quad_Block,
				       Rte2D_Input_Parameters>::
Coarse_Grid_Correction(const int top_Level,
		       const int Current_Level) {

  int N_Smooths = 0;
  int N_MCycle = 0;
  double dTime;

  /* Store current top level */
  if (Current_Level == top_Level) {
    Top_Level = top_Level;
  }

#ifdef _MG_COUT_CYCLE_STAGE_
  for (int l = 0; l<Current_Level; l++) cout << "   ";
  cout << "Entering Level " << Current_Level << endl;  
#endif

  /* Check if on the coarsest level */
  if (Current_Level == Input_Parameters->Multigrid_IP.Levels-1) {

    /* Set flag_for_first_time_restriction to false, since we've reached
       the bottom at least once */
    flag_for_First_Time_Restriction = false;

    /* Store the unsmoothed solution for later use. */

    Store_Current_Solution_in_uo_MG(Current_Level);

#ifdef _MG_COUT_CYCLE_STAGE_
    for (int l = 0; l<Current_Level; l++) cout << "   ";
    cout << "Coarse Smoothing" << endl;
#endif
    
    /* Smooth on the coarsest grid */
    for (int nu = 1; nu <= Input_Parameters->Multigrid_IP.Number_of_Smooths_on_Coarsest_Level; nu++) {
      error_flag = Smooth(Current_Level);
      
      if (error_flag) return (error_flag);
    } /* endfor */
    
  } else { /* Not on the coarsest level */
    
    /* Check to see which type of cycle we're performing */
    if (Input_Parameters->Multigrid_IP.i_Cycle == MULTIGRID_V_CYCLE || Current_Level == 0) {
      N_MCycle = 1;
    } else if (Input_Parameters->Multigrid_IP.i_Cycle == MULTIGRID_W_CYCLE) {
      N_MCycle = 2;
    } else {
      return error_flag == 1;
    } /* end if */

    /* Implement V-type (N_MCycle = 1) or W-type (N_MCycle = 2) Cycle */
    for (int m = 1; m <= N_MCycle; m++) {

#ifdef _MG_COUT_CYCLE_STAGE_
      for (int l = 0; l<Current_Level; l++) cout << "   ";
      cout << "m = " << m << endl;
#endif

      /* Check to see if we're on the top level, if so, reset 
	 flag_for_First_Time_Restriction to true */
      if (Current_Level == Top_Level) {
	flag_for_First_Time_Restriction = true;
      }

      /* If this is the first time restrict from this level, store
	 the unsmoothed solution for later use. */
      /* and calculate time step sizes for next coarser level */
      if (flag_for_First_Time_Restriction == true) {
#ifdef _MG_COUT_CYCLE_STAGE_
	for (int l = 0; l<Current_Level; l++) cout << "   ";
	cout << "Storing uo_MG" << endl;
#endif
	Store_Current_Solution_in_uo_MG(Current_Level);	
      }

      /*****************************
       * PreSmooth N_Smooths times *
       *****************************/

      /* Check if on the top level and set the number of iterations */
      N_Smooths = Current_Level == Top_Level ?
	Input_Parameters->Multigrid_IP.Number_of_Smooths_on_Finest_Level :
	Input_Parameters->Multigrid_IP.Number_of_Pre_Smooths;

#ifdef _MG_COUT_CYCLE_STAGE_
      for (int l = 0; l<Current_Level; l++) cout << "   ";
      cout << "Presmoothing " << N_Smooths << " times" << endl;
#endif

      for (int nu1 = 1; nu1 <= N_Smooths; nu1++) {
	//	cout << nu1 << endl;
	error_flag = Smooth(Current_Level);

	if (error_flag) return (error_flag);
      } /* endfor */

      /************************************
       * Restrict solution to coarse grid *
       ************************************/
#ifdef _MG_COUT_CYCLE_STAGE_
      for (int l = 0; l<Current_Level; l++) cout << "   ";
      cout << "Restrict Solution " << Current_Level << endl;
#endif
      Restrict_Solution_Blocks(Current_Level);
      Update_Primitive_Variables(Current_Level+1);

      // Exchange solution information between neighbouring blocks.
#ifdef _MG_COUT_CYCLE_STAGE_
      for (int l = 0; l<Current_Level; l++) cout << "   ";
      cout << "Message-Passing, post Restriction, on level " << Current_Level+1 << endl;
#endif
      error_flag = Exchange_solution_information(Current_Level+1,
						 OFF);

      if (error_flag) return (error_flag);

      // Apply boundary conditions for stage.
#ifdef _MG_COUT_CYCLE_STAGE_
      for (int l = 0; l<Current_Level; l++) cout << "   ";
      cout << "Apply BCs on level " << Current_Level+1 << endl;
#endif
      BCs(Local_SolnBlk[Current_Level+1], 
	  List_of_Local_Solution_Blocks[Current_Level+1]);
      Prescribe_NonSol(Local_SolnBlk[Current_Level+1], 
		       List_of_Local_Solution_Blocks[Current_Level+1],
		       *Input_Parameters);

      /************************************************************************
       * Evaluate Forcing Term P by first restricting the fine grid residuals *
       * from dUdt[][][0] on Level_Fine to P on Level_Coarse.                 *
       * P=Restrict(Fine Grid Residual)-Residual(Restrict(Fine Grid Solution))*
       ************************************************************************/
       if (Input_Parameters->Multigrid_IP.INCLUDE_FINE_TO_COARSE_DEFECT_CORRECTIONS) {
#ifdef _MG_COUT_CYCLE_STAGE_
 	  for (int l = 0; l<Current_Level; l++) cout << "   ";
	  cout << "Form defect corrections, P " << endl;
#endif
          // Restrict residual from fine to coarse grid.
	  Restrict_Residuals(Current_Level);

	  /*******************************************************
	   * Evaluate residuals of restricted fine grid solution *
	   *******************************************************/
	  error_flag = Residual_Evaluation_for_Multigrid(Current_Level+1,
		  				         0);	
	  if (error_flag) return (error_flag);

  	  /*************************************************************************
	   * Finish Creating Forcing Term P by subtracting the residual evaluated  *
	   * from the restricted fine grid solution with P, which already contains *
	   * the restricted fine grid residuals.                                   *
	   * P=Restrict(Fine Grid Residual)-Residual(Restrict(Fine Grid Solution)) *
	   *************************************************************************/
	  Subtract_dUdt_from_P(Current_Level+1);
       } /* endif */

      /*************************************************************************
       * Call Solver Routine recursively to solve the coarse grid level problem*
       *************************************************************************/
#ifdef _MG_COUT_CYCLE_STAGE_
      for (int l = 0; l<Current_Level; l++) cout << "   ";
      cout << "Heading Down from Level " << Current_Level << endl;
#endif      

      error_flag = Coarse_Grid_Correction(Top_Level,
					  Current_Level+1);

      if (error_flag) return (error_flag);

#ifdef _MG_COUT_CYCLE_STAGE_
       for (int l = 0; l<Current_Level; l++) cout << "   ";
       cout << "Coming back up to Level " << Current_Level << endl;
#endif
      
      /************************************************************************
       * Prolong Solution changes, stored in U[][].u from Current_Level+1     *
       * back to Current_Level, and update solution U[][].u on this level     *
       ************************************************************************/
#ifdef _MG_COUT_CYCLE_STAGE_
      for (int l = 0; l<Current_Level; l++) cout << "   ";
      cout << "Prolong & Update " << endl;
#endif

      error_flag = Prolong_and_Update_Solution_Blocks(Current_Level+1);

      if (error_flag) return (error_flag);

      /* Update Primitive variables for current level */

      Update_Primitive_Variables(Current_Level);

      /* Exchange messages and apply BCs */

      error_flag = Exchange_solution_information(Current_Level,
						 OFF);

      if (error_flag) return (error_flag);
      
      // Apply boundary conditions for level.
      BCs(Local_SolnBlk[Current_Level], 
	  List_of_Local_Solution_Blocks[Current_Level]);
      Prescribe_NonSol(Local_SolnBlk[Current_Level], 
	  List_of_Local_Solution_Blocks[Current_Level],
	  *Input_Parameters);

      /*****************************
       * PostSmooth N_Smooths times *
       *****************************/

      /* Check if on the top level and set the number of iterations */
      N_Smooths = Current_Level == Top_Level ? 0 :
	Input_Parameters->Multigrid_IP.Number_of_Post_Smooths;

#ifdef _MG_COUT_CYCLE_STAGE_
      for (int l = 0; l<Current_Level; l++) cout << "   ";
      cout << "Postsmooth " << N_Smooths << " times" << endl;
#endif

      for (int nu2 = 1; nu2 <= N_Smooths; nu2++) {
	error_flag = Smooth(Current_Level);

	if (error_flag) return (error_flag);
      } /* endfor */

    } /* endfor */
  } /* end if */

  /* if not yet back on the top level, and on the way up, calculate 
     solution changes */
  if (Current_Level != Top_Level) {
    /*************************************************************
     *  Calculate Solution changes = current U[][] - uo_MG[][] *
     *************************************************************/

#ifdef _MG_COUT_CYCLE_STAGE_
    for (int l = 0; l<Current_Level; l++) cout << "   ";
    cout << "Evaluate Solution changes " << endl;
#endif

    Evaluate_Solution_Changes(Current_Level);
    
    /* Solution changes ready to be prolonged to the finer level... */

  } /* end if */
  
#ifdef _MG_COUT_CYCLE_STAGE_
    for (int l = 0; l<Current_Level; l++) cout << "   ";
    cout << "Leaving Level " << Current_Level << endl;
#endif

  return (error_flag);
}


/********************************************************
 * Routine: Execute                                     *
 *                                                      *
 * Commence computation using FAS multigrid method      *
 *                                                      *
 ********************************************************/
template <> int FAS_Multigrid2D_Solver<Rte2D_State,
				       Rte2D_Quad_Block,
				       Rte2D_Input_Parameters>::
Execute(int& number_of_time_steps,
	double& Time,
	CPUTime& processor_cpu_time,
	CPUTime& total_cpu_time,
	ofstream& residual_file) {

  /* Other local solution variables. */

  int first_step, line_number, command_flag, level, limiter_freezing;

  /********************************************************  
   * Initialize solution variables.                       *
   ********************************************************/
      
  /* Initialize all coarser levels of solution blocks with ICs */
  
  for (level = 1; level < Input_Parameters->Multigrid_IP.Levels; level++) {
    ICs(Local_SolnBlk[level], 
	List_of_Local_Solution_Blocks[level], 
	*Input_Parameters);
    Restrict_Boundary_Ref_States(level-1);
  }
  
  /* Send solution information between neighbouring blocks to complete
     prescription of initial data. */
  
  for (level = 1; level < Input_Parameters->Multigrid_IP.Levels; level++) {
    
    error_flag = Exchange_solution_information(level,ON);
    
    if (error_flag) return (error_flag);
  }

  /********************************************************  
   * Solve IBVP or BVP on multi-block solution-adaptive   *
   * quadrilateral mesh using Multigrid for time-marching *
   ********************************************************/
  
  CFDkit_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

  /* Open residual file and reset the CPU time. */
  
  first_step = 1;  
  limiter_freezing = OFF;
  
  if (CFDkit_Primary_MPI_Processor()) {
    error_flag = Open_Progress_File(residual_file,
				    Input_Parameters->Output_File_Name,
				    number_of_time_steps);
    if (error_flag) {
      cout << "\n FASMultigrid ERROR: Unable to open residual file for calculation.\n";
      cout.flush();
    } /* endif */
  } /* endif */
  
  CFDkit_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
  CFDkit_Broadcast_MPI(&error_flag, 1);
  if (error_flag) return (error_flag);
  
  processor_cpu_time.reset();
  
  /* Perform required number of iterations (time steps). */
  
  if (!Input_Parameters->Time_Accurate &&
      Input_Parameters->Maximum_Number_of_Time_Steps > 0) {

    if (!batch_flag) cout << "\n Beginning FAS Multigrid computations on "
			  << Date_And_Time() << ".";
    
    /* Perform the required number of Full & Regular Multigrid cycles */
    
    /* Determine if Full Multigrid is required */
    int initial_top_level = FINEST_LEVEL;
    if (Input_Parameters->Multigrid_IP.Number_of_Cycles_per_Stage_for_Full_Multigrid > 0) {
      if (!batch_flag) cout << "\n\n Perform Full multigrid cycles\n\n";
      initial_top_level = Input_Parameters->Multigrid_IP.Levels-2;                          // WHY ALWAYS -2 ???
    }
    
    /* Loop through each top_level (for FMG) */
    for (int top_level = initial_top_level; top_level >= FINEST_LEVEL; top_level--) {

      if (!batch_flag && top_level == FINEST_LEVEL) 
	cout << "\n\n Perform Regular multigrid cycles\n";
      
      /* Unfreeze limiters if frozen */
      if (!first_step &&
	  Input_Parameters->Freeze_Limiter &&
	  limiter_freezing == ON &&
	  Input_Parameters->Limiter_Type != LIMITER_ZERO) {
	for (level = top_level; level < Input_Parameters->Multigrid_IP.Levels; level++) {
	  Evaluate_Limiters(Local_SolnBlk[level], 
			    List_of_Local_Solution_Blocks[level]);
	} /* end for */
	limiter_freezing = OFF;
      }

      /* Perform the appropriate number of cycles for this level */

      /* Set the appropriate number of cycles for this level */
      /* Note:  cycles_for_this_level is set to an arbitrarily
	 large number if the current_top_level is the finest_level.
	 A separate exit criterion based on number_of_time_steps is
	 used to ensure finishing at the right time */

#ifdef _GNU_GCC_V3      
      unsigned long cycles_for_this_level =
         numeric_limits<unsigned long>::max();
#else
      unsigned long cycles_for_this_level = 10000000;
#endif
      if (Input_Parameters->Multigrid_IP.Number_of_Cycles_per_Stage_for_Full_Multigrid > 0) { // FMG on?
	if (top_level != FINEST_LEVEL) {  // FMG - Not on the finest level
	  cycles_for_this_level = Input_Parameters->Multigrid_IP.Number_of_Cycles_per_Stage_for_Full_Multigrid;
	} else { // FMG - On the finest level
#ifdef _GNU_GCC_V3
	  cycles_for_this_level = numeric_limits<unsigned long>::max();
#else
	  cycles_for_this_level = 10000000;
#endif 
	} /* end if */
      } /* end if */

      for (int cycles = 1; cycles <= cycles_for_this_level; cycles++) {
		  
	/* Determine the L1, L2, and max norms of the solution residual. */
	residual_l1_norm = L1_Norm_Residual(Local_SolnBlk[top_level], 
					    List_of_Local_Solution_Blocks[top_level]);
	residual_l1_norm = CFDkit_Summation_MPI(residual_l1_norm); // L1 norm for all processors.
	
	residual_l2_norm = L2_Norm_Residual(Local_SolnBlk[top_level], 
					    List_of_Local_Solution_Blocks[top_level]);
	residual_l2_norm = sqr(residual_l2_norm);
	residual_l2_norm = CFDkit_Summation_MPI(residual_l2_norm); // L2 norm for all processors.
	residual_l2_norm = sqrt(residual_l2_norm);
	
	residual_max_norm = Max_Norm_Residual(Local_SolnBlk[top_level], 
					      List_of_Local_Solution_Blocks[top_level]);
	residual_max_norm = CFDkit_Maximum_MPI(residual_max_norm); // Max norm for all processors.
	
	/* Update CPU time used for the calculation so far. */
	processor_cpu_time.update();
	// Total CPU time for all processors. 
	total_cpu_time.cput = CFDkit_Summation_MPI(processor_cpu_time.cput); 
	/* Periodically save restart solution files. */
	if (!first_step &&
	    top_level == FINEST_LEVEL &&
	    number_of_time_steps-Input_Parameters->Restart_Solution_Save_Frequency*
	    (number_of_time_steps/Input_Parameters->Restart_Solution_Save_Frequency) == 0 ) {
	  if (!batch_flag) cout << "\n\n  Saving solution to restart data file(s) after"
				<< " n = " << number_of_time_steps << " steps (iterations).";
	  error_flag = Write_Restart_Solution(Local_SolnBlk[top_level], 
					      List_of_Local_Solution_Blocks[top_level], 
					      *Input_Parameters,
					      number_of_time_steps,
					      Time,
					      processor_cpu_time);
	  if (error_flag) {
	    cout << "\n FASMultigrid ERROR: Unable to open restart output data file(s) "
		 << "on processor "
		 << List_of_Local_Solution_Blocks[top_level].ThisCPU
		 << ".\n";
	    cout.flush();
	  } /* endif */
	  error_flag = CFDkit_OR_MPI(error_flag);
	  if (error_flag) return (error_flag);
	  // cout << "\n";
	  cout.flush();
	} /* endif */
	
	/* Output progress information for the calculation. */
	
	if (!batch_flag) Output_Progress_L2norm(number_of_time_steps,
						Time,
						total_cpu_time,
						residual_l2_norm,
						first_step,
						50);
	// 	  if (!batch_flag) Output_Progress(number_of_time_steps,
	// 					   Time,
	// 					   total_cpu_time,
	// 					   residual_l1_norm,
	// 					   first_step,
	// 					   50);
	if (CFDkit_Primary_MPI_Processor() && !first_step) {
	  Output_Progress_to_File(residual_file,
				  number_of_time_steps,
				  Time,
				  total_cpu_time,
				  residual_l1_norm,
				  residual_l2_norm,
				  residual_max_norm);
	} /* endif */

	/* Check if the maximum number of time steps has been reached 
	   or if the residual has dropped below the prescribed level for
	   a FMG cycle */
	if (number_of_time_steps >= Input_Parameters->Maximum_Number_of_Time_Steps ||
	    (residual_l2_norm < Input_Parameters->Multigrid_IP.Convergence_Residual_Level && 
	     cycles != 1 && 
	     !first_step &&
	     top_level != FINEST_LEVEL)) {
	  
	  /* Output final progress information for the calculation. */

	  if (!batch_flag) Output_Progress_L2norm(number_of_time_steps,
						  Time,
						  total_cpu_time,
						  residual_l2_norm,
						  first_step,
						  number_of_time_steps);
	  
	  //  	if (!batch_flag) Output_Progress(number_of_time_steps,
	  //  					 Time,
	  //  					 total_cpu_time,
	  //  					 residual_l1_norm,
	  //  					 first_step,
	  //  					 number_of_time_steps);
	  
	  /* quit */
	  break;                    // break out of cycles_for_this_level loop
	}

	/* If residual is lower than specified, freeze limiters */
	if (!first_step &&
	    cycles != 1 &&
	    Input_Parameters->Freeze_Limiter &&
	    limiter_freezing == OFF &&
	    Input_Parameters->Limiter_Type != LIMITER_ZERO &&
	    residual_l2_norm <= Input_Parameters->Freeze_Limiter_Residual_Level) {
	  for (level = top_level; level < Input_Parameters->Multigrid_IP.Levels; level++) {
	    Freeze_Limiters(Local_SolnBlk[level], 
			    List_of_Local_Solution_Blocks[level]);
	  } /* end for */

	  limiter_freezing = ON;	
	}
	
	/* Update solution for next time step using multigrid
	   time stepping scheme. */

	error_flag = Coarse_Grid_Correction(top_level,
					    top_level);
	
	if (error_flag) return (error_flag);
	
	/* update step count */
	if (first_step) first_step = 0;
	number_of_time_steps++;
      } /* end for */                              // end of cycles_for_this_level loop

      /* Prolong solution up one level if not on the finest level yet */
      if (top_level != FINEST_LEVEL) {

	/* Prolong */
	
	error_flag = Prolong_Solution_Blocks(top_level);
	if (error_flag) return (error_flag);
     	
	/* Enforce BCs and pass messages */
	
	error_flag = Exchange_solution_information(top_level-1,
						   OFF);	
	if (error_flag) return (error_flag);
	
	BCs(Local_SolnBlk[top_level-1], 
	    List_of_Local_Solution_Blocks[top_level-1]);

	// need to set non-solution vector components
	Prescribe_NonSol(Local_SolnBlk[top_level-1],
			 List_of_Local_Solution_Blocks[top_level-1],
			 *Input_Parameters);
	

	
      } /* end if */
    } /* end for */

    if (!batch_flag) cout << "\n\n FAS Multigrid computations complete on " 
			  << Date_And_Time() << ".\n";
    
  } /* endif */

  /* Update ghostcell information and prescribe boundary conditions to ensure
     that the solution is consistent on each block. */
  
  error_flag = Exchange_solution_information(FINEST_LEVEL,OFF);
  
  if (error_flag) return (error_flag);
  
  BCs(Local_SolnBlk[FINEST_LEVEL], 
      List_of_Local_Solution_Blocks[FINEST_LEVEL]);

  // need to set non-solution vector components
  Prescribe_NonSol(Local_SolnBlk[FINEST_LEVEL],
		   List_of_Local_Solution_Blocks[FINEST_LEVEL],
		   *Input_Parameters);
  
  /* Close residual file. */
  
  error_flag = Close_Progress_File(residual_file);
  
  /********************************************************
   * Solution calculations using multigrid complete.      *
   ********************************************************/

  // DEBUG
  // plot solution on specific grid level, then exit
//   Output_Cells_Tecplot(Local_SolnBlk[0],
// 		 List_of_Local_Solution_Blocks[0],
// 		 *Input_Parameters,
// 		 number_of_time_steps,
// 		 Time);
//   exit(-1);

  return (0);  
}



#endif // _RTE2D_MULTIGRID_SPECIALIZATION_INCLUDED 
