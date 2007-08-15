/* AdvectDiffuse2DQuadMultigrid.cc:
   FAS Multigrid 2D Explicit Specialization for
   AdvectDiffuse2D - function definitions */

/* Include required C++ libraries. */

#ifdef _GNU_GCC_V3
#include <limits>
#endif

/* Include 2D advection diffusion equation quadrilateral mesh solution header file. */

#ifndef _ADVECTDIFFUSE2D_QUAD_INCLUDED
#include "AdvectDiffuse2DQuad.h"
#endif // _ADVECTDIFFUSE2D_QUAD_INCLUDED

/* Include 2D advection diffusion equation quadrilateral mesh solution header file. */

#ifndef _FASMULTIGRID2D_INCLUDED
#include "../FASMultigrid2D/FASMultigrid2D.h"
#endif // _FASMULTIGRID2D_INCLUDED

/**********************************************************************
 * FAS_Multigrid_Quad_Block::allocate --                              *
 *                                                                    *
 * Explicit specialization of the allocate routine for the FAS        *
 * Multigrid quadrilateral blocks.                                    *
 *                                                                    *
 **********************************************************************/
template <> void FAS_Multigrid_Quad_Block<double>::
allocate(const int &Ni, const int &Nj) {
  assert(Ni > 1 && Nj > 1);
  NCi = Ni; NCj = Nj;
  P = new double*[NCi];
  uo_MG = new double*[NCi];
  for (int i = 0; i < NCi; i++) {
    P[i] = new double[NCj];
    uo_MG[i] = new double[NCj];
    for (int j = 0; j < NCj; j++) {
      P[i][j] = ZERO;
      uo_MG[i][j] = ZERO;
    } 
  }
}
/**********************************************************************
 * DTS_Multigrid_Quad_Block::allocate --                              *
 *                                                                    *
 * Explicit specialization of the allocate routinie for the dual      *
 * time-stepping FAS multigrid quadrilateral blocks.                  *
 *                                                                    *
 **********************************************************************/
template <> void DTS_Multigrid_Quad_Block<double>::
allocate(const int &Ni, const int &Nj) {
  NCi = Ni; NCj = Nj;
  Un = new double*[NCi];
  Uo = new double*[NCi];
  for (int i = 0; i < NCi; i++) {
    Un[i] = new double[NCj];
    Uo[i] = new double[NCj];
    for (int j = 0; j < NCj; j++) {
      Un[i][j] = ZERO;
      Uo[i][j] = ZERO;
    }
  }
}

/**********************************************************************
 * Routine: Output_Multigrid                                          *
 *                                                                    *
 * This routine writes out the node information for each multigrid    *
 * grid level (finest to coarest) for a 1D array of 2D quadrilateral  *
 * multi-block solution blocks.                                       *
 *                                                                    *
 **********************************************************************/
template <> int FAS_Multigrid2D_Solver<double,
				       AdvectDiffuse2D_Quad_Block,
				       AdvectDiffuse2D_Input_Parameters>::
Output_Multigrid(int &number_of_time_steps,
		 double &Time) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;    

  // Determine main prefix of output data file names.
  i = 0;
  while (1) {
    if (IP->Output_File_Name[i] == ' ' ||
	IP->Output_File_Name[i] == '.') break;
    prefix[i] = IP->Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP->Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_multigrid_");

  // Output to a seperate file for each multigrid level.
  for (int level = 0; level < IP->Multigrid_IP.Levels; level++) {

    // Determine prefix of output data file names for current grid level.
    strcpy(output_file_name,prefix);
    sprintf(extension,"%.2d",level);
    strcat(output_file_name,extension);
    strcat(output_file_name,"_cpu");

    // Determine output data file name for this processor.
    sprintf(extension,"%.6d",List_of_Local_Solution_Blocks[level].ThisCPU);
    strcat(extension,".dat");
    strcat(output_file_name,extension);
    output_file_name_ptr = output_file_name;

    // Open the output data file.
    output_file.open(output_file_name_ptr,ios::out);
    if (output_file.bad()) return 1;

    // Write the solution data for each solution block.
    i_output_title = 1;
    for (int nb = 0; nb < List_of_Local_Solution_Blocks[level].Nblk; nb++) {
      if (List_of_Local_Solution_Blocks[level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	Output_Tecplot(Local_SolnBlk[level][nb],
 		       number_of_time_steps,
 		       Time,
 		       List_of_Local_Solution_Blocks[level].Block[nb].gblknum,
 		       i_output_title,
 		       output_file);
 	if (i_output_title) i_output_title = 0;
      }
    }

    // Close the output data file.
    output_file.close();

  }

  // Writing of multigrid grids data files complete.
  return 0;

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
template <> void FAS_Multigrid2D_Solver<double,
					AdvectDiffuse2D_Quad_Block,
					AdvectDiffuse2D_Input_Parameters>::
Restrict_Residuals(const int &Level_Fine) {
  
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
	  
	  MG_SolnBlk[Level_Coarse][b].P[i_coarse][j_coarse] = 
	    (Local_SolnBlk[Level_Fine][b].dudt[i_fine][j_fine][0] *
	     Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].A +
	     Local_SolnBlk[Level_Fine][b].dudt[i_fine+1][j_fine][0] *
	     Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][j_fine].A +
	     Local_SolnBlk[Level_Fine][b].dudt[i_fine][j_fine+1][0] *
	     Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine+1].A +
	     Local_SolnBlk[Level_Fine][b].dudt[i_fine+1][j_fine+1][0] *
	     Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][j_fine+1].A) /
	    Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].A;
			       
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
template <> void FAS_Multigrid2D_Solver<double,
					AdvectDiffuse2D_Quad_Block,
					AdvectDiffuse2D_Input_Parameters>::
Restrict_Boundary_Ref_States(const int &Level_Fine) {
  
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
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][JCu].A +
	   Local_SolnBlk[Level_Fine][b].UoN[i_fine+1] *
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][JCu].A) /
	  (Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][JCu].A +
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][JCu].A);

	Local_SolnBlk[Level_Coarse][b].UoS[i_coarse] = 
	  (Local_SolnBlk[Level_Fine][b].UoS[i_fine] *
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][JCl].A +
	   Local_SolnBlk[Level_Fine][b].UoS[i_fine+1] *
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][JCl].A) /
	  (Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][JCl].A +
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine+1][JCl].A);

      } /* endfor */

      /* Get lower and upper indices */
      int ICl = Local_SolnBlk[Level_Fine][b].ICl;
      int ICu = Local_SolnBlk[Level_Fine][b].ICu;

      /* Loop through j-direction cells */
      for (int j_coarse = Local_SolnBlk[Level_Coarse][b].JCl; j_coarse <= Local_SolnBlk[Level_Coarse][b].JCu; j_coarse++) {
	j_fine = 2*(j_coarse-Nghost)+Nghost;

	Local_SolnBlk[Level_Coarse][b].UoW[j_coarse] = 
	  (Local_SolnBlk[Level_Fine][b].UoW[j_fine] *
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[ICl][j_fine].A +
	   Local_SolnBlk[Level_Fine][b].UoW[j_fine+1] *
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[ICl][j_fine+1].A) /
	  (Local_SolnBlk[Level_Fine][b].Grid.Cell[ICl][j_fine].A +
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[ICl][j_fine+1].A);

	Local_SolnBlk[Level_Coarse][b].UoE[j_coarse] = 
	  (Local_SolnBlk[Level_Fine][b].UoE[j_fine] *
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[ICu][j_fine].A +
	   Local_SolnBlk[Level_Fine][b].UoE[j_fine+1] *
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[ICu][j_fine+1].A) /
	  (Local_SolnBlk[Level_Fine][b].Grid.Cell[ICu][j_fine].A +
	   Local_SolnBlk[Level_Fine][b].Grid.Cell[ICu][j_fine+1].A);

      } /* endfor */
      
    } /* end if */     
  } /* endfor */
}

/********************************************************************
 * Routine: Prolong_Solution_Blocks (for Multigrid)                 *
 *          for AdvectDiffuse2D                                     *
 * Prolong solution states stored in U[][].u from Level_Coarse      *
 * to Level_Fine for all blocks on local block list **Soln_ptr      *
 *                                                                  *
 ********************************************************************/
template <> int FAS_Multigrid2D_Solver<double,
				       AdvectDiffuse2D_Quad_Block,
				       AdvectDiffuse2D_Input_Parameters>::
Prolong_Solution_Blocks(const int &Level_Coarse) {

  int error_flag, i_fine,j_fine;
  double solution;
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

	  // Set injection flag if prolongation by injection is desired.
	  if (IP->Multigrid_IP.Prolong_Using_Injection) injection_for_this_coarse_cell = true;

	restart: ;
	    
	  /*******************************************************************/
	  /*** Calculate cell-centred solution at SW Corner of coarse cell ***/
	  /*******************************************************************/
	  i_fine = 2*(i_coarse-Nghost)+Nghost; // SW Corner i index (fine)
	  j_fine = 2*(j_coarse-Nghost)+Nghost; // SW Corner j index (fine)

	  if (!injection_for_this_coarse_cell) {
	    error_flag = Bilinear_Interpolation(
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1].u,
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse].u,
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1].u,
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
	       Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc,
	       solution);
	    
	    /* If fine grid cell-centre falls outside of quad in SW dir
	       formed by coarse grid cell-centres, try NW one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
	    }
	    /* If fine grid cell-centre falls outside of quad in NW dir
	       formed by coarse grid cell-centres, try SE one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
	    }
	    /* If fine grid cell-centre falls outside of quad in SE dir
	       formed by coarse grid cell-centres, try NE one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
	    }
	    if (error_flag != 0) {
	      injection_for_this_coarse_cell = true;
	      goto restart;
	    }
	  } else {
	    // use simple injection instead of bilinear interpolation
	    solution = Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u;
	    error_flag = 0;
	  }
	  
	  Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine].u = solution;
	    
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
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1].u,
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse].u,
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1].u,
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
	       Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
	       solution);
	    
	    /* If fine grid cell-centre falls outside of quad in SE dir
	       formed by coarse grid cell-centres, try NE one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
	    }
	    /* If fine grid cell-centre falls outside of quad in NE dir
	       formed by coarse grid cell-centres, try SW one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
	    }
	    /* If fine grid cell-centre falls outside of quad in SW dir
	       formed by coarse grid cell-centres, try NW one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
	    }
	    if (error_flag != 0) {
	      injection_for_this_coarse_cell = true;
	      goto restart;
	    }
	  } else {
	    // use simple injection instead of bilinear interpolation
	    solution = Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u;
	    error_flag = 0;
	  }
	  
	  Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine].u = solution;
	  
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
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1].u,
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1].u,
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse].u,
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
	       Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
	       solution);
	    
	    /* If fine grid cell-centre falls outside of quad in NE dir
	       formed by coarse grid cell-centres, try SE one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
	    }
	    /* If fine grid cell-centre falls outside of quad in SE dir
	       formed by coarse grid cell-centres, try NW one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
	    }
	    /* If fine grid cell-centre falls outside of quad in NW dir
	       formed by coarse grid cell-centres, try SW one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
	    }
	    if (error_flag != 0) {
	      injection_for_this_coarse_cell = true;
	      goto restart;
	    }
	  } else {
	    // use simple injection instead of bilinear interpolation
	    solution = Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u;
	    error_flag = 0;
	  }
	  
	  Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine].u = solution;
	  
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
               Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse].u,
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1].u,
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1].u,
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
	       Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
	       Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
	       Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
	       solution);
	    
	    /* If fine grid cell-centre falls outside of quad in NW dir
	       formed by coarse grid cell-centres, try NE one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
	    }
	    /* If fine grid cell-centre falls outside of quad in NE dir
	       formed by coarse grid cell-centres, try SW one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
	    }
	    /* If fine grid cell-centre falls outside of quad in SW dir
	       formed by coarse grid cell-centres, try SE one */
	    if (error_flag != 0) {
	      error_flag = Bilinear_Interpolation(
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 solution);
	    }
	    if (error_flag != 0) {
	      injection_for_this_coarse_cell = true;
	      goto restart;
	    }
	  } else {
	    // use simple injection instead of bilinear interpolation
	    solution = Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u;
	    error_flag = 0;
	  }
	
	  Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine].u = solution;

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
  return 0;
}

/********************************************************************
 * Routine: Prolong_and_Update_Solution_Blocks (for Multigrid)      *
 *          for AdvectDiffuse2D                                     *
 * Prolong solution changes stored in U[][].u from Level_Coarse     *
 * to Level_Fine for all blocks on local block list **Soln_ptr and  *
 * add the changes to the solution U[][].u on Level_Fine blocks     *
 * (modifies any value present in U[][] on Level_Fine)              *
 *                                                                  *
 ********************************************************************/
template <> int FAS_Multigrid2D_Solver<double,
			   AdvectDiffuse2D_Quad_Block,
			   AdvectDiffuse2D_Input_Parameters>::
Prolong_and_Update_Solution_Blocks(const int &Level_Coarse) {

  int error_flag, i_fine,j_fine;
  double sol_change;
  int Level_Fine = Level_Coarse - 1;
  bool injection_for_this_coarse_cell;
  injection_for_this_coarse_cell = false;

  // loop through used blocks only
  for (int b = 0 ; b <= List_of_Local_Solution_Blocks[Level_Coarse].Nblk-1 ; ++b ) {
    if (List_of_Local_Solution_Blocks[Level_Coarse].Block[b].used == ADAPTIVEBLOCK2D_USED) { 

      /* Get number of ghost cells */
      int Nghost = Local_SolnBlk[Level_Coarse][b].Nghost;

      /* Enforce Neumann boundary conditions as appropriate */
      if (!IP->Multigrid_IP.Prolong_Using_Injection) {
	for (int j_coarse = Local_SolnBlk[Level_Coarse][b].JCl; j_coarse <= Local_SolnBlk[Level_Coarse][b].JCu; j_coarse++) {

	  // West Boundary
	  if (Local_SolnBlk[Level_Coarse][b].Grid.BCtypeW[j_coarse] == BC_NEUMANN) {

	    Local_SolnBlk[Level_Coarse][b].U[Local_SolnBlk[Level_Coarse][b].ICl-1][j_coarse].u =
	      Local_SolnBlk[Level_Coarse][b].U[Local_SolnBlk[Level_Coarse][b].ICl][j_coarse].u;
	  } /* end if */
	  // East Boundary
	  if (Local_SolnBlk[Level_Coarse][b].Grid.BCtypeE[j_coarse] == BC_NEUMANN) {
	    Local_SolnBlk[Level_Coarse][b].U[Local_SolnBlk[Level_Coarse][b].ICu+1][j_coarse].u =
	      Local_SolnBlk[Level_Coarse][b].U[Local_SolnBlk[Level_Coarse][b].ICu][j_coarse].u;
	  } /* end if */
	} /* end for */
	for (int i_coarse = Local_SolnBlk[Level_Coarse][b].ICl; i_coarse <= Local_SolnBlk[Level_Coarse][b].ICu; i_coarse++) {
	  // South Boundary
	  if (Local_SolnBlk[Level_Coarse][b].Grid.BCtypeS[i_coarse] == BC_NEUMANN) {
	    
	    Local_SolnBlk[Level_Coarse][b].U[i_coarse][Local_SolnBlk[Level_Coarse][b].JCl-1].u =
	      Local_SolnBlk[Level_Coarse][b].U[i_coarse][Local_SolnBlk[Level_Coarse][b].JCl].u;
	  } /* end if */
	  // North Boundary
	  if (Local_SolnBlk[Level_Coarse][b].Grid.BCtypeN[i_coarse] == BC_NEUMANN) {
	    Local_SolnBlk[Level_Coarse][b].U[i_coarse][Local_SolnBlk[Level_Coarse][b].JCu+1].u =
	      Local_SolnBlk[Level_Coarse][b].U[i_coarse][Local_SolnBlk[Level_Coarse][b].JCu].u;
	  } /* end if */
	} /* end for */
      } /* end if */

      /* Loop through coarse cells */
      for (int i_coarse = Local_SolnBlk[Level_Coarse][b].ICl; i_coarse <= Local_SolnBlk[Level_Coarse][b].ICu; i_coarse++) {
	for (int j_coarse = Local_SolnBlk[Level_Coarse][b].JCl; j_coarse <= Local_SolnBlk[Level_Coarse][b].JCu; j_coarse++) {

	  /* Check if coarse grid cell is on a Dirichlet Boundary.  If so, use
	     simple injection instead of bilinear interpolation for all
	     4 related fine cells. */
	  if ((i_coarse == Local_SolnBlk[Level_Coarse][b].ICl &&
	       Local_SolnBlk[Level_Coarse][b].Grid.BCtypeW[j_coarse] == BC_DIRICHLET) ||
	      (i_coarse == Local_SolnBlk[Level_Coarse][b].ICu &&
	       Local_SolnBlk[Level_Coarse][b].Grid.BCtypeE[j_coarse] == BC_DIRICHLET) ||
	      (j_coarse == Local_SolnBlk[Level_Coarse][b].JCl &&
	       Local_SolnBlk[Level_Coarse][b].Grid.BCtypeS[i_coarse] == BC_DIRICHLET) ||
	      (j_coarse == Local_SolnBlk[Level_Coarse][b].JCu &&
	       Local_SolnBlk[Level_Coarse][b].Grid.BCtypeN[i_coarse] == BC_DIRICHLET)) {
	    
	    /* Find i,j indices of the SW Corner of coarse cell */
	    i_fine = 2*(i_coarse-Nghost)+Nghost; // SW Corner i index (fine)
	    j_fine = 2*(j_coarse-Nghost)+Nghost; // SW Corner j index (fine)
	    
	    sol_change = Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u;
	    Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine].u += sol_change;
	    Local_SolnBlk[Level_Fine][b].U[i_fine+1][j_fine].u += sol_change;
	    Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine+1].u += sol_change;
	    Local_SolnBlk[Level_Fine][b].U[i_fine+1][j_fine+1].u += sol_change;
	    error_flag = 0;	    	    

	  } else {

	    // Set injection flag if prolongation by injection is desired.
	    if (IP->Multigrid_IP.Prolong_Using_Injection) injection_for_this_coarse_cell = true;

	  restart: ;

	    /*******************************************************************/
	    /*** Calculate cell-centred solution at SW Corner of coarse cell ***/
	    /*******************************************************************/
	    i_fine = 2*(i_coarse-Nghost)+Nghost; // SW Corner i index (fine)
	    j_fine = 2*(j_coarse-Nghost)+Nghost; // SW Corner j index (fine)
	    
	    if (!injection_for_this_coarse_cell) {
	      error_flag = Bilinear_Interpolation(
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc,
		 sol_change);
	      
	      /* If fine grid cell-centre falls outside of quad in SW dir
		 formed by coarse grid cell-centres, try NW one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
	      }
	      /* If fine grid cell-centre falls outside of quad in NW dir
		 formed by coarse grid cell-centres, try SE one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
	      }
	      /* If fine grid cell-centre falls outside of quad in SE dir
		 formed by coarse grid cell-centres, try NE one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
	      }
	      if (error_flag != 0) {
		injection_for_this_coarse_cell = true;
		goto restart;
	      }
	    } else {
	      // use simple injection instead of bilinear interpolation
	      sol_change = Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u;
	      error_flag = 0;
	    }
	    Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine].u += sol_change;
	    
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
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 sol_change);
	      
	      /* If fine grid cell-centre falls outside of quad in SE dir
		 formed by coarse grid cell-centres, try NE one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
	      }
	      /* If fine grid cell-centre falls outside of quad in NE dir
		 formed by coarse grid cell-centres, try SW one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
	      }
	      /* If fine grid cell-centre falls outside of quad in SW dir
		 formed by coarse grid cell-centres, try NW one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
	      }
	      if (error_flag != 0) {
		injection_for_this_coarse_cell = true;
		goto restart;
	      }
	    } else {
	      // use simple injection instead of bilinear interpolation
	      sol_change = Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u;
	      error_flag = 0;
	    }
	    Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine].u += sol_change;
	    
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
	         Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 sol_change);
	      
	      /* If fine grid cell-centre falls outside of quad in NE dir
		 formed by coarse grid cell-centres, try SE one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
	      }
	      /* If fine grid cell-centre falls outside of quad in SE dir
		 formed by coarse grid cell-centres, try NW one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
	      }
	      /* If fine grid cell-centre falls outside of quad in NW dir
		 formed by coarse grid cell-centres, try SW one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
	      }
	      if (error_flag != 0) {
		injection_for_this_coarse_cell = true;
		goto restart;
	      }
	    } else {
	      // use simple injection instead of bilinear interpolation
	      sol_change = Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u;
	      error_flag = 0;
	    }
	    Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine].u += sol_change;
	    
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
                 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse+1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		 Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		 Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		 Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		 sol_change);
	      
	      /* If fine grid cell-centre falls outside of quad in NW dir
		 formed by coarse grid cell-centres, try NE one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse+1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse+1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
	      }
	      /* If fine grid cell-centre falls outside of quad in NE dir
		 formed by coarse grid cell-centres, try SW one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse-1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse-1][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse-1][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
	      }
	      /* If fine grid cell-centre falls outside of quad in SW dir
		 formed by coarse grid cell-centres, try SE one */
	      if (error_flag != 0) {
		error_flag = Bilinear_Interpolation(
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse-1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse].Xc,
		   Local_SolnBlk[Level_Coarse][b].U[i_coarse+1][j_coarse-1].u,
		   Local_SolnBlk[Level_Coarse][b].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
		   Local_SolnBlk[Level_Fine][b].Grid.Cell[i_fine][j_fine].Xc, 
		   sol_change);
	      }
	      if (error_flag != 0) {
		injection_for_this_coarse_cell = true;
		goto restart;
	      }
	    } else {
	      // use simple injection instead of bilinear interpolation
	      sol_change = Local_SolnBlk[Level_Coarse][b].U[i_coarse][j_coarse].u;
	      error_flag = 0;
	    }
	    Local_SolnBlk[Level_Fine][b].U[i_fine][j_fine].u += sol_change;
	    
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
  return 0;
}

/********************************************************
 * Routine: Store_Current_Solution_in_uo_MG             *
 *          for AdvectDiffuse2D                         *
 * This routine copies the solution state in each block *
 * into the array uo_MG                                 *
 *                                                      *
 ********************************************************/
template <> void FAS_Multigrid2D_Solver<double,
					AdvectDiffuse2D_Quad_Block,
					AdvectDiffuse2D_Input_Parameters>::
Store_Current_Solution_in_uo_MG(const int &Level) {

  /* Copy for each block */

  for (int b = 0 ; b <= List_of_Local_Solution_Blocks[Level].Nblk-1 ; ++b ) {
    if (List_of_Local_Solution_Blocks[Level].Block[b].used == ADAPTIVEBLOCK2D_USED) {
      int Ng = Local_SolnBlk[Level][b].Nghost;
      for (int i = Local_SolnBlk[Level][b].ICl-Ng; i <= Local_SolnBlk[Level][b].ICu+Ng; i++) {
	for (int j = Local_SolnBlk[Level][b].JCl-Ng; j <= Local_SolnBlk[Level][b].JCu+Ng; j++) {
	  MG_SolnBlk[Level][b].uo_MG[i][j] = Local_SolnBlk[Level][b].U[i][j].u;
	}
      }
    } /* endif */
  }  /* endfor */
}

/*********************************************************
 * Routine: Subtract_dUdt_from_P                         *
 *                                                       *
 * This routine subtracts dUdt from P.  P is expected to *
 * contain the restricted fine grid residuals, where     *
 * dUdt is expected to contain the residuals evaluated   *
 * from the restricted fine grid solution, such that     *
 * after this operation the forcing term P is formed.    *
 *********************************************************/
template <> void FAS_Multigrid2D_Solver<double,
					AdvectDiffuse2D_Quad_Block,
					AdvectDiffuse2D_Input_Parameters>::
Subtract_dUdt_from_P(const int &Level) {

  /* Loop through each block */

  for (int b = 0 ; b <= List_of_Local_Solution_Blocks[Level].Nblk-1 ; ++b ) {
    if (List_of_Local_Solution_Blocks[Level].Block[b].used == ADAPTIVEBLOCK2D_USED) {
      for (int i = Local_SolnBlk[Level][b].ICl; i <= Local_SolnBlk[Level][b].ICu; i++) {
	for (int j = Local_SolnBlk[Level][b].JCl; j <= Local_SolnBlk[Level][b].JCu; j++) {
	  MG_SolnBlk[Level][b].P[i][j] -= Local_SolnBlk[Level][b].dudt[i][j][0];
	} /* endfor */
      } /* endfor */
    } /* endif */
  }  /* endfor */
}

/**********************************************************
 * Routine: Update_Primitive_Variables                    *
 *                                                        *
 * This routine updates all primitive variables W using   *
 * the conserved variables U, for all solution blocks on  *
 * a level.                                               *
 **********************************************************/
template <> void FAS_Multigrid2D_Solver<double,
                                        AdvectDiffuse2D_Quad_Block,
                                        AdvectDiffuse2D_Input_Parameters>::
Update_Primitive_Variables(const int &Level) {
  /* Doesn't Apply for AdvectDiffuse2D */
  /* Leave body of function blank */
}

/**********************************************************
 * Routine: Evaluate_Solution_Changes                     *
 *                                                        *
 * This routine subtracts uo_MG[][] from U[][].u          *
 * and store it into U[][].u for later prolongation       *
 * to the finer level.                                    *
 **********************************************************/
template <> void FAS_Multigrid2D_Solver<double,
					AdvectDiffuse2D_Quad_Block,
					AdvectDiffuse2D_Input_Parameters>::
Evaluate_Solution_Changes(const int &Level) {

  /* Loop through each block */

  for (int b = 0 ; b <= List_of_Local_Solution_Blocks[Level].Nblk-1 ; ++b ) {
    if (List_of_Local_Solution_Blocks[Level].Block[b].used == ADAPTIVEBLOCK2D_USED) {
      int Ng = Local_SolnBlk[Level][b].Nghost;
      for (int i = Local_SolnBlk[Level][b].ICl-Ng; i <= Local_SolnBlk[Level][b].ICu+Ng; i++) {
	for (int j = Local_SolnBlk[Level][b].JCl-Ng; j <= Local_SolnBlk[Level][b].JCu+Ng; j++) {
	  Local_SolnBlk[Level][b].U[i][j].u -= MG_SolnBlk[Level][b].uo_MG[i][j];
	}
      }
    } /* endif */
  }  /* endfor */
}

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
template <> int FAS_Multigrid2D_Solver<double,
				       AdvectDiffuse2D_Quad_Block,
				       AdvectDiffuse2D_Input_Parameters>::
dUdt_Multistage_Explicit_for_Multigrid(const int &Level,
				       const int &Top_Level,
				       const int &i_stage) {

  int error_flag, limiter_type, k_residual;
  
  /* Force use of first-order spatial reconstruction on coarse levels. */

  limiter_type = IP->i_Limiter;
  if (Level > Top_Level && IP->Multigrid_IP.First_Order_Coarse_Mesh_Reconstruction) {
    IP->i_Limiter = LIMITER_ZERO;
  }

  /* Temporarily overwrite IP->i_Time_Integration
     for using Euler2DQuadMultiBlock function */
  IP->i_Time_Integration = IP->Multigrid_IP.i_Smoothing;
  
  /* Evaluate the residual for each solution block. */

  error_flag = dUdt_Multistage_Explicit(Local_SolnBlk[Level],
					*List_of_Global_Solution_Blocks,
					List_of_Local_Solution_Blocks[Level],
					*IP,
					i_stage);
  if (error_flag) return (error_flag);

  // Apply the defect correction forcing term, if necessary.
  //error_flag = Apply_the_FAS_Multigrid_Forcing_Term(Level,
  //					    Top_Level,
  //					    i_stage);
  //if (error_flag) return error_flag;

  /* Reset limiter, and time integration parameters. */

  if (Level > Top_Level && IP->Multigrid_IP.First_Order_Coarse_Mesh_Reconstruction) {
    IP->i_Limiter = limiter_type;
  }
  IP->i_Time_Integration = TIME_STEPPING_MULTIGRID;

  /* Quadrilateral multi-block solution residuals
     successfully calculated.  Return. */

  return 0;

}

/**********************************************************************
 * Routine: Apply_the_FAS_Multigrid_Forcing_Term                      *
 *                                                                    *
 * This routine applies the FAS multigrid defect correction forcing   *
 * term, P, for a 1D array of 2D quadrilateral solution blocks.       *
 *                                                                    *
 **********************************************************************/
template <> int FAS_Multigrid2D_Solver<double,
				       AdvectDiffuse2D_Quad_Block,
				       AdvectDiffuse2D_Input_Parameters>::
Apply_the_FAS_Multigrid_Forcing_Term(const int &Level,
				     const int &Top_Level,
				     const int &i_stage) {

  int k_residual;

  // Do not apply the forcing term if the current level is finer or
  // equal to the current top level.  Exit immediately.
  if (Level <= Top_Level) return 0;

  // Temporarily overwrite the time integration type to the multistage
  // optimally smoothing scheme.
  IP->i_Time_Integration = IP->Multigrid_IP.i_Smoothing;

  /* Add the defect correction forcing term as necessary. */

  if (Level > Top_Level) {
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
	     Local_SolnBlk[Level][b].dudt[i][j][k_residual] += 
	       (IP->CFL_Number*Local_SolnBlk[Level][b].dt[i][j])*
	       MG_SolnBlk[Level][b].P[i][j];
	   } /* endfor */
         } /* endfor */
       } /* endif */
     }  /* endfor */
  } /* endif */

  // Reset the time-integration type.
  IP->i_Time_Integration = TIME_STEPPING_MULTIGRID;

  // The forcing term has been applied successfully.
  return 0;

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
template <> int FAS_Multigrid2D_Solver<double,
				       AdvectDiffuse2D_Quad_Block,
				       AdvectDiffuse2D_Input_Parameters>::
Residual_Evaluation_for_Multigrid(const int &Level,
				  const int &Top_Level,
				  const double &dt,
				  const int &apply_forcing_term_flag) {
  
  int error_flag, limiter_type;
  
  /* Force use of first-order spatial reconstruction on coarse levels. */

  limiter_type = IP->i_Limiter;
  if (Level > Top_Level && IP->Multigrid_IP.First_Order_Coarse_Mesh_Reconstruction) {
    IP->i_Limiter = LIMITER_ZERO;
  }

  /* Evaluate the residual for each solution block. */
  
  for (int b = 0 ; b <= List_of_Local_Solution_Blocks[Level].Nblk-1 ; ++b ) {
    if (List_of_Local_Solution_Blocks[Level].Block[b].used == ADAPTIVEBLOCK2D_USED) {
      error_flag = dUdt_Residual_Evaluation(Local_SolnBlk[Level][b],
					    *IP);
      if (error_flag) return (error_flag);
    } /* endif */
  }  /* endfor */
  
  /* Reset limiter. */

  if (Level > Top_Level && IP->Multigrid_IP.First_Order_Coarse_Mesh_Reconstruction) {
    IP->i_Limiter = limiter_type;
  }

  // Add dual-time-stepping physical-time source term if required.
  if (IP->i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
    for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
      if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	for (int j = Local_SolnBlk[Level][nb].JCl; j <= Local_SolnBlk[Level][nb].JCu; j++) {
	  for (int i = Local_SolnBlk[Level][nb].ICl; i <= Local_SolnBlk[Level][nb].ICu; i++) {
	    if (IP->Multigrid_IP.i_Physical_Time_Integration == TIME_STEPPING_IMPLICIT_EULER) {
	      Local_SolnBlk[Level][nb].dudt[i][j][0] -= (Local_SolnBlk[Level][nb].U[i][j].u -
							 DTS_SolnBlk[Level][nb].Un[i][j])/
		                                        (IP->Multigrid_IP.Physical_Time_CFL_Number*dt);
	    } else if (IP->Multigrid_IP.i_Physical_Time_Integration == TIME_STEPPING_IMPLICIT_SECOND_ORDER_BACKWARD) {
	      Local_SolnBlk[Level][nb].dudt[i][j][0] -= (3.0*Local_SolnBlk[Level][nb].U[i][j].u -
							 4.0*DTS_SolnBlk[Level][nb].Un[i][j] +
							 DTS_SolnBlk[Level][nb].Uo[i][j])/
		                                        (TWO*IP->Multigrid_IP.Physical_Time_CFL_Number*dt);
	    }
	  }
	}
      }
    }
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
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return (error_flag);

  // Apply boundary flux corrections to ensure that method is conservative.
  Apply_Boundary_Flux_Corrections(Local_SolnBlk[Level], 
				  List_of_Local_Solution_Blocks[Level]);

  // Add the forcing term if required.
  if (apply_forcing_term_flag) {
    for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
      if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	for (int j = Local_SolnBlk[Level][nb].JCl; j <= Local_SolnBlk[Level][nb].JCu; j++) {
	  for (int i = Local_SolnBlk[Level][nb].ICl; i <= Local_SolnBlk[Level][nb].ICu; i++) {
	    Local_SolnBlk[Level][nb].dudt[i][j][0] += MG_SolnBlk[Level][nb].P[i][j];
	  }
	}
      }
    }
  }

  /* Residuals successfully calculated for quadrilateral multi-block solution blocks
     Return. */

  return 0;

}

/********************************************************
 * Routine: Execute                                     *
 *                                                      *
 * Commence computation using FAS multigrid method      *
 *                                                      *
 ********************************************************/
template <> int FAS_Multigrid2D_Solver<double,
				       AdvectDiffuse2D_Quad_Block,
				       AdvectDiffuse2D_Input_Parameters>::
Execute(const int &batch_flag,
	int &number_of_time_steps,
	double &Time,
	CPUTime &processor_cpu_time,
	CPUTime &total_cpu_time,
	ofstream &residual_file) {

  /* Other local solution variables. */

  int error_flag, first_step, line_number, command_flag, level, limiter_freezing;
  double residual_l1_norm, residual_l2_norm, residual_max_norm;
  unsigned long cycles_for_this_level, max_cycles_for_this_level;
#ifdef _GNU_GCC_V3
  max_cycles_for_this_level = numeric_limits<unsigned long>::max();
#else
  max_cycles_for_this_level = 10000000;
#endif

  /********************************************************  
   * Initialize solution variables.                       *
   ********************************************************/

  // Send solution information between neighbouring blocks to complete
  // prescription of initial data.
  for (int level = 1; level < IP->Multigrid_IP.Levels; level++) {
    error_flag = Exchange_Solution_Information(level,ON);
    if (error_flag) return error_flag;
  }

  /* Initialize all coarser levels of solution blocks with ICs */

  for (level = 1; level < IP->Multigrid_IP.Levels; level++) {
    ICs(Local_SolnBlk[level], 
	List_of_Local_Solution_Blocks[level], 
	*IP,
	*QuadTree);

    Set_Advection_Velocity_Field(Local_SolnBlk[level],
				 List_of_Local_Solution_Blocks[level],
				 *IP);

    Restrict_Boundary_Ref_States(level-1);
  }

  /* Send solution information between neighbouring blocks to complete
     prescription of initial data. */

  for (level = 1; level < IP->Multigrid_IP.Levels; level++) {

    error_flag = Exchange_Solution_Information(level,OFF);

    if (error_flag) return (error_flag);
  } /* end for */

  /********************************************************  
   * Solve IBVP or BVP on multi-block solution-adaptive   *
   * quadrilateral mesh using Multigrid for time-marching *
   ********************************************************/

  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

  /* Open residual file and reset the CPU time. */

  first_step = 1;
  limiter_freezing = OFF;

  if (CFFC_Primary_MPI_Processor()) {
    error_flag = Open_Progress_File(residual_file,
				    IP->Output_File_Name,
				    number_of_time_steps);
    if (error_flag) {
      cout << "\n FASMultigrid ERROR: Unable to open residual file for calculation.\n";
      cout.flush();
    } /* endif */
  } /* endif */

  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
  CFFC_Broadcast_MPI(&error_flag, 1);
  if (error_flag) return (error_flag);

  processor_cpu_time.reset();

  /* Perform required number of iterations (time steps). */

  if (!IP->Time_Accurate &&
      IP->Maximum_Number_of_Time_Steps > 0) {

    if (!batch_flag) cout << "\n\n Beginning FAS Multigrid computations on "
			  << Date_And_Time() << ".\n\n";

    /* Perform the required number of Full & Regular Multigrid cycles */

    /* Determine if Full Multigrid is required */
    int initial_top_level = FINEST_LEVEL;
    if (IP->Multigrid_IP.Number_of_Cycles_per_Stage_for_Full_Multigrid > 0) {
      if (!batch_flag) cout << " Perform Full multigrid cycles\n\n";
      initial_top_level = IP->Multigrid_IP.Levels-2;
    }

    /* Loop through each top_level (for FMG) */
    for (int top_level = initial_top_level; top_level >= FINEST_LEVEL; top_level--) {

      if (!batch_flag && top_level == FINEST_LEVEL) 
	cout << "\n\n Perform Regular multigrid cycles \n  ";

      /* Unfreeze limiters if frozen */
      if (!first_step &&
	  IP->Freeze_Limiter &&
	  limiter_freezing == ON &&
	  IP->Limiter_Type != LIMITER_ZERO) {
	for (int level = top_level; level < IP->Multigrid_IP.Levels; level++) {
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

      if (IP->Multigrid_IP.Number_of_Cycles_per_Stage_for_Full_Multigrid > 0) {
	if (top_level != FINEST_LEVEL) {
	  // The number of full multigrid cycles is specified in the input
	  // parameters.
	  cycles_for_this_level = IP->Multigrid_IP.Number_of_Cycles_per_Stage_for_Full_Multigrid;
	} else {
	  // Set the number of cycles to an arbitrarily large number for
	  // the regualar FAS multigrid calculations.
	  cycles_for_this_level = max_cycles_for_this_level;
	} /* end if */
      } /* end if */

      for (int cycles = 1; cycles <= cycles_for_this_level; cycles++) {

	/* Determine the L1, L2, and max norms of the solution residual. */
	residual_l1_norm = L1_Norm_Residual(Local_SolnBlk[top_level], 
					    List_of_Local_Solution_Blocks[top_level]);
	residual_l1_norm = CFFC_Summation_MPI(residual_l1_norm); // L1 norm for all processors.
	
	residual_l2_norm = L2_Norm_Residual(Local_SolnBlk[top_level], 
					    List_of_Local_Solution_Blocks[top_level]);
	residual_l2_norm = sqr(residual_l2_norm);
	residual_l2_norm = CFFC_Summation_MPI(residual_l2_norm); // L2 norm for all processors.
	residual_l2_norm = sqrt(residual_l2_norm);
	
	residual_max_norm = Max_Norm_Residual(Local_SolnBlk[top_level], 
					      List_of_Local_Solution_Blocks[top_level]);
	residual_max_norm = CFFC_Maximum_MPI(residual_max_norm); // Max norm for all processors.
	
	/* Update CPU time used for the calculation so far. */
	processor_cpu_time.update();
	// Total CPU time for all processors. 
	total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput); 
	/* Periodically save restart solution files. */
	if (!first_step &&
	    top_level == FINEST_LEVEL &&
	    number_of_time_steps-IP->Restart_Solution_Save_Frequency*
	    (number_of_time_steps/IP->Restart_Solution_Save_Frequency) == 0 ) {
	  if (!batch_flag) cout << "\n\n  Saving solution to restart data file(s) after"
				<< " n = " << number_of_time_steps << " steps (iterations).";
	  error_flag = Write_Restart_Solution(Local_SolnBlk[top_level], 
					      List_of_Local_Solution_Blocks[top_level], 
					      *IP,
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
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return (error_flag);
	  cout << "\n";
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
	if (CFFC_Primary_MPI_Processor() && !first_step) {
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
	if (number_of_time_steps >= IP->Maximum_Number_of_Time_Steps ||
	    (residual_l2_norm < IP->Multigrid_IP.Convergence_Residual_Level && 
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
	  break;
	}

	/* If residual is lower than specified, freeze limiters */
	if (!first_step &&
	    cycles != 1 &&
	    IP->Freeze_Limiter &&
	    limiter_freezing == OFF &&
	    IP->Limiter_Type != LIMITER_ZERO &&
	    residual_l2_norm <= IP->Freeze_Limiter_Residual_Level) {
	  for (int level = top_level; level < IP->Multigrid_IP.Levels; level++) {
	    Freeze_Limiters(Local_SolnBlk[level], 
			    List_of_Local_Solution_Blocks[level]);
	  } /* end for */

	  limiter_freezing = ON;	
	}
	
	/* Update solution for next time step using multigrid
	   time stepping scheme. */
	
	error_flag = Coarse_Grid_Correction(top_level,
					    top_level,
					    ZERO);
	
	if (error_flag) return (error_flag);
	
	/* update step count */
	if (first_step) first_step = 0;
	number_of_time_steps++;
      } /* end for */

      /* Prolong solution up one level if not on the finest level yet */
      if (top_level != FINEST_LEVEL) {

	/* Prolong */
	
	error_flag = Prolong_Solution_Blocks(top_level);
	if (error_flag) return (error_flag);
	
	/* Enforce BCs and pass messages */
	
	error_flag = Exchange_Solution_Information(top_level-1,
						   OFF);	
	if (error_flag) return (error_flag);
	
	BCs(Local_SolnBlk[top_level-1], 
	    List_of_Local_Solution_Blocks[top_level-1],
	    *IP);

      } /* end if */
    } /* end for */
    
    if (!batch_flag) cout << "\n\n FAS Multigrid computations complete on " 
			  << Date_And_Time() << ".\n";
    
  } /* endif */

  /* Update ghostcell information and prescribe boundary conditions to ensure
     that the solution is consistent on each block. */
  
  error_flag = Exchange_Solution_Information(FINEST_LEVEL,OFF);
  
  if (error_flag) return (error_flag);
  
  BCs(Local_SolnBlk[FINEST_LEVEL], 
      List_of_Local_Solution_Blocks[FINEST_LEVEL],
      *IP);

  /* Close residual file. */
  
  error_flag = Close_Progress_File(residual_file);
  
  /********************************************************
   * Solution calculations using multigrid complete.      *
   ********************************************************/
  
  return 0;

}

/**********************************************************************
 **********************************************************************
 *****       Routines required for dual-time-stepping.       **********
 **********************************************************************
 **********************************************************************/

/**********************************************************************
 * FAS_Multigrid2D_Solver::dUdtau_Multistage_Explicit --              *
 *                                                                    *
 **********************************************************************/
template <> int FAS_Multigrid2D_Solver<double,
				       AdvectDiffuse2D_Quad_Block,
				       AdvectDiffuse2D_Input_Parameters>::
dUdtau_Multistage_Explicit(const int &Level,
			   const int &i_stage,
			   const double &dt) {

  // Exit immediately if dual-time-stepping is not required.
  if (IP->i_Time_Integration != TIME_STEPPING_DUAL_TIME_STEPPING) return 0;

  int error_flag, k_residual, time_integration;

  // Temporarily overwrite the time integration type to the multistage
  // optimally smoothing scheme.
  time_integration = IP->i_Time_Integration;
  IP->i_Time_Integration = IP->Multigrid_IP.i_Smoothing;

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
    break;
  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
    k_residual = 0;
    break;
  default:
    k_residual = 0;
    break;
  }

  // Evaluate the solution residual for each block.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = Local_SolnBlk[Level][nb].JCl; j <= Local_SolnBlk[Level][nb].JCu; j++) {
	for (int i = Local_SolnBlk[Level][nb].ICl; i <= Local_SolnBlk[Level][nb].ICu; i++) {
	  if (IP->Multigrid_IP.i_Physical_Time_Integration == TIME_STEPPING_IMPLICIT_EULER) {
	    Local_SolnBlk[Level][nb].dudt[i][j][k_residual] -= (IP->CFL_Number*Local_SolnBlk[Level][nb].dt[i][j])*(Local_SolnBlk[Level][nb].U[i][j].u -
														   DTS_SolnBlk[Level][nb].Un[i][j])/
		                                                                                                  (IP->Multigrid_IP.Physical_Time_CFL_Number*dt);
	  } else if (IP->Multigrid_IP.i_Physical_Time_Integration == TIME_STEPPING_IMPLICIT_SECOND_ORDER_BACKWARD) {
	    Local_SolnBlk[Level][nb].dudt[i][j][k_residual] -= (IP->CFL_Number*Local_SolnBlk[Level][nb].dt[i][j])*(3.0*Local_SolnBlk[Level][nb].U[i][j].u -
														   4.0*DTS_SolnBlk[Level][nb].Un[i][j] +
														   DTS_SolnBlk[Level][nb].Uo[i][j])/
		                                                                                                  (TWO*IP->Multigrid_IP.Physical_Time_CFL_Number*dt);
	  }
	}
      }
    }
  }

  // Reset the time-integration type.
  IP->i_Time_Integration = time_integration;

  // Solution residual successfully determined.
  return 0;

}

/**********************************************************************
 * FAS_Multigrid2D_Solver::Apply_Melson_Time_Step --                  *
 *                                                                    *
 **********************************************************************/
template <> int FAS_Multigrid2D_Solver<double,
				       AdvectDiffuse2D_Quad_Block,
				       AdvectDiffuse2D_Input_Parameters>::
Apply_Melson_Time_Step(const int &Level,
		       const int &i_stage,
		       const double &dt) {

  // Exit immediately if dual-time-stepping is not required.
  if (IP->i_Time_Integration != TIME_STEPPING_DUAL_TIME_STEPPING) return 0;

  int error_flag, k_residual, time_integration;

  // Temporarily overwrite the time integration type to the multistage
  // optimally smoothing scheme.
  time_integration = IP->i_Time_Integration;
  IP->i_Time_Integration = IP->Multigrid_IP.i_Smoothing;

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
    break;
  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
    k_residual = 0;
    break;
  default:
    k_residual = 0;
    break;
  }

  // Evaluate the solution residual for each block.
  for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = Local_SolnBlk[Level][nb].JCl; j <= Local_SolnBlk[Level][nb].JCu; j++) {
	for (int i = Local_SolnBlk[Level][nb].ICl; i <= Local_SolnBlk[Level][nb].ICu; i++) {
	  if (IP->Multigrid_IP.i_Physical_Time_Integration == TIME_STEPPING_IMPLICIT_EULER) {
 	    Local_SolnBlk[Level][nb].dudt[i][j][k_residual] /= ONE + ONE*IP->CFL_Number*Local_SolnBlk[Level][nb].dt[i][j]/(IP->Multigrid_IP.Physical_Time_CFL_Number*dt);
	  } else if (IP->Multigrid_IP.i_Physical_Time_Integration == TIME_STEPPING_IMPLICIT_SECOND_ORDER_BACKWARD) {
 	    Local_SolnBlk[Level][nb].dudt[i][j][k_residual] /= ONE + 1.5*IP->CFL_Number*Local_SolnBlk[Level][nb].dt[i][j]/(IP->Multigrid_IP.Physical_Time_CFL_Number*dt);
	  }
	}
      }
    }
  }

  // Reset the time-integration type.
  IP->i_Time_Integration = time_integration;

  // Solution residual successfully determined.
  return 0;

}

/**********************************************************************
 * FAS_Multigrid2D_Solver::Store_Previous_Solution --                 *
 *                                                                    *
 * Store the previous solution and restrict to the coarse levels.     *
 *                                                                    *
 **********************************************************************/
template <> int FAS_Multigrid2D_Solver<double,
				       AdvectDiffuse2D_Quad_Block,
				       AdvectDiffuse2D_Input_Parameters>::
Store_Previous_Solution(void) {

  int i_fine, j_fine;

  for (int nb = 0; nb < List_of_Local_Solution_Blocks[FINEST_LEVEL].Nblk; nb++) {
    if (List_of_Local_Solution_Blocks[FINEST_LEVEL].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      for (int j = Local_SolnBlk[FINEST_LEVEL][nb].JCl; j <= Local_SolnBlk[FINEST_LEVEL][nb].JCu; j++) {
	for (int i = Local_SolnBlk[FINEST_LEVEL][nb].ICl; i <= Local_SolnBlk[FINEST_LEVEL][nb].ICu; i++) {
 	  DTS_SolnBlk[FINEST_LEVEL][nb].Uo[i][j] = DTS_SolnBlk[FINEST_LEVEL][nb].Un[i][j];
 	  DTS_SolnBlk[FINEST_LEVEL][nb].Un[i][j] = Local_SolnBlk[FINEST_LEVEL][nb].U[i][j].u;
	}
      }
    }
  }

  for (int level = 0; level < IP->Multigrid_IP.Levels-1; level++) {
    for (int nb = 0; nb < List_of_Local_Solution_Blocks[level+1].Nblk; nb++) {
      if (List_of_Local_Solution_Blocks[level+1].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
 	for (int i_coarse = Local_SolnBlk[level+1][nb].ICl; i_coarse <= Local_SolnBlk[level+1][nb].ICu; i_coarse++) {
 	  for (int j_coarse = Local_SolnBlk[level+1][nb].JCl; j_coarse <= Local_SolnBlk[level+1][nb].JCu; j_coarse++) {
 	    // Determine the (i,j) index of the SW corner fine cell.
 	    i_fine = 2*(i_coarse-Local_SolnBlk[level][nb].Nghost)+Local_SolnBlk[level][nb].Nghost;
 	    j_fine = 2*(j_coarse-Local_SolnBlk[level][nb].Nghost)+Local_SolnBlk[level][nb].Nghost;
 	    // Determine the solution state of the coarse grid cell by
 	    // a area-weighted average of the associated fine grid cells.
 	    // Determine the solution state of the coarse grid cell by
 	    // a area-weighted average of the associated fine grid cells.
	    DTS_SolnBlk[level+1][nb].Un[i_coarse][j_coarse] = 
	      (DTS_SolnBlk[level][nb].Un[i_fine  ][j_fine  ]*Local_SolnBlk[level][nb].Grid.Cell[i_fine  ][j_fine  ].A +
	       DTS_SolnBlk[level][nb].Un[i_fine+1][j_fine  ]*Local_SolnBlk[level][nb].Grid.Cell[i_fine+1][j_fine  ].A +
	       DTS_SolnBlk[level][nb].Un[i_fine  ][j_fine+1]*Local_SolnBlk[level][nb].Grid.Cell[i_fine  ][j_fine+1].A +
	       DTS_SolnBlk[level][nb].Un[i_fine+1][j_fine+1]*Local_SolnBlk[level][nb].Grid.Cell[i_fine+1][j_fine+1].A)/
	      Local_SolnBlk[level+1][nb].Grid.Cell[i_coarse][j_coarse].A;
	    DTS_SolnBlk[level+1][nb].Uo[i_coarse][j_coarse] = 
	      (DTS_SolnBlk[level][nb].Uo[i_fine  ][j_fine  ]*Local_SolnBlk[level][nb].Grid.Cell[i_fine  ][j_fine  ].A +
	       DTS_SolnBlk[level][nb].Uo[i_fine+1][j_fine  ]*Local_SolnBlk[level][nb].Grid.Cell[i_fine+1][j_fine  ].A +
	       DTS_SolnBlk[level][nb].Uo[i_fine  ][j_fine+1]*Local_SolnBlk[level][nb].Grid.Cell[i_fine  ][j_fine+1].A +
	       DTS_SolnBlk[level][nb].Uo[i_fine+1][j_fine+1]*Local_SolnBlk[level][nb].Grid.Cell[i_fine+1][j_fine+1].A)/
	      Local_SolnBlk[level+1][nb].Grid.Cell[i_coarse][j_coarse].A;
 	  }
 	}
      }
    }
  }

  return 0;

}
