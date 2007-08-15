/**********************************************************************
 * Dusty2DQuadMultigrid.h: Dusty2D versions of some multigrid         *
 *                         functions that could not use templated     *
 *                         version.                                   *
 **********************************************************************/

#ifndef _DUSTY2D_QUAD_MULTIGRID_INCLUDED
#define _DUSTY2D_QUAD_MULTIGRID_INCLUDED

// Include 2D Dusty quadrilateral mesh solution header file.

#ifndef _DUSTY2D_QUAD_INCLUDED
#include "Dusty2DQuad.h"
#endif // _DUSTY2D_QUAD_INCLUDED

// Include 2D Dusty quadrilateral mesh solution header file.

#ifndef _FASMULTIGRID2D_INCLUDED
#include "../FASMultigrid2D/FASMultigrid2D.h"
#endif // _FASMULTIGRID2D_INCLUDED

/**********************************************************************
 * Dusty2D_Quad_Block -- Dusty2D Multigrid Subroutines.               *
 **********************************************************************/

/**********************************************************************
 * Routine: Output_Multigrid_Cells                                    *
 *                                                                    *
 * This routine writes out the cell-centred information for each      *
 * multigrid grid level (finest to coarest) for a 1D array of 2D      *
 * quadrilateral multi-block solution blocks.                         *
 *                                                                    *
 **********************************************************************/
// template <> int FAS_Multigrid2D_Solver<Dusty2D_cState, 
// 				       Dusty2D_Quad_Block,
// 				       Dusty2D_Input_Parameters>::
// Output_Multigrid_Cells(int &number_of_time_steps,
// 		       double &Time) {

//   int i, i_output_title;
//   char prefix[256], extension[256], output_file_name[256];
//   char *output_file_name_ptr;
//   ofstream output_file;

//   // Determine main prefix of output data file names.
//   i = 0;
//   while (1) {
//     if (IP->Output_File_Name[i] == ' ' ||
// 	IP->Output_File_Name[i] == '.') break;
//     prefix[i] = IP->Output_File_Name[i];
//     i = i + 1;
//     if (i > strlen(IP->Output_File_Name)) break;
//   }
//   prefix[i] = '\0';
//   strcat(prefix,"_multigrid_cells_");

//   // Output to a seperate file for each multigrid level.
//   for (int level = 0; level < IP->Multigrid_IP.Levels; level++) {

//     // Determine prefix of output data file names for current grid level.
//     strcpy(output_file_name,prefix);
//     sprintf(extension,"%.2d",level);
//     strcat(output_file_name,extension);
//     strcat(output_file_name,"_cpu");

//     // Determine output data file name for this processor.
//     sprintf(extension,"%.6d",List_of_Local_Solution_Blocks[level].ThisCPU);
//     strcat(extension,".dat");
//     strcat(output_file_name,extension);
//     output_file_name_ptr = output_file_name;

//     // Open the output data file.
//     output_file.open(output_file_name_ptr,ios::out);
//     if (output_file.bad()) return 1;

//     // Write the solution data for each solution block.
//     i_output_title = 1;
//     for (int nb = 0; nb < List_of_Local_Solution_Blocks[level].Nblk; nb++) {
//       if (List_of_Local_Solution_Blocks[level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
// 	Output_Cells_Tecplot(Local_SolnBlk[level][nb],
// 			     *IP,
// 			     number_of_time_steps,
// 			     Time,
// 			     List_of_Local_Solution_Blocks[level].Block[nb].gblknum,
// 			     i_output_title,
// 			     output_file);
// 	if (i_output_title) i_output_title = 0;
//       }
//     }

//     // Close the output data file.
//     output_file.close();

//   }

//   // Writing of multigrid grids data files complete.
//   return 0;

// }

/**********************************************************************
 * Routine: Restrict_Solution_Blocks                                  *
 *                                                                    *
 * Restrict solution from Level_Fine to Level_Coarse for all blocks   *
 * on local block list **Soln_ptr (overwrites any solution present on *
 * Level_Coarse)                                                      *
 *                                                                    *
 **********************************************************************/
// template <> void FAS_Multigrid2D_Solver<Dusty2D_cState, 
//                                         Dusty2D_Quad_Block,
//                                         Dusty2D_Input_Parameters>::
// Restrict_Solution_Blocks(const int Level_Fine) {

//   int i_fine, j_fine, Nghost, nghost;
//   int Level_Coarse = Level_Fine + 1;
//   double A, total_area;
//   Polygon Pc, Pf;
//   Dusty2D_cState Dusty2D_U_STDATM; Dusty2D_U_STDATM.Standard_Atmosphere();
//   Dusty2D_cState Dusty2D_U_VACUUM; Dusty2D_U_VACUUM.Vacuum();

//   // Determine if the restriction includes the ghost cells.
//   //if (IP->Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions) nghost = 0;
//   if (APPLY_BOUNDARY_CONDITIONS_ON_COARSE_MESH == true) nghost = 0;
//   else nghost = 1;

//   // Loop through each solution block.
//   for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level_Fine].Nblk; nb++) {
//     if (List_of_Local_Solution_Blocks[Level_Fine].Block[nb].used == ADAPTIVEBLOCK2D_USED) { 

//       // Get the number of ghost cells on the fine grid.
//       Nghost = Local_SolnBlk[Level_Fine][nb].Nghost;

//       // Loop through the coarse grid cells.
//       for (int i_coarse = Local_SolnBlk[Level_Coarse][nb].ICl-nghost; i_coarse <= Local_SolnBlk[Level_Coarse][nb].ICu+nghost; i_coarse++) {
// 	for (int j_coarse = Local_SolnBlk[Level_Coarse][nb].JCl-nghost; j_coarse <= Local_SolnBlk[Level_Coarse][nb].JCu+nghost; j_coarse++) {

// 	  // Determine the (i,j) index of the SW corner fine cell.
// 	  i_fine = 2*(i_coarse-Nghost)+Nghost;
// 	  j_fine = 2*(j_coarse-Nghost)+Nghost;

// 	  if (Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[i_coarse][j_coarse] != CELL_STATUS_ACTIVE) {
// 	    // The cell is inactive.  Set to standard atmosphere.
// 	    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse] = Dusty2D_U_STDATM;

// 	  } else if (!Local_SolnBlk[Level_Coarse][nb].Interface_Union_List.Ni) {
// 	    // No embedded boundaries exist.  Restrict the solution
// 	    // information from the fine cells to the coarse cell using an
// 	    // area-weighted average of the fine cell solution information.
// 	    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse] = (Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].A*
// 								     Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] +
// 								     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine].A*
// 								     Local_SolnBlk[Level_Fine][nb].U[i_fine+1][j_fine] +
// 								     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine+1].A*
// 								     Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine+1] +
// 								     Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine+1].A*
// 								     Local_SolnBlk[Level_Fine][nb].U[i_fine+1][j_fine+1])/
// 	                                                            (Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].A);

// 	  } else {

// 	    if (Local_SolnBlk[Level_Coarse][nb].Grid.cell_type[i_coarse][j_coarse] == CELL_TYPE_QUADRILATERAL &&
// 		Local_SolnBlk[Level_Fine][nb].Grid.cell_type[i_fine  ][j_fine  ] == CELL_TYPE_QUADRILATERAL &&
// 		Local_SolnBlk[Level_Fine][nb].Grid.cell_type[i_fine+1][j_fine  ] == CELL_TYPE_QUADRILATERAL &&
// 		Local_SolnBlk[Level_Fine][nb].Grid.cell_type[i_fine  ][j_fine+1] == CELL_TYPE_QUADRILATERAL &&
// 		Local_SolnBlk[Level_Fine][nb].Grid.cell_type[i_fine+1][j_fine+1] == CELL_TYPE_QUADRILATERAL &&
// 		Local_SolnBlk[Level_Fine][nb].Grid.cell_status[i_fine  ][j_fine  ] == CELL_STATUS_ACTIVE &&
// 		Local_SolnBlk[Level_Fine][nb].Grid.cell_status[i_fine+1][j_fine  ] == CELL_STATUS_ACTIVE &&
// 		Local_SolnBlk[Level_Fine][nb].Grid.cell_status[i_fine  ][j_fine+1] == CELL_STATUS_ACTIVE &&
// 		Local_SolnBlk[Level_Fine][nb].Grid.cell_status[i_fine+1][j_fine+1] == CELL_STATUS_ACTIVE) {
// 	      // None of the associated cells have been adjusted.
// 	      // Restrict the solution information using an area-weighted
// 	      // average of the cells from the stored adjusted fine mesh
// 	      // and the coarse mesh.
// 	      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse] = (Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].A*
// 								       Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] +
// 								       Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine].A*
// 								       Local_SolnBlk[Level_Fine][nb].U[i_fine+1][j_fine] +
// 								       Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine+1].A*
// 								       Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine+1] +
// 								       Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine+1].A*
// 								       Local_SolnBlk[Level_Fine][nb].U[i_fine+1][j_fine+1])/
// 	                                                              (Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].A);

// 	    } else {
// 	      // The coarse and/or fine cells have been adjusted.
// 	      // Restrict the solution information using an area-weight
// 	      // average based on the area of intersection of the coarse
// 	      // grid and fine grid cells.
// 	      total_area = ZERO;
// 	      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse] = Dusty2D_U_VACUUM;
// 	      // Create the coarse cell polygon.
// 	      Pc.convert(Local_SolnBlk[Level_Coarse][nb].Grid.nodeSW(i_coarse,j_coarse).X,
// 			 Local_SolnBlk[Level_Coarse][nb].Grid.nodeSE(i_coarse,j_coarse).X,
// 			 Local_SolnBlk[Level_Coarse][nb].Grid.nodeNE(i_coarse,j_coarse).X,
// 			 Local_SolnBlk[Level_Coarse][nb].Grid.nodeNW(i_coarse,j_coarse).X);
// 	      // Search the fine cells for intersecting cells.
// 	      for (int j = j_fine-2; j < j_fine+4; j++) {
// 		for (int i = i_fine-2; i < i_fine+4; i++) {
// 		  if (Local_SolnBlk[Level_Fine][nb].Grid.cell_status[i][j] == CELL_STATUS_ACTIVE) {
// 		    // Create the fine cell polygon.
// 		    Pf.convert(Local_SolnBlk[Level_Fine][nb].Grid.nodeSW(i,j).X,
// 			       Local_SolnBlk[Level_Fine][nb].Grid.nodeSE(i,j).X,
// 			       Local_SolnBlk[Level_Fine][nb].Grid.nodeNE(i,j).X,
// 			       Local_SolnBlk[Level_Fine][nb].Grid.nodeNW(i,j).X);
// 		    // Determine the area of intersection between the fine
// 		    // and coarse cell polygons.
// 		    A = Polygon_Intersection_Area(Pc,Pf);
// 		    Pf.deallocate();
// 		    // Add the area of intersection to the total area of
// 		    // intersection.
// 		    total_area += A;
// 		    // Restrict fine cell solution to the coarse cell
// 		    // based on the area of interesection.
// 		    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse] += Local_SolnBlk[Level_Fine][nb].U[i][j]*A;
// 		  }
// 		}
// 	      }
// 	      Pc.deallocate();
// 	      // Complete solution restriction by dividing the restricted
// 	      // solution by the total area of intersection.
// 	      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse] /= total_area;

// 	    }

// 	  }

// 	}
//       }

//     }
//   }

// }

/**********************************************************************
 * Routine: Restrict_Residuals                                        *
 *                                                                    *
 * Restrict residual dUdt from Level_Fine to Level_Coarse (stored in  *
 * P) for all blocks on local block list **Soln_ptr (overwrites any   *
 * value stored in P on Level_Coarse)                                 *
 *                                                                    *
 **********************************************************************/
// template <> void FAS_Multigrid2D_Solver<Dusty2D_cState, 
// 		  	                Dusty2D_Quad_Block,
// 			                Dusty2D_Input_Parameters>::
// Restrict_Residuals(const int Level_Fine) {

//   int i_fine, j_fine, Nghost, nghost;
//   int Level_Coarse = Level_Fine + 1;
//   double A, total_area;
//   Polygon Pc, Pf;
//   Dusty2D_cState Dusty2D_U_STDATM; Dusty2D_U_STDATM.Standard_Atmosphere();
//   Dusty2D_cState Dusty2D_U_VACUUM; Dusty2D_U_VACUUM.Vacuum();

//   // Determine if the restriction includes the ghost cells.
//   //if (IP->Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions) nghost = 0;
//   if (APPLY_BOUNDARY_CONDITIONS_ON_COARSE_MESH == true) nghost = 0;
//   else nghost = 1;

//   // Loop through each local solution block.
//   for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level_Fine].Nblk; nb++) {
//     if (List_of_Local_Solution_Blocks[Level_Fine].Block[nb].used == ADAPTIVEBLOCK2D_USED) { 

//       // Get the number of ghost cells.
//       Nghost = Local_SolnBlk[Level_Fine][nb].Nghost;

//       // Loop through the coarse cells.
//       for (int i_coarse = Local_SolnBlk[Level_Coarse][nb].ICl-nghost; i_coarse <= Local_SolnBlk[Level_Coarse][nb].ICu+nghost; i_coarse++) {
// 	for (int j_coarse = Local_SolnBlk[Level_Coarse][nb].JCl-nghost; j_coarse <= Local_SolnBlk[Level_Coarse][nb].JCu+nghost; j_coarse++) {
// 	  // Determine the (i,j) index of the corresponding SW corner
// 	  // fine cell.
// 	  i_fine = 2*(i_coarse-Nghost)+Nghost;
// 	  j_fine = 2*(j_coarse-Nghost)+Nghost;

// 	  if (Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[i_coarse][j_coarse] != CELL_STATUS_ACTIVE) {
// 	    // The cell is inactive.  Set to standard atmosphere.
// 	    MG_SolnBlk[Level_Coarse][nb].P[i_coarse][j_coarse] = Dusty2D_U_STDATM;

// 	  } else if (!Local_SolnBlk[Level_Coarse][nb].Interface_Union_List.Ni) {
// 	    // No embedded boundaries exist.  Restrict the solution
// 	    // information from the fine cells to the coarse cell using an
// 	    // area-weighted average of the fine cell solution information.
// 	    MG_SolnBlk[Level_Coarse][nb].P[i_coarse][j_coarse] = (Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].A*
// 								  Local_SolnBlk[Level_Fine][nb].dUdt[i_fine][j_fine][0] +
// 								  Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine].A*
// 								  Local_SolnBlk[Level_Fine][nb].dUdt[i_fine+1][j_fine][0] +
// 								  Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine+1].A*
// 								  Local_SolnBlk[Level_Fine][nb].dUdt[i_fine][j_fine+1][0] +
// 								  Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine+1].A*
// 								  Local_SolnBlk[Level_Fine][nb].dUdt[i_fine+1][j_fine+1][0])/
// 	                                                         (Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].A);

// 	  } else {

// 	    if (Local_SolnBlk[Level_Coarse][nb].Grid.cell_type[i_coarse][j_coarse] == CELL_TYPE_QUADRILATERAL &&
// 		Local_SolnBlk[Level_Fine][nb].Grid.cell_type[i_fine  ][j_fine  ] == CELL_TYPE_QUADRILATERAL &&
// 		Local_SolnBlk[Level_Fine][nb].Grid.cell_type[i_fine+1][j_fine  ] == CELL_TYPE_QUADRILATERAL &&
// 		Local_SolnBlk[Level_Fine][nb].Grid.cell_type[i_fine  ][j_fine+1] == CELL_TYPE_QUADRILATERAL &&
// 		Local_SolnBlk[Level_Fine][nb].Grid.cell_type[i_fine+1][j_fine+1] == CELL_TYPE_QUADRILATERAL &&
// 		Local_SolnBlk[Level_Fine][nb].Grid.cell_status[i_fine  ][j_fine  ] == CELL_STATUS_ACTIVE &&
// 		Local_SolnBlk[Level_Fine][nb].Grid.cell_status[i_fine+1][j_fine  ] == CELL_STATUS_ACTIVE &&
// 		Local_SolnBlk[Level_Fine][nb].Grid.cell_status[i_fine  ][j_fine+1] == CELL_STATUS_ACTIVE &&
// 		Local_SolnBlk[Level_Fine][nb].Grid.cell_status[i_fine+1][j_fine+1] == CELL_STATUS_ACTIVE) {
// 	      // None of the associated cells have been adjusted.
// 	      // Restrict the solution information using an area-weighted
// 	      // average of the cells from the stored adjusted fine mesh
// 	      // and the coarse mesh.
// 	      MG_SolnBlk[Level_Coarse][nb].P[i_coarse][j_coarse] = (Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].A*
// 								    Local_SolnBlk[Level_Fine][nb].dUdt[i_fine][j_fine][0] +
// 								    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine].A*
// 								    Local_SolnBlk[Level_Fine][nb].dUdt[i_fine+1][j_fine][0] +
// 								    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine+1].A*
// 								    Local_SolnBlk[Level_Fine][nb].dUdt[i_fine][j_fine+1][0] +
// 								    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine+1][j_fine+1].A*
// 								    Local_SolnBlk[Level_Fine][nb].dUdt[i_fine+1][j_fine+1][0])/
// 	                                                           (Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].A);

// 	    } else {
// 	      // The coarse and/or fine cells have been adjusted.
// 	      // Restrict the solution information using an area-weight
// 	      // average based on the area of intersection of the coarse
// 	      // grid and fine grid cells.
// 	      total_area = ZERO;
// 	      MG_SolnBlk[Level_Coarse][nb].P[i_coarse][j_coarse] = Dusty2D_U_VACUUM;
// 	      // Create the coarse cell polygon.
// 	      Pc.convert(Local_SolnBlk[Level_Coarse][nb].Grid.nodeSW(i_coarse,j_coarse).X,
// 			 Local_SolnBlk[Level_Coarse][nb].Grid.nodeSE(i_coarse,j_coarse).X,
// 			 Local_SolnBlk[Level_Coarse][nb].Grid.nodeNE(i_coarse,j_coarse).X,
// 			 Local_SolnBlk[Level_Coarse][nb].Grid.nodeNW(i_coarse,j_coarse).X);
// 	      // Search the fine cells for intersecting cells.
// 	      for (int j = j_fine-2; j < j_fine+4; j++) {
// 		for (int i = i_fine-2; i < i_fine+4; i++) {
//   		  if (Local_SolnBlk[Level_Fine][nb].Grid.cell_status[i][j] == CELL_STATUS_ACTIVE) {
// 		    // Create the fine cell polygon.
// 		    Pf.convert(Local_SolnBlk[Level_Fine][nb].Grid.nodeSW(i,j).X,
// 			       Local_SolnBlk[Level_Fine][nb].Grid.nodeSE(i,j).X,
// 			       Local_SolnBlk[Level_Fine][nb].Grid.nodeNE(i,j).X,
// 			       Local_SolnBlk[Level_Fine][nb].Grid.nodeNW(i,j).X);
// 		    // Determine the area of intersection between the fine
// 		    // and coarse cell polygons.
// 		    A = Polygon_Intersection_Area(Pc,Pf);
// 		    Pf.deallocate();
// 		    // Add the area of intersection to the total area of
// 		    // intersection.
// 		    total_area += A;
// 		    // Restrict fine cell solution to the coarse cell
// 		    // based on the area of interesection.
// 		    MG_SolnBlk[Level_Coarse][nb].P[i_coarse][j_coarse] += Local_SolnBlk[Level_Fine][nb].dUdt[i][j][0]*A;
// 		  }
// 		}
// 	      }
// 	      Pc.deallocate();
// 	      // Complete solution restriction by dividing the restricted
// 	      // solution by the total area of intersection.
// 	      MG_SolnBlk[Level_Coarse][nb].P[i_coarse][j_coarse] /= total_area;

// 	    }

// 	  }

// 	}
//       }

//     }
//   }

// }

/**********************************************************************
 * Routine: Prolong_Solution_Blocks (for Multigrid)                   *
 *                                                                    *
 * Prolong solution states stored in U[][] from Level_Coarse to       *
 * Level_Fine for all blocks on local block list **Soln_ptr.          *
 *                                                                    *
 **********************************************************************/
// template <> int FAS_Multigrid2D_Solver<Dusty2D_cState, 
// 			               Dusty2D_Quad_Block,
// 			               Dusty2D_Input_Parameters>::
// Prolong_Solution_Blocks(const int Level_Coarse) {

//   int error_flag, i_fine, j_fine, Nghost;
//   int Level_Fine = Level_Coarse - 1;
//   Dusty2D_cState Ufine;
//   bool injection_for_this_coarse_cell;
//   injection_for_this_coarse_cell = false;
//   int i_coarse_min, j_coarse_min, coarse_cell_found;
//   double distance;
//   Dusty2D_cState Dusty2D_U_STDATM;
//   Dusty2D_pState Dusty2D_W_STDATM;

//   // Loop through each local solution block.
//   for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level_Coarse].Nblk; nb++) {
//     if (List_of_Local_Solution_Blocks[Level_Coarse].Block[nb].used == ADAPTIVEBLOCK2D_USED) { 

//       // Get the number of ghost cells on the coarse grid.
//       Nghost = Local_SolnBlk[Level_Coarse][nb].Nghost;

//       // Loop through the coarse grid cells.
//       for (int i_coarse = Local_SolnBlk[Level_Coarse][nb].ICl; i_coarse <= Local_SolnBlk[Level_Coarse][nb].ICu; i_coarse++) {
// 	for (int j_coarse = Local_SolnBlk[Level_Coarse][nb].JCl; j_coarse <= Local_SolnBlk[Level_Coarse][nb].JCu; j_coarse++) {

//  	  if (Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[i_coarse-1][j_coarse-1] != CELL_STATUS_ACTIVE ||
//  	      Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[i_coarse  ][j_coarse-1] != CELL_STATUS_ACTIVE ||
//  	      Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[i_coarse+1][j_coarse-1] != CELL_STATUS_ACTIVE ||
//  	      Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[i_coarse-1][j_coarse  ] != CELL_STATUS_ACTIVE ||
//  	      Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[i_coarse  ][j_coarse  ] != CELL_STATUS_ACTIVE ||
//  	      Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[i_coarse+1][j_coarse  ] != CELL_STATUS_ACTIVE ||
//  	      Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[i_coarse-1][j_coarse+1] != CELL_STATUS_ACTIVE ||
//  	      Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[i_coarse  ][j_coarse+1] != CELL_STATUS_ACTIVE ||
//  	      Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[i_coarse+1][j_coarse+1] != CELL_STATUS_ACTIVE) {
//  	    // The fine cell is near an embedded boundary.  Prolong by 
//  	    // injecting the solution information from the nearest
//  	    // active cell.
//  	    i_fine = 2*(i_coarse-Nghost)+Nghost;
//  	    j_fine = 2*(j_coarse-Nghost)+Nghost;
//  	    for (int jj = j_fine; jj <= j_fine+1; jj++) {
//  	      for (int ii = i_fine; ii <= i_fine+1; ii++) {
//  		if (Local_SolnBlk[Level_Fine][nb].Grid.cell_status[ii][jj] != CELL_STATUS_ACTIVE) {
// 		  // The fine grid cell is inactive.  Set to standard
// 		  // atmosphere.
//  		  Local_SolnBlk[Level_Fine][nb].U[ii][jj] = Dusty2D_U_STDATM;
//  		} else {
// 		  // The fine grid cell is inactive.  Inject the
// 		  // solution information from the nearest coarse grid
// 		  // cell.
//  		  distance = MILLION;
//  		  coarse_cell_found = OFF;
//  		  for (int jc = j_coarse-1; jc <= j_coarse+1; jc++) {
//  		    for (int ic = i_coarse-1; ic <= i_coarse+1; ic++) {
//  		      if (Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[ic][jc] == CELL_STATUS_ACTIVE &&
//  			  abs(Local_SolnBlk[Level_Coarse][nb].Grid.Cell[ic][jc].Xc -
//  			      Local_SolnBlk[Level_Fine][nb].Grid.Cell[ii][jj].Xc) < distance) {
//  			i_coarse_min = ic; j_coarse_min = jc;
//  			distance = abs(Local_SolnBlk[Level_Coarse][nb].Grid.Cell[ic][jc].Xc -
//  				       Local_SolnBlk[Level_Fine][nb].Grid.Cell[ii][jj].Xc);
//  			coarse_cell_found = ON;
//  		      }
//  		    }
//  		  }
//  		  if (coarse_cell_found) {
//  		    Local_SolnBlk[Level_Fine][nb].U[ii][jj] = Local_SolnBlk[Level_Coarse][nb].U[i_coarse_min][j_coarse_min];
//  		  } else {
// 		    return 27081;
//  		  }
//  		}
//  	      }
//  	    }

//  	  } else {

//  	    // The fine cell is not near an embedded boundary.  Prolong
// 	    // by either by injection or a bilinear interpolation.

// 	    // Set injection flag if prolongation by injection is desired.
// 	    if (IP->Multigrid_IP.Prolong_Using_Injection) injection_for_this_coarse_cell = true;

// 	    // If the bilinear interpolation is used as the prolongation
// 	    // operator and a valid set of interpolation points can not
// 	    // be found for one or more of the fine grid cells then
// 	    // switch the prolongation operator to injection for all
// 	    // four fine grid cells.  This is facilitated by the use of
// 	    // a goto statement where 'restart' is used as the goto
// 	    // label.
// 	  restart: ;

// 	    /**********************************************************
// 	     * Determine the cell-centred solution state for the      *
// 	     * south-west fine grid cell.                             *
// 	     **********************************************************/

// 	    // Determine the (i,j) index of the SW fine grid cell.
// 	    i_fine = 2*(i_coarse-Nghost)+Nghost;
// 	    j_fine = 2*(j_coarse-Nghost)+Nghost;

// 	    // Attempt to use bilinear interpolation if the injection
// 	    // flag has not been set true.
// 	    if (!injection_for_this_coarse_cell) {

// 	      // Attempt to use the SW coarse grid solution states as
// 	      // the interpolants ([i-1,i],[j-1,j]).
// 	      error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse-1],
// 						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
// 						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
// 						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
// 						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
// 						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
// 						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 						  Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
// 						  Ufine);
// 	      // If fine grid cell-centre falls outside of the SW
// 	      // coarse grid interpolants then try the NW coarse grid
// 	      // solution states as the interpolants ([i-1,i],[j,j+1]).
// 	      if (error_flag) {
// 		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse+1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
// 						    Ufine);
// 	      }
// 	      // If fine grid cell-centre falls outside of the NW
// 	      // coarse grid interpolants then try the SE coarse grid
// 	      // solution states as the interpolants ([i,i+1],[j-1,j]).
// 	      if (error_flag) {
// 		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse-1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
// 						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
// 						    Ufine);
// 	      }
// 	      // If fine grid cell-centre falls outside of the SE
// 	      // coarse grid interpolants then try the NE coarse grid
// 	      // solution states as the interpolants ([i,i+1],[j,j+1]).
// 	      if (error_flag != 0) {
// 		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse+1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
// 						    Ufine);
// 	      }
// 	      // If none of the four sets of possible interpolants 
// 	      // provide a suitable set of interpolants then restart
// 	      // the prolongation associated with the current coarse
// 	      // grid cell using injection.
// 	      if (error_flag) {
// 		injection_for_this_coarse_cell = true;
// 		goto restart;
// 	      }

// 	    } else {

// 	      // Perform the prolongation using injection.
// 	      Ufine = Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse];

// 	    }

// 	    // Assign the SW fine grid cell solution state.
// 	    Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] = Ufine;

// 	    /**********************************************************
// 	     * Determine the cell-centred solution state for the      *
// 	     * south-east fine grid cell.                             *
// 	     **********************************************************/

// 	    // Determine the (i,j) index of the SE fine grid cell.
// 	    i_fine = 2*(i_coarse-Nghost)+Nghost+1;
// 	    j_fine = 2*(j_coarse-Nghost)+Nghost;

// 	    // Attempt to use bilinear interpolation if the injection
// 	    // flag has not been set true.
// 	    if (!injection_for_this_coarse_cell) {

// 	      // Attempt to use the SE coarse grid solution states as
// 	      // the interpolants ([i,i+1],[j-1,j]).
// 	      error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
// 						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
// 						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
// 						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
// 						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse-1],
// 						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
// 						  Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
// 						  Ufine);
// 	      // If fine grid cell-centre falls outside of the SE coarse
// 	      // grid interpolants then try the NE coarse grid solution
// 	      // states as the interpolants ([i,i+1],[j,j+1]).
// 	      if (error_flag) {
// 		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse+1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
// 						    Ufine);
// 	      }
// 	      // If fine grid cell-centre falls outside of the NE coarse
// 	      // grid interpolants then try the SW coarse grid solution
// 	      // states as the interpolants ([i-1,i],[j-1,j]).
// 	      if (error_flag) {
// 		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse-1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
// 						    Ufine);
// 	      }
// 	      // If fine grid cell-centre falls outside of the SW coarse
// 	      // grid interpolants then try the N coarse grid solution
// 	      // states as the interpolants ([i-1,i],[j,j+1]).
// 	      if (error_flag) {
// 		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse+1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
// 						    Ufine);
// 	      }
// 	      // If none of the four sets of possible interpolants 
// 	      // provide a suitable set of interpolants then restart
// 	      // the prolongation associated with the current coarse
// 	      // grid cell using injection.
// 	      if (error_flag != 0) {
// 		injection_for_this_coarse_cell = true;
// 		goto restart;
// 	      }

// 	    } else {

// 	      // Perform the prolongation using injection.
// 	      Ufine = Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse];

// 	    }

// 	    // Assign the SE fine grid cell solution state.
// 	    Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] = Ufine;

// 	    /************************************************************
// 	     * Determine the cell-centred solution state for the north- *
// 	     * east fine grid cell.                                     *
// 	     ************************************************************/

// 	    // Determine the (i,j) index of the NE fine grid cell.
// 	    i_fine = 2*(i_coarse-Nghost)+Nghost+1;
// 	    j_fine = 2*(j_coarse-Nghost)+Nghost+1;

// 	    // Attempt to use bilinear interpolation if the injection flag
// 	    // has not been set true.
// 	    if (!injection_for_this_coarse_cell) {

// 	      // Attempt to use the NE coarse grid solution states as the
// 	      // interpolants ([i,i+1],[j,j+1]).
// 	      error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
// 						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
// 						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
// 						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse+1],
// 						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
// 						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
// 						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 						  Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
// 						  Ufine);
// 	      // If fine grid cell-centre falls outside of the NE coarse
// 	      // grid interpolants then try the SE coarse grid solution
// 	      // states as the interpolants ([i,i+1],[j-1,j]).
// 	      if (error_flag) {
// 		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse-1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
// 						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
// 						    Ufine);
// 	      }
// 	      // If fine grid cell-centre falls outside of the SE coarse
// 	      // grid interpolants then try the NW coarse grid solution
// 	      // states as the interpolants ([i-1,i],[j,j+1]).
// 	      if (error_flag) {
// 		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse+1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
// 						    Ufine);
// 	      }
// 	      // If fine grid cell-centre falls outside of the NW coarse
// 	      // grid interpolants then try the SW coarse grid solution
// 	      // states as the interpolants ([i-1,i],[j-1,j]).
// 	      if (error_flag) {
// 		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse-1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
// 						    Ufine);
// 	      }
// 	      // If none of the four sets of possible interpolants 
// 	      // provide a suitable set of interpolants then restart
// 	      // the prolongation associated with the current coarse
// 	      // grid cell using injection.
// 	      if (error_flag != 0) {
// 		injection_for_this_coarse_cell = true;
// 		goto restart;
// 	      }

// 	    } else {

// 	      // Perform the prolongation using injection.
// 	      Ufine = Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse];

// 	    }

// 	    // Assign the NE fine grid cell solution state.
// 	    Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] = Ufine;

// 	    /************************************************************
// 	     * Determine the cell-centred solution state for the north- *
// 	     * west fine grid cell.                                     *
// 	     ************************************************************/

// 	    // Determine the (i,j) index of the NW fine grid cell.
// 	    i_fine = 2*(i_coarse-Nghost)+Nghost;
// 	    j_fine = 2*(j_coarse-Nghost)+Nghost+1;

// 	    // Attempt to use bilinear interpolation if the injection flag
// 	    // has not been set true.
// 	    if (!injection_for_this_coarse_cell) {

// 	      // Attempt to use the NW coarse grid solution states as the
// 	      // interpolants ([i-1,i],[j,j+1]).
// 	      error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
// 						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse+1],
// 						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
// 						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
// 						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 						  Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
// 						  Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
// 						  Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
// 						  Ufine);
// 	      // If fine grid cell-centre falls outside of the NW coarse
// 	      // grid interpolants then try the NE coarse grid solution
// 	      // states as the interpolants ([i,i+1],[j,j+1]).
// 	      if (error_flag) {
// 		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse+1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
// 						    Ufine);
// 	      }
// 	      // If fine grid cell-centre falls outside of the NE coarse
// 	      // grid interpolants then try the SW coarse grid solution
// 	      // states as the interpolants ([i-1,i],[j-1,j]).
// 	      if (error_flag) {
// 		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse-1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
// 						    Ufine);
// 	      }
// 	      // If fine grid cell-centre falls outside of the SW coarse
// 	      // grid interpolants then try the SE coarse grid solution
// 	      // states as the interpolants ([i,i+1],[j-1,j]).
// 	      if (error_flag) {
// 		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
// 						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse-1],
// 						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
// 						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
// 						    Ufine);
// 	      }
// 	      // If none of the four sets of possible interpolants 
// 	      // provide a suitable set of interpolants then restart
// 	      // the prolongation associated with the current coarse
// 	      // grid cell using injection.
// 	      if (error_flag != 0) {
// 		injection_for_this_coarse_cell = true;
// 		goto restart;
// 	      }

// 	    } else {

// 	      // Perform the prolongation using injection.
// 	      Ufine = Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse];

// 	    }

// 	    // Assign the NW fine grid cell solution state.
// 	    Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] = Ufine;

// 	    // Reset the injection flag to false.
// 	    injection_for_this_coarse_cell = false;

// 	  }

// 	}
//       }

//     }
//   }

//   // Prolongation of the solution block was successful.
//   return 0;

// }

/**********************************************************************
 * Routine: Prolong_and_Update_Solution_Blocks                        *
 *                                                                    *
 * Prolong solution changes stored in U[][].u from Level_Coarse to    *
 * Level_Fine for all blocks on local block list **Soln_ptr and add   *
 * the changes to the solution U[][].u on Level_Fine blocks (modifies *
 * any value present in U[][] on Level_Fine).                         *
 *                                                                    *
 **********************************************************************/
// template <> int FAS_Multigrid2D_Solver<Dusty2D_cState, 
//                                        Dusty2D_Quad_Block,
//         		               Dusty2D_Input_Parameters>::
// Prolong_and_Update_Solution_Blocks(const int Level_Coarse) {

//   int error_flag, i_fine, j_fine, Nghost, residual_reduction_number;
//   int Level_Fine = Level_Coarse - 1;
//   Dusty2D_cState Uchange, Unew;
//   bool injection_for_this_coarse_cell;
//   injection_for_this_coarse_cell = false;
//   int i_coarse_min, j_coarse_min, coarse_cell_found;
//   double distance;
//   Dusty2D_cState Dusty2D_U_STDATM;

//   // Loop through each local solution block.
//   for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level_Coarse].Nblk; nb++) {
//     if (List_of_Local_Solution_Blocks[Level_Coarse].Block[nb].used == ADAPTIVEBLOCK2D_USED) { 

//       // Get the number of ghost cells on the coarse grid.
//       Nghost = Local_SolnBlk[Level_Coarse][nb].Nghost;

//       // Enforce Neumann boundary conditions as appropriate on the 
//       // coarse grid if bilinear interpolation is used as the 
//       // prolongation operator.
//       if (!IP->Multigrid_IP.Prolong_Using_Injection) {

// 	// West and east face Neumann boundary conditions.
// 	for (int j_coarse = Local_SolnBlk[Level_Coarse][nb].JCl; j_coarse <= Local_SolnBlk[Level_Coarse][nb].JCu; j_coarse++) {
// 	  // East boundary.
// 	  if (Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeE[j_coarse] == BC_CONSTANT_EXTRAPOLATION) {
// 	    Local_SolnBlk[Level_Coarse][nb].U[Local_SolnBlk[Level_Coarse][nb].ICu+1][j_coarse] =
// 	      Local_SolnBlk[Level_Coarse][nb].U[Local_SolnBlk[Level_Coarse][nb].ICu][j_coarse];
// 	  }
// 	  // West boundary.
// 	  if (Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeW[j_coarse] == BC_CONSTANT_EXTRAPOLATION) {
// 	    Local_SolnBlk[Level_Coarse][nb].U[Local_SolnBlk[Level_Coarse][nb].ICl-1][j_coarse] =
// 	      Local_SolnBlk[Level_Coarse][nb].U[Local_SolnBlk[Level_Coarse][nb].ICl][j_coarse];
// 	  }
// 	}

// 	// North and south face Neumann boundary coditions.
// 	for (int i_coarse = Local_SolnBlk[Level_Coarse][nb].ICl; i_coarse <= Local_SolnBlk[Level_Coarse][nb].ICu; i_coarse++) {
// 	  // South boundary.
// 	  if (Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeS[i_coarse] == BC_CONSTANT_EXTRAPOLATION) {
// 	    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][Local_SolnBlk[Level_Coarse][nb].JCl-1] =
// 	      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][Local_SolnBlk[Level_Coarse][nb].JCl];
// 	  }
// 	  // North boundary.
// 	  if (Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeN[i_coarse] == BC_CONSTANT_EXTRAPOLATION) {
// 	    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][Local_SolnBlk[Level_Coarse][nb].JCu+1] =
// 	      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][Local_SolnBlk[Level_Coarse][nb].JCu];
// 	  }
// 	}
//       }

//       // Loop through the coarse grid cells.
//       for (int i_coarse = Local_SolnBlk[Level_Coarse][nb].ICl; i_coarse <= Local_SolnBlk[Level_Coarse][nb].ICu; i_coarse++) {
// 	for (int j_coarse = Local_SolnBlk[Level_Coarse][nb].JCl; j_coarse <= Local_SolnBlk[Level_Coarse][nb].JCu; j_coarse++) {

// 	  if (Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[i_coarse-1][j_coarse-1] != CELL_STATUS_ACTIVE ||
// 	      Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[i_coarse  ][j_coarse-1] != CELL_STATUS_ACTIVE ||
// 	      Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[i_coarse+1][j_coarse-1] != CELL_STATUS_ACTIVE ||
// 	      Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[i_coarse-1][j_coarse  ] != CELL_STATUS_ACTIVE ||
// 	      Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[i_coarse  ][j_coarse  ] != CELL_STATUS_ACTIVE ||
// 	      Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[i_coarse+1][j_coarse  ] != CELL_STATUS_ACTIVE ||
// 	      Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[i_coarse-1][j_coarse+1] != CELL_STATUS_ACTIVE ||
// 	      Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[i_coarse  ][j_coarse+1] != CELL_STATUS_ACTIVE ||
// 	      Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[i_coarse+1][j_coarse+1] != CELL_STATUS_ACTIVE) {
// 	    // The fine cell is near an embedded boundary.  Prolong by 
// 	    // injecting the solution information from the nearest
// 	    // active cell.
// 	    i_fine = 2*(i_coarse-Nghost)+Nghost;
// 	    j_fine = 2*(j_coarse-Nghost)+Nghost;
// 	    for (int jj = j_fine; jj <= j_fine+1; jj++) {
// 	      for (int ii = i_fine; ii <= i_fine+1; ii++) {
// 		if (Local_SolnBlk[Level_Fine][nb].Grid.cell_status[ii][jj] == CELL_STATUS_INACTIVE) {
// 		  // The fine grid cell is inactive.  Set to standard
// 		  // atmosphere.
// 		  Uchange = Dusty2D_U_STDATM;

// 		} else {
// 		  // The fine grid cell is inactive.  Inject the
// 		  // solution information from the nearest coarse grid
// 		  // cell.
// 		  distance = MILLION;
// 		  coarse_cell_found = OFF;
// 		  for (int jc = j_coarse-1; jc <= j_coarse+1; jc++) {
// 		    for (int ic = i_coarse-1; ic <= i_coarse+1; ic++) {
// 		      if (Local_SolnBlk[Level_Coarse][nb].Grid.cell_status[ic][jc] == CELL_STATUS_ACTIVE &&
// 			  abs(Local_SolnBlk[Level_Coarse][nb].Grid.Cell[ic][jc].Xc -
// 			      Local_SolnBlk[Level_Fine][nb].Grid.Cell[ii][jj].Xc) < distance) {
// 			i_coarse_min = ic; j_coarse_min = jc;
// 			distance = abs(Local_SolnBlk[Level_Coarse][nb].Grid.Cell[ic][jc].Xc -
// 				       Local_SolnBlk[Level_Fine][nb].Grid.Cell[ii][jj].Xc);
// 			coarse_cell_found = ON;
// 		      }
// 		    }
// 		  }
// 		  if (coarse_cell_found) {
// 		    Uchange = Local_SolnBlk[Level_Coarse][nb].U[i_coarse_min][j_coarse_min];
// 		  } else {
// 		    return 27082;
// 		  }
// 		}

// 		// Update the conservative solution state in the fine
// 		// grid cell.  If required, use residual reduction if
// 		// the update creates a unphysical state (negative 
// 		// density or pressure).  Residual reduction is
// 		// performed by halving the residual update.
// 		if (UPDATE_STABILITY_SWITCH) {
// 		  Unew = Local_SolnBlk[Level_Fine][nb].U[ii][jj] + Uchange;
// 		  if (Unew[1] < 0 || Unew.p() < 0) {
// 		    residual_reduction_number = 1;
// 		    while (residual_reduction_number <= NUMBER_OF_UPDATE_REDUCTIONS) {
// 		      Uchange *= HALF;
// 		      Unew = Local_SolnBlk[Level_Fine][nb].U[ii][jj] + Uchange;
// 		      if (Unew[1] > ZERO && Unew.p() > ZERO) break;
// 		      residual_reduction_number++;
// 		    }
// 		    if (residual_reduction_number > NUMBER_OF_UPDATE_REDUCTIONS) Uchange *= ZERO;
// 		  }
// 		}
// 		Local_SolnBlk[Level_Fine][nb].U[ii][jj] += Uchange;
// 	      }
// 	    }

// 	  } else {

// 	    // If the coarse grid cell is next to a Dirichlet (fixed) 
// 	    // boundary condition then prolong using injection instead of 
// 	    // interpolation for all four of the fine cells.  Interpolation
// 	    // for these cells can lead to a skewed (unphysical) result.
// 	    if ((i_coarse == Local_SolnBlk[Level_Coarse][nb].ICl &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeW[j_coarse] == BC_FIXED) ||
// 		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICu &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeE[j_coarse] == BC_FIXED) ||
// 		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCl &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeS[i_coarse] == BC_FIXED) ||
// 		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCu &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeN[i_coarse] == BC_FIXED) ||
// 		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICl &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeW[j_coarse] == BC_INFLOW_SUBSONIC) ||
// 		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICu &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeE[j_coarse] == BC_INFLOW_SUBSONIC) ||
// 		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCl &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeS[i_coarse] == BC_INFLOW_SUBSONIC) ||
// 		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCu &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeN[i_coarse] == BC_INFLOW_SUBSONIC) ||
// 		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICl &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeW[j_coarse] == BC_OUTFLOW_SUBSONIC) ||
// 		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICu &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeE[j_coarse] == BC_OUTFLOW_SUBSONIC) ||
// 		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCl &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeS[i_coarse] == BC_OUTFLOW_SUBSONIC) ||
// 		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCu &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeN[i_coarse] == BC_OUTFLOW_SUBSONIC) ||
// 		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICl &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeW[j_coarse] == BC_FIXED_PRESSURE) ||
// 		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICu &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeE[j_coarse] == BC_FIXED_PRESSURE) ||
// 		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCl &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeS[i_coarse] == BC_FIXED_PRESSURE) ||
// 		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCu &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeN[i_coarse] == BC_FIXED_PRESSURE) ||
// 		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICl &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeW[j_coarse] == BC_FIXED_TEMP_WALL) ||
// 		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICu &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeE[j_coarse] == BC_FIXED_TEMP_WALL) ||
// 		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCl &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeS[i_coarse] == BC_FIXED_TEMP_WALL) ||
// 		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCu &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeN[i_coarse] == BC_FIXED_TEMP_WALL) ||
// 		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICl &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeW[j_coarse] == BC_ADIABATIC_WALL) ||
// 		(i_coarse == Local_SolnBlk[Level_Coarse][nb].ICu &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeE[j_coarse] == BC_ADIABATIC_WALL) ||
// 		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCl &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeS[i_coarse] == BC_ADIABATIC_WALL) ||
// 		(j_coarse == Local_SolnBlk[Level_Coarse][nb].JCu &&
// 		 Local_SolnBlk[Level_Coarse][nb].Grid.BCtypeN[i_coarse] == BC_ADIABATIC_WALL)) {

// 	      // Determine the (i,j) index of the SW fine grid cell.
// 	      i_fine = 2*(i_coarse-Nghost)+Nghost;
// 	      j_fine = 2*(j_coarse-Nghost)+Nghost;

// 	      // Perform the prolongation using injection.
// 	      Uchange = Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse];
// 	      Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] += Uchange;
// 	      Local_SolnBlk[Level_Fine][nb].U[i_fine+1][j_fine] += Uchange;
// 	      Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine+1] += Uchange;
// 	      Local_SolnBlk[Level_Fine][nb].U[i_fine+1][j_fine+1] += Uchange;

// 	    } else {

//  	      // Set injection flag if prolongation by injection is desired.
//  	      if (IP->Multigrid_IP.Prolong_Using_Injection) injection_for_this_coarse_cell = true;

//  	      // If the bilinear interpolation is used as the prolongation
//  	      // operator and a valid set of interpolation points can not be
//  	      // found for one or more of the fine grid cells then switch
//  	      // the prolongation operator to injection for all four fine
//  	      // grid cells.  This is facilitated by the use of a goto
//  	      // statement where 'restart' is used as the goto label.
//  	    restart: ;

//  	      /**********************************************************
//  	       * Update the cell-centred solution state for the south-  *
//  	       * west fine grid cell.                                   *
//  	       **********************************************************/

//  	      // Determine the (i,j) index of the SW fine grid cell.
//  	      i_fine = 2*(i_coarse-Nghost)+Nghost;
//  	      j_fine = 2*(j_coarse-Nghost)+Nghost;

// 	      // Attempt to use bilinear interpolation if the injection
//  	      // flag has not been set true.
//  	      if (!injection_for_this_coarse_cell) {

//  		// Attempt to use the SW coarse grid solution states as
//  		// the interpolants ([i-1,i],[j-1,j]).
//  		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse-1],
//  						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
//  						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
//  						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
//  						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
//  						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
//  						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
//  						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
//  						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
//  						    Uchange);
//  		// If fine grid cell-centre falls outside of the SW coarse
//  		// grid interpolants then try the NW coarse grid solution
//  		// states as the interpolants ([i-1,i],[j,j+1]).
//  		if (error_flag) {
//  		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse+1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
//  						      Uchange);
//  		}
//  		// If fine grid cell-centre falls outside of the NW coarse
//  		// grid interpolants then try the SE coarse grid solution
//  		// states as the interpolants ([i,i+1],[j-1,j]).
//  		if (error_flag) {
//  		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse-1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
//  						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
//  						      Uchange);
//  		}
//  		// If fine grid cell-centre falls outside of the SE coarse
//  		// grid interpolants then try the NE coarse grid solution
//  		// states as the interpolants ([i,i+1],[j,j+1]).
//  		if (error_flag) {
//  		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse+1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
//  						      Uchange);
//  		}
//  		// If none of the four sets of possible interpolants 
//  		// provide a suitable set of interpolants then restart
//  		// the prolongation associated with the current coarse
//  		// grid cell using injection.
//  		if (error_flag) {
//  		  injection_for_this_coarse_cell = true;
//  		  goto restart;
//  		}

//  	      } else {

//  		// Perform the prolongation using injection.
//  		Uchange = Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse];

//  	      }

//  	      // Update the conservative solution state in the SW fine
//  	      // grid cell.  If required, use residual reduction if the
//  	      // update creates a unphysical state (negative density or 
//  	      // pressure).  Residual reduction is performed by halving
//  	      // the residual update.
//  	      if (UPDATE_STABILITY_SWITCH) {
//  		Unew = Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] + Uchange;
//  		if (Unew[1] < ZERO || Unew.p() < ZERO) {
//  		  residual_reduction_number = 1;
//  		  while (residual_reduction_number <= NUMBER_OF_UPDATE_REDUCTIONS) {
//  		    Uchange *= HALF;
//  		    Unew = Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] + Uchange;
//  		    if (Unew[1] > ZERO && Unew.p() > ZERO) break;
//  		    residual_reduction_number++;
//  		  }
//  		  if (residual_reduction_number > NUMBER_OF_UPDATE_REDUCTIONS) Uchange *= ZERO;
//  		}
//  	      }
//  	      Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] += Uchange;

//  	      /**********************************************************
//  	       * Update the cell-centred solution state for the south-  *
//  	       * east fine grid cell.                                   *
//  	       **********************************************************/

//  	      // Determine the (i,j) index of the SE fine grid cell.
//  	      i_fine = 2*(i_coarse-Nghost)+Nghost+1;
//  	      j_fine = 2*(j_coarse-Nghost)+Nghost;

//  	      // Attempt to use bilinear interpolation if the injection
//  	      // flag has not been set true.
//  	      if (!injection_for_this_coarse_cell) {

//  		// Attempt to use the SE coarse grid solution states as
//  		// the interpolants ([i,i+1],[j-1,j]).
//  		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
//  						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
//  						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
//  						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
//  						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
//  						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
//  						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse-1],
//  						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
//  						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
//  						    Uchange);
//  		// If fine grid cell-centre falls outside of the SE coarse
//  		// grid interpolants then try the NE coarse grid solution
//  		// states as the interpolants ([i,i+1],[j,j+1]).
// 		if (error_flag) {
//  		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse+1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
//  						      Uchange);
//  		}
//  		// If fine grid cell-centre falls outside of the NE coarse
//  		// grid interpolants then try the SW coarse grid solution
//  		// states as the interpolants ([i-1,i],[j-1,j]).
//  		if (error_flag) {
//  		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse-1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
//  						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
//  						      Uchange);
//  		}
//  		// If fine grid cell-centre falls outside of the SW coarse
//  		// grid interpolants then try the N coarse grid solution
//  		// states as the interpolants ([i-1,i],[j,j+1]).
//  		if (error_flag) {
//  		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse+1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
//  						      Uchange);
//  		}
//  		// If none of the four sets of possible interpolants 
//  		// provide a suitable set of interpolants then restart
//  		// the prolongation associated with the current coarse
//  		// grid cell using injection.
//  		if (error_flag) {
//  		  injection_for_this_coarse_cell = true;
//  		  goto restart;
//  		}

//  	      } else {

//  		// Perform the prolongation using injection.
//  		Uchange = Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse];

//  	      }

//  	      // Update the conservative solution state in the SE fine
//  	      // grid cell.  If required, use residual reduction if the
//  	      // update creates a unphysical state (negative density or 
//  	      // pressure).  Residual reduction is performed by halving
//  	      // the residual update.
//  	      if (UPDATE_STABILITY_SWITCH) {
//  		Unew = Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] + Uchange;
//  		if (Unew[1] < ZERO || Unew.p() < ZERO) {
//  		  residual_reduction_number = 1;
//  		  while (residual_reduction_number <= NUMBER_OF_UPDATE_REDUCTIONS) {
//  		    Uchange *= HALF;
//  		    Unew = Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] + Uchange;
//  		    if (Unew[1] > ZERO && Unew.p() > ZERO) break;
//  		    residual_reduction_number++;
//  		  }
//  		  if (residual_reduction_number > NUMBER_OF_UPDATE_REDUCTIONS) Uchange *= ZERO;
//  		}
//  	      }
//  	      Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] += Uchange;

//  	      /**********************************************************
//  	       * Update the cell-centred solution state for the north-  *
//  	       * east fine grid cell.                                   *
//  	       **********************************************************/

//  	      // Determine the (i,j) index of the NE fine grid cell.
//  	      i_fine = 2*(i_coarse-Nghost)+Nghost+1;
//  	      j_fine = 2*(j_coarse-Nghost)+Nghost+1;


//  	      // Attempt to use bilinear interpolation if the injection flag
//  	      // has not been set true.
//  	      if (!injection_for_this_coarse_cell) {

// 		// Attempt to use the NE coarse grid solution states as the
//  		// interpolants ([i,i+1],[j,j+1]).
//  		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
//  						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
//  						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
//  						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
//  						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse+1],
//  						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
//  						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
//  						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
//  						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
//  						    Uchange);
//  		// If fine grid cell-centre falls outside of the NE coarse
//  		// grid interpolants then try the SE coarse grid solution
//  		// states as the interpolants ([i,i+1],[j-1,j]).
//  		if (error_flag) {
//  		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse-1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
//  						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
//  						      Uchange);
//  		}
//  		// If fine grid cell-centre falls outside of the SE coarse
//  		// grid interpolants then try the NW coarse grid solution
//  		// states as the interpolants ([i-1,i],[j,j+1]).
//  		if (error_flag) {
//  		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse+1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
//  						      Uchange);
//  		}
// 		// If fine grid cell-centre falls outside of the NW coarse
//  		// grid interpolants then try the SW coarse grid solution
//  		// states as the interpolants ([i-1,i],[j-1,j]).
//  		if (error_flag) {
//  		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse-1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
//  						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
//  						      Uchange);
//  		}
//  		// If none of the four sets of possible interpolants 
//  		// provide a suitable set of interpolants then restart
//  		// the prolongation associated with the current coarse
//  		// grid cell using injection.
//  		if (error_flag) {
//  		  injection_for_this_coarse_cell = true;
//  		  goto restart;
//  		}

//  	      } else {

//  		// Perform the prolongation using injection.
//  		Uchange = Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse];

//  	      }

//  	      // Update the conservative solution state in the NE fine
//  	      // grid cell.  If required, use residual reduction if the
// 	      // update creates a unphysical state (negative density or 
//  	      // pressure).  Residual reduction is performed by halving
//  	      // the residual update.
//  	      if (UPDATE_STABILITY_SWITCH) {
//  		Unew = Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] + Uchange;
//  		if (Unew[1] < ZERO || Unew.p() < ZERO) {
//  		  residual_reduction_number = 1;
//  		  while (residual_reduction_number <= NUMBER_OF_UPDATE_REDUCTIONS) {
//  		    Uchange *= HALF;
//  		    Unew = Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] + Uchange;
//  		    if (Unew[1] > ZERO && Unew.p() > ZERO) break;
//  		    residual_reduction_number++;
//  		  }
//  		  if (residual_reduction_number > NUMBER_OF_UPDATE_REDUCTIONS) Uchange *= ZERO;
//  		}
//  	      }
//  	      Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] += Uchange;

//  	      /**********************************************************
//  	       * Update the cell-centred solution state for the north-  *
//  	       * west fine grid cell.                                   *
//  	       **********************************************************/

//  	      // Determine the (i,j) index of the NW fine grid cell.
//  	      i_fine = 2*(i_coarse-Nghost)+Nghost;
//  	      j_fine = 2*(j_coarse-Nghost)+Nghost+1;

//  	      // Attempt to use bilinear interpolation if the injection flag
//  	      // has not been set true.
//  	      if (!injection_for_this_coarse_cell) {

// 		// Attempt to use the NW coarse grid solution states as the
//  		// interpolants ([i-1,i],[j,j+1]).
//  		error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
//  						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
//  						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse+1],
//  						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse+1].Xc,
//  						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
//  						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
//  						    Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
//  						    Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
//  						    Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
//  						    Uchange);
//  		// If fine grid cell-centre falls outside of the NW coarse
//  		// grid interpolants then try the NE coarse grid solution
//  		// states as the interpolants ([i,i+1],[j,j+1]).
//  		if (error_flag) {
//  		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse+1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse+1].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse+1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse+1].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
//  						      Uchange);
//  		}
//  		// If fine grid cell-centre falls outside of the NE coarse
//  		// grid interpolants then try the SW coarse grid solution
// 		// states as the interpolants ([i-1,i],[j-1,j]).
//  		if (error_flag) {
//  		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse-1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse-1].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse-1][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse-1][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
//  						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
//  						      Uchange);
//  		}
//  		// If fine grid cell-centre falls outside of the SW coarse
//  		// grid interpolants then try the SE coarse grid solution
//  		// states as the interpolants ([i,i+1],[j-1,j]).
//  		if (error_flag) {
//  		  error_flag = Bilinear_Interpolation(Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse-1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse-1].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse].Xc,
//  						      Local_SolnBlk[Level_Coarse][nb].U[i_coarse+1][j_coarse-1],
//  						      Local_SolnBlk[Level_Coarse][nb].Grid.Cell[i_coarse+1][j_coarse-1].Xc,
//  						      Local_SolnBlk[Level_Fine][nb].Grid.Cell[i_fine][j_fine].Xc,
//  						      Uchange);
//  		}
//  		// If none of the four sets of possible interpolants 
//  		// provide a suitable set of interpolants then restart
//  		// the prolongation associated with the current coarse
//  		// grid cell using injection.
//  		if (error_flag) {
//  		  injection_for_this_coarse_cell = true;
//  		  goto restart;
//  		}

//  	      } else {

//  		// Perform the prolongation using injection.
//  		Uchange = Local_SolnBlk[Level_Coarse][nb].U[i_coarse][j_coarse];

//  	      }

//  	      // Update the conservative solution state in the NW fine
//  	      // grid cell.  If required, use residual reduction if the
//  	      // update creates a unphysical state (negative density or 
//  	      // pressure).  Residual reduction is performed by halving
//  	      // the residual update.
//  	      if (UPDATE_STABILITY_SWITCH) {
//  		Unew = Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] + Uchange;
//  		if (Unew[1] < ZERO || Unew.p() < ZERO) {
//  		  residual_reduction_number = 1;
//  		  while (residual_reduction_number <= NUMBER_OF_UPDATE_REDUCTIONS) {
//  		    Uchange *= HALF;
//  		    Unew = Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] + Uchange;
//  		    if (Unew[1] > ZERO && Unew.p() > ZERO) break;
//  		    residual_reduction_number++;
//  		  }
//  		  if (residual_reduction_number > NUMBER_OF_UPDATE_REDUCTIONS) Uchange *= ZERO;
//  		}
//  	      }
//  	      Local_SolnBlk[Level_Fine][nb].U[i_fine][j_fine] += Uchange;

//  	      // Reset the injection flag to false.
//  	      injection_for_this_coarse_cell = false;

// 	    }

// 	  }

// 	}
//       }

//     }
//   }

//   // Prolongation and update of the solution block was successful.
//   return 0;

// }

/**********************************************************************
 * Routine: Store_Current_Solution_in_uo_MG                           *
 *                                                                    *
 * This routine copies the solution state in each block into the      *
 * array uo_MG                                                        *
 *                                                                    *
 **********************************************************************/
// template <> void FAS_Multigrid2D_Solver<Dusty2D_cState,
// 			                Dusty2D_Quad_Block,
//                                         Dusty2D_Input_Parameters>::
// Store_Current_Solution_in_uo_MG(const int Level) {

//   int Nghost;

//   // Loop through each local solution block.
//   for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
//     if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
//       // Get the number of ghost cells on the current grid level.
//       Nghost = Local_SolnBlk[Level][nb].Nghost;
//       // Loop through and copy the solution content for each cell of the
//       // current grid level.
//       for (int i = Local_SolnBlk[Level][nb].ICl-Nghost; i <= Local_SolnBlk[Level][nb].ICu+Nghost; i++) {
// 	for (int j = Local_SolnBlk[Level][nb].JCl-Nghost; j <= Local_SolnBlk[Level][nb].JCu+Nghost; j++) {
// 	  if (Local_SolnBlk[Level][nb].Grid.cell_status[i][j] == CELL_STATUS_ACTIVE) {
// 	    MG_SolnBlk[Level][nb].uo_MG[i][j] = Local_SolnBlk[Level][nb].U[i][j];
// 	  }
// 	}
//       }
//     }
//   }

// }

/**********************************************************************
 * Routine: Subtract_dUdt_from_P                                      *
 *                                                                    *
 * This routine subtracts dUdt from P.  P is expected to contain the  *
 * restricted fine grid residuals, where dUdt is expected to contain  *
 * the residuals evaluated from the restricted fine grid solution,    *
 * such that after this operation the forcing term P is formed.       *
 *                                                                    *
 **********************************************************************/
// template <> void FAS_Multigrid2D_Solver<Dusty2D_cState,
// 			                Dusty2D_Quad_Block,
//                                         Dusty2D_Input_Parameters>::
// Subtract_dUdt_from_P(const int Level) {

//   // Loop through each local solution block.
//   for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
//     if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
//       // Loop through and subtract residuals evaluated from the
//       // restricted fine grid solution, dUdt, from the forcing term, P,
//       // for each cell of the current grid level.
//       for (int i = Local_SolnBlk[Level][nb].ICl; i <= Local_SolnBlk[Level][nb].ICu; i++) {
// 	for (int j = Local_SolnBlk[Level][nb].JCl; j <= Local_SolnBlk[Level][nb].JCu; j++) {
// 	  if (Local_SolnBlk[Level][nb].Grid.cell_status[i][j] == CELL_STATUS_ACTIVE) {
// 	    MG_SolnBlk[Level][nb].P[i][j] -= Local_SolnBlk[Level][nb].dUdt[i][j][0];
// 	  }
// 	}
//       }
//     }
//   }

// }

/**********************************************************************
 * Routine: Update_Primitive_Variables                                *
 *                                                                    *
 * This routine updates all primitive variables W using the conserved *
 * variables U, for all solution blocks on  a level.                  *
 *                                                                    *
 **********************************************************************/
// template <> void FAS_Multigrid2D_Solver<Dusty2D_cState,
// 			                Dusty2D_Quad_Block,
//                                         Dusty2D_Input_Parameters>::
// Update_Primitive_Variables(const int Level) {

//   // Loop through each local solution block.
//   for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
//     if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
//       // Loop through and determine the primitive solution state for 
//       // each (interior) cell of the current grid level.
//       for (int j = Local_SolnBlk[Level][nb].JCl; j <= Local_SolnBlk[Level][nb].JCu; j++) {
// 	for (int i = Local_SolnBlk[Level][nb].ICl; i <= Local_SolnBlk[Level][nb].ICu; i++) {
// 	  if (Local_SolnBlk[Level][nb].Grid.cell_status[i][j] == CELL_STATUS_ACTIVE) {
// 	    Local_SolnBlk[Level][nb].W[i][j] = W(Local_SolnBlk[Level][nb].U[i][j]);
// 	  }
// 	}
//       }
//     }
//   }

// }

/**********************************************************************
 * Routine: Evaluate_Solution_Changes                                 *
 *                                                                    *
 * This routine subtracts uo_MG[][] from U[][] and store it into      *
 * U[][] for later prolongation to the finer level.                   *
 *                                                                    *
 **********************************************************************/
// template <> void FAS_Multigrid2D_Solver<Dusty2D_cState,
//                                         Dusty2D_Quad_Block,
// 			                Dusty2D_Input_Parameters>::
// Evaluate_Solution_Changes(const int Level) {

//   int Nghost;

//   // Loop through each local solution block.
//   for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
//     if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
//       // Get the number of ghost cells on the current grid level.
//       Nghost = Local_SolnBlk[Level][nb].Nghost;
//       // Determine the solution change in every cell if the current grid.
//       for (int i = Local_SolnBlk[Level][nb].ICl-Nghost; i <= Local_SolnBlk[Level][nb].ICu+Nghost; i++) {
// 	for (int j = Local_SolnBlk[Level][nb].JCl-Nghost; j <= Local_SolnBlk[Level][nb].JCu+Nghost; j++) {
// 	  if (Local_SolnBlk[Level][nb].Grid.cell_status[i][j] == CELL_STATUS_ACTIVE) {
// 	    Local_SolnBlk[Level][nb].U[i][j] -= MG_SolnBlk[Level][nb].uo_MG[i][j];
// 	  }
// 	}
//       }
//     }
//   }

// }

/**********************************************************************
 * Routine: dUdt_Multistage_Explicit_for_Multigrid                    *
 *                                                                    *
 * This routine updates the solution for a 1D array of 2D             *
 * quadrilateral multi-block solution blocks.  A variety of           *
 * multistage explicit time integration and a 2nd-ororder limited     *
 * upwind finite-volume spatial discretization scheme for the         *
 * convective flux coupled with a centrally-weighted finite-volume    *
 * discretization for the diffused flux can be used depending on the  *
 * specified input values.                                            *
 *                                                                    *
 **********************************************************************/
// template <> int FAS_Multigrid2D_Solver<Dusty2D_cState,
//                                        Dusty2D_Quad_Block,
//                                        Dusty2D_Input_Parameters>::
// dUdt_Multistage_Explicit_for_Multigrid(const int Level,
// 				       const int Top_Level,
// 				       const int i_stage) {

//   int error_flag, flux_function_type, limiter_type, k_residual;

//   // Force use of first-order spatial reconstruction and HLLE flux
//   // function on coarse levels.
//   flux_function_type = IP->i_Flux_Function;
//   limiter_type = IP->i_Limiter;
//   if (Level > Top_Level &&
//       FORCE_FIRST_ORDER_SPATIAL_RECONSTRUCTION == true) {
//     //if (Level > Top_Level &&
//     //IP->Multigrid_IP.First_Order_Coarse_Mesh_Reconstruction) {
//     if (IP->Local_Time_Stepping != LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER ||
//         ((IP->Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) &&
//          (IP->i_Flux_Function != FLUX_FUNCTION_ROE_PRECON_WS)) &&
//         ((IP->Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) &&
//          (IP->i_Flux_Function != FLUX_FUNCTION_HLLE_PRECON_WS)))
//       IP->i_Flux_Function = FLUX_FUNCTION_HLLE;
//     IP->i_Limiter = LIMITER_ZERO;
//   }

//   // Temporarily overwrite the time integration type to the multistage
//   // optimally smoothing scheme.
//   IP->i_Time_Integration = IP->Multigrid_IP.i_Smoothing;

//   // Evaluate the residual for each solution block.
//   error_flag = dUdt_Multistage_Explicit(Local_SolnBlk[Level],
// 					List_of_Local_Solution_Blocks[Level],
// 					*IP,
// 					i_stage);
//   if (error_flag) return error_flag;

//   // Apply the defect correction forcing term, if necessary.
//   //error_flag = Apply_the_FAS_Multigrid_Forcing_Term(Level,
//   //					    Top_Level,
//   //					    i_stage);
//   //if (error_flag) return error_flag;

//   // Reset flux function, limiter, and time integration parameters.
//   if (Level > Top_Level &&
//       FORCE_FIRST_ORDER_SPATIAL_RECONSTRUCTION == true) {
//     //if (Level > Top_Level &&
//     //IP->Multigrid_IP.First_Order_Coarse_Mesh_Reconstruction) {
//     IP->i_Flux_Function = flux_function_type;
//     IP->i_Limiter = limiter_type;
//   }
// //   if (!IP->Multigrid_IP.i_Dual_Time_Stepping) {
//   IP->i_Time_Integration = TIME_STEPPING_MULTIGRID;
// //   } else {
// //     IP->i_Time_Integration = TIME_STEPPING_DUAL_TIME_STEPPING;
// //   }

//   // Residuals for each quadrilateral multi-block solution has been
//   // successfully calculated.
//   return 0;

// }

/**********************************************************************
 * Routine: Residual_Evaluation_for_Multigrid                         *
 *                                                                    *
 * This routine evaluates the residual for a 1D array of solution     *
 * blocks given the solution U using an upwind finite-volume spatial  *
 * discretization procedure for hyperbolic fluxes and a centrally-    *
 * weighted finite-volume spatial discretization procedure for the    *
 * elliptic fluxes.                                                   *
 *                                                                    *
 **********************************************************************/
// template <> int FAS_Multigrid2D_Solver<Dusty2D_cState,
// 			               Dusty2D_Quad_Block,
//                                        Dusty2D_Input_Parameters>::
// Residual_Evaluation_for_Multigrid(const int Level,
// 				  const int Top_Level,
// 				  const int apply_forcing_term_flag) {

//   int error_flag, flux_function_type, limiter_type;

//   // Force use of first-order spatial reconstruction and HLLE flux
//   // function flux function on coarse levels (better damping).
//   flux_function_type = IP->i_Flux_Function;
//   limiter_type = IP->i_Limiter;
//   if (Level > Top_Level &&
//       FORCE_FIRST_ORDER_SPATIAL_RECONSTRUCTION == true) {
//      //  if (Level > Top_Level &&
//      //IP->Multigrid_IP.First_Order_Coarse_Mesh_Reconstruction) {
//     if (IP->Local_Time_Stepping != LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER ||
//         ((IP->Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) &&
//          (IP->i_Flux_Function != FLUX_FUNCTION_ROE_PRECON_WS)) &&
//         ((IP->Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) &&
//          (IP->i_Flux_Function != FLUX_FUNCTION_HLLE_PRECON_WS)))
//        IP->i_Flux_Function = FLUX_FUNCTION_HLLE;
//     IP->i_Limiter = LIMITER_ZERO;
//   }

//   // Evaluate the residual for each solution block.
//   for (int nb = 0 ; nb <= List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
//     if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
//       // Evaluate residual.
//       error_flag = dUdt_Residual_Evaluation(Local_SolnBlk[Level][nb],
// 					    *(IP));
//       if (error_flag) return error_flag;
// //       if (IP->Multigrid_IP.i_Dual_Time_Stepping) {
// // 	error_flag = dUdtau_Residual_Evaluation(Level);
// // 	if (error_flag) return error_flag;
// //       }
//     }
//   }

//   // Reset flux function and limiter.
//   if (Level > Top_Level &&
//       FORCE_FIRST_ORDER_SPATIAL_RECONSTRUCTION == true) {
//     //if (Level > Top_Level &&
//     //IP->Multigrid_IP.First_Order_Coarse_Mesh_Reconstruction) {
//     IP->i_Flux_Function = flux_function_type;
//     IP->i_Limiter = limiter_type;
//   }

//   // Send boundary flux corrections at block interfaces with resolution changes.
//   error_flag = Send_Conservative_Flux_Corrections(Local_SolnBlk[Level],
// 						  List_of_Local_Solution_Blocks[Level],
// 						  Local_SolnBlk[Level][0].NumVar());
//   if (error_flag) return error_flag;

//   // Apply the boundary flux corrections to ensure that method is
//   // conservative.
//   Apply_Boundary_Flux_Corrections(Local_SolnBlk[Level],
// 				  List_of_Local_Solution_Blocks[Level]);

//   for (int nb = 0; nb < List_of_Local_Solution_Blocks[Level].Nblk; nb++) {
//     if (List_of_Local_Solution_Blocks[Level].Block[nb].used == ADAPTIVEBLOCK2D_USED) {
//       // Add the forcing term if required.
//       if (apply_forcing_term_flag) {
// 	for (int j = Local_SolnBlk[Level][nb].JCl; j <= Local_SolnBlk[Level][nb].JCu; j++) {
// 	  for (int i = Local_SolnBlk[Level][nb].ICl; i <= Local_SolnBlk[Level][nb].ICu; i++) {
// 	    Local_SolnBlk[Level][nb].dUdt[i][j][0] += MG_SolnBlk[Level][nb].P[i][j];
// 	  }
// 	}
//       }
//     }
//   }

//   // Residuals successfully calculated for quadrilateral multi-block solution blocks.
//   return 0;

// }

/**********************************************************************
 * Routine: Smooth                                                    *
 *                                                                    *
 * This routine implements the multigrid smoothing for a 1D array of  *
 * 2D quadrilateral multi-block solution blocks at the specified grid *
 * level.                                                             *
 *                                                                    *
 **********************************************************************/
// template <> int FAS_Multigrid2D_Solver<Dusty2D_cState,
// 			               Dusty2D_Quad_Block,
//                                        Dusty2D_Input_Parameters>::
// Smooth(const int Level,
//        const int Top_Level) {

//   int error_flag;
//   double dTime;

//   // Determine the local (and minimum) time-step for each cell of the 
//   // current grid level.
//   dTime = CFL(Local_SolnBlk[Level],
// 	      List_of_Local_Solution_Blocks[Level],
//               *(IP));
	
//   // Set each cell to the global minimum time-step if required.
//   if (!IP->Local_Time_Stepping) {
//      dTime = CFFC_Minimum_MPI(dTime);
//      Set_Global_TimeStep(Local_SolnBlk[Level],
// 			 List_of_Local_Solution_Blocks[Level],
// 			 dTime);
//   }

//   // If the current grid level is not the top grid level then ensure
//   // that the coarse grid time-step is less than or equal to the 
//   // time-steps of the associated finer grid level.
//   if (Level != Top_Level) CFL_Multigrid(Level);

//   // Smooth the solution using N stage multistage time stepping scheme.
//   for (int i_stage = 1; i_stage <= IP->N_Stage; i_stage++) {

//     // 1. Exchange solution information between neighbouring blocks.
//     error_flag = Exchange_Solution_Information(Level,
// 					       OFF);
//     if (error_flag) {
//       cout << "\n FASMultigrid2D ERROR: Message passing error on processor "
// 	   << List_of_Local_Solution_Blocks[Level].ThisCPU
// 	   << ".\n";
//       cout.flush();
//     }
//     error_flag = CFFC_OR_MPI(error_flag);
//     if (error_flag) return error_flag;

//     // 2. Apply boundary conditions for stage.
//     if (Level == Top_Level || APPLY_BOUNDARY_CONDITIONS_ON_COARSE_MESH == true) {
//       //if (Level == Top_Level ||
//       //IP->Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions) {
//       BCs(Local_SolnBlk[Level],
// 	  List_of_Local_Solution_Blocks[Level],
// 	  *IP);
//     }

//     // 3. Determine solution residuals for stage.
//     error_flag = dUdt_Multistage_Explicit_for_Multigrid(Level,
// 							Top_Level,
// 							i_stage);
//     if (error_flag) {
//       cout << "\n FASMultigrid2D ERROR: Residual calculation error on processor "
// 	   << List_of_Local_Solution_Blocks[Level].ThisCPU
// 	   << ".\n";
//       cout.flush();
//     }
//     error_flag = CFFC_OR_MPI(error_flag);
//     if (error_flag) return error_flag;

//     // 4. Send boundary flux corrections at block interfaces with
//     // resolution changes.
//     error_flag = Send_Conservative_Flux_Corrections(Local_SolnBlk[Level],
// 						    List_of_Local_Solution_Blocks[Level],
// 						    Local_SolnBlk[Level][0].NumVar());
//     if (error_flag) {
//       cout << "\n FASMultigrid2D ERROR: Flux correction message passing error on processor "
// 	   << List_of_Local_Solution_Blocks[Level].ThisCPU
// 	   << ".\n";
//       cout.flush();
//     }
//     error_flag = CFFC_OR_MPI(error_flag);
//     if (error_flag) return error_flag;

//     // 5. Apply boundary flux corrections to ensure that method is
//     // conservative.
//     Apply_Boundary_Flux_Corrections_Multistage_Explicit_for_Multigrid(Level,
// 								      Top_Level,
// 								      i_stage);

//     // 6. Apply the defect correction forcing term, if necessary.
//     error_flag = Apply_the_FAS_Multigrid_Forcing_Term(Level,
// 						      Top_Level,
// 						      i_stage);
//     if (error_flag) {
//       cout << "\n FASMultigrid2D ERROR: Forcing term application error on processor "
// 	   << List_of_Local_Solution_Blocks[Level].ThisCPU
// 	   << ".\n";
//       cout.flush();
//     }
//     error_flag = CFFC_OR_MPI(error_flag);
//     if (error_flag) return error_flag;

//     // 7. Smooth the solution residual using implicit residual smoothing.
//     if (IP->Residual_Smoothing) {
//        Residual_Smoothing(Local_SolnBlk[Level],
// 	     	          List_of_Local_Solution_Blocks[Level],
// 		          *IP,
// 		          i_stage);
//     }

//     // 8. Update solution for the current stage.
//     error_flag = Update_Solution_Multistage_Explicit_for_Multigrid(Level,
// 								   Top_Level,
// 								   i_stage);
//     if (error_flag) {
//       cout << "\n FASMultigrid2D ERROR: Solution update error on processor "
// 	   << List_of_Local_Solution_Blocks[Level].ThisCPU
// 	   << ".\n";
//       cout.flush();
//     }
//     error_flag = CFFC_OR_MPI(error_flag);
//     if (error_flag) return error_flag;

//   }

//   // Exchange solution information between neighbouring blocks
//   // after smoothing is complete.
//   error_flag = Exchange_Solution_Information(Level,
//          				     OFF);
//   if (error_flag) {
//     cout << "\n FASMultigrid2D ERROR: Message passing error on processor "
// 	 << List_of_Local_Solution_Blocks[Level].ThisCPU
// 	 << ".\n";
//     cout.flush();
//   }
//   error_flag = CFFC_OR_MPI(error_flag);
//   if (error_flag) return error_flag;

//   // Apply boundary conditions after smoothing is complete.
//   if (Level == Top_Level || APPLY_BOUNDARY_CONDITIONS_ON_COARSE_MESH == true) {
//     //if (Level == Top_Level ||
//     //IP->Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions) {
//     BCs(Local_SolnBlk[Level],
// 	List_of_Local_Solution_Blocks[Level],
// 	*IP);
//   }

//   // Determine solution residuals for the entire stage.
//   error_flag = Residual_Evaluation_for_Multigrid(Level,
// 						 Top_Level,
//         		 	 		 Level-Top_Level);
//   if (error_flag) {
//     cout << "\n FASMultigrid2D ERROR: Residual calculation error on processor "
// 	 << List_of_Local_Solution_Blocks[Level].ThisCPU
// 	 << ".\n";
//     cout.flush();
//   }
//   error_flag = CFFC_OR_MPI(error_flag);
//   if (error_flag) return error_flag;

//   // Multigrid smoothing applied successfully.
//   return 0;

// }

#endif // _DUSTY2D_QUAD_MULTIGRID_INCLUDED
