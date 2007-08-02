/**********************************************************************
 * LevelSet2DQuadIOMultiBlock.cc                                      *
 *                                                                    *
 * Multi-block versions of the input and output subroutines for the   *
 * 2D Level Set multi-block quadrilateral mesh solution class.        *
 *                                                                    *
 **********************************************************************/

// Include 2D LevelSet quadrilateral mesh solution header file.

#ifndef _LEVELSET2D_QUAD_INCLUDED
#include "LevelSet2DQuad.h"
#endif // _LEVELSET2D_QUAD_INCLUDED

/**********************************************************************
 * LevelSet2D_Quad_Block -- Input and Output Multiple Block External  *
 *                          Subroutines.                              *
 **********************************************************************/

/**********************************************************************
 * Routine: Read_Restart_Solution                                     *
 *                                                                    *
 * Reads restart solution file(s) and assigns values to the solution  *
 * variables of a 1D array of 2D quadrilateral multi-block solution   *
 * blocks.  Returns a non-zero value if cannot read any of the        *
 * restart solution files.                                            *
 *                                                                    *
 **********************************************************************/
int Read_Restart_Solution(LevelSet2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          LevelSet2D_Input_Parameters &Input_Parameters,
                          int &Number_of_Time_Steps,
                          double &Time,
                          CPUTime &CPU_Time) {

  int i, i_new_time_set, nsteps;
  char prefix[256], extension[256], restart_file_name[256];
  char *restart_file_name_ptr;
  ifstream restart_file;
  double time0;
  CPUTime cpu_time0;
  
  // Determine prefix of restart file names.
  i = 0;
  while (1) {
    if (Input_Parameters.Restart_File_Name[i] == ' ' ||
	Input_Parameters.Restart_File_Name[i] == '.') break;
    prefix[i]=Input_Parameters.Restart_File_Name[i];
    i = i + 1;
    if (i > strlen(Input_Parameters.Restart_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_blk");

  // Read the initial data for each solution block.
  i_new_time_set = 0;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      // Restart file name base on global block number.
      sprintf(extension,"%.6d",Soln_Block_List.Block[nb].gblknum);
      strcat(extension,".soln");
      strcpy(restart_file_name,prefix);
      strcat(restart_file_name,extension);
      restart_file_name_ptr = restart_file_name;

      // Open restart file.
      restart_file.open(restart_file_name_ptr,ios::in);
      if (restart_file.bad()) return 1;

      // Read solution block data.
      restart_file.setf(ios::skipws);
      restart_file >> nsteps >> time0 >> cpu_time0;
      restart_file.unsetf(ios::skipws);
      if (!i_new_time_set) {
	Number_of_Time_Steps = nsteps;
	Input_Parameters.Maximum_Number_of_Time_Steps += Number_of_Time_Steps;
	Time = time0;
	CPU_Time.cput = cpu_time0.cput;
	i_new_time_set = 1;
      }
      restart_file >> Soln_ptr[nb];

      // Close restart file.
      restart_file.close();
    }
  }

  // Reading of restart files complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Write_Restart_Solution                                    *
 *                                                                    *
 * Writes restart solution file(s) for a 1D array of 2D quadrilateral *
 * multi-block solution blocks. Returns a non-zero value if cannot    *
 * write any of the restart solution files.                           *
 *                                                                    *
 **********************************************************************/
int Write_Restart_Solution(LevelSet2D_Quad_Block *Soln_ptr,
                           AdaptiveBlock2D_List &Soln_Block_List,
                           LevelSet2D_Input_Parameters &Input_Parameters,
                           const int Number_of_Time_Steps,
                           const double &Time,
                           const CPUTime &CPU_Time) {

  int i;
  char prefix[256], extension[256], restart_file_name[256];
  char *restart_file_name_ptr;
  ofstream restart_file;

  // Determine prefix of restart file names.
  i = 0;
  while (1) {
    if (Input_Parameters.Restart_File_Name[i] == ' ' ||
	Input_Parameters.Restart_File_Name[i] == '.') break;
    prefix[i]=Input_Parameters.Restart_File_Name[i];
    i = i + 1;
    if (i > strlen(Input_Parameters.Restart_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_blk");

  // Write the solution data for each solution block.
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      // Restart file name base on global block number.
      sprintf(extension,"%.6d",Soln_Block_List.Block[nb].gblknum);
      strcat(extension,".soln");
      strcpy(restart_file_name,prefix);
      strcat(restart_file_name,extension);
      restart_file_name_ptr = restart_file_name;

      // Open restart file.
      restart_file.open(restart_file_name_ptr,ios::out);
      if (restart_file.bad()) return 1;

      // Write solution block data.
      restart_file.setf(ios::scientific);
      restart_file << setprecision(14) << Number_of_Time_Steps 
		   << " " << Time << " " << CPU_Time << "\n";
      restart_file.unsetf(ios::scientific);
      restart_file << setprecision(14) << Soln_ptr[nb];

      // Close restart file.
      restart_file.close();
    }
  }

  // Writing of restart files complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Output_Tecplot                                            *
 *                                                                    *
 * Writes the solution values at the nodes for a 1D array of 2D       *
 * quadrilateral multi-block solution blocks to the specified output  *
 * data file(s) in a format suitable for plotting with TECPLOT.       *
 * Returns a non-zero value if cannot write any of the TECPLOT        *
 * solution files.                                                    *
 *                                                                    *
 **********************************************************************/
int Output_Tecplot(LevelSet2D_Quad_Block *Soln_ptr,
                   AdaptiveBlock2D_List &Soln_Block_List,
                   LevelSet2D_Input_Parameters &Input_Parameters,
                   const int Number_of_Time_Steps,
                   const double &Time) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;    

  // Only create output file if the number of blocks used by the 
  // current processor is greater than zero.
  if (Soln_Block_List.Nused() == 0) return 0;

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (Input_Parameters.Output_File_Name[i] == ' ' ||
	Input_Parameters.Output_File_Name[i] == '.') break;
    prefix[i]=Input_Parameters.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(Input_Parameters.Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_cpu");

  // Determine output data file name for this processor.
  sprintf(extension,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;

  // Open the output data file.
  output_file.open(output_file_name_ptr,ios::out);
  if (output_file.bad()) return 1;

  // Write the solution data for each solution block.
  i_output_title = 1;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Output_Tecplot(Soln_ptr[nb],
		     Number_of_Time_Steps,
		     Time,
		     Soln_Block_List.Block[nb].gblknum,
		     i_output_title,
		     output_file);
      if (i_output_title) i_output_title = 0;
    }
  }

  // Close the output data file.
  output_file.close();

  // Writing of output data files complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Output_Cells_Tecplot                                      *
 *                                                                    *
 * Writes the cell centred solution values for a 1D array of 2D       *
 * quadrilateral multi-block solution blocks to the specified output  *
 * data file(s) in a format suitable for plotting with TECPLOT.       *
 * Returns a non-zero value if cannot write any of the TECPLOT        *
 * solution files.                                                    *
 *                                                                    *
 **********************************************************************/
int Output_Cells_Tecplot(LevelSet2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List,
                         LevelSet2D_Input_Parameters &Input_Parameters,
                         const int Number_of_Time_Steps,
                         const double &Time) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;    

  // Only create output file if the number of blocks used by the 
  // current processor is greater than zero.
  if (Soln_Block_List.Nused() == 0) return 0;

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (Input_Parameters.Output_File_Name[i] == ' ' ||
	Input_Parameters.Output_File_Name[i] == '.') break;
    prefix[i]=Input_Parameters.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(Input_Parameters.Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_cells_cpu");

  // Determine output data file name for this processor.
  sprintf(extension,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;

  // Open the output data file.
  output_file.open(output_file_name_ptr,ios::out);
  if (output_file.bad()) return 1;

  // Write the solution data for each solution block.
  i_output_title = 1;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Output_Cells_Tecplot(Soln_ptr[nb],
			   Number_of_Time_Steps,
			   Time,
			   Soln_Block_List.ThisCPU,
			   Soln_Block_List.Block[nb].gblknum,
			   i_output_title,
			   output_file);
      if (i_output_title) i_output_title = 0;
    }
  }

  // Close the output data file.
  output_file.close();

  // Writing of output data files complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Output_Nodes_Tecplot                                      *
 *                                                                    *
 * Writes the node location values for a 1D array of 2D quadrilateral *
 * multi-block solution blocks to the specified output data file(s)   *
 * in a format suitable for plotting with TECPLOT.  Returns a non-    *
 * zero value if cannot write any of the TECPLOT solution files.      *
 *                                                                    *
 **********************************************************************/
int Output_Nodes_Tecplot(LevelSet2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List,
                         LevelSet2D_Input_Parameters &Input_Parameters) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;    

  // Only create output file if the number of blocks used by the 
  // current processor is greater than zero.
  if (Soln_Block_List.Nused() == 0) return 0;

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (Input_Parameters.Output_File_Name[i] == ' ' ||
	Input_Parameters.Output_File_Name[i] == '.') break;
    prefix[i]=Input_Parameters.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(Input_Parameters.Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_nodes_cpu");

  // Determine output data file name for this processor.
  sprintf(extension,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;

  // Open the output data file.
  output_file.open(output_file_name_ptr,ios::out);
  if (output_file.bad()) return 1;

  // Write the solution data for each solution block.
  i_output_title = 1;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Output_Nodes_Tecplot(Soln_ptr[nb],
			   Soln_Block_List.Block[nb].gblknum,
			   i_output_title,
			   output_file);
      if (i_output_title) i_output_title = 0;
    }
  }

  // Close the output data file.
  output_file.close();

  // Writing of output data files complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Output_Mesh_Tecplot                                       *
 *                                                                    *
 * Writes the nodes of the mesh for a 1D array of 2D quadrilateral    *
 * multi-block solution blocks to the specified output data file(s)   *
 * in a format suitable for plotting with TECPLOT. Returns a non-zero *
 * value if cannot write any of the TECPLOT solution files.           *
 *                                                                    *
 **********************************************************************/
int Output_Mesh_Tecplot(LevelSet2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        LevelSet2D_Input_Parameters &Input_Parameters,
                        const int Number_of_Time_Steps,
                        const double &Time) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;    

  // Only create output file if the number of blocks used by the 
  // current processor is greater than zero.
  if (Soln_Block_List.Nused() == 0) return 0;

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (Input_Parameters.Output_File_Name[i] == ' ' ||
	Input_Parameters.Output_File_Name[i] == '.') break;
    prefix[i]=Input_Parameters.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(Input_Parameters.Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_mesh_cpu");

  // Determine output data file name for this processor.
  sprintf(extension,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;

  // Open the output data file.
  output_file.open(output_file_name_ptr,ios::out);
  if (output_file.bad()) return 1;

  // Write the solution data for each solution block.
  i_output_title = 1;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Output_Tecplot(Soln_ptr[nb].Grid,
		     Soln_Block_List.Block[nb].gblknum,
		     i_output_title,
		     output_file);
      if (i_output_title) i_output_title = 0;
    }
  }
  
  // Close the output data file.
  output_file.close();

  // Writing of output data files complete.  Return zero value.
  return 0;
  
}

/**********************************************************************
 * Routine: Output_Mesh_Gnuplot                                       *
 *                                                                    *
 * Writes the nodes of the mesh for a 1D array of 2D quadrilateral    *
 * multi-block solution blocks to the specified output data file(s)   *
 * in a format suitable for plotting with GNUPLOT. Returns a non-zero *
 * value if cannot write any of the GNUPLOT solution files.           *
 *                                                                    *
 **********************************************************************/
int Output_Mesh_Gnuplot(LevelSet2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        LevelSet2D_Input_Parameters &Input_Parameters,
                        const int Number_of_Time_Steps,
                        const double &Time) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;    

  // Only create output file if the number of blocks used by the 
  // current processor is greater than zero.
  if (Soln_Block_List.Nused() == 0) return 0;

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (Input_Parameters.Output_File_Name[i] == ' ' ||
	Input_Parameters.Output_File_Name[i] == '.') break;
    prefix[i]=Input_Parameters.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(Input_Parameters.Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_mesh_cpu");

  // Determine output data file name for this processor.
  sprintf(extension,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;

  // Open the output data file.
  output_file.open(output_file_name_ptr,ios::out);
  if (output_file.bad()) return 1;

  // Write the solution data for each solution block.
  i_output_title = 1;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Output_Gnuplot(Soln_ptr[nb].Grid,
		     Soln_Block_List.Block[nb].gblknum,
		     i_output_title,
		     output_file);
      if (i_output_title) i_output_title = 0;
    }
  }

  // Close the output data file.
  output_file.close();

  // Writing of output data files complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Output_Interface_Tecplot                                  *
 *                                                                    *
 * Output the interface spline in a format suitabile for plotting     *
 * with TECPLOT.  Returns a non-zero value if cannot write any of the *
 * TECPLOT solution files.                                            *
 *                                                                    *
 **********************************************************************/
int Output_Interface_Tecplot(LevelSet2D_Quad_Block *Soln_ptr,
			     AdaptiveBlock2D_List &Soln_Block_List,
			     LevelSet2D_Input_Parameters &Input_Parameters,
			     const int Number_of_Time_Steps,
			     const double &Time) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;    
  
  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (Input_Parameters.Output_File_Name[i] == ' ' ||
	Input_Parameters.Output_File_Name[i] == '.') break;
    prefix[i]=Input_Parameters.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(Input_Parameters.Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_interface_list_cpu");
  
  // Determine output data file name for this processor.
  sprintf(extension,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;
  
  // Open the output data file.
  output_file.open(output_file_name_ptr,ios::out);
  if (output_file.bad()) return 1;

  // Write the solution data for each solution block.
  i_output_title = 1;
  //for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
  for (int nb = 0; nb < 1; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Output_Interface_Tecplot(Soln_ptr[nb],
			       Number_of_Time_Steps,
			       Time,
			       Soln_Block_List.Block[nb].gblknum,
			       i_output_title,
			       output_file);
      if (i_output_title) i_output_title = 0;
    }
  }

  // Close the output data file.
  output_file.close();

  // Writing of output data files complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Output_Circle_Tecplot                                     *
 *                                                                    *
 * Output a comparison of the exact and computed solutions for a      *
 * circle interface in a format suitabile for plotting with TECPLOT.  *
 *                                                                    *
 **********************************************************************/
int Output_Circle_Tecplot(LevelSet2D_Quad_Block *Soln_ptr,
			  AdaptiveBlock2D_List &Soln_Block_List,
			  LevelSet2D_Input_Parameters &Input_Parameters,
			  const int Number_of_Time_Steps,
			  const double &Time) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;
  int n_cells;
  double l1_norm, l2_norm, max_norm;
  double l1_temp, l2_temp, max_temp;

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (Input_Parameters.Output_File_Name[i] == ' ' ||
	Input_Parameters.Output_File_Name[i] == '.') break;
    prefix[i]=Input_Parameters.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(Input_Parameters.Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_circle_cpu");
  
  // Determine output data file name for this processor.
  sprintf(extension,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;
  
  // Open the output data file.
  output_file.open(output_file_name_ptr,ios::out);
  if (output_file.bad()) return 1;

  // Write the solution data for each solution block.
  i_output_title = 1;
  l1_norm = ZERO; l2_norm = ZERO; max_norm = ZERO; n_cells = 0;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO;
      Output_Circle_Tecplot(Soln_ptr[nb],
			    Input_Parameters,
			    Time,
			    Soln_Block_List.Block[nb].gblknum,
			    i_output_title,
			    output_file,
			    l1_temp,l2_temp,max_temp);
      l1_norm += l1_temp;
      l2_norm += l2_temp;
      max_norm = max(max_norm,max_temp);
      n_cells += (Soln_ptr[nb].JCu-Soln_ptr[nb].JCl+1)*(Soln_ptr[nb].ICu-Soln_ptr[nb].ICl+1);
      if (i_output_title) i_output_title = 0;
    }
  }

  // Close the output data file.
  output_file.close();

#ifdef _MPI_VERSION
  l1_norm = CFDkit_Summation_MPI(l1_norm);
  l2_norm = CFDkit_Summation_MPI(l2_norm);
  max_norm = CFDkit_Maximum_MPI(max_norm);
  n_cells = CFDkit_Summation_MPI(n_cells);
#endif

  // Calculate the L1-norm and L2-norm for all blocks.
  l1_norm /= double(n_cells);
  l2_norm = sqrt(l2_norm/double(n_cells));

  if (CFDkit_Primary_MPI_Processor()) {
    cout << endl
	 << endl
	 << " ==================================================================== "
	 << endl
	 << " Error norms for the level set solution of a circle interface:"
	 << endl
	 << "   L1_Norm = " << setprecision(14) << l1_norm
	 << endl
	 << "   L2_Norm = " << setprecision(14) << l2_norm
	 << endl
	 << "   Max_Norm = " << setprecision(14) << max_norm
	 << endl
	 << " ==================================================================== "
	 << endl;
  }

  // Writing of output data files complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Output_Ellipse_Tecplot                                    *
 *                                                                    *
 * Output a comparison of the exact and computed solutions for an     *
 * ellipse interface in a format suitabile for plotting with TECPLOT. *
 *                                                                    *
 **********************************************************************/
int Output_Ellipse_Tecplot(LevelSet2D_Quad_Block *Soln_ptr,
			   AdaptiveBlock2D_List &Soln_Block_List,
			   LevelSet2D_Input_Parameters &Input_Parameters,
			   const int Number_of_Time_Steps,
			   const double &Time) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;
  int n_cells;
  double l1_norm, l2_norm, max_norm;
  double l1_temp, l2_temp, max_temp;

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (Input_Parameters.Output_File_Name[i] == ' ' ||
	Input_Parameters.Output_File_Name[i] == '.') break;
    prefix[i]=Input_Parameters.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(Input_Parameters.Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_ellipse_cpu");
  
  // Determine output data file name for this processor.
  sprintf(extension,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;
  
  // Open the output data file.
  output_file.open(output_file_name_ptr,ios::out);
  if (output_file.bad()) return 1;

  // Write the solution data for each solution block.
  i_output_title = 1;
  l1_norm = ZERO; l2_norm = ZERO; max_norm = ZERO; n_cells = 0;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO;
      Output_Ellipse_Tecplot(Soln_ptr[nb],
			     Input_Parameters,
			     Time,
			     Soln_Block_List.Block[nb].gblknum,
			     i_output_title,
			     output_file,
			     l1_temp,l2_temp,max_temp);
      l1_norm += l1_temp;
      l2_norm += l2_temp;
      max_norm = max(max_norm,max_temp);
      n_cells += (Soln_ptr[nb].JCu-Soln_ptr[nb].JCl+1)*(Soln_ptr[nb].ICu-Soln_ptr[nb].ICl+1);
      if (i_output_title) i_output_title = 0;
    }
  }

  // Close the output data file.
  output_file.close();

#ifdef _MPI_VERSION
  l1_norm = CFDkit_Summation_MPI(l1_norm);
  l2_norm = CFDkit_Summation_MPI(l2_norm);
  max_norm = CFDkit_Maximum_MPI(max_norm);
  n_cells = CFDkit_Summation_MPI(n_cells);
#endif

  // Calculate the L1-norm and L2-norm for all blocks.
  l1_norm /= double(n_cells);
  l2_norm = sqrt(l2_norm/double(n_cells));

  if (CFDkit_Primary_MPI_Processor()) {
    cout << endl
	 << endl
	 << " ==================================================================== "
	 << endl
	 << " Error norms for the level set solution of an ellipse interface:"
	 << endl
	 << "   L1_Norm = " << l1_norm
	 << endl
	 << "   L2_Norm = " << l2_norm
	 << endl
	 << "   Max_Norm = " << max_norm
	 << endl
	 << " ==================================================================== "
	 << endl;
  }

  // Writing of output data files complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Output_Zalesaks_Disk_Tecplot                              *
 *                                                                    *
 * Output a comparison of the exact and computed solutions for        *
 * Zalesak's disk in a format suitabile for plotting with TECPLOT.    *
 *                                                                    *
 **********************************************************************/
int Output_Zalesaks_Disk_Tecplot(LevelSet2D_Quad_Block *Soln_ptr,
				 AdaptiveBlock2D_List &Soln_Block_List,
				 LevelSet2D_Input_Parameters &Input_Parameters,
				 const int Number_of_Time_Steps,
				 const double &Time) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;
  int n_cells;
  double l1_norm, l2_norm, max_norm;
  double l1_temp, l2_temp, max_temp;

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (Input_Parameters.Output_File_Name[i] == ' ' ||
	Input_Parameters.Output_File_Name[i] == '.') break;
    prefix[i]=Input_Parameters.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(Input_Parameters.Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_zalesak_cpu");
  
  // Determine output data file name for this processor.
  sprintf(extension,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;
  
  // Open the output data file.
  output_file.open(output_file_name_ptr,ios::out);
  if (output_file.bad()) return 1;

  // Write the solution data for each solution block.
  i_output_title = 1;
  l1_norm = ZERO; l2_norm = ZERO; max_norm = ZERO; n_cells = 0;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO;
      Output_Zalesaks_Disk_Tecplot(Soln_ptr[nb],
				   Input_Parameters,
				   Time,
				   Soln_Block_List.Block[nb].gblknum,
				   i_output_title,
				   output_file,
				   l1_temp,l2_temp,max_temp);
      l1_norm += l1_temp;
      l2_norm += l2_temp;
      max_norm = max(max_norm,max_temp);
      n_cells += (Soln_ptr[nb].JCu-Soln_ptr[nb].JCl+1)*(Soln_ptr[nb].ICu-Soln_ptr[nb].ICl+1);
      if (i_output_title) i_output_title = 0;
    }
  }

  // Close the output data file.
  output_file.close();

#ifdef _MPI_VERSION
  l1_norm = CFDkit_Summation_MPI(l1_norm);
  l2_norm = CFDkit_Summation_MPI(l2_norm);
  max_norm = CFDkit_Maximum_MPI(max_norm);
  n_cells = CFDkit_Summation_MPI(n_cells);
#endif

  // Calculate the L1-norm and L2-norm for all blocks.
  l1_norm /= double(n_cells);
  l2_norm = sqrt(l2_norm/double(n_cells));

  if (CFDkit_Primary_MPI_Processor()) {
    cout << endl
	 << endl
	 << " ==================================================================== "
	 << endl
	 << " Error norms for the level set solution of Zalesak's disk:"
	 << endl
	 << "   L1_Norm = " << l1_norm
	 << endl
	 << "   L2_Norm = " << l2_norm
	 << endl
	 << "   Max_Norm = " << max_norm
	 << endl
	 << " ==================================================================== "
	 << endl;
  }

  // Writing of output data files complete.  Return zero value.
  return 0;

}
