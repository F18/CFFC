/**********************************************************************
 * Electrostatic2DQuadIOMultiBlock.cc                                 *
 *                                                                    *
 * Multi-block versions of input and output subroutines for 2D        *
 * Electrostatic multi-block quadrilateral mesh solution classes.     *
 *                                                                    *
 **********************************************************************/

// Include 2D Electrostatic quadrilateral mesh solution header file.

#ifndef _ELECTROSTATIC2D_QUAD_INCLUDED
#include "Electrostatic2DQuad.h"
#endif // _ELECTROSTATIC2D_QUAD_INCLUDED

/**********************************************************************
 * Electrostatic2D_Quad_Block -- IO Multiple Block External           *
 *                               Subroutines.                         *
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
int Read_Restart_Solution(Electrostatic2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          Electrostatic2D_Input_Parameters &IP,
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
    if (IP.Restart_File_Name[i] == ' ' ||
	IP.Restart_File_Name[i] == '.') break;
    prefix[i] = IP.Restart_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Restart_File_Name)) break;
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
	IP.Maximum_Number_of_Time_Steps += Number_of_Time_Steps;
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
int Write_Restart_Solution(Electrostatic2D_Quad_Block *Soln_ptr,
                           AdaptiveBlock2D_List &Soln_Block_List,
                           Electrostatic2D_Input_Parameters &IP,
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
    if (IP.Restart_File_Name[i] == ' ' ||
	IP.Restart_File_Name[i] == '.') break;
    prefix[i] = IP.Restart_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Restart_File_Name)) break;
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
      //restart_file << IP.Gas_Type << "\n";
      //restart_file << IP.Solid_Type << "\n";
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
int Output_Tecplot(Electrostatic2D_Quad_Block *Soln_ptr,
                   AdaptiveBlock2D_List &Soln_Block_List,
                   Electrostatic2D_Input_Parameters &IP,
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
    if (IP.Output_File_Name[i] == ' ' || IP.Output_File_Name[i] == '.') break;
    prefix[i] = IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name)) break;
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
		     IP,
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
int Output_Cells_Tecplot(Electrostatic2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List,
                         Electrostatic2D_Input_Parameters &IP,
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
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i] = IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name)) break;
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
			   IP,
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
 * Routine: Output_Quasi3D_Tecplot                                    *
 *                                                                    *
 * Writes the solution values at nodes of the positive cells for a 1D *
 * array of 2D quadrilateral multi-block solution blocks to the       *
 * specified output data file(s) in a quasi-3D format suitable for    *
 * plotting a with TECPLOT.  Returns a non-zero value if cannot write *
 * any of the TECPLOT solution files.                                 *
 *                                                                    *
 **********************************************************************/
int Output_Quasi3D_Tecplot(Electrostatic2D_Quad_Block *Soln_ptr,
			   AdaptiveBlock2D_List &Soln_Block_List,
			   Electrostatic2D_Input_Parameters &IP,
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
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i] = IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_quasi3D_cpu");

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
      Output_Quasi3D_Tecplot(Soln_ptr[nb],
			     IP,
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
 * Routine: Output_Mesh_Tecplot                                       *
 *                                                                    *
 * Writes the nodes of the mesh for a 1D array of 2D quadrilateral    *
 * multi-block solution blocks to the specified output data file(s)   *
 * in a format suitable for plotting with TECPLOT. Returns a non-zero *
 * value if cannot write any of the TECPLOT solution files.           *
 *                                                                    *
 **********************************************************************/
int Output_Mesh_Tecplot(Electrostatic2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        Electrostatic2D_Input_Parameters &IP,
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
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i] = IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name)) break;
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
int Output_Mesh_Gnuplot(Electrostatic2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        Electrostatic2D_Input_Parameters &IP,
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
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i] = IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name)) break;
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
