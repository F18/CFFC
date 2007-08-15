/**********************************************************************
 * HighTemp2DQuadIOMultiBlock.cc: Multi-block versions of input   *
 *                                    and output subroutines for 2D   *
 *                                    High-Temp multi-block       *
 *                                    quadrilateral mesh solution     *
 *                                    classes.                        *
 **********************************************************************/

#include "HighTemp2DQuad.h"

/**********************************************************************
 * HighTemp2D_Quad_Block -- IO Multiple Block External            *
 *                              Subroutines.                          *
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
int Read_Restart_Solution(HighTemp2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          HighTemp2D_Input_Parameters &IP,
                          int &Number_of_Time_Steps,
                          double &Time,
                          CPUTime &CPU_Time) {

  int i, i_new_time_set, nsteps;
  char prefix[256], extension[256], restart_file_name[256];
  //char gas_type[256], particle_type[256];
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
      //restart_file.getline(gas_type,sizeof(gas_type));
      //restart_file.getline(gas_type,sizeof(gas_type));
      if (!i_new_time_set) {

	Number_of_Time_Steps = nsteps;
	IP.Maximum_Number_of_Time_Steps += Number_of_Time_Steps;
	Time = time0;
	CPU_Time.cput = cpu_time0.cput;
	//if ((strcmp(gas_type,IP.Gas_Type) != 0) &&
	//  (strcmp(solid_type,IP.Solid_Type) != 0)) {
	//strcpy(IP.Gas_Type,gas_type);
	//strcpy(IP.Particle_Type,particle_type);
	//IP.Wo.Set_Static_Variables(IP.Gas_Type,IP.FlowType,IP.Electrostatic);
	//IP.Uo.setgas(IP.Gas_Type);
	//IP.Uo = U(IP.Wo);
	//}
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
int Write_Restart_Solution(HighTemp2D_Quad_Block *Soln_ptr,
                           AdaptiveBlock2D_List &Soln_Block_List,
                           HighTemp2D_Input_Parameters &IP,
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
int Output_Tecplot(HighTemp2D_Quad_Block *Soln_ptr,
                   const AdaptiveBlock2D_List &Soln_Block_List,
                   const HighTemp2D_Input_Parameters &IP,
                   int Number_of_Time_Steps,
                   double Time) 
{

	bool writing_intermediate_soln = false;
	double l2_norm = -1.0, l2_norm_rel = -1.0;

	return Output_Tecplot(Soln_ptr,
                   Soln_Block_List,
                   IP,
                   Number_of_Time_Steps,
                   Time,
									 writing_intermediate_soln,
									 l2_norm,
									 l2_norm_rel);
}

int Output_Tecplot(HighTemp2D_Quad_Block *Soln_ptr,
                   const AdaptiveBlock2D_List &Soln_Block_List,
                   const HighTemp2D_Input_Parameters &IP,
                   int Number_of_Time_Steps,
                   double Time,
									 bool writing_intermediate_soln,
									 double l2_norm,
									 double l2_norm_rel)
{

  // Only create output file if the number of blocks used by the 
  // current processor is greater than zero.
  if (Soln_Block_List.Nused() == 0) return 0;
  
  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;    

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' || IP.Output_File_Name[i] == '.') break;
    prefix[i] = IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name)) break;
  }
  prefix[i] = '\0';
  
  // Determine output data file name for this processor.
	if (writing_intermediate_soln) {
		sprintf(extension,"_n1%.4d",Number_of_Time_Steps);
		strcat(prefix,extension);
	}

  strcat(prefix,"_cpu");
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
				 l2_norm,
				 l2_norm_rel,
		     output_file);
      	i_output_title = 0;
    }
  }

  // Close the output data file.
  output_file.close();

  // Writing of output data files complete.
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
int Output_Cells_Tecplot(HighTemp2D_Quad_Block *Soln_ptr,
                         const AdaptiveBlock2D_List &Soln_Block_List,
                         const HighTemp2D_Input_Parameters &IP,
                         int Number_of_Time_Steps,
                         double Time) 
{

	bool writing_intermediate_soln = false;
	double l2_norm = -1.0, l2_norm_rel = -1.0;

	return Output_Cells_Tecplot(Soln_ptr,
                         Soln_Block_List,
                         IP,
                         Number_of_Time_Steps,
                         Time,
												 writing_intermediate_soln,
												 l2_norm,
												 l2_norm_rel);
}

int Output_Cells_Tecplot(HighTemp2D_Quad_Block *Soln_ptr,
                         const AdaptiveBlock2D_List &Soln_Block_List,
                         const HighTemp2D_Input_Parameters &IP,
                         int Number_of_Time_Steps,
                         double Time,
												 bool writing_intermediate_soln,
												 double l2_norm,
												 double l2_norm_rel)
{
  
  // Only create output file if the number of blocks used by the 
  // current processor is greater than zero.
  if (Soln_Block_List.Nused() == 0) return 0;
  
  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;    
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
  
  // Determine output data file name for this processor.
	if (writing_intermediate_soln) {
		sprintf(extension,"_n1%.4d",Number_of_Time_Steps);
		strcat(prefix,extension);
	}

  strcat(prefix,"_cells_cpu");
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
			   Soln_Block_List.Block[nb].gblknum,
			   i_output_title,
				 l2_norm, l2_norm_rel,
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
 * Writes the node values for a 1D array of 2D quadrilateral multi-   *
 * block solution blocks to the specified output data file(s) in a    *
 * format suitable for plotting with TECPLOT.  Returns a non-zero     *
 * value if cannot write any of the TECPLOT solution files.           *
 *                                                                    *
 **********************************************************************/
int Output_Nodes_Tecplot(HighTemp2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List,
                         HighTemp2D_Input_Parameters &IP,
                         const int Number_of_Time_Steps,
                         const double &Time) {
  
  // Only create output file if the number of blocks used by the 
  // current processor is greater than zero.
  if (Soln_Block_List.Nused() == 0) return 0;
  
  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;    
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
 * Routine: Output_Gradients_Tecplot                                  *
 *                                                                    *
 * Writes the gradients of the primitive solution states (as well as  *
 * the slope limiters) for a 1D array of 2D quadrilateral multi-block *
 * solution blocks to the specified output data file(s) in a format   *
 * suitable for plotting with TECPLOT.  Returns a non-zero value if   *
 * cannot write any of the TECPLOT solution files.                    *
 *                                                                    *
 **********************************************************************/
int Output_Gradients_Tecplot(HighTemp2D_Quad_Block *Soln_ptr,
			     AdaptiveBlock2D_List &Soln_Block_List,
			     HighTemp2D_Input_Parameters &IP,
			     const int Number_of_Time_Steps,
			     const double &Time) {
  
  // Only create output file if the number of blocks used by the 
  // current processor is greater than zero.
  if (Soln_Block_List.Nused() == 0) return 0;
  
  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;    
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
  strcat(prefix,"_gradients_cpu");
  
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
      Output_Gradients_Tecplot(Soln_ptr[nb],
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
 * Routine: Output_Quasi3D_Tecplot                                    *
 *                                                                    *
 * Writes the solution values at nodes of the positive cells for a 1D *
 * array of 2D quadrilateral multi-block solution blocks to the       *
 * specified output data file(s) in a quasi-3D format suitable for    *
 * plotting a with TECPLOT.  Returns a non-zero value if cannot write *
 * any of the TECPLOT solution files.                                 *
 *                                                                    *
 **********************************************************************/
int Output_Quasi3D_Tecplot(HighTemp2D_Quad_Block *Soln_ptr,
			   AdaptiveBlock2D_List &Soln_Block_List,
			   HighTemp2D_Input_Parameters &IP,
			   const int Number_of_Time_Steps,
			   const double &Time) {

  // Only create output file if the number of blocks used by the 
  // current processor is greater than zero.
  if (Soln_Block_List.Nused() == 0) return 0;
  
  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;

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
int Output_Mesh_Tecplot(HighTemp2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        HighTemp2D_Input_Parameters &IP,
                        const int Number_of_Time_Steps,
                        const double &Time) {

  // Only create output file if the number of blocks used by the 
  // current processor is greater than zero.
  if (Soln_Block_List.Nused() == 0) return 0;
  
  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;    
  
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
int Output_Mesh_Gnuplot(HighTemp2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        HighTemp2D_Input_Parameters &IP,
                        const int Number_of_Time_Steps,
                        const double &Time) {
  
  // Only create output file if the number of blocks used by the 
  // current processor is greater than zero.
  if (Soln_Block_List.Nused() == 0) return 0;
  
  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;    
  
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

/**********************************************************************
 * Routine: Output_Ringleb_Tecplot                                    *
 *                                                                    *
 * Writes the exact and computed Ringleb's flow solution values at    *
 * the nodes for a 1D array of 2D quadrilateral multi-block solution  *
 * blocks to the specified output data file(s) in a format suitable   *
 * for plotting with TECPLOT.  Returns a non-zero value if cannot     *
 * write any of the TECPLOT solution files.  This routine will also   *
 * compute the error norms of the computed solution (L1, L2, and max  *
 * norms).                                                            *
 *                                                                    *
 **********************************************************************/
int Output_Ringleb_Tecplot(HighTemp2D_Quad_Block *Soln_ptr,
			   AdaptiveBlock2D_List &Soln_Block_List,
			   HighTemp2D_Input_Parameters &IP) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;    
  double l1_norm, l2_norm, max_norm, area;
  double l1_temp, l2_temp, max_temp, area_temp;

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' || IP.Output_File_Name[i] == '.') break;
    prefix[i] = IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_ringleb_cpu");

  // Determine output data file name for this processor.
  sprintf(extension,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;

  // Open the output data file.
  output_file.open(output_file_name_ptr,ios::out);
  if (output_file.bad()) return 1;

  // Write the solution data for each solution block and determine the
  // L1, L2, and Max norms of the error and the number of cells.
  i_output_title = 1;
  l1_norm = ZERO; l2_norm = ZERO; max_norm = ZERO; area = ZERO;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO; area_temp = ZERO;
      Output_Ringleb_Tecplot(Soln_ptr[nb],nb,i_output_title,output_file,
			     l1_temp,l2_temp,max_temp,area_temp);
      l1_norm += l1_temp;
      l2_norm += l2_temp;
      max_norm = max(max_norm,max_temp);
      area += area_temp;
      if (i_output_title) i_output_title = 0;
    }
  }

  // Close the output data file.
  output_file.close();

#ifdef _MPI_VERSION
  l1_norm = CFFC_Summation_MPI(l1_norm);
  l2_norm = CFFC_Summation_MPI(l2_norm);
  max_norm = CFFC_Maximum_MPI(max_norm);
  area = CFFC_Summation_MPI(area);
#endif

  // Calculate the L1-norm and L2-norm for all blocks.
  l1_norm /= area;
  l2_norm = sqrt(l2_norm/area);

  if (CFFC_Primary_MPI_Processor()) {
    cout << endl
	 << endl
	 << " ==================================================================== "
	 << endl
	 << " Error norms for Ringleb's flow:"
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
 * Routine: Output_Viscous_Channel_Tecplot                            *
 *                                                                    *
 * Writes the exact and computed laminar channel flow solution        *
 * (Couette/Poiseuille solutions) at the nodes for a 1D array of 2D   *
 * quadrilateral multi-block solution blocks to the specified output  *
 * data file(s) in a format suitable for plotting with TECPLOT.       *
 * Returns a non-zero value if cannot write any of the TECPLOT        *
 * solution files.  This routine will also compute the error norms of *
 * the u-velocity component of the computed solution (L1, L2, and max *
 * norms).                                                            *
 *                                                                    *
 **********************************************************************/
int Output_Viscous_Channel_Tecplot(HighTemp2D_Quad_Block *Soln_ptr,
				   AdaptiveBlock2D_List &Soln_Block_List,
				   HighTemp2D_Input_Parameters &IP) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;    
  double l1_norm, l2_norm, max_norm;
  double l1_temp, l2_temp, max_temp;
  int numberofcells, numberofcells_temp;

  // Initialize error variables.
  l1_norm = ZERO; l2_norm = ZERO; max_norm = ZERO;
  l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO; 
  numberofcells = 0; numberofcells_temp = 0;

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
  strcat(prefix,"_viscous_channel_cpu");
  
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
      l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO;
      Output_Viscous_Channel_Tecplot(Soln_ptr[nb],nb,i_output_title,output_file,
				     l1_temp,l2_temp,max_temp,IP.Vwall,IP.dp,
				     IP.Box_Width,IP.Box_Height);
      l1_norm += l1_temp;
      l2_norm += sqr(l2_temp);
      max_norm = max(max_norm,max_temp);
      numberofcells_temp = (Soln_ptr[nb].Grid.NCi-2*Soln_ptr[nb].Nghost)*
	                   (Soln_ptr[nb].Grid.NCj-2*Soln_ptr[nb].Nghost);
      numberofcells += numberofcells_temp;
      if (i_output_title) i_output_title = 0;
    }
  }
  
  // Close the output data file.
  output_file.close();
  
#ifdef _MPI_VERSION
  l1_norm = CFFC_Summation_MPI(l1_norm);
  l2_norm = CFFC_Summation_MPI(l2_norm);
  max_norm = CFFC_Maximum_MPI(max_norm);
  numberofcells = CFFC_Summation_MPI(numberofcells);
#endif

  // Calculate the L1-norm and L2-norm for all blocks.
  if (Soln_Block_List.Block[0].used == ADAPTIVEBLOCK2D_USED) {
    l1_norm = l1_norm/double(numberofcells);
    l2_norm = sqrt(l2_norm/double(numberofcells));
  }

  if (CFFC_Primary_MPI_Processor()) {
    cout << endl
	 << endl
	 << " ==================================================================== "
	 << endl
	 << " Error norms for the viscous channel flow:"
	 << endl
	 << "   L1_Norm = " << l1_norm
	 << endl
	 << "   L2_Norm = " << l2_norm
	 << endl
	 << "   Max_Norm = " << max_norm
	 << endl
	 << "   Number of cells = " << numberofcells
	 << endl
	 << " ==================================================================== "
	 << endl;
  }

  // Writing of output data files complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Output_Viscous_Pipe_Tecplot                               *
 *                                                                    *
 * Writes the exact and computed laminar viscous pipe flow solution   *
 * (Hagen-Poiseuille solution) at the nodes for a 1D array of 2D      *
 * quadrilateral multi-block solution blocks to the specified output  *
 * data file(s) in a format suitable for plotting with TECPLOT.       *
 * Returns a non-zero value if cannot write any of the TECPLOT        *
 * solution files.  This routine will also compute the error norms of *
 * the u-velocity component of the computed solution (L1, L2, and max *
 * norms).                                                            *
 *                                                                    *
 **********************************************************************/
int Output_Viscous_Pipe_Tecplot(HighTemp2D_Quad_Block *Soln_ptr,
				AdaptiveBlock2D_List &Soln_Block_List,
				HighTemp2D_Input_Parameters &IP) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;    
  double l1_norm, l2_norm, max_norm;
  double l1_temp, l2_temp, max_temp;
  int numberofcells, numberofcells_temp;

  // Initialize error variables.
  l1_norm = ZERO; l2_norm = ZERO; max_norm = ZERO;
  l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO; 
  numberofcells = 0; numberofcells_temp = 0;

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
  strcat(prefix,"_viscous_pipe_cpu");
  
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
      l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO;
      Output_Viscous_Pipe_Tecplot(Soln_ptr[nb],nb,i_output_title,output_file,
				  l1_temp,l2_temp,max_temp,
				  IP.dp,IP.Pipe_Length,IP.Pipe_Radius);
      l1_norm += l1_temp;
      l2_norm += sqr(l2_temp);
      max_norm = max(max_norm,max_temp);
      numberofcells_temp = (Soln_ptr[nb].Grid.NCi-2*Soln_ptr[nb].Nghost)*
	                   (Soln_ptr[nb].Grid.NCj-2*Soln_ptr[nb].Nghost);
      numberofcells += numberofcells_temp;
      if (i_output_title) i_output_title = 0;
    }
  }
  
  // Close the output data file.
  output_file.close();
  
#ifdef _MPI_VERSION
  l1_norm = CFFC_Summation_MPI(l1_norm);
  l2_norm = CFFC_Summation_MPI(l2_norm);
  max_norm = CFFC_Maximum_MPI(max_norm);
  numberofcells = CFFC_Summation_MPI(numberofcells);
#endif

  // Calculate the L1-norm and L2-norm for all blocks.
  if (Soln_Block_List.Block[0].used == ADAPTIVEBLOCK2D_USED) {
    l1_norm = l1_norm/double(numberofcells);
    l2_norm = sqrt(l2_norm/double(numberofcells));
  }

  if (CFFC_Primary_MPI_Processor()) {
    cout << endl
	 << endl
	 << " ==================================================================== "
	 << endl
	 << " Error norms for the laminar pipe flow:"
	 << endl
	 << "   L1_Norm = " << l1_norm
	 << endl
	 << "   L2_Norm = " << l2_norm
	 << endl
	 << "   Max_Norm = " << max_norm
	 << endl
	 << "   Number of cells = " << numberofcells
	 << endl
	 << " ==================================================================== "
	 << endl;
  }

  // Writing of output data files complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Output_Turbulent_Pipe_Tecplot                             *
 *                                                                    *
 * Writes the experiment and computed solution data for a turbulent   *
 * pipe flow for a 1D array of 2D quadrilateral multi-block solution  *
 * blocks to the specified output data file(s) in a format suitable   *
 * for plotting with TECPLOT.  Returns a non-zero value if cannot     *
 * write any of the TECPLOT solution files.  The experimental data    *
 * was reported by Laufer (19??).                                     *
 *                                                                    *
 **********************************************************************/
int Output_Turbulent_Pipe_Tecplot(HighTemp2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  HighTemp2D_Input_Parameters &IP) {

  int i, i_output_title, i_output_data;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;

  for (int n = 1; n <= 3; n++) {

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
    if (n == 1) strcat(prefix,"_turbulent_pipe_u_cpu");
    else if (n == 2) strcat(prefix,"_turbulent_pipe_k_cpu");
    else if (n == 3) strcat(prefix,"_turbulent_pipe_uv_cpu");
  
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
    if (CFFC_Primary_MPI_Processor()) i_output_data = 1;
    else i_output_data = 0;
    for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
      if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	Output_Turbulent_Pipe_Tecplot(Soln_ptr[nb],nb,i_output_title,i_output_data,output_file,
				      IP.Reynolds_Number,IP.Pipe_Radius,n);
	if (i_output_title) i_output_title = 0;
	if (i_output_data) i_output_data = 0;
      }
    }

    // Close the output data file.
    output_file.close();

  }

  // Writing of output data files complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Output_Flat_Plate_Tecplot                                 *
 *                                                                    *
 * This routine outputs the non-dimensionalized computed flat plate   *
 * solution and the corresponding Blasius solution.                   *
 *                                                                    *
 **********************************************************************/
int Output_Flat_Plate_Tecplot(HighTemp2D_Quad_Block *Soln_ptr,
			      AdaptiveBlock2D_List &Soln_Block_List,
			      HighTemp2D_Input_Parameters &IP) {

  int i, i_output_title_soln, i_output_title_skin;
  char prefix_soln[256], prefix_skin[256], extension_soln[256], extension_skin[256];
  char output_file_name_soln[256], output_file_name_skin[256];
  char *output_file_name_soln_ptr, *output_file_name_skin_ptr;
  ofstream output_file_soln, output_file_skin;
  double l1_norm, l2_norm, max_norm, area;
  double l1_temp, l2_temp, max_temp, area_temp;
  double l1_norm_cf, l2_norm_cf, max_norm_cf, area_cf;
  double l1_cf_temp, l2_cf_temp, max_cf_temp, area_cf_temp;
  int numberofcells, numberofcells_temp;
  int numberofcells_cf, numberofcells_cf_temp;

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' || IP.Output_File_Name[i] == '.') break;
    prefix_soln[i] = IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name)) break;
  }
  prefix_soln[i] = '\0';
  strcat(prefix_soln,"_flatplate_soln_cpu");
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' || IP.Output_File_Name[i] == '.') break;
    prefix_skin[i] = IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name)) break;
  }
  prefix_skin[i] = '\0';
  strcat(prefix_skin,"_flatplate_skin_friction_cpu");
  
  // Determine output data file name for this processor.
  sprintf(extension_soln,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension_soln,".dat");
  strcpy(output_file_name_soln,prefix_soln);
  strcat(output_file_name_soln,extension_soln);
  output_file_name_soln_ptr = output_file_name_soln;
  sprintf(extension_skin,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension_skin,".dat");
  strcpy(output_file_name_skin,prefix_skin);
  strcat(output_file_name_skin,extension_skin);
  output_file_name_skin_ptr = output_file_name_skin;
  
  // Open the output data files.
  output_file_soln.open(output_file_name_soln_ptr,ios::out);
  if (output_file_soln.bad()) return 1;
  output_file_skin.open(output_file_name_skin_ptr,ios::out);
  if (output_file_skin.bad()) return 1;
  
  // Write the solution data for each solution block.
  i_output_title_soln = 1;
  i_output_title_skin = 1;
  l1_norm = ZERO; l2_norm = ZERO; max_norm = ZERO; area = ZERO;
  l1_norm_cf = ZERO; l2_norm_cf = ZERO; max_norm_cf = ZERO; area_cf = ZERO;
  numberofcells = 0; numberofcells_cf = 0;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO; area_temp = 0; numberofcells_temp = 0;
      l1_cf_temp = ZERO; l2_cf_temp = ZERO; max_cf_temp = ZERO; area_cf_temp = 0; numberofcells_cf_temp = 0;
      Output_Flat_Plate_Tecplot(Soln_ptr[nb],nb,
				i_output_title_soln,output_file_soln,
				i_output_title_skin,output_file_skin,
				IP.Wo,IP.Plate_Length,
				l1_temp,l2_temp,max_temp,area_temp,numberofcells_temp,
				l1_cf_temp,l2_cf_temp,max_cf_temp,area_cf_temp,numberofcells_cf_temp);
      l1_norm += l1_temp;
      l2_norm += l2_temp;
      max_norm = max(max_norm,max_temp);
      area += area_temp;
      numberofcells += numberofcells_temp;
      l1_norm_cf += l1_cf_temp;
      l2_norm_cf += l2_cf_temp;
      max_norm_cf = max(max_norm_cf,max_cf_temp);
      area_cf += area_cf_temp;
      numberofcells_cf += numberofcells_cf_temp;
      if (i_output_title_soln) i_output_title_soln = 0;
      if (i_output_title_skin) i_output_title_skin = 0;
    }
  }
  
  // Close the output data file.
  output_file_soln.close();
  output_file_skin.close();

  // Calculate the L1-norm and L2-norm for all blocks.
#ifdef _MPI_VERSION
  l1_norm = CFFC_Summation_MPI(l1_norm);
  l2_norm = CFFC_Summation_MPI(l2_norm);
  max_norm = CFFC_Maximum_MPI(max_norm);
  area = CFFC_Summation_MPI(area);
  numberofcells = CFFC_Summation_MPI(numberofcells);
  l1_norm_cf = CFFC_Summation_MPI(l1_norm_cf);
  l2_norm_cf = CFFC_Summation_MPI(l2_norm_cf);
  max_norm_cf = CFFC_Maximum_MPI(max_norm_cf);
  area_cf = CFFC_Summation_MPI(area_cf);
  numberofcells_cf = CFFC_Summation_MPI(numberofcells_cf);
#endif
  l1_norm /= area;//= l1_norm/numberofcells;
  l2_norm = sqrt(l2_norm/area);//sqrt(l2_norm/numberofcells);
  l1_norm_cf /= area_cf;//= l1_norm_cf/numberofcells_cf;
  l2_norm_cf = sqrt(l2_norm_cf/area_cf);//sqrt(l2_norm_cf/numberofcells_cf);

  if (CFFC_Primary_MPI_Processor()) {
    cout << endl
	 << endl
	 << " ==================================================================== "
	 << endl
	 << " Error norms for the u-velocity component:"
	 << endl
	 << "   L1_Norm = " << l1_norm
	 << endl
	 << "   L2_Norm = " << l2_norm
	 << endl
	 << "   Max_Norm = " << max_norm
	 << endl
	 << endl
	 << " Error norms for the skin friction coefficient:"
	 << endl
	 << "   L1_Norm = " << l1_norm_cf
	 << endl
	 << "   L2_Norm = " << l2_norm_cf
	 << endl
	 << "   Max_Norm = " << max_norm_cf
	 << endl
	 << " ==================================================================== "
	 << endl;
  }

  // Writing of output data files complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Output_Driven_Cavity_Flow_Tecplot                         *
 *                                                                    *
 * This routine outputs a comparison of the computed solution for the *
 * driven cavity flow with the computations done by Ghia et al. (J.   *
 * Comp. Phys. Vol. 48 1982) for a 1D array of 2D quadrilateral       *
 * multi-block solution blocks to the specified output data file(s)   *
 * in a format suitable for plotting with tecplot.                    *
 *                                                                    *
 **********************************************************************/
int Output_Driven_Cavity_Flow_Tecplot(HighTemp2D_Quad_Block *Soln_ptr,
				      AdaptiveBlock2D_List &Soln_Block_List,
				      HighTemp2D_Input_Parameters &IP) {

  int i, i_output_title;
  char prefix_u[256], extension_u[256], output_file_name_u[256];
  char *output_file_name_ptr_u;
  ofstream output_file_u;
  char prefix_v[256], extension_v[256], output_file_name_v[256];
  char *output_file_name_ptr_v;
  ofstream output_file_v;

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' || IP.Output_File_Name[i] == '.') break;
    prefix_u[i] = IP.Output_File_Name[i];
    prefix_v[i] = IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name)) break;
  }
  prefix_u[i] = '\0';
  prefix_v[i] = '\0';
  strcat(prefix_u,"_driven_cavity_flow_u_cpu");
  strcat(prefix_v,"_driven_cavity_flow_v_cpu");

  // Determine output data file name for this processor.
  sprintf(extension_u,"%.6d",Soln_Block_List.ThisCPU);
  sprintf(extension_v,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension_u,".dat");
  strcat(extension_v,".dat");
  strcpy(output_file_name_u,prefix_u);
  strcpy(output_file_name_v,prefix_v);
  strcat(output_file_name_u,extension_u);
  strcat(output_file_name_v,extension_v);
  output_file_name_ptr_u = output_file_name_u;
  output_file_name_ptr_v = output_file_name_v;

  // Open the output data files.
  output_file_u.open(output_file_name_ptr_u,ios::out);
  if (output_file_u.bad()) return 1;
  output_file_v.open(output_file_name_ptr_v,ios::out);
  if (output_file_v.bad()) return 2;

  // Write the solution data for each solution block.
  i_output_title = 1;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Output_Driven_Cavity_Flow_Tecplot(Soln_ptr[nb],nb,
					i_output_title,
					output_file_u,
					output_file_v,
					IP.Re_lid,
					IP.Vwall,
					IP.Box_Width);
      if (i_output_title) i_output_title = 0;
    }
  }
  
  // Close the output data files.
  output_file_u.close();
  output_file_v.close();

  // Writing of output data files complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Output_Backward_Facing_Step_Tecplot                       *
 *                                                                    *
 * This routine outputs a comparison of the computed solution for the *
 * backward facing step flow with the experimental data published by  *
 * Driver and Seepmiller (AIAA J. Vol. 23 No. 2 1985) for a 1D array  *
 * of 2D quadrilateral multi-block solution blocks to the specified   *
 * output data file(s) in a format suitable for plotting with         *
 * tecplot.                                                           *
 *                                                                    *
 **********************************************************************/
int Output_Backward_Facing_Step_Tecplot(HighTemp2D_Quad_Block *Soln_ptr,
					AdaptiveBlock2D_List &Soln_Block_List,
					HighTemp2D_Input_Parameters &IP) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' || IP.Output_File_Name[i] == '.') break;
    prefix[i] = IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_backstep_cpu");

  // Determine output data file name for this processor.
  sprintf(extension,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;

  // Open the output data files.
  output_file.open(output_file_name_ptr,ios::out);
  if (output_file.bad()) return 1;

  // Write the solution data for each solution block.
  i_output_title = 1;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Output_Backward_Facing_Step_Tecplot(Soln_ptr[nb],nb,
					  i_output_title,
					  output_file,
					  IP.Step_Height,
					  IP.Top_Wall_Deflection);
      if (i_output_title) i_output_title = 0;
    }
  }
  
  // Close the output data files.
  output_file.close();

  // Writing of output data files complete.  Return zero value.
  return 0;

}


/**********************************************************************
 * Routine: Output_Forward_Facing_Step_Tecplot                       *
 *                                                                    *
 * This routine outputs a comparison of the computed solution for the *
 * forward facing step flow with CFD data published by  *
 * Glaister (Computational Physics J. No. 77 1988) for a 1D array  *
 * of 2D quadrilateral multi-block solution blocks to the specified   *
 * output data file(s) in a format suitable for plotting with         *
 * tecplot.                                                           *
 *                                                                    *
 **********************************************************************/
int Output_Forward_Facing_Step_Tecplot(HighTemp2D_Quad_Block *Soln_ptr,
					AdaptiveBlock2D_List &Soln_Block_List,
					HighTemp2D_Input_Parameters &IP) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' || IP.Output_File_Name[i] == '.') break;
    prefix[i] = IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_frontstep_cpu");

  // Determine output data file name for this processor.
  sprintf(extension,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;

  // Open the output data files.
  output_file.open(output_file_name_ptr,ios::out);
  if (output_file.bad()) return 1;

  // Write the solution data for each solution block.
  i_output_title = 1;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Output_Forward_Facing_Step_Tecplot(Soln_ptr[nb],nb,
					  i_output_title,
					  output_file,
					  IP.Step_Height,
					  IP.Top_Wall_Deflection);
      if (i_output_title) i_output_title = 0;
    }
  }
  
  // Close the output data files.
  output_file.close();

  // Writing of output data files complete.  Return zero value.
  return 0;

}
