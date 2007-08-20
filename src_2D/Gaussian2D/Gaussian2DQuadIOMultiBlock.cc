/* Gaussian2DQuadMultiBlock.cc:  Multi-Block Versions of IO Subroutines for 2D Gaussian
                                 Multi-Block Quadrilateral Mesh 
                                 Solution Classes. */

/* Include 2D Gaussian quadrilateral mesh solution header file. */

#ifndef _GAUSSIAN2D_QUAD_INCLUDED
#include "Gaussian2DQuad.h"
#endif // _GAUSSIAN2D_QUAD_INCLUDED


/********************************************************
 * Routine: Read_Restart_Solution                       *
 *                                                      *
 * Reads restart solution file(s) and assigns values to *
 * the solution variables of a 1D array of 2D           *
 * quadrilateral multi-block solution blocks.           *
 * Returns a non-zero value if cannot read any of the   *
 * restart solution files.                              *
 *                                                      *
 ********************************************************/
int Read_Restart_Solution(Gaussian2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          Gaussian2D_Input_Parameters &Input_Parameters,
                          int &Number_of_Time_Steps,
                          double &Time,
                          CPUTime &CPU_Time) {

    int i, i_new_time_set, nsteps;
    char prefix[256], extension[256], restart_file_name[256], gas_type[256];
    char *restart_file_name_ptr;
    ifstream restart_file;
    double time0, alpha;
    CPUTime cpu_time0;

    /* Determine prefix of restart file names. */

    i = 0;
    while (1) {
       if (Input_Parameters.Restart_File_Name[i] == ' ' ||
           Input_Parameters.Restart_File_Name[i] == '.') break;
       prefix[i]=Input_Parameters.Restart_File_Name[i];
       i = i + 1;
       if (i > strlen(Input_Parameters.Restart_File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';
    strcat(prefix, "_blk");
    
    /* Read the initial data for each solution block. */

    i_new_time_set = 0;
    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
	  // Restart file name base on global block number.
          sprintf(extension, "%.6d", Soln_Block_List.Block[i].gblknum);
          strcat(extension, ".soln");
          strcpy(restart_file_name, prefix);
          strcat(restart_file_name, extension);
          restart_file_name_ptr = restart_file_name;

          // Open restart file.
          restart_file.open(restart_file_name_ptr, ios::in);
          if (restart_file.bad()) return (1);

          // Read solution block data.
          restart_file.setf(ios::skipws);
          restart_file >> nsteps >> time0 >> cpu_time0;
          restart_file.unsetf(ios::skipws);
          restart_file.getline(gas_type, sizeof(gas_type));
          restart_file.getline(gas_type, sizeof(gas_type));
	  restart_file >> alpha;
          if (!i_new_time_set) {
             Number_of_Time_Steps = nsteps;
             Input_Parameters.Maximum_Number_of_Time_Steps += Number_of_Time_Steps;
	     Time = time0;
             CPU_Time.cput = cpu_time0.cput;
	     Input_Parameters.alpha    = alpha;
	     Input_Parameters.Wo.alpha = alpha;
	     Input_Parameters.Uo.alpha = alpha;
             if (strcmp(gas_type, Input_Parameters.Gas_Type) != 0) {
                strcpy(Input_Parameters.Gas_Type, 
                       gas_type);
                Input_Parameters.Wo.setgas(Input_Parameters.Gas_Type);
                Input_Parameters.Uo.setgas(Input_Parameters.Gas_Type);
                Input_Parameters.Uo = U(Input_Parameters.Wo);
             } /* endif */
             i_new_time_set = 1;
          } /* endif */

          restart_file >> Soln_ptr[i];

          // Close restart file.
          restart_file.close();
       } /* endif */
    }  /* endfor */

    /* Reading of restart files complete.  Return zero value. */

    return(0);

}

/********************************************************
 * Routine: Write_Restart_Solution                      *
 *                                                      *
 * Writes restart solution file(s) for a 1D array of 2D *
 * quadrilateral multi-block solution blocks.           *
 * Returns a non-zero value if cannot write any of the  *
 * restart solution files.                              *
 *                                                      *
 ********************************************************/
int Write_Restart_Solution(Gaussian2D_Quad_Block *Soln_ptr,
                           AdaptiveBlock2D_List &Soln_Block_List,
                           Gaussian2D_Input_Parameters &Input_Parameters,
                           const int Number_of_Time_Steps,
                           const double &Time,
                           const CPUTime &CPU_Time) {

    int i;
    char prefix[256], extension[256], restart_file_name[256];
    char *restart_file_name_ptr;
    ofstream restart_file;

//    //  Save and delete old restart files in compressed archive (just in case)
//    if (CFFC_Primary_MPI_Processor()) {
//      cout << "\n  Creating compressed archive of (and deleting) old restarts.";
//      System::Compress_Restart();
//      cout << "\n  Writing new restart files.";
//      cout.flush();
//    }
//
//    CFFC_Barrier_MPI(); // MPI barrier so that other processors do
//                          // not start over writing restarts
//
//    if (CFFC_Primary_MPI_Processor()) {
//      System::Set_Restart_Flag();  //Set flag to indicate a restart is being saved
//    }

    /* Determine prefix of restart file names. */

    i = 0;
    while (1) {
       if (Input_Parameters.Restart_File_Name[i] == ' ' ||
           Input_Parameters.Restart_File_Name[i] == '.') break;
       prefix[i]=Input_Parameters.Restart_File_Name[i];
       i = i + 1;
       if (i > strlen(Input_Parameters.Restart_File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';
    strcat(prefix, "_blk");

    /* Write the solution data for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
	  // Restart file name base on global block number.
          sprintf(extension, "%.6d", Soln_Block_List.Block[i].gblknum);
          strcat(extension, ".soln");
          strcpy(restart_file_name, prefix);
          strcat(restart_file_name, extension);
          restart_file_name_ptr = restart_file_name;

          // Open restart file.
          restart_file.open(restart_file_name_ptr, ios::out);
          if (restart_file.bad()) return (1);

          // Write solution block data.
          restart_file.setf(ios::scientific);
          restart_file << setprecision(14) << Number_of_Time_Steps 
                       << " " << Time << " " << CPU_Time << "\n";
          restart_file.unsetf(ios::scientific);
          restart_file << Input_Parameters.Gas_Type << "\n";
	  restart_file << Input_Parameters.Wo.alpha << "\n";
          restart_file << setprecision(14) << Soln_ptr[i];

          // Close restart file.
          restart_file.close();
       } /* endif */
    }  /* endfor */

    CFFC_Barrier_MPI(); // MPI barrier so that all processors have completed
                          // writing.

//    if (CFFC_Primary_MPI_Processor()) {
//      System::Remove_Restart_Flag();  //Remove flag to indicate the restart is finished
//    }
//
    /* Writing of restart files complete.  Return zero value. */

    return(0);

}

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the solution values at the nodes for a 1D     *
 * array of 2D quadrilateral multi-block solution       *
 * blocks to the specified output data file(s) in a     *
 * format suitable for plotting with TECPLOT.           *
 * Returns a non-zero value if cannot write any of the  *
 * TECPLOT solution files.                              *
 *                                                      *
 ********************************************************/
int Output_Tecplot(Gaussian2D_Quad_Block *Soln_ptr,
                   AdaptiveBlock2D_List &Soln_Block_List,
                   Gaussian2D_Input_Parameters &Input_Parameters,
                   const int Number_of_Time_Steps,
                   const double &Time) {

    int i, i_output_title;
    char prefix[256], extension[256], output_file_name[256];
    char *output_file_name_ptr;
    ofstream output_file;    

    /* Determine prefix of output data file names. */

    i = 0;
    while (1) {
       if (Input_Parameters.Output_File_Name[i] == ' ' ||
           Input_Parameters.Output_File_Name[i] == '.') break;
       prefix[i]=Input_Parameters.Output_File_Name[i];
       i = i + 1;
       if (i > strlen(Input_Parameters.Output_File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';
    strcat(prefix, "_cpu");

    /* Determine output data file name for this processor. */

    sprintf(extension, "%.6d", Soln_Block_List.ThisCPU);
    strcat(extension, ".dat");
    strcpy(output_file_name, prefix);
    strcat(output_file_name, extension);
    output_file_name_ptr = output_file_name;

    /* Open the output data file. */

    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) return (1);

    /* Write the solution data for each solution block. */

    i_output_title = 1;
    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          Output_Tecplot(Soln_ptr[i],
                         Input_Parameters,
                         Number_of_Time_Steps, 
                         Time,
                         Soln_Block_List.Block[i].gblknum,
                         i_output_title,
                         output_file);
	  if (i_output_title) i_output_title = 0;
       } /* endif */
    }  /* endfor */

    /* Close the output data file. */

    output_file.close();

    /* Writing of output data files complete.  Return zero value. */

    return(0);

}

/********************************************************
 * Routine: Output_Cells_Tecplot                        *
 *                                                      *
 * Writes the cell centred solution values for a 1D     *
 * array of 2D quadrilateral multi-block solution       *
 * blocks to the specified output data file(s) in a     *
 * format suitable for plotting with TECPLOT.           *
 * Returns a non-zero value if cannot write any of the  *
 * TECPLOT solution files.                              *
 *                                                      *
 ********************************************************/
int Output_Cells_Tecplot(Gaussian2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List,
                         Gaussian2D_Input_Parameters &Input_Parameters,
                         const int Number_of_Time_Steps,
                         const double &Time) {

    int i, i_output_title;
    char prefix[256], extension[256], output_file_name[256];
    char *output_file_name_ptr;
    ofstream output_file;    

    /* Determine prefix of output data file names. */

    i = 0;
    while (1) {
       if (Input_Parameters.Output_File_Name[i] == ' ' ||
           Input_Parameters.Output_File_Name[i] == '.') break;
       prefix[i]=Input_Parameters.Output_File_Name[i];
       i = i + 1;
       if (i > strlen(Input_Parameters.Output_File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';
    strcat(prefix, "_cells_cpu");

    /* Determine output data file name for this processor. */

    sprintf(extension, "%.6d", Soln_Block_List.ThisCPU);
    strcat(extension, ".dat");
    strcpy(output_file_name, prefix);
    strcat(output_file_name, extension);
    output_file_name_ptr = output_file_name;

    /* Open the output data file. */

    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) return (1);

    /* Write the solution data for each solution block. */

    i_output_title = 1;
    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          Output_Cells_Tecplot(Soln_ptr[i], 
			       Input_Parameters,
                               Number_of_Time_Steps, 
                               Time,
                               Soln_Block_List.Block[i].gblknum,
                               i_output_title,
                               output_file);
	  if (i_output_title) i_output_title = 0;
       } /* endif */
    }  /* endfor */

    /* Close the output data file. */

    output_file.close();

    /* Writing of output data files complete.  Return zero value. */

    return(0);

}

/********************************************************
 * Routine: Output_Mesh_Tecplot                         *
 *                                                      *
 * Writes the nodes of the mesh for a 1D array of 2D    *
 * quadrilateral multi-block solution blocks to the     *
 * specified output data file(s) in a format suitable   *
 * for plotting with TECPLOT. Returns a non-zero value  *
 * if cannot write any of the TECPLOT solution files.   *
 *                                                      *
 ********************************************************/
int Output_Mesh_Tecplot(Gaussian2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        Gaussian2D_Input_Parameters &Input_Parameters,
                        const int Number_of_Time_Steps,
                        const double &Time) {

    int i, i_output_title;
    char prefix[256], extension[256], output_file_name[256];
    char *output_file_name_ptr;
    ofstream output_file;    

    /* Determine prefix of output data file names. */

    i = 0;
    while (1) {
       if (Input_Parameters.Output_File_Name[i] == ' ' ||
           Input_Parameters.Output_File_Name[i] == '.') break;
       prefix[i]=Input_Parameters.Output_File_Name[i];
       i = i + 1;
       if (i > strlen(Input_Parameters.Output_File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';
    strcat(prefix, "_mesh_cpu");

    /* Determine output data file name for this processor. */

    sprintf(extension, "%.6d", Soln_Block_List.ThisCPU);
    strcat(extension, ".dat");
    strcpy(output_file_name, prefix);
    strcat(output_file_name, extension);
    output_file_name_ptr = output_file_name;

    /* Open the output data file. */

    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) return (1);

    /* Write the solution data for each solution block. */

    i_output_title = 1;
    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          Output_Tecplot(Soln_ptr[i].Grid, 
                         Soln_Block_List.Block[i].gblknum,
                         i_output_title,
                         output_file);
//           Output_Nodes_Tecplot(Soln_ptr[i].Grid, 
//                                Soln_Block_List.Block[i].gblknum,
//                                i_output_title,
//                                output_file);
	  if (i_output_title) i_output_title = 0;
       } /* endif */
    }  /* endfor */

    /* Close the output data file. */

    output_file.close();

    /* Writing of output data files complete.  Return zero value. */

    return(0);

}

/********************************************************
 * Routine: Output_Mesh_Gnuplot                         *
 *                                                      *
 * Writes the nodes of the mesh for a 1D array of 2D    *
 * quadrilateral multi-block solution blocks to the     *
 * specified output data file(s) in a format suitable   *
 * for plotting with GNUPLOT. Returns a non-zero value  *
 * if cannot write any of the GNUPLOT solution files.   *
 *                                                      *
 ********************************************************/
int Output_Mesh_Gnuplot(Gaussian2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        Gaussian2D_Input_Parameters &Input_Parameters,
                        const int Number_of_Time_Steps,
                        const double &Time) {

    int i, i_output_title;
    char prefix[256], extension[256], output_file_name[256];
    char *output_file_name_ptr;
    ofstream output_file;    

    /* Determine prefix of output data file names. */

    i = 0;
    while (1) {
       if (Input_Parameters.Output_File_Name[i] == ' ' ||
           Input_Parameters.Output_File_Name[i] == '.') break;
       prefix[i]=Input_Parameters.Output_File_Name[i];
       i = i + 1;
       if (i > strlen(Input_Parameters.Output_File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';
    strcat(prefix, "_mesh_cpu");

    /* Determine output data file name for this processor. */

    sprintf(extension, "%.6d", Soln_Block_List.ThisCPU);
    strcat(extension, ".dat");
    strcpy(output_file_name, prefix);
    strcat(output_file_name, extension);
    output_file_name_ptr = output_file_name;

    /* Open the output data file. */

    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) return (1);

    /* Write the solution data for each solution block. */

    i_output_title = 1;
    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          Output_Gnuplot(Soln_ptr[i].Grid, 
                         Soln_Block_List.Block[i].gblknum,
                         i_output_title,
                         output_file);
	  if (i_output_title) i_output_title = 0;
       } /* endif */
    }  /* endfor */

    /* Close the output data file. */

    output_file.close();

    /* Writing of output data files complete.  Return zero value. */

    return(0);

}

/**********************************************************************
 * Routine: Output_Flat_Plate                                         *
 *                                                                    *
 * This routine outputs the non-dimensionalized computed flat plate   *
 * solution and the corresponding Blasius solution.                   *
 *                                                                    *
 **********************************************************************/
int Output_Flat_Plate(Gaussian2D_Quad_Block *Soln_ptr,
		      AdaptiveBlock2D_List &Soln_Block_List,
		      Gaussian2D_Input_Parameters &IP,
		      const int Number_of_Time_Steps,
		      const double &Time) {

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
  strcat(prefix,"_flatplate_soln_cpu");

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
  for (i = 0; i <= Soln_Block_List.Nblk-1; i++) {
    if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
      Output_Flat_Plate(Soln_ptr[i],
			Number_of_Time_Steps,
			Time,
			Soln_Block_List.Block[i].gblknum,
			i_output_title,
			output_file,
			IP.Wo);
      if (i_output_title) i_output_title = 0;
    }
  }
  
  // Close the output data file.
  output_file.close();

  // Writing of output data files complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Output_Cylinder_Free_Molecular                            *
 *                                                                    *
 * This routine outputs the boudary pressure and shear profile for a  *
 * cylinder in free molecular flow.                                   *
 *                                                                    *
 **********************************************************************/
int Output_Cylinder_Free_Molecular(Gaussian2D_Quad_Block *Soln_ptr,
				   AdaptiveBlock2D_List &Soln_Block_List,
				   Gaussian2D_Input_Parameters &IP,
				   const int Number_of_Time_Steps,
				   const double &Time) {

  int i, j, i_output_title, count, number_sent(0), offset(0), start(0),
      number_of_boundary_cells_this_cpu(0), total_number_of_boundary_cells(0),
      max_number_of_boundary_cells_per_cpu(0), sorted(0);
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;
  double *buffer1, *buffer2, temp_angle1, temp_angle2;
  Vector2D *nodes, temp_vector;

  //Determine number of cells touching cylinder on this cpu, all cpus and max on any cpu

  for (i = 0; i <= Soln_Block_List.Nblk-1; i++) {
    if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
      number_of_boundary_cells_this_cpu += Count_south_adiabatic_wall_cells(Soln_ptr[i]);
    }
  }
  
  total_number_of_boundary_cells = CFFC_Summation_MPI(number_of_boundary_cells_this_cpu);
  max_number_of_boundary_cells_per_cpu = CFFC_Maximum_MPI(number_of_boundary_cells_this_cpu);

  //Create buffers and array of nodes

  buffer1 = new double[max_number_of_boundary_cells_per_cpu];
  buffer2 = new double[max_number_of_boundary_cells_per_cpu];
  nodes = new Vector2D[total_number_of_boundary_cells];

  //combine list of nodes from each processor

  for(i = 0; i < CFFC_MPI::Number_of_Processors; i++) {

    if(CFFC_MPI::This_Processor_Number == i) {
      //load send buffers
      start = 0;
      number_sent = number_of_boundary_cells_this_cpu;
      for(j = 0; j <= Soln_Block_List.Nblk-1; j++) {
	if (Soln_Block_List.Block[j].used == ADAPTIVEBLOCK2D_USED) {
	  Append_nodes_to_send_buffer(Soln_ptr[j], buffer1, buffer2, start);
	}
      }
    } //send buffer loaded


#ifdef _MPI_VERSION
    MPI::COMM_WORLD.Bcast(&number_sent, 1, MPI::INT, i);
    MPI::COMM_WORLD.Bcast(buffer1, number_sent, MPI::DOUBLE, i);
    MPI::COMM_WORLD.Bcast(buffer2, number_sent, MPI::DOUBLE, i);
#endif

    //Unload sent buffer

    for(j = 0; j < number_sent; j++) {
      assert(j+offset < total_number_of_boundary_cells);  //check that we're still in the array!
      nodes[j + offset ] = Vector2D(buffer1[j], buffer2[j]);
    }
    offset += number_sent;
  }

  //Bubble Sort nodes array

  while(!sorted) {

    for(i = 0 ; i < total_number_of_boundary_cells; i++) {

      if(i == total_number_of_boundary_cells-1){
	sorted = 1;
	break;
      }

      temp_angle1 = atan2(nodes[i].y, nodes[i].x);
      temp_angle2 = atan2(nodes[i+1].y, nodes[i+1].x);
      if(temp_angle1 < 0.0) temp_angle1 += 2.0*PI;
      if(temp_angle2 < 0.0) temp_angle2 += 2.0*PI;

      if(temp_angle1 > temp_angle2){ 
	temp_vector = nodes[i+1];
	nodes[i+1] = nodes[i];
	nodes[i] = temp_vector;
	break;
      }

    }

  }

  // Determine prefix of output data file names.
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' || IP.Output_File_Name[i] == '.') break;
    prefix[i] = IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_free_molecular_soln_cpu");

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
  for (i = 0; i <= Soln_Block_List.Nblk-1; i++) {
    if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
      Output_Cylinder_Free_Molecular(Soln_ptr[i],
				     Number_of_Time_Steps,
				     Time,
				     Soln_Block_List.Block[i].gblknum,
				     i_output_title,
				     output_file,
				     IP.Wo,
				     nodes,
				     total_number_of_boundary_cells);

      if (i_output_title) i_output_title = 0;
    }
  }
  
  // Close the output data file.
  output_file.close();

  //delete dynamic buffers

  delete []buffer1;
  delete []buffer2;
  delete []nodes;

  // Writing of output data files complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Output_Drag                                               *
 *                                                                    *
 * This routine outputs the Drag on a body                            *
 *                                                                    *
 **********************************************************************/
void Output_Drag(Gaussian2D_Quad_Block *Soln_ptr,
		 AdaptiveBlock2D_List &Soln_Block_List,
		 double &drag, double &lift) {
  int i;

  for(i = 0; i < Soln_Block_List.Nblk-1;i++) {
    if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
      Output_Drag(Soln_ptr[i],
		  drag, lift);
    }
  }

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
int Output_Gradients_Tecplot(Gaussian2D_Quad_Block *Soln_ptr,
			     AdaptiveBlock2D_List &Soln_Block_List,
			     Gaussian2D_Input_Parameters &IP,
			     const int Number_of_Time_Steps,
			     const double &Time) {

  cout << "This has not been implemented." << endl;

  //this is here only for compalibility with embedded boundaries
  return 0;
}
