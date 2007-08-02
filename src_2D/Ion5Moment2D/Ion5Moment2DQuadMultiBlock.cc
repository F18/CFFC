/* Ion5Moment2DQuadMultiBlock.cc:  Multi-Block Versions of Subroutines for
                                   2D 5-Moment Ion Transport Model
                                   Multi-Block Quadrilateral Mesh 
                                   Solution Classes. */

/* Include 2D 5-moment ion transport model quadrilateral mesh solution header file. */

#ifndef _ION5MOMENT2D_QUAD_INCLUDED
#include "Ion5Moment2DQuad.h"
#endif // _ION5MOMENT2D_QUAD_INCLUDED

/**************************************************************************
 * Ion5Moment2D_Quad_Block -- Multiple Block External Subroutines.        *
 **************************************************************************/

/********************************************************
 * Routine: Allocate                                    *
 *                                                      *
 * Allocate memory for 1D array of 2D quadrilateral     *
 * multi-block solution blocks.                         *
 *                                                      *
 ********************************************************/
Ion5Moment2D_Quad_Block* Allocate(Ion5Moment2D_Quad_Block *Soln_ptr,
                                  Ion5Moment2D_Input_Parameters &Input_Parameters) {

    /* Allocate memory. */

    Soln_ptr = new Ion5Moment2D_Quad_Block[Input_Parameters.Number_of_Blocks_Per_Processor];

    /* Return memory location. */

    return(Soln_ptr);

}

/********************************************************
 * Routine: Deallocate                                  *
 *                                                      *
 * Deallocate memory for 1D array of 2D quadrilateral   *
 * multi-block solution blocks.                         *
 *                                                      *
 ********************************************************/
Ion5Moment2D_Quad_Block* Deallocate(Ion5Moment2D_Quad_Block *Soln_ptr,
                                    Ion5Moment2D_Input_Parameters &Input_Parameters) {

    int i;
 
    /* Deallocate memory. */

    for ( i = 0 ; i <= Input_Parameters.Number_of_Blocks_Per_Processor-1 ; ++i ) {
       if (Soln_ptr[i].W != NULL) Soln_ptr[i].deallocate();
    }  /* endfor */
    delete []Soln_ptr;
    Soln_ptr = NULL;

    /* Return memory location. */

    return(Soln_ptr);

}

/********************************************************
 * Routine: ICs                                         *
 *                                                      *
 * Assigns initial conditions and data to the           *
 * solution variables of a 1D array of 2D quadrilateral *
 * multi-block solution blocks.                         *
 *                                                      *
 ********************************************************/
void ICs(Ion5Moment2D_Quad_Block *Soln_ptr,
         AdaptiveBlock2D_List &Soln_Block_List,
         Ion5Moment2D_Input_Parameters &Input_Parameters) {

    int i;
    Ion5Moment2D_pState Wo[5];

    /* Assign the ion and neutral gas constants for the 
       ions and neutral gas of interest. */

    Input_Parameters.Wo.setion(Input_Parameters.Ion_Type);
    Input_Parameters.Uo.setion(Input_Parameters.Ion_Type);
    Input_Parameters.Uo = U(Input_Parameters.Wo);

    Input_Parameters.Wno.setgas(Input_Parameters.Neutral_Gas_Type); // for neutral gas solution states
    Input_Parameters.Wo.setgas(Input_Parameters.Neutral_Gas_Type); // for ion solution states
    Input_Parameters.Uo.setgas(Input_Parameters.Neutral_Gas_Type);

    /* Define various reference flow states. */

    Wo[0] = Input_Parameters.Wo;
    Wo[1] = Input_Parameters.W1;
    Wo[2] = Input_Parameters.W2;

    /* Assign initial data for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
	  // Set flow geometry indicator (planar/axisymmetric) flow block.
 	  Soln_ptr[i].Axisymmetric = Input_Parameters.Axisymmetric;

          // Set initial data for the ions.
          ICs(Soln_ptr[i], 
              Input_Parameters.i_ICs,
              Wo);

          // Set initial data for the neutral gas.
          ICs_Neutral_Gas(Soln_ptr[i], 
                          Input_Parameters.i_ICs,
                          Input_Parameters.Wno);

          // Set initial data for the electric field.
          ICs_Electric_Field(Soln_ptr[i], 
                             Input_Parameters.i_Electric_Field,
                             Input_Parameters.Electric_Field_Strength,
                             Input_Parameters.Electric_Field_Angle);
       } /* endif */
    }  /* endfor */

}

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
int Read_Restart_Solution(Ion5Moment2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          Ion5Moment2D_Input_Parameters &Input_Parameters,
                          int &Number_of_Time_Steps,
                          double &Time,
                          CPUTime &CPU_Time) {

    int i, i_new_time_set, nsteps;
    char prefix[256], extension[256], restart_file_name[256], 
         ion_type[256], neutral_gas_type[256];
    char *restart_file_name_ptr;
    ifstream restart_file;
    double time0;
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
          restart_file.getline(ion_type, sizeof(ion_type));
          restart_file.getline(ion_type, sizeof(ion_type));
          restart_file.getline(neutral_gas_type, sizeof(neutral_gas_type));
          if (!i_new_time_set) {
             Number_of_Time_Steps = nsteps;
             Input_Parameters.Maximum_Number_of_Time_Steps += Number_of_Time_Steps;
	     Time = time0;
             CPU_Time.cput = cpu_time0.cput;
             if (strcmp(ion_type, Input_Parameters.Ion_Type) != 0) {
                strcpy(Input_Parameters.Ion_Type, 
                       ion_type);
                Input_Parameters.Wo.setion(Input_Parameters.Ion_Type);
                Input_Parameters.Uo.setion(Input_Parameters.Ion_Type);
                Input_Parameters.Uo = U(Input_Parameters.Wo);
             } /* endif */
             if (strcmp(neutral_gas_type, Input_Parameters.Neutral_Gas_Type) != 0) {
                strcpy(Input_Parameters.Neutral_Gas_Type, 
                       neutral_gas_type);
                Input_Parameters.Wno.setgas(Input_Parameters.Neutral_Gas_Type);
                Input_Parameters.Wo.setgas(Input_Parameters.Neutral_Gas_Type);
                Input_Parameters.Uo.setgas(Input_Parameters.Neutral_Gas_Type);
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
int Write_Restart_Solution(Ion5Moment2D_Quad_Block *Soln_ptr,
                           AdaptiveBlock2D_List &Soln_Block_List,
                           Ion5Moment2D_Input_Parameters &Input_Parameters,
                           const int Number_of_Time_Steps,
                           const double &Time,
                           const CPUTime &CPU_Time) {

    int i;
    char prefix[256], extension[256], restart_file_name[256];
    char *restart_file_name_ptr;
    ofstream restart_file;

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
          restart_file << Input_Parameters.Ion_Type << "\n";
          restart_file << Input_Parameters.Neutral_Gas_Type << "\n";
          restart_file << setprecision(14) << Soln_ptr[i];

          // Close restart file.
          restart_file.close();
       } /* endif */
    }  /* endfor */

    /* Writing of restart files complete.  Return zero value. */

    return(0);

}

/********************************************************
 * Routine: Read_Neutral_Gas_Solution                   *
 *                                                      *
 * Reads the neutral gas solution file and assigns      *
 * values to the neutral gas flow solution variables    *
 * for a 1D array of 2D quadrilateral multi-block       *
 * solution blocks.  The netural flow solution file is  *
 * expected to be in ASCII TECPLOT point data format.   *
 * Returns a non-zero value if cannot read the input    *
 * solution file.                                       *
 *                                                      *
 ********************************************************/
int Read_Neutral_Gas_Solution(Ion5Moment2D_Quad_Block *Soln_ptr,
                              AdaptiveBlock2D_List &Soln_Block_List,
                              Ion5Moment2D_Input_Parameters &Input_Parameters) {

    int i, line, blk_num, ii, jj, kk;
    double x, y, z;
    char temp[256];
    char *neutral_gas_solution_file_name_ptr, *line_ptr;
    ifstream neutral_gas_solution_file;

    /* Read the neutral gas solution data for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
	  // Open neutral gas input solution file.
          neutral_gas_solution_file_name_ptr = Input_Parameters.Neutral_Gas_Solution_File_Name;
          neutral_gas_solution_file.open(neutral_gas_solution_file_name_ptr, ios::in);
          if (neutral_gas_solution_file.bad()) return (1);

	  // Read the TECPLOT title and variable information.
          for (line = 1; line <= 11; ++line) {
             neutral_gas_solution_file.getline(temp, sizeof(temp));
          } /* endfor */

          // Skip over unread solution blocks in TECPLOT file.
	  for (blk_num = 0; blk_num < Soln_Block_List.Block[i].gblknum; ++blk_num) {
  	     // Read the TECPLOT zone information and get zone size (ii, jj, kk).
             neutral_gas_solution_file.getline(temp, sizeof(temp));
             neutral_gas_solution_file.getline(temp, sizeof(temp));
             line_ptr = index(temp, '=')+1;
             sscanf(line_ptr, "%d", &ii);
             line_ptr = index(line_ptr, '=')+1;
             sscanf(line_ptr, "%d", &jj);
             line_ptr = index(line_ptr, '=')+1;
             sscanf(line_ptr, "%d", &kk);
             neutral_gas_solution_file.getline(temp, sizeof(temp));

             // Skip over data points.
             for (line = 1; line <= ii*jj; ++line) {
                neutral_gas_solution_file.getline(temp, sizeof(temp));
             } /* endfor */
          } /* endfor */

  	  // Read the TECPLOT zone information for zone of interest.
          for (line = 1; line <= 3; ++line) {
             neutral_gas_solution_file.getline(temp, sizeof(temp));
          } /* endfor */

          // Read neutral gas solution data for zone of interest.
          Read_Neutral_Gas_Solution_Block(Soln_ptr[i], 
                                          neutral_gas_solution_file);

          // Close neutral gas input solution file.
          neutral_gas_solution_file.close();
       } /* endif */
    }  /* endfor */

    /* Reading of neutral gas input solution file complete.  
       Return zero value. */

    return(0);

}

/********************************************************
 * Routine: Read_Electric_Field_Solution                *
 *                                                      *
 * Reads the electric field solution file and assigns   *
 * values to the electric field solution variables      *
 * for a 1D array of 2D quadrilateral multi-block       *
 * solution blocks.  The electric field solution file   *
 * is expected to be in ASCII TECPLOT point data format.*
 * Returns a non-zero value if cannot read the input    *
 * solution file.                                       *
 *                                                      *
 ********************************************************/
int Read_Electric_Field_Solution(Ion5Moment2D_Quad_Block *Soln_ptr,
                                 AdaptiveBlock2D_List &Soln_Block_List,
                                 Ion5Moment2D_Input_Parameters &Input_Parameters) {

    int i, line, blk_num, ii, jj, kk;
    double x, y, z;
    char temp[256];
    char *electric_field_solution_file_name_ptr, *line_ptr;
    ifstream electric_field_solution_file;

    /* Read the electric field solution data for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
	  // Open electric field input solution file.
          electric_field_solution_file_name_ptr = Input_Parameters.Electric_Field_Solution_File_Name;
          electric_field_solution_file.open(electric_field_solution_file_name_ptr, ios::in);
          if (electric_field_solution_file.bad()) return (1);

	  // Read the TECPLOT title and variable information.
          for (line = 1; line <= 8; ++line) {
             electric_field_solution_file.getline(temp, sizeof(temp));
          } /* endfor */

          // Skip over unread solution blocks in TECPLOT file.
	  for (blk_num = 0; blk_num < Soln_Block_List.Block[i].gblknum; ++blk_num) {
  	     // Read the TECPLOT zone information and get zone size (ii, jj, kk).
             electric_field_solution_file.getline(temp, sizeof(temp));
             electric_field_solution_file.getline(temp, sizeof(temp));
             line_ptr = index(temp, '=')+1;
             sscanf(line_ptr, "%d", &ii);
             line_ptr = index(line_ptr, '=')+1;
             sscanf(line_ptr, "%d", &jj);
             line_ptr = index(line_ptr, '=')+1;
             sscanf(line_ptr, "%d", &kk);
             electric_field_solution_file.getline(temp, sizeof(temp));

             // Skip over data points.
             for (line = 1; line <= ii*jj; ++line) {
                electric_field_solution_file.getline(temp, sizeof(temp));
             } /* endfor */
          } /* endfor */

  	  // Read the TECPLOT zone information for zone of interest.
          for (line = 1; line <= 3; ++line) {
             electric_field_solution_file.getline(temp, sizeof(temp));
          } /* endfor */

          // Read electric field solution data for zone of interest.
          Read_Electric_Field_Solution_Block(Soln_ptr[i], 
                                             electric_field_solution_file,
                                             Input_Parameters.Add_Initial_and_Solution_File_Electric_Fields);

          // Close electric field input solution file.
          electric_field_solution_file.close();
       } /* endif */
    }  /* endfor */

    /* Reading of electric field input solution file complete.  
       Return zero value. */

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
int Output_Tecplot(Ion5Moment2D_Quad_Block *Soln_ptr,
                   AdaptiveBlock2D_List &Soln_Block_List,
                   Ion5Moment2D_Input_Parameters &Input_Parameters,
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
int Output_Cells_Tecplot(Ion5Moment2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List,
                         Ion5Moment2D_Input_Parameters &Input_Parameters,
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
int Output_Mesh_Tecplot(Ion5Moment2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        Ion5Moment2D_Input_Parameters &Input_Parameters,
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
int Output_Mesh_Gnuplot(Ion5Moment2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        Ion5Moment2D_Input_Parameters &Input_Parameters,
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

/********************************************************
 * Routine: BCs                                         *
 *                                                      *
 * Apply boundary conditions at boundaries of a 1D      *
 * array of 2D quadrilateral multi-block solution       *
 * blocks.                                              *
 *                                                      *
 ********************************************************/
void BCs(Ion5Moment2D_Quad_Block *Soln_ptr,
         AdaptiveBlock2D_List &Soln_Block_List,
	 Ion5Moment2D_Input_Parameters &Input_Parameters) {

    int i;

    /* Prescribe ion flow boundary data for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          BCs(Soln_ptr[i]);
       } /* endif */
    }  /* endfor */

}

/********************************************************
 * Routine: BCs_Neutral_Gas                             *
 *                                                      *
 * Apply boundary conditions for the neutral gas flow   *
 * at boundaries of a 1D array of 2D quadrilateral      *
 * multi-block solution blocks.                         *
 *                                                      *
 ********************************************************/
void BCs_Neutral_Gas(Ion5Moment2D_Quad_Block *Soln_ptr,
                     AdaptiveBlock2D_List &Soln_Block_List) {

    int i;

    /* Prescribe neutral gas flow boundary data for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          BCs_Neutral_Gas(Soln_ptr[i]);
       } /* endif */
    }  /* endfor */

}

/********************************************************
 * Routine: BCs_Electric_Field                          *
 *                                                      *
 * Apply boundary conditions for the electric field     *
 * solution at boundaries of a 1D array of 2D           *
 * quadrilateral multi-block solution blocks.           *
 *                                                      *
 ********************************************************/
void BCs_Electric_Field(Ion5Moment2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List) {

    int i;

    /* Prescribe electric field boundary data for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          BCs_Electric_Field(Soln_ptr[i]);
       } /* endif */
    }  /* endfor */

}

/********************************************************
 * Routine: CFL                                         *
 *                                                      *
 * Determines the allowable global and local time steps *
 * (for explicit Euler time stepping scheme) for a 1D   *
 * array of 2D quadrilateral multi-block solution       *
 * blocks according to the Courant-Friedrichs-Lewy      *
 * condition for the wave motion terms and a            *
 * semi-impirical criteria for the ion-neutral          *
 * collisions source terms.                             *
 *                                                      *
 ********************************************************/
double CFL(Ion5Moment2D_Quad_Block *Soln_ptr,
           AdaptiveBlock2D_List &Soln_Block_List,
           Ion5Moment2D_Input_Parameters &Input_Parameters) {

    int i;
    double dtMin;

    dtMin = MILLION;

    /* Determine the allowable time step for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          dtMin = min(dtMin, CFL(Soln_ptr[i], Input_Parameters));
       } /* endif */
    }  /* endfor */

    /* Return the global time step. */

    return (dtMin);

}

/********************************************************
 * Routine: Set_Global_TimeStep                         *
 *                                                      *
 * Assigns global time step to a 1D array of 2D         *
 * quadrilateral multi-block solution blocks for        *
 * time-accurate calculations.                          *
 *                                                      *
 ********************************************************/
void Set_Global_TimeStep(Ion5Moment2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List,
                         const double &Dt_min) {

    int i;

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          Set_Global_TimeStep(Soln_ptr[i], Dt_min);
       } /* endif */
    }  /* endfor */

}

/********************************************************
 * Routine: L1_Norm_Residual                            *
 *                                                      *
 * Determines the L1-norm of the solution residual for  *
 * a 1D array of 2D quadrilateral multi-block solution  *
 * blocks.  Useful for monitoring convergence of the    *
 * solution for steady state problems.                  *
 *                                                      *
 ********************************************************/
double L1_Norm_Residual(Ion5Moment2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
			Ion5Moment2D_Input_Parameters &Input_Parameters) {

    int i;
    double l1_norm;

    l1_norm = ZERO;

    /* Calculate the L1-norm.
       Sum the L1-norm for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          l1_norm += L1_Norm_Residual(Soln_ptr[i]);
       } /* endif */
    }  /* endfor */

    /* Return the L1-norm. */

    return (l1_norm);

}

/********************************************************
 * Routine: L2_Norm_Residual                            *
 *                                                      *
 * Determines the L2-norm of the solution residual for  *
 * a 1D array of 2D quadrilateral multi-block solution  *
 * blocks.  Useful for monitoring convergence of the    *
 * solution for steady state problems.                  *
 *                                                      *
 ********************************************************/
double L2_Norm_Residual(Ion5Moment2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
			Ion5Moment2D_Input_Parameters &Input_Parameters) {

    int i;
    double l2_norm;

    l2_norm = ZERO;

    /* Sum the square of the L2-norm for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          l2_norm += sqr(L2_Norm_Residual(Soln_ptr[i]));
       } /* endif */
    }  /* endfor */

    /* Calculate the L2-norm for all blocks. */

    l2_norm = sqrt(l2_norm);

    /* Return the L2-norm. */

    return (l2_norm);

}

/********************************************************
 * Routine: Max_Norm_Residual                           *
 *                                                      *
 * Determines the maximum norm of the solution residual *
 * for a 1D array of 2D quadrilateral multi-block       *
 * solution blocks.  Useful for monitoring convergence  *
 * of the solution for steady state problems.           *
 *                                                      *
 ********************************************************/
double Max_Norm_Residual(Ion5Moment2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List,
			 Ion5Moment2D_Input_Parameters &Input_Parameters) {

    int i;
    double max_norm;

    max_norm = ZERO;

    /* Find the maximum norm for all solution blocks. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          max_norm = max(max_norm, Max_Norm_Residual(Soln_ptr[i]));
       } /* endif */
    }  /* endfor */

    /* Return the maximum norm. */

    return (max_norm);

}

/********************************************************
 * Routine:  Residual_Smoothing                         *
 *                                                      *
 * Applies implicit residual smoothing to a 1D array of *
 * 2D quadrilateral multi-block  solution blocks.       *
 *                                                      *
 ********************************************************/
void Residual_Smoothing(Ion5Moment2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        Ion5Moment2D_Input_Parameters &Input_Parameters,
                        const int I_Stage) {

    int i, k_residual;

    switch(Input_Parameters.i_Time_Integration) {
      case TIME_STEPPING_EXPLICIT_EULER :
        k_residual = 0;
        break;
      case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
        k_residual = 0;
        break;
      case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
        k_residual = 0;
	if (Input_Parameters.N_Stage == 4) {
	   if (I_Stage == 4) {
	      k_residual = 0;
           } else {
	      k_residual = I_Stage - 1;
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

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
	 Residual_Smoothing(Soln_ptr[i],
                            k_residual,
			    Input_Parameters.Residual_Smoothing_Epsilon, 
                            Input_Parameters.Residual_Smoothing_Gauss_Seidel_Iterations);
       } /* endif */
    }  /* endfor */

}

/****************************************************************
 * Routine: Apply_Boundary_Flux_Correction                      *
 *                                                              *
 * Apply flux corrections at boundaries of a 1D array           *
 * of 2D quadrilateral multi-block solution blocks to           *
 * ensure that the scheme is conservative at boundaries         *
 * with resolution mesh changes.                                *
 *                                                              *
 ****************************************************************/
void Apply_Boundary_Flux_Corrections(Ion5Moment2D_Quad_Block *Soln_ptr,
                                     AdaptiveBlock2D_List &Soln_Block_List,
                                     Ion5Moment2D_Input_Parameters &Input_Parameters) {

    int i;

    /* Prescribe boundary data for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          Apply_Boundary_Flux_Corrections(Soln_ptr[i],
                                          Soln_Block_List.Block[i].nN,
                                          Soln_Block_List.Block[i].nS,
                                          Soln_Block_List.Block[i].nE,
                                          Soln_Block_List.Block[i].nW);
       } /* endif */
    }  /* endfor */

}

/****************************************************************
 * Routine: Apply_Boundary_Flux_Corrections_Multistage_Explicit *
 *                                                              *
 * Apply flux corrections at boundaries of a 1D array           *
 * of 2D quadrilateral multi-block solution blocks to           *
 * ensure that the scheme is conservative at boundaries         *
 * with resolution mesh changes.                                *
 *                                                              *
 ****************************************************************/
void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Ion5Moment2D_Quad_Block *Soln_ptr,
                                                         AdaptiveBlock2D_List &Soln_Block_List,
                                                         Ion5Moment2D_Input_Parameters &Input_Parameters,
   	                                                 const int I_Stage) {

    int i;

    /* Prescribe boundary data for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          Apply_Boundary_Flux_Corrections_Multistage_Explicit(Soln_ptr[i],
                                                              I_Stage,
                                                              Input_Parameters.N_Stage,
                                                              Input_Parameters.CFL_Number,
                                                              Input_Parameters.i_Time_Integration,
                                                              Input_Parameters.Local_Time_Stepping,
                                                              Input_Parameters.i_Reconstruction,
                                                              Input_Parameters.i_Limiter, 
                                                              Input_Parameters.i_Flux_Function,
                                                              Soln_Block_List.Block[i].nN,
                                                              Soln_Block_List.Block[i].nS,
                                                              Soln_Block_List.Block[i].nE,
                                                              Soln_Block_List.Block[i].nW);
       } /* endif */
    }  /* endfor */

}

/********************************************************
 * Routine: dUdt_Multistage_Explicit                    *
 *                                                      *
 * This routine evaluates the stage solution residual   *
 * for a 1D array of 2D quadrilateral multi-block       *
 * solution blocks.  A variety of multistage explicit   *
 * time integration and upwind finite-volume spatial    *
 * discretization procedures can be used depending on   *
 * the specified input values.                          *
 *                                                      *
 ********************************************************/
int dUdt_Multistage_Explicit(Ion5Moment2D_Quad_Block *Soln_ptr,
			     AdaptiveBlockResourceList &Global_Soln_Block_List,
                             AdaptiveBlock2D_List &Local_Soln_Block_List,
                             Ion5Moment2D_Input_Parameters &Input_Parameters,
   	                     const int I_Stage) {

    int i, error_flag;

    error_flag = 0;

    /* Update the solution for each solution block. */

    for ( i = 0 ; i <= Local_Soln_Block_List.Nblk-1 ; ++i ) {
       if (Local_Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          error_flag = dUdt_Multistage_Explicit(Soln_ptr[i],
                                                I_Stage,
                                                Input_Parameters);
          if (error_flag) return (error_flag);
       } /* endif */
    }  /* endfor */

    /* Quadrilateral multi-block solution blocks
       successfully updated.  Return. */

    return(error_flag);

}

/********************************************************
 * Routine: Update_Solution_Multistage_Explicit         *
 *                                                      *
 * This routine updates the solution for a 1D array of  *
 * 2D quadrilateral multi-block solution blocks.  A     *
 * variety of multistage explicit time integration      *
 * and upwind finite-volume spatial discretization      *
 * procedures can be used depending on the specified    *
 * input values.                                        *
 *                                                      *
 ********************************************************/
int Update_Solution_Multistage_Explicit(Ion5Moment2D_Quad_Block *Soln_ptr,
                                        AdaptiveBlock2D_List &Soln_Block_List,
                                        Ion5Moment2D_Input_Parameters &Input_Parameters,
   	                                const int I_Stage) {

    int i, error_flag;

    error_flag = 0;

    /* Update the solution for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          error_flag = Update_Solution_Multistage_Explicit(Soln_ptr[i],
                                                           I_Stage,
                                                           Input_Parameters);
          if (error_flag) return (error_flag);
       } /* endif */
    }  /* endfor */

    /* Quadrilateral multi-block solution blocks
       successfully updated.  Return. */

    return(error_flag);

}

/********************************************************
 * Routine: dUdt_Output_Cells_Tecplot                   *
 *                                                      *
 * Writes the cell centred solution changes for a 1D    *
 * array of 2D quadrilateral multi-block solution       *
 * blocks to the specified output data file(s) in a     *
 * format suitable for plotting with TECPLOT.           *
 * Returns a non-zero value if cannot write any of the  *
 * TECPLOT solution files.                              *
 *                                                      *
 ********************************************************/
int dUdt_Output_Cells_Tecplot(Ion5Moment2D_Quad_Block *Soln_ptr,
                              AdaptiveBlock2D_List &Soln_Block_List,
                              Ion5Moment2D_Input_Parameters &Input_Parameters,
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
    strcat(prefix, "_dudt_cells_cpu");

    /* Determine output data file name for this processor. */

    sprintf(extension, "%.6d", Soln_Block_List.ThisCPU);
    strcat(extension, ".dat");
    strcpy(output_file_name, prefix);
    strcat(output_file_name, extension);
    output_file_name_ptr = output_file_name;

    /* Open the output data file. */

    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) return (1);

    /* Write the solution change data for each solution block. */

    i_output_title = 1;
    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          dUdt_Output_Cells_Tecplot(Soln_ptr[i], 
                                    Number_of_Time_Steps, 
                                    Time,
                                    Soln_Block_List.Block[i].gblknum,
                                    i_output_title,
                                    Input_Parameters.i_Reconstruction,
                                    Input_Parameters.i_Limiter,
                                    Input_Parameters.i_Flux_Function,
                                    output_file);
	  if (i_output_title) i_output_title = 0;
       } /* endif */
    }  /* endfor */

    /* Close the output data file. */

    output_file.close();

    /* Writing of output data files complete.  Return zero value. */

    return(0);

}
