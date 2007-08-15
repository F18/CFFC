/* AdvectDiffuse2DQuadMultiBlock.cc:  Multi-Block Versions of Subroutines 
                                      for 2D Advection Diffusion Equation
                                      Multi-Block Quadrilateral Mesh 
                                      Solution Classes. */

/* Include 2D advection diffusion equation quadrilateral mesh solution header file. */

#ifndef _ADVECTDIFFUSE2D_QUAD_INCLUDED
#include "AdvectDiffuse2DQuad.h"
#endif // _ADVECTDIFFUSE2D_QUAD_INCLUDED

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block -- Multiple Block External Subroutines.     *
 **************************************************************************/

/********************************************************
 * Routine: Allocate                                    *
 *                                                      *
 * Allocate memory for 1D array of 2D quadrilateral     *
 * multi-block solution blocks.                         *
 *                                                      *
 ********************************************************/
AdvectDiffuse2D_Quad_Block* Allocate(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                                     AdvectDiffuse2D_Input_Parameters &Input_Parameters) {

    /* Allocate memory. */

    Soln_ptr = new AdvectDiffuse2D_Quad_Block[Input_Parameters.Number_of_Blocks_Per_Processor];

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
AdvectDiffuse2D_Quad_Block* Deallocate(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                                       AdvectDiffuse2D_Input_Parameters &Input_Parameters) {

    int i;
 
    /* Deallocate memory. */

    for ( i = 0 ; i <= Input_Parameters.Number_of_Blocks_Per_Processor-1 ; ++i ) {
       if (Soln_ptr[i].U != NULL) Soln_ptr[i].deallocate();
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
void ICs(AdvectDiffuse2D_Quad_Block *Soln_ptr,
         AdaptiveBlock2D_List &Soln_Block_List,
         AdvectDiffuse2D_Input_Parameters &Input_Parameters,
	 QuadTreeBlock_DataStructure &QuadTree) {

    int i;
    AdvectDiffuse2D_State Uo[5];
    int I = Input_Parameters.Number_of_Blocks_Idir;
    int J = Input_Parameters.Number_of_Blocks_Jdir;

    /* Define various reference flow states. */

    Uo[0] = Input_Parameters.Uo;
    Uo[1] = Input_Parameters.U1;
    Uo[2] = Input_Parameters.U2;

    /* Assign initial data for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
	  // Set flow geometry indicator (planar/axisymmetric) flow block.
          Soln_ptr[i].Axisymmetric = Input_Parameters.Axisymmetric;

          // Set initial data.
          ICs(Soln_ptr[i], 
              Input_Parameters.i_ICs,
              Uo);

 //          // Zero the gradients and residuals. 
//           Zero_Gradients_and_Residuals(Soln_ptr[i]);

	  // Copy Initial data and other overridden boundary data into
	  // solution block boundary reference states
	  
	  Set_Boundary_Ref_State(Soln_ptr[i], 
				 Soln_Block_List.Block[i].gblknum,
				 QuadTree,
				 Input_Parameters);
       } /* endif */
    }  /* endfor */
}

/********************************************************
 * Routine: Set_Advection_Velocity_Field                *
 *                                                      *
 * Specify the advection velocity field for a 1D array  *
 * of 2D quadrilateral multi-block solution blocks.     *
 *                                                      *
 ********************************************************/
void Set_Advection_Velocity_Field(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                                  AdaptiveBlock2D_List &Soln_Block_List,
                                  AdvectDiffuse2D_Input_Parameters &Input_Parameters) {

    int i;

    /* Specify the advection velocity field for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          Set_Advection_Velocity_Field(Soln_ptr[i], 
                                       Input_Parameters.i_Velocity_Field,
                                       Input_Parameters.a,
                                       Input_Parameters.b);
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
int Read_Restart_Solution(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          AdvectDiffuse2D_Input_Parameters &Input_Parameters,
                          int &Number_of_Time_Steps,
                          double &Time,
                          CPUTime &CPU_Time) {

    int i, i_new_time_set, nsteps;
    char prefix[256], extension[256], restart_file_name[256], gas_type[256];
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
          if (!i_new_time_set) {
             Number_of_Time_Steps = nsteps;
             Input_Parameters.Maximum_Number_of_Time_Steps += Number_of_Time_Steps;
	     Time = time0;
             CPU_Time.cput = cpu_time0.cput;
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
int Write_Restart_Solution(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                           AdaptiveBlock2D_List &Soln_Block_List,
                           AdvectDiffuse2D_Input_Parameters &Input_Parameters,
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
          restart_file << setprecision(14) << Soln_ptr[i];

          // Close restart file.
          restart_file.close();
       } /* endif */
    }  /* endfor */

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
int Output_Tecplot(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                   AdaptiveBlock2D_List &Soln_Block_List,
                   AdvectDiffuse2D_Input_Parameters &Input_Parameters,
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
int Output_Cells_Tecplot(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List,
                         AdvectDiffuse2D_Input_Parameters &Input_Parameters,
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
 * Routine: Output_Nodes_Tecplot                        *
 *                                                      *
 * Writes the node values for a 1D array of 2D          *
 * quadrilateral multi-block solution blocks to the     *
 * specified output data file(s) in a format suitable   *
 * for plotting with TECPLOT.  Returns a non-zero value *
 * if cannot write any of the TECPLOT solution files.   *
 *                                                      *
 ********************************************************/
int Output_Nodes_Tecplot(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List,
                         AdvectDiffuse2D_Input_Parameters &Input_Parameters,
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
    strcat(prefix, "_nodes_cpu");

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
          Output_Nodes_Tecplot(Soln_ptr[i], 
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
int Output_Mesh_Tecplot(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        AdvectDiffuse2D_Input_Parameters &Input_Parameters,
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
int Output_Mesh_Gnuplot(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        AdvectDiffuse2D_Input_Parameters &Input_Parameters,
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
void BCs(AdvectDiffuse2D_Quad_Block *Soln_ptr,
         AdaptiveBlock2D_List &Soln_Block_List,
	 AdvectDiffuse2D_Input_Parameters &Input_Parameters) {

    int i;

    /* Prescribe boundary data for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          BCs(Soln_ptr[i]);
       } /* endif */
    }  /* endfor */

}

/********************************************************
 * Routine: Set_Analytical_Solution                     *
 *                                                      *
 * Set Analytical Solution for each solution block from *
 * the list of local solution blocks.                   *
 ********************************************************/
void Set_Analytical_Solution(AdvectDiffuse2D_Quad_Block *Soln_ptr,
			     const AdaptiveBlock2D_List &Soln_Block_List,
			     const AdvectDiffuse2D_Input_Parameters &Input_Parameters) {
    int i;

    /* Prescribe boundary data for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
	  // Set Analytical Solution
	  Set_Analytical_Solution(Soln_ptr[i],
				  Input_Parameters);
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
 * condition for the advection terms a semi-impirical   *
 * criteria for the diffusion and source terms.         *
 *                                                      *
 ********************************************************/
double CFL(AdvectDiffuse2D_Quad_Block *Soln_ptr,
           AdaptiveBlock2D_List &Soln_Block_List,
           AdvectDiffuse2D_Input_Parameters &Input_Parameters) {

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
void Set_Global_TimeStep(AdvectDiffuse2D_Quad_Block *Soln_ptr,
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
double L1_Norm_Residual(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List) {

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
double L2_Norm_Residual(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List) {

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
double Max_Norm_Residual(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List) {

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
 * Routine: Evaluate_Limiters                           *
 *                                                      *
 * Set conditions to evaluate the limiters for a        *
 * 1D array of 2D quadrilateral multi-block solution    *
 * blocks.                                              *
 *                                                      *
 ********************************************************/
void Evaluate_Limiters(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                       AdaptiveBlock2D_List &Soln_Block_List) {

    int i;

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) Soln_ptr[i].evaluate_limiters();
    }  /* endfor */

}

/********************************************************
 * Routine: Freeze_Limiters                             *
 *                                                      *
 * Set conditions to freeze the limiters for a          *
 * 1D array of 2D quadrilateral multi-block solution    *
 * blocks.                                              *
 *                                                      *
 ********************************************************/
void Freeze_Limiters(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                     AdaptiveBlock2D_List &Soln_Block_List) {

    int i;

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) Soln_ptr[i].freeze_limiters();
    }  /* endfor */

}


/********************************************************
 * Routine:  Residual_Smoothing                         *
 *                                                      *
 * Applies implicit residual smoothing to a 1D array of *
 * 2D quadrilateral multi-block  solution blocks.       *
 *                                                      *
 ********************************************************/
void Residual_Smoothing(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        AdvectDiffuse2D_Input_Parameters &Input_Parameters,
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
void Apply_Boundary_Flux_Corrections(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                                     AdaptiveBlock2D_List &Soln_Block_List) {

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
void Apply_Boundary_Flux_Corrections_Multistage_Explicit(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                                                         AdaptiveBlock2D_List &Soln_Block_List,
                                                         AdvectDiffuse2D_Input_Parameters &Input_Parameters,
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
                                                              Input_Parameters.i_Reconstruction,
                                                              Input_Parameters.i_Limiter, 
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
 * time integration and a 2nd-ororder limited upwind    *
 * finite-volume spatial discretization scheme for the  *
 * convective flux coupled with a centrally-weighted    *
 * finite-volume discretization for the diffused flux   *
 * can be used depending on the specified input values. *
 *                                                      *
 ********************************************************/
int dUdt_Multistage_Explicit(AdvectDiffuse2D_Quad_Block *Soln_ptr,
			     AdaptiveBlockResourceList &Global_Soln_Block_List,
                             AdaptiveBlock2D_List &Local_Soln_Block_List,
                             AdvectDiffuse2D_Input_Parameters &Input_Parameters,
   	                     const int I_Stage) {

    int i, error_flag;

    error_flag = 0;

    /* Evaluate the stage solution residual for each solution block. */

    for ( i = 0 ; i <= Local_Soln_Block_List.Nblk-1 ; ++i ) {
       if (Local_Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
	
	 error_flag = dUdt_Multistage_Explicit(Soln_ptr[i],
					       I_Stage,
					       Input_Parameters);
	 if (error_flag) return (error_flag);
       } /* endif */
    }  /* endfor */

    /* Residuals for each quadrilateral multi-block solution block
       successfully calculated.  Return. */

    return(error_flag);

}

/********************************************************
 * Routine: Update_Solution_Multistage_Explicit         *
 *                                                      *
 * This routine updates the solution for a 1D array of  *
 * 2D quadrilateral multi-block solution blocks.        *
 * A variety of multistage explicit                     *
 * time integration and a 2nd-ororder limited upwind    *
 * finite-volume spatial discretization scheme for the  *
 * convective flux coupled with a centrally-weighted    *
 * finite-volume discretization for the diffused flux   *
 * can be used depending on the specified input values. *
 *                                                      *
 ********************************************************/
int Update_Solution_Multistage_Explicit(AdvectDiffuse2D_Quad_Block *Soln_ptr,
                                        AdaptiveBlock2D_List &Soln_Block_List,
                                        AdvectDiffuse2D_Input_Parameters &Input_Parameters,
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

/**********************************************************
 * Routine: Initial_AMR2                                  *
 *                                                        *
 * Performs initial refinement of the adaptive solution   *
 * block mesh based on initial data.  Returns a zero      *
 * value if no error in the AMR precedure has occurred.   *
 *                                                        *
 **********************************************************/
int Initial_AMR2(AdvectDiffuse2D_Quad_Block       *Soln_ptr,
                 AdvectDiffuse2D_Input_Parameters &InputParameters,
                 QuadTreeBlock_DataStructure      &QuadTree,
                 AdaptiveBlockResourceList        &GlobalSolnBlockList,
                 AdaptiveBlock2D_List             &LocalSolnBlockList) {

    int error_flag, number_of_initial_mesh_refinements;

    error_flag = 0;

    if (InputParameters.Number_of_Initial_Mesh_Refinements == 0) return(error_flag);

    for (number_of_initial_mesh_refinements = 1; 
         number_of_initial_mesh_refinements <= InputParameters.Number_of_Initial_Mesh_Refinements; 
         ++number_of_initial_mesh_refinements) {

       error_flag = AMR(Soln_ptr,
			InputParameters,
   	                QuadTree,
 	                GlobalSolnBlockList,
	                LocalSolnBlockList,
	                ON,ON);
       if (error_flag) return (error_flag);

       ICs(Soln_ptr,
           LocalSolnBlockList,
           InputParameters,
	   QuadTree);

       Set_Advection_Velocity_Field(Soln_ptr,
                                    LocalSolnBlockList,
                                    InputParameters);

       BCs(Soln_ptr,
           LocalSolnBlockList,
	   InputParameters);

       error_flag = Send_All_Messages(Soln_ptr,
                                      LocalSolnBlockList,
                                      Soln_ptr[0].NumVar(),
                                      OFF);
       if (error_flag) return (error_flag);

       if (CFFC_Primary_MPI_Processor()) {
          cout << "\n Refinement Level #" << number_of_initial_mesh_refinements
               << " : Number of Blocks = " << QuadTree.countUsedBlocks()
               << ", Number of Cells = " << QuadTree.countUsedCells()
               << ", Refinement Efficiency = " << QuadTree.efficiencyRefinement();
       } /* endif */

    } /* endfor */

    return(error_flag);

}

/**********************************************************
 * Routine: Uniform_AMR2                                  *
 *                                                        *
 * Performs uniform refinement of the adaptive solution   *
 * block mesh.  Returns a zero value if no error in the   *
 * AMR precedure has occurred.                            *
 *                                                        *
 **********************************************************/
int Uniform_AMR2(AdvectDiffuse2D_Quad_Block       *Soln_ptr,
                 AdvectDiffuse2D_Input_Parameters &InputParameters,
                 QuadTreeBlock_DataStructure      &QuadTree,
                 AdaptiveBlockResourceList        &GlobalSolnBlockList,
                 AdaptiveBlock2D_List             &LocalSolnBlockList) {

    int error_flag, number_of_uniform_mesh_refinements;

    error_flag = 0;

    if (InputParameters.Number_of_Uniform_Mesh_Refinements == 0) return(error_flag);

    for (number_of_uniform_mesh_refinements = 1; 
         number_of_uniform_mesh_refinements <= InputParameters.Number_of_Uniform_Mesh_Refinements; 
         ++number_of_uniform_mesh_refinements) {

       // Set refinement flags to all.
       QuadTree.nochangeAll();
       LocalSolnBlockList.refineAll();

       error_flag = AMR(Soln_ptr,
			InputParameters,
   	                QuadTree,
 	                GlobalSolnBlockList,
	                LocalSolnBlockList,
	                OFF,ON);
       if (error_flag) return (error_flag);

       ICs(Soln_ptr,
           LocalSolnBlockList,
           InputParameters,
	   QuadTree);

       Set_Advection_Velocity_Field(Soln_ptr,
                                    LocalSolnBlockList,
                                    InputParameters);

       BCs(Soln_ptr,
           LocalSolnBlockList,
	   InputParameters);

       error_flag = Send_All_Messages(Soln_ptr,
                                      LocalSolnBlockList,
                                      Soln_ptr[0].NumVar(),
                                      OFF);
       if (error_flag) return (error_flag);

       if (CFFC_Primary_MPI_Processor()) {
          cout << "\n Refinement Level #" << number_of_uniform_mesh_refinements
               << " : Number of Blocks = " << QuadTree.countUsedBlocks()
               << ", Number of Cells = " << QuadTree.countUsedCells()
               << ", Refinement Efficiency = " << QuadTree.efficiencyRefinement();
       } /* endif */

    } /* endfor */

    return(error_flag);

}

/**********************************************************
 * Routine: Boundary_AMR2                                 *
 *                                                        *
 * Performs refinement of the adaptive solution block     *
 * mesh based on boundary conditions data.  Returns a     *
 * zero value if no error in the AMR precedure has        *
 * occurred.                                              *
 *                                                        *
 **********************************************************/
int Boundary_AMR2(AdvectDiffuse2D_Quad_Block       *Soln_ptr,
                  AdvectDiffuse2D_Input_Parameters &InputParameters,
                  QuadTreeBlock_DataStructure      &QuadTree,
                  AdaptiveBlockResourceList        &GlobalSolnBlockList,
                  AdaptiveBlock2D_List             &LocalSolnBlockList) {

    int i, j, nb, error_flag, number_of_boundary_mesh_refinements;

    error_flag = 0;

    if (InputParameters.Number_of_Boundary_Mesh_Refinements == 0) return(error_flag);

    for (number_of_boundary_mesh_refinements = 1; 
         number_of_boundary_mesh_refinements <= InputParameters.Number_of_Boundary_Mesh_Refinements; 
         ++number_of_boundary_mesh_refinements) {

       // Set refinement flags to all.
       QuadTree.nochangeAll();
       for (nb = 0; nb < LocalSolnBlockList.Nblk; ++nb) {
          LocalSolnBlockList.RefineFlag[nb] = ADAPTIVEBLOCK2D_NOCHANGE;
	  if (LocalSolnBlockList.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	     for (i = Soln_ptr[nb].ICl-Soln_ptr[nb].Nghost; i <= Soln_ptr[nb].ICu+Soln_ptr[nb].Nghost; i++) {
	        if (Soln_ptr[nb].Grid.BCtypeN[i] == BC_DIRICHLET ||
		    Soln_ptr[nb].Grid.BCtypeN[i] == BC_NEUMANN) {
		   LocalSolnBlockList.RefineFlag[nb] = ADAPTIVEBLOCK2D_REFINE;
		} /* endif */
	     } /* endfor */
	     for (j = Soln_ptr[nb].JCl-Soln_ptr[nb].Nghost; j <= Soln_ptr[nb].JCu+Soln_ptr[nb].Nghost; j++) {
	        if (Soln_ptr[nb].Grid.BCtypeE[j] == BC_DIRICHLET ||
		    Soln_ptr[nb].Grid.BCtypeE[j] == BC_NEUMANN) {
		   LocalSolnBlockList.RefineFlag[nb] = ADAPTIVEBLOCK2D_REFINE;
		} /* endif */
	     } /* endfor */
	  } /* endif */
       } /* endfor */

       error_flag = AMR(Soln_ptr,
			InputParameters,
   	                QuadTree,
 	                GlobalSolnBlockList,
	                LocalSolnBlockList,
	                OFF,ON);
       if (error_flag) return (error_flag);

       ICs(Soln_ptr,
           LocalSolnBlockList,
           InputParameters,
	   QuadTree);

       Set_Advection_Velocity_Field(Soln_ptr,
                                    LocalSolnBlockList,
                                    InputParameters);

       BCs(Soln_ptr,
           LocalSolnBlockList,
	   InputParameters);

       error_flag = Send_All_Messages(Soln_ptr,
                                      LocalSolnBlockList,
                                      Soln_ptr[0].NumVar(),
                                      OFF);
       if (error_flag) return (error_flag);

       if (CFFC_Primary_MPI_Processor()) {
          cout << "\n Refinement Level #" << number_of_boundary_mesh_refinements
               << " : Number of Blocks = " << QuadTree.countUsedBlocks()
               << ", Number of Cells = " << QuadTree.countUsedCells()
               << ", Refinement Efficiency = " << QuadTree.efficiencyRefinement();
       } /* endif */

    } /* endfor */

    return(error_flag);

}
