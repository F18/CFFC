/* Rte2DQuadMultiBlock.cc:  Multi-Block Versions of Subroutines for 2D Rte 
                              Multi-Block Quadrilateral Mesh 
                              Solution Classes. */

/* Include 2D Rte quadrilateral mesh solution header file. */
#include "Rte2DQuad.h"


/**************************************************************************
 * Rte2D_Quad_Block -- Multiple Block External Subroutines.             *
 **************************************************************************/

/********************************************************
 * Routine: Allocate                                    *
 *                                                      *
 * Allocate memory for 1D array of 2D quadrilateral     *
 * multi-block solution blocks.                         *
 *                                                      *
 ********************************************************/
Rte2D_Quad_Block* Allocate(Rte2D_Quad_Block *Soln_ptr,
                             Rte2D_Input_Parameters &Input_Parameters) {

    /* Allocate memory. */

    Soln_ptr = new Rte2D_Quad_Block[Input_Parameters.Number_of_Blocks_Per_Processor];

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
Rte2D_Quad_Block* Deallocate(Rte2D_Quad_Block *Soln_ptr,
                               Rte2D_Input_Parameters &Input_Parameters) {

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

/**********************************************************************
 * Routine: CreateInitialSolutionBlocks                               *
 *                                                                    *
 * This routine calls the Create_Initial_Solution_Blocks routine      *
 * found in AMR.h.                                                    *
 *                                                                    *
 **********************************************************************/
Rte2D_Quad_Block* CreateInitialSolutionBlocks(Grid2D_Quad_Block **InitMeshBlk,
					      Rte2D_Quad_Block *Soln_ptr,
					      Rte2D_Input_Parameters &Input_Parameters,
					      QuadTreeBlock_DataStructure &QuadTree,
					      AdaptiveBlockResourceList &GlobalSolnBlockList,
					      AdaptiveBlock2D_List &LocalSolnBlockList) {
  
  int error_flag;

  // Create the initial solution blocks.
  Soln_ptr = Create_Initial_Solution_Blocks(InitMeshBlk,
					    Soln_ptr,
					    Input_Parameters,
					    QuadTree,
					    GlobalSolnBlockList,
					    LocalSolnBlockList);

  // Return the solution blocks.
  return Soln_ptr;

}

/********************************************************
 * Routine: ICs                                         *
 *                                                      *
 * Assigns initial conditions and data to the           *
 * solution variables of a 1D array of 2D quadrilateral *
 * multi-block solution blocks.                         *
 *                                                      *
 ********************************************************/
void ICs(Rte2D_Quad_Block *Soln_ptr,
         AdaptiveBlock2D_List &Soln_Block_List,
         Rte2D_Input_Parameters &Input_Parameters) {

    int i;


    /* Assign initial data for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
	  // Set flow geometry indicator (planar/axisymmetric) flow block.
 	  Soln_ptr[i].Axisymmetric = Input_Parameters.Axisymmetric;
	  Soln_ptr[i].Medium_Field_Type = Input_Parameters.Medium_Field_Type;
	  Soln_ptr[i].NorthWallTemp  = Input_Parameters.NorthWallTemp;   
	  Soln_ptr[i].SouthWallTemp  = Input_Parameters.SouthWallTemp;  
	  Soln_ptr[i].EastWallTemp   = Input_Parameters.EastWallTemp;   
	  Soln_ptr[i].WestWallTemp   = Input_Parameters.WestWallTemp;  
	  Soln_ptr[i].NorthWallEmiss = Input_Parameters.NorthWallEmiss;   
	  Soln_ptr[i].SouthWallEmiss = Input_Parameters.SouthWallEmiss;  
	  Soln_ptr[i].EastWallEmiss  = Input_Parameters.EastWallEmiss;   
	  Soln_ptr[i].WestWallEmiss  = Input_Parameters.WestWallEmiss;

          // Set initial data.
          ICs(Soln_ptr[i], Input_Parameters);

	  //
          // Set initial data for medium
	  //
	  // For an analytically prescribed field
	  if (Soln_ptr[i].Medium_Field_Type == MEDIUM2D_FIELD_ANALYTIC)
	    PrescribeFields(Soln_ptr[i]);
	  //
	  // for a discrete field specified by the user
	  // do nothing
	  
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
int Read_Restart_Solution(Rte2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          Rte2D_Input_Parameters &Input_Parameters,
                          int &Number_of_Time_Steps,
                          double &Time,
                          CPUTime &CPU_Time) {

    int i, i_new_time_set, nsteps;
    char prefix[256], extension[256], restart_file_name[256], gas_type[256];
    char *restart_file_name_ptr;
    ifstream restart_file;
    double time0;
    CPUTime cpu_time0;
    double temp;

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
          if (restart_file.fail()) return (1);

          // Read solution block data.
          restart_file.setf(ios::skipws);
          restart_file >> nsteps >> time0 >> cpu_time0;
          restart_file.unsetf(ios::skipws);

	  //------------------------ Rte2D Specific -------------------------//
          restart_file.setf(ios::skipws);
	  restart_file >> Input_Parameters.i_RTE_Solver;	    
	  restart_file >> Input_Parameters.Number_of_Angles_Mdir;	    
	  restart_file >> Input_Parameters.Number_of_Angles_Ldir;	    
	  restart_file >> Input_Parameters.i_DOM_Quadrature;	    
	  restart_file >> Input_Parameters.i_AbsorptionModel;	    
	  restart_file >> Input_Parameters.SNBCK_IP;	    
	  restart_file >> Input_Parameters.Axisymmetric;	    
          restart_file.unsetf(ios::skipws);

	  // Setup conserved and medium state
	  Input_Parameters.SetupInputState();
	  //---------------------- End Rte2D Specific -----------------------//

          if (!i_new_time_set) {
             Number_of_Time_Steps = nsteps;
             Input_Parameters.Maximum_Number_of_Time_Steps += Number_of_Time_Steps;
	     Time = time0;
             CPU_Time.cput = cpu_time0.cput;

             i_new_time_set = 1;
          } /* endif */

	  /************** Quad Single Block Data from file *********************/
          restart_file >> Soln_ptr[i];

          // Close restart file.
          restart_file.close();

	  //------------------------ Rte2D Specific -------------------------//
	  //OVERIDE Restart File settings with Input File Settings (may have changed)
	  // Soln_ptr[i].Flow_Type = Input_Parameters.FlowType; // <-- not used (static)
	  Soln_ptr[i].Axisymmetric = Input_Parameters.Axisymmetric;
	  Soln_ptr[i].Medium_Field_Type = Input_Parameters.Medium_Field_Type;
	  //---------------------- End Rte2D Specific -----------------------//

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
int Write_Restart_Solution(Rte2D_Quad_Block *Soln_ptr,
                           AdaptiveBlock2D_List &Soln_Block_List,
                           Rte2D_Input_Parameters &Input_Parameters,
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
          if (restart_file.fail()) return (1);

          // Write solution block data.
          restart_file.setf(ios::scientific);
          restart_file << setprecision(14) << Number_of_Time_Steps 
                       << " " << Time << " " << CPU_Time << "\n";

	  //------------------------ Rte2D Specific -------------------------//

	  restart_file << Input_Parameters.i_RTE_Solver << "\n";	    
	  restart_file << Input_Parameters.Number_of_Angles_Mdir << "\n";	    
	  restart_file << Input_Parameters.Number_of_Angles_Ldir << "\n";	    
	  restart_file << Input_Parameters.DOM_Quadrature << "\n";	    
	  restart_file << Input_Parameters.i_AbsorptionModel << "\n";
	  restart_file << Input_Parameters.SNBCK_IP << "\n";	    
	  restart_file << Input_Parameters.Axisymmetric << "\n";	    

	  //---------------------- End Rte2D Specific -----------------------//
	  
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
int Output_Tecplot(Rte2D_Quad_Block *Soln_ptr,
                   AdaptiveBlock2D_List &Soln_Block_List,
                   Rte2D_Input_Parameters &Input_Parameters,
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
    if (output_file.fail()) return (1);

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
int Output_Cells_Tecplot(Rte2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List,
                         Rte2D_Input_Parameters &Input_Parameters,
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
    if (output_file.fail()) return (1);

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
 * Routine: Output_Nodes_Tecplot                        *
 *                                                      *
 * Writes the node values for a 1D array of 2D          *
 * quadrilateral multi-block solution blocks to the     *
 * specified output data file(s) in a format suitable   *
 * for plotting with TECPLOT.  Returns a non-zero value *
 * if cannot write any of the TECPLOT solution files.   *
 *                                                      *
 ********************************************************/
int Output_Nodes_Tecplot(Rte2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List,
                         Rte2D_Input_Parameters &Input_Parameters,
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
    if (output_file.fail()) return (1);

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
int Output_Gradients_Tecplot(Rte2D_Quad_Block *Soln_ptr,
			     AdaptiveBlock2D_List &Soln_Block_List,
			     Rte2D_Input_Parameters &IP,
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
  if (output_file.fail()) return 1;

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
int Output_Mesh_Tecplot(Rte2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        Rte2D_Input_Parameters &Input_Parameters,
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
    if (output_file.fail()) return (1);

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
int Output_Mesh_Gnuplot(Rte2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        Rte2D_Input_Parameters &Input_Parameters,
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
    if (output_file.fail()) return (1);

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
void BCs(Rte2D_Quad_Block *Soln_ptr,
         AdaptiveBlock2D_List &Soln_Block_List,
	 Rte2D_Input_Parameters &IP) {

    int i;

    /* Prescribe boundary data for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
	 BCs(Soln_ptr[i], IP);
       } /* endif */
    }  /* endfor */

}

void BCs_Space_March(Rte2D_Quad_Block *Soln_ptr,
		     AdaptiveBlock2D_List &Soln_Block_List,
		     Rte2D_Input_Parameters &IP) {

    int i;

    /* Prescribe boundary data for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
	 BCs_Space_March(Soln_ptr[i], IP);
       } /* endif */
    }  /* endfor */

}


/********************************************************
 * Routine: CFL                                         *
 *                                                      *
 * Determines the allowable global and local time steps *
 * (for explicit Rte time stepping scheme) for a 1D   *
 * array of 2D quadrilateral multi-block solution       *
 * blocks according to the Courant-Friedrichs-Lewy      *
 * condition.                                           *
 *                                                      *
 ********************************************************/
double CFL(Rte2D_Quad_Block *Soln_ptr,
           AdaptiveBlock2D_List &Soln_Block_List,
           Rte2D_Input_Parameters &Input_Parameters) {

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
void Set_Global_TimeStep(Rte2D_Quad_Block *Soln_ptr,
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
double L1_Norm_Residual(Rte2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List) {

    int i;
    double l1_norm;

    l1_norm = ZERO;

    /* Calculate the L1-norm.
       Sum the L1-norm for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          l1_norm += L1_Norm_Residual(Soln_ptr[i],Soln_ptr[i].residual_variable);
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
double L2_Norm_Residual(Rte2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List) {

    int i;
    double l2_norm;

    l2_norm = ZERO;

    /* Sum the square of the L2-norm for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          l2_norm += sqr(L2_Norm_Residual(Soln_ptr[i],Soln_ptr[i].residual_variable));
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
double Max_Norm_Residual(Rte2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List) {

    int i;
    double max_norm;

    max_norm = ZERO;

    /* Find the maximum norm for all solution blocks. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          max_norm = max(max_norm, Max_Norm_Residual(Soln_ptr[i],Soln_ptr[i].residual_variable));
       } /* endif */
    }  /* endfor */

    /* Return the maximum norm. */

    return (max_norm);

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
void L1_Norm_Residual(Rte2D_Quad_Block *Soln_ptr,
		      AdaptiveBlock2D_List &Soln_Block_List,
		      double *l1_norm) {

  /* Calculate the L1-norm. Sum the L1-norm for each solution block. */
  for( int k = 0; k< Soln_ptr[0].Number_of_Residual_Norms; k++){    
    l1_norm[k] = ZERO;
    for ( int i = 0 ; i < Soln_Block_List.Nblk; ++i ) {
      if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {     
	l1_norm[k] += L1_Norm_Residual(Soln_ptr[i],k+1);       
      }
    }
  }  


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
void L2_Norm_Residual(Rte2D_Quad_Block *Soln_ptr,
		      AdaptiveBlock2D_List &Soln_Block_List,
		      double *l2_norm) {

  /* Sum the square of the L2-norm for each solution block. */
  for( int k = 0; k< Soln_ptr[0].Number_of_Residual_Norms; k++){
    l2_norm[k] =ZERO;
    for ( int i = 0 ; i < Soln_Block_List.Nblk; ++i ) {
      if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
	l2_norm[k] += sqr(L2_Norm_Residual(Soln_ptr[i],k+1));     
      } 
    }  
    l2_norm[k] = sqrt(l2_norm[k]);
  }

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
void Max_Norm_Residual(Rte2D_Quad_Block *Soln_ptr,
		       AdaptiveBlock2D_List &Soln_Block_List,
		       double *max_norm) {
  
  /* Find the maximum norm for all solution blocks. */
  for( int k = 0; k< Soln_ptr[0].Number_of_Residual_Norms; k++){
    max_norm[k] = ZERO;
    for (int i = 0 ; i < Soln_Block_List.Nblk ; ++i ) {
      if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
	max_norm[k] = max(max_norm[k], Max_Norm_Residual(Soln_ptr[i],k+1));
      } 
    } 
  }

}


/********************************************************
 * Routine: Evaluate_Limiters                           *
 *                                                      *
 * Set conditions to evaluate the limiters for a        *
 * 1D array of 2D quadrilateral multi-block solution    *
 * blocks.                                              *
 *                                                      *
 ********************************************************/
void Evaluate_Limiters(Rte2D_Quad_Block *Soln_ptr,
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
void Freeze_Limiters(Rte2D_Quad_Block *Soln_ptr,
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
void Residual_Smoothing(Rte2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        Rte2D_Input_Parameters &Input_Parameters,
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
void Apply_Boundary_Flux_Corrections(Rte2D_Quad_Block *Soln_ptr,
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
void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Rte2D_Quad_Block *Soln_ptr,
                                                         AdaptiveBlock2D_List &Soln_Block_List,
                                                         Rte2D_Input_Parameters &Input_Parameters,
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
int dUdt_Multistage_Explicit(Rte2D_Quad_Block *Soln_ptr,
			     AdaptiveBlockResourceList &Global_Soln_Block_List,
                             AdaptiveBlock2D_List &Soln_Block_List,
                             Rte2D_Input_Parameters &Input_Parameters,
   	                     const int I_Stage) {

    int i, error_flag;

    error_flag = 0;

    /* Evaluate the stage solution residual for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          error_flag = dUdt_Multistage_Explicit(Soln_ptr[i],
                                                I_Stage,
                                                Input_Parameters);
          if (error_flag) return (error_flag);
       } /* endif */
    }  /* endfor */

    /* Residuals for each quadrilateral multi-block solution block
       successfully calcualted.  Return. */

    return(error_flag);

}

/********************************************************
 * Routine: dUdt_Space_March                            *
 *                                                      *
 * This routine evaluates the stage solution residual   *
 * for a 1D array of 2D quadrilateral multi-block       *
 * solution blocks using a space-marching procedure.  A *
 * variety of spatial                                   *
 * discretization procedures can be used depending on   *
 * the specified input values.                          *
 *                                                      *
 ********************************************************/
int dUdt_Space_March(Rte2D_Quad_Block *Soln_ptr,
		     AdaptiveBlockResourceList &Global_Soln_Block_List,
		     AdaptiveBlock2D_List &Soln_Block_List,
		     Rte2D_Input_Parameters &Input_Parameters) {
  
    int i, error_flag;

    error_flag = 0;

    /* Evaluate the stage solution residual for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          error_flag = dUdt_Space_March(Soln_ptr[i],
					Input_Parameters);
          if (error_flag) return (error_flag);
       } /* endif */
    }  /* endfor */

    /* Residuals for each quadrilateral multi-block solution block
       successfully calcualted.  Return. */

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
int Update_Solution_Multistage_Explicit(Rte2D_Quad_Block *Soln_ptr,
                                        AdaptiveBlock2D_List &Soln_Block_List,
                                        Rte2D_Input_Parameters &Input_Parameters,
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



/**********************************************************************
 * Routine: Output_Exact                                              *
 *                                                                    *
 * Writes the exact and computed  solution values                     *
 * at the nodes for a 1D array of 2D quadrilateral multi-block        *
 * solution blocks to the specified output data file(s) in a format   *
 * suitable for plotting with TECPLOT.  Returns a non-zero value if   *
 * cannot write any of the TECPLOT solution files.  This routine will *
 * also compute the error norms of the computed solution (L1, L2, and *
 * max norms).                                                        *
 *                                                                    *
 * The exact solution for multi-dimensional radiative heat flux in a  *
 * cylindrical or rectangular emitting-absorbing isothermal medium    *
 * isothermal bounding walls.  See Dua and Cheng (1975) for           *
 * cylindrical and Cheng (1972) for cartesian coordinate systems.     *
 *                                                                    *
 **********************************************************************/
int Output_Exact(Rte2D_Quad_Block *Soln_ptr,
		 AdaptiveBlock2D_List &Soln_Block_List,
		 Rte2D_Input_Parameters &IP) {

  int i, i_output_title;
  char prefix[256], extension[256], output_file_name[256];
  char *output_file_name_ptr;
  ofstream output_file;    
  double l1_norm, l2_norm, max_norm;
  double l1_temp, l2_temp, max_temp;
  int numberofactivecells, numberofactivecells_temp;

  // Initialize error variables.
  l1_norm = ZERO; l2_norm = ZERO; max_norm = ZERO;
  l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO; 
  numberofactivecells = 0; numberofactivecells_temp = 0;

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
  strcat(prefix,"_exact_cpu");
  
  // Determine output data file name for this processor.
  sprintf(extension,"%.6d",Soln_Block_List.ThisCPU);
  strcat(extension,".dat");
  strcpy(output_file_name,prefix);
  strcat(output_file_name,extension);
  output_file_name_ptr = output_file_name;
  
  // Open the output data file.
  output_file.open(output_file_name_ptr,ios::out);
  if (output_file.fail()) return 1;
  
  // Write the solution data for each solution block.
  i_output_title = 1;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO;
      Output_Exact(Soln_ptr[nb],nb,i_output_title,output_file,
		   l1_temp,l2_temp,max_temp,IP);
      l1_norm += l1_temp;
      l2_norm += sqr(l2_temp);
      max_norm = max(max_norm,max_temp);
      numberofactivecells_temp = (Soln_ptr[nb].Grid.NCi-2*Soln_ptr[nb].Nghost)*(Soln_ptr[nb].Grid.NCj-2*Soln_ptr[nb].Nghost);
      numberofactivecells += numberofactivecells_temp;
      if (i_output_title) i_output_title = 0;
    }
  }
  
  // Close the output data file.
  output_file.close();
  
#ifdef _MPI_VERSION
  l1_norm = CFFC_Summation_MPI(l1_norm);
  l2_norm = CFFC_Summation_MPI(l2_norm);
  max_norm = CFFC_Maximum_MPI(max_norm);
  numberofactivecells = CFFC_Summation_MPI(numberofactivecells);
#endif

  // Calculate the L1-norm and L2-norm for all blocks.
  if (Soln_Block_List.Block[0].used == ADAPTIVEBLOCK2D_USED) {
    l1_norm = l1_norm/double(numberofactivecells);
    l2_norm = sqrt(l2_norm/double(numberofactivecells));
  }

  if (CFFC_Primary_MPI_Processor()) {
    cout << endl
	 << endl
	 << " ==================================================================== "
	 << endl
	 << " Error norms for the black enclosure:"
	 << endl
	 << "   L1_Norm = " << l1_norm
	 << endl
	 << "   L2_Norm = " << l2_norm
	 << endl
	 << "   Max_Norm = " << max_norm
	 << endl
	 << "   Number of cells = " << numberofactivecells
	 << endl
	 << " ==================================================================== "
	 << endl;
  }

  // Writing of output data files complete.  Return zero value.
  return 0;

}


/********************************************************
 * Routine: ScaleGridTo3D                               *
 *                                                      *
 * This is routine scales the 2D multiblock mesh to a   *
 * quasi-3D mesh.  Basically, it calculates scaling     *
 * parameters define a "thickness".  This is a          *
 * discusting hack - T-Bone doesn't like.               *
 *                                                      *
 ********************************************************/
extern void ScaleGridTo3D( Rte2D_Quad_Block *Soln_ptr,
			   AdaptiveBlock2D_List &Soln_Block_List ) {

  // loop along the blocks
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) 
      Soln_ptr[nb].ScaleGridTo3D();
  } /* endfor */
}
