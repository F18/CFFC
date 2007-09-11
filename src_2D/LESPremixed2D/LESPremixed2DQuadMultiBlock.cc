/****************** LESPremixed2DMultiBlock.cc **********************************
    Multi-Block Versions of Subroutines for 2D Multi-Block 
    Quadrilateral Mesh Solution Classes. 

    NOTES: - In IC's I currently only use one of the "5" 
             pass in funcitons
                        
************************************************************************/

/* Include 2D LES Premixed quadrilateral mesh solution header file. */

#ifndef _LESPREMIXED2D_QUAD_INCLUDED
#include "LESPremixed2DQuad.h"
#endif // _LESPREMIXED2D_QUAD_INCLUDED

/**************************************************************************
 * LESPremixed2D_Quad_Block -- Multiple Block External Subroutines.       *
 **************************************************************************/

/********************************************************
 * Routine: Allocate                                    *
 *                                                      *
 * Allocate memory for 1D array of 2D quadrilateral     *
 * multi-block solution blocks.                         *
 *                                                      *
 ********************************************************/
LESPremixed2D_Quad_Block* Allocate(LESPremixed2D_Quad_Block *Soln_ptr,
				   LESPremixed2D_Input_Parameters &Input_Parameters) {

    /* Allocate memory. */
    Soln_ptr = new LESPremixed2D_Quad_Block[Input_Parameters.Number_of_Blocks_Per_Processor];

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
LESPremixed2D_Quad_Block* Deallocate(LESPremixed2D_Quad_Block *Soln_ptr,
				     LESPremixed2D_Input_Parameters &Input_Parameters) {
 
    /* Deallocate memory. */
    for (int i = 0 ; i <= Input_Parameters.Number_of_Blocks_Per_Processor-1 ; ++i ) {
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
void ICs(LESPremixed2D_Quad_Block *Soln_ptr,
         AdaptiveBlock2D_List &Soln_Block_List, 
         LESPremixed2D_Input_Parameters &Input_Parameters) {

  LESPremixed2D_pState Wo[5];
  
  /* Define various reference flow states. */
  Wo[0] = Input_Parameters.Wo;
    
  //Assign initial data for each solution block.
  for (int i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
    if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
      // Set flow geometry indicator (planar/axisymmetric) flow block.
      Soln_ptr[i].Flow_Type = Input_Parameters.FlowType;
      Soln_ptr[i].Axisymmetric = Input_Parameters.Axisymmetric;
      //      Soln_ptr[i].Wall_Functions = Input_Parameters.Wall_Functions;
      Soln_ptr[i].Gravity = Input_Parameters.Gravity;
      Soln_ptr[i].debug_level = Input_Parameters.debug_level;
      Soln_ptr[i].Moving_wall_velocity = Input_Parameters.Moving_wall_velocity;

      // Set initial data. (LESPremixed2DQuadSingleBlock.cc)
      ICs(Soln_ptr[i], Input_Parameters.i_ICs, Wo, Input_Parameters, 
	  Soln_Block_List.Block[i].gblknum);
    }
  }  /* endfor */


  //cout<<"\n Initial conditons assigned \n"; cout.flush();
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
int Read_Restart_Solution(LESPremixed2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          LESPremixed2D_Input_Parameters &Input_Parameters,
                          int &Number_of_Time_Steps,
                          double &Time,
                          CPUTime &CPU_Time) {

    int i, i_new_time_set, nsteps;
    char prefix[256], extension[256], restart_file_name[256], line[256];
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
	  if(restart_file.fail()){ 
	    cerr<<"\nError opening file: "<<restart_file_name<<endl;
	    return 1; 
	  }	 
          
          // Read solution block data.
          restart_file.setf(ios::skipws);
          restart_file >> nsteps >> time0 >> cpu_time0;
          restart_file.unsetf(ios::skipws);

          // This is for the initialization of the 2D wrinkled premixed flame
          if (Input_Parameters.i_ICs == IC_2DWRINKLED_FLAME) {
            nsteps = 0; 
            time0 = ZERO; 
            cpu_time0.zero();
	  } 
    
	  /********** LESPREMIXED2D SPECIFIC ****************************************/
	  restart_file.getline(line,sizeof(line)); 
	  // get reaction set name
	  restart_file >>Input_Parameters.react_name;
	  if( Input_Parameters.react_name != "CANTERA")
	    Input_Parameters.Wo.React.set_reactions(Input_Parameters.react_name);
	  else
	    Input_Parameters.Wo.React.ct_load_mechanism(Input_Parameters.ct_mech_file, 
							Input_Parameters.ct_mech_name);

	  // multispecies but no reactions
	  if( Input_Parameters.react_name == "NO_REACTIONS"){
	    restart_file.setf(ios::skipws);
	    int num_species;
	    //get number of species
	    restart_file >> num_species; 
	    string *species = new string[num_species];
	    //get species names 
	    for(int k=0; k<num_species; k++){
	      restart_file >> species[k];
	    } 
	    restart_file.unsetf(ios::skipws);  
	    Input_Parameters.Wo.React.set_species(species,num_species);
	    delete[] species; 
	  }

	  //OVERRIDE 2-step CH4 mechanism name for 2D turbulent flame with mesh refinement
	  if (Input_Parameters.i_ICs == IC_2DWRINKLED_FLAME) {
            Input_Parameters.react_name = "CH4_2STEP";
	    Input_Parameters.Wo.React.set_reactions(Input_Parameters.react_name);
	    //cout << "\n" << Input_Parameters.react_name << flush;
	  }

	  
	  // scalar system
	  restart_file.getline(line,sizeof(line)); 
	  restart_file >> Input_Parameters.scalar_system_name;
	  //cout << "\n Scalar_system_name read: " << Input_Parameters.scalar_system_name << flush << endl;
          Input_Parameters.Wo.Scal_sys.scalar_set(Input_Parameters.scalar_system_name);
	  

	  //Set Data Path
	  Input_Parameters.get_cffc_path();

	  //setup properties 
	  Input_Parameters.Wo.set_species_data(Input_Parameters.Wo.React.num_species,
					       Input_Parameters.Wo.Scal_sys.num_scalars,
					       Input_Parameters.Wo.React.species,
					       Input_Parameters.CFFC_Path,
					       Input_Parameters.Mach_Number_Reference,
					       Input_Parameters.Schmidt,
	    			               Input_Parameters.i_trans_type);   
	  Input_Parameters.Uo.set_species_data(Input_Parameters.Wo.React.num_species,
					       Input_Parameters.Wo.Scal_sys.num_scalars,
					       Input_Parameters.Wo.React.species,
					       Input_Parameters.CFFC_Path,
					       Input_Parameters.Mach_Number_Reference,
					       Input_Parameters.Schmidt,
	  				       Input_Parameters.i_trans_type);    
          Input_Parameters.Wo.set_premixed_flame_variables(Input_Parameters.Wo.laminar_thickness,
				                           Input_Parameters.Wo.laminar_speed,
				                           Input_Parameters.Wo.TFactor,
                                                           Input_Parameters.Wo.adiabatic_temp,
                                                           Input_Parameters.Wo.equivalence_ratio,
                                                           Input_Parameters.Wo.reactants_den);
          Input_Parameters.Uo.set_premixed_flame_variables(Input_Parameters.Wo.laminar_thickness,
				                           Input_Parameters.Wo.laminar_speed,
				                           Input_Parameters.Wo.TFactor,
                                                           Input_Parameters.Wo.adiabatic_temp,
                                                           Input_Parameters.Wo.equivalence_ratio,
                                                           Input_Parameters.Wo.reactants_den);
	  Input_Parameters.Uo = U(Input_Parameters.Wo);

	  //Reset  Wo species mass fractions
          Input_Parameters.Wo.set_initial_values(Input_Parameters.mass_fractions);  
	  
          // Update Uo from Wo
          Input_Parameters.Uo = U(Input_Parameters.Wo);


          /********* Turbulent Premixed Reacting Flow ***********************************/
          if (Input_Parameters.FlowType == FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY ||
	      Input_Parameters.FlowType == FLOWTYPE_TURBULENT_LES_TF_K ||
              Input_Parameters.FlowType == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY ||
	      Input_Parameters.FlowType == FLOWTYPE_TURBULENT_LES_C_FSD_K) {
            restart_file.setf(ios::skipws); 
            restart_file >> Input_Parameters.Fresh_Fuel_Mass_Fraction
                         >> Input_Parameters.Burnt_Fuel_Mass_Fraction
                         >> Input_Parameters.Fresh_Density;
            restart_file.unsetf(ios::skipws);    
	  }  

	  if (Input_Parameters.i_Grid == GRID_PERIODIC_BOX ||
              Input_Parameters.i_Grid == GRID_2DTURBULENT_PREMIXED_FLAME) {
            restart_file.setf(ios::skipws);
	    restart_file >>  Input_Parameters.TKEo;
            restart_file.unsetf(ios::skipws);
	  }

	  /********** END LESPREMIXED2D SPECIFIC ****************************************/
	  
          if (!i_new_time_set) {
             Number_of_Time_Steps = nsteps;
             Input_Parameters.Maximum_Number_of_Time_Steps += Number_of_Time_Steps;
	     Time = time0;
             CPU_Time.cput = cpu_time0.cput;
	 
             i_new_time_set = 1;
          } /* endif */

	  /************** Quad Single Block Data from file *********************/
	  restart_file >> Soln_ptr[i];

	  /*********************************************************************/
	  //OVERIDE Restart File settings with Input File Settings (may have changed)
	  Soln_ptr[i].Flow_Type = Input_Parameters.FlowType;
	  Soln_ptr[i].Axisymmetric = Input_Parameters.Axisymmetric;
	  Soln_ptr[i].Gravity = Input_Parameters.Gravity;
	  Soln_ptr[i].debug_level = Input_Parameters.debug_level;
	  Soln_ptr[i].Moving_wall_velocity = Input_Parameters.Moving_wall_velocity;
	  /*********************************************************************/
	  	  
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
int Write_Restart_Solution(LESPremixed2D_Quad_Block *Soln_ptr,
                           AdaptiveBlock2D_List &Soln_Block_List,
                           LESPremixed2D_Input_Parameters &Input_Parameters,
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
	  /********* LESPREMIXED2D SPECIFIC ***********************************/
	  restart_file << Input_Parameters.react_name << "\n";
	  if(Input_Parameters.react_name == "NO_REACTIONS"){
	    restart_file << Input_Parameters.Wo.ns <<" ";
	    for(int k=0; k< Input_Parameters.Wo.ns; k++){ 
	      restart_file << Input_Parameters.multispecies[k] <<" ";
	    }
	    restart_file<<endl;
	  }


	  // scalar system name
	  restart_file << Input_Parameters.scalar_system_name << "\n";
 
          /********* Turbulent Premixed Reacting Flow ***********************************/
          if (Input_Parameters.FlowType == FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY ||
	     Input_Parameters.FlowType == FLOWTYPE_TURBULENT_LES_TF_K ||
              Input_Parameters.FlowType == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY ||
	      Input_Parameters.FlowType == FLOWTYPE_TURBULENT_LES_C_FSD_K) {
            restart_file << Input_Parameters.Fresh_Fuel_Mass_Fraction << " "
                         << Input_Parameters.Burnt_Fuel_Mass_Fraction << " "
                         << Input_Parameters.Fresh_Density << endl;
	  }

	  if (Input_Parameters.i_Grid == GRID_PERIODIC_BOX ||
              Input_Parameters.i_Grid == GRID_2DTURBULENT_PREMIXED_FLAME) {
	    restart_file <<  Input_Parameters.TKEo << endl;
	  }
 
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
int Output_Tecplot(LESPremixed2D_Quad_Block *Soln_ptr,
                   AdaptiveBlock2D_List &Soln_Block_List,
                   LESPremixed2D_Input_Parameters &Input_Parameters,
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
      
	Output_Tecplot(Soln_ptr[i], Input_Parameters,
		       Number_of_Time_Steps, 
		       Time,
		       Soln_Block_List.Block[i].gblknum,
		       i_output_title,
		       output_file);
	//cout<<" \n Returned to Multi \n";
          
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
int Output_Cells_Tecplot(LESPremixed2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List,
                         LESPremixed2D_Input_Parameters &Input_Parameters,
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
	 Output_Cells_Tecplot(Soln_ptr[i],Input_Parameters, 
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
int Output_Nodes_Tecplot(LESPremixed2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List,
                         LESPremixed2D_Input_Parameters &Input_Parameters,
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
int Output_Mesh_Tecplot(LESPremixed2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        LESPremixed2D_Input_Parameters &Input_Parameters,
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
 * Routine: Output_RHS                                  *
 * For testing Jacobians                                *
 *                                                      *
 *                                                      *
 ********************************************************/
int Output_RHS(LESPremixed2D_Quad_Block *Soln_ptr,
	       AdaptiveBlock2D_List &Soln_Block_List,
	       LESPremixed2D_Input_Parameters &Input_Parameters,
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
    strcat(prefix, "_RHS_cpu");

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
	 Output_RHS(Soln_ptr[i],
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
 * Routine: Output_PERTURB                              *
 * For testing Jacobians                                *
 *                                                      *
 *                                                      *
 ********************************************************/
int Output_PERTURB(LESPremixed2D_Quad_Block *Soln_ptr,
		   AdaptiveBlock2D_List &Soln_Block_List,
		   LESPremixed2D_Input_Parameters &Input_Parameters,
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
	  /********* LESPREMIXED2D SPECIFIC ***********************************/
	  restart_file << Input_Parameters.react_name << "\n";
	  if(Input_Parameters.react_name == "NO_REACTIONS"){
	    restart_file << Input_Parameters.Wo.ns <<" ";
	    for(int k=0; k< Input_Parameters.Wo.ns; k++){ 
	      restart_file << Input_Parameters.multispecies[k] <<" ";
	    }
	    restart_file<<endl;
	  }
	  /**************************************************************/

	  for(int jj= Soln_ptr[i].JCl; jj<=Soln_ptr[i].JCu; jj++)
	    for  (int ii= Soln_ptr[i].ICl; ii<=Soln_ptr[i].ICu ; ii++)
	      {	
		if((ii==16)&&(jj==4)){
		  Soln_ptr[i].U[ii][jj][7] =Soln_ptr[i].U[ii][jj][7]*(1.0+1e-8);
		  
		}
	      }
	  
	  
      	  restart_file << setprecision(14) << Soln_ptr[i];

    
          // Close restart file.
          restart_file.close();
       } /* endif */
    }  /* endfor */
    /* Writing of restart files complete.  Return zero value. */
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
int Output_Mesh_Gnuplot(LESPremixed2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        LESPremixed2D_Input_Parameters &Input_Parameters,
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
void BCs(LESPremixed2D_Quad_Block *Soln_ptr,
         AdaptiveBlock2D_List &Soln_Block_List,
	 LESPremixed2D_Input_Parameters &Input_Parameters) {

    int i;

    /* Prescribe boundary data for each solution block. */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
	 BCs(Soln_ptr[i],Input_Parameters);
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
 * condition.                                           *
 *                                                      *
 ********************************************************/
double CFL(LESPremixed2D_Quad_Block *Soln_ptr,
           AdaptiveBlock2D_List &Soln_Block_List,
           LESPremixed2D_Input_Parameters &Input_Parameters) {

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
void Set_Global_TimeStep(LESPremixed2D_Quad_Block *Soln_ptr,
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
double L1_Norm_Residual(LESPremixed2D_Quad_Block *Soln_ptr,
		      AdaptiveBlock2D_List &Soln_Block_List) {

  /* Calculate the L1-norm. Sum the L1-norm for each solution block. */
  double l1_norm(ZERO);
  for ( int i = 0 ; i < Soln_Block_List.Nblk; ++i ) {
    if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {     
      l1_norm += L1_Norm_Residual(Soln_ptr[i],Soln_ptr[i].residual_variable);       
    }
  }    

  return l1_norm;
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
double L2_Norm_Residual(LESPremixed2D_Quad_Block *Soln_ptr,
			AdaptiveBlock2D_List &Soln_Block_List) {
	
  /* Sum the square of the L2-norm for each solution block. */
  double l2_norm(ZERO);
  for ( int i = 0 ; i < Soln_Block_List.Nblk; ++i ) {
    if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
      l2_norm += sqr(L2_Norm_Residual(Soln_ptr[i],Soln_ptr[i].residual_variable));       
    } 
  }  
  return sqrt(l2_norm);
  

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
double Max_Norm_Residual(LESPremixed2D_Quad_Block *Soln_ptr,
			 AdaptiveBlock2D_List &Soln_Block_List){	
  
  /* Find the maximum norm for all solution blocks. */
  double max_norm(ZERO);
  for (int i = 0 ; i < Soln_Block_List.Nblk ; ++i ) {
    if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
      max_norm = max(max_norm, Max_Norm_Residual(Soln_ptr[i],Soln_ptr[i].residual_variable));
    } 
  } 

  return max_norm;
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
void L1_Norm_Residual(LESPremixed2D_Quad_Block *Soln_ptr,
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
void L2_Norm_Residual(LESPremixed2D_Quad_Block *Soln_ptr,
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
void Max_Norm_Residual(LESPremixed2D_Quad_Block *Soln_ptr,
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
void Evaluate_Limiters(LESPremixed2D_Quad_Block *Soln_ptr,
                       AdaptiveBlock2D_List &Soln_Block_List) {

    for (int i = 0 ; i < Soln_Block_List.Nblk ; ++i ) {
      if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) Soln_ptr[i].evaluate_limiters();
    }
}

/********************************************************
 * Routine: Freeze_Limiters                             *
 *                                                      *
 * Set conditions to freeze the limiters for a          *
 * 1D array of 2D quadrilateral multi-block solution    *
 * blocks.                                              *
 *                                                      *
 ********************************************************/
void Freeze_Limiters(LESPremixed2D_Quad_Block *Soln_ptr,
                     AdaptiveBlock2D_List &Soln_Block_List) {

    for (int i = 0 ; i < Soln_Block_List.Nblk ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) Soln_ptr[i].freeze_limiters();
    }  
}

/********************************************************
 * Routine: Change Mref                                 *
 *                                                      *
 ********************************************************/
void Change_Mref(LESPremixed2D_Quad_Block *Soln_ptr,
		 AdaptiveBlock2D_List &Soln_Block_List,
		 const double &Mr) {

    for ( int i = 0 ; i < Soln_Block_List.Nblk ; ++i ) {
      if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED){
	Soln_ptr[i].W[0][0].Mref = Mr;  //static so only need to change once per block.
	Soln_ptr[i].U[0][0].Mref = Mr;
	Soln_ptr[i].Uo[0][0].Mref = Mr;
      }
    }  
}

/********************************************************
 * Routine:  Residual_Smoothing                         *
 *                                                      *
 * Applies implicit residual smoothing to a 1D array of *
 * 2D quadrilateral multi-block  solution blocks.       *
 *                                                      *
 ********************************************************/
void Residual_Smoothing(LESPremixed2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        LESPremixed2D_Input_Parameters &Input_Parameters,
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

/********************************************************
 * Routine: Distance_to_Wall                            *
 *                                                      *
 * Determines the distance to the nearest wall for      *
 * given cell centered location.                        *
 *                                                      *
 ********************************************************/
double Distance_to_Wall(LESPremixed2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        const Vector2D X_cell) {

    int i;
    double y_wall;

    y_wall = 1e70;

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
	 y_wall = min(y_wall, Distance_to_Wall(Soln_ptr[i], X_cell));
       } /* endif */
    }  /* endfor */

    return (y_wall);

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
void Apply_Boundary_Flux_Corrections(LESPremixed2D_Quad_Block *Soln_ptr,
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
void Apply_Boundary_Flux_Corrections_Multistage_Explicit(LESPremixed2D_Quad_Block *Soln_ptr,
                                                         AdaptiveBlock2D_List &Soln_Block_List,
                                                         LESPremixed2D_Input_Parameters &Input_Parameters,
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
int dUdt_Multistage_Explicit(LESPremixed2D_Quad_Block *Soln_ptr,
			     AdaptiveBlockResourceList &Global_Soln_Block_List,
                             AdaptiveBlock2D_List &Local_Soln_Block_List,
                             LESPremixed2D_Input_Parameters &Input_Parameters,
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
int Update_Solution_Multistage_Explicit(LESPremixed2D_Quad_Block *Soln_ptr,
                                        AdaptiveBlock2D_List &Soln_Block_List,
                                        LESPremixed2D_Input_Parameters &Input_Parameters,
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
 * Routine: Update_Dual_Solution_States                 *
 *                                                      *
 * This routine updates solution states of the given    *
 * solution block corresponding to different times,     *
 * required in the dual time stepping.                  * 
 *                                                      *
 ********************************************************/
int Update_Dual_Solution_States(LESPremixed2D_Quad_Block *Soln_ptr,
                                 AdaptiveBlock2D_List &Soln_Block_List) {
    int i, error_flag;
 
    error_flag = 0;

    /* Update the solution states required by the dual time stepping for each * 
     *  solution block.                                                       */

    for ( i = 0 ; i <= Soln_Block_List.Nblk-1 ; ++i ) {
       if (Soln_Block_List.Block[i].used == ADAPTIVEBLOCK2D_USED) {
          error_flag =  Update_Dual_Solution_States(Soln_ptr[i]);
          if (error_flag) return (error_flag);
       }
    }  
    
    return(error_flag); 
}


//FROM DUSTY2D

/**********************************************************************
 * Routine: Output_Ringleb                                            *
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
int Output_Ringleb(LESPremixed2D_Quad_Block *Soln_ptr,
		   AdaptiveBlock2D_List &Soln_Block_List,
		   LESPremixed2D_Input_Parameters &IP) {

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
  
  // Write the solution data for each solution block.
  i_output_title = 1;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      Output_Ringleb_Solution(Soln_ptr[nb],nb,i_output_title,output_file);
      if (i_output_title) i_output_title = 0;
    }
  }
  
  // Close the output data file.
  output_file.close();

  // Now determine and output the L1, L2, and Max norms of the error and
  // the number of active cells.
  int numberofactivecells = 0;
  double l1_norm = ZERO; 
  double l2_norm = ZERO;
  double max_norm = ZERO;
  int numberofactivecells_temp = 0;
  double l1_temp = ZERO;
  double l2_temp = ZERO;
  double max_temp = ZERO;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO; numberofactivecells_temp = 0;
      Output_Ringleb_Error(Soln_ptr[nb],l1_temp,l2_temp,max_temp,numberofactivecells_temp);
      numberofactivecells += numberofactivecells_temp;
      l1_norm += l1_temp;
      l2_norm += sqr(l2_temp);
      max_norm = max(max_norm,max_temp);
    }
  }

#ifdef _MPI_VERSION
  l1_norm = CFFC_Summation_MPI(l1_norm);
  l2_norm = CFFC_Summation_MPI(l2_norm);
  max_norm = CFFC_Maximum_MPI(max_norm);
  numberofactivecells = CFFC_Summation_MPI(numberofactivecells);
#endif

  // Calculate the L1-norm and L2-norm for all blocks.
  l1_norm = l1_norm/numberofactivecells;
  l2_norm = sqrt(l2_norm/numberofactivecells);

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
	 << "   Number of cells = " << numberofactivecells
	 << endl
	 << " ==================================================================== "
	 << endl;
  }

  // Writing of output data files complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Output_Viscous_Channel                                    *
 *                                                                    *
 * Writes the exact and computed viscous channel flow solution values *
 * at the nodes for a 1D array of 2D quadrilateral multi-block        *
 * solution blocks to the specified output data file(s) in a format   *
 * suitable for plotting with TECPLOT.  Returns a non-zero value if   *
 * cannot write any of the TECPLOT solution files.  This routine will *
 * also compute the error norms of the computed solution (L1, L2, and *
 * max norms).                                                        *
 *                                                                    *
 **********************************************************************/
int Output_Viscous_Channel(LESPremixed2D_Quad_Block *Soln_ptr,
			   AdaptiveBlock2D_List &Soln_Block_List,
			   LESPremixed2D_Input_Parameters &IP) {

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
  
  double dpdx = -3177.7;  //should be from dp/dx = Input_Parameter.delta_pres/length

  // Write the solution data for each solution block.
  i_output_title = 1;
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO;
      Output_Viscous_Channel(Soln_ptr[nb],nb,i_output_title,output_file,
			     l1_temp,l2_temp,max_temp,IP.Moving_wall_velocity,dpdx);
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
	 << " Error norms for the viscous channel flow:"
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

/**********************************************************************
 * Routine: Output_Flat_Plate                                         *
 *                                                                    *
 * This routine outputs the non-dimensionalized computed flat plate   *
 * solution and the corresponding Blasius solution.                   *
 *                                                                    *
 **********************************************************************/
int Output_Flat_Plate(LESPremixed2D_Quad_Block *Soln_ptr,
		      AdaptiveBlock2D_List &Soln_Block_List,
		      LESPremixed2D_Input_Parameters &IP) {

  int i, i_output_title_soln, i_output_title_skin;
  char prefix_soln[256], prefix_skin[256], extension_soln[256], extension_skin[256];
  char output_file_name_soln[256], output_file_name_skin[256];
  char *output_file_name_soln_ptr, *output_file_name_skin_ptr;
  ofstream output_file_soln, output_file_skin;
  double l1_norm, l2_norm, max_norm;
  double l1_temp, l2_temp, max_temp;
  
  // Initialize error variables.
  l1_norm = ZERO; l2_norm = ZERO; max_norm = ZERO;
  l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO; 

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
  for (int nb = 0; nb < Soln_Block_List.Nblk; nb++) {
    if (Soln_Block_List.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
      l1_temp = ZERO; l2_temp = ZERO; max_temp = ZERO; 
      Output_Flat_Plate(Soln_ptr[nb],nb,
			i_output_title_soln,output_file_soln,
			i_output_title_skin,output_file_skin,
			IP.Wo,
			l1_temp,l2_temp,max_temp);
      l1_norm += l1_temp;
      l2_norm += sqr(l2_temp);
      max_norm = max(max_norm,max_temp);
      if (i_output_title_soln) i_output_title_soln = 0;
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
#endif
 
  l2_norm = sqrt(l2_norm);

  if (CFFC_Primary_MPI_Processor()) {
    cout << endl
	 << endl
	 << " ==================================================================== "
	 << endl
	 << " Error norms for the skin friction coefficient for the laminar flow"
	 << " over a flat plate:"
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
 * Routine: Output_Driven_Cavity_Flow                                 *
 *                                                                    *
 * This routine outputs a comparison of the computed solution for the *
 * driven cavity flow with the computations done by Ghia et al. (J.   *
 * Comp. Phys. Vol. 48 1982) for a 1D array of 2D quadrilateral       *
 * multi-block solution blocks to the specified output data file(s)   *
 * in a format suitable for plotting with tecplot.                    *
 *                                                                    *
 **********************************************************************/
int Output_Driven_Cavity_Flow(LESPremixed2D_Quad_Block *Soln_ptr,
			      AdaptiveBlock2D_List &Soln_Block_List,
			      LESPremixed2D_Input_Parameters &IP) {

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
      Output_Driven_Cavity_Flow(Soln_ptr[nb],nb,
				i_output_title,
				output_file_u,
				output_file_v,
				IP.Re_lid,
				IP.Moving_wall_velocity,
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
