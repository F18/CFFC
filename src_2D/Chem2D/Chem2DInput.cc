/****************** Chem2DInput.cc ************************************
  Constructors for the Chem2DInput class.

  associated header file:  Chem2DInput.h
***********************************************************************/

// CHEM2D Header file
#ifndef _CHEM2D_INPUT_INCLUDED
#include "Chem2DInput.h"
#endif // _CHEM2D_INPUT_INCLUDED

/*************************************************************
 * Chem2D_Input_Parameters -- External subroutines.         *
 *************************************************************/

/********************************************************
 * Routine: Open_Input_File                             *
 *                                                      *
 * Opens the appropriate input data file.               *
 *                                                      *
 ********************************************************/
void Open_Input_File(Chem2D_Input_Parameters &IP) {

    IP.Input_File.open(IP.Input_File_Name, ios::in);
    if (!IP.Input_File.bad()) {
       IP.Line_Number = 0;
       IP.Input_File.setf(ios::skipws);
    } /* endif */

}

/********************************************************
 * Routine: Close_Input_File                            *
 *                                                      *
 * Closes the appropriate input data file.              *
 *                                                      *
 ********************************************************/
void Close_Input_File(Chem2D_Input_Parameters &IP) {

    IP.Input_File.unsetf(ios::skipws);
    IP.Input_File.close();

}

/********************************************************
 * Routine: Set_Default_Input_Parameters                *
 *                                                      *
 * Assigns default values to the input parameters.      *
 *                                                      *
 ********************************************************/
void Set_Default_Input_Parameters(Chem2D_Input_Parameters &IP) {

    int i;
    char *string_ptr;

    string_ptr = "Chem2D.in";
    strcpy(IP.Input_File_Name, string_ptr);

    // Time-stepping parameters:
    string_ptr = "Explicit_Euler";
    strcpy(IP.Time_Integration_Type, string_ptr);
    IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
    IP.Time_Accurate = 0;
    IP.Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
    IP.Preconditioning = 0; //default off
    IP.Maximum_Number_of_Time_Steps = 100;
    IP.N_Stage = 1;
    IP.CFL_Number = 0.5;
    IP.Time_Max = ZERO;

    // Residual variable:
    IP.i_Residual_Variable = 1;

    // Residual smoothing:
    IP.Residual_Smoothing = 0;
    IP.Residual_Smoothing_Epsilon = ZERO;
    IP.Residual_Smoothing_Gauss_Seidel_Iterations = 2;

    // Reconstruction type:
    string_ptr = "Least_Squares";
    strcpy(IP.Reconstruction_Type, string_ptr);
    IP.i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;

    // Limiter type:
    string_ptr = "Barth_Jespersen";
    strcpy(IP.Limiter_Type, string_ptr);
    IP.i_Limiter = LIMITER_BARTH_JESPERSEN;

    // Flux function:
    string_ptr = "ROE";
    strcpy(IP.Flux_Function_Type, string_ptr);
    IP.i_Flux_Function = FLUX_FUNCTION_ROE;
 
    // Viscous gradient reconstruction type:
    string_ptr = "Arithmetic";
    strcpy(IP.Viscous_Flux_Evaluation_Type, string_ptr);
    IP.i_Viscous_Flux_Evaluation = VISCOUS_RECONSTRUCTION_DIAMOND_PATH;

    // Initial conditions:
    string_ptr = "Uniform";
    strcpy(IP.ICs_Type, string_ptr);
    IP.i_ICs = IC_UNIFORM;

    /* Directory Path */
    IP.get_cfdkit_path();

    /* Debug Level */
    IP.debug_level = 0;  //default no debug information

    IP.Mach_Number = ZERO;
    IP.Mach_Number_Reference = ONE;  //for precondtioning
    IP.Flow_Angle = ZERO; 
    IP.Global_Schmidt = ONE;

    /******************************************/
    /********** CHEM2D SPECIFIC ***************/
    //define multispecies with no reactions.
    IP.react_name ="NO_REACTIONS";   
    //strcpy(IP.React_Name,"NO_REACTIONS"); 
    //Use air with 79% N2, and 21% 02 by volume.(ie. mol)
    IP.num_species = 2;
    IP.Allocate();
    IP.multispecies[0] = "N2"; 
    IP.multispecies[1] = "O2"; 
    IP.mass_fractions[0] = 0.765; 
    IP.mass_fractions[1] = 0.235;
    IP.Schmidt[0] = IP.Global_Schmidt;
    IP.Schmidt[1] = IP.Global_Schmidt;

    IP.Wo.React.set_reactions(IP.react_name); 
    IP.Wo.React.set_species(IP.multispecies,IP.num_species);
   
    //Get Species parameters and set default initial values
    IP.Wo.set_species_data(IP.num_species,IP.multispecies,IP.CFDkit_Path,
			   IP.debug_level,IP.Mach_Number_Reference,IP.Schmidt); 
    IP.Uo.set_species_data(IP.num_species,IP.multispecies,IP.CFDkit_Path,
			   IP.debug_level,IP.Mach_Number_Reference,IP.Schmidt);

    //Air at STD_ATM
    IP.Pressure = IP.Wo.p;
    IP.Temperature = IP.Wo.T(); 
    IP.Wo.v.x =IP.Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Flow_Angle/360.00);
    IP.Wo.v.y =IP.Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Flow_Angle/360.00);
    IP.Wo.set_initial_values(IP.mass_fractions);
    IP.Uo.set_initial_values(IP.mass_fractions);

    IP.Uo = U(IP.Wo);   
    /***** END CHEM2D SPECFIC *****************/
    /******************************************/
    //BC

    IP.Moving_wall_velocity = ZERO; 

    /* Flow type */
    string_ptr = "Inviscid";
    strcpy(IP.Flow_Type, string_ptr);
    IP.FlowType = FLOWTYPE_LAMINAR;
    IP.Wo.set_flow_type(IP.FlowType); 
    IP.Uo.set_flow_type(IP.FlowType); 

    /* Flow geometry type */
    string_ptr = "Planar";
    strcpy(IP.Flow_Geometry_Type, string_ptr);
    IP.Axisymmetric = 0;
    IP.Gravity = 0;  //default sans gravity

    // Turbulence parameters:
    string_ptr = "Direct_Integration";
    strcpy(IP.Turbulence_BC_Type,string_ptr);
    IP.i_Turbulence_BCs = TURBULENT_BC_DIRECT_INTEGRATION;
    string_ptr = "Wall";
    strcpy(IP.Friction_Velocity_Type,string_ptr);
    IP.i_Friction_Velocity = FRICTION_VELOCITY_WALL_SHEAR_STRESS;
    IP.C_constant = 5.0;
    IP.von_Karman_Constant = 0.41;
    IP.yplus_sublayer = 5.0;
    IP.yplus_buffer_layer = 30.0;
    IP.yplus_outer_layer = 100.0;
    IP.Wo.set_turbulence_variables(IP.C_constant,
				   IP.von_Karman_Constant,
				   IP.yplus_sublayer,
				   IP.yplus_buffer_layer,
				   IP.yplus_outer_layer); 
    IP.Uo.set_turbulence_variables(IP.C_constant,
				   IP.von_Karman_Constant,
				   IP.yplus_sublayer,
				   IP.yplus_buffer_layer,
				   IP.yplus_outer_layer);

    /* Grid Parameters */
    string_ptr = "Square";
    strcpy(IP.Grid_Type, string_ptr);
    IP.i_Grid = GRID_SQUARE;
    IP.Box_Width = ONE;
    IP.Box_Height = ONE;
    IP.Number_of_Cells_Idir = 100;
    IP.Number_of_Cells_Jdir = 100;
    IP.Number_of_Ghost_Cells = 2;
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Blocks_Jdir = 1;
 
    IP.Plate_Length = ONE;
    IP.Pipe_Length = ONE;
    IP.Pipe_Radius = 0.1;
    IP.Blunt_Body_Radius = ONE;
    IP.Blunt_Body_Mach_Number = TWO;
    IP.Chamber_Length = 0.835;
    IP.Chamber_Radius = 0.020;
    IP.Chamber_To_Throat_Length = 0.05;
    IP.Nozzle_Length = 0.150;
    IP.Nozzle_Radius_Exit = 0.030;
    IP.Nozzle_Radius_Throat = 0.010;
    IP.Nozzle_Type = NOZZLE_GOTTLIEB_FUNCTION;
    IP.Grain_Radius = 0.0;
    IP.Cylinder_Radius = ONE;
    IP.Ellipse_Length_X_Axis = TWO;
    IP.Ellipse_Length_Y_Axis = HALF;
    IP.Chord_Length = ONE;
    IP.Orifice_Radius = ONE;
    IP.Wedge_Angle = 25.0;
    IP.Wedge_Length = HALF;
    IP.Smooth_Bump = OFF;

    IP.X_Shift = Vector2D_ZERO;
    IP.X_Scale = ONE;
    IP.X_Rotate = ZERO;

    // Boundary conditions:
    string_ptr = "Off";
    strcpy(IP.Boundary_Conditions_Specified,string_ptr);
    IP.BCs_Specified = OFF;
    string_ptr = "None";
    strcpy(IP.BC_North_Type,string_ptr);
    strcpy(IP.BC_South_Type,string_ptr);
    strcpy(IP.BC_East_Type,string_ptr);
    strcpy(IP.BC_West_Type,string_ptr);
    IP.BC_North = BC_NONE;
    IP.BC_South = BC_NONE;
    IP.BC_East  = BC_NONE;
    IP.BC_West  = BC_NONE;

    // Mesh stretching factor:
    IP.i_Mesh_Stretching = OFF;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_LINEAR;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_LINEAR;
    IP.Mesh_Stretching_Factor_Idir = 1.01;
    IP.Mesh_Stretching_Factor_Jdir = 1.01;

    // ICEM:
    IP.ICEMCFD_FileNames = ICEMCFD_get_filenames();

    // AMR:
    IP.AMR = 0;
    IP.AMR_Frequency = 100;
    IP.Number_of_Initial_Mesh_Refinements = 0;
    IP.Number_of_Uniform_Mesh_Refinements = 0;
    IP.Number_of_Boundary_Mesh_Refinements = 0;
    IP.Maximum_Refinement_Level = 100;
    IP.Minimum_Refinement_Level = 1;
    IP.Threshold_for_Refinement = 0.50;
    IP.Threshold_for_Coarsening = 0.10;
    IP.Number_of_Refinement_Criteria = 6;
    IP.Refinement_Criteria_Gradient_Density = ON;
    IP.Refinement_Criteria_Divergence_Velocity = ON;
    IP.Refinement_Criteria_Curl_Velocity = ON;
    IP.Refinement_Criteria_Gradient_Temperature = ON;
    IP.Refinement_Criteria_Gradient_CH4 = ON;
    IP.Refinement_Criteria_Gradient_CO2 = ON;

    // Smooth quad block flag:
    IP.i_Smooth_Quad_Block = ON;

    string_ptr = "outputfile.dat";
    strcpy(IP.Output_File_Name, string_ptr);

    string_ptr = "gridfile.grid";
    strcpy(IP.Grid_File_Name, string_ptr);
    string_ptr = "gridfile.griddef";
    strcpy(IP.Grid_Definition_File_Name, string_ptr);

    string_ptr = "restartfile.soln";
    strcpy(IP.Restart_File_Name, string_ptr);

    string_ptr = "gnuplotfile.gplt";
    strcpy(IP.Gnuplot_File_Name, string_ptr);

    string_ptr = "Tecplot";
    strcpy(IP.Output_Format_Type, string_ptr);
    IP.i_Output_Format = IO_TECPLOT;
    IP.Restart_Solution_Save_Frequency = 1000;

    // Default output progress frequency:
    IP.Output_Progress_Frequency = 50;

    //CHEM2D
    IP.Time_Accurate_Plot_Frequency = 0;

    string_ptr = " ";
    strcpy(IP.Next_Control_Parameter, string_ptr);

    IP.Line_Number = 0;

    IP.Number_of_Processors = CFDkit_MPI::Number_of_Processors;
    IP.Number_of_Blocks_Per_Processor = 10;      

    // Freezing_Limiter
    IP.Freeze_Limiter = 0;

    // Limiter_switch
    IP.Freeze_Limiter_Residual_Level = 1e-4;
}

/********************************************************
 * Routine: Broadcast_Input_Parameters                  *
 *                                                      *
 * Broadcast the input parameters variables to all      *
 * processors involved in the calculation from the      *
 * primary processor using the MPI broadcast routine.   *
 *                                                      *
 ********************************************************/
void Broadcast_Input_Parameters(Chem2D_Input_Parameters &IP) {

#ifdef _MPI_VERSION

    MPI::COMM_WORLD.Bcast(IP.Input_File_Name, 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Line_Number), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Time_Integration_Type, 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Time_Integration), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Time_Accurate), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Local_Time_Stepping), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Maximum_Number_of_Time_Steps), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.N_Stage), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.CFL_Number), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Time_Max), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Residual_Variable),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Residual_Smoothing), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Residual_Smoothing_Epsilon), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Residual_Smoothing_Gauss_Seidel_Iterations), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Reconstruction_Type, 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Reconstruction), 
                          1, 
                          MPI::INT, 0);   
    MPI::COMM_WORLD.Bcast(IP.Viscous_Flux_Evaluation_Type, 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Viscous_Flux_Evaluation), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Limiter_Type, 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Limiter), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Flux_Function_Type, 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Flux_Function), 
                          1, 
                          MPI::INT, 0);
    // Turbulence parameters:
    MPI::COMM_WORLD.Bcast(IP.Turbulence_BC_Type,
			  INPUT_PARAMETER_LENGTH_CHEM2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Turbulence_BCs),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(IP.Friction_Velocity_Type,
			  INPUT_PARAMETER_LENGTH_CHEM2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Friction_Velocity),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.C_constant),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.von_Karman_Constant),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.yplus_sublayer),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.yplus_buffer_layer),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.yplus_outer_layer),
			  1,
			  MPI::DOUBLE,0);
    // Initial conditions:
    MPI::COMM_WORLD.Bcast(IP.ICs_Type,
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_ICs), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Pressure), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Temperature), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Mach_Number), 
			  1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Mach_Number_Reference), 
			  1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Flow_Angle), 
                          1, 
                          MPI::DOUBLE, 0);   

    /*********************************************************
     ******************* CHEM2D SPECIFIC *********************
     *********************************************************/
    MPI::COMM_WORLD.Bcast(IP.CFDkit_Path, 
			  INPUT_PARAMETER_LENGTH_CHEM2D, 
			  MPI::CHAR, 0);
    //reaction name
    MPI::COMM_WORLD.Bcast(IP.React_Name, 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
			  MPI::CHAR, 0);
    
    //delete current dynamic memory before changing num_species
    if(!CFDkit_Primary_MPI_Processor()) {   
      IP.Deallocate();
    } 
    //number of species
    MPI::COMM_WORLD.Bcast(&(IP.num_species), 
                          1, 
                          MPI::INT, 0);
    //set up new dynamic memory
    if(!CFDkit_Primary_MPI_Processor()) {   
      IP.Allocate();
    }

    //species names & mass fractions
    for(int i =0; i < IP.num_species; i++){
      MPI::COMM_WORLD.Bcast(&(IP.mass_fractions[i]), 
			    1, 
			    MPI::DOUBLE, 0);
      MPI::COMM_WORLD.Bcast(&(IP.Schmidt[i]), 
			    1, 
			    MPI::DOUBLE, 0);
      MPI::COMM_WORLD.Bcast(IP.Multispecies[i], 
			    INPUT_PARAMETER_LENGTH_CHEM2D, 
			    MPI::CHAR, 0);
    }
    //set recaction and species parameters
    if (!CFDkit_Primary_MPI_Processor()) {      
      IP.react_name = IP.React_Name;
      for (int i = 0; i < IP.num_species; i++) {
	IP.multispecies[i] = IP.Multispecies[i];  
      }    

      //load reaction names
      IP.Wo.React.set_reactions(IP.react_name);
      
      //Set species if non-reacting
      if( IP.Wo.React.reactset_flag == NO_REACTIONS){
	IP.Wo.React.set_species(IP.multispecies,IP.num_species);
      }  
          
      //set the data for each
      IP.get_cfdkit_path();
      IP.Wo.set_species_data(IP.num_species,IP.multispecies,IP.CFDkit_Path,
			     IP.debug_level,IP.Mach_Number_Reference,IP.Schmidt); 
      IP.Uo.set_species_data(IP.num_species,IP.multispecies,IP.CFDkit_Path,
			     IP.debug_level,IP.Mach_Number_Reference,IP.Schmidt);
      IP.Wo.set_initial_values(IP.mass_fractions);
      IP.Uo.set_initial_values(IP.mass_fractions);
   
      //set proper Temp & Pressure instead of defaults
      IP.Wo.rho = IP.Pressure/(IP.Wo.Rtot()*IP.Temperature); 
      IP.Wo.p = IP.Pressure;	
      IP.Wo.v.zero();

      IP.Uo = U(IP.Wo);
    } 
    /*********************************************************
     ******************* CHEM2D END **************************
     *********************************************************/

   if(!CFDkit_Primary_MPI_Processor()) {   
      IP.Wo.v.x = IP.Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Flow_Angle/360.00);
      IP.Wo.v.y = IP.Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Flow_Angle/360.00);
    }

    /***************************************/
    MPI::COMM_WORLD.Bcast(IP.Flow_Type, 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.FlowType), 
                          1, 
                          MPI::INT, 0);

   /***************************************/
    MPI::COMM_WORLD.Bcast(IP.Flow_Geometry_Type, 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Axisymmetric), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Global_Schmidt), 
                          1, 
                          MPI::INT, 0);    
    MPI::COMM_WORLD.Bcast(&(IP.Gravity), 
                          1, 
			  MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.debug_level), 
                          1, 
			  MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Preconditioning), 
                          1, 
			  MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Moving_wall_velocity), 
                          1, 
			  MPI::DOUBLE, 0);

    /**************************************/
    MPI::COMM_WORLD.Bcast(IP.Grid_Type, 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.NACA_Aerofoil_Type, 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Grid), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Cells_Idir), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Cells_Jdir),
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Ghost_Cells),
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Idir), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Jdir), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Box_Width), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Box_Height), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Plate_Length), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Pipe_Length), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Pipe_Radius), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Blunt_Body_Radius), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Blunt_Body_Mach_Number), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Chamber_Length),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Chamber_Radius),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Chamber_To_Throat_Length),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Nozzle_Length),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Nozzle_Radius_Exit),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Nozzle_Radius_Throat),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Nozzle_Type),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Grain_Radius),
			  1,
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Cylinder_Radius), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Ellipse_Length_X_Axis), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Ellipse_Length_Y_Axis), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Chord_Length), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Orifice_Radius), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Wedge_Angle), 
			  1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Wedge_Length), 
			  1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Smooth_Bump),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.X_Shift.x), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.X_Shift.y), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.X_Scale), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.X_Rotate), 
                          1, 
                          MPI::DOUBLE, 0);
    // Boundary Conditions:
    MPI::COMM_WORLD.Bcast(IP.Boundary_Conditions_Specified,
			  INPUT_PARAMETER_LENGTH_CHEM2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(&(IP.BCs_Specified),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(IP.BC_North_Type,
			  INPUT_PARAMETER_LENGTH_CHEM2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(IP.BC_South_Type,
			  INPUT_PARAMETER_LENGTH_CHEM2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(IP.BC_East_Type,
			  INPUT_PARAMETER_LENGTH_CHEM2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(IP.BC_West_Type,
			  INPUT_PARAMETER_LENGTH_CHEM2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(&(IP.BC_North),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.BC_South),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.BC_East),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.BC_West),
			  1,
			  MPI::INT,0);
    // ICEM:
    if (!CFDkit_Primary_MPI_Processor()) {
      IP.ICEMCFD_FileNames = new char*[3];
       for (int i = 0; i < 3; i++) {
          IP.ICEMCFD_FileNames[i] = new char[INPUT_PARAMETER_LENGTH_CHEM2D];
       } /* endfor */
    } /* endif */
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[0], 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[1], 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[2], 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
                          MPI::CHAR, 0);

    // AMR & Refinement Parameters
    MPI::COMM_WORLD.Bcast(&(IP.AMR), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.AMR_Frequency),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Initial_Mesh_Refinements), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Uniform_Mesh_Refinements),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Boundary_Mesh_Refinements),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Maximum_Refinement_Level),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Minimum_Refinement_Level),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Threshold_for_Refinement), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Threshold_for_Coarsening), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Refinement_Criteria),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Refinement_Criteria_Gradient_Density),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Refinement_Criteria_Divergence_Velocity),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Refinement_Criteria_Curl_Velocity),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Refinement_Criteria_Gradient_Temperature),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Refinement_Criteria_Gradient_CH4),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Refinement_Criteria_Gradient_CO2),
			  1,
			  MPI::DOUBLE,0);
    // Mesh stretching flag:
    MPI::COMM_WORLD.Bcast(&(IP.i_Mesh_Stretching),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Mesh_Stretching_Type_Idir),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Mesh_Stretching_Type_Jdir),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Mesh_Stretching_Factor_Idir),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Mesh_Stretching_Factor_Jdir),
			  1,
			  MPI::DOUBLE,0);
    // Smooth quad block flag:
    MPI::COMM_WORLD.Bcast(&(IP.i_Smooth_Quad_Block),
			  1,
			  MPI::INT,0);

    // File Names
    MPI::COMM_WORLD.Bcast(IP.Output_File_Name, 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Grid_File_Name, 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Grid_Definition_File_Name, 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Restart_File_Name, 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Gnuplot_File_Name, 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Output_Format_Type, 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Output_Format), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Restart_Solution_Save_Frequency), 
                          1, 
                          MPI::INT, 0); 

    // Output progress frequency:
    MPI::COMM_WORLD.Bcast(&(IP.Output_Progress_Frequency),
			  1,
			  MPI::INT,0);
    // Multigrid Related Parameters
    IP.Multigrid_IP.Broadcast_Input_Parameters();
    
    if (!CFDkit_Primary_MPI_Processor()) {
       IP.Number_of_Processors = CFDkit_MPI::Number_of_Processors;
    } /* endif */
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Per_Processor), 
                          1, 
                          MPI::INT, 0); 
    
    // Freeze_Limiter
    MPI::COMM_WORLD.Bcast(&(IP.Freeze_Limiter), 
			  1, 
			  MPI::INT, 0);

    // Freeze_Limiter_Residual_Level
    MPI::COMM_WORLD.Bcast(&(IP.Freeze_Limiter_Residual_Level),
                          1, 
                          MPI::DOUBLE, 0);

    // Reset the static variables.
    IP.Wo.set_flow_type(IP.FlowType);
    IP.Uo.set_flow_type(IP.FlowType);
    IP.Wo.set_turbulence_variables(IP.C_constant,
				   IP.von_Karman_Constant,
				   IP.yplus_sublayer,
				   IP.yplus_buffer_layer,
				   IP.yplus_outer_layer); 
    IP.Uo.set_turbulence_variables(IP.C_constant,
				   IP.von_Karman_Constant,
				   IP.yplus_sublayer,
				   IP.yplus_buffer_layer,
				   IP.yplus_outer_layer);

#endif
}

#ifdef _MPI_VERSION
/********************************************************
 * Routine: Broadcast_Input_Parameters                  *
 *                                                      *
 * Broadcast the input parameters variables to all      *
 * processors associated with the specified communicator*
 * from the specified processor using the MPI broadcast *
 *routine.                                              *
 *                                                      *
 ********************************************************/
void Broadcast_Input_Parameters(Chem2D_Input_Parameters &IP,
                                MPI::Intracomm &Communicator, 
                                const int Source_CPU) {
 
    int Source_Rank = 0;
    int i;
   
    Communicator.Bcast(IP.Input_File_Name, 
                       INPUT_PARAMETER_LENGTH_CHEM2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.Line_Number), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Time_Integration_Type, 
                       INPUT_PARAMETER_LENGTH_CHEM2D, 
                       MPI::CHAR, Source_Rank);
    // Time integration:
    Communicator.Bcast(&(IP.i_Time_Integration), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Time_Accurate), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Local_Time_Stepping), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Maximum_Number_of_Time_Steps), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.N_Stage), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.CFL_Number), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Time_Max), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    // Residual variable:
    Communicator.Bcast(&(IP.i_Residual_Variable),
		       1,
		       MPI::DOUBLE,Source_Rank);
    // Residual Smoothing:
    Communicator.Bcast(&(IP.Residual_Smoothing), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Residual_Smoothing_Epsilon), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Residual_Smoothing_Gauss_Seidel_Iterations), 
                       1, 
                       MPI::INT, Source_Rank);
    // Reconstruction:
    Communicator.Bcast(IP.Reconstruction_Type, 
                       INPUT_PARAMETER_LENGTH_CHEM2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Reconstruction), 
                       1, 
                       MPI::INT, Source_Rank);
    // Flux functions:
    Communicator.Bcast(IP.Viscous_Flux_Evaluation_Type, 
                       INPUT_PARAMETER_LENGTH_CHEM2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Viscous_Flux_Evaluation), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Limiter_Type, 
                       INPUT_PARAMETER_LENGTH_CHEM2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Limiter), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Flux_Function_Type, 
                       INPUT_PARAMETER_LENGTH_CHEM2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Flux_Function), 
                       1, 
                       MPI::INT, Source_Rank);
    // Turbulence parameters:
    Communicator.Bcast(IP.Turbulence_BC_Type,
		       INPUT_PARAMETER_LENGTH_CHEM2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(&(IP.i_Turbulence_BCs),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(IP.Friction_Velocity_Type,
		       INPUT_PARAMETER_LENGTH_CHEM2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(&(IP.i_Friction_Velocity),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.C_constant),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.von_Karman_Constant),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.yplus_sublayer),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.yplus_buffer_layer),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.yplus_outer_layer),
		       1,
		       MPI::DOUBLE,Source_Rank);
    // Initial conditions:
    Communicator.Bcast(IP.ICs_Type, 
                       INPUT_PARAMETER_LENGTH_CHEM2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_ICs), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Pressure), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Temperature), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Mach_Number), 
                       1, 
                       MPI::DOUBLE, Source_Rank); 
    Communicator.Bcast(&(IP.Mach_Number_Reference), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Flow_Angle), 
                       1, 
                       MPI::DOUBLE, Source_Rank);

    /*********************************************************
     ******************* CHEM2D SPECIFIC *********************
     *********************************************************/
    Communicator.Bcast(IP.CFDkit_Path, 
			  INPUT_PARAMETER_LENGTH_CHEM2D, 
			  MPI::CHAR, Source_Rank);
    //reaction name
    Communicator.Bcast(IP.React_Name, 
                          INPUT_PARAMETER_LENGTH_CHEM2D, 
			  MPI::CHAR, Source_Rank);
    //delete orginal dynamic memory before changing num_species
    if(!CFDkit_Primary_MPI_Processor()) {  
      IP.Deallocate();  
    } 
    //number of species
    Communicator.Bcast(&(IP.num_species), 
                          1, 
                          MPI::INT,Source_Rank );
    //set up dynamic memory
    if(!CFDkit_Primary_MPI_Processor()) {  
      IP.Allocate();
    } 
    //species names & mass fractions
    for(int i =0; i < IP.num_species; i++){
      Communicator.Bcast(&(IP.mass_fractions[i]), 
			    1, 
			    MPI::DOUBLE, Source_Rank);
      Communicator.Bcast(&(IP.Schmidt[i]), 
			 1, 
			 MPI::DOUBLE, Source_Rank);
      Communicator.Bcast(IP.Multispecies[i], 
			 INPUT_PARAMETER_LENGTH_CHEM2D, 
			 MPI::CHAR, Source_Rank);
    }
    //set recaction and species parameters
    if (!CFDkit_Primary_MPI_Processor()) {      
      IP.react_name = IP.React_Name;
      for (int i = 0; i < IP.num_species; i++) {
	IP.multispecies[i] = IP.Multispecies[i];  
      }     
      //load reaction names
      IP.Wo.React.set_reactions(IP.react_name);
      
      //Set species if non-reacting
      if( IP.Wo.React.reactset_flag == NO_REACTIONS){
	IP.Wo.React.set_species(IP.multispecies,IP.num_species);
      }  

      //set the data for each
      IP.get_cfdkit_path();
      IP.Wo.set_species_data(IP.num_species,IP.multispecies,IP.CFDkit_Path,
			     IP.debug_level,IP.Mach_Number_Reference,IP.Schmidt); 
      IP.Uo.set_species_data(IP.num_species,IP.multispecies,IP.CFDkit_Path,
			     IP.debug_level,IP.Mach_Number_Reference,IP.Schmidt);
      IP.Wo.set_initial_values(IP.mass_fractions);
      IP.Uo.set_initial_values(IP.mass_fractions);

      //set proper Temp & Pressure instead of defaults
      IP.Wo.rho = IP.Pressure/(IP.Wo.Rtot()*IP.Temperature); 
      IP.Wo.p = IP.Pressure;	
      IP.Wo.v.zero();

      IP.Uo = U(IP.Wo);
    } 
    /*********************************************************
     ******************* CHEM2D END **************************
     *********************************************************/
   
    if(!CFDkit_Primary_MPI_Processor()) {   
      IP.Wo.v.x = IP.Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Flow_Angle/360.00);
      IP.Wo.v.y = IP.Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Flow_Angle/360.00);
    }

    /********************************************/
    Communicator.Bcast(IP.Flow_Type, 
                       INPUT_PARAMETER_LENGTH_CHEM2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.FlowType), 
                       1, 
                       MPI::INT, Source_Rank);

    /********************************************/
    Communicator.Bcast(IP.Flow_Geometry_Type, 
                       INPUT_PARAMETER_LENGTH_CHEM2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.Axisymmetric), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Gravity), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Global_Schmidt),
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.debug_level), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Preconditioning), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Moving_wall_velocity), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    /********************************************/
    Communicator.Bcast(IP.Grid_Type, 
                       INPUT_PARAMETER_LENGTH_CHEM2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.NACA_Aerofoil_Type, 
                       INPUT_PARAMETER_LENGTH_CHEM2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Grid), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Cells_Idir), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Cells_Jdir), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Ghost_Cells), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Blocks_Idir), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Blocks_Jdir), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Box_Width), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Box_Height), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Plate_Length), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Pipe_Length), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Pipe_Radius), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Blunt_Body_Radius), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Blunt_Body_Mach_Number), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Chamber_Length),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Chamber_Radius),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Chamber_To_Throat_Length),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Nozzle_Length),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Nozzle_Radius_Exit),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Nozzle_Radius_Throat),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Nozzle_Type),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Grain_Radius),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Cylinder_Radius), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Ellipse_Length_X_Axis), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Ellipse_Length_Y_Axis), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Chord_Length), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Orifice_Radius), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Wedge_Angle), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Wedge_Length), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Smooth_Bump),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.X_Shift.x), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.X_Shift.y), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.X_Scale), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.X_Rotate), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    // Boundary Conditions:
    Communicator.Bcast(IP.Boundary_Conditions_Specified,
		       INPUT_PARAMETER_LENGTH_CHEM2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(&(IP.BCs_Specified),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(IP.BC_North_Type,
		       INPUT_PARAMETER_LENGTH_CHEM2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(IP.BC_South_Type,
		       INPUT_PARAMETER_LENGTH_CHEM2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(IP.BC_East_Type,
		       INPUT_PARAMETER_LENGTH_CHEM2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(IP.BC_West_Type,
		       INPUT_PARAMETER_LENGTH_CHEM2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(&(IP.BC_North),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.BC_South),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.BC_East),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.BC_West),
		       1,
		       MPI::INT,Source_Rank);
    // ICEM:
    if (!(CFDkit_MPI::This_Processor_Number == Source_Rank)) {
       IP.ICEMCFD_FileNames = new char*[3];
       for (i = 0; i < 3; i++) {
          IP.ICEMCFD_FileNames[i] = new char[INPUT_PARAMETER_LENGTH_CHEM2D];
       } /* endfor */
    } /* endif */
    Communicator.Bcast(IP.ICEMCFD_FileNames[0], 
                       INPUT_PARAMETER_LENGTH_CHEM2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.ICEMCFD_FileNames[1], 
                       INPUT_PARAMETER_LENGTH_CHEM2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.ICEMCFD_FileNames[2], 
                       INPUT_PARAMETER_LENGTH_CHEM2D, 
                       MPI::CHAR, Source_Rank);

    // AMR & Refinement Parameters
    Communicator.Bcast(&(IP.AMR), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.AMR_Frequency),
                       1,
                       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Initial_Mesh_Refinements), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Uniform_Mesh_Refinements),
                       1,
                       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Boundary_Mesh_Refinements),
                       1,
                       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Maximum_Refinement_Level),
                       1,
                       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Minimum_Refinement_Level),
                       1,
                       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Threshold_for_Refinement), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Threshold_for_Coarsening), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Refinement_Criteria),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Refinement_Criteria_Gradient_Density),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Refinement_Criteria_Divergence_Velocity),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Refinement_Criteria_Curl_Velocity),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Refinement_Criteria_Gradient_Temperature),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Refinement_Criteria_Gradient_CH4),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Refinement_Criteria_Gradient_CO2),
		       1,
		       MPI::DOUBLE,Source_Rank);

    // Mesh stretching flag:
    Communicator.Bcast(&(IP.i_Mesh_Stretching),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Mesh_Stretching_Type_Idir),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Mesh_Stretching_Type_Jdir),		     
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Mesh_Stretching_Factor_Idir),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Mesh_Stretching_Factor_Jdir),
		       1,
		       MPI::DOUBLE,Source_Rank);
    // Smooth quad block flag:
    Communicator.Bcast(&(IP.i_Smooth_Quad_Block),
		       1,
		       MPI::INT,Source_Rank);

    // File Names
    Communicator.Bcast(IP.Output_File_Name, 
                       INPUT_PARAMETER_LENGTH_CHEM2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Grid_File_Name, 
                       INPUT_PARAMETER_LENGTH_CHEM2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Grid_Definition_File_Name, 
                       INPUT_PARAMETER_LENGTH_CHEM2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Restart_File_Name, 
                       INPUT_PARAMETER_LENGTH_CHEM2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Gnuplot_File_Name, 
                       INPUT_PARAMETER_LENGTH_CHEM2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Output_Format_Type, 
                       INPUT_PARAMETER_LENGTH_CHEM2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Output_Format), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Restart_Solution_Save_Frequency), 
                       1, 
                       MPI::INT, Source_Rank);
    // Output progress frequency:
    Communicator.Bcast(&(IP.Output_Progress_Frequency),
		       1,
		       MPI::INT,Source_Rank);
    // Multigrid Related Parameters
    IP.Multigrid_IP.Broadcast_Input_Parameters(Communicator,
					       Source_Rank);

    if (!(CFDkit_MPI::This_Processor_Number == Source_Rank)) {
       IP.Number_of_Processors = CFDkit_MPI::Number_of_Processors;
    } /* endif */
    Communicator.Bcast(&(IP.Number_of_Blocks_Per_Processor), 
                       1, 
                       MPI::INT, Source_Rank);
   // Freeze_Limiter
    Communicator.Bcast(&(IP.Freeze_Limiter), 
                       1, 
                       MPI::INT, Source_Rank);

    // Freeze_Limiter_Residual_Level
    Communicator.Bcast(&(IP.Freeze_Limiter_Residual_Level), 
                       1, 
                       MPI::DOUBLE, Source_Rank);

    // Reset the static variables.
    IP.Wo.set_flow_type(IP.FlowType);
    IP.Uo.set_flow_type(IP.FlowType);
    IP.Wo.set_turbulence_variables(IP.C_constant,
				   IP.von_Karman_Constant,
				   IP.yplus_sublayer,
				   IP.yplus_buffer_layer,
				   IP.yplus_outer_layer); 
    IP.Uo.set_turbulence_variables(IP.C_constant,
				   IP.von_Karman_Constant,
				   IP.yplus_sublayer,
				   IP.yplus_buffer_layer,
				   IP.yplus_outer_layer);

}
#endif

/********************************************************
 * Routine: Get_Next_Input_Control_Parameter            *
 *                                                      *
 * Get the next input control parameter from the input  *
 * file.                                                *
 *                                                      *
 ********************************************************/
void Get_Next_Input_Control_Parameter(Chem2D_Input_Parameters &IP) {

    int i;
    char buffer[256];

    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File.getline(buffer, sizeof(buffer));
    i = 0;
    if (buffer[0] != '#') {
 
       while (1) {
          if (buffer[i] == ' ' || buffer[i] == '=' ) break;
          i = i + 1;
          if (i > strlen(buffer) ) break;
       } /* endwhile */
       buffer[i] = '\0';
    } /* endif */
    strcpy(IP.Next_Control_Parameter, buffer);

    //    cout<<"\n "<<IP.Next_Control_Parameter<<endl; cout.flush();

}

/********************************************************
 * Routine: Parse_Next_Input_Control_Parameter          *
 *                                                      *
 * Parses and executes the next input control parameter *
 * from the input file.                                 *
 *                                                      *
 ********************************************************/
int Parse_Next_Input_Control_Parameter(Chem2D_Input_Parameters &IP) {

    int i_command;
    char buffer[256];

    i_command = 0;

    if (strcmp(IP.Next_Control_Parameter, "Time_Integration_Type") == 0) {
       i_command = 1;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Time_Integration_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Time_Integration_Type, "Explicit_Euler") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "Explicit_Predictor_Corrector") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR;
           IP.N_Stage = 2;
       } else if (strcmp(IP.Time_Integration_Type, "Explicit_Runge_Kutta") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_RUNGE_KUTTA;
           IP.N_Stage = 4;
       } else if (strcmp(IP.Time_Integration_Type, "Multistage_Optimal_Smoothing") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING;
           IP.N_Stage = 4;
       	   /* Multigrid */
       } else if (strcmp(IP.Time_Integration_Type, "Multigrid") == 0) {
	   IP.i_Time_Integration = TIME_STEPPING_MULTIGRID;
       } else {
           IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
           IP.N_Stage = 1;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Reconstruction_Type") == 0) {
       i_command = 2;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Reconstruction_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Reconstruction_Type, "Green_Gauss") == 0) {
          IP.i_Reconstruction = RECONSTRUCTION_GREEN_GAUSS;
       } else if (strcmp(IP.Reconstruction_Type, "Least_Squares") == 0) {
          IP.i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;
       } else {
          IP.i_Reconstruction = RECONSTRUCTION_GREEN_GAUSS;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Limiter_Type") == 0) {
       i_command = 3;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Limiter_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Limiter_Type, "One") == 0) {
          IP.i_Limiter = LIMITER_ONE;
       } else if (strcmp(IP.Limiter_Type, "Zero") == 0) {
          IP.i_Limiter = LIMITER_ZERO;
       } else if (strcmp(IP.Limiter_Type, "VanLeer") == 0) {
          IP.i_Limiter = LIMITER_VANLEER;
       } else if (strcmp(IP.Limiter_Type, "VanAlbada") == 0) {
          IP.i_Limiter = LIMITER_VANALBADA;
       } else if (strcmp(IP.Limiter_Type, "Barth_Jespersen") == 0) {
          IP.i_Limiter = LIMITER_BARTH_JESPERSEN;
       } else if (strcmp(IP.Limiter_Type, "Venkatakrishnan") == 0) {
          IP.i_Limiter = LIMITER_VENKATAKRISHNAN;
       } else {
          IP.i_Limiter = LIMITER_VANLEER ;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Flux_Function_Type") == 0) {
       i_command = 4;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Flux_Function_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Flux_Function_Type, "Godunov") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV;
       } else if (strcmp(IP.Flux_Function_Type, "Roe") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_ROE;
       } else if (strcmp(IP.Flux_Function_Type, "Rusanov") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_RUSANOV;
       } else if (strcmp(IP.Flux_Function_Type, "HLLE") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_HLLE;
       } else if (strcmp(IP.Flux_Function_Type, "Linde") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_LINDE;
       } else if (strcmp(IP.Flux_Function_Type, "HLLC") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_HLLC;
       } else {
          IP.i_Flux_Function = FLUX_FUNCTION_ROE;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "ICs_Type") == 0) {
       i_command = 5;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.ICs_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.ICs_Type, "Constant") == 0) {
          IP.i_ICs = IC_CONSTANT;
       } else if (strcmp(IP.ICs_Type, "Uniform") == 0) {
          IP.i_ICs = IC_UNIFORM;
       } else if (strcmp(IP.ICs_Type, "Sod") == 0) {
          IP.i_ICs = IC_SOD;
       } else if (strcmp(IP.ICs_Type, "Sod_Xdir") == 0) {
          IP.i_ICs = IC_SOD_XDIR;
       } else if (strcmp(IP.ICs_Type, "Sod_Ydir") == 0) {
          IP.i_ICs = IC_SOD_YDIR;
       } else if (strcmp(IP.ICs_Type, "Groth") == 0) {
          IP.i_ICs = IC_GROTH;
       } else if (strcmp(IP.ICs_Type, "Groth_Xdir") == 0) {
          IP.i_ICs = IC_GROTH_XDIR;
       } else if (strcmp(IP.ICs_Type, "Groth_Ydir") == 0) {
          IP.i_ICs = IC_GROTH_YDIR;
       } else if (strcmp(IP.ICs_Type, "Einfeldt") == 0) {
          IP.i_ICs = IC_EINFELDT;
       } else if (strcmp(IP.ICs_Type, "Einfeldt_Xdir") == 0) {
          IP.i_ICs = IC_EINFELDT_XDIR;
       } else if (strcmp(IP.ICs_Type, "Einfeldt_Ydir") == 0) {
          IP.i_ICs = IC_EINFELDT_YDIR;
       } else if (strcmp(IP.ICs_Type, "Shock_Wave_Xdir") == 0) {
          IP.i_ICs = IC_SHOCK_WAVE_XDIR;
       } else if (strcmp(IP.ICs_Type, "Shock_Wave_Ydir") == 0) {
          IP.i_ICs = IC_SHOCK_WAVE_YDIR;
       } else if (strcmp(IP.ICs_Type, "Contact_Surface_Xdir") == 0) {
          IP.i_ICs = IC_CONTACT_SURFACE_XDIR;
       } else if (strcmp(IP.ICs_Type, "Contact_Surface_Ydir") == 0) {
          IP.i_ICs = IC_CONTACT_SURFACE_YDIR;
       } else if (strcmp(IP.ICs_Type, "Rarefaction_Wave_Xdir") == 0) {
          IP.i_ICs = IC_RAREFACTION_WAVE_XDIR;
       } else if (strcmp(IP.ICs_Type, "Rarefaction_Wave_Ydir") == 0) {
          IP.i_ICs = IC_RAREFACTION_WAVE_YDIR;
       } else if (strcmp(IP.ICs_Type, "ShockBox") == 0) {
          IP.i_ICs = IC_SHOCK_BOX;
       } else if (strcmp(IP.ICs_Type, "High_Pressure_Reservoir") == 0) {
          IP.i_ICs = IC_HIGH_PRESSURE_RESERVOIR;
       } else if (strcmp(IP.ICs_Type, "Low_Pressure_Reservoir") == 0) {
          IP.i_ICs = IC_LOW_PRESSURE_RESERVOIR;
       } else if (strcmp(IP.ICs_Type, "Riemann") == 0) {
          IP.i_ICs = IC_RIEMANN;
       } else if (strcmp(IP.ICs_Type, "Riemann_Xdir") == 0) {
          IP.i_ICs = IC_RIEMANN_XDIR;
       } else if (strcmp(IP.ICs_Type, "Riemann_Ydir") == 0) {
	 IP.i_ICs = IC_RIEMANN_YDIR;  
       } else if (strcmp(IP.ICs_Type, "Wedge_Flow") == 0) {
	 IP.i_ICs = IC_WEDGE_FLOW;   
       } else if (strcmp(IP.ICs_Type,"Ringleb_Flow") == 0) {
	 IP.i_ICs = IC_RINGLEB_FLOW;
       } else if (strcmp(IP.ICs_Type,"Flat_Plate") == 0) {
	 IP.i_ICs = IC_VISCOUS_FLAT_PLATE;  
	 IP.BC_South = BC_WALL_VISCOUS_HEATFLUX;
       } else if (strcmp(IP.ICs_Type, "Couette") == 0 ){
	 IP.i_ICs = IC_VISCOUS_COUETTE; 
       } else if (strcmp(IP.ICs_Type, "Couette_with_Pressure_Gradient") == 0 ){
	 IP.i_ICs = IC_VISCOUS_COUETTE_PRESSURE_GRADIENT;
       /****************** CHEMD2D ******************************/
       } else if (strcmp(IP.ICs_Type, "Mix") == 0) {
	 IP.i_ICs = IC_GAS_MIX;
       } else if (strcmp(IP.ICs_Type, "Core_Flame") == 0 ){
	 IP.i_ICs = IC_CHEM_CORE_FLAME;
       } else if (strcmp(IP.ICs_Type, "Inverse_Flame") == 0 ){
	 IP.i_ICs = IC_CHEM_INVERSE_FLAME ; 
       } else if (strcmp(IP.ICs_Type, "Pressure_Gradient_x") == 0 ){
	 IP.i_ICs = IC_PRESSURE_GRADIENT_X;
       } else if (strcmp(IP.ICs_Type, "Pressure_Gradient_y") == 0 ){
	 IP.i_ICs = IC_PRESSURE_GRADIENT_Y;
       } else if (strcmp(IP.ICs_Type, "1DPremixedFlame") == 0 ){
	 IP.i_ICs = IC_CHEM_1DFLAME;  
       } else if (strcmp(IP.ICs_Type, "Pipe_Flow") == 0) {
          IP.i_ICs = IC_VISCOUS_PIPE_FLOW;
       /**************** END - CHEMD2D **************************/
       } else if (strcmp(IP.ICs_Type, "Restart") == 0) {
          IP.i_ICs = IC_RESTART;
       } else {
          IP.i_ICs = IC_UNIFORM;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Grid_Type") == 0) {
       i_command = 6;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Grid_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Grid_Type, "Cartesian") == 0) {
          IP.i_Grid = GRID_CARTESIAN_UNIFORM;
          IP.Box_Width = ONE;
          IP.Box_Height = ONE;
       } else if (strcmp(IP.Grid_Type, "Square") == 0) {
          IP.i_Grid = GRID_SQUARE;
          IP.Box_Width = ONE;
          IP.Box_Height = ONE;
       } else if (strcmp(IP.Grid_Type, "Rectangular_Box") == 0) {
          IP.i_Grid = GRID_RECTANGULAR_BOX;
          IP.Box_Width = ONE;
          IP.Box_Height = ONE;
      } else if (strcmp(IP.Grid_Type, "Couette") == 0) {
          IP.i_Grid = GRID_COUETTE;
          IP.Box_Width = ONE;
          IP.Box_Height = 0.001;
      } else if (strcmp(IP.Grid_Type, "1DFlame") == 0) {
          IP.i_Grid = GRID_1DFLAME;
          IP.Box_Width = ONE;
          IP.Box_Height = 0.1;
      } else if (strcmp(IP.Grid_Type, "Laminar_Flame") == 0) {
          IP.i_Grid = GRID_LAMINAR_FLAME;
          IP.Pipe_Length = 0.1;
          IP.Pipe_Radius = 0.2;
      } else if (strcmp(IP.Grid_Type,"Flat_Plate") == 0 ||
		 strcmp(IP.Grid_Type,"Adiabatic_Flat_Plate") == 0) {
	  IP.i_Grid = GRID_FLAT_PLATE;
	  IP.Plate_Length = ONE;
	  IP.BC_South = BC_WALL_VISCOUS_HEATFLUX;
      } else if (strcmp(IP.Grid_Type,"Isothermal_Flat_Plate") == 0) {
	  IP.i_Grid = GRID_FLAT_PLATE;
	  IP.Plate_Length = ONE;
	  IP.BC_South = BC_WALL_VISCOUS_ISOTHERMAL;
       } else if (strcmp(IP.Grid_Type, "Pipe") == 0) {
          IP.i_Grid = GRID_PIPE;
          IP.Pipe_Length = ONE;
          IP.Pipe_Radius = 0.1;
       } else if (strcmp(IP.Grid_Type, "Blunt_Body") == 0) {
          IP.i_Grid = GRID_BLUNT_BODY;
          IP.Blunt_Body_Radius = ONE;
          IP.Blunt_Body_Mach_Number = TWO;
       } else if (strcmp(IP.Grid_Type, "Rocket_Motor") == 0) {
	  IP.i_Grid = GRID_ROCKET_MOTOR;
	  IP.Chamber_Length = 0.835;
	  IP.Chamber_Radius = 0.020;
	  IP.Chamber_To_Throat_Length = 0.05;
	  IP.Nozzle_Length = 0.150;
	  IP.Nozzle_Radius_Exit = 0.030;
	  IP.Nozzle_Radius_Throat = 0.010;
	  IP.Nozzle_Type = NOZZLE_GOTTLIEB_FUNCTION;
	  IP.Grain_Radius = 0.0;
	  IP.BC_North = BC_WALL_VISCOUS_ISOTHERMAL;
       } else if (strcmp(IP.Grid_Type, "Circular_Cylinder") == 0) {
          IP.i_Grid = GRID_CIRCULAR_CYLINDER;
          IP.Cylinder_Radius = ONE;
	  IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
	  IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
	  IP.Mesh_Stretching_Factor_Idir = 1.025;
	  IP.Mesh_Stretching_Factor_Jdir = 1.001;
       } else if (strcmp(IP.Grid_Type, "Ellipse") == 0) {
          IP.i_Grid = GRID_ELLIPSE;
          IP.Ellipse_Length_X_Axis = TWO;
	  IP.Ellipse_Length_Y_Axis = HALF;
       } else if (strcmp(IP.Grid_Type, "NACA_Aerofoil") == 0) {
          IP.i_Grid = GRID_NACA_AEROFOIL;
          IP.Chord_Length = ONE;
          strcpy(IP.NACA_Aerofoil_Type, "0012");
       } else if (strcmp(IP.Grid_Type, "Free_Jet") == 0) {
          IP.i_Grid = GRID_FREE_JET;
          IP.Orifice_Radius = ONE;
       } else if (strcmp(IP.Grid_Type, "Wedge") == 0) {
	  IP.i_Grid = GRID_WEDGE;
	  IP.Wedge_Angle = 25.0;
	  IP.Wedge_Length = HALF;
	  IP.BC_South = BC_WALL_VISCOUS_ISOTHERMAL;
       } else if (strcmp(IP.Grid_Type,"Ringleb_Flow") == 0) {
	 IP.i_Grid = GRID_RINGLEB_FLOW;
       } else if (strcmp(IP.Grid_Type,"Bump_Channel_Flow") == 0) {
	 IP.i_Grid = GRID_BUMP_CHANNEL_FLOW;
       } else if (strcmp(IP.Grid_Type,"Non_Smooth_Bump_Channel_Flow") == 0) {
	 IP.i_Grid = GRID_BUMP_CHANNEL_FLOW;
	 IP.Smooth_Bump = OFF;
       } else if (strcmp(IP.Grid_Type,"Smooth_Bump_Channel_Flow") == 0) {
	 IP.i_Grid = GRID_BUMP_CHANNEL_FLOW;
	 IP.Smooth_Bump = ON;
       } else if (strcmp(IP.Grid_Type, "ICEMCFD") == 0) {
          IP.i_Grid = GRID_ICEMCFD;
       } else if (strcmp(IP.Grid_Type, "Read_From_Definition_File") == 0) {
          IP.i_Grid = GRID_READ_FROM_DEFINITION_FILE;
       } else if (strcmp(IP.Grid_Type, "Read_From_Data_File") == 0) {
          IP.i_Grid = GRID_READ_FROM_GRID_DATA_FILE;
       } else {
          IP.i_Grid = GRID_SQUARE;
          IP.Box_Width = ONE;
          IP.Box_Height = ONE;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Output_File_Name") == 0) {
       i_command = 7;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Output_File_Name, 
              IP.Next_Control_Parameter);
       strcat(IP.Output_File_Name, ".dat");
       strcpy(IP.Grid_File_Name, 
              IP.Next_Control_Parameter);
       strcat(IP.Grid_File_Name, ".grid");
       strcpy(IP.Grid_Definition_File_Name, 
              IP.Next_Control_Parameter);
       strcat(IP.Grid_Definition_File_Name, ".griddef");
       strcpy(IP.Restart_File_Name, 
              IP.Next_Control_Parameter);
       strcat(IP.Restart_File_Name, ".soln");
       strcpy(IP.Gnuplot_File_Name, 
              IP.Next_Control_Parameter);
       strcat(IP.Gnuplot_File_Name, ".gplt");

    } else if (strcmp(IP.Next_Control_Parameter, "Grid_File_Name") == 0) {
       i_command = 8;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Grid_File_Name, 
              IP.Next_Control_Parameter);
       strcat(IP.Grid_File_Name, ".grid");
       strcpy(IP.Grid_Definition_File_Name, 
              IP.Next_Control_Parameter);
       strcat(IP.Grid_Definition_File_Name, ".griddef");

    } else if (strcmp(IP.Next_Control_Parameter, "Restart_File_Name") == 0) {
       i_command = 9;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Restart_File_Name, 
              IP.Next_Control_Parameter);
       strcat(IP.Restart_File_Name, ".soln");

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Cells_Idir") == 0) {
       i_command = 10;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Cells_Idir;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Cells_Idir < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Cells_Jdir") == 0) {
       i_command = 11;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Cells_Jdir;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Cells_Jdir < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Ghost_Cells") == 0) {
       i_command = 12;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Ghost_Cells;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Ghost_Cells < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Blocks_Idir") == 0) {
       i_command = 12;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Blocks_Idir;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Blocks_Idir < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Blocks_Jdir") == 0) {
       i_command = 13;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Blocks_Jdir;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Blocks_Jdir < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Time_Accurate") == 0) {
       i_command = 14;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Time_Accurate;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Time_Accurate != 0 &&
           IP.Time_Accurate != 1) IP.Time_Accurate = 0;
       if (IP.Time_Accurate) {
          IP.Local_Time_Stepping = GLOBAL_TIME_STEPPING;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Local_Time_Stepping") == 0) {
       i_command = 15;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Local_Time_Stepping;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Local_Time_Stepping != GLOBAL_TIME_STEPPING &&
           IP.Local_Time_Stepping != SCALAR_LOCAL_TIME_STEPPING &&
	   IP.Local_Time_Stepping != MATRIX_LOCAL_TIME_STEPPING &&	   
           IP.Local_Time_Stepping != LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER &&
	   IP.Local_Time_Stepping != SEMI_IMPLICIT_LOCAL_TIME_STEPPING &&
	   IP.Local_Time_Stepping != SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER) {
         IP.Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
	}
       if(IP.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER ||
       	  IP.Local_Time_Stepping == SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER){
	 IP.Preconditioning = 1;
       }       

    } else if (strcmp(IP.Next_Control_Parameter, "Maximum_Number_of_Time_Steps") == 0) {
       i_command = 16;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Maximum_Number_of_Time_Steps;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Maximum_Number_of_Time_Steps < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "N_Stage") == 0) {
       i_command = 17;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.N_Stage;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.N_Stage < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "CFL_Number") == 0) {
       i_command = 18;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.CFL_Number;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.CFL_Number <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Box_Width") == 0) {
       i_command = 19;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Box_Width;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Box_Width <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Box_Height") == 0) {
       i_command = 20;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Box_Height;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Box_Height <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Plate_Length") == 0) {
       i_command = 21;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Plate_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Plate_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Pipe_Length") == 0) {
       i_command = 22;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Pipe_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Pipe_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Pipe_Radius") == 0) {
       i_command = 23;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Pipe_Radius;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Pipe_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Blunt_Body_Radius") == 0) {
       i_command = 24;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Blunt_Body_Radius;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Blunt_Body_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Blunt_Body_Mach_Number") == 0) {
       i_command = 25;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Blunt_Body_Mach_Number;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Blunt_Body_Mach_Number <= ONE) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Cylinder_Radius") == 0) {
       i_command = 26;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Cylinder_Radius;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Cylinder_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Ellipse_Length_X_Axis") == 0) {
       i_command = 27;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Ellipse_Length_X_Axis;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Ellipse_Length_X_Axis <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Ellipse_Length_Y_Axis") == 0) {
       i_command = 28;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Ellipse_Length_Y_Axis;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Ellipse_Length_Y_Axis <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Chord_Length") == 0) {
       i_command = 29;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Chord_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Chord_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "NACA_Aerofoil_Type") == 0) {
       i_command = 30;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.NACA_Aerofoil_Type, 
              IP.Next_Control_Parameter);
       if (strlen(IP.NACA_Aerofoil_Type) != 4 &&
           strlen(IP.NACA_Aerofoil_Type) != 5) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Orifice_Radius") == 0) {
       i_command = 31;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Orifice_Radius;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Orifice_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Wedge_Angle") == 0) {
       i_command = 31;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Wedge_Angle;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Wedge_Angle <= ZERO) i_command = INVALID_INPUT_VALUE;
    
    } else if (strcmp(IP.Next_Control_Parameter, "Wedge_Length") == 0) {
       i_command = 31;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Wedge_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Wedge_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Smooth_Bump") == 0) {
       i_command = 32;
       Get_Next_Input_Control_Parameter(IP);
       if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	 IP.Smooth_Bump = ON;
       } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	 IP.Smooth_Bump = OFF;
       } else {
	 i_command = INVALID_INPUT_VALUE;
       }

    } else if (strcmp(IP.Next_Control_Parameter,"Chamber_Length") == 0) {
       i_command = 33;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Chamber_Length;
       IP.Input_File.getline(buffer,sizeof(buffer));
       if (IP.Chamber_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Chamber_Radius") == 0) {
       i_command = 34;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Chamber_Radius;
       IP.Input_File.getline(buffer,sizeof(buffer));
       if (IP.Chamber_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Chamber_To_Throat_Length") == 0) {
       i_command = 35;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Chamber_To_Throat_Length;
       IP.Input_File.getline(buffer,sizeof(buffer));
       if (IP.Chamber_To_Throat_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Nozzle_Length") == 0) {
       i_command = 36;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Nozzle_Length;
       IP.Input_File.getline(buffer,sizeof(buffer));
       if (IP.Nozzle_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Nozzle_Radius_Exit") == 0) {
       i_command = 37;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Nozzle_Radius_Exit;
       IP.Input_File.getline(buffer,sizeof(buffer));
       if (IP.Nozzle_Radius_Exit <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Nozzle_Radius_Throat") == 0) {
       i_command = 38;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Nozzle_Radius_Throat;
       IP.Input_File.getline(buffer,sizeof(buffer));
       if (IP.Nozzle_Radius_Throat <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Nozzle_Type") == 0) {
       i_command = 39;
       Get_Next_Input_Control_Parameter(IP);
       if (strcmp(IP.Next_Control_Parameter,"Conical") == 0) {
	 IP.Nozzle_Type = NOZZLE_CONICAL;
       } else if (strcmp(IP.Next_Control_Parameter,"Gottlieb") == 0) {
	 IP.Nozzle_Type = NOZZLE_GOTTLIEB_FUNCTION;
       } else if (strcmp(IP.Next_Control_Parameter,"Quartic") == 0) {
	 IP.Nozzle_Type = NOZZLE_QUARTIC_FUNCTION;
       } else if (strcmp(IP.Next_Control_Parameter,"Hybrid_Conical") == 0) {
	 IP.Nozzle_Type = NOZZLE_HYBRID_CONICAL_FUNCTION;
       } else {
	 i_command = INVALID_INPUT_VALUE;
       }

    } else if (strcmp(IP.Next_Control_Parameter,"Grain_Radius") == 0) {
       i_command = 40;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Grain_Radius;
       IP.Input_File.getline(buffer,sizeof(buffer));

       /***********************************************************************
        ************************ CHEM2D SPECIFIC ******************************
        For Chem2D the following input parameters are added:
             "Reaction_Mechanism" is 

             "User_Reaciton_Mechanism"
 
             "Species"

       *************************************************************************
       *************************************************************************/
       
       /*************************************/
       /**** REACTIONS SET FOR HARDCODED ****/
       /*************************************/
    } else if (strcmp(IP.Next_Control_Parameter, "Reaction_Mechanism") == 0) {
       i_command = 200;
       Get_Next_Input_Control_Parameter(IP);
       IP.Deallocate();  //DEALLOCATE BEFORE CHANGING num_species
       int flag =0;

       //convert IP to string & define setup which Reaction Mechanism
       IP.react_name = IP.Next_Control_Parameter;
       IP.Wo.React.set_reactions(IP.react_name);
   
       IP.num_species = IP.Wo.React.num_species;      
       IP.Allocate();

       //Get species and load appropriate data
       for(int i=0; i<IP.num_species; i++){
	 IP.multispecies[i] = IP.Wo.React.species[i];
	 IP.Schmidt[i] = IP.Global_Schmidt;
       }   

       //Get next line and read in Schmidt numbers else will use defaults
       Get_Next_Input_Control_Parameter(IP);
       if (strcmp(IP.Next_Control_Parameter, "Schmidt_Numbers") == 0){
	 for(int i=0; i<IP.num_species; i++){
	   IP.Input_File >> IP.Schmidt[i];	 	 
	 }	 
	 //fudge the line number and istream counters
	 IP.Input_File.getline(buffer, sizeof(buffer)); 
	 IP.Line_Number = IP.Line_Number + 1 ;
	 flag = 1;
       } else { //Set to one value ~1
	 for(int i=0; i<IP.num_species; i++){
	   IP.Schmidt[i] = IP.Global_Schmidt;	 	 
	 }
	 //To fix up line numbers
	 IP.Line_Number = IP.Line_Number - 1 ;
       }

       //Set appropriate species data
       IP.Wo.set_species_data(IP.num_species,IP.multispecies,IP.CFDkit_Path,
			      IP.debug_level,IP.Mach_Number_Reference,IP.Schmidt); 
       IP.Uo.set_species_data(IP.num_species,IP.multispecies,IP.CFDkit_Path,
			      IP.debug_level,IP.Mach_Number_Reference,IP.Schmidt);
  
       //Get next line and read in mass fractions or set defaults
       if(flag){
	 Get_Next_Input_Control_Parameter(IP);
       }
       if (strcmp(IP.Next_Control_Parameter, "Mass_Fractions") == 0){
	 //Get Initial Mass Fractions from user	
	 double temp=0.0;
	 for(int i=0; i<IP.num_species; i++){
	   IP.Input_File >> IP.mass_fractions[i];
	   temp += IP.mass_fractions[i];
	 }
	 //check to make sure it adds to 1
	 if(temp < ONE-MICRO || temp > ONE+MICRO){ 
	   cout<<"\n Mass Fractions summed to "<<temp<<". Should sum to 1\n";
	   i_command = INVALID_INPUT_VALUE;
	 }
	 //Set inital Values; 
	 IP.Wo.set_initial_values(IP.mass_fractions);  
	 IP.Uo.set_initial_values(IP.mass_fractions);  
	 IP.Uo = U(IP.Wo);
	 
	 //fudge the line number and istream counters
	 IP.Input_File.getline(buffer, sizeof(buffer));  
	 IP.Line_Number = IP.Line_Number + 1; 
        
       //Spit out appropriate mass fractions and exit
       } else if (strcmp(IP.Next_Control_Parameter, "Equivalence_Ratio") == 0){        
	 double phi;
	 IP.Input_File >> phi;
	 Equivalence_Ratio(phi);

	 //If no mass fraction data is set to defaults (all equal to 1/num_species) 
       } else {
	 IP.Uo = U(IP.Wo);
	 IP.Line_Number = IP.Line_Number - 1 ;
       }
         
       
       /***************************************/
       /**** REACTIONS FOR USER DEFINED *******/ 
       /***************************************/
    } else if (strcmp(IP.Next_Control_Parameter, "User_Reaction_Mechanism") == 0) { 
      // this will be added but its not quite yet
      i_command=201;
      cout<<endl<<IP.Next_Control_Parameter<<"\n not currently available in freeware version :)\n";
      i_command = INVALID_INPUT_VALUE;   
      
      /******************************************/
      /**** NON REACTING, BUT MULTIPLE GASES ****/
      /******************************************/
    } else if (strcmp(IP.Next_Control_Parameter, "Species") == 0) { 
      i_command = 203;
  
      IP.Deallocate();  //DEALLOCATE BEFORE CHANGING num_species
 
      // Non Reaction case so set NO_REACTIONS flag in reactions class
      IP.react_name ="NO_REACTIONS";
      IP.Wo.React.set_reactions(IP.react_name);
     
      //read in the number of species (should be first in line) 
      IP.Input_File>>IP.num_species;     
      IP.Allocate();

      //read in species names
      for(int i=0; i<IP.num_species; i++){
	IP.Input_File >> IP.multispecies[i];
      }
      
      //copy names into Reaction class for storage
      IP.Wo.React.set_species(IP.multispecies,IP.num_species);
     
      //Setup State class data and find species thermo and transport properties
      IP.Wo.set_species_data(IP.num_species,IP.multispecies,IP.CFDkit_Path,
			     IP.debug_level,IP.Mach_Number_Reference,IP.Schmidt); 
      IP.Uo.set_species_data(IP.num_species,IP.multispecies,IP.CFDkit_Path,
			     IP.debug_level,IP.Mach_Number_Reference,IP.Schmidt);
    
      //More Fudging of lines 
      IP.Line_Number = IP.Line_Number + 1 ;
      IP.Input_File.getline(buffer, sizeof(buffer));  
      
      //Get next line and read in mass fractions or set defaults
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter, "Mass_Fractions") == 0){
	 //Get Initial Mass Fractions from user 
	 double temp=0.0;
	 for(int i=0; i<IP.num_species; i++){
	   IP.Input_File >> IP.mass_fractions[i];
	   temp += IP.mass_fractions[i];
	 }
	 //check to make sure it adds to 1
	 if(temp < ONE-MICRO || temp > ONE+MICRO){ 
	   cout<<"\n Mass Fractions summed to "<<temp<<". Should be sum to 1\n";
	   i_command = INVALID_INPUT_VALUE;
	 }
	 //Set inital Values; 
	 IP.Wo.set_initial_values(IP.mass_fractions);  
	 IP.Uo.set_initial_values(IP.mass_fractions);
	 IP.Uo = U(IP.Wo);
	 
	 //fudge the line number and istream counters
	 IP.Input_File.getline(buffer, sizeof(buffer));  
	 IP.Line_Number = IP.Line_Number + 1; 
       } 
       //If no mass fraction data is set to defaults (all equal to 1/num_species)
       else{        
	 IP.Uo = U(IP.Wo);
	 IP.Line_Number = IP.Line_Number - 1 ;
       }
    
      /************* TEMPERATURE *************/
    } else if (strcmp(IP.Next_Control_Parameter, "Temperature") == 0) {
      i_command = 42;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Temperature;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Temperature <= ZERO) {
	i_command = INVALID_INPUT_VALUE;
      } else {
	IP.Wo.rho = IP.Pressure/(IP.Wo.Rtot()*IP.Temperature); 	
	//IP.Wo.v.zero();
      } /* endif */
   
      /************* PRESSURE ****************/
    } else if (strcmp(IP.Next_Control_Parameter, "Pressure") == 0) {
       i_command = 41;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Pressure;
       IP.Input_File.getline(buffer, sizeof(buffer));
       IP.Pressure = IP.Pressure*THOUSAND;
       if (IP.Pressure <= ZERO) {
	 i_command = INVALID_INPUT_VALUE;
       } else {
	 IP.Wo.rho = IP.Pressure/(IP.Wo.Rtot()*IP.Temperature); 
	 IP.Wo.p = IP.Pressure;	
	 //IP.Wo.v.zero();
       } /* endif */
      
       /***********************************************************************
	**************** END CHEM2D MODIFICATIONS *****************************
	***********************************************************************/

      /********** MACH NUMBER ***********/
    } else if (strcmp(IP.Next_Control_Parameter, "Mach_Number") == 0) {
      i_command = 39;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Mach_Number;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Mach_Number < ZERO) {
	i_command = INVALID_INPUT_VALUE;
      } else {
	IP.Mach_Number_Reference = IP.Mach_Number;
	IP.Wo.Mref = IP.Mach_Number;
	IP.Uo.Mref = IP.Mach_Number;
	IP.Wo.v.x = IP.Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Flow_Angle/360.00);
	IP.Wo.v.y = IP.Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Flow_Angle/360.00);
      } /* endif */
      
      /********** MACH NUMBER REFERENCE (for Low Mach Number Preconditioning ********/

    } else if (strcmp(IP.Next_Control_Parameter, "Mach_Number_Reference") == 0) {
      i_command = 40;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Mach_Number_Reference;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Mach_Number_Reference < ZERO) {
	i_command = INVALID_INPUT_VALUE;
      } else {
	IP.Wo.Mref = IP.Mach_Number_Reference;
	IP.Uo.Mref = IP.Mach_Number_Reference;
      } /* endif */
  
    } else if (strcmp(IP.Next_Control_Parameter, "Flow_Angle") == 0) {
       i_command = 41;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Flow_Angle;
       IP.Input_File.getline(buffer, sizeof(buffer));
       IP.Wo.v.x = IP.Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Flow_Angle/360.00);
       IP.Wo.v.y = IP.Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Flow_Angle/360.00);
       
    } else if (strcmp(IP.Next_Control_Parameter, "Time_Max") == 0) {
       i_command = 43;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Time_Max;
       IP.Input_File.getline(buffer, sizeof(buffer));
       IP.Time_Max = IP.Time_Max/THOUSAND;
       if (IP.Time_Max < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Blocks_Per_Processor") == 0) {
       i_command = 44;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Blocks_Per_Processor;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Blocks_Per_Processor < 1) i_command = INVALID_INPUT_VALUE;
 
    } else if (strcmp(IP.Next_Control_Parameter, "Output_Format_Type") == 0) {
       i_command = 45;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Output_Format_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Output_Format_Type, "Gnuplot") == 0  ||
           strcmp(IP.Output_Format_Type, "gnuplot") == 0  ||
           strcmp(IP.Output_Format_Type, "GNUPLOT") == 0) {
          IP.i_Output_Format = IO_GNUPLOT;
       } else if (strcmp(IP.Output_Format_Type, "Tecplot") == 0  ||
                  strcmp(IP.Output_Format_Type, "tecplot") == 0  ||
                  strcmp(IP.Output_Format_Type, "TECPLOT") == 0) {
          IP.i_Output_Format = IO_TECPLOT;
       } else if (strcmp(IP.Output_Format_Type, "Matlab") == 0  ||
                  strcmp(IP.Output_Format_Type, "matlab") == 0  ||
                  strcmp(IP.Output_Format_Type, "MATLAB") == 0) {
          IP.i_Output_Format = IO_MATLAB;
       } else if (strcmp(IP.Output_Format_Type, "Octave") == 0  ||
                  strcmp(IP.Output_Format_Type, "octave") == 0  ||
                  strcmp(IP.Output_Format_Type, "OCTAVE") == 0) {
          IP.i_Output_Format = IO_OCTAVE;
       } else {
          IP.i_Output_Format = IO_TECPLOT;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Flow_Type") == 0) {
       i_command = 46;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Flow_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Flow_Geometry_Type, "Inviscid") == 0) {
          IP.FlowType = FLOWTYPE_INVISCID;
       } else if (strcmp(IP.Flow_Type, "Laminar") == 0) {
	 IP.FlowType = FLOWTYPE_LAMINAR;
       } else if (strcmp(IP.Flow_Type, "Turbulent-k-epsilon") == 0) {
	 IP.FlowType = FLOWTYPE_TURBULENT_RANS_K_EPSILON;
       } else if (strcmp(IP.Flow_Type, "Turbulent-k-omega") == 0) {
	 IP.FlowType = FLOWTYPE_TURBULENT_RANS_K_OMEGA;
       } else if (strcmp(IP.Flow_Type, "Turbulent-LES") == 0) {
	 IP.FlowType = FLOWTYPE_TURBULENT_LES;
       } else if (strcmp(IP.Flow_Type, "Turbulent-DES-k-omega") == 0) {
	 IP.FlowType = FLOWTYPE_TURBULENT_DES_K_OMEGA;
       } else if (strcmp(IP.Flow_Type, "Turbulent-DNS") == 0) {
	 IP.FlowType = FLOWTYPE_TURBULENT_DNS;
       } else {
         IP.FlowType = FLOWTYPE_INVISCID;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Flow_Geometry_Type") == 0) {
       i_command = 50;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Flow_Geometry_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Flow_Geometry_Type, "Planar") == 0) {
          IP.Axisymmetric = 0;
       } else if (strcmp(IP.Flow_Geometry_Type, "Axisymmetric") == 0) {
	 IP.Axisymmetric = 1;
       } else if (strcmp(IP.Flow_Geometry_Type, "Axisymmetric-x") == 0) {
	 IP.Axisymmetric = 2;
       } else if (strcmp(IP.Flow_Geometry_Type, "Axisymmetric-y") == 0) {
	 IP.Axisymmetric = 1;
       } else {
          IP.Axisymmetric = 0;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Restart_Solution_Save_Frequency") == 0) {
       i_command = 51;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Restart_Solution_Save_Frequency;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Restart_Solution_Save_Frequency < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Output_Progress_Frequency") == 0) {
       i_command = 51;
       IP.Line_Number++;
       IP.Input_File >> IP.Output_Progress_Frequency;
       IP.Input_File.getline(buffer,sizeof(buffer));
       if (IP.Output_Progress_Frequency < 1) i_command = INVALID_INPUT_VALUE;

       /******* CHEM2D *********/
    } else if (strcmp(IP.Next_Control_Parameter, "Viscous_Flux_Evaluation_Type") == 0) {
      i_command = 510; 
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Viscous_Flux_Evaluation_Type,IP.Next_Control_Parameter);
      if (strcmp(IP.Viscous_Flux_Evaluation_Type, "Arithmetic_Mean") == 0) {
	IP.i_Viscous_Flux_Evaluation = VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE; 
      } else if (strcmp(IP.Viscous_Flux_Evaluation_Type, "Cartesian") == 0) {
	IP.i_Viscous_Flux_Evaluation = VISCOUS_RECONSTRUCTION_CARTESIAN; 
      } else if (strcmp(IP.Viscous_Flux_Evaluation_Type,"Diamond_Path") == 0) {
	IP.i_Viscous_Flux_Evaluation = VISCOUS_RECONSTRUCTION_DIAMOND_PATH;
      } else {
	i_command = INVALID_INPUT_VALUE;
	IP.i_Viscous_Flux_Evaluation = VISCOUS_RECONSTRUCTION_DIAMOND_PATH; 	
      } 

    } else if (strcmp(IP.Next_Control_Parameter, "Gravity") == 0) {
      i_command = 511; 
      IP.Gravity = 1;
     
    } else if (strcmp(IP.Next_Control_Parameter, "Schmidt") == 0) {
      i_command = 512; 
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Global_Schmidt;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Global_Schmidt < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Debug_Level") == 0) {
      i_command = 513; 
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.debug_level;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.debug_level != 0 && IP.debug_level != 1 ) i_command = INVALID_INPUT_VALUE; 

    } else if (strcmp(IP.Next_Control_Parameter, "Moving_Wall_Velocity") == 0) {
      i_command = 514;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Moving_wall_velocity;
      IP.Input_File.getline(buffer, sizeof(buffer));
    } else if (strcmp(IP.Next_Control_Parameter, "Time_Accurate_Plot_Frequency") == 0) {
       i_command = 515; 
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Time_Accurate_Plot_Frequency;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Time_Accurate_Plot_Frequency < 0) i_command = INVALID_INPUT_VALUE;
     
      /***** END CHEM2D *******/

    } else if (strcmp(IP.Next_Control_Parameter, "ICEMCFD_Topology_File") == 0) {
      i_command = 52;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.ICEMCFD_FileNames[0], IP.Next_Control_Parameter);
      
    } else if (strcmp(IP.Next_Control_Parameter, "ICEMCFD_Family_Boco_File") == 0) {
       i_command = 53;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.ICEMCFD_FileNames[1], IP.Next_Control_Parameter);

    } else if (strcmp(IP.Next_Control_Parameter, "ICEMCFD_Family_Topo_File") == 0) {
       i_command = 54;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.ICEMCFD_FileNames[2], IP.Next_Control_Parameter);

    } else if (strcmp(IP.Next_Control_Parameter, "X_Shift") == 0) {
       i_command = 55;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.X_Shift;
       IP.Input_File.setf(ios::skipws);
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "X_Scale") == 0) {
       i_command = 56;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.X_Scale;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "X_Rotate") == 0) {
       i_command = 57;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.X_Rotate;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Freeze_Limiter") == 0) {
      // Freeze_Limiter:
      i_command = 68;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Freeze_Limiter;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Freeze_Limiter < ZERO) i_command = INVALID_INPUT_VALUE;
      
    } else if (strcmp(IP.Next_Control_Parameter, "Freeze_Limiter_Residual_Level") == 0) {
      // Freeze_Limiter_Residual_Level:
      i_command = 69;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Freeze_Limiter_Residual_Level;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Freeze_Limiter_Residual_Level < 0) i_command = INVALID_INPUT_VALUE;
      
   } else if (strcmp(IP.Next_Control_Parameter, "Residual_Smoothing_Epsilon") == 0) {
      i_command = 70;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Residual_Smoothing_Epsilon;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Residual_Smoothing_Epsilon <= ZERO) {
         IP.Residual_Smoothing = 0;
         IP.Residual_Smoothing_Epsilon = ZERO;
      } else {
         IP.Residual_Smoothing = 1;
      } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Residual_Smoothing_Gauss_Seidel_Iterations") == 0) {
      i_command = 71;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Residual_Smoothing_Gauss_Seidel_Iterations;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Residual_Smoothing_Gauss_Seidel_Iterations < 0) {
         IP.Residual_Smoothing_Gauss_Seidel_Iterations = 0;
      } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter,"Mesh_Stretching") == 0) {
      i_command = 101;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	IP.i_Mesh_Stretching = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	IP.i_Mesh_Stretching = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Mesh_Stretching_Type_Idir") == 0) {
      i_command = 102;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"Linear") == 0) {
	IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_LINEAR;
      } else if (strcmp(IP.Next_Control_Parameter,"Min_Clustering") == 0) {
	IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MIN_CLUSTERING;
      } else if (strcmp(IP.Next_Control_Parameter,"Max_Clustering") == 0) {
	IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MAX_CLUSTERING;
      } else if (strcmp(IP.Next_Control_Parameter,"MinMax_Clustering") == 0) {
	IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
      } else if (strcmp(IP.Next_Control_Parameter,"Midpt_Clustering") == 0) {
	IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MIDPT_CLUSTERING;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Mesh_Stretching_Type_Jdir") == 0) {
      i_command = 103;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"Linear") == 0) {
	IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_LINEAR;
      } else if (strcmp(IP.Next_Control_Parameter,"Min_Clustering") == 0) {
	IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
      } else if (strcmp(IP.Next_Control_Parameter,"Max_Clustering") == 0) {
	IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MAX_CLUSTERING;
      } else if (strcmp(IP.Next_Control_Parameter,"MinMax_Clustering") == 0) {
	IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MINMAX_CLUSTERING;
      } else if (strcmp(IP.Next_Control_Parameter,"Midpt_Clustering") == 0) {
	IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIDPT_CLUSTERING;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Mesh_Stretching_Factor_Idir") == 0) {
      i_command = 104;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Mesh_Stretching_Factor_Idir;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Mesh_Stretching_Factor_Idir < ONE) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Mesh_Stretching_Factor_Jdir") == 0) {
      i_command = 105;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Mesh_Stretching_Factor_Jdir;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Mesh_Stretching_Factor_Jdir < ONE) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Smooth_Quad_Block") == 0) {
      i_command = 106;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	IP.i_Smooth_Quad_Block = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	IP.i_Smooth_Quad_Block = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter, "AMR") == 0) {
      i_command = 72;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	IP.AMR = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	IP.AMR = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter, "AMR_Frequency") == 0) {
      i_command = 73;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.AMR_Frequency;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.AMR_Frequency < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Initial_Mesh_Refinements") == 0) {
      i_command = 74;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Number_of_Initial_Mesh_Refinements;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Number_of_Initial_Mesh_Refinements < 0) IP.Number_of_Initial_Mesh_Refinements = 0;

    } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Uniform_Mesh_Refinements") == 0) {
      i_command = 75;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Number_of_Uniform_Mesh_Refinements;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Number_of_Uniform_Mesh_Refinements < 0) IP.Number_of_Uniform_Mesh_Refinements = 0;

    } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Boundary_Mesh_Refinements") == 0) {
      i_command = 76;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Number_of_Boundary_Mesh_Refinements;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Number_of_Boundary_Mesh_Refinements < 0) IP.Number_of_Boundary_Mesh_Refinements = 0;

    } else if (strcmp(IP.Next_Control_Parameter,"Maximum_Refinement_Level") == 0) {
      i_command = 77;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Maximum_Refinement_Level;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Maximum_Refinement_Level < 1) IP.Maximum_Refinement_Level = 1;

    } else if (strcmp(IP.Next_Control_Parameter,"Minimum_Refinement_Level") == 0) {
      i_command = 78;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Minimum_Refinement_Level;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Minimum_Refinement_Level < 1) IP.Minimum_Refinement_Level = 1;

    } else if (strcmp(IP.Next_Control_Parameter, "Threshold_for_Refinement") == 0) {
       i_command = 79;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Threshold_for_Refinement;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Threshold_for_Refinement <= ZERO ||
           IP.Threshold_for_Refinement > ONE) IP.Threshold_for_Refinement = 0.50;

    } else if (strcmp(IP.Next_Control_Parameter, "Threshold_for_Coarsening") == 0) {
      i_command = 80;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Threshold_for_Coarsening;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Threshold_for_Coarsening < ZERO ||
	  IP.Threshold_for_Coarsening >= ONE) IP.Threshold_for_Coarsening = 0.10;
      
    } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Refinement_Criteria") == 0) {
      i_command = 81;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Number_of_Refinement_Criteria;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Number_of_Refinement_Criteria < 1 || IP.Number_of_Refinement_Criteria > 6) {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Gradient_Density") == 0) {
      i_command = 82;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	IP.Refinement_Criteria_Gradient_Density = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	IP.Refinement_Criteria_Gradient_Density = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Divergence_Velocity") == 0) {
      i_command = 83;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	IP.Refinement_Criteria_Divergence_Velocity = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	IP.Refinement_Criteria_Divergence_Velocity = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Curl_Velocity") == 0) {
      i_command = 84;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	IP.Refinement_Criteria_Curl_Velocity = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	IP.Refinement_Criteria_Curl_Velocity = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Gradient_Temperature") == 0) {
      i_command = 85;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	IP.Refinement_Criteria_Gradient_Temperature = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	IP.Refinement_Criteria_Gradient_Temperature = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Gradient_CH4") == 0) {
      i_command = 86;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	IP.Refinement_Criteria_Gradient_CH4 = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	IP.Refinement_Criteria_Gradient_CH4 = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Gradient_CO2") == 0) {
      i_command = 87;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	IP.Refinement_Criteria_Gradient_CO2 = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	IP.Refinement_Criteria_Gradient_CO2 = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Residual_Variable") == 0) {
      i_command = 88;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.i_Residual_Variable;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.i_Residual_Variable < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Residual_Smoothing_Epsilon") == 0) {
      i_command = 89;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Residual_Smoothing_Epsilon;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Residual_Smoothing_Epsilon <= ZERO) {
         IP.Residual_Smoothing = 0;
         IP.Residual_Smoothing_Epsilon = ZERO;
      } else {
         IP.Residual_Smoothing = 1;
      } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter,"Mesh_Stretching") == 0) {
      i_command = 101;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	IP.i_Mesh_Stretching = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	IP.i_Mesh_Stretching = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Mesh_Stretching_Type_Idir") == 0) {
      i_command = 102;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"Linear") == 0) {
	IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_LINEAR;
      } else if (strcmp(IP.Next_Control_Parameter,"Min_Clustering") == 0) {
	IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MIN_CLUSTERING;
      } else if (strcmp(IP.Next_Control_Parameter,"Max_Clustering") == 0) {
	IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MAX_CLUSTERING;
      } else if (strcmp(IP.Next_Control_Parameter,"MinMax_Clustering") == 0) {
	IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
      } else if (strcmp(IP.Next_Control_Parameter,"Midpt_Clustering") == 0) {
	IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MIDPT_CLUSTERING;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Mesh_Stretching_Type_Jdir") == 0) {
      i_command = 103;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"Linear") == 0) {
	IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_LINEAR;
      } else if (strcmp(IP.Next_Control_Parameter,"Min_Clustering") == 0) {
	IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
      } else if (strcmp(IP.Next_Control_Parameter,"Max_Clustering") == 0) {
	IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MAX_CLUSTERING;
      } else if (strcmp(IP.Next_Control_Parameter,"MinMax_Clustering") == 0) {
	IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MINMAX_CLUSTERING;
      } else if (strcmp(IP.Next_Control_Parameter,"Midpt_Clustering") == 0) {
	IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIDPT_CLUSTERING;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Mesh_Stretching_Factor_Idir") == 0) {
      i_command = 104;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Mesh_Stretching_Factor_Idir;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Mesh_Stretching_Factor_Idir < ONE) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Mesh_Stretching_Factor_Jdir") == 0) {
      i_command = 105;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Mesh_Stretching_Factor_Jdir;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Mesh_Stretching_Factor_Jdir < ONE) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Smooth_Quad_Block") == 0) {
      i_command = 106;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	IP.i_Smooth_Quad_Block = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	IP.i_Smooth_Quad_Block = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Turbulence_BC_Type") == 0) {
      i_command = 161;
      Get_Next_Input_Control_Parameter(IP);;
      strcpy(IP.Turbulence_BC_Type,IP.Next_Control_Parameter);
      if (strcmp(IP.Turbulence_BC_Type,"Direct_Integration") == 0) {
	IP.i_Turbulence_BCs = TURBULENT_BC_DIRECT_INTEGRATION;
      } else if (strcmp(IP.Turbulence_BC_Type,"Standard_Wall_Function") == 0) {
	IP.i_Turbulence_BCs = TURBULENT_BC_STANDARD_WALL_FUNCTION;
      } else if (strcmp(IP.Turbulence_BC_Type,"Two_Layer_Wall_Function") == 0) {
	IP.i_Turbulence_BCs = TURBULENT_BC_TWO_LAYER_WALL_FUNCTION;
      } else if (strcmp(IP.Turbulence_BC_Type,"Three_Layer_Wall_Function") == 0) {
	IP.i_Turbulence_BCs = TURBULENT_BC_THREE_LAYER_WALL_FUNCTION;
      } else if (strcmp(IP.Turbulence_BC_Type,"Automatic_Wall_Treatment") == 0) {
	IP.i_Turbulence_BCs = TURBULENT_BC_AUTOMATIC_WALL_TREATMENT;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Friction_Velocity_Type") == 0) {
      i_command = 162;
      Get_Next_Input_Control_Parameter(IP);;
      strcpy(IP.Friction_Velocity_Type,IP.Next_Control_Parameter);
      if (strcmp(IP.Friction_Velocity_Type,"Local") == 0) {
	IP.i_Friction_Velocity = FRICTION_VELOCITY_LOCAL_SHEAR_STRESS;
      } else if (strcmp(IP.Friction_Velocity_Type,"Wall") == 0) {
	IP.i_Friction_Velocity = FRICTION_VELOCITY_WALL_SHEAR_STRESS;
      } else if (strcmp(IP.Friction_Velocity_Type,"Iterative") == 0) {
	IP.i_Friction_Velocity = FRICTION_VELOCITY_ITERATIVE;
      } else if (strcmp(IP.Friction_Velocity_Type,"Pipe") == 0) {
	IP.i_Friction_Velocity = FRICTION_VELOCITY_PIPE;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"C_constant") == 0) {
      i_command = 163;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.C_constant;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.C_constant < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"von_Karman_Constant") == 0) {
      i_command = 164;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.von_Karman_Constant;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.von_Karman_Constant < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"yplus_sublayer") == 0) {
      i_command = 165;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.yplus_sublayer;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.yplus_sublayer < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"yplus_buffer_layer") == 0) {
      i_command = 166;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.yplus_buffer_layer;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.yplus_buffer_layer < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"yplus_outer_layer") == 0) {
      i_command = 167;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.yplus_outer_layer;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.yplus_outer_layer < ZERO) i_command = INVALID_INPUT_VALUE;

      ////////////////////////////////////////////////////////////////////
      // FAS MULTIGRID PARAMETERS                                       //
      ////////////////////////////////////////////////////////////////////

    } else if (strcmp(IP.Next_Control_Parameter,"Multigrid_Levels") == 0) {
      i_command = 201;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Multigrid_IP.Levels;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Multigrid_IP.Levels <= ONE) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Multigrid_Cycle_Type") == 0) {
      i_command = 202;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Multigrid_IP.Cycle_Type,IP.Next_Control_Parameter);
      if (strcmp(IP.Multigrid_IP.Cycle_Type,"V") == 0 ||
	  strcmp(IP.Multigrid_IP.Cycle_Type,"v") == 0) {
	IP.Multigrid_IP.i_Cycle = MULTIGRID_V_CYCLE;
      } else if (strcmp(IP.Multigrid_IP.Cycle_Type,"W") == 0 ||
		 strcmp(IP.Multigrid_IP.Cycle_Type,"w") == 0) {
	IP.Multigrid_IP.i_Cycle = MULTIGRID_W_CYCLE;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Full_Multigrid") == 0) {
      i_command = 203;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Multigrid_IP.Number_of_Cycles_per_Stage_for_Full_Multigrid;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Multigrid_IP.Number_of_Cycles_per_Stage_for_Full_Multigrid < 0) 
	IP.Multigrid_IP.Number_of_Cycles_per_Stage_for_Full_Multigrid = 0;

    } else if (strcmp(IP.Next_Control_Parameter,"Defect_Correction") == 0) {
      i_command = 204;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	IP.Multigrid_IP.Defect_Correction = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	IP.Multigrid_IP.Defect_Correction = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"First_Order_Coarse_Mesh_Reconstruction") == 0) {
      i_command = 205;
      IP.Line_Number = IP.Line_Number + 1;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	IP.Multigrid_IP.First_Order_Coarse_Mesh_Reconstruction = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	IP.Multigrid_IP.First_Order_Coarse_Mesh_Reconstruction = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Prolong_Using_Injection") == 0) {
      i_command = 206;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	IP.Multigrid_IP.Prolong_Using_Injection = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	IP.Multigrid_IP.Prolong_Using_Injection = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Apply_Coarse_Mesh_Boundary_Conditions") == 0) {
      i_command = 207;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	IP.Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	IP.Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Injection_at_Dirichlet_Boundary_Conditions") == 0) {
      i_command = 208;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	IP.Multigrid_IP.Injection_at_Dirichlet_Boundary_Conditions = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	IP.Multigrid_IP.Injection_at_Dirichlet_Boundary_Conditions = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Update_Stability_Switch") == 0) {
      i_command = 209;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	IP.Multigrid_IP.Update_Stability_Switch = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	IP.Multigrid_IP.Update_Stability_Switch = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Maximum_Number_of_Update_Reductions") == 0) {
      i_command = 210;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Multigrid_IP.Maximum_Number_of_Update_Reductions;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Multigrid_IP.Maximum_Number_of_Update_Reductions < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Multigrid_Number_of_Smooths_on_Finest_Level") == 0) {
      i_command = 211;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Multigrid_IP.Number_of_Smooths_on_Finest_Level;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Multigrid_IP.Number_of_Smooths_on_Finest_Level < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Multigrid_Number_of_Pre_Smooths") == 0) {
      i_command = 212;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Multigrid_IP.Number_of_Pre_Smooths;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Multigrid_IP.Number_of_Pre_Smooths < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Multigrid_Number_of_Post_Smooths") == 0) {
      i_command = 213;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Multigrid_IP.Number_of_Post_Smooths;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Multigrid_IP.Number_of_Post_Smooths < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Multigrid_Number_of_Smooths_on_Coarsest_Level") == 0) {
      i_command = 214;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Multigrid_IP.Number_of_Smooths_on_Coarsest_Level;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Multigrid_IP.Number_of_Smooths_on_Coarsest_Level < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Convergence_Residual_Level") == 0) {
      i_command = 215;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Multigrid_IP.Convergence_Residual_Level;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Multigrid_IP.Convergence_Residual_Level < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Multigrid_Smoothing_Type") == 0) {
      i_command = 216;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Multigrid_IP.Smoothing_Type,IP.Next_Control_Parameter);
      if (strcmp(IP.Multigrid_IP.Smoothing_Type,"Explicit_Euler") == 0) {
	IP.Multigrid_IP.i_Smoothing = TIME_STEPPING_EXPLICIT_EULER;
	IP.N_Stage = 1;
      } else if (strcmp(IP.Multigrid_IP.Smoothing_Type,"Explicit_Predictor_Corrector") == 0) {
	IP.Multigrid_IP.i_Smoothing = TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR;
	IP.N_Stage = 2;
      } else if (strcmp(IP.Multigrid_IP.Smoothing_Type,"Explicit_Runge_Kutta") == 0) {
	IP.Multigrid_IP.i_Smoothing = TIME_STEPPING_EXPLICIT_RUNGE_KUTTA;
	IP.N_Stage = 5;
      } else if (strcmp(IP.Multigrid_IP.Smoothing_Type,"Multistage_Optimal_Smoothing") == 0) {
	IP.Multigrid_IP.i_Smoothing = TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING;
	IP.N_Stage = 4;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Ncycles_Regular_Multigrid") == 0) {
      i_command = 217;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Multigrid_IP.Ncycles_Regular_Multigrid;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Multigrid_IP.Ncycles_Regular_Multigrid < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Ncycles_Full_Multigrid") == 0) {
      i_command = 218;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Multigrid_IP.Ncycles_Full_Multigrid;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Multigrid_IP.Ncycles_Full_Multigrid < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Physical_Time_Integration_Type") == 0) {
      i_command = 219;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Multigrid_IP.Physical_Time_Integration_Type,
	     IP.Next_Control_Parameter);
      if (strcmp(IP.Multigrid_IP.Physical_Time_Integration_Type,"Implicit_Euler") == 0) {
	IP.Multigrid_IP.i_Physical_Time_Integration = TIME_STEPPING_IMPLICIT_EULER;
      } else if (strcmp(IP.Multigrid_IP.Physical_Time_Integration_Type,"Implicit_Trapezoidal") == 0) {
	IP.Multigrid_IP.i_Physical_Time_Integration = TIME_STEPPING_IMPLICIT_TRAPEZOIDAL;
	i_command = INVALID_INPUT_VALUE;
      } else if (strcmp(IP.Multigrid_IP.Physical_Time_Integration_Type,"Second_Order_Backwards") == 0) {
	IP.Multigrid_IP.i_Physical_Time_Integration = TIME_STEPPING_IMPLICIT_SECOND_ORDER_BACKWARD;
      } else if (strcmp(IP.Multigrid_IP.Physical_Time_Integration_Type,"Adams_Type") == 0) {
	IP.Multigrid_IP.i_Physical_Time_Integration = TIME_STEPPING_IMPLICIT_ADAMS_TYPE;
	i_command = INVALID_INPUT_VALUE;
      } else if (strcmp(IP.Multigrid_IP.Physical_Time_Integration_Type,"Lees") == 0) {
	IP.Multigrid_IP.i_Physical_Time_Integration = TIME_STEPPING_IMPLICIT_LEES;
	i_command = INVALID_INPUT_VALUE;
      } else if (strcmp(IP.Multigrid_IP.Physical_Time_Integration_Type,"Two_Step_Trapezoidal") == 0) {
	IP.Multigrid_IP.i_Physical_Time_Integration = TIME_STEPPING_IMPLICIT_TWO_STEP_TRAPEZOIDAL;
	i_command = INVALID_INPUT_VALUE;
      } else if (strcmp(IP.Multigrid_IP.Physical_Time_Integration_Type,"A_Contractive") == 0) {
	IP.Multigrid_IP.i_Physical_Time_Integration = TIME_STEPPING_IMPLICIT_A_CONTRACTIVE;
	i_command = INVALID_INPUT_VALUE;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Physical_Time_CFL_Number") == 0) {
      i_command = 220;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Multigrid_IP.Physical_Time_CFL_Number;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Multigrid_IP.Physical_Time_CFL_Number <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Dual_Time_Convergence_Residual_Level") == 0) {
      i_command = 221;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Multigrid_IP.Dual_Time_Convergence_Residual_Level;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Multigrid_IP.Dual_Time_Convergence_Residual_Level < ZERO) i_command = INVALID_INPUT_VALUE;

       ////////////////////////////////////////////////////////////////////
       // SPECIFIED BOUNDARY CONDITIONS                                  //
       ////////////////////////////////////////////////////////////////////

    } else if (strcmp(IP.Next_Control_Parameter,"Boundary_Conditions_Specified") == 0) {
      i_command = 500;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Boundary_Conditions_Specified,IP.Next_Control_Parameter);
      if (strcmp(IP.Boundary_Conditions_Specified,"On") == 0) {
	IP.BCs_Specified = ON;
      } else if (strcmp(IP.Boundary_Conditions_Specified,"Off") == 0) {
	IP.BCs_Specified = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"BC_North") == 0) {
      i_command = 501;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.BC_North_Type,IP.Next_Control_Parameter);
      if (strcmp(IP.BC_North_Type,"Reflection") == 0) {
	IP.BC_North = BC_REFLECTION;
      } else if (strcmp(IP.BC_North_Type,"Burning_Surface") == 0) {
	IP.BC_North = BC_BURNING_SURFACE;
      } else if (strcmp(IP.BC_North_Type,"Wall_Viscous_Isothermal") == 0) {
	IP.BC_North = BC_WALL_VISCOUS_ISOTHERMAL;
      } else if (strcmp(IP.BC_North_Type,"Wall_Viscous_Heatflux") == 0) {
	IP.BC_North = BC_WALL_VISCOUS_HEATFLUX;
      } else if (strcmp(IP.BC_North_Type,"Moving_Wall_Isothermal") == 0) {
	IP.BC_North = BC_MOVING_WALL_ISOTHERMAL;
      } else if (strcmp(IP.BC_North_Type,"Moving_Wall_Heatflux") == 0) {
	IP.BC_North = BC_MOVING_WALL_HEATFLUX;
      } else if (strcmp(IP.BC_North_Type,"Fixed") == 0) {
	IP.BC_North = BC_FIXED;
      } else if (strcmp(IP.BC_North_Type,"Inflow_Subsonic") == 0) {
	IP.BC_North = BC_INFLOW_SUBSONIC;
      } else if (strcmp(IP.BC_North_Type,"Outflow_Subsonic") == 0) {
	IP.BC_North = BC_OUTFLOW_SUBSONIC;
      } else if (strcmp(IP.BC_North_Type,"Fixed_Pressure") == 0) {
	IP.BC_North = BC_FIXED_PRESSURE;
      } else if (strcmp(IP.BC_North_Type,"Constant_Extrapolation") == 0) {
	IP.BC_North = BC_CONSTANT_EXTRAPOLATION;
      } else if (strcmp(IP.BC_North_Type,"Linear_Extrapolation") == 0) {
	IP.BC_North = BC_LINEAR_EXTRAPOLATION;
      } else if (strcmp(IP.BC_North_Type,"Characteristic") == 0) {
	IP.BC_North = BC_CHARACTERISTIC;
      } else if (strcmp(IP.BC_North_Type,"None") == 0) {
	IP.BC_North = BC_NONE;
      } else if (strcmp(IP.BC_North_Type,"Ringleb") == 0) {
	IP.BC_North = BC_RINGLEB_FLOW;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"BC_South") == 0) {
      i_command = 502;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.BC_South_Type,IP.Next_Control_Parameter);
      if (strcmp(IP.BC_South_Type,"Reflection") == 0) {
	IP.BC_South = BC_REFLECTION;
      } else if (strcmp(IP.BC_South_Type,"Burning_Surface") == 0) {
	IP.BC_South = BC_BURNING_SURFACE;
      } else if (strcmp(IP.BC_South_Type,"Wall_Viscous_Isothermal") == 0) {
	IP.BC_South = BC_WALL_VISCOUS_ISOTHERMAL;
      } else if (strcmp(IP.BC_South_Type,"Wall_Viscous_Heatflux") == 0) {
	IP.BC_South = BC_WALL_VISCOUS_HEATFLUX;
      } else if (strcmp(IP.BC_South_Type,"Moving_Wall_Isothermal") == 0) {
	IP.BC_South = BC_MOVING_WALL_ISOTHERMAL;
      } else if (strcmp(IP.BC_South_Type,"Moving_Wall_Heatflux") == 0) {
	IP.BC_South = BC_MOVING_WALL_HEATFLUX;
      } else if (strcmp(IP.BC_South_Type,"Fixed") == 0) {
	IP.BC_South = BC_FIXED;
      } else if (strcmp(IP.BC_South_Type,"Inflow_Subsonic") == 0) {
	IP.BC_South = BC_INFLOW_SUBSONIC;
      } else if (strcmp(IP.BC_South_Type,"Outflow_Subsonic") == 0) {
	IP.BC_South = BC_OUTFLOW_SUBSONIC;
      } else if (strcmp(IP.BC_South_Type,"Fixed_Pressure") == 0) {
	IP.BC_South = BC_FIXED_PRESSURE;
      } else if (strcmp(IP.BC_South_Type,"Constant_Extrapolation") == 0) {
	IP.BC_South = BC_CONSTANT_EXTRAPOLATION;
      } else if (strcmp(IP.BC_South_Type,"Linear_Extrapolation") == 0) {
	IP.BC_South = BC_LINEAR_EXTRAPOLATION;
      } else if (strcmp(IP.BC_South_Type,"Characteristic") == 0) {
	IP.BC_South = BC_CHARACTERISTIC;
      } else if (strcmp(IP.BC_South_Type,"None") == 0) {
	IP.BC_South = BC_NONE;
      } else if (strcmp(IP.BC_South_Type,"Ringleb") == 0) {
	IP.BC_South = BC_RINGLEB_FLOW;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"BC_East") == 0) {
      i_command = 503;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.BC_East_Type,IP.Next_Control_Parameter);
      if (strcmp(IP.BC_East_Type,"Reflection") == 0) {
	IP.BC_East = BC_REFLECTION;
      } else if (strcmp(IP.BC_East_Type,"Burning_Surface") == 0) {
	IP.BC_East = BC_BURNING_SURFACE;
      } else if (strcmp(IP.BC_East_Type,"Wall_Viscous_Isothermal") == 0) {
	IP.BC_East = BC_WALL_VISCOUS_ISOTHERMAL;
      } else if (strcmp(IP.BC_East_Type,"Wall_Viscous_Heatflux") == 0) {
	IP.BC_East = BC_WALL_VISCOUS_HEATFLUX;
      } else if (strcmp(IP.BC_East_Type,"Moving_Wall_Isothermal") == 0) {
	IP.BC_East = BC_MOVING_WALL_ISOTHERMAL;
      } else if (strcmp(IP.BC_East_Type,"Moving_Wall_Heatflux") == 0) {
	IP.BC_East = BC_MOVING_WALL_HEATFLUX;
      } else if (strcmp(IP.BC_East_Type,"Fixed") == 0) {
	IP.BC_East = BC_FIXED;
      } else if (strcmp(IP.BC_East_Type,"Inflow_Subsonic") == 0) {
	IP.BC_East = BC_INFLOW_SUBSONIC;
      } else if (strcmp(IP.BC_East_Type,"Outflow_Subsonic") == 0) {
	IP.BC_East = BC_OUTFLOW_SUBSONIC;
      } else if (strcmp(IP.BC_East_Type,"Fixed_Pressure") == 0) {
	IP.BC_East = BC_FIXED_PRESSURE;
      } else if (strcmp(IP.BC_East_Type,"Constant_Extrapolation") == 0) {
	IP.BC_East = BC_CONSTANT_EXTRAPOLATION;
      } else if (strcmp(IP.BC_East_Type,"Linear_Extrapolation") == 0) {
	IP.BC_East = BC_LINEAR_EXTRAPOLATION;
      } else if (strcmp(IP.BC_East_Type,"Characteristic") == 0) {
	IP.BC_East = BC_CHARACTERISTIC;
      } else if (strcmp(IP.BC_East_Type,"None") == 0) {
	IP.BC_East = BC_NONE;
      } else if (strcmp(IP.BC_East_Type,"Ringleb") == 0) {
	IP.BC_East = BC_RINGLEB_FLOW;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"BC_West") == 0) {
      i_command = 504;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.BC_West_Type,IP.Next_Control_Parameter);
      if (strcmp(IP.BC_West_Type,"Reflection") == 0) {
	IP.BC_West = BC_REFLECTION;
      } else if (strcmp(IP.BC_West_Type,"Burning_Surface") == 0) {
	IP.BC_West = BC_BURNING_SURFACE;
      } else if (strcmp(IP.BC_West_Type,"Wall_Viscous_Isothermal") == 0) {
	IP.BC_West = BC_WALL_VISCOUS_ISOTHERMAL;
      } else if (strcmp(IP.BC_West_Type,"Wall_Viscous_Heatflux") == 0) {
	IP.BC_West = BC_WALL_VISCOUS_HEATFLUX;
      } else if (strcmp(IP.BC_West_Type,"Moving_Wall_Isothermal") == 0) {
	IP.BC_West = BC_MOVING_WALL_ISOTHERMAL;
      } else if (strcmp(IP.BC_West_Type,"Moving_Wall_Heatflux") == 0) {
	IP.BC_West = BC_MOVING_WALL_HEATFLUX;
      } else if (strcmp(IP.BC_West_Type,"Fixed") == 0) {
	IP.BC_West = BC_FIXED;
      } else if (strcmp(IP.BC_West_Type,"Inflow_Subsonic") == 0) {
	IP.BC_West = BC_INFLOW_SUBSONIC;
      } else if (strcmp(IP.BC_West_Type,"Outflow_Subsonic") == 0) {
	IP.BC_West = BC_OUTFLOW_SUBSONIC;
      } else if (strcmp(IP.BC_West_Type,"Fixed_Pressure") == 0) {
	IP.BC_West = BC_FIXED_PRESSURE;
      } else if (strcmp(IP.BC_West_Type,"Constant_Extrapolation") == 0) {
	IP.BC_West = BC_CONSTANT_EXTRAPOLATION;
      } else if (strcmp(IP.BC_West_Type,"Linear_Extrapolation") == 0) {
	IP.BC_West = BC_LINEAR_EXTRAPOLATION;
      } else if (strcmp(IP.BC_West_Type,"Characteristic") == 0) {
	IP.BC_West = BC_CHARACTERISTIC;
      } else if (strcmp(IP.BC_West_Type,"None") == 0) {
	IP.BC_West = BC_NONE;
      } else if (strcmp(IP.BC_West_Type,"Ringleb") == 0) {
	IP.BC_West = BC_RINGLEB_FLOW;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter, "Execute") == 0) {
       i_command = EXECUTE_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Terminate") == 0) {
       i_command = TERMINATE_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Continue") == 0) {
       i_command = CONTINUE_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output") == 0) {
       i_command = WRITE_OUTPUT_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Cells") == 0) {
       i_command = WRITE_OUTPUT_CELLS_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Nodes") == 0) {
       i_command = WRITE_OUTPUT_NODES_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Restart") == 0) {
       i_command = WRITE_RESTART_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_RHS") == 0) {
       i_command = WRITE_OUTPUT_RHS_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Mesh") == 0) {
       i_command = WRITE_OUTPUT_GRID_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Mesh_Definition") == 0) {
       i_command = WRITE_GRID_DEFINITION_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Mesh_Nodes") == 0) {
       i_command = WRITE_OUTPUT_GRID_NODES_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Mesh_Cells") == 0) {
       i_command = WRITE_OUTPUT_GRID_CELLS_CODE;

    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Ringleb") == 0) {
      i_command = WRITE_OUTPUT_RINGLEB_CODE;
    
    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Viscous_Channel") == 0) {
      i_command = WRITE_OUTPUT_VISCOUS_CHANNEL_CODE;
        
    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Flat_Plate") == 0) {
      i_command = WRITE_OUTPUT_FLAT_PLATE_CODE;
    
    } else if (strcmp(IP.Next_Control_Parameter, "Refine_Grid") == 0) {
       i_command = REFINE_GRID_CODE;

    } else if (IP.Next_Control_Parameter[0] == '#') {
       i_command = COMMENT_CODE;
    } else {
       i_command = INVALID_INPUT_CODE;
       cout<<"\n ERROR ERROR WILL ROBINSON.... INVALID_INPUT_CODE \n";
    } /* endif */

   
    /* Return the parser command type indicator. */
    return (i_command);

}

/********************************************************
 * Routine: Process_Input_Control_Parameter_File        *
 *                                                      *
 * Reads, parses, and executes the list of input        *
 * control parameters from the standard input file.     *
 *                                                      *
 ********************************************************/
int Process_Input_Control_Parameter_File(Chem2D_Input_Parameters &Input_Parameters,
                                         char *Input_File_Name_ptr,
                                         int &Command_Flag) {

    int error_flag, line_number;

    /* Assign initial value for error indicator flag. */
    error_flag = 0;

    /* Assign default values to the input parameters. */
    Set_Default_Input_Parameters(Input_Parameters);

   /* Copy input file name (a string) to appropriate input parameter variable. */
    if (Input_File_Name_ptr != NULL) strcpy(Input_Parameters.Input_File_Name, Input_File_Name_ptr);

    /* Open the input file containing the input parameters. */
    Open_Input_File(Input_Parameters);
    error_flag = Input_Parameters.Input_File.bad();

    if (error_flag) {
       cout << "\n Chem2D ERROR: Unable to open Chem2D input data file.\n";
       return (error_flag);
    } /* endif */

    /* Read and parse control parameters contained in
       the input file. */

    while (1) {
       Get_Next_Input_Control_Parameter(Input_Parameters);
    
       Command_Flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
       line_number = Input_Parameters.Line_Number;
       
   
       if (Command_Flag == EXECUTE_CODE) {
          break;
       } else if (Command_Flag == TERMINATE_CODE) {
          break;
       } else if (Command_Flag == INVALID_INPUT_CODE ||
                  Command_Flag == INVALID_INPUT_VALUE) {
          line_number = -line_number;
          cout << "\n Chem2D ERROR: Error reading Chem2D data at line #"
               << -line_number  << " of input data file.\n";
          error_flag = line_number;
          break;
       } /* endif */
    } /* endwhile */

    /* Initial processing of input control parameters complete.  
       Return the error indicator flag. */
   
    //Load the C-Type strings from the C++ strings 
    strcpy(Input_Parameters.React_Name,Input_Parameters.react_name.c_str());
    
    for (int i = 0; i < Input_Parameters.num_species; i++) {
      strcpy(Input_Parameters.Multispecies[i],Input_Parameters.multispecies[i].c_str());
    }

    // Proper temperature for display
    Input_Parameters.Temperature = Input_Parameters.Wo.T();
 
    // Reset the static variables.
    Input_Parameters.Wo.set_flow_type(Input_Parameters.FlowType);
    Input_Parameters.Uo.set_flow_type(Input_Parameters.FlowType);
    Input_Parameters.Wo.set_turbulence_variables(Input_Parameters.C_constant,
						 Input_Parameters.von_Karman_Constant,
						 Input_Parameters.yplus_sublayer,
						 Input_Parameters.yplus_buffer_layer,
						 Input_Parameters.yplus_outer_layer);
    Input_Parameters.Uo.set_turbulence_variables(Input_Parameters.C_constant,
						 Input_Parameters.von_Karman_Constant,
						 Input_Parameters.yplus_sublayer,
						 Input_Parameters.yplus_buffer_layer,
						 Input_Parameters.yplus_outer_layer);

    // Perform consitency checks on the refinement criteria.
    Input_Parameters.Number_of_Refinement_Criteria = 0;
    if (Input_Parameters.Refinement_Criteria_Gradient_Density) Input_Parameters.Number_of_Refinement_Criteria++;
    if (Input_Parameters.Refinement_Criteria_Divergence_Velocity) Input_Parameters.Number_of_Refinement_Criteria++;
    if (Input_Parameters.Refinement_Criteria_Curl_Velocity) Input_Parameters.Number_of_Refinement_Criteria++;
    if (Input_Parameters.Refinement_Criteria_Gradient_Temperature) Input_Parameters.Number_of_Refinement_Criteria++;
    if (Input_Parameters.Refinement_Criteria_Gradient_CH4) Input_Parameters.Number_of_Refinement_Criteria++;
    if (Input_Parameters.Refinement_Criteria_Gradient_CO2) Input_Parameters.Number_of_Refinement_Criteria++;
    if (Input_Parameters.Number_of_Refinement_Criteria < 1 || Input_Parameters.Number_of_Refinement_Criteria > 6) return 1011;

    // Initial processing of input control parameters complete.  
    return (error_flag);

}

/********************************************************
 * Routine: Equivalence_Ratio                           *
 * 
 ********************************************************/
void Equivalence_Ratio(const double &phi){

  // Methane & Air
  // CH4 + 2(O2 +3.76N2) -> CO2 + 2H2O 
  double n = 9.52*(1.0 - phi)/(phi +9.52);

  double X_CH4 = ( 1.0 - n)/10.52;
  double X_O2  = (2.0 + 0.21*n)/10.52;
  double X_N2  = (7.52 +0.79*n)/10.52;

  double Mtot = X_CH4*(16.0) + X_O2*(32.0) + X_N2*(28.0);

  double c_CH4 =  X_CH4*(16.0)/Mtot;
  double c_O2  =  X_O2*(32.0)/Mtot;
  double c_N2  =  X_N2*(28.0)/Mtot;

  cout<<"\n For phi = "<<phi;
  cout<<"\n c_CH4 "<<c_CH4;
  cout<<"\n c_O2  "<<c_O2;
  cout<<"\n c_N2  "<<c_N2;
  cout<<"\n X_CH4 "<<X_CH4;
  cout<<"\n X_O2  "<<X_O2;
  cout<<"\n X_N2  "<<X_N2;
  cout<<"\n Mtot  "<<Mtot;

  cout<<endl;
  cout<<"\n Should pass in equation type, and pass back the appropriate mass fractions ";
  cout<<endl;
  
  exit(0); 

}
