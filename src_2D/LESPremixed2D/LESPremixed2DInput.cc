/*****************************************************//**
 * \file LESPremixed2DInput.cc
 *
 * Subroutines for the LESPremixed2DInput class.
 *
 **********************************************************/

// LESPREMIXED2D Header file
#ifndef _LESPREMIXED2D_INPUT_INCLUDED
#include "LESPremixed2DInput.h"
#endif // _LESPREMIXED2D_INPUT_INCLUDED

/*************************************************************
 * LESPremixed2D_Input_Parameters -- External subroutines.   *
 *************************************************************/

/********************************************************
 * Routine: Open_Input_File                             *
 *                                                      *
 * Opens the appropriate input data file.               *
 *                                                      *
 ********************************************************/
void Open_Input_File(LESPremixed2D_Input_Parameters &IP) {

    IP.Input_File.open(IP.Input_File_Name, ios::in);
    if (!IP.Input_File.fail()) {
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
void Close_Input_File(LESPremixed2D_Input_Parameters &IP) {

    IP.Input_File.unsetf(ios::skipws);
    IP.Input_File.close();

}

/********************************************************
 * Routine: Set_Default_Input_Parameters                *
 *                                                      *
 * Assigns default values to the input parameters.      *
 *                                                      *
 ********************************************************/
void Set_Default_Input_Parameters(LESPremixed2D_Input_Parameters &IP) {

    int i;
    char *string_ptr;

    string_ptr = "LESPremixed2D.in";
    strcpy(IP.Input_File_Name, string_ptr);

    // Time-stepping parameters:
    string_ptr = "Explicit_Euler";
    strcpy(IP.Time_Integration_Type, string_ptr);
    IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
    IP.Time_Accurate = 0;
    IP.Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
    IP.Preconditioning = 0; //default off
    IP.Dual_Time_Stepping = 0; //default off
    IP.Maximum_Number_of_Time_Steps = 100;
    IP.N_Stage = 1;
    IP.CFL_Number = 0.5;
    IP.Time_Max = ZERO;
    
    /* Dual time stepping */
    IP.Physical_CFL_Number = 0.9;
    IP.dTime = ZERO;
    IP.Max_Inner_Steps = 1;
    IP.first_step = 0;
    
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
    IP.get_cffc_path();

    /* Debug Level */
    IP.debug_level = 0;  //default no debug information

    IP.Mach_Number = ZERO;
    IP.Mach_Number_Reference = ONE;  //for precondtioning
    IP.Mach_Number_Reference_target = ONE;
    IP.Flow_Angle = ZERO; 
    IP.Global_Schmidt = ONE;
    IP.Reynolds_Number = 500000.0;
    IP.Kinematic_Viscosity_Wall = 1.4590E-5;
    IP.i_Grid_Level = 0;
    
    /*************************************************/
    /********** LESPREMIXED2D SPECIFIC ***************/
    //define multispecies with no reactions.
    IP.react_name ="NO_REACTIONS";   
    //define multispecies with no scalars
    IP.scalar_system_name = "LAMINAR";   
    
    IP.Fresh_Fuel_Mass_Fraction = 0.0551;
    IP.Burnt_Fuel_Mass_Fraction = ZERO;
    IP.Fresh_Density = 1.1301;
    IP.TKEo = ZERO;
    //Use air with 79% N2, and 21% 02 by volume.(ie. mol)
    IP.num_species = 2;
    IP.num_scalars = 0;
    // Energy spectrum for turbulence
    string_ptr = "von-Karman-Pao";
    strcpy(IP.Spectrum_Type, string_ptr);
    IP.i_Spectrum = VON_KARMAN_PAO;
    IP.Rescale_Spectrum = 0;
    IP.Read_Fluctuations_From_File = 0;

    // transport data 
    string_ptr = "Transport-NASA";
    strcpy(IP.trans_type, string_ptr);
    IP.i_trans_type = TRANSPORT_NASA;

    IP.Allocate();
    IP.multispecies[0] = "N2"; 
    IP.multispecies[1] = "O2"; 
    IP.mass_fractions[0] = 0.765; 
    IP.mass_fractions[1] = 0.235;
    IP.Schmidt[0] = IP.Global_Schmidt;
    IP.Schmidt[1] = IP.Global_Schmidt;

    IP.Wo.React.set_reactions(IP.react_name); 
    IP.Wo.React.set_species(IP.multispecies,IP.num_species);
    IP.Wo.Scal_sys.scalar_set(IP.scalar_system_name);
   
    //Get Species parameters and set default initial values
    IP.Wo.set_species_data(IP.num_species, IP.num_scalars, IP.multispecies, 
			   IP.CFFC_Path,
			   IP.Mach_Number_Reference,
			   IP.Schmidt,
      			   IP.i_trans_type); 
    IP.Uo.set_species_data(IP.num_species, IP.num_scalars, IP.multispecies,
			   IP.CFFC_Path,
			   IP.Mach_Number_Reference,
			   IP.Schmidt,
      			   IP.i_trans_type);

    //Air at STD_ATM
    IP.Pressure = IP.Wo.p;
    IP.Temperature = IP.Wo.T(); 
    IP.Wo.v.x =IP.Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Flow_Angle/360.00);
    IP.Wo.v.y =IP.Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Flow_Angle/360.00);
    IP.Wo.set_initial_values(IP.mass_fractions);
    IP.Uo.set_initial_values(IP.mass_fractions);
    IP.Uo = U(IP.Wo);   
    IP.Heat_Source = ZERO;

    // cantera parameters
    IP.ct_mech_name = "none";
    IP.ct_mech_file = "none";

    IP.Smagorinsky_Constant = 0.18;
    IP.Yoshizawa_coefficient = 0.005;
    IP.filter_width = ZERO;

    IP.Wo.set_SFSmodel_variables(IP.filter_width,
				 IP.Smagorinsky_Constant,
				 IP.Yoshizawa_coefficient);
    IP.Uo.set_SFSmodel_variables(IP.filter_width,
				 IP.Smagorinsky_Constant,
				 IP.Yoshizawa_coefficient);

    IP.laminar_flame_thickness = 0.446E-3; // m
    IP.laminar_flame_speed = 0.38; // m/s
    IP.TFactor = ONE;
    IP.adiabatic_temp = 2000.0;
    IP.equivalence_ratio = ONE;
    IP.reactants_den = 1.13;

    IP.Wo.set_premixed_flame_variables(IP.laminar_flame_thickness,
				       IP.laminar_flame_speed,
				       IP.TFactor,
                                       IP.adiabatic_temp,
                                       IP.equivalence_ratio,
                                       IP.reactants_den);
    IP.Uo.set_premixed_flame_variables(IP.laminar_flame_thickness,
				       IP.laminar_flame_speed,
				       IP.TFactor,
                                       IP.adiabatic_temp,
                                       IP.equivalence_ratio,
                                       IP.reactants_den);

    /***** END LESPREMIXED2D SPECFIC *****************/
    /*************************************************/

    //BC
    IP.Moving_wall_velocity = ZERO; 
    IP.Re_lid = 100.0;
    IP.Pressure_Gradient = ZERO; 
    /* Flow type */
    string_ptr = "Inviscid";
    strcpy(IP.Flow_Type, string_ptr);
    IP.FlowType = FLOWTYPE_INVISCID;

    /* Flow geometry type */
    string_ptr = "Planar";
    strcpy(IP.Flow_Geometry_Type, string_ptr);
    IP.Axisymmetric = 0;
    IP.Gravity = 0;  //default sans gravity

    IP.BluffBody_Data_Usage = 0; 
    IP.Wall_Boundary_Treatments = 0;


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
    IP.Smooth_Bump = ON;

    IP.X_Shift = Vector2D_ZERO;
    IP.X_Scale = ONE;
    IP.X_Rotate = ZERO;
    IP.Length_Shroud = 0.1;
    IP.Radius_Shroud = 0.043;
    IP.Length_BluffBody = 0.04;
    IP.Radius_BluffBody = 0.02 ;
    IP.Radius_Orifice = 0.001;
    IP.Radius_Inlet_Pipe = 0.0508;
    IP.Radius_Combustor_Tube = 0.0762;
    IP.Length_Inlet_Pipe = 0.127;
    IP.Length_Combustor_Tube = 0.508;
    
    IP.BluffBody_Coflow_Air_Velocity = 20.0;
    IP.BluffBody_Coflow_Fuel_Velocity = 61.0; 
  
    // Boundary conditions:
    string_ptr = "OFF";
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
    IP.Number_of_Uniform_Mesh_Coarsenings = 0;
    IP.Number_of_Initial_Mesh_Refinements = 0;
    IP.Number_of_Uniform_Mesh_Refinements = 0;
    IP.Number_of_Boundary_Mesh_Refinements = 0;
    IP.Maximum_Refinement_Level = 100;
    IP.Minimum_Refinement_Level = 1;
    IP.Threshold_for_Refinement = 0.50;
    IP.Threshold_for_Coarsening = 0.10;
    IP.Number_of_Refinement_Criteria = 6;
    IP.Refinement_Criteria_Gradient_Density = ON;
    IP.Refinement_Criteria_Divergence_Velocity = OFF;
    IP.Refinement_Criteria_Curl_Velocity = OFF;
    IP.Refinement_Criteria_Gradient_Temperature = OFF;
    IP.Refinement_Criteria_Gradient_CH4 = OFF;
    IP.Refinement_Criteria_Gradient_CO2 = OFF;

    // Smooth quad block flag:
    IP.i_Smooth_Quad_Block = ON;

    // Default Solver Type
    IP.Solver_Type = EXPLICIT;
  
    IP.Morton = 0;
    IP.Morton_Reordering_Frequency = 1000000;
 
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

    //LESPREMIXED2D
    IP.Time_Accurate_Plot_Frequency = 0;

    string_ptr = " ";
    strcpy(IP.Next_Control_Parameter, string_ptr);

    IP.Line_Number = 0;
   
    IP.Number_of_Processors = CFFC_MPI::Number_of_Processors;
    IP.Number_of_Blocks_Per_Processor = 10;      

    // Freezing_Limiter
    IP.Freeze_Limiter = 0;
   
    IP.i_Residual_Variable = 1; //density

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
void Broadcast_Input_Parameters(LESPremixed2D_Input_Parameters &IP) {

#ifdef _MPI_VERSION

    MPI::COMM_WORLD.Bcast(IP.Input_File_Name, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Line_Number), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Time_Integration_Type, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
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
    MPI::COMM_WORLD.Bcast(&(IP.dTime), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Physical_CFL_Number), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Max_Inner_Steps), 
                          1, 
                          MPI::INT, 0);
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
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Reconstruction), 
                          1, 
                          MPI::INT, 0);   
    MPI::COMM_WORLD.Bcast(IP.Viscous_Flux_Evaluation_Type, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Viscous_Flux_Evaluation), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Limiter_Type, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Limiter), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Flux_Function_Type, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Flux_Function), 
                          1, 
                          MPI::INT, 0);
    // Initial conditions:
    MPI::COMM_WORLD.Bcast(IP.ICs_Type,
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_ICs), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Grid_Level), 
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
    MPI::COMM_WORLD.Bcast(&(IP.Mach_Number_Reference_target), 
			  1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Flow_Angle), 
                          1, 
                          MPI::DOUBLE, 0);   
    MPI::COMM_WORLD.Bcast(&(IP.Re_lid),
			  1,
			  MPI::DOUBLE,0);
    /*********************************************************
     ******************* LESPREMIXED2D SPECIFIC **************
     *********************************************************/
    MPI::COMM_WORLD.Bcast(IP.CFFC_Path, 
			  INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
			  MPI::CHAR, 0);
    //fresh and burnt fuel mass fraction
    MPI::COMM_WORLD.Bcast(&(IP.Fresh_Fuel_Mass_Fraction), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Burnt_Fuel_Mass_Fraction), 
                          1, 
                          MPI::DOUBLE, 0);
    //fresh density
    MPI::COMM_WORLD.Bcast(&(IP.Fresh_Density), 
                          1, 
                          MPI::DOUBLE, 0);
    //Initial turbulence kinetic energy
    MPI::COMM_WORLD.Bcast(&(IP.TKEo), 
                          1, 
                          MPI::DOUBLE, 0);
    //Spectrum flag
    MPI::COMM_WORLD.Bcast(IP.Spectrum_Type, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Spectrum), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Rescale_Spectrum), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Read_Fluctuations_From_File), 
                          1, 
                          MPI::INT, 0);
    //reaction name
    MPI::COMM_WORLD.Bcast(IP.React_Name, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
			  MPI::CHAR, 0);
    // cantera parameters
    MPI::COMM_WORLD.Bcast(IP.ct_Mech_Name,
                          INPUT_PARAMETER_LENGTH_CHEM2D,
			  MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.ct_Mech_File,
                          INPUT_PARAMETER_LENGTH_CHEM2D,
			  MPI::CHAR, 0);
    //scalar system name
    MPI::COMM_WORLD.Bcast(IP.Scalar_system_name, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
			  MPI::CHAR, 0);
    //delete current dynamic memory before changing num_species
    if(!CFFC_Primary_MPI_Processor()) {   
      IP.Deallocate();
    } 
    //number of species
    MPI::COMM_WORLD.Bcast(&(IP.num_species), 
                          1, 
                          MPI::INT, 0);
    //number of scalars
    MPI::COMM_WORLD.Bcast(&(IP.num_scalars), 
                          1, 
                          MPI::INT, 0);
    //transport data
    MPI::COMM_WORLD.Bcast(IP.trans_type, 
			  INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
			  MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_trans_type), 
			  1, 
			  MPI::INT, 0);


    //set up new dynamic memory
    if(!CFFC_Primary_MPI_Processor()) {   
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
			    INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
			    MPI::CHAR, 0);
    }
    //scalar names & values
    if(IP.num_scalars > 0){
      for(int i=0; i < IP.num_scalars; i++){
        MPI::COMM_WORLD.Bcast(&(IP.scalar_values[i]), 
			    1, 
			    MPI::DOUBLE, 0);
        MPI::COMM_WORLD.Bcast(IP.Scalars[i], 
			    INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
			    MPI::CHAR, 0);
      }
    }
    //set reaction and species parameters
    if (!CFFC_Primary_MPI_Processor()) {      
      IP.react_name = IP.React_Name;
      for (int i = 0; i < IP.num_species; i++) {
	IP.multispecies[i] = IP.Multispecies[i];  
      }
      IP.scalar_system_name = IP.Scalar_system_name;
      if(IP.num_scalars > 0){
        for(int i=0; i<IP.num_scalars; i++){
          IP.scalars[i] = IP.Scalars[i];
        }
      }

      //load reaction and scalar names
      if (IP.Wo.React.reactset_flag != CANTERA)
	IP.Wo.React.set_reactions(IP.react_name);
      else
	 IP.Wo.React.ct_load_mechanism(IP.ct_mech_file, IP.ct_mech_name);
      IP.Wo.Scal_sys.scalar_set(IP.scalar_system_name);
      
      //Set species if non-reacting
      if( IP.Wo.React.reactset_flag == NO_REACTIONS){
	IP.Wo.React.set_species(IP.multispecies,IP.num_species);
      }  
          
      //set the data for each
      IP.get_cffc_path();
      IP.Wo.set_species_data(IP.num_species, IP.num_scalars, IP.multispecies,
			     IP.CFFC_Path,
			     IP.Mach_Number_Reference,
			     IP.Schmidt,
			     IP.i_trans_type); 
      IP.Uo.set_species_data(IP.num_species, IP.num_scalars, IP.multispecies,
			     IP.CFFC_Path,
			     IP.Mach_Number_Reference,
			     IP.Schmidt,
			     IP.i_trans_type);
      IP.Wo.set_initial_values(IP.mass_fractions);
      IP.Uo.set_initial_values(IP.mass_fractions);
      if(IP.num_scalars > 0){
        IP.Wo.set_initial_values_scal(IP.scalar_values);   
        IP.Uo.set_initial_values_scal(IP.scalar_values);
      }
   
      //set proper Temp & Pressure instead of defaults
      IP.Wo.rho = IP.Pressure/(IP.Wo.Rtot()*IP.Temperature); 
      IP.Wo.p = IP.Pressure;	
      IP.Wo.v.zero();

      IP.Uo = U(IP.Wo);
    }
    
    // SFS modelling parameters
    MPI::COMM_WORLD.Bcast(&(IP.Smagorinsky_Constant), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Yoshizawa_coefficient), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.filter_width), 
                          1, 
                          MPI::DOUBLE, 0);
    // premixed flame parameters
    MPI::COMM_WORLD.Bcast(&(IP.laminar_flame_thickness), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.laminar_flame_speed), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.TFactor), 
                          1, 
                          MPI::DOUBLE, 0);  
    MPI::COMM_WORLD.Bcast(&(IP.adiabatic_temp), 
                          1, 
                          MPI::DOUBLE, 0);  
    MPI::COMM_WORLD.Bcast(&(IP.equivalence_ratio), 
                          1, 
                          MPI::DOUBLE, 0);  
    MPI::COMM_WORLD.Bcast(&(IP.reactants_den), 
                          1, 
                          MPI::DOUBLE, 0);  

    // Reset the static variables.
    IP.Wo.set_SFSmodel_variables(IP.filter_width,
				 IP.Smagorinsky_Constant,
				 IP.Yoshizawa_coefficient);
    IP.Uo.set_SFSmodel_variables(IP.filter_width,
				 IP.Smagorinsky_Constant,
				 IP.Yoshizawa_coefficient);
    IP.Wo.set_premixed_flame_variables(IP.laminar_flame_thickness,
				       IP.laminar_flame_speed,
				       IP.TFactor,
                                       IP.adiabatic_temp,
                                       IP.equivalence_ratio,
                                       IP.reactants_den);
    IP.Uo.set_premixed_flame_variables(IP.laminar_flame_thickness,
				       IP.laminar_flame_speed,
				       IP.TFactor,
                                       IP.adiabatic_temp,
                                       IP.equivalence_ratio,
                                       IP.reactants_den);
 
    /*********************************************************
     ******************* LESPREMIXED2D END *******************
     *********************************************************/


   if(!CFFC_Primary_MPI_Processor()) {   
      IP.Wo.v.x = IP.Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Flow_Angle/360.00);
      IP.Wo.v.y = IP.Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Flow_Angle/360.00);
    }

    /***************************************/
    MPI::COMM_WORLD.Bcast(IP.Flow_Type, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.FlowType), 
                          1, 
                          MPI::INT, 0);

   /***************************************/
    MPI::COMM_WORLD.Bcast(IP.Flow_Geometry_Type, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Axisymmetric), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Wall_Boundary_Treatments), 
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
    MPI::COMM_WORLD.Bcast(&(IP.Dual_Time_Stepping), 
                          1, 
			  MPI::INT, 0); 
    MPI::COMM_WORLD.Bcast(&(IP.BluffBody_Data_Usage), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Moving_wall_velocity), 
                          1, 
			  MPI::DOUBLE, 0);

    MPI::COMM_WORLD.Bcast(&(IP.Pressure_Gradient), 
                          1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Reynolds_Number), 
                          1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Kinematic_Viscosity_Wall), 
                          1, 
			  MPI::DOUBLE, 0);
    /**************************************/
    MPI::COMM_WORLD.Bcast(IP.Grid_Type, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.NACA_Aerofoil_Type, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
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
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Idir), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Jdir), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Ghost_Cells), 
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
    MPI::COMM_WORLD.Bcast(&(IP.Length_Shroud), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Radius_Shroud), 
                          1, 
                          MPI::DOUBLE, 0);

    MPI::COMM_WORLD.Bcast(&(IP.Length_BluffBody), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Radius_BluffBody), 
                          1, 
                          MPI::DOUBLE, 0);

    MPI::COMM_WORLD.Bcast(&(IP.Radius_Orifice), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.BluffBody_Coflow_Air_Velocity), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.BluffBody_Coflow_Fuel_Velocity), 
                          1, 
                          MPI::DOUBLE, 0);
    
    MPI::COMM_WORLD.Bcast(&(IP.Radius_Inlet_Pipe), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Radius_Combustor_Tube), 
                          1, 
                          MPI::DOUBLE, 0);
    
    MPI::COMM_WORLD.Bcast(&(IP.Length_Inlet_Pipe), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Length_Combustor_Tube), 
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
			  INPUT_PARAMETER_LENGTH_LESPREMIXED2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(&(IP.BCs_Specified),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(IP.BC_North_Type,
			  INPUT_PARAMETER_LENGTH_LESPREMIXED2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(IP.BC_South_Type,
			  INPUT_PARAMETER_LENGTH_LESPREMIXED2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(IP.BC_East_Type,
			  INPUT_PARAMETER_LENGTH_LESPREMIXED2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(IP.BC_West_Type,
			  INPUT_PARAMETER_LENGTH_LESPREMIXED2D,
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
    if (!CFFC_Primary_MPI_Processor()) {
      IP.ICEMCFD_FileNames = new char*[3];
       for (int i = 0; i < 3; i++) {
          IP.ICEMCFD_FileNames[i] = new char[INPUT_PARAMETER_LENGTH_LESPREMIXED2D];
       } /* endfor */
    } /* endif */
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[0], 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[1], 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[2], 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                          MPI::CHAR, 0);

    // AMR & Refinement Parameters
    MPI::COMM_WORLD.Bcast(&(IP.AMR), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.AMR_Frequency),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Uniform_Mesh_Coarsenings), 
                          1, 
                          MPI::INT, 0);
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
    // Morton Ordering Parameters
    MPI::COMM_WORLD.Bcast(&(IP.Morton), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Morton_Reordering_Frequency),
                          1,
                          MPI::INT,0);
    // File Names
    MPI::COMM_WORLD.Bcast(IP.Output_File_Name, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Grid_File_Name, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Grid_Definition_File_Name, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Restart_File_Name, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Gnuplot_File_Name, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Output_Format_Type, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
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
    
    // NKS Parameters
    IP.NKS_IP.Broadcast_Input_Parameters();

    if (!CFFC_Primary_MPI_Processor()) {
       IP.Number_of_Processors = CFFC_MPI::Number_of_Processors;
    } /* endif */
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Per_Processor), 
                          1, 
                          MPI::INT, 0); 
    
    // Freeze_Limiter
    MPI::COMM_WORLD.Bcast(&(IP.Freeze_Limiter), 
			  1, 
			  MPI::INT, 0);

    MPI::COMM_WORLD.Bcast(&(IP.i_Residual_Variable), 
			  1, 
			  MPI::INT, 0);
    // Freeze_Limiter_Residual_Level
    MPI::COMM_WORLD.Bcast(&(IP.Freeze_Limiter_Residual_Level),
                          1, 
                          MPI::DOUBLE, 0);

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
void Broadcast_Input_Parameters(LESPremixed2D_Input_Parameters &IP,
                                MPI::Intracomm &Communicator, 
                                const int Source_CPU) {
 
    int Source_Rank = 0;
    int i;
   
    Communicator.Bcast(IP.Input_File_Name, 
                       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.Line_Number), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Time_Integration_Type, 
                       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
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
    Communicator.Bcast(&(IP.dTime), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Physical_CFL_Number), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Max_Inner_Steps), 
                       1, 
                       MPI::INT, Source_Rank);
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
                       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Reconstruction), 
                       1, 
                       MPI::INT, Source_Rank);
    // Flux functions:
    Communicator.Bcast(IP.Viscous_Flux_Evaluation_Type, 
                       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Viscous_Flux_Evaluation), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Limiter_Type, 
                       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Limiter), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Flux_Function_Type, 
                       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Flux_Function), 
                       1, 
                       MPI::INT, Source_Rank);
    // Initial conditions:
    Communicator.Bcast(IP.ICs_Type, 
                       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_ICs), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.i_Grid_Level), 
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
    Communicator.Bcast(&(IP.Mach_Number_Reference_target), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Flow_Angle), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Re_lid), 
                       1, 
		       MPI::DOUBLE,Source_Rank);
    /*********************************************************
     ******************* LESPREMIXED2D SPECIFIC **************
     *********************************************************/
    Communicator.Bcast(IP.CFFC_Path, 
		       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
		       MPI::CHAR, Source_Rank);
    //fresh and burnt fuel mass fraction
    Communicator.Bcast(&(IP.Fresh_Fuel_Mass_Fraction), 
                          1, 
                          MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Burnt_Fuel_Mass_Fraction), 
                          1, 
                          MPI::DOUBLE, Source_Rank);
    //fresh density
    Communicator.Bcast(&(IP.Fresh_Density), 
                          1, 
                          MPI::DOUBLE, Source_Rank);
    //Initial turbulence kinetic energy
    Communicator.Bcast(&(IP.TKEo), 
                          1, 
                          MPI::DOUBLE, Source_Rank);
    //Spectrum flag
    Communicator.Bcast(IP.Spectrum_Type, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                          MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Spectrum), 
                          1, 
                          MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Rescale_Spectrum), 
                          1, 
                          MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Read_Fluctuations_From_File), 
                          1, 
                          MPI::INT, Source_Rank);
    //reaction name
    Communicator.Bcast(IP.React_Name, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
			  MPI::CHAR, Source_Rank);
    // cantera parameters
    Communicator.Bcast(IP.ct_Mech_Name,
		       INPUT_PARAMETER_LENGTH_CHEM2D,
		       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.ct_Mech_File,
		       INPUT_PARAMETER_LENGTH_CHEM2D,
		       MPI::CHAR, Source_Rank);
    //scalar system name
    Communicator.Bcast(IP.Scalar_system_name, 
                          INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
			  MPI::CHAR, Source_Rank);

    //delete orginal dynamic memory before changing num_species
    if(!CFFC_Primary_MPI_Processor()) {  
      IP.Deallocate();  
    } 
    //number of species
    Communicator.Bcast(&(IP.num_species), 
                          1, 
                          MPI::INT, Source_Rank);
    //number of scalars
    Communicator.Bcast(&(IP.num_scalars), 
                          1, 
                          MPI::INT, Source_Rank);
    //transport data
    Communicator.Bcast(IP.trans_type, 
			  INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
			  MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_trans_type), 
			  1, 
			  MPI::INT, Source_Rank);

    //set up dynamic memory
    if(!CFFC_Primary_MPI_Processor()) {  
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
			 INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
			 MPI::CHAR, Source_Rank);
    }
    //scalar names & values
    if(IP.num_scalars > 0){
      for(int i=0; i < IP.num_scalars; i++){
        Communicator.Bcast(&(IP.scalar_values[i]), 
			    1, 
			    MPI::DOUBLE, Source_Rank);
        Communicator.Bcast(IP.Scalars[i], 
			    INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
			    MPI::CHAR, Source_Rank);
      }
    }
    //set reaction and species parameters
    if (!CFFC_Primary_MPI_Processor()) {      
      IP.react_name = IP.React_Name;
      for (int i = 0; i < IP.num_species; i++) {
	IP.multispecies[i] = IP.Multispecies[i];  
      }
      IP.scalar_system_name = IP.Scalar_system_name;
      if(IP.num_scalars > 0){
        for(int i=0; i<IP.num_scalars; i++){
          IP.scalars[i] = IP.Scalars[i];
        }
      }
     
      //load reaction and scalar names
      if (IP.Wo.React.reactset_flag != CANTERA)
	IP.Wo.React.set_reactions(IP.react_name);
      else
	 IP.Wo.React.ct_load_mechanism(IP.ct_mech_file, IP.ct_mech_name);
      IP.Wo.Scal_sys.scalar_set(IP.scalar_system_name);
      
      //Set species if non-reacting
      if( IP.Wo.React.reactset_flag == NO_REACTIONS){
	IP.Wo.React.set_species(IP.multispecies,IP.num_species);
      }  

      //set the data for each
      IP.get_cffc_path();
      IP.Wo.set_species_data(IP.num_species, IP.num_scalars, IP.multispecies,
			     IP.CFFC_Path,
			     IP.Mach_Number_Reference,
			     IP.Schmidt,
			     IP.i_trans_type); 
      IP.Uo.set_species_data(IP.num_species, IP.num_scalars, IP.multispecies,
			     IP.CFFC_Path,
			     IP.Mach_Number_Reference,
			     IP.Schmidt,
			     IP.i_trans_type);
      IP.Wo.set_initial_values(IP.mass_fractions);
      IP.Uo.set_initial_values(IP.mass_fractions);
      if(IP.num_scalars > 0){
        IP.Wo.set_initial_values_scal(IP.scalar_values);   
        IP.Uo.set_initial_values_scal(IP.scalar_values);
      }

      //set proper Temp & Pressure instead of defaults
      IP.Wo.rho = IP.Pressure/(IP.Wo.Rtot()*IP.Temperature); 
      IP.Wo.p = IP.Pressure;	
      IP.Wo.v.zero();

      IP.Uo = U(IP.Wo);
    }

    // SFS modelling parameters
    Communicator.Bcast(&(IP.Smagorinsky_Constant), 
                          1, 
                          MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Yoshizawa_coefficient), 
                          1, 
                          MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.filter_width), 
                          1, 
                          MPI::DOUBLE, Source_Rank);
    // premixed flame parameters
    Communicator.Bcast(&(IP.laminar_flame_thickness), 
                          1, 
                          MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.laminar_flame_speed), 
                          1, 
                          MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.TFactor), 
                          1, 
                          MPI::DOUBLE, Source_Rank);  
    Communicator.Bcast(&(IP.adiabatic_temp), 
                          1, 
                          MPI::DOUBLE, Source_Rank);  
    Communicator.Bcast(&(IP.equivalence_ratio), 
                          1, 
                          MPI::DOUBLE, Source_Rank);  
    Communicator.Bcast(&(IP.reactants_den), 
                          1, 
                          MPI::DOUBLE, Source_Rank);  

    // Reset the static variables.
    IP.Wo.set_SFSmodel_variables(IP.filter_width,
				 IP.Smagorinsky_Constant,
				 IP.Yoshizawa_coefficient);
    IP.Uo.set_SFSmodel_variables(IP.filter_width,
				 IP.Smagorinsky_Constant,
				 IP.Yoshizawa_coefficient);
    IP.Wo.set_premixed_flame_variables(IP.laminar_flame_thickness,
				       IP.laminar_flame_speed,
				       IP.TFactor,
                                       IP.adiabatic_temp,
                                       IP.equivalence_ratio,
                                       IP.reactants_den);
    IP.Uo.set_premixed_flame_variables(IP.laminar_flame_thickness,
				       IP.laminar_flame_speed,
				       IP.TFactor,
                                       IP.adiabatic_temp,
                                       IP.equivalence_ratio,
                                       IP.reactants_den);

    /*********************************************************
     ******************* LESPREMIXED2D END *******************
     *********************************************************/
   
    if(!CFFC_Primary_MPI_Processor()) {   
      IP.Wo.v.x = IP.Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Flow_Angle/360.00);
      IP.Wo.v.y = IP.Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Flow_Angle/360.00);
    }

    /********************************************/
    Communicator.Bcast(IP.Flow_Type, 
                       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.FlowType), 
                       1, 
                       MPI::INT, Source_Rank);

    /********************************************/
    Communicator.Bcast(IP.Flow_Geometry_Type, 
                       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.Axisymmetric), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Wall_Boundary_Treatments), 
                       1, 
                       MPI::INT, Source_Rank);
//     Communicator.Bcast(&(IP.Stretch_Level), 
//                        1, 
//                        MPI::INT, Source_Rank);
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
    Communicator.Bcast(&(IP.Dual_Time_Stepping), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.BluffBody_Data_Usage), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Moving_wall_velocity), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Pressure_Gradient), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Reynolds_Number), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Kinematic_Viscosity_Wall), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    /********************************************/
    Communicator.Bcast(IP.Grid_Type, 
                       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.NACA_Aerofoil_Type, 
                       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
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
    Communicator.Bcast(&(IP.Length_Shroud), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Radius_Shroud), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Length_BluffBody), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Radius_BluffBody), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Radius_Orifice), 
                       1, 
                       MPI::DOUBLE, Source_Rank);

    Communicator.Bcast(&(IP.BluffBody_Coflow_Air_Velocity), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    
    Communicator.Bcast(&(IP.BluffBody_Coflow_Fuel_Velocity), 
                       1, 
                       MPI::DOUBLE, Source_Rank);

    Communicator.Bcast(&(IP.Radius_Inlet_Pipe), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Radius_Combustor_Tube), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Length_Inlet_Pipe), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Length_Combustor_Tube), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
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
		       INPUT_PARAMETER_LENGTH_LESPREMIXED2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(&(IP.BCs_Specified),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(IP.BC_North_Type,
		       INPUT_PARAMETER_LENGTH_LESPREMIXED2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(IP.BC_South_Type,
		       INPUT_PARAMETER_LENGTH_LESPREMIXED2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(IP.BC_East_Type,
		       INPUT_PARAMETER_LENGTH_LESPREMIXED2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(IP.BC_West_Type,
		       INPUT_PARAMETER_LENGTH_LESPREMIXED2D,
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
    if (!(CFFC_MPI::This_Processor_Number == Source_Rank)) {
       IP.ICEMCFD_FileNames = new char*[3];
       for (i = 0; i < 3; i++) {
          IP.ICEMCFD_FileNames[i] = new char[INPUT_PARAMETER_LENGTH_LESPREMIXED2D];
       } /* endfor */
    } /* endif */
    Communicator.Bcast(IP.ICEMCFD_FileNames[0], 
                       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.ICEMCFD_FileNames[1], 
                       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.ICEMCFD_FileNames[2], 
                       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                       MPI::CHAR, Source_Rank);

    // AMR & Refinement Parameters
    Communicator.Bcast(&(IP.AMR), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.AMR_Frequency),
                       1,
                       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Uniform_Mesh_Coarsenings), 
                       1, 
                       MPI::INT, Source_Rank);
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
    // Morton Ordering Parameters
    Communicator.Bcast(&(IP.Morton), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Morton_Reordering_Frequency),
                       1,
                       MPI::INT,Source_Rank);
    // File Names
    Communicator.Bcast(IP.Output_File_Name, 
                       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Grid_File_Name, 
                       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Grid_Definition_File_Name, 
                       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Restart_File_Name, 
                       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Gnuplot_File_Name, 
                       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Output_Format_Type, 
                       INPUT_PARAMETER_LENGTH_LESPREMIXED2D, 
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
    // Multigrid Related Parameters
    IP.NKS_IP.Broadcast_Input_Parameters(Communicator,
					 Source_CPU);

    if (!(CFFC_MPI::This_Processor_Number == Source_Rank)) {
       IP.Number_of_Processors = CFFC_MPI::Number_of_Processors;
    } /* endif */
    Communicator.Bcast(&(IP.Number_of_Blocks_Per_Processor), 
                       1, 
                       MPI::INT, Source_Rank);
   // Freeze_Limiter
    Communicator.Bcast(&(IP.Freeze_Limiter), 
                       1, 
                       MPI::INT, Source_Rank);

    Communicator.Bcast(&(IP.i_Residual_Variable), 
                       1, 
                       MPI::INT, Source_Rank);

    // Freeze_Limiter_Residual_Level
    Communicator.Bcast(&(IP.Freeze_Limiter_Residual_Level), 
                       1, 
                       MPI::DOUBLE, Source_Rank);

}
#endif

/********************************************************
 * Routine: Get_Next_Input_Control_Parameter            *
 *                                                      *
 * Get the next input control parameter from the input  *
 * file.                                                *
 *                                                      *
 ********************************************************/
void Get_Next_Input_Control_Parameter(LESPremixed2D_Input_Parameters &IP) {

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
int Parse_Next_Input_Control_Parameter(LESPremixed2D_Input_Parameters &IP) {

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
	   // Jai's Dual Time Stepping
//        } else if (strcmp(IP.Time_Integration_Type, "Dual_Time_Stepping") == 0) {
// 	  IP.i_Time_Integration = TIME_STEPPING_DUAL_TIME_STEPPING;
// 	  IP.Multigrid_IP.i_Dual_Time_Stepping = ON; 
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
       } else if (strcmp(IP.Flux_Function_Type, "AUSM_plus_up") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_AUSM_PLUS_UP;
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
      /****************** LESPREMIXEDD2D ******************************/
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
       } else if (strcmp(IP.ICs_Type, "2DPremixedFlame") == 0 ){
	 IP.i_ICs = IC_LESPREMIXED_2DFLAME;
       } else if (strcmp(IP.ICs_Type, "2DWrinkledFlame") == 0 ){
	 IP.i_ICs = IC_2DWRINKLED_FLAME;
       } else if (strcmp(IP.ICs_Type, "HomogeneousTurbulence") == 0 ){
	 IP.i_ICs = IC_HOMOGENEOUS_TURBULENCE;
       } else if (strcmp(IP.ICs_Type, "Mixing_Layer") == 0 ){
	 IP.i_ICs = IC_MIXING_LAYER;
       } else if (strcmp(IP.ICs_Type, "Acoustic_Wave") == 0 ){
	 IP.i_ICs = IC_ACOUSTIC_WAVE;
       } else if (strcmp(IP.ICs_Type, "Impacting_Flow_Xdir") == 0 ){
	 IP.i_ICs = IC_SLOWLY_IMPACTING_XDIR ;
       } else if (strcmp(IP.ICs_Type, "Vortex_Xdir") == 0 ){
	 IP.i_ICs = IC_VORTEX_XDIR;
       } else if (strcmp(IP.ICs_Type, "Pipe_Flow") == 0) {
          IP.i_ICs = IC_TURBULENT_PIPE_FLOW;
       } else if (strcmp(IP.ICs_Type, "Coflow") == 0) {
          IP.i_ICs = IC_TURBULENT_COFLOW;
       }else if (strcmp(IP.ICs_Type, "Driven_Cavity_Flow") == 0) {
          IP.i_ICs = IC_VISCOUS_DRIVEN_CAVITY_FLOW;
       } else if (strcmp(IP.ICs_Type, "Turbulent_Dump_Combustor") == 0) {
          IP.i_ICs = IC_TURBULENT_DUMP_COMBUSTOR;
       }else if (strcmp(IP.ICs_Type, "Turbulent_Diffusion_Flame") == 0) {
          IP.i_ICs = IC_TURBULENT_DIFFUSION_FLAME;
       }else if (strcmp(IP.ICs_Type, "Turbulent_Free_Jet_Flame") == 0) {
          IP.i_ICs = IC_FREE_JET_FLAME;
       }else if (strcmp(IP.ICs_Type, "Restart") == 0) {
          IP.i_ICs = IC_RESTART;
       } else {
	 cerr<<"\n Not a vaild Initial Condition for LESPremixed2D, check LESPremixed2DInput.cc :  "<<IP.ICs_Type; exit(1);
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
       } else if (strcmp(IP.Grid_Type, "Testing") == 0) {
          IP.i_Grid = GRID_TEST;
          IP.Box_Width = ONE;
          IP.Box_Height = ONE;
      } else if (strcmp(IP.Grid_Type, "1DFlame") == 0) {
          IP.i_Grid = GRID_1DFLAME;
          IP.Box_Width = ONE;
          IP.Box_Height = 0.1;
      } else if (strcmp(IP.Grid_Type, "Laminar_Flame") == 0) {
          IP.i_Grid = GRID_LAMINAR_FLAME;
          IP.Pipe_Length = 0.1;
          IP.Pipe_Radius = 0.2;
      } else if (strcmp(IP.Grid_Type, "Periodic_Box") == 0) {
          IP.i_Grid = GRID_PERIODIC_BOX;
          IP.Box_Width = 6.283185307;
          IP.Box_Height = 6.283185307;
      } else if (strcmp(IP.Grid_Type, "Mixing_Layer_Box") == 0) {
          IP.i_Grid = GRID_MIXING_LAYER_BOX;
          IP.Box_Width = TWO;
          IP.Box_Height = TWO;
      } else if (strcmp(IP.Grid_Type, "Vortex_Box") == 0) {
          IP.i_Grid = GRID_VORTEX_BOX;
          IP.Box_Width = TWO;
          IP.Box_Height = TWO;
      } else if (strcmp(IP.Grid_Type, "2DTurbulent_Flame") == 0) {
          IP.i_Grid = GRID_2DTURBULENT_PREMIXED_FLAME;
          IP.Box_Width = TWO;
          IP.Box_Height = TWO;
       } else if (strcmp(IP.Grid_Type, "Flat_Plate") == 0) {
          IP.i_Grid = GRID_FLAT_PLATE;
          IP.Plate_Length = ONE;
       } else if (strcmp(IP.Grid_Type, "Pipe") == 0) {
          IP.i_Grid = GRID_PIPE;
	  //channel (laufer 1953)
	  //IP.Pipe_Length = 1.524 ;
          //IP.Pipe_Radius = 0.0635;
	  //pipe (laufer 1953)
          IP.Kinematic_Viscosity_Wall = 1.4590E-5;
	  IP.Pipe_Length = 4.877;
	  IP.Pipe_Radius = 0.123;
       }else if (strcmp(IP.Grid_Type, "Driven_Cavity_Flow") == 0) {
	 IP.i_Grid = GRID_DRIVEN_CAVITY_FLOW;
	 IP.Box_Width = 0.00044;
	 IP.Box_Height = 0.00044;	  
       }else if (strcmp(IP.Grid_Type, "Bluff_Body") == 0) {
          IP.i_Grid = GRID_BLUFF_BODY;
          IP.Radius_BluffBody = 0.025;
          IP.Radius_Shroud = 0.07;
          IP.Radius_Orifice = 0.0018;
          IP.Length_BluffBody = 0.1;
          IP.Length_Shroud = 0.3;
          
       }else if (strcmp(IP.Grid_Type, "Free_Jet_Flame") == 0) {
          IP.i_Grid = GRID_FREE_JET_FLAME;
          IP.Radius_Inlet_Pipe = 0.00387;
          IP.Radius_Shroud = 0.48;
          IP.Length_Shroud = 1.6;
       //    IP.Radius_Inlet_Pipe = 0.00305;
//           IP.Radius_Shroud = 0.96;
//           IP.Length_Shroud = 3.2;
          
       }else if (strcmp(IP.Grid_Type, "Dump_Combustor") == 0) {
          IP.i_Grid = GRID_DUMP_COMBUSTOR;
          IP.Radius_Inlet_Pipe = 0.0508;
          IP.Radius_Combustor_Tube = 0.0762;
          IP.Length_Inlet_Pipe = 0.127;
          IP.Length_Combustor_Tube = 0.508;
          
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

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Ghost_Cells") == 0) {
       i_command = 100;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Ghost_Cells;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Ghost_Cells < 2) i_command = INVALID_INPUT_VALUE;

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
	   IP.Local_Time_Stepping != SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER &&
           IP.Local_Time_Stepping != DUAL_TIME_STEPPING &&
           IP.Local_Time_Stepping != DUAL_LOW_MACH_NUMBER_PRECONDITIONER &&
           IP.Local_Time_Stepping != DUAL_SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER) {
         IP.Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
	}
       if (IP.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER ||
       	 IP.Local_Time_Stepping == SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER ||
         IP.Local_Time_Stepping == DUAL_LOW_MACH_NUMBER_PRECONDITIONER ||
         IP.Local_Time_Stepping == DUAL_SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER) {
	 IP.Preconditioning = 1;
       }
       if (IP.Local_Time_Stepping == DUAL_TIME_STEPPING ||
           IP.Local_Time_Stepping == DUAL_LOW_MACH_NUMBER_PRECONDITIONER ||
	   IP.Local_Time_Stepping == DUAL_SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER) {
         IP.Dual_Time_Stepping = 1;
       }       

    } else if (strcmp(IP.Next_Control_Parameter, "Maximum_Number_of_Time_Steps") == 0) {
       i_command = 16;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Maximum_Number_of_Time_Steps;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Maximum_Number_of_Time_Steps < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Maximum_Number_of_Inner_Time_Steps") == 0) {
       i_command = 1600;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Max_Inner_Steps;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Max_Inner_Steps < 1) i_command = INVALID_INPUT_VALUE;

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
       if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	 IP.Smooth_Bump = ON;
       } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
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

    } else if (strcmp(IP.Next_Control_Parameter, "Length_BluffBody") == 0) {
       i_command =38 ;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Length_BluffBody;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Length_BluffBody <ZERO) i_command = INVALID_INPUT_VALUE;
    } else if (strcmp(IP.Next_Control_Parameter, "Radius_Orifice") == 0) {
       i_command =39 ;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Radius_Orifice;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Radius_Orifice <ZERO) i_command = INVALID_INPUT_VALUE;   
    }else if (strcmp(IP.Next_Control_Parameter, "BluffBody_Coflow_Air_Velocity") == 0) {
       i_command =40 ;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.BluffBody_Coflow_Air_Velocity;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.BluffBody_Coflow_Air_Velocity <ZERO) i_command = INVALID_INPUT_VALUE;
       
    }else if (strcmp(IP.Next_Control_Parameter, "BluffBody_Coflow_Fuel_Velocity") == 0) {
       i_command =41 ;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.BluffBody_Coflow_Fuel_Velocity;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.BluffBody_Coflow_Fuel_Velocity <ZERO) i_command = INVALID_INPUT_VALUE;       


       /***********************************************************************
        ************************ LESPREMIXED2D SPECIFIC ***********************
        For LESPremixed2D the following input parameters are added:

             "Scalar_system"

             "Reaction_Mechanism"

             "User_Reaction_Mechanism"
 
             "Species"

       *************************************************************************
       *************************************************************************/
       
       /*************************************/
       /******   SYSTEM OF SCALARS     ******/
       /*************************************/
    } else if (strcmp(IP.Next_Control_Parameter, "Scalar_system") == 0) {
       i_command = 200;
       Get_Next_Input_Control_Parameter(IP);
       IP.Deallocate_scalars();  //DEALLOCATE BEFORE CHANGING num_scalars
     
       //convert IP to string & define setup for Scalar System
       IP.scalar_system_name = IP.Next_Control_Parameter;

       IP.Wo.Scal_sys.scalar_set(IP.scalar_system_name);  
       IP.num_scalars = IP.Wo.Scal_sys.num_scalars;
      
       IP.Allocate_scalars();

       //Get scalars
       for(int i=0; i<IP.num_scalars; i++){
	 IP.scalars[i] = IP.Wo.Scal_sys.scalars[i];
       }   
       
       //Get next line and read in scalar values or set defaults
       Get_Next_Input_Control_Parameter(IP);
       if (strcmp(IP.Next_Control_Parameter, "Scalar_values") == 0){
	 //Get Initial scalar values from user	
	 for(int i=0; i<IP.num_scalars; i++){
	   IP.Input_File >> IP.scalar_values[i];
	 }
	
	 //Set inital Values;
         IP.Wo.nscal = IP.num_scalars;
         IP.Uo.nscal = IP.num_scalars;  
	 IP.Wo.set_initial_values_scal(IP.scalar_values);  
	 IP.Uo.set_initial_values_scal(IP.scalar_values);  
	 IP.Uo = U(IP.Wo);        
	 
	 //fudge the line number and istream counters
	 IP.Input_File.getline(buffer, sizeof(buffer));  
	 IP.Line_Number = IP.Line_Number + 1; 
       
	 //If no scalar values is set to defaults (ZERO) 
       } else {
         IP.Wo.nscal = IP.num_scalars;
         IP.Uo.nscal = IP.num_scalars;  
         IP.Wo.set_initial_values_scal();  
	 IP.Uo.set_initial_values_scal();
	 IP.Uo = U(IP.Wo);
         for(int i=0; i<IP.num_scalars; i++){
	   IP.scalar_values[i] = IP.Wo.scalar[i];
	 }
	 IP.Line_Number = IP.Line_Number - 1 ;
       }

       /*************************************/
       /**** REACTIONS SET FOR HARDCODED ****/
       /*************************************/
    } else if (strcmp(IP.Next_Control_Parameter, "Reaction_Mechanism") == 0) {
       i_command = 201;
       Get_Next_Input_Control_Parameter(IP);
       IP.Deallocate_species();  //DEALLOCATE BEFORE CHANGING num_species
       int flag =0;

       //convert IP to string & define setup which Reaction Mechanism
       IP.react_name = IP.Next_Control_Parameter;
       IP.Wo.React.set_reactions(IP.react_name);
   
       IP.num_species = IP.Wo.React.num_species;      
       IP.Allocate_species();

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
       IP.Wo.set_species_data(IP.num_species, IP.num_scalars, IP.multispecies,
			      IP.CFFC_Path,
			      IP.Mach_Number_Reference,
			      IP.Schmidt,
			      IP.i_trans_type); 
       IP.Uo.set_species_data(IP.num_species, IP.num_scalars, IP.multispecies,
			      IP.CFFC_Path,
			      IP.Mach_Number_Reference,
			      IP.Schmidt,
			      IP.i_trans_type);
  
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

	 // assign initial fuel mass fraction
	 IP.Fresh_Fuel_Mass_Fraction = IP.mass_fractions[0];

	 //Set inital Values; 
	 IP.Wo.set_initial_values(IP.mass_fractions);  
	 IP.Uo.set_initial_values(IP.mass_fractions);
	 if (IP.Wo.nscal > 0 && IP.scalar_values != NULL) {
           IP.Wo.set_initial_values_scal(IP.scalar_values);  
	   IP.Uo.set_initial_values_scal(IP.scalar_values); 
	 }
	 IP.Uo = U(IP.Wo);
	 
	 //fudge the line number and istream counters
	 IP.Input_File.getline(buffer, sizeof(buffer));  
	 IP.Line_Number = IP.Line_Number + 1; 
        
       //Spit out appropriate mass fractions and exit
       // } else if (strcmp(IP.Next_Control_Parameter, "Equivalence_Ratio") == 0){        
// 	 double phi;
// 	 IP.Input_File >> phi;
// 	 Equivalence_Ratio(phi);

	 //If no mass fraction data is set to defaults (all equal to 1/num_species) 
       } else {
	 if (IP.Wo.nscal > 0  && IP.scalar_values != NULL) {
           IP.Wo.set_initial_values_scal(IP.scalar_values);  
	   IP.Uo.set_initial_values_scal(IP.scalar_values); 
	 }
	 IP.Uo = U(IP.Wo);
	 IP.Line_Number = IP.Line_Number - 1 ;
       }
         
       
       /***************************************/
       /**** REACTIONS FOR USER DEFINED *******/ 
       /***************************************/
    } else if (strcmp(IP.Next_Control_Parameter, "User_Reaction_Mechanism") == 0) { 
      // this will be added but its not quite yet
      i_command = 202;
      cout<<endl<<IP.Next_Control_Parameter<<"\n not currently available in freeware version :)\n";
      i_command = INVALID_INPUT_VALUE;   
      
       /*************************************/
       /***** REACTIONS SET FOR CANTERA *****/
       /*************************************/
    } else if (strcmp(IP.Next_Control_Parameter, "Cantera_Reaction_Mechanism") == 0) {
       i_command = 203;

       Get_Next_Input_Control_Parameter(IP);
       IP.Deallocate();  //DEALLOCATE BEFORE CHANGING num_species
       int flag =0;

       //convert IP to string & define Reaction Mechanism
       IP.react_name = "CANTERA";
       IP.ct_mech_name = IP.Next_Control_Parameter;

       //get the mechanism file name and load the mechanism
       Get_Next_Input_Control_Parameter(IP);
       if (strcmp(IP.Next_Control_Parameter, "Mechanism_File") == 0){
	 Get_Next_Input_Control_Parameter(IP);
	 IP.ct_mech_file = IP.Next_Control_Parameter;
	 IP.Wo.React.ct_load_mechanism(IP.ct_mech_file, IP.ct_mech_name);   
	 IP.num_species = IP.Wo.React.num_species;      
       // if no reaction file was given, return in error
       } else {
	 IP.num_species = 0;
	 i_command = INVALID_INPUT_CODE;
       }

       // allocate storage
       IP.Allocate();

       //Get species and load appropriate data
       for(int i=0; i<IP.num_species; i++){
	 IP.multispecies[i] = IP.Wo.React.species[i];
	 IP.Schmidt[i] = IP.Global_Schmidt;
       }   

       //
       //Get next line and read in Schmidt numbers else will use defaults
       //
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

       //
       //Set appropriate species data
       //
       IP.Wo.set_species_data(IP.num_species, IP.num_scalars, IP.multispecies,
	 		      IP.CFFC_Path,
			      IP.Mach_Number_Reference,
			      IP.Schmidt,
			      IP.i_trans_type); 
       IP.Uo.set_species_data(IP.num_species, IP.num_scalars, IP.multispecies,
			      IP.CFFC_Path,
			      IP.Mach_Number_Reference,
			      IP.Schmidt,
			      IP.i_trans_type);
  
       //
       //Get next line and read in mass fractions or set defaults
       //
       if(flag){
	 Get_Next_Input_Control_Parameter(IP);
       }
       if (strcmp(IP.Next_Control_Parameter, "Mass_Fractions") == 0){

	 // Get Initial Mass Fractions from user 
	 // Here we use Cantera to parse the string of the form:
	 //       CH4:0.5, O2:0.5
	 // All other species will be assumed to have 0 mass fractions.  
	 // Cantera also normalizes the mass fractions to sum to unity.  
	 // Returns them in an array.
	 IP.Input_File.getline(buffer, sizeof(buffer)); 
	 IP.Line_Number = IP.Line_Number + 1 ;
	 IP.Wo.React.ct_parse_mass_string(buffer, IP.mass_fractions);

	 //Set inital Values; 
	 IP.Wo.set_initial_values(IP.mass_fractions);  
	 IP.Uo.set_initial_values(IP.mass_fractions);  
	 IP.Uo = U(IP.Wo);
	 
	 //fudge the line number and istream counters
	 IP.Input_File.getline(buffer, sizeof(buffer));  
	 IP.Line_Number = IP.Line_Number + 1; 
        
       //If no mass fraction data is set to defaults (all equal to 1/num_species) 
       } else {
	 IP.Uo = U(IP.Wo);
	 IP.Line_Number = IP.Line_Number - 1 ;
       }

      /******************************************/
      /**** NON REACTING, BUT MULTIPLE GASES ****/
      /******************************************/
    } else if (strcmp(IP.Next_Control_Parameter, "Species") == 0) { 
      i_command = 204;
  
      IP.Deallocate_species();  //DEALLOCATE BEFORE CHANGING num_species
 
      // Non Reaction case so set NO_REACTIONS flag in reactions class
      IP.react_name ="NO_REACTIONS";
      IP.Wo.React.set_reactions(IP.react_name);
     
      //read in the number of species (should be first in line) 
      IP.Input_File>>IP.num_species;     
      IP.Allocate_species();

      //read in species names
      for(int i=0; i<IP.num_species; i++){
	IP.Input_File >> IP.multispecies[i];
        IP.Schmidt[i] = IP.Global_Schmidt;
      }

      //copy names into Reaction class for storage
      IP.Wo.React.set_species(IP.multispecies,IP.num_species);
       
      //Setup State class data and find species thermo and transport properties
      IP.Wo.set_species_data(IP.num_species, IP.num_scalars, IP.multispecies,
			     IP.CFFC_Path,
			     IP.Mach_Number_Reference,
			     IP.Schmidt,
			     IP.i_trans_type); 
      IP.Uo.set_species_data(IP.num_species, IP.num_scalars, IP.multispecies,
			     IP.CFFC_Path,
			     IP.Mach_Number_Reference,
			     IP.Schmidt,
			     IP.i_trans_type);
    
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
	 if (IP.Wo.nscal > 0 && IP.scalar_values != NULL) {
           IP.Wo.set_initial_values_scal(IP.scalar_values);  
	   IP.Uo.set_initial_values_scal(IP.scalar_values); 
	 }
	 IP.Uo = U(IP.Wo);
	 
	 //fudge the line number and istream counters
	 IP.Input_File.getline(buffer, sizeof(buffer));  
	 IP.Line_Number = IP.Line_Number + 1; 
       } 
       //If no mass fraction data is set to defaults (all equal to 1/num_species)
       else{
	 if (IP.Wo.nscal > 0 && IP.scalar_values != NULL) {
           IP.Wo.set_initial_values_scal(IP.scalar_values);  
	   IP.Uo.set_initial_values_scal(IP.scalar_values); 
	 }
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

       /************* HEAT SOURCE ****************/
    } else if (strcmp(IP.Next_Control_Parameter, "Temperature_Rise") == 0) {
       i_command = 43;
       IP.Line_Number = IP.Line_Number + 1;
       double temp;
       IP.Input_File >> temp;
       IP.Heat_Source = IP.Wo.rho*(IP.Wo.h(temp+IP.Temperature) - IP.Wo.Rtot()*IP.Temperature);

       IP.Input_File.getline(buffer, sizeof(buffer));
         
       /***********************************************************************
	**************** END LESPREMIXED2D MODIFICATIONS **********************
	***********************************************************************/

       /********** TRANSPORT DATA ***********/
    } else if (strcmp(IP.Next_Control_Parameter, "Transport_Data_Type") == 0) {
      i_command = 390;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.trans_type, 
	     IP.Next_Control_Parameter);
      if (strcmp(IP.trans_type, "Transport-NASA") == 0) {
	IP.i_trans_type = TRANSPORT_NASA;
      } else if (strcmp(IP.trans_type, "Transport-Lennard-Jones") == 0) {
	IP.i_trans_type = TRANSPORT_LENNARD_JONES;
      } // end if


      /********** SFS PARAMETERS ***********/
    } else if (strcmp(IP.Next_Control_Parameter, "Smagorinsky_Coefficient") == 0) {
      i_command = 391;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Smagorinsky_Constant;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Smagorinsky_Constant < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Yoshizawa_Coefficient") == 0) {
      i_command = 392;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Yoshizawa_coefficient;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Yoshizawa_coefficient < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Filter_Width") == 0) {
      i_command = 393;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.filter_width;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.filter_width < ZERO) i_command = INVALID_INPUT_VALUE;

       /********** PREMIXED FLAME PARAMETERS  ***********/
    } else if (strcmp(IP.Next_Control_Parameter, "Laminar_Flame_Thickness") == 0) {
      i_command = 394;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.laminar_flame_thickness;
      IP.Input_File.getline(buffer, sizeof(buffer));   	 	 
      if (IP.laminar_flame_thickness <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Laminar_Flame_Speed") == 0) {
      i_command = 395;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.laminar_flame_speed;
      IP.Input_File.getline(buffer, sizeof(buffer));   	 	 
      if (IP.laminar_flame_speed <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Thickening_Factor") == 0) {
      i_command = 396;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.TFactor;
      IP.Input_File.getline(buffer, sizeof(buffer));   	 	 
      if (IP.TFactor < 1.0) i_command = INVALID_INPUT_VALUE;
    } else if (strcmp(IP.Next_Control_Parameter, "Adiabatic_Temperature") == 0) {
      i_command = 400;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.adiabatic_temp;
      IP.Input_File.getline(buffer, sizeof(buffer));   	 	 
      if (IP.adiabatic_temp < 0.0) i_command = INVALID_INPUT_VALUE;
    } else if (strcmp(IP.Next_Control_Parameter, "Equivalence_Ratio") == 0) {
      i_command = 401;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.equivalence_ratio;
      IP.Input_File.getline(buffer, sizeof(buffer));   	 	 
      if (IP.equivalence_ratio < 0.0) i_command = INVALID_INPUT_VALUE;
    } else if (strcmp(IP.Next_Control_Parameter, "Reactants_Density") == 0) {
      i_command = 402;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.reactants_den;
      IP.Input_File.getline(buffer, sizeof(buffer));   	 	 
      if (IP.reactants_den < 1.0) i_command = INVALID_INPUT_VALUE;

      /********** ENERGY SPECTRUM ***********/
    } else if (strcmp(IP.Next_Control_Parameter, "Energy_Spectrum") == 0) {
      i_command = 397;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Spectrum_Type, 
             IP.Next_Control_Parameter);
       if (strcmp(IP.Spectrum_Type, "Lee-Reynolds") == 0) {
          IP.i_Spectrum = LEE_REYNOLDS;
       } else if (strcmp(IP.Spectrum_Type, "Laval-Nazarenko") == 0) {
          IP.i_Spectrum = LAVAL_NAZARENKO;
       } else if (strcmp(IP.Spectrum_Type, "von-Karman-Pao") == 0) {
          IP.i_Spectrum = VON_KARMAN_PAO;
       } else if (strcmp(IP.Spectrum_Type, "Haworth-Poinsot") == 0) {
          IP.i_Spectrum = HAWORTH_POINSOT;
       } else if (strcmp(IP.Spectrum_Type, "Chasnov") == 0) {
          IP.i_Spectrum = CHASNOV;
       } else if (strcmp(IP.Spectrum_Type, "Bell-Day") == 0) {
          IP.i_Spectrum = BELL_DAY;
       } else {
          IP.i_Spectrum = VON_KARMAN_PAO; 
       } /* end if */
    } else if (strcmp(IP.Next_Control_Parameter, "Rescale_Spectrum") == 0) {
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Rescale_Spectrum;
      IP.Input_File.getline(buffer, sizeof(buffer));   	 	 
      if (IP.Rescale_Spectrum < 0  ||  IP.Rescale_Spectrum > 1 ) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Read_Fluctuations_From_File") == 0) {
      i_command = 398;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Read_Fluctuations_From_File;
      IP.Input_File.getline(buffer, sizeof(buffer));   	 	 
      if (IP.Read_Fluctuations_From_File < 0  ||  IP.Read_Fluctuations_From_File > 1 ) {
	i_command = INVALID_INPUT_VALUE;
      }


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
      IP.Input_File >> IP.Mach_Number_Reference; // >> IP.Mach_Number_Reference_target;
      //if( IP.Mach_Number_Reference_target > IP.Mach_Number_Reference) 
      IP.Mach_Number_Reference_target = IP.Mach_Number_Reference; 
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
       
    } else if (strcmp(IP.Next_Control_Parameter,"Re_lid") == 0) {
      i_command = 42;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Re_lid;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Re_lid < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Time_Max") == 0) {
       i_command = 43;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Time_Max;
       IP.Input_File.getline(buffer, sizeof(buffer));
       IP.Time_Max = IP.Time_Max/THOUSAND;
       if (IP.Time_Max < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Physical_Time_Step") == 0) {
       i_command = 4300;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.dTime;
       IP.Input_File.getline(buffer, sizeof(buffer));
       IP.dTime = IP.dTime/THOUSAND;
       if (IP.dTime < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Physical_CFL_Number") == 0) {
       i_command = 4301;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Physical_CFL_Number;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Physical_CFL_Number <= ZERO) i_command = INVALID_INPUT_VALUE;

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
       } else if (strcmp(IP.Flow_Type, "Turbulent-LES") == 0) {
	 IP.FlowType = FLOWTYPE_TURBULENT_LES;
       } else if (strcmp(IP.Flow_Type, "Turbulent-LES-No-Model") == 0) {
	 IP.FlowType = FLOWTYPE_TURBULENT_LES_NO_MODEL;
       } else if (strcmp(IP.Flow_Type, "Turbulent-LES-TF-Smagorinsky") == 0) {
         IP.FlowType = FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY;
       } else if (strcmp(IP.Flow_Type, "Turbulent-LES-TF-k") == 0) {
         IP.FlowType = FLOWTYPE_TURBULENT_LES_TF_K;
       } else if (strcmp(IP.Flow_Type, "Laminar-C") == 0) {
         IP.FlowType = FLOWTYPE_LAMINAR_C;
       } else if (strcmp(IP.Flow_Type, "Laminar-C-Algebraic") == 0) {
         IP.FlowType = FLOWTYPE_LAMINAR_C_ALGEBRAIC;
       } else if (strcmp(IP.Flow_Type, "Laminar-C-FSD") == 0) {
         IP.FlowType = FLOWTYPE_LAMINAR_C_FSD;
       } else if (strcmp(IP.Flow_Type, "Laminar-NGT-C-FSD") == 0) {
         IP.FlowType = FLOWTYPE_LAMINAR_NGT_C_FSD;
       } else if (strcmp(IP.Flow_Type, "Turbulent-LES-C") == 0) {
         IP.FlowType = FLOWTYPE_TURBULENT_LES_C;
       } else if (strcmp(IP.Flow_Type, "Turbulent-LES-C-Algebraic") == 0) {
         IP.FlowType = FLOWTYPE_TURBULENT_LES_C_ALGEBRAIC;
       } else if (strcmp(IP.Flow_Type, "Turbulent-LES-C-FSD-Smagorinsky") == 0) {
         IP.FlowType = FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY;
       } else if (strcmp(IP.Flow_Type, "Turbulent-LES-C-FSD-Charlette") == 0) {
         IP.FlowType = FLOWTYPE_TURBULENT_LES_C_FSD_CHARLETTE;
       } else if (strcmp(IP.Flow_Type, "Turbulent-LES-NGT-C-FSD-Smagorinsky") == 0) {
         IP.FlowType = FLOWTYPE_TURBULENT_LES_NGT_C_FSD_SMAGORINSKY;
       } else if (strcmp(IP.Flow_Type, "Turbulent-LES-C-FSD-K") == 0) {
         IP.FlowType = FLOWTYPE_TURBULENT_LES_C_FSD_K;
       } else if (strcmp(IP.Flow_Type, "Frozen-Turbulent-LES-C-FSD") == 0) {
         IP.FlowType = FLOWTYPE_FROZEN_TURBULENT_LES_C_FSD;
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
          IP.Axisymmetric = PLANAR;
       } else if (strcmp(IP.Flow_Geometry_Type, "Axisymmetric") == 0) {
	 IP.Axisymmetric = AXISYMMETRIC_X;
       } else if (strcmp(IP.Flow_Geometry_Type, "Axisymmetric-x") == 0) {
	 IP.Axisymmetric = AXISYMMETRIC_X;
       } else if (strcmp(IP.Flow_Geometry_Type, "Axisymmetric-y") == 0) {
	 IP.Axisymmetric = AXISYMMETRIC_Y;
       } else {
          IP.Axisymmetric = PLANAR;
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

       /******* LESPREMIXED2D *********/
    } else if (strcmp(IP.Next_Control_Parameter, "Viscous_Flux_Evaluation_Type") == 0) {
      i_command = 510; 
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Viscous_Flux_Evaluation_Type,IP.Next_Control_Parameter);
      if (strcmp(IP.Viscous_Flux_Evaluation_Type, "Arithmetic_Mean") == 0) {
	IP.i_Viscous_Flux_Evaluation = VISCOUS_RECONSTRUCTION_ARITHMETIC; 
      } else if (strcmp(IP.Viscous_Flux_Evaluation_Type, "Mean_Gradient") == 0) {
	IP.i_Viscous_Flux_Evaluation = VISCOUS_RECONSTRUCTION_MEAN_GRADIENT;
      } else if (strcmp(IP.Viscous_Flux_Evaluation_Type, "Cartesian") == 0) {
	IP.i_Viscous_Flux_Evaluation = VISCOUS_RECONSTRUCTION_CARTESIAN; 
      } else if (strcmp(IP.Viscous_Flux_Evaluation_Type,"Diamond_Path_Least_Squares") == 0) {
	IP.i_Viscous_Flux_Evaluation = VISCOUS_RECONSTRUCTION_DIAMONDPATH_LEAST_SQUARES;
	//Must have regular reconstruction set to least squares as well.
	IP.i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES; 
      } else if (strcmp(IP.Viscous_Flux_Evaluation_Type,"Diamond_Path_Green_Gauss") == 0) {
	IP.i_Viscous_Flux_Evaluation = VISCOUS_RECONSTRUCTION_DIAMONDPATH_GREEN_GAUSS;
	//Must have regular reconstruction set to green_gauss as well.
	IP.i_Reconstruction = RECONSTRUCTION_GREEN_GAUSS;
      } else {
	cout<<"\n Not a valid Viscous Flux Evaluation Type, Using Diamond_Path_Green_Gauss\n";
	i_command = INVALID_INPUT_VALUE;	
	IP.i_Viscous_Flux_Evaluation = VISCOUS_RECONSTRUCTION_DIAMONDPATH_GREEN_GAUSS; 	
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
     
    }else if (strcmp(IP.Next_Control_Parameter, "Wall_Boundary_Treatments") == 0) {
       i_command = 516;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Wall_Boundary_Treatments;
       IP.Input_File.getline(buffer, sizeof(buffer));
        
//     }else if (strcmp(IP.Next_Control_Parameter, "Stretch_Level") == 0) {
//        i_command = 517;
//        IP.Line_Number = IP.Line_Number + 1;
//        IP.Input_File >> IP.Stretch_Level;
//        IP.Input_File.getline(buffer, sizeof(buffer));
       
    }else if (strcmp(IP.Next_Control_Parameter, "Pressure_Gradient") == 0) {
       i_command = 518;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Pressure_Gradient;
       IP.Input_File.getline(buffer, sizeof(buffer));
    }else if (strcmp(IP.Next_Control_Parameter, "Reynolds_Number") == 0) {
       i_command = 519;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Reynolds_Number;
       IP.Input_File.getline(buffer, sizeof(buffer));
    }else if (strcmp(IP.Next_Control_Parameter, "BluffBody_Data_Usage") == 0) {
       i_command = 520;
       IP.BluffBody_Data_Usage = 1;
       
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

    } else if (strcmp(IP.Next_Control_Parameter, "Residual_Variable") == 0) {
      i_command = 666;
      IP.Line_Number = IP.Line_Number + 1;                   //should change to read (rho, momentum, energy)
      IP.Input_File >> IP.i_Residual_Variable;                     //as apposed to 1, 2, or 4
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.i_Residual_Variable < ZERO) i_command = INVALID_INPUT_VALUE;

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
   
    } else if (strcmp(IP.Next_Control_Parameter, "Morton") == 0) {
      i_command = 83;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Morton;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Morton < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Morton_Reordering_Frequency") == 0) {
      i_command = 84;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Morton_Reordering_Frequency;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Morton_Reordering_Frequency < 0) i_command = INVALID_INPUT_VALUE;

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
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.i_Mesh_Stretching = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
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
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.i_Smooth_Quad_Block = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.i_Smooth_Quad_Block = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter, "AMR") == 0) {
      i_command = 72;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.AMR = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
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

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Uniform_Mesh_Coarsenings") == 0) {
      i_command = 740;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Number_of_Uniform_Mesh_Coarsenings;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Number_of_Uniform_Mesh_Coarsenings < 0) IP.Number_of_Uniform_Mesh_Coarsenings = 0;

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
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Refinement_Criteria_Gradient_Density = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.Refinement_Criteria_Gradient_Density = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Divergence_Velocity") == 0) {
      i_command = 83;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Refinement_Criteria_Divergence_Velocity = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.Refinement_Criteria_Divergence_Velocity = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Curl_Velocity") == 0) {
      i_command = 84;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Refinement_Criteria_Curl_Velocity = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.Refinement_Criteria_Curl_Velocity = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Gradient_Temperature") == 0) {
      i_command = 85;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Refinement_Criteria_Gradient_Temperature = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.Refinement_Criteria_Gradient_Temperature = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Gradient_CH4") == 0) {
      i_command = 86;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Refinement_Criteria_Gradient_CH4 = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.Refinement_Criteria_Gradient_CH4 = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Gradient_CO2") == 0) {
      i_command = 87;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Refinement_Criteria_Gradient_CO2 = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
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
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.i_Mesh_Stretching = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
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
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.i_Smooth_Quad_Block = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.i_Smooth_Quad_Block = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

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
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Multigrid_IP.Defect_Correction = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.Multigrid_IP.Defect_Correction = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"First_Order_Coarse_Mesh_Reconstruction") == 0) {
      i_command = 205;
      IP.Line_Number = IP.Line_Number + 1;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Multigrid_IP.First_Order_Coarse_Mesh_Reconstruction = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.Multigrid_IP.First_Order_Coarse_Mesh_Reconstruction = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Prolong_Using_Injection") == 0) {
      i_command = 206;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Multigrid_IP.Prolong_Using_Injection = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.Multigrid_IP.Prolong_Using_Injection = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Apply_Coarse_Mesh_Boundary_Conditions") == 0) {
      i_command = 207;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.Multigrid_IP.Apply_Coarse_Mesh_Boundary_Conditions = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Injection_at_Dirichlet_Boundary_Conditions") == 0) {
      i_command = 208;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Multigrid_IP.Injection_at_Dirichlet_Boundary_Conditions = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.Multigrid_IP.Injection_at_Dirichlet_Boundary_Conditions = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Update_Stability_Switch") == 0) {
      i_command = 209;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Multigrid_IP.Update_Stability_Switch = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
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
  
    } else if (strcmp(IP.Next_Control_Parameter, "Multigrid_Absolute_Convergence_Tolerance") == 0) {
       i_command = 215;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Multigrid_IP.Absolute_Convergence_Tolerance;
       IP.Input_File.getline(buffer, sizeof(buffer));
  
    } else if (strcmp(IP.Next_Control_Parameter, "Multigrid_Relative_Convergence_Tolerance") == 0) {
       i_command = 216;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Multigrid_IP.Relative_Convergence_Tolerance;
       IP.Input_File.getline(buffer, sizeof(buffer));
  
    } else if (strcmp(IP.Next_Control_Parameter, "FMG_Absolute_Convergence_Tolerance") == 0) {
       i_command = 217;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Multigrid_IP.FMG_Absolute_Convergence_Tolerance;
       IP.Input_File.getline(buffer, sizeof(buffer));
  
    } else if (strcmp(IP.Next_Control_Parameter, "FMG_Relative_Convergence_Tolerance") == 0) {
       i_command = 218;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Multigrid_IP.FMG_Relative_Convergence_Tolerance;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter,"Multigrid_Smoothing_Type") == 0) {
      i_command = 219;
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
      i_command = 220;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Multigrid_IP.Ncycles_Regular_Multigrid;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Multigrid_IP.Ncycles_Regular_Multigrid < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Ncycles_Full_Multigrid") == 0) {
      i_command = 221;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Multigrid_IP.Ncycles_Full_Multigrid;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Multigrid_IP.Ncycles_Full_Multigrid < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Physical_Time_Integration_Type") == 0) {
      i_command = 222;
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
      i_command = 223;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Multigrid_IP.Physical_Time_CFL_Number;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Multigrid_IP.Physical_Time_CFL_Number <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Dual_Time_Convergence_Residual_Level") == 0) {
      i_command = 224;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Multigrid_IP.Dual_Time_Convergence_Residual_Level;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Multigrid_IP.Dual_Time_Convergence_Residual_Level < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Multigrid_Write_Output_Cells_Frequency") == 0) {
      i_command = 225;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Multigrid_IP.Write_Output_Cells_Frequency;
      IP.Input_File.getline(buffer,sizeof(buffer));

       ////////////////////////////////////////////////////////////////////
       // SPECIFIED BOUNDARY CONDITIONS                                  //
       ////////////////////////////////////////////////////////////////////

    } else if (strcmp(IP.Next_Control_Parameter,"Boundary_Conditions_Specified") == 0) {
      i_command = 500;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Boundary_Conditions_Specified,IP.Next_Control_Parameter);
      if (strcmp(IP.Boundary_Conditions_Specified,"ON") == 0) {
	IP.BCs_Specified = ON;
      } else if (strcmp(IP.Boundary_Conditions_Specified,"OFF") == 0) {
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

 //    }else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Multigrid") == 0) {
//        i_command = WRITE_OUTPUT_MULTIGRID_CODE;
       
//     } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Multigrid_Cells") == 0) {
//        i_command = WRITE_OUTPUT_MULTIGRID_CELLS_CODE;
       
    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_RHS") == 0) {
       i_command = WRITE_OUTPUT_RHS_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Mesh") == 0) {
       i_command = WRITE_OUTPUT_GRID_CODE;

    }else if (strcmp(IP.Next_Control_Parameter, "Perturbation") == 0) {
       i_command = WRITE_OUTPUT_PERTURB_CODE;

    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Driven_Cavity_Flow") == 0) {
       i_command = WRITE_OUTPUT_DRIVEN_CAVITY_FLOW_CODE;
       
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

    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Driven_Cavity_Flow") == 0) {
      i_command = WRITE_OUTPUT_DRIVEN_CAVITY_FLOW_CODE;

    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Quasi3D") == 0) {
      i_command = WRITE_OUTPUT_QUASI3D_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Refine_Grid") == 0) {
       i_command = REFINE_GRID_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Morton_Ordering") == 0) {
       i_command = MORTON_ORDERING_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Fix_BCs") == 0) {
      i_command = SWITCH_BCS_TO_FIXED;

    } else if (strcmp(IP.Next_Control_Parameter, "Postprocess_1Dflame") == 0) {
       i_command = POSTPROCESS_1DFLAME;

    } else if (strcmp(IP.Next_Control_Parameter, "Postprocess_Turbulence") == 0) {
       i_command = POSTPROCESS_TURBULENCE;

    } else if (IP.Next_Control_Parameter[0] == '#') {
       i_command = COMMENT_CODE;
    } else {
       i_command = INVALID_INPUT_CODE;    
    } /* endif */

    if (i_command == INVALID_INPUT_CODE) {
      
      strcpy(buffer, IP.Next_Control_Parameter);
      Get_Next_Input_Control_Parameter(IP);
      i_command = IP.NKS_IP.Parse_Next_Input_Control_Parameter(buffer, IP.Next_Control_Parameter);
      
      // If it's still unknown then ignore it. 
      // This could be a bad idea if it was an unknown command 
      // as opposed to an unknown code.
      if (i_command == INVALID_INPUT_CODE) {
	cout << "\n***\n\nWarning: input file line " << IP.Line_Number << ": ";
	cout << "ignoring unknown input code:\n";
	cout << "code: " << buffer;
	cout << "\nvalue: " << IP.Next_Control_Parameter;
	cout << "\n\n***\n";
      }
      i_command = COMMENT_CODE; // sure why not
    }
    
    if (!IP.Input_File.good()) { i_command = INVALID_INPUT_VALUE; }
       
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
int Process_Input_Control_Parameter_File(LESPremixed2D_Input_Parameters &Input_Parameters,
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
    error_flag = Input_Parameters.Input_File.fail();

    if (error_flag) {
       cout << "\n LESPremixed2D ERROR: Unable to open LESPremixed2D input data file.\n";
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
          cout << "\n LESPremixed2D ERROR: Error reading LESPremixed2D data at line #"
               << -line_number  << " of input data file.\n";
          error_flag = line_number;
          break;
       } /* endif */
    } /* endwhile */

    /* Initial processing of input control parameters complete.  
       Return the error indicator flag. */
   
    //Load the C-Type strings from the C++ strings 
    strcpy(Input_Parameters.React_Name,Input_Parameters.react_name.c_str());
    strcpy(Input_Parameters.ct_Mech_Name,Input_Parameters.ct_mech_name.c_str());
    strcpy(Input_Parameters.ct_Mech_File,Input_Parameters.ct_mech_file.c_str());
    strcpy(Input_Parameters.Scalar_system_name,Input_Parameters.scalar_system_name.c_str());
    
    for (int i = 0; i < Input_Parameters.num_species; i++) {
      strcpy(Input_Parameters.Multispecies[i],Input_Parameters.multispecies[i].c_str());
    }
    if (Input_Parameters.num_scalars > 0) {
      for (int i = 0; i < Input_Parameters.num_scalars; i++) {
        strcpy(Input_Parameters.Scalars[i],Input_Parameters.scalars[i].c_str());
      }
    }

    // Proper temperature for display
    Input_Parameters.Temperature = Input_Parameters.Wo.T();
 
    // Reset the static variables.
    Input_Parameters.Wo.set_SFSmodel_variables(Input_Parameters.filter_width,
					       Input_Parameters.Smagorinsky_Constant,
					       Input_Parameters.Yoshizawa_coefficient);
    Input_Parameters.Uo.set_SFSmodel_variables(Input_Parameters.filter_width,
					       Input_Parameters.Smagorinsky_Constant,
					       Input_Parameters.Yoshizawa_coefficient);

    Input_Parameters.Wo.set_premixed_flame_variables(Input_Parameters.laminar_flame_thickness,
						     Input_Parameters.laminar_flame_speed,
						     Input_Parameters.TFactor,
                                                     Input_Parameters.adiabatic_temp,
                                                     Input_Parameters.equivalence_ratio,
                                                     Input_Parameters.reactants_den);
    Input_Parameters.Uo.set_premixed_flame_variables(Input_Parameters.laminar_flame_thickness,
						     Input_Parameters.laminar_flame_speed,
						     Input_Parameters.TFactor,
                                                     Input_Parameters.adiabatic_temp,
                                                     Input_Parameters.equivalence_ratio,
                                                     Input_Parameters.reactants_den);

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
 *                                                      *
 ********************************************************/
// void Equivalence_Ratio(const double &phi){

//   // Methane & Air
//   // CH4 + 2(O2 +3.76N2) -> CO2 + 2H2O 
//   double n = 9.52*(1.0 - phi)/(phi +9.52);

//   double X_CH4 = ( 1.0 - n)/10.52;
//   double X_O2  = (2.0 + 0.21*n)/10.52;
//   double X_N2  = (7.52 +0.79*n)/10.52;

//   double Mtot = X_CH4*(16.0) + X_O2*(32.0) + X_N2*(28.0);

//   double c_CH4 =  X_CH4*(16.0)/Mtot;
//   double c_O2  =  X_O2*(32.0)/Mtot;
//   double c_N2  =  X_N2*(28.0)/Mtot;

//   cout<<"\n For phi = "<<phi;
//   cout<<"\n c_CH4 "<<c_CH4;
//   cout<<"\n c_O2  "<<c_O2;
//   cout<<"\n c_N2  "<<c_N2;
//   cout<<"\n X_CH4 "<<X_CH4;
//   cout<<"\n X_O2  "<<X_O2;
//   cout<<"\n X_N2  "<<X_N2;
//   cout<<"\n Mtot  "<<Mtot;

//   cout<<endl;
//   cout<<"\n Should pass in equation type, and pass back the appropriate mass fractions ";
//   cout<<endl;
  
//   exit(0); 

// }
