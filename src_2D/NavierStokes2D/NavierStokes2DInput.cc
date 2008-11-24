/*!\file NavierStokes2DInput.cc
  \brief Subroutines for the 2D Navier-Stokes Input Class. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "NavierStokes2DInput.h"   // Include 2D NavierStokes input parameter header file.
#include "NavierStokes2DQuad.h"    /* Include NavierStokes2D_Quad_Block header file. */
#include "../Grid/HO_Grid2DQuad_ExecutionMode.h" // Include high-order 2D grid execution mode header file
#include "../Grid/Tecplot_ExecutionMode.h" // Include Tecplot execution mode header file

/*********************************************************
 * NavierStokes2D_Input_Parameters -- Member functions.  *
 ********************************************************/

/******************************************************//**
 * Parse the input file
 ********************************************************/
int NavierStokes2D_Input_Parameters::Parse_Input_File(char *Input_File_Name_ptr){

  ostringstream msg;
  int command_flag, error_flag;

  /* Assign initial value for error indicator flag. */
  error_flag = 0;

  strcpy(Input_File_Name, Input_File_Name_ptr);
  Open_Input_File(*this);
  if (Input_File.fail()) {
    msg << "NavierStokes2D_Input_Parameters::Parse_Input_File() ERROR: Unable to open "
	<<string(Input_File_Name_ptr) 
	<< " input data file.";
    if (Verbose()) {
      cerr << msg.str() << endl;
    }
    throw runtime_error(msg.str());
  } /* endif */

  if (Verbose()) {
    cout << "\n Reading input data file `"
	 << Input_File_Name << "'." << endl;
    cout.flush();
  }
  while (1) {
    Get_Next_Input_Control_Parameter();
    command_flag = Parse_Next_Input_Control_Parameter(*this);
    if (command_flag == EXECUTE_CODE) {
      break;
      
    } else if (command_flag == TERMINATE_CODE) {
      break;
      
    } else if (command_flag == INVALID_INPUT_CODE ||
	       command_flag == INVALID_INPUT_VALUE) {
      Line_Number = -Line_Number;
      
      msg << "NavierStokes2D_Input_Parameters::Parse_Input_File() ERROR: Error reading data at line # " 
	  << -Line_Number
	  << " of input data file.";
      if (Verbose()){
	cerr << msg.str() << endl;
      }
      
      throw runtime_error(msg.str());
    } /* endif */
  } /* endwhile */

  /* Perform consistency checks and internal parameter setup */
  doInternalSetupAndConsistencyChecks(error_flag);

  /* Initial processing of input control parameters complete.  
     Return the error indicator flag. */
  return (error_flag);

}

/******************************************************//**
 * Get the next input control parameter from the input file.                                                
 ********************************************************/
void NavierStokes2D_Input_Parameters::Get_Next_Input_Control_Parameter(void){

  int i, index, LineSize, IndexFirstChar(0);
  char buffer[256], ControlParameter[256];

  // Initialize ControlParameter and Next_Control_Parameter to end of string
  ControlParameter[0] = '\0';
  strcpy(Next_Control_Parameter, ControlParameter);

  // While the input stream is 'good' for reading and the end of file is not attained
  while ( Input_File.good() && !Input_File.getline(buffer, sizeof(buffer)).eof() ){

    // Process the line 
    Line_Number = Line_Number + 1;
    LineSize = Input_File.gcount(); // Get the size of the line. Last character is "\0"!

    // Determine the index of the first character different than 'space' and 'tab'
    for (i=0; i<LineSize; ++i){
      if (buffer[i] != ' ' && buffer[i] != '\t'){
	IndexFirstChar = i;
	break;
      }
    }

    /* Parse the line if the first character different than 'space' 
       is also different than '#' or end of string ('\0').
       Otherwise skip the line because it is either a comment or an empty line. */
    if ( buffer[IndexFirstChar] != '#' && buffer[IndexFirstChar] != '\0'){

      // Get the ControlParameter
      for(i=IndexFirstChar, index=0;  i<LineSize;  ++i, ++index){
	if (buffer[i] == ' ' || buffer[i] == '=' || buffer[i] == '\t'){
	  ControlParameter[index] = '\0';
	  break;
	} else {
	  ControlParameter[index] = buffer[i];
	}
      }

      // Set the Next_Control_Parameter
      strcpy(Next_Control_Parameter, ControlParameter);
      break;
    }

  }//endwhile
}

/******************************************************//**
 * Perform setup of internal and related external parameters.
 * Do also check and validation of input parameters.
 ********************************************************/
void NavierStokes2D_Input_Parameters::doInternalSetupAndConsistencyChecks(int & error_flag){

  /* Perform consistency checks on input parameters for multigrid */
  if (!error_flag &&
      i_Time_Integration == TIME_STEPPING_MULTIGRID) {
    error_flag = Check_Input_Parameters<NavierStokes2D_Input_Parameters>(*this);
    if (error_flag) {
      cout << "\n NavierStokes2D ERROR: Input Parameters consistency check failure\n";
    }
  }

  // Perform consitency checks on the refinement criteria.
  Number_of_Refinement_Criteria = 0;
  if (Refinement_Criteria_Gradient_Density) Number_of_Refinement_Criteria++;
  if (Refinement_Criteria_Divergence_Velocity) Number_of_Refinement_Criteria++;
  if (Refinement_Criteria_Curl_Velocity) Number_of_Refinement_Criteria++;
  if (Refinement_Criteria_Gradient_Turbulence_Kinetic_Energy) Number_of_Refinement_Criteria++;
  if (Number_of_Refinement_Criteria < 1 || Number_of_Refinement_Criteria > 4){
    error_flag = 1011;
  }

  // Perform consitency checks on the time marching parameters.
  if (Time_Accurate == 1 && Local_Time_Stepping != GLOBAL_TIME_STEPPING){
    Local_Time_Stepping = GLOBAL_TIME_STEPPING;
  }

  // Enforce NO mesh stretching is required
  if (!i_Mesh_Stretching){
    // Mesh stretching is not ON
    Mesh_Stretching_Type_Idir = STRETCHING_FCN_LINEAR;
    Mesh_Stretching_Type_Jdir = STRETCHING_FCN_LINEAR;
    Mesh_Stretching_Factor_Idir = 1.0;
    Mesh_Stretching_Factor_Jdir = 1.0;
  }

  // Perform update of the internal variables of the exact solution
  ExactSoln->Set_ParticularSolution_Parameters(*this);

  // Perform update of the internal variables of the high-order input parameters
  HighOrder2D_Input::Set_Final_Parameters(*this);

  // Set reference state in the NavierStokes2D_Quad_Block class
  NavierStokes2D_Quad_Block::Set_Normalization_Reference_State(RefW);
  
  // Set limiter in CENO class
  CENO_Execution_Mode::Limiter = i_Limiter;

}

/**********************************************************************
 * NavierStokes2D_Input_Parameters -- External subroutines.           *
 **********************************************************************/

/**********************************************************************
 * Routine: Open_Input_File                                           *
 *                                                                    *
 * Opens the appropriate input data file.                             *
 *                                                                    *
 **********************************************************************/
void Open_Input_File(NavierStokes2D_Input_Parameters &IP) {

  IP.Input_File.open(IP.Input_File_Name,ios::in);
  if (!IP.Input_File.fail()) {
    IP.Line_Number = 0;
    IP.Input_File.setf(ios::skipws);
  }

}

/**********************************************************************
 * Routine: Close_Input_File                                          *
 *                                                                    *
 * Closes the appropriate input data file.                            *
 *                                                                    *
 **********************************************************************/
void Close_Input_File(NavierStokes2D_Input_Parameters &IP) {

  IP.Input_File.unsetf(ios::skipws);
  IP.Input_File.close();

}

/**********************************************************************
 * Routine: Set_Default_Input_Parameters                              *
 *                                                                    *
 * Assigns default values to the input parameters.                    *
 *                                                                    *
 **********************************************************************/
void Set_Default_Input_Parameters(NavierStokes2D_Input_Parameters &IP) {

  char *string_ptr;
  double cos_angle, sin_angle;
  
  // CFFC root directory path:
  IP.get_cffc_path();

  string_ptr = "NavierStokes2D.in";
  strcpy(IP.Input_File_Name,string_ptr);

  // Time-stepping parameters:
  string_ptr = "Explicit_Euler";
  strcpy(IP.Time_Integration_Type,string_ptr);
  IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
  IP.Time_Accurate = 0;
  IP.Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
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
  strcpy(IP.Reconstruction_Type,string_ptr);
  IP.i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;

  // Limiter type:
  string_ptr = "Barth_Jespersen";
  strcpy(IP.Limiter_Type,string_ptr);
  IP.i_Limiter = LIMITER_BARTH_JESPERSEN;
  IP.Freeze_Limiter = 0;
  IP.Freeze_Limiter_Residual_Level = 1e-4;

  // Gas-phase flux function:
  string_ptr = "Roe";
  strcpy(IP.Flux_Function_Type,string_ptr);
  IP.i_Flux_Function = FLUX_FUNCTION_ROE;

  // Viscous gradient reconstruction type:
  string_ptr = "Diamond_Path";
  strcpy(IP.Viscous_Reconstruction_Type,string_ptr);
  IP.i_Viscous_Reconstruction = VISCOUS_RECONSTRUCTION_DIAMOND_PATH;

  // Initial conditions:
  string_ptr = "Uniform";
  strcpy(IP.ICs_Type,string_ptr);
  IP.i_ICs = IC_UNIFORM;

  // Flow-type switch:
  string_ptr = "Inviscid";
  strcpy(IP.Flow_Type,string_ptr);
  IP.FlowType = FLOWTYPE_LAMINAR;

  // Geometry switch:
  string_ptr = "Planar";
  strcpy(IP.Flow_Geometry_Type,string_ptr);
  IP.Axisymmetric = OFF;

  // Turbulence parameters:
  string_ptr = "Direct_Integration";
  strcpy(IP.Turbulence_BC_Type,string_ptr);
  IP.i_Turbulence_BCs = TURBULENT_BC_DIRECT_INTEGRATION;
  string_ptr = "Wall";
  strcpy(IP.Friction_Velocity_Type,string_ptr);
  IP.i_Friction_Velocity = FRICTION_VELOCITY_WALL_SHEAR_STRESS;
  IP.C_constant = 5.0;
  IP.von_Karman_Constant = 0.41;
  IP.yplus_sublayer = 2.5;
  IP.yplus_buffer_layer = 30.0;
  IP.yplus_outer_layer = 250.0;//If this comment exists, you yet need to verify this value (6th sept, 2007)
  IP.n_cells_sublayer_support = 5;
  IP.i_Turbulent_Wall_Injection = OFF;
  IP.sigmav = 0.035;
  IP.lw = 0.000200;
  
  // Compressibility effect switch:
  string_ptr = "Off";
  strcpy(IP.Compressibility_Effect,string_ptr);
  IP.i_Compressibility_Effect = OFF;

  // Transition model:
  string_ptr = "Off";
  strcpy(IP.Transition_Model,string_ptr);
  IP.i_Transition_Model = OFF;

  // Variable Prandtl number model:
  string_ptr = "Off";
  strcpy(IP.Variable_Prandtl,string_ptr);
  IP.i_Variable_Prandtl = OFF;
  IP.C_lambda = 0.14; IP.Cd1 = 2.0; IP.Cd4 = 2.2; IP.Cd5 = 0.8; 
 
  // Propellant type:
  string_ptr = "AP_HTPB";
  strcpy(IP.Propellant_Type,string_ptr);

  // Gas-phase conditions:
  string_ptr = "AIR";
  strcpy(IP.Gas_Type,string_ptr);
  IP.Wo.set_static_variables(IP.Gas_Type,IP.FlowType,IP.C_constant,IP.von_Karman_Constant,
			     IP.yplus_sublayer,IP.yplus_buffer_layer,
			     IP.yplus_outer_layer,IP.Propellant_Type,
			     IP.i_Compressibility_Effect, IP.i_Transition_Model, IP.i_Variable_Prandtl,
			     IP.sigmav,IP.lw,IP.C_lambda,IP.Cd1,IP.Cd4,IP.Cd5);
  IP.Pressure = IP.Wo.p;
  IP.Temperature = IP.Wo.T();
  IP.Mach_Number = 0.80;
  IP.Mach_Number2 = ZERO;
  IP.Flow_Angle = ZERO;
  cos_angle = cos(TWO*PI*IP.Flow_Angle/360.00); if (fabs(cos_angle) < TOLER*TOLER) cos_angle = ZERO;
  sin_angle = sin(TWO*PI*IP.Flow_Angle/360.00); if (fabs(sin_angle) < TOLER*TOLER) sin_angle = ZERO;
  IP.Wo.v.x = IP.Mach_Number*IP.Wo.a()*cos_angle;
  IP.Wo.v.y = IP.Mach_Number*IP.Wo.a()*sin_angle;
  IP.Reynolds_Number = ZERO;
  IP.dp = ZERO;
  IP.Re_lid = 100.0;
  IP.Wave_Position = Vector2D_ZERO;
  IP.Wave_Width = ZERO;

  // State conditions:
  IP.Uo.set_static_variables(IP.Gas_Type,IP.FlowType,IP.C_constant,IP.von_Karman_Constant,
			     IP.yplus_sublayer,IP.yplus_buffer_layer,
			     IP.yplus_outer_layer,IP.Propellant_Type,
			     IP.i_Compressibility_Effect, IP.i_Transition_Model, IP.i_Variable_Prandtl,
			     IP.sigmav,IP.lw,IP.C_lambda,IP.Cd1,IP.Cd4,IP.Cd5);
  IP.Uo = U(IP.Wo);
  IP.W1 = IP.Wo;
  IP.W2 = IP.Wo;

  // Grid parameters:
  string_ptr = "Square";
  strcpy(IP.Grid_Type,string_ptr);
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
  IP.Pipe_Radius = HALF;
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
  IP.Inner_Streamline_Number = 0.80;
  IP.Outer_Streamline_Number = 0.40;
  IP.Isotach_Line = 0.30;
  IP.Wedge_Angle = 25.0;
  IP.Wedge_Length = HALF;
  IP.Step_Height = 0.0127;
  IP.Top_Wall_Deflection = ZERO;
  IP.Smooth_Bump = OFF;

  IP.VertexSW = Vector2D(-0.5,-0.5);
  IP.VertexSE = Vector2D( 0.5,-0.5);
  IP.VertexNE = Vector2D( 0.5, 0.5);
  IP.VertexNW = Vector2D(-0.5, 0.5);

  IP.X_Shift = Vector2D_ZERO;
  IP.X_Scale = ONE;
  IP.X_Rotate = ZERO;

  IP.IterationsOfInteriorNodesDisturbances = 0;     /* Number of iterations of disturbing the mesh 
						       (create an unsmooth interior mesh) */
  IP.Num_Of_Spline_Control_Points = 361; /* Number of control points on the 2D spline (used for some grids) */

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

  // NASA rotor input variables:
  strcpy(IP.NASA_Rotor37_Data_Directory, IP.CFFC_Path);
  strcpy(IP.NASA_Rotor37_Data_Directory, "/data/NASA_Rotors/R37/");
  strcpy(IP.NASA_Rotor67_Data_Directory, IP.CFFC_Path);
  strcpy(IP.NASA_Rotor67_Data_Directory, "/data/NASA_Rotors/R67/");
  IP.Rotor_Flow_Type = PEAK_FLOW;
  IP.Rotor_Percent_Span = 50.00;

  // Viscous channel flow input variables:
  IP.Vwall = Vector2D(0.0,0.0);
  IP.Twall = TEMPERATURE_STDATM;

  // ICEM:
  IP.ICEMCFD_FileNames = ICEMCFD_get_filenames();

  // AMR:
  IP.AMR = 0;
  IP.AMR_Frequency = 100;
  IP.Number_of_Initial_Mesh_Refinements = 0;
  IP.Number_of_Uniform_Mesh_Refinements = 0;
  IP.Number_of_Boundary_Mesh_Refinements = 0;
  IP.Number_of_Interface_Mesh_Refinements = 0;
  IP.Number_of_Bounding_Box_Mesh_Refinements = 0;
  IP.Number_of_Flat_Plate_Mesh_Refinements = 0;
  IP.Interface_Refinement_Condition = OFF;
  IP.Maximum_Refinement_Level = 100;
  IP.Minimum_Refinement_Level = 1;
  IP.Threshold_for_Refinement = 0.50;
  IP.Threshold_for_Coarsening = 0.10;
  IP.Number_of_Refinement_Criteria = 3;
  IP.Refinement_Criteria_Gradient_Density = ON;
  IP.Refinement_Criteria_Divergence_Velocity = ON;
  IP.Refinement_Criteria_Curl_Velocity = ON;
  IP.Refinement_Criteria_Gradient_Turbulence_Kinetic_Energy = OFF;
  IP.AMR_Xmin = Vector2D_ZERO;
  IP.AMR_Xmax = Vector2D_ZERO;
  IP.Morton = 0;
  IP.Morton_Reordering_Frequency = 0;

  // Smooth quad block indicator:
  IP.i_Smooth_Quad_Block = ON;

  IP.Solver_Type = EXPLICIT;

  // Embedded boundary input parameters:
  IP.Reset_Interface_Motion_Type = OFF;

  // Default output file names and parameters:
  string_ptr = "outputfile.dat";
  strcpy(IP.Output_File_Name,string_ptr);
  string_ptr = "gridfile.grid";
  strcpy(IP.Grid_File_Name,string_ptr);
  string_ptr = "gridfile.griddef";
  strcpy(IP.Grid_Definition_File_Name,string_ptr);
  string_ptr = "restartfile.soln";
  strcpy(IP.Restart_File_Name,string_ptr);
  string_ptr = "gnuplotfile.gplt";
  strcpy(IP.Gnuplot_File_Name,string_ptr);
  string_ptr = "Tecplot";
  strcpy(IP.Output_Format_Type,string_ptr);
  IP.i_Output_Format = IO_TECPLOT;
  IP.Restart_Solution_Save_Frequency = 1000;

  // Default output progress frequency:
  IP.Output_Progress_Frequency = 50;

  // Input_file parameters:
  string_ptr = " ";
  strcpy(IP.Next_Control_Parameter,string_ptr);
  IP.Line_Number = 0;
  IP.Number_of_Processors = CFFC_MPI::Number_of_Processors;
  IP.Number_of_Blocks_Per_Processor = 10;

}

/**********************************************************************
 * Routine: Broadcast_Input_Parameters                                *
 *                                                                    *
 * Broadcast the input parameters variables to all processors in the  *
 * calculation from the primary processor using the MPI broadcast     *
 * routine.                                                           *
 *                                                                    *
 **********************************************************************/
void Broadcast_Input_Parameters(NavierStokes2D_Input_Parameters &IP) {

#ifdef _MPI_VERSION

  double cos_angle, sin_angle;

  // CFFC path:
  MPI::COMM_WORLD.Bcast(IP.CFFC_Path, 
 			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR, 0);
  // Input file name and line number:
  MPI::COMM_WORLD.Bcast(IP.Input_File_Name,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.Line_Number),
			1,
			MPI::INT,0);
  // Flow type:
  MPI::COMM_WORLD.Bcast(IP.Flow_Type,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.FlowType),
			1,
			MPI::INT,0);
  // Flow geometry type:
  MPI::COMM_WORLD.Bcast(IP.Flow_Geometry_Type,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.Axisymmetric),
			1,
			MPI::INT,0);
  // Time integration:
  MPI::COMM_WORLD.Bcast(IP.Time_Integration_Type,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Time_Integration),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Time_Accurate),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Local_Time_Stepping),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Maximum_Number_of_Time_Steps),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.N_Stage),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.CFL_Number),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Time_Max),
			1,
			MPI::DOUBLE,0);
  // Residual variable:
  MPI::COMM_WORLD.Bcast(&(IP.i_Residual_Variable),
			1,
			MPI::DOUBLE,0);
  // Residual smoothing:
  MPI::COMM_WORLD.Bcast(&(IP.Residual_Smoothing),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Residual_Smoothing_Epsilon),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Residual_Smoothing_Gauss_Seidel_Iterations),
			1,
			MPI::INT,0);
  // Multigrid related parameters:
  IP.Multigrid_IP.Broadcast_Input_Parameters();
  IP.NKS_IP.Broadcast_Input_Parameters();
  // Reconstruction:
  MPI::COMM_WORLD.Bcast(IP.Reconstruction_Type,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Reconstruction),
			1,
			MPI::INT,0);
  // Limiters:
  MPI::COMM_WORLD.Bcast(IP.Limiter_Type,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Limiter),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Freeze_Limiter),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Freeze_Limiter_Residual_Level),
			1,
			MPI::DOUBLE,0);
  // Flux functions:
  MPI::COMM_WORLD.Bcast(IP.Flux_Function_Type,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Flux_Function),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(IP.Viscous_Reconstruction_Type,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Viscous_Reconstruction),
			1,
			MPI::INT,0);
  // Turbulence parameters:
  MPI::COMM_WORLD.Bcast(IP.Turbulence_BC_Type,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Turbulence_BCs),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(IP.Friction_Velocity_Type,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
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
  MPI::COMM_WORLD.Bcast(&(IP.n_cells_sublayer_support),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Turbulent_Wall_Injection),
			1,
			MPI::INT,0);
  // Propellant type:
  MPI::COMM_WORLD.Bcast(IP.Propellant_Type,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.sigmav),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.lw),
			1,
			MPI::DOUBLE,0);
  // Compressibility Effects:
  MPI::COMM_WORLD.Bcast(IP.Compressibility_Effect,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Compressibility_Effect),
			1,
			MPI::INT,0);
  // Transition Model
  MPI::COMM_WORLD.Bcast(IP.Transition_Model,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.Transition_Model),
			1,
			MPI::INT,0);
  // Variable Prandtl number
  MPI::COMM_WORLD.Bcast(IP.Variable_Prandtl,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Variable_Prandtl),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.C_lambda),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Cd1),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Cd4),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Cd5),
			1,
			MPI::DOUBLE,0);
  // Initial conditions:
  MPI::COMM_WORLD.Bcast(IP.ICs_Type,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.Gas_Type,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_ICs),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Pressure),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Temperature),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Mach_Number),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Mach_Number2),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Flow_Angle),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Reynolds_Number),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.dp),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Re_lid),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Wave_Position.x),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Wave_Position.y),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Wave_Width),
			1,
			MPI::DOUBLE,0);
  if (!CFFC_Primary_MPI_Processor()) {
    Initialize_Reference_State(IP);
  }
  for (int nv = 1; nv <= NUM_VAR_NAVIERSTOKES2D; nv++) {
    MPI::COMM_WORLD.Bcast(&(IP.W1[nv]),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.W2[nv]),
			  1,
			  MPI::DOUBLE,0);
  }
  // Grid variables:
  MPI::COMM_WORLD.Bcast(IP.Grid_Type,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.NACA_Aerofoil_Type,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Grid),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Cells_Idir),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Cells_Jdir),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Ghost_Cells),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Idir),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Jdir),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Box_Width),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Box_Height),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Plate_Length),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Pipe_Length),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Pipe_Radius),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Blunt_Body_Radius),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Blunt_Body_Mach_Number),
			1,
			MPI::DOUBLE,0);
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
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Cylinder_Radius),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Ellipse_Length_X_Axis),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Ellipse_Length_Y_Axis),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Chord_Length),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Orifice_Radius),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Inner_Streamline_Number),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Outer_Streamline_Number),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Isotach_Line),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Wedge_Angle),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Wedge_Length),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Step_Height),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Top_Wall_Deflection),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Smooth_Bump),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.X_Shift.x),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.X_Shift.y),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.X_Scale),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.X_Rotate),
			1,
			MPI::DOUBLE,0);
  // Boundary Conditions:
  MPI::COMM_WORLD.Bcast(IP.Boundary_Conditions_Specified,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.BCs_Specified),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(IP.BC_North_Type,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.BC_South_Type,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.BC_East_Type,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.BC_West_Type,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
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
  // NASA rotors:
  MPI::COMM_WORLD.Bcast(IP.NASA_Rotor37_Data_Directory,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.NASA_Rotor67_Data_Directory,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
  		MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.Rotor_Flow_Type),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Rotor_Percent_Span),
			1,
			MPI::DOUBLE,0);
  // Viscous channel flow input variables:
  MPI::COMM_WORLD.Bcast(&(IP.Vwall.x),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Vwall.y),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Twall),
			1,
			MPI::DOUBLE,0);
  // ICEM:
  if (!CFFC_Primary_MPI_Processor()) {
    IP.ICEMCFD_FileNames = new char*[3];
    for (int i = 0; i < 3; i++) {
      IP.ICEMCFD_FileNames[i] = new char[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
    }
  }
  MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[0],
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[1],
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[2],
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  // AMR & Refinement Parameters:
  MPI::COMM_WORLD.Bcast(&(IP.AMR),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.AMR_Frequency),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Initial_Mesh_Refinements),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Uniform_Mesh_Refinements),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Boundary_Mesh_Refinements),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Interface_Mesh_Refinements),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Bounding_Box_Mesh_Refinements),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Flat_Plate_Mesh_Refinements),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Interface_Refinement_Condition),
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
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Threshold_for_Coarsening),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Refinement_Criteria),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Refinement_Criteria_Gradient_Density),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Refinement_Criteria_Divergence_Velocity),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Refinement_Criteria_Curl_Velocity),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Refinement_Criteria_Gradient_Turbulence_Kinetic_Energy),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.AMR_Xmin.x),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.AMR_Xmin.y),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.AMR_Xmax.x),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.AMR_Xmax.y),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Morton),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Morton_Reordering_Frequency),
			1,
			MPI::INT,0);
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
  // Interface input parameters:
  IP.Interface_IP.Broadcast_Input_Parameters();
  MPI::COMM_WORLD.Bcast(&(IP.Reset_Interface_Motion_Type),
			1,
			MPI::INT,0);
  // File Names:
  MPI::COMM_WORLD.Bcast(IP.Output_File_Name,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.Grid_File_Name,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.Grid_Definition_File_Name,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.Restart_File_Name,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.Gnuplot_File_Name,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.Output_Format_Type,
			INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Output_Format),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Restart_Solution_Save_Frequency),
			1,
			MPI::INT,0);
  // Output progress frequency:
  MPI::COMM_WORLD.Bcast(&(IP.Output_Progress_Frequency),
			1,
			MPI::INT,0);
  // Number of processors:
  if (!CFFC_Primary_MPI_Processor()) {
    IP.Number_of_Processors = CFFC_MPI::Number_of_Processors;
  }
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Per_Processor),
			1,
			MPI::INT,0);

  // Reinitialize the reference state.
  Reinitialize_Reference_State(IP);

#endif

}

#ifdef _MPI_VERSION
/**********************************************************************
 * Routine: Broadcast_Input_Parameters                                *
 *                                                                    *
 * Broadcast the input parameters variables to all processors         *
 * associated with the specified communicator from the specified      *
 * processor using the MPI broadcast routine.                         *
 *                                                                    *
 **********************************************************************/
void Broadcast_Input_Parameters(NavierStokes2D_Input_Parameters &IP,
                                MPI::Intracomm &Communicator,
                                const int Source_CPU) {

  int Source_Rank = 0;
  double cos_angle, sin_angle;

  // CFFC path:
  Communicator.Bcast(IP.CFFC_Path, 
 		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR, Source_Rank);
  Communicator.Bcast(IP.Input_File_Name,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.Line_Number),
		     1,
		     MPI::INT,Source_Rank);
  // Flow type:
  Communicator.Bcast(IP.Flow_Type,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.FlowType),
		     1,
		     MPI::INT,Source_Rank);
  // Flow geometry type:
  Communicator.Bcast(IP.Flow_Geometry_Type,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.Axisymmetric),
		     1,
		     MPI::INT,Source_Rank);
  // Time integration:
  Communicator.Bcast(IP.Time_Integration_Type,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Time_Integration),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Time_Accurate),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Local_Time_Stepping),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Maximum_Number_of_Time_Steps),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.N_Stage),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.CFL_Number),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Time_Max),
		     1,
		     MPI::DOUBLE,Source_Rank);
  // Residual variable:
  Communicator.Bcast(&(IP.i_Residual_Variable),
		     1,
		     MPI::DOUBLE,Source_Rank);
  // Residual Smoothing:
  Communicator.Bcast(&(IP.Residual_Smoothing),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Residual_Smoothing_Epsilon),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Residual_Smoothing_Gauss_Seidel_Iterations),
		     1,
		     MPI::INT,Source_Rank);
  // Multigrid related parameters:
  IP.Multigrid_IP.Broadcast_Input_Parameters(Communicator,
					     Source_CPU);
  IP.NKS_IP.Broadcast_Input_Parameters(Communicator, Source_CPU);
  // Reconstruction:
  Communicator.Bcast(IP.Reconstruction_Type,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Reconstruction),
		     1,
		     MPI::INT,Source_Rank);
  // Limiters:
  Communicator.Bcast(IP.Limiter_Type,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Limiter),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Freeze_Limiter),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Freeze_Limiter_Residual_Level),
		     1,
		     MPI::DOUBLE,Source_Rank);
  // Flux functions:
  Communicator.Bcast(IP.Flux_Function_Type,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Flux_Function),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(IP.Viscous_Reconstruction_Type,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Viscous_Reconstruction),
		     1,
		     MPI::INT,Source_Rank);
  // Turbulence parameters:
  Communicator.Bcast(IP.Turbulence_BC_Type,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Turbulence_BCs),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(IP.Friction_Velocity_Type,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
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
  Communicator.Bcast(&(IP.n_cells_sublayer_support),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.i_Turbulent_Wall_Injection),
		     1,
		     MPI::INT,Source_Rank);
  // Propellant type:
  Communicator.Bcast(IP.Propellant_Type,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.sigmav),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.lw),
		     1,
		     MPI::DOUBLE,Source_Rank);
  //Compressibility Effects
  Communicator.Bcast(IP.Compressibility_Effect,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Compressibility_Effect),
		     1,
		     MPI::INT,Source_Rank);
  //Transition Effects
  Communicator.Bcast(IP.Transition_Model,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Transition_Model),
		     1,
		     MPI::INT,Source_Rank);
  //Variable Prandtl number
  Communicator.Bcast(IP.Variable_Prandtl,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Variable_Prandtl),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.C_lambda),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Cd1),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Cd4),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Cd5),
		     1,
		     MPI::DOUBLE,Source_Rank);
  // Initial conditions:
  Communicator.Bcast(IP.ICs_Type,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.Gas_Type,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_ICs),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Pressure),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Temperature),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Mach_Number),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Mach_Number2),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Flow_Angle),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Reynolds_Number),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.dp),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Re_lid),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Wave_Position.x),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Wave_Position.y),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Wave_Width),
		     1,
		     MPI::DOUBLE,Source_Rank);
  if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
    Initialize_Reference_State(IP);
  }
  for (int nv = 1; nv <= NUM_VAR_NAVIERSTOKES2D; nv++) {
    Communicator.Bcast(&(IP.W1[nv]),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.W2[nv]),
		       1,
		       MPI::DOUBLE,Source_Rank);
  }
  // Grid variables:
  Communicator.Bcast(IP.Grid_Type,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.NACA_Aerofoil_Type,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Grid),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Cells_Idir),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Cells_Jdir),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Ghost_Cells),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Blocks_Idir),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Blocks_Jdir),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Box_Width),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Box_Height),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Plate_Length),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Pipe_Length),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Pipe_Radius),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Blunt_Body_Radius),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Blunt_Body_Mach_Number),
		     1,
		     MPI::DOUBLE,Source_Rank);
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
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Ellipse_Length_X_Axis),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Ellipse_Length_Y_Axis),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Chord_Length),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Orifice_Radius),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Inner_Streamline_Number),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Outer_Streamline_Number),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Isotach_Line),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Wedge_Angle),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Wedge_Length),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Step_Height),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Top_Wall_Deflection),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Smooth_Bump),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.X_Shift.x),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.X_Shift.y),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.X_Scale),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.X_Rotate),
		     1,
		     MPI::DOUBLE,Source_Rank);
  // Boundary Conditions:
  Communicator.Bcast(IP.Boundary_Conditions_Specified,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.BCs_Specified),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(IP.BC_North_Type,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.BC_South_Type,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.BC_East_Type,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.BC_West_Type,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
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
  // NASA rotors:
  Communicator.Bcast(IP.NASA_Rotor37_Data_Directory,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.NASA_Rotor67_Data_Directory,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.Rotor_Flow_Type),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Rotor_Percent_Span),
		     1,
		     MPI::DOUBLE,Source_Rank);
  // Viscous channel flow input variables:
  Communicator.Bcast(&(IP.Vwall.x),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Vwall.y),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Twall),
		     1,
		     MPI::DOUBLE,Source_Rank);
  // ICEM:
  if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
    IP.ICEMCFD_FileNames = new char*[3];
    for (int i = 0; i < 3; i++) {
      IP.ICEMCFD_FileNames[i] = new char[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
    }
  }
  Communicator.Bcast(IP.ICEMCFD_FileNames[0],
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.ICEMCFD_FileNames[1],
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
                       MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.ICEMCFD_FileNames[2],
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  // AMR & Refinement Parameters:
  Communicator.Bcast(&(IP.AMR),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.AMR_Frequency),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Initial_Mesh_Refinements),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Uniform_Mesh_Refinements),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Boundary_Mesh_Refinements),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Flat_Plate_Mesh_Refinements),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Interface_Mesh_Refinements),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Bounding_Box_Mesh_Refinements),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Interface_Refinement_Condition),
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
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Threshold_for_Coarsening),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Refinement_Criteria),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Refinement_Criteria_Gradient_Density),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Refinement_Criteria_Divergence_Velocity),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Refinement_Criteria_Curl_Velocity),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Refinement_Criteria_Gradient_Turbulence_Kinetic_Energy),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.AMR_Xmin.x),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.AMR_Xmin.y),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.AMR_Xmax.x),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.AMR_Xmax.y),
		     1,
		     MPI::DOUBLE,Source_Rank);
  MPI::COMM_WORLD.Bcast(&(IP.Morton),
		     1,
		     MPI::INT,Source_Rank);
  MPI::COMM_WORLD.Bcast(&(IP.Morton_Reordering_Frequency),
		     1,
		     MPI::INT,Source_Rank);
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
  // Interface input parameters:
  IP.Interface_IP.Broadcast_Input_Parameters(Communicator,
					     Source_CPU);
  Communicator.Bcast(&(IP.Reset_Interface_Motion_Type),
		     1,
		     MPI::INT,Source_Rank);
  // File Names:
  Communicator.Bcast(IP.Output_File_Name,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.Grid_File_Name,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.Grid_Definition_File_Name,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.Restart_File_Name,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.Gnuplot_File_Name,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.Output_Format_Type,
		     INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Output_Format),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Restart_Solution_Save_Frequency),
		     1,
		     MPI::INT,Source_Rank);
  // Output progress frequency:
  Communicator.Bcast(&(IP.Output_Progress_Frequency),
		     1,
		     MPI::INT,Source_Rank);
  // Number of blocks per processor:
  if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
    IP.Number_of_Processors = CFFC_MPI::Number_of_Processors;
  }
  Communicator.Bcast(&(IP.Number_of_Blocks_Per_Processor),
		     1,
		     MPI::INT,Source_Rank);

  // Reinitialize the reference state.
  Reinitialize_Reference_State(IP);

}
#endif

/**********************************************************************
 * Routine: Get_Next_Input_Control_Parameter                          *
 *                                                                    *
 * Get the next input control parameter from the input file.          *
 *                                                                    *
 **********************************************************************/
void Get_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters &IP) {

  int i;
  char buffer[256];

  IP.Line_Number = IP.Line_Number + 1;
  IP.Input_File.getline(buffer,sizeof(buffer));
  i = 0;
  if (buffer[0] != '#') {
    while (1) {
      if (buffer[i] == ' ' || buffer[i] == '=') break;
      i = i + 1;
      if (i > strlen(buffer)) break;
    }
    buffer[i] = '\0';
  }
  strcpy(IP.Next_Control_Parameter,buffer);

}

/**********************************************************************
 * Routine: Parse_Next_Input_Control_Parameter                        *
 *                                                                    *
 * Parses and executes the next input control parameter from the      *
 * input file.                                                        *
 *                                                                    *
 **********************************************************************/
int Parse_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters &IP) {

  int i_command = 0;
  char buffer[256];
  int tpt, bct;
  Vector2D Xt;

  if (strcmp(IP.Next_Control_Parameter, "CFFC_Path") == 0) {
    i_command = 1111;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.CFFC_Path, IP.Next_Control_Parameter);
    strcpy(IP.NASA_Rotor37_Data_Directory, IP.CFFC_Path);
    strcpy(IP.NASA_Rotor37_Data_Directory, "/data/NASA_Rotors/R37/");
    strcpy(IP.NASA_Rotor67_Data_Directory, IP.CFFC_Path);
    strcpy(IP.NASA_Rotor67_Data_Directory, "/data/NASA_Rotors/R67/");

  } else if (strcmp(IP.Next_Control_Parameter,"Time_Integration_Type") == 0) {
    i_command = 1;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Time_Integration_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Time_Integration_Type,"Explicit_Euler") == 0) {
      IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
      IP.N_Stage = 1;
    } else if (strcmp(IP.Time_Integration_Type,"Explicit_Predictor_Corrector") == 0) {
      IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR;
      IP.N_Stage = 2;
    } else if (strcmp(IP.Time_Integration_Type,"Explicit_Runge_Kutta") == 0) {
      IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_RUNGE_KUTTA;
      IP.N_Stage = 4;
    } else if (strcmp(IP.Time_Integration_Type,"Multistage_Optimal_Smoothing") == 0) {
      IP.i_Time_Integration = TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING;
      IP.N_Stage = 3;
    } else if (strcmp(IP.Time_Integration_Type,"Multigrid") == 0) {
      IP.i_Time_Integration = TIME_STEPPING_MULTIGRID;
    } else if (strcmp(IP.Time_Integration_Type,"Dual_Time_Stepping") == 0) {
      IP.i_Time_Integration = TIME_STEPPING_DUAL_TIME_STEPPING;
      IP.Multigrid_IP.i_Dual_Time_Stepping = ON;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Reconstruction_Type") == 0) {
    i_command = 2;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Reconstruction_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Reconstruction_Type,"Green_Gauss") == 0) {
      IP.i_Reconstruction = RECONSTRUCTION_GREEN_GAUSS;
    } else if (strcmp(IP.Reconstruction_Type,"Least_Squares") == 0 ||
	       strcmp(IP.Reconstruction_Type,"Linear_Least_Squares") == 0) {
      IP.i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;
    } else if (strcmp(IP.Reconstruction_Type,"Diamond_Path") == 0) {
      IP.i_Reconstruction = RECONSTRUCTION_DIAMOND_PATH;
    } else if (strcmp(IP.Reconstruction_Type,"Quadratic_Least_Squares") == 0) {
      IP.i_Reconstruction = RECONSTRUCTION_QUADRATIC_LEAST_SQUARES;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Limiter_Type") == 0) {
    i_command = 3;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Limiter_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Limiter_Type,"One") == 0) {
      IP.i_Limiter = LIMITER_ONE;
    } else if (strcmp(IP.Limiter_Type,"Zero") == 0) {
      IP.i_Limiter = LIMITER_ZERO;
    } else if (strcmp(IP.Limiter_Type,"VanLeer") == 0) {
      IP.i_Limiter = LIMITER_VANLEER;
    } else if (strcmp(IP.Limiter_Type,"VanAlbada") == 0) {
      IP.i_Limiter = LIMITER_VANALBADA;
    } else if (strcmp(IP.Limiter_Type,"Barth_Jespersen") == 0) {
      IP.i_Limiter = LIMITER_BARTH_JESPERSEN;
    } else if (strcmp(IP.Limiter_Type,"Venkatakrishnan") == 0) {
      IP.i_Limiter = LIMITER_VENKATAKRISHNAN;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Freeze_Limiter") == 0) {
    i_command = 3;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Freeze_Limiter;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Freeze_Limiter < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Freeze_Limiter_Residual_Level") == 0) {
    i_command = 3;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Freeze_Limiter_Residual_Level;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Freeze_Limiter_Residual_Level < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Flux_Function_Type") == 0) {
    i_command = 4;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Flux_Function_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Flux_Function_Type,"Godunov") == 0) {
      IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV;
    } else if (strcmp(IP.Flux_Function_Type,"Roe") == 0) {
      IP.i_Flux_Function = FLUX_FUNCTION_ROE;
    } else if (strcmp(IP.Flux_Function_Type,"Rusanov") == 0) {
      IP.i_Flux_Function = FLUX_FUNCTION_RUSANOV;
    } else if (strcmp(IP.Flux_Function_Type,"HLLE") == 0) {
      IP.i_Flux_Function = FLUX_FUNCTION_HLLE;
    } else if (strcmp(IP.Flux_Function_Type,"HLLL") == 0 ||
	       strcmp(IP.Flux_Function_Type,"Linde") == 0) {
      IP.i_Flux_Function = FLUX_FUNCTION_HLLL;
    } else if (strcmp(IP.Flux_Function_Type,"HLLC") == 0) {
      IP.i_Flux_Function = FLUX_FUNCTION_HLLC;
    } else if (strcmp(IP.Flux_Function_Type,"VanLeer") == 0) {
      IP.i_Flux_Function = FLUX_FUNCTION_VANLEER;
    } else if (strcmp(IP.Flux_Function_Type,"AUSM") == 0) {
      IP.i_Flux_Function = FLUX_FUNCTION_AUSM;
    } else if (strcmp(IP.Flux_Function_Type,"AUSM+") == 0) {
      IP.i_Flux_Function = FLUX_FUNCTION_AUSMplus;
    } else if (strcmp(IP.Flux_Function_Type,"Godunov_MB") == 0) {
      IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV_MB;
    } else if (strcmp(IP.Flux_Function_Type,"Roe_MB") == 0) {
      IP.i_Flux_Function = FLUX_FUNCTION_ROE_MB;
    } else if (strcmp(IP.Flux_Function_Type,"HLLE_MB") == 0) {
      IP.i_Flux_Function = FLUX_FUNCTION_HLLE_MB;
    } else if (strcmp(IP.Flux_Function_Type,"VanLeer_MB") == 0) {
      IP.i_Flux_Function = FLUX_FUNCTION_VANLEER_MB;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Viscous_Reconstruction_Type") == 0) {
    i_command = 4;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Viscous_Reconstruction_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Viscous_Reconstruction_Type,"Arithmetic_Average") == 0) {
      IP.i_Viscous_Reconstruction = VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE;
    } else if (strcmp(IP.Viscous_Reconstruction_Type,"Diamond_Path") == 0) {
      IP.i_Viscous_Reconstruction = VISCOUS_RECONSTRUCTION_DIAMOND_PATH;
    } else if (strcmp(IP.Viscous_Reconstruction_Type,"Hybrid") == 0) {
      IP.i_Viscous_Reconstruction = VISCOUS_RECONSTRUCTION_HYBRID;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"ICs_Type") == 0) {
    i_command = 5;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.ICs_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.ICs_Type,"Constant") == 0) {
      IP.i_ICs = IC_CONSTANT;
    } else if (strcmp(IP.ICs_Type,"Uniform") == 0) {
      IP.i_ICs = IC_UNIFORM;
    } else if (strcmp(IP.ICs_Type,"Sod") == 0) {
      IP.i_ICs = IC_SOD;
    } else if (strcmp(IP.ICs_Type,"Sod_Xdir") == 0) {
      IP.i_ICs = IC_SOD_XDIR;
    } else if (strcmp(IP.ICs_Type,"Sod_Ydir") == 0) {
      IP.i_ICs = IC_SOD_YDIR;
    } else if (strcmp(IP.ICs_Type,"Groth") == 0) {
      IP.i_ICs = IC_GROTH;
    } else if (strcmp(IP.ICs_Type,"Groth_Xdir") == 0) {
      IP.i_ICs = IC_GROTH_XDIR;
    } else if (strcmp(IP.ICs_Type,"Groth_Ydir") == 0) {
      IP.i_ICs = IC_GROTH_YDIR;
    } else if (strcmp(IP.ICs_Type,"Einfeldt") == 0) {
      IP.i_ICs = IC_EINFELDT;
    } else if (strcmp(IP.ICs_Type,"Einfeldt_Xdir") == 0) {
      IP.i_ICs = IC_EINFELDT_XDIR;
    } else if (strcmp(IP.ICs_Type,"Einfeldt_Ydir") == 0) {
      IP.i_ICs = IC_EINFELDT_YDIR;
    } else if (strcmp(IP.ICs_Type,"Shock_Wave_Xdir") == 0) {
      IP.i_ICs = IC_SHOCK_WAVE_XDIR;
    } else if (strcmp(IP.ICs_Type,"Shock_Wave_Ydir") == 0) {
      IP.i_ICs = IC_SHOCK_WAVE_YDIR;
    } else if (strcmp(IP.ICs_Type,"Contact_Surface_Xdir") == 0) {
      IP.i_ICs = IC_CONTACT_SURFACE_XDIR;
    } else if (strcmp(IP.ICs_Type,"Contact_Surface_Ydir") == 0) {
      IP.i_ICs = IC_CONTACT_SURFACE_YDIR;
    } else if (strcmp(IP.ICs_Type,"Rarefaction_Wave_Xdir") == 0) {
      IP.i_ICs = IC_RAREFACTION_WAVE_XDIR;
    } else if (strcmp(IP.ICs_Type,"Rarefaction_Wave_Ydir") == 0) {
      IP.i_ICs = IC_RAREFACTION_WAVE_YDIR;
    } else if (strcmp(IP.ICs_Type,"ShockBox") == 0) {
      IP.i_ICs = IC_SHOCK_BOX;
    } else if (strcmp(IP.ICs_Type,"High_Pressure_Reservoir") == 0) {
      IP.i_ICs = IC_HIGH_PRESSURE_RESERVOIR;
    } else if (strcmp(IP.ICs_Type,"Low_Pressure_Reservoir") == 0) {
      IP.i_ICs = IC_LOW_PRESSURE_RESERVOIR;
    } else if (strcmp(IP.ICs_Type,"Riemann") == 0) {
      IP.i_ICs = IC_RIEMANN;
    } else if (strcmp(IP.ICs_Type,"Riemann_Xdir") == 0) {
      IP.i_ICs = IC_RIEMANN_XDIR;
    } else if (strcmp(IP.ICs_Type,"Riemann_Ydir") == 0) {
      IP.i_ICs = IC_RIEMANN_YDIR;
    } else if (strcmp(IP.ICs_Type,"Square_Wave_Xdir") == 0) {
      IP.i_ICs = IC_SQUARE_WAVE_XDIR;
    } else if (strcmp(IP.ICs_Type,"Square_Wave_Ydir") == 0) {
      IP.i_ICs = IC_SQUARE_WAVE_YDIR;
    } else if (strcmp(IP.ICs_Type,"Sine_Wave_Xdir") == 0) {
      IP.i_ICs = IC_SINE_WAVE_XDIR;
    } else if (strcmp(IP.ICs_Type,"Sine_Wave_Ydir") == 0) {
      IP.i_ICs = IC_SINE_WAVE_YDIR;
    } else if (strcmp(IP.ICs_Type,"Viscous_Channel_Flow") == 0) {
      IP.i_ICs = IC_VISCOUS_CHANNEL_FLOW;
    } else if (strcmp(IP.ICs_Type,"Viscous_Pipe_Flow") == 0 ||
	       strcmp(IP.ICs_Type,"Turbulent_Pipe_Flow") == 0) {
      IP.i_ICs = IC_VISCOUS_PIPE_FLOW;
    } else if (strcmp(IP.ICs_Type,"Flat_Plate") == 0) {
      IP.i_ICs = IC_VISCOUS_FLAT_PLATE;
      IP.BC_South = BC_WALL_VISCOUS_HEATFLUX;
    } else if (strcmp(IP.ICs_Type,"Stokes_Flow") == 0) {
      IP.i_ICs = IC_VISCOUS_STOKES_FLOW;
    } else if (strcmp(IP.ICs_Type,"Driven_Cavity_Flow") == 0) {
      IP.i_ICs = IC_VISCOUS_DRIVEN_CAVITY_FLOW;
    } else if (strcmp(IP.ICs_Type,"Backward_Facing_Step") == 0) {
      IP.i_ICs = IC_VISCOUS_BACKWARD_FACING_STEP;
    } else if (strcmp(IP.ICs_Type,"Jet_Flow") == 0) {
      IP.i_ICs = IC_JET_FLOW; 
    } else if (strcmp(IP.ICs_Type,"Mixing_Layer") == 0) {
      IP.i_ICs = IC_MIXING_LAYER;  
    } else if (strcmp(IP.ICs_Type,"Branched_Duct") == 0) {
      IP.i_ICs = IC_VISCOUS_BRANCHED_DUCT;
    } else if (strcmp(IP.ICs_Type,"Cylindrical_Explosion") == 0) {
      IP.i_ICs = IC_CYLINDRICAL_EXPLOSION;
    } else if (strcmp(IP.ICs_Type,"Cylindrical_Implosion") == 0) {
      IP.i_ICs = IC_CYLINDRICAL_IMPLOSION;
    } else if (strcmp(IP.ICs_Type,"Ringleb_Flow") == 0) {
      IP.i_ICs = IC_RINGLEB_FLOW;
    } else if (strcmp(IP.ICs_Type,"Electrostatic_Channel") == 0) {
      IP.i_ICs = IC_ELECTROSTATIC_CHANNEL;
    } else if (strcmp(IP.ICs_Type,"Restart") == 0) {
      IP.i_ICs = IC_RESTART;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Grid_Type") == 0) {
    i_command = 6;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Grid_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Grid_Type,"Cartesian") == 0) {
      IP.i_Grid = GRID_CARTESIAN_UNIFORM;
      IP.Box_Width = ONE;
      IP.Box_Height = ONE;
    } else if (strcmp(IP.Grid_Type,"Square") == 0) {
      IP.i_Grid = GRID_SQUARE;
      IP.Box_Width = ONE;
      IP.Box_Height = ONE;
    } else if (strcmp(IP.Grid_Type,"Rectangular_Box") == 0) {
      IP.i_Grid = GRID_RECTANGULAR_BOX;
      IP.Box_Width = ONE;
      IP.Box_Height = ONE;
    } else if (strcmp(IP.Grid_Type,"Flat_Plate") == 0 ||
	       strcmp(IP.Grid_Type,"Adiabatic_Flat_Plate") == 0) {
      IP.i_Grid = GRID_FLAT_PLATE;
      IP.Plate_Length = ONE;
      IP.BC_South = BC_WALL_VISCOUS_HEATFLUX;
    } else if (strcmp(IP.Grid_Type,"Isothermal_Flat_Plate") == 0) {
      IP.i_Grid = GRID_FLAT_PLATE;
      IP.Plate_Length = ONE;
      IP.BC_South = BC_WALL_VISCOUS_ISOTHERMAL;
    } else if (strcmp(IP.Grid_Type,"Burning_Surface_Flat_Plate") == 0) {
      IP.i_Grid = GRID_FLAT_PLATE;
      IP.Plate_Length = ONE;
      IP.BC_South = BC_BURNING_SURFACE;
    } else if (strcmp(IP.Grid_Type,"Pipe") == 0) {
      IP.i_Grid = GRID_PIPE;
      IP.Pipe_Length = ONE;
      IP.Pipe_Radius = HALF;
    } else if (strcmp(IP.Grid_Type,"Blunt_Body") == 0) {
      IP.i_Grid = GRID_BLUNT_BODY;
      IP.Blunt_Body_Radius = ONE;
      IP.Blunt_Body_Mach_Number = TWO;
    } else if (strcmp(IP.Grid_Type,"Rocket_Motor") == 0) {
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
      IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_LINEAR;
      IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MAX_CLUSTERING;
      IP.Mesh_Stretching_Factor_Idir = 1.10;
      IP.Mesh_Stretching_Factor_Jdir = 1.01;
    } else if (strcmp(IP.Grid_Type,"Nozzleless_Rocket_Motor") == 0) {
      IP.i_Grid = GRID_NOZZLELESS_ROCKET_MOTOR;
      IP.Chamber_Length = 0.48;
      IP.Chamber_Radius = 0.01;
      IP.Nozzle_Length = 0.032;
      IP.Wedge_Angle = 15.0*PI/180.0;
      IP.Nozzle_Radius_Exit = IP.Chamber_Radius + IP.Nozzle_Length*tan(IP.Wedge_Angle);
      IP.BC_North = BC_WALL_VISCOUS_ISOTHERMAL;
      IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_LINEAR;
      IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MAX_CLUSTERING;
      IP.Mesh_Stretching_Factor_Idir = 1.99;
      IP.Mesh_Stretching_Factor_Jdir = 1.01;
    } else if (strcmp(IP.Grid_Type,"Circular_Cylinder") == 0) {
      IP.i_Grid = GRID_CIRCULAR_CYLINDER;
      IP.Cylinder_Radius = ONE;
      IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
      IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
      IP.Mesh_Stretching_Factor_Idir = 1.025;
      IP.Mesh_Stretching_Factor_Jdir = 1.001;
    } else if (strcmp(IP.Grid_Type,"Ellipse") == 0) {
      IP.i_Grid = GRID_ELLIPSE;
      IP.Ellipse_Length_X_Axis = TWO;
      IP.Ellipse_Length_Y_Axis = HALF;
    } else if (strcmp(IP.Grid_Type,"NACA_Aerofoil") == 0) {
      IP.i_Grid = GRID_NACA_AEROFOIL;
      IP.Chord_Length = ONE;
      strcpy(IP.NACA_Aerofoil_Type,"0012");
    } else if (strcmp(IP.Grid_Type,"Free_Jet") == 0) {
      IP.i_Grid = GRID_FREE_JET;
      IP.Orifice_Radius = ONE;
    } else if (strcmp(IP.Grid_Type,"Mixing_Layer") ==0) {
      IP.i_Grid = GRID_MIXING_LAYER;
      IP.Box_Width = 2.0;
      IP.Mesh_Stretching_Factor_Idir = 1.0001;
      IP.Mesh_Stretching_Factor_Jdir = 1.0010;
    } else if (strcmp(IP.Grid_Type,"Jet_Flow") == 0) {
      IP.i_Grid = GRID_JET_FLOW;
    } else if (strcmp(IP.Grid_Type,"Ringleb_Flow") == 0) {
      IP.i_Grid = GRID_RINGLEB_FLOW;
      IP.Inner_Streamline_Number = 0.80;
      IP.Outer_Streamline_Number = 0.40;
      IP.Isotach_Line = 0.30;
    } else if (strcmp(IP.Grid_Type,"Wedge") == 0) {
      IP.i_Grid = GRID_WEDGE;
      IP.Wedge_Angle = 25.0;
      IP.Wedge_Length = HALF;
      IP.BC_South = BC_WALL_VISCOUS_ISOTHERMAL;
    } else if (strcmp(IP.Grid_Type,"NASA_Rotor_37") == 0) {
      IP.i_Grid = GRID_NASA_ROTOR_37;
      IP.Rotor_Flow_Type = PEAK_FLOW;
      IP.Rotor_Percent_Span = 50.00;
    } else if (strcmp(IP.Grid_Type,"NASA_Rotor_67") == 0) {
      IP.i_Grid = GRID_NASA_ROTOR_67;
      IP.Rotor_Flow_Type = PEAK_FLOW;
      IP.Rotor_Percent_Span = 50.00;
    } else if (strcmp(IP.Grid_Type,"Ringleb_Flow") == 0) {
      IP.i_Grid = GRID_RINGLEB_FLOW;
    } else if (strcmp(IP.Grid_Type,"Unsteady_Blunt_Body") == 0) {
      IP.i_Grid = GRID_UNSTEADY_BLUNT_BODY;
      IP.Blunt_Body_Radius = ONE;
      IP.Blunt_Body_Mach_Number = TWO;
    } else if (strcmp(IP.Grid_Type,"Bump_Channel_Flow") == 0) {
      IP.i_Grid = GRID_BUMP_CHANNEL_FLOW;
    } else if (strcmp(IP.Grid_Type,"Non_Smooth_Bump_Channel_Flow") == 0) {
      IP.i_Grid = GRID_BUMP_CHANNEL_FLOW;
      IP.Smooth_Bump = OFF;
    } else if (strcmp(IP.Grid_Type,"Smooth_Bump_Channel_Flow") == 0) {
      IP.i_Grid = GRID_BUMP_CHANNEL_FLOW;
      IP.Smooth_Bump = ON;
    } else if (strcmp(IP.Grid_Type,"Driven_Cavity_Flow") == 0) {
      IP.i_Grid = GRID_DRIVEN_CAVITY_FLOW;
    } else if (strcmp(IP.Grid_Type,"Backward_Facing_Step") == 0) {
      IP.i_Grid = GRID_BACKWARD_FACING_STEP;
    } else if (strcmp(IP.Grid_Type,"Nozzle") == 0) {
      IP.i_Grid = GRID_NOZZLE;
    } else if (strcmp(IP.Grid_Type,"Desolvation_Chamber") == 0) {
      IP.i_Grid = GRID_DESOLVATION_CHAMBER;
    } else if (strcmp(IP.Grid_Type,"ICEMCFD") == 0) {
      IP.i_Grid = GRID_ICEMCFD;
      strcpy(IP.ICEMCFD_FileNames[0],"topo_mulcad_out");
      strcpy(IP.ICEMCFD_FileNames[1],"family_boco");
      strcpy(IP.ICEMCFD_FileNames[2],"family_topo");
    } else if (strcmp(IP.Grid_Type,"Read_From_Definition_File") == 0) {
      IP.i_Grid = GRID_READ_FROM_DEFINITION_FILE;
    } else if (strcmp(IP.Grid_Type,"Read_From_Data_File") == 0) {
      IP.i_Grid = GRID_READ_FROM_GRID_DATA_FILE;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Output_File_Name") == 0) {
    i_command = 7;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Output_File_Name,IP.Next_Control_Parameter);
    strcat(IP.Output_File_Name,".dat");
    strcpy(IP.Grid_File_Name,IP.Next_Control_Parameter);
    strcat(IP.Grid_File_Name,".grid");
    strcpy(IP.Grid_Definition_File_Name,IP.Next_Control_Parameter);
    strcat(IP.Grid_Definition_File_Name,".griddef");
    strcpy(IP.Restart_File_Name,IP.Next_Control_Parameter);
    strcat(IP.Restart_File_Name,".soln");
    strcpy(IP.Gnuplot_File_Name,IP.Next_Control_Parameter);
    strcat(IP.Gnuplot_File_Name,".gplt");

  } else if (strcmp(IP.Next_Control_Parameter,"Grid_File_Name") == 0) {
    i_command = 8;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Grid_File_Name,IP.Next_Control_Parameter);
    strcat(IP.Grid_File_Name,".grid");
    strcpy(IP.Grid_Definition_File_Name,IP.Next_Control_Parameter);
    strcat(IP.Grid_Definition_File_Name,".griddef");

  } else if (strcmp(IP.Next_Control_Parameter,"Restart_File_Name") == 0) {
    i_command = 9;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Restart_File_Name,IP.Next_Control_Parameter);
    strcat(IP.Restart_File_Name,".soln");

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Cells_Idir") == 0) {
    i_command = 10;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Cells_Idir;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Cells_Idir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Cells_Jdir") == 0) {
    i_command = 11;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Cells_Jdir;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Cells_Jdir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Ghost_Cells") == 0) {
    i_command = 12;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Ghost_Cells;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Ghost_Cells < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Blocks_Idir") == 0) {
    i_command = 12;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Blocks_Idir;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Blocks_Idir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Blocks_Jdir") == 0) {
    i_command = 13;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Blocks_Jdir;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Blocks_Jdir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Time_Accurate") == 0) {
    i_command = 14;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Time_Accurate;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Time_Accurate != 0 && IP.Time_Accurate != 1) IP.Time_Accurate = 0;
    if (IP.Time_Accurate) IP.Local_Time_Stepping = GLOBAL_TIME_STEPPING;

  } else if (strcmp(IP.Next_Control_Parameter,"Local_Time_Stepping") == 0) {
    i_command = 15;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Local_Time_Stepping;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Local_Time_Stepping != GLOBAL_TIME_STEPPING &&
	IP.Local_Time_Stepping != SCALAR_LOCAL_TIME_STEPPING &&
	IP.Local_Time_Stepping != SEMI_IMPLICIT_LOCAL_TIME_STEPPING) {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Maximum_Number_of_Time_Steps") == 0) {
    i_command = 16;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Maximum_Number_of_Time_Steps;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Maximum_Number_of_Time_Steps < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"N_Stage") == 0) {
    i_command = 17;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.N_Stage;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.i_Time_Integration == TIME_STEPPING_EXPLICIT_EULER && IP.N_Stage != 1) i_command = INVALID_INPUT_VALUE;
    if (IP.i_Time_Integration == TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR && IP.N_Stage != 2) i_command = INVALID_INPUT_VALUE;
    if (IP.i_Time_Integration == TIME_STEPPING_EXPLICIT_RUNGE_KUTTA && IP.N_Stage != 4) i_command = INVALID_INPUT_VALUE;
    if (IP.N_Stage < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"CFL_Number") == 0) {
    i_command = 18;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.CFL_Number;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.CFL_Number <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Box_Width") == 0) {
    i_command = 19;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Box_Width;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Box_Width <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Box_Height") == 0) {
    i_command = 20;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Box_Height;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Box_Height <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Plate_Length") == 0) {
    i_command = 21;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Plate_Length;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Plate_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Pipe_Length") == 0) {
    i_command = 22;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Pipe_Length;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Pipe_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Pipe_Radius") == 0) {
    i_command = 23;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Pipe_Radius;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Pipe_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Blunt_Body_Radius") == 0) {
    i_command = 24;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Blunt_Body_Radius;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Blunt_Body_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Blunt_Body_Mach_Number") == 0) {
    i_command = 25;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Blunt_Body_Mach_Number;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Blunt_Body_Mach_Number <= ONE) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Cylinder_Radius") == 0) {
    i_command = 26;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Cylinder_Radius;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Cylinder_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Ellipse_Length_X_Axis") == 0) {
    i_command = 27;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Ellipse_Length_X_Axis;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Ellipse_Length_X_Axis <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Ellipse_Length_Y_Axis") == 0) {
    i_command = 28;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Ellipse_Length_Y_Axis;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Ellipse_Length_Y_Axis <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Chord_Length") == 0) {
    i_command = 29;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Chord_Length;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Chord_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"NACA_Aerofoil_Type") == 0) {
    i_command = 30;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.NACA_Aerofoil_Type,IP.Next_Control_Parameter);
    if (strlen(IP.NACA_Aerofoil_Type) != 4 &&
	strlen(IP.NACA_Aerofoil_Type) != 5) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Orifice_Radius") == 0) {
    i_command = 31;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Orifice_Radius;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Orifice_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Inner_Streamline_Number") == 0) {
    i_command = 32;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Inner_Streamline_Number;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Inner_Streamline_Number <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Outer_Streamline_Number") == 0) {
    i_command = 33;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Outer_Streamline_Number;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Outer_Streamline_Number <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Isotach_Line") == 0) {
    i_command = 34;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Isotach_Line;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Isotach_Line <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Wedge_Angle") == 0) {
    i_command = 35;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Wedge_Angle;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Wedge_Angle <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Wedge_Length") == 0) {
    i_command = 36;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Wedge_Length;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Wedge_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Smooth_Bump") == 0) {
    i_command = 37;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
      IP.Smooth_Bump = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
      IP.Smooth_Bump = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Step_Height") == 0) {
    i_command = 38;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Step_Height;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Step_Height <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Top_Wall_Deflection") == 0) {
    i_command = 39;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Top_Wall_Deflection;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Top_Wall_Deflection < -TWO || IP.Top_Wall_Deflection > TEN) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Chamber_Length") == 0) {
    i_command = 40;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Chamber_Length;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Chamber_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Chamber_Radius") == 0) {
    i_command = 41;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Chamber_Radius;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Chamber_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Chamber_To_Throat_Length") == 0) {
    i_command = 42;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Chamber_To_Throat_Length;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Chamber_To_Throat_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Nozzle_Length") == 0) {
    i_command = 43;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Nozzle_Length;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Nozzle_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Nozzle_Radius_Exit") == 0) {
    i_command = 44;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Nozzle_Radius_Exit;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Nozzle_Radius_Exit <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Nozzle_Radius_Throat") == 0) {
    i_command = 45;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Nozzle_Radius_Throat;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Nozzle_Radius_Throat <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Nozzle_Type") == 0) {
    i_command = 46;
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
    i_command = 47;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Grain_Radius;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"Gas_Type") == 0) {
    i_command = 38;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Gas_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Gas_Type,"AIR") != 0 &&
	strcmp(IP.Gas_Type,"H2") != 0 &&
	strcmp(IP.Gas_Type,"HE") != 0 &&
	strcmp(IP.Gas_Type,"N2") != 0 &&
	strcmp(IP.Gas_Type,"O2") != 0 &&
	strcmp(IP.Gas_Type,"AIR-viscous") != 0 &&
	strcmp(IP.Gas_Type,"AP_HTPB") != 0) {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Mach_Number") == 0) {
    i_command = 39;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Mach_Number;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Mach_Number < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Mach_Number2") == 0) {
    i_command = 39;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Mach_Number2;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Mach_Number2 < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Flow_Angle") == 0) {
    i_command = 40;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Flow_Angle;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"Pressure") == 0) {
    i_command = 41;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Pressure;
    IP.Input_File.getline(buffer,sizeof(buffer));
    IP.Pressure = IP.Pressure*THOUSAND;
    if (IP.Pressure <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Temperature") == 0) {
    i_command = 42;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Temperature;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Temperature <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Reynolds_Number") == 0) {
    i_command = 43;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Reynolds_Number;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Reynolds_Number < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Pressure_Gradient") == 0) {
    i_command = 44;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.dp;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"Re_lid") == 0) {
    i_command = 44;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Re_lid;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Re_lid < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Wave_Position") == 0) {
    i_command = 44;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Wave_Position;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"Wave_Width") == 0) {
    i_command = 44;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Wave_Width;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Wave_Width < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Time_Max") == 0) {
    i_command = 45;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Time_Max;
    IP.Input_File.getline(buffer,sizeof(buffer));
    IP.Time_Max = IP.Time_Max/THOUSAND;
    if (IP.Time_Max < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Blocks_Per_Processor") == 0) {
    i_command = 46;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Blocks_Per_Processor;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Blocks_Per_Processor < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Output_Format_Type") == 0) {
    i_command = 47;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Output_Format_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Output_Format_Type,"Gnuplot") == 0  ||
	strcmp(IP.Output_Format_Type,"gnuplot") == 0  ||
	strcmp(IP.Output_Format_Type,"GNUPLOT") == 0) {
      IP.i_Output_Format = IO_GNUPLOT;
    } else if (strcmp(IP.Output_Format_Type,"Tecplot") == 0  ||
	       strcmp(IP.Output_Format_Type,"tecplot") == 0  ||
	       strcmp(IP.Output_Format_Type,"TECPLOT") == 0) {
      IP.i_Output_Format = IO_TECPLOT;
    } else if (strcmp(IP.Output_Format_Type,"Matlab") == 0  ||
	       strcmp(IP.Output_Format_Type,"matlab") == 0  ||
	       strcmp(IP.Output_Format_Type,"MATLAB") == 0) {
      IP.i_Output_Format = IO_MATLAB;
    } else if (strcmp(IP.Output_Format_Type,"Octave") == 0  ||
	       strcmp(IP.Output_Format_Type,"octave") == 0  ||
	       strcmp(IP.Output_Format_Type,"OCTAVE") == 0) {
      IP.i_Output_Format = IO_OCTAVE;
    } else {
      IP.i_Output_Format = IO_TECPLOT;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Flow_Type") == 0) {
    i_command = 48;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Flow_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Flow_Type,"Inviscid") == 0) {
      IP.FlowType = FLOWTYPE_INVISCID;
    } else if (strcmp(IP.Flow_Type,"Laminar") == 0) {
      IP.FlowType = FLOWTYPE_LAMINAR;
    } else if (strcmp(IP.Flow_Type,"Turbulent") == 0 ||
	       strcmp(IP.Flow_Type,"Turbulent_k_omega") == 0) {
      IP.FlowType = FLOWTYPE_TURBULENT_RANS_K_OMEGA;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Flow_Geometry_Type") == 0) {
    i_command = 49;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Flow_Geometry_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Flow_Geometry_Type,"Planar") == 0) {
      IP.Axisymmetric = OFF;
    } else if (strcmp(IP.Flow_Geometry_Type,"Axisymmetric") == 0) {
      IP.Axisymmetric = ON;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Wall_Speed") == 0) {
    i_command = 37;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Vwall.x >> IP.Vwall.y;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"Wall_Temperature") == 0) {
    i_command = 37;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Twall;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"NASA_Rotor37_Data_Directory") == 0) {
    i_command = 46;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.NASA_Rotor37_Data_Directory,IP.Next_Control_Parameter);

  } else if (strcmp(IP.Next_Control_Parameter,"NASA_Rotor67_Data_Directory") == 0) {
    i_command = 47;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.NASA_Rotor67_Data_Directory,IP.Next_Control_Parameter);

  } else if (strcmp(IP.Next_Control_Parameter,"Rotor_Flow_Type") == 0) {
    i_command = 48;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Rotor_Flow_Type;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Rotor_Flow_Type != PEAK_FLOW && IP.Rotor_Flow_Type != STALL_FLOW)
      IP.Rotor_Flow_Type = PEAK_FLOW;

  } else if (strcmp(IP.Next_Control_Parameter,"Rotor_Percent_Span") == 0) {
    i_command = 49;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Rotor_Percent_Span;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Rotor_Percent_Span < ZERO || IP.Rotor_Percent_Span > HUNDRED)
      i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Restart_Solution_Save_Frequency") == 0) {
    i_command = 50;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Restart_Solution_Save_Frequency;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Restart_Solution_Save_Frequency < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Output_Progress_Frequency") == 0) {
    i_command = 51;
    IP.Line_Number++;
    IP.Input_File >> IP.Output_Progress_Frequency;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Output_Progress_Frequency < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"ICEMCFD_Topology_File") == 0) {
    i_command = 53;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.ICEMCFD_FileNames[0],IP.Next_Control_Parameter);

  } else if (strcmp(IP.Next_Control_Parameter,"ICEMCFD_Family_Boco_File") == 0) {
    i_command = 54;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.ICEMCFD_FileNames[1],IP.Next_Control_Parameter);

  } else if (strcmp(IP.Next_Control_Parameter,"ICEMCFD_Family_Topo_File") == 0) {
    i_command = 55;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.ICEMCFD_FileNames[2],IP.Next_Control_Parameter);

  } else if (strcmp(IP.Next_Control_Parameter,"X_Shift") == 0) {
    i_command = 56;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.X_Shift;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"X_Scale") == 0) {
    i_command = 57;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.X_Scale;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"X_Rotate") == 0) {
    i_command = 58;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.X_Rotate;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"AMR") == 0) {
    i_command = 71;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
      IP.AMR = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
      IP.AMR = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"AMR_Frequency") == 0) {
    i_command = 72;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.AMR_Frequency;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.AMR_Frequency < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Initial_Mesh_Refinements") == 0) {
    i_command = 73;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Initial_Mesh_Refinements;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Initial_Mesh_Refinements < 0) IP.Number_of_Initial_Mesh_Refinements = 0;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Uniform_Mesh_Refinements") == 0) {
    i_command = 74;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Uniform_Mesh_Refinements;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Uniform_Mesh_Refinements < 0) IP.Number_of_Uniform_Mesh_Refinements = 0;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Boundary_Mesh_Refinements") == 0) {
    i_command = 75;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Boundary_Mesh_Refinements;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Boundary_Mesh_Refinements < 0) IP.Number_of_Boundary_Mesh_Refinements = 0;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Interface_Mesh_Refinements") == 0) {
    i_command = 76;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Interface_Mesh_Refinements;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Interface_Mesh_Refinements < 0) IP.Number_of_Interface_Mesh_Refinements = 0;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Bounding_Box_Mesh_Refinements") == 0) {
    i_command = 77;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Bounding_Box_Mesh_Refinements;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Bounding_Box_Mesh_Refinements < 0) IP.Number_of_Bounding_Box_Mesh_Refinements = 0;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Flat_Plate_Mesh_Refinements") == 0) {
    i_command = 78;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Flat_Plate_Mesh_Refinements;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Flat_Plate_Mesh_Refinements < 0) IP.Number_of_Flat_Plate_Mesh_Refinements = 0;

  } else if (strcmp(IP.Next_Control_Parameter,"Interface_Refinement_Condition") == 0) {
    i_command = 79;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
      IP.Interface_Refinement_Condition = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
      IP.Interface_Refinement_Condition = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Maximum_Refinement_Level") == 0) {
    i_command = 80;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Maximum_Refinement_Level;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Maximum_Refinement_Level < 1) IP.Maximum_Refinement_Level = 1;

  } else if (strcmp(IP.Next_Control_Parameter,"Minimum_Refinement_Level") == 0) {
    i_command = 81;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Minimum_Refinement_Level;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Minimum_Refinement_Level < 1) IP.Minimum_Refinement_Level = 1;

  } else if (strcmp(IP.Next_Control_Parameter,"Threshold_for_Refinement") == 0) {
    i_command = 82;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Threshold_for_Refinement;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Threshold_for_Refinement <= ZERO ||
	IP.Threshold_for_Refinement > ONE) IP.Threshold_for_Refinement = 0.50;

  } else if (strcmp(IP.Next_Control_Parameter,"Threshold_for_Coarsening") == 0) {
    i_command = 83;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Threshold_for_Coarsening;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Threshold_for_Coarsening < ZERO ||
	IP.Threshold_for_Coarsening >= ONE) IP.Threshold_for_Coarsening = 0.10;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Refinement_Criteria") == 0) {
    i_command = 84;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Refinement_Criteria;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Refinement_Criteria < 1 || IP.Number_of_Refinement_Criteria > 6) {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Gradient_Density") == 0) {
    i_command = 85;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
      IP.Refinement_Criteria_Gradient_Density = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
      IP.Refinement_Criteria_Gradient_Density = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Divergence_Velocity") == 0) {
    i_command = 86;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
      IP.Refinement_Criteria_Divergence_Velocity = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
      IP.Refinement_Criteria_Divergence_Velocity = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Curl_Velocity") == 0) {
    i_command = 87;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
      IP.Refinement_Criteria_Curl_Velocity = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
      IP.Refinement_Criteria_Curl_Velocity = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Gradient_Turbulence_Kinetic_Energy") == 0) {
    i_command = 88;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
      IP.Refinement_Criteria_Gradient_Turbulence_Kinetic_Energy = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
      IP.Refinement_Criteria_Gradient_Turbulence_Kinetic_Energy = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"AMR_Xmin") == 0) {
    i_command = 89;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.AMR_Xmin;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"AMR_Xmax") == 0) {
    i_command = 90;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.AMR_Xmax;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Morton") == 0) {
    i_command = 90;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
      IP.Morton = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
      IP.Morton = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter, "Morton_Reordering_Frequency") == 0) {
    i_command = 90;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Morton_Reordering_Frequency;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"Residual_Variable") == 0) {
    i_command = 91;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.i_Residual_Variable;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.i_Residual_Variable < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Residual_Smoothing_Epsilon") == 0) {
    i_command = 92;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Residual_Smoothing_Epsilon;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Residual_Smoothing_Epsilon <= ZERO) {
      IP.Residual_Smoothing = OFF;
      IP.Residual_Smoothing_Epsilon = ZERO;
    } else {
      IP.Residual_Smoothing = ON;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Residual_Smoothing_Gauss_Seidel_Iterations") == 0) {
    i_command = 93;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Residual_Smoothing_Gauss_Seidel_Iterations;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Residual_Smoothing_Gauss_Seidel_Iterations < 0) i_command = INVALID_INPUT_VALUE;
    if (IP.Time_Accurate && IP.i_Time_Integration != TIME_STEPPING_DUAL_TIME_STEPPING) i_command = INVALID_INPUT_VALUE;

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
    // TURBULENCE PARAMETERS                                          //
    ////////////////////////////////////////////////////////////////////

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

  } else if (strcmp(IP.Next_Control_Parameter,"Compressibility_Effect") == 0) {
    i_command = 70;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Compressibility_Effect,IP.Next_Control_Parameter);
    if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
      IP.i_Compressibility_Effect = OFF;
    } else if (strcmp(IP.Compressibility_Effect,"Sarkar_Model") == 0) {
      IP.i_Compressibility_Effect = COMPRESSIBILITY_CORRECTION_SARKAR;
    } else if (strcmp(IP.Compressibility_Effect,"Zeman_Model") == 0) {
      IP.i_Compressibility_Effect = COMPRESSIBILITY_CORRECTION_ZEMAN;
    } else if (strcmp(IP.Compressibility_Effect,"Wilcox_Model") == 0) {
      IP.i_Compressibility_Effect = COMPRESSIBILITY_CORRECTION_WILCOX;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Variable_Prandtl") == 0) {
    i_command = 68;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Variable_Prandtl,IP.Next_Control_Parameter);
    if (strcmp(IP.Variable_Prandtl,"OFF") == 0) {
      IP.i_Variable_Prandtl = OFF;
    } else if (strcmp(IP.Variable_Prandtl,"ON") == 0) {
      IP.i_Variable_Prandtl = ON;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Transition_Model") == 0) {
    i_command = 69;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Transition_Model,IP.Next_Control_Parameter);
    if (strcmp(IP.Transition_Model,"OFF") == 0) {
      IP.i_Transition_Model = OFF;
    } else if (strcmp(IP.Transition_Model,"Wilcox") == 0) {
      IP.i_Transition_Model = TRANSITION_WILCOX;
    } else if (strcmp(IP.Transition_Model,"Menter") == 0) {
      IP.i_Transition_Model = TRANSITION_MENTER;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"C_lambda") == 0) {
    i_command = 70;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.C_lambda;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.C_lambda < ZERO) i_command = INVALID_INPUT_VALUE;
    
  } else if (strcmp(IP.Next_Control_Parameter,"Cd1") == 0) {
    i_command = 71;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Cd1;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Cd1 < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Cd4") == 0) {
    i_command = 72;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Cd4;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Cd4 < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Cd5") == 0) {
    i_command = 73;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Cd5;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Cd5 < ZERO) i_command = INVALID_INPUT_VALUE;
    
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

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Cells_Sublayer_Support") == 0) {
    i_command = 168;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.n_cells_sublayer_support;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.n_cells_sublayer_support < ZERO) i_command = INVALID_INPUT_VALUE;
    if (IP.n_cells_sublayer_support > IP.Number_of_Cells_Idir) i_command = INVALID_INPUT_VALUE;
    if (IP.n_cells_sublayer_support > IP.Number_of_Cells_Jdir) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Turbulent_Wall_Injection") == 0) {
    i_command = 169;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
      IP.i_Turbulent_Wall_Injection = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
      IP.i_Turbulent_Wall_Injection = OFF;
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
    // SOLID PROPELLANT TYPE                                          //
    ////////////////////////////////////////////////////////////////////

  } else if (strcmp(IP.Next_Control_Parameter,"Propellant_Type") == 0) {
    i_command = 450;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Propellant_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Propellant_Type,"AP_HTPB") == 0) {
    } else if (strcmp(IP.Propellant_Type,"QUICK_AP_HTPB") == 0) {
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Surface_Roughness") == 0) {
    i_command = 451;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.sigmav;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.sigmav < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Turbulence_Length_Scale") == 0) {
    i_command = 452;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.lw;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.lw < ZERO) i_command = INVALID_INPUT_VALUE;

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
    } else if (strcmp(IP.BC_North_Type,"Mass_Injection") == 0) {
      IP.BC_North = BC_MASS_INJECTION;
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
    } else if (strcmp(IP.BC_South_Type,"Mass_Injection") == 0) {
      IP.BC_South = BC_MASS_INJECTION;
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
    } else if (strcmp(IP.BC_East_Type,"Mass_Injection") == 0) {
      IP.BC_East = BC_MASS_INJECTION;
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
    } else if (strcmp(IP.BC_West_Type,"Mass_Injection") == 0) {
      IP.BC_West = BC_MASS_INJECTION;
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

    ////////////////////////////////////////////////////////////////////
    // INTERFACE PARAMETERS                                           //
    ////////////////////////////////////////////////////////////////////

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Interface_Components") == 0) {
    i_command = 600;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Interface_IP.Number_of_Components;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Interface_IP.Number_of_Components < 0) {
      i_command = INVALID_INPUT_VALUE;
    } else if (IP.Interface_IP.Number_of_Components > 0) {
      IP.Interface_IP.Component_List.allocate(IP.Interface_IP.Number_of_Components);
      IP.Interface_IP.ns = new int[IP.Interface_IP.Number_of_Components];
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Interface_Spline_Type") == 0) {
    i_command = 601;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Interface_IP.Type,IP.Next_Control_Parameter);
    IP.Interface_IP.ci++; // increment the interface counter.
    if (strcmp(IP.Interface_IP.Type,"Circle") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_CIRCLE;
    } else if (strcmp(IP.Interface_IP.Type,"Ellipse") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_ELLIPSE;
    } else if (strcmp(IP.Interface_IP.Type,"Square") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_SQUARE;
    } else if (strcmp(IP.Interface_IP.Type,"Rectangle") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_RECTANGLE;
    } else if (strcmp(IP.Interface_IP.Type,"Rocket_Propellant_Grain") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_ROCKET_PROPELLANT_GRAIN;
    } else if (strcmp(IP.Interface_IP.Type,"NACA0012_Aerofoil") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_NACA0012_AEROFOIL;
    } else if (strcmp(IP.Interface_IP.Type,"NACA0015_Aerofoil") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_NACA0015_AEROFOIL;
    } else if (strcmp(IP.Interface_IP.Type,"NASA_Rotor_37") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_NASA_ROTOR_37;
    } else if (strcmp(IP.Interface_IP.Type,"NASA_Rotor_67") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_NASA_ROTOR_67;
    } else if (strcmp(IP.Interface_IP.Type,"Zalesaks_Disk") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_ZALESAK;
    } else if (strcmp(IP.Interface_IP.Type,"Ringleb_Flow") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_RINGLEB;
    } else if (strcmp(IP.Interface_IP.Type,"Bump_Channel_Flow") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_BUMP;
    } else if (strcmp(IP.Interface_IP.Type,"Flat_Plate") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_FLAT_PLATE;
    } else if (strcmp(IP.Interface_IP.Type,"User_Specified") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_USER_SPECIFIED;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }
    IP.Interface_IP.Component_List[IP.Interface_IP.ci].allocate(2);

  } else if (strcmp(IP.Next_Control_Parameter,"Interface_Characteristic_Length") == 0 ||
	     strcmp(IP.Next_Control_Parameter,"Interface_Characteristic_Length_1") == 0) {
    i_command = 602;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Interface_IP.Component_List[IP.Interface_IP.ci].Length1;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Interface_IP.Component_List[IP.Interface_IP.ci].Length1 <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Interface_Characteristic_Length_2") == 0) {
    i_command = 603;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Interface_IP.Component_List[IP.Interface_IP.ci].Length2;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Interface_IP.Component_List[IP.Interface_IP.ci].Length2 <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Interface_BC_Type") == 0) {
    i_command = 604;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Interface_IP.BC_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Interface_IP.BC_Type,"Reflection") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].BC_Type = INTERFACE_BC_REFLECTION;
    } else if (strcmp(IP.Interface_IP.BC_Type,"Absorption") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].BC_Type = INTERFACE_BC_ABSORPTION;
    } else if (strcmp(IP.Interface_IP.BC_Type,"Wall_Viscous_Heatflux") == 0 ||
	       strcmp(IP.Interface_IP.BC_Type,"Wall_Viscous_Adiabatic") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].BC_Type = INTERFACE_BC_WALL_VISCOUS_HEATFLUX;
    } else if (strcmp(IP.Interface_IP.BC_Type,"Wall_Viscous_Isothermal") == 0 ||
	       strcmp(IP.Interface_IP.BC_Type,"Wall_Viscous_Fixed_Temperature") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].BC_Type = INTERFACE_BC_WALL_VISCOUS_ISOTHERMAL; 
    } else if (strcmp(IP.Interface_IP.BC_Type,"Burning_Surface") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].BC_Type = INTERFACE_BC_BURNING_SURFACE;
    } else if (strcmp(IP.Interface_IP.BC_Type,"Ringleb_Flow") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].BC_Type = INTERFACE_BC_RINGLEB;
    } else if (strcmp(IP.Interface_IP.BC_Type,"Flat_Plate") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].BC_Type = INTERFACE_BC_FLAT_PLATE;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Interface_Reference_Point") == 0) {
    i_command = 605;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Interface_IP.Component_List[IP.Interface_IP.ci].Xref.x
		  >> IP.Interface_IP.Component_List[IP.Interface_IP.ci].Xref.y;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"Interface_Motion") == 0) {
    i_command = 606;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Interface_IP.Motion,IP.Next_Control_Parameter);
    if (strcmp(IP.Interface_IP.Motion,"Stationary") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_STATIONARY;
    } else if (strcmp(IP.Interface_IP.Motion,"Constant") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_CONSTANT;
    } else if (strcmp(IP.Interface_IP.Motion,"Expand") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_EXPAND;
    } else if (strcmp(IP.Interface_IP.Motion,"Uniform") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_UNIFORM;
    } else if (strcmp(IP.Interface_IP.Motion,"Translate") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_TRANSLATE;
    } else if (strcmp(IP.Interface_IP.Motion,"Stretch") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_STRETCH;
    } else if (strcmp(IP.Interface_IP.Motion,"Rotate") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_ROTATE;
    } else if (strcmp(IP.Interface_IP.Motion,"Burning_Surface") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_BURNING_SURFACE;
    } else if (strcmp(IP.Interface_IP.Motion,"Momentum_Transfer") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_MOMENTUM_TRANSFER;
    } else if (strcmp(IP.Interface_IP.Motion,"LevelSet_Stationary") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_LEVELSET_STATIONARY;
    } else if (strcmp(IP.Interface_IP.Motion,"LevelSet_Expand") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_LEVELSET_EXPAND;
    } else if (strcmp(IP.Interface_IP.Motion,"LevelSet_Stretch") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_LEVELSET_STRETCH;
    } else if (strcmp(IP.Interface_IP.Motion,"LevelSet_BulkFlow") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_LEVELSET_BULKFLOW;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Interface_Characteristic_Velocity") == 0) {
    i_command = 607;
    IP.Line_Number = IP.Line_Number + 1;
    //IP.Input_File >> IP.Interface_IP.Component_List[IP.Interface_IP.ci].Speed.x
    //	  >> IP.Interface_IP.Component_List[IP.Interface_IP.ci].Speed.y;
    IP.Input_File >> IP.Interface_IP.Component_List[IP.Interface_IP.ci].Speed;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Specified_Interface_Points") == 0) {
    i_command = 608;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Interface_IP.ns[IP.Interface_IP.ci-1];
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Interface_IP.ns[IP.Interface_IP.ci-1] <  0) i_command = INVALID_INPUT_VALUE;
    else IP.Interface_IP.Component_List[IP.Interface_IP.ci].allocate(max(IP.Interface_IP.ns[IP.Interface_IP.ci-1],2));

  } else if (strcmp(IP.Next_Control_Parameter,"Interface_Points") == 0) {
    i_command = 609;
    for (int np = 0; np < IP.Interface_IP.ns[IP.Interface_IP.ci-1]; np++) {
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> Xt.x >> Xt.y >> tpt >> bct;
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Spline.Xp[np] = Xt;
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Spline.tp[np] = tpt;
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Spline.bc[np] = bct;
      IP.Input_File.getline(buffer,sizeof(buffer));
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Solution_Redistribution_Type") == 0) {
    i_command = 610;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Interface_IP.Solution_Redistribution_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Interface_IP.Solution_Redistribution_Type,"None") == 0) {
      IP.Interface_IP.Redistribution_Type = SOLUTION_REDISTRIBUTION_NONE;
    } else if (strcmp(IP.Interface_IP.Solution_Redistribution_Type,"Injection") == 0) {
      IP.Interface_IP.Redistribution_Type = SOLUTION_REDISTRIBUTION_INJECTION;
    } else if (strcmp(IP.Interface_IP.Solution_Redistribution_Type,"Weighted_Injection") == 0) {
      IP.Interface_IP.Redistribution_Type = SOLUTION_REDISTIBUTION_WEIGHTED_INJECTION;
    } else if (strcmp(IP.Interface_IP.Solution_Redistribution_Type,"Neighbour_Average") == 0) {
      IP.Interface_IP.Redistribution_Type = SOLUTION_REDISTIBUTION_NEIGHBOUR_AVERAGE;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Evolution_Frequency") == 0) {
    i_command = 611;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Interface_IP.Evolution_Frequency;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"Reset_Interface_Motion_Type") == 0) {
    i_command = 612;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Interface_IP.BC_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Interface_IP.BC_Type,"ON") == 0) {
      IP.Reset_Interface_Motion_Type = ON;
    } else if (strcmp(IP.Interface_IP.BC_Type,"OFF") == 0) {
      IP.Reset_Interface_Motion_Type = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Execute") == 0) {
    i_command = EXECUTE_CODE;
 
  } else if (strcmp(IP.Next_Control_Parameter,"Execute_Zero_Steps") == 0) {
    i_command = EXECUTE_ZERO_STEPS_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Terminate") == 0) {
    i_command = TERMINATE_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Continue") == 0) {
    i_command = CONTINUE_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output") == 0) {
    i_command = WRITE_OUTPUT_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Cells") == 0) {
    i_command = WRITE_OUTPUT_CELLS_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Nodes") == 0) {
    i_command = WRITE_OUTPUT_NODES_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Gradients") == 0) {
    i_command = WRITE_OUTPUT_GRADIENTS_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Quasi3D") == 0) {
    i_command = WRITE_OUTPUT_QUASI3D_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Elements") == 0) {
    i_command = WRITE_OUTPUT_ELEMENTS_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Cell_Status") == 0) {
    i_command = WRITE_OUTPUT_CELL_STATUS_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Interface_Component_List") == 0) {
    i_command = WRITE_OUTPUT_INTERFACE_COMPONENT_LIST_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Interface_Union_List") == 0) {
    i_command = WRITE_OUTPUT_INTERFACE_UNION_LIST_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Level_Set") == 0) {
    i_command = WRITE_OUTPUT_LEVEL_SET_CODE;
    
  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Level_Set_Cells") == 0) {
    i_command = WRITE_OUTPUT_LEVEL_SET_CELLS_CODE;
    
  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Level_Set_Interface_List") == 0) {
    i_command = WRITE_OUTPUT_LEVEL_SET_INTERFACE_LIST_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Restart") == 0) {
    i_command = WRITE_RESTART_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Mesh") == 0) {
    i_command = WRITE_OUTPUT_GRID_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Mesh_Definition") == 0) {
    i_command = WRITE_GRID_DEFINITION_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Mesh_Nodes") == 0) {
    i_command = WRITE_OUTPUT_GRID_NODES_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Mesh_Cells") == 0) {
    i_command = WRITE_OUTPUT_GRID_CELLS_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Ringleb") == 0) {
    i_command = WRITE_OUTPUT_RINGLEB_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Viscous_Channel") == 0) {
    i_command = WRITE_OUTPUT_VISCOUS_CHANNEL_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Viscous_Pipe") == 0) {
    i_command = WRITE_OUTPUT_VISCOUS_PIPE_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Flat_Plate") == 0) {
    i_command = WRITE_OUTPUT_FLAT_PLATE_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Driven_Cavity_Flow") == 0) {
    i_command = WRITE_OUTPUT_DRIVEN_CAVITY_FLOW_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Turbulent_Pipe") == 0) {
    i_command = WRITE_OUTPUT_TURBULENT_PIPE_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Backward_Facing_Step") == 0) {
    i_command = WRITE_OUTPUT_BACKWARD_FACING_STEP_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Aerodynamic_Coefficients") == 0) {
    i_command = WRITE_OUTPUT_AERODYNAMIC_COEFFICIENTS_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Elliptic_Operator_Analysis") == 0) {
    i_command = WRITE_OUTPUT_ELLIPTIC_OPERATOR_ANALYSIS;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Supersonic_Hot_Jet") == 0) {
    i_command = WRITE_OUTPUT_SUPERSONIC_HOT_JET_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Subsonic_Hot_Jet") == 0) {
    i_command = WRITE_OUTPUT_SUBSONIC_HOT_JET_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Refine_Grid") == 0) {
    i_command = REFINE_GRID_CODE;

  } else if (strcmp(IP.Next_Control_Parameter, "Morton_Ordering") == 0) {
     i_command = MORTON_ORDERING_CODE;

  } else if (strcmp(IP.Next_Control_Parameter, "Bounding_Box_Refine_Grid") == 0) {
    IP.Number_of_Bounding_Box_Mesh_Refinements = 1;
    i_command = BOUNDING_BOX_REFINE_GRID_CODE;

  } else if (IP.Next_Control_Parameter[0] == '#') {
    i_command = COMMENT_CODE;

  } else {
    i_command = INVALID_INPUT_CODE;

  }

  if (i_command == INVALID_INPUT_CODE) {
    strcpy(buffer, IP.Next_Control_Parameter);
    Get_Next_Input_Control_Parameter(IP);
    i_command = IP.NKS_IP.Parse_Next_Input_Control_Parameter(buffer, 
                                                             IP.Next_Control_Parameter);

//     if (i_command == INVALID_INPUT_CODE) {
//       cout << "\n***\n\nWarning: input file line " << IP.Line_Number << ": ";
//       cout << "ignoring unknown input code:\n";
//       cout << "code: " << buffer;
//       cout << "\nvalue: " << IP.Next_Control_Parameter;
//       cout << "\n\n***\n";
//     }
//     i_command = COMMENT_CODE;
  }

  if (!IP.Input_File.good()) { i_command = INVALID_INPUT_VALUE; }

  // Return the parser command type indicator.
  return i_command;

}

/**********************************************************************
 * Routine: Process_Input_Control_Parameter_File                      *
 *                                                                    *
 * Reads, parses, and executes the list of input control parameters   *
 * from the standard input file.                                      *
 *                                                                    *
 **********************************************************************/
int Process_Input_Control_Parameter_File(NavierStokes2D_Input_Parameters &IP,
                                         char *Input_File_Name_ptr,
                                         int &Command_Flag) {

  int error_flag, line_number;

  // Assign default values to the input parameters.
  Set_Default_Input_Parameters(IP);

  // Copy input file name (a string) to appropriate input parameter variable.
  if (Input_File_Name_ptr != NULL) strcpy(IP.Input_File_Name,Input_File_Name_ptr);

  // Open the input file containing the input parameters.
  Open_Input_File(IP);
  error_flag = IP.Input_File.fail();

  if (error_flag) {
    cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D input data file.\n";
    return error_flag;
  }

  // Read and parse control parameters contained in the input file.
  while (1) {
    Get_Next_Input_Control_Parameter(IP);
    Command_Flag = Parse_Next_Input_Control_Parameter(IP);
    line_number = IP.Line_Number;
    if (Command_Flag == EXECUTE_CODE || Command_Flag == EXECUTE_ZERO_STEPS_CODE) {
      break;
    } else if (Command_Flag == TERMINATE_CODE) {
      break;
    } else if (Command_Flag == INVALID_INPUT_CODE ||
	       Command_Flag == INVALID_INPUT_VALUE) {
      line_number = -line_number;
      cout << "\n NavierStokes2D ERROR: Error reading NavierStokes2D data at line #"
	   << -line_number  << " of input data file.\n";
      error_flag = line_number;
      break;
    }
  }
  if (error_flag) return error_flag;

  // Set static variables and initialize reference state.
  Initialize_Reference_State(IP);

  // Perform consistency checks on the input parameters.
  if (IP.i_Time_Integration == TIME_STEPPING_MULTIGRID) {
    error_flag = Check_Input_Parameters<NavierStokes2D_Input_Parameters>(IP);
    if (error_flag) {
      cout << "\n NavierStokes2D ERROR: Input Parameters consistency check failure\n";
      return error_flag;
    }
  }

  // Perform consitency checks on the refinement criteria.
  IP.Number_of_Refinement_Criteria = 0;
  if (IP.Refinement_Criteria_Gradient_Density) IP.Number_of_Refinement_Criteria++;
  if (IP.Refinement_Criteria_Divergence_Velocity) IP.Number_of_Refinement_Criteria++;
  if (IP.Refinement_Criteria_Curl_Velocity) IP.Number_of_Refinement_Criteria++;
  if (IP.Refinement_Criteria_Gradient_Turbulence_Kinetic_Energy) IP.Number_of_Refinement_Criteria++;
  if (IP.Number_of_Refinement_Criteria < 1 || IP.Number_of_Refinement_Criteria > 4) return 1011;

  // Initial processing of input control parameters complete.  
  return 0;

}

/**********************************************************************
 * Routine: Initialize_Reference_State                                *
 *                                                                    *
 * This function sets the static variables and reference states.      *
 *                                                                    *
 **********************************************************************/
void Initialize_Reference_State(NavierStokes2D_Input_Parameters &IP) {

  double cos_angle, sin_angle;

  // Set static variables.
  IP.Wo.set_static_variables(IP.Gas_Type,
			     IP.FlowType,
			     IP.C_constant,
			     IP.von_Karman_Constant,
			     IP.yplus_sublayer,
			     IP.yplus_buffer_layer,
			     IP.yplus_outer_layer,
			     IP.Propellant_Type,
			     IP.i_Compressibility_Effect,
			     IP.i_Transition_Model, 
			     IP.i_Variable_Prandtl,
			     IP.sigmav,IP.lw,
			     IP.C_lambda,
			     IP.Cd1,IP.Cd4,
			     IP.Cd5);
  IP.Uo.set_static_variables(IP.Gas_Type,
			     IP.FlowType,
			     IP.C_constant,
			     IP.von_Karman_Constant,
			     IP.yplus_sublayer,
			     IP.yplus_buffer_layer,
			     IP.yplus_outer_layer,
			     IP.Propellant_Type,
			     IP.i_Compressibility_Effect,
			     IP.i_Transition_Model, 
			     IP.i_Variable_Prandtl,
			     IP.sigmav,IP.lw,
			     IP.C_lambda,
			     IP.Cd1,IP.Cd4,
			     IP.Cd5);

  // Initialize gas-phase reference state.
  IP.Wo.rho = IP.Pressure/(IP.Wo.R*IP.Temperature);
  IP.Wo.p   = IP.Pressure;
  cos_angle = cos(TWO*PI*IP.Flow_Angle/360.00); if (fabs(cos_angle) < TOLER*TOLER) cos_angle = ZERO;
  sin_angle = sin(TWO*PI*IP.Flow_Angle/360.00); if (fabs(sin_angle) < TOLER*TOLER) sin_angle = ZERO;
  IP.Wo.v.x = IP.Mach_Number*IP.Wo.a()*cos_angle;
  IP.Wo.v.y = IP.Mach_Number*IP.Wo.a()*sin_angle;
  IP.Uo = U(IP.Wo);

  // NASA rotor variables.
  if (strcmp(IP.Grid_Type,"NASA_Rotor_37") == 0) {
    IP.NASA_Rotor37.init(IP.Rotor_Flow_Type,IP.NASA_Rotor37_Data_Directory);
    IP.NASA_Rotor37.getPstateREL_up(IP.Wo,IP.Rotor_Percent_Span);
    IP.Pressure = IP.Wo.p;
    IP.Temperature = IP.Wo.T();
    IP.Mach_Number = IP.NASA_Rotor37.getMachREL_up(IP.Rotor_Percent_Span);
    IP.Flow_Angle = atan2(IP.Wo.v.y,IP.Wo.v.x); 
    if (IP.Flow_Angle < ZERO) IP.Flow_Angle += TWO*PI;
    IP.Flow_Angle = 180.00*IP.Flow_Angle/PI;
    IP.NASA_Rotor37.getPstateREL_down(IP.W1,IP.Rotor_Percent_Span);
  } else if (strcmp(IP.Grid_Type,"NASA_Rotor_67") == 0) {
    IP.NASA_Rotor67.init(IP.Rotor_Flow_Type,IP.NASA_Rotor67_Data_Directory);
    IP.NASA_Rotor67.getPstateREL_up(IP.Wo,IP.Rotor_Percent_Span);
    IP.Pressure = IP.Wo.p;
    IP.Temperature = IP.Wo.T();
    IP.Mach_Number = IP.NASA_Rotor67.getMachREL_up(IP.Rotor_Percent_Span);
    IP.Flow_Angle = atan2(IP.Wo.v.y,IP.Wo.v.x); 
    if (IP.Flow_Angle < ZERO) IP.Flow_Angle += TWO*PI;
    IP.Flow_Angle = 180.00*IP.Flow_Angle/PI;
    IP.NASA_Rotor67.getPstateREL_down(IP.W1,IP.Rotor_Percent_Span);
  }

  // Set the state for the laminar flat-plate boundary layer flow.
  if (IP.i_Grid == GRID_FLAT_PLATE && (IP.FlowType == FLOWTYPE_LAMINAR ||
				       IP.FlowType == FLOWTYPE_TURBULENT_RANS_K_OMEGA)) {
    IP.Wo.rho = (IP.Reynolds_Number/IP.Plate_Length)*IP.Wo.mu()/(sqrt(IP.Wo.g*IP.Wo.R*IP.Temperature)*IP.Mach_Number);
    IP.Wo.p = IP.Wo.rho*IP.Wo.R*IP.Temperature;
    IP.Wo.v.x = sqrt(IP.Wo.g*IP.Wo.R*IP.Temperature)*IP.Mach_Number;
    IP.Wo.v.y = ZERO;
    IP.Uo = U(IP.Wo);
    //IP.Twall = 1.2*IP.Temperature;
  }

  if (IP.Interface_IP.Number_of_Components) {
    if (IP.Interface_IP.Component_List[1].Type == INTERFACE_FLAT_PLATE) {
      IP.Plate_Length = IP.Interface_IP.Component_List[1].Length1;
      IP.Wo.rho = (IP.Reynolds_Number/IP.Plate_Length)*IP.Wo.mu()/(sqrt(IP.Wo.g*IP.Wo.R*IP.Temperature)*IP.Mach_Number);
      IP.Wo.p = IP.Wo.rho*IP.Wo.R*IP.Temperature;
      IP.Wo.v.x = sqrt(IP.Wo.g*IP.Wo.R*IP.Temperature)*IP.Mach_Number;
      IP.Wo.v.y = ZERO;
      IP.Uo = U(IP.Wo);
      strcpy(IP.Grid_Type,"Rectangular_Box");
      IP.i_Grid = GRID_RECTANGULAR_BOX;
      IP.Box_Width = TWO*(IP.Interface_IP.Component_List[1].Length2*
			  IP.Interface_IP.Component_List[1].Length1*
			  cos(TWO*PI*IP.X_Rotate/360.00));
      IP.Box_Height = IP.Box_Width;
      IP.X_Shift = Vector2D(ZERO,
			    HALF*IP.Box_Height +
			    IP.Interface_IP.Component_List[1].Length2*
			    IP.Interface_IP.Component_List[1].Length1*
			    sin(TWO*PI*IP.X_Rotate/360.00));
    }
  }

  // Set the driven cavity flow lid speed.
  if (IP.i_ICs == IC_VISCOUS_DRIVEN_CAVITY_FLOW || IP.i_Grid == GRID_DRIVEN_CAVITY_FLOW) {
    IP.Box_Width = 0.001; IP.Box_Height = 0.001;
    IP.Wo = DrivenCavityFlow(IP.Wo,IP.Box_Width,IP.Re_lid);
    IP.Vwall.x = IP.Re_lid*IP.Wo.nu()/IP.Box_Width;
    IP.Vwall.y = ZERO;
    IP.Uo = U(IP.Wo);
  }

  // Set the maximum mean flow velocity for the turbulent pipe flow.
  if (IP.i_Grid == GRID_PIPE && IP.FlowType > FLOWTYPE_LAMINAR) {
    IP.Wo.v.x = IP.Reynolds_Number*IP.Wo.nu()/(TWO*IP.Pipe_Radius);
    IP.Wo.v.y = ZERO;
    IP.Uo = U(IP.Wo);
  }

  // Set the state for the mixing layer case.
  if (IP.i_ICs == IC_MIXING_LAYER) {
    IP.Wo.rho = (IP.Wo.p/(TEMPERATURE_STDATM*IP.Wo.R))*sqrt(ONE + HALF*IP.Wo.gm1*sqr(IP.Mach_Number));
    IP.Box_Width = (IP.Reynolds_Number/IP.Wo.rho)*IP.Wo.mu()/(IP.Wo.a()*IP.Mach_Number);
    IP.Wo.v.x = IP.Wo.a()*IP.Mach_Number;
    IP.Wo.v.y = sqrt(IP.Wo.g*IP.Wo.R*TEMPERATURE_STDATM)*IP.Mach_Number2;
    IP.Uo = U(IP.Wo);
  }

}

/**********************************************************************
 * Routine: Reinitialize_Reference_State                              *
 *                                                                    *
 * This function sets the static variables and reference states.      *
 *                                                                    *
 **********************************************************************/
void Reinitialize_Reference_State(NavierStokes2D_Input_Parameters &IP) {

  // NASA rotor variables.
  if (strcmp(IP.Grid_Type,"NASA_Rotor_37") == 0) {
    IP.NASA_Rotor37.init(IP.Rotor_Flow_Type,IP.NASA_Rotor37_Data_Directory);
    IP.NASA_Rotor37.getPstateREL_up(IP.Wo,IP.Rotor_Percent_Span);
    IP.Pressure = IP.Wo.p;
    IP.Temperature = IP.Wo.T();
    IP.Mach_Number = IP.NASA_Rotor37.getMachREL_up(IP.Rotor_Percent_Span);
    IP.Flow_Angle = atan2(IP.Wo.v.y,IP.Wo.v.x); 
    if (IP.Flow_Angle < ZERO) IP.Flow_Angle += TWO*PI;
    IP.Flow_Angle = 180.00*IP.Flow_Angle/PI;
    IP.NASA_Rotor37.getPstateREL_down(IP.W1,IP.Rotor_Percent_Span);
  } else if (strcmp(IP.Grid_Type,"NASA_Rotor_67") == 0) {
    IP.NASA_Rotor67.init(IP.Rotor_Flow_Type,IP.NASA_Rotor67_Data_Directory);
    IP.NASA_Rotor67.getPstateREL_up(IP.Wo,IP.Rotor_Percent_Span);
    IP.Pressure = IP.Wo.p;
    IP.Temperature = IP.Wo.T();
    IP.Mach_Number = IP.NASA_Rotor67.getMachREL_up(IP.Rotor_Percent_Span);
    IP.Flow_Angle = atan2(IP.Wo.v.y,IP.Wo.v.x); 
    if (IP.Flow_Angle < ZERO) IP.Flow_Angle += TWO*PI;
    IP.Flow_Angle = 180.00*IP.Flow_Angle/PI;
    IP.NASA_Rotor67.getPstateREL_down(IP.W1,IP.Rotor_Percent_Span);
  }

  // Set the state for the laminar flat-plate boundary layer flow.
  if (IP.i_Grid == GRID_FLAT_PLATE && (IP.FlowType == FLOWTYPE_LAMINAR ||
				       IP.FlowType == FLOWTYPE_TURBULENT_RANS_K_OMEGA)) {
    IP.Wo.rho = (IP.Reynolds_Number/IP.Plate_Length)*IP.Wo.mu()/(sqrt(IP.Wo.g*IP.Wo.R*IP.Temperature)*IP.Mach_Number);
    IP.Wo.p = IP.Wo.rho*IP.Wo.R*IP.Temperature;
    IP.Wo.v.x = sqrt(IP.Wo.g*IP.Wo.R*IP.Temperature)*IP.Mach_Number;
    IP.Wo.v.y = ZERO;
    IP.Uo = U(IP.Wo);
    //IP.Twall = 1.2*IP.Temperature;
  }

  if (IP.Interface_IP.Number_of_Components) {
    if (IP.Interface_IP.Component_List[1].Type == INTERFACE_FLAT_PLATE) {
      IP.Plate_Length = IP.Interface_IP.Component_List[1].Length1;
      IP.Wo.rho = (IP.Reynolds_Number/IP.Plate_Length)*IP.Wo.mu()/(sqrt(IP.Wo.g*IP.Wo.R*IP.Temperature)*IP.Mach_Number);
      IP.Wo.p = IP.Wo.rho*IP.Wo.R*IP.Temperature;
      IP.Wo.v.x = sqrt(IP.Wo.g*IP.Wo.R*IP.Temperature)*IP.Mach_Number;
      IP.Wo.v.y = ZERO;
      IP.Uo = U(IP.Wo);
      strcpy(IP.Grid_Type,"Rectangular_Box");
      IP.i_Grid = GRID_RECTANGULAR_BOX;
      IP.Box_Width = TWO*(IP.Interface_IP.Component_List[1].Length2*
			  IP.Interface_IP.Component_List[1].Length1*
			  cos(TWO*PI*IP.X_Rotate/360.00));
      IP.Box_Height = IP.Box_Width;
      IP.X_Shift = Vector2D(ZERO,
			    HALF*IP.Box_Height +
			    IP.Interface_IP.Component_List[1].Length2*
			    IP.Interface_IP.Component_List[1].Length1*
			    sin(TWO*PI*IP.X_Rotate/360.00));
    }
  }

  // Set the driven cavity flow lid speed.
  if (IP.i_ICs == IC_VISCOUS_DRIVEN_CAVITY_FLOW || IP.i_Grid == GRID_DRIVEN_CAVITY_FLOW) {
    IP.Box_Width = 0.001; IP.Box_Height = 0.001;
    IP.Wo = DrivenCavityFlow(IP.Wo,IP.Box_Width,IP.Re_lid);
    IP.Vwall.x = IP.Re_lid*IP.Wo.nu()/IP.Box_Width;
    IP.Vwall.y = ZERO;
    IP.Uo = U(IP.Wo);
  }

  // Set the maximum mean flow velocity for the turbulent pipe flow.
  if (IP.i_Grid == GRID_PIPE && IP.FlowType > FLOWTYPE_LAMINAR) {
    IP.Wo.v.x = IP.Reynolds_Number*IP.Wo.nu()/(TWO*IP.Pipe_Radius);
    IP.Wo.v.y = ZERO;
    IP.Uo = U(IP.Wo);
  }

  // Set the state for the mixing layer case.
  if (IP.i_ICs == IC_MIXING_LAYER) {
    IP.Wo.rho = (IP.Wo.p/(TEMPERATURE_STDATM*IP.Wo.R))*sqrt(ONE + HALF*IP.Wo.gm1*sqr(IP.Mach_Number));
    IP.Box_Width = (IP.Reynolds_Number/IP.Wo.rho)*IP.Wo.mu()/(IP.Wo.a()*IP.Mach_Number);
    IP.Wo.v.x = IP.Wo.a()*IP.Mach_Number;
    IP.Wo.v.y = sqrt(IP.Wo.g*IP.Wo.R*TEMPERATURE_STDATM)*IP.Mach_Number2;
    IP.Uo = U(IP.Wo);
  }

}
