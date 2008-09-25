/*!\file Euler2DInput.cc
  \brief Subroutines for 2D Euler Input Classes. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "Euler2DInput.h" /* Include 2D Euler input parameter header file. */
#include "Euler2DQuad.h"  /* Include Euler2D_Quad_Block header file. */
#include "../Grid/HO_Grid2DQuad_ExecutionMode.h" // Include high-order 2D grid execution mode header file
#include "../Grid/Tecplot_ExecutionMode.h" // Include Tecplot execution mode header file

/*****************************************************
 * Euler2D_Input_Parameters -- Member functions.     *
 ****************************************************/

/*!
 * Decide whether to output or not the boundary reference state for a particular boundary condition.
 * To get output for a particular BCtype, just add it to the list.
 */
bool Euler2D_Input_Parameters::OutputBoundaryReferenceState(const int & BCtype) const{
  if (BCtype == BC_DIRICHLET ||
      BCtype == BC_NEUMANN ||
      BCtype == BC_FARFIELD){
    // Get output
    return true;
  } else {
    // No output
    return false;
  }
}

/******************************************************//**
 * Parse the input file
 ********************************************************/
int Euler2D_Input_Parameters::Parse_Input_File(char *Input_File_Name_ptr){

  ostringstream msg;
  int command_flag, error_flag;

  /* Assign initial value for error indicator flag. */
  error_flag = 0;

  strcpy(Input_File_Name, Input_File_Name_ptr);
  Open_Input_File(*this);
  if (Input_File.fail()) {
    msg << "Euler2D_Input_Parameters::Parse_Input_File() ERROR: Unable to open "
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
      
      msg << "Euler2D_Input_Parameters::Parse_Input_File() ERROR: Error reading data at line # " 
	  << -Line_Number
	  << " of input data file.";
      if (Verbose()){
	cerr << msg.str() << endl;
      }
      
      throw runtime_error(msg.str());
    } /* endif */
  } /* endwhile */

  /* Perform consistency checks on Input_Parameters */
  if (!error_flag && 
      i_Time_Integration == TIME_STEPPING_MULTIGRID) {
    error_flag = Check_Input_Parameters<Euler2D_Input_Parameters>(*this);
    if (error_flag) {
      cout << "\n Euler2D ERROR: Input Parameters consistency check failure\n";
      return (error_flag);
    }
  }

  // Perform consitency checks on the refinement criteria.
  Number_of_Refinement_Criteria = 0;
  if (Refinement_Criteria_Gradient_Density) Number_of_Refinement_Criteria++;
  if (Refinement_Criteria_Divergence_Velocity) Number_of_Refinement_Criteria++;
  if (Refinement_Criteria_Curl_Velocity) Number_of_Refinement_Criteria++;

  // Perform consitency checks on the time marching parameters.
  if (Time_Accurate == 1 && Local_Time_Stepping != GLOBAL_TIME_STEPPING){
    Local_Time_Stepping = GLOBAL_TIME_STEPPING;
  }
  
  // Set Wo state based on the final input parameters
  Wo = Euler2D_pState(Pressure/(Wo.R*Temperature), 
		      ZERO, 
		      ZERO, 
		      Pressure);
  Wo.v.x = Mach_Number*Wo.a()*cos(TWO*PI* Flow_Angle/360.00);
  Wo.v.y = Mach_Number*Wo.a()*sin(TWO*PI* Flow_Angle/360.00);
    
  Uo.setgas(Gas_Type);
  Uo = U(Wo);

  Wo.Mr_min = Mr_Min_Factor*Mach_Number;
  Uo.Mr_min = Mr_Min_Factor*Mach_Number;
  
  // Perform update of the internal variables of the exact solution
  ExactSoln->Set_ParticularSolution_Parameters(*this);
  
  // Perform update of the internal variables of the high-order input parameters
  HighOrder2D_Input::Set_Final_Parameters(*this);
  
  // Set reference state in the Euler2D_Quad_Block class
  Euler2D_Quad_Block::Set_Normalization_Reference_State(RefU);
  
  // Set limiter in CENO class
  CENO_Execution_Mode::Limiter = i_Limiter;

  /* Initial processing of input control parameters complete.  
     Return the error indicator flag. */
  return (error_flag);

}

/******************************************************//**
 * Get the next input control parameter from the input file.                                                
 ********************************************************/
void Euler2D_Input_Parameters::Get_Next_Input_Control_Parameter(void){

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


/*************************************************************
 * Euler2D_Input_Parameters -- External subroutines.         *
 *************************************************************/

/******************************************************//**
 * Routine: Open_Input_File                             
 *                                                      
 * Opens the appropriate input data file.               
 *                                                      
 ********************************************************/
void Open_Input_File(Euler2D_Input_Parameters &IP) {

    IP.Input_File.open(IP.Input_File_Name, ios::in);
    if (!IP.Input_File.fail()) {
       IP.Line_Number = 0;
       IP.Input_File.setf(ios::skipws);
    } /* endif */

}

/******************************************************//**
 * Routine: Close_Input_File                            
 *                                                      
 * Closes the appropriate input data file.              
 *                                                      
 ********************************************************/
void Close_Input_File(Euler2D_Input_Parameters &IP) {

    IP.Input_File.unsetf(ios::skipws);
    IP.Input_File.close();

}

/******************************************************//**
 * Routine: Set_Default_Input_Parameters               
 *                                                     
 * Assigns default values to the input parameters.     
 *                                                     
 ********************************************************/
void Set_Default_Input_Parameters(Euler2D_Input_Parameters &IP) {

    int i;
    char *string_ptr;

    // CFFC root directory path:
    IP.get_cffc_path();

    string_ptr = "Euler2D.in";
    strcpy(IP.Input_File_Name, string_ptr);

    string_ptr = "Explicit_Euler";
    strcpy(IP.Time_Integration_Type, string_ptr);
    IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
    IP.Time_Accurate = 0;
    IP.Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
    IP.Maximum_Number_of_Time_Steps = 100;
    IP.N_Stage = 1;
    IP.Time_Max = ZERO;
    IP.Mr_Min_Factor = ONE;

    // Residual variable:
    IP.i_Residual_Variable = 1;

    IP.Residual_Smoothing = 0;
    IP.Residual_Smoothing_Epsilon = ZERO;
    IP.Residual_Smoothing_Gauss_Seidel_Iterations = 2;

    // Reconstruction type:
    string_ptr = "Least_Squares";
    strcpy(IP.Reconstruction_Type, string_ptr);
    IP.i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;
    IP.Space_Accuracy = 1;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_ReconstructionMethod = RECONSTRUCTION_LEAST_SQUARES;
    CENO_Execution_Mode::USE_CENO_ALGORITHM = OFF;

    // Limiter type:
    string_ptr = "Barth_Jespersen";
    strcpy(IP.Limiter_Type, string_ptr);
    IP.i_Limiter = LIMITER_BARTH_JESPERSEN;

    // Flux function:
    string_ptr = "HLLE";
    strcpy(IP.Flux_Function_Type, string_ptr);
    IP.i_Flux_Function = FLUX_FUNCTION_HLLE;

    // Initial conditions:
    string_ptr = "Uniform";
    strcpy(IP.ICs_Type, string_ptr);
    IP.i_ICs = IC_UNIFORM;
    IP.Exact_Integration_Digits = 9;

    string_ptr = "AIR";
    strcpy(IP.Gas_Type, string_ptr);
    IP.Wo.setgas(IP.Gas_Type);
    IP.Wo = Euler2D_W_STDATM;
    IP.Pressure = IP.Wo.p;
    IP.Temperature = IP.Wo.T();
    IP.Mach_Number = 0.80;
    IP.Flow_Angle = ZERO;
    IP.Wo.v.x = IP.Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Flow_Angle/360.00);
    IP.Wo.v.y = IP.Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Flow_Angle/360.00);
    IP.Uo.setgas(IP.Gas_Type);
    IP.Uo = U(IP.Wo);
    IP.W1 = IP.Wo;
    IP.W2 = IP.Wo;
    IP.Wo.Mr_min = IP.Mr_Min_Factor*IP.Mach_Number;
    IP.Uo.Mr_min = IP.Mr_Min_Factor*IP.Mach_Number;
    IP.Wave_Position = Vector2D_ZERO;
    IP.Wave_Width = ZERO;

    // Geometry switch:
    string_ptr = "Planar";
    strcpy(IP.Flow_Geometry_Type, string_ptr);
    IP.Axisymmetric = 0;

    // Grid parameters:
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
    IP.Cylinder_Radius2 = 32.00;
    IP.Annulus_Theta_Start = 0.0;
    IP.Annulus_Theta_End = 90.0;
    IP.Ellipse_Length_X_Axis = TWO;
    IP.Ellipse_Length_Y_Axis = HALF;
    IP.Chord_Length = ONE;
    IP.Orifice_Radius = ONE;
    IP.Inner_Streamline_Number = 0.80;
    IP.Outer_Streamline_Number = 0.40;
    IP.Isotach_Line = 0.30;
    IP.Wedge_Angle = 25.0;
    IP.Wedge_Length = HALF;
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

    // Mesh stretching factor:
    IP.i_Mesh_Stretching = OFF;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_LINEAR;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_LINEAR;
    IP.Mesh_Stretching_Factor_Idir = 1.01;
    IP.Mesh_Stretching_Factor_Jdir = 1.01;

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
    IP.Ref_State_BC_North = Euler2D_pState(0.0);
    IP.Ref_State_BC_South = Euler2D_pState(0.0);
    IP.Ref_State_BC_East = Euler2D_pState(0.0);
    IP.Ref_State_BC_West = Euler2D_pState(0.0);

    // NASA rotor input variables:
    strcpy(IP.NASA_Rotor37_Data_Directory, IP.CFFC_Path);
    strcpy(IP.NASA_Rotor37_Data_Directory, "/data/NASA_Rotors/R37/");
    strcpy(IP.NASA_Rotor67_Data_Directory, IP.CFFC_Path);
    strcpy(IP.NASA_Rotor67_Data_Directory, "/data/NASA_Rotors/R67/");
    IP.Rotor_Flow_Type = PEAK_FLOW;
    IP.Rotor_Percent_Span = 50.00;

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
    IP.Interface_Refinement_Condition = OFF;
    IP.Maximum_Refinement_Level = 100;
    IP.Minimum_Refinement_Level = 1;
    IP.Threshold_for_Refinement = 0.50;
    IP.Threshold_for_Coarsening = 0.10;
    IP.Number_of_Refinement_Criteria = 3;
    IP.Refinement_Criteria_Gradient_Density = ON;
    IP.Refinement_Criteria_Divergence_Velocity = ON;
    IP.Refinement_Criteria_Curl_Velocity = ON;
    IP.AMR_Xmin = Vector2D_ZERO;
    IP.AMR_Xmax = Vector2D_ZERO;

    IP.Solver_Type = EXPLICIT;

    IP.Morton = 0;
    IP.Morton_Reordering_Frequency = 1000;

    // Smooth quad block flag:
    IP.i_Smooth_Quad_Block = ON;

    // Embedded boundary parameters:
    IP.Reset_Interface_Motion_Type = OFF;

    // Default output file names and parameters:
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

    // Input_file parameters:
    string_ptr = " ";
    strcpy(IP.Next_Control_Parameter, string_ptr);

    IP.Line_Number = 0;

    IP.Number_of_Processors = CFFC_MPI::Number_of_Processors;
    IP.Number_of_Blocks_Per_Processor = 10;

    // Freezing_Limiter
    IP.Freeze_Limiter = 0;

    // Limiter_switch
    IP.Freeze_Limiter_Residual_Level = 1e-4;

    // Accuracy assessment parameters:
    IP.Accuracy_Assessment_Mode = ACCURACY_ASSESSMENT_BASED_ON_EXACT_SOLUTION;
    IP.Accuracy_Assessment_Exact_Digits = 10;
    IP.Accuracy_Assessment_Parameter = 1;

    // High-order parameters:
    HighOrder2D_Input::SetDefaults();
}

/******************************************************//**
 * Routine: Broadcast_Input_Parameters                  
 *                                                      
 * Broadcast the input parameters variables to all      
 * processors involved in the calculation from the      
 * primary processor using the MPI broadcast routine.   
 *                                                      
 ********************************************************/
void Broadcast_Input_Parameters(Euler2D_Input_Parameters &IP) {

#ifdef _MPI_VERSION
    int i;

    // CFFC path:
    MPI::COMM_WORLD.Bcast(IP.CFFC_Path, 
 			  INPUT_PARAMETER_LENGTH_EULER2D, 
			  MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Input_File_Name, 
                          INPUT_PARAMETER_LENGTH_EULER2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Line_Number), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Time_Integration_Type, 
                          INPUT_PARAMETER_LENGTH_EULER2D, 
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
    MPI::COMM_WORLD.Bcast(&(IP.Mr_Min_Factor), 
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
                          INPUT_PARAMETER_LENGTH_EULER2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Reconstruction), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_ReconstructionMethod), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Space_Accuracy), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.IncludeHighOrderBoundariesRepresentation), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Limiter_Type, 
                          INPUT_PARAMETER_LENGTH_EULER2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Limiter), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Flux_Function_Type, 
                          INPUT_PARAMETER_LENGTH_EULER2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Flux_Function), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.ICs_Type, 
                          INPUT_PARAMETER_LENGTH_EULER2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Gas_Type, 
                          INPUT_PARAMETER_LENGTH_EULER2D, 
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
    MPI::COMM_WORLD.Bcast(&(IP.Flow_Angle), 
                          1, 
                          MPI::DOUBLE, 0);
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
       IP.Wo.setgas(IP.Gas_Type);
       IP.Wo = Euler2D_pState(IP.Pressure/(IP.Wo.R*IP.Temperature), 
                              ZERO, 
                              ZERO, 
                              IP.Pressure);
       IP.Wo.v.x = IP.Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Flow_Angle/360.00);
       IP.Wo.v.y = IP.Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Flow_Angle/360.00);
       IP.Uo.setgas(IP.Gas_Type);
       IP.Uo = U(IP.Wo);
       IP.Wo.Mr_min = IP.Mr_Min_Factor*IP.Mach_Number;
       IP.Uo.Mr_min = IP.Mr_Min_Factor*IP.Mach_Number;
    } /* endif */

    // Reference states
    MPI::COMM_WORLD.Bcast(&(IP.W1.d), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W1.v.x), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W1.v.y), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W1.p), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W2.d), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W2.v.x), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W2.v.y), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W2.p), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.RefU.d), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.RefU.v.x), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.RefU.v.y), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.RefU.p), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Exact_Integration_Digits), 
                          1, 
                          MPI::DOUBLE, 0);

    MPI::COMM_WORLD.Bcast(IP.Flow_Geometry_Type, 
                          INPUT_PARAMETER_LENGTH_EULER2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Axisymmetric), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Grid_Type, 
                          INPUT_PARAMETER_LENGTH_EULER2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.NACA_Aerofoil_Type, 
                          INPUT_PARAMETER_LENGTH_EULER2D, 
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
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Cylinder_Radius), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Cylinder_Radius2), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Annulus_Theta_Start), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Annulus_Theta_End), 
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
    MPI::COMM_WORLD.Bcast(&(IP.Inner_Streamline_Number),
			  1,
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Outer_Streamline_Number),
			  1,
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Isotach_Line),
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
    MPI::COMM_WORLD.Bcast(&(IP.VertexSW.x), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.VertexSW.y), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.VertexSE.x), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.VertexSE.y), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.VertexNW.x), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.VertexNW.y), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.VertexNE.x), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.VertexNE.y), 
                          1, 
                          MPI::DOUBLE, 0);
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
			  INPUT_PARAMETER_LENGTH_EULER2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(&(IP.BCs_Specified),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(IP.BC_North_Type,
			  INPUT_PARAMETER_LENGTH_EULER2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(IP.BC_South_Type,
			  INPUT_PARAMETER_LENGTH_EULER2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(IP.BC_East_Type,
			  INPUT_PARAMETER_LENGTH_EULER2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(IP.BC_West_Type,
			  INPUT_PARAMETER_LENGTH_EULER2D,
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
    MPI::COMM_WORLD.Bcast(&(IP.Ref_State_BC_North.d),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Ref_State_BC_North.v.x),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Ref_State_BC_North.v.y),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Ref_State_BC_North.p),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Ref_State_BC_South.d),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Ref_State_BC_South.v.x),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Ref_State_BC_South.v.y),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Ref_State_BC_South.p),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Ref_State_BC_East.d),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Ref_State_BC_East.v.x),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Ref_State_BC_East.v.y),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Ref_State_BC_East.p),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Ref_State_BC_West.d),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Ref_State_BC_West.v.x),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Ref_State_BC_West.v.y),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Ref_State_BC_West.p),
			  1,
			  MPI::DOUBLE,0);

    // NASA rotors:
    MPI::COMM_WORLD.Bcast(IP.NASA_Rotor37_Data_Directory, 
                          INPUT_PARAMETER_LENGTH_EULER2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.NASA_Rotor67_Data_Directory, 
                          INPUT_PARAMETER_LENGTH_EULER2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Rotor_Flow_Type), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Rotor_Percent_Span), 
                          1, 
                          MPI::DOUBLE, 0);

    // ICEM:
    if (!CFFC_Primary_MPI_Processor()) {
       IP.ICEMCFD_FileNames = new char*[3];
       for (i = 0; i < 3; i++) {
          IP.ICEMCFD_FileNames[i] = new char[INPUT_PARAMETER_LENGTH_EULER2D];
       } /* endfor */
    } /* endif */
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[0], 
                          INPUT_PARAMETER_LENGTH_EULER2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[1], 
                          INPUT_PARAMETER_LENGTH_EULER2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[2], 
                          INPUT_PARAMETER_LENGTH_EULER2D, 
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
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Interface_Mesh_Refinements),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Bounding_Box_Mesh_Refinements),
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

    MPI::COMM_WORLD.Bcast(&(IP.IterationsOfInteriorNodesDisturbances),
			  1,
			  MPI::INT,0);

    // Interface input parameters:
    IP.Interface_IP.Broadcast_Input_Parameters();
    MPI::COMM_WORLD.Bcast(&(IP.Reset_Interface_Motion_Type),
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
                          INPUT_PARAMETER_LENGTH_EULER2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Grid_File_Name, 
                          INPUT_PARAMETER_LENGTH_EULER2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Grid_Definition_File_Name, 
                          INPUT_PARAMETER_LENGTH_EULER2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Restart_File_Name, 
                          INPUT_PARAMETER_LENGTH_EULER2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Gnuplot_File_Name, 
                          INPUT_PARAMETER_LENGTH_EULER2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Output_Format_Type, 
                          INPUT_PARAMETER_LENGTH_EULER2D, 
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

    // Freeze_Limiter_Residual_Level
    MPI::COMM_WORLD.Bcast(&(IP.Freeze_Limiter_Residual_Level),
                          1, 
                          MPI::DOUBLE, 0);

    // Accuracy assessment parameters:
    MPI::COMM_WORLD.Bcast(&(IP.Accuracy_Assessment_Mode), 
			  1, 
			  MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Accuracy_Assessment_Exact_Digits), 
			  1, 
			  MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Accuracy_Assessment_Parameter), 
			  1, 
			  MPI::INT, 0);

    // CENO_Execution_Mode variables
    CENO_Execution_Mode::Broadcast();
    
    // CENO_Tolerances variables
    CENO_Tolerances::Broadcast();

    // Exact solution variables
    IP.ExactSoln->Broadcast();

    // HO_Grid2D_Execution_Mode variables
    HO_Grid2D_Execution_Mode::Broadcast();

    // Tecplot_Execution_Mode variables
    Tecplot_Execution_Mode::Broadcast();

    // HighOrder2D_Input variables
    HighOrder2D_Input::Broadcast();    

    // Update all dependent variables
    if (!CFFC_Primary_MPI_Processor()) {

      // Set reference state in the Euler2D_Quad_Block class
      Euler2D_Quad_Block::Set_Normalization_Reference_State(IP.RefU);
    }

#endif

}

#ifdef _MPI_VERSION
/******************************************************//**
 * Routine: Broadcast_Input_Parameters                  
 *                                                      
 * Broadcast the input parameters variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 *routine.                                              
 *                                                      
 ********************************************************/
void Broadcast_Input_Parameters(Euler2D_Input_Parameters &IP,
                                MPI::Intracomm &Communicator, 
                                const int Source_CPU) {

    int Source_Rank = 0;
    int i;

    // CFFC path:
    Communicator.Bcast(IP.CFFC_Path, 
 		       INPUT_PARAMETER_LENGTH_EULER2D, 
		       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Input_File_Name, 
                       INPUT_PARAMETER_LENGTH_EULER2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.Line_Number), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Time_Integration_Type, 
                       INPUT_PARAMETER_LENGTH_EULER2D, 
                       MPI::CHAR, Source_Rank);
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
    Communicator.Bcast(&(IP.Mr_Min_Factor), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.i_Residual_Variable),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Residual_Smoothing), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Residual_Smoothing_Epsilon), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Residual_Smoothing_Gauss_Seidel_Iterations), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Reconstruction_Type, 
                       INPUT_PARAMETER_LENGTH_EULER2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Reconstruction), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.i_ReconstructionMethod), 
		       1, 
		       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Space_Accuracy), 
		       1, 
		       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.IncludeHighOrderBoundariesRepresentation), 
		       1, 
		       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Limiter_Type, 
                       INPUT_PARAMETER_LENGTH_EULER2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Limiter), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Flux_Function_Type, 
                       INPUT_PARAMETER_LENGTH_EULER2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Flux_Function), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.ICs_Type, 
                       INPUT_PARAMETER_LENGTH_EULER2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Gas_Type, 
                       INPUT_PARAMETER_LENGTH_EULER2D, 
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
    Communicator.Bcast(&(IP.Flow_Angle), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
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
       IP.Wo.setgas(IP.Gas_Type);
       IP.Wo = Euler2D_pState(IP.Pressure/(IP.Wo.R*IP.Temperature), 
                              ZERO, 
                              ZERO, 
                              IP.Pressure);
       IP.Wo.v.x = IP.Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Flow_Angle/360.00);
       IP.Wo.v.y = IP.Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Flow_Angle/360.00);
       IP.Uo.setgas(IP.Gas_Type);
       IP.Uo = U(IP.Wo);
       IP.Wo.Mr_min = IP.Mr_Min_Factor*IP.Mach_Number;
       IP.Uo.Mr_min = IP.Mr_Min_Factor*IP.Mach_Number;
    } /* endif */

    // Reference state
    Communicator.Bcast(&(IP.W1.d), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W1.v.x), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W1.v.y), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W1.p), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W2.d), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W2.v.x), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W2.v.y), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W2.p), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.RefU.d), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.RefU.v.x), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.RefU.v.y), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.RefU.p), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Exact_Integration_Digits), 
		       1, 
		       MPI::DOUBLE, Source_Rank);

    Communicator.Bcast(IP.Flow_Geometry_Type, 
                       INPUT_PARAMETER_LENGTH_EULER2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.Axisymmetric), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Grid_Type, 
                       INPUT_PARAMETER_LENGTH_EULER2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.NACA_Aerofoil_Type, 
                       INPUT_PARAMETER_LENGTH_EULER2D, 
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
    Communicator.Bcast(&(IP.Cylinder_Radius2), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Annulus_Theta_Start), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Annulus_Theta_End), 
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
    Communicator.Bcast(&(IP.VertexSW.x), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.VertexSW.y), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.VertexSE.x), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.VertexSE.y), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.VertexNE.x), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.VertexNE.y), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.VertexNW.x), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.VertexNW.y), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Inner_Streamline_Number),
		       1,
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Outer_Streamline_Number),
		       1,
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Isotach_Line),
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
		       INPUT_PARAMETER_LENGTH_EULER2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(&(IP.BCs_Specified),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(IP.BC_North_Type,
		       INPUT_PARAMETER_LENGTH_EULER2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(IP.BC_South_Type,
		       INPUT_PARAMETER_LENGTH_EULER2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(IP.BC_East_Type,
		       INPUT_PARAMETER_LENGTH_EULER2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(IP.BC_West_Type,
		       INPUT_PARAMETER_LENGTH_EULER2D,
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
    Communicator.Bcast(&(IP.Ref_State_BC_North.d),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Ref_State_BC_North.v.x),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Ref_State_BC_North.v.y),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Ref_State_BC_North.p),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Ref_State_BC_South.d),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Ref_State_BC_South.v.x),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Ref_State_BC_South.v.y),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Ref_State_BC_South.p),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Ref_State_BC_East.d),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Ref_State_BC_East.v.x),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Ref_State_BC_East.v.y),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Ref_State_BC_East.p),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Ref_State_BC_West.d),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Ref_State_BC_West.v.x),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Ref_State_BC_West.v.y),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Ref_State_BC_West.p),
		       1,
		       MPI::DOUBLE,Source_Rank);


    // NASA rotors:
    Communicator.Bcast(IP.NASA_Rotor37_Data_Directory, 
                       INPUT_PARAMETER_LENGTH_EULER2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.NASA_Rotor67_Data_Directory, 
                       INPUT_PARAMETER_LENGTH_EULER2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.Rotor_Flow_Type), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Rotor_Percent_Span), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
       IP.ICEMCFD_FileNames = new char*[3];
       for (i = 0; i < 3; i++) {
          IP.ICEMCFD_FileNames[i] = new char[INPUT_PARAMETER_LENGTH_EULER2D];
       } /* endfor */
    } /* endif */
    Communicator.Bcast(IP.ICEMCFD_FileNames[0], 
                       INPUT_PARAMETER_LENGTH_EULER2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.ICEMCFD_FileNames[1], 
                       INPUT_PARAMETER_LENGTH_EULER2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.ICEMCFD_FileNames[2], 
                       INPUT_PARAMETER_LENGTH_EULER2D, 
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

    Communicator.Bcast(&(IP.IterationsOfInteriorNodesDisturbances),
		       1,
		       MPI::INT,Source_Rank);

    // Interface input parameters:
    IP.Interface_IP.Broadcast_Input_Parameters(Communicator,
					       Source_CPU);
    Communicator.Bcast(&(IP.Reset_Interface_Motion_Type),
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
                       INPUT_PARAMETER_LENGTH_EULER2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Grid_File_Name, 
                       INPUT_PARAMETER_LENGTH_EULER2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Grid_Definition_File_Name, 
                       INPUT_PARAMETER_LENGTH_EULER2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Restart_File_Name, 
                       INPUT_PARAMETER_LENGTH_EULER2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Gnuplot_File_Name, 
                       INPUT_PARAMETER_LENGTH_EULER2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Output_Format_Type, 
                       INPUT_PARAMETER_LENGTH_EULER2D, 
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
					       Source_CPU);

    // NKS Parameters
    IP.NKS_IP.Broadcast_Input_Parameters(Communicator, Source_CPU);

    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
       IP.Number_of_Processors = CFFC_MPI::Number_of_Processors;
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
    // Accuracy assessment parameters:
    Communicator.Bcast(&(IP.Accuracy_Assessment_Mode), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Accuracy_Assessment_Exact_Digits), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Accuracy_Assessment_Parameter), 
                       1, 
                       MPI::INT, Source_Rank);

    // Update all dependent variables
    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {    
      // Set reference state in the Euler2D_Quad_Block class
      Euler2D_Quad_Block::Set_Normalization_Reference_State(IP.RefU);
    }

}
#endif

/******************************************************//**
 * Routine: Get_Next_Input_Control_Parameter            
 *                                                      
 * Get the next input control parameter from the input  
 * file.                                                
 *                                                      
 ********************************************************/
void Get_Next_Input_Control_Parameter(Euler2D_Input_Parameters &IP) {

  return IP.Get_Next_Input_Control_Parameter();
}

/******************************************************//**
 * Routine: Parse_Next_Input_Control_Parameter          
 *                                                      
 * Parses and executes the next input control parameter 
 * from the input file.                                 
 *                                                      
 ********************************************************/
int Parse_Next_Input_Control_Parameter(Euler2D_Input_Parameters &IP) {

    int i_command;
    char buffer[256];
    int tpt, bct;
    Vector2D Xt;

    i_command = 0;

    if (strcmp(IP.Next_Control_Parameter, "CFFC_Path") == 0) {
       i_command = 1111;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.CFFC_Path, IP.Next_Control_Parameter);
       strcpy(IP.NASA_Rotor37_Data_Directory, IP.CFFC_Path);
       strcpy(IP.NASA_Rotor37_Data_Directory, "/data/NASA_Rotors/R37/");
       strcpy(IP.NASA_Rotor67_Data_Directory, IP.CFFC_Path);
       strcpy(IP.NASA_Rotor67_Data_Directory, "/data/NASA_Rotors/R67/");

    } else if (strcmp(IP.Next_Control_Parameter, "Time_Integration_Type") == 0) {
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
       } else if (strcmp(IP.Time_Integration_Type, "Multigrid") == 0) {
	  IP.i_Time_Integration = TIME_STEPPING_MULTIGRID;
       } else if (strcmp(IP.Time_Integration_Type, "Dual_Time_Stepping") == 0) {
	  IP.i_Time_Integration = TIME_STEPPING_DUAL_TIME_STEPPING;
	  IP.Multigrid_IP.i_Dual_Time_Stepping = ON;
       } else {
	 std::cout << "\n ==> Unknown time marching scheme!";
	 i_command = INVALID_INPUT_VALUE;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Reconstruction_Type") == 0) {
       i_command = 2;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Reconstruction_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Reconstruction_Type, "Green_Gauss") == 0) {
          IP.i_Reconstruction = RECONSTRUCTION_GREEN_GAUSS;
	  IP.i_ReconstructionMethod = RECONSTRUCTION_GREEN_GAUSS;
	  CENO_Execution_Mode::USE_CENO_ALGORITHM = OFF;
       } else if (strcmp(IP.Reconstruction_Type, "Least_Squares") == 0) {
          IP.i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;
	  IP.i_ReconstructionMethod = RECONSTRUCTION_LEAST_SQUARES;
	  CENO_Execution_Mode::USE_CENO_ALGORITHM = OFF;
       } else if (strcmp(IP.Reconstruction_Type, "CENO") == 0) {
	 IP.i_Reconstruction = RECONSTRUCTION_HIGH_ORDER;
	 IP.i_ReconstructionMethod = RECONSTRUCTION_CENO;
	 CENO_Execution_Mode::USE_CENO_ALGORITHM = ON;
       } else {
	 std::cout << "\n ==> Unknown reconstruction method!";
	 i_command = INVALID_INPUT_VALUE;
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
       } else if (strcmp(IP.Limiter_Type, "Venkatakrishnan_Modified") == 0) {
	 IP.i_Limiter = LIMITER_VENKATAKRISHNAN_CORRECTED;
       } else {
	 std::cout << "\n ==> Unknown limiter type!";
	 i_command = INVALID_INPUT_VALUE;
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
       } else if (strcmp(IP.Flux_Function_Type, "VanLeer") == 0) {
	  IP.i_Flux_Function = FLUX_FUNCTION_VANLEER;
       } else if (strcmp(IP.Flux_Function_Type, "AUSM") == 0) {
	  IP.i_Flux_Function = FLUX_FUNCTION_AUSM;
       } else if (strcmp(IP.Flux_Function_Type, "AUSM+") == 0) {
	  IP.i_Flux_Function = FLUX_FUNCTION_AUSMplus;
       } else if (strcmp(IP.Flux_Function_Type, "Roe_Precon_WS") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_ROE_PRECON_WS;
       } else if (strcmp(IP.Flux_Function_Type, "HLLE_Precon_WS") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_HLLE_PRECON_WS;
       } else if (strcmp(IP.Flux_Function_Type, "Godunov_MB") == 0) {
	  IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV_MB;
       } else if (strcmp(IP.Flux_Function_Type, "Roe_MB") == 0) {
	  IP.i_Flux_Function = FLUX_FUNCTION_ROE_MB;
       } else if (strcmp(IP.Flux_Function_Type, "HLLE_MB") == 0) {
	  IP.i_Flux_Function = FLUX_FUNCTION_HLLE_MB;
       } else if (strcmp(IP.Flux_Function_Type, "VanLeer_MB") == 0) {
	  IP.i_Flux_Function = FLUX_FUNCTION_VANLEER_MB;
       } else {
	 std::cout << "\n ==> Unknown flux function type!";
	  i_command = INVALID_INPUT_VALUE;
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
       } else if (strcmp(IP.ICs_Type,"Ringleb_Flow") == 0) {
          IP.i_ICs = IC_RINGLEB_FLOW;
       } else if (strcmp(IP.ICs_Type,"Cylindrical_Explosion") == 0) {
	  IP.i_ICs = IC_CYLINDRICAL_EXPLOSION;
       } else if (strcmp(IP.ICs_Type,"Cylindrical_Implosion") == 0) {
          IP.i_ICs = IC_CYLINDRICAL_IMPLOSION;
       } else if (strcmp(IP.ICs_Type, "Wedge_Flow") == 0) {
	  IP.i_ICs = IC_WEDGE_FLOW;
       } else if (strcmp(IP.ICs_Type, "Unsteady_Blunt_Body") == 0) {
          IP.i_ICs = IC_UNSTEADY_BLUNT_BODY;
       } else if (strcmp(IP.ICs_Type,"Sin_X") == 0) {
          IP.i_ICs = IC_PERIODIC_SINX_WAVE;
       } else if (strcmp(IP.ICs_Type,"Sin_Y") == 0) {
          IP.i_ICs = IC_PERIODIC_SINY_WAVE;
       } else if (strcmp(IP.ICs_Type,"Multiple_Sin_X") == 0) {
          IP.i_ICs = IC_PERIODIC_SINX_MULTIWAVE;
       } else if (strcmp(IP.ICs_Type,"Multiple_Sin_Y") == 0) {
          IP.i_ICs = IC_PERIODIC_SINY_MULTIWAVE;
       } else if (strcmp(IP.ICs_Type,"Complex_Waves") == 0) {
          IP.i_ICs = IC_PERIODIC_COMPLEX_MULTIWAVE;
       } else if (strcmp(IP.ICs_Type,"2D_Abgrall") == 0) {
          IP.i_ICs = IC_ABGRALL_FUNCTION;
       } else if (strcmp(IP.ICs_Type,"SinExp_X") == 0) {
          IP.i_ICs = IC_SIN_EXP_X_WAVE;
       } else if (strcmp(IP.ICs_Type,"SinExp_Y") == 0) {
          IP.i_ICs = IC_SIN_EXP_Y_WAVE;
       } else if (strcmp(IP.ICs_Type,"SinExp_Rotated") == 0) {
          IP.i_ICs = IC_SIN_EXP_ROTATED_WAVE;
       } else if (strcmp(IP.ICs_Type,"Cosine_Hill") == 0) {
          IP.i_ICs = IC_COSINE_HILL;
       } else if (strcmp(IP.ICs_Type,"Hyper_Tangent") == 0) {
          IP.i_ICs = IC_HYPER_TANGENT;
       } else if (strcmp(IP.ICs_Type, "Blast_Wave_IVP") == 0) {
	  IP.i_ICs = IC_BLAST_WAVE_INTERACTION;
       } else if (strcmp(IP.ICs_Type, "Acoustic_Shock_IVP") == 0) {
	  IP.i_ICs = IC_SHOCK_ACOUSTIC_INTERACTION;
       } else if (strcmp(IP.ICs_Type,"Exact_Solution") == 0) {
	 IP.i_ICs = IC_EXACT_SOLUTION;
       } else if (strcmp(IP.ICs_Type,"Uniform_Interior_Exact_Ghost_Cells") == 0) {
	 IP.i_ICs = IC_INTERIOR_UNIFORM_GHOSTCELLS_EXACT;
       } else if (strcmp(IP.ICs_Type, "Restart") == 0) {
	 IP.i_ICs = IC_RESTART;
       } else {
	 std::cout << "\n ==> Unknown initial condition!";
	 i_command = INVALID_INPUT_VALUE;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Exact_Integration_Digits") == 0) {
      i_command = 0;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Exact_Integration_Digits;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Exact_Integration_Digits < 0) i_command = INVALID_INPUT_VALUE;

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
       } else if (strcmp(IP.Grid_Type, "Deformed_Box") == 0) {
	 IP.i_Grid = GRID_DEFORMED_BOX;
       } else if (strcmp(IP.Grid_Type, "Periodic_Box") == 0) {
	 IP.i_Grid = GRID_PERIODIC_BOX;
       } else if (strcmp(IP.Grid_Type, "Flat_Plate") == 0) {
          IP.i_Grid = GRID_FLAT_PLATE;
          IP.Plate_Length = ONE;
	  IP.BC_South = BC_REFLECTION;
       } else if (strcmp(IP.Grid_Type, "Pipe") == 0) {
          IP.i_Grid = GRID_PIPE;
          IP.Pipe_Length = ONE;
          IP.Pipe_Radius = HALF;
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
	  IP.BC_North = BC_REFLECTION;
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
       } else if (strcmp(IP.Grid_Type, "Circular_Cylinder") == 0) {
          IP.i_Grid = GRID_CIRCULAR_CYLINDER;
          IP.Cylinder_Radius = ONE;
          IP.Cylinder_Radius2 = 32.00;
	  IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
	  IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
	  IP.Mesh_Stretching_Factor_Idir = 1.025;
	  IP.Mesh_Stretching_Factor_Jdir = 1.001;
	  IP.BC_South = BC_REFLECTION;
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
	  IP.BC_South = BC_REFLECTION;
       } else if (strcmp(IP.Grid_Type, "Unsteady_Blunt_Body") == 0) {
	  IP.i_Grid = GRID_UNSTEADY_BLUNT_BODY;
	  IP.Blunt_Body_Radius = ONE;
	  IP.Blunt_Body_Mach_Number = TWO;
       } else if (strcmp(IP.Grid_Type, "NASA_Rotor37") == 0) {
          IP.i_Grid = GRID_NASA_ROTOR_37;
          IP.Rotor_Flow_Type = PEAK_FLOW;
          IP.Rotor_Percent_Span = 50.00;
          IP.NASA_Rotor37.init(IP.Rotor_Flow_Type,
                               IP.NASA_Rotor37_Data_Directory);
          IP.Wo.setgas("AIR");
	  IP.NASA_Rotor37.getPstateREL_up(IP.Wo,IP.Rotor_Percent_Span);
          IP.Pressure = IP.Wo.p;
          IP.Temperature = IP.Wo.T();
          IP.Mach_Number = IP.NASA_Rotor37.getMachREL_up(IP.Rotor_Percent_Span);
          IP.Flow_Angle = atan2(IP.Wo.v.y, IP.Wo.v.x); 
          if (IP.Flow_Angle < ZERO) IP.Flow_Angle = TWO*PI + IP.Flow_Angle;
          IP.Flow_Angle = 180.00*IP.Flow_Angle/PI;
	  IP.NASA_Rotor37.getPstateREL_down(IP.W1,IP.Rotor_Percent_Span);
          IP.Wo.Mr_min = IP.Mr_Min_Factor*IP.Mach_Number;
          IP.Uo.Mr_min = IP.Mr_Min_Factor*IP.Mach_Number;
       } else if (strcmp(IP.Grid_Type, "NASA_Rotor67") == 0) {
          IP.i_Grid = GRID_NASA_ROTOR_67;
          IP.Rotor_Flow_Type = PEAK_FLOW;
          IP.Rotor_Percent_Span = 50.00;
          IP.NASA_Rotor67.init(IP.Rotor_Flow_Type,
                               IP.NASA_Rotor67_Data_Directory);
          IP.Wo.setgas("AIR");
          IP.NASA_Rotor67.getPstateREL_up(IP.Wo,IP.Rotor_Percent_Span);
          IP.Pressure = IP.Wo.p;
          IP.Temperature = IP.Wo.T();
          IP.Mach_Number = IP.NASA_Rotor67.getMachREL_up(IP.Rotor_Percent_Span);
          IP.Flow_Angle = atan2(IP.Wo.v.y, IP.Wo.v.x); 
          if (IP.Flow_Angle < ZERO) IP.Flow_Angle = TWO*PI + IP.Flow_Angle;
          IP.Flow_Angle = 180.00*IP.Flow_Angle/PI;
	  IP.NASA_Rotor67.getPstateREL_down(IP.W1,IP.Rotor_Percent_Span);
          IP.Wo.Mr_min = IP.Mr_Min_Factor*IP.Mach_Number;
          IP.Uo.Mr_min = IP.Mr_Min_Factor*IP.Mach_Number;
       } else if (strcmp(IP.Grid_Type,"Ringleb_Flow") == 0) {
          IP.i_Grid = GRID_RINGLEB_FLOW;
	  IP.Inner_Streamline_Number = 0.80;
	  IP.Outer_Streamline_Number = 0.40;
	  IP.Isotach_Line = 0.30;
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
	 std::cout << "\n ==> Unknown grid type!";
	 i_command = INVALID_INPUT_VALUE;
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
       } else {
	 IP.Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Local_Time_Stepping") == 0) {
       i_command = 15;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Local_Time_Stepping;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Local_Time_Stepping != GLOBAL_TIME_STEPPING &&
           IP.Local_Time_Stepping != SCALAR_LOCAL_TIME_STEPPING &&
	   IP.Local_Time_Stepping != MATRIX_LOCAL_TIME_STEPPING &&
           IP.Local_Time_Stepping != LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER)
         IP.Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;

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

    } else if (strcmp(IP.Next_Control_Parameter, "Cylinder_Radius2") == 0) {
       i_command = 26;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Cylinder_Radius2;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Cylinder_Radius2 <= ZERO) i_command = INVALID_INPUT_VALUE;
       if (IP.Cylinder_Radius2 <= IP.Cylinder_Radius) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Annulus_Start_Angle") == 0) {
      i_command = 26;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Annulus_Theta_Start;
      IP.Input_File.getline(buffer, sizeof(buffer));
      
    } else if (strcmp(IP.Next_Control_Parameter, "Annulus_End_Angle") == 0) {
      i_command = 26;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Annulus_Theta_End;
      IP.Input_File.getline(buffer, sizeof(buffer));
      
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

    } else if (strcmp(IP.Next_Control_Parameter,"Inner_Streamline_Number") == 0) {
       i_command = 31;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Inner_Streamline_Number;
       IP.Input_File.getline(buffer,sizeof(buffer));
       if (IP.Inner_Streamline_Number <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Outer_Streamline_Number") == 0) {
       i_command = 31;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Outer_Streamline_Number;
       IP.Input_File.getline(buffer,sizeof(buffer));
       if (IP.Outer_Streamline_Number <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Isotach_Line") == 0) {
       i_command = 31;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Isotach_Line;
       IP.Input_File.getline(buffer,sizeof(buffer));
       if (IP.Isotach_Line <= ZERO) i_command = INVALID_INPUT_VALUE;

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

    } else if (strcmp(IP.Next_Control_Parameter, "VertexSW") == 0) {
      i_command = 0;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.VertexSW;
      IP.Input_File.setf(ios::skipws);
      IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "VertexSE") == 0) {
      i_command = 0;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.VertexSE;
      IP.Input_File.setf(ios::skipws);
      IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "VertexNE") == 0) {
      i_command = 0;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.VertexNE;
      IP.Input_File.setf(ios::skipws);
      IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "VertexNW") == 0) {
      i_command = 0;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.VertexNW;
      IP.Input_File.setf(ios::skipws);
      IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Gas_Type") == 0) {
       i_command = 41;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Gas_Type, 
              IP.Next_Control_Parameter);
       IP.Wo.setgas(IP.Gas_Type);
       IP.Wo = Euler2D_pState(IP.Pressure/(IP.Wo.R*IP.Temperature), 
                              ZERO, 
                              ZERO, 
                              IP.Pressure);
       IP.Wo.v.x = IP.Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Flow_Angle/360.00);
       IP.Wo.v.y = IP.Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Flow_Angle/360.00);
       IP.Uo.setgas(IP.Gas_Type);
       IP.Uo = U(IP.Wo);

    } else if (strcmp(IP.Next_Control_Parameter, "Mach_Number") == 0) {
       i_command = 42;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Mach_Number;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Mach_Number < ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } else {
          IP.Wo.v.x = IP.Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Flow_Angle/360.00);
          IP.Wo.v.y = IP.Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Flow_Angle/360.00);
       } /* endif */
       IP.Wo.Mr_min = IP.Mr_Min_Factor*IP.Mach_Number;
       IP.Uo.Mr_min = IP.Mr_Min_Factor*IP.Mach_Number;

    } else if (strcmp(IP.Next_Control_Parameter, "Flow_Angle") == 0) {
       i_command = 43;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Flow_Angle;
       IP.Input_File.getline(buffer, sizeof(buffer));
       IP.Wo.v.x = IP.Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Flow_Angle/360.00);
       IP.Wo.v.y = IP.Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Flow_Angle/360.00);

    } else if (strcmp(IP.Next_Control_Parameter, "Pressure") == 0) {
       i_command = 44;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Pressure;
       IP.Input_File.getline(buffer, sizeof(buffer));
       IP.Pressure = IP.Pressure*THOUSAND;
       if (IP.Pressure <= ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } else {
          IP.Wo = Euler2D_pState(IP.Pressure/(IP.Wo.R*IP.Temperature), 
                                 ZERO, 
                                 ZERO, 
                                 IP.Pressure);
          IP.Wo.v.x = IP.Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Flow_Angle/360.00);
          IP.Wo.v.y = IP.Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Flow_Angle/360.00);
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Temperature") == 0) {
       i_command = 45;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Temperature;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Temperature <= ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } else {
          IP.Wo = Euler2D_pState(IP.Pressure/(IP.Wo.R*IP.Temperature), 
                                 ZERO, 
                                 ZERO, 
                                 IP.Pressure);
          IP.Wo.v.x = IP.Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Flow_Angle/360.00);
          IP.Wo.v.y = IP.Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Flow_Angle/360.00);
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter,"Wave_Position") == 0) {
       i_command = 46;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Wave_Position.x >> IP.Wave_Position.y;
       IP.Input_File.getline(buffer,sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter,"Wave_Width") == 0) {
       i_command = 47;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Wave_Width;
       IP.Input_File.getline(buffer,sizeof(buffer));
       if (IP.Wave_Width < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Time_Max") == 0) {
       i_command = 48;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Time_Max;
       IP.Input_File.getline(buffer, sizeof(buffer));
       IP.Time_Max = IP.Time_Max/THOUSAND;
       if (IP.Time_Max < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Blocks_Per_Processor") == 0) {
       i_command = 49;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Blocks_Per_Processor;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Blocks_Per_Processor < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Output_Format_Type") == 0) {
       i_command = 50;
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

    } else if (strcmp(IP.Next_Control_Parameter, "NASA_Rotor37_Data_Directory") == 0) {
       i_command = 51;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.NASA_Rotor37_Data_Directory, IP.Next_Control_Parameter);

    } else if (strcmp(IP.Next_Control_Parameter, "NASA_Rotor67_Data_Directory") == 0) {
       i_command = 52;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.NASA_Rotor67_Data_Directory, IP.Next_Control_Parameter);

    } else if (strcmp(IP.Next_Control_Parameter, "Rotor_Flow_Type") == 0) {
       i_command = 53;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Rotor_Flow_Type;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Rotor_Flow_Type != PEAK_FLOW &&
           IP.Rotor_Flow_Type != STALL_FLOW) IP.Rotor_Flow_Type = PEAK_FLOW;
       if (IP.i_Grid == GRID_NASA_ROTOR_37) {
          IP.NASA_Rotor37.init(IP.Rotor_Flow_Type,
                               IP.NASA_Rotor37_Data_Directory);
       } else if (IP.i_Grid == GRID_NASA_ROTOR_67) {
          IP.NASA_Rotor67.init(IP.Rotor_Flow_Type,
                               IP.NASA_Rotor67_Data_Directory);
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Rotor_Percent_Span") == 0) {
       i_command = 54;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Rotor_Percent_Span;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Rotor_Percent_Span < ZERO ||
           IP.Rotor_Percent_Span > HUNDRED) {
          i_command = INVALID_INPUT_VALUE;
       } else if (IP.i_Grid == GRID_NASA_ROTOR_37) {
          IP.Wo.setgas("AIR");
	  IP.NASA_Rotor37.getPstateREL_up(IP.Wo,IP.Rotor_Percent_Span);
          IP.Pressure = IP.Wo.p;
          IP.Temperature = IP.Wo.T();
          IP.Mach_Number = IP.NASA_Rotor37.getMachREL_up(IP.Rotor_Percent_Span);
          IP.Flow_Angle = atan2(IP.Wo.v.y, IP.Wo.v.x); 
          if (IP.Flow_Angle < ZERO) IP.Flow_Angle = TWO*PI + IP.Flow_Angle;
          IP.Flow_Angle = 180.00*IP.Flow_Angle/PI;
	  IP.NASA_Rotor37.getPstateREL_down(IP.W1,IP.Rotor_Percent_Span);
          IP.Wo.Mr_min = IP.Mr_Min_Factor*IP.Mach_Number;
          IP.Uo.Mr_min = IP.Mr_Min_Factor*IP.Mach_Number;
       } else if (IP.i_Grid == GRID_NASA_ROTOR_67) {
          IP.Wo.setgas("AIR");
          IP.NASA_Rotor67.getPstateREL_up(IP.Wo,IP.Rotor_Percent_Span);
          IP.Pressure = IP.Wo.p;
          IP.Temperature = IP.Wo.T();
          IP.Mach_Number = IP.NASA_Rotor67.getMachREL_up(IP.Rotor_Percent_Span);
          IP.Flow_Angle = atan2(IP.Wo.v.y, IP.Wo.v.x); 
          if (IP.Flow_Angle < ZERO) IP.Flow_Angle = TWO*PI + IP.Flow_Angle;
          IP.Flow_Angle = 180.00*IP.Flow_Angle/PI;
          IP.NASA_Rotor67.getPstateREL_down(IP.W1,IP.Rotor_Percent_Span);
          IP.Wo.Mr_min = IP.Mr_Min_Factor*IP.Mach_Number;
          IP.Uo.Mr_min = IP.Mr_Min_Factor*IP.Mach_Number;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Flow_Geometry_Type") == 0) {
       i_command = 55;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Flow_Geometry_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Flow_Geometry_Type, "Planar") == 0) {
          IP.Axisymmetric = 0;
       } else if (strcmp(IP.Flow_Geometry_Type, "Axisymmetric") == 0) {
          IP.Axisymmetric = 1;
       } else {
          IP.Axisymmetric = 0;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Restart_Solution_Save_Frequency") == 0) {
       i_command = 56;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Restart_Solution_Save_Frequency;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Restart_Solution_Save_Frequency < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Output_Progress_Frequency") == 0) {
       i_command = 57;
       IP.Line_Number++;
       IP.Input_File >> IP.Output_Progress_Frequency;
       IP.Input_File.getline(buffer,sizeof(buffer));
       if (IP.Output_Progress_Frequency < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "ICEMCFD_Topology_File") == 0) {
       i_command = 58;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.ICEMCFD_FileNames[0], IP.Next_Control_Parameter);

    } else if (strcmp(IP.Next_Control_Parameter, "ICEMCFD_Family_Boco_File") == 0) {
       i_command = 59;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.ICEMCFD_FileNames[1], IP.Next_Control_Parameter);

    } else if (strcmp(IP.Next_Control_Parameter, "ICEMCFD_Family_Topo_File") == 0) {
       i_command = 60;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.ICEMCFD_FileNames[2], IP.Next_Control_Parameter);

    } else if (strcmp(IP.Next_Control_Parameter, "X_Shift") == 0) {
       i_command = 61;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.X_Shift;
       IP.Input_File.setf(ios::skipws);
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "X_Scale") == 0) {
       i_command = 62;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.X_Scale;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "X_Rotate") == 0) {
       i_command = 63;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.X_Rotate;
       IP.Input_File.getline(buffer, sizeof(buffer));
    } else if (strcmp(IP.Next_Control_Parameter, "Freeze_Limiter") == 0) {
      // Freeze_Limiter:
      i_command = 74;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Freeze_Limiter;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Freeze_Limiter < 0) i_command = INVALID_INPUT_VALUE;
      
    } else if (strcmp(IP.Next_Control_Parameter, "Freeze_Limiter_Residual_Level") == 0) {
      // Freeze_Limiter_Residual_Level:
      i_command = 75;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Freeze_Limiter_Residual_Level;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Freeze_Limiter_Residual_Level < ZERO) i_command = INVALID_INPUT_VALUE;
      
    } else if (strcmp(IP.Next_Control_Parameter, "AMR") == 0) {
      i_command = 77;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.AMR = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.AMR = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter, "AMR_Frequency") == 0) {
      i_command = 78;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.AMR_Frequency;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.AMR_Frequency < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Initial_Mesh_Refinements") == 0) {
      i_command = 79;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Number_of_Initial_Mesh_Refinements;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Number_of_Initial_Mesh_Refinements < 0) IP.Number_of_Initial_Mesh_Refinements = 0;

    } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Uniform_Mesh_Refinements") == 0) {
      i_command = 80;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Number_of_Uniform_Mesh_Refinements;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Number_of_Uniform_Mesh_Refinements < 0) IP.Number_of_Uniform_Mesh_Refinements = 0;

    } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Boundary_Mesh_Refinements") == 0) {
      i_command = 81;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Number_of_Boundary_Mesh_Refinements;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Number_of_Boundary_Mesh_Refinements < 0) IP.Number_of_Boundary_Mesh_Refinements = 0;

    } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Interface_Mesh_Refinements") == 0) {
      i_command = 82;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Number_of_Interface_Mesh_Refinements;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Number_of_Interface_Mesh_Refinements < 0) IP.Number_of_Interface_Mesh_Refinements = 0;

    } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Bounding_Box_Mesh_Refinements") == 0) {
      i_command = 83;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Number_of_Bounding_Box_Mesh_Refinements;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Number_of_Bounding_Box_Mesh_Refinements < 0) IP.Number_of_Bounding_Box_Mesh_Refinements = 0;

    } else if (strcmp(IP.Next_Control_Parameter,"Interface_Refinement_Condition") == 0) {
      i_command = 84;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Interface_Refinement_Condition = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.Interface_Refinement_Condition = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Maximum_Refinement_Level") == 0) {
      i_command = 85;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Maximum_Refinement_Level;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Maximum_Refinement_Level < 1) IP.Maximum_Refinement_Level = 1;

    } else if (strcmp(IP.Next_Control_Parameter,"Minimum_Refinement_Level") == 0) {
      i_command = 86;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Minimum_Refinement_Level;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Minimum_Refinement_Level < 1) IP.Minimum_Refinement_Level = 1;

    } else if (strcmp(IP.Next_Control_Parameter, "Threshold_for_Refinement") == 0) {
       i_command = 87;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Threshold_for_Refinement;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Threshold_for_Refinement <= ZERO ||
           IP.Threshold_for_Refinement > ONE) IP.Threshold_for_Refinement = 0.50;

    } else if (strcmp(IP.Next_Control_Parameter, "Threshold_for_Coarsening") == 0) {
       i_command = 88;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Threshold_for_Coarsening;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Threshold_for_Coarsening < ZERO ||
           IP.Threshold_for_Coarsening >= ONE) IP.Threshold_for_Coarsening = 0.10;

    } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Refinement_Criteria") == 0) {
      i_command = 89;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Number_of_Refinement_Criteria;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Number_of_Refinement_Criteria < 1 || IP.Number_of_Refinement_Criteria > 6) {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Gradient_Density") == 0) {
      i_command = 90;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Refinement_Criteria_Gradient_Density = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.Refinement_Criteria_Gradient_Density = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Divergence_Velocity") == 0) {
      i_command = 91;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Refinement_Criteria_Divergence_Velocity = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.Refinement_Criteria_Divergence_Velocity = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Curl_Velocity") == 0) {
      i_command = 92;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Refinement_Criteria_Curl_Velocity = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.Refinement_Criteria_Curl_Velocity = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"AMR_Xmin") == 0) {
      i_command = 93;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.AMR_Xmin;
      IP.Input_File.setf(ios::skipws);
      IP.Input_File.getline(buffer,sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter,"AMR_Xmax") == 0) {
      i_command = 94;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.AMR_Xmax;
      IP.Input_File.setf(ios::skipws);
      IP.Input_File.getline(buffer,sizeof(buffer));
      
    } else if (strcmp(IP.Next_Control_Parameter,"Residual_Variable") == 0) {
      i_command = 95;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.i_Residual_Variable;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.i_Residual_Variable < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Residual_Smoothing_Epsilon") == 0) {
      i_command = 96;
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
      i_command = 97;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Residual_Smoothing_Gauss_Seidel_Iterations;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Residual_Smoothing_Gauss_Seidel_Iterations < 0) {
         IP.Residual_Smoothing_Gauss_Seidel_Iterations = 0;
      } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Mr_Min_Factor") == 0) {
      i_command = 98;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Mr_Min_Factor;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Mr_Min_Factor <= ZERO) {
         IP.Mr_Min_Factor = ZERO;
      } /* endif */
      IP.Wo.Mr_min = IP.Mr_Min_Factor*IP.Mach_Number;
      IP.Uo.Mr_min = IP.Mr_Min_Factor*IP.Mach_Number;

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

    } else if (strcmp(IP.Next_Control_Parameter,"Iteration_Disturb_Mesh") == 0) {
      i_command = 0;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.IterationsOfInteriorNodesDisturbances;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.IterationsOfInteriorNodesDisturbances < 0 ) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Number_Spline_Points") == 0) {
      i_command = 101;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Num_Of_Spline_Control_Points;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Num_Of_Spline_Control_Points <= TWO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Morton") == 0) {
      i_command = 107;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Morton;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Morton < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Morton_Reordering_Frequency") == 0) {
      i_command = 108;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Morton_Reordering_Frequency;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Morton_Reordering_Frequency < 0) i_command = INVALID_INPUT_VALUE;

    /**************************************
     * Multigrid Related Input Parameters *
     **************************************/

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

    /******************************************
     * End of Multigrid related Input parsing *
     ******************************************/

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
      } else if (strcmp(IP.BC_North_Type,"Periodic") == 0) {
	IP.BC_North = BC_PERIODIC;
      } else if (strcmp(IP.BC_North_Type,"Burning_Surface") == 0) {
	IP.BC_North = BC_BURNING_SURFACE;
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
      } else if (strcmp(IP.BC_North_Type,"Characteristic_Velocity") == 0) {
	IP.BC_North = BC_CHARACTERISTIC_VELOCITY;
      } else if (strcmp(IP.BC_North_Type,"None") == 0) {
	IP.BC_North = BC_NONE;
      } else if (strcmp(IP.BC_North_Type,"Ringleb") == 0) {
	IP.BC_North = BC_RINGLEB_FLOW;
      } else if (strcmp(IP.BC_North_Type,"Dirichlet") == 0) {
	IP.BC_North = BC_DIRICHLET;
      } else if (strcmp(IP.BC_North_Type,"Neumann") == 0) {
	IP.BC_North = BC_NEUMANN;
      } else if (strcmp(IP.BC_North_Type,"Robin") == 0) {
	IP.BC_North = BC_ROBIN;
      } else if (strcmp(IP.BC_North_Type,"Farfield") == 0) {
	IP.BC_North = BC_FARFIELD;
      } else if (strcmp(IP.BC_North_Type,"Frozen") == 0) {
	IP.BC_North = BC_FROZEN;
      } else if (strcmp(IP.BC_North_Type,"Exact_Solution") == 0) {
	IP.BC_North = BC_EXACT_SOLUTION;
      } else if (strcmp(IP.BC_North_Type,"Streamline") == 0) {
	IP.BC_North = BC_STREAMLINE;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter, "Ref_State_North") == 0) {
      i_command = 0;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Ref_State_BC_North;
      IP.Input_File.setf(ios::skipws);
      IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter,"BC_South") == 0) {
      i_command = 502;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.BC_South_Type,IP.Next_Control_Parameter);
      if (strcmp(IP.BC_South_Type,"Reflection") == 0) {
	IP.BC_South = BC_REFLECTION;
      } else if (strcmp(IP.BC_South_Type,"Periodic") == 0) {
	IP.BC_South = BC_PERIODIC;
      } else if (strcmp(IP.BC_South_Type,"Burning_Surface") == 0) {
	IP.BC_South = BC_BURNING_SURFACE;
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
      } else if (strcmp(IP.BC_South_Type,"Characteristic_Velocity") == 0) {
	IP.BC_South = BC_CHARACTERISTIC_VELOCITY;
      } else if (strcmp(IP.BC_South_Type,"None") == 0) {
	IP.BC_South = BC_NONE;
      } else if (strcmp(IP.BC_South_Type,"Ringleb") == 0) {
	IP.BC_South = BC_RINGLEB_FLOW;
      } else if (strcmp(IP.BC_South_Type,"Dirichlet") == 0) {
	IP.BC_South = BC_DIRICHLET;
      } else if (strcmp(IP.BC_South_Type,"Neumann") == 0) {
	IP.BC_South = BC_NEUMANN;
      } else if (strcmp(IP.BC_South_Type,"Robin") == 0) {
	IP.BC_South = BC_ROBIN;
      } else if (strcmp(IP.BC_South_Type,"Farfield") == 0) {
	IP.BC_South = BC_FARFIELD;
      } else if (strcmp(IP.BC_South_Type,"Frozen") == 0) {
	IP.BC_South = BC_FROZEN;
      } else if (strcmp(IP.BC_South_Type,"Exact_Solution") == 0) {
	IP.BC_South = BC_EXACT_SOLUTION;
      } else if (strcmp(IP.BC_South_Type,"Streamline") == 0) {
	IP.BC_South = BC_STREAMLINE;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter, "Ref_State_South") == 0) {
      i_command = 0;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Ref_State_BC_South;
      IP.Input_File.setf(ios::skipws);
      IP.Input_File.getline(buffer, sizeof(buffer));
      
    } else if (strcmp(IP.Next_Control_Parameter,"BC_East") == 0) {
      i_command = 503;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.BC_East_Type,IP.Next_Control_Parameter);
      if (strcmp(IP.BC_East_Type,"Reflection") == 0) {
	IP.BC_East = BC_REFLECTION;
      } else if (strcmp(IP.BC_East_Type,"Periodic") == 0) {
	IP.BC_East = BC_PERIODIC;
      } else if (strcmp(IP.BC_East_Type,"Burning_Surface") == 0) {
	IP.BC_East = BC_BURNING_SURFACE;
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
      } else if (strcmp(IP.BC_East_Type,"Characteristic_Velocity") == 0) {
	IP.BC_East = BC_CHARACTERISTIC_VELOCITY;
      } else if (strcmp(IP.BC_East_Type,"None") == 0) {
	IP.BC_East = BC_NONE;
      } else if (strcmp(IP.BC_East_Type,"Ringleb") == 0) {
	IP.BC_East = BC_RINGLEB_FLOW;
      } else if (strcmp(IP.BC_East_Type,"Dirichlet") == 0) {
	IP.BC_East = BC_DIRICHLET;
      } else if (strcmp(IP.BC_East_Type,"Neumann") == 0) {
	IP.BC_East = BC_NEUMANN;
      } else if (strcmp(IP.BC_East_Type,"Robin") == 0) {
	IP.BC_East = BC_ROBIN;
      } else if (strcmp(IP.BC_East_Type,"Farfield") == 0) {
	IP.BC_East = BC_FARFIELD;
      } else if (strcmp(IP.BC_East_Type,"Frozen") == 0) {
	IP.BC_East = BC_FROZEN;
      } else if (strcmp(IP.BC_East_Type,"Exact_Solution") == 0) {
	IP.BC_East = BC_EXACT_SOLUTION;
      } else if (strcmp(IP.BC_East_Type,"Streamline") == 0) {
	IP.BC_East = BC_STREAMLINE;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter, "Ref_State_East") == 0) {
      i_command = 0;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Ref_State_BC_East;
      IP.Input_File.setf(ios::skipws);
      IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter,"BC_West") == 0) {
      i_command = 504;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.BC_West_Type,IP.Next_Control_Parameter);
      if (strcmp(IP.BC_West_Type,"Reflection") == 0) {
	IP.BC_West = BC_REFLECTION;
      } else if (strcmp(IP.BC_West_Type,"Periodic") == 0) {
	IP.BC_West = BC_PERIODIC;
      } else if (strcmp(IP.BC_West_Type,"Burning_Surface") == 0) {
	IP.BC_West = BC_BURNING_SURFACE;
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
      } else if (strcmp(IP.BC_West_Type,"Characteristic_Velocity") == 0) {
	IP.BC_West = BC_CHARACTERISTIC_VELOCITY;
      } else if (strcmp(IP.BC_West_Type,"None") == 0) {
	IP.BC_West = BC_NONE;
      } else if (strcmp(IP.BC_West_Type,"Ringleb") == 0) {
	IP.BC_West = BC_RINGLEB_FLOW;
      } else if (strcmp(IP.BC_West_Type,"Dirichlet") == 0) {
	IP.BC_West = BC_DIRICHLET;
      } else if (strcmp(IP.BC_West_Type,"Neumann") == 0) {
	IP.BC_West = BC_NEUMANN;
      } else if (strcmp(IP.BC_West_Type,"Robin") == 0) {
	IP.BC_West = BC_ROBIN;
      } else if (strcmp(IP.BC_West_Type,"Farfield") == 0) {
	IP.BC_West = BC_FARFIELD;
      } else if (strcmp(IP.BC_West_Type,"Frozen") == 0) {
	IP.BC_West = BC_FROZEN;
      } else if (strcmp(IP.BC_West_Type,"Exact_Solution") == 0) {
	IP.BC_West = BC_EXACT_SOLUTION;
      } else if (strcmp(IP.BC_West_Type,"Streamline") == 0) {
	IP.BC_West = BC_STREAMLINE;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter, "Ref_State_West") == 0) {
      i_command = 0;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Ref_State_BC_West;
      IP.Input_File.setf(ios::skipws);
      IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Ref_State_Normalization") == 0) {
      i_command = 0;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.RefU;
      IP.Input_File.setf(ios::skipws);
      IP.Input_File.getline(buffer, sizeof(buffer));
      
    } else if (strcmp(IP.Next_Control_Parameter, "Space_Accuracy") == 0) {
      i_command = 210;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Space_Accuracy;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Space_Accuracy <= 0 && IP.Space_Accuracy >= 7){
	IP.Space_Accuracy = 1;
	cout << "\n Space Accuracy should be between 1 and 6 \n"
	     << "Space Accuracy set to 1" << endl;
      }/* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Accuracy_Assessment_Exact_Digits") == 0) {
      i_command = 0;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Accuracy_Assessment_Exact_Digits;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Accuracy_Assessment_Exact_Digits < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Accuracy_Assessment_Parameter") == 0) {
      i_command = 0;
      ++IP.Line_Number;
      IP.Input_File >> IP.Accuracy_Assessment_Parameter;
      IP.Input_File.getline(buffer, sizeof(buffer));
      // check that index makes sense for 2D Euler state
      if (IP.Accuracy_Assessment_Parameter < 1 || IP.Accuracy_Assessment_Parameter > 4)
	i_command = INVALID_INPUT_VALUE;

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
      } else if (strcmp(IP.Interface_IP.BC_Type,"Burning_Surface") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].BC_Type = INTERFACE_BC_BURNING_SURFACE;
      } else if (strcmp(IP.Interface_IP.BC_Type,"RinglebFlow") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].BC_Type = INTERFACE_BC_RINGLEB;
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
      } else if (strcmp(IP.Interface_IP.Motion,"BurningSurface") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_BURNING_SURFACE;
      } else if (strcmp(IP.Interface_IP.Motion,"MomentumTransfer") == 0) {
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
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Speed = Vector2D_ZERO;

    } else if (strcmp(IP.Next_Control_Parameter,"Interface_Characteristic_Velocity") == 0) {
      i_command = 607;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Interface_IP.Component_List[IP.Interface_IP.ci].Speed;
      //IP.Input_File >> IP.Interface_IP.Component_List[IP.Interface_IP.ci].Speed.x
      //	    >> IP.Interface_IP.Component_List[IP.Interface_IP.ci].Speed.y;
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

    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Nodes") == 0) {
      i_command = WRITE_OUTPUT_NODES_CODE;

    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Gradients") == 0) {
      i_command = WRITE_OUTPUT_GRADIENTS_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Quasi3D") == 0) {
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

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Restart") == 0) {
       i_command = WRITE_RESTART_CODE;

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

    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Aerodynamic_Coefficients") == 0) {
       i_command = WRITE_OUTPUT_AERODYNAMIC_COEFFICIENTS_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Refine_Grid") == 0) {
       i_command = REFINE_GRID_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Morton_Ordering") == 0) {
       i_command = MORTON_ORDERING_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Bounding_Box_Refine_Grid") == 0) {
       IP.Number_of_Bounding_Box_Mesh_Refinements = 1;
       i_command = BOUNDING_BOX_REFINE_GRID_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Determine_Maximum_Surface_Pressure") == 0) {
       i_command = DETERMINE_MAXIMUM_SURFACE_PRESSURE_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Determine_Mach_Stem_Height") == 0) {
       i_command = DETERMINE_MACH_STEM_HEIGHT_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Wedge_Solution_Distribution") == 0) {
       i_command = WRITE_OUTPUT_WEDGE_SOLUTION_DISTRIBUTION_CODE;

    } else if (strcmp(IP.Next_Control_Parameter,"Print_Accuracy") == 0) {
      i_command = WRITE_ERROR_NORMS_TO_SCREEN;
      
    } else if (strcmp(IP.Next_Control_Parameter,"Write_Accuracy_To_File") == 0) {
      i_command = WRITE_ERROR_NORMS_TO_FILE;
      
    } else if (strcmp(IP.Next_Control_Parameter,"Append_Accuracy_To_File") == 0) {
      i_command = APPEND_ERROR_NORMS_TO_FILE;

    } else if (IP.Next_Control_Parameter[0] == '#') {
       i_command = COMMENT_CODE;

    } else {
       i_command = INVALID_INPUT_CODE;

    } /* endif */


    /* Parse next control parameter with ExactSoln parser */
    IP.ExactSoln->Parse_Next_Input_Control_Parameter(IP,i_command);

    /* Parse next control parameter with CENO_Execution_Mode parser */
    CENO_Execution_Mode::Parse_Next_Input_Control_Parameter(IP,i_command);
  
    /* Parse next control parameter with CENO_Tolerances parser */
    CENO_Tolerances::Parse_Next_Input_Control_Parameter(IP,i_command);

    /* Parse next control parameter with HO_Grid2D_Execution_Mode parser */
    HO_Grid2D_Execution_Mode::Parse_Next_Input_Control_Parameter(IP,i_command);

    /* Parse next control parameter with Tecplot_Execution_Mode parser */
    Tecplot_Execution_Mode::Parse_Next_Input_Control_Parameter(IP,i_command);

    /* Parse next control parameter with HighOrder2D_Input parser */
    HighOrder2D_Input::Parse_Next_Input_Control_Parameter(IP,i_command);

    if (i_command == INVALID_INPUT_CODE) {
       // that is, we have an input line which:
       //  - is not a comment (that's COMMENT_CODE), and,
       //  - is not a valid code with an invalid value (that's INVALID_INPUT_VALUE), 
       // and so is an unknown option. Maybe it's an NKS option:
       strcpy(buffer, IP.Next_Control_Parameter);
       Get_Next_Input_Control_Parameter(IP);
       i_command = IP.NKS_IP.Parse_Next_Input_Control_Parameter(buffer, 
                                                                IP.Next_Control_Parameter);

      // If it's still unknown then ignore it. 
      // This could be a bad idea if it was an unknown command 
      // as opposed to an unknown code.
//       if (i_command == INVALID_INPUT_CODE) {
// 	cout << "\n***\n\nWarning: input file line " << IP.Line_Number << ": ";
// 	cout << "ignoring unknown input code:\n";
// 	cout << "code: " << buffer;
// 	cout << "\nvalue: " << IP.Next_Control_Parameter;
// 	cout << "\n\n***\n";
//       }
//       i_command = COMMENT_CODE; // sure why not
    }

    if (!IP.Input_File.good()) { i_command = INVALID_INPUT_VALUE; }


    /* Return the parser command type indicator. */

    return (i_command);

}

/******************************************************//**
 * Routine: Process_Input_Control_Parameter_File        
 *                                                      
 * Reads, parses, and executes the list of input        
 * control parameters from the standard input file.     
 *                                                      
 ********************************************************/
int Process_Input_Control_Parameter_File(Euler2D_Input_Parameters &Input_Parameters,
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
       cout << "\n Euler2D ERROR: Unable to open Euler2D input data file.\n";
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
          cout << "\n Euler2D ERROR: Error reading Euler2D data at line #"
               << -line_number  << " of input data file.\n";
          error_flag = line_number;
          break;
       } /* endif */
    } /* endwhile */

    /* Perform consistency checks on Input_Parameters */

    if (!error_flag &&
	Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {
      error_flag = Check_Input_Parameters<Euler2D_Input_Parameters>(Input_Parameters);
      if (error_flag) {
	cout << "\n Euler2D ERROR: Input Parameters consistency check failure\n";
	return (error_flag);
      }
    }

    // Perform consitency checks on the refinement criteria.
    Input_Parameters.Number_of_Refinement_Criteria = 0;
    if (Input_Parameters.Refinement_Criteria_Gradient_Density) Input_Parameters.Number_of_Refinement_Criteria++;
    if (Input_Parameters.Refinement_Criteria_Divergence_Velocity) Input_Parameters.Number_of_Refinement_Criteria++;
    if (Input_Parameters.Refinement_Criteria_Curl_Velocity) Input_Parameters.Number_of_Refinement_Criteria++;

    // Perform consitency checks on the time marching parameters.
    if (Input_Parameters.Time_Accurate == 1 && Input_Parameters.Local_Time_Stepping != GLOBAL_TIME_STEPPING){
      Input_Parameters.Local_Time_Stepping = GLOBAL_TIME_STEPPING;
    }

    // Set Wo state based on the final input parameters
    Input_Parameters.Wo = Euler2D_pState(Input_Parameters.Pressure/(Input_Parameters.Wo.R*Input_Parameters.Temperature), 
					 ZERO, 
					 ZERO, 
					 Input_Parameters.Pressure);
    Input_Parameters.Wo.v.x = Input_Parameters.Mach_Number*Input_Parameters.Wo.a()*cos(TWO*PI*
										       Input_Parameters.Flow_Angle/360.00);
    Input_Parameters.Wo.v.y = Input_Parameters.Mach_Number*Input_Parameters.Wo.a()*sin(TWO*PI*
										       Input_Parameters.Flow_Angle/360.00);
    
    Input_Parameters.Uo.setgas(Input_Parameters.Gas_Type);
    Input_Parameters.Uo = U(Input_Parameters.Wo);

    Input_Parameters.Wo.Mr_min = Input_Parameters.Mr_Min_Factor*Input_Parameters.Mach_Number;
    Input_Parameters.Uo.Mr_min = Input_Parameters.Mr_Min_Factor*Input_Parameters.Mach_Number;

    // Perform update of the internal variables of the exact solution
    Input_Parameters.ExactSoln->Set_ParticularSolution_Parameters(Input_Parameters);

    // Perform update of the internal variables of the high-order input parameters
    HighOrder2D_Input::Set_Final_Parameters(Input_Parameters);

    // Set reference state in the Euler2D_Quad_Block class
    Euler2D_Quad_Block::Set_Normalization_Reference_State(Input_Parameters.RefU);

    // Set limiter in CENO class
    CENO_Execution_Mode::Limiter = Input_Parameters.i_Limiter;

    /* Initial processing of input control parameters complete.  
       Return the error indicator flag. */

    return (error_flag);

}
