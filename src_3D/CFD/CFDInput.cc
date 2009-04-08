/* CFDInput.cc: Definition of CFD_Input_Parameters class member functions. */

/* Include the CFDInput header file. */

#ifndef _CFD_INPUT_INCLUDED
#include "CFDInput.h"
#endif // _CFD_INPUT_INCLUDED

/* Define member functions. */

/********************************************************
 * Routine: Get_CFFC_Path                               *
 *                                                      *
 * Obtains the location of CFFC source directory        *
 *                                                      *
 ********************************************************/
void CFD_Input_Parameters::Get_CFFC_Path(void){

  char *string_ptr;
  // Check to see if environment varible exists.
  if (getenv(PATHVAR) == NULL) {
    //Set default path
     string_ptr = "CFFC";
     strcpy(CFFC_Path, string_ptr);
  } else {
     //Set path specified by environment variable
     strcpy(CFFC_Path, getenv(PATHVAR));
  } /* endif */

}

/********************************************************
 * Routine: Open_Input_File                             *
 *                                                      *
 * Opens the appropriate input data file.               *
 *                                                      *
 ********************************************************/
void CFD_Input_Parameters::Open_Input_File(void) {

    Input_File.open(Input_File_Name, ios::in);
    if (!Input_File.fail()) {
       Line_Number = 0;
       Input_File.setf(ios::skipws);
    } /* endif */

}

/********************************************************
 * Routine: Close_Input_File                            *
 *                                                      *
 * Closes the appropriate input data file.              *
 *                                                      *
 ********************************************************/
void CFD_Input_Parameters::Close_Input_File(void) {

    Input_File.unsetf(ios::skipws);
    Input_File.close();

}

/********************************************************
 * Routine: Get_Next_Input_Control_Parameter            *
 *                                                      *
 * Get the next input control parameter from the input  *
 * file.  Returns non-zero value if a read error is     *
 * encountered.                                         *
 *                                                      *
 ********************************************************/
int CFD_Input_Parameters::Get_Next_Input_Control_Parameter(const bool read_control_parameter) {

    int i;
    char buffer[256], buffer_no_leading_blanks[256];

    Line_Number = Line_Number + 1;
    Input_File.getline(buffer, sizeof(buffer));

    // Remove leading blanks and tabs.
    i = 0;
    if (buffer[0] == ' ' || buffer[0] == '\t') {
       while (1) {
	  if (buffer[i] != ' ' && buffer[i] != '\t') break;
          i = i + 1;
          if (i > strlen(buffer) ) break;
       } /* endwhile */
       for (int k = i; k <= strlen(buffer); k++) {
	  buffer_no_leading_blanks[k-i] = buffer[k];
       } /* endfor */
       buffer_no_leading_blanks[strlen(buffer)-1] = '\0';
       strcpy(buffer, buffer_no_leading_blanks);
    } /* endif */

    // Obtain next control parameter.
    i = 0;
    if (buffer[0] != '#') {
       while (1) {
	  if (read_control_parameter && 
              (buffer[i] == ' ' || buffer[i] == '=' )) break;
          i = i + 1;
          if (i > strlen(buffer) ) break;
       } /* endwhile */
       buffer[i] = '\0';
    } /* endif */
    strcpy(Next_Control_Parameter, buffer);
    
    // Check for long input lines that have not been completely read in.
    if (!Input_File.good()) { 
       return (1);
    } else {
       return (0);
    } /* endif */

}

/********************************************************
 * Routine: Parse_Next_Input_Control_Parameter          *
 *                                                      *
 * Parses and executes the next input control parameter *
 * from the input file.                                 *
 *                                                      *
 ********************************************************/
int CFD_Input_Parameters::Parse_Next_Input_Control_Parameter(void) {

    int i_command, error_flag;
    char code[256];
    string value_string;
    stringstream value_stream;

    i_command = 0;
    strcpy(code, Next_Control_Parameter);
   
    //
    //  Parse execution control parameters.
    //       
    if (strcmp(code, "Execute") == 0) {
       i_command = EXECUTE_CODE;
            
    } else if (strcmp(code, "Terminate") == 0) {
       i_command = TERMINATE_CODE;
       
    } else if (strcmp(code, "Continue") == 0) {
       i_command = CONTINUE_CODE;
       
    } else if (strcmp(code, "Write_Output") == 0) {
       i_command = WRITE_OUTPUT_CODE;
       
    } else if (strcmp(code, "Write_Output_Cells") == 0) {
       i_command = WRITE_OUTPUT_CELLS_CODE;
       
    } else if (strcmp(code, "Write_Output_Nodes") == 0) {
       i_command = WRITE_OUTPUT_NODES_CODE;

    } else if (strcmp(code, "Write_Restart") == 0) {
       i_command = WRITE_RESTART_CODE;
       
    } else if (strcmp(code, "Write_Output_RHS") == 0) {
       i_command = WRITE_OUTPUT_RHS_CODE;
       
    } else if (strcmp(code, "Write_Output_Mesh") == 0) {
       i_command = WRITE_OUTPUT_GRID_CODE;

    } else if (strcmp(code, "Perturbation") == 0) {
       i_command = WRITE_OUTPUT_PERTURB_CODE;

    } else if (strcmp(code, "Write_Mesh_Definition") == 0) {
       i_command = WRITE_GRID_DEFINITION_CODE;
       
    } else if (strcmp(code, "Write_Output_Mesh_Nodes") == 0) {
       i_command = WRITE_OUTPUT_GRID_NODES_CODE;

    } else if (strcmp(code, "Write_Output_Mesh_Cells") == 0) {
       i_command = WRITE_OUTPUT_GRID_CELLS_CODE;

    }  else if (strcmp(code, "Morton_Ordering") == 0) {
       i_command = MORTON_ORDERING_CODE;

    } else if (code[0] == '#') {
       i_command = COMMENT_CODE;

    } else if (strcmp(code, "Refine_Grid") == 0) {
       i_command = REFINE_GRID_CODE;

    } else {
       i_command = INVALID_INPUT_CODE;

    } /* endif */

    // If valid execution control command, then return.
    if (i_command != INVALID_INPUT_CODE) {
       return (i_command);
    } /* endif */

    //
    // Else continue parsing.
    //
    error_flag = Get_Next_Input_Control_Parameter(false);
    while (!error_flag &&
           Next_Control_Parameter[0] == '#') {
       error_flag = Get_Next_Input_Control_Parameter(false);
    } /* endwhile */
    if (error_flag) {
       i_command = INVALID_INPUT_VALUE;
       return (i_command);
    } else {
       value_stream << Next_Control_Parameter;
    } /* endif */

    //
    // Input file parameters:
    //
    if (strcmp(code, "CFFC_Path") == 0) {
       i_command = 1;
       value_stream >> value_string;
       strcpy(CFFC_Path, value_string.c_str());

    //
    // Output parameters:
    //
    } else if (strcmp(code, "Output_File_Name") == 0) {
       i_command = 10;
       value_stream >> value_string;
       strcpy(Output_File_Name, value_string.c_str());
       strcat(Output_File_Name, ".dat");
       strcpy(Restart_File_Name, value_string.c_str());
       strcat(Restart_File_Name, ".soln");
      
    } else if (strcmp(code, "Restart_File_Name") == 0) {
       i_command = 11;
       value_stream >> value_string;
       strcpy(Restart_File_Name, value_string.c_str());
       strcat(Restart_File_Name, ".soln");

    } else if (strcmp(code, "Output_Format_Type") == 0) {
       i_command = 12;
       value_stream >> value_string;
       strcpy(Output_Format_Type, value_string.c_str());
       if (strcmp(Output_Format_Type, "Gnuplot") == 0  ||
           strcmp(Output_Format_Type, "gnuplot") == 0  ||
           strcmp(Output_Format_Type, "GNUPLOT") == 0) {
          i_Output_Format = IO_GNUPLOT;
       } else if (strcmp(Output_Format_Type, "Tecplot") == 0  ||
                  strcmp(Output_Format_Type, "tecplot") == 0  ||
                  strcmp(Output_Format_Type, "TECPLOT") == 0) {
          i_Output_Format = IO_TECPLOT;
       } else if (strcmp(Output_Format_Type, "Matlab") == 0  ||
                  strcmp(Output_Format_Type, "matlab") == 0  ||
                  strcmp(Output_Format_Type, "MATLAB") == 0) {
          i_Output_Format = IO_MATLAB;
       } else if (strcmp(Output_Format_Type, "Octave") == 0  ||
                  strcmp(Output_Format_Type, "octave") == 0  ||
                  strcmp(Output_Format_Type, "OCTAVE") == 0) {
          i_Output_Format = IO_OCTAVE;
       } else {
          i_command = INVALID_INPUT_VALUE;
       } /* endif */

    } else if (strcmp(code, "Restart_Solution_Save_Frequency") == 0) {
       i_command = 13;
       value_stream >> Restart_Solution_Save_Frequency;
       if (Restart_Solution_Save_Frequency < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(code, "Time_Accurate_Output_Frequency") == 0) {
       i_command = 14; 
       value_stream >> Time_Accurate_Output_Frequency;
       if (Time_Accurate_Output_Frequency < 0) i_command = INVALID_INPUT_VALUE;

    //
    // Debugging parameters:
    //
    } else if (strcmp(code, "Debug_Level") == 0) {
       i_command = 20;
       value_stream >> Debug_Level;
       if (Debug_Level != 0 && Debug_Level != 1 ) i_command = INVALID_INPUT_VALUE; 

    //
    // Flow type indicator and related input parameters:
    //
    } else if (strcmp(code, "Flow_Type") == 0) {
       i_command = 30;
       value_stream >> value_string;
       strcpy(Flow_Type, value_string.c_str());
       if (strcmp(Flow_Type, "Inviscid") == 0) {
          i_Flow_Type = FLOWTYPE_INVISCID;
       } else if (strcmp(Flow_Type, "Laminar") == 0) {
          i_Flow_Type = FLOWTYPE_LAMINAR;
       } else if (strcmp(Flow_Type, "Turbulent-k-epsilon") == 0) {
          i_Flow_Type = FLOWTYPE_TURBULENT_RANS_K_EPSILON;
       } else if (strcmp(Flow_Type, "Turbulent-k-omega") == 0) {
          i_Flow_Type = FLOWTYPE_TURBULENT_RANS_K_OMEGA;
       } else if (strcmp(Flow_Type, "Turbulent-LES") == 0) {
          i_Flow_Type = FLOWTYPE_TURBULENT_LES;
       } else if (strcmp(Flow_Type, "Turbulent-LES-C-Fsd-k") == 0) {
	 i_Flow_Type = FLOWTYPE_TURBULENT_LES_C_FSD_K;
       } else if (strcmp(Flow_Type, "Turbulent-LES-C-Fsd-Smagorinsky") == 0) {
	 i_Flow_Type = FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY;
       } else if (strcmp(Flow_Type, "Turbulent-LES-TF-k") == 0) {
	 i_Flow_Type = FLOWTYPE_TURBULENT_LES_TF_K;
       } else if (strcmp(Flow_Type, "Turbulent-LES-TF-Smagorinsky") == 0) {
	 i_Flow_Type = FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY;
       } else if (strcmp(Flow_Type, "Turbulent-DES-k-omega") == 0) {
          i_Flow_Type = FLOWTYPE_TURBULENT_DES_K_OMEGA;
       } else if (strcmp(Flow_Type, "Turbulent-DNS") == 0) {
          i_Flow_Type = FLOWTYPE_TURBULENT_DNS;
       } else {
          i_command = INVALID_INPUT_VALUE;
       } /* endif */

    } else if (strcmp(code, "Axisymmetric") == 0) {
       i_command = 31;
       value_stream >> value_string;
       if (value_string == "ON") {
         Axisymmetric = ON;
       } else if (value_string == "OFF") {
         Axisymmetric = OFF;
       } else {
         i_command = INVALID_INPUT_VALUE;
       } /* endif */

    //
    // Time integration type indicator and related input parameters:
    //
    } else if (strcmp(code, "Time_Integration_Type") == 0) {
       i_command = 40;
       value_stream >> value_string;
       strcpy(Time_Integration_Type, value_string.c_str());
       if (strcmp(Time_Integration_Type, "Explicit_Euler") == 0) {
          i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
          N_Stage = 1;
       } else if (strcmp(Time_Integration_Type, "Explicit_Predictor_Corrector") == 0) {
          i_Time_Integration = TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR;
          N_Stage = 2;
       } else if (strcmp(Time_Integration_Type, "Explicit_Runge_Kutta") == 0) {
          i_Time_Integration = TIME_STEPPING_EXPLICIT_RUNGE_KUTTA;
          N_Stage = 4;
       } else if (strcmp(Time_Integration_Type, "Multistage_Optimal_Smoothing") == 0) {
          i_Time_Integration = TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING;
          N_Stage = 4;
       } else {
          i_command = INVALID_INPUT_VALUE;
       } /* endif */

    } else if (strcmp(code, "Time_Accurate") == 0) {
       i_command = 41;
       value_stream >> value_string;
       if (value_string == "ON") {
         Time_Accurate = ON;
       } else if (value_string == "OFF") {
         Time_Accurate = OFF;
       } else {
         i_command = INVALID_INPUT_VALUE;
       } /* endif */
       if (Time_Accurate) {
          Local_Time_Stepping = GLOBAL_TIME_STEPPING;
       } /* endif */

    } else if (strcmp(code, "Local_Time_Stepping") == 0) {
       i_command = 42;
       value_stream >> Local_Time_Stepping;
       Preconditioning = 0;
       if (Local_Time_Stepping != GLOBAL_TIME_STEPPING &&
           Local_Time_Stepping != SCALAR_LOCAL_TIME_STEPPING &&
	   Local_Time_Stepping != MATRIX_LOCAL_TIME_STEPPING &&	   
           Local_Time_Stepping != LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER &&
	   Local_Time_Stepping != SEMI_IMPLICIT_LOCAL_TIME_STEPPING &&
	   Local_Time_Stepping != SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER) {
         Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
       }
       if(Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER ||
       	  Local_Time_Stepping == SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER){
	 Preconditioning = 1;
       }       

    } else if (strcmp(code, "Maximum_Number_of_Time_Steps") == 0) {
       i_command = 43;
       value_stream >> Maximum_Number_of_Time_Steps;
       if (Maximum_Number_of_Time_Steps < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(code, "N_Stage") == 0) {
       i_command = 44;
       value_stream >> N_Stage;
       if (N_Stage < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(code, "CFL_Number") == 0) {
       i_command = 45;
       value_stream >> CFL_Number;
       if (CFL_Number <= ZERO) i_command = INVALID_INPUT_VALUE;

    }  else if (strcmp(code, "Time_Max") == 0) {
       i_command = 46;
       value_stream >> Time_Max;
       Time_Max = Time_Max/THOUSAND;
       if (Time_Max < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(code, "Residual_Norm") == 0) {
       i_command = 47;
       value_stream >> Residual_Norm;
       if (Residual_Norm < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(code, "Number_of_Residual_Norms") == 0) {
       i_command = 48;
       value_stream >> Number_of_Residual_Norms;
       if (Number_of_Residual_Norms < 0) i_command = INVALID_INPUT_VALUE;
    //
    // Spatial Order of Accuracy type indicator and related input parameters:
    //
    } else if (strcmp(code, "Spatial_Order_of_Accuracy") == 0) {
      i_command = 49;
      value_stream >> value_string;
      strcpy(Spatial_Accuracy, value_string.c_str());
      if (strcmp(Spatial_Accuracy, "1") == 0) {
	Reconstruction_Order = 0;
      } else if (strcmp(Spatial_Accuracy, "2") == 0) {
	Reconstruction_Order = 1;
      } else if (strcmp(Spatial_Accuracy, "3") == 0) {
	Reconstruction_Order = 2;
      } else if (strcmp(Spatial_Accuracy, "4") == 0) {
	Reconstruction_Order = 3;
      } else {
	i_command = INVALID_INPUT_VALUE;
      } /* endif */
      Grid3D_HO_Execution_Mode::RECONSTRUCTION_ORDER = Reconstruction_Order;
    //
    // Reconstruction type indicator and related input parameters:
    //
    } else if (strcmp(code, "Reconstruction_Type") == 0) {
       i_command = 50;
       value_stream >> value_string;
       strcpy(Reconstruction_Type, value_string.c_str());
       if (strcmp(Reconstruction_Type, "Green_Gauss") == 0) {
          i_Reconstruction = RECONSTRUCTION_GREEN_GAUSS;
       } else if (strcmp(Reconstruction_Type, "Least_Squares") == 0) {
          i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;
       } else if (strcmp(Reconstruction_Type, "CENO") == 0) {
	  i_Reconstruction = RECONSTRUCTION_HIGH_ORDER;
	  Grid3D_HO_Execution_Mode::USE_HO_CENO_GRID = ON;
       } else {
          i_command = INVALID_INPUT_VALUE;
       } /* endif */

    //
    // Limiter type indicator and related input parameters:
    //
    } else if (strcmp(code, "Limiter_Type") == 0) {
       i_command = 60;
       value_stream >> value_string;
       strcpy(Limiter_Type, value_string.c_str());
       if (strcmp(Limiter_Type, "One") == 0) {
          i_Limiter = LIMITER_ONE;
       } else if (strcmp(Limiter_Type, "Zero") == 0) {
          i_Limiter = LIMITER_ZERO;
       } else if (strcmp(Limiter_Type, "VanLeer") == 0) {
          i_Limiter = LIMITER_VANLEER;
       } else if (strcmp(Limiter_Type, "VanAlbada") == 0) {
          i_Limiter = LIMITER_VANALBADA;
       } else if (strcmp(Limiter_Type, "Barth_Jespersen") == 0) {
          i_Limiter = LIMITER_BARTH_JESPERSEN;
       } else if (strcmp(Limiter_Type, "Venkatakrishnan") == 0) {
          i_Limiter = LIMITER_VENKATAKRISHNAN;
       } else {
          i_command = INVALID_INPUT_VALUE;
       } /* endif */

    } else if (strcmp(code, "Freeze_Limiter") == 0) {
      i_command = 61;
      value_stream >> Freeze_Limiter;
      if (Freeze_Limiter < ZERO) i_command = INVALID_INPUT_VALUE;
      
    } else if (strcmp(code, "Freeze_Limiter_Residual_Level") == 0) {
      i_command = 62;
      value_stream >> Freeze_Limiter_Residual_Level;
      if (Freeze_Limiter_Residual_Level < 0) i_command = INVALID_INPUT_VALUE;

    //
    // Inviscid flux function type and related input parameters:
    //
    } else if (strcmp(code, "Flux_Function_Type") == 0) {
       i_command = 70;
       value_stream >> value_string;
       strcpy(Flux_Function_Type, value_string.c_str());
       if (strcmp(Flux_Function_Type, "Godunov") == 0) {
          i_Flux_Function = FLUX_FUNCTION_GODUNOV;
       } else if (strcmp(Flux_Function_Type, "Roe") == 0) {
          i_Flux_Function = FLUX_FUNCTION_ROE;
       } else if (strcmp(Flux_Function_Type, "Rusanov") == 0) {
          i_Flux_Function = FLUX_FUNCTION_RUSANOV;
       } else if (strcmp(Flux_Function_Type, "HLLE") == 0) {
          i_Flux_Function = FLUX_FUNCTION_HLLE;
       } else if (strcmp(Flux_Function_Type, "Linde") == 0) {
          i_Flux_Function = FLUX_FUNCTION_LINDE;
       } else if (strcmp(Flux_Function_Type, "HLLC") == 0) {
          i_Flux_Function = FLUX_FUNCTION_HLLC;
       } else if (strcmp(Flux_Function_Type, "AUSM_plus_up") == 0) {
          i_Flux_Function = FLUX_FUNCTION_AUSM_PLUS_UP;
       } else {
          i_command = INVALID_INPUT_VALUE;
       } /* endif */

    //
    // Preconditioner indicator and related input parameters:
    //
    } else if (strcmp(code, "Mach_Number_Reference") == 0) {
       i_command = 80;
       value_stream >> Mach_Number_Reference; 
       if (Mach_Number_Reference < ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } /* endif */

    //
    // Implicit residual smoothing control parameters:
    //
    } else if (strcmp(code, "Residual_Smoothing") == 0) {
       i_command = 90;
       value_stream >> value_string;
       if (value_string == "ON") {
         Residual_Smoothing = ON;
       } else if (value_string == "OFF") {
         Residual_Smoothing = OFF;
       } else {
         i_command = INVALID_INPUT_VALUE;
       } /* endif */

    } else if (strcmp(code, "Residual_Smoothing_Gauss_Seidel_Iterations") == 0) {
       i_command = 91;
       value_stream >> Residual_Smoothing_Gauss_Seidel_Iterations;
       if (Residual_Smoothing_Gauss_Seidel_Iterations < 0) {
          Residual_Smoothing_Gauss_Seidel_Iterations = 0;
       } /* endif */
   
    } else if (strcmp(code, "Residual_Smoothing_Epsilon") == 0) {
       i_command = 92;
       value_stream >> Residual_Smoothing_Epsilon;
       if (Residual_Smoothing_Epsilon <= ZERO) {
          Residual_Smoothing = OFF;
          Residual_Smoothing_Epsilon = ZERO;
       } else {
          Residual_Smoothing = ON;
       } /* endif */

    //
    // Initial condition type indicator and related input parameters:
    //
    } else if (strcmp(code, "ICs_Type") == 0) {
       i_command = 100;
       value_stream >> value_string;
       strcpy(ICs_Type, value_string.c_str());
       if (strcmp(ICs_Type, "Constant") == 0) {
          i_ICs = IC_CONSTANT;
       } else if (strcmp(ICs_Type, "Uniform") == 0) {
          i_ICs = IC_UNIFORM;
       } else if (strcmp(ICs_Type, "Sod") == 0) {
          i_ICs = IC_SOD;
       } else if (strcmp(ICs_Type, "Sod_Xdir") == 0) {
          i_ICs = IC_SOD_XDIR;
       } else if (strcmp(ICs_Type, "Sod_Ydir") == 0) {
          i_ICs = IC_SOD_YDIR;
       } else if (strcmp(ICs_Type, "Sod_Zdir") == 0) {
          i_ICs = IC_SOD_ZDIR;
       } else if (strcmp(ICs_Type, "Groth") == 0) {
          i_ICs = IC_GROTH;
       } else if (strcmp(ICs_Type, "Groth_Xdir") == 0) {
          i_ICs = IC_GROTH_XDIR;
       } else if (strcmp(ICs_Type, "Groth_Ydir") == 0) {
          i_ICs = IC_GROTH_YDIR;
       } else if (strcmp(ICs_Type, "Einfeldt") == 0) {
          i_ICs = IC_EINFELDT;
       } else if (strcmp(ICs_Type, "Einfeldt_Xdir") == 0) {
          i_ICs = IC_EINFELDT_XDIR;
       } else if (strcmp(ICs_Type, "Einfeldt_Ydir") == 0) {
          i_ICs = IC_EINFELDT_YDIR;
       } else if (strcmp(ICs_Type, "Shock_Wave_Xdir") == 0) {
          i_ICs = IC_SHOCK_WAVE_XDIR;
       } else if (strcmp(ICs_Type, "Shock_Wave_Ydir") == 0) {
          i_ICs = IC_SHOCK_WAVE_YDIR;
       } else if (strcmp(ICs_Type, "Contact_Surface_Xdir") == 0) {
          i_ICs = IC_CONTACT_SURFACE_XDIR;
       } else if (strcmp(ICs_Type, "Contact_Surface_Ydir") == 0) {
          i_ICs = IC_CONTACT_SURFACE_YDIR;
       } else if (strcmp(ICs_Type, "Rarefaction_Wave_Xdir") == 0) {
          i_ICs = IC_RAREFACTION_WAVE_XDIR;
       } else if (strcmp(ICs_Type, "Rarefaction_Wave_Ydir") == 0) {
          i_ICs = IC_RAREFACTION_WAVE_YDIR;
       } else if (strcmp(ICs_Type, "ShockBox") == 0) {
          i_ICs = IC_SHOCK_BOX;
       } else if (strcmp(ICs_Type, "ShockBox_XY") == 0) {
           i_ICs = IC_SHOCK_BOX_XY;
       } else if (strcmp(ICs_Type, "ShockBox_XZ") == 0) {
           i_ICs = IC_SHOCK_BOX_XZ;
       } else if (strcmp(ICs_Type, "ShockBox_YZ") == 0) {
           i_ICs = IC_SHOCK_BOX_YZ;
       } else if (strcmp(ICs_Type, "High_Pressure_Reservoir") == 0) {
          i_ICs = IC_HIGH_PRESSURE_RESERVOIR;
       } else if (strcmp(ICs_Type, "Low_Pressure_Reservoir") == 0) {
          i_ICs = IC_LOW_PRESSURE_RESERVOIR;
       } else if (strcmp(ICs_Type, "Riemann") == 0) {
          i_ICs = IC_RIEMANN;
       } else if (strcmp(ICs_Type, "Riemann_Xdir") == 0) {
          i_ICs = IC_RIEMANN_XDIR;
       } else if (strcmp(ICs_Type, "Riemann_Ydir") == 0) {
          i_ICs = IC_RIEMANN_YDIR;  
       } else if (strcmp(ICs_Type, "Wedge_Flow") == 0) {
          i_ICs = IC_WEDGE_FLOW;	 
       } else if (strcmp(ICs_Type, "Mix") == 0) {
          i_ICs = IC_GAS_MIX;
       } else if (strcmp(ICs_Type, "Core_Flame") == 0 ) {
          i_ICs = IC_CHEM_CORE_FLAME ;
       } else if (strcmp(ICs_Type, "Inverse_Flame") == 0 ) {
          i_ICs = IC_CHEM_INVERSE_FLAME ; 
       } else if (strcmp(ICs_Type, "Pressure_Gradient_x") == 0 ) {
          i_ICs = IC_PRESSURE_GRADIENT_X;
       } else if (strcmp(ICs_Type, "Pressure_Gradient_y") == 0 ) {
          i_ICs = IC_PRESSURE_GRADIENT_Y;
       } else if (strcmp(ICs_Type, "Pressure_Gradient_z") == 0 ) {
          i_ICs = IC_PRESSURE_GRADIENT_Z;
       } else if (strcmp(ICs_Type, "Couette") == 0 ) {
          i_ICs = IC_VISCOUS_COUETTE; 
       } else if (strcmp(ICs_Type, "Couette_Pressure_Gradient_x") == 0 ) {
          i_ICs = IC_VISCOUS_COUETTE_PRESSURE_GRADIENT_X;
       } else if (strcmp(ICs_Type, "Couette_Pressure_Gradient_y") == 0 ) {
          i_ICs = IC_VISCOUS_COUETTE_PRESSURE_GRADIENT_Y;
       } else if (strcmp(ICs_Type, "Couette_Pressure_Gradient_z") == 0 ) {
          i_ICs = IC_VISCOUS_COUETTE_PRESSURE_GRADIENT_Z;
       } else if (strcmp(ICs_Type, "1DPremixedFlame") == 0 ) {
          i_ICs = IC_CHEM_1DFLAME; 
       } else if (strcmp(ICs_Type, "Flat_Plate") == 0) {
          i_ICs = IC_VISCOUS_FLAT_PLATE;
       } else if (strcmp(ICs_Type, "Pipe_Flow") == 0) {
          i_ICs = IC_TURBULENT_PIPE_FLOW;
       } else if (strcmp(ICs_Type, "Coflow") == 0) {
          i_ICs = IC_TURBULENT_COFLOW;
       } else if (strcmp(ICs_Type, "Driven_Cavity_Flow") == 0) {
          i_ICs = IC_VISCOUS_DRIVEN_CAVITY_FLOW;
       } else if (strcmp(ICs_Type, "Turbulent_Dump_Combustor") == 0) {
          i_ICs = IC_TURBULENT_DUMP_COMBUSTOR;
       } else if (strcmp(ICs_Type, "Channel_Flow") == 0) {
          i_ICs = IC_CHANNEL_FLOW;
       } else if (strcmp(ICs_Type, "Laminar_Channel_Flow") == 0) {
          i_ICs = IC_CHANNEL_FLOW;
       } else if (strcmp(ICs_Type, "Turbulent_Channel_Flow") == 0) {
          i_ICs = IC_CHANNEL_FLOW;
       } else if (strcmp(ICs_Type, "Turbulent_Diffusion_Flame") == 0) {
          i_ICs = IC_TURBULENT_DIFFUSION_FLAME;
       } else if (strcmp(ICs_Type, "Turbulent_Free_Jet_Flame") == 0) {
          i_ICs = IC_FREE_JET_FLAME;
       } else if (strcmp(ICs_Type, "Turbulent_Premixed_Flame") == 0) {
          i_ICs = IC_TURBULENT_PREMIXED_FLAME;
       }else if (strcmp(ICs_Type, "Turbulent_Bunsen_Flame") == 0) {
          i_ICs = IC_TURBULENT_BUNSEN_FLAME;
       }else if (strcmp(ICs_Type, "Turbulent_Bunsen_Box") == 0) {
          i_ICs = IC_TURBULENT_BUNSEN_BOX;
       }else if (strcmp(ICs_Type, "Turbulent_Box") == 0) {
          i_ICs = IC_TURBULENT_BOX;
       } else if (strcmp(ICs_Type, "Restart") == 0) {
          i_ICs = IC_RESTART;
       } else {
          i_command = INVALID_INPUT_VALUE;
       } /* endif */
       if (i_ICs != IC_RESTART) {
           i_Original_ICs = i_ICs;
           strcpy(Original_ICs_Type, ICs_Type);
       }
        
    } else if (strcmp(code, "Original_ICs_Type") == 0) {
        i_command = 100;
        value_stream >> value_string;
        strcpy(Original_ICs_Type, value_string.c_str());
        if (strcmp(Original_ICs_Type, "Constant") == 0) {
            i_Original_ICs = IC_CONSTANT;
        } else if (strcmp(Original_ICs_Type, "Uniform") == 0) {
            i_Original_ICs = IC_UNIFORM;
        } else if (strcmp(Original_ICs_Type, "Sod") == 0) {
            i_Original_ICs = IC_SOD;
        } else if (strcmp(Original_ICs_Type, "Sod_Xdir") == 0) {
            i_Original_ICs = IC_SOD_XDIR;
        } else if (strcmp(Original_ICs_Type, "Sod_Ydir") == 0) {
            i_Original_ICs = IC_SOD_YDIR;
        } else if (strcmp(Original_ICs_Type, "Sod_Zdir") == 0) {
            i_Original_ICs = IC_SOD_ZDIR;
        } else if (strcmp(Original_ICs_Type, "Groth") == 0) {
            i_Original_ICs = IC_GROTH;
        } else if (strcmp(Original_ICs_Type, "Groth_Xdir") == 0) {
            i_Original_ICs = IC_GROTH_XDIR;
        } else if (strcmp(Original_ICs_Type, "Groth_Ydir") == 0) {
            i_Original_ICs = IC_GROTH_YDIR;
        } else if (strcmp(Original_ICs_Type, "Einfeldt") == 0) {
            i_Original_ICs = IC_EINFELDT;
        } else if (strcmp(Original_ICs_Type, "Einfeldt_Xdir") == 0) {
            i_Original_ICs = IC_EINFELDT_XDIR;
        } else if (strcmp(Original_ICs_Type, "Einfeldt_Ydir") == 0) {
            i_Original_ICs = IC_EINFELDT_YDIR;
        } else if (strcmp(Original_ICs_Type, "Shock_Wave_Xdir") == 0) {
            i_Original_ICs = IC_SHOCK_WAVE_XDIR;
        } else if (strcmp(Original_ICs_Type, "Shock_Wave_Ydir") == 0) {
            i_Original_ICs = IC_SHOCK_WAVE_YDIR;
        } else if (strcmp(Original_ICs_Type, "Contact_Surface_Xdir") == 0) {
            i_Original_ICs = IC_CONTACT_SURFACE_XDIR;
        } else if (strcmp(Original_ICs_Type, "Contact_Surface_Ydir") == 0) {
            i_Original_ICs = IC_CONTACT_SURFACE_YDIR;
        } else if (strcmp(Original_ICs_Type, "Rarefaction_Wave_Xdir") == 0) {
            i_Original_ICs = IC_RAREFACTION_WAVE_XDIR;
        } else if (strcmp(Original_ICs_Type, "Rarefaction_Wave_Ydir") == 0) {
            i_Original_ICs = IC_RAREFACTION_WAVE_YDIR;
        } else if (strcmp(Original_ICs_Type, "ShockBox") == 0) {
            i_Original_ICs = IC_SHOCK_BOX;
        } else if (strcmp(Original_ICs_Type, "ShockBox_XY") == 0) {
            i_Original_ICs = IC_SHOCK_BOX_XY;
        } else if (strcmp(Original_ICs_Type, "ShockBox_XZ") == 0) {
            i_Original_ICs = IC_SHOCK_BOX_XZ;
        } else if (strcmp(Original_ICs_Type, "ShockBox_YZ") == 0) {
            i_Original_ICs = IC_SHOCK_BOX_YZ;
        } else if (strcmp(Original_ICs_Type, "High_Pressure_Reservoir") == 0) {
            i_Original_ICs = IC_HIGH_PRESSURE_RESERVOIR;
        } else if (strcmp(Original_ICs_Type, "Low_Pressure_Reservoir") == 0) {
            i_Original_ICs = IC_LOW_PRESSURE_RESERVOIR;
        } else if (strcmp(Original_ICs_Type, "Riemann") == 0) {
            i_Original_ICs = IC_RIEMANN;
        } else if (strcmp(Original_ICs_Type, "Riemann_Xdir") == 0) {
            i_Original_ICs = IC_RIEMANN_XDIR;
        } else if (strcmp(Original_ICs_Type, "Riemann_Ydir") == 0) {
            i_Original_ICs = IC_RIEMANN_YDIR;  
        } else if (strcmp(Original_ICs_Type, "Wedge_Flow") == 0) {
            i_Original_ICs = IC_WEDGE_FLOW;	 
        } else if (strcmp(Original_ICs_Type, "Mix") == 0) {
            i_Original_ICs = IC_GAS_MIX;
        } else if (strcmp(Original_ICs_Type, "Core_Flame") == 0 ) {
            i_Original_ICs = IC_CHEM_CORE_FLAME ;
        } else if (strcmp(Original_ICs_Type, "Inverse_Flame") == 0 ) {
            i_Original_ICs = IC_CHEM_INVERSE_FLAME ; 
        } else if (strcmp(Original_ICs_Type, "Pressure_Gradient_x") == 0 ) {
            i_Original_ICs = IC_PRESSURE_GRADIENT_X;
        } else if (strcmp(Original_ICs_Type, "Pressure_Gradient_y") == 0 ) {
            i_Original_ICs = IC_PRESSURE_GRADIENT_Y;
        } else if (strcmp(Original_ICs_Type, "Pressure_Gradient_z") == 0 ) {
            i_Original_ICs = IC_PRESSURE_GRADIENT_Z;
        } else if (strcmp(Original_ICs_Type, "Couette") == 0 ) {
            i_Original_ICs = IC_VISCOUS_COUETTE; 
        } else if (strcmp(Original_ICs_Type, "Couette_Pressure_Gradient_x") == 0 ) {
            i_Original_ICs = IC_VISCOUS_COUETTE_PRESSURE_GRADIENT_X;
        } else if (strcmp(Original_ICs_Type, "Couette_Pressure_Gradient_y") == 0 ) {
            i_Original_ICs = IC_VISCOUS_COUETTE_PRESSURE_GRADIENT_Y;
        } else if (strcmp(Original_ICs_Type, "Couette_Pressure_Gradient_z") == 0 ) {
            i_Original_ICs = IC_VISCOUS_COUETTE_PRESSURE_GRADIENT_Z;
        } else if (strcmp(Original_ICs_Type, "1DPremixedFlame") == 0 ) {
            i_Original_ICs = IC_CHEM_1DFLAME; 
        } else if (strcmp(Original_ICs_Type, "Flat_Plate") == 0) {
            i_Original_ICs = IC_VISCOUS_FLAT_PLATE;
        } else if (strcmp(Original_ICs_Type, "Pipe_Flow") == 0) {
            i_Original_ICs = IC_TURBULENT_PIPE_FLOW;
        } else if (strcmp(Original_ICs_Type, "Coflow") == 0) {
            i_Original_ICs = IC_TURBULENT_COFLOW;
        } else if (strcmp(Original_ICs_Type, "Driven_Cavity_Flow") == 0) {
            i_Original_ICs = IC_VISCOUS_DRIVEN_CAVITY_FLOW;
        } else if (strcmp(Original_ICs_Type, "Turbulent_Dump_Combustor") == 0) {
            i_Original_ICs = IC_TURBULENT_DUMP_COMBUSTOR;
        } else if (strcmp(Original_ICs_Type, "Channel_Flow") == 0) {
            i_Original_ICs = IC_CHANNEL_FLOW;
        } else if (strcmp(Original_ICs_Type, "Laminar_Channel_Flow") == 0) {
            i_Original_ICs = IC_CHANNEL_FLOW;
        } else if (strcmp(Original_ICs_Type, "Turbulent_Channel_Flow") == 0) {
            i_Original_ICs = IC_CHANNEL_FLOW;
        } else if (strcmp(Original_ICs_Type, "Turbulent_Diffusion_Flame") == 0) {
            i_Original_ICs = IC_TURBULENT_DIFFUSION_FLAME;
        } else if (strcmp(Original_ICs_Type, "Turbulent_Free_Jet_Flame") == 0) {
            i_Original_ICs = IC_FREE_JET_FLAME;
        } else if (strcmp(Original_ICs_Type, "Turbulent_Premixed_Flame") == 0) {
            i_Original_ICs = IC_TURBULENT_PREMIXED_FLAME;
        } else if (strcmp(Original_ICs_Type, "Restart") == 0) {
            i_Original_ICs = IC_RESTART;
        } else {
            i_command = INVALID_INPUT_VALUE;
        } /* endif */
        

    } else if (strcmp(code, "Gas_Type") == 0) {
       i_command = 101;
       value_stream >> value_string;
       strcpy(Gas_Type, value_string.c_str());

    } else if (strcmp(code, "Mach_Number") == 0) {
       i_command = 102;
       value_stream >> Mach_Number;
       if (Mach_Number < ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } /* endif */

    } else if (strcmp(code, "Reynolds_Number") == 0) {
       i_command = 103;
       value_stream >> Reynolds_Number;
       if (Reynolds_Number < ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } /* endif */

    } else if (strcmp(code, "Pressure") == 0) {
       i_command = 104;
       value_stream >> Pressure;
       Pressure = Pressure*THOUSAND;
       if (Pressure <= ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } /* endif */

    } else if (strcmp(code, "Temperature") == 0) {
       i_command = 105;
       value_stream >> Temperature;
       if (Temperature <= ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } /* endif */
            
    } else if (strcmp(code, "Flow_Angle") == 0) {
       i_command = 106;
       value_stream >> Flow_Angle;

    //
    // Other initial and boundary condition related input parameters:
    //
    } else if (strcmp(code, "Moving_Wall_Velocity") == 0) {
       i_command = 110;
       value_stream >> Moving_Wall_Velocity;

    } else if (strcmp(code, "Pressure_Gradient") == 0) {
       i_command = 111;
       value_stream >> Pressure_Gradient;

    } else if (strcmp(code, "Mean_Velocity") == 0) {
       i_command = 112;
       value_stream >> Mean_Velocity;

    } else if (strcmp(code, "Fresh_Gas_Height") == 0) {
       i_command = 113;
       value_stream >> Fresh_Gas_Height;

    //
    // Invalid input parameter code:
    //
    } else {
       i_command = INVALID_INPUT_CODE;

    } /* endif */

    // Parse other input parameters.      
    if (i_command == INVALID_INPUT_CODE) {
       // Multigrid
       i_command = Multigrid_IP.Parse_Next_Input_Control_Parameter(code, 
                                                                   value_stream);

       // NKS
       if (i_command == INVALID_INPUT_CODE) {
          i_command = NKS_IP.Parse_Next_Input_Control_Parameter(code, 
                                                                value_stream);
       } /* endif */

       // Grid
       if (i_command == INVALID_INPUT_CODE) {
          i_command = Grid_IP.Parse_Next_Input_Control_Parameter(code, 
                                                                 value_stream);
       } /* endif */

       // AMR
       if (i_command == INVALID_INPUT_CODE) {
          i_command = AMR_IP.Parse_Next_Input_Control_Parameter(code, 
                                                                value_stream);
       } /* endif */

       // Multi-Species
       if (i_command == INVALID_INPUT_CODE) {
          i_command = Species_IP.Parse_Next_Input_Control_Parameter(code, 
                                                                    value_stream);
       } /* endif */

       // Turbulence modelling
       if (i_command == INVALID_INPUT_CODE) {
          i_command = Turbulence_IP.Parse_Next_Input_Control_Parameter(code, 
                                                                       value_stream);
       } /* endif */

       // High Order
       if (i_command == INVALID_INPUT_CODE) {
	 i_command = HighOrder_IP.Parse_Next_Input_Control_Parameter(code,
								     value_stream);
       } /* endif */
       
    } /* endif */

    // Check for long input lines that have not been completely read in.
    if (!Input_File.good()) { 
       i_command = INVALID_INPUT_VALUE;
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
int CFD_Input_Parameters::Process_Input_Control_Parameter_File(char *Input_File_Name_ptr,
                                                               int &Command_Flag) {
   
   int error_flag, line_number;
   
   /* Assign initial value for error indicator flag. */
   error_flag = 0;
   
   /* Copy input file name (a string) to appropriate input parameter variable. */
   if (Input_File_Name_ptr != NULL) strcpy(Input_File_Name, Input_File_Name_ptr);
   
   /* Open the input file containing the input parameters. */
   Open_Input_File();
   error_flag = Input_File.fail();

   if (error_flag) {
      cout << "\n ERROR: Unable to open input data file.\n";
      return (error_flag);
   } /* endif */
    
   /* Read and parse control parameters contained in
      the input file. */
   while (1) {
      error_flag = Get_Next_Input_Control_Parameter(true);
      if (error_flag) {
         line_number = -line_number;
         cout << "\n ERROR: Error reading data at line #"
              << -line_number  << " of input data file.\n";
         error_flag = line_number;
         return (error_flag);
      } /* endif */
     
      Command_Flag = Parse_Next_Input_Control_Parameter();

      line_number = Line_Number;

      if (Command_Flag == EXECUTE_CODE) {
         break;
      } else if (Command_Flag == TERMINATE_CODE) {
         break;
      } else if (Command_Flag == INVALID_INPUT_CODE ||
                 Command_Flag == INVALID_INPUT_VALUE) {
         line_number = -line_number;
         cout << "\n ERROR: Error reading data at line #"
              << -line_number  << " of input data file.\n";
         error_flag = line_number;
         break;
      } /* endif */
   } /* endwhile */

   /* Perform consistency checks on all input parameters. */

   if (!error_flag) {
      error_flag = Check_Inputs();
      if (error_flag) {
	 cout << "\n ERROR: Input parameters consistency check failed.\n";
	 return (error_flag);
      } /* endif */
   } /* endif */

   /* Processing of input control parameters complete.  
      Return the error indicator flag. */
    
   return (error_flag);

}

/********************************************************
 * Routine: Broadcast                                   *
 *                                                      *
 * Broadcast the input parameters variables to all      *
 * processors involved in the calculation from the      *
 * primary processor using the MPI broadcast routine.   *
 *                                                      *
*********************************************************/
void CFD_Input_Parameters::Broadcast(void) {
   
#ifdef _MPI_VERSION
    // Input file parameters:
    MPI::COMM_WORLD.Bcast(CFFC_Path,
                          INPUT_PARAMETER_LENGTH,
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(Input_File_Name,
                          INPUT_PARAMETER_LENGTH,
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(Line_Number),
                          1,
                          MPI::INT, 0);

    // Output parameters:
    MPI::COMM_WORLD.Bcast(Output_File_Name,
                          INPUT_PARAMETER_LENGTH,
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(Restart_File_Name,
                          INPUT_PARAMETER_LENGTH,
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(Output_Format_Type,
                          INPUT_PARAMETER_LENGTH,
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(i_Output_Format),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Restart_Solution_Save_Frequency),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Time_Accurate_Output_Frequency),
                          1,
                          MPI::INT, 0);

    // Debugging parameters:
    MPI::COMM_WORLD.Bcast(&(Debug_Level),
                          1,
			  MPI::INT, 0);

    // Flow type indicator and related input parameters:
    MPI::COMM_WORLD.Bcast(Flow_Type,
                          INPUT_PARAMETER_LENGTH,
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(i_Flow_Type),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Axisymmetric),
                          1,
                          MPI::INT, 0);

    // Time integration type indicator and related input parameters:
    MPI::COMM_WORLD.Bcast(Time_Integration_Type,
                          INPUT_PARAMETER_LENGTH,
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(i_Time_Integration),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Time_Accurate),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Local_Time_Stepping),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Maximum_Number_of_Time_Steps),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(N_Stage),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(CFL_Number),
                           1,
                           MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Time_Max),
                          1,
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Residual_Norm),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Number_of_Residual_Norms),
                          1,
                          MPI::INT, 0);
    // Reconstruction type indicator and related input parameters:
    MPI::COMM_WORLD.Bcast(Reconstruction_Type,
                          INPUT_PARAMETER_LENGTH,
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(i_Reconstruction),
                          1,
                          MPI::INT, 0);

    // Limiter type indicator and related input parameters:
    MPI::COMM_WORLD.Bcast(Limiter_Type,
                          INPUT_PARAMETER_LENGTH,
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(i_Limiter),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Freeze_Limiter),
			  1,
			  MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Freeze_Limiter_Residual_Level),
                          1,
                          MPI::DOUBLE, 0);

    // Inviscid flux function type and related input parameters:
    MPI::COMM_WORLD.Bcast(Flux_Function_Type,
                          INPUT_PARAMETER_LENGTH,
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(i_Flux_Function),
                          1,
                          MPI::INT, 0);

    // Preconditioner indicator and related input parameters:
    MPI::COMM_WORLD.Bcast(&(Preconditioning),
                          1,
			  MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Mach_Number_Reference),
			  1,
			  MPI::DOUBLE, 0);

    // Implicit residual smoothing control parameters:
    MPI::COMM_WORLD.Bcast(&(Residual_Smoothing),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Residual_Smoothing_Epsilon),
                          1,
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Residual_Smoothing_Gauss_Seidel_Iterations),
                          1,
                          MPI::INT, 0);

    // Initial condition type indicator and related input parameters:
    MPI::COMM_WORLD.Bcast(ICs_Type,
                          INPUT_PARAMETER_LENGTH,
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(i_ICs),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(Original_ICs_Type,
                          INPUT_PARAMETER_LENGTH,
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(i_Original_ICs),
                          1,
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(Gas_Type, 
                          INPUT_PARAMETER_LENGTH, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(Mach_Number),
			  1,
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Reynolds_Number),
                          1,
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Pressure),
                          1,
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Temperature),
                          1,
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Flow_Angle),
                          1,
                          MPI::DOUBLE, 0);

    // Other initial and boundary condition related input parameters:
    MPI::COMM_WORLD.Bcast(&(Moving_Wall_Velocity.x),
                          1,
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Moving_Wall_Velocity.y),
                          1,
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Moving_Wall_Velocity.z),
                          1,
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Mean_Velocity.x),
                          1,
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Mean_Velocity.y),
                          1,
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Mean_Velocity.z),
                          1,
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Fresh_Gas_Height),
                          1,
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Pressure_Gradient.x),
                          1,
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Pressure_Gradient.y),
                          1,
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Pressure_Gradient.z),
                          1,
			  MPI::DOUBLE, 0);

    // Multigrid Input Parameters:
    Multigrid_IP.Broadcast();
    
    // NKS Input Parameters:
    NKS_IP.Broadcast();

    // Grid Input Parameters:
    Grid_IP.Broadcast();

    // AMR Input Parameters:
    AMR_IP.Broadcast();

    // Multi-Species Input Parameters:
    Species_IP.Broadcast();

    // Turbulence modelling
    Turbulence_IP.Broadcast();

    // High Order Grid Execution Mode:
    // Grid_HO_Execution_Mode.Broadcast(); --> RR: broadcast later

    // CENO Execution Mode:
    // CENO_Execution_Mode.Broadcast();  --> RR: broadcast later

#endif

}

/********************************************************
 * Routine: Check_Inputs                                *
 *                                                      *
 * Checks the validity and consistency of all input     *
 * parameters.                                          *
 *                                                      *
*********************************************************/
int CFD_Input_Parameters::Check_Inputs(void) {

    int error_flag = 0;

    // CFD Input Parameters:
    // ---------------------

    // Make sure that the spatial order of accuracy is in agreement
    // with a limiter type of ZERO for the piecewise constant case
    if (i_Limiter == LIMITER_ZERO && Reconstruction_Order == 1){
      Reconstruction_Order = 0;
      Grid3D_HO_Execution_Mode::RECONSTRUCTION_ORDER = Reconstruction_Order;
      cout << "\n CFD::Check_Inputs: Note: Spatial Order of Accuracy has been reduced from second to first-order due to the Limiter_Type having been set to Zero." << endl;
      cout.flush();
    } else if (Reconstruction_Order == 0 && i_Limiter != LIMITER_ZERO){
      cout << "CFD::Check_Inputs: Limiter_Type must be set to zero for a Spatial_Order_of_Accuracy = 1." << endl;
      cout.flush();
      return 1;
    } else if (i_Limiter == LIMITER_ZERO && Reconstruction_Order > 1){
      cout << "CFD::Check_Inputs: Caution! Limiter_Type has been set to zero which conflicts with the desired Spatial_Order_of_Accracy." <<endl;
      cout.flush();
      return 1;
    }/* endif */

    // Make sure that the reconstruction type for high-order reconstruction is set to CENO
    if (i_Reconstruction != RECONSTRUCTION_HIGH_ORDER && Reconstruction_Order > 1){
      cout << "CFD::Check_Inputs: Cannot perform high-order reconstruction with the given Reconstruction_Type.\n";
      cout << "For high-order reconstruction please use the CENO Reconstruction_Type." << endl;
      cout.flush();
      return 1;
    }

    // Multigrid Input Parameters:
    if (!error_flag &&
        i_Time_Integration == TIME_STEPPING_MULTIGRID) {
       error_flag = Multigrid_IP.Check_Inputs();
    } /* endif */
    
    // NKS Input Parameters:
    if (!error_flag) {
       error_flag = NKS_IP.Check_Inputs();
    } /* endif */

    // Grid Input Parameters:
    if (!error_flag) {
       error_flag = Grid_IP.Check_Inputs();
    } /* endif */

    // AMR Input Parameters:
    if (!error_flag) {
       error_flag = AMR_IP.Check_Inputs();
    } /* endif */

    // Multi-Species Input Parameters:
    if (!error_flag) {
       error_flag = Species_IP.Check_Inputs();
    } /* endif */

    // Turbulence modelling
    if (!error_flag) {
       error_flag = Turbulence_IP.Check_Inputs();
    } /* endif */

    // Input parameters are consistent.  Exit successfully.
    return (error_flag);

}

/*************************************************************
 * Input_Parameters -- Input-output operators.               *
 *************************************************************/

ostream &operator << (ostream &out_file,
                      const CFD_Input_Parameters &IP) {

   IP.Output(out_file);   
   return (out_file);

}

istream  &operator >> (istream &in_file,
                       CFD_Input_Parameters &IP) {

   return (in_file);

}

void CFD_Input_Parameters::Output(ostream &out_file) const {

   out_file << setprecision(6);
   Output_Problem_Type(out_file);
   Output_Solver_Type(out_file);
   Output_ICsBCs_Types(out_file);
   Output_Solution_Type(out_file);
   Output_IO_Types(out_file);
   Output_GridAMR_Types(out_file);

}

void CFD_Input_Parameters::Output_Problem_Type(ostream &out_file) const {

   out_file << "\n  -> CFFC Path: " << CFFC_Path;

   out_file << "\n\n Solving 3D ";
   if (i_Flow_Type ==  FLOWTYPE_INVISCID){
      out_file<<"Euler (Inviscid) ";
   } else {
      out_file<<"Navier-Stokes (Viscous) ";
   } /* endif */
   out_file << "Equations (IBVP/BVP)"; 
   if (i_Flow_Type ==  FLOWTYPE_INVISCID) {
   } else if (i_Flow_Type == FLOWTYPE_LAMINAR) {
      out_file << "\n  -> Laminar flow";
   } else if (i_Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
      out_file << "\n  -> Turbulent flow: RANS with k-epsilon turbulence model";
   } else if (i_Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      out_file << "\n  -> Turbulent flow: RANS with k-oemga turbulence model";
   } else if (i_Flow_Type == FLOWTYPE_TURBULENT_LES) {
      out_file << "\n  -> Turbulent flow: LES ";
   } else if (i_Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K) {
      out_file << "\n  -> Turbulent flow: LES with C-Fsd-k model";
   } else if (i_Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY) {
      out_file << "\n  -> Turbulent flow: LES with C-Fsd-Smagorinsky model";
   } else if (i_Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {
      out_file << "\n  -> Turbulent flow: LES with TF-k model";
   } else if (i_Flow_Type == FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY) {
      out_file << "\n  -> Turbulent flow: LES with TF-Smagorinsky model";	 
   } else if (i_Flow_Type == FLOWTYPE_TURBULENT_DES_K_OMEGA) {
      out_file << "\n  -> Turbulent flow: DES with k-omega SGS turbulence model ";
   } else if (i_Flow_Type == FLOWTYPE_TURBULENT_DNS) {
      out_file << "\n  -> Turbulent flow: DNS ";
   } /* endif */

   // Turbulence modelling
   if (i_Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON ||
       i_Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      out_file << "\n\n Turbulence Modelling";
      out_file << Turbulence_IP;
   } /* endif */

   if (Time_Accurate) { 
      out_file << "\n  -> Time Accurate (Unsteady) Solution";
   } else {
      out_file << "\n  -> Time Invariant (Steady-State) Solution";
   } /* endif */
   
   if (Debug_Level) {
      out_file << "\n  -> Debug level: "
	       << Debug_Level;
   } /* endif */

}

void CFD_Input_Parameters::Output_Solver_Type(ostream &out_file) const {

   out_file << "\n\n Time Marching Scheme";
   out_file << "\n  -> Time Integration: " 
            << Time_Integration_Type;
   out_file << "\n  -> Number of Stages in Multi-Stage Scheme: " 
            << N_Stage;
   if (Local_Time_Stepping == GLOBAL_TIME_STEPPING) {
      out_file << "\n  -> Global Time Stepping";
   } else if (Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
      out_file << "\n  -> Scalar Local Time Stepping";
   } else if (Local_Time_Stepping == MATRIX_LOCAL_TIME_STEPPING) {
      out_file << "\n  -> Matrix Local Time Stepping";
   } else if (Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) {
      out_file << "\n  -> Low-Mach-Number Local Preconditioning (Weiss-Smith)";
   } else if (Local_Time_Stepping == SEMI_IMPLICIT_LOCAL_TIME_STEPPING) {
      out_file << "\n  -> Semi-Implicit Local Time Stepping";
   } else if (Local_Time_Stepping == SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER) {
      out_file << "\n  -> Semi-Implicit Low-Mach-Number Local Preconditioned Time Stepping";
   } /* endif */
   out_file << "\n  -> CFL Number: " 
            << CFL_Number;
   if (Preconditioning == 1) {
      out_file <<"\n  -> Mach Number Reference "<< Mach_Number_Reference;
   } /* endif */
   
   out_file << "\n\n Spatial Discretization";
   out_file << "\n  -> Reconstruction Type: " 
            << Reconstruction_Type;
   out_file << "\n  -> Order of Accuracy: "
	    << Reconstruction_Order + 1;
   out_file << "\n  -> Limiter: " 
            << Limiter_Type;   
   if (Limiter_Type != LIMITER_ZERO && Freeze_Limiter) {
      out_file << "\n  -> Freeze Limiter when L2-norm of residual is < "
	       << Freeze_Limiter_Residual_Level;
   } /* endif */

   out_file << "\n  -> Flux Function: " 
            << Flux_Function_Type;
   
   if (Residual_Smoothing) {
      out_file << "\n  -> Residual Smoothing";
      out_file << "\n  -> Epsilon: " 
               << Residual_Smoothing_Epsilon;
      out_file << "\n  -> Gauss_Seidel_Iterations: " 
               << Residual_Smoothing_Gauss_Seidel_Iterations;
   } /* endif */

   out_file << "\n  -> Residual_Norm : " << Residual_Norm;

}

void CFD_Input_Parameters::Output_ICsBCs_Types(ostream &out_file) const {

   out_file << "\n\n Initial and Boundary Data";
   out_file << "\n  -> Initial Conditions: " 
            << ICs_Type;

   out_file << "\n\n Free Stream Conditions"; 
   out_file << "\n  -> Mach Number: " 
            << Mach_Number;
   out_file << "\n  -> Reynolds Number: " 
            << Reynolds_Number;
   out_file << "\n  -> Pressure (kPa): " 
            << Pressure/THOUSAND;
   out_file << "\n  -> Temperature (K): " 
            << Temperature;
   out_file << "\n  -> Flow Angle: " 
            << Flow_Angle;

   if (i_Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY ||
       i_Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K ) {
      out_file << "\n\n Laminar Flame Parameters";
      out_file << "\n  -> Fuel Equivalence Ratio: " 
               << Turbulence_IP.Fuel_Equivalence_Ratio;
      out_file << "\n  -> Unburnt Fuel Mass Fraction: " 
               << Turbulence_IP.Unburnt_Fuel_Mass_Fraction;
      out_file << "\n  -> Reactants Density: " 
               << Turbulence_IP.Reactants_Density;
      out_file << "\n  -> Laminar Flame Speed: " 
               << Turbulence_IP.Laminar_Flame_Speed;
      out_file << "\n  -> Laminar Flame Thickness: " 
               << Turbulence_IP.Laminar_Flame_Thickness;
      out_file << "\n  -> Adiabatic Flame Temperature: " 
               << Turbulence_IP.Adiabatic_Flame_Temperature;
      out_file << "\n  -> Filter Width: " 
               << Turbulence_IP.Filter_Width;
   }
   if (i_ICs == IC_TURBULENT_BUNSEN_FLAME) {
      out_file << "\n\n Flame Parameters";
      out_file << "\n  -> Mean Flow Velocity -- X Direction: " 
               << Mean_Velocity.x;
      out_file << "\n  -> Mean Flow Velocity -- Y Direction: " 
               << Mean_Velocity.y;
      out_file << "\n  -> Mean Flow Velocity -- Z Direction: " 
               << Mean_Velocity.z;
      out_file << "\n  -> Initial Flame Height: " 
               << Fresh_Gas_Height;
   }
}

void CFD_Input_Parameters::Output_Solution_Type(ostream &out_file) const {

   out_file << "\n\n Polytropic Gas"; 
   out_file << "\n  -> Gas Type: " 
            << Gas_Type;

}

void CFD_Input_Parameters::Output_IO_Types(ostream &out_file) const {

   out_file << "\n\n Time Step Control and I/O";
   out_file << "\n  -> Maximum Time (ms): " 
            << Time_Max*THOUSAND;
   out_file << "\n  -> Maximum Number of Time Steps (Iterations): " 
            << Maximum_Number_of_Time_Steps;
   out_file << "\n  -> Input File Name: " 
            << Input_File_Name;
   out_file << "\n  -> Output File Name: " 
            << Output_File_Name;
   out_file << "\n  -> Output Format: " 
            << Output_Format_Type;
   out_file << "\n  -> Restart Solution Save Frequency: "
            << Restart_Solution_Save_Frequency
            << " steps (iterations)"; 
   if (Time_Accurate_Output_Frequency !=0 && Time_Accurate){
      out_file << "\n  -> Time Accurate Solution Plot Frequency: "
	       << Time_Accurate_Output_Frequency
	       << " steps (iterations)"; 
   } /* endif */

}

void CFD_Input_Parameters::Output_GridAMR_Types(ostream &out_file) const {

    // Grid Input Parameters:
   out_file << "\n\n Using Hexhedral Multi-Block Grid";
   out_file << Grid_IP;

   // AMR Input Parameters:
   out_file << "\n\n Parallel Block-Based Adaptive Mesh Refinement (AMR)";
   out_file << AMR_IP;

}
