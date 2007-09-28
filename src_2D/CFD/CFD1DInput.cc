/*!file CFD1DInput.cc
  \brief Implementation of subroutines prototyped in CFD1DInput.h. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "CFD1DInput.h"
#include "../Euler1D/ExactSolutions/ExactSolutions.h"

/*************************************************************
 * CFD1D_Input_Parameters -- External subroutines.           *
 *************************************************************/

/********************************************************
 * Routine: Open_Input_File                             *
 *                                                      *
 * Opens the appropriate input data file.               *
 *                                                      *
 ********************************************************/
void Open_Input_File(CFD1D_Input_Parameters &IP) {

    IP.Input_File.open(IP.Input_File_Name, ios::in);
    if (! IP.Input_File.bad()) {
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
void Close_Input_File(CFD1D_Input_Parameters &IP) {

    IP.Input_File.unsetf(ios::skipws);
    IP.Input_File.close();

}

/********************************************************
 * Routine: Set_Default_Input_Parameters                *
 *                                                      *
 * Assigns default values to the input parameters.      *
 *                                                      *
 ********************************************************/
void Set_Default_Input_Parameters(CFD1D_Input_Parameters &IP) {

    char *string_ptr;

    string_ptr = "CFD1D.in";
    strcpy(IP.Input_File_Name, string_ptr);

    string_ptr = "Explicit_Euler";
    strcpy(IP.Time_Integration_Type, string_ptr);
    IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
    IP.Time_Accurate = 1;
    IP.Local_Time_Stepping = 0;
    IP.Maximum_Number_of_Time_Steps = 100;
    IP.N_Stage = 1;
    IP.Time_Max = ZERO;

    string_ptr = "MUSCL";
    strcpy(IP.Reconstruction_Type, string_ptr);
    IP.i_Reconstruction = RECONSTRUCTION_MUSCL;
    IP.i_ReconstructionMethod = RECONSTRUCTION_MUSCL;
    IP.Space_Accuracy = 2;
    IP.CENO_Cutoff = 30;
    IP.ExactFunction = NULL;

    string_ptr = "VanLeer";
    strcpy(IP.Limiter_Type, string_ptr);
    IP.i_Limiter = LIMITER_VANLEER;

    string_ptr = "Roe";
    strcpy(IP.Flux_Function_Type, string_ptr);
    IP.i_Flux_Function = FLUX_FUNCTION_ROE;

    string_ptr = "Uniform";
    strcpy(IP.ICs_Type, string_ptr);
    IP.i_ICs = IC_UNIFORM;

    IP.Kappa = 0.02;
    IP.a = ONE;
    IP.Tau = ONE;

    string_ptr = "Uniform";
    strcpy(IP.Grid_Type, string_ptr);
    IP.i_Grid = GRID_CARTESIAN_UNIFORM;
    IP.X_Min = -ONE;
    IP.X_Max = ONE;
    IP.Number_of_Cells = 50;
    IP.Number_of_Nodes = 51;

    string_ptr = "outputfile.dat";
    strcpy(IP.Output_File_Name, string_ptr);

    string_ptr = "gnuplotfile.gplt";
    strcpy(IP.Gnuplot_File_Name, string_ptr);

    string_ptr = "Gnuplot";
    strcpy(IP.Output_Format_Type, string_ptr);
    IP.i_Output_Format = IO_GNUPLOT;

    string_ptr = " ";
    strcpy(IP.Next_Control_Parameter, string_ptr);

    IP.Line_Number = 0;

    // set output to verbose
    IP.Verbose() = ON;
}

/********************************************************
 * Routine: Get_Next_Input_Control_Parameter            *
 *                                                      *
 * Get the next input control parameter from the input  *
 * file.                                                *
 *                                                      *
 ********************************************************/
void Get_Next_Input_Control_Parameter(CFD1D_Input_Parameters &IP) {

  int i, index, LineSize, IndexFirstChar(0);
  char buffer[256], ControlParameter[256];

  // Initialize ControlParameter and IP.Next_Control_Parameter to end of string
  ControlParameter[0] = '\0';
  strcpy(IP.Next_Control_Parameter, ControlParameter);

  // While the input stream is 'good' for reading and the end of file is not attained
  while ( IP.Input_File.good() && !IP.Input_File.getline(buffer, sizeof(buffer)).eof() ){

    // Process the line 
    IP.Line_Number = IP.Line_Number + 1;
    LineSize = IP.Input_File.gcount(); // Get the size of the line. Last character is "\0"!

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
	if (buffer[i] == ' ' || buffer[i] == '='){
	  ControlParameter[index] = '\0';
	  break;
	} else {
	  ControlParameter[index] = buffer[i];
	}
      }

      // Set the Next_Control_Parameter
      strcpy(IP.Next_Control_Parameter, ControlParameter);
      break;
    }

  }//endwhile

}


/********************************************************
 * Routine: Parse_Next_Input_Control_Parameter          *
 *                                                      *
 * Parses and executes the next input control parameter *
 * from the input file.                                 *
 *                                                      *
 ********************************************************/
int Parse_Next_Input_Control_Parameter(CFD1D_Input_Parameters &IP) {

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
       } else if (strcmp(IP.Time_Integration_Type, "Implicit_Euler") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_IMPLICIT_EULER;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "Explicit_Euler_HO") == 0){
           IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER_HIGH_ORDER;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "Implicit_Trapezoidal") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_IMPLICIT_TRAPEZOIDAL;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "Explicit_Predictor_Corrector") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR;
           IP.N_Stage = 2;
       } else if (strcmp(IP.Time_Integration_Type, "Explicit_Predictor_Corrector_HO") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR_HIGH_ORDER;
           IP.N_Stage = 2;
       } else if (strcmp(IP.Time_Integration_Type, "Explicit_Predictor_Corrector4_HO") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_RUNGE_KUTTA_4_HIGH_ORDER;
           IP.N_Stage = 4;
       } else if (strcmp(IP.Time_Integration_Type, "Semi_Implicit_Euler") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_SEMI_IMPLICIT_EULER;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "Semi_Implicit_Trapezoidal") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_SEMI_IMPLICIT_TRAPEZOIDAL;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "Semi_Implicit_Predictor_Corrector") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_SEMI_IMPLICIT_PREDICTOR_CORRECTOR;
           IP.N_Stage = 2;
       } else if (strcmp(IP.Time_Integration_Type, "Multistage_Optimal_Smoothing") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING;
           IP.N_Stage = 4;
       } else if (strcmp(IP.Time_Integration_Type, "Lax_Friedrichs") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_LAX_FRIEDRICHS;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "Lax_Wendroff") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_LAX_WENDROFF;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "MacCormack") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_MACCORMACK;
           IP.N_Stage = 2;
       } else if (strcmp(IP.Time_Integration_Type, "Hancock") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_HANCOCK;
           IP.N_Stage = 2;
       } else if (strcmp(IP.Time_Integration_Type, "Simple_Explicit") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_SIMPLE_EXPLICIT;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "Simple_Implicit") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_SIMPLE_IMPLICIT;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "Crank_Nicolson") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_CRANK_NICOLSON;
           IP.N_Stage = 1;
       } else if (strcmp(IP.Time_Integration_Type, "ADE") == 0) {
           IP.i_Time_Integration = TIME_STEPPING_ADE;
           IP.N_Stage = 2;
       } else {
	 i_command = INVALID_INPUT_CODE;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Reconstruction_Type") == 0) {
       i_command = 2;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Reconstruction_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Reconstruction_Type, "MUSCL") == 0) {
          IP.i_Reconstruction = RECONSTRUCTION_MUSCL;
	  IP.i_ReconstructionMethod = RECONSTRUCTION_MUSCL;
       } else if (strcmp(IP.Reconstruction_Type, "Green_Gauss") == 0) {
          IP.i_Reconstruction = RECONSTRUCTION_GREEN_GAUSS;
          IP.i_ReconstructionMethod = RECONSTRUCTION_GREEN_GAUSS;
       } else if (strcmp(IP.Reconstruction_Type, "Least_Squares") == 0) {
          IP.i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;
          IP.i_ReconstructionMethod = RECONSTRUCTION_LEAST_SQUARES;
       } else if (strcmp(IP.Reconstruction_Type, "Characteristic") == 0) {
          IP.i_Reconstruction = RECONSTRUCTION_CHARACTERISTIC;
          IP.i_ReconstructionMethod = RECONSTRUCTION_CHARACTERISTIC;
       } else if (strcmp(IP.Reconstruction_Type, "ENO") == 0) {
          IP.i_Reconstruction = RECONSTRUCTION_HIGH_ORDER;
          IP.i_ReconstructionMethod = RECONSTRUCTION_ENO;
       } else if (strcmp(IP.Reconstruction_Type, "ENO_Characteristic") == 0) {
          IP.i_Reconstruction = RECONSTRUCTION_HIGH_ORDER;
          IP.i_ReconstructionMethod = RECONSTRUCTION_ENO_CHARACTERISTIC;
       } else if (strcmp(IP.Reconstruction_Type, "CENO") == 0) {
          IP.i_Reconstruction = RECONSTRUCTION_HIGH_ORDER;
          IP.i_ReconstructionMethod = RECONSTRUCTION_CENO;
       } else {
	 i_command = INVALID_INPUT_CODE;
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
       } else if (strcmp(IP.Limiter_Type, "Minmod") == 0) {
          IP.i_Limiter = LIMITER_MINMOD;
       } else if (strcmp(IP.Limiter_Type, "UMIST") == 0) {
          IP.i_Limiter = LIMITER_UMIST;
       } else if (strcmp(IP.Limiter_Type, "Double_Minmod") == 0) {
          IP.i_Limiter = LIMITER_DOUBLE_MINMOD;
       } else if (strcmp(IP.Limiter_Type, "Superbee") == 0) {
          IP.i_Limiter = LIMITER_SUPERBEE;
       } else if (strcmp(IP.Limiter_Type, "Phi") == 0) {
          IP.i_Limiter = LIMITER_PHI;
       } else if (strcmp(IP.Limiter_Type, "VanLeer") == 0) {
          IP.i_Limiter = LIMITER_VANLEER;
       } else if (strcmp(IP.Limiter_Type, "VanAlbada") == 0) {
          IP.i_Limiter = LIMITER_VANALBADA;
       } else if (strcmp(IP.Limiter_Type, "Barth_Jespersen") == 0) {
          IP.i_Limiter = LIMITER_BARTH_JESPERSEN;
       } else if (strcmp(IP.Limiter_Type, "Venkatakrishnan") == 0) {
          IP.i_Limiter = LIMITER_VENKATAKRISHNAN;
       } else {
	 i_command = INVALID_INPUT_CODE;
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
       } else if (strcmp(IP.Flux_Function_Type, "Osher") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_OSHER;
       } else {
	 i_command = INVALID_INPUT_CODE;
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
       } else if (strcmp(IP.ICs_Type, "Groth") == 0) {
          IP.i_ICs = IC_GROTH;
       } else if (strcmp(IP.ICs_Type, "Einfeldt") == 0) {
          IP.i_ICs = IC_EINFELDT;
       } else if (strcmp(IP.ICs_Type, "Shock_Wave") == 0) {
          IP.i_ICs = IC_SHOCK_WAVE;
       } else if (strcmp(IP.ICs_Type, "Contact_Surface") == 0) {
          IP.i_ICs = IC_CONTACT_SURFACE;
       } else if (strcmp(IP.ICs_Type, "Rarefaction_Wave") == 0) {
          IP.i_ICs = IC_RAREFACTION_WAVE;
       } else if (strcmp(IP.ICs_Type, "Brio_Wu") == 0) {
          IP.i_ICs = IC_BRIO_WU;
       } else if (strcmp(IP.ICs_Type, "Riemann") == 0) {
          IP.i_ICs = IC_RIEMANN;
       } else if (strcmp(IP.ICs_Type, "Square_Wave") == 0) {
          IP.i_ICs = IC_SQUARE_WAVE;
       } else if (strcmp(IP.ICs_Type, "Sinx2_Wave") == 0) {
          IP.i_ICs = IC_SINX2_WAVE;
       } else if (strcmp(IP.ICs_Type, "Impulsive_Rod") == 0) {
          IP.i_ICs = IC_IMPULSIVE_ROD;
          IP.X_Min = ZERO;
          IP.X_Max = TEN;
       } else if (strcmp(IP.ICs_Type, "Sinusoidal_Rod1") == 0) {
          IP.i_ICs = IC_SINUSOIDAL_ROD1;
          IP.X_Min = ZERO;
          IP.X_Max = TEN;
       } else if (strcmp(IP.ICs_Type, "Sinusoidal_Rod4") == 0) {
          IP.i_ICs = IC_SINUSOIDAL_ROD4;
          IP.X_Min = ZERO;
          IP.X_Max = TEN;
       } else if (strcmp(IP.ICs_Type, "Riemann_IVP_qx=0") == 0) {
          IP.i_ICs = IC_RIEMANN_IVP_QX0;
          IP.X_Min = -ONE;
          IP.X_Max = ONE;
       } else if (strcmp(IP.ICs_Type, "Riemann_IVP_T=0") == 0) {
          IP.i_ICs = IC_RIEMANN_IVP_T0;
          IP.X_Min = -ONE;
          IP.X_Max = ONE;
       } else if (strcmp(IP.ICs_Type, "Riemann_IVP") == 0) {
          IP.i_ICs = IC_RIEMANN_IVP;
          IP.X_Min = -ONE;
          IP.X_Max = ONE;
       } else if (strcmp(IP.ICs_Type, "Sin_IVP") == 0) {
          IP.i_ICs = IC_SIN_WAVE;
	  IP.ExactFunction = SIN_WAVE_Solution;
       } else if (strcmp(IP.ICs_Type, "Jiang_IVP") == 0) {
          IP.i_ICs = IC_JIANG_WAVE;
	  IP.ExactFunction = JIANG_IVP_Solution;
       } else if (strcmp(IP.ICs_Type, "Acoustic_Shock_IVP") == 0) {
          IP.i_ICs = IC_SHOCK_ACOUSTIC_INTERACTION;
       } else if (strcmp(IP.ICs_Type, "Blast_Wave_IVP") == 0) {
          IP.i_ICs = IC_BLAST_WAVE_INTERACTION;
       } else if (strcmp(IP.ICs_Type, "Convection_Shapes_IVP") == 0) {
          IP.i_ICs = IC_CONVECTION_OF_DIFFERENT_SHAPES;
	  IP.ExactFunction = ConvectionShapes;
       } else if (strcmp(IP.ICs_Type, "Density_Step_IVP") == 0) {
          IP.i_ICs = IC_DENSITY_STEP_WAVE;
       } else {
	 i_command = INVALID_INPUT_CODE;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Grid_Type") == 0) {
       i_command = 6;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Grid_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Grid_Type, "Uniform") == 0) {
          IP.i_Grid = GRID_CARTESIAN_UNIFORM;
       } else {
	 i_command = INVALID_INPUT_CODE;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Output_File_Name") == 0) {
       i_command = 7;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Output_File_Name, 
              IP.Next_Control_Parameter);
       strcat(IP.Output_File_Name, ".dat");
       strcpy(IP.Gnuplot_File_Name, 
              IP.Next_Control_Parameter);
       strcat(IP.Gnuplot_File_Name, ".gplt");

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Cells") == 0) {
       i_command = 8;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Cells;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Cells < 1) i_command = INVALID_INPUT_VALUE;
       IP.Number_of_Nodes = IP.Number_of_Cells + 1;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Nodes") == 0) {
       i_command = 9;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Nodes;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Nodes < 2) i_command = INVALID_INPUT_VALUE;
       IP.Number_of_Cells = IP.Number_of_Nodes - 1;

    } else if (strcmp(IP.Next_Control_Parameter, "Time_Accurate") == 0) {
       i_command = 10;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Time_Accurate;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Time_Accurate != 0 &&
           IP.Time_Accurate != 1) IP.Time_Accurate = 0;
       if (IP.Time_Accurate) {
          IP.Local_Time_Stepping = 0;
       } else {
          IP.Local_Time_Stepping = 1;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Local_Time_Stepping") == 0) {
       i_command = 11;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Local_Time_Stepping;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Local_Time_Stepping != 0 &&
           IP.Local_Time_Stepping != 1) IP.Local_Time_Stepping = 1;

    } else if (strcmp(IP.Next_Control_Parameter, "Maximum_Number_of_Time_Steps") == 0) {
       i_command = 12;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Maximum_Number_of_Time_Steps;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Maximum_Number_of_Time_Steps < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "N_Stage") == 0) {
       i_command = 13;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.N_Stage;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.N_Stage < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "CFL_Number") == 0) {
       i_command = 14;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.CFL_Number;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.CFL_Number <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "X_Min") == 0) {
       i_command = 15;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.X_Min;
       IP.Input_File.getline(buffer, sizeof(buffer));
       IP.X_ExactSolution_Min = IP.X_Min; // assume that the exact solution is given on the real domain.
       // See ICs for domains that are different

    } else if (strcmp(IP.Next_Control_Parameter, "X_Max") == 0) {
       i_command = 16;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.X_Max;
       IP.Input_File.getline(buffer, sizeof(buffer));
       IP.X_ExactSolution_Max = IP.X_Max; // assume that the exact solution is given on the real domain.
       // See ICs for domains that are different

    } else if (strcmp(IP.Next_Control_Parameter, "Time_Max") == 0) {
       i_command = 17;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Time_Max;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Time_Max < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Kappa") == 0) {
       i_command = 18;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Kappa;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Kappa < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "a") == 0) {
       i_command = 19;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.a;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Tau") == 0) {
       i_command = 20;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Tau;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Tau < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Output_Format_Type") == 0) {
       i_command = 21;
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

    } else if (strcmp(IP.Next_Control_Parameter, "Space_Accuracy") == 0) {
       i_command = 22;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Space_Accuracy;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Space_Accuracy <= 0 && IP.Space_Accuracy >= 7){
	 IP.Space_Accuracy = 1;
	 cout << "\n Space Accuracy should be between 1 and 6 \n"
	      << "Space Accuracy set to 1" << endl;
       }/* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "CENO_Tolerance") == 0) {
       i_command = 24;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.FitTolerance();
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Execute") == 0) {
       i_command = EXECUTE_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Terminate") == 0) {
       i_command = TERMINATE_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Continue") == 0) {
       i_command = CONTINUE_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output") == 0) {
       i_command = WRITE_OUTPUT_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Accuracy") == 0) {
      i_command = WRITE_OUTPUT_ACCURACY_CODE;
      
    } else if (strcmp(IP.Next_Control_Parameter, "Print_Norms") == 0) {
      i_command = WRITE_NORM_ON_SCREEN;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Restart") == 0) {
       i_command = WRITE_RESTART_CODE;

    } else if (IP.Next_Control_Parameter[0] == '#') {
       i_command = COMMENT_CODE;

    } else {
       i_command = INVALID_INPUT_CODE;

    } /* endif */

    /* Return the parser command type indicator. */

    return (i_command);

}

