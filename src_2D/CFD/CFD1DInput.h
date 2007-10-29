/*!file CFD1DInput.h
  \brief Header file defining CFD1D input subroutines and macros. */

#ifndef _CFD1DINPUT_INCLUDED
#define _CFD1DINPUT_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "CFD.h"
#include "../HighOrderReconstruction/CENO_ExecutionMode.h"
#include "../HighOrderReconstruction/CENO_Tolerances.h"

#define	INPUT_PARAMETER_LENGTH_CFD1D    80
#define INPUT_PARAMETER_LENGTH_CFD1D_FILE_NAME 350

/********************************************************
 * Class:  CFD1D_Input_Parameters                       *
 ********************************************************/
class CFD1D_Input_Parameters{
  private:
  public:
  // Input file name:
  char Input_File_Name[INPUT_PARAMETER_LENGTH_CFD1D_FILE_NAME];

  // Input file stream:
  ifstream Input_File;

  // Time integration type indicator and related input parameters:
  char Time_Integration_Type[INPUT_PARAMETER_LENGTH_CFD1D];
  int i_Time_Integration;
  bool Reconstruction_In_Each_Stage;
  int Time_Accurate, Local_Time_Stepping, 
      Maximum_Number_of_Time_Steps, N_Stage;
  double CFL_Number, Time_Max;

  // Reconstruction type indicator and related input parameters:
  char Reconstruction_Type[INPUT_PARAMETER_LENGTH_CFD1D];
  int i_Reconstruction;
  int i_ReconstructionMethod;
  int Space_Accuracy;

  // Limiter type indicator and related input parameters:
  char Limiter_Type[INPUT_PARAMETER_LENGTH_CFD1D];
  int i_Limiter;

  // Flux_Function_Type and related input parameters
  char Flux_Function_Type[INPUT_PARAMETER_LENGTH_CFD1D];
  int i_Flux_Function;

  // Initial condition type indicator and related input parameters:
  char ICs_Type[INPUT_PARAMETER_LENGTH_CFD1D];
  int i_ICs;
  FunctionType1D ExactFunction; // pointer to the exact function

  // Grid type indicator and related input parameters:
  char Grid_Type[INPUT_PARAMETER_LENGTH_CFD1D];
  int i_Grid;
  int Number_of_Cells, Number_of_Nodes;
  double X_Min, X_Max;
  double X_ExactSolution_Min, X_ExactSolution_Max; /* determine the definition domain of the exact solution */

  // Diffusion coefficient, wave speed, and relaxation time:
  double Kappa, a, Tau;

  // Output file name:
  char Output_File_Name[INPUT_PARAMETER_LENGTH_CFD1D];

  // Gnuplot file name:
  char Gnuplot_File_Name[INPUT_PARAMETER_LENGTH_CFD1D];

  // Next_Control_Parameter:
  char Next_Control_Parameter[INPUT_PARAMETER_LENGTH_CFD1D];

  // Output format type indicator:
  char Output_Format_Type[INPUT_PARAMETER_LENGTH_CFD1D];
  int i_Output_Format;

  // Input file line number:
  int Line_Number;

  // Locally used variable
  char **charPtr;

  // Batch mode or verbose
  short verbose_flag;

  unsigned ErrorParameter; //!< Error assessment parameter (the parameter used to compute the solution error)

  /* Access fields */
  short & Verbose(void) {return verbose_flag;}
  const short & Verbose(void) const {return verbose_flag;}

  int ReconstructionOrder(void) const {return (Space_Accuracy-1);}
  int & Limiter(void) {return i_Limiter;}
  const int & Limiter(void) const {return i_Limiter;} 

  int Parse_Input_File(char *Input_File_Name_ptr);

  // Return the number of required ghost cells for the specified ReconstructionMethod
  int Nghost(void) const;  

  // Read the next input control parameter
  void Get_Next_Input_Control_Parameter(void);

  /* Input-output operators. */

  friend ostream &operator << (ostream &out_file,
		               const CFD1D_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
			       CFD1D_Input_Parameters &IP);

};

/**********************************************************************
 * CFD1D_Input_Parameters -- Input-output operators.         *
 **********************************************************************/
inline ostream &operator << (ostream &out_file,
			     const CFD1D_Input_Parameters &IP) {
    if (IP.i_Grid == GRID_CARTESIAN_UNIFORM) {
       out_file << "\n\n Solving PDE(s) on 1D uniform Cartesian mesh.";
    } else {
       out_file << "\n\n Solving PDE(s) on 1D mesh.";
    } /* endif */
    out_file << "\n  -> Input File Name: " 
             << IP.Input_File_Name;

    if (IP.Local_Time_Stepping) {
      out_file << "\n  -> Time Invariant (Steady-State) Solution";
      out_file << "\n  -> Local Time Stepping";
      out_file << "\n  -> Maximum Numer of Time Steps: " 
	       << IP.Maximum_Number_of_Time_Steps;
    } else if (IP.Time_Accurate) { 
      out_file << "\n  -> Time Accurate (Unsteady) Solution";
      out_file << "\n  -> Maximum Time: " 
	       << IP.Time_Max;
    } else {
      out_file << "\n  -> ERROR: Time accuracy not defined!";
      
    }

    out_file << "\n  -> Time Integration: " 
             << IP.Time_Integration_Type;
    out_file << "\n  -> Number of Stages in Multi-Stage Scheme: " 
             << IP.N_Stage;
    out_file << "\n  -> CFL Number: " 
             << IP.CFL_Number;

    out_file << "\n  -> Space Accuracy Solution: ";
    switch(IP.Space_Accuracy){
    case 1: 
      out_file << "1st order";
      break;
    case 2:
      out_file << "2nd order";
      break;
    case 3:
      out_file << "3rd order";
      break;
    case 4:
      out_file << "4th order";
      break;
    case 5:
      out_file << "5th order";
      break;
    case 6:
      out_file << "6th order";
      break;
    default:
      out_file << "higher than 6th order";
   }

    out_file << "\n  -> Reconstruction: " 
             << IP.Reconstruction_Type;
    if ( IP.i_Reconstruction != RECONSTRUCTION_HIGH_ORDER ){
      out_file << "\n  -> Limiter: " << IP.Limiter_Type;
    }
    if (IP.i_ReconstructionMethod == RECONSTRUCTION_CENO){
      CENO_Execution_Mode::Print_Info(out_file);
      CENO_Tolerances::Print_Info(out_file);
      out_file << "\n  -> Limiter: " << IP.Limiter_Type;
    }

    out_file << "\n  -> Flux Function: " 
             << IP.Flux_Function_Type;

    out_file << "\n  -> Initial Conditions: " 
             << IP.ICs_Type;
    switch(IP.i_ICs) {
      case IC_CONSTANT :
        break;
      case IC_UNIFORM :
        break;
      case IC_SOD :
        break;
      case IC_GROTH :
        break;
      case IC_EINFELDT :
        break;
      default:
        break;
    } /* endswitch */
    out_file << "\n  -> Diffusion Coefficient : " 
             << IP.Kappa;
    out_file << "\n  -> Wave Speed : " 
             << IP.a;
    out_file << "\n  -> Relaxation Time : " 
             << IP.Tau;

    out_file << "\n  -> Grid: " 
             << IP.Grid_Type;
    switch(IP.i_Grid) {
      case GRID_CARTESIAN_UNIFORM :
        out_file << "\n  -> Minimum X Location : " 
                 << IP.X_Min;
        out_file << "\n  -> Maximum X Location : " 
                 << IP.X_Max;
        out_file << "\n  -> Width of Solution Domain : " 
                 << IP.X_Max-IP.X_Min;
        break;
      default:
        out_file << "\n  -> Minimum X Location : " 
                 << IP.X_Min;
        out_file << "\n  -> Maximum X Location : " 
                 << IP.X_Max;
        out_file << "\n  -> Width of Solution Domain : " 
                 << IP.X_Max-IP.X_Min;
        break;
    } /* endswitch */
    out_file << "\n  -> Number of Cells: "
             << IP.Number_of_Cells;
    out_file << "\n  -> Number of Nodes: " 
             << IP.Number_of_Nodes;

    out_file << "\n  -> Output File Name: " 
             << IP.Output_File_Name;
    out_file << "\n  -> Output Format: " 
             << IP.Output_Format_Type;

    return (out_file);
}

inline istream &operator >> (istream &in_file,
			     CFD1D_Input_Parameters &IP) {
    return (in_file);
}

/**********************************************************************
 * CFD1D_Input_Parameters -- External subroutines.           *
 **********************************************************************/

extern void Open_Input_File(CFD1D_Input_Parameters &IP);

extern void Close_Input_File(CFD1D_Input_Parameters &IP);

extern void Set_Default_Input_Parameters(CFD1D_Input_Parameters &IP);

extern void Get_Next_Input_Control_Parameter(CFD1D_Input_Parameters &IP);

extern int Parse_Next_Input_Control_Parameter(CFD1D_Input_Parameters &IP);

#endif
