/* CFDInput.h:  Header file defining general CFD 
                solution input parameter class. */

#ifndef _CFD_INPUT_INCLUDED
#define _CFD_INPUT_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

/* Include required CFFC header files. */

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

#ifndef _GAS_CONSTANTS_INCLUDED
#include "../Physics/GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED

#ifndef _FASMULTIGRIDINPUT_INCLUDED
#include "../FASMultigrid/FASMultigridInput.h"
#endif // _FASMULTIGRIDINPUT_INCLUDED

#ifndef _NKSINPUT_INCLUDED
#include "../NewtonKrylovSchwarz/NKSInput.h"
#endif // _NKSINPUT_INCLUDED

#ifndef _GRID3D_INPUT_INCLUDED
#include "../Grid/Grid3DInput.h"
#endif // _GRID3D_INPUT_INCLUDED

#ifndef _AMRINPUT_INCLUDED
#include "../AMR/AMRInput.h"
#endif // _AMRINPUT_INCLUDED

#ifndef _SPECIESINPUT_INCLUDED
#include "../Physics/SpeciesInput.h"
#endif // _SPECIESINPUT_INCLUDED

#ifndef _TURBULENCEMODEL_INPUT_INCLUDED
#include "../TurbulenceModelling/TurbulenceModellingInput.h"
#endif // _TURBULENCEMODEL_INPUT_INCLUDED

/* Define the class. */

#define	INPUT_PARAMETER_LENGTH    256

// Enviroment flag for CFFC root directory path
#define PATHVAR "CFFC_Path"

/*!
 * Class: CFD_Input_Parameters
 *
 * @brief Input parameters for general CFFC CFD solutions.
 *
 * This class defines and handles the input variables related to
 * CFFC CFD solutions.
 *
 */
class CFD_Input_Parameters{
  private:
  public:
  //@{ @name Input file parameters:
  //! CFFC root directory path
  char CFFC_Path[INPUT_PARAMETER_LENGTH];
  //! Input file name
  char Input_File_Name[INPUT_PARAMETER_LENGTH];
  //! Input file stream
  ifstream Input_File;
  //! Input file line number
  int Line_Number;
  //! Next_Control_Parameter
  char Next_Control_Parameter[INPUT_PARAMETER_LENGTH];
  //@}

  //@{ @name Output parameters:
  //! Output file name
  char Output_File_Name[INPUT_PARAMETER_LENGTH];
  //! Restart file name
  char Restart_File_Name[INPUT_PARAMETER_LENGTH];
  //! Output format type indicator
  char Output_Format_Type[INPUT_PARAMETER_LENGTH];
  int i_Output_Format;
  //! Frequency for saving restart solutions
  int Restart_Solution_Save_Frequency;
  //! Frequency of outputting solution progress information
  int Output_Progress_Frequency;
  //! Frequency of outputting solution for time accurate calculations
  int Time_Accurate_Output_Frequency;
  //@}

  //@{ @name Debugging parameters:
  //! Debug level for calculation (0 for none, 1,2,3... level of verboseness)
  int Debug_Level;
  //@}

  //@{ @name Flow type indicator and related input parameters:
  char Flow_Type[INPUT_PARAMETER_LENGTH];
  int i_Flow_Type;
  //! Axisymmetric flow indicator (1=axisymmetric flow, 0=2D planar or 3D flow)
  int Axisymmetric;
  //@}

  //@{ @name Time integration type indicator and related input parameters:
  char Time_Integration_Type[INPUT_PARAMETER_LENGTH];
  int i_Time_Integration;
  int Time_Accurate, Local_Time_Stepping, 
      Maximum_Number_of_Time_Steps, N_Stage;
  //! p-norm calculation indicator (0 density, 1, momentum, 2 k,...)
  int p_Norm_Indicator;
  int Number_of_Residual_Norms;
  double CFL_Number, Time_Max;
  //@}

  //@{ @name Reconstruction type indicator and related input parameters:
  char Reconstruction_Type[INPUT_PARAMETER_LENGTH];
  int i_Reconstruction;
  //@}

  //@{ @name Limiter type indicator and related input parameters:
  char Limiter_Type[INPUT_PARAMETER_LENGTH];
  int i_Limiter;
  int Freeze_Limiter;
  double Freeze_Limiter_Residual_Level;
  //@}

  //@{ @name Inviscid flux function type and related input parameters:
  char Flux_Function_Type[INPUT_PARAMETER_LENGTH];
  int i_Flux_Function;
  //@}
  
  //@{ @name Preconditioner type and related input parameters:
  int Preconditioning;
  double Mach_Number_Reference;
  //@}

  //@{ @name Implicit residual smoothing control parameters:
  int Residual_Smoothing;
  double Residual_Smoothing_Epsilon;
  int Residual_Smoothing_Gauss_Seidel_Iterations;
  //@}

  //@{ @name Multigrid related input parameters:
  Multigrid_Input_Parameters Multigrid_IP;
  //@}

  //@{ @name Newton-Krylov-Schwarz related input parameters:
  NKS_Input_Parameters NKS_IP;
  //@}

  //@{ @name Computational grid/mesh related input parameters:
  Grid3D_Input_Parameters Grid_IP;
  //@}

  //@{ @name AMR related input parameters:
  AMR_Input_Parameters AMR_IP;
  //@}

  //@{ @name Multi-species gaseous mixture related input parameters:
  Species_Input_Parameters Species_IP;
  //@}

  //@{ @name Turbulence modelling related input parameters:
  Turbulence_Modelling_Input_Parameters Turbulence_IP;
  //@}

  //@{ @name Initial condition type indicator and related input parameters:
  char ICs_Type[INPUT_PARAMETER_LENGTH];
  int i_ICs;
  //! Gas Type
  char Gas_Type[INPUT_PARAMETER_LENGTH];
  //! Freestream Mach number
  double Mach_Number;
  //! Freestream Reynolds number
  double Reynolds_Number;
  //! Freestream pressure
  double Pressure;
  //! Freestream temperature
  double Temperature;
  //! Freestream flow angle
  double Flow_Angle;
  //@}

  //@{ @name Other initial and boundary condition related input parameters:
  //! Wall velocity 
  Vector3D Moving_Wall_Velocity;
  //! Pressure gradient 
  Vector3D Pressure_Gradient;
  //@}
 
  //@{ @name Constructors and desctructors:
  //! Constructor (assign default values).
  CFD_Input_Parameters(void) {
    // Input file parameters:
    Get_CFFC_Path();
    strcpy(Input_File_Name, "CFFC.in");
    Line_Number = 0;
    strcpy(Next_Control_Parameter, " ");
    // Output parameters:
    strcpy(Output_File_Name, "outputfile.dat");
    strcpy(Restart_File_Name, "restartfile.soln");
    strcpy(Output_Format_Type, "Tecplot");
    i_Output_Format = IO_TECPLOT;
    Restart_Solution_Save_Frequency = 1000;
    Time_Accurate_Output_Frequency = 0;
    // Debugging parameters:
    Debug_Level = 0;  //default no debug information
    // Flow type indicator and related input parameters:
    strcpy(Flow_Type, "Inviscid");
    i_Flow_Type = FLOWTYPE_INVISCID;
    Axisymmetric = 0;
    // Time integration type indicator and related input parameters:
    strcpy(Time_Integration_Type, "Explicit_Euler");
    i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
    Time_Accurate = 0;
    Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
    Maximum_Number_of_Time_Steps = 100;
    N_Stage = 1;
    CFL_Number = HALF;
    Time_Max = ZERO;
    p_Norm_Indicator = 1;              //This is a stupid NAME!
    Number_of_Residual_Norms =1;
    // Reconstruction type indicator and related input parameters:
    strcpy(Reconstruction_Type, "Least_Squares");
    i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;
    // Limiter type indicator and related input parameters:
    strcpy(Limiter_Type, "Barth_Jespersen");
    i_Limiter = LIMITER_BARTH_JESPERSEN;
    Freeze_Limiter = 0;
    Freeze_Limiter_Residual_Level = 1e-4;
    // Inviscid flux function type and related input parameters:
    strcpy(Flux_Function_Type, "HLLE");
    i_Flux_Function = FLUX_FUNCTION_HLLE;
    // Preconditioner indicator and related input parameters:
    Preconditioning = 0; //default off
    Mach_Number_Reference = ONE;
    // Implicit residual smoothing control parameters:
    Residual_Smoothing = 0;
    Residual_Smoothing_Epsilon = ZERO;
    Residual_Smoothing_Gauss_Seidel_Iterations = 2;
    // Initial condition type indicator and related input parameters:
    strcpy(ICs_Type, "Uniform");
    i_ICs = IC_UNIFORM;
    strcpy(Gas_Type, "AIR");
    Mach_Number = ZERO;
    Reynolds_Number = 500000.00;
    Pressure = PRESSURE_STDATM;
    Temperature = TEMPERATURE_STDATM;
    Flow_Angle = ZERO;
    // Other initial and boundary condition related input parameters:
    Moving_Wall_Velocity.zero(); 
    Pressure_Gradient.zero(); 
  }

  //! Destructor
  ~CFD_Input_Parameters(void){ }
  //@}

  //@{ @name Other Member functions:
  //! Obtain location/path of CFFC source directory from environment variable
  void Get_CFFC_Path(void);
  //! Open input file
  void Open_Input_File(void);
  //! Close input file
  void Close_Input_File(void);
  //! Broadcast input parameters to all processors
  void Broadcast(void);
  //! Read next input line from input file
  int Get_Next_Input_Control_Parameter(const bool read_control_parameter);
  //! Parse next input line
  int Parse_Next_Input_Control_Parameter(void);
  //! Process/parse all input from input file
  int Process_Input_Control_Parameter_File(char *Input_File_Name_ptr,
                                           int &Command_Flag);
  //! Check validity of specified input parameters
  int Check_Inputs(void);
  //@}

  //@{ @name Input-output operators:
  friend ostream &operator << (ostream &out_file,
                               const CFD_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
                               CFD_Input_Parameters &IP);
  void Output(ostream &out_file) const;
  void Output_Problem_Type(ostream &out_file) const;
  void Output_Solver_Type(ostream &out_file) const;
  void Output_ICsBCs_Types(ostream &out_file) const;
  void Output_Solution_Type(ostream &out_file) const;
  void Output_IO_Types(ostream &out_file) const;
  void Output_GridAMR_Types(ostream &out_file) const;
  //@}

};

#endif // _CFD_INPUT_INCLUDED
