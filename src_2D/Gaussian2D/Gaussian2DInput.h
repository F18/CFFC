/* Gaussian2DInput.h:  Header file defining 
                    2D Gaussian Input Parameter Class. */

#ifndef _GAUSSIAN2D_INPUT_INCLUDED
#define _GAUSSIAN2D_INPUT_INCLUDED

/* Include 2D Gaussian state, 2D cell, and 2D quadrilateral 
   multiblock grid */

#ifndef _GAUSSIAN2D_STATE_INCLUDED
#include "Gaussian2DState.h"
#endif // _GAUSSIAN2D_STATE_INCLUDED

#ifndef _CELL2D_INCLUDED
#include "../Grid/Cell2D.h"
#endif // _CELL2D_INCLUDED

#ifndef _GRID2D_QUAD_BLOCK_INCLUDED
#include "../Grid/Grid2DQuad.h"
#endif // _GRID2D_QUAD_BLOCK_INCLUDED

#ifndef _NASA_ROTOR37_INCLUDED
#include "../Grid/NASARotor37.h"
#endif // _NASA_ROTOR37_INCLUDED

#ifndef _NASA_ROTOR67_INCLUDED
#include "../Grid/NASARotor67.h"
#endif // _NASA_ROTOR67_INCLUDED

/* Include multigrid input header file. */

#ifndef _FASMULTIGRID2DINPUT_INCLUDED
#include "../FASMultigrid2D/FASMultigrid2DInput.h"
#endif // _FASMULTIGRID2DINPUT_INCLUDED

// Include embedded boundary input header file.

#ifndef _EMBEDDEDBOUNDARIES2DINPUT_INCLUDED
#include "../Interface2D/EmbeddedBoundaries2D_Input.h"
#endif // _EMBEDDEDBOUNDARIES2DINPUT_INCLUDED

/* Include ICEMCFD input header file. */

#ifndef _ICEMCFD_INCLUDED
#include "../ICEM/ICEMCFD.h"
#endif // _ICEMCFD_INCLUDED

/* Define the structures and classes. */

#define	INPUT_PARAMETER_LENGTH_GAUSSIAN2D    128

// Enviroment flag for CFFC root directory path
#define PATHVAR_GAUSSIAN2D "CFFC_Path"

/********************************************************
 * Class:  Gaussian2D_Input_Parameters                  *
 ********************************************************/
class Gaussian2D_Input_Parameters{
  private:
  public:
  // CFFC root directory path:
  char CFFC_Path[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];

  // GMRES restart:
  //int GMRES_Restart;

  // GMRES tolerance:
  //double GMRES_Toler;

  // Overall tolerance:
  double Overall_Toler;

  // GMRES overlap:
  //int GMRES_Overlap;

  // GMRES P_Switch:
  //int GMRES_P_Switch;

  // ILU(k) - level of fill
  //int GMRES_ILUK_Level_of_Fill;

  // Finite_Time_Step:
  int Finite_Time_Step;
  double Finite_Time_Step_Initial_CFL;

  // Normalization:
  int Normalization;

  // Input file name:
  char Input_File_Name[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];

  // Input file stream:
  ifstream Input_File;

  // Input file line number:
  int Line_Number;

  // Time integration type indicator and related input parameters:
  char Time_Integration_Type[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];
  int i_Time_Integration;
  int Time_Accurate, Local_Time_Stepping, 
    Maximum_Number_of_Time_Steps, N_Stage;
  //Maximum_Number_of_NKS_Iterations,
  //Maximum_Number_of_GMRES_Iterations;
  double CFL_Number, Time_Max;

  // Residual variable:
  int i_Residual_Variable;

  // Implicit residual smoothing control parameters:
  int Residual_Smoothing;
  double Residual_Smoothing_Epsilon;
  int Residual_Smoothing_Gauss_Seidel_Iterations;

  /***************************************
   * Multigrid related input parameters: */
  Multigrid_Input_Parameters Multigrid_IP;
  /***************************************/  

  // Reconstruction type indicator and related input parameters:
  char Reconstruction_Type[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];
  int i_Reconstruction;

  // Limiter type indicator and related input parameters:
  char Limiter_Type[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];
  int i_Limiter;
  int  Freeze_Limiter;
  double Freeze_Limiter_Residual_Level;

  // Flux_Function_Type and related input parameters
  char Flux_Function_Type[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];
  int i_Flux_Function;

  char Heat_Reconstruction_Type[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];
  int i_Heat_Reconstruction;

  // Initial condition type indicator and related input parameters:
  char ICs_Type[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];
  char Gas_Type[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];
  int i_ICs;
  Gaussian2D_pState Wo, W1, W2;
  Gaussian2D_cState Uo;
  double Pressure, Temperature, Mach_Number, Flow_Angle;

  // Flow geometry (planar or axisymmetric):
  char Flow_Geometry_Type[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];
  int Axisymmetric;

  // Heat Transfer (on or off)
  int Heat_Transfer;

  //Prandtl number
  double pr;  

  // Grid type indicator and related input parameters:
  char Grid_Type[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];
  char NACA_Aerofoil_Type[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];
  int i_Grid;
  int Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Ghost_Cells,
      Number_of_Blocks_Idir, Number_of_Blocks_Jdir; 
  double Box_Width, Box_Height, Plate_Length, 
         Pipe_Length, Pipe_Radius, 
         Blunt_Body_Radius, Blunt_Body_Mach_Number,
         Grain_Length, Grain_Radius, Grain_To_Throat_Length,
         Nozzle_Length, Nozzle_Radius_Exit, Nozzle_Radius_Throat,
         Cylinder_Radius, Cylinder_Radius2, Ellipse_Length_X_Axis, 
         Ellipse_Length_Y_Axis, Chord_Length, Orifice_Radius,
         Inner_Streamline_Number, Outer_Streamline_Number, Isotach_Line,
         Wedge_Angle, Wedge_Length, Couette_Plate_Separation,
         Couette_Plate_Velocity, Pressure_Drop;
  int Smooth_Bump;
  double X_Scale, X_Rotate;
  Vector2D X_Shift;

  //Boundary information
  double alpha_m, alpha_t;
  double T_damping;
  double Ramp_by_Mach_Number;
  int Number_of_Time_Steps_to_Ramp;
  char Boundary_Conditions_Specified[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];
  int BCs_Specified;
  char BC_North_Type[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];
  char BC_South_Type[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];
  char BC_East_Type[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];
  char BC_West_Type[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];
  int BC_North, BC_South, BC_East, BC_West;
  double Temperature_North_BC, Temperature_South_BC, Temperature_East_BC, Temperature_West_BC;

  char NASA_Rotor37_Data_Directory[INPUT_PARAMETER_LENGTH_GAUSSIAN2D],
       NASA_Rotor67_Data_Directory[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];
  int Rotor_Flow_Type;
  double Rotor_Percent_Span;
  NASARotor37 NASA_Rotor37;
  NASARotor67 NASA_Rotor67;

  char **ICEMCFD_FileNames;

  // AMR
  int AMR, AMR_Frequency;

  // Number of initial mesh refinements.
  int Number_of_Initial_Mesh_Refinements;

  //Number of Uniform mesh refinements.
  int Number_of_Uniform_Mesh_Refinements;

  //Number of boundary mesh refinements.
  int Number_of_Boundary_Mesh_Refinements;

  // Number of interface mesh refinements.
  int Number_of_Interface_Mesh_Refinements;

  // Number of bounding-box mesh refinements:
  int Number_of_Bounding_Box_Mesh_Refinements;

  // Interface refinement condition.
  int Interface_Refinement_Condition;

  //Maximum and minimum number of refinement levels.
  int Maximum_Refinement_Level, Minimum_Refinement_Level;

  //Thresholds for refinement and coarsening of mesh.
  double Threshold_for_Refinement, Threshold_for_Coarsening;

  // Number of refinement criteria.
  int Number_of_Refinement_Criteria;
  // Gradient of the density field refinement criteria.
  int Refinement_Criteria_Gradient_Density;
  // Divergence of the velocity field refinement criteria.
  int Refinement_Criteria_Divergence_Velocity;
  // Curl of the velocity field refinement criteria.
  int Refinement_Criteria_Curl_Velocity;

  // Bounding-box for bounding-box mesh refinement.
  Vector2D AMR_Xmin, AMR_Xmax;

  // Mesh stretching factor.
  int i_Mesh_Stretching;
  int Mesh_Stretching_Type_Idir;
  int Mesh_Stretching_Type_Jdir;
  double Mesh_Stretching_Factor_Idir;
  double Mesh_Stretching_Factor_Jdir;


  // Smooth quad block flag:
  int i_Smooth_Quad_Block;

  // Embedded boundary input parameters:
  EmbeddedBoundaries2D_Input_Parameters Interface_IP;
  int Reset_Interface_Motion_Type;

  // Output file name:
  char Output_File_Name[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];

  // Multi-block mesh definition input file names:
  char Grid_File_Name[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];
  char Grid_Definition_File_Name[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];

  // Restart file name:
  char Restart_File_Name[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];

  // Gnuplot file name:
  char Gnuplot_File_Name[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];

  // Next_Control_Parameter:
  char Next_Control_Parameter[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];

  // Output format type indicator:
  char Output_Format_Type[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];
  int i_Output_Format;
  int Restart_Solution_Save_Frequency;

  // Output progress frequency:
  int Output_Progress_Frequency;

  // Multi-block solution-adaption and parallel domain
  // decomposition input parameters

  int Number_of_Processors, Number_of_Blocks_Per_Processor;

  // Obtain the CFFC root directory path:
  void get_cffc_path();

  /* Input-output operators. */

  friend ostream &operator << (ostream &out_file,
		               const Gaussian2D_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
			       Gaussian2D_Input_Parameters &IP);

};

/****************************************************************
 * Gaussian2D_Input_Parameters::get_cffc_path -- Get CFFC path. *
 ****************************************************************/
inline void Gaussian2D_Input_Parameters::get_cffc_path(){
  char *string_ptr;
  // Check to see if environment varible exists.
  if (getenv(PATHVAR_GAUSSIAN2D) == NULL) {
    //Set default path
     string_ptr = "CFFC";
     strcpy(CFFC_Path, string_ptr);
  } else {
     //Set path specified by environment variable
     strcpy(CFFC_Path, getenv(PATHVAR_GAUSSIAN2D));
  }
}

/****************************************************************
 * Gaussian2D_Input_Parameters -- Input-output operators.       *
 ****************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Gaussian2D_Input_Parameters &IP) {
    out_file << setprecision(6);
    out_file << "\n  -> CFFC Path: " 
	     << IP.CFFC_Path;
    if (IP.i_Grid == GRID_CARTESIAN_UNIFORM) {
       out_file << "\n\n Solving 2D Gaussian equations (IBVP/BVP) on uniform Cartesian mesh.";
    } else {
       out_file << "\n\n Solving 2D Gaussian equations (IBVP/BVP) on multi-block solution-adaptive quadrilateral mesh.";
    } /* endif */
    out_file << "\n  -> Input File Name: " 
             << IP.Input_File_Name;
    if (IP.Time_Accurate) { 
       out_file << "\n  -> Time Accurate (Unsteady) Solution";
    } else {
       out_file << "\n  -> Time Invariant (Steady-State) Solution";
    }
    if (IP.Axisymmetric) { 
       out_file << "\n  -> 2D Axisymmetric Flow";
    } else {
       out_file << "\n  -> 2D Planar Flow";
    }
    if(IP.Heat_Transfer) {
       out_file << "\n  -> Heat Transfer added";
       out_file << "\n  -> Heat flux evaluation: " << IP.Heat_Reconstruction_Type;
       out_file << "\n  -> Prandtl number: " << IP.pr;
       out_file << "\n  -> Temperature_North_BC: " << IP.Temperature_North_BC;
       out_file << "\n  -> Temperature_South_BC: " << IP.Temperature_South_BC;
       out_file << "\n  -> Temperature_East_BC: " << IP.Temperature_East_BC;
       out_file << "\n  -> Temperature_West_BC: " << IP.Temperature_West_BC;
    }
    out_file << "\n  -> Time Integration: " 
             << IP.Time_Integration_Type;
    out_file << "\n  -> Number of Stages in Multi-Stage Scheme: " 
             << IP.N_Stage;
    /*****************************************************
     *                     Multigrid                     */
    if (IP.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	IP.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
      out_file << IP.Multigrid_IP;
    }
    /*****************************************************/
    if (IP.Local_Time_Stepping == GLOBAL_TIME_STEPPING) {
      out_file << "\n  -> Global Time Stepping";
    } else if (IP.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
      out_file << "\n  -> Scalar Local Time Stepping";
    } /* endif */
    out_file << "\n  -> L1-, L2-, and max-norms computed on residual variable: " << IP.i_Residual_Variable;
    if (IP.Residual_Smoothing) {
      out_file << "\n  -> Residual Smoothing";
      out_file << "\n  -> Epsilon: " 
               << IP.Residual_Smoothing_Epsilon;
      out_file << "\n  -> Gauss_Seidel_Iterations: " 
               << IP.Residual_Smoothing_Gauss_Seidel_Iterations;
    } /* endif */
    out_file << "\n  -> Reconstruction: " 
             << IP.Reconstruction_Type;
    out_file << "\n  -> Limiter: " 
             << IP.Limiter_Type;
    if (IP.Limiter_Type != LIMITER_ZERO && IP.Freeze_Limiter) {
      out_file << "\n  -> Freeze Limiter when L2-norm of residual is < "
	       << IP.Freeze_Limiter_Residual_Level;
    } /* endif */
    out_file << "\n  -> Flux Function: " 
             << IP.Flux_Function_Type;
    out_file << "\n  -> Initial Conditions: " 
             << IP.ICs_Type;
    out_file << "\n  -> Gas Type: " 
             << IP.Gas_Type;
    switch(IP.i_ICs) {
      case IC_CONSTANT :
        out_file << "\n  -> Pressure (kPa): " 
                 << IP.Pressure/THOUSAND;
        out_file << "\n  -> Temperature (K): " 
                 << IP.Temperature;
        out_file << "\n  -> Mach Number: " 
                 << IP.Mach_Number;
        out_file << "\n  -> Flow Angle: " 
                 << IP.Flow_Angle;
        break;
      case IC_UNIFORM :
        out_file << "\n  -> Pressure (kPa): " 
                 << IP.Pressure/THOUSAND;
        out_file << "\n  -> Temperature (K): " 
                 << IP.Temperature;
        out_file << "\n  -> Mach Number: " 
                 << IP.Mach_Number;
        out_file << "\n  -> Flow Angle (degrees): " 
                 << IP.Flow_Angle;
        break;
      case IC_SOD_XDIR :
        break;
      case IC_SOD_YDIR :
        break;
      case IC_GROTH_XDIR :
        break;
      case IC_GROTH_YDIR :
        break;
      case IC_EINFELDT_XDIR :
        break;
      case IC_EINFELDT_YDIR :
        break;
      case IC_SHOCK_BOX :
        break;
      case IC_HIGH_PRESSURE_RESERVOIR :
        break;
      case IC_LOW_PRESSURE_RESERVOIR :
        break;
      case IC_RIEMANN :
      case IC_RIEMANN_XDIR :
      case IC_RIEMANN_YDIR :
        out_file << "\n  -> Pressure (kPa): " 
                 << IP.Pressure/THOUSAND;
        out_file << "\n  -> Temperature (K): " 
                 << IP.Temperature;
        out_file << "\n  -> Mach Number: " 
                 << IP.Mach_Number;
        out_file << "\n  -> Flow Angle (degrees): " 
                 << IP.Flow_Angle;
        break;
      default:
        break;
    } /* endswitch */
    out_file << "\n  -> Grid: " 
             << IP.Grid_Type;
    switch(IP.i_Grid) {
      case GRID_CARTESIAN_UNIFORM :
        out_file << "\n  -> Width of Solution Domain (m): " 
                 << IP.Box_Width;
        out_file << "\n  -> Height of Solution Domain (m): " 
                 << IP.Box_Height;
        break;
      case GRID_SQUARE :
        out_file << "\n  -> Size of Solution Domain (m): " 
                 << IP.Box_Width;
        break;
      case GRID_RECTANGULAR_BOX :
        out_file << "\n  -> Width of Solution Domain (m): " 
                 << IP.Box_Width;
        out_file << "\n  -> Height of Solution Domain (m): " 
                 << IP.Box_Height;
        break;
      case GRID_FLAT_PLATE :
        out_file << "\n  -> Plate Length (m): " 
                 << IP.Plate_Length;
        break;
      case GRID_PIPE :
        out_file << "\n  -> Pipe Length (m): " 
                 << IP.Pipe_Length;
        out_file << "\n  -> Pipe Radius (m): " 
                 << IP.Pipe_Radius;
        break;
      case GRID_BLUNT_BODY :
        out_file << "\n  -> Cylinder Radius (m): " 
                 << IP.Blunt_Body_Radius;
        break;
      case GRID_ROCKET_MOTOR :
        out_file << "\n  -> Length of Grain (m): " 
                 << IP.Grain_Length;
        out_file << "\n  -> Radius of Grain (m): " 
                 << IP.Grain_Radius;
        out_file << "\n  -> Distance from Grain to Nozzle Throat (m): " 
                 << IP.Grain_To_Throat_Length;
        out_file << "\n  -> Length of the Nozzle (m): " 
                 << IP.Nozzle_Length;
        out_file << "\n  -> Radius of the Nozzle at Throat (m): " 
                 << IP.Nozzle_Radius_Throat;
        out_file << "\n  -> Radius of the Nozzle at Exit (m): " 
                 << IP.Nozzle_Radius_Exit;
        break;
      case GRID_CIRCULAR_CYLINDER :
        out_file << "\n  -> Cylinder Radius (m): " 
                 << IP.Cylinder_Radius;
        out_file << "\n  -> Cylinder Outer Radius (m): " 
                 << IP.Cylinder_Radius2;
        break;
      case GRID_ELLIPSE :
        out_file << "\n  -> Width of Ellipse along x-axis (m): " 
                 << IP.Ellipse_Length_X_Axis;
        out_file << "\n  -> Height of Ellipse along y-axis (m): " 
                 << IP.Ellipse_Length_Y_Axis;
        break;
      case GRID_NACA_AEROFOIL :
        out_file << "\n  -> NACA " 
                 << IP.NACA_Aerofoil_Type;
        out_file << "\n  -> Chord Length (m): " 
                 << IP.Chord_Length;
        break;
      case GRID_FREE_JET :
        out_file << "\n  -> Orifice Radius (m): " 
                 << IP.Orifice_Radius;
        break;
      case GRID_WEDGE :
        out_file << "\n  -> Wedge Angle (degrees): " << IP.Wedge_Angle;
	out_file << "\n  -> Wedge Length (m): " << IP.Wedge_Length;
	break;
      case GRID_UNSTEADY_BLUNT_BODY :
	out_file << "\n  -> Cylinder Radius (m): " 
		 << IP.Blunt_Body_Radius;
        break;
        //case GRID_NASA_ROTOR_37 :
        //out_file << "\n  -> Percent Span: " 
        //         << IP.Rotor_Percent_Span;
        //break;
        //case GRID_NASA_ROTOR_67 :
        //out_file << "\n  -> Percent Span: " 
        //         << IP.Rotor_Percent_Span;
        //break;
      case GRID_BUMP_CHANNEL_FLOW :
        if (strcmp(IP.Grid_Type,"Bump_Channel_Flow") == 0) {
	  if (IP.Smooth_Bump) out_file << " (smooth bump)";
	  else out_file << " (non-smooth bump)";
	}
	break;

      case GRID_ADIABATIC_FLAT_PLATE :
        out_file << "\n  -> Plate Length (m): " 
                 << IP.Plate_Length;
        break;
      case GRID_ADIABATIC_PIPE :
        out_file << "\n  -> Pipe Length (m): " 
                 << IP.Pipe_Length;
        out_file << "\n  -> Pipe Radius (m): " 
                 << IP.Pipe_Radius;
        out_file << "\n  -> Pipe Pressure Drop (Kpa): " 
                 << IP.Pressure_Drop;
        break;
      case GRID_ADIABATIC_CIRCULAR_CYLINDER :
        out_file << "\n  -> Cylinder Radius (m): " 
                 << IP.Cylinder_Radius;
        break;
      case GRID_ADIABATIC_COUETTE :
        out_file << "\n  -> Plate Separation (m): " 
                 << IP.Couette_Plate_Separation;
        out_file << "\n  -> Plate Velocity (m/s): " 
                 << IP.Couette_Plate_Velocity;
        out_file << "\n  -> Pipe Pressure Drop (Kpa): " 
                 << IP.Pressure_Drop;
        break;
      case GRID_ICEMCFD :
        break;
      case GRID_READ_FROM_DEFINITION_FILE :
        break;
      case GRID_READ_FROM_GRID_DATA_FILE :
        break;
      default:
        out_file << "\n  -> Width of Solution Domain (m): " 
                 << IP.Box_Width;
        out_file << "\n  -> Height of Solution Domain (m): " 
                 << IP.Box_Height;
        break;
    } /* endswitch */
    out_file << "\n  -> Mesh shift, scale, and rotate: " 
             << IP.X_Shift << " " << IP.X_Scale << " " << IP.X_Rotate;
    out_file << "\n  -> Number of Blocks i-direction: "
             << IP.Number_of_Blocks_Idir;
    out_file << "\n  -> Number of Blocks j-direction: " 
             << IP.Number_of_Blocks_Jdir;
    out_file << "\n  -> Number of Cells i-direction: "
             << IP.Number_of_Cells_Idir;
    out_file << "\n  -> Number of Cells j-direction: " 
             << IP.Number_of_Cells_Jdir;
    out_file << "\n  -> Number of Ghost Cells: "
	     << IP.Number_of_Ghost_Cells;
    if (IP.Interface_IP.Component_List.Ni) out_file << IP.Interface_IP;
    if(IP.i_ICs != IC_RESTART)
    out_file << "\n  -> Momentum Accomodation Coefficient: "
	     << IP.alpha_m;
    out_file << "\n  -> Thermal Accomodation Coefficient: "
	     << IP.alpha_t;
    out_file << "\n  -> Damping for Slip_T Boundary conditions: "
	     << IP.T_damping;
    if(IP.Number_of_Initial_Mesh_Refinements >0)
    out_file << "\n  -> Number of Initial Mesh Refinements : " 
             << IP.Number_of_Initial_Mesh_Refinements;
    if (IP.Number_of_Uniform_Mesh_Refinements > 0)
    out_file << "\n  -> Number of Uniform Mesh Refinements : " 
	     << IP.Number_of_Uniform_Mesh_Refinements;
    if (IP.Number_of_Boundary_Mesh_Refinements > 0)
    out_file << "\n  -> Number of Boundary Mesh Refinements : " 
	     << IP.Number_of_Boundary_Mesh_Refinements;
    out_file << "\n  -> Number of Interface Mesh Refinements: " 
	     << IP.Number_of_Interface_Mesh_Refinements;
    out_file << "\n  -> Number of Bounding-Box Mesh Refinements: " 
	     << IP.Number_of_Bounding_Box_Mesh_Refinements;
    if (IP.Number_of_Bounding_Box_Mesh_Refinements)
      out_file << "\n     -> Bounding-box for bounding-box AMR:" << IP.AMR_Xmin << IP.AMR_Xmax;
    out_file << "\n  -> CFL Number: " 
             << IP.CFL_Number;
    out_file << "\n  -> Maximum Time (ms): " 
             << IP.Time_Max*THOUSAND;
    out_file << "\n  -> Maximum Number of Time Steps (Iterations): " 
             << IP.Maximum_Number_of_Time_Steps;
    //out_file << "\n  -> Maximum Number of NKS Iterations: " 
    //         << IP.Maximum_Number_of_NKS_Iterations;
    out_file << "\n  -> Number of Processors: " 
             << IP.Number_of_Processors;
    out_file << "\n  -> Number of Blocks Per Processor: " 
             << IP.Number_of_Blocks_Per_Processor;
    out_file << "\n  -> Output File Name: " 
             << IP.Output_File_Name;
    out_file << "\n  -> Output Format: " 
             << IP.Output_Format_Type;
    out_file << "\n  -> Restart Solution Save Frequency: "
             << IP.Restart_Solution_Save_Frequency
             << " steps (iterations)"; 
    out_file << "\n  -> Output Progress Frequency: "
	     << IP.Output_Progress_Frequency
	     << " steps (iterations)";
    return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Gaussian2D_Input_Parameters &IP) {
    return (in_file);
}

/*************************************************************
 * Gaussian2D_Input_Parameters -- External subroutines.         *
 *************************************************************/

extern void Open_Input_File(Gaussian2D_Input_Parameters &IP);

extern void Close_Input_File(Gaussian2D_Input_Parameters &IP);

extern void Set_Default_Input_Parameters(Gaussian2D_Input_Parameters &IP);

extern void Broadcast_Input_Parameters(Gaussian2D_Input_Parameters &IP);

#ifdef _MPI_VERSION
extern void Broadcast_Input_Parameters(Gaussian2D_Input_Parameters &IP,
                                       MPI::Intracomm &Communicator, 
                                       const int Source_CPU);
#endif

extern void Get_Next_Input_Control_Parameter(Gaussian2D_Input_Parameters &IP);

extern int Parse_Next_Input_Control_Parameter(Gaussian2D_Input_Parameters &IP);

extern int Process_Input_Control_Parameter_File(Gaussian2D_Input_Parameters &Input_Parameters,
                                                char *Input_File_Name_ptr,
                                                int &Command_Flag);

#endif /* _GAUSSIAN2D_INPUT_INCLUDED  */
