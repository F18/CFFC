/* Rte2DInput.h:  Header file defining 
                    2D Rte Input Parameter Class. */

#ifndef _RTE2D_INPUT_INCLUDED
#define _RTE2D_INPUT_INCLUDED

/* Include 2D Rte state, 2D cell, 2D quadrilateral multiblock 
   grid, and NASA rotor header files. */

#ifndef _RTE2D_STATE_INCLUDED
#include "Rte2DState.h"
#endif // _RTE2D_STATE_INCLUDED

#ifndef _CELL2D_INCLUDED
#include "../Grid/Cell2D.h"
#endif // _CELL2D_INCLUDED

#ifndef _GRID2D_QUAD_BLOCK_INCLUDED
#include "../Grid/Grid2DQuad.h"
#endif // _GRID2D_QUAD_BLOCK_INCLUDED

/* Also include multigrid input header file. */

#ifndef _FASMULTIGRID2DINPUT_INCLUDED
#include "../FASMultigrid2D/FASMultigrid2DInput.h"
#endif // _FASMULTIGRID2DINPUT_INCLUDED

#ifndef _NKS2DINPUT_INCLUDED
#include "../NewtonKrylovSchwarz2D/NKSInput2D.h"
#endif // _NKS2DINPUT_INCLUDED

// Include ICEMCFD input header file.

#ifndef _ICEMCFD_INCLUDED
#include "../ICEM/ICEMCFD.h"
#endif // _ICEMCFD_INCLUDED


/* Define the structures and classes. */

#define	INPUT_PARAMETER_LENGTH_RTE2D    128

/********************************************************
 * Class:  Rte2D_Input_Parameters                     *
 ********************************************************/
class Rte2D_Input_Parameters{
  private:
  public:
  
  // Input file name:
  char Input_File_Name[INPUT_PARAMETER_LENGTH_RTE2D];

  // Input file stream:
  ifstream Input_File;

  // Input file line number:
  int Line_Number;

  // Time integration type indicator and related input parameters:
  char Time_Integration_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  int i_Time_Integration;
  int Time_Accurate, Local_Time_Stepping, 
    Maximum_Number_of_Time_Steps, N_Stage;
  double CFL_Number, Time_Max;

  // Implicit residual smoothing control parameters:
  int Residual_Smoothing;
  double Residual_Smoothing_Epsilon;
  int Residual_Smoothing_Gauss_Seidel_Iterations;

  /***************************************
   * Multigrid related input parameters: */
  Multigrid_Input_Parameters Multigrid_IP;
  /***************************************/  

  NKS_Input_Parameters  NKS_IP;

  // Reconstruction type indicator and related input parameters:
  char Reconstruction_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  int i_Reconstruction;

  // Limiter type indicator and related input parameters:
  char Limiter_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  int i_Limiter;
  int  Freeze_Limiter;
  double Freeze_Limiter_Residual_Level;

  // space marching related parameters
  char SpaceMarch_Scheme[INPUT_PARAMETER_LENGTH_RTE2D];
  int i_SpaceMarch_Scheme;


  // Initial condition type indicator and related input parameters:
  char ICs_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  int i_ICs;
  Rte2D_State Uo;
  double AbsorbsionCoef, Temperature, ScatteringCoef;
  int i_ScatteringFunc;
  char ScatteringFunc[INPUT_PARAMETER_LENGTH_RTE2D];
  int Number_of_Angles_Mdir, Number_of_Angles_Ldir;
  double NorthWallTemp, SouthWallTemp, EastWallTemp, WestWallTemp;
  double NorthWallEmiss, SouthWallEmiss, EastWallEmiss, WestWallEmiss;
  int i_RTE_Solver;
  char RTE_Solver[INPUT_PARAMETER_LENGTH_RTE2D];
  int i_DOM_Quadrature;
  char DOM_Quadrature[INPUT_PARAMETER_LENGTH_RTE2D];

  // Flow geometry (planar or axisymmetric):
  char Flow_Geometry_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  int Axisymmetric;

  // Grid type indicator and related input parameters:
  char Grid_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  char NACA_Aerofoil_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  int i_Grid;
  int Number_of_Cells_Idir, Number_of_Cells_Jdir,
      Number_of_Blocks_Idir, Number_of_Blocks_Jdir; 
  double Box_Width, Box_Height, Plate_Length, 
         Pipe_Length, Pipe_Radius, 
         Blunt_Body_Radius, Blunt_Body_Mach_Number,
         Grain_Length, Grain_Radius, Grain_To_Throat_Length,
         Nozzle_Length, Nozzle_Radius_Exit, Nozzle_Radius_Throat,
         Cylinder_Radius, Ellipse_Length_X_Axis, 
         Ellipse_Length_Y_Axis, Chord_Length, Orifice_Radius,
         Wedge_Angle, Wedge_Length;
  double X_Scale, X_Rotate;
  Vector2D X_Shift;

  char **ICEMCFD_FileNames;

  // AMR
  int AMR, AMR_Frequency;

  // Number of initial mesh refinements.
  int Number_of_Initial_Mesh_Refinements;

  // Number of uniform mesh refinements.
  int Number_of_Uniform_Mesh_Refinements;

  // Number of boundary mesh refinements.
  int Number_of_Boundary_Mesh_Refinements;

  // Maximum and minimum number of refinement levels.
  int Maximum_Refinement_Level, Minimum_Refinement_Level;

  // Thresholds for refinement and coarsening of mesh.
  double Threshold_for_Refinement, Threshold_for_Coarsening;

  // Morton Re-Ordering
  int Morton, Morton_Reordering_Frequency;

  // Output file name:
  char Output_File_Name[INPUT_PARAMETER_LENGTH_RTE2D];

  // Multi-block mesh definition input file names:
  char Grid_File_Name[INPUT_PARAMETER_LENGTH_RTE2D];
  char Grid_Definition_File_Name[INPUT_PARAMETER_LENGTH_RTE2D];

  // Restart file name:
  char Restart_File_Name[INPUT_PARAMETER_LENGTH_RTE2D];

  // Gnuplot file name:
  char Gnuplot_File_Name[INPUT_PARAMETER_LENGTH_RTE2D];

  // Next_Control_Parameter:
  char Next_Control_Parameter[INPUT_PARAMETER_LENGTH_RTE2D];

  // Output format type indicator:
  char Output_Format_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  int i_Output_Format;
  int Restart_Solution_Save_Frequency;

  // Multi-block solution-adaption and parallel domain
  // decomposition input parameters

  int Number_of_Processors, Number_of_Blocks_Per_Processor, Number_of_Ghost_Cells;

  //Destructor
  ~Rte2D_Input_Parameters(){
    for (int i = 0; i < 3; i++) {
      delete[] ICEMCFD_FileNames[i];
    } /* endfor */
    delete[]  ICEMCFD_FileNames; ICEMCFD_FileNames=NULL;
  }

  /* Input-output operators. */

  friend ostream &operator << (ostream &out_file,
		               const Rte2D_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
			       Rte2D_Input_Parameters &IP);

};

/*************************************************************
 * Rte2D_Input_Parameters -- Input-output operators.       *
 *************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Rte2D_Input_Parameters &IP) {
    out_file << setprecision(6);
    if (IP.i_Grid == GRID_CARTESIAN_UNIFORM) {
       out_file << "\n\n Solving 2D Rte equations (IBVP/BVP) on uniform Cartesian mesh.";
    } else {
       out_file << "\n\n Solving 2D Rte equations (IBVP/BVP) on multi-block solution-adaptive quadrilateral mesh.";
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
    out_file << "\n  -> RTE Solver: " 
             << IP.RTE_Solver;
    out_file << "\n  -> Time Integration: " 
             << IP.Time_Integration_Type;
    out_file << "\n  -> Number of Stages in Multi-Stage Scheme: " 
             << IP.N_Stage;
    /*****************************************************
     *                     Multigrid                     */
    if (IP.i_Time_Integration == TIME_STEPPING_MULTIGRID)
      {
	out_file << IP.Multigrid_IP;
      }
    /*****************************************************/
    else if (IP.i_Time_Integration != TIME_STEPPING_SPACE_MARCH) {
      if (IP.Local_Time_Stepping == GLOBAL_TIME_STEPPING) {
	out_file << "\n  -> Global Time Stepping";
      } else if (IP.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
	out_file << "\n  -> Scalar Local Time Stepping";
      } else if (IP.Local_Time_Stepping == MATRIX_LOCAL_TIME_STEPPING) {
	out_file << "\n  -> Matrix Local Time Stepping";
      } else if (IP.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) {
	out_file << "\n  -> Low-Mach-Number Local Preconditioning (Weiss-Smith)";
      } /* endif */
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

    /*****************************************************/
    // space marching specific
    } else {
      out_file << "\n  -> Differencing scheme: "
	       << IP.SpaceMarch_Scheme;
    }

    //--------------------- Rte2D Specific ----------------------//      
    out_file << "\n  -> Initial Conditions: " 
             << IP.ICs_Type;
    switch(IP.i_ICs) {
      case IC_CONSTANT :
        out_file << "\n  -> Absorbsion Cefficient (m^-1): " 
                 << IP.AbsorbsionCoef;
        out_file << "\n  -> Scattering Cefficient (m^-1): " 
                 << IP.ScatteringCoef;
        out_file << "\n  -> Scattering Function: " 
                 << IP.ScatteringFunc;
        out_file << "\n  -> Temperature (K): " 
                 << IP.Temperature;
        break;
    case IC_UNIFORM :
        out_file << "\n  -> Absorbsion Cefficient (m^-1): " 
                 << IP.AbsorbsionCoef;
        out_file << "\n  -> Scattering Cefficient (m^-1): " 
                 << IP.ScatteringCoef;
        out_file << "\n  -> Scattering Function: " 
                 << IP.ScatteringFunc;
        out_file << "\n  -> Temperature (K): " 
                 << IP.Temperature;
        break;
      default:
        break;
    } /* endswitch */
    out_file << "\n  -> Wall Conditions: ";
    out_file << "\n  -> Wall Temperature (K): " 
	     << IP.NorthWallTemp << "  "
	     << IP.SouthWallTemp << "  "
	     << IP.EastWallTemp << "  "
	     << IP.WestWallTemp;
    out_file << "\n  -> Wall Emissivity: " 
	     << IP.NorthWallEmiss << "  "
	     << IP.SouthWallEmiss << "  "
	     << IP.EastWallEmiss << "  "
	     << IP.WestWallEmiss;

    if (IP.i_RTE_Solver == RTE2D_SOLVER_DOM) {
      out_file << "\n  -> Quadrature: "
	       << IP.DOM_Quadrature;      
    } else {
      out_file << "\n  -> Number of CA polar-direction: "
	       << IP.Number_of_Angles_Mdir;
      out_file << "\n  -> Number of CA azim-direction: " 
	       << IP.Number_of_Angles_Ldir;
    }
    //------------------- End Rte2D Specific --------------------//

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
    if (IP.Number_of_Initial_Mesh_Refinements > 0)
    out_file << "\n  -> Number of Initial Mesh Refinements : " 
             << IP.Number_of_Initial_Mesh_Refinements;
    if (IP.Number_of_Uniform_Mesh_Refinements > 0)
    out_file << "\n  -> Number of Uniform Mesh Refinements : " 
	     << IP.Number_of_Uniform_Mesh_Refinements;
    if (IP.Number_of_Boundary_Mesh_Refinements > 0)
    out_file << "\n  -> Number of Boundary Mesh Refinements : " 
	     << IP.Number_of_Boundary_Mesh_Refinements;
    out_file << "\n  -> CFL Number: " 
             << IP.CFL_Number;
    out_file << "\n  -> Maximum Time (ms): " 
             << IP.Time_Max*THOUSAND;
    out_file << "\n  -> Maximum Number of Time Steps (Iterations): " 
             << IP.Maximum_Number_of_Time_Steps;
    out_file << "\n  -> Maximum Number of NKS Iterations: " 
             << IP.NKS_IP.Maximum_Number_of_NKS_Iterations;
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
    if (IP.AMR) out_file << "\n  -> AMR Frequency: "
                         << IP.AMR_Frequency
                         << " steps (iterations)";
    if (IP.Morton) out_file << "\n  -> Morton Re-Ordering Frequency: "
                            << IP.Morton_Reordering_Frequency
                            << " steps (iterations)"; 
    return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Rte2D_Input_Parameters &IP) {
    return (in_file);
}

/*************************************************************
 * Rte2D_Input_Parameters -- External subroutines.         *
 *************************************************************/

extern void Open_Input_File(Rte2D_Input_Parameters &IP);

extern void Close_Input_File(Rte2D_Input_Parameters &IP);

extern void Set_Default_Input_Parameters(Rte2D_Input_Parameters &IP);

extern void Broadcast_Input_Parameters(Rte2D_Input_Parameters &IP);

#ifdef _MPI_VERSION
extern void Broadcast_Input_Parameters(Rte2D_Input_Parameters &IP,
                                       MPI::Intracomm &Communicator, 
                                       const int Source_CPU);
#endif

extern void Get_Next_Input_Control_Parameter(Rte2D_Input_Parameters &IP);

extern int Parse_Next_Input_Control_Parameter(Rte2D_Input_Parameters &IP);

extern int Process_Input_Control_Parameter_File(Rte2D_Input_Parameters &Input_Parameters,
                                                char *Input_File_Name_ptr,
                                                int &Command_Flag);

#endif /* _RTE2D_INPUT_INCLUDED  */
