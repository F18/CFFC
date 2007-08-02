/**********************************************************************
 * Electrostatic2DInput.h: Header file defining 2D Electrostatic      *
 *                         input parameter class.                     *
 **********************************************************************/

#ifndef _ELECTROSTATIC2D_INPUT_INCLUDED
#define _ELECTROSTATIC2D_INPUT_INCLUDED

// Include 2D Electrostatic state header file.

#ifndef _ELECTROSTATIC2D_STATE_INCLUDED
#include "Electrostatic2DState.h"
#endif // _ELECTROSTATIC2D_STATE_INCLUDED

// Include the math header file.

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

// Include the 2D cell header file.

#ifndef _CELL2D_INCLUDED
#include "../Grid/Cell2D.h"
#endif // _CELL2D_INCLUDED

// Include the linked list header file.

#ifndef _LINKEDLIST_INCLUDED
#include "../Math/LinkedList.h"
#endif // _LINKEDLIST_INCLUDED

// Include the 2D quadrilateral multiblock grid header file.

#ifndef _GRID2D_QUAD_BLOCK_INCLUDED
#include "../Grid/Grid2DQuad.h"
#endif // _GRID2D_QUAD_BLOCK_INCLUDED

// Include ICEMCFD input header file.

#ifndef _ICEMCFD_INCLUDED
#include "../ICEM/ICEMCFD.h"
#endif // _ICEMCFD_INCLUDED

// Define the structures and classes.

#define	INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D    128

/*!
 * Class:  Electrostatic2D_Input_Parameters
 *
 * @brief Definition and manipulation of electrostatic2D input variables.
 *
 */
class Electrostatic2D_Input_Parameters{
private:
public:
  //@{ @name Input file parameters.
  //! Input file name:
  char Input_File_Name[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
  //! Input file stream:
  ifstream Input_File;
  //! Input file line number:
  int Line_Number;
  //@}

  //@{ @name Time integration type indicator and related input parameters:
  char Time_Integration_Type[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
  int i_Time_Integration;
  int Time_Accurate, Local_Time_Stepping, 
      Maximum_Number_of_Time_Steps, N_Stage;
  double CFL_Number, Time_Max;
  //@}

  //@{ @name Implicit residual smoothing control parameters:
  int Residual_Smoothing;
  double Residual_Smoothing_Epsilon;
  int Residual_Smoothing_Gauss_Seidel_Iterations;
  //@}

  //@{ @name Reconstruction type indicator and related input parameters:
  char Reconstruction_Type[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
  int i_Reconstruction;
  //@}

  //@{ @name Flux function type and related input parameters:
  char Flux_Reconstruction_Type[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
  int i_Flux_Reconstruction;
  //@}

  //@{ @name Initial condition type indicator and related input parameters:
  char ICs_Type[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
  int i_ICs;
  Electrostatic2DState Uo, U1, U2;
  double Potential, Electric_Field_Strength, Electric_Field_Angle;
  //@}

  //@{ @name Flow geometry (planar or axisymmetric):
  char Flow_Geometry_Type[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
  int Axisymmetric;
  //@}

  //@{ @name Grid type indicator and related input parameters:
  char Grid_Type[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
  int i_Grid;
  int Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Ghost_Cells,
      Number_of_Blocks_Idir, Number_of_Blocks_Jdir;
  double Box_Width, Box_Height,
         Cylinder_Radius, Orifice_Radius,
         Wedge_Angle, Wedge_Length,
         Chamber_Length, Chamber_Radius, Chamber_To_Throat_Length,
         Nozzle_Length, Nozzle_Radius_Exit, Nozzle_Radius_Throat, Grain_Radius;
  int Smooth_Bump, Nozzle_Type;
  double X_Scale, X_Rotate;
  Vector2D X_Shift;
  char **ICEMCFD_FileNames;
  //@}

  //@{ @name Mesh stretching factor.
  int i_Mesh_Stretching;
  int Mesh_Stretching_Type_Idir;
  int Mesh_Stretching_Type_Jdir;
  double Mesh_Stretching_Factor_Idir;
  double Mesh_Stretching_Factor_Jdir;
  //@}

  //@{ @name Boundary conditions:
  char Boundary_Conditions_Specified[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
  int BCs_Specified;
  char BC_North_Type[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
  char BC_South_Type[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
  char BC_East_Type[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
  char BC_West_Type[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
  int BC_North, BC_South, BC_East, BC_West;
  //@}

  //@{ @name AMR input parameters:
  //! Unsteady AMR flag.
  int AMR;
  //! Unsteady AMR frequency.
  int AMR_Frequency;
  //! Number of initial mesh refinements.
  int Number_of_Initial_Mesh_Refinements;
  //! Number of uniform mesh refinements.
  int Number_of_Uniform_Mesh_Refinements;
  //! Number of boundary mesh refinements.
  int Number_of_Boundary_Mesh_Refinements;
  //! Maximum number of refinement levels.
  int Maximum_Refinement_Level;
  //! Minimum number of refinement levels.
  int Minimum_Refinement_Level;
  //! Refinement threshold.
  double Threshold_for_Refinement;
  //! Coarsening threshhold.
  double Threshold_for_Coarsening;
  //! Number of refinement criteria.
  int Number_of_Refinement_Criteria;
  //! Gradient of the potential field refinement criteria.
  int Refinement_Criteria_Gradient_Potential_Field;
  //! Divergence of the electric field refinement criteria.
  int Refinement_Criteria_Divergence_Electric_Field;
  //! Curl of the electric field refinement criteria.
  int Refinement_Criteria_Curl_Electric_Field;
  //! Smooth quad block flag:
  int i_Smooth_Quad_Block;
  //@}

  //@{ @name Output parameters:
  //! Output file name:
  char Output_File_Name[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
  //! Multi-block mesh definition input file names:
  char Grid_File_Name[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
  char Grid_Definition_File_Name[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
  //! Restart file name:
  char Restart_File_Name[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
  //! Gnuplot file name:
  char Gnuplot_File_Name[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
  //! Next_Control_Parameter:
  char Next_Control_Parameter[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
  //! Output format type indicator:
  char Output_Format_Type[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
  int i_Output_Format;
  int Restart_Solution_Save_Frequency;
  //! Output progress frequency:
  int Output_Progress_Frequency;
  //@}

  //@{ @name Multi-block solution-adaption and parallel domain decomposition input parameters:
  int Number_of_Processors, Number_of_Blocks_Per_Processor;
  //@}

  //@{ @name Input-output operators:
  friend ostream &operator << (ostream &out_file, const Electrostatic2D_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file, Electrostatic2D_Input_Parameters &IP);
  //@}

};

/**********************************************************************
 * Electrostatic2D_Input_Parameters -- Input-output operators.        *
 **********************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Electrostatic2D_Input_Parameters &IP) {
  out_file << setprecision(6);
  out_file << "\n\n Solving 2D Electrostatic equations (IBVP/BVP) on multi-block solution-adaptive quadrilateral mesh.";
  out_file << "\n  -> Input File Name: " << IP.Input_File_Name;
  if (IP.Time_Accurate) { 
    out_file << "\n  -> Time Accurate (Unsteady) Solution";
  } else {
    out_file << "\n  -> Time Invariant (Steady-State) Solution";
  }
  if (!IP.Axisymmetric) out_file << "\n  -> 2D Planar Flow";
  else out_file << "\n  -> 2D Axisymmetric Flow";
  out_file << "\n  -> Time Integration: " << IP.Time_Integration_Type;
  out_file << "\n  -> Number of Stages in Multi-Stage Scheme: " << IP.N_Stage;
  if (IP.Local_Time_Stepping == GLOBAL_TIME_STEPPING) {
    out_file << "\n  -> Global Time Stepping";
  } else if (IP.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
    out_file << "\n  -> Scalar Local Time Stepping";
  } else if (IP.Local_Time_Stepping == SEMI_IMPLICIT_LOCAL_TIME_STEPPING) {
    out_file << "\n  -> Semi-Implicit Local Time Stepping";
  }
  if (IP.Residual_Smoothing) {
    out_file << "\n  -> Residual Smoothing:";
    out_file << "\n     -> Epsilon: " << IP.Residual_Smoothing_Epsilon;
    out_file << "\n     -> Gauss_Seidel_Iterations: " 
	     << IP.Residual_Smoothing_Gauss_Seidel_Iterations;
  }
  out_file << "\n  -> Reconstruction: " << IP.Reconstruction_Type;
  out_file << "\n  -> Flux reconstruction: " << IP.Flux_Reconstruction_Type;
  out_file << "\n  -> Initial Conditions: " << IP.ICs_Type;
  out_file << "\n  -> Potential: " << IP.Potential;
  out_file << "\n  -> Electric Field Strength (V/m): " << IP.Electric_Field_Strength;
  out_file << "\n  -> Electric Field Angle: " << IP.Electric_Field_Angle;
  out_file << "\n  -> Grid: " << IP.Grid_Type;
  switch(IP.i_Grid) {
  case GRID_CARTESIAN_UNIFORM :
    out_file << "\n  -> Width of Solution Domain (m): " << IP.Box_Width;
    out_file << "\n  -> Height of Solution Domain (m): " << IP.Box_Height;
        break;
  case GRID_SQUARE :
    out_file << "\n  -> Size of Solution Domain (m): " << IP.Box_Width;
    break;
  case GRID_RECTANGULAR_BOX :
    out_file << "\n  -> Width of Solution Domain (m): " << IP.Box_Width;
    out_file << "\n  -> Height of Solution Domain (m): " << IP.Box_Height;
    break;
  case GRID_CIRCULAR_CYLINDER :
    out_file << "\n  -> Cylinder Radius (m): " << IP.Cylinder_Radius;
    break;
  case GRID_FREE_JET :
    out_file << "\n  -> Orifice Radius (m): " << IP.Orifice_Radius;
    break;
  case GRID_WEDGE :
    out_file << "\n  -> Wedge Angle (degrees): " << IP.Wedge_Angle;
    out_file << "\n  -> Wedge Length (m): " << IP.Wedge_Length;
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
  }
  if (IP.BCs_Specified) {
    out_file << "\n  -> Boundary conditions specified as: "
	     << "\n     -> BC_North = " << IP.BC_North_Type
	     << "\n     -> BC_South = " << IP.BC_South_Type
	     << "\n     -> BC_East = " << IP.BC_East_Type
	     << "\n     -> BC_West = " << IP.BC_West_Type;
  }
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
  out_file << "\n  -> Number of Initial Mesh Refinements: " 
	   << IP.Number_of_Initial_Mesh_Refinements;
  out_file << "\n  -> Number of Uniform Mesh Refinements: " 
	   << IP.Number_of_Uniform_Mesh_Refinements;
  out_file << "\n  -> Number of Boundary Mesh Refinements: " 
	   << IP.Number_of_Boundary_Mesh_Refinements;
  out_file << "\n  -> Refinement Criteria: ";
  if (IP.Refinement_Criteria_Gradient_Potential_Field)
  out_file << "\n     -> Gradient of the potential field refinement criteria";
  if (IP.Refinement_Criteria_Divergence_Electric_Field)
  out_file << "\n     -> Divergence of the electric field refinement criteria";
  if (IP.Refinement_Criteria_Curl_Electric_Field)
  out_file << "\n     -> Curl of the electric field refinement criteria";
  out_file << "\n  -> Smooth Quad Block: "
	   << IP.i_Smooth_Quad_Block;
  out_file << "\n  -> CFL Number: "
	   << IP.CFL_Number;
  out_file << "\n  -> Maximum Time (ms): " 
	   << IP.Time_Max*THOUSAND;
  out_file << "\n  -> Maximum Number of Time Steps (Iterations): " 
	   << IP.Maximum_Number_of_Time_Steps;
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
  return out_file;
}

inline istream &operator >> (istream &in_file, Electrostatic2D_Input_Parameters &IP) {
  return in_file;
}

/**********************************************************************
 * Electrostatic2D_Input_Parameters -- External subroutines.          *
 **********************************************************************/

extern void Open_Input_File(Electrostatic2D_Input_Parameters &IP);

extern void Close_Input_File(Electrostatic2D_Input_Parameters &IP);

extern void Set_Default_Input_Parameters(Electrostatic2D_Input_Parameters &IP);

extern void Broadcast_Input_Parameters(Electrostatic2D_Input_Parameters &IP);

#ifdef _MPI_VERSION
extern void Broadcast_Input_Parameters(Electrostatic2D_Input_Parameters &IP,
                                       MPI::Intracomm &Communicator,
                                       const int Source_CPU);
#endif

extern void Get_Next_Input_Control_Parameter(Electrostatic2D_Input_Parameters &IP);

extern int Parse_Next_Input_Control_Parameter(Electrostatic2D_Input_Parameters &IP);

extern int Process_Input_Control_Parameter_File(Electrostatic2D_Input_Parameters &IP,
                                                char *Input_File_Name_ptr,
                                                int &Command_Flag);

extern void Initialize_Reference_State(Electrostatic2D_Input_Parameters &IP);

extern void Reinitialize_Reference_State(Electrostatic2D_Input_Parameters &IP);

#endif // _ELECTROSTATIC2D_INPUT_INCLUDED
