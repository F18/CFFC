/**********************************************************************
 * LevelSet2DInput.h: Header file defining 2D Level Set input         *
 *                    parameter class.                                *
 **********************************************************************/

#ifndef _LEVELSET2D_INPUT_INCLUDED
#define _LEVELSET2D_INPUT_INCLUDED

// Include 2D Level Set state header file.

#ifndef _LEVELSET2D_STATE_INCLUDED
#include "LevelSet2DState.h"
#endif // _LEVELSET2D_STATE_INCLUDED

// Include the 2D cell header file.

#ifndef _CELL2D_INCLUDED
#include "../Grid/Cell2D.h"
#endif // _CELL2D_INCLUDED

// Include the 2D cell header file.

#ifndef _GRID2D_QUAD_BLOCK_INCLUDED
#include "../Grid/Grid2DQuad.h"
#endif // _GRID2D_QUAD_BLOCK_INCLUDED

// Include embedded boundary input header file.

#ifndef _EMBEDDEDBOUNDARIES2DINPUT_INCLUDED
#include "../Interface2D/EmbeddedBoundaries2D_Input.h"
#endif // _EMBEDDEDBOUNDARIES2DINPUT_INCLUDED

// Define the structures and classes.

#define	INPUT_PARAMETER_LENGTH_LEVELSET2D  128

// Enviroment flag for CFFC root directory path
#define PATHVAR_LEVELSET2D "CFFC_Path"

/*!
 * Class: LevelSet2D_Input_Parameters
 *
 * @brief Definition and manipulation of levelset2D input variables.
 *
 */
class LevelSet2D_Input_Parameters{
  private:
public:
  //@{ @name Input file parameters.
  //! Input file name:
  //! CFFC root directory path:
  char CFFC_Path[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  char Input_File_Name[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  //! Input file stream:
  ifstream Input_File;
  //! Input file line number:
  int Line_Number;
  //@}

  //@{ @name Time integration type indicator and related input parameters:
  char Time_Integration_Type[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  int i_Time_Integration;
  int Time_Accurate, Local_Time_Stepping, 
      Maximum_Number_of_Time_Steps, N_Stage;
  double Hamilton_Jacobi_CFL_Number, Time_Max;
  // Residual variable:
  int i_Residual_Variable;
  //@}

  //@{ @name Initial distance function:
  char Initial_Distance_Type[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  int i_Initial_Distance_Type;
  int Perturb_Distance_Function;
  double Extension_Distance;
  //@}

  //@{ @name Redistance type indicator and related input parameters:
  char Redistance_Criteria[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  int i_Redistance_Criteria;
  double Eikonal_Threshold;
  int Redistance_Frequency;
  double Redistance_Tolerance;
  int Number_of_Initial_Redistance_Iterations, Number_of_Redistance_Iterations;
  double Eikonal_CFL_Number;
  char Eikonal_Scheme[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  int i_Eikonal_Scheme;
  char Eikonal_Sign_Function[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  int i_Eikonal_Sign_Function;
  //@}

  //@{ @name Scalar (front speed) extension problem input parameters:
  int Number_of_Scalar_Extension_Iterations;
  double Scalar_Extension_CFL_Number;
  //@}

  //@{ @name Reconstruction type indicator and related input parameters:
  char Reconstruction_Type[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  int i_Reconstruction;
  //@}

  //@{ @name Limiter type indicator and related input parameters:
  char Limiter_Type[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  int i_Limiter;
  //@}

  //@{ @name Limiter type indicator and related input parameters:
  char BC_Type[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  int i_BC_Type;
  //@}

  //@{ @name Embedded boundary input parameters:
  EmbeddedBoundaries2D_Input_Parameters Interface_IP;
  //@}

  //@{ @name Curvature driven flow parameters:
  double Curvature_Speed;
  char Curvature_Scheme[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  int i_Curvature_Scheme;
  //@}

  //@{ @name Bullk flowfield.
  char BulkFlowField_Type[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  int i_BulkFlowField_Type;
  Vector2D V;
  //@}

  //@{ @name Grid type indicator and related input parameters:
  char Grid_Type[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  char NACA_Aerofoil_Type[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  int i_Grid;
  int Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Ghost_Cells,
      Number_of_Blocks_Idir, Number_of_Blocks_Jdir; 
  double Box_Width, Box_Height, Plate_Length, 
         Pipe_Length, Pipe_Radius, 
         Blunt_Body_Radius, Blunt_Body_Mach_Number,
         Grain_Length, Grain_Radius, Grain_To_Throat_Length,
         Nozzle_Length, Nozzle_Radius_Exit, Nozzle_Radius_Throat,
         Cylinder_Radius, Ellipse_Length_X_Axis, 
         Ellipse_Length_Y_Axis, Chord_Length, Orifice_Radius;
  double X_Scale, X_Rotate;
  Vector2D X_Shift;
  //@}

  //@{ @name Mesh stretching factor.
  int i_Mesh_Stretching;
  int Mesh_Stretching_Type_Idir;
  int Mesh_Stretching_Type_Jdir;
  double Mesh_Stretching_Factor_Idir;
  double Mesh_Stretching_Factor_Jdir;
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
  //! Curvature of the level set function refinement criteria.
  int Refinement_Criteria_Curvature;
  //! Zero contour of the level set function refinement criteria.
  int Refinement_Criteria_Zero_Level_Set;
  //! Smooth quad block flag:
  int i_Smooth_Quad_Block;
  //@}

  //@{ @name Output parameters:
  //! Output file name:
  char Output_File_Name[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  //! Multi-block mesh definition input file names:
  char Grid_File_Name[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  char Grid_Definition_File_Name[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  //! Restart file name:
  char Restart_File_Name[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  //! Gnuplot file name:
  char Gnuplot_File_Name[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  //! Next_Control_Parameter:
  char Next_Control_Parameter[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  //! Output format type indicator:
  char Output_Format_Type[INPUT_PARAMETER_LENGTH_LEVELSET2D];
  int i_Output_Format;
  int Restart_Solution_Save_Frequency;
  //! Output progress frequency:
  int Output_Progress_Frequency;
  //@}

  //@{ @name Multi-block solution-adaption and parallel domain decomposition input parameters:
  int Number_of_Processors, Number_of_Blocks_Per_Processor;
  //@}

  //@{ @name Obtain the CFFC root directory path:
  void get_cffc_path();
  //@}

  //@{ @name Input-output operators:
  friend ostream &operator << (ostream &out_file, const LevelSet2D_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file, LevelSet2D_Input_Parameters &IP);
  //@}

};

/*****************************************************************
 *  LevelSet2D_Input_Parameters::get_cffc_path -- Get CFFC path. *
 *****************************************************************/
inline void LevelSet2D_Input_Parameters::get_cffc_path(){
  char *string_ptr;
  // Check to see if environment varible exists.
  if (getenv(PATHVAR_LEVELSET2D) == NULL) {
    //Set default path
     string_ptr = "CFFC";
     strcpy(CFFC_Path, string_ptr);
  } else {
     //Set path specified by environment variable
     strcpy(CFFC_Path, getenv(PATHVAR_LEVELSET2D));
  }
}

/**********************************************************************
 * LevelSet2D_Input_Parameters -- Input-output operators.             *
 **********************************************************************/
inline ostream &operator << (ostream &out_file,
			     const LevelSet2D_Input_Parameters &IP) {
  out_file << setprecision(6);
  out_file << "\n  -> CFFC Path: " 
           << IP.CFFC_Path;
  out_file << "\n\n Solving 2D Level Set equations on multi-block solution-adaptive quadrilateral mesh.";
  out_file << "\n  -> Input File Name: " << IP.Input_File_Name;
  out_file << "\n  -> Time Accurate (Unsteady) Solution";
  out_file << "\n  -> Time Integration Type: " << IP.Time_Integration_Type;
  out_file << "\n  -> Reconstruction type: "  << IP.Reconstruction_Type;
  //out_file << "\n  -> Limiter: " << IP.Limiter_Type;
  out_file << "\n  -> Boundary Condition: " << IP.BC_Type;
  // Display interface variables:
  if (IP.Interface_IP.Component_List.Ni) out_file << IP.Interface_IP;
  out_file << "\n  -> Grid: " << IP.Grid_Type;
  switch(IP.i_Grid) {
  case GRID_SQUARE :
    out_file << "\n  -> Size of Solution Domain (m): " << IP.Box_Width;
    break;
  case GRID_RECTANGULAR_BOX :
    out_file << "\n  -> Width of Solution Domain (m): " << IP.Box_Width;
    out_file << "\n  -> Height of Solution Domain (m): " << IP.Box_Height;
    break;
  case GRID_ROCKET_MOTOR :
    out_file << "\n  -> Length of Grain (m): " << IP.Grain_Length;
    out_file << "\n  -> Radius of Grain (m): " << IP.Grain_Radius;
    out_file << "\n  -> Distance from Grain to Nozzle Throat (m): " 
	     << IP.Grain_To_Throat_Length;
    out_file << "\n  -> Length of the Nozzle (m): " << IP.Nozzle_Length;
    out_file << "\n  -> Radius of the Nozzle at Throat (m): " 
	     << IP.Nozzle_Radius_Throat;
    out_file << "\n  -> Radius of the Nozzle at Exit(m): " 
	     << IP.Nozzle_Radius_Exit;
    break;
  default:
    out_file << "\n  -> Size of Solution Domain (m): " << IP.Box_Width;
    break;
  };
  out_file << "\n  -> Mesh shift, scale, and rotate: " 
	   << IP.X_Shift << " " << IP.X_Scale << " " << IP.X_Rotate;
  out_file << "\n  -> Number of Blocks i-direction: " << IP.Number_of_Blocks_Idir;
  out_file << "\n  -> Number of Blocks j-direction: " << IP.Number_of_Blocks_Jdir;
  out_file << "\n  -> Number of Cells i-direction: " << IP.Number_of_Cells_Idir;
  out_file << "\n  -> Number of Cells j-direction: " << IP.Number_of_Cells_Jdir;
  out_file << "\n  -> Number of Ghost Cells: " << IP.Number_of_Ghost_Cells;
  out_file << "\n  -> Number of Initial Mesh Refinements: " 
	   << IP.Number_of_Initial_Mesh_Refinements;
  out_file << "\n  -> Number of Uniform Mesh Refinements: " 
	   << IP.Number_of_Uniform_Mesh_Refinements;
  out_file << "\n  -> Refinement Criteria: ";
  if (IP.Refinement_Criteria_Curvature)
  out_file << "\n     -> Curvature of the level set function refinement criteria.";
  if (IP.Refinement_Criteria_Zero_Level_Set)
  out_file << "\n     -> Zero contour of the level set function refinement criteria.";
  // Hamilton-Jacobi-type equation parameters:
  out_file << "\n  -> Hamilton-Jacobi-type equation parameters: ";
  out_file << "\n     -> CFL Number: " << IP.Hamilton_Jacobi_CFL_Number;
  if (IP.Time_Accurate) out_file << "\n     -> Maximum Time (ms): " << IP.Time_Max*THOUSAND;
  if (!IP.Time_Accurate) out_file << "\n     -> Maximum Number of Time Steps (iterations): " 
				  << IP.Maximum_Number_of_Time_Steps;
  // Initial extension problem parameters:
  out_file << "\n  -> Initial extension problem parameters: ";
  out_file << "\n     -> Extension type: " << IP.Initial_Distance_Type;
  if (IP.i_Initial_Distance_Type == LEVELSET_INITIAL_EXTENSION_GEOMETRIC)
    out_file << "\n  -> Extension distance: " << IP.Extension_Distance;
  out_file << "\n  -> Perturb initial distance function: " << IP.Perturb_Distance_Function;  
  // Eikonal equation parameters:
  out_file << "\n  -> Eikonal equation parameters: ";
  out_file << "\n     -> CFL Number: " << IP.Eikonal_CFL_Number;
  out_file << "\n     -> Redistance Tolerance: " << IP.Redistance_Tolerance;
  out_file << "\n     -> Number of Initial Redistance Iterations: " << IP.Number_of_Initial_Redistance_Iterations;
  out_file << "\n     -> Number of Redistance Iterations: " << IP.Number_of_Redistance_Iterations;
  out_file << "\n     -> Redistance Criteria: " << IP.Redistance_Criteria;
  if (IP.i_Redistance_Criteria == EIKONAL_CRITERIA_THRESHOLD) {
    out_file << "\n     -> Redistance Threshold: " << IP.Eikonal_Threshold;
  } else {
    out_file << "\n     -> Redistance Frequency: " << IP.Redistance_Frequency;
  }
  out_file << "\n     -> Eikonal scheme: " << IP.Eikonal_Scheme;
  out_file << "\n     -> Eikonal sign function: " << IP.Eikonal_Sign_Function;
  // Scalar extension parameters:
  out_file << "\n  -> Scalar (front speed) extension parameters: ";
  out_file << "\n     -> CFL Number: " << IP.Scalar_Extension_CFL_Number;
  out_file << "\n     -> Number of Iterations: " << IP.Number_of_Scalar_Extension_Iterations;
  // Curvature driven flow parameter:
  out_file << "\n  -> Curvature driven flow parameters: ";
  out_file << "\n     -> Curvature_Speed: " << IP.Curvature_Speed;
  out_file << "\n     -> Curvature_Scheme: " << IP.Curvature_Scheme;
  // Bulk flow-field parameters:
  out_file << "\n  -> Bulk flow-field parameters: ";
  switch(IP.i_BulkFlowField_Type) {
  case INTERFACE_BULKFLOWFIELD_NONE :
    out_file << "\n     -> No bulk flow-field.";
    break;
  case INTERFACE_BULKFLOWFIELD_UNIFORM :
    out_file << "\n     -> Uniform bulk flow-field with V = (" << IP.V.x << "," << IP.V.y << ")";
    break;
  case INTERFACE_BULKFLOWFIELD_SWIRL :
    out_file << "\n     -> Swirling bulk flow-field with omega = " << IP.V.x;
  default:
    break;
  };
  //
  out_file << "\n  -> Number of Processors: " << IP.Number_of_Processors;
  out_file << "\n  -> Number of Blocks Per Processor: " 
	   << IP.Number_of_Blocks_Per_Processor;
  out_file << "\n  -> Output File Name: " << IP.Output_File_Name;
  out_file << "\n  -> Output Format: " << IP.Output_Format_Type;
  out_file << "\n  -> Restart Solution Save Frequency: "
	   << IP.Restart_Solution_Save_Frequency
	   << " steps (iterations)"; 
  out_file << "\n  -> Output Progress Frequency: "
	   << IP.Output_Progress_Frequency
	   << " steps (iterations)";
  return out_file;
}

inline istream &operator >> (istream &in_file, LevelSet2D_Input_Parameters &IP) {
  return in_file;
}

/**********************************************************************
 * LevelSet2D_Input_Parameters -- External subroutines.               *
 **********************************************************************/

extern void Open_Input_File(LevelSet2D_Input_Parameters &IP);

extern void Close_Input_File(LevelSet2D_Input_Parameters &IP);

extern void Set_Default_Input_Parameters(LevelSet2D_Input_Parameters &IP);

extern void Broadcast_Input_Parameters(LevelSet2D_Input_Parameters &IP);

#ifdef _MPI_VERSION
extern void Broadcast_Input_Parameters(LevelSet2D_Input_Parameters &IP,
                                       MPI::Intracomm &Communicator, 
                                       const int Source_CPU);
#endif

extern void Get_Next_Input_Control_Parameter(LevelSet2D_Input_Parameters &IP);

extern int Parse_Next_Input_Control_Parameter(LevelSet2D_Input_Parameters &IP);

extern int Process_Input_Control_Parameter_File(LevelSet2D_Input_Parameters &Input_Parameters,
                                                char *Input_File_Name_ptr,
                                                int &Command_Flag);

#endif // _LEVELSET2D_INPUT_INCLUDED
