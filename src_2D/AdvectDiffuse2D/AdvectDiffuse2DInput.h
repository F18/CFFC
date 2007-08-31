/* AdvectDiffuse2DInput.h:  Header file defining 
                            2D Advection Diffusion Equation Input Parameter Class. */

#ifndef _ADVECTDIFFUSE2D_INPUT_INCLUDED
#define _ADVECTDIFFUSE2D_INPUT_INCLUDED

/* Include 2D advection diffusion equation solution state, 
   2D cell, and 2D quadrilateral multiblock grid header files. */

#ifndef _ADVECTDIFFUSE2D_STATE_INCLUDED
#include "AdvectDiffuse2DState.h"
#endif // _ADVECTDIFFUSE2D_STATE_INCLUDED

#ifndef _CELL2D_INCLUDED
#include "../Grid/Cell2D.h"
#endif // _CELL2D_INCLUDED

#ifndef _GRID2D_QUAD_BLOCK_INCLUDED
#include "../Grid/Grid2DQuad.h"
#endif // _GRID2D_QUAD_BLOCK_INCLUDED

/* Include multigrid input header file. */

#ifndef _FASMULTIGRID2DINPUT_INCLUDED
#include "../FASMultigrid2D/FASMultigrid2DInput.h"
#endif // _FASMULTIGRID2DINPUT_INCLUDED

/* Include ICEMCFD input header file. */

#ifndef _ICEMCFD_INCLUDED
#include "../ICEM/ICEMCFD.h"
#endif // _ICEMCFD_INCLUDED

/* Define the structures and classes. */

#define	INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D    128

// Enviroment flag for CFFC root directory path
#define PATHVAR_ADVECTDIFFUSE2D "CFFC_Path"

#define VELOCITY_FIELD_ZERO        0
#define VELOCITY_FIELD_UNIFORM     1
#define VELOCITY_FIELD_ROTATING    2

/*!
 * Class: AdvectDiffuse2D_Input_Parameters
 *
 * @brief Definition and manipulation of 2D Advection-Diffusion input variables.
 *
 */
class AdvectDiffuse2D_Input_Parameters{
private:
public:
  //@{ @name Input file parameters.
  //! CFFC root directory path:
  char CFFC_Path[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  //! Input file name:
  char Input_File_Name[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  //! Input file stream:
  ifstream Input_File;
  //! Input file line number:
  int Line_Number;
  //@}

  //@{ @name Time integration type indicator and related input parameters:
  char Time_Integration_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  int i_Time_Integration;
  int Time_Accurate, Local_Time_Stepping, 
      Maximum_Number_of_Time_Steps, N_Stage;
  double CFL_Number, Time_Max;
  //! Residual variable:
  int i_Residual_Variable;
  //@}

  //@{ @name Newton-Krylov-Schwarz input parameters.
  //! Maximum number of Newton-Krylov-Schwarz iterations.
  int Maximum_Number_of_NKS_Iterations;
  //! Maximum number of GMRES iterations.
  int Maximum_Number_of_GMRES_Iterations;
  //! GMRES restart:
  int GMRES_Restart;
  //! GMRES tolerance:
  double GMRES_Toler;
  //! Overall tolerance:
  double Overall_Toler;
  //! GMRES overlap:
  int GMRES_Overlap;
  //! GMRES P_Switch:
  int GMRES_P_Switch;
  //! ILU(k) - level of fill
  int GMRES_ILUK_Level_of_Fill;
  //! Finite_Time_Step:
  int Finite_Time_Step;
  double Finite_Time_Step_Initial_CFL;
  //! Normalization:
  int Normalization;
  //@}

  //@{ @name Implicit residual smoothing control parameters:
  int Residual_Smoothing;
  double Residual_Smoothing_Epsilon;
  int Residual_Smoothing_Gauss_Seidel_Iterations;
  //@}

  //@{ @name Multigrid related input parameters:
  Multigrid_Input_Parameters Multigrid_IP;
  //@}

  //@{ @name Reconstruction type indicator and related input parameters:
  char Reconstruction_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D]; /*!< Spatial reconstruction type. */
  int i_Reconstruction;		/*!< Index to store the reconstruction type. */
  int i_ReconstructionMethod;	/*!< Index to store the reconstruction method. */
  int Space_Accuracy;		/*!< Parameter to show the order of accuracy in space. */
  double CENO_Cutoff;		/*!< Value used by the CENO reconstruction to differentiate 
				 between smooth and non-smooth solution content.*/
  double CENO_RefinementUnits;	/*!< Parameter used in the AMR process when the CENO approach is employed. */
  int IncludeHighOrderBoundariesRepresentation;	/*!< Flag for including or excluding high-order BCs. */
  //@}

  //@{ @name Limiter type indicator and related input parameters:
  char Limiter_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  int i_Limiter;
  int  Freeze_Limiter;
  double Freeze_Limiter_Residual_Level;
  //@}

  //@{ @name Flux function type:
  int i_Flux_Function;
  //@}

  //@{ @name Initial condition type indicator and related input parameters:
  char ICs_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  int i_ICs;
  AdvectDiffuse2D_State Uo, U1, U2;
  FunctionType2D ExactSolution;
  int ExactSolution_Flag; 	        /*!< Flag used by MPI and RESTART to set the ExactSolution pointer. */
  AdvectDiffuse2D_State RefU;		/*!< Reference state, used by CENO to normalize the
					   variables in the computation of the smoothness indicator. */
  //@}

  //@{ @name Diffusion coefficient, advection speeds, and relaxation time:
  double Kappa, a, b, Tau;
  FunctionType2D KappaVar;          /*!< Function pointer which is set to the diffusion coefficient variation. */
  FunctionType2D SourceTermVar;	    /*!< Function pointer which is set to the source term variation. */
  //@}

  //@{ @name Convection velocity field type parameters:
  char Velocity_Field_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  int i_Velocity_Field;
  FunctionType2D VelocityField;	    /*!< Function pointer which is set to the velocity field. */
  //@}

  //@{ @name Flow geometry (planar or axisymmetric):
  char Flow_Geometry_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  int Axisymmetric;
  //@}

  //@{ @name Grid type indicator and related input parameters:
  char Grid_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  char NACA_Aerofoil_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  int i_Grid;
  int Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Ghost_Cells,
      Number_of_Blocks_Idir, Number_of_Blocks_Jdir;
  double Box_Width, Box_Height, Plate_Length, 
         Pipe_Length, Pipe_Radius, 
         Blunt_Body_Radius, Blunt_Body_Mach_Number,
         Chamber_Length, Chamber_Radius, Chamber_To_Throat_Length,
         Nozzle_Length, Nozzle_Radius_Exit, Nozzle_Radius_Throat, Grain_Radius,
         Cylinder_Radius, Ellipse_Length_X_Axis, 
         Ellipse_Length_Y_Axis, Chord_Length, Orifice_Radius;
  int Nozzle_Type;
  double X_Scale, X_Rotate;
  Vector2D X_Shift;
  char **ICEMCFD_FileNames;
  int IterationsOfInteriorNodesDisturbancies; /* <! Number of iterations of disturbing the mesh
						 (create an unsmooth interior mesh). */
  //@}

  //@{ @name Mesh stretching factor.
  int i_Mesh_Stretching;
  int Mesh_Stretching_Type_Idir;
  int Mesh_Stretching_Type_Jdir;
  double Mesh_Stretching_Factor_Idir;
  double Mesh_Stretching_Factor_Jdir;
  //@}

  //@{ @name Boundary conditions:
  char Boundary_Conditions_Specified[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  int BCs_Specified;
  char BC_North_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  char BC_South_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  char BC_East_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  char BC_West_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
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
  //! Threshold for refinement of mesh.
  double Threshold_for_Refinement;
  //! Threshold for coarsening of mesh.
  double Threshold_for_Coarsening;
  //! Smooth quad block flag:
  int i_Smooth_Quad_Block;
  //@}

  //@{ @name Output parameters:
  //! Output file name:
  char Output_File_Name[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  //! Multi-block mesh definition input file names:
  char Grid_File_Name[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  char Grid_Definition_File_Name[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  //! Restart file name:
  char Restart_File_Name[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  //! Gnuplot file name:
  char Gnuplot_File_Name[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  //! Next_Control_Parameter:
  char Next_Control_Parameter[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
  //! Output format type indicator:
  char Output_Format_Type[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
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
  friend ostream &operator << (ostream &out_file,
		               const AdvectDiffuse2D_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
			       AdvectDiffuse2D_Input_Parameters &IP);
  //@}

};

/*********************************************************************
 * AdvectDiffuse2D_Input_Parameters::get_cffc_path -- Get CFFC path. *
 *********************************************************************/
inline void AdvectDiffuse2D_Input_Parameters::get_cffc_path(){
  char *string_ptr;
  // Check to see if environment varible exists.
  if (getenv(PATHVAR_ADVECTDIFFUSE2D) == NULL) {
    //Set default path
     string_ptr = "CFFC";
     strcpy(CFFC_Path, string_ptr);
  } else {
     //Set path specified by environment variable
     strcpy(CFFC_Path, getenv(PATHVAR_ADVECTDIFFUSE2D));
  }
}

/***************************************************************
 * AdvectDiffuse2D_Input_Parameters -- Input-output operators. *
 ***************************************************************/
inline ostream &operator << (ostream &out_file,
			     const AdvectDiffuse2D_Input_Parameters &IP) {
    out_file << setprecision(6);
    out_file << "\n  -> CFFC Path: " 
	     << IP.CFFC_Path;
    out_file << "\n\n Solving 2D advection diffusion equations (IBVP/BVP) on multi-block solution-adaptive quadrilateral mesh.";
    out_file << "\n  -> Input File Name: " 
             << IP.Input_File_Name;
    if (IP.Time_Accurate) { 
       out_file << "\n  -> Time Accurate (Unsteady) Solution";
    } else {
       out_file << "\n  -> Time Invariant (Steady-State) Solution";
    }
    if (IP.Axisymmetric) { 
       out_file << "\n  -> 2D Axisymmetric Geometry";
    } else {
       out_file << "\n  -> 2D Planar Geometry";
    }
    out_file << "\n  -> Time Integration: " 
             << IP.Time_Integration_Type;
    out_file << "\n  -> Number of Stages in Multi-Stage Scheme: " 
             << IP.N_Stage;
    /*****************************************************
     *                     Multigrid                     */
    if (IP.i_Time_Integration == TIME_STEPPING_MULTIGRID) {
      out_file << IP.Multigrid_IP;
    }
    /*****************************************************/
    if (IP.Local_Time_Stepping == GLOBAL_TIME_STEPPING) {
      out_file << "\n  -> Global Time Stepping";
    } else if (IP.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
      out_file << "\n  -> Local Time Stepping";
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
    out_file << "\n  -> Initial Conditions: " 
             << IP.ICs_Type;
    out_file << "\n  -> Flow Velocity Field: " 
             << IP.Velocity_Field_Type; 
    out_file << "\n  -> Advection velocity x-direction: " 
             << IP.a;
    out_file << "\n  -> Advection velocity y-direction: " 
             << IP.b;
    out_file << "\n  -> Diffusion Coefficient : " 
             << IP.Kappa;
    out_file << "\n  -> Relaxation Time : " 
             << IP.Tau;
    out_file << "\n  -> Grid: " 
             << IP.Grid_Type;
    switch(IP.i_Grid) {
      case GRID_SQUARE :
        out_file << "\n  -> Size of Solution Domain: " 
                 << IP.Box_Width;
        break;
      case GRID_RECTANGULAR_BOX :
        out_file << "\n  -> Width of Solution Domain: " 
                 << IP.Box_Width;
        out_file << "\n  -> Height of Solution Domain: " 
                 << IP.Box_Height;
        break;
      case GRID_FLAT_PLATE :
        out_file << "\n  -> Plate Length: " 
                 << IP.Plate_Length;
        break;
      case GRID_PIPE :
        out_file << "\n  -> Pipe Length: " 
                 << IP.Pipe_Length;
        out_file << "\n  -> Pipe Radius: " 
                 << IP.Pipe_Radius;
        break;
      case GRID_BLUNT_BODY :
        out_file << "\n  -> Cylinder Radius: " 
                 << IP.Blunt_Body_Radius;
        break;
      case GRID_ROCKET_MOTOR :
	out_file << "\n  -> Length of Chamber (m): "
		 << IP.Chamber_Length;
	out_file << "\n  -> Radius of Chamber (m): "
		 << IP.Chamber_Radius;
	out_file << "\n  -> Distance from Chamber to Nozzle Throat (m): "
		 << IP.Chamber_To_Throat_Length;
	out_file << "\n  -> Length of the Nozzle (m): "
		 << IP.Nozzle_Length;
	out_file << "\n  -> Radius of the Nozzle at Throat (m): "
		 << IP.Nozzle_Radius_Throat;
	out_file << "\n  -> Radius of the Nozzle at Exit(m): "
		 << IP.Nozzle_Radius_Exit;
	out_file << "\n  -> Radius of the Propellant Grain (m): "
		 << IP.Grain_Radius;
	out_file << "\n  -> Nozzle type: "
		 << IP.Nozzle_Type;
        break;
      case GRID_CIRCULAR_CYLINDER :
        out_file << "\n  -> Cylinder Radius: " 
                 << IP.Cylinder_Radius;
        break;
      case GRID_ELLIPSE :
        out_file << "\n  -> Width of Ellipse along x-axis: " 
                 << IP.Ellipse_Length_X_Axis;
        out_file << "\n  -> Height of Ellipse along y-axis: " 
                 << IP.Ellipse_Length_Y_Axis;
        break;
      case GRID_NACA_AEROFOIL :
        out_file << "\n  -> NACA " 
                 << IP.NACA_Aerofoil_Type;
        out_file << "\n  -> Chord Length: " 
                 << IP.Chord_Length;
        break;
      case GRID_FREE_JET :
        out_file << "\n  -> Orifice Radius: " 
                 << IP.Orifice_Radius;
        break;
      case GRID_ICEMCFD :
        break;
      case GRID_READ_FROM_DEFINITION_FILE :
        break;
      case GRID_READ_FROM_GRID_DATA_FILE :
        break;
      default:
        out_file << "\n  -> Width of Solution Domain: " 
                 << IP.Box_Width;
        out_file << "\n  -> Height of Solution Domain: " 
                 << IP.Box_Height;
        break;
    } /* endswitch */
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
    out_file << "\n  -> Maximum Time: " 
             << IP.Time_Max;
    out_file << "\n  -> Maximum Number of Time Steps (Iterations): " 
             << IP.Maximum_Number_of_Time_Steps;
    out_file << "\n  -> Maximum Number of NKS Iterations: " 
             << IP.Maximum_Number_of_NKS_Iterations;
    out_file << "\n  -> Maximum Number of GMRES Iterations: " 
             << IP.Maximum_Number_of_GMRES_Iterations;
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
    return (out_file);
}

inline istream &operator >> (istream &in_file,
			     AdvectDiffuse2D_Input_Parameters &IP) {
    return (in_file);
}

/*************************************************************
 * AdvectDiffuse2D_Input_Parameters -- External subroutines. *
 *************************************************************/

extern void Open_Input_File(AdvectDiffuse2D_Input_Parameters &IP);

extern void Close_Input_File(AdvectDiffuse2D_Input_Parameters &IP);

extern void Set_Default_Input_Parameters(AdvectDiffuse2D_Input_Parameters &IP);

extern void Broadcast_Input_Parameters(AdvectDiffuse2D_Input_Parameters &IP);

#ifdef _MPI_VERSION
extern void Broadcast_Input_Parameters(AdvectDiffuse2D_Input_Parameters &IP,
                                       MPI::Intracomm &Communicator, 
                                       const int Source_CPU);
#endif

extern void Get_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters &IP);

extern int Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters &IP);

extern int Check_Input_Parameters(AdvectDiffuse2D_Input_Parameters &IP);

extern int Process_Input_Control_Parameter_File(AdvectDiffuse2D_Input_Parameters &Input_Parameters,
                                                char *Input_File_Name_ptr,
                                                int &Command_Flag);

#endif /* _ADVECTDIFFUSE2D_INPUT_INCLUDED  */
