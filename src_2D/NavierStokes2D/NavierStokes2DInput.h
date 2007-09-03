/**********************************************************************
 * NavierStokes2DInput.h: Header file defining 2D Navier-Stokes input *
 *                        parameter class.                            *
 **********************************************************************/

#ifndef _NAVIERSTOKES2D_INPUT_INCLUDED
#define _NAVIERSTOKES2D_INPUT_INCLUDED

// Include 2D NavierStokes state header file.

#ifndef _NAVIERSTOKES2D_STATE_INCLUDED
#include "NavierStokes2DState.h"
#endif // _NAVIERSTOKES2D_STATE_INCLUDED

// Include the math header file.

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

// Include the 2D cell header file.

#ifndef _CELL2D_INCLUDED
#include "../Grid/Cell2D.h"
#endif // _CELL2D_INCLUDED

// Include the 2D quadrilateral multiblock grid header file.

#ifndef _GRID2D_QUAD_BLOCK_INCLUDED
#include "../Grid/Grid2DQuad.h"
#endif // _GRID2D_QUAD_BLOCK_INCLUDED

// Include NASA rotor 37 and 67 header files.

#ifndef _NASA_ROTOR37_INCLUDED
#include "../Grid/NASARotor37.h"
#endif // _NASA_ROTOR37_INCLUDED

#ifndef _NASA_ROTOR67_INCLUDED
#include "../Grid/NASARotor67.h"
#endif // _NASA_ROTOR67_INCLUDED

// Include ICEMCFD input header file.

#ifndef _ICEMCFD_INCLUDED
#include "../ICEM/ICEMCFD.h"
#endif // _ICEMCFD_INCLUDED

// Include multigrid input header file.

#ifndef _FASMULTIGRID2DINPUT_INCLUDED
#include "../FASMultigrid2D/FASMultigrid2DInput.h"
#endif // _FASMULTIGRID2DINPUT_INCLUDED

// Include embedded boundary input header file.

#ifndef _EMBEDDEDBOUNDARIES2DINPUT_INCLUDED
#include "../Interface2D/EmbeddedBoundaries2D_Input.h"
#endif // _EMBEDDEDBOUNDARIES2DINPUT_INCLUDED
#include "../NewtonKrylovSchwarz2D/NKSInput2D.h"

// Define the structures and classes.

#define	INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D    128

// Enviroment flag for CFFC root directory path
#define PATHVAR_NAVIERSTOKES2D "CFFC_Path"

/*!
 * Class:  NavierStokes2D_Input_Parameters
 *
 * @brief Definition and manipulation of 2D Navier-Stokes input variables.
 *
 */
class NavierStokes2D_Input_Parameters{
private:
public:
  //@{ @name Input file parameters.
  //! CFFC root directory path:
  char CFFC_Path[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  //! Input file name:
  char Input_File_Name[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  //! Input file stream:
  ifstream Input_File;
  //! Input file line number:
  int Line_Number;
  //@}

  //@{ @name Time integration type indicator and related input parameters:
  char Time_Integration_Type[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  int i_Time_Integration;
  int Time_Accurate, Local_Time_Stepping, 
      Maximum_Number_of_Time_Steps, N_Stage;
  double CFL_Number, Time_Max;
  //! Residual variable:
  int i_Residual_Variable;
  //@}

  //@{ @name Implicit residual smoothing control parameters:
  int Residual_Smoothing;
  double Residual_Smoothing_Epsilon;
  int Residual_Smoothing_Gauss_Seidel_Iterations;
  //@}

  //@{ @name Multigrid related input parameters:
  Multigrid_Input_Parameters Multigrid_IP;
  //@}

  //@{ Newton-Krylov-Schwarz related input parameters:
  NKS_Input_Parameters NKS_IP;
  int Solver_Type; 
  //@}

  //@{ @name Reconstruction type indicator and related input parameters:
  char Reconstruction_Type[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  int i_Reconstruction;
  //@}

  //@{ @name Limiter type indicator and related input parameters:
  char Limiter_Type[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  int i_Limiter;
  int  Freeze_Limiter;
  double Freeze_Limiter_Residual_Level;
  //@}

  //@{ @name Flux function type and related input parameters:
  char Flux_Function_Type[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  int i_Flux_Function;
  char NavierStokes_Flux_Function_Type[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  int i_NavierStokes_Flux_Function;
  char Viscous_Reconstruction_Type[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  int i_Viscous_Reconstruction;
  //@}

  //@{ @name Initial condition type indicator and related input parameters:
  char ICs_Type[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  char Gas_Type[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  int i_ICs;
  NavierStokes2D_pState Wo, W1, W2;
  NavierStokes2D_cState Uo;
  double Pressure, Temperature, Mach_Number, Mach_Number2, Flow_Angle, Reynolds_Number, dp, Re_lid;
  Vector2D Wave_Position;
  double Wave_Width;
  //@}

  //@{ @name Propellant parameters:
  char Propellant_Type[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  //@}

  //@{ @name Flow type flag:
  char Flow_Type[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  int FlowType;
  //@}

  //@{ @name Flow geometry (planar or axisymmetric):
  char Flow_Geometry_Type[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  int Axisymmetric;
  //@}

  //@{ @name Turbulence parameters:
  char Turbulence_BC_Type[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  int i_Turbulence_BCs;
  char Friction_Velocity_Type[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  int i_Friction_Velocity;
  double C_constant, von_Karman_Constant;
  double yplus_sublayer, yplus_buffer_layer, yplus_outer_layer;
  int i_Turbulent_Wall_Injection, n_cells_sublayer_support;
  double sigmav, lw;
  //@}

  //@{ @name Grid type indicator and related input parameters:
  char Grid_Type[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  char NACA_Aerofoil_Type[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  int i_Grid;
  int Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Ghost_Cells,
      Number_of_Blocks_Idir, Number_of_Blocks_Jdir;
  double Box_Width, Box_Height, Plate_Length,
         Pipe_Length, Pipe_Radius,
         Blunt_Body_Radius, Blunt_Body_Mach_Number,
         Chamber_Length, Chamber_Radius, Chamber_To_Throat_Length,
         Nozzle_Length, Nozzle_Radius_Exit, Nozzle_Radius_Throat, Grain_Radius,
         Cylinder_Radius, Ellipse_Length_X_Axis, 
         Ellipse_Length_Y_Axis, Chord_Length, Orifice_Radius,
         Inner_Streamline_Number, Outer_Streamline_Number, Isotach_Line,
         Wedge_Angle, Wedge_Length,
         Step_Height, Top_Wall_Deflection;
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
  char Boundary_Conditions_Specified[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  int BCs_Specified;
  char BC_North_Type[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  char BC_South_Type[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  char BC_East_Type[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  char BC_West_Type[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  int BC_North, BC_South, BC_East, BC_West;
  //@}

  //@{ @name NASA rotor input variables:
  char NASA_Rotor37_Data_Directory[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D],
       NASA_Rotor67_Data_Directory[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  int Rotor_Flow_Type;
  double Rotor_Percent_Span;
  NASARotor37 NASA_Rotor37;
  NASARotor67 NASA_Rotor67;
  //@}

  //@{ @name Flow constants:
  Vector2D Vwall;
  double Twall;
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
  //! Number of interface mesh refinements.
  int Number_of_Interface_Mesh_Refinements;
  //! Number of bounding-box mesh refinements.
  int Number_of_Bounding_Box_Mesh_Refinements;
  //! Number of flat plate mesh refinements.
  int Number_of_Flat_Plate_Mesh_Refinements;
  //! Interface refinement condition.
  int Interface_Refinement_Condition;
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
  //! Gradient of the gas density field refinement criteria.
  int Refinement_Criteria_Gradient_Density;
  //! Divergence of the gas velocity field refinement criteria.
  int Refinement_Criteria_Divergence_Velocity;
  //! Curl of the gas velocity field refinement criteria.
  int Refinement_Criteria_Curl_Velocity;
  //! Gradient of the turbulence kinetic energy refinement criteria.
  int Refinement_Criteria_Gradient_Turbulence_Kinetic_Energy;
  //! Bounding-box for bounding-box mesh refinement.
  Vector2D AMR_Xmin, AMR_Xmax;
  //! Morton Re-Ordering
  int Morton, Morton_Reordering_Frequency;
  //! Smooth quad block flag:
  int i_Smooth_Quad_Block;
  //@}

  //@{ @name Embedded boundary input parameters:
  EmbeddedBoundaries2D_Input_Parameters Interface_IP;
  int Reset_Interface_Motion_Type;
  //@}

  //@{ @name Output parameters:
  //! Output file name:
  char Output_File_Name[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  //! Multi-block mesh definition input file names:
  char Grid_File_Name[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  char Grid_Definition_File_Name[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  //! Restart file name:
  char Restart_File_Name[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  //! Gnuplot file name:
  char Gnuplot_File_Name[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  //! Next_Control_Parameter:
  char Next_Control_Parameter[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
  //! Output format type indicator:
  char Output_Format_Type[INPUT_PARAMETER_LENGTH_NAVIERSTOKES2D];
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
  friend ostream &operator << (ostream &out_file, const NavierStokes2D_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file, NavierStokes2D_Input_Parameters &IP);
  //@}

};

/*************************************************************
 * NavierStokes2D_Input_Parameters::get_cffc_path -- Get CFFC path. *
 *************************************************************/
inline void NavierStokes2D_Input_Parameters::get_cffc_path(){
  char *string_ptr;
  // Check to see if environment varible exists.
  if (getenv(PATHVAR_NAVIERSTOKES2D) == NULL) {
    //Set default path
     string_ptr = "CFFC";
     strcpy(CFFC_Path, string_ptr);
  } else {
     //Set path specified by environment variable
     strcpy(CFFC_Path, getenv(PATHVAR_NAVIERSTOKES2D));
  }
}

/**********************************************************************
 * NavierStokes2D_Input_Parameters -- Input-output operators.         *
 **********************************************************************/
inline ostream &operator << (ostream &out_file,
			     const NavierStokes2D_Input_Parameters &IP) {
  out_file << setprecision(6);
  out_file << "\n  -> CFFC Path: " 
	   << IP.CFFC_Path;
  if (IP.i_Grid == GRID_CARTESIAN_UNIFORM) {
    out_file << "\n\n Solving 2D Navier-Stokes equations (IBVP/BVP) on uniform Cartesian mesh.";
  } else {
    out_file << "\n\n Solving 2D Navier-Stokes equations (IBVP/BVP) on multi-block solution-adaptive quadrilateral mesh.";
  }
  out_file << "\n  -> Input File Name: " << IP.Input_File_Name;
  if (IP.Time_Accurate) { 
    out_file << "\n  -> Time Accurate (Unsteady) Solution";
  } else {
    out_file << "\n  -> Time Invariant (Steady-State) Solution";
  }
  if (IP.FlowType == FLOWTYPE_INVISCID) {
    out_file << "\n  -> 2D Euler equations (inviscid flow)";   
  } else if (IP.FlowType == FLOWTYPE_LAMINAR) {
    out_file << "\n  -> 2D Laminar Navier-Stokes equations";
  } else if (IP.FlowType == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    out_file << "\n  -> 2D Turbulent Navier-Stokes equations with the k-omega turbulence model";
  } else {
    out_file << "\n  -> 2D Inviscid Flow";   
  }
  if (IP.FlowType == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    out_file << "\n  -> Turbulence boundary condition: " << IP.Turbulence_BC_Type;
    out_file << "\n  -> Friction velocity evaluation type: " << IP.Friction_Velocity_Type;
    out_file << "\n  -> Wall surface constant, C: " << IP.C_constant
	     << "\n  -> von Karman constant, kappa: " << IP.von_Karman_Constant;
    out_file << "\n  -> yplus sublayer: " << IP.yplus_sublayer
	     << ", yplus buffer layer: " << IP.yplus_buffer_layer
	     << ", yplus outer layer: " << IP.yplus_outer_layer;
    if (IP.i_Turbulent_Wall_Injection)
      out_file << "\n  -> Turbulence wall injection is on.";
  }
  if (!IP.Axisymmetric) out_file << "\n  -> 2D Planar Flow";
  else out_file << "\n  -> 2D Axisymmetric Flow";
  out_file << "\n  -> Time Integration: " << IP.Time_Integration_Type;
  out_file << "\n  -> Number of Stages in Multi-Stage Scheme: " << IP.N_Stage;
  if (IP.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
      IP.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
    out_file << IP.Multigrid_IP;
  }
  if (IP.Local_Time_Stepping == GLOBAL_TIME_STEPPING) {
    out_file << "\n  -> Global Time Stepping";
  } else if (IP.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
    out_file << "\n  -> Scalar Local Time Stepping";
  } else if (IP.Local_Time_Stepping == SEMI_IMPLICIT_LOCAL_TIME_STEPPING) {
    out_file << "\n  -> Semi-Implicit Local Time Stepping";
  }
  out_file << "\n  -> L1-, L2-, and max-norms computed on residual variable: " << IP.i_Residual_Variable;
  if (IP.Residual_Smoothing) {
    out_file << "\n  -> Residual Smoothing:";
    out_file << "\n     -> Epsilon: " << IP.Residual_Smoothing_Epsilon;
    out_file << "\n     -> Gauss_Seidel_Iterations: " 
	     << IP.Residual_Smoothing_Gauss_Seidel_Iterations;
  }
  out_file << "\n  -> Reconstruction: " << IP.Reconstruction_Type;
  out_file << "\n  -> Limiter: " << IP.Limiter_Type;
  if (IP.Limiter_Type != LIMITER_ZERO && IP.Freeze_Limiter)
    out_file << "\n  -> Freeze Limiter when L2-norm of residual is < "
	     << IP.Freeze_Limiter_Residual_Level;
  out_file << "\n  -> Hyperbolic Flux Function: " << IP.Flux_Function_Type;
  if (IP.FlowType)
    out_file << "\n  -> Elliptic Flux Evaluation: " << IP.Viscous_Reconstruction_Type;
  out_file << "\n  -> Initial Conditions: " << IP.ICs_Type;
  out_file << "\n  -> Gas Type: " << IP.Gas_Type;
  if (IP.i_Grid == GRID_ROCKET_MOTOR) {
    out_file << "\n  -> Solid propellant type: " << IP.Propellant_Type;
    out_file << "\n     -> Propellant density = " << IP.Wo.rhos;
    out_file << "\n     -> Propellant burning constant = " << IP.Wo.n;
    out_file << "\n     -> Propellant burning coefficient = " << IP.Wo.beta;
    out_file << "\n     -> Propellant flame temperature = " << IP.Wo.Tf;
    out_file << "\n     -> Propellant surface temperature = " << IP.Wo.Ts;
  }
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
    out_file << "\n  -> Pressure (kPa): " << IP.Pressure/THOUSAND;
    out_file << "\n  -> Temperature (K): " << IP.Temperature;
    out_file << "\n  -> Mach Number: " << IP.Mach_Number;
    out_file << "\n  -> Flow Angle (degrees): " << IP.Flow_Angle;
    break;
  case IC_VISCOUS_CHANNEL_FLOW :
    out_file << "\n  -> Pressure Change (Pa): " << IP.dp;
    out_file << "\n  -> Upper Wall Speed (m/s): " << IP.Vwall.x;
    break;
  case IC_VISCOUS_PIPE_FLOW :
    out_file << "\n  -> Pressure Change (Pa): " << IP.dp;
    break;
  case IC_VISCOUS_STOKES_FLOW :
    out_file << "\n  -> Stokes flow Reynolds number: " << IP.Reynolds_Number;
    break;
  case IC_VISCOUS_FLAT_PLATE :
    out_file << "\n  -> Flat plate free stream Mach number: " << IP.Mach_Number;
    out_file << "\n  -> Flat plate free stream Reynolds number: " << IP.Reynolds_Number;
    break;
  case IC_VISCOUS_DRIVEN_CAVITY_FLOW :
    out_file << "\n  -> Driven cavity flow with lid speed: " << IP.Vwall.x
	     << " and Reynolds number: " << IP.Re_lid;
    break;
  default:
    break;
  }
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
  case GRID_FLAT_PLATE :
    out_file << "\n  -> Plate Length (m): " << IP.Plate_Length;
    break;
  case GRID_PIPE :
    out_file << "\n  -> Pipe Length (m): " << IP.Pipe_Length;
    out_file << "\n  -> Pipe Radius (m): " << IP.Pipe_Radius;
    break;
  case GRID_BLUNT_BODY :
    out_file << "\n  -> Cylinder Radius (m): " << IP.Blunt_Body_Radius;
    out_file << "\n  -> Blunt Body Mach Number: " << IP.Blunt_Body_Mach_Number;
    break;
  case GRID_ROCKET_MOTOR :
    out_file << "\n  -> Length of Chamber (m): " << IP.Chamber_Length;
    out_file << "\n  -> Radius of Chamber (m): " << IP.Chamber_Radius;
    out_file << "\n  -> Distance from Chamber to Nozzle Throat (m): " << IP.Chamber_To_Throat_Length;
    out_file << "\n  -> Length of the Nozzle (m): " << IP.Nozzle_Length;
    out_file << "\n  -> Radius of the Nozzle at Throat (m): " << IP.Nozzle_Radius_Throat;
    out_file << "\n  -> Radius of the Nozzle at Exit(m): " << IP.Nozzle_Radius_Exit;
    out_file << "\n  -> Radius of the Propellant Grain (m): " << IP.Grain_Radius;
    out_file << "\n  -> Nozzle type: " << IP.Nozzle_Type;
    break;
  case GRID_CIRCULAR_CYLINDER :
    out_file << "\n  -> Cylinder Radius (m): " << IP.Cylinder_Radius;
    break;
  case GRID_ELLIPSE :
    out_file << "\n  -> Width of Ellipse along x-axis (m): " 
	     << IP.Ellipse_Length_X_Axis;
    out_file << "\n  -> Height of Ellipse along y-axis (m): " 
	     << IP.Ellipse_Length_Y_Axis;
    break;
  case GRID_NACA_AEROFOIL :
    out_file << "\n  -> NACA " << IP.NACA_Aerofoil_Type;
    out_file << "\n  -> Chord Length (m): " << IP.Chord_Length;
    break;
  case GRID_FREE_JET :
    out_file << "\n  -> Orifice Radius (m): " << IP.Orifice_Radius;
    break;
  case GRID_RINGLEB_FLOW :
    out_file << "\n  -> Inner streamline number: " << IP.Inner_Streamline_Number;
    out_file << "\n  -> Outer streamline number: " << IP.Outer_Streamline_Number;
    out_file << "\n  -> Isotach line: " << IP.Isotach_Line;
    break;
  case GRID_WEDGE :
    out_file << "\n  -> Wedge Angle (degrees): " << IP.Wedge_Angle;
    out_file << "\n  -> Wedge Length (m): " << IP.Wedge_Length;
    break;
  case GRID_UNSTEADY_BLUNT_BODY :
    out_file << "\n  -> Cylinder Radius (m): " << IP.Blunt_Body_Radius;
    break;
  case GRID_BUMP_CHANNEL_FLOW :
    if (strcmp(IP.Grid_Type,"Bump_Channel_Flow") == 0) {
      if (IP.Smooth_Bump) out_file << " (smooth bump)";
      else out_file << " (non-smooth bump)";
    }
    break;
  case GRID_BACKWARD_FACING_STEP :
    out_file << "\n  -> Step height: " << IP.Step_Height;
    out_file << "\n  -> Reynolds number: " << IP.Reynolds_Number;
    out_file << "\n  -> Inflow height: " << IP.Top_Wall_Deflection;
    out_file << "\n  -> Maximum inflow velocity x-component: " << IP.Wo.v.x;
    break;
  case GRID_MIXING_LAYER:
    out_file << "\n  -> Reynolds Number: " << IP.Reynolds_Number;
    out_file << "\n  -> Mach Number: "<<IP.Mach_Number;
    out_file << "\n  -> Mach Number2: "<<IP.Mach_Number2;
    break;
  case GRID_NASA_ROTOR_37 :
  case GRID_NASA_ROTOR_67 :
    out_file << "\n  -> Percent Span: " << IP.Rotor_Percent_Span;
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
  out_file << "\n  -> Mesh Stretching: "
	   << IP.i_Mesh_Stretching;
  out_file << "\n  -> Mesh Stretching Type Idir: "
	   << IP.Mesh_Stretching_Type_Idir;
  out_file << "\n  -> Mesh Stretching Type Jdir: "
	   << IP.Mesh_Stretching_Type_Jdir;
  out_file << "\n  -> Mesh Stretching Factor Idir: "
	   << IP.Mesh_Stretching_Factor_Idir;
  out_file << "\n  -> Mesh Stretching Factor Jdir: "
	   << IP.Mesh_Stretching_Factor_Jdir;
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
  out_file << "\n  -> Number of Initial Mesh Refinements: " 
	   << IP.Number_of_Initial_Mesh_Refinements;
  out_file << "\n  -> Number of Uniform Mesh Refinements: " 
	   << IP.Number_of_Uniform_Mesh_Refinements;
  out_file << "\n  -> Number of Boundary Mesh Refinements: " 
	   << IP.Number_of_Boundary_Mesh_Refinements;
  out_file << "\n  -> Number of Interface Mesh Refinements: " 
	   << IP.Number_of_Interface_Mesh_Refinements;
  out_file << "\n  -> Number of Bounding-Box Mesh Refinements: " 
	   << IP.Number_of_Bounding_Box_Mesh_Refinements;
  out_file << "\n  -> Number of Flat-Plate Mesh Refinements: " 
	   << IP.Number_of_Flat_Plate_Mesh_Refinements;
  out_file << "\n  -> Refinement Criteria: ";
  if (IP.Refinement_Criteria_Gradient_Density)
  out_file << "\n     -> Gradient of the density field refinement criteria";
  if (IP.Refinement_Criteria_Divergence_Velocity)
  out_file << "\n     -> Divergence of the velocity field refinement criteria";
  if (IP.Refinement_Criteria_Curl_Velocity)
  out_file << "\n     -> Curl of the velocity field refinement criteria";
  if (IP.Refinement_Criteria_Gradient_Turbulence_Kinetic_Energy)
  out_file << "\n     -> Gradient of the turbulence kinetic energy.";
  if (IP.Number_of_Bounding_Box_Mesh_Refinements)
    out_file << "\n     -> Bounding-box for bounding-box AMR:" << IP.AMR_Xmin << IP.AMR_Xmax;
  out_file << "\n  -> Smooth Quad Block: "
	   << IP.i_Smooth_Quad_Block;
  out_file << "\n  -> CFL Number: "
	   << IP.CFL_Number;
  out_file << "\n  -> Maximum Time (ms): " 
	   << IP.Time_Max*THOUSAND;
  out_file << "\n  -> Maximum Number of Time Steps (Iterations): " 
	   << IP.Maximum_Number_of_Time_Steps;
  if (IP.NKS_IP.Maximum_Number_of_NKS_Iterations > 0) {
  out_file << "\n  -> Maximum Number of NKS Iterations: " 
             << IP.NKS_IP.Maximum_Number_of_NKS_Iterations;
  }
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

  if (IP.Morton) {
    out_file << "\n  -> Morton Re-Ordering Frequency: "
      << IP.Morton_Reordering_Frequency
      << " time steps (iterations)";
  }
  
  return out_file;
}

inline istream &operator >> (istream &in_file,NavierStokes2D_Input_Parameters &IP) {
  return in_file;
}

/**********************************************************************
 * NavierStokes2D_Input_Parameters -- External subroutines.           *
 **********************************************************************/

extern void Open_Input_File(NavierStokes2D_Input_Parameters &IP);

extern void Close_Input_File(NavierStokes2D_Input_Parameters &IP);

extern void Set_Default_Input_Parameters(NavierStokes2D_Input_Parameters &IP);

extern void Broadcast_Input_Parameters(NavierStokes2D_Input_Parameters &IP);

#ifdef _MPI_VERSION
extern void Broadcast_Input_Parameters(NavierStokes2D_Input_Parameters &IP,
                                       MPI::Intracomm &Communicator,
                                       const int Source_CPU);
#endif

extern void Get_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters &IP);

extern int Parse_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters &IP);

extern int Process_Input_Control_Parameter_File(NavierStokes2D_Input_Parameters &IP,
                                                char *Input_File_Name_ptr,
                                                int &Command_Flag);

extern void Initialize_Reference_State(NavierStokes2D_Input_Parameters &IP);

extern void Reinitialize_Reference_State(NavierStokes2D_Input_Parameters &IP);

#endif // _NAVIERSTOKES2D_INPUT_INCLUDED
