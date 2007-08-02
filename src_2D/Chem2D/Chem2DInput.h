/****************** Chem2DInput.h ************************************
  This class defines the all the input parameters read from the 
  standard input and the input parameter file.  It is based on 
  the Euler2DInput class and extended for the Chem2D flow
  solver.  The main changes are the addition of:
           
   -> Mulitple species
         -User Defined
   -> Reaction Mechanisms
         -User Defined
   -> Thermally Perfect Thermodynamic Properties
          -P, T, rho, Vel, etc... ie functions of Temperature
   -> Viscosity
   -> Gravity
   -> Low Mach Number Precondtioning 
   -> BCs NO_SLIP, MOVING_WALL, SUBSONIC_INFLOW, SUBSONIC_OUTFLOW 
 
  - based on Euler2DInput.h 
***********************************************************************/


#ifndef _CHEM2D_INPUT_INCLUDED
#define _CHEM2D_INPUT_INCLUDED

/* Include 2D Euler state, 2D cell, 2D quadrilateral multiblock 
   grid, and NASA rotor header files. */

#include <cstdlib> 

using namespace std;

#ifndef _CHEM2D_STATE_INCLUDED
#include "Chem2DState.h"
#endif // _CHEM2D_STATE_INCLUDED

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

// Include ICEMCFD input header file.

#ifndef _ICEMCFD_INCLUDED
#include "../ICEM/ICEMCFD.h"
#endif // _ICEMCFD_INCLUDED

/* Define the structures and classes. */
#define	INPUT_PARAMETER_LENGTH_CHEM2D    128

//Enviroment Flag 
#define PATHVAR "CFDkit_Path"

/*!
 * Class:  Chem2D_Input_Parameters
 *
 * @brief Definition and manipulation of 2D chemical reacting gas input
 *        variables.
 *
 */
class Chem2D_Input_Parameters{
  private:
  public:
  //@{ @name Input file parameters.
  //! Input file name:
  char Input_File_Name[INPUT_PARAMETER_LENGTH_CHEM2D];
  //! Input file stream:
  ifstream Input_File;
  //! Input file line number:
  int Line_Number;
  //@}

  //@{ @name Time integration type indicator and related input parameters:
  char Time_Integration_Type[INPUT_PARAMETER_LENGTH_CHEM2D];
  int i_Time_Integration;
  int Time_Accurate, Local_Time_Stepping, 
      Maximum_Number_of_Time_Steps, N_Stage;
  double CFL_Number, Time_Max;
  // Residual variable:
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

  //@{ @name Reconstruction type indicator and related input parameters:
  char Reconstruction_Type[INPUT_PARAMETER_LENGTH_CHEM2D];
  int i_Reconstruction;
  //@}

  //@{ @name Limiter type indicator and related input parameters:
  char Limiter_Type[INPUT_PARAMETER_LENGTH_CHEM2D];
  int i_Limiter;
  int  Freeze_Limiter;
  double Freeze_Limiter_Residual_Level;
  //@}

  //@{ @name Flux function type and related input parameters:
  char Flux_Function_Type[INPUT_PARAMETER_LENGTH_CHEM2D];
  int i_Flux_Function;
  char Viscous_Flux_Evaluation_Type[INPUT_PARAMETER_LENGTH_CHEM2D];
  int i_Viscous_Flux_Evaluation;
  //@}

  //@{ @name Initial condition type indicator and related input parameters:
  char ICs_Type[INPUT_PARAMETER_LENGTH_CHEM2D];
  int i_ICs;
  //Pa, K, m/s, degress from horizontal
  double Pressure, Temperature, Flow_Angle;
  double Mach_Number_Reference,Mach_Number;
  int Preconditioning;
  //@}

  //@{ @name Chemical reacting flow inpput parameters.
  string react_name;
  char React_Name[INPUT_PARAMETER_LENGTH_CHEM2D];
  char **Multispecies;
  //! Array of species names
  string *multispecies;
  //! Array of mass fractions
  double *mass_fractions;
  int num_species;
  Chem2D_pState Wo;
  Chem2D_cState Uo;
  double Global_Schmidt; // Depricated, use each individual Schmidt's
  //! Individual Schmidt number for each species.
  double *Schmidt;

  //! Root path of CFDkit+caboodle 
  char CFDkit_Path[INPUT_PARAMETER_LENGTH_CHEM2D];
  void get_cfdkit_path();
  //@}

  //@{ @name Flow type indicator and related input parameters:
  char Flow_Type[INPUT_PARAMETER_LENGTH_CHEM2D];
  int FlowType;
  //@}

  //@{ @name Flow geometry (planar or axisymmetric):
  char Flow_Geometry_Type[INPUT_PARAMETER_LENGTH_CHEM2D];
  int Axisymmetric; //0 no, 1 yes
  //@}

  //@{ @name Gravity indicator (yes/no) = (1,0).
  int Gravity;
  //@}

  //@{ @name Turbulence parameters:
  char Turbulence_BC_Type[INPUT_PARAMETER_LENGTH_CHEM2D];
  int i_Turbulence_BCs;
  char Friction_Velocity_Type[INPUT_PARAMETER_LENGTH_CHEM2D];
  int i_Friction_Velocity;
  double C_constant, von_Karman_Constant;
  double yplus_sublayer, yplus_buffer_layer, yplus_outer_layer;
  //@}

  //@{ @name Grid type indicator and related input parameters:
  char Grid_Type[INPUT_PARAMETER_LENGTH_CHEM2D];
  char NACA_Aerofoil_Type[INPUT_PARAMETER_LENGTH_CHEM2D];
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
         Wedge_Angle, Wedge_Length;
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
  char Boundary_Conditions_Specified[INPUT_PARAMETER_LENGTH_CHEM2D];
  int BCs_Specified;
  char BC_North_Type[INPUT_PARAMETER_LENGTH_CHEM2D];
  char BC_South_Type[INPUT_PARAMETER_LENGTH_CHEM2D];
  char BC_East_Type[INPUT_PARAMETER_LENGTH_CHEM2D];
  char BC_West_Type[INPUT_PARAMETER_LENGTH_CHEM2D];
  int BC_North, BC_South, BC_East, BC_West;
  //@}

  //@{ @name Flow constants:
  double Moving_wall_velocity;
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
  //! Gradient of the density field refinement criteria.
  int Refinement_Criteria_Gradient_Density;
  //! Divergence of the velocity field refinement criteria.
  int Refinement_Criteria_Divergence_Velocity;
  //! Curl of the velocity field refinement criteria.
  int Refinement_Criteria_Curl_Velocity;
  //! Gradient of the temperature field refinement criteria.
  int Refinement_Criteria_Gradient_Temperature;
  //! Gradient of the CH4 mass fraction field refinement criteria.
  int Refinement_Criteria_Gradient_CH4;
  //! Gradient of the CO2 mass fraction field refinement criteria.
  int Refinement_Criteria_Gradient_CO2;
  //! Smooth quad block flag:
  int i_Smooth_Quad_Block;
  //@}

  //@{ @name Debug Level 0 for none, 1,2,3... level of verboseness
  int debug_level;
  //@}

  //@{ @name Output parameters:
  //! Output file name:
  char Output_File_Name[INPUT_PARAMETER_LENGTH_CHEM2D];
  //! Multi-block mesh definition input file names:
  char Grid_File_Name[INPUT_PARAMETER_LENGTH_CHEM2D];
  char Grid_Definition_File_Name[INPUT_PARAMETER_LENGTH_CHEM2D];
  //! Restart file name:
  char Restart_File_Name[INPUT_PARAMETER_LENGTH_CHEM2D];
  //! Gnuplot file name:
  char Gnuplot_File_Name[INPUT_PARAMETER_LENGTH_CHEM2D];
  //! Next_Control_Parameter:
  char Next_Control_Parameter[INPUT_PARAMETER_LENGTH_CHEM2D];
  //! Output format type indicator:
  char Output_Format_Type[INPUT_PARAMETER_LENGTH_CHEM2D];
  int i_Output_Format;
  int Restart_Solution_Save_Frequency;
  int Time_Accurate_Plot_Frequency;
  //! Output progress frequency:
  int Output_Progress_Frequency;
  //@}

  //@{ @name Multi-block solution-adaption and parallel domain decomposition input parameters:
  int Number_of_Processors, Number_of_Blocks_Per_Processor;
  //@}

  //@{ @name Constructor and destructor.
  //! Default constructor.
  Chem2D_Input_Parameters(){
    Multispecies = NULL;
    multispecies = NULL; 
    mass_fractions = NULL;
    ICEMCFD_FileNames = NULL;
  }

  //! Default constructor.
  ~Chem2D_Input_Parameters(){
    for (int i = 0; i < 3; i++) {
      delete[] ICEMCFD_FileNames[i]; 
    } /* endfor */
    delete[]  ICEMCFD_FileNames; ICEMCFD_FileNames=NULL;
    //State class memory
    Deallocate(); 
  }
  //@}

  //@{ @name Allocation and deallocation.
  //! Allocate Memory.
  void Allocate();
  //!  Deallocate Memory.
  void Deallocate();
  //@}

  //@{ @name Input-output operators:
  friend ostream &operator << (ostream &out_file, const Chem2D_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file, Chem2D_Input_Parameters &IP);
  //@}

};

/*************************************************************
 * Chem2D_Input_Parameters -- Memory Management              *
 *************************************************************/
inline void Chem2D_Input_Parameters::Allocate() {
  Multispecies = new char*[num_species];
  for (int i = 0; i < num_species; i++) {
    Multispecies[i] = new char[INPUT_PARAMETER_LENGTH_CHEM2D];
  } 
  multispecies = new string[num_species]; 
  mass_fractions = new double[num_species];
  Schmidt = new double[num_species];
}

inline void Chem2D_Input_Parameters::Deallocate() {
  if(Multispecies != NULL ){
    for (int i = 0; i < num_species; i++) {
      delete[] Multispecies[i]; Multispecies[i]=NULL;
    }
    delete[] Multispecies; Multispecies=NULL;
  }
  
  if(multispecies != NULL){
    delete[] multispecies; multispecies=NULL;
  }
  
  if(mass_fractions != NULL){
    delete[] mass_fractions; mass_fractions=NULL;
    if( Schmidt != NULL) delete[] Schmidt; Schmidt = NULL;
    Wo.Deallocate_static(); 
    Uo.Deallocate_static();
  }
 
} 

/*************************************************************
 * Chem2D_Input_Parameters -- Input-output operators.        *
 *************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Chem2D_Input_Parameters &IP) {
    out_file << setprecision(6);
  
    /*********************************************************/
    out_file << "\n\n Solving 2D MulitSpecies";
    if (IP.FlowType ==  FLOWTYPE_INVISCID){
      out_file<<" Euler (Inviscid) ";
    } else {
      out_file<<" Navier-Stokes (Viscous) ";
    }
    out_file<<"equations (IBVP/BVP) "; 
    if (IP.i_Grid == GRID_CARTESIAN_UNIFORM) {
      out_file<<"\n on uniform Cartesian mesh.";
    } else { 
      out_file <<"\n on multi-block solution-adaptive quadrilateral mesh.";
    } 

    /*********************************************************/
    if (IP.FlowType ==  FLOWTYPE_INVISCID) {
    } else if (IP.FlowType == FLOWTYPE_LAMINAR) {
      out_file << "\n  -> Laminar flow";
    } else if (IP.FlowType == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
      out_file << "\n  -> Turbulent flow: RANS with k-epsilon turbulence model";
      if (IP.i_Turbulence_BCs){
         out_file <<" (wall function formulation)";
      } else {
         out_file <<" (low-Reynolds-number formulation)";
      }
    } else if (IP.FlowType == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      out_file << "\n  -> Turbulent flow: RANS with k-oemga turbulence model";
      if (IP.i_Turbulence_BCs) {
         out_file <<" (wall function formulation)";
      } else {
         out_file <<" (low-Reynolds-number formulation)";
      }
    } else if (IP.FlowType == FLOWTYPE_TURBULENT_LES) {
      out_file << "\n  -> Turbulent flow: LES ";
    } else if (IP.FlowType == FLOWTYPE_TURBULENT_DES_K_OMEGA) {
      out_file << "\n  -> Turbulent flow: DES with k-omega SGS turbulence model ";
    } else if (IP.FlowType == FLOWTYPE_TURBULENT_DNS) {
      out_file << "\n  -> Turbulent flow: DNS ";
    }
    if (IP.FlowType == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      out_file << "\n  -> Turbulence boundary condition: " << IP.Turbulence_BC_Type;
      out_file << "\n  -> Friction velocity evaluation type: " << IP.Friction_Velocity_Type;
      out_file << "\n  -> Wall surface constant, C: " << IP.C_constant
	     << "\n  -> von Karman constant, kappa: " << IP.von_Karman_Constant;
      out_file << "\n  -> yplus sublayer: " << IP.yplus_sublayer
	       << ", yplus buffer layer: " << IP.yplus_buffer_layer
	       << ", yplus outer layer: " << IP.yplus_outer_layer;
    }

    /*********************************************************/
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
    /********************************************************/
    if(IP.Gravity) {
      out_file << "\n  -> With Gravity (-z)";
    }
    if(IP.FlowType != FLOWTYPE_INVISCID) {
      out_file << "\n  -> Viscous Reconstruction Method: "
	       << IP.Viscous_Flux_Evaluation_Type;
    }
    if(IP.debug_level){
      out_file << "\n  -> Debug level: "
	       << IP.debug_level;
    }
    /********************************************************/
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
    if (IP.Local_Time_Stepping == GLOBAL_TIME_STEPPING) {
      out_file << "\n  -> Global Time Stepping";
    } else if (IP.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
      out_file << "\n  -> Scalar Local Time Stepping";
    } else if (IP.Local_Time_Stepping == MATRIX_LOCAL_TIME_STEPPING) {
      out_file << "\n  -> Matrix Local Time Stepping";
    } else if (IP.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) {
      out_file << "\n  -> Low-Mach-Number Local Preconditioning (Weiss-Smith)";
    } else if (IP.Local_Time_Stepping == SEMI_IMPLICIT_LOCAL_TIME_STEPPING) {
      out_file << "\n  -> Semi-Implicit Local Time Stepping";
    } else if (IP.Local_Time_Stepping == SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER) {
      out_file << "\n  -> Semi-Implicit Low-Mach-Number Local Preconditioned Time Stepping";
    }

    /*********************************************************/
    if (IP.Preconditioning == 1){
	out_file <<"\n  -> Mach Number Reference "<<IP.Mach_Number_Reference;
    }

    /*********************************************************/
    out_file << "\n  -> L1-, L2-, and max-norms computed on residual variable: " << IP.i_Residual_Variable;

    /*********************************************************/
    if (IP.Residual_Smoothing) {
      out_file << "\n  -> Residual Smoothing";
      out_file << "\n  -> Epsilon: " 
               << IP.Residual_Smoothing_Epsilon;
      out_file << "\n  -> Gauss_Seidel_Iterations: " 
               << IP.Residual_Smoothing_Gauss_Seidel_Iterations;
    }

    /*****************************************************/
    out_file << "\n  -> Reconstruction: " 
             << IP.Reconstruction_Type;
    out_file << "\n  -> Limiter: " 
             << IP.Limiter_Type;   
    if (IP.Limiter_Type != LIMITER_ZERO && IP.Freeze_Limiter) {
      out_file << "\n  -> Freeze Limiter when L2-norm of residual is < "
	       << IP.Freeze_Limiter_Residual_Level;
    } 
    out_file << "\n  -> Flux Function: " 
             << IP.Flux_Function_Type;
    /********** CHEM2D ****************************/
    out_file << "\n  -> Reaction Mechanism: " 
	     << IP.react_name;
    out_file << "\n  -> Species: "<<IP.Wo.ns
	     << "\n  -> Initial mass fractions: ";
    for(int i=0; i<IP.Wo.ns; i++){
      out_file  <<"c["<<IP.multispecies[i]<<"]= ";
      out_file  << IP.Wo.spec[i].c<<", ";
    } 
    if(IP.FlowType != FLOWTYPE_INVISCID){
      out_file << "\n  -> Schmidt Numbers for Viscous flow: ";
      for(int i=0; i <IP.Wo.ns; i++){
	out_file  <<"Sc["<<IP.multispecies[i]<<"]= ";
	out_file << IP.Schmidt[i]<<", ";
      }
    }
    out_file << "\n  -> CFDkit+caboodle Path: " 
	     << IP.CFDkit_Path;
    /*********************************************/ 
    out_file << "\n  -> Initial Conditions: " 
             << IP.ICs_Type;
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
	break;
      case IC_RIEMANN_XDIR :
	break;
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
	/******** CHEM2D ********/
    case IC_GAS_MIX :
      break;
    case IC_CHEM_CORE_FLAME:
      break;
    case IC_CHEM_INVERSE_FLAME:
      break;
    case IC_CHEM_1DFLAME:
      break;
    case IC_PRESSURE_GRADIENT_X:
      break;
    case IC_PRESSURE_GRADIENT_Y:
      break;
    case IC_VISCOUS_CHANNEL_FLOW:
      break;
      /********** CHEM2D *******/
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
      case GRID_COUETTE :
	out_file << "\n  -> Width of Solution Domain (m): " 
		 << IP.Box_Width;
	out_file << "\n  -> Height of Solution Domain (m): " 
		 << IP.Box_Height;
	out_file << "\n  -> Moving Wall Velocity (m/s): "
		 << IP.Moving_wall_velocity;
	break;
      case GRID_1DFLAME :
	out_file << "\n  -> Width of Solution Domain (m): " 
		 << IP.Box_Width;
	out_file << "\n  -> Height of Solution Domain (m): " 
		 << IP.Box_Height;
	break;
      case GRID_LAMINAR_FLAME :
	out_file << "\n  -> Width of Solution Domain (m): " 
		 << IP.Pipe_Length;
	out_file << "\n  -> Height of Solution Domain (m): " 
		 << IP.Pipe_Radius;
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
      case GRID_BUMP_CHANNEL_FLOW :
	if (strcmp(IP.Grid_Type,"Bump_Channel_Flow") == 0) {
	  if (IP.Smooth_Bump) out_file << " (smooth bump)";
	  else out_file << " (non-smooth bump)";
	}
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
    if (IP.Refinement_Criteria_Gradient_Density)
    out_file << "\n     -> Gradient of the density field refinement criteria";
    if (IP.Refinement_Criteria_Divergence_Velocity)
    out_file << "\n     -> Divergence of the velocity field refinement criteria";
    if (IP.Refinement_Criteria_Curl_Velocity)
    out_file << "\n     -> Curl of the velocity field refinement criteria";
    if (IP.Refinement_Criteria_Gradient_Temperature)
    out_file << "\n     -> Gradient of the temperature field refinement criteria";
    if (IP.Refinement_Criteria_Gradient_CH4)
    out_file << "\n     -> Gradient of the CH4 mass fraction field refinement criteria";
    if (IP.Refinement_Criteria_Gradient_CO2)
    out_file << "\n     -> Gradient of the CO2 mass fraction field refinement criteria";
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
    out_file << "\n  -> Output Progress Frequency: "
	     << IP.Output_Progress_Frequency
	     << " steps (iterations)";
    //CHEM2D
    if (IP.Time_Accurate_Plot_Frequency !=0 && IP.Time_Accurate){
      out_file << "\n  -> Time Accurate Solution Plot Frequency: "
	       << IP.Time_Accurate_Plot_Frequency
	       << " steps (iterations)"; 
    }
    return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Chem2D_Input_Parameters &IP) {
    return (in_file);
}

/****************** Get CFD kit Path *********/
inline void Chem2D_Input_Parameters::get_cfdkit_path(){
 
  // Check to see if environment varible exists.
  if(getenv(PATHVAR) == NULL){
    cerr<<"Missing Environment Variable: "<<PATHVAR<<"\n\n"; 
    exit(1);
  }
  
  //Set Path
  strcpy(CFDkit_Path,getenv(PATHVAR));
}

/**************** DESTRUCTOR ******************/
//inline Chem2D_Input_Parameters::Chem2D_Input_Parameters(){
//  delete[] mass_fractions;
//}

/*************************************************************
 * Chem2D_Input_Parameters -- External subroutines.         *
 *************************************************************/

extern void Open_Input_File(Chem2D_Input_Parameters &IP);

extern void Close_Input_File(Chem2D_Input_Parameters &IP);

extern void Set_Default_Input_Parameters(Chem2D_Input_Parameters &IP);

extern void Broadcast_Input_Parameters(Chem2D_Input_Parameters &IP);

#ifdef _MPI_VERSION
extern void Broadcast_Input_Parameters(Chem2D_Input_Parameters &IP,
                                       MPI::Intracomm &Communicator, 
                                       const int Source_CPU);
#endif

extern void Get_Next_Input_Control_Parameter(Chem2D_Input_Parameters &IP);

extern int Parse_Next_Input_Control_Parameter(Chem2D_Input_Parameters &IP);

extern int Process_Input_Control_Parameter_File(Chem2D_Input_Parameters &Input_Parameters,
                                                char *Input_File_Name_ptr,
                                                int &Command_Flag);

extern void Equivalence_Ratio(const double &phi);

#endif /* _CHEM2D_INPUT_INCLUDED  */
