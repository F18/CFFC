/****************** Flame2DInput.h ************************************
  This class defines the all the input parameters read from the 
  standard input and the input parameter file.  It is based on 
  the Euler2DInput class and extended for the Flame2D flow
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
 
NEW

  - based on Euler2DInput.h 
***********************************************************************/

#ifndef _FLAME2D_INPUT_INCLUDED
#define _FLAME2D_INPUT_INCLUDED

/* Include 2D Euler state, 2D cell, 2D quadrilateral multiblock 
   grid, and NASA rotor header files. */
#include "../Grid/Cell2D.h"
#include "../Grid/Grid2DQuad.h"

/* Also include multigrid input header file. */
#include "../FASMultigrid2D/FASMultigrid2DInput.h"

/* Also include NKS  input header file. */
#include "../NewtonKrylovSchwarz2D/NKSInput2D.h"

#include "Flame2DState.h"

/////////////////////////////////////////////////////////////////////
/// Defines
/////////////////////////////////////////////////////////////////////

/* Define the structures and classes. */
#undef INPUT_PARAMETER_LENGTH_FLAME2D
#define	INPUT_PARAMETER_LENGTH_FLAME2D    128

//Enviroment Flag 
#undef PATHVAR
#define PATHVAR "CFFC_Path"

//! Composition flag:
//! This flag species what method was used to specify 
//! the reference composition.
enum CompType { FLAME2D_INPUT_MASS_FRACTIONS,
		FLAME2D_INPUT_MOLE_FRACTIONS,
		FLAME2D_INPUT_EQUIVALENCE_RATIO };


/////////////////////////////////////////////////////////////////////
/// Function Definitions
/////////////////////////////////////////////////////////////////////

/*!
 * Class:  Flame2D_Input_Parameters
 *
 * @brief Definition and manipulation of 2D chemical reacting gas input
 *        variables.
 *
 */
class Flame2D_Input_Parameters{
private:
public:
  //@{ @name Input file parameters.
  //! Input file name:
  char Input_File_Name[INPUT_PARAMETER_LENGTH_FLAME2D];
  //! Input file stream:
  ifstream Input_File;
  //! Input file line number:
  int Line_Number;
  //@}

  // Input file line number, the index of solution parameters:
  int Solution_Parameters_Index;
  
  // Time integration type indicator and related input parameters:
  char Time_Integration_Type[INPUT_PARAMETER_LENGTH_FLAME2D];
  int i_Time_Integration;
  int Time_Accurate, Local_Time_Stepping, 
    Maximum_Number_of_Time_Steps, N_Stage;
  double CFL_Number, Time_Max, Time_Step;
  bool Fixed_Time_Step;

  //! source term time-step size multiplyer
  double Source_Term_Multiplyer;
  
  // Additional input parameters for dual time stepping
  int     Max_Inner_Steps, first_step;
  double  dTime, Physical_CFL_Number;

  //@{ @name Implicit residual smoothing control parameters:
  int Residual_Smoothing;
  double Residual_Smoothing_Epsilon;
  int Residual_Smoothing_Gauss_Seidel_Iterations;
  //@}

  //@{ @name Multigrid related input parameters:
  Multigrid_Input_Parameters Multigrid_IP;
  //@}

  //@{ @name NKS related input parametrs:
  NKS_Input_Parameters  NKS_IP;
  int Solver_Type; 
  //@}

  //@{ @name Reconstruction type indicator and related input parameters:
  char Reconstruction_Type[INPUT_PARAMETER_LENGTH_FLAME2D];
  int i_Reconstruction;
  //@}

  //@{ @name Limiter type indicator and related input parameters:
  char Limiter_Type[INPUT_PARAMETER_LENGTH_FLAME2D];
  int i_Limiter;
  int  Freeze_Limiter;
  int i_Residual_Variable;
  int Number_of_Residual_Norms;
  double Freeze_Limiter_Residual_Level;
  //@}

  //@{ @name Flux function type and related input parameters:
  char Flux_Function_Type[INPUT_PARAMETER_LENGTH_FLAME2D];
  int i_Flux_Function;
  char Viscous_Flux_Evaluation_Type[INPUT_PARAMETER_LENGTH_FLAME2D];
  int i_Viscous_Flux_Evaluation;
  //@}

  //@{ @name Initial condition type indicator and related input parameters:
  char ICs_Type[INPUT_PARAMETER_LENGTH_FLAME2D];
  int i_ICs;
  int i_Grid_Level;
  
  //Pa, K, m/s, degress from horizontal
  double Pressure, Temperature, Flow_Angle;
  double Mach_Number_Reference,Mach_Number_Reference_target,Mach_Number,Re_lid;
  int Preconditioning,Dual_Time_Stepping;
  //@}

  //@{ @name Chemical reacting flow inpput parameters.
  Flame2D_pState Wo;
  //! Number of species
  int num_species;
  //! Array of mass fractions
  double *mass_fractions;
  //! Fuel species
  char Fuel_Species[INPUT_PARAMETER_LENGTH_FLAME2D];
  //! Fuel species
  double equivalence_ratio;
  //!laminar flame speed [m/s]
  double flame_speed;
  //! flags for flow properties
  bool reacting, constant_schmidt;
  //! temporary storage
  string schmidt_string, composition_string;
  int i_specified_composition;
  

  //! BC Pressure Gradient 
  double Pressure_Gradient;

  //! Root path of CFFC 
  char CFFC_Path[INPUT_PARAMETER_LENGTH_FLAME2D];
  void get_cffc_path();
  //@}

  //@{ @name Cantera input parameters.
  //! Mechanism name
  string ct_mech_name;
  char ct_Mech_Name[INPUT_PARAMETER_LENGTH_FLAME2D];
  //! Mechanism file
  string ct_mech_file;
  char ct_Mech_File[INPUT_PARAMETER_LENGTH_FLAME2D];
  //@}

  //@{ @name Flow type indicator and related input parameters:
  char Flow_Type[INPUT_PARAMETER_LENGTH_FLAME2D];
  int FlowType;
  //@}

  //@{ @name Flow geometry (planar or axisymmetric):
  char Flow_Geometry_Type[INPUT_PARAMETER_LENGTH_FLAME2D];
  int Axisymmetric; //0 no, 1 yes
  //@}

  //@{ @name Flow Conditions:
  double Global_Schmidt;  //depricated, use each individual Schmidt's
  double *Schmidt;        //individual for each species
  //@}

  //@{ @Debug Level 0 for none, 1,2,3... level of verboseness  
  int debug_level;
  //@}

  //@{ @name Gravity indicator (yes/no) = (1,0).
  int Gravity;
  double gravity_z; // acceration due to gravity [m/s] (negative => acts downward)
  //@}


  //@{ @name Radiation related input parameters.
  //! Input file name:
  char Rte_Input_File_Name[INPUT_PARAMETER_LENGTH_FLAME2D];
  //! Radiation Flag (0 -> no radiation, 1 -> radiation modelled)
  int Radiation;
  //! Number of sequential solves
  int Max_Number_Sequential_Solves;
  //@}

  //@{ @name Grid type indicator and related input parameters:
  char Grid_Type[INPUT_PARAMETER_LENGTH_FLAME2D];
  char NACA_Aerofoil_Type[INPUT_PARAMETER_LENGTH_FLAME2D];
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
    Wedge_Angle, Wedge_Length, Heat_Source; 

  //! Bluff Body
  double Length_Shroud,Radius_Shroud,Length_BluffBody,Radius_BluffBody,
    Radius_Orifice, Length_Inlet_Pipe, Radius_Inlet_Pipe, Length_Combustor_Tube, 
    Radius_Combustor_Tube;
  double BluffBody_Coflow_Air_Velocity, BluffBody_Coflow_Fuel_Velocity;
  int BluffBody_Data_Usage; // 0 no, 1 yes, 

  int Smooth_Bump, Nozzle_Type;
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

  //@{ @name Boundary conditions:
  char Boundary_Conditions_Specified[INPUT_PARAMETER_LENGTH_FLAME2D];
  int BCs_Specified;
  char BC_North_Type[INPUT_PARAMETER_LENGTH_FLAME2D];
  char BC_South_Type[INPUT_PARAMETER_LENGTH_FLAME2D];
  char BC_East_Type[INPUT_PARAMETER_LENGTH_FLAME2D];
  char BC_West_Type[INPUT_PARAMETER_LENGTH_FLAME2D];
  int BC_North, BC_South, BC_East, BC_West;
  //@}

  //@{ @name Flow constants:
  double Moving_wall_velocity;
  //@}

  //@{ @name Morton Re-Ordering
  int Morton, Morton_Reordering_Frequency;
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

  //@{ @name Output parameters:
  //! Output file name:
  char Output_File_Name[INPUT_PARAMETER_LENGTH_FLAME2D];
  //! Multi-block mesh definition input file names:
  char Grid_File_Name[INPUT_PARAMETER_LENGTH_FLAME2D];
  char Grid_Definition_File_Name[INPUT_PARAMETER_LENGTH_FLAME2D];
  //! Restart file name:
  char Restart_File_Name[INPUT_PARAMETER_LENGTH_FLAME2D];
  //! Gnuplot file name:
  char Gnuplot_File_Name[INPUT_PARAMETER_LENGTH_FLAME2D];
  //! Next_Control_Parameter:
  char Next_Control_Parameter[INPUT_PARAMETER_LENGTH_FLAME2D];
  //! Output format type indicator:
  char Output_Format_Type[INPUT_PARAMETER_LENGTH_FLAME2D];
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
  Flame2D_Input_Parameters(){
    Schmidt = NULL; 
    mass_fractions = NULL;
  }

  //! Default constructor.
  ~Flame2D_Input_Parameters(){
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

  //! set reference state
  void setRefSolutionState(void);

  //@{ @name Input-output operators:
  friend ostream &operator << (ostream &out_file, const Flame2D_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file, Flame2D_Input_Parameters &IP);
  //@}

};

/*************************************************************
 * Flame2D_Input_Parameters -- Memory Management              *
 *************************************************************/
inline void Flame2D_Input_Parameters::Allocate() {
  mass_fractions = new double[num_species];
  Schmidt = new double[num_species];
}

inline void Flame2D_Input_Parameters::Deallocate() {
  if(mass_fractions != NULL){
    delete[] mass_fractions; mass_fractions=NULL;
    if( Schmidt != NULL) delete[] Schmidt; Schmidt = NULL;
    Wo.DeallocateStatic(); 
  }
 
} 


/*************************************************************
 * Flame2D_Input_Parameters -- Setup ref solution state.      *
 *************************************************************/
inline void Flame2D_Input_Parameters::setRefSolutionState(void) {

  // load the mechanism
  Flame2D_pState::setMixture(ct_mech_name, ct_mech_file);
  
  // setup static
  Flame2D_pState::set_Mref(Mach_Number_Reference );
  Flame2D_pState::set_gravity(gravity_z);
  if (!reacting) Flame2D_pState::setNonReacting();

  //-----------------------------------------------------------------
  // for constant schmidt
  //-----------------------------------------------------------------
  if (constant_schmidt) {
    if (!schmidt_string.empty()) 
      Mixture::parse_schmidt_string(schmidt_string, Schmidt);
    Flame2D_pState::setConstantSchmidt(Schmidt);
  }

  //-----------------------------------------------------------------
  // Compute the species mass fractions, depending upon
  // what was specified
  //-----------------------------------------------------------------
  // -> Mass fractions
  if (i_specified_composition == FLAME2D_INPUT_MASS_FRACTIONS) {
    
    if (!composition_string.empty()) 
      Mixture::parse_mass_string(composition_string, mass_fractions);      

    // -> Mole fractions
  } else if (i_specified_composition == FLAME2D_INPUT_MOLE_FRACTIONS) {

    // store mole fractions in mass_fractions array
    Mixture::parse_mole_string(composition_string, mass_fractions);

    // convert molar fractions to mass fractions
    double sum(0.0);
    for(int i=0; i<num_species; i++) 
      sum += mass_fractions[i]*Mixture::molarMass(i);
    for(int i=0; i<num_species; i++) 
      mass_fractions[i] = mass_fractions[i]*Mixture::molarMass(i)/sum;

    // -> Equivalence Ratio
  } else if (i_specified_composition == FLAME2D_INPUT_EQUIVALENCE_RATIO) {
    Mixture::composition(Fuel_Species, 
			 equivalence_ratio, 
			 mass_fractions);
  } // endif - composition
    

    //-----------------------------------------------------------------
    // set the solution state
    //-----------------------------------------------------------------
  Wo.setState_TPY(Temperature, Pressure, mass_fractions);
  Wo.setVelocity( Mach_Number*Wo.a()*cos(TWO*PI*Flow_Angle/360.00),
		  Mach_Number*Wo.a()*sin(TWO*PI*Flow_Angle/360.00) );

  // Now that the mass_fractions and schmidt numbers have been set,
  // force the flags
  i_specified_composition = FLAME2D_INPUT_MASS_FRACTIONS;
  schmidt_string = "";
  composition_string = "";

}


/*************************************************************
 * Flame2D_Input_Parameters -- Input-output operators.        *
 *************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Flame2D_Input_Parameters &IP) {
  out_file << setprecision(6);
  out_file << "\n  -> CFFC Path: " 
	   << IP.CFFC_Path;
  
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
  if (IP.Radiation)
    out_file<<"\n Radiation heat transfer modelled.";
  else
    out_file<<"\n Radiation Neglected.";

  /*********************************************************/
  if (IP.FlowType ==  FLOWTYPE_INVISCID) {
  } else if (IP.FlowType == FLOWTYPE_LAMINAR) {
    out_file << "\n  -> Laminar flow";
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
    out_file << "\n  -> With Gravity (-z): g_z="<<IP.gravity_z<<" [m/s]";
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
  } else if (IP.Local_Time_Stepping == MATRIX_WITH_LOW_MACH_NUMBER_PRECONDITIONER) {
    out_file << "\n  -> Matrix Local Time Stepping with Low-Mach-Number Local Preconditioning";
  } else if (IP.Local_Time_Stepping == DUAL_TIME_STEPPING) {
    out_file << "\n  -> Dual Time Stepping";
  } else if (IP.Local_Time_Stepping == DUAL_LOW_MACH_NUMBER_PRECONDITIONER) {
    out_file << "\n  -> Dual Low-Mach-Number Local Preconditioned Time Stepping";
  } else if (IP.Local_Time_Stepping == DUAL_SEMI_IMPLICIT_LOW_MACH_NUMBER_PRECONDITIONER) {
    out_file << "\n  -> Dual Semi-Implicit Low-Mach-Number Local Preconditioned Time Stepping";
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
  out_file << "\n  -> Residual L2 Norm of: "; 
  if (IP.i_Residual_Variable == 1) out_file <<" density ";
  if (IP.i_Residual_Variable == 2 || IP.i_Residual_Variable == 3) out_file <<" momentum ";
  if (IP.i_Residual_Variable == 4) out_file <<" energy ";
  out_file << "\n  -> Flux Function: " 
	   << IP.Flux_Function_Type;
  /********** FLAME2D ****************************/
  out_file << "\n  -> Reaction Mechanism: " 
	   << Flame2D_pState::mechName();
  out_file << "\n  -> Species: "<<IP.Wo.NumSpecies()
	   << "\n  -> Initial mass fractions: ";
  for(int i=0; i<IP.Wo.NumSpecies(); i++){
    out_file  <<"c["<<Flame2D_pState::speciesName(i)<<"]= ";
    out_file  << IP.Wo.c(i)<<", ";
  } 
  if(IP.FlowType != FLOWTYPE_INVISCID){
    out_file << "\n  -> Schmidt Numbers for Viscous flow: ";
    for(int i=0; i <IP.Wo.NumSpecies(); i++){
      out_file  <<"Sc["<<Flame2D_pState::speciesName(i)<<"]= ";
      out_file << IP.Schmidt[i]<<", ";
    }
  }
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
    out_file << "\n  -> Heat_source : " 
	     << IP.Heat_Source;
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
  case IC_VISCOUS_DRIVEN_CAVITY_FLOW :
    out_file << "\n  -> Driven cavity flow with lid speed: " << IP.Moving_wall_velocity
	     << " and Reynolds number: " << IP.Re_lid;
    break;
    /******** FLAME2D ********/
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
  case IC_VISCOUS_COUETTE:
    break;
  case IC_VISCOUS_FLAT_PLATE :
    out_file << "\n  -> Viscous flat plate free stream Mach number: " << IP.Mach_Number;
    out_file << "\n  -> Viscous flat plate Reynolds number: "
	     <<IP.Wo.Re(IP.Plate_Length);
    break;
    /********** FLAME2D *******/
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
  case GRID_TEST :
    out_file << "\n  -> Width of Solution Domain (m): " 
	     << IP.Box_Width;
    out_file << "\n  -> Height of Solution Domain (m): " 
	     << IP.Box_Height;
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
  case GRID_DRIVEN_CAVITY_FLOW:
    out_file << "\n  -> Cavity Depth (m): " 
	     << IP.Box_Width; 
    //        out_file <<" \n  -> Cavity Grid Stretch level "<<IP.Stretch_Level;
    break;
       
  case GRID_BLUNT_BODY :
    out_file << "\n  -> Cylinder Radius (m): " 
	     << IP.Blunt_Body_Radius;
    break;
  case GRID_BLUFF_BODY :
    out_file << "\n  -> Shroud Length (m): " 
	     << IP.Length_Shroud;
    out_file << "\n  -> Shroud Radius (m): " 
	     << IP.Radius_Shroud;
    out_file << "\n  -> Bluff Body Length (m): " 
	     << IP.Length_BluffBody;
    out_file << "\n  -> Bluff Body Radius (m): " 
	     << IP.Radius_BluffBody;
    out_file << "\n  -> Fuel Orifice Radius (m): " 
	     << IP.Radius_Orifice;
    break;
  case GRID_DUMP_COMBUSTOR :
    out_file << "\n  -> Inlet Pipe Length (m): " 
	     << IP.Length_Inlet_Pipe;
    out_file << "\n  -> Inlet Pipe Radius (m): " 
	     << IP.Radius_Inlet_Pipe;
    out_file << "\n  -> Combustor Tube Length (m): " 
	     << IP.Length_Combustor_Tube;
    out_file << "\n  -> Combustor Tube Radius(m): " 
	     << IP.Radius_Combustor_Tube;
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
  if (IP.Dual_Time_Stepping)
    out_file << "\n  -> Physical CFL Number: "
	     << IP.Physical_CFL_Number; 
  out_file << "\n  -> Maximum Time (ms): " 
	   << IP.Time_Max*THOUSAND;
  out_file << "\n  -> Maximum Number of Explicit Time Steps (Iterations): " 
	   << IP.Maximum_Number_of_Time_Steps;
  out_file << "\n  -> Maximum Number of Implicit (NKS) Steps (Iterations): " 
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
  if (IP.Time_Accurate_Plot_Frequency !=0 && IP.Time_Accurate){
    out_file << "\n  -> Time Accurate Solution Plot Frequency: "
	     << IP.Time_Accurate_Plot_Frequency
	     << " steps (iterations)"; 
  }
  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Flame2D_Input_Parameters &IP) {
  return (in_file);
}


/****************** Get CFFC Path *********/
inline void Flame2D_Input_Parameters::get_cffc_path(){
  // Check to see if environment varible exists.
  if(getenv(PATHVAR) == NULL){
    cerr<<"Missing Environment Variable: "<<PATHVAR<<"\n\n"; 
    exit(1);
  }
  //Set Path
  strcpy(CFFC_Path,getenv(PATHVAR));
}


/**************** DESTRUCTOR ******************/
//inline Flame2D_Input_Parameters::Flame2D_Input_Parameters(){
//  delete[] mass_fractions;
//}

/*************************************************************
 * Flame2D_Input_Parameters -- External subroutines.         *
 *************************************************************/

extern void Open_Input_File(Flame2D_Input_Parameters &IP);

extern void Close_Input_File(Flame2D_Input_Parameters &IP);

extern void Set_Default_Input_Parameters(Flame2D_Input_Parameters &IP);

extern void Broadcast_Input_Parameters(Flame2D_Input_Parameters &IP);

extern void Get_Next_Input_Control_Parameter(Flame2D_Input_Parameters &IP);

extern int Parse_Next_Input_Control_Parameter(Flame2D_Input_Parameters &IP);

extern int Process_Input_Control_Parameter_File(Flame2D_Input_Parameters &Input_Parameters,
						char *Input_File_Name_ptr,
						int &Command_Flag);

extern void Equivalence_Ratio(const double &phi);


#endif /* _FLAME2D_INPUT_INCLUDED  */
