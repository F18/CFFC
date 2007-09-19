/******************* Rte2DInput.h ************************************
  This class defines the all the input parameters read from the 
  standard input and the input parameter file.  It is based on 
  the Euler2DInput class and extended for the Rte2D 
  solver.  The main changes are the addition of:
           
   -> Type of radiating gas
         -User Defined
   -> Angular discretization 
         -User Defined
   -> Absorbsion coefficient
   -> Scattering model

***********************************************************************/


#ifndef _RTE2D_INPUT_INCLUDED
#define _RTE2D_INPUT_INCLUDED

/* Class Declaration */
class Rte2D_Input_Parameters;

/* Include 2D Rte state, 2D cell, 2D quadrilateral multiblock 
   grid, and NASA rotor header files. */

#include "Rte2DState.h"
#include "Rte2DTools.h"
#include "../Grid/Cell2D.h"
#include "../Grid/Grid2DQuad.h"
#include "../Grid/NASARotor37.h"
#include "../Grid/NASARotor67.h"

/* Also include multigrid input header file. */
#include "../FASMultigrid2D/FASMultigrid2DInput.h"

/* Also include NKS  input header file. */
#include "../NewtonKrylovSchwarz2D/NKSInput2D.h"

/* Include ICEMCFD input header file.*/
#include "../ICEM/ICEMCFD.h"

/* Define the structures and classes. */

#define	INPUT_PARAMETER_LENGTH_RTE2D    128

//Enviroment Flag 
#define PATHVAR "CFFC_Path"


/*!
 * Class: Rte2D_Input_Parameters
 *
 * @brief Definition and manipulation of 2D Rte input variables.
 *
 */
class Rte2D_Input_Parameters{
  private:
  public:
  //@{ @name Input file parameters.
  //! Input file name:
  char Input_File_Name[INPUT_PARAMETER_LENGTH_RTE2D];
  //! Input file stream:
  ifstream Input_File;
  //! Input file line number:
  int Line_Number;
  //@}

  //@{ @name CFFC path parameters/functions.
  //! Root path of CFFC
  char CFFC_Path[INPUT_PARAMETER_LENGTH_RTE2D];
  //! Get root path of CFFC
  void get_cffc_path();
  //@}

  //@{ @name Time integration type indicator and related input parameters:
  char Time_Integration_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  int i_Time_Integration;
  int Time_Accurate, Local_Time_Stepping, 
      Maximum_Number_of_Time_Steps, N_Stage;
  double CFL_Number, Time_Max;
  double Mr_Min_Factor;
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

  //@{ @name Newton-Krylov-Schwarz related input parameters:
  NKS_Input_Parameters  NKS_IP;
  int Solver_Type;
  //@}

  //@{ @name SNBCK related input parameters
  SNBCK_Input_Parameters SNBCK_IP;
  //@}


  //@{ @name RTE equation solver and space-marching related parameters::
  int i_RTE_Solver;
  char RTE_Solver[INPUT_PARAMETER_LENGTH_RTE2D];
  char SpaceMarch_Scheme[INPUT_PARAMETER_LENGTH_RTE2D];
  int i_SpaceMarch_Scheme;
  int i_DOM_Quadrature;
  char DOM_Quadrature[INPUT_PARAMETER_LENGTH_RTE2D];
  //@}

  //@{ @name Reconstruction type indicator and related input parameters:
  char Reconstruction_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  int i_Reconstruction;
  //@}

  //@{ @name Limiter type indicator and related input parameters:
  char Limiter_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  int i_Limiter;
  int  Freeze_Limiter;
  double Freeze_Limiter_Residual_Level;
  //@}

  //@{ @name Flux function type and related input parameters:
  char Flux_Function_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  int i_Flux_Function;
  //@}

  //@{ @name Initial condition type indicator and related input parameters:
  char ICs_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  char Gas_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  int i_ICs;
  Rte2D_State Uo;
  Medium2D_State Mo;
  double Intensity;
  double Temperature, Pressure;
  double xco, xco2, xh2o, xo2, fsoot;
  double AbsorptionCoef, ScatteringCoef;
  int Number_of_Angles_Mdir, Number_of_Angles_Ldir;
  double NorthWallTemp, SouthWallTemp, EastWallTemp, WestWallTemp;
  double NorthWallEmiss, SouthWallEmiss, EastWallEmiss, WestWallEmiss;
  int i_ScatteringFunc;
  char ScatteringFunc[INPUT_PARAMETER_LENGTH_RTE2D];
  int i_AbsorptionModel;
  char AbsorptionModel[INPUT_PARAMETER_LENGTH_RTE2D];
  char ICs_Medium[INPUT_PARAMETER_LENGTH_RTE2D];
  int i_ICs_Medium;
  int Medium_Field_Type;
  //@}

  //@{ @name Flow geometry (planar or axisymmetric):
  char Flow_Geometry_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  int Axisymmetric;
  //@}

  //@{ @name Grid type indicator and related input parameters:
  char Grid_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  char NACA_Aerofoil_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  int i_Grid;
  int Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Ghost_Cells,
      Number_of_Blocks_Idir, Number_of_Blocks_Jdir; 
  double Box_Width, Box_Height, Plate_Length, 
         Pipe_Length, Pipe_Radius, 
         Blunt_Body_Radius, Blunt_Body_Mach_Number,
         Chamber_Length, Chamber_Radius, Chamber_To_Throat_Length,
         Nozzle_Length, Nozzle_Radius_Exit, Nozzle_Radius_Throat, Grain_Radius,
         Cylinder_Radius, Cylinder_Radius2, Ellipse_Length_X_Axis, 
         Ellipse_Length_Y_Axis, Chord_Length, Orifice_Radius,
         Inner_Streamline_Number, Outer_Streamline_Number, Isotach_Line,
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
  char Boundary_Conditions_Specified[INPUT_PARAMETER_LENGTH_RTE2D];
  int BCs_Specified;
  char BC_North_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  char BC_South_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  char BC_East_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  char BC_West_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  int BC_North, BC_South, BC_East, BC_West;
  //@}

  //@{ @name NASA rotor input variables:
  char NASA_Rotor37_Data_Directory[INPUT_PARAMETER_LENGTH_RTE2D],
       NASA_Rotor67_Data_Directory[INPUT_PARAMETER_LENGTH_RTE2D];
  int Rotor_Flow_Type;
  double Rotor_Percent_Span;
  NASARotor37 NASA_Rotor37;
  NASARotor67 NASA_Rotor67;
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
  //! Number of bounding-box mesh refinements:
  int Number_of_Bounding_Box_Mesh_Refinements;
  //! Interface refinement condition.
  int Interface_Refinement_Condition;
  //! Maximum number of refinement levels.
  int Maximum_Refinement_Level;
  //! Minimum number of refinement levels.
  int Minimum_Refinement_Level;
  //! Threshold for refinement of mesh.
  double Threshold_for_Refinement;
  //! Threshold for coarsening of mesh.
  double Threshold_for_Coarsening;
  //! Number of refinement criteria.
  int Number_of_Refinement_Criteria;
  //! Gradient of the density field refinement criteria.
  int Refinement_Criteria_Gradient_Density;
  //! Divergence of the velocity field refinement criteria.
  int Refinement_Criteria_Divergence_Velocity;
  //! Curl of the velocity field refinement criteria.
  int Refinement_Criteria_Curl_Velocity;
  //! Bounding-box for bounding-box mesh refinement.
  Vector2D AMR_Xmin, AMR_Xmax;
  //! Smooth quad block flag:
  int i_Smooth_Quad_Block;
  //@}

  //@{ @name Morton re-ordering parameters:
  int Morton, Morton_Reordering_Frequency;
  //@}

  //@{ @name Output parameters:
  //! Output file name:
  char Output_File_Name[INPUT_PARAMETER_LENGTH_RTE2D];
  //! Multi-block mesh definition input file names:
  char Grid_File_Name[INPUT_PARAMETER_LENGTH_RTE2D];
  char Grid_Definition_File_Name[INPUT_PARAMETER_LENGTH_RTE2D];
  //! Restart file name:
  char Restart_File_Name[INPUT_PARAMETER_LENGTH_RTE2D];
  //! Gnuplot file name:
  char Gnuplot_File_Name[INPUT_PARAMETER_LENGTH_RTE2D];
  //! Next_Control_Parameter:
  char Next_Control_Parameter[INPUT_PARAMETER_LENGTH_RTE2D];
  //! Output format type indicator:
  char Output_Format_Type[INPUT_PARAMETER_LENGTH_RTE2D];
  int i_Output_Format;
  int Restart_Solution_Save_Frequency;
  //! Output progress frequency:
  int Output_Progress_Frequency;
  //@}

  //@{ @name Multi-block solution-adaption and parallel domain decomposition input parameters:
  int Number_of_Processors, Number_of_Blocks_Per_Processor;
  //@}


  //@{ @name Member function to setup conserved state Uo and medium state Mo
  void SetupInputState();


  //@{ @name Input-output operators:
  friend ostream &operator << (ostream &out_file, const Rte2D_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file, Rte2D_Input_Parameters &IP);
  //@}

};

/*******************************************************
 * Main function to setup the Rte2D_State and          *
 * Medium2D_State static parameters.  This function is *
 * called everytime input parameters are changed.      *
 *******************************************************/  
inline void Rte2D_Input_Parameters::SetupInputState() 
{

  // deallocate
  Mo.Deallocate();
  Uo.Deallocate();
  
  // Setup static state variables
  Medium2D_State:: SetupStatic( i_AbsorptionModel, 
				i_ScatteringFunc,
				SNBCK_IP,
				CFFC_Path );
  Rte2D_State::SetupStatic( i_RTE_Solver, 
			    Medium2D_State::Nband,
			    Number_of_Angles_Mdir,
			    Number_of_Angles_Ldir,
			    i_DOM_Quadrature,
			    Axisymmetric,
			    i_ScatteringFunc,
			    CFFC_Path );

  // allocate
  Mo.Allocate();
  Uo.Allocate();
  
  // initialize
  Mo.SetInitialValues(Pressure,
		      Temperature,
		      xco,
		      xh2o,
		      xco2,
		      xo2,
		      fsoot,
		      AbsorptionCoef,
		      ScatteringCoef);
  Uo.SetInitialValues(Intensity);

  //
  // setup Medium2D fields
  //
  // CONSTANT
  if (i_ICs_Medium == IC_UNIFORM  || i_ICs_Medium == IC_CONSTANT) {
    Medium2D_State::SetConstantField( Mo );
    
  // DISCONTINUOUS
  } else if (i_ICs_Medium == IC_DISCONTINUOUS) {
    Medium2D_State::SetDiscontinuousField( Mo, 
					   Mo/FIVE, 
					   Vector2D(-0.5, -0.5), 
					   Vector2D(0.5, 0.5) );
    
  // ERROR
  } else {
    cerr << "Rte2D_Input_Parameters::SetupInputState - Invalid flag for field type\n";
    exit(-1);
  } // endif

}


/*************************************************************
 * Rte2D_Input_Parameters -- Input-output operators.         *
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
    out_file << "\n  -> Time Integration: " 
             << IP.Time_Integration_Type;
    out_file << "\n  -> Number of Stages in Multi-Stage Scheme: " 
             << IP.N_Stage;
    /*****************************************************
     *                     Multigrid                     */
    if (IP.i_Time_Integration == TIME_STEPPING_MULTIGRID) {
      out_file << IP.Multigrid_IP;
    }
    /*****************************************************
     *                   Space March                     */
    else if (IP.i_Time_Integration == TIME_STEPPING_SPACE_MARCH) {
      out_file << "\n  -> Space Marching";
      out_file << "\n  -> Differencing scheme: "
	       << IP.SpaceMarch_Scheme;
    }
    /*****************************************************/
    else if (IP.Local_Time_Stepping == GLOBAL_TIME_STEPPING) {
      out_file << "\n  -> Global Time Stepping";
    } else if (IP.Local_Time_Stepping == SCALAR_LOCAL_TIME_STEPPING) {
      out_file << "\n  -> Scalar Local Time Stepping";
    } else if (IP.Local_Time_Stepping == MATRIX_LOCAL_TIME_STEPPING) {
      out_file << "\n  -> Matrix Local Time Stepping";
    } else if (IP.Local_Time_Stepping == LOW_MACH_NUMBER_WEISS_SMITH_PRECONDITIONER) {
      out_file << "\n  -> Low-Mach-Number Local Preconditioning (Weiss-Smith)";
    } /* endif */
    // out_file << "\n  -> L1-, L2-, and max-norms computed on residual variable: " << IP.i_Residual_Variable;
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
    

    /***********************************************************************
     *************************** RTE SPECIFIC ******************************/
    // SNBCK parameters
    if (IP.i_AbsorptionModel == MEDIUM2D_ABSORB_SNBCK) {
      out_file << endl;
      IP.SNBCK_IP.Output(out_file);
    } else {
      out_file << endl << endl << string(75,'*') << endl;
    }

    out_file << string(29,'*') << "   RTE SPECIFIC  " << string(29,'*') << endl;
    out_file << string(75,'*');
   
    // space marching ?
    if (IP.i_Time_Integration == TIME_STEPPING_SPACE_MARCH) {
      out_file << "\n  -> Space Marching";
      out_file << "\n      -> Differencing scheme: "
	       << IP.SpaceMarch_Scheme;
    }

    // solver type
    out_file << "\n  -> RTE Solver: " 
             << IP.RTE_Solver;

    // DOM quadrature
    if (IP.i_RTE_Solver == RTE2D_SOLVER_DOM) {
      out_file << "\n      -> Quadrature: "
	       << IP.DOM_Quadrature;      
    // FVM angular discretization
    } else {
      out_file << "\n      -> Number of CA polar-direction: "
	       << IP.Number_of_Angles_Mdir;
      out_file << "\n      -> Number of CA azim-direction: " 
	       << IP.Number_of_Angles_Ldir;
    }

    // absorbsion model
    out_file << "\n  -> Absorption Model: " 
	     << IP.AbsorptionModel;
      
    // scattering function
    out_file << "\n  -> Scattering Function: " 
	     << IP.ScatteringFunc;

    // ICs
    out_file << "\n  -> Initial Conditions: " 
             << IP.ICs_Type;
    if (IP.i_ICs == IC_CONSTANT || IP.i_ICs == IC_UNIFORM ){
      out_file << "\n      -> Intensity (K): " 
	       << IP.Intensity;
      out_file << "\n      -> Temperature (K): " 
	       << IP.Temperature;
      out_file << "\n      -> Pressure (Pa): " 
	       << IP.Pressure;
      out_file << "\n      -> Mixture: " 
	       << "xco = " << IP.xco << ",  "
	       << "xh2o = " << IP.xh2o << ",  "
	       << "xco2 = " << IP.xco2 << ",  "
	       << "xo2 = " << IP.xo2 << ",  "
	       << "fsoot = " << IP.fsoot;
      out_file << "\n      -> Absorption Cefficient (m^-1): " 
	       << IP.AbsorptionCoef;
      out_file << "\n      -> Scattering Cefficient (m^-1): " 
	       << IP.ScatteringCoef;
    } /* endif */

    // Medium ICs
    out_file << "\n  -> Medium Initial Conditions: " 
             << IP.ICs_Medium;
    /***********************************************************************
     ***********************************************************************/

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
	out_file << "\n  -> Cylinder Radius (m): " 
		 << IP.Blunt_Body_Radius;
        break;
      case GRID_BUMP_CHANNEL_FLOW :
	if (strcmp(IP.Grid_Type,"Bump_Channel_Flow") == 0) {
	  if (IP.Smooth_Bump) out_file << " (smooth bump)";
	  else out_file << " (non-smooth bump)";
	}
	break;
      case GRID_NASA_ROTOR_37 :
        out_file << "\n  -> Percent Span: " 
                 << IP.Rotor_Percent_Span;
        break;
      case GRID_NASA_ROTOR_67 :
        out_file << "\n  -> Percent Span: " 
                 << IP.Rotor_Percent_Span;
        break;
	/***********************************************************************
	 *************************** RTE SPECIFIC ******************************/
      case GRID_RECTANGULAR_ENCLOSURE :
        out_file << "\n      -> Width of Solution Domain (m): " 
                 << IP.Box_Width;
        out_file << "\n      -> Height of Solution Domain (m): " 
                 << IP.Box_Height;

	// Wall conditions
	out_file << "\n      -> Wall Conditions: ";
	out_file << "\n          -> Wall Temperature (K) [N,S,E,W]: " 
		 << IP.NorthWallTemp << "  "
		 << IP.SouthWallTemp << "  "
		 << IP.EastWallTemp << "  "
		 << IP.WestWallTemp;
	out_file << "\n          -> Wall Emissivity [N,S,E,W]: " 
		 << IP.NorthWallEmiss << "  "
		 << IP.SouthWallEmiss << "  "
		 << IP.EastWallEmiss << "  "
		 << IP.WestWallEmiss;
        break;
      case GRID_CYLINDRICAL_ENCLOSURE :
        out_file << "\n      -> Pipe Length (m): " 
                 << IP.Pipe_Length;
        out_file << "\n      -> Pipe Radius (m): " 
                 << IP.Pipe_Radius;

	// Wall conditions
	out_file << "\n      -> Wall Conditions: ";
	out_file << "\n          -> Wall Temperature (K) [N,S,E,W]: " 
		 << IP.NorthWallTemp << "  "
		 << IP.SouthWallTemp << "  "
		 << IP.EastWallTemp << "  "
		 << IP.WestWallTemp;
	out_file << "\n          -> Wall Emissivity [N,S,E,W]: " 
		 << IP.NorthWallEmiss << "  "
		 << IP.SouthWallEmiss << "  "
		 << IP.EastWallEmiss << "  "
		 << IP.WestWallEmiss;
        break;
	/***********************************************************************
	 ***********************************************************************/
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

    out_file << endl << string(75,'*');
    out_file << endl << string(75,'*') << endl;


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
	     << IP.Number_of_Interface_Mesh_Refinements;
    if (IP.Number_of_Boundary_Mesh_Refinements > 0)
    out_file << "\n  -> Number of Interface Mesh Refinements: " 
	     << IP.Number_of_Interface_Mesh_Refinements;
    out_file << "\n  -> Number of Bounding-Box Mesh Refinements: " 
	     << IP.Number_of_Bounding_Box_Mesh_Refinements;
    out_file << "\n  -> Refinement Criteria: ";
    if (IP.Refinement_Criteria_Gradient_Density)
    out_file << "\n     -> Gradient of the density field refinement criteria";
    if (IP.Refinement_Criteria_Divergence_Velocity)
    out_file << "\n     -> Divergence of the velocity field refinement criteria";
    if (IP.Refinement_Criteria_Curl_Velocity)
    out_file << "\n     -> Curl of the velocity field refinement criteria";
    if (IP.Number_of_Bounding_Box_Mesh_Refinements)
      out_file << "\n     -> Bounding-box for bounding-box AMR:" << IP.AMR_Xmin << IP.AMR_Xmax;
    out_file << "\n  -> CFL Number: " 
             << IP.CFL_Number;
    if (IP.Time_Accurate) { 
      out_file << "\n  -> Maximum Time (ms): " << IP.Time_Max*THOUSAND;
    }
    out_file << "\n  -> Maximum Number of ";
    if (IP.Time_Accurate) { 
      out_file << "Time Steps: ";
    } else { 
      out_file << "Relaxation (Pseudo-Time) Iterations: ";
    }
    out_file << IP.Maximum_Number_of_Time_Steps;
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
 * Get CFD kit Path.                                         *
 *************************************************************/
inline void Rte2D_Input_Parameters::get_cffc_path(){
 
  // Check to see if environment varible exists.
  if(getenv(PATHVAR) == NULL){
    cerr<<"Missing Environment Variable: "<<PATHVAR<<"\n\n"; 
    exit(1);
  }
  
  //Set Path
  strcpy(CFFC_Path,getenv(PATHVAR));
}



/*************************************************************
 * Rte2D_Input_Parameters -- External subroutines.           *
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
