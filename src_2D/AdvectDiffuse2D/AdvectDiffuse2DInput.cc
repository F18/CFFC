/* AdvectDiffuse2DInput.cc:  Subroutines for 
                             2D Advection Diffusion Equation Input Classes. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "AdvectDiffuse2DInput.h" /* Include 2D advection diffusion equation input parameter header file. */
#include "AdvectDiffuse2DQuad.h"  /* Include AdvectDiffuse2D_Quad_Block header file. */
#include "../HighOrderReconstruction/CENO_ExecutionMode.h" // Include high-order CENO execution mode header file
#include "../HighOrderReconstruction/CENO_Tolerances.h"	   // Include high-order CENO tolerances header file
#include "../Grid/HO_Grid2DQuad_ExecutionMode.h" // Include high-order 2D grid execution mode header file
#include "../Grid/Tecplot_ExecutionMode.h" // Include Tecplot execution mode header file

/*************************************************************
 * AdvectDiffuse2D_Input_Parameters -- Member functions.     *
 *************************************************************/

/*!
 * \todo Set CENO ghost cells based on the HighOrder2D<int>::Nghost(ReconstructionOrder())
 */
int AdvectDiffuse2D_Input_Parameters::Nghost(void) const{

  int Number_Of_Ghost_Cells(0);
  
  switch(i_ReconstructionMethod){

  case RECONSTRUCTION_CENO:
    // This number if function of the specified reconstruction order
    // Use 'int' type to get the number of ghost cells (it doesn't matter which type is used!)
    Number_Of_Ghost_Cells = 5;
    break;

  default:
    Number_Of_Ghost_Cells = 2;
  }

  return Number_Of_Ghost_Cells;
}

/*!
 * Decide whether to output or not the boundary reference state for a particular boundary condition.
 * To get output for a particular BCtype, just add it to the list.
 */
bool AdvectDiffuse2D_Input_Parameters::OutputBoundaryReferenceState(const int & BCtype) const{
  if (BCtype == BC_DIRICHLET ||
      BCtype == BC_NEUMANN ||
      BCtype == BC_FARFIELD){
    // Get output
    return true;
  } else {
    // No output
    return false;
  }
}

/******************************************************//**
 * Parse the input file
 ********************************************************/
int AdvectDiffuse2D_Input_Parameters::Parse_Input_File(char *Input_File_Name_ptr){

  ostringstream msg;
  int command_flag;

  strcpy(Input_File_Name, Input_File_Name_ptr);
  Open_Input_File(*this);
  if (Input_File.fail()) {
    msg << "AdvectDiffuse2D_Input_Parameters::Parse_Input_File() ERROR: Unable to open "
	<<string(Input_File_Name_ptr) 
	<< " input data file.";
    if (Verbose()) {
      cerr << msg.str() << endl;
    }
    throw runtime_error(msg.str());
  } /* endif */

  if (Verbose()) {
    cout << "\n Reading input data file `"
	 << Input_File_Name << "'." << endl;
    cout.flush();
  }
  while (1) {
    Get_Next_Input_Control_Parameter();
    command_flag = Parse_Next_Input_Control_Parameter(*this);
    if (command_flag == EXECUTE_CODE) {
      break;
      
    } else if (command_flag == TERMINATE_CODE) {
      return (0);
      
    } else if (command_flag == INVALID_INPUT_CODE ||
	       command_flag == INVALID_INPUT_VALUE) {
      Line_Number = -Line_Number;
      
      msg << "AdvectDiffuse2D_Input_Parameters::Parse_Input_File() ERROR: Error reading data at line # " 
	  << -Line_Number
	  << " of input data file.";
      if (Verbose()){
	cerr << msg.str() << endl;
      }
      
      throw runtime_error(msg.str());
    } /* endif */
  } /* endwhile */

  // Perform update of the internal variables of the exact solution
  ExactSoln->Set_ParticularSolution_Parameters();

  // Perform update of the internal variables of the inflow field
  Inflow->Set_InflowField_Parameters();

  // Perform update of the internal variables of the high-order input parameters
  HighOrder2D_Input::Set_Final_Parameters(*this);

  // Set reference state in the AdvectDiffuse2D_Quad_Block class
  AdvectDiffuse2D_Quad_Block::Set_Normalization_Reference_State(RefU);

  // Set flag for including/excluding source term in the model equation
  AdvectDiffuse2D_Quad_Block::Include_Source_Term = Include_Source_Term;
}

/******************************************************//**
 * Get the next input control parameter from the input  
 * file.                                                
 ********************************************************/
void AdvectDiffuse2D_Input_Parameters::Get_Next_Input_Control_Parameter(void){

  int i, index, LineSize, IndexFirstChar(0);
  char buffer[256], ControlParameter[256];

  // Initialize ControlParameter and Next_Control_Parameter to end of string
  ControlParameter[0] = '\0';
  strcpy(Next_Control_Parameter, ControlParameter);

  // While the input stream is 'good' for reading and the end of file is not attained
  while ( Input_File.good() && !Input_File.getline(buffer, sizeof(buffer)).eof() ){

    // Process the line 
    Line_Number = Line_Number + 1;
    LineSize = Input_File.gcount(); // Get the size of the line. Last character is "\0"!

    // Determine the index of the first character different than 'space' and 'tab'
    for (i=0; i<LineSize; ++i){
      if (buffer[i] != ' ' && buffer[i] != '\t'){
	IndexFirstChar = i;
	break;
      }
    }

    /* Parse the line if the first character different than 'space' 
       is also different than '#' or end of string ('\0').
       Otherwise skip the line because it is either a comment or an empty line. */
    if ( buffer[IndexFirstChar] != '#' && buffer[IndexFirstChar] != '\0'){

      // Get the ControlParameter
      for(i=IndexFirstChar, index=0;  i<LineSize;  ++i, ++index){
	if (buffer[i] == ' ' || buffer[i] == '=' || buffer[i] == '\t'){
	  ControlParameter[index] = '\0';
	  break;
	} else {
	  ControlParameter[index] = buffer[i];
	}
      }

      // Set the Next_Control_Parameter
      strcpy(Next_Control_Parameter, ControlParameter);
      break;
    }

  }//endwhile
}

/***************************************************************
 * AdvectDiffuse2D_Input_Parameters -- Input-output operators. *
 ***************************************************************/
ostream &operator << (ostream &out_file,
		      const AdvectDiffuse2D_Input_Parameters &IP) {
    out_file << setprecision(6);
    out_file << "\n  -> CFFC Path: " 
	     << IP.CFFC_Path;
    out_file << "\n\n Solving 2D advection diffusion equations (IBVP/BVP) on multi-block solution-adaptive quadrilateral mesh.";
    out_file << "\n  -> Input File Name: " 
             << IP.Input_File_Name;

    // ==== Initial condition parameters ====
    out_file << "\n  -> Initial Conditions: " 
             << IP.ICs_Type;
    if (IP.i_ICs == IC_EXACT_SOLUTION || IP.i_ICs == IC_INTERIOR_UNIFORM_GHOSTCELLS_EXACT){
      out_file << "\n     -> Exact integration digits = " << IP.Exact_Integration_Digits;
    }

    // ====    Velocity field parameters ====
    VelocityFields::Print_Info(out_file);

    // ====    Diffusion field parameters ====
    DiffusionFields::Print_Info(out_file);

    // ====    Source field parameters ====
    if (IP.Include_Source_Term) {
      out_file << "\n  -> Source Term Included in Equation: Yes ";
    } else {
      out_file << "\n  -> Source Term Included in Equation: No ";
    }
    IP.SourceTerm->Print_Info(out_file);

    // ====    Boundary conditions ====
    if (IP.BCs_Specified) {
      out_file << "\n  -> Boundary conditions specified as: ";

      // North
      out_file << "\n     -> BC_North = " << IP.BC_North_Type;
      if (IP.OutputBoundaryReferenceState(IP.BC_North)){
	out_file << " , Ref_State = " << IP.Ref_State_BC_North;
      }
      // South
      out_file << "\n     -> BC_South = " << IP.BC_South_Type;
      if (IP.OutputBoundaryReferenceState(IP.BC_South)){
	out_file << " , Ref_State = " << IP.Ref_State_BC_South;
      }
      // East
      out_file << "\n     -> BC_East  = " << IP.BC_East_Type;
      if (IP.OutputBoundaryReferenceState(IP.BC_East)){
	out_file << " , Ref_State = " << IP.Ref_State_BC_East;
      }
      // West
      out_file << "\n     -> BC_West  = " << IP.BC_West_Type;
      if (IP.OutputBoundaryReferenceState(IP.BC_West)){
	out_file << " , Ref_State = " << IP.Ref_State_BC_West;
      }
    }

    // ====    Inflow field ====
    if (IP.Inflow->IsInflowFieldSet()){
      IP.Inflow->Print_Info(out_file);
    }

    // ====    Exact solution parameters ====
    IP.ExactSoln->Print_Info(out_file);

    // ==== Time marching scheme parameters ====
    if (IP.Time_Accurate) { 
       out_file << "\n  -> Time Accurate (Unsteady) Solution";
    } else {
       out_file << "\n  -> Time Invariant (Steady-State) Solution";
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
    out_file << "\n  -> CFL Number: " 
             << IP.CFL_Number;
    out_file << "\n  -> Maximum Time: " 
             << IP.Time_Max;
    out_file << "\n  -> Maximum Number of Time Steps (Iterations): " 
             << IP.Maximum_Number_of_Time_Steps;
    if (IP.NKS_IP.Maximum_Number_of_NKS_Iterations > 0) {
      out_file << "\n  -> Maximum Number of NKS Iterations: " 
	       << IP.NKS_IP.Maximum_Number_of_NKS_Iterations;
    }

    // ==== Spatial approximation parameters ====
    if (IP.Axisymmetric) { 
       out_file << "\n  -> 2D Axisymmetric Geometry";
    } else {
       out_file << "\n  -> 2D Planar Geometry";
    }
    if (IP.i_Reconstruction == RECONSTRUCTION_HIGH_ORDER){
      out_file << "\n  -> Space Accuracy : ";
      switch(IP.Space_Accuracy){
      case 1: 
	out_file << "1st-order";
	break;
      case 2:
	out_file << "2nd-order";
	break;		
      case 3:		
	out_file << "3rd-order";
	break;		
      case 4:		
	out_file << "4th-order";
	break;		
      case 5:		
	out_file << "5th-order";
	break;		
      case 6:		
	out_file << "6th-order";
	break;
      default:
	out_file << "bigger than 6th-order";
      }
    } else {
      out_file << "\n  -> Space Accuracy : 2nd-order";
    }
    out_file << "\n  -> Reconstruction: " 
             << IP.Reconstruction_Type;
    if (IP.i_ReconstructionMethod == RECONSTRUCTION_CENO){
      CENO_Execution_Mode::Print_Info(out_file);
      CENO_Tolerances::Print_Info(out_file);
      out_file << "\n     -> Reference State: "
	       << IP.RefU;
    }

    // output information related to auxiliary reconstructions.
    HighOrder2D_Input::Print_Info(out_file);

    out_file << "\n  -> Limiter: " 
             << IP.Limiter_Type;
    if (IP.Limiter_Type != LIMITER_ZERO && IP.Freeze_Limiter) {
      out_file << "\n  -> Freeze Limiter when L2-norm of residual is < "
	       << IP.Freeze_Limiter_Residual_Level;
    } /* endif */
    if (IP.i_ReconstructionMethod != RECONSTRUCTION_CENO ){
      out_file << "\n  -> Elliptic Flux Evaluation: " << IP.Viscous_Reconstruction_Type;
    }
    // output information related to the treatment of curved boundaries.
    HO_Grid2D_Execution_Mode::Print_Info(out_file);

    // ==== Grid parameters ====
    out_file << "\n  -> Grid: " 
             << IP.Grid_Type;
    switch(IP.i_Grid) {
      case GRID_SQUARE :
        out_file << "\n     -> Size of Solution Domain: " 
                 << IP.Box_Width;
        break;
      case GRID_RECTANGULAR_BOX :
        out_file << "\n     -> Width of Solution Domain: " 
                 << IP.Box_Width;
        out_file << "\n     -> Height of Solution Domain: " 
                 << IP.Box_Height;
        break;
      case GRID_DEFORMED_BOX :
        out_file << "\n     -> SW Corner: " 
                 << IP.VertexSW;
        out_file << "\n     -> SE Corner: " 
                 << IP.VertexSE;
        out_file << "\n     -> NE Corner: " 
                 << IP.VertexNE;
        out_file << "\n     -> NW Corner: " 
                 << IP.VertexNW;
        break;
      case GRID_FLAT_PLATE :
        out_file << "\n     -> Plate Length: " 
                 << IP.Plate_Length;
        break;
      case GRID_PIPE :
        out_file << "\n     -> Pipe Length: " 
                 << IP.Pipe_Length;
        out_file << "\n     -> Pipe Radius: " 
                 << IP.Pipe_Radius;
        break;
      case GRID_BLUNT_BODY :
        out_file << "\n     -> Cylinder Radius: " 
                 << IP.Blunt_Body_Radius;
        break;
      case GRID_ROCKET_MOTOR :
	out_file << "\n     -> Length of Chamber (m): "
		 << IP.Chamber_Length;
	out_file << "\n     -> Radius of Chamber (m): "
		 << IP.Chamber_Radius;
	out_file << "\n     -> Distance from Chamber to Nozzle Throat (m): "
		 << IP.Chamber_To_Throat_Length;
	out_file << "\n     -> Length of the Nozzle (m): "
		 << IP.Nozzle_Length;
	out_file << "\n     -> Radius of the Nozzle at Throat (m): "
		 << IP.Nozzle_Radius_Throat;
	out_file << "\n     -> Radius of the Nozzle at Exit(m): "
		 << IP.Nozzle_Radius_Exit;
	out_file << "\n     -> Radius of the Propellant Grain (m): "
		 << IP.Grain_Radius;
	out_file << "\n     -> Nozzle type: "
		 << IP.Nozzle_Type;
        break;
      case GRID_CIRCULAR_CYLINDER :
        out_file << "\n     -> Inner Cylinder Radius (m): " 
                 << IP.Cylinder_Radius
		 << "\n     -> Outer Cylinder Radius (m): " 
                 << IP.Cylinder_Radius2;
        break;
      case GRID_ANNULUS :
        out_file << "\n     -> Inner Cylinder Radius (m): " 
                 << IP.Cylinder_Radius
		 << "\n     -> Outer Cylinder Radius (m): " 
                 << IP.Cylinder_Radius2
		 << "\n     -> Start Theta (degrees): " 
		 << IP.Annulus_Theta_Start
		 << "\n     -> End Theta (degrees): " 
		 << IP.Annulus_Theta_End;
        break;
      case GRID_ELLIPSE :
        out_file << "\n     -> Width of Ellipse along x-axis: " 
                 << IP.Ellipse_Length_X_Axis;
        out_file << "\n     -> Height of Ellipse along y-axis: " 
                 << IP.Ellipse_Length_Y_Axis;
        break;
      case GRID_NACA_AEROFOIL :
        out_file << "\n     -> NACA " 
                 << IP.NACA_Aerofoil_Type;
        out_file << "\n     -> Chord Length: " 
                 << IP.Chord_Length;
        break;
      case GRID_FREE_JET :
        out_file << "\n     -> Orifice Radius: " 
                 << IP.Orifice_Radius;
        break;
      case GRID_ICEMCFD :
        break;
      case GRID_READ_FROM_DEFINITION_FILE :
        break;
      case GRID_READ_FROM_GRID_DATA_FILE :
        break;
      default:
        out_file << "\n     -> Width of Solution Domain: " 
                 << IP.Box_Width;
        out_file << "\n     -> Height of Solution Domain: " 
                 << IP.Box_Height;
        break;
    } /* endswitch */
    out_file << "\n     -> Smooth Quad Block: ";
    if (IP.i_Smooth_Quad_Block){
      out_file << "Yes";
    } else {
      out_file << "No";
    }
    if (IP.IterationsOfInteriorNodesDisturbances > 0){
      out_file << "\n     -> Disturbed Interior Quad Block Nodes: "
	       << IP.IterationsOfInteriorNodesDisturbances << " iterations.";
    }
    out_file << "\n     -> Mesh shift, scale, and rotate: " 
             << IP.X_Shift << " " << IP.X_Scale << " " << IP.X_Rotate;
    out_file << "\n     -> Number of Blocks i-direction: "
             << IP.Number_of_Blocks_Idir;
    out_file << "\n     -> Number of Blocks j-direction: " 
             << IP.Number_of_Blocks_Jdir;
    out_file << "\n     -> Number of Cells i-direction: "
             << IP.Number_of_Cells_Idir;
    out_file << "\n     -> Number of Cells j-direction: " 
             << IP.Number_of_Cells_Jdir;
    out_file << "\n     -> Number of Ghost Cells: "
	     << IP.Number_of_Ghost_Cells;

    // ==== AMR parameters ====
    if (IP.Number_of_Initial_Mesh_Refinements > 0)
    out_file << "\n  -> Number of Initial Mesh Refinements : " 
             << IP.Number_of_Initial_Mesh_Refinements;
    if (IP.Number_of_Uniform_Mesh_Refinements > 0)
    out_file << "\n  -> Number of Uniform Mesh Refinements : " 
	     << IP.Number_of_Uniform_Mesh_Refinements;
    if (IP.Number_of_Boundary_Mesh_Refinements > 0)
    out_file << "\n  -> Number of Boundary Mesh Refinements : " 
	     << IP.Number_of_Boundary_Mesh_Refinements;
    out_file << "\n  -> Number of Processors: " 
             << IP.Number_of_Processors;
    out_file << "\n  -> Number of Blocks Per Processor: " 
             << IP.Number_of_Blocks_Per_Processor;
    if (IP.AMR) {
      out_file << "\n  -> AMR Frequency: "
	       << IP.AMR_Frequency
	       << " steps (iterations)";
      if(IP.AMR_Global_Reference == ON){
	out_file << "\n  -> AMR Reference: Absolute Simulation Start";
      } else {
	out_file << "\n  -> AMR Reference: Current Simulation Start";
      }
      
      out_file << "\n  -> Threshold for Refinement: "
	       << IP.Threshold_for_Refinement
	       << "\n  -> Threshold for Coarsening: "
	       << IP.Threshold_for_Coarsening
	       << "\n  -> Maximum Refinement Level: "
	       << IP.Maximum_Refinement_Level
	       << "\n  -> Minimum Refinement Level: "
	       << IP.Minimum_Refinement_Level;
    }

    // ==== Output parameters ====
    out_file << "\n  -> Output File Name: " 
             << IP.Output_File_Name;
    out_file << "\n  -> Output Format: " 
             << IP.Output_Format_Type;
    if (IP.i_Output_Format == IO_TECPLOT){
      // output information related to Tecplot output
      Tecplot_Execution_Mode::Print_Info(out_file);
    }
    out_file << "\n  -> Restart Solution Save Frequency: "
             << IP.Restart_Solution_Save_Frequency
             << " steps (iterations)"; 
    out_file.flush();
    return (out_file);
}

istream &operator >> (istream &in_file,
		      AdvectDiffuse2D_Input_Parameters &IP) {
  return (in_file);
}

/*************************************************************
 * AdvectDiffuse2D_Input_Parameters -- External subroutines. *
 *************************************************************/

/******************************************************//**
 * Routine: Open_Input_File                             
 *                                                      
 * Opens the appropriate input data file.               
 *                                                      
 ********************************************************/
void Open_Input_File(AdvectDiffuse2D_Input_Parameters &IP) {

    IP.Input_File.open(IP.Input_File_Name, ios::in);
    if (!IP.Input_File.fail()) {
       IP.Line_Number = 0;
       IP.Input_File.setf(ios::skipws);
    } /* endif */

}

/******************************************************//**
 * Routine: Close_Input_File                            
 *                                                      
 * Closes the appropriate input data file.              
 *                                                      
 ********************************************************/
void Close_Input_File(AdvectDiffuse2D_Input_Parameters &IP) {

    IP.Input_File.unsetf(ios::skipws);
    IP.Input_File.close();

}

/******************************************************//**
 * Routine: Set_Default_Input_Parameters                
 *                                                      
 * Assigns default values to the input parameters.      
 *                                                      
 ********************************************************/
void Set_Default_Input_Parameters(AdvectDiffuse2D_Input_Parameters &IP) {

    int i;
    char *string_ptr;

    // CFFC root directory path:
    IP.get_cffc_path();

    string_ptr = "AdvectDiffuse2D.in";
    strcpy(IP.Input_File_Name, string_ptr);

    // Time-stepping parameters:
    string_ptr = "Explicit_Euler";
    strcpy(IP.Time_Integration_Type, string_ptr);
    IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
    IP.Time_Accurate = 0;
    IP.CFL_Number = 0.5;
    IP.Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
    IP.Maximum_Number_of_Time_Steps = 100;
    IP.N_Stage = 1;
    IP.Time_Max = ZERO;

    // Residual variable:
    IP.i_Residual_Variable = 1;

    // Residual smoothing:
    IP.Residual_Smoothing = 0;
    IP.Residual_Smoothing_Epsilon = ZERO;
    IP.Residual_Smoothing_Gauss_Seidel_Iterations = 2;

    // Reconstruction type:
    string_ptr = "Least_Squares";
    strcpy(IP.Reconstruction_Type, string_ptr);
    IP.i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;
    IP.Space_Accuracy = 1;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_ReconstructionMethod = RECONSTRUCTION_LEAST_SQUARES;
    CENO_Execution_Mode::USE_CENO_ALGORITHM = OFF;

    // Viscous gradient reconstruction type:
    string_ptr = "Diamond_Path";
    strcpy(IP.Viscous_Reconstruction_Type,string_ptr);
    IP.i_Viscous_Reconstruction = VISCOUS_RECONSTRUCTION_DIAMOND_PATH;

    // Limiter type:
    string_ptr = "Barth_Jespersen";
    strcpy(IP.Limiter_Type, string_ptr);
    IP.i_Limiter = LIMITER_BARTH_JESPERSEN;
    IP.Freeze_Limiter = 0;
    IP.Freeze_Limiter_Residual_Level = 1e-4;

    // Initial conditions:
    string_ptr = "Uniform";
    strcpy(IP.ICs_Type, string_ptr);
    IP.i_ICs = IC_UNIFORM;
    IP.Exact_Integration_Digits = 9;

    // Source term:
    IP.SourceTerm->SetSourceField(SOURCE_FIELD_ZERO); // set the default source term (no source field)

    // State conditions:
    IP.Uo = AdvectDiffuse2D_State(ONE);
    IP.U1 = IP.Uo; IP.U1.u = ZERO;
    IP.U2 = IP.Uo; IP.U2.u = -ONE;
    IP.RefU = IP.Uo;

    // Geometry switch:
    string_ptr = "Planar";
    strcpy(IP.Flow_Geometry_Type, string_ptr);
    IP.Axisymmetric = 0;
    IP.Include_Source_Term = ON;

    // Grid parameters:
    string_ptr = "Square";
    strcpy(IP.Grid_Type, string_ptr);
    IP.i_Grid = GRID_SQUARE;
    IP.Box_Width = ONE;
    IP.Box_Height = ONE;
    IP.Number_of_Cells_Idir = 100;
    IP.Number_of_Cells_Jdir = 100;
    IP.Number_of_Ghost_Cells = 2;
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Plate_Length = ONE;
    IP.Pipe_Length = ONE;
    IP.Pipe_Radius = HALF;
    IP.Blunt_Body_Radius = ONE;
    IP.Blunt_Body_Mach_Number = TWO;
    IP.Chamber_Length = 0.835;
    IP.Chamber_Radius = 0.020;
    IP.Chamber_To_Throat_Length = 0.05;
    IP.Nozzle_Length = 0.150;
    IP.Nozzle_Radius_Exit = 0.030;
    IP.Nozzle_Radius_Throat = 0.010;
    IP.Nozzle_Type = NOZZLE_GOTTLIEB_FUNCTION;
    IP.Grain_Radius = 0.0;
    IP.Cylinder_Radius = ONE;
    IP.Cylinder_Radius2 = 32.00;
    IP.Annulus_Theta_Start = 0.0;
    IP.Annulus_Theta_End = 90.0;
    IP.Ellipse_Length_X_Axis = TWO;
    IP.Ellipse_Length_Y_Axis = HALF;
    IP.Chord_Length = ONE;
    IP.Orifice_Radius = ONE;
    IP.VertexSW = Vector2D(-0.5,-0.5);
    IP.VertexSE = Vector2D( 0.5,-0.5);
    IP.VertexNE = Vector2D( 0.5, 0.5);
    IP.VertexNW = Vector2D(-0.5, 0.5);
    IP.X_Shift = Vector2D_ZERO;
    IP.X_Scale = ONE;
    IP.X_Rotate = ZERO;
    IP.i_Smooth_Quad_Block = ON;        // Smooth quad block flag:
    IP.IterationsOfInteriorNodesDisturbances = 0;     /* Number of iterations of disturbing the mesh 
							 (create an unsmooth interior mesh) */
    IP.Num_Of_Spline_Control_Points = 361; /* Number of control points on the 2D spline (used for some grids) */

    // ICEM:
    IP.ICEMCFD_FileNames = ICEMCFD_get_filenames();

    // Mesh stretching factor:
    IP.i_Mesh_Stretching = OFF;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_LINEAR;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_LINEAR;
    IP.Mesh_Stretching_Factor_Idir = 1.01;
    IP.Mesh_Stretching_Factor_Jdir = 1.01;

    // Boundary conditions:
    string_ptr = "OFF";
    strcpy(IP.Boundary_Conditions_Specified,string_ptr);
    IP.BCs_Specified = OFF;
    string_ptr = "None";
    strcpy(IP.BC_North_Type,string_ptr);
    strcpy(IP.BC_South_Type,string_ptr);
    strcpy(IP.BC_East_Type,string_ptr);
    strcpy(IP.BC_West_Type,string_ptr);
    IP.BC_North = BC_NONE;
    IP.BC_South = BC_NONE;
    IP.BC_East  = BC_NONE;
    IP.BC_West  = BC_NONE;
    IP.Ref_State_BC_North = 0.0;
    IP.Ref_State_BC_South = 0.0;
    IP.Ref_State_BC_East = 0.0;
    IP.Ref_State_BC_West = 0.0;

    // AMR:
    IP.AMR = 0;
    IP.AMR_Frequency = 100;
    IP.Number_of_Initial_Mesh_Refinements = 0;
    IP.Number_of_Uniform_Mesh_Refinements = 0;
    IP.Number_of_Boundary_Mesh_Refinements = 0;
    IP.Maximum_Refinement_Level = 100;
    IP.Minimum_Refinement_Level = 1;
    IP.Threshold_for_Refinement = 0.50;
    IP.Threshold_for_Coarsening = 0.10;

    // Default output file names and parameters:
    string_ptr = "outputfile.dat";
    strcpy(IP.Output_File_Name, string_ptr);
    string_ptr = "gridfile.grid";
    strcpy(IP.Grid_File_Name, string_ptr);
    string_ptr = "gridfile.griddef";
    strcpy(IP.Grid_Definition_File_Name, string_ptr);
    string_ptr = "restartfile.soln";
    strcpy(IP.Restart_File_Name, string_ptr);
    string_ptr = "gnuplotfile.gplt";
    strcpy(IP.Gnuplot_File_Name, string_ptr);
    string_ptr = "Tecplot";
    strcpy(IP.Output_Format_Type, string_ptr);
    IP.i_Output_Format = IO_TECPLOT;
    IP.Restart_Solution_Save_Frequency = 1000;

    // Default output progress frequency:
    IP.Output_Progress_Frequency = 50;

    // Input_file parameters:
    string_ptr = " ";
    strcpy(IP.Next_Control_Parameter, string_ptr);
    IP.Line_Number = 0;
    IP.Number_of_Processors = CFFC_MPI::Number_of_Processors;
    IP.Number_of_Blocks_Per_Processor = 10;

    // Accuracy assessment parameters:
    IP.Accuracy_Assessment_Mode = ACCURACY_ASSESSMENT_BASED_ON_EXACT_SOLUTION;
    IP.Accuracy_Assessment_Exact_Digits = 10;
    IP.Accuracy_Assessment_Parameter = 1;

    // High-order parameters:
    HighOrder2D_Input::SetDefaults();
}

/******************************************************//**
 * Routine: Broadcast_Input_Parameters                  
 *                                                      
 * Broadcast the input parameters variables to all      
 * processors involved in the calculation from the      
 * primary processor using the MPI broadcast routine.   
 *                                                      
 ********************************************************/
void Broadcast_Input_Parameters(AdvectDiffuse2D_Input_Parameters &IP) {

#ifdef _MPI_VERSION
    int i;

    MPI::COMM_WORLD.Bcast(IP.CFFC_Path, 
 			  INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
			  MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Input_File_Name, 
                          INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Line_Number), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Time_Integration_Type, 
                          INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Time_Integration), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Time_Accurate), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Local_Time_Stepping), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Maximum_Number_of_Time_Steps), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.N_Stage), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.CFL_Number), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Time_Max), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Residual_Smoothing), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Residual_Smoothing_Epsilon), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Residual_Smoothing_Gauss_Seidel_Iterations), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Reconstruction_Type, 
                          INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Reconstruction), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_ReconstructionMethod), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Space_Accuracy), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.IncludeHighOrderBoundariesRepresentation), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Viscous_Reconstruction_Type,
			  INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Viscous_Reconstruction),
			  1,
			  MPI::INT,0);

    MPI::COMM_WORLD.Bcast(IP.Limiter_Type, 
                          INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Limiter), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.ICs_Type, 
                          INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_ICs), 
                          1, 
                          MPI::INT, 0);
    if (!CFFC_Primary_MPI_Processor()) {
      IP.Uo = AdvectDiffuse2D_State(ONE);
      IP.U1 = IP.Uo; IP.U1.u = ZERO;
      IP.U2 = IP.Uo; IP.U2.u = -ONE;
    } /* endif */

    // Reference state
    MPI::COMM_WORLD.Bcast(&(IP.RefU.u), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Exact_Integration_Digits), 
                          1, 
                          MPI::DOUBLE, 0);

    MPI::COMM_WORLD.Bcast(IP.Flow_Geometry_Type, 
                          INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Axisymmetric), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Include_Source_Term), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Grid_Type, 
                          INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.NACA_Aerofoil_Type, 
                          INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Grid), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Cells_Idir), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Cells_Jdir), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Ghost_Cells), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Idir), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Jdir), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Box_Width), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Box_Height), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Plate_Length), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Pipe_Length), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Pipe_Radius), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Blunt_Body_Radius), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Blunt_Body_Mach_Number), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Chamber_Length),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Chamber_Radius),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Chamber_To_Throat_Length),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Nozzle_Length),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Nozzle_Radius_Exit),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Nozzle_Radius_Throat),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Nozzle_Type),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Grain_Radius),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Cylinder_Radius), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Cylinder_Radius2), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Annulus_Theta_Start), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Annulus_Theta_End), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Ellipse_Length_X_Axis), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Ellipse_Length_Y_Axis), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Chord_Length), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Orifice_Radius), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.VertexSW.x), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.VertexSW.y), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.VertexSE.x), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.VertexSE.y), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.VertexNW.x), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.VertexNW.y), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.VertexNE.x), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.VertexNE.y), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.X_Shift.x), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.X_Shift.y), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.X_Scale), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.X_Rotate), 
                          1, 
                          MPI::DOUBLE, 0);
    // Boundary Conditions:
    MPI::COMM_WORLD.Bcast(IP.Boundary_Conditions_Specified,
			  INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(&(IP.BCs_Specified),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(IP.BC_North_Type,
			  INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(IP.BC_South_Type,
			  INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(IP.BC_East_Type,
			  INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(IP.BC_West_Type,
			  INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(&(IP.BC_North),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.BC_South),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.BC_East),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.BC_West),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Ref_State_BC_North.u),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Ref_State_BC_South.u),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Ref_State_BC_East.u),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Ref_State_BC_West.u),
			  1,
			  MPI::DOUBLE,0);

    // ICEM:
    if (!CFFC_Primary_MPI_Processor()) {
       IP.ICEMCFD_FileNames = new char*[3];
       for (i = 0; i < 3; i++) {
          IP.ICEMCFD_FileNames[i] = new char[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
       } /* endfor */
    } /* endif */
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[0], 
                          INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[1], 
                          INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[2], 
                          INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                          MPI::CHAR, 0);

    // AMR & Refinement Parameters
    MPI::COMM_WORLD.Bcast(&(IP.AMR), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.AMR_Frequency),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.AMR_Global_Reference),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Initial_Mesh_Refinements), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Uniform_Mesh_Refinements),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Boundary_Mesh_Refinements),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Maximum_Refinement_Level),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Minimum_Refinement_Level),
                          1,
                          MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Threshold_for_Refinement), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Threshold_for_Coarsening), 
                          1, 
                          MPI::DOUBLE, 0);

    // Mesh stretching flag:
    MPI::COMM_WORLD.Bcast(&(IP.i_Mesh_Stretching),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Mesh_Stretching_Type_Idir),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Mesh_Stretching_Type_Jdir),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Mesh_Stretching_Factor_Idir),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Mesh_Stretching_Factor_Jdir),
			  1,
			  MPI::DOUBLE,0);
    // Smooth quad block flag:
    MPI::COMM_WORLD.Bcast(&(IP.i_Smooth_Quad_Block),
			  1,
			  MPI::INT,0);

    MPI::COMM_WORLD.Bcast(&(IP.IterationsOfInteriorNodesDisturbances),
			  1,
			  MPI::INT,0);

    // File Names
    MPI::COMM_WORLD.Bcast(IP.Output_File_Name, 
                          INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Grid_File_Name, 
                          INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Grid_Definition_File_Name, 
                          INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Restart_File_Name, 
                          INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Gnuplot_File_Name, 
                          INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Output_Format_Type, 
                          INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Output_Format), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Restart_Solution_Save_Frequency), 
                          1, 
                          MPI::INT, 0);

    // Output progress frequency:
    MPI::COMM_WORLD.Bcast(&(IP.Output_Progress_Frequency),
			  1,
			  MPI::INT,0);

    // Multigrid Related Parameters
    IP.Multigrid_IP.Broadcast_Input_Parameters();

    // NKS Parameters
    IP.NKS_IP.Broadcast_Input_Parameters();

    if (!CFFC_Primary_MPI_Processor()) {
       IP.Number_of_Processors = CFFC_MPI::Number_of_Processors;
    } /* endif */
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Per_Processor), 
                          1, 
                          MPI::INT, 0);

    // Freeze_Limiter
    MPI::COMM_WORLD.Bcast(&(IP.Freeze_Limiter), 
			  1, 
			  MPI::INT, 0);
    // Freeze_Limiter_Residual_Level
    MPI::COMM_WORLD.Bcast(&(IP.Freeze_Limiter_Residual_Level),
                          1, 
                          MPI::DOUBLE, 0);

    // Accuracy assessment parameters:
    MPI::COMM_WORLD.Bcast(&(IP.Accuracy_Assessment_Mode), 
			  1, 
			  MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Accuracy_Assessment_Exact_Digits), 
			  1, 
			  MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Accuracy_Assessment_Parameter), 
			  1, 
			  MPI::INT, 0);

    // CENO_Execution_Mode variables
    CENO_Execution_Mode::Broadcast();
    
    // CENO_Tolerances variables
    CENO_Tolerances::Broadcast();

    // VelocityFields variables
    VelocityFields::Broadcast();

    // DiffusionFields variables
    DiffusionFields::Broadcast();

    // SourceTermFields variables
    IP.SourceTerm->Broadcast();

    // Exact solution variables
    IP.ExactSoln->Broadcast();

    // Inflow field variables
    IP.Inflow->Broadcast();

    // HO_Grid2D_Execution_Mode variables
    HO_Grid2D_Execution_Mode::Broadcast();

    // Tecplot_Execution_Mode variables
    Tecplot_Execution_Mode::Broadcast();

    // HighOrder2D_Input variables
    HighOrder2D_Input::Broadcast();    

    // Update all dependent variables
    if (!CFFC_Primary_MPI_Processor()) {
      // Perform update of the internal variables of the exact solution
      IP.ExactSoln->Set_ParticularSolution_Parameters();

      // Perform update of the internal variables of the inflow field
      IP.Inflow->Set_InflowField_Parameters();

      // Perform update of the internal variables of the high-order input parameters
      HighOrder2D_Input::Set_Final_Parameters(IP);

      // Set reference state in the AdvectDiffuse2D_Quad_Block class
      AdvectDiffuse2D_Quad_Block::Set_Normalization_Reference_State(IP.RefU);

      // Set flag for including/excluding source term in the model equation
      AdvectDiffuse2D_Quad_Block::Include_Source_Term = IP.Include_Source_Term;
    }

#endif

}

#ifdef _MPI_VERSION
/******************************************************//**
 * Routine: Broadcast_Input_Parameters                  
 *                                                      
 * Broadcast the input parameters variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.                                              
 *                                                      
 ********************************************************/
void Broadcast_Input_Parameters(AdvectDiffuse2D_Input_Parameters &IP,
                                MPI::Intracomm &Communicator, 
                                const int Source_CPU) {

    int Source_Rank = 0;
    int i;

    Communicator.Bcast(IP.CFFC_Path, 
 		       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
		       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Input_File_Name, 
                       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.Line_Number), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Time_Integration_Type, 
                       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Time_Integration), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Time_Accurate), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Local_Time_Stepping), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Maximum_Number_of_Time_Steps), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.N_Stage), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.CFL_Number), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Time_Max), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Residual_Smoothing), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Residual_Smoothing_Epsilon), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Residual_Smoothing_Gauss_Seidel_Iterations), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Reconstruction_Type, 
                       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Reconstruction), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.i_ReconstructionMethod), 
		       1, 
		       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Space_Accuracy), 
		       1, 
		       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.IncludeHighOrderBoundariesRepresentation), 
		       1, 
		       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Viscous_Reconstruction_Type,
		       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(&(IP.i_Viscous_Reconstruction),
		       1,
		       MPI::INT,Source_Rank);

    Communicator.Bcast(IP.Limiter_Type, 
                       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Limiter), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.ICs_Type, 
                       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_ICs), 
                       1, 
                       MPI::INT, Source_Rank);
    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
      IP.Uo = AdvectDiffuse2D_State(ONE);
      IP.U1 = IP.Uo; IP.U1.u = ZERO;
      IP.U2 = IP.Uo; IP.U2.u = -ONE;
    } /* endif */

    // Reference state
    Communicator.Bcast(&(IP.RefU.u), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Exact_Integration_Digits), 
		       1, 
		       MPI::DOUBLE, Source_Rank);

    Communicator.Bcast(IP.Flow_Geometry_Type, 
                       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.Axisymmetric), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Include_Source_Term), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Grid_Type, 
                       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.NACA_Aerofoil_Type, 
                       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Grid), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Cells_Idir), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Cells_Jdir), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Ghost_Cells), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Blocks_Idir), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Blocks_Jdir), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Box_Width), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Box_Height), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Plate_Length), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Pipe_Length), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Pipe_Radius), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Blunt_Body_Radius), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Blunt_Body_Mach_Number), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Chamber_Length),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Chamber_Radius),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Chamber_To_Throat_Length),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Nozzle_Length),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Nozzle_Radius_Exit),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Nozzle_Radius_Throat),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Nozzle_Type),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Grain_Radius),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Cylinder_Radius), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Cylinder_Radius2), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Annulus_Theta_Start), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Annulus_Theta_End), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Ellipse_Length_X_Axis), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Ellipse_Length_Y_Axis), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Chord_Length), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Orifice_Radius), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.VertexSW.x), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.VertexSW.y), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.VertexSE.x), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.VertexSE.y), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.VertexNE.x), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.VertexNE.y), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.VertexNW.x), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.VertexNW.y), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.X_Shift.x), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.X_Shift.y), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.X_Scale), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.X_Rotate), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    // Boundary Conditions:
    Communicator.Bcast(IP.Boundary_Conditions_Specified,
		       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(&(IP.BCs_Specified),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(IP.BC_North_Type,
		       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(IP.BC_South_Type,
		       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(IP.BC_East_Type,
		       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(IP.BC_West_Type,
		       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(&(IP.BC_North),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.BC_South),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.BC_East),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.BC_West),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Ref_State_BC_North.u),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Ref_State_BC_South.u),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Ref_State_BC_East.u),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Ref_State_BC_West.u),
		       1,
		       MPI::DOUBLE,Source_Rank);

    // ICEM :
    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
       IP.ICEMCFD_FileNames = new char*[3];
       for (i = 0; i < 3; i++) {
          IP.ICEMCFD_FileNames[i] = new char[INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D];
       } /* endfor */
    } /* endif */
    Communicator.Bcast(IP.ICEMCFD_FileNames[0], 
                       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.ICEMCFD_FileNames[1], 
                       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.ICEMCFD_FileNames[2], 
                       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                       MPI::CHAR, Source_Rank);
    // AMR & Refinement Parameters
    Communicator.Bcast(&(IP.AMR), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.AMR_Frequency),
                       1,
                       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.AMR_Global_Reference),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Initial_Mesh_Refinements), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Uniform_Mesh_Refinements),
                       1,
                       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Boundary_Mesh_Refinements),
                       1,
                       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Maximum_Refinement_Level),
                       1,
                       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Minimum_Refinement_Level),
                       1,
                       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Threshold_for_Refinement), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Threshold_for_Coarsening), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    // Mesh stretching flag:
    Communicator.Bcast(&(IP.i_Mesh_Stretching),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Mesh_Stretching_Type_Idir),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Mesh_Stretching_Type_Jdir),		     
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Mesh_Stretching_Factor_Idir),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Mesh_Stretching_Factor_Jdir),
		       1,
		       MPI::DOUBLE,Source_Rank);
    // Smooth quad block flag:
    Communicator.Bcast(&(IP.i_Smooth_Quad_Block),
		       1,
		       MPI::INT,Source_Rank);

    Communicator.Bcast(&(IP.IterationsOfInteriorNodesDisturbances),
		       1,
		       MPI::INT,Source_Rank);

    // File Names
    Communicator.Bcast(IP.Output_File_Name, 
                       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Grid_File_Name, 
                       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Grid_Definition_File_Name, 
                       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Restart_File_Name, 
                       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Gnuplot_File_Name, 
                       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Output_Format_Type, 
                       INPUT_PARAMETER_LENGTH_ADVECTDIFFUSE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Output_Format), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Restart_Solution_Save_Frequency), 
                       1, 
                       MPI::INT, Source_Rank);

    // Output progress frequency:
    Communicator.Bcast(&(IP.Output_Progress_Frequency),
		       1,
		       MPI::INT,Source_Rank);

    // Multigrid Related Parameters
    IP.Multigrid_IP.Broadcast_Input_Parameters(Communicator,
					       Source_CPU);

    // NKS Parameters
    IP.NKS_IP.Broadcast_Input_Parameters(Communicator, Source_CPU);

    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
       IP.Number_of_Processors = CFFC_MPI::Number_of_Processors;
    } /* endif */
    Communicator.Bcast(&(IP.Number_of_Blocks_Per_Processor), 
                       1, 
                       MPI::INT, Source_Rank);

   // Freeze_Limiter
    Communicator.Bcast(&(IP.Freeze_Limiter), 
                       1, 
                       MPI::INT, Source_Rank);
    // Freeze_Limiter_Residual_Level
    Communicator.Bcast(&(IP.Freeze_Limiter_Residual_Level), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    // Accuracy assessment parameters:
    Communicator.Bcast(&(IP.Accuracy_Assessment_Mode), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Accuracy_Assessment_Exact_Digits), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Accuracy_Assessment_Parameter), 
                       1, 
                       MPI::INT, Source_Rank);
}
#endif

/******************************************************//**
 * Routine: Get_Next_Input_Control_Parameter            
 *                                                      
 * Get the next input control parameter from the input  
 * file.                                                
 *                                                      
 ********************************************************/
void Get_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters &IP) {
  return IP.Get_Next_Input_Control_Parameter();
}

/******************************************************//**
 * Routine: Parse_Next_Input_Control_Parameter          
 *                                                      
 * Parses and executes the next input control parameter 
 * from the input file.                                 
 *                                                      
 ********************************************************/
int Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters &IP) {

  int i_command;
  char buffer[256];

  i_command = 0;

  if (strcmp(IP.Next_Control_Parameter, "CFFC_Path") == 0) {
    i_command = 1111;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.CFFC_Path, IP.Next_Control_Parameter);

  } else if (strcmp(IP.Next_Control_Parameter, "Time_Integration_Type") == 0) {
    i_command = 1;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Time_Integration_Type, 
	   IP.Next_Control_Parameter);
    if (strcmp(IP.Time_Integration_Type, "Explicit_Euler") == 0) {
      IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
      IP.N_Stage = 1;
    } else if (strcmp(IP.Time_Integration_Type, "Explicit_Predictor_Corrector") == 0) {
      IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR;
      IP.N_Stage = 2;
    } else if (strcmp(IP.Time_Integration_Type, "Explicit_Runge_Kutta") == 0) {
      IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_RUNGE_KUTTA;
      IP.N_Stage = 4;
    } else if (strcmp(IP.Time_Integration_Type, "Multistage_Optimal_Smoothing") == 0) {
      IP.i_Time_Integration = TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING;
      IP.N_Stage = 4;
      /* Multigrid */
    } else if (strcmp(IP.Time_Integration_Type, "Multigrid") == 0) {
      IP.i_Time_Integration = TIME_STEPPING_MULTIGRID;
    } else {
      std::cout << "\n ==> Unknown time marching scheme!";
      i_command = INVALID_INPUT_VALUE;
    } /* endif */

  } else if (strcmp(IP.Next_Control_Parameter, "Reconstruction_Type") == 0) {
    i_command = 2;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Reconstruction_Type, 
	   IP.Next_Control_Parameter);
    if (strcmp(IP.Reconstruction_Type, "Green_Gauss") == 0) {
      IP.i_Reconstruction = RECONSTRUCTION_GREEN_GAUSS;
      IP.i_ReconstructionMethod = RECONSTRUCTION_GREEN_GAUSS;
      CENO_Execution_Mode::USE_CENO_ALGORITHM = OFF;
    } else if (strcmp(IP.Reconstruction_Type, "Least_Squares") == 0) {
      IP.i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;
      IP.i_ReconstructionMethod = RECONSTRUCTION_LEAST_SQUARES;
      CENO_Execution_Mode::USE_CENO_ALGORITHM = OFF;
    } else if (strcmp(IP.Reconstruction_Type, "CENO") == 0) {
      IP.i_Reconstruction = RECONSTRUCTION_HIGH_ORDER;
      IP.i_ReconstructionMethod = RECONSTRUCTION_CENO;
      CENO_Execution_Mode::USE_CENO_ALGORITHM = ON;
    } else {
      std::cout << "\n ==> Unknown reconstruction method!";
      i_command = INVALID_INPUT_VALUE;
    } /* endif */

  } else if (strcmp(IP.Next_Control_Parameter,"Viscous_Reconstruction_Type") == 0) {
    i_command = 0;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Viscous_Reconstruction_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Viscous_Reconstruction_Type,"Arithmetic_Average") == 0) {
      IP.i_Viscous_Reconstruction = VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE;
    } else if (strcmp(IP.Viscous_Reconstruction_Type,"Diamond_Path") == 0) {
      IP.i_Viscous_Reconstruction = VISCOUS_RECONSTRUCTION_DIAMOND_PATH;
    } else if (strcmp(IP.Viscous_Reconstruction_Type,"Hybrid") == 0) {
      IP.i_Viscous_Reconstruction = VISCOUS_RECONSTRUCTION_HYBRID;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter, "Limiter_Type") == 0) {
    i_command = 3;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Limiter_Type, 
	   IP.Next_Control_Parameter);
    if (strcmp(IP.Limiter_Type, "One") == 0) {
      IP.i_Limiter = LIMITER_ONE;
    } else if (strcmp(IP.Limiter_Type, "Zero") == 0) {
      IP.i_Limiter = LIMITER_ZERO;
    } else if (strcmp(IP.Limiter_Type, "VanLeer") == 0) {
      IP.i_Limiter = LIMITER_VANLEER;
    } else if (strcmp(IP.Limiter_Type, "VanAlbada") == 0) {
      IP.i_Limiter = LIMITER_VANALBADA;
    } else if (strcmp(IP.Limiter_Type, "Barth_Jespersen") == 0) {
      IP.i_Limiter = LIMITER_BARTH_JESPERSEN;
    } else if (strcmp(IP.Limiter_Type, "Venkatakrishnan") == 0) {
      IP.i_Limiter = LIMITER_VENKATAKRISHNAN;
    } else {
      std::cout << "\n ==> Unknown limiter type!";
      i_command = INVALID_INPUT_VALUE;
    } /* endif */

  } else if (strcmp(IP.Next_Control_Parameter, "ICs_Type") == 0) {
    i_command = 4;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.ICs_Type, 
	   IP.Next_Control_Parameter);
    if (strcmp(IP.ICs_Type, "Constant") == 0) {
      IP.i_ICs = IC_CONSTANT;
    } else if (strcmp(IP.ICs_Type, "Uniform") == 0) {
      IP.i_ICs = IC_UNIFORM;
    } else if (strcmp(IP.ICs_Type, "Riemann") == 0) {
      IP.i_ICs = IC_RIEMANN;
    } else if (strcmp(IP.ICs_Type, "Riemann_Xdir") == 0) {
      IP.i_ICs = IC_RIEMANN_XDIR;
    } else if (strcmp(IP.ICs_Type, "Riemann_Ydir") == 0) {
      IP.i_ICs = IC_RIEMANN_YDIR;
    } else if (strcmp(IP.ICs_Type, "Square_Box_IVP") == 0) {
      IP.i_ICs = IC_SQUARE_BOX_IVP;
    } else if (strcmp(IP.ICs_Type, "Circular_Box_IVP") == 0) {
      IP.i_ICs = IC_CIRCULAR_BOX_IVP;
    } else if (strcmp(IP.ICs_Type,"Exact_Solution") == 0) {
      IP.i_ICs = IC_EXACT_SOLUTION;
    } else if (strcmp(IP.ICs_Type,"Uniform_Interior_Exact_Ghost_Cells") == 0) {
      IP.i_ICs = IC_INTERIOR_UNIFORM_GHOSTCELLS_EXACT;
    } else if (strcmp(IP.ICs_Type, "Restart") == 0) {
      IP.i_ICs = IC_RESTART;
    } else {
      std::cout << "\n ==> Unknown initial condition!";
      i_command = INVALID_INPUT_VALUE;
    } /* endif */

  } else if (strcmp(IP.Next_Control_Parameter, "Exact_Integration_Digits") == 0) {
    i_command = 0;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Exact_Integration_Digits;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Exact_Integration_Digits < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Grid_Type") == 0) {
    i_command = 5;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Grid_Type, 
	   IP.Next_Control_Parameter);
    if (strcmp(IP.Grid_Type, "Square") == 0) {
      IP.i_Grid = GRID_SQUARE;
      IP.Box_Width = ONE;
      IP.Box_Height = ONE;
    } else if (strcmp(IP.Grid_Type, "Rectangular_Box") == 0) {
      IP.i_Grid = GRID_RECTANGULAR_BOX;
      IP.Box_Width = ONE;
      IP.Box_Height = ONE;
    } else if (strcmp(IP.Grid_Type, "Deformed_Box") == 0) {
      IP.i_Grid = GRID_DEFORMED_BOX;
    } else if (strcmp(IP.Grid_Type, "Periodic_Box") == 0) {
      IP.i_Grid = GRID_PERIODIC_BOX;
    } else if (strcmp(IP.Grid_Type, "Interior_Inflow_Outflow_Box") == 0) {
      IP.i_Grid = GRID_INTERIOR_INFLOW_OUTFLOW_BOX;
    } else if (strcmp(IP.Grid_Type, "Flat_Plate") == 0) {
      IP.i_Grid = GRID_FLAT_PLATE;
      IP.Plate_Length = ONE;
      IP.BC_South = BC_DIRICHLET;
    } else if (strcmp(IP.Grid_Type, "Pipe") == 0) {
      IP.i_Grid = GRID_PIPE;
      IP.Pipe_Length = ONE;
      IP.Pipe_Radius = HALF;
    } else if (strcmp(IP.Grid_Type, "Blunt_Body") == 0) {
      IP.i_Grid = GRID_BLUNT_BODY;
      IP.Blunt_Body_Radius = ONE;
      IP.Blunt_Body_Mach_Number = TWO;
    } else if (strcmp(IP.Grid_Type, "Rocket_Motor") == 0) {
      IP.i_Grid = GRID_ROCKET_MOTOR;
      IP.Chamber_Length = 0.835;
      IP.Chamber_Radius = 0.020;
      IP.Chamber_To_Throat_Length = 0.05;
      IP.Nozzle_Length = 0.150;
      IP.Nozzle_Radius_Exit = 0.030;
      IP.Nozzle_Radius_Throat = 0.010;
      IP.Nozzle_Type = NOZZLE_GOTTLIEB_FUNCTION;
      IP.Grain_Radius = 0.0;
      IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
      IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
      IP.Mesh_Stretching_Factor_Idir = 1.025;
      IP.Mesh_Stretching_Factor_Jdir = 1.001;
      IP.BC_North = BC_DIRICHLET;
    } else if (strcmp(IP.Grid_Type, "Circular_Cylinder") == 0) {
      IP.i_Grid = GRID_CIRCULAR_CYLINDER;
      IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
      IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
      IP.Mesh_Stretching_Factor_Idir = 1.025;
      IP.Mesh_Stretching_Factor_Jdir = 1.001;
      IP.BC_South = BC_DIRICHLET;
    } else if (strcmp(IP.Grid_Type, "Annulus") == 0) {
      IP.i_Grid = GRID_ANNULUS;
      IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MINMAX_CLUSTERING;
      IP.Mesh_Stretching_Factor_Jdir = 1.01;
    } else if (strcmp(IP.Grid_Type, "Ellipse") == 0) {
      IP.i_Grid = GRID_ELLIPSE;
      IP.Ellipse_Length_X_Axis = TWO;
      IP.Ellipse_Length_Y_Axis = HALF;
    } else if (strcmp(IP.Grid_Type, "NACA_Aerofoil") == 0) {
      IP.i_Grid = GRID_NACA_AEROFOIL;
      IP.Chord_Length = ONE;
      strcpy(IP.NACA_Aerofoil_Type, "0012");
    } else if (strcmp(IP.Grid_Type, "Free_Jet") == 0) {
      IP.i_Grid = GRID_FREE_JET;
      IP.Orifice_Radius = ONE;
    } else if (strcmp(IP.Grid_Type, "ICEMCFD") == 0) {
      IP.i_Grid = GRID_ICEMCFD;
    } else if (strcmp(IP.Grid_Type, "Read_From_Definition_File") == 0) {
      IP.i_Grid = GRID_READ_FROM_DEFINITION_FILE;
    } else if (strcmp(IP.Grid_Type, "Read_From_Data_File") == 0) {
      IP.i_Grid = GRID_READ_FROM_GRID_DATA_FILE;
    } else {
      std::cout << "\n ==> Unknown grid type!";
      i_command = INVALID_INPUT_VALUE;
    } /* endif */

  } else if (strcmp(IP.Next_Control_Parameter, "Output_File_Name") == 0) {
    i_command = 6;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Output_File_Name, 
	   IP.Next_Control_Parameter);
    strcat(IP.Output_File_Name, ".dat");
    strcpy(IP.Grid_File_Name, 
	   IP.Next_Control_Parameter);
    strcat(IP.Grid_File_Name, ".grid");
    strcpy(IP.Grid_Definition_File_Name, 
	   IP.Next_Control_Parameter);
    strcat(IP.Grid_Definition_File_Name, ".griddef");
    strcpy(IP.Restart_File_Name, 
	   IP.Next_Control_Parameter);
    strcat(IP.Restart_File_Name, ".soln");
    strcpy(IP.Gnuplot_File_Name, 
	   IP.Next_Control_Parameter);
    strcat(IP.Gnuplot_File_Name, ".gplt");

  } else if (strcmp(IP.Next_Control_Parameter, "Grid_File_Name") == 0) {
    i_command = 7;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Grid_File_Name, 
	   IP.Next_Control_Parameter);
    strcat(IP.Grid_File_Name, ".grid");
    strcpy(IP.Grid_Definition_File_Name, 
	   IP.Next_Control_Parameter);
    strcat(IP.Grid_Definition_File_Name, ".griddef");

  } else if (strcmp(IP.Next_Control_Parameter, "Restart_File_Name") == 0) {
    i_command = 8;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Restart_File_Name, 
	   IP.Next_Control_Parameter);
    strcat(IP.Restart_File_Name, ".soln");

  } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Cells_Idir") == 0) {
    i_command = 9;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Cells_Idir;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Number_of_Cells_Idir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Cells_Jdir") == 0) {
    i_command = 10;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Cells_Jdir;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Number_of_Cells_Jdir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Ghost_Cells") == 0) {
    i_command = 11;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Ghost_Cells;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Number_of_Ghost_Cells < 2) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Blocks_Idir") == 0) {
    i_command = 11;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Blocks_Idir;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Number_of_Blocks_Idir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Blocks_Jdir") == 0) {
    i_command = 12;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Blocks_Jdir;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Number_of_Blocks_Jdir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Time_Accurate") == 0) {
    i_command = 13;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Time_Accurate;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Time_Accurate != 0 &&
	IP.Time_Accurate != 1) IP.Time_Accurate = 0;
    if (IP.Time_Accurate) {
      IP.Local_Time_Stepping = GLOBAL_TIME_STEPPING;
    } else {
      IP.Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
    } /* endif */

  } else if (strcmp(IP.Next_Control_Parameter, "Local_Time_Stepping") == 0) {
    i_command = 14;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Local_Time_Stepping;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Local_Time_Stepping != GLOBAL_TIME_STEPPING &&
	IP.Local_Time_Stepping != SCALAR_LOCAL_TIME_STEPPING) 
      IP.Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;

  } else if (strcmp(IP.Next_Control_Parameter, "Maximum_Number_of_Time_Steps") == 0) {
    i_command = 15;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Maximum_Number_of_Time_Steps;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Maximum_Number_of_Time_Steps < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "N_Stage") == 0) {
    i_command = 16;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.N_Stage;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.N_Stage < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "CFL_Number") == 0) {
    i_command = 17;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.CFL_Number;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.CFL_Number <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Box_Width") == 0) {
    i_command = 18;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Box_Width;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Box_Width <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Box_Height") == 0) {
    i_command = 19;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Box_Height;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Box_Height <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Plate_Length") == 0) {
    i_command = 20;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Plate_Length;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Plate_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Pipe_Length") == 0) {
    i_command = 21;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Pipe_Length;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Pipe_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Pipe_Radius") == 0) {
    i_command = 22;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Pipe_Radius;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Pipe_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Blunt_Body_Radius") == 0) {
    i_command = 23;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Blunt_Body_Radius;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Blunt_Body_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Blunt_Body_Mach_Number") == 0) {
    i_command = 24;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Blunt_Body_Mach_Number;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Blunt_Body_Mach_Number <= ONE) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Cylinder_Radius") == 0) {
    i_command = 25;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Cylinder_Radius;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Cylinder_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Cylinder_Radius2") == 0) {
    i_command = 26;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Cylinder_Radius2;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Cylinder_Radius2 <= ZERO) i_command = INVALID_INPUT_VALUE;
    if (IP.Cylinder_Radius2 <= IP.Cylinder_Radius) i_command = INVALID_INPUT_VALUE;
    
  } else if (strcmp(IP.Next_Control_Parameter, "Annulus_Start_Angle") == 0) {
    i_command = 26;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Annulus_Theta_Start;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Annulus_End_Angle") == 0) {
    i_command = 26;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Annulus_Theta_End;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Ellipse_Length_X_Axis") == 0) {
    i_command = 26;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Ellipse_Length_X_Axis;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Ellipse_Length_X_Axis <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Ellipse_Length_Y_Axis") == 0) {
    i_command = 27;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Ellipse_Length_Y_Axis;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Ellipse_Length_Y_Axis <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Chord_Length") == 0) {
    i_command = 28;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Chord_Length;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Chord_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "NACA_Aerofoil_Type") == 0) {
    i_command = 29;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.NACA_Aerofoil_Type, 
	   IP.Next_Control_Parameter);
    if (strlen(IP.NACA_Aerofoil_Type) != 4 &&
	strlen(IP.NACA_Aerofoil_Type) != 5) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Orifice_Radius") == 0) {
    i_command = 30;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Orifice_Radius;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Orifice_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Chamber_Length") == 0) {
    i_command = 31;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Chamber_Length;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Chamber_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Chamber_Radius") == 0) {
    i_command = 32;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Chamber_Radius;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Chamber_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Chamber_To_Throat_Length") == 0) {
    i_command = 33;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Chamber_To_Throat_Length;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Chamber_To_Throat_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Nozzle_Length") == 0) {
    i_command = 34;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Nozzle_Length;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Nozzle_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Nozzle_Radius_Exit") == 0) {
    i_command = 35;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Nozzle_Radius_Exit;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Nozzle_Radius_Exit <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Nozzle_Radius_Throat") == 0) {
    i_command = 36;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Nozzle_Radius_Throat;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Nozzle_Radius_Throat <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Nozzle_Type") == 0) {
    i_command = 37;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"Conical") == 0) {
      IP.Nozzle_Type = NOZZLE_CONICAL;
    } else if (strcmp(IP.Next_Control_Parameter,"Gottlieb") == 0) {
      IP.Nozzle_Type = NOZZLE_GOTTLIEB_FUNCTION;
    } else if (strcmp(IP.Next_Control_Parameter,"Quartic") == 0) {
      IP.Nozzle_Type = NOZZLE_QUARTIC_FUNCTION;
    } else if (strcmp(IP.Next_Control_Parameter,"Hybrid_Conical") == 0) {
      IP.Nozzle_Type = NOZZLE_HYBRID_CONICAL_FUNCTION;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Grain_Radius") == 0) {
    i_command = 38;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Grain_Radius;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "VertexSW") == 0) {
    i_command = 0;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.VertexSW;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "VertexSE") == 0) {
    i_command = 0;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.VertexSE;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "VertexNE") == 0) {
    i_command = 0;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.VertexNE;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "VertexNW") == 0) {
    i_command = 0;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.VertexNW;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Time_Max") == 0) {
    i_command = 39;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Time_Max;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Time_Max < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Blocks_Per_Processor") == 0) {
    i_command = 44;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Blocks_Per_Processor;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Number_of_Blocks_Per_Processor < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Output_Format_Type") == 0) {
    i_command = 45;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Output_Format_Type, 
	   IP.Next_Control_Parameter);
    if (strcmp(IP.Output_Format_Type, "Gnuplot") == 0  ||
	strcmp(IP.Output_Format_Type, "gnuplot") == 0  ||
	strcmp(IP.Output_Format_Type, "GNUPLOT") == 0) {
      IP.i_Output_Format = IO_GNUPLOT;
    } else if (strcmp(IP.Output_Format_Type, "Tecplot") == 0  ||
	       strcmp(IP.Output_Format_Type, "tecplot") == 0  ||
	       strcmp(IP.Output_Format_Type, "TECPLOT") == 0) {
      IP.i_Output_Format = IO_TECPLOT;
    } else if (strcmp(IP.Output_Format_Type, "Matlab") == 0  ||
	       strcmp(IP.Output_Format_Type, "matlab") == 0  ||
	       strcmp(IP.Output_Format_Type, "MATLAB") == 0) {
      IP.i_Output_Format = IO_MATLAB;
    } else if (strcmp(IP.Output_Format_Type, "Octave") == 0  ||
	       strcmp(IP.Output_Format_Type, "octave") == 0  ||
	       strcmp(IP.Output_Format_Type, "OCTAVE") == 0) {
      IP.i_Output_Format = IO_OCTAVE;
    } else {
      IP.i_Output_Format = IO_TECPLOT;
    } /* endif */

  } else if (strcmp(IP.Next_Control_Parameter, "Flow_Geometry_Type") == 0) {
    i_command = 46;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Flow_Geometry_Type, 
	   IP.Next_Control_Parameter);
    if (strcmp(IP.Flow_Geometry_Type, "Planar") == 0) {
      IP.Axisymmetric = 0;
    } else if (strcmp(IP.Flow_Geometry_Type, "Axisymmetric") == 0) {
      IP.Axisymmetric = 1;
    } else {
      IP.Axisymmetric = 0;
    } /* endif */

  } else if (strcmp(IP.Next_Control_Parameter, "Include_Source_Term_In_Equation") == 0) {
    i_command = 0;
    IP.Get_Next_Input_Control_Parameter();
    if (strcmp(IP.Next_Control_Parameter, "Yes") == 0) {
      IP.Include_Source_Term = ON;
    } else {
      IP.Include_Source_Term = OFF;
    } /* endif */

  } else if (strcmp(IP.Next_Control_Parameter, "Restart_Solution_Save_Frequency") == 0) {
    i_command = 47;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Restart_Solution_Save_Frequency;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Restart_Solution_Save_Frequency < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "ICEMCFD_Topology_File") == 0) {
    i_command = 48;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.ICEMCFD_FileNames[0], IP.Next_Control_Parameter);

  } else if (strcmp(IP.Next_Control_Parameter, "ICEMCFD_Family_Boco_File") == 0) {
    i_command = 49;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.ICEMCFD_FileNames[1], IP.Next_Control_Parameter);

  } else if (strcmp(IP.Next_Control_Parameter, "ICEMCFD_Family_Topo_File") == 0) {
    i_command = 50;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.ICEMCFD_FileNames[2], IP.Next_Control_Parameter);

  } else if (strcmp(IP.Next_Control_Parameter, "X_Shift") == 0) {
    i_command = 51;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.X_Shift;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "X_Scale") == 0) {
    i_command = 52;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.X_Scale;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "X_Rotate") == 0) {
    i_command = 53;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.X_Rotate;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Freeze_Limiter") == 0) {
    // Freeze_Limiter:
    i_command = 61;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Freeze_Limiter;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Freeze_Limiter < ZERO) i_command = INVALID_INPUT_VALUE;
      
  } else if (strcmp(IP.Next_Control_Parameter, "Freeze_Limiter_Residual_Level") == 0) {
    // Freeze_Limiter_Residual_Level:
    i_command = 62;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Freeze_Limiter_Residual_Level;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Freeze_Limiter_Residual_Level < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "AMR") == 0) {
    i_command = 65;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0 || strcmp(IP.Next_Control_Parameter,"On") == 0) {
      IP.AMR = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0 || strcmp(IP.Next_Control_Parameter,"Off") == 0 ){
      IP.AMR = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter, "AMR_Frequency") == 0) {
    i_command = 66;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.AMR_Frequency;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.AMR_Frequency < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Initial_Mesh_Refinements") == 0) {
    i_command = 67;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Initial_Mesh_Refinements;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Number_of_Initial_Mesh_Refinements < 0) IP.Number_of_Initial_Mesh_Refinements = 0;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Uniform_Mesh_Refinements") == 0) {
    i_command = 68;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Uniform_Mesh_Refinements;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Uniform_Mesh_Refinements < 0) IP.Number_of_Uniform_Mesh_Refinements = 0;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Boundary_Mesh_Refinements") == 0) {
    i_command = 69;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Boundary_Mesh_Refinements;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Boundary_Mesh_Refinements < 0) IP.Number_of_Boundary_Mesh_Refinements = 0;

  } else if (strcmp(IP.Next_Control_Parameter,"Maximum_Refinement_Level") == 0) {
    i_command = 70;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Maximum_Refinement_Level;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Maximum_Refinement_Level < 1) IP.Maximum_Refinement_Level = 1;

  } else if (strcmp(IP.Next_Control_Parameter,"Minimum_Refinement_Level") == 0) {
    i_command = 71;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Minimum_Refinement_Level;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Minimum_Refinement_Level < 1) IP.Minimum_Refinement_Level = 1;

  } else if (strcmp(IP.Next_Control_Parameter, "Threshold_for_Refinement") == 0) {
    i_command = 72;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Threshold_for_Refinement;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Threshold_for_Refinement <= ZERO ||
	IP.Threshold_for_Refinement > ONE) IP.Threshold_for_Refinement = 0.50;

  } else if (strcmp(IP.Next_Control_Parameter, "Threshold_for_Coarsening") == 0) {
    i_command = 71;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Threshold_for_Coarsening;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Threshold_for_Coarsening < ZERO ||
	IP.Threshold_for_Coarsening >= ONE) IP.Threshold_for_Coarsening = 0.10;

  } else if (strcmp(IP.Next_Control_Parameter, "Residual_Smoothing_Epsilon") == 0) {
    i_command = 73;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Residual_Smoothing_Epsilon;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Residual_Smoothing_Epsilon <= ZERO) {
      IP.Residual_Smoothing = 0;
      IP.Residual_Smoothing_Epsilon = ZERO;
    } else {
      IP.Residual_Smoothing = 1;
    } /* endif */

  } else if (strcmp(IP.Next_Control_Parameter, "Residual_Smoothing_Gauss_Seidel_Iterations") == 0) {
    i_command = 74;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Residual_Smoothing_Gauss_Seidel_Iterations;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Residual_Smoothing_Gauss_Seidel_Iterations < 0) {
      IP.Residual_Smoothing_Gauss_Seidel_Iterations = 0;
    } /* endif */

  } else if (strcmp(IP.Next_Control_Parameter,"Mesh_Stretching") == 0) {
    i_command = 101;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0 || strcmp(IP.Next_Control_Parameter,"On") == 0 ) {
      IP.i_Mesh_Stretching = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0 || strcmp(IP.Next_Control_Parameter,"Off") == 0) {
      IP.i_Mesh_Stretching = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Mesh_Stretching_Type_Idir") == 0) {
    i_command = 102;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"Linear") == 0) {
      IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_LINEAR;
    } else if (strcmp(IP.Next_Control_Parameter,"Min_Clustering") == 0) {
      IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MIN_CLUSTERING;
    } else if (strcmp(IP.Next_Control_Parameter,"Max_Clustering") == 0) {
      IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MAX_CLUSTERING;
    } else if (strcmp(IP.Next_Control_Parameter,"MinMax_Clustering") == 0) {
      IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
    } else if (strcmp(IP.Next_Control_Parameter,"Midpt_Clustering") == 0) {
      IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MIDPT_CLUSTERING;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Mesh_Stretching_Type_Jdir") == 0) {
    i_command = 103;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"Linear") == 0) {
      IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_LINEAR;
    } else if (strcmp(IP.Next_Control_Parameter,"Min_Clustering") == 0) {
      IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
    } else if (strcmp(IP.Next_Control_Parameter,"Max_Clustering") == 0) {
      IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MAX_CLUSTERING;
    } else if (strcmp(IP.Next_Control_Parameter,"MinMax_Clustering") == 0) {
      IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MINMAX_CLUSTERING;
    } else if (strcmp(IP.Next_Control_Parameter,"Midpt_Clustering") == 0) {
      IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIDPT_CLUSTERING;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Mesh_Stretching_Factor_Idir") == 0) {
    i_command = 104;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Mesh_Stretching_Factor_Idir;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Mesh_Stretching_Factor_Idir < ONE) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Mesh_Stretching_Factor_Jdir") == 0) {
    i_command = 105;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Mesh_Stretching_Factor_Jdir;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Mesh_Stretching_Factor_Jdir < ONE) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Smooth_Quad_Block") == 0) {
    i_command = 0;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0 || strcmp(IP.Next_Control_Parameter,"On") == 0) {
      IP.i_Smooth_Quad_Block = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0 || strcmp(IP.Next_Control_Parameter,"Off") == 0) {
      IP.i_Smooth_Quad_Block = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Iteration_Disturb_Mesh") == 0) {
    i_command = 0;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.IterationsOfInteriorNodesDisturbances;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.IterationsOfInteriorNodesDisturbances < 0 ) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_Spline_Points") == 0) {
    i_command = 101;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Num_Of_Spline_Control_Points;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Num_Of_Spline_Control_Points <= TWO) i_command = INVALID_INPUT_VALUE;

    /**************************************
     * Multigrid Related Input Parameters */
  } else if (strcmp(IP.Next_Control_Parameter, "Multigrid_Levels") == 0) {
    i_command = 201;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Multigrid_IP.Levels;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Multigrid_IP.Levels <= ONE) i_command = INVALID_INPUT_VALUE;
       
  } else if (strcmp(IP.Next_Control_Parameter, "Multigrid_Cycle_Type") == 0) {
    i_command = 202;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Multigrid_IP.Cycle_Type,
	   IP.Next_Control_Parameter);
    if (strcmp(IP.Multigrid_IP.Cycle_Type, "V") == 0 ||
	strcmp(IP.Multigrid_IP.Cycle_Type, "v") == 0) {
      IP.Multigrid_IP.i_Cycle = MULTIGRID_V_CYCLE;
    } else if (strcmp(IP.Multigrid_IP.Cycle_Type, "W") == 0 ||
	       strcmp(IP.Multigrid_IP.Cycle_Type, "w") == 0) {
      IP.Multigrid_IP.i_Cycle = MULTIGRID_W_CYCLE;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }
  } else if (strcmp(IP.Next_Control_Parameter, "Full_Multigrid") == 0) {
    i_command = 203;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Multigrid_IP.Number_of_Cycles_per_Stage_for_Full_Multigrid;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Multigrid_IP.Number_of_Cycles_per_Stage_for_Full_Multigrid < 0) 
      IP.Multigrid_IP.Number_of_Cycles_per_Stage_for_Full_Multigrid = 0;

  } else if (strcmp(IP.Next_Control_Parameter, "Multigrid_Number_of_Smooths_on_Finest_Level") == 0) {
    i_command = 204;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Multigrid_IP.Number_of_Smooths_on_Finest_Level;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Multigrid_IP.Number_of_Smooths_on_Finest_Level < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Multigrid_Number_of_Pre_Smooths") == 0) {
    i_command = 205;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Multigrid_IP.Number_of_Pre_Smooths;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Multigrid_IP.Number_of_Pre_Smooths < ZERO) i_command = INVALID_INPUT_VALUE;
       
  } else if (strcmp(IP.Next_Control_Parameter, "Multigrid_Number_of_Post_Smooths") == 0) {
    i_command = 206;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Multigrid_IP.Number_of_Post_Smooths;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Multigrid_IP.Number_of_Post_Smooths < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Multigrid_Number_of_Smooths_on_Coarsest_Level") == 0) {
    i_command = 207;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Multigrid_IP.Number_of_Smooths_on_Coarsest_Level;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Multigrid_IP.Number_of_Smooths_on_Coarsest_Level < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Multigrid_Absolute_Convergence_Tolerance") == 0) {
    i_command = 208;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Multigrid_IP.Absolute_Convergence_Tolerance;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Multigrid_Relative_Convergence_Tolerance") == 0) {
    i_command = 209;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Multigrid_IP.Relative_Convergence_Tolerance;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "FMG_Absolute_Convergence_Tolerance") == 0) {
    i_command = 210;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Multigrid_IP.FMG_Absolute_Convergence_Tolerance;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "FMG_Relative_Convergence_Tolerance") == 0) {
    i_command = 211;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Multigrid_IP.FMG_Relative_Convergence_Tolerance;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Multigrid_Smoothing_Type") == 0) {
    i_command = 212;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Multigrid_IP.Smoothing_Type, 
	   IP.Next_Control_Parameter);
    if (strcmp(IP.Multigrid_IP.Smoothing_Type, "Explicit_Euler") == 0) {
      IP.Multigrid_IP.i_Smoothing = TIME_STEPPING_EXPLICIT_EULER;
      IP.N_Stage = 1;
    } else if (strcmp(IP.Multigrid_IP.Smoothing_Type, "Explicit_Predictor_Corrector") == 0) {
      IP.Multigrid_IP.i_Smoothing = TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR;
      IP.N_Stage = 2;
    } else if (strcmp(IP.Multigrid_IP.Smoothing_Type, "Explicit_Runge_Kutta") == 0) {
      IP.Multigrid_IP.i_Smoothing = TIME_STEPPING_EXPLICIT_RUNGE_KUTTA;
      IP.N_Stage = 5;
    } else if (strcmp(IP.Multigrid_IP.Smoothing_Type, "Multistage_Optimal_Smoothing") == 0) {
      IP.Multigrid_IP.i_Smoothing = TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING;
      IP.N_Stage = 4;
    } else {
      i_command = INVALID_INPUT_VALUE;
    } /* endif */

    /* End of Multigrid related Input parsing *
******************************************/

  } else if (strcmp(IP.Next_Control_Parameter,"Boundary_Conditions_Specified") == 0) {
    i_command = 500;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Boundary_Conditions_Specified,IP.Next_Control_Parameter);
    if (strcmp(IP.Boundary_Conditions_Specified,"ON") == 0 || strcmp(IP.Boundary_Conditions_Specified,"On") == 0) {
      IP.BCs_Specified = ON;
    } else if (strcmp(IP.Boundary_Conditions_Specified,"OFF") == 0 || strcmp(IP.Boundary_Conditions_Specified,"Off") == 0) {
      IP.BCs_Specified = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"BC_North") == 0) {
    i_command = 501;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.BC_North_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.BC_North_Type,"Reflection") == 0) {
      IP.BC_North = BC_REFLECTION;
    } else if (strcmp(IP.BC_North_Type,"Burning_Surface") == 0) {
      IP.BC_North = BC_BURNING_SURFACE;
    } else if (strcmp(IP.BC_North_Type,"Wall_Viscous_Isothermal") == 0) {
      IP.BC_North = BC_WALL_VISCOUS_ISOTHERMAL;
    } else if (strcmp(IP.BC_North_Type,"Wall_Viscous_Heatflux") == 0) {
      IP.BC_North = BC_WALL_VISCOUS_HEATFLUX;
    } else if (strcmp(IP.BC_North_Type,"Moving_Wall_Isothermal") == 0) {
      IP.BC_North = BC_MOVING_WALL_ISOTHERMAL;
    } else if (strcmp(IP.BC_North_Type,"Moving_Wall_Heatflux") == 0) {
      IP.BC_North = BC_MOVING_WALL_HEATFLUX;
    } else if (strcmp(IP.BC_North_Type,"Fixed") == 0) {
      IP.BC_North = BC_FIXED;
    } else if (strcmp(IP.BC_North_Type,"Dirichlet") == 0) {
      IP.BC_North = BC_DIRICHLET;
    } else if (strcmp(IP.BC_North_Type,"Neumann") == 0) {
      IP.BC_North = BC_NEUMANN;
    } else if (strcmp(IP.BC_North_Type,"Robin") == 0) {
      IP.BC_North = BC_ROBIN;
    } else if (strcmp(IP.BC_North_Type,"Inflow") == 0) {
      IP.BC_North = BC_INFLOW;
    } else if (strcmp(IP.BC_North_Type,"Outflow") == 0) {
      IP.BC_North = BC_OUTFLOW;
    } else if (strcmp(IP.BC_North_Type,"Farfield") == 0) {
      IP.BC_North = BC_FARFIELD;
    } else if (strcmp(IP.BC_North_Type,"Constant_Extrapolation") == 0) {
      IP.BC_North = BC_CONSTANT_EXTRAPOLATION;
    } else if (strcmp(IP.BC_North_Type,"Linear_Extrapolation") == 0) {
      IP.BC_North = BC_LINEAR_EXTRAPOLATION;
    } else if (strcmp(IP.BC_North_Type,"Frozen") == 0) {
      IP.BC_North = BC_FROZEN;
    } else if (strcmp(IP.BC_North_Type,"Exact_Solution") == 0) {
      IP.BC_North = BC_EXACT_SOLUTION;
    } else if (strcmp(IP.BC_North_Type,"None") == 0) {
      IP.BC_North = BC_NONE;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter, "Ref_State_North") == 0) {
    i_command = 0;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Ref_State_BC_North;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"BC_South") == 0) {
    i_command = 502;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.BC_South_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.BC_South_Type,"Reflection") == 0) {
      IP.BC_South = BC_REFLECTION;
    } else if (strcmp(IP.BC_South_Type,"Burning_Surface") == 0) {
      IP.BC_South = BC_BURNING_SURFACE;
    } else if (strcmp(IP.BC_South_Type,"Wall_Viscous_Isothermal") == 0) {
      IP.BC_South = BC_WALL_VISCOUS_ISOTHERMAL;
    } else if (strcmp(IP.BC_South_Type,"Wall_Viscous_Heatflux") == 0) {
      IP.BC_South = BC_WALL_VISCOUS_HEATFLUX;
    } else if (strcmp(IP.BC_South_Type,"Moving_Wall_Isothermal") == 0) {
      IP.BC_South = BC_MOVING_WALL_ISOTHERMAL;
    } else if (strcmp(IP.BC_South_Type,"Moving_Wall_Heatflux") == 0) {
      IP.BC_South = BC_MOVING_WALL_HEATFLUX;
    } else if (strcmp(IP.BC_South_Type,"Fixed") == 0) {
      IP.BC_South = BC_FIXED;
    } else if (strcmp(IP.BC_South_Type,"Dirichlet") == 0) {
      IP.BC_South = BC_DIRICHLET;
    } else if (strcmp(IP.BC_South_Type,"Neumann") == 0) {
      IP.BC_South = BC_NEUMANN;
    } else if (strcmp(IP.BC_South_Type,"Robin") == 0) {
      IP.BC_South = BC_ROBIN;
    } else if (strcmp(IP.BC_South_Type,"Inflow") == 0) {
      IP.BC_South = BC_INFLOW;
    } else if (strcmp(IP.BC_South_Type,"Outflow") == 0) {
      IP.BC_South = BC_OUTFLOW;
    } else if (strcmp(IP.BC_South_Type,"Farfield") == 0) {
      IP.BC_South = BC_FARFIELD;
    } else if (strcmp(IP.BC_South_Type,"Constant_Extrapolation") == 0) {
      IP.BC_South = BC_CONSTANT_EXTRAPOLATION;
    } else if (strcmp(IP.BC_South_Type,"Linear_Extrapolation") == 0) {
      IP.BC_South = BC_LINEAR_EXTRAPOLATION;
    } else if (strcmp(IP.BC_South_Type,"Frozen") == 0) {
      IP.BC_South = BC_FROZEN;
    } else if (strcmp(IP.BC_South_Type,"Exact_Solution") == 0) {
      IP.BC_South = BC_EXACT_SOLUTION;
    } else if (strcmp(IP.BC_South_Type,"None") == 0) {
      IP.BC_South = BC_NONE;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter, "Ref_State_South") == 0) {
    i_command = 0;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Ref_State_BC_South;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"BC_East") == 0) {
    i_command = 503;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.BC_East_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.BC_East_Type,"Reflection") == 0) {
      IP.BC_East = BC_REFLECTION;
    } else if (strcmp(IP.BC_East_Type,"Burning_Surface") == 0) {
      IP.BC_East = BC_BURNING_SURFACE;
    } else if (strcmp(IP.BC_East_Type,"Wall_Viscous_Isothermal") == 0) {
      IP.BC_East = BC_WALL_VISCOUS_ISOTHERMAL;
    } else if (strcmp(IP.BC_East_Type,"Wall_Viscous_Heatflux") == 0) {
      IP.BC_East = BC_WALL_VISCOUS_HEATFLUX;
    } else if (strcmp(IP.BC_East_Type,"Moving_Wall_Isothermal") == 0) {
      IP.BC_East = BC_MOVING_WALL_ISOTHERMAL;
    } else if (strcmp(IP.BC_East_Type,"Moving_Wall_Heatflux") == 0) {
      IP.BC_East = BC_MOVING_WALL_HEATFLUX;
    } else if (strcmp(IP.BC_East_Type,"Fixed") == 0) {
      IP.BC_East = BC_FIXED;
    } else if (strcmp(IP.BC_East_Type,"Dirichlet") == 0) {
      IP.BC_East = BC_DIRICHLET;
    } else if (strcmp(IP.BC_East_Type,"Neumann") == 0) {
      IP.BC_East = BC_NEUMANN;
    } else if (strcmp(IP.BC_East_Type,"Robin") == 0) {
      IP.BC_East = BC_ROBIN;
    } else if (strcmp(IP.BC_East_Type,"Inflow") == 0) {
      IP.BC_East = BC_INFLOW;
    } else if (strcmp(IP.BC_East_Type,"Outflow") == 0) {
      IP.BC_East = BC_OUTFLOW;
    } else if (strcmp(IP.BC_East_Type,"Farfield") == 0) {
      IP.BC_East = BC_FARFIELD;
    } else if (strcmp(IP.BC_East_Type,"Constant_Extrapolation") == 0) {
      IP.BC_East = BC_CONSTANT_EXTRAPOLATION;
    } else if (strcmp(IP.BC_East_Type,"Linear_Extrapolation") == 0) {
      IP.BC_East = BC_LINEAR_EXTRAPOLATION;
    } else if (strcmp(IP.BC_East_Type,"Frozen") == 0) {
      IP.BC_East = BC_FROZEN;
    } else if (strcmp(IP.BC_East_Type,"Exact_Solution") == 0) {
      IP.BC_East = BC_EXACT_SOLUTION;
    } else if (strcmp(IP.BC_East_Type,"None") == 0) {
      IP.BC_East = BC_NONE;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter, "Ref_State_East") == 0) {
    i_command = 0;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Ref_State_BC_East;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"BC_West") == 0) {
    i_command = 504;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.BC_West_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.BC_West_Type,"Reflection") == 0) {
      IP.BC_West = BC_REFLECTION;
    } else if (strcmp(IP.BC_West_Type,"Burning_Surface") == 0) {
      IP.BC_West = BC_BURNING_SURFACE;
    } else if (strcmp(IP.BC_West_Type,"Wall_Viscous_Isothermal") == 0) {
      IP.BC_West = BC_WALL_VISCOUS_ISOTHERMAL;
    } else if (strcmp(IP.BC_West_Type,"Wall_Viscous_Heatflux") == 0) {
      IP.BC_West = BC_WALL_VISCOUS_HEATFLUX;
    } else if (strcmp(IP.BC_West_Type,"Moving_Wall_Isothermal") == 0) {
      IP.BC_West = BC_MOVING_WALL_ISOTHERMAL;
    } else if (strcmp(IP.BC_West_Type,"Moving_Wall_Heatflux") == 0) {
      IP.BC_West = BC_MOVING_WALL_HEATFLUX;
    } else if (strcmp(IP.BC_West_Type,"Fixed") == 0) {
      IP.BC_West = BC_FIXED;
    } else if (strcmp(IP.BC_West_Type,"Dirichlet") == 0) {
      IP.BC_West = BC_DIRICHLET;
    } else if (strcmp(IP.BC_West_Type,"Neumann") == 0) {
      IP.BC_West = BC_NEUMANN;
    } else if (strcmp(IP.BC_West_Type,"Robin") == 0) {
      IP.BC_West = BC_ROBIN;
    } else if (strcmp(IP.BC_West_Type,"Inflow") == 0) {
      IP.BC_West = BC_INFLOW;
    } else if (strcmp(IP.BC_West_Type,"Outflow") == 0) {
      IP.BC_West = BC_OUTFLOW;
    } else if (strcmp(IP.BC_West_Type,"Farfield") == 0) {
      IP.BC_West = BC_FARFIELD;
    } else if (strcmp(IP.BC_West_Type,"Constant_Extrapolation") == 0) {
      IP.BC_West = BC_CONSTANT_EXTRAPOLATION;
    } else if (strcmp(IP.BC_West_Type,"Linear_Extrapolation") == 0) {
      IP.BC_West = BC_LINEAR_EXTRAPOLATION;
    } else if (strcmp(IP.BC_West_Type,"Frozen") == 0) {
      IP.BC_West = BC_FROZEN;
    } else if (strcmp(IP.BC_West_Type,"Exact_Solution") == 0) {
      IP.BC_West = BC_EXACT_SOLUTION;
    } else if (strcmp(IP.BC_West_Type,"None") == 0) {
      IP.BC_West = BC_NONE;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }
    
  } else if (strcmp(IP.Next_Control_Parameter, "Ref_State_West") == 0) {
    i_command = 0;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Ref_State_BC_West;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Ref_State_Normalization") == 0) {
    i_command = 0;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.RefU;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Space_Accuracy") == 0) {
    i_command = 210;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Space_Accuracy;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Space_Accuracy <= 0 && IP.Space_Accuracy >= 7){
      IP.Space_Accuracy = 1;
      cout << "\n Space Accuracy should be between 1 and 6 \n"
	   << "Space Accuracy set to 1" << endl;
    }/* endif */

  } else if (strcmp(IP.Next_Control_Parameter, "Accuracy_Assessment_Exact_Digits") == 0) {
    i_command = 0;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Accuracy_Assessment_Exact_Digits;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Accuracy_Assessment_Exact_Digits < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter, "Execute") == 0) {
    i_command = EXECUTE_CODE;

  } else if (strcmp(IP.Next_Control_Parameter, "Terminate") == 0) {
    i_command = TERMINATE_CODE;

  } else if (strcmp(IP.Next_Control_Parameter, "Continue") == 0) {
    i_command = CONTINUE_CODE;

  } else if (strcmp(IP.Next_Control_Parameter, "Write_Output") == 0) {
    i_command = WRITE_OUTPUT_CODE;

  } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Cells") == 0) {
    i_command = WRITE_OUTPUT_CELLS_CODE;

  } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Nodes") == 0) {
    i_command = WRITE_OUTPUT_NODES_CODE;

  } else if (strcmp(IP.Next_Control_Parameter, "Write_Restart") == 0) {
    i_command = WRITE_RESTART_CODE;

  } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Mesh") == 0) {
    i_command = WRITE_OUTPUT_GRID_CODE;

  } else if (strcmp(IP.Next_Control_Parameter, "Write_Mesh_Definition") == 0) {
    i_command = WRITE_GRID_DEFINITION_CODE;

  } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Mesh_Nodes") == 0) {
    i_command = WRITE_OUTPUT_GRID_NODES_CODE;

  } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Mesh_Cells") == 0) {
    i_command = WRITE_OUTPUT_GRID_CELLS_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Print_Accuracy") == 0) {
    i_command = WRITE_ERROR_NORMS_TO_SCREEN;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Accuracy_To_File") == 0) {
    i_command = WRITE_ERROR_NORMS_TO_FILE;

  } else if (strcmp(IP.Next_Control_Parameter,"Append_Accuracy_To_File") == 0) {
    i_command = APPEND_ERROR_NORMS_TO_FILE;

  } else if (strcmp(IP.Next_Control_Parameter, "Refine_Grid") == 0) {
    i_command = REFINE_GRID_CODE;

  } else if (IP.Next_Control_Parameter[0] == '#') {
    i_command = COMMENT_CODE;

  } else {
    i_command = INVALID_INPUT_CODE;

  } /* endif */

  

  /* Parse next control parameter with VelocityFields parser */
  VelocityFields::Parse_Next_Input_Control_Parameter(IP,i_command);

  /* Parse next control parameter with DiffusionFields parser */
  DiffusionFields::Parse_Next_Input_Control_Parameter(IP,i_command);

  /* Parse next control parameter with SourceTerm parser */
  IP.SourceTerm->Parse_Next_Input_Control_Parameter(IP,i_command);

  /* Parse next control parameter with ExactSoln parser */
  IP.ExactSoln->Parse_Next_Input_Control_Parameter(IP,i_command);

  /* Parse next control parameter with Inflow parser */
  IP.Inflow->Parse_Next_Input_Control_Parameter(IP,i_command);

  /* Parse next control parameter with CENO_Execution_Mode parser */
  CENO_Execution_Mode::Parse_Next_Input_Control_Parameter(IP,i_command);
  
  /* Parse next control parameter with CENO_Tolerances parser */
  CENO_Tolerances::Parse_Next_Input_Control_Parameter(IP,i_command);

  /* Parse next control parameter with HO_Grid2D_Execution_Mode parser */
  HO_Grid2D_Execution_Mode::Parse_Next_Input_Control_Parameter(IP,i_command);

  /* Parse next control parameter with Tecplot_Execution_Mode parser */
  Tecplot_Execution_Mode::Parse_Next_Input_Control_Parameter(IP,i_command);

  /* Parse next control parameter with HighOrder2D_Input parser */
  HighOrder2D_Input::Parse_Next_Input_Control_Parameter(IP,i_command);

  if (i_command == INVALID_INPUT_CODE){
    // that is, we have an input line which:
    //  - is not a comment (that's COMMENT_CODE), and,
    //  - is not a valid code with an invalid value (that's INVALID_INPUT_VALUE), 
    // and so is an unknown option. Maybe it's an NKS option:
    
    strcpy(buffer, IP.Next_Control_Parameter);
    Get_Next_Input_Control_Parameter(IP);
    i_command = IP.NKS_IP.Parse_Next_Input_Control_Parameter(buffer, 
							     IP.Next_Control_Parameter);
  }

  /* Return the parser command type indicator. */

  return (i_command);

}

/******************************************************//**
 * Routine: Process_Input_Control_Parameter_File        
 *                                                      
 * Reads, parses, and executes the list of input        
 * control parameters from the standard input file.     
 *                                                      
 ********************************************************/
int Process_Input_Control_Parameter_File(AdvectDiffuse2D_Input_Parameters &Input_Parameters,
                                         char *Input_File_Name_ptr,
                                         int &Command_Flag) {

    int error_flag, line_number;

    /* Assign initial value for error indicator flag. */

    error_flag = 0;

    /* Assign default values to the input parameters. */

    Set_Default_Input_Parameters(Input_Parameters);

    /* Copy input file name (a string) to appropriate input parameter variable. */

    if (Input_File_Name_ptr != NULL) strcpy(Input_Parameters.Input_File_Name, Input_File_Name_ptr);

    /* Open the input file containing the input parameters. */

    Open_Input_File(Input_Parameters);
    error_flag = !Input_Parameters.Input_File.is_open();

    if (error_flag) {
       cout << "\n AdvectDiffuse2D ERROR: Unable to open AdvectDiffuse2D input data file.\n";
       return (error_flag);
    } /* endif */

    /* Read and parse control parameters contained in
       the input file. */
  
    while (1) {
       Get_Next_Input_Control_Parameter(Input_Parameters);
       Command_Flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
       line_number = Input_Parameters.Line_Number;
       if (Command_Flag == EXECUTE_CODE) {
          break;
       } else if (Command_Flag == TERMINATE_CODE) {
	 return(0);
       } else if (Command_Flag == INVALID_INPUT_CODE ||
                  Command_Flag == INVALID_INPUT_VALUE) {
          line_number = -line_number;
          cout << "\n AdvectDiffuse2D ERROR: Error reading AdvectDiffuse2D data at line #"
               << -line_number  << " of input data file.\n";
          error_flag = line_number;
          break;
       } /* endif */
    } /* endwhile */

    /* Perform consistency checks on Input_Parameters */

    if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {
      error_flag = Check_Input_Parameters<AdvectDiffuse2D_Input_Parameters>(Input_Parameters);
      if (error_flag) {
	cout << "\n AdvectDiffuse2D ERROR: Input Parameters consistency check failure\n";
	return (error_flag);
      }
    }

    // Perform update of the internal variables of the exact solution
    Input_Parameters.ExactSoln->Set_ParticularSolution_Parameters();

    // Perform update of the internal variables of the inflow field
    Input_Parameters.Inflow->Set_InflowField_Parameters();

    // Perform update of the internal variables of the high-order input parameters
    HighOrder2D_Input::Set_Final_Parameters(Input_Parameters);

    // Set reference states
    // Uo state
    Input_Parameters.Uo = AdvectDiffuse2D_State(ONE);
    // U1 state
    Input_Parameters.U1 = Input_Parameters.Uo;
    Input_Parameters.U1.u = ZERO;

    // U2 state
    Input_Parameters.U2 = Input_Parameters.Uo;
    Input_Parameters.U2.u = -ONE;

    // Set reference state in the AdvectDiffuse2D_Quad_Block class
    AdvectDiffuse2D_Quad_Block::Set_Normalization_Reference_State(Input_Parameters.RefU);

    // Set flag for including/excluding source term in the model equation
    AdvectDiffuse2D_Quad_Block::Include_Source_Term = Input_Parameters.Include_Source_Term;

    /* Initial processing of input control parameters complete.  
       Return the error indicator flag. */

    return (error_flag);

}
