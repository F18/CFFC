/*!\file HO_Grid2DQuadMultiBlock_InputForTesting.cc
  \brief Subroutines for 2D Multi-block Quadrilateral Grid Input for Testing Purposes Class. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "HO_Grid2DQuadMultiBlock_InputForTesting.h"	      // Include definition class header file

/*************************************************************
 * AdvectDiffuse2D_Input_Parameters -- Member functions.     *
 *************************************************************/

/*!
 * \todo Set CENO ghost cells based on the HighOrder2D<int>::Nghost(ReconstructionOrder())
 */
int Grid2DTesting_Input_Parameters::Nghost(void) const{

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
bool Grid2DTesting_Input_Parameters::OutputBoundaryReferenceState(const int & BCtype) const{
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
int Grid2DTesting_Input_Parameters::Parse_Input_File(char *Input_File_Name_ptr){

  ostringstream msg;
  int command_flag;

  strcpy(Input_File_Name, Input_File_Name_ptr);
  Open_Input_File();
  if (Input_File.fail()) {
    msg << "Grid2DTesting_Input_Parameters::Parse_Input_File() ERROR: Unable to open "
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
    command_flag = Parse_Next_Input_Control_Parameter();
    if (command_flag == EXECUTE_CODE) {
      break;
      
    } else if (command_flag == TERMINATE_CODE) {
      return (0);
      
    } else if (command_flag == INVALID_INPUT_CODE ||
	       command_flag == INVALID_INPUT_VALUE) {
      Line_Number = -Line_Number;
      
      msg << "Grid2DTesting_Input_Parameters::Parse_Input_File() ERROR: Error reading data at line # " 
	  << -Line_Number
	  << " of input data file.";
      if (Verbose()){
	cerr << msg.str() << endl;
      }
      
      throw runtime_error(msg.str());
    } /* endif */
  } /* endwhile */

}

/******************************************************//**
 * Get the next input control parameter from the input  
 * file.                                                
 ********************************************************/
void Grid2DTesting_Input_Parameters::Get_Next_Input_Control_Parameter(void){

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
 * Grid2DTesting_Input_Parameters -- Input-output operators. *
 ***************************************************************/
ostream &operator << (ostream &out_file,
		      const Grid2DTesting_Input_Parameters &IP) {
    out_file << setprecision(6);
    out_file << "\n  -> CFFC Path: " 
	     << IP.CFFC_Path;
    out_file << "\n  -> Input File Name: " 
             << IP.Input_File_Name;

    // ====    Boundary conditions ====
    if (IP.BCs_Specified) {
      out_file << "\n  -> Boundary conditions specified as: ";

      // North
      out_file << "\n     -> BC_North = " << IP.BC_North_Type;
      // South
      out_file << "\n     -> BC_South = " << IP.BC_South_Type;
      // East
      out_file << "\n     -> BC_East  = " << IP.BC_East_Type;
      // West
      out_file << "\n     -> BC_West  = " << IP.BC_West_Type;
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
    if (IP.i_ReconstructionMethod != RECONSTRUCTION_CENO ){
      out_file << "\n  -> Elliptic Flux Evaluation: " << IP.Viscous_Reconstruction_Type;
    }
    if (IP.IncludeHighOrderBoundariesRepresentation == OFF){
      out_file << "\n  -> Boundary Accuracy: " << "2nd-Order";
    } else {
      out_file << "\n  -> Boundary Accuracy: " << "high-order";
    }

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
    out_file << "\n  -> Restart Solution Save Frequency: "
             << IP.Restart_Solution_Save_Frequency
             << " steps (iterations)"; 
    out_file.flush();
    return (out_file);
}

istream &operator >> (istream &in_file,
		      Grid2DTesting_Input_Parameters &IP) {
  return (in_file);
}

/*************************************************************
 * Grid2DTesting_Input_Parameters -- External subroutines. *
 *************************************************************/

/******************************************************//**
 * Routine: Open_Input_File                             
 *                                                      
 * Opens the appropriate input data file.               
 *                                                      
 ********************************************************/
void Grid2DTesting_Input_Parameters::Open_Input_File(void) {

    Input_File.open(Input_File_Name, ios::in);
    if (!Input_File.fail()) {
       Line_Number = 0;
       Input_File.setf(ios::skipws);
    } /* endif */

}

/******************************************************//**
 * Routine: Close_Input_File                            
 *                                                      
 * Closes the appropriate input data file.              
 *                                                      
 ********************************************************/
void Grid2DTesting_Input_Parameters::Close_Input_File(void) {

    Input_File.unsetf(ios::skipws);
    Input_File.close();

}

/******************************************************//**
 * Routine: Set_Default_Input_Parameters                
 *                                                      
 * Assigns default values to the input parameters.      
 *                                                      
 ********************************************************/
void Grid2DTesting_Input_Parameters::Set_Default_Input_Parameters(void) {

    int i;
    char *string_ptr;

    // CFFC root directory path:
    get_cffc_path();

    string_ptr = "AdvectDiffuse2D.in";
    strcpy(Input_File_Name, string_ptr);

    // Reconstruction type:
    string_ptr = "Least_Squares";
    strcpy(Reconstruction_Type, string_ptr);
    i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;
    Space_Accuracy = 1;
    IncludeHighOrderBoundariesRepresentation = OFF;
    i_ReconstructionMethod = RECONSTRUCTION_LEAST_SQUARES;

    // Viscous gradient reconstruction type:
    string_ptr = "Diamond_Path";
    strcpy(Viscous_Reconstruction_Type,string_ptr);
    i_Viscous_Reconstruction = VISCOUS_RECONSTRUCTION_DIAMOND_PATH;

    // Geometry switch:
    string_ptr = "Planar";
    strcpy(Flow_Geometry_Type, string_ptr);
    Axisymmetric = 0;

    // Grid parameters:
    string_ptr = "Square";
    strcpy(Grid_Type, string_ptr);
    i_Grid = GRID_SQUARE;
    Box_Width = ONE;
    Box_Height = ONE;
    Number_of_Cells_Idir = 100;
    Number_of_Cells_Jdir = 100;
    Number_of_Ghost_Cells = 2;
    Number_of_Blocks_Idir = 1;
    Number_of_Blocks_Jdir = 1;
    Plate_Length = ONE;
    Pipe_Length = ONE;
    Pipe_Radius = HALF;
    Blunt_Body_Radius = ONE;
    Blunt_Body_Mach_Number = TWO;
    Chamber_Length = 0.835;
    Chamber_Radius = 0.020;
    Chamber_To_Throat_Length = 0.05;
    Nozzle_Length = 0.150;
    Nozzle_Radius_Exit = 0.030;
    Nozzle_Radius_Throat = 0.010;
    Nozzle_Type = NOZZLE_GOTTLIEB_FUNCTION;
    Grain_Radius = 0.0;
    Cylinder_Radius = ONE;
    Cylinder_Radius2 = 32.00;
    Annulus_Theta_Start = 0.0;
    Annulus_Theta_End = 90.0;
    Ellipse_Length_X_Axis = TWO;
    Ellipse_Length_Y_Axis = HALF;
    Chord_Length = ONE;
    Orifice_Radius = ONE;
    VertexSW = Vector2D(-0.5,-0.5);
    VertexSE = Vector2D( 0.5,-0.5);
    VertexNE = Vector2D( 0.5, 0.5);
    VertexNW = Vector2D(-0.5, 0.5);
    X_Shift = Vector2D_ZERO;
    X_Scale = ONE;
    X_Rotate = ZERO;
    i_Smooth_Quad_Block = ON;        // Smooth quad block flag:
    IterationsOfInteriorNodesDisturbances = 0;     /* Number of iterations of disturbing the mesh 
							 (create an unsmooth interior mesh) */
    Num_Of_Spline_Control_Points = 361; /* Number of control points on the 2D spline (used for some grids) */

    // ICEM:
    ICEMCFD_FileNames = ICEMCFD_get_filenames();

    // Mesh stretching factor:
    i_Mesh_Stretching = OFF;
    Mesh_Stretching_Type_Idir = STRETCHING_FCN_LINEAR;
    Mesh_Stretching_Type_Jdir = STRETCHING_FCN_LINEAR;
    Mesh_Stretching_Factor_Idir = 1.01;
    Mesh_Stretching_Factor_Jdir = 1.01;

    // Boundary conditions:
    string_ptr = "OFF";
    strcpy(Boundary_Conditions_Specified,string_ptr);
    BCs_Specified = OFF;
    string_ptr = "None";
    strcpy(BC_North_Type,string_ptr);
    strcpy(BC_South_Type,string_ptr);
    strcpy(BC_East_Type,string_ptr);
    strcpy(BC_West_Type,string_ptr);
    BC_North = BC_NONE;
    BC_South = BC_NONE;
    BC_East  = BC_NONE;
    BC_West  = BC_NONE;

    // AMR:
    AMR = 0;
    AMR_Frequency = 100;
    Number_of_Initial_Mesh_Refinements = 0;
    Number_of_Uniform_Mesh_Refinements = 0;
    Number_of_Boundary_Mesh_Refinements = 0;
    Maximum_Refinement_Level = 100;
    Minimum_Refinement_Level = 1;
    Threshold_for_Refinement = 0.50;
    Threshold_for_Coarsening = 0.10;

    // Default output file names and parameters:
    string_ptr = "outputfile.dat";
    strcpy(Output_File_Name, string_ptr);
    string_ptr = "gridfile.grid";
    strcpy(Grid_File_Name, string_ptr);
    string_ptr = "gridfile.griddef";
    strcpy(Grid_Definition_File_Name, string_ptr);
    string_ptr = "restartfile.soln";
    strcpy(Restart_File_Name, string_ptr);
    string_ptr = "gnuplotfile.gplt";
    strcpy(Gnuplot_File_Name, string_ptr);
    string_ptr = "Tecplot";
    strcpy(Output_Format_Type, string_ptr);
    i_Output_Format = IO_TECPLOT;
    Restart_Solution_Save_Frequency = 1000;

    // Default output progress frequency:
    Output_Progress_Frequency = 50;

    // Input_file parameters:
    string_ptr = " ";
    strcpy(Next_Control_Parameter, string_ptr);
    Line_Number = 0;
    Number_of_Processors = CFFC_MPI::Number_of_Processors;
    Number_of_Blocks_Per_Processor = 10;

}

/******************************************************//**
 * Routine: Parse_Next_Input_Control_Parameter          
 *                                                      
 * Parses and executes the next input control parameter 
 * from the input file.                                 
 *                                                      
 ********************************************************/
int Grid2DTesting_Input_Parameters::Parse_Next_Input_Control_Parameter(void) {

  int i_command;
  char buffer[256];

  i_command = 0;

  if (strcmp(Next_Control_Parameter, "CFFC_Path") == 0) {
    i_command = 1111;
    Get_Next_Input_Control_Parameter();
    strcpy(CFFC_Path, Next_Control_Parameter);

  } else if (strcmp(Next_Control_Parameter, "Reconstruction_Type") == 0) {
    i_command = 2;
    Get_Next_Input_Control_Parameter();
    strcpy(Reconstruction_Type, 
	   Next_Control_Parameter);
    if (strcmp(Reconstruction_Type, "Green_Gauss") == 0) {
      i_Reconstruction = RECONSTRUCTION_GREEN_GAUSS;
      i_ReconstructionMethod = RECONSTRUCTION_GREEN_GAUSS;
    } else if (strcmp(Reconstruction_Type, "Least_Squares") == 0) {
      i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;
      i_ReconstructionMethod = RECONSTRUCTION_LEAST_SQUARES;
    } else if (strcmp(Reconstruction_Type, "CENO") == 0) {
      i_Reconstruction = RECONSTRUCTION_HIGH_ORDER;
      i_ReconstructionMethod = RECONSTRUCTION_CENO;
    } else {
      std::cout << "\n ==> Unknown reconstruction method!";
      i_command = INVALID_INPUT_VALUE;
    } /* endif */

  } else if (strcmp(Next_Control_Parameter,"Viscous_Reconstruction_Type") == 0) {
    i_command = 0;
    Get_Next_Input_Control_Parameter();
    strcpy(Viscous_Reconstruction_Type,Next_Control_Parameter);
    if (strcmp(Viscous_Reconstruction_Type,"Arithmetic_Average") == 0) {
      i_Viscous_Reconstruction = VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE;
    } else if (strcmp(Viscous_Reconstruction_Type,"Diamond_Path") == 0) {
      i_Viscous_Reconstruction = VISCOUS_RECONSTRUCTION_DIAMOND_PATH;
    } else if (strcmp(Viscous_Reconstruction_Type,"Hybrid") == 0) {
      i_Viscous_Reconstruction = VISCOUS_RECONSTRUCTION_HYBRID;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(Next_Control_Parameter, "Grid_Type") == 0) {
    i_command = 5;
    Get_Next_Input_Control_Parameter();
    strcpy(Grid_Type, 
	   Next_Control_Parameter);
    if (strcmp(Grid_Type, "Square") == 0) {
      i_Grid = GRID_SQUARE;
      Box_Width = ONE;
      Box_Height = ONE;
    } else if (strcmp(Grid_Type, "Rectangular_Box") == 0) {
      i_Grid = GRID_RECTANGULAR_BOX;
      Box_Width = ONE;
      Box_Height = ONE;
    } else if (strcmp(Grid_Type, "Deformed_Box") == 0) {
      i_Grid = GRID_DEFORMED_BOX;
    } else if (strcmp(Grid_Type, "Periodic_Box") == 0) {
      i_Grid = GRID_PERIODIC_BOX;
    } else if (strcmp(Grid_Type, "Interior_Inflow_Outflow_Box") == 0) {
      i_Grid = GRID_INTERIOR_INFLOW_OUTFLOW_BOX;
    } else if (strcmp(Grid_Type, "Flat_Plate") == 0) {
      i_Grid = GRID_FLAT_PLATE;
      Plate_Length = ONE;
      BC_South = BC_DIRICHLET;
    } else if (strcmp(Grid_Type, "Pipe") == 0) {
      i_Grid = GRID_PIPE;
      Pipe_Length = ONE;
      Pipe_Radius = HALF;
    } else if (strcmp(Grid_Type, "Blunt_Body") == 0) {
      i_Grid = GRID_BLUNT_BODY;
      Blunt_Body_Radius = ONE;
      Blunt_Body_Mach_Number = TWO;
    } else if (strcmp(Grid_Type, "Rocket_Motor") == 0) {
      i_Grid = GRID_ROCKET_MOTOR;
      Chamber_Length = 0.835;
      Chamber_Radius = 0.020;
      Chamber_To_Throat_Length = 0.05;
      Nozzle_Length = 0.150;
      Nozzle_Radius_Exit = 0.030;
      Nozzle_Radius_Throat = 0.010;
      Nozzle_Type = NOZZLE_GOTTLIEB_FUNCTION;
      Grain_Radius = 0.0;
      Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
      Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
      Mesh_Stretching_Factor_Idir = 1.025;
      Mesh_Stretching_Factor_Jdir = 1.001;
      BC_North = BC_DIRICHLET;
    } else if (strcmp(Grid_Type, "Circular_Cylinder") == 0) {
      i_Grid = GRID_CIRCULAR_CYLINDER;
      Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
      Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
      Mesh_Stretching_Factor_Idir = 1.025;
      Mesh_Stretching_Factor_Jdir = 1.001;
      BC_South = BC_DIRICHLET;
    } else if (strcmp(Grid_Type, "Annulus") == 0) {
      i_Grid = GRID_ANNULUS;
      Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MINMAX_CLUSTERING;
      Mesh_Stretching_Factor_Jdir = 1.01;
    } else if (strcmp(Grid_Type, "Ellipse") == 0) {
      i_Grid = GRID_ELLIPSE;
      Ellipse_Length_X_Axis = TWO;
      Ellipse_Length_Y_Axis = HALF;
    } else if (strcmp(Grid_Type, "NACA_Aerofoil") == 0) {
      i_Grid = GRID_NACA_AEROFOIL;
      Chord_Length = ONE;
      strcpy(NACA_Aerofoil_Type, "0012");
    } else if (strcmp(Grid_Type, "Free_Jet") == 0) {
      i_Grid = GRID_FREE_JET;
      Orifice_Radius = ONE;
    } else if (strcmp(Grid_Type, "ICEMCFD") == 0) {
      i_Grid = GRID_ICEMCFD;
    } else if (strcmp(Grid_Type, "Read_From_Definition_File") == 0) {
      i_Grid = GRID_READ_FROM_DEFINITION_FILE;
    } else if (strcmp(Grid_Type, "Read_From_Data_File") == 0) {
      i_Grid = GRID_READ_FROM_GRID_DATA_FILE;
    } else {
      std::cout << "\n ==> Unknown grid type!";
      i_command = INVALID_INPUT_VALUE;
    } /* endif */

  } else if (strcmp(Next_Control_Parameter, "Output_File_Name") == 0) {
    i_command = 6;
    Get_Next_Input_Control_Parameter();
    strcpy(Output_File_Name, 
	   Next_Control_Parameter);
    strcat(Output_File_Name, ".dat");
    strcpy(Grid_File_Name, 
	   Next_Control_Parameter);
    strcat(Grid_File_Name, ".grid");
    strcpy(Grid_Definition_File_Name, 
	   Next_Control_Parameter);
    strcat(Grid_Definition_File_Name, ".griddef");
    strcpy(Restart_File_Name, 
	   Next_Control_Parameter);
    strcat(Restart_File_Name, ".soln");
    strcpy(Gnuplot_File_Name, 
	   Next_Control_Parameter);
    strcat(Gnuplot_File_Name, ".gplt");

  } else if (strcmp(Next_Control_Parameter, "Grid_File_Name") == 0) {
    i_command = 7;
    Get_Next_Input_Control_Parameter();
    strcpy(Grid_File_Name, 
	   Next_Control_Parameter);
    strcat(Grid_File_Name, ".grid");
    strcpy(Grid_Definition_File_Name, 
	   Next_Control_Parameter);
    strcat(Grid_Definition_File_Name, ".griddef");

  } else if (strcmp(Next_Control_Parameter, "Restart_File_Name") == 0) {
    i_command = 8;
    Get_Next_Input_Control_Parameter();
    strcpy(Restart_File_Name, 
	   Next_Control_Parameter);
    strcat(Restart_File_Name, ".soln");

  } else if (strcmp(Next_Control_Parameter, "Number_of_Cells_Idir") == 0) {
    i_command = 9;
    Line_Number = Line_Number + 1;
    Input_File >> Number_of_Cells_Idir;
    Input_File.getline(buffer, sizeof(buffer));
    if (Number_of_Cells_Idir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter, "Number_of_Cells_Jdir") == 0) {
    i_command = 10;
    Line_Number = Line_Number + 1;
    Input_File >> Number_of_Cells_Jdir;
    Input_File.getline(buffer, sizeof(buffer));
    if (Number_of_Cells_Jdir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter, "Number_of_Ghost_Cells") == 0) {
    i_command = 11;
    Line_Number = Line_Number + 1;
    Input_File >> Number_of_Ghost_Cells;
    Input_File.getline(buffer, sizeof(buffer));
    if (Number_of_Ghost_Cells < 2) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter, "Number_of_Blocks_Idir") == 0) {
    i_command = 11;
    Line_Number = Line_Number + 1;
    Input_File >> Number_of_Blocks_Idir;
    Input_File.getline(buffer, sizeof(buffer));
    if (Number_of_Blocks_Idir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter, "Number_of_Blocks_Jdir") == 0) {
    i_command = 12;
    Line_Number = Line_Number + 1;
    Input_File >> Number_of_Blocks_Jdir;
    Input_File.getline(buffer, sizeof(buffer));
    if (Number_of_Blocks_Jdir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter, "Box_Width") == 0) {
    i_command = 18;
    Line_Number = Line_Number + 1;
    Input_File >> Box_Width;
    Input_File.getline(buffer, sizeof(buffer));
    if (Box_Width <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter, "Box_Height") == 0) {
    i_command = 19;
    Line_Number = Line_Number + 1;
    Input_File >> Box_Height;
    Input_File.getline(buffer, sizeof(buffer));
    if (Box_Height <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter, "Plate_Length") == 0) {
    i_command = 20;
    Line_Number = Line_Number + 1;
    Input_File >> Plate_Length;
    Input_File.getline(buffer, sizeof(buffer));
    if (Plate_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter, "Pipe_Length") == 0) {
    i_command = 21;
    Line_Number = Line_Number + 1;
    Input_File >> Pipe_Length;
    Input_File.getline(buffer, sizeof(buffer));
    if (Pipe_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter, "Pipe_Radius") == 0) {
    i_command = 22;
    Line_Number = Line_Number + 1;
    Input_File >> Pipe_Radius;
    Input_File.getline(buffer, sizeof(buffer));
    if (Pipe_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter, "Blunt_Body_Radius") == 0) {
    i_command = 23;
    Line_Number = Line_Number + 1;
    Input_File >> Blunt_Body_Radius;
    Input_File.getline(buffer, sizeof(buffer));
    if (Blunt_Body_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter, "Blunt_Body_Mach_Number") == 0) {
    i_command = 24;
    Line_Number = Line_Number + 1;
    Input_File >> Blunt_Body_Mach_Number;
    Input_File.getline(buffer, sizeof(buffer));
    if (Blunt_Body_Mach_Number <= ONE) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter, "Cylinder_Radius") == 0) {
    i_command = 25;
    Line_Number = Line_Number + 1;
    Input_File >> Cylinder_Radius;
    Input_File.getline(buffer, sizeof(buffer));
    if (Cylinder_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter, "Cylinder_Radius2") == 0) {
    i_command = 26;
    Line_Number = Line_Number + 1;
    Input_File >> Cylinder_Radius2;
    Input_File.getline(buffer, sizeof(buffer));
    if (Cylinder_Radius2 <= ZERO) i_command = INVALID_INPUT_VALUE;
    if (Cylinder_Radius2 <= Cylinder_Radius) i_command = INVALID_INPUT_VALUE;
    
  } else if (strcmp(Next_Control_Parameter, "Annulus_Start_Angle") == 0) {
    i_command = 26;
    Line_Number = Line_Number + 1;
    Input_File >> Annulus_Theta_Start;
    Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(Next_Control_Parameter, "Annulus_End_Angle") == 0) {
    i_command = 26;
    Line_Number = Line_Number + 1;
    Input_File >> Annulus_Theta_End;
    Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(Next_Control_Parameter, "Ellipse_Length_X_Axis") == 0) {
    i_command = 26;
    Line_Number = Line_Number + 1;
    Input_File >> Ellipse_Length_X_Axis;
    Input_File.getline(buffer, sizeof(buffer));
    if (Ellipse_Length_X_Axis <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter, "Ellipse_Length_Y_Axis") == 0) {
    i_command = 27;
    Line_Number = Line_Number + 1;
    Input_File >> Ellipse_Length_Y_Axis;
    Input_File.getline(buffer, sizeof(buffer));
    if (Ellipse_Length_Y_Axis <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter, "Chord_Length") == 0) {
    i_command = 28;
    Line_Number = Line_Number + 1;
    Input_File >> Chord_Length;
    Input_File.getline(buffer, sizeof(buffer));
    if (Chord_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter, "NACA_Aerofoil_Type") == 0) {
    i_command = 29;
    Get_Next_Input_Control_Parameter();
    strcpy(NACA_Aerofoil_Type, 
	   Next_Control_Parameter);
    if (strlen(NACA_Aerofoil_Type) != 4 &&
	strlen(NACA_Aerofoil_Type) != 5) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter, "Orifice_Radius") == 0) {
    i_command = 30;
    Line_Number = Line_Number + 1;
    Input_File >> Orifice_Radius;
    Input_File.getline(buffer, sizeof(buffer));
    if (Orifice_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter,"Chamber_Length") == 0) {
    i_command = 31;
    Line_Number = Line_Number + 1;
    Input_File >> Chamber_Length;
    Input_File.getline(buffer,sizeof(buffer));
    if (Chamber_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter,"Chamber_Radius") == 0) {
    i_command = 32;
    Line_Number = Line_Number + 1;
    Input_File >> Chamber_Radius;
    Input_File.getline(buffer,sizeof(buffer));
    if (Chamber_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter,"Chamber_To_Throat_Length") == 0) {
    i_command = 33;
    Line_Number = Line_Number + 1;
    Input_File >> Chamber_To_Throat_Length;
    Input_File.getline(buffer,sizeof(buffer));
    if (Chamber_To_Throat_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter,"Nozzle_Length") == 0) {
    i_command = 34;
    Line_Number = Line_Number + 1;
    Input_File >> Nozzle_Length;
    Input_File.getline(buffer,sizeof(buffer));
    if (Nozzle_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter,"Nozzle_Radius_Exit") == 0) {
    i_command = 35;
    Line_Number = Line_Number + 1;
    Input_File >> Nozzle_Radius_Exit;
    Input_File.getline(buffer,sizeof(buffer));
    if (Nozzle_Radius_Exit <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter,"Nozzle_Radius_Throat") == 0) {
    i_command = 36;
    Line_Number = Line_Number + 1;
    Input_File >> Nozzle_Radius_Throat;
    Input_File.getline(buffer,sizeof(buffer));
    if (Nozzle_Radius_Throat <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter,"Nozzle_Type") == 0) {
    i_command = 37;
    Get_Next_Input_Control_Parameter();
    if (strcmp(Next_Control_Parameter,"Conical") == 0) {
      Nozzle_Type = NOZZLE_CONICAL;
    } else if (strcmp(Next_Control_Parameter,"Gottlieb") == 0) {
      Nozzle_Type = NOZZLE_GOTTLIEB_FUNCTION;
    } else if (strcmp(Next_Control_Parameter,"Quartic") == 0) {
      Nozzle_Type = NOZZLE_QUARTIC_FUNCTION;
    } else if (strcmp(Next_Control_Parameter,"Hybrid_Conical") == 0) {
      Nozzle_Type = NOZZLE_HYBRID_CONICAL_FUNCTION;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(Next_Control_Parameter,"Grain_Radius") == 0) {
    i_command = 38;
    Line_Number = Line_Number + 1;
    Input_File >> Grain_Radius;
    Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(Next_Control_Parameter, "VertexSW") == 0) {
    i_command = 0;
    Line_Number = Line_Number + 1;
    Input_File >> VertexSW;
    Input_File.setf(ios::skipws);
    Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(Next_Control_Parameter, "VertexSE") == 0) {
    i_command = 0;
    Line_Number = Line_Number + 1;
    Input_File >> VertexSE;
    Input_File.setf(ios::skipws);
    Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(Next_Control_Parameter, "VertexNE") == 0) {
    i_command = 0;
    Line_Number = Line_Number + 1;
    Input_File >> VertexNE;
    Input_File.setf(ios::skipws);
    Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(Next_Control_Parameter, "VertexNW") == 0) {
    i_command = 0;
    Line_Number = Line_Number + 1;
    Input_File >> VertexNW;
    Input_File.setf(ios::skipws);
    Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(Next_Control_Parameter, "Number_of_Blocks_Per_Processor") == 0) {
    i_command = 44;
    Line_Number = Line_Number + 1;
    Input_File >> Number_of_Blocks_Per_Processor;
    Input_File.getline(buffer, sizeof(buffer));
    if (Number_of_Blocks_Per_Processor < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter, "Output_Format_Type") == 0) {
    i_command = 45;
    Get_Next_Input_Control_Parameter();
    strcpy(Output_Format_Type, 
	   Next_Control_Parameter);
    if (strcmp(Output_Format_Type, "Gnuplot") == 0  ||
	strcmp(Output_Format_Type, "gnuplot") == 0  ||
	strcmp(Output_Format_Type, "GNUPLOT") == 0) {
      i_Output_Format = IO_GNUPLOT;
    } else if (strcmp(Output_Format_Type, "Tecplot") == 0  ||
	       strcmp(Output_Format_Type, "tecplot") == 0  ||
	       strcmp(Output_Format_Type, "TECPLOT") == 0) {
      i_Output_Format = IO_TECPLOT;
    } else if (strcmp(Output_Format_Type, "Matlab") == 0  ||
	       strcmp(Output_Format_Type, "matlab") == 0  ||
	       strcmp(Output_Format_Type, "MATLAB") == 0) {
      i_Output_Format = IO_MATLAB;
    } else if (strcmp(Output_Format_Type, "Octave") == 0  ||
	       strcmp(Output_Format_Type, "octave") == 0  ||
	       strcmp(Output_Format_Type, "OCTAVE") == 0) {
      i_Output_Format = IO_OCTAVE;
    } else {
      i_Output_Format = IO_TECPLOT;
    } /* endif */

  } else if (strcmp(Next_Control_Parameter, "Flow_Geometry_Type") == 0) {
    i_command = 46;
    Get_Next_Input_Control_Parameter();
    strcpy(Flow_Geometry_Type, 
	   Next_Control_Parameter);
    if (strcmp(Flow_Geometry_Type, "Planar") == 0) {
      Axisymmetric = 0;
    } else if (strcmp(Flow_Geometry_Type, "Axisymmetric") == 0) {
      Axisymmetric = 1;
    } else {
      Axisymmetric = 0;
    } /* endif */

  } else if (strcmp(Next_Control_Parameter, "Restart_Solution_Save_Frequency") == 0) {
    i_command = 47;
    Line_Number = Line_Number + 1;
    Input_File >> Restart_Solution_Save_Frequency;
    Input_File.getline(buffer, sizeof(buffer));
    if (Restart_Solution_Save_Frequency < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter, "ICEMCFD_Topology_File") == 0) {
    i_command = 48;
    Get_Next_Input_Control_Parameter();
    strcpy(ICEMCFD_FileNames[0], Next_Control_Parameter);

  } else if (strcmp(Next_Control_Parameter, "ICEMCFD_Family_Boco_File") == 0) {
    i_command = 49;
    Get_Next_Input_Control_Parameter();
    strcpy(ICEMCFD_FileNames[1], Next_Control_Parameter);

  } else if (strcmp(Next_Control_Parameter, "ICEMCFD_Family_Topo_File") == 0) {
    i_command = 50;
    Get_Next_Input_Control_Parameter();
    strcpy(ICEMCFD_FileNames[2], Next_Control_Parameter);

  } else if (strcmp(Next_Control_Parameter, "X_Shift") == 0) {
    i_command = 51;
    Line_Number = Line_Number + 1;
    Input_File >> X_Shift;
    Input_File.setf(ios::skipws);
    Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(Next_Control_Parameter, "X_Scale") == 0) {
    i_command = 52;
    Line_Number = Line_Number + 1;
    Input_File >> X_Scale;
    Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(Next_Control_Parameter, "X_Rotate") == 0) {
    i_command = 53;
    Line_Number = Line_Number + 1;
    Input_File >> X_Rotate;
    Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(Next_Control_Parameter, "AMR") == 0) {
    i_command = 65;
    Get_Next_Input_Control_Parameter();
    if (strcmp(Next_Control_Parameter,"ON") == 0 || strcmp(Next_Control_Parameter,"On") == 0) {
      AMR = ON;
    } else if (strcmp(Next_Control_Parameter,"OFF") == 0 || strcmp(Next_Control_Parameter,"Off") == 0 ){
      AMR = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(Next_Control_Parameter, "AMR_Frequency") == 0) {
    i_command = 66;
    Line_Number = Line_Number + 1;
    Input_File >> AMR_Frequency;
    Input_File.getline(buffer, sizeof(buffer));
    if (AMR_Frequency < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter, "Number_of_Initial_Mesh_Refinements") == 0) {
    i_command = 67;
    Line_Number = Line_Number + 1;
    Input_File >> Number_of_Initial_Mesh_Refinements;
    Input_File.getline(buffer, sizeof(buffer));
    if (Number_of_Initial_Mesh_Refinements < 0) Number_of_Initial_Mesh_Refinements = 0;

  } else if (strcmp(Next_Control_Parameter,"Number_of_Uniform_Mesh_Refinements") == 0) {
    i_command = 68;
    Line_Number = Line_Number + 1;
    Input_File >> Number_of_Uniform_Mesh_Refinements;
    Input_File.getline(buffer,sizeof(buffer));
    if (Number_of_Uniform_Mesh_Refinements < 0) Number_of_Uniform_Mesh_Refinements = 0;

  } else if (strcmp(Next_Control_Parameter,"Number_of_Boundary_Mesh_Refinements") == 0) {
    i_command = 69;
    Line_Number = Line_Number + 1;
    Input_File >> Number_of_Boundary_Mesh_Refinements;
    Input_File.getline(buffer,sizeof(buffer));
    if (Number_of_Boundary_Mesh_Refinements < 0) Number_of_Boundary_Mesh_Refinements = 0;

  } else if (strcmp(Next_Control_Parameter,"Maximum_Refinement_Level") == 0) {
    i_command = 70;
    Line_Number = Line_Number + 1;
    Input_File >> Maximum_Refinement_Level;
    Input_File.getline(buffer,sizeof(buffer));
    if (Maximum_Refinement_Level < 1) Maximum_Refinement_Level = 1;

  } else if (strcmp(Next_Control_Parameter,"Minimum_Refinement_Level") == 0) {
    i_command = 71;
    Line_Number = Line_Number + 1;
    Input_File >> Minimum_Refinement_Level;
    Input_File.getline(buffer,sizeof(buffer));
    if (Minimum_Refinement_Level < 1) Minimum_Refinement_Level = 1;

  } else if (strcmp(Next_Control_Parameter, "Threshold_for_Refinement") == 0) {
    i_command = 72;
    Line_Number = Line_Number + 1;
    Input_File >> Threshold_for_Refinement;
    Input_File.getline(buffer, sizeof(buffer));
    if (Threshold_for_Refinement <= ZERO ||
	Threshold_for_Refinement > ONE) Threshold_for_Refinement = 0.50;

  } else if (strcmp(Next_Control_Parameter, "Threshold_for_Coarsening") == 0) {
    i_command = 71;
    Line_Number = Line_Number + 1;
    Input_File >> Threshold_for_Coarsening;
    Input_File.getline(buffer, sizeof(buffer));
    if (Threshold_for_Coarsening < ZERO ||
	Threshold_for_Coarsening >= ONE) Threshold_for_Coarsening = 0.10;

  } else if (strcmp(Next_Control_Parameter,"Mesh_Stretching") == 0) {
    i_command = 101;
    Get_Next_Input_Control_Parameter();
    if (strcmp(Next_Control_Parameter,"ON") == 0 || strcmp(Next_Control_Parameter,"On") == 0 ) {
      i_Mesh_Stretching = ON;
    } else if (strcmp(Next_Control_Parameter,"OFF") == 0 || strcmp(Next_Control_Parameter,"Off") == 0) {
      i_Mesh_Stretching = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(Next_Control_Parameter,"Mesh_Stretching_Type_Idir") == 0) {
    i_command = 102;
    Get_Next_Input_Control_Parameter();
    if (strcmp(Next_Control_Parameter,"Linear") == 0) {
      Mesh_Stretching_Type_Idir = STRETCHING_FCN_LINEAR;
    } else if (strcmp(Next_Control_Parameter,"Min_Clustering") == 0) {
      Mesh_Stretching_Type_Idir = STRETCHING_FCN_MIN_CLUSTERING;
    } else if (strcmp(Next_Control_Parameter,"Max_Clustering") == 0) {
      Mesh_Stretching_Type_Idir = STRETCHING_FCN_MAX_CLUSTERING;
    } else if (strcmp(Next_Control_Parameter,"MinMax_Clustering") == 0) {
      Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
    } else if (strcmp(Next_Control_Parameter,"Midpt_Clustering") == 0) {
      Mesh_Stretching_Type_Idir = STRETCHING_FCN_MIDPT_CLUSTERING;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(Next_Control_Parameter,"Mesh_Stretching_Type_Jdir") == 0) {
    i_command = 103;
    Get_Next_Input_Control_Parameter();
    if (strcmp(Next_Control_Parameter,"Linear") == 0) {
      Mesh_Stretching_Type_Jdir = STRETCHING_FCN_LINEAR;
    } else if (strcmp(Next_Control_Parameter,"Min_Clustering") == 0) {
      Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
    } else if (strcmp(Next_Control_Parameter,"Max_Clustering") == 0) {
      Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MAX_CLUSTERING;
    } else if (strcmp(Next_Control_Parameter,"MinMax_Clustering") == 0) {
      Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MINMAX_CLUSTERING;
    } else if (strcmp(Next_Control_Parameter,"Midpt_Clustering") == 0) {
      Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIDPT_CLUSTERING;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(Next_Control_Parameter,"Mesh_Stretching_Factor_Idir") == 0) {
    i_command = 104;
    Line_Number = Line_Number + 1;
    Input_File >> Mesh_Stretching_Factor_Idir;
    Input_File.getline(buffer,sizeof(buffer));
    if (Mesh_Stretching_Factor_Idir < ONE) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter,"Mesh_Stretching_Factor_Jdir") == 0) {
    i_command = 105;
    Line_Number = Line_Number + 1;
    Input_File >> Mesh_Stretching_Factor_Jdir;
    Input_File.getline(buffer,sizeof(buffer));
    if (Mesh_Stretching_Factor_Jdir < ONE) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter,"Smooth_Quad_Block") == 0) {
    i_command = 0;
    Get_Next_Input_Control_Parameter();
    if (strcmp(Next_Control_Parameter,"ON") == 0 || strcmp(Next_Control_Parameter,"On") == 0) {
      i_Smooth_Quad_Block = ON;
    } else if (strcmp(Next_Control_Parameter,"OFF") == 0 || strcmp(Next_Control_Parameter,"Off") == 0) {
      i_Smooth_Quad_Block = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(Next_Control_Parameter,"Iteration_Disturb_Mesh") == 0) {
    i_command = 0;
    Line_Number = Line_Number + 1;
    Input_File >> IterationsOfInteriorNodesDisturbances;
    Input_File.getline(buffer, sizeof(buffer));
    if (IterationsOfInteriorNodesDisturbances < 0 ) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter,"Number_Spline_Points") == 0) {
    i_command = 101;
    Line_Number = Line_Number + 1;
    Input_File >> Num_Of_Spline_Control_Points;
    Input_File.getline(buffer, sizeof(buffer));
    if (Num_Of_Spline_Control_Points <= TWO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(Next_Control_Parameter,"Boundary_Conditions_Specified") == 0) {
    i_command = 500;
    Get_Next_Input_Control_Parameter();
    strcpy(Boundary_Conditions_Specified,Next_Control_Parameter);
    if (strcmp(Boundary_Conditions_Specified,"ON") == 0 || strcmp(Boundary_Conditions_Specified,"On") == 0) {
      BCs_Specified = ON;
    } else if (strcmp(Boundary_Conditions_Specified,"OFF") == 0 || strcmp(Boundary_Conditions_Specified,"Off") == 0) {
      BCs_Specified = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(Next_Control_Parameter,"BC_North") == 0) {
    i_command = 501;
    Get_Next_Input_Control_Parameter();
    strcpy(BC_North_Type,Next_Control_Parameter);
    if (strcmp(BC_North_Type,"Reflection") == 0) {
      BC_North = BC_REFLECTION;
    } else if (strcmp(BC_North_Type,"Burning_Surface") == 0) {
      BC_North = BC_BURNING_SURFACE;
    } else if (strcmp(BC_North_Type,"Wall_Viscous_Isothermal") == 0) {
      BC_North = BC_WALL_VISCOUS_ISOTHERMAL;
    } else if (strcmp(BC_North_Type,"Wall_Viscous_Heatflux") == 0) {
      BC_North = BC_WALL_VISCOUS_HEATFLUX;
    } else if (strcmp(BC_North_Type,"Moving_Wall_Isothermal") == 0) {
      BC_North = BC_MOVING_WALL_ISOTHERMAL;
    } else if (strcmp(BC_North_Type,"Moving_Wall_Heatflux") == 0) {
      BC_North = BC_MOVING_WALL_HEATFLUX;
    } else if (strcmp(BC_North_Type,"Fixed") == 0) {
      BC_North = BC_FIXED;
    } else if (strcmp(BC_North_Type,"Dirichlet") == 0) {
      BC_North = BC_DIRICHLET;
    } else if (strcmp(BC_North_Type,"Neumann") == 0) {
      BC_North = BC_NEUMANN;
    } else if (strcmp(BC_North_Type,"Robin") == 0) {
      BC_North = BC_ROBIN;
    } else if (strcmp(BC_North_Type,"Inflow") == 0) {
      BC_North = BC_INFLOW;
    } else if (strcmp(BC_North_Type,"Outflow") == 0) {
      BC_North = BC_OUTFLOW;
    } else if (strcmp(BC_North_Type,"Farfield") == 0) {
      BC_North = BC_FARFIELD;
    } else if (strcmp(BC_North_Type,"Constant_Extrapolation") == 0) {
      BC_North = BC_CONSTANT_EXTRAPOLATION;
    } else if (strcmp(BC_North_Type,"Linear_Extrapolation") == 0) {
      BC_North = BC_LINEAR_EXTRAPOLATION;
    } else if (strcmp(BC_North_Type,"Frozen") == 0) {
      BC_North = BC_FROZEN;
    } else if (strcmp(BC_North_Type,"Exact_Solution") == 0) {
      BC_North = BC_EXACT_SOLUTION;
    } else if (strcmp(BC_North_Type,"None") == 0) {
      BC_North = BC_NONE;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(Next_Control_Parameter,"BC_South") == 0) {
    i_command = 502;
    Get_Next_Input_Control_Parameter();
    strcpy(BC_South_Type,Next_Control_Parameter);
    if (strcmp(BC_South_Type,"Reflection") == 0) {
      BC_South = BC_REFLECTION;
    } else if (strcmp(BC_South_Type,"Burning_Surface") == 0) {
      BC_South = BC_BURNING_SURFACE;
    } else if (strcmp(BC_South_Type,"Wall_Viscous_Isothermal") == 0) {
      BC_South = BC_WALL_VISCOUS_ISOTHERMAL;
    } else if (strcmp(BC_South_Type,"Wall_Viscous_Heatflux") == 0) {
      BC_South = BC_WALL_VISCOUS_HEATFLUX;
    } else if (strcmp(BC_South_Type,"Moving_Wall_Isothermal") == 0) {
      BC_South = BC_MOVING_WALL_ISOTHERMAL;
    } else if (strcmp(BC_South_Type,"Moving_Wall_Heatflux") == 0) {
      BC_South = BC_MOVING_WALL_HEATFLUX;
    } else if (strcmp(BC_South_Type,"Fixed") == 0) {
      BC_South = BC_FIXED;
    } else if (strcmp(BC_South_Type,"Dirichlet") == 0) {
      BC_South = BC_DIRICHLET;
    } else if (strcmp(BC_South_Type,"Neumann") == 0) {
      BC_South = BC_NEUMANN;
    } else if (strcmp(BC_South_Type,"Robin") == 0) {
      BC_South = BC_ROBIN;
    } else if (strcmp(BC_South_Type,"Inflow") == 0) {
      BC_South = BC_INFLOW;
    } else if (strcmp(BC_South_Type,"Outflow") == 0) {
      BC_South = BC_OUTFLOW;
    } else if (strcmp(BC_South_Type,"Farfield") == 0) {
      BC_South = BC_FARFIELD;
    } else if (strcmp(BC_South_Type,"Constant_Extrapolation") == 0) {
      BC_South = BC_CONSTANT_EXTRAPOLATION;
    } else if (strcmp(BC_South_Type,"Linear_Extrapolation") == 0) {
      BC_South = BC_LINEAR_EXTRAPOLATION;
    } else if (strcmp(BC_South_Type,"Frozen") == 0) {
      BC_South = BC_FROZEN;
    } else if (strcmp(BC_South_Type,"Exact_Solution") == 0) {
      BC_South = BC_EXACT_SOLUTION;
    } else if (strcmp(BC_South_Type,"None") == 0) {
      BC_South = BC_NONE;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(Next_Control_Parameter,"BC_East") == 0) {
    i_command = 503;
    Get_Next_Input_Control_Parameter();
    strcpy(BC_East_Type,Next_Control_Parameter);
    if (strcmp(BC_East_Type,"Reflection") == 0) {
      BC_East = BC_REFLECTION;
    } else if (strcmp(BC_East_Type,"Burning_Surface") == 0) {
      BC_East = BC_BURNING_SURFACE;
    } else if (strcmp(BC_East_Type,"Wall_Viscous_Isothermal") == 0) {
      BC_East = BC_WALL_VISCOUS_ISOTHERMAL;
    } else if (strcmp(BC_East_Type,"Wall_Viscous_Heatflux") == 0) {
      BC_East = BC_WALL_VISCOUS_HEATFLUX;
    } else if (strcmp(BC_East_Type,"Moving_Wall_Isothermal") == 0) {
      BC_East = BC_MOVING_WALL_ISOTHERMAL;
    } else if (strcmp(BC_East_Type,"Moving_Wall_Heatflux") == 0) {
      BC_East = BC_MOVING_WALL_HEATFLUX;
    } else if (strcmp(BC_East_Type,"Fixed") == 0) {
      BC_East = BC_FIXED;
    } else if (strcmp(BC_East_Type,"Dirichlet") == 0) {
      BC_East = BC_DIRICHLET;
    } else if (strcmp(BC_East_Type,"Neumann") == 0) {
      BC_East = BC_NEUMANN;
    } else if (strcmp(BC_East_Type,"Robin") == 0) {
      BC_East = BC_ROBIN;
    } else if (strcmp(BC_East_Type,"Inflow") == 0) {
      BC_East = BC_INFLOW;
    } else if (strcmp(BC_East_Type,"Outflow") == 0) {
      BC_East = BC_OUTFLOW;
    } else if (strcmp(BC_East_Type,"Farfield") == 0) {
      BC_East = BC_FARFIELD;
    } else if (strcmp(BC_East_Type,"Constant_Extrapolation") == 0) {
      BC_East = BC_CONSTANT_EXTRAPOLATION;
    } else if (strcmp(BC_East_Type,"Linear_Extrapolation") == 0) {
      BC_East = BC_LINEAR_EXTRAPOLATION;
    } else if (strcmp(BC_East_Type,"Frozen") == 0) {
      BC_East = BC_FROZEN;
    } else if (strcmp(BC_East_Type,"Exact_Solution") == 0) {
      BC_East = BC_EXACT_SOLUTION;
    } else if (strcmp(BC_East_Type,"None") == 0) {
      BC_East = BC_NONE;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(Next_Control_Parameter,"BC_West") == 0) {
    i_command = 504;
    Get_Next_Input_Control_Parameter();
    strcpy(BC_West_Type,Next_Control_Parameter);
    if (strcmp(BC_West_Type,"Reflection") == 0) {
      BC_West = BC_REFLECTION;
    } else if (strcmp(BC_West_Type,"Burning_Surface") == 0) {
      BC_West = BC_BURNING_SURFACE;
    } else if (strcmp(BC_West_Type,"Wall_Viscous_Isothermal") == 0) {
      BC_West = BC_WALL_VISCOUS_ISOTHERMAL;
    } else if (strcmp(BC_West_Type,"Wall_Viscous_Heatflux") == 0) {
      BC_West = BC_WALL_VISCOUS_HEATFLUX;
    } else if (strcmp(BC_West_Type,"Moving_Wall_Isothermal") == 0) {
      BC_West = BC_MOVING_WALL_ISOTHERMAL;
    } else if (strcmp(BC_West_Type,"Moving_Wall_Heatflux") == 0) {
      BC_West = BC_MOVING_WALL_HEATFLUX;
    } else if (strcmp(BC_West_Type,"Fixed") == 0) {
      BC_West = BC_FIXED;
    } else if (strcmp(BC_West_Type,"Dirichlet") == 0) {
      BC_West = BC_DIRICHLET;
    } else if (strcmp(BC_West_Type,"Neumann") == 0) {
      BC_West = BC_NEUMANN;
    } else if (strcmp(BC_West_Type,"Robin") == 0) {
      BC_West = BC_ROBIN;
    } else if (strcmp(BC_West_Type,"Inflow") == 0) {
      BC_West = BC_INFLOW;
    } else if (strcmp(BC_West_Type,"Outflow") == 0) {
      BC_West = BC_OUTFLOW;
    } else if (strcmp(BC_West_Type,"Farfield") == 0) {
      BC_West = BC_FARFIELD;
    } else if (strcmp(BC_West_Type,"Constant_Extrapolation") == 0) {
      BC_West = BC_CONSTANT_EXTRAPOLATION;
    } else if (strcmp(BC_West_Type,"Linear_Extrapolation") == 0) {
      BC_West = BC_LINEAR_EXTRAPOLATION;
    } else if (strcmp(BC_West_Type,"Frozen") == 0) {
      BC_West = BC_FROZEN;
    } else if (strcmp(BC_West_Type,"Exact_Solution") == 0) {
      BC_West = BC_EXACT_SOLUTION;
    } else if (strcmp(BC_West_Type,"None") == 0) {
      BC_West = BC_NONE;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }
    
  } else if (strcmp(Next_Control_Parameter, "Space_Accuracy") == 0) {
    i_command = 210;
    Line_Number = Line_Number + 1;
    Input_File >> Space_Accuracy;
    Input_File.getline(buffer, sizeof(buffer));
    if (Space_Accuracy <= 0 && Space_Accuracy >= 7){
      Space_Accuracy = 1;
      cout << "\n Space Accuracy should be between 1 and 6 \n"
	   << "Space Accuracy set to 1" << endl;
    }/* endif */

  } else if (strcmp(Next_Control_Parameter, "High_Order_Boundary") == 0) {
    i_command = 0;
    Get_Next_Input_Control_Parameter();
    if (strcmp(Next_Control_Parameter,"ON") == 0 || strcmp(Next_Control_Parameter,"On") == 0) {
      IncludeHighOrderBoundariesRepresentation = ON;
    } else if (strcmp(Next_Control_Parameter,"OFF") == 0 || strcmp(Next_Control_Parameter,"Off") == 0) {
      IncludeHighOrderBoundariesRepresentation = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(Next_Control_Parameter, "Execute") == 0) {
    i_command = EXECUTE_CODE;

  } else if (strcmp(Next_Control_Parameter, "Terminate") == 0) {
    i_command = TERMINATE_CODE;

  } else if (strcmp(Next_Control_Parameter, "Continue") == 0) {
    i_command = CONTINUE_CODE;

  } else if (strcmp(Next_Control_Parameter, "Write_Output") == 0) {
    i_command = WRITE_OUTPUT_CODE;

  } else if (strcmp(Next_Control_Parameter, "Write_Output_Cells") == 0) {
    i_command = WRITE_OUTPUT_CELLS_CODE;

  } else if (strcmp(Next_Control_Parameter, "Write_Output_Nodes") == 0) {
    i_command = WRITE_OUTPUT_NODES_CODE;

  } else if (strcmp(Next_Control_Parameter, "Write_Restart") == 0) {
    i_command = WRITE_RESTART_CODE;

  } else if (strcmp(Next_Control_Parameter, "Write_Output_Mesh") == 0) {
    i_command = WRITE_OUTPUT_GRID_CODE;

  } else if (strcmp(Next_Control_Parameter, "Write_Mesh_Definition") == 0) {
    i_command = WRITE_GRID_DEFINITION_CODE;

  } else if (strcmp(Next_Control_Parameter, "Write_Output_Mesh_Nodes") == 0) {
    i_command = WRITE_OUTPUT_GRID_NODES_CODE;

  } else if (strcmp(Next_Control_Parameter, "Write_Output_Mesh_Cells") == 0) {
    i_command = WRITE_OUTPUT_GRID_CELLS_CODE;

  } else if (strcmp(Next_Control_Parameter,"Print_Accuracy") == 0) {
    i_command = WRITE_ERROR_NORMS_TO_SCREEN;

  } else if (strcmp(Next_Control_Parameter,"Write_Accuracy_To_File") == 0) {
    i_command = WRITE_ERROR_NORMS_TO_FILE;

  } else if (strcmp(Next_Control_Parameter,"Append_Accuracy_To_File") == 0) {
    i_command = APPEND_ERROR_NORMS_TO_FILE;

  } else if (strcmp(Next_Control_Parameter, "Refine_Grid") == 0) {
    i_command = REFINE_GRID_CODE;

  } else if (Next_Control_Parameter[0] == '#') {
    i_command = COMMENT_CODE;

  } else {
    i_command = INVALID_INPUT_CODE;

  } /* endif */

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
int Grid2DTesting_Input_Parameters::Process_Input_Control_Parameter_File(char *Input_File_Name_ptr,
									 int &Command_Flag) {

    int error_flag, line_number;

    /* Assign initial value for error indicator flag. */

    error_flag = 0;

    /* Assign default values to the input parameters. */

    Set_Default_Input_Parameters();

    /* Copy input file name (a string) to appropriate input parameter variable. */

    if (Input_File_Name_ptr != NULL) strcpy(Input_File_Name, Input_File_Name_ptr);

    /* Open the input file containing the input parameters. */

    Open_Input_File();
    error_flag = !Input_File.is_open();

    if (error_flag) {
       cout << "\n AdvectDiffuse2D ERROR: Unable to open AdvectDiffuse2D input data file.\n";
       return (error_flag);
    } /* endif */

    /* Read and parse control parameters contained in
       the input file. */
  
    while (1) {
       Get_Next_Input_Control_Parameter();
       Command_Flag = Parse_Next_Input_Control_Parameter();
       line_number = Line_Number;
       if (Command_Flag == EXECUTE_CODE) {
          break;
       } else if (Command_Flag == TERMINATE_CODE) {
	 return(0);
       } else if (Command_Flag == INVALID_INPUT_CODE ||
                  Command_Flag == INVALID_INPUT_VALUE) {
          line_number = -line_number;
          cout << "\n Grid2DQuadMultiBlock Testing ERROR: Error reading Grid2DQuadMultiBlock data at line #"
               << -line_number  << " of input data file.\n";
          error_flag = line_number;
          break;
       } /* endif */
    } /* endwhile */

    /* Initial processing of input control parameters complete.  
       Return the error indicator flag. */

    return (error_flag);

}
