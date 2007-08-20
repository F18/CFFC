/**********************************************************************
 * Electrostatic2DInput.cc: Subroutines for 2D Electrostatic input    *
 *                          classes.                                  *
 **********************************************************************/

// Include 2D Electrostatic input parameter header file.

#ifndef _ELECTROSTATIC2D_INPUT_INCLUDED
#include "Electrostatic2DInput.h"
#endif // _ELECTROSTATIC2D_INPUT_INCLUDED

/**********************************************************************
 * Electrostatic2D_Input_Parameters -- External subroutines.          *
 **********************************************************************/

/**********************************************************************
 * Routine: Open_Input_File                                           *
 *                                                                    *
 * Opens the appropriate input data file.                             *
 *                                                                    *
 **********************************************************************/
void Open_Input_File(Electrostatic2D_Input_Parameters &IP) {

  IP.Input_File.open(IP.Input_File_Name,ios::in);
  if (!IP.Input_File.bad()) {
    IP.Line_Number = 0;
    IP.Input_File.setf(ios::skipws);
  }

}

/**********************************************************************
 * Routine: Close_Input_File                                          *
 *                                                                    *
 * Closes the appropriate input data file.                            *
 *                                                                    *
 **********************************************************************/
void Close_Input_File(Electrostatic2D_Input_Parameters &IP) {

  // Close input file.
  IP.Input_File.unsetf(ios::skipws);
  IP.Input_File.close();

}

/**********************************************************************
 * Routine: Set_Default_Input_Parameters                              *
 *                                                                    *
 * Assigns default values to the input parameters.                    *
 *                                                                    *
 **********************************************************************/
void Set_Default_Input_Parameters(Electrostatic2D_Input_Parameters &IP) {

  char *string_ptr;
  
  // CFFC root directory path:
  IP.get_cffc_path();

  string_ptr = "Electrostatic2D.in";
  strcpy(IP.Input_File_Name,string_ptr);

  // Time-stepping parameters:
  string_ptr = "Explicit_Euler";
  strcpy(IP.Time_Integration_Type,string_ptr);
  IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
  IP.Time_Accurate = 0;
  IP.Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
  IP.Maximum_Number_of_Time_Steps = 100;
  IP.N_Stage = 1;
  IP.CFL_Number = 0.5;
  IP.Time_Max = ZERO;

  // Residual smoothing:
  IP.Residual_Smoothing = 0;
  IP.Residual_Smoothing_Epsilon = ZERO;
  IP.Residual_Smoothing_Gauss_Seidel_Iterations = 2;

  // Reconstruction type:
  string_ptr = "Linear_Least_Squares";
  strcpy(IP.Reconstruction_Type,string_ptr);
  IP.i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;

  // Flux Reconstruction type:
  string_ptr = "Diamond_Path";
  strcpy(IP.Flux_Reconstruction_Type,string_ptr);
  IP.i_Flux_Reconstruction = VISCOUS_RECONSTRUCTION_HYBRID;

  // Initial conditions:
  string_ptr = "Uniform";
  strcpy(IP.ICs_Type,string_ptr);
  IP.i_ICs = IC_UNIFORM;

  // Geometry switch:
  string_ptr = "Planar";
  strcpy(IP.Flow_Geometry_Type,string_ptr);
  IP.Axisymmetric = OFF;

  // Electric field conditions:
  IP.Potential = ZERO;
  IP.Electric_Field_Strength = ZERO;
  IP.Electric_Field_Angle = ZERO;

  // Grid parameters:
  string_ptr = "Square";
  strcpy(IP.Grid_Type,string_ptr);
  IP.i_Grid = GRID_SQUARE;
  IP.Box_Width = ONE;
  IP.Box_Height = ONE;
  IP.Number_of_Cells_Idir = 100;
  IP.Number_of_Cells_Jdir = 100;
  IP.Number_of_Ghost_Cells = 2;
  IP.Number_of_Blocks_Idir = 1;
  IP.Number_of_Blocks_Jdir = 1;
  IP.Cylinder_Radius = ONE;
  IP.Orifice_Radius = ONE;
  IP.Wedge_Angle = 25.0;
  IP.Wedge_Length = HALF;
  IP.Chamber_Length = 0.835;
  IP.Chamber_Radius = 0.020;
  IP.Chamber_To_Throat_Length = 0.05;
  IP.Nozzle_Length = 0.150;
  IP.Nozzle_Radius_Exit = 0.030;
  IP.Nozzle_Radius_Throat = 0.010;
  IP.Nozzle_Type = NOZZLE_GOTTLIEB_FUNCTION;
  IP.X_Shift = Vector2D_ZERO;
  IP.X_Scale = ONE;
  IP.X_Rotate = ZERO;

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

  // ICEM:
  IP.ICEMCFD_FileNames = ICEMCFD_get_filenames();

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
  IP.Number_of_Refinement_Criteria = 3;
  IP.Refinement_Criteria_Gradient_Potential_Field = ON;
  IP.Refinement_Criteria_Divergence_Electric_Field = ON;
  IP.Refinement_Criteria_Curl_Electric_Field = ON;

  // Mesh stretching factor:
  IP.i_Mesh_Stretching = OFF;
  IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_LINEAR;
  IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_LINEAR;
  IP.Mesh_Stretching_Factor_Idir = 1.01;
  IP.Mesh_Stretching_Factor_Jdir = 1.01;

  // Smooth quad block indicator:
  IP.i_Smooth_Quad_Block = ON;

  // Default output file names and parameters:
  string_ptr = "outputfile.dat";
  strcpy(IP.Output_File_Name,string_ptr);
  string_ptr = "gridfile.grid";
  strcpy(IP.Grid_File_Name,string_ptr);
  string_ptr = "gridfile.griddef";
  strcpy(IP.Grid_Definition_File_Name,string_ptr);
  string_ptr = "restartfile.soln";
  strcpy(IP.Restart_File_Name,string_ptr);
  string_ptr = "gnuplotfile.gplt";
  strcpy(IP.Gnuplot_File_Name,string_ptr);
  string_ptr = "Tecplot";
  strcpy(IP.Output_Format_Type,string_ptr);
  IP.i_Output_Format = IO_TECPLOT;
  IP.Restart_Solution_Save_Frequency = 1000;

  // Default output progress frequency:
  IP.Output_Progress_Frequency = 50;

  // Input_file parameters:
  string_ptr = " ";
  strcpy(IP.Next_Control_Parameter,string_ptr);
  IP.Line_Number = 0;
  IP.Number_of_Processors = CFFC_MPI::Number_of_Processors;
  IP.Number_of_Blocks_Per_Processor = 10;

}

/**********************************************************************
 * Routine: Broadcast_Input_Parameters                                *
 *                                                                    *
 * Broadcast the input parameters variables to all processors in the  *
 * calculation from the primary processor using the MPI broadcast     *
 * routine.                                                           *
 *                                                                    *
 **********************************************************************/
void Broadcast_Input_Parameters(Electrostatic2D_Input_Parameters &IP) {

#ifdef _MPI_VERSION

  // CFFC path:
  MPI::COMM_WORLD.Bcast(IP.CFFC_Path, 
 			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D, 
			MPI::CHAR, 0);
  // Input file name and line number:
  MPI::COMM_WORLD.Bcast(IP.Input_File_Name,
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.Line_Number),
			1,
			MPI::INT,0);
  // Flow geometry type:
  MPI::COMM_WORLD.Bcast(IP.Flow_Geometry_Type,
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.Axisymmetric),
			1,
			MPI::INT,0);
  // Electrostatic force:
  MPI::COMM_WORLD.Bcast(&(IP.Potential),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Electric_Field_Strength),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Electric_Field_Angle),
			1,
			MPI::DOUBLE,0);
  // Time integration:
  MPI::COMM_WORLD.Bcast(IP.Time_Integration_Type,
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Time_Integration),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Time_Accurate),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Local_Time_Stepping),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Maximum_Number_of_Time_Steps),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.N_Stage),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.CFL_Number),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Time_Max),
			1,
			MPI::DOUBLE,0);
  // Residual smoothing:
  MPI::COMM_WORLD.Bcast(&(IP.Residual_Smoothing),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Residual_Smoothing_Epsilon),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Residual_Smoothing_Gauss_Seidel_Iterations),
			1,
			MPI::INT,0);
  // Reconstruction:
  MPI::COMM_WORLD.Bcast(IP.Reconstruction_Type,
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Reconstruction),
			1,
			MPI::INT,0);
  // Flux Reconstruction:
  MPI::COMM_WORLD.Bcast(IP.Flux_Reconstruction_Type,
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Flux_Reconstruction),
			1,
			MPI::INT,0);
  // Initial conditions:
  MPI::COMM_WORLD.Bcast(IP.ICs_Type,
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_ICs),
			1,
			MPI::INT,0);
  if (!CFFC_Primary_MPI_Processor()) {
    Initialize_Reference_State(IP);
  }
  for (int nv = 1; nv <= NUM_VAR_ELECTROSTATIC2D; nv++) {
    MPI::COMM_WORLD.Bcast(&(IP.U1[nv]),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.U2[nv]),
			  1,
			  MPI::DOUBLE,0);
  }
  // Grid variables:
  MPI::COMM_WORLD.Bcast(IP.Grid_Type,
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Grid),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Cells_Idir),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Cells_Jdir),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Ghost_Cells),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Idir),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Jdir),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Box_Width),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Box_Height),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Cylinder_Radius),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Orifice_Radius),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Wedge_Angle),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Wedge_Length),
			1,
			MPI::DOUBLE,0);
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
  MPI::COMM_WORLD.Bcast(&(IP.X_Shift.x),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.X_Shift.y),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.X_Scale),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.X_Rotate),
			1,
			MPI::DOUBLE,0);
  // Boundary Conditions:
  MPI::COMM_WORLD.Bcast(IP.Boundary_Conditions_Specified,
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.BCs_Specified),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(IP.BC_North_Type,
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.BC_South_Type,
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.BC_East_Type,
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.BC_West_Type,
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
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
  // ICEM:
  if (!CFFC_Primary_MPI_Processor()) {
    IP.ICEMCFD_FileNames = new char*[3];
    for (int i = 0; i < 3; i++) {
      IP.ICEMCFD_FileNames[i] = new char[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
    }
  }
  MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[0],
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[1],
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[2],
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
			MPI::CHAR,0);
  // AMR & Refinement Parameters:
  MPI::COMM_WORLD.Bcast(&(IP.AMR),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.AMR_Frequency),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Initial_Mesh_Refinements),
			1,
			MPI::INT,0);
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
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Threshold_for_Coarsening),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Refinement_Criteria),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Refinement_Criteria_Gradient_Potential_Field),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Refinement_Criteria_Divergence_Electric_Field),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Refinement_Criteria_Curl_Electric_Field),
			1,
			MPI::DOUBLE,0);
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
  // File Names:
  MPI::COMM_WORLD.Bcast(IP.Output_File_Name,
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.Grid_File_Name,
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.Grid_Definition_File_Name,
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.Restart_File_Name,
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.Gnuplot_File_Name,
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.Output_Format_Type,
			INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Output_Format),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Restart_Solution_Save_Frequency),
			1,
			MPI::INT,0);
  // Output progress frequency:
  MPI::COMM_WORLD.Bcast(&(IP.Output_Progress_Frequency),
			1,
			MPI::INT,0);
  // Number of processors:
  if (!CFFC_Primary_MPI_Processor()) {
    IP.Number_of_Processors = CFFC_MPI::Number_of_Processors;
  }
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Per_Processor),
			1,
			MPI::INT,0);

  // Reinitialize the reference state.
  Reinitialize_Reference_State(IP);

#endif

}

#ifdef _MPI_VERSION
/**********************************************************************
 * Routine: Broadcast_Input_Parameters                                *
 *                                                                    *
 * Broadcast the input parameters variables to all processors         *
 * associated with the specified communicator from the specified      *
 * processor using the MPI broadcast routine.                         *
 *                                                                    *
 **********************************************************************/
void Broadcast_Input_Parameters(Electrostatic2D_Input_Parameters &IP,
                                MPI::Intracomm &Communicator,
                                const int Source_CPU) {

  int Source_Rank = 0;

  // CFFC path:
  Communicator.Bcast(IP.CFFC_Path, 
 		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D, 
		     MPI::CHAR, Source_Rank);
  Communicator.Bcast(IP.Input_File_Name,
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.Line_Number),
		     1,
		     MPI::INT,Source_Rank);
  // Flow geometry type:
  Communicator.Bcast(IP.Flow_Geometry_Type,
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.Axisymmetric),
		     1,
		     MPI::INT,Source_Rank);
  // Electrostatic force:
  Communicator.Bcast(&(IP.Potential),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Electric_Field_Strength),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Electric_Field_Angle),
		     1,
		     MPI::DOUBLE,Source_Rank);
  // Time integration:
  Communicator.Bcast(IP.Time_Integration_Type,
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Time_Integration),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Time_Accurate),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Local_Time_Stepping),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Maximum_Number_of_Time_Steps),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.N_Stage),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.CFL_Number),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Time_Max),
		     1,
		     MPI::DOUBLE,Source_Rank);
  // Residual Smoothing:
  Communicator.Bcast(&(IP.Residual_Smoothing),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Residual_Smoothing_Epsilon),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Residual_Smoothing_Gauss_Seidel_Iterations),
		     1,
		     MPI::INT,Source_Rank);
  // Reconstruction:
  Communicator.Bcast(IP.Reconstruction_Type,
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Reconstruction),
		     1,
		     MPI::INT,Source_Rank);
  // Flux Reconstruction:
  Communicator.Bcast(IP.Flux_Reconstruction_Type,
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Flux_Reconstruction),
		     1,
		     MPI::INT,Source_Rank);
  // Initial conditions:
  Communicator.Bcast(IP.ICs_Type,
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_ICs),
		     1,
		     MPI::INT,Source_Rank);
  if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
    Initialize_Reference_State(IP);
  }
  for (int nv = 1; nv <= NUM_VAR_ELECTROSTATIC2D; nv++) {
    Communicator.Bcast(&(IP.U1[nv]),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.U2[nv]),
		       1,
		       MPI::DOUBLE,Source_Rank);
  }
  // Grid variables:
  Communicator.Bcast(IP.Grid_Type,
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Grid),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Cells_Idir),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Cells_Jdir),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Ghost_Cells),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Blocks_Idir),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Blocks_Jdir),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Box_Width),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Box_Height),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Cylinder_Radius),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Orifice_Radius),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Wedge_Angle),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Wedge_Length),
		     1,
		     MPI::DOUBLE,Source_Rank);
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
  Communicator.Bcast(&(IP.X_Shift.x),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.X_Shift.y),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.X_Scale),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.X_Rotate),
		     1,
		     MPI::DOUBLE,Source_Rank);
  // Boundary Conditions:
  Communicator.Bcast(IP.Boundary_Conditions_Specified,
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.BCs_Specified),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(IP.BC_North_Type,
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.BC_South_Type,
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.BC_East_Type,
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.BC_West_Type,
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
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
  // ICEM:
  if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
    IP.ICEMCFD_FileNames = new char*[3];
    for (int i = 0; i < 3; i++) {
      IP.ICEMCFD_FileNames[i] = new char[INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D];
    }
  }
  Communicator.Bcast(IP.ICEMCFD_FileNames[0],
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.ICEMCFD_FileNames[1],
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
                       MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.ICEMCFD_FileNames[2],
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
		     MPI::CHAR,Source_Rank);
  // AMR & Refinement Parameters:
  Communicator.Bcast(&(IP.AMR),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.AMR_Frequency),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Initial_Mesh_Refinements),
		     1,
		     MPI::INT,Source_Rank);
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
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Threshold_for_Coarsening),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Refinement_Criteria),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Refinement_Criteria_Gradient_Potential_Field),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Refinement_Criteria_Divergence_Electric_Field),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Refinement_Criteria_Curl_Electric_Field),
		     1,
		     MPI::DOUBLE,Source_Rank);
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
  // File Names:
  Communicator.Bcast(IP.Output_File_Name,
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.Grid_File_Name,
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.Grid_Definition_File_Name,
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.Restart_File_Name,
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.Gnuplot_File_Name,
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.Output_Format_Type,
		     INPUT_PARAMETER_LENGTH_ELECTROSTATIC2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Output_Format),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Restart_Solution_Save_Frequency),
		     1,
		     MPI::INT,Source_Rank);
  // Number of blocks per processor:
  if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
    IP.Number_of_Processors = CFFC_MPI::Number_of_Processors;
  }
  Communicator.Bcast(&(IP.Number_of_Blocks_Per_Processor),
		     1,
		     MPI::INT,Source_Rank);
  // Output progress frequency:
  Communicator.Bcast(&(IP.Output_Progress_Frequency),
		     1,
		     MPI::INT,Source_Rank);

  // Reinitialize the reference state.
  Reinitialize_Reference_State(IP);

}
#endif

/**********************************************************************
 * Routine: Get_Next_Input_Control_Parameter                          *
 *                                                                    *
 * Get the next input control parameter from the input file.          *
 *                                                                    *
 **********************************************************************/
void Get_Next_Input_Control_Parameter(Electrostatic2D_Input_Parameters &IP) {

  int i;
  char buffer[256];

  IP.Line_Number++;
  IP.Input_File.getline(buffer,sizeof(buffer));
  i = 0;
  if (buffer[0] != '#') {
    while (1) {
      if (buffer[i] == ' ' || buffer[i] == '=') break;
      i = i + 1;
      if (i > strlen(buffer)) break;
    }
    buffer[i] = '\0';
  }
  strcpy(IP.Next_Control_Parameter,buffer);

}

/**********************************************************************
 * Routine: Parse_Next_Input_Control_Parameter                        *
 *                                                                    *
 * Parses and executes the next input control parameter from the      *
 * input file.                                                        *
 *                                                                    *
 **********************************************************************/
int Parse_Next_Input_Control_Parameter(Electrostatic2D_Input_Parameters &IP) {

  int i_command = 0;
  char buffer[256];
  int tpt, bct;
  Vector2D Xt;

  if (strcmp(IP.Next_Control_Parameter, "CFFC_Path") == 0) {
     i_command = 1111;
     Get_Next_Input_Control_Parameter(IP);
     strcpy(IP.CFFC_Path, IP.Next_Control_Parameter);

  } else if (strcmp(IP.Next_Control_Parameter,"Time_Integration_Type") == 0) {
    i_command = 1;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Time_Integration_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Time_Integration_Type,"Explicit_Euler") == 0) {
      IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
      IP.N_Stage = 1;
    } else if (strcmp(IP.Time_Integration_Type,"Explicit_Predictor_Corrector") == 0) {
      IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR;
      IP.N_Stage = 2;
    } else if (strcmp(IP.Time_Integration_Type,"Explicit_Runge_Kutta") == 0) {
      IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_RUNGE_KUTTA;
      IP.N_Stage = 4;
    } else if (strcmp(IP.Time_Integration_Type,"Multistage_Optimal_Smoothing") == 0) {
      IP.i_Time_Integration = TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING;
      IP.N_Stage = 3;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Reconstruction_Type") == 0) {
    i_command = 2;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Reconstruction_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Reconstruction_Type,"Green_Gauss") == 0) {
      IP.i_Reconstruction = RECONSTRUCTION_GREEN_GAUSS;
    } else if (strcmp(IP.Reconstruction_Type,"Least_Squares") == 0 ||
	       strcmp(IP.Reconstruction_Type,"Linear_Least_Squares") == 0) {
      IP.i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;
    } else if (strcmp(IP.Reconstruction_Type,"Quadratic_Least_Squares") == 0) {
      i_command = INVALID_INPUT_VALUE;
      IP.i_Reconstruction = RECONSTRUCTION_QUADRATIC_LEAST_SQUARES;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Flux_Reconstruction_Type") == 0) {
    i_command = 3;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Flux_Reconstruction_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Flux_Reconstruction_Type,"Diamond_Path") == 0) {
      IP.i_Flux_Reconstruction = VISCOUS_RECONSTRUCTION_DIAMOND_PATH;
    } else if (strcmp(IP.Flux_Reconstruction_Type,"Arithmetic_Average") == 0) {
      IP.i_Flux_Reconstruction = VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE;
    } else if (strcmp(IP.Flux_Reconstruction_Type,"Hybrid") == 0) {
      IP.i_Flux_Reconstruction = VISCOUS_RECONSTRUCTION_HYBRID;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"ICs_Type") == 0) {
    i_command = 5;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.ICs_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.ICs_Type,"Constant") == 0) {
      IP.i_ICs = IC_CONSTANT;
    } else if (strcmp(IP.ICs_Type,"Uniform") == 0) {
      IP.i_ICs = IC_UNIFORM;
    } else if (strcmp(IP.ICs_Type,"Electrostatic_Channel") == 0) {
      IP.i_ICs = IC_ELECTROSTATIC_CHANNEL;
    } else if (strcmp(IP.ICs_Type,"Desolvation_Chamber") == 0) {
      IP.i_ICs = IC_DESOLVATION_CHAMBER;
    } else if (strcmp(IP.ICs_Type,"Restart") == 0) {
      IP.i_ICs = IC_RESTART;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Grid_Type") == 0) {
    i_command = 6;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Grid_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Grid_Type,"Cartesian") == 0) {
      IP.i_Grid = GRID_CARTESIAN_UNIFORM;
      IP.Box_Width = ONE;
      IP.Box_Height = ONE;
    } else if (strcmp(IP.Grid_Type,"Square") == 0) {
      IP.i_Grid = GRID_SQUARE;
      IP.Box_Width = ONE;
      IP.Box_Height = ONE;
    } else if (strcmp(IP.Grid_Type,"Rectangular_Box") == 0) {
      IP.i_Grid = GRID_RECTANGULAR_BOX;
      IP.Box_Width = ONE;
      IP.Box_Height = ONE;
    } else if (strcmp(IP.Grid_Type,"Circular_Cylinder") == 0) {
      IP.i_Grid = GRID_CIRCULAR_CYLINDER;
      IP.Cylinder_Radius = ONE;
      IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
      IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
      IP.Mesh_Stretching_Factor_Idir = 1.025;
      IP.Mesh_Stretching_Factor_Jdir = 1.001;
    } else if (strcmp(IP.Grid_Type,"Desolvation_Chamber") == 0) {
      IP.i_Grid = GRID_DESOLVATION_CHAMBER;
    } else if (strcmp(IP.Grid_Type,"ICEMCFD") == 0) {
      IP.i_Grid = GRID_ICEMCFD;
      strcpy(IP.ICEMCFD_FileNames[0],"topo_mulcad_out");
      strcpy(IP.ICEMCFD_FileNames[1],"family_boco");
      strcpy(IP.ICEMCFD_FileNames[2],"family_topo");
    } else if (strcmp(IP.Grid_Type,"Read_From_Definition_File") == 0) {
      IP.i_Grid = GRID_READ_FROM_DEFINITION_FILE;
    } else if (strcmp(IP.Grid_Type,"Read_From_Data_File") == 0) {
      IP.i_Grid = GRID_READ_FROM_GRID_DATA_FILE;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Output_File_Name") == 0) {
    i_command = 7;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Output_File_Name,IP.Next_Control_Parameter);
    strcat(IP.Output_File_Name,".dat");
    strcpy(IP.Grid_File_Name,IP.Next_Control_Parameter);
    strcat(IP.Grid_File_Name,".grid");
    strcpy(IP.Grid_Definition_File_Name,IP.Next_Control_Parameter);
    strcat(IP.Grid_Definition_File_Name,".griddef");
    strcpy(IP.Restart_File_Name,IP.Next_Control_Parameter);
    strcat(IP.Restart_File_Name,".soln");
    strcpy(IP.Gnuplot_File_Name,IP.Next_Control_Parameter);
    strcat(IP.Gnuplot_File_Name,".gplt");

  } else if (strcmp(IP.Next_Control_Parameter,"Grid_File_Name") == 0) {
    i_command = 8;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Grid_File_Name,IP.Next_Control_Parameter);
    strcat(IP.Grid_File_Name,".grid");
    strcpy(IP.Grid_Definition_File_Name,IP.Next_Control_Parameter);
    strcat(IP.Grid_Definition_File_Name,".griddef");

  } else if (strcmp(IP.Next_Control_Parameter,"Restart_File_Name") == 0) {
    i_command = 9;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Restart_File_Name,IP.Next_Control_Parameter);
    strcat(IP.Restart_File_Name,".soln");

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Cells_Idir") == 0) {
    i_command = 10;
    IP.Line_Number++;
    IP.Input_File >> IP.Number_of_Cells_Idir;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Cells_Idir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Cells_Jdir") == 0) {
    i_command = 11;
    IP.Line_Number++;
    IP.Input_File >> IP.Number_of_Cells_Jdir;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Cells_Jdir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Ghost_Cells") == 0) {
    i_command = 12;
    IP.Line_Number++;
    IP.Input_File >> IP.Number_of_Ghost_Cells;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Ghost_Cells < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Blocks_Idir") == 0) {
    i_command = 12;
    IP.Line_Number++;
    IP.Input_File >> IP.Number_of_Blocks_Idir;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Blocks_Idir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Blocks_Jdir") == 0) {
    i_command = 13;
    IP.Line_Number++;
    IP.Input_File >> IP.Number_of_Blocks_Jdir;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Blocks_Jdir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Time_Accurate") == 0) {
    i_command = 14;
    IP.Line_Number++;
    IP.Input_File >> IP.Time_Accurate;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Time_Accurate != 0 && IP.Time_Accurate != 1) IP.Time_Accurate = 0;
    if (IP.Time_Accurate) IP.Local_Time_Stepping = GLOBAL_TIME_STEPPING;

  } else if (strcmp(IP.Next_Control_Parameter,"Local_Time_Stepping") == 0) {
    i_command = 15;
    IP.Line_Number++;
    IP.Input_File >> IP.Local_Time_Stepping;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Local_Time_Stepping != GLOBAL_TIME_STEPPING &&
	IP.Local_Time_Stepping != SCALAR_LOCAL_TIME_STEPPING) {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Maximum_Number_of_Time_Steps") == 0) {
    i_command = 16;
    IP.Line_Number++;
    IP.Input_File >> IP.Maximum_Number_of_Time_Steps;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Maximum_Number_of_Time_Steps < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"N_Stage") == 0) {
    i_command = 17;
    IP.Line_Number++;
    IP.Input_File >> IP.N_Stage;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.i_Time_Integration == TIME_STEPPING_EXPLICIT_EULER && IP.N_Stage > 1) i_command = INVALID_INPUT_VALUE;
    if (IP.N_Stage < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"CFL_Number") == 0) {
    i_command = 18;
    IP.Line_Number++;
    IP.Input_File >> IP.CFL_Number;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.CFL_Number <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Box_Width") == 0) {
    i_command = 19;
    IP.Line_Number++;
    IP.Input_File >> IP.Box_Width;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Box_Width <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Box_Height") == 0) {
    i_command = 20;
    IP.Line_Number++;
    IP.Input_File >> IP.Box_Height;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Box_Height <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Cylinder_Radius") == 0) {
    i_command = 26;
    IP.Line_Number++;
    IP.Input_File >> IP.Cylinder_Radius;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Cylinder_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Orifice_Radius") == 0) {
    i_command = 31;
    IP.Line_Number++;
    IP.Input_File >> IP.Orifice_Radius;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Orifice_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Wedge_Angle") == 0) {
    i_command = 31;
    IP.Line_Number++;
    IP.Input_File >> IP.Wedge_Angle;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Wedge_Angle <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Wedge_Length") == 0) {
    i_command = 31;
    IP.Line_Number++;
    IP.Input_File >> IP.Wedge_Length;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Wedge_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Time_Max") == 0) {
    i_command = 45;
    IP.Line_Number++;
    IP.Input_File >> IP.Time_Max;
    IP.Input_File.getline(buffer,sizeof(buffer));
    IP.Time_Max = IP.Time_Max/THOUSAND;
    if (IP.Time_Max < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Blocks_Per_Processor") == 0) {
    i_command = 46;
    IP.Line_Number++;
    IP.Input_File >> IP.Number_of_Blocks_Per_Processor;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Blocks_Per_Processor < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Output_Format_Type") == 0) {
    i_command = 47;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Output_Format_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Output_Format_Type,"Gnuplot") == 0  ||
	strcmp(IP.Output_Format_Type,"gnuplot") == 0  ||
	strcmp(IP.Output_Format_Type,"GNUPLOT") == 0) {
      IP.i_Output_Format = IO_GNUPLOT;
    } else if (strcmp(IP.Output_Format_Type,"Tecplot") == 0  ||
	       strcmp(IP.Output_Format_Type,"tecplot") == 0  ||
	       strcmp(IP.Output_Format_Type,"TECPLOT") == 0) {
      IP.i_Output_Format = IO_TECPLOT;
    } else if (strcmp(IP.Output_Format_Type,"Matlab") == 0  ||
	       strcmp(IP.Output_Format_Type,"matlab") == 0  ||
	       strcmp(IP.Output_Format_Type,"MATLAB") == 0) {
      IP.i_Output_Format = IO_MATLAB;
    } else if (strcmp(IP.Output_Format_Type,"Octave") == 0  ||
	       strcmp(IP.Output_Format_Type,"octave") == 0  ||
	       strcmp(IP.Output_Format_Type,"OCTAVE") == 0) {
      IP.i_Output_Format = IO_OCTAVE;
    } else {
      IP.i_Output_Format = IO_TECPLOT;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Flow_Geometry_Type") == 0) {
    i_command = 49;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Flow_Geometry_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Flow_Geometry_Type,"Planar") == 0) {
      IP.Axisymmetric = OFF;
    } else if (strcmp(IP.Flow_Geometry_Type,"Axisymmetric") == 0) {
      IP.Axisymmetric = ON;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Restart_Solution_Save_Frequency") == 0) {
    i_command = 49;
    IP.Line_Number++;
    IP.Input_File >> IP.Restart_Solution_Save_Frequency;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Restart_Solution_Save_Frequency < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"ICEMCFD_Topology_File") == 0) {
    i_command = 53;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.ICEMCFD_FileNames[0],IP.Next_Control_Parameter);

  } else if (strcmp(IP.Next_Control_Parameter,"ICEMCFD_Family_Boco_File") == 0) {
    i_command = 54;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.ICEMCFD_FileNames[1],IP.Next_Control_Parameter);

  } else if (strcmp(IP.Next_Control_Parameter,"ICEMCFD_Family_Topo_File") == 0) {
    i_command = 55;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.ICEMCFD_FileNames[2],IP.Next_Control_Parameter);

  } else if (strcmp(IP.Next_Control_Parameter,"X_Shift") == 0) {
    i_command = 56;
    IP.Line_Number++;
    IP.Input_File >> IP.X_Shift;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"X_Scale") == 0) {
    i_command = 57;
    IP.Line_Number++;
    IP.Input_File >> IP.X_Scale;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"X_Rotate") == 0) {
    i_command = 58;
    IP.Line_Number++;
    IP.Input_File >> IP.X_Rotate;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"Potential") == 0) {
    i_command = 61;
    IP.Line_Number++;
    IP.Input_File >> IP.Potential;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"Electric_Field_Strength") == 0) {
    i_command = 62;
    IP.Line_Number++;
    IP.Input_File >> IP.Electric_Field_Strength;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"Electric_Field_Angle") == 0) {
    i_command = 63;
    IP.Line_Number++;
    IP.Input_File >> IP.Electric_Field_Angle;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"AMR") == 0) {
    i_command = 71;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
      IP.AMR = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
      IP.AMR = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"AMR_Frequency") == 0) {
    i_command = 72;
    IP.Line_Number++;
    IP.Input_File >> IP.AMR_Frequency;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.AMR_Frequency < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Initial_Mesh_Refinements") == 0) {
    i_command = 73;
    IP.Line_Number++;
    IP.Input_File >> IP.Number_of_Initial_Mesh_Refinements;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Initial_Mesh_Refinements < 0) IP.Number_of_Initial_Mesh_Refinements = 0;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Uniform_Mesh_Refinements") == 0) {
    i_command = 74;
    IP.Line_Number++;
    IP.Input_File >> IP.Number_of_Uniform_Mesh_Refinements;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Uniform_Mesh_Refinements < 0) IP.Number_of_Uniform_Mesh_Refinements = 0;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Boundary_Mesh_Refinements") == 0) {
    i_command = 75;
    IP.Line_Number++;
    IP.Input_File >> IP.Number_of_Boundary_Mesh_Refinements;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Boundary_Mesh_Refinements < 0) IP.Number_of_Boundary_Mesh_Refinements = 0;

  } else if (strcmp(IP.Next_Control_Parameter,"Maximum_Refinement_Level") == 0) {
    i_command = 76;
    IP.Line_Number++;
    IP.Input_File >> IP.Maximum_Refinement_Level;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Maximum_Refinement_Level < 1) IP.Maximum_Refinement_Level = 1;

  } else if (strcmp(IP.Next_Control_Parameter,"Minimum_Refinement_Level") == 0) {
    i_command = 77;
    IP.Line_Number++;
    IP.Input_File >> IP.Minimum_Refinement_Level;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Minimum_Refinement_Level < 1) IP.Minimum_Refinement_Level = 1;

  } else if (strcmp(IP.Next_Control_Parameter,"Threshold_for_Refinement") == 0) {
    i_command = 78;
    IP.Line_Number++;
    IP.Input_File >> IP.Threshold_for_Refinement;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Threshold_for_Refinement <= ZERO ||
	IP.Threshold_for_Refinement > ONE) IP.Threshold_for_Refinement = 0.50;

  } else if (strcmp(IP.Next_Control_Parameter,"Threshold_for_Coarsening") == 0) {
    i_command = 79;
    IP.Line_Number++;
    IP.Input_File >> IP.Threshold_for_Coarsening;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Threshold_for_Coarsening < ZERO ||
	IP.Threshold_for_Coarsening >= ONE) IP.Threshold_for_Coarsening = 0.10;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Refinement_Criteria") == 0) {
    i_command = 80;
    IP.Line_Number++;
    IP.Input_File >> IP.Number_of_Refinement_Criteria;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Refinement_Criteria < 1 || IP.Number_of_Refinement_Criteria > 6) {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Gradient_Potential_Field") == 0) {
    i_command = 81;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
      IP.Refinement_Criteria_Gradient_Potential_Field = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
      IP.Refinement_Criteria_Gradient_Potential_Field = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Divergence_Electric_Field") == 0) {
    i_command = 82;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
      IP.Refinement_Criteria_Divergence_Electric_Field = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
      IP.Refinement_Criteria_Divergence_Electric_Field = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Curl_Electric_Field") == 0) {
    i_command = 83;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
      IP.Refinement_Criteria_Curl_Electric_Field = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
      IP.Refinement_Criteria_Curl_Electric_Field = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Residual_Smoothing_Epsilon") == 0) {
    i_command = 91;
    IP.Line_Number++;
    IP.Input_File >> IP.Residual_Smoothing_Epsilon;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Residual_Smoothing_Epsilon <= ZERO) {
      IP.Residual_Smoothing = OFF;
      IP.Residual_Smoothing_Epsilon = ZERO;
    } else {
      IP.Residual_Smoothing = ON;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Residual_Smoothing_Gauss_Seidel_Iterations") == 0) {
    i_command = 92;
    IP.Line_Number++;
    IP.Input_File >> IP.Residual_Smoothing_Gauss_Seidel_Iterations;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Residual_Smoothing_Gauss_Seidel_Iterations < 0)
      IP.Residual_Smoothing_Gauss_Seidel_Iterations = 0;

  } else if (strcmp(IP.Next_Control_Parameter,"Mesh_Stretching") == 0) {
    i_command = 101;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
      IP.i_Mesh_Stretching = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
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
    i_command = 103;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Mesh_Stretching_Factor_Idir;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Mesh_Stretching_Factor_Idir < ONE) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Mesh_Stretching_Factor_Jdir") == 0) {
    i_command = 104;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Mesh_Stretching_Factor_Jdir;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Mesh_Stretching_Factor_Jdir < ONE) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Smooth_Quad_Block") == 0) {
    i_command = 105;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
      IP.i_Smooth_Quad_Block = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
      IP.i_Smooth_Quad_Block = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Boundary_Conditions_Specified") == 0) {
    i_command = 500;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Boundary_Conditions_Specified,IP.Next_Control_Parameter);
    if (strcmp(IP.Boundary_Conditions_Specified,"ON") == 0) {
      IP.BCs_Specified = ON;
    } else if (strcmp(IP.Boundary_Conditions_Specified,"OFF") == 0) {
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
    } else if (strcmp(IP.BC_North_Type,"Wall_Viscous_Isothermal") == 0) {
      IP.BC_North = BC_WALL_VISCOUS_ISOTHERMAL;
    } else if (strcmp(IP.BC_North_Type,"Wall_Viscous_Heatflux") == 0) {
      IP.BC_North = BC_WALL_VISCOUS_HEATFLUX;
    } else if (strcmp(IP.BC_North_Type,"Fixed") == 0) {
      IP.BC_North = BC_FIXED;
    } else if (strcmp(IP.BC_North_Type,"Constant_Extrapolation") == 0) {
      IP.BC_North = BC_CONSTANT_EXTRAPOLATION;
    } else if (strcmp(IP.BC_North_Type,"Linear_Extrapolation") == 0) {
      IP.BC_North = BC_LINEAR_EXTRAPOLATION;
    } else if (strcmp(IP.BC_North_Type,"None") == 0) {
      IP.BC_North = BC_NONE;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"BC_South") == 0) {
    i_command = 502;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.BC_South_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.BC_South_Type,"Reflection") == 0) {
      IP.BC_South = BC_REFLECTION;
    } else if (strcmp(IP.BC_South_Type,"Wall_Viscous_Isothermal") == 0) {
      IP.BC_South = BC_WALL_VISCOUS_ISOTHERMAL;
    } else if (strcmp(IP.BC_South_Type,"Wall_Viscous_Heatflux") == 0) {
      IP.BC_South = BC_WALL_VISCOUS_HEATFLUX;
    } else if (strcmp(IP.BC_South_Type,"Fixed") == 0) {
      IP.BC_South = BC_FIXED;
    } else if (strcmp(IP.BC_South_Type,"Constant_Extrapolation") == 0) {
      IP.BC_South = BC_CONSTANT_EXTRAPOLATION;
    } else if (strcmp(IP.BC_South_Type,"Linear_Extrapolation") == 0) {
      IP.BC_South = BC_LINEAR_EXTRAPOLATION;
    } else if (strcmp(IP.BC_South_Type,"None") == 0) {
      IP.BC_South = BC_NONE;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"BC_East") == 0) {
    i_command = 503;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.BC_East_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.BC_East_Type,"Reflection") == 0) {
      IP.BC_East = BC_REFLECTION;
    } else if (strcmp(IP.BC_East_Type,"Wall_Viscous_Isothermal") == 0) {
      IP.BC_East = BC_WALL_VISCOUS_ISOTHERMAL;
    } else if (strcmp(IP.BC_East_Type,"Wall_Viscous_Heatflux") == 0) {
      IP.BC_East = BC_WALL_VISCOUS_HEATFLUX;
    } else if (strcmp(IP.BC_East_Type,"Fixed") == 0) {
      IP.BC_East = BC_FIXED;
    } else if (strcmp(IP.BC_East_Type,"Constant_Extrapolation") == 0) {
      IP.BC_East = BC_CONSTANT_EXTRAPOLATION;
    } else if (strcmp(IP.BC_East_Type,"Linear_Extrapolation") == 0) {
      IP.BC_East = BC_LINEAR_EXTRAPOLATION;
    } else if (strcmp(IP.BC_East_Type,"None") == 0) {
      IP.BC_East = BC_NONE;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"BC_West") == 0) {
    i_command = 504;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.BC_West_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.BC_West_Type,"Reflection") == 0) {
      IP.BC_West = BC_REFLECTION;
    } else if (strcmp(IP.BC_West_Type,"Wall_Viscous_Isothermal") == 0) {
      IP.BC_West = BC_WALL_VISCOUS_ISOTHERMAL;
    } else if (strcmp(IP.BC_West_Type,"Wall_Viscous_Heatflux") == 0) {
      IP.BC_West = BC_WALL_VISCOUS_HEATFLUX;
    } else if (strcmp(IP.BC_West_Type,"Fixed") == 0) {
      IP.BC_West = BC_FIXED;
    } else if (strcmp(IP.BC_West_Type,"Constant_Extrapolation") == 0) {
      IP.BC_West = BC_CONSTANT_EXTRAPOLATION;
    } else if (strcmp(IP.BC_West_Type,"Linear_Extrapolation") == 0) {
      IP.BC_West = BC_LINEAR_EXTRAPOLATION;
    } else if (strcmp(IP.BC_West_Type,"None") == 0) {
      IP.BC_West = BC_NONE;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Execute") == 0) {
    i_command = EXECUTE_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Terminate") == 0) {
    i_command = TERMINATE_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Continue") == 0) {
    i_command = CONTINUE_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output") == 0) {
    i_command = WRITE_OUTPUT_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Cells") == 0) {
    i_command = WRITE_OUTPUT_CELLS_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Quasi3D") == 0) {
    i_command = WRITE_OUTPUT_QUASI3D_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Restart") == 0) {
    i_command = WRITE_RESTART_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Mesh") == 0) {
    i_command = WRITE_OUTPUT_GRID_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Mesh_Definition") == 0) {
    i_command = WRITE_GRID_DEFINITION_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Mesh_Nodes") == 0) {
    i_command = WRITE_OUTPUT_GRID_NODES_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Mesh_Cells") == 0) {
    i_command = WRITE_OUTPUT_GRID_CELLS_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Refine_Grid") == 0) {
    i_command = REFINE_GRID_CODE;

  } else if (IP.Next_Control_Parameter[0] == '#') {
    i_command = COMMENT_CODE;

  } else {
    i_command = INVALID_INPUT_CODE;

  }

  // Return the parser command type indicator.
  return i_command;

}

/**********************************************************************
 * Routine: Process_Input_Control_Parameter_File                      *
 *                                                                    *
 * Reads, parses, and executes the list of input control parameters   *
 * from the standard input file.                                      *
 *                                                                    *
 **********************************************************************/
int Process_Input_Control_Parameter_File(Electrostatic2D_Input_Parameters &IP,
                                         char *Input_File_Name_ptr,
                                         int &Command_Flag) {

  int error_flag, line_number;

  // Assign default values to the input parameters.
  Set_Default_Input_Parameters(IP);

  // Copy input file name (a string) to appropriate input parameter variable.
  if (Input_File_Name_ptr != NULL) strcpy(IP.Input_File_Name,Input_File_Name_ptr);

  // Open the input file containing the input parameters.
  Open_Input_File(IP);
  error_flag = IP.Input_File.bad();

  if (error_flag) {
    cout << "\n Electrostatic2D ERROR: Unable to open Electrostatic2D input data file.\n";
    return error_flag;
  }

  // Read and parse control parameters contained in the input file.
  while (1) {
    Get_Next_Input_Control_Parameter(IP);
    Command_Flag = Parse_Next_Input_Control_Parameter(IP);
    line_number = IP.Line_Number;
    if (Command_Flag == EXECUTE_CODE) {
      break;
    } else if (Command_Flag == TERMINATE_CODE) {
      break;
    } else if (Command_Flag == INVALID_INPUT_CODE ||
	       Command_Flag == INVALID_INPUT_VALUE) {
      line_number = -line_number;
      cout << "\n Electrostatic2D ERROR: Error reading Electrostatic2D data at line #"
	   << -line_number  << " of input data file.\n";
      error_flag = line_number;
      break;
    }
  }
  if (error_flag) return error_flag;

  // Set static variables and initialize reference state.
  Initialize_Reference_State(IP);

  // Perform consitency checks on the refinement criteria.
  IP.Number_of_Refinement_Criteria = 0;
  if (IP.Refinement_Criteria_Gradient_Potential_Field) IP.Number_of_Refinement_Criteria++;
  if (IP.Refinement_Criteria_Divergence_Electric_Field) IP.Number_of_Refinement_Criteria++;
  if (IP.Refinement_Criteria_Curl_Electric_Field) IP.Number_of_Refinement_Criteria++;
  if (IP.Number_of_Refinement_Criteria < 1 || IP.Number_of_Refinement_Criteria > 6) return 1011;

  // Initial processing of input control parameters complete.  
  return 0;

}

/**********************************************************************
 * Routine: Initialize_Reference_State                                *
 *                                                                    *
 * This function sets the static variables and reference states.      *
 *                                                                    *
 **********************************************************************/
void Initialize_Reference_State(Electrostatic2D_Input_Parameters &IP) {

  double cos_angle, sin_angle;

  IP.Uo.V = IP.Potential;
  cos_angle = cos(TWO*PI*IP.Electric_Field_Angle/360.00);
  sin_angle = sin(TWO*PI*IP.Electric_Field_Angle/360.00);
  IP.Uo.E.x = IP.Electric_Field_Strength*cos_angle;
  IP.Uo.E.y = IP.Electric_Field_Strength*sin_angle;

}

/**********************************************************************
 * Routine: Reinitialize_Reference_State                              *
 *                                                                    *
 * This function sets the static variables and reference states.      *
 *                                                                    *
 **********************************************************************/
void Reinitialize_Reference_State(Electrostatic2D_Input_Parameters &IP) {

  double cos_angle, sin_angle;

  IP.Uo.V = IP.Potential;
  cos_angle = cos(TWO*PI*IP.Electric_Field_Angle/360.00);
  sin_angle = sin(TWO*PI*IP.Electric_Field_Angle/360.00);
  IP.Uo.E.x = IP.Electric_Field_Strength*cos_angle;
  IP.Uo.E.y = IP.Electric_Field_Strength*sin_angle;

}
