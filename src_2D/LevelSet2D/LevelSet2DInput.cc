/**********************************************************************
 * LevelSet2DInput.cc                                                 *
 *                                                                    *
 * Subroutines for 2D Level Set input classes.                        *
 *                                                                    *
 **********************************************************************/

// Include 2D LevelSet input parameter header file.
#ifndef _LEVELSET2D_INPUT_INCLUDED
#include "LevelSet2DInput.h"
#endif // _LEVELSET2D_INPUT_INCLUDED

/**********************************************************************
 * LevelSet2D_Input_Parameters -- External subroutines.               *
 **********************************************************************/

/**********************************************************************
 * Routine: Open_Input_File                                           *
 *                                                                    *
 * Opens the appropriate input data file.                             *
 *                                                                    *
 **********************************************************************/
void Open_Input_File(LevelSet2D_Input_Parameters &IP) {

  IP.Input_File.open(IP.Input_File_Name,ios::in);
  if (!IP.Input_File.fail()) {
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
void Close_Input_File(LevelSet2D_Input_Parameters &IP) {

  IP.Input_File.unsetf(ios::skipws);
  IP.Input_File.close();

}

/**********************************************************************
 * Routine: Set_Default_Input_Parameters                              *
 *                                                                    *
 * Assigns default values to the input parameters.                    *
 *                                                                    *
 **********************************************************************/
void Set_Default_Input_Parameters(LevelSet2D_Input_Parameters &IP) {

  char *string_ptr;

  // CFFC root directory path:
  IP.get_cffc_path();

  string_ptr = "LevelSet2D.in";
  strcpy(IP.Input_File_Name,string_ptr);

  // Hamilton-Jacobi time-stepping parameters:
  string_ptr = "Explicit_Predictor_Corrector";
  strcpy(IP.Time_Integration_Type,string_ptr);
  IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR;
  IP.Time_Accurate = 1;
  IP.Local_Time_Stepping = 0;
  IP.Maximum_Number_of_Time_Steps = 100;
  IP.N_Stage = 2;
  IP.Time_Max = ZERO;

  // Residual variable:
  IP.i_Residual_Variable = 1;

  // Initial distance function parameters:
  string_ptr = "Geometric";
  strcpy(IP.Initial_Distance_Type,string_ptr);
  IP.i_Initial_Distance_Type = LEVELSET_INITIAL_EXTENSION_GEOMETRIC;
  IP.Perturb_Distance_Function = OFF;
  IP.Extension_Distance = MILLION;

  // Eikonal equation parameters:
  string_ptr = "Frequency";
  strcpy(IP.Redistance_Criteria,string_ptr);
  IP.i_Redistance_Criteria = EIKONAL_CRITERIA_FREQUENCY;
  IP.Redistance_Frequency = 1;
  IP.Redistance_Tolerance = 0.001;
  IP.Number_of_Initial_Redistance_Iterations = 200;
  IP.Number_of_Redistance_Iterations = 50;
  string_ptr = "Sussman";
  strcpy(IP.Eikonal_Scheme,string_ptr);
  IP.i_Eikonal_Scheme = EIKONAL_SCHEME_SUSSMAN;
  string_ptr = "Derivative";
  strcpy(IP.Eikonal_Sign_Function,string_ptr);
  IP.i_Eikonal_Sign_Function = EIKONAL_SIGN_FUNCTION_DERIVATIVE;
  IP.Eikonal_Threshold = HALF;

  // Scalar extension equation parameters:
  IP.Scalar_Extension_CFL_Number = HALF;
  IP.Number_of_Scalar_Extension_Iterations = 10;

  // Reconstruction type:
  string_ptr = "Linear_Essentially_Non_Oscillatory";
  strcpy(IP.Reconstruction_Type,string_ptr);
  IP.i_Reconstruction = RECONSTRUCTION_LINEAR_ESSENTIALLY_NON_OSCILLATORY;

  // Limiter type:
  string_ptr = "Unlimited";
  strcpy(IP.Limiter_Type,string_ptr);
  IP.i_Limiter = LIMITER_UNLIMITED;

  // Boundary condition type:
  string_ptr = "Constant_Extrapolation";
  strcpy(IP.BC_Type,string_ptr);
  IP.i_BC_Type = BC_CONSTANT_EXTRAPOLATION;

  // Curvature driven flow default value.
  IP.Curvature_Speed = ZERO;
  string_ptr = "Laplacian";
  strcpy(IP.Curvature_Scheme,string_ptr);
  IP.i_Curvature_Scheme = CURVATURE_SCHEME_LAPLACIAN;

  // Bulk flow field default values.
  string_ptr = "None";
  strcpy(IP.BulkFlowField_Type,string_ptr);
  IP.i_BulkFlowField_Type = INTERFACE_BULKFLOWFIELD_NONE;
  IP.V = Vector2D_ZERO;

  // Grid parameters:
  string_ptr = "Square";
  strcpy(IP.Grid_Type,string_ptr);
  IP.i_Grid = GRID_SQUARE;
  IP.Box_Width = ONE;
  IP.Box_Height = ONE;
  IP.Grain_Length = 0.835;
  IP.Grain_Radius = 0.020;
  IP.Grain_To_Throat_Length = 0.05;
  IP.Nozzle_Length = 0.150;
  IP.Nozzle_Radius_Exit = 0.030;
  IP.Nozzle_Radius_Throat = 0.010;
  IP.Number_of_Cells_Idir = 100;
  IP.Number_of_Cells_Jdir = 100;
  IP.Number_of_Ghost_Cells = 2;
  IP.Number_of_Blocks_Idir = 1;
  IP.Number_of_Blocks_Jdir = 1;
  IP.X_Shift = Vector2D_ZERO;
  IP.X_Scale = ONE;
  IP.X_Rotate = ZERO;

  // Mesh stretching factor:
  IP.i_Mesh_Stretching = OFF;
  IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_LINEAR;
  IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_LINEAR;
  IP.Mesh_Stretching_Factor_Idir = 1.01;
  IP.Mesh_Stretching_Factor_Jdir = 1.01;

  // AMR:
  IP.AMR = 0;
  IP.AMR_Frequency = 100;
  IP.Number_of_Initial_Mesh_Refinements = 0;
  IP.Number_of_Uniform_Mesh_Refinements = 0;
  IP.Maximum_Refinement_Level = 100;
  IP.Minimum_Refinement_Level = 1;
  IP.Threshold_for_Refinement = 0.50;
  IP.Threshold_for_Coarsening = 0.10;
  IP.Number_of_Refinement_Criteria = 1;
  IP.Refinement_Criteria_Curvature = ON;
  IP.Refinement_Criteria_Zero_Level_Set = OFF;

  // Smooth quad block flag:
  IP.i_Smooth_Quad_Block = OFF;

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

  string_ptr = " ";
  strcpy(IP.Next_Control_Parameter,string_ptr);

  IP.Line_Number = 0;
  
  IP.Number_of_Processors = CFFC_MPI::Number_of_Processors;
  IP.Number_of_Blocks_Per_Processor = 10;

}

/**********************************************************************
 * Routine: Broadcast_Input_Parameters                                *
 *                                                                    *
 * Broadcast the input parameters variables to all processors         *
 * involved in the calculation from the primary processor using the   *
 * MPI broadcast routine.                                             *
 *                                                                    *
 **********************************************************************/
void Broadcast_Input_Parameters(LevelSet2D_Input_Parameters &IP) {

#ifdef _MPI_VERSION

  // CFFC path:
  MPI::COMM_WORLD.Bcast(IP.CFFC_Path, 
 			INPUT_PARAMETER_LENGTH_LEVELSET2D, 
			MPI::CHAR, 0);
  MPI::COMM_WORLD.Bcast(IP.Input_File_Name,
			INPUT_PARAMETER_LENGTH_LEVELSET2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.Line_Number),
                        1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(IP.Time_Integration_Type,
			INPUT_PARAMETER_LENGTH_LEVELSET2D,
			MPI::CHAR,0);
  // Time integration:
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
  MPI::COMM_WORLD.Bcast(&(IP.Hamilton_Jacobi_CFL_Number),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Time_Max),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(IP.Reconstruction_Type,
			INPUT_PARAMETER_LENGTH_LEVELSET2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Reconstruction),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(IP.Limiter_Type,
			INPUT_PARAMETER_LENGTH_LEVELSET2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Limiter),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(IP.BC_Type,
			INPUT_PARAMETER_LENGTH_LEVELSET2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_BC_Type),
			1,
			MPI::INT,0);
  // Residual variable:
  MPI::COMM_WORLD.Bcast(&(IP.i_Residual_Variable),
			1,
			MPI::DOUBLE,0);
  // Interface input parameters:
  IP.Interface_IP.Broadcast_Input_Parameters();
  // Initial distance function parameters:
  MPI::COMM_WORLD.Bcast(IP.Initial_Distance_Type,
			INPUT_PARAMETER_LENGTH_LEVELSET2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Initial_Distance_Type),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Perturb_Distance_Function),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Extension_Distance),
			1,
			MPI::DOUBLE,0);
  // Eikonal equation parameters.
  MPI::COMM_WORLD.Bcast(&(IP.i_Redistance_Criteria),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Redistance_Frequency),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Redistance_Tolerance),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Initial_Redistance_Iterations),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Redistance_Iterations),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Eikonal_CFL_Number),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(IP.Eikonal_Scheme,
			INPUT_PARAMETER_LENGTH_LEVELSET2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Eikonal_Scheme),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(IP.Eikonal_Sign_Function,
			INPUT_PARAMETER_LENGTH_LEVELSET2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Eikonal_Sign_Function),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Eikonal_Threshold),
			1,
			MPI::DOUBLE,0);
  // Scalar (front speed) extension equation parameters.
  MPI::COMM_WORLD.Bcast(&(IP.Scalar_Extension_CFL_Number),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Scalar_Extension_Iterations),
			1,
			MPI::INT,0);
  // Pass the curvature motion parameters.
  MPI::COMM_WORLD.Bcast(&(IP.Curvature_Speed),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(IP.Curvature_Scheme,
			INPUT_PARAMETER_LENGTH_LEVELSET2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_Curvature_Scheme),
			1,
			MPI::INT,0);
  // Pass the bulk flowfield type.
  MPI::COMM_WORLD.Bcast(IP.BulkFlowField_Type,
			INPUT_PARAMETER_LENGTH_LEVELSET2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(&(IP.i_BulkFlowField_Type),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.V.x),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.V.y),
			1,
			MPI::DOUBLE,0);
  // Pass the grid variables.
  MPI::COMM_WORLD.Bcast(IP.Grid_Type,
			INPUT_PARAMETER_LENGTH_LEVELSET2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.NACA_Aerofoil_Type,
			INPUT_PARAMETER_LENGTH_LEVELSET2D,
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
			MPI::INT, 0);
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
  MPI::COMM_WORLD.Bcast(&(IP.Plate_Length),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Pipe_Length),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Pipe_Radius),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Blunt_Body_Radius),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Blunt_Body_Mach_Number),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Grain_Length),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Grain_Radius),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Grain_To_Throat_Length),
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
  MPI::COMM_WORLD.Bcast(&(IP.Cylinder_Radius),
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
  MPI::COMM_WORLD.Bcast(&(IP.Ellipse_Length_X_Axis),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Ellipse_Length_Y_Axis),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Chord_Length),
			1,
			MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Orifice_Radius),
			1,
			MPI::DOUBLE,0);
  // AMR parameters:
  MPI::COMM_WORLD.Bcast(&(IP.AMR), 
                        1, 
                        MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(IP.AMR_Frequency),
                        1,
                        MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Initial_Mesh_Refinements), 
                        1, 
                        MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Uniform_Mesh_Refinements),
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
                        MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Refinement_Criteria),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Refinement_Criteria_Curvature),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Refinement_Criteria_Zero_Level_Set),
			1,
			MPI::INT,0);
  // Smooth quad block flag:
  MPI::COMM_WORLD.Bcast(&(IP.i_Smooth_Quad_Block),
			1,
			MPI::INT,0);
  // File names.
  MPI::COMM_WORLD.Bcast(IP.Output_File_Name,
			INPUT_PARAMETER_LENGTH_LEVELSET2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.Grid_File_Name,
			INPUT_PARAMETER_LENGTH_LEVELSET2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.Grid_Definition_File_Name,
			INPUT_PARAMETER_LENGTH_LEVELSET2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.Restart_File_Name,
			INPUT_PARAMETER_LENGTH_LEVELSET2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.Gnuplot_File_Name,
			INPUT_PARAMETER_LENGTH_LEVELSET2D,
			MPI::CHAR,0);
  MPI::COMM_WORLD.Bcast(IP.Output_Format_Type,
			INPUT_PARAMETER_LENGTH_LEVELSET2D,
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
void Broadcast_Input_Parameters(LevelSet2D_Input_Parameters &IP,
                                MPI::Intracomm &Communicator,
                                const int Source_CPU) {

  int Source_Rank = 0;

  // CFFC path:
  Communicator.Bcast(IP.CFFC_Path, 
 		     INPUT_PARAMETER_LENGTH_LEVELSET2D, 
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.Input_File_Name,
		     INPUT_PARAMETER_LENGTH_LEVELSET2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.Line_Number),
		     1,
		     MPI::INT,Source_Rank);
  // Time integration parameters:
  Communicator.Bcast(IP.Time_Integration_Type,
		     INPUT_PARAMETER_LENGTH_LEVELSET2D,
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
  Communicator.Bcast(&(IP.Hamilton_Jacobi_CFL_Number),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Time_Max),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(IP.Reconstruction_Type,
		     INPUT_PARAMETER_LENGTH_LEVELSET2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Reconstruction),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(IP.Limiter_Type,
		     INPUT_PARAMETER_LENGTH_LEVELSET2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Limiter),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(IP.BC_Type,
		     INPUT_PARAMETER_LENGTH_LEVELSET2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_BC_Type),
		     1,
		     MPI::INT,Source_Rank);
  // Residual variable:
  Communicator.Bcast(&(IP.i_Residual_Variable),
		     1,
		     MPI::DOUBLE,Source_Rank);
  // Interface input parameters:
  IP.Interface_IP.Broadcast_Input_Parameters(Communicator,
					     Source_CPU);
  // Initial distance function parameters:
  Communicator.Bcast(IP.Initial_Distance_Type,
		     INPUT_PARAMETER_LENGTH_LEVELSET2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Initial_Distance_Type),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Perturb_Distance_Function),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Extension_Distance),
		     1,
		     MPI::DOUBLE,Source_Rank);
  // Eikonal equation parameters:
  Communicator.Bcast(&(IP.i_Redistance_Criteria),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Redistance_Frequency),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Redistance_Tolerance),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Initial_Redistance_Iterations),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Redistance_Iterations),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Eikonal_CFL_Number),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(IP.Eikonal_Scheme,
		     INPUT_PARAMETER_LENGTH_LEVELSET2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Eikonal_Scheme),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(IP.Eikonal_Sign_Function,
		     INPUT_PARAMETER_LENGTH_LEVELSET2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Eikonal_Sign_Function),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Eikonal_Threshold),
		     1,
		     MPI::DOUBLE,Source_Rank);
  // Scalar (front speed) extension equation parameters:
  Communicator.Bcast(&(IP.Scalar_Extension_CFL_Number),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Scalar_Extension_Iterations),
		     1,
		     MPI::INT,Source_Rank);
  // Pass the curvature motion parameters:
  Communicator.Bcast(&(IP.Curvature_Speed),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(IP.Curvature_Scheme,
		     INPUT_PARAMETER_LENGTH_LEVELSET2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Curvature_Scheme),
		     1,
		     MPI::INT,Source_Rank);
  // Pass the bulk flowfield type:
  Communicator.Bcast(IP.BulkFlowField_Type,
		     INPUT_PARAMETER_LENGTH_LEVELSET2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_BulkFlowField_Type),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.V.x),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.V.y),
		     1,
		     MPI::DOUBLE,Source_Rank);
  // Pass grid variables.
  Communicator.Bcast(IP.Grid_Type,
		     INPUT_PARAMETER_LENGTH_LEVELSET2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.NACA_Aerofoil_Type,
		     INPUT_PARAMETER_LENGTH_LEVELSET2D,
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
		     MPI::INT, Source_Rank);
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
  Communicator.Bcast(&(IP.Plate_Length),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Pipe_Length),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Pipe_Radius),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Blunt_Body_Radius),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Blunt_Body_Mach_Number),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Grain_Length),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Grain_Radius),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Grain_To_Throat_Length),
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
  Communicator.Bcast(&(IP.Cylinder_Radius),
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
  Communicator.Bcast(&(IP.Ellipse_Length_X_Axis),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Ellipse_Length_Y_Axis),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Chord_Length),
		     1,
		     MPI::DOUBLE,Source_Rank);
  Communicator.Bcast(&(IP.Orifice_Radius),
		     1,
		     MPI::DOUBLE,Source_Rank);
  // AMR parameters:
  Communicator.Bcast(&(IP.AMR), 
                     1, 
                     MPI::INT, Source_Rank);
  Communicator.Bcast(&(IP.AMR_Frequency),
                     1,
                     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Initial_Mesh_Refinements), 
                     1, 
                     MPI::INT, Source_Rank);
  Communicator.Bcast(&(IP.Number_of_Uniform_Mesh_Refinements),
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
  Communicator.Bcast(&(IP.Number_of_Refinement_Criteria),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Refinement_Criteria_Curvature),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Refinement_Criteria_Zero_Level_Set),
		     1,
		     MPI::INT,Source_Rank);
  // Smooth quad block flag:
  Communicator.Bcast(&(IP.i_Smooth_Quad_Block),
		     1,
		     MPI::INT,Source_Rank);
  // File Names:
  Communicator.Bcast(IP.Output_File_Name,
		     INPUT_PARAMETER_LENGTH_LEVELSET2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.Grid_File_Name,
		     INPUT_PARAMETER_LENGTH_LEVELSET2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.Grid_Definition_File_Name,
		     INPUT_PARAMETER_LENGTH_LEVELSET2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.Restart_File_Name,
		     INPUT_PARAMETER_LENGTH_LEVELSET2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.Gnuplot_File_Name,
		     INPUT_PARAMETER_LENGTH_LEVELSET2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(IP.Output_Format_Type,
		     INPUT_PARAMETER_LENGTH_LEVELSET2D,
		     MPI::CHAR,Source_Rank);
  Communicator.Bcast(&(IP.i_Output_Format),
		     1,
		     MPI::INT,Source_Rank);
  Communicator.Bcast(&(IP.Restart_Solution_Save_Frequency),
		     1,
		     MPI::INT,Source_Rank);
  // Output progress frequency:
  Communicator.Bcast(&(IP.Output_Progress_Frequency),
		     1,
		     MPI::INT,Source_Rank);
  // Number of blocks per processor:
  if (!CFFC_Primary_MPI_Processor()) {
    IP.Number_of_Processors = CFFC_MPI::Number_of_Processors;
  }
  Communicator.Bcast(&(IP.Number_of_Blocks_Per_Processor),
		     1,
		     MPI::INT,Source_Rank);

}
#endif

/**********************************************************************
 * Routine: Get_Next_Input_Control_Parameter                          *
 *                                                                    *
 * Get the next input control parameter from the input file.          *
 *                                                                    *
 **********************************************************************/
void Get_Next_Input_Control_Parameter(LevelSet2D_Input_Parameters &IP) {

  int i;
  char buffer[256];

  IP.Line_Number = IP.Line_Number + 1;
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
int Parse_Next_Input_Control_Parameter(LevelSet2D_Input_Parameters &IP) {

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
      IP.N_Stage = 4;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Reconstruction_Type") == 0) {
    i_command = 2;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Reconstruction_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Reconstruction_Type,"ENO") == 0 ||
	strcmp(IP.Reconstruction_Type,"Essentially_Non_Oscillatory") == 0 ||
	strcmp(IP.Reconstruction_Type,"Linear_ENO") == 0 ||
	strcmp(IP.Reconstruction_Type,"Linear_Essentially_Non_Oscillatory") == 0) {
      IP.i_Reconstruction = RECONSTRUCTION_LINEAR_ESSENTIALLY_NON_OSCILLATORY;
    } else if (strcmp(IP.Reconstruction_Type,"Quadratic_ENO") == 0 ||
	       strcmp(IP.Reconstruction_Type,"Quadratic_Essentially_Non_Oscillatory") == 0) {
      IP.i_Reconstruction = RECONSTRUCTION_QUADRATIC_ESSENTIALLY_NON_OSCILLATORY;
    } else if (strcmp(IP.Reconstruction_Type,"Cubic_ENO") == 0 ||
	       strcmp(IP.Reconstruction_Type,"Cubic_Essentially_Non_Oscillatory") == 0) {
      IP.i_Reconstruction = RECONSTRUCTION_CUBIC_ESSENTIALLY_NON_OSCILLATORY;
    } else if (strcmp(IP.Reconstruction_Type,"WENO") == 0 ||
	       strcmp(IP.Reconstruction_Type,"Weighted_Essentially_Non_Oscillatory") == 0) {
      IP.i_Reconstruction = RECONSTRUCTION_WEIGHTED_ESSENTIALLY_NON_OSCILLATORY;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Limiter_Type") == 0) {
    i_command = 3;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Limiter_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Limiter_Type,"One") == 0) {
      IP.i_Limiter = LIMITER_ONE;
    } else if (strcmp(IP.Limiter_Type,"Zero") == 0) {
      IP.i_Limiter = LIMITER_ZERO;
    } else if (strcmp(IP.Limiter_Type,"VanLeer") == 0) {
      IP.i_Limiter = LIMITER_VANLEER;
    } else if (strcmp(IP.Limiter_Type,"VanAlbada") == 0) {
      IP.i_Limiter = LIMITER_VANALBADA;
    } else if (strcmp(IP.Limiter_Type,"Barth_Jespersen") == 0) {
      IP.i_Limiter = LIMITER_BARTH_JESPERSEN;
    } else if (strcmp(IP.Limiter_Type,"Venkatakrishnan") == 0) {
      IP.i_Limiter = LIMITER_VENKATAKRISHNAN;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"BC_Type") == 0) {
    i_command = 4;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.BC_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.BC_Type,"Constant_Extrapolation") == 0) {
      IP.i_BC_Type = BC_CONSTANT_EXTRAPOLATION;
    } else if (strcmp(IP.BC_Type,"Linear_Extrapolation") == 0) {
      IP.i_BC_Type = BC_LINEAR_EXTRAPOLATION;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Initial_Distance_Type") == 0) {
    i_command = 5;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Initial_Distance_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.Initial_Distance_Type,"Geometric") == 0) {
      IP.i_Initial_Distance_Type = LEVELSET_INITIAL_EXTENSION_GEOMETRIC;
    } else if (strcmp(IP.Initial_Distance_Type,"Exact") == 0) {
      IP.i_Initial_Distance_Type = LEVELSET_INITIAL_EXTENSION_EXACT;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Perturb_Distance_Function") == 0) {
    i_command = 5;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
      IP.Perturb_Distance_Function = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
      IP.Perturb_Distance_Function = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Extension_Distance") == 0) {
    i_command = 5;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Extension_Distance;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Extension_Distance <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Redistance_Criteria") == 0) {
    i_command = 5;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Redistance_Criteria;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (strcmp(IP.Redistance_Criteria,"Threshold") == 0) {
      IP.i_Redistance_Criteria = EIKONAL_CRITERIA_THRESHOLD;
    } else if (strcmp(IP.Redistance_Criteria,"Frequency") == 0) {
      IP.i_Redistance_Criteria = EIKONAL_CRITERIA_FREQUENCY;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Redistance_Frequency") == 0) {
    i_command = 5;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Redistance_Frequency;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Redistance_Frequency < 0) i_command = INVALID_INPUT_VALUE;
    //if (IP.Redistance_Frequency > 1) IP.Redistance_Frequency = 1;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Initial_Redistance_Iterations") == 0) {
    i_command = 5;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Initial_Redistance_Iterations;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Initial_Redistance_Iterations < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Redistance_Iterations") == 0) {
    i_command = 5;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Redistance_Iterations;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Redistance_Iterations < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Redistance_Tolerance") == 0) {
    i_command = 5;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Redistance_Tolerance;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Redistance_Tolerance <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Eikonal_Solution_Type") == 0) {
    i_command = 5;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Eikonal_Scheme,IP.Next_Control_Parameter);
    if (strcmp(IP.Eikonal_Scheme,"Sussman") == 0) {
      IP.i_Eikonal_Scheme = EIKONAL_SCHEME_SUSSMAN;
    } else if (strcmp(IP.Eikonal_Scheme,"Russo_Smereka") == 0) {
      IP.i_Eikonal_Scheme = EIKONAL_SCHEME_RUSSO_SMEREKA;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Eikonal_Sign_Function") == 0) {
    i_command = 5;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Eikonal_Sign_Function,IP.Next_Control_Parameter);
    if (strcmp(IP.Eikonal_Sign_Function,"Discrete") == 0) {
      IP.i_Eikonal_Sign_Function = EIKONAL_SIGN_FUNCTION_DISCRETE;
    } else if (strcmp(IP.Eikonal_Sign_Function,"Smeared") == 0) {
      IP.i_Eikonal_Sign_Function = EIKONAL_SIGN_FUNCTION_SMEARED;
    } else if (strcmp(IP.Eikonal_Sign_Function,"Smeared_New") == 0) {
      IP.i_Eikonal_Sign_Function = EIKONAL_SIGN_FUNCTION_SMEARED_RUUTH;
    } else if (strcmp(IP.Eikonal_Sign_Function,"Derivative") == 0) {
      IP.i_Eikonal_Sign_Function = EIKONAL_SIGN_FUNCTION_DERIVATIVE;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Eikonal_Threshold") == 0) {
    i_command = 5;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Eikonal_Threshold;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Eikonal_Threshold < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Scalar_Extension_Iterations") == 0) {
    i_command = 5;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Scalar_Extension_Iterations;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Scalar_Extension_Iterations < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Curvature_Speed") == 0) {
    i_command = 5;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Curvature_Speed;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Curvature_Speed < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Curvature_Scheme") == 0) {
    i_command = 5;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Curvature_Scheme,IP.Next_Control_Parameter);
    if (strcmp(IP.Curvature_Scheme,"Laplacian") == 0) {
      IP.i_Curvature_Scheme = CURVATURE_SCHEME_LAPLACIAN;
    } else if (strcmp(IP.Curvature_Scheme,"Regular") == 0) {
      IP.i_Curvature_Scheme = CURVATURE_SCHEME_REGULAR;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Bulk_Flow_Field_Type") == 0) {
    i_command = 6;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.BulkFlowField_Type,IP.Next_Control_Parameter);
    if (strcmp(IP.BulkFlowField_Type,"None") == 0) {
      IP.i_BulkFlowField_Type = INTERFACE_BULKFLOWFIELD_NONE;
      IP.V = Vector2D_ZERO;
    } else if (strcmp(IP.BulkFlowField_Type,"Uniform") == 0) {
      IP.i_BulkFlowField_Type = INTERFACE_BULKFLOWFIELD_UNIFORM;
    } else if (strcmp(IP.BulkFlowField_Type,"Swirl") == 0) {
      IP.i_BulkFlowField_Type = INTERFACE_BULKFLOWFIELD_SWIRL;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Bulk_Flow_Characteristic_Speed") == 0) {
    i_command = 6;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.V.x >> IP.V.y;
    IP.Input_File.getline(buffer,sizeof(buffer));

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
    } else if (strcmp(IP.Grid_Type,"Rocket_Motor") == 0) {
      IP.i_Grid = GRID_ROCKET_MOTOR;
      IP.Grain_Length = 0.835;
      IP.Grain_Radius = 0.020;
      IP.Grain_To_Throat_Length = 0.05;
      IP.Nozzle_Length = 0.150;
      IP.Nozzle_Radius_Exit = 0.030;
      IP.Nozzle_Radius_Throat = 0.010;
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
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Cells_Idir;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Cells_Idir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Cells_Jdir") == 0) {
    i_command = 11;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Cells_Jdir;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Cells_Jdir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Ghost_Cells") == 0) {
    i_command = 12;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Ghost_Cells;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Ghost_Cells < 2) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Blocks_Idir") == 0) {
    i_command = 12;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Blocks_Idir;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Blocks_Idir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Blocks_Jdir") == 0) {
    i_command = 13;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Blocks_Jdir;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Blocks_Jdir < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Time_Accurate") == 0) {
    i_command = 14;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Time_Accurate;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Time_Accurate != 0 && IP.Time_Accurate != 1) IP.Time_Accurate = 0;
    if (IP.Time_Accurate) {
      IP.Local_Time_Stepping = 0;
    } else {
      IP.Local_Time_Stepping = 0;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Local_Time_Stepping") == 0) {
    i_command = 15;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Local_Time_Stepping;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Local_Time_Stepping != 0 && IP.Local_Time_Stepping != 1) IP.Local_Time_Stepping = 1;
    cout << endl << " Local_Time_Stepping is not allowed." << endl;
    i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Maximum_Number_of_Time_Steps") == 0) {
    i_command = 16;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Maximum_Number_of_Time_Steps;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Maximum_Number_of_Time_Steps < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"N_Stage") == 0) {
    i_command = 17;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.N_Stage;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.N_Stage < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Hamilton_Jacobi_CFL_Number") == 0) {
    i_command = 18;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Hamilton_Jacobi_CFL_Number;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Hamilton_Jacobi_CFL_Number <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Eikonal_CFL_Number") == 0) {
    i_command = 18;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Eikonal_CFL_Number;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Eikonal_CFL_Number <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Scalar_Extension_CFL_Number") == 0) {
    i_command = 18;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Scalar_Extension_CFL_Number;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Scalar_Extension_CFL_Number <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Box_Width") == 0) {
    i_command = 19;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Box_Width;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Box_Width <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Box_Height") == 0) {
    i_command = 20;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Box_Height;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Box_Height <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Plate_Length") == 0) {
    i_command = 21;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Plate_Length;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Plate_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Pipe_Length") == 0) {
    i_command = 22;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Pipe_Length;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Pipe_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Pipe_Radius") == 0) {
    i_command = 23;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Pipe_Radius;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Pipe_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Blunt_Body_Radius") == 0) {
    i_command = 24;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Blunt_Body_Radius;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Blunt_Body_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Blunt_Body_Mach_Number") == 0) {
    i_command = 25;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Blunt_Body_Mach_Number;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Blunt_Body_Mach_Number <= ONE) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Cylinder_Radius") == 0) {
    i_command = 26;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Cylinder_Radius;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Cylinder_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Ellipse_Length_X_Axis") == 0) {
    i_command = 27;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Ellipse_Length_X_Axis;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Ellipse_Length_X_Axis <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Ellipse_Length_Y_Axis") == 0) {
    i_command = 28;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Ellipse_Length_Y_Axis;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Ellipse_Length_Y_Axis <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Chord_Length") == 0) {
    i_command = 29;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Chord_Length;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Chord_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"NACA_Aerofoil_Type") == 0) {
    i_command = 30;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.NACA_Aerofoil_Type,IP.Next_Control_Parameter);
    if (strlen(IP.NACA_Aerofoil_Type) != 4 &&
	strlen(IP.NACA_Aerofoil_Type) != 5) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Orifice_Radius") == 0) {
    i_command = 31;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Orifice_Radius;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Orifice_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Grain_Length") == 0) {
    i_command = 32;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Grain_Length;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Grain_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Grain_Radius") == 0) {
    i_command = 33;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Grain_Radius;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Grain_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Grain_To_Throat_Length") == 0) {
    i_command = 34;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Grain_To_Throat_Length;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Grain_To_Throat_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Nozzle_Length") == 0) {
    i_command = 35;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Nozzle_Length;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Nozzle_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Nozzle_Radius_Exit") == 0) {
    i_command = 36;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Nozzle_Radius_Exit;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Nozzle_Radius_Exit <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Nozzle_Radius_Throat") == 0) {
    i_command = 37;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Nozzle_Radius_Throat;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Nozzle_Radius_Throat <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Time_Max") == 0) {
    i_command = 45;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Time_Max;
    IP.Input_File.getline(buffer,sizeof(buffer));
    IP.Time_Max = IP.Time_Max/THOUSAND;
    if (IP.Time_Max < ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Blocks_Per_Processor") == 0) {
    i_command = 46;
    IP.Line_Number = IP.Line_Number + 1;
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

  } else if (strcmp(IP.Next_Control_Parameter,"Restart_Solution_Save_Frequency") == 0) {
    i_command = 49;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Restart_Solution_Save_Frequency;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Restart_Solution_Save_Frequency < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Output_Progress_Frequency") == 0) {
    i_command = 50;
    IP.Line_Number++;
    IP.Input_File >> IP.Output_Progress_Frequency;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Output_Progress_Frequency < 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"X_Shift") == 0) {
    i_command = 56;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.X_Shift;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"X_Scale") == 0) {
    i_command = 57;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.X_Scale;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"X_Rotate") == 0) {
    i_command = 58;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.X_Rotate;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "AMR") == 0) {
    i_command = 59;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
      IP.AMR = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
      IP.AMR = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter, "AMR_Frequency") == 0) {
    i_command = 60;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.AMR_Frequency;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.AMR_Frequency < 0) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Initial_Mesh_Refinements") == 0) {
    i_command = 61;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Initial_Mesh_Refinements;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Initial_Mesh_Refinements < 0) IP.Number_of_Initial_Mesh_Refinements = 0;
    if (IP.Number_of_Initial_Mesh_Refinements > IP.Maximum_Refinement_Level)
      IP.Maximum_Refinement_Level = IP.Number_of_Initial_Mesh_Refinements;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Uniform_Mesh_Refinements") == 0) {
    i_command = 62;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Uniform_Mesh_Refinements;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Uniform_Mesh_Refinements < 0) IP.Number_of_Uniform_Mesh_Refinements = 0;
    if (IP.Number_of_Uniform_Mesh_Refinements > IP.Maximum_Refinement_Level)
      IP.Maximum_Refinement_Level = IP.Number_of_Uniform_Mesh_Refinements;

  } else if (strcmp(IP.Next_Control_Parameter,"Maximum_Refinement_Level") == 0) {
    i_command = 65;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Maximum_Refinement_Level;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Maximum_Refinement_Level < 1) IP.Maximum_Refinement_Level = 1;

  } else if (strcmp(IP.Next_Control_Parameter,"Minimum_Refinement_Level") == 0) {
    i_command = 66;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Minimum_Refinement_Level;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Minimum_Refinement_Level < 1) IP.Minimum_Refinement_Level = 1;

  } else if (strcmp(IP.Next_Control_Parameter, "Threshold_for_Refinement") == 0) {
    i_command = 67;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Threshold_for_Refinement;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Threshold_for_Refinement <= ZERO ||
	IP.Threshold_for_Refinement > ONE) IP.Threshold_for_Refinement = 0.50;

  } else if (strcmp(IP.Next_Control_Parameter, "Threshold_for_Coarsening") == 0) {
    i_command = 68;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Threshold_for_Coarsening;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (IP.Threshold_for_Coarsening < ZERO ||
	IP.Threshold_for_Coarsening >= ONE) IP.Threshold_for_Coarsening = 0.10;

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Refinement_Criteria") == 0) {
    i_command = 69;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Number_of_Refinement_Criteria;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Number_of_Refinement_Criteria < 1 || IP.Number_of_Refinement_Criteria > 6) {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Curvature") == 0) {
    i_command = 70;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
      IP.Refinement_Criteria_Curvature = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
      IP.Refinement_Criteria_Curvature = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Zero_Level_Set") == 0) {
    i_command = 72;
    Get_Next_Input_Control_Parameter(IP);
    if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
      IP.Refinement_Criteria_Zero_Level_Set = ON;
    } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
      IP.Refinement_Criteria_Zero_Level_Set = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Smooth_Quad_Block") == 0) {
    i_command = 73;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.i_Smooth_Quad_Block;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.i_Smooth_Quad_Block < 0 || IP.i_Smooth_Quad_Block > 1) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Residual_Variable") == 0) {
    i_command = 81;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.i_Residual_Variable;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.i_Residual_Variable != 1) i_command = INVALID_INPUT_VALUE;

    ////////////////////////////////////////////////////////////////////
    // INTERFACE PARAMETERS                                           //
    ////////////////////////////////////////////////////////////////////

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Interface_Components") == 0) {
    i_command = 600;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Interface_IP.Number_of_Components;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Interface_IP.Number_of_Components < 0) {
      i_command = INVALID_INPUT_VALUE;
    } else if (IP.Interface_IP.Number_of_Components > 0) {
      IP.Interface_IP.Component_List.allocate(IP.Interface_IP.Number_of_Components);
      IP.Interface_IP.ns = new int[IP.Interface_IP.Number_of_Components];
    }

  } else if (strcmp(IP.Next_Control_Parameter,"Interface_Spline_Type") == 0) {
    i_command = 601;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Interface_IP.Type,IP.Next_Control_Parameter);
    IP.Interface_IP.ci++; // increment the interface counter.
    if (strcmp(IP.Interface_IP.Type,"Circle") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_CIRCLE;
    } else if (strcmp(IP.Interface_IP.Type,"Ellipse") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_ELLIPSE;
    } else if (strcmp(IP.Interface_IP.Type,"Square") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_SQUARE;
    } else if (strcmp(IP.Interface_IP.Type,"Rectangle") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_RECTANGLE;
    } else if (strcmp(IP.Interface_IP.Type,"Rocket_Propellant_Grain") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_ROCKET_PROPELLANT_GRAIN;
    } else if (strcmp(IP.Interface_IP.Type,"NACA0012_Aerofoil") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_NACA0012_AEROFOIL;
    } else if (strcmp(IP.Interface_IP.Type,"NACA0015_Aerofoil") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_NACA0015_AEROFOIL;
    } else if (strcmp(IP.Interface_IP.Type,"Zalesaks_Disk") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_ZALESAK;
    } else if (strcmp(IP.Interface_IP.Type,"User_Specified") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_USER_SPECIFIED;
    } else if (strcmp(IP.Interface_IP.Type,"Restart") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_RESTART;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }
    IP.Interface_IP.Component_List[IP.Interface_IP.ci].allocate(2);

  } else if (strcmp(IP.Next_Control_Parameter,"Interface_Characteristic_Length") == 0 ||
	     strcmp(IP.Next_Control_Parameter,"Interface_Characteristic_Length_1") == 0) {
    i_command = 602;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Interface_IP.Component_List[IP.Interface_IP.ci].Length1;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Interface_IP.Component_List[IP.Interface_IP.ci].Length1 <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Interface_Characteristic_Length_2") == 0) {
    i_command = 603;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Interface_IP.Component_List[IP.Interface_IP.ci].Length2;
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Interface_IP.Component_List[IP.Interface_IP.ci].Length2 <= ZERO) i_command = INVALID_INPUT_VALUE;

  } else if (strcmp(IP.Next_Control_Parameter,"Interface_Reference_Point") == 0) {
    i_command = 605;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Interface_IP.Component_List[IP.Interface_IP.ci].Xref.x
		  >> IP.Interface_IP.Component_List[IP.Interface_IP.ci].Xref.y;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"Interface_Motion") == 0) {
    i_command = 606;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Interface_IP.Motion,IP.Next_Control_Parameter);
    if (strcmp(IP.Interface_IP.Motion,"Stationary") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_LEVELSET_STATIONARY;
    } else if (strcmp(IP.Interface_IP.Motion,"Constant") == 0 ||
	       strcmp(IP.Interface_IP.Motion,"Expand") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_LEVELSET_EXPAND;
    } else if (strcmp(IP.Interface_IP.Motion,"Stretch") == 0) {
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_LEVELSET_STRETCH;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }
    IP.Interface_IP.Component_List[IP.Interface_IP.ci].Speed = Vector2D_ZERO;

  } else if (strcmp(IP.Next_Control_Parameter,"Interface_Characteristic_Velocity") == 0) {
    i_command = 607;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Interface_IP.Component_List[IP.Interface_IP.ci].Speed.x
		  >> IP.Interface_IP.Component_List[IP.Interface_IP.ci].Speed.y;
    IP.Input_File.getline(buffer,sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Specified_Interface_Points") == 0) {
    i_command = 608;
    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File >> IP.Interface_IP.ns[IP.Interface_IP.ci-1];
    IP.Input_File.getline(buffer,sizeof(buffer));
    if (IP.Interface_IP.ns[IP.Interface_IP.ci-1] <  0) i_command = INVALID_INPUT_VALUE;
    else IP.Interface_IP.Component_List[IP.Interface_IP.ci].allocate(max(IP.Interface_IP.ns[IP.Interface_IP.ci-1],2));

  } else if (strcmp(IP.Next_Control_Parameter,"Interface_Points") == 0) {
    i_command = 609;
    for (int np = 0; np < IP.Interface_IP.ns[IP.Interface_IP.ci-1]; np++) {
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> Xt.x >> Xt.y >> tpt >> bct;
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Spline.Xp[np] = Xt;
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Spline.tp[np] = tpt;
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Spline.bc[np] = bct;
      IP.Input_File.getline(buffer,sizeof(buffer));
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

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Nodes") == 0) {
    i_command = WRITE_OUTPUT_NODES_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Interface_List") == 0) {
    i_command = WRITE_OUTPUT_INTERFACE_COMPONENT_LIST_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Circle") == 0) {
    i_command = WRITE_OUTPUT_LEVEL_SET_CIRCLE_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Ellipse") == 0) {
    i_command = WRITE_OUTPUT_LEVEL_SET_ELLIPSE_CODE;
    i_command = INVALID_INPUT_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Zalesak") == 0) {
    i_command = WRITE_OUTPUT_LEVEL_SET_ZALESAK_CODE;
    i_command = INVALID_INPUT_CODE;

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

  } else if (IP.Next_Control_Parameter[0] == '#') {
    i_command = COMMENT_CODE;

  } else if (strcmp(IP.Next_Control_Parameter,"Refine_Grid") == 0) {
    i_command = REFINE_GRID_CODE;

  } else {
    i_command = INVALID_INPUT_CODE;

  }

  // Return the parser command type indicator.
  return i_command;

}

/**********************************************************************
 * Routine: Process_Input_Control_Parameter_File                      *
 *                                                                    *
 * Reads,parses,and executes the list of input control parameters     *
 * from the standard input file.                                      *
 *                                                                    *
 **********************************************************************/
int Process_Input_Control_Parameter_File(LevelSet2D_Input_Parameters &Input_Parameters,
                                         char *Input_File_Name_ptr,
                                         int  &Command_Flag) {

  int error_flag,line_number;

  // Assign initial value for error indicator flag.
  error_flag = 0;

  // Assign default values to the input parameters.
  Set_Default_Input_Parameters(Input_Parameters);

  // Copy input file name (a string) to appropriate input parameter variable.
  if (Input_File_Name_ptr != NULL) 
    strcpy(Input_Parameters.Input_File_Name,Input_File_Name_ptr);

  // Open the input file containing the input parameters.
  Open_Input_File(Input_Parameters);
  error_flag = Input_Parameters.Input_File.fail();

  if (error_flag) {
    cout << "\n LevelSet2D ERROR: Unable to open LevelSet2D input data file.\n";
    return error_flag;
  }

  // Read and parse control parameters contained in the input file.
  while (1) {
    Get_Next_Input_Control_Parameter(Input_Parameters);
    Command_Flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
    line_number = Input_Parameters.Line_Number;
    if (Command_Flag == EXECUTE_CODE) {
      break;
    } else if (Command_Flag == TERMINATE_CODE) {
      break;
    } else if (Command_Flag == INVALID_INPUT_CODE ||
	       Command_Flag == INVALID_INPUT_VALUE) {
      line_number = -line_number;
      cout << "\n LevelSet2D ERROR: Error reading LevelSet2D data at line #"
	   << -line_number  << " of input data file.\n";
      error_flag = line_number;
      break;
    }
  }

  // Perform consistency checks on the refinement criteria.
  Input_Parameters.Number_of_Refinement_Criteria = 0;
  if (Input_Parameters.Refinement_Criteria_Curvature) Input_Parameters.Number_of_Refinement_Criteria++;
  if (Input_Parameters.Refinement_Criteria_Zero_Level_Set) {
    Input_Parameters.Number_of_Refinement_Criteria = 1;
    Input_Parameters.Refinement_Criteria_Curvature = OFF;
  }

  // Initial processing of input control parameters complete.
  // Return the error indicator flag.
  return error_flag;

}
