/* Rte2DInput.cc:  Subroutines for 2D Rte Input Classes. */

/* Include 2D Rte input parameter header file. */

#ifndef _RTE2D_INPUT_INCLUDED
#include "Rte2DInput.h"
#endif // _RTE2D_INPUT_INCLUDED

/*************************************************************
 * Rte2D_Input_Parameters -- External subroutines.         *
 *************************************************************/

/********************************************************
 * Routine: Open_Input_File                             *
 *                                                      *
 * Opens the appropriate input data file.               *
 *                                                      *
 ********************************************************/
void Open_Input_File(Rte2D_Input_Parameters &IP) {

    IP.Input_File.open(IP.Input_File_Name, ios::in);
    if (!IP.Input_File.bad()) {
       IP.Line_Number = 0;
       IP.Input_File.setf(ios::skipws);
    } /* endif */

}

/********************************************************
 * Routine: Close_Input_File                            *
 *                                                      *
 * Closes the appropriate input data file.              *
 *                                                      *
 ********************************************************/
void Close_Input_File(Rte2D_Input_Parameters &IP) {

    IP.Input_File.unsetf(ios::skipws);
    IP.Input_File.close();

}

/********************************************************
 * Routine: Set_Default_Input_Parameters                *
 *                                                      *
 * Assigns default values to the input parameters.      *
 *                                                      *
 ********************************************************/
void Set_Default_Input_Parameters(Rte2D_Input_Parameters &IP) {

    int i;
    char *string_ptr;

    string_ptr = "Rte2D.in";
    strcpy(IP.Input_File_Name, string_ptr);

    string_ptr = "Explicit_Euler";
    strcpy(IP.Time_Integration_Type, string_ptr);
    IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
    IP.Time_Accurate = 0;
    IP.Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
    IP.Maximum_Number_of_Time_Steps = 100;
    IP.N_Stage = 1;
    IP.Time_Max = ZERO;

    IP.Residual_Smoothing = 0;
    IP.Residual_Smoothing_Epsilon = ZERO;
    IP.Residual_Smoothing_Gauss_Seidel_Iterations = 2;

    string_ptr = "Least_Squares";
    strcpy(IP.Reconstruction_Type, string_ptr);
    IP.i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;

    string_ptr = "Barth_Jespersen";
    strcpy(IP.Limiter_Type, string_ptr);
    IP.i_Limiter = LIMITER_BARTH_JESPERSEN;

    string_ptr = "Upwind_Difference";
    strcpy(IP.SpaceMarch_Scheme, string_ptr);
    IP.i_Limiter = SPACE_MARCH_UPWIND;

    string_ptr = "Uniform";
    strcpy(IP.ICs_Type, string_ptr);
    IP.i_ICs = IC_UNIFORM;

    string_ptr = "Planar";
    strcpy(IP.Flow_Geometry_Type, string_ptr);
    IP.Axisymmetric = 0;

    //------------------------ Rte2D Specific -------------------------//
    string_ptr = "FVM";
    strcpy(IP.RTE_Solver, string_ptr);
    IP.i_RTE_Solver = RTE2D_SOLVER_FVM;
    IP.Uo.RTE_Type = IP.i_RTE_Solver;
    
    string_ptr = "S2";
    strcpy(IP.DOM_Quadrature, string_ptr);
    IP.i_DOM_Quadrature = DOM_S2;

    Rte2D_State::SetGas( );
    IP.Number_of_Angles_Mdir = 4;
    IP.Number_of_Angles_Ldir = 4;
    Rte2D_State::SetDirs( IP.Number_of_Angles_Mdir, 
			  IP.Number_of_Angles_Ldir, 
			  IP.i_DOM_Quadrature,
			  IP.Axisymmetric );
    IP.Uo.Allocate();
    IP.Uo.Zero();

    IP.i_ScatteringFunc = SCATTER_ISO;
    string_ptr = "Isotropic";
    strcpy(IP.ScatteringFunc, string_ptr);

    IP.AbsorbsionCoef = ONE;
    IP.ScatteringCoef = ZERO;
    IP.Temperature = THOUSAND;
    IP.Uo.SetAbsorbsion( IP.AbsorbsionCoef );
    IP.Uo.SetScattering( IP.ScatteringCoef );
    IP.Uo.SetBlackbody( Ib(IP.Temperature) );

    Rte2D_State :: SetupPhase( IP.i_ScatteringFunc );

    IP.NorthWallTemp = ZERO;
    IP.SouthWallTemp = ZERO;
    IP.EastWallTemp = ZERO;
    IP.WestWallTemp = ZERO;
    IP.NorthWallEmiss = ZERO;      
    IP.SouthWallEmiss = ZERO;      
    IP.EastWallEmiss = ZERO;      
    IP.WestWallEmiss = ZERO;     

    //---------------------- End Rte2D Specific -----------------------//


    string_ptr = "Square";
    strcpy(IP.Grid_Type, string_ptr);
    IP.i_Grid = GRID_SQUARE;
    IP.Box_Width = ONE;
    IP.Box_Height = ONE;
    IP.Number_of_Cells_Idir = 100;
    IP.Number_of_Cells_Jdir = 100;
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Ghost_Cells = 2;

    IP.Plate_Length = ONE;
    IP.Pipe_Length = ONE;
    IP.Pipe_Radius = HALF;
    IP.Blunt_Body_Radius = ONE;
    IP.Blunt_Body_Mach_Number = TWO;
    IP.Grain_Length = 0.835;
    IP.Grain_Radius = 0.020;
    IP.Grain_To_Throat_Length = 0.05;
    IP.Nozzle_Length = 0.150;
    IP.Nozzle_Radius_Exit = 0.030;
    IP.Nozzle_Radius_Throat = 0.010;
    IP.Cylinder_Radius = ONE;
    IP.Ellipse_Length_X_Axis = TWO;
    IP.Ellipse_Length_Y_Axis = HALF;
    IP.Chord_Length = ONE;
    IP.Orifice_Radius = ONE;
    IP.Wedge_Angle = 25.0;
    IP.Wedge_Length = HALF;

    IP.X_Shift = Vector2D_ZERO;
    IP.X_Scale = ONE;
    IP.X_Rotate = ZERO;

    IP.ICEMCFD_FileNames = ICEMCFD_get_filenames();

    IP.AMR = 0;
    IP.AMR_Frequency = 100;
    IP.Number_of_Initial_Mesh_Refinements = 0;
    IP.Number_of_Uniform_Mesh_Refinements = 0;
    IP.Number_of_Boundary_Mesh_Refinements = 0;
    IP.Maximum_Refinement_Level = 100;
    IP.Minimum_Refinement_Level = 1;

    IP.Threshold_for_Refinement = 0.50;
    IP.Threshold_for_Coarsening = 0.10;

    IP.Morton = 0;
    IP.Morton_Reordering_Frequency = 1000;

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

    string_ptr = " ";
    strcpy(IP.Next_Control_Parameter, string_ptr);

    IP.Line_Number = 0;

    IP.Number_of_Processors = CFDkit_MPI::Number_of_Processors;
    IP.Number_of_Blocks_Per_Processor = 10;
    
    // Freezing_Limiter
    IP.Freeze_Limiter = 0;

    // Limiter_switch
    IP.Freeze_Limiter_Residual_Level = 1e-4;

}

/********************************************************
 * Routine: Broadcast_Input_Parameters                  *
 *                                                      *
 * Broadcast the input parameters variables to all      *
 * processors involved in the calculation from the      *
 * primary processor using the MPI broadcast routine.   *
 *                                                      *
 ********************************************************/
void Broadcast_Input_Parameters(Rte2D_Input_Parameters &IP) {

#ifdef _MPI_VERSION
    int i;

    MPI::COMM_WORLD.Bcast(IP.Input_File_Name, 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Line_Number), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Time_Integration_Type, 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
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
                          INPUT_PARAMETER_LENGTH_RTE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Reconstruction), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Limiter_Type, 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Limiter), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.SpaceMarch_Scheme, 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_SpaceMarch_Scheme), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.ICs_Type, 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_ICs), 
                          1, 
                          MPI::INT, 0);
    //------------------------ Rte2D Specific -------------------------//
    MPI::COMM_WORLD.Bcast(IP.RTE_Solver, 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_RTE_Solver), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.DOM_Quadrature, 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_DOM_Quadrature), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Temperature), 
                          1, 
                          MPI::DOUBLE, 0);
     MPI::COMM_WORLD.Bcast(&(IP.AbsorbsionCoef), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.ScatteringCoef), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(IP.ScatteringFunc, 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_ScatteringFunc), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Angles_Mdir), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Angles_Ldir), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.NorthWallTemp), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.SouthWallTemp), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.EastWallTemp), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.WestWallTemp), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.NorthWallEmiss), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.SouthWallEmiss), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.EastWallEmiss), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.WestWallEmiss), 
                          1, 
                          MPI::DOUBLE, 0);
   if (!CFDkit_Primary_MPI_Processor()) {
     IP.Uo.RTE_Type = IP.i_RTE_Solver;
     Rte2D_State::SetGas( );
     Rte2D_State::SetDirs( IP.Number_of_Angles_Mdir, 
			   IP.Number_of_Angles_Ldir, 
			   IP.i_DOM_Quadrature,
			   IP.Axisymmetric);
    IP.Uo.Allocate();
    IP.Uo.Zero();
    IP.Uo.SetAbsorbsion( IP.AbsorbsionCoef );
    IP.Uo.SetScattering( IP.ScatteringCoef );
    IP.Uo.SetBlackbody( Ib(IP.Temperature) );
    Rte2D_State :: SetupPhase( IP.i_ScatteringFunc );
    } /* endif */
    //---------------------- End Rte2D Specific -----------------------//
    MPI::COMM_WORLD.Bcast(IP.Flow_Geometry_Type, 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Axisymmetric), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Grid_Type, 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.NACA_Aerofoil_Type, 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
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
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Idir), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Jdir), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Ghost_Cells), 
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
    MPI::COMM_WORLD.Bcast(&(IP.Grain_Length), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Grain_Radius), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Grain_To_Throat_Length), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Nozzle_Length), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Nozzle_Radius_Exit), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Nozzle_Radius_Throat), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Cylinder_Radius), 
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
    MPI::COMM_WORLD.Bcast(&(IP.Wedge_Angle), 
			  1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Wedge_Length), 
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
    if (!CFDkit_Primary_MPI_Processor()) {
       IP.ICEMCFD_FileNames = new char*[3];
       for (i = 0; i < 3; i++) {
          IP.ICEMCFD_FileNames[i] = new char[INPUT_PARAMETER_LENGTH_RTE2D];
       } /* endfor */
    } /* endif */
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[0], 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[1], 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[2], 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
                          MPI::CHAR, 0);

    // AMR & Refinement Parameters
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

    // Morton Ordering Parameters
    MPI::COMM_WORLD.Bcast(&(IP.Morton), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Morton_Reordering_Frequency),
                          1,
                          MPI::INT,0);

    // File Names
    MPI::COMM_WORLD.Bcast(IP.Output_File_Name, 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Grid_File_Name, 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Grid_Definition_File_Name, 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Restart_File_Name, 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Gnuplot_File_Name, 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Output_Format_Type, 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Output_Format), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Restart_Solution_Save_Frequency), 
                          1, 
                          MPI::INT, 0);

    // Multigrid Related Parameters
    IP.Multigrid_IP.Broadcast_Input_Parameters();

    // NKS Parameters
    IP.NKS_IP.Broadcast_Input_Parameters();

    if (!CFDkit_Primary_MPI_Processor()) {
       IP.Number_of_Processors = CFDkit_MPI::Number_of_Processors;
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

#endif

}

#ifdef _MPI_VERSION
/********************************************************
 * Routine: Broadcast_Input_Parameters                  *
 *                                                      *
 * Broadcast the input parameters variables to all      *
 * processors associated with the specified communicator*
 * from the specified processor using the MPI broadcast *
 *routine.                                              *
 *                                                      *
 ********************************************************/
void Broadcast_Input_Parameters(Rte2D_Input_Parameters &IP,
                                MPI::Intracomm &Communicator, 
                                const int Source_CPU) {

    int Source_Rank = 0;
    int i;

    Communicator.Bcast(IP.Input_File_Name, 
                       INPUT_PARAMETER_LENGTH_RTE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.Line_Number), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Time_Integration_Type, 
                       INPUT_PARAMETER_LENGTH_RTE2D, 
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
                       INPUT_PARAMETER_LENGTH_RTE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Reconstruction), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Limiter_Type, 
                       INPUT_PARAMETER_LENGTH_RTE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Limiter), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.SpaceMarch_Scheme, 
                       INPUT_PARAMETER_LENGTH_RTE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_SpaceMarch_Scheme), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.ICs_Type, 
                       INPUT_PARAMETER_LENGTH_RTE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_ICs), 
                       1, 
                       MPI::INT, Source_Rank);
    //------------------------ Rte2D Specific -------------------------//
    Communicator.Bcast(IP.RTE_Solver, 
                       INPUT_PARAMETER_LENGTH_RTE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_RTE_Solver), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.DOM_Quadrature, 
                       INPUT_PARAMETER_LENGTH_RTE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_DOM_Quadrature), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Temperature), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.AbsorbsionCoef), 
                          1, 
                          MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.ScatteringCoef), 
                          1, 
                          MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(IP.ScatteringFunc, 
                          INPUT_PARAMETER_LENGTH_RTE2D, 
                          MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_ScatteringFunc), 
                          1, 
                          MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Angles_Mdir), 
                          1, 
                          MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Angles_Ldir), 
                          1, 
                          MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.NorthWallTemp), 
                          1, 
                          MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.SouthWallTemp), 
                          1, 
                          MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.EastWallTemp), 
                          1, 
                          MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.WestWallTemp), 
                          1, 
                          MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.NorthWallEmiss), 
                          1, 
                          MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.SouthWallEmiss), 
                          1, 
                          MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.EastWallEmiss), 
                          1, 
                          MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.WestWallEmiss), 
                          1, 
                          MPI::DOUBLE, Source_Rank);
    if (!(CFDkit_MPI::This_Processor_Number == Source_CPU)) {
      IP.Uo.RTE_Type = IP.i_RTE_Solver;
      Rte2D_State::SetGas( );
      Rte2D_State::SetDirs( IP.Number_of_Angles_Mdir, 
			    IP.Number_of_Angles_Ldir, 
			    IP.i_DOM_Quadrature,
			    IP.Axisymmetric);
      IP.Uo.Allocate();
      IP.Uo.Zero();
      IP.Uo.SetAbsorbsion( IP.AbsorbsionCoef );
      IP.Uo.SetScattering( IP.ScatteringCoef );
      IP.Uo.SetBlackbody( Ib(IP.Temperature) );
      Rte2D_State :: SetupPhase( IP.i_ScatteringFunc );
    } /* endif */
    //---------------------- End Rte2D Specific -----------------------//
    Communicator.Bcast(IP.Flow_Geometry_Type, 
                       INPUT_PARAMETER_LENGTH_RTE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.Axisymmetric), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Grid_Type, 
                       INPUT_PARAMETER_LENGTH_RTE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.NACA_Aerofoil_Type, 
                       INPUT_PARAMETER_LENGTH_RTE2D, 
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
    Communicator.Bcast(&(IP.Number_of_Blocks_Idir), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Blocks_Jdir), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Ghost_Cells), 
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
    Communicator.Bcast(&(IP.Grain_Length), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Grain_Radius), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Grain_To_Throat_Length), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Nozzle_Length), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Nozzle_Radius_Exit), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Nozzle_Radius_Throat), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Cylinder_Radius), 
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
    Communicator.Bcast(&(IP.Wedge_Angle), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Wedge_Length), 
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
    if (!(CFDkit_MPI::This_Processor_Number == Source_CPU)) {
       IP.ICEMCFD_FileNames = new char*[3];
       for (i = 0; i < 3; i++) {
          IP.ICEMCFD_FileNames[i] = new char[INPUT_PARAMETER_LENGTH_RTE2D];
       } /* endfor */
    } /* endif */
    Communicator.Bcast(IP.ICEMCFD_FileNames[0], 
                       INPUT_PARAMETER_LENGTH_RTE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.ICEMCFD_FileNames[1], 
                       INPUT_PARAMETER_LENGTH_RTE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.ICEMCFD_FileNames[2], 
                       INPUT_PARAMETER_LENGTH_RTE2D, 
                       MPI::CHAR, Source_Rank);

    // AMR & Refinement Parameters
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

   // Morton Ordering Parameters
    Communicator.Bcast(&(IP.Morton), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Morton_Reordering_Frequency),
                       1,
                       MPI::INT,Source_Rank);

    // File Names
    Communicator.Bcast(IP.Output_File_Name, 
                       INPUT_PARAMETER_LENGTH_RTE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Grid_File_Name, 
                       INPUT_PARAMETER_LENGTH_RTE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Grid_Definition_File_Name, 
                       INPUT_PARAMETER_LENGTH_RTE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Restart_File_Name, 
                       INPUT_PARAMETER_LENGTH_RTE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Gnuplot_File_Name, 
                       INPUT_PARAMETER_LENGTH_RTE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Output_Format_Type, 
                       INPUT_PARAMETER_LENGTH_RTE2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Output_Format), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Restart_Solution_Save_Frequency), 
                       1, 
                       MPI::INT, Source_Rank);

    // Multigrid Related Parameters
    IP.Multigrid_IP.Broadcast_Input_Parameters(Communicator,
					       Source_CPU);   
    
    // Multigrid Related Parameters
    IP.NKS_IP.Broadcast_Input_Parameters(Communicator,
					 Source_CPU);

    if (!(CFDkit_MPI::This_Processor_Number == Source_CPU)) {
       IP.Number_of_Processors = CFDkit_MPI::Number_of_Processors;
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

}
#endif

/********************************************************
 * Routine: Get_Next_Input_Control_Parameter            *
 *                                                      *
 * Get the next input control parameter from the input  *
 * file.                                                *
 *                                                      *
 ********************************************************/
void Get_Next_Input_Control_Parameter(Rte2D_Input_Parameters &IP) {

    int i;
    char buffer[256];

    IP.Line_Number = IP.Line_Number + 1;
    IP.Input_File.getline(buffer, sizeof(buffer));
    i = 0;
    if (buffer[0] != '#') {
       while (1) {
          if (buffer[i] == ' ' || buffer[i] == '=' ) break;
          i = i + 1;
          if (i > strlen(buffer) ) break;
       } /* endwhile */
       buffer[i] = '\0';
    } /* endif */
    strcpy(IP.Next_Control_Parameter, buffer);

}

/********************************************************
 * Routine: Parse_Next_Input_Control_Parameter          *
 *                                                      *
 * Parses and executes the next input control parameter *
 * from the input file.                                 *
 *                                                      *
 ********************************************************/
int Parse_Next_Input_Control_Parameter(Rte2D_Input_Parameters &IP) {

    int i_command;
    char buffer[256];

    i_command = 0;

    if (strcmp(IP.Next_Control_Parameter, "Time_Integration_Type") == 0) {
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
       } else if (strcmp(IP.Time_Integration_Type, "Space_March") == 0) {
	 IP.i_Time_Integration = TIME_STEPPING_SPACE_MARCH;
       } else {
           IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
           IP.N_Stage = 1;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Reconstruction_Type") == 0) {
       i_command = 2;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Reconstruction_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Reconstruction_Type, "Green_Gauss") == 0) {
          IP.i_Reconstruction = RECONSTRUCTION_GREEN_GAUSS;
       } else if (strcmp(IP.Reconstruction_Type, "Least_Squares") == 0) {
          IP.i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;
       } else {
          IP.i_Reconstruction = RECONSTRUCTION_GREEN_GAUSS;
       } /* endif */

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
          IP.i_Limiter = LIMITER_VANLEER ;
       } /* endif */


    } else if (strcmp(IP.Next_Control_Parameter, "ICs_Type") == 0) {
       i_command = 5;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.ICs_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.ICs_Type, "Constant") == 0) {
          IP.i_ICs = IC_CONSTANT;
       } else if (strcmp(IP.ICs_Type, "Uniform") == 0) {
          IP.i_ICs = IC_UNIFORM;
       } else if (strcmp(IP.ICs_Type, "Restart") == 0) {
          IP.i_ICs = IC_RESTART;
       } else {
          IP.i_ICs = IC_UNIFORM;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Grid_Type") == 0) {
       i_command = 6;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Grid_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Grid_Type, "Cartesian") == 0) {
          IP.i_Grid = GRID_CARTESIAN_UNIFORM;
          IP.Box_Width = ONE;
          IP.Box_Height = ONE;
       } else if (strcmp(IP.Grid_Type, "Square") == 0) {
          IP.i_Grid = GRID_SQUARE;
          IP.Box_Width = ONE;
          IP.Box_Height = ONE;
       } else if (strcmp(IP.Grid_Type, "Rectangular_Box") == 0) {
          IP.i_Grid = GRID_RECTANGULAR_BOX;
          IP.Box_Width = ONE;
          IP.Box_Height = ONE;
       } else if (strcmp(IP.Grid_Type, "Flat_Plate") == 0) {
          IP.i_Grid = GRID_FLAT_PLATE;
          IP.Plate_Length = ONE;
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
          IP.Grain_Length = 0.835;
          IP.Grain_Radius = 0.020;
          IP.Grain_To_Throat_Length = 0.05;
          IP.Nozzle_Length = 0.150;
          IP.Nozzle_Radius_Exit = 0.030;
          IP.Nozzle_Radius_Throat = 0.010;
       } else if (strcmp(IP.Grid_Type, "Circular_Cylinder") == 0) {
          IP.i_Grid = GRID_CIRCULAR_CYLINDER;
          IP.Cylinder_Radius = ONE;
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
       } else if (strcmp(IP.Grid_Type, "Wedge") == 0) {
	  IP.i_Grid = GRID_WEDGE;
	  IP.Wedge_Angle = 25.0;
	  IP.Wedge_Length = HALF;
       } else if (strcmp(IP.Grid_Type, "Unsteady_Blunt_Body") == 0) {
	  IP.i_Grid = GRID_UNSTEADY_BLUNT_BODY;
	  IP.Blunt_Body_Radius = ONE;
	  IP.Blunt_Body_Mach_Number = TWO;
       } else if (strcmp(IP.Grid_Type,"Bump_Channel_Flow") == 0) {
          IP.i_Grid = GRID_BUMP_CHANNEL_FLOW;
       } else if (strcmp(IP.Grid_Type, "ICEMCFD") == 0) {
          IP.i_Grid = GRID_ICEMCFD;
       } else if (strcmp(IP.Grid_Type, "Read_From_Definition_File") == 0) {
          IP.i_Grid = GRID_READ_FROM_DEFINITION_FILE;
       } else if (strcmp(IP.Grid_Type, "Read_From_Data_File") == 0) {
          IP.i_Grid = GRID_READ_FROM_GRID_DATA_FILE;
       } else {
          IP.i_Grid = GRID_SQUARE;
          IP.Box_Width = ONE;
          IP.Box_Height = ONE;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Output_File_Name") == 0) {
       i_command = 7;
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
       i_command = 8;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Grid_File_Name, 
              IP.Next_Control_Parameter);
       strcat(IP.Grid_File_Name, ".grid");
       strcpy(IP.Grid_Definition_File_Name, 
              IP.Next_Control_Parameter);
       strcat(IP.Grid_Definition_File_Name, ".griddef");

    } else if (strcmp(IP.Next_Control_Parameter, "Restart_File_Name") == 0) {
       i_command = 9;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Restart_File_Name, 
              IP.Next_Control_Parameter);
       strcat(IP.Restart_File_Name, ".soln");

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Cells_Idir") == 0) {
       i_command = 10;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Cells_Idir;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Cells_Idir < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Cells_Jdir") == 0) {
       i_command = 11;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Cells_Jdir;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Cells_Jdir < 1) i_command = INVALID_INPUT_VALUE;

   } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Ghost_Cells") == 0) {
       i_command = 100;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Ghost_Cells;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Ghost_Cells < 2) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Blocks_Idir") == 0) {
       i_command = 12;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Blocks_Idir;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Blocks_Idir < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Blocks_Jdir") == 0) {
       i_command = 13;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Blocks_Jdir;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Blocks_Jdir < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Time_Accurate") == 0) {
       i_command = 14;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Time_Accurate;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Time_Accurate != 0 &&
           IP.Time_Accurate != 1) IP.Time_Accurate = 0;
       if (IP.Time_Accurate) {
          IP.Local_Time_Stepping = GLOBAL_TIME_STEPPING;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Local_Time_Stepping") == 0) {
       i_command = 15;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Local_Time_Stepping;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Local_Time_Stepping != GLOBAL_TIME_STEPPING &&
           IP.Local_Time_Stepping != SCALAR_LOCAL_TIME_STEPPING &&
	   IP.Local_Time_Stepping != MATRIX_LOCAL_TIME_STEPPING ) 
         IP.Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;

    } else if (strcmp(IP.Next_Control_Parameter, "Maximum_Number_of_Time_Steps") == 0) {
       i_command = 16;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Maximum_Number_of_Time_Steps;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Maximum_Number_of_Time_Steps < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "N_Stage") == 0) {
       i_command = 17;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.N_Stage;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.N_Stage < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "CFL_Number") == 0) {
       i_command = 18;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.CFL_Number;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.CFL_Number <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Box_Width") == 0) {
       i_command = 19;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Box_Width;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Box_Width <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Box_Height") == 0) {
       i_command = 20;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Box_Height;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Box_Height <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Plate_Length") == 0) {
       i_command = 21;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Plate_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Plate_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Pipe_Length") == 0) {
       i_command = 22;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Pipe_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Pipe_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Pipe_Radius") == 0) {
       i_command = 23;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Pipe_Radius;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Pipe_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Blunt_Body_Radius") == 0) {
       i_command = 24;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Blunt_Body_Radius;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Blunt_Body_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Blunt_Body_Mach_Number") == 0) {
       i_command = 25;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Blunt_Body_Mach_Number;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Blunt_Body_Mach_Number <= ONE) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Cylinder_Radius") == 0) {
       i_command = 26;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Cylinder_Radius;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Cylinder_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Ellipse_Length_X_Axis") == 0) {
       i_command = 27;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Ellipse_Length_X_Axis;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Ellipse_Length_X_Axis <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Ellipse_Length_Y_Axis") == 0) {
       i_command = 28;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Ellipse_Length_Y_Axis;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Ellipse_Length_Y_Axis <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Chord_Length") == 0) {
       i_command = 29;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Chord_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Chord_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "NACA_Aerofoil_Type") == 0) {
       i_command = 30;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.NACA_Aerofoil_Type, 
              IP.Next_Control_Parameter);
       if (strlen(IP.NACA_Aerofoil_Type) != 4 &&
           strlen(IP.NACA_Aerofoil_Type) != 5) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Orifice_Radius") == 0) {
       i_command = 31;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Orifice_Radius;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Orifice_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Wedge_Angle") == 0) {
       i_command = 31;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Wedge_Angle;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Wedge_Angle <= ZERO) i_command = INVALID_INPUT_VALUE;
    
    } else if (strcmp(IP.Next_Control_Parameter, "Wedge_Length") == 0) {
       i_command = 31;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Wedge_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Wedge_Length <= ZERO) i_command = INVALID_INPUT_VALUE;
    
    } else if (strcmp(IP.Next_Control_Parameter, "Grain_Length") == 0) {
       i_command = 32;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Grain_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Grain_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Grain_Radius") == 0) {
       i_command = 33;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Grain_Radius;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Grain_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Grain_To_Throat_Length") == 0) {
       i_command = 34;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Grain_To_Throat_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Grain_To_Throat_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Nozzle_Length") == 0) {
       i_command = 35;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Nozzle_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Nozzle_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Nozzle_Radius_Exit") == 0) {
       i_command = 36;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Nozzle_Radius_Exit;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Nozzle_Radius_Exit <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Nozzle_Radius_Throat") == 0) {
       i_command = 37;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Nozzle_Radius_Throat;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Nozzle_Radius_Throat <= ZERO) i_command = INVALID_INPUT_VALUE;



       //------------------------ Rte2D Specific -------------------------//

    } else if (strcmp(IP.Next_Control_Parameter, "Difference_Scheme") == 0) {
       i_command = 4;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.SpaceMarch_Scheme, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.SpaceMarch_Scheme, "Upwind") == 0) {
          IP.i_SpaceMarch_Scheme = SPACE_MARCH_UPWIND;
       } else if (strcmp(IP.SpaceMarch_Scheme, "CLAM") == 0) {
          IP.i_SpaceMarch_Scheme = SPACE_MARCH_CLAM;
       } else if (strcmp(IP.SpaceMarch_Scheme, "GM") == 0) {
          IP.i_SpaceMarch_Scheme = SPACE_MARCH_GM;
       } else if (strcmp(IP.SpaceMarch_Scheme, "Central") == 0) {
          IP.i_SpaceMarch_Scheme = SPACE_MARCH_CENTRAL;
       } else if (strcmp(IP.SpaceMarch_Scheme, "Exponential") == 0) {
          IP.i_SpaceMarch_Scheme = SPACE_MARCH_EXPONENTIAL;
       } else {
          IP.i_SpaceMarch_Scheme = SPACE_MARCH_UPWIND;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "RTE_Solver") == 0) {
       i_command = 999;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.RTE_Solver, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.RTE_Solver, "FVM") == 0) {
          IP.i_RTE_Solver = RTE2D_SOLVER_FVM;
       } else if (strcmp(IP.RTE_Solver, "DOM") == 0) {
          IP.i_RTE_Solver = RTE2D_SOLVER_DOM;
       } else {
          IP.i_RTE_Solver = RTE2D_SOLVER_FVM;
       } /* endif */
       IP.Uo.RTE_Type = IP.i_RTE_Solver;

    } else if (strcmp(IP.Next_Control_Parameter, "DOM_Quadrature") == 0) {
       i_command = 999;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.DOM_Quadrature, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.DOM_Quadrature, "S2") == 0) {
          IP.i_DOM_Quadrature = DOM_S2;
       } else if (strcmp(IP.DOM_Quadrature, "S4") == 0) {
          IP.i_DOM_Quadrature = DOM_S4;
       } else if (strcmp(IP.DOM_Quadrature, "S6") == 0) {
          IP.i_DOM_Quadrature = DOM_S6;
       } else if (strcmp(IP.DOM_Quadrature, "S8") == 0) {
          IP.i_DOM_Quadrature = DOM_S8;
       } else if (strcmp(IP.DOM_Quadrature, "S12") == 0) {
          IP.i_DOM_Quadrature = DOM_S12;
       } else if (strcmp(IP.DOM_Quadrature, "T3") == 0) {
          IP.i_DOM_Quadrature = DOM_T3;
       } else {
          IP.i_DOM_Quadrature = DOM_S2;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Absorbsion_Coefficient") == 0) {

       i_command = 38;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.AbsorbsionCoef;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.AbsorbsionCoef < ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } else {
	 IP.Uo.SetAbsorbsion( IP.AbsorbsionCoef );
       } 

    } else if (strcmp(IP.Next_Control_Parameter, "Scattering_Coefficient") == 0) {

       i_command = 39;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.ScatteringCoef;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.ScatteringCoef < ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } else {
	 IP.Uo.SetScattering( IP.ScatteringCoef );
       } 

    } else if (strcmp(IP.Next_Control_Parameter, "Scattering_Function") == 0) {

      i_command = 40;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.ScatteringFunc, IP.Next_Control_Parameter);
       if (strcmp(IP.ScatteringFunc, "Isotropic") == 0) {
	 IP.i_ScatteringFunc = SCATTER_ISO;
       } else if (strcmp(IP.ScatteringFunc, "F1") == 0) {
	 IP.i_ScatteringFunc = SCATTER_F1;
       } else if (strcmp(IP.ScatteringFunc, "F2") == 0) {
	 IP.i_ScatteringFunc = SCATTER_F2;
       } else if (strcmp(IP.ScatteringFunc, "F3") == 0) {
	 IP.i_ScatteringFunc = SCATTER_F3;
       } else if (strcmp(IP.ScatteringFunc, "B1") == 0) {
	 IP.i_ScatteringFunc = SCATTER_B1;
       } else if (strcmp(IP.ScatteringFunc, "B2") == 0) {
	 IP.i_ScatteringFunc = SCATTER_B2;
       } else {
	 i_command = INVALID_INPUT_VALUE;
       } /* endif */


    } else if (strcmp(IP.Next_Control_Parameter, "Gas_Temperature") == 0) {

       i_command = 41;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Temperature;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Temperature <= ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } else {
	 IP.Uo.SetBlackbody( Ib(IP.Temperature) );
       } 

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Angles_Mdir") == 0) {

       i_command = 42;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Angles_Mdir;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Angles_Mdir <= 0) {
          i_command = INVALID_INPUT_VALUE;
       } else {
	 Rte2D_State::SetGas( );
	 Rte2D_State::SetDirs( IP.Number_of_Angles_Mdir, 
			       IP.Number_of_Angles_Ldir, 
			       IP.i_DOM_Quadrature,
			       IP.Axisymmetric );
	 IP.Uo.Allocate();
	 IP.Uo.Zero();
	 IP.Uo.SetBlackbody( Ib(IP.Temperature) );
	 IP.Uo.SetAbsorbsion( IP.AbsorbsionCoef );
	 IP.Uo.SetScattering( IP.ScatteringCoef);
	 Rte2D_State :: SetupPhase( IP.i_ScatteringFunc );
      } 

       
    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Angles_Ldir") == 0) {

       i_command = 43;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Angles_Ldir;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Angles_Ldir <= 0) {
          i_command = INVALID_INPUT_VALUE;
       } else {
	 Rte2D_State::SetGas( );
	 Rte2D_State::SetDirs( IP.Number_of_Angles_Mdir, 
			       IP.Number_of_Angles_Ldir, 
			       IP.i_DOM_Quadrature,
			       IP.Axisymmetric );
	 IP.Uo.Allocate();
	 IP.Uo.Zero();
	 IP.Uo.SetBlackbody( Ib(IP.Temperature) );
	 IP.Uo.SetAbsorbsion( IP.AbsorbsionCoef );
	 IP.Uo.SetScattering( IP.ScatteringCoef );
 	 Rte2D_State :: SetupPhase( IP.i_ScatteringFunc );
      } 

    } else if (strcmp(IP.Next_Control_Parameter, "Wall_Temperature") == 0) {
      
      i_command = 44;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.NorthWallTemp; if (IP.NorthWallTemp <= ZERO) i_command = INVALID_INPUT_VALUE;
      IP.Input_File >> IP.SouthWallTemp; if (IP.SouthWallTemp <= ZERO) i_command = INVALID_INPUT_VALUE;
      IP.Input_File >> IP.EastWallTemp; if (IP.EastWallTemp <= ZERO) i_command = INVALID_INPUT_VALUE;
      IP.Input_File >> IP.WestWallTemp; if (IP.WestWallTemp <= ZERO) i_command = INVALID_INPUT_VALUE;
      IP.Input_File.getline(buffer, sizeof(buffer));


    } else if (strcmp(IP.Next_Control_Parameter, "Wall_Emissivity") == 0) {
      
      i_command = 45;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.NorthWallEmiss; if (IP.NorthWallEmiss <= ZERO) i_command = INVALID_INPUT_VALUE;
      IP.Input_File >> IP.SouthWallEmiss; if (IP.SouthWallEmiss <= ZERO) i_command = INVALID_INPUT_VALUE;
      IP.Input_File >> IP.EastWallEmiss; if (IP.EastWallEmiss <= ZERO) i_command = INVALID_INPUT_VALUE;
      IP.Input_File >> IP.WestWallEmiss; if (IP.WestWallEmiss <= ZERO) i_command = INVALID_INPUT_VALUE;
      IP.Input_File.getline(buffer, sizeof(buffer));


      //---------------------- End Rte2D Specific -----------------------//




    } else if (strcmp(IP.Next_Control_Parameter, "Time_Max") == 0) {
       i_command = 43;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Time_Max;
       IP.Input_File.getline(buffer, sizeof(buffer));
       IP.Time_Max = IP.Time_Max/THOUSAND;
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
       i_command = 50;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Flow_Geometry_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Flow_Geometry_Type, "Planar") == 0) {
          IP.Axisymmetric = 0;
       } else if (strcmp(IP.Flow_Geometry_Type, "Axisymmetric") == 0) {
	 IP.Axisymmetric = AXISYMMETRIC_Y;
       } else if (strcmp(IP.Flow_Geometry_Type, "Axisymmetric-x") == 0) {
	 IP.Axisymmetric = AXISYMMETRIC_X;
       } else if (strcmp(IP.Flow_Geometry_Type, "Axisymmetric-y") == 0) {
	 IP.Axisymmetric = AXISYMMETRIC_Y;
       } else {
          IP.Axisymmetric = PLANAR;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Restart_Solution_Save_Frequency") == 0) {
       i_command = 51;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Restart_Solution_Save_Frequency;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Restart_Solution_Save_Frequency < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "ICEMCFD_Topology_File") == 0) {
       i_command = 52;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.ICEMCFD_FileNames[0], IP.Next_Control_Parameter);

    } else if (strcmp(IP.Next_Control_Parameter, "ICEMCFD_Family_Boco_File") == 0) {
       i_command = 53;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.ICEMCFD_FileNames[1], IP.Next_Control_Parameter);

    } else if (strcmp(IP.Next_Control_Parameter, "ICEMCFD_Family_Topo_File") == 0) {
       i_command = 54;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.ICEMCFD_FileNames[2], IP.Next_Control_Parameter);

    } else if (strcmp(IP.Next_Control_Parameter, "X_Shift") == 0) {
       i_command = 55;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.X_Shift;
       IP.Input_File.setf(ios::skipws);
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "X_Scale") == 0) {
       i_command = 56;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.X_Scale;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "X_Rotate") == 0) {
       i_command = 57;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.X_Rotate;
       IP.Input_File.getline(buffer, sizeof(buffer));

       /****************************************************************/
       /***************** NKS PARAMETERS *******************************/  
       /****************************************************************/

    } else if (strcmp(IP.Next_Control_Parameter, "GMRES_Restart") == 0) {
      // GMRES restart:
      i_command = 58;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.NKS_IP.GMRES_Restart;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.NKS_IP.GMRES_Restart < 0) i_command = INVALID_INPUT_VALUE;
      
    } else if (strcmp(IP.Next_Control_Parameter, "GMRES_Overlap") == 0) {
      // GMRES overlap:
      i_command = 59;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.NKS_IP.GMRES_Overlap;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.NKS_IP.GMRES_Overlap < 0) i_command = INVALID_INPUT_VALUE;
      
    } else if (strcmp(IP.Next_Control_Parameter, "GMRES_Tolerance") == 0) {
      // GMRES tolerance:
      i_command = 60;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.NKS_IP.GMRES_Tolerance;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.NKS_IP.GMRES_Tolerance <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Overall_Tolerance") == 0) {
      // GMRES tolerance:
      i_command = 61;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.NKS_IP.Overall_Tolerance;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.NKS_IP.Overall_Tolerance <= ZERO) i_command = INVALID_INPUT_VALUE;
      
    } else if (strcmp(IP.Next_Control_Parameter, "GMRES_Block_Preconditioner") == 0) {
      // GMRES Precondtioner
      i_command = 62;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter, "ILUK") == 0) {
	IP.NKS_IP.GMRES_Block_Preconditioner = Block_ILUK;
      } else if (strcmp(IP.Next_Control_Parameter, "Diagonal") == 0) {
	IP.NKS_IP.GMRES_Block_Preconditioner = Block_Jacobi;
      } else {
	IP.NKS_IP.GMRES_Block_Preconditioner = Block_ILUK;
      }

    } else if (strcmp(IP.Next_Control_Parameter, "GMRES_ILUK_Level_of_Fill") == 0) {
      i_command = 63;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.NKS_IP.GMRES_ILUK_Level_of_Fill;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.NKS_IP.GMRES_ILUK_Level_of_Fill < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Jacobian_Order") == 0) {
      i_command = 631;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter, "Source_Terms_Only") == 0) {
	IP.NKS_IP.Jacobian_Order = SOURCE_TERMS_ONLY;
      } else if (strcmp(IP.Next_Control_Parameter, "First_Order_Inviscid_HLLE") == 0) {
	IP.NKS_IP.Jacobian_Order = FIRST_ORDER_INVISCID_HLLE;
      } else if (strcmp(IP.Next_Control_Parameter, "First_Order_Inviscid_Roe") == 0) {
	IP.NKS_IP.Jacobian_Order = FIRST_ORDER_INVISCID_ROE;
      } else if (strcmp(IP.Next_Control_Parameter, "Second_Order_Diamond_Path_with_HLLE") == 0) {
	IP.NKS_IP.Jacobian_Order = SECOND_ORDER_DIAMOND_WITH_HLLE;
      } else if (strcmp(IP.Next_Control_Parameter, "Second_Order_Diamond_Path_with_Roe") == 0) {
	IP.NKS_IP.Jacobian_Order = SECOND_ORDER_DIAMOND_WITH_ROE;	
      } else {
	IP.NKS_IP.Jacobian_Order = FIRST_ORDER_INVISCID_HLLE;
      }

    } else if (strcmp(IP.Next_Control_Parameter, "Maximum_Number_of_NKS_Iterations") == 0) {
      i_command = 64;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.NKS_IP.Maximum_Number_of_NKS_Iterations;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.NKS_IP.Maximum_Number_of_NKS_Iterations < 0) i_command = INVALID_INPUT_VALUE;
     
    } else if (strcmp(IP.Next_Control_Parameter, "Maximum_Number_of_GMRES_Iterations") == 0) {
      i_command = 65;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.NKS_IP.Maximum_Number_of_GMRES_Iterations;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.NKS_IP.Maximum_Number_of_GMRES_Iterations < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Finite_Time_Step") == 0) {
      i_command = 66; 
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter, "ON") == 0) {
	IP.NKS_IP.Finite_Time_Step = true;
      } else if (strcmp(IP.Next_Control_Parameter, "OFF") == 0) {
	IP.NKS_IP.Finite_Time_Step = false;
      } else {
	IP.NKS_IP.Finite_Time_Step = true;
      }

    } else if (strcmp(IP.Next_Control_Parameter, "Normalization") == 0) {
      i_command = 67;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter, "ON") == 0) {
	IP.NKS_IP.Normalization = true;
      } else if (strcmp(IP.Next_Control_Parameter, "OFF") == 0) {
	IP.NKS_IP.Normalization = false;
      } else {
	IP.NKS_IP.Normalization = true;
      }

    } else if (strcmp(IP.Next_Control_Parameter, "Finite_Time_Step_Initial_CFL") == 0) {
      i_command = 70;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.NKS_IP.Finite_Time_Step_Initial_CFL;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.NKS_IP.Finite_Time_Step_Initial_CFL < ZERO) i_command = INVALID_INPUT_VALUE;
     
      /****************************************************************/
      /********************* END OF NKS SPECIFIC **********************/
      /****************************************************************/

    } else if (strcmp(IP.Next_Control_Parameter, "Freeze_Limiter") == 0) {
      // Freeze_Limiter:
      i_command = 68;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Freeze_Limiter;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Freeze_Limiter < 0) i_command = INVALID_INPUT_VALUE;
      
    } else if (strcmp(IP.Next_Control_Parameter, "Freeze_Limiter_Residual_Level") == 0) {
      // Freeze_Limiter_Residual_Level:
      i_command = 69;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Freeze_Limiter_Residual_Level;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Freeze_Limiter_Residual_Level < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "AMR") == 0) {
      i_command = 71;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.AMR;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.AMR < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "AMR_Frequency") == 0) {
      i_command = 72;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.AMR_Frequency;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.AMR_Frequency < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Initial_Mesh_Refinements") == 0) {
      i_command = 73;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Number_of_Initial_Mesh_Refinements;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Number_of_Initial_Mesh_Refinements < 0) IP.Number_of_Initial_Mesh_Refinements = 0;

    } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Uniform_Mesh_Refinements") == 0) {
      i_command = 74;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Number_of_Uniform_Mesh_Refinements;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Number_of_Uniform_Mesh_Refinements < 0) IP.Number_of_Uniform_Mesh_Refinements = 0;

    } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Boundary_Mesh_Refinements") == 0) {
      i_command = 75;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Number_of_Boundary_Mesh_Refinements;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Number_of_Boundary_Mesh_Refinements < 0) IP.Number_of_Boundary_Mesh_Refinements = 0;

    } else if (strcmp(IP.Next_Control_Parameter,"Maximum_Refinement_Level") == 0) {
      i_command = 76;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Maximum_Refinement_Level;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Maximum_Refinement_Level < 1) IP.Maximum_Refinement_Level = 1;

    } else if (strcmp(IP.Next_Control_Parameter,"Minimum_Refinement_Level") == 0) {
      i_command = 77;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Minimum_Refinement_Level;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Minimum_Refinement_Level < 1) IP.Minimum_Refinement_Level = 1;

    } else if (strcmp(IP.Next_Control_Parameter, "Threshold_for_Refinement") == 0) {
       i_command = 78;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Threshold_for_Refinement;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Threshold_for_Refinement <= ZERO ||
           IP.Threshold_for_Refinement > ONE) IP.Threshold_for_Refinement = 0.50;

    } else if (strcmp(IP.Next_Control_Parameter, "Threshold_for_Coarsening") == 0) {
       i_command = 79;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Threshold_for_Coarsening;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Threshold_for_Coarsening < ZERO ||
           IP.Threshold_for_Coarsening >= ONE) IP.Threshold_for_Coarsening = 0.10;

    } else if (strcmp(IP.Next_Control_Parameter, "Residual_Smoothing_Epsilon") == 0) {
      i_command = 80;
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
      i_command = 81;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Residual_Smoothing_Gauss_Seidel_Iterations;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Residual_Smoothing_Gauss_Seidel_Iterations < 0) {
         IP.Residual_Smoothing_Gauss_Seidel_Iterations = 0;
      } /* endif */

//     } else if (strcmp(IP.Next_Control_Parameter, "Mr_Min_Factor") == 0) {
//       i_command = 82;
//       IP.Line_Number = IP.Line_Number + 1;
//       IP.Input_File >> IP.Mr_Min_Factor;
//       IP.Input_File.getline(buffer, sizeof(buffer));
//       if (IP.Mr_Min_Factor <= ZERO) {
//          IP.Mr_Min_Factor = ZERO;
//       } /* endif */
//       IP.Wo.Mr_min = IP.Mr_Min_Factor*IP.Mach_Number;
//       IP.Uo.Mr_min = IP.Mr_Min_Factor*IP.Mach_Number;

//     } else if (strcmp(IP.Next_Control_Parameter, "Morton") == 0) {
//       i_command = 83;
//       IP.Line_Number = IP.Line_Number + 1;
//       IP.Input_File >> IP.Morton;
//       IP.Input_File.getline(buffer, sizeof(buffer));
//       if (IP.Morton < 0) i_command = INVALID_INPUT_VALUE;

//     } else if (strcmp(IP.Next_Control_Parameter, "Morton_Reordering_Frequency") == 0) {
//       i_command = 84;
//       IP.Line_Number = IP.Line_Number + 1;
//       IP.Input_File >> IP.Morton_Reordering_Frequency;
//       IP.Input_File.getline(buffer, sizeof(buffer));
//       if (IP.Morton_Reordering_Frequency < 0) i_command = INVALID_INPUT_VALUE;

    /**************************************
     * Multigrid Related Input Parameters *
     **************************************/

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

    } else if (strcmp(IP.Next_Control_Parameter, "Convergence_Residual_Level") == 0) {
       i_command = 208;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Multigrid_IP.Convergence_Residual_Level;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Multigrid_IP.Convergence_Residual_Level < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Multigrid_Smoothing_Type") == 0) {
       i_command = 209;
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
           IP.Multigrid_IP.i_Smoothing = TIME_STEPPING_EXPLICIT_EULER;
           IP.N_Stage = 1;
       } /* endif */

    /******************************************
     * End of Multigrid related Input parsing *
     ******************************************/

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

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Quasi3D") == 0) {
       i_command = WRITE_OUTPUT_QUASI3D_CODE;
    
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

    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Black_Enclosure") == 0) {
       i_command = WRITE_OUTPUT_BLACK_ENCLOSURE_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Refine_Grid") == 0) {
       i_command = REFINE_GRID_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Morton_Ordering") == 0) {
       i_command = MORTON_ORDERING_CODE;

    } else if (IP.Next_Control_Parameter[0] == '#') {
       i_command = COMMENT_CODE;

    } else {
       i_command = INVALID_INPUT_CODE;

    } /* endif */

    /* Return the parser command type indicator. */

    return (i_command);

}

/********************************************************
 * Routine: Process_Input_Control_Parameter_File        *
 *                                                      *
 * Reads, parses, and executes the list of input        *
 * control parameters from the standard input file.     *
 *                                                      *
 ********************************************************/
int Process_Input_Control_Parameter_File(Rte2D_Input_Parameters &Input_Parameters,
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
    error_flag = Input_Parameters.Input_File.bad();

    if (error_flag) {
       cout << "\n Rte2D ERROR: Unable to open Rte2D input data file.\n";
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
          break;
       } else if (Command_Flag == INVALID_INPUT_CODE ||
                  Command_Flag == INVALID_INPUT_VALUE) {
          line_number = -line_number;
          cout << "\n Rte2D ERROR: Error reading Rte2D data at line #"
               << -line_number  << " of input data file.\n";
          error_flag = line_number;
          break;
       } /* endif */
    } /* endwhile */

    /* Perform consistency checks on Input_Parameters */

    if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {
      error_flag = Check_Input_Parameters<Rte2D_Input_Parameters>(Input_Parameters);
      if (error_flag) {
	cout << "\n AdvectDiffuse2D ERROR: Input Parameters consistency check failure\n";
	return (error_flag);
      }
    }

    /* Initial processing of input control parameters complete.  
       Return the error indicator flag. */

    return (error_flag);

}
