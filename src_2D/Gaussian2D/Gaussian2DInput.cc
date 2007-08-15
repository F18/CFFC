/* Gaussian2DInput.cc:  Subroutines for 2D Gaussian Input Classes. */

/* Include 2D Gaussian input parameter header file. */

#ifndef _GAUSSIAN2D_INPUT_INCLUDED
#include "Gaussian2DInput.h"
#endif // _GAUSSIAN2D_INPUT_INCLUDED

/*************************************************************
 * GAUSSIAN2D_Input_Parameters -- External subroutines.      *
 *************************************************************/

/********************************************************
 * Routine: Open_Input_File                             *
 *                                                      *
 * Opens the appropriate input data file.               *
 *                                                      *
 ********************************************************/
void Open_Input_File(Gaussian2D_Input_Parameters &IP) {

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
void Close_Input_File(Gaussian2D_Input_Parameters &IP) {

    IP.Input_File.unsetf(ios::skipws);
    IP.Input_File.close();

}

/********************************************************
 * Routine: Set_Default_Input_Parameters                *
 *                                                      *
 * Assigns default values to the input parameters.      *
 *                                                      *
 ********************************************************/
void Set_Default_Input_Parameters(Gaussian2D_Input_Parameters &IP) {

    int i;
    char *string_ptr;

    string_ptr = "Gaussian2D.in";
    strcpy(IP.Input_File_Name, string_ptr);

    string_ptr = "Explicit_Euler";
    strcpy(IP.Time_Integration_Type, string_ptr);
    IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
    IP.Time_Accurate = 0;
    IP.Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
    IP.Maximum_Number_of_Time_Steps = 100;
    //IP.Maximum_Number_of_NKS_Iterations = 0;
    //IP.Maximum_Number_of_GMRES_Iterations = 0;
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

    string_ptr = "Roe";
    strcpy(IP.Flux_Function_Type, string_ptr);
    IP.i_Flux_Function = FLUX_FUNCTION_ROE;

    string_ptr = "Uniform";
    strcpy(IP.ICs_Type, string_ptr);
    IP.i_ICs = IC_UNIFORM;

    string_ptr = "AIR";
    strcpy(IP.Gas_Type, string_ptr);
    IP.Wo.setgas(IP.Gas_Type);
    IP.Wo = Gaussian2D_W_STDATM;
    IP.Pressure = IP.Wo.pressure();
    IP.Temperature = IP.Wo.T();
    IP.Mach_Number = 0.80;
    IP.Flow_Angle = ZERO;
    IP.Wo.v.x = IP.Mach_Number*IP.Wo.sound()*cos(TWO*PI*IP.Flow_Angle/360.00);
    IP.Wo.v.y = IP.Mach_Number*IP.Wo.sound()*sin(TWO*PI*IP.Flow_Angle/360.00);
    IP.Uo.setgas(IP.Gas_Type);
    IP.Uo = U(IP.Wo);
    IP.W1 = IP.Wo;
    IP.W2 = IP.Wo;

    string_ptr = "Planar";
    strcpy(IP.Flow_Geometry_Type, string_ptr);
    IP.Axisymmetric = 0;

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

    // Boundary Conditions
    IP.alpha = 1.0;
    IP.Ramp_by_Mach_Number = 0.0;
    IP.Number_of_Time_Steps_to_Ramp = 0; 
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
    IP.Cylinder_Radius2 = 32.0;
    IP.Ellipse_Length_X_Axis = TWO;
    IP.Ellipse_Length_Y_Axis = HALF;
    IP.Chord_Length = ONE;
    IP.Orifice_Radius = ONE;
    IP.Inner_Streamline_Number = 0.80;
    IP.Outer_Streamline_Number = 0.40;
    IP.Isotach_Line = 0.30;
    IP.Wedge_Angle = 25.0;
    IP.Wedge_Length = HALF;
    IP.Couette_Plate_Separation = ONE;
    IP.Couette_Plate_Velocity = ZERO;
    IP.Pressure_Drop = ZERO;

    IP.Smooth_Bump = OFF;

    IP.X_Shift = Vector2D_ZERO;
    IP.X_Scale = ONE;
    IP.X_Rotate = ZERO;

    string_ptr = "/home/groth/CFDkit+caboodle/data/NASA_Rotors/R37/";
    strcpy(IP.NASA_Rotor37_Data_Directory, string_ptr);
    string_ptr = "/home/groth/CFDkit+caboodle/data/NASA_Rotors/R67/";
    strcpy(IP.NASA_Rotor67_Data_Directory, string_ptr);
    IP.Rotor_Flow_Type = PEAK_FLOW;
    IP.Rotor_Percent_Span = 50.00;

    IP.ICEMCFD_FileNames = ICEMCFD_get_filenames();

    IP.AMR = 0;
    IP.AMR_Frequency = 100;
    IP.Number_of_Initial_Mesh_Refinements = 0;
    IP.Number_of_Uniform_Mesh_Refinements = 0;
    IP.Number_of_Boundary_Mesh_Refinements = 0;
    IP.Number_of_Interface_Mesh_Refinements = 0;
    IP.Number_of_Bounding_Box_Mesh_Refinements = 0;
    IP.Interface_Refinement_Condition = OFF;
    IP.Maximum_Refinement_Level = 100;
    IP.Minimum_Refinement_Level = 1;

    IP.Threshold_for_Refinement = 0.50;
    IP.Threshold_for_Coarsening = 0.10;

    IP.Number_of_Refinement_Criteria = 3;
    IP.Refinement_Criteria_Gradient_Density = ON;
    IP.Refinement_Criteria_Divergence_Velocity = ON;
    IP.Refinement_Criteria_Curl_Velocity = ON;
    IP.AMR_Xmin = Vector2D_ZERO;
    IP.AMR_Xmax = Vector2D_ZERO;

    // Mesh stretching factor:
    IP.i_Mesh_Stretching = OFF;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_LINEAR;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_LINEAR;
    IP.Mesh_Stretching_Factor_Idir = 1.01;
    IP.Mesh_Stretching_Factor_Jdir = 1.01;

    // Smooth quad block flag:
    IP.i_Smooth_Quad_Block = ON;

    // Embedded boundary parameters:
    IP.Reset_Interface_Motion_Type = OFF;

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

    string_ptr = " ";
    strcpy(IP.Next_Control_Parameter, string_ptr);

    IP.Line_Number = 0;

    IP.Number_of_Processors = CFDkit_MPI::Number_of_Processors;
    IP.Number_of_Blocks_Per_Processor = 10;

    // GMRES restart:
    //IP.GMRES_Restart = 30;
    
    // GMRES tolerance:
    //IP.GMRES_Toler = 1e-5;
    
    // Overall tolerance:
    IP.Overall_Toler = 1e-5;
    
    // GMRES overlap:
    //IP.GMRES_Overlap = 0;

    // GMRES P_Switch:
    //IP.GMRES_P_Switch = 1;

    // GMRES ILUK Level of Fill:
    //IP.GMRES_ILUK_Level_of_Fill = 0;

    // Finite_Time_Step:
    IP.Finite_Time_Step = 0;
    IP.Finite_Time_Step_Initial_CFL = 1;

    // Normalization:
    IP.Normalization = 0;

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
void Broadcast_Input_Parameters(Gaussian2D_Input_Parameters &IP) {

#ifdef _MPI_VERSION
    int i;
    MPI::COMM_WORLD.Bcast(IP.Input_File_Name, 
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Line_Number), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Time_Integration_Type, 
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
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
    //MPI::COMM_WORLD.Bcast(&(IP.Maximum_Number_of_NKS_Iterations), 
    //                      1, 
    //                      MPI::INT, 0);
    //MPI::COMM_WORLD.Bcast(&(IP.Maximum_Number_of_GMRES_Iterations), 
    //                      1, 
    //                      MPI::INT, 0);
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
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Reconstruction), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Limiter_Type, 
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Limiter), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Flux_Function_Type, 
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Flux_Function), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.ICs_Type, 
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Gas_Type, 
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_ICs), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Pressure), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Temperature), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Mach_Number), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Flow_Angle), 
                          1, 
                          MPI::DOUBLE, 0);
    if (!CFDkit_Primary_MPI_Processor()) {
       IP.Wo.setgas(IP.Gas_Type);
       IP.Wo = Gaussian2D_pState(IP.Pressure*IP.Wo.M/(AVOGADRO*BOLTZMANN*IP.Temperature)/THOUSAND, 
                              ZERO, 
                              ZERO, 
                              IP.Pressure);
       IP.Wo.v.x = IP.Mach_Number*IP.Wo.sound()*cos(TWO*PI*IP.Flow_Angle/360.00);
       IP.Wo.v.y = IP.Mach_Number*IP.Wo.sound()*sin(TWO*PI*IP.Flow_Angle/360.00);
       IP.Uo.setgas(IP.Gas_Type);
       IP.Uo = U(IP.Wo);
    } // endif 
    MPI::COMM_WORLD.Bcast(&(IP.W1.d), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W1.v.x), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W1.v.y), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W1.p.xx), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W1.p.xy), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W1.p.yy), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W1.p.zz), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W1.erot), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W2.d), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W2.v.x), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W2.v.y), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W2.p.xx), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W2.p.xy), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W2.p.yy), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W2.p.zz), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W2.erot), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(IP.Flow_Geometry_Type, 
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Axisymmetric), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Grid_Type, 
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.NACA_Aerofoil_Type, 
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
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
    MPI::COMM_WORLD.Bcast(&(IP.Pressure_Drop), 
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
    MPI::COMM_WORLD.Bcast(&(IP.Cylinder_Radius2), 
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
    MPI::COMM_WORLD.Bcast(&(IP.Inner_Streamline_Number),
			  1,
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Outer_Streamline_Number),
			  1,
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Isotach_Line),
			  1,
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Wedge_Angle), 
			  1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Wedge_Length), 
			  1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Couette_Plate_Separation), 
			  1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Couette_Plate_Velocity), 
			  1, 
			  MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Smooth_Bump),
			  1,
			  MPI::INT,0);
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

    //Boundary infromation

    MPI::COMM_WORLD.Bcast(&(IP.alpha), 
                          1, 
                          MPI::DOUBLE, 0);
    if (!CFDkit_Primary_MPI_Processor()) {
      IP.Wo.alpha = IP.alpha;
      IP.Uo.alpha = IP.alpha;
    } // endif 
    MPI::COMM_WORLD.Bcast(&(IP.Ramp_by_Mach_Number), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Time_Steps_to_Ramp), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Boundary_Conditions_Specified,
			  INPUT_PARAMETER_LENGTH_GAUSSIAN2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(&(IP.BCs_Specified),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(IP.BC_North_Type,
			  INPUT_PARAMETER_LENGTH_GAUSSIAN2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(IP.BC_South_Type,
			  INPUT_PARAMETER_LENGTH_GAUSSIAN2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(IP.BC_East_Type,
			  INPUT_PARAMETER_LENGTH_GAUSSIAN2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(IP.BC_West_Type,
			  INPUT_PARAMETER_LENGTH_GAUSSIAN2D,
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
    MPI::COMM_WORLD.Bcast(IP.NASA_Rotor37_Data_Directory, 
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.NASA_Rotor67_Data_Directory, 
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Rotor_Flow_Type), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Rotor_Percent_Span), 
                          1, 
                          MPI::DOUBLE, 0);
    if (!CFDkit_Primary_MPI_Processor()) {
       IP.ICEMCFD_FileNames = new char*[3];
       for (i = 0; i < 3; i++) {
          IP.ICEMCFD_FileNames[i] = new char[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];
       } // endfor 
    } // endif 
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[0], 
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[1], 
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[2], 
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
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
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Interface_Mesh_Refinements),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Bounding_Box_Mesh_Refinements),
			  1,
			  MPI::INT,0);
    MPI::COMM_WORLD.Bcast(&(IP.Interface_Refinement_Condition),
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
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Refinement_Criteria),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Refinement_Criteria_Gradient_Density),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Refinement_Criteria_Divergence_Velocity),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.Refinement_Criteria_Curl_Velocity),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.AMR_Xmin.x),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.AMR_Xmin.y),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.AMR_Xmax.x),
			  1,
			  MPI::DOUBLE,0);
    MPI::COMM_WORLD.Bcast(&(IP.AMR_Xmax.y),
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
    // Interface input parameters:
    IP.Interface_IP.Broadcast_Input_Parameters();
    MPI::COMM_WORLD.Bcast(&(IP.Reset_Interface_Motion_Type),
			  1,
			  MPI::INT,0);
    // File names
    MPI::COMM_WORLD.Bcast(IP.Output_File_Name, 
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Grid_File_Name, 
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Grid_Definition_File_Name, 
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Restart_File_Name, 
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Gnuplot_File_Name, 
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Output_Format_Type, 
                          INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
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

    if (!CFDkit_Primary_MPI_Processor()) {
       IP.Number_of_Processors = CFDkit_MPI::Number_of_Processors;
    } // endif
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Per_Processor), 
                          1, 
                          MPI::INT, 0);
    // GMRES restart:
    //MPI::COMM_WORLD.Bcast(&(IP.GMRES_Restart), 
    //                      1, 
    //                      MPI::INT, 0);
    // GMRES overlap:
    //MPI::COMM_WORLD.Bcast(&(IP.GMRES_Overlap), 
    //			  1, 
    //		  MPI::INT, 0);
    // GMRES tolerance:
    //MPI::COMM_WORLD.Bcast(&(IP.GMRES_Toler), 
    //			  1, 
    //                    MPI::DOUBLE, 0);
    // Overall tolerance:
    MPI::COMM_WORLD.Bcast(&(IP.Overall_Toler), 
			  1, 
                          MPI::DOUBLE, 0);
    // GMRES P_Switch:
    //MPI::COMM_WORLD.Bcast(&(IP.GMRES_P_Switch), 
    //			  1, 
    //			  MPI::INT, 0);

    // GMRES ILUK_Level_of_Fill:
    //MPI::COMM_WORLD.Bcast(&(IP.GMRES_ILUK_Level_of_Fill), 
    //		  1, 
    //		  MPI::INT, 0);

    // Finite Time Step:
    MPI::COMM_WORLD.Bcast(&(IP.Finite_Time_Step), 
			  1, 
			  MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Finite_Time_Step_Initial_CFL), 
			  1, 
			  MPI::DOUBLE, 0);

    // Normalization:
    MPI::COMM_WORLD.Bcast(&(IP.Normalization), 
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
void Broadcast_Input_Parameters(Gaussian2D_Input_Parameters &IP,
                                MPI::Intracomm &Communicator, 
                                const int Source_CPU) {

    int Source_Rank = 0;
    int i;

    Communicator.Bcast(IP.Input_File_Name, 
                       INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.Line_Number), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Time_Integration_Type, 
                       INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
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
    //Communicator.Bcast(&(IP.Maximum_Number_of_NKS_Iterations), 
    //                   1, 
    //                   MPI::INT, Source_Rank);
    //Communicator.Bcast(&(IP.Maximum_Number_of_GMRES_Iterations), 
    //                 1, 
    //                 MPI::INT, Source_Rank);
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
                       INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Reconstruction), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Limiter_Type, 
                       INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Limiter), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Flux_Function_Type, 
                       INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Flux_Function), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.ICs_Type, 
                       INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Gas_Type, 
                       INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_ICs), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Pressure), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Temperature), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Mach_Number), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Flow_Angle), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    if (!(CFDkit_MPI::This_Processor_Number == Source_CPU)) {
       IP.Wo.setgas(IP.Gas_Type);
       IP.Wo = Gaussian2D_pState(IP.Pressure*IP.Wo.M/(AVOGADRO*BOLTZMANN*IP.Temperature)/THOUSAND, 
                              ZERO, 
                              ZERO, 
                              IP.Pressure);
       IP.Wo.v.x = IP.Mach_Number*IP.Wo.sound()*cos(TWO*PI*IP.Flow_Angle/360.00);
       IP.Wo.v.y = IP.Mach_Number*IP.Wo.sound()*sin(TWO*PI*IP.Flow_Angle/360.00);
       IP.Uo.setgas(IP.Gas_Type);
       IP.Uo = U(IP.Wo);
    } // endif
    Communicator.Bcast(&(IP.W1.d), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W1.v.x), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W1.v.y), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W1.p.xx), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W1.p.xy), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W1.p.yy), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W1.p.zz), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W1.erot), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W2.d), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W2.v.x), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W2.v.y), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W2.p.xx), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W2.p.xy), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W2.p.yy), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W2.p.zz), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W2.erot), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(IP.Flow_Geometry_Type, 
                       INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.Axisymmetric), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Grid_Type, 
                       INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.NACA_Aerofoil_Type, 
                       INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
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
    Communicator.Bcast(&(IP.Pressure_Drop), 
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
    Communicator.Bcast(&(IP.Cylinder_Radius2), 
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
    Communicator.Bcast(&(IP.Inner_Streamline_Number),
		       1,
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Outer_Streamline_Number),
		       1,
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Isotach_Line),
		       1,
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Wedge_Angle), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Wedge_Length), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Couette_Plate_Separation), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Couette_Plate_Velocity), 
		       1, 
		       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Smooth_Bump),
		       1,
		       MPI::INT,Source_Rank);
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

    //Boundary information

    Communicator.Bcast(&(IP.alpha), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    if (!(CFDkit_MPI::This_Processor_Number == Source_CPU)) {
      IP.Wo.alpha = IP.alpha;
      IP.Uo.alpha = IP.alpha;
    } // endif
    Communicator.Bcast(&(IP.Ramp_by_Mach_Number), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Time_Steps_to_Ramp), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Boundary_Conditions_Specified,
		       INPUT_PARAMETER_LENGTH_GAUSSIAN2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(&(IP.BCs_Specified),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(IP.BC_North_Type,
		       INPUT_PARAMETER_LENGTH_GAUSSIAN2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(IP.BC_South_Type,
		       INPUT_PARAMETER_LENGTH_GAUSSIAN2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(IP.BC_East_Type,
		       INPUT_PARAMETER_LENGTH_GAUSSIAN2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(IP.BC_West_Type,
		       INPUT_PARAMETER_LENGTH_GAUSSIAN2D,
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
    Communicator.Bcast(IP.NASA_Rotor37_Data_Directory, 
                     INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                     MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.NASA_Rotor67_Data_Directory, 
                     INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                     MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.Rotor_Flow_Type), 
                     1, 
                     MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Rotor_Percent_Span), 
                     1, 
                     MPI::DOUBLE, Source_Rank);
    if (!(CFDkit_MPI::This_Processor_Number == Source_CPU)) {
       IP.ICEMCFD_FileNames = new char*[3];
       for (i = 0; i < 3; i++) {
          IP.ICEMCFD_FileNames[i] = new char[INPUT_PARAMETER_LENGTH_GAUSSIAN2D];
       } // endfor 
    } // endif 
    Communicator.Bcast(IP.ICEMCFD_FileNames[0], 
                       INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.ICEMCFD_FileNames[1], 
                       INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.ICEMCFD_FileNames[2], 
                       INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
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
    Communicator.Bcast(&(IP.Number_of_Interface_Mesh_Refinements),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Number_of_Bounding_Box_Mesh_Refinements),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(&(IP.Interface_Refinement_Condition),
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
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Refinement_Criteria_Gradient_Density),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Refinement_Criteria_Divergence_Velocity),
		     1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.Refinement_Criteria_Curl_Velocity),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.AMR_Xmin.x),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.AMR_Xmin.y),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.AMR_Xmax.x),
		       1,
		       MPI::DOUBLE,Source_Rank);
    Communicator.Bcast(&(IP.AMR_Xmax.y),
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
    // Interface input parameters:
    IP.Interface_IP.Broadcast_Input_Parameters(Communicator,
					       Source_CPU);
    Communicator.Bcast(&(IP.Reset_Interface_Motion_Type),
		       1,
		       MPI::INT,Source_Rank);
    // File Names
    Communicator.Bcast(IP.Output_File_Name, 
                       INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Grid_File_Name, 
                       INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Grid_Definition_File_Name, 
                       INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Restart_File_Name, 
                       INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Gnuplot_File_Name, 
                       INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Output_Format_Type, 
                       INPUT_PARAMETER_LENGTH_GAUSSIAN2D, 
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
    
    if (!(CFDkit_MPI::This_Processor_Number == Source_CPU)) {
       IP.Number_of_Processors = CFDkit_MPI::Number_of_Processors;
    } // endif
    Communicator.Bcast(&(IP.Number_of_Blocks_Per_Processor), 
                       1, 
                       MPI::INT, Source_Rank);
    // GMRES restart:
    //Communicator.Bcast(&(IP.GMRES_Restart), 
    //                  1, 
    //                 MPI::INT, Source_Rank);
    // GMRES overlap:
    //Communicator.Bcast(&(IP.GMRES_Overlap), 
    //                  1, 
    //                 MPI::INT, Source_Rank);
    // GMRES tolerance:
    //Communicator.Bcast(&(IP.GMRES_Toler), 
    //                  1, 
    //                  MPI::DOUBLE, Source_Rank);
    // Overll tolerance:
    Communicator.Bcast(&(IP.Overall_Toler), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    // GMRES P_Switch:
    //Communicator.Bcast(&(IP.GMRES_P_Switch), 
    //                  1, 
    //                 MPI::INT, Source_Rank);
    // GMRES ILUK_Level_of_Fill:
    //Communicator.Bcast(&(IP.GMRES_ILUK_Level_of_Fill), 
    //                  1, 
    //                 MPI::INT, Source_Rank);

    // GMRES Finite_Time_Step:
    Communicator.Bcast(&(IP.Finite_Time_Step), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Finite_Time_Step_Initial_CFL), 
                       1, 
                       MPI::DOUBLE, Source_Rank);

    // GMRES Normalization:
    Communicator.Bcast(&(IP.Normalization), 
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
void Get_Next_Input_Control_Parameter(Gaussian2D_Input_Parameters &IP) {

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
int Parse_Next_Input_Control_Parameter(Gaussian2D_Input_Parameters &IP) {

    int i_command;
    char *string_ptr;
    char buffer[256];
    int tpt, bct;
    Vector2D Xt;

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
	 cout << "Warning multigrid has not been tested with Gaussian2D and will probably not work!" << endl;
	 IP.i_Time_Integration = TIME_STEPPING_MULTIGRID;
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
	  i_command = INVALID_INPUT_VALUE;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Flux_Function_Type") == 0) {
       i_command = 4;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Flux_Function_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Flux_Function_Type, "Roe") == 0) {
	  IP.i_Flux_Function = FLUX_FUNCTION_ROE;
       } else if (strcmp(IP.Flux_Function_Type, "HLLE") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_HLLE;
       } else if (strcmp(IP.Flux_Function_Type, "Roe_MB") == 0) {
	  IP.i_Flux_Function = FLUX_FUNCTION_ROE_MB;
       } else {
	  i_command = INVALID_INPUT_VALUE;
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
       } else if (strcmp(IP.ICs_Type, "Sod") == 0) {
          IP.i_ICs = IC_SOD;
       } else if (strcmp(IP.ICs_Type, "Sod_Xdir") == 0) {
          IP.i_ICs = IC_SOD_XDIR;
       } else if (strcmp(IP.ICs_Type, "Sod_Ydir") == 0) {
          IP.i_ICs = IC_SOD_YDIR;
       } else if (strcmp(IP.ICs_Type, "Groth") == 0) {
          IP.i_ICs = IC_GROTH;
       } else if (strcmp(IP.ICs_Type, "Groth_Xdir") == 0) {
          IP.i_ICs = IC_GROTH_XDIR;
       } else if (strcmp(IP.ICs_Type, "Groth_Ydir") == 0) {
          IP.i_ICs = IC_GROTH_YDIR;
       } else if (strcmp(IP.ICs_Type, "Einfeldt") == 0) {
          IP.i_ICs = IC_EINFELDT;
       } else if (strcmp(IP.ICs_Type, "Einfeldt_Xdir") == 0) {
          IP.i_ICs = IC_EINFELDT_XDIR;
       } else if (strcmp(IP.ICs_Type, "Einfeldt_Ydir") == 0) {
          IP.i_ICs = IC_EINFELDT_YDIR;
       } else if (strcmp(IP.ICs_Type, "Shock_Wave_Xdir") == 0) {
          IP.i_ICs = IC_SHOCK_WAVE_XDIR;
       } else if (strcmp(IP.ICs_Type, "Shock_Wave_Ydir") == 0) {
          IP.i_ICs = IC_SHOCK_WAVE_YDIR;
       } else if (strcmp(IP.ICs_Type, "Shock_Structure_M1_1") == 0) {
          IP.i_ICs = IC_SHOCK_STRUCTURE_M1_1;
	  string_ptr = "A";
	  strcpy(IP.Gas_Type, string_ptr);
       } else if (strcmp(IP.ICs_Type, "Shock_Structure_M1_3") == 0) {
          IP.i_ICs = IC_SHOCK_STRUCTURE_M1_3;
	  string_ptr = "A";
	  strcpy(IP.Gas_Type, string_ptr);
       } else if (strcmp(IP.ICs_Type, "Shock_Structure_M1_5") == 0) {
          IP.i_ICs = IC_SHOCK_STRUCTURE_M1_5;
	  string_ptr = "A";
	  strcpy(IP.Gas_Type, string_ptr);
       } else if (strcmp(IP.ICs_Type, "Shock_Structure_M2_0") == 0) {
          IP.i_ICs = IC_SHOCK_STRUCTURE_M2_0;
	  string_ptr = "A";
	  strcpy(IP.Gas_Type, string_ptr);
       } else if (strcmp(IP.ICs_Type, "Shock_Structure_M10_0") == 0) {
          IP.i_ICs = IC_SHOCK_STRUCTURE_M10_0;
	  string_ptr = "A";
	  strcpy(IP.Gas_Type, string_ptr);
       } else if (strcmp(IP.ICs_Type, "Contact_Surface_Xdir") == 0) {
          IP.i_ICs = IC_CONTACT_SURFACE_XDIR;
       } else if (strcmp(IP.ICs_Type, "Contact_Surface_Ydir") == 0) {
          IP.i_ICs = IC_CONTACT_SURFACE_YDIR;
       } else if (strcmp(IP.ICs_Type, "Rarefaction_Wave_Xdir") == 0) {
          IP.i_ICs = IC_RAREFACTION_WAVE_XDIR;
       } else if (strcmp(IP.ICs_Type, "Rarefaction_Wave_Ydir") == 0) {
          IP.i_ICs = IC_RAREFACTION_WAVE_YDIR;
       } else if (strcmp(IP.ICs_Type, "ShockBox") == 0) {
          IP.i_ICs = IC_SHOCK_BOX;
       } else if (strcmp(IP.ICs_Type, "High_Pressure_Reservoir") == 0) {
          IP.i_ICs = IC_HIGH_PRESSURE_RESERVOIR;
       } else if (strcmp(IP.ICs_Type, "Low_Pressure_Reservoir") == 0) {
          IP.i_ICs = IC_LOW_PRESSURE_RESERVOIR;
       } else if (strcmp(IP.ICs_Type, "Riemann") == 0) {
          IP.i_ICs = IC_RIEMANN;
       } else if (strcmp(IP.ICs_Type, "Riemann_Xdir") == 0) {
          IP.i_ICs = IC_RIEMANN_XDIR;
       } else if (strcmp(IP.ICs_Type, "Riemann_Ydir") == 0) {
          IP.i_ICs = IC_RIEMANN_YDIR;
       } else if (strcmp(IP.ICs_Type, "Wedge_Flow") == 0) {
	  IP.i_ICs = IC_WEDGE_FLOW;
       } else if (strcmp(IP.ICs_Type, "Unsteady_Blunt_Body") == 0) {
          IP.i_ICs = IC_UNSTEADY_BLUNT_BODY;
       } else if (strcmp(IP.ICs_Type, "Pipe") == 0) {
          IP.i_ICs = IC_PIPE;
       } else if (strcmp(IP.ICs_Type, "Couette") == 0) {
          IP.i_ICs = IC_COUETTE;
       } else if (strcmp(IP.ICs_Type, "Flat_Plate") == 0) {
          IP.i_ICs = IC_VISCOUS_FLAT_PLATE;
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
       } else if (strcmp(IP.Grid_Type, "Micro_Cartesian") == 0) {
          IP.i_Grid = GRID_SQUARE;
          IP.Box_Width = 2.0e-6;
          IP.Box_Height = 2.0e-6;
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
       } else if (strcmp(IP.Grid_Type, "Circular_Cylinder") == 0) {
          IP.i_Grid = GRID_CIRCULAR_CYLINDER;
          IP.Cylinder_Radius = ONE;
          IP.Cylinder_Radius2 = 32.0;
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
//       } else if (strcmp(IP.Grid_Type, "NASA_Rotor37") == 0) {
//          IP.i_Grid = GRID_NASA_ROTOR_37;
//          IP.Rotor_Flow_Type = PEAK_FLOW;
//          IP.Rotor_Percent_Span = 50.00;
//          IP.NASA_Rotor37.init(IP.Rotor_Flow_Type,
//                               IP.NASA_Rotor37_Data_Directory);
//          IP.Wo.setgas("AIR");
//          IP.Wo = IP.NASA_Rotor37.getPstateREL_up(IP.Rotor_Percent_Span);
//          IP.Pressure = IP.Wo.p;
//          IP.Temperature = IP.Wo.T();
//          IP.Mach_Number = IP.NASA_Rotor37.getMachREL_up(IP.Rotor_Percent_Span);
//          IP.Flow_Angle = atan2(IP.Wo.v.y, IP.Wo.v.x); 
//          if (IP.Flow_Angle < ZERO) IP.Flow_Angle = TWO*PI + IP.Flow_Angle;
//          IP.Flow_Angle = 180.00*IP.Flow_Angle/PI;
//          IP.W1 = IP.NASA_Rotor37.getPstateREL_down(IP.Rotor_Percent_Span);
//       } else if (strcmp(IP.Grid_Type, "NASA_Rotor67") == 0) {
//	 IP.i_Grid = GRID_NASA_ROTOR_67;
//	 IP.Rotor_Flow_Type = PEAK_FLOW;
//	 IP.Rotor_Percent_Span = 50.00;
//	 IP.NASA_Rotor67.init(IP.Rotor_Flow_Type,
//			      IP.NASA_Rotor67_Data_Directory);
//	 IP.Wo.setgas("AIR");
//	 IP.Wo = IP.NASA_Rotor67.getPstateREL_up(IP.Rotor_Percent_Span);
//	 IP.Pressure = IP.Wo.p;
//	 IP.Temperature = IP.Wo.T();
//	 IP.Mach_Number = IP.NASA_Rotor67.getMachREL_up(IP.Rotor_Percent_Span);
//	 IP.Flow_Angle = atan2(IP.Wo.v.y, IP.Wo.v.x); 
//	 if (IP.Flow_Angle < ZERO) IP.Flow_Angle = TWO*PI + IP.Flow_Angle;
//	 IP.Flow_Angle = 180.00*IP.Flow_Angle/PI;
//	 IP.W1 = IP.NASA_Rotor67.getPstateREL_down(IP.Rotor_Percent_Span);
       } else if (strcmp(IP.Grid_Type, "Bump_Channel_Flow") == 0) {
	 IP.i_Grid = GRID_BUMP_CHANNEL_FLOW;
       } else if (strcmp(IP.Grid_Type,"Non_Smooth_Bump_Channel_Flow") == 0) {
	 IP.i_Grid = GRID_BUMP_CHANNEL_FLOW;
	 IP.Smooth_Bump = OFF;
       } else if (strcmp(IP.Grid_Type,"Smooth_Bump_Channel_Flow") == 0) {
	 IP.i_Grid = GRID_BUMP_CHANNEL_FLOW;
	 IP.Smooth_Bump = ON;
       } else if (strcmp(IP.Grid_Type, "Adiabatic_Flat_Plate") == 0) {
	  IP.i_Grid = GRID_ADIABATIC_FLAT_PLATE;
          IP.Plate_Length = ONE;
       } else if (strcmp(IP.Grid_Type, "Adiabatic_Pipe") == 0) {
	  IP.i_Grid = GRID_ADIABATIC_PIPE;
          IP.Pipe_Length = ONE;
          IP.Pipe_Radius = HALF;
       } else if (strcmp(IP.Grid_Type, "Adiabatic_Circular_Cylinder") == 0) {
	  IP.i_Grid = GRID_ADIABATIC_CIRCULAR_CYLINDER;
          IP.Cylinder_Radius = ONE;
          IP.Cylinder_Radius2 = 32.0;
       } else if (strcmp(IP.Grid_Type, "Adiabatic_Couette") == 0) {
	  IP.i_Grid = GRID_ADIABATIC_COUETTE;
          IP.Couette_Plate_Separation = ONE;
	  IP.Couette_Plate_Velocity = ZERO;
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

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Ghost_Cells") == 0) {
       i_command = 12;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Ghost_Cells;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Ghost_Cells < 1) i_command = INVALID_INPUT_VALUE;

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
	   IP.Local_Time_Stepping != MATRIX_LOCAL_TIME_STEPPING) 
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

    } else if (strcmp(IP.Next_Control_Parameter, "Pressure_Drop") == 0) {
       i_command = 24;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Pressure_Drop;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Pipe_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Blunt_Body_Radius") == 0) {
       i_command = 25;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Blunt_Body_Radius;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Blunt_Body_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Blunt_Body_Mach_Number") == 0) {
       i_command = 26;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Blunt_Body_Mach_Number;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Blunt_Body_Mach_Number <= ONE) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Cylinder_Radius") == 0) {
       i_command = 27;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Cylinder_Radius;
       IP.Cylinder_Radius2 = 32.0*IP.Cylinder_Radius;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Cylinder_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Cylinder_Radius2") == 0) {
       i_command = 27;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Cylinder_Radius2;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Cylinder_Radius2 <= ZERO) i_command = INVALID_INPUT_VALUE;
       if (IP.Cylinder_Radius2 <= IP.Cylinder_Radius) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Ellipse_Length_X_Axis") == 0) {
       i_command = 28;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Ellipse_Length_X_Axis;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Ellipse_Length_X_Axis <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Ellipse_Length_Y_Axis") == 0) {
       i_command = 29;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Ellipse_Length_Y_Axis;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Ellipse_Length_Y_Axis <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Chord_Length") == 0) {
       i_command = 30;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Chord_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Chord_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "NACA_Aerofoil_Type") == 0) {
       i_command = 31;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.NACA_Aerofoil_Type, 
              IP.Next_Control_Parameter);
       if (strlen(IP.NACA_Aerofoil_Type) != 4 &&
           strlen(IP.NACA_Aerofoil_Type) != 5) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Orifice_Radius") == 0) {
       i_command = 32;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Orifice_Radius;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Orifice_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Wedge_Angle") == 0) {
       i_command = 33;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Wedge_Angle;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Wedge_Angle <= ZERO) i_command = INVALID_INPUT_VALUE;
    
    } else if (strcmp(IP.Next_Control_Parameter, "Wedge_Length") == 0) {
       i_command = 34;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Wedge_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Wedge_Length <= ZERO) i_command = INVALID_INPUT_VALUE;
    
    } else if (strcmp(IP.Next_Control_Parameter,"Smooth_Bump") == 0) {
      i_command = 37;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Smooth_Bump = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.Smooth_Bump = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter, "Couette_Plate_Separation") == 0) {
       i_command = 35;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Couette_Plate_Separation;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Couette_Plate_Separation <= ZERO) i_command = INVALID_INPUT_VALUE;
    
    } else if (strcmp(IP.Next_Control_Parameter, "Couette_Plate_Velocity") == 0) {
       i_command = 36;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Couette_Plate_Velocity;
       IP.Input_File.getline(buffer, sizeof(buffer));
    
    } else if (strcmp(IP.Next_Control_Parameter, "Grain_Length") == 0) {
       i_command = 37;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Grain_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Grain_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Grain_Radius") == 0) {
       i_command = 38;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Grain_Radius;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Grain_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Grain_To_Throat_Length") == 0) {
       i_command = 39;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Grain_To_Throat_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Grain_To_Throat_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Nozzle_Length") == 0) {
       i_command = 40;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Nozzle_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Nozzle_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Nozzle_Radius_Exit") == 0) {
       i_command = 41;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Nozzle_Radius_Exit;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Nozzle_Radius_Exit <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Nozzle_Radius_Throat") == 0) {
       i_command = 42;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Nozzle_Radius_Throat;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Nozzle_Radius_Throat <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Alpha") == 0) {
       i_command = 43;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.alpha;
       IP.Wo.alpha = IP.alpha;
       IP.Uo.alpha = IP.alpha;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.alpha < ZERO || IP.alpha > ONE) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Ramp_by_Mach_Number") == 0) {
       i_command = 43;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Ramp_by_Mach_Number;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Time_Steps_to_Ramp") == 0) {
       i_command = 43;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Time_Steps_to_Ramp;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Time_Steps_to_Ramp < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Gas_Type") == 0) {
       i_command = 44;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Gas_Type, 
              IP.Next_Control_Parameter);
       IP.Wo.setgas(IP.Gas_Type);
       IP.Wo = Gaussian2D_pState(IP.Pressure*IP.Wo.M/(AVOGADRO*BOLTZMANN*IP.Temperature)/THOUSAND, 
                              ZERO, 
                              ZERO, 
                              IP.Pressure);
       IP.Wo.v.x = IP.Mach_Number*IP.Wo.sound()*cos(TWO*PI*IP.Flow_Angle/360.00);
       IP.Wo.v.y = IP.Mach_Number*IP.Wo.sound()*sin(TWO*PI*IP.Flow_Angle/360.00);
       IP.Uo.setgas(IP.Gas_Type);
       IP.Uo = U(IP.Wo);

    } else if (strcmp(IP.Next_Control_Parameter, "Mach_Number") == 0) {
       i_command = 45;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Mach_Number;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Mach_Number < ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } else {
          IP.Wo.v.x = IP.Mach_Number*IP.Wo.sound()*cos(TWO*PI*IP.Flow_Angle/360.00);
          IP.Wo.v.y = IP.Mach_Number*IP.Wo.sound()*sin(TWO*PI*IP.Flow_Angle/360.00);
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Flow_Angle") == 0) {
       i_command = 46;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Flow_Angle;
       IP.Input_File.getline(buffer, sizeof(buffer));
       IP.Wo.v.x = IP.Mach_Number*IP.Wo.sound()*cos(TWO*PI*IP.Flow_Angle/360.00);
       IP.Wo.v.y = IP.Mach_Number*IP.Wo.sound()*sin(TWO*PI*IP.Flow_Angle/360.00);

    } else if (strcmp(IP.Next_Control_Parameter, "Pressure") == 0) {
       i_command = 47;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Pressure;
       IP.Input_File.getline(buffer, sizeof(buffer));
       IP.Pressure = IP.Pressure*THOUSAND;
       if (IP.Pressure <= ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } else {
          IP.Wo = Gaussian2D_pState(IP.Pressure*IP.Wo.M/(AVOGADRO*BOLTZMANN*IP.Temperature)/THOUSAND, 
                                 ZERO, 
                                 ZERO, 
                                 IP.Pressure);
          IP.Wo.v.x = IP.Mach_Number*IP.Wo.sound()*cos(TWO*PI*IP.Flow_Angle/360.00);
          IP.Wo.v.y = IP.Mach_Number*IP.Wo.sound()*sin(TWO*PI*IP.Flow_Angle/360.00);
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Temperature") == 0) {
       i_command = 48;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Temperature;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Temperature <= ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } else {
          IP.Wo = Gaussian2D_pState(IP.Pressure*IP.Wo.M/(AVOGADRO*BOLTZMANN*IP.Temperature)/THOUSAND, 
                                 ZERO, 
                                 ZERO, 
                                 IP.Pressure);
          IP.Wo.v.x = IP.Mach_Number*IP.Wo.sound()*cos(TWO*PI*IP.Flow_Angle/360.00);
          IP.Wo.v.y = IP.Mach_Number*IP.Wo.sound()*sin(TWO*PI*IP.Flow_Angle/360.00);
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Time_Max") == 0) {
       i_command = 49;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Time_Max;
       IP.Input_File.getline(buffer, sizeof(buffer));
       IP.Time_Max = IP.Time_Max/THOUSAND;
       if (IP.Time_Max < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Blocks_Per_Processor") == 0) {
       i_command = 50;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Blocks_Per_Processor;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Blocks_Per_Processor < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Output_Format_Type") == 0) {
       i_command = 51;
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

       //} else if (strcmp(IP.Next_Control_Parameter, "NASA_Rotor37_Data_Directory") == 0) {
       //i_command = 46;
       //Get_Next_Input_Control_Parameter(IP);
       //strcpy(IP.NASA_Rotor37_Data_Directory, IP.Next_Control_Parameter);

       //} else if (strcmp(IP.Next_Control_Parameter, "NASA_Rotor67_Data_Directory") == 0) {
       //i_command = 47;
       //Get_Next_Input_Control_Parameter(IP);
       //strcpy(IP.NASA_Rotor67_Data_Directory, IP.Next_Control_Parameter);

       //} else if (strcmp(IP.Next_Control_Parameter, "Rotor_Flow_Type") == 0) {
       //i_command = 48;
       //IP.Line_Number = IP.Line_Number + 1;
       //IP.Input_File >> IP.Rotor_Flow_Type;
       //IP.Input_File.getline(buffer, sizeof(buffer));
       //if (IP.Rotor_Flow_Type != PEAK_FLOW &&
       //    IP.Rotor_Flow_Type != STALL_FLOW) IP.Rotor_Flow_Type = PEAK_FLOW;
       //if (IP.i_Grid == GRID_NASA_ROTOR_37) {
       //   IP.NASA_Rotor37.init(IP.Rotor_Flow_Type,
       //                        IP.NASA_Rotor37_Data_Directory);
       //} else if (IP.i_Grid == GRID_NASA_ROTOR_67) {
       //   IP.NASA_Rotor67.init(IP.Rotor_Flow_Type,
       //                        IP.NASA_Rotor67_Data_Directory);
       //} /* endif */

       //} else if (strcmp(IP.Next_Control_Parameter, "Rotor_Percent_Span") == 0) {
       //i_command = 49;
       //IP.Line_Number = IP.Line_Number + 1;
       //IP.Input_File >> IP.Rotor_Percent_Span;
       //IP.Input_File.getline(buffer, sizeof(buffer));
       //if (IP.Rotor_Percent_Span < ZERO ||
       //    IP.Rotor_Percent_Span > HUNDRED) {
       //   i_command = INVALID_INPUT_VALUE;
       //} else if (IP.i_Grid == GRID_NASA_ROTOR_37) {
       //   IP.Wo.setgas("AIR");
       //   IP.Wo = IP.NASA_Rotor37.getPstateREL_up(IP.Rotor_Percent_Span);
       //   IP.Pressure = IP.Wo.p;
       //   IP.Temperature = IP.Wo.T();
       //   IP.Mach_Number = IP.NASA_Rotor37.getMachREL_up(IP.Rotor_Percent_Span);
       //   IP.Flow_Angle = atan2(IP.Wo.v.y, IP.Wo.v.x); 
       //   if (IP.Flow_Angle < ZERO) IP.Flow_Angle = TWO*PI + IP.Flow_Angle;
       //   IP.Flow_Angle = 180.00*IP.Flow_Angle/PI;
       //   IP.W1 = IP.NASA_Rotor37.getPstateREL_down(IP.Rotor_Percent_Span);
       //} else if (IP.i_Grid == GRID_NASA_ROTOR_67) {
       //   IP.Wo.setgas("AIR");
       //   IP.Wo = IP.NASA_Rotor67.getPstateREL_up(IP.Rotor_Percent_Span);
       //   IP.Pressure = IP.Wo.p;
       //   IP.Temperature = IP.Wo.T();
       //   IP.Mach_Number = IP.NASA_Rotor67.getMachREL_up(IP.Rotor_Percent_Span);
       //   IP.Flow_Angle = atan2(IP.Wo.v.y, IP.Wo.v.x); 
       //   if (IP.Flow_Angle < ZERO) IP.Flow_Angle = TWO*PI + IP.Flow_Angle;
       //   IP.Flow_Angle = 180.00*IP.Flow_Angle/PI;
       //   IP.W1 = IP.NASA_Rotor67.getPstateREL_down(IP.Rotor_Percent_Span);
       //} /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Flow_Geometry_Type") == 0) {
       i_command = 52;
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

    } else if (strcmp(IP.Next_Control_Parameter, "Restart_Solution_Save_Frequency") == 0) {
       i_command = 53;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Restart_Solution_Save_Frequency;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Restart_Solution_Save_Frequency < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Output_Progress_Frequency") == 0) {
       i_command = 51;
       IP.Line_Number++;
       IP.Input_File >> IP.Output_Progress_Frequency;
       IP.Input_File.getline(buffer,sizeof(buffer));
       if (IP.Output_Progress_Frequency < 1) i_command = INVALID_INPUT_VALUE;

       //} else if (strcmp(IP.Next_Control_Parameter, "ICEMCFD_Topology_File") == 0) {
       //i_command = 52;
       //Get_Next_Input_Control_Parameter(IP);
       //strcpy(IP.ICEMCFD_FileNames[0], IP.Next_Control_Parameter);

       //} else if (strcmp(IP.Next_Control_Parameter, "ICEMCFD_Family_Boco_File") == 0) {
       //i_command = 53;
       //Get_Next_Input_Control_Parameter(IP);
       //strcpy(IP.ICEMCFD_FileNames[1], IP.Next_Control_Parameter);

       //} else if (strcmp(IP.Next_Control_Parameter, "ICEMCFD_Family_Topo_File") == 0) {
       //i_command = 54;
       //Get_Next_Input_Control_Parameter(IP);
       //strcpy(IP.ICEMCFD_FileNames[2], IP.Next_Control_Parameter);

    } else if (strcmp(IP.Next_Control_Parameter, "X_Shift") == 0) {
       i_command = 54;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.X_Shift;
       IP.Input_File.setf(ios::skipws);
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "X_Scale") == 0) {
       i_command = 55;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.X_Scale;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "X_Rotate") == 0) {
       i_command = 56;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.X_Rotate;
       IP.Input_File.getline(buffer, sizeof(buffer));

       //} else if (strcmp(IP.Next_Control_Parameter, "GMRES_Restart") == 0) {
      // GMRES restart:
       //i_command = 58;
       //IP.Line_Number = IP.Line_Number + 1;
       //IP.Input_File >> IP.GMRES_Restart;
       //IP.Input_File.getline(buffer, sizeof(buffer));
       //if (IP.GMRES_Restart < 0) i_command = INVALID_INPUT_VALUE;
      
       //} else if (strcmp(IP.Next_Control_Parameter, "GMRES_Overlap") == 0) {
      // GMRES overlap:
       //i_command = 59;
       //IP.Line_Number = IP.Line_Number + 1;
       //IP.Input_File >> IP.GMRES_Overlap;
       //IP.Input_File.getline(buffer, sizeof(buffer));
       //if (IP.GMRES_Overlap < 0) i_command = INVALID_INPUT_VALUE;
      
       //} else if (strcmp(IP.Next_Control_Parameter, "GMRES_Tolerance") == 0) {
      // GMRES tolerance:
       //i_command = 60;
       //IP.Line_Number = IP.Line_Number + 1;
       //IP.Input_File >> IP.GMRES_Toler;
       //IP.Input_File.getline(buffer, sizeof(buffer));
       //if (IP.GMRES_Toler <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Overall_Tolerance") == 0) {
      // GMRES tolerance:
      i_command = 57;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Overall_Toler;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Overall_Toler <= ZERO) i_command = INVALID_INPUT_VALUE;
      
      //} else if (strcmp(IP.Next_Control_Parameter, "GMRES_P_Switch") == 0) {
      // GMRES P_Switch:
      //i_command = 62;
      //IP.Line_Number = IP.Line_Number + 1;
      //IP.Input_File >> IP.GMRES_P_Switch;
      //IP.Input_File.getline(buffer, sizeof(buffer));
      //if (IP.GMRES_P_Switch < 0) i_command = INVALID_INPUT_VALUE;

      //} else if (strcmp(IP.Next_Control_Parameter, "GMRES_ILUK_Level_of_Fill") == 0) {
      // GMRES ILUK_Level_of_Fill:
      //i_command = 63;
      //IP.Line_Number = IP.Line_Number + 1;
      //IP.Input_File >> IP.GMRES_ILUK_Level_of_Fill;
      //IP.Input_File.getline(buffer, sizeof(buffer));
      //if (IP.GMRES_ILUK_Level_of_Fill < 0) i_command = INVALID_INPUT_VALUE;

      //} else if (strcmp(IP.Next_Control_Parameter, "Maximum_Number_of_NKS_Iterations") == 0) {
      //i_command = 64;
      //IP.Line_Number = IP.Line_Number + 1;
      //IP.Input_File >> IP.Maximum_Number_of_NKS_Iterations;
      //IP.Input_File.getline(buffer, sizeof(buffer));
      //if (IP.Maximum_Number_of_NKS_Iterations < 0) i_command = INVALID_INPUT_VALUE;
     
      //} else if (strcmp(IP.Next_Control_Parameter, "Maximum_Number_of_GMRES_Iterations") == 0) {
      //i_command = 65;
      //IP.Line_Number = IP.Line_Number + 1;
      //IP.Input_File >> IP.Maximum_Number_of_GMRES_Iterations;
      //IP.Input_File.getline(buffer, sizeof(buffer));
      //if (IP.Maximum_Number_of_GMRES_Iterations < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Finite_Time_Step") == 0) {
      i_command = 58;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Finite_Time_Step;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Finite_Time_Step < 0) i_command = INVALID_INPUT_VALUE;
      
    } else if (strcmp(IP.Next_Control_Parameter, "Normalization") == 0) {
      i_command = 59;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Normalization;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Normalization < 0) i_command = INVALID_INPUT_VALUE;
      
    } else if (strcmp(IP.Next_Control_Parameter, "Freeze_Limiter") == 0) {
      // Freeze_Limiter:
      i_command = 60;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Freeze_Limiter;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Freeze_Limiter < 0) i_command = INVALID_INPUT_VALUE;
      
    } else if (strcmp(IP.Next_Control_Parameter, "Freeze_Limiter_Residual_Level") == 0) {
      // Freeze_Limiter_Residual_Level:
      i_command = 61;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Freeze_Limiter_Residual_Level;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Freeze_Limiter_Residual_Level < ZERO) i_command = INVALID_INPUT_VALUE;
      
    } else if (strcmp(IP.Next_Control_Parameter, "Finite_Time_Step_Initial_CFL") == 0) {
      i_command = 62;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Finite_Time_Step_Initial_CFL;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Finite_Time_Step_Initial_CFL < ZERO) i_command = INVALID_INPUT_VALUE;

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

    } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Interface_Mesh_Refinements") == 0) {
      i_command = 76;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Number_of_Interface_Mesh_Refinements;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Number_of_Interface_Mesh_Refinements < 0) IP.Number_of_Interface_Mesh_Refinements = 0;

    } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Bounding_Box_Mesh_Refinements") == 0) {
      i_command = 77;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Number_of_Bounding_Box_Mesh_Refinements;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Number_of_Bounding_Box_Mesh_Refinements < 0) IP.Number_of_Bounding_Box_Mesh_Refinements = 0;

    } else if (strcmp(IP.Next_Control_Parameter,"Interface_Refinement_Condition") == 0) {
      i_command = 78;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Interface_Refinement_Condition = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.Interface_Refinement_Condition = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

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

    } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Refinement_Criteria") == 0) {
      i_command = 80;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Number_of_Refinement_Criteria;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Number_of_Refinement_Criteria < 1 || IP.Number_of_Refinement_Criteria > 6) {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Gradient_Density") == 0) {
      i_command = 81;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Refinement_Criteria_Gradient_Density = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.Refinement_Criteria_Gradient_Density = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Divergence_Velocity") == 0) {
      i_command = 82;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Refinement_Criteria_Divergence_Velocity = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.Refinement_Criteria_Divergence_Velocity = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Curl_Velocity") == 0) {
      i_command = 83;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"ON") == 0) {
	IP.Refinement_Criteria_Curl_Velocity = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"OFF") == 0) {
	IP.Refinement_Criteria_Curl_Velocity = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"AMR_Xmin") == 0) {
      i_command = 87;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.AMR_Xmin;
      IP.Input_File.setf(ios::skipws);
      IP.Input_File.getline(buffer,sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter,"AMR_Xmax") == 0) {
      i_command = 88;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.AMR_Xmax;
      IP.Input_File.setf(ios::skipws);
      IP.Input_File.getline(buffer,sizeof(buffer));
      
    } else if (strcmp(IP.Next_Control_Parameter, "Residual_Smoothing_Epsilon") == 0) {
      i_command = 85;
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
      i_command = 86;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Residual_Smoothing_Gauss_Seidel_Iterations;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Residual_Smoothing_Gauss_Seidel_Iterations < 0) {
         IP.Residual_Smoothing_Gauss_Seidel_Iterations = 0;
      } /* endif */

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
      i_command = 91;
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
      } else if (strcmp(IP.BC_North_Type,"Adiabatic_Wall") == 0) {
	IP.BC_North = BC_ADIABATIC_WALL;
      } else if (strcmp(IP.BC_North_Type,"Fixed") == 0) {
	IP.BC_North = BC_FIXED;
      } else if (strcmp(IP.BC_North_Type,"Constant_Extrapolation") == 0) {
	IP.BC_North = BC_CONSTANT_EXTRAPOLATION;
      } else if (strcmp(IP.BC_North_Type,"Linear_Extrapolation") == 0) {
	IP.BC_North = BC_LINEAR_EXTRAPOLATION;
      } else if (strcmp(IP.BC_North_Type,"Characteristic") == 0) {
	IP.BC_North = BC_CHARACTERISTIC;
      } else if (strcmp(IP.BC_North_Type,"Characteristic_Velocity") == 0) {
	IP.BC_North = BC_CHARACTERISTIC_VELOCITY;
      } else if (strcmp(IP.BC_North_Type,"Periodic") == 0) {
	IP.BC_North = BC_PERIODIC;
      } else if (strcmp(IP.BC_North_Type,"Couette") == 0) {
	IP.BC_North = BC_COUETTE;
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
      } else if (strcmp(IP.BC_South_Type,"Adiabatic_Wall") == 0) {
	IP.BC_South = BC_ADIABATIC_WALL;
      } else if (strcmp(IP.BC_South_Type,"Fixed") == 0) {
	IP.BC_South = BC_FIXED;
      } else if (strcmp(IP.BC_South_Type,"Constant_Extrapolation") == 0) {
	IP.BC_South = BC_CONSTANT_EXTRAPOLATION;
      } else if (strcmp(IP.BC_South_Type,"Linear_Extrapolation") == 0) {
	IP.BC_South = BC_LINEAR_EXTRAPOLATION;
      } else if (strcmp(IP.BC_South_Type,"Characteristic") == 0) {
	IP.BC_South = BC_CHARACTERISTIC;
      } else if (strcmp(IP.BC_South_Type,"Characteristic_Velocity") == 0) {
	IP.BC_South = BC_CHARACTERISTIC_VELOCITY;
      } else if (strcmp(IP.BC_South_Type,"Periodic") == 0) {
	IP.BC_South = BC_PERIODIC;
      } else if (strcmp(IP.BC_South_Type,"Couette") == 0) {
	IP.BC_South = BC_COUETTE;
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
      } else if (strcmp(IP.BC_East_Type,"Adiabatic_Wall") == 0) {
	IP.BC_East = BC_ADIABATIC_WALL;
      } else if (strcmp(IP.BC_East_Type,"Fixed") == 0) {
	IP.BC_East = BC_FIXED;
      } else if (strcmp(IP.BC_East_Type,"Constant_Extrapolation") == 0) {
	IP.BC_East = BC_CONSTANT_EXTRAPOLATION;
      } else if (strcmp(IP.BC_East_Type,"Linear_Extrapolation") == 0) {
	IP.BC_East = BC_LINEAR_EXTRAPOLATION;
      } else if (strcmp(IP.BC_East_Type,"Characteristic") == 0) {
	IP.BC_East = BC_CHARACTERISTIC;
      } else if (strcmp(IP.BC_East_Type,"Characteristic_Velocity") == 0) {
	IP.BC_East = BC_CHARACTERISTIC_VELOCITY;
      } else if (strcmp(IP.BC_East_Type,"Periodic") == 0) {
	IP.BC_East = BC_PERIODIC;
      } else if (strcmp(IP.BC_East_Type,"Developed_Channel") == 0) {
	IP.BC_East = BC_DEVELOPED_CHANNEL;
      } else if (strcmp(IP.BC_East_Type,"Couette") == 0) {
	IP.BC_East = BC_COUETTE;
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
      } else if (strcmp(IP.BC_West_Type,"Adiabatic_Wall") == 0) {
	IP.BC_West = BC_ADIABATIC_WALL;
      } else if (strcmp(IP.BC_West_Type,"Fixed") == 0) {
	IP.BC_West = BC_FIXED;
      } else if (strcmp(IP.BC_West_Type,"Constant_Extrapolation") == 0) {
	IP.BC_West = BC_CONSTANT_EXTRAPOLATION;
      } else if (strcmp(IP.BC_West_Type,"Linear_Extrapolation") == 0) {
	IP.BC_West = BC_LINEAR_EXTRAPOLATION;
      } else if (strcmp(IP.BC_West_Type,"Characteristic") == 0) {
	IP.BC_West = BC_CHARACTERISTIC;
      } else if (strcmp(IP.BC_West_Type,"Characteristic_Velocity") == 0) {
	IP.BC_West = BC_CHARACTERISTIC_VELOCITY;
      } else if (strcmp(IP.BC_West_Type,"Periodic") == 0) {
	IP.BC_West = BC_PERIODIC;
      } else if (strcmp(IP.BC_West_Type,"Developed_Channel") == 0) {
	IP.BC_West = BC_DEVELOPED_CHANNEL;
      } else if (strcmp(IP.BC_West_Type,"Couette") == 0) {
	IP.BC_West = BC_COUETTE;
      } else if (strcmp(IP.BC_West_Type,"None") == 0) {
	IP.BC_West = BC_NONE;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

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
      } else if (strcmp(IP.Interface_IP.Type,"Line") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_LINE;
      } else if (strcmp(IP.Interface_IP.Type,"Square") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_SQUARE;
      } else if (strcmp(IP.Interface_IP.Type,"Rectangle") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_RECTANGLE;
      } else if (strcmp(IP.Interface_IP.Type,"Rocket_Propellant_Grain") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_ROCKET_PROPELLANT_GRAIN;
      } else if (strcmp(IP.Interface_IP.Type,"NACA_Aerofoil") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_NACA0012_AEROFOIL;
      } else if (strcmp(IP.Interface_IP.Type,"NASA_Rotor_37") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_NASA_ROTOR_37;
      } else if (strcmp(IP.Interface_IP.Type,"NASA_Rotor_67") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_NASA_ROTOR_67;
      } else if (strcmp(IP.Interface_IP.Type,"Zalesaks_Disk") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_ZALESAK;
      } else if (strcmp(IP.Interface_IP.Type,"Ringleb_Flow") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_RINGLEB;
      } else if (strcmp(IP.Interface_IP.Type,"Bump_Channel_Flow") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_BUMP;
      } else if (strcmp(IP.Interface_IP.Type,"Flat_Plate") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_FLAT_PLATE;
      } else if (strcmp(IP.Interface_IP.Type,"User_Specified") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Type = INTERFACE_USER_SPECIFIED;
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

    } else if (strcmp(IP.Next_Control_Parameter,"Interface_BC_Type") == 0) {
      i_command = 604;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Interface_IP.BC_Type,IP.Next_Control_Parameter);
      if (strcmp(IP.Interface_IP.BC_Type,"Reflection") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].BC_Type = INTERFACE_BC_REFLECTION;
      } else if (strcmp(IP.Interface_IP.BC_Type,"Adiabatic_Wall") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].BC_Type = INTERFACE_BC_ADIABATIC_WALL;
      } else if (strcmp(IP.Interface_IP.BC_Type,"Flat_Plate") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].BC_Type = INTERFACE_BC_FLAT_PLATE;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

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
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_STATIONARY;
      } else if (strcmp(IP.Interface_IP.Motion,"Constant") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_CONSTANT;
      } else if (strcmp(IP.Interface_IP.Motion,"Expand") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_EXPAND;
      } else if (strcmp(IP.Interface_IP.Motion,"Uniform") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_UNIFORM;
      } else if (strcmp(IP.Interface_IP.Motion,"Translate") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_TRANSLATE;
      } else if (strcmp(IP.Interface_IP.Motion,"Stretch") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_STRETCH;
      } else if (strcmp(IP.Interface_IP.Motion,"Rotate") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_ROTATE;
      } else if (strcmp(IP.Interface_IP.Motion,"BurningSurface") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_BURNING_SURFACE;
      } else if (strcmp(IP.Interface_IP.Motion,"MomentumTransfer") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_MOMENTUM_TRANSFER;
      } else if (strcmp(IP.Interface_IP.Motion,"LevelSet_Stationary") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_LEVELSET_STATIONARY;
      } else if (strcmp(IP.Interface_IP.Motion,"LevelSet_Expand") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_LEVELSET_EXPAND;
      } else if (strcmp(IP.Interface_IP.Motion,"LevelSet_Stretch") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_LEVELSET_STRETCH;
      } else if (strcmp(IP.Interface_IP.Motion,"LevelSet_BulkFlow") == 0) {
	IP.Interface_IP.Component_List[IP.Interface_IP.ci].Motion = INTERFACE_MOTION_LEVELSET_BULKFLOW;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }
      IP.Interface_IP.Component_List[IP.Interface_IP.ci].Speed = Vector2D_ZERO;

    } else if (strcmp(IP.Next_Control_Parameter,"Interface_Characteristic_Velocity") == 0) {
      i_command = 607;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Interface_IP.Component_List[IP.Interface_IP.ci].Speed;
      //IP.Input_File >> IP.Interface_IP.Component_List[IP.Interface_IP.ci].Speed.x
      //		    >> IP.Interface_IP.Component_List[IP.Interface_IP.ci].Speed.y;
      IP.Input_File.setf(ios::skipws);
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

    } else if (strcmp(IP.Next_Control_Parameter,"Solution_Redistribution_Type") == 0) {
      i_command = 610;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Interface_IP.Solution_Redistribution_Type,IP.Next_Control_Parameter);
      if (strcmp(IP.Interface_IP.Solution_Redistribution_Type,"None") == 0) {
	IP.Interface_IP.Redistribution_Type = SOLUTION_REDISTRIBUTION_NONE;
      } else if (strcmp(IP.Interface_IP.Solution_Redistribution_Type,"Injection") == 0) {
	IP.Interface_IP.Redistribution_Type = SOLUTION_REDISTRIBUTION_INJECTION;
      } else if (strcmp(IP.Interface_IP.Solution_Redistribution_Type,"Weighted_Injection") == 0) {
	IP.Interface_IP.Redistribution_Type = SOLUTION_REDISTIBUTION_WEIGHTED_INJECTION;
      } else if (strcmp(IP.Interface_IP.Solution_Redistribution_Type,"Neighbour_Average") == 0) {
	IP.Interface_IP.Redistribution_Type = SOLUTION_REDISTIBUTION_NEIGHBOUR_AVERAGE;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter,"Reset_Interface_Motion_Type") == 0) {
      i_command = 604;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Interface_IP.BC_Type,IP.Next_Control_Parameter);
      if (strcmp(IP.Interface_IP.BC_Type,"ON") == 0) {
	IP.Reset_Interface_Motion_Type = ON;
      } else if (strcmp(IP.Interface_IP.BC_Type,"OFF") == 0) {
	IP.Reset_Interface_Motion_Type = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    /**************************************
     * Multigrid Related Input Parameters */
      //} else if (strcmp(IP.Next_Control_Parameter, "Multigrid_Levels") == 0) {
      //i_command = 201;
      //IP.Line_Number = IP.Line_Number + 1;
      //IP.Input_File >> IP.Multigrid_IP.Levels;
      //IP.Input_File.getline(buffer, sizeof(buffer));
      //if (IP.Multigrid_IP.Levels <= ONE) i_command = INVALID_INPUT_VALUE;
       
      // } else if (strcmp(IP.Next_Control_Parameter, "Multigrid_Cycle_Type") == 0) {
      //i_command = 202;
      //Get_Next_Input_Control_Parameter(IP);
      //strcpy(IP.Multigrid_IP.Cycle_Type,
      //	      IP.Next_Control_Parameter);
      //if (strcmp(IP.Multigrid_IP.Cycle_Type, "V") == 0 ||
      //   strcmp(IP.Multigrid_IP.Cycle_Type, "v") == 0) {
      // IP.Multigrid_IP.i_Cycle = MULTIGRID_V_CYCLE;
      //} else if (strcmp(IP.Multigrid_IP.Cycle_Type, "W") == 0 ||
      //	  strcmp(IP.Multigrid_IP.Cycle_Type, "w") == 0) {
      // IP.Multigrid_IP.i_Cycle = MULTIGRID_W_CYCLE;
      //} else {
      // i_command = INVALID_INPUT_VALUE;
      //}
      //} else if (strcmp(IP.Next_Control_Parameter, "Full_Multigrid") == 0) {
      //i_command = 203;
      //IP.Line_Number = IP.Line_Number + 1;
      //IP.Input_File >> IP.Multigrid_IP.Number_of_Cycles_per_Stage_for_Full_Multigrid;
      //IP.Input_File.getline(buffer, sizeof(buffer));
      //if (IP.Multigrid_IP.Number_of_Cycles_per_Stage_for_Full_Multigrid < 0) 
      //  IP.Multigrid_IP.Number_of_Cycles_per_Stage_for_Full_Multigrid = 0;

      //} else if (strcmp(IP.Next_Control_Parameter, "Multigrid_Number_of_Smooths_on_Finest_Level") == 0) {
      //i_command = 204;
      //IP.Line_Number = IP.Line_Number + 1;
      //IP.Input_File >> IP.Multigrid_IP.Number_of_Smooths_on_Finest_Level;
      //IP.Input_File.getline(buffer, sizeof(buffer));
      //if (IP.Multigrid_IP.Number_of_Smooths_on_Finest_Level < ZERO) i_command = INVALID_INPUT_VALUE;

      //} else if (strcmp(IP.Next_Control_Parameter, "Multigrid_Number_of_Pre_Smooths") == 0) {
      //i_command = 205;
      //IP.Line_Number = IP.Line_Number + 1;
      //IP.Input_File >> IP.Multigrid_IP.Number_of_Pre_Smooths;
      //IP.Input_File.getline(buffer, sizeof(buffer));
      //if (IP.Multigrid_IP.Number_of_Pre_Smooths < ZERO) i_command = INVALID_INPUT_VALUE;
       
      //} else if (strcmp(IP.Next_Control_Parameter, "Multigrid_Number_of_Post_Smooths") == 0) {
      //i_command = 206;
      //IP.Line_Number = IP.Line_Number + 1;
      //IP.Input_File >> IP.Multigrid_IP.Number_of_Post_Smooths;
      //IP.Input_File.getline(buffer, sizeof(buffer));
      //if (IP.Multigrid_IP.Number_of_Post_Smooths < ZERO) i_command = INVALID_INPUT_VALUE;

      //} else if (strcmp(IP.Next_Control_Parameter, "Multigrid_Number_of_Smooths_on_Coarsest_Level") == 0) {
      //i_command = 207;
      //IP.Line_Number = IP.Line_Number + 1;
      //IP.Input_File >> IP.Multigrid_IP.Number_of_Smooths_on_Coarsest_Level;
      //IP.Input_File.getline(buffer, sizeof(buffer));
      //if (IP.Multigrid_IP.Number_of_Smooths_on_Coarsest_Level < ZERO) i_command = INVALID_INPUT_VALUE;

      //} else if (strcmp(IP.Next_Control_Parameter, "Convergence_Residual_Level") == 0) {
      //i_command = 208;
      //IP.Line_Number = IP.Line_Number + 1;
      //IP.Input_File >> IP.Multigrid_IP.Convergence_Residual_Level;
      //IP.Input_File.getline(buffer, sizeof(buffer));
      //if (IP.Multigrid_IP.Convergence_Residual_Level < ZERO) i_command = INVALID_INPUT_VALUE;

      //} else if (strcmp(IP.Next_Control_Parameter, "Multigrid_Smoothing_Type") == 0) {
      //i_command = 209;
      //Get_Next_Input_Control_Parameter(IP);
      //strcpy(IP.Multigrid_IP.Smoothing_Type, 
      //       IP.Next_Control_Parameter);
      //if (strcmp(IP.Multigrid_IP.Smoothing_Type, "Explicit_Euler") == 0) {
      //    IP.Multigrid_IP.i_Smoothing = TIME_STEPPING_EXPLICIT_EULER;
      //     IP.N_Stage = 1;
      //} else if (strcmp(IP.Multigrid_IP.Smoothing_Type, "Explicit_Predictor_Corrector") == 0) {
      //    IP.Multigrid_IP.i_Smoothing = TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR;
      //    IP.N_Stage = 2;
      //} else if (strcmp(IP.Multigrid_IP.Smoothing_Type, "Explicit_Runge_Kutta") == 0) {
      //    IP.Multigrid_IP.i_Smoothing = TIME_STEPPING_EXPLICIT_RUNGE_KUTTA;
      //    IP.N_Stage = 5;
      //} else if (strcmp(IP.Multigrid_IP.Smoothing_Type, "Multistage_Optimal_Smoothing") == 0) {
      //    IP.Multigrid_IP.i_Smoothing = TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING;
      //    IP.N_Stage = 4;
      //} else {
      //    IP.Multigrid_IP.i_Smoothing = TIME_STEPPING_EXPLICIT_EULER;
      //    IP.N_Stage = 1;
      //} /* endif */

       /* End of Multigrid related Input parsing *
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
    
    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Elements") == 0) {
      i_command = WRITE_OUTPUT_ELEMENTS_CODE;

    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Nodes") == 0) {
      i_command = WRITE_OUTPUT_NODES_CODE;

    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Cell_Status") == 0) {
      i_command = WRITE_OUTPUT_CELL_STATUS_CODE;

    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Interface_Component_List") == 0) {
      i_command = WRITE_OUTPUT_INTERFACE_COMPONENT_LIST_CODE;

    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Interface_Union_List") == 0) {
      i_command = WRITE_OUTPUT_INTERFACE_UNION_LIST_CODE;

    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Level_Set") == 0) {
      i_command = WRITE_OUTPUT_LEVEL_SET_CODE;
    
    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Level_Set_Cells") == 0) {
      i_command = WRITE_OUTPUT_LEVEL_SET_CELLS_CODE;
    
    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Level_Set_Interface_List") == 0) {
      i_command = WRITE_OUTPUT_LEVEL_SET_INTERFACE_LIST_CODE;

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

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Flat_Plate") == 0) {
       i_command = WRITE_OUTPUT_FLAT_PLATE_CODE;

    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Cell_Status") == 0) {
      i_command = WRITE_OUTPUT_CELL_STATUS_CODE;

    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Interface_Component_List") == 0) {
      i_command = WRITE_OUTPUT_INTERFACE_COMPONENT_LIST_CODE;

    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Interface_Union_List") == 0) {
      i_command = WRITE_OUTPUT_INTERFACE_UNION_LIST_CODE;

    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Level_Set") == 0) {
      i_command = WRITE_OUTPUT_LEVEL_SET_CODE;
    
    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Level_Set_Cells") == 0) {
      i_command = WRITE_OUTPUT_LEVEL_SET_CELLS_CODE;
    
    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Level_Set_Interface_List") == 0) {
      i_command = WRITE_OUTPUT_LEVEL_SET_INTERFACE_LIST_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Cylinder_Free_Molecular") == 0) {
       i_command = WRITE_OUTPUT_CYLINDER_FREE_MOLECULAR_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Drag") == 0) {
       i_command = WRITE_OUTPUT_DRAG_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Cylinder_Drag") == 0) {
       i_command = WRITE_OUTPUT_CYLINDER_DRAG_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Couette") == 0) {
       i_command = WRITE_OUTPUT_COUETTE;

    } else if (strcmp(IP.Next_Control_Parameter,"Write_Output_Aerodynamic_Coefficients") == 0) {
       i_command = WRITE_OUTPUT_AERODYNAMIC_COEFFICIENTS_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Bounding_Box_Refine_Grid") == 0) {
       IP.Number_of_Bounding_Box_Mesh_Refinements = 1;
       i_command = BOUNDING_BOX_REFINE_GRID_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Refine_Grid") == 0) {
       i_command = REFINE_GRID_CODE;

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
int Process_Input_Control_Parameter_File(Gaussian2D_Input_Parameters &Input_Parameters,
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
       cout << "\n Gaussian2D ERROR: Unable to open Gaussian2D input data file.\n";
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
          cout << "\n Gaussian2D ERROR: Error reading Gaussian2D data at line #"
               << -line_number  << " of input data file.\n";
          error_flag = line_number;
          break;
       } /* endif */
    } /* endwhile */

    /* Perform consistency checks on Input_Parameters */

    //if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {
    //  error_flag = Check_Input_Parameters<Euler2D_Input_Parameters>(Input_Parameters);
    //  if (error_flag) {
    //cout << "\n AdvectDiffuse2D ERROR: Input Parameters consistency check failure\n";
    //	return (error_flag);
    //  }
    //}

    /* Initial processing of input control parameters complete.  
       Return the error indicator flag. */

    return (error_flag);

}
