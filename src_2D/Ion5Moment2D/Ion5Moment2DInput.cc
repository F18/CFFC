/* Ion5Moment2DInput.cc:  Subroutines for 
                          2D 5-Moment Ion Transport Model Input Classes. */

/* Include 2D 5-moment ion transport model input parameter header file. */

#ifndef _ION5MOMENT2D_INPUT_INCLUDED
#include "Ion5Moment2DInput.h"
#endif // _ION5MOMENT2D_INPUT_INCLUDED

/*************************************************************
 * Ion5Moment2D_Input_Parameters -- External subroutines.    *
 *************************************************************/

/********************************************************
 * Routine: Open_Input_File                             *
 *                                                      *
 * Opens the appropriate input data file.               *
 *                                                      *
 ********************************************************/
void Open_Input_File(Ion5Moment2D_Input_Parameters &IP) {

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
void Close_Input_File(Ion5Moment2D_Input_Parameters &IP) {

    IP.Input_File.unsetf(ios::skipws);
    IP.Input_File.close();

}

/********************************************************
 * Routine: Set_Default_Input_Parameters                *
 *                                                      *
 * Assigns default values to the input parameters.      *
 *                                                      *
 ********************************************************/
void Set_Default_Input_Parameters(Ion5Moment2D_Input_Parameters &IP) {

    int i;
    char *string_ptr;

    string_ptr = "Ion5Moment2D.in";
    strcpy(IP.Input_File_Name, string_ptr);

    string_ptr = "Explicit_Euler";
    strcpy(IP.Time_Integration_Type, string_ptr);
    IP.i_Time_Integration = TIME_STEPPING_EXPLICIT_EULER;
    IP.Time_Accurate = 0;
    IP.Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
    IP.Maximum_Number_of_Time_Steps = 100;
    IP.N_Stage = 1;
    IP.Time_Max = ZERO;

    // Residual variable:
    IP.i_Residual_Variable = 1;

    IP.Residual_Smoothing = 0;
    IP.Residual_Smoothing_Epsilon = ZERO;
    IP.Residual_Smoothing_Gauss_Seidel_Iterations = 2;

    string_ptr = "Least_Squares";
    strcpy(IP.Reconstruction_Type, string_ptr);
    IP.i_Reconstruction = RECONSTRUCTION_LEAST_SQUARES;

    string_ptr = "Barth_Jespersen";
    strcpy(IP.Limiter_Type, string_ptr);
    IP.i_Limiter = LIMITER_BARTH_JESPERSEN;

    string_ptr = "HLLE";
    strcpy(IP.Flux_Function_Type, string_ptr);
    IP.i_Flux_Function = FLUX_FUNCTION_HLLE;

    string_ptr = "Uniform";
    strcpy(IP.ICs_Type, string_ptr);
    IP.i_ICs = IC_UNIFORM;

    // Ion solution states.
    string_ptr = "H+";  
    strcpy(IP.Ion_Type, string_ptr);
    IP.Wo.setion(IP.Ion_Type);
    IP.Wo = Ion5Moment2D_W_STDATM;
    IP.Ion_Pressure = IP.Wo.p;
    IP.Ion_Temperature = IP.Wo.T();
    IP.Ion_Mach_Number = 0.80;
    IP.Ion_Flow_Angle = ZERO;
    IP.Wo.v.x = IP.Ion_Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Ion_Flow_Angle/360.00);
    IP.Wo.v.y = IP.Ion_Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Ion_Flow_Angle/360.00);
    IP.Uo.setion(IP.Ion_Type);
    IP.Uo = U(IP.Wo);
    IP.W1 = IP.Wo;
    IP.W2 = IP.Wo;

    // Neutral gas solution states.
    string_ptr = "AIR";
    strcpy(IP.Neutral_Gas_Type, string_ptr);
    IP.Wno.setgas(IP.Neutral_Gas_Type);
    IP.Wno = Euler2D_W_STDATM;
    IP.Neutral_Gas_Pressure = IP.Wno.p;
    IP.Neutral_Gas_Temperature = IP.Wno.T();
    IP.Neutral_Gas_Mach_Number = 0.80;
    IP.Neutral_Gas_Flow_Angle = ZERO;
    IP.Wno.v.x = IP.Neutral_Gas_Mach_Number*IP.Wno.a()*cos(TWO*PI*IP.Neutral_Gas_Flow_Angle/360.00);
    IP.Wno.v.y = IP.Neutral_Gas_Mach_Number*IP.Wno.a()*sin(TWO*PI*IP.Neutral_Gas_Flow_Angle/360.00);

    // Set neutral gas type for ions.
    IP.Wo.setgas(IP.Neutral_Gas_Type);
    IP.Uo.setgas(IP.Neutral_Gas_Type);

    string_ptr = "none";
    strcpy(IP.Neutral_Gas_Solution_File_Name, string_ptr);

    // Set electric field parameters.
    string_ptr = "none";
    strcpy(IP.Electric_Field_Solution_File_Name, string_ptr);
    IP.Electric_Field_Strength = ZERO;
    IP.Electric_Field_Angle = ZERO;
    string_ptr = "Uniform";
    strcpy(IP.Electric_Field_Type, string_ptr);
    IP.i_Electric_Field = IC_ELECTRIC_FIELD_UNIFORM;
    IP.Add_Initial_and_Solution_File_Electric_Fields = 0;

    // Set geometry and mesh parameters.
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
    IP.Ellipse_Length_X_Axis = TWO;
    IP.Ellipse_Length_Y_Axis = HALF;
    IP.Chord_Length = ONE;
    IP.Orifice_Radius = ONE;

    IP.X_Shift = Vector2D_ZERO;
    IP.X_Scale = ONE;
    IP.X_Rotate = ZERO;

    // Mesh stretching factor:
    IP.i_Mesh_Stretching = OFF;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_LINEAR;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_LINEAR;
    IP.Mesh_Stretching_Factor_Idir = 1.01;
    IP.Mesh_Stretching_Factor_Jdir = 1.01;

    // Boundary conditions:
    string_ptr = "Off";
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
    IP.Refinement_Criteria_Gradient_Density = ON;
    IP.Refinement_Criteria_Divergence_Velocity = ON;
    IP.Refinement_Criteria_Curl_Velocity = ON;

    // Smooth quad block indicator:
    IP.i_Smooth_Quad_Block = ON;

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

}

/********************************************************
 * Routine: Broadcast_Input_Parameters                  *
 *                                                      *
 * Broadcast the input parameters variables to all      *
 * processors involved in the calculation from the      *
 * primary processor using the MPI broadcast routine.   *
 *                                                      *
 ********************************************************/
void Broadcast_Input_Parameters(Ion5Moment2D_Input_Parameters &IP) {

#ifdef _MPI_VERSION
    int i;

    MPI::COMM_WORLD.Bcast(IP.Input_File_Name, 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Line_Number), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Time_Integration_Type, 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
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
    // Residual variable:
    MPI::COMM_WORLD.Bcast(&(IP.i_Residual_Variable),
			  1,
			  MPI::DOUBLE,0);
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
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Reconstruction), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Limiter_Type, 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Limiter), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Flux_Function_Type, 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Flux_Function), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.ICs_Type, 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_ICs), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Ion_Type, 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Neutral_Gas_Type, 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Ion_Pressure), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Ion_Temperature), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Ion_Mach_Number), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Ion_Flow_Angle), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Neutral_Gas_Pressure), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Neutral_Gas_Temperature), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Neutral_Gas_Mach_Number), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Neutral_Gas_Flow_Angle), 
                          1, 
                          MPI::DOUBLE, 0);
    if (!CFDkit_Primary_MPI_Processor()) {
       IP.Wo.setion(IP.Ion_Type);
       IP.Wo = Ion5Moment2D_pState(IP.Ion_Pressure/(IP.Wo.R*IP.Ion_Temperature), 
                                   ZERO, 
                                   ZERO, 
                                   IP.Ion_Pressure);
       IP.Wo.v.x = IP.Ion_Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Ion_Flow_Angle/360.00);
       IP.Wo.v.y = IP.Ion_Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Ion_Flow_Angle/360.00);
       IP.Uo.setion(IP.Ion_Type);
       IP.Uo = U(IP.Wo);
       IP.Wno.setgas(IP.Neutral_Gas_Type); // for neutral gas solution states
       IP.Wo.setgas(IP.Neutral_Gas_Type); // for ion solution states
       IP.Uo.setgas(IP.Neutral_Gas_Type);
       IP.Wno = Euler2D_pState(IP.Neutral_Gas_Pressure/(IP.Wno.R*IP.Neutral_Gas_Temperature), 
                               ZERO, 
                               ZERO, 
                               IP.Neutral_Gas_Pressure);
       IP.Wno.v.x = IP.Neutral_Gas_Mach_Number*IP.Wno.a()*cos(TWO*PI*IP.Neutral_Gas_Flow_Angle/360.00);
       IP.Wno.v.y = IP.Neutral_Gas_Mach_Number*IP.Wno.a()*sin(TWO*PI*IP.Neutral_Gas_Flow_Angle/360.00);
    } /* endif */
    MPI::COMM_WORLD.Bcast(&(IP.W1.d), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W1.v.x), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W1.v.y), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.W1.p), 
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
    MPI::COMM_WORLD.Bcast(&(IP.W2.p), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(IP.Neutral_Gas_Solution_File_Name, 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                          MPI::CHAR, 0);

    // Electric field solution parameters.
    MPI::COMM_WORLD.Bcast(IP.Electric_Field_Solution_File_Name, 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Electric_Field_Strength), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Electric_Field_Angle), 
                          1, 
                          MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(IP.Electric_Field_Type, 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.i_Electric_Field), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Add_Initial_and_Solution_File_Electric_Fields), 
                          1, 
                          MPI::INT, 0);

    // Grids
    MPI::COMM_WORLD.Bcast(IP.Flow_Geometry_Type, 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(&(IP.Axisymmetric), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(IP.Grid_Type, 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.NACA_Aerofoil_Type, 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
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
			  INPUT_PARAMETER_LENGTH_ION5MOMENT2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(&(IP.BCs_Specified),
			  1,
			MPI::INT,0);
    MPI::COMM_WORLD.Bcast(IP.BC_North_Type,
			  INPUT_PARAMETER_LENGTH_ION5MOMENT2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(IP.BC_South_Type,
			  INPUT_PARAMETER_LENGTH_ION5MOMENT2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(IP.BC_East_Type,
			  INPUT_PARAMETER_LENGTH_ION5MOMENT2D,
			  MPI::CHAR,0);
    MPI::COMM_WORLD.Bcast(IP.BC_West_Type,
			  INPUT_PARAMETER_LENGTH_ION5MOMENT2D,
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
    if (!CFDkit_Primary_MPI_Processor()) {
       IP.ICEMCFD_FileNames = new char*[3];
       for (i = 0; i < 3; i++) {
          IP.ICEMCFD_FileNames[i] = new char[INPUT_PARAMETER_LENGTH_ION5MOMENT2D];
       } /* endfor */
    } /* endif */
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[0], 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[1], 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.ICEMCFD_FileNames[2], 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
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
    // File Names
    MPI::COMM_WORLD.Bcast(IP.Output_File_Name, 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Grid_File_Name, 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Grid_Definition_File_Name, 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Restart_File_Name, 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Gnuplot_File_Name, 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                          MPI::CHAR, 0);
    MPI::COMM_WORLD.Bcast(IP.Output_Format_Type, 
                          INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
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
    } /* endif */
    MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Per_Processor), 
                          1, 
                          MPI::INT, 0);
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
void Broadcast_Input_Parameters(Ion5Moment2D_Input_Parameters &IP,
                                MPI::Intracomm &Communicator, 
                                const int Source_CPU) {

    int Source_Rank = 0;
    int i;

    Communicator.Bcast(IP.Input_File_Name, 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.Line_Number), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Time_Integration_Type, 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
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
    // Residual variable:
    Communicator.Bcast(&(IP.i_Residual_Variable),
		       1,
		       MPI::DOUBLE,Source_Rank);
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
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Reconstruction), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Limiter_Type, 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Limiter), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Flux_Function_Type, 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Flux_Function), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.ICs_Type, 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_ICs), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Ion_Type, 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Neutral_Gas_Type, 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.Ion_Pressure), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Ion_Temperature), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Ion_Mach_Number), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Ion_Flow_Angle), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Neutral_Gas_Pressure), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Neutral_Gas_Temperature), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Neutral_Gas_Mach_Number), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Neutral_Gas_Flow_Angle), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    if (!(CFDkit_MPI::This_Processor_Number == Source_CPU)) {
       IP.Wo.setion(IP.Ion_Type);
       IP.Wo = Ion5Moment2D_pState(IP.Ion_Pressure/(IP.Wo.R*IP.Ion_Temperature), 
                                   ZERO, 
                                   ZERO, 
                                   IP.Ion_Pressure);
       IP.Wo.v.x = IP.Ion_Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Ion_Flow_Angle/360.00);
       IP.Wo.v.y = IP.Ion_Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Ion_Flow_Angle/360.00);
       IP.Uo.setion(IP.Ion_Type);
       IP.Uo = U(IP.Wo);
       IP.Wno.setgas(IP.Neutral_Gas_Type); // for neutral gas solution states
       IP.Wo.setgas(IP.Neutral_Gas_Type); // for ion solution states
       IP.Uo.setgas(IP.Neutral_Gas_Type);
       IP.Wno = Euler2D_pState(IP.Neutral_Gas_Pressure/(IP.Wno.R*IP.Neutral_Gas_Temperature), 
                               ZERO, 
                               ZERO, 
                               IP.Neutral_Gas_Pressure);
       IP.Wno.v.x = IP.Neutral_Gas_Mach_Number*IP.Wno.a()*cos(TWO*PI*IP.Neutral_Gas_Flow_Angle/360.00);
       IP.Wno.v.y = IP.Neutral_Gas_Mach_Number*IP.Wno.a()*sin(TWO*PI*IP.Neutral_Gas_Flow_Angle/360.00);
    } /* endif */
    Communicator.Bcast(&(IP.W1.d), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W1.v.x), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W1.v.y), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.W1.p), 
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
    Communicator.Bcast(&(IP.W2.p), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(IP.Neutral_Gas_Solution_File_Name, 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                       MPI::CHAR, Source_Rank);
    // Electric field solution parameters.
    Communicator.Bcast(IP.Electric_Field_Solution_File_Name, 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.Electric_Field_Strength), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(IP.Electric_Field_Angle), 
                       1, 
                       MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(IP.Electric_Field_Type, 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.i_Electric_Field), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(&(IP.Add_Initial_and_Solution_File_Electric_Fields), 
                       1, 
                       MPI::INT, Source_Rank);
    // Grids.
    Communicator.Bcast(IP.Flow_Geometry_Type, 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(&(IP.Axisymmetric), 
                       1, 
                       MPI::INT, Source_Rank);
    Communicator.Bcast(IP.Grid_Type, 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.NACA_Aerofoil_Type, 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
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
		       INPUT_PARAMETER_LENGTH_ION5MOMENT2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(&(IP.BCs_Specified),
		       1,
		       MPI::INT,Source_Rank);
    Communicator.Bcast(IP.BC_North_Type,
		       INPUT_PARAMETER_LENGTH_ION5MOMENT2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(IP.BC_South_Type,
		       INPUT_PARAMETER_LENGTH_ION5MOMENT2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(IP.BC_East_Type,
		       INPUT_PARAMETER_LENGTH_ION5MOMENT2D,
		       MPI::CHAR,Source_Rank);
    Communicator.Bcast(IP.BC_West_Type,
		       INPUT_PARAMETER_LENGTH_ION5MOMENT2D,
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
    if (!(CFDkit_MPI::This_Processor_Number == Source_CPU)) {
       IP.ICEMCFD_FileNames = new char*[3];
       for (i = 0; i < 3; i++) {
          IP.ICEMCFD_FileNames[i] = new char[INPUT_PARAMETER_LENGTH_ION5MOMENT2D];
       } /* endfor */
    } /* endif */
    Communicator.Bcast(IP.ICEMCFD_FileNames[0], 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.ICEMCFD_FileNames[1], 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.ICEMCFD_FileNames[2], 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
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
    // File Names
    Communicator.Bcast(IP.Output_File_Name, 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Grid_File_Name, 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Grid_Definition_File_Name, 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Restart_File_Name, 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Gnuplot_File_Name, 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
                       MPI::CHAR, Source_Rank);
    Communicator.Bcast(IP.Output_Format_Type, 
                       INPUT_PARAMETER_LENGTH_ION5MOMENT2D, 
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
    } /* endif */
    Communicator.Bcast(&(IP.Number_of_Blocks_Per_Processor), 
                       1, 
                       MPI::INT, Source_Rank);

}
#endif

/********************************************************
 * Routine: Get_Next_Input_Control_Parameter            *
 *                                                      *
 * Get the next input control parameter from the input  *
 * file.                                                *
 *                                                      *
 ********************************************************/
void Get_Next_Input_Control_Parameter(Ion5Moment2D_Input_Parameters &IP) {

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
int Parse_Next_Input_Control_Parameter(Ion5Moment2D_Input_Parameters &IP) {

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

    } else if (strcmp(IP.Next_Control_Parameter, "Flux_Function_Type") == 0) {
       i_command = 4;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Flux_Function_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Flux_Function_Type, "Godunov") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV;
       } else if (strcmp(IP.Flux_Function_Type, "Roe") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_ROE;
       } else if (strcmp(IP.Flux_Function_Type, "Rusanov") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_RUSANOV;
       } else if (strcmp(IP.Flux_Function_Type, "HLLE") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_HLLE;
       } else if (strcmp(IP.Flux_Function_Type, "Linde") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_LINDE;
       } else if (strcmp(IP.Flux_Function_Type, "HLLC") == 0) {
          IP.i_Flux_Function = FLUX_FUNCTION_HLLC;
       } else {
          IP.i_Flux_Function = FLUX_FUNCTION_ROE;
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
       if (strcmp(IP.Grid_Type, "Square") == 0) {
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
	  IP.BC_South = BC_REFLECTION;
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
	  IP.BC_North = BC_REFLECTION;
       } else if (strcmp(IP.Grid_Type, "Circular_Cylinder") == 0) {
          IP.i_Grid = GRID_CIRCULAR_CYLINDER;
          IP.Cylinder_Radius = ONE;
	  IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
	  IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
	  IP.Mesh_Stretching_Factor_Idir = 1.025;
	  IP.Mesh_Stretching_Factor_Jdir = 1.001;
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
       i_command = 12;
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
       } else {
          IP.Local_Time_Stepping = SCALAR_LOCAL_TIME_STEPPING;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Local_Time_Stepping") == 0) {
       i_command = 15;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Local_Time_Stepping;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Local_Time_Stepping != GLOBAL_TIME_STEPPING &&
           IP.Local_Time_Stepping != SCALAR_LOCAL_TIME_STEPPING) 
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

    } else if (strcmp(IP.Next_Control_Parameter,"Chamber_Length") == 0) {
      i_command = 32;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Chamber_Length;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Chamber_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Chamber_Radius") == 0) {
      i_command = 33;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Chamber_Radius;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Chamber_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Chamber_To_Throat_Length") == 0) {
      i_command = 34;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Chamber_To_Throat_Length;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Chamber_To_Throat_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

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

    } else if (strcmp(IP.Next_Control_Parameter,"Nozzle_Type") == 0) {
      i_command = 38;
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
      i_command = 39;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Grain_Radius;
      IP.Input_File.getline(buffer,sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Ion_Type") == 0) {
       i_command = 40;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Ion_Type, 
              IP.Next_Control_Parameter);
       IP.Wo.setion(IP.Ion_Type);
       IP.Wo = Ion5Moment2D_pState(IP.Ion_Pressure/(IP.Wo.R*IP.Ion_Temperature), 
                                   ZERO, 
                                   ZERO, 
                                   IP.Ion_Pressure);
       IP.Wo.v.x = IP.Ion_Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Ion_Flow_Angle/360.00);
       IP.Wo.v.y = IP.Ion_Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Ion_Flow_Angle/360.00);
       IP.Uo.setion(IP.Ion_Type);
       IP.Uo = U(IP.Wo);

    } else if (strcmp(IP.Next_Control_Parameter, "Ion_Mach_Number") == 0) {
       i_command = 41;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Ion_Mach_Number;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Ion_Mach_Number < ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } else {
          IP.Wo.v.x = IP.Ion_Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Ion_Flow_Angle/360.00);
          IP.Wo.v.y = IP.Ion_Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Ion_Flow_Angle/360.00);
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Ion_Flow_Angle") == 0) {
       i_command = 42;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Ion_Flow_Angle;
       IP.Input_File.getline(buffer, sizeof(buffer));
       IP.Wo.v.x = IP.Ion_Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Ion_Flow_Angle/360.00);
       IP.Wo.v.y = IP.Ion_Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Ion_Flow_Angle/360.00);

    } else if (strcmp(IP.Next_Control_Parameter, "Ion_Pressure") == 0) {
       i_command = 43;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Ion_Pressure;
       IP.Input_File.getline(buffer, sizeof(buffer));
       IP.Ion_Pressure = IP.Ion_Pressure*THOUSAND;
       if (IP.Ion_Pressure <= ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } else {
          IP.Wo = Ion5Moment2D_pState(IP.Ion_Pressure/(IP.Wo.R*IP.Ion_Temperature), 
                                      ZERO, 
                                      ZERO, 
                                      IP.Ion_Pressure);
          IP.Wo.v.x = IP.Ion_Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Ion_Flow_Angle/360.00);
          IP.Wo.v.y = IP.Ion_Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Ion_Flow_Angle/360.00);
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Ion_Temperature") == 0) {
       i_command = 44;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Ion_Temperature;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Ion_Temperature <= ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } else {
          IP.Wo = Ion5Moment2D_pState(IP.Ion_Pressure/(IP.Wo.R*IP.Ion_Temperature), 
                                      ZERO, 
                                      ZERO, 
                                      IP.Ion_Pressure);
          IP.Wo.v.x = IP.Ion_Mach_Number*IP.Wo.a()*cos(TWO*PI*IP.Ion_Flow_Angle/360.00);
          IP.Wo.v.y = IP.Ion_Mach_Number*IP.Wo.a()*sin(TWO*PI*IP.Ion_Flow_Angle/360.00);
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Neutral_Gas_Type") == 0) {
       i_command = 45;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Neutral_Gas_Type, 
              IP.Next_Control_Parameter);
       IP.Wno.setgas(IP.Neutral_Gas_Type); // for neutral gas solution states
       IP.Wo.setgas(IP.Neutral_Gas_Type); // for ion solution states
       IP.Uo.setgas(IP.Neutral_Gas_Type);
       IP.Wno = Euler2D_pState(IP.Neutral_Gas_Pressure/(IP.Wno.R*IP.Neutral_Gas_Temperature), 
                               ZERO, 
                               ZERO, 
                               IP.Neutral_Gas_Pressure);
       IP.Wno.v.x = IP.Neutral_Gas_Mach_Number*IP.Wno.a()*cos(TWO*PI*IP.Neutral_Gas_Flow_Angle/360.00);
       IP.Wno.v.y = IP.Neutral_Gas_Mach_Number*IP.Wno.a()*sin(TWO*PI*IP.Neutral_Gas_Flow_Angle/360.00);

    } else if (strcmp(IP.Next_Control_Parameter, "Neutral_Gas_Mach_Number") == 0) {
       i_command = 46;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Neutral_Gas_Mach_Number;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Neutral_Gas_Mach_Number < ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } else {
          IP.Wno.v.x = IP.Neutral_Gas_Mach_Number*IP.Wno.a()*cos(TWO*PI*IP.Neutral_Gas_Flow_Angle/360.00);
          IP.Wno.v.y = IP.Neutral_Gas_Mach_Number*IP.Wno.a()*sin(TWO*PI*IP.Neutral_Gas_Flow_Angle/360.00);
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Neutral_Gas_Flow_Angle") == 0) {
       i_command = 47;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Neutral_Gas_Flow_Angle;
       IP.Input_File.getline(buffer, sizeof(buffer));
       IP.Wno.v.x = IP.Neutral_Gas_Mach_Number*IP.Wno.a()*cos(TWO*PI*IP.Neutral_Gas_Flow_Angle/360.00);
       IP.Wno.v.y = IP.Neutral_Gas_Mach_Number*IP.Wno.a()*sin(TWO*PI*IP.Neutral_Gas_Flow_Angle/360.00);

    } else if (strcmp(IP.Next_Control_Parameter, "Neutral_Gas_Pressure") == 0) {
       i_command = 48;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Neutral_Gas_Pressure;
       IP.Input_File.getline(buffer, sizeof(buffer));
       IP.Neutral_Gas_Pressure = IP.Neutral_Gas_Pressure*THOUSAND;
       if (IP.Neutral_Gas_Pressure <= ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } else {
          IP.Wno = Euler2D_pState(IP.Neutral_Gas_Pressure/(IP.Wno.R*IP.Neutral_Gas_Temperature), 
                                  ZERO, 
                                  ZERO, 
                                  IP.Neutral_Gas_Pressure);
          IP.Wno.v.x = IP.Neutral_Gas_Mach_Number*IP.Wno.a()*cos(TWO*PI*IP.Neutral_Gas_Flow_Angle/360.00);
          IP.Wno.v.y = IP.Neutral_Gas_Mach_Number*IP.Wno.a()*sin(TWO*PI*IP.Neutral_Gas_Flow_Angle/360.00);
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Neutral_Gas_Temperature") == 0) {
       i_command = 49;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Neutral_Gas_Temperature;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Neutral_Gas_Temperature <= ZERO) {
          i_command = INVALID_INPUT_VALUE;
       } else {
          IP.Wno = Euler2D_pState(IP.Neutral_Gas_Pressure/(IP.Wno.R*IP.Neutral_Gas_Temperature), 
                                  ZERO, 
                                  ZERO, 
                                  IP.Neutral_Gas_Pressure);
          IP.Wno.v.x = IP.Neutral_Gas_Mach_Number*IP.Wno.a()*cos(TWO*PI*IP.Neutral_Gas_Flow_Angle/360.00);
          IP.Wno.v.y = IP.Neutral_Gas_Mach_Number*IP.Wno.a()*sin(TWO*PI*IP.Neutral_Gas_Flow_Angle/360.00);
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Time_Max") == 0) {
       i_command = 50;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Time_Max;
       IP.Input_File.getline(buffer, sizeof(buffer));
       IP.Time_Max = IP.Time_Max/THOUSAND;
       if (IP.Time_Max < ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Blocks_Per_Processor") == 0) {
       i_command = 51;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Blocks_Per_Processor;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Blocks_Per_Processor < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Output_Format_Type") == 0) {
       i_command = 52;
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
       i_command = 53;
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
       i_command = 54;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Restart_Solution_Save_Frequency;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Restart_Solution_Save_Frequency < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter,"Output_Progress_Frequency") == 0) {
       i_command = 55;
       IP.Line_Number++;
       IP.Input_File >> IP.Output_Progress_Frequency;
       IP.Input_File.getline(buffer,sizeof(buffer));
       if (IP.Output_Progress_Frequency < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "ICEMCFD_Topology_File") == 0) {
       i_command = 56;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.ICEMCFD_FileNames[0], IP.Next_Control_Parameter);

    } else if (strcmp(IP.Next_Control_Parameter, "ICEMCFD_Family_Boco_File") == 0) {
       i_command = 57;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.ICEMCFD_FileNames[1], IP.Next_Control_Parameter);

    } else if (strcmp(IP.Next_Control_Parameter, "ICEMCFD_Family_Topo_File") == 0) {
       i_command = 58;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.ICEMCFD_FileNames[2], IP.Next_Control_Parameter);

    } else if (strcmp(IP.Next_Control_Parameter, "X_Shift") == 0) {
       i_command = 59;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.X_Shift;
       IP.Input_File.setf(ios::skipws);
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "X_Scale") == 0) {
       i_command = 60;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.X_Scale;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "X_Rotate") == 0) {
       i_command = 61;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.X_Rotate;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Neutral_Gas_Solution_File_Name") == 0) {
       i_command = 62;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Neutral_Gas_Solution_File_Name, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Neutral_Gas_Solution_File_Name, "none") != 0 &&
           strcmp(IP.Neutral_Gas_Solution_File_Name, "NONE") != 0) {
          strcat(IP.Neutral_Gas_Solution_File_Name, ".dat");
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Electric_Field_Solution_File_Name") == 0) {
       i_command = 63;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Electric_Field_Solution_File_Name, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Electric_Field_Solution_File_Name, "none") != 0 &&
           strcmp(IP.Electric_Field_Solution_File_Name, "NONE") != 0) {
          strcat(IP.Electric_Field_Solution_File_Name, ".dat");
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Electric_Field_Strength") == 0) {
       i_command = 64;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Electric_Field_Strength;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Electric_Field_Angle") == 0) {
       i_command = 65;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Electric_Field_Angle;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Electric_Field_Type") == 0) {
       i_command = 66;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Electric_Field_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Electric_Field_Type, "Uniform") == 0) {
          IP.i_Electric_Field = IC_ELECTRIC_FIELD_UNIFORM;
       } else if (strcmp(IP.Electric_Field_Type, "Quadrupole") == 0) {
          IP.i_Electric_Field = IC_ELECTRIC_FIELD_QUADRUPOLE;
       } else if (strcmp(IP.Electric_Field_Type, "Octapole") == 0) {
          IP.i_Electric_Field = IC_ELECTRIC_FIELD_OCTAPOLE;
       } else {
          IP.i_Electric_Field = IC_ELECTRIC_FIELD_UNIFORM;
       } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter, "Add_Initial_and_Solution_File_Electric_Fields") == 0) {
      i_command = 67;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Add_Initial_and_Solution_File_Electric_Fields;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Add_Initial_and_Solution_File_Electric_Fields < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "AMR") == 0) {
      i_command = 68;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	IP.AMR = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	IP.AMR = OFF;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

    } else if (strcmp(IP.Next_Control_Parameter, "AMR_Frequency") == 0) {
      i_command = 69;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.AMR_Frequency;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.AMR_Frequency < 0) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Initial_Mesh_Refinements") == 0) {
      i_command = 70;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Number_of_Initial_Mesh_Refinements;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Number_of_Initial_Mesh_Refinements < 0) IP.Number_of_Initial_Mesh_Refinements = 0;

    } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Uniform_Mesh_Refinements") == 0) {
      i_command = 71;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Number_of_Uniform_Mesh_Refinements;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Number_of_Uniform_Mesh_Refinements < 0) IP.Number_of_Uniform_Mesh_Refinements = 0;

    } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Boundary_Mesh_Refinements") == 0) {
      i_command = 72;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Number_of_Boundary_Mesh_Refinements;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Number_of_Boundary_Mesh_Refinements < 0) IP.Number_of_Boundary_Mesh_Refinements = 0;

    } else if (strcmp(IP.Next_Control_Parameter,"Maximum_Refinement_Level") == 0) {
      i_command = 73;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Maximum_Refinement_Level;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Maximum_Refinement_Level < 1) IP.Maximum_Refinement_Level = 1;

    } else if (strcmp(IP.Next_Control_Parameter,"Minimum_Refinement_Level") == 0) {
      i_command = 74;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Minimum_Refinement_Level;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.Minimum_Refinement_Level < 1) IP.Minimum_Refinement_Level = 1;

    } else if (strcmp(IP.Next_Control_Parameter, "Threshold_for_Refinement") == 0) {
       i_command = 75;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Threshold_for_Refinement;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Threshold_for_Refinement <= ZERO ||
           IP.Threshold_for_Refinement > ONE) IP.Threshold_for_Refinement = 0.50;

    } else if (strcmp(IP.Next_Control_Parameter, "Threshold_for_Coarsening") == 0) {
       i_command = 76;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Threshold_for_Coarsening;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Threshold_for_Coarsening < ZERO ||
           IP.Threshold_for_Coarsening >= ONE) IP.Threshold_for_Coarsening = 0.10;

    } else if (strcmp(IP.Next_Control_Parameter,"Number_of_Refinement_Criteria") == 0) {
       i_command = 77;
       IP.Line_Number = IP.Line_Number + 1;
       IP.Input_File >> IP.Number_of_Refinement_Criteria;
       IP.Input_File.getline(buffer,sizeof(buffer));
       if (IP.Number_of_Refinement_Criteria < 1 || IP.Number_of_Refinement_Criteria > 6) {
	 i_command = INVALID_INPUT_VALUE;
       }

    } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Gradient_Density") == 0) {
       i_command = 78;
       Get_Next_Input_Control_Parameter(IP);
       if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	 IP.Refinement_Criteria_Gradient_Density = ON;
       } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	 IP.Refinement_Criteria_Gradient_Density = OFF;
       } else {
	 i_command = INVALID_INPUT_VALUE;
       }

    } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Divergence_Velocity") == 0) {
       i_command = 79;
       Get_Next_Input_Control_Parameter(IP);
       if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	 IP.Refinement_Criteria_Divergence_Velocity = ON;
       } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	 IP.Refinement_Criteria_Divergence_Velocity = OFF;
       } else {
	i_command = INVALID_INPUT_VALUE;
       }

    } else if (strcmp(IP.Next_Control_Parameter,"Refinement_Criteria_Curl_Velocity") == 0) {
       i_command = 80;
       Get_Next_Input_Control_Parameter(IP);
       if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	 IP.Refinement_Criteria_Curl_Velocity = ON;
       } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
	 IP.Refinement_Criteria_Curl_Velocity = OFF;
       } else {
	 i_command = INVALID_INPUT_VALUE;
       }

    } else if (strcmp(IP.Next_Control_Parameter,"Residual_Variable") == 0) {
      i_command = 81;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.i_Residual_Variable;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.i_Residual_Variable < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Residual_Smoothing_Epsilon") == 0) {
      i_command = 82;
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
      i_command = 83;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.Residual_Smoothing_Gauss_Seidel_Iterations;
      IP.Input_File.getline(buffer, sizeof(buffer));
      if (IP.Residual_Smoothing_Gauss_Seidel_Iterations < 0) {
         IP.Residual_Smoothing_Gauss_Seidel_Iterations = 0;
      } /* endif */

    } else if (strcmp(IP.Next_Control_Parameter,"Mesh_Stretching") == 0) {
      i_command = 101;
      Get_Next_Input_Control_Parameter(IP);
      if (strcmp(IP.Next_Control_Parameter,"On") == 0) {
	IP.i_Mesh_Stretching = ON;
      } else if (strcmp(IP.Next_Control_Parameter,"Off") == 0) {
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
      i_command = 106;
      IP.Line_Number = IP.Line_Number + 1;
      IP.Input_File >> IP.i_Smooth_Quad_Block;
      IP.Input_File.getline(buffer,sizeof(buffer));
      if (IP.i_Smooth_Quad_Block < 0 || IP.i_Smooth_Quad_Block > 1) i_command = INVALID_INPUT_VALUE;

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

    } else if (strcmp(IP.Next_Control_Parameter, "Multigrid_Smoothing_Type") == 0) {
       i_command = 208;
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

       /* End of Multigrid related Input parsing *
        ******************************************/

       ////////////////////////////////////////////////////////////////////
       // SPECIFIED BOUNDARY CONDITIONS                                  //
       ////////////////////////////////////////////////////////////////////

    } else if (strcmp(IP.Next_Control_Parameter,"Boundary_Conditions_Specified") == 0) {
      i_command = 500;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Boundary_Conditions_Specified,IP.Next_Control_Parameter);
      if (strcmp(IP.Boundary_Conditions_Specified,"On") == 0) {
	IP.BCs_Specified = ON;
      } else if (strcmp(IP.Boundary_Conditions_Specified,"Off") == 0) {
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
      } else if (strcmp(IP.BC_North_Type,"Fixed") == 0) {
	IP.BC_North = BC_FIXED;
      } else if (strcmp(IP.BC_North_Type,"Inflow_Subsonic") == 0) {
	IP.BC_North = BC_INFLOW_SUBSONIC;
      } else if (strcmp(IP.BC_North_Type,"Outflow_Subsonic") == 0) {
	IP.BC_North = BC_OUTFLOW_SUBSONIC;
      } else if (strcmp(IP.BC_North_Type,"Fixed_Pressure") == 0) {
	IP.BC_North = BC_FIXED_PRESSURE;
      } else if (strcmp(IP.BC_North_Type,"Constant_Extrapolation") == 0) {
	IP.BC_North = BC_CONSTANT_EXTRAPOLATION;
      } else if (strcmp(IP.BC_North_Type,"Linear_Extrapolation") == 0) {
	IP.BC_North = BC_LINEAR_EXTRAPOLATION;
      } else if (strcmp(IP.BC_North_Type,"Characteristic") == 0) {
	IP.BC_North = BC_CHARACTERISTIC;
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
      } else if (strcmp(IP.BC_South_Type,"Fixed") == 0) {
	IP.BC_South = BC_FIXED;
      } else if (strcmp(IP.BC_South_Type,"Inflow_Subsonic") == 0) {
	IP.BC_South = BC_INFLOW_SUBSONIC;
      } else if (strcmp(IP.BC_South_Type,"Outflow_Subsonic") == 0) {
	IP.BC_South = BC_OUTFLOW_SUBSONIC;
      } else if (strcmp(IP.BC_South_Type,"Fixed_Pressure") == 0) {
	IP.BC_South = BC_FIXED_PRESSURE;
      } else if (strcmp(IP.BC_South_Type,"Constant_Extrapolation") == 0) {
	IP.BC_South = BC_CONSTANT_EXTRAPOLATION;
      } else if (strcmp(IP.BC_South_Type,"Linear_Extrapolation") == 0) {
	IP.BC_South = BC_LINEAR_EXTRAPOLATION;
      } else if (strcmp(IP.BC_South_Type,"Characteristic") == 0) {
	IP.BC_South = BC_CHARACTERISTIC;
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
      } else if (strcmp(IP.BC_East_Type,"Fixed") == 0) {
	IP.BC_East = BC_FIXED;
      } else if (strcmp(IP.BC_East_Type,"Inflow_Subsonic") == 0) {
	IP.BC_East = BC_INFLOW_SUBSONIC;
      } else if (strcmp(IP.BC_East_Type,"Outflow_Subsonic") == 0) {
	IP.BC_East = BC_OUTFLOW_SUBSONIC;
      } else if (strcmp(IP.BC_East_Type,"Fixed_Pressure") == 0) {
	IP.BC_East = BC_FIXED_PRESSURE;
      } else if (strcmp(IP.BC_East_Type,"Constant_Extrapolation") == 0) {
	IP.BC_East = BC_CONSTANT_EXTRAPOLATION;
      } else if (strcmp(IP.BC_East_Type,"Linear_Extrapolation") == 0) {
	IP.BC_East = BC_LINEAR_EXTRAPOLATION;
      } else if (strcmp(IP.BC_East_Type,"Characteristic") == 0) {
	IP.BC_East = BC_CHARACTERISTIC;
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
      } else if (strcmp(IP.BC_West_Type,"Fixed") == 0) {
	IP.BC_West = BC_FIXED;
      } else if (strcmp(IP.BC_West_Type,"Inflow_Subsonic") == 0) {
	IP.BC_West = BC_INFLOW_SUBSONIC;
      } else if (strcmp(IP.BC_West_Type,"Outflow_Subsonic") == 0) {
	IP.BC_West = BC_OUTFLOW_SUBSONIC;
      } else if (strcmp(IP.BC_West_Type,"Fixed_Pressure") == 0) {
	IP.BC_West = BC_FIXED_PRESSURE;
      } else if (strcmp(IP.BC_West_Type,"Constant_Extrapolation") == 0) {
	IP.BC_West = BC_CONSTANT_EXTRAPOLATION;
      } else if (strcmp(IP.BC_West_Type,"Linear_Extrapolation") == 0) {
	IP.BC_West = BC_LINEAR_EXTRAPOLATION;
      } else if (strcmp(IP.BC_West_Type,"Characteristic") == 0) {
	IP.BC_West = BC_CHARACTERISTIC;
      } else if (strcmp(IP.BC_West_Type,"None") == 0) {
	IP.BC_West = BC_NONE;
      } else {
	i_command = INVALID_INPUT_VALUE;
      }

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

    } else if (IP.Next_Control_Parameter[0] == '#') {
       i_command = COMMENT_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Refine_Grid") == 0) {
       i_command = REFINE_GRID_CODE;

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
int Process_Input_Control_Parameter_File(Ion5Moment2D_Input_Parameters &Input_Parameters,
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
       cout << "\n Ion5Moment2D ERROR: Unable to open Ion5Moment2D input data file.\n";
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
          cout << "\n Ion5Moment2D ERROR: Error reading Ion5Moment2D data at line #"
               << -line_number  << " of input data file.\n";
          error_flag = line_number;
          break;
       } /* endif */
    } /* endwhile */

    /* Perform consistency checks on Input_Parameters */

    if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {
      error_flag = Check_Input_Parameters<Ion5Moment2D_Input_Parameters>(Input_Parameters);
      if (error_flag) {
	cout << "\n Ion5Moment2D ERROR: Input Parameters consistency check failure\n";
	return (error_flag);
      }
    }

    // Perform consitency checks on the refinement criteria.
    Input_Parameters.Number_of_Refinement_Criteria = 0;
    if (Input_Parameters.Refinement_Criteria_Gradient_Density)
      Input_Parameters.Number_of_Refinement_Criteria++;
    if (Input_Parameters.Refinement_Criteria_Divergence_Velocity)
      Input_Parameters.Number_of_Refinement_Criteria++;
    if (Input_Parameters.Refinement_Criteria_Curl_Velocity)
      Input_Parameters.Number_of_Refinement_Criteria++;
    if (Input_Parameters.Number_of_Refinement_Criteria < 1 ||
	Input_Parameters.Number_of_Refinement_Criteria > 3) return 1011;

    /* Initial processing of input control parameters complete.  
       Return the error indicator flag. */

    return (error_flag);

}
