/* Reconstruct3DInput.cc:  Subroutines for 
                           3D Reconstuct Input Classes. */

/* Include Reconstruct 3D input parameter header file. */

#ifndef _RECONSTRUCT3D_INPUT_INCLUDED
#include "Reconstruct3DInput.h"
#include "../CENO_Tolerances.h" // Include high-order CENO tolerances header file
#endif // _RECONSTRUCT3D_INPUT_INCLUDED

/***************************************************************
 * Reconstruct3D_Input_Parameters -- Input-output operators. *
 ***************************************************************/
ostream &operator << (ostream &out_file,
			     const Reconstruct3D_Input_Parameters &IP) {

  if (IP.TestF != NULL){
    out_file << setprecision(6);
    out_file << "\n\n Solving 3D Reconstruction for the function \""
	     << IP.Function_Type << "\"";
    out_file << "\n  -> Method: "
	     << IP.Method_Used;

    switch (IP.Method){
    case DD_ENO:
      out_file << "\n  -> Cutoff Knob = "
	       << IP.CutoffKnob();
      break;
    case ENO:
      // comments
      break;
    case WENO:
      // comments
      break;
    case SpectralDiff:
      // comments
      break;

    case CENO:
      out_file << "\n  -> Limiter: "
	       << IP.Limiter_Type; 
      out_file << "\n  -> Fit Tolerance: "
	       << CENO_Tolerances::Fit_Tolerance;
      break;
    } /* endswitch */

    out_file << "\n  -> Reconstruction Order: "
	     << IP.Reconstruction_Order;
    out_file << "\n  -> Input File Name: " 
             << IP.Input_File_Name;
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
        out_file << "\n  -> Length of Grain (m): " 
                 << IP.Grain_Length;
        out_file << "\n  -> Radius of Grain (m): " 
                 << IP.Grain_Radius;
        out_file << "\n  -> Distance from Grain to Nozzle Throat (m): " 
                 << IP.Grain_To_Throat_Length;
        out_file << "\n  -> Length of the Nozzle (m): " 
                 << IP.Nozzle_Length;
        out_file << "\n  -> Radius of the Nozzle at Throat (m): " 
                 << IP.Nozzle_Radius_Throat;
        out_file << "\n  -> Radius of the Nozzle at Exit (m): " 
                 << IP.Nozzle_Radius_Exit;
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
      case GRID_UNSTEADY_BLUNT_BODY :
	out_file << "\n  -> Cylinder Radius (m): " 
		 << IP.Blunt_Body_Radius;
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
    out_file << "\n  -> Mesh shift, scale, and rotate: " 
             << IP.X_Shift << " " << IP.X_Scale << " " << IP.X_Rotate;
    out_file << "\n  -> Mesh Stretching Function I_Direction: "
	     << IP.Stretching_Function_I
	     << "\n  -> Mesh Stretching Function J_Direction: "
	     << IP.Stretching_Function_J
	     << "\n  -> Stretching Parameters: Beta_I, Beta_J, Tau_I, Tau_J: "
	     << IP.Beta_I << " " << IP.Beta_J << " " << IP.Tau_I << " " << IP.Tau_J;
    out_file << "\n  -> Number of iteration for unsmoothing the mesh: " << IP.NumOfIter_UnsmoothMesh;
    out_file << "\n  -> Number of Blocks i-direction: "
             << IP.Number_of_Blocks_Idir;
    out_file << "\n  -> Number of Blocks j-direction: " 
             << IP.Number_of_Blocks_Jdir;
    out_file << "\n  -> Number of Cells i-direction: "
             << IP.Number_of_Cells_Idir;
    out_file << "\n  -> Number of Cells j-direction: " 
             << IP.Number_of_Cells_Jdir;
    out_file << "\n  -> Number of SubGrid Points in X-dir: "
	     << IP.Number_of_SubGrid_Points_Idir;
    out_file << "\n  -> Number of SubGrid Points in Y-dir: "
	     << IP.Number_of_SubGrid_Points_Jdir;
    out_file << "\n  -> Output File Name: " 
             << IP.Output_File_Name;

    out_file << "\n";
  } else {



  }

  return (out_file);
}

istream &operator >> (istream &in_file,
			     Reconstruct3D_Input_Parameters &IP) {
    return (in_file);
}

/*************************************************************
 * _Reconstruct3D_Input_Parameters -- External subroutines. *
 *************************************************************/

/********************************************************
 * Routine: Open_Input_File                             *
 *                                                      *
 * Opens the appropriate input data file.               *
 *                                                      *
 ********************************************************/
void Open_Input_File(Reconstruct3D_Input_Parameters &IP) {

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
void Close_Input_File(Reconstruct3D_Input_Parameters &IP) {

    IP.Input_File.unsetf(ios::skipws);
    IP.Input_File.close();
}

/********************************************************
 * Routine: Set_Default_Input_Parameters                *
 *                                                      *
 * Assigns default values to the input parameters.      *
 *                                                      *
 ********************************************************/
void Set_Default_Input_Parameters(Reconstruct3D_Input_Parameters &IP) {

    int i;

    strcpy(IP.Input_File_Name, "reconstruct3D.in");
    // Method settings
    strcpy(IP.Method_Used, "Normal_Equation");
    IP.Method = DD_ENO;
    strcpy(IP.Geometric_Weighting, "No");
    IP.geom_weighting = OFF;
    strcpy(IP.Data_Dependent_Weighting, "No");
    IP.data_depend_weighting = OFF;

    strcpy(IP.Limiter_Type, "Barth-Jespersen");
    IP.i_Limiter = LIMITER_BARTH_JESPERSEN;
    IP.CENO_Cutoff = 430;

    // Function settings
    strcpy(IP.Function_Type, "Function_Default");
    IP.i_Function = FUNCTION_DEFAULT;
    IP.TestF = Test_Default3D;
    IP.IntTestF = Test_Default3D_Integral;
    IP.Reconstruction_Order = 1;

    // Grid settings
    strcpy(IP.Grid_Type, "Cube");
    IP.i_Grid = GRID_CUBE;
    IP.Box_Width = ONE;
    IP.Box_Height = ONE;
    IP.Box_Length = ONE;
    IP.Number_of_Cells_Idir = 10;
    IP.Number_of_Cells_Jdir = 10;
    IP.Number_of_Cells_Kdir = 10;
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Kdir = 1;
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

    IP.X_Shift = Vector3D_ZERO;
    IP.X_Scale = ONE;
    IP.X_Rotate = ZERO;
    IP.CharacteristicLength = ONE;
    IP.CutoffKnob() = ONE;
    strcpy(IP.Stretching_Function_I, "Linear_Function");
    strcpy(IP.Stretching_Function_J, "Linear_Function"); 
    strcpy(IP.Stretching_Function_K, "Linear_Function"); 
    IP.Stretch_I = STRETCHING_FCN_LINEAR;
    IP.Beta_I = ZERO; 
    IP.Tau_I = ZERO;
    IP.Stretch_J = STRETCHING_FCN_LINEAR;
    IP.Beta_J = ZERO;
    IP.Tau_J = ZERO;
    IP.Stretch_K = STRETCHING_FCN_LINEAR;
    IP.Beta_K = ZERO;
    IP.Tau_K = ZERO;
    IP.NumOfIter_UnsmoothMesh = 0;


    IP.Number_of_SubGrid_Points_Idir = 5;
    IP.Number_of_SubGrid_Points_Jdir = 5;
    IP.Number_of_SubGrid_Points_Kdir = 5;

    strcpy(IP.Output_File_Name, "outputfile.dat");

    strcpy(IP.Grid_File_Name, "gridfile.grid");
    strcpy(IP.Grid_Definition_File_Name, "gridfile.griddef");

    strcpy(IP.Gnuplot_File_Name, "gnuplotfile.gplt");

    strcpy(IP.Output_Format_Type, "Tecplot");
    IP.i_Output_Format = IO_TECPLOT;

    strcpy(IP.Next_Control_Parameter, " ");

    IP.Line_Number = 0;
    IP.Message_Number = 1;

}

/********************************************************
 * Routine: Get_Next_Input_Control_Parameter            *
 *                                                      *
 * Get the next input control parameter from the input  *
 * file.                                                *
 *                                                      *
 ********************************************************/
void Get_Next_Input_Control_Parameter(Reconstruct3D_Input_Parameters &IP) {

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
int Parse_Next_Input_Control_Parameter(Reconstruct3D_Input_Parameters &IP) {

    int i_command;
    char buffer[256];

    i_command = 0;

    if (strcmp(IP.Next_Control_Parameter, "Function_Definition") == 0) {
       i_command = 1;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Function_Type, 
              IP.Next_Control_Parameter);

       // Add the parse word for the new function //
       if (strcmp(IP.Function_Type, "Read_Data_From_File") == 0){
	 IP.TestF = NULL;
	 IP.IntTestF = NULL;
       }
       else if (strcmp(IP.Function_Type, "Default") == 0){
	 IP.TestF = Test_Default3D;
	 IP.IntTestF = Test_Default3D_Integral;
       }
       else if (strcmp(IP.Function_Type, "Example1") == 0){
	 IP.TestF = Test_Example1;
	 IP.IntTestF = Test_Example1_Integral;
       }
       else if (strcmp(IP.Function_Type, "Example2") == 0){
	 IP.TestF = Test_Example2;
	 IP.IntTestF = Test_Example2_Integral;
       }
       else if (strcmp(IP.Function_Type, "Example3") == 0){
	 IP.TestF = Test_Example3;
	 IP.IntTestF = Test_Example3_Integral;
       }
       else if (strcmp(IP.Function_Type, "Example4") == 0){
	 IP.TestF = Test_Example4;
	 IP.IntTestF = Test_Example4_Integral;
       }
       else if (strcmp(IP.Function_Type, "Example5") == 0){
	 IP.TestF = Test_Example5;
	 IP.IntTestF = Test_Example5_Integral;
       }
       else if (strcmp(IP.Function_Type, "Example6") == 0){
	 IP.TestF = Test_Example6;
	 IP.IntTestF = Test_Example6_Integral;
       }
       else if (strcmp(IP.Function_Type, "Example7") == 0){
	 IP.TestF = Test_Example7;
	 IP.IntTestF = Test_Example7_Integral;
       }
       else if (strcmp(IP.Function_Type, "Example8") == 0){
	 IP.TestF = Test_Example8;
	 IP.IntTestF = Test_Example8_Integral;
       }
       // endif 

    } else if (strcmp(IP.Next_Control_Parameter, "Reconstruction_Order") == 0){
	i_command = 2;
	++IP.Line_Number;

	if (!(IP.Input_File >> IP.Reconstruction_Order))
	  i_command = INVALID_INPUT_VALUE;
	IP.Input_File.getline(buffer, sizeof(buffer));
	if (IP.Reconstruction_Order < ZERO) {
	    IP.Reconstruction_Order = 0;
	    cout << "\n MESSAGE "<< IP.Message_Number
		 << ": Reconstruction Order must be at least 0.\a";
	    IP.Message_Number++;
	};
	if (IP.Reconstruction_Order > SEVEN){
	    IP.Reconstruction_Order = 7;
	    cout << "\n MESSAGE "<< IP.Message_Number
		 << ": This program doesn't perform reconstruction"
		 << " higher than 7.\a";
	    IP.Message_Number++;
	} // endif

	// Set the required number of Ghost cells
	IP.Number_of_Ghost_Cells = 5;

    } else if (strcmp(IP.Next_Control_Parameter, "Grid_Type") == 0) {
       i_command = 3;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Grid_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Grid_Type, "Cube") == 0) {
          IP.i_Grid = GRID_CUBE;
          IP.Box_Width  = ONE;
          IP.Box_Height = ONE;
	  IP.Box_Length = ONE;
       } else {
	 IP.i_Grid = GRID_CUBE;
          IP.Box_Width  = ONE;
          IP.Box_Height = ONE;
	  IP.Box_Length = ONE;
       } /* endif */

       /* ********** 2D Grid Types ******************************

       else if (strcmp(IP.Grid_Type, "Square") == 0) {
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
       } else if (strcmp(IP.Grid_Type, "Rocket_Cold_Flow") == 0) {
	  IP.i_Grid = GRID_ROCKET_MOTOR_COLD_FLOW ;
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
       } else if (strcmp(IP.Grid_Type,"Ringleb_Flow") == 0) {
          IP.i_Grid = GRID_RINGLEB_FLOW;
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
       } 
       **************** End of 2D Grid Types ************************* */

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Cells_Idir") == 0) {
       i_command = 4;
       ++IP.Line_Number;
       IP.Input_File >> IP.Number_of_Cells_Idir;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Cells_Idir < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Cells_Jdir") == 0) {
       i_command = 5;
       ++IP.Line_Number;
       IP.Input_File >> IP.Number_of_Cells_Jdir;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Cells_Jdir < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Cells_Kdir") == 0) {
       i_command = 5;
       ++IP.Line_Number;
       IP.Input_File >> IP.Number_of_Cells_Kdir;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Cells_Kdir < 1) i_command = INVALID_INPUT_VALUE;

    }else if (strcmp(IP.Next_Control_Parameter, "Number_of_Blocks_Idir") == 0) {
       i_command = 6;
       ++IP.Line_Number;
       IP.Input_File >> IP.Number_of_Blocks_Idir;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Blocks_Idir < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Blocks_Jdir") == 0) {
       i_command = 7;
       ++IP.Line_Number;
       IP.Input_File >> IP.Number_of_Blocks_Jdir;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Blocks_Jdir < 1) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Blocks_Kdir") == 0) {
       i_command = 7;
       ++IP.Line_Number;
       IP.Input_File >> IP.Number_of_Blocks_Kdir;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Number_of_Blocks_Kdir < 1) i_command = INVALID_INPUT_VALUE;

    }else if (strcmp(IP.Next_Control_Parameter, "Box_Width") == 0) {
       i_command = 8;
       ++IP.Line_Number;
       IP.Input_File >> IP.Box_Width;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Box_Width <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Box_Height") == 0) {
       i_command = 9;
       ++IP.Line_Number;
       IP.Input_File >> IP.Box_Height;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Box_Height <= ZERO) i_command = INVALID_INPUT_VALUE;

    }else if (strcmp(IP.Next_Control_Parameter, "Box_Length") == 0) {
       i_command = 9;
       ++IP.Line_Number;
       IP.Input_File >> IP.Box_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Box_Height <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Plate_Length") == 0) {
       i_command = 10;
       ++IP.Line_Number;
       IP.Input_File >> IP.Plate_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Plate_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Pipe_Length") == 0) {
       i_command = 11;
       ++IP.Line_Number;
       IP.Input_File >> IP.Pipe_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Pipe_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Pipe_Radius") == 0) {
       i_command = 12;
       ++IP.Line_Number;
       IP.Input_File >> IP.Pipe_Radius;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Pipe_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Blunt_Body_Radius") == 0) {
       i_command = 13;
       ++IP.Line_Number;
       IP.Input_File >> IP.Blunt_Body_Radius;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Blunt_Body_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Blunt_Body_Mach_Number") == 0) {
       i_command = 14;
       ++IP.Line_Number;
       IP.Input_File >> IP.Blunt_Body_Mach_Number;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Blunt_Body_Mach_Number <= ONE) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Cylinder_Radius") == 0) {
       i_command = 15;
       ++IP.Line_Number;
       IP.Input_File >> IP.Cylinder_Radius;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Cylinder_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Ellipse_Length_X_Axis") == 0) {
       i_command = 16;
       ++IP.Line_Number;
       IP.Input_File >> IP.Ellipse_Length_X_Axis;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Ellipse_Length_X_Axis <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Ellipse_Length_Y_Axis") == 0) {
       i_command = 17;
       ++IP.Line_Number;
       IP.Input_File >> IP.Ellipse_Length_Y_Axis;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Ellipse_Length_Y_Axis <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Chord_Length") == 0) {
       i_command = 18;
       ++IP.Line_Number;
       IP.Input_File >> IP.Chord_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Chord_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "NACA_Aerofoil_Type") == 0) {
       i_command = 19;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.NACA_Aerofoil_Type, 
              IP.Next_Control_Parameter);
       if (strlen(IP.NACA_Aerofoil_Type) != 4 &&
           strlen(IP.NACA_Aerofoil_Type) != 5) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Orifice_Radius") == 0) {
       i_command = 20;
       ++IP.Line_Number;
       IP.Input_File >> IP.Orifice_Radius;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Orifice_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Wedge_Angle") == 0) {
       i_command = 21;
       ++IP.Line_Number;
       IP.Input_File >> IP.Wedge_Angle;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Wedge_Angle <= ZERO) i_command = INVALID_INPUT_VALUE;
    
    } else if (strcmp(IP.Next_Control_Parameter, "Wedge_Length") == 0) {
       i_command = 22;
       ++IP.Line_Number;
       IP.Input_File >> IP.Wedge_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Wedge_Length <= ZERO) i_command = INVALID_INPUT_VALUE;
    
    } else if (strcmp(IP.Next_Control_Parameter, "Grain_Length") == 0) {
       i_command = 23;
       ++IP.Line_Number;
       IP.Input_File >> IP.Grain_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Grain_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Grain_Radius") == 0) {
       i_command = 24;
       ++IP.Line_Number;
       IP.Input_File >> IP.Grain_Radius;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Grain_Radius <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Grain_To_Throat_Length") == 0) {
       i_command = 25;
       ++IP.Line_Number;
       IP.Input_File >> IP.Grain_To_Throat_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Grain_To_Throat_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Nozzle_Length") == 0) {
       i_command = 26;
       ++IP.Line_Number;
       IP.Input_File >> IP.Nozzle_Length;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Nozzle_Length <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Nozzle_Radius_Exit") == 0) {
       i_command = 27;
       ++IP.Line_Number;
       IP.Input_File >> IP.Nozzle_Radius_Exit;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Nozzle_Radius_Exit <= ZERO) i_command = INVALID_INPUT_VALUE;

    } else if (strcmp(IP.Next_Control_Parameter, "Nozzle_Radius_Throat") == 0) {
       i_command = 28;
       ++IP.Line_Number;
       IP.Input_File >> IP.Nozzle_Radius_Throat;
       IP.Input_File.getline(buffer, sizeof(buffer));
       if (IP.Nozzle_Radius_Throat <= ZERO) i_command = INVALID_INPUT_VALUE;
  
    } else if (strcmp(IP.Next_Control_Parameter, "Output_File_Name") == 0) {
      i_command = 29;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Output_File_Name, 
	     IP.Next_Control_Parameter);
      strcpy(IP.Grid_File_Name, 
	     IP.Next_Control_Parameter);
      strcpy(IP.Grid_Definition_File_Name, 
	     IP.Next_Control_Parameter);
      strcpy(IP.Gnuplot_File_Name, 
	     IP.Next_Control_Parameter);

    } else if (strcmp(IP.Next_Control_Parameter, "Input_File_Name") == 0) {
      i_command = 29;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Input_File_Name, 
	     IP.Next_Control_Parameter);
      
    } else if (strcmp(IP.Next_Control_Parameter, "Grid_File_Name") == 0) {
      i_command = 30;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Grid_File_Name, 
	     IP.Next_Control_Parameter);
      strcat(IP.Grid_File_Name, ".grid");
      strcpy(IP.Grid_Definition_File_Name, 
	     IP.Next_Control_Parameter);
      strcat(IP.Grid_Definition_File_Name, ".griddef");
      
    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Subgrid_Points_Idir") == 0){
      i_command = 31;
      ++IP.Line_Number;
      if ((IP.Input_File >> IP.Number_of_SubGrid_Points_Idir) == 0 ||
	  (IP.Number_of_SubGrid_Points_Idir < 1))
	i_command = INVALID_INPUT_VALUE;
      IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Subgrid_Points_Jdir") == 0){
      i_command = 32;
      ++IP.Line_Number;
      if ((IP.Input_File >> IP.Number_of_SubGrid_Points_Jdir) == 0 ||
	  (IP.Number_of_SubGrid_Points_Jdir < 1))
	i_command = INVALID_INPUT_VALUE;
      IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Subgrid_Points_Kdir") == 0){
      i_command = 32;
      ++IP.Line_Number;
      if ((IP.Input_File >> IP.Number_of_SubGrid_Points_Kdir) == 0 ||
	  (IP.Number_of_SubGrid_Points_Kdir < 1))
	i_command = INVALID_INPUT_VALUE;
      IP.Input_File.getline(buffer, sizeof(buffer));

    }else if (strcmp(IP.Next_Control_Parameter, "Method_Used") == 0) {
      i_command = 33;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Method_Used, 
	     IP.Next_Control_Parameter);

       // Add the parse word for the new function //
       if (strcmp(IP.Method_Used, "ENO") == 0)
	 IP.Method = ENO;
       else if (strcmp(IP.Method_Used, "DD_ENO") == 0)
	 IP.Method = DD_ENO;
       else if (strcmp(IP.Method_Used, "Spectral") == 0)
	 IP.Method = SpectralDiff;
       else if (strcmp(IP.Method_Used, "WENO") == 0)
	 IP.Method = WENO;
       else if (strcmp(IP.Method_Used, "ENO_LS") == 0)
	 IP.Method = ENO_LS;
       else if (strcmp(IP.Method_Used, "CENO") == 0)
	 IP.Method = CENO;

    } else if (strcmp(IP.Next_Control_Parameter, "Geometric_Weighting") == 0){
      i_command = 34;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Geometric_Weighting, 
	     IP.Next_Control_Parameter);

       // Add the parse word for the new function //
       if (strcmp(IP.Geometric_Weighting, "Yes") == 0)
	 IP.geom_weighting = ON;
       else if (strcmp(IP.Geometric_Weighting, "No") == 0)
	 IP.geom_weighting = OFF;
      
    } else if (strcmp(IP.Next_Control_Parameter, 
		      "Data_Dependent_Weighting") == 0){
      i_command = 35;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Data_Dependent_Weighting, 
	     IP.Next_Control_Parameter);

       // Add the parse word for the new function //
       if (strcmp(IP.Data_Dependent_Weighting, "No") == 0)
	 IP.data_depend_weighting = OFF;
       else if ((strcmp(IP.Data_Dependent_Weighting, "Yes") == 0) &&
		(IP.geom_weighting))
	 IP.data_depend_weighting = ON;
     
      
    } else if (strcmp(IP.Next_Control_Parameter, "Output_Format_Type") == 0) {
      i_command = 36;
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
      } // endif 

    } else if (strcmp(IP.Next_Control_Parameter, "X_Shift") == 0) {
       i_command = 37;
       ++IP.Line_Number;
       IP.Input_File >> IP.X_Shift;
       IP.Input_File.setf(ios::skipws);
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "X_Scale") == 0) {
       i_command = 38;
       ++IP.Line_Number;
       IP.Input_File >> IP.X_Scale;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "X_Rotate") == 0) {
       i_command = 39;
       ++IP.Line_Number;
       IP.Input_File >> IP.X_Rotate;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Characteristic_Length") == 0) {
       i_command = 40;
       ++IP.Line_Number;
       IP.Input_File >> IP.CharacteristicLength;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Cutoff_Knob") == 0) {
       i_command = 41;
       ++IP.Line_Number;
       IP.Input_File >> IP.CutoffKnob();
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Integration_Type") == 0){
      i_command = 42;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Integration_Type,
	     IP.Next_Control_Parameter);

    } else if (strcmp(IP.Next_Control_Parameter, "Limiter") == 0){
      i_command = 43;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Limiter_Type, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Limiter_Type, "One") == 0) {
          IP.i_Limiter = LIMITER_ONE;
       } else if (strcmp(IP.Limiter_Type, "Zero") == 0) {
          IP.i_Limiter = LIMITER_ZERO;
       } else if (strcmp(IP.Limiter_Type, "Minmod") == 0) {
          IP.i_Limiter = LIMITER_MINMOD;
       } else if (strcmp(IP.Limiter_Type, "UMIST") == 0) {
          IP.i_Limiter = LIMITER_UMIST;
       } else if (strcmp(IP.Limiter_Type, "Double_Minmod") == 0) {
          IP.i_Limiter = LIMITER_DOUBLE_MINMOD;
       } else if (strcmp(IP.Limiter_Type, "Superbee") == 0) {
          IP.i_Limiter = LIMITER_SUPERBEE;
       } else if (strcmp(IP.Limiter_Type, "Phi") == 0) {
          IP.i_Limiter = LIMITER_PHI;
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

    } else if (strcmp(IP.Next_Control_Parameter, "Stretch_I") == 0) {
       i_command = 45;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Stretching_Function_I, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Stretching_Function_I, "Linear") == 0) {
	 IP.Stretch_I = STRETCHING_FCN_LINEAR;
       } else if (strcmp(IP.Stretching_Function_I, "Min_Clustering") == 0) {
	 IP.Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
       } else if (strcmp(IP.Stretching_Function_I, "Max_Clustering") == 0) {
	 IP.Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
       } else if (strcmp(IP.Stretching_Function_I, "MinMax_Clustering") == 0) {
	 IP.Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
       } else if (strcmp(IP.Stretching_Function_I, "Midpt_Clustering") == 0) {
	 IP.Stretch_I = STRETCHING_FCN_MIDPT_CLUSTERING;
       } else if (strcmp(IP.Stretching_Function_I, "Sine") == 0) {
	 IP.Stretch_I = STRETCHING_FCN_SINE;
       } else if (strcmp(IP.Stretching_Function_I, "Cosine") == 0) {
	 IP.Stretch_I = STRETCHING_FCN_COSINE;
       } else {
       	 IP.Stretch_I = STRETCHING_FCN_LINEAR;
       }

    } else if (strcmp(IP.Next_Control_Parameter, "Stretch_J") == 0) {
       i_command = 45;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Stretching_Function_J, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Stretching_Function_J, "Linear") == 0) {
	 IP.Stretch_J = STRETCHING_FCN_LINEAR;
       } else if (strcmp(IP.Stretching_Function_J, "Min_Clustering") == 0) {
	 IP.Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
       } else if (strcmp(IP.Stretching_Function_J, "Max_Clustering") == 0) {
	 IP.Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
       } else if (strcmp(IP.Stretching_Function_J, "MinMax_Clustering") == 0) {
	 IP.Stretch_J = STRETCHING_FCN_MINMAX_CLUSTERING;
       } else if (strcmp(IP.Stretching_Function_J, "Midpt_Clustering") == 0) {
	 IP.Stretch_J = STRETCHING_FCN_MIDPT_CLUSTERING;
       } else if (strcmp(IP.Stretching_Function_J, "Sine") == 0) {
	 IP.Stretch_J = STRETCHING_FCN_SINE;
       } else if (strcmp(IP.Stretching_Function_J, "Cosine") == 0) {
	 IP.Stretch_J = STRETCHING_FCN_COSINE;
       } else {
       	 IP.Stretch_J = STRETCHING_FCN_LINEAR;
       }

    } else if (strcmp(IP.Next_Control_Parameter, "Stretch_K") == 0) {
       i_command = 45;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Stretching_Function_K, 
              IP.Next_Control_Parameter);
       if (strcmp(IP.Stretching_Function_K, "Linear") == 0) {
	 IP.Stretch_K = STRETCHING_FCN_LINEAR;
       } else if (strcmp(IP.Stretching_Function_K, "Min_Clustering") == 0) {
	 IP.Stretch_K = STRETCHING_FCN_MIN_CLUSTERING;
       } else if (strcmp(IP.Stretching_Function_K, "Max_Clustering") == 0) {
	 IP.Stretch_K = STRETCHING_FCN_MAX_CLUSTERING;
       } else if (strcmp(IP.Stretching_Function_K, "MinMax_Clustering") == 0) {
	 IP.Stretch_K = STRETCHING_FCN_MINMAX_CLUSTERING;
       } else if (strcmp(IP.Stretching_Function_K, "Midpt_Clustering") == 0) {
	 IP.Stretch_K = STRETCHING_FCN_MIDPT_CLUSTERING;
       } else if (strcmp(IP.Stretching_Function_K, "Sine") == 0) {
	 IP.Stretch_K = STRETCHING_FCN_SINE;
       } else if (strcmp(IP.Stretching_Function_K, "Cosine") == 0) {
	 IP.Stretch_K = STRETCHING_FCN_COSINE;
       } else {
       	 IP.Stretch_K = STRETCHING_FCN_LINEAR;
       }

    }else if (strcmp(IP.Next_Control_Parameter, "Beta_I") == 0) {
       i_command = 46;
       ++IP.Line_Number;
       IP.Input_File >> IP.Beta_I;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Beta_J") == 0) {
       i_command = 46;
       ++IP.Line_Number;
       IP.Input_File >> IP.Beta_J;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Beta_K") == 0) {
       i_command = 46;
       ++IP.Line_Number;
       IP.Input_File >> IP.Beta_K;
       IP.Input_File.getline(buffer, sizeof(buffer));

    }else if (strcmp(IP.Next_Control_Parameter, "Tau_I") == 0) {
       i_command = 46;
       ++IP.Line_Number;
       IP.Input_File >> IP.Tau_I;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Tau_J") == 0) {
       i_command = 46;
       ++IP.Line_Number;
       IP.Input_File >> IP.Tau_J;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Tau_K") == 0) {
       i_command = 46;
       ++IP.Line_Number;
       IP.Input_File >> IP.Tau_K;
       IP.Input_File.getline(buffer, sizeof(buffer));

    }else if (strcmp(IP.Next_Control_Parameter, "Iter_UnsmoothMesh") == 0) {
       i_command = 47;
       ++IP.Line_Number;
       IP.Input_File >> IP.NumOfIter_UnsmoothMesh;
       IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Execute") == 0) {
      i_command = EXECUTE_CODE;
      
    } else if (strcmp(IP.Next_Control_Parameter, "Terminate") == 0) {
      i_command = TERMINATE_CODE;
      
    } else if (strcmp(IP.Next_Control_Parameter, "Continue") == 0) {
      i_command = CONTINUE_CODE;
      
    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output") == 0) {
      i_command = WRITE_OUTPUT_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_OneZone") == 0) {
       i_command = WRITE_UNIFIED_OUTPUT_CODE;
      
    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Cells") == 0) {
      i_command = WRITE_OUTPUT_CELLS_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Stencil_Reconstruction") == 0) {
      i_command = WRITE_RECONSTRUCTED_FUNCTIONS;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Mesh") == 0) {
      i_command = WRITE_OUTPUT_GRID_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Mesh_Definition") == 0) {
      i_command = WRITE_GRID_DEFINITION_CODE;
      
    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Solution_Nodes") == 0) {
      i_command = WRITE_OUTPUT_GRID_NODES_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Full_Solution_Nodes") == 0) {
      i_command = WRITE_OUTPUT_FULL_GRID_NODES_CODE;
      
    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Mesh_Cells") == 0) {
      i_command = WRITE_OUTPUT_GRID_CELLS_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Accuracy") == 0) {
      i_command = WRITE_OUTPUT_ACCURACY_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Print_Norms") == 0) {
      i_command = WRITE_NORM_ON_SCREEN;

    } else if (strcmp(IP.Next_Control_Parameter, "Assess_Accuracy") == 0) {
      i_command = ASSESS_ACCURACY_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Function_Graph") == 0) {
      i_command = WRITE_OUTPUT_FUNCTION_GRAPH_CODE;

    } else if (IP.Next_Control_Parameter[0] == '#') {
      i_command = COMMENT_CODE;
      
    } else {
      i_command = INVALID_INPUT_CODE;
      
    } /* endif */

    /* Parse next control parameter with CENO_Tolerances parser */
    CENO_Tolerances::Parse_Next_Input_Control_Parameter(IP,i_command);

    
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
int Process_Input_Control_Parameter_File(Reconstruct3D_Input_Parameters
					 &Input_Parameters,
                                         char *Input_File_Name_ptr,
                                         int &Command_Flag) {

    int error_flag, line_number;

    /* Assign initial value for error indicator flag. */

    error_flag = 0;
 
    /* Assign default values to the input parameters. */

    Set_Default_Input_Parameters(Input_Parameters);

    /* Copy input file name (a string) to appropriate input parameter variable. */

    if (Input_File_Name_ptr != NULL) strcpy(Input_Parameters.Input_File_Name,
					    Input_File_Name_ptr);

    /* Open the input file containing the input parameters. */

    Open_Input_File(Input_Parameters);
    error_flag = Input_Parameters.Input_File.bad();

    if (error_flag) {
       cout << "\n Reconstruction3D ERROR: Unable to open Reconstruction3D input"
	    << " data file.\n";
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
          cout << "\n Reconstruction3D ERROR: Error reading Reconstruction3D"
	       << " data at line #"
               << -line_number  << " of input data file.\n\a";
          error_flag = line_number;
          break;
       } /* endif */
    } /* endwhile */

    /* Initial processing of input control parameters complete.  
       Return the error indicator flag. */

    return (error_flag);
}
