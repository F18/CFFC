/* Reconstruct1DInput.cc:  Subroutines for 
                           1D Reconstuct Input Classes. */

/* Include Reconstruct 1D input parameter header file. */

#ifndef _RECONSTRUCT1D_INPUT_INCLUDED
#include "Reconstruct1DInput.h"
#endif // _RECONSTRUCT1D_INPUT_INCLUDED 

/***************************************************************
 * Reconstruct1D_Input_Parameters -- Input-output operators. *
 ***************************************************************/
ostream &operator << (ostream &out_file,
		      const Reconstruct1D_Input_Parameters &IP) {
  out_file << setprecision(6);
  out_file << "\n\n Solving 1D Reconstruction for the function \""
	   << IP.Function_Type << "\"";
  out_file << "\n  -> Method: "
	   << IP.Method_Used;

  switch (IP.Method){
  case DD_ENO:
    out_file << "\n  -> Solving the Linear System with: Householder "
	     << "transformation";
    out_file << "\n  -> Geometric weighting: "
	     << IP.Geometric_Weighting;
    out_file << "\n  -> Data-dependent weighting: "
	     << IP.Data_Dependent_Weighting;
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
	     << IP.FitTolerance();
    break;
  } /* endswitch */

  out_file << "\n  -> Reconstruction Order: "
	   << IP.Reconstruction_Order;
  out_file << "\n  -> Input File Name: " 
	   << IP.Input_File_Name;
  out_file << "\n  -> Grid: " 
	   << IP.Grid_Type;

  switch (IP.i_Grid){
  case GRID_READ_FROM_GRID_DATA_FILE:
    break;
  default:
    out_file << "\n  -> Starting Point of Solution Domain: " 
	     << IP.X_min;
    out_file << "\n  -> Size of Solution Domain: " 
	     << IP.LengthX();
    out_file << "\n  -> Characteristic Length of the Domain: " 
	     << IP.CharacteristicLength;
  } /* endswitch */

  out_file << "\n  -> Mesh shift and scale: " 
	   << IP.X_Shift << "  ,    " << IP.X_Scale;
  out_file << "\n  -> Number of Cells i-direction: "
	   << IP.Number_of_Cells_Idir;
  out_file << "\n  -> Number of SubGrid Points: "
	   << IP.Number_of_SubGrid_Points;
  out_file << "\n  -> Output File Name: " 
	   << IP.Output_File_Name;
  out_file << "\n  -> Output Format: " 
	   << IP.Output_Format_Type;
  out_file << "\n  -> Method for Computing Average Cell Solution: "
	   << IP.Integration_Type << " Integration";
    
  out_file << "\n";
  return (out_file);
}

istream &operator >> (istream &in_file,
		      Reconstruct1D_Input_Parameters &IP) {
  return (in_file);
}

/*************************************************************
 * _Reconstruct1D_Input_Parameters -- External subroutines. *
 *************************************************************/

/********************************************************
 * Routine: Open_Input_File                             *
 *                                                      *
 * Opens the appropriate input data file.               *
 *                                                      *
 ********************************************************/
void Open_Input_File(Reconstruct1D_Input_Parameters &IP) {

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
void Close_Input_File(Reconstruct1D_Input_Parameters &IP) {

    IP.Input_File.unsetf(ios::skipws);
    IP.Input_File.close();
}

/********************************************************
 * Routine: Set_Default_Input_Parameters                *
 *                                                      *
 * Assigns default values to the input parameters.      *
 *                                                      *
 ********************************************************/
void Set_Default_Input_Parameters(Reconstruct1D_Input_Parameters &IP) {

  strcpy(IP.Input_File_Name, "reconstruct.in");
  // Method settings
  strcpy(IP.Method_Used,"Normal_Equation");
  IP.Method = DD_ENO;
  strcpy(IP.Geometric_Weighting, "No");
  IP.geom_weighting = OFF;
  strcpy(IP.Data_Dependent_Weighting,"No");
  IP.data_depend_weighting = OFF;
  strcpy(IP.Use_Residual_DI,"No");
  IP.use_residual_DI = OFF;

  strcpy(IP.Limiter_Type, "Barth-Jespersen");
  IP.i_Limiter = LIMITER_BARTH_JESPERSEN;
  IP.CENO_Cutoff = 30;

  // Function settings
  strcpy(IP.Function_Type, "Function_Default");
  IP.i_Function = FUNCTION_DEFAULT;
  IP.TestF = Test_Default1D;
  IP.IntTestF = Test_Default1D_Integral;
  strcpy(IP.Integration_Type, "Theoretic");
  IP.Reconstruction_Order = 1;
  // Grid settings
  IP.X_min = -ONE;
  IP.X_max = ONE;
  IP.X_Shift = ZERO;
  IP.X_Scale = ONE;
  IP.CharacteristicLength = 1.0;

  strcpy(IP.Grid_Type, "Grid_Uniform" );
  IP.i_Grid = GRID_UNIFORM;
  strcpy(IP.CellNumber_or_DeltaCell, "N");
  IP.Number_of_Cells_Idir = 100;
  IP.Number_of_Ghost_Cells = 2;

  IP.Delta_Cell = IP.LengthX()/IP.Number_of_Cells_Idir;
  strcpy(IP.SubGridPoints_or_DeltaSubGrid, "N");
  IP.Number_of_SubGrid_Points = 5;
  IP.Delta_of_SubGrid = IP.Delta_Cell/IP.Number_of_SubGrid_Points;

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
void Get_Next_Input_Control_Parameter(Reconstruct1D_Input_Parameters &IP) {

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
int Parse_Next_Input_Control_Parameter(Reconstruct1D_Input_Parameters &IP) {

    int i_command;
    char buffer[256];

    i_command = 0;

    if (strcmp(IP.Next_Control_Parameter, "Function_Definition") == 0) {
       i_command = 1;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Function_Type, 
              IP.Next_Control_Parameter);

       // Add the parse word for the new function //
       if (strcmp(IP.Function_Type, "Default") == 0){
	 IP.i_Function = FUNCTION_DEFAULT;
	 IP.TestF = Test_Default1D;
	 IP.IntTestF = Test_Default1D_Integral;
       }
       else if (strcmp(IP.Function_Type, "Example1") == 0){
	 IP.i_Function = FUNCTION_EXAMPLE_1;
	 IP.TestF = Test_Example1;
	 IP.IntTestF = Test_Example1_Integral;
       }
       else if (strcmp(IP.Function_Type, "Example2") == 0){
	 IP.i_Function = FUNCTION_EXAMPLE_2;
	 IP.TestF = Test_Example2;
	 IP.IntTestF = Test_Example2_Integral;
       }
       else if (strcmp(IP.Function_Type, "Example3") == 0){
	 IP.i_Function = FUNCTION_EXAMPLE_3;
	 IP.TestF = Test_Example3;
	 IP.IntTestF = Test_Example3_Integral;
       }
       else if (strcmp(IP.Function_Type, "Example4") == 0){
	 IP.i_Function = FUNCTION_EXAMPLE_4;
	 IP.TestF = Test_Example4;
	 IP.IntTestF = Test_Example4_Integral;
       }
       else if (strcmp(IP.Function_Type, "Example5") == 0){
	 IP.i_Function = FUNCTION_EXAMPLE_5;
	 IP.TestF = Test_Example5;
	 IP.IntTestF = Test_Example5_Integral;
       }
       else if (strcmp(IP.Function_Type, "Example6") == 0){
	 IP.i_Function = FUNCTION_EXAMPLE_6;
	 IP.TestF = Test_Example6;
	 IP.IntTestF = Test_Example6_Integral;
       }
       else if (strcmp(IP.Function_Type, "Example7") == 0){
	 IP.i_Function = FUNCTION_EXAMPLE_7;
	 IP.TestF = Test_Example7;
	 IP.IntTestF = Test_Example7_Integral;
       }
       else if (strcmp(IP.Function_Type, "Example8") == 0){
	 IP.i_Function = FUNCTION_EXAMPLE_8;
	 IP.TestF = Test_Example8;
	 IP.IntTestF = Test_Example8_Integral;
       }

       else if (strcmp(IP.Function_Type, "Example9") == 0){
	 IP.i_Function = FUNCTION_EXAMPLE_9;
	 IP.TestF = Test_Example9;
	 IP.IntTestF = Test_Example9_Integral;
       }

       else if (strcmp(IP.Function_Type, "Example10") == 0){
	 IP.i_Function = FUNCTION_EXAMPLE_10;
	 IP.TestF = Test_Example10;
	 IP.IntTestF = Test_Example10_Integral;
       }

       // endif 

    } else if (strcmp(IP.Next_Control_Parameter, "Reconstruction_Order") == 0){
	i_command = 2;
	IP.Line_Number ++;

	if (!(IP.Input_File >> IP.Reconstruction_Order))
	  i_command = INVALID_INPUT_VALUE;
	IP.Input_File.getline(buffer, sizeof(buffer));
	if (IP.Reconstruction_Order < ZERO) {
	    IP.Reconstruction_Order = 0;
	    cout << "\n MESSAGE "<< IP.Message_Number
		 << ": Reconstruction Order must be at least 1.\a";
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
	IP.Number_of_Ghost_Cells = 6;

    } else if (strcmp(IP.Next_Control_Parameter, "Grid_Type") == 0) {
       i_command = 3;
       Get_Next_Input_Control_Parameter(IP);
       strcpy(IP.Grid_Type, 
              IP.Next_Control_Parameter);

       if (strcmp(IP.Grid_Type, "Grid_Uniform") == 0) {
	 IP.i_Grid = GRID_UNIFORM;
       } else if (strcmp(IP.Grid_Type, "Grid_NonUniform") == 0) {
	 IP.i_Grid = GRID_NONUNIFORM; //!!!!! It should be specified the definition file
       } else if (strcmp(IP.Grid_Type, "Grid_Predefined1") == 0) {
	 IP.i_Grid = GRID_PREDEFINED1;
	 IP.X_min = -1.5;
	 IP.X_max = 1.5;
       } else if (strcmp(IP.Grid_Type, "Grid_Predefined2") == 0) {
	 IP.i_Grid = GRID_PREDEFINED2;
	 IP.X_min = -10.0;
	 IP.X_max = 10.0;
       } else if (strcmp(IP.Grid_Type, "Read_From_Definition_File") == 0) {
	 IP.i_Grid = GRID_READ_FROM_DEFINITION_FILE;
       } else if (strcmp(IP.Grid_Type, "Read_From_Data_File") == 0) {
	 IP.i_Grid = GRID_READ_FROM_GRID_DATA_FILE;
       } else {
	 IP.i_Grid = GRID_UNIFORM ;
       } // endif 
       
    } else if (strcmp(IP.Next_Control_Parameter, "Output_File_Name") == 0) {
      i_command = 4;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Output_File_Name, 
	     IP.Next_Control_Parameter);
      strcpy(IP.Grid_File_Name, 
	     IP.Next_Control_Parameter);
      strcpy(IP.Grid_Definition_File_Name, 
	     IP.Next_Control_Parameter);
      strcpy(IP.Gnuplot_File_Name, 
	     IP.Next_Control_Parameter);
      
    } else if (strcmp(IP.Next_Control_Parameter, "Grid_File_Name") == 0) {
      i_command = 5;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Grid_File_Name, 
	     IP.Next_Control_Parameter);
      strcat(IP.Grid_File_Name, ".grid");
      strcpy(IP.Grid_Definition_File_Name, 
	     IP.Next_Control_Parameter);
      strcat(IP.Grid_Definition_File_Name, ".griddef");
      
    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Cells_Idir") == 0) {
      i_command = 6;
      IP.Line_Number ++;
      if ((IP.Input_File >> IP.Number_of_Cells_Idir) == 0 || 
	  (IP.Number_of_Cells_Idir < ONE)){
	i_command = INVALID_INPUT_VALUE;
      }
      IP.Input_File.getline(buffer, sizeof(buffer));
      
    } else if (strcmp(IP.Next_Control_Parameter, "Number_of_Subgrid_Points") == 0){
      i_command = 7;
      IP.Line_Number ++;
      if ((IP.Input_File >> IP.Number_of_SubGrid_Points) == 0 ||
	  (IP.Number_of_SubGrid_Points < 1))
	i_command = INVALID_INPUT_VALUE;

      IP.Input_File.getline(buffer, sizeof(buffer));
      
    } else if (strcmp(IP.Next_Control_Parameter, "X_min") == 0) {
      i_command = 8;
      IP.Line_Number ++;
      if ((IP.i_Grid != GRID_PREDEFINED1) && (IP.i_Grid != GRID_PREDEFINED2)){
	if (!(IP.Input_File >> IP.X_min))
	  i_command = INVALID_INPUT_VALUE;
      } //endif
      IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Characteristic_Length") == 0) {
      i_command = 9;
      IP.Line_Number ++;
      if (!(IP.Input_File >> IP.CharacteristicLength))
	i_command = INVALID_INPUT_VALUE;
      IP.Input_File.getline(buffer, sizeof(buffer));
      
    } else if (strcmp(IP.Next_Control_Parameter, "X_max") == 0) {
      i_command = 10;
      IP.Line_Number ++;
      if ((IP.i_Grid != GRID_PREDEFINED1) && (IP.i_Grid != GRID_PREDEFINED2)){
	if (!(IP.Input_File >> IP.X_max) || (IP.MinX() >= IP.MaxX() ))
	  i_command = INVALID_INPUT_VALUE;
      } // endif
      IP.Input_File.getline(buffer, sizeof(buffer));
      
    } else if (strcmp(IP.Next_Control_Parameter, "Output_Format_Type") == 0) {
      i_command = 11;
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
      i_command = 12;
      IP.Line_Number ++;
      if(!(IP.Input_File >> IP.X_Shift))
	i_command = INVALID_INPUT_VALUE;
      IP.Input_File.setf(ios::skipws);
      IP.Input_File.getline(buffer, sizeof(buffer));
      
    } else if (strcmp(IP.Next_Control_Parameter, "X_Scale") == 0) {
      i_command = 13;
      IP.Line_Number = IP.Line_Number + 1;
      if (!(IP.Input_File >> IP.X_Scale))
	i_command = INVALID_INPUT_VALUE;
      IP.Input_File.getline(buffer, sizeof(buffer));

    } else if (strcmp(IP.Next_Control_Parameter, "Method_Used") == 0) {
      i_command = 14;
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
      i_command = 15;
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
      i_command = 16;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Data_Dependent_Weighting, 
	     IP.Next_Control_Parameter);

       // Add the parse word for the new function //
       if (strcmp(IP.Data_Dependent_Weighting, "No") == 0)
	 IP.data_depend_weighting = OFF;
       else if (strcmp(IP.Data_Dependent_Weighting, "Yes") == 0)
	 IP.data_depend_weighting = ON;

    } else if (strcmp(IP.Next_Control_Parameter, 
		      "Use_Residual_DI") == 0){
      i_command = 17;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Use_Residual_DI, 
	     IP.Next_Control_Parameter);

       // Add the parse word for the new function //
       if (strcmp(IP.Use_Residual_DI, "No") == 0)
	 IP.use_residual_DI = OFF;
       else if (strcmp(IP.Use_Residual_DI, "Yes") == 0)
	 IP.use_residual_DI = ON;

    } else if (strcmp(IP.Next_Control_Parameter, "Integration_Type") == 0){
      i_command = 18;
      Get_Next_Input_Control_Parameter(IP);
      strcpy(IP.Integration_Type,
	     IP.Next_Control_Parameter);

    } else if (strcmp(IP.Next_Control_Parameter, "Limiter") == 0){
      i_command = 19;
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

    } else if (strcmp(IP.Next_Control_Parameter, "CENO_Tolerance") == 0){
      i_command = 20;
      IP.Line_Number ++;
      if(!(IP.Input_File >> IP.FitTolerance()))
	i_command = INVALID_INPUT_VALUE;
      IP.Input_File.setf(ios::skipws);
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
      
    } else if (strcmp(IP.Next_Control_Parameter, "Write_Stencil_Reconstruction") == 0) {
      i_command = WRITE_RECONSTRUCTED_FUNCTIONS;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Cells") == 0) {
      i_command = WRITE_OUTPUT_CELLS_CODE;
      
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

    } else if (strcmp(IP.Next_Control_Parameter, "Write_L1_Norm") == 0) {
      i_command = WRITE_L1_NORM;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_L2_Norm") == 0) {
      i_command = WRITE_L2_NORM;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Norm_On_Screen") == 0) {
      i_command = WRITE_NORM_ON_SCREEN;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Mesh_Nodes") == 0) {
      i_command = WRITE_OUTPUT_GRID_NODES_CODE;
      
    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Mesh_Cells") == 0) {
      i_command = WRITE_OUTPUT_GRID_CELLS_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Derivatives") == 0) {
      i_command = WRITE_OUTPUT_DERIVATIVES_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Cell_Comparison") == 0) {
      i_command = WRITE_OUTPUT_COMPARISON_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Function_Graph") == 0) {
      i_command = WRITE_OUTPUT_FUNCTION_GRAPH_CODE;

    } else if (strcmp(IP.Next_Control_Parameter, "Write_Output_Accuracy") == 0) {
      i_command = WRITE_OUTPUT_ACCURACY_CODE;
      
    } else if (strcmp(IP.Next_Control_Parameter, "Print_Norms") == 0) {
      i_command = WRITE_NORM_ON_SCREEN;

    } else if (strcmp(IP.Next_Control_Parameter, "Assess_Accuracy") == 0) {
      i_command = ASSESS_ACCURACY_CODE;
      
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
int Process_Input_Control_Parameter_File(Reconstruct1D_Input_Parameters &Input_Parameters,
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
       cout << "\n Reconstruction1D ERROR: Unable to open Reconstruction1D input data file.\n";
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
          cout << "\n Reconstruction1D ERROR: Error reading Reconstruction1D data at line #"
               << -line_number  << " of input data file.\n\a";
          error_flag = line_number;
          break;
       } /* endif */
    } /* endwhile */

    /* Initial processing of input control parameters complete.  
       Return the error indicator flag. */

    return (error_flag);

}
