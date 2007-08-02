/* Heat1DSolvers.cc:  1D Heat Equation Solvers. */

/* Include 1D Heat equation solution header file. */

#ifndef _HEAT1D_INCLUDED
#include "Heat1D.h"
#endif // _HEAT1D_INCLUDED

/********************************************************
 * Routine: Heat1DSolver                                *
 *                                                      *
 * Computes solutions to 1D Heat equation.              *
 *                                                      *
 ********************************************************/
int Heat1DSolver(char *Input_File_Name_ptr,
	         int batch_flag) {

  /********************************************************  
   * Local variable declarations.                         *
   ********************************************************/

  // 1D CFD input variables and parameters:
  CFD1D_Input_Parameters Input_Parameters;

  // Output file name pointer:
  char *Output_File_Name_ptr;

  // Gnuplot file name pointer:
  char *Gnuplot_File_Name_ptr;

  // Output file stream:
  ofstream Output_File;

  // Gnuplot macro file stream:
  ofstream Gnuplot_File;

  // Exact Solution Output file name:
  char *Exact_Output_File_Name_ptr;
  char Exact_Output_File_Name[256];

  /* Solution variables. */

  Heat1D_UniformMesh *Soln_ptr, *Exact_Soln_ptr;

  /* Variables for tridiagonal system of linear equations. */

  TriDiagonalSystemLinEqs Ax_equal_b;

  /* Other local solution variables. */

  int number_of_time_steps,
      command_flag, error_flag, line_number, i;

  double time, dtime;
  
  /********************************************************  
   * Set default values for the input solution parameters *
   * and then read user specified input values from the   *
   * specified input file.                                *
   ********************************************************/

  Set_Default_Input_Parameters(Input_Parameters);

  strcpy(Input_Parameters.Input_File_Name, Input_File_Name_ptr);
  Open_Input_File(Input_Parameters);
  if (Input_Parameters.Input_File.bad()) {
     if (batch_flag) {
        cout << "\nPDES++ ERROR: Unable to open CFD1D input data file.\n\n";
     } else {
        cout << "\n PDES++ ERROR: Unable to open CFD1D input data file.";
        cout << "\n\nPDES++: Execution terminated.\n";
     } /* endif */
     return (-1);
  } /* endif */

  if (! batch_flag) {
    cout << "\n Reading CFD1D input data file `"
	 << Input_Parameters.Input_File_Name << "'.";
  }
  while (1) {
     Get_Next_Input_Control_Parameter(Input_Parameters);
     command_flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
     line_number = Input_Parameters.Line_Number;
     if (command_flag == EXECUTE_CODE) {
         break;

     } else if (command_flag == TERMINATE_CODE) {
         return (0);

     } else if (command_flag == INVALID_INPUT_CODE ||
                command_flag == INVALID_INPUT_VALUE) {
         line_number = -line_number;
         if (batch_flag) {
             cout << "\nPDES++ ERROR: Error reading CFD1D data at line #"
	          << -line_number << " of input data file.\n\n";
         } else {
             cout << "\n PDES++ ERROR: Error reading CFD1D data at line #"
	          << -line_number  << " of input data file.";
             cout << "\n\nPDES++: Execution terminated.\n";
         } /* end if */
         return (line_number);
     } /* endif */
  } /* endwhile */

  if (! batch_flag)  {
     cout << Input_Parameters << "\n";
  } /* endif */

  execute_new_calculation: ;

  i = 0;
  while (Input_Parameters.Output_File_Name[i] != '.' ) {
     Exact_Output_File_Name[i] = Input_Parameters.Output_File_Name[i];
     i = i + 1;
     if (i > strlen(Input_Parameters.Output_File_Name) ) break;
  } /* endwhile */
  Exact_Output_File_Name[i] = '\0';
  strcat(Exact_Output_File_Name, ".exact");

  /*********************************************************  
   * Create mesh and allocate Heat1D solution variables.   *
   *********************************************************/

  /* Allocate memory for 1D heat equation solution on
     uniform mesh. */

  if (! batch_flag) cout << "\n Creating memory for Heat1D solution variables.";
  Soln_ptr=Allocate(Soln_ptr,
                    Input_Parameters.Number_of_Nodes);
  Exact_Soln_ptr=Allocate(Exact_Soln_ptr,
                          Input_Parameters.Number_of_Nodes);

  /* Allocate memory for tridiagonal system of linear equations
     to be solved in simple implicit and Crank-Nicolson methods. */

  Ax_equal_b.allocate(Input_Parameters.Number_of_Nodes);

  /* Create uniform mesh. */

  if (! batch_flag) cout << "\n Creating uniform mesh.";
  Grid(Soln_ptr, 
       Input_Parameters.X_Min, 
       Input_Parameters.X_Max, 
       Input_Parameters.Number_of_Nodes);
  Grid(Exact_Soln_ptr, 
       Input_Parameters.X_Min, 
       Input_Parameters.X_Max,
       Input_Parameters.Number_of_Nodes);

  /********************************************************  
   * Initialize Heat1D solution variables.               *
   ********************************************************/

  /* Set the initial time level. */

  time = ZERO;
  number_of_time_steps = 0;
  
  /* Initialize solution variables. */

  if (! batch_flag) cout << "\n Prescribing Heat1D initial data.";
  ICs(Soln_ptr, 
      Input_Parameters.i_ICs,
      Input_Parameters.Kappa, 
      Input_Parameters.Number_of_Nodes);

  /* Determine the heat flux. */

  HeatFlux(Soln_ptr, 
           Input_Parameters.Number_of_Nodes);

  /* Determine the exact analytical solution at initial time level. */

  if (! batch_flag) cout << "\n Computing exact initial solution.";
  AnalyticSoln(Exact_Soln_ptr, 
               Input_Parameters.i_ICs,
               Input_Parameters.Kappa, 
               time, 
               Input_Parameters.Number_of_Nodes);
 
  /********************************************************  
   * Solve IBVP or BVP for heat equation on a             *
   * uniform (constant cell size) 1D mesh.                *
   ********************************************************/

  continue_existing_calculation: ;

  if ((Input_Parameters.Local_Time_Stepping &&
       Input_Parameters.Maximum_Number_of_Time_Steps > 0) ||
      (!Input_Parameters.Local_Time_Stepping &&
       Input_Parameters.Time_Max > time)) {
     if (! batch_flag) cout << "\n\n Beginning Heat1D computations.\n\n";
     while (1) {
         /* Determine local and global time steps. */
         dtime = Rstability(Soln_ptr, 
                            Input_Parameters.Number_of_Nodes);
         if (time + Input_Parameters.CFL_Number*dtime > 
             Input_Parameters.Time_Max) {
            dtime = (Input_Parameters.Time_Max-time)/Input_Parameters.CFL_Number;
         } /* endif */

         /* Output progress information for the calculation. */
         if (! batch_flag && number_of_time_steps == 0 ) {
             cout << " Time Step = " << number_of_time_steps
	          << " Time = " << time << "\n .";
         } else if (! batch_flag && 
                    number_of_time_steps-100*(number_of_time_steps/100) == 0 ) {
	     cout << "\n" << " Time Step = " << number_of_time_steps
	          << " Time = " << time << "\n .";
         } else if (! batch_flag && 
                    number_of_time_steps-50*(number_of_time_steps/50) == 0 ) {
	     cout << "\n .";
         } else if (! batch_flag) {
	     cout << ".";
         } /* endif */

         /* Check to see if calculations are complete. */
         if (Input_Parameters.Local_Time_Stepping &&
             number_of_time_steps >= 
             Input_Parameters.Maximum_Number_of_Time_Steps) break;
         if (!Input_Parameters.Local_Time_Stepping &&
             time >= Input_Parameters.Time_Max) break;

         /* Update solution for next time step. */
         switch (Input_Parameters.i_Time_Integration) {
           case TIME_STEPPING_SIMPLE_EXPLICIT:
             error_flag = dTdt_SimpleExplicit(Soln_ptr,
                           Input_Parameters.Number_of_Nodes,
          		   dtime,
                           Input_Parameters.CFL_Number,
       			   Input_Parameters.Local_Time_Stepping);
             break;
           case TIME_STEPPING_SIMPLE_IMPLICIT:
             error_flag = dTdt_SimpleImplicit(Soln_ptr,
                           Input_Parameters.Number_of_Nodes,
        		   dtime,
                           Input_Parameters.CFL_Number,
       			   Input_Parameters.Local_Time_Stepping,
                           Ax_equal_b);
             break;
           case TIME_STEPPING_CRANK_NICOLSON:
             error_flag = dTdt_CrankNicolson(Soln_ptr,
                           Input_Parameters.Number_of_Nodes,
        		   dtime,
                           Input_Parameters.CFL_Number,
       			   Input_Parameters.Local_Time_Stepping,
                           Ax_equal_b);
             break;
           case TIME_STEPPING_ADE:
             error_flag = dTdt_ADE(Soln_ptr,
                           Input_Parameters.Number_of_Nodes,
        		   dtime,
                           Input_Parameters.CFL_Number,
       			   Input_Parameters.Local_Time_Stepping);
             break;
           default:
             error_flag = dTdt_SimpleExplicit(Soln_ptr,
                           Input_Parameters.Number_of_Nodes,
          		   dtime,
                           Input_Parameters.CFL_Number,
       			   Input_Parameters.Local_Time_Stepping);
             break;
         } /* endswitch */

         if (error_flag) {
             if (batch_flag) {
                 cout << "\nPDES++ ERROR: Heat1D solution error.\n\n";
             } else {
                 cout << "\n PDES++ ERROR: Heat1D solution error.";
                 cout << "\n\nPDES++: Execution terminated.\n";
             } /* endif */
             return (error_flag);
         } /* endif */

         /* Compute the heat flux. */
         HeatFlux(Soln_ptr, 
                  Input_Parameters.Number_of_Nodes);
      
         /* Update time and time step counter. */
         number_of_time_steps = number_of_time_steps + 1;
         time = time + Input_Parameters.CFL_Number*dtime;

     } /* endwhile */
     if (! batch_flag) cout << "\n\n Heat1D computations complete.\n";
  } /* endif */

  /* Determine the exact analytical solution at final time level. */

  if (! batch_flag) cout << "\n Computing exact final solution.";
  AnalyticSoln(Exact_Soln_ptr, 
               Input_Parameters.i_ICs,
               Input_Parameters.Kappa, 
               time, 
               Input_Parameters.Number_of_Nodes);

  /********************************************************  
   * Solution calculations complete.                      *
   * Write 1D Heat equation solution to output file as    *
   * required, reset solution parameters, and run other   *
   * cases as specified by input data.                    *
   ********************************************************/

  while (1) {
     Get_Next_Input_Control_Parameter(Input_Parameters);
     command_flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
     line_number = Input_Parameters.Line_Number;
     if (command_flag == EXECUTE_CODE) {
         // Deallocate memory for 1D heat equation solution.
         if (! batch_flag) cout << "\n Deallocating Heat1D solution variables.";
         Ax_equal_b.deallocate();
         Soln_ptr=Deallocate(Soln_ptr,
                             Input_Parameters.Number_of_Nodes);
         Exact_Soln_ptr=Deallocate(Exact_Soln_ptr,
                                   Input_Parameters.Number_of_Nodes);
         // Execute new calculation.
         if (! batch_flag)  {
            cout << "\n\n Starting a new calculation.";
            cout << Input_Parameters << "\n";
         } /* endif */
         goto execute_new_calculation;

     } else if (command_flag == TERMINATE_CODE) {
         // Deallocate memory for 1D heat equation solution.
         if (! batch_flag) cout << "\n Deallocating Heat1D solution variables.";
         Ax_equal_b.deallocate();
         Soln_ptr=Deallocate(Soln_ptr,
                             Input_Parameters.Number_of_Nodes);
         Exact_Soln_ptr=Deallocate(Exact_Soln_ptr,
                                   Input_Parameters.Number_of_Nodes);
         // Close input data file.
         if (! batch_flag) cout << "\n\n Closing CFD1D input data file.";
         Close_Input_File(Input_Parameters);
         // Terminate calculation.
         return (0);

     } else if (command_flag == CONTINUE_CODE) {
         // Reset maximum time step counter.
         Input_Parameters.Maximum_Number_of_Time_Steps += number_of_time_steps;
         // Continue existing calculation.
         if (! batch_flag)  {
            cout << "\n\n Continuing existing calculation.";
            cout << Input_Parameters << "\n";
         } /* endif */
         goto continue_existing_calculation;

     } else if (command_flag == WRITE_OUTPUT_CODE) {
         // Output solution data.
         Output_File_Name_ptr = Input_Parameters.Output_File_Name;
         if (! batch_flag)
           cout << "\n Writing Heat1D solution to output data file `"
	        << Output_File_Name_ptr << "'.";
         Output_File.open(Output_File_Name_ptr, ios::out);
         if (Output_File.bad()) {
            if (batch_flag) {
               cout << "\nPDES++ ERROR: Unable to open Heat1D output data file.\n\n";
            } else {
               cout << "\n PDES++ ERROR: Unable to open Heat1D output data file.";
               cout << "\n\nPDES++: Execution terminated.\n";
            } /* endif */
            return (1);
         } /* endif */

         if (Input_Parameters.i_Output_Format == IO_GNUPLOT) {
            Output_Gnuplot(Soln_ptr,
	  	           Input_Parameters.Number_of_Nodes,
		           number_of_time_steps,
		           time,
		           Output_File);
         } else if (Input_Parameters.i_Output_Format == IO_TECPLOT) {
            Output_Tecplot(Soln_ptr,
	  	           Input_Parameters.Number_of_Nodes,
		           number_of_time_steps,
		           time,
		           Output_File);
         } /* endif */

         Output_File.close();

         Exact_Output_File_Name_ptr = Exact_Output_File_Name;
         if (! batch_flag)
           cout << "\n Writing exact analytical Heat1D solution to output data file `"
	        << Exact_Output_File_Name_ptr << "'.";
         Output_File.open(Exact_Output_File_Name_ptr, ios::out);
         if (Output_File.bad()) {
            if (batch_flag) {
               cout << "\nPDES++ ERROR: Unable to open exact Heat1D output data file.\n\n";
            } else {
               cout << "\n PDES++ ERROR: Unable to open exact Heat1D output data file.";
               cout << "\n\nPDES++: Execution terminated.\n";
            } /* endif */
            return (2);
         } /* endif */

         if (Input_Parameters.i_Output_Format == IO_GNUPLOT) {
            Output_Gnuplot(Exact_Soln_ptr,
	         	    Input_Parameters.Number_of_Nodes,
		           number_of_time_steps,
		           time,
		           Output_File);
         } else if (Input_Parameters.i_Output_Format == IO_TECPLOT) {
            Output_Tecplot(Exact_Soln_ptr,
	  	           Input_Parameters.Number_of_Nodes,
		           number_of_time_steps,
		           time,
		           Output_File);
         } /* endif */

         Output_File.close();

         if (Input_Parameters.i_Output_Format == IO_GNUPLOT) {
            Gnuplot_File_Name_ptr = Input_Parameters.Gnuplot_File_Name;
            Gnuplot_File.open(Gnuplot_File_Name_ptr, ios::out);
            if (Gnuplot_File.bad()) {
               if (batch_flag) {
                  cout << "\nPDES++ ERROR: Unable to open Heat1D Gnuplot macro file.\n\n";
               } else {
                  cout << "\n PDES++ ERROR: Unable to open Heat1D Gnuplot macro file.";
                  cout << "\n\nPDES++: Execution terminated.\n";
               } /* endif */
               return (3);
            } /* endif */

            Gnuplot_File << "set title \"PDES++: 1D Heat Equation Solution\"\n"
                         << "set xlabel \"x\"\n"
	                 << "plot \"" << Output_File_Name_ptr << "\""
                         << " using 1:2 \"%*lf%lf%*lf%lf%*lf%*lf\" \\\n"
	                 << "     title \"T\" with points, \\\n"
                         << "\"" << Exact_Output_File_Name_ptr << "\""
                         << " using 1:2 \"%*lf%lf%*lf%lf%*lf%*lf\" \\\n"
	                 << "     title \"Te\" with linespoints\n"
                         << "pause -1  \"Hit return to continue\"\n"
	                 << "plot \"" << Output_File_Name_ptr << "\""
                         << " using 1:2 \"%*lf%lf%*lf%*lf%lf%*lf\" \\\n"
	                 << "     title \"qx\" with points, \\\n"
                         << "\"" << Exact_Output_File_Name_ptr << "\""
                         << " using 1:2 \"%*lf%lf%*lf%*lf%lf%*lf\" \\\n"
	                 << "     title \"qxe\" with linespoints\n"
                         << "pause -1  \"Hit return to continue\"\n";

            Gnuplot_File.close();
         } /* endif */

     } else if (command_flag == INVALID_INPUT_CODE ||
                command_flag == INVALID_INPUT_VALUE) {
         line_number = -line_number;
         if (batch_flag) {
             cout << "\nPDES++ ERROR: Error reading CFD1D data at line #"
	          << -line_number << " of input data file.\n\n";
         } else {
             cout << "\n PDES++ ERROR: Error reading CFD1D data at line #"
	          << -line_number  << " of input data file.";
             cout << "\n\nPDES++: Execution terminated.\n";
         } /* end if */
         return (line_number);
     } /* endif */
  } /* endwhile */

  /********************************************************  
   * End of all Heat1DSolver computations and I/O.        *
   ********************************************************/

  return (0);
  
}
