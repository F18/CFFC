/* MHD1DSolvers.cc:  1D MHD Equation Solvers. */

/* Include 1D MHD solution header file. */

#ifndef _MHD1D_INCLUDED
#include "MHD1D.h"
#endif // _MHD1D_INCLUDED

/********************************************************
 * Routine: MHD1DSolver                                 *
 *                                                      *
 * Computes solutions to 1D MHD equations.              *
 *                                                      *
 ********************************************************/
int MHD1DSolver(char *Input_File_Name_ptr,
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

  /* Solution variables. */

  MHD1D_UniformMesh *Soln_ptr;

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

  /*********************************************************  
   * Create mesh and allocate MHD1D solution variables.    *
   *********************************************************/

  execute_new_calculation: ;

  /* Allocate memory for 1D MHD equation solution on
     uniform mesh. */

  if (! batch_flag) cout << "\n Creating memory for MHD1D solution variables.";
  Soln_ptr=Allocate(Soln_ptr,
                    Input_Parameters.Number_of_Cells);

  /* Create uniform mesh. */

  if (! batch_flag) cout << "\n Creating uniform mesh.";
  Grid(Soln_ptr, 
       Input_Parameters.X_Min, 
       Input_Parameters.X_Max, 
       Input_Parameters.Number_of_Cells);

  /* Set the initial time level. */

  /********************************************************  
   * Initialize MHD1D solution variables.                 *
   ********************************************************/

  /* Set the initial time level. */

  time = ZERO;
  number_of_time_steps = 0;
  
  /* Initialize the conserved and primitive state
     solution variables. */
  
  if (! batch_flag) cout << "\n Prescribing MHD1D initial data.";
  if (Input_Parameters.i_ICs == IC_BRIO_WU) {
      ICs(Soln_ptr, 
          "BRIOWU", 
          Input_Parameters.i_ICs, 
          Input_Parameters.Number_of_Cells);
  } else if (Input_Parameters.i_ICs == IC_DAI_WOODWARD) {
      ICs(Soln_ptr, 
          "MONATOMIC", 
          Input_Parameters.i_ICs, 
          Input_Parameters.Number_of_Cells);
  } else {
      ICs(Soln_ptr, 
          "MONATOMIC", 
          Input_Parameters.i_ICs, 
          Input_Parameters.Number_of_Cells);
  } /* endif */
  
  /********************************************************  
   * Solve conservation form of 1D MHD equations for      *
   * specified IBVP or BVP on uniform mesh.               *
   ********************************************************/

  continue_existing_calculation: ;

  if ((Input_Parameters.Local_Time_Stepping &&
       Input_Parameters.Maximum_Number_of_Time_Steps > 0) ||
      (!Input_Parameters.Local_Time_Stepping &&
       Input_Parameters.Time_Max > time)) {
     if (! batch_flag) cout << "\n\n Beginning MHD1D computations.\n\n";
     while (1) {
         /* Determine local and global time steps. */
         dtime = CFL(Soln_ptr, 
                     Input_Parameters.Number_of_Cells);
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
           case TIME_STEPPING_EXPLICIT_EULER:
             error_flag = dUdt_explicitEuler_upwind(Soln_ptr,
                           Input_Parameters.Number_of_Cells,
        		   dtime,
                           Input_Parameters.CFL_Number,
                           Input_Parameters.i_Flux_Function,
       			   Input_Parameters.Local_Time_Stepping);
             break;
           case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR:

             break;
           case TIME_STEPPING_LAX_FRIEDRICHS:
             error_flag = dUdt_LaxFriedrichs(Soln_ptr,
                           Input_Parameters.Number_of_Cells,
        		   dtime,
                           Input_Parameters.CFL_Number,
       			   Input_Parameters.Local_Time_Stepping);
             break;
           case TIME_STEPPING_LAX_WENDROFF:

             break;
           case TIME_STEPPING_MACCORMACK:

             break;
           case TIME_STEPPING_HANCOCK:

             break;
           default:
             error_flag = dUdt_explicitEuler_upwind(Soln_ptr,
                           Input_Parameters.Number_of_Cells,
        		   dtime,
                           Input_Parameters.CFL_Number,
                           Input_Parameters.i_Flux_Function,
       			   Input_Parameters.Local_Time_Stepping);
             break;
         } /* endswitch */

         if (error_flag) {
             if (batch_flag) {
                 cout << "\nPDES++ ERROR: MHD1D solution error.\n\n";
             } else {
                 cout << "\n PDES++ ERROR: MHD1D solution error.";
                 cout << "\n\nPDES++: Execution terminated.\n";
             } /* endif */
             return (error_flag);
         } /* endif */
      
         /* Update time and time step counter. */
         number_of_time_steps = number_of_time_steps + 1;
         time = time + Input_Parameters.CFL_Number*dtime;

     } /* endwhile */
     if (! batch_flag) cout << "\n\n MHD1D computations complete.\n";
  } /* endif */

  /********************************************************  
   * Solution calculations complete.                      *
   * Write 1D MHD solution to output file as              *
   * required, reset solution parameters, and run other   *
   * cases as specified by input data.                    *
   ********************************************************/

  while (1) {
     Get_Next_Input_Control_Parameter(Input_Parameters);
     command_flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
     line_number = Input_Parameters.Line_Number;
     if (command_flag == EXECUTE_CODE) {
         // Deallocate memory for 1D MHD equation solution.
         if (! batch_flag) cout << "\n Deallocating MHD1D solution variables.";
         Soln_ptr=Deallocate(Soln_ptr,
                             Input_Parameters.Number_of_Cells);
         // Execute new calculation.
         if (! batch_flag)  {
            cout << "\n\n Starting a new calculation.";
            cout << Input_Parameters << "\n";
         } /* endif */
         goto execute_new_calculation;

     } else if (command_flag == TERMINATE_CODE) {
         // Deallocate memory for 1D MHD equation solution.
         if (! batch_flag) cout << "\n Deallocating MHD1D solution variables.";
         Soln_ptr=Deallocate(Soln_ptr,
                             Input_Parameters.Number_of_Cells);
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
           cout << "\n Writing MHD1D solution to output data file `"
	        << Output_File_Name_ptr << "'.";
         Output_File.open(Output_File_Name_ptr, ios::out);
         if (Output_File.bad()) {
            if (batch_flag) {
               cout << "\nPDES++ ERROR: Unable to open MHD1D output data file.\n\n";
            } else {
               cout << "\n PDES++ ERROR: Unable to open MHD1D output data file.";
               cout << "\n\nPDES++: Execution terminated.\n";
            } /* endif */
            return (1);
         } /* endif */

         if (Input_Parameters.i_Output_Format == IO_GNUPLOT) {
            Output_Gnuplot(Soln_ptr,
	   	           Input_Parameters.Number_of_Cells,
		           number_of_time_steps,
		           time,
		           Output_File);
         } else if (Input_Parameters.i_Output_Format == IO_TECPLOT) {
            Output_Tecplot(Soln_ptr,
	   	           Input_Parameters.Number_of_Cells,
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
                  cout << "\nPDES++ ERROR: Unable to open MHD1D Gnuplot macro file.\n\n";
               } else {
                  cout << "\n PDES++ ERROR: Unable to open MHD1D Gnuplot macro file.";
                  cout << "\n\nPDES++: Execution terminated.\n";
               } /* endif */
               return (2);
            } /* endif */

            Gnuplot_File << "set title \"PDES++: 1D MHD Solution\"\n"
                         << "set xlabel \"x\"\n"
	                 << "plot \"" << Output_File_Name_ptr << "\""
                         << " using 1:2 \"%*lf%lf%*lf%lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf\" \\\n"
	                 << "     title \"density\" with points\n"
                         << "pause -1  \"Hit return to continue\"\n"
                         << "plot \"" << Output_File_Name_ptr << "\""
                         << " using 1:2 \"%*lf%lf%*lf%*lf%lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf\" \\\n"
	                 << "     title \"ux\" with points\n"
                         << "pause -1  \"Hit return to continue\"\n"
                         << "plot \"" << Output_File_Name_ptr << "\""
                         << " using 1:2 \"%*lf%lf%*lf%*lf%*lf%lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf\" \\\n"
	                 << "     title \"uy\" with points\n"
                         << "pause -1  \"Hit return to continue\"\n"
                         << "plot \"" << Output_File_Name_ptr << "\""
                         << " using 1:2 \"%*lf%lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%lf%*lf\" \\\n"
	                 << "     title \"pressure\" with points\n"
                         << "pause -1  \"Hit return to continue\"\n"
                         << "plot \"" << Output_File_Name_ptr << "\""
                         << " using 1:2 \"%*lf%lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%lf\" \\\n"
	                 << "     title \"temperature\" with points\n"
                         << "pause -1  \"Hit return to continue\"\n"
                         << "plot \"" << Output_File_Name_ptr << "\""
                         << " using 1:2 \"%*lf%lf%*lf%*lf%*lf%*lf%*lf%*lf%lf%*lf%*lf%*lf%*lf%*lf%*lf\" \\\n"
	                 << "     title \"By\" with points\n"
                         << "pause -1  \"Hit return to continue\"\n"
                         << "plot \"" << Output_File_Name_ptr << "\""
                         << " using 1:2 \"%*lf%lf%*lf%*lf%*lf%*lf%*lf%*lf%*lf%lf%*lf%*lf%*lf%*lf%*lf\" \\\n"
	                 << "     title \"Bz\" with points\n"
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
   * End of all MHD1DSolver computations and I/O.         *
   ********************************************************/

  return (0);
  
}
