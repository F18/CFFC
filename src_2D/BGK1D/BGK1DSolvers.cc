/*!\file BGK1DSolvers.cc
  \brief 1D BGK Equation Solvers. */

/* Include 1D BGK solution header file. */

#ifndef _BGK1D_INCLUDED
#include "BGK1D.h"
#endif // _BGK1D_INCLUDED

/********************************************************
 * Routine: BGK1DSolver                               *
 *                                                      *
 * Computes solutions to 1D BGK equations.            *
 *                                                      *
 ********************************************************/
int BGK1DSolver(char *Input_File_Name_ptr,
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

  /* Solution variables. */

  BGK1D_UniformMesh *Soln_ptr = NULL;

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

  // set batch_mode if required
  if (batch_flag){
    Input_Parameters.Verbose() = OFF;
  }

  strcpy(Input_Parameters.Input_File_Name, Input_File_Name_ptr);
  Open_Input_File(Input_Parameters);
  if (Input_Parameters.Input_File.fail()) {
     if (batch_flag) {
        cout << "\n BGK1D ERROR: Unable to open BGK1D input data file.\n\n";
     } else {
        cout << "\n BGK1D ERROR: Unable to open BGK1D input data file.";
        cout << "\n\nBGK1D: Execution terminated.\n";
     } /* endif */
     return (-1);
  } /* endif */

  if (! batch_flag) {
    cout << "\n Reading BGK1D input data file `"
	 << Input_Parameters.Input_File_Name << "'.";
    cout.flush();
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
             cout << "\n BGK1D ERROR: Error reading BGK1D data at line #"
	          << -line_number << " of input data file.\n\n";
         } else {
             cout << "\n BGK1D ERROR: Error reading BGK1D data at line #"
	          << -line_number  << " of input data file.";
             cout << "\n\nBGK1D: Execution terminated.\n";
         } /* end if */
         return (line_number);
     } /* endif */
  } /* endwhile */

  if (! batch_flag)  {
     cout << Input_Parameters << "\n";
     cout.flush();
  } /* endif */

  /*********************************************************
   * Create mesh and allocate BGK1D solution variables.  *
   *********************************************************/

  execute_new_calculation: ;

  /* Setup distribution function discretization used in calculation. */
  BGK1D_Vector::setup(Input_Parameters.bgk_v_number,
		      Input_Parameters.bgk_v_min,
		      Input_Parameters.bgk_v_max);

  /* Allocate memory for 1D BGK equation solution on
     uniform mesh. */

  if (! batch_flag){
    cout << "\n Creating memory for BGK1D solution variables.";
    cout.flush();
  }
  Soln_ptr=Allocate(Soln_ptr,
                    Input_Parameters);

  if (Soln_ptr == NULL){
    cout << "\n BGK1DSolvers::Allocate() Error! Probably not enough memory!";
    cout.flush();
    exit(1);
  }

  /* Create uniform mesh. */

  if (! batch_flag){
    cout << "\n Creating uniform mesh.";
    cout.flush();
  }
  Grid(Soln_ptr,
       Input_Parameters.X_Min,
       Input_Parameters.X_Max,
       Input_Parameters.Number_of_Cells);

  /********************************************************
   * Initialize BGK1D solution variables.               *
   ********************************************************/

  /* Set the initial time level. */

  time = ZERO;
  number_of_time_steps = 0;

  /* Initialize the conserved and primitive state
     solution variables. */

  if (! batch_flag){
    cout << "\n Prescribing BGK1D initial data.";
    cout.flush();
  }
  error_flag = ICs(Soln_ptr,
		   "ZB",
		   Input_Parameters.i_ICs,
		   Input_Parameters.Number_of_Cells,
		   Input_Parameters);

  if(error_flag) {
    cout << "Unable to set desired ICs, exiting." << endl;
    return 1;
  }

  /********************************************************
   * Solve conservation form of 1D BGK equations for    *
   * specified IBVP or BVP on uniform mesh.               *
   ********************************************************/

  continue_existing_calculation: ;

  if ((Input_Parameters.Local_Time_Stepping &&
       Input_Parameters.Maximum_Number_of_Time_Steps > 0) ||
      (!Input_Parameters.Local_Time_Stepping &&
       Input_Parameters.Time_Max > time)) {
    if (! batch_flag){
      cout << "\n\n Beginning BGK1D computations.\n\n";
      cout.flush();
    }
     while (1) {
         /* Determine local and global time steps. */
         dtime = CFL(Soln_ptr, Input_Parameters.Number_of_Cells);

         if (time + Input_Parameters.CFL_Number*dtime >
             Input_Parameters.Time_Max) {
            dtime = (Input_Parameters.Time_Max-time)/Input_Parameters.CFL_Number;
         } /* endif */

         /* Output progress information for the calculation. */
         if (! batch_flag && number_of_time_steps == 0 ) {
             cout << " Time Step = " << number_of_time_steps
	          << " Time = " << time*THOUSAND << "\n .";
	     cout.flush();
         } else if (! batch_flag && number_of_time_steps%100 == 0 ) {
	     cout << "\n" << " Time Step = " << number_of_time_steps
	          << " Time = " << time*THOUSAND << "\n .";
	     cout.flush();
         } else if (! batch_flag && number_of_time_steps%50 == 0 ) {
	     cout << "\n .";
	     cout.flush();
         } else if (! batch_flag) {
	     cout << ".";
	     cout.flush();
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
             error_flag = dVdt_explicitEuler_upwind(Soln_ptr,
						    Input_Parameters.Number_of_Cells,
						    dtime,
						    Input_Parameters.CFL_Number,
						    Input_Parameters.i_Flux_Function,
						    Input_Parameters.Local_Time_Stepping);
             break;
           case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR:
             error_flag = dVdt_2stage_2ndOrder_upwind(Soln_ptr,
						      Input_Parameters,
						      dtime,
						      Input_Parameters.CFL_Number,
						      Input_Parameters.i_Reconstruction,
						      Input_Parameters.i_Limiter,
						      Input_Parameters.i_Flux_Function,
						      Input_Parameters.Local_Time_Stepping);
             break;
           default:
	     cout << "Invalid time marching chosen!" << endl;
	     return (1);
             break;
         } /* endswitch */

         if (error_flag) {
             if (batch_flag) {
                 cout << "\nCFFC ERROR: BGK1D solution error.\n\n";
		 cout.flush();
             } else {
                 cout << "\n CFFC ERROR: BGK1D solution error.";
                 cout << "\n\nCFFC: Execution terminated.\n";
		 cout.flush();
             } /* endif */
             return (error_flag);
         } /* endif */

         /* Update time and time step counter. */
         number_of_time_steps = number_of_time_steps + 1;
         time = time + Input_Parameters.CFL_Number*dtime;

     } /* endwhile */
     if (! batch_flag){
       cout << "\n\n BGK1D computations complete.\n";
       cout.flush();
     }
  } /* endif */

  /********************************************************
   * Write 1D BGK solution to output file.              *
   ********************************************************/

  while (1) {
     Get_Next_Input_Control_Parameter(Input_Parameters);
     command_flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
     line_number = Input_Parameters.Line_Number;
     if (command_flag == EXECUTE_CODE) {
         // Deallocate memory for 1D BGK equation solution.
         if (! batch_flag) cout << "\n Deallocating BGK1D solution variables.";
         Soln_ptr=Deallocate(Soln_ptr);
         // Execute new calculation.
         if (! batch_flag)  {
            cout << "\n\n Starting a new calculation.";
            cout << Input_Parameters << "\n";
         } /* endif */
         goto execute_new_calculation;

     } else if (command_flag == TERMINATE_CODE) {
         // Deallocate memory for 1D BGK equation solution.
         if (! batch_flag) cout << "\n Deallocating BGK1D solution variables.";
         Soln_ptr=Deallocate(Soln_ptr);
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
           cout << "\n Writing BGK1D solution to output data file `"
	        << Output_File_Name_ptr << "'.";
         Output_File.open(Output_File_Name_ptr, ios::out);
         if (Output_File.fail()) {
            if (batch_flag) {
               cout << "\nCFFC ERROR: Unable to open BGK1D output data file.\n\n";
            } else {
               cout << "\n CFFC ERROR: Unable to open BGK1D output data file.";
               cout << "\n\nCFFC: Execution terminated.\n";
            } /* endif */
            return (1);
         } /* endif */

	 Output_Tecplot(Soln_ptr,
			Input_Parameters,
			number_of_time_steps,
			time,
			Output_File);
         Output_File.close();

     } else if (command_flag == INVALID_INPUT_CODE ||
                command_flag == INVALID_INPUT_VALUE) {
         line_number = -line_number;
         if (batch_flag) {
             cout << "\nBGK1D ERROR: Error reading BGK1D data at line #"
	          << -line_number << " of input data file.\n\n";
         } else {
             cout << "\nBGK1D ERROR: Error reading BGK1D data at line #"
	          << -line_number  << " of input data file.";
             cout << "\n\nCFFC: Execution terminated.\n";
         } /* end if */
         return (line_number);
     } /* endif */
  } /* endwhile */

  /********************************************************
   * End of all BGK1DSolver computations and I/O.         *
   ********************************************************/

  return (0);

}
