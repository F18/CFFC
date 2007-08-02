/* HyperHeat1DSolvers.cc:  1D Hyperbolic Heat Equations Solvers. */

/* Include 1D hyperbolic heat equations solution header file. */

#ifndef _HYPERHEAT1D_INCLUDED
#include "HyperHeat1D.h"
#endif // _HYPERHEAT1D_INCLUDED

/********************************************************
 * Routine: HyperHeat1DSolver                           *
 *                                                      *
 * Computes solutions to 1D Hyperbolic heat equations.  *
 *                                                      *
 ********************************************************/
int HyperHeat1DSolver(char *Input_File_Name_ptr,
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

  HyperHeat1D_UniformMesh *Soln_ptr;
  HyperHeat1D_State U_left, U_right;

  /* Other local solution variables. */

  int number_of_time_steps,
      command_flag, error_flag, line_number, i;
  int bc_type_left, bc_type_right;

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
   * Create mesh and allocate HyperHeart1D solution        *
   * variables.                                            *
   *********************************************************/

  execute_new_calculation: ;

  /* Allocate memory for 1D hyperbolic heat equations solution on
     uniform mesh. */

  if (! batch_flag) cout << "\n Creating memory for HyperHeat1D solution variables.";
  Soln_ptr=Allocate(Soln_ptr,
                    Input_Parameters.Number_of_Cells);

  /* Create uniform mesh. */

  if (! batch_flag) cout << "\n Creating uniform mesh.";
  Grid(Soln_ptr,
       Input_Parameters.X_Min,
       Input_Parameters.X_Max, 
       Input_Parameters.Number_of_Cells);  

  /********************************************************  
   * Initialize HyperHeat1D solution variables.           *
   ********************************************************/

  /* Set the initial time level. */

  time = ZERO;
  number_of_time_steps = 0;
  
  /* Initialize solution variables. */

  if (! batch_flag) cout << "\n Prescribing HyperHeat1D initial data.";
  ICs(Soln_ptr, 
      Input_Parameters.i_ICs,
      Input_Parameters.Kappa, 
      Input_Parameters.Tau, 
      Input_Parameters.Number_of_Cells);

  /********************************************************  
   * Set HyperHeat1D boundary conditions.                 *
   ********************************************************/

  if (! batch_flag) cout << "\n Prescribing HyperHeat1D boundary data.";

  /* Set the type of boundary conditions to apply at either end
     of the numerical grid. */

  if (Input_Parameters.i_ICs == IC_UNIFORM ||
      Input_Parameters.i_ICs == IC_RIEMANN_IVP_QX0 ||
      Input_Parameters.i_ICs == IC_RIEMANN_IVP_T0 ||
      Input_Parameters.i_ICs == IC_RIEMANN_IVP) {
      bc_type_left = BC_CONSTANT_EXTRAPOLATION;
      bc_type_right = BC_CONSTANT_EXTRAPOLATION;
  } else if (Input_Parameters.i_ICs == IC_IMPULSIVE_ROD ||
             Input_Parameters.i_ICs == IC_SINUSOIDAL_ROD1 ||
             Input_Parameters.i_ICs == IC_SINUSOIDAL_ROD4) {
      bc_type_left = BC_FIXED_TEMP;
      bc_type_right = BC_FIXED_TEMP;
  } /* endif */
  
  /* Set left and right boundary solution states. */

  U_left = Soln_ptr[0].U;
  U_right = Soln_ptr[Input_Parameters.Number_of_Cells+1].U;

  /********************************************************
   * Solve weak conservation form of 1D hyperbolic heat   *
   * equations for specified IBVP or BVP on uniform mesh. *
   ********************************************************/

  continue_existing_calculation: ;

  if ((Input_Parameters.Local_Time_Stepping &&
       Input_Parameters.Maximum_Number_of_Time_Steps > 0) ||
      (!Input_Parameters.Local_Time_Stepping &&
       Input_Parameters.Time_Max > time)) {
     if (! batch_flag) cout << "\n\n Beginning HyperHeat1D computations.\n\n";
     while (1) {
         /* Determine local and global time steps. */
         dtime = CFL(Soln_ptr, 
                     Input_Parameters.Number_of_Cells);
         if (time + Input_Parameters.CFL_Number*dtime > Input_Parameters.Time_Max) {
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
                           bc_type_left, 
                           bc_type_right,
                           U_left, 
                           U_right,
          		   dtime,
                           Input_Parameters.CFL_Number,
                           Input_Parameters.i_Flux_Function,
       			   Input_Parameters.Local_Time_Stepping);
             break;
          case TIME_STEPPING_SEMI_IMPLICIT_EULER:
             error_flag = dUdt_semiimplicitEuler_upwind(Soln_ptr,
                           Input_Parameters.Number_of_Cells,
                           bc_type_left, 
                           bc_type_right,
                           U_left, 
                           U_right,
        		   dtime,
                           Input_Parameters.CFL_Number,
                           Input_Parameters.i_Flux_Function,
       			   Input_Parameters.Local_Time_Stepping);
             break;
           case TIME_STEPPING_SEMI_IMPLICIT_PREDICTOR_CORRECTOR:
             error_flag = dUdt_2stage_2ndOrder_upwind(Soln_ptr,
                           Input_Parameters.Number_of_Cells,
                           bc_type_left, 
                           bc_type_right,
                           U_left, 
                           U_right,
        		   dtime,
                           Input_Parameters.CFL_Number,
                           Input_Parameters.i_Reconstruction,
                           Input_Parameters.i_Limiter,
                           Input_Parameters.i_Flux_Function,
       			   Input_Parameters.Local_Time_Stepping);
             break;
           default:
             error_flag = dUdt_explicitEuler_upwind(Soln_ptr,
                           Input_Parameters.Number_of_Cells,
                           bc_type_left, 
                           bc_type_right,
                           U_left, 
                           U_right,
          		   dtime,
                           Input_Parameters.CFL_Number,
                           Input_Parameters.i_Flux_Function,
       			   Input_Parameters.Local_Time_Stepping);
             break;
         } /* endswitch */

         if (error_flag) {
             if (batch_flag) {
                 cout << "\nPDES++ ERROR: HyperHeat1D solution error.\n\n";
             } else {
                 cout << "\n PDES++ ERROR: HyperHeat1D solution error.";
                 cout << "\n\nPDES++: Execution terminated.\n";
             } /* endif */
             return (error_flag);
         } /* endif */
      
         /* Update time and time step counter. */
         number_of_time_steps = number_of_time_steps + 1;
         time = time + Input_Parameters.CFL_Number*dtime;

     } /* endwhile */
     if (! batch_flag) cout << "\n\n HyperHeat1D computations complete.\n";
  } /* endif */

  /********************************************************  
   * Solution calculations complete.                      *
   * Write 1D hyperbolic heat equations solution          *
   * to output file as required, reset solution           *
   * parameters, and run other cases as specified by      *
   * input data.                                          *
   ********************************************************/

  while (1) {
     Get_Next_Input_Control_Parameter(Input_Parameters);
     command_flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
     line_number = Input_Parameters.Line_Number;
     if (command_flag == EXECUTE_CODE) {
         // Deallocate memory for 1D hyperbolic heat equations solution.
         if (! batch_flag) cout << "\n Deallocating HyperHeat1D solution variables.";
         Soln_ptr=Deallocate(Soln_ptr,
                             Input_Parameters.Number_of_Cells);
         // Execute new calculation.
         if (! batch_flag)  {
            cout << "\n\n Starting a new calculation.";
            cout << Input_Parameters << "\n";
         } /* endif */
         goto execute_new_calculation;

     } else if (command_flag == TERMINATE_CODE) {
         // Deallocate memory for 1D hyperbolic heat equations solution.
         if (! batch_flag) cout << "\n Deallocating HyperHeat1D solution variables.";
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
           cout << "\n Writing HyperHeat1D solution to output data file `"
	        << Output_File_Name_ptr << "'.";
         Output_File.open(Output_File_Name_ptr, ios::out);
         if (Output_File.bad()) {
            if (batch_flag) {
               cout << "\nPDES++ ERROR: Unable to open HyperHeat1D output data file.\n\n";
            } else {
               cout << "\n PDES++ ERROR: Unable to open HyperHeat1D output data file.";
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
                  cout << "\nPDES++ ERROR: Unable to open HyperHeat1D Gnuplot macro file.\n\n";
               } else {
                  cout << "\n PDES++ ERROR: Unable to open HyperHeat1D Gnuplot macro file.";
                  cout << "\n\nPDES++: Execution terminated.\n";
               } /* endif */
               return (2);
            } /* endif */

            Gnuplot_File << "set title \"PDES++: 1D Hyperbolic Heat Equations Solution\"\n"
                         << "set xlabel \"x\"\n"
	                 << "plot \"" << Output_File_Name_ptr << "\""
                         << " using 1:2 \"%*lf%lf%*lf%lf%*lf\" \\\n"
	                 << "     title \"T\" with points\n"
                         << "pause -1  \"Hit return to continue\"\n"
                         << "plot \"" << Output_File_Name_ptr << "\""
                         << " using 1:2 \"%*lf%lf%*lf%*lf%lf\" \\\n"
	                 << "     title \"qx\" with points\n"
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
   * End of all HyperHeat1DSolver computations and I/O.   *
   ********************************************************/

  return (0);
  
}
