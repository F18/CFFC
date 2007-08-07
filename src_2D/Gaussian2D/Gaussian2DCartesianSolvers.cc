/* Gaussian2DCartesianSolvers.cc:  2D Gaussian Equation Cartesian Mesh Solvers. */

/* Include 2D Gaussian Cartesian mesh solution header file. */

#ifndef _GAUSSIAN2D_CARTESIAN_INCLUDED
#include "Gaussian2DCartesian.h"
#endif // _GAUSSIAN2D_CARTESIAN_INCLUDED

/********************************************************
 * Routine: Euler2DCartesianSolver                      *
 *                                                      *
 * Computes solutions to 2D Gaussian equations on 2D    *
 * uniform Cartesian mesh.                              *
 *                                                      *
 ********************************************************/
int Gaussian2DCartesianSolver(char *Input_File_Name_ptr,
	                   int batch_flag) {

  /********************************************************  
   * Local variable declarations.                         *
   ********************************************************/

  // Gaussian2D input variables and parameters:
  Gaussian2D_Input_Parameters Input_Parameters;

  // Output file name pointer:
  char *Output_File_Name_ptr;

  // Gnuplot file name pointer:
  char *Gnuplot_File_Name_ptr;

  // Output file stream:
  ofstream Output_File;

  // Gnuplot macro file stream:
  ofstream Gnuplot_File;

  /* Uniform Cartesian mesh solution variables. */

  Gaussian2D_Cartesian_Cell **Soln_ptr;
  Gaussian2D_pState W_left, W_right, W_bottom, W_top;

  /* Other local solution variables. */

  int number_of_time_steps,
      command_flag, error_flag, line_number, i;
  int bc_type_left, bc_type_right, 
      bc_type_bottom, bc_type_top;

  double time, dtime;
  
  /********************************************************  
   * Set default values for the input solution parameters *
   * and then read user specified input values from the   *
   * specified input file.                                *
   ********************************************************/

  Set_Default_Input_Parameters(Input_Parameters);
  // Set time accurate mode with no local time stepping.
  Input_Parameters.Time_Accurate = 1;
  Input_Parameters.Local_Time_Stepping = 0;

  strcpy(Input_Parameters.Input_File_Name, Input_File_Name_ptr);
  Open_Input_File(Input_Parameters);
  if (Input_Parameters.Input_File.bad()) {
     if (batch_flag) {
        cout << "\nPDES++ ERROR: Unable to open Gaussian2D input data file.\n\n";
     } else {
        cout << "\n PDES++ ERROR: Unable to open Gaussian2D input data file.";
        cout << "\n\nPDES++: Execution terminated.\n";
     } /* endif */
     return (-1);
  } /* endif */

  if (! batch_flag) {
    cout << "\n Reading Gaussian2D input data file `"
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
             cout << "\nPDES++ ERROR: Error reading Gaussian2D data at line #"
	          << -line_number << " of input data file.\n\n";
         } else {
             cout << "\n PDES++ ERROR: Error reading Gaussian2D data at line #"
	          << -line_number  << " of input data file.";
             cout << "\n\nPDES++: Execution terminated.\n";
         } /* end if */
         return (line_number);
     } /* endif */
  } /* endwhile */

  // Force Cartesian grid type.
  strcpy(Input_Parameters.Grid_Type, "Cartesian");
  Input_Parameters.i_Grid = GRID_CARTESIAN_UNIFORM;

  if (! batch_flag)  {
     cout << Input_Parameters << "\n";
  } /* endif */

  /***********************************************************  
   * Create mesh and allocate Gaussian2D solution variables. *
   ***********************************************************/

  execute_new_calculation: ;

  /* Allocate memory for 2D Gaussian equation solution for
     uniform Cartesian mesh. */

  if (! batch_flag) cout << "\n Creating memory for Gaussian2D solution variables.";
  Soln_ptr=Allocate(Soln_ptr, 
                    Input_Parameters.Number_of_Cells_Idir, 
                    Input_Parameters.Number_of_Cells_Jdir);

  /* Create uniform Cartesian mesh. */

  if (! batch_flag) cout << "\n Creating uniform Cartesian mesh.";
  if (Input_Parameters.i_ICs == IC_SHOCK_BOX) {
      Grid(Soln_ptr, 
           ZERO, 
           Input_Parameters.Box_Width, 
           ZERO, 
           Input_Parameters.Box_Height, 
           Input_Parameters.Number_of_Cells_Idir, 
           Input_Parameters.Number_of_Cells_Jdir);
  } else {
      Grid(Soln_ptr, 
           -HALF*Input_Parameters.Box_Width, 
           HALF*Input_Parameters.Box_Width, 
           -HALF*Input_Parameters.Box_Height, 
           HALF*Input_Parameters.Box_Height, 
           Input_Parameters.Number_of_Cells_Idir, 
           Input_Parameters.Number_of_Cells_Jdir);
  } /* endif */

  /********************************************************  
   * Initialize Gaussian2D solution variables.            *
   ********************************************************/

  /* Set the initial time level. */

  time = ZERO;
  number_of_time_steps = 0;
  
  /* Initialize the conserved and primitive state
     solution variables. */
  
  if (! batch_flag) cout << "\n Prescribing Gaussian2D initial data.";
  ICs(Soln_ptr, 
      Input_Parameters.Gas_Type, 
      Input_Parameters.i_ICs,
      Input_Parameters.Number_of_Cells_Idir, 
      Input_Parameters.Number_of_Cells_Jdir);
  
  /********************************************************  
   * Set Gaussian2D boundary conditions.                  *
   ********************************************************/

  /* Set the type of boundary conditions to apply at the boundaries
     of the numerical grid. */

  if (! batch_flag) cout << "\n Prescribing Gaussian2D boundary data.";
  if (Input_Parameters.i_ICs == IC_SHOCK_BOX) {
      bc_type_left = BC_REFLECTION;
      bc_type_right = BC_REFLECTION;
      bc_type_bottom = BC_REFLECTION;
      bc_type_top = BC_REFLECTION;
  } else {
      bc_type_left = BC_CONSTANT_EXTRAPOLATION;
      bc_type_right = BC_CONSTANT_EXTRAPOLATION;
      bc_type_bottom = BC_CONSTANT_EXTRAPOLATION;
      bc_type_top = BC_CONSTANT_EXTRAPOLATION;
  } /* endif */

  /* Set boundary solution states. */

  W_left = Soln_ptr[0][1].W;
  W_right = Soln_ptr[Input_Parameters.Number_of_Cells_Idir+1][1].W;
  W_bottom = Soln_ptr[1][0].W;
  W_top = Soln_ptr[1][Input_Parameters.Number_of_Cells_Jdir+1].W;

  /***********************************************************  
   * Solve IBVP or BVP for conservation form of 2D Gaussian  *
   * equations on uniform Cartesian mesh.                    *
   ***********************************************************/

  continue_existing_calculation: ;

  if ((Input_Parameters.Local_Time_Stepping &&
       Input_Parameters.Maximum_Number_of_Time_Steps > 0) ||
      (!Input_Parameters.Local_Time_Stepping &&
       Input_Parameters.Time_Max > time)) {
     if (! batch_flag) cout << "\n\n Beginning Gaussian2D computations.\n\n";
     while (1) {
         /* Determine local and global time steps. */
         dtime = CFL(Soln_ptr, 
                     Input_Parameters.Number_of_Cells_Idir, 
                     Input_Parameters.Number_of_Cells_Jdir);
         if (time + Input_Parameters.CFL_Number*dtime > 
             Input_Parameters.Time_Max) {
            dtime = (Input_Parameters.Time_Max-time)/Input_Parameters.CFL_Number;
         } /* endif */

         /* Determine local and global time steps. */
         if (! batch_flag && number_of_time_steps ==0 ) {
             cout << " Time Step = " << number_of_time_steps
	          << " Time = " << time*THOUSAND << "\n .";
         } else if (! batch_flag && 
                    number_of_time_steps-100*(number_of_time_steps/100)==0 ) {
	     cout << "\n" << " Time Step = " << number_of_time_steps
	          << " Time = " << time*THOUSAND << "\n .";
         } else if (! batch_flag && 
                    number_of_time_steps-50*(number_of_time_steps/50)==0 ) {
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
                           Input_Parameters.Number_of_Cells_Idir,
                           Input_Parameters.Number_of_Cells_Jdir,
                           bc_type_left, bc_type_right,
                           bc_type_bottom, bc_type_top,
                           W_left, W_right,
                           W_bottom, W_top,
          		   dtime,
                           Input_Parameters.CFL_Number,
                           Input_Parameters.i_Flux_Function,
       			   Input_Parameters.Local_Time_Stepping);
             break;
           case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR:
             error_flag = dUdt_2stage_2ndOrder_upwind(Soln_ptr,
                           Input_Parameters.Number_of_Cells_Idir,
                           Input_Parameters.Number_of_Cells_Jdir,
                           bc_type_left, bc_type_right,
                           bc_type_bottom, bc_type_top,
                           W_left, W_right,
                           W_bottom, W_top,
        		   dtime,
                           Input_Parameters.CFL_Number,
                           Input_Parameters.i_Reconstruction,
                           Input_Parameters.i_Limiter, 
                           Input_Parameters.i_Flux_Function,
       			   Input_Parameters.Local_Time_Stepping);
             break;
           default:
             error_flag = dUdt_explicitEuler_upwind(Soln_ptr,
                           Input_Parameters.Number_of_Cells_Idir,
                           Input_Parameters.Number_of_Cells_Jdir,
                           bc_type_left, bc_type_right,
                           bc_type_bottom, bc_type_top,
                           W_left, W_right,
                           W_bottom, W_top,
          		   dtime,
                           Input_Parameters.CFL_Number,
                           Input_Parameters.i_Flux_Function,
       			   Input_Parameters.Local_Time_Stepping);
             break;
         } /* endswitch */

         if (error_flag) {
             if (batch_flag) {
                 cout << "\nPDES++ ERROR: Gaussian2D solution error.\n\n";
             } else {
                 cout << "\n PDES++ ERROR: Gaussian2D solution error.";
                 cout << "\n\nPDES++: Execution terminated.\n";
             } /* endif */
             return (error_flag);
         } /* endif */
      
         /* Update time and time step counter. */
         number_of_time_steps = number_of_time_steps + 1;
         time = time + Input_Parameters.CFL_Number*dtime;

     } /* endwhile */
     if (! batch_flag) cout << "\n\n Gaussian2D computations complete.\n";
  } /* endif */

  /********************************************************  
   * Solution calculations complete.                      *
   * Write 2D Gaussian solution to output file as         *
   * required, reset solution parameters, and run other   *
   * cases as specified by input data.                    *
   ********************************************************/

  while (1) {
     Get_Next_Input_Control_Parameter(Input_Parameters);
     command_flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
     line_number = Input_Parameters.Line_Number;
     if (command_flag == EXECUTE_CODE) {
         // Deallocate memory for 2D Euler equation solution.
         if (! batch_flag) cout << "\n Deallocating Gaussian2D solution variables.";
         Soln_ptr=Deallocate(Soln_ptr, 
                             Input_Parameters.Number_of_Cells_Idir, 
                             Input_Parameters.Number_of_Cells_Jdir);
         // Execute new calculation.
         if (! batch_flag)  {
            cout << "\n\n Starting a new calculation.";
            cout << Input_Parameters << "\n";
         } /* endif */
         goto execute_new_calculation;

     } else if (command_flag == TERMINATE_CODE) {
         // Deallocate memory for 2D Euler equation solution.
         if (! batch_flag) cout << "\n Deallocating Gaussian2D solution variables.";
         Soln_ptr=Deallocate(Soln_ptr, 
                             Input_Parameters.Number_of_Cells_Idir, 
                             Input_Parameters.Number_of_Cells_Jdir);
         // Close input data file.
         if (! batch_flag) cout << "\n\n Closing Gaussian2D input data file.";
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
           cout << "\n Writing Gaussian2D solution to output data file `"
	        << Output_File_Name_ptr << "'.";
         Output_File.open(Output_File_Name_ptr, ios::out);
         if (Output_File.bad()) {
            if (batch_flag) {
               cout << "\nPDES++ ERROR: Unable to open Gaussian2D output data file.\n\n";
            } else {
               cout << "\n PDES++ ERROR: Unable to open Gaussian2D output data file.";
               cout << "\n\nPDES++: Execution terminated.\n";
            } /* endif */
            return (1);
         } /* endif */

         if (Input_Parameters.i_Output_Format == IO_GNUPLOT) {
            Output_Gnuplot(Soln_ptr,
                           Input_Parameters.Number_of_Cells_Idir,
                           Input_Parameters.Number_of_Cells_Jdir,
		           number_of_time_steps,
		           time,
		           Output_File);
         } else if (Input_Parameters.i_Output_Format == IO_TECPLOT) {
            Output_Tecplot(Soln_ptr,
                           Input_Parameters.Number_of_Cells_Idir,
                           Input_Parameters.Number_of_Cells_Jdir,
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
                  cout << "\nPDES++ ERROR: Unable to open Gaussian2D Gnuplot macro file.\n\n";
               } else {
                  cout << "\n PDES++ ERROR: Unable to open Gaussian2D Gnuplot macro file.";
                  cout << "\n\nPDES++: Execution terminated.\n";
               } /* endif */
               return (2);
            } /* endif */

            Gnuplot_File << "set title \"PDES++: 2D Gaussian Solution\"\n"
                         << "set xlabel \"x\"\n"
                         << "set ylabel \"y\"\n"
	                 << "splot \"" << Output_File_Name_ptr << "\""
                         << " using 1:2:3 \"%*lf%*lf%lf%lf%lf%*lf%*lf%*lf%*lf\" \\\n"
	                 << "     title \"density\" with points\n"
                         << "pause -1  \"Hit return to continue\"\n"
                         << "splot \"" << Output_File_Name_ptr << "\""
                         << " using 1:2:3 \"%*lf%*lf%lf%lf%*lf%lf%*lf%*lf%*lf\" \\\n"
	                 << "     title \"velocity x-direction\" with points\n"
                         << "pause -1  \"Hit return to continue\"\n"
                         << "splot \"" << Output_File_Name_ptr << "\""
                         << " using 1:2:3 \"%*lf%*lf%lf%lf%*lf%*lf%lf%*lf%*lf\" \\\n"
	                 << "     title \"velocity y-direction\" with points\n"
                         << "pause -1  \"Hit return to continue\"\n"
                         << "splot \"" << Output_File_Name_ptr << "\""
                         << " using 1:2:3 \"%*lf%*lf%lf%lf%*lf%*lf%*lf%lf%*lf\" \\\n"
	                 << "     title \"pressure\" with points\n"
                         << "pause -1  \"Hit return to continue\"\n"
                         << "splot \"" << Output_File_Name_ptr << "\""
                         << " using 1:2:3 \"%*lf%*lf%lf%lf%*lf%*lf%*lf%*lf%lf\" \\\n"
	                 << "     title \"temperature\" with points\n"
                         << "pause -1  \"Hit return to continue\"\n";

            Gnuplot_File.close();
         } /* endif */

     } else if (command_flag == INVALID_INPUT_CODE ||
                command_flag == INVALID_INPUT_VALUE) {
         line_number = -line_number;
         if (batch_flag) {
             cout << "\nPDES++ ERROR: Error reading Gaussian2D data at line #"
	          << -line_number << " of input data file.\n\n";
         } else {
             cout << "\n PDES++ ERROR: Error reading Gaussian2D data at line #"
	          << -line_number  << " of input data file.";
             cout << "\n\nPDES++: Execution terminated.\n";
         } /* end if */
         return (line_number);
     } /* endif */
  } /* endwhile */

  /***********************************************************  
   * End of all Gaussian2DSolver_Cartesian computations and  *
   * I/O.                                                    *
   ***********************************************************/

  return (0);
  
}

