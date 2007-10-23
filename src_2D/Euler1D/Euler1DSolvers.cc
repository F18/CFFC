/*!\file Euler1DSolvers.cc
  \brief 1D Euler Equation Solvers. */

/* Include 1D Euler solution header file. */

#ifndef _EULER1D_INCLUDED
#include "Euler1D.h"
#endif // _EULER1D_INCLUDED

#include "Euler1D_HighOrder.h"	// High-order 1D Euler header file 

/********************************************************
 * Routine: Euler1DSolver                               *
 *                                                      *
 * Computes solutions to 1D Euler equations.            *
 *                                                      *
 ********************************************************/
int Euler1DSolver(char *Input_File_Name_ptr,
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

  Euler1D_UniformMesh *Soln_ptr = NULL;

 /* Other local solution variables. */

  int number_of_time_steps,
      command_flag, error_flag, line_number, i;

  /* Accuracy variables */
  vector<double> LNorms; LNorms.reserve(3);
  bool AccuracyAssessed_flag = false;
  bool Title_Error_Norms = true;


  /* Create memory storage variable for the ENO subroutines */
  MemoryStorageENO_1D MemoryENO;

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
        cout << "\n Euler1D ERROR: Unable to open Euler1D input data file.\n\n";
     } else {
        cout << "\n Euler1D ERROR: Unable to open Euler1D input data file.";
        cout << "\n\nEuler1D: Execution terminated.\n";
     } /* endif */
     return (-1);
  } /* endif */

  if (! batch_flag) {
    cout << "\n Reading Euler1D input data file `"
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
             cout << "\n Euler1D ERROR: Error reading Euler1D data at line #"
	          << -line_number << " of input data file.\n\n";
         } else {
             cout << "\n Euler1D ERROR: Error reading Euler1D data at line #"
	          << -line_number  << " of input data file.";
             cout << "\n\nEuler1D: Execution terminated.\n";
         } /* end if */
         return (line_number);
     } /* endif */
  } /* endwhile */

  if (! batch_flag)  {
     cout << Input_Parameters << "\n";
  } /* endif */

  /*********************************************************  
   * Create mesh and allocate Euler1D solution variables.  *
   *********************************************************/

  execute_new_calculation: ;

  /* Allocate memory for 1D Euler equation solution on
     uniform mesh. */

  if (! batch_flag) cout << "\n Creating memory for Euler1D solution variables.";
  Soln_ptr=Allocate(Soln_ptr,
                    Input_Parameters);

  if (Soln_ptr == NULL){
    cout << "\n Euler1DSolvers::Allocate() Error! Probably not enough memory!";
    exit(1);
  }

  /* Create uniform mesh. */

  if (! batch_flag) cout << "\n Creating uniform mesh.";
  Grid(Soln_ptr, 
       Input_Parameters.X_Min, 
       Input_Parameters.X_Max, 
       Input_Parameters.Number_of_Cells);

  /********************************************************  
   * Initialize Euler1D solution variables.               *
   ********************************************************/

  /* Set the initial time level. */

  time = ZERO;
  number_of_time_steps = 0;
  
  /* Initialize the conserved and primitive state
     solution variables. */
  
  if (! batch_flag) cout << "\n Prescribing Euler1D initial data.";
  ICs(Soln_ptr, 
      "AIR", 
      Input_Parameters.i_ICs, 
      Input_Parameters.Number_of_Cells,
      Input_Parameters);
  

  /********************************************************  
   * Solve conservation form of 1D Euler equations for    *
   * specified IBVP or BVP on uniform mesh.               *
   ********************************************************/

  continue_existing_calculation: ;

  if ((Input_Parameters.Local_Time_Stepping &&
       Input_Parameters.Maximum_Number_of_Time_Steps > 0) ||
      (!Input_Parameters.Local_Time_Stepping &&
       Input_Parameters.Time_Max > time)) {
     if (! batch_flag) cout << "\n\n Beginning Euler1D computations.\n\n";
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
	          << " Time = " << time*THOUSAND << "\n .";
         } else if (! batch_flag && 
                    number_of_time_steps-100*(number_of_time_steps/100) == 0 ) {
	     cout << "\n" << " Time Step = " << number_of_time_steps
	          << " Time = " << time*THOUSAND << "\n .";
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
	   case TIME_STEPPING_EXPLICIT_EULER_HIGH_ORDER:
	     error_flag = dUdt_explicitEuler_upwind(Soln_ptr,
						    Input_Parameters,
						    dtime,
						    Input_Parameters.Local_Time_Stepping,
						    &Euler1D_UniformMesh::CellHighOrder);
	     break;
           case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR:
             error_flag = dUdt_2stage_2ndOrder_upwind(Soln_ptr,
						      Input_Parameters,
						      dtime,
						      Input_Parameters.CFL_Number,
						      Input_Parameters.i_Reconstruction,
						      Input_Parameters.i_Limiter, 
						      Input_Parameters.i_Flux_Function,
						      Input_Parameters.Local_Time_Stepping);
             break;
	   case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR_HIGH_ORDER:
	     error_flag = dUdt_2stage_HighOrder_upwind(Soln_ptr,
						       Input_Parameters,
						       dtime,
						       Input_Parameters.Local_Time_Stepping,
						       &Euler1D_UniformMesh::CellHighOrder);
	     break;
	   case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA_4_HIGH_ORDER:
	     error_flag = dUdt_4stage_HighOrder_upwind(Soln_ptr,
						       Input_Parameters,
						       dtime,
						       Input_Parameters.Local_Time_Stepping,
						       &Euler1D_UniformMesh::CellHighOrder);
	     break;
           case TIME_STEPPING_LAX_FRIEDRICHS:
             error_flag = dUdt_LaxFriedrichs(Soln_ptr,
					     Input_Parameters.Number_of_Cells,
					     dtime,
					     Input_Parameters.CFL_Number,
					     Input_Parameters.Local_Time_Stepping);
             break;
           case TIME_STEPPING_LAX_WENDROFF:
             error_flag = dUdt_LaxWendroff(Soln_ptr,
					   Input_Parameters.Number_of_Cells,
					   dtime,
					   Input_Parameters.CFL_Number,
					   Input_Parameters.Local_Time_Stepping);
             break;
           case TIME_STEPPING_MACCORMACK:
             error_flag = dUdt_MacCormack(Soln_ptr,
					  Input_Parameters.Number_of_Cells,
					  dtime,
					  Input_Parameters.CFL_Number,
					  Input_Parameters.Local_Time_Stepping);
             break;
           case TIME_STEPPING_HANCOCK:
             error_flag = dUdt_Hancock(Soln_ptr,
				       Input_Parameters.Number_of_Cells,
				       dtime,
				       Input_Parameters.CFL_Number,
				       Input_Parameters.i_Limiter, 
				       Input_Parameters.i_Flux_Function,
				       Input_Parameters.Local_Time_Stepping);
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
                 cout << "\nPDES++ ERROR: Euler1D solution error.\n\n";
             } else {
                 cout << "\n PDES++ ERROR: Euler1D solution error.";
                 cout << "\n\nPDES++: Execution terminated.\n";
             } /* endif */
             return (error_flag);
         } /* endif */
      
         /* Update time and time step counter. */
         number_of_time_steps = number_of_time_steps + 1;
         time = time + Input_Parameters.CFL_Number*dtime;

     } /* endwhile */
     if (! batch_flag) cout << "\n\n Euler1D computations complete.\n";
  } /* endif */

  /**************************************************  
   * Reconstruct the final solution (high-order or  *
   * limited linear reconstruction)                 *
   *************************************************/
  if (Input_Parameters.Verbose()){
    cout << '\n'
	 << " Reconstruct the final solution. "
	 << "\n";
  }
  if(Input_Parameters.i_Reconstruction == RECONSTRUCTION_HIGH_ORDER){
    HighOrderSolutionReconstructionOverDomain(Soln_ptr,Input_Parameters,&Euler1D_UniformMesh::CellHighOrder);
  } else {
    LimitedLinearReconstructionOverDomain(Soln_ptr,Input_Parameters);
  }

  /********************************************************  
   * Write 1D Euler solution to output file.              *
   ********************************************************/

  while (1) {
     Get_Next_Input_Control_Parameter(Input_Parameters);
     command_flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
     line_number = Input_Parameters.Line_Number;
     if (command_flag == EXECUTE_CODE) {
         // Deallocate memory for 1D Euler equation solution.
         if (! batch_flag) cout << "\n Deallocating Euler1D solution variables.";
         Soln_ptr=Deallocate(Soln_ptr);
         // Execute new calculation.
         if (! batch_flag)  {
            cout << "\n\n Starting a new calculation.";
            cout << Input_Parameters << "\n";
         } /* endif */
         goto execute_new_calculation;

     } else if (command_flag == TERMINATE_CODE) {
         // Deallocate memory for 1D Euler equation solution.
         if (! batch_flag) cout << "\n Deallocating Euler1D solution variables.";
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
           cout << "\n Writing Euler1D solution to output data file `"
	        << Output_File_Name_ptr << "'.";
         Output_File.open(Output_File_Name_ptr, ios::out);
         if (Output_File.fail()) {
            if (batch_flag) {
               cout << "\nPDES++ ERROR: Unable to open Euler1D output data file.\n\n";
            } else {
               cout << "\n PDES++ ERROR: Unable to open Euler1D output data file.";
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
	   	           Input_Parameters,
		           number_of_time_steps,
		           time,
		           Output_File);
         } /* endif */

         Output_File.close();

         if (Input_Parameters.i_Output_Format == IO_GNUPLOT) {
            Gnuplot_File_Name_ptr = Input_Parameters.Gnuplot_File_Name;
            Gnuplot_File.open(Gnuplot_File_Name_ptr, ios::out);
            if (Gnuplot_File.fail()) {
               if (batch_flag) {
                  cout << "\nPDES++ ERROR: Unable to open Euler1D Gnuplot macro file.\n\n";
               } else {
                  cout << "\n PDES++ ERROR: Unable to open Euler1D Gnuplot macro file.";
                  cout << "\n\nPDES++: Execution terminated.\n";
               } /* endif */
               return (2);
            } /* endif */

            Gnuplot_File << "set title \"PDES++: 1D Euler Solution\"\n"
                         << "set xlabel \"x\"\n"
	                 << "plot \"" << Output_File_Name_ptr << "\""
                         << " using 1:2 \"%*lf%lf%*lf%lf%*lf%*lf%*lf\" \\\n"
	                 << "     title \"density\" with points\n"
                         << "pause -1  \"Hit return to continue\"\n"
                         << "plot \"" << Output_File_Name_ptr << "\""
                         << " using 1:2 \"%*lf%lf%*lf%*lf%lf%*lf%*lf\" \\\n"
	                 << "     title \"velocity\" with points\n"
                         << "pause -1  \"Hit return to continue\"\n"
                         << "plot \"" << Output_File_Name_ptr << "\""
                         << " using 1:2 \"%*lf%lf%*lf%*lf%*lf%lf%*lf\" \\\n"
	                 << "     title \"pressure\" with points\n"
                         << "pause -1  \"Hit return to continue\"\n"
                         << "plot \"" << Output_File_Name_ptr << "\""
                         << " using 1:2 \"%*lf%lf%*lf%*lf%*lf%*lf%lf\" \\\n"
	                 << "     title \"temperature\" with points\n"
                         << "pause -1  \"Hit return to continue\"\n";

            Gnuplot_File.close();
         } /* endif */

     } else if (command_flag == INVALID_INPUT_CODE ||
                command_flag == INVALID_INPUT_VALUE) {
         line_number = -line_number;
         if (batch_flag) {
             cout << "\nEuler1D ERROR: Error reading Euler1D data at line #"
	          << -line_number << " of input data file.\n\n";
         } else {
             cout << "\nEuler1D ERROR: Error reading Euler1D data at line #"
	          << -line_number  << " of input data file.";
             cout << "\n\nPDES++: Execution terminated.\n";
         } /* end if */
         return (line_number);
     } /* endif */
  } /* endwhile */

  /********************************************************  
   * End of all Euler1DSolver computations and I/O.       *
   ********************************************************/

  return (0);
  
}
