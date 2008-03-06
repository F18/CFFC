/*!\file Levermore1DSolvers.cc
  \brief 1D Levermore Equation Solvers. */

/* Include 1D Levermore solution header file. */

#ifndef _LEVERMORE1D_INCLUDED
#include "Levermore1D.h"
#endif // _LEVERMORE1D_INCLUDED

/********************************************************
 * Routine: Levermore1DSolver                               *
 *                                                      *
 * Computes solutions to 1D Levermore equations.            *
 *                                                      *
 ********************************************************/
int Levermore1DSolver(char *Input_File_Name_ptr,
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

  Levermore1D_UniformMesh *Soln_ptr = NULL;

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
        cout << "\n Levermore1D ERROR: Unable to open Levermore1D input data file.\n\n";
     } else {
        cout << "\n Levermore1D ERROR: Unable to open Levermore1D input data file.";
        cout << "\n\nLevermore1D: Execution terminated.\n";
     } /* endif */
     return (-1);
  } /* endif */

  if (! batch_flag) {
    cout << "\n Reading Levermore1D input data file `"
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
             cout << "\n Levermore1D ERROR: Error reading Levermore1D data at line #"
	          << -line_number << " of input data file.\n\n";
         } else {
             cout << "\n Levermore1D ERROR: Error reading Levermore1D data at line #"
	          << -line_number  << " of input data file.";
             cout << "\n\nLevermore1D: Execution terminated.\n";
         } /* end if */
         return (line_number);
     } /* endif */
  } /* endwhile */

  if (! batch_flag)  {
     cout << Input_Parameters << "\n";
     cout.flush();
  } /* endif */

  /*********************************************************
   * Create mesh and allocate Levermore1D solution variables.  *
   *********************************************************/

  execute_new_calculation: ;

  /* Set number of moments used in calculation. */
  if (! batch_flag){
    cout << "\n Setting number of moments to " << Input_Parameters.number_of_moments << ".";
    cout.flush();
  }

  Levermore1D_Vector::set_length(Input_Parameters.number_of_moments);

  /* Allocate memory for 1D Levermore equation solution on
     uniform mesh. */

  if (! batch_flag){
    cout << "\n Creating memory for Levermore1D solution variables.";
    cout.flush();
  }
  Soln_ptr=Allocate(Soln_ptr,
                    Input_Parameters);

  if (Soln_ptr == NULL){
    cout << "\n Levermore1DSolvers::Allocate() Error! Probably not enough memory!";
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
   * Initialize Levermore1D solution variables.               *
   ********************************************************/

  /* Set the initial time level. */

  time = ZERO;
  number_of_time_steps = 0;

  /* Initialize the conserved and primitive state
     solution variables. */

  if (! batch_flag){
    cout << "\n Prescribing Levermore1D initial data.";
    cout.flush();
  }
  ICs(Soln_ptr,
      "ZB",
      Input_Parameters.i_ICs,
      Input_Parameters.Number_of_Cells,
      Input_Parameters);

  /********************************************************
   * Solve conservation form of 1D Levermore equations for    *
   * specified IBVP or BVP on uniform mesh.               *
   ********************************************************/

  continue_existing_calculation: ;

  if ((Input_Parameters.Local_Time_Stepping &&
       Input_Parameters.Maximum_Number_of_Time_Steps > 0) ||
      (!Input_Parameters.Local_Time_Stepping &&
       Input_Parameters.Time_Max > time)) {
    if (! batch_flag){
      cout << "\n\n Beginning Levermore1D computations.\n\n";
      cout.flush();
    }
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
             error_flag = dUdt_explicitEuler_upwind(Soln_ptr,
						    Input_Parameters.Number_of_Cells,
						    dtime,
						    Input_Parameters.CFL_Number,
						    Input_Parameters.i_Flux_Function,
						    Input_Parameters.Local_Time_Stepping);
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
//           case TIME_STEPPING_HANCOCK:
//             error_flag = dUdt_Hancock(Soln_ptr,
//				       Input_Parameters.Number_of_Cells,
//				       dtime,
//				       Input_Parameters.CFL_Number,
//				       Input_Parameters.i_Limiter,
//				       Input_Parameters.i_Flux_Function,
//				       Input_Parameters.Local_Time_Stepping);
//             break;
           default:
	     cout << "Invalid time marching chosen!" << endl;
	     return (1);
             break;
         } /* endswitch */

         if (error_flag) {
             if (batch_flag) {
                 cout << "\nCFFC ERROR: Levermore1D solution error.\n\n";
		 cout.flush();
             } else {
                 cout << "\n CFFC ERROR: Levermore1D solution error.";
                 cout << "\n\nCFFC: Execution terminated.\n";
		 cout.flush();
             } /* endif */
             return (error_flag);
         } /* endif */

         /* Update time and time step counter. */
         number_of_time_steps = number_of_time_steps + 1;
         time = time + Input_Parameters.CFL_Number*dtime;
	 cout << "[" << time << "]";cout.flush();
     } /* endwhile */
     if (! batch_flag){
       cout << "\n\n Levermore1D computations complete.\n";
       cout.flush();
     }
  } /* endif */

  /**************************************************
   * Reconstruct the final solution (high-order or  *
   * limited linear reconstruction)                 *
   *************************************************/
//  if (Input_Parameters.Verbose()){
//    cout << '\n'
//	 << " Reconstruct the final solution. "
//	 << "\n";
//    cout.flush();
//  }
//  if(Input_Parameters.i_Reconstruction == RECONSTRUCTION_HIGH_ORDER){
//    HighOrderSolutionReconstructionOverDomain(Soln_ptr,Input_Parameters,&Levermore1D_UniformMesh::CellHighOrder);
//  } else {
//    LimitedLinearReconstructionOverDomain(Soln_ptr,Input_Parameters);
//  }

  /********************************************************
   * Write 1D Levermore solution to output file.              *
   ********************************************************/

  while (1) {
     Get_Next_Input_Control_Parameter(Input_Parameters);
     command_flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
     line_number = Input_Parameters.Line_Number;
     if (command_flag == EXECUTE_CODE) {
         // Deallocate memory for 1D Levermore equation solution.
         if (! batch_flag) cout << "\n Deallocating Levermore1D solution variables.";
         Soln_ptr=Deallocate(Soln_ptr);
         // Execute new calculation.
         if (! batch_flag)  {
            cout << "\n\n Starting a new calculation.";
            cout << Input_Parameters << "\n";
         } /* endif */
//	 // Reset internal flags AccuracyAssessment1D
//	 AccuracyAssessment1D::ResetForNewCalculation();
         goto execute_new_calculation;

     } else if (command_flag == TERMINATE_CODE) {
         // Deallocate memory for 1D Levermore equation solution.
         if (! batch_flag) cout << "\n Deallocating Levermore1D solution variables.";
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
//	 // Reset internal flags AccuracyAssessment1D
//	 AccuracyAssessment1D::ResetForNewCalculation();
         goto continue_existing_calculation;

     } else if (command_flag == WRITE_OUTPUT_CODE) {
         // Output solution data.
         Output_File_Name_ptr = Input_Parameters.Output_File_Name;
         if (! batch_flag)
           cout << "\n Writing Levermore1D solution to output data file `"
	        << Output_File_Name_ptr << "'.";
         Output_File.open(Output_File_Name_ptr, ios::out);
         if (Output_File.fail()) {
            if (batch_flag) {
               cout << "\nCFFC ERROR: Unable to open Levermore1D output data file.\n\n";
            } else {
               cout << "\n CFFC ERROR: Unable to open Levermore1D output data file.";
               cout << "\n\nCFFC: Execution terminated.\n";
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
                  cout << "\nCFFC ERROR: Unable to open Levermore1D Gnuplot macro file.\n\n";
               } else {
                  cout << "\n CFFC ERROR: Unable to open Levermore1D Gnuplot macro file.";
                  cout << "\n\nCFFC: Execution terminated.\n";
               } /* endif */
               return (2);
            } /* endif */

            Gnuplot_File << "set title \"CFFC: 1D Levermore Solution\"\n"
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

//     } else if (command_flag == WRITE_OUTPUT_ACCURACY_CODE) {
//       // Output error norms to tecplot file
//       AccuracyAssessment1D::OutputErrorNormsTecplot(Soln_ptr,Input_Parameters);
//
//    } else if (command_flag == WRITE_NORM_ON_SCREEN) {
//       // Print solution
//       AccuracyAssessment1D::PrintErrorNorms(Soln_ptr,Input_Parameters,cout);

     } else if (command_flag == INVALID_INPUT_CODE ||
                command_flag == INVALID_INPUT_VALUE) {
         line_number = -line_number;
         if (batch_flag) {
             cout << "\nLevermore1D ERROR: Error reading Levermore1D data at line #"
	          << -line_number << " of input data file.\n\n";
         } else {
             cout << "\nLevermore1D ERROR: Error reading Levermore1D data at line #"
	          << -line_number  << " of input data file.";
             cout << "\n\nCFFC: Execution terminated.\n";
         } /* end if */
         return (line_number);
     } /* endif */
  } /* endwhile */

  /********************************************************
   * End of all Levermore1DSolver computations and I/O.   *
   ********************************************************/

  return (0);

}
