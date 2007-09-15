/* Reconstruction1DSolverFramework.cc:  Source file defining 
   1D ReconstructionSolver function */

/* Include header files. */

#include "CFD/CFD.h"
#include "Reconstruction/ComputationalDomain.h"
#include "Common/SourceRevisionData.h"

/********************************************************
 * Routine: Reconstruction1DSolver                      *
 *                                                      *
 * Computes numerical reconstruction to 1D functions    *
 *                                                      *
 ********************************************************/
int Reconstruction1DSolver(char *Input_File_Name_ptr){

  // Reconstruction1D input variables and parameters:
  Reconstruct1D_Input_Parameters Input_Parameters;

  // Other local solution variables.
  int command_flag, error_flag, line_number;
  bool AccuracyAssessed_flag = false;
  bool Title_Error_Norms = true;
  int CurrentOrderOfReconstruction(-1);

  /**********************************************************
   *               Pre-processing Step                      *
   *********************************************************/

  cout << '\n'
       << "****************************************************" << "\n"
       << "*\t Running " << SourceCode::ProgramName() << " for 1D " << "\n"
       << "****************************************************" << "\n"
       << " Executable information:" << "\n"
       << "    --> repository version: rev. " << SourceCode::LastCommitted_Revision() << '\n'
       << "    --> committed on: " << SourceCode::LastCommitted_Date() << '\n'
       << "    --> committed by: " << SourceCode::LastCommitted_Author() << '\n'
       << "    --> local version: rev. " << SourceCode::RevisionAtCompileTime() << '\n'
       << "    --> compiled on: " << SourceCode::TimeAtCompilation() << '\n'
       << "####################################################" << "\n";
  cout.flush();
  
  /********************************************************  
   * Set default values for the input solution parameters *
   * and then read user specified input values from the   *
   * specified input parameter file.                      *
   ********************************************************/

  // Read and parse the input file.
  cout << "\n Reading Reconstruction1D input data file `"
       << Input_File_Name_ptr << "'.";

  error_flag = Process_Input_Control_Parameter_File(Input_Parameters,
	 					    Input_File_Name_ptr,
		 				    command_flag);
  if (error_flag == 0) {
    cout << Input_Parameters << "\n";
    cout.flush();
  } else {
    // Error in the input file
    exit(1);
  }//endif 
  
  // update the CurrentOrderOfReconstruction
  CurrentOrderOfReconstruction = Input_Parameters.RecOrder();

  /***************************************************************  
   *           Set-Up the Computational Domain                   *
   **************************************************************/

  // Solution variables.
  ComputationalDomain<OneD,Cell1D_NonUniform,double> SolnBlkDouble;
  //ComputationalDomain<OneD,Cell1D_NonUniform,Gaussian1D_pState> SolnBlkGauss;

 execute_new_calculation: ;

  SetUpDomain(SolnBlkDouble,Input_Parameters);
  
  
  /**********************************************************
   *               Solver Step                              *
   *********************************************************/

  SolveTask(SolnBlkDouble,Input_Parameters);

  /***************************************************
   * Solution calculations complete.                 *
   * Write 1D reconstruction solution to output and  *
   * restart files as required, reset solution       *
   * parameters, and run other cases as specified by *
   * input parameters                                *
   **************************************************/

  /**********************************************************
   *               Post-processing Step                     *
   *********************************************************/
  
  // postprocess_current_calculation: ;

  while (1) {
    Get_Next_Input_Control_Parameter(Input_Parameters);
    command_flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
    line_number = Input_Parameters.Line_Number;

    if (command_flag == EXECUTE_CODE) {
      // Output input parameters for new caluculation.
      cout << "\n \n Starting a new calculation.";
      cout << Input_Parameters << "\n";
      cout.flush();
      
      // set accuracy flag
      AccuracyAssessed_flag = false;

      // Execute new calculation.
      goto execute_new_calculation;
      
    } else if (command_flag == TERMINATE_CODE) {
      // Close input data file.
      cout << "\n \n Closing Reconstruction1D input" 
	   << " data file.";
      Close_Input_File(Input_Parameters);
      // Terminate calculation.
      return (0);
      
    } else if (command_flag == 2) {
      // The reconstruction order has been read
      if (CurrentOrderOfReconstruction != Input_Parameters.RecOrder()){
	// set the title flag
	Title_Error_Norms = true;

	// update the CurrentOrderOfReconstruction
	CurrentOrderOfReconstruction = Input_Parameters.RecOrder();
      }

    } else if (command_flag == WRITE_OUTPUT_CODE) {
      // Output solution data.
      if (!SolnBlkDouble.null()){
	cout << "\n Writing Reconstruction1D" 
	     << " solution to output data file(s): ";
	Output_Tecplot(SolnBlkDouble,Input_Parameters);
      }
      
    } else if (command_flag == WRITE_UNIFIED_OUTPUT_CODE) {
      // Output solution data.
      if (!SolnBlkDouble.null()){
	cout << "\n Writing Reconstruction1D" 
	     << " unified solution to output data file(s): ";
	Output_Tecplot(SolnBlkDouble,Input_Parameters,false);
      }

    } else if (command_flag == WRITE_RECONSTRUCTED_FUNCTIONS) {
      // Output solution data.
      if (!SolnBlkDouble.null()){
	cout << "\n Writing Reconstruction1D" 
	     << " solution in every stencil to output data file(s): ";
	Output_Tecplot_Stencil_Reconstruction(SolnBlkDouble,Input_Parameters);
      }

    } else if (command_flag == WRITE_OUTPUT_FUNCTION_GRAPH_CODE) {
      // Output the data defining the reconstructed function
      // using the its mathematical definition
      // on the specified domain.
      cout << "\n Writing exact function data" 
	   << " to output data file(s): ";

      Output_Function_Graph_Tecplot(SolnBlkDouble,Input_Parameters);

      if (error_flag) {
	cout << "\n  Reconstruct1D ERROR: Unable to open Reconstruct1D function"
	     << " definition output data file(s)" << ".\n";
	cout.flush();
      } // endif
      
    } else if (command_flag == WRITE_OUTPUT_ACCURACY_CODE) {
      if(AccuracyAssessed_flag){
      	if (!SolnBlkDouble.null()){
	  // output to file
	  Output_Error_Norms_Tecplot(SolnBlkDouble,Input_Parameters,Title_Error_Norms);
	}
      }
      else {
	// Analyze reconstruction to assess the accuracy
	std::cout << "\nAssess the accuracy:";

	SolnBlkDouble.AssessReconstructionAccuracy(Input_Parameters);
	
	std::cout << "\nAccuracy assessed";
	// change flag
	AccuracyAssessed_flag = true;
	
	// output to file
	cout << "\n Writing Reconstruction1D" 
	     << " error norms to output data file(s): ";
	Output_Error_Norms_Tecplot(SolnBlkDouble,Input_Parameters,Title_Error_Norms);
      }
	 
    } else if (command_flag == WRITE_NORM_ON_SCREEN) {
      // Print solution
      if(AccuracyAssessed_flag){
	std::cout << " The accuracy norms of the reconstruction are:\n";
	SolnBlkDouble.PrintErrorNorms();
      }
      else{
	std::cout << "\nAssess the accuracy";
	// Analyze reconstruction to assess the accuracy
	SolnBlkDouble.AssessReconstructionAccuracy(Input_Parameters);
	
	std::cout << "\nAccuracy assessed";
	// change flag
	AccuracyAssessed_flag = true;
	// print solution
	std::cout << "The accuracy norms of the reconstruction are:\n";
	SolnBlkDouble.PrintErrorNorms();
      }
      
    } else if (command_flag == WRITE_OUTPUT_GRID_NODES_CODE) {
      // Output mesh node data file.
      if (!SolnBlkDouble.null()){
	cout << "\n Writing Reconstruction1D solution at nodes to"
	     << " output data file(s).";
	Output_Solution_Nodes_Tecplot(SolnBlkDouble, Input_Parameters);
      }

    } else if (command_flag == WRITE_OUTPUT_FULL_GRID_NODES_CODE) {
      // Output mesh node data file.
      if (!SolnBlkDouble.null()){
	cout << "\n Writing Reconstruction1D solution at nodes to"
	     << " output data file(s): ";
	Output_Full_Solution_Nodes_Tecplot(SolnBlkDouble, Input_Parameters);
      }
      
    } else if (command_flag == INVALID_INPUT_CODE ||
	       command_flag == INVALID_INPUT_VALUE) {
      line_number = -line_number;
      cout << "\n Reconstruct1D ERROR: Error reading Reconstruct1D" 
	   << " data at line #" << -line_number << " of input data file.\n";
      cout.flush();
      return (line_number);
    } // endif
    
  } // endwhile
  
  /**********************************************************  
   * End of all Reconstruction1DSolver computations and I/O. *
   **********************************************************/
  
  return (0);
}
