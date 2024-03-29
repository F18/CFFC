/* Reconstruction3DSolver.cc:  Source file defining 
   3D ReconstructionSolver function */

/* Include header files. */
#include "CFD/CFD.h"
#include "Reconstruction/ComputationalDomain.h"
#include "Common/SourceRevisionData.h"

/********************************************************
 * Routine: Reconstruction3DSolver                      *
 *                                                      *
 * Computes numerical reconstruction to 3D functions    *
 *                                                      *
 ********************************************************/
int Reconstruction3DSolver(char *Input_File_Name_ptr){

  // Reconstruction3D input variables and parameters:
  Reconstruct3D_Input_Parameters Input_Parameters;

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
       << "*\t Running " << SourceCode::ProgramName() << " for 3D " << "\n"
       << "****************************************************" << "\n"
       << " Executable information:" << "\n"
       << "    --> GIT repository version: hash = " << SourceCode::LastCommitted_HashID() << '\n'
       << "    --> committed on: " << SourceCode::LastCommitted_Date() << '\n'
       << "    --> committed by: " << SourceCode::LastCommitted_Author() << '\n'
       << "    --> compiled on: " << SourceCode::TimeAtCompilation() << '\n'
       << "####################################################" << "\n";
  cout.flush();
  
  /********************************************************  
   * Set default values for the input solution parameters *
   * and then read user specified input values from the   *
   * specified input parameter file.                      *
   ********************************************************/

  // Read and parse the input file.
  cout << "\n Reading Reconstruction3D input data file `"
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
  
  /***************************************************************  
   * Create the mesh variable and call for the proper solver:    *
   * - Reconstruction3D_HexahedralMesh                           *
   **************************************************************/

  // Solution variables.
  ComputationalDomain<ThreeD,Cell3D_Hexa,double> SolnBlkDouble;           //

  // Mesh variables
  Grid3D_Hexa_Block   Grid;

 execute_new_calculation: ;

  if (Input_Parameters.TestF == NULL){ // the grid and the solution must be read from a file
    std::cout << "\n Read grid and solution data from file: " << Input_Parameters.Input_File_Name << "\n";

    // Set the Computation Domain
    SolnBlkDouble.SetDomain(Input_Parameters);

    // Reconstruct the solution
    SolnBlkDouble.ReconstructSolution(Input_Parameters);


  } 
   else {
    std::cout << "\n Create Grid:\n";
    // Set Grid
    Grid.Create_Block(Input_Parameters.Box_Length,
		      Input_Parameters.Box_Width,
		      Input_Parameters.Box_Height,
		      ZERO,    // x-orig
		      ZERO,    // y-orig
		      ZERO,    // z-orig
		      ZERO,    // alpha
		      ZERO,    // beta
		      ZERO,    // gamma
		      BC_NONE, // top
		      BC_NONE, // bottom
		      BC_NONE, // north
		      BC_NONE, // south
		      BC_NONE, // west
		      BC_NONE, // east
		      Input_Parameters.NCells_Idir,
		      Input_Parameters.NCells_Jdir,
		      Input_Parameters.NCells_Kdir,
		      Input_Parameters.Nghost_Cells);
    
    
    std::cout << " Set the Computational Domain:\n";
    // Set the Computation Domain
    SolnBlkDouble.SetDomain(Grid,Input_Parameters);
  
    std::cout << " Compute the intial solution:\n";
    // Compute the initial solution
    SolnBlkDouble.SetInitialData(Input_Parameters);

    /**********************************************************
     *               Solver Step                              *
     *********************************************************/
  
    SolveTask(SolnBlkDouble,Input_Parameters);

  }

  /***************************************************
   * Solution calculations complete.                 *
   * Write 3D reconstruction solution to output and  *
   * restart files as required, reset solution       *
   * parameters, and run other cases as specified by *
   * input parameters                                *
   **************************************************/
  
  // postprocess_current_calculation: ;

  while (1) {
    Get_Next_Input_Control_Parameter(Input_Parameters);
    command_flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
    line_number = Input_Parameters.Line_Number;

    if (command_flag == EXECUTE_CODE) {

      // Output input parameters for new caluculation.
      cout << "\n\n Starting a new calculation.";
      cout << Input_Parameters << "\n";
      cout.flush();
      
      // set accuracy flag
      AccuracyAssessed_flag = false;

      // Execute new calculation.
      goto execute_new_calculation;
      
    } else if (command_flag == TERMINATE_CODE) {
      // Deallocate memory for 3D reconstruction solution.
      cout << "\n\n Deallocating Reconstruction3D"
	   <<" solution variables.";

      Grid.deallocate();

      // Close input data file.
      cout << "\n\n Closing Reconstruction3D input" 
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
      //      if (MeshBlk != NULL){
	cout << "\n Writing Reconstruction3D" 
	     << " solution to output data file(s): ";
	Output_Tecplot(SolnBlkDouble,Input_Parameters);
	//      }
      
    } else if (command_flag == WRITE_UNIFIED_OUTPUT_CODE) {
      // Output solution data.
      //      if (MeshBlk != NULL){
	cout << "\n Writing Reconstruction3D" 
	     << " solution to output data file(s): ";
	Output_Tecplot(SolnBlkDouble,Input_Parameters,false);
	//      }

    } else if (command_flag == WRITE_OUTPUT_FUNCTION_GRAPH_CODE) {
      // Output the data defining the reconstructed function
      // using the its mathematical definition
      // on the specified domain.
      cout << "\n Writing exact function data" 
	   << " to output data file(s): ";

      Output_Function_Graph_Tecplot(SolnBlkDouble,Input_Parameters);
      
    } else if (command_flag == WRITE_OUTPUT_ACCURACY_CODE) {
      if(AccuracyAssessed_flag){
	//      	if (MeshBlk != NULL){
	  // output to file
	  Output_Error_Norms_Tecplot(SolnBlkDouble,Input_Parameters,Title_Error_Norms);
	  //	}
      }
      else {
	// Analyze reconstruction to assess the accuracy
	std::cout << "\nAssess the accuracy:\n";

	SolnBlkDouble.AssessReconstructionAccuracy(Input_Parameters);
	
	std::cout << "\nAccuracy assessed\n";
	// change flag
	AccuracyAssessed_flag = true;
	
	// output to file
	Output_Error_Norms_Tecplot(SolnBlkDouble,Input_Parameters,Title_Error_Norms);
      }
	 
    } else if (command_flag == WRITE_NORM_ON_SCREEN) {
      // Print solution
      if(AccuracyAssessed_flag){
	std::cout << " The accuracy norms of the reconstruction are:\n";
	SolnBlkDouble.PrintErrorNorms();
      }
      else{
	std::cout << "\nAssess the accuracy\n";
	// Analyze reconstruction to assess the accuracy
	SolnBlkDouble.AssessReconstructionAccuracy(Input_Parameters);
	
	std::cout << "\nAccuracy assessed\n";
	// change flag
	AccuracyAssessed_flag = true;
	// print solution
	std::cout << "The accuracy norms of the reconstruction are:\n";
	SolnBlkDouble.PrintErrorNorms();
      }

    } else if (command_flag == WRITE_OUTPUT_GRID_CODE) {
      // Output mesh node data file.
      //      if (MeshBlk != NULL){
	Output_Mesh_Nodes_Tecplot(SolnBlkDouble, Input_Parameters);
	//      }
      
    } else if (command_flag == WRITE_OUTPUT_GRID_NODES_CODE) {
      // Output node solution data file.
      //      if (MeshBlk != NULL){
	cout << "\n Writing grid at cells to"
	     << " output data file(s): ";
	Output_Mesh_Nodes_Tecplot(SolnBlkDouble, Input_Parameters);
	//      }
    } else if (command_flag == WRITE_OUTPUT_GRID_CELLS_CODE) {
      // Output node solution data file.
      //      if (MeshBlk != NULL){
	cout << "\n Writing Reconstruction3D solution at nodes to"
	     << " output data file(s): ";
	Output_Mesh_Cells_Tecplot(SolnBlkDouble, Input_Parameters);
	//      }

    } else if (command_flag == WRITE_OUTPUT_FULL_GRID_NODES_CODE) {
      // Output mesh node data file.
      //      if (MeshBlk != NULL){
	cout << "\n Writing Reconstruction3D solution at nodes to"
	     << " output data file(s): ";
	Output_Full_Solution_Nodes_Tecplot(SolnBlkDouble, Input_Parameters);
	//      }
    } else if (command_flag == INVALID_INPUT_CODE ||
	       command_flag == INVALID_INPUT_VALUE) {
      line_number = -line_number;
      cout << "\n Reconstruct3D ERROR: Error reading Reconstruct3D" 
	   << " data at line #" << -line_number << " of input data file.";
      cout.flush();
      return (line_number);
    } // endif
    
  } // endwhile

  /**********************************************************  
   * End of all Reconstruction3DSolver computations and I/O. *
   **********************************************************/
  
  return (0);
}
