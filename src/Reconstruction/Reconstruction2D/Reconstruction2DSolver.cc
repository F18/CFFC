/* Reconstruction2DSolver.cc:  Source file defining 
   2D ReconstructionSolver function */

/* Include header files. */
#include "CFD/CFD.h"
#include "Reconstruction/ComputationalDomain.h"

/********************************************************
 * Routine: Reconstruction2DSolver                      *
 *                                                      *
 * Computes numerical reconstruction to 2D functions    *
 *                                                      *
 ********************************************************/
int Reconstruction2DSolver(char *Input_File_Name_ptr){

  // Reconstruction2D input variables and parameters:
  Reconstruct2D_Input_Parameters Input_Parameters;

  // Other local solution variables.
  int command_flag, error_flag, line_number;
  bool AccuracyAssessed_flag = false;
  bool Title_Error_Norms = true;
  int CurrentOrderOfReconstruction(-1);

  /**********************************************************
   *               Pre-processing Step                      *
   *********************************************************/
  
  /********************************************************  
   * Set default values for the input solution parameters *
   * and then read user specified input values from the   *
   * specified input parameter file.                      *
   ********************************************************/

  // Read and parse the input file.
  cout << "\n Reading Reconstruction2D input data file `"
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
   * - Reconstruction2D_CartesianMesh                            *
   * - Reconstruction2D_QuadrilateralMesh                        *
   **************************************************************/

  // Solution variables.
  ComputationalDomain<TwoD,Cell2D_Quad,double> SolnBlkDouble;           //

  // Mesh variables
  Grid2D_Quad_Block   Grid;
  Grid2D_Quad_Block   **MeshBlk(NULL); // the mesh

 execute_new_calculation: ;

  if (Input_Parameters.TestF == NULL){ // the grid and the solution must be read from a file
    std::cout << "\n Read grid and solution data from file: " << Input_Parameters.Input_File_Name << "\n";

    // Set the Computation Domain
    SolnBlkDouble.SetDomain(Input_Parameters);

    // Reconstruct the solution
    SolnBlkDouble.ReconstructZonalSolution(Input_Parameters);


  } else {
    std::cout << "\n Create Grid:\n";
    // Set Grid
    MeshBlk = Multi_Block_Grid(MeshBlk,Input_Parameters);
    if (MeshBlk != NULL) {
      if (Check_Multi_Block_Grid(MeshBlk,
				 Input_Parameters.Number_of_Blocks_Idir,
				 Input_Parameters.Number_of_Blocks_Jdir)) {
	std::cout << "Check_Multi_Block_Grid ERROR\n";
      } 
    } /* endif */

    // copy the mesh to the computational domain
    Copy_Quad_Block(Grid,MeshBlk[0][0]);
  
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
   * Write 2D reconstruction solution to output and  *
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
      if (MeshBlk != NULL){
	  // Deallocate memory for 2D reconstruction solution.
	  cout << "\n Deallocating Reconstruction2D "
	       << "solution variables.";
	  MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
						Input_Parameters.Number_of_Blocks_Idir, 
						Input_Parameters.Number_of_Blocks_Jdir);
      }

      // Output input parameters for new caluculation.
      cout << "\n\n Starting a new calculation.";
      cout << Input_Parameters << "\n";
      cout.flush();
      
      // set accuracy flag
      AccuracyAssessed_flag = false;

      // Execute new calculation.
      goto execute_new_calculation;
      
    } else if (command_flag == TERMINATE_CODE) {
      // Deallocate memory for 2D reconstruction solution.
      cout << "\n\n Deallocating Reconstruction2D"
	   <<" solution variables.";

      MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
					    Input_Parameters.Number_of_Blocks_Idir, 
					    Input_Parameters.Number_of_Blocks_Jdir);

      Grid.deallocate();

      // Close input data file.
      cout << "\n\n Closing Reconstruction2D input" 
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
	cout << "\n Writing Reconstruction2D" 
	     << " solution to output data file(s): ";
	Output_Tecplot(SolnBlkDouble,Input_Parameters);
	//      }
      
    } else if (command_flag == WRITE_UNIFIED_OUTPUT_CODE) {
      // Output solution data.
      //      if (MeshBlk != NULL){
	cout << "\n Writing Reconstruction2D" 
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
	cout << "\n Writing Reconstruction2D solution at nodes to"
	     << " output data file(s): ";
	Output_Solution_Nodes_Tecplot(SolnBlkDouble, Input_Parameters);
	//      }

    } else if (command_flag == WRITE_OUTPUT_FULL_GRID_NODES_CODE) {
      // Output mesh node data file.
      //      if (MeshBlk != NULL){
	cout << "\n Writing Reconstruction2D solution at nodes to"
	     << " output data file(s): ";
	Output_Full_Solution_Nodes_Tecplot(SolnBlkDouble, Input_Parameters);
	//      }

    } else if (command_flag == INVALID_INPUT_CODE ||
	       command_flag == INVALID_INPUT_VALUE) {
      line_number = -line_number;
      cout << "\n Reconstruct2D ERROR: Error reading Reconstruct2D" 
	   << " data at line #" << -line_number << " of input data file.";
      cout.flush();
      return (line_number);
    } // endif
    
  } // endwhile
  
  /**********************************************************  
   * End of all Reconstruction2DSolver computations and I/O. *
   **********************************************************/
  
  return (0);
}
