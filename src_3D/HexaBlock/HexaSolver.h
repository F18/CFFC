#ifndef _HEXA_SOLVER_INCLUDED
#define _HEXA_SOLVER_INCLUDED

/* Include required CFFC header files. */

#ifndef _HEXA_MULTIBLOCK_INCLUDED
#include "HexaMultiBlock.h"
#endif //_HEXA_MULTIBLOCK_INCLUDED

#ifndef _AMR_INCLUDED
#include "../AMR/AMR.h"
#endif // _AMR_INCLUDED

#ifndef _ADAPTIVEBLOCK3D_MESSAGEPASSING_INCLUDED
#include "../AMR/AdaptiveBlock3D_MessagePassing.h"
#endif // _ADAPTIVEBLOCK3D_MESSAGEPASSING_INCLUDED

#ifndef _HEXA_SOLVER_CLASSES_INCLUDED
#include  "HexaSolverClasses.h"
#endif  //_HEXA_SOLVER_CLASSES_INCLUDED

#ifndef _EXPLICIT_FILTER_COMMANDS_INCLUDED
#include "../ExplicitFilters/Explicit_Filter_Commands.h"
#endif // _EXPLICIT_FILTER_COMMANDS_INCLUDED

#ifndef _HEXA_PRE_PROCESSING_INCLUDED
#include "HexaPreProcessing.h"
#endif //_HEXA_PRE_PROCESSING_INCLUDED

#ifndef _HEXA_POST_PROCESSING_INCLUDED
#include "HexaPostProcessing.h"
#endif //_HEXA_POST_PROCESSING_INCLUDED

#ifndef _HEXA_EXPLICIT_SOLVER
#include "HexaExplicitSolver.h"
#endif //_HEXA_EXPLICIT_SOLVER

#ifndef _NKS_INCLUDED
#include "../NewtonKrylovSchwarz/NKS.h"
#endif //_NKS_INCLUDED


/********************************************************
 * Routine: HexaSolver                                  *
 *                                                      *
 * Computes solutions to the governing PDEs on          *
 * multi-block AMR mesh composed of body-fitted         *
 * hexahedral solution blocks.                          *
 *                                                      *
 ********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int HexaSolver(char *Input_File_Name_ptr, int batch_flag) {
   
  int error_flag(0);

  // Non Solution Specific
  HexaSolver_Data                                     Data(batch_flag);  
  // Solution Specific 
  HexaSolver_Solution_Data<SOLN_pSTATE, SOLN_cSTATE>  Solution_Data;   

  /******************* INPUT PARAMETERS  **********************************
    Set default values for the input solution parameters and then read user 
    specified input values from the specified input parameter file.               
  *************************************************************************/    
  error_flag = Solution_Data.Get_Input_Parameters(Input_File_Name_ptr, batch_flag);  
  if(error_flag) return(error_flag);
  
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
  
  /******************* SOLVER LOOP ****************************************
    Loop until command_flag is set to TERMINATE_CODE, most likely by
    Input Parameters read in Post_Processing.
  *************************************************************************/  
  while (Solution_Data.command_flag != TERMINATE_CODE) {
   
    /******************* INITIAL GRID & SOLUTION BLOCKS *********************************
      Create initial mesh and allocate solution variables for specified IBVP/BVP problem. 
    *************************************************************************************/    
    
    // New Calculation (  != CONTINUE_CODE )
    if (Solution_Data.command_flag == EXECUTE_CODE) { 
      CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization. 
      
      error_flag = Initialize_Solution_Blocks(Data,
                                              Solution_Data);
      
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);

      error_flag = Initial_Conditions(Data,
                                      Solution_Data);      
      if (error_flag) return (error_flag);
    } /* endif */
    
    /********************** MAIN SOLVER ****************************************
      Solve IBVP or BVP for conservation form of equations on multi-block 
      solution-adaptive quadrilateral mesh.                                  
    ****************************************************************************/    
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.    
    if (!Data.batch_flag) cout << "\n\n Initialization of mesh and solution data complete on " 
                               << Date_And_Time() << ".";
    
    /****************** OPEN SOLUTION PROGRESS FILES ******************/  
    // Open residual progress file 
    if (CFFC_Primary_MPI_Processor()) {    
       error_flag = Open_Progress_File(Data.residual_file,
	                               Solution_Data.Input.Output_File_Name,
				       Data.number_of_explicit_time_steps);
       if (error_flag && !Data.batch_flag) {
	 cout << "\n ERROR: Unable to open residual progress file for the calculation.\n";
	 cout.flush(); 
       } /* endif */
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);

    // Open other solution progress file(s) 
    if (CFFC_Primary_MPI_Processor()) {    
       error_flag = Open_Other_Solution_Progress_Specialization_Files(Data,
                                                                      Solution_Data);
       if (error_flag && !Data.batch_flag) {
	 cout << "\n ERROR: Unable to open other solution progress file(s) for the calculation.\n";
	 cout.flush(); 
       } /* endif */
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);

    /********************** EXPLICIT **********************************/  
    if ((Data.number_of_explicit_time_steps < Solution_Data.Input.Maximum_Number_of_Time_Steps) ||
	(Solution_Data.Input.Time_Accurate && Solution_Data.Input.Time_Max > Data.Time)) {    
      if (Solution_Data.Input.i_Time_Integration == TIME_STEPPING_MULTIGRID) {
	cerr << "\n MULTIGRID would be here, but not yet. \n"; return error_flag;
      } else {
	error_flag = Hexa_MultiStage_Explicit_Solver<SOLN_pSTATE, SOLN_cSTATE>(Data,
                                                                               Solution_Data);
      } /* endif */
    } /* endif */

    /********************** IMPLICIT **********************************/  
    if (Data.number_of_implicit_time_steps < Solution_Data.Input.NKS_IP.Maximum_Number_of_NKS_Iterations) {
      Hexa_Newton_Krylov_Schwarz_Solver<SOLN_pSTATE, SOLN_cSTATE> NKS(Data, Solution_Data);
      error_flag = NKS.Solve();
    } /* endif */
    
    /****************** CLOSE SOLUTION PROGRESS FILES *****************/  
    // Close residual progress file 
    if (CFFC_Primary_MPI_Processor()) {
       error_flag = Close_Progress_File(Data.residual_file);
       if (error_flag && !Data.batch_flag) {
	 cout << "\n ERROR: Unable to close residual progress file for the calculation.\n";
	 cout.flush(); 
       } /* endif */
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);

    // Close other solution progress file(s) 
    if (CFFC_Primary_MPI_Processor()) {    
       error_flag = Close_Other_Solution_Progress_Specialization_Files(Data,
                                                                       Solution_Data);
       if (error_flag && !Data.batch_flag) {
	 cout << "\n ERROR: Unable to close other solution progress file(s) for the calculation.\n";
	 cout.flush(); 
       } /* endif */
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);


    /***************************************************************
     * Perform solution reconstruction with the final average      *
     * states in order to use the true piecewise representation    *
     * of the solution for post-processing steps, such as solution *
     * plotting or accuracy assessment.                            *
     **************************************************************/

    if (CFFC_Primary_MPI_Processor() && (!batch_flag)) {
      std::cout << "\n\n---------------------------------------\n"
		<< " Reconstruct final solution.\n";
    }

    if ( Solution_Data.Input.i_Reconstruction == RECONSTRUCTION_HIGH_ORDER){
      // Use high-order reconstruction
      HighOrder_Multi_Block::HighOrder_Reconstruction(Solution_Data.Local_Solution_Blocks.Soln_Blks,
						      Data.Local_Adaptive_Block_List,
						      0); 
    } else {
      // Use low-order reconstruction
      HighOrder_Multi_Block::Linear_Reconstruction(Solution_Data.Local_Solution_Blocks.Soln_Blks,
						   Data.Local_Adaptive_Block_List,
						   Solution_Data.Input.i_Limiter);
    } // endif
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);

    if (CFFC_Primary_MPI_Processor() && (!batch_flag)) {
      std::cout << " Solution reconstruction done.\n" 
		<< " ---------------------------------------\n";
    }

    /***************************** POST PROCESSSING *******************************
      Solution calculations complete. Write 3D solution to output and restart files  
      as required, reset solution parameters, and run other cases as specified 
      by input parameters.        
    *******************************************************************************/ 
    
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.    
    if (!Data.batch_flag) cout << "\n";

    error_flag = Hexa_Post_Processing(Data,
                                      Solution_Data);   
    if (error_flag) return (error_flag);

  } //END while

  return (error_flag);
  
}

#endif // _HEXA_SOLVER_INCLUDED
