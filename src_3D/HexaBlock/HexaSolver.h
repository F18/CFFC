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

#ifndef _HEXA_POST_PROCESSING_INCLUDED
#include "Hexa_PostProcessing.h"
#endif //_HEXA_POST_PROCESSING_INCLUDED

#ifndef _HEXA_EXPLICIT_SOLVER
#include "Hexa_Explicit_Solver.h"
#endif //_HEXA_EXPLICIT_SOLVER

// #ifndef _NKS_INCLUDED
// #include "../NewtonKrylovSchwarz/NKS.h"
// #endif //_NKS_INCLUDED

/*! *****************************************************
 * Routine: Initialize_Grid                             *
 *                                                      *
 * Create initial mesh.  Read mesh from grid definition *
 * or data files as specified by input parameters.      *
 *                                                      *
 ********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Initialize_Grid(HexaSolver_Data &Data,
		    HexaSolver_Solution_Data<SOLN_pSTATE, SOLN_cSTATE> &Solution_Data) {
  
  //NOTE: None of thes "GRID" functions pass back an error_flag...
  int error_flag(0);

  // The primary MPI processor creates the initial mesh.
  if (CFFC_Primary_MPI_Processor()) {
    Data.Initial_Mesh.Create_Grid(Solution_Data.Input.Grid_IP);
    //Outputting solution input parameters
    if (!Data.batch_flag) {
      cout << Solution_Data.Input << "\n";
      cout.flush(); 
    }
  } 

  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

  //Broadcast the mesh to other MPI processors.
  Data.Initial_Mesh.Broadcast();                    
  
  /* Create (allocate) list of hexahedral solution blocks on each processor. */
  if (!Data.batch_flag) {
    cout << "\n Creating multi-block octree data structure and assigning"
	 << "\n  solution blocks corresponding to initial mesh.";
    cout.flush();
  } 
  
  // FROM AMR 
  Create_Initial_Solution_Blocks<SOLN_pSTATE, SOLN_cSTATE>(Data.Initial_Mesh,
							   Solution_Data.Local_Solution_Blocks,
							   Solution_Data.Input,
							   Data.Octree,
							   Data.Global_Adaptive_Block_List,
							   Data.Local_Adaptive_Block_List);
  return error_flag;
} 


/*! *****************************************************
 * Routine: Initial_Conditions                          *
 *                                                      *
 *                                                      *
 ********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Initial_Conditions(HexaSolver_Data &Data,
		       HexaSolver_Solution_Data<SOLN_pSTATE, SOLN_cSTATE> &Solution_Data) {

  int error_flag(0);
  
  /* Set the initial time level. */
  Data.Time = ZERO;
  Data.number_of_explicit_time_steps = 0;
  Data.number_of_implicit_time_steps = 0;

  /* Set the CPU time to zero. */
  Data.processor_cpu_time.zero();
  Data.total_cpu_time.zero();

  /* Initialize the conserved and primitive state solution variables. */         
  if (!Data.batch_flag) {
    cout << "\n Prescribing initial data.";  cout.flush();
  } 
  
  //Restart 
  if (Solution_Data.Input.i_ICs == IC_RESTART) {
    if (!Data.batch_flag){ 
      cout << "\n Reading solution from restart data files."; 
      cout.flush();
    }

    // Read Restart Octree
    // error_flag = Read_Octree(Octree,Input);
    if (!Data.batch_flag && error_flag) {
      cout << "\n ERROR: Unable to open Octree data file on processor "
	   << Data.Local_Adaptive_Block_List.ThisCPU<< ".\n";
         cout.flush();
    } 
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    
    // Read Restart Solution Files
    error_flag = Solution_Data.Local_Solution_Blocks.Read_Restart_Solution(Solution_Data.Input,  
									   Data.Local_Adaptive_Block_List,
									   Data.number_of_explicit_time_steps,  
									   Data.Time,
									   Data.processor_cpu_time);
    if (!Data.batch_flag && error_flag) {
      cout << "\n  ERROR: Unable to open restart input data file(s) "
	   << "on processor "<< CFFC_MPI::This_Processor_Number<< ".\n";
      cout.flush();
    } 
    
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    
    // Ensure each processor has the correct number of steps and time!!!
    Data.number_of_explicit_time_steps = CFFC_Maximum_MPI(Data.number_of_explicit_time_steps); 
    Data.number_of_implicit_time_steps = CFFC_Maximum_MPI(Data.number_of_implicit_time_steps);
    Data.Time = CFFC_Maximum_MPI(Data.Time);
    Data.processor_cpu_time.cput = CFFC_Maximum_MPI(Data.processor_cpu_time.cput);

    Solution_Data.Input.Maximum_Number_of_Time_Steps = 
      CFFC_Maximum_MPI(Solution_Data.Input.Maximum_Number_of_Time_Steps);
   
    // NON RESTART 
  } else {
    error_flag = Wall_Distance(Solution_Data.Local_Solution_Blocks.Soln_Blks,  // Turbulence function in GENERIC TYPE????
			       Data.Octree, 
			       Data.Local_Adaptive_Block_List);
    if (!Data.batch_flag && error_flag) {
      cout << "\n  ERROR: Difficulty determining the wall distance "
	   << "on processor "<< CFFC_MPI::This_Processor_Number
	   << ".\n";
      cout.flush();
    } 
    
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    
    Solution_Data.Local_Solution_Blocks.ICs(Solution_Data.Input);
  } 

   /* Send solution information between neighbouring blocks to complete
      prescription of initial data. */   
  Send_All_Messages<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >(Solution_Data.Local_Solution_Blocks.Soln_Blks,
							   Data.Local_Adaptive_Block_List,
							   Solution_Data.Local_Solution_Blocks.Soln_Blks[0].NumVar(),
							   OFF);
  
  /* Prescribe boundary data consistent with initial data. */
  Solution_Data.Local_Solution_Blocks.BCs(Solution_Data.Input);
  
  return error_flag;

}


/*! *****************************************************
 * Routine: HexaSolver                                  *
 *                                                      *
 * Computes solutions to the governing PDEs on          *
 * multi-block AMR mesh composed of body-fitted         *
 * hexahedral solution blocks.                          *
 *                                                      *
 ********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int HexaSolver(char *Input_File_Name_ptr,int batch_flag){
   
  int error_flag(0);

  //Non Solution Specific
  HexaSolver_Data                                     Data(batch_flag);  
  //Solution Specific 
  HexaSolver_Solution_Data<SOLN_pSTATE, SOLN_cSTATE>  Solution_Data;   

  /*! *************** INPUT PARAMETERS  **********************************
    Set default values for the input solution parameters and then read user 
    specified input values from the specified input parameter file.               
  *************************************************************************/    
  error_flag = Solution_Data.Get_Input_Parameters(Input_File_Name_ptr,batch_flag);  
  error_flag = CFFC_OR_MPI(error_flag);
  if(error_flag) return(error_flag);
  
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
  
  /*! **************** SOLVER LOOP ****************************************
    Loop until command_flag is set to TERMINATE_CODE, most likely by
    Input Parameters read in Post_Processing.
  *************************************************************************/  
  while(Solution_Data.command_flag != TERMINATE_CODE){
   
    /*! **************** INITIAL GRID & SOLUTION SPACE **********************
      Create initial mesh and allocate Chem2D solution variables for 
      specified IBVP/BVP problem. 
    *************************************************************************/    
    
    // New Calculation (  != CONTINUE_CODE )
    if(Solution_Data.command_flag == EXECUTE_CODE) { 

      CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization. 

      error_flag = Initialize_Grid(Data,Solution_Data);      
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);

      error_flag = Initial_Conditions(Data,Solution_Data);      
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
    }
    
    /*! *********************** MAIN SOLVER ****************************************
      Solve IBVP or BVP for conservation form of equations on multi-block 
      solution-adaptive quadrilateral mesh.                                  
    ****************************************************************************/  
    
    /* 
      Need logic for explicit, implicit, multigrid, & unsteady, with combinations
      of both...    

      Input.Time_Accurate                  (true/false)
      Input.NKS_IP.Dual_Time_Stepping      (true/false) 

      Input.Time_Max                       (0 -> n) seconds
      Input.i_Time_Integration             (Explicit_Euler,Explicit_Runge_Kutta,Multistage_Optimal_Smoothing, etc.)

      Input.NKS_IP.Maximum_Number_of_NKS_Iterations   (0 -> n) its
      Input.Maximum_Number_of_Time_Steps              (0 -> n) its  explicit
   */
  
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.    
    
    /********************** EXPLICIT **********************************/  
    if( (Data.number_of_explicit_time_steps < Solution_Data.Input.Maximum_Number_of_Time_Steps) ||
	(Solution_Data.Input.Time_Accurate && Solution_Data.Input.Time_Max > Data.Time) ) {    
      
      if(Solution_Data.Input.i_Time_Integration == TIME_STEPPING_MULTIGRID){
	cerr<< "\n MULTIGRID would be here, but not yet. \n"; return error_flag;
      } else {
	error_flag = Hexa_MultiStage_Explicit_Solver<SOLN_pSTATE, SOLN_cSTATE>
	  (Data,Solution_Data);
      }
    }

    /********************** IMPLICIT **********************************/  
    if(Data.number_of_implicit_time_steps < Solution_Data.Input.NKS_IP.Maximum_Number_of_NKS_Iterations){
      cerr<< "\n NKS would be here, but not yet. \n"; return error_flag;
      //   	error_flag = Hexa_Newton_Krylov_Schwarz_Solver<SOLN_pSTATE, SOLN_cSTATE> 
      //       	  (Data,Solution_Data);
    }
     
    /*! ************************** POST PROCESSSING *******************************
      Solution calculations complete. Write 3D solution to output and restart files  
      as required, reset solution parameters, and run other cases as specified 
      by input parameters.        
    *******************************************************************************/     
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.    

    error_flag = Hexa_PostProcessing(Data,Solution_Data);   
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);

  } //END while

  return (error_flag);
  
}


#endif // _HEXA_SOLVER_INCLUDED


