#ifndef _HEXA_PRE_PROCESSING_INCLUDED
#define _HEXA_PRE_PROCESSING_INCLUDED


#ifndef _HIGHORDER_MULTIBLOCK_INCLUDED
#include "../HighOrderReconstruction/HighOrderMultiBlock.h"
#endif

/********************************************************
 * Routine: Initialize_Solution_Blocks                  *
 *                                                      *
 * Create initial mesh and solution blocks.  The        *
 * computational mesh is either generated directly by   *
 * the program or read in from data files as specified  *
 * by input parameters.                                 *
 *                                                      *
 ********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Initialize_Solution_Blocks(HexaSolver_Data &Data,
		               HexaSolver_Solution_Data<SOLN_pSTATE, SOLN_cSTATE> &Solution_Data) {
  
  //NOTE: None of the initial "GRID" functions pass back an error_flag...
  int error_flag(0);

  /* Create the initial mesh on the primary MPI processor. */

  if (CFFC_Primary_MPI_Processor()) {
    //Outputting solution input parameters
    if (!Data.batch_flag) {
      cout << Solution_Data.Input << "\n";
      cout << "\n Generating the initial mesh.";
      cout.flush(); 
    } /* endif */
    Data.Initial_Mesh.Create_Grid(Solution_Data.Input.Grid_IP);
  } /* endif */

  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
  
  /* Broadcast the mesh to other MPI processors. */

  Data.Initial_Mesh.Broadcast();
  
  /* Initialize solution blocks specializations. */

  error_flag = Initialize_Solution_Blocks_Specializations(Data, Solution_Data);
  if (error_flag) return (error_flag);

  /* Create (allocate) list of hexahedral solution blocks on each processor. */

  if (!Data.batch_flag) {
    cout << "\n Creating multi-block octree data structure, assigning"
      << "\n  solution blocks corresponding to initial mesh, and computing"
      << "\n  geometric coefficients (if required).";
    cout.flush();
  } /* endif */
  
  // FROM AMR 
  error_flag = Create_Initial_Solution_Blocks<SOLN_pSTATE, SOLN_cSTATE>(Data.Initial_Mesh,
		  					                Solution_Data.Local_Solution_Blocks,
							                Solution_Data.Input,
							                Data.Octree,
							                Data.Global_Adaptive_Block_List,
							                Data.Local_Adaptive_Block_List);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return (error_flag);

  /* Create (allocate) the high-order variables in each of the
     local 3D solution blocks, if using high order reconstruction */
  if (Solution_Data.Input.i_Reconstruction == RECONSTRUCTION_HIGH_ORDER){  
    error_flag = HighOrder_Multi_Block::Create_Initial_HighOrder_Variables(Solution_Data.Local_Solution_Blocks.Soln_Blks,
									   Data.Local_Adaptive_Block_List,
									   Solution_Data.Input.Reconstruction_Order);
  }

  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return (error_flag);

  /* Initialization of solution blocks complete, return. */

  return error_flag;

} 

/********************************************************
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
  } /* endif */
  
  //======================================================
  // Read solution from restart data files.
  //======================================================
  if (Solution_Data.Input.i_ICs == IC_RESTART) {
    if (!Data.batch_flag){ 
      cout << "\n Reading solution from restart data files."; 
      cout.flush();
    } /* endif */

    // Read octree
    error_flag = Read_Octree(Data.Octree,
                             Data.Global_Adaptive_Block_List,
                             Data.Local_Adaptive_Block_List,
                             Solution_Data.Local_Solution_Blocks,
                             Solution_Data.Input);
    if (!Data.batch_flag && error_flag) {
      cout << "\n ERROR: Unable to open octree data file on processor "
	   << Data.Local_Adaptive_Block_List.ThisCPU << ".\n";
         cout.flush();
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    
    // Read restart solution files
    error_flag = Solution_Data.Local_Solution_Blocks.Read_Restart_Solution(Solution_Data.Input,  
									   Data.Local_Adaptive_Block_List,
									   Data.number_of_explicit_time_steps,  
									   Data.Time,
									   Data.processor_cpu_time);
    if (!Data.batch_flag && error_flag) {
      cout << "\n ERROR: Unable to open restart input data file(s) "
	   << "on processor "<< CFFC_MPI::This_Processor_Number << ".\n";
      cout.flush();
    } /* endif */ 
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    
    // Calculate y+:
    error_flag = Wall_Distance(Solution_Data.Local_Solution_Blocks.Soln_Blks, // Turbulence function in GENERIC TYPE????
			       Data.Octree, 
			       Data.Local_Adaptive_Block_List);
    if (!Data.batch_flag && error_flag) {
      cout << "\n ERROR: Difficulty determining the wall distance "
	   << "on processor "<< CFFC_MPI::This_Processor_Number
	   << ".\n";
      cout.flush();
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);

    // Calculate wall shear:
    error_flag = Solution_Data.Local_Solution_Blocks.Wall_Shear();
    if (!Data.batch_flag && error_flag) {
      cout << "\n ERROR: Difficulty determining the wall shear "
	   << "on processor "<< CFFC_MPI::This_Processor_Number
	   << ".\n";
      cout.flush();
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);

    // Call specializations:
    error_flag = Hexa_Pre_Processing_Specializations(Data,
                                                     Solution_Data);
    if (!Data.batch_flag && error_flag) {
      cout << "\n ERROR: Issues with pre-processing specializations "
           << "on processor "<< CFFC_MPI::This_Processor_Number
           << ".\n";
      cout.flush();
    } /* endif */

    // Ensure each processor has the correct number of steps and time!!!
    Data.number_of_explicit_time_steps = CFFC_Maximum_MPI(Data.number_of_explicit_time_steps); 
    Data.number_of_implicit_time_steps = CFFC_Maximum_MPI(Data.number_of_implicit_time_steps);
    Data.Time = CFFC_Maximum_MPI(Data.Time);
    Data.processor_cpu_time.cput = CFFC_Maximum_MPI(Data.processor_cpu_time.cput);

    Solution_Data.Input.Maximum_Number_of_Time_Steps = 
      CFFC_Maximum_MPI(Solution_Data.Input.Maximum_Number_of_Time_Steps);
   
  //======================================================
  // Generate initial solution data to begin calculation. 
  //======================================================
  } else {
    // Calculate y+:
    error_flag = Wall_Distance(Solution_Data.Local_Solution_Blocks.Soln_Blks, // Turbulence function in GENERIC TYPE????
			       Data.Octree, 
			       Data.Local_Adaptive_Block_List);
    if (!Data.batch_flag && error_flag) {
      cout << "\n ERROR: Difficulty determining the wall distance "
	   << "on processor "<< CFFC_MPI::This_Processor_Number
	   << ".\n";
      cout.flush();
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);

    // Call ICs:
    error_flag = Solution_Data.Local_Solution_Blocks.ICs(Solution_Data.Input);
    if (!Data.batch_flag && error_flag) {
      cout << "\n ERROR: Issues with initialization of solution data "
	   << "on processor "<< CFFC_MPI::This_Processor_Number
	   << ".\n";
      cout.flush();
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);

    // Call specializations:
    error_flag = Hexa_Pre_Processing_Specializations(Data,
                                                     Solution_Data);
    if (!Data.batch_flag && error_flag) {
      cout << "\n ERROR: Issues with pre-processing specializations "
           << "on processor "<< CFFC_MPI::This_Processor_Number
           << ".\n";
      cout.flush();
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);

  } /* endif */

  /* Send solution information between neighbouring blocks to complete
     prescription of initial data. */

  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization. 

  // First send mesh and geometry information.
  if (Solution_Data.Input.Grid_IP.i_Grid != GRID_PERIODIC_BOX && // temporary check for periodic grids (to be removed, CPTG)
      Solution_Data.Input.Grid_IP.i_Grid != GRID_PERIODIC_BOX_WITH_INFLOW) { 
     error_flag = Send_Messages_Mesh_Geometry_Only<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >
                     (Solution_Data.Local_Solution_Blocks.Soln_Blks,
                      Data.Local_Adaptive_Block_List);
     if (!Data.batch_flag && error_flag) {
         cout    << "\n ERROR: Message passing error during geometry intialization "
              << "on processor "
              << CFFC_MPI::This_Processor_Number
              << ".\n";
         cout.flush();
     } /* endif */
     error_flag = CFFC_OR_MPI(error_flag);
     if (error_flag) return (error_flag);
     // Correct exterior nodes to match with message passed geometry information.
     Solution_Data.Local_Solution_Blocks.Correct_Grid_Exterior_Nodes(Data.Local_Adaptive_Block_List);

     // Fix corner ghost nodes for three-block abutting cases.
     Solution_Data.Local_Solution_Blocks.Fix_Corner_Cells_for_3_Blks_Abutting(Data.Local_Adaptive_Block_List);     
  } /* endif */

  // Now send solution information and data.
  error_flag = Send_Messages<Hexa_Block<SOLN_pSTATE, SOLN_cSTATE> >
                  (Solution_Data.Local_Solution_Blocks.Soln_Blks, 
                   Data.Local_Adaptive_Block_List);
  if (!Data.batch_flag && error_flag) {
      cout << "\n ERROR: Message passing error during solution intialization "
           << "on processor "
           << CFFC_MPI::This_Processor_Number
           << ".\n";
      cout.flush();
   } /* endif */
   error_flag = CFFC_OR_MPI(error_flag);
   if (error_flag) return (error_flag);

  /* Prescribe boundary data consistent with initial data. */

  Solution_Data.Local_Solution_Blocks.BCs(Solution_Data.Input);

  error_flag = Hexa_BCs_Specializations(Data, Solution_Data);
  if (error_flag) return (error_flag);

  /* Output multi-block solution-adaptive quadrilateral mesh statistics. */

  if (!Data.batch_flag) {
     cout << "\n\n Multi-block solution-adaptive hexahedral mesh statistics: "; 
     cout << "\n  -> Number of Root Blocks: "
          << Data.Octree.NR;
     cout << "\n  -> Number of Root Blocks i-direction: "
          << Data.Octree.NRi;
     cout << "\n  -> Number of Root Blocks j-direction: " 
          << Data.Octree.NRj;
     cout << "\n  -> Number of Root Blocks k-direction: " 
          << Data.Octree.NRk;
     cout << "\n  -> Total Number of Used Blocks: " 
          << Data.Octree.countUsedBlocks();
     cout << "\n  -> Total Number of Computational Cells: " 
          << Data.Octree.countUsedCells();
     cout << "\n  -> Number of Mesh Refinement Levels: " 
	  << Data.Octree.highestRefinementLevel()+1;
     cout << "\n  -> Refinement Efficiency: " 
          << Data.Octree.efficiencyRefinement(); 
     cout.flush();
  } /* endif */

  /* End of prepocessing. */  

  return error_flag;

}

/********************************************************
 * Routine: Pre_Processing_Specializations              *
 *                                                      *
 *                                                      *
 ********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Hexa_Pre_Processing_Specializations(HexaSolver_Data &Data,
		                        HexaSolver_Solution_Data<SOLN_pSTATE, SOLN_cSTATE> &Solution_Data) {

  int error_flag(0);
  
  // Call ICs_Specializations:
  error_flag = Solution_Data.Local_Solution_Blocks.ICs_Specializations(Solution_Data.Input);
  if (error_flag) return (error_flag);

  return error_flag;

}

/********************************************************
 * Routine: Initialize_Solution_Blocks_Specializations  *
 *                                                      *
 *                                                      *
 ********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Initialize_Solution_Blocks_Specializations(HexaSolver_Data &Data,
					       HexaSolver_Solution_Data<SOLN_pSTATE, SOLN_cSTATE> &Solution_Data) {
  int error_flag(0);

  return error_flag;

}

/********************************************************
 * Routine: BCs_Specializations                         *
 *                                                      *
 *                                                      *
 ********************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Hexa_BCs_Specializations(HexaSolver_Data &Data,
			     HexaSolver_Solution_Data<SOLN_pSTATE, SOLN_cSTATE> &Solution_Data) {

  int error_flag(0);

  return error_flag;

}

/*****************************************************************
 * Routine: Open_Other_Solution_Progress_Specialization_Files    *
 *                                                               *
 *                                                               *
 *****************************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Open_Other_Solution_Progress_Specialization_Files(HexaSolver_Data &Data,
		                                      HexaSolver_Solution_Data<SOLN_pSTATE, SOLN_cSTATE> &Solution_Data) {

  int error_flag(0);
  
  return error_flag;

}

/*****************************************************************
 * Routine: Close_Other_Solution_Progress_Specialization_Files   *
 *                                                               *
 *                                                               *
 *****************************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Close_Other_Solution_Progress_Specialization_Files(HexaSolver_Data &Data,
                                                       HexaSolver_Solution_Data<SOLN_pSTATE, SOLN_cSTATE> &Solution_Data) {

  int error_flag(0);
  
  return error_flag;

}

/*****************************************************************
 * Routine: Output_Other_Solution_Progress_Specialization_Data   *
 *                                                               *
 *                                                               *
 *****************************************************************/
template<typename SOLN_pSTATE, typename SOLN_cSTATE>
int Output_Other_Solution_Progress_Specialization_Data(HexaSolver_Data &Data,
		                                       HexaSolver_Solution_Data<SOLN_pSTATE, SOLN_cSTATE> &Solution_Data) {

  int error_flag(0);
  
  return error_flag;

}

#endif // _HEXA_PRE_PROCESSING_INCLUDED
