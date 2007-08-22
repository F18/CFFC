/**********************************************************************
 * HighTemp2DQuadSolvers.cc: 2D High-Temp equation multi-block        *
 *                           quadrilateral mesh solvers.              *
 **********************************************************************/

#include "HighTemp2DQuad.h"
#include "../FASMultigrid2D/FASMultigrid2D.h"
#include "HighTemp2DQuadMultigrid.h"
#include "../NewtonKrylovSchwarz2D/NKS2D.h"
#include "HighTemp2DQuadNKS.h"
#include "../CFD/EllipticOperatorAnalysis2D.h"

/**********************************************************************
 * Routine: HighTemp2DQuadSolver                                      *
 *                                                                    *
 * Computes solutions to 2D Navier-Stokes HT  equations on 2D         *
 * quadrilateral multi-block solution-adaptive mesh.                  *
 *                                                                    *
 **********************************************************************/
int HighTemp2DQuadSolver(char *Input_File_Name_ptr, int batch_flag) {

  /********************************************************************
   * Local variable declarations.                                     *
   ********************************************************************/

  // HighTemp2D input variables and parameters.
  HighTemp2D_Input_Parameters Input_Parameters;

  // Multi-block solution-adaptive quadrilateral mesh solution variables.
  Grid2D_Quad_Block             **MeshBlk;
  QuadTreeBlock_DataStructure     QuadTree;
  AdaptiveBlockResourceList       List_of_Global_Solution_Blocks;
  AdaptiveBlock2D_List            List_of_Local_Solution_Blocks;
  HighTemp2D_Quad_Block      *Local_SolnBlk;

  // Multigrid declaration.
  FAS_Multigrid2D_Solver<HighTemp2D_cState, 
                         HighTemp2D_Quad_Block, 
                         HighTemp2D_Input_Parameters> MGSolver;

  // Define residual file and cpu time variables.
  ofstream residual_file;
 
  // processor_cpu_time is the real deal whereas total_cpu_time is used
  // only to print. Whenever total_cpu_time is to be printed, it is
  // calculated from processor_cpu_time. processor_cpu_time is passed from
  // solver to solver (FAS, plain time stepping, NKS, etc) and is the only
  // thing passed. Those solvers who wish to track their own CPU time
  // should record the value of processor_cpu_time before and after their
  // execution but should only ever update processor_cpu_time. This is 
  // especially important for cases with mesh refinements.
  CPUTime processor_cpu_time, total_cpu_time;

  // NKS time needs to be allocated here since 
  // it needs to persist after a mesh refinement. 
  CPUTime NKS_total_cpu_time;

  // These are for the current mesh only. Used only to verify CPU times.
  time_t start_explicit = 0, end_explicit = 0;
  time_t start_NKS = 0, end_NKS = 0;

  // Other local solution variables.
  int number_of_time_steps, start_number_of_time_steps, first_step,
      command_flag, error_flag, line_number,
      perform_explicit_time_marching, limiter_freezing_flag;
  double Time, dTime;
  double residual_l1_norm, residual_l2_norm, residual_max_norm;
  double first_residual_l2_norm = -1.0, ratio_residual_l2_norm = 1.0;
  bool mgsolver_is_allocated = false;

  /********************************************************************  
   * Set default values for the input solution parameters and then    *
   * read user specified input values from the specified input        *
   * parameter file.                                                  *
   ********************************************************************/

  // The primary MPI processor processes the input parameter file.
  if (CFFC_Primary_MPI_Processor()) {
    if (!batch_flag)
      cout << "\n Reading HighTemp2D input data file `"
	   << Input_File_Name_ptr << "'." << endl;
    error_flag = Process_Input_Control_Parameter_File(Input_Parameters,
						      Input_File_Name_ptr,
						      command_flag);
    if (!batch_flag && !error_flag && command_flag != EXECUTE_ZERO_STEPS_CODE) {
      cout << Input_Parameters << endl;
      cout.flush();
    } else if (error_flag) {
      cout << endl << "HighTemp2D ERROR: During processing of input parameters." << endl << endl;
    }
  } else {
    error_flag = 0;
  }
  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Broadcast input solution parameters to other MPI processors.
  CFFC_Broadcast_MPI(&error_flag,1);
  if (error_flag) return error_flag;
  CFFC_Broadcast_MPI(&command_flag,1);
  if (command_flag == TERMINATE_CODE) return 0;
  Broadcast_Input_Parameters(Input_Parameters);

  /********************************************************************
   * Create initial mesh and allocate HighTemp2D solution             *
   * variables for the specified IBVP/BVP problem.                    *
   ********************************************************************/

 execute_new_calculation: ;

  // Synchronize processors.
  CFFC_Barrier_MPI();

  // Create initial mesh.  Read mesh from grid definition or data  
  // files when specified by input parameters.

  // The primary MPI processor creates the initial mesh.
  if (CFFC_Primary_MPI_Processor()) {

    if (!batch_flag) 
      cout << "\n Creating (or reading) initial quadrilateral multi-block mesh.";

    MeshBlk = NULL;
    MeshBlk = Multi_Block_Grid(MeshBlk,Input_Parameters);

    if (MeshBlk == NULL) {
      error_flag = 1;
    } else if (Check_Multi_Block_Grid(MeshBlk,
				      Input_Parameters.Number_of_Blocks_Idir,
				      Input_Parameters.Number_of_Blocks_Jdir)) {
      error_flag = 1;
    } else {
      error_flag = 0;
    }

    if (error_flag) {
      cout << "\n HighTemp2D ERROR: Unable to create valid HighTemp2D multi-block mesh.\n";
      cout.flush();
    } else {
//        if (Input_Parameters.i_Grid == GRID_NASA_ROTOR_37) {
//  	if (!batch_flag) cout << "\n Writing geometry and flow field data files for NASA Rotor 37.";
//  	Input_Parameters.NASA_Rotor37.outputTP_UpstreamFlowConditions("NASARotor37_UpstreamProfile.dat");
//  	Input_Parameters.NASA_Rotor37.outputTP_DownstreamFlowConditions("NASARotor37_DownstreamProfile.dat");
//  	Input_Parameters.NASA_Rotor37.outputTP_Geometry("NASARotor37_Geometry.dat");
//        } else if (Input_Parameters.i_Grid == GRID_NASA_ROTOR_67) {
//  	if (!batch_flag) cout << "\n Writing geometry and flow field data files for NASA Rotor 67.";
//  	Input_Parameters.NASA_Rotor67.outputTP_UpstreamFlowConditions("NASARotor67_UpstreamProfile.dat");
//  	Input_Parameters.NASA_Rotor67.outputTP_DownstreamFlowConditions("NASARotor67_DownstreamProfile.dat");
//  	Input_Parameters.NASA_Rotor67.outputTP_Geometry("NASARotor67_Geometry.dat");
//  	Input_Parameters.NASA_Rotor67.outputTP_FlowField2D(Input_Parameters.Rotor_Percent_Span,
//  							   "NASARotor67_InterbladeFlow.dat");
//        }
    }
  } else {
    MeshBlk = NULL;
  }

  // Synchronize processors.
  CFFC_Barrier_MPI();

  // Broadcast the mesh to other MPI processors.
  CFFC_Broadcast_MPI(&error_flag,1);
  if (error_flag) return error_flag;
  MeshBlk = Broadcast_Multi_Block_Grid(MeshBlk,Input_Parameters);

  // Create (allocate) multi-block quadtree data structure, create
  // (allocate) array of local 2D HighTemp equation solution blocks, 
  // assign and create (allocate) 2D HighTemp equation solution blocks
  // corresponding to the initial mesh.
  if (!batch_flag) cout << "\n Creating multi-block quadtree data structure and assigning"
                        << "\n HighTemp2D solution blocks corresponding to initial mesh.";
  Local_SolnBlk = CreateInitialSolutionBlocks(MeshBlk,
					      Local_SolnBlk,
					      Input_Parameters,
					      QuadTree,
					      List_of_Global_Solution_Blocks,
					      List_of_Local_Solution_Blocks);
  if (Local_SolnBlk == NULL) return 1;

  /********************************************************************
   * Initialize HighTemp2D solution variables.                        *
   ********************************************************************/

  // Set the initial time level.
  Time = ZERO;
  number_of_time_steps = 0;

  // Set the CPU time to zero.
  processor_cpu_time.zero();
  total_cpu_time.zero();
	NKS_total_cpu_time.zero();

  // Initialize the conserved and primitive state solution variables.
  if (!batch_flag) cout << "\n Prescribing HighTemp2D initial data.";
  if (Input_Parameters.i_ICs == IC_RESTART) {
    if (!batch_flag) cout << "\n Reading HighTemp2D solution from restart data files.";
    // Read the quadtree restart file.
    error_flag = Read_QuadTree(QuadTree,
			       List_of_Global_Solution_Blocks,
			       List_of_Local_Solution_Blocks,
			       Input_Parameters);
    if (error_flag) {
      cout << "\n HighTemp2D ERROR: Unable to open HighTemp2D quadtree data file "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;
    // Allocate the message buffers.
    Allocate_Message_Buffers(List_of_Local_Solution_Blocks,
			     NUM_VAR_HIGHTEMP2D+NUM_COMP_VECTOR2D);
    // Read the solution block restart files.
    error_flag = Read_Restart_Solution(Local_SolnBlk,
				       List_of_Local_Solution_Blocks,
				       Input_Parameters,
				       number_of_time_steps,
				       Time,
				       processor_cpu_time);
    if (error_flag) {
      cout << endl << " HighTemp2D ERROR: Unable to open HighTemp2D restart "
	   << "input data file(s) on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Determine the distance to the nearest wall distance.
    error_flag = Determine_Wall_Distance(Local_SolnBlk,
					 QuadTree,
					 List_of_Local_Solution_Blocks,
					 Input_Parameters);
    if (error_flag) {
      cout << " HighTemp2D ERROR: During wall distance calculation."
	   << "  Error #" << error_flag << "." << endl;
    }
    CFFC_Broadcast_MPI(&error_flag,1);
    if (error_flag) return error_flag;

    // Ensure each processor has the correct time and time.
    number_of_time_steps = CFFC_Maximum_MPI(number_of_time_steps);
    Time = CFFC_Maximum_MPI(Time);
    processor_cpu_time.cput = CFFC_Maximum_MPI(processor_cpu_time.cput);
    Input_Parameters.Maximum_Number_of_Time_Steps = CFFC_Maximum_MPI(Input_Parameters.Maximum_Number_of_Time_Steps);

    // Synchronize processors.
    CFFC_Barrier_MPI();
    // Broadcast input solution parameters to other MPI processors.
    CFFC_Broadcast_MPI(&error_flag,1);
    if (error_flag != 0) return error_flag;
    CFFC_Broadcast_MPI(&command_flag,1);
    if (command_flag == TERMINATE_CODE) return 0;
    Broadcast_Input_Parameters(Input_Parameters);

  } else {

    // Apply initial conditions.
    ICs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

  }

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Send solution information between neighbouring blocks to complete
  // prescription of initial data.
  error_flag = Send_All_Messages(Local_SolnBlk,
				 List_of_Local_Solution_Blocks,
				 NUM_COMP_VECTOR2D,
				 ON);
  if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk,
						  List_of_Local_Solution_Blocks,
						  NUM_VAR_HIGHTEMP2D,
						  OFF);
  if (error_flag) {
    cout << "\n HighTemp2D ERROR: Message passing error during HighTemp2D "
	 << "solution intialization on processor "
	 << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Prescribe boundary data consistent with initial data.
  BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

  // Perform uniform, boundary, and initial mesh refinement.
  if (Input_Parameters.i_ICs != IC_RESTART) {

    // Perform uniform mesh refinement.
    if (!batch_flag) cout << "\n Performing HighTemp2D uniform mesh refinement.";
    error_flag = Uniform_Adaptive_Mesh_Refinement(Local_SolnBlk,
						  Input_Parameters,
						  QuadTree,
						  List_of_Global_Solution_Blocks,
						  List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n HighTemp2D ERROR: Uniform AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Perform boundary mesh refinement.
    if (!batch_flag) cout << "\n Performing HighTemp2D boundary mesh refinement.";
    error_flag = Boundary_Adaptive_Mesh_Refinement(Local_SolnBlk,
						   Input_Parameters,
						   QuadTree,
						   List_of_Global_Solution_Blocks,
						   List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n HighTemp2D ERROR: Boundary AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Perform flat-plate mesh refinement.
    if (!batch_flag) cout << "\n Performing HighTemp2D flat-plate mesh refinement.";
    error_flag = Flat_Plate_Adaptive_Mesh_Refinement(Local_SolnBlk,
						     Input_Parameters,
						     QuadTree,
						     List_of_Global_Solution_Blocks,
						     List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n HighTemp2D ERROR: Flat-plate AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Perform initial mesh refinement.
    if (!batch_flag) cout << "\n Performing HighTemp2D initial mesh refinement.";
    error_flag = Initial_Adaptive_Mesh_Refinement(Local_SolnBlk,
						  Input_Parameters,
						  QuadTree,
						  List_of_Global_Solution_Blocks,
						  List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n HighTemp2D ERROR: Initial AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // MPI barrier to ensure processor synchronization.
    CFFC_Barrier_MPI();

    // Send solution information between neighbouring blocks to complete
    // prescription of initial data.
    error_flag = Send_All_Messages(Local_SolnBlk,
 				   List_of_Local_Solution_Blocks,
				   NUM_COMP_VECTOR2D,
				   ON);
    if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk,
						    List_of_Local_Solution_Blocks,
						    NUM_VAR_HIGHTEMP2D,
						    OFF);
    if (error_flag) {
      cout << "\n HighTemp2D ERROR: Message passing error during HighTemp2D solution intialization "
	   << "on processor " << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

  }
  
   // Compute wall injection speed.
  error_flag = Wall_Injection_Speed(Local_SolnBlk,
				    List_of_Global_Solution_Blocks,
				    List_of_Local_Solution_Blocks,
				    Input_Parameters);
  if (error_flag) {
    cout << "\n HighTemp2D ERROR: HT2D error during vw-plus calculation "
	 << "on processor " << List_of_Local_Solution_Blocks.ThisCPU << endl;
  }

  // Compute dimensionless wall distance.
  error_flag = Dimensionless_Wall_Distance(Local_SolnBlk,
					   //QuadTree,
					   List_of_Local_Solution_Blocks,
					   Input_Parameters);
  if (error_flag) {
    cout << "\n HighTemp2D ERROR: HT2D error during y-plus calculation "
	 << "on processor " << List_of_Local_Solution_Blocks.ThisCPU << endl;
  }

  // Compute dimensionless wall injection speed.
  error_flag = Dimensionless_Wall_Injection_Speed(Local_SolnBlk,
						  List_of_Local_Solution_Blocks,
						  Input_Parameters);
  if (error_flag) {
    cout << "\n HighTemp2D ERROR: HT2D error during vw-plus calculation "
	 << "on processor " << List_of_Local_Solution_Blocks.ThisCPU << endl;
  }

  // Output multi-block solution-adaptive quadrilateral mesh statistics.
  if (!batch_flag) {
    cout << "\n\n Multi-block solution-adaptive quadrilateral mesh statistics: ";
    cout << "\n  -> Number of Root Blocks i-direction: "
	 << QuadTree.NRi;
    cout << "\n  -> Number of Root Blocks j-direction: "
	 << QuadTree.NRj;
    cout << "\n  -> Total Number of Used Blocks: "
	 << QuadTree.countUsedBlocks();
    cout << "\n  -> Total Number of Computational Cells: "
	 << QuadTree.countUsedCells();
    cout << "\n  -> Number of Mesh Refinement Levels: "
	 << QuadTree.highestRefinementLevel()+1;
    cout << "\n  -> Refinement Efficiency: "
	 << QuadTree.efficiencyRefinement() << endl;
  }

//       if (CFFC_Primary_MPI_Processor()) {
//          for (int j_blk = 0; j_blk < QuadTree.Nblk; j_blk++) {
//             for (int i_blk = 0; i_blk < QuadTree.Ncpu; i_blk++) {
//                if (QuadTree.Blocks[i_blk][j_blk] != NULL) {
//                   cout << "\n cpu = " 
//                        << i_blk
//                        << " blk = "
//                        << j_blk
//                        << " blk = "
//                        << QuadTree.Blocks[i_blk][j_blk]->block;
//                } else {
//                   cout << "\n cpu = " 
//                        << i_blk
//                        << " blk = "
//                        << j_blk;
//                }
//             }
//          }
//       }

	if (Input_Parameters.Morton) {
		if (!batch_flag) { cout << "\n\n Applying Morton re-ordering algorithm to initial solution blocks. "; }
		error_flag = Morton_ReOrdering_of_Solution_Blocks(QuadTree, 
				List_of_Global_Solution_Blocks, 
				List_of_Local_Solution_Blocks, 
				Local_SolnBlk, 
				Input_Parameters, 
				number_of_time_steps, 
				Time, 
				processor_cpu_time); 
		if (error_flag) {
			cout <<"\n HighTemp2D ERROR: Morton re-ordering error on processor "
				<< List_of_Local_Solution_Blocks.ThisCPU << ".\n";
		} 
		error_flag = CFFC_OR_MPI(error_flag);
		if (error_flag) { return error_flag; }

		if (!batch_flag) { cout << "\n Outputting space filling curve showing block loading for CPUs."; }
		Morton_SFC_Output_Tecplot(Local_SolnBlk, 
				Input_Parameters, 
				List_of_Local_Solution_Blocks);
	} // End of Morton Re-ordering

	if (command_flag == EXECUTE_ZERO_STEPS_CODE) {
		// Maybe I could do Initial_AMR instead of this. 
		// Maybe I could not.
		//   -- Alistair Wood Wed Jun 27 2007 
		start_number_of_time_steps = number_of_time_steps;
		goto postprocess_current_calculation;
	}

  /********************************************************************
   * Solve IBVP or BVP for conservation form of 2D Navier-Stokes      *
   * equations on multi-block solution-adaptive quadrilateral mesh.   *
   ********************************************************************/

 continue_existing_calculation: ;

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();
  time(&start_explicit);
  start_number_of_time_steps = number_of_time_steps;

  if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {
    // Allocate memory for multigrid solver.
    error_flag = MGSolver.allocate(Local_SolnBlk,
				   &QuadTree,
				   &List_of_Global_Solution_Blocks,
				   &List_of_Local_Solution_Blocks,
				   &Input_Parameters);
    mgsolver_is_allocated = true;
    if (error_flag) {
      cout << "\n HighTemp2D ERROR: Unable to allocate memory for multigrid solver.\n";
      cout.flush();
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Execute multigrid solver.
    error_flag = MGSolver.Execute(batch_flag,
				  number_of_time_steps,
				  Time,
				  processor_cpu_time,
				  total_cpu_time,
				  residual_file);
    if (error_flag) {
      cout << "\n HighTemp2D ERROR: Multigrid error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
      cout.flush();
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

  } else if (Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {

    // Allocate memory for the multigrid solver.
    error_flag = MGSolver.allocate(Local_SolnBlk,
				   &QuadTree,
				   &List_of_Global_Solution_Blocks,
				   &List_of_Local_Solution_Blocks,
				   &Input_Parameters);
		mgsolver_is_allocated = true;
    if (error_flag) {
      cout << "\n HighTemp2D ERROR: Unable to allocate memory for DTS multigrid"
	   << "\n solver.  Error number = " << error_flag << endl;
    }
    CFFC_Broadcast_MPI(&error_flag,1);
    if (error_flag) return error_flag;

    // Execute DTS FAS multigrid solver.
    error_flag = MGSolver.DTS_Multigrid_Solution(batch_flag,
						 number_of_time_steps,
						 Time,
						 processor_cpu_time,
						 total_cpu_time,
						 residual_file);
    if (error_flag) {
      cout << "\n HighTemp2D ERROR: Error during DTS multigrid solution.  Error"
	   << "\n number = " << error_flag << endl;
    }
    CFFC_Broadcast_MPI(&error_flag,1);
    if (error_flag) return error_flag;
    if (error_flag) {
      cout << "\n HighTemp2D ERROR: Error during DTS multigrid solution.\n";
      cout.flush();
    }
    CFFC_Broadcast_MPI(&error_flag,1);
    if (error_flag) return error_flag;

  } else {

    // Open residual file.
    first_step = 1;
    limiter_freezing_flag = OFF;

    if (CFFC_Primary_MPI_Processor()) {
      error_flag = Open_Progress_File(residual_file,
				      Input_Parameters.Output_File_Name,
				      number_of_time_steps);
      if (error_flag) {
	cout << "\n HighTemp2D ERROR: Unable to open residual file for HighTemp2D calculation.\n";
	cout.flush();
      }
    }
    // MPI barrier to ensure processor synchronization.  
    CFFC_Barrier_MPI();
    CFFC_Broadcast_MPI(&error_flag,1);
    if (error_flag) return error_flag;
    // Reset the CPU time.
    processor_cpu_time.reset();

    // Perform required number of iterations (time steps).
    if ((!Input_Parameters.Time_Accurate &&
	 Input_Parameters.Maximum_Number_of_Time_Steps > 0 &&
	 number_of_time_steps < Input_Parameters.Maximum_Number_of_Time_Steps) ||
	(Input_Parameters.Time_Accurate &&
	 Input_Parameters.Time_Max > Time)) {

      if (!batch_flag) cout << "\n Beginning HighTemp2D computations on "
			    << Date_And_Time() << ".\n\n";

      if ((!Input_Parameters.Time_Accurate &&
	   Input_Parameters.Maximum_Number_of_Time_Steps > 0 &&
	   number_of_time_steps < Input_Parameters.Maximum_Number_of_Time_Steps) ||
  	  (Input_Parameters.Time_Accurate &&
	   Input_Parameters.Time_Max > Time)) {
	perform_explicit_time_marching = ON;
      } else {
	perform_explicit_time_marching = OFF;
      }

      while (perform_explicit_time_marching) {

	// Periodically refine the mesh (AMR).
        if (Input_Parameters.AMR) {
	  if (!first_step &&
	      number_of_time_steps-Input_Parameters.AMR_Frequency*
	      (number_of_time_steps/Input_Parameters.AMR_Frequency) == 0 ) {
	    if (!batch_flag) cout << "\n\n Refining Grid.  Performing adaptive mesh refinement at n = "
				  << number_of_time_steps << ".";
	    Evaluate_Limiters(Local_SolnBlk,List_of_Local_Solution_Blocks);
	    error_flag = Adaptive_Mesh_Refinement(Local_SolnBlk,
						  Input_Parameters,
						  QuadTree,
						  List_of_Global_Solution_Blocks,
						  List_of_Local_Solution_Blocks);
	    if (error_flag) {
	      cout << "\n HighTemp2D ERROR: HighTemp2D AMR error number "
		   << error_flag << " on processor "
		   << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
	      cout.flush();
	    }
	    error_flag = CFFC_OR_MPI(error_flag);
	    if (error_flag) {
	      command_flag = Output_Tecplot(Local_SolnBlk,
					    List_of_Local_Solution_Blocks,
					    Input_Parameters,
					    number_of_time_steps,
					    Time);
	      return error_flag;
	    }
	    if (!batch_flag) {
	      cout << "\n New multi-block solution-adaptive quadrilateral mesh statistics: ";
	      cout << "\n  -> Number of Root Blocks i-direction: "
		   << QuadTree.NRi;
	      cout << "\n  -> Number of Root Blocks j-direction: "
		   << QuadTree.NRj;
	      cout << "\n  -> Total Number of Used Blocks: "
		   << QuadTree.countUsedBlocks();
	      cout << "\n  -> Total Number of Computational Cells: "
		   << QuadTree.countUsedCells();
	      cout << "\n  -> Number of Mesh Refinement Levels: "
		   << QuadTree.highestRefinementLevel()+1;
	      cout << "\n  -> Refinement Efficiency: "
		   << QuadTree.efficiencyRefinement() << endl;
	      cout.flush();
	    }
	  }
        }

	// Periodically re-order the solution blocks on the processors
	// by applying the Morton ordering space filling curve algorithm.
	//
	// Probably this should be coupled with AMR Frequency.
	if (Input_Parameters.Morton &&
  	    !first_step &&
	    number_of_time_steps % Input_Parameters.Morton_Reordering_Frequency == 0) {
	   if (!batch_flag) {
	     cout << "\n\n Applying Morton re-ordering algorithm to solution blocks at n = "
		  << number_of_time_steps << ".";
	   }
	   error_flag = Morton_ReOrdering_of_Solution_Blocks(QuadTree, 
		 			                     List_of_Global_Solution_Blocks, 
					                     List_of_Local_Solution_Blocks, 
					                     Local_SolnBlk, 
					                     Input_Parameters, 
					                     number_of_time_steps, 
					                     Time, 
					                     processor_cpu_time); 
	   if (error_flag) {
	      cout <<"\n HighTemp2D ERROR: Morton re-ordering error on processor "
		   << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
	   } 
	   error_flag = CFFC_OR_MPI(error_flag);
	   if (error_flag) { return error_flag; }

	   if (!batch_flag) { cout << "\n Outputting space filling curve showing block loading for CPUs."; }
	     Morton_SFC_Output_Tecplot(Local_SolnBlk, 
		   	               Input_Parameters, 
				       List_of_Local_Solution_Blocks);
	} // End of Morton Re-ordering

	// Determine local and global time steps.
	dTime = CFL(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);
	// Find global minimum time step for all processors.
	dTime = CFFC_Minimum_MPI(dTime);
	if (Input_Parameters.Time_Accurate) {
	  if ((Input_Parameters.i_Time_Integration != 
	       TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) &&
	      (Time + Input_Parameters.CFL_Number*dTime > 
	       Input_Parameters.Time_Max)) {
	    dTime = (Input_Parameters.Time_Max-Time)/Input_Parameters.CFL_Number;
	  } else if (Time + Input_Parameters.CFL_Number*dTime*
		     MultiStage_Optimally_Smoothing(Input_Parameters.N_Stage,
						    Input_Parameters.N_Stage,
						    Input_Parameters.i_Limiter) > 
		     Input_Parameters.Time_Max) {
	    dTime = (Input_Parameters.Time_Max-Time)/
	            (Input_Parameters.CFL_Number*
		     MultiStage_Optimally_Smoothing(Input_Parameters.N_Stage,
						    Input_Parameters.N_Stage,
						    Input_Parameters.i_Limiter));
	  }
	}
	// Set global time step.
	if (!Input_Parameters.Local_Time_Stepping) {
	  Set_Global_TimeStep(Local_SolnBlk,
			      List_of_Local_Solution_Blocks,
			      dTime);
	}

	// Determine the L1, L2, and max norms of the solution residual.
	residual_l1_norm = L1_Norm_Residual(Local_SolnBlk,List_of_Local_Solution_Blocks);
	residual_l1_norm = CFFC_Summation_MPI(residual_l1_norm);
      
	residual_l2_norm = L2_Norm_Residual(Local_SolnBlk,List_of_Local_Solution_Blocks);
	residual_l2_norm = sqr(residual_l2_norm);
	residual_l2_norm = CFFC_Summation_MPI(residual_l2_norm);
	residual_l2_norm = sqrt(residual_l2_norm);

	residual_max_norm = Max_Norm_Residual(Local_SolnBlk,List_of_Local_Solution_Blocks);
	residual_max_norm = CFFC_Maximum_MPI(residual_max_norm);

	//  Why ignore the first two steps and not just the first step?
	//  Because when the time integration type is set to Explicit Euler
	//  then the L2 norm on the first step is many orders of magnitude
	//  greater than during the rest of the simulation. Why? I don't know
	//  and I don't care. So why not ignore the first two steps only for
	//  Explicit Euler? Keep it simple, kiddo, keep it simple. Besides,
	//  we only use the ratio norm to print and to stop - neither of
	//  which require absolute accuracy.
	//    -- Alistair Wood. Tue Jul 17 2007.
	if (number_of_time_steps < 2) {
		ratio_residual_l2_norm = 1.0;
	} else {
		if (first_residual_l2_norm <= 0.0) { first_residual_l2_norm = residual_l2_norm; }
		ratio_residual_l2_norm = residual_l2_norm / first_residual_l2_norm;
	}

	// Update CPU time used for the calculation so far.
	processor_cpu_time.update();
	total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput);

	// Periodically save restart solution files.
	if (!first_step &&
	    number_of_time_steps-Input_Parameters.Restart_Solution_Save_Frequency*
	    (number_of_time_steps/Input_Parameters.Restart_Solution_Save_Frequency) == 0) {
	  if (!batch_flag) 
	    cout << "\n\n  Saving HighTemp2D solution to restart data file(s) after"
		 << " n = " << number_of_time_steps << " steps (iterations).";
	  // Write the quadtree restart file.
	  error_flag = Write_QuadTree(QuadTree,Input_Parameters);
	  if (error_flag) {
	    cout << "\n HighTemp2D ERROR: Unable to open HighTemp2D quadtree data file "
		 << "on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  }
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return error_flag;
	  // Write the solution block restart files.
	  error_flag = Write_Restart_Solution(Local_SolnBlk,
					      List_of_Local_Solution_Blocks,
					      Input_Parameters,
					      number_of_time_steps,
					      Time,
					      processor_cpu_time);
	  if (error_flag) {
	    cout << "\n HighTemp2D ERROR: Unable to open HighTemp2D restart output data file(s) "
		 << "on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  }
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return error_flag;
	  if (!batch_flag) cout << endl;
	}

	// Output progress information for the calculation.
#define MINUS_SIGN_CHARACTER 2
	if (!batch_flag) Output_Progress_L2norm(number_of_time_steps,
						Time*THOUSAND,
						total_cpu_time,
						residual_l2_norm,
						ratio_residual_l2_norm,
						first_step,
						Input_Parameters.Output_Progress_Frequency,
						MINUS_SIGN_CHARACTER);

	if (CFFC_Primary_MPI_Processor() && !first_step)
	  Output_Progress_to_File(residual_file,
				  number_of_time_steps,
				  Time*THOUSAND,
				  total_cpu_time,
				  residual_l1_norm,
				  residual_l2_norm,
				  residual_max_norm);

	if (Input_Parameters.Explicit_Write_Output_Freq > 0 &&
	    (number_of_time_steps % Input_Parameters.Explicit_Write_Output_Freq == 0)) {
	  if (!batch_flag) { 
			cout << endl << " Writing out solution in Tecplot format at time step ";
			cout << number_of_time_steps << "." << endl; 
		}
		Output_Tecplot(Local_SolnBlk, List_of_Local_Solution_Blocks,
				Input_Parameters, number_of_time_steps, Time, 
				true, residual_l2_norm, ratio_residual_l2_norm);
		Output_Cells_Tecplot(Local_SolnBlk, List_of_Local_Solution_Blocks,
				Input_Parameters, number_of_time_steps, Time, 
				true, residual_l2_norm, ratio_residual_l2_norm);
	}



	

	// Check to see if calculations are complete.
	if (!Input_Parameters.Time_Accurate &&
	    number_of_time_steps >= 
	    Input_Parameters.Maximum_Number_of_Time_Steps) break;
	if (Input_Parameters.Time_Accurate &&
	    Time >= Input_Parameters.Time_Max) break;
	if (!Input_Parameters.Time_Accurate &&
			ratio_residual_l2_norm < Input_Parameters.Explicit_Rel_Tol) {
		if (!batch_flag) { 
			int tempp = cout.precision(); cout.precision(2);
			cout.setf(ios::scientific);
			cout << "\nRelaxation Explicit Time Stepping met relative tolerance\n";
			cout << " (" << ratio_residual_l2_norm << " < ";
			cout << Input_Parameters.Explicit_Rel_Tol << ").\n";
			cout.unsetf(ios::scientific);
			cout.precision(tempp);
		}
		break;
	}

        // Freeze limiters as necessary.
        if (!first_step &&
            Input_Parameters.Freeze_Limiter &&
            !limiter_freezing_flag &&           
            residual_l2_norm <= Input_Parameters.Freeze_Limiter_Residual_Level) {
	  Freeze_Limiters(Local_SolnBlk,List_of_Local_Solution_Blocks);
	  limiter_freezing_flag = ON;
	}

	// Update solution for next time step using a multistage
	// time stepping scheme.
	for (int i_stage = 1; i_stage <= Input_Parameters.N_Stage; i_stage++) {

	  // Step 1. Exchange solution information between neighbouring blocks.
	  if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk,
							  List_of_Local_Solution_Blocks,
							  NUM_VAR_HIGHTEMP2D,
							  OFF);
	  if (error_flag) {
	    cout << "\n HighTemp2D ERROR: HighTemp2D message passing error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
	  }
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return error_flag;

	  // Step 2. Apply boundary conditions for stage.
  	  BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

	  // Step 3. Determine solution residuals for stage.
	  error_flag = dUdt_Multistage_Explicit(Local_SolnBlk,
						List_of_Global_Solution_Blocks,
						List_of_Local_Solution_Blocks,
						Input_Parameters,
						i_stage);
	  if (error_flag) {
	    cout << "\n HighTemp2D ERROR: HighTemp2D solution error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
	  }
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return error_flag;

	  // Step 4. Send boundary flux corrections at block interfaces with resolution changes.
	  error_flag = Send_Conservative_Flux_Corrections(Local_SolnBlk,
							  List_of_Local_Solution_Blocks,
							  NUM_VAR_HIGHTEMP2D);
	  if (error_flag) {
	    cout << "\n HighTemp2D ERROR: HighTemp2D flux correction message passing error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
	  }
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return error_flag;

	  // Step 5. Apply boundary flux corrections to ensure that method is conservative.
	  Apply_Boundary_Flux_Corrections_Multistage_Explicit(Local_SolnBlk,
							      List_of_Local_Solution_Blocks,
							      Input_Parameters,
							      i_stage);

          // Step 6. Smooth the solution residual using implicit residual smoothing.
          if (Input_Parameters.Residual_Smoothing) {
	    Residual_Smoothing(Local_SolnBlk,
			       List_of_Local_Solution_Blocks,
			       Input_Parameters,
			       i_stage);
	  }

	  // Step 7. Update solution for stage.
	  error_flag = Update_Solution_Multistage_Explicit(Local_SolnBlk,
							   List_of_Local_Solution_Blocks,
							   Input_Parameters,
							   i_stage);
	  if (error_flag) {
	    cout << "\n HighTemp2D ERROR: HighTemp2D solution update error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
	  }
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return error_flag;

	   // Step 8. Apply turbulence boundary conditions.
	  //Turbulence --> Turbulent
	  error_flag = Turbulent_BCs(Local_SolnBlk,
				      List_of_Local_Solution_Blocks,
				      Input_Parameters);
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return error_flag;
	  
	}

	first_step = 0;
	// Update time and time step counter.
	number_of_time_steps = number_of_time_steps + 1;
	if (Input_Parameters.i_Time_Integration != 
	    TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) {
	  Time = Time + Input_Parameters.CFL_Number*dTime;
	} else {
	  Time = Time + Input_Parameters.CFL_Number*dTime*
	    MultiStage_Optimally_Smoothing(Input_Parameters.N_Stage,
					   Input_Parameters.N_Stage,
					   Input_Parameters.i_Limiter);
	}

      } // while (perform_explicit_time_marching)

      if (!batch_flag) cout << "\n\n HighTemp2D computations complete on " 
			    << Date_And_Time() << "." << endl;
    } // if we should continue explicit scheme

    // MPI barrier to ensure processor synchronization.
    CFFC_Barrier_MPI();

    // Update ghostcell information and prescribe boundary conditions to 
    // ensure that the solution is consistent on each block.
    error_flag = Send_All_Messages(Local_SolnBlk,
				   List_of_Local_Solution_Blocks,
				   NUM_COMP_VECTOR2D,
				   ON);
    if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk,
						    List_of_Local_Solution_Blocks,
						    NUM_VAR_HIGHTEMP2D,
						    OFF);
    if (error_flag) {
      cout << "\n HighTemp2D ERROR: HighTemp2D message passing error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Apply boundary conditions.
    BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

    // Close residual file.
    if (CFFC_Primary_MPI_Processor()) error_flag = Close_Progress_File(residual_file);

  } // endif - multigrid or not 

  time(&end_explicit);

  // Start APPLY Newton_Krylov_Schwarz

  time(&start_NKS); 

  if (Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0) {
     processor_cpu_time.update();
     total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput);  
     double temp_t = total_cpu_time.cput;

     if (Input_Parameters.FlowType) {
	//  Where should we allocate the face gradient arrays
	//  dWd[xy]_face[NE]? Well:
	//    - they are not necessary for all equation types (Euler NKS
	//      doesn't need them for example), and
	//    - they are not necessary for High-Temp when using an
	//      explicit time stepping.
	//  So this is the best place to allocate them. Note that
	//  they are deallocated in the general deallocate() call.
	for (int nb = 0; nb < List_of_Local_Solution_Blocks.Nblk; nb++) {
	   if (List_of_Local_Solution_Blocks.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
		Local_SolnBlk[nb].allocate_face_grad_arrays();
	   }
	}
     } else {
	bool say_it = false;
	switch (Input_Parameters.NKS_IP.Jacobian_Order) {
  	  case SECOND_ORDER_DIAMOND_WITH_HLLE:
	    Input_Parameters.NKS_IP.Jacobian_Order = FIRST_ORDER_INVISCID_HLLE; 
	    say_it = true;
	    break;
	  case SECOND_ORDER_DIAMOND_WITH_ROE:
	    Input_Parameters.NKS_IP.Jacobian_Order = FIRST_ORDER_INVISCID_ROE;
	    say_it = true;
	  break;
	}
	if (say_it && CFFC_Primary_MPI_Processor()) {
	   cout << "\n Warning: viscous preconditioner Jacobian requested even though flow is inviscid.";
	   cout << "\n Changing Jacobian to only evaluate inviscid terms.";
	}
     }

     if (CFFC_Primary_MPI_Processor()) {
	error_flag = Open_Progress_File(residual_file,
				        Input_Parameters.Output_File_Name,
					number_of_time_steps);
        if (error_flag) {
  	   cout << "\n HighTemp2D ERROR: Unable to open residual file "
		<< "for HighTemp2D calculation.\n";
	   cout.flush();
	} 
     }

     CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
     CFFC_Broadcast_MPI(&error_flag, 1);
     if (error_flag) return (error_flag);

     // Turn off limiter freezing.
     Evaluate_Limiters(Local_SolnBlk, List_of_Local_Solution_Blocks);

     if (!batch_flag) {
	cout << "\n\n Beginning HighTemp2D NKS computations on ";
	cout << Date_And_Time() << ".\n\n";
     }

     error_flag = Newton_Krylov_Schwarz_Solver<HighTemp2D_pState,
					       HighTemp2D_Quad_Block,
					       HighTemp2D_Input_Parameters>(
			processor_cpu_time,
			residual_file,
			number_of_time_steps, // For printing. Set to total steps by NKS.
			Local_SolnBlk, 
			List_of_Local_Solution_Blocks,
			Input_Parameters);

     processor_cpu_time.update();
     total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput);  
     NKS_total_cpu_time.cput += total_cpu_time.cput - temp_t;

     if (error_flag) {
        if (CFFC_Primary_MPI_Processor()) { 
  	   cout << "\n HighTemp2D_NKS ERROR: HighTemp2D solution error on processor " 
		<< List_of_Local_Solution_Blocks.ThisCPU << ".\n";
	   cout.flush();
	} 
     }

     CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
     CFFC_Broadcast_MPI(&error_flag, 1);
     if (error_flag) return (error_flag);

     /***********************************************************************/
     if (!batch_flag) { 
	cout << "\n\n HighTemp2D NKS computations complete on " 
 	     << Date_And_Time() << ".\n"; 
     }
     if (CFFC_Primary_MPI_Processor()) error_flag = Close_Progress_File(residual_file);
  } 

  time(&end_NKS); 

  // End of the Newton-Krylov-Schwarz solver. 
  // Perhaps we should print the memory usage of the NKS code?

  if (!batch_flag) { 
     int tmpp = cout.precision();
     cout.precision(5);

     cout << "\n |-----------------------------------------------------------";
     cout << "\n |  ";
     cout << "\n |   Solution Computational Timings:";
     cout << "\n |  ";
     cout << "\n |   The CPU times are summed from all the processors";
     cout << "\n |   and include all the time spent on this test case";
     cout << "\n |   including time on pre-refined grids for cases with";
     cout << "\n |   AMR and time spent before a restart for cases";
     cout << "\n |   started from a restart file. ";
     cout << "\n |  ";
     cout << "\n |   The clock times are for the current mesh and";
     cout << "\n |   current restart only.";
     cout << "\n |  ";
     if (Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0) {
        cout << "\n |   The NKS CPU time is the total time spent doing NKS";
        cout << "\n |   calculations on all the meshes for cases with AMR. ";
        cout << "\n |  ";
        cout << "\n |   To do: incorporate NKS CPU time into the restart file.";
        cout << "\n |  ";
     }
     cout << "\n |-----------------------------------------------------------";

     if (Input_Parameters.i_ICs != IC_RESTART &&
	 Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0) {
	cout <<"\n Startup CPU Time  = " << setw(8) 
             << total_cpu_time.min() - NKS_total_cpu_time.min() << "  minutes";
	cout <<"\n NKS CPU Time      = " << setw(8) 
             << NKS_total_cpu_time.min() << "  minutes";
     }
     cout <<"\n Total CPU Time    = " << setw(8) << total_cpu_time.min() << "  minutes";

     if (Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0) {
	cout <<"\n Startup Clock Time= " << setw(8) 
             << difftime(end_explicit,start_explicit)/60.0 << "  minutes";
	cout <<"\n NKS Clock Time    = " << setw(8) 
             << difftime(end_NKS,start_NKS)/60.0 << "  minutes";
     }
     cout <<"\n Total Clock Time  = " << setw(8) 
          << difftime(end_NKS,start_explicit)/60.0 << "  minutes";

     cout<<"\n ----------------------------------------------------------------";
     cout<<"\n ----------------------------------------------------------------";
     cout<<"\n ----------------------------------------------------------------\n";
     cout.precision(tmpp);
  } 

  /********************************************************************
   * Solution calculations complete.  Write 2D High-Temp solution *
   * to output and restart files as required, reset solution          *
   * parameters, and run other cases specified by input parameters.   *
   ********************************************************************/

postprocess_current_calculation: ;

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // the dummy_index is to elminate compiler warnings of "unreachable code".
  for (int dummy_index = 0, a_big_number = 10000; dummy_index < a_big_number; dummy_index++) {
    if (CFFC_Primary_MPI_Processor()) {
      Get_Next_Input_Control_Parameter(Input_Parameters);
      command_flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
      line_number = Input_Parameters.Line_Number;
      if (command_flag == INVALID_INPUT_CODE ||
	  command_flag == INVALID_INPUT_VALUE) {
	line_number = -line_number;
	cout << "\n HighTemp2D ERROR: Error reading HighTemp2D data at line #"
	     << -line_number << " of input data file.\n";
	return 1;
      }
      Reinitialize_Reference_State(Input_Parameters);
    }
    // MPI barrier to ensure processor synchronization.
    CFFC_Barrier_MPI(); 
    Broadcast_Input_Parameters(Input_Parameters);
    CFFC_Broadcast_MPI(&command_flag,1);
    for (int nb = 0; nb < List_of_Local_Solution_Blocks.Nblk; nb++) {
      if (List_of_Local_Solution_Blocks.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	// Set flow type indicator (inviscid/viscous).
	Local_SolnBlk[nb].Flow_Type = Input_Parameters.FlowType;
	// Set flow geometry indicator (planar/axisymmetric) flow block.
	Local_SolnBlk[nb].Axisymmetric = Input_Parameters.Axisymmetric;
      }
    }

    if (command_flag == EXECUTE_CODE) {
      // Deallocate memory for 2D HighTemp equation solution.
      if (!batch_flag) cout << "\n Deallocating HighTemp2D solution variables.";

			if (mgsolver_is_allocated) {
      if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	  Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
	MGSolver.deallocate();
		mgsolver_is_allocated = false;
      }
			}

      Local_SolnBlk = Deallocate(Local_SolnBlk,Input_Parameters);
      //Deallocate_Message_Buffers(List_of_Local_Solution_Blocks); // Not necessary here!
      List_of_Local_Solution_Blocks.deallocate();
      List_of_Global_Solution_Blocks.deallocate();
      QuadTree.deallocate();
      MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk,
					    Input_Parameters.Number_of_Blocks_Idir,
					    Input_Parameters.Number_of_Blocks_Jdir);
      // Output input parameters for new caluculation.
      if (!batch_flag) cout << "\n\n Starting a new calculation." << Input_Parameters << "\n";
      // Execute new calculation.
      goto execute_new_calculation;

    } else if (command_flag == TERMINATE_CODE) {

      CFFC_Barrier_MPI();
      // Deallocate memory for 2D HighTemp equation solution.
      if (!batch_flag) cout << "\n Deallocating HighTemp2D solution variables.";
      if (mgsolver_is_allocated) {
         if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	     Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
	    MGSolver.deallocate();
	    mgsolver_is_allocated = false;
         }
      }
      Local_SolnBlk = Deallocate(Local_SolnBlk,Input_Parameters);
      //Deallocate_Message_Buffers(List_of_Local_Solution_Blocks); // Not necessary here!
      List_of_Local_Solution_Blocks.deallocate();
      List_of_Global_Solution_Blocks.deallocate();
      QuadTree.deallocate();
      MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk,
					    Input_Parameters.Number_of_Blocks_Idir,
					    Input_Parameters.Number_of_Blocks_Jdir);
      // Close input data file.
      if (!batch_flag) cout << "\n\n Closing HighTemp2D input data file.";
      if (CFFC_Primary_MPI_Processor()) Close_Input_File(Input_Parameters);
      // Terminate calculation.
      return 0;

    } else if (command_flag == CONTINUE_CODE) {
      //  Bump the Max Time Steps by the number of time steps run on the
      //  previous grid. Think about it. Before starting the previous grid, we
      //  had run "x" steps and so Maximum_Number_of_Time_Steps was "x" +
      //  "m0", m0 the max on the first grid. Now suppose that on this
      //  previous grid we actually ran "y" steps so that now
      //  number_of_time_steps is "x" + "y". Clearly for the next grid we want
      //  Maximum_Number_of_Time_Steps to be "x" + "y" + "m0".
      //
      //    -- Alistair Wood Thu Mar 29 2007 
      Input_Parameters.Maximum_Number_of_Time_Steps += 
      number_of_time_steps - start_number_of_time_steps;
      // Output input parameters for continuing calculation.
      if (!batch_flag) cout << "\n\n Continuing existing calculation."
			    << Input_Parameters << "\n";
      // Deallocate multigrid if necessary.
      if (mgsolver_is_allocated) {
         if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	     Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
	    MGSolver.deallocate();
	    mgsolver_is_allocated = false;
	 }
      }
      // Reinitialize the initial conditions for a wedge flow.
      error_flag = Reinitialize_Wedge_Initial_Conditions(Local_SolnBlk,
							 List_of_Local_Solution_Blocks,
							 Input_Parameters);
      if (error_flag) {
	cout << "\n HighTemp2D ERROR: HighTemp2D reinitializing wedge flow initial"
	     << "\n conditions on processor " << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;
      // Continue existing calculation.
      goto continue_existing_calculation;

    } else if (command_flag == REFINE_GRID_CODE) {
      // Refine mesh using block based adaptive mesh refinement algorithm.
      if (!batch_flag) cout << "\n Refining Grid.  Performing adaptive mesh refinement.";
      Evaluate_Limiters(Local_SolnBlk,List_of_Local_Solution_Blocks);
      error_flag = Adaptive_Mesh_Refinement(Local_SolnBlk,
					    Input_Parameters,
					    QuadTree,
					    List_of_Global_Solution_Blocks,
					    List_of_Local_Solution_Blocks);
      if (error_flag) {
	cout << "\n HighTemp2D ERROR: HighTemp2D AMR error on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".  Error number = " << error_flag << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) {
	command_flag = Output_Tecplot(Local_SolnBlk,
				      List_of_Local_Solution_Blocks,
				      Input_Parameters,
				      number_of_time_steps,
				      Time);
	return error_flag;
      }

      // Output multi-block solution-adaptive quadrilateral mesh statistics.
      if (!batch_flag) {
	cout << "\n\n New multi-block solution-adaptive quadrilateral mesh statistics: ";
	cout << "\n  -> Number of Root Blocks i-direction: "
	     << QuadTree.NRi;
	cout << "\n  -> Number of Root Blocks j-direction: "
	     << QuadTree.NRj;
	cout << "\n  -> Total Number of Used Blocks: "
	     << QuadTree.countUsedBlocks();
	cout << "\n  -> Total Number of Computational Cells: "
	     << QuadTree.countUsedCells();
	cout << "\n  -> Number of Mesh Refinement Levels: "
	     << QuadTree.highestRefinementLevel()+1;
	cout << "\n  -> Refinement Efficiency: "
	     << QuadTree.efficiencyRefinement() << endl;
      }
      // Apply boundary conditions.
      BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

    } else if (command_flag == MORTON_ORDERING_CODE) {
      if (!batch_flag) { cout << "\n\n Applying Morton re-ordering algorithm. "; }
      error_flag = Morton_ReOrdering_of_Solution_Blocks(QuadTree, 
					                List_of_Global_Solution_Blocks, 
					                List_of_Local_Solution_Blocks, 
					                Local_SolnBlk, 
					                Input_Parameters, 
					                number_of_time_steps, 
					                Time, 
					                processor_cpu_time); 
      if (error_flag) {
	cout <<"\n HighTemp2D ERROR: Morton re-ordering error on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
      } 
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) { return error_flag; }
      //  The Morton Re-ordering writes out and then reads back in the
      //  solution using restart files. The function which reads the
      //  restart files makes the assumption that we are starting a new
      //  calculation with restart initial conditions and so does:
      //  
      //    restart_file >> Number_of_Time_Steps;
      //    IP.Maximum_Number_of_Time_Steps += Number_of_Time_Steps;
      //  
      //  which we need to correct here.
      //  
      //  Is this maximum number of time steps the best way to do it? Any
      //  better way would have to keep in mind Morton, steady-state
      //  problems on one grid with multiple restarts from restart files,
      //  steady-state problems with successively refined grids, NKS,
      //  time-accurate problems and so on. Maybe someone else can think
      //  about this.
      //     -- Alistair Wood Fri Jun 15 2007 
      Input_Parameters.Maximum_Number_of_Time_Steps -= number_of_time_steps;

      if (!batch_flag) { cout << "\n Outputting space filling curve showing block loading for CPUs."; }
      Morton_SFC_Output_Tecplot(Local_SolnBlk, 
				Input_Parameters, 
				List_of_Local_Solution_Blocks);

    } else if (command_flag == WRITE_OUTPUT_CODE) {
	bool output_multigrid = false;
	if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	    Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
	   output_multigrid = true;
	}
	// but if we ran NKS then multigrid was just a startup 
	// and so we don't care about it.
	if (Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0) {
	   output_multigrid = false;
	}
	if (!batch_flag) {
	   cout << "\n Writing HighTemp2D solution (at the cell nodes) ";
	   cout << "to output data file(s).";
	} 
	if (output_multigrid) {
	   error_flag = MGSolver.Output_Multigrid(number_of_time_steps, Time);
	} else {
	   error_flag = Output_Tecplot(Local_SolnBlk,
		  		       List_of_Local_Solution_Blocks,
				       Input_Parameters,
				       number_of_time_steps,
				       Time);
        }
        if (error_flag) {
	   cout << "\n HighTemp2D ERROR: Unable to open HighTemp2D "
		<< "output data file(s) on processor "
	        << List_of_Local_Solution_Blocks.ThisCPU
	        << "." << endl;
        }
        error_flag = CFFC_OR_MPI(error_flag);
        if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_CELLS_CODE) {
	bool output_multigrid = false;
	if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	    Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
	   output_multigrid = true;
	}
	// but if we ran NKS then multigrid was just a startup 
	// and so we don't care about it.
	if (Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0) {
	   output_multigrid = false;
	}
	if (!batch_flag) {
	   cout << "\n Writing cell-centered HighTemp2D solution ";
	   cout << "to output data file(s).";
	} 
	if (output_multigrid) {
	   error_flag = MGSolver.Output_Multigrid_Cells(number_of_time_steps, Time);
	} else {
	error_flag = Output_Cells_Tecplot(Local_SolnBlk,
					  List_of_Local_Solution_Blocks,
					  Input_Parameters,
					  number_of_time_steps,
					  Time);
        }
        if (error_flag) {
	   cout << "\n HighTemp2D ERROR: Unable to open HighTemp2D "
		<< "output data file(s) on processor "
	        << List_of_Local_Solution_Blocks.ThisCPU
	        << "." << endl;
        }
        error_flag = CFFC_OR_MPI(error_flag);
        if (error_flag) return error_flag;

      /*  
    } else if (command_flag == WRITE_OUTPUT_MULTIGRID_CODE) {
      // Output multi-block solution-adaptive multigrid nodes.
      if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	  Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
	if (!batch_flag) cout << "\n Writing HighTemp2D multigrid nodes to output file.";
	error_flag = MGSolver.Output_Multigrid(number_of_time_steps,
					       Time);
	if (error_flag) {
	  cout << "\n HighTemp2D ERROR: Unable to open HighTemp2D multigrid nodes output file.\n";
	  cout.flush();
	}
	CFFC_Broadcast_MPI(&error_flag, 1);
	if (error_flag) return (error_flag);
      }
      
    } else if (command_flag == WRITE_OUTPUT_MULTIGRID_CELLS_CODE) {
      // Output multi-block solution-adaptive multigrid cells.
      if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	  Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
	if (!batch_flag) cout << "\n Writing HighTemp2D multigrid cells to output file.";
	error_flag = MGSolver.Output_Multigrid_Cells(number_of_time_steps,
						     Time);
	if (error_flag) {
	  cout << "\n HighTemp2D ERROR: Unable to open HighTemp2D multigrid cells output file.\n";
	  cout.flush();
	}
	CFFC_Broadcast_MPI(&error_flag, 1);
	if (error_flag) return (error_flag);
      }
      */
    } else if (command_flag == WRITE_OUTPUT_NODES_CODE) {
	bool output_multigrid = false;
	if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	    Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
	   output_multigrid = true;
	}
	// but if we ran NKS then multigrid was just a startup 
	// and so we don't care about it.
	if (Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0) {
	   output_multigrid = false;
	}
	if (!batch_flag) {
	   cout << "\n Writing HighTemp2D node locations ";
	   cout << "to output data file(s).";
	} 
	if (output_multigrid) {
	   error_flag = MGSolver.Output_Multigrid_Nodes(number_of_time_steps, Time);
	} else {
	   error_flag = Output_Nodes_Tecplot(Local_SolnBlk,
		  			     List_of_Local_Solution_Blocks,
					     Input_Parameters,
					     number_of_time_steps,
					     Time);
        }
        if (error_flag) {
	   cout << "\n HighTemp2D ERROR: Unable to open HighTemp2D "
		<< "output data file(s) on processor "
	        << List_of_Local_Solution_Blocks.ThisCPU
	        << "." << endl;
        }
        error_flag = CFFC_OR_MPI(error_flag);
        if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_GRADIENTS_CODE) {
      if (!batch_flag) cout << "\n Writing HighTemp2D node locations to output data file(s).";
      error_flag = Output_Gradients_Tecplot(Local_SolnBlk,
					    List_of_Local_Solution_Blocks,
					    Input_Parameters,
					    number_of_time_steps,
					    Time);
      if (error_flag) {
	cout << "\n HighTemp2D ERROR: Unable to open HighTemp2D nodes output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;  

    } else if (command_flag == WRITE_OUTPUT_QUASI3D_CODE) {
      // Output solution data.
      if (!batch_flag) cout << "\n Writing HighTemp2D quasi3D solution to output data file(s).";
      error_flag = Output_Quasi3D_Tecplot(Local_SolnBlk,
					  List_of_Local_Solution_Blocks,
					  Input_Parameters,
					  number_of_time_steps,
					  Time);
      if (error_flag) {
	cout << "\n HighTemp2D ERROR: Unable to open HighTemp2D quasi3D output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_RESTART_CODE) {
      // Write restart files.
      if (!batch_flag) cout << "\n Writing HighTemp2D solution to restart data file(s).";
      // Write the quadtree restart file.
      error_flag = Write_QuadTree(QuadTree,Input_Parameters);
      if (error_flag) {
	cout << "\n HighTemp2D ERROR: Unable to open HighTemp2D quadtree data file "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;
      // Write the solution block restart files.
      error_flag = Write_Restart_Solution(Local_SolnBlk,
					  List_of_Local_Solution_Blocks,
					  Input_Parameters,
					  number_of_time_steps,
					  Time,
					  processor_cpu_time);
      if (error_flag) {
	cout << "\n HighTemp2D ERROR: Unable to open HighTemp2D restart output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_GRID_CODE) {
      // Output multi-block solution-adaptive mesh data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing HighTemp2D multi-block mesh to grid data output file.";
	error_flag = Output_Tecplot(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << "\n HighTemp2D ERROR: Unable to open HighTemp2D mesh data output file.\n";
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      
    } else if (command_flag == WRITE_GRID_DEFINITION_CODE) {
      // Write multi-block solution-adaptive mesh definition files.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing HighTemp2D multi-block mesh to grid definition files.";
	error_flag = Write_Multi_Block_Grid_Definition(MeshBlk,
						       Input_Parameters);
	error_flag = Write_Multi_Block_Grid(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << "\n HighTemp2D ERROR: Unable to open HighTemp2D multi-block mesh definition files.\n";
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      
    } else if (command_flag == WRITE_OUTPUT_GRID_NODES_CODE) {
      // Output multi-block solution-adaptive mesh node data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing HighTemp2D multi-block mesh to node data output file.";
	error_flag = Output_Nodes_Tecplot(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << "\n HighTemp2D ERROR: Unable to open HighTemp2D mesh node data output file.\n";
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      
    } else if (command_flag == WRITE_OUTPUT_GRID_CELLS_CODE) {
      // Output multi-block solution-adaptive mesh cell data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing HighTemp2D multi-block mesh to cell data output file.";
	error_flag = Output_Cells_Tecplot(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << "\n HighTemp2D ERROR: Unable to open HighTemp2D mesh cell data output file.\n";
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_RINGLEB_CODE) {
      if (!batch_flag) cout << endl << " Writing exact solution and error norms for Ringleb's flow.";
      error_flag = Output_Ringleb_Tecplot(Local_SolnBlk,
					  List_of_Local_Solution_Blocks,
					  Input_Parameters);
      if (error_flag) {
	cout << endl << "\n HighTemp2D ERROR: Unable to open HighTemp2D Ringleb's flow output file." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_VISCOUS_CHANNEL_CODE) {
      if (!batch_flag) cout << endl << " Writing exact solution and error norms for the viscous channel flow.";
      error_flag = Output_Viscous_Channel_Tecplot(Local_SolnBlk,
						  List_of_Local_Solution_Blocks,
						  Input_Parameters);
      if (error_flag) {
	cout << endl << "\n HighTemp2D ERROR: Unable to open HighTemp2D viscous channel flow output file." << endl;
	cout.flush();
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_VISCOUS_PIPE_CODE) {
      if (!batch_flag) cout << endl << " Writing exact solution and error norms for the viscous pipe flow.";
      error_flag = Output_Viscous_Pipe_Tecplot(Local_SolnBlk,
					       List_of_Local_Solution_Blocks,
					       Input_Parameters);
      if (error_flag) {
	cout << endl << "\n HighTemp2D ERROR: Unable to open HighTemp2D viscous pipe flow output file." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_TURBULENT_PIPE_CODE) {
      if (!batch_flag) cout << endl 
                            << " Writing experimental and computed solution data for the turbulent pipe flow.";
      error_flag = Output_Turbulent_Pipe_Tecplot(Local_SolnBlk,
						 List_of_Local_Solution_Blocks,
						 Input_Parameters);
      if (error_flag) {
	cout << endl << "\n HighTemp2D ERROR: Unable to open HighTemp2D viscous pipe flow output file." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_FLAT_PLATE_CODE) {
      if (!batch_flag) cout << endl << " Writing exact solution and error norms for the flat plate flow (Blasius solution).";
      error_flag = Output_Flat_Plate_Tecplot(Local_SolnBlk,
					     List_of_Local_Solution_Blocks,
					     Input_Parameters);
      if (error_flag) {
	cout << endl << "\n HighTemp2D ERROR: Unable to open HighTemp2D flat plate output file." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_DRIVEN_CAVITY_FLOW_CODE) {
      if (!batch_flag) cout << endl << " Writing the driven cavity flow output file.";
      error_flag = Output_Driven_Cavity_Flow_Tecplot(Local_SolnBlk,
						     List_of_Local_Solution_Blocks,
						     Input_Parameters);
      if (error_flag) {
	cout << endl << "\n HighTemp2D ERROR: Unable to open HighTemp2D driven cavity flow output file." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_BACKWARD_FACING_STEP_CODE) {
      if (!batch_flag) cout << endl << " Writing the backward facing step output file.";
      error_flag = Output_Backward_Facing_Step_Tecplot(Local_SolnBlk,
						       List_of_Local_Solution_Blocks,
						       Input_Parameters);
      if (error_flag) {
	cout << endl << "\n HighTemp2D ERROR: Unable to open HighTemp2D backward facing step output file." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_FORWARD_FACING_STEP_CODE) {
      if (!batch_flag) cout << endl << " Writing the forward facing step output file.";
      error_flag = Output_Forward_Facing_Step_Tecplot(Local_SolnBlk,
						       List_of_Local_Solution_Blocks,
						       Input_Parameters);
      if (error_flag) {
	cout << endl << "\n HighTemp2D ERROR: Unable to open HighTemp2D forward facing step output file." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_ELLIPTIC_OPERATOR_ANALYSIS) {
      if (!batch_flag) cout << endl << " Performing the elliptic operator analyses:";
      if (!batch_flag) cout << endl << "  -> Face gradient reconstruction using the average cell-centred"
			    << endl << "     gradient constructed with linear least-squares reconstruction.";
      error_flag = Elliptic_Operator_Analysis_Average_Gradient_Linear_Least_Squares(Local_SolnBlk,
										    List_of_Local_Solution_Blocks,
										    Input_Parameters,
										    OFF);
      if (error_flag) {
	cout << endl << "\n HighTemp2D ERROR: Unable to perform the elliptic operator analysis." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      if (!batch_flag) cout << endl << "  -> Face gradient reconstruction using a directional derivative.";
      error_flag = Elliptic_Operator_Analysis_Directional_Derivative(Local_SolnBlk,
 								     List_of_Local_Solution_Blocks,
 								     Input_Parameters);
      if (error_flag) {
	cout << endl << "\n HighTemp2D ERROR: Unable to perform the elliptic operator analysis." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      if (!batch_flag) cout << endl << "  -> Face gradient reconstruction using Green-Gauss reconstruction over"
			    << endl << "     the diamond-path with simple vertex weighting.";
      error_flag = Elliptic_Operator_Analysis_Diamond_Path(Local_SolnBlk,
 							   List_of_Local_Solution_Blocks,
 							   Input_Parameters,
 							   0);
      if (error_flag) {
	cout << endl << "\n HighTemp2D ERROR: Unable to perform the elliptic operator analysis." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      if (!batch_flag) cout << endl << "  -> Face gradient reconstruction using Green-Gauss reconstruction over"
			    << endl << "     the diamond-path with the linearity-preserving weighting function"
			    << endl << "     of Holmes and Connell.";
      error_flag = Elliptic_Operator_Analysis_Diamond_Path(Local_SolnBlk,
 							   List_of_Local_Solution_Blocks,
 							   Input_Parameters,
 							   1);
      if (error_flag) {
	cout << endl << "\n HighTemp2D ERROR: Unable to perform the elliptic operator analysis." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      if (!batch_flag) cout << endl << "  -> Face gradient reconstruction using Green-Gauss reconstruction over"
			    << endl << "     the diamond-path with the linearity-preserving weighting function"
			    << endl << "     of Zingg and Yarrow.";
      error_flag = Elliptic_Operator_Analysis_Diamond_Path(Local_SolnBlk,
 							   List_of_Local_Solution_Blocks,
 							   Input_Parameters,
 							   2);
      if (error_flag) {
	cout << endl << "\n HighTemp2D ERROR: Unable to perform the elliptic operator analysis." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      if (!batch_flag) cout << endl << "  -> Face gradient reconstruction using the corrected average cell-"
			    << endl << "     centred gradient constructed with linear least-squares reconstruction";
      error_flag = Elliptic_Operator_Analysis_Average_Gradient_Linear_Least_Squares(Local_SolnBlk,
										    List_of_Local_Solution_Blocks,
										    Input_Parameters,
										    ON);
      if (error_flag) {
	cout << endl << "\n HighTemp2D ERROR: Unable to perform the elliptic operator analysis." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == INVALID_INPUT_CODE ||
	       command_flag == INVALID_INPUT_VALUE) {
      line_number = -line_number;
      cout << "\n HighTemp2D ERROR: Error reading HighTemp2D data at line #"
	   << -line_number  << " of input data file.\n";
      return line_number;
    }
    
  }
  
  /********************************************************************
   * End of all HighTemp2DSolver computations and I/O.                *
   ********************************************************************/
  return 0;
  
}
