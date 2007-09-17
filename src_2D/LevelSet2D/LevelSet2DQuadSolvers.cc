/**********************************************************************
 * LevelSet2DQuadSolvers.cc                                           *
 *                                                                    *
 * 2D Level Set equation multi-block quadrilateral mesh solvers.      *
 *                                                                    *
 **********************************************************************/

// Include 2D LevelSet quadrilateral mesh solution header file.

#ifndef _LEVELSET2D_QUAD_INCLUDED
#include "LevelSet2DQuad.h"
#endif // _LEVELSET2D_QUAD_INCLUDED

/**********************************************************************
 * Routine: LevelSet2DQuadSolver                                      *
 *                                                                    *
 * Computes solutions to 2D LevelSet equations on 2D quadrilateral    *
 * multi-block solution-adaptive mesh.                                *
 *                                                                    *
 **********************************************************************/
int LevelSet2DQuadSolver(char *Input_File_Name_ptr, int batch_flag) {

  /********************************************************************
   * Local variable declarations.                                     *
   ********************************************************************/
  
  // LevelSet2D input variables and parameters:
  LevelSet2D_Input_Parameters Input_Parameters;
  
  // Multi-block solution-adaptive quadrilateral mesh solution variables.
  Grid2D_Quad_Block           **MeshBlk;
  QuadTreeBlock_DataStructure   QuadTree;
  AdaptiveBlockResourceList     List_of_Global_Solution_Blocks;
  AdaptiveBlock2D_List          List_of_Local_Solution_Blocks;
  LevelSet2D_Quad_Block        *Local_SolnBlk;

  // Define residual file and cpu time variables.
  ofstream residual_file;
  CPUTime processor_cpu_time, total_cpu_time;
  
  // Other local solution variables.
  int number_of_time_steps, first_step,
      command_flag, error_flag, line_number, 
      perform_explicit_time_marching, limiter_freezing_off,
      redistance_count;

  double global_error, global_area, weighted_global_error;
  global_error = ZERO;
  global_area = ZERO;
  weighted_global_error = ZERO;

  double Time, dTime;

  double residual_l2norm_first, residual_ratio;
  double residual_l1_norm, residual_l2_norm, residual_max_norm;

  /********************************************************************
   * Set default values for the input solution parameters and then    *
   * read user specified input values from the specified input        *
   * parameter file.                                                  *
   ********************************************************************/

  // The primary MPI processor processes the input parameter file.
  if (CFFC_Primary_MPI_Processor()) {
    if (!batch_flag) cout << "\n Reading LevelSet2D input data file `"
			  << Input_File_Name_ptr << "'.";
    // Process input file.
    error_flag = Process_Input_Control_Parameter_File(Input_Parameters,
						      Input_File_Name_ptr,
						      command_flag);
    // Output input data to screen.
    if (!batch_flag && error_flag == 0) cout << Input_Parameters << "\n";
  } else {
    error_flag = 0;
  }

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Broadcast input solution parameters to other MPI processors.
  CFFC_Broadcast_MPI(&error_flag,1);
  if (error_flag != 0) return error_flag;
  CFFC_Broadcast_MPI(&command_flag,1);
  if (command_flag == TERMINATE_CODE) return 0;
  Broadcast_Input_Parameters(Input_Parameters);

  /********************************************************************
   * Create initial mesh and allocate LevelSet2D solution variables   *
   * specified IBVP/BVP problem.                                      *
   ********************************************************************/

 execute_new_calculation: ;
  
  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Create initial mesh.  Read mesh from grid definition or data files 
  // when specified by input parameters.

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
      cout << "\n LevelSet2D ERROR: Unable to create valid LevelSet2D multi-block mesh." << endl;
    }
  } else {
    MeshBlk = NULL;
  }

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Broadcast the mesh to other MPI processors.
  CFFC_Broadcast_MPI(&error_flag,1);
  if (error_flag) return error_flag;
  MeshBlk = Broadcast_Multi_Block_Grid(MeshBlk,Input_Parameters);

  // Create (allocate) multi-block quadtree data structure, create
  // (allocate) array of local 2D LevelSet equation solution blocks,
  // assign and create (allocate) 2D LevelSet equation solution blocks
  // corresponding to the initial mesh.
  if (!batch_flag) cout << "\n Creating multi-block quadtree data structure and assigning"
                        << "\n  LevelSet2D solution blocks corresponding to initial mesh.";
  Local_SolnBlk = Create_Initial_Solution_Blocks(MeshBlk,
						 Local_SolnBlk,
						 Input_Parameters,
						 QuadTree,
						 List_of_Global_Solution_Blocks,
						 List_of_Local_Solution_Blocks);
  if (Local_SolnBlk == NULL) return 1;

  /********************************************************************
   * Initialize LevelSet2D solution variables.                        *
   ********************************************************************/

  // Set the initial time level.
  Time = ZERO;
  number_of_time_steps = 0;

  // Set the CPU time to zero.
  processor_cpu_time.zero();
  total_cpu_time.zero();

  // Initialize the conserved and primitive state solution variables.
  if (!batch_flag) cout << "\n Prescribing LevelSet2D initial data.";
  if (Input_Parameters.Interface_IP.Component_List[1].Type == INTERFACE_RESTART) {
    if (!batch_flag) cout << "\n Reading LevelSet2D solution from restart data files.";
    // Read the quadtree restart file.
    error_flag = Read_QuadTree(QuadTree,
			       List_of_Global_Solution_Blocks,
			       List_of_Local_Solution_Blocks,
			       Input_Parameters);
    if (error_flag) {
      cout << "\n LevelSet2D ERROR: Unable to open LevelSet2D quadtree data file "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;
    // Allocate the message buffers.
    Allocate_Message_Buffers(List_of_Local_Solution_Blocks,
			     Local_SolnBlk[0].NumVar()+NUM_COMP_VECTOR2D);
    // Read the solution block restart files.
    error_flag = Read_Restart_Solution(Local_SolnBlk,
				       List_of_Local_Solution_Blocks,
				       Input_Parameters,
				       number_of_time_steps,
				       Time,
				       processor_cpu_time);
    if (error_flag) {
      cout << "\n LevelSet2D ERROR: Unable to open LevelSet2D restart input data file(s) "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;
    // Ensure each processor has the correct time and time.
    number_of_time_steps = CFFC_Maximum_MPI(number_of_time_steps);
    Time = CFFC_Maximum_MPI(Time);
    processor_cpu_time.cput = CFFC_Maximum_MPI(processor_cpu_time.cput);
    // MPI barrier to ensure processor synchronization.
    CFFC_Barrier_MPI();
    // Broadcast input solution parameters to other MPI processors.
    CFFC_Broadcast_MPI(&error_flag,1);
    if (error_flag != 0) return error_flag;
    CFFC_Broadcast_MPI(&command_flag,1);
    if (command_flag == TERMINATE_CODE) return 0;
    Broadcast_Input_Parameters(Input_Parameters);

  } else {

    // Initialize interface(s).
    error_flag = Initialize_Interfaces(Local_SolnBlk,
				       List_of_Local_Solution_Blocks,
				       Input_Parameters);
    if (error_flag) {
      cout << "\n LevelSet2D ERROR: Message passing error during LevelSet2D interface"
	   << "\n                   intialization on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Construct the bulk flow-field.
    error_flag = Construct_Bulk_Flow_Field(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);
    if (error_flag) {
      cout << "\n LevelSet2D ERROR: Message passing error during LevelSet2D bulk"
	   << "\n                   flow-field intialization on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

  }

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Send solution information between neighbouring blocks to complete
  // prescription of initial data.
  error_flag = Send_All_Messages(Local_SolnBlk,
                                 List_of_Local_Solution_Blocks,
                                 NUM_COMP_VECTOR2D,
                                 ON);
  if (error_flag) {
    cout << "\n LevelSet2D ERROR: Message passing error during LevelSet2D solution"
	 << "\n intialization on processor " << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Solve extension problems if not restart.
  if (Input_Parameters.Interface_IP.Component_List[1].Type != INTERFACE_RESTART) {

    // Solve the initial extension problem.
    if (!batch_flag) cout << "\n Solving the initial extension problem ";
    switch(Input_Parameters.i_Initial_Distance_Type) {
    case LEVELSET_INITIAL_EXTENSION_GEOMETRIC :
      error_flag = Geometric_Extension_Problem(Local_SolnBlk,
					       List_of_Local_Solution_Blocks,
					       Input_Parameters);
      if (error_flag) {
	cout << "\n LevelSet2D ERROR: LevelSet2D error " << error_flag
	     << " on processor " << List_of_Local_Solution_Blocks.ThisCPU 
	     << " while computing geometric extension problem." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;
      break;
    case LEVELSET_INITIAL_EXTENSION_EXACT :
      error_flag = Exact_Initial_Extension(Local_SolnBlk,
					   List_of_Local_Solution_Blocks,
					   Input_Parameters);
      if (error_flag) {
	cout << "\n LevelSet2D ERROR: LevelSet2D error " << error_flag
	     << " on processor " << List_of_Local_Solution_Blocks.ThisCPU 
	     << " while computing exact extension problem." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;
      break;
    };

    // Redistance level set function to be a signed distance function.
    if (!batch_flag) cout << "\n Performing initial solution of the Eikonal equation.";
    error_flag = Explicit_Eikonal_Equation(Local_SolnBlk,
 					   Input_Parameters,
 					   QuadTree,
 					   List_of_Global_Solution_Blocks,
 					   List_of_Local_Solution_Blocks,
 					   OFF);
    if (error_flag) {
      cout << "\n LevelSet2D ERROR: LevelSet2D error during redistancing on processor "
 	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Solve scalar (front speed) extension equation.
    if (!batch_flag) cout << "\n Performing solution of the scalar (front speed) extension equation.";
    error_flag = Explicit_Scalar_Extension_Equation(Local_SolnBlk,
						    Input_Parameters,
						    QuadTree,
						    List_of_Global_Solution_Blocks,
						    List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n LevelSet2D ERROR: LevelSet2D error during the solution of the scalar extension euqtion on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Update ghostcell information and prescribe boundary conditions 
    // to ensure that the solution is consistent on each block.
    error_flag = Send_All_Messages(Local_SolnBlk,
                                   List_of_Local_Solution_Blocks,
                                   NUM_VAR_LEVELSET2D,
                                   OFF);
    if (error_flag) {
      cout << "\n LevelSet2D ERROR: LevelSet2D message passing error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Prescribe boundary data consistent with initial data.
    BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

    // Perform uniform mesh refinement.
    if (!batch_flag) cout << "\n Performing LevelSet2D uniform mesh refinement.";
    error_flag = Uniform_Adaptive_Mesh_Refinement(Local_SolnBlk,
						  Input_Parameters,
						  QuadTree,
						  List_of_Global_Solution_Blocks,
						  List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n LevelSet2D ERROR: Uniform AMR error #" << error_flag
	   << " on processor " << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Perform initial mesh refinement.
    if (!batch_flag) cout << "\n Performing LevelSet2D interface mesh refinement.";
    error_flag = Initial_Adaptive_Mesh_Refinement(Local_SolnBlk,
						  Input_Parameters,
						  QuadTree,
						  List_of_Global_Solution_Blocks,
						  List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n LevelSet2D ERROR: Initial AMR error #" << error_flag
	   << " on processor " << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

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

  /********************************************************************
   * Solve IBVP or BVP for conservation form of 2D Level Set          *
   * equations on multi-block solution-adaptive quadrilateral mesh.   *
   ********************************************************************/

 continue_existing_calculation: ;
  
  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  redistance_count = 0;

  // Open residual file.
  first_step = 1;
  if (CFFC_Primary_MPI_Processor()) {
    error_flag = Open_Progress_File(residual_file,
				    Input_Parameters.Output_File_Name,
				    number_of_time_steps);
    if (error_flag) {
      cout << "\n LevelSet2D ERROR: Unable to open residual file for LevelSet2D calculation.\n";
      cout.flush();
    }
  }
  CFFC_Barrier_MPI();
  CFFC_Broadcast_MPI(&error_flag,1);
  if (error_flag) return error_flag;
  // Reset the CPU time.
  processor_cpu_time.reset();
  
  // Perform required number of iterations (time steps).
  if ((!Input_Parameters.Time_Accurate &&
       Input_Parameters.Maximum_Number_of_Time_Steps > 0) ||
      (Input_Parameters.Time_Accurate &&
       Input_Parameters.Time_Max > Time)) {

    if (!batch_flag) cout << "\n Beginning LevelSet2D computations on "
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
	    (number_of_time_steps/Input_Parameters.AMR_Frequency) == 0) {
	  if (!batch_flag) cout << "\n\n Refining Grid.  Performing adaptive mesh refinement at n = "
				<< number_of_time_steps << ".";
	  error_flag = AMR(Local_SolnBlk,
			   Input_Parameters,
			   QuadTree,
			   List_of_Global_Solution_Blocks,
			   List_of_Local_Solution_Blocks,
			   ON,ON);
	  if (error_flag) {
	    cout << "\n LevelSet2D ERROR: LevelSet2D AMR error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
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
	  }
	}
      }

      // Determine local and global time steps.
      dTime = CFL_Hamilton_Jacobi(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);
      // Find global minimum time step for all processors.
      dTime = CFFC_Minimum_MPI(dTime);
      if (Input_Parameters.Time_Accurate) {
	if ((Input_Parameters.i_Time_Integration != TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) &&
	    (Time + Input_Parameters.Hamilton_Jacobi_CFL_Number*dTime > Input_Parameters.Time_Max)) {
	  dTime = (Input_Parameters.Time_Max-Time)/Input_Parameters.Hamilton_Jacobi_CFL_Number;
	} else if (Time + Input_Parameters.Hamilton_Jacobi_CFL_Number*dTime*
		   MultiStage_Optimally_Smoothing(Input_Parameters.N_Stage,
						  Input_Parameters.N_Stage,
						  Input_Parameters.i_Limiter) > 
		   Input_Parameters.Time_Max) {
	  dTime = (Input_Parameters.Time_Max-Time)/
	          (Input_Parameters.Hamilton_Jacobi_CFL_Number*MultiStage_Optimally_Smoothing(Input_Parameters.N_Stage,
									      Input_Parameters.N_Stage,
									      Input_Parameters.i_Limiter));
	}
      }
      // Set global time step.
      Set_Global_TimeStep(Local_SolnBlk,List_of_Local_Solution_Blocks,dTime);

      // Determine the L1, L2, and max norms of the solution residual.
      residual_l1_norm = L1_Norm_Residual(Local_SolnBlk,List_of_Local_Solution_Blocks,1);
      residual_l1_norm = CFFC_Summation_MPI(residual_l1_norm);

      residual_l2_norm = L2_Norm_Residual(Local_SolnBlk,List_of_Local_Solution_Blocks,1);
      residual_l2_norm = sqr(residual_l2_norm);
      residual_l2_norm = CFFC_Summation_MPI(residual_l2_norm);
      residual_l2_norm = sqrt(residual_l2_norm);

      residual_max_norm = Max_Norm_Residual(Local_SolnBlk,List_of_Local_Solution_Blocks,1);
      residual_max_norm = CFFC_Maximum_MPI(residual_max_norm);

      // Update CPU time used for the calculation so far.
      processor_cpu_time.update();
      total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput);

      // Periodically save restart solution files.
      if (!first_step &&
	  number_of_time_steps-Input_Parameters.Restart_Solution_Save_Frequency*
	  (number_of_time_steps/Input_Parameters.Restart_Solution_Save_Frequency) == 0) {
	if (!batch_flag)
	  cout << "\n\n  Saving LevelSet2D solution to restart data file(s) after"
	       << " n = " << number_of_time_steps << " steps (iterations).";
	// Write the quadtree restart file.
	error_flag = Write_QuadTree(QuadTree,Input_Parameters);
	if (error_flag) {
	  cout << "\n LevelSet2D ERROR: Unable to open LevelSet2D quadtree data file on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
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
	  cout << "\n LevelSet2D ERROR: Unable to open LevelSet2D restart output data file(s) "
	       << "on processor " << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;
	if (!batch_flag) cout << "\n";
      }

      // Output progress information for the calculation.
      if (!batch_flag) Output_Progress_L2norm(number_of_time_steps,
					      Time*THOUSAND,
					      total_cpu_time,
					      residual_l2_norm,
					      first_step,
					      Input_Parameters.Output_Progress_Frequency);
//  	if (!batch_flag) Output_Progress(number_of_time_steps,
//  					 Time*THOUSAND,
//  					 total_cpu_time,
//  					 residual_l1_norm,
//  					 first_step,
//  					 Input_Parameters.Output_Progress_Frequency);
      if (CFFC_Primary_MPI_Processor() && !first_step)
	Output_Progress_to_File(residual_file,
				number_of_time_steps,
				Time*THOUSAND,
				total_cpu_time,
				residual_l1_norm,
				residual_l2_norm,
				residual_max_norm);

      // Check to see if calculations are complete.
      if (!Input_Parameters.Time_Accurate &&
	  number_of_time_steps >= Input_Parameters.Maximum_Number_of_Time_Steps) break;
      if (Input_Parameters.Time_Accurate && Time >= Input_Parameters.Time_Max) break;

      // Update solution for next time step using a multistage time
      // stepping scheme.
      for (int i_stage = 1; i_stage <= Input_Parameters.N_Stage; i_stage++) {

	// Step 1. Exchange solution information between neighbouring blocks.
	error_flag = Send_All_Messages(Local_SolnBlk,
				       List_of_Local_Solution_Blocks,
				       NUM_VAR_LEVELSET2D,
				       OFF);
	if (error_flag) {
	  cout << "\n LevelSet2D ERROR: LevelSet2D message passing error on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;
	
	// Step 2. Apply boundary conditions for stage.
	BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

	// Step 3. Determine solution residuals for stage.
	error_flag = dUdt_Multistage_Hamilton_Jacobi(Local_SolnBlk,
						     List_of_Local_Solution_Blocks,
						     Input_Parameters,
						     i_stage);
	if (error_flag) {
	  cout << "\n LevelSet2D ERROR: LevelSet2D solution error on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;

	// Step 4. Update solution for stage.
	error_flag = Update_Solution_Multistage_Hamilton_Jacobi(Local_SolnBlk,
								List_of_Local_Solution_Blocks,
								Input_Parameters,
								i_stage);
	if (error_flag) {
	  cout << "\n LevelSet2D ERROR: LevelSet2D solution update error on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;

      }
      
      // Update time and time step counter.
      if (first_step) first_step = 0;
      number_of_time_steps++;
      if (Input_Parameters.i_Time_Integration != 
	  TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) {
	Time += Input_Parameters.Hamilton_Jacobi_CFL_Number*dTime;
      } else {
	Time += Input_Parameters.Hamilton_Jacobi_CFL_Number*dTime*
	        MultiStage_Optimally_Smoothing(Input_Parameters.N_Stage,
					       Input_Parameters.N_Stage,
					       Input_Parameters.i_Limiter);
      }

      // Exchange solution information between neighbouring blocks.
      error_flag = Send_All_Messages(Local_SolnBlk,
				     List_of_Local_Solution_Blocks,
				     NUM_VAR_LEVELSET2D,
				     OFF);
      if (error_flag) {
	cout << "\n LevelSet2D ERROR: LevelSet2D message passing error on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

      // Apply boundary conditions before redistancing.
      BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

      // Reinitialization procedure.
      if (Input_Parameters.i_Redistance_Criteria == EIKONAL_CRITERIA_THRESHOLD) {
	// Use a threshold-based reinitilization criteria.

	error_flag = Eikonal_Error(Local_SolnBlk,
				   Input_Parameters,
				   List_of_Local_Solution_Blocks,
				   global_error,
				   global_area);
	if (error_flag) {
	  cout << "\n LevelSet2D ERROR: LevelSet2D error during redistance error-checking on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;
	CFFC_Barrier_MPI();

	// Compile the errors from each processor.
	global_error = CFFC_Summation_MPI(global_error);
	global_area = CFFC_Summation_MPI(global_area);
	weighted_global_error = global_error/global_area;
	CFFC_Barrier_MPI();      

	// Perform Eikonal redistancing if error exceeds threshold.
	if (weighted_global_error > Input_Parameters.Eikonal_Threshold) {

	  // Reset.
	  global_error = ZERO;
	  global_area = ZERO;
	  weighted_global_error = ZERO;

	  // Notify the user that the function is being reinitialized.
	  if (!batch_flag) {
	    cout << "\n\n  Redistancing to a signed distance function after n = " << number_of_time_steps << " steps (iterations)." << endl;
	  }

	  // Redistance level set function to be a signed distance function.
	  error_flag = Explicit_Eikonal_Equation(Local_SolnBlk,
						 Input_Parameters,
						 QuadTree,
						 List_of_Global_Solution_Blocks,
						 List_of_Local_Solution_Blocks,
						 ON);
	  if (error_flag) {
	    cout << "\n LevelSet2D ERROR: LevelSet2D error during redistancing on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
	  }
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return error_flag;

	}

      } else {
	// Use a frequency-based reinitialization criteria.
	if (redistance_count == Input_Parameters.Redistance_Frequency-1) {
	  // Reset redistance count.
	  redistance_count = 0;
	  // Redistance level set function to be a signed distance function.
	  error_flag = Explicit_Eikonal_Equation(Local_SolnBlk,
						 Input_Parameters,
						 QuadTree,
						 List_of_Global_Solution_Blocks,
						 List_of_Local_Solution_Blocks,
						 ON);
	  if (error_flag) {
	    cout << "\n LevelSet2D ERROR: LevelSet2D error during redistancing on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
	  }
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return error_flag;
	} else {
	  // Increment redistance count.
	  redistance_count++;
	}
      }

    }

    if (!batch_flag) cout << "\n\n LevelSet2D computations complete on " 
			  << Date_And_Time() << ".\n";

  }

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Update ghostcell information and prescribe boundary conditions to 
  // ensure that the solution is consistent on each block.
  error_flag = Send_All_Messages(Local_SolnBlk,
				 List_of_Local_Solution_Blocks,
				 NUM_VAR_LEVELSET2D,
				 OFF);
  if (error_flag) {
    cout << "\n LevelSet2D ERROR: LevelSet2D message passing error on processor "
	 << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);
  
  // Close residual file.
  error_flag = Close_Progress_File(residual_file);


  /********************************************************************
   * Solution calculations complete.                                  *
   * Write 2D LevelSet solution to output and restart files as        *
   * required, reset solution parameters, and run other cases as      *
   * specified by input parameters.                                   *
   ********************************************************************/
  
 postprocess_current_calculation: ;

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  while (1) {
    if (CFFC_Primary_MPI_Processor()) {
      Get_Next_Input_Control_Parameter(Input_Parameters);
      command_flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
      line_number = Input_Parameters.Line_Number;
    }
    CFFC_Barrier_MPI();
    Broadcast_Input_Parameters(Input_Parameters);
    CFFC_Broadcast_MPI(&command_flag,1);
    
    if (command_flag == EXECUTE_CODE) {
      // Deallocate memory for 2D LevelSet equation solution.
      if (!batch_flag) cout << "\n Deallocating LevelSet2D solution variables.";
      Local_SolnBlk = Deallocate(Local_SolnBlk,Input_Parameters);
      //Deallocate_Message_Buffers(List_of_Local_Solution_Blocks); // Not necessary here!
      List_of_Local_Solution_Blocks.deallocate();
      List_of_Global_Solution_Blocks.deallocate();
      QuadTree.deallocate();
      MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk,
					    Input_Parameters.Number_of_Blocks_Idir,
					    Input_Parameters.Number_of_Blocks_Jdir);
      // Output input parameters for new caluculation.
      if (!batch_flag)
	cout << "\n\n Starting a new calculation." << Input_Parameters << "\n";
      // Execute new calculation.
      goto execute_new_calculation;
      
    } else if (command_flag == TERMINATE_CODE) {
      // Deallocate memory for 2D LevelSet equation solution.
      if (!batch_flag) cout << "\n Deallocating LevelSet2D solution variables.";
      Local_SolnBlk = Deallocate(Local_SolnBlk,Input_Parameters);
      //Deallocate_Message_Buffers(List_of_Local_Solution_Blocks); // Not necessary here!
      List_of_Local_Solution_Blocks.deallocate();
      List_of_Global_Solution_Blocks.deallocate();
      QuadTree.deallocate();
      MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk,
					    Input_Parameters.Number_of_Blocks_Idir,
					    Input_Parameters.Number_of_Blocks_Jdir);
      // Close input data file.
      if (!batch_flag) cout << "\n\n Closing LevelSet2D input data file.";
      if (CFFC_Primary_MPI_Processor()) Close_Input_File(Input_Parameters);
      // Terminate calculation.
      return 0;
      
    } else if (command_flag == CONTINUE_CODE) {
      // Reset maximum time step counter.
      Input_Parameters.Maximum_Number_of_Time_Steps += number_of_time_steps;
      // Output input parameters for continuing calculation.
      if (!batch_flag)
	cout << "\n\n Continuing existing calculation."
	     << Input_Parameters << "\n";
      // Continue existing calculation.
      goto continue_existing_calculation;

    } else if (command_flag == REFINE_GRID_CODE) {
      // Refine mesh.  Must first retrieve the interface and broadcast 
      // it to all of the processors.  The extension problems must be 
      // resolved on the new mesh.

      // Refine mesh using block based adaptive mesh refinement algorithm.
      if (!batch_flag) cout << "\n\n Refining Grid.  Performing adaptive mesh refinement.";

      // MPI barrier to ensure processor synchronization.
      CFFC_Barrier_MPI();

      // Update ghostcell information and prescribe boundary conditions 
      // to ensure that the solution is consistent on each block.
      error_flag = Send_All_Messages(Local_SolnBlk,
				     List_of_Local_Solution_Blocks,
				     NUM_VAR_LEVELSET2D,
				     OFF);
      if (error_flag) {
	cout << "\n LevelSet2D ERROR: LevelSet2D message passing error on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

      BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

      // Refine mesh.
      error_flag = AMR(Local_SolnBlk,
		       Input_Parameters,
                       QuadTree,
                       List_of_Global_Solution_Blocks,
                       List_of_Local_Solution_Blocks,
                       ON,ON);
      if (error_flag) {
         cout << "\n LevelSet2D ERROR: LevelSet2D AMR error on processor "
              << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;
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
	     << QuadTree.efficiencyRefinement() << "\n";
      }

      // Explicit solution of the Eikonal equation for iteratively forcing the
      // level set function, psi, to be a signed distance fuction.
      if (!batch_flag) cout << "\n Performing the solution of the Eikonal equation after refinement.";
      error_flag = Explicit_Eikonal_Equation(Local_SolnBlk,
					     Input_Parameters,
					     QuadTree,
					     List_of_Global_Solution_Blocks,
					     List_of_Local_Solution_Blocks,
					     ON);
      if (error_flag) {
	cout << "\n LevelSet2D ERROR: LevelSet2D error during"
	     << " solution of the Eikonal equation on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

      // MPI barrier to ensure processor synchronization.
      CFFC_Barrier_MPI();

      // Update ghostcell information and prescribe boundary conditions 
      // to ensure that the solution is consistent on each block.
      error_flag = Send_All_Messages(Local_SolnBlk,
				     List_of_Local_Solution_Blocks,
				     NUM_VAR_LEVELSET2D,
				     OFF);
      if (error_flag) {
	cout << "\n LevelSet2D ERROR: LevelSet2D message passing error on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
	cout.flush();
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;
      
      BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);


    } else if (command_flag == WRITE_OUTPUT_CODE) {
      // Output solution data.
      if (!batch_flag) cout << "\n Writing LevelSet2D solution to output data file(s).";
      error_flag = Output_Tecplot(Local_SolnBlk,
				  List_of_Local_Solution_Blocks,
				  Input_Parameters,
				  number_of_time_steps,
				  Time);
      if (error_flag) {
	cout << "\n LevelSet2D ERROR: Unable to open LevelSet2D output data file(s) "
	     << "on processor " << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
	cout.flush();
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_CELLS_CODE) {
      // Output solution data.
      if (!batch_flag) cout << "\n Writing cell-centered LevelSet2D solution to output data file(s).";
      error_flag = Output_Cells_Tecplot(Local_SolnBlk,
					List_of_Local_Solution_Blocks,
					Input_Parameters,
					number_of_time_steps,
					Time);
      if (error_flag) {
	cout << "\n LevelSet2D ERROR: Unable to open LevelSet2D cell output data file(s) "
	     << "on processor " << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
	cout.flush();
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_NODES_CODE) {
      // Output solution data.
      if (!batch_flag) cout << "\n Writing LevelSet2D node locations to output data file(s).";
      error_flag = Output_Nodes_Tecplot(Local_SolnBlk,
					List_of_Local_Solution_Blocks,
					Input_Parameters);
      if (error_flag) {
	cout << "\n LevelSet2D ERROR: Unable to open LevelSet2D node output data file(s) "
	     << "on processor " << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
	cout.flush();
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_INTERFACE_COMPONENT_LIST_CODE) {

      if (!batch_flag) {
	cout << "\n Writing LevelSet2D interface nodes to output data file(s).";
	cout.flush();
      }

      // MPI barrier to ensure processor synchronization.
      CFFC_Barrier_MPI();

      // Retrieve spline(s).
      error_flag = Retrieve_Interface_Spline(Local_SolnBlk,
					     List_of_Local_Solution_Blocks);
      if (error_flag) {
	cout << "\n LevelSet2D ERROR: Error #" << error_flag
	     << " occured when recapturing front location on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
	cout.flush();
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

      // Broadcast spline(s).
      error_flag = Share_Interface_Information(Local_SolnBlk,
					       QuadTree,
					       List_of_Global_Solution_Blocks,
					       List_of_Local_Solution_Blocks,
					       Input_Parameters);
      if (error_flag) {
	cout << "\n LevelSet2D ERROR: LevelSet2D sharing interface error on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
	cout.flush();
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

      // Output solution data.
      //if (CFFC_Primary_MPI_Processor())
      error_flag = Output_Interface_Tecplot(Local_SolnBlk,
					    List_of_Local_Solution_Blocks,
					    Input_Parameters,
					    number_of_time_steps,
					    Time);
      if (error_flag) {
	cout << "\n LevelSet2D ERROR: Unable to open LevelSet2D output "
	     << "data file(s) on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
	cout.flush();
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_LEVEL_SET_CIRCLE_CODE) {
      if (!batch_flag) cout << "\n Writing comparison of the LevelSet2D solution with the exact solution"
			    << "\n for a circle to the output data file(s).";
      error_flag = Output_Circle_Tecplot(Local_SolnBlk,
					 List_of_Local_Solution_Blocks,
					 Input_Parameters,
					 number_of_time_steps,
					 Time);
      if (error_flag) {
	cout << "\n LevelSet2D ERROR: Unable to open LevelSet2D circle output data file(s) "
	     << "on processor " << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
	cout.flush();
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_LEVEL_SET_ELLIPSE_CODE) {
      if (!batch_flag) cout << "\n Writing comparison of the LevelSet2D solution with the exact solution"
			    << "\n for an ellipse to the output data file(s).";
      error_flag = Output_Ellipse_Tecplot(Local_SolnBlk,
					  List_of_Local_Solution_Blocks,
					  Input_Parameters,
					  number_of_time_steps,
					  Time);
      if (error_flag) {
	cout << "\n LevelSet2D ERROR: Unable to open LevelSet2D ellipse output data file(s) "
	     << "on processor " << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
	cout.flush();
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_LEVEL_SET_ZALESAK_CODE) {
      if (!batch_flag) cout << "\n Writing comparison of the LevelSet2D solution with the exact solution"
			    << "\n for a Zalesak's Disk to the output data file(s).";
      error_flag = Output_Zalesaks_Disk_Tecplot(Local_SolnBlk,
						List_of_Local_Solution_Blocks,
						Input_Parameters,
						number_of_time_steps,
						Time);
      if (error_flag) {
	cout << "\n LevelSet2D ERROR: Unable to open LevelSet2D Zalesak's Disk output data file(s) "
	     << "on processor " << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
	cout.flush();
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_RESTART_CODE) {

      // Write restart files.
      if (!batch_flag) cout << "\n Writing LevelSet2D solution to restart data file(s).";

      // MPI barrier to ensure processor synchronization.
      CFFC_Barrier_MPI();

      // Write the quadtree restart file.
      error_flag = Write_QuadTree(QuadTree,Input_Parameters);
      if (error_flag) {
	cout << "\n LevelSet2D ERROR: Unable to open LevelSet2D quadtree data file "
	     << "on processor " << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
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
	cout << "\n LevelSet2D ERROR: Unable to open LevelSet2D restart output data file(s) "
	     << "on processor " << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
	cout.flush();
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

      // Retrieve spline(s).
      error_flag = Retrieve_Interface_Spline(Local_SolnBlk,List_of_Local_Solution_Blocks);
      if (error_flag) {
	cout << "\n LevelSet2D ERROR: An error occured when recapturing"
	     << " front location on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
	cout.flush();
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

      // Broadcast spline(s).
      error_flag = Share_Interface_Information(Local_SolnBlk,
					       QuadTree,
					       List_of_Global_Solution_Blocks,
					       List_of_Local_Solution_Blocks,
					       Input_Parameters);
      if (error_flag) {
	cout << "\n LevelSet2D ERROR: LevelSet2D sharing interface error on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
	cout.flush();
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_GRID_CODE) {
      // Output multi-block solution-adaptive mesh data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag)cout << "\n Writing LevelSet2D multi-block mesh to grid data output file.";
	error_flag = Output_Tecplot(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << "\n LevelSet2D ERROR: Unable to open LevelSet2D mesh data output file.\n";
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_GRID_DEFINITION_CODE) {
      // Write multi-block solution-adaptive mesh definition files.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing LevelSet2D multi-block mesh to grid definition files.";
	error_flag = Write_Multi_Block_Grid_Definition(MeshBlk,Input_Parameters);
	error_flag = Write_Multi_Block_Grid(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << "\n LevelSet2D ERROR: Unable to open LevelSet2D multi-block mesh definition files.\n";
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_GRID_NODES_CODE) {
      // Output multi-block solution-adaptive mesh node data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing LevelSet2D multi-block mesh to node data output file.";
	error_flag = Output_Nodes_Tecplot(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << "\n LevelSet2D ERROR: Unable to open LevelSet2D mesh node data output file.\n";
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_GRID_CELLS_CODE) {
      // Output multi-block solution-adaptive mesh cell data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing LevelSet2D multi-block mesh to cell data output file.";
	error_flag = Output_Cells_Tecplot(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << "\n LevelSet2D ERROR: Unable to open LevelSet2D mesh cell data output file.\n";
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      
    } else if (command_flag == INVALID_INPUT_CODE ||
	       command_flag == INVALID_INPUT_VALUE) {
      line_number = -line_number;
      cout << "\n LevelSet2D ERROR: Error reading LevelSet2D data at line #"
	   << -line_number  << " of input data file.\n";
      cout.flush();
      return line_number;
    }
    
  }

  /********************************************************************  
   * End of all LevelSet2DSolver computations and I/O.                *
   ********************************************************************/
  return 0;

}

/**********************************************************************
 * Routine: Initialize_Level_Set_Solution                             *
 *                                                                    *
 *                                                                    *
 *                                                                    *
 **********************************************************************/
LevelSet2D_Quad_Block* Initialize_Level_Set_Solution(char *Input_File_Name,
						     const int &batch_flag,
						     int &error_flag,
						     LevelSet2D_Quad_Block *Local_SolnBlk,
						     LevelSet2D_Input_Parameters &Input_Parameters,
						     QuadTreeBlock_DataStructure &QuadTree,
						     AdaptiveBlockResourceList &List_of_Global_Solution_Blocks,
						     AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
						     const Interface2D_List &Interface_List) {

  int command_flag;

  // Multi-block solution-adaptive quadrilateral mesh solution variables.
  Grid2D_Quad_Block **MeshBlk;

  // LevelSet2D input file name.
  char Input_File_Name_ptr[256];

  // Determine LevelSet2D input file name.
  strcpy(Input_File_Name_ptr,"levelset_");
  strcat(Input_File_Name_ptr,Input_File_Name);

  // The primary MPI processor processes the input parameter file.
  if (CFFC_Primary_MPI_Processor()) {
    if (!batch_flag) cout << "\n\n Reading LevelSet2D input data file `"
			  << Input_File_Name_ptr << "'.";
    // Process input file.
    error_flag = Process_Input_Control_Parameter_File(Input_Parameters,
						      Input_File_Name_ptr,
						      command_flag);
    // Output input data to screen.
    if (!batch_flag && error_flag == 0) cout << Input_Parameters << "\n";
  } else {
    error_flag = 0;
  }

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Broadcast input solution parameters to other MPI processors.
  CFFC_Broadcast_MPI(&error_flag,1);
  if (error_flag != 0) return NULL;
  CFFC_Broadcast_MPI(&command_flag,1);
  if (command_flag == TERMINATE_CODE) return NULL;
  Broadcast_Input_Parameters(Input_Parameters);

  /********************************************************************
   * Create initial mesh and allocate LevelSet2D solution variables   *
   * specified IBVP/BVP problem.                                      *
   ********************************************************************/
  
  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Create initial mesh.  Read mesh from grid definition or data files 
  // when specified by input parameters.

  // The primary MPI processor creates the initial mesh.
  if (CFFC_Primary_MPI_Processor()) {

    if (!batch_flag)
      cout << "\n Creating (or reading) initial levelset2D quadrilateral multi-block mesh.";

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
      cout << "\n LevelSet2D ERROR: Unable to create valid LevelSet2D multi-block mesh." << "\n";
      cout.flush();
    }
  } else {
    MeshBlk = NULL;
  }

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Broadcast the mesh to other MPI processors.
  CFFC_Broadcast_MPI(&error_flag,1);
  if (error_flag) return NULL;
  MeshBlk = Broadcast_Multi_Block_Grid(MeshBlk,Input_Parameters);

  // Create (allocate) multi-block quadtree data structure, create
  // (allocate) array of local 2D LevelSet equation solution blocks,
  // assign and create (allocate) 2D LevelSet equation solution blocks
  // corresponding to the initial mesh.
  if (!batch_flag) cout << "\n Creating multi-block quadtree data structure and assigning"
                        << "\n  LevelSet2D solution blocks corresponding to initial mesh.";
  Local_SolnBlk = Create_Initial_Solution_Blocks(MeshBlk,
						 Local_SolnBlk,
						 Input_Parameters,
						 QuadTree,
						 List_of_Global_Solution_Blocks,
						 List_of_Local_Solution_Blocks);
  if (Local_SolnBlk == NULL) return NULL;

  /********************************************************************
   * Initialize LevelSet2D solution variables.                        *
   ********************************************************************/

  // Initialize the conserved and primitive state solution variables.
  if (!batch_flag) cout << "\n Prescribing LevelSet2D initial data.";

  // Initialize interface(s).
  error_flag = Initialize_Interfaces(Local_SolnBlk,
				     List_of_Local_Solution_Blocks,
				     Input_Parameters,
				     Interface_List);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return NULL;

  // Construct the bulk flow-field.
  error_flag = Construct_Bulk_Flow_Field(Local_SolnBlk,
					 List_of_Local_Solution_Blocks,
					 Input_Parameters);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return NULL;

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Send solution information between neighbouring blocks to complete
  // prescription of initial data.
  error_flag = Send_All_Messages(Local_SolnBlk,
                                 List_of_Local_Solution_Blocks,
                                 NUM_COMP_VECTOR2D,
                                 ON);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return NULL;

  // Solve geometric extension problem.
  if (!batch_flag) cout << "\n Solving initial geometric extension problem.";
  error_flag = Geometric_Extension_Problem(Local_SolnBlk,
					   List_of_Local_Solution_Blocks,
					   Input_Parameters);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return NULL;

  // Redistance level set function to be a signed distance function.
  if (!batch_flag) cout << "\n Performing initial solution of the Eikonal equation.";
  error_flag = Explicit_Eikonal_Equation(Local_SolnBlk,
					 Input_Parameters,
					 QuadTree,
					 List_of_Global_Solution_Blocks,
					 List_of_Local_Solution_Blocks,
					 OFF);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return NULL;

  // Solve scalar (front speed) extension equation.
  if (!batch_flag) cout << "\n Performing solution of the scalar (front speed) extension equation.";
  error_flag = Explicit_Scalar_Extension_Equation(Local_SolnBlk,
						  Input_Parameters,
						  QuadTree,
						  List_of_Global_Solution_Blocks,
						  List_of_Local_Solution_Blocks);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return NULL;

  // Update ghostcell information and prescribe boundary conditions 
  // to ensure that the solution is consistent on each block.
  error_flag = Send_All_Messages(Local_SolnBlk,
				 List_of_Local_Solution_Blocks,
				 NUM_VAR_LEVELSET2D,
				 OFF);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return NULL;

  // Prescribe boundary data consistent with initial data.
  BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

  // Perform uniform mesh refinement.
  if (!batch_flag) cout << "\n Performing LevelSet2D interface mesh refinement.";
  error_flag = Uniform_Adaptive_Mesh_Refinement(Local_SolnBlk,
						Input_Parameters,
						QuadTree,
						List_of_Global_Solution_Blocks,
						List_of_Local_Solution_Blocks,
						Interface_List);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return NULL;

  // Perform interface mesh refinement.
  if (!batch_flag) cout << "\n Performing LevelSet2D initial mesh refinement.";
  error_flag = Initial_Adaptive_Mesh_Refinement(Local_SolnBlk,
						Input_Parameters,
						QuadTree,
						List_of_Global_Solution_Blocks,
						List_of_Local_Solution_Blocks,
 						Interface_List);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return NULL;

  // Prescribe boundary data consistent with initial data.
  BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

  // Capture zero-level set contour(s) (interface spline(s)).
  if (!batch_flag) cout << "\n Retrieving interface spline(s).";
  error_flag = Retrieve_Interface_Spline(Local_SolnBlk,
					 List_of_Local_Solution_Blocks);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return NULL;

  // Synchronize processors.
  CFFC_Barrier_MPI();

  // Broadcast spline(s).
  if (!batch_flag) cout << "\n Sharing interface information across processors.";
  error_flag = Share_Interface_Information(Local_SolnBlk,
					   QuadTree,
					   List_of_Global_Solution_Blocks,
					   List_of_Local_Solution_Blocks,
					   Input_Parameters);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return NULL;

  // Output multi-block solution-adaptive quadrilateral mesh statistics.
  if (!batch_flag) {
    cout << "\n\n LevelSet2D multi-block solution-adaptive quadrilateral mesh statistics: "; 
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
	 << QuadTree.efficiencyRefinement() << "\n"; 
  }

  // Deallocate mesh block.
  MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk,
					Input_Parameters.Number_of_Blocks_Idir, 
					Input_Parameters.Number_of_Blocks_Jdir);

  // Level set problem successfully initialized.
  return Local_SolnBlk;

}

/**********************************************************************
 * Routine: Restart_Level_Set_Solution                                *
 *                                                                    *
 *                                                                    *
 *                                                                    *
 **********************************************************************/
LevelSet2D_Quad_Block* Restart_Level_Set_Solution(char *Input_File_Name,
						  const int &batch_flag,
						  int &error_flag,
						  int &number_of_time_steps,
						  double &Time,
						  LevelSet2D_Quad_Block *Local_SolnBlk,
						  LevelSet2D_Input_Parameters &Input_Parameters,
						  QuadTreeBlock_DataStructure &QuadTree,
						  AdaptiveBlockResourceList &List_of_Global_Solution_Blocks,
						  AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
						  const Interface2D_List &Interface_List) {

  CPUTime processor_cpu_time, total_cpu_time;
  int command_flag;

  // Multi-block solution-adaptive quadrilateral mesh solution variables.
  Grid2D_Quad_Block **MeshBlk;

  // LevelSet2D input file name.
  char Input_File_Name_ptr[256];

  // Determine LevelSet2D input file name.
  strcpy(Input_File_Name_ptr,"levelset_");
  strcat(Input_File_Name_ptr,Input_File_Name);

  // The primary MPI processor processes the input parameter file.
  if (CFFC_Primary_MPI_Processor()) {
    if (!batch_flag) cout << "\n\n Reading LevelSet2D input data file `"
			  << Input_File_Name_ptr << "'.";
    // Process input file.
    error_flag = Process_Input_Control_Parameter_File(Input_Parameters,
						      Input_File_Name_ptr,
						      command_flag);
    // Output input data to screen.
    if (!batch_flag && error_flag == 0) cout << Input_Parameters << "\n";
  } else {
    error_flag = 0;
  }

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Broadcast input solution parameters to other MPI processors.
  CFFC_Broadcast_MPI(&error_flag,1);
  if (error_flag != 0) return NULL;
  CFFC_Broadcast_MPI(&command_flag,1);
  if (command_flag == TERMINATE_CODE) return NULL;
  Broadcast_Input_Parameters(Input_Parameters);

  /********************************************************************
   * Create initial mesh and allocate LevelSet2D solution variables   *
   * specified IBVP/BVP problem.                                      *
   ********************************************************************/
  
  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Create initial mesh.  Read mesh from grid definition or data files 
  // when specified by input parameters.

  // The primary MPI processor creates the initial mesh.
  if (CFFC_Primary_MPI_Processor()) {

    if (!batch_flag)
      cout << "\n Creating (or reading) initial levelset2D quadrilateral multi-block mesh.";

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
      cout << "\n LevelSet2D ERROR: Unable to create valid LevelSet2D multi-block mesh." << "\n";
      cout.flush();
    }
  } else {
    MeshBlk = NULL;
  }

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Broadcast the mesh to other MPI processors.
  CFFC_Broadcast_MPI(&error_flag,1);
  if (error_flag) return NULL;
  MeshBlk = Broadcast_Multi_Block_Grid(MeshBlk,Input_Parameters);

  // Create (allocate) multi-block quadtree data structure, create
  // (allocate) array of local 2D LevelSet equation solution blocks,
  // assign and create (allocate) 2D LevelSet equation solution blocks
  // corresponding to the initial mesh.
  if (!batch_flag) cout << "\n Creating multi-block quadtree data structure and assigning"
                        << "\n  LevelSet2D solution blocks corresponding to initial mesh.";
  Local_SolnBlk = Create_Initial_Solution_Blocks(MeshBlk,
						 Local_SolnBlk,
						 Input_Parameters,
						 QuadTree,
						 List_of_Global_Solution_Blocks,
						 List_of_Local_Solution_Blocks);
  if (Local_SolnBlk == NULL) return NULL;

  /********************************************************************
   * Initialize LevelSet2D solution variables.                        *
   ********************************************************************/

  // Set the initial time level.
  Time = ZERO;
  number_of_time_steps = 0;

  // Set the CPU time to zero.
  processor_cpu_time.zero();
  total_cpu_time.zero();

  // Initialize interface(s).
  if (!batch_flag) cout << "\n Creating LevelSet2D interface list.";
  error_flag = Initialize_Interfaces(Local_SolnBlk,
				     List_of_Local_Solution_Blocks,
				     Input_Parameters,
				     Interface_List);
  if (error_flag) return NULL;

  // Initialize the conserved and primitive state solution variables.
  if (!batch_flag) cout << "\n Reading LevelSet2D solution from restart data files.";
  // Read the quadtree restart file.
  error_flag = Read_QuadTree(QuadTree,
			     List_of_Global_Solution_Blocks,
			     List_of_Local_Solution_Blocks,
			     Input_Parameters);
  if (error_flag) {
    cout << "\n LevelSet2D ERROR: Unable to open LevelSet2D quadtree data file "
	 << "on processor "
	 << List_of_Local_Solution_Blocks.ThisCPU
	 << "." << "\n";
    cout.flush();
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag)  return NULL;
  // Allocate the message buffers.
  Allocate_Message_Buffers(List_of_Local_Solution_Blocks,
			   Local_SolnBlk[0].NumVar()+NUM_COMP_VECTOR2D);
  // Read the solution block restart files.
  error_flag = Read_Restart_Solution(Local_SolnBlk,
				     List_of_Local_Solution_Blocks,
				     Input_Parameters,
				     number_of_time_steps,
				     Time,
				     processor_cpu_time);
  if (error_flag) {
    cout << "\n LevelSet2D ERROR: Unable to open LevelSet2D restart input data file(s) "
	 << "on processor "
	 << List_of_Local_Solution_Blocks.ThisCPU
	 << "." << "\n";
    cout.flush();
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return NULL;
  // Ensure each processor has the correct time and time.
  number_of_time_steps = CFFC_Maximum_MPI(number_of_time_steps);
  Time = CFFC_Maximum_MPI(Time);
  processor_cpu_time.cput = CFFC_Maximum_MPI(processor_cpu_time.cput);
  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();
  // Broadcast input solution parameters to other MPI processors.
  CFFC_Broadcast_MPI(&error_flag,1);
  if (error_flag != 0) return NULL;
  CFFC_Broadcast_MPI(&command_flag,1);
  if (command_flag == TERMINATE_CODE) return NULL;
  Broadcast_Input_Parameters(Input_Parameters);

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Send solution information between neighbouring blocks to complete
  // prescription of initial data.
  error_flag = Send_All_Messages(Local_SolnBlk,
                                 List_of_Local_Solution_Blocks,
                                 NUM_COMP_VECTOR2D,
                                 ON);
  if (error_flag) {
    cout << "\n LevelSet2D ERROR: Message passing error during LevelSet2D solution intialization "
	 << "on processor " << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
    cout.flush();
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return NULL;

  // Output multi-block solution-adaptive quadrilateral mesh statistics.
  if (!batch_flag) {
    cout << "\n\n LevelSet2D multi-block solution-adaptive quadrilateral mesh statistics: "; 
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
	 << QuadTree.efficiencyRefinement() << "\n"; 
  }

  // Deallocate mesh block.
  MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
					Input_Parameters.Number_of_Blocks_Idir, 
					Input_Parameters.Number_of_Blocks_Jdir);

  // Level set problem successfully initialized.
  return Local_SolnBlk;

}

/**********************************************************************
 * Routine: Evolve_Level_Set_Solution                                 *
 **********************************************************************/
int Evolve_Level_Set_Solution(const int &batch_flag,
			      LevelSet2D_Quad_Block *Local_SolnBlk,
			      LevelSet2D_Input_Parameters &Input_Parameters,
			      QuadTreeBlock_DataStructure &QuadTree,
			      AdaptiveBlockResourceList &List_of_Global_Solution_Blocks,
			      AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
			      double &Time,
			      const double &Time_Max,
			      int &number_of_time_steps) {

  int error_flag, command_flag;
  int redistance_count = 0, first_step = 1;
  double dTime;

  double global_error, global_area, weighted_global_error;
  global_error = ZERO;
  global_area = ZERO;
  weighted_global_error = ZERO;

//   char extension[256], output_file_name[256];
//   char *output_file_name_ptr;
//   ofstream dout;
//   strcpy(output_file_name,"evolve_cpu");
//   sprintf(extension,"%.6d",List_of_Local_Solution_Blocks.ThisCPU);
//   strcat(extension,".txt");
//   strcat(output_file_name,extension);
//   output_file_name_ptr = output_file_name;
//   dout.open(output_file_name_ptr,ios::out);
//   if (dout.bad()) return ;

  // Perform required number of iterations (time steps).
  if (Time_Max > Time) {

    while (1) {

      // Periodically refine the mesh (AMR).
      if (Input_Parameters.AMR) {
	if (!first_step &&
	    number_of_time_steps-Input_Parameters.AMR_Frequency*
	    (number_of_time_steps/Input_Parameters.AMR_Frequency) == 0) {
	  if (!batch_flag)
	    cout << "\n\n Refining Level Set Grid.  Performing adaptive mesh refinement at n = "
		 << number_of_time_steps << ".";
	  error_flag = AMR(Local_SolnBlk,
			   Input_Parameters,
			   QuadTree,
			   List_of_Global_Solution_Blocks,
			   List_of_Local_Solution_Blocks,
			   ON,ON);
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
	  }
	}
      }

      // Determine local and global time steps.
      dTime = CFL_Hamilton_Jacobi(Local_SolnBlk,
				  List_of_Local_Solution_Blocks,
				  Input_Parameters);
      // Find global minimum time step for all processors.
      dTime = CFFC_Minimum_MPI(dTime);
      if (Time + Input_Parameters.Hamilton_Jacobi_CFL_Number*dTime > Time_Max)
	dTime = (Time_Max-Time)/Input_Parameters.Hamilton_Jacobi_CFL_Number;
      // Set global time step.
      Set_Global_TimeStep(Local_SolnBlk,List_of_Local_Solution_Blocks,dTime);

      // Check to see if calculations are complete.
      if (Time >= Time_Max) break;

      if (!batch_flag) cout << "o";

      // Update solution for next time step using a multistage time
      // stepping scheme.
      for (int i_stage = 1; i_stage <= Input_Parameters.N_Stage; i_stage++) {

	// Step 1. Exchange solution information between neighbouring blocks.
	error_flag = Send_All_Messages(Local_SolnBlk,
				       List_of_Local_Solution_Blocks,
				       NUM_VAR_LEVELSET2D,
				       OFF);
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;
	
	// Step 2. Apply boundary conditions for stage.
	BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

	// Step 3. Determine solution residuals for stage.
	error_flag = dUdt_Multistage_Hamilton_Jacobi(Local_SolnBlk,
						     List_of_Local_Solution_Blocks,
						     Input_Parameters,
						     i_stage);
	if (error_flag) return error_flag;

	// Step 4. Update solution for stage.
	error_flag = Update_Solution_Multistage_Hamilton_Jacobi(Local_SolnBlk,
								List_of_Local_Solution_Blocks,
								Input_Parameters,
								i_stage);
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;

      }
      
      // Update time and time step counter.
      if (first_step) first_step = 0;
      number_of_time_steps++;
      Time += Input_Parameters.Hamilton_Jacobi_CFL_Number*dTime;

      // Exchange solution information between neighbouring blocks.
      error_flag = Send_All_Messages(Local_SolnBlk,
				     List_of_Local_Solution_Blocks,
				     NUM_VAR_LEVELSET2D,
				     OFF);
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

      // Apply boundary conditions before redistancing.
      BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

      // Samuel: Calculate the error in the Eikonal equation solution.
      error_flag = Eikonal_Error(Local_SolnBlk,
				 Input_Parameters,
				 List_of_Local_Solution_Blocks,
				 global_error,
				 global_area);
      if (error_flag) {
	cout << "\n LevelSet2D ERROR: LevelSet2D error during redistance error-checking on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

      // MPI barrier to ensure processor synchronization.
      CFFC_Barrier_MPI();

      global_error = CFFC_Summation_MPI(global_error);
      global_area = CFFC_Summation_MPI(global_area);
      weighted_global_error = global_error/global_area;

      // MPI barrier to ensure processor synchronization.
      CFFC_Barrier_MPI();

      // Print out the weighted global error.
      if (!batch_flag) {
	cout << "\n Weighted global error = " << weighted_global_error << endl;
      }

      // Perform Eikonal redistancing if necessary.
      //  if (weighted_global_error > Input_Parameters.Eikonal_Threshold) {
      if (redistance_count == Input_Parameters.Redistance_Frequency-1) {
	
	// Reset.
	redistance_count = 0;
	global_error = ZERO;
	global_area = ZERO;
	weighted_global_error = ZERO;

	if (!batch_flag) {
	  cout << "\n\n  Redistancing to a signed distance function after n = " << number_of_time_steps << " steps (iterations).\n";
	}
	error_flag = Explicit_Eikonal_Equation(Local_SolnBlk,
					       Input_Parameters,
					       QuadTree,
					       List_of_Global_Solution_Blocks,
					       List_of_Local_Solution_Blocks,
					       ON);
	if (error_flag) {
	  cout << "\n LevelSet2D ERROR: LevelSet2D error during redistancing on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;
      } else {
	// Increment redistance count.
	redistance_count++;
      }

//       // Redistance level set function to be a signed distance function.
//       if (redistance_count == Input_Parameters.Redistance_Frequency-1) {
// 	redistance_count = 0;
// 	error_flag = Explicit_Eikonal_Equation(Local_SolnBlk,
// 					       Input_Parameters,
// 					       QuadTree,
// 					       List_of_Global_Solution_Blocks,
// 					       List_of_Local_Solution_Blocks,
// 					       ON);
// 	error_flag = CFFC_OR_MPI(error_flag);
// 	if (error_flag) return error_flag;
//       } else {
// 	// Increment redistance count.
// 	redistance_count++;
//       }

    }

  }

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Update ghostcell information and prescribe boundary conditions to 
  // ensure that the solution is consistent on each block.
  error_flag = Send_All_Messages(Local_SolnBlk,
				 List_of_Local_Solution_Blocks,
				 NUM_VAR_LEVELSET2D,
				 OFF);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Apply boundary conditions.
  BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

  // Capture zero-level set contour(s) (interface spline(s)).
  error_flag = Retrieve_Interface_Spline(Local_SolnBlk,
					 List_of_Local_Solution_Blocks);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) {
    error_flag = Output_Tecplot(Local_SolnBlk,
				List_of_Local_Solution_Blocks,
				Input_Parameters,
				0,ZERO);
    error_flag = Output_Cells_Tecplot(Local_SolnBlk,
				      List_of_Local_Solution_Blocks,
				      Input_Parameters,
				      0,ZERO);
    return 88888;
  }

  // Broadcast spline(s).
  error_flag = Share_Interface_Information(Local_SolnBlk,
					   QuadTree,
					   List_of_Global_Solution_Blocks,
					   List_of_Local_Solution_Blocks,
					   Input_Parameters);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) {
    error_flag = Output_Tecplot(Local_SolnBlk,
				List_of_Local_Solution_Blocks,
				Input_Parameters,
				0,ZERO);
    error_flag = Output_Cells_Tecplot(Local_SolnBlk,
				      List_of_Local_Solution_Blocks,
				      Input_Parameters,
				      0,ZERO);
    error_flag = Output_Interface_Tecplot(Local_SolnBlk,
					  List_of_Local_Solution_Blocks,
					  Input_Parameters,
					  0,ZERO);
    return 88888;
  }
  //   if (error_flag) return error_flag;

  //dout.close();

  // Level set evolution successfully computed.
  return 0;

}
