/**********************************************************************
 * Electrostatic2DQuadSolvers.cc: 2D Electrostatic equation multi-    *
 *                                block quadrilateral mesh solvers.   *
 **********************************************************************/

// Include 2D Electrostatic quadrilateral mesh solution header file.

#ifndef _ELECTROSTATIC2D_QUAD_INCLUDED
#include "Electrostatic2DQuad.h"
#endif // _ELECTROSTATIC2D_QUAD_INCLUDED

/**********************************************************************
 * Routine: Electrostatic2DQuadSolver                                 *
 *                                                                    *
 * Computes solutions to 2D Electrostatic equations on 2D             *
 * quadrilateral multi-block solution-adaptive mesh.                  *
 *                                                                    *
 **********************************************************************/
int Electrostatic2DQuadSolver(char *Input_File_Name_ptr, int batch_flag) {

  /********************************************************************
   * Local variable declarations.                                     *
   ********************************************************************/

  // Electrostatic2D input variables and parameters.
  Electrostatic2D_Input_Parameters      Input_Parameters;

  // Multi-block solution-adaptive quadrilateral mesh solution variables.
  Grid2D_Quad_Block           **MeshBlk;
  QuadTreeBlock_DataStructure   QuadTree;
  AdaptiveBlockResourceList     List_of_Global_Solution_Blocks;
  AdaptiveBlock2D_List          List_of_Local_Solution_Blocks;
  Electrostatic2D_Quad_Block   *Local_SolnBlk;

  // Define residual file and cpu time variables.
  ofstream residual_file;
  CPUTime processor_cpu_time, total_cpu_time;

  // Other local solution variables.
  int number_of_time_steps, first_step,
      command_flag, error_flag, line_number,
      perform_explicit_time_marching;
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
    if (!batch_flag)
      cout << endl << " Reading Electrostatic2D input data file `"
	   << Input_File_Name_ptr << "'.";
    error_flag = Process_Input_Control_Parameter_File(Input_Parameters,
						      Input_File_Name_ptr,
						      command_flag);
    if (!batch_flag && !error_flag) {
      cout << Input_Parameters << endl;
    } else if (error_flag) {
      cout << endl << "Electrostatic2D ERROR: During processing of input parameters." << endl << endl;
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
   * Create initial mesh and allocate Electrostatic2D solution        *
   * variables for specified IBVP/BVP problem.                        *
   ********************************************************************/

 execute_new_calculation: ;

  // Synchronize processors.
  CFFC_Barrier_MPI();

  // Create initial mesh.  Read mesh from grid definition or data  
  // files when specified by input parameters.

  // The primary MPI processor creates the initial mesh.
  if (CFFC_Primary_MPI_Processor()) {

    if (!batch_flag) 
      cout << endl << " Creating (or reading) initial quadrilateral multi-block mesh.";

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
      cout << endl << " Electrostatic2D ERROR: Unable to create valid Electrostatic2D multi-block mesh." << endl;
    } else {
//        if (Input_Parameters.i_Grid == GRID_NASA_ROTOR_37) {
//  	if (!batch_flag) cout << endl << " Writing geometry and flow field data files for NASA Rotor 37.";
//  	Input_Parameters.NASA_Rotor37.outputTP_UpstreamFlowConditions("NASARotor37_UpstreamProfile.dat");
//  	Input_Parameters.NASA_Rotor37.outputTP_DownstreamFlowConditions("NASARotor37_DownstreamProfile.dat");
//  	Input_Parameters.NASA_Rotor37.outputTP_Geometry("NASARotor37_Geometry.dat");
//        } else if (Input_Parameters.i_Grid == GRID_NASA_ROTOR_67) {
//  	if (!batch_flag) cout << endl << " Writing geometry and flow field data files for NASA Rotor 67.";
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
  // (allocate) array of local 2D Electrostatic equation solution blocks, 
  // assign and create (allocate) 2D Electrostatic equation solution blocks
  // corresponding to the initial mesh.
  if (!batch_flag) cout << endl << " Creating multi-block quadtree data structure and assigning"
                        << endl << " Electrostatic2D solution blocks corresponding to initial mesh.";
  Local_SolnBlk = CreateInitialSolutionBlocks(MeshBlk,
					      Local_SolnBlk,
					      Input_Parameters,
					      QuadTree,
					      List_of_Global_Solution_Blocks,
					      List_of_Local_Solution_Blocks);
  if (Local_SolnBlk == NULL) return 1;

  /********************************************************************
   * Initialize Electrostatic2D solution variables.                   *
   ********************************************************************/

  // Set the initial time level.
  Time = ZERO;
  number_of_time_steps = 0;

  // Set the CPU time to zero.
  processor_cpu_time.zero();
  total_cpu_time.zero();

  // Initialize the conserved and primitive state solution variables.
  if (!batch_flag) cout << endl << " Prescribing Electrostatic2D initial data.";
  if (Input_Parameters.i_ICs == IC_RESTART) {
    if (!batch_flag) cout << endl << " Reading Electrostatic2D solution from restart data files.";
    // Read the quadtree restart file.
    error_flag = Read_QuadTree(QuadTree,
			       List_of_Global_Solution_Blocks,
			       List_of_Local_Solution_Blocks,
			       Input_Parameters);
    if (error_flag) {
      cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D quadtree data file "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;
    // Allocate the message buffers.
    Allocate_Message_Buffers(List_of_Local_Solution_Blocks,
			     NUM_VAR_ELECTROSTATIC2D+NUM_COMP_VECTOR2D);
    // Read the solution block restart files.
    error_flag = Read_Restart_Solution(Local_SolnBlk,
				       List_of_Local_Solution_Blocks,
				       Input_Parameters,
				       number_of_time_steps,
				       Time,
				       processor_cpu_time);
    if (error_flag) {
      cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D restart "
	   << "input data file(s) on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
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
						  NUM_VAR_ELECTROSTATIC2D,
						  OFF);
  if (error_flag) {
    cout << endl << " Electrostatic2D ERROR: Message passing error during Electrostatic2D "
	 << "solution intialization on processor "
	 << List_of_Local_Solution_Blocks.ThisCPU
	 << "." << endl;
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Prescribe boundary data consistent with initial data.
  BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

  // Perform uniform, boundary, and initial mesh refinement.
  if (Input_Parameters.i_ICs != IC_RESTART) {

    // Perform uniform mesh refinement.
    if (!batch_flag) cout << endl << " Performing Electrostatic2D uniform mesh refinement.";
    error_flag = Uniform_AMR(Local_SolnBlk,
			     Input_Parameters,
			     QuadTree,
			     List_of_Global_Solution_Blocks,
			     List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << endl << " Electrostatic2D ERROR: Uniform AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Perform boundary mesh refinement.
    if (!batch_flag) cout << endl << " Performing Electrostatic2D boundary mesh refinement.";
    error_flag = Boundary_AMR(Local_SolnBlk,
			      Input_Parameters,
			      QuadTree,
			      List_of_Global_Solution_Blocks,
			      List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << endl << " Electrostatic2D ERROR: Boundary AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Perform initial mesh refinement.
    if (!batch_flag) cout << endl << " Performing Electrostatic2D initial mesh refinement.";
    error_flag = Initial_AMR(Local_SolnBlk,
			     Input_Parameters,
			     QuadTree,
			     List_of_Global_Solution_Blocks,
			     List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << endl << " Electrostatic2D ERROR: Initial AMR error on processor "
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
						    NUM_VAR_ELECTROSTATIC2D,
						    OFF);
    if (error_flag) {
      cout << endl << " Electrostatic2D ERROR: Message passing error "
	   << "during Electrostatic2D solution intialization on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;
  }

  // Output multi-block solution-adaptive quadrilateral mesh statistics.
  if (!batch_flag) {
    cout << endl
	 << endl << " Multi-block solution-adaptive quadrilateral mesh statistics: ";
    cout << endl << "  -> Number of Root Blocks i-direction: "
	 << QuadTree.NRi;
    cout << endl << "  -> Number of Root Blocks j-direction: "
	 << QuadTree.NRj;
    cout << endl << "  -> Total Number of Used Blocks: "
	 << QuadTree.countUsedBlocks();
    cout << endl << "  -> Total Number of Computational Cells: "
	 << QuadTree.countUsedCells();
    cout << endl << "  -> Number of Mesh Refinement Levels: "
	 << QuadTree.highestRefinementLevel()+1;
    cout << endl << "  -> Refinement Efficiency: "
	 << QuadTree.efficiencyRefinement() << endl; 
  }

//       if (CFFC_Primary_MPI_Processor()) {
//          for (int j_blk = 0; j_blk < QuadTree.Nblk; j_blk++) {
//             for (int i_blk = 0; i_blk < QuadTree.Ncpu; i_blk++) {
//                if (QuadTree.Blocks[i_blk][j_blk] != NULL) {
//                   cout << endl << " cpu = " 
//                        << i_blk
//                        << " blk = "
//                        << j_blk
//                        << " blk = "
//                        << QuadTree.Blocks[i_blk][j_blk]->block;
//                } else {
//                   cout << endl << " cpu = " 
//                        << i_blk
//                        << " blk = "
//                        << j_blk;
//                }
//             }
//          }
//       }

  /********************************************************************
   * Solve IBVP or BVP for conservation form of 2D Electrostatic      *
   * equations on multi-block solution-adaptive quadrilateral mesh.   *
   ********************************************************************/

 continue_existing_calculation: ;

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Open residual file.
  first_step = 1;

  if (CFFC_Primary_MPI_Processor()) {
    error_flag = Open_Progress_File(residual_file,
				    Input_Parameters.Output_File_Name,
				    number_of_time_steps);
    if (error_flag) {
      cout << endl << " Electrostatic2D ERROR: Unable to open residual file for Electrostatic2D calculation." << endl;
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

    if (!batch_flag) cout << endl << " Beginning Electrostatic2D computations on "
			  << Date_And_Time() << "." << endl << endl;

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
	  if (!batch_flag) cout << endl << endl << " Refining Grid.  Performing adaptive mesh refinement at n = "
				<< number_of_time_steps << ".";
	  error_flag = AMR(Local_SolnBlk,
			   Input_Parameters,
			   QuadTree,
			   List_of_Global_Solution_Blocks,
			   List_of_Local_Solution_Blocks,
			   ON,ON);
	  if (error_flag) {
	    cout << endl << " Electrostatic2D ERROR: Electrostatic2D AMR error number "
		 << error_flag << " on processor "
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
	    cout << endl << " New multi-block solution-adaptive quadrilateral mesh statistics: ";
	    cout << endl << "  -> Number of Root Blocks i-direction: "
		 << QuadTree.NRi;
	    cout << endl << "  -> Number of Root Blocks j-direction: "
		 << QuadTree.NRj;
	    cout << endl << "  -> Total Number of Used Blocks: "
		 << QuadTree.countUsedBlocks();
	    cout << endl << "  -> Total Number of Computational Cells: "
		 << QuadTree.countUsedCells();
	    cout << endl << "  -> Number of Mesh Refinement Levels: "
		 << QuadTree.highestRefinementLevel()+1;
	    cout << endl << "  -> Refinement Efficiency: "
		 << QuadTree.efficiencyRefinement() << endl;
	  }
	}
      }

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
						  LIMITER_ONE) > 
		   Input_Parameters.Time_Max) {
	  dTime = (Input_Parameters.Time_Max-Time)/
 	          (Input_Parameters.CFL_Number*
	          MultiStage_Optimally_Smoothing(Input_Parameters.N_Stage,
						 Input_Parameters.N_Stage,
						 LIMITER_ONE));
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

      // Update CPU time used for the calculation so far.
      processor_cpu_time.update();
      total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput);

      // Periodically save restart solution files.
      if (!first_step &&
	  number_of_time_steps-Input_Parameters.Restart_Solution_Save_Frequency*
	  (number_of_time_steps/Input_Parameters.Restart_Solution_Save_Frequency) == 0) {
	if (!batch_flag) 
	  cout << endl << endl << "  Saving Electrostatic2D solution to restart data file(s) after"
	       << " n = " << number_of_time_steps << " steps (iterations).";
	// Write the quadtree restart file.
	error_flag = Write_QuadTree(QuadTree,Input_Parameters);
	if (error_flag) {
	  cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D quadtree data file "
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
	  cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D restart output data file(s) "
	       << "on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU
	       << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;
	if (!batch_flag) cout << endl;
      }

      // Output progress information for the calculation.
      if (!batch_flag) Output_Progress_L2norm(number_of_time_steps,
					      Time*THOUSAND,
					      total_cpu_time,
					      residual_l2_norm,
					      first_step,
					      50);
//	if (!batch_flag) Output_Progress(number_of_time_steps,
//					 Time*THOUSAND,
//					 total_cpu_time,
//					 residual_l1_norm,
//					 first_step,
//					 50);
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
	  number_of_time_steps >= 
	  Input_Parameters.Maximum_Number_of_Time_Steps) break;
      if (Input_Parameters.Time_Accurate &&
	  Time >= Input_Parameters.Time_Max) break;

      // Update solution for next time step using a multistage
      // time stepping scheme.
      for (int i_stage = 1; i_stage <= Input_Parameters.N_Stage; i_stage++) {

	// Step 1. Exchange solution information between neighbouring blocks.
	if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk,
							List_of_Local_Solution_Blocks,
							NUM_VAR_ELECTROSTATIC2D,
							OFF);
	if (error_flag) {
	  cout << endl << " Electrostatic2D ERROR: Electrostatic2D message passing error on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU
	       << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;

	// Step 2. Apply boundary conditions for stage.
	BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

	// Step 3. Determine solution residuals for stage.
	error_flag = dUdt_Multistage_Explicit(Local_SolnBlk,
					      List_of_Local_Solution_Blocks,
					      Input_Parameters,
					      i_stage);
	if (error_flag) {
	  cout << endl << " Electrostatic2D ERROR: Electrostatic2D solution error on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU
	       << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;

	// Step 4. Send boundary flux corrections at block interfaces with resolution changes.
	error_flag = Send_Conservative_Flux_Corrections(Local_SolnBlk,
							List_of_Local_Solution_Blocks,
							NUM_VAR_ELECTROSTATIC2D);
	if (error_flag) {
	  cout << endl << " Electrostatic2D ERROR: Electrostatic2D flux correction message passing error on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU
	       << "." << endl;
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
	  cout << endl << " Electrostatic2D ERROR: Electrostatic2D solution update error on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU
	       << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return error_flag;

      }

      // Update time and time step counter.
      if (first_step) first_step = 0;
      number_of_time_steps++;
      if (Input_Parameters.i_Time_Integration != 
	  TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) {
	Time += Input_Parameters.CFL_Number*dTime;
      } else {
	Time += Input_Parameters.CFL_Number*dTime*
	        MultiStage_Optimally_Smoothing(Input_Parameters.N_Stage,
					       Input_Parameters.N_Stage,
					       LIMITER_ONE);
      }

    }

    if (!batch_flag) cout << endl << endl << " Electrostatic2D computations complete on " 
			  << Date_And_Time() << "." << endl;
  }

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
						  NUM_VAR_ELECTROSTATIC2D,
						  OFF);
  if (error_flag) {
    cout << endl << " Electrostatic2D ERROR: Electrostatic2D message passing error on processor "
	 << List_of_Local_Solution_Blocks.ThisCPU
	 << "." << endl;
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Apply boundary conditions.
  BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

  error_flag = Determine_Electric_Field(Local_SolnBlk,
					List_of_Local_Solution_Blocks,
					Input_Parameters);
  if (error_flag) {
    cout << endl << " Electrostatic2D ERROR: Electrostatic2D error during calculation of electric field on processor "
	 << List_of_Local_Solution_Blocks.ThisCPU
	 << "." << endl;
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Close residual file.
  if (CFFC_Primary_MPI_Processor()) error_flag = Close_Progress_File(residual_file);

  /********************************************************************
   * Solution calculations complete.  Write the 2D Electrostatic      *
   * solution to output and restart files as required, reset solution *
   * parameters, and run other cases as specified by input            *
   * parameters.                                                      *
   ********************************************************************/

 postprocess_current_calculation: ;

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  while (1) {
    if (CFFC_Primary_MPI_Processor()) {
      Get_Next_Input_Control_Parameter(Input_Parameters);
      command_flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
      line_number = Input_Parameters.Line_Number;
      if (command_flag == INVALID_INPUT_CODE ||
	  command_flag == INVALID_INPUT_VALUE) {
	line_number = -line_number;
	cout << endl << " Electrostatic2D ERROR: Error reading Electrostatic2D data at line #"
	     << -line_number << " of input data file." << endl;
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
	// Set flow geometry indicator (planar/axisymmetric).
	Local_SolnBlk[nb].Axisymmetric = Input_Parameters.Axisymmetric;
      }
    }

    if (command_flag == EXECUTE_CODE) {
      // Deallocate memory for 2D Electrostatic equation solution.
      if (!batch_flag) cout << endl << " Deallocating Electrostatic2D solution variables.";
      Local_SolnBlk = Deallocate(Local_SolnBlk,Input_Parameters);
      //Deallocate_Message_Buffers(List_of_Local_Solution_Blocks); // Not necessary here!
      List_of_Local_Solution_Blocks.deallocate();
      List_of_Global_Solution_Blocks.deallocate();
      QuadTree.deallocate();
      MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk,
					    Input_Parameters.Number_of_Blocks_Idir,
					    Input_Parameters.Number_of_Blocks_Jdir);
      // Output input parameters for new caluculation.
      if (!batch_flag) cout << endl << endl << " Starting a new calculation." << Input_Parameters << endl;
      // Execute new calculation.
      goto execute_new_calculation;

    } else if (command_flag == TERMINATE_CODE) {

      CFFC_Barrier_MPI();

      // Deallocate memory for 2D Electrostatic equation solution.
      if (!batch_flag) cout << endl << " Deallocating Electrostatic2D solution variables.";
      Local_SolnBlk = Deallocate(Local_SolnBlk,Input_Parameters);
      //Deallocate_Message_Buffers(List_of_Local_Solution_Blocks); // Not necessary here!
      List_of_Local_Solution_Blocks.deallocate();
      List_of_Global_Solution_Blocks.deallocate();
      QuadTree.deallocate();
      MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk,
					    Input_Parameters.Number_of_Blocks_Idir,
					    Input_Parameters.Number_of_Blocks_Jdir);
      // Close input data file.
      if (!batch_flag) cout << endl << endl << " Closing Electrostatic2D input data file.";
      if (CFFC_Primary_MPI_Processor()) Close_Input_File(Input_Parameters);
      // Terminate calculation.
      return 0;

    } else if (command_flag == CONTINUE_CODE) {
      // Reset maximum time step counter.
      Input_Parameters.Maximum_Number_of_Time_Steps += number_of_time_steps;
      // Output input parameters for continuing calculation.
      if (!batch_flag) cout << endl << endl << " Continuing existing calculation."
			    << Input_Parameters << endl;
      // Continue existing calculation.
      goto continue_existing_calculation;

    } else if (command_flag == REFINE_GRID_CODE) {
      // Refine mesh using block based adaptive mesh refinement algorithm.
      if (!batch_flag) cout << endl << " Refining Grid.  Performing adaptive mesh refinement.";
      error_flag = AMR(Local_SolnBlk,
		       Input_Parameters,
		       QuadTree,
		       List_of_Global_Solution_Blocks,
		       List_of_Local_Solution_Blocks,
		       ON,ON);
      if (error_flag) {
	cout << endl << " Electrostatic2D ERROR: Electrostatic2D AMR error on processor "
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
	cout << endl
	     << endl << " New multi-block solution-adaptive quadrilateral mesh statistics: ";
	cout << endl << "  -> Number of Root Blocks i-direction: "
	     << QuadTree.NRi;
	cout << endl << "  -> Number of Root Blocks j-direction: "
	     << QuadTree.NRj;
	cout << endl << "  -> Total Number of Used Blocks: "
	     << QuadTree.countUsedBlocks();
	cout << endl << "  -> Total Number of Computational Cells: "
	     << QuadTree.countUsedCells();
	cout << endl << "  -> Number of Mesh Refinement Levels: "
	     << QuadTree.highestRefinementLevel()+1;
	cout << endl << "  -> Refinement Efficiency: "
	     << QuadTree.efficiencyRefinement() << endl;
      }

//       if (CFFC_Primary_MPI_Processor()) {
//          for (int j_blk = 0; j_blk < QuadTree.Nblk; j_blk++) {
//             for (int i_blk = 0; i_blk < QuadTree.Ncpu; i_blk++) {
//                if (QuadTree.Blocks[i_blk][j_blk] != NULL) {
//                   cout << endl << " cpu = " 
//                        << i_blk
//                        << " blk = "
//                        << j_blk
//                        << " blk = "
//                        << QuadTree.Blocks[i_blk][j_blk]->block;
//                } else {
//                   cout << endl << " cpu = " 
//                        << i_blk
//                        << " blk = "
//                        << j_blk;
//                }
//             }
//          }
//       }

      // Apply boundary conditions.
      BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

    } else if (command_flag == WRITE_OUTPUT_CODE) {
      // Output solution data.
      if (!batch_flag) cout << endl << " Writing Electrostatic2D solution to output data file(s).";
      error_flag = Output_Tecplot(Local_SolnBlk,
				  List_of_Local_Solution_Blocks,
				  Input_Parameters,
				  number_of_time_steps,
				  Time);
      if (error_flag) {
	cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_CELLS_CODE) {
      // Output solution data.
      if (!batch_flag) cout << endl << " Writing cell-centered Electrostatic2D solution to output data file(s).";
      error_flag = Output_Cells_Tecplot(Local_SolnBlk,
					List_of_Local_Solution_Blocks,
					Input_Parameters,
					number_of_time_steps,
					Time);
      if (error_flag) {
	cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D cell output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;
      
    } else if (command_flag == WRITE_OUTPUT_QUASI3D_CODE) {
      // Output solution data.
      if (!batch_flag) cout << endl << " Writing Electrostatic2D quasi3D solution to output data file(s).";
      error_flag = Output_Quasi3D_Tecplot(Local_SolnBlk,
					  List_of_Local_Solution_Blocks,
					  Input_Parameters,
					  number_of_time_steps,
					  Time);
      if (error_flag) {
	cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D quasi3D output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_RESTART_CODE) {
      // Write restart files.
      if (!batch_flag) cout << endl << " Writing Electrostatic2D solution to restart data file(s).";
      // Write the quadtree restart file.
      error_flag = Write_QuadTree(QuadTree,Input_Parameters);
      if (error_flag) {
	cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D quadtree data file "
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
	cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D restart output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_GRID_CODE) {
      // Output multi-block solution-adaptive mesh data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << endl << " Writing Electrostatic2D multi-block mesh to grid data output file.";
	error_flag = Output_Tecplot(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D mesh data output file." << endl;
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      
    } else if (command_flag == WRITE_GRID_DEFINITION_CODE) {
      // Write multi-block solution-adaptive mesh definition files.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << endl << " Writing Electrostatic2D multi-block mesh to grid definition files.";
	error_flag = Write_Multi_Block_Grid_Definition(MeshBlk,
						       Input_Parameters);
	error_flag = Write_Multi_Block_Grid(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D multi-block mesh definition files." << endl;
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      
    } else if (command_flag == WRITE_OUTPUT_GRID_NODES_CODE) {
      // Output multi-block solution-adaptive mesh node data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << endl << " Writing Electrostatic2D multi-block mesh to node data output file.";
	error_flag = Output_Nodes_Tecplot(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D mesh node data output file." << endl;
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      
    } else if (command_flag == WRITE_OUTPUT_GRID_CELLS_CODE) {
      // Output multi-block solution-adaptive mesh cell data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << endl << " Writing Electrostatic2D multi-block mesh to cell data output file.";
	error_flag = Output_Cells_Tecplot(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D mesh cell data output file." << endl;
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == INVALID_INPUT_CODE ||
	       command_flag == INVALID_INPUT_VALUE) {
      line_number = -line_number;
      cout << endl << " Electrostatic2D ERROR: Error reading Electrostatic2D data at line #"
	   << -line_number  << " of input data file." << endl;
      return line_number;
    }
    
  }
  
  /********************************************************************
   * End of all Electrostatic2DSolver computations and I/O.           *
   ********************************************************************/
  return 0;
  
}

/**********************************************************************
 * Routine: Electrostatic2DQuadSolver                                 *
 *                                                                    *
 * Computes solutions to 2D Electrostatic equations on 2D             *
 * quadrilateral multi-block solution-adaptive mesh.                  *
 *                                                                    *
 **********************************************************************/
Electrostatic2D_Quad_Block* Electrostatic2DQuadSolver(char *Input_File_Name,
						      int batch_flag,
						      Electrostatic2D_Quad_Block *Local_SolnBlk,
						      Electrostatic2D_Input_Parameters &Input_Parameters,
						      QuadTreeBlock_DataStructure &QuadTree,
						      AdaptiveBlockResourceList &List_of_Global_Solution_Blocks,
						      AdaptiveBlock2D_List &List_of_Local_Solution_Blocks) {

  /********************************************************************
   * Local variable declarations.                                     *
   ********************************************************************/

  // 2D multi-block grid blocks for mesh generation:
  Grid2D_Quad_Block **MeshBlk;

  // Define residual file and cpu time variables.
  ofstream residual_file;
  CPUTime processor_cpu_time, total_cpu_time;

  // Electrostatic2D input file name.
  char Input_File_Name_ptr[256];

  // Other local solution variables.
  int number_of_time_steps, first_step,
      command_flag, error_flag, line_number,
      perform_explicit_time_marching;
  double Time, dTime;
  double residual_l2norm_first, residual_ratio;
  double residual_l1_norm, residual_l2_norm, residual_max_norm;

  /********************************************************************  
   * Set default values for the input solution parameters and then    *
   * read user specified input values from the specified input        *
   * parameter file.                                                  *
   ********************************************************************/

  // Determine LevelSet2D input file name.
  strcpy(Input_File_Name_ptr,"electrostatic_");
  strcat(Input_File_Name_ptr,Input_File_Name);

  // The primary MPI processor processes the input parameter file.
  if (CFFC_Primary_MPI_Processor()) {
    if (!batch_flag)
      cout << endl << " Reading Electrostatic2D input data file `"
	   << Input_File_Name_ptr << "'.";
    error_flag = Process_Input_Control_Parameter_File(Input_Parameters,
						      Input_File_Name_ptr,
						      command_flag);
    if (!batch_flag && !error_flag) {
      cout << Input_Parameters << endl;
    } else if (error_flag) {
      cout << endl << "Electrostatic2D ERROR: During processing of input parameters." << endl << endl;
    }
  } else {
    error_flag = 0;
  }
  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Broadcast input solution parameters to other MPI processors.
  CFFC_Broadcast_MPI(&error_flag,1);
  if (error_flag) return NULL;
  CFFC_Broadcast_MPI(&command_flag,1);
  if (command_flag == TERMINATE_CODE) return NULL;
  Broadcast_Input_Parameters(Input_Parameters);

  /********************************************************************
   * Create initial mesh and allocate Electrostatic2D solution        *
   * variables for specified IBVP/BVP problem.                        *
   ********************************************************************/

 execute_new_calculation: ;

  // Synchronize processors.
  CFFC_Barrier_MPI();

  // Create initial mesh.  Read mesh from grid definition or data  
  // files when specified by input parameters.

  // The primary MPI processor creates the initial mesh.
  if (CFFC_Primary_MPI_Processor()) {

    if (!batch_flag) 
      cout << endl << " Creating (or reading) initial quadrilateral multi-block mesh.";

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
      cout << endl << " Electrostatic2D ERROR: Unable to create valid Electrostatic2D multi-block mesh." << endl;
      cout.flush();
    } else {
//        if (Input_Parameters.i_Grid == GRID_NASA_ROTOR_37) {
//  	if (!batch_flag) cout << endl << " Writing geometry and flow field data files for NASA Rotor 37.";
//  	Input_Parameters.NASA_Rotor37.outputTP_UpstreamFlowConditions("NASARotor37_UpstreamProfile.dat");
//  	Input_Parameters.NASA_Rotor37.outputTP_DownstreamFlowConditions("NASARotor37_DownstreamProfile.dat");
//  	Input_Parameters.NASA_Rotor37.outputTP_Geometry("NASARotor37_Geometry.dat");
//        } else if (Input_Parameters.i_Grid == GRID_NASA_ROTOR_67) {
//  	if (!batch_flag) cout << endl << " Writing geometry and flow field data files for NASA Rotor 67.";
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
  if (error_flag) return NULL;
  MeshBlk = Broadcast_Multi_Block_Grid(MeshBlk,Input_Parameters);

  // Create (allocate) multi-block quadtree data structure, create
  // (allocate) array of local 2D Electrostatic equation solution blocks, 
  // assign and create (allocate) 2D Electrostatic equation solution blocks
  // corresponding to the initial mesh.
  if (!batch_flag) cout << endl << " Creating multi-block quadtree data structure and assigning"
                        << endl << " Electrostatic2D solution blocks corresponding to initial mesh.";
  Local_SolnBlk = CreateInitialSolutionBlocks(MeshBlk,
					      Local_SolnBlk,
					      Input_Parameters,
					      QuadTree,
					      List_of_Global_Solution_Blocks,
					      List_of_Local_Solution_Blocks);
  if (Local_SolnBlk == NULL) return NULL;

  /********************************************************************
   * Initialize Electrostatic2D solution variables.                   *
   ********************************************************************/

  // Set the initial time level.
  Time = ZERO;
  number_of_time_steps = 0;

  // Set the CPU time to zero.
  processor_cpu_time.zero();
  total_cpu_time.zero();

  // Initialize the conserved and primitive state solution variables.
  if (!batch_flag) cout << endl << " Prescribing Electrostatic2D initial data.";
  if (Input_Parameters.i_ICs == IC_RESTART) {
    if (!batch_flag) cout << endl << " Reading Electrostatic2D solution from restart data files.";
    // Read the quadtree restart file.
    error_flag = Read_QuadTree(QuadTree,
			       List_of_Global_Solution_Blocks,
			       List_of_Local_Solution_Blocks,
			       Input_Parameters);
    if (error_flag) {
      cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D quadtree data file "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return NULL;
    // Allocate the message buffers.
    Allocate_Message_Buffers(List_of_Local_Solution_Blocks,
			     NUM_VAR_ELECTROSTATIC2D+NUM_COMP_VECTOR2D);
    // Read the solution block restart files.
    error_flag = Read_Restart_Solution(Local_SolnBlk,
				       List_of_Local_Solution_Blocks,
				       Input_Parameters,
				       number_of_time_steps,
				       Time,
				       processor_cpu_time);
    if (error_flag) {
      cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D restart "
	   << "input data file(s) on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return NULL;
    // Ensure each processor has the correct time and time.
    number_of_time_steps = CFFC_Maximum_MPI(number_of_time_steps);
    Time = CFFC_Maximum_MPI(Time);
    processor_cpu_time.cput = CFFC_Maximum_MPI(processor_cpu_time.cput);
    Input_Parameters.Maximum_Number_of_Time_Steps = CFFC_Maximum_MPI(Input_Parameters.Maximum_Number_of_Time_Steps);

    // Synchronize processors.
    CFFC_Barrier_MPI();
    // Broadcast input solution parameters to other MPI processors.
    CFFC_Broadcast_MPI(&error_flag,1);
    if (error_flag != 0) return NULL;
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
						  NUM_VAR_ELECTROSTATIC2D,
						  OFF);
  if (error_flag) {
    cout << endl << " Electrostatic2D ERROR: Message passing error during Electrostatic2D "
	 << "solution intialization on processor "
	 << List_of_Local_Solution_Blocks.ThisCPU
	 << "." << endl;
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return NULL;

  // Prescribe boundary data consistent with initial data.
  BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

  // Perform uniform, boundary, and initial mesh refinement.
  if (Input_Parameters.i_ICs != IC_RESTART) {

    // Perform uniform mesh refinement.
    if (!batch_flag) cout << endl << " Performing Electrostatic2D uniform mesh refinement.";
    error_flag = Uniform_AMR(Local_SolnBlk,
			     Input_Parameters,
			     QuadTree,
			     List_of_Global_Solution_Blocks,
			     List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << endl << " Electrostatic2D ERROR: Uniform AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return NULL;

    // Perform boundary mesh refinement.
    if (!batch_flag) cout << endl << " Performing Electrostatic2D boundary mesh refinement.";
    error_flag = Boundary_AMR(Local_SolnBlk,
			      Input_Parameters,
			      QuadTree,
			      List_of_Global_Solution_Blocks,
			      List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << endl << " Electrostatic2D ERROR: Boundary AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return NULL;

    // Perform initial mesh refinement.
    if (!batch_flag) cout << endl << " Performing Electrostatic2D initial mesh refinement.";
    error_flag = Initial_AMR(Local_SolnBlk,
			     Input_Parameters,
			     QuadTree,
			     List_of_Global_Solution_Blocks,
			     List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << endl << " Electrostatic2D ERROR: Initial AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
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
    if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk,
						    List_of_Local_Solution_Blocks,
						    NUM_VAR_ELECTROSTATIC2D,
						    OFF);
    if (error_flag) {
      cout << endl << " Electrostatic2D ERROR: Message passing error during Electrostatic2D solution intialization "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << "." << endl;
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return NULL;
  }

  // Output multi-block solution-adaptive quadrilateral mesh statistics.
  if (!batch_flag) {
    cout << endl << endl << " Multi-block solution-adaptive quadrilateral mesh statistics: ";
    cout << endl << "  -> Number of Root Blocks i-direction: "
	 << QuadTree.NRi;
    cout << endl << "  -> Number of Root Blocks j-direction: "
	 << QuadTree.NRj;
    cout << endl << "  -> Total Number of Used Blocks: "
	 << QuadTree.countUsedBlocks();
    cout << endl << "  -> Total Number of Computational Cells: "
	 << QuadTree.countUsedCells();
    cout << endl << "  -> Number of Mesh Refinement Levels: "
	 << QuadTree.highestRefinementLevel()+1;
    cout << endl << "  -> Refinement Efficiency: "
	 << QuadTree.efficiencyRefinement() << endl; 
  }

//       if (CFFC_Primary_MPI_Processor()) {
//          for (int j_blk = 0; j_blk < QuadTree.Nblk; j_blk++) {
//             for (int i_blk = 0; i_blk < QuadTree.Ncpu; i_blk++) {
//                if (QuadTree.Blocks[i_blk][j_blk] != NULL) {
//                   cout << endl << " cpu = " 
//                        << i_blk
//                        << " blk = "
//                        << j_blk
//                        << " blk = "
//                        << QuadTree.Blocks[i_blk][j_blk]->block;
//                } else {
//                   cout << endl << " cpu = " 
//                        << i_blk
//                        << " blk = "
//                        << j_blk;
//                }
//             }
//          }
//       }

  /********************************************************************
   * Solve IBVP or BVP for conservation form of 2D Electrostatic      *
   * equations on multi-block solution-adaptive quadrilateral mesh.   *
   ********************************************************************/

 continue_existing_calculation: ;

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Open residual file.
  first_step = 1;

  if (CFFC_Primary_MPI_Processor()) {
    error_flag = Open_Progress_File(residual_file,
				    Input_Parameters.Output_File_Name,
				    number_of_time_steps);
    if (error_flag) {
      cout << endl << " Electrostatic2D ERROR: Unable to open residual file for Electrostatic2D calculation." << endl;
      cout.flush();
    }
  }
  // MPI barrier to ensure processor synchronization.  
  CFFC_Barrier_MPI();
  CFFC_Broadcast_MPI(&error_flag,1);
  if (error_flag) return NULL;
  // Reset the CPU time.
  processor_cpu_time.reset();

  // Perform required number of iterations (time steps).
  if ((!Input_Parameters.Time_Accurate &&
       Input_Parameters.Maximum_Number_of_Time_Steps > 0 &&
       number_of_time_steps < Input_Parameters.Maximum_Number_of_Time_Steps) ||
      (Input_Parameters.Time_Accurate &&
       Input_Parameters.Time_Max > Time)) {

    if (!batch_flag) cout << endl << " Beginning Electrostatic2D computations on "
			  << Date_And_Time() << "." << endl << endl;

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
	  if (!batch_flag) cout << endl << endl << " Refining Grid.  Performing adaptive mesh refinement at n = "
				<< number_of_time_steps << ".";
	  error_flag = AMR(Local_SolnBlk,
			   Input_Parameters,
			   QuadTree,
			   List_of_Global_Solution_Blocks,
			   List_of_Local_Solution_Blocks,
			   ON,ON);
	  if (error_flag) {
	    cout << endl << " Electrostatic2D ERROR: Electrostatic2D AMR error number "
		 << error_flag << " on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
	  }
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) {
	    command_flag = Output_Tecplot(Local_SolnBlk,
					  List_of_Local_Solution_Blocks,
					  Input_Parameters,
					  number_of_time_steps,
					  Time);
	    return NULL;
	  }
	  if (!batch_flag) {
	    cout << endl << " New multi-block solution-adaptive quadrilateral mesh statistics: ";
	    cout << endl << "  -> Number of Root Blocks i-direction: "
		 << QuadTree.NRi;
	    cout << endl << "  -> Number of Root Blocks j-direction: "
		 << QuadTree.NRj;
	    cout << endl << "  -> Total Number of Used Blocks: "
		 << QuadTree.countUsedBlocks();
	    cout << endl << "  -> Total Number of Computational Cells: "
		 << QuadTree.countUsedCells();
	    cout << endl << "  -> Number of Mesh Refinement Levels: "
		 << QuadTree.highestRefinementLevel()+1;
	    cout << endl << "  -> Refinement Efficiency: "
		 << QuadTree.efficiencyRefinement() << endl;
	  }
	}
      }

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
						  LIMITER_ONE) > 
		   Input_Parameters.Time_Max) {
	  dTime = (Input_Parameters.Time_Max-Time)/
 	          (Input_Parameters.CFL_Number*
	          MultiStage_Optimally_Smoothing(Input_Parameters.N_Stage,
						 Input_Parameters.N_Stage,
						 LIMITER_ONE));
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

      // Update CPU time used for the calculation so far.
      processor_cpu_time.update();
      total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput);

      // Periodically save restart solution files.
      if (!first_step &&
	  number_of_time_steps-Input_Parameters.Restart_Solution_Save_Frequency*
	  (number_of_time_steps/Input_Parameters.Restart_Solution_Save_Frequency) == 0) {
	if (!batch_flag) 
	  cout << endl << endl << "  Saving Electrostatic2D solution to restart data file(s) after"
	       << " n = " << number_of_time_steps << " steps (iterations).";
	// Write the quadtree restart file.
	error_flag = Write_QuadTree(QuadTree,Input_Parameters);
	if (error_flag) {
	  cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D quadtree data file "
	       << "on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU
	       << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return NULL;
	// Write the solution block restart files.
	error_flag = Write_Restart_Solution(Local_SolnBlk,
					    List_of_Local_Solution_Blocks,
					    Input_Parameters,
					    number_of_time_steps,
					    Time,
					    processor_cpu_time);
	if (error_flag) {
	  cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D restart output data file(s) "
	       << "on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU
	       << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return NULL;
	if (!batch_flag) cout << endl;
      }

      // Output progress information for the calculation.
      if (!batch_flag) Output_Progress_L2norm(number_of_time_steps,
					      Time*THOUSAND,
					      total_cpu_time,
					      residual_l2_norm,
					      first_step,
					      50);
//	if (!batch_flag) Output_Progress(number_of_time_steps,
//					 Time*THOUSAND,
//					 total_cpu_time,
//					 residual_l1_norm,
//					 first_step,
//					 50);
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
	  number_of_time_steps >= 
	  Input_Parameters.Maximum_Number_of_Time_Steps) break;
      if (Input_Parameters.Time_Accurate &&
	  Time >= Input_Parameters.Time_Max) break;

      // Update solution for next time step using a multistage
      // time stepping scheme.
      for (int i_stage = 1; i_stage <= Input_Parameters.N_Stage; i_stage++) {

	// Step 1. Exchange solution information between neighbouring blocks.
	if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk,
							List_of_Local_Solution_Blocks,
							NUM_VAR_ELECTROSTATIC2D,
							OFF);
	if (error_flag) {
	  cout << endl << " Electrostatic2D ERROR: Electrostatic2D message passing error on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU
	       << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return NULL;

	// Step 2. Apply boundary conditions for stage.
	BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

	// Step 3. Determine solution residuals for stage.
	error_flag = dUdt_Multistage_Explicit(Local_SolnBlk,
					      List_of_Local_Solution_Blocks,
					      Input_Parameters,
					      i_stage);
	if (error_flag) {
	  cout << endl << " Electrostatic2D ERROR: Electrostatic2D solution error on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU
	       << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return NULL;

	// Step 4. Send boundary flux corrections at block interfaces with resolution changes.
	error_flag = Send_Conservative_Flux_Corrections(Local_SolnBlk,
							List_of_Local_Solution_Blocks,
							NUM_VAR_ELECTROSTATIC2D);
	if (error_flag) {
	  cout << endl << " Electrostatic2D ERROR: Electrostatic2D flux correction message passing error on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU
	       << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return NULL;

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
	  cout << endl << " Electrostatic2D ERROR: Electrostatic2D solution update error on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU
	       << "." << endl;
	}
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return NULL;

      }

      // Update time and time step counter.
      if (first_step) first_step = 0;
      number_of_time_steps++;
      if (Input_Parameters.i_Time_Integration != 
	  TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) {
	Time += Input_Parameters.CFL_Number*dTime;
      } else {
	Time += Input_Parameters.CFL_Number*dTime*
	        MultiStage_Optimally_Smoothing(Input_Parameters.N_Stage,
					       Input_Parameters.N_Stage,
					       LIMITER_ONE);
      }

    }

    if (!batch_flag) cout << endl << endl << " Electrostatic2D computations complete on " 
			  << Date_And_Time() << "." << endl;
  }

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
						  NUM_VAR_ELECTROSTATIC2D,
						  OFF);
  if (error_flag) {
    cout << endl << " Electrostatic2D ERROR: Electrostatic2D message passing error on processor "
	 << List_of_Local_Solution_Blocks.ThisCPU
	 << "." << endl;
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return NULL;

  // Apply boundary conditions.
  BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

  error_flag = Determine_Electric_Field(Local_SolnBlk,
					List_of_Local_Solution_Blocks,
					Input_Parameters);
  if (error_flag) {
    cout << endl << " Electrostatic2D ERROR: Electrostatic2D error during calculation of electric field on processor "
	 << List_of_Local_Solution_Blocks.ThisCPU
	 << "." << endl;
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return NULL;

  // Close residual file.
  if (CFFC_Primary_MPI_Processor()) error_flag = Close_Progress_File(residual_file);

  /********************************************************************
   * Solution calculations complete.  Write the 2D Electrostatic      *
   * solution to output and restart files as required, reset solution *
   * parameters, and run other cases as specified by input            *
   * parameters.                                                      *
   ********************************************************************/

 postprocess_current_calculation: ;

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  while (1) {
    if (CFFC_Primary_MPI_Processor()) {
      Get_Next_Input_Control_Parameter(Input_Parameters);
      command_flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
      line_number = Input_Parameters.Line_Number;
      if (command_flag == INVALID_INPUT_CODE ||
	  command_flag == INVALID_INPUT_VALUE) {
	line_number = -line_number;
	cout << endl << " Electrostatic2D ERROR: Error reading Electrostatic2D data at line #"
	     << -line_number << " of input data file." << endl;
	return NULL;
      }
      Reinitialize_Reference_State(Input_Parameters);
    }
    // MPI barrier to ensure processor synchronization.
    CFFC_Barrier_MPI(); 
    Broadcast_Input_Parameters(Input_Parameters);
    CFFC_Broadcast_MPI(&command_flag,1);
    for (int nb = 0; nb < List_of_Local_Solution_Blocks.Nblk; nb++) {
      if (List_of_Local_Solution_Blocks.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	// Set flow geometry indicator (planar/axisymmetric).
	Local_SolnBlk[nb].Axisymmetric = Input_Parameters.Axisymmetric;
      }
    }

    if (command_flag == EXECUTE_CODE) {
      // Deallocate memory for 2D Electrostatic equation solution.
      if (!batch_flag) cout << endl << " Deallocating Electrostatic2D solution variables.";
      Local_SolnBlk = Deallocate(Local_SolnBlk,Input_Parameters);
      //Deallocate_Message_Buffers(List_of_Local_Solution_Blocks); // Not necessary here!
      List_of_Local_Solution_Blocks.deallocate();
      List_of_Global_Solution_Blocks.deallocate();
      QuadTree.deallocate();
      MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk,
					    Input_Parameters.Number_of_Blocks_Idir,
					    Input_Parameters.Number_of_Blocks_Jdir);
      // Output input parameters for new caluculation.
      if (!batch_flag) cout << endl << endl << " Starting a new calculation." << Input_Parameters << endl;
      // Execute new calculation.
      goto execute_new_calculation;

    } else if (command_flag == TERMINATE_CODE) {

      CFFC_Barrier_MPI();

      // Deallocate memory for 2D Electrostatic equation solution.
      //if (!batch_flag) cout << endl << " Deallocating Electrostatic2D solution variables.";
      //Local_SolnBlk = Deallocate(Local_SolnBlk,Input_Parameters);
      //Deallocate_Message_Buffers(List_of_Local_Solution_Blocks); // Not necessary here!
      //List_of_Local_Solution_Blocks.deallocate();
      //List_of_Global_Solution_Blocks.deallocate();
      //QuadTree.deallocate();
      MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk,
					    Input_Parameters.Number_of_Blocks_Idir,
					    Input_Parameters.Number_of_Blocks_Jdir);
      // Close input data file.
      if (!batch_flag) cout << endl << endl << " Closing Electrostatic2D input data file." << endl;
      if (CFFC_Primary_MPI_Processor()) Close_Input_File(Input_Parameters);
      // Terminate calculation.
      return Local_SolnBlk;

    } else if (command_flag == CONTINUE_CODE) {
      // Reset maximum time step counter.
      Input_Parameters.Maximum_Number_of_Time_Steps += number_of_time_steps;
      // Output input parameters for continuing calculation.
      if (!batch_flag) cout << endl << endl << " Continuing existing calculation."
			    << Input_Parameters << endl;
      // Continue existing calculation.
      goto continue_existing_calculation;

    } else if (command_flag == REFINE_GRID_CODE) {
      // Refine mesh using block based adaptive mesh refinement algorithm.
      if (!batch_flag) cout << endl << " Refining Grid.  Performing adaptive mesh refinement.";
      error_flag = AMR(Local_SolnBlk,
		       Input_Parameters,
		       QuadTree,
		       List_of_Global_Solution_Blocks,
		       List_of_Local_Solution_Blocks,
		       ON,ON);
      if (error_flag) {
	cout << endl << " Electrostatic2D ERROR: Electrostatic2D AMR error on processor "
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
	return NULL;
      }

      // Output multi-block solution-adaptive quadrilateral mesh statistics.
      if (!batch_flag) {
	cout << endl
	     << endl << " New multi-block solution-adaptive quadrilateral mesh statistics: ";
	cout << endl << "  -> Number of Root Blocks i-direction: "
	     << QuadTree.NRi;
	cout << endl << "  -> Number of Root Blocks j-direction: "
	     << QuadTree.NRj;
	cout << endl << "  -> Total Number of Used Blocks: "
	     << QuadTree.countUsedBlocks();
	cout << endl << "  -> Total Number of Computational Cells: "
	     << QuadTree.countUsedCells();
	cout << endl << "  -> Number of Mesh Refinement Levels: "
	     << QuadTree.highestRefinementLevel()+1;
	cout << endl << "  -> Refinement Efficiency: "
	     << QuadTree.efficiencyRefinement() << endl;
      }

//       if (CFFC_Primary_MPI_Processor()) {
//          for (int j_blk = 0; j_blk < QuadTree.Nblk; j_blk++) {
//             for (int i_blk = 0; i_blk < QuadTree.Ncpu; i_blk++) {
//                if (QuadTree.Blocks[i_blk][j_blk] != NULL) {
//                   cout << endl << " cpu = " 
//                        << i_blk
//                        << " blk = "
//                        << j_blk
//                        << " blk = "
//                        << QuadTree.Blocks[i_blk][j_blk]->block;
//                } else {
//                   cout << endl << " cpu = " 
//                        << i_blk
//                        << " blk = "
//                        << j_blk;
//                }
//             }
//          }
//       }

      // Apply boundary conditions.
      BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

    } else if (command_flag == WRITE_OUTPUT_CODE) {
      // Output solution data.
      if (!batch_flag) cout << endl << " Writing Electrostatic2D solution to output data file(s).";
      error_flag = Output_Tecplot(Local_SolnBlk,
				  List_of_Local_Solution_Blocks,
				  Input_Parameters,
				  number_of_time_steps,
				  Time);
      if (error_flag) {
	cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return NULL;

    } else if (command_flag == WRITE_OUTPUT_CELLS_CODE) {
      // Output solution data.
      if (!batch_flag) cout << endl << " Writing cell-centered Electrostatic2D solution to output data file(s).";
      error_flag = Output_Cells_Tecplot(Local_SolnBlk,
					List_of_Local_Solution_Blocks,
					Input_Parameters,
					number_of_time_steps,
					Time);
      if (error_flag) {
	cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D cell output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return NULL;
      
    } else if (command_flag == WRITE_OUTPUT_QUASI3D_CODE) {
      // Output solution data.
      if (!batch_flag) cout << endl << " Writing Electrostatic2D quasi3D solution to output data file(s).";
      error_flag = Output_Quasi3D_Tecplot(Local_SolnBlk,
					  List_of_Local_Solution_Blocks,
					  Input_Parameters,
					  number_of_time_steps,
					  Time);
      if (error_flag) {
	cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D quasi3D output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return NULL;

    } else if (command_flag == WRITE_RESTART_CODE) {
      // Write restart files.
      if (!batch_flag) cout << endl << " Writing Electrostatic2D solution to restart data file(s).";
      // Write the quadtree restart file.
      error_flag = Write_QuadTree(QuadTree,Input_Parameters);
      if (error_flag) {
	cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D quadtree data file "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return NULL;
      // Write the solution block restart files.
      error_flag = Write_Restart_Solution(Local_SolnBlk,
					  List_of_Local_Solution_Blocks,
					  Input_Parameters,
					  number_of_time_steps,
					  Time,
					  processor_cpu_time);
      if (error_flag) {
	cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D restart output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return NULL;

    } else if (command_flag == WRITE_OUTPUT_GRID_CODE) {
      // Output multi-block solution-adaptive mesh data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << endl << " Writing Electrostatic2D multi-block mesh to grid data output file.";
	error_flag = Output_Tecplot(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D mesh data output file." << endl;
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return NULL;
      
    } else if (command_flag == WRITE_GRID_DEFINITION_CODE) {
      // Write multi-block solution-adaptive mesh definition files.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << endl << " Writing Electrostatic2D multi-block mesh to grid definition files.";
	error_flag = Write_Multi_Block_Grid_Definition(MeshBlk,
						       Input_Parameters);
	error_flag = Write_Multi_Block_Grid(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D multi-block mesh definition files." << endl;
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return NULL;
      
    } else if (command_flag == WRITE_OUTPUT_GRID_NODES_CODE) {
      // Output multi-block solution-adaptive mesh node data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << endl << " Writing Electrostatic2D multi-block mesh to node data output file.";
	error_flag = Output_Nodes_Tecplot(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D mesh node data output file." << endl;
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return NULL;
      
    } else if (command_flag == WRITE_OUTPUT_GRID_CELLS_CODE) {
      // Output multi-block solution-adaptive mesh cell data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << endl << " Writing Electrostatic2D multi-block mesh to cell data output file.";
	error_flag = Output_Cells_Tecplot(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << endl << " Electrostatic2D ERROR: Unable to open Electrostatic2D mesh cell data output file." << endl;
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return NULL;

    } else if (command_flag == INVALID_INPUT_CODE ||
	       command_flag == INVALID_INPUT_VALUE) {
      line_number = -line_number;
      cout << endl << " Electrostatic2D ERROR: Error reading Electrostatic2D data at line #"
	   << -line_number  << " of input data file." << endl;
      return NULL;
    }
    
  }
  
  /********************************************************************
   * End of all Electrostatic2DSolver computations and I/O.           *
   ********************************************************************/
  return Local_SolnBlk;
  
}
