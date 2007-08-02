/**********************************************************************
 * NavierStokes2DQuadSolvers.cc: 2D Navier-Stokes equation multi-     *
 *                               block quadrilateral mesh solvers.    *
 **********************************************************************/

// Include 2D Navier-Stokes quadrilateral mesh solution header file.

#ifndef _NAVIERSTOKES2D_QUAD_INCLUDED
#include "NavierStokes2DQuad.h"
#endif // _NAVIERSTOKES2D_QUAD_INCLUDED

// Include the multigrid header file.

#ifndef _FASMULTIGRID2D_INCLUDED
#include "../FASMultigrid2D/FASMultigrid2D.h"
#endif // _FASMULTIGRID2D_INCLUDED

// Include 2D Navier-Stokes multigrid specializations header file.

#ifndef _NAVIERSTOKES2D_QUAD_MULTIGRID_INCLUDED
#include "NavierStokes2DQuadMultigrid.h"
#endif // _NAVIERSTOKES2D_QUAD_MULTIGRID_INCLUDED

// Include the elliptic operator analysis header file.

#ifndef _ELLIPTIC2D_INCLUDED
#include "../CFD/EllipticOperatorAnalysis2D.h"
#endif // _ELLIPTIC2D_INCLUDED

/**********************************************************************
 * Routine: NavierStokes2DQuadSolver                                  *
 *                                                                    *
 * Computes solutions to 2D Navier-Stokes equations on 2D             *
 * quadrilateral multi-block solution-adaptive mesh.                  *
 *                                                                    *
 **********************************************************************/
int NavierStokes2DQuadSolver(char *Input_File_Name_ptr, int batch_flag) {

  /********************************************************************
   * Local variable declarations.                                     *
   ********************************************************************/

  // NavierStokes2D input variables and parameters.
  NavierStokes2D_Input_Parameters Input_Parameters;

  // Multi-block solution-adaptive quadrilateral mesh solution variables.
  Grid2D_Quad_Block             **MeshBlk;
  QuadTreeBlock_DataStructure     QuadTree;
  AdaptiveBlockResourceList       List_of_Global_Solution_Blocks;
  AdaptiveBlock2D_List            List_of_Local_Solution_Blocks;
  NavierStokes2D_Quad_Block      *Local_SolnBlk;

  // Multigrid declaration.
  FAS_Multigrid2D_Solver<NavierStokes2D_cState, 
                         NavierStokes2D_Quad_Block, 
                         NavierStokes2D_Input_Parameters> MGSolver;

  // Define residual file and cpu time variables.
  ofstream residual_file;
  CPUTime processor_cpu_time, total_cpu_time;

  // Other local solution variables.
  int number_of_time_steps, first_step,
      command_flag, error_flag, line_number,
      perform_explicit_time_marching, limiter_freezing_flag;
  double Time, dTime;
  double residual_l2norm_first, residual_ratio;
  double residual_l1_norm, residual_l2_norm, residual_max_norm;

  /********************************************************************  
   * Set default values for the input solution parameters and then    *
   * read user specified input values from the specified input        *
   * parameter file.                                                  *
   ********************************************************************/

  // The primary MPI processor processes the input parameter file.
  if (CFDkit_Primary_MPI_Processor()) {
    if (!batch_flag)
      cout << "\n Reading NavierStokes2D input data file `"
	   << Input_File_Name_ptr << "'.";
    error_flag = Process_Input_Control_Parameter_File(Input_Parameters,
						      Input_File_Name_ptr,
						      command_flag);
    if (!batch_flag && !error_flag) {
      cout << Input_Parameters << endl;
      cout.flush();
    } else if (error_flag) {
      cout << endl << "NavierStokes2D ERROR: During processing of input parameters." << endl << endl;
    }
  } else {
    error_flag = 0;
  }
  // MPI barrier to ensure processor synchronization.
  CFDkit_Barrier_MPI();

  // Broadcast input solution parameters to other MPI processors.
  CFDkit_Broadcast_MPI(&error_flag,1);
  if (error_flag) return error_flag;
  CFDkit_Broadcast_MPI(&command_flag,1);
  if (command_flag == TERMINATE_CODE) return 0;
  Broadcast_Input_Parameters(Input_Parameters);

  /********************************************************************
   * Create initial mesh and allocate NavierStokes2D solution         *
   * variables for the specified IBVP/BVP problem.                    *
   ********************************************************************/

 execute_new_calculation: ;

  // Synchronize processors.
  CFDkit_Barrier_MPI();

  // Create initial mesh.  Read mesh from grid definition or data  
  // files when specified by input parameters.

  // The primary MPI processor creates the initial mesh.
  if (CFDkit_Primary_MPI_Processor()) {

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
      cout << "\n NavierStokes2D ERROR: Unable to create valid NavierStokes2D multi-block mesh.\n";
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
  CFDkit_Barrier_MPI();

  // Broadcast the mesh to other MPI processors.
  CFDkit_Broadcast_MPI(&error_flag,1);
  if (error_flag) return error_flag;
  MeshBlk = Broadcast_Multi_Block_Grid(MeshBlk,Input_Parameters);

  // Create (allocate) multi-block quadtree data structure, create
  // (allocate) array of local 2D NavierStokes equation solution blocks, 
  // assign and create (allocate) 2D NavierStokes equation solution blocks
  // corresponding to the initial mesh.
  if (!batch_flag) cout << "\n Creating multi-block quadtree data structure and assigning"
                        << "\n NavierStokes2D solution blocks corresponding to initial mesh.";
  Local_SolnBlk = CreateInitialSolutionBlocks(MeshBlk,
					      Local_SolnBlk,
					      Input_Parameters,
					      QuadTree,
					      List_of_Global_Solution_Blocks,
					      List_of_Local_Solution_Blocks);
  if (Local_SolnBlk == NULL) return 1;

#ifdef _NS_PARALLEL_DEBUG_
  error_flag = Local_SolnBlk[0].Open_Debug_Output_File(List_of_Local_Solution_Blocks.ThisCPU);
  error_flag = CFDkit_OR_MPI(error_flag);
  if (error_flag) return error_flag;
#endif

  /********************************************************************
   * Initialize NavierStokes2D solution variables.                    *
   ********************************************************************/

  // Set the initial time level.
  Time = ZERO;
  number_of_time_steps = 0;

  // Set the CPU time to zero.
  processor_cpu_time.zero();
  total_cpu_time.zero();

  // Initialize the conserved and primitive state solution variables.
  if (!batch_flag) cout << "\n Prescribing NavierStokes2D initial data.";
  if (Input_Parameters.i_ICs == IC_RESTART) {
    if (!batch_flag) cout << "\n Reading NavierStokes2D solution from restart data files.";
    // Read the quadtree restart file.
    error_flag = Read_QuadTree(QuadTree,
			       List_of_Global_Solution_Blocks,
			       List_of_Local_Solution_Blocks,
			       Input_Parameters);
    if (error_flag) {
      cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D quadtree data file "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;
    // Allocate the message buffers.
    Allocate_Message_Buffers(List_of_Local_Solution_Blocks,
			     NUM_VAR_NAVIERSTOKES2D+NUM_COMP_VECTOR2D);
    // Read the solution block restart files.
    error_flag = Read_Restart_Solution(Local_SolnBlk,
				       List_of_Local_Solution_Blocks,
				       Input_Parameters,
				       number_of_time_steps,
				       Time,
				       processor_cpu_time);
    if (error_flag) {
      cout << endl << " NavierStokes2D ERROR: Unable to open NavierStokes2D restart "
	   << "input data file(s) on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;
    // Determine the distance to the nearest wall distance.
    error_flag = Determine_Wall_Distance(Local_SolnBlk,
					 QuadTree,
					 List_of_Local_Solution_Blocks,
					 Input_Parameters);
    if (error_flag) {
      cout << " NavierStokes2D ERROR: During wall distance calculation."
	   << "  Error #" << error_flag << "." << endl;
    }
    CFDkit_Broadcast_MPI(&error_flag,1);
    if (error_flag) return error_flag;
    // Ensure each processor has the correct time and time.
    number_of_time_steps = CFDkit_Maximum_MPI(number_of_time_steps);
    Time = CFDkit_Maximum_MPI(Time);
    processor_cpu_time.cput = CFDkit_Maximum_MPI(processor_cpu_time.cput);
    Input_Parameters.Maximum_Number_of_Time_Steps = CFDkit_Maximum_MPI(Input_Parameters.Maximum_Number_of_Time_Steps);

    // Synchronize processors.
    CFDkit_Barrier_MPI();
    // Broadcast input solution parameters to other MPI processors.
    CFDkit_Broadcast_MPI(&error_flag,1);
    if (error_flag != 0) return error_flag;
    CFDkit_Broadcast_MPI(&command_flag,1);
    if (command_flag == TERMINATE_CODE) return 0;
    Broadcast_Input_Parameters(Input_Parameters);

  } else {

    // Apply initial conditions.
    ICs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

  }

  // MPI barrier to ensure processor synchronization.
  CFDkit_Barrier_MPI();

  // Send solution information between neighbouring blocks to complete
  // prescription of initial data.
  error_flag = Send_All_Messages(Local_SolnBlk,
				 List_of_Local_Solution_Blocks,
				 NUM_COMP_VECTOR2D,
				 ON);
  if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk,
						  List_of_Local_Solution_Blocks,
						  NUM_VAR_NAVIERSTOKES2D,
						  OFF);
  if (error_flag) {
    cout << "\n NavierStokes2D ERROR: Message passing error during NavierStokes2D "
	 << "solution intialization on processor "
	 << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
  }
  error_flag = CFDkit_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Prescribe boundary data consistent with initial data.
  BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

  // Perform uniform, boundary, and initial mesh refinement.
  if (Input_Parameters.i_ICs != IC_RESTART) {

    // Perform uniform mesh refinement.
    if (!batch_flag) cout << "\n Performing NavierStokes2D uniform mesh refinement.";
    error_flag = Uniform_Adaptive_Mesh_Refinement(Local_SolnBlk,
						  Input_Parameters,
						  QuadTree,
						  List_of_Global_Solution_Blocks,
						  List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n NavierStokes2D ERROR: Uniform AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Perform boundary mesh refinement.
    if (!batch_flag) cout << "\n Performing NavierStokes2D boundary mesh refinement.";
    error_flag = Boundary_Adaptive_Mesh_Refinement(Local_SolnBlk,
						   Input_Parameters,
						   QuadTree,
						   List_of_Global_Solution_Blocks,
						   List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n NavierStokes2D ERROR: Boundary AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Perform flat-plate mesh refinement.
    if (!batch_flag) cout << "\n Performing NavierStokes2D flat-plate mesh refinement.";
    error_flag = Flat_Plate_Adaptive_Mesh_Refinement(Local_SolnBlk,
						     Input_Parameters,
						     QuadTree,
						     List_of_Global_Solution_Blocks,
						     List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n NavierStokes2D ERROR: Flat-plate AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Perform initial mesh refinement.
    if (!batch_flag) cout << "\n Performing NavierStokes2D initial mesh refinement.";
    error_flag = Initial_Adaptive_Mesh_Refinement(Local_SolnBlk,
						  Input_Parameters,
						  QuadTree,
						  List_of_Global_Solution_Blocks,
						  List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n NavierStokes2D ERROR: Initial AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // MPI barrier to ensure processor synchronization.
    CFDkit_Barrier_MPI();

    // Send solution information between neighbouring blocks to complete
    // prescription of initial data.
    error_flag = Send_All_Messages(Local_SolnBlk,
 				   List_of_Local_Solution_Blocks,
				   NUM_COMP_VECTOR2D,
				   ON);
    if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk,
						    List_of_Local_Solution_Blocks,
						    NUM_VAR_NAVIERSTOKES2D,
						    OFF);
    if (error_flag) {
      cout << "\n NavierStokes2D ERROR: Message passing error during NavierStokes2D solution intialization "
	   << "on processor " << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;

  }

  // Compute wall injection speed.
  error_flag = Wall_Injection_Speed(Local_SolnBlk,
				    List_of_Global_Solution_Blocks,
				    List_of_Local_Solution_Blocks,
				    Input_Parameters);
  if (error_flag) {
    cout << "\n NavierStokes2D ERROR: NavierStokes2D error during vw-plus calculation "
	 << "on processor " << List_of_Local_Solution_Blocks.ThisCPU << endl;
  }

  // Compute dimensionless wall distance.
  error_flag = Dimensionless_Wall_Distance(Local_SolnBlk,
					   //QuadTree,
					   List_of_Local_Solution_Blocks,
					   Input_Parameters);
  if (error_flag) {
    cout << "\n NavierStokes2D ERROR: NavierStokes2D error during y-plus calculation "
	 << "on processor " << List_of_Local_Solution_Blocks.ThisCPU << endl;
  }

  // Compute dimensionless wall injection speed.
  error_flag = Dimensionless_Wall_Injection_Speed(Local_SolnBlk,
						  List_of_Local_Solution_Blocks,
						  Input_Parameters);
  if (error_flag) {
    cout << "\n NavierStokes2D ERROR: NavierStokes2D error during vw-plus calculation "
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

//       if (CFDkit_Primary_MPI_Processor()) {
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

  /********************************************************************
   * Solve IBVP or BVP for conservation form of 2D Navier-Stokes      *
   * equations on multi-block solution-adaptive quadrilateral mesh.   *
   ********************************************************************/

 continue_existing_calculation: ;

  // MPI barrier to ensure processor synchronization.
  CFDkit_Barrier_MPI();

  if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {

    // Allocate memory for multigrid solver.
    error_flag = MGSolver.allocate(Local_SolnBlk,
				   &QuadTree,
				   &List_of_Global_Solution_Blocks,
				   &List_of_Local_Solution_Blocks,
				   &Input_Parameters);
    if (error_flag) {
      cout << "\n NavierStokes2D ERROR: Unable to allocate memory for multigrid solver.\n";
      cout.flush();
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Execute multigrid solver.
    error_flag = MGSolver.Execute(batch_flag,
				  number_of_time_steps,
				  Time,
				  processor_cpu_time,
				  total_cpu_time,
				  residual_file);
    if (error_flag) {
      cout << "\n NavierStokes2D ERROR: Multigrid error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
      cout.flush();
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;

  } else if (Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {

    // Allocate memory for the multigrid solver.
    error_flag = MGSolver.allocate(Local_SolnBlk,
				   &QuadTree,
				   &List_of_Global_Solution_Blocks,
				   &List_of_Local_Solution_Blocks,
				   &Input_Parameters);
    if (error_flag) {
      cout << "\n NavierStokes2D ERROR: Unable to allocate memory for DTS multigrid"
	   << "\n solver.  Error number = " << error_flag << endl;
    }
    CFDkit_Broadcast_MPI(&error_flag,1);
    if (error_flag) return error_flag;

    // Execute DTS FAS multigrid solver.
    error_flag = MGSolver.DTS_Multigrid_Solution(batch_flag,
						 number_of_time_steps,
						 Time,
						 processor_cpu_time,
						 total_cpu_time,
						 residual_file);
    if (error_flag) {
      cout << "\n NavierStokes2D ERROR: Error during DTS multigrid solution.  Error"
	   << "\n number = " << error_flag << endl;
    }
    CFDkit_Broadcast_MPI(&error_flag,1);
    if (error_flag) return error_flag;
    if (error_flag) {
      cout << "\n NavierStokes2D ERROR: Error during DTS multigrid solution.\n";
      cout.flush();
    }
    CFDkit_Broadcast_MPI(&error_flag,1);
    if (error_flag) return error_flag;

  } else {

    // Open residual file.
    first_step = 1;
    limiter_freezing_flag = OFF;

    if (CFDkit_Primary_MPI_Processor()) {
      error_flag = Open_Progress_File(residual_file,
				      Input_Parameters.Output_File_Name,
				      number_of_time_steps);
      if (error_flag) {
	cout << "\n NavierStokes2D ERROR: Unable to open residual file for NavierStokes2D calculation.\n";
	cout.flush();
      }
    }
    // MPI barrier to ensure processor synchronization.  
    CFDkit_Barrier_MPI();
    CFDkit_Broadcast_MPI(&error_flag,1);
    if (error_flag) return error_flag;
    // Reset the CPU time.
    processor_cpu_time.reset();

    // Perform required number of iterations (time steps).
    if ((!Input_Parameters.Time_Accurate &&
	 Input_Parameters.Maximum_Number_of_Time_Steps > 0 &&
	 number_of_time_steps < Input_Parameters.Maximum_Number_of_Time_Steps) ||
	(Input_Parameters.Time_Accurate &&
	 Input_Parameters.Time_Max > Time)) {

      if (!batch_flag) cout << "\n Beginning NavierStokes2D computations on "
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
	      cout << "\n NavierStokes2D ERROR: NavierStokes2D AMR error number "
		   << error_flag << " on processor "
		   << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
	      cout.flush();
	    }
	    error_flag = CFDkit_OR_MPI(error_flag);
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

	// Determine local and global time steps.
	dTime = CFL(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);
	// Find global minimum time step for all processors.
	dTime = CFDkit_Minimum_MPI(dTime);
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
	residual_l1_norm = CFDkit_Summation_MPI(residual_l1_norm);
      
	residual_l2_norm = L2_Norm_Residual(Local_SolnBlk,List_of_Local_Solution_Blocks);
	residual_l2_norm = sqr(residual_l2_norm);
	residual_l2_norm = CFDkit_Summation_MPI(residual_l2_norm);
	residual_l2_norm = sqrt(residual_l2_norm);
      
	residual_max_norm = Max_Norm_Residual(Local_SolnBlk,List_of_Local_Solution_Blocks);
	residual_max_norm = CFDkit_Maximum_MPI(residual_max_norm);

	// Update CPU time used for the calculation so far.
	processor_cpu_time.update();
	total_cpu_time.cput = CFDkit_Summation_MPI(processor_cpu_time.cput);

	// Periodically save restart solution files.
	if (!first_step &&
	    number_of_time_steps-Input_Parameters.Restart_Solution_Save_Frequency*
	    (number_of_time_steps/Input_Parameters.Restart_Solution_Save_Frequency) == 0) {
	  if (!batch_flag) 
	    cout << "\n\n  Saving NavierStokes2D solution to restart data file(s) after"
		 << " n = " << number_of_time_steps << " steps (iterations).";
	  // Write the quadtree restart file.
	  error_flag = Write_QuadTree(QuadTree,Input_Parameters);
	  if (error_flag) {
	    cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D quadtree data file "
		 << "on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  }
	  error_flag = CFDkit_OR_MPI(error_flag);
	  if (error_flag) return error_flag;
	  // Write the solution block restart files.
	  error_flag = Write_Restart_Solution(Local_SolnBlk,
					      List_of_Local_Solution_Blocks,
					      Input_Parameters,
					      number_of_time_steps,
					      Time,
					      processor_cpu_time);
	  if (error_flag) {
	    cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D restart output data file(s) "
		 << "on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  }
	  error_flag = CFDkit_OR_MPI(error_flag);
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
	if (CFDkit_Primary_MPI_Processor() && !first_step)
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
							  NUM_VAR_NAVIERSTOKES2D,
							  OFF);
	  if (error_flag) {
	    cout << "\n NavierStokes2D ERROR: NavierStokes2D message passing error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
	  }
	  error_flag = CFDkit_OR_MPI(error_flag);
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
	    cout << "\n NavierStokes2D ERROR: NavierStokes2D solution error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
	  }
	  error_flag = CFDkit_OR_MPI(error_flag);
	  if (error_flag) return error_flag;

	  // Step 4. Send boundary flux corrections at block interfaces with resolution changes.
	  error_flag = Send_Conservative_Flux_Corrections(Local_SolnBlk,
							  List_of_Local_Solution_Blocks,
							  NUM_VAR_NAVIERSTOKES2D);
	  if (error_flag) {
	    cout << "\n NavierStokes2D ERROR: NavierStokes2D flux correction message passing error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
	  }
	  error_flag = CFDkit_OR_MPI(error_flag);
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
	    cout << "\n NavierStokes2D ERROR: NavierStokes2D solution update error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
	  }
	  error_flag = CFDkit_OR_MPI(error_flag);
	  if (error_flag) return error_flag;

	  // Step 8. Apply turbulence boundary conditions.
	  error_flag = Turbulence_BCs(Local_SolnBlk,
				      List_of_Local_Solution_Blocks,
				      Input_Parameters);
	  error_flag = CFDkit_OR_MPI(error_flag);
	  if (error_flag) return error_flag;

	}

	// Update time and time step counter.
	if (first_step) first_step = 0;
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

      }

      if (!batch_flag) cout << "\n\n NavierStokes2D computations complete on " 
			    << Date_And_Time() << "." << endl;
    }

    // MPI barrier to ensure processor synchronization.
    CFDkit_Barrier_MPI();

    // Update ghostcell information and prescribe boundary conditions to 
    // ensure that the solution is consistent on each block.
    error_flag = Send_All_Messages(Local_SolnBlk,
				   List_of_Local_Solution_Blocks,
				   NUM_COMP_VECTOR2D,
				   ON);
    if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk,
						    List_of_Local_Solution_Blocks,
						    NUM_VAR_NAVIERSTOKES2D,
						    OFF);
    if (error_flag) {
      cout << "\n NavierStokes2D ERROR: NavierStokes2D message passing error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Apply boundary conditions.
    BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

    // Close residual file.
    if (CFDkit_Primary_MPI_Processor()) error_flag = Close_Progress_File(residual_file);

  }

  /********************************************************************
   * Solution calculations complete.  Write 2D Navier-Stokes solution *
   * to output and restart files as required, reset solution          *
   * parameters, and run other cases specified by input parameters.   *
   ********************************************************************/

 postprocess_current_calculation: ;

  // MPI barrier to ensure processor synchronization.
  CFDkit_Barrier_MPI();

  while (1) {
    if (CFDkit_Primary_MPI_Processor()) {
      Get_Next_Input_Control_Parameter(Input_Parameters);
      command_flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
      line_number = Input_Parameters.Line_Number;
      if (command_flag == INVALID_INPUT_CODE ||
	  command_flag == INVALID_INPUT_VALUE) {
	line_number = -line_number;
	cout << "\n NavierStokes2D ERROR: Error reading NavierStokes2D data at line #"
	     << -line_number << " of input data file.\n";
	return 1;
      }
      Reinitialize_Reference_State(Input_Parameters);
    }
    // MPI barrier to ensure processor synchronization.
    CFDkit_Barrier_MPI(); 
    Broadcast_Input_Parameters(Input_Parameters);
    CFDkit_Broadcast_MPI(&command_flag,1);
    for (int nb = 0; nb < List_of_Local_Solution_Blocks.Nblk; nb++) {
      if (List_of_Local_Solution_Blocks.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
	// Set flow type indicator (inviscid/viscous).
	Local_SolnBlk[nb].Flow_Type = Input_Parameters.FlowType;
	// Set flow geometry indicator (planar/axisymmetric) flow block.
	Local_SolnBlk[nb].Axisymmetric = Input_Parameters.Axisymmetric;
      }
    }

    if (command_flag == EXECUTE_CODE) {
      // Deallocate memory for 2D NavierStokes equation solution.
      if (!batch_flag) cout << "\n Deallocating NavierStokes2D solution variables.";
      if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	  Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
	MGSolver.deallocate();
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

      CFDkit_Barrier_MPI();

      // Deallocate memory for 2D NavierStokes equation solution.
      if (!batch_flag) cout << "\n Deallocating NavierStokes2D solution variables.";
      if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	  Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
	MGSolver.deallocate();
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
      if (!batch_flag) cout << "\n\n Closing NavierStokes2D input data file.";
      if (CFDkit_Primary_MPI_Processor()) Close_Input_File(Input_Parameters);

#ifdef _NS_PARALLEL_DEBUG_
      Local_SolnBlk[0].Close_Debug_Output_File();
#endif

      // Terminate calculation.
      return 0;

    } else if (command_flag == CONTINUE_CODE) {
      // Reset maximum time step counter.
      Input_Parameters.Maximum_Number_of_Time_Steps += number_of_time_steps;
      // Output input parameters for continuing calculation.
      if (!batch_flag) cout << "\n\n Continuing existing calculation."
			    << Input_Parameters << "\n";
      // Deallocate multigrid if necessary.
      if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	  Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
	MGSolver.deallocate();
      }
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
	cout << "\n NavierStokes2D ERROR: NavierStokes2D AMR error on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".  Error number = " << error_flag << "." << endl;
      }
      error_flag = CFDkit_OR_MPI(error_flag);
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

//       if (CFDkit_Primary_MPI_Processor()) {
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

      // Apply boundary conditions.
      BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

    } else if (command_flag == WRITE_OUTPUT_CODE) {
      if (!batch_flag) cout << "\n Writing NavierStokes2D solution to output data file(s).";
      if (!(Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	    Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING)) {
	error_flag = Output_Tecplot(Local_SolnBlk,
				    List_of_Local_Solution_Blocks,
				    Input_Parameters,
				    number_of_time_steps,
				    Time);
      } else {
	error_flag = MGSolver.Output_Multigrid(number_of_time_steps,
					       Time);
      }
      if (error_flag) {
	cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_CELLS_CODE) {
      if (!batch_flag) cout << "\n Writing cell-centered NavierStokes2D solution to output data file(s).";
      if (!(Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	    Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING)) {
	error_flag = Output_Cells_Tecplot(Local_SolnBlk,
					  List_of_Local_Solution_Blocks,
					  Input_Parameters,
					  number_of_time_steps,
					  Time);
      } else {
	error_flag = MGSolver.Output_Multigrid_Cells(number_of_time_steps,
						     Time);
      }
      if (error_flag) {
	cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D cell output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_NODES_CODE) {
      if (!batch_flag) cout << "\n Writing NavierStokes2D node locations to output data file(s).";
      if (!(Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	    Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING)) {
	error_flag = Output_Nodes_Tecplot(Local_SolnBlk,
					  List_of_Local_Solution_Blocks,
					  Input_Parameters,
					  number_of_time_steps,
					  Time);
      } else {
	error_flag = MGSolver.Output_Multigrid_Nodes(number_of_time_steps,
						     Time);
      }
      if (error_flag) {
	cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D nodes output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_GRADIENTS_CODE) {
      if (!batch_flag) cout << "\n Writing NavierStokes2D node locations to output data file(s).";
      error_flag = Output_Gradients_Tecplot(Local_SolnBlk,
					    List_of_Local_Solution_Blocks,
					    Input_Parameters,
					    number_of_time_steps,
					    Time);
      if (error_flag) {
	cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D nodes output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_QUASI3D_CODE) {
      // Output solution data.
      if (!batch_flag) cout << "\n Writing NavierStokes2D quasi3D solution to output data file(s).";
      error_flag = Output_Quasi3D_Tecplot(Local_SolnBlk,
					  List_of_Local_Solution_Blocks,
					  Input_Parameters,
					  number_of_time_steps,
					  Time);
      if (error_flag) {
	cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D quasi3D output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_RESTART_CODE) {
      // Write restart files.
      if (!batch_flag) cout << "\n Writing NavierStokes2D solution to restart data file(s).";
      // Write the quadtree restart file.
      error_flag = Write_QuadTree(QuadTree,Input_Parameters);
      if (error_flag) {
	cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D quadtree data file "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return error_flag;
      // Write the solution block restart files.
      error_flag = Write_Restart_Solution(Local_SolnBlk,
					  List_of_Local_Solution_Blocks,
					  Input_Parameters,
					  number_of_time_steps,
					  Time,
					  processor_cpu_time);
      if (error_flag) {
	cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D restart output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_GRID_CODE) {
      // Output multi-block solution-adaptive mesh data file.
      if (CFDkit_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing NavierStokes2D multi-block mesh to grid data output file.";
	error_flag = Output_Tecplot(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D mesh data output file.\n";
	  cout.flush();
	}
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      
    } else if (command_flag == WRITE_GRID_DEFINITION_CODE) {
      // Write multi-block solution-adaptive mesh definition files.
      if (CFDkit_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing NavierStokes2D multi-block mesh to grid definition files.";
	error_flag = Write_Multi_Block_Grid_Definition(MeshBlk,
						       Input_Parameters);
	error_flag = Write_Multi_Block_Grid(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D multi-block mesh definition files.\n";
	  cout.flush();
	}
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      
    } else if (command_flag == WRITE_OUTPUT_GRID_NODES_CODE) {
      // Output multi-block solution-adaptive mesh node data file.
      if (CFDkit_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing NavierStokes2D multi-block mesh to node data output file.";
	error_flag = Output_Nodes_Tecplot(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D mesh node data output file.\n";
	  cout.flush();
	}
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      
    } else if (command_flag == WRITE_OUTPUT_GRID_CELLS_CODE) {
      // Output multi-block solution-adaptive mesh cell data file.
      if (CFDkit_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing NavierStokes2D multi-block mesh to cell data output file.";
	error_flag = Output_Cells_Tecplot(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D mesh cell data output file.\n";
	  cout.flush();
	}
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_RINGLEB_CODE) {
      if (!batch_flag) cout << endl << " Writing exact solution and error norms for Ringleb's flow.";
      error_flag = Output_Ringleb_Tecplot(Local_SolnBlk,
					  List_of_Local_Solution_Blocks,
					  Input_Parameters);
      if (error_flag) {
	cout << endl << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D Ringleb's flow output file." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_VISCOUS_CHANNEL_CODE) {
      if (!batch_flag) cout << endl << " Writing exact solution and error norms for the viscous channel flow.";
      error_flag = Output_Viscous_Channel_Tecplot(Local_SolnBlk,
						  List_of_Local_Solution_Blocks,
						  Input_Parameters);
      if (error_flag) {
	cout << endl << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D viscous channel flow output file." << endl;
	cout.flush();
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_VISCOUS_PIPE_CODE) {
      if (!batch_flag) cout << endl << " Writing exact solution and error norms for the viscous pipe flow.";
      error_flag = Output_Viscous_Pipe_Tecplot(Local_SolnBlk,
					       List_of_Local_Solution_Blocks,
					       Input_Parameters);
      if (error_flag) {
	cout << endl << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D viscous pipe flow output file." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_TURBULENT_PIPE_CODE) {
      if (!batch_flag) cout << endl << " Writing experimental and computed solution data for the turbulent pipe flow.";
      error_flag = Output_Turbulent_Pipe_Tecplot(Local_SolnBlk,
						 List_of_Local_Solution_Blocks,
						 Input_Parameters);
      if (error_flag) {
	cout << endl << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D viscous pipe flow output file." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_FLAT_PLATE_CODE) {
      if (!batch_flag) cout << endl << " Writing exact solution and error norms for the flat plate flow (Blasius solution).";
      error_flag = Output_Flat_Plate_Tecplot(Local_SolnBlk,
					     List_of_Local_Solution_Blocks,
					     Input_Parameters);
      if (error_flag) {
	cout << endl << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D flat plate output file." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_DRIVEN_CAVITY_FLOW_CODE) {
      if (!batch_flag) cout << endl << " Writing the driven cavity flow output file.";
      error_flag = Output_Driven_Cavity_Flow_Tecplot(Local_SolnBlk,
						     List_of_Local_Solution_Blocks,
						     Input_Parameters);
      if (error_flag) {
	cout << endl << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D driven cavity flow output file." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_BACKWARD_FACING_STEP_CODE) {
      if (!batch_flag) cout << endl << " Writing the backward facing step output file.";
      error_flag = Output_Backward_Facing_Step_Tecplot(Local_SolnBlk,
						       List_of_Local_Solution_Blocks,
						       Input_Parameters);
      if (error_flag) {
	cout << endl << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D backward facing step output file." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
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
	cout << endl << "\n NavierStokes2D ERROR: Unable to perform the elliptic operator analysis." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      if (!batch_flag) cout << endl << "  -> Face gradient reconstruction using a directional derivative.";
      error_flag = Elliptic_Operator_Analysis_Directional_Derivative(Local_SolnBlk,
 								     List_of_Local_Solution_Blocks,
 								     Input_Parameters);
      if (error_flag) {
	cout << endl << "\n NavierStokes2D ERROR: Unable to perform the elliptic operator analysis." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      if (!batch_flag) cout << endl << "  -> Face gradient reconstruction using Green-Gauss reconstruction over"
			    << endl << "     the diamond-path with simple vertex weighting.";
      error_flag = Elliptic_Operator_Analysis_Diamond_Path(Local_SolnBlk,
 							   List_of_Local_Solution_Blocks,
 							   Input_Parameters,
 							   0);
      if (error_flag) {
	cout << endl << "\n NavierStokes2D ERROR: Unable to perform the elliptic operator analysis." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      if (!batch_flag) cout << endl << "  -> Face gradient reconstruction using Green-Gauss reconstruction over"
			    << endl << "     the diamond-path with the linearity-preserving weighting function"
			    << endl << "     of Holmes and Connell.";
      error_flag = Elliptic_Operator_Analysis_Diamond_Path(Local_SolnBlk,
 							   List_of_Local_Solution_Blocks,
 							   Input_Parameters,
 							   1);
      if (error_flag) {
	cout << endl << "\n NavierStokes2D ERROR: Unable to perform the elliptic operator analysis." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      if (!batch_flag) cout << endl << "  -> Face gradient reconstruction using Green-Gauss reconstruction over"
			    << endl << "     the diamond-path with the linearity-preserving weighting function"
			    << endl << "     of Zingg and Yarrow.";
      error_flag = Elliptic_Operator_Analysis_Diamond_Path(Local_SolnBlk,
 							   List_of_Local_Solution_Blocks,
 							   Input_Parameters,
 							   2);
      if (error_flag) {
	cout << endl << "\n NavierStokes2D ERROR: Unable to perform the elliptic operator analysis." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      if (!batch_flag) cout << endl << "  -> Face gradient reconstruction using the corrected average cell-"
			    << endl << "     centred gradient constructed with linear least-squares reconstruction";
      error_flag = Elliptic_Operator_Analysis_Average_Gradient_Linear_Least_Squares(Local_SolnBlk,
										    List_of_Local_Solution_Blocks,
										    Input_Parameters,
										    ON);
      if (error_flag) {
	cout << endl << "\n NavierStokes2D ERROR: Unable to perform the elliptic operator analysis." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == INVALID_INPUT_CODE ||
	       command_flag == INVALID_INPUT_VALUE) {
      line_number = -line_number;
      cout << "\n NavierStokes2D ERROR: Error reading NavierStokes2D data at line #"
	   << -line_number  << " of input data file.\n";
      return line_number;
    }
    
  }
  
  /********************************************************************
   * End of all NavierStokes2DSolver computations and I/O.            *
   ********************************************************************/
  return 0;
  
}
