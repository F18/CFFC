/**********************************************************************
 * Dusty2DQuadSolvers.cc: 2D Dusty equation multi-block quadrilateral *
 *                        mesh solvers.                               *
 **********************************************************************/

// Include 2D Dusty quadrilateral mesh solution header file.

#ifndef _DUSTY2D_QUAD_INCLUDED
#include "Dusty2DQuad.h"
#endif // _DUSTY2D_QUAD_INCLUDED

// Include 2D Dusty multigrid specializations header file.

#ifndef _DUSTY2D_QUAD_MULTIGRID_INCLUDED
#include "Dusty2DQuadMultigrid.h"
#endif // _DUSTY2D_QUAD_MULTIGRID_INCLUDED

// Include 2D Electrostatic quadrilateral mesh solution header file.

#ifndef _ELECTROSTATIC2D_QUAD_INCLUDED
#include "../Electrostatic2D/Electrostatic2DQuad.h"
#endif // _ELECTROSTATIC2D_QUAD_INCLUDED

/**********************************************************************
 * Routine: Dusty2DQuadSolver                                         *
 *                                                                    *
 * Computes solutions to 2D Dusty equations on 2D quadrilateral       *
 * multi-block solution-adaptive mesh.                                *
 *                                                                    *
 **********************************************************************/
int Dusty2DQuadSolver(char *Input_File_Name_ptr, int batch_flag) {

  /********************************************************************
   * Local variable declarations.                                     *
   ********************************************************************/

  // 2D multi-block grid blocks for mesh generation:
  Grid2D_Quad_Block  **MeshBlk;

  // Dusty2D input parameters and multi-block solution-adaptive
  // quadrilateral mesh solution variables:
  Dusty2D_Input_Parameters      Input_Parameters;
  QuadTreeBlock_DataStructure   QuadTree;
  AdaptiveBlockResourceList     List_of_Global_Solution_Blocks;
  AdaptiveBlock2D_List          List_of_Local_Solution_Blocks;
  Dusty2D_Quad_Block           *Local_SolnBlk;

  // Electrostatic2D input parameters and multi-block solution-adaptive
  // quadrilateral mesh solution variables:
  Electrostatic2D_Input_Parameters   ES_Input_Parameters;
  Electrostatic2D_Quad_Block        *ES_Local_SolnBlk;

  // Multigrid declaration.
  FAS_Multigrid2D_Solver<Dusty2D_cState, 
                         Dusty2D_Quad_Block, 
                         Dusty2D_Input_Parameters> MGSolver;

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

  ofstream conservation_file;
  Dusty2D_cState Umass;

  // Level set solution variables.
  int levelset_iterations;

  /********************************************************************  
   * Set default values for the input solution parameters and then    *
   * read user specified input values from the specified input        *
   * parameter file.                                                  *
   ********************************************************************/

  // The primary MPI processor processes the input parameter file.
  if (CFFC_Primary_MPI_Processor()) {
    if (!batch_flag)
      cout << "\n Reading Dusty2D input data file `"
	   << Input_File_Name_ptr << "'.";
    error_flag = Process_Input_Control_Parameter_File(Input_Parameters,
						      Input_File_Name_ptr,
						      command_flag);
    if (!batch_flag && !error_flag) {
      cout << Input_Parameters << endl;
      cout.flush();
    } else if (error_flag) {
      cout << "\n Dusty2D ERROR: During processing of input parameters."
	   << " " << command_flag << " " << error_flag << endl << endl;
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

  const int NUM_VAR_DUSTY2D = Input_Parameters.Wo.NUM_VAR_DUSTY2D;

  /********************************************************************
   * Create initial mesh and allocate Dusty2D solution variables for  *
   * specified IBVP/BVP problem.                                      *
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
      cout << "\n Dusty2D ERROR: Unable to create valid Dusty2D multi-block mesh.\n";
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
  // (allocate) array of local 2D Dusty equation solution blocks, 
  // assign and create (allocate) 2D Dusty equation solution blocks
  // corresponding to the initial mesh.
  if (!batch_flag) cout << "\n Creating multi-block quadtree data structure and assigning"
                        << "\n Dusty2D solution blocks corresponding to initial mesh.";
  Local_SolnBlk = CreateInitialSolutionBlocks(MeshBlk,
					      Local_SolnBlk,
					      Input_Parameters,
					      QuadTree,
					      List_of_Global_Solution_Blocks,
					      List_of_Local_Solution_Blocks);
  if (Local_SolnBlk == NULL) return 1;

  /********************************************************************
   * Initialize Dusty2D solution variables.                           *
   ********************************************************************/

  // Set the initial time level.
  Time = ZERO;
  number_of_time_steps = 0;

  // Set the CPU time to zero.
  processor_cpu_time.zero();
  total_cpu_time.zero();

  // Initialize and solve the 2D electrostatic problem if required.
  if (Input_Parameters.Electrostatic) {
    if (!batch_flag) cout << endl << endl << " Solving 2D Electrostatic potential and electric fields.";
    ES_Local_SolnBlk = Electrostatic2DQuadSolver(Input_File_Name_ptr,
						 batch_flag,
						 ES_Local_SolnBlk,
						 ES_Input_Parameters,
						 QuadTree,
						 List_of_Global_Solution_Blocks,
						 List_of_Local_Solution_Blocks);
    error_flag = (ES_Local_SolnBlk == NULL) ? 1 : 0;
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;
    error_flag = Copy_Electrostatic_Field_Variables(Local_SolnBlk,
						    ES_Local_SolnBlk,
						    List_of_Local_Solution_Blocks);
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;
    Allocate_Message_Buffers(List_of_Local_Solution_Blocks,
			     Local_SolnBlk[0].NumVar()+NUM_COMP_VECTOR2D);
  }

  // Initialize the conserved and primitive state solution variables.
  if (!batch_flag) cout << "\n Prescribing Dusty2D initial data.";
  if (Input_Parameters.i_ICs == IC_RESTART) {
    if (!batch_flag) cout << "\n Reading Dusty2D solution from restart data files.";
    // Read the quadtree restart file.
    error_flag = Read_QuadTree(QuadTree,
			       List_of_Global_Solution_Blocks,
			       List_of_Local_Solution_Blocks,
			       Input_Parameters);
    if (error_flag) {
      cout << "\n Dusty2D ERROR: Unable to open Dusty2D quadtree data file "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;
    // Allocate the message buffers.
    Allocate_Message_Buffers(List_of_Local_Solution_Blocks,
			     NUM_VAR_DUSTY2D+NUM_COMP_VECTOR2D);
    // Read the solution block restart files.
    error_flag = Read_Restart_Solution(Local_SolnBlk,
				       List_of_Local_Solution_Blocks,
				       Input_Parameters,
				       number_of_time_steps,
				       Time,
				       processor_cpu_time);
    if (error_flag) {
      cout << "\n Dusty2D ERROR: Unable to open Dusty2D restart "
	   << "input data file(s) on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
      cout.flush();
    }
    error_flag = CFFC_OR_MPI(error_flag);
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
						  NUM_VAR_DUSTY2D,
						  OFF);
  if (error_flag) {
    cout << "\n Dusty2D ERROR: Message passing error during Dusty2D "
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
    if (!batch_flag) cout << "\n Performing Dusty2D uniform mesh refinement.";
    error_flag = Uniform_Adaptive_Mesh_Refinement(Local_SolnBlk,
						  Input_Parameters,
						  QuadTree,
						  List_of_Global_Solution_Blocks,
						  List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n Dusty2D ERROR: Uniform AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
      cout.flush();
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Perform boundary mesh refinement.
    if (!batch_flag) cout << "\n Performing Dusty2D boundary mesh refinement.";
    error_flag = Boundary_Adaptive_Mesh_Refinement(Local_SolnBlk,
						   Input_Parameters,
						   QuadTree,
						   List_of_Global_Solution_Blocks,
						   List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n Dusty2D ERROR: Boundary AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
      cout.flush();
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Perform initial mesh refinement.
    if (!batch_flag) cout << "\n Performing Dusty2D initial mesh refinement.";
    error_flag = Initial_Adaptive_Mesh_Refinement(Local_SolnBlk,
						  Input_Parameters,
						  QuadTree,
						  List_of_Global_Solution_Blocks,
						  List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n Dusty2D ERROR: Initial AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
      cout.flush();
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
						    NUM_VAR_DUSTY2D,
						    OFF);
    if (error_flag) {
      cout << "\n Dusty2D ERROR: Message passing error during Dusty2D solution intialization "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
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
    cout.flush();
  }

//   if (CFFC_Primary_MPI_Processor()) {
//     for (int nblk = 0; nblk < QuadTree.Nblk; nblk++) {
//       for (int ncpu = 0; ncpu < QuadTree.Ncpu; ncpu++) {
// 	if (QuadTree.Blocks[ncpu][nblk] != NULL) {
// 	  cout << endl
// 	       << " cpu = " << ncpu << " blk = " << nblk
// 	       << " blk = " << QuadTree.Blocks[ncpu][nblk]->block;
// 	} else {
// 	  cout << endl << " cpu = " << ncpu << " blk = " << nblk;
// 	}
//       }
//     }
//   }

  /********************************************************************
   * Solve IBVP or BVP for conservation form of 2D Dusty equations on *
   * multi-block solution-adaptive quadrilateral mesh.                *
   ********************************************************************/

 continue_existing_calculation: ;

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {

    // Allocate memory for multigrid solver.
    error_flag = MGSolver.allocate(Local_SolnBlk,
				   &QuadTree,
				   &List_of_Global_Solution_Blocks,
				   &List_of_Local_Solution_Blocks,
				   &Input_Parameters);
    if (error_flag) {
      cout << "\n Dusty2D ERROR: Unable to allocate memory for multigrid solver.\n";
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
      cout << "\n Dusty2D ERROR: Multigrid error on processor "
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
    if (error_flag) {
      cout << "\n Euler2D ERROR: Unable to allocate memory for DTS multigrid solver.\n";
      cout.flush();
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
      cout << "\n Euler2D ERROR: Error during DTS multigrid solution.\n";
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
	cout << "\n Dusty2D ERROR: Unable to open residual file for Dusty2D calculation.\n";
	cout.flush();
      }
    }
    // MPI barrier to ensure processor synchronization.  
    CFFC_Barrier_MPI();
    CFFC_Broadcast_MPI(&error_flag,1);
    if (error_flag) return error_flag;
    // Reset the CPU time.
    processor_cpu_time.reset();

    // Open conservation file.
    if (CFFC_Primary_MPI_Processor() && Input_Parameters.MeasureConservation) {
      error_flag = Open_Conservation_File(conservation_file,
					  Input_Parameters.Output_File_Name,
					  number_of_time_steps,
					  Umass);
      if (error_flag) {
	cout << "\n Dusty2D ERROR: Unable to open conservation file for Dusty2D calculation.\n";
	cout.flush();
      }
    }
    CFFC_Barrier_MPI();
    CFFC_Broadcast_MPI(&error_flag,1);
    if (error_flag) return error_flag;

    // Perform required number of iterations (time steps).
    if ((!Input_Parameters.Time_Accurate &&
	 Input_Parameters.Maximum_Number_of_Time_Steps > 0 &&
	 number_of_time_steps < Input_Parameters.Maximum_Number_of_Time_Steps) ||
	(Input_Parameters.Time_Accurate &&
	 Input_Parameters.Time_Max > Time)) {

      if (!batch_flag) cout << "\n Beginning Dusty2D computations on "
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
						  List_of_Local_Solution_Blocks,
						  ON);
	    if (error_flag) {
	      cout << "\n Dusty2D ERROR: Dusty2D AMR error number "
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

	// Determine the total mass, momentum, and energy for the gas
	// and particle-phases.
	if (Input_Parameters.MeasureConservation) {
	  error_flag = Determine_Conservation_Properties(Local_SolnBlk,List_of_Local_Solution_Blocks,Umass);
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return error_flag;
	  Umass = CFFC_Summation_MPI(Umass);
	}

	// Update CPU time used for the calculation so far.
	processor_cpu_time.update();
	total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput);

	// Periodically save restart solution files.
	if (!first_step &&
	    number_of_time_steps-Input_Parameters.Restart_Solution_Save_Frequency*
	    (number_of_time_steps/Input_Parameters.Restart_Solution_Save_Frequency) == 0) {
	  if (!batch_flag) 
	    cout << "\n\n  Saving Dusty2D solution to restart data file(s) after"
		 << " n = " << number_of_time_steps << " steps (iterations).";
	  // Write the quadtree restart file.
	  error_flag = Write_QuadTree(QuadTree,Input_Parameters);
	  if (error_flag) {
	    cout << "\n Dusty2D ERROR: Unable to open Dusty2D quadtree data file "
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
	    cout << "\n Dusty2D ERROR: Unable to open Dusty2D restart output data file(s) "
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
	if (!batch_flag) Output_Progress_L2norm(number_of_time_steps,
						Time*THOUSAND,
						total_cpu_time,
						residual_l2_norm,
						first_step,
						Input_Parameters.Output_Progress_Frequency);
//	if (!batch_flag) Output_Progress(number_of_time_steps,
//					 Time*THOUSAND,
//					 total_cpu_time,
//					 residual_l1_norm,
//					 first_step,
//					 Input_Parameters.Output_Progress_Frequency);
	if (CFFC_Primary_MPI_Processor() && !first_step)
	  Output_Progress_to_File(residual_file,
				  number_of_time_steps,
				  Time*THOUSAND,
				  total_cpu_time,
				  residual_l1_norm,
				  residual_l2_norm,
				  residual_max_norm);

	if (CFFC_Primary_MPI_Processor() && !first_step && Input_Parameters.MeasureConservation)
	  Output_Conservation_to_File(conservation_file,
				      number_of_time_steps,
				      Time*THOUSAND,
				      total_cpu_time,
				      Umass);

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
							  NUM_VAR_DUSTY2D,
							  OFF);
	  if (error_flag) {
	    cout << "\n Dusty2D ERROR: Dusty2D message passing error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
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
	    cout << "\n Dusty2D ERROR: Dusty2D solution error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  }
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return error_flag;

	  // Step 4. Send boundary flux corrections at block interfaces with resolution changes.
	  error_flag = Send_Conservative_Flux_Corrections(Local_SolnBlk,
							  List_of_Local_Solution_Blocks,
							  NUM_VAR_DUSTY2D);
	  if (error_flag) {
	    cout << "\n Dusty2D ERROR: Dusty2D flux correction message passing error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
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
	    cout << "\n Dusty2D ERROR: Dusty2D solution update error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << "." << endl;
	    cout.flush();
	  }
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return error_flag;

	  // Step 8. Apply turbulence boundary conditions.
	  error_flag = Turbulence_BCs(Local_SolnBlk,
				      List_of_Local_Solution_Blocks,
				      Input_Parameters);
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return error_flag;

	}

	// Conduct particle-phase component switching if required.
	error_flag = Multi_Velocity_Component_Particle_Phase_Switch(Local_SolnBlk,Input_Parameters,
								    List_of_Local_Solution_Blocks);
	if (error_flag) {
	  cout << "\n Dusty2D ERROR: Dusty2D multicomponent particle-phase "
	       << "switching error on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU << endl;
	  cout.flush();
	}

	// Update time and time step counter.
	if (first_step) first_step = 0;
	number_of_time_steps = number_of_time_steps + 1;
	if (Input_Parameters.i_Time_Integration != 
	    TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) {
	  Time += Input_Parameters.CFL_Number*dTime;
	} else {
	  Time += Input_Parameters.CFL_Number*dTime*
	    MultiStage_Optimally_Smoothing(Input_Parameters.N_Stage,
					   Input_Parameters.N_Stage,
					   Input_Parameters.i_Limiter);
	}

      }

      if (!batch_flag) cout << "\n\n Dusty2D computations complete on " 
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
						    NUM_VAR_DUSTY2D,
						    OFF);
    if (error_flag) {
      cout << "\n Dusty2D ERROR: Dusty2D message passing error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Apply boundary conditions.
    BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

    // Close residual file.
    if (CFFC_Primary_MPI_Processor()) error_flag = Close_Progress_File(residual_file);

    // Close conservation file.
    if (CFFC_Primary_MPI_Processor() && Input_Parameters.MeasureConservation)
      error_flag = Close_Conservation_File(conservation_file);

  }

  /********************************************************************
   * Solution calculations complete.                                  *
   * Write 2D Dusty solution to output and restart files as required, *
   * reset solution parameters, and run other cases as  specified by  *
   * input parameters.                                                *
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
	cout << "\n Dusty2D ERROR: Error reading Dusty2D data at line #"
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
	// Set particle flag (on or off).
	Local_SolnBlk[nb].Particles = Input_Parameters.Particles;
	// Set flow type indicator (inviscid/viscous).
	Local_SolnBlk[nb].Flow_Type = Input_Parameters.FlowType;
	// Set flow geometry indicator (planar/axisymmetric) flow block.
	Local_SolnBlk[nb].Axisymmetric = Input_Parameters.Axisymmetric;
	// Set electrostatic flow indicator flow block.
	Local_SolnBlk[nb].Electrostatic = Input_Parameters.Electrostatic;
      }
    }

    if (command_flag == EXECUTE_CODE) {
      // Deallocate memory for 2D Dusty equation solution.
      if (!batch_flag) cout << "\n Deallocating Dusty2D solution variables.";
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

      CFFC_Barrier_MPI();

      // Deallocate memory for 2D Dusty equation solution.
      if (!batch_flag) cout << "\n Deallocating Dusty2D solution variables.";
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
      if (!batch_flag) cout << "\n\n Closing Dusty2D input data file.";
      if (CFFC_Primary_MPI_Processor()) Close_Input_File(Input_Parameters);
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
					    List_of_Local_Solution_Blocks,
					    ON);
      if (error_flag) {
	cout << "\n Dusty2D ERROR: Dusty2D AMR error on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".  Error number = " << error_flag << ".\n";
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
	cout.flush();
      }

//       if (CFFC_Primary_MPI_Processor()) {
// 	for (int nblk = 0; nblk < QuadTree.Nblk; nblk++) {
// 	  for (int ncpu = 0; ncpu < QuadTree.Ncpu; ncpu++) {
// 	    if (QuadTree.Blocks[ncpu][nblk] != NULL) {
// 	      cout << endl
// 		   << " cpu = " << ncpu << " blk = " << nblk
// 		   << " blk = " << QuadTree.Blocks[ncpu][nblk]->block;
// 	    } else {
// 	      cout << endl << " cpu = " << ncpu << " blk = " << nblk;
// 	    }
// 	  }
// 	}
//       }

      // Apply boundary conditions.
      BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

    } else if (command_flag == WRITE_OUTPUT_CODE) {
      if (!batch_flag) cout << "\n Writing Dusty2D solution to output data file(s).";
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
	cout << "\n Dusty2D ERROR: Unable to open Dusty2D output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_CELLS_CODE) {
      if (!batch_flag) cout << "\n Writing cell-centered Dusty2D solution to output data file(s).";
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
	cout << "\n Dusty2D ERROR: Unable to open Dusty2D cell output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_NODES_CODE) {
      if (!batch_flag) cout << "\n Writing Dusty2D node locations to output data file(s).";
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
	cout << "\n Dusty2D ERROR: Unable to open Dusty2D nodes output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_QUASI3D_CODE) {
      if (!batch_flag) cout << "\n Writing Dusty2D quasi3D solution to output data file(s).";
      error_flag = Output_Quasi3D_Tecplot(Local_SolnBlk,
					  List_of_Local_Solution_Blocks,
					  Input_Parameters,
					  number_of_time_steps,
					  Time);
      if (error_flag) {
	cout << "\n Dusty2D ERROR: Unable to open Dusty2D quasi3D output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_RESTART_CODE) {
      // Write restart files.
      if (!batch_flag) cout << "\n Writing Dusty2D solution to restart data file(s).";
      // Write the quadtree restart file.
      error_flag = Write_QuadTree(QuadTree,Input_Parameters);
      if (error_flag) {
	cout << "\n Dusty2D ERROR: Unable to open Dusty2D quadtree data file "
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
	cout << "\n Dusty2D ERROR: Unable to open Dusty2D restart output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_GRID_CODE) {
      // Output multi-block solution-adaptive mesh data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing Dusty2D multi-block mesh to grid data output file.";
	error_flag = Output_Tecplot(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << "\n Dusty2D ERROR: Unable to open Dusty2D mesh data output file.\n";
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      
    } else if (command_flag == WRITE_GRID_DEFINITION_CODE) {
      // Write multi-block solution-adaptive mesh definition files.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing Dusty2D multi-block mesh to grid definition files.";
	error_flag = Write_Multi_Block_Grid_Definition(MeshBlk,
						       Input_Parameters);
	error_flag = Write_Multi_Block_Grid(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << "\n Dusty2D ERROR: Unable to open Dusty2D multi-block mesh definition files.\n";
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_GRID_CELLS_CODE) {
      // Output multi-block solution-adaptive mesh cell data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing Dusty2D multi-block mesh to cell data output file.";
	error_flag = Output_Cells_Tecplot(MeshBlk,Input_Parameters);
	if (error_flag) {
	  cout << "\n Dusty2D ERROR: Unable to open Dusty2D mesh cell data output file.\n";
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_RINGLEB_CODE) {
      if (!batch_flag) cout << endl << " Writing exact solution and error norms for Ringleb's flow.";
      error_flag = Output_Ringleb(Local_SolnBlk,
				  List_of_Local_Solution_Blocks,
				  Input_Parameters);
      if (error_flag) {
	cout << endl << "\n Dusty2D ERROR: Unable to open Dusty2D Ringleb's flow output file." << endl;
	cout.flush();
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_VISCOUS_CHANNEL_CODE) {
      if (!batch_flag) cout << endl << " Writing exact solution and error norms for the viscous channel flow.";
      error_flag = Output_Viscous_Channel(Local_SolnBlk,
					  List_of_Local_Solution_Blocks,
					  Input_Parameters);
      if (error_flag) {
	cout << endl << "\n Dusty2D ERROR: Unable to open Dusty2D viscous channel flow output file." << endl;
	cout.flush();
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_VISCOUS_PIPE_CODE) {
      if (!batch_flag) cout << endl << " Writing exact solution and error norms for the viscous pipe flow.";
      error_flag = Output_Viscous_Pipe(Local_SolnBlk,
				       List_of_Local_Solution_Blocks,
				       Input_Parameters);
      if (error_flag) {
	cout << endl << "\n Dusty2D ERROR: Unable to open Dusty2D viscous pipe flow output file." << endl;
	cout.flush();
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_TURBULENT_PIPE_CODE) {
      if (!batch_flag) cout << endl << " Writing experimental and computed solution data for the turbulent pipe flow.";
      error_flag = Output_Turbulent_Pipe_Tecplot(Local_SolnBlk,
						 List_of_Local_Solution_Blocks,
						 Input_Parameters);
      if (error_flag) {
	cout << endl << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D viscous pipe flow output file." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_FLAT_PLATE_CODE) {
      if (!batch_flag) cout << endl << " Writing exact solution and error norms for the flat plate flow (Blasius solution).";
      error_flag = Output_Flat_Plate(Local_SolnBlk,
				     List_of_Local_Solution_Blocks,
				     Input_Parameters);
      if (error_flag) {
	cout << endl << "\n Dusty2D ERROR: Unable to open Dusty2D flat plate output file." << endl;
	cout.flush();
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_DRIVEN_CAVITY_FLOW_CODE) {
      if (!batch_flag) cout << endl << " Writing the driven cavity flow output file.";
      error_flag = Output_Driven_Cavity_Flow(Local_SolnBlk,
					     List_of_Local_Solution_Blocks,
					     Input_Parameters);
      if (error_flag) {
	cout << endl << "\n Dusty2D ERROR: Unable to open Dusty2D driven cavity flow." << endl;
	cout.flush();
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == INVALID_INPUT_CODE ||
	       command_flag == INVALID_INPUT_VALUE) {
      line_number = -line_number;
      cout << "\n Dusty2D ERROR: Error reading Dusty2D data at line #"
	   << -line_number  << " of input data file.\n";
      return line_number;
    }
    
  }
  
  /********************************************************************
   * End of all Dusty2DSolver computations and I/O.                   *
   ********************************************************************/
  return 0;
  
}

/**********************************************************************
 **********************************************************************/
int Open_Conservation_File(ofstream &Conservation_File,
			   char *File_Name,
			   const int Append_to_File,
			   const Dusty2D_cState &U) {

  int i;
  char prefix[256], extension[256], 
       conservation_file_name[256], gnuplot_file_name[256];
  char *conservation_file_name_ptr, *gnuplot_file_name_ptr;
  ofstream gnuplot_file;

  // Determine the name of the conservation file.
  i = 0;
  while (1) {
    if (File_Name[i] == ' ' || File_Name[i] == '.') break;
    prefix[i] = File_Name[i];
    i++;
    if (i > strlen(File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_conservation");

  strcpy(extension,".dat");
  strcpy(conservation_file_name,prefix);
  strcat(conservation_file_name,extension);

  conservation_file_name_ptr = conservation_file_name;

  // Open the conservation file.
  if (Append_to_File) {
    Conservation_File.open(conservation_file_name_ptr,ios::out|ios::app);
  } else {
    Conservation_File.open(conservation_file_name_ptr,ios::out);
  }
  if (Conservation_File.bad()) return 1;

  // Write the appropriate GNUPLOT command file for plotting
  // file information.

  strcpy(extension,".gplt");
  strcpy(gnuplot_file_name,prefix);
  strcat(gnuplot_file_name,extension);

  gnuplot_file_name_ptr = gnuplot_file_name;

  gnuplot_file.open(gnuplot_file_name_ptr,ios::out);
  if (gnuplot_file.bad()) return(1);

  gnuplot_file << "set title \"Solution Convergence\"\n"
	       << "set xlabel \"N (iterations)\"\n"
	       << "set ylabel \"conservation\"\n" 
	       << "set logscale y\n";
  gnuplot_file << "plot ";
  for (int n = 1; n <= U.NumVar(); n++) {
    gnuplot_file << "\"" << conservation_file_name_ptr << "\" using 1:2 \"%lf%*lf%*lf";
    for (int m = 1; m <= 4; m++) {
      gnuplot_file << "%";
      if (m != n) gnuplot_file << "*";
      gnuplot_file << "lf";
    }
    gnuplot_file << "\" \\\n     title \"Property " << n << "\" with lines";
    if (n == 4) gnuplot_file << "\n";
    else gnuplot_file << ", \\\n";
  }
  gnuplot_file << "pause -1  \"Hit return to continue\"\n";

  gnuplot_file.close();

  // Preparation of conservation file complete. Return zero value.
  return 0;

}

int Close_Conservation_File(ofstream &Conservation_File) {
  Conservation_File.close();
  return 0;
}

void Output_Conservation_to_File(ostream &Conservation_File,
				 const int Number_of_Time_Steps,
				 const double &Time,
				 const CPUTime &CPU_Time,
				 const Dusty2D_cState &U) {
  Conservation_File << setprecision(6);
  Conservation_File << Number_of_Time_Steps
		    << " " << Time
		    << " " << CPU_Time.min();
  Conservation_File.setf(ios::scientific);
  for (int n = 1; n <= 4; n++) Conservation_File << " " << U[n];
  Conservation_File << endl;
  Conservation_File.unsetf(ios::scientific);
  Conservation_File.flush();
}
