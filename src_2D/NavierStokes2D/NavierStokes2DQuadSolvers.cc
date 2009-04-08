/*!\file NavierStokes2DQuadSolvers.cc:
  \brief 2D Navier-Stokes equation multi-block quadrilateral mesh solvers. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "NavierStokes2DQuad.h"    // Include 2D Navier-Stokes quadrilateral mesh solution header file.
#include "../FASMultigrid2D/FASMultigrid2D.h" // Include the multigrid header file.
#include "HO_NavierStokes2DQuadGrid.h" /* Include 2D quadrilateral multiblock grid header file for Navier-Stokes */
#include "NavierStokes2DQuadMultigrid.h" // Include 2D Navier-Stokes multigrid specializations header file.
#include "../NewtonKrylovSchwarz2D/NKS2D.h"
#include "NavierStokes2DQuadNKS.h"
#include "../CFD/EllipticOperatorAnalysis2D.h"  // Include the elliptic operator analysis header file.
#include "NavierStokes2DAccuracyAssessmentMultiBlock.h" /* Include 2D accuracy assessment for multi-block level. */
#include "../HighOrderReconstruction/HighOrder2D_MultiBlock.h" /* Include 2D high-order header file for multi-block level. */


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
  Grid2D_Quad_MultiBlock_HO       MeshBlk;
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
  CPUTime NKS_total_cpu_time;
  time_t start_explicit = 0, end_explicit = 0;
  time_t start_NKS = 0, end_NKS = 0;

  // Other local solution variables.
  int number_of_time_steps, first_step,
      command_flag, error_flag, line_number,
      perform_explicit_time_marching, limiter_freezing_flag;
  int start_number_of_time_steps;
  double Time, dTime;
  double residual_l2norm_first, residual_ratio;
  double residual_l1_norm, residual_l2_norm, residual_max_norm;
  bool mgsolver_is_allocated = false;

  /********************************************************************  
   * Set default values for the input solution parameters and then    *
   * read user specified input values from the specified input        *
   * parameter file.                                                  *
   ********************************************************************/

  // The primary MPI processor processes the input parameter file.
  if (CFFC_Primary_MPI_Processor()) {
    if (!batch_flag)
      cout << "\n Reading NavierStokes2D input data file `"
	   << Input_File_Name_ptr << "'.";
    error_flag = Process_Input_Control_Parameter_File(Input_Parameters,
						      Input_File_Name_ptr,
						      command_flag);
    if (!batch_flag && !error_flag && command_flag != EXECUTE_ZERO_STEPS_CODE) {
      cout << Input_Parameters << endl;
      cout.flush();
    } else if (error_flag) {
      cout << endl << "NavierStokes2D ERROR: During processing of input parameters." << endl << endl;
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
  Input_Parameters.Verbose(batch_flag);    //< Set Input_Parameters to batch_mode if required

  /********************************************************************
   * Create initial mesh and allocate NavierStokes2D solution         *
   * variables for the specified IBVP/BVP problem.                    *
   ********************************************************************/

 execute_new_calculation: ;

  // Synchronize processors.
  CFFC_Barrier_MPI();

  // Create initial mesh.  Read mesh from grid definition or data  
  // files when specified by input parameters.

  if (Input_Parameters.i_ICs != IC_RESTART && Input_Parameters.i_ICs != IC_GIVEN_STARTUP) {
    // Generate the mesh only if the current run is NOT a restart!

    // The primary MPI processor creates the initial mesh.
    if (CFFC_Primary_MPI_Processor()) {
      if (!batch_flag){ 
	cout << "\n Creating (or reading) initial quadrilateral multi-block mesh.";
	cout.flush();
      }
      error_flag = MeshBlk.Multi_Block_Grid(Input_Parameters);

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
    }

    // Synchronize processors.
    CFFC_Barrier_MPI();

    // Broadcast the mesh to other MPI processors.
    CFFC_Broadcast_MPI(&error_flag,1);
    if (error_flag) return error_flag;
    MeshBlk.Broadcast_Multi_Block_Grid();

    /* Set the number of blocks in I-dir and J-dir in Input_Parameters based
       on what resulted after the mesh has been created. */
    Input_Parameters.Number_of_Blocks_Idir = MeshBlk.Blocks_Idir();
    Input_Parameters.Number_of_Blocks_Jdir = MeshBlk.Blocks_Jdir();


    // Create (allocate) multi-block quadtree data structure, create
    // (allocate) array of local 2D NavierStokes equation solution blocks, 
    // assign and create (allocate) 2D NavierStokes equation solution blocks
    // corresponding to the initial mesh.
    if (!batch_flag) cout << "\n Creating multi-block quadtree data structure and assigning"
			  << "\n NavierStokes2D solution blocks corresponding to initial mesh.";
    Local_SolnBlk = CreateInitialSolutionBlocks(MeshBlk.Grid_ptr,
						Local_SolnBlk,
						Input_Parameters,
						QuadTree,
						List_of_Global_Solution_Blocks,
						List_of_Local_Solution_Blocks);

    /* Create (allocate) the high-order variables in each of the
       local 2D Euler solution blocks */
    HighOrder2D_MultiBlock::Create_Initial_HighOrder_Variables(Local_SolnBlk,
							       List_of_Local_Solution_Blocks);
  } else {
    // Allocate the minimum information related to the solution blocks. (i.e. use the default constructors)
    Local_SolnBlk = Allocate(Local_SolnBlk,Input_Parameters);
  }
  if (Local_SolnBlk == NULL) return 1;

#ifdef _NS_PARALLEL_DEBUG_
  error_flag = Local_SolnBlk[0].Open_Debug_Output_File(List_of_Local_Solution_Blocks.ThisCPU);
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;
#endif

  if (Input_Parameters.i_ICs != IC_RESTART && Input_Parameters.i_ICs != IC_GIVEN_STARTUP) {
    /* Output multi-block solution-adaptive quadrilateral mesh statistics. */

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
      cout << "\n  -> Refinement Efficiency: " 
	   << QuadTree.efficiencyRefinement() << "\n";
      cout.flush();
    } /* endif */
  }
 
  /********************************************************************
   * Initialize NavierStokes2D solution variables.                    *
   ********************************************************************/

  // Set the initial time level.
  Time = ZERO;
  number_of_time_steps = 0;

  // Set the CPU time to zero.
  processor_cpu_time.zero();
  total_cpu_time.zero();
  NKS_total_cpu_time.zero();

  // Initialize the conserved and primitive state solution variables.
  if (Input_Parameters.i_ICs == IC_RESTART  || Input_Parameters.i_ICs == IC_GIVEN_STARTUP) {
    if (!batch_flag){ cout << "\n Reading NavierStokes2D solution from restart data files."; cout.flush(); }

    //Check that restart files are probably not corrupt.
    if (CFFC_Primary_MPI_Processor()) {
      if(System::Restart_In_Progress()) {
	cout << "\n  Restart-in-progress flag detected, assuming data is corrupt."
	     << "\n  Uncompressing backups.";
	System::Uncompress_Restart();
	System::Remove_Restart_Flag();
	cout << "\n  Backup successfully uncompressed; reading.";
      }
    }
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

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
    error_flag = CFFC_OR_MPI(error_flag);
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
      cout.flush();
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Ensure that grid and solution can be used with the new scheme parameters
    if (Input_Parameters.i_ICs == IC_GIVEN_STARTUP){
      
      // Adjust the grid properties based on the new input parameters
      Grid2D_Quad_MultiBlock_HO::Adjust_Grid_To_New_InputParameters(Local_SolnBlk,
								    List_of_Local_Solution_Blocks, 
								    Input_Parameters,
								    HighOrder2D_Input::MaximumReconstructionOrder());

      /* Create (allocate) the high-order variables in each of the
	 local 2D Navier-Stokes solution blocks, if necessary. */
      HighOrder2D_MultiBlock::Create_Initial_HighOrder_Variables(Local_SolnBlk,
								 List_of_Local_Solution_Blocks);

      //!\todo Set BCs values if possible.
    }

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
    if (!batch_flag){ cout << "\n Prescribing NavierStokes2D initial data."; cout.flush(); }
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
						  NUM_VAR_NAVIERSTOKES2D,
						  OFF);
  if (error_flag) {
    cout << "\n NavierStokes2D ERROR: Message passing error during NavierStokes2D "
	 << "solution intialization on processor "
	 << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
    cout.flush();
  }
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Prescribe boundary data consistent with initial data.
  BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

  // Perform uniform, boundary, and initial mesh refinement.
  if (Input_Parameters.i_ICs != IC_RESTART && Input_Parameters.i_ICs != IC_GIVEN_STARTUP) {

    // Perform uniform mesh refinement.
    if (!batch_flag){ cout << "\n Performing NavierStokes2D uniform mesh refinement."; cout.flush(); }
    error_flag = Uniform_Adaptive_Mesh_Refinement(Local_SolnBlk,
						  Input_Parameters,
						  QuadTree,
						  List_of_Global_Solution_Blocks,
						  List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n NavierStokes2D ERROR: Uniform AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
      cout.flush();
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Perform boundary mesh refinement.
    if (!batch_flag){ cout << "\n Performing NavierStokes2D boundary mesh refinement."; cout.flush(); }
    error_flag = Boundary_Adaptive_Mesh_Refinement(Local_SolnBlk,
						   Input_Parameters,
						   QuadTree,
						   List_of_Global_Solution_Blocks,
						   List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n NavierStokes2D ERROR: Boundary AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
      cout.flush();
    }
    error_flag = CFFC_OR_MPI(error_flag);
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
    error_flag = CFFC_OR_MPI(error_flag);
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
						    NUM_VAR_NAVIERSTOKES2D,
						    OFF);
    if (error_flag) {
      cout << "\n NavierStokes2D ERROR: Message passing error during NavierStokes2D solution intialization "
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

  /********************************************************************
   * Solve IBVP or BVP for conservation form of 2D Navier-Stokes      *
   * equations on multi-block solution-adaptive quadrilateral mesh.   *
   ********************************************************************/

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
      cout <<"\n NavierStokes2D ERROR: Morton re-ordering error on processor "
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
    // Initial_AMR does not coarsen.
    start_number_of_time_steps = number_of_time_steps;
    if (!batch_flag) { 
      cout << "\n NavierStokes2DQuadSolver: EXECUTE_ZERO_STEPS_CODE: jumping to postprocess."; 
    }
    goto postprocess_current_calculation;
  }

 continue_existing_calculation: ;

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Set the same number of maximum solid body objects on all CPUs
  Spline2D_HO::Broadcast_Maximum_Number_Of_SolidBodies();

  // Reset accuracy assessment
  AccuracyAssessment2D_MultiBlock::ResetForNewCalculation(Local_SolnBlk,
							  List_of_Local_Solution_Blocks);

  time(&start_explicit);
  start_number_of_time_steps = number_of_time_steps;

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
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;
    mgsolver_is_allocated = true;

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
      cout << "\n NavierStokes2D ERROR: Unable to allocate memory for DTS multigrid"
	   << "\n solver.  Error number = " << error_flag << endl;
    }
    CFFC_Broadcast_MPI(&error_flag,1);
    if (error_flag) return error_flag;
    mgsolver_is_allocated = true;

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
    CFFC_Broadcast_MPI(&error_flag,1);
    if (error_flag) return error_flag;
    if (error_flag) {
      cout << "\n NavierStokes2D ERROR: Error during DTS multigrid solution.\n";
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
	cout << "\n NavierStokes2D ERROR: Unable to open residual file for NavierStokes2D calculation.\n";
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
	    
	    /* Update ghostcell information and prescribe boundary conditions to ensure
	       that the solution is consistent on each block. */
    
	    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

	    error_flag = Send_All_Messages(Local_SolnBlk, 
					   List_of_Local_Solution_Blocks,
					   NUM_VAR_NAVIERSTOKES2D,
					   OFF);
	    if (error_flag) {
	      cout << "\n NavierStokes2D ERROR: NavierStokes2D message passing error on processor "
		   << List_of_Local_Solution_Blocks.ThisCPU
		   << ".\n";
	      cout.flush();
	    } /* endif */
	    error_flag = CFFC_OR_MPI(error_flag);
	    if (error_flag) return (error_flag);
	      
	    BCs(Local_SolnBlk, 
		List_of_Local_Solution_Blocks,
		Input_Parameters);

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
	    error_flag = CFFC_OR_MPI(error_flag);
	    if (error_flag) {
	      command_flag = Output_Tecplot(Local_SolnBlk,
					    List_of_Local_Solution_Blocks,
					    Input_Parameters,
					    number_of_time_steps,
					    Time);
	      return error_flag;
	    }
	    // Set the same number of maximum solid body objects on all CPUs after AMR
	    Spline2D_HO::Broadcast_Maximum_Number_Of_SolidBodies();

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
            cout <<"\n NavierStokes2D ERROR: Morton re-ordering error on processor "
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

	// Update CPU time used for the calculation so far.
	processor_cpu_time.update();
	total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput);

	// Periodically save restart solution files.
	if (!first_step &&
	    number_of_time_steps-Input_Parameters.Restart_Solution_Save_Frequency*
	    (number_of_time_steps/Input_Parameters.Restart_Solution_Save_Frequency) == 0) {
	  
	  //  Save and delete old restart files in compressed archive (just in case)
	  if (CFFC_Primary_MPI_Processor()) {
	    cout << "\n  Creating compressed archive of (and deleting) old restarts.";
	    System::Compress_Restart();
	    cout << "\n  Writing new restart files.";
	    cout.flush();
	  }
	  CFFC_Barrier_MPI(); // MPI barrier so that other processors do
	                      // not start over writing restarts
	  
	  if (CFFC_Primary_MPI_Processor()) {
	    System::Set_Restart_Flag();  //Set flag to indicate a restart is being saved
	  }
	  
	  if (!batch_flag) {
	    cout << "\n\n  Saving NavierStokes2D solution to restart data file(s) after"
		 << " n = " << number_of_time_steps << " steps (iterations).";
	    cout.flush();
	  }
	  // Write the quadtree restart file.
	  error_flag = Write_QuadTree(QuadTree,Input_Parameters);
	  if (error_flag) {
	    cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D quadtree data file "
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
	    cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D restart output data file(s) "
		 << "on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  }
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return error_flag;
	  if (!batch_flag) { cout << endl; cout.flush(); }

	  if (CFFC_Primary_MPI_Processor()) {
	    System::Remove_Restart_Flag();  //Remove flag to indicate the restart is finished
	  }

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
	    cout << "\n NavierStokes2D ERROR: NavierStokes2D solution error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
	    cout.flush();
	  }
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return error_flag;

	  // Step 4. Send boundary flux corrections at block interfaces with resolution changes.
	  error_flag = Send_Conservative_Flux_Corrections(Local_SolnBlk,
							  List_of_Local_Solution_Blocks,
							  NUM_VAR_NAVIERSTOKES2D);
	  if (error_flag) {
	    cout << "\n NavierStokes2D ERROR: NavierStokes2D flux correction message passing error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
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
	    cout << "\n NavierStokes2D ERROR: NavierStokes2D solution update error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
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
    CFFC_Barrier_MPI();

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
      cout.flush();
    }
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Apply boundary conditions.
    BCs(Local_SolnBlk,List_of_Local_Solution_Blocks,Input_Parameters);

    // Close residual file.
    if (CFFC_Primary_MPI_Processor()) error_flag = Close_Progress_File(residual_file);

  }

  time(&end_explicit);

  /*************************************************************************************************************************/
  /************************ APPLY Newton_Krylov_Schwarz ********************************************************************/
  /*************************************************************************************************************************/\
  time(&start_NKS); 

  if (Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0) {

    processor_cpu_time.update();
    total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput);  
    double temp_t = total_cpu_time.cput;

    if (Input_Parameters.FlowType) {
      for (int nb = 0; nb < List_of_Local_Solution_Blocks.Nblk; nb++) {
        if (List_of_Local_Solution_Blocks.Block[nb].used == ADAPTIVEBLOCK2D_USED) {
          Local_SolnBlk[nb].allocate_face_grad_arrays();
        }
      }
    } else {

      switch (Input_Parameters.NKS_IP.Jacobian_Order) {
        case SECOND_ORDER_DIAMOND_WITH_HLLE:
          Input_Parameters.NKS_IP.Jacobian_Order = FIRST_ORDER_INVISCID_HLLE; 
          break;
        case SECOND_ORDER_DIAMOND_WITH_ROE:
          Input_Parameters.NKS_IP.Jacobian_Order = FIRST_ORDER_INVISCID_ROE;
          break;
        case SECOND_ORDER_DIAMOND_WITH_AUSMPLUSUP:
          Input_Parameters.NKS_IP.Jacobian_Order = FIRST_ORDER_INVISCID_AUSMPLUSUP;
          break;
      }

    }

    if (CFFC_Primary_MPI_Processor()) {
      error_flag = Open_Progress_File(residual_file,
          Input_Parameters.Output_File_Name,
          number_of_time_steps);
      if (error_flag) {
        cout << "\n NavierStokes2D ERROR: Unable to open residual file "
          << "for NavierStokes2D calculation.\n";
        cout.flush();
      } 
    }

    CFFC_Broadcast_MPI(&error_flag, 1);
    if (error_flag) { return error_flag; }

    // Turn off limiter freezing.
    Evaluate_Limiters(Local_SolnBlk, List_of_Local_Solution_Blocks);

    if (!batch_flag) {
      cout << "\n\n Beginning NavierStokes2D NKS computations on ";
      cout << Date_And_Time() << ".\n\n";
    }

    error_flag = Newton_Krylov_Schwarz_Solver<NavierStokes2D_pState,
               NavierStokes2D_Quad_Block,
               NavierStokes2D_Input_Parameters>(
                   processor_cpu_time,
                   residual_file,
                   number_of_time_steps,
		   Time,
                   Local_SolnBlk,
		   QuadTree,
		   List_of_Global_Solution_Blocks, 
                   List_of_Local_Solution_Blocks,
                   Input_Parameters);

    if (error_flag) {
      cout << "\n NavierStokes2D_NKS ERROR: NavierStokes2D solution error on processor " 
        << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
      cout.flush();
    }

    CFFC_Broadcast_MPI(&error_flag, 1);
    if (error_flag) { return error_flag; }

    processor_cpu_time.update();
    total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput);  
    NKS_total_cpu_time.cput += total_cpu_time.cput - temp_t;

    if (!batch_flag) { 
      cout << "\n\n NavierStokes2D NKS computations complete on " 
        << Date_And_Time() << ".\n"; 
    }
    if (CFFC_Primary_MPI_Processor()) {
      error_flag = Close_Progress_File(residual_file);
    }
  } // if (Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0) 

  time(&end_NKS); 

  if (!batch_flag) { 
    int tmpp = cout.precision(); cout.precision(5);
    cout << "\n |-----------------------------------------------------------";
    cout << "\n |  ";
    cout << "\n |   Solution Computational Timings:";
    cout << "\n |  ";
    cout << "\n |   The CPU times are summed from all the processors and";
    cout << "\n |   for cases with AMR includes time on all the meshes.";
    cout << "\n |  ";
    cout << "\n |   For cases with AMR the clock times are for the";
    cout << "\n |   current mesh only.";
    cout << "\n |  ";
    if (Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0) {
      cout << "\n |   For cases with AMR the NKS CPU time is the total time";
      cout << "\n |   spent doing NKS on all the meshes.";
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
    cout <<"\n Total CPU Time    = " << setw(8) 
         << total_cpu_time.min() << "  minutes";

    if (Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0) {
        cout <<"\n Startup Clock Time= " << setw(8) 
             << difftime(end_explicit,start_explicit)/60.0 << "  minutes";
        cout <<"\n NKS Clock Time    = " << setw(8) 
             << difftime(end_NKS,start_NKS)/60.0 << "  minutes";
    }
    cout <<"\n Total Clock Time  = " << setw(8) 
         << difftime(end_NKS,start_explicit)/60.0 << "  minutes";

    cout <<"\n ----------------------------------------------------------------";
    cout <<"\n ----------------------------------------------------------------";
    cout <<"\n ----------------------------------------------------------------\n";
    cout.precision(tmpp);
  } 

  /*************************************************************************************************************************/
  /*************************************************************************************************************************/
  /*************************************************************************************************************************/


  /***************************************************************
   * Perform solution reconstruction with the final average      *
   * states in order to use the true piecewise representation    *
   * of the solution for post-processing steps, such as solution *
   * plotting or accuracy assessment.                            *
   **************************************************************/
  if (CFFC_Primary_MPI_Processor() && (!batch_flag)) {
    std::cout << "\n\n ---------------------------------------\n"
	      << " Reconstruct final solution.\n";
  }

  if ( Input_Parameters.i_Reconstruction == RECONSTRUCTION_HIGH_ORDER){
    // Use high-order reconstruction
    HighOrder2D_MultiBlock::HighOrder_Reconstruction(Local_SolnBlk,
						     List_of_Local_Solution_Blocks,
						     Input_Parameters,
						     0);
  } else {
    // Use low-order reconstruction
    Linear_Reconstruction(Local_SolnBlk, 
			  List_of_Local_Solution_Blocks,
			  Input_Parameters);
  } // endif
  
  if (CFFC_Primary_MPI_Processor() && (!batch_flag)) {
    std::cout << " Solution reconstruction done.\n" << " ---------------------------------------\n";
  }
  
  /*************************************************************************************************************************/
  /*************************************************************************************************************************/


  /********************************************************************
   * Solution calculations complete.  Write 2D Navier-Stokes solution *
   * to output and restart files as required, reset solution          *
   * parameters, and run other cases specified by input parameters.   *
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
	cout << "\n NavierStokes2D ERROR: Error reading NavierStokes2D data at line #"
	     << -line_number << " of input data file.\n";
	return 1;
      }
      Input_Parameters.doInternalSetupAndConsistencyChecks(error_flag);
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
      // Deallocate memory for 2D NavierStokes equation solution.
      if (!batch_flag) cout << "\n Deallocating NavierStokes2D solution variables.";
      if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	  Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
			if (mgsolver_is_allocated) {
	MGSolver.deallocate();
		mgsolver_is_allocated = false;
			}
      }
      Local_SolnBlk = Deallocate(Local_SolnBlk,Input_Parameters);
      //Deallocate_Message_Buffers(List_of_Local_Solution_Blocks); // Not necessary here!
      List_of_Local_Solution_Blocks.deallocate();
      List_of_Global_Solution_Blocks.deallocate();
      QuadTree.deallocate();
      Spline2D_HO::ResetCounter(); //< reset the counter for the number of track solid bodies.
      // Output input parameters for new caluculation.
      if (!batch_flag) {
	cout << "\n\n Starting a new calculation." 
	     << Input_Parameters << "\n";
	cout.flush();
      } /* endif */
      // Execute new calculation.
      goto execute_new_calculation;

    } else if (command_flag == TERMINATE_CODE) {

      CFFC_Barrier_MPI();

      // Deallocate memory for 2D NavierStokes equation solution.
      if (!batch_flag) cout << "\n Deallocating NavierStokes2D solution variables.";
      if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	  Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
			if (mgsolver_is_allocated) {
	MGSolver.deallocate();
		mgsolver_is_allocated = false;
			}
      }
      Local_SolnBlk = Deallocate(Local_SolnBlk,Input_Parameters);
      //Deallocate_Message_Buffers(List_of_Local_Solution_Blocks); // Not necessary here!
      List_of_Local_Solution_Blocks.deallocate();
      List_of_Global_Solution_Blocks.deallocate();
      QuadTree.deallocate();
      // Close input data file.
      if (!batch_flag) { cout << "\n\n Closing NavierStokes2D input data file."; cout.flush(); }
      if (CFFC_Primary_MPI_Processor()) Close_Input_File(Input_Parameters);

#ifdef _NS_PARALLEL_DEBUG_
      Local_SolnBlk[0].Close_Debug_Output_File();
#endif

      // Terminate calculation.
      return 0;

    } else if (command_flag == CONTINUE_CODE) {
      // Reset maximum time step counter.
        Input_Parameters.Maximum_Number_of_Time_Steps += 
        number_of_time_steps - start_number_of_time_steps;
	// Output input parameters for continuing calculation.
	if (!batch_flag){ 
	  cout << "\n\n Continuing existing calculation."
	       << Input_Parameters << "\n";
	  cout.flush();
	}
      // Deallocate multigrid if necessary.
      if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	  Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
			if (mgsolver_is_allocated) {
	MGSolver.deallocate();
		mgsolver_is_allocated = false;
			}
      }
      // Continue existing calculation.
      goto continue_existing_calculation;

    } else if (command_flag == REFINE_GRID_CODE) {
      // Refine mesh using block based adaptive mesh refinement algorithm.
      if (!batch_flag){ cout << "\n Refining Grid.  Performing adaptive mesh refinement."; cout.flush(); }
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
      // Set the same number of maximum solid body objects on all CPUs after AMR
      Spline2D_HO::Broadcast_Maximum_Number_Of_SolidBodies();

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
        cout <<"\n NavierStokes2D ERROR: Morton re-ordering error on processor "
             << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
      } 
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) { return error_flag; }

      // Fix time step due to write/read of restart file in Morton code.
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
      if (Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0) {
        output_multigrid = false;
      }
      if (!batch_flag){ cout << "\n Writing NavierStokes2D solution to output data file(s)."; cout.flush(); }
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
        cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D output data file(s) "
             << "on processor "
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
      if (Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0) {
        output_multigrid = false;
      }
      if (!batch_flag){ cout << "\n Writing cell-centered NavierStokes2D solution to output data file(s)."; cout.flush();}
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
        cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D cell output data file(s) "
             << "on processor "
             << List_of_Local_Solution_Blocks.ThisCPU
             << "." << endl;
      }
      error_flag = CFFC_OR_MPI(error_flag);
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
      error_flag = CFFC_OR_MPI(error_flag);
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
      error_flag = CFFC_OR_MPI(error_flag);
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
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_RESTART_CODE) {
      // Write restart files.

      //  Save and delete old restart files in compressed archive (just in case)
      if (CFFC_Primary_MPI_Processor()) {
	cout << "\n  Creating compressed archive of (and deleting) old restarts.";
	System::Compress_Restart();
	cout << "\n  Writing new restart files.";
	cout.flush();
      }
      CFFC_Barrier_MPI(); // MPI barrier so that other processors do
                          // not start over writing restarts

      if (CFFC_Primary_MPI_Processor()) {
	System::Set_Restart_Flag();  //Set flag to indicate a restart is being saved
      }

      if (!batch_flag){ cout << "\n Writing NavierStokes2D solution to restart data file(s)."; cout.flush(); }
      // Write the quadtree restart file.
      error_flag = Write_QuadTree(QuadTree,Input_Parameters);
      if (error_flag) {
	cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D quadtree data file "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
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
	cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D restart output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
	cout.flush();
      }
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;
      if (CFFC_Primary_MPI_Processor()) {
	System::Remove_Restart_Flag();  //Remove flag to indicate the restart is finished
      }

    } else if (command_flag == WRITE_OUTPUT_GRID_CODE) {
      // Output multi-block solution-adaptive mesh data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag){ cout << "\n Writing NavierStokes2D multi-block mesh to grid data output file."; cout.flush(); }
	error_flag = MeshBlk.Output_Tecplot_Using_IP(Input_Parameters);
	if (error_flag) {
	  cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D mesh data output file.\n";
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      
    } else if (command_flag == WRITE_GRID_DEFINITION_CODE) {
      // Write multi-block solution-adaptive mesh definition files.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag){ cout << "\n Writing NavierStokes2D multi-block mesh to grid definition files."; cout.flush(); }
	error_flag = MeshBlk.Write_Multi_Block_Grid_Definition_Using_IP(Input_Parameters);
	error_flag = MeshBlk.Write_Multi_Block_Grid_Using_IP(Input_Parameters);
	if (error_flag) {
	  cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D multi-block mesh definition files.\n";
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      
    } else if (command_flag == WRITE_OUTPUT_GRID_NODES_CODE) {
      // Output multi-block solution-adaptive mesh node data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag){ cout << "\n Writing NavierStokes2D multi-block mesh to node data output file."; cout.flush(); }
	error_flag = MeshBlk.Output_Nodes_Tecplot_Using_IP(Input_Parameters);
	if (error_flag) {
	  cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D mesh node data output file.\n";
	  cout.flush();
	}
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      
    } else if (command_flag == WRITE_OUTPUT_GRID_CELLS_CODE) {
      // Output multi-block solution-adaptive mesh cell data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag){ cout << "\n Writing NavierStokes2D multi-block mesh to cell data output file."; cout.flush();}
	error_flag = MeshBlk.Output_Cells_Tecplot_Using_IP(Input_Parameters);
	if (error_flag) {
	  cout << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D mesh cell data output file.\n";
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
	cout << endl << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D Ringleb's flow output file." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
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
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_VISCOUS_PIPE_CODE) {
      if (!batch_flag) cout << endl << " Writing exact solution and error norms for the viscous pipe flow.";
      error_flag = Output_Viscous_Pipe_Tecplot(Local_SolnBlk,
					       List_of_Local_Solution_Blocks,
					       Input_Parameters);
      if (error_flag) {
	cout << endl << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D viscous pipe flow output file." << endl;
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

    } else if (command_flag == WRITE_OUTPUT_SUPERSONIC_HOT_JET_CODE) {
      if (!batch_flag) cout << endl << " Writing experimental solution data for the supersonic hot jet.";
      error_flag = Output_Supersonic_Hot_Jet_Tecplot(Local_SolnBlk,
						     List_of_Local_Solution_Blocks,
						     Input_Parameters);
      if (error_flag) {
	cout << endl << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D supersonic hot jet file." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_SUBSONIC_HOT_JET_CODE) {
      if (!batch_flag) cout << endl << " Writing experimental solution data for the subsonic hot jet.";
      error_flag = Output_Subsonic_Hot_Jet_Tecplot(Local_SolnBlk,
						   List_of_Local_Solution_Blocks,
						   Input_Parameters);
      if (error_flag) {
	cout << endl << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D subsonic hot jet file." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_FLAT_PLATE_CODE) {
      if (!batch_flag) cout << endl << " Writing exact solution and error norms for the flat plate flow (Blasius solution).";
      error_flag = Output_Flat_Plate_Tecplot(Local_SolnBlk,
					     List_of_Local_Solution_Blocks,
					     Input_Parameters);
      if (error_flag) {
	cout << endl << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D flat plate output file." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_DRIVEN_CAVITY_FLOW_CODE) {
      if (!batch_flag) cout << endl << " Writing the driven cavity flow output file.";
      error_flag = Output_Driven_Cavity_Flow_Tecplot(Local_SolnBlk,
						     List_of_Local_Solution_Blocks,
						     Input_Parameters);
      if (error_flag) {
	cout << endl << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D driven cavity flow output file." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_BACKWARD_FACING_STEP_CODE) {
      if (!batch_flag) cout << endl << " Writing the backward facing step output file.";
      error_flag = Output_Backward_Facing_Step_Tecplot(Local_SolnBlk,
						       List_of_Local_Solution_Blocks,
						       Input_Parameters);
      if (error_flag) {
	cout << endl << "\n NavierStokes2D ERROR: Unable to open NavierStokes2D backward facing step output file." << endl;
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
	cout << endl << "\n NavierStokes2D ERROR: Unable to perform the elliptic operator analysis." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      if (!batch_flag) cout << endl << "  -> Face gradient reconstruction using a directional derivative.";
      error_flag = Elliptic_Operator_Analysis_Directional_Derivative(Local_SolnBlk,
 								     List_of_Local_Solution_Blocks,
 								     Input_Parameters);
      if (error_flag) {
	cout << endl << "\n NavierStokes2D ERROR: Unable to perform the elliptic operator analysis." << endl;
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
	cout << endl << "\n NavierStokes2D ERROR: Unable to perform the elliptic operator analysis." << endl;
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
	cout << endl << "\n NavierStokes2D ERROR: Unable to perform the elliptic operator analysis." << endl;
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
	cout << endl << "\n NavierStokes2D ERROR: Unable to perform the elliptic operator analysis." << endl;
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
	cout << endl << "\n NavierStokes2D ERROR: Unable to perform the elliptic operator analysis." << endl;
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_ERROR_NORMS_TO_SCREEN) {
      if (CFFC_Primary_MPI_Processor() && (!batch_flag)) {
	cout << "\n\n ---------------------------------------\n" 
	     << " Writing error norms to screen ...\n";
	cout.flush();
      }

      error_flag = AccuracyAssessment2D_MultiBlock::PrintErrorNorms(Local_SolnBlk, 
								    List_of_Local_Solution_Blocks, 
								    Input_Parameters,
								    std::cout);
       
      if (CFFC_Primary_MPI_Processor() && error_flag) {
	cout << "\n NavierStokes2D ERROR: Unable to write NavierStokes2D error norms data.\n"; cout.flush();
      } // endif

      CFFC_Broadcast_MPI(&error_flag, 1);
      if (error_flag) return (error_flag);

      if (CFFC_Primary_MPI_Processor() && (!batch_flag)) {
	cout << "\n ---------------------------------------\n";       
	cout.flush();
      }

    } else if (command_flag == WRITE_ERROR_NORMS_TO_FILE) {
      if (CFFC_Primary_MPI_Processor() && (!batch_flag)) {
	cout << "\n\n ---------------------------------------\n" 
	     << " Writing error norms to output file ...\n";
	cout.flush();
      }

      error_flag = AccuracyAssessment2D_MultiBlock::WriteErrorNormsToOutputFile(Local_SolnBlk, 
										List_of_Local_Solution_Blocks, 
										Input_Parameters);

      if (CFFC_Primary_MPI_Processor() && error_flag) {
	cout << "\n NavierStokes2D ERROR: Unable to write NavierStokes2D error norms data.\n"; cout.flush();
      } // endif

      CFFC_Broadcast_MPI(&error_flag, 1);
      if (error_flag) return (error_flag);

      if (CFFC_Primary_MPI_Processor() && (!batch_flag)) {
	cout << "\n ---------------------------------------\n";       
	cout.flush();
      }

    } else if (command_flag == APPEND_ERROR_NORMS_TO_FILE) {
      if (CFFC_Primary_MPI_Processor() && (!batch_flag)) {
	cout << "\n\n ---------------------------------------\n" 
	     << " Appending error norms to output file ...\n";
	cout.flush();
      }

      error_flag = AccuracyAssessment2D_MultiBlock::AppendErrorNormsToOutputFile(Local_SolnBlk, 
										 List_of_Local_Solution_Blocks, 
										 Input_Parameters);

      if (CFFC_Primary_MPI_Processor() && error_flag) {
	cout << "\n NavierStokes2D ERROR: Unable to write NavierStokes2D error norms data.\n"; cout.flush();
      } // endif

      CFFC_Broadcast_MPI(&error_flag, 1);
      if (error_flag) return (error_flag);

      if (CFFC_Primary_MPI_Processor() && (!batch_flag)) {
	cout << "\n ---------------------------------------\n";       
	cout.flush();
      }

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
