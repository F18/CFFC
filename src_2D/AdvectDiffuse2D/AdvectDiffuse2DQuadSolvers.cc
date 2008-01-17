/*! \file AdvectDiffuse2DQuadSolvers.cc
  @brief 2D Advection Diffusion Equation Multi-Block Quadrilateral Mesh Solvers. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "AdvectDiffuse2DQuad.h" /* Include 2D advection diffusion equation quadrilateral mesh solution header file. */
#include "../Grid/HO_Grid2DQuadMultiBlock.h" /* Include 2D quadrilateral multiblock grid header file */
#include "AdvectDiffuse2DQuadFASMultigrid.h" /* Include the multigrid header file. */
#include "AdvectDiffuse2DQuadNKS.h"          /* Include 2D Newton-Krylov-Schwarz solver header file for advection-diffusion. */
#include "../HighOrderReconstruction/AccuracyAssessment2DMultiBlock.h" /* Include 2D accuracy assessment for multi-block level. */

/******************************************************//**
 * Routine: AdvectDiffuse2DQuadSolver                   
 *                                                      
 * Computes solutions to 2D advection diffusion         
 * equations on 2D quadrilateral multi-block            
 * solution-adaptive mesh.                              
 *                                                      
 ********************************************************/
int AdvectDiffuse2DQuadSolver(char *Input_File_Name_ptr,
			      int batch_flag) {
  
  /********************************************************  
   * Local variable declarations.                         *
   ********************************************************/
  
  // AdvectDiffuse2D input variables and parameters:
  AdvectDiffuse2D_Input_Parameters Input_Parameters;

  /* Multi-block solution-adaptive quadrilateral mesh 
     solution variables. */

#ifdef USE_HIGH_ORDER_GRID 
  Grid2D_Quad_MultiBlock_HO     MeshBlk;
#else
  Grid2D_Quad_Block           **MeshBlk;
#endif

  QuadTreeBlock_DataStructure   QuadTree;
  AdaptiveBlockResourceList     List_of_Global_Solution_Blocks;
  AdaptiveBlock2D_List          List_of_Local_Solution_Blocks;
  AdvectDiffuse2D_Quad_Block   *Local_SolnBlk;

  FAS_Multigrid2D_Solver<AdvectDiffuse2D_State, 
                         AdvectDiffuse2D_Quad_Block, 
                         AdvectDiffuse2D_Input_Parameters> MGSolver;

  /* Define residual file and cpu time variables. */

  ofstream residual_file;
  CPUTime processor_cpu_time, total_cpu_time, NKS_total_cpu_time;
  time_t start_explicit = 0, end_explicit = 0;

  /* Other local solution variables. */

  int number_of_time_steps, first_step,
    command_flag, error_flag, line_number, 
    i_stage, perform_explicit_time_marching, limiter_freezing_off;

  double Time, dTime;

  double residual_l1_norm, residual_l2_norm, residual_max_norm;
  double residual_l2norm_first, residual_ratio;

  /********************************************************  
   * Set default values for the input solution parameters *
   * and then read user specified input values from the   *
   * specified input parameter file.                      *
   ********************************************************/

  // The primary MPI processor reads and parses the input file.
  if (CFFC_Primary_MPI_Processor()) {
    if (!batch_flag) {
      cout << "\n Reading AdvectDiffuse2D input data file `"
	   << Input_File_Name_ptr << "'.";
      cout.flush();
    } /* endif */
    error_flag = Process_Input_Control_Parameter_File(Input_Parameters,
						      Input_File_Name_ptr,
						      command_flag);
    if (!batch_flag && error_flag == 0) {
      cout << Input_Parameters << "\n";
      cout.flush();
    } /* endif */
  } else {
    error_flag = 0;
  } /* endif */

  // Broadcast input solution parameters to other MPI processors.
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
  CFFC_Broadcast_MPI(&error_flag, 1);
  if (error_flag != 0) return (error_flag);
  CFFC_Broadcast_MPI(&command_flag, 1);
  if (command_flag == TERMINATE_CODE) return (0);
  Broadcast_Input_Parameters(Input_Parameters);
  Input_Parameters.Verbose(batch_flag);    //< Set Input_Parameters to batch_mode if required

  /*************************************************************  
   * Create initial mesh and allocate AdvectDiffuse2D solution *
   * variables for specified IBVP/BVP problem.                 *
   *************************************************************/

 execute_new_calculation: ;
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

  /* Create initial mesh.  Read mesh from grid definition or data files 
     when specified by input parameters. */

#ifdef USE_HIGH_ORDER_GRID
  // The primary MPI processor creates the mesh.
  if (CFFC_Primary_MPI_Processor()) {
    if (!batch_flag) cout << "\n Creating (or reading) initial quadrilateral multi-block mesh.";
    error_flag = MeshBlk.Multi_Block_Grid(Input_Parameters);

    if (error_flag) {
      cout << "\n AdvectDiffuse2D ERROR: Unable to create valid AdvectDiffuse2D multi-block mesh.\n";
      cout.flush();
    } /* endif */
  }

  // Broadcast the mesh to other MPI processors.
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
  CFFC_Broadcast_MPI(&error_flag, 1); // Broadcast mesh error flag.
  if (error_flag) return (error_flag);
  MeshBlk.Broadcast_Multi_Block_Grid();

  /* Create (allocate) multi-block quadtree data structure, create
     (allocate) array of local 2D advection diffusion equation solution blocks, 
     assign and create (allocate) 2D advection diffusion equation solution blocks
     corresponding to the initial mesh. */

  if (!batch_flag) cout << "\n Creating multi-block quadtree data structure and assigning"
                        << "\n  AdvectDiffuse2D solution blocks corresponding to initial mesh.";
  Local_SolnBlk = Create_Initial_Solution_Blocks(MeshBlk.Grid_ptr,
                                                 Local_SolnBlk,
                                                 Input_Parameters,
                                                 QuadTree,
                                                 List_of_Global_Solution_Blocks,
                                                 List_of_Local_Solution_Blocks);
  if (Local_SolnBlk == NULL) return (1);

#else

  // The primary MPI processor creates the mesh.
  if (CFFC_Primary_MPI_Processor()) {
    if (!batch_flag) cout << "\n Creating (or reading) initial quadrilateral multi-block mesh.";
    MeshBlk = NULL;
    MeshBlk = Multi_Block_Grid(MeshBlk, 
			       Input_Parameters);

    if (MeshBlk == NULL) {
      error_flag = 1;
    } else if (Check_Multi_Block_Grid(MeshBlk,
				      Input_Parameters.Number_of_Blocks_Idir,
				      Input_Parameters.Number_of_Blocks_Jdir)) {
      error_flag = 1;
    } else {
      error_flag = 0;
    } /* endif */

    if (error_flag) {
      cout << "\n AdvectDiffuse2D ERROR: Unable to create valid AdvectDiffuse2D multi-block mesh.\n";
      cout.flush();
    } /* endif */
  } else {
    MeshBlk = NULL;
  } /* endif */

  // Broadcast the mesh to other MPI processors.
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
  CFFC_Broadcast_MPI(&error_flag, 1); // Broadcast mesh error flag.
  if (error_flag) return (error_flag);
  MeshBlk = Broadcast_Multi_Block_Grid(MeshBlk, 
                                       Input_Parameters);

  /* Create (allocate) multi-block quadtree data structure, create
     (allocate) array of local 2D advection diffusion equation solution blocks, 
     assign and create (allocate) 2D advection diffusion equation solution blocks
     corresponding to the initial mesh. */

  if (!batch_flag) cout << "\n Creating multi-block quadtree data structure and assigning"
                        << "\n  AdvectDiffuse2D solution blocks corresponding to initial mesh.";
  Local_SolnBlk = Create_Initial_Solution_Blocks(MeshBlk,
                                                 Local_SolnBlk,
                                                 Input_Parameters,
                                                 QuadTree,
                                                 List_of_Global_Solution_Blocks,
                                                 List_of_Local_Solution_Blocks);
  if (Local_SolnBlk == NULL) return (1);
#endif

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

  /********************************************************  
   * Initialize AdvectDiffuse2D solution variables.       *
   ********************************************************/

  /* Set the initial time level. */

  Time = ZERO;
  number_of_time_steps = 0;

  /* Set the CPU time to zero. */

  processor_cpu_time.zero();
  total_cpu_time.zero();

  /* Initialize the state solution variables. */
  if (!batch_flag) cout << "\n Prescribing AdvectDiffuse2D initial data.";
  if (Input_Parameters.i_ICs == IC_RESTART) {
    if (!batch_flag) cout << "\n Reading AdvectDiffuse2D solution from restart data files.";

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

    error_flag = Read_QuadTree(QuadTree,
			       List_of_Global_Solution_Blocks,
			       List_of_Local_Solution_Blocks, 
			       Input_Parameters);
    if (error_flag) {
      cout << "\n AdvectDiffuse2D ERROR: Unable to open AdvectDiffuse2D quadtree data file "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    Allocate_Message_Buffers(List_of_Local_Solution_Blocks,
			     Local_SolnBlk[0].NumVar()+NUM_COMP_VECTOR2D);

    error_flag = Read_Restart_Solution(Local_SolnBlk, 
				       List_of_Local_Solution_Blocks, 
				       Input_Parameters,
				       number_of_time_steps,
				       Time,
				       processor_cpu_time);
    if (error_flag) {
      cout << "\n AdvectDiffuse2D ERROR: Unable to open AdvectDiffuse2D restart input data file(s) "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } /* endif */

    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);

    // Ensure each processor has the correct time and time!!!
    number_of_time_steps = CFFC_Maximum_MPI(number_of_time_steps);
    Time = CFFC_Maximum_MPI(Time);
    processor_cpu_time.cput = CFFC_Maximum_MPI(processor_cpu_time.cput);
    Input_Parameters.Maximum_Number_of_Time_Steps =
      CFFC_Maximum_MPI(Input_Parameters.Maximum_Number_of_Time_Steps);
  } else {
    ICs(Local_SolnBlk, 
	List_of_Local_Solution_Blocks, 
	Input_Parameters);
  } /* endif */

  /* Send solution information between neighbouring blocks to complete
     prescription of initial data. */

  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

  error_flag = Send_All_Messages(Local_SolnBlk, 
                                 List_of_Local_Solution_Blocks,
                                 NUM_COMP_VECTOR2D,
                                 ON);

  if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk, 
                                                  List_of_Local_Solution_Blocks,
                                                  NUM_VAR_ADVECTDIFFUSE2D,
                                                  OFF);

  if (error_flag) {
    cout << "\n AdvectDiffuse2D ERROR: Message passing error during AdvectDiffuse2D solution intialization "
	 << "on processor "
	 << List_of_Local_Solution_Blocks.ThisCPU
	 << ".\n";
    cout.flush();
  } /* endif */

  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return (error_flag);

  /* Prescribe boundary data consistent with initial data. */

  BCs(Local_SolnBlk, 
      List_of_Local_Solution_Blocks,
      Input_Parameters);

  /* Perform initial mesh refinement. */

  if (Input_Parameters.i_ICs != IC_RESTART) {
    if (!batch_flag) cout << "\n Performing AdvectDiffuse2D uniform mesh refinement.";
    error_flag = Uniform_AMR(Local_SolnBlk,
			     Input_Parameters,
			     QuadTree,
			     List_of_Global_Solution_Blocks,
			     List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n AdvectDiffuse2D ERROR: Uniform AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
      cout.flush();
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    if (!batch_flag) cout << "\n Performing AdvectDiffuse2D boundary mesh refinement.";
    error_flag = Boundary_AMR(Local_SolnBlk,
			      Input_Parameters,
			      QuadTree,
			      List_of_Global_Solution_Blocks,
			      List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n AdvectDiffuse2D ERROR: Boundary AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
      cout.flush();
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;
     
    if (!batch_flag) cout << "\n Performing AdvectDiffuse2D initial mesh refinement.";
    error_flag = Initial_AMR(Local_SolnBlk,
			     Input_Parameters,
			     QuadTree,
			     List_of_Global_Solution_Blocks,
			     List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n AdvectDiffuse2D ERROR: Initial AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
  } /* endif */

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
     cout << "\n  -> Number of Mesh Refinement Levels: " 
	  << QuadTree.highestRefinementLevel()+1;
     cout << "\n  -> Refinement Efficiency: " 
          << QuadTree.efficiencyRefinement() << "\n"; 
     cout.flush();
  } /* endif */

  /********************************************************  
   * Solve IBVP or BVP for 2D advection diffusion         *
   * equations on multi-block solution-adaptive           *
   * quadrilateral mesh.                                  *
   ********************************************************/

  continue_existing_calculation: ;
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

  // Reset accuracy assessment
  AccuracyAssessment2D_MultiBlock::ResetForNewCalculation(Local_SolnBlk,
							  List_of_Local_Solution_Blocks);

  if(!batch_flag) { time(&start_explicit); }

  if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {

    // Allocate memory for multigrid solver.
    error_flag = MGSolver.allocate(Local_SolnBlk,
				   &QuadTree,
				   &List_of_Global_Solution_Blocks,
				   &List_of_Local_Solution_Blocks,
				   &Input_Parameters);
    if (error_flag) {
      cout << "\n AdvectDiffuse2D ERROR: Unable to allocate memory for multigrid solver.\n";
      cout.flush();
    }
    CFFC_Broadcast_MPI(&error_flag,1);
    if (error_flag) return error_flag;

    // Execute multigrid solver.
    error_flag = MGSolver.Execute(batch_flag,
				  number_of_time_steps,
				  Time,
				  processor_cpu_time,
				  total_cpu_time,
				  residual_file);
    if (error_flag) {
      cout << "\n AdvectDiffuse2D ERROR: Error during multigrid solution.\n";
      cout.flush();
    }
    CFFC_Broadcast_MPI(&error_flag,1);
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
    
    /* Open residual file and reset the CPU time. */
    
    first_step = 1;
    limiter_freezing_off = ON;
    
    if (CFFC_Primary_MPI_Processor()) {
      error_flag = Open_Progress_File(residual_file,
				      Input_Parameters.Output_File_Name,
				      number_of_time_steps);
      if (error_flag) {
        cout << "\n AdvectDiffuse2D ERROR: Unable to open residual file for AdvectDiffuse2D calculation.\n";
        cout.flush();
      } /* endif */
    } /* endif */
    
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
    CFFC_Broadcast_MPI(&error_flag, 1);
    if (error_flag) return (error_flag);

    processor_cpu_time.reset();
    
    /* Perform required number of iterations (time steps). */
    
    if ((!Input_Parameters.Time_Accurate &&
	 Input_Parameters.Maximum_Number_of_Time_Steps > 0 &&
	 number_of_time_steps < Input_Parameters.Maximum_Number_of_Time_Steps) ||
        (!Input_Parameters.Time_Accurate &&
	 Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0) ||
	(Input_Parameters.Time_Accurate &&
	 Input_Parameters.Time_Max > Time)) {
      if (!batch_flag){ cout << "\n\n Beginning AdvectDiffuse2D computations on "
			     << Date_And_Time() << ".\n\n"; start_explicit = clock(); }

      if ((!Input_Parameters.Time_Accurate &&
	   Input_Parameters.Maximum_Number_of_Time_Steps > 0 &&
	   number_of_time_steps < Input_Parameters.Maximum_Number_of_Time_Steps) ||
  	  (Input_Parameters.Time_Accurate &&
	   Input_Parameters.Time_Max > Time)) {
	 perform_explicit_time_marching = ON;
      } else {
	 perform_explicit_time_marching = OFF;
      } /* endif */

      while (perform_explicit_time_marching) {
	/* Periodically refine the mesh (AMR). */
        if (Input_Parameters.AMR) {
           if (!first_step &&
	       number_of_time_steps-Input_Parameters.AMR_Frequency*
	       (number_of_time_steps/Input_Parameters.AMR_Frequency) == 0 ) {
              if (!batch_flag) cout << "\n\n Refining Grid.  Performing adaptive mesh refinement at n = "
                                    << number_of_time_steps << ".";
              Evaluate_Limiters(Local_SolnBlk, 
                                List_of_Local_Solution_Blocks);
              error_flag = AMR(Local_SolnBlk,
			       Input_Parameters,
                               QuadTree,
                               List_of_Global_Solution_Blocks,
                               List_of_Local_Solution_Blocks,
                               ON,ON);
              if (error_flag) {
                 cout << "\n AdvectDiffuse2D ERROR: AdvectDiffuse2D AMR error on processor "
	              << List_of_Local_Solution_Blocks.ThisCPU
	              << ".\n";
	         cout.flush();
              } /* endif */
              error_flag = CFFC_OR_MPI(error_flag);
              if (error_flag) {
                 command_flag = Output_Tecplot(Local_SolnBlk,
                                               List_of_Local_Solution_Blocks,
                                               Input_Parameters,
                                               number_of_time_steps,
                                               Time);
                 return (error_flag);
              } /* endif */
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
                      << QuadTree.efficiencyRefinement() << "\n";
                 cout.flush();
              } /* endif */
           } /* endif */
        } /* endif */

	/* Determine local and global time steps. */
	dTime = CFL(Local_SolnBlk, 
		    List_of_Local_Solution_Blocks,
                    Input_Parameters);
	dTime = CFFC_Minimum_MPI(dTime); // Find global minimum time step for all processors.

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
	  } /* endif */
	} /* endif */
	if (!Input_Parameters.Local_Time_Stepping) { // Set global time step.
	  Set_Global_TimeStep(Local_SolnBlk, 
			      List_of_Local_Solution_Blocks,
			      dTime);
	} /* endif */
	
	/* Determine the L1, L2, and max norms of the solution residual. */
	residual_l1_norm = L1_Norm_Residual(Local_SolnBlk, List_of_Local_Solution_Blocks);
	residual_l1_norm = CFFC_Summation_MPI(residual_l1_norm); // L1 norm for all processors.
	
	residual_l2_norm = L2_Norm_Residual(Local_SolnBlk, List_of_Local_Solution_Blocks);
	residual_l2_norm = sqr(residual_l2_norm);
	residual_l2_norm = CFFC_Summation_MPI(residual_l2_norm); // L2 norm for all processors.
	residual_l2_norm = sqrt(residual_l2_norm);
	
	residual_max_norm = Max_Norm_Residual(Local_SolnBlk, List_of_Local_Solution_Blocks);
	residual_max_norm = CFFC_Maximum_MPI(residual_max_norm); // Max norm for all processors.
	
	/* Update CPU time used for the calculation so far. */
	processor_cpu_time.update();
	total_cpu_time.cput = 
	  CFFC_Summation_MPI(processor_cpu_time.cput); // Total CPU time for all processors.
	
	/* Periodically save restart solution files. */
	if (!first_step &&
	    number_of_time_steps-Input_Parameters.Restart_Solution_Save_Frequency*
	    (number_of_time_steps/Input_Parameters.Restart_Solution_Save_Frequency) == 0 ) {

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

	  if (!batch_flag) cout << "\n\n  Saving AdvectDiffuse2D solution to restart data file(s) after"
				<< " n = " << number_of_time_steps << " steps (iterations).";
	  error_flag = Write_Restart_Solution(Local_SolnBlk, 
					      List_of_Local_Solution_Blocks, 
					      Input_Parameters,
					      number_of_time_steps,
					      Time,
					      processor_cpu_time);
	  if (error_flag) {
	    cout << "\n AdvectDiffuse2D ERROR: Unable to open AdvectDiffuse2D restart output data file(s) "
		 << "on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  } /* endif */

	  if (CFFC_Primary_MPI_Processor()) {
	    if (!batch_flag) cout << "\n  Saving QuadTree data file.";
	    error_flag = Write_QuadTree(QuadTree, Input_Parameters);
	    if (error_flag) {
	      cout << "\n AdvectDiffuse2D ERROR: Unable to open QuadTree data file(s).\n";
	      cout.flush();
	    } /* endif */
	  } /* endif */

	  CFFC_Broadcast_MPI(&error_flag, 1);
	  if (error_flag) return (error_flag);

	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return (error_flag);
	  cout << "\n";
	  cout.flush();
	} /* endif */
	
	/* Output progress information for the calculation. */
	 if (!batch_flag) Output_Progress_L2norm(number_of_time_steps,
						 Time*THOUSAND,
						 total_cpu_time,
						 residual_l2_norm,
						 first_step,
						 50);

	 //  	if (!batch_flag) Output_Progress(number_of_time_steps,
	 //  					 Time*THOUSAND,
	 //  					 total_cpu_time,
	 //  					 residual_l1_norm,
	 //  					 first_step,
	 //  					 50);

	 if (CFFC_Primary_MPI_Processor() && !first_step) {
	   Output_Progress_to_File(residual_file,
				   number_of_time_steps,
				   Time*THOUSAND,
				   total_cpu_time,
				   residual_l1_norm,
				   residual_l2_norm,
				   residual_max_norm);
	 } /* endif */
	
	/* Check to see if calculations are complete. */
	if (!Input_Parameters.Time_Accurate &&
	    number_of_time_steps >= 
	    Input_Parameters.Maximum_Number_of_Time_Steps) break;
	if (Input_Parameters.Time_Accurate &&
	    Time >= Input_Parameters.Time_Max) break;
	
        /* Freeze limiters as necessary. */

        if (!first_step &&
            Input_Parameters.Freeze_Limiter &&
            limiter_freezing_off &&           
            residual_l2_norm <= Input_Parameters.Freeze_Limiter_Residual_Level) {
  	   Freeze_Limiters(Local_SolnBlk, 
			   List_of_Local_Solution_Blocks);
	   limiter_freezing_off = ON;	
        } /* endif */

	/* Update solution for next time step using a multistage
	   time stepping scheme. */
	for ( i_stage  = 1 ; i_stage <= Input_Parameters.N_Stage ; ++i_stage ) {
	  // 1. Exchange solution information between neighbouring blocks.
	  error_flag = Send_All_Messages(Local_SolnBlk, 
					 List_of_Local_Solution_Blocks,
					 NUM_VAR_ADVECTDIFFUSE2D,
					 OFF);
	  if (error_flag) {
	    cout << "\n AdvectDiffuse2D ERROR: AdvectDiffuse2D message passing error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  } /* endif */
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return (error_flag);
	  
	  // 2. Apply boundary conditions for stage.
	  BCs(Local_SolnBlk, 
	      List_of_Local_Solution_Blocks,
	      Input_Parameters);
	  
	  // 3. Determine solution residuals for stage.
	  error_flag = dUdt_Multistage_Explicit(Local_SolnBlk, 
						List_of_Global_Solution_Blocks,
						List_of_Local_Solution_Blocks,
						Input_Parameters,
						i_stage);
	  if (error_flag) {
	    cout << "\n AdvectDiffuse2D ERROR: AdvectDiffuse2D solution residual error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  } /* endif */
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return (error_flag);
	  
	  // 4. Send boundary flux corrections at block interfaces with resolution changes.
	  error_flag = Send_Conservative_Flux_Corrections(Local_SolnBlk, 
							  List_of_Local_Solution_Blocks,
							  NUM_VAR_ADVECTDIFFUSE2D);
	  if (error_flag) {
	    cout << "\n AdvectDiffuse2D ERROR: AdvectDiffuse2D flux correction message passing error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  } /* endif */
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return (error_flag);
	  
	  // 5. Apply boundary flux corrections to ensure that method is conservative.
	  Apply_Boundary_Flux_Corrections_Multistage_Explicit(Local_SolnBlk, 
							      List_of_Local_Solution_Blocks,
							      Input_Parameters,
							      i_stage);
	  
          // 6. Smooth the solution residual using implicit residual smoothing. */
          if (Input_Parameters.Residual_Smoothing) {
             Residual_Smoothing(Local_SolnBlk,
                                List_of_Local_Solution_Blocks,
                                Input_Parameters,
                                i_stage);
          } /* endif */

	  // 7. Update solution for stage.
	  error_flag = Update_Solution_Multistage_Explicit(Local_SolnBlk, 
							   List_of_Local_Solution_Blocks,
							   Input_Parameters,
							   i_stage);
	  if (error_flag) {
	    cout << "\n AdvectDiffuse2D ERROR: AdvectDiffuse2D solution update error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  } /* endif */
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return (error_flag);
	  
	} /* endfor */
	
	/* Update time and time step counter. */
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
	} /* endif */
	
      } /* endwhile */

      if (!batch_flag){ cout << "\n\n AdvectDiffuse2D computations complete on " 
			     << Date_And_Time() << ".\n"; end_explicit = clock(); }
    } /* endif */
    
    /* Update ghostcell information and prescribe boundary conditions to ensure
       that the solution is consistent on each block. */
    
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
    
    error_flag = Send_All_Messages(Local_SolnBlk, 
				   List_of_Local_Solution_Blocks,
				   NUM_VAR_ADVECTDIFFUSE2D,
				   OFF);
    if (error_flag) {
      cout << "\n AdvectDiffuse2D ERROR: AdvectDiffuse2D message passing error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    
    BCs(Local_SolnBlk, 
	List_of_Local_Solution_Blocks,
	Input_Parameters);
    
    /* Close residual file. */
    
    if (CFFC_Primary_MPI_Processor()) error_flag = Close_Progress_File(residual_file);

  } /* endif - Multigrid or Not */
  
  if(!batch_flag) { time(&end_explicit); }


  /*************************************************************************************************************************/
  /************************ APPLY Newton_Krylov_Schwarz ********************************************************************/
  /*************************************************************************************************************************/\
  if (Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0) {
     time_t start_NKS, end_NKS;

     if (CFFC_Primary_MPI_Processor()) {
        error_flag = Open_Progress_File(residual_file,
	 			        Input_Parameters.Output_File_Name,
				        number_of_time_steps);
        if (error_flag) {
           cout << "\n AdvectDiffuse2D ERROR: Unable to open residual file for AdvectDiffuse2D calculation.\n";
           cout.flush();
        } /* endif */ 
     } /* endif */

     CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
     CFFC_Broadcast_MPI(&error_flag, 1);
     if (error_flag) return (error_flag);

     //Turn Limiter Freezing OFF for startup 
     Evaluate_Limiters(Local_SolnBlk, List_of_Local_Solution_Blocks);

     if (!batch_flag){ cout << "\n\n Beginning AdvectDiffuse2D NKS computations on " 
			    << Date_And_Time() << ".\n\n"; time(&start_NKS); }

     //Store Explicit times for output
     CPUTime Explicit_processor_cpu_time = processor_cpu_time;
     CPUTime Explicit_total_cpu_time =  total_cpu_time;
    
     //Perform NKS Iterations 
     error_flag = Newton_Krylov_Schwarz_Solver<AdvectDiffuse2D_State,                    //pass in Time for DTS ??
                                               AdvectDiffuse2D_Quad_Block,                                               
                                               AdvectDiffuse2D_Input_Parameters>(processor_cpu_time,
										 residual_file,
										 number_of_time_steps, // explicit time steps
										 Time,
										 Local_SolnBlk, 
										 List_of_Local_Solution_Blocks,
										 Input_Parameters);
      
     processor_cpu_time.update();
     total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput); 

     if (error_flag) {
        if (CFFC_Primary_MPI_Processor()) { 
   	   cout << "\n AdvectDiffuse2D_NKS ERROR: AdvectDiffuse2D solution error on processor " 
                << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
   	   cout.flush();
   	} /* endif */
     } /* endif */
   
     CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
     CFFC_Broadcast_MPI(&error_flag, 1);
     if (error_flag) return (error_flag);

     /***********************************************************************/
     if (!batch_flag) { cout << "\n\n AdvectDiffuse2D NKS computations complete on " 
			     << Date_And_Time() << ".\n"; time(&end_NKS); }

     if (!batch_flag) {
       cout<<"\n ----------------------------------------------------------------";
       cout<<"\n -------- Solution Computations Summary in minutes --------------";
       cout<<"\n ----------------------------------------------------------------";
       cout<<"\n Total Startup CPU Time\t\t = "<<Explicit_total_cpu_time.min();
       cout<<"\n Total NKS CPU Time \t\t = "<<total_cpu_time.min()-Explicit_total_cpu_time.min();
       cout<<"\n Total CPU Time \t\t = "<<total_cpu_time.min(); 
       cout<<"\n Total Startup Clock Time\t = "<<difftime(end_explicit,start_explicit)/60.0;
       cout<<"\n Total NKS Clock Time\t\t = "<<difftime(end_NKS,start_NKS)/60.0;
       cout<<"\n Total Clock Time\t\t = "<<difftime(end_NKS,start_explicit)/60.0;   //if no explicit start_explicit not defined
       cout<<"\n ----------------------------------------------------------------";
       cout<<"\n ----------------------------------------------------------------";
       cout<<"\n ----------------------------------------------------------------\n";
     } 
     //Also want to output total GMRES & NKS Iterations, and maybe max memory usage possibly??

     if (CFFC_Primary_MPI_Processor()) error_flag = Close_Progress_File(residual_file); 
     
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
    
    // \todo Add high-order reconstruction
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


  /********************************************************
   * Solution calculations complete.                      *
   * Write 2D advection diffusion solution to output and  *
   * restart files as required, reset solution            *
   * parameters, and run other cases as specified by      *
   * input parameters.                                    *
   ********************************************************/
 
  postprocess_current_calculation: ;
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

  while (1) {
     if (CFFC_Primary_MPI_Processor()) {
        Get_Next_Input_Control_Parameter(Input_Parameters);
        command_flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
        line_number = Input_Parameters.Line_Number;
     } /* endif */
     CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
     Broadcast_Input_Parameters(Input_Parameters);
     CFFC_Broadcast_MPI(&command_flag, 1);

     if (command_flag == EXECUTE_CODE) {
         // Deallocate memory for 2D advection diffusion equation solution.
         if (!batch_flag) cout << "\n Deallocating AdvectDiffuse2D solution variables.";
	 if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {
 	   MGSolver.deallocate();
	 }
         Local_SolnBlk = Deallocate(Local_SolnBlk, 
                                    Input_Parameters);
         List_of_Local_Solution_Blocks.deallocate();
         List_of_Global_Solution_Blocks.deallocate();
         QuadTree.deallocate();
#ifndef USE_HIGH_ORDER_GRID
         MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
                                               Input_Parameters.Number_of_Blocks_Idir, 
                                               Input_Parameters.Number_of_Blocks_Jdir);
#endif
         // Output input parameters for new caluculation.
         if (!batch_flag)  {
            cout << "\n\n Starting a new calculation.";
            cout << Input_Parameters << "\n";
            cout.flush();
         } /* endif */
         // Execute new calculation.
         goto execute_new_calculation;

     } else if (command_flag == TERMINATE_CODE) {
         // Deallocate memory for 2D advection diffusion equation solution.
         if (!batch_flag) cout << "\n Deallocating AdvectDiffuse2D solution variables.";
	 if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {
 	   MGSolver.deallocate();
	 }
         Local_SolnBlk = Deallocate(Local_SolnBlk, 
                                    Input_Parameters);
         List_of_Local_Solution_Blocks.deallocate();
         List_of_Global_Solution_Blocks.deallocate();
         QuadTree.deallocate();
#ifndef USE_HIGH_ORDER_GRID
         MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
                                               Input_Parameters.Number_of_Blocks_Idir, 
                                               Input_Parameters.Number_of_Blocks_Jdir);
#endif
         // Close input data file.
         if (!batch_flag) cout << "\n\n Closing AdvectDiffuse2D input data file.";
         if (CFFC_Primary_MPI_Processor()) Close_Input_File(Input_Parameters);
         // Terminate calculation.
         return (0);

     } else if (command_flag == CONTINUE_CODE) {
         // Reset maximum time step counter.
         Input_Parameters.Maximum_Number_of_Time_Steps += number_of_time_steps;
         // Output input parameters for new caluculation.
         if (!batch_flag)  {
            cout << "\n\n Continuing existing calculation.";
            cout << Input_Parameters << "\n";
            cout.flush();
         } /* endif */
         // Continue existing calculation.
         goto continue_existing_calculation;

     } else if (command_flag == REFINE_GRID_CODE) {
         // Refine mesh using block based adaptive mesh refinement algorithm.
         if (!batch_flag) cout << "\n\n Refining Grid.  Performing adaptive mesh refinement.";
         Evaluate_Limiters(Local_SolnBlk, 
                           List_of_Local_Solution_Blocks);
         error_flag = AMR(Local_SolnBlk,
			  Input_Parameters,
                          QuadTree,
                          List_of_Global_Solution_Blocks,
                          List_of_Local_Solution_Blocks,
                          ON,ON);
         if (error_flag) {
            cout << "\n AdvectDiffuse2D ERROR: AdvectDiffuse2D AMR error on processor "
                 << List_of_Local_Solution_Blocks.ThisCPU
                 << ".\n";
            cout.flush();
         } /* endif */
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);
         // Output multi-block solution-adaptive quadrilateral mesh statistics.
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
                 << QuadTree.efficiencyRefinement() << "\n";
            cout.flush();
         } /* endif */
	 //          if (CFFC_Primary_MPI_Processor()) {
	 //             for ( int j_blk = 0 ; j_blk <= QuadTree.Nblk-1 ; ++j_blk ) {
	 //                for ( int i_blk = 0 ; i_blk <= QuadTree.Ncpu-1 ; ++i_blk ) {
	 // 	          if (QuadTree.Blocks[i_blk][j_blk] != NULL) {
	 //                      cout << "\n cpu = " 
	 //                           << i_blk
	 //                           << " blk = "
	 //                           << j_blk
	 //                           << " blk = "
	 //                           << QuadTree.Blocks[i_blk][j_blk]->block;
	 //                   } else {
	 //                      cout << "\n cpu = " 
	 //                           << i_blk
	 //                           << " blk = "
	 //                           << j_blk;
	 //                   } /* endif */
	 //                } /* endfor */
	 //             } /* endfor */
	 //          } /* endif */

     } else if (command_flag == WRITE_OUTPUT_CODE) {
         // Output solution data.
         if (!batch_flag) cout << "\n Writing AdvectDiffuse2D solution to output data file(s).";
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
            cout << "\n AdvectDiffuse2D ERROR: Unable to open AdvectDiffuse2D output data file(s) "
                 << "on processor "
                 << List_of_Local_Solution_Blocks.ThisCPU
                 << ".\n";
            cout.flush();
         } /* endif */
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);

     } else if (command_flag == WRITE_OUTPUT_CELLS_CODE) {
         // Output solution data.
         if (!batch_flag) cout << "\n Writing cell-centered AdvectDiffuse2D solution to output data file(s).";
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
            cout << "\n AdvectDiffuse2D ERROR: Unable to open AdvectDiffuse2D cell output data file(s) "
                 << "on processor "
                 << List_of_Local_Solution_Blocks.ThisCPU
                 << ".\n";
            cout.flush();
         } /* endif */
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);

     } else if (command_flag == WRITE_OUTPUT_NODES_CODE) {
       if (!batch_flag) cout << "\n Writing AdvectDiffuse2D node locations to output data file(s).";
       if (Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0 ||
	   !(Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
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
	 cout << "\n AdvectDiffuse2D ERROR: Unable to open AdvectDiffuse2D nodes output data file(s) "
	      << "on processor "
	      << List_of_Local_Solution_Blocks.ThisCPU
	      << "." << endl;
       } /* endif */
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

       if (!batch_flag) cout << "\n Writing AdvectDiffuse2D solution to restart data file(s).";
       error_flag = Write_Restart_Solution(Local_SolnBlk, 
					   List_of_Local_Solution_Blocks, 
					   Input_Parameters,
					   number_of_time_steps,
					   Time,
					   processor_cpu_time);
       if (error_flag) {
	 cout << "\n AdvectDiffuse2D ERROR: Unable to open AdvectDiffuse2D restart output data file(s) "
	      << "on processor "
	      << List_of_Local_Solution_Blocks.ThisCPU
	      << ".\n";
	 cout.flush();
       } /* endif */
       error_flag = CFFC_OR_MPI(error_flag);
       if (error_flag) return (error_flag);
       error_flag = Write_QuadTree(QuadTree,
				   Input_Parameters);
       if (error_flag) {
	 cout << "\n AdvectDiffuse2D ERROR: Unable to open AdvectDiffuse2D quadtree data file "
	      << "on processor "
	      << List_of_Local_Solution_Blocks.ThisCPU
	      << ".\n";
	 cout.flush();
       } /* endif */
       error_flag = CFFC_OR_MPI(error_flag);
       if (error_flag) return (error_flag);
       if (CFFC_Primary_MPI_Processor()) {
	 System::Remove_Restart_Flag();  //Remove flag to indicate the restart is finished
       }

     } else if (command_flag == WRITE_OUTPUT_GRID_CODE) {
         // Output multi-block solution-adaptive mesh data file.
         if (CFFC_Primary_MPI_Processor()) {
            if (!batch_flag) cout << "\n Writing AdvectDiffuse2D multi-block mesh to grid data output file.";
#ifdef USE_HIGH_ORDER_GRID
	    //	    error_flag = MeshBlk.Output_Tecplot();
#else
            error_flag = Output_Tecplot(MeshBlk,
                                        Input_Parameters);
#endif
            if (error_flag) {
               cout << "\n AdvectDiffuse2D ERROR: Unable to open AdvectDiffuse2D mesh data output file.\n";
               cout.flush();
            } /* endif */
         } /* endif */
         CFFC_Broadcast_MPI(&error_flag, 1);
         if (error_flag) return (error_flag);

     } else if (command_flag == WRITE_GRID_DEFINITION_CODE) {
         // Write multi-block solution-adaptive mesh definition files.
         if (CFFC_Primary_MPI_Processor()) {
            if (!batch_flag) cout << "\n Writing AdvectDiffuse2D multi-block mesh to grid definition files.";
#ifdef USE_HIGH_ORDER_GRID

#else
            error_flag = Write_Multi_Block_Grid_Definition(MeshBlk,
                                                           Input_Parameters);
            error_flag = Write_Multi_Block_Grid(MeshBlk,
                                                Input_Parameters);
#endif
            if (error_flag) {
               cout << "\n AdvectDiffuse2D ERROR: Unable to open AdvectDiffuse2D multi-block mesh definition files.\n";
               cout.flush();
            } /* endif */
         } /* endif */
         CFFC_Broadcast_MPI(&error_flag, 1);
         if (error_flag) return (error_flag);

     } else if (command_flag == WRITE_OUTPUT_GRID_NODES_CODE) {
         // Output multi-block solution-adaptive mesh node data file.
         if (CFFC_Primary_MPI_Processor()) {
            if (!batch_flag) cout << "\n Writing AdvectDiffuse2D multi-block mesh to node data output file.";
#ifdef USE_HIGH_ORDER_GRID
            error_flag = MeshBlk.Output_Nodes_Tecplot_Using_IP(Input_Parameters);
#else
            error_flag = Output_Nodes_Tecplot(MeshBlk,
                                              Input_Parameters);
#endif
            if (error_flag) {
               cout << "\n AdvectDiffuse2D ERROR: Unable to open AdvectDiffuse2D mesh node data output file.\n";
               cout.flush();
            } /* endif */
         } /* endif */
         CFFC_Broadcast_MPI(&error_flag, 1);
         if (error_flag) return (error_flag);

     } else if (command_flag == WRITE_OUTPUT_GRID_CELLS_CODE) {
         // Output multi-block solution-adaptive mesh cell data file.
         if (CFFC_Primary_MPI_Processor()) {
            if (!batch_flag) cout << "\n Writing AdvectDiffuse2D multi-block mesh to cell data output file.";
#ifdef USE_HIGH_ORDER_GRID
            error_flag = MeshBlk.Output_Cells_Tecplot_Using_IP(Input_Parameters);
#else
            error_flag = Output_Cells_Tecplot(MeshBlk,
                                              Input_Parameters);
#endif
            if (error_flag) {
               cout << "\n AdvectDiffuse2D ERROR: Unable to open AdvectDiffuse2D mesh cell data output file.\n";
               cout.flush();
            } /* endif */
         } /* endif */
         CFFC_Broadcast_MPI(&error_flag, 1);
         if (error_flag) return (error_flag);

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
	 cout << "\n AdvectDiffuse2D ERROR: Unable to write AdvectDiffuse2D error norms data.\n"; cout.flush();
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
	 cout << "\n AdvectDiffuse2D ERROR: Unable to write AdvectDiffuse2D error norms data.\n"; cout.flush();
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
	 cout << "\n AdvectDiffuse2D ERROR: Unable to write AdvectDiffuse2D error norms data.\n"; cout.flush();
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
         cout << "\n AdvectDiffuse2D ERROR: Error reading AdvectDiffuse2D data at line #"
	      << -line_number  << " of input data file.\n";
         cout.flush();
         return (line_number);
     } /* endif */

  } /* endwhile */


  /**********************************************************  
   * End of all AdvectDiffuse2DSolver computations and I/O. *
   **********************************************************/

  return (0);
  
}
