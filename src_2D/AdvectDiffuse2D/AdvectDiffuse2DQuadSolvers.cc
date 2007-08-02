/* AdvectDiffuse2DQuadSolvers.cc:  2D Advection Diffusion Equation 
                                   Multi-Block Quadrilateral Mesh Solvers. */

/* Include 2D advection diffusion equation quadrilateral mesh solution header file. */

#ifndef _ADVECTDIFFUSE2D_QUAD_INCLUDED
#include "AdvectDiffuse2DQuad.h"
#endif // _ADVECTDIFFUSE2D_QUAD_INCLUDED

/* Include the multigrid header file. */

#ifndef _FASMULTIGRID2D_INCLUDED
#include "../FASMultigrid2D/FASMultigrid2D.h"
#endif // _FASMULTIGRID2D_INCLUDED

/* Include 2D advection diffusion multigrid specializations header file. */

#ifndef _ADVECTDIFFUSE2D_QUAD_MULTIGRID_INCLUDED
#include "AdvectDiffuse2DQuadMultigrid.h"
#endif // _ADVECTDIFFUSE2D_QUAD_MULTIGRID_INCLUDED

/* Include 2D advection diffusion equation Newton-Krylov-Schwarz solver header file. */

#ifndef _ADVECTDIFFUSE2D_QUAD_NKS_INCLUDED
#include "AdvectDiffuse2DQuad_NKS.h"
#endif // _ADVECTDIFFUSE2D_QUAD_NKS_INCLUDED

/********************************************************
 * Routine: AdvectDiffuse2DQuadSolver                   *
 *                                                      *
 * Computes solutions to 2D advection diffusion         *
 * equations on 2D quadrilateral multi-block            *
 * solution-adaptive mesh.                              *
 *                                                      *
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
 
  Grid2D_Quad_Block           **MeshBlk;
  QuadTreeBlock_DataStructure   QuadTree;
  AdaptiveBlockResourceList     List_of_Global_Solution_Blocks;
  AdaptiveBlock2D_List          List_of_Local_Solution_Blocks;
  AdvectDiffuse2D_Quad_Block   *Local_SolnBlk;

  FAS_Multigrid2D_Solver<double, 
                         AdvectDiffuse2D_Quad_Block, 
                         AdvectDiffuse2D_Input_Parameters> MGSolver;

  /* Define residual file and cpu time variables. */

  ofstream residual_file;
  CPUTime processor_cpu_time, total_cpu_time, NKS_time;

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
  if (CFDkit_Primary_MPI_Processor()) {
     if (!batch_flag) {
        cout << "\n Reading AdvectDiffuse2D input data file `"
             << Input_File_Name_ptr << "'.";
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
  CFDkit_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
  CFDkit_Broadcast_MPI(&error_flag, 1);
  if (error_flag != 0) return (error_flag);
  CFDkit_Broadcast_MPI(&command_flag, 1);
  if (command_flag == TERMINATE_CODE) return (0);
  Broadcast_Input_Parameters(Input_Parameters);

  /*************************************************************  
   * Create initial mesh and allocate AdvectDiffuse2D solution *
   * variables for specified IBVP/BVP problem.                 *
   *************************************************************/

  execute_new_calculation: ;
  CFDkit_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

  /* Create initial mesh.  Read mesh from grid definition or data files 
     when specified by input parameters. */

  // The primary MPI processor creates the mesh.
  if (CFDkit_Primary_MPI_Processor()) {
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
  CFDkit_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
  CFDkit_Broadcast_MPI(&error_flag, 1); // Broadcast mesh error flag.
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
          << QuadTree.countUsedCells() << "\n"; 
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
  
  /* Initialize the conserved and primitive state
     solution variables. */
  if (!batch_flag) cout << "\n Prescribing AdvectDiffuse2D initial data.";
  if (Input_Parameters.i_ICs == IC_RESTART) {
     if (!batch_flag) cout << "\n Reading AdvectDiffuse2D solution from restart data files.";
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
     error_flag = CFDkit_OR_MPI(error_flag);
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
     error_flag = CFDkit_OR_MPI(error_flag);
     if (error_flag) return (error_flag);
     // Ensure each processor has the correct time and time!!!
     number_of_time_steps = CFDkit_Maximum_MPI(number_of_time_steps);
     Time = CFDkit_Maximum_MPI(Time);
     processor_cpu_time.cput = CFDkit_Maximum_MPI(processor_cpu_time.cput);
  } else {
     ICs(Local_SolnBlk, 
         List_of_Local_Solution_Blocks, 
         Input_Parameters,
	 QuadTree);
     Set_Advection_Velocity_Field(Local_SolnBlk, 
                                  List_of_Local_Solution_Blocks, 
                                  Input_Parameters);
  } /* endif */

  /* Send solution information between neighbouring blocks to complete
     prescription of initial data. */

  CFDkit_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

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
  error_flag = CFDkit_OR_MPI(error_flag);
  if (error_flag) return (error_flag);

  /* Prescribe boundary data consistent with initial data. */

  BCs(Local_SolnBlk, 
      List_of_Local_Solution_Blocks,
      Input_Parameters);

  /* Perform initial mesh refinement. */

  if (Input_Parameters.i_ICs != IC_RESTART) {
     if (!batch_flag) cout << "\n Performing AdvectDiffuse2D uniform mesh refinement.";
     error_flag = Uniform_AMR2(Local_SolnBlk,
                               Input_Parameters,
                               QuadTree,
                               List_of_Global_Solution_Blocks,
                               List_of_Local_Solution_Blocks);
     if (error_flag) {
       cout << "\n AdvectDiffuse2D ERROR: Uniform AMR error on processor "
	    << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
       cout.flush();
     } /* endif */
     error_flag = CFDkit_OR_MPI(error_flag);
     if (error_flag) return error_flag;

     if (!batch_flag) cout << "\n Performing AdvectDiffuse2D boundary mesh refinement.";
     error_flag = Boundary_AMR2(Local_SolnBlk,
                                Input_Parameters,
                                QuadTree,
                                List_of_Global_Solution_Blocks,
                                List_of_Local_Solution_Blocks);
     if (error_flag) {
       cout << "\n AdvectDiffuse2D ERROR: Boundary AMR error on processor "
	    << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
       cout.flush();
     } /* endif */
     error_flag = CFDkit_OR_MPI(error_flag);
     if (error_flag) return error_flag;

     if (!batch_flag) cout << "\n Performing AdvectDiffuse2D initial mesh refinement.";
     error_flag = Initial_AMR2(Local_SolnBlk,
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
     error_flag = CFDkit_OR_MPI(error_flag);
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

  /* Set Analytical Solution if problem is 1D-Diffusion or Circular Advection */
  Set_Analytical_Solution(Local_SolnBlk,
			  List_of_Local_Solution_Blocks,
			  Input_Parameters);

  /********************************************************  
   * Solve IBVP or BVP for 2D advection diffusion         *
   * equations on multi-block solution-adaptive           *
   * quadrilateral mesh.                                  *
   ********************************************************/

  continue_existing_calculation: ;
  CFDkit_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

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
    CFDkit_Broadcast_MPI(&error_flag,1);
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
    CFDkit_Broadcast_MPI(&error_flag,1);
    if (error_flag) return error_flag;
  
  } else {
    
    /* Open residual file and reset the CPU time. */
    
    first_step = 1;
    limiter_freezing_off = ON;
    
    if (CFDkit_Primary_MPI_Processor()) {
      error_flag = Open_Progress_File(residual_file,
				      Input_Parameters.Output_File_Name,
				      number_of_time_steps);
      if (error_flag) {
        cout << "\n AdvectDiffuse2D ERROR: Unable to open residual file for AdvectDiffuse2D calculation.\n";
        cout.flush();
      } /* endif */
    } /* endif */
    
    CFDkit_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
    CFDkit_Broadcast_MPI(&error_flag, 1);
    if (error_flag) return (error_flag);

    processor_cpu_time.reset();
    
    /* Perform required number of iterations (time steps). */
    
    if ((!Input_Parameters.Time_Accurate &&
	 Input_Parameters.Maximum_Number_of_Time_Steps > 0 &&
	 number_of_time_steps < Input_Parameters.Maximum_Number_of_Time_Steps) ||
        (!Input_Parameters.Time_Accurate &&
	 Input_Parameters.Maximum_Number_of_NKS_Iterations > 0) ||
	(Input_Parameters.Time_Accurate &&
	 Input_Parameters.Time_Max > Time)) {
      if (!batch_flag) cout << "\n\n Beginning AdvectDiffuse2D computations on "
			    << Date_And_Time() << ".\n\n";

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
              error_flag = CFDkit_OR_MPI(error_flag);
              if (error_flag) return (error_flag);
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
	dTime = CFDkit_Minimum_MPI(dTime); // Find global minimum time step for all processors.
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
	residual_l1_norm = CFDkit_Summation_MPI(residual_l1_norm); // L1 norm for all processors.
	
	residual_l2_norm = L2_Norm_Residual(Local_SolnBlk, List_of_Local_Solution_Blocks);
	residual_l2_norm = sqr(residual_l2_norm);
	residual_l2_norm = CFDkit_Summation_MPI(residual_l2_norm); // L2 norm for all processors.
	residual_l2_norm = sqrt(residual_l2_norm);
	
	residual_max_norm = Max_Norm_Residual(Local_SolnBlk, List_of_Local_Solution_Blocks);
	residual_max_norm = CFDkit_Maximum_MPI(residual_max_norm); // Max norm for all processors.
	
	/* Update CPU time used for the calculation so far. */
	processor_cpu_time.update();
	total_cpu_time.cput = 
	  CFDkit_Summation_MPI(processor_cpu_time.cput); // Total CPU time for all processors.
	
	/* Periodically save restart solution files. */
	if (!first_step &&
	    number_of_time_steps-Input_Parameters.Restart_Solution_Save_Frequency*
	    (number_of_time_steps/Input_Parameters.Restart_Solution_Save_Frequency) == 0 ) {
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
	  error_flag = CFDkit_OR_MPI(error_flag);
	  if (error_flag) return (error_flag);
	  cout << "\n";
	  cout.flush();
	} /* endif */
	
	/* Output progress information for the calculation. */
	 if (!batch_flag) Output_Progress_L2norm(number_of_time_steps,
						 Time,
						 total_cpu_time,
						 residual_l2_norm,
						 first_step,
						 50);
//  	if (!batch_flag) Output_Progress(number_of_time_steps,
//  					 Time,
//  					 total_cpu_time,
//  					 residual_l1_norm,
//  					 first_step,
//  					 50);
	if (CFDkit_Primary_MPI_Processor() && !first_step) {
	  Output_Progress_to_File(residual_file,
				  number_of_time_steps,
				  Time,
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
	  error_flag = CFDkit_OR_MPI(error_flag);
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
	  error_flag = CFDkit_OR_MPI(error_flag);
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
	  error_flag = CFDkit_OR_MPI(error_flag);
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
	  error_flag = CFDkit_OR_MPI(error_flag);
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



      if (!batch_flag) cout << "\n\n AdvectDiffuse2D computations complete on " 
			    << Date_And_Time() << ".\n";
    } /* endif */
    
    /* Update ghostcell information and prescribe boundary conditions to ensure
       that the solution is consistent on each block. */
    
    CFDkit_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
    
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
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    
    BCs(Local_SolnBlk, 
	List_of_Local_Solution_Blocks,
	Input_Parameters);
    
    /* Close residual file. */
    
    if (CFDkit_Primary_MPI_Processor()) error_flag = Close_Progress_File(residual_file);

  } /* endif - Multigrid or Not */
  
  /*******************************/
  /* APPLY Newton_Krylov_Schwarz */
  /*******************************/

  if (Input_Parameters.Maximum_Number_of_NKS_Iterations > 0) {
     if (CFDkit_Primary_MPI_Processor()) {
         error_flag = Open_Progress_File(residual_file,
	  			         Input_Parameters.Output_File_Name,
				         number_of_time_steps);
         if (error_flag) {
            cout << "\n AdvectDiffuse2D ERROR: Unable to open residual file for AdvectDiffuse2D calculation.\n";
            cout.flush();
         } /* endif */
     } /* endif */
    
     CFDkit_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
     CFDkit_Broadcast_MPI(&error_flag, 1);
     if (error_flag) return (error_flag);

     Set_Analytical_Solution(Local_SolnBlk, 
     			     List_of_Local_Solution_Blocks,
     			     Input_Parameters);
     Evaluate_Limiters(Local_SolnBlk, 
                       List_of_Local_Solution_Blocks);
     
     if (!batch_flag) cout << "\n\n Beginning AdvectDiffuse2D NKS computations on "
	 		   << Date_And_Time() << ".\n\n";

     NKS_time.zero();

     error_flag = Newton_Krylov_Schwarz_Solver(residual_file,
                                               number_of_time_steps,
     					       Local_SolnBlk, 
     					       List_of_Local_Solution_Blocks,
     					       Input_Parameters);

     NKS_time.update();
     
     if (CFDkit_Primary_MPI_Processor()) { 
        cout << "\nNKS_time in minutes = " << NKS_time.min() << endl;
     } /* endif */
     
     if (error_flag) {
        if (CFDkit_Primary_MPI_Processor()) { 
     	   cout << "\n AdvectDiffuse2D_NKS ERROR: AdvectDiffuse2D solution error on processor " 
                << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
     	   cout.flush();
     	} /* endif */
     } /* endif */

     CFDkit_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
     CFDkit_Broadcast_MPI(&error_flag, 1);
     if (error_flag) return (error_flag);

     if (!batch_flag) cout << "\n\n AdvectDiffuse2D NKS computations complete on " 
                           << Date_And_Time() << ".\n";

     if (CFDkit_Primary_MPI_Processor()) error_flag = Close_Progress_File(residual_file);
  } /* endif */

  /********************************************************
   * Solution calculations complete.                      *
   * Write 2D advection diffusion solution to output and  *
   * restart files as required, reset solution            *
   * parameters, and run other cases as specified by      *
   * input parameters.                                    *
   ********************************************************/
 
  postprocess_current_calculation: ;
  CFDkit_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

  while (1) {
     if (CFDkit_Primary_MPI_Processor()) {
        Get_Next_Input_Control_Parameter(Input_Parameters);
        command_flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
        line_number = Input_Parameters.Line_Number;
     } /* endif */
     CFDkit_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
     Broadcast_Input_Parameters(Input_Parameters);
     CFDkit_Broadcast_MPI(&command_flag, 1);

     if (command_flag == EXECUTE_CODE) {
         // Deallocate memory for 2D advection diffusion equation solution.
         if (!batch_flag) cout << "\n Deallocating AdvectDiffuse2D solution variables.";
	 if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {
 	   MGSolver.deallocate();
	 }
         Local_SolnBlk = Deallocate(Local_SolnBlk, 
                                    Input_Parameters);
         //Deallocate_Message_Buffers(List_of_Local_Solution_Blocks); // Not necessary here!
         List_of_Local_Solution_Blocks.deallocate();
         List_of_Global_Solution_Blocks.deallocate();
         QuadTree.deallocate();
         MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
                                               Input_Parameters.Number_of_Blocks_Idir, 
                                               Input_Parameters.Number_of_Blocks_Jdir);
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
         //Deallocate_Message_Buffers(List_of_Local_Solution_Blocks); // Not necessary here!
         List_of_Local_Solution_Blocks.deallocate();
         List_of_Global_Solution_Blocks.deallocate();
         QuadTree.deallocate();
         MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
                                               Input_Parameters.Number_of_Blocks_Idir, 
                                               Input_Parameters.Number_of_Blocks_Jdir);
         // Close input data file.
         if (!batch_flag) cout << "\n\n Closing AdvectDiffuse2D input data file.";
         if (CFDkit_Primary_MPI_Processor()) Close_Input_File(Input_Parameters);
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
         error_flag = CFDkit_OR_MPI(error_flag);
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
//          if (CFDkit_Primary_MPI_Processor()) {
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
         // Reset advection velocity field.
         Set_Advection_Velocity_Field(Local_SolnBlk, 
                                      List_of_Local_Solution_Blocks, 
                                      Input_Parameters);
         // Reset analytic solution if problem is 1D-Diffusion or Circular Advection.
         Set_Analytical_Solution(Local_SolnBlk,
                                 List_of_Local_Solution_Blocks,
                                 Input_Parameters);

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
         error_flag = CFDkit_OR_MPI(error_flag);
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
         error_flag = CFDkit_OR_MPI(error_flag);
         if (error_flag) return (error_flag);

     } else if (command_flag == WRITE_OUTPUT_NODES_CODE) {
       if (!batch_flag) cout << "\n Writing AdvectDiffuse2D node locations to output data file(s).";
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
	 cout << "\n AdvectDiffuse2D ERROR: Unable to open AdvectDiffuse2D nodes output data file(s) "
	      << "on processor "
	      << List_of_Local_Solution_Blocks.ThisCPU
	      << "." << endl;
       } /* endif */
       error_flag = CFDkit_OR_MPI(error_flag);
       if (error_flag) return error_flag;

     } else if (command_flag == WRITE_RESTART_CODE) {
         // Write restart files.
         if (!batch_flag) cout << "\n Writing AdvectDiffuse2D solution to restart data file(s).";
         error_flag = Write_QuadTree(QuadTree,
                                     Input_Parameters);
         if (error_flag) {
            cout << "\n AdvectDiffuse2D ERROR: Unable to open AdvectDiffuse2D quadtree data file "
                 << "on processor "
                 << List_of_Local_Solution_Blocks.ThisCPU
                 << ".\n";
            cout.flush();
         } /* endif */
         error_flag = CFDkit_OR_MPI(error_flag);
         if (error_flag) return (error_flag);
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
         error_flag = CFDkit_OR_MPI(error_flag);
         if (error_flag) return (error_flag);

     } else if (command_flag == WRITE_OUTPUT_GRID_CODE) {
         // Output multi-block solution-adaptive mesh data file.
         if (CFDkit_Primary_MPI_Processor()) {
            if (!batch_flag) cout << "\n Writing AdvectDiffuse2D multi-block mesh to grid data output file.";
            error_flag = Output_Tecplot(MeshBlk,
                                        Input_Parameters);
            if (error_flag) {
               cout << "\n AdvectDiffuse2D ERROR: Unable to open AdvectDiffuse2D mesh data output file.\n";
               cout.flush();
            } /* endif */
         } /* endif */
         CFDkit_Broadcast_MPI(&error_flag, 1);
         if (error_flag) return (error_flag);

     } else if (command_flag == WRITE_GRID_DEFINITION_CODE) {
         // Write multi-block solution-adaptive mesh definition files.
         if (CFDkit_Primary_MPI_Processor()) {
            if (!batch_flag) cout << "\n Writing AdvectDiffuse2D multi-block mesh to grid definition files.";
            error_flag = Write_Multi_Block_Grid_Definition(MeshBlk,
                                                           Input_Parameters);
            error_flag = Write_Multi_Block_Grid(MeshBlk,
                                                Input_Parameters);
            if (error_flag) {
               cout << "\n AdvectDiffuse2D ERROR: Unable to open AdvectDiffuse2D multi-block mesh definition files.\n";
               cout.flush();
            } /* endif */
         } /* endif */
         CFDkit_Broadcast_MPI(&error_flag, 1);
         if (error_flag) return (error_flag);

     } else if (command_flag == WRITE_OUTPUT_GRID_NODES_CODE) {
         // Output multi-block solution-adaptive mesh node data file.
         if (CFDkit_Primary_MPI_Processor()) {
            if (!batch_flag) cout << "\n Writing AdvectDiffuse2D multi-block mesh to node data output file.";
            error_flag = Output_Nodes_Tecplot(MeshBlk,
                                              Input_Parameters);
            if (error_flag) {
               cout << "\n AdvectDiffuse2D ERROR: Unable to open AdvectDiffuse2D mesh node data output file.\n";
               cout.flush();
            } /* endif */
         } /* endif */
         CFDkit_Broadcast_MPI(&error_flag, 1);
         if (error_flag) return (error_flag);

     } else if (command_flag == WRITE_OUTPUT_GRID_CELLS_CODE) {
         // Output multi-block solution-adaptive mesh cell data file.
         if (CFDkit_Primary_MPI_Processor()) {
            if (!batch_flag) cout << "\n Writing AdvectDiffuse2D multi-block mesh to cell data output file.";
            error_flag = Output_Cells_Tecplot(MeshBlk,
                                              Input_Parameters);
            if (error_flag) {
               cout << "\n AdvectDiffuse2D ERROR: Unable to open AdvectDiffuse2D mesh cell data output file.\n";
               cout.flush();
            } /* endif */
         } /* endif */
         CFDkit_Broadcast_MPI(&error_flag, 1);
         if (error_flag) return (error_flag);

     } else if (command_flag == INVALID_INPUT_CODE ||
                command_flag == INVALID_INPUT_VALUE) {
         line_number = -line_number;
         cout << "\n AdvectDiffuse2D ERROR: Error reading AdvectDiffuse2D data at line #"
	      << -line_number  << " of input data file.";
         cout.flush();
         return (line_number);
     } /* endif */

  } /* endwhile */

  /**********************************************************  
   * End of all AdvectDiffuse2DSolver computations and I/O. *
   **********************************************************/

  return (0);
  
}
