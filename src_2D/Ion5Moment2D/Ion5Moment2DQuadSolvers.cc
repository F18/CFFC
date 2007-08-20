/* Ion5Moment2DQuadSolvers.cc:  2D 5-Moment Ion Transport Model
                                Multi-Block Quadrilateral Mesh Solvers. */

/* Include 2D 5-moment ion transport model quadrilateral mesh solution header file. */

#ifndef _ION5MOMENT2D_QUAD_INCLUDED
#include "Ion5Moment2DQuad.h"
#endif // _ION5MOMENT2D_QUAD_INCLUDED

/********************************************************
 * Routine: Ion5Moment2DQuadSolver                      *
 *                                                      *
 * Computes solutions to 2D 5-moment ion transport 2D   *
 * model on 2D quadrilateral multi-block solution-      *
 * adaptive mesh.                                       *
 *                                                      *
 ********************************************************/
int Ion5Moment2DQuadSolver(char *Input_File_Name_ptr,
                           int batch_flag) {

  /********************************************************  
   * Local variable declarations.                         *
   ********************************************************/

  // Ion5Moment2D input variables and parameters:
  Ion5Moment2D_Input_Parameters Input_Parameters;

  /* Multi-block solution-adaptive quadrilateral mesh 
     solution variables. */
 
  Grid2D_Quad_Block             **MeshBlk;
  QuadTreeBlock_DataStructure     QuadTree;
  AdaptiveBlockResourceList       List_of_Global_Solution_Blocks;
  AdaptiveBlock2D_List            List_of_Local_Solution_Blocks;
  Ion5Moment2D_Quad_Block        *Local_SolnBlk;

//   FAS_Multigrid2D_Solver<Ion5Moment2D_cState, 
//                          Ion5Moment2D_Quad_Block, 
//                          Ion5Moment2D_Input_Parameters> MGSolver;

  /* Define residual file and cpu time variables. */

  ofstream residual_file;
  CPUTime processor_cpu_time, total_cpu_time;

  /* Other local solution variables. */

  int number_of_time_steps, first_step,
      command_flag, error_flag, line_number, 
      i_stage, i_blk;

  double Time, dTime;

  double residual_l1_norm, residual_l2_norm, residual_max_norm;

  /********************************************************  
   * Set default values for the input solution parameters *
   * and then read user specified input values from the   *
   * specified input parameter file.                      *
   ********************************************************/

  // The primary MPI processor processes the input parameter file.
  if (CFFC_Primary_MPI_Processor()) {
     if (!batch_flag) {
        cout << "\n Reading Ion5Moment2D input data file `"
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
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
  CFFC_Broadcast_MPI(&error_flag, 1);
  if (error_flag != 0) return (error_flag);
  CFFC_Broadcast_MPI(&command_flag, 1);
  if (command_flag == TERMINATE_CODE) return (0);
  Broadcast_Input_Parameters(Input_Parameters);

  /********************************************************  
   * Create initial mesh and allocate Ion5Moment2D        *
   * solution variables for specified IBVP/BVP problem.   *
   ********************************************************/

  execute_new_calculation: ;
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

  /* Create initial mesh.  Read mesh from grid definition or data files 
     when specified by input parameters. */

  // The primary MPI processor creates the initial mesh.
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
        cout << "\n Ion5Moment2D ERROR: Unable to create valid Ion5Moment2D multi-block mesh.\n";
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
     (allocate) array of local 2D 5-moment ion transport model solution blocks, 
     assign and create (allocate) 5-moment ion transport model solution blocks
     corresponding to the initial mesh. */

  if (!batch_flag) cout << "\n Creating multi-block quadtree data structure and assigning"
                        << "\n  Ion5Moment2D solution blocks corresponding to initial mesh.";
  Local_SolnBlk = Create_Initial_Solution_Blocks(MeshBlk,
                                                 Local_SolnBlk,
                                                 Input_Parameters,
                                                 QuadTree,
                                                 List_of_Global_Solution_Blocks,
                                                 List_of_Local_Solution_Blocks);
  if (Local_SolnBlk == NULL) return (1);

  /********************************************************  
   * Initialize Ion5Moment2D solution variables.          *
   ********************************************************/

  /* Set the initial time level. */

  Time = ZERO;
  number_of_time_steps = 0;

 /* Set the CPU time to zero. */

  processor_cpu_time.zero();
  total_cpu_time.zero();
  
  /* Initialize the conserved and primitive state
     solution variables. */
  
  if (!batch_flag) cout << "\n Prescribing Ion5Moment2D initial data.";
  if (Input_Parameters.i_ICs == IC_RESTART) {
     if (!batch_flag) cout << "\n Reading Ion5Moment2D solution from restart data files.";
     error_flag = Read_Restart_Solution(Local_SolnBlk, 
                                        List_of_Local_Solution_Blocks, 
                                        Input_Parameters,
                                        number_of_time_steps,
                                        Time,
                                        processor_cpu_time);
     if (error_flag) {
        cout << "\n Ion5Moment2D ERROR: Unable to open Ion5Moment2D restart input data file(s) "
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
  } else {
     // Set initial data using subroutine ICs.
     ICs(Local_SolnBlk, 
         List_of_Local_Solution_Blocks, 
         Input_Parameters);

     // Read the neutral gas flow solution from appropriate input file as required.
     if (strcmp(Input_Parameters.Neutral_Gas_Solution_File_Name, "none") != 0 &&
         strcmp(Input_Parameters.Neutral_Gas_Solution_File_Name, "NONE") != 0) {
        if (!batch_flag) cout << "\n Reading Ion5Moment2D neutral gas solution file.";
        for (i_blk = 0; i_blk <= Input_Parameters.Number_of_Processors-1; ++i_blk) {
           CFFC_Barrier_MPI(); // Put a barrier here to ensure that only one processor is accessing input file at a time.
           if (i_blk == List_of_Local_Solution_Blocks.ThisCPU) {
              error_flag = Read_Neutral_Gas_Solution(Local_SolnBlk, 
                                                     List_of_Local_Solution_Blocks, 
                                                     Input_Parameters);
              BCs_Neutral_Gas(Local_SolnBlk, 
                              List_of_Local_Solution_Blocks); // Apply boundary conditions for neutral gas flow.
           } /* endif */
        } /* endif */
     } /* endif */

     // Read the electric field solution from appropriate input file as required.
     if (strcmp(Input_Parameters.Electric_Field_Solution_File_Name, "none") != 0 &&
         strcmp(Input_Parameters.Electric_Field_Solution_File_Name, "NONE") != 0) {
        if (!batch_flag) cout << "\n Reading Ion5Moment2D electric field solution file.";
        for (i_blk = 0; i_blk <= Input_Parameters.Number_of_Processors-1; ++i_blk) {
           CFFC_Barrier_MPI(); // Put a barrier here to ensure that only one processor is accessing input file at a time.
           if (i_blk == List_of_Local_Solution_Blocks.ThisCPU) {
              error_flag = Read_Electric_Field_Solution(Local_SolnBlk, 
                                                        List_of_Local_Solution_Blocks, 
                                                        Input_Parameters);
              BCs_Electric_Field(Local_SolnBlk, 
                                 List_of_Local_Solution_Blocks); // Apply boundary conditions for electric field.
           } /* endif */
        } /* endif */
     } /* endif */
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
                                                  NUM_VAR_ION5MOMENT2D + NUM_VAR_EULER2D + NUM_COMP_VECTOR2D + 1,
                                                  OFF);

  if (error_flag) {
     cout << "\n Ion5Moment2D ERROR: Message passing error during Ion5Moment2D solution intialization "
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
     cout << "\n  -> Number of Refinement Levels: " 
	  << QuadTree.highestRefinementLevel();
     cout << "\n  -> Refinement Efficiency: " 
          << QuadTree.efficiencyRefinement() << "\n"; 
     cout.flush();
  } /* endif */

  /********************************************************  
   * Solve IBVP or BVP for conservation form of 2D        *
   * 5-moment ion transport model equations on            *
   * multi-block solution-adaptive quadrilateral mesh.    *
   ********************************************************/

  continue_existing_calculation: ;
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

  if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {

//     MGSolver.allocate(batch_flag,
// 		      Local_SolnBlk,
// 		      &QuadTree,
// 		      &List_of_Local_Solution_Blocks,
// 		      &Input_Parameters,
// 		      2*NUM_COMP_VECTOR2D + NUM_VAR_EULER2D + 1);

//     error_flag = MGSolver.Execute(number_of_time_steps,
// 				  Time,
// 				  processor_cpu_time,
// 				  total_cpu_time,
// 				  residual_file);

  } else {

    /* Open residual file and reset the CPU time. */
    
    first_step = 1;
    
    if (CFFC_Primary_MPI_Processor()) {
      error_flag = Open_Progress_File(residual_file,
				      Input_Parameters.Output_File_Name,
				      number_of_time_steps);
      if (error_flag) {
        cout << "\n Ion5Moment2D ERROR: Unable to open residual file for Ion5Moment2D calculation.\n";
        cout.flush();
      } /* endif */
    } /* endif */
    
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
    CFFC_Broadcast_MPI(&error_flag, 1);
    if (error_flag) return (error_flag);
    
    processor_cpu_time.reset();
    
    /* Perform required number of iterations (time steps). */
    
    if ((!Input_Parameters.Time_Accurate &&
	 Input_Parameters.Maximum_Number_of_Time_Steps > 0) ||
	(Input_Parameters.Time_Accurate &&
	 Input_Parameters.Time_Max > Time)) {
      if (!batch_flag) cout << "\n\n Beginning Ion5Moment2D computations on "
			    << Date_And_Time() << ".\n\n";
      while (1) {
	/* Periodically refine the mesh (AMR). */
        if (Input_Parameters.AMR) {
           if (!first_step &&
	       number_of_time_steps-Input_Parameters.AMR_Frequency*
	       (number_of_time_steps/Input_Parameters.AMR_Frequency) == 0 ) {
              if (!batch_flag) cout << "\n\n Refining Grid.  Performing adaptive mesh refinement at n = "
                                    << number_of_time_steps << ".";
//               Evaluate_Limiters(Local_SolnBlk, 
//                                 List_of_Local_Solution_Blocks);
              error_flag = AMR(Local_SolnBlk,
			       Input_Parameters,
                               QuadTree,
                               List_of_Global_Solution_Blocks,
                               List_of_Local_Solution_Blocks,
                               ON,ON);
              if (error_flag) {
                 cout << "\n Ion5Moment2D ERROR: Ion5Moment2D AMR error on processor "
	              << List_of_Local_Solution_Blocks.ThisCPU
	              << ".\n";
	         cout.flush();
              } /* endif */
              error_flag = CFFC_OR_MPI(error_flag);
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
	residual_l1_norm = L1_Norm_Residual(Local_SolnBlk, List_of_Local_Solution_Blocks,Input_Parameters);
	residual_l1_norm = CFFC_Summation_MPI(residual_l1_norm); // L1 norm for all processors.
	
	residual_l2_norm = L2_Norm_Residual(Local_SolnBlk, List_of_Local_Solution_Blocks,Input_Parameters);
	residual_l2_norm = sqr(residual_l2_norm);
	residual_l2_norm = CFFC_Summation_MPI(residual_l2_norm); // L2 norm for all processors.
	residual_l2_norm = sqrt(residual_l2_norm);
	
	residual_max_norm = Max_Norm_Residual(Local_SolnBlk, List_of_Local_Solution_Blocks,Input_Parameters);
	residual_max_norm = CFFC_Maximum_MPI(residual_max_norm); // Max norm for all processors.
	
	/* Update CPU time used for the calculation so far. */
	processor_cpu_time.update();
	total_cpu_time.cput = 
	  CFFC_Summation_MPI(processor_cpu_time.cput); // Total CPU time for all processors.
	
	/* Periodically save restart solution files. */
	if (!first_step &&
	    number_of_time_steps-Input_Parameters.Restart_Solution_Save_Frequency*
	    (number_of_time_steps/Input_Parameters.Restart_Solution_Save_Frequency) == 0 ) {
	  if (!batch_flag) cout << "\n\n  Saving Ion5Moment2D solution to restart data file(s) after"
				<< " n = " << number_of_time_steps << " steps (iterations).";
	  error_flag = Write_Restart_Solution(Local_SolnBlk, 
					      List_of_Local_Solution_Blocks, 
					      Input_Parameters,
					      number_of_time_steps,
					      Time,
					      processor_cpu_time);
	  if (error_flag) {
	    cout << "\n Ion5Moment2D ERROR: Unable to open Ion5Moment2D restart output data file(s) "
		 << "on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  } /* endif */
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return (error_flag);
	  cout << "\n";
	  cout.flush();
	} /* endif */
	
	/* Output progress information for the calculation. */
	if (!batch_flag) Output_Progress(number_of_time_steps,
					 Time*THOUSAND,
					 total_cpu_time,
					 residual_l1_norm,
					 first_step,
					 50);
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
	
	/* Update solution for next time step using a multistage
	   time stepping scheme. */
	for ( i_stage  = 1 ; i_stage <= Input_Parameters.N_Stage ; ++i_stage ) {
	  // 1. Exchange solution information between neighbouring blocks.
	  error_flag = Send_All_Messages(Local_SolnBlk, 
					 List_of_Local_Solution_Blocks,
					 NUM_VAR_ION5MOMENT2D + NUM_VAR_EULER2D + NUM_COMP_VECTOR2D + 1, 
					 OFF);
	  if (error_flag) {
	    cout << "\n Ion5Moment2D ERROR: Ion5Moment2D message passing error on processor "
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
	    cout << "\n Ion5Moment2D ERROR: Ion5Moment2D solution residual error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  } /* endif */
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return (error_flag);
	  
	  // 4. Send boundary flux corrections at block interfaces with resolution changes.
	  error_flag = Send_Conservative_Flux_Corrections(Local_SolnBlk, 
							  List_of_Local_Solution_Blocks,
							  NUM_VAR_ION5MOMENT2D + NUM_VAR_EULER2D + NUM_COMP_VECTOR2D + 1);
	  if (error_flag) {
	    cout << "\n Ion5Moment2D ERROR: Ion5Moment2D flux correction message passing error on processor "
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
	    cout << "\n Ion5Moment2D ERROR: Ion5Moment2D solution update error on processor "
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
      if (!batch_flag) cout << "\n\n Ion5Moment2D computations complete on " 
			    << Date_And_Time() << ".\n";
    } /* endif */
    
    /* Update ghostcell information and prescribe boundary conditions to ensure
       that the solution is consistent on each block. */
    
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
    
    error_flag = Send_All_Messages(Local_SolnBlk, 
				   List_of_Local_Solution_Blocks,
				   NUM_VAR_ION5MOMENT2D + NUM_VAR_EULER2D + NUM_COMP_VECTOR2D + 1,
				   OFF);
    
    if (error_flag) {
      cout << "\n Ion5Moment2D ERROR: Ion5Moment2D message passing error on processor "
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

  /********************************************************
   * Solution calculations complete.                      *
   * Write 2D 5-moment ion transport model solution to    *
   * output and restart files as required, reset solution *
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
         // Deallocate memory for 2D 5-moment ion transport model solution.
         if (!batch_flag) cout << "\n Deallocating Ion5Moment2D solution variables.";
	 if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {
	   //MGSolver.deallocate();
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
         // Deallocate memory for 2D 5-moment ion transport model solution.
         if (!batch_flag) cout << "\n Deallocating Ion5Moment2D solution variables.";
	 if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {
	   //MGSolver.deallocate();
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
         if (!batch_flag) cout << "\n\n Closing Ion5Moment2D input data file.";
         if (CFFC_Primary_MPI_Processor()) Close_Input_File(Input_Parameters);
         // Terminate calculation.
         return (0);

     } else if (command_flag == CONTINUE_CODE) {
         // Reset maximum time step counter.
         Input_Parameters.Maximum_Number_of_Time_Steps += number_of_time_steps;
         // Output input parameters for continuing calculation.
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
         error_flag = AMR(Local_SolnBlk,
			  Input_Parameters,
                          QuadTree,
                          List_of_Global_Solution_Blocks,
                          List_of_Local_Solution_Blocks,
                          ON,ON);
         if (error_flag) {
            cout << "\n Ion5Moment2D ERROR: Ion5Moment2D AMR error on processor "
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

     } else if (command_flag == WRITE_OUTPUT_CODE) {
         // Output solution data.
         if (!batch_flag) cout << "\n Writing Ion5Moment2D solution to output data file(s).";
         error_flag = Output_Tecplot(Local_SolnBlk, 
                                     List_of_Local_Solution_Blocks, 
                                     Input_Parameters,
                                     number_of_time_steps,
                                     Time);
         if (error_flag) {
            cout << "\n Ion5Moment2D ERROR: Unable to open Ion5Moment2D output data file(s) "
                 << "on processor "
                 << List_of_Local_Solution_Blocks.ThisCPU
                 << ".\n";
            cout.flush();
         } /* endif */
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);

     } else if (command_flag == WRITE_OUTPUT_CELLS_CODE) {
         // Output solution data.
         if (!batch_flag) cout << "\n Writing cell-centered Ion5Moment2D solution to output data file(s).";
         error_flag = Output_Cells_Tecplot(Local_SolnBlk, 
                                           List_of_Local_Solution_Blocks, 
                                           Input_Parameters,
                                           number_of_time_steps,
                                           Time);
         if (!error_flag) error_flag = dUdt_Output_Cells_Tecplot(Local_SolnBlk,
                                                                 List_of_Local_Solution_Blocks,
                                                                 Input_Parameters,
                                                                 number_of_time_steps,
                                                                 Time);
         if (error_flag) {
            cout << "\n Ion5Moment2D ERROR: Unable to open Ion5Moment2D cell output data file(s) "
                 << "on processor "
                 << List_of_Local_Solution_Blocks.ThisCPU
                 << ".\n";
            cout.flush();
         } /* endif */
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);

     } else if (command_flag == WRITE_RESTART_CODE) {
         // Write restart files.
         if (!batch_flag) cout << "\n Writing Ion5Moment2D solution to restart data file(s).";
         error_flag = Write_Restart_Solution(Local_SolnBlk, 
                                             List_of_Local_Solution_Blocks, 
                                             Input_Parameters,
                                             number_of_time_steps,
                                             Time,
                                             processor_cpu_time);
         if (error_flag) {
            cout << "\n Ion5Moment2D ERROR: Unable to open Ion5Moment2D restart output data file(s) "
                 << "on processor "
                 << List_of_Local_Solution_Blocks.ThisCPU
                 << ".\n";
            cout.flush();
         } /* endif */
         error_flag = CFFC_OR_MPI(error_flag);
         if (error_flag) return (error_flag);

     } else if (command_flag == WRITE_OUTPUT_GRID_CODE) {
         // Output multi-block solution-adaptive mesh data file.
         if (CFFC_Primary_MPI_Processor()) {
            if (!batch_flag) cout << "\n Writing Ion5Moment2D multi-block mesh to grid data output file.";
            error_flag = Output_Tecplot(MeshBlk,
                                        Input_Parameters);
            if (error_flag) {
               cout << "\n Ion5Moment2D ERROR: Unable to open Ion5Moment2D mesh data output file.\n";
               cout.flush();
            } /* endif */
         } /* endif */
         CFFC_Broadcast_MPI(&error_flag, 1);
         if (error_flag) return (error_flag);

     } else if (command_flag == WRITE_GRID_DEFINITION_CODE) {
         // Write multi-block solution-adaptive mesh definition files.
         if (CFFC_Primary_MPI_Processor()) {
            if (!batch_flag) cout << "\n Writing Ion5Moment2D multi-block mesh to grid definition files.";
            error_flag = Write_Multi_Block_Grid_Definition(MeshBlk,
                                                           Input_Parameters);
            error_flag = Write_Multi_Block_Grid(MeshBlk,
                                                Input_Parameters);
            if (error_flag) {
               cout << "\n Ion5Moment2D ERROR: Unable to open Ion5Moment2D multi-block mesh definition files.\n";
               cout.flush();
            } /* endif */
         } /* endif */
         CFFC_Broadcast_MPI(&error_flag, 1);
         if (error_flag) return (error_flag);

     } else if (command_flag == WRITE_OUTPUT_GRID_NODES_CODE) {
         // Output multi-block solution-adaptive mesh node data file.
         if (CFFC_Primary_MPI_Processor()) {
            if (!batch_flag) cout << "\n Writing Ion5Moment2D multi-block mesh to node data output file.";
            error_flag = Output_Nodes_Tecplot(MeshBlk,
                                              Input_Parameters);
            if (error_flag) {
               cout << "\n Ion5Moment2D ERROR: Unable to open Ion5Moment2D mesh node data output file.\n";
               cout.flush();
            } /* endif */
         } /* endif */
         CFFC_Broadcast_MPI(&error_flag, 1);
         if (error_flag) return (error_flag);

     } else if (command_flag == WRITE_OUTPUT_GRID_CELLS_CODE) {
         // Output multi-block solution-adaptive mesh cell data file.
         if (CFFC_Primary_MPI_Processor()) {
            if (!batch_flag) cout << "\n Writing Ion5Moment2D multi-block mesh to cell data output file.";
            error_flag = Output_Cells_Tecplot(MeshBlk,
                                              Input_Parameters);
            if (error_flag) {
               cout << "\n Ion5Moment2D ERROR: Unable to open Ion5Moment2D mesh cell data output file.\n";
               cout.flush();
            } /* endif */
         } /* endif */
         CFFC_Broadcast_MPI(&error_flag, 1);
         if (error_flag) return (error_flag);

     } else if (command_flag == INVALID_INPUT_CODE ||
                command_flag == INVALID_INPUT_VALUE) {
         line_number = -line_number;
         cout << "\n Ion5Moment2D ERROR: Error reading Ion5Moment2D data at line #"
              << -line_number  << " of input data file.\n";
         cout.flush();
         return (line_number);
     } /* endif */

  } /* endwhile */

  /********************************************************  
   * End of all Ion5Moment2DSolver computations and I/O.  *
   ********************************************************/

  return (0);
  
}
