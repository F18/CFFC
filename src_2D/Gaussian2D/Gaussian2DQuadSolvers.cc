/* Gaussian2DQuadSolvers.cc:  2D Gaussian Equation 
                              Multi-Block Quadrilateral Mesh Solvers. */

/* Include 2D Gaussian quadrilateral mesh solution header file. */

#ifndef _GAUSSIAN2D_QUAD_INCLUDED
#include "Gaussian2DQuad.h"
#endif // _GAUSSIAN2D_QUAD_INCLUDED

/********************************************************
 * Routine: Gaussian2DQuadSolver                        *
 *                                                      *
 * Computes solutions to 2D Gaussian equations on 2D    *
 * quadrilateral multi-block solution-adaptive mesh.    *
 *                                                      *
 ********************************************************/
int Gaussian2DQuadSolver(char *Input_File_Name_ptr,
                         int batch_flag) {

  /********************************************************  
   * Local variable declarations.                         *
   ********************************************************/

  // Gaussian2D input variables and parameters:
  Gaussian2D_Input_Parameters Input_Parameters;

  /* Multi-block solution-adaptive quadrilateral mesh 
     solution variables. */
 
  Grid2D_Quad_Block          **MeshBlk;
  QuadTreeBlock_DataStructure  QuadTree;
  AdaptiveBlockResourceList    List_of_Global_Solution_Blocks;
  AdaptiveBlock2D_List         List_of_Local_Solution_Blocks;
  Gaussian2D_Quad_Block       *Local_SolnBlk;

  /* Define residual file and cpu time variables. */

  ofstream residual_file, drag_file;
  CPUTime processor_cpu_time, total_cpu_time;

  /* Other local solution variables. */

  int number_of_time_steps, first_step,
      command_flag, error_flag, line_number, 
      i_stage, perform_explicit_time_marching, limiter_freezing_off;

  double Time, dTime;

  double residual_l2norm_first, residual_ratio;
  double residual_l1_norm, residual_l2_norm, residual_max_norm;
  double drag, Cd, S, Kn, Re, CdC, CdF, CdP;
  double lift, Cl, speed;

  /********************************************************  
   * Set default values for the input solution parameters *
   * and then read user specified input values from the   *
   * specified input parameter file.                      *
   ********************************************************/

  // The primary MPI processor processes the input parameter file.
  if (CFFC_Primary_MPI_Processor()) {
     if (!batch_flag) {
        cout << "\n Reading Gaussian2D input data file `"
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
   * Create initial mesh and allocate Gaussian2D solution *
   * variables for specified IBVP/BVP problem.            *
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
        cout << "\n Gaussian2D ERROR: Unable to create valid Gaussian2D multi-block mesh.\n";
        cout.flush();
     }
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
     (allocate) array of local 2D Gaussian equation solution blocks, 
     assign and create (allocate) 2D Gaussian equation solution blocks
     corresponding to the initial mesh. */

  if (!batch_flag) cout << "\n Creating multi-block quadtree data structure and assigning"
                        << "\n  Gaussian2D solution blocks corresponding to initial mesh.";
  Local_SolnBlk = Create_Initial_Solution_Blocks(MeshBlk,
                                                 Local_SolnBlk,
                                                 Input_Parameters,
                                                 QuadTree,
                                                 List_of_Global_Solution_Blocks,
                                                 List_of_Local_Solution_Blocks);
  if (Local_SolnBlk == NULL) return (1);

  /********************************************************  
   * Initialize Gaussian2D solution variables.            *
   ********************************************************/

  /* Set the initial time level. */

  Time = ZERO;
  number_of_time_steps = 0;

 /* Set the CPU time to zero. */

  processor_cpu_time.zero();
  total_cpu_time.zero();
  
  /* Initialize the conserved and primitive state
     solution variables. */
  
  if (!batch_flag) cout << "\n Prescribing Gaussian2D initial data.";
  if (Input_Parameters.i_ICs == IC_RESTART) {
     if (!batch_flag) cout << "\n Reading Gaussian2D solution from restart data files.";

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
        cout << "\n Gaussian2D ERROR: Unable to open Gaussian2D quadtree data file "
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
        cout << "\n Gaussian2D ERROR: Unable to open Gaussian2D restart input data file(s) "
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
                                                  NUM_VAR_GAUSSIAN2D,
                                                  OFF);

  if (error_flag) {
     cout << "\n Gaussian2D ERROR: Message passing error during Gaussian2D solution intialization "
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
     if (!batch_flag) cout << "\n Performing Gaussian2D uniform mesh refinement.";
     error_flag = Uniform_AMR(Local_SolnBlk,
                              Input_Parameters,
                              QuadTree,
                              List_of_Global_Solution_Blocks,
                              List_of_Local_Solution_Blocks);
     if (error_flag) {
       cout << "\n Gaussian2D ERROR: Uniform AMR error on processor "
	    << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
       cout.flush();
     } /* endif */
     error_flag = CFFC_OR_MPI(error_flag);
     if (error_flag) return error_flag;

     if (!batch_flag) cout << "\n Performing Gaussian2D boundary mesh refinement.";
     error_flag = Boundary_AMR(Local_SolnBlk,
                               Input_Parameters,
                               QuadTree,
                               List_of_Global_Solution_Blocks,
                               List_of_Local_Solution_Blocks);
     if (error_flag) {
       cout << "\n Gaussian2D ERROR: Boundary AMR error on processor "
	    << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
       cout.flush();
     } /* endif */
     error_flag = CFFC_OR_MPI(error_flag);
     if (error_flag) return error_flag;


     if (!batch_flag) cout << "\n Performing Gaussian2D initial mesh refinement.";
     error_flag = Initial_AMR(Local_SolnBlk,
                              Input_Parameters,
	   	              QuadTree,
		              List_of_Global_Solution_Blocks,
		              List_of_Local_Solution_Blocks);
     if (error_flag) {
        cout << "\n Gaussian2D ERROR: Initial AMR error on processor "
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

  /***********************************************************  
   * Solve IBVP or BVP for conservation form of 2D Gaussian  *
   * equations on multi-block solution-adaptive              *
   * quadrilateral mesh.                                     *
   ***********************************************************/

  continue_existing_calculation: ;
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

  /* Open residual file and reset the CPU time. */
    
  first_step = 1;
  limiter_freezing_off = ON;
    
  if (CFFC_Primary_MPI_Processor()) {
    error_flag = Open_Progress_File(residual_file,
				    Input_Parameters.Output_File_Name,
				    number_of_time_steps);
    drag_file.open("lift_and_drag.dat", ios::out|ios::app);              //DRAG_FILE STUFF
    error_flag = drag_file.fail();                                       //DRAG_FILE STUFF
    if (error_flag) {
      cout << "\n Gaussian2D ERROR: Unable to open residual file for Gaussian2D calculation.\n";
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
      (!Input_Parameters.Time_Accurate) ||
      (Input_Parameters.Time_Accurate &&
       Input_Parameters.Time_Max > Time)) {
    if (!batch_flag) cout << "\n\n Beginning Gaussian2D computations on "
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
                 cout << "\n Gaussian2D ERROR: Gaussian2D AMR error on processor "
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

	if (!batch_flag) cout << "\n\n  Saving Gaussian2D solution to restart data file(s) after"
			      << " n = " << number_of_time_steps << " steps (iterations).";
	error_flag = Write_Restart_Solution(Local_SolnBlk, 
					    List_of_Local_Solution_Blocks, 
					    Input_Parameters,
					    number_of_time_steps,
					    Time,
					    processor_cpu_time);
	if (error_flag) {
	  cout << "\n Gaussian2D ERROR: Unable to open Gaussian2D restart output data file(s) "
	       << "on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU
	       << ".\n";
	  cout.flush();
	} /* endif */

        if (CFFC_Primary_MPI_Processor()) {
	  if (!batch_flag) cout << "\n  Saving QuadTree data file.";
 	  error_flag = Write_QuadTree(QuadTree, Input_Parameters);
          if (error_flag) {
	    cout << "\n Gaussian2D ERROR: Unable to open QuadTree data file(s).\n";
	    cout.flush();
	  } /* endif */
	} /* endif */

//       if (CFFC_Primary_MPI_Processor()) {
//          if (!batch_flag) cout << "\n  Writing Gaussian2D multi-block mesh to grid definition files.";
//          error_flag = Write_Multi_Block_Grid_Definition(MeshBlk,
//                                                         Input_Parameters);
//          error_flag = Write_Multi_Block_Grid(MeshBlk,
//                                              Input_Parameters);
//          if (error_flag) {
//             cout << "\n Gaussian2D ERROR: Unable to open Gaussian2D multi-block mesh definition files.\n";
//             cout.flush();
//          } /* endif */
//       } /* endif */

       CFFC_Broadcast_MPI(&error_flag, 1);
       if (error_flag) return (error_flag);

	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return (error_flag);
	if (!batch_flag) cout << "\n";
	cout.flush();

	if (CFFC_Primary_MPI_Processor()) {
	  System::Remove_Restart_Flag();  //Remove flag to indicate the restart is finished
	}

      } /* endif */  //end save_restart
	
      /* Output progress information for the calculation. */
      if (!batch_flag) Output_Progress_L2norm(number_of_time_steps,
					      Time*THOUSAND,
					      total_cpu_time,
					      residual_l2_norm,
					      first_step,
					      50);
      //Output Cd and Cl for now
      
      drag = 0.0;                                                                                  //Drag File Stuff!!
      lift = 0.0;
      Output_Drag(Local_SolnBlk,
		  List_of_Local_Solution_Blocks,
		  drag,lift);
      drag = CFFC_Summation_MPI(drag);
      lift = CFFC_Summation_MPI(lift);
      speed = sqrt(sqr(Input_Parameters.Wo.v.x)+sqr(Input_Parameters.Wo.v.y))
	+Input_Parameters.Ramp_by_Mach_Number*Input_Parameters.Wo.sound();
      if(speed < TOLER){speed = TOLER;}
      Cd = drag/(0.5*Input_Parameters.Wo.d*sqr(speed)*2.0*Input_Parameters.Cylinder_Radius);
      Cl = lift/(0.5*Input_Parameters.Wo.d*sqr(speed)*2.0*Input_Parameters.Cylinder_Radius);
      if(!batch_flag){
	drag_file << Time << "     " << Cl << "     " << Cd << endl;
      }
      
      if ((number_of_time_steps+1)%50 == 0 ) {
//	drag = 0.0;
//	lift = 0.0;
//	Output_Drag(Local_SolnBlk,
//		    List_of_Local_Solution_Blocks,
//		    drag,lift);
//	drag = CFFC_Summation_MPI(drag);
//	lift = CFFC_Summation_MPI(lift);
//	speed = sqrt(sqr(Input_Parameters.Wo.v.x)+sqr(Input_Parameters.Wo.v.y))
//	        +Input_Parameters.Ramp_by_Mach_Number*Input_Parameters.Wo.sound();
//	Cd = drag/(0.5*Input_Parameters.Wo.d*sqr(speed)*2.0*Input_Parameters.Cylinder_Radius);
//	Cl = lift/(0.5*Input_Parameters.Wo.d*sqr(speed)*2.0*Input_Parameters.Cylinder_Radius);
	if(!batch_flag) {
	  cout << endl << "  Cd = " << Cd << "          Cl = " << Cl;
	}
      }

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
				       NUM_VAR_GAUSSIAN2D, 
				       OFF); 
	if (error_flag) {
	  cout << "\n Gaussian2D ERROR: Gaussian2D message passing error on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU
	       << ".\n";
	  cout.flush();
	} /* endif */
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return (error_flag);

	// 2. Apply boundary conditions for stage.

	if(number_of_time_steps < Input_Parameters.Number_of_Time_Steps_to_Ramp &&
	   i_stage == 1) {
	  Ramp_up_Reference_Mach_Number(Local_SolnBlk,
					List_of_Local_Solution_Blocks,
					Input_Parameters,
					number_of_time_steps);
	}

	BCs(Local_SolnBlk, 
	    List_of_Local_Solution_Blocks,
	    Input_Parameters);

	// 3. Determine solution residuals for stage.
	error_flag = dUdt_Multistage_Explicit(Local_SolnBlk, 
					      List_of_Local_Solution_Blocks,
					      Input_Parameters,
					      i_stage);
	if (error_flag) {
	  cout << "\n Gaussian2D ERROR: Gaussian2D solution residual error on processor "
	       << List_of_Local_Solution_Blocks.ThisCPU
	       << ".\n";
	  cout.flush();
	} /* endif */
	error_flag = CFFC_OR_MPI(error_flag);
	if (error_flag) return (error_flag);

	// 4. Send boundary flux corrections at block interfaces with resolution changes.
	error_flag = Send_Conservative_Flux_Corrections(Local_SolnBlk, 
							List_of_Local_Solution_Blocks,
							NUM_VAR_GAUSSIAN2D);
	if (error_flag) {
	  cout << "\n Gaussian2D ERROR: Gaussian2D flux correction message passing error on processor "
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

	// 6. Update solution for stage.
	error_flag = Update_Solution_Multistage_Explicit(Local_SolnBlk, 
							 List_of_Local_Solution_Blocks,
							 Input_Parameters,
							 i_stage);

	if (error_flag) {
	  cout << "\n Gaussian2D ERROR: Gaussian2D solution update error on processor "
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

    if (!batch_flag) cout << "\n\n Gaussian2D computations complete on " 
			  << Date_And_Time() << ".\n";
  } /* endif */
    
  /* Update ghostcell information and prescribe boundary conditions to ensure
       that the solution is consistent on each block. */
    
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
    
  error_flag = Send_All_Messages(Local_SolnBlk, 
				 List_of_Local_Solution_Blocks,
				 NUM_VAR_GAUSSIAN2D,
				 OFF);
  if (error_flag) {
    cout << "\n Gaussian2D ERROR: Gaussian2D message passing error on processor "
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
  if (CFFC_Primary_MPI_Processor()) drag_file.close();                      //Drag File


  /***********************************************************
   * Solution calculations complete.                         *
   * Write 2D Gaussian solution to output and restart files  *
   * as required, reset solution parameters, and run         *
   * other cases as specified by input parameters.           *
   ***********************************************************/
  
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
      // Deallocate memory for 2D Gaussian equation solution.
      if (!batch_flag) cout << "\n Deallocating Gaussian2D solution variables.";
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
      // Deallocate memory for 2D Gaussian equation solution.
      if (!batch_flag) cout << "\n Deallocating Gaussian2D solution variables.";
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
      if (!batch_flag) cout << "\n\n Closing Gaussian2D input data file.";
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
      Evaluate_Limiters(Local_SolnBlk, 
                        List_of_Local_Solution_Blocks);
      error_flag = AMR(Local_SolnBlk,
		       Input_Parameters,
		       QuadTree,
		       List_of_Global_Solution_Blocks,
		       List_of_Local_Solution_Blocks,
		       ON,ON);
      if (error_flag) {
	cout << "\n Gaussian2D ERROR: Gaussian2D AMR error on processor "
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
//       if (CFFC_Primary_MPI_Processor()) {
//          for ( int j_blk = 0 ; j_blk <= QuadTree.Nblk-1 ; ++j_blk ) {
//             for ( int i_blk = 0 ; i_blk <= QuadTree.Ncpu-1 ; ++i_blk ) {
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
//                } /* endif */
//             } /* endfor */
//          } /* endfor */
//       } /* endif */

    } else if (command_flag == WRITE_OUTPUT_CODE) {
       // Output solution data.
       if (!batch_flag) cout << "\n Writing Gaussian2D solution to output data file(s).";
       error_flag = Output_Tecplot(Local_SolnBlk, 
                                   List_of_Local_Solution_Blocks, 
                                   Input_Parameters,
                                   number_of_time_steps,
                                   Time);
       if (error_flag) {
          cout << "\n Gaussian2D ERROR: Unable to open Gaussian2D output data file(s) "
               << "on processor "
               << List_of_Local_Solution_Blocks.ThisCPU
               << ".\n";
          cout.flush();
       } /* endif */
       error_flag = CFFC_OR_MPI(error_flag);
       if (error_flag) return (error_flag);

   } else if (command_flag == WRITE_OUTPUT_CELLS_CODE) {
       // Output solution data.
       if (!batch_flag) cout << "\n Writing cell-centered Gaussian2D solution to output data file(s).";
       error_flag = Output_Cells_Tecplot(Local_SolnBlk, 
                                         List_of_Local_Solution_Blocks, 
                                         Input_Parameters,
                                         number_of_time_steps,
                                         Time);
       if (error_flag) {
          cout << "\n Gaussian2D ERROR: Unable to open Gaussian2D cell output data file(s) "
               << "on processor "
               << List_of_Local_Solution_Blocks.ThisCPU
               << ".\n";
          cout.flush();
       } /* endif */
       error_flag = CFFC_OR_MPI(error_flag);
       if (error_flag) return (error_flag);

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

       if (!batch_flag) cout << "\n Writing Gaussian2D solution to restart data file(s).";
       error_flag = Write_Restart_Solution(Local_SolnBlk, 
                                           List_of_Local_Solution_Blocks, 
                                           Input_Parameters,
                                           number_of_time_steps,
                                           Time,
                                           processor_cpu_time);
       if (error_flag) {
          cout << "\n Gaussian2D ERROR: Unable to open Gaussian2D restart output data file(s) "
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
          cout << "\n Gaussian2D ERROR: Unable to open Gaussian2D quadtree data file "
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
          if (!batch_flag) cout << "\n Writing Gaussian2D multi-block mesh to grid data output file.";
          error_flag = Output_Tecplot(MeshBlk,
                                      Input_Parameters);
          if (error_flag) {
             cout << "\n Gaussian2D ERROR: Unable to open Gaussian2D mesh data output file.\n";
             cout.flush();
          } /* endif */
       } /* endif */
       CFFC_Broadcast_MPI(&error_flag, 1);
       if (error_flag) return (error_flag);

    } else if (command_flag == WRITE_GRID_DEFINITION_CODE) {
       // Write multi-block solution-adaptive mesh definition files.
       if (CFFC_Primary_MPI_Processor()) {
          if (!batch_flag) cout << "\n Writing Gaussian2D multi-block mesh to grid definition files.";
          error_flag = Write_Multi_Block_Grid_Definition(MeshBlk,
                                                         Input_Parameters);
          error_flag = Write_Multi_Block_Grid(MeshBlk,
                                              Input_Parameters);
          if (error_flag) {
             cout << "\n Gaussian2D ERROR: Unable to open Gaussian2D multi-block mesh definition files.\n";
             cout.flush();
          } /* endif */
       } /* endif */
       CFFC_Broadcast_MPI(&error_flag, 1);
       if (error_flag) return (error_flag);

    } else if (command_flag == WRITE_OUTPUT_GRID_NODES_CODE) {
       // Output multi-block solution-adaptive mesh node data file.
       if (CFFC_Primary_MPI_Processor()) {
          if (!batch_flag) cout << "\n Writing Gaussian2D multi-block mesh to node data output file.";
          error_flag = Output_Nodes_Tecplot(MeshBlk,
                                            Input_Parameters);
          if (error_flag) {
             cout << "\n Gaussian2D ERROR: Unable to open Gaussian2D mesh node data output file.\n";
             cout.flush();
          } /* endif */
       } /* endif */
       CFFC_Broadcast_MPI(&error_flag, 1);
       if (error_flag) return (error_flag);

    } else if (command_flag == WRITE_OUTPUT_GRID_CELLS_CODE) {
       // Output multi-block solution-adaptive mesh cell data file.
       if (CFFC_Primary_MPI_Processor()) {
          if (!batch_flag) cout << "\n Writing Gaussian2D multi-block mesh to cell data output file.";
          error_flag = Output_Cells_Tecplot(MeshBlk,
                                            Input_Parameters);
          if (error_flag) {
             cout << "\n Gaussian2D ERROR: Unable to open Gaussian2D mesh cell data output file.\n";
             cout.flush();
          } /* endif */
       } /* endif */
       CFFC_Broadcast_MPI(&error_flag, 1);
       if (error_flag) return (error_flag);

    } else if (command_flag == WRITE_OUTPUT_FLAT_PLATE_CODE) {
       // Output Flat Plate Computed and Blasius solution
          if (!batch_flag) cout << "\n Writing Gaussian2D multi-block Flat Plate Solution....This may take a while.";
          error_flag = Output_Flat_Plate(Local_SolnBlk,
					 List_of_Local_Solution_Blocks,
					 Input_Parameters,
					 number_of_time_steps,
					 Time);
          if (error_flag) {
             cout << "\n Gaussian2D ERROR: Unable to open Gaussian2D flat plate data output file.\n";
             cout.flush();
       } /* endif */
       CFFC_Broadcast_MPI(&error_flag, 1);
       if (error_flag) return (error_flag);

    } else if (command_flag == WRITE_OUTPUT_SHOCK_STRUCTURE_CODE) {
       // Output shock structure solution
          if (!batch_flag) cout << "\n Writing Gaussian2D multi-block shock-structure solution.";
          error_flag = Output_Shock_Structure(Local_SolnBlk,
					      List_of_Local_Solution_Blocks,
					      Input_Parameters,
					      number_of_time_steps,
					      Time);
          if (error_flag) {
             cout << "\n Gaussian2D ERROR: Unable to open Gaussian2D shock-structure data output file.\n";
             cout.flush();
       } /* endif */
       CFFC_Broadcast_MPI(&error_flag, 1);
       if (error_flag) return (error_flag);

    } else if (command_flag == WRITE_OUTPUT_CYLINDER_FREE_MOLECULAR_CODE) {
       // Output Free Molecular Cylinder Solution
      if (!batch_flag) cout << "\n Writing Gaussian2D multi-block Free Molecular Cylinder Solution...This may take a while.";
      error_flag = Output_Cylinder_Free_Molecular(Local_SolnBlk,
						  List_of_Local_Solution_Blocks,
						  Input_Parameters,
						  number_of_time_steps,
						  Time);
      if (error_flag) {
	cout << "\n Gaussian2D ERROR: Unable to open Gaussian2D Free Molecular Cylinder data output file.\n";
	cout.flush();
      } /* endif */

      CFFC_Broadcast_MPI(&error_flag, 1);
      if (error_flag) return (error_flag);

    } else if (command_flag == WRITE_OUTPUT_DRAG_CODE) {
       // Output Drag
      drag = 0.0;
      lift = 0.0;
      if (!batch_flag) cout << "\n Writing Gaussian2D Drag.";
      Output_Drag(Local_SolnBlk,
		  List_of_Local_Solution_Blocks,
		  drag,lift);
      drag = CFFC_Summation_MPI(drag);
      lift = CFFC_Summation_MPI(lift);
      speed = sqrt(sqr(Input_Parameters.Wo.v.x)+sqr(Input_Parameters.Wo.v.y))
	      +Input_Parameters.Ramp_by_Mach_Number*Input_Parameters.Wo.sound();
      Cd = drag/(0.5*Input_Parameters.Wo.d*sqr(speed)*2.0*Input_Parameters.Cylinder_Radius);
      Cl = lift/(0.5*Input_Parameters.Wo.d*sqr(speed)*2.0*Input_Parameters.Cylinder_Radius);
      S = speed/sqrt(2.0*AVOGADRO*BOLTZMANN*Input_Parameters.Temperature/Input_Parameters.Wo.M*THOUSAND);
      Kn = Input_Parameters.Wo.mfp()/(2.0*Input_Parameters.Cylinder_Radius);
      Re = Input_Parameters.Wo.d*speed*2.0*Input_Parameters.Cylinder_Radius/Input_Parameters.Wo.viscosity();
      CdF = sqrt(PI)/S*(PI/4.0+3.0/2.0);
      CdC = 8.0*PI/Re/(log(8.0/Re)-0.0772);
      CdP = 8.0*PI/Re/(log(8.0/Re)-0.0772+1.0161*(4.0*S/Re)+0.0166*sqr(4.0*S/Re));
      if (!batch_flag) cout << "\n\n Drag = " << drag << " N/m (valid only for \"Adiabatic Walls\")        Cd = " 
			    << Cd <<" (valid only for cylinder)\n"
			    << " Lift = " << lift << " N/m (valid only for \"Adiabatic Walls\")        Cl = " 
			    << Cl <<" (valid only for cylinder)\n"
			    << " Speed Ratio = " << S << "                                       Kn = "
			    << Kn <<" (valid only for cylinder)\n           Re = " << Re
			    << "\n CdF = " << CdF << "                CdC = " << CdC
			    << " Valid only for Re<0.5\n Cd/CdF = " << Cd/CdF << "          CdC/CdF = " << CdC/CdF 
			    << "\n CdP = " << CdP << "                CdP/CdF = " << CdP/CdF << "\n";
    } else if (command_flag == INVALID_INPUT_CODE ||
               command_flag == INVALID_INPUT_VALUE) {
        line_number = -line_number;
        cout << "\n Gaussian2D ERROR: Error reading Gaussian2D data at line #"
             << -line_number  << " of input data file.\n";
        cout.flush();
        return (line_number);
    } /* endif */

  } /* endwhile */

  /********************************************************  
   * End of all Gaussian2DSolver computations and I/O.    *
   ********************************************************/

  return (0);
  
}









