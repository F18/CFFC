/*!\file Euler2DQuadSolvers.cc
  \brief 2D Euler Equation Multi-Block Quadrilateral Mesh Solvers. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "Euler2DQuad.h" /* Include 2D Euler quadrilateral mesh solution header file. */
#include "../FASMultigrid2D/FASMultigrid2D.h" /* Include 2D FAS Multigrid solver header file. */
#include "HO_Euler2DQuadGrid.h" /* Include 2D quadrilateral multiblock grid header file for Euler */
#include "../NewtonKrylovSchwarz2D/NKS2D.h" /* Include 2D NKS solver header file. */
#include "Euler2DQuadNKS.h"  /* Inlcude Euler Specializaitons for NKS */


/******************************************************//**
 * Routine: Euler2DQuadSolver                           
 *                                                      
 * Computes solutions to 2D Euler equations on 2D       
 * quadrilateral multi-block solution-adaptive mesh.    
 *                                                      
 ********************************************************/
int Euler2DQuadSolver(char *Input_File_Name_ptr,
                      int batch_flag) {

  /********************************************************  
   * Local variable declarations.                         *
   ********************************************************/

  // Euler2D input variables and parameters:
  Euler2D_Input_Parameters Input_Parameters;

  /* Multi-block solution-adaptive quadrilateral mesh 
     solution variables. */
 
  Grid2D_Quad_MultiBlock_HO    MeshBlk;
  QuadTreeBlock_DataStructure  QuadTree;
  AdaptiveBlockResourceList    List_of_Global_Solution_Blocks;
  AdaptiveBlock2D_List         List_of_Local_Solution_Blocks;
  Euler2D_Quad_Block          *Local_SolnBlk;

  FAS_Multigrid2D_Solver<Euler2D_cState, 
                         Euler2D_Quad_Block, 
                         Euler2D_Input_Parameters> MGSolver;

  /* Define residual file and cpu time variables. */

  ofstream residual_file;
  CPUTime processor_cpu_time, total_cpu_time, NKS_total_cpu_time;
  time_t start_explicit = 0, end_explicit = 0;

  /* Other local solution variables. */

  int number_of_time_steps, first_step,
      command_flag, error_flag, line_number, 
      i_stage, perform_explicit_time_marching, limiter_freezing_off;

  double Time, dTime;

  double residual_l2norm_first, residual_ratio;
  double residual_l1_norm, residual_l2_norm, residual_max_norm;

  /********************************************************  
   * Set default values for the input solution parameters *
   * and then read user specified input values from the   *
   * specified input parameter file.                      *
   ********************************************************/

  // The primary MPI processor processes the input parameter file.
  if (CFFC_Primary_MPI_Processor()) {
     if (!batch_flag) {
        cout << "\n Reading Euler2D input data file `"
             << Input_File_Name_ptr << "'." << endl;
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
   * Create initial mesh and allocate Euler2D solution    *
   * variables for specified IBVP/BVP problem.            *
   ********************************************************/

  execute_new_calculation: ;
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

  /* Create initial mesh.  Read mesh from grid definition or data files 
     when specified by input parameters. */

  if (Input_Parameters.i_ICs != IC_RESTART) {
    // Generate the mesh only if the current run is NOT a restart!

    // The primary MPI processor creates the initial mesh.
    if (CFFC_Primary_MPI_Processor()) {
      if (!batch_flag){ cout << "\n Creating (or reading) initial quadrilateral multi-block mesh."; cout.flush(); }
      error_flag = MeshBlk.Multi_Block_Grid(Input_Parameters);

      if (error_flag) {
        cout << "\n Euler2D ERROR: Unable to create valid Euler2D multi-block mesh.\n";
        cout.flush();
      } else {
        if (Input_Parameters.i_Grid == GRID_NASA_ROTOR_37) {
	  if (!batch_flag) cout << "\n Writing geometry and flow field data files for NASA Rotor 37.";
	  Input_Parameters.NASA_Rotor37.outputTP_UpstreamFlowConditions("NASARotor37_UpstreamProfile.dat");
	  Input_Parameters.NASA_Rotor37.outputTP_DownstreamFlowConditions("NASARotor37_DownstreamProfile.dat");
	  Input_Parameters.NASA_Rotor37.outputTP_Geometry("NASARotor37_Geometry.dat");
	  //             Input_Parameters.NASA_Rotor37.outputTP_Mesh3D(48, 24, 12,
	  //                                                           Input_Parameters.NASA_Rotor37.z_up, 
	  //                                                           Input_Parameters.NASA_Rotor37.z_d,
	  //                                                           "NASARotor37_3DMesh_Top.dat",
	  //                                                           "NASARotor37_3DMesh_Bottom.dat");
	  //             Input_Parameters.NASA_Rotor37.output_Tetin("NASARotor37.tetin", 
	  //                                                        2);
        } else if (Input_Parameters.i_Grid == GRID_NASA_ROTOR_67 ){
	  if (!batch_flag) cout << "\n Writing geometry and flow field data files for NASA Rotor 67.";
	  Input_Parameters.NASA_Rotor67.outputTP_UpstreamFlowConditions("NASARotor67_UpstreamProfile.dat");
	  Input_Parameters.NASA_Rotor67.outputTP_DownstreamFlowConditions("NASARotor67_DownstreamProfile.dat");
	  Input_Parameters.NASA_Rotor67.outputTP_Geometry("NASARotor67_Geometry.dat");
	  Input_Parameters.NASA_Rotor67.outputTP_FlowField2D(Input_Parameters.Rotor_Percent_Span,
							     "NASARotor67_InterbladeFlow.dat");
	  //             Input_Parameters.NASA_Rotor67.outputTP_FlowField2D_All("NASARotor67_InterbladeFlow.dat");
	  //             Input_Parameters.NASA_Rotor67.outputTP_Mesh3D(48, 24, 12,
	  //                                                           Input_Parameters.NASA_Rotor67.z_up, 
	  //                                                           Input_Parameters.NASA_Rotor67.z_d,
	  //                                                           "NASARotor67_3DMesh_Top.dat",
	  //                                                           "NASARotor67_3DMesh_Bottom.dat");
	  //             Input_Parameters.NASA_Rotor67.output_Tetin("NASARotor67.tetin", 
	  //                                                        2);
        } /* endif */
      } /* endif */
    }

    // Broadcast the mesh to other MPI processors.
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
    CFFC_Broadcast_MPI(&error_flag, 1); // Broadcast mesh error flag.
    if (error_flag) return (error_flag);
    MeshBlk.Broadcast_Multi_Block_Grid();

    /* Set the number of blocks in I-dir and J-dir in Input_Parameters based
       on what resulted after the mesh has been created. */
    Input_Parameters.Number_of_Blocks_Idir = MeshBlk.Blocks_Idir();
    Input_Parameters.Number_of_Blocks_Jdir = MeshBlk.Blocks_Jdir();

    /* Create (allocate) multi-block quadtree data structure, create
       (allocate) array of local 2D Euler equation solution blocks, 
       assign and create (allocate) 2D Euler equation solution blocks
       corresponding to the initial mesh. */
    
    if (!batch_flag) cout << "\n Creating multi-block quadtree data structure and assigning"
			  << "\n  Euler2D solution blocks corresponding to initial mesh.";
    Local_SolnBlk = Create_Initial_Solution_Blocks(MeshBlk.Grid_ptr,
						   Local_SolnBlk,
						   Input_Parameters,
						   QuadTree,
						   List_of_Global_Solution_Blocks,
						   List_of_Local_Solution_Blocks);

  } else {
    // Allocate the minimum information related to the solution blocks. (i.e. use the default constructors)
    Local_SolnBlk = Allocate(Local_SolnBlk,Input_Parameters);
  }
  if (Local_SolnBlk == NULL) return (1);

  if (Input_Parameters.i_ICs != IC_RESTART) {
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

  /********************************************************  
   * Initialize Euler2D solution variables.               *
   ********************************************************/

  /* Set the initial time level. */

  Time = ZERO;
  number_of_time_steps = 0;
  
  /* Set the CPU time to zero. */
  
  processor_cpu_time.zero();
  total_cpu_time.zero();
  NKS_total_cpu_time.zero();
  
  /* Initialize the conserved and primitive state
     solution variables. */
  if (Input_Parameters.i_ICs == IC_RESTART) {
    if (!batch_flag) { cout << "\n Reading Euler2D solution from restart data files."; cout.flush(); }

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
        cout << "\n Euler2D ERROR: Unable to open Euler2D quadtree data file "
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
        cout << "\n Euler2D ERROR: Unable to open Euler2D restart input data file(s) "
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
    if (!batch_flag){ cout << "\n Prescribing Euler2D initial data."; cout.flush(); }
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
                                                  NUM_VAR_EULER2D,
                                                  OFF);
  if (error_flag) {
     cout << "\n Euler2D ERROR: Message passing error during Euler2D solution intialization "
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

  /* Perform uniform, boundary, and, initial mesh refinement. */

  if (Input_Parameters.i_ICs != IC_RESTART) {
     if (!batch_flag) cout << "\n Performing Euler2D uniform mesh refinement.";
     error_flag = Uniform_AMR(Local_SolnBlk,
                              Input_Parameters,
                              QuadTree,
                              List_of_Global_Solution_Blocks,
                              List_of_Local_Solution_Blocks);
     if (error_flag) {
       cout << "\n Euler2D ERROR: Uniform AMR error on processor "
	    << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
       cout.flush();
     } /* endif */
     error_flag = CFFC_OR_MPI(error_flag);
     if (error_flag) return error_flag;

     if (!batch_flag) cout << "\n Performing Euler2D boundary mesh refinement.";
     error_flag = Boundary_AMR(Local_SolnBlk,
                               Input_Parameters,
                               QuadTree,
                               List_of_Global_Solution_Blocks,
                               List_of_Local_Solution_Blocks);
     if (error_flag) {
       cout << "\n Euler2D ERROR: Boundary AMR error on processor "
	    << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
       cout.flush();
     } /* endif */
     error_flag = CFFC_OR_MPI(error_flag);
     if (error_flag) return error_flag;

     if (!batch_flag){ cout << "\n Performing Euler2D initial mesh refinement."; cout.flush(); }
     error_flag = Initial_AMR(Local_SolnBlk,
                              Input_Parameters,
	   	              QuadTree,
		              List_of_Global_Solution_Blocks,
		              List_of_Local_Solution_Blocks);
     if (error_flag) {
        cout << "\n Euler2D ERROR: Initial AMR error on processor "
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

  /* Re-order the solution blocks on the processors by applying 
     the Morton ordering space filling curve algorithm. */
  
  if (Input_Parameters.Morton) {
     if (!batch_flag) cout << "\n\n Applying Morton re-ordering algorithm to initial solution blocks.";
     error_flag = Morton_ReOrdering_of_Solution_Blocks(QuadTree, 
                                                       List_of_Global_Solution_Blocks, 
                                                       List_of_Local_Solution_Blocks, 
                                                       Local_SolnBlk, 
                                                       Input_Parameters, 
                                                       number_of_time_steps, 
                                                       Time, 
                                                       processor_cpu_time); 
     if (error_flag) {
        cout <<"\n Euler2D ERROR: Morton re-ordering error on processor "
	           << List_of_Local_Solution_Blocks.ThisCPU
	           << ".\n";
	      cout.flush();
	      return (error_flag);
     } /* endif */
     error_flag = CFFC_OR_MPI(error_flag);
     if (error_flag) return (error_flag);
     //Output space filling curve in Tecplot format
     if (!batch_flag) {
        cout << "\n Outputting space filling curve showing block loading for CPUs.\n";
        cout.flush();
     } /* endif */
     Morton_SFC_Output_Tecplot(Local_SolnBlk, 
                               Input_Parameters, 
                               List_of_Local_Solution_Blocks);
  } /* endif */

  /********************************************************  
   * Solve IBVP or BVP for conservation form of 2D Euler  *
   * equations on multi-block solution-adaptive           *
   * quadrilateral mesh.                                  *
   ********************************************************/

  continue_existing_calculation: ;
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

  if(!batch_flag) { time(&start_explicit); }

  if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {

    // Allocate memory for multigrid solver.
    error_flag = MGSolver.allocate(Local_SolnBlk,
				   &QuadTree,
				   &List_of_Global_Solution_Blocks,
				   &List_of_Local_Solution_Blocks,
				   &Input_Parameters);
    if (error_flag) {
      cout << "\n Euler2D ERROR: Unable to allocate memory for multigrid solver.\n";
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
      cout << "\n Euler2D ERROR: Error during multigrid solution.\n";
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
        cout << "\n Euler2D ERROR: Unable to open residual file for Euler2D calculation.\n";
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
      if (!batch_flag){ cout << "\n\n Beginning Euler2D computations on "
			     << Date_And_Time() << ".\n\n"; start_explicit = clock();  }

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
                 cout << "\n Euler2D ERROR: Euler2D AMR error on processor "
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

	/* Periodically re-order the solution blocks on the processors
           by applying the Morton ordering space filling curve algorithm. */
	if (Input_Parameters.Morton &&
            !first_step &&
            number_of_time_steps%Input_Parameters.Morton_Reordering_Frequency == 0) {
	  if (!batch_flag) { cout << "\n\n Applying Morton re-ordering algorithm to solution blocks at n = "
				  << number_of_time_steps << "."; cout.flush(); }
	   error_flag = Morton_ReOrdering_of_Solution_Blocks(QuadTree, 
                                                             List_of_Global_Solution_Blocks, 
                                                             List_of_Local_Solution_Blocks, 
                                                             Local_SolnBlk, 
                                                             Input_Parameters, 
                                                             number_of_time_steps, 
                                                             Time, 
                                                             processor_cpu_time); 
	   if (error_flag) {
              cout <<"\n Euler2D ERROR: Morton re-ordering error on processor "
	           << List_of_Local_Solution_Blocks.ThisCPU
	           << ".\n";
	      cout.flush();
 	      return (error_flag);
           } /* endif */
           error_flag = CFFC_OR_MPI(error_flag);
           if (error_flag) return (error_flag);
           //Output space filling curve in Tecplot format
           if (!batch_flag) cout << "\n Outputting space filling curve showing block loading for CPUs.";
	   Morton_SFC_Output_Tecplot(Local_SolnBlk, 
                                     Input_Parameters, 
                                     List_of_Local_Solution_Blocks);
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

	  if (!batch_flag){ cout << "\n\n  Saving Euler2D solution to restart data file(s) after"
				 << " n = " << number_of_time_steps << " steps (iterations)."; cout.flush(); }
          error_flag = Write_QuadTree(QuadTree,
                                      Input_Parameters);
          if (error_flag) {
             cout << "\n Euler2D ERROR: Unable to open Euler2D quadtree data file "
                  << "on processor "
                  << List_of_Local_Solution_Blocks.ThisCPU
                  << ".\n";
             cout.flush();
          } /* endif */
          error_flag = CFFC_OR_MPI(error_flag);
          if (error_flag) return (error_flag);

	  error_flag = Write_Restart_Solution(Local_SolnBlk, 
					      List_of_Local_Solution_Blocks, 
					      Input_Parameters,
					      number_of_time_steps,
					      Time,
					      processor_cpu_time);
	  if (error_flag) {
	    cout << "\n Euler2D ERROR: Unable to open Euler2D restart output data file(s) "
		 << "on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  } /* endif */
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return (error_flag);
	  if (!batch_flag) {
	    cout << "\n";
	    cout.flush();
	  }

	  if (CFFC_Primary_MPI_Processor()) {
	    System::Remove_Restart_Flag();  //Remove flag to indicate the restart is finished
	  }

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
					 NUM_VAR_EULER2D, 
					 OFF);
	  if (error_flag) {
	    cout << "\n Euler2D ERROR: Euler2D message passing error on processor "
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
	    cout << "\n Euler2D ERROR: Euler2D solution residual error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  } /* endif */
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return (error_flag);

	  // 4. Send boundary flux corrections at block interfaces with resolution changes.
	  error_flag = Send_Conservative_Flux_Corrections(Local_SolnBlk, 
							  List_of_Local_Solution_Blocks,
							  NUM_VAR_EULER2D);
	  if (error_flag) {
	    cout << "\n Euler2D ERROR: Euler2D flux correction message passing error on processor "
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

          // 6. Smooth the solution residual using implicit residual smoothing.
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
	    cout << "\n Euler2D ERROR: Euler2D solution update error on processor "
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

      if (!batch_flag){ cout << "\n\n Euler2D computations complete on " 
			     << Date_And_Time() << ".\n"; end_explicit = clock(); }

    } /* endif */
    
    /* Update ghostcell information and prescribe boundary conditions to ensure
       that the solution is consistent on each block. */
    
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
    
    error_flag = Send_All_Messages(Local_SolnBlk, 
				   List_of_Local_Solution_Blocks,
				   NUM_VAR_EULER2D,
				   OFF);
    if (error_flag) {
      cout << "\n Euler2D ERROR: Euler2D message passing error on processor "
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
           cout << "\n Euler2D ERROR: Unable to open residual file for Euler2D calculation.\n";
           cout.flush();
        } /* endif */ 
     } /* endif */

     CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
     CFFC_Broadcast_MPI(&error_flag, 1);
     if (error_flag) return (error_flag);

     //Turn Limiter Freezing OFF for startup 
     Evaluate_Limiters(Local_SolnBlk, List_of_Local_Solution_Blocks);

     if (!batch_flag){ cout << "\n\n Beginning Euler2D NKS computations on " << Date_And_Time() << ".\n\n"; time(&start_NKS); }

     //Store Explicit times for output
     CPUTime Explicit_processor_cpu_time = processor_cpu_time;
     CPUTime Explicit_total_cpu_time =  total_cpu_time;
    
     //Perform NKS Iterations 
     error_flag = Newton_Krylov_Schwarz_Solver<Euler2D_pState,                    //pass in Time for DTS ??
                                               Euler2D_Quad_Block,                                               
                                               Euler2D_Input_Parameters>(processor_cpu_time,
									 residual_file,
									 number_of_time_steps, // explicit time steps
									 Time,	
									 Local_SolnBlk, 
									 QuadTree,
									 List_of_Global_Solution_Blocks,	
						       			 List_of_Local_Solution_Blocks,
									 Input_Parameters);
      
     processor_cpu_time.update();
     total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput); 

     if (error_flag) {
        if (CFFC_Primary_MPI_Processor()) { 
   	   cout << "\n Euler2D_NKS ERROR: Euler2D solution error on processor " 
                << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
   	   cout.flush();
   	} /* endif */
     } /* endif */
   
     CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
     CFFC_Broadcast_MPI(&error_flag, 1);
     if (error_flag) return (error_flag);

     /***********************************************************************/
     if (!batch_flag) { cout << "\n\n Euler2D NKS computations complete on " << Date_And_Time() << ".\n"; time(&end_NKS); }

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

  /********************************************************
   * Solution calculations complete.                      *
   * Write 2D Euler solution to output and restart files  *
   * as required, reset solution parameters, and run      *
   * other cases as specified by input parameters.        *
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
      // Deallocate memory for 2D Euler equation solution.
      if (!batch_flag) cout << "\n Deallocating Euler2D solution variables.";
      if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {
	MGSolver.deallocate();
      }
      Local_SolnBlk = Deallocate(Local_SolnBlk, 
				 Input_Parameters);
      //Deallocate_Message_Buffers(List_of_Local_Solution_Blocks); // Not necessary here!
      List_of_Local_Solution_Blocks.deallocate();
      List_of_Global_Solution_Blocks.deallocate();
      QuadTree.deallocate();
      // Output input parameters for new caluculation.
      if (!batch_flag)  {
	cout << "\n\n Starting a new calculation.";
	cout << Input_Parameters << "\n";
	cout.flush();
      } /* endif */
      // Execute new calculation.
      goto execute_new_calculation;
      
    } else if (command_flag == TERMINATE_CODE) {
      // Deallocate memory for 2D Euler equation solution.
      if (!batch_flag) cout << "\n Deallocating Euler2D solution variables.";
      if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	  Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
	MGSolver.deallocate();
      }
      Local_SolnBlk = Deallocate(Local_SolnBlk, 
				 Input_Parameters);
      //Deallocate_Message_Buffers(List_of_Local_Solution_Blocks); // Not necessary here!
      List_of_Local_Solution_Blocks.deallocate();
      List_of_Global_Solution_Blocks.deallocate();
      QuadTree.deallocate();
      // Close input data file.
      if (!batch_flag) cout << "\n\n Closing Euler2D input data file.";
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
	cout << "\n Euler2D ERROR: Euler2D AMR error on processor "
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

   } else if (command_flag == MORTON_ORDERING_CODE) {
      if (!batch_flag) cout << "\n\n Applying Morton re-ordering algorithm.";
        error_flag = Morton_ReOrdering_of_Solution_Blocks(QuadTree, 
                                                          List_of_Global_Solution_Blocks, 
                                                          List_of_Local_Solution_Blocks, 
                                                          Local_SolnBlk, 
                                                          Input_Parameters, 
                                                          number_of_time_steps, 
                                                          Time, 
                                                          processor_cpu_time); 
        if (error_flag) {
           cout <<"\n Euler2D ERROR: Morton re-ordering error on processor "
                << List_of_Local_Solution_Blocks.ThisCPU
                << ".\n";
           cout.flush();
           return (error_flag);
        } /* endif */
        error_flag = CFFC_OR_MPI(error_flag);
        if (error_flag) return (error_flag);
        //Output space filling curve in Tecplot format
        if (!batch_flag) cout << "\n Outputting space filling curve showing block loading for CPUs.";
        Morton_SFC_Output_Tecplot(Local_SolnBlk, 
                                  Input_Parameters,
                                  List_of_Local_Solution_Blocks);

    } else if (command_flag == WRITE_OUTPUT_CODE) {
       // Output solution data.
       if (!batch_flag) cout << "\n Writing Euler2D solution to output data file(s).";

       if (Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0 ||
           !(Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
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
          cout << "\n Euler2D ERROR: Unable to open Euler2D output data file(s) "
               << "on processor "
               << List_of_Local_Solution_Blocks.ThisCPU
               << ".\n";
          cout.flush();
       } /* endif */
       error_flag = CFFC_OR_MPI(error_flag);
       if (error_flag) return (error_flag);

   } else if (command_flag == WRITE_OUTPUT_CELLS_CODE) {
       // Output solution data.
       if (!batch_flag) cout << "\n Writing cell-centered Euler2D solution to output data file(s).";
       if (Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0 ||
	   !(Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
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
          cout << "\n Euler2D ERROR: Unable to open Euler2D cell output data file(s) "
               << "on processor "
               << List_of_Local_Solution_Blocks.ThisCPU
               << ".\n";
          cout.flush();
       } /* endif */
       error_flag = CFFC_OR_MPI(error_flag);
       if (error_flag) return (error_flag);

    } else if (command_flag == WRITE_OUTPUT_NODES_CODE) {
      if (!batch_flag) cout << "\n Writing Euler2D node locations to output data file(s).";
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
	cout << "\n Euler2D ERROR: Unable to open Euler2D nodes output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      } /* endif */
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_GRADIENTS_CODE) {
      if (!batch_flag) cout << "\n Writing Euler2D primitive state gradients to output data file(s).";
      error_flag = Output_Gradients_Tecplot(Local_SolnBlk,
					    List_of_Local_Solution_Blocks,
					    Input_Parameters,
					    number_of_time_steps,
					    Time);
      if (error_flag) {
	cout << "\n Euler2D ERROR: Unable to open Euler2D gradient output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      } /* endif */
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_QUASI3D_CODE) {
       // Output solution data.
       if (!batch_flag) cout << "\n Writing Euler2D quasi3D solution to output data file(s).";
       error_flag = Output_Tecplot_Quasi3D(Local_SolnBlk, 
					   List_of_Local_Solution_Blocks, 
					   Input_Parameters,
					   number_of_time_steps,
					   Time);
       if (error_flag) {
          cout << "\n Euler2D ERROR: Unable to open Euler2D quasi3D output data file(s) "
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

       if (!batch_flag) cout << "\n Writing Euler2D solution to restart data file(s).";
       error_flag = Write_QuadTree(QuadTree,
                                   Input_Parameters);
       if (error_flag) {
          cout << "\n Euler2D ERROR: Unable to open Euler2D quadtree data file "
               << "on processor "
               << List_of_Local_Solution_Blocks.ThisCPU
               << ".\n";
          cout.flush();
       } /* endif */
       error_flag = CFFC_OR_MPI(error_flag);
       if (error_flag) return (error_flag);
       error_flag = Write_Restart_Solution(Local_SolnBlk, 
                                           List_of_Local_Solution_Blocks, 
                                           Input_Parameters,
                                           number_of_time_steps,
                                           Time,
                                           processor_cpu_time);
       if (error_flag) {
          cout << "\n Euler2D ERROR: Unable to open Euler2D restart output data file(s) "
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
          if (!batch_flag) cout << "\n Writing Euler2D multi-block mesh to grid data output file.";
          error_flag = MeshBlk.Output_Tecplot_Using_IP(Input_Parameters);
          if (error_flag) {
             cout << "\n Euler2D ERROR: Unable to open Euler2D mesh data output file.\n";
             cout.flush();
          } /* endif */
       } /* endif */
       CFFC_Broadcast_MPI(&error_flag, 1);
       if (error_flag) return (error_flag);

    } else if (command_flag == WRITE_GRID_DEFINITION_CODE) {
       // Write multi-block solution-adaptive mesh definition files.
       if (CFFC_Primary_MPI_Processor()) {
          if (!batch_flag) cout << "\n Writing Euler2D multi-block mesh to grid definition files.";
          error_flag = MeshBlk.Write_Multi_Block_Grid_Definition_Using_IP(Input_Parameters);
          if (error_flag) {
             cout << "\n Euler2D ERROR: Unable to open Euler2D multi-block mesh definition files.\n";
             cout.flush();
          } /* endif */
       } /* endif */
       CFFC_Broadcast_MPI(&error_flag, 1);
       if (error_flag) return (error_flag);

    } else if (command_flag == WRITE_OUTPUT_GRID_NODES_CODE) {
       // Output multi-block solution-adaptive mesh node data file.
       if (CFFC_Primary_MPI_Processor()) {
          if (!batch_flag) cout << "\n Writing Euler2D multi-block mesh to node data output file.";
          error_flag = MeshBlk.Output_Nodes_Tecplot_Using_IP(Input_Parameters);
          if (error_flag) {
             cout << "\n Euler2D ERROR: Unable to open Euler2D mesh node data output file.\n";
             cout.flush();
          } /* endif */
       } /* endif */
       CFFC_Broadcast_MPI(&error_flag, 1);
       if (error_flag) return (error_flag);

    } else if (command_flag == WRITE_OUTPUT_GRID_CELLS_CODE) {
       // Output multi-block solution-adaptive mesh cell data file.
       if (CFFC_Primary_MPI_Processor()) {
          if (!batch_flag) cout << "\n Writing Euler2D multi-block mesh to cell data output file.";
          error_flag = MeshBlk.Output_Cells_Tecplot_Using_IP(Input_Parameters);
          if (error_flag) {
             cout << "\n Euler2D ERROR: Unable to open Euler2D mesh cell data output file.\n";
             cout.flush();
          } /* endif */
       } /* endif */
       CFFC_Broadcast_MPI(&error_flag, 1);
       if (error_flag) return (error_flag);

    } else if (command_flag == WRITE_OUTPUT_RINGLEB_CODE) {
       // Output Ringleb flow solution information.
       if (CFFC_Primary_MPI_Processor()) {
          if (!batch_flag) cout << "\n Writing Euler2D Ringleb flow solution data.";
          error_flag = Output_Ringleb_Flow(Local_SolnBlk, 
					   List_of_Local_Solution_Blocks, 
                                           Input_Parameters);
          if (error_flag) {
             cout << "\n Euler2D ERROR: Unable to write Euler2D Ringleb flow data.\n";
             cout.flush();
          } /* endif */
       } /* endif */
       CFFC_Broadcast_MPI(&error_flag, 1);
       if (error_flag) return (error_flag);

    } else if (command_flag == DETERMINE_MACH_STEM_HEIGHT_CODE) {
      if (!batch_flag) cout << "\n Determine the height of the Mach stem.";
      double maximum_pressure;
      error_flag = Determine_Maximum_Surface_Pressure(Local_SolnBlk,
						      List_of_Local_Solution_Blocks,
						      Input_Parameters,
						      maximum_pressure);
      if (error_flag) {
	cout << "\n Euler2D ERROR: Unable to determine the maximum pressure on the surface of the wedge.\n";
	cout.flush();
      } /* endif */
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
      error_flag = Determine_Mach_Stem_Height(Local_SolnBlk,
					      List_of_Local_Solution_Blocks,
					      Input_Parameters,
					      HALF*(maximum_pressure + 101325.0),
					      number_of_time_steps,
					      Time);
      if (error_flag) {
	cout << "\n Euler2D ERROR: Unable to determine the Mach stem height.\n";
	cout.flush();
      } /* endif */
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);

    } else if (command_flag == WRITE_OUTPUT_WEDGE_SOLUTION_DISTRIBUTION_CODE) {
      if (!batch_flag) cout << "\n Writing Euler2D solution on the wedge.";
      error_flag = Output_Wedge_Solution_Distribution_Tecplot(Local_SolnBlk,
							      List_of_Local_Solution_Blocks,
							      Input_Parameters,
							      number_of_time_steps,
							      Time);
      if (error_flag) {
	cout << "\n Euler2D ERROR: Unable to write the solution on the wedge.\n";
	cout.flush();
      } /* endif */
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);

    } else if (command_flag == INVALID_INPUT_CODE ||
               command_flag == INVALID_INPUT_VALUE) {
        line_number = -line_number;
        cout << "\n Euler2D ERROR: Error reading Euler2D data at line #"
             << -line_number  << " of input data file.\n";
        cout.flush();
        return (line_number);
    } /* endif */

  } /* endwhile */

  /********************************************************  
   * End of all Euler2DSolver computations and I/O.       *
   ********************************************************/

  return (0);
  
}
