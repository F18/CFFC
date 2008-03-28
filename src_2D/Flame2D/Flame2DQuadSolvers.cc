/////////////////////////////////////////////////////////////////////
///
/// \file Flame2DQuadSolvers.cc
/// 
/// \author Marc R.J. Charest
/// 
/// \brief 2D Axisymmetric Solver for....
///
/////////////////////////////////////////////////////////////////////
// Include main solution block header file.
#include "Flame2DQuad.h"

// Include some tools
#include "Flame2DTools.h"

// Include FLAME2D NKS Specializations header file.
#include "Flame2DQuadNKS.h"

// Include Rte2D solver Sepcialization header file.
#include "Flame2DQuadRte.h"


/////////////////////////////////////////////////////////////////////
///
/// Routine: Flame2DQuadSolver
///
/// Computes solutions to 2D .....  equations on 2D
/// quadrilateral multi-block solution-adaptive mesh.
///
/////////////////////////////////////////////////////////////////////
int Flame2DQuadSolver(char *Input_File_Name_ptr,  int batch_flag) {
  
  /********************************************************  
   * Local variable declarations.                         *
   ********************************************************/
  
  // Flame2D input variables and parameters:
  Flame2D_Input_Parameters Input_Parameters;
  
  // Multi-block solution-adaptive quadrilateral mesh 
  // solution variables.
  Grid2D_Quad_Block          **MeshBlk;
  QuadTreeBlock_DataStructure  QuadTree;
  AdaptiveBlockResourceList    List_of_Global_Solution_Blocks;
  AdaptiveBlock2D_List         List_of_Local_Solution_Blocks;
  Flame2D_Quad_Block           *Local_SolnBlk;
    
  // Radiation solver object pointer and parameters
#ifdef _FLAME2D_WITH_RTE
  Rte2DSolver *RteSolver(NULL);     // object pointer
#endif
  bool Rte_PostProcess(false);      // don't post process radiation solution after solve
  int Rte_batch_flag(1);            // don't print out information
  int number_sequential_solves(0);  // current sequential solve number

  // Define residual file and cpu time variables.
  ofstream residual_file;
  ofstream time_accurate_data_file;
  CPUTime processor_cpu_time, total_cpu_time;
  time_t start_explicit, end_explicit;

  // Other local solution variables. 
  int number_of_time_steps, first_step, last_step,
    command_flag, error_flag, line_number, 
    i_stage, limiter_freezing_off;
  int Max_Number_Global_Steps;
  double Time, dTime;


  /*************************************************************************
   ******************** INPUT PARAMETERS  **********************************
     Set default values for the input solution parameters and then read user 
     specified input values from the specified input parameter file.            
     *************************************************************************
     *************************************************************************/  
  //The primary MPI processor processes the input parameter file.
  if (CFFC_Primary_MPI_Processor()) {
    if (!batch_flag) {
      cout << "\n Reading Flame2D input data file `"
	   << Input_File_Name_ptr << "'."; 
    } // endif
    error_flag = Process_Input_Control_Parameter_File(Input_Parameters,
						      Input_File_Name_ptr,
						      command_flag);
    if (!batch_flag && error_flag == 0) {
      cout << Input_Parameters << "\n";
      cout.flush();
    } // endif
  } else {
    error_flag = 0;
  } // endif


    // Broadcast input solution parameters to other MPI processors.
    // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();
  CFFC_Broadcast_MPI(&error_flag, 1);
  if (error_flag != 0) return (error_flag); 
  CFFC_Broadcast_MPI(&command_flag, 1);
  if (command_flag == TERMINATE_CODE) return (0);

  Broadcast_Input_Parameters(Input_Parameters);
 
  // Common Use Variables
  const int NUM_VAR_FLAME2D = Input_Parameters.Wo.NumVar();
  
  /*************************************************************************
   ******************** INITIAL GRID & SOLUTION SPACE **********************
     Create initial mesh and allocate Flame2D solution variables for 
     specified IBVP/BVP problem. 
     *************************************************************************
     *************************************************************************/

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 execute_new_calculation: ;
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();
  
  // Create initial mesh.  Read mesh from grid definition or data files 
  // when specified by input parameters.
  
  // The primary MPI processor creates the initial mesh.
  if (CFFC_Primary_MPI_Processor()) {
    if (!batch_flag) cout << "\n Creating (or reading) initial quadrilateral multi-block mesh.";
    MeshBlk = NULL;
    MeshBlk = Multi_Block_Grid(MeshBlk, Input_Parameters);
    
    if (MeshBlk == NULL) {
      error_flag = 1;
    } else if (Check_Multi_Block_Grid(MeshBlk,
				      Input_Parameters.Number_of_Blocks_Idir,
				      Input_Parameters.Number_of_Blocks_Jdir)) {
      error_flag = 1;
    } else {
      error_flag = 0;
    } // endif
    
    if (error_flag) {
      cout << "\n Flame2D ERROR: Unable to create valid Flame2D multi-block mesh.\n";
      cout.flush();
    } // endif
  } else {
    MeshBlk = NULL;
  } // endif

    // Broadcast the mesh to other MPI processors.
    // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();
  // Broadcast mesh error flag.
  CFFC_Broadcast_MPI(&error_flag, 1); 

  if (error_flag) return (error_flag);
  MeshBlk = Broadcast_Multi_Block_Grid(MeshBlk, Input_Parameters);

  /************************************************************************
   ******************** MULTIBLOCK QUADTREE *******************************
       Create (allocate) multi-block quadtree data structure, create
       (allocate) array of local  2D equation solution blocks, 
       assign and create (allocate) 2D equation solution blocks
       corresponding to the initial mesh. 
       ************************************************************************
       ************************************************************************/
  if (!batch_flag){
    cout << "\n Creating multi-block quadtree data structure and assigning"
	 << "\n Flame2D solution blocks corresponding to initial mesh.";
  }
  //needs to be run by each processor
  Local_SolnBlk = Create_Initial_Solution_Blocks(MeshBlk,
						 Local_SolnBlk,
						 Input_Parameters,
						 QuadTree,
						 List_of_Global_Solution_Blocks,
						 List_of_Local_Solution_Blocks);

 
  
  if (Local_SolnBlk == NULL){ return (1); }

  /************************************************************************  
   * Initialize Flame2D solution variables, Initial conditions             *  
   ************************************************************************/

  // Set the initial time level.
  Time = ZERO;
  number_of_time_steps = 0;
  Max_Number_Global_Steps = 0;

  // Set the CPU time to zero.
  processor_cpu_time.zero();
  total_cpu_time.zero(); 

  // Initialize the conserved and primitive state solution variables.
  if (!batch_flag){ 
    cout << "\n Prescribing Flame2D initial data.";
  }
 
  /*************************************************************************
   **************** RESTART SOLUTION or INITIAL CONDITIONS *****************
   *************************************************************************/
  //
  // Read restart data
  //
  if (Input_Parameters.i_ICs == IC_RESTART) {
    if (!batch_flag) cout << "\n Reading Flame2D solution from restart data files.";
    error_flag = Read_QuadTree(QuadTree,
			       List_of_Global_Solution_Blocks,
			       List_of_Local_Solution_Blocks, 
			       Input_Parameters);
    if (error_flag) {
      cout << "\n Flame2D ERROR: Unable to open Flame2D quadtree data file "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } // endif
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
      cout << "\n Flame2D ERROR: Unable to open Flame2D restart input data file(s) "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } // endif
    error_flag = CFFC_OR_MPI(error_flag); 
    if (error_flag){ 
      return (error_flag);
    } 
    // Ensure each processor has the correct time and time!!!
    number_of_time_steps = CFFC_Maximum_MPI(number_of_time_steps);
    Time = CFFC_Maximum_MPI(Time);
    processor_cpu_time.cput = CFFC_Maximum_MPI(processor_cpu_time.cput);
    Input_Parameters.Maximum_Number_of_Time_Steps = CFFC_Maximum_MPI(Input_Parameters.Maximum_Number_of_Time_Steps);
    Input_Parameters.Time_Max = CFFC_Maximum_MPI(Input_Parameters.Time_Max);

    // Send grid information between neighbouring blocks
    CFFC_Barrier_MPI(); 
    error_flag = Send_All_Messages(Local_SolnBlk, 
				   List_of_Local_Solution_Blocks,
				   NUM_COMP_VECTOR2D,
				   ON);     

    //
    // Else apply initial conditions from input parameters
    //
  } else {   

    // Send grid information between neighbouring blocks BEFORE applying ICs
    CFFC_Barrier_MPI(); 
    error_flag = Send_All_Messages(Local_SolnBlk, 
				   List_of_Local_Solution_Blocks,
				   NUM_COMP_VECTOR2D,
				   ON);     

    ICs(Local_SolnBlk, List_of_Local_Solution_Blocks, Input_Parameters);
  } // endif
  
  
    /******************************************************************************  
      Send solution information between neighbouring blocks to complete
      prescription of initial data. 
    *******************************************************************************/
    // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();
 
  // OLD grid send_all_messages, moved before ICs to help with values used in BCs 
  // set from ICs, ie ghost cell node locations, changed before ICs 

  if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk, 
						  List_of_Local_Solution_Blocks,
						  NUM_VAR_FLAME2D,
						  OFF);
  if (error_flag) {
    cout << "\n Flame2D ERROR: Message passing error during Flame2D solution intialization "
	 << "on processor "
	 << List_of_Local_Solution_Blocks.ThisCPU
	 << ".\n";
    cout.flush();
    exit(1);
  } // endif 

  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return (error_flag);
    
  /*******************************************************************************
   ************************* BOUNDARY CONDITIONS *********************************
   *******************************************************************************/
  // Prescribe boundary data consistent with initial data. 
  BCs(Local_SolnBlk, List_of_Local_Solution_Blocks, Input_Parameters);

  /*******************************************************************************
   ******************* ADAPTIVE MESH REFINEMENT (AMR) ****************************
   *******************************************************************************/  
  // Perform uniform, boundary, and, initial mesh refinement.
  if (Input_Parameters.i_ICs != IC_RESTART) {
    if (!batch_flag) cout << "\n Performing Flame2D uniform mesh refinement.";
    error_flag = Uniform_AMR(Local_SolnBlk,
			     Input_Parameters,
			     QuadTree,
			     List_of_Global_Solution_Blocks,
			     List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n Flame2D ERROR: Uniform AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
      cout.flush();
    } // endif
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    ///////////////////////////////////////////////////////////////////////////
    if (!batch_flag) cout << "\n Performing Flame2D boundary mesh refinement.";
    error_flag = Boundary_AMR(Local_SolnBlk,
			      Input_Parameters,
			      QuadTree,
			      List_of_Global_Solution_Blocks,
			      List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n Flame2D ERROR: Boundary AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
      cout.flush();
    } // endif
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    ///////////////////////////////////////////////////////////////////////////
    if (!batch_flag) cout << "\n Performing Flame2D initial mesh refinement.";
    error_flag = Initial_AMR(Local_SolnBlk,
			     Input_Parameters,
			     QuadTree,
			     List_of_Global_Solution_Blocks,
			     List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n Flame2D ERROR: Initial AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } // endif
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
  } // endif

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
	 << QuadTree.efficiencyRefinement() << "\n"; 
    cout.flush();
  } // endif
 
    /***********************************************************************	
      MORTON ORDERING of initial solution blocks 
      (should be meshed with AMR, ie when Refine_Grid is done call the ordering)
    ************************************************************************/
  if (Input_Parameters.Morton){
    if (!batch_flag) cout << "\n\n Applying Morton re-ordering algorithm to initial solution blocks. ";

    //NOTES: Issue here with Input_Parameters.Maximum_Number_of_Time_Steps related to reading restart files

    error_flag = Morton_ReOrdering_of_Solution_Blocks(QuadTree, 
						      List_of_Global_Solution_Blocks, 
						      List_of_Local_Solution_Blocks, 
						      Local_SolnBlk, 
						      Input_Parameters, 
						      number_of_time_steps, 
						      Time, 
						      processor_cpu_time); 
    if (error_flag) {
      cout <<"\n Flame2D ERROR: Morton re-ordering error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
      return (error_flag);
    } 
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    //Output space filling curve in Tecplot format
    if (!batch_flag) cout << "\n Outputting space filling curve showing block loading for CPUs.";
    Morton_SFC_Output_Tecplot(Local_SolnBlk, 
			      Input_Parameters, 
			      List_of_Local_Solution_Blocks);
  } 
  
  /***********************************************************************	
      RADIATION SOLVER setup of initial grid, solution blocks, quadtree,
      initial conditions.                                                      
  ************************************************************************/
  if (Input_Parameters.Radiation == RADIATION_RTE){

#ifdef _FLAME2D_WITH_RTE

    if (!batch_flag) cout << "\n Setting up radiation solver.";

    // deallocate radiation solver just in case
    if (RteSolver != NULL) delete RteSolver;

    // allocate a new radiation object
    RteSolver = new Rte2DSolver;

    // Pass the chem2D input parameters, quadtree data structure, and mesh
    // to initialize the radiation solver.  Note that the radiation solver
    // uses a copy of the same grid and quadtree.  All AMR is performed first
    // on the Flame2D mesh (by Flame2D), and the resulting refinement flags 
    // are passed to the radation solver and applied.  This keeps the meshes 
    // in sync.
    RteSolver->SetupSequentialSolve(Rte_batch_flag, 
				    Input_Parameters.Rte_Input_File_Name,
				    MeshBlk, 
				    Input_Parameters);

    // We need to initialize the radiation field. There are several cases 
    // here:
    //  1 -> Flame2D is being restarted but Rte2D is not (ie. changing 
    //       radiation parameters or turning on radiation for first time).
    //       |-> Read Flame2D solution from restart, compute initial Rte2D solution,
    //           variables to Flame2D
    //  2 -> Flame2D and Rte2D are both being restarted.
    //       |-> Read both solutions from restarts and copy over Rte2D solution
    //           variables to Flame2D
    //  3 -> Completely new simulation (Neither are being restarted).
    //       |-> No initialization required
    if (Input_Parameters.i_ICs == IC_RESTART && 
	RteSolver->Input_Parameters.i_ICs == IC_RESTART) {

      // Copy over computed Rte2D solution variables only
      //RteSolver->Copy_SRC_Solution_Vars(Local_SolnBlk); // set medium state using Flame2D data
      RteSolver->Copy_Rte2D_Solution_Vars(Local_SolnBlk); // compute Flame2D rad source using Rte2D data

    } // endif - restart

      // set the current number of sequential solves and the update frequency
    number_sequential_solves = 0;

#else

    cerr << " \nCoupled Flame2D / Rte2D solver requested."
	 << " Recompile code with _FLAME2D_WITH_RTE #define.\n";
    exit(-1);

#endif //_FLAME2D_WITH_RTE

  } // endif

    /****************************************************************************
     *********************** MAIN SOLVER ****************************************
     Solve IBVP or BVP for conservation form of 2D Axisymmetric multispecies 
     chemically reacting thermally perfect equations on multi-block 
     solution-adaptive quadrilateral mesh.                                  
     ****************************************************************************
     ****************************************************************************/  
  
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 continue_existing_calculation: ;
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI(); 
 

  // print out a header for sequential solver
  if (!batch_flag && Input_Parameters.Radiation == RADIATION_RTE &&
      Input_Parameters.Max_Number_Sequential_Solves > 0) {
    cout << endl << endl << string(70,'=');
    cout << "\n Sequential Solver Step = " << number_sequential_solves+1;
    cout << endl << string(70,'=');
  }

  /////////////////////////////////////////////////////////////////////////////
  /// NON-MULTIGRID EXPLICT
  /////////////////////////////////////////////////////////////////////////////

  if(Input_Parameters.Maximum_Number_of_Time_Steps > 0){ 
    double *residual_l1_norm = new double[Local_SolnBlk[0].Number_of_Residual_Norms]; 
    double *residual_l2_norm = new double[Local_SolnBlk[0].Number_of_Residual_Norms]; 
    double *residual_max_norm = new double[Local_SolnBlk[0].Number_of_Residual_Norms];  

    /* Open residual file and reset the CPU time. */
    first_step = 1;
    limiter_freezing_off = ON;
    if (Input_Parameters.i_ICs != IC_RESTART) {
      Input_Parameters.first_step = first_step;
    }

    if (CFFC_Primary_MPI_Processor()) {
      error_flag = Open_Progress_File(residual_file,
				      Input_Parameters.Output_File_Name,
				      number_of_time_steps,
				      Local_SolnBlk[0].residual_variable,
				      Local_SolnBlk[0].Number_of_Residual_Norms);
      //for unsteady plotting
      if( Input_Parameters.Time_Accurate_Plot_Frequency != 0){
	error_flag = Open_Time_Accurate_File(time_accurate_data_file,
					     Input_Parameters.Output_File_Name,
					     number_of_time_steps,
					     Local_SolnBlk[0].W[2][2]);
      }
      if (error_flag) {
	cout << "\n Flame2D ERROR: Unable to open residual file for Flame2D calculation.\n";
	cout.flush();
      } // endif
    } // endif
    
      // MPI barrier to ensure processor synchronization.
    CFFC_Barrier_MPI();
    CFFC_Broadcast_MPI(&error_flag, 1);
    if (error_flag) return (error_flag);
    
    processor_cpu_time.reset();
    
    /**************************************************************************
        Perform required number of iterations (time steps). 
    **************************************************************************/ 
    if ((!Input_Parameters.Time_Accurate && 
	 Input_Parameters.Maximum_Number_of_Time_Steps > 0) ||
	(Input_Parameters.Time_Accurate && Input_Parameters.Time_Max > Time)) {
     
      if (!batch_flag) { cout << "\n\n Beginning Explicit Flame2D computations on "
			      << Date_And_Time() << ".\n\n"; time(&start_explicit); /*start_explicit = clock();*/ }

      last_step = 0;
      int i=1;
      double number_explicit_steps = 0;

      while ((!Input_Parameters.Time_Accurate &&
	      number_explicit_steps < Input_Parameters.Maximum_Number_of_Time_Steps) ||
	     (Input_Parameters.Time_Accurate && Time < Input_Parameters.Time_Max)) {

	/***********************************************************************	
	   MORTON ORDERING of solution blocks during solution ever "n" steps
           ??? Should this be coupled with AMR Frequency ????
	************************************************************************/
	if (Input_Parameters.Morton && !first_step &&
	    number_of_time_steps%Input_Parameters.Morton_Reordering_Frequency == 0) {
	  if (!batch_flag) cout << "\n\n Applying Morton re-ordering algorithm to solution blocks at n = "
				<< number_of_time_steps << ".";
	  error_flag = Morton_ReOrdering_of_Solution_Blocks(QuadTree, 
							    List_of_Global_Solution_Blocks, 
							    List_of_Local_Solution_Blocks, 
							    Local_SolnBlk, 
							    Input_Parameters, 
							    number_of_time_steps, 
							    Time, 
							    processor_cpu_time); 
	  if (error_flag) {
	    cout <<"\n Flame2D ERROR: Morton re-ordering error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	    return (error_flag);
	  } 
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return (error_flag);
	  //Output space filling curve in Tecplot format
	  if (!batch_flag) cout << "\n Outputting space filling curve showing block loading for CPUs.";
	  Morton_SFC_Output_Tecplot(Local_SolnBlk, 
				    Input_Parameters, 
				    List_of_Local_Solution_Blocks);
	} 
       
	/***********************************************************************	
	   MESH REFINEMENT: Periodically refine the mesh (AMR). 
	************************************************************************/
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
			     ON, 
			     ON);
	    if (error_flag) {
	      cout << "\n Flame2D ERROR: Flame2D AMR error on processor "
		   << List_of_Local_Solution_Blocks.ThisCPU
		   << ".\n";
	      cout.flush();
	    } // endif
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
	    } // endif
	  } // endif
	} // endif

	  /********************** TIME STEPS **************************************
           Determine local and global time steps. 
	  *************************************************************************/

	  // Find global minimum time step for all processors.
	dTime = CFL(Local_SolnBlk, List_of_Local_Solution_Blocks, Input_Parameters);
	dTime = CFFC_Minimum_MPI(dTime);
	    
	if (Input_Parameters.Time_Accurate) {
	  if ((Input_Parameters.i_Time_Integration != 
	       TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) &&
	      (Time + Input_Parameters.CFL_Number*dTime > Input_Parameters.Time_Max)) {
	    dTime = (Input_Parameters.Time_Max-Time)/Input_Parameters.CFL_Number;
	    last_step = 1;
	       
	  } else if (Time + Input_Parameters.CFL_Number*dTime*
		     MultiStage_Optimally_Smoothing(Input_Parameters.N_Stage, 
						    Input_Parameters.N_Stage,
						    Input_Parameters.i_Limiter) > Input_Parameters.Time_Max) {
	      
	    dTime = (Input_Parameters.Time_Max-Time)/
	      (Input_Parameters.CFL_Number*
	       MultiStage_Optimally_Smoothing(Input_Parameters.N_Stage, 
					      Input_Parameters.N_Stage,
					      Input_Parameters.i_Limiter)); 
	    last_step = 1;
	  }
	} 	 
	  
	if (!Input_Parameters.Local_Time_Stepping) { 
	  // Set global time step.
	  if (Input_Parameters.Fixed_Time_Step) dTime = Input_Parameters.Time_Step;
	  Set_Global_TimeStep(Local_SolnBlk, List_of_Local_Solution_Blocks,dTime);
	} 

	/************************ NORMS *****************************************
           Determine the L1, L2, and max norms of the solution residual. 
	*************************************************************************/
	L1_Norm_Residual(Local_SolnBlk, List_of_Local_Solution_Blocks,residual_l1_norm);	  
	L2_Norm_Residual(Local_SolnBlk, List_of_Local_Solution_Blocks,residual_l2_norm);       	  
	Max_Norm_Residual(Local_SolnBlk, List_of_Local_Solution_Blocks,residual_max_norm);
	  
	for(int q=0; q < Local_SolnBlk[0].Number_of_Residual_Norms; q++){
	  residual_l1_norm[q] = CFFC_Summation_MPI(residual_l1_norm[q]); // L1 norm for all processors.
	  residual_l2_norm[q] =residual_l2_norm[q]*residual_l2_norm[q];
	  residual_l2_norm[q] = sqrt(CFFC_Summation_MPI(residual_l2_norm[q])); // L2 norm for all processors.
	  residual_max_norm[q] = CFFC_Maximum_MPI(residual_max_norm[q]); // Max norm for all processors.
	}
	  
	/* Update CPU time used for the calculation so far. */
	processor_cpu_time.update();
	// Total CPU time for all processors.
	total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput); 
	/************************ RESTART *****************************************
           Periodically save restart solution files. 
	***************************************************************************/ 
	    
	if (!first_step &&
	    number_of_time_steps%Input_Parameters.Restart_Solution_Save_Frequency == 0 ) {
	  if (!batch_flag){
	    cout << "\n\n  Saving Flame2D solution to restart data file(s) after"
		 << " n = " << number_of_time_steps << " steps (iterations).";
	  }
	  if (!batch_flag) cout << "\n  Writing Flame2D solution to restart data file(s). \n";
	  error_flag = Write_QuadTree(QuadTree,
				      Input_Parameters);
	  if (error_flag) {
	    cout << "\n Flame2D ERROR: Unable to open Flame2D quadtree data file "
		 << "on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  } 
	  error_flag = Write_Restart_Solution(Local_SolnBlk, 
					      List_of_Local_Solution_Blocks, 
					      Input_Parameters,
					      number_of_time_steps,
					      Time,
					      processor_cpu_time);
	  if (error_flag) {
	    cout << "\n Flame2D ERROR: Unable to open Flame2D restart output data file(s) "
		 << "on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  } 
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return (error_flag);
	  cout.flush();
	} 
	  
	/************************ PROGRESS *****************************************
           Output progress information for the calculation. 
	***************************************************************************/
	//screen 	    
	if (!batch_flag) Output_Progress_L2norm(number_of_time_steps,
						Time*THOUSAND, //time in ms
						total_cpu_time,
						residual_l2_norm[Local_SolnBlk[0].residual_variable-1],
						first_step,
						50);
	//residual to file
	if (CFFC_Primary_MPI_Processor() && !first_step) {
	  Output_Progress_to_File(residual_file,
				  number_of_time_steps,
				  Time*THOUSAND,
				  total_cpu_time,
				  residual_l1_norm,
				  residual_l2_norm,
				  residual_max_norm,
				  Local_SolnBlk[0].residual_variable,
				  Local_SolnBlk[0].Number_of_Residual_Norms);
	}
	//time vs mass frac (time accurate) for debugging
	//chemical solutions with no flow ie. just nonequilibrium chemistry
	if (Input_Parameters.Time_Accurate && Input_Parameters.Time_Accurate_Plot_Frequency != 0){
	  if ( (number_of_time_steps%Input_Parameters.Time_Accurate_Plot_Frequency) == 0 ){ 
	    Output_to_Time_Accurate_File(time_accurate_data_file,
					 Time,
					 Local_SolnBlk[0].W[2][2]);
	  }
	}

	/******************* CALCULATION CHECK ************************************
           Check to see if calculations are complete and if so jump of out of 
           this infinite loop.   

           Possibly replace while(1) so that these if's are the conditions 
           used or some sort of resonable facisimile.
	***************************************************************************/
 
	if (!Input_Parameters.Time_Accurate &&
	    number_explicit_steps >= 
	    Input_Parameters.Maximum_Number_of_Time_Steps) break;
	if (Input_Parameters.Time_Accurate && 
	    Time >= Input_Parameters.Time_Max) break;

	/******************* LIMITER FREEZE ***************************************	
	   Freeze limiters as necessary
	***************************************************************************/
	//         if (!first_step &&
	//             Input_Parameters.Freeze_Limiter &&
	//             limiter_freezing_off &&           
	//             residual_l2_norm[Local_SolnBlk[0].residual_variable-1] <= Input_Parameters.Freeze_Limiter_Residual_Level) {
	//   	   Freeze_Limiters(Local_SolnBlk, 
	// 			   List_of_Local_Solution_Blocks);
	// 	   limiter_freezing_off = ON;	
	//         } 

	/******************* BLOCK SOLUTION UPDATE ********************************
            Update solution for next time step using a multistage
            time stepping scheme. 
	***************************************************************************/
	for ( i_stage  = 1 ; i_stage <= Input_Parameters.N_Stage ; ++i_stage ) {
	  // 1. Exchange solution information between neighbouring blocks.
	  error_flag = Send_All_Messages(Local_SolnBlk, 
					 List_of_Local_Solution_Blocks,
					 NUM_VAR_FLAME2D, 
					 OFF);
	  if (error_flag) {
	    cout << "\n Flame2D ERROR: Flame2D message passing error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  } // endif
	    
	    // Reduce message passing error flag to other MPI processors.
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return (error_flag);
	    
	    
	  /************* BOUNDARY CONDITIONS *********************************/
	  // 2. Apply boundary conditions for stage.
	  BCs(Local_SolnBlk, List_of_Local_Solution_Blocks,Input_Parameters);
	    
	  /*************** UPDATE SOLUTION ************************************/
	  // 3. Determine solution residuals for stage.
	    
	  error_flag = dUdt_Multistage_Explicit(Local_SolnBlk,
						List_of_Global_Solution_Blocks,
						List_of_Local_Solution_Blocks,
						Input_Parameters,
						i_stage);
	  if (error_flag) {
	    cout << "\n Flame2D ERROR: Flame2D solution residual error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  } // endif
	    
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return (error_flag);
	    
	  // 4. Send boundary flux corrections at block interfaces with resolution changes.
	  error_flag = Send_Conservative_Flux_Corrections(Local_SolnBlk, 
							  List_of_Local_Solution_Blocks,
							  NUM_VAR_FLAME2D);
	  if (error_flag) {
	    cout << "\n Flame2D ERROR: Flame2D flux correction message passing error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  } // endif
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
	  } 
	    
	  // 7. Update solution for stage.
	  error_flag = Update_Solution_Multistage_Explicit(Local_SolnBlk, 
							   List_of_Local_Solution_Blocks,
							   Input_Parameters,
							   i_stage);
	    
	  if (error_flag) {
	    cout << "\n Flame2D ERROR: Flame2D solution update error on processor "
		 << List_of_Local_Solution_Blocks.ThisCPU
		 << ".\n";
	    cout.flush();
	  } 
	    
	  error_flag = CFFC_OR_MPI(error_flag);
	  if (error_flag) return (error_flag);
	    
	} /* endfor */
	
	  /******************* UPDATE TIMER & COUNTER *******************************
           Update time and time step counter. 
	  ***************************************************************************/
	if (first_step) {
	  first_step = 0;
	  Input_Parameters.first_step = 0;
	}

	number_of_time_steps++;
	number_explicit_steps++;

	// check for last step
	if (!Input_Parameters.Time_Accurate &&
	    number_explicit_steps == Input_Parameters.Maximum_Number_of_Time_Steps) {
	  last_step = 1;
	}
		
	if (Input_Parameters.i_Time_Integration != 
	    TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) {
	  Time = Time + Input_Parameters.CFL_Number*dTime;
 
	} else if (Input_Parameters.i_Time_Integration == 
		   TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) {
	  Time = Time + Input_Parameters.CFL_Number*dTime*
	    MultiStage_Optimally_Smoothing(Input_Parameters.N_Stage, 
					   Input_Parameters.N_Stage,
					   Input_Parameters.i_Limiter);
	} // endif
        
	i++; // is this necessary???
	
      } /* endwhile */
      
      if (!batch_flag) { cout << "\n\n Explicit Flame2D computations complete on " 
			      << Date_And_Time() << ".\n"; time(&end_explicit); /*end_explicit = clock();*/ }
          
    } // endif


      /************************************************************************************  
       Update ghostcell information and prescribe boundary conditions to ensure
       that the solution is consistent on each block. 
      *************************************************************************************/
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
    
    error_flag = Send_All_Messages(Local_SolnBlk, 
				   List_of_Local_Solution_Blocks,
				   NUM_VAR_FLAME2D,
				   OFF);
    
    if (error_flag) {
      cout << "\n Flame2D ERROR: Flame2D message passing error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } // endif
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    
    BCs(Local_SolnBlk,List_of_Local_Solution_Blocks, Input_Parameters);
    /* Close residual file. */
    
    if (CFFC_Primary_MPI_Processor()) {
      error_flag = Close_Progress_File(residual_file);
      error_flag = Close_Time_Accurate_File(time_accurate_data_file);
    } // endif
      
      //housecleaning
    delete[] residual_l1_norm;  delete[] residual_l2_norm;  delete[] residual_max_norm;

  }  //END EXPLICT NON-MULTIGRID

  /////////////////////////////////////////////////////////////////////////////
  /// APPLY Newton_Krylov_Schwarz
  /////////////////////////////////////////////////////////////////////////////

  if (Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0) {
    time_t start_NKS, end_NKS;

    if (CFFC_Primary_MPI_Processor()) {
      error_flag = Open_Progress_File(residual_file,
				      Input_Parameters.Output_File_Name,
				      number_of_time_steps,
				      Local_SolnBlk[0].residual_variable,
				      Local_SolnBlk[0].Number_of_Residual_Norms);
      if (error_flag) {
	cout << "\n Flame2D ERROR: Unable to open residual file for Chem2D calculation.\n";
	cout.flush();
      } /* endif */
    } /* endif */

    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
    CFFC_Broadcast_MPI(&error_flag, 1);
    if (error_flag) return (error_flag);

    //Turn Limiter Freezing OFF for startup
    Evaluate_Limiters(Local_SolnBlk, List_of_Local_Solution_Blocks);

    if (!batch_flag){ cout << "\n\n Beginning Flame2D NKS computations on " << Date_And_Time() << ".\n\n"; time(&start_NKS); }

    //Store Explicit times for output
    CPUTime Explicit_processor_cpu_time = processor_cpu_time;
    CPUTime Explicit_total_cpu_time =  total_cpu_time;
    
    //Perform NKS Iterations 
    error_flag = Newton_Krylov_Schwarz_Solver<Flame2D_pState,
      Flame2D_Quad_Block,
      Flame2D_Input_Parameters>(processor_cpu_time,
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
	cout << "\n Flame2D_NKS ERROR: Chem2D solution error on processor " 
	     << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
	cout.flush();
      } /* endif */
    } /* endif */

    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
    CFFC_Broadcast_MPI(&error_flag, 1);
    if (error_flag) return (error_flag);
    
    /***********************************************************************/
    if (!batch_flag) { cout << "\n\n Chem2D NKS computations complete on " << Date_And_Time() << ".\n"; time(&end_NKS); }

    if (!batch_flag) {
      cout<<"\n ----------------------------------------------------------------";
      cout<<"\n -------- Solution Computations Summary in minutes --------------";
      cout<<"\n ----------------------------------------------------------------";
      cout<<"\n Total Startup CPU Time\t\t = "<<Explicit_total_cpu_time.min();
      cout<<"\n Total NKS CPU Time \t\t = "<<total_cpu_time.min()-Explicit_total_cpu_time.min();
      cout<<"\n Total CPU Time \t\t = "<<total_cpu_time.min(); 
      cout<<"\n Total Startup Clock Time\t = "<<difftime(end_explicit,start_explicit)/60.0;
      cout<<"\n Total NKS Clock Time\t\t = "<<difftime(end_NKS,start_NKS)/60.0;
      cout<<"\n Total Clock Time\t\t = "<<difftime(end_NKS,start_explicit)/60.0;    //if no explicit start_eplicit not defined...
      cout<<"\n ----------------------------------------------------------------";
      cout<<"\n ----------------------------------------------------------------";
      cout<<"\n ----------------------------------------------------------------\n";
    } 
    //Also want to output total GMRES & NKS Iterations, and maybe max memory usage possibly??

    if (CFFC_Primary_MPI_Processor()) error_flag = Close_Progress_File(residual_file); 
     
  } 



  /***********************************************************************	
      RADIATION SOLVER solution of equation of radiation transfer to obtain
      the radiation heat flux term.
  ************************************************************************/
  if (Input_Parameters.Radiation == RADIATION_RTE &&
      Input_Parameters.Max_Number_Sequential_Solves > 0) {

#ifdef _FLAME2D_WITH_RTE

    if (!batch_flag) cout << "\n Solving radiation transfer equation...";

    // perform sequential solve
    RteSolver->SequentialSolve(Local_SolnBlk, Rte_PostProcess);
    
    if (!batch_flag) cout << "\n ...Done";

    // output some radiation solver stats
    if (!batch_flag) RteSolver->OutputSolverStats(cout);


    // increment sequential solve number
    number_sequential_solves++;

    //
    // Check for last step.
    // If it is not, continue an existing calculation.
    // This is not the nicest way of doing this... should be a for loop.
    //
    if (number_sequential_solves < Input_Parameters.Max_Number_Sequential_Solves) {

      // Continue existing calculation.     
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      goto continue_existing_calculation;
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      //
      // print out finalstats for sequential solver
    } else if (!batch_flag) {
      cout << endl << endl << string(70,'=');
      cout << "\n Sequential Solver Finished on " << Date_And_Time();
      cout << "\n Total Sequential Solves = " << number_sequential_solves;
      cout << "\n Total CPU Time          = " << total_cpu_time.min() << " min";
      cout << endl << string(70,'=') << endl;
    } // endif - continue

#endif // _FLAME2D_WITH_RTE

  } // endif - radiation


    /***************************************************************************
     ************************** POST PROCESSSING *******************************
     Solution calculations complete. Write 2D solution to output and restart files  
     as required, reset solution parameters, and run other cases as specified 
     by input parameters.        
     *************************************************************************** 
     ****************************************************************************/ 

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 postprocess_current_calculation: ;
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI(); 

  // To infinity and beyond....
  while (1) {
    
    if (CFFC_Primary_MPI_Processor()) {
      Get_Next_Input_Control_Parameter(Input_Parameters);
      command_flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
      line_number = Input_Parameters.Line_Number;
    } 

    //Debugging -- gives command directives from the end of the input file
    //cout<<" \n HERE "<<command_flag;
    //cout.flush();

    // MPI barrier to ensure processor synchronization.
    CFFC_Barrier_MPI();
    Broadcast_Input_Parameters(Input_Parameters);
    CFFC_Broadcast_MPI(&command_flag, 1);
    
    /************************************************************************
     **************** EXECUTE CODE ******************************************
     ************************************************************************/
    if (command_flag == EXECUTE_CODE) {
      // Deallocate memory for 2D Chem equation solution.
      if (!batch_flag) cout << "\n Deallocating Flame2D solution variables."; 
      Local_SolnBlk = Deallocate(Local_SolnBlk, 
				 Input_Parameters);
      //Deallocate_Message_Buffers(List_of_Local_Solution_Blocks); // Not necessary here!
      List_of_Local_Solution_Blocks.deallocate();
      List_of_Global_Solution_Blocks.deallocate();
      QuadTree.deallocate();
      MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
					    Input_Parameters.Number_of_Blocks_Idir, 
					    Input_Parameters.Number_of_Blocks_Jdir);
      Flame2D_Quad_Block::deallocate_static();

      // Output input parameters for new caluculation.
      if (!batch_flag)  {
	cout << "\n\n Starting a new calculation.";
	cout << Input_Parameters << "\n";
	cout.flush();
      } // endif
	// Execute new calculation.
      
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      goto execute_new_calculation;
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }
    /************************************************************************
     ******************** TERMINATE CODE ************************************ 
     ************************************************************************/
    else if (command_flag == TERMINATE_CODE) {
      // Deallocate memory for 2D Chem equation solution.
      if (!batch_flag) cout << "\n Deallocating Flame2D solution variables."; 
      Local_SolnBlk = Deallocate(Local_SolnBlk, 
				 Input_Parameters);

      //Deallocate_Message_Buffers(List_of_Local_Solution_Blocks); // Not necessary here!
      List_of_Local_Solution_Blocks.deallocate(); 
      List_of_Global_Solution_Blocks.deallocate();
      QuadTree.deallocate();
      MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
					    Input_Parameters.Number_of_Blocks_Idir, 
					    Input_Parameters.Number_of_Blocks_Jdir);
      Flame2D_Quad_Block::deallocate_static();

      // deallocate radiation solver
#ifdef _FLAME2D_WITH_RTE
      if (RteSolver != NULL) {
	delete RteSolver;
	RteSolver = NULL;
      } //endif
#endif // _FLAME2D_WITH_RTE

      // Close input data file.
      if (!batch_flag) cout << "\n\n Closing Flame2D input data file.";
      if (CFFC_Primary_MPI_Processor()) Close_Input_File(Input_Parameters);
      // Terminate calculation.
      return (0);
    }
    /************************************************************************
     ******************** CONTINUE CODE ie RESTART ***************************
     *************************************************************************/
    else if (command_flag == CONTINUE_CODE) {
      // Reset maximum time step counter.
      // Input_Parameters.Maximum_Number_of_Time_Steps += number_of_time_steps;

      // Output input parameters for continuing calculation.
      if (!batch_flag)  {
	cout << "\n\n Continuing existing calculation.";
	cout << Input_Parameters << "\n";
	cout.flush();
      } // endif
	// Continue existing calculation.     
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      goto continue_existing_calculation;
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }
    
    /************************************************************************
     ******************** REFINE GRID (AMR) *********************************
     *************************************************************************/
    else if (command_flag == REFINE_GRID_CODE) {
      // Refine mesh using block based adaptive mesh refinement algorithm.
      if (!batch_flag) cout << "\n\n Refining Grid.  Performing adaptive mesh refinement.";
      error_flag = AMR(Local_SolnBlk,
		       Input_Parameters,
		       QuadTree,
		       List_of_Global_Solution_Blocks,
		       List_of_Local_Solution_Blocks,
		       ON,
		       ON);
      if (error_flag) {
	cout << "\n Flame2D ERROR: Flame2D AMR error on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
      } // endif
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
      } // endif
      
    }
    /************************************************************************
     ********************* PERFORM MORTON ORDERING **************************
     *************************************************************************/
    else if (command_flag == MORTON_ORDERING_CODE) {
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
	cout <<"\n Flame2D ERROR: Morton re-ordering error on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
	return (error_flag);
      } // endif
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
      //Output space filling curve in Tecplot format
      if (!batch_flag) cout << "\n Outputting space filling curve showing block loading for CPUs.";
      Morton_SFC_Output_Tecplot(Local_SolnBlk, 
				Input_Parameters,
				List_of_Local_Solution_Blocks);

    }
    /************************************************************************
     ********************* WRITE OUTPUT @ NODES ******************************
     *************************************************************************/
    else if (command_flag == WRITE_OUTPUT_CODE) {
      // Output solution data.
      if (!batch_flag) cout << "\n Writing Flame2D solution to output data file(s).";
      
      error_flag = Output_Tecplot(Local_SolnBlk, 
				  List_of_Local_Solution_Blocks, 
				  Input_Parameters,
				  number_of_time_steps,
				  Time);
      if (error_flag) {
	cout << "\n Flame2D ERROR: Unable to open Flame2D output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
      } // endif
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
    }
 
    /*************************************************************************
     ********************** WRITE OUTPUT CELL-CENTERED ***********************
     *************************************************************************/
    
    else if (command_flag == WRITE_OUTPUT_CELLS_CODE) {
      // Output solution data.
      if (!batch_flag) cout << "\n Writing cell-centered Flame2D solution to output data file(s).";
      error_flag = Output_Cells_Tecplot(Local_SolnBlk, 
					List_of_Local_Solution_Blocks, 
					Input_Parameters,
					number_of_time_steps,
					Time);
      if (error_flag) {
	cout << "\n Flame2D ERROR: Unable to open Flame2D cell output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
      } // endif
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
      
    }
    /*************************************************************************
     ******************** WRITE RESTART FILE *********************************
     *************************************************************************/
    else if (command_flag == WRITE_RESTART_CODE) {
      // Write restart files.
      if (!batch_flag) cout << "\n Writing Flame2D solution to restart data file(s).";
      error_flag = Write_QuadTree(QuadTree,
				  Input_Parameters);
      if (error_flag) {
	cout << "\n Flame2D ERROR: Unable to open Flame2D quadtree data file "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
      } // endif
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
      
      error_flag = Write_Restart_Solution(Local_SolnBlk, 
					  List_of_Local_Solution_Blocks, 
					  Input_Parameters,
					  number_of_time_steps,
					  Time,
					  processor_cpu_time);
      if (error_flag) {
	cout << "\n Flame2D ERROR: Unable to open Flame2D restart output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
      } // endif
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
    }
  
    /*************************************************************************
     ******************** WRITE OUTPUT GRID **********************************
     *************************************************************************/
    else if (command_flag == WRITE_OUTPUT_GRID_CODE) {
      // Output multi-block solution-adaptive mesh data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing Flame2D multi-block mesh to grid data output file.";
	error_flag = Output_Tecplot(MeshBlk,
				    Input_Parameters);
	if (error_flag) {
	  cout << "\n Flame2D ERROR: Unable to open Flame2D mesh data output file.\n";
	  cout.flush();
	} // endif
      } // endif
      CFFC_Broadcast_MPI(&error_flag, 1);
      if (error_flag) return (error_flag);
    }

    /*************************************************************************
     **************** WRITE GRID DEFINITION **********************************
     *************************************************************************/
    else if (command_flag == WRITE_GRID_DEFINITION_CODE) {
      // Write multi-block solution-adaptive mesh definition files.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing Flame2D multi-block mesh to grid definition files.";
	error_flag = Write_Multi_Block_Grid_Definition(MeshBlk,
						       Input_Parameters);
	error_flag = Write_Multi_Block_Grid(MeshBlk,
					    Input_Parameters);
	if (error_flag) {
	  cout << "\n Flame2D ERROR: Unable to open Flame2D multi-block mesh definition files.\n";
	  cout.flush();
	} // endif
      } // endif
      CFFC_Broadcast_MPI(&error_flag, 1);
      if (error_flag) return (error_flag);
    }
    /*************************************************************************
     ******************** WRITE OUTPUT GRID NODES ****************************
     *************************************************************************/
    else if (command_flag == WRITE_OUTPUT_GRID_NODES_CODE) {
      // Output multi-block solution-adaptive mesh node data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing Flame2D multi-block mesh to node data output file.";
	error_flag = Output_Nodes_Tecplot(MeshBlk,
					  Input_Parameters);
	if (error_flag) {
	  cout << "\n Flame2D ERROR: Unable to open Flame2D mesh node data output file.\n";
	  cout.flush();
	} // endif
      } // endif
      CFFC_Broadcast_MPI(&error_flag, 1);
      if (error_flag) return (error_flag);
    }
    /*************************************************************************
     **************** WRITE OUTPUT GRID CELLS ********************************
     *************************************************************************/
    else if (command_flag == WRITE_OUTPUT_GRID_CELLS_CODE) {
      // Output multi-block solution-adaptive mesh cell data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing Flame2D multi-block mesh to cell data output file.";
	error_flag = Output_Cells_Tecplot(MeshBlk,
					  Input_Parameters);
	if (error_flag) {
	  cout << "\n Flame2D ERROR: Unable to open Flame2D mesh cell data output file.\n";
	  cout.flush();
	} // endif
      } // endif
      CFFC_Broadcast_MPI(&error_flag, 1);
      if (error_flag) return (error_flag);
    }

    /*************************************************************************
     **************** WRITE VICOUS CHANNEL ***********************************
     *************************************************************************/ 
    else if (command_flag == WRITE_OUTPUT_VISCOUS_CHANNEL_CODE) {
      if (!batch_flag) cout << endl << " Writing exact solution and error norms for the viscous channel flow.";
      error_flag = Output_Viscous_Channel(Local_SolnBlk,
					  List_of_Local_Solution_Blocks,
					  Input_Parameters);
      if (error_flag) {
	cout << endl << "\n Chem2D ERROR: Unable to open Chem2D viscous channel flow output file." << endl;
	cout.flush();
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
    }
    /*************************************************************************
     **************** WRITE FLAT PLATE ***************************************
     *************************************************************************/ 
    else if (command_flag == WRITE_OUTPUT_FLAT_PLATE_CODE) {
      if (!batch_flag) cout << endl << " Writing exact solution and error norms for the flat plate flow (Blasius solution).";
      error_flag = Output_Flat_Plate(Local_SolnBlk,
				     List_of_Local_Solution_Blocks,
				     Input_Parameters);
      if (error_flag) {
	cout << endl << "\n Chem2D ERROR: Unable to open Chem2D flat plate output file." << endl;
	cout.flush();
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
    } 
    /*************************************************************************
     **************** WRITE CAVITY DRIVEN FLOW *******************************
     *************************************************************************/ 
    else if (command_flag == WRITE_OUTPUT_DRIVEN_CAVITY_FLOW_CODE) {
      if (!batch_flag) cout << endl << " Writing the driven cavity flow output file.";
      error_flag = Output_Driven_Cavity_Flow(Local_SolnBlk,
					     List_of_Local_Solution_Blocks,
					     Input_Parameters);
      if (error_flag) {
	cout << endl << "\n Chem2D ERROR: Unable to open Chem2D driven cavity flow." << endl;
	cout.flush();
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
      
    } 

    /*************************************************************************
     ***************** POSTPROCESS 2D FLAME CODE *****************************
     *************************************************************************/ 
    else if (command_flag == POSTPROCESS_2DFLAME_CODE) {

      // Output solution data.
      if (!batch_flag) cout << "\n Postprocessing 2D flame.";
      

      error_flag = Output_2D_Flame(Local_SolnBlk,
				   List_of_Local_Solution_Blocks,
				   Input_Parameters);
      if (error_flag) {
	cout << endl << "\n Chem2D ERROR: Unable to process 2D flame." << endl;
	cout.flush();
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    }

    /*************************************************************************
     **************** POSTPROCESS RADIATION CODE *****************************
     *************************************************************************/ 
    else if (command_flag == POSTPROCESS_RADIATION_CODE) {

#ifdef _FLAME2D_WITH_RTE

      // Output solution data.
      if (!batch_flag) cout << "\n Postprocessing Rte2D solution.";
      
      error_flag = RteSolver->PostProcess(command_flag);
      if (error_flag) {
	cout << "\n Flame2D ERROR: Unable to postprocess Rte2D data."
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
      } // endif
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);

#endif //_FLAME2D_WITH_RTE

    }
 

    /*************************************************************************
     **************** NOT A VALID INPUT_CODE  ********************************
     *************************************************************************/  
    else if (command_flag == INVALID_INPUT_CODE ||
	     command_flag == INVALID_INPUT_VALUE) {
      line_number = -line_number;
      cout << "\n Flame2D ERROR: Error reading Flame2D data at line #"
	   << -line_number  << " of input data file.\n";
      cout.flush();
      return (line_number);
    } // endif
    
  } /* endwhile */

    /********************************************************  
     * End of all Flame2DSolver computations and I/O.       *
     ********************************************************/    
    //should use terminate code ie. should never get here?
  cout<<"\nEND OF FLAME2DSOLVER BUT DIDN'T USE TERMINATE!!\n"; cout.flush();
  
  return (0);
 
} //end Flame2DSolver
