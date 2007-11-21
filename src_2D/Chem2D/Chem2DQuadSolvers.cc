/************* Chem2DQuadSolvers.cc ***************************** 
   Chem2DQuadSolvers.cc:  2D Axisymmetric Solver for.... 
     
    NOTES:      

   - based on Euler2DQuadSolvers.cc
*****************************************************************/
#ifndef _CHEM2D_QUAD_INCLUDED
#include "Chem2DQuad.h"
#endif // _CHEM2D_QUAD_INCLUDED

/* Include CHEM2D Multigrid Specializations header file. */
#ifndef _CHEM_QUAD_MULTIGRID_INCLUDED
#include "Chem2DQuadMultigrid.h"
#endif //_CHEM_QUAD_MULTIGRID_INCLUDED

/* Include CHEM2D NKS Specializations header file. */
#ifndef _CHEM2D_NKS_INCLUDED 
#include "Chem2DQuadNKS.h"
#endif // _CHEM2D_NKS_INCLUDED 

/********************************************************
 * Routine: Chem2DQuadSolver                            *
 *                                                      *
 * Computes solutions to 2D .....  equations on 2D      *
 * quadrilateral multi-block solution-adaptive mesh.    *
 *                                                      *
 ********************************************************/
int Chem2DQuadSolver(char *Input_File_Name_ptr,  int batch_flag) {

  /********************************************************  
   * Local variable declarations.                         *
   ********************************************************/
  
  // Chem2D input variables and parameters:

  Chem2D_Input_Parameters Input_Parameters;
  
  /* Multi-block solution-adaptive quadrilateral mesh 
     solution variables. */
  Grid2D_Quad_Block          **MeshBlk;
  QuadTreeBlock_DataStructure  QuadTree;
  AdaptiveBlockResourceList    List_of_Global_Solution_Blocks;
  AdaptiveBlock2D_List         List_of_Local_Solution_Blocks;
  Chem2D_Quad_Block           *Local_SolnBlk;
  
  //Multigrid instantiation
  FAS_Multigrid2D_Solver<Chem2D_cState, 
                         Chem2D_Quad_Block, 
                         Chem2D_Input_Parameters> MGSolver;


  /* Define residual file and cpu time variables. */
  ofstream residual_file;
  ofstream time_accurate_data_file;
  CPUTime processor_cpu_time, total_cpu_time;
  time_t start_explicit, end_explicit;

  /* Other local solution variables. */
  int number_of_time_steps, first_step, last_step,
    command_flag, error_flag, line_number, 
    i_stage, limiter_freezing_off;
  
  double Time, dTime, initial_residual_l2_norm;

  /* Variables used for dual time stepping. */
  int n_inner = 0, n_inner_temp = 0;
  const double dual_eps = 1.0E-3;  // 5.0E-4;

  /*************************************************************************
   ******************** INPUT PARAMETERS  **********************************
     Set default values for the input solution parameters and then read user 
     specified input values from the specified input parameter file.            
   *************************************************************************
   *************************************************************************/  
  //The primary MPI processor processes the input parameter file.
  if (CFFC_Primary_MPI_Processor()) {
    if (!batch_flag) {
        cout << "\n Reading Chem2D input data file `"
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
 
  /* Common Use Variables */
  const Chem2D_pState Chem2D_W_STDATM;
  const Chem2D_pState Chem2D_W_VACUUM(ZERO, Vector2D_ZERO, ZERO,ZERO,ZERO);
  const Chem2D_cState Chem2D_U_STDATM; //default
  const Chem2D_cState Chem2D_U_VACUUM(ZERO, Vector2D_ZERO, ZERO,ZERO, ZERO);
  const int NUM_VAR_CHEM2D = Input_Parameters.Wo.NUM_VAR_CHEM2D;
  
  /*************************************************************************
   ******************** INITIAL GRID & SOLUTION SPACE **********************
   Create initial mesh and allocate Chem2D solution variables for 
   specified IBVP/BVP problem. 
   *************************************************************************
   *************************************************************************/

  execute_new_calculation: ;
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
  
  /* Create initial mesh.  Read mesh from grid definition or data files 
     when specified by input parameters. */
  
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
    } /* endif */
    
     if (error_flag) {
       cout << "\n Chem2D ERROR: Unable to create valid Chem2D multi-block mesh.\n";
       cout.flush();
     } /* endif */
  } else {
    MeshBlk = NULL;
  } /* endif */

  // Broadcast the mesh to other MPI processors.
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
  CFFC_Broadcast_MPI(&error_flag, 1); // Broadcast mesh error flag.

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
	 << "\n Chem2D solution blocks corresponding to initial mesh.";
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
   * Initialize Chem2D solution variables, Initial conditions             *  
   ************************************************************************/

  /* Set the initial time level. */
  Time = ZERO;
  number_of_time_steps = 0;

  /* Set the CPU time to zero. */
  processor_cpu_time.zero();
  total_cpu_time.zero(); 

  /* Initialize the conserved and primitive state solution variables. */
  if (!batch_flag){ 
    cout << "\n Prescribing Chem2D initial data.";
  }
 
  /*************************************************************************
   **************** RESTART SOLUTION or INITIAL CONDITIONS *****************
   *************************************************************************/
  if (Input_Parameters.i_ICs == IC_RESTART) {
    if (!batch_flag) cout << "\n Reading Chem2D solution from restart data files.";
    error_flag = Read_QuadTree(QuadTree,
			       List_of_Global_Solution_Blocks,
			       List_of_Local_Solution_Blocks, 
			       Input_Parameters);
    if (error_flag) {
      cout << "\n Chem2D ERROR: Unable to open Chem2D quadtree data file "
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
         cout << "\n Chem2D ERROR: Unable to open Chem2D restart input data file(s) "
              << "on processor "
              << List_of_Local_Solution_Blocks.ThisCPU
              << ".\n";
         cout.flush();
    } /* endif */
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

    /****** Else apply initial conditions from input parameters *******/
  } else {   

    // Send grid information between neighbouring blocks BEFORE applying ICs
    CFFC_Barrier_MPI(); 
    error_flag = Send_All_Messages(Local_SolnBlk, 
				   List_of_Local_Solution_Blocks,
				   NUM_COMP_VECTOR2D,
				   ON);     

    ICs(Local_SolnBlk, List_of_Local_Solution_Blocks, Input_Parameters);
  } /* endif */
  
  
  /******************************************************************************  
    Send solution information between neighbouring blocks to complete
    prescription of initial data. 
  *******************************************************************************/
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
 
  // OLD grid send_all_messages, moved before ICs to help with values used in BCs 
  // set from ICs, ie ghost cell node locations, changed before ICs 

  if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk, 
                                                  List_of_Local_Solution_Blocks,
                                                  NUM_VAR_CHEM2D,
                                                  OFF);
  if (error_flag) {
    cout << "\n Chem2D ERROR: Message passing error during Chem2D solution intialization "
	 << "on processor "
	 << List_of_Local_Solution_Blocks.ThisCPU
	 << ".\n";
    cout.flush();
    exit(1);
  } /* endif */ 

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
  /* Perform uniform, boundary, and, initial mesh refinement. */
  if (Input_Parameters.i_ICs != IC_RESTART) {
     if (!batch_flag) cout << "\n Performing Chem2D uniform mesh refinement.";
     error_flag = Uniform_AMR(Local_SolnBlk,
                              Input_Parameters,
                              QuadTree,
                              List_of_Global_Solution_Blocks,
                              List_of_Local_Solution_Blocks);
     if (error_flag) {
       cout << "\n Chem2D ERROR: Uniform AMR error on processor "
	    << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
       cout.flush();
     } /* endif */
     error_flag = CFFC_OR_MPI(error_flag);
     if (error_flag) return error_flag;

     ///////////////////////////////////////////////////////////////////////////
     if (!batch_flag) cout << "\n Performing Chem2D boundary mesh refinement.";
     error_flag = Boundary_AMR(Local_SolnBlk,
                               Input_Parameters,
                               QuadTree,
                               List_of_Global_Solution_Blocks,
                               List_of_Local_Solution_Blocks);
     if (error_flag) {
       cout << "\n Chem2D ERROR: Boundary AMR error on processor "
	    << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
       cout.flush();
     } /* endif */
     error_flag = CFFC_OR_MPI(error_flag);
     if (error_flag) return error_flag;

     ///////////////////////////////////////////////////////////////////////////
     if (!batch_flag) cout << "\n Performing Chem2D initial mesh refinement.";
     error_flag = Initial_AMR(Local_SolnBlk,
                              Input_Parameters,
	   	              QuadTree,
		              List_of_Global_Solution_Blocks,
		              List_of_Local_Solution_Blocks);
     if (error_flag) {
        cout << "\n Chem2D ERROR: Initial AMR error on processor "
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
      cout <<"\n Chem2D ERROR: Morton re-ordering error on processor "
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
  
  /****************************************************************************
   *********************** MAIN SOLVER ****************************************
   Solve IBVP or BVP for conservation form of 2D Axisymmetric multispecies 
   chemically reacting thermally perfect equations on multi-block 
   solution-adaptive quadrilateral mesh.                                  
   ****************************************************************************
   ****************************************************************************/  
  continue_existing_calculation: ;
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
 
  /******************* MULTIGRID SETUP ****************************************/
  if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {
 
    MGSolver.allocate(Local_SolnBlk,
		      &QuadTree,
		      &List_of_Global_Solution_Blocks,
		      &List_of_Local_Solution_Blocks,
		      &Input_Parameters); 

    if(!batch_flag) time(&start_explicit);

    error_flag = MGSolver.Execute(batch_flag,
				  number_of_time_steps,
				  Time,
				  processor_cpu_time,
				  total_cpu_time,
				  residual_file);

    if(!batch_flag) time(&end_explicit);
    
  /*********************** NON-MULTIGRID EXPLICT ***********************************/
  } else if(Input_Parameters.Maximum_Number_of_Time_Steps > 0){ 
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
        cout << "\n Chem2D ERROR: Unable to open residual file for Chem2D calculation.\n";
        cout.flush();
      } /* endif */
    } /* endif */
    
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
    CFFC_Broadcast_MPI(&error_flag, 1);
    if (error_flag) return (error_flag);
    
    processor_cpu_time.reset();
    
    /**************************************************************************
     Perform required number of iterations (time steps). 
    **************************************************************************/ 
    if ((!Input_Parameters.Time_Accurate && 
	 Input_Parameters.Maximum_Number_of_Time_Steps > 0 &&
	 number_of_time_steps < Input_Parameters.Maximum_Number_of_Time_Steps) ||
	(Input_Parameters.Time_Accurate && Input_Parameters.Time_Max > Time)) {
     
      if (!batch_flag) { cout << "\n\n Beginning Explicit Chem2D computations on "
			      << Date_And_Time() << ".\n\n"; time(&start_explicit); /*start_explicit = clock();*/ }

     last_step = 0;
     int i=1; //markthis

     while ((!Input_Parameters.Time_Accurate &&
	     number_of_time_steps < Input_Parameters.Maximum_Number_of_Time_Steps) ||
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
	   cout <<"\n Chem2D ERROR: Morton re-ordering error on processor "
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
                 cout << "\n Chem2D ERROR: Chem2D AMR error on processor "
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

// 	n_inner_temp = n_inner;
//         n_inner = 0;

// 	do{
	  /********************** TIME STEPS **************************************
           Determine local and global time steps. 
	  *************************************************************************/

	  dTime = CFL(Local_SolnBlk, List_of_Local_Solution_Blocks, Input_Parameters);
	  if (!Input_Parameters.Dual_Time_Stepping) {
	    // Find global minimum time step for all processors.
	    dTime = CFFC_Minimum_MPI(dTime);
	  } else {
            // Assign physical time step for dual time stepping.
            if (n_inner == 0) { 
	      dTime = Input_Parameters.Physical_CFL_Number*CFFC_Minimum_MPI(dTime);
              Input_Parameters.dTime = dTime;              
	    }
            dTime = Input_Parameters.dTime;
	  }

	  if (Input_Parameters.Time_Accurate) {
	    if ((Input_Parameters.i_Time_Integration != 
		 TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) &&
		(Time + Input_Parameters.CFL_Number*dTime > Input_Parameters.Time_Max)) {
	      dTime = (Input_Parameters.Time_Max-Time)/Input_Parameters.CFL_Number;
	      last_step = 1;
	       
	    } else if (Time + Input_Parameters.CFL_Number*dTime*
		       MultiStage_Optimally_Smoothing(Input_Parameters.N_Stage, 
						      Input_Parameters.N_Stage,
						      Input_Parameters.i_Limiter) > 
		       Input_Parameters.Time_Max && (!Input_Parameters.Dual_Time_Stepping)) {
	      
	      dTime = (Input_Parameters.Time_Max-Time)/
		(Input_Parameters.CFL_Number*
		 MultiStage_Optimally_Smoothing(Input_Parameters.N_Stage, 
						Input_Parameters.N_Stage,
						Input_Parameters.i_Limiter)); 
	      last_step = 1;
	   
	    } else if (Time + dTime > Input_Parameters.Time_Max && Input_Parameters.Dual_Time_Stepping) {
	      dTime = Input_Parameters.Time_Max - Time;
              Input_Parameters.dTime = dTime;
              last_step = 1;
            } 
	  } 	 
	  
	  if (!Input_Parameters.Local_Time_Stepping) { 
	    // Set global time step.	
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
	  if (n_inner == 1) initial_residual_l2_norm = residual_l2_norm[Local_SolnBlk[0].residual_variable-1];
	  
	  /* Update CPU time used for the calculation so far. */
	  processor_cpu_time.update();
	  // Total CPU time for all processors.
	  total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput); 
	  /************************ RESTART *****************************************
          Periodically save restart solution files. 
	  ***************************************************************************/ 
	  if (!Input_Parameters.Dual_Time_Stepping  ||  
	      (Input_Parameters.Dual_Time_Stepping && n_inner == 0)) {
	    
	    if (!first_step &&
		number_of_time_steps%Input_Parameters.Restart_Solution_Save_Frequency == 0 ) {
	      if (!batch_flag){
		cout << "\n\n  Saving Chem2D solution to restart data file(s) after"
		     << " n = " << number_of_time_steps << " steps (iterations).";
	      }
	      if (!batch_flag) cout << "\n  Writing Chem2D solution to restart data file(s). \n";
	      error_flag = Write_QuadTree(QuadTree,
					  Input_Parameters);
	      if (error_flag) {
		cout << "\n Chem2D ERROR: Unable to open Chem2D quadtree data file "
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
		cout << "\n Chem2D ERROR: Unable to open Chem2D restart output data file(s) "
		     << "on processor "
		     << List_of_Local_Solution_Blocks.ThisCPU
		     << ".\n";
		cout.flush();
	      } 
	      error_flag = CFFC_OR_MPI(error_flag);
	      if (error_flag) return (error_flag);
	      cout.flush();
	    } 
	  }
	  
	  /************************ PROGRESS *****************************************
           Output progress information for the calculation. 
	  ***************************************************************************/
	  //screen 
	  if (!Input_Parameters.Dual_Time_Stepping || 
	      (Input_Parameters.Dual_Time_Stepping && n_inner == 0)) {
	    
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
	  }
	  /******************* CALCULATION CHECK ************************************
           Check to see if calculations are complete and if so jump of out of 
           this infinite loop.   

           Possibly replace while(1) so that these if's are the conditions 
           used or some sort of resonable facisimile.
	  ***************************************************************************/
 
	if (!Input_Parameters.Time_Accurate &&
 	    number_of_time_steps >= 
 	    Input_Parameters.Maximum_Number_of_Time_Steps) break;
 	if (Input_Parameters.Time_Accurate && 
 	    Time >= Input_Parameters.Time_Max) break;
	if (Input_Parameters.Dual_Time_Stepping && 
	    (n_inner == Input_Parameters.Max_Inner_Steps ||
	     (n_inner > 0  && (residual_l2_norm[Local_SolnBlk[0].residual_variable-1] < dual_eps  || 
			       residual_l2_norm[Local_SolnBlk[0].residual_variable-1]/
			       initial_residual_l2_norm < dual_eps)))) break;

// 	/******************* LIMITER FREEZE ***************************************	
// 	 Freeze limiters as necessary
// 	***************************************************************************/
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
					   NUM_VAR_CHEM2D, 
					   OFF);
	    if (error_flag) {
	      cout << "\n Chem2D ERROR: Chem2D message passing error on processor "
		   << List_of_Local_Solution_Blocks.ThisCPU
		   << ".\n";
	      cout.flush();
	    } /* endif */
	    
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
	      cout << "\n Chem2D ERROR: Chem2D solution residual error on processor "
		   << List_of_Local_Solution_Blocks.ThisCPU
		   << ".\n";
	      cout.flush();
	    } /* endif */
	    
	    error_flag = CFFC_OR_MPI(error_flag);
	    if (error_flag) return (error_flag);
	    
	    // 4. Send boundary flux corrections at block interfaces with resolution changes.
	    error_flag = Send_Conservative_Flux_Corrections(Local_SolnBlk, 
							    List_of_Local_Solution_Blocks,
							    NUM_VAR_CHEM2D);
	    if (error_flag) {
	      cout << "\n Chem2D ERROR: Chem2D flux correction message passing error on processor "
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
	    } 
	    
	    // 7. Update solution for stage.
	    error_flag = Update_Solution_Multistage_Explicit(Local_SolnBlk, 
							     List_of_Local_Solution_Blocks,
							     Input_Parameters,
							     i_stage);
	    
	    if (error_flag) {
	      cout << "\n Chem2D ERROR: Chem2D solution update error on processor "
		   << List_of_Local_Solution_Blocks.ThisCPU
		   << ".\n";
	      cout.flush();
	    } 
	    
	    error_flag = CFFC_OR_MPI(error_flag);
	    if (error_flag) return (error_flag);
	    
	  } /* endfor */

	  n_inner++;
// 	} while (Input_Parameters.Dual_Time_Stepping && 
//                  (n_inner < Input_Parameters.Max_Inner_Steps+1  || first_step)); 

// 	/********** UPDATE STATES FOR DUAL TIME STEPPING ***************************  
//          Update solution states, at different times, required by dual time stepping.
//         ****************************************************************************/
//         if (Input_Parameters.Dual_Time_Stepping) {
// 	  error_flag = Update_Dual_Solution_States(Local_SolnBlk, 
// 						 List_of_Local_Solution_Blocks);
// 	  if (error_flag) {
// 	    cout << "\n Chem2D ERROR: Chem2D solution states update error on processor "
// 		 << List_of_Local_Solution_Blocks.ThisCPU
// 		 << ".\n";
// 	    cout.flush();
// 	  }
// 	}
	
	/******************* UPDATE TIMER & COUNTER *******************************
          Update time and time step counter. 
	***************************************************************************/
	if (first_step) {
	  first_step = 0;
	  Input_Parameters.first_step = 0;
	}

	number_of_time_steps = number_of_time_steps + 1;

	// check for last step
        if (!Input_Parameters.Time_Accurate &&
	    number_of_time_steps == Input_Parameters.Maximum_Number_of_Time_Steps) {
          last_step = 1;
	}
		
	if (Input_Parameters.i_Time_Integration != 
	  TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING  &&
          !Input_Parameters.Dual_Time_Stepping) {
	  Time = Time + Input_Parameters.CFL_Number*dTime;
 
	} else if (Input_Parameters.i_Time_Integration == 
	  TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING  &&
	  !Input_Parameters.Dual_Time_Stepping) {
	  Time = Time + Input_Parameters.CFL_Number*dTime*
	    MultiStage_Optimally_Smoothing(Input_Parameters.N_Stage, 
					   Input_Parameters.N_Stage,
					   Input_Parameters.i_Limiter);
	} else if (Input_Parameters.Dual_Time_Stepping) {
          Time = Time + dTime;
	} /* endif */
        
	i++; // is this necessary???
	
      } /* endwhile */
      
      if (!batch_flag) { cout << "\n\n Explicit Chem2D computations complete on " 
			      << Date_And_Time() << ".\n"; time(&end_explicit); /*end_explicit = clock();*/ }
          
    } /* endif */


    /************************************************************************************  
    Update ghostcell information and prescribe boundary conditions to ensure
    that the solution is consistent on each block. 
    *************************************************************************************/
    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
    
    error_flag = Send_All_Messages(Local_SolnBlk, 
				   List_of_Local_Solution_Blocks,
				   NUM_VAR_CHEM2D,
				   OFF);
    
    if (error_flag) {
      cout << "\n Chem2D ERROR: Chem2D message passing error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    
    BCs(Local_SolnBlk,List_of_Local_Solution_Blocks, Input_Parameters);
    /* Close residual file. */
    
    if (CFFC_Primary_MPI_Processor()) {
      error_flag = Close_Progress_File(residual_file);
      error_flag = Close_Time_Accurate_File(time_accurate_data_file);
    } /* endif */
      
    //housecleaning
    delete[] residual_l1_norm;  delete[] residual_l2_norm;  delete[] residual_max_norm;

  }  //END EXPLICT NON-MULTIGRID

  /*************************************************************************************************************************/
  /************************ APPLY Newton_Krylov_Schwarz ********************************************************************/
  /*************************************************************************************************************************/
  if (Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0) {
     time_t start_NKS, end_NKS;

     if (CFFC_Primary_MPI_Processor()) {
        error_flag = Open_Progress_File(residual_file,
	  			        Input_Parameters.Output_File_Name,
				        number_of_time_steps,
					Local_SolnBlk[0].residual_variable,
					Local_SolnBlk[0].Number_of_Residual_Norms);
        if (error_flag) {
           cout << "\n Chem2D ERROR: Unable to open residual file for Chem2D calculation.\n";
           cout.flush();
        } /* endif */
     } /* endif */

     CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
     CFFC_Broadcast_MPI(&error_flag, 1);
     if (error_flag) return (error_flag);

     //Turn Limiter Freezing OFF for startup
     Evaluate_Limiters(Local_SolnBlk, List_of_Local_Solution_Blocks);

     if (!batch_flag){ cout << "\n\n Beginning Chem2D NKS computations on " << Date_And_Time() << ".\n\n"; time(&start_NKS); }

     //Store Explicit times for output
     CPUTime Explicit_processor_cpu_time = processor_cpu_time;
     CPUTime Explicit_total_cpu_time =  total_cpu_time;
    
     //Perform NKS Iterations 
     error_flag = Newton_Krylov_Schwarz_Solver<Chem2D_pState,
                                               Chem2D_Quad_Block,                                               
                                               Chem2D_Input_Parameters>(processor_cpu_time,
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
   	   cout << "\n Chem2D_NKS ERROR: Chem2D solution error on processor " 
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
  /*************************************************************************************************************************/
  /*************************************************************************************************************************/
  /*************************************************************************************************************************/

  /***************************************************************************
   ************************** POST PROCESSSING *******************************
    Solution calculations complete. Write 2D solution to output and restart files  
    as required, reset solution parameters, and run other cases as specified 
    by input parameters.        
   *************************************************************************** 
   ****************************************************************************/ 
  
  postprocess_current_calculation: ;
  CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

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

    CFFC_Barrier_MPI(); // MPI barrier to ensure processor synchronization.
    Broadcast_Input_Parameters(Input_Parameters);
    CFFC_Broadcast_MPI(&command_flag, 1);
    
    /************************************************************************
     **************** EXECUTE CODE ******************************************
     ************************************************************************/
    if (command_flag == EXECUTE_CODE) {
      // Deallocate memory for 2D Chem equation solution.
      if (!batch_flag) cout << "\n Deallocating Chem2D solution variables."; 
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
      
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      goto execute_new_calculation;
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }
    /************************************************************************
     ******************** TERMINATE CODE ************************************ 
     ************************************************************************/
    else if (command_flag == TERMINATE_CODE) {
      // Deallocate memory for 2D Chem equation solution.
      if (!batch_flag) cout << "\n Deallocating Chem2D solution variables."; 
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
      if (!batch_flag) cout << "\n\n Closing Chem2D input data file.";
      if (CFFC_Primary_MPI_Processor()) Close_Input_File(Input_Parameters);
      // Terminate calculation.
      return (0);
    }
    /************************************************************************
    ******************** CONTINUE CODE ie RESTART ***************************
    *************************************************************************/
    else if (command_flag == CONTINUE_CODE) {
      // Reset maximum time step counter.
      Input_Parameters.Maximum_Number_of_Time_Steps += number_of_time_steps;
      // Output input parameters for continuing calculation.
      if (!batch_flag)  {
	cout << "\n\n Continuing existing calculation.";
	cout << Input_Parameters << "\n";
	cout.flush();
      } /* endif */
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
	cout << "\n Chem2D ERROR: Chem2D AMR error on processor "
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
           cout <<"\n Chem2D ERROR: Morton re-ordering error on processor "
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

    }
    /************************************************************************
     ********************* WRITE OUTPUT @ NODES ******************************
     *************************************************************************/
    else if (command_flag == WRITE_OUTPUT_CODE) {
      // Output solution data.
      if (!batch_flag) cout << "\n Writing Chem2D solution to output data file(s).";
      
      error_flag = Output_Tecplot(Local_SolnBlk, 
				  List_of_Local_Solution_Blocks, 
				  Input_Parameters,
				  number_of_time_steps,
				  Time);
      if (error_flag) {
	cout << "\n Chem2D ERROR: Unable to open Chem2D output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
      } /* endif */
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
    }
 
    /*************************************************************************
     ********************** WRITE OUTPUT CELL-CENTERED ***********************
     *************************************************************************/
    
    else if (command_flag == WRITE_OUTPUT_CELLS_CODE) {
      // Output solution data.
      if (!batch_flag) cout << "\n Writing cell-centered Chem2D solution to output data file(s).";
      error_flag = Output_Cells_Tecplot(Local_SolnBlk, 
					List_of_Local_Solution_Blocks, 
					Input_Parameters,
					number_of_time_steps,
					Time);
      if (error_flag) {
	cout << "\n Chem2D ERROR: Unable to open Chem2D cell output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
      } /* endif */
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
      
    }
    /************************************************************************
     ********************* WRITE OUTPUT RightHandSide ***********************
     *************************************************************************/
	 else if (command_flag == WRITE_OUTPUT_RHS_CODE) {
      // Output the RHS
      if (!batch_flag) cout << "\n Writing right hand side to output data file(s).";
      
      error_flag = Output_RHS(Local_SolnBlk, 
				  List_of_Local_Solution_Blocks, 
				  Input_Parameters,
				  number_of_time_steps,
				  Time);
      if (error_flag) {
	cout << "\n Chem2D ERROR: Unable to open Chem2D output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
      } /* endif */
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
    }
    /************************************************************************
     ********************* WRITE OUTPUT RightHandSide after Perturbation******
     *************************************************************************/
    else if (command_flag == WRITE_OUTPUT_PERTURB_CODE) {
      // Output the RHS
      if (!batch_flag) cout << "\n Writing right hand side to output data file(s).";
      
      error_flag = Output_PERTURB(Local_SolnBlk, 
					  List_of_Local_Solution_Blocks, 
					  Input_Parameters,
					  number_of_time_steps,
					  Time,
					  processor_cpu_time);
 
      if (error_flag) {
	cout << "\n Chem2D ERROR: Unable to open Chem2D output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
      } /* endif */
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
    }
    /*************************************************************************
     ******************** WRITE RESTART FILE *********************************
     *************************************************************************/
    else if (command_flag == WRITE_RESTART_CODE) {
      // Write restart files.
      if (!batch_flag) cout << "\n Writing Chem2D solution to restart data file(s).";
      error_flag = Write_QuadTree(QuadTree,
				  Input_Parameters);
      if (error_flag) {
	cout << "\n Chem2D ERROR: Unable to open Chem2D quadtree data file "
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
	cout << "\n Chem2D ERROR: Unable to open Chem2D restart output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
      } /* endif */
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
    }
  
    /*************************************************************************
     ******************** WRITE OUTPUT GRID **********************************
     *************************************************************************/
    else if (command_flag == WRITE_OUTPUT_GRID_CODE) {
      // Output multi-block solution-adaptive mesh data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing Chem2D multi-block mesh to grid data output file.";
	error_flag = Output_Tecplot(MeshBlk,
				    Input_Parameters);
	if (error_flag) {
	  cout << "\n Chem2D ERROR: Unable to open Chem2D mesh data output file.\n";
	  cout.flush();
	} /* endif */
      } /* endif */
      CFFC_Broadcast_MPI(&error_flag, 1);
      if (error_flag) return (error_flag);
    }

    /*************************************************************************
     **************** WRITE GRID DEFINITION **********************************
     *************************************************************************/
    else if (command_flag == WRITE_GRID_DEFINITION_CODE) {
      // Write multi-block solution-adaptive mesh definition files.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing Chem2D multi-block mesh to grid definition files.";
	error_flag = Write_Multi_Block_Grid_Definition(MeshBlk,
						       Input_Parameters);
	error_flag = Write_Multi_Block_Grid(MeshBlk,
					    Input_Parameters);
	if (error_flag) {
	  cout << "\n Chem2D ERROR: Unable to open Chem2D multi-block mesh definition files.\n";
	  cout.flush();
	} /* endif */
      } /* endif */
      CFFC_Broadcast_MPI(&error_flag, 1);
      if (error_flag) return (error_flag);
    }
    /*************************************************************************
     ******************** WRITE OUTPUT GRID NODES ****************************
     *************************************************************************/
    else if (command_flag == WRITE_OUTPUT_GRID_NODES_CODE) {
         // Output multi-block solution-adaptive mesh node data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing Chem2D multi-block mesh to node data output file.";
	error_flag = Output_Nodes_Tecplot(MeshBlk,
					  Input_Parameters);
	if (error_flag) {
	  cout << "\n Chem2D ERROR: Unable to open Chem2D mesh node data output file.\n";
	  cout.flush();
	} /* endif */
      } /* endif */
      CFFC_Broadcast_MPI(&error_flag, 1);
      if (error_flag) return (error_flag);
    }
    /*************************************************************************
     **************** WRITE OUTPUT GRID CELLS ********************************
     *************************************************************************/
    else if (command_flag == WRITE_OUTPUT_GRID_CELLS_CODE) {
      // Output multi-block solution-adaptive mesh cell data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing Chem2D multi-block mesh to cell data output file.";
	error_flag = Output_Cells_Tecplot(MeshBlk,
					  Input_Parameters);
	if (error_flag) {
	  cout << "\n Chem2D ERROR: Unable to open Chem2D mesh cell data output file.\n";
	  cout.flush();
	} /* endif */
      } /* endif */
      CFFC_Broadcast_MPI(&error_flag, 1);
      if (error_flag) return (error_flag);
    }
    /*************************************************************************
     **************** WRITE RINGLEB ******************************************
     *************************************************************************/
    else if (command_flag == WRITE_OUTPUT_RINGLEB_CODE) {
      if (!batch_flag) cout << endl << " Writing exact solution and error norms for Ringleb's flow.";
      error_flag = Output_Ringleb(Local_SolnBlk,
				  List_of_Local_Solution_Blocks,
				  Input_Parameters);
      if (error_flag) {
	cout << endl << "\n Chem2D ERROR: Unable to open Chem2D Ringleb's flow output file." << endl;
	cout.flush();
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
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
     **************** WRITE QUASI 3D for AXISYMMETRIC ************************
     *************************************************************************/ 
    else if (command_flag == WRITE_OUTPUT_QUASI3D_CODE) {
      // Output solution data.
      if (!batch_flag) cout << "\n Writing Chem2D quasi3D solution to output data file(s).";
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
 
    }
    /*************************************************************************
     **************** SWITCH BCs TO FIXED  ***********************************
     *************************************************************************/ 
    else if (command_flag == SWITCH_BCS_TO_FIXED) {      
      if (!batch_flag) cout << endl << " Setting EAST & WEST BC's Fixed. -> WARNING: ONLY VALID FOR 1D FLAME CASE!!!! \n";

      //Should be calling a quadmultiblock function that would do this per block,
      //but since the 1D FLAME case only has one block this should work, even though it is a hack.
      
      Local_SolnBlk[0].Grid.set_BCs(WEST,BC_CONSTANT_EXTRAPOLATION);      
      Local_SolnBlk[0].Grid.set_BCs(EAST,BC_CONSTANT_EXTRAPOLATION); 
      
      //SET all V-velocity to ZERO
      Local_SolnBlk[0].set_v_zero();  //????

      //Switch to const_extrap then run BCS to set ghost cell W values to internal values
      BCs(Local_SolnBlk, List_of_Local_Solution_Blocks,Input_Parameters);

      //Copy adjusted ghost cell values to overwrite previous WoS, WoN, etc..
      Reset_Wo(Local_SolnBlk[0]);

      //Swich to Fixed and now BCs will use updated WoS, WoN, etc...
      Local_SolnBlk[0].Grid.set_BCs(WEST,BC_FIXED);      
      Local_SolnBlk[0].Grid.set_BCs(EAST,BC_FIXED); 

      if (error_flag) {
	cout << endl << "\n Chem2D ERROR: Problem in Chem2D SWITCH_BCS_TO_FIXED. " << endl;
	cout.flush();
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;      
    }
    /*************************************************************************
     **************** NOT A VALID INPUT_CODE  ********************************
     *************************************************************************/  
    else if (command_flag == INVALID_INPUT_CODE ||
	     command_flag == INVALID_INPUT_VALUE) {
      line_number = -line_number;
      cout << "\n Chem2D ERROR: Error reading Chem2D data at line #"
	   << -line_number  << " of input data file.\n";
      cout.flush();
      return (line_number);
     } /* endif */
    
  } /* endwhile */

  /********************************************************  
   * End of all Chem2DSolver computations and I/O.       *
   ********************************************************/    
  //should use terminate code ie. should never get here?
  cout<<"\nEND OF CHEM2DSOLVER BUT DIDN'T USE TERMINATE!!\n"; cout.flush();
  
  return (0);
 
} //end Chem2DSolver



