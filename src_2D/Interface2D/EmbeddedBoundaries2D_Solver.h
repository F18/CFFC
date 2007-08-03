/**********************************************************************
 * EmbeddedBoundaries2D_Solver.h: 2D Embedded Boundaries multi-block  *
 *                                quadrilateral mesh solvers.         *
 **********************************************************************/

// Include the embedded boundaries header file.

#ifndef _EMBEDDEDBOUNDARIES2D_INCLUDED
#include "EmbeddedBoundaries2D.h"
#endif // _EMBEDDEDBOUNDARIES2D_INCLUDED

// Include 2D Level Set quadrilateral mesh solution header file.

#ifndef _LEVELSET2D_QUAD_INCLUDED
#include "../LevelSet2D/LevelSet2DQuad.h"
#endif // _LEVELSET2D_QUAD_INCLUDED

// Include the embedded boundaries multigrid header file.

#ifndef _EMBEDDEDBOUNDARIES2D_FASMULTIGRID2D_INCLUDED
#include "EmbeddedBoundaries2D_FASMultigrid.h"
#endif // _EMBEDDEDBOUNDARIES2D_FASMULTIGRID2D_INCLUDED

// Include the embedded boundaries solution header files.

#ifndef _EMBEDDEDBOUNDARIES2D_EULER_INCLUDED
#include "EmbeddedBoundaries2D_Euler.h"
#endif // _EMBEDDEDBOUNDARIES2D_EULER_INCLUDED

#ifndef _EMBEDDEDBOUNDARIES2D_NAVIERSTOKES_INCLUDED
#include "EmbeddedBoundaries2D_NavierStokes.h"
#endif // _EMBEDDEDBOUNDARIES2D_NAVIERSTOKES_INCLUDED

#ifndef _EMBEDDEDBOUNDARIES2D_DUSTY_INCLUDED
#include "EmbeddedBoundaries2D_Dusty.h"
#endif // _EMBEDDEDBOUNDARIES2D_DUSTY_INCLUDED

// #ifndef _EMBEDDEDBOUNDARIES2D_GAUSSIAN_INCLUDED
 #include "EmbeddedBoundaries2D_Gaussian.h"
// #endif // _EMBEDDEDBOUNDARIES2D_GAUSSIAN_INCLUDED

/**********************************************************************
 * Routine:EmbeddedBoundaries2D_Solver                                *
 *                                                                    *
 * Computes solutions to the templated equations on 2D quadrilateral  *
 * multi-block solution-adaptive mesh with embedded boundaries.       *
 *                                                                    *
 **********************************************************************/
template <class cState, class pState, class Quad_Soln_Block, class Quad_Soln_Input_Parameters>
int EmbeddedBoundaries2D_Solver(char *Input_File_Name_ptr, int batch_flag) {

  /********************************************************************
   * Local variable declarations.                                     *
   ********************************************************************/

  // Input variables and parameters.
  Quad_Soln_Input_Parameters    Input_Parameters;

  // Multi-block solution-adaptive quadrilateral mesh solution variables.
  Grid2D_Quad_Block           **MeshBlk;
  QuadTreeBlock_DataStructure   QuadTree;
  AdaptiveBlockResourceList     Global_Solution_Block_List;
  AdaptiveBlock2D_List          Local_Solution_Block_List;
  Quad_Soln_Block              *Local_SolnBlk;

  // LevelSet2D input parameters and multi-block solution-adaptive
  // quadrilateral mesh solution variables.
  LevelSet2D_Input_Parameters   LS_Input_Parameters;
  QuadTreeBlock_DataStructure   LS_QuadTree;
  AdaptiveBlockResourceList     LS_Global_Solution_Block_List;
  AdaptiveBlock2D_List          LS_Local_Solution_Block_List;
  LevelSet2D_Quad_Block        *LS_Local_SolnBlk;

  // Embedded boundaries solver declaration.
  EmbeddedBoundaries2D<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters> EBSolver;

  // Embedded boundaries FAS-multigrid declaration.
  EB_FAS_Multigrid2D_Solver<cState, pState, Quad_Soln_Block, Quad_Soln_Input_Parameters> MGSolver;

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

  // Level set solution variables.
  int levelset_iterations, evolution_counter;
  double levelset_Time;

  /********************************************************************  
   * Set default values for the input solution parameters and then    *
   * read user specified input values from the specified input        *
   * parameter file.                                                  *
   ********************************************************************/

  // The primary MPI processor processes the input parameter file.
  if (CFDkit_Primary_MPI_Processor()) {
    if (!batch_flag)
      cout << "\n Reading input data file `"
	   << Input_File_Name_ptr << "'.";
    error_flag = Process_Input_Control_Parameter_File(Input_Parameters,
						      Input_File_Name_ptr,
						      command_flag);
    if (!batch_flag && !error_flag) {
      cout << Input_Parameters << endl;
      cout.flush();
    } else if (error_flag) {
      cout << endl << " ERROR: During processing of input parameters."
	   << endl << endl;
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
   * Create initial mesh and allocate solution variables for the      *
   * specified set of PDES and IBVP/BVP problem .                     *
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
      error_flag = Output_Tecplot(MeshBlk,Input_Parameters);
      error_flag = Output_Nodes_Tecplot(MeshBlk,Input_Parameters);
      error_flag = 1;
      cout << "\n ERROR: Unable to create valid multi-block mesh.\n";
      cout.flush();
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
  // (allocate) array of local solution blocks, assign and create 
  // (allocate) solution blocks corresponding to the initial mesh.
  if (!batch_flag) cout << "\n Creating multi-block quadtree data structure and assigning"
                        << "\n solution blocks corresponding to initial mesh.";
  Local_SolnBlk = CreateInitialSolutionBlocks(MeshBlk,
					      Local_SolnBlk,
					      Input_Parameters,
					      QuadTree,
					      Global_Solution_Block_List,
					      Local_Solution_Block_List);
  if (Local_SolnBlk == NULL) return 1;

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //char extension[256], pout_name[256];
  //char *pout_name_ptr;
  //ofstream pout;
  //sprintf(extension,"%.6d",Local_Solution_Block_List.ThisCPU);
  //strcat(extension,".in");
  //strcpy(pout_name,"ip_cpu");
  //strcat(pout_name,extension);
  //pout_name_ptr = pout_name;
  //pout.open(pout_name_ptr,ios::out);
  //if (pout.bad()) return 1;
  //pout << Input_Parameters;
  //pout << Input_Parameters.Interface_IP.Component_List;
  //pout.close();
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////

  /********************************************************************
   * Initialize solution variables.                                   *
   ********************************************************************/

  // Set the initial time level.
  Time = ZERO;
  number_of_time_steps = 0;
  levelset_iterations = 0;
  evolution_counter = 1;

  // Set the CPU time to zero.
  processor_cpu_time.zero();
  total_cpu_time.zero();

  // Initialize the conserved and primitive state solution variables.
  if (!batch_flag) cout << "\n Prescribing initial data.";
  if (Input_Parameters.i_ICs == IC_RESTART) {
    if (!batch_flag) cout << "\n Reading solution from restart data files.";

    //Check that restart files are probably not corrupt. //added by ~james
    if (CFDkit_Primary_MPI_Processor()) {
      if(System::Restart_In_Progress()) {
	cout << "\n  Restart-in-progress flag detected, assuming data is corrupt."
	     << "\n  Uncompressing backups.";
	System::Uncompress_Restart();
	System::Remove_Restart_Flag();
	cout << "\n  Backup successfully uncompressed; reading.";
      }
    }
    CFDkit_Barrier_MPI(); // MPI barrier to ensure processor synchronization.

    // Read the quadtree restart file.
    error_flag = Read_QuadTree(QuadTree,
			       Global_Solution_Block_List,
			       Local_Solution_Block_List,
			       Input_Parameters);
    if (error_flag) {
      cout << "\n ERROR: Unable to open quadtree data file on processor "
	   << Local_Solution_Block_List.ThisCPU
	   << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;
    // Allocate the message buffers.
    Allocate_Message_Buffers(Local_Solution_Block_List,
			     Local_SolnBlk[0].NumVar()+NUM_COMP_VECTOR2D);
    // Read the solution block restart files.
    error_flag = Read_Restart_Solution(Local_SolnBlk,
				       Local_Solution_Block_List,
				       Input_Parameters,
				       number_of_time_steps,
				       Time,
				       processor_cpu_time);
    if (error_flag) {
      cout << endl << " ERROR: Unable to open restart input data "
	   << "file(s) on processor "
	   << Local_Solution_Block_List.ThisCPU << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
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
    if (error_flag) return error_flag;
    CFDkit_Broadcast_MPI(&command_flag,1);
    if (command_flag == TERMINATE_CODE) return 0;
    Broadcast_Input_Parameters(Input_Parameters);

  }

  // Send solution information between neighbouring blocks to complete
  // prescription of initial data.
  error_flag = Send_All_Messages(Local_SolnBlk,
				 Local_Solution_Block_List,
				 NUM_COMP_VECTOR2D,
				 ON);
  if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk,
						  Local_Solution_Block_List,
						  Local_SolnBlk[0].NumVar(),
						  OFF);
  if (error_flag) {
    cout << "\n ERROR: Message passing error during solution initialization on "
	 << "processor " << Local_Solution_Block_List.ThisCPU << "." << endl;
  }
  error_flag = CFDkit_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Allocate memory for the embedded boundary solution.
  EBSolver.allocate(Local_SolnBlk,
		    &QuadTree,
		    &Global_Solution_Block_List,
		    &Local_Solution_Block_List,
		    &Input_Parameters,
		    LS_Local_SolnBlk,
		    &LS_QuadTree,
		    &LS_Global_Solution_Block_List,
		    &LS_Local_Solution_Block_List,
		    &LS_Input_Parameters);

  // Complete prescription of initial data.
  if (Input_Parameters.i_ICs == IC_RESTART) {

    // Read embedded boundary restart files.
    error_flag = EBSolver.Read_Restart_Files(batch_flag,
					     levelset_iterations,
					     levelset_Time);
    if (error_flag) {
      cout << endl << " ERROR: Unable to process embedded boundary process restart input "
	   << endl << " data file(s) on processor " << Local_Solution_Block_List.ThisCPU << "."
	   << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    evolution_counter = number_of_time_steps%Input_Parameters.Interface_IP.Evolution_Frequency;
    if (!evolution_counter) evolution_counter = 1;

  } else {

    // Construct the interface component list.
    if (!batch_flag) cout << endl << " Constructing interface component list.";
    error_flag = EBSolver.Construct_Interface_Component_List();
    if (error_flag) {
      cout << endl << " ERROR: Unable to construct interface component list on processor "
	   << Local_Solution_Block_List.ThisCPU << ".  Error number = "
	   << error_flag << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;
    // Construct the interface level set list.
    error_flag = EBSolver.Construct_Interface_Level_Set_List(batch_flag);
    if (error_flag) {
      cout << endl << " ERROR: Unable to construct level set interface list on processor "
	   << Local_Solution_Block_List.ThisCPU << ".  Error number = "
	   << error_flag << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;
    // Construct the interface union list.
    if (!batch_flag) cout << endl << " Constructing interface union list.";
    error_flag = EBSolver.Construct_Interface_Union_List();
    if (error_flag) {
      cout << endl << " ERROR: Unable to construct interface union list on processor "
	   << Local_Solution_Block_List.ThisCPU << ".  Error number = "
	   << error_flag << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;
    // Synchronize processors.
    CFDkit_Barrier_MPI();
    // Initialize grids and mesh data required for the adjustment
    // algorithm.
    EBSolver.Initialize_Adjustment_Grids();
    // Perform mesh adjustment according to interface location(s).
    if (!batch_flag) cout << "\n Performing initial mesh adjustment.";
    error_flag = EBSolver.Mesh_Adjustment(ON,ON);
    if (error_flag) {
      cout << endl << " ERROR: Unable to adjust grid on processor "
	   << Local_Solution_Block_List.ThisCPU << ".  Error number = " 
	   << error_flag << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;
    // Apply initial conditions.
    if (!batch_flag) cout << "\n Applying initial conditions.";
    ICs(Local_SolnBlk,Local_Solution_Block_List,Input_Parameters);
    // Determine the interface velocity function on both interface lists.
    if (!batch_flag) cout << "\n Determining the interface velocity function.";
    error_flag = EBSolver.Compute_Interface_Velocity_Function(Time);
    if (error_flag) {
      cout << endl << " ERROR: Unable to determine the interface velocity function on processor "
	   << Local_Solution_Block_List.ThisCPU << ".  Error number = " 
	   << error_flag << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;
  }

  // MPI barrier to ensure processor synchronization.
  CFDkit_Barrier_MPI();

  // Reset the interface motion type if required (used for oscillaring
  // NACA0012 aerofoils).
  if (Input_Parameters.Reset_Interface_Motion_Type) {
    error_flag = EBSolver.Reset_Interface_Motion_Type();
    Time = ZERO;
    if (!error_flag) error_flag = EBSolver.Compute_Interface_Velocity_Function(Time);
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;
  }

  // Send solution information between neighbouring blocks to complete
  // prescription of initial data.
  error_flag = Send_All_Messages(Local_SolnBlk,
				 Local_Solution_Block_List,
				 NUM_COMP_VECTOR2D,
				 ON);
  if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk,
						  Local_Solution_Block_List,
						  Local_SolnBlk[0].NumVar(),
						  OFF);
  if (error_flag) {
    cout << "\n ERROR: Message passing error during solution initialization on "
	 << "processor " << Local_Solution_Block_List.ThisCPU << "." << endl;
  }
  error_flag = CFDkit_OR_MPI(error_flag);
  if (error_flag) return error_flag;

  // Prescribe boundary data consistent with initial data.
  EBSolver.Boundary_Conditions(Time);

  // Perform uniform, boundary, and initial mesh refinement.
  if (Input_Parameters.i_ICs != IC_RESTART) {

    // Perform uniform mesh refinement.
    if (!batch_flag) cout << "\n Performing uniform mesh refinement.";
    error_flag = EBSolver.Uniform_Adaptive_Mesh_Refinement();
    if (error_flag) {
      cout << "\n ERROR: Uniform AMR error on processor "
	   << Local_Solution_Block_List.ThisCPU << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Perform boundary mesh refinement.
    if (!batch_flag) cout << "\n Performing boundary mesh refinement.";
    error_flag = EBSolver.Boundary_Adaptive_Mesh_Refinement();
    if (error_flag) {
      cout << "\n ERROR: Boundary AMR error on processor "
	   << Local_Solution_Block_List.ThisCPU << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Perform interface mesh refinement.
    if (!batch_flag) cout << "\n Performing interface mesh refinement.";
    error_flag = EBSolver.Interface_Adaptive_Mesh_Refinement();
    if (error_flag) {
      cout << "\n ERROR: Interface AMR error #" << error_flag << " on processor "
	   << Local_Solution_Block_List.ThisCPU << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Perform bounding-box mesh refinement.
    if (!batch_flag) cout << "\n Performing bounding-box mesh refinement.";
    error_flag = EBSolver.Bounding_Box_Adaptive_Mesh_Refinement(ON);
    if (error_flag) {
      cout << "\n ERROR: Bounding-box AMR error #" << error_flag << " on processor "
	   << Local_Solution_Block_List.ThisCPU << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Perform initial mesh refinement.
    if (!batch_flag) cout << "\n Performing initial mesh refinement.";
    error_flag = EBSolver.Initial_Adaptive_Mesh_Refinement();
    if (error_flag) {
      cout << "\n ERROR: Initial AMR error on processor "
	   << Local_Solution_Block_List.ThisCPU << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    //re-apply ICs after refinement
    if (Input_Parameters.Number_of_Uniform_Mesh_Refinements ||
	Input_Parameters.Number_of_Boundary_Mesh_Refinements ||
	Input_Parameters.Number_of_Interface_Mesh_Refinements ||
	Input_Parameters.Number_of_Initial_Mesh_Refinements ||
	Input_Parameters.Number_of_Bounding_Box_Mesh_Refinements) {
      if (!batch_flag) cout << "\n Re-applying ICs after mesh refinement.";
      // Apply initial conditions.
      ICs(Local_SolnBlk,Local_Solution_Block_List,Input_Parameters);
    }

    // MPI barrier to ensure processor synchronization.
    CFDkit_Barrier_MPI();

    // Send solution information between neighbouring blocks to complete
    // prescription of initial data.
    error_flag = Send_All_Messages(Local_SolnBlk,
				   Local_Solution_Block_List,
				   NUM_COMP_VECTOR2D,
				   ON);
    if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk,
						    Local_Solution_Block_List,
						    Local_SolnBlk[0].NumVar(),
						    OFF);
    if (error_flag) {
      cout << "\n ERROR: Message passing error during solution intialization "
	   << "on processor " << Local_Solution_Block_List.ThisCPU << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;

    // Determine the interface velocity function on both interface lists.
    if (!batch_flag) cout << "\n Redetermining the interface velocity function.";
    error_flag = EBSolver.Compute_Interface_Velocity_Function(Time);
    if (error_flag) {
      cout << endl << " ERROR: Unable to determine the interface velocity function on processor "
	   << Local_Solution_Block_List.ThisCPU << ".  Error number = " 
	   << error_flag << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;

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

//   if (CFDkit_Primary_MPI_Processor()) {
//     for (int j_blk = 0; j_blk < QuadTree.Nblk; j_blk++) {
//       for (int i_blk = 0; i_blk < QuadTree.Ncpu; i_blk++) {
// 	if (QuadTree.Blocks[i_blk][j_blk] != NULL) {
// 	  cout << "\n cpu = " << i_blk
// 	       << " blk = " << j_blk
// 	       << " blk = " << QuadTree.Blocks[i_blk][j_blk]->block;
// 	} else {
// 	  cout << "\n cpu = " << i_blk
// 	       << " blk = " << j_blk;
// 	}
//       }
//     }
//   }

  /********************************************************************
   * Solve IBVP or BVP for conservation form of the equations on      *
   * multi-block solution-adaptive quadrilateral mesh.                *
   ********************************************************************/

 continue_existing_calculation: ;

  // MPI barrier to ensure processor synchronization.
  CFDkit_Barrier_MPI();

  // Allocate memory for multigrid solver if required.

        //Note: this has never been tested for Gaussian2D!!!!
  if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
      Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
    error_flag = MGSolver.allocate(Local_SolnBlk,
				   &QuadTree,
				   &Global_Solution_Block_List,
				   &Local_Solution_Block_List,
				   &Input_Parameters,
				   LS_Local_SolnBlk,
				   &LS_QuadTree,
				   &LS_Global_Solution_Block_List,
				   &LS_Local_Solution_Block_List,
				   &LS_Input_Parameters,
				   &EBSolver);
    if (error_flag) {
      cout << "\n ERROR: Unable to allocate memory for multigrid solver."
	   << "\n Error number = " << error_flag << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;
  }

  // Execute required solution algorithm.
  if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {
    // Perform FAS multigrid computations.
    error_flag = MGSolver.Execute(batch_flag,
				  number_of_time_steps,
				  evolution_counter,
				  levelset_iterations,
				  Time,
				  levelset_Time,
				  processor_cpu_time,
				  total_cpu_time,
				  residual_file);
    if (error_flag) {
      cout << "\n ERROR: Multigrid error on processor "
	   << Local_Solution_Block_List.ThisCPU << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;

  } else if (Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
    // Perform dual-time computations with FAS multigrid.

       //Note: This has never been tested for Gaussian2D!!!
    error_flag = MGSolver.DTS_Multigrid_Solution(batch_flag,
						 number_of_time_steps,
						 evolution_counter,
						 levelset_iterations,
						 Time,
						 levelset_Time,
						 processor_cpu_time,
						 total_cpu_time,
						 residual_file);
    if (error_flag) {
      cout << "\n ERROR: Error during DTS multigrid solution on processor "
	   << Local_Solution_Block_List.ThisCPU << ".  Error number = "
	   << error_flag << "." << endl;
    }
    CFDkit_Broadcast_MPI(&error_flag,1);
    if (error_flag) return error_flag;

  } else {
    // Perform explicit time-stepping computations.
    error_flag = EBSolver.Execute(batch_flag,
				  number_of_time_steps,
				  evolution_counter,
				  levelset_iterations,
				  Time,
				  levelset_Time,
				  processor_cpu_time,
				  total_cpu_time,
				  residual_file);
    if (error_flag) {
      cout << "\n ERROR: Embedded boundaries error on processor "
	   << Local_Solution_Block_List.ThisCPU << "." << endl;
    }
    error_flag = CFDkit_OR_MPI(error_flag);
    if (error_flag) return error_flag;

  }

  /********************************************************************
   * Solution calculations complete.  Write the solution to output    *
   * and restart files as required, reset solution parameters, and    *
   * run other cases specified by input parameters.                   *
   ********************************************************************/

 postprocess_current_calculation: ;

  // MPI barrier to ensure processor synchronization.
  CFDkit_Barrier_MPI();

  while (1) {
    if (CFDkit_Primary_MPI_Processor()) {
      Get_Next_Input_Control_Parameter(Input_Parameters);
      command_flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
      line_number = Input_Parameters.Line_Number;
//       if (command_flag == INVALID_INPUT_CODE ||
// 	  command_flag == INVALID_INPUT_VALUE) {
// 	line_number = -line_number;
// 	cout << "\n ERROR: Error reading data at line #"
// 	     << -line_number << " of input data file.\n";
// 	return 1;
//       }
      //Reinitialize_Reference_State(Input_Parameters);
    }
    // MPI barrier to ensure processor synchronization.
    CFDkit_Barrier_MPI();
    Broadcast_Input_Parameters(Input_Parameters);
    CFDkit_Broadcast_MPI(&command_flag,1);

    if (command_flag == EXECUTE_CODE) {
      // Deallocate memory for the solution variables.
      if (!batch_flag) cout << "\n Deallocating solution variables.";
      if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	  Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
      }
      Local_SolnBlk = Deallocate(Local_SolnBlk,Input_Parameters);
      //Deallocate_Message_Buffers(Local_Solution_Block_List); // Not necessary here!
      Local_Solution_Block_List.deallocate();
      Global_Solution_Block_List.deallocate();
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

      // Deallocate memory for equation solution.
      if (!batch_flag) cout << "\n Deallocating solution variables.";
      if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	  Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING) {
      }
      if (Input_Parameters.Interface_IP.Number_of_Components) EBSolver.deallocate();
      Local_SolnBlk = Deallocate(Local_SolnBlk,Input_Parameters);
      //Deallocate_Message_Buffers(Local_Solution_Block_List); // Not necessary here!
      Local_Solution_Block_List.deallocate();
      Global_Solution_Block_List.deallocate();
      QuadTree.deallocate();
      MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk,
					    Input_Parameters.Number_of_Blocks_Idir,
					    Input_Parameters.Number_of_Blocks_Jdir);
      // Close input data file.
      if (!batch_flag) cout << "\n\n Closing input data file.";
      if (CFDkit_Primary_MPI_Processor()) Close_Input_File(Input_Parameters);
      // Terminate calculation.
      return 0;

    } else if (command_flag == CONTINUE_CODE) {
      // Reset maximum time step counter.
      Input_Parameters.Maximum_Number_of_Time_Steps += number_of_time_steps;
      // Output input parameters for continuing calculation.
      if (!batch_flag) cout << "\n\n Continuing existing calculation."
			    << Input_Parameters << endl;
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
      Evaluate_Limiters(Local_SolnBlk,Local_Solution_Block_List);
      error_flag = EBSolver.Adaptive_Mesh_Refinement(ON,ON);
      if (error_flag) {
	cout << "\n ERROR: AMR error on processor "
	     << Local_Solution_Block_List.ThisCPU
	     << ".  Error number = " << error_flag << "." << endl;
      }
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) {
	command_flag = Output_Tecplot(Local_SolnBlk,
				      Local_Solution_Block_List,
				      Input_Parameters,
				      number_of_time_steps,
				      Time);
	return error_flag;
      }
      error_flag = EBSolver.Boundary_Conditions(Time);
      if (error_flag) {
	cout << "\n ERROR: Boundary conditions error on processor "
	     << Local_Solution_Block_List.ThisCPU << "." << endl;
      }
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return error_flag;

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

    } else if (command_flag == BOUNDING_BOX_REFINE_GRID_CODE) {
      // Divide the mesh blocks that fall within the specified bounding box.
      if (!batch_flag) cout << "\n Refining Grid.  Performing bounding-box adaptive mesh refinement.";
      Evaluate_Limiters(Local_SolnBlk,Local_Solution_Block_List);
      error_flag = EBSolver.Bounding_Box_Adaptive_Mesh_Refinement(OFF);
      if (error_flag) {
	cout << "\n ERROR: Bounding-Box AMR error on processor "
	     << Local_Solution_Block_List.ThisCPU
	     << ".  Error number = " << error_flag << "." << endl;
      }
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) {
	command_flag = Output_Tecplot(Local_SolnBlk,
				      Local_Solution_Block_List,
				      Input_Parameters,
				      number_of_time_steps,
				      Time);
	return error_flag;
      }
      error_flag = EBSolver.Boundary_Conditions(Time);
      if (error_flag) {
	cout << "\n ERROR: Boundary conditions error on processor "
	     << Local_Solution_Block_List.ThisCPU << "." << endl;
      }
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return error_flag;

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

    } else if (command_flag == WRITE_OUTPUT_CODE) {
      // Output solution data.
      if (!batch_flag) cout << "\n Writing solution to output data file(s).";
      if (!(Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	    Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING)) {
	error_flag = Output_Tecplot(Local_SolnBlk,
				    Local_Solution_Block_List,
				    Input_Parameters,
				    number_of_time_steps,
				    Time);
      } else {
  	error_flag = MGSolver.Output_Multigrid(number_of_time_steps,
  					       Time);
      }
      if (error_flag) {
	cout << "\n ERROR: Unable to open output data file(s) on processor "
	     << Local_Solution_Block_List.ThisCPU << "." << endl;
      }
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_CELLS_CODE) {
      // Output solution data.
      if (!batch_flag) cout << "\n Writing cell-centered solution to output data file(s).";
      if (!(Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	    Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING)) {
	error_flag = Output_Cells_Tecplot(Local_SolnBlk,
					  Local_Solution_Block_List,
					  Input_Parameters,
					  number_of_time_steps,
					  Time);
      } else {
  	error_flag = MGSolver.Output_Multigrid_Cells(number_of_time_steps,
  						     Time);
      }
      if (error_flag) {
	cout << "\n ERROR: Unable to open cell output data file(s) on processor "
	     << Local_Solution_Block_List.ThisCPU << "." << endl;
      }
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_ELEMENTS_CODE) {
      if (!batch_flag) cout << "\n Writing solution elements to output data files.";
      if (!(Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	    Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING)) {
	error_flag = EBSolver.Output_Elements_Tecplot(number_of_time_steps,
						      Time);
      } else {
  	error_flag = MGSolver.Output_Multigrid_Elements(number_of_time_steps,
  							Time);
      }
      if (error_flag) {
	cout << "\n ERROR: Unable to open solution element output file.\n";
	cout.flush();
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_NODES_CODE) {
      if (!batch_flag) cout << "\n Writing solution nodes to output data files.";
      if (!(Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	    Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING)) {
	error_flag = EBSolver.Output_Nodes_Tecplot();
      } else {
  	error_flag = MGSolver.Output_Multigrid_Nodes(number_of_time_steps,
  						     Time);
      }
      if (error_flag) {
	cout << "\n ERROR: Unable to open solution nodes output file." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_GRADIENTS_CODE) {
      if (!batch_flag) cout << "\n Writing solution gradients to output data files.";
      error_flag = Output_Gradients_Tecplot(Local_SolnBlk,
					    Local_Solution_Block_List,
					    Input_Parameters,
					    number_of_time_steps,
					    Time);
      if (error_flag) {
	cout << "\n ERROR: Unable to open solution gradients output file." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_CELL_STATUS_CODE) {
      if (!batch_flag) cout << "\n Writing cell status data to output data files.";
      if (!(Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	    Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING)) {
	error_flag = EBSolver.Output_Cell_Status_Tecplot();
      } else {
  	error_flag = MGSolver.Output_Multigrid_Cell_Status(number_of_time_steps,
  							   Time);
      }
      if (error_flag) {
	cout << "\n ERROR: Unable to open cell status data output file." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_INTERFACE_COMPONENT_LIST_CODE) {
      if (!batch_flag) cout << "\n Writing interface component list to output file.";
      if (!(Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	    Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING)) {
	error_flag = EBSolver.Output_Interface_Component_List_Tecplot();
      } else {
 	error_flag = MGSolver.EBSolver[0].Output_Interface_Component_List_Tecplot();
      }
      if (error_flag) {
	cout << "\n ERROR: Unable to open interface component list output file." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_INTERFACE_UNION_LIST_CODE) {
      if (!batch_flag) cout << "\n Writing interface union list to output file.";
      if (!(Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	    Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING)) {
	error_flag = EBSolver.Output_Interface_Union_List_Tecplot();
      } else {
 	error_flag = MGSolver.EBSolver[0].Output_Interface_Union_List_Tecplot();
      }
      if (error_flag) {
	cout << "\n ERROR: Unable to open interface union list output file." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_LEVEL_SET_CODE) { 
      if (!batch_flag) cout << "\n Writing LevelSet2D solution to output data file(s).";
      error_flag = EBSolver.Output_Level_Set_Tecplot(number_of_time_steps,
						     Time);
      if (error_flag) {
	cout << "\n ERROR: Unable to open LevelSet2D output data file(s)." << endl;
      }
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_LEVEL_SET_CELLS_CODE) {
      if (!batch_flag) cout << "\n Writing cell-centered LevelSet2D solution to output data file(s).";
      error_flag = EBSolver.Output_Level_Set_Cells_Tecplot(number_of_time_steps,
							   Time);
      if (error_flag) {
	cout << "\n ERROR: Unable to open LevelSet2D cell output data file(s)." << endl;
      }
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_LEVEL_SET_INTERFACE_LIST_CODE) {
      if (!batch_flag) cout << "\n Writing LevelSet2D interface nodes to output data file(s).";
      error_flag = EBSolver.Output_Level_Set_Interface_Tecplot(number_of_time_steps,
							       Time);
      if (error_flag) {
	cout << "\n ERROR: Unable to open LevelSet2D interface nodes output data file(s)." << endl;
      }
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_RINGLEB_CODE) {
      if (!batch_flag) cout << endl << " Writing exact solution and error norms for Ringleb's flow.";
      if (!(Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	    Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING)) {
	error_flag = EBSolver.Output_Ringleb_Tecplot();
      } else {
 	error_flag = MGSolver.EBSolver[0].Output_Ringleb_Tecplot();
      }
      if (error_flag) {
	cout << endl << "\n ERROR: Unable to open Ringleb's flow output file." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_FLAT_PLATE_CODE) {
      if (!batch_flag) cout << endl << " Writing exact solution and error norms for the flat plate flow (Blasius solution).";
      if (!(Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	    Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING)) {
	error_flag = EBSolver.Output_Flat_Plate_Tecplot();
      } else {
 	error_flag = MGSolver.EBSolver[0].Output_Flat_Plate_Tecplot();
      }
      if (error_flag) {
	cout << endl << "\n ERROR: Unable to open flat plate output file." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_CYLINDER_DRAG_CODE) {
      if (!batch_flag) cout << endl << " Calculating Cylinder Cd and Cl.";
      if (!(Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	    Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING)) {
	error_flag = EBSolver.Output_Cylinder_Drag();
      } else {
	// 	error_flag = MGSolver.EBSolver[0].Output_Flat_Plate_Tecplot();
      }
      if (error_flag) {
	cout << endl << "\n ERROR: Problem calculating Cd and Cl." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_COUETTE) {
      if (!batch_flag) cout << endl << " Calculating Couette info.";
      if (!(Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	    Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING)) {
	error_flag = EBSolver.Output_Couette();
      } else {
	// 	error_flag = MGSolver.EBSolver[0].Output_Flat_Plate_Tecplot();
      }
      if (error_flag) {
	cout << endl << "\n ERROR: Problem calculating Couette info." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_OUTPUT_AERODYNAMIC_COEFFICIENTS_CODE) {
      if (!batch_flag) cout << endl << " Writing the aerodynamic coefficients to an output data file.";
      if (!(Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	    Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING)) {
	error_flag = EBSolver.Output_Aerodynamic_Coefficients_Tecplot(number_of_time_steps,
								      Time);
      } else {
 	error_flag = MGSolver.EBSolver[0].Output_Aerodynamic_Coefficients_Tecplot(number_of_time_steps,
										  Time);
      }
      if (error_flag) {
	cout << endl << "\n ERROR: Unable to open aerodynamic coefficients output file." << endl;
      }
      CFDkit_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;

    } else if (command_flag == WRITE_RESTART_CODE) {
      // Write restart files.
      if (!batch_flag) cout << "\n Writing solution to restart data file(s).";
      //  Save and delete old restart files in compressed archive (just in case)
      if (CFDkit_Primary_MPI_Processor()) {
	cout << "\n  Creating compressed archive of (and deleting) old restarts.";
	System::Compress_Restart();
	cout << "\n  Writing new restart files.";
	cout.flush();
      }
      CFDkit_Barrier_MPI(); // MPI barrier so that other processors do
                            // not start over writing restarts

      if (CFDkit_Primary_MPI_Processor()) {
	System::Set_Restart_Flag();  //Set flag to indicate a restart is being saved
      }

      if (!batch_flag) cout << "\n Writing solution to restart data file(s).";
      // Write the quadtree restart file.
      error_flag = Write_QuadTree(QuadTree,Input_Parameters);
      if (error_flag) {
	cout << "\n ERROR: Unable to open quadtree data file on processor "
	     << Local_Solution_Block_List.ThisCPU << "." << endl;
      }
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return error_flag;
      // Write the solution block restart files.
      error_flag = Write_Restart_Solution(Local_SolnBlk,
					  Local_Solution_Block_List,
					  Input_Parameters,
					  number_of_time_steps,
					  Time,
					  processor_cpu_time);
      if (error_flag) {
	cout << "\n ERROR: Unable to open restart output data file(s) on processor "
	     << Local_Solution_Block_List.ThisCPU << "." << endl;
      }
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return error_flag;
      // Write out data files for the embedded boundaries.
      if (!(Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID ||
	    Input_Parameters.i_Time_Integration == TIME_STEPPING_DUAL_TIME_STEPPING)) {
	error_flag = EBSolver.Write_Restart_Files(number_of_time_steps,
						  levelset_iterations,
						  Time,
						  levelset_Time,
						  processor_cpu_time);
      } else {
 	error_flag = MGSolver.EBSolver[0].Write_Restart_Files(number_of_time_steps,
 							      levelset_iterations,
 							      Time,
							      levelset_Time,
 							      processor_cpu_time);
      }
      if (error_flag) {
	cout << "\n ERROR: Unable to open restart output data file(s) on processor "
	     << Local_Solution_Block_List.ThisCPU << "." << endl;
      }
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return error_flag;

       if (CFDkit_Primary_MPI_Processor()) {
	 System::Remove_Restart_Flag();  //Remove flag to indicate the restart is finished
       }

    } else if (command_flag == INVALID_INPUT_CODE ||
	       command_flag == INVALID_INPUT_VALUE) {
      line_number = -line_number;
      cout << "\n ERROR: Error reading data at line #"
	   << -line_number  << " of input data file." << endl;
      return line_number;
    }
    
  }
  
  /********************************************************************
   * End of all EmbeddedBoundaries2D_Solver computations and I/O.     *
   ********************************************************************/
  return 0;
  
}
