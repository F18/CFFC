/* Rte2DQuadSolvers.h:  2D Rte Equation 
                        Multi-Block Quadrilateral Mesh Solvers. */

#ifndef _RTE2D_QUADSOLVER__INCLUDED
#define _RTE2D_QUADSOLVER__INCLUDED


/* Include 2D Rte quadrilateral mesh solution header file. */
#include "Rte2DQuad.h"

/* Include Rte2D Multigrid Specialization header file. */
#include "Rte2DQuadMultigrid.h"

/* Include Rte2D Newton-Krylov-Schwarz Specialization header file. */
#include "Rte2DQuadNKS.h"

/* Include Rte2D AMR Specialization header file. */
#include "Rte2DQuadAMR.h"


/****************************************************************************
 ****************************************************************************
  Class: Rte2DQuadSolver

  The old-style 2DQuadSolver() function was modularized 
  into several parts so that one can easily play with the
  order the modules are called.  This also helps make
  the overall process more 'readable'.

  The class itself holds the same local variables as 
  before, only now they are members.  This way, you
  don't have to worry about passing everything everywhere.
 ****************************************************************************
 ****************************************************************************/
class Rte2DSolver {
  
 public:

  /********************************************************  
   * Object Declarations                                  *
   ********************************************************/   
  // Rte2D input variables and parameters:
  Rte2D_Input_Parameters Input_Parameters;


  //
  // Multi-block solution-adaptive quadrilateral mesh 
  // solution variables.
  //

  // For Rte2D State (beam transport)
  Grid2D_Quad_Block          **MeshBlk;
  QuadTreeBlock_DataStructure  QuadTree;
  AdaptiveBlockResourceList    List_of_Global_Solution_Blocks;
  AdaptiveBlock2D_List         List_of_Local_Solution_Blocks;
  Rte2D_Quad_Block            *Local_SolnBlk;

  // For Medium2D Field (participating medium description)
  Grid2D_Quad_Block          **MF_MeshBlk;
  QuadTreeBlock_DataStructure  MF_QuadTree;
  AdaptiveBlockResourceList    MF_List_of_Global_Solution_Blocks;
  AdaptiveBlock2D_List         MF_List_of_Local_Solution_Blocks;
  Rte2D_Quad_Block            *MF_Local_SolnBlk;


  // Multigrid declaration.
  FAS_Multigrid2D_Solver<Rte2D_State, 
                         Rte2D_Quad_Block, 
                         Rte2D_Input_Parameters> MGSolver;


  // cpu time variables.
  CPUTime processor_cpu_time, total_cpu_time;
  time_t start_explicit, end_explicit;


  // Other local solution variables.
  int number_of_time_steps;
  bool first_step;
  double Time;
  int NUM_VAR_RTE2D;

  /********************************************************  
   * Constructors/Destructors                             *
   ********************************************************/
  Rte2DSolver () : Input_Parameters(), 
		   MeshBlk(NULL), QuadTree(),
                   List_of_Global_Solution_Blocks(),
		   List_of_Local_Solution_Blocks(),
		   Local_SolnBlk(NULL), 
		   MF_MeshBlk(NULL), MF_QuadTree(),
                   MF_List_of_Global_Solution_Blocks(),
		   MF_List_of_Local_Solution_Blocks(),
		   MF_Local_SolnBlk(NULL), 
		   MGSolver(),
		   processor_cpu_time(), total_cpu_time(),
		   start_explicit(), end_explicit(), 
		   number_of_time_steps(0), Time(0),
		   first_step(true), NUM_VAR_RTE2D(0) {}
  //~Rte2DSolver() { DeallocateSoln(); }

  /********************************************************  
   * Member functions                                     *
   *                                                      *
   * The overal solution process has been broken into     *
   * several modular components.  This allows one to      *
   * enter/exit the solution process whenever they want.  *
   ********************************************************/
  void DeallocateSoln();
  int ReadInputFile(char *Input_File_Name_ptr, 
		    const int batch_flag,
		    int &command_flag);
  int SetupMesh(const int batch_flag);
  int SetICs(const int batch_flag);
  int BroadcastSolution();
  int PerformAMR(const int batch_flag);
  int PerformPeriodicAMR(const int batch_flag);
  int PerformMortonOrdering(const int batch_flag);
  int PerformMortonReOrdering(const int batch_flag);
  void ComputeNorms(double &residual_l1_norm, 
		    double &residual_l2_norm, 
		    double &residual_max_norm);
  int WritePeriodicRestart(const int batch_flag);
  int OpenResidualFile(ofstream &);
  void OutputProgress(const int batch_flag, 
		      const double residual_l1_norm, 
		      const double residual_l2_norm, 
		      const double residual_max_norm,
		      ofstream &);
  int SolveMultigrid(const int batch_flag);
  int SolveExplicitTimeStep(const int batch_flag);
  int SolveSpaceMarch(const int batch_flag);
  int SolveNKS(const int batch_flag);
  int PostProcess(const int batch_flag,int &command_flag);
  int StandAloneSolve(const int batch_flag, char *Input_File_Name_ptr);

};


/********************************************************  
 * Rte2DSolver :: Deallocate                            *
 *                                                      *
 * Deallocate all solution variables and objects.       *
 ********************************************************/
inline void Rte2DSolver::DeallocateSoln() {

  // multigrid
  if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {
    MGSolver.deallocate();
  }

  // local solution blocks
  Local_SolnBlk = Deallocate(Local_SolnBlk, 
			     Input_Parameters);
  
  // message buffers (Not necessary here!)
  //Deallocate_Message_Buffers(List_of_Local_Solution_Blocks);
  
  // solution block lists, quadtree, meshblock
  List_of_Local_Solution_Blocks.deallocate();
  List_of_Global_Solution_Blocks.deallocate();
  QuadTree.deallocate();
  MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
					Input_Parameters.Number_of_Blocks_Idir, 
					Input_Parameters.Number_of_Blocks_Jdir);
  
}
  
/********************************************************  
 * Rte2DSolver :: StandAloneSolve                       *
 *                                                      *
 * The main control function which performs the stand-  *
 * alone version of the solve.                          *
 ********************************************************/
inline int Rte2DSolver::StandAloneSolve(const int batch_flag, 
					char *Input_File_Name_ptr) {

  // Declares
  int error_flag;
  int command_flag(EXECUTE_CODE);

  //--------------------------------------------------
  // Setup default input parameters
  //--------------------------------------------------
  error_flag = ReadInputFile(Input_File_Name_ptr,
			     batch_flag,command_flag);
  if (error_flag) {
    cout << "\n Rte2D ERROR: Unable to set default inputs.\n"
	 << flush;
    return error_flag;
  } // endif

  // if we want to quit, then exit
  if (command_flag == TERMINATE_CODE) return (error_flag);

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // We are redirected here when we want to execute a
  // new calculation.
  //
  execute_new_calculation:
  //
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  //--------------------------------------------------
  // Setup mesh
  //--------------------------------------------------
  error_flag = SetupMesh(batch_flag);
  if (error_flag) {
    cout << "\n Rte2D ERROR: Unable to setup mesh.\n"
	 << flush;
    return error_flag;
  } // endif


  //--------------------------------------------------
  // Initialize Rte2D solution variables.
  //--------------------------------------------------

  // Set the initial time level.
  Time = ZERO;
  number_of_time_steps = 0;

  // Set the CPU time to zero.
  processor_cpu_time.zero();
  total_cpu_time.zero();


  //--------------------------------------------------
  // Setup ICs
  //--------------------------------------------------
  error_flag = SetICs(batch_flag);
  if (error_flag) {
    cout << "\n Rte2D ERROR: Unable to setup ICs.\n"
	 << flush;
    return error_flag;
  } // endif


  //--------------------------------------------------
  // Broadcast Solution Information
  //--------------------------------------------------
  error_flag = BroadcastSolution();
  if (error_flag) {
    cout << "\n Rte2D ERROR: Unable to broadcast solution.\n"
	 << flush;
    return error_flag;
  } // endif


  //--------------------------------------------------
  // Boundary Conditions
  //--------------------------------------------------
  //Prescribe boundary data consistent with initial data.
  BCs(Local_SolnBlk, List_of_Local_Solution_Blocks, Input_Parameters);


  //--------------------------------------------------
  // AMR
  //--------------------------------------------------
  error_flag = PerformAMR(batch_flag);
  if (error_flag) {
    cout << "\n Rte2D ERROR: Unable to perform AMR.\n"
	 << flush;
    return error_flag;
  } // endif


  //--------------------------------------------------
  // Morton Ordering
  //--------------------------------------------------
  error_flag = PerformMortonOrdering(batch_flag);
  if (error_flag) {
    cout << "\n Rte2D ERROR: Unable to perform Morton Ordering.\n"
	 << flush;
    return error_flag;
  } // endif
  

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // We are redirected here when we want to continue an
  // exisiting calculation.
  //
  continue_existing_calculation: 
  //
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();


  //--------------------------------------------------
  // Solve
  //--------------------------------------------------
  //
  // Multigrid
  //
  if (Input_Parameters.i_Time_Integration == TIME_STEPPING_MULTIGRID) {
    error_flag = SolveMultigrid(batch_flag);
    if (error_flag) {
      cout << "\n Rte2D ERROR: Unable to perform multigrid computations.\n"
	   << flush;
      return error_flag;
    } // endif

  //
  // Space marching
  //
  } else if (Input_Parameters.i_Time_Integration == TIME_STEPPING_SPACE_MARCH) {
    error_flag = SolveSpaceMarch(batch_flag);
    if (error_flag) {
      cout << "\n Rte2D ERROR: Unable to perform space-marching computations.\n"
	   << flush;
      return error_flag;
    } // endif

  //
  // Non-multigrid, explicit time marching
  //
  } else {
    error_flag = SolveExplicitTimeStep(batch_flag);
    if (error_flag) {
      cout << "\n Rte2D ERROR: Unable to perform explicit time-stepping computations.\n"
	   << flush;
      return error_flag;
    } // endif

  } // endif
  
  //
  // NKS time marching
  //
  if (Input_Parameters.NKS_IP.Maximum_Number_of_NKS_Iterations > 0) {
    error_flag = SolveNKS(batch_flag);
    if (error_flag) {
      cout << "\n Rte2D ERROR: Unable to perform implicit NKS computations.\n"
	   << flush;
      return error_flag;
    } // endif
  } // endif


  //--------------------------------------------------
  // Postprocess
  //--------------------------------------------------
  error_flag = PostProcess(batch_flag,command_flag);
  if (error_flag) {
    cout << "\n Rte2D ERROR: Unable to perform postprocessing.\n"
	 << flush;
    return error_flag;
  } // endif

  //--------------------------------------------------
  // Redirect
  //--------------------------------------------------
  if (command_flag == EXECUTE_CODE) 
    goto execute_new_calculation;
  else if (command_flag == CONTINUE_CODE) 
    goto continue_existing_calculation;

  // return error flag
  return error_flag;

}

/********************************************************  
 * Rte2DSolver :: ReadInputFile                         *
 *                                                      *
 * Set default values for the input solution parameters *
 * and then read user specified input values from the   *
 * specified input parameter file.                      *
 ********************************************************/
inline int Rte2DSolver::ReadInputFile(char *Input_File_Name_ptr, 
				      const int batch_flag,
				      int &command_flag) {
  
  // declares
  int error_flag(0);

  // The primary MPI processor processes the input parameter file.
  if (CFFC_Primary_MPI_Processor()) {
    if (!batch_flag) {
      cout << "\n Reading Rte2D input data file `"
	   << Input_File_Name_ptr << "'." << endl;
    } // endif
    error_flag = Process_Input_Control_Parameter_File(Input_Parameters,
						      Input_File_Name_ptr,
						      command_flag);
    if (!batch_flag && error_flag == 0) {
      cout << Input_Parameters << "\n";
      cout.flush();
    } // endif
  } // endif 
  
  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Broadcast input solution parameters to other MPI processors.
  CFFC_Broadcast_MPI(&error_flag, 1);
  if (error_flag != 0) return (error_flag);
  CFFC_Broadcast_MPI(&command_flag, 1);
  Broadcast_Input_Parameters(Input_Parameters);

  // common use variables
  NUM_VAR_RTE2D = Input_Parameters.Uo.NUM_VAR_RTE2D;

  // return error flag
  return error_flag;

}


/********************************************************  
 * Rte2DSolver :: SetupMesh                             *
 *                                                      *
 * Create initial mesh and allocate Rte2D solution      *
 * avariables for specified IBVP/BVP problem.           *
 * Additionally, create (allocate) multi-block quadtree *
 * data structure, create (allocate) array of local 2D  *
 * solution blocks, assign and create (allocate) 2D     *
 * equation solution blocks corresponding to the initial*
 * mesh.                                                *
 ********************************************************/
inline int Rte2DSolver::SetupMesh(const int batch_flag) {

  // declares
  int error_flag(0);
  
  //--------------------------------------------------
  // INITIAL GRID & SOLUTION SPACE
  // Create initial mesh and allocate Chem2D solution variables for 
  // specified IBVP/BVP problem. 
  //--------------------------------------------------
    
  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI(); 
  
  // Create initial mesh.  Read mesh from grid definition or data files 
  // when specified by input parameters.
  
  // The primary MPI processor creates the initial mesh.
  MeshBlk = NULL;
  if (CFFC_Primary_MPI_Processor()) {
    if (!batch_flag) cout << "\n Creating (or reading) initial quadrilateral multi-block mesh.";
    MeshBlk = Multi_Block_Grid(MeshBlk, 
			       Input_Parameters);
    
    if (MeshBlk == NULL) {
      error_flag = 1;
    } else if (Check_Multi_Block_Grid(MeshBlk,
				      Input_Parameters.Number_of_Blocks_Idir,
				      Input_Parameters.Number_of_Blocks_Jdir)) {
      error_flag = 1;
    } // endif 
    
    if (error_flag) {
      cout << "\n Rte2D ERROR: Unable to create valid Rte2D multi-block mesh.\n";
      cout.flush();
    }  // endif
  } // endif
  
  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // Broadcast the mesh to other MPI processors.
  CFFC_Broadcast_MPI(&error_flag, 1); // Broadcast mesh error flag.
  if (error_flag) return (error_flag);
  MeshBlk = Broadcast_Multi_Block_Grid(MeshBlk, 
                                       Input_Parameters);

  //--------------------------------------------------
  // MULTIBLOCK QUADTREE
  // Create (allocate) multi-block quadtree data structure, create
  // (allocate) array of local  2D equation solution blocks, 
  // assign and create (allocate) 2D equation solution blocks
  // corresponding to the initial mesh. 
  //--------------------------------------------------

  if (!batch_flag) cout << "\n Creating multi-block quadtree data structure and assigning"
                        << "\n  Rte2D solution blocks corresponding to initial mesh.";
  Local_SolnBlk = Create_Initial_Solution_Blocks(MeshBlk,
                                                 Local_SolnBlk,
                                                 Input_Parameters,
                                                 QuadTree,
                                                 List_of_Global_Solution_Blocks,
                                                 List_of_Local_Solution_Blocks);
  if (Local_SolnBlk == NULL) error_flag = 1;

  // return error flag
  return error_flag;

}

/********************************************************  
 * Rte2DSolver :: SetICs                                *
 *                                                      *
 * Prescribe the inial Rte2D solution data.  We either  *
 * read from a RESTART file or set the ICs manually.    *
 ********************************************************/
inline int Rte2DSolver::SetICs(const int batch_flag) {

  // declares
  int error_flag(0);

 
  // Initialize the conserved and primitive state
  // solution variables.
  
  if (!batch_flag) cout << "\n Prescribing Rte2D initial data.";

  //--------------------------------------------------
  // RESTART SOLUTION
  //--------------------------------------------------
  if (Input_Parameters.i_ICs == IC_RESTART) {

    // read quad tree data structure
    if (!batch_flag) cout << "\n Reading Rte2D solution from restart data files.";
    error_flag = Read_QuadTree(QuadTree,
			       List_of_Global_Solution_Blocks,
			       List_of_Local_Solution_Blocks, 
			       Input_Parameters);
    if (error_flag) {
      cout << "\n Rte2D ERROR: Unable to open Rte2D quadtree data file "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } /* endif */
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);

    // allocate message buffers
    Allocate_Message_Buffers(List_of_Local_Solution_Blocks,
			     Local_SolnBlk[0].NumVar()+NUM_COMP_VECTOR2D);

    // read restart solution file
    error_flag = Read_Restart_Solution(Local_SolnBlk, 
				       List_of_Local_Solution_Blocks, 
				       Input_Parameters,
				       number_of_time_steps,
				       Time,
				       processor_cpu_time);
    if (error_flag) {
      cout << "\n Rte2D ERROR: Unable to open Rte2D restart input data file(s) "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } // endif
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);

    // Ensure each processor has the correct time and time!!!
    number_of_time_steps = CFFC_Maximum_MPI(number_of_time_steps);
    Time = CFFC_Maximum_MPI(Time);
    processor_cpu_time.cput = CFFC_Maximum_MPI(processor_cpu_time.cput);
    
  //--------------------------------------------------
  // RESTART SOLUTION
  //--------------------------------------------------
  } else {
     ICs(Local_SolnBlk, 
         List_of_Local_Solution_Blocks, 
         Input_Parameters);

  } /* endif */

  //--------------------------------------------------
  // ADDITIONAL SETUP
  //--------------------------------------------------
  // Compute the scaling parameter required to treat the 2D grids as 3D. 
  ScaleGridTo3D(Local_SolnBlk, List_of_Local_Solution_Blocks);

  // return error flag
  return error_flag;

}


/********************************************************  
 * Rte2DSolver :: BroadcastSolution                     *
 *                                                      *
 * Send solution information between neighboring blocks *
 * to complete prescription of initial data.            *
 ********************************************************/
inline int Rte2DSolver::BroadcastSolution() {

  // declares
  int error_flag(0);

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();

  // send solution info
  error_flag = Send_All_Messages(Local_SolnBlk, 
                                 List_of_Local_Solution_Blocks,
                                 NUM_COMP_VECTOR2D,
                                 ON);

  if (!error_flag) error_flag = Send_All_Messages(Local_SolnBlk, 
                                                  List_of_Local_Solution_Blocks,
                                                  NUM_VAR_RTE2D,
                                                  OFF);

  if (error_flag) {
     cout << "\n Rte2D ERROR: Message passing error during Rte2D solution intialization "
          << "on processor "
          << List_of_Local_Solution_Blocks.ThisCPU
          << ".\n";
     cout.flush();
  } // endif
  error_flag = CFFC_OR_MPI(error_flag);

  // return error flag
  return error_flag;

}


/********************************************************  
 * Rte2DSolver :: PerformAMR                            *
 *                                                      *
 * Perform uniform, boundary, and initial mesh          *
 * refinement.                                          *
 ********************************************************/
inline int Rte2DSolver::PerformAMR(const int batch_flag) {

  // declares
  int error_flag(0);

  //--------------------------------------------------
  // Perform uniform, boundary, and, initial mesh refinement. 
  //--------------------------------------------------
  if (Input_Parameters.i_ICs != IC_RESTART) {
    
    //
    // uniform AMR
    //
    if (!batch_flag) cout << "\n Performing Rte2D uniform mesh refinement.";
    error_flag = Uniform_AMR(Local_SolnBlk,
			     Input_Parameters,
			     QuadTree,
			     List_of_Global_Solution_Blocks,
			     List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n Rte2D ERROR: Uniform AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
      cout.flush();
    } // endif
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;
    
    //
    // boundary AMR
    //
    if (!batch_flag) cout << "\n Performing Rte2D boundary mesh refinement.";
    error_flag = Boundary_AMR(Local_SolnBlk,
			      Input_Parameters,
			      QuadTree,
			      List_of_Global_Solution_Blocks,
			      List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n Rte2D ERROR: Boundary AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU << "." << endl;
      cout.flush();
    } // endif
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return error_flag;
    
    //
    // initial mesh refinement
    //
    if (!batch_flag) cout << "\n Performing Rte2D initial mesh refinement.";
    error_flag = Initial_AMR(Local_SolnBlk,
			     Input_Parameters,
			     QuadTree,
			     List_of_Global_Solution_Blocks,
			     List_of_Local_Solution_Blocks);
    if (error_flag) {
      cout << "\n Rte2D ERROR: Initial AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } // endif
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    
  } // endif
  
  //--------------------------------------------------
  // Output multi-block solution-adaptive quadrilateral mesh
  // statistics.
  //--------------------------------------------------

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


  // return error flag
  return error_flag;

}

/********************************************************  
 * Rte2DSolver :: PerformPeriodicAMR                    *
 *                                                      *
 * Perform uniform, boundary, and initial mesh          *
 * refinement.                                          *
 ********************************************************/
inline int Rte2DSolver::PerformPeriodicAMR(const int batch_flag) {

  // declares
  int error_flag(0);

  if (Input_Parameters.AMR && !first_step && 
      number_of_time_steps % Input_Parameters.AMR_Frequency == 0 ) {
    
    // evaluate limiters & peform AMR
    if (!batch_flag) 
      cout << "\n\n Refining Grid.  Performing adaptive mesh refinement at n = "
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
      cout << "\n Rte2D ERROR: Rte2D AMR error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } // endif
    error_flag = CFFC_OR_MPI(error_flag);
    
    // output to tecplot if there is an error
    if (error_flag) {
      Output_Tecplot(Local_SolnBlk,
		     List_of_Local_Solution_Blocks,
		     Input_Parameters,
		     number_of_time_steps,
		     Time);
      return (error_flag);
    } // endif
    
    // display info
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
  
  // return error flag
  return error_flag;

}

/********************************************************  
 * Rte2DSolver :: PerformMortonOrdering                 *
 *                                                      *
 * Morton ordering of initial solution blocks.  This    *
 * should be meshed with AMR, ie. when Refine_Grid      *
 * is done, call the ordering process.                  *
 ********************************************************/
inline int Rte2DSolver::PerformMortonOrdering(const int batch_flag) {

  // declares
  int error_flag(0);

  if (Input_Parameters.Morton) {
    
    // NOTES: Issue here with Input_Parameters.Maximum_Number_of_Time_Steps 
    // related to reading restart files

    // Perorm ordering
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
    } // endif
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    
    //Output space filling curve in Tecplot format
    if (!batch_flag) {
      cout << "\n Outputting space filling curve showing block loading for CPUs.\n";
      cout.flush();
    } // endif
    Morton_SFC_Output_Tecplot(Local_SolnBlk, 
			      Input_Parameters, 
			      List_of_Local_Solution_Blocks);
  } // endif

  // return error flag
  return error_flag;

}

/********************************************************  
 * Rte2DSolver :: PerformMortonReOrdering               *
 *                                                      *
 * Periodically re-order the solution blocks on the     *
 * processors.                                          *
 ********************************************************/
inline int Rte2DSolver::PerformMortonReOrdering(const int batch_flag) {

  // declares
  int error_flag(0);

  if (Input_Parameters.Morton && !first_step &&
      number_of_time_steps%Input_Parameters.Morton_Reordering_Frequency == 0) {
    
    // peroform re-ordering
    if (!batch_flag) 
      cout << "\n\n Applying Morton re-ordering algorithm to solution blocks at n = "
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
      cout <<"\n Rte2D ERROR: Morton re-ordering error on processor "
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
  } // endif
 
  // return error flag
  return error_flag;
 
}

/********************************************************  
 * Rte2DSolver :: OpenResidualFile                      *
 *                                                      *
 * Open the residual file for writing.                  *
 ********************************************************/
inline int Rte2DSolver::OpenResidualFile(ofstream &residual_file) {

  // declares
  int error_flag(0);

  if (CFFC_Primary_MPI_Processor()) {
    error_flag = Open_Progress_File(residual_file,
				    Input_Parameters.Output_File_Name,
				    number_of_time_steps);
    if (error_flag) {
      cout << "\n Rte2D ERROR: Unable to open residual file for Rte2D calculation.\n";
      cout.flush();
    } // endif
  } // endif
  
  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();
  CFFC_Broadcast_MPI(&error_flag, 1);
 
  // return error flag
  return (error_flag);
}

/********************************************************  
 * Rte2DSolver :: OutputProgress()                      *
 *                                                      *
 * Output progress to screen and file.                  *
 ********************************************************/
inline void Rte2DSolver::OutputProgress(const int batch_flag, 
					const double residual_l1_norm, 
					const double residual_l2_norm, 
					const double residual_max_norm,
					ofstream &residual_file) {

    //screen
    if (!batch_flag) Output_Progress_L2norm(number_of_time_steps,
					    Time*THOUSAND,
					    total_cpu_time,
					    residual_l2_norm,
					    first_step,
					    50);
//     if (!batch_flag) Output_Progress(number_of_time_steps,
// 				     Time*THOUSAND,
// 				     total_cpu_time,
// 				     residual_l1_norm,
// 				     first_step,
// 				     50);

    //residual to file
    if (CFFC_Primary_MPI_Processor() && !first_step)
      Output_Progress_to_File(residual_file,
			      number_of_time_steps,
			      Time*THOUSAND,
			      total_cpu_time,
			      residual_l1_norm,
			      residual_l2_norm,
			      residual_max_norm);
}

/********************************************************  
 * Rte2DSolver :: ComputeNorms                          *
 *                                                      *
 * Determine the L1, L2, and max norms of the solution  *
 * residual and update the CPU time for all processesors*
 ********************************************************/
inline void Rte2DSolver::ComputeNorms(double &residual_l1_norm, 
				      double &residual_l2_norm, 
				      double &residual_max_norm) {

  // L1 norm for all processors.
  residual_l1_norm = L1_Norm_Residual(Local_SolnBlk, 
				      List_of_Local_Solution_Blocks);
  residual_l1_norm = CFFC_Summation_MPI(residual_l1_norm);
  
  // L2 norm for all processors.
  residual_l2_norm = L2_Norm_Residual(Local_SolnBlk, 
				      List_of_Local_Solution_Blocks);
  residual_l2_norm = sqr(residual_l2_norm);
  residual_l2_norm = CFFC_Summation_MPI(residual_l2_norm);
  residual_l2_norm = sqrt(residual_l2_norm);
  
  // Max norm for all processors.
  residual_max_norm = Max_Norm_Residual(Local_SolnBlk, 
					List_of_Local_Solution_Blocks);
  residual_max_norm = CFFC_Maximum_MPI(residual_max_norm);
   
}

/********************************************************  
 * Rte2DSolver :: WritePeriodicRestart                  *
 *                                                      *
 * Periodically save restart solution files.            *
 ********************************************************/
inline int Rte2DSolver::WritePeriodicRestart(const int batch_flag) {

  // declares
  int error_flag(0);

  if (!first_step && 
      number_of_time_steps%Input_Parameters.Restart_Solution_Save_Frequency==0) {
    
    // write quadtree
    if (!batch_flag) 
      cout << "\n\n  Saving Rte2D solution to restart data file(s) after"
	   << " n = " << number_of_time_steps << " steps (iterations).";
    error_flag = Write_QuadTree(QuadTree,
				Input_Parameters);
    if (error_flag) {
      cout << "\n Rte2D ERROR: Unable to open Rte2D quadtree data file "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } // endif
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    
    // write restart solution
    error_flag = Write_Restart_Solution(Local_SolnBlk, 
					List_of_Local_Solution_Blocks, 
					Input_Parameters,
					number_of_time_steps,
					Time,
					processor_cpu_time);
    if (error_flag) {
      cout << "\n Rte2D ERROR: Unable to open Rte2D restart output data file(s) "
	   << "on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } // endif
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
    cout << "\n";
    cout.flush();

  } // endif

  // return error flag
  return error_flag;
}

/****************************************************************************
 *********************** MULTIGRID SOLVER ***********************************
   Rte2DSolver :: SolveMultigrid

   Solve IBVP or BVP for conservation form of 2D planar/axisymmetric 
   radiation transport through absorbing, anisotropically scattering, 
   non-gray media on multi-block solution-adaptive quadrilateral mesh.                                  
 ****************************************************************************
 ****************************************************************************/  
inline int Rte2DSolver::SolveMultigrid(const int batch_flag) {

  // declares
  int error_flag(0);
  ofstream residual_file;

  // Allocate memory for multigrid solver.
  error_flag = MGSolver.allocate(Local_SolnBlk,
				 &QuadTree,
				 &List_of_Global_Solution_Blocks,
				 &List_of_Local_Solution_Blocks,
				 &Input_Parameters);
  if (error_flag) {
    cout << "\n Rte2D ERROR: Unable to allocate memory for multigrid solver.\n";
    cout.flush();
  }
  CFFC_Broadcast_MPI(&error_flag,1);
  if (error_flag) return error_flag;

  // Execute multigrid solver.
  start_explicit = clock();  
  error_flag = MGSolver.Execute(batch_flag,
				number_of_time_steps,
				Time,
				processor_cpu_time,
				total_cpu_time,
				residual_file);
  end_explicit = clock();  
  if (error_flag) {
    cout << "\n Euler2D ERROR: Error during multigrid solution.\n";
    cout.flush();
  }
  CFFC_Broadcast_MPI(&error_flag,1);

  // return error flag
  return error_flag;

}

/****************************************************************************
 ********************** SPACE-MARCHING SOLVER *******************************
   Rte2DSolver :: SolveSpaceMarch

   Solve IBVP or BVP for conservation form of 2D planar/axisymmetric 
   radiation transport through absorbing, anisotropically scattering, 
   non-gray media on multi-block solution-adaptive quadrilateral mesh.                                  
 ****************************************************************************
 ****************************************************************************/  
inline int Rte2DSolver::SolveSpaceMarch(const int batch_flag) {

  // declares
  int error_flag(0);
  int i_SpaceMarch_Scheme;
  double residual_l1_norm, residual_l2_norm, residual_max_norm;
  ofstream residual_file;

  //--------------------------------------------------
  // Open residual file and reset the CPU time
  //--------------------------------------------------
  error_flag = OpenResidualFile(residual_file);
  if (error_flag) return (error_flag);
  first_step = true;
  processor_cpu_time.reset();
  
  //--------------------------------------------------
  // Start explicit clock
  //--------------------------------------------------
  if (!batch_flag){ 
    cout << "\n\n Beginning Rte2D computations on "
	 << Date_And_Time() << ".\n\n"; 
    start_explicit = clock();  
  } // endif
  
  
  //==================================================
  //
  // Perform required number of iterations (space marchs).
  //
  //==================================================
  while (Input_Parameters.Maximum_Number_of_Time_Steps > 0 &&
	 number_of_time_steps <= Input_Parameters.Maximum_Number_of_Time_Steps) {
    
    //--------------------------------------------------
    // MORTON ORDERING: Periodically re-order the solution 
    // blocks on the processors.
    //--------------------------------------------------
    PerformMortonReOrdering(batch_flag);
    if (error_flag) return (error_flag);

    //--------------------------------------------------
    // MESH REFINEMENT: Periodically refine the mesh (AMR). 
    //--------------------------------------------------
    PerformPeriodicAMR(batch_flag);
    if (error_flag) return (error_flag);

    //--------------------------------------------------
    // NORMS: Determine the L1, L2, and max norms of the
    // solution residual. 
    //--------------------------------------------------
    ComputeNorms(residual_l1_norm, residual_l2_norm, residual_max_norm);

    //--------------------------------------------------
    // UPDATE CPU
    //--------------------------------------------------
    // Update CPU time used for the calculation so far.
    processor_cpu_time.update();
  
    // Total CPU time for all processors.
    total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput);

    //--------------------------------------------------
    // RESTART: Periodically save restart solution files.
    //--------------------------------------------------
    WritePeriodicRestart(batch_flag);
    if (error_flag) return (error_flag);

    //--------------------------------------------------
    // PROGRESS: Output progress information for the calculation.
    //--------------------------------------------------
    OutputProgress(batch_flag, residual_l1_norm, residual_l2_norm, 
		   residual_max_norm, residual_file);
      
    //--------------------------------------------------
    // CALCULATION CHECK: Check to see if calculations are 
    // complete and if so jump of out of this infinite loop.   
    //--------------------------------------------------
    if (number_of_time_steps >= 
	Input_Parameters.Maximum_Number_of_Time_Steps) break;	

    //--------------------------------------------------
    // BLOCK SOLUTION UPDATE: Update solution for next 
    // time step using a space marching scheme.
    //--------------------------------------------------
    // 1a. Exchange solution information between neighbouring blocks.
    error_flag = Send_All_Messages(Local_SolnBlk, 
				   List_of_Local_Solution_Blocks,
				   NUM_VAR_RTE2D, 
				   OFF);
    if (error_flag) {
      cout << "\n Rte2D ERROR: Rte2D message passing error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } // endif
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
      
    // 1b. Send boundary ref states at block interfaces.
    error_flag = Send_Boundary_Ref_States(Local_SolnBlk, 
					  List_of_Local_Solution_Blocks,
					  NUM_VAR_RTE2D);
    if (error_flag) {
      cout << "\n Rte2D ERROR: Rte2D boundary ref state message passing error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } // endif
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
      
    //--------------------------------------------------
    // BOUNDARY CONDITIONS
    //--------------------------------------------------
    // 2. Apply boundary conditions for stage.
    BCs_Space_March(Local_SolnBlk, 
		    List_of_Local_Solution_Blocks,
		    Input_Parameters);			
	
    //--------------------------------------------------
    // UPDATE SOLUTION
    //--------------------------------------------------
    // use upwind for first step of CLAM
    if ( first_step && ( Input_Parameters.i_SpaceMarch_Scheme==SPACE_MARCH_CLAM || 
			 Input_Parameters.i_SpaceMarch_Scheme==SPACE_MARCH_GM ) ) {
      i_SpaceMarch_Scheme = Input_Parameters.i_SpaceMarch_Scheme;
      Input_Parameters.i_SpaceMarch_Scheme = SPACE_MARCH_UPWIND;
    } else {
      i_SpaceMarch_Scheme = Input_Parameters.i_SpaceMarch_Scheme;	  
    }

    // 3. Determine solution residuals for stage.
    error_flag = dUdt_Space_March(Local_SolnBlk, 
				  List_of_Global_Solution_Blocks,
				  List_of_Local_Solution_Blocks,
				  Input_Parameters);
    if (error_flag) {
      cout << "\n Rte2D ERROR: Rte2D solution residual error on processor "
	   << List_of_Local_Solution_Blocks.ThisCPU
	   << ".\n";
      cout.flush();
    } // endif
    error_flag = CFFC_OR_MPI(error_flag);
    if (error_flag) return (error_flag);
      
    // reset differencing scheme
    if (first_step)
      Input_Parameters.i_SpaceMarch_Scheme = i_SpaceMarch_Scheme;
      
    //--------------------------------------------------
    // UPDATE TIMER & COUNTER: Update time and time step 
    // counter.
    //--------------------------------------------------
      
    if (first_step) first_step = false;
    number_of_time_steps = number_of_time_steps + 1;
      

  //==================================================
  } // endwhile
  //==================================================
   
  if (!batch_flag) { 
    cout << "\n\n Rte2D computations complete on " 
	 << Date_And_Time() << ".\n"; 
    end_explicit = clock(); 
  }// endif

  //--------------------------------------------------
  // UPDATE GHOSTCELLS: Update ghostcell information and 
  // prescribe boundary conditions to ensure that the 
  // solution is consistent on each block. 
  //--------------------------------------------------
  
  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI(); 

  // broadcast solution
  error_flag = Send_All_Messages(Local_SolnBlk, 
				 List_of_Local_Solution_Blocks,
				 NUM_VAR_RTE2D,
				 OFF);
  if (error_flag) {
    cout << "\n Rte2D ERROR: Rte2D message passing error on processor "
	 << List_of_Local_Solution_Blocks.ThisCPU
	 << ".\n";
    cout.flush();
  } // endif
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return (error_flag);
    
  // broadcast boundary ref states
  error_flag = Send_Boundary_Ref_States(Local_SolnBlk, 
					List_of_Local_Solution_Blocks,
					NUM_VAR_RTE2D);
  if (error_flag) {
    cout << "\n Rte2D ERROR: Rte2D boundary ref state message passing error on processor "
	 << List_of_Local_Solution_Blocks.ThisCPU
	 << ".\n";
    cout.flush();
  } // endif
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return (error_flag);

  // use time-marching BCs to fudge the ghost cells properly
  BCs(Local_SolnBlk, List_of_Local_Solution_Blocks, Input_Parameters);
    
  // Close residual file.
  if (CFFC_Primary_MPI_Processor()) 
    error_flag = Close_Progress_File(residual_file);

  // return error flag
  return error_flag;
  
}

/****************************************************************************
 ********************** TIME-MARCHING SOLVER ********************************
   Rte2DSolver :: SolveTimeMarch

   Solve IBVP or BVP for conservation form of 2D planar/axisymmetric 
   radiation transport through absorbing, anisotropically scattering, 
   non-gray media on multi-block solution-adaptive quadrilateral mesh.                                  
 ****************************************************************************
 ****************************************************************************/  
inline int Rte2DSolver::SolveExplicitTimeStep(const int batch_flag) {
  
  // declares
  int error_flag(0);
  double dTime;
  double residual_l1_norm, residual_l2_norm, residual_max_norm;
  ofstream residual_file;

  //--------------------------------------------------
  // Open residual file and reset the CPU time
  //--------------------------------------------------
  error_flag = OpenResidualFile(residual_file);
  if (error_flag) return (error_flag);
  first_step = true;
  processor_cpu_time.reset();
     
  //--------------------------------------------------
  // Start explicit clock
  //--------------------------------------------------
  if (!batch_flag) { 
    cout << "\n\n Beginning Rte2D computations on "
	 << Date_And_Time() << ".\n\n"; 
    start_explicit = clock();  
  } // endif


  //==================================================
  //
  // Perform required number of iterations (time-steps).
  //
  //==================================================
  while ( (!Input_Parameters.Time_Accurate &&
	   Input_Parameters.Maximum_Number_of_Time_Steps > 0 &&
	   number_of_time_steps <= Input_Parameters.Maximum_Number_of_Time_Steps) ||
	  (Input_Parameters.Time_Accurate &&
	   Input_Parameters.Time_Max >= Time) ) {

    //--------------------------------------------------
    // MORTON ORDERING: Periodically re-order the solution 
    // blocks on the processors.
    //--------------------------------------------------
    PerformMortonReOrdering(batch_flag);
    if (error_flag) return (error_flag);

    //--------------------------------------------------
    // MESH REFINEMENT: Periodically refine the mesh (AMR). 
    //--------------------------------------------------
    PerformPeriodicAMR(batch_flag);
    if (error_flag) return (error_flag);

    //--------------------------------------------------
    // TIME STEPS: Determine local and global time steps.
    //--------------------------------------------------
    // determine max time step
    dTime = CFL(Local_SolnBlk, 
		List_of_Local_Solution_Blocks,
		Input_Parameters);
    
    // Find global minimum time step for all processors.
    dTime = CFFC_Minimum_MPI(dTime); 
    
    //
    // Compute time step for time accurate case
    //
    if (Input_Parameters.Time_Accurate) {
      
      // non-multistage optimally smoothing time integration
      if ((Input_Parameters.i_Time_Integration != TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) &&
	  (Time + Input_Parameters.CFL_Number*dTime > Input_Parameters.Time_Max)) {
	dTime = (Input_Parameters.Time_Max-Time)/Input_Parameters.CFL_Number;

      // multistage optimally smoothing time integration
      } else if (Time + 
		 Input_Parameters.CFL_Number*dTime*
		 MultiStage_Optimally_Smoothing(Input_Parameters.N_Stage, 
						Input_Parameters.N_Stage,
						Input_Parameters.i_Limiter) > Input_Parameters.Time_Max) {
	dTime = (Input_Parameters.Time_Max-Time) /
	        (Input_Parameters.CFL_Number*
		 MultiStage_Optimally_Smoothing(Input_Parameters.N_Stage, 
						Input_Parameters.N_Stage,
						Input_Parameters.i_Limiter));
      } // endif
    } // endif


    //
    // if not local time stepping, set global time step.
    //
    if (!Input_Parameters.Local_Time_Stepping)
      Set_Global_TimeStep(Local_SolnBlk, 
			  List_of_Local_Solution_Blocks,
			  dTime);
	
	
    //--------------------------------------------------
    // NORMS: Determine the L1, L2, and max norms of the
    // solution residual. 
    //--------------------------------------------------
    ComputeNorms(residual_l1_norm, residual_l2_norm, residual_max_norm);

    //--------------------------------------------------
    // UPDATE CPU
    //--------------------------------------------------
    // Update CPU time used for the calculation so far.
    processor_cpu_time.update();
  
    // Total CPU time for all processors.
    total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput);

    //--------------------------------------------------
    // RESTART: Periodically save restart solution files.
    //--------------------------------------------------
    WritePeriodicRestart(batch_flag);
    if (error_flag) return (error_flag);
	
    //--------------------------------------------------
    // PROGRESS: Output progress information for the calculation.
    //--------------------------------------------------
    OutputProgress(batch_flag, residual_l1_norm, residual_l2_norm, 
		   residual_max_norm, residual_file);
	
    //--------------------------------------------------
    // CALCULATION CHECK: Check to see if calculations are 
    // complete and if so jump of out of this infinite loop.   
    //--------------------------------------------------
    if (!Input_Parameters.Time_Accurate && 
	number_of_time_steps >= Input_Parameters.Maximum_Number_of_Time_Steps) 
      break;
    if (Input_Parameters.Time_Accurate &&
	Time >= Input_Parameters.Time_Max) 
      break;
	
    //--------------------------------------------------
    // LIMITER FREEZE: Freeze limiters as necessary
    //--------------------------------------------------
    if (!first_step && Input_Parameters.Freeze_Limiter &&
	residual_l2_norm <= Input_Parameters.Freeze_Limiter_Residual_Level) {
      Freeze_Limiters(Local_SolnBlk, 
		      List_of_Local_Solution_Blocks);
    } // endif

    //--------------------------------------------------
    // BLOCK SOLUTION UPDATE: Update solution for next 
    // time step using a space marching scheme.
    //--------------------------------------------------
    for ( int i_stage = 1 ; i_stage <= Input_Parameters.N_Stage ; ++i_stage ) {

      // 1. Exchange solution information between neighbouring blocks.
      error_flag = Send_All_Messages(Local_SolnBlk, 
				     List_of_Local_Solution_Blocks,
				     NUM_VAR_RTE2D, 
				     OFF);
      if (error_flag) {
	cout << "\n Rte2D ERROR: Rte2D message passing error on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
      } // endif
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
	  
      //--------------------------------------------------
      // BOUNDARY CONDITIONS
      //--------------------------------------------------
      // 2. Apply boundary conditions for stage.
      BCs(Local_SolnBlk, 
	  List_of_Local_Solution_Blocks,
	  Input_Parameters);
	  
      //--------------------------------------------------
      // UPDATE SOLUTION
      //--------------------------------------------------
      // 3. Determine solution residuals for stage.
      error_flag = dUdt_Multistage_Explicit(Local_SolnBlk, 
					    List_of_Global_Solution_Blocks,
					    List_of_Local_Solution_Blocks,
					    Input_Parameters,
					    i_stage);
      if (error_flag) {
	cout << "\n Rte2D ERROR: Rte2D solution residual error on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
      } // endif
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
      
      // 4. Send boundary flux corrections at block interfaces with resolution changes.
      error_flag = Send_Conservative_Flux_Corrections(Local_SolnBlk, 
						      List_of_Local_Solution_Blocks,
						      NUM_VAR_RTE2D);
      if (error_flag) {
	cout << "\n Rte2D ERROR: Rte2D flux correction message passing error on processor "
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
      
      // 6. Smooth the solution residual using implicit residual smoothing.
      if (Input_Parameters.Residual_Smoothing)
	Residual_Smoothing(Local_SolnBlk,
			   List_of_Local_Solution_Blocks,
			   Input_Parameters,
			   i_stage);

      // 7. Update solution for stage.
      error_flag = Update_Solution_Multistage_Explicit(Local_SolnBlk, 
						       List_of_Local_Solution_Blocks,
						       Input_Parameters,
						       i_stage);
      if (error_flag) {
	cout << "\n Rte2D ERROR: Rte2D solution update error on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
      } // endif
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
      
    } // endfor
	
    //--------------------------------------------------
    // UPDATE TIMER & COUNTER: Update time and time step 
    // counter.
    //--------------------------------------------------
    if (first_step) first_step = false;
    number_of_time_steps = number_of_time_steps + 1;
    
    //
    // Time accurate case
    //
    // non-multistage optimally smoothing
    if (Input_Parameters.i_Time_Integration != 
	TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING) {
      Time = Time + Input_Parameters.CFL_Number*dTime;

    // multistage optimally smoothing
    } else {
      Time = Time + 
	     Input_Parameters.CFL_Number*dTime*
	     MultiStage_Optimally_Smoothing(Input_Parameters.N_Stage, 
					    Input_Parameters.N_Stage,
					    Input_Parameters.i_Limiter);
    } // endif
    
  //==================================================
  } // endwhile
  //==================================================

  if (!batch_flag){ 
    cout << "\n\n Rte2D computations complete on " 
	 << Date_And_Time() << ".\n"; 
    end_explicit = clock(); 
  } // endif

  //--------------------------------------------------
  // UPDATE GHOSTCELLS: Update ghostcell information and 
  // prescribe boundary conditions to ensure that the 
  // solution is consistent on each block. 
  //--------------------------------------------------
  
  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();
    
  // broadcast solution
  error_flag = Send_All_Messages(Local_SolnBlk, 
				 List_of_Local_Solution_Blocks,
				 NUM_VAR_RTE2D,
				 OFF);
  if (error_flag) {
    cout << "\n Rte2D ERROR: Rte2D message passing error on processor "
	 << List_of_Local_Solution_Blocks.ThisCPU
	 << ".\n";
    cout.flush();
  } // endif
  error_flag = CFFC_OR_MPI(error_flag);
  if (error_flag) return (error_flag);
  
  // recompute BCs
  BCs(Local_SolnBlk, 
      List_of_Local_Solution_Blocks,
      Input_Parameters);
  
  // Close residual file.
  if (CFFC_Primary_MPI_Processor()) 
    error_flag = Close_Progress_File(residual_file);


  // return error flag
  return error_flag;

}

/****************************************************************************
 *************************** NKS SOLVER *************************************
   Rte2DSolver :: SolveTimeMarch

   Solve IBVP or BVP for conservation form of 2D planar/axisymmetric 
   radiation transport through absorbing, anisotropically scattering, 
   non-gray media on multi-block solution-adaptive quadrilateral mesh.                                  
 ****************************************************************************
 ****************************************************************************/  
inline int Rte2DSolver::SolveNKS(const int batch_flag) {

  // declares
  int error_flag(0);
  time_t start_NKS, end_NKS;
  ofstream residual_file;

  //--------------------------------------------------
  // Open residual file and reset the CPU time
  //--------------------------------------------------
  OpenResidualFile(residual_file);
  if (error_flag) return (error_flag);

  //Turn Limiter Freezing OFF for startup 
  Evaluate_Limiters(Local_SolnBlk, List_of_Local_Solution_Blocks);

  //--------------------------------------------------
  // Perform NKS computations.
  //--------------------------------------------------
  if (!batch_flag){ cout << "\n\n Beginning Rte2D NKS computations on " << Date_And_Time() << ".\n\n"; time(&start_NKS); }

  //Store Explicit times for output
  CPUTime Explicit_processor_cpu_time = processor_cpu_time;
  CPUTime Explicit_total_cpu_time =  total_cpu_time;
    
  //Perform NKS Iterations 
  error_flag = Newton_Krylov_Schwarz_Solver<Rte2D_State,
                                            Rte2D_Quad_Block,
                                            Rte2D_Input_Parameters>(processor_cpu_time,   //NKS_processor_cpu_time,
								    residual_file,
								    number_of_time_steps, // explicit time steps
								    Time,
								    Local_SolnBlk, 
								    List_of_Local_Solution_Blocks,
								    Input_Parameters);
  
  // update time
  processor_cpu_time.update();
  total_cpu_time.cput = CFFC_Summation_MPI(processor_cpu_time.cput);  

  if (error_flag) {
    if (CFFC_Primary_MPI_Processor()) { 
      cout << "\n Rte2D_NKS ERROR: Rte2D solution error on processor " 
	   << List_of_Local_Solution_Blocks.ThisCPU << ".\n";
      cout.flush();
    } // endif
  } // endif

  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();
  CFFC_Broadcast_MPI(&error_flag, 1);
  if (error_flag) return (error_flag);
  
  //--------------------------------------------------
  // Output NKS stats
  //--------------------------------------------------
  if (!batch_flag) { 
    cout << "\n\n Rte2D NKS computations complete on " 
	 << Date_And_Time() 
	 << ".\n"; 
    time(&end_NKS); }
  
  if (!batch_flag) {
    cout<<"\n ----------------------------------------------------------------";
    cout<<"\n -------- Solution Computations Summary in minutes --------------";
    cout<<"\n ----------------------------------------------------------------";
    cout<<"\n Total Startup CPU Time\t\t = "<<Explicit_total_cpu_time.min();
    cout<<"\n Total NKS CPU Time \t\t = "<<total_cpu_time.min()-Explicit_total_cpu_time.min();
    cout<<"\n Total CPU Time \t\t = "<<total_cpu_time.min(); 
    cout<<"\n Total Startup Clock Time\t = "<<difftime(end_explicit,start_explicit)/60.0;
    cout<<"\n Total NKS Clock Time\t\t = "<<difftime(end_NKS,start_NKS)/60.0;
    cout<<"\n Total Clock Time\t\t = "<<difftime(end_NKS,start_explicit)/60.0;
    cout<<"\n ----------------------------------------------------------------";
    cout<<"\n ----------------------------------------------------------------";
    cout<<"\n ----------------------------------------------------------------\n";
  } 
  //Also want to output total GMRES & NKS Iterations, and maybe max memory usage possibly??

  // Close residual file.
  if (CFFC_Primary_MPI_Processor()) 
    error_flag = Close_Progress_File(residual_file); 
     
  //add implicit and explicit times 
  //total_cpu_time.cput += NKS_total_cpu_time.cput;

  // return error flag
  return error_flag;

}

/***************************************************************************
 ************************** POST PROCESSSING *******************************
   Rte2DSolver :: PostProcess

   Solution calculations complete. Write 2D solution to output and restart files  
   as required, reset solution parameters, and run other cases as specified 
   by input parameters.        
 *************************************************************************** 
 ****************************************************************************/ 
inline int Rte2DSolver::PostProcess(const int batch_flag,int &command_flag) {

  // declares
  int error_flag(0);
  int line_number;
  
  // MPI barrier to ensure processor synchronization.
  CFFC_Barrier_MPI();
  
  //==================================================
  //
  // Continue processing input file
  //
  //==================================================
  while (1) {

    
    //--------------------------------------------------
    // get next input value
    //--------------------------------------------------
    if (CFFC_Primary_MPI_Processor()) {
      Get_Next_Input_Control_Parameter(Input_Parameters);
      command_flag = Parse_Next_Input_Control_Parameter(Input_Parameters);
      line_number = Input_Parameters.Line_Number;
    } // endif

    // MPI barrier to ensure processor synchronization.
    CFFC_Barrier_MPI();

    //broadcast input parameters
    Broadcast_Input_Parameters(Input_Parameters);
    CFFC_Broadcast_MPI(&command_flag, 1);
   
    //--------------------------------------------------
    // EXECUTE CODE
    //--------------------------------------------------
    if (command_flag == EXECUTE_CODE) {

      //
      // Deallocate memory for 2D Rte equation solution.
      //
      if (!batch_flag) cout << "\n Deallocating Rte2D solution variables.";
      DeallocateSoln();

      // Output input parameters for new caluculation.
      if (!batch_flag)  {
	cout << "\n\n Starting a new calculation.";
	cout << Input_Parameters << "\n";
	cout.flush();
      } /* endif */
      // Execute new calculation.

      //return error_flag
      return error_flag;
      
    //--------------------------------------------------
    // TERMINATE CODE
    //--------------------------------------------------
    } else if (command_flag == TERMINATE_CODE) {

      //
      // Deallocate memory for 2D Rte equation solution.
      //
      if (!batch_flag) cout << "\n Deallocating Rte2D solution variables.";
      DeallocateSoln();
      

      // Close input data file.
      if (!batch_flag) cout << "\n\n Closing Rte2D input data file.";
      if (CFFC_Primary_MPI_Processor()) Close_Input_File(Input_Parameters);

      // Terminate calculation.
      return (error_flag);
      

    //--------------------------------------------------
    // CONTINUE CODE ie RESTART
    //--------------------------------------------------
    } else if (command_flag == CONTINUE_CODE) {

      // Reset maximum time step counter.
      Input_Parameters.Maximum_Number_of_Time_Steps += number_of_time_steps;

      // Output input parameters for continuing calculation.
      if (!batch_flag)  {
	cout << "\n\n Continuing existing calculation.";
	cout << Input_Parameters << "\n";
	cout.flush();
      } // endif

      // Continue existing calculation.
      return (error_flag);
      

    //--------------------------------------------------
    // REFINE GRID (AMR)
    //--------------------------------------------------
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
		       ON,
		       ON);
      if (error_flag) {
	cout << "\n Rte2D ERROR: Rte2D AMR error on processor "
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

    //--------------------------------------------------
    // MORTON ORDERING
    //--------------------------------------------------
    } else if (command_flag == MORTON_ORDERING_CODE) {

      // perform reordering
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
	cout <<"\n Rte2D ERROR: Morton re-ordering error on processor "
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
      
    //--------------------------------------------------
    // WRITE OUTPUT @ NODES
    //--------------------------------------------------
    } else if (command_flag == WRITE_OUTPUT_CODE) {

      // Output solution data.
      if (!batch_flag) cout << "\n Writing Rte2D solution to output data file(s).";
      error_flag = Output_Tecplot(Local_SolnBlk, 
				  List_of_Local_Solution_Blocks, 
				  Input_Parameters,
				  number_of_time_steps,
				  Time);
      if (error_flag) {
	cout << "\n Rte2D ERROR: Unable to open Rte2D output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
      } // endif
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
      

    //--------------------------------------------------
    // WRITE OUTPUT CELL-CENTERED
    //--------------------------------------------------
    } else if (command_flag == WRITE_OUTPUT_CELLS_CODE) {

      // Output solution data.
      if (!batch_flag) cout << "\n Writing cell-centered Rte2D solution to output data file(s).";
      error_flag = Output_Cells_Tecplot(Local_SolnBlk, 
					List_of_Local_Solution_Blocks, 
					Input_Parameters,
					number_of_time_steps,
					Time);
      if (error_flag) {
	cout << "\n Rte2D ERROR: Unable to open Rte2D cell output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
      } // endif
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
      

    //--------------------------------------------------
    // WRITE OUTPUT NODES
    //--------------------------------------------------
    } else if (command_flag == WRITE_OUTPUT_NODES_CODE) {

      // Output node data.
      if (!batch_flag) cout << "\n Writing Euler2D node locations to output data file(s).";
      error_flag = Output_Nodes_Tecplot(Local_SolnBlk,
					List_of_Local_Solution_Blocks,
					Input_Parameters,
					number_of_time_steps,
					Time);
      if (error_flag) {
	cout << "\n Euler2D ERROR: Unable to open Euler2D nodes output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << "." << endl;
      } // endif
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;
      
    //--------------------------------------------------
    // WRITE OUTPUT GRADIENTS
    //--------------------------------------------------
    } else if (command_flag == WRITE_OUTPUT_GRADIENTS_CODE) {
      
      // output solution gradient info
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
      } // endif
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return error_flag;

    //--------------------------------------------------
    // WRITE RESTART FILE
    //--------------------------------------------------
    } else if (command_flag == WRITE_RESTART_CODE) {
      
      // Write quadtree.
      if (!batch_flag) cout << "\n Writing Rte2D solution to restart data file(s).";
      error_flag = Write_QuadTree(QuadTree,
				  Input_Parameters);
      if (error_flag) {
	cout << "\n Rte2D ERROR: Unable to open Rte2D quadtree data file "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
      } // endif
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);

      // write restart file
      error_flag = Write_Restart_Solution(Local_SolnBlk, 
					  List_of_Local_Solution_Blocks, 
					  Input_Parameters,
					  number_of_time_steps,
					  Time,
					  processor_cpu_time);
      if (error_flag) {
	cout << "\n Rte2D ERROR: Unable to open Rte2D restart output data file(s) "
	     << "on processor "
	     << List_of_Local_Solution_Blocks.ThisCPU
	     << ".\n";
	cout.flush();
      } // endif
      error_flag = CFFC_OR_MPI(error_flag);
      if (error_flag) return (error_flag);

      
    //--------------------------------------------------
    // WRITE OUTPUT GRID
    //--------------------------------------------------
    } else if (command_flag == WRITE_OUTPUT_GRID_CODE) {

      // Output multi-block solution-adaptive mesh data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing Rte2D multi-block mesh to grid data output file.";
	error_flag = Output_Tecplot(MeshBlk,
				    Input_Parameters);
	if (error_flag) {
	  cout << "\n Rte2D ERROR: Unable to open Rte2D mesh data output file.\n";
	  cout.flush();
	} // endif
      } // endif
      CFFC_Broadcast_MPI(&error_flag, 1);
      if (error_flag) return (error_flag);
      
    //--------------------------------------------------
    // WRITE GRID DEFINITION
    //--------------------------------------------------
    } else if (command_flag == WRITE_GRID_DEFINITION_CODE) {
      
      // Write multi-block solution-adaptive mesh definition files.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing Rte2D multi-block mesh to grid definition files.";
	error_flag = Write_Multi_Block_Grid_Definition(MeshBlk,
						       Input_Parameters);
	error_flag = Write_Multi_Block_Grid(MeshBlk,
					    Input_Parameters);
	if (error_flag) {
	  cout << "\n Rte2D ERROR: Unable to open Rte2D multi-block mesh definition files.\n";
	  cout.flush();
	} // endif
      } // endif
      CFFC_Broadcast_MPI(&error_flag, 1);
      if (error_flag) return (error_flag);
      

    //--------------------------------------------------
    // WRITE OUTPUT GRID NODES
    //--------------------------------------------------
    } else if (command_flag == WRITE_OUTPUT_GRID_NODES_CODE) {
      
      // Output multi-block solution-adaptive mesh node data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing Rte2D multi-block mesh to node data output file.";
	error_flag = Output_Nodes_Tecplot(MeshBlk,
					  Input_Parameters);
	if (error_flag) {
	  cout << "\n Rte2D ERROR: Unable to open Rte2D mesh node data output file.\n";
	  cout.flush();
	} // endif
      } // endif
      CFFC_Broadcast_MPI(&error_flag, 1);
      if (error_flag) return (error_flag);
      

    //--------------------------------------------------
    // WRITE OUTPUT GRID CELLS
    //--------------------------------------------------
    } else if (command_flag == WRITE_OUTPUT_GRID_CELLS_CODE) {
      
      // Output multi-block solution-adaptive mesh cell data file.
      if (CFFC_Primary_MPI_Processor()) {
	if (!batch_flag) cout << "\n Writing Rte2D multi-block mesh to cell data output file.";
	error_flag = Output_Cells_Tecplot(MeshBlk,
					  Input_Parameters);
	if (error_flag) {
	  cout << "\n Rte2D ERROR: Unable to open Rte2D mesh cell data output file.\n";
	  cout.flush();
	} // endif
      } // endif
      CFFC_Broadcast_MPI(&error_flag, 1);
      if (error_flag) return (error_flag);
      

    //--------------------------------------------------
    // WRITE BLACK_ENCLOSURE
    //--------------------------------------------------
    } else if (command_flag == WRITE_OUTPUT_BLACK_ENCLOSURE_CODE) {

      // write out exact solution for uniformly absorbing, 
      // isothermal media in a black, cold enclosure.
      if (!batch_flag) cout << endl << " Writing exact solution and error norms for Black Enclosure.";
      error_flag = Output_Exact(Local_SolnBlk,
				List_of_Local_Solution_Blocks,
				Input_Parameters);
      if (error_flag) {
	cout << endl << "\n Rte2D ERROR: Unable to open Rte2D Exact Solution output file." << endl;
	cout.flush();
      }
      CFFC_Broadcast_MPI(&error_flag,1);
      if (error_flag) return error_flag;
   

    //--------------------------------------------------
    // NOT A VALID INPUT_CODE
    //--------------------------------------------------
    } else if (command_flag == INVALID_INPUT_CODE ||
               command_flag == INVALID_INPUT_VALUE) {
      line_number = -line_number;
      cout << "\n Rte2D ERROR: Error reading Rte2D data at line #"
	   << -line_number  << " of input data file.\n";
      cout.flush();
      return (line_number);
    } // endif
    
  //==================================================
  } // endwhile
  //==================================================
  
  // return the error flag
  return error_flag;
  
}


#endif // _RTE2D_QUADSOLVER__INCLUDED
