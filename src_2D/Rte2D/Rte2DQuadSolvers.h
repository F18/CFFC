/**************************************************************** 
   Rte2DQuadSolvers.h

   2D Axisymmetric Solver for.... This header defines a class 
   'Rte2DSolver' which encapsulates the main Rte2D solver.  It
   can be called in full standalone mode (ie. by itself), or 
   in sequential solve mode (ie. say from Chem2D).
     
    NOTES:      

   - based on Euler2DQuadSolvers.cc
*****************************************************************/

#ifndef _RTE2D_QUADSOLVER__INCLUDED
#define _RTE2D_QUADSOLVER__INCLUDED


// Include 2D Rte quadrilateral mesh solution header file.
#include "Rte2DQuad.h"

// Include Rte2D Multigrid Specialization header file.
#include "Rte2DQuadMultigrid.h"

// Include Rte2D Newton-Krylov-Schwarz Specialization header file.
#include "Rte2DQuadNKS.h"

// Include Rte2D AMR Specialization header file.
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

  //------------------------------------------------------
  // Object Declarations                                 
  //------------------------------------------------------

  // batch flag
  int batch_flag;


  // Rte2D input variables and parameters:
  Rte2D_Input_Parameters Input_Parameters;


  // Multi-block solution-adaptive quadrilateral mesh 
  // solution variables.
  //
  // For Rte2D State (beam transport)
  Grid2D_Quad_Block          **MeshBlk;
  QuadTreeBlock_DataStructure  QuadTree;
  AdaptiveBlockResourceList    List_of_Global_Solution_Blocks;
  AdaptiveBlock2D_List         List_of_Local_Solution_Blocks;
  Rte2D_Quad_Block            *Local_SolnBlk;


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

  //------------------------------------------------------
  // Constructors/Destructors                             *
  //------------------------------------------------------
  Rte2DSolver () : batch_flag(0),
                   Input_Parameters(), 
		   MeshBlk(NULL), QuadTree(),
                   List_of_Global_Solution_Blocks(),
		   List_of_Local_Solution_Blocks(),
		   Local_SolnBlk(NULL), 
		   MGSolver(),
		   processor_cpu_time(), total_cpu_time(),
		   start_explicit(), end_explicit(), 
		   number_of_time_steps(0), Time(0),
		   first_step(true), NUM_VAR_RTE2D(0) {}
  ~Rte2DSolver() { DeallocateSoln(); }

  //------------------------------------------------------
  // Member functions
  // 
  // The overal solution process has been broken into
  // several modular components.  This allows one to
  // enter/exit the solution process whenever they want.
  //------------------------------------------------------

  // Allocate / deallocate solution objects
  void DeallocateSoln();

  // read input file and setup inputs
  int ReadInputFile(char *Input_File_Name_ptr, 
		    int &command_flag);

  // override input parameters with those of another object
  template <class Quad_Soln_Input_Parameters>
  void OverrideInputs(const Quad_Soln_Input_Parameters &SRC_Input_Parameters);

  // setup mesh
  int SetupMesh();
  int CopyMesh(Grid2D_Quad_Block **SRC_MeshBlk);

  // setup ics or read restart
  int SetICs();

  // broadcast solution to all processors
  int BroadcastSolution();

  // peroform initial/periodic amr
  int PerformAMR();
  int PerformPeriodicAMR();

  // perform initial/periodic morton ordering
  int PerformMortonOrdering();
  int PerformMortonReOrdering();

  // compute the solution norms
  void ComputeNorms(double &residual_l1_norm, 
		    double &residual_l2_norm, 
		    double &residual_max_norm);

  // open / write to residual file
  int WritePeriodicRestart();
  int OpenResidualFile(ofstream &);

  // output progress
  void OutputProgress(const double residual_l1_norm, 
		      const double residual_l2_norm, 
		      const double residual_max_norm,
		      ofstream &);

  // output some stats
  void OutputSolverStats(ostream &) const;

  // solvers
  int SolveMultigrid();
  int SolveExplicitTimeStep();
  int SolveSpaceMarch();
  int SolveNKS();

  // postprocess solution
  int PostProcess(int &command_flag);

  //
  // the standalone solver
  //
  int StandAloneSolve(const int batch_flag, char *Input_File_Name_ptr);

  //
  // The sequential solver
  //
  // Main setup function
  template <class Quad_Soln_Input_Parameters>
  int SetupSequentialSolve(const int batch_flag, char *Input_File_Name_ptr,
			   Grid2D_Quad_Block **SRC_MeshBlk,
			   const Quad_Soln_Input_Parameters &SRC_Input_Parameters);

  // copy src solution data over to Rte2D solution block
  template <class Quad_Soln_Block>
  void Copy_SRC_Solution_Vars(Quad_Soln_Block *SRC_Local_SolnBlk) 
  { 
    cerr<<"\n EXPLICIT SPECIALIZATION OF Copy_SRC_Solution_Vars for Rte2DQuadSolvers.h requried \n";
    exit(-1);
  }

  // copy Rte2D solution data back to src solution block
  template <class Quad_Soln_Block>
  void Copy_Rte2D_Solution_Vars(Quad_Soln_Block *SRC_Local_SolnBlk)
  { 
    cerr<<"\n EXPLICIT SPECIALIZATION OF Copy_SRC_Solution_Vars for Rte2DQuadSolvers.h requried \n";
    exit(-1);
  }
  
  // solver
  template <class Quad_Soln_Block>
  int SequentialSolve(Quad_Soln_Block *SRC_Local_SolnBlk,
		      const bool postProcess);

};

//
// TEMPLATED FUNCTION DEFINITIONS
//

/*****************************************************************************************
 *********************** MAIN CONTROL LOOPS **********************************************
 *****************************************************************************************/

/********************************************************  
 * Rte2DSolver :: SetupSequentialSolve                  *
 *                                                      *
 * The main control function which sets up the          *
 * repeated sequential solves.                          *
 ********************************************************/
template <class Quad_Soln_Input_Parameters>
inline int Rte2DSolver::SetupSequentialSolve( const int batch_flag_, 
					      char *Input_File_Name_ptr,
					      Grid2D_Quad_Block **SRC_MeshBlk,
					      const Quad_Soln_Input_Parameters &SRC_Input_Parameters ) {

  //
  // Declares
  //
  int error_flag(0);
  int command_flag(0);

  //
  // Initialization
  //
  batch_flag = batch_flag_;


  //--------------------------------------------------
  // Setup default input parameters
  //--------------------------------------------------
  error_flag = ReadInputFile(Input_File_Name_ptr, command_flag);
  if (error_flag) {
    cout << "\n Rte2D ERROR: Unable to set default inputs.\n"
	 << flush;
    return error_flag;
  } // endif

  // if we want to quit, then exit
  if (command_flag == TERMINATE_CODE) return (error_flag);

  //--------------------------------------------------
  // Override input parameters
  //--------------------------------------------------
  OverrideInputs(SRC_Input_Parameters);  

  //--------------------------------------------------
  // Copy Quadtree data structure and mesh
  //--------------------------------------------------
  error_flag = CopyMesh(SRC_MeshBlk);
  if (error_flag) {
    cout << "\n Rte2D ERROR: Unable to copy mesh.\n"
	 << flush;
    return error_flag;
  } // endif


  //--------------------------------------------------
  // Setup ICs
  //--------------------------------------------------
  error_flag = SetICs();
  if (error_flag) {
    cout << "\n Rte2D ERROR: Unable to setup ICs.\n"
	 << flush;
    return error_flag;
  } // endif


  // return the error flag
  return error_flag;

}


/********************************************************  
 * Rte2DSolver :: SequentialSolve                       *
 *                                                      *
 * The main control function which performs the         *
 * repeated sequential solves.                          *
 ********************************************************/
template <class Quad_Soln_Block>
inline int Rte2DSolver::SequentialSolve( Quad_Soln_Block *SRC_Local_SolnBlk,
					 const bool postProcess ) {

  //
  // Declares
  //
  int error_flag(0);
  int command_flag(0);


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
  // Copy Rte2D solution variables.
  //--------------------------------------------------
  // only do this if we are not prescribing an analytic field
  if (Input_Parameters.Medium_Field_Type == MEDIUM2D_FIELD_DISCRETE)
    Copy_SRC_Solution_Vars(SRC_Local_SolnBlk);


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
  error_flag = PerformAMR();
  if (error_flag) {
    cout << "\n Rte2D ERROR: Unable to perform AMR.\n"
	 << flush;
    return error_flag;
  } // endif


  //--------------------------------------------------
  // Morton Ordering
  //--------------------------------------------------
  error_flag = PerformMortonOrdering();
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
    error_flag = SolveMultigrid();
    if (error_flag) {
      cout << "\n Rte2D ERROR: Unable to perform multigrid computations.\n"
	   << flush;
      return error_flag;
    } // endif

  //
  // Space marching
  //
  } else if (Input_Parameters.i_Time_Integration == TIME_STEPPING_SPACE_MARCH) {
    error_flag = SolveSpaceMarch();
    if (error_flag) {
      cout << "\n Rte2D ERROR: Unable to perform space-marching computations.\n"
	   << flush;
      return error_flag;
    } // endif

  //
  // Non-multigrid, explicit time marching
  //
  } else {
    error_flag = SolveExplicitTimeStep();
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
    error_flag = SolveNKS();
    if (error_flag) {
      cout << "\n Rte2D ERROR: Unable to perform implicit NKS computations.\n"
	   << flush;
      return error_flag;
    } // endif
  } // endif

  //--------------------------------------------------
  // Copy over computed Rte2D solution variables.
  //--------------------------------------------------
  Copy_Rte2D_Solution_Vars(SRC_Local_SolnBlk);

  //--------------------------------------------------
  // Postprocess
  //--------------------------------------------------
  // Only postprocess if the user wants to.  
  // This could save some time at the end of each solve
  if (postProcess) {
    
    error_flag = PostProcess(command_flag);
    if (error_flag) {
      cout << "\n Rte2D ERROR: Unable to perform postprocessing.\n"
	   << flush;
      return error_flag;
    } // endif
    
    //
    // Redirect
    //
    if (command_flag == CONTINUE_CODE || 
	command_flag == EXECUTE_CODE) 
      goto continue_existing_calculation;
    
  }// endif - postProcess
  
  // return error flag
  return error_flag;


}


/*****************************************************************************************
 *********************** SETUP FUNCTIONS *************************************************
 *****************************************************************************************/

/********************************************************  
 * Rte2DSolver :: OverrideInputs                        *
 *                                                      *
 * When performing a sequential solve, override some    *
 * input parameters to ensure both the Rte2D solver     *
 * and the calling solver are in sync.                  *
 *                                                      *
 ********************************************************/
template <class Quad_Soln_Input_Parameters>
inline void Rte2DSolver::OverrideInputs( const Quad_Soln_Input_Parameters &SRC_IP ) {

  // grid related
  Input_Parameters.Number_of_Cells_Idir           = SRC_IP.Number_of_Cells_Idir;
  Input_Parameters.Number_of_Cells_Jdir           = SRC_IP.Number_of_Cells_Jdir;
  Input_Parameters.Number_of_Ghost_Cells          = SRC_IP.Number_of_Ghost_Cells;
  Input_Parameters.Number_of_Blocks_Idir          = SRC_IP.Number_of_Blocks_Idir;
  Input_Parameters.Number_of_Blocks_Jdir          = SRC_IP.Number_of_Blocks_Jdir;

  // Quadtree data structure / AMR related
  Input_Parameters.Number_of_Processors           = SRC_IP.Number_of_Processors;
  Input_Parameters.Number_of_Blocks_Per_Processor = SRC_IP.Number_of_Blocks_Per_Processor;
  Input_Parameters.Maximum_Refinement_Level       = SRC_IP.Maximum_Refinement_Level;
  Input_Parameters.Threshold_for_Refinement       = SRC_IP.Threshold_for_Refinement;
  Input_Parameters.Threshold_for_Coarsening       = SRC_IP.Threshold_for_Coarsening;

  // turn off AMR, Morton ordering
  Input_Parameters.Morton = OFF;
  Input_Parameters.AMR = OFF;
  
}


#endif // _RTE2D_QUADSOLVER__INCLUDED
