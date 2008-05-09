/******************* test_Rte2DSolver.cc ******************************

  Test file for problems solved with Rte2D

***********************************************************************/

// Include required C++ libraries.
// None

// Using std namespace functions
// None

// Include CFFC header files
#include "../../UnitTesting/TestData.h"
#include "../Rte2D/Rte2DQuadSolvers.h"

using namespace std;

namespace tut
{

  /*************************************************************
   * Class Definition:  Data_Rte2DSolver                       *
   *                                                           *
   * Defines the test data class for Rte2D solvers.            *
   *                                                           *
   *************************************************************/
  class Data_Rte2DSolver: public TestData {
    // Public data and member functions
  public:

    //
    // Member functions
    //
    int error_flag;
    int verbose;
    
    //
    // Member functions
    //

    // Default Constructor
    Data_Rte2DSolver(void);

    // Destructor
    ~Data_Rte2DSolver(void);

    // Solve Problem --> run Rte2DQuadSolver() subroutine with the proper input_file_name
    //                   and analyze the return error_flag
    void Solve_Problem(int batch_flag = 1);

  };

  /*************************************************************
   * Constructor / Destructor                                  *
   *************************************************************/
  Data_Euler2DSolver::Data_Euler2DSolver(void){

    set_local_output_path("test_Rte2D_Solver");
    set_local_input_path("test_Rte2D_Solver");

    verbose = 0; // ZERO value is verbose
  }

  Data_Euler2DSolver::~Data_Euler2DSolver(void){
  }

  /*************************************************************
   * Main call to solve a problem                              *
   *************************************************************/
  void Data_Euler2DSolver::Solve_Problem(int batch_flag){

    //
    // Call solver, throw any errors that are caught
    //
    try{

      // Call solver
      error_flag = Euler2DQuadSolver(input_file_name,
				     batch_flag);
      
      // Check error_flag
      if (error_flag) {
	CFDkit_Finalize_MPI();
	throw runtime_error("Rte2DQuadSolver() ERROR: Runtime error! For details, run verbose the regression test.");
      } // endif

    //
    // Catch errors thro wn
    //
    } catch(...){
      CFDkit_Finalize_MPI();
      throw runtime_error("Rte2DQuadSolver() ERROR: Runtime error! For details, run verbose the regression test.");
    }
  }


  /*************************************************************
   * Typedefs for Test Group                                   *
   *          -- tg = Rte2DSolver_TestGroup                    *
   *          -- object = Rte2DSolver_TestObject               *
   *************************************************************/
  typedef test_group<Data_Rte2DSolver> tg;
  typedef tg::object object;
}

/******************************************************************************************************
 ******************************************************************************************************
 *                      ********    *******        ***       ********        ***                      *
 *                         **       **           **   **        **         **   **                    *
 *                         **       **          **              **        **                          *
 *                         **       *******      *****          **          *****                     *
 *                         **       **                **        **               **                   *
 *                         **       **            *    **       **          *   **                    *
 *                         **       *******        ****         **           ****                     *
 ******************************************************************************************************
 ******************************************************************************************************/

namespace tut
{
  /*************************************************************
   * Test 1                                                    *
   *************************************************************/
  template<>
  template<>
  void object::test<1>()
  {
    // simulation settings
    set_test_name("Simulation_001");
    RunRegression = ON;
 
    // create a new solver object
    Rte2DSolver SolverObj;

    // Setup default input parameters
    Set_Default_Input_Parameters(SolverObj.Input_Parameters);

    // call solver
    Solve_Problem();

  }




}

/*************************************************************
 * Register the test                                         *
 *************************************************************/
tut::tg Rte2DSolver_Test("Rte2D_Simulations");
