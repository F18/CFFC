/*!\file test_Euler2D_Simulations.cc
  \brief Regression tests for Euler2D solver . */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../Euler2D/Euler2DQuad.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_Euler2DSolver : public TestData {

    // Local variables
  public:
    int error_flag;

    // Constructor
    Data_Euler2DSolver(){
      /* Set the global path for this test suite. 
	 It's a good practice to put it in the constructor of the data class in order to be set
	 automatically for each individual test. Declare it relative to the /src_2D directory,
	 otherwise the framework might not find the input and output files. */
      
      set_test_suite_path("Euler2D/UnitTests/SimulationsData");
    }

    // Solve Problem --> run Euler2DQuadSolver() subroutine with the proper input_file_name
    //                   and analyze the returned error_flag
    void Solve_Problem(int batch_flag = 1);

  private:
    
  };

  void Data_Euler2DSolver::Solve_Problem(int batch_flag){
    // Call solver
    error_flag = Euler2DQuadSolver(input_file_name,
				   batch_flag);
      
    // Check error_flag
    if (error_flag) {
      throw runtime_error("Euler2DQuadSolver() ERROR: Runtime error! For details, run verbose the regression test.");
    } /* endif */
  }


  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. 'NAME'_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_Euler2DSolver> Euler2DSolver_TestSuite;
  typedef Euler2DSolver_TestSuite::object Euler2DSolver_object;


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


  /*********************************************
     TESTING OPERATIONS:
        -->  ensure("ConditionName", test_condition);
        -->  ensure_not("ConditionName", test_condition);
        -->  ensure_equals("ConditionName", expected_value, got_value);
                 // operators '<<' and '!=' must be defined!!!
        -->  ensure_distance("ConditionName", expected_value, got_value, tolerance);
	         // operators '<<', '>=' and '<=' must be defined!!!
        -->  fail("Message");

     Obs: "ConditionName" is optional
  */



  /* Test 1:*/
  template<>
  template<>
  void Euler2DSolver_object::test<1>()
  {

    set_test_name("Better solution start-up for high-order");

    set_local_input_path("RinglebFlow/");
    set_local_output_path("RinglebFlow/");

    RunRegression = ON;
 
    // Set input file name
    this->Open_Input_File("Ringleb.in");

    // call solver
    Solve_Problem();

    if (RunRegression){

      //===== Check nodal solution
      MasterFile  = "Master_Ringleb_Flow_40x40_cpu000000.dat";
      CurrentFile = "Current_Ringleb_Flow_40x40_cpu000000.dat";
      // check
      RunRegressionTest("Nodal Solution", CurrentFile, MasterFile, 3.0e-5, 5.0e-9);

      //===== Check cell solution
      MasterFile  = "Master_Ringleb_Flow_40x40_cells_cpu000000.dat";
      CurrentFile = "Current_Ringleb_Flow_40x40_cells_cpu000000.dat";
      // check
      RunRegressionTest("Cell Solution", CurrentFile, MasterFile, 3.0e-5, 5.0e-9);
    }
  }


}



// Test suite constructor
tut::Euler2DSolver_TestSuite Euler2DSolverTestSuite("Solver:Euler2D");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

