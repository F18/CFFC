/*!\file test_Euler1D_Solver.cc
  \brief Regression tests for Euler1D solver . */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../Euler1D.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_Euler1DSolver : public TestData {

    // Local variables
  public:

    int error_flag;

    // Member functions

    // Constructor
    Data_Euler1DSolver(){
      /* Set the global path for this test suite. 
	 It's a good practice to put it in the constructor of the data class in order to be set
	 automatically for each individual test. Declare it relative to the /src_2D directory,
	 otherwise the framework might not find the input and output files. */
      
      set_test_suite_path("Euler1D/UnitTests");
    }

    ~Data_Euler1DSolver(){};

    // Solve Problem --> run Euler1DSolver() subroutine with the proper input_file_name
    //                   and analyze the returned error_flag
    void Solve_Problem(int batch_flag = 1);

  private:
    
  };

  void Data_Euler1DSolver::Solve_Problem(int batch_flag){
    // Call solver
    error_flag = Euler1DSolver(input_file_name,
			       batch_flag);
      
    // Check error_flag
    if (error_flag) {
      throw runtime_error("Euler1DSolver() ERROR: Runtime error! For details, run verbose the regression test.");
    } /* endif */
  }


  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. Euler1DSolver_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_Euler1DSolver> Euler1DSolver_TestSuite;
  typedef Euler1DSolver_TestSuite::object Euler1DSolver_object;


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
  void Euler1DSolver_object::test<1>()
  {
    set_test_name("Sod Problem 2nd-order");
    set_local_output_path("SodProblem");
    set_local_input_path("SodProblem");

    RunRegression = ON;
    
    // Set input file name
    Open_Input_File("sod.in");

    // call solver
    Solve_Problem();

    if (RunRegression){

      //===== Check solution
      MasterFile  = "sod.dat";
      CurrentFile = "Current_sod.dat";

      // check
      RunRegressionTest("Sod Solution", CurrentFile, MasterFile, 5.0e-12, 5.0e-12);
    }
  }

  /* Test 2:*/
  template<>
  template<>
  void Euler1DSolver_object::test<2>()
  {
    set_test_name("Sod Problem 4th-order with CENO");
    set_local_output_path("SodProblem");
    set_local_input_path("SodProblem");

    RunRegression = ON;
    
    // Set input file name
    Open_Input_File("sod_CENO_HighOrder_Simulation.in");

    // call solver
    Solve_Problem();
    
    if (RunRegression){

      //===== Check solution
      MasterFile  = "sod_CENO_HighOrder_Simulation.dat";
      CurrentFile = "Current_sod_CENO_HighOrder_Simulation.dat";

      // check
      RunRegressionTest("CENO high-order Sod solution", CurrentFile, MasterFile, 5.0e-12, 5.0e-12);
    }
  }

  /* Test 3:*/
  template<>
  template<>
  void Euler1DSolver_object::test<3>()
  {
    set_test_name("Sod Problem 4th-order with ENO characteristics");
    set_local_output_path("SodProblem");
    set_local_input_path("SodProblem");

    RunRegression = ON;
    
    // Set input file name
    Open_Input_File("sod_HighOrder_ENO.in");

    // call solver
    Solve_Problem();
    
    if (RunRegression){

      //===== Check solution
      MasterFile  = "sod_HighOrder_ENO.dat";
      CurrentFile = "Current_sod_HighOrder_ENO.dat";

      // check
      RunRegressionTest("ENO with charcteristics high-order Sod solution", CurrentFile, MasterFile, 5.0e-12, 5.0e-12);
    }
  }


}



// Test suite constructor
tut::Euler1DSolver_TestSuite Euler1DSolverTestSuite("Solver: Euler1D");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

