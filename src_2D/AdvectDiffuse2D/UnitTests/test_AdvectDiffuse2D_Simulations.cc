/*!\file test_AdvectDiffuse2D_Simulations.cc
  \brief Regression tests for 2D advection-diffusion solver. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_AdvectDiffuse2DSolver : public TestData {

    // Local variables
  public:

    // Constructor
    Data_AdvectDiffuse2DSolver(){
      /* Set the global path for this test suite. 
	 It's a good practice to put it in the constructor of the data class in order to be set
	 automatically for each individual test. Declare it relative to the /src_2D directory,
	 otherwise the framework might not find the input and output files. */
      
      set_test_suite_path("AdvectDiffuse2D/UnitTests/");
    }

  private:
    
  };

  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. 'NAME'_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_AdvectDiffuse2DSolver> AdvectDiffuse2DSolver_TestSuite;
  typedef AdvectDiffuse2DSolver_TestSuite::object AdvectDiffuse2DSolver_object;


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
  void AdvectDiffuse2DSolver_object::test<1>()
  {

    set_test_name("Name of the test");

  }


}



// Test suite constructor
tut::AdvectDiffuse2DSolver_TestSuite AdvectDiffuse2DSolverTestSuite("Solver: AdvectDiffuse2D");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

