/*!\file test_HO_Spline2D.cc
  \brief Regression tests for class Spline2D_HO. */

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
  class Data_Spline2D_HO : public TestData {

    // Local variables
  public:

    // Constructor
    Data_Spline2D_HO(){
      /* Set the global path for this test suite. 
	 It's a good practice to put it in the constructor of the data class in order to be set
	 automatically for each individual test. Declare it relative to the /src_2D directory,
	 otherwise the framework might not find the input and output files. */
      
      set_test_suite_path("Grid/UnitTests");
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
  typedef test_group<Data_Spline2D_HO> Spline2D_HO_TestSuite;
  typedef Spline2D_HO_TestSuite::object Spline2D_HO_object;


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
  void Spline2D_HO_object::test<1>()
  {

    set_test_name("Name of the test");

  }


}



// Test suite constructor
tut::Spline2D_HO_TestSuite Spline2D_HOTestSuite("Class:Spline2D_HO");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

