/*!\file test_Euler2DExactSolutions.cc
  \brief Regression tests for Euler2D exact solution classes. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../Euler2DExactSolutions.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_Euler2DExactSolutions : public TestData {

    // Local variables
  public:

    // Constructor
    Data_Euler2DExactSolutions(){
      /* Set the global path for this test suite. 
	 It's a good practice to put it in the constructor of the data class in order to be set
	 automatically for each individual test. Declare it relative to the /src_2D directory,
	 otherwise the framework might not find the input and output files. */
      
      set_test_suite_path("Euler2D/UnitTests");
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
  typedef test_group<Data_Euler2DExactSolutions> Euler2DExactSolutions_TestSuite;
  typedef Euler2DExactSolutions_TestSuite::object Euler2DExactSolutions_object;


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
  void Euler2DExactSolutions_object::test<1>()
  {

    set_test_name("Ringleb Flow");

    Ringleb_Flow_ExactSolution_Euler RinglebSoln;
    Euler2D_pState Solution, Result;

    // First point
    Result.d   =  1.1259056146519972e+00;
    Result.v.x = -4.4645197801413531e+01;
    Result.v.y = -1.3121387860770773e+02;
    Result.p   = 9.0038627671024893e+04;

    Solution = RinglebSoln.EvaluateSolutionAt(3.0, -2.0);
    ensure_distance("Solution I", Solution, Result, AcceptedError(Result, 1.0e-12) );

    // Second point
    Result.d   = 1.1447543917607794e+00;
    Result.v.x = 8.6439529910015267e+01;
    Result.v.y = -8.9490610133323813e+01;
    Result.p   = 9.2155940123422304e+04;

    Solution = RinglebSoln.EvaluateSolutionAt(0.6, 4.0);
    ensure_distance("Solution II", Solution, Result, AcceptedError(Result, 1.0e-12) );

    // Third point
    Result.d   = 8.1141619533406362e-01; 
    Result.v.x = 6.0985898787481808e+01;
    Result.v.y = -2.9023331442404145e+02; 
    Result.p   = 5.6920306319364616e+04;

    Solution = RinglebSoln.EvaluateSolutionAt(0.8, 0.4);
    ensure_distance("Solution III", Solution, Result, AcceptedError(Result, 1.0e-12) );
  }


}



// Test suite constructor
tut::Euler2DExactSolutions_TestSuite Euler2DExactSolutionsTestSuite("Class:Euler2D_ExactSolutions");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

