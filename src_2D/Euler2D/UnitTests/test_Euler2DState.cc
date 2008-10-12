/*!\file test_Euler2DState.cc
  \brief Regression tests for Euler2D_pState and Euler2D_cState classes. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../Euler2DState.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_Euler2D_pState : public TestData {

    // Local variables
  public:

    // Constructor
    Data_Euler2D_pState(){
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
  typedef test_group<Data_Euler2D_pState> Euler2D_pState_TestSuite;
  typedef Euler2D_pState_TestSuite::object Euler2D_pState_object;


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
  void Euler2D_pState_object::test<1>()
  {

    set_test_name("Flux in normal direction");

    Euler2D_cState Flux_n, Flux;
    Euler2D_pState W(1.5, 2.5, 3.5, 4.5);
    Vector2D normal(0.41,0.91);
    
    // Compute flux with Fx and Fy
    Flux_n = Fx(W)*normal.x + Fy(W)*normal.y;
    // Compute flux directly
    Flux = Fn(W,normal);

    ensure_distance("Flux",Flux, Flux_n, Euler2D_cState(tol));
  }

  /* Test 2:*/
  template<>
  template<>
  void Euler2D_pState_object::test<2>()
  {

    set_test_name("Min member function");

    // Set initial data
    Euler2D_pState W1(1.5, 2.5, 3.5, 4.5);
    Euler2D_pState W2(-1.5, -2.5, -3.5, -4.5);
    Euler2D_pState Wm(-23.434, 2.500001, 0.0, -4.5);

    // === check max
    ensure_equals("Min I"  , min(W1,W2), W2);
    ensure_equals("Min II" , min(W1,Wm), Euler2D_pState(-23.434, 2.5, 0.0, -4.5) );
    ensure_equals("Min III", min(W2,Wm), Euler2D_pState(-23.434, -2.5, -3.5, -4.5) );
  }

  /* Test 3:*/
  template<>
  template<>
  void Euler2D_pState_object::test<3>()
  {

    set_test_name("Max member function");

    // Set initial data
    Euler2D_pState W1(1.5, 2.5, 3.5, 4.5);
    Euler2D_pState W2(-1.5, -2.5, -3.5, -4.5);
    Euler2D_pState Wm(-23.434, 2.500001, 0.0, -4.5);

    // === check max
    ensure_equals("Max I"  , max(W1,W2), W1);
    ensure_equals("Max II" , max(W1,Wm), Euler2D_pState(1.5, 2.500001, 3.5, 4.5) );
    ensure_equals("Max III", max(W2,Wm), Euler2D_pState(-1.5, 2.500001, 0.0, -4.5) );
  }
}



// Test suite constructor
tut::Euler2D_pState_TestSuite Euler2D_pStateTestSuite("Class:Euler2D_pState");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

