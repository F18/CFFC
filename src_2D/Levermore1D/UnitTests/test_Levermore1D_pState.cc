/*!\file test_Levermore1D_pState.cc
  \brief Regression tests for template class Levermore1D_pState datatype. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../Levermore1DState.h"

/* define useful constants for tests */
#define LEVERMORE1D_VECTOR_LENGTH    5

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_Levermore1D_pState : public TestData {

    // Local variables
  public:

    // Constructor
    Data_Levermore1D_pState(){ }

  private:
    
  };

  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. 'NAME'_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_Levermore1D_pState> Levermore1D_pState_TestSuite;
  typedef Levermore1D_pState_TestSuite::object Levermore1D_pState_object;


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
  void Levermore1D_pState_object::test<1>()
  {
    set_test_name("Constructors");

    int i(0);
    Levermore1D_pState<LEVERMORE1D_VECTOR_LENGTH> P1;
    
    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      ensure_distance("Default Constructor, set it to zero.",P1[i],0.0,tol);
      //set to a value for copy constructor test
      P1[i] = pow((double)i,1.23456) / 98.765;
    }

    Levermore1D_pState<LEVERMORE1D_VECTOR_LENGTH> P2(P1);

    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      ensure_distance("P2=P1",P2[i],P1[i],fabs(P1[i])*tol);
    }

  }


  //end tests
}



// Test suite constructor
tut::Levermore1D_pState_TestSuite Levermore1D_pStateTestSuite("Template Class:Levermore1D_pState");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

