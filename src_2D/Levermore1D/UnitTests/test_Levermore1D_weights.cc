/*!\file test_Levermore1D_weights.cc
  \brief Regression tests for template class Levermore1D_weights datatype. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../Levermore1DState.h"
#include "test_Levermore1D_defines.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_Levermore1D_weights : public TestData {

    // Local variables
  public:

    // Constructor
    Data_Levermore1D_weights(){
      if(!Levermore1D_Vector::length_is_set())
	Levermore1D_Vector::set_length(LEVERMORE1D_VECTOR_LENGTH);
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
  typedef test_group<Data_Levermore1D_weights> Levermore1D_weights_TestSuite;
  typedef Levermore1D_weights_TestSuite::object Levermore1D_weights_object;


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
  void Levermore1D_weights_object::test<1>()
  {
    set_test_name("Copy Constructors");

    int i(0);
    Levermore1D_weights U1;
    
    for(i=1;i<=Levermore1D_Vector::get_length();++i) {
      //set to a value for copy constructor test
      U1[i] = pow((double)i,3.54321) / 18.765;
    }

    Levermore1D_weights U2(U1);

    for(i=1;i<=Levermore1D_Vector::get_length();++i) {
      ensure_distance("U2=U1",U2[i],U1[i],fabs(U1[i])*tol);
    }

  }

  //end tests
}



// Test suite constructor
tut::Levermore1D_weights_TestSuite Levermore1D_weightsTestSuite("Class:Levermore1D_weights");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

