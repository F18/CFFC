/*!\file test_BGK1DVector.cc
  \brief Regression tests for template class BGK1D_Vector datatype. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../BGK1DVector.h"
#include "test_BGK1D_defines.h"

namespace tut
{

  /* Define the test-specific data class and add data members
     when tests have complex or repeating creation phase. */
  class Data_BGK1DVector : public TestData {

    // Local variables
  public:

    // Constructor
    Data_BGK1DVector(){
      if(!BGK1D_Vector::length_is_set())
	BGK1D_Vector::set_length(BGK1D_TEST_VECTOR_LENGTH);
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
  typedef test_group<Data_BGK1DVector> BGK1DVector_TestSuite;
  typedef BGK1DVector_TestSuite::object BGK1DVector_object;


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



  //end tests
}



// Test suite constructor
tut::BGK1DVector_TestSuite BGK1DVectorTestSuite("Class:BGK1D_Vector");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

