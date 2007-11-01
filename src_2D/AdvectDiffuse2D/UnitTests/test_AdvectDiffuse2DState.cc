/*!\file test_AdvectDiffuse2DState.cc
  \brief Regression tests for AdvectDiffuse2D_State class. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../New_AdvectDiffuse2DState.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_AdvectDiffuse2DState : public TestData {

    // Local variables
  public:

    // Constructor
    Data_AdvectDiffuse2DState(){
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
  typedef test_group<Data_AdvectDiffuse2DState> AdvectDiffuse2DState_TestSuite;
  typedef AdvectDiffuse2DState_TestSuite::object AdvectDiffuse2DState_object;


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
  void AdvectDiffuse2DState_object::test<1>()
  {

    set_test_name("Velocity field 1");

    VelocityFields::Set_UniformFlow_Velocities(10.0,0);
    //     Print_(VelocityFields::Uniform_Flow_Xdir(1.0,2.0));
    //     Print_(VelocityFields::Uniform_Flow_Ydir(1.0,2.0));
    //     Print_(VelocityFields::Uniform_Flow(1.0,2.0));
    //     Print_(VelocityFields::Uniform_Flow(2.0,1.0));

  }


}



// Test suite constructor
tut::AdvectDiffuse2DState_TestSuite AdvectDiffuse2DStateTestSuite("Class:AdvectDiffuse2D_State");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

