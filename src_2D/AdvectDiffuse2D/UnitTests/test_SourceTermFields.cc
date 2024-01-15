/*!\file test_SourceTermFields.cc
  \brief Regression tests for class SourceTermFields. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../SourceTermFields.h"
#include "../AdvectDiffuse2DInput.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_SourceTermFields : public TestData {

    // Local variables
  public:

    // Constructor
    Data_SourceTermFields(){
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
  typedef test_group<Data_SourceTermFields> SourceTermFields_TestSuite;
  typedef SourceTermFields_TestSuite::object SourceTermFields_object;


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
  void SourceTermFields_object::test<1>()
  {
    set_test_name("Set source term in input parameters");
    set_local_input_path("SourceTermFieldsData");
    set_local_output_path("SourceTermFieldsData");

    AdvectDiffuse2D_Input_Parameters IP;

    IP.Verbose() = OFF;
    Set_Default_Input_Parameters(IP);

    // Set input file name
    this->Open_Input_File("AdvectDiffuse2D_Test.in");

    // Parse the input file
    IP.Parse_Input_File(input_file_name);

    // == check that the SourceTerm object was generated correctly
    ensure_equals("FieldRequireIntegration()", IP.SourceTerm->FieldRequireIntegration(), true);
    ensure_distance("Computed value", IP.SourceTerm->operator()(2,4,1.2), 0.03888717100122214, tol);
  }

}



// Test suite constructor
tut::SourceTermFields_TestSuite SourceTermFieldsTestSuite("Class:SourceTermFields");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

