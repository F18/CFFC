/*!\file test_SourceTermFields.cc
  \brief Regression tests for class SourceTermFields. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../SourceTermFields.h"

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
    set_test_name("Test Basic field type");

    SourceFieldBasicType * FieldS;

    FieldS = new ZERO_SourceField;

    Print_(FieldS->FieldSoln(10.0));
    Print_(FieldS->FieldRequireIntegration());

    delete FieldS;

    FieldS = new Linear_SourceField;

    Print_(FieldS->FieldSoln(12.0));
    Print_(FieldS->FieldRequireIntegration());

    delete FieldS;

    FieldS = new Exponential_SourceField;

    Print_(FieldS->FieldSoln(2.0,4.0,1.2));
    Print_(FieldS->FieldRequireIntegration());

    delete FieldS;


    Print_("Here")

    SourceTermFields& myField = SourceTermFields::getInstance();
    myField.SetSource(3);
    if(myField.Source->FieldRequireIntegration()){
      Print_(myField.Source->FieldSoln(2,4,1.2));
    } else {
      Print_(myField.Source->FieldSoln(12.0));
    }
  }

  /* Test 2:*/
  template<>
  template<>
  void SourceTermFields_object::test<2>()
  {
    set_test_name("Input param");

    AdvectDiffuse2D_Input_Parameters IP;
    strcpy(IP.Next_Control_Parameter, "Source_Linear_Tau_Coeff");
    int i_command = INVALID_INPUT_CODE;

    SourceTermFields& myField = SourceTermFields::getInstance();
    myField.SetSource(2);
    myField.Source->Parse_Next_Input_Control_Parameter(IP,i_command);
    if(myField.Source->FieldRequireIntegration()){
      Print_(myField.Source->FieldSoln(2,4,1.2));
    } else {
      Print_(myField.Source->FieldSoln(12.0));
    }

  }

  /* Test 3:*/
  template<>
  template<>
  void SourceTermFields_object::test<3>()
  {

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

