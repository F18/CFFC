/*!\file test_HO_Cell2D.cc
  \brief Regression tests for Cell2D_HO class. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../HO_Cell2D.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_Cell2D_HO : public TestData {

    // Local variables
  public:

    // Constructor
    Data_Cell2D_HO(){
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
  typedef test_group<Data_Cell2D_HO> Cell2D_HO_TestSuite;
  typedef Cell2D_HO_TestSuite::object Cell2D_HO_object;


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
  void Cell2D_HO_object::test<1>()
  {

    set_test_name("Cell constructor");

    Cell2D_HO Cell(3);
    double Result;

    // Initialize Cell
    Cell.I = 1; Cell.J = 2; Cell.A = 12.1212;
    Cell.Xc.x = 0.1212; Cell.Xc.y = 1.2323;
    for(int i = 0; i<=Cell.GeomCoeff().LastElem(); ++i){
      Cell.GeomCoeffValue(i) = (i+3)*0.34;
    }

    // Check geometric coefficients
    for(int i = 0; i<=Cell.GeomCoeff().LastElem(); ++i){
      Result = (i+3)*0.34;
      ensure_equals("Check geometric coeff. value", Cell.GeomCoeffValue(i), Result);
    }

    // Check different field access
    ensure_equals("P1", Cell.GeomCoeff(9).P1(), 3);
    ensure_equals("P2", Cell.GeomCoeff(9).P2(), 0);
    ensure_equals("Value", Cell.GeomCoeff(9).D(), Result);
    ensure_equals("Value with powers", Cell.GeomCoeffValue(3,0), Result);
  }


  /* Test 2:*/
  template<>
  template<>
  void Cell2D_HO_object::test<2>()
  {

    set_test_name("Output operator <<");
    set_local_output_path("HO_Cell2D");
    set_local_input_path("HO_Cell2D");

    RunRegression = ON;

    Cell2D_HO Cell(3);

    // Initialize Cell
    Cell.I = 1; Cell.J = 2; Cell.A = 12.1212;
    Cell.Xc.x = 0.1212; Cell.Xc.y = 1.2323;
    for(int i = 0; i<=Cell.GeomCoeff().LastElem(); ++i){
      Cell.GeomCoeffValue(i) = (i+3)*0.34;
    }

    MasterFile = "HighOrderCell.dat";
    CurrentFile ="Current_HighOrderCell.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // write new data
      out() << Cell;

      // check
      RunRegressionTest("operator <<", CurrentFile, MasterFile, 5.0e-9, 5.0e-9);

    } else {
      Open_Output_File(MasterFile);

      // write data
      out() << Cell;
    }
  }

  /* Test 3:*/
  template<>
  template<>
  void Cell2D_HO_object::test<3>()
  {

    set_test_name("Output operator <<");

    Cell2D_HO Cell(3);

    // Initialize Cell
    Cell.I = 1; Cell.J = 2; Cell.A = 12.1212;
    Cell.Xc.x = 0.1212; Cell.Xc.y = 1.2323;
    for(int i = 0; i<=Cell.GeomCoeff().LastElem(); ++i){
      Cell.GeomCoeffValue(i) = (i+3)*0.34;
    }

    Check_Input_Output_Operator(Cell);
  }



}



// Test suite constructor
tut::Cell2D_HO_TestSuite Cell2D_HOTestSuite("Class:Cell2D_HO");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

