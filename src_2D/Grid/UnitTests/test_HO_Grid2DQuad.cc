/*!\file test_HO_Grid2DQuad.cc
  \brief Regression tests for 2D high-order quadrilateral block grid. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../HO_Grid2DQuad.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_Grid2DQuad_HO : public TestData {

    // Local variables
  public:

    // Constructor
    Data_Grid2DQuad_HO(){
      /* Set the global path for this test suite. 
	 It's a good practice to put it in the constructor of the data class in order to be set
	 automatically for each individual test. Declare it relative to the /src_2D directory,
	 otherwise the framework might not find the input and output files. */
      
      set_test_suite_path("Grid/UnitTests/");
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
  typedef test_group<Data_Grid2DQuad_HO> Grid2DQuad_HO_TestSuite;
  typedef Grid2DQuad_HO_TestSuite::object Grid2DQuad_HO_object;


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
  void Grid2DQuad_HO_object::test<1>()
  {

    set_test_name("Constructor");

    // Create block grid
    Grid2D_Quad_Block_HO Grid;
    Grid.allocate(10,15,2);
    
    // == check 
    ensure_equals("NNi", Grid.NNi, 15);
    ensure_equals("NNj", Grid.NNj, 20);
    ensure_equals("NCi", Grid.NCi, 14);
    ensure_equals("NCj", Grid.NCj, 19);
    ensure_equals("Nghost", Grid.Nghost, 2);
    ensure_equals("ICl", Grid.ICl, 2);
    ensure_equals("ICu", Grid.ICu, 11);
    ensure_equals("JCl", Grid.JCl, 2);
    ensure_equals("JCu", Grid.JCu, 16);
    ensure_equals("INl", Grid.INl, 2);
    ensure_equals("INu", Grid.INu, 12);
    ensure_equals("JNl", Grid.JNl, 2);
    ensure_equals("JNu", Grid.JNu, 17);
  }


}



// Test suite constructor
tut::Grid2DQuad_HO_TestSuite Grid2DQuad_HOTestSuite("Class:Grid2DQuad_HO");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

