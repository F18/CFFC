/*!\file test_Grid2DQuadIntegration.cc
  \brief Regression tests for class Grid2DQuadIntegration. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../Grid2DQuad.h"
#include "../../Math/UnitTests/TestFunctions/TestFunctions_2D.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_Grid2DQuadIntegration : public TestData {

    // Local variables
  public:

    // Constructor
    Data_Grid2DQuadIntegration(){
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
  typedef test_group<Data_Grid2DQuadIntegration> Grid2DQuadIntegration_TestSuite;
  typedef Grid2DQuadIntegration_TestSuite::object Grid2DQuadIntegration_object;


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
  void Grid2DQuadIntegration_object::test<1>()
  {

    set_test_name("Constructor");

    Grid2D_Quad_Block Grid;
    
    // Create integration object
    Grid2DQuadIntegration<Grid2D_Quad_Block> Integration(&Grid);

    // == check that the pointer is set correctly
    ensure_equals("Geometry block", Integration.getGrid(), &Grid);
  }

  /* Test 2:*/
  template<>
  template<>
  void Grid2DQuadIntegration_object::test<2>()
  {

    set_test_name("Integrate over square domain");

    Grid2D_Quad_Block **Mesh;
    
    int Blocks_Idir, Blocks_Jdir, IdirCells, JdirCells, Nghost;
    double Box_Width;
    double Result(18.916666666666667);

    Blocks_Idir = 1; Blocks_Jdir = 1;
    Box_Width = 10.0;
    IdirCells = 10; JdirCells = 10;
    Nghost = 3;

    // Set the square grid
    Mesh = Grid_Rectangular_Box(Mesh,
				Blocks_Idir, Blocks_Jdir,
				Box_Width, Box_Width,
				IdirCells, JdirCells, Nghost);

    // Create integration object
    Grid2DQuadIntegration<Grid2D_Quad_Block> Integration(&Mesh[0][0]);

    // == check integration
    ensure_distance("IntegrateFunctionOverCell()",
		    Integration.IntegrateFunctionOverCell(5,5,Test_Example1,14,Result), Result,
		    AcceptedError(Result));
  }

  /* Test 3:*/
  template<>
  template<>
  void Grid2DQuadIntegration_object::test<3>()
  {

    set_test_name("Integrate over square domain");

    Grid2D_Quad_Block **Mesh;
    
    int Blocks_Idir, Blocks_Jdir, IdirCells, JdirCells, Nghost;
    double Box_Width;
    double Result(18.916666666666667);

    Blocks_Idir = 1; Blocks_Jdir = 1;
    Box_Width = 10.0;
    IdirCells = 10; JdirCells = 10;
    Nghost = 3;

    // Set the square grid
    Mesh = Grid_Rectangular_Box(Mesh,
				Blocks_Idir, Blocks_Jdir,
				Box_Width, Box_Width,
				IdirCells, JdirCells, Nghost);

    // Create integration object
    Grid2DQuadIntegration<Grid2D_Quad_Block> Integration(&Mesh[0][0]);

    // == check integration
    ensure_distance("IntegratePolynomialOverCell()",
		    Integration.IntegratePolynomialOverCell(5,5,Test_Example1,Result), Result,
		    AcceptedError(Result));
  }

  /* Test 4:*/
  template<>
  template<>
  void Grid2DQuadIntegration_object::test<4>()
  {

    set_test_name("Use mesh to integrate over square domain");

    Grid2D_Quad_Block **Mesh;
    
    int Blocks_Idir, Blocks_Jdir, IdirCells, JdirCells, Nghost;
    double Box_Width;
    double Result(18.916666666666667);

    Blocks_Idir = 1; Blocks_Jdir = 1;
    Box_Width = 10.0;
    IdirCells = 10; JdirCells = 10;
    Nghost = 3;

    // Set the square grid
    Mesh = Grid_Rectangular_Box(Mesh,
				Blocks_Idir, Blocks_Jdir,
				Box_Width, Box_Width,
				IdirCells, JdirCells, Nghost);

    // == check integration
    ensure_distance("IntegrateFunctionOverCell()",
		    Mesh[0][0].Integration.IntegrateFunctionOverCell(5,5,Test_Example1,14,Result), Result,
		    AcceptedError(Result));
  }


}



// Test suite constructor
tut::Grid2DQuadIntegration_TestSuite Grid2DQuadIntegrationTestSuite("Template Class:Grid2DQuadIntegration");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

