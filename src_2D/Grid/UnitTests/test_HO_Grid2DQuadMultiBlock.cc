/*!\file test_HO_MultiBlock_Grid2DQuad.cc
  \brief Regression tests for 2D high-order quadrilateral multi-block grid. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../HO_Grid2DQuadMultiBlock.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_Grid2DQuadMultiBlock_HO : public TestData {

    // Local variables
  public:

    // Constructor
    Data_Grid2DQuadMultiBlock_HO(){
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
  typedef test_group<Data_Grid2DQuadMultiBlock_HO> Grid2DQuadMultiBlock_HO_TestSuite;
  typedef Grid2DQuadMultiBlock_HO_TestSuite::object Grid2DQuadMultiBlock_HO_object;


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
  void Grid2DQuadMultiBlock_HO_object::test<1>()
  {

    set_test_name("Constructor");

    // Create block grid
    Grid2D_Quad_MultiBlock_HO MultiBlockGrid;
    MultiBlockGrid.allocate(10,15);
    
    // == check 
    ensure_equals("Number of Blocks in Idir", MultiBlockGrid.Blocks_Idir(), 10);
    ensure_equals("Number of Blocks in Jdir", MultiBlockGrid.Blocks_Jdir(), 15);

    ensure_equals("Index of last block in Idir", MultiBlockGrid.Last_iBlock(), 9);
    ensure_equals("Index of last block in Jdir", MultiBlockGrid.Last_jBlock(), 14);

    ensure_equals("operator ()", &MultiBlockGrid.Grid_ptr[9][14], &MultiBlockGrid(9,14) );
  }

  /* Test 2:*/
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<2>()
  {

    set_test_name("operator = (assignment)");

    // Create block grid
    Grid2D_Quad_MultiBlock_HO MultiBlockGrid, MultiBlockGrid_Copy;
    MultiBlockGrid.allocate(1,2);
    
    // allocate single blocks
    MultiBlockGrid(0,0).allocate(10,15,2);
    MultiBlockGrid(0,1).allocate( 6,10,3);

    // operation
    MultiBlockGrid_Copy = MultiBlockGrid;
    
    // == check 
    ensure_equals("Number of Blocks in Idir", MultiBlockGrid_Copy.Blocks_Idir(), 1);
    ensure_equals("Number of Blocks in Jdir", MultiBlockGrid_Copy.Blocks_Jdir(), 2);

    ensure_equals("Index of last block in Idir", MultiBlockGrid_Copy.Last_iBlock(), 0);
    ensure_equals("Index of last block in Jdir", MultiBlockGrid_Copy.Last_jBlock(), 1);

    ensure_equals("Block(0,0) ICl", MultiBlockGrid_Copy(0,0).ICl, 2);
    ensure_equals("Block(0,0) ICu", MultiBlockGrid_Copy(0,0).ICu, 11);
    ensure_equals("Block(0,0) JCl", MultiBlockGrid_Copy(0,0).JCl, 2);
    ensure_equals("Block(0,0) JCu", MultiBlockGrid_Copy(0,0).JCu, 16);
    ensure_equals("Block(0,1) ICl", MultiBlockGrid_Copy(0,1).ICl, 3);
    ensure_equals("Block(0,1) ICu", MultiBlockGrid_Copy(0,1).ICu, 8);
    ensure_equals("Block(0,1) JCl", MultiBlockGrid_Copy(0,1).JCl, 3);
    ensure_equals("Block(0,1) JCu", MultiBlockGrid_Copy(0,1).JCu, 12);
  }

  /* Test 3:*/
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<3>()
  {

    set_test_name("operator = (self-assignment)");

    // Create block grid
    Grid2D_Quad_MultiBlock_HO MultiBlockGrid;
    MultiBlockGrid.allocate(1,2);
    
    // allocate single blocks
    MultiBlockGrid(0,0).allocate(10,15,2);
    MultiBlockGrid(0,1).allocate( 6,10,3);

    // operation
    MultiBlockGrid = MultiBlockGrid;
    
    // == check 
    ensure_equals("Number of Blocks in Idir", MultiBlockGrid.Blocks_Idir(), 1);
    ensure_equals("Number of Blocks in Jdir", MultiBlockGrid.Blocks_Jdir(), 2);

    ensure_equals("Index of last block in Idir", MultiBlockGrid.Last_iBlock(), 0);
    ensure_equals("Index of last block in Jdir", MultiBlockGrid.Last_jBlock(), 1);

    ensure_equals("Block(0,0) ICl", MultiBlockGrid(0,0).ICl, 2);
    ensure_equals("Block(0,0) ICu", MultiBlockGrid(0,0).ICu, 11);
    ensure_equals("Block(0,0) JCl", MultiBlockGrid(0,0).JCl, 2);
    ensure_equals("Block(0,0) JCu", MultiBlockGrid(0,0).JCu, 16);
    ensure_equals("Block(0,1) ICl", MultiBlockGrid(0,1).ICl, 3);
    ensure_equals("Block(0,1) ICu", MultiBlockGrid(0,1).ICu, 8);
    ensure_equals("Block(0,1) JCl", MultiBlockGrid(0,1).JCl, 3);
    ensure_equals("Block(0,1) JCu", MultiBlockGrid(0,1).JCu, 12);
  }

  /* Test 4:*/
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<4>()
  {

    set_test_name("operator = (re-assignment)");

    // Create block grid
    Grid2D_Quad_MultiBlock_HO MultiBlockGrid_I, MultiBlockGrid_II, TestedGrid;
    
    // allocate MultiBlockGrid_I
    MultiBlockGrid_I.allocate(1,2);
    MultiBlockGrid_I(0,0).allocate(10,15,2);
    MultiBlockGrid_I(0,1).allocate( 6,10,3);

    // allocate MultiBlockGrid_II
    MultiBlockGrid_II.allocate(2,1);
    MultiBlockGrid_II(0,0).allocate(10,15,2);
    MultiBlockGrid_II(1,0).allocate( 6,10,3);

    // operation I
    TestedGrid = MultiBlockGrid_I;

    // == check operation I
    ensure_equals("Number of Blocks in Idir", TestedGrid.Blocks_Idir(), 1);
    ensure_equals("Number of Blocks in Jdir", TestedGrid.Blocks_Jdir(), 2);

    ensure_equals("Index of last block in Idir", TestedGrid.Last_iBlock(), 0);
    ensure_equals("Index of last block in Jdir", TestedGrid.Last_jBlock(), 1);

    ensure_equals("Block(0,0) ICl", TestedGrid(0,0).ICl, 2);
    ensure_equals("Block(0,0) ICu", TestedGrid(0,0).ICu, 11);
    ensure_equals("Block(0,0) JCl", TestedGrid(0,0).JCl, 2);
    ensure_equals("Block(0,0) JCu", TestedGrid(0,0).JCu, 16);
    ensure_equals("Block(0,1) ICl", TestedGrid(0,1).ICl, 3);
    ensure_equals("Block(0,1) ICu", TestedGrid(0,1).ICu, 8);
    ensure_equals("Block(0,1) JCl", TestedGrid(0,1).JCl, 3);
    ensure_equals("Block(0,1) JCu", TestedGrid(0,1).JCu, 12);

    // operation II
    TestedGrid = MultiBlockGrid_II;

    // == check operation II
    ensure_equals("Number of Blocks in Idir", TestedGrid.Blocks_Idir(), 2);
    ensure_equals("Number of Blocks in Jdir", TestedGrid.Blocks_Jdir(), 1);

    ensure_equals("Index of last block in Idir", TestedGrid.Last_iBlock(), 1);
    ensure_equals("Index of last block in Jdir", TestedGrid.Last_jBlock(), 0);

    ensure_equals("Block(0,0) ICl", TestedGrid(0,0).ICl, 2);
    ensure_equals("Block(0,0) ICu", TestedGrid(0,0).ICu, 11);
    ensure_equals("Block(0,0) JCl", TestedGrid(0,0).JCl, 2);
    ensure_equals("Block(0,0) JCu", TestedGrid(0,0).JCu, 16);
    ensure_equals("Block(1,0) ICl", TestedGrid(1,0).ICl, 3);
    ensure_equals("Block(1,0) ICu", TestedGrid(1,0).ICu, 8);
    ensure_equals("Block(1,0) JCl", TestedGrid(1,0).JCl, 3);
    ensure_equals("Block(1,0) JCu", TestedGrid(1,0).JCu, 12); 
  }

  /* Test 5:*/
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<5>()
  {

    set_test_name("Check input-output operators");

    // Create block grid
    Grid2D_Quad_MultiBlock_HO MultiBlockGrid;
    MultiBlockGrid.allocate(10,15);
    
    // == check 
    Check_Input_Output_Operator("MultiBlockGrid variable", MultiBlockGrid);
  }


}



// Test suite constructor
tut::Grid2DQuadMultiBlock_HO_TestSuite Grid2DQuadMultiBlock_HOTestSuite("Class:Grid2DQuadMultiBlock_HO");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

