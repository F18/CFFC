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
    //Initialize grid and set variables to certain predefined values in order to be verified.
    void InitializeGrid_Algorithm1(Grid2D_Quad_Block_HO &Grid, int Ni, int Nj, int Ng);

  private:
    
  };

  //Initialize grid and set variables to certain predefined values in order to be verified.
  void Data_Grid2DQuad_HO::InitializeGrid_Algorithm1(Grid2D_Quad_Block_HO &Grid, int Ni, int Nj, int Ng){

    // Initialize grid
    Grid.allocate(Ni,Nj,Ng);
           
    Grid.SminN = 10; 
    Grid.SmaxN = 9; 
    Grid.SminS = 8; 
    Grid.SmaxS = 7; 
    Grid.SminW = 6; 
    Grid.SmaxW = 5; 
    Grid.SminE = 4; 
    Grid.SmaxE = 3; 

    Grid.StretchI = 10;
    Grid.StretchJ = 10;
    Grid.BetaI = 10;
    Grid.TauI = 10;
    Grid.BetaJ = 10;
    Grid.TauJ = 10;
    Grid.OrthogonalN = 10;
    Grid.OrthogonalS = 10;
    Grid.OrthogonalE = 10;
    Grid.OrthogonalW = 10;

    //Initialize all cells including ghost cells.
    //The area for the [i,j] cell is i*10+j and the 
    //location of the cell center is (i,j)
    Vector2D CellLocation;
    for (int i = 0; i < Ni+2*Ng; ++i){
      for (int j = 0; j < Nj+2*Ng; ++j){
	CellLocation = Vector2D(i,j);
	Grid.Cell[i][j].Xc = CellLocation; 
	Grid.Cell[i][j].A = (i*10+j);
	Grid.Cell[i][j].I = i;
	Grid.Cell[i][j].J = j;
      }
    }
 
    //Initialize values for each node.
    //The location of the node is (i,j).
    Vector2D NodeLocation;
    for (int i = 0; i < Ni+2*Ng+1; ++i){
      for (int j = 0; j < Nj+2*Ng+1; ++j){
	NodeLocation = Vector2D(i,j);
	Grid.Node[i][j].X = NodeLocation; 
      }
    }
    
  }


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

  /* Test 2:*/
  template<>
  template<>
  void Grid2DQuad_HO_object::test<2>()
  {
    
    
    set_test_name("operator = (assignment)");

    // Create block grid
    Grid2D_Quad_Block_HO Grid, Grid_Copy;
    
    int Ni = 10;
    int Nj = 15;
    int Nghost = 2;
    InitializeGrid_Algorithm1(Grid, Ni, Nj, Nghost);


    // operation
    Grid_Copy = Grid;

    // == check 
    ensure_equals("NNi", Grid_Copy.NNi, 15);
    ensure_equals("NNj", Grid_Copy.NNj, 20);
    ensure_equals("NCi", Grid_Copy.NCi, 14);
    ensure_equals("NCj", Grid_Copy.NCj, 19);
    ensure_equals("Nghost", Grid_Copy.Nghost, 2);
    ensure_equals("ICl", Grid_Copy.ICl, 2);
    ensure_equals("ICu", Grid_Copy.ICu, 11);
    ensure_equals("JCl", Grid_Copy.JCl, 2);
    ensure_equals("JCu", Grid_Copy.JCu, 16);
    ensure_equals("INl", Grid_Copy.INl, 2);
    ensure_equals("INu", Grid_Copy.INu, 12);
    ensure_equals("JNl", Grid_Copy.JNl, 2);
    ensure_equals("JNu", Grid_Copy.JNu, 17);
    
    ensure_equals("SminN",Grid_Copy.SminN, 10);
    ensure_equals("SmaxN",Grid_Copy.SmaxN, 9);
    ensure_equals("SminS",Grid_Copy.SminS, 8);
    ensure_equals("SmaxS",Grid_Copy.SmaxS, 7);
    ensure_equals("SminW",Grid_Copy.SminW, 6);
    ensure_equals("SmaxW",Grid_Copy.SmaxW, 5);
    ensure_equals("SminE",Grid_Copy.SminE, 4);
    ensure_equals("SmaxE",Grid_Copy.SmaxE, 3);

    ensure_equals("StretchI",Grid_Copy.StretchI, 10);
    ensure_equals("StretchJ",Grid_Copy.StretchJ, 10);
    ensure_equals("BetaI",Grid_Copy.BetaI, 10);
    ensure_equals("TauI",Grid_Copy.TauI, 10);
    ensure_equals("BetaJ",Grid_Copy.BetaJ, 10);
    ensure_equals("TauJ",Grid_Copy.TauJ, 10);
    ensure_equals("OrthogonalN",Grid_Copy.OrthogonalN, 10);
    ensure_equals("OrthogonalS",Grid_Copy.OrthogonalS, 10);
    ensure_equals("OrthogonalE",Grid_Copy.OrthogonalE, 10);
    ensure_equals("OrthogonalW",Grid_Copy.OrthogonalW, 10);

    //Loop through all cells including the ghost cells to ensure
    //that the values are copied correctly.
    //The expected area for the [i,j] cell is i*10+j and the 
    //location of the cell center is (i,j)
    Vector2D CellLocation;
    for (int i = 0; i < Ni+2*Nghost; ++i){
      for (int j = 0; j < Nj+2*Nghost; ++j){
	CellLocation = Vector2D(i,j);
	ensure_equals("Cell.Xc", Grid_Copy.Cell[i][j].Xc, CellLocation);
	ensure_equals("Cell.A", Grid_Copy.Cell[i][j].A, i*10+j);
	ensure_equals("Cell.I", Grid_Copy.Cell[i][j].I, i);
	ensure_equals("Cell.J", Grid_Copy.Cell[i][j].J, j);
      }
    }

    //Loop through all nodes to ensure that the locations are 
    //copied correctly.
    //The expected location of the node is (i,j).
    Vector2D NodeLocation;
    for (int i = 0; i < Ni+2*Nghost+1; ++i){
      for (int j = 0; j < Nj+2*Nghost+1; ++j){
	NodeLocation = Vector2D(i,j);
	ensure_equals("Node.X", Grid_Copy.Node[i][j].X, NodeLocation);
      }
    }

  }

  /* Test 3:*/
  template<>
  template<>
  void Grid2DQuad_HO_object::test<3>()
  {

    set_test_name("Check input-output operators");

    // Create block grid
    Grid2D_Quad_Block_HO Grid;
    InitializeGrid_Algorithm1(Grid, 10, 15, 2);

    // == check 
    Check_Input_Output_Operator("Grid variable", Grid);
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

