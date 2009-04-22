/*!\file test_Grid2DQuadIntegration.cc
  \brief Regression tests for class Grid2DQuadIntegration. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../Grid2DQuad.h"
#include "../HO_Grid2DQuadMultiBlock.h"
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

  /* Test 5:*/
  template<>
  template<>
  void Grid2DQuadIntegration_object::test<5>()
  {

    set_test_name("Use mesh to integrate over square domain, high-order too");

    Grid2D_Quad_MultiBlock_HO Mesh;
    
    int Blocks_Idir, Blocks_Jdir, IdirCells, JdirCells, Nghost, RecOrder;
    int i,j;
    double Box_Width;
    double Result;

    Blocks_Idir = 1; Blocks_Jdir = 1;
    Box_Width = 10.0;
    IdirCells = 10; JdirCells = 10;
    Nghost = 3;
    RecOrder = 1;

    // Require high-order boundaries.
    Grid2D_Quad_Block_HO::setHighOrderBoundaryRepresentation();

    // Set the square grid
    Mesh.Grid_Rectangular_Box(Blocks_Idir, Blocks_Jdir,
			      Box_Width, Box_Width,
			      IdirCells, JdirCells, Nghost,
			      RecOrder);

    // == check integration
    for (i = 0; i<=Mesh(0,0).ICu + Mesh(0,0).Nghost; ++i){
      for (j = 0; j<=Mesh(0,0).JCu + Mesh(0,0).Nghost; ++j){
	// Calculate analytic result
	Result = Test_Example1_Integral(Mesh(0,0).nodeSW(i,j).x(),
					Mesh(0,0).nodeSE(i,j).x(),
					Mesh(0,0).nodeSE(i,j).y(),
					Mesh(0,0).nodeNE(i,j).y());

	// == check numerical integration
	ensure_distance("IntegrateFunctionOverCell()",
			Mesh(0,0).Integration.IntegrateFunctionOverCell(i,j,
									Test_Example1,
									Test_Example1_XDependencyIntegrated,
									14,Result),
			Result,
			AcceptedError(Result));
      }	// endfor
    } // endfor

  }

  /* Test 6:*/
  template<>
  template<>
  void Grid2DQuadIntegration_object::test<6>()
  {

    set_test_name("CalculateFunctionIntegralWithGaussQuadratures()");

    Grid2D_Quad_MultiBlock_HO Mesh;
    
    int i,j;
    double Result;

    int Blocks_Idir = 1;
    int Blocks_Jdir = 1;
    double Box_Width = 10.0;
    int IdirCells = 10;
    int JdirCells = 10;
    int Nghost = 2;
    int RecOrder = 1;
    int Stretching_Flag = 0;
    int Stretching_Type_Idir = 1;
    int Stretching_Type_Jdir = 1;
    double Stretching_Factor_Idir = 1.0;
    double Stretching_Factor_Jdir = 1.0;
    Vector2D SW = Vector2D(-2.0,-2.0);
    Vector2D SE = Vector2D(5.0, -0.5);
    Vector2D NE = Vector2D(3.5, 2.7);
    Vector2D NW = Vector2D(-1.0, 2.5);

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    // Set the square grid
    Mesh.Grid_Deformed_Box(Blocks_Idir, Blocks_Jdir,
			   SW, SE, NE, NW,
			   Stretching_Flag,
			   Stretching_Type_Idir,
			   Stretching_Type_Jdir,
			   Stretching_Factor_Idir,
			   Stretching_Factor_Jdir,
			   IdirCells, JdirCells, Nghost,
			   RecOrder);

    int NumGQP = Spline2DInterval_HO::get_NumGQPoints_ContourIntegral();
    double Val(0);
    Vector2D *GaussQuadPoints = new Vector2D [NumGQP];   // the GQPs used to calculate the integral along a line segment
    double * GaussQuadWeights = new double [NumGQP];     // the Gauss integration weights for each Gauss quadrature
    double DeltaY;


    // === Check with East face of cell (1,2)

    // get the Gauss quadrature points
    Mesh(0,0).getGaussQuadPointsFaceE(1,2,GaussQuadPoints,NumGQP);

    /* set the GaussQuadWeights. */
    GaussQuadratureData::getGaussQuadWeights(GaussQuadWeights, NumGQP);

    // get DeltaY
    DeltaY = (Mesh(0,0).nodeNE(1,2).y() - Mesh(0,0).nodeSE(1,2).y());

    // integrate along the line segment
    Val = Mesh(0,0).Integration.CalculateFunctionIntegralWithGaussQuadratures(Test_Example1_XDependencyIntegrated,
									      GaussQuadPoints, GaussQuadWeights,
									      NumGQP, DeltaY, Val);

    Result = -5.420681250000002;

    // == check integration
    ensure_distance("Integration along line segment",Val, Result, AcceptedError(Result));
  }

  /* Test 7:*/
  template<>
  template<>
  void Grid2DQuadIntegration_object::test<7>()
  {

    set_test_name("Use mesh to integrate over quadrilateral domain, high-order too");

    Grid2D_Quad_MultiBlock_HO Mesh;
    
    int i,j;
    double Result(0), IntVal;

    int Blocks_Idir = 1;
    int Blocks_Jdir = 1;
    double Box_Width = 10.0;
    int IdirCells = 10;
    int JdirCells = 10;
    int Nghost = 2;
    int RecOrder = 1;
    int Stretching_Flag = 0;
    int Stretching_Type_Idir = 1;
    int Stretching_Type_Jdir = 1;
    double Stretching_Factor_Idir = 1.0;
    double Stretching_Factor_Jdir = 1.0;
    Vector2D SW = Vector2D(-2.0,-2.0);
    Vector2D SE = Vector2D(5.0, -0.5);
    Vector2D NE = Vector2D(3.5, 2.7);
    Vector2D NW = Vector2D(-1.0, 2.5);

    // Require high-order boundaries.
    Grid2D_Quad_Block_HO::setHighOrderBoundaryRepresentation();

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    // Set the square grid
    Mesh.Grid_Deformed_Box(Blocks_Idir, Blocks_Jdir,
			   SW, SE, NE, NW,
			   Stretching_Flag,
			   Stretching_Type_Idir,
			   Stretching_Type_Jdir,
			   Stretching_Factor_Idir,
			   Stretching_Factor_Jdir,
			   IdirCells, JdirCells, Nghost,
			   RecOrder);

    // == check integration
    for (i = 1; i<=Mesh(0,0).ICu + Mesh(0,0).Nghost - 1; ++i){
      for (j = 1; j<=Mesh(0,0).JCu + Mesh(0,0).Nghost - 1; ++j){

	// Calculate integration with low-order boundaries
	Result = Mesh(0,0).Integration.IntegrateFunctionOverCell(i,j,
								 Test_Example1,
								 14,Result);

	
	IntVal = Mesh(0,0).Integration.IntegrateFunctionOverCell(i,j,
								 Test_Example1,
								 Test_Example1_XDependencyIntegrated,
								 14,Result);

	// == check numerical integration with high-order boundaries
	ensure_distance("IntegrateFunctionOverCell()", IntVal, Result, AcceptedError(Result, 1.0e-12));

      }	// endfor
    } // endfor

  }


  /* Test 8:*/
  template<>
  template<>
  void Grid2DQuadIntegration_object::test<8>()
  {

    set_test_name("Integrate over high-order quadrilateral domain with polygonal adaptive integration");

    Grid2D_Quad_MultiBlock_HO Mesh;
    
    int i,j;
    double Result(0), IntVal;

    int Blocks_Idir = 1;
    int Blocks_Jdir = 1;
    double Box_Width = 10.0;
    int IdirCells = 10;
    int JdirCells = 10;
    int Nghost = 3;
    int RecOrder = 1;
    int Stretching_Flag = 0;
    int Stretching_Type_Idir = 1;
    int Stretching_Type_Jdir = 1;
    double Stretching_Factor_Idir = 1.0;
    double Stretching_Factor_Jdir = 1.0;
    Vector2D SW = Vector2D(-2.0,-2.0);
    Vector2D SE = Vector2D(5.0, -0.5);
    Vector2D NE = Vector2D(4.5, 2.7);
    Vector2D NW = Vector2D(-1.0, 3.3);

    // Require high-order boundaries.
    Grid2D_Quad_Block_HO::setHighOrderBoundaryRepresentation();
    
    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    Grid2D_Quad_Block_HO::setPolygonalAdaptiveQuadratureIntegrationON();
    NumericalLibrary_Execution_Mode::Output_Error_Messages = false;

    // Set the square grid
    Mesh.Grid_Deformed_Box(Blocks_Idir, Blocks_Jdir,
			   SW, SE, NE, NW,
			   Stretching_Flag,
			   Stretching_Type_Idir,
			   Stretching_Type_Jdir,
			   Stretching_Factor_Idir,
			   Stretching_Factor_Jdir,
			   IdirCells, JdirCells, Nghost,
			   RecOrder);

    Mesh.Check_Multi_Block_Grid_Completely();

    // Activate the extension splines
    int INl, INu, JNl, JNu, NNi, NNj;
    INl = Mesh(0,0).INl;
    INu = Mesh(0,0).INu;
    JNl = Mesh(0,0).JNl;
    JNu = Mesh(0,0).JNu;
    NNi = Mesh(0,0).NNi;
    NNj = Mesh(0,0).NNj;
   
    Mesh(0,0).ExtendWest_BndNorthSpline.Create_Spline_Line(Mesh(0,0).Node[0][JNu].X,Mesh(0,0).Node[INl][JNu].X,2);
    Mesh(0,0).ExtendEast_BndNorthSpline.Create_Spline_Line(Mesh(0,0).Node[INu][JNu].X,Mesh(0,0).Node[NNi-1][JNu].X,2); 
    Mesh(0,0).ExtendWest_BndSouthSpline.Create_Spline_Line(Mesh(0,0).Node[0][JNl].X,Mesh(0,0).Node[INl][JNl].X,2); 
    Mesh(0,0).ExtendEast_BndSouthSpline.Create_Spline_Line(Mesh(0,0).Node[INu][JNl].X,Mesh(0,0).Node[NNi-1][JNl].X,2); 
    Mesh(0,0).ExtendNorth_BndEastSpline.Create_Spline_Line(Mesh(0,0).Node[INu][JNu].X, Mesh(0,0).Node[INu][NNj-1].X, 2); 
    Mesh(0,0).ExtendSouth_BndEastSpline.Create_Spline_Line(Mesh(0,0).Node[INu][0].X, Mesh(0,0).Node[INu][JNl].X, 2); 
    Mesh(0,0).ExtendNorth_BndWestSpline.Create_Spline_Line(Mesh(0,0).Node[INl][JNu].X, Mesh(0,0).Node[INl][NNj-1].X, 2); 
    Mesh(0,0).ExtendSouth_BndWestSpline.Create_Spline_Line(Mesh(0,0).Node[INl][0].X, Mesh(0,0).Node[INl][JNl].X, 2); 

    Mesh(0,0).ExtendWest_BndNorthSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendEast_BndNorthSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendWest_BndSouthSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendEast_BndSouthSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendNorth_BndEastSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendSouth_BndEastSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendNorth_BndWestSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendSouth_BndWestSpline.setBCtype(BC_DIRICHLET);

    Mesh(0,0).Schedule_Ghost_Cells_Update();
    Mesh(0,0).Update_Cells();

    // == check integration
    for (i = 0; i<Mesh(0,0).NCi; ++i){
      for (j = 0; j<Mesh(0,0).NCj; ++j){

	// Calculate integration with low-order boundaries
	Result = Mesh(0,0).Integration.IntegrateFunctionOverCell(i,j,
								 Test_Example1,
								 14,Result);

	IntVal = Mesh(0,0).Integration.IntegrateFunctionOverCell(i,j,
								 Test_Example1,
								 Test_Example1_XDependencyIntegrated,
								 14,Result);

	// == check numerical integration with high-order boundaries
	ostmClear();
	ostm() << "IntegrateFunctionOverCell(), Cell (" << i << "," << j << ")";
	ensure_distance(ostm().str(), IntVal, Result, AcceptedError(Result, 1.0e-4));

      }	// endfor
    } // endfor

  }

  /* Test 9:*/
  template<>
  template<>
  void Grid2DQuadIntegration_object::test<9>()
  {

    set_test_name("Integrate over high-order quadrilateral domain with Monte Carlo integration");

    Grid2D_Quad_MultiBlock_HO Mesh;
    
    int i,j;
    double Result(0), IntVal;

    int Blocks_Idir = 1;
    int Blocks_Jdir = 1;
    double Box_Width = 10.0;
    int IdirCells = 10;
    int JdirCells = 10;
    int Nghost = 3;
    int RecOrder = 1;
    int Stretching_Flag = 0;
    int Stretching_Type_Idir = 1;
    int Stretching_Type_Jdir = 1;
    double Stretching_Factor_Idir = 1.0;
    double Stretching_Factor_Jdir = 1.0;
    Vector2D SW = Vector2D(-2.0,-2.0);
    Vector2D SE = Vector2D(5.0, -0.5);
    Vector2D NE = Vector2D(4.5, 2.7);
    Vector2D NW = Vector2D(-1.0, 3.3);

    // Require high-order boundaries.
    Grid2D_Quad_Block_HO::setHighOrderBoundaryRepresentation();
    
    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    Grid2D_Quad_Block_HO::setPolygonalAdaptiveQuadratureIntegrationOFF();
    Grid2D_Quad_Block_HO::setMonteCarloIntegrationON();
    NumericalLibrary_Execution_Mode::Number_Monte_Carlo_Samples = 800000;

    // Set the square grid
    Mesh.Grid_Deformed_Box(Blocks_Idir, Blocks_Jdir,
			   SW, SE, NE, NW,
			   Stretching_Flag,
			   Stretching_Type_Idir,
			   Stretching_Type_Jdir,
			   Stretching_Factor_Idir,
			   Stretching_Factor_Jdir,
			   IdirCells, JdirCells, Nghost,
			   RecOrder);

    Mesh.Check_Multi_Block_Grid_Completely();

    // Activate the extension splines
    int INl, INu, JNl, JNu, NNi, NNj;
    INl = Mesh(0,0).INl;
    INu = Mesh(0,0).INu;
    JNl = Mesh(0,0).JNl;
    JNu = Mesh(0,0).JNu;
    NNi = Mesh(0,0).NNi;
    NNj = Mesh(0,0).NNj;
   
    Mesh(0,0).ExtendWest_BndNorthSpline.Create_Spline_Line(Mesh(0,0).Node[0][JNu].X,Mesh(0,0).Node[INl][JNu].X,2);
    Mesh(0,0).ExtendEast_BndNorthSpline.Create_Spline_Line(Mesh(0,0).Node[INu][JNu].X,Mesh(0,0).Node[NNi-1][JNu].X,2); 
    Mesh(0,0).ExtendWest_BndSouthSpline.Create_Spline_Line(Mesh(0,0).Node[0][JNl].X,Mesh(0,0).Node[INl][JNl].X,2); 
    Mesh(0,0).ExtendEast_BndSouthSpline.Create_Spline_Line(Mesh(0,0).Node[INu][JNl].X,Mesh(0,0).Node[NNi-1][JNl].X,2); 
    Mesh(0,0).ExtendNorth_BndEastSpline.Create_Spline_Line(Mesh(0,0).Node[INu][JNu].X, Mesh(0,0).Node[INu][NNj-1].X, 2); 
    Mesh(0,0).ExtendSouth_BndEastSpline.Create_Spline_Line(Mesh(0,0).Node[INu][0].X, Mesh(0,0).Node[INu][JNl].X, 2); 
    Mesh(0,0).ExtendNorth_BndWestSpline.Create_Spline_Line(Mesh(0,0).Node[INl][JNu].X, Mesh(0,0).Node[INl][NNj-1].X, 2); 
    Mesh(0,0).ExtendSouth_BndWestSpline.Create_Spline_Line(Mesh(0,0).Node[INl][0].X, Mesh(0,0).Node[INl][JNl].X, 2); 

    Mesh(0,0).ExtendWest_BndNorthSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendEast_BndNorthSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendWest_BndSouthSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendEast_BndSouthSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendNorth_BndEastSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendSouth_BndEastSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendNorth_BndWestSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendSouth_BndWestSpline.setBCtype(BC_DIRICHLET);

    Mesh(0,0).Schedule_Ghost_Cells_Update();
    Mesh(0,0).Update_Cells();

    // == check integration
    for (i = 0; i<Mesh(0,0).NCi; ++i){
      for (j = 0; j<Mesh(0,0).NCj; ++j){

	// Calculate integration with low-order boundaries
	Result = Mesh(0,0).Integration.IntegrateFunctionOverCell(i,j,
								 Test_Example1,
								 14,Result);

	IntVal = Mesh(0,0).Integration.IntegrateFunctionOverCell(i,j,
								 Test_Example1,
								 Test_Example1_XDependencyIntegrated,
								 14,Result);

	// == check numerical integration with high-order boundaries
	ostmClear();
	ostm() << "IntegrateFunctionOverCell(), Cell (" << i << "," << j << ")";
	ensure_distance(ostm().str(), IntVal, Result, AcceptedError(Result, 1.0e-2));

      }	// endfor
    } // endfor

  }


  /* Test 10:*/
  template<>
  template<>
  void Grid2DQuadIntegration_object::test<10>()
  {

    set_test_name("Integrate over high-order quadrilateral domain with Monte Carlo integration");

    Grid2D_Quad_MultiBlock_HO Mesh;
    
    int i,j;
    double Result(0), IntVal;

    int Blocks_Idir = 1;
    int Blocks_Jdir = 1;
    double Box_Width = 10.0;
    int IdirCells = 10;
    int JdirCells = 10;
    int Nghost = 3;
    int RecOrder = 1;
    int Stretching_Flag = 0;
    int Stretching_Type_Idir = 1;
    int Stretching_Type_Jdir = 1;
    double Stretching_Factor_Idir = 1.0;
    double Stretching_Factor_Jdir = 1.0;
    Vector2D SW = Vector2D(0.0,0.0);
    Vector2D SE = Vector2D(5.0, 0.0);
    Vector2D NE = Vector2D(5.0, 5.0);
    Vector2D NW = Vector2D(0.0, 5.0);

    // Require high-order boundaries.
    Grid2D_Quad_Block_HO::setHighOrderBoundaryRepresentation();

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    Grid2D_Quad_Block_HO::setPolygonalAdaptiveQuadratureIntegrationON();
    NumericalLibrary_Execution_Mode::Output_Error_Messages = false;

    // Set the square grid
    Mesh.Grid_Deformed_Box(Blocks_Idir, Blocks_Jdir,
			   SW, SE, NE, NW,
			   Stretching_Flag,
			   Stretching_Type_Idir,
			   Stretching_Type_Jdir,
			   Stretching_Factor_Idir,
			   Stretching_Factor_Jdir,
			   IdirCells, JdirCells, Nghost,
			   RecOrder);

    // Activate the extension splines
    int INl, INu, JNl, JNu, NNi, NNj;
    INl = Mesh(0,0).INl;
    INu = Mesh(0,0).INu;
    JNl = Mesh(0,0).JNl;
    JNu = Mesh(0,0).JNu;
    NNi = Mesh(0,0).NNi;
    NNj = Mesh(0,0).NNj;
   
    Mesh(0,0).ExtendWest_BndNorthSpline.Create_Spline_Line(Mesh(0,0).Node[0][JNu].X,Mesh(0,0).Node[INl][JNu].X,2);
    Mesh(0,0).ExtendEast_BndNorthSpline.Create_Spline_Line(Mesh(0,0).Node[INu][JNu].X,Mesh(0,0).Node[NNi-1][JNu].X,2); 
    Mesh(0,0).ExtendWest_BndSouthSpline.Create_Spline_Line(Mesh(0,0).Node[0][JNl].X,Mesh(0,0).Node[INl][JNl].X,2); 
    Mesh(0,0).ExtendEast_BndSouthSpline.Create_Spline_Line(Mesh(0,0).Node[INu][JNl].X,Mesh(0,0).Node[NNi-1][JNl].X,2); 
    Mesh(0,0).ExtendNorth_BndEastSpline.Create_Spline_Line(Mesh(0,0).Node[INu][JNu].X, Mesh(0,0).Node[INu][NNj-1].X, 2); 
    Mesh(0,0).ExtendSouth_BndEastSpline.Create_Spline_Line(Mesh(0,0).Node[INu][0].X, Mesh(0,0).Node[INu][JNl].X, 2); 
    Mesh(0,0).ExtendNorth_BndWestSpline.Create_Spline_Line(Mesh(0,0).Node[INl][JNu].X, Mesh(0,0).Node[INl][NNj-1].X, 2); 
    Mesh(0,0).ExtendSouth_BndWestSpline.Create_Spline_Line(Mesh(0,0).Node[INl][0].X, Mesh(0,0).Node[INl][JNl].X, 2); 

    Mesh(0,0).ExtendWest_BndNorthSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendEast_BndNorthSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendWest_BndSouthSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendEast_BndSouthSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendNorth_BndEastSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendSouth_BndEastSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendNorth_BndWestSpline.setBCtype(BC_DIRICHLET);
    Mesh(0,0).ExtendSouth_BndWestSpline.setBCtype(BC_DIRICHLET);

    Mesh(0,0).Schedule_Ghost_Cells_Update();
    Mesh(0,0).Update_Cells();

    // == check integration
    for (i = 0; i<Mesh(0,0).NCi; ++i){
      for (j = 0; j<Mesh(0,0).NCj; ++j){

	// Calculate integration with low-order boundaries
	Result = Mesh(0,0).Integration.IntegrateFunctionOverCell(i,j,
								 Test_Example1,
								 14,Result);

	
	IntVal = Mesh(0,0).Integration.IntegrateFunctionOverCell(i,j,
								 Test_Example1,
								 Test_Example1_XDependencyIntegrated,
								 14,Result);

	// == check numerical integration with high-order boundaries
	ostmClear();
	ostm() << "IntegrateFunctionOverCell(), Cell (" << i << "," << j << ")";
	ensure_distance(ostm().str(), IntVal, Result, AcceptedError(Result, 1.0e-12));

      }	// endfor
    } // endfor

  }


  /* Test 11:*/
  template<>
  template<>
  void Grid2DQuadIntegration_object::test<11>()
  {

    set_test_name("Test QuadrilateralQuadrature routine for integration on rotated rectangle");

    Node2D_HO SW, NW, NE, SE;
    double NumericResult(0), AnalyticResult;

    // Set the integration domain (i.e. a rectangle rotated counterclockwise at 30 degrees) 
    double cos30, sin30, alpha1, alpha2;
    cos30 = cos(PI/6);
    sin30 = sin(PI/6);
    alpha1 = arctan(2.0,1.0);
    alpha2 = arctan(2.0,2.0);

    SW.x() = cos30;
    SW.y() = sin30;
    
    SE.x() = 2*cos30;
    SE.y() = 2*sin30;

    NE.x() = 2*sqrt(2.0)*cos(alpha2 + PI/6);
    NE.y() = 2*sqrt(2.0)*sin(alpha2 + PI/6);

    NW.x() = sqrt(5.0)*cos(alpha1 + PI/6);
    NW.y() = sqrt(5.0)*sin(alpha1 + PI/6);
    
    AnalyticResult = 14./3.*(1. + cos30*sin30) + 8./3.*(1. - cos30*sin30) + 3.*(cos30*cos30 - sin30*sin30);

    NumericResult = QuadrilateralQuadrature(Test_Example1,
					    SW, NW, NE, SE,
					    14, NumericResult);
   
    ensure_distance("Numeric result", NumericResult, AnalyticResult, AcceptedError(AnalyticResult));
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

