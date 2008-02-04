/*!\file test_HO_MultiBlock_Grid2DQuad.cc
  \brief Regression tests for 2D high-order quadrilateral multi-block grid. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../HO_Grid2DQuadMultiBlock.h"
#include "HO_Grid2DQuadMultiBlock_InputForTesting.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_Grid2DQuadMultiBlock_HO : public TestData {

    // Local variables
  public:
    Grid2D_Quad_MultiBlock_HO   MeshBlk, MeshBlk_Fine; // the mesh
    Grid2D_Quad_Block_HO   Grid;

    Grid2DTesting_Input_Parameters IP;

    int error_flag;
    Vector2D GQPoint1, GQPoint2; // Gauss Quadrature Points
    int iCell,jCell;		// cell coordinates

    Spline2DInterval_HO * SInfoNULL;
    
    // Member functions

    // Constructor
    Data_Grid2DQuadMultiBlock_HO(void);

    // Destructor
    ~Data_Grid2DQuadMultiBlock_HO(void);

    // Create Mesh
    template<class Input_Parameters>
    void CreateMesh(Grid2D_Quad_MultiBlock_HO & _MeshBlk_, Input_Parameters & IP) throw(std::runtime_error);

    // Set flux calculation method in all the mesh blocks
    template<class Input_Parameters>
    void SetFluxCalculationMethod(Grid2D_Quad_MultiBlock_HO & _MeshBlk_,
				  int FluxMethod, Input_Parameters & IP);

    // Verify the number of constrained Gauss quadrature points per edge for all the interior mesh cells
    // Ghost cells are not used for constrained reconstruction and therefore they are not checked.
    template<class Input_Parameters>
    void CheckNumberOfConstrainedGQP(Grid2D_Quad_MultiBlock_HO & MeshBlk,
				     Input_Parameters & IP,
				     int Result);

  private:
    
  };

  // Constructor
  Data_Grid2DQuadMultiBlock_HO::Data_Grid2DQuadMultiBlock_HO(void){
    /* Set the global path for this test suite. 
       It's a good practice to put it in the constructor of the data class in order to be set
       automatically for each individual test. Declare it relative to the /src_2D directory,
       otherwise the framework might not find the input and output files. */
    
    set_test_suite_path("Grid/UnitTests/");

    SInfoNULL = NULL;

    IP.Set_Default_Input_Parameters();
    set_local_output_path("HO_MultiBlockQuadGrids");
    set_local_input_path("HO_MultiBlockQuadGrids");
  }

  Data_Grid2DQuadMultiBlock_HO::~Data_Grid2DQuadMultiBlock_HO(void){
    // reset to default value
    Grid2D_Quad_Block_HO::setDefaultBoundaryRepresentation();
  }

  template<class Input_Parameters>
  void Data_Grid2DQuadMultiBlock_HO::CreateMesh(Grid2D_Quad_MultiBlock_HO & _MeshBlk_,
						Input_Parameters & IP) throw(std::runtime_error){

    /* Initialize all static variables within the class */
    if (IP.IncludeHighOrderBoundariesRepresentation == ON){
      Grid2D_Quad_Block_HO::setHighOrderBoundaryRepresentation();
    } else {
      Grid2D_Quad_Block_HO::setLowOrderBoundaryRepresentation();
    }

    error_flag = _MeshBlk_.Multi_Block_Grid(IP);
    
    if (error_flag) {
      throw runtime_error("CreateMesh() ERROR: Unable to create valid Euler2D multi-block mesh.");
    }
   
  }

  template<class Input_Parameters>
  void Data_Grid2DQuadMultiBlock_HO::SetFluxCalculationMethod(Grid2D_Quad_MultiBlock_HO & MeshBlk,
							      int FluxMethod, Input_Parameters & IP){

    int iBlock, jBlock;

    for (iBlock = 0; iBlock < IP.Number_of_Blocks_Idir ; ++iBlock){
      for (jBlock = 0; jBlock < IP.Number_of_Blocks_Jdir ; ++jBlock){

	if(MeshBlk(iBlock,jBlock).BndNorthSpline.bc[0] != BC_NONE){
	  MeshBlk(iBlock,jBlock).BndNorthSpline.setFluxCalcMethod(FluxMethod);
	}
	if(MeshBlk(iBlock,jBlock).BndSouthSpline.bc[0] != BC_NONE){
	  MeshBlk(iBlock,jBlock).BndSouthSpline.setFluxCalcMethod(FluxMethod);
	}
	if(MeshBlk(iBlock,jBlock).BndEastSpline.bc[0] != BC_NONE){
	  MeshBlk(iBlock,jBlock).BndEastSpline.setFluxCalcMethod(FluxMethod);
	}
	if(MeshBlk(iBlock,jBlock).BndWestSpline.bc[0] != BC_NONE){
	  MeshBlk(iBlock,jBlock).BndWestSpline.setFluxCalcMethod(FluxMethod);
	}
      }	// endfor
    } // endfor

  }

  template<class Input_Parameters>
  void Data_Grid2DQuadMultiBlock_HO::CheckNumberOfConstrainedGQP(Grid2D_Quad_MultiBlock_HO & MeshBlk,
								 Input_Parameters & IP,
								 int Result){

    int iBlock, jBlock;


    for (iBlock = 0; iBlock < IP.Number_of_Blocks_Idir ; ++iBlock){
      for (jBlock = 0; jBlock < IP.Number_of_Blocks_Jdir ; ++jBlock){

	for (iCell=MeshBlk(iBlock,jBlock).ICl; iCell<=MeshBlk(iBlock,jBlock).ICu; ++iCell){
	  for (jCell=MeshBlk(iBlock,jBlock).JCl; jCell<=MeshBlk(iBlock,jBlock).JCu; ++jCell){


#if 0
	    if (jCell == MeshBlk(iBlock,jBlock).JCu){
	      // check cells with boundary at North

	      ostmClear();
	      ostm() << "NorthBnd Check, Block [" << iBlock << "," << jBlock << "], cell[" << iCell << "," << jCell << "]";
	    
	      if (MeshBlk(iBlock,jBlock).BndNorthSpline.getFluxCalcMethod() == ReconstructionBasedFlux){
		ensure_equals(ostm().str(), MeshBlk(iBlock,jBlock).NumOfConstrainedGaussQuadPoints_North(iCell,jCell) ,
			      Result);
	      } else {
		// zero constrained GQP
		ensure_equals(ostm().str(), MeshBlk(iBlock,jBlock).NumOfConstrainedGaussQuadPoints_North(iCell,jCell) ,
			      0);
	      }
	    } 
	    
	    if (jCell == MeshBlk(iBlock,jBlock).JCl){
	      // check cells with boundary at South

	      ostmClear();
	      ostm() << "SouthBnd Check, Block [" << iBlock << "," << jBlock << "], cell[" << iCell << "," << jCell << "]";
	      
	      if (MeshBlk(iBlock,jBlock).BndSouthSpline.getFluxCalcMethod() == ReconstructionBasedFlux){
		ensure_equals(ostm().str(), MeshBlk(iBlock,jBlock).NumOfConstrainedGaussQuadPoints_South(iCell,jCell) ,
			      Result);
	      } else {
		// zero constrained GQP
		ensure_equals(ostm().str(), MeshBlk(iBlock,jBlock).NumOfConstrainedGaussQuadPoints_South(iCell,jCell) ,
			      0);
	      }
	    } 

	    if (iCell == MeshBlk(iBlock,jBlock).ICu){
	      // check cells with boundary at East

	      ostmClear();
	      ostm() << "EastBnd Check, Block [" << iBlock << "," << jBlock << "], cell[" << iCell << "," << jCell << "]";

	      if (MeshBlk(iBlock,jBlock).BndEastSpline.getFluxCalcMethod() == ReconstructionBasedFlux){
		ensure_equals(ostm().str(), MeshBlk(iBlock,jBlock).NumOfConstrainedGaussQuadPoints_East(iCell,jCell) ,
			      Result);
	      } else {
		// zero constrained GQP
		ensure_equals(ostm().str(), MeshBlk(iBlock,jBlock).NumOfConstrainedGaussQuadPoints_East(iCell,jCell) ,
			      0);
	      }
	    } 
	    
	    if (iCell == MeshBlk(iBlock,jBlock).ICl){
	      // check cells with boundary at West
	      
	      ostmClear();
	      ostm() << "WestBnd Check, Block [" << iBlock << "," << jBlock << "], cell[" << iCell << "," << jCell << "]";
	      
	      if (MeshBlk(iBlock,jBlock).BndWestSpline.getFluxCalcMethod() == ReconstructionBasedFlux){
		ensure_equals(ostm().str(), MeshBlk(iBlock,jBlock).NumOfConstrainedGaussQuadPoints_West(iCell,jCell) ,
			      Result);
	      } else {
		// zero constrained GQP
		ensure_equals(ostm().str(), MeshBlk(iBlock,jBlock).NumOfConstrainedGaussQuadPoints_West(iCell,jCell) ,
			      0);
	      }
	    }

	    if (iCell != MeshBlk(iBlock,jBlock).ICl && iCell != MeshBlk(iBlock,jBlock).ICu &&
		jCell != MeshBlk(iBlock,jBlock).JCl && jCell != MeshBlk(iBlock,jBlock).JCu){
	      // check interior cells 

	      ostmClear();
	      ostm() << "InteriorCell North, Block [" << iBlock << "," << jBlock << "], cell[" << iCell << "," << jCell << "]";
	      ensure_equals(ostm().str(), MeshBlk(iBlock,jBlock).NumOfConstrainedGaussQuadPoints_North(iCell,jCell), 0);

	      ostmClear();
	      ostm() << "InteriorCell South, Block [" << iBlock << "," << jBlock << "], cell[" << iCell << "," << jCell << "]";
	      ensure_equals(ostm().str(), MeshBlk(iBlock,jBlock).NumOfConstrainedGaussQuadPoints_South(iCell,jCell), 0);

	      ostmClear();
	      ostm() << "InteriorCell East, Block [" << iBlock << "," << jBlock << "], cell[" << iCell << "," << jCell << "]";
	      ensure_equals(ostm().str(), MeshBlk(iBlock,jBlock).NumOfConstrainedGaussQuadPoints_East(iCell,jCell), 0);

	      ostmClear();
	      ostm() << "InteriorCell West, Block [" << iBlock << "," << jBlock << "], cell[" << iCell << "," << jCell << "]";
	      ensure_equals(ostm().str(), MeshBlk(iBlock,jBlock).NumOfConstrainedGaussQuadPoints_West(iCell,jCell), 0);
	    }
#endif

	  } // endfor
	}// endfor

      }// endfor
    }// endfor
  }

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

  /* Test 6:*/
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<6>()
  {

    set_test_name("Grid rectangular box");

    RunRegression = ON;
    
    // Create block grid
    Grid2D_Quad_MultiBlock_HO MultiBlockGrid;

    // Input parameters
    int Number_of_Blocks_Idir = 5;
    int Number_of_Blocks_Jdir = 3;
    double Width = 2.0;
    double Height = 3.0;
    int Number_of_Cells_Idir = 40;
    int Number_of_Cells_Jdir = 30;
    int Number_of_Ghost_Cells = 2;

    // Build grid
    MultiBlockGrid.Grid_Rectangular_Box(Number_of_Blocks_Idir,
                                        Number_of_Blocks_Jdir,
                                        Width,
                                        Height,
                                        Number_of_Cells_Idir,
					Number_of_Cells_Jdir,
					Number_of_Ghost_Cells,
					0);

    MultiBlockGrid.Translate_Multi_Block_Grid(Vector2D(1.0,2.0));
    MultiBlockGrid.Scale_Multi_Block_Grid(2.0);
    MultiBlockGrid.Reflect_Multi_Block_Grid();
    MultiBlockGrid.Rotate_Multi_Block_Grid(0.785398163);
  
    if (RunRegression){

      //open file for output
      CurrentFile = "Current_RectangularBox.dat";
      MasterFile  = "RectangularBox.dat";
      Open_Output_File(CurrentFile);
      MultiBlockGrid.Output_Tecplot(out());
      RunRegressionTest("Rectangular Box", CurrentFile, MasterFile, 5.0e-12, 5.0e-12);
      
      //open file for all mesh node output
      CurrentFile  = "Current_RectangularBox_node.dat";
      MasterFile  = "RectangularBox_node.dat";
      Open_Output_File(CurrentFile);
      MultiBlockGrid.Output_Nodes_Tecplot(out());
      RunRegressionTest("Rectangular Box node", CurrentFile, MasterFile, 5.0e-12, 5.0e-12);

      //open file for all mesh cell output
      CurrentFile  = "Current_RectangularBox_cell.dat";
      MasterFile  = "RectangularBox_cell.dat";
      Open_Output_File(CurrentFile);
      MultiBlockGrid.Output_Cells_Tecplot(out());
      RunRegressionTest("Rectangular Box cell", CurrentFile, MasterFile, 5.0e-12, 5.0e-12);

    } else {
      //open file for interior node output
      MasterFile  = "RectangularBox.dat";
      Open_Output_File(MasterFile);
      MultiBlockGrid.Output_Tecplot(out());

      //open file for all mesh node output
      MasterFile  = "RectangularBox_node.dat";
      Open_Output_File(MasterFile);
      MultiBlockGrid.Output_Nodes_Tecplot(out());

      //open file for all mesh cell output
      MasterFile  = "RectangularBox_cell.dat";
      Open_Output_File(MasterFile);
      MultiBlockGrid.Output_Cells_Tecplot(out());
    }

  }

  /* Test 7:*/
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<7>()
  {

    set_test_name("Grid flat plate");

    RunRegression = ON;
    
    // Create block grid
    Grid2D_Quad_MultiBlock_HO MultiBlockGrid;

    // Input parameters
    int Number_of_Blocks_Idir = 5;
    int Number_of_Blocks_Jdir = 3;
    double Length = 2.5;
    int Flat_Plate_Type = 1;
    int Stretching_Flag = 1;
    double Stretching_Factor_Idir = 2.0;
    double Stretching_Factor_Jdir = 3.0;
    int Number_of_Cells_Idir = 40;
    int Number_of_Cells_Jdir = 30;
    int Number_of_Ghost_Cells = 2;

    // Build grid
    MultiBlockGrid.Grid_Flat_Plate(Number_of_Blocks_Idir,
				   Number_of_Blocks_Jdir,
				   Length,
				   Flat_Plate_Type,
				   Stretching_Flag,
				   Stretching_Factor_Idir,
				   Stretching_Factor_Jdir,
				   Number_of_Cells_Idir,
				   Number_of_Cells_Jdir,
				   Number_of_Ghost_Cells,
				   0);

    MultiBlockGrid.Translate_Multi_Block_Grid(Vector2D(1.0,2.0));
    MultiBlockGrid.Scale_Multi_Block_Grid(2.0);
    MultiBlockGrid.Reflect_Multi_Block_Grid();
    MultiBlockGrid.Rotate_Multi_Block_Grid(0.785398163);
  
    if (RunRegression){

      //open file for output
      CurrentFile = "Current_FlatPlate.dat";
      MasterFile  = "FlatPlate.dat";
      Open_Output_File(CurrentFile);
      MultiBlockGrid.Output_Tecplot(out());
      RunRegressionTest("Flat Plate", CurrentFile, MasterFile, 5.0e-12, 5.0e-12);
      
      //open file for all mesh node output
      CurrentFile  = "Current_FlatPlate_node.dat";
      MasterFile  = "FlatPlate_node.dat";
      Open_Output_File(CurrentFile);
      MultiBlockGrid.Output_Nodes_Tecplot(out());
      RunRegressionTest("Flat Plate node", CurrentFile, MasterFile, 5.0e-12, 5.0e-12);

      //open file for all mesh cell output
      CurrentFile  = "Current_FlatPlate_cell.dat";
      MasterFile  = "FlatPlate_cell.dat";
      Open_Output_File(CurrentFile);
      MultiBlockGrid.Output_Cells_Tecplot(out());
      RunRegressionTest("Flat Plate cell", CurrentFile, MasterFile, 5.0e-12, 5.0e-12);

    } else {
      //open file for interior node output
      MasterFile  = "FlatPlate.dat";
      Open_Output_File(MasterFile);
      MultiBlockGrid.Output_Tecplot(out());

      //open file for all mesh node output
      MasterFile  = "FlatPlate_node.dat";
      Open_Output_File(MasterFile);
      MultiBlockGrid.Output_Nodes_Tecplot(out());

      //open file for all mesh cell output
      MasterFile  = "FlatPlate_cell.dat";
      Open_Output_File(MasterFile);
      MultiBlockGrid.Output_Cells_Tecplot(out());
    }

  }

  /* Test 8:*/
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<8>()
  {

    set_test_name("Grid 2D Laminar Flame");

    RunRegression = ON;
    
    // Create block grid
    Grid2D_Quad_MultiBlock_HO MultiBlockGrid;

    // Input parameters
    int Number_of_Blocks_Idir = 12;
    int Number_of_Blocks_Jdir = 3;
    double Height = 2.0;
    double Length = 3.0;
    int Number_of_Cells_Idir = 40;
    int Number_of_Cells_Jdir = 30;
    int Number_of_Ghost_Cells = 2;
    int Flame_Type_Flag = IC_RESTART;

    // Build grid
    MultiBlockGrid.Grid_2D_Laminar_Flame(Number_of_Blocks_Idir,
					 Number_of_Blocks_Jdir,
					 Length,
					 Height,
					 Number_of_Cells_Idir,
					 Number_of_Cells_Jdir,
					 Number_of_Ghost_Cells,
					 0,
					 Flame_Type_Flag);

    MultiBlockGrid.Translate_Multi_Block_Grid(Vector2D(1.0,2.0));
    MultiBlockGrid.Scale_Multi_Block_Grid(2.0);
    MultiBlockGrid.Reflect_Multi_Block_Grid();
    MultiBlockGrid.Rotate_Multi_Block_Grid(0.785398163);
  
    if (RunRegression){

      //open file for output
      CurrentFile = "Current_2D_Laminar_Flame.dat";
      MasterFile  = "2D_Laminar_Flame.dat";
      Open_Output_File(CurrentFile);
      MultiBlockGrid.Output_Tecplot(out());
      RunRegressionTest("2D Laminar Flame", CurrentFile, MasterFile, 5.0e-12, 5.0e-12);
      
      //open file for all mesh node output
      CurrentFile  = "Current_2D_Laminar_Flame_node.dat";
      MasterFile  = "2D_Laminar_Flame_node.dat";
      Open_Output_File(CurrentFile);
      MultiBlockGrid.Output_Nodes_Tecplot(out());
      RunRegressionTest("2D Laminar Flame node", CurrentFile, MasterFile, 5.0e-12, 5.0e-12);

      //open file for all mesh cell output
      CurrentFile  = "Current_2D_Laminar_Flame_cell.dat";
      MasterFile  = "2D_Laminar_Flame_cell.dat";
      Open_Output_File(CurrentFile);
      MultiBlockGrid.Output_Cells_Tecplot(out());
      RunRegressionTest("2D Laminar Flame cell", CurrentFile, MasterFile, 5.0e-12, 5.0e-12);

    } else {
      //open file for interior node output
      MasterFile  = "2D_Laminar_Flame.dat";
      Open_Output_File(MasterFile);
      MultiBlockGrid.Output_Tecplot(out());

      //open file for all mesh node output
      MasterFile  = "2D_Laminar_Flame_node.dat";
      Open_Output_File(MasterFile);
      MultiBlockGrid.Output_Nodes_Tecplot(out());

      //open file for all mesh cell output
      MasterFile  = "2D_Laminar_Flame_cell.dat";
      Open_Output_File(MasterFile);
      MultiBlockGrid.Output_Cells_Tecplot(out());
    }

  }

  /* Test 9:*/
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<9>()
  {

    set_test_name("Grid NACA AirFoil");

    RunRegression = ON;
    
    // Create block grid
    Grid2D_Quad_MultiBlock_HO MultiBlockGrid;

    // Input parameters
    int Number_of_Blocks_Idir = 12;
    int Number_of_Blocks_Jdir = 3;
    char NACA_Aerofoil_Type_ptr[6] = "23015";
    double Chord_Length = 2.0;
    int Number_of_Cells_Idir = 40;
    int Number_of_Cells_Jdir = 30;
    int Number_of_Ghost_Cells = 5;
    
    // Build grid
    MultiBlockGrid.Grid_NACA_Aerofoil(Number_of_Blocks_Idir,
				      Number_of_Blocks_Jdir,
				      NACA_Aerofoil_Type_ptr,
				      Chord_Length,
				      Number_of_Cells_Idir,
				      Number_of_Cells_Jdir,
				      Number_of_Ghost_Cells,
				      0);

    MultiBlockGrid.Translate_Multi_Block_Grid(Vector2D(1.0,2.0));
    MultiBlockGrid.Scale_Multi_Block_Grid(2.0);
    MultiBlockGrid.Reflect_Multi_Block_Grid();
  
    if (RunRegression){

      //open file for output
      CurrentFile = "Current_NACA_Aerofoil.dat";
      MasterFile  = "NACA_Aerofoil.dat";
      Open_Output_File(CurrentFile);
      MultiBlockGrid.Output_Tecplot(out());
      RunRegressionTest("NACA Aerofoil", CurrentFile, MasterFile, 5.0e-12, 5.0e-12);
      
      //open file for all mesh node output
      CurrentFile  = "Current_NACA_Aerofoil_node.dat";
      MasterFile  = "NACA_Aerofoil_node.dat";
      Open_Output_File(CurrentFile);
      MultiBlockGrid.Output_Nodes_Tecplot(out());
      RunRegressionTest("NACA Aerofoil node", CurrentFile, MasterFile, 5.0e-12, 5.0e-12);

      //open file for all mesh cell output
      CurrentFile  = "Current_NACA_Aerofoil_cell.dat";
      MasterFile  = "NACA_Aerofoil_cell.dat";
      Open_Output_File(CurrentFile);
      MultiBlockGrid.Output_Cells_Tecplot(out());
      RunRegressionTest("NACA Aerofoil cell", CurrentFile, MasterFile, 5.0e-12, 5.0e-12);

    } else {
      //open file for interior node output
      MasterFile  = "NACA_Aerofoil.dat";
      Open_Output_File(MasterFile);
      MultiBlockGrid.Output_Tecplot(out());

      //open file for all mesh node output
      MasterFile  = "NACA_Aerofoil_node.dat";
      Open_Output_File(MasterFile);
      MultiBlockGrid.Output_Nodes_Tecplot(out());

      //open file for all mesh cell output
      MasterFile  = "NACA_Aerofoil_cell.dat";
      Open_Output_File(MasterFile);
      MultiBlockGrid.Output_Cells_Tecplot(out());
    }

  }

  /* Test 10:*/
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<10>()
  {

    set_test_name("Grid Nozzle");

    RunRegression = ON;
    
    // Create block grid
    Grid2D_Quad_MultiBlock_HO MultiBlockGrid;

    // Input parameters
    int Number_of_Blocks_Idir = 12;
    int Number_of_Blocks_Jdir = 3;
    double Length_Nozzle = 10.0;
    double Radius_Chamber = 2;
    double Radius_Nozzle_Exit = 4;
    double Radius_Nozzle_Throat = 1;
    int Nozzle_Type = 1;
    int Stretching_Flag = 1;
    int Stretching_Type_Idir = 1;
    int Stretching_Type_Jdir = 1;
    double Stretching_Factor_Idir = 1.9;
    double Stretching_Factor_Jdir = 2.1;
    int Number_of_Cells_Idir = 40;
    int Number_of_Cells_Jdir = 30;
    int Number_of_Ghost_Cells = 5;
    
    // Build grid
    MultiBlockGrid.Grid_Nozzle(Number_of_Blocks_Idir,
			       Number_of_Blocks_Jdir,
			       Length_Nozzle,
			       Radius_Chamber,
			       Radius_Nozzle_Exit,
			       Radius_Nozzle_Throat,
			       Nozzle_Type,
			       Stretching_Flag,
			       Stretching_Type_Idir,
			       Stretching_Type_Jdir,
			       Stretching_Factor_Idir,
			       Stretching_Factor_Jdir,
			       Number_of_Cells_Idir,
			       Number_of_Cells_Jdir,
			       Number_of_Ghost_Cells,
			       0);

    MultiBlockGrid.Translate_Multi_Block_Grid(Vector2D(1.0,2.0));
    MultiBlockGrid.Scale_Multi_Block_Grid(2.0);
    MultiBlockGrid.Reflect_Multi_Block_Grid();
  
    if (RunRegression){

      //open file for output
      CurrentFile = "Current_Nozzle.dat";
      MasterFile  = "Nozzle.dat";
      Open_Output_File(CurrentFile);
      MultiBlockGrid.Output_Tecplot(out());
      RunRegressionTest("Nozzle", CurrentFile, MasterFile, 5.0e-12, 5.0e-12);
      
      //open file for all mesh node output
      CurrentFile  = "Current_Nozzle_node.dat";
      MasterFile  = "Nozzle_node.dat";
      Open_Output_File(CurrentFile);
      MultiBlockGrid.Output_Nodes_Tecplot(out());
      RunRegressionTest("Nozzle node", CurrentFile, MasterFile, 5.0e-12, 5.0e-12);

      //open file for all mesh cell output
      CurrentFile  = "Current_Nozzle_cell.dat";
      MasterFile  = "Nozzle_cell.dat";
      Open_Output_File(CurrentFile);
      MultiBlockGrid.Output_Cells_Tecplot(out());
      RunRegressionTest("Nozzle cell", CurrentFile, MasterFile, 5.0e-12, 5.0e-12);

    } else {
      //open file for interior node output
      MasterFile  = "Nozzle.dat";
      Open_Output_File(MasterFile);
      MultiBlockGrid.Output_Tecplot(out());

      //open file for all mesh node output
      MasterFile  = "Nozzle_node.dat";
      Open_Output_File(MasterFile);
      MultiBlockGrid.Output_Nodes_Tecplot(out());

      //open file for all mesh cell output
      MasterFile  = "Nozzle_cell.dat";
      Open_Output_File(MasterFile);
      MultiBlockGrid.Output_Cells_Tecplot(out());
    }

  }

  /* Test 11:*/
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<11>()
  {
    set_test_name("Grid_Rectangular_Box()");
    RunRegression = ON;
 
    // Add test particular input parameters
    IP.i_Grid = GRID_RECTANGULAR_BOX;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 10;
    IP.Number_of_Cells_Jdir = 10;
    IP.Number_of_Ghost_Cells = 2;
    IP.Space_Accuracy = 4;
    IP.Box_Width = 2;
    IP.Box_Height = 2;
    IP.X_Shift.x = 6;
    IP.X_Shift.y = 6;
    IP.X_Scale = 2.0;
    IP.X_Rotate = 10.0;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_Smooth_Quad_Block = ON;
    strcpy(IP.BC_North_Type, "Reflection");
    strcpy(IP.BC_East_Type, "Reflection");
    strcpy(IP.BC_South_Type, "Reflection");
    strcpy(IP.BC_West_Type, "Reflection");
    IP.BCs_Specified = ON;
    IP.BC_North = BC_REFLECTION;
    IP.BC_South = BC_REFLECTION;
    IP.BC_East = BC_REFLECTION;
    IP.BC_West = BC_REFLECTION;

    // Build the mesh
    CreateMesh(MeshBlk,IP);

    MasterFile = "GridRectangularBox_TecplotCells.dat";
    CurrentFile = "Current_GridRectangularBox_TecplotCells.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // OutputMeshTecplot
      MeshBlk.Output_Cells_Tecplot(out());

      // check
      RunRegressionTest("Cells_Tecplot", CurrentFile, MasterFile, 1.0e-9, 1.0e-9);

    } else {
      Open_Output_File(MasterFile);
      
      // write data
      MeshBlk.Output_Cells_Tecplot(out());
    }

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

