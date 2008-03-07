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

  double Function_XCentroid(const double &x, const double &y){
    return x;
  }

  double Function_YCentroid(const double &x, const double &y){
    return y;
  }

  double Function_area(const double &x, const double &y){
    return 1;
  }

  /******************************************************************************************
   * Function for generating the geom coeff. for a cartesian cell
   ******************************************************************************************/
  
  double GeomCoeffCartesian(int p1, int p2, double deltaX, double deltaY, double deltaXC, double deltaYC){

    /* p1 -> the first power coefficient
       p2 -> the second power coefficient
       deltaX -> the grid size in the X direction
       deltaY -> the grid size in the Y direction
       deltaXC -> the X distance between the center of the reconstructed cell and that of the cell used in the reconstruction
       deltaYC -> the Y distance between the center of the reconstructed cell and that of the cell used in the reconstruction

       Obs. To compute the coefficient of the reconstructed cell, deltaXC and deltaYC must be ZERO. 
    */

    double val1, val2;
    double coef_x1, coef_x2, coef_y1, coef_y2;

    val1 = val2 = 0.0;
    coef_x1 = deltaX/2  + deltaXC;
    coef_x2 = -deltaX/2 + deltaXC;
    coef_y1 = deltaY/2  + deltaYC;
    coef_y2 = -deltaY/2 + deltaYC;

    for (int m=1; m<=p1+1; m++){
      val1 += pow(coef_x1,p1+1-m)*pow(coef_x2,m-1);
    }
    for (int l=1; l<=p2+1; l++){
      val2 += pow(coef_y1,p2+1-l)*pow(coef_y2,l-1);
    }

    return val1*val2/((p1+1)*(p2+1));
  }


  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_Grid2DQuadMultiBlock_HO : public TestData {

    // Local variables
  public:
    Grid2D_Quad_MultiBlock_HO   MeshBlk, MeshBlk_Fine; // the mesh
    Grid2D_Quad_Block_HO   Grid;

    Grid2DTesting_Input_Parameters IP;

    int error_flag;
    Vector2D* GQPoints; // Gauss Quadrature Points
    int iCell,jCell;   // cell coordinates

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

    // allocate memory for the GQPs
    GQPoints = new Vector2D [3];

    Grid2D_Quad_Block_HO::setContourIntegrationBasedOnGaussQuadratures();
    Grid2D_Quad_Block_HO::setNoSpecialTreatmentForNumericalError();
    Grid2D_Quad_Block_HO::setDefaultPrecisionTecplotPlotting();
  }

  Data_Grid2DQuadMultiBlock_HO::~Data_Grid2DQuadMultiBlock_HO(void){
    // reset to default value
    Grid2D_Quad_Block_HO::setDefaultBoundaryRepresentation();
    Grid2D_Quad_Block_HO::setNoSpecialTreatmentForNumericalError();
    Grid2D_Quad_Block_HO::setDefaultPrecisionTecplotPlotting();
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
      // try to output the nodes
      _MeshBlk_.Output_Nodes_Tecplot_Using_IP(IP);
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

	    if (jCell == MeshBlk(iBlock,jBlock).JCu){
	      // check cells with boundary at North

	      ostmClear();
	      ostm() << "NorthBnd Check, Block (" << iBlock << "," << jBlock << "), cell[" << iCell << "," << jCell << "]";
	    
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
	      ostm() << "SouthBnd Check, Block (" << iBlock << "," << jBlock << "), cell[" << iCell << "," << jCell << "]";
	      
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
	      ostm() << "EastBnd Check, Block (" << iBlock << "," << jBlock << "), cell[" << iCell << "," << jCell << "]";

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
	      ostm() << "WestBnd Check, Block (" << iBlock << "," << jBlock << "), cell[" << iCell << "," << jCell << "]";
	      
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
	      ostm() << "InteriorCell North, Block (" << iBlock << "," << jBlock << "), cell[" << iCell << "," << jCell << "]";
	      ensure_equals(ostm().str(), MeshBlk(iBlock,jBlock).NumOfConstrainedGaussQuadPoints_North(iCell,jCell), 0);

	      ostmClear();
	      ostm() << "InteriorCell South, Block (" << iBlock << "," << jBlock << "), cell[" << iCell << "," << jCell << "]";
	      ensure_equals(ostm().str(), MeshBlk(iBlock,jBlock).NumOfConstrainedGaussQuadPoints_South(iCell,jCell), 0);

	      ostmClear();
	      ostm() << "InteriorCell East, Block (" << iBlock << "," << jBlock << "), cell[" << iCell << "," << jCell << "]";
	      ensure_equals(ostm().str(), MeshBlk(iBlock,jBlock).NumOfConstrainedGaussQuadPoints_East(iCell,jCell), 0);

	      ostmClear();
	      ostm() << "InteriorCell West, Block (" << iBlock << "," << jBlock << "), cell[" << iCell << "," << jCell << "]";
	      ensure_equals(ostm().str(), MeshBlk(iBlock,jBlock).NumOfConstrainedGaussQuadPoints_West(iCell,jCell), 0);
	    }

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
  typedef test_group<Data_Grid2DQuadMultiBlock_HO,60> Grid2DQuadMultiBlock_HO_TestSuite;
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
      RunRegressionTest("2D Laminar Flame cell", CurrentFile, MasterFile, 5.0e-11, 5.0e-11);

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
      RunRegressionTest("NACA Aerofoil", CurrentFile, MasterFile, 5.0e-10, 5.0e-10);
      
      //open file for all mesh node output
      CurrentFile  = "Current_NACA_Aerofoil_node.dat";
      MasterFile  = "NACA_Aerofoil_node.dat";
      Open_Output_File(CurrentFile);
      MultiBlockGrid.Output_Nodes_Tecplot(out());
      RunRegressionTest("NACA Aerofoil node", CurrentFile, MasterFile, 5.0e-10, 5.0e-10);

      //open file for all mesh cell output
      CurrentFile  = "Current_NACA_Aerofoil_cell.dat";
      MasterFile  = "NACA_Aerofoil_cell.dat";
      Open_Output_File(CurrentFile);
      MultiBlockGrid.Output_Cells_Tecplot(out());
      RunRegressionTest("NACA Aerofoil cell", CurrentFile, MasterFile, 5.0e-10, 5.0e-10);

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

  // Test 12:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<12>()
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

    MasterFile = "GridRectangularBox_TecplotNodes.dat";
    CurrentFile = "Current_GridRectangularBox_TecplotNodes.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // OutputMeshTecplot
      MeshBlk.Output_Nodes_Tecplot(out());

      // check
      RunRegressionTest("Nodes_Tecplot", CurrentFile, MasterFile, 1.0e-9, 1.0e-9);

    } else {
      Open_Output_File(MasterFile);
      
      // write data
      MeshBlk.Output_Nodes_Tecplot(out());
    }
  }

  // Test 13:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<13>()
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

    MasterFile = "GridRectangularBox_Tecplot.dat";
    CurrentFile = "Current_GridRectangularBox_Tecplot.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // OutputMeshTecplot
      MeshBlk.Output_Tecplot(out());

      // check
      RunRegressionTest("Output_Tecplot", CurrentFile, MasterFile, 1.0e-9, 1.0e-9);

    } else {
      Open_Output_File(MasterFile);
      
      // write data
      MeshBlk.Output_Tecplot(out());
    }
  }

  // Test 14:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<14>()
  {
    set_test_name("Grid_Rectangular_Box(), Cell properties");
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

    MasterFile = "GridRectangularBox_GeomProperties.dat";
    CurrentFile = "Current_GridRectangularBox_GeomProperties.dat";

    // Checked cell --> Block(0,0), Cell (ICl,JCu)
    iCell = MeshBlk(0,0).ICl; jCell = MeshBlk(0,0).JCl;

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // Cell Nodes
      Print_File(MeshBlk(0,0).nodeSW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeSE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNE(iCell,jCell), out());

      // Centroid
      Print_File(MeshBlk(0,0).CellCentroid(iCell,jCell), out());

      // Centroid, area, cell indexes, geometric moments
      Print_File(MeshBlk(0,0).Cell[iCell][jCell], out());

      // Cell face midpoints
      Print_File(MeshBlk(0,0).xfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceW(iCell,jCell), out());

      // Cell Gauss Points
      MeshBlk(0,0).getGaussQuadPointsFaceN(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceS(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceE(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceW(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());

      MeshBlk(0,0).getGaussQuadPointsFaceN(iCell,jCell, GQPoints, 2); 
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceS(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceE(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceW(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      // Face Lengths
      Print_File(MeshBlk(0,0).lfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceW(iCell,jCell), out());

      // Normals
      Print_File(MeshBlk(0,0).nfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceW(iCell,jCell), out());


      // check geometry
      RunRegressionTest("Output_Tecplot", CurrentFile, MasterFile, 1.0e-9, 1.0e-9);


      // Check boundaries
      ensure_equals("BndNorthSplineInfo",MeshBlk(0,0).BndNorthSplineInfo, SInfoNULL);
      ensure_equals("BndSouthSplineInfo",MeshBlk(0,0).BndSouthSplineInfo, SInfoNULL);
      ensure_equals("BndEastSplineInfo",MeshBlk(0,0).BndEastSplineInfo, SInfoNULL);
      ensure_equals("BndWestSplineInfo",MeshBlk(0,0).BndWestSplineInfo, SInfoNULL);

      ensure_equals("Constaints North", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_North(iCell,jCell), 0);
      ensure_equals("Constaints South", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_South(iCell,jCell), 0);
      ensure_equals("Constaints East", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_East(iCell,jCell), 0);
      ensure_equals("Constaints West", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_West(iCell,jCell), 0);
      ensure_equals("Total Constaints", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints(iCell,jCell), 0);

      // Check boundary conditions
      ensure_equals("BCtypeS", MeshBlk(0,0).BCtypeS[iCell], BC_REFLECTION);
      ensure_equals("BCtypeW", MeshBlk(0,0).BCtypeW[jCell], BC_REFLECTION);

      // Check update state
      ensure_equals("InteriorMesh", MeshBlk(0,0).Value_InteriorMeshUpdate_Flag(), OFF);
      ensure_equals("GhostCells", MeshBlk(0,0).Value_GhostCellsUpdate_Flag(), OFF);
      ensure_equals("CornerGhostCells", MeshBlk(0,0).Value_CornerGhostCellsUpdate_Flag(), OFF);

    } else {
      Open_Output_File(MasterFile);
      
      // write data

      // Cell Nodes
      Print_File(MeshBlk(0,0).nodeSW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeSE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNE(iCell,jCell), out());

      // Centroid
      Print_File(MeshBlk(0,0).CellCentroid(iCell,jCell), out());

      // Centroid, area, cell indexes, geometric moments
      Print_File(MeshBlk(0,0).Cell[iCell][jCell], out());

      // Cell face midpoints
      Print_File(MeshBlk(0,0).xfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceW(iCell,jCell), out());

      // Cell Gauss Points
      MeshBlk(0,0).getGaussQuadPointsFaceN(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceS(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceE(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceW(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());

      MeshBlk(0,0).getGaussQuadPointsFaceN(iCell,jCell, GQPoints, 2); 
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceS(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceE(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceW(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      // Face Lengths
      Print_File(MeshBlk(0,0).lfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceW(iCell,jCell), out());

      // Normals
      Print_File(MeshBlk(0,0).nfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceW(iCell,jCell), out());

    }
  }

  // Test 15:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<15>()
  {
    set_test_name("Grid_Circular_Cylinder()");
    RunRegression = ON;

    // Add test particular input parameters
    IP.i_Grid = GRID_CIRCULAR_CYLINDER;
    IP.Cylinder_Radius = 1;
    IP.Cylinder_Radius2 = 32;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 20;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_Smooth_Quad_Block = OFF;

    IP.i_Mesh_Stretching = ON;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
    IP.Mesh_Stretching_Factor_Idir = 1.025;
    IP.Mesh_Stretching_Factor_Jdir = 1.001;
    
    // Build the mesh
    CreateMesh(MeshBlk,IP);
    
    MasterFile = "Grid_CircularCylinder_Tecplot.dat";
    CurrentFile = "Current_Grid_CircularCylinder_Tecplot.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // OutputMeshTecplot
      MeshBlk.Output_Tecplot(out());

      // check
      RunRegressionTest("Output_Tecplot", CurrentFile, MasterFile, 1.0e-9, 1.0e-9);

    } else {
      Open_Output_File(MasterFile);
      
      // write data
      MeshBlk.Output_Tecplot(out());
    }
  }

  // Test 16:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<16>()
  {
    set_test_name("Grid_Circular_Cylinder()");
    RunRegression = ON;
 
    // Add test particular input parameters
    IP.i_Grid = GRID_CIRCULAR_CYLINDER;
    IP.Cylinder_Radius = 1;
    IP.Cylinder_Radius2 = 32;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 20;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_Smooth_Quad_Block = OFF;

    IP.i_Mesh_Stretching = ON;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
    IP.Mesh_Stretching_Factor_Idir = 1.025;
    IP.Mesh_Stretching_Factor_Jdir = 1.001;

    // Build the mesh
    CreateMesh(MeshBlk,IP);

    MasterFile = "Grid_CircularCylinder_CellsTecplot.dat";
    CurrentFile = "Current_Grid_CircularCylinder_CellsTecplot.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // OutputMeshTecplot
      MeshBlk.Output_Cells_Tecplot(out());

      // check
      RunRegressionTest("Output_Cells_Tecplot", CurrentFile, MasterFile, 1.0e-9, 1.0e-9);

    } else {
      Open_Output_File(MasterFile);
      
      // write data
      MeshBlk.Output_Cells_Tecplot(out());
    }
  }

  // Test 17:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<17>()
  {
    set_test_name("Grid_Circular_Cylinder()");
    RunRegression = ON;
 
    // Add test particular input parameters
    IP.i_Grid = GRID_CIRCULAR_CYLINDER;
    IP.Cylinder_Radius = 1;
    IP.Cylinder_Radius2 = 32;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 20;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_Smooth_Quad_Block = OFF;

    IP.i_Mesh_Stretching = ON;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
    IP.Mesh_Stretching_Factor_Idir = 1.025;
    IP.Mesh_Stretching_Factor_Jdir = 1.001;

    // Build the mesh
    CreateMesh(MeshBlk,IP);

    MasterFile = "Grid_CircularCylinder_NodesTecplot.dat";
    CurrentFile = "Current_Grid_CircularCylinder_NodesTecplot.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // OutputMeshTecplot
      MeshBlk.Output_Nodes_Tecplot(out());

      // check
      RunRegressionTest("Output_Nodes_Tecplot", CurrentFile, MasterFile, 1.0e-9, 1.0e-9);

    } else {
      Open_Output_File(MasterFile);
      
      // write data
      MeshBlk.Output_Nodes_Tecplot(out());
    }
  }

  // Test 18:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<18>()
  {
    set_test_name("Grid_Circular_Cylinder(), Interior cell properties");
    RunRegression = ON;
 
    // Add test particular input parameters
    IP.i_Grid = GRID_CIRCULAR_CYLINDER;
    IP.Cylinder_Radius = 1;
    IP.Cylinder_Radius2 = 32;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 20;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_Smooth_Quad_Block = OFF;

    IP.i_Mesh_Stretching = ON;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
    IP.Mesh_Stretching_Factor_Idir = 1.025;
    IP.Mesh_Stretching_Factor_Jdir = 1.001;

    // Build the mesh
    CreateMesh(MeshBlk,IP);

    MasterFile = "GridCircularCylinder_GeomProperties_InteriorCell.dat";
    CurrentFile = "Current_GridCircularCylinder_GeomProperties_InteriorCell.dat";

    // Checked cell --> Block(0,0), Cell (16,25)
    iCell = 16; jCell = 25;

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // Cell Nodes
      Print_File(MeshBlk(0,0).nodeSW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeSE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNE(iCell,jCell), out());

      // Centroid
      Print_File(MeshBlk(0,0).CellCentroid(iCell,jCell), out());

      // Centroid, area, cell indexes, geometric moments
      Print_File(MeshBlk(0,0).Cell[iCell][jCell], out());

      // Cell face midpoints
      Print_File(MeshBlk(0,0).xfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceW(iCell,jCell), out());

      // Cell Gauss Points
      MeshBlk(0,0).getGaussQuadPointsFaceN(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceS(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceE(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceW(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());

      MeshBlk(0,0).getGaussQuadPointsFaceN(iCell,jCell, GQPoints, 2); 
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceS(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceE(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceW(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      // Face Lengths
      Print_File(MeshBlk(0,0).lfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceW(iCell,jCell), out());

      // Normals
      Print_File(MeshBlk(0,0).nfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceW(iCell,jCell), out());


      // check geometry
      RunRegressionTest("Cell Geom Properties", CurrentFile, MasterFile, 1.0e-9, 1.0e-9);

      // Check boundaries
      ensure_equals("BndNorthSplineInfo",MeshBlk(0,0).BndNorthSplineInfo, SInfoNULL);
      ensure_equals("BndSouthSplineInfo",MeshBlk(0,0).BndSouthSplineInfo, SInfoNULL);
      ensure_equals("BndEastSplineInfo",MeshBlk(0,0).BndEastSplineInfo, SInfoNULL);
      ensure_equals("BndWestSplineInfo",MeshBlk(0,0).BndWestSplineInfo, SInfoNULL);

      ensure_equals("Constaints North", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_North(iCell,jCell), 0);
      ensure_equals("Constaints South", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_South(iCell,jCell), 0);
      ensure_equals("Constaints East", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_East(iCell,jCell), 0);
      ensure_equals("Constaints West", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_West(iCell,jCell), 0);
      ensure_equals("Total Constaints", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints(iCell,jCell), 0);

      // Check update state
      ensure_equals("InteriorMesh", MeshBlk(0,0).Value_InteriorMeshUpdate_Flag(), OFF);
      ensure_equals("GhostCells", MeshBlk(0,0).Value_GhostCellsUpdate_Flag(), OFF);
      ensure_equals("CornerGhostCells", MeshBlk(0,0).Value_CornerGhostCellsUpdate_Flag(), OFF);

    } else {
      Open_Output_File(MasterFile);
      
      // write data

      // Cell Nodes
      Print_File(MeshBlk(0,0).nodeSW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeSE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNE(iCell,jCell), out());

      // Centroid
      Print_File(MeshBlk(0,0).CellCentroid(iCell,jCell), out());

      // Centroid, area, cell indexes, geometric moments
      Print_File(MeshBlk(0,0).Cell[iCell][jCell], out());

      // Cell face midpoints
      Print_File(MeshBlk(0,0).xfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceW(iCell,jCell), out());

      // Cell Gauss Points
      MeshBlk(0,0).getGaussQuadPointsFaceN(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceS(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceE(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceW(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());

      MeshBlk(0,0).getGaussQuadPointsFaceN(iCell,jCell, GQPoints, 2); 
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceS(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceE(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceW(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      // Face Lengths
      Print_File(MeshBlk(0,0).lfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceW(iCell,jCell), out());

      // Normals
      Print_File(MeshBlk(0,0).nfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceW(iCell,jCell), out());

    }
  }

  // Test 19:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<19>()
  {
    set_test_name("Grid_Circular_Cylinder(), Boundary cell properties");
    RunRegression = ON;
 
    // Add test particular input parameters
    IP.i_Grid = GRID_CIRCULAR_CYLINDER;
    IP.Cylinder_Radius = 1;
    IP.Cylinder_Radius2 = 32;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 20;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.IncludeHighOrderBoundariesRepresentation = ON;
    IP.i_Smooth_Quad_Block = OFF;
    Grid2D_Quad_Block_HO::setContourIntegrationBasedOnLinearSegments();

    IP.i_Mesh_Stretching = ON;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
    IP.Mesh_Stretching_Factor_Idir = 1.025;
    IP.Mesh_Stretching_Factor_Jdir = 1.001;

    // Build the mesh
    CreateMesh(MeshBlk,IP);

    MeshBlk(0,0).BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux);

    MasterFile = "GridCircularCylinder_GeomProperties_BoundaryCell.dat";
    CurrentFile = "Current_GridCircularCylinder_GeomProperties_BoundaryCell.dat";

    // Checked cell --> Block(0,0), Cell (16,5)
    iCell = 16; jCell = 5;

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // Cell Nodes
      Print_File(MeshBlk(0,0).nodeSW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeSE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNE(iCell,jCell), out());

      // Centroid
      Print_File(MeshBlk(0,0).CellCentroid(iCell,jCell), out());

      // Centroid, area, cell indexes, geometric moments
      Print_File(MeshBlk(0,0).Cell[iCell][jCell], out());

      // Cell Gauss Points
      if (MeshBlk(0,0).BndSouthSplineInfo != NULL){
	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].NumOfSubIntervals(), out());

	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].GQPoint(1), out());
	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].NormalGQPoint(1), out());

	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].GQPoint(2), out());
	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].NormalGQPoint(2), out());

	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].IntLength(1), out());
      }

      if (MeshBlk(0,0).BndNorthSplineInfo != NULL){
	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].NumOfSubIntervals(), out());

	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].GQPoint(1), out());
	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].NormalGQPoint(1), out());

	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].GQPoint(2), out());
	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].NormalGQPoint(2), out());

	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].IntLength(1), out());
      }

      if (MeshBlk(0,0).BndEastSplineInfo != NULL){
	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].NumOfSubIntervals(), out());

	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].GQPoint(1), out());
	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].NormalGQPoint(1), out());

	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].GQPoint(2), out());
	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].NormalGQPoint(2), out());

	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].IntLength(1), out());
      }

      if (MeshBlk(0,0).BndWestSplineInfo != NULL){
	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].NumOfSubIntervals(), out());

	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].GQPoint(1), out());
	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].NormalGQPoint(1), out());

	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].GQPoint(2), out());
	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].NormalGQPoint(2), out());

	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].IntLength(1), out());
      }


      // check geometry
      RunRegressionTest("Cell Geom Properties", CurrentFile, MasterFile, 1.0e-9, 1.0e-9);

      // Check boundaries
      ensure_equals("BndNorthSplineInfo",MeshBlk(0,0).BndNorthSplineInfo == SInfoNULL, 0);
      ensure_equals("BndSouthSplineInfo",MeshBlk(0,0).BndSouthSplineInfo == SInfoNULL, 0);
      ensure_equals("BndEastSplineInfo",MeshBlk(0,0).BndEastSplineInfo == SInfoNULL, 1);
      ensure_equals("BndWestSplineInfo",MeshBlk(0,0).BndWestSplineInfo == SInfoNULL, 1);

      ensure_equals("Constaints North", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_North(iCell,jCell), 0);
      ensure_equals("Constaints South", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_South(iCell,jCell), 2);
      ensure_equals("Constaints East", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_East(iCell,jCell), 0);
      ensure_equals("Constaints West", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_West(iCell,jCell), 0);
      ensure_equals("Total Constaints", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints(iCell,jCell), 2);

      // Check update state
      ensure_equals("InteriorMesh", MeshBlk(0,0).Value_InteriorMeshUpdate_Flag(), OFF);
      ensure_equals("GhostCells", MeshBlk(0,0).Value_GhostCellsUpdate_Flag(), OFF);
      ensure_equals("CornerGhostCells", MeshBlk(0,0).Value_CornerGhostCellsUpdate_Flag(), OFF);

    } else {
      Open_Output_File(MasterFile);
      
      // write data

      // Cell Nodes
      Print_File(MeshBlk(0,0).nodeSW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeSE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNE(iCell,jCell), out());

      // Centroid
      Print_File(MeshBlk(0,0).CellCentroid(iCell,jCell), out());

      // Centroid, area, cell indexes, geometric moments
      Print_File(MeshBlk(0,0).Cell[iCell][jCell], out());

      // Cell Gauss Points
      if (MeshBlk(0,0).BndSouthSplineInfo != NULL){
	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].NumOfSubIntervals(), out());

	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].GQPoint(1), out());
	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].NormalGQPoint(1), out());

	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].GQPoint(2), out());
	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].NormalGQPoint(2), out());

	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].IntLength(1), out());
      }

      if (MeshBlk(0,0).BndNorthSplineInfo != NULL){
	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].NumOfSubIntervals(), out());

	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].GQPoint(1), out());
	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].NormalGQPoint(1), out());

	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].GQPoint(2), out());
	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].NormalGQPoint(2), out());

	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].IntLength(1), out());
      }

      if (MeshBlk(0,0).BndEastSplineInfo != NULL){
	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].NumOfSubIntervals(), out());

	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].GQPoint(1), out());
	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].NormalGQPoint(1), out());

	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].GQPoint(2), out());
	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].NormalGQPoint(2), out());

	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].IntLength(1), out());
      }

      if (MeshBlk(0,0).BndWestSplineInfo != NULL){
	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].NumOfSubIntervals(), out());

	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].GQPoint(1), out());
	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].NormalGQPoint(1), out());

	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].GQPoint(2), out());
	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].NormalGQPoint(2), out());

	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].IntLength(1), out());
      }
    }
  }

  // Test 20:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<20>()
  {
    set_test_name("Check output operator");
    RunRegression = ON;
 
    // Add test particular input parameters
    IP.i_Grid = GRID_CIRCULAR_CYLINDER;
    IP.Cylinder_Radius = 1;
    IP.Cylinder_Radius2 = 32;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 20;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.IncludeHighOrderBoundariesRepresentation = ON;
    IP.i_Smooth_Quad_Block = OFF;
    Grid2D_Quad_Block_HO::setContourIntegrationBasedOnLinearSegments();

    IP.i_Mesh_Stretching = ON;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
    IP.Mesh_Stretching_Factor_Idir = 1.025;
    IP.Mesh_Stretching_Factor_Jdir = 1.001;

    // Build the mesh
    CreateMesh(MeshBlk,IP);

    MasterFile = "Grid_Output_Operator.dat";
    CurrentFile = "Current_Grid_Output_Operator.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);

      out() << MeshBlk(0,0) << endl << endl;
      out() << MeshBlk(1,0) << endl << endl;

      // check operator <<
      RunRegressionTest("operator <<", CurrentFile, MasterFile, 5.0e-10, 5.0e-10);

    } else {
      Open_Output_File(MasterFile);

      out() << MeshBlk(0,0) << endl << endl;
      out() << MeshBlk(1,0) << endl << endl;
    }
  }

  // Test 21:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<21>()
  {
    set_test_name("Check operator << & >>, with curved boundaries");
 
    // Add test particular input parameters
    IP.i_Grid = GRID_CIRCULAR_CYLINDER;
    IP.Cylinder_Radius = 1;
    IP.Cylinder_Radius2 = 32;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 20;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.IncludeHighOrderBoundariesRepresentation = ON;
    IP.i_Smooth_Quad_Block = OFF;

    IP.i_Mesh_Stretching = ON;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
    IP.Mesh_Stretching_Factor_Idir = 1.025;
    IP.Mesh_Stretching_Factor_Jdir = 1.001;

    // Build the mesh
    CreateMesh(MeshBlk,IP);

    // Check output input operation
    Check_Input_Output_Operator("Grid, Block 1", MeshBlk(0,0));
    Check_Input_Output_Operator("Grid, Block 2", MeshBlk(1,0));
  }

  // Test 22:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<22>()
  {
    set_test_name("Check operator << & >>, without curved boundaries");
 
    // Add test particular input parameters
    IP.i_Grid = GRID_CIRCULAR_CYLINDER;
    IP.Cylinder_Radius = 1;
    IP.Cylinder_Radius2 = 32;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 20;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_Smooth_Quad_Block = OFF;

    IP.i_Mesh_Stretching = ON;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
    IP.Mesh_Stretching_Factor_Idir = 1.025;
    IP.Mesh_Stretching_Factor_Jdir = 1.001;

    // Build the mesh
    CreateMesh(MeshBlk,IP);

    // Check output input operation
    Check_Input_Output_Operator("Grid, Block 1", MeshBlk(0,0));
    Check_Input_Output_Operator("Grid, Block 2", MeshBlk(1,0));
  }


  // Test 23:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<23>()
  {
    set_test_name("Translate Rotate Grid, -V");
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
    IP.X_Shift.x = -6;
    IP.X_Shift.y = -6;
    IP.X_Scale = 1.0;
    IP.X_Rotate = 45.0;
    IP.IncludeHighOrderBoundariesRepresentation = ON;
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

    MasterFile = "GridRectangularBox_Translate.dat";
    CurrentFile = "Current_GridRectangularBox_Translate.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // OutputMeshTecplot
      MeshBlk.Output_Nodes_Tecplot(out());

      // check
      RunRegressionTest("Translate Grid", CurrentFile, MasterFile, 1.0e-9, 1.0e-9);

    } else {
      Open_Output_File(MasterFile);
      
      // write data
      MeshBlk.Output_Nodes_Tecplot(out());
    }

  }

  // Test 24:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<24>()
  {
    set_test_name("Copy_Quad_Block()");
    RunRegression = ON;
 
    // Add test particular input parameters
    IP.i_Grid = GRID_NACA_AEROFOIL;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 40;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 3;
    strcpy(IP.NACA_Aerofoil_Type, "0012");
    IP.Chord_Length = ONE;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_Smooth_Quad_Block = ON;

    // Build the mesh
    CreateMesh(MeshBlk,IP);

    MasterFile = "CopyGrid_NACA0012.dat";
    CurrentFile = "Current_CopyGrid_NACA0012.dat";
    
    if (RunRegression){
      Open_Output_File(CurrentFile);

      // Copy block 0
      Copy_Quad_Block(Grid, MeshBlk(2,0));

      out() << Grid << endl;

      // check
      RunRegressionTest("Copy Grid", CurrentFile, MasterFile, 1.0e-10, 1.0e-10);

    } else {
      Open_Output_File(MasterFile);
      
      // write data
      out() << MeshBlk(2,0) << endl;
    }

  }

  // Test 25:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<25>()
  {
    set_test_name("Copy_Quad_Block(), HighOrder boundaries");
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
    IP.X_Shift.x = -6;
    IP.X_Shift.y = -6;
    IP.X_Scale = 1.0;
    IP.X_Rotate = 45.0;
    IP.IncludeHighOrderBoundariesRepresentation = ON;
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

    MasterFile = "CopyGrid_HighOrder_Boundaries.dat";
    CurrentFile = "Current_CopyGrid_HighOrder_Boundaries.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // Copy block 0
      Copy_Quad_Block(Grid, MeshBlk(1,0));

      out() << Grid << endl;

      // check
      RunRegressionTest("Copy Grid", CurrentFile, MasterFile, 1.0e-9, 1.0e-9);

    } else {
      Open_Output_File(MasterFile);
      
      // write data
      out() << MeshBlk(1,0) << endl;
    }

  }


  // Test 26:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<26>()
  {
    set_test_name("Check_Quad_Block()");
 
    // Add test particular input parameters
    IP.i_Grid = GRID_CIRCULAR_CYLINDER;
    IP.Cylinder_Radius = 1;
    IP.Cylinder_Radius2 = 32;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 20;
    IP.Number_of_Cells_Jdir = 20;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_Smooth_Quad_Block = ON;

    IP.i_Mesh_Stretching = ON;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
    IP.Mesh_Stretching_Factor_Idir = 1.025;
    IP.Mesh_Stretching_Factor_Jdir = 1.001;

    // Build the mesh
    CreateMesh(MeshBlk,IP);

    int error_flag;
    error_flag = MeshBlk.Check_Multi_Block_Grid();

    ensure_equals("Check_Quad_Block", error_flag, 0);

  }

  // Test 27:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<27>()
  {
    set_test_name("Number of constrained GQPs, 'SolveRiemannProblem' boundaries, Low-order Boundaries");
 
    // Add test particular input parameters
    IP.i_Grid = GRID_RECTANGULAR_BOX;
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Blocks_Jdir = 2;
    IP.Number_of_Cells_Idir = 10;
    IP.Number_of_Cells_Jdir = 10;
    IP.Number_of_Ghost_Cells = 5;
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

    // Reset the method of computing the fluxes for all the splines
    SetFluxCalculationMethod(MeshBlk,SolveRiemannProblem,IP);

    // check --> The Result is ZERO GQP/Edge
    CheckNumberOfConstrainedGQP(MeshBlk,IP,
				0);
  }

  // Test 28:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<28>()
  {
    set_test_name("Number of constrained GQPs, 'ReconstructionBasedFlux' && Low-order boundaries, 2 GQPs/edge");
 
    // Add test particular input parameters
    IP.i_Grid = GRID_RECTANGULAR_BOX;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Cells_Idir = 10;
    IP.Number_of_Cells_Jdir = 10;
    IP.Number_of_Ghost_Cells = 5;
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

    // Reset the method of computing the fluxes for all the splines
    SetFluxCalculationMethod(MeshBlk,ReconstructionBasedFlux,IP);

    // check --> The Result is TWO GQP/Edge
    CheckNumberOfConstrainedGQP(MeshBlk,IP,
				2);
  }

  // Test 29:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<29>()
  {
    set_test_name("Number of constrained GQPs, 'ReconstructionBasedFlux' && Low-order boundaries, 1 GQP/edge");
 
    // Add test particular input parameters
    IP.i_Grid = GRID_RECTANGULAR_BOX;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Cells_Idir = 10;
    IP.Number_of_Cells_Jdir = 10;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 2;
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

    // Reset the method of computing the fluxes for all the splines
    SetFluxCalculationMethod(MeshBlk,ReconstructionBasedFlux,IP);

    // check --> The Result is ONE GQP/Edge
    CheckNumberOfConstrainedGQP(MeshBlk,IP,
				1);
  }

  // Test 30:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<30>()
  {
    set_test_name("Number of constrained GQPs, 'ReconstructionBasedFlux' && High-order boundaries, 2 GQPs/edge");
 
    // Add test particular input parameters
    IP.i_Grid = GRID_RECTANGULAR_BOX;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Cells_Idir = 10;
    IP.Number_of_Cells_Jdir = 10;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.Box_Width = 2;
    IP.Box_Height = 2;
    IP.X_Shift.x = 6;
    IP.X_Shift.y = 6;
    IP.X_Scale = 2.0;
    IP.X_Rotate = 10.0;
    IP.IncludeHighOrderBoundariesRepresentation = ON;
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

    // Reset the method of computing the fluxes for all the splines
    SetFluxCalculationMethod(MeshBlk,ReconstructionBasedFlux,IP);

    // check --> The Result is TWO GQP/Edge
    CheckNumberOfConstrainedGQP(MeshBlk,IP,
				2);
  }

  // Test 31:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<31>()
  {
    set_test_name("Number of constrained GQPs, 'ReconstructionBasedFlux' && High-order boundaries, 1 GQP/edge");
 
    // Add test particular input parameters
    IP.i_Grid = GRID_RECTANGULAR_BOX;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Cells_Idir = 10;
    IP.Number_of_Cells_Jdir = 10;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 2;
    IP.Box_Width = 2;
    IP.Box_Height = 2;
    IP.X_Shift.x = 6;
    IP.X_Shift.y = 6;
    IP.X_Scale = 2.0;
    IP.X_Rotate = 10.0;
    IP.IncludeHighOrderBoundariesRepresentation = ON;
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

    // Reset the method of computing the fluxes for all the splines
    SetFluxCalculationMethod(MeshBlk,ReconstructionBasedFlux,IP);

    // check --> The Result is ONE GQP/Edge
    CheckNumberOfConstrainedGQP(MeshBlk,IP,
				1);

  }

  // Test 32:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<32>()
  {
    set_test_name("Number of constrained GQPs, 'SolveRiemannProblem' boundaries, High-order Boundaries");
 
    // Add test particular input parameters
    IP.i_Grid = GRID_RECTANGULAR_BOX;
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Blocks_Jdir = 2;
    IP.Number_of_Cells_Idir = 10;
    IP.Number_of_Cells_Jdir = 10;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.Box_Width = 2;
    IP.Box_Height = 2;
    IP.X_Shift.x = 6;
    IP.X_Shift.y = 6;
    IP.X_Scale = 2.0;
    IP.X_Rotate = 10.0;
    IP.IncludeHighOrderBoundariesRepresentation = ON;
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

    // Reset the method of computing the fluxes for all the splines
    SetFluxCalculationMethod(MeshBlk,SolveRiemannProblem,IP);

    // check --> The Result is ZERO GQP
    CheckNumberOfConstrainedGQP(MeshBlk,IP,
				0);
  }

  // Test 33:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<33>()
  {
    set_test_name("Number of constrained GQPs, Mixed methods in corners , Low-order Boundaries");
 
    // Add test particular input parameters
    IP.i_Grid = GRID_RECTANGULAR_BOX;
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Cells_Idir = 10;
    IP.Number_of_Cells_Jdir = 10;
    IP.Number_of_Ghost_Cells = 5;
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

    // Reset the method of computing the fluxes for all the splines
    SetFluxCalculationMethod(MeshBlk,SolveRiemannProblem,IP);

    // Change the method of computing the fluxes for North and South splines
    MeshBlk(0,0).BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    MeshBlk(0,0).BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux);

    // check corners
    ostmClear();
    ostm() << "SW Corner, South";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_South(5,5) , 2);
    ostmClear();
    ostm() << "SW Corner, West";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_West(5,5) , 0);
    ostmClear();
    ostm() << "SW Corner, Total";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints(5,5) , 2);

    ostmClear();
    ostm() << "SE Corner, South";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_South(14,5) , 2);
    ostmClear();
    ostm() << "SE Corner, East";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_East(14,5) , 0);
    ostmClear();
    ostm() << "SE Corner, Total";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints(14,5) , 2);

    ostmClear();
    ostm() << "NW Corner, North";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_North(5,14) , 2);
    ostmClear();
    ostm() << "NW Corner, West";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_West(5,14) , 0);
    ostmClear();
    ostm() << "NW Corner, Total";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints(5,14) , 2);

    ostmClear();
    ostm() << "NE Corner, North";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_North(14,14) , 2);
    ostmClear();
    ostm() << "NE Corner, East";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_East(14,14) , 0);
    ostmClear();
    ostm() << "NE Corner, Total";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints(14,14) , 2);

  }

  // Test 34:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<34>()
  {
    set_test_name("Number of constrained GQPs, Mixed methods in corners , High-order Boundaries");
 
    // Add test particular input parameters
    IP.i_Grid = GRID_RECTANGULAR_BOX;
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Cells_Idir = 10;
    IP.Number_of_Cells_Jdir = 10;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 3;
    IP.Box_Width = 2;
    IP.Box_Height = 2;
    IP.X_Shift.x = 6;
    IP.X_Shift.y = 6;
    IP.X_Scale = 2.0;
    IP.X_Rotate = 10.0;
    IP.IncludeHighOrderBoundariesRepresentation = ON;
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

    // Reset the method of computing the fluxes for all the splines
    SetFluxCalculationMethod(MeshBlk,SolveRiemannProblem,IP);

    // Change the method of computing the fluxes for North and South splines
    MeshBlk(0,0).BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    MeshBlk(0,0).BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux);

    // check corners
    ostmClear();
    ostm() << "SW Corner, South";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_South(5,5) , 2);
    ostmClear();
    ostm() << "SW Corner, West";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_West(5,5) , 0);
    ostmClear();
    ostm() << "SW Corner, Total";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints(5,5) , 2);

    ostmClear();
    ostm() << "SE Corner, South";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_South(14,5) , 2);
    ostmClear();
    ostm() << "SE Corner, East";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_East(14,5) , 0);
    ostmClear();
    ostm() << "SE Corner, Total";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints(14,5) , 2);

    ostmClear();
    ostm() << "NW Corner, North";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_North(5,14) , 2);
    ostmClear();
    ostm() << "NW Corner, West";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_West(5,14) , 0);
    ostmClear();
    ostm() << "NW Corner, Total";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints(5,14) , 2);

    ostmClear();
    ostm() << "NE Corner, North";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_North(14,14) , 2);
    ostmClear();
    ostm() << "NE Corner, East";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_East(14,14) , 0);
    ostmClear();
    ostm() << "NE Corner, Total";
    ensure_equals(ostm().str(), MeshBlk(0,0).NumOfConstrainedGaussQuadPoints(14,14) , 2);
  }

  // Test 35:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<35>()
  {
    set_test_name("Grid_Rectangular_Box(), Cell properties for 5th-order");
    RunRegression = ON;
 
    // Add test particular input parameters
    IP.i_Grid = GRID_RECTANGULAR_BOX;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 10;
    IP.Number_of_Cells_Jdir = 10;
    IP.Number_of_Ghost_Cells = 2;
    IP.Space_Accuracy = 5;
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

    MasterFile = "GridRectangularBox_GeomProperties_5thOrder.dat";
    CurrentFile = "Current_GridRectangularBox_GeomProperties_5thOrder.dat";

    // Checked cell --> Block(0,0), Cell (ICl,JCu)
    iCell = MeshBlk(0,0).ICl; jCell = MeshBlk(0,0).JCl;


    if (RunRegression){
      Open_Output_File(CurrentFile);

      // Cell Nodes
      Print_File(MeshBlk(0,0).nodeSW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeSE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNE(iCell,jCell), out());

      // Centroid
      Print_File(MeshBlk(0,0).CellCentroid(iCell,jCell), out());

      // Centroid, area, cell indexes, geometric moments
      Print_File(MeshBlk(0,0).Cell[iCell][jCell], out());

      // Cell face midpoints
      Print_File(MeshBlk(0,0).xfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceW(iCell,jCell), out());

      // Cell Gauss Points
      MeshBlk(0,0).getGaussQuadPointsFaceN(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceS(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceE(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceW(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());

      MeshBlk(0,0).getGaussQuadPointsFaceN(iCell,jCell, GQPoints, 2); 
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceS(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceE(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceW(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      // Face Lengths
      Print_File(MeshBlk(0,0).lfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceW(iCell,jCell), out());

      // Normals
      Print_File(MeshBlk(0,0).nfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceW(iCell,jCell), out());


      // check geometry
      RunRegressionTest("Output_Tecplot", CurrentFile, MasterFile, 1.0e-9, 1.0e-9);


      // Check boundaries
      ensure_equals("BndNorthSplineInfo",MeshBlk(0,0).BndNorthSplineInfo, SInfoNULL);
      ensure_equals("BndSouthSplineInfo",MeshBlk(0,0).BndSouthSplineInfo, SInfoNULL);
      ensure_equals("BndEastSplineInfo",MeshBlk(0,0).BndEastSplineInfo, SInfoNULL);
      ensure_equals("BndWestSplineInfo",MeshBlk(0,0).BndWestSplineInfo, SInfoNULL);

      ensure_equals("Constaints North", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_North(iCell,jCell), 0);
      ensure_equals("Constaints South", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_South(iCell,jCell), 0);
      ensure_equals("Constaints East", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_East(iCell,jCell), 0);
      ensure_equals("Constaints West", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_West(iCell,jCell), 0);
      ensure_equals("Total Constaints", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints(iCell,jCell), 0);

      // Check boundary conditions
      ensure_equals("BCtypeS", MeshBlk(0,0).BCtypeS[iCell], BC_REFLECTION);
      ensure_equals("BCtypeW", MeshBlk(0,0).BCtypeW[jCell], BC_REFLECTION);

      // Check update state
      ensure_equals("InteriorMesh", MeshBlk(0,0).Value_InteriorMeshUpdate_Flag(), OFF);
      ensure_equals("GhostCells", MeshBlk(0,0).Value_GhostCellsUpdate_Flag(), OFF);
      ensure_equals("CornerGhostCells", MeshBlk(0,0).Value_CornerGhostCellsUpdate_Flag(), OFF);

    } else {
      Open_Output_File(MasterFile);
      
      // write data

      // Cell Nodes
      Print_File(MeshBlk(0,0).nodeSW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeSE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNE(iCell,jCell), out());

      // Centroid
      Print_File(MeshBlk(0,0).CellCentroid(iCell,jCell), out());

      // Centroid, area, cell indexes, geometric moments
      Print_File(MeshBlk(0,0).Cell[iCell][jCell], out());

      // Cell face midpoints
      Print_File(MeshBlk(0,0).xfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceW(iCell,jCell), out());

      // Cell Gauss Points
      MeshBlk(0,0).getGaussQuadPointsFaceN(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceS(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceE(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceW(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());

      MeshBlk(0,0).getGaussQuadPointsFaceN(iCell,jCell, GQPoints, 2); 
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceS(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceE(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceW(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      // Face Lengths
      Print_File(MeshBlk(0,0).lfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceW(iCell,jCell), out());

      // Normals
      Print_File(MeshBlk(0,0).nfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceW(iCell,jCell), out());

    }
  }

  // Test 36:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<36>()
  {
    set_test_name("Centroid and Area for concave quadrilateral");

    // Add test particular input parameters
    IP.i_Grid = GRID_NACA_AEROFOIL;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 40;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 3;
    strcpy(IP.NACA_Aerofoil_Type, "0012");
    IP.Chord_Length = ONE;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_Smooth_Quad_Block = ON;

    Vector2D X[4], Centroid, CentroidResult;
    double area, areaResult;
    int Info;
    int iCell, jCell; 
    iCell = 11;
    jCell = 46;

    // Define results
    CentroidResult = Vector2D(2.333333333333333333,2.25);
    areaResult = 0.375;

    // Build the mesh
    CreateMesh(MeshBlk,IP);

    // Copy block 0
    Copy_Quad_Block(Grid, MeshBlk(2,0));
    
    Grid.Node[iCell  ][jCell  ].X = Vector2D(1.0,1.0);
    Grid.Node[iCell+1][jCell  ].X = Vector2D(2.0,2.75);
    Grid.Node[iCell+1][jCell+1].X = Vector2D(4.0,1.0);
    Grid.Node[iCell  ][jCell+1].X = Vector2D(2.0,3.0);

    // Use the grid functions
    Centroid = Grid.centroid(iCell,jCell);
    area = Grid.area(iCell,jCell);

    // == check
    ensure_distance("Grid Centroid function", Centroid, CentroidResult, AcceptedError(CentroidResult));
    ensure_distance("Grid Area function", area, areaResult, AcceptedError(areaResult));
    
    
    // Use the numeric integration subroutine
    area = Grid.Integration.IntegrateFunctionOverCell(iCell,jCell,Function_area,14,area);
    Centroid = Vector2D( Grid.Integration.IntegrateFunctionOverCell(iCell,jCell,Function_XCentroid,14,area),
 			 Grid.Integration.IntegrateFunctionOverCell(iCell,jCell,Function_YCentroid,14,area) )/area;

    // == check
    ensure_distance("Integration Centroid function", Centroid, CentroidResult, AcceptedError(CentroidResult));
    ensure_distance("Integration Area function", area, areaResult, AcceptedError(areaResult));

    // Use the polyCentroid subroutine directly
    X[0] = Grid.nodeSW(iCell,jCell).X;
    X[1] = Grid.nodeSE(iCell,jCell).X;
    X[2] = Grid.nodeNE(iCell,jCell).X;
    X[3] = Grid.nodeNW(iCell,jCell).X;

    Info = polyCentroid(X, 4, Centroid, area);

    // == check
    ensure_equals("Error check", Info, 0); //< Zero value corresponds to no error
    ensure_distance("polyCentroid() function", Centroid, CentroidResult, AcceptedError(CentroidResult));
    ensure_distance("Area with polyCentroid() function", area, areaResult, AcceptedError(areaResult));
  }
  
  // Test 37:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<37>()
  {
    set_test_name("Centroid and Area for concave quadrilateral with quadAreaAndCentroid");

    // Add test particular input parameters
    IP.i_Grid = GRID_NACA_AEROFOIL;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 40;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 3;
    strcpy(IP.NACA_Aerofoil_Type, "0012");
    IP.Chord_Length = ONE;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_Smooth_Quad_Block = ON;

    Vector2D X[4], Centroid, CentroidResult;
    double area, areaResult;
    int Info;
    int iCell, jCell; 
    iCell = 11;
    jCell = 46;

    // Define results
    CentroidResult = Vector2D(2.333333333333333333,2.25);
    areaResult = 0.375;

    // Build the mesh
    CreateMesh(MeshBlk,IP);

    // Copy block 0
    Copy_Quad_Block(Grid, MeshBlk(2,0));
    
    // Create a concave quadrilateral
    Grid.Node[iCell  ][jCell  ].X = Vector2D(1.0,1.0);
    Grid.Node[iCell+1][jCell  ].X = Vector2D(2.0,2.75);
    Grid.Node[iCell+1][jCell+1].X = Vector2D(4.0,1.0);
    Grid.Node[iCell  ][jCell+1].X = Vector2D(2.0,3.0);

    // Update the cell centroid and area
    Grid.Update_Cell(iCell,jCell);

    // == check 
    ensure_distance("Centroid with quadAreaAndCentroid()",
		    Grid.Cell[iCell][jCell].Xc,
		    CentroidResult,
		    AcceptedError(CentroidResult));
    ensure_distance("Area with quadAreaAndCentroid()",
		    Grid.Cell[iCell][jCell].A,
		    areaResult,
		    AcceptedError(areaResult));
  }

  // Test 38:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<38>()
  {
    set_test_name("Detect convex quadrilateral");

    // Add test particular input parameters
    Vector2D X[4];
    int Info;

    X[0] = Vector2D(1,1);    
    X[1] = Vector2D(1.5,0.0);
    X[2] = Vector2D(4,1);
    X[3] = Vector2D(2,3);

    // Determine quadrilateral type
    ensure_equals("Find_Quadrilateral_Type() for convex quad I" , Find_Quadrilateral_Type(X[0],X[1],X[2],X[3]), 1);
    ensure_equals("Find_Quadrilateral_Type() for convex quad II", Find_Quadrilateral_Type(X[1],X[2],X[3],X[0]), 1);
  }

  // Test 39:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<39>()
  {
    set_test_name("Detect concave quadrilateral");

    // Add test particular input parameters
    Vector2D X[4];
    int Info;

    X[0] = Vector2D(1.0,1.0);
    X[1] = Vector2D(2.0,2.75);
    X[2] = Vector2D(4.0,1.0);
    X[3] = Vector2D(2.0,3.0);

    // Determine quadrilateral type
    ensure_equals("Find_Quadrilateral_Type() for concave quad I" , Find_Quadrilateral_Type(X[0],X[1],X[2],X[3]), 2);
    ensure_equals("Find_Quadrilateral_Type() for concave quad II", Find_Quadrilateral_Type(X[1],X[2],X[3],X[0]), 3);
  }

  // Test 40:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<40>()
  {
    set_test_name("Detect crossed quadrilateral");

    // Add test particular input parameters
    Vector2D X[4];
    int Info;

    X[0] = Vector2D(1.0,1.0);
    X[1] = Vector2D(3.0,2.75);
    X[2] = Vector2D(4.0,1.0);
    X[3] = Vector2D(2.0,3.0);

    // Determine quadrilateral type
    ensure_equals("Find_Quadrilateral_Type() for crossed quad I" , Find_Quadrilateral_Type(X[0],X[1],X[2],X[3]), 4);
    ensure_equals("Find_Quadrilateral_Type() for crossed quad II", Find_Quadrilateral_Type(X[1],X[2],X[3],X[0]), 4);
  }

  // Test 41:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<41>()
  {
    set_test_name("Detect quadrilateral degenerated into a triangle");

    // Add test particular input parameters
    Vector2D X[4];
    int Info;

    X[0] = Vector2D(1.0,1.0);
    X[1] = Vector2D(3.0,2.75);
    X[2] = Vector2D(3.0,2.75);
    X[3] = Vector2D(2.0,3.0);

    // Determine quadrilateral type
    ensure_equals("Find_Quadrilateral_Type() for triangle I" , Find_Quadrilateral_Type(X[0],X[1],X[2],X[3]), 1);
    ensure_equals("Find_Quadrilateral_Type() for triangle II", Find_Quadrilateral_Type(X[1],X[2],X[3],X[0]), 1);
  }

  // Test 43:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<43>()
  {
    set_test_name("Rectangular Box");

    RunRegression = ON;

    // Add test particular input parameters
    IP.i_Grid = GRID_RECTANGULAR_BOX;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 40;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_Smooth_Quad_Block = OFF;
    IP.BCs_Specified = ON;
    IP.BC_North = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_South = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_East = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_West = BC_CONSTANT_EXTRAPOLATION;

    MasterFile = "CurvedBoundaries_For_RectangularBox.dat";
    CurrentFile = "Current_CurvedBoundaries_For_RectangularBox.dat";
    
    if (RunRegression){
      IP.IncludeHighOrderBoundariesRepresentation = ON;
      Grid2D_Quad_Block_HO::setContourIntegrationBasedOnGaussQuadratures();
      Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

      // Build the mesh
      CreateMesh(MeshBlk,IP);

      // open CurrentFile
      Open_Output_File(CurrentFile);

      // output cell data
      MeshBlk.Output_Cells_Data(out());

      // == check geometric properties
      RunRegressionTest("Rectangular Box with Curved Bnds", CurrentFile, MasterFile, 5.0e-11, 5.0e-11);

    } else {
      // Build the mesh
      CreateMesh(MeshBlk,IP);

      // open MasterFile
      Open_Output_File(MasterFile);

      // output cell data
      MeshBlk.Output_Cells_Data(out());
    }
  }

  // Test 44:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<44>()
  {
    set_test_name("Grid_Circular_Cylinder(), Boundary cell properties, 5-point Gauss quad");
    RunRegression = ON;
 
    // Add test particular input parameters
    IP.i_Grid = GRID_CIRCULAR_CYLINDER;
    IP.Cylinder_Radius = 1;
    IP.Cylinder_Radius2 = 32;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 20;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.IncludeHighOrderBoundariesRepresentation = ON;
    IP.i_Smooth_Quad_Block = OFF;
    Grid2D_Quad_Block_HO::setContourIntegrationBasedOnLinearSegments();
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    IP.i_Mesh_Stretching = ON;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
    IP.Mesh_Stretching_Factor_Idir = 1.025;
    IP.Mesh_Stretching_Factor_Jdir = 1.001;

    MasterFile = "GaussQuad_For_GridCircularCylinder_GeomProperties_BoundaryCell.dat";
    CurrentFile = "Current_GaussQuad_For_GridCircularCylinder_GeomProperties_BoundaryCell.dat";

    // Checked cell --> Block(0,0), Cell (16,5)
    iCell = 16; jCell = 5;

    if (RunRegression){
      Grid2D_Quad_Block_HO::setContourIntegrationBasedOnGaussQuadratures();

      // Build the mesh
      CreateMesh(MeshBlk,IP);
      
      MeshBlk(0,0).BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
      
      Open_Output_File(CurrentFile);

      // Cell Nodes
      Print_File(MeshBlk(0,0).nodeSW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeSE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNE(iCell,jCell), out());

      // Centroid
      Print_File(MeshBlk(0,0).CellCentroid(iCell,jCell), out());

      // Centroid, area, cell indexes, geometric moments
      Print_File(MeshBlk(0,0).Cell[iCell][jCell], out());

      // Cell Gauss Points
      if (MeshBlk(0,0).BndSouthSplineInfo != NULL){
	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].NumOfSubIntervals(), out());

	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].GQPoint(1), out());
	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].NormalGQPoint(1), out());

	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].GQPoint(2), out());
	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].NormalGQPoint(2), out());

	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].IntLength(1), out());
      }

      if (MeshBlk(0,0).BndNorthSplineInfo != NULL){
	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].NumOfSubIntervals(), out());

	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].GQPoint(1), out());
	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].NormalGQPoint(1), out());

	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].GQPoint(2), out());
	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].NormalGQPoint(2), out());

	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].IntLength(1), out());
      }

      if (MeshBlk(0,0).BndEastSplineInfo != NULL){
	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].NumOfSubIntervals(), out());

	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].GQPoint(1), out());
	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].NormalGQPoint(1), out());

	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].GQPoint(2), out());
	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].NormalGQPoint(2), out());

	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].IntLength(1), out());
      }

      if (MeshBlk(0,0).BndWestSplineInfo != NULL){
	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].NumOfSubIntervals(), out());

	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].GQPoint(1), out());
	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].NormalGQPoint(1), out());

	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].GQPoint(2), out());
	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].NormalGQPoint(2), out());

	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].IntLength(1), out());
      }


      // check geometry
      RunRegressionTest("Cell Geom Properties", CurrentFile, MasterFile, 1.0e-4, 1.0e-4);

      // Check boundaries
      ensure_equals("BndNorthSplineInfo",MeshBlk(0,0).BndNorthSplineInfo == SInfoNULL, 0);
      ensure_equals("BndSouthSplineInfo",MeshBlk(0,0).BndSouthSplineInfo == SInfoNULL, 0);
      ensure_equals("BndEastSplineInfo",MeshBlk(0,0).BndEastSplineInfo == SInfoNULL, 1);
      ensure_equals("BndWestSplineInfo",MeshBlk(0,0).BndWestSplineInfo == SInfoNULL, 1);

      ensure_equals("Constaints North", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_North(iCell,jCell), 0);
      ensure_equals("Constaints South", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_South(iCell,jCell), 2);
      ensure_equals("Constaints East", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_East(iCell,jCell), 0);
      ensure_equals("Constaints West", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_West(iCell,jCell), 0);
      ensure_equals("Total Constaints", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints(iCell,jCell), 2);

      // Check update state
      ensure_equals("InteriorMesh", MeshBlk(0,0).Value_InteriorMeshUpdate_Flag(), OFF);
      ensure_equals("GhostCells", MeshBlk(0,0).Value_GhostCellsUpdate_Flag(), OFF);
      ensure_equals("CornerGhostCells", MeshBlk(0,0).Value_CornerGhostCellsUpdate_Flag(), OFF);

    } else {

      // Build the mesh
      CreateMesh(MeshBlk,IP);

      MeshBlk(0,0).BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux);

      Open_Output_File(MasterFile);
      
      // write data

      // Cell Nodes
      Print_File(MeshBlk(0,0).nodeSW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeSE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNE(iCell,jCell), out());

      // Centroid
      Print_File(MeshBlk(0,0).CellCentroid(iCell,jCell), out());

      // Centroid, area, cell indexes, geometric moments
      Print_File(MeshBlk(0,0).Cell[iCell][jCell], out());

      // Cell Gauss Points
      if (MeshBlk(0,0).BndSouthSplineInfo != NULL){
	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].NumOfSubIntervals(), out());

	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].GQPoint(1), out());
	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].NormalGQPoint(1), out());

	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].GQPoint(2), out());
	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].NormalGQPoint(2), out());

	Print_File(MeshBlk(0,0).BndSouthSplineInfo[iCell].IntLength(1), out());
      }

      if (MeshBlk(0,0).BndNorthSplineInfo != NULL){
	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].NumOfSubIntervals(), out());

	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].GQPoint(1), out());
	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].NormalGQPoint(1), out());

	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].GQPoint(2), out());
	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].NormalGQPoint(2), out());

	Print_File(MeshBlk(0,0).BndNorthSplineInfo[iCell].IntLength(1), out());
      }

      if (MeshBlk(0,0).BndEastSplineInfo != NULL){
	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].NumOfSubIntervals(), out());

	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].GQPoint(1), out());
	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].NormalGQPoint(1), out());

	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].GQPoint(2), out());
	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].NormalGQPoint(2), out());

	Print_File(MeshBlk(0,0).BndEastSplineInfo[jCell].IntLength(1), out());
      }

      if (MeshBlk(0,0).BndWestSplineInfo != NULL){
	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].NumOfSubIntervals(), out());

	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].GQPoint(1), out());
	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].NormalGQPoint(1), out());

	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].GQPoint(2), out());
	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].NormalGQPoint(2), out());

	Print_File(MeshBlk(0,0).BndWestSplineInfo[jCell].IntLength(1), out());
      }
    }
  }

  // Test 45:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<45>()
  {
    set_test_name("Check output operator, 5-point Gauss quad");
    RunRegression = ON;
 
    // Add test particular input parameters
    IP.i_Grid = GRID_CIRCULAR_CYLINDER;
    IP.Cylinder_Radius = 1;
    IP.Cylinder_Radius2 = 32;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 20;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.IncludeHighOrderBoundariesRepresentation = ON;
    IP.i_Smooth_Quad_Block = OFF;
    Grid2D_Quad_Block_HO::setContourIntegrationBasedOnLinearSegments();
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    IP.i_Mesh_Stretching = ON;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
    IP.Mesh_Stretching_Factor_Idir = 1.025;
    IP.Mesh_Stretching_Factor_Jdir = 1.001;

    MasterFile = "GaussQuad_For_Grid_Output_Operator.dat";
    CurrentFile = "Current_GaussQuad_For_Grid_Output_Operator.dat";

    if (RunRegression){
      Grid2D_Quad_Block_HO::setContourIntegrationBasedOnGaussQuadratures();

      // Build the mesh
      CreateMesh(MeshBlk,IP);

      Open_Output_File(CurrentFile);

      out() << MeshBlk(0,0) << endl << endl;
      out() << MeshBlk(1,0) << endl << endl;

      // check operator <<
      RunRegressionTest("operator <<", CurrentFile, MasterFile, 5.0e-3, 5.0e-3);

    } else {
      // Build the mesh
      CreateMesh(MeshBlk,IP);

      Open_Output_File(MasterFile);

      out() << MeshBlk(0,0) << endl << endl;
      out() << MeshBlk(1,0) << endl << endl;
    }
  }

  // Test 46:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<46>()
  {
    set_test_name("Check output operator, 3-point Gauss quad");
    RunRegression = ON;
 
    // Add test particular input parameters
    IP.i_Grid = GRID_CIRCULAR_CYLINDER;
    IP.Cylinder_Radius = 1;
    IP.Cylinder_Radius2 = 32;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 20;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.IncludeHighOrderBoundariesRepresentation = ON;
    IP.i_Smooth_Quad_Block = OFF;
    Grid2D_Quad_Block_HO::setContourIntegrationBasedOnLinearSegments();
    Spline2DInterval_HO::setThreePointGaussQuadContourIntegration();

    IP.i_Mesh_Stretching = ON;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
    IP.Mesh_Stretching_Factor_Idir = 1.025;
    IP.Mesh_Stretching_Factor_Jdir = 1.001;

    MasterFile = "GaussQuad_3Points_For_Grid_Output_Operator.dat";
    CurrentFile = "Current_GaussQuad_3Points_For_Grid_Output_Operator.dat";

    if (RunRegression){
      Grid2D_Quad_Block_HO::setContourIntegrationBasedOnGaussQuadratures();

      // Build the mesh
      CreateMesh(MeshBlk,IP);

      Open_Output_File(CurrentFile);

      out() << MeshBlk(0,0) << endl << endl;
      out() << MeshBlk(1,0) << endl << endl;

      // check operator <<
      RunRegressionTest("operator <<", CurrentFile, MasterFile, 5.0e-2, 5.0e-2);

    } else {
      // Build the mesh
      CreateMesh(MeshBlk,IP);

      Open_Output_File(MasterFile);

      out() << MeshBlk(0,0) << endl << endl;
      out() << MeshBlk(1,0) << endl << endl;
    }
  }

  // Test 47:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<47>()
  {
    set_test_name("Check output operator, 5-point Gauss quad");
    RunRegression = ON;
 
    // Add test particular input parameters
    IP.i_Grid = GRID_CIRCULAR_CYLINDER;
    IP.Cylinder_Radius = 1;
    IP.Cylinder_Radius2 = 32;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 160;
    IP.Number_of_Cells_Jdir = 80;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.IncludeHighOrderBoundariesRepresentation = ON;
    IP.i_Smooth_Quad_Block = OFF;
    Grid2D_Quad_Block_HO::setContourIntegrationBasedOnLinearSegments();
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    IP.i_Mesh_Stretching = ON;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
    IP.Mesh_Stretching_Factor_Idir = 1.025;
    IP.Mesh_Stretching_Factor_Jdir = 1.001;

    MasterFile = "GaussQuad_5Points_For_Grid_Output_Operator_LargeMesh.dat";
    CurrentFile = "Current_GaussQuad_5Points_For_Grid_Output_Operator_LargeMesh.dat";

    if (RunRegression){
      Grid2D_Quad_Block_HO::setContourIntegrationBasedOnGaussQuadratures();

      // Build the mesh
      CreateMesh(MeshBlk,IP);

      Open_Output_File(CurrentFile);

      out() << MeshBlk(0,0) << endl << endl;
      out() << MeshBlk(1,0) << endl << endl;

      // check operator <<
      RunRegressionTest("operator <<", CurrentFile, MasterFile, 5.0e-3, 5.0e-3);

    } else {
      // Build the mesh
      CreateMesh(MeshBlk,IP);

      Open_Output_File(MasterFile);

      out() << MeshBlk(0,0) << endl << endl;
      out() << MeshBlk(1,0) << endl << endl;
    }
  }

  // Test 48:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<48>()
  {
    set_test_name("Large Deformed Box");

    RunRegression = ON;

    // Add test particular input parameters
    IP.i_Grid = GRID_DEFORMED_BOX;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 40;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_Smooth_Quad_Block = OFF;
    IP.BCs_Specified = ON;
    IP.BC_North = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_South = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_East = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_West = BC_CONSTANT_EXTRAPOLATION;
    IP.IterationsOfInteriorNodesDisturbances = 300;

    IP.VertexSW = Vector2D(0.0,0.0);
    IP.VertexSE = Vector2D(4.0,1.0);
    IP.VertexNE = Vector2D(2.5,4.0);
    IP.VertexNW = Vector2D(0.5,5.0);

    IP.X_Scale = 50;

    MasterFile = "Large_Deformed_Box.dat";
    CurrentFile = "Current_Large_Deformed_Box.dat";
    char * MeshFile = "Large_Deformed_Box_Mesh.dat";
    
    if (RunRegression){
      // read the mesh
      Open_Input_File(MeshFile);
      in() >> MeshBlk;

      // check mesh
      MeshBlk.Check_Multi_Block_Grid_Completely();

      // Set high-order flags
      Grid2D_Quad_Block_HO::setHighOrderBoundaryRepresentation();
      Grid2D_Quad_Block_HO::setContourIntegrationBasedOnGaussQuadratures();
      Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

      // schedule update of all cells
      MeshBlk(0,0).Schedule_Interior_Mesh_Update();
      MeshBlk(0,0).Schedule_Ghost_Cells_Update();

      // recompute the geoemtric properties with the current method
      MeshBlk.Update_All_Cells();

      // open CurrentFile
      Open_Output_File(CurrentFile);

      // output cell data
      MeshBlk.Output_Cells_Data(out());

      // == check geometric properties
      RunRegressionTest("Large Deformed Box", CurrentFile, MasterFile, 9.0e-7, 9.0e-7);

    } else {
      // Build the mesh
      CreateMesh(MeshBlk,IP);
      MeshBlk.Check_Multi_Block_Grid_Completely();

      // open MasterFile
      Open_Output_File(MasterFile);

      // output cell data
      MeshBlk.Output_Cells_Data(out());

      // open MeshFile
      Open_Output_File(MeshFile);
      
      // output mesh data
      out() << MeshBlk;
    }
  }

  // Test 49:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<49>()
  {
    set_test_name("Large Deformed Box");

    RunRegression = ON;

    // Add test particular input parameters
    IP.i_Grid = GRID_DEFORMED_BOX;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 40;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_Smooth_Quad_Block = OFF;
    IP.BCs_Specified = ON;
    IP.BC_North = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_South = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_East = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_West = BC_CONSTANT_EXTRAPOLATION;
    IP.IterationsOfInteriorNodesDisturbances = 300;

    IP.VertexSW = Vector2D(0.0,0.0);
    IP.VertexSE = Vector2D(4.0,1.0);
    IP.VertexNE = Vector2D(2.5,4.0);
    IP.VertexNW = Vector2D(0.5,5.0);

    IP.X_Scale = 50;

    MasterFile = "Large_Deformed_Box.dat";
    CurrentFile = "Current_Large_Deformed_Box_LineBasedIntegration.dat";
    char * MeshFile = "Large_Deformed_Box_Mesh.dat";
    
    if (RunRegression){
      // read the mesh (which was generated in test 48)
      Open_Input_File(MeshFile);
      in() >> MeshBlk;

      // check mesh
      MeshBlk.Check_Multi_Block_Grid_Completely();

      // Set high-order flags
      Grid2D_Quad_Block_HO::setHighOrderBoundaryRepresentation();
      Grid2D_Quad_Block_HO::setContourIntegrationBasedOnLinearSegments();

      // schedule update of all cells
      MeshBlk(0,0).Schedule_Interior_Mesh_Update();
      MeshBlk(0,0).Schedule_Ghost_Cells_Update();

      // recompute the geoemtric properties with the current method
      MeshBlk.Update_All_Cells();

      // open CurrentFile
      Open_Output_File(CurrentFile);

      // output cell data
      MeshBlk.Output_Cells_Data(out());

      // == check geometric properties
      RunRegressionTest("Large Deformed Box", CurrentFile, MasterFile, 1.0e-6, 1.0e-6);

    } else {
      // Build the mesh
      CreateMesh(MeshBlk,IP);
      MeshBlk.Check_Multi_Block_Grid_Completely();

      // open MasterFile
      Open_Output_File(MasterFile);

      // output cell data
      MeshBlk.Output_Cells_Data(out());
    }
  }

  // Test 50:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<50>()
  {
    set_test_name("Cartesian Box with up to 4th-order moments");

    RunRegression = ON;

    // Add test particular input parameters
    IP.i_Grid = GRID_DEFORMED_BOX;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 40;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 5;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_Smooth_Quad_Block = OFF;
    IP.BCs_Specified = ON;
    IP.BC_North = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_South = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_East = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_West = BC_CONSTANT_EXTRAPOLATION;
    IP.IterationsOfInteriorNodesDisturbances = 0;

    IP.VertexSW = Vector2D(0.0,0.0);
    IP.VertexSE = Vector2D(4.0,0.0);
    IP.VertexNE = Vector2D(4.0,2.0);
    IP.VertexNW = Vector2D(0.0,2.0);

    IP.X_Scale = 50;

    MasterFile = "Cartesian_Box_4thOrder_Geom_Moments.dat";
    CurrentFile = "Current_Cartesian_Box_4thOrder_Geom_Moments.dat";
    
    if (RunRegression){
      // Build the low-order mesh (i.e. straight boundaries)
      CreateMesh(MeshBlk,IP);
      MeshBlk.Check_Multi_Block_Grid_Completely();

      // Set high-order flags (i.e. treat the mesh as if it had curved boundaries)
      Grid2D_Quad_Block_HO::setHighOrderBoundaryRepresentation();
      Grid2D_Quad_Block_HO::setContourIntegrationBasedOnGaussQuadratures(); // treat the curved boundaries with Gauss points.
      Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

      // schedule update of all cells
      MeshBlk(0,0).Schedule_Interior_Mesh_Update();
      MeshBlk(0,0).Schedule_Ghost_Cells_Update();

      // recompute the geoemtric properties with curved boundaries
      MeshBlk.Update_All_Cells();

      // open CurrentFile
      Open_Output_File(CurrentFile);

      // output cell data
      MeshBlk.Output_Cells_Data(out());

      // == check geometric properties
      RunRegressionTest("Large Rectangular Box 4th-order moments", CurrentFile, MasterFile, 1.0e-10, 1.0e-10);

    } else {
      // Build the mesh
      CreateMesh(MeshBlk,IP);
      MeshBlk.Check_Multi_Block_Grid_Completely();

      // open MasterFile
      Open_Output_File(MasterFile);

      // output cell data
      MeshBlk.Output_Cells_Data(out());
    }
  }

  // Test 51:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<51>()
  {
    set_test_name("Cartesian Box with up to 4th-order moments");

    RunRegression = ON;

    // Add test particular input parameters
    IP.i_Grid = GRID_DEFORMED_BOX;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 40;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 5;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_Smooth_Quad_Block = OFF;
    IP.BCs_Specified = ON;
    IP.BC_North = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_South = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_East = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_West = BC_CONSTANT_EXTRAPOLATION;
    IP.IterationsOfInteriorNodesDisturbances = 0;

    IP.VertexSW = Vector2D(0.0,0.0);
    IP.VertexSE = Vector2D(4.0,0.0);
    IP.VertexNE = Vector2D(4.0,2.0);
    IP.VertexNW = Vector2D(0.0,2.0);

    IP.X_Scale = 50;

    MasterFile = "Cartesian_Box_4thOrder_Geom_Moments.dat";
    CurrentFile = "Current_Cartesian_Box_4thOrder_Geom_Moments_LineBasedIntegration.dat";
    
    if (RunRegression){
      // Build the low-order mesh (i.e. straight boundaries)
      CreateMesh(MeshBlk,IP);
      MeshBlk.Check_Multi_Block_Grid_Completely();

      // Set high-order flags (i.e. treat the mesh as if it had curved boundaries)
      Grid2D_Quad_Block_HO::setHighOrderBoundaryRepresentation();
      Grid2D_Quad_Block_HO::setContourIntegrationBasedOnLinearSegments(); // treat the curved boundaries with line segments.

      // schedule update of all cells
      MeshBlk(0,0).Schedule_Interior_Mesh_Update();
      MeshBlk(0,0).Schedule_Ghost_Cells_Update();

      // recompute the geoemtric properties with curved boundaries
      MeshBlk.Update_All_Cells();

      // open CurrentFile
      Open_Output_File(CurrentFile);

      // output cell data
      MeshBlk.Output_Cells_Data(out());

      // == check geometric properties
      RunRegressionTest("Large Cartesian Box 4th-order moments", CurrentFile, MasterFile, 1.0e-10, 1.0e-10);

    } else {
      // Build the mesh
      CreateMesh(MeshBlk,IP);
      MeshBlk.Check_Multi_Block_Grid_Completely();

      // open MasterFile
      Open_Output_File(MasterFile);

      // output cell data
      MeshBlk.Output_Cells_Data(out());
    }
  }

  // Test 52:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<52>()
  {
    set_test_name("Large Deformed Box with up to 4th-order moments");

    RunRegression = ON;

    // Add test particular input parameters
    IP.i_Grid = GRID_DEFORMED_BOX;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 40;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 5;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_Smooth_Quad_Block = OFF;
    IP.BCs_Specified = ON;
    IP.BC_North = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_South = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_East = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_West = BC_CONSTANT_EXTRAPOLATION;
    IP.IterationsOfInteriorNodesDisturbances = 300;

    IP.VertexSW = Vector2D(0.0,0.0);
    IP.VertexSE = Vector2D(4.0,1.0);
    IP.VertexNE = Vector2D(2.5,4.0);
    IP.VertexNW = Vector2D(0.5,5.0);
    
    IP.X_Scale = 50;

    MasterFile = "Large_Deformed_Box_4thOrder_Geom_Moments.dat";
    CurrentFile = "Current_Large_Deformed_Box_4thOrder_Geom_Moments.dat";
    char * MeshFile = "Large_Deformed_Box_Mesh_4thOrder_Geom_Moments.dat";
    
    if (RunRegression){
      // read the mesh
      Open_Input_File(MeshFile);
      in() >> MeshBlk;

      // check mesh
      MeshBlk.Check_Multi_Block_Grid_Completely();

      // Set high-order flags
      Grid2D_Quad_Block_HO::setHighOrderBoundaryRepresentation();
      Grid2D_Quad_Block_HO::setContourIntegrationBasedOnGaussQuadratures();
      Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

      // schedule update of all cells
      MeshBlk(0,0).Schedule_Interior_Mesh_Update();
      MeshBlk(0,0).Schedule_Ghost_Cells_Update();

      // recompute the geoemtric properties with the current method
      MeshBlk.Update_All_Cells();

      // open CurrentFile
      Open_Output_File(CurrentFile);

      // output cell data
      MeshBlk.Output_Cells_Data(out());

      // == check geometric properties
      RunRegressionTest("Large Deformed Box 4th-order moments", CurrentFile, MasterFile, 1.0e-8, 1.0e-8);

    } else {
      // Build the mesh
      CreateMesh(MeshBlk,IP);
      MeshBlk.Check_Multi_Block_Grid_Completely();

      // open MasterFile
      Open_Output_File(MasterFile);

      // output cell data
      MeshBlk.Output_Cells_Data(out());

      // open MeshFile
      Open_Output_File(MeshFile);
      
      // output mesh data
      out() << MeshBlk;
    }
  }

  // Test 53:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<53>()
  {
    set_test_name("Large Deformed Box with up to 4th-order moments");

    RunRegression = ON;

    // Add test particular input parameters
    IP.i_Grid = GRID_DEFORMED_BOX;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 40;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 5;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_Smooth_Quad_Block = OFF;
    IP.BCs_Specified = ON;
    IP.BC_North = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_South = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_East = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_West = BC_CONSTANT_EXTRAPOLATION;
    IP.IterationsOfInteriorNodesDisturbances = 300;

    IP.VertexSW = Vector2D(0.0,0.0);
    IP.VertexSE = Vector2D(4.0,1.0);
    IP.VertexNE = Vector2D(2.5,4.0);
    IP.VertexNW = Vector2D(0.5,5.0);

    IP.X_Scale = 50;

    MasterFile = "Large_Deformed_Box_4thOrder_Geom_Moments.dat";
    CurrentFile = "Current_Large_Deformed_Box_4thOrder_Geom_Moments_LineBasedIntegration.dat";
    char * MeshFile = "Large_Deformed_Box_Mesh_4thOrder_Geom_Moments.dat";
    
    if (RunRegression){
      // read the mesh (which was generated in test 48)
      Open_Input_File(MeshFile);
      in() >> MeshBlk;

      // check mesh
      MeshBlk.Check_Multi_Block_Grid_Completely();

      // Set high-order flags
      Grid2D_Quad_Block_HO::setHighOrderBoundaryRepresentation();
      Grid2D_Quad_Block_HO::setContourIntegrationBasedOnLinearSegments();

      // schedule update of all cells
      MeshBlk(0,0).Schedule_Interior_Mesh_Update();
      MeshBlk(0,0).Schedule_Ghost_Cells_Update();

      // recompute the geoemtric properties with the current method
      MeshBlk.Update_All_Cells();

      // open CurrentFile
      Open_Output_File(CurrentFile);

      // output cell data
      MeshBlk.Output_Cells_Data(out());

      // == check geometric properties
      RunRegressionTest("Large Deformed Box 4th-order moments", CurrentFile, MasterFile, 1.0e-8, 1.0e-8);

    } else {
      // The mesh and the data are generated at test<52> 
    }
  }

  // Test 54:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<54>()
  {
    set_test_name("Translated Cartesian Box with up to 4th-order moments");

    RunRegression = ON;

    // Add test particular input parameters
    IP.i_Grid = GRID_DEFORMED_BOX;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 40;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 5;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_Smooth_Quad_Block = OFF;
    IP.BCs_Specified = ON;
    IP.BC_North = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_South = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_East = BC_CONSTANT_EXTRAPOLATION;
    IP.BC_West = BC_CONSTANT_EXTRAPOLATION;
    IP.IterationsOfInteriorNodesDisturbances = 0;
    IP.X_Scale = 1;

    IP.VertexSW = Vector2D(0.0,0.0);
    IP.VertexSE = Vector2D(4.0,0.0);
    IP.VertexNE = Vector2D(4.0,2.0);
    IP.VertexNW = Vector2D(0.0,2.0);

    MasterFile = "Cartesian_Box_4thOrder_Geom_Moments_Without_Translation.dat";
    CurrentFile = "Current_Translated_Cartesian_Box_4thOrder_Geom_Moments.dat";
    
    if (RunRegression){
      // Translate the box and update the mesh without curved boundaries.
      IP.X_Shift = Vector2D(5.5e6,0.01);

      // Set special treatment for this mesh
      Grid2D_Quad_Block_HO::setTreatMeshWithExtraCareForNumericalError();
      Grid2D_Quad_Block_HO::setDoublePrecisionTecplotPlotting();

      // Build the mesh
      CreateMesh(MeshBlk,IP);
      MeshBlk.Check_Multi_Block_Grid_Completely();

      // open CurrentFile
      Open_Output_File(CurrentFile);

      // output cell data
      MeshBlk.Output_Cells_Translation_Rotation_Invariant_Properties(out());
      
      // == check geometric properties
      RunRegressionTest("Translated Rectangular Box 4th-order moments", CurrentFile, MasterFile, 5.0e-10, 5.0e-10);

    } else {
      // Build the mesh
      CreateMesh(MeshBlk,IP);
      MeshBlk.Check_Multi_Block_Grid_Completely();

      // open MasterFile
      Open_Output_File(MasterFile);

      // output cell data
      MeshBlk.Output_Cells_Translation_Rotation_Invariant_Properties(out());
    }
  }

  // Test 55:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<55>()
  {
    set_test_name("Treat mesh with extra care for numerical errors");
    RunRegression = ON;
 
    // Add test particular input parameters
    IP.i_Grid = GRID_CIRCULAR_CYLINDER;
    IP.Cylinder_Radius = 1;
    IP.Cylinder_Radius2 = 32;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 40;
    IP.Number_of_Cells_Jdir = 20;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.IncludeHighOrderBoundariesRepresentation = OFF;
    IP.i_Smooth_Quad_Block = OFF;

    IP.i_Mesh_Stretching = ON;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
    IP.Mesh_Stretching_Factor_Idir = 1.025;
    IP.Mesh_Stretching_Factor_Jdir = 1.001;

    // Set special treatment for this mesh
    Grid2D_Quad_Block_HO::setTreatMeshWithExtraCareForNumericalError();

    // Build the low-order mesh
    CreateMesh(MeshBlk,IP);

    // Set the file names
    MasterFile = "GridCircularCylinder_GeomProperties_InteriorCell.dat";
    CurrentFile = "Current_GridCircularCylinder_GeomProperties_InteriorCell_SpecialTreatmentMesh.dat";

    // Checked cell --> Block(0,0), Cell (16,25)
    iCell = 16; jCell = 25;

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // Cell Nodes
      Print_File(MeshBlk(0,0).nodeSW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNW(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeSE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nodeNE(iCell,jCell), out());

      // Centroid
      Print_File(MeshBlk(0,0).CellCentroid(iCell,jCell), out());

      // Centroid, area, cell indexes, geometric moments
      Print_File(MeshBlk(0,0).Cell[iCell][jCell], out());

      // Cell face midpoints
      Print_File(MeshBlk(0,0).xfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).xfaceW(iCell,jCell), out());

      // Cell Gauss Points
      MeshBlk(0,0).getGaussQuadPointsFaceN(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceS(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceE(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());
      MeshBlk(0,0).getGaussQuadPointsFaceW(iCell,jCell, GQPoints, 1);   Print_File(GQPoints[0], out());

      MeshBlk(0,0).getGaussQuadPointsFaceN(iCell,jCell, GQPoints, 2); 
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceS(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceE(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      MeshBlk(0,0).getGaussQuadPointsFaceW(iCell,jCell, GQPoints, 2);
      Print_File(GQPoints[0], out()); Print_File(GQPoints[1], out());

      // Face Lengths
      Print_File(MeshBlk(0,0).lfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).lfaceW(iCell,jCell), out());

      // Normals
      Print_File(MeshBlk(0,0).nfaceN(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceS(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceE(iCell,jCell), out());
      Print_File(MeshBlk(0,0).nfaceW(iCell,jCell), out());

      // check geometry
      RunRegressionTest("Cell Geom Properties", CurrentFile, MasterFile, 1.0e-10, 1.0e-10);

      // Check boundaries
      ensure_equals("BndNorthSplineInfo",MeshBlk(0,0).BndNorthSplineInfo, SInfoNULL);
      ensure_equals("BndSouthSplineInfo",MeshBlk(0,0).BndSouthSplineInfo, SInfoNULL);
      ensure_equals("BndEastSplineInfo",MeshBlk(0,0).BndEastSplineInfo, SInfoNULL);
      ensure_equals("BndWestSplineInfo",MeshBlk(0,0).BndWestSplineInfo, SInfoNULL);

      ensure_equals("Constaints North", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_North(iCell,jCell), 0);
      ensure_equals("Constaints South", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_South(iCell,jCell), 0);
      ensure_equals("Constaints East", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_East(iCell,jCell), 0);
      ensure_equals("Constaints West", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints_West(iCell,jCell), 0);
      ensure_equals("Total Constaints", MeshBlk(0,0).NumOfConstrainedGaussQuadPoints(iCell,jCell), 0);

      // Check update state
      ensure_equals("InteriorMesh", MeshBlk(0,0).Value_InteriorMeshUpdate_Flag(), OFF);
      ensure_equals("GhostCells", MeshBlk(0,0).Value_GhostCellsUpdate_Flag(), OFF);
      ensure_equals("CornerGhostCells", MeshBlk(0,0).Value_CornerGhostCellsUpdate_Flag(), OFF);
    }
  }

  // Test 56:
  template<>
  template<>
  void Grid2DQuadMultiBlock_HO_object::test<56>()
  {
    set_test_name("Treat mesh with extra care for numerical errors");
    RunRegression = OFF;
 
    // Add test particular input parameters
    IP.i_Grid = GRID_CIRCULAR_CYLINDER;
    IP.Cylinder_Radius = 1;
    IP.Cylinder_Radius2 = 32;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 1440;
    IP.Number_of_Cells_Jdir = 720;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.IncludeHighOrderBoundariesRepresentation = ON;
    IP.i_Smooth_Quad_Block = OFF;

    IP.i_Mesh_Stretching = ON;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
    IP.Mesh_Stretching_Factor_Idir = 1.025;
    IP.Mesh_Stretching_Factor_Jdir = 1.001;

    // Set high-order flags
    Grid2D_Quad_Block_HO::setContourIntegrationBasedOnGaussQuadratures();
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    // Set special treatment for this mesh
    Grid2D_Quad_Block_HO::setTreatMeshWithExtraCareForNumericalError();

    // Build the low-order mesh
    CreateMesh(MeshBlk,IP);

#if 0
    // Set the file names
    MasterFile = "GridCircularCylinder_GeomProperties_InteriorCell.dat";
    CurrentFile = "Current_GridCircularCylinder_GeomProperties_InteriorCell_SpecialTreatmentMesh.dat";

    // Checked cell --> Block(0,0), Cell (16,25)
    iCell = 16; jCell = 25;
#endif

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

