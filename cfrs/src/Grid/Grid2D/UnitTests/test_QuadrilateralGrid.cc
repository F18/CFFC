/*!\file test_PointWiseSolution.cc
   \brief Regression tests for PointWiseSolution class.*/

/* Include required C++ libraries. */

/* Using std namespace functions */

/* Include CFFC header files */
#include "TestData.h"
#include "Grid/Grid2D/QuadrilateralGrid.h"

using namespace std;

namespace tut
{

  // **********************************************
  //     TEST SUITE: QuadrilateralGrid_TestGroup
  // **********************************************
  // Data used for testing
  // **********************
  class Data_QuadrilateralGrid: public TestData{
  public:
    Grid2D_Quad_Block          **MeshBlk; // the mesh
    Reconstruct2D_Input_Parameters IP; // input parameters
    char path[256];

    int error_flag;

    Data_QuadrilateralGrid(void){
      set_test_suite_path("Grid/Grid2D/UnitTests");
    }

    ~Data_QuadrilateralGrid(){
      // Deallocate the mesh if previously allocated 
      MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
					    IP.Number_of_Blocks_Idir, 
					    IP.Number_of_Blocks_Jdir);
    };

    // Create Mesh
    void CreateMesh(Grid2D_Quad_Block **& _MeshBlk_,
		    Reconstruct2D_Input_Parameters & IP) throw(std::runtime_error);

  };


  void Data_QuadrilateralGrid::CreateMesh(Grid2D_Quad_Block **& _MeshBlk_,
					  Reconstruct2D_Input_Parameters & IP) throw(std::runtime_error){

    // Generate the mesh
    _MeshBlk_ = Multi_Block_Grid(_MeshBlk_, 
				 IP);
    
    if (_MeshBlk_ == NULL) {
      error_flag = 1;
    } else if (Check_Multi_Block_Grid(_MeshBlk_,
				      IP.Number_of_Blocks_Idir,
				      IP.Number_of_Blocks_Jdir)) {
      error_flag = 1;
    } else {
      error_flag = 0;
    } /* endif */
    
    if (error_flag) {
      throw runtime_error("CreateMesh() ERROR: Unable to create valid Euler2D multi-block mesh.");
    }
  }



  // Define Test Group && Test Object
  typedef test_group<Data_QuadrilateralGrid> QuadrilateralGrid_TestGroup;
  typedef QuadrilateralGrid_TestGroup::object QuadrilateralGrid_TestObject;


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


  /*Test1: */
  template<>
  template<>
  void QuadrilateralGrid_TestObject::test<1>()
  {

    set_test_name("GRID_RECTANGULAR_BOX");
    RunRegression = ON;

    // Set Input Parameters
    IP.i_Grid = GRID_RECTANGULAR_BOX;
    strcpy(IP.Grid_Type, "Grid Rectangular Box");
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Box_Width = 2.0;
    IP.Box_Height = 2.0;
    IP.Number_of_Cells_Idir = 20;
    IP.Number_of_Cells_Jdir = 20;
    Vector2D Shift(0.0,0.0);
    IP.X_Shift = Shift;
    IP.X_Scale = 1.0;
    IP.X_Rotate = 0.0;
    IP.Number_of_Ghost_Cells = 2;
    IP.NumOfIter_UnsmoothMesh = 100;

    // Build the mesh
    CreateMesh(MeshBlk,IP);

    MasterFile = "GridRectangularBox_TecplotCells_Master.dat";
    CurrentFile = "Current.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // OutputMeshTecplot
      Output_Cells_Tecplot(MeshBlk,IP.Number_of_Blocks_Idir,IP.Number_of_Blocks_Jdir,out());

      // check
      RunRegressionTest("Cells_Tecplot", CurrentFile, MasterFile, 1.0e-13, 1.0e-13);

    } else {
      Open_Output_File(MasterFile);
      
      // write data
      Output_Cells_Tecplot(MeshBlk,IP.Number_of_Blocks_Idir,IP.Number_of_Blocks_Jdir,out());
    }

  }


  /*Test3: DistanceFromPointToLine */
  template<>
  template<>
  void QuadrilateralGrid_TestObject::test<7>()
  {
    double Result, AnalyticResult, tol;
    int digits;
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);


    Vector2D Point1(0.0, 0.0), Point2(1.0,0.0), Point3(0.0,1.0);
    // Test 1
    AnalyticResult = sqrt(2.0)/2.0 ; 
    Result = DistanceFromPointToLine(Point1,Point2,Point3);

    ensure("DistanceFromPointToLine::Test1", fabs(Result - AnalyticResult)/AnalyticResult <= tol);

    Result = DistanceFromPointToLine(Point1,Point3,Point2);
    ensure("DistanceFromPointToLine::Test1 reversed", fabs(Result - AnalyticResult)/AnalyticResult <= tol);


    // Test 2
    Point1.x = 2.5;
    Point1.y = 4.1;

    Point2.x = -1.1;
    Point2.y = -3.0;

    Point3.x = -1.1;
    Point3.y = 3.0;

    AnalyticResult = 3.6;
    Result = DistanceFromPointToLine(Point1,Point2,Point3);
    ensure("DistanceFromPointToLine::Test2", fabs(Result - AnalyticResult)/AnalyticResult <= tol);

    Result = DistanceFromPointToLine(Point1,Point3,Point2);
    ensure("DistanceFromPointToLine::Test2 reversed", fabs(Result - AnalyticResult)/AnalyticResult <= tol);

  }


}

tut::QuadrilateralGrid_TestGroup QuadrilateralGrid_Test("Grid2D:QuadrilateralGrid");

