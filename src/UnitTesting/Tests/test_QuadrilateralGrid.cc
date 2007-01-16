#include <tut.h>
#include "Grid/Grid2D/QuadrilateralGrid.h"
#include <limits>
#include <cstring>

using namespace std;

namespace tut
{
  struct Data_QuadrilateralGrid{
    Grid2D_Quad_Block          **MeshBlk; // the mesh
    Reconstruct2D_Input_Parameters IP; // input parameters
    char path[256];
  };

  typedef test_group<Data_QuadrilateralGrid> QuadrilateralGrid_TestGroup;
  typedef QuadrilateralGrid_TestGroup::object QuadrilateralGrid_TestObject;
}

namespace tut
{

  /*Test1: GRID_RECTANGULAR_BOX */
  template<>
  template<>
  void QuadrilateralGrid_TestObject::test<1>()
  {
    // Set Path
    //     strcpy(path,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    //     strcat(path,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcpy(path,"TUT_Tests/TestResults/");
    strcpy(IP.Output_File_Name,path);

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
    strcat(IP.Output_File_Name, "RectangularBox_Test");

    // Set Mesh
    MeshBlk = Multi_Block_Grid(MeshBlk,IP);

    if (MeshBlk != NULL) {
      if (Check_Multi_Block_Grid(MeshBlk,
				 IP.Number_of_Blocks_Idir,
				 IP.Number_of_Blocks_Jdir)) {
	std::cout << "Check_Multi_Block_Grid ERROR\n";
      } 
    } /* endif */

    // Output Mesh
    Output_Tecplot(MeshBlk, IP);

    MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
					  IP.Number_of_Blocks_Idir, 
					  IP.Number_of_Blocks_Jdir);

  }

  /*Test2: GRID_RINGLEB_FLOW */
  template<>
  template<>
  void QuadrilateralGrid_TestObject::test<2>()
  {
    // Set Path
    //     strcpy(path,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    //     strcat(path,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcpy(path,"TUT_Tests/TestResults/");
    strcpy(IP.Output_File_Name,path);

    // Set Input Parameters
    IP.i_Grid = GRID_RINGLEB_FLOW;
    strcpy(IP.Grid_Type, "Grid Ringleb Flow");
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Cells_Idir = 160;
    IP.Number_of_Cells_Jdir = 160;
    Vector2D Shift(0.0,0.0);
    IP.X_Shift = Shift;
    IP.X_Scale = 1.5;
    IP.X_Rotate = 0.0;
    IP.Number_of_Ghost_Cells = 2;
    strcat(IP.Output_File_Name, "RinglebFlow");

    // Set Mesh
 
    MeshBlk = Multi_Block_Grid(MeshBlk,IP);

    if (MeshBlk != NULL) {
      if (Check_Multi_Block_Grid(MeshBlk,
				 IP.Number_of_Blocks_Idir,
				 IP.Number_of_Blocks_Jdir)) {
	std::cout << "Check_Multi_Block_Grid ERROR\n";
      } 
    } /* endif */

    // Output Mesh
    Output_Tecplot(MeshBlk, IP);

    MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
					  IP.Number_of_Blocks_Idir, 
					  IP.Number_of_Blocks_Jdir);
  }

  /*Test3: GRID_RECTANGULAR_BOX */
  template<>
  template<>
  void QuadrilateralGrid_TestObject::test<3>()
  {
    // Set Path
    strcpy(path,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(path,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcpy(IP.Output_File_Name,path);

    // Set Input Parameters
    IP.i_Grid = GRID_RECTANGULAR_BOX;
    strcpy(IP.Grid_Type, "Grid Rectangular Box");
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Box_Width = 2.0;
    IP.Box_Height = 3.0;
    IP.Number_of_Cells_Idir = 9;
    IP.Number_of_Cells_Jdir = 15;
    Vector2D Shift(1.0,0.0);
    IP.X_Shift = Shift;
    IP.X_Scale = 1.0;
    IP.X_Rotate = 0.0;
    IP.Number_of_Ghost_Cells = 2;
    strcat(IP.Output_File_Name, "New_Test3_RectangularBox");

    // Set Mesh
 
    MeshBlk = Multi_Block_Grid(MeshBlk,IP);

    if (MeshBlk != NULL) {
      if (Check_Multi_Block_Grid(MeshBlk,
				 IP.Number_of_Blocks_Idir,
				 IP.Number_of_Blocks_Jdir)) {
	std::cout << "Check_Multi_Block_Grid ERROR\n";
      } 
    } /* endif */

    // Output Mesh
    Output_Tecplot(MeshBlk,IP);
    Output_Nodes_Tecplot(MeshBlk,IP);
    Output_Cells_Tecplot(MeshBlk,IP);

   MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
					  IP.Number_of_Blocks_Idir, 
					  IP.Number_of_Blocks_Jdir);
  }


  /*Test4: GRID_NACA_AEROFOIL */
  template<>
  template<>
  void QuadrilateralGrid_TestObject::test<4>()
  {
    // Set Path
    strcpy(path,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(path,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcpy(IP.Output_File_Name,path);

    // Set Input Parameters
    IP.i_Grid = GRID_RECTANGULAR_BOX;
    strcpy(IP.Grid_Type, "NACA Aerofoil");
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Blocks_Jdir = 1;
    strcpy(IP.NACA_Aerofoil_Type, "0012");
    IP.i_Grid = GRID_NACA_AEROFOIL;
    IP.Chord_Length = ONE;
    IP.Number_of_Cells_Idir = 100;
    IP.Number_of_Cells_Jdir = 100;
    IP.X_Shift = Vector2D_ZERO;
    IP.X_Scale = 1.0;
    IP.X_Rotate = 0.0;
    IP.Number_of_Ghost_Cells = 2;
    strcat(IP.Output_File_Name, "NACA0012_Test");

    // Set Mesh
 
    MeshBlk = Multi_Block_Grid(MeshBlk,IP);

    if (MeshBlk != NULL) {
      if (Check_Multi_Block_Grid(MeshBlk,
				 IP.Number_of_Blocks_Idir,
				 IP.Number_of_Blocks_Jdir)) {
	std::cout << "Check_Multi_Block_Grid ERROR\n";
      } 
    } /* endif */

    // Output Mesh
    Output_Tecplot(MeshBlk,IP);
    Output_Nodes_Tecplot(MeshBlk,IP);
    Output_Cells_Tecplot(MeshBlk,IP);

    MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
					  IP.Number_of_Blocks_Idir, 
					  IP.Number_of_Blocks_Jdir);

  }


  /*Test5: GRID_RINGLEB_FLOW */
  template<>
  template<>
  void QuadrilateralGrid_TestObject::test<5>()
  {
    // Set Path
    strcpy(path,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(path,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcpy(IP.Output_File_Name,path);

    // Set Input Parameters
    IP.i_Grid = GRID_RINGLEB_FLOW;
    strcpy(IP.Grid_Type, "Grid Ringleb Flow");
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Cells_Idir = 16;
    IP.Number_of_Cells_Jdir = 16;
    IP.X_Shift = Vector2D_ZERO;
    IP.X_Scale = 1.5;
    IP.X_Rotate = 0.0;
    IP.Number_of_Ghost_Cells = 4;
    strcat(IP.Output_File_Name, "RinglebFlow_4GhostCells");

    // Set Mesh
 
    MeshBlk = Multi_Block_Grid(MeshBlk,IP);

    if (MeshBlk != NULL) {
      if (Check_Multi_Block_Grid(MeshBlk,
				 IP.Number_of_Blocks_Idir,
				 IP.Number_of_Blocks_Jdir)) {
	std::cout << "Check_Multi_Block_Grid ERROR\n";
      } 
    } /* endif */

    Grid2D_Quad_Block DoubleGrid, HalfGrid;
    Grid2D_Quad_Block RefineGrid0, RefineGrid1, RefineGrid2, RefineGrid3;
    Grid2D_Quad_Block CoarsenGrid;

    Grid2D_Quad_Block MeshBlkRef;
    Copy_Quad_Block(MeshBlkRef,MeshBlk[0][0]);

    ofstream mesh_file;
    char *mesh_file_name_ptr;
    char mesh_file_name[256];
    mesh_file_name_ptr = mesh_file_name;

    // Double Grid
    strcpy(mesh_file_name, IP.Output_File_Name);
    strcat(mesh_file_name,"_DoubleMesh_nodes.dat");
    mesh_file.open(mesh_file_name_ptr, ios::out);
    if (mesh_file.bad()) std::cout << "Problems\n";

    Double_Mesh_Resolution(DoubleGrid,MeshBlkRef);

    // Output Mesh
    Output_Nodes_Tecplot(DoubleGrid,0,1, mesh_file);
    mesh_file.close();

    // Half Grid
    strcpy(mesh_file_name, IP.Output_File_Name);
    strcat(mesh_file_name,"_HalfMesh_nodes.dat");
    mesh_file.open(mesh_file_name_ptr, ios::out);
    if (mesh_file.bad()) std::cout << "Problems\n";

    Half_Mesh_Resolution(HalfGrid,MeshBlkRef);

    // Output Mesh
    Output_Nodes_Tecplot(HalfGrid,0,1, mesh_file);
    mesh_file.close();

    // Refine Grid
    strcpy(mesh_file_name, IP.Output_File_Name);
    strcat(mesh_file_name,"_RefineMesh1_nodes.dat");
    mesh_file.open(mesh_file_name_ptr, ios::out);
    if (mesh_file.bad()) std::cout << "Problems\n";

    Refine_Mesh(RefineGrid1,MeshBlkRef,1);

    // Output Mesh
    Output_Nodes_Tecplot(RefineGrid1,0,1, mesh_file);
    mesh_file.close();

    // Refine Grid
    strcpy(mesh_file_name, IP.Output_File_Name);
    strcat(mesh_file_name,"_RefineMesh0_nodes.dat");
    mesh_file.open(mesh_file_name_ptr, ios::out);
    if (mesh_file.bad()) std::cout << "Problems\n";

    Refine_Mesh(RefineGrid0,MeshBlkRef,0);

    // Output Mesh
    Output_Nodes_Tecplot(RefineGrid0,0,1, mesh_file);
    mesh_file.close();
    
    strcpy(mesh_file_name, IP.Output_File_Name);
    strcat(mesh_file_name,"_RefineMesh.dat");
    mesh_file.open(mesh_file_name_ptr, ios::out);
    if (mesh_file.bad()) std::cout << "Problems\n";
    Output_Tecplot(RefineGrid1,0,1, mesh_file);
    mesh_file.close();


    // Coarsen Grid
    strcpy(mesh_file_name, IP.Output_File_Name);
    strcat(mesh_file_name,"_CoarsenMesh_nodes.dat");
    mesh_file.open(mesh_file_name_ptr, ios::out);
    if (mesh_file.bad()) std::cout << "Problems\n";

    Refine_Mesh(RefineGrid2,MeshBlkRef,2);
    Refine_Mesh(RefineGrid3,MeshBlkRef,3);

    Coarsen_Mesh(CoarsenGrid,RefineGrid0,RefineGrid1,RefineGrid2,RefineGrid3);

    Output_Nodes_Tecplot(CoarsenGrid,0,1, mesh_file);
    mesh_file.close();

    Output_Nodes_Tecplot(MeshBlk,IP);
    Output_Tecplot(MeshBlk,IP);

   MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
					  IP.Number_of_Blocks_Idir, 
					  IP.Number_of_Blocks_Jdir);
  }


  /*Test6: Update_Corner_Ghost_Nodes(Grid2D_Quad_Block &Grid) */
  template<>
  template<>
  void QuadrilateralGrid_TestObject::test<6>()
  {
    // Set Path
    strcpy(path,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(path,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcpy(IP.Output_File_Name,path);

    // Set Input Parameters
    IP.i_Grid = GRID_RINGLEB_FLOW;
    strcpy(IP.Grid_Type, "Grid Ringleb Flow");
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Cells_Idir = 16;
    IP.Number_of_Cells_Jdir = 16;
    IP.X_Shift = Vector2D_ZERO;
    IP.X_Scale = 1.5;
    IP.X_Rotate = 0.0;
    IP.Number_of_Ghost_Cells = 2;
    strcat(IP.Output_File_Name, "RinglebFlow");

    // Set Mesh
 
    MeshBlk = Multi_Block_Grid(MeshBlk,IP);

    if (MeshBlk != NULL) {
      if (Check_Multi_Block_Grid(MeshBlk,
				 IP.Number_of_Blocks_Idir,
				 IP.Number_of_Blocks_Jdir)) {
	std::cout << "Check_Multi_Block_Grid ERROR\n";
      } 
    } /* endif */

    Grid2D_Quad_Block MeshBlkRef;
    Copy_Quad_Block(MeshBlkRef,MeshBlk[0][0]);

    ofstream mesh_file;
    char *mesh_file_name_ptr;
    char mesh_file_name[256];
    mesh_file_name_ptr = mesh_file_name;

    strcpy(mesh_file_name, IP.Output_File_Name);
    strcat(mesh_file_name,"UpdateCornerGhostNodes_nodes.dat");
    mesh_file.open(mesh_file_name_ptr, ios::out);
    if (mesh_file.bad()) std::cout << "Problems\n";

    Update_Corner_Ghost_Nodes(MeshBlkRef);

    // Output Mesh
    Output_Nodes_Tecplot(MeshBlkRef,0,1,mesh_file);
    mesh_file.close();

   MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
					  IP.Number_of_Blocks_Idir, 
					  IP.Number_of_Blocks_Jdir);
  }


  /*Test7: DistanceFromPointToLine */
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

namespace tut
{
  QuadrilateralGrid_TestGroup QuadrilateralGrid_Test("QuadrilateralGrid_Test");
}
