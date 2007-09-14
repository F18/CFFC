#include <tut.h>
#include "Reconstruction/ReconstructionFunction.h"
#include "TestFunctions/TestFunctions_2D.h"
#include "CFD/Gaussian2DState.h"
#include <limits>
#include "Reconstruction/ComputationalDomain.h"

using namespace std;

namespace tut
{
  struct Data_Reconstruction{

    typedef ComputationalDomain<OneD,Cell1D_NonUniform,double> CompDomainType1D;

    vector<int> i_index, j_index, k_index;
    int iCell,jCell,kCell;
    ComputationalDomain<> SolnBlk;
    ComputationalDomain<TwoD,Cell2D_Quad,Gaussian2D_pState> SolnBlkGauss;
    CompDomainType1D Domain1D;
    Grid2D_Quad_Block   Grid;
    Grid2D_Quad_Block   **MeshBlk; // the mesh
    Reconstruct2D_Input_Parameters IP; // input parameters
    Reconstruct1D_Input_Parameters IP1D;
  };

  typedef test_group<Data_Reconstruction> Reconstruction_TestGroup;
  typedef Reconstruction_TestGroup::object Reconstruction_TestObject;
}

namespace tut
{

  /*Test1: */
  template<>
  template<>
  void Reconstruction_TestObject::test<1>()
  {

    // Set the domain

    // Set Path
    char path[256];
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
    Vector2D Shift(0.0,0.0);
    IP.X_Shift = Shift;
    IP.X_Scale = 1.5;
    IP.X_Rotate = 0.0;
    IP.Number_of_Ghost_Cells = 3;
    strcat(IP.Output_File_Name, "RinglebFlow");
    IP.Number_of_SubGrid_Points_Idir = 4;
    IP.Number_of_SubGrid_Points_Jdir = 4;
    IP.Reconstruction_Order = 3;
    IP.ReconstructionMethod() = DD_ENO ;
    strcat(IP.Output_File_Name, "ReconstructionRinglebFlow");

    // Set Grid
    MeshBlk = Multi_Block_Grid(MeshBlk,IP);
    if (MeshBlk != NULL) {
      if (Check_Multi_Block_Grid(MeshBlk,
				 IP.Number_of_Blocks_Idir,
				 IP.Number_of_Blocks_Jdir)) {
	std::cout << "Check_Multi_Block_Grid ERROR\n";
      } 
    } /* endif */

    Copy_Quad_Block(Grid,MeshBlk[0][0]);

    // Set ComputationalDomain
    SolnBlk.SetDomain(Grid,IP);

    // Set the stencil -> Reconstruction with one ring around cell (3,3)
    iCell = 3;
    jCell = 3;

    MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
					  IP.Number_of_Blocks_Idir, 
					  IP.Number_of_Blocks_Jdir);

  }

  /*Test2: Check the LHS of the linear system in the reconstruction function */
  template<>
  template<>
  void Reconstruction_TestObject::test<2>()
  {
    // Set the domain

    // Set Path
    char path[256];
    strcpy(path,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(path,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcpy(IP.Output_File_Name,path);

    // Set Input Parameters
    IP.i_Grid = GRID_RECTANGULAR_BOX;
    strcpy(IP.Grid_Type, "Grid Rectangular Box");
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Cells_Idir = 3;
    IP.Number_of_Cells_Jdir = 5;
    IP.Box_Width = 2.0;
    IP.Box_Height = 3.0;
    Vector2D Shift(0.0,0.0);
    IP.X_Shift = Shift;
    IP.X_Scale = 1.0;
    IP.X_Rotate = 0.0;
    IP.Number_of_Ghost_Cells = 3;
    strcat(IP.Output_File_Name, "RectangularBox");
    IP.Number_of_SubGrid_Points_Idir = 4;
    IP.Number_of_SubGrid_Points_Jdir = 4;
    IP.Reconstruction_Order = 3;
    strcat(IP.Output_File_Name, "ReconstructionOverRectangularBox");

    // Set Grid
    MeshBlk = Multi_Block_Grid(MeshBlk,IP);
    if (MeshBlk != NULL) {
      if (Check_Multi_Block_Grid(MeshBlk,
				 IP.Number_of_Blocks_Idir,
				 IP.Number_of_Blocks_Jdir)) {
	std::cout << "Check_Multi_Block_Grid ERROR\n";
      } 
    } /* endif */

    Copy_Quad_Block(Grid,MeshBlk[0][0]);

    // Set ComputationalDomain
    SolnBlk.SetDomain(Grid,IP);

    // Set the stencil -> Reconstruction with one ring around cell (3,3)
    iCell = 3;
    jCell = 3;

    i_index.reserve(9);
    j_index.reserve(9);
    for(int i=1; i<=3; ++i){
      // Set i_index
      i_index.push_back(2);
      i_index.push_back(3);
      i_index.push_back(4);
      // Set j_index
      j_index.push_back(i+1);
      j_index.push_back(i+1);
      j_index.push_back(i+1);
    }


    // Compute the theoretical coefficients

    TaylorDerivativesContainer<TwoD,double> & NCell = SolnBlk( iCell, jCell ).CellGeomCoeff(); 
    /*coeff. of neighbour cell */

    double deltaX = IP.Box_Width/IP.Number_of_Cells_Idir;
    double deltaY = IP.Box_Height/IP.Number_of_Cells_Jdir;

    vector<double> CartesianCoeff;
    // cell (4,4)
    int Memb = (IP.Reconstruction_Order+1)*(IP.Reconstruction_Order+2)/2;
    CartesianCoeff.reserve(Memb);

    int p2_limit = IP.Reconstruction_Order;

    for (int p1=0; p1<=IP.Reconstruction_Order; ++p1, p2_limit--){
      for (int p2=0; p2<=p2_limit; ++p2){
	CartesianCoeff.push_back(GeomCoeffCartesian(p1,p2,deltaX,deltaY,deltaX,deltaY) - NCell(p1,p2));
      }
    }

    for (int i=0; i<=Memb-1; ++i){
      std::cout << CartesianCoeff[i] << std::endl;
    }

    // This result is compared visually. Print the matrix A in the Weighted_ENO_Reconstruction() subroutine and
    // check the last line. It should be equal to the CartesianCoeff except for the first one (p1=0,p2=0).

    MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
					  IP.Number_of_Blocks_Idir, 
					  IP.Number_of_Blocks_Jdir);
  }

#if 0
  /*Test3: */
  template<>
  template<>
  void Reconstruction_TestObject::test<3>()
  {
    // Set the domain

    // Set Path
    char path[256];
    strcpy(path,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(path,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcpy(IP.Output_File_Name,path);

    // Set Input Parameters
    IP.i_Grid = GRID_RECTANGULAR_BOX;
    strcpy(IP.Grid_Type, "Grid Rectangular Box");
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Cells_Idir = 3;
    IP.Number_of_Cells_Jdir = 5;
    IP.Box_Width = 2.0;
    IP.Box_Height = 3.0;
    Vector2D Shift(0.0,0.0);
    IP.X_Shift = Shift;
    IP.X_Scale = 1.0;
    IP.X_Rotate = 0.0;
    IP.Number_of_Ghost_Cells = 3;
    strcat(IP.Output_File_Name, "RectangularBox");
    IP.Number_of_SubGrid_Points_Idir = 4;
    IP.Number_of_SubGrid_Points_Jdir = 4;
    IP.Reconstruction_Order = 2;
    IP.ReconstructionMethod() = DD_ENO;
    strcat(IP.Output_File_Name, "Gaussian2DReconstructionOverRectangularBox");

    // Set Grid
    MeshBlk = Multi_Block_Grid(MeshBlk,IP);
    if (MeshBlk != NULL) {
      if (Check_Multi_Block_Grid(MeshBlk,
				 IP.Number_of_Blocks_Idir,
				 IP.Number_of_Blocks_Jdir)) {
	std::cout << "Check_Multi_Block_Grid ERROR\n";
      } 
    } /* endif */

    Copy_Quad_Block(Grid,MeshBlk[0][0]);

    // Set ComputationalDomain
    SolnBlkGauss.SetDomain(Grid,IP);

    // Set the stencil -> Reconstruction with one ring around cell (3,3)
    iCell = 3;
    jCell = 3;

    i_index.reserve(9);
    j_index.reserve(9);
    for(int i=1; i<=3; ++i){
      // Set i_index
      i_index.push_back(2);
      i_index.push_back(3);
      i_index.push_back(4);
      // Set j_index
      j_index.push_back(i+1);
      j_index.push_back(i+1);
      j_index.push_back(i+1);
    }

    // Assign solution data
    for(int i=0; i<=SolnBlkGauss.iLastCell(); ++i)
      for(int j=0; j<=SolnBlkGauss.jLastCell(); ++j){
	for(int param=1; param<=8; ++param)
	  SolnBlkGauss(i,j).CellSolution(param) = i+j+param;
      }

    // Compute the theoretical coefficients

    MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
					  IP.Number_of_Blocks_Idir, 
					  IP.Number_of_Blocks_Jdir);
  }

  /*Test4: */
  template<>
  template<>
  void Reconstruction_TestObject::test<4>()
  {
    // Set the domain
    try{
    // Set Path
    char path[256];
    strcpy(path,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(path,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcpy(IP.Output_File_Name,path);

    // Set Input Parameters
    IP.i_Grid = GRID_RECTANGULAR_BOX;
    strcpy(IP.Grid_Type, "Grid Rectangular Box");
    IP.Number_of_Blocks_Idir = 1;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Cells_Idir = 6;
    IP.Number_of_Cells_Jdir = 10;
    IP.Box_Width = 2.0;
    IP.Box_Height = 3.0;
    Vector2D Shift(0.0,0.0);
    IP.X_Shift = Shift;
    IP.X_Scale = 1.0;
    IP.X_Rotate = 0.0;
    IP.Number_of_Ghost_Cells = 3;
    strcat(IP.Output_File_Name, "RectangularBox");
    IP.Number_of_SubGrid_Points_Idir = 4;
    IP.Number_of_SubGrid_Points_Jdir = 4;
    IP.Reconstruction_Order = 6;
    IP.TestF = Test_Default2D;
    IP.ReconstructionMethod() = DD_ENO;
    strcat(IP.Output_File_Name, "ReconstructionOverRectangularBox");

    // Set Grid
    MeshBlk = Multi_Block_Grid(MeshBlk,IP);
    if (MeshBlk != NULL) {
      if (Check_Multi_Block_Grid(MeshBlk,
				 IP.Number_of_Blocks_Idir,
				 IP.Number_of_Blocks_Jdir)) {
	std::cout << "Check_Multi_Block_Grid ERROR\n";
      } 
    } /* endif */

    Copy_Quad_Block(Grid,MeshBlk[0][0]);

    // Set ComputationalDomain
    SolnBlk.SetDomain(Grid,IP);
    // Compute the initial solution
    SolnBlk.SetInitialData(IP);
    SolnBlk.UpdateSubgridSolution();

    // Print solution at the cell center
    // Set Path
    ofstream output_file;
    char output_file_name[256];
    strcpy(output_file_name,IP.Output_File_Name);
    strcat(output_file_name,".dat");

    output_file.open(output_file_name, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    SolnBlk.OutputSolutionCellTecplot(output_file);
    output_file.close();

    // Output solution in the subgrid
    strcpy(output_file_name,IP.Output_File_Name);
    strcat(output_file_name,"_Subgrid.dat");

    output_file.open(output_file_name, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    SolnBlk.OutputSolutionNodesTecplot(output_file);

    output_file.close();

    // Output solution of the reconstructed function in the subgrid
    strcpy(output_file_name,IP.Output_File_Name);
    strcat(output_file_name,"_ReconstructSubgrid.dat");

    output_file.open(output_file_name, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    // reconstruct solution
    SolnBlk.ReconstructSolution(IP);
    SolnBlk.UpdateSubgridSolution();
    SolnBlk.AssessReconstructionAccuracy(IP);
    SolnBlk.PrintErrorNorms();
    //    Print(SolnBlk);
    SolnBlk.OutputSolutionNodesTecplot(output_file);
    output_file.close();

    MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk, 
					  IP.Number_of_Blocks_Idir, 
					  IP.Number_of_Blocks_Jdir);
    }
    catch(...){
      Print("Ex")
    }

  }

  /*Test5: DivDifference() */
  template<>
  template<>
  void Reconstruction_TestObject::test<5>()
  {
    // Set the domain

    // Set Path
    char path[256];
    strcpy(path,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(path,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcpy(IP.Output_File_Name,path);

    // Set Input Parameters
    IP1D.i_Grid = GRID_UNIFORM;
    strcpy(IP1D.Grid_Type, "Grid Uniform");
    IP1D.Number_of_Cells_Idir = 5;
    IP1D.X_Shift = 0.0;
    IP1D.X_Scale = 1.0;
    IP1D.Number_of_Ghost_Cells = 3;
    strcat(IP1D.Output_File_Name, "Flow1D");
    IP1D.Number_of_SubGrid_Points = 4;
    IP1D.Reconstruction_Order = 3;
    IP1D.X_min = 1.0;
    IP1D.X_max = 4.0;
    IP1D.TestF = Test_Default1D;
    IP.ReconstructionMethod() = ENO ;

    // Set ComputationalDomain
    Domain1D.SetDomain(IP1D);
    Domain1D.SetInitialData(IP1D);

    // Set new values for the geometry

    Cell1D_NonUniform NewCell;

    NewCell.dx = 1.0;   
    // Cell 0
    NewCell.x = -4.0;
    Domain1D(3).SetCell(NewCell);
    Domain1D(3).CellSolution() = -240.0;
    // Cell 1
    NewCell.x = -1.0;
    Domain1D(4).SetCell(NewCell);
    Domain1D(4).CellSolution() = -30.0;
    // Cell 2
    NewCell.x = 1.0;
    Domain1D(5).SetCell(NewCell);
    Domain1D(5).CellSolution() = 0.0;
    // Cell 3
    NewCell.x = 2.0;
    Domain1D(6).SetCell(NewCell);
    Domain1D(6).CellSolution() = 0.0;
    // Cell 4
    NewCell.x = 5.0;
    Domain1D(7).SetCell(NewCell);
    Domain1D(7).CellSolution() = 12.0;

    //    Print(Domain1D);

    ensure("Order1", DivDifference(Domain1D,1,1,Domain1D.iStart()) == 70);
    ensure("Order1", DivDifference(Domain1D,1,1,Domain1D.iStart()+3) == 4);

    ensure("Order2", DivDifference(Domain1D,1,2,Domain1D.iStart()) == -11);
    ensure("Order2", DivDifference(Domain1D,1,2,Domain1D.iStart()+1) == -5);

    ensure("Order3", DivDifference(Domain1D,1,3,Domain1D.iStart()) == 1);
    ensure("Order3", DivDifference(Domain1D,1,3,Domain1D.iStart()+1) == 1);

    ensure("Order4", DivDifference(Domain1D,1,4,Domain1D.iStart()) == 0.0);

  }

  /*Test6: ENO_Reconstruction() */
  template<>
  template<>
  void Reconstruction_TestObject::test<6>()
  {
    // Set the domain

    // Set Path
    char path[256];
    strcpy(path,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(path,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcpy(IP.Output_File_Name,path);

    // Set Input Parameters
    IP1D.i_Grid = GRID_UNIFORM;
    strcpy(IP1D.Grid_Type, "Grid Uniform");
    IP1D.Number_of_Cells_Idir = 10;
    IP1D.X_Shift = 0.0;
    IP1D.X_Scale = 1.0;
    IP1D.Number_of_Ghost_Cells = 4;
    strcat(IP1D.Output_File_Name, "Flow1D");
    IP1D.Number_of_SubGrid_Points = 4;
    IP1D.Reconstruction_Order = 4;
    IP1D.X_min = 1.0;
    IP1D.X_max = 4.0;
    IP1D.TestF = Test_Example3;
    IP1D.ReconstructionMethod() = ENO ;

    // Set ComputationalDomain
    Domain1D.SetDomain(IP1D);
    Domain1D.SetInitialData(IP1D);

    Domain1D.ReconstructSolution(IP1D);
  }

#endif

}

namespace tut
{
  Reconstruction_TestGroup Reconstruction_Test("Reconstruction_Test");
}
