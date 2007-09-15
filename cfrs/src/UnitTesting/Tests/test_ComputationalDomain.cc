#include <tut.h>
#include "Reconstruction/ComputationalDomain.h"
#include "TestFunctions/TestFunctions_2D.h"
#include "CFD/Gaussian2DState.h"
#include "Grid/Grid2D/QuadrilateralGrid.h"
#include <limits>

using namespace std;

namespace tut
{
  struct Data_ComputationalDomain{
    // data types
    typedef ComputationalDomain<TwoD,Cell2D_Quad,Gaussian2D_pState> GaussianCompDomainType;
    typedef ComputationalDomain<TwoD,Cell2D_Quad,double> CompDomainType;
    typedef ComputationalDomain<OneD,Cell1D_NonUniform,double> CompDomainType1D;
    typedef ComputationalDomain<ThreeD,Cell3D,Gaussian2D_pState> Gaussian3DCompDomainType;

    // variable declaration
    CompDomainType1D Domain1D;
    CompDomainType A;
    GaussianCompDomainType GDom;
    Gaussian3DCompDomainType Domain3D;
    int NCx, NCy, NCz, Nghost;
    Grid2D_Quad_Block   **MeshBlk; // the mesh
    Grid2D_Quad_Block   Grid;
    Reconstruct2D_Input_Parameters IP; // input parameters
    Reconstruct1D_Input_Parameters IP1D;
  };

  typedef test_group<Data_ComputationalDomain> ComputationalDomain_TestGroup;
  typedef ComputationalDomain_TestGroup::object ComputationalDomain_TestObject;
}

namespace tut
{

  /*Test1: ComputationalDomain::SetDomain() */
  template<>
  template<>
  void ComputationalDomain_TestObject::test<1>()
  {
    // initial data
    NCx = 2; NCy = 1; NCz = 2;
    Nghost = 2;

    // check 1D domain
    Domain1D.SetDomain(NCx,Nghost);

    //    Print(Domain1D);
    //    Print(Domain1D(1));
    // Copy Constructor
    CompDomainType1D NewDomain1D = Domain1D;
    ensure("Copy Constructor", NewDomain1D == Domain1D);

        // Assignment operator
    CompDomainType1D NewNewDomain1D;
    NewNewDomain1D = Domain1D;
    ensure("Assignment operator", NewNewDomain1D == Domain1D);

    // check 2D domain
    A.SetDomain(NCx,NCy,Nghost);

    //    Print(A(1,1));

        // Copy Constructor
    CompDomainType B = A;
    ensure("Copy Constructor", B == A);

        // Assignment operator
    CompDomainType C;
    C = A;
    ensure("Assignment operator", C == A);

    //    check 3D domain
    Domain3D.SetDomain(NCx,NCy,NCz,Nghost);

        // Copy Constructor
    Gaussian3DCompDomainType NewDomain3D = Domain3D;
    ensure("Copy Constructor", NewDomain3D == Domain3D);

        // Assignment operator
    Gaussian3DCompDomainType NewNewDomain3D;
    NewNewDomain3D = Domain3D;
    ensure("Assignment operator", NewNewDomain3D == Domain3D);

    //    Print(Domain3D(1,1,1));
    
    fail("This test is dezactivated!");

  }

  /*Test2: ComputationalDomain::*/
  template<>
  template<>
  void ComputationalDomain_TestObject::test<2>()
  {

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
    IP.Box_Width = 2.0;
    IP.Box_Height = 3.0;
    IP.Number_of_Cells_Idir = 10;
    IP.Number_of_Cells_Jdir = 5;
    Vector2D Shift(2.0,3.0);
    IP.X_Shift = Shift;
    IP.X_Scale = 1.5;
    IP.X_Rotate = 0.0;
    IP.Number_of_Ghost_Cells = 2;
    IP.Number_of_SubGrid_Points_Idir = 3;
    IP.Number_of_SubGrid_Points_Jdir = 3;
    IP.Reconstruction_Order = 3;
    strcat(IP.Output_File_Name, "RectangularBox");

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

    //    cout << Grid << endl;

    // Set ComputationalDomain
    A.SetDomain(Grid,IP);

    // Plot the geometry
    // Set Path
    ofstream output_file;
    char *output_file_name_ptr;
    char output_file_name[256];
    output_file_name_ptr = output_file_name;

    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"ComputationalDomain_Nodes_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    // OutputNodesTecplot
    A.OutputNodesTecplot(output_file);
    output_file.close();

    // Output Mesh
    Output_Nodes_Tecplot(MeshBlk, IP);
  }

#if 0
  /*Test3: ComputationalDomain:: */
  template<>
  template<>
  void ComputationalDomain_TestObject::test<3>()
  {

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
    strcat(IP.Output_File_Name, "OriginalRinglebFlow");

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
    A.SetDomain(Grid,IP);

    ensure("NumberOfDerivatives",A.NumberOfTaylorDerivatives() == 10);

    // Plot the geometry
    // Set Path
    ofstream output_file;
    char *output_file_name_ptr;
    char output_file_name[256];
    output_file_name_ptr = output_file_name;

    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"ComputationalDomainRinglebFlow.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    A.OutputNodesTecplot(output_file);
    output_file.close();

    // OutputSolutionNodesTecplot
    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"CompDomainRF_SolutionAtNodes_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    A.OutputSolutionNodesTecplot(output_file);
    output_file.close();

    // OutputFullSolutionNodesTecplot
    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"CompDomainRF_FullSolutionAtNodes_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";
    
    A.OutputFullSolutionNodesTecplot(output_file);
    output_file.close();

    // OutputSolutionCellTecplot
    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"CompDomainRF_SolutionCell_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    A.OutputSolutionCellTecplot(output_file);
    output_file.close();

    // OutputFullSolutionCellTecplot
    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"CompDomainRF_FullSolutionCell_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    A.OutputFullSolutionCellTecplot(output_file);
    output_file.close();

    //    Output Mesh
    Output_Nodes_Tecplot(MeshBlk, IP);

  }

  /*Test4: ComputationalDomain::*/
  template<>
  template<>
  void ComputationalDomain_TestObject::test<4>()
  {
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
    strcat(IP.Output_File_Name, "GaussianRinglebFlow");

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
    GDom.SetDomain(Grid,IP);

    HeaderData Gaussian;
    Gaussian.add("Var1"); Gaussian.add("Var2"); Gaussian.add("Var3");
    Gaussian.add("Var4"); Gaussian.add("Var5"); Gaussian.add("Var6");
    Gaussian.add("Var7"); Gaussian.add("Var8");

    GDom.DefineHeader(Gaussian);

    // Plot the geometry
    // Set Path
    ofstream output_file;
    char *output_file_name_ptr;
    char output_file_name[256];
    output_file_name_ptr = output_file_name;

    // OutputSolutionNodesTecplot
    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"GaussianRinblebFlow_SolutionAtNodes_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    GDom.OutputSolutionNodesTecplot(output_file);
    output_file.close();

    // OutputFullSolutionNodesTecplot
    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"GaussianRinblebFlow_FullSolutionAtNodes_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";
    
    GDom.OutputFullSolutionNodesTecplot(output_file);
    output_file.close();

    // OutputSolutionCellTecplot
    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"GaussianRinblebFlow_SolutionCell_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    GDom.OutputSolutionCellTecplot(output_file);
    output_file.close();

    // OutputFullSolutionCellTecplot
    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"GaussianRinblebFlow_FullSolutionCell_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    GDom.OutputFullSolutionCellTecplot(output_file);
    output_file.close();
  }

  /*Test5: ComputationalDomain::SetMaxDeltaSolutionOverDomain()*/
  template<>
  template<>
  void ComputationalDomain_TestObject::test<5>()
  {
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
    IP.Number_of_Cells_Idir = 4;
    IP.Number_of_Cells_Jdir = 3;
    Vector2D Shift(0.0,0.0);
    IP.X_Shift = Shift;
    IP.X_Scale = 1.5;
    IP.X_Rotate = 0.0;
    IP.Number_of_Ghost_Cells = 3;
    strcat(IP.Output_File_Name, "RinglebFlow");
    IP.Number_of_SubGrid_Points_Idir = 4;
    IP.Number_of_SubGrid_Points_Jdir = 4;
    IP.Reconstruction_Order = 3;
    strcat(IP.Output_File_Name, "GaussianRinglebFlow");

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
    GDom.SetDomain(Grid,IP);
    
    // Initialize Gaussian State
    
    Gaussian2D_pState GState;

    for(int i = GDom.iStart(); i<= GDom.iEnd(); ++i )
      for(int j = GDom.jStart(); j<= GDom.jEnd(); ++j ){
	for(int param=1; param<=8; ++param){
	  GState[param] = 3.14*param + i+j;
	}
	GDom(i,j).CellSolution() = GState;
      }


    // Compute the MaxDeltaSolutionOverDomain

    GDom.SetMaxDeltaSolutionOverDomain();

    // Result
    vector<double> Result;
    Result.reserve(8);

    Result.push_back(14.14); Result.push_back(17.28); Result.push_back(20.42); Result.push_back(23.56);
    Result.push_back(26.7); Result.push_back(29.84); Result.push_back(32.98); Result.push_back(36.12);

    for(int i=1; i<=8; ++i){
      ensure("MaxDeltaSolutionOverDomain check result", (GDom.MaxDeltaSolutionDomain(i-1) - Result[i-1])< 1.0e-14);
    }
  }

  /*Test6: ComputationalDomain1D::SetDomain()*/
  template<>
  template<>
  void ComputationalDomain_TestObject::test<6>()
  {
    //Domain1D    
    
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

    // Set ComputationalDomain
    Domain1D.SetDomain(IP1D);
    Domain1D.SetInitialData(IP1D);

    //    Print(Domain1D);

    // Plot the geometry
    // Set Path
    ofstream output_file;
    char *output_file_name_ptr;
    char output_file_name[256];
    output_file_name_ptr = output_file_name;

    // OutputSolutionNodesTecplot
    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"Domain1D_SolutionAtNodes_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    Domain1D.OutputSolutionNodesTecplot(output_file);
    output_file.close();

    // OutputFullSolutionNodesTecplot
    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"Domain1D_FullSolutionAtNodes_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";
    
    Domain1D.OutputFullSolutionNodesTecplot(output_file);
    output_file.close();

    // OutputSolutionCellTecplot
    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"Domain1D_SolutionCell_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    Domain1D.OutputSolutionCellTecplot(output_file);
    output_file.close();

    // OutputFullSolutionCellTecplot
    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"Domain1D_FullSolutionCell_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    Domain1D.OutputFullSolutionCellTecplot(output_file);
    output_file.close();

  }

  /*Test7: Reflection */
  template<>
  template<>
  void ComputationalDomain_TestObject::test<7>()
  {
    //Domain1D    
    
    // Set Input Parameters
    IP1D.i_Grid = GRID_UNIFORM;
    strcpy(IP1D.Grid_Type, "Grid Uniform");
    IP1D.Number_of_Cells_Idir = 7;
    IP1D.X_Shift = 0.0;
    IP1D.X_Scale = 1.0;
    IP1D.Number_of_Ghost_Cells = 3;
    strcat(IP1D.Output_File_Name, "Flow1D");
    IP1D.Number_of_SubGrid_Points = 4;
    IP1D.Reconstruction_Order = 3;
    IP1D.X_min = 1.0;
    IP1D.X_max = 4.0;
    IP1D.TestF = Test_Default1D;

    // Set ComputationalDomain
    Domain1D.SetDomain(IP1D);
    Domain1D.SetInitialData(IP1D);

    // Plot the geometry
    // Set Path
    ofstream output_file;
    char *output_file_name_ptr;
    char output_file_name[256];
    output_file_name_ptr = output_file_name;

    // OutputSolutionNodesTecplot
    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"Domain1D_SolutionAtNodes_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    Domain1D.OutputSolutionNodesTecplot(output_file);
    output_file.close();

    // OutputFullSolutionNodesTecplot
    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"Domain1D_FullSolutionAtNodes_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";
    
    Domain1D.OutputFullSolutionNodesTecplot(output_file);
    output_file.close();

    // OutputSolutionCellTecplot
    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"Domain1D_SolutionCell_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    Domain1D.OutputSolutionCellTecplot(output_file);
    output_file.close();

    // OutputFullSolutionCellTecplot
    strcpy(output_file_name,"/nfs/carv/d1/people/lucian/Documents/Fluids/");
    strcat(output_file_name,"Thesis/Reconstruction/src/TUT_Tests/TestResults/");
    strcat(output_file_name,"Domain1D_FullSolutionCell_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    Domain1D.OutputFullSolutionCellTecplot(output_file);
    output_file.close();
  }

  /*Test8: ComputationalDomain1D::ComputeMultipleCorrelationCoefficient */
  template<>
  template<>
  void ComputationalDomain_TestObject::test<8>()
  {
    // Set Input Parameters
    IP1D.i_Grid = GRID_UNIFORM;
    strcpy(IP1D.Grid_Type, "Grid Uniform");
    IP1D.Number_of_Cells_Idir = 32;
    IP1D.X_Shift = 0.0;
    IP1D.X_Scale = 1.0;
    IP1D.Number_of_Ghost_Cells = 3;
    strcat(IP1D.Output_File_Name, "Flow1D_MCC");
    IP1D.Number_of_SubGrid_Points = 4;
    IP1D.Reconstruction_Order = 3;
    IP1D.X_min = -1.0;
    IP1D.X_max = 4.0;
    // IP1D.TestF = Test_Default1D;
    IP1D.TestF = Test_Example5;
    //    IP1D.TestF = Test_Example3;
    // IP1D.TestF = Test_Example2;
    IP1D.Method = DD_ENO;

    // Set ComputationalDomain
    Domain1D.SetDomain(IP1D);
    Domain1D.SetInitialData(IP1D);

    // Reconstruct solution
    Domain1D.ReconstructSolution(IP1D);

    // Update subgrid solution
    Domain1D.UpdateSubgridSolution();


    // Compute the Multiple-Correlation Coefficient
    Domain1D.ComputeMultipleCorrelationCoefficient();
  }

  /*Test9: ComputationalDomain::*/
  template<>
  template<>
  void ComputationalDomain_TestObject::test<9>()
  {
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
    IP.Number_of_Cells_Idir = 64;
    IP.Number_of_Cells_Jdir = 64;
    Vector2D Shift(0.0,0.0);
    IP.X_Shift = Shift;
    IP.X_Scale = 1.5;
    IP.X_Rotate = 0.0;
    IP.Number_of_Ghost_Cells = 3;
    strcat(IP.Output_File_Name, "RinglebFlow");
    IP.Number_of_SubGrid_Points_Idir = 4;
    IP.Number_of_SubGrid_Points_Jdir = 4;
    IP.Reconstruction_Order = 3;
    strcat(IP.Output_File_Name, "RinglebFlow");
    IP.Method = DD_ENO;
    //    IP.TestF = Test_Default2D;
    IP.TestF = Test_Example5;

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
    A.SetDomain(Grid,IP);

    A.SetInitialData(IP);

    // Reconstruct solution
    A.ReconstructSolution(IP);

    // Update subgrid solution
    A.UpdateSubgridSolution();

#if 0
    vector<int> i_index,j_index;
    int CellsInOneDirection;
    
    CellsInOneDirection = 3 + 2*(A(0,0).CellRings() - 1);
    i_index.reserve(CellsInOneDirection);
    j_index.reserve(CellsInOneDirection);
    
    for (int i=A.iStart(); i<=A.iEnd(); ++i)
      for (int j=A.jStart(); j<=A.jEnd(); ++j){
    	// Make Stencil			
    	MakeReconstructionStencil(A(0,0).CellRings(),i,j,i_index,j_index);
	
    	// estimate the correlation coefficient for the current cell
    	MultipleCorrelationCoefficient(A,i_index,j_index,i,j);
      }
#endif

    for (int i=A.iStart(); i<=A.iEnd(); ++i)
      for (int j=A.jStart(); j<=A.jEnd(); ++j){
	if ( A(i,j).CellMCC() < 430 ){
	  // piecewise constant
	  A(i,j).CellDeriv(0,true,true,true) = A(i,j).CellSolution();
	  for (int TD = 1; TD< A(i,j).NumberOfTaylorDerivatives(); ++TD){
	    A(i,j).CellDeriv(TD,true,true,true) = 0.0;
	  }
	}
      }

    // Update subgrid solution
    A.UpdateSubgridSolution();

    // Plot the geometry
    // Set Path
    ofstream output_file;
    char *output_file_name_ptr;
    char output_file_name[256];
    output_file_name_ptr = output_file_name;

    // OutputSolutionNodesTecplot
    strcpy(output_file_name,"TUT_Tests/TestResults/");
    strcat(output_file_name,"CompDomainRF_SolutionAtNodes_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    A.OutputSolutionNodesTecplot(output_file);
    output_file.close();

    // OutputFullSolutionNodesTecplot
    strcpy(output_file_name,"TUT_Tests/TestResults/");
    strcat(output_file_name,"CompDomainRF_FullSolutionAtNodes_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";
    
    A.OutputFullSolutionNodesTecplot(output_file);
    output_file.close();

    // OutputSolutionCellTecplot
    strcpy(output_file_name,"TUT_Tests/TestResults/");
    strcat(output_file_name,"CompDomainRF_SolutionCell_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    A.OutputSolutionCellTecplot(output_file);
    output_file.close();

    // OutputFullSolutionCellTecplot
    strcpy(output_file_name,"TUT_Tests/TestResults/");
    strcat(output_file_name,"CompDomainRF_FullSolutionCell_Tecplot.dat");
    output_file.open(output_file_name_ptr, ios::out);
    if (output_file.bad()) std::cout << "Problems\n";

    A.OutputFullSolutionCellTecplot(output_file);
    output_file.close();


    //    Output Mesh
    Output_Nodes_Tecplot(MeshBlk, IP);
  }
#endif

}

namespace tut
{
  ComputationalDomain_TestGroup ComputationalDomain_Test("ComputationalDomain_Test");
}
