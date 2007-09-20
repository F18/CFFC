/*!\file test_ComputationalCell.cc
   \brief Regression tests for templated class ComputationalCell. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */

/* Include CFFC header files */
#include "TestData.h"
#include "include/ComputationalCell.h"
#include "TestFunctions/TestFunctions_2D.h"
#include "CFD/Gaussian2DState.h"
#include "../../../../src_2D/HighOrderReconstruction/ReconstructionHelpers.h"


using namespace std;

namespace tut
{

  // **********************************************
  //                 TEST SUITE: 
  // **********************************************
  // Data used for testing
  // **********************
  struct Data_ComputationalCell{

    // data type definition
    typedef ComputationalCell<OneD,Cell1D_NonUniform,double> ComputationalCell1DType;
    typedef ComputationalCell<TwoD,Cell2D_Cartesian,double> ComputationalCell2DType;
    typedef ComputationalCell<TwoD,Cell2D_Cartesian,Gaussian2D_pState> GaussianComputationalCell2DType;
    typedef ComputationalCell<ThreeD,Cell3D,double> ComputationalCell3DType;

    ComputationalCell1DType CompCell1D;
    ComputationalCell2DType CompCell2D;
    GaussianComputationalCell2DType G_CompCell2D;
    ComputationalCell3DType CompCell3D;
    int Nx,Ny,Nz;
    int RecOrder;
    double tol;
    int digits;
  };

  typedef test_group<Data_ComputationalCell> ComputationalCell_TestGroup;
  typedef ComputationalCell_TestGroup::object ComputationalCell_TestObject;


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

  /*Test1: ComputationalCell::Constructor & Copy constructor & Assignment operator */
  template<>
  template<>
  void ComputationalCell_TestObject::test<1>()
  {

    GaussianComputationalCell2DType G_CompCell2D_Test = G_CompCell2D;
    ensure("Copy constructor", G_CompCell2D_Test == G_CompCell2D);

    GaussianComputationalCell2DType G_CompCell2D_Test2;
    G_CompCell2D_Test2 = G_CompCell2D;
    ensure("Assignment operator", G_CompCell2D_Test2 == G_CompCell2D);

    // Change the parameters and verify that the computational cells are different
    Nx = 3; Ny = 4; Nz = 2;
    RecOrder = 3;
    std::vector<int> SubgridPoints;
    SubgridPoints.push_back(Nx);
    SubgridPoints.push_back(Ny);
    SubgridPoints.push_back(Nz);
   
    double Xloc, Yloc, Xsize, Ysize;

    Cell2D_Cartesian Cell;
    Xloc = 1.2; Yloc = 3.4;
    Cell.setloc(Xloc,Yloc);
    G_CompCell2D_Test.SetCell(Cell,SubgridPoints,RecOrder);
    ensure("Equal operator", G_CompCell2D_Test != G_CompCell2D);
  }

  /*Test3: ComputationalCell::SetCell() */
  template<>
  template<>
  void ComputationalCell_TestObject::test<3>()
  {
    Nx = 3; Ny = 4; Nz = 2;
    RecOrder = 2;
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);

    std::vector<int> SubgridPoints;
    SubgridPoints.push_back(Nx);
    SubgridPoints.push_back(Ny);
    SubgridPoints.push_back(Nz);
    
    double Xloc, Yloc, Xsize, Ysize;

    ComputationalCell<TwoD,Cell2D_Quad,Gaussian2D_pState> G_CompCell2D_Test;

    G_CompCell2D_Test.SetCell(SubgridPoints,RecOrder);

    Vector2D Node1, Node2, Node3, Node4;

    Node1.x = 1.2; Node1.y = 0.3;
    Node2.x = 1.5; Node2.y = 1.3;
    Node3.x = 2.1; Node3.y = 1.5;
    Node4.x = 2.5; Node4.y = 1.2;

    G_CompCell2D_Test.CellGeometry().setnodes(Node1,Node2,Node3,Node4);
    G_CompCell2D_Test.SetSubGridGeometry();

#ifdef __Use_Iterator__ 
    TaylorDerivativesContainer<TwoD,Gaussian2D_pState>::iterator Iter = G_CompCell2D_Test.CellDeriv();  
    for(Iter=G_CompCell2D_Test.BeginD(); Iter!=G_CompCell2D_Test.EndD(); Iter++){
      Iter->D()[1] = 2;
      Iter->D()[2] = 1.2;
      Iter->D()[3] = 2.3;
      Iter->D()[4] = 0.2323;
      Iter->D()[5] = 0.2;
      Iter->D()[6] = 0.1;
      Iter->D()[7] = 0.2;
      Iter->D()[8] = 0.02;
    }
#else
    for(int i=G_CompCell2D_Test.CellDeriv().FirstElem(); i<=G_CompCell2D_Test.CellDeriv().LastElem(); ++i){

      G_CompCell2D_Test.CellDeriv(i,true,true,true).D()[1] = 2;
      G_CompCell2D_Test.CellDeriv(i,true,true,true).D()[2] = 1,2;
      G_CompCell2D_Test.CellDeriv(i,true,true,true).D()[3] = 2.3;
      G_CompCell2D_Test.CellDeriv(i,true,true,true).D()[4] = 0.2323;
      G_CompCell2D_Test.CellDeriv(i,true,true,true).D()[5] = 0.2;
      G_CompCell2D_Test.CellDeriv(i,true,true,true).D()[6] = 0.1;
      G_CompCell2D_Test.CellDeriv(i,true,true,true).D()[7] = 0.2;
      G_CompCell2D_Test.CellDeriv(i,true,true,true).D()[8] = 0.02;
    }
#endif

    ensure("Derivatives", G_CompCell2D_Test.CellDeriv(0,0)[1] == 2);
    ensure("Derivatives", G_CompCell2D_Test.CellDeriv(0,0)[8] == 0.02);

#if 0
    Gaussian2D_pState G1,G2,G3;
    double deltaX, deltaY;
    deltaX = 0.1; deltaY = 0.3;

    // G1 = G_CompCell2D_Test.CellDeriv().ComputeSolutionFor(deltaX,deltaY);

    //    ensure("SolutionaAt", (G1[1] -3.060e00) < tol );
    G2 = G_CompCell2D_Test.SolutionAt(deltaX,deltaY);


    ensure("SolutionAt", (G2[1]- G1[1]) < tol);

    Node2D node(1.925,1.375);
    G3 = G_CompCell2D_Test.SolutionAt(node);
    ensure("SolutionAt", (G3[1] - G1[1]) < tol);

    G_CompCell2D_Test.UpdateSubgridSolution();

    for (int i=1; i<=G_CompCell2D_Test.NumberOfVariables; ++i){
      G_CompCell2D_Test.CellSolution()[i] = i;
    }

    // Obs. There is no function to access the subgrid and therefore I cannot use the "ensure" method
    //      Visual inspection of the coordinates was used instead 
    // cout <<  "The Computational Cell" << endl << G_CompCell2D_Test << endl;

    for (int i=1; i<=G_CompCell2D_Test.NumberOfVariables; ++i){
      ensure("CellSolution(i)", G_CompCell2D_Test.CellSolution(i) == i);
    }
#endif

  }

  /*Test4: ComputationalCell::GenerateCell( )*/
  template<>
  template<>
  void ComputationalCell_TestObject::test<4>()
  {
    Nx = 3; Ny = 4; Nz = 2;
    RecOrder = 3;
    std::vector<int> SubgridPoints;
    SubgridPoints.push_back(Nx);
    SubgridPoints.push_back(Ny);
    SubgridPoints.push_back(Nz);

    double Xloc, Yloc, Xsize, Ysize;
    CompCell1D.GenerateCell(Nx,RecOrder);

    Cell2D_Cartesian Cell;
    Xloc = 1.2; Yloc = 3.4; Xsize = 0.43; Ysize = 0.54;
    Cell.setloc(Xloc,Yloc);
    Cell.setsize(Xsize,Ysize);
    GaussianComputationalCell2DType G_CompCell2D_Test;
    G_CompCell2D_Test.SetCell(Cell,SubgridPoints,RecOrder);
    CompCell3D.GenerateCell(Nx,Ny,Nz,RecOrder);
  }

  /*Test5: ComputationalCell::SetCell( )*/
  template<>
  template<>
  void ComputationalCell_TestObject::test<5>()
  {
    Nx = 3; Ny = 4; Nz = 2;
    RecOrder = 3;

    // Construct the cell
    Cell1D_NonUniform Cell;
    Cell.x = 2.3;
    Cell.dx = 0.45;
    std::vector<int> SubgridPoints;
    SubgridPoints.push_back(Nx);
    SubgridPoints.push_back(Ny);
    SubgridPoints.push_back(Nz);

    ComputationalCell<OneD,Cell1D_NonUniform,Gaussian2D_pState> CompCell1D_Test;
    CompCell1D_Test.SetCell(Cell,SubgridPoints,RecOrder);

    ComputationalCell<OneD,Cell1D_NonUniform,Gaussian2D_pState> CompCell1D_Test2;
    CompCell1D_Test2 = CompCell1D_Test;
    ensure("Equal 1D", CompCell1D_Test2 == CompCell1D_Test);

    CompCell1D.SetCell(Cell,SubgridPoints,RecOrder);

    ensure("Geometry", CompCell1D.CellGeometry() == Cell);
    ensure("CellCenter", CompCell1D.CellCenter() == Cell.x);
    ensure("CellRings", CompCell1D.CellRings() == 2);
    ensure("RecOrder", CompCell1D.CellRecOrder() == RecOrder);
    ensure("Final Order", CompCell1D.CellFOrder() == 3);
    ensure("Solution", CompCell1D.CellSolution() == 0.0);
    ensure("Derivatives", CompCell1D.CellDeriv(3) == 0.0);
    ensure("Error", CompCell1D.CellErrorL1() == 0.0);

    int p1 = 3;
    double Deriv = 10.0;
    CompCell1D.CellDeriv(p1) = Deriv;
    CompCell1D.CellDeriv(0) = 1.0;
    CompCell1D.CellDeriv(1) = 2.0;
    CompCell1D.CellDeriv(2) = 3.0;

    ensure("Derivatives", CompCell1D.CellDeriv(p1) == Deriv);

    double DeltaPoint = +Cell.dx/2;
    double Point = Cell.x + Cell.dx/2;
    ensure("SolutionAt", CompCell1D.SolutionAt(DeltaPoint) - 1.71578125 < 1.0e-8);
    ensure("SolutionAtCoordinate", CompCell1D.SolutionAtCoordinates(Point) - 1.71578125 < 1.0e-8);

  }

  /*Test6: ComputationalCell::SetCell() */
  template<>
  template<>
  void ComputationalCell_TestObject::test<6>()
  {
    Nx = 4; Ny = 6;
    RecOrder = 3;
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);

    std::vector<int> SubgridPoints;
    SubgridPoints.push_back(Nx);
    SubgridPoints.push_back(Ny);
    
    double Xloc, Yloc, Xsize, Ysize;
    Xloc = 5.3; Yloc = 6.7;
    Xsize = 0.0234; Ysize = 0.000123;

    GaussianComputationalCell2DType G_CompCell2D_Test;

    G_CompCell2D_Test.SetCell(SubgridPoints,RecOrder);

    G_CompCell2D_Test.CellGeometry().setloc(Xloc,Yloc);
    G_CompCell2D_Test.CellGeometry().setsize(Xsize,Ysize);
    G_CompCell2D_Test.ComputeGeometricCoefficients();
    G_CompCell2D_Test.SetSubGridGeometry();
  
#ifdef __Use_Iterator__
    TaylorDerivativesContainer<TwoD,Gaussian2D_pState>::iterator Iter = G_CompCell2D_Test.CellDeriv();
    for(Iter=G_CompCell2D_Test.BeginD(); Iter!=G_CompCell2D_Test.EndD(); Iter++){
      Iter->D()[1] = 2;
      Iter->D()[2] = 1.2;
      Iter->D()[3] = 2.3;
      Iter->D()[4] = 0.2323;
      Iter->D()[5] = 0.2;
      Iter->D()[6] = 0.1;
      Iter->D()[7] = 0.2;
      Iter->D()[8] = 0.02;
    }
#else
    for(int i=G_CompCell2D_Test.CellDeriv().FirstElem(); i<=G_CompCell2D_Test.CellDeriv().LastElem(); ++i){

      G_CompCell2D_Test.CellDeriv(i,true,true,true).D()[1] = 2;
      G_CompCell2D_Test.CellDeriv(i,true,true,true).D()[2] = 1,2;
      G_CompCell2D_Test.CellDeriv(i,true,true,true).D()[3] = 2.3;
      G_CompCell2D_Test.CellDeriv(i,true,true,true).D()[4] = 0.2323;
      G_CompCell2D_Test.CellDeriv(i,true,true,true).D()[5] = 0.2;
      G_CompCell2D_Test.CellDeriv(i,true,true,true).D()[6] = 0.1;
      G_CompCell2D_Test.CellDeriv(i,true,true,true).D()[7] = 0.2;
      G_CompCell2D_Test.CellDeriv(i,true,true,true).D()[8] = 0.02;
    }
#endif

    G_CompCell2D_Test.UpdateSubgridSolution();
    // Obs. There is no function to access the subgrid and therefore I cannot use the "ensure" method
    //      Visual inspection of the coordinates was used instead 
    //    cout <<  "The Computational Cell" << endl << G_CompCell2D_Test << endl;

    int rings, RecOrder;
    rings = G_CompCell2D_Test.CellRings();
    RecOrder = G_CompCell2D_Test.CellRecOrder();

    // Check the geometric coefficients
    vector<double> CartesianCoeff;
    int Memb = (RecOrder+1)*(RecOrder+2)/2;
    CartesianCoeff.reserve(Memb);

    int p2_limit = RecOrder;

    for (int p1=0; p1<=RecOrder; ++p1, p2_limit--){
      for (int p2=0; p2<=p2_limit; ++p2){
	CartesianCoeff.push_back(GeomCoeffCartesian(p1,p2,Xsize,Ysize,0,0));
      }
    }

  }

}

tut::ComputationalCell_TestGroup ComputationalCell_Test("Template Class:ComputationalCell");

