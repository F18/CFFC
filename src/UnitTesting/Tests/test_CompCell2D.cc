#include "../tut.h"
#include "../CompCell2D.h"
#include "../TestFunctions_2D.h"


using namespace std;

namespace tut

{
  // Data used for testing
  struct Data_CompCell2D{
//     CompCell2D_Cartesian TestCell1;
//     CompCell2D_Cartesian TestCell2;
//     FunctionType2D ExactSolution;
//     int Nx, Ny, OrderOfReconstruction;
//     double L1Norm;
//     double MatlabResult;
//     double tol;
//     double EpsMachine;
//     double RelErr;
 };

   /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   */
  typedef test_group<Data_CompCell2D> CompCell2D_TestGroup;
  typedef CompCell2D_TestGroup::object  CompCell2D_TestObject;
}

/**********************************************
 *                THE TESTS                   *
 *********************************************/

namespace tut
{

  /*Test1: Copy Constructor */
  template<>
  template<>
  void CompCell2D_TestObject::test<1>()
  {
//     // set number of points in the subdomain
//     Nx = 5; Ny = 10; OrderOfReconstruction = 3;

//     // create subdomain and 
//     TestCell1.allocate(Nx,Ny,OrderOfReconstruction);
//     CompCell2D_Cartesian NewCell = TestCell1;
//     ensure("Copy Constructor", NewCell == TestCell1);
  }

  /*Test2: Assignment Operator */
  template<>
  template<>
  void CompCell2D_TestObject::test<2>()
  {
    // set number of points in the subdomain
//     Nx = 5; Ny = 10; OrderOfReconstruction = 3;

//     // create subdomain and 
//     TestCell1.allocate(Nx,Ny,OrderOfReconstruction);
//     CompCell2D_Cartesian NewCell;
//     NewCell = TestCell1;
//     ensure("Assignment operator", NewCell == TestCell1);
  }



}

namespace tut
{
  CompCell2D_TestGroup CompCell2D_Test("CompCell2D_Test");
}
