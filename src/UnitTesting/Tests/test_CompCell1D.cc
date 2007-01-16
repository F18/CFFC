#include <tut.h>
#include "Grid/Grid1D/Cell1D.h"
#include "Reconstruction/Reconstruction1D/CompCell1D.h"
#include "TestFunctions/TestFunctions.h"


using namespace std;

namespace tut

{
  // Data used for testing
  struct Data_CompCell1D{
    CompCell1D_NonUniform TestCell;
    FunctionType1D ExactSolution;
    double L1Norm;
    double MatlabResult;
    double tol;
    double EpsMachine;
    double RelErr;
 };

   /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   */
  typedef test_group<Data_CompCell1D> CompCell1D_TestGroup;
  typedef CompCell1D_TestGroup::object  CompCell1D_TestObject;
}

/**********************************************
 *                THE TESTS                   *
 *********************************************/

namespace tut
{
  CompCell1D_TestGroup CompCell1D_Test("CompCell1D_Test");
}
