// Test file for class Cell1D_Uniform
#include <tut.h>
#include "Grid/Grid2D/Cell2D.h"

namespace tut

{
  // Data used for testing
  struct Data_Cell2DCartesian{
    Cell2D_Cartesian VarA;
    Cell2D_Cartesian VarB;
    double LocationX, LocationY, SizeX, SizeY;
  };

   /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   */
  typedef test_group<Data_Cell2DCartesian> Cell2DCartesian_TestGroup;
  typedef Cell2DCartesian_TestGroup::object Cell2DCartesian_TestObject;
}

/**********************************************
 *                THE TESTS                   *
 *********************************************/

namespace tut
{

  /*Test1: Copy constructor */
  template<>
  template<>
  void Cell2DCartesian_TestObject::test<1>()
  {
    LocationX = 3.45;
    LocationY = 8.543434;
    SizeX = 23.3434;
    SizeY = 0.2323;
    VarA.setloc(LocationX, LocationY);
    VarA.setsize(SizeX, SizeY);
    Cell2D_Cartesian TestCell = VarA;
    ensure("Cell2DCartesian: Copy Constructor", TestCell == VarA );
  }

  /* Test2: Assignment operator */
  template<>
  template<>
  void Cell2DCartesian_TestObject::test<2>()
  {
    LocationX = 3.45;
    LocationY = 8.543434;
    SizeX = 23.3434;
    SizeY = 0.2323;
    VarA.setloc(LocationX, LocationY);
    VarA.setsize(SizeX, SizeY);
    Cell2D_Cartesian TestCell;

    TestCell = VarA;
    ensure("Cell2DCartesian: Assignment operator", TestCell == VarA );
  }

  /* Test3: Different*/
  template<>
  template<>
  void Cell2DCartesian_TestObject::test<3>()
  {
    LocationX = 3.45;
    LocationY = 8.543434;
    SizeX = 23.3434;
    SizeY = 0.2323;
    VarA.setloc(LocationX, LocationY);
    VarA.setsize(SizeX, SizeY);
    Cell2D_Cartesian TestCell;

    TestCell.setloc(LocationX + 0.34, LocationY/0.23);
    TestCell.setsize(SizeX,SizeY);

    ensure("Cell2DCartesian: != operator", TestCell != VarA);
  }

}

namespace tut
{
  Cell2DCartesian_TestGroup Cell2DCartesian_Test("Cell2DCartesian_Test");
}
