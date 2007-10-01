/*!\file test_Cell2D.cc
   \brief Regression tests for Cell2D class.*/

/* Include required C++ libraries. */
// None

/* Using std namespace functions */

/* Include CFFC header files */
#include "TestData.h"
#include "Grid/Grid2D/Cell2D.h"

namespace tut

{
  // Data used for testing
  class Data_Cell2DCartesian: public TestData{
  public:
    Cell2D_Cartesian VarA;
    Cell2D_Cartesian VarB;
    double LocationX, LocationY, SizeX, SizeY;

    Data_Cell2DCartesian(void){
      set_test_suite_path("Grid/Grid2D/UnitTests");
    }
  };

  typedef test_group<Data_Cell2DCartesian> Cell2DCartesian_TestGroup;
  typedef Cell2DCartesian_TestGroup::object Cell2DCartesian_TestObject;

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
    ensure_equals("Cell2DCartesian: Copy Constructor", TestCell, VarA );
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
    ensure_equals("Cell2DCartesian: Assignment operator", TestCell, VarA );
  }

  /* Test3: */
  template<>
  template<>
  void Cell2DCartesian_TestObject::test<3>()
  {

    set_test_name("operator !=");

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

tut::Cell2DCartesian_TestGroup Cell2DCartesian_Test("Class:Cell2D_Cartesian");

