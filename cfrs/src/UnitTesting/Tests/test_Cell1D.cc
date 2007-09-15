// Test file for class Cell1D_Uniform
#include <tut.h>
#include "Grid/Grid1D/Cell1D.h"

namespace tut

{
  // Data used for testing
  struct Data_UniformCell1D{
    Cell1D_Uniform VarA;
    Cell1D_Uniform VarB;
  };

   /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (UniformCell_TestGroup) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_UniformCell1D> UniformCell1D_TestGroup;
  typedef UniformCell1D_TestGroup::object UniformCell1D_TestObject;
}

/**********************************************
 *                THE TESTS                   *
 *********************************************/

namespace tut
{

  Cell1D_Uniform VarC(11.2); 

  template<>
  template<>
  void UniformCell1D_TestObject::test<1>()
  {
    VarB = VarC;
    ensure ("Cell Constructor with Length",VarB.x==11.2);
  }

  template<>
  template<>
  void UniformCell1D_TestObject::test<2>()
  {
    VarB = VarC;
    ensure ("equality of two cells",!(VarA==VarB));
  }

  template<>
  template<>
  void UniformCell1D_TestObject::test<3>()
  {
    ensure("Cell Constructor",VarA.x == 1);
  }
}

namespace tut
{
  UniformCell1D_TestGroup Cell1D_Uniform_Test("Cell1D_Uniform_Test");
}
