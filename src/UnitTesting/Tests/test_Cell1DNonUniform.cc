#include <tut.h>
#include "Grid/Grid1D/Cell1D.h"

namespace tut
{
  struct shared_data_NonUniformCell{ };

  typedef test_group<shared_data_NonUniformCell> testgroup;
  typedef testgroup::object testobject;
}

namespace tut
{

  Cell1D_NonUniform CellVarA;
  Cell1D_NonUniform CellVarB(1.2, 3.0);

  template<>
  template<>
  void testobject::test<1>()
  {
    ensure ("Cell Constructor",CellVarB.x==1.2 && CellVarB.dx==3.0);
  }

}

namespace
{
  tut::testgroup shared_data_testgroup("Cell1D_NonUniform_Test");
}
