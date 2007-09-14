// Test file for class Cell1D_Uniform
#include <tut.h>
#include "Common/SourceRevisionData.h"

namespace tut

{
  // Data used for testing
  struct Data_SourceCode{ };

   /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (UniformCell_TestGroup) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_SourceCode> SourceCode_TestGroup;
  typedef SourceCode_TestGroup::object SourceCode_TestObject;
}

/**********************************************
 *                THE TESTS                   *
 *********************************************/

namespace tut
{

  template<>
  template<>
  void SourceCode_TestObject::test<1>()
  {
    std::cout << SourceCode::ProgramName() << std::endl;
    std::cout << SourceCode::RevisionAtCompileTime() << std::endl;
    std::cout << SourceCode::TimeAtCompilation() << std::endl;
    std::cout << SourceCode::LastCommitted_Author() << std::endl;
    std::cout << SourceCode::LastCommitted_Revision() << std::endl;
    std::cout << SourceCode::LastCommitted_Date() << std::endl;

  }

}

namespace tut
{
  SourceCode_TestGroup SourceCode_Test("SourceCode_Test");
}
