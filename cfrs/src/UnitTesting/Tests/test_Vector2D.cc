#include <tut.h>
#include "Math/Vector2D.h"

using namespace std;

namespace tut

{
  // Data used for testing
  struct Data_Vector2D{
    Vector2D A;
 };

   /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   */
  typedef test_group<Data_Vector2D> Vector2D_TestGroup;
  typedef Vector2D_TestGroup::object  Vector2D_TestObject;
}

/**********************************************
 *                THE TESTS                   *
 *********************************************/

namespace tut
{

  /* Test 1: Assignment operator*/
  template<>
  template<>
  void Vector2D_TestObject::test<1>()
  {
    Vector2D B(2.0,3.0);
    ensure("Check X", A.x ==  0);
    ensure("Check Y", A.y ==  0);
    A = B;
    ensure("Check X", A.x ==  2.0);
    ensure("Check Y", A.y == 3.0);
  }

  /* Test 2: Copy Constructor*/
  template<>
  template<>
  void Vector2D_TestObject::test<2>()
  {
    A.x = 2.0; A.y = 3.0;
    Vector2D B = A;
    ensure("Copy Constructor", A ==  B);
  }


}

namespace tut
{
  Vector2D_TestGroup Vector2D_Test("Vector2D_Test");
}
