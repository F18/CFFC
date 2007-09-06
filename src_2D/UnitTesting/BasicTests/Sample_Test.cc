#include "TestData.h"

namespace tut
{

  class Data_Test : public TestData {

    // Add more local variables
  public:
    double varA;

    Data_Test(){
      // Set the global path for this test suit
      // It's a good practice to put it in the 
      // constructor of the data class.
      // Then it is automatically set for each individual test.
      // Declare it relative to the /src_2D directory,
      // otherwise the framework might not find the input
      // and output files.
      set_test_suit_path("UnitTesting/BasicTests");
    }
  };

  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. TestGroup) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_Test> TestGroup;
  typedef TestGroup::object TestObject;
}

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


/***********************************************
 *  THE TESTS THAT ARE PART OF THE TEST SUITE  *
 ***********************************************/

namespace tut
{

  /* Test 1:*/
  template<>
  template<>
  void TestObject::test<1>()
  {

    set_test_name("Example 1");

    // *************** Test code ***************** 
    ensure ("NameOfTheCondition", 1 <= 1); // test condition

  }

  /*Test 2:*/
  template<>
  template<>
  void TestObject::test<2>()
  {

    set_test_name("Example 2");

    // *************** Test code ***************** 
    varA = -1.0;
    ensure_equals("Ensure_equal example", varA, -1.0); //test condition

  }

  /*Test 3:*/
  template<>
  template<>
  void TestObject::test<3>()
  {

    set_test_name("Example 3: how to use local path for output and input");
    set_local_output_path("test_MyTest"); // <-- this path is relative to the global one
    set_local_input_path("test_MyTest"); // <-- this path is relative to the global one


    // *************** Test code ***************** 
    double a;
    Open_Input_File("TestFile");

    in() >> a;			// the in() function can be used for reading from the input stream

    ensure_distance("Ensure_distance example", a, 1231313.131, 0.001 ); //test condition
  }


}


namespace tut
{
  TestGroup TestSuitObject("Simple_Example");
}
