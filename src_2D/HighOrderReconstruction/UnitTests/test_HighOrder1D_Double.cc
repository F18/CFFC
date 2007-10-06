/*!\file test_HighOrder1DState.cc
  \brief Regression tests for template class HighOrder1DState. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../HighOrder1D.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_HighOrder1D : public TestData {

    // Local variables
  public:

    // generate a geometry set
    Cell1D_Uniform Cell1D;

    // Constructor
    Data_HighOrder1D(){

      set_test_suite_path("HighOrderReconstruction/UnitTests/");

      // set geometry
      Cell1D.setloc(2.3);
      Cell1D.setsize(4.34);
    
    }

  private:
    
  };

  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. 'NAME'_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_HighOrder1D> HighOrder1D_TestSuite;
  typedef HighOrder1D_TestSuite::object HighOrder1D_object;


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


  /*********************************************
     TESTING OPERATIONS:
        -->  ensure("ConditionName", test_condition);
        -->  ensure_not("ConditionName", test_condition);
        -->  ensure_equals("ConditionName", expected_value, got_value);
                 // operators '<<' and '!=' must be defined!!!
        -->  ensure_distance("ConditionName", expected_value, got_value, tolerance);
	         // operators '<<', '>=' and '<=' must be defined!!!
        -->  fail("Message");

     Obs: "ConditionName" is optional
  */


  /* Test 1:*/
  template<>
  template<>
  void HighOrder1D_object::test<1>()
  {
    set_test_name("Defalt Constructor");

    HighOrder1D<double> HO;
    Cell1D_Uniform * Geom(NULL);

    ensure("Geometry Pointer", HO.Geometry() == NULL);
  }

  /* Test 2:*/
  template<>
  template<>
  void HighOrder1D_object::test<2>()
  {
    set_test_name("Copy Constructor");

    // set HighOrder1D variable
    HighOrder1D<double> HO;

    // set the geometry pointer
    HO.SetGeometryPointer(Cell1D);

    // call copy constructor
    HighOrder1D<double> HO_2(HO);

    // == check if both geometry pointers point to the same object
    ensure_equals("Copy constructor (geometry pointer)", HO_2.CellGeometry(), HO.CellGeometry());
  }

  /* Test 3:*/
  template<>
  template<>
  void HighOrder1D_object::test<3>()
  {
    set_test_name("Assignment Operator");

    // set HighOrder1D variable
    HighOrder1D<double> HO, HO_2;

    // set the geometry pointer
    HO.SetGeometryPointer(Cell1D);

    // call assignment operator
    HO_2 = HO;

    // == check if both geometry pointers point to the same object
    ensure_equals("Assignment operator (geometry pointer)", HO_2.CellGeometry(), HO.CellGeometry());
  }

  /* Test 4:*/
  template<>
  template<>
  void HighOrder1D_object::test<4>()
  {
    set_test_name("SetGeometryPointer()");

    // set HighOrder1D variable
    HighOrder1D<double> HO;

    // set the geometry pointer
    HO.SetGeometryPointer(Cell1D);

    // == check
    ensure_equals("geometry pointer", HO.CellGeometry(), Cell1D);
  }

  /* Test 5:*/
  template<>
  template<>
  void HighOrder1D_object::test<5>()
  {
    set_test_name("operator ==");

    // set HighOrder1D variable
    HighOrder1D<double> HO, HO_2;

    // set the geometry pointer
    HO.SetGeometryPointer(Cell1D);

    HO_2 = HO;

    // == check
    ensure("==", (HO_2 == HO) == true);
  }

  /* Test 6:*/
  template<>
  template<>
  void HighOrder1D_object::test<6>()
  {
    set_test_name("operator <<");

    // set HighOrder1D variable
    HighOrder1D<double> HO, HO_2;

    // set the geometry pointer
    HO.SetGeometryPointer(Cell1D);

    // == check

    //    cout << HO << endl;

  }

  /* Test 7:*/
  template<>
  template<>
  void HighOrder1D_object::test<7>()
  {
    set_test_name("operator >>");

    // set HighOrder1D variable
    HighOrder1D<double> HO, HO_2;

    // set the geometry pointer
    HO.SetGeometryPointer(Cell1D);

    // == check

    //    cout << HO << endl;

  }

  /* Test 8:*/
  template<>
  template<>
  void HighOrder1D_object::test<8>()
  {
    set_test_name("Constructor(RecOrder=3, Geometry)");

    // set HighOrder1D variable
    HighOrder1D<double> HO(3,Cell1D);

    // == check

    cout << HO << endl;

  }




}



// Test suite constructor
tut::HighOrder1D_TestSuite HighOrder1DTestSuite("Template Class:HighOrder1D");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

