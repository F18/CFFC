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

    // generate a pseudo-inverse situation
    DenseMatrix A;
    ColumnVector B;

    // Constructor
    Data_HighOrder1D(): A(2,3), B(3){

      set_test_suite_path("HighOrderReconstruction/UnitTests/");
      set_local_input_path("test_TaylorDerivatives1D");
      set_local_output_path("test_TaylorDerivatives1D");

      // set geometry
      Cell1D.setloc(2.3);
      Cell1D.setsize(4.34);

      // set the pseudo-inverse situation
      A(0,0) = 1.0; A(0,1) = 2.0; A(0,2) = 3.0;
      A(1,0) = 1.0; A(1,1) = 4.0; A(1,2) = 3.0;
      
      B(0) = 1.0; B(1) = 2.34; B(2) = 34.0;
    
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
    set_test_name("Defalt Constructor (for piecewise constant)");

    HighOrder1D<double> HO;
    Cell1D_Uniform * Geom(NULL);

    ensure_equals("TD size", HO.NumberOfTaylorDerivatives(), 1);
    ensure_equals("TD power", HO.CellDeriv_InPosition(0).P1(), 0);
    ensure_equals("TD value I", HO.CellDeriv_InPosition(0).D(), 0.0);
    ensure_equals("TD value II", HO.CellDeriv(0), 0.0);
    ensure_equals("TD value III", HO.CellDeriv(0,1), 0.0);

    ensure_equals("GeomCoeff size", HO.CellGeomCoeff().size(), 1);
    ensure_equals("GeomCoeff power", HO.CellGeomCoeff().operator()(0,true,true,true).P1(), 0);
    ensure_equals("GeomCoeff value I", HO.CellGeomCoeff_InPosition(0), 1.0);
    ensure_equals("GeomCoeff value II", HO.CellGeomCoeff(0), 1.0);
    ensure_equals("Rings", HO.CellRings(), 0);
    ensure_equals("Monotonicity flag", HO.CellInadequateFit(1), OFF);
    ensure_distance("Smoothness indicator I", HO.CellSmoothnessIndicator(1), 0.0, 1.0e-14);
    ensure_distance("Smoothness indicator II", HO.CellSmoothnessIndicator(), 0.0, 1.0e-14);

    ensure("Geometry Pointer", HO.Geometry() == NULL);
  }

  /* Test 2:*/
  template<>
  template<>
  void HighOrder1D_object::test<2>()
  {
    set_test_name("Main Constructor(RecOrder=3, Cell1D)");

    Open_Output_File("HO.dat");

    // reset the geometry
    Cell1D.setloc(23.3);
    Cell1D.setsize(2.34);
    
    // set HighOrder1D variable
    HighOrder1D<double> HO(3,Cell1D);

    // set values of TD
    HO.CellDeriv(0) = 1.03; HO.CellDeriv(1) = 2.03;
    HO.CellDeriv(2) = 3.03; HO.CellDeriv(3) = 4.03;

    // == check
    out() << HO;

    RunRegressionTest("Main Constructor", "HO.dat", "TD_1D_MainConstructor.dat", 1.0e-14);
  }

  /* Test 3:*/
  template<>
  template<>
  void HighOrder1D_object::test<3>()
  {
    set_test_name("SetGeometryPointer()");

    // set HighOrder1D variable
    HighOrder1D<double> HO;

    // set the geometry pointer
    HO.SetGeometryPointer(Cell1D);

    // == check
    ensure_equals("geometry pointer", HO.CellGeometry(), Cell1D);
  }

  /* Test 4:*/
  template<>
  template<>
  void HighOrder1D_object::test<4>()
  {
    set_test_name("operator ==");

    // set HighOrder1D variables
    HighOrder1D<double> HO(3,Cell1D), HO_2(3,Cell1D);

    // set values of TD
    HO.CellDeriv(0) = 1.03; HO.CellDeriv(1) = 2.03;
    HO.CellDeriv(2) = 3.03; HO.CellDeriv(3) = 4.03;

    HO_2.CellDeriv(0) = 1.03; HO_2.CellDeriv(1) = 2.03;
    HO_2.CellDeriv(2) = 3.03; HO_2.CellDeriv(3) = 4.03;
    
    // == check
    ensure("==", (HO_2 == HO) == true);
  }

  /* Test 5:*/
  template<>
  template<>
  void HighOrder1D_object::test<5>()
  {
    set_test_name("operator ==, with LHS and Weights");

    // set HighOrder1D variables
    HighOrder1D<double> HO(3,Cell1D), HO_2(3,Cell1D);

    // set values of TD
    HO.CellDeriv(0) = 1.03; HO.CellDeriv(1) = 2.03;
    HO.CellDeriv(2) = 3.03; HO.CellDeriv(3) = 4.03;

    HO_2.CellDeriv(0) = 1.03; HO_2.CellDeriv(1) = 2.03;
    HO_2.CellDeriv(2) = 3.03; HO_2.CellDeriv(3) = 4.03;

    // set LHS and GeomWeights
    HO.LHS() = A; HO_2.LHS() = A;
    HO.GeomWeights() = B; HO_2.GeomWeights() = B;
    
    // == check
    ensure("==", (HO_2 == HO) == true);
  }

  /* Test 6:*/
  template<>
  template<>
  void HighOrder1D_object::test<6>()
  {
    set_test_name("operator !=");

    // set HighOrder1D variables
    HighOrder1D<double> HO(3,Cell1D), HO_2(3,Cell1D);

    // set values of TD
    HO.CellDeriv(0) = 1.03; HO.CellDeriv(1) = 2.03;
    HO.CellDeriv(2) = 3.03; HO.CellDeriv(3) = 4.03;

    HO_2.CellDeriv(0) = 1.03; HO_2.CellDeriv(1) = 2.03;
    HO_2.CellDeriv(2) = 3.034; HO_2.CellDeriv(3) = 4.03;

    // == check
    ensure("==", (HO_2 == HO) == false);
    ensure("!=", (HO_2 != HO) == true);
  }

  /* Test 7:*/
  template<>
  template<>
  void HighOrder1D_object::test<7>()
  {
    set_test_name("operator !=, with LHS and Weights");

    // set HighOrder1D variables
    HighOrder1D<double> HO(3,Cell1D), HO_2(3,Cell1D);

    // set values of TD
    HO.CellDeriv(0) = 1.03; HO.CellDeriv(1) = 2.03;
    HO.CellDeriv(2) = 3.03; HO.CellDeriv(3) = 4.03;

    HO_2.CellDeriv(0) = 1.03; HO_2.CellDeriv(1) = 2.03;
    HO_2.CellDeriv(2) = 3.03; HO_2.CellDeriv(3) = 4.03;

    // set LHS and GeomWeights
    HO.LHS() = A; 
    A(1,0) = -12.232;		// make a modification
    HO_2.LHS() = A;
    HO.GeomWeights() = B; HO_2.GeomWeights() = B;
    
    // == check
    ensure("==", (HO_2 == HO) == false);
    ensure("!=", (HO_2 != HO) == true);
  }

  /* Test 8:*/
  template<>
  template<>
  void HighOrder1D_object::test<8>()
  {
    set_test_name("Copy Constructor");

    // set HighOrder1D variable
    HighOrder1D<double> HO(3,Cell1D);

    // set values of TD
    HO.CellDeriv(0) = 1.03; HO.CellDeriv(1) = 2.03;
    HO.CellDeriv(2) = 3.03; HO.CellDeriv(3) = 4.03;

    // call copy constructor
    HighOrder1D<double> HO_2(HO);

    // == check if the objects are equal
    ensure_equals("Copy constructor", HO_2, HO);
  }

  /* Test 9:*/
  template<>
  template<>
  void HighOrder1D_object::test<9>()
  {
    set_test_name("Assignment Operator");

    // set HighOrder1D variable
    HighOrder1D<double> HO(3,Cell1D), HO_2;

    // set values of TD
    HO.CellDeriv(0) = 1.03; HO.CellDeriv(1) = 2.03;
    HO.CellDeriv(2) = 3.03; HO.CellDeriv(3) = 4.03;

    // call assignment operator
    HO_2 = HO;

    // == check if the objects are equal
    ensure_equals("Assignment operator", HO_2, HO);
  }

  /* Test 10:*/
  template<>
  template<>
  void HighOrder1D_object::test<10>()
  {
    set_test_name("operator >>");

    // set HighOrder1D variable
    HighOrder1D<double> HO(3,Cell1D);

    // set values of TD
    HO.CellDeriv(0) = 1.03; HO.CellDeriv(1) = 2.03;
    HO.CellDeriv(2) = 3.03; HO.CellDeriv(3) = 4.03;

    // == check
    Check_Input_Output_Operator(HO);
  }

  
  /* Test 11:*/
  template<>
  template<>
  void HighOrder1D_object::test<11>()
  {
    set_test_name("operator !=");

    // set HighOrder1D variable
    HighOrder1D<double> HO(3,Cell1D), HO_2(2,Cell1D);

    // set values of TD
    HO.CellDeriv(0) = 1.03; HO.CellDeriv(1) = 2.03;
    HO.CellDeriv(2) = 3.03; HO.CellDeriv(3) = 4.03;

    HO_2.CellDeriv(0) = 1.03; HO_2.CellDeriv(1) = 2.03;
    HO_2.CellDeriv(2) = 3.03;

    // == check
    ensure_not("==", HO == HO_2);
  }

  /* Test 12:*/
  template<>
  template<>
  void HighOrder1D_object::test<12>()
  {
    set_test_name("SolutionAtCoordinates()");

    // set geometry
    Cell1D.setloc(2.3);
    Cell1D.setsize(4.34);

    // set HighOrder1D variable
    HighOrder1D<double> HO(3,Cell1D);

    // set values of TD
    HO.CellDeriv(0) = 1.03; HO.CellDeriv(1) = 2.03;
    HO.CellDeriv(2) = 3.03; HO.CellDeriv(3) = 4.03;

    // == check
    ensure_distance("solution", HO.SolutionAtCoordinates(0.13), -30.28693439, 1.0e-14);
    ensure_distance("solution", HO.SolutionAtCoordinates(0.13,1), -30.28693439, 1.0e-14);
  }

  /* Test 13:*/
  template<>
  template<>
  void HighOrder1D_object::test<13>()
  {
    set_test_name("left_state()");

    // set geometry
    Cell1D.setloc(2.3);
    Cell1D.setsize(4.34);

    // set HighOrder1D variable
    HighOrder1D<double> HO(3,Cell1D);

    // set values of TD
    HO.CellDeriv(0) = 1.03; HO.CellDeriv(1) = 2.03;
    HO.CellDeriv(2) = 3.03; HO.CellDeriv(3) = 4.03;

    // == check
    ensure_distance("left interface solution", HO.left_state(), -30.28693439, 1.0e-14);
  }

  /* Test 14:*/
  template<>
  template<>
  void HighOrder1D_object::test<14>()
  {
    set_test_name("right_state()");

    // set geometry
    Cell1D.setloc(2.3);
    Cell1D.setsize(4.34);

    // set HighOrder1D variable
    HighOrder1D<double> HO(3,Cell1D);

    // set values of TD
    HO.CellDeriv(0) = 1.03; HO.CellDeriv(1) = 2.03;
    HO.CellDeriv(2) = 3.03; HO.CellDeriv(3) = 4.03;

    // == check
    ensure_distance("right interface solution", HO.right_state(), 60.88286839, 1.0e-14);
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

