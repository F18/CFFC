/*!\file test_HighOrder1D_Double.cc
  \brief Regression tests for template class HighOrder1D with double datatype. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../HighOrder1D.h"
#include "../../Math/UnitTests/TestFunctions/TestFunctions_1D.h"

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
    Data_HighOrder1D();

  private:
    
  };

  // Constructor
  Data_HighOrder1D::Data_HighOrder1D(): A(2,3), B(3){

    set_test_suite_path("HighOrderReconstruction/UnitTests/");
    set_local_input_path("test_TaylorDerivatives1D");
    set_local_output_path("test_TaylorDerivatives1D");
    
    // set CENO_Execution_Mode to default values
    CENO_Execution_Mode::SetDefaults();

    // set geometry
    Cell1D.setloc(2.3);
    Cell1D.setsize(4.34);
    
    // set the pseudo-inverse situation
    A(0,0) = 1.0; A(0,1) = 2.0; A(0,2) = 3.0;
    A(1,0) = 1.0; A(1,1) = 4.0; A(1,2) = 3.0;
    
    B(0) = 1.0; B(1) = 2.34; B(2) = 34.0;
    
  }

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
    ensure_equals("Smoothness indicator I", HO.CellSmoothnessIndicator(1), 0.0);
    ensure_equals("Smoothness indicator II", HO.CellSmoothnessIndicator(), 0.0);

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

    CENO_Execution_Mode::CENO_SPEED_EFFICIENT = ON;

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

  /* Test 15:*/
  template<>
  template<>
  void HighOrder1D_object::test<15>()
  {
    set_test_name("Associate geometry");

    // reset the geometry
    Cell1D.setloc(23.3);
    Cell1D.setsize(2.34);
    
    // set HighOrder1D variable
    HighOrder1D<double> HO(3,Cell1D);

    // set values of TD
    HO.CellDeriv(0) = 1.03; HO.CellDeriv(1) = 2.03;
    HO.CellDeriv(2) = 3.03; HO.CellDeriv(3) = 4.03;

    Cell1D_Uniform Cell2;
    Cell2.setloc(21.3);
    Cell2.setsize(2.1);

    HO.AssociateGeometry(Cell2);

    // == check
    ensure_distance("Second-order geometric moment", HO.CellGeomCoeff(2), 0.3675, 1.0e-14);
  }

  /* Test 16:*/
  template<>
  template<>
  void HighOrder1D_object::test<16>()
  {
    set_test_name("InitializeVariable() && RECONSTRUCTION_CENO");
    Open_Output_File("CENO_Current.dat");

    CENO_Execution_Mode::CENO_SPEED_EFFICIENT = ON;

    // set HighOrder1D variable
    HighOrder1D<double> HO;

    HO.InitializeVariable(3,RECONSTRUCTION_CENO);

    // set values of TD
    HO.CellDeriv(0) = 1.03; HO.CellDeriv(1) = 2.03;
    HO.CellDeriv(2) = 3.03; HO.CellDeriv(3) = 4.03;

    // == check
    out() << HO;

    RunRegressionTest("initialized variable", "CENO_Current.dat", "InitializeVariable_CENO_Master.dat", 1.0e-14);
  }

  /* Test 17:*/
  template<>
  template<>
  void HighOrder1D_object::test<17>()
  {
    set_test_name("InitializeVariable() && RECONSTRUCTION_ENO");
    Open_Output_File("ENO_Current.dat");

    // set HighOrder1D variable
    HighOrder1D<double> HO;

    HO.InitializeVariable(3,RECONSTRUCTION_ENO);

    // set values of TD
    HO.CellDeriv(0) = 1.03; HO.CellDeriv(1) = 2.03;
    HO.CellDeriv(2) = 3.03; HO.CellDeriv(3) = 4.03;

    // == check
    out() << HO;

    RunRegressionTest("initialized variable", "ENO_Current.dat", "InitializeVariable_ENO_Master.dat", 1.0e-14);
  }

  /* Test 18:*/
  template<>
  template<>
  void HighOrder1D_object::test<18>()
  {
    set_test_name("InitializeVariable() && RECONSTRUCTION_LEAST_SQUARES");
    Open_Output_File("LEAST_SQUARES_Current.dat");

    // set HighOrder1D variable
    HighOrder1D<double> HO;

    HO.InitializeVariable(3,RECONSTRUCTION_LEAST_SQUARES);

    // == check
    out() << HO;

    RunRegressionTest("initialized variable", "LEAST_SQUARES_Current.dat",
		      "InitializeVariable_LEAST_SQUARES_Master.dat", 1.0e-14);
  }

  /* Test 19:*/
  template<>
  template<>
  void HighOrder1D_object::test<19>()
  {
    set_test_name("operator << >> && RECONSTRUCTION_LEAST_SQUARES");

    // set HighOrder1D variable
    HighOrder1D<double> HO;

    HO.InitializeVariable(3,RECONSTRUCTION_LEAST_SQUARES);

    // == check
    Check_Input_Output_Operator("initialized variable for RECONSTRUCTION_LEAST_SQUARES", HO);
  }

  /* Test 20:*/
  template<>
  template<>
  void HighOrder1D_object::test<20>()
  {
    set_test_name("operator << >> && ENO");

    // set HighOrder1D variable
    HighOrder1D<double> HO;

    HO.InitializeVariable(3,RECONSTRUCTION_ENO);

    // == check
    Check_Input_Output_Operator("initialized variable for RECONSTRUCTION_ENO", HO);
  }

  /* Test 21:*/
  template<>
  template<>
  void HighOrder1D_object::test<21>()
  {
    set_test_name("ComputeSolutionError()");

    // set HighOrder1D variable
    HighOrder1D<double> HO, HO2;

    HO.InitializeVariable(3,RECONSTRUCTION_ENO);
    HO2.InitializeVariable(2,RECONSTRUCTION_CENO);

    Cell1D_Uniform Cell;
    Cell.setloc(1.3);
    Cell.setsize(2.1);

    HO.AssociateGeometry(Cell);
    HO2.AssociateGeometry(Cell);

    // set values of TD
    HO.CellDeriv(0) = 1.03; HO.CellDeriv(1) = 2.03;
    HO.CellDeriv(2) = 3.03; HO.CellDeriv(3) = 4.03;

    HO2.CellDeriv(0) = 1.03; HO2.CellDeriv(1) = 2.03;
    HO2.CellDeriv(2) = 3.03;

    // compute error
    double ErrorL1, ErrorL2;
    ErrorL1 = HO.ComputeSolutionErrorL1(mapped_function(Test_Example9, ErrorL1, -5.0, 5.0, -5.0, 5.0), 1);
    ErrorL2 = HO.ComputeSolutionErrorL2(mapped_function(Test_Example9, ErrorL2, -5.0, 5.0, -5.0, 5.0), 1);

    // == check
    ensure_distance("L1 error norm", ErrorL1, 4.40320869669138, 1.0e-10);
    ensure_distance("L2 error norm", ErrorL2, 21.90907364279464, 1.0e-10);

    ErrorL1 = HO.ComputeSolutionErrorL1(HO2,1);
    ensure_distance("Relative error norm", ErrorL1, 2.44924509375, 1.0e-12);
  }

}



// Test suite constructor
tut::HighOrder1D_TestSuite HighOrder1DTestSuite("Template Class:HighOrder1D && double");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

