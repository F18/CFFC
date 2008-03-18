/*!\file test_HighOrder2D_Double.cc
  \brief Regression tests for template class HighOrder2D with double datatype. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../HighOrder2D.h"
#include "../../Math/UnitTests/TestFunctions/TestFunctions_2D.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_HighOrder2D : public TestData {

    // Local variables
  public:

    // generate geometry 

    // generate a pseudo-inverse situation
    DenseMatrix A;
    ColumnVector B;

    // Constructor
    Data_HighOrder2D();

  private:
    
  };

  // Constructor
  Data_HighOrder2D::Data_HighOrder2D(): A(2,3), B(3){

    set_test_suite_path("HighOrderReconstruction/UnitTests/");
    set_local_input_path("test_TaylorDerivatives1D");
    set_local_output_path("test_TaylorDerivatives1D");
    
    // set CENO_Execution_Mode to default values
    CENO_Execution_Mode::SetDefaults();

    // set geometry
    
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
  typedef test_group<Data_HighOrder2D> HighOrder2D_TestSuite;
  typedef HighOrder2D_TestSuite::object HighOrder2D_object;


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
  void HighOrder2D_object::test<1>()
  {
    set_test_name("Defalt Constructor");

    HighOrder2D<double> HO;

    ensure("TD", HO.TaylorDeriv() == NULL);
    ensure("SI", HO.SmoothnessIndicator() == NULL);
    ensure("LimitedCell", HO.InadequateFit() == NULL);
    ensure("LHS", HO.LHS_Inv() == NULL);
    ensure("CENO_Geometric_Weights", HO.GeomWeights() == NULL);
    ensure_equals("Order of reconstruction", HO.RecOrder(), -1);
    ensure_equals("Rings", HO.Rings(), 0);
    ensure_equals("PseudoInverseFlag", HO.IsPseudoInversePreComputed(), false);    
    ensure("Geometry Pointer", HO.Geometry() == NULL);
  }

  /* Test 2:*/
  template<>
  template<>
  void HighOrder2D_object::test<2>()
  {
    set_test_name("TaylorDerivativesSize()");

    HighOrder2D<double> HO;

    ensure_equals("k=-1", HO.TaylorDerivativesSize(-1), 0);
    ensure_equals("k=0", HO.TaylorDerivativesSize(0), 1);
    ensure_equals("k=1", HO.TaylorDerivativesSize(1), 3);
    ensure_equals("k=2", HO.TaylorDerivativesSize(2), 6);
    ensure_equals("k=3", HO.TaylorDerivativesSize(3), 10);
    ensure_equals("k=4", HO.TaylorDerivativesSize(4), 15);
  }

  /* Test 3:*/
  template<>
  template<>
  void HighOrder2D_object::test<3>()
  {
    set_test_name("MinimumNghost()");

    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;

    HighOrder2D<double> HO;

    ensure_equals("k=-1", HO.MinimumNghost(-1), 0);
    ensure_equals("k=0", HO.MinimumNghost(0), 1);
    ensure_equals("k=1", HO.MinimumNghost(1), 3);
    ensure_equals("k=2", HO.MinimumNghost(2), 5);
    ensure_equals("k=3", HO.MinimumNghost(3), 5);
    ensure_equals("k=4", HO.MinimumNghost(4), 5);
    ensure_equals("k=10", HO.MinimumNghost(10), 20);
  }

  /* Test 4:*/
  template<>
  template<>
  void HighOrder2D_object::test<4>()
  {
    set_test_name("MinimumNghost()");

    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = ON;

    HighOrder2D<double> HO;

    ensure_equals("k=-1", HO.MinimumNghost(-1), 0);
    ensure_equals("k=0", HO.MinimumNghost(0), 1);
    ensure_equals("k=1", HO.MinimumNghost(1), 3);
    ensure_equals("k=2", HO.MinimumNghost(2), 4);
    ensure_equals("k=3", HO.MinimumNghost(3), 4);
    ensure_equals("k=4", HO.MinimumNghost(4), 4);
    ensure_equals("k=10", HO.MinimumNghost(10), 20);
  }

  /* Test 5:*/
  template<>
  template<>
  void HighOrder2D_object::test<5>()
  {
    set_test_name("NumberOfRings()");

    HighOrder2D<double> HO;

    ensure_equals("k=-1",HO.NumberOfRings(0), 0);
    ensure_equals("k=0", HO.NumberOfRings(1), 0);
    ensure_equals("k=1", HO.NumberOfRings(3), 1);
    ensure_equals("k=2", HO.NumberOfRings(6), 2);
    ensure_equals("k=3", HO.NumberOfRings(10), 2);
    ensure_equals("k=4", HO.NumberOfRings(15), 2);
    ensure_equals("k=5", HO.NumberOfRings(21), 3);
  }

  /* Test 6:*/
  template<>
  template<>
  void HighOrder2D_object::test<6>()
  {
    set_test_name("StencilSize()");

    HighOrder2D<double> HO;

    ensure_equals("StencilSize, k=-1",HO.StencilSize(0), 1);
    ensure_equals("StencilSize, k=0", HO.StencilSize(0), 1);
    ensure_equals("StencilSize, k=1", HO.StencilSize(1), 9);
    ensure_equals("StencilSize, k=2", HO.StencilSize(2), 25);
    ensure_equals("StencilSize, k=3", HO.StencilSize(3), 25);
    ensure_equals("StencilSize, k=4", HO.StencilSize(4), 25);
    ensure_equals("StencilSize, k=5", HO.StencilSize(5), 49);
  }

  /* Test 7:*/
  template<>
  template<>
  void HighOrder2D_object::test<7>()
  {
    set_test_name("allocate()");

    HighOrder2D<double> HO;

    //    HO.allocate(2,3,2,true,1);

    //    ensure_equals("TD value", HO.CellTaylorDerivState(1,1,0,0), 0.0);

#if 0
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
#endif
  }

}



// Test suite constructor
tut::HighOrder2D_TestSuite HighOrder2DTestSuite("Template Class:HighOrder2D && double");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

