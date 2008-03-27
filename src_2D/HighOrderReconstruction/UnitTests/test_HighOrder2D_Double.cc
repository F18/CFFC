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
#include "../HighOrder2D_Input.h"

namespace tut
{

  double Function_Test(const double &x, const double &y){
    return 1;
  }

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
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");
    
    // set CENO_Execution_Mode to default values
    CENO_Execution_Mode::SetDefaults();

    // set High-order parameters to default values
    HighOrder2D_Input::SetDefaults();

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
    set_test_name("getTaylorDerivativesSize()");

    HighOrder2D<double> HO;

    ensure_equals("k=-1", HO.getTaylorDerivativesSize(-1), 0);
    ensure_equals("k=0", HO.getTaylorDerivativesSize(0), 1);
    ensure_equals("k=1", HO.getTaylorDerivativesSize(1), 3);
    ensure_equals("k=2", HO.getTaylorDerivativesSize(2), 6);
    ensure_equals("k=3", HO.getTaylorDerivativesSize(3), 10);
    ensure_equals("k=4", HO.getTaylorDerivativesSize(4), 15);
  }

  /* Test 3:*/
  template<>
  template<>
  void HighOrder2D_object::test<3>()
  {
    set_test_name("getMinimumNghost()");

    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;

    HighOrder2D<double> HO;

    ensure_equals("k=-1", HO.getMinimumNghost(-1), 0);
    ensure_equals("k=0", HO.getMinimumNghost(0), 1);
    ensure_equals("k=1", HO.getMinimumNghost(1), 3);
    ensure_equals("k=2", HO.getMinimumNghost(2), 5);
    ensure_equals("k=3", HO.getMinimumNghost(3), 5);
    ensure_equals("k=4", HO.getMinimumNghost(4), 5);
    ensure_equals("k=10", HO.getMinimumNghost(10), 20);
  }

  /* Test 4:*/
  template<>
  template<>
  void HighOrder2D_object::test<4>()
  {
    set_test_name("getMinimumNghost()");

    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = ON;

    HighOrder2D<double> HO;

    ensure_equals("k=-1", HO.getMinimumNghost(-1), 0);
    ensure_equals("k=0", HO.getMinimumNghost(0), 1);
    ensure_equals("k=1", HO.getMinimumNghost(1), 3);
    ensure_equals("k=2", HO.getMinimumNghost(2), 4);
    ensure_equals("k=3", HO.getMinimumNghost(3), 4);
    ensure_equals("k=4", HO.getMinimumNghost(4), 4);
    ensure_equals("k=10", HO.getMinimumNghost(10), 20);
  }

  /* Test 5:*/
  template<>
  template<>
  void HighOrder2D_object::test<5>()
  {
    set_test_name("getNumberOfRings()");

    HighOrder2D<double> HO;

    ensure_equals("k=-1",HO.getNumberOfRings(0), 0);
    ensure_equals("k=0", HO.getNumberOfRings(1), 0);
    ensure_equals("k=1", HO.getNumberOfRings(3), 1);
    ensure_equals("k=2", HO.getNumberOfRings(6), 2);
    ensure_equals("k=3", HO.getNumberOfRings(10), 2);
    ensure_equals("k=4", HO.getNumberOfRings(15), 2);
    ensure_equals("k=5", HO.getNumberOfRings(21), 3);
  }

  /* Test 6:*/
  template<>
  template<>
  void HighOrder2D_object::test<6>()
  {
    set_test_name("getStencilSize()");

    HighOrder2D<double> HO;

    ensure_equals("getStencilSize, k=-1",HO.getStencilSize(0), 1);
    ensure_equals("getStencilSize, k=0", HO.getStencilSize(0), 1);
    ensure_equals("getStencilSize, k=1", HO.getStencilSize(1), 9);
    ensure_equals("getStencilSize, k=2", HO.getStencilSize(2), 25);
    ensure_equals("getStencilSize, k=3", HO.getStencilSize(3), 25);
    ensure_equals("getStencilSize, k=4", HO.getStencilSize(4), 25);
    ensure_equals("getStencilSize, k=5", HO.getStencilSize(5), 49);
  }

  /* Test 7:*/
  template<>
  template<>
  void HighOrder2D_object::test<7>()
  {
    set_test_name("allocate()");

    HighOrder2D<double> HO;

    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = ON;

    try {

      // allocate wrong number of cells
      HO.allocate(2,1,3,true,1);

      fail("Fail to detect inconsistent dimensions");

    } catch (runtime_error){
      // test successful
    }
  }

  /* Test 8:*/
  template<>
  template<>
  void HighOrder2D_object::test<8>()
  {
    set_test_name("allocate()");

    HighOrder2D<double> HO;

    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = ON;

    try {

      // allocate wrong number of ghost cells
      HO.allocate(2,2,2,true,1);

      fail("Fail to detect inconsistent dimensions");

    } catch (runtime_error){
      // test successful
    }
  }

  /* Test 9:*/
  template<>
  template<>
  void HighOrder2D_object::test<9>()
  {
    set_test_name("allocate()");

    HighOrder2D<double> HO;

    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = ON;

    try {

      // allocate wrong reconstruction order
      HO.allocate(2,1,3,true,-2);

      fail("Fail to detect inconsistent dimensions");

    } catch (runtime_error){
      // test successful
    }
  }

  /* Test 10:*/
  template<>
  template<>
  void HighOrder2D_object::test<10>()
  {
    set_test_name("deallocate() with memory allocation");

    HighOrder2D<double> HO;

    // call
    HO.deallocate();

    // check
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

  /* Test 11:*/
  template<>
  template<>
  void HighOrder2D_object::test<11>()
  {
    set_test_name("allocate() without pseudo-inverse");

    HighOrder2D<double> HO;
    int NCi, NCj, Nghost, RecOrder;
    NCi = 2; NCj = 3; Nghost = 3; RecOrder = 1;
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = ON;

    // Allocate memory for the high-order WITHOUT pseudo-inverse
    HO.allocate(NCi,NCj,Nghost,false,RecOrder);

    // == check 
    ensure_equals("Stencil size", HO.getStencilSize(), 9);
    ensure_equals("Taylor deriv. size", HO.getTaylorDerivativesSize(), 3);

    // == check TD containers
    ensure("TD", HO.TaylorDeriv() != NULL);
    ensure_equals("TD size", HO.NumberOfTaylorDerivatives(), 3);

    ensure_equals("Pos. 0, TD power", HO.CellTaylorDeriv(1,1,0).P1(), 0);
    ensure_equals("Pos. 0, TD power", HO.CellTaylorDeriv(1,1,0).P2(), 0);
    ensure_equals("Pos. 0, TD value I", HO.CellTaylorDeriv(1,1,0).D(), 0.0);
    ensure_equals("Pos. 1, TD power", HO.CellTaylorDeriv(1,1,1).P1(), 0);
    ensure_equals("Pos. 1, TD power", HO.CellTaylorDeriv(1,1,1).P2(), 1);
    ensure_equals("Pos. 1, TD value I", HO.CellTaylorDeriv(1,1,1).D(), 0.0);
    ensure_equals("Pos. 2, TD power", HO.CellTaylorDeriv(1,1,2).P1(), 1);
    ensure_equals("Pos. 2, TD power", HO.CellTaylorDeriv(1,1,2).P2(), 0);
    ensure_equals("Pos. 2, TD value I", HO.CellTaylorDeriv(1,1,2).D(), 0.0);

    ensure_equals("TD value II", HO.CellTaylorDerivValue(1,1,1,0,0), 0.0);
    ensure_equals("TD value III", HO.CellTaylorDerivState(1,1,0,1), 0.0);

    // check Reconstruction Order and number of rings
    ensure_equals("Rings", HO.Rings(), 1);
    ensure_equals("RecOrder", HO.RecOrder(), 1);

    // check monotonicity containers
    std::vector<short int> MonotonicityContainer(1); MonotonicityContainer[0] = OFF;
    std::vector<double> MonotonicityValue(1); MonotonicityValue[0] = 0.0;
    ensure("LimitedCell", HO.InadequateFit() != NULL);
    ensure("Monotonicity container", HO.CellInadequateFit(1,1) == MonotonicityContainer);
    ensure_equals("Monotonicity value", HO.CellInadequateFitValue(1,1,1), OFF);

    ensure("SI", HO.SmoothnessIndicator() != NULL);
    ensure("SI container", HO.CellSmoothnessIndicator(1,1) == MonotonicityValue);
    ensure_equals("SI value", HO.CellSmoothnessIndicatorValue(1,1,1), 0.0);
    
    // check pseudo-inverse data
    ensure("CENO_LHS", HO.LHS_Inv() == NULL);
    ensure("CENO_Geometric_Weights", HO.GeomWeights() == NULL);

    ensure("Geometry Pointer", HO.Geometry() == NULL);

    // == check containers for the last cell with high-order containers
    ensure_equals("Last Cell, Pos. 0, TD power", HO.CellTaylorDeriv(6,7,0).P1(), 0);
    ensure_equals("Last Cell, Pos. 0, TD power", HO.CellTaylorDeriv(6,7,0).P2(), 0);
    ensure_equals("Last Cell, Pos. 0, TD value I", HO.CellTaylorDeriv(6,7,0).D(), 0.0);
    ensure_equals("Last Cell, Pos. 1, TD power", HO.CellTaylorDeriv(6,7,1).P1(), 0);
    ensure_equals("Last Cell, Pos. 1, TD power", HO.CellTaylorDeriv(6,7,1).P2(), 1);
    ensure_equals("Last Cell, Pos. 1, TD value I", HO.CellTaylorDeriv(6,7,1).D(), 0.0);
    ensure_equals("Last Cell, Pos. 2, TD power", HO.CellTaylorDeriv(6,7,2).P1(), 1);
    ensure_equals("Last Cell, Pos. 2, TD power", HO.CellTaylorDeriv(6,7,2).P2(), 0);
    ensure_equals("Last Cell, Pos. 2, TD value I", HO.CellTaylorDeriv(6,7,2).D(), 0.0);

    ensure_equals("Last Cell, TD value II", HO.CellTaylorDerivValue(6,7,1,0,0), 0.0);
    ensure_equals("Last Cell, TD value III", HO.CellTaylorDerivState(6,7,0,1), 0.0);

    // check monotonicity containers
    ensure("Last Cell, Monotonicity container", HO.CellInadequateFit(6,7) == MonotonicityContainer);
    ensure_equals("Last Cell, Monotonicity value", HO.CellInadequateFitValue(6,7,1), OFF);
    ensure("Last Cell, SI container", HO.CellSmoothnessIndicator(6,7) == MonotonicityValue);
    ensure_equals("Last Cell, SI value", HO.CellSmoothnessIndicatorValue(6,7,1), 0.0);
  }

  /* Test 12:*/
  template<>
  template<>
  void HighOrder2D_object::test<12>()
  {
    set_test_name("allocate() with pseudo-inverse");

    HighOrder2D<double> HO;
    int NCi, NCj, Nghost, RecOrder;
    NCi = 2; NCj = 3; Nghost = 3; RecOrder = 1;
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = ON;

    // Allocate memory for the high-order WITH pseudo-inverse
    HO.allocate(NCi,NCj,Nghost,true,RecOrder);

    // == check 
    ensure_equals("Stencil size", HO.getStencilSize(), 9);
    ensure_equals("Taylor deriv. size", HO.getTaylorDerivativesSize(), 3);

    // == check TD containers
    ensure("TD", HO.TaylorDeriv() != NULL);
    ensure_equals("TD size", HO.NumberOfTaylorDerivatives(), 3);

    ensure_equals("Pos. 0, TD power", HO.CellTaylorDeriv(1,1,0).P1(), 0);
    ensure_equals("Pos. 0, TD power", HO.CellTaylorDeriv(1,1,0).P2(), 0);
    ensure_equals("Pos. 0, TD value I", HO.CellTaylorDeriv(1,1,0).D(), 0.0);
    ensure_equals("Pos. 1, TD power", HO.CellTaylorDeriv(1,1,1).P1(), 0);
    ensure_equals("Pos. 1, TD power", HO.CellTaylorDeriv(1,1,1).P2(), 1);
    ensure_equals("Pos. 1, TD value I", HO.CellTaylorDeriv(1,1,1).D(), 0.0);
    ensure_equals("Pos. 2, TD power", HO.CellTaylorDeriv(1,1,2).P1(), 1);
    ensure_equals("Pos. 2, TD power", HO.CellTaylorDeriv(1,1,2).P2(), 0);
    ensure_equals("Pos. 2, TD value I", HO.CellTaylorDeriv(1,1,2).D(), 0.0);

    ensure_equals("TD value II", HO.CellTaylorDerivValue(1,1,1,0,0), 0.0);
    ensure_equals("TD value III", HO.CellTaylorDerivState(1,1,0,1), 0.0);

    // check Reconstruction Order and number of rings
    ensure_equals("Rings", HO.Rings(), 1);
    ensure_equals("RecOrder", HO.RecOrder(), 1);

    // check monotonicity containers
    std::vector<short int> MonotonicityContainer(1); MonotonicityContainer[0] = OFF;
    std::vector<double> MonotonicityValue(1); MonotonicityValue[0] = 0.0;
    ensure("LimitedCell", HO.InadequateFit() != NULL);
    ensure("Monotonicity container", HO.CellInadequateFit(1,1) == MonotonicityContainer);
    ensure_equals("Monotonicity value", HO.CellInadequateFitValue(1,1,1), OFF);

    ensure("SI", HO.SmoothnessIndicator() != NULL);
    ensure("SI container", HO.CellSmoothnessIndicator(1,1) == MonotonicityValue);
    ensure_equals("SI value", HO.CellSmoothnessIndicatorValue(1,1,1), 0.0);
    
    // check pseudo-inverse data
    ensure("CENO_LHS", HO.LHS_Inv() != NULL);
    ensure_equals("LHS value", HO.Cell_LHS_Inv(1,1).size(0), 8);
    ensure_equals("LHS value", HO.Cell_LHS_Inv(1,1).size(1), 2);
    ensure("CENO_Geometric_Weights", HO.GeomWeights() != NULL);
    ensure_equals("CENO_Geometric_Weights", HO.GeomWeights(1,1).size(), 9);
    ensure_equals("CENO_Geometric_Weights", HO.GeomWeightValue(1,1,8), 0.0);

    ensure("Geometry Pointer", HO.Geometry() == NULL);

    // == check containers for the last cell with high-order containers
    ensure_equals("Last Cell, Pos. 0, TD power", HO.CellTaylorDeriv(6,7,0).P1(), 0);
    ensure_equals("Last Cell, Pos. 0, TD power", HO.CellTaylorDeriv(6,7,0).P2(), 0);
    ensure_equals("Last Cell, Pos. 0, TD value I", HO.CellTaylorDeriv(6,7,0).D(), 0.0);
    ensure_equals("Last Cell, Pos. 1, TD power", HO.CellTaylorDeriv(6,7,1).P1(), 0);
    ensure_equals("Last Cell, Pos. 1, TD power", HO.CellTaylorDeriv(6,7,1).P2(), 1);
    ensure_equals("Last Cell, Pos. 1, TD value I", HO.CellTaylorDeriv(6,7,1).D(), 0.0);
    ensure_equals("Last Cell, Pos. 2, TD power", HO.CellTaylorDeriv(6,7,2).P1(), 1);
    ensure_equals("Last Cell, Pos. 2, TD power", HO.CellTaylorDeriv(6,7,2).P2(), 0);
    ensure_equals("Last Cell, Pos. 2, TD value I", HO.CellTaylorDeriv(6,7,2).D(), 0.0);

    ensure_equals("Last Cell, TD value II", HO.CellTaylorDerivValue(6,7,1,0,0), 0.0);
    ensure_equals("Last Cell, TD value III", HO.CellTaylorDerivState(6,7,0,1), 0.0);

    // check monotonicity containers
    ensure("Last Cell, Monotonicity container", HO.CellInadequateFit(6,7) == MonotonicityContainer);
    ensure_equals("Last Cell, Monotonicity value", HO.CellInadequateFitValue(6,7,1), OFF);
    ensure("Last Cell, SI container", HO.CellSmoothnessIndicator(6,7) == MonotonicityValue);
    ensure_equals("Last Cell, SI value", HO.CellSmoothnessIndicatorValue(6,7,1), 0.0);
  }

  /* Test 13:*/
  template<>
  template<>
  void HighOrder2D_object::test<13>()
  {
    set_test_name("allocate(), SI with all neighbours");

    HighOrder2D<double> HO;
    int NCi, NCj, Nghost, RecOrder;
    NCi = 2; NCj = 3; Nghost = 5; RecOrder = 4;
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;

    // Allocate memory for the high-order WITH pseudo-inverse
    HO.allocate(NCi,NCj,Nghost,true,RecOrder);

    // == check 
    ensure_equals("Stencil size", HO.getStencilSize(), 25);
    ensure_equals("Taylor deriv. size", HO.getTaylorDerivativesSize(), 15);

    // == check TD containers
    ensure("TD", HO.TaylorDeriv() != NULL);
    ensure_equals("TD size", HO.NumberOfTaylorDerivatives(), 15);

    ensure_equals("Pos. 0, TD power", HO.CellTaylorDeriv(2,2,0).P1(), 0);
    ensure_equals("Pos. 0, TD power", HO.CellTaylorDeriv(2,2,0).P2(), 0);
    ensure_equals("Pos. 0, TD value I", HO.CellTaylorDeriv(2,2,0).D(), 0.0);
    ensure_equals("Pos. 1, TD power", HO.CellTaylorDeriv(2,2,1).P1(), 0);
    ensure_equals("Pos. 1, TD power", HO.CellTaylorDeriv(2,2,1).P2(), 1);
    ensure_equals("Pos. 1, TD value I", HO.CellTaylorDeriv(2,2,1).D(), 0.0);
    ensure_equals("Pos. 14, TD power", HO.CellTaylorDeriv(2,2,14).P1(), 4);
    ensure_equals("Pos. 14, TD power", HO.CellTaylorDeriv(2,2,14).P2(), 0);
    ensure_equals("Pos. 14, TD value I", HO.CellTaylorDeriv(2,2,14).D(), 0.0);

    ensure_equals("TD value II", HO.CellTaylorDerivValue(2,2,1,0,0), 0.0);
    ensure_equals("TD value III", HO.CellTaylorDerivState(2,2,0,1), 0.0);

    // check Reconstruction Order and number of rings
    ensure_equals("Rings", HO.Rings(), 2);
    ensure_equals("RecOrder", HO.RecOrder(), 4);

    // check monotonicity containers
    std::vector<short int> MonotonicityContainer(1); MonotonicityContainer[0] = OFF;
    std::vector<double> MonotonicityValue(1); MonotonicityValue[0] = 0.0;
    ensure("LimitedCell", HO.InadequateFit() != NULL);
    ensure("Monotonicity container", HO.CellInadequateFit(2,2) == MonotonicityContainer);
    ensure_equals("Monotonicity value", HO.CellInadequateFitValue(2,2,1), OFF);

    ensure("SI", HO.SmoothnessIndicator() != NULL);
    ensure("SI container", HO.CellSmoothnessIndicator(2,2) == MonotonicityValue);
    ensure_equals("SI value", HO.CellSmoothnessIndicatorValue(2,2,1), 0.0);
    
    // check pseudo-inverse data
    ensure("CENO_LHS", HO.LHS_Inv() != NULL);
    ensure_equals("LHS value", HO.Cell_LHS_Inv(2,2).size(0), 24);
    ensure_equals("LHS value", HO.Cell_LHS_Inv(2,2).size(1), 14);
    ensure("CENO_Geometric_Weights", HO.GeomWeights() != NULL);
    ensure_equals("CENO_Geometric_Weights", HO.GeomWeights(2,2).size(), 25);
    ensure_equals("CENO_Geometric_Weights", HO.GeomWeightValue(2,2,24), 0.0);

    ensure("Geometry Pointer", HO.Geometry() == NULL);

    // == check containers for the last cell with high-order containers
    ensure_equals("Last Cell, Pos. 0, TD power", HO.CellTaylorDeriv(9,10,0).P1(), 0);
    ensure_equals("Last Cell, Pos. 0, TD power", HO.CellTaylorDeriv(9,10,0).P2(), 0);
    ensure_equals("Last Cell, Pos. 0, TD value I", HO.CellTaylorDeriv(9,10,0).D(), 0.0);
    ensure_equals("Last Cell, Pos. 1, TD power", HO.CellTaylorDeriv(9,10,1).P1(), 0);
    ensure_equals("Last Cell, Pos. 1, TD power", HO.CellTaylorDeriv(9,10,1).P2(), 1);
    ensure_equals("Last Cell, Pos. 1, TD value I", HO.CellTaylorDeriv(9,10,1).D(), 0.0);
    ensure_equals("Last Cell, Pos. 14, TD power", HO.CellTaylorDeriv(9,10,14).P1(), 4);
    ensure_equals("Last Cell, Pos. 14, TD power", HO.CellTaylorDeriv(9,10,14).P2(), 0);
    ensure_equals("Last Cell, Pos. 14, TD value I", HO.CellTaylorDeriv(9,10,14).D(), 0.0);

    ensure_equals("Last Cell, TD value II", HO.CellTaylorDerivValue(9,10,1,0,0), 0.0);
    ensure_equals("Last Cell, TD value III", HO.CellTaylorDerivState(9,10,0,1), 0.0);

    // check monotonicity containers
    ensure("Last Cell, Monotonicity container", HO.CellInadequateFit(9,10) == MonotonicityContainer);
    ensure_equals("Last Cell, Monotonicity value", HO.CellInadequateFitValue(9,10,1), OFF);
    ensure("Last Cell, SI container", HO.CellSmoothnessIndicator(9,10) == MonotonicityValue);
    ensure_equals("Last Cell, SI value", HO.CellSmoothnessIndicatorValue(9,10,1), 0.0);
  }

  /* Test 14:*/
  template<>
  template<>
  void HighOrder2D_object::test<14>()
  {
    set_test_name("allocate(), Change execution flag");

    HighOrder2D<double> HO;
    int NCi, NCj, Nghost, RecOrder;
    NCi = 2; NCj = 3; Nghost = 5; RecOrder = 4;
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;

    // Allocate memory for the high-order WITH pseudo-inverse
    HO.allocate(NCi,NCj,Nghost,true,RecOrder);
  
    ensure_equals("Nghost high-order", HO.NghostHO(), 3);

    // Change the execution flag
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = ON;

    // Re-allocate memory for the high-order WITH pseudo-inverse
    HO.allocate(NCi,NCj,Nghost,true,RecOrder);

    // == check 
    ensure_equals("Nghost high-order", HO.NghostHO(), 2);
    ensure_equals("Stencil size", HO.getStencilSize(), 25);
    ensure_equals("Taylor deriv. size", HO.getTaylorDerivativesSize(), 15);

  
    // == check TD containers
    ensure("TD", HO.TaylorDeriv() != NULL);
    ensure_equals("TD size", HO.NumberOfTaylorDerivatives(), 15);

    ensure_equals("Pos. 0, TD power", HO.CellTaylorDeriv(3,3,0).P1(), 0);
    ensure_equals("Pos. 0, TD power", HO.CellTaylorDeriv(3,3,0).P2(), 0);
    ensure_equals("Pos. 0, TD value I", HO.CellTaylorDeriv(3,3,0).D(), 0.0);
    ensure_equals("Pos. 1, TD power", HO.CellTaylorDeriv(3,3,1).P1(), 0);
    ensure_equals("Pos. 1, TD power", HO.CellTaylorDeriv(3,3,1).P2(), 1);
    ensure_equals("Pos. 1, TD value I", HO.CellTaylorDeriv(3,3,1).D(), 0.0);
    ensure_equals("Pos. 14, TD power", HO.CellTaylorDeriv(3,3,14).P1(), 4);
    ensure_equals("Pos. 14, TD power", HO.CellTaylorDeriv(3,3,14).P2(), 0);
    ensure_equals("Pos. 14, TD value I", HO.CellTaylorDeriv(3,3,14).D(), 0.0);

    ensure_equals("TD value II", HO.CellTaylorDerivValue(3,3,1,0,0), 0.0);
    ensure_equals("TD value III", HO.CellTaylorDerivState(3,3,0,1), 0.0);

    // check Reconstruction Order and number of rings
    ensure_equals("Rings", HO.Rings(), 2);
    ensure_equals("RecOrder", HO.RecOrder(), 4);

    // check monotonicity containers
    std::vector<short int> MonotonicityContainer(1); MonotonicityContainer[0] = OFF;
    std::vector<double> MonotonicityValue(1); MonotonicityValue[0] = 0.0;
    ensure("LimitedCell", HO.InadequateFit() != NULL);
    ensure("Monotonicity container", HO.CellInadequateFit(3,3) == MonotonicityContainer);
    ensure_equals("Monotonicity value", HO.CellInadequateFitValue(3,3,1), OFF);

    ensure("SI", HO.SmoothnessIndicator() != NULL);
    ensure("SI container", HO.CellSmoothnessIndicator(3,3) == MonotonicityValue);
    ensure_equals("SI value", HO.CellSmoothnessIndicatorValue(3,3,1), 0.0);
    
    // check pseudo-inverse data
    ensure("CENO_LHS", HO.LHS_Inv() != NULL);
    ensure_equals("LHS value", HO.Cell_LHS_Inv(3,3).size(0), 24);
    ensure_equals("LHS value", HO.Cell_LHS_Inv(3,3).size(1), 14);
    ensure("CENO_Geometric_Weights", HO.GeomWeights() != NULL);
    ensure_equals("CENO_Geometric_Weights", HO.GeomWeights(3,3).size(), 25);
    ensure_equals("CENO_Geometric_Weights", HO.GeomWeightValue(3,3,24), 0.0);
  
    ensure("Geometry Pointer", HO.Geometry() == NULL);
  
    // == check containers for the last cell with high-order containers
    ensure_equals("Last Cell, Pos. 0, TD power", HO.CellTaylorDeriv(8,9,0).P1(), 0);
    ensure_equals("Last Cell, Pos. 0, TD power", HO.CellTaylorDeriv(8,9,0).P2(), 0); 
    ensure_equals("Last Cell, Pos. 0, TD value I", HO.CellTaylorDeriv(8,9,0).D(), 0.0);
    ensure_equals("Last Cell, Pos. 1, TD power", HO.CellTaylorDeriv(8,9,1).P1(), 0);
    ensure_equals("Last Cell, Pos. 1, TD power", HO.CellTaylorDeriv(8,9,1).P2(), 1);
    ensure_equals("Last Cell, Pos. 1, TD value I", HO.CellTaylorDeriv(8,9,1).D(), 0.0);
    ensure_equals("Last Cell, Pos. 14, TD power", HO.CellTaylorDeriv(8,9,14).P1(), 4);
    ensure_equals("Last Cell, Pos. 14, TD power", HO.CellTaylorDeriv(8,9,14).P2(), 0);
    ensure_equals("Last Cell, Pos. 14, TD value I", HO.CellTaylorDeriv(8,9,14).D(), 0.0);
  
    ensure_equals("Last Cell, TD value II", HO.CellTaylorDerivValue(8,9,1,0,0), 0.0);
    ensure_equals("Last Cell, TD value III", HO.CellTaylorDerivState(8,9,0,1), 0.0);

    // check monotonicity containers
    ensure("Last Cell, Monotonicity container", HO.CellInadequateFit(8,9) == MonotonicityContainer);
    ensure_equals("Last Cell, Monotonicity value", HO.CellInadequateFitValue(8,9,1), OFF);
    ensure("Last Cell, SI container", HO.CellSmoothnessIndicator(8,9) == MonotonicityValue);
    ensure_equals("Last Cell, SI value", HO.CellSmoothnessIndicatorValue(8,9,1), 0.0);
  }

  /* Test 15:*/
  template<>
  template<>
  void HighOrder2D_object::test<15>()
  {
    set_test_name("allocate(), Change reconstruction order");

    HighOrder2D<double> HO;
    int NCi, NCj, Nghost, RecOrder;
    NCi = 2; NCj = 3; Nghost = 5; RecOrder = 4;
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = ON;

    // Allocate memory for the high-order WITH pseudo-inverse
    HO.allocate(NCi,NCj,Nghost,true,RecOrder);

    ensure_equals("Nghost high-order", HO.NghostHO(), 2);

    // Change the reconstruction order
    RecOrder = 2;

    // Re-allocate memory for the high-order WITH pseudo-inverse
    HO.allocate(NCi,NCj,Nghost,true,RecOrder);

    // == check 
    ensure_equals("Nghost high-order", HO.NghostHO(), 2);
    ensure_equals("Stencil size", HO.getStencilSize(), 25);
    ensure_equals("Taylor deriv. size", HO.getTaylorDerivativesSize(), 6);


    // == check TD containers
    ensure("TD", HO.TaylorDeriv() != NULL);
    ensure_equals("TD size", HO.NumberOfTaylorDerivatives(), 6);

    ensure_equals("Pos. 0, TD power", HO.CellTaylorDeriv(3,3,0).P1(), 0);
    ensure_equals("Pos. 0, TD power", HO.CellTaylorDeriv(3,3,0).P2(), 0);
    ensure_equals("Pos. 0, TD value I", HO.CellTaylorDeriv(3,3,0).D(), 0.0);
    ensure_equals("Pos. 1, TD power", HO.CellTaylorDeriv(3,3,1).P1(), 0);
    ensure_equals("Pos. 1, TD power", HO.CellTaylorDeriv(3,3,1).P2(), 1);
    ensure_equals("Pos. 1, TD value I", HO.CellTaylorDeriv(3,3,1).D(), 0.0);
    ensure_equals("Pos. 5, TD power", HO.CellTaylorDeriv(3,3,5).P1(), 2);
    ensure_equals("Pos. 5, TD power", HO.CellTaylorDeriv(3,3,5).P2(), 0);
    ensure_equals("Pos. 5, TD value I", HO.CellTaylorDeriv(3,3,5).D(), 0.0);

    ensure_equals("TD value II", HO.CellTaylorDerivValue(3,3,1,0,0), 0.0);
    ensure_equals("TD value III", HO.CellTaylorDerivState(3,3,0,1), 0.0);

    // check Reconstruction Order and number of rings
    ensure_equals("Rings", HO.Rings(), 2);
    ensure_equals("RecOrder", HO.RecOrder(), 2);

    // check monotonicity containers
    std::vector<short int> MonotonicityContainer(1); MonotonicityContainer[0] = OFF;
    std::vector<double> MonotonicityValue(1); MonotonicityValue[0] = 0.0;
    ensure("LimitedCell", HO.InadequateFit() != NULL);
    ensure("Monotonicity container", HO.CellInadequateFit(3,3) == MonotonicityContainer);
    ensure_equals("Monotonicity value", HO.CellInadequateFitValue(3,3,1), OFF);

    ensure("SI", HO.SmoothnessIndicator() != NULL);
    ensure("SI container", HO.CellSmoothnessIndicator(3,3) == MonotonicityValue);
    ensure_equals("SI value", HO.CellSmoothnessIndicatorValue(3,3,1), 0.0);
    
    // check pseudo-inverse data
    ensure("CENO_LHS", HO.LHS_Inv() != NULL);
    ensure_equals("LHS value", HO.Cell_LHS_Inv(3,3).size(0), 24);
    ensure_equals("LHS value", HO.Cell_LHS_Inv(3,3).size(1), 5);
    ensure("CENO_Geometric_Weights", HO.GeomWeights() != NULL);
    ensure_equals("CENO_Geometric_Weights", HO.GeomWeights(3,3).size(), 25);
    ensure_equals("CENO_Geometric_Weights", HO.GeomWeightValue(3,3,24), 0.0);

    ensure("Geometry Pointer", HO.Geometry() == NULL);

    // == check containers for the last cell with high-order containers
    ensure_equals("Last Cell, Pos. 0, TD power", HO.CellTaylorDeriv(8,9,0).P1(), 0);
    ensure_equals("Last Cell, Pos. 0, TD power", HO.CellTaylorDeriv(8,9,0).P2(), 0);
    ensure_equals("Last Cell, Pos. 0, TD value I", HO.CellTaylorDeriv(8,9,0).D(), 0.0);
    ensure_equals("Last Cell, Pos. 1, TD power", HO.CellTaylorDeriv(8,9,1).P1(), 0);
    ensure_equals("Last Cell, Pos. 1, TD power", HO.CellTaylorDeriv(8,9,1).P2(), 1);
    ensure_equals("Last Cell, Pos. 1, TD value I", HO.CellTaylorDeriv(8,9,1).D(), 0.0);
    ensure_equals("Last Cell, Pos. 2, TD power", HO.CellTaylorDeriv(8,9,5).P1(), 2);
    ensure_equals("Last Cell, Pos. 2, TD power", HO.CellTaylorDeriv(8,9,5).P2(), 0);
    ensure_equals("Last Cell, Pos. 2, TD value I", HO.CellTaylorDeriv(8,9,5).D(), 0.0);

    ensure_equals("Last Cell, TD value II", HO.CellTaylorDerivValue(8,9,1,0,0), 0.0);
    ensure_equals("Last Cell, TD value III", HO.CellTaylorDerivState(8,9,0,1), 0.0);

    // check monotonicity containers
    ensure("Last Cell, Monotonicity container", HO.CellInadequateFit(8,9) == MonotonicityContainer);
    ensure_equals("Last Cell, Monotonicity value", HO.CellInadequateFitValue(8,9,1), OFF);
    ensure("Last Cell, SI container", HO.CellSmoothnessIndicator(8,9) == MonotonicityValue);
    ensure_equals("Last Cell, SI value", HO.CellSmoothnessIndicatorValue(8,9,1), 0.0);
  }

  /* Test 16:*/
  template<>
  template<>
  void HighOrder2D_object::test<16>()
  {
    set_test_name("Output operator");

    RunRegression = ON;

    HighOrder2D<double> HO;
    int NCi, NCj, Nghost, RecOrder;
    NCi = 2; NCj = 3; Nghost = 5; RecOrder = 4;
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;

    // Allocate memory for the high-order WITH pseudo-inverse
    HO.allocate(NCi,NCj,Nghost,true,RecOrder);

    MasterFile = "Output_Operator.dat";
    CurrentFile = "Current_Output_Operator.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);
      
      out() << HO;

      RunRegressionTest("operator <<", CurrentFile, MasterFile, 1.0e-12);

    } else {
      // Generate the master file
      Open_Output_File(MasterFile);

      out() << HO;
    }
      
  }

  /* Test 17:*/
  template<>
  template<>
  void HighOrder2D_object::test<17>()
  {
    set_test_name("Input-output operators");

    HighOrder2D<double> HO;
    int NCi, NCj, Nghost, RecOrder;
    NCi = 2; NCj = 3; Nghost = 5; RecOrder = 4;
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = ON;

    // Allocate memory for the high-order WITH pseudo-inverse
    HO.allocate(NCi,NCj,Nghost,true,RecOrder);

    Check_Input_Output_Operator(HO);
  }

  /* Test 18:*/
  template<>
  template<>
  void HighOrder2D_object::test<18>()
  {
    set_test_name("Copy constructor I");

    HighOrder2D<double> HO;
    int NCi, NCj, Nghost, RecOrder;
    NCi = 2; NCj = 3; Nghost = 5; RecOrder = 4;
    int i,j;
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = ON;

    // Generate a geometry
    Grid2D_Quad_Block_HO MyGrid;

    // Allocate memory for the high-order WITH pseudo-inverse
    HO.allocate(NCi,NCj,Nghost,true,RecOrder);
    HO.SetGeometryPointer(MyGrid);

    // Use copy constructor
    HighOrder2D<double> HO_Copy(HO);

    // check
    ensure_equals("Geometry", HO_Copy.Geometry(), HO.Geometry());
    ensure_equals("RecOrder", HO_Copy.RecOrder(), HO.RecOrder());
    ensure_equals("Rings", HO_Copy.Rings(), HO.Rings());

    for (i = 3; i<= 8; ++i) {
      for (j = 3; j<= 9; ++j) {
	ensure_equals("TD", HO_Copy.CellTaylorDeriv(i,j), HO.CellTaylorDeriv(i,j));
	ensure("LimitedCell", HO_Copy.CellInadequateFit(i,j) == HO.CellInadequateFit(i,j));
	ensure("SI", HO_Copy.CellSmoothnessIndicator(i,j) == HO.CellSmoothnessIndicator(i,j));
	ensure_equals("LHS", HO_Copy.Cell_LHS_Inv(i,j), HO.Cell_LHS_Inv(i,j));
	ensure("GeomWeights", HO_Copy.GeomWeights(i,j) == HO.GeomWeights(i,j));
      }
    }
  }

  /* Test 19:*/
  template<>
  template<>
  void HighOrder2D_object::test<19>()
  {
    set_test_name("Copy constructor II");

    HighOrder2D<double> HO;
    int NCi, NCj, Nghost, RecOrder;
    NCi = 2; NCj = 3; Nghost = 5; RecOrder = 4;
    int i,j;
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;

    // Generate a geometry
    Grid2D_Quad_Block_HO MyGrid;

    // Allocate memory for the high-order WITH pseudo-inverse
    HO.allocate(NCi,NCj,Nghost,false,RecOrder);
    HO.SetGeometryPointer(MyGrid);

    // Use copy constructor
    HighOrder2D<double> HO_Copy(HO);

    // check
    ensure_equals("Geometry", HO_Copy.Geometry(), HO.Geometry());
    ensure_equals("RecOrder", HO_Copy.RecOrder(), HO.RecOrder());
    ensure_equals("Rings", HO_Copy.Rings(), HO.Rings());

    for (i = 2; i<= 9; ++i) {
      for (j = 2; j<= 10; ++j) {
	ensure_equals("TD", HO_Copy.CellTaylorDeriv(i,j), HO.CellTaylorDeriv(i,j));
	ensure("LimitedCell", HO_Copy.CellInadequateFit(i,j) == HO.CellInadequateFit(i,j));
	ensure("SI", HO_Copy.CellSmoothnessIndicator(i,j) == HO.CellSmoothnessIndicator(i,j));
      }
    }
  }

  /* Test 20:*/
  template<>
  template<>
  void HighOrder2D_object::test<20>()
  {
    set_test_name("Assignment operator I");

    HighOrder2D<double> HO;
    int NCi, NCj, Nghost, RecOrder;
    NCi = 2; NCj = 3; Nghost = 5; RecOrder = 4;
    int i,j;
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = ON;

    // Generate a geometry
    Grid2D_Quad_Block_HO MyGrid;

    // Allocate memory for the high-order WITH pseudo-inverse
    HO.allocate(NCi,NCj,Nghost,true,RecOrder);
    HO.SetGeometryPointer(MyGrid);

    // Use assignment operator
    HighOrder2D<double> HO_Copy;
    HO_Copy = HO;
    HO_Copy.SetGeometryPointer(MyGrid);

    // check
    ensure_equals("Geometry", HO_Copy.Geometry(), HO.Geometry());
    ensure_equals("RecOrder", HO_Copy.RecOrder(), HO.RecOrder());
    ensure_equals("Rings", HO_Copy.Rings(), HO.Rings());

    for (i = 3; i<= 8; ++i) {
      for (j = 3; j<= 9; ++j) {
	ensure_equals("TD", HO_Copy.CellTaylorDeriv(i,j), HO.CellTaylorDeriv(i,j));
	ensure("LimitedCell", HO_Copy.CellInadequateFit(i,j) == HO.CellInadequateFit(i,j));
	ensure("SI", HO_Copy.CellSmoothnessIndicator(i,j) == HO.CellSmoothnessIndicator(i,j));
	ensure_equals("LHS", HO_Copy.Cell_LHS_Inv(i,j), HO.Cell_LHS_Inv(i,j));
	ensure("GeomWeights", HO_Copy.GeomWeights(i,j) == HO.GeomWeights(i,j));
      }
    }
  }

  /* Test 21:*/
  template<>
  template<>
  void HighOrder2D_object::test<21>()
  {
    set_test_name("Assignment operator II");

    try {
      HighOrder2D<double> HO;
      int NCi, NCj, Nghost, RecOrder;
      NCi = 2; NCj = 3; Nghost = 5; RecOrder = 4;
      int i,j;
      
      // Set execution mode
      CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
      CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = ON;
      
      // Generate a geometry
      Grid2D_Quad_Block_HO MyGrid;

      // Allocate memory for the high-order WITH pseudo-inverse
      HO.allocate(NCi,NCj,Nghost,true,RecOrder);
      HO.SetGeometryPointer(MyGrid);

      // Change Execution_Mode class setting
      CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;

      // Use assignment operator
      HighOrder2D<double> HO_Copy;
      HO_Copy = HO;

      // fail
      fail("Fail to detect different settings in Execution_Mode!");

    } catch (runtime_error){
      // test successful
    }
  }

  /* Test 22:*/
  template<>
  template<>
  void HighOrder2D_object::test<22>()
  {
    set_test_name("Evaluate the polynomial interpolant");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    HighOrder2D<double> HO;
    int RecOrder(4);
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = ON;

    // Generate a geometry
    Grid2D_Quad_Block_HO Grid;

    // Read the geometry from input file
    Open_Input_File("CartesianMesh.dat");
    in() >> Grid;

    // Initialize high-order variable
    HO.InitializeVariable(RecOrder,Grid,true);

    // Set the Taylor derivatives
    HO.CellTaylorDerivState(3,3,0,0) = 1.0;
    HO.CellTaylorDerivState(3,3,0,1) = 2.0;
    HO.CellTaylorDerivState(3,3,0,2) = 3.0;
    HO.CellTaylorDerivState(3,3,0,3) = 4.0;
    HO.CellTaylorDerivState(3,3,0,4) = 5.0;
    HO.CellTaylorDerivState(3,3,1,0) = 1.0;
    HO.CellTaylorDerivState(3,3,1,1) = 2.0;
    HO.CellTaylorDerivState(3,3,1,2) = 3.0;
    HO.CellTaylorDerivState(3,3,1,3) = 4.0;
    HO.CellTaylorDerivState(3,3,2,0) = 1.0;
    HO.CellTaylorDerivState(3,3,2,1) = 2.0;
    HO.CellTaylorDerivState(3,3,2,2) = 3.0;
    HO.CellTaylorDerivState(3,3,3,0) = 1.0;
    HO.CellTaylorDerivState(3,3,3,1) = 2.0;
    HO.CellTaylorDerivState(3,3,4,0) = 1.0;

    // == check that constrained reconstruction is not required
    ensure_equals("Block Flag", HO.IsConstrainedReconstructionRequired(), false);
    ensure_equals("West Bnd Flag", HO.IsWestConstrainedReconstructionRequired(), false);
    ensure_equals("East Bnd Flag", HO.IsEastConstrainedReconstructionRequired(), false);
    ensure_equals("North Bnd Flag", HO.IsNorthConstrainedReconstructionRequired(), false);
    ensure_equals("South Bnd Flag", HO.IsSouthConstrainedReconstructionRequired(), false);

    // == check solution
    double Result = 0.59509999999999985;
    Vector2D Point(-1.2, -1.8);	// This corresponds to DeltaX = 0.1 and DeltaY = 0.5

    ensure_distance("Point value with cell (3,3)",
		    HO.SolutionAtCoordinates(3,3,Point.x, Point.y), Result, AcceptedError(Result));

    ensure_distance("Point value with cell (3,3)",
		    HO.SolutionAtCoordinates(3,3,Point.x, Point.y, 1), Result, AcceptedError(Result));
  }

  /* Test 23:*/
  template<>
  template<>
  void HighOrder2D_object::test<23>()
  {
    set_test_name("Integrate over the cell domain");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    HighOrder2D<double> HO;
    int RecOrder(4);
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;

    // Generate a geometry
    Grid2D_Quad_Block_HO Grid;

    // Read the geometry from input file
    Open_Input_File("CartesianMesh.dat");
    in() >> Grid;

    // Initialize high-order variable
    HO.InitializeVariable(RecOrder,Grid,true);

    // Integrate over the cell (3,3)
    double Result = Grid.CellArea(3,3);

    // == check integration
    ensure_distance("Area", HO.IntegrateOverTheCell(3,3,Function_Test,10, Result), Result, AcceptedError(Result));

    // == check that constrained reconstruction is not required
    ensure_equals("Block Flag", HO.IsConstrainedReconstructionRequired(), false);
    ensure_equals("West Bnd Flag", HO.IsWestConstrainedReconstructionRequired(), false);
    ensure_equals("East Bnd Flag", HO.IsEastConstrainedReconstructionRequired(), false);
    ensure_equals("North Bnd Flag", HO.IsNorthConstrainedReconstructionRequired(), false);
    ensure_equals("South Bnd Flag", HO.IsSouthConstrainedReconstructionRequired(), false);

    // Change spline type
    Grid.BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.BndEastSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    
    // Update
    HO.AssociateGeometry(Grid);

    // == check that constrained reconstruction is required
    ensure_equals("Block Flag", HO.IsConstrainedReconstructionRequired(), true);
    ensure_equals("West Bnd Flag", HO.IsWestConstrainedReconstructionRequired(), false);
    ensure_equals("East Bnd Flag", HO.IsEastConstrainedReconstructionRequired(), true);
    ensure_equals("North Bnd Flag", HO.IsNorthConstrainedReconstructionRequired(), true);
    ensure_equals("South Bnd Flag", HO.IsSouthConstrainedReconstructionRequired(), false);

  }

  /* Test 24:*/
  template<>
  template<>
  void HighOrder2D_object::test<24>()
  {
    set_test_name("Calculate pseudo-inverse");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    RunRegression = ON;

    HighOrder2D<double> HO;
    int RecOrder(4);
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;

    // Generate a geometry
    Grid2D_Quad_Block_HO Grid;

    // Read the geometry from input file
    Open_Input_File("CartesianMesh.dat");
    in() >> Grid;

    // Initialize high-order variable
    HO.InitializeVariable(RecOrder,Grid,true);

    MasterFile  = "PseudoInverse.dat";
    CurrentFile = "Current_PseudoInverse.dat";

    if (RunRegression){

      // open output file
      Open_Output_File(CurrentFile);

      ensure_equals("Pseudo-inv Flag I", HO.IsPseudoInversePreComputed(), true);

      // output the pseudo-inverse of cell (5,7)
      Print_File(HO.Cell_LHS_Inv(5,7), out());

      // switch order
      HO.SetReconstructionOrder(1);

      ensure_equals("Pseudo-inv Flag II", HO.IsPseudoInversePreComputed(), true);

      // output the pseudo-inverse of cell (5,7)
      Print_File(HO.Cell_LHS_Inv(5,7), out());

      RunRegressionTest("Pseudo-inverse cell (5,7)", CurrentFile, MasterFile, 1.0e-10, 1.0e-10);

    } else {

      // open output file
      Open_Output_File(MasterFile);

      // output the pseudo-inverse of cell (5,7)
      Print_File(HO.Cell_LHS_Inv(5,7), out());

      // switch order
      HO.SetReconstructionOrder(1);

      // output the pseudo-inverse of cell (5,7)
      Print_File(HO.Cell_LHS_Inv(5,7), out());
    }
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

