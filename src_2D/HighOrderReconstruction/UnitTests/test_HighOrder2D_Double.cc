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
#include "../../Grid/UnitTests/HO_Grid2DQuadMultiBlock_InputForTesting.h"

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
    typedef HighOrder2D<double> HighOrderVariable;

    // generate geometry 

    // generate a pseudo-inverse situation
    DenseMatrix A;
    ColumnVector B;
    int rings;

    Grid2DTesting_Input_Parameters IP;

    int error_flag;

    // Constructor
    Data_HighOrder2D();

    // Check consistency of reconstruction stencils for a multi-block mesh
    void CheckHighOrderReconstructionStencilConsistency(const HighOrderVariable &CheckedHO,
							const int &i_Start, const int &i_End,
							const int &j_Start, const int &j_End,
							const HighOrderVariable &MasterHO,
							const int &i_Master, const int &j_Master,
							const std::string & BaseMsg = "");

    // Create Mesh
    void CreateMesh(Grid2D_Quad_MultiBlock_HO & _MeshBlk_,
		    Grid2DTesting_Input_Parameters & IP) throw(std::runtime_error);

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
    Spline2D_HO::ResetCounter();
    Spline2DInterval_HO::Set_Default_Parameters();

    Grid2D_Quad_Block_HO::setLowOrderBoundaryRepresentation();
    
    // set the pseudo-inverse situation
    A(0,0) = 1.0; A(0,1) = 2.0; A(0,2) = 3.0;
    A(1,0) = 1.0; A(1,1) = 4.0; A(1,2) = 3.0;
    
    B(0) = 1.0; B(1) = 2.34; B(2) = 34.0;
    
  }

  // === CheckHighOrderReconstructionStencilConsistency()
  void Data_HighOrder2D::CheckHighOrderReconstructionStencilConsistency(const HighOrderVariable &CheckedHO,
									const int &i_Start, const int &i_End,
									const int &j_Start, const int &j_End,
									const HighOrderVariable &MasterHO,
									const int &i_Master, const int &j_Master,
									const std::string & BaseMsg){

    int iCell,jCell;		// cell indexes for the checked block
    int iMast,jMast;		// cell indexes for the master block that corresponds to the CheckedBlock cell
    int iShift, jShift;		// i- and j-shift between the indexes of the two blocks
    bool ICond, JCond;		// indicators for how to loop over indexes

    IndexType i_index_Master, j_index_Master, i_index_Checked, j_index_Checked;
    bool EqualStencils(true);
    int Counter;


    // Determine looping conditions
    ICond = (i_End - i_Start) > 0;
    JCond = (j_End - j_Start) > 0;

    // Determine index shifts
    iShift = i_Master - i_Start;
    jShift = j_Master - j_Start;
   

    // Compare geometric properties
    for (iCell = i_Start; ICond? (iCell<=i_End): (iCell>=i_End); ICond? (++iCell): (--iCell)){
      for (jCell = j_Start; JCond? (jCell<=j_End): (jCell>=j_End); JCond? (++jCell): (--jCell)){

	iMast = iCell + iShift;
	jMast = jCell + jShift;

	// Determine stencil for the master cell
	MasterHO.SetDeviatedReconstructionStencil(iMast, jMast,
						  i_index_Master, j_index_Master,
						  rings);

	// Determine stencil for the checked cell
	CheckedHO.SetDeviatedReconstructionStencil(iCell, jCell,
						   i_index_Checked, j_index_Checked,
						   rings);

	// === Check stencils
	for (Counter = 0; Counter<i_index_Master.size() ; ++Counter){
	  if ( ((i_index_Master[Counter]-iMast) != (i_index_Checked[Counter]-iCell)) || 
	       ((j_index_Master[Counter]-jMast) != (j_index_Checked[Counter]-jCell)) ){
	    EqualStencils = false;
	  }
	}

	if (EqualStencils == false){
	  ostm() << "Equal Stencils, " << BaseMsg << "\n"
		 << "Checked Cell (" << iCell << "," << jCell << ")"; 
	  CheckedHO.displayDeviatedReconstructionStencil(ostm(), iCell, jCell, rings);
	  
	  ostm() << "Master Cell (" << iMast << "," << jMast << ")"; 
	  MasterHO.displayDeviatedReconstructionStencil(ostm(), iMast, jMast, rings);
	  ostm() << "\n";

	  ensure_equals(ostm().str(), EqualStencils, true);
	}

      }
    }
  }


  // === CreateMesh()
  void Data_HighOrder2D::CreateMesh(Grid2D_Quad_MultiBlock_HO & _MeshBlk_,
				    Grid2DTesting_Input_Parameters & IP) throw(std::runtime_error){
    
    // Ensure that the highest reconstruction order was set correctly
    HighOrder2D_Input::Set_Final_Parameters(IP);

    /* Initialize all static variables within the class */
    if (IP.IncludeHighOrderBoundariesRepresentation == ON){
      Grid2D_Quad_Block_HO::setHighOrderBoundaryRepresentation();
    } else {
      Grid2D_Quad_Block_HO::setLowOrderBoundaryRepresentation();
    }

    error_flag = _MeshBlk_.Multi_Block_Grid(IP);
    
    if (error_flag) {
      // try to output the nodes
      _MeshBlk_.Output_Nodes_Tecplot_Using_IP(IP);
      throw runtime_error("CreateMesh() ERROR: Unable to create valid multi-block mesh.");
    }
   
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
    ensure_equals("Rings SI", HO.RingsSI(), 0);
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
    ensure_equals("Rings SI", HO.RingsSI(), 1);
    ensure_equals("RecOrder", HO.RecOrder(), 1);

    // check monotonicity containers
    std::vector<int> MonotonicityContainer(1); MonotonicityContainer[0] = OFF;
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
    ensure_equals("Rings SI", HO.RingsSI(), 1);
    ensure_equals("RecOrder", HO.RecOrder(), 1);

    // check monotonicity containers
    std::vector<int> MonotonicityContainer(1); MonotonicityContainer[0] = OFF;
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
    ensure_equals("Rings SI", HO.RingsSI(), 2);
    ensure_equals("RecOrder", HO.RecOrder(), 4);

    // check monotonicity containers
    std::vector<int> MonotonicityContainer(1); MonotonicityContainer[0] = OFF;
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
    ensure_equals("Rings SI", HO.RingsSI(), 1);
    ensure_equals("RecOrder", HO.RecOrder(), 4);

    // check monotonicity containers
    std::vector<int> MonotonicityContainer(1); MonotonicityContainer[0] = OFF;
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
    ensure_equals("Rings SI", HO.RingsSI(), 1);
    ensure_equals("RecOrder", HO.RecOrder(), 2);

    // check monotonicity containers
    std::vector<int> MonotonicityContainer(1); MonotonicityContainer[0] = OFF;
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
    ensure_equals("Rings SI", HO_Copy.RingsSI(), HO.RingsSI());

    for (i = 3; i<= 8; ++i) {
      for (j = 3; j<= 9; ++j) {
	ensure_equals("TD", HO_Copy.CellTaylorDeriv(i,j), HO.CellTaylorDeriv(i,j));
	ensure("LimitedCell", HO_Copy.CellInadequateFit(i,j) == HO.CellInadequateFit(i,j));
	ensure("SI", HO_Copy.CellSmoothnessIndicator(i,j) == HO.CellSmoothnessIndicator(i,j));
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
    ensure_equals("Rings SI", HO_Copy.RingsSI(), HO.RingsSI());

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
    ensure_equals("Rings SI", HO_Copy.RingsSI(), HO.RingsSI());

    for (i = 3; i<= 8; ++i) {
      for (j = 3; j<= 9; ++j) {
	ensure_equals("TD", HO_Copy.CellTaylorDeriv(i,j), HO.CellTaylorDeriv(i,j));
	ensure("LimitedCell", HO_Copy.CellInadequateFit(i,j) == HO.CellInadequateFit(i,j));
	ensure("SI", HO_Copy.CellSmoothnessIndicator(i,j) == HO.CellSmoothnessIndicator(i,j));
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
    Grid.BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux); Grid.BndNorthSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.BndEastSpline.setFluxCalcMethod(ReconstructionBasedFlux); Grid.BndEastSpline.setBCtype(BC_SYMMETRY_PLANE);
    
    // Update
    HO.AssociateGeometry(Grid);

    // == check that constrained reconstruction is required
    ensure_equals("Block Flag", HO.IsConstrainedReconstructionRequired(), true);
    ensure_equals("West Bnd Flag", HO.IsWestConstrainedReconstructionRequired(), false);
    ensure_equals("East Bnd Flag", HO.IsEastConstrainedReconstructionRequired(), true);
    ensure_equals("North Bnd Flag", HO.IsNorthConstrainedReconstructionRequired(), true);
    ensure_equals("South Bnd Flag", HO.IsSouthConstrainedReconstructionRequired(), false);
    // == check indexes for smoothness indicator
    ensure_equals("Rings SI", HO.RingsSI(), 2);

    // Change spline types again
    Grid.BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    
    // Update
    HO.AssociateGeometry(Grid);

    // == check that constrained reconstruction is required
    ensure_equals("Block Flag, II", HO.IsConstrainedReconstructionRequired(), true);
    ensure_equals("West Bnd Flag, II", HO.IsWestConstrainedReconstructionRequired(), true);
    ensure_equals("East Bnd Flag, II", HO.IsEastConstrainedReconstructionRequired(), false);
    ensure_equals("North Bnd Flag, II", HO.IsNorthConstrainedReconstructionRequired(), false);
    ensure_equals("South Bnd Flag, II", HO.IsSouthConstrainedReconstructionRequired(), true);
    // == check indexes for smoothness indicator
    ensure_equals("Rings SI, II", HO.RingsSI(), 2);
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

  /* Test 25:*/
  template<>
  template<>
  void HighOrder2D_object::test<25>()
  {
    set_test_name("Check loop indexes influenced by constraints, all sides constrained");
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

    // Change spline type
    Grid.BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux); Grid.BndNorthSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.BndEastSpline.setFluxCalcMethod(ReconstructionBasedFlux);  Grid.BndEastSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux); Grid.BndSouthSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux);  Grid.BndWestSpline.setBCtype(BC_SYMMETRY_PLANE);
    
    // Update
    HO.AssociateGeometry(Grid);

    // == check that constrained reconstruction is required
    ensure_equals("Block Flag", HO.IsConstrainedReconstructionRequired(), true);
    ensure_equals("West Bnd Flag", HO.IsWestConstrainedReconstructionRequired(), true);
    ensure_equals("East Bnd Flag", HO.IsEastConstrainedReconstructionRequired(), true);
    ensure_equals("North Bnd Flag", HO.IsNorthConstrainedReconstructionRequired(), true);
    ensure_equals("South Bnd Flag", HO.IsSouthConstrainedReconstructionRequired(), true);

    // == check indexes for smoothness indicator
    ensure_equals("Rings SI", HO.RingsSI(), 2);
  }

  /* Test 26:*/
  template<>
  template<>
  void HighOrder2D_object::test<26>()
  {
    set_test_name("Check loop indexes influenced by constraints, W and E sides constrained");
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

    // Change spline type
    Grid.BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.BndEastSpline.setFluxCalcMethod(ReconstructionBasedFlux); Grid.BndEastSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux); Grid.BndWestSpline.setBCtype(BC_SYMMETRY_PLANE);
    
    // Update
    HO.AssociateGeometry(Grid);

    // == check that constrained reconstruction is required
    ensure_equals("Block Flag", HO.IsConstrainedReconstructionRequired(), true);
    ensure_equals("West Bnd Flag", HO.IsWestConstrainedReconstructionRequired(), true);
    ensure_equals("East Bnd Flag", HO.IsEastConstrainedReconstructionRequired(), true);
    ensure_equals("North Bnd Flag", HO.IsNorthConstrainedReconstructionRequired(), false);
    ensure_equals("South Bnd Flag", HO.IsSouthConstrainedReconstructionRequired(), false);

    // == check indexes for smoothness indicator
    ensure_equals("Rings SI", HO.RingsSI(), 2);
  }

  /* Test 27:*/
  template<>
  template<>
  void HighOrder2D_object::test<27>()
  {
    set_test_name("Check loop indexes influenced by constraints, N and S sides constrained");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    HighOrder2D<double> HO;
    int RecOrder(3);
    
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

    // Change spline type
    Grid.BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    
    // Update
    HO.AssociateGeometry(Grid);

    // == check that constrained reconstruction is required
    ensure_equals("Block Flag", HO.IsConstrainedReconstructionRequired(), true);
    ensure_equals("West Bnd Flag", HO.IsWestConstrainedReconstructionRequired(), false);
    ensure_equals("East Bnd Flag", HO.IsEastConstrainedReconstructionRequired(), false);
    ensure_equals("North Bnd Flag", HO.IsNorthConstrainedReconstructionRequired(), true);
    ensure_equals("South Bnd Flag", HO.IsSouthConstrainedReconstructionRequired(), true);

    // == check indexes for smoothness indicator
    ensure_equals("Rings SI", HO.RingsSI(), 2);
  }

  /* Test 28:*/
  template<>
  template<>
  void HighOrder2D_object::test<28>()
  {
    set_test_name("Check loop indexes influenced by constraints, W and E sides constrained, SI with first neighbours");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    HighOrder2D<double> HO;
    int RecOrder(3);
    
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

    // Change spline type
    Grid.BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.BndEastSpline.setFluxCalcMethod(ReconstructionBasedFlux); Grid.BndEastSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux); Grid.BndWestSpline.setBCtype(BC_SYMMETRY_PLANE);
    
    // Update
    HO.AssociateGeometry(Grid);

    // == check that constrained reconstruction is required
    ensure_equals("Block Flag", HO.IsConstrainedReconstructionRequired(), true);
    ensure_equals("West Bnd Flag", HO.IsWestConstrainedReconstructionRequired(), true);
    ensure_equals("East Bnd Flag", HO.IsEastConstrainedReconstructionRequired(), true);
    ensure_equals("North Bnd Flag", HO.IsNorthConstrainedReconstructionRequired(), false);
    ensure_equals("South Bnd Flag", HO.IsSouthConstrainedReconstructionRequired(), false);

    // == check indexes for smoothness indicator
    ensure_equals("Rings SI", HO.RingsSI(), 1);
  }

  /* Test 29:*/
  template<>
  template<>
  void HighOrder2D_object::test<29>()
  {
    set_test_name("Check loop indexes influenced by constraints, N and S sides constrained, SI with first neighbours");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    HighOrder2D<double> HO;
    int RecOrder(3);
    
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

    // Change spline type
    Grid.BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux);  Grid.BndNorthSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux);  Grid.BndSouthSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    
    // Update
    HO.AssociateGeometry(Grid);

    // == check that constrained reconstruction is required
    ensure_equals("Block Flag", HO.IsConstrainedReconstructionRequired(), true);
    ensure_equals("West Bnd Flag", HO.IsWestConstrainedReconstructionRequired(), false);
    ensure_equals("East Bnd Flag", HO.IsEastConstrainedReconstructionRequired(), false);
    ensure_equals("North Bnd Flag", HO.IsNorthConstrainedReconstructionRequired(), true);
    ensure_equals("South Bnd Flag", HO.IsSouthConstrainedReconstructionRequired(), true);

    
    ensure_equals("West Bnd Reconstruction Flag", HO.getWestBnd().IsReconstructionConstrained(), false);
    ensure_equals("S_West Bnd Reconstruction Flag", HO.get_South_WestBnd().IsReconstructionConstrained(), false);
    ensure_equals("N_West Bnd Reconstruction Flag", HO.get_North_WestBnd().IsReconstructionConstrained(), false);

    ensure_equals("East Bnd Reconstruction Flag", HO.getEastBnd().IsReconstructionConstrained(), false);
    ensure_equals("S_East Bnd Reconstruction Flag", HO.get_South_EastBnd().IsReconstructionConstrained(), false);
    ensure_equals("N_East Bnd Reconstruction Flag", HO.get_North_EastBnd().IsReconstructionConstrained(), false);

    ensure_equals("North Bnd Reconstruction Flag", HO.getNorthBnd().IsReconstructionConstrained(), true);
    ensure_equals("E_North Bnd Reconstruction Flag", HO.get_East_NorthBnd().IsReconstructionConstrained(), false);
    ensure_equals("W_North Bnd Reconstruction Flag", HO.get_West_NorthBnd().IsReconstructionConstrained(), false);

    ensure_equals("South Bnd Reconstruction Flag", HO.getSouthBnd().IsReconstructionConstrained(), true);
    ensure_equals("E_South Bnd Reconstruction Flag", HO.get_East_SouthBnd().IsReconstructionConstrained(), false);
    ensure_equals("W_South Bnd Reconstruction Flag", HO.get_West_SouthBnd().IsReconstructionConstrained(), false);

    // == check indexes for smoothness indicator
    ensure_equals("Rings SI", HO.RingsSI(), 1);
  }

  /* Test 30:*/
  template<>
  template<>
  void HighOrder2D_object::test<30>()
  {
    set_test_name("Check consistency rule violation for constrained reconstruction");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    HighOrder2D<double> HO;
    int RecOrder(3);
    
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

    // Change spline type
    Grid.BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.BndEastSpline.setFluxCalcMethod(ReconstructionBasedFlux); Grid.BndEastSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux); Grid.BndWestSpline.setBCtype(BC_SYMMETRY_PLANE);

    // Force rule violation
    Grid.ExtendEast_BndNorthSpline = Grid.BndEastSpline;
    
    try {
      // Update
      HO.AssociateGeometry(Grid);

      // Ensure test failure
      fail("Check Rule #1 violation for East boundary");

    } catch (runtime_error){
      // Successful
    }
  }

  /* Test 31:*/
  template<>
  template<>
  void HighOrder2D_object::test<31>()
  {
    set_test_name("Check correct setup of block boundaries I");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    HighOrder2D<double> HO;
    int RecOrder(3);
    
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

    // Change spline type (i.e. create scenario)
    Grid.BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendEast_BndNorthSpline = Grid.BndNorthSpline;
    Grid.ExtendEast_BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendWest_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);

    Grid.BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndSouthSpline = Grid.BndSouthSpline;
    Grid.ExtendEast_BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux); 

    Grid.BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndWestSpline = Grid.BndWestSpline;
    Grid.ExtendSouth_BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux); 
    Grid.ExtendNorth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem); 
    
    // Update
    HO.AssociateGeometry(Grid);

    // == check that the right setup is obtained
    ensure_equals("Block Flag", HO.IsConstrainedReconstructionRequired(), true);

    // === West Bnd 
    ensure_equals("West Bnd Reconstruction Flag", HO.getWestBnd().IsReconstructionConstrained(), false);
    ensure_equals("West Bnd Stencil Flag", HO.getWestBnd().IsReconstructionStencilAffected(), false);
    
    ensure_equals("S_West Bnd Reconstruction Flag", HO.get_South_WestBnd().IsReconstructionConstrained(), true);
    ensure_equals("S_West Bnd Stencil Flag", HO.get_South_WestBnd().IsReconstructionStencilAffected(), true);

    ensure_equals("N_West Bnd Reconstruction Flag", HO.get_North_WestBnd().IsReconstructionConstrained(), false);
    ensure_equals("N_West Bnd Stencil Flag", HO.get_North_WestBnd().IsReconstructionStencilAffected(), true);

    // === East Bnd
    ensure_equals("East Bnd Reconstruction Flag", HO.getEastBnd().IsReconstructionConstrained(), false);
    ensure_equals("East Bnd Stencil Flag", HO.getEastBnd().IsReconstructionStencilAffected(), false);

    ensure_equals("S_East Bnd Reconstruction Flag", HO.get_South_EastBnd().IsReconstructionConstrained(), false);
    ensure_equals("S_East Bnd Stencil Flag", HO.get_South_EastBnd().IsReconstructionStencilAffected(), true);

    ensure_equals("N_East Bnd Reconstruction Flag", HO.get_North_EastBnd().IsReconstructionConstrained(), false);
    ensure_equals("N_East Bnd Stencil Flag", HO.get_North_EastBnd().IsReconstructionStencilAffected(), true);

    // === North Bnd
    ensure_equals("North Bnd Reconstruction Flag", HO.getNorthBnd().IsReconstructionConstrained(), true);
    ensure_equals("North Bnd Stencil Flag", HO.getNorthBnd().IsReconstructionStencilAffected(), true);

    ensure_equals("E_North Bnd Reconstruction Flag", HO.get_East_NorthBnd().IsReconstructionConstrained(), true);
    ensure_equals("E_North Bnd Stencil Flag", HO.get_East_NorthBnd().IsReconstructionStencilAffected(), true);

    ensure_equals("W_North Bnd Reconstruction Flag", HO.get_West_NorthBnd().IsReconstructionConstrained(), false);
    ensure_equals("W_North Bnd Stencil Flag", HO.get_West_NorthBnd().IsReconstructionStencilAffected(), false);

    // === South Bnd
    ensure_equals("South Bnd Reconstruction Flag", HO.getSouthBnd().IsReconstructionConstrained(), false);
    ensure_equals("South Bnd Stencil Flag", HO.getSouthBnd().IsReconstructionStencilAffected(), false);

    ensure_equals("E_South Bnd Reconstruction Flag", HO.get_East_SouthBnd().IsReconstructionConstrained(), true);
    ensure_equals("E_South Bnd Stencil Flag", HO.get_East_SouthBnd().IsReconstructionStencilAffected(), true);

    ensure_equals("W_South Bnd Reconstruction Flag", HO.get_West_SouthBnd().IsReconstructionConstrained(), false);
    ensure_equals("W_South Bnd Stencil Flag", HO.get_West_SouthBnd().IsReconstructionStencilAffected(), true);
  }

  /* Test 32:*/
  template<>
  template<>
  void HighOrder2D_object::test<32>()
  {
    set_test_name("Check correct setup of block boundaries II");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    HighOrder2D<double> HO;
    int RecOrder(3);
    
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

    // Change spline type (i.e. create scenario)
    Grid.BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendEast_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndEastSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.BndEastSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.ExtendNorth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);

    Grid.BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendWest_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendSouth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem); 
    Grid.ExtendNorth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem); 
    
    // Update
    HO.AssociateGeometry(Grid);

    // == check that the right setup is obtained
    ensure_equals("Block Flag", HO.IsConstrainedReconstructionRequired(), true);

    // === West Bnd 
    ensure_equals("West Bnd Reconstruction Flag", HO.getWestBnd().IsReconstructionConstrained(), true);
    ensure_equals("West Bnd Stencil Flag", HO.getWestBnd().IsReconstructionStencilAffected(), true);
    
    ensure_equals("S_West Bnd Reconstruction Flag", HO.get_South_WestBnd().IsReconstructionConstrained(), false);
    ensure_equals("S_West Bnd Stencil Flag", HO.get_South_WestBnd().IsReconstructionStencilAffected(), true);

    ensure_equals("N_West Bnd Reconstruction Flag", HO.get_North_WestBnd().IsReconstructionConstrained(), false);
    ensure_equals("N_West Bnd Stencil Flag", HO.get_North_WestBnd().IsReconstructionStencilAffected(), true);

    // === East Bnd
    ensure_equals("East Bnd Reconstruction Flag", HO.getEastBnd().IsReconstructionConstrained(), true);
    ensure_equals("East Bnd Stencil Flag", HO.getEastBnd().IsReconstructionStencilAffected(), true);

    ensure_equals("S_East Bnd Reconstruction Flag", HO.get_South_EastBnd().IsReconstructionConstrained(), false);
    ensure_equals("S_East Bnd Stencil Flag", HO.get_South_EastBnd().IsReconstructionStencilAffected(), true);

    ensure_equals("N_East Bnd Reconstruction Flag", HO.get_North_EastBnd().IsReconstructionConstrained(), false);
    ensure_equals("N_East Bnd Stencil Flag", HO.get_North_EastBnd().IsReconstructionStencilAffected(), true);

    // === North Bnd
    ensure_equals("North Bnd Reconstruction Flag", HO.getNorthBnd().IsReconstructionConstrained(), true);
    ensure_equals("North Bnd Stencil Flag", HO.getNorthBnd().IsReconstructionStencilAffected(), true);

    ensure_equals("E_North Bnd Reconstruction Flag", HO.get_East_NorthBnd().IsReconstructionConstrained(), false);
    ensure_equals("E_North Bnd Stencil Flag", HO.get_East_NorthBnd().IsReconstructionStencilAffected(), true);

    ensure_equals("W_North Bnd Reconstruction Flag", HO.get_West_NorthBnd().IsReconstructionConstrained(), false);
    ensure_equals("W_North Bnd Stencil Flag", HO.get_West_NorthBnd().IsReconstructionStencilAffected(), true);

    // === South Bnd
    ensure_equals("South Bnd Reconstruction Flag", HO.getSouthBnd().IsReconstructionConstrained(), true);
    ensure_equals("South Bnd Stencil Flag", HO.getSouthBnd().IsReconstructionStencilAffected(), true);

    ensure_equals("E_South Bnd Reconstruction Flag", HO.get_East_SouthBnd().IsReconstructionConstrained(), false);
    ensure_equals("E_South Bnd Stencil Flag", HO.get_East_SouthBnd().IsReconstructionStencilAffected(), true);

    ensure_equals("W_South Bnd Reconstruction Flag", HO.get_West_SouthBnd().IsReconstructionConstrained(), false);
    ensure_equals("W_South Bnd Stencil Flag", HO.get_West_SouthBnd().IsReconstructionStencilAffected(), true);
  }

  /* Test 33:*/
  template<>
  template<>
  void HighOrder2D_object::test<33>()
  {
    set_test_name("Check stencil setting for complex opaque boundaries");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    RunRegression = ON;

    HighOrder2D<double> HO;
    int RecOrder(3);
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;
    CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL = ON;

    // Generate a geometry
    Grid2D_Quad_Block_HO Grid;

    // Read the geometry from input file
    Open_Input_File("CartesianMesh.dat");
    in() >> Grid;

    // Initialize high-order variable
    HO.InitializeVariable(RecOrder,Grid,true);

    // Change spline type (i.e. create scenario)
    Grid.BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendEast_BndNorthSpline = Grid.BndNorthSpline;
    Grid.ExtendEast_BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendWest_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);

    Grid.BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndSouthSpline = Grid.BndSouthSpline;
    Grid.ExtendEast_BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux); 

    Grid.BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndWestSpline = Grid.BndWestSpline;
    Grid.ExtendSouth_BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux); 
    Grid.ExtendNorth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem); 
    
    // Update
    HO.AssociateGeometry(Grid);

    int iCell, jCell, i, j;
    IndexType i_index, j_index;
    rings = HO.Rings();

    MasterFile = "ReconstructionStencilSetupForComplexOpaqueBoundaryConfiguration.dat";
    CurrentFile = "Current_ReconstructionStencilSetupForComplexOpaqueBoundaryConfiguration.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);      
      out() << "Reconstruction type map for complex configuration\n";
      HO.outputReconstructionTypeMap(out());

      out() << endl 
	    << "Stencils for each cell between ICl-NghostHO, ICu+NghostHO, JCl-NghostHO, JCu+NghostHO\n"
	    << "Warning! In the cells with 'n' reconstruction the stencils are probably wrong! That's okay!\n";

      for (iCell = Grid.ICl-HO.NghostHO(); iCell<=Grid.ICu+HO.NghostHO(); ++iCell){
	for (jCell = Grid.JCl-HO.NghostHO(); jCell<=Grid.JCu+HO.NghostHO(); ++jCell){
	  HO.displayDeviatedReconstructionStencil(out(), iCell, jCell, rings);
	}
      }

      // == Compare current file against master
      RunRegressionTest("Reconstruction Map + Cell Stencils", CurrentFile, MasterFile, 1.0e-12);

    } else {
      // Generate the master file
      Open_Output_File(MasterFile);

      out() << "Reconstruction type map for complex configuration\n";
      HO.outputReconstructionTypeMap(out());

      out() << endl 
	    << "Stencils for each cell between ICl-NghostHO, ICu+NghostHO, JCl-NghostHO, JCu+NghostHO\n"
	    << "Warning! In the cells with 'n' reconstruction the stencils are probably wrong! That's okay!\n";

      for (iCell = Grid.ICl-HO.NghostHO(); iCell<=Grid.ICu+HO.NghostHO(); ++iCell){
	for (jCell = Grid.JCl-HO.NghostHO(); jCell<=Grid.JCu+HO.NghostHO(); ++jCell){
	  HO.displayDeviatedReconstructionStencil(out(), iCell, jCell, rings);
	}
      }
    }
  }

  /* Test 34:*/
  template<>
  template<>
  void HighOrder2D_object::test<34>()
  {
    set_test_name("Check reconstruction map for unconstrained block");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    HighOrder2D<double> HO;
    int RecOrder(3);
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;
    CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL = ON;

    // Generate a geometry
    Grid2D_Quad_Block_HO Grid;

    // Read the geometry from input file
    Open_Input_File("CartesianMesh.dat");
    in() >> Grid;

    // Initialize high-order variable
    HO.InitializeVariable(RecOrder,Grid,true);

    // Change spline type (i.e. create scenario)
    Grid.BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.BndEastSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.ExtendSouth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);

    Grid.BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem); 
    
    // Update
    HO.AssociateGeometry(Grid);

    ensure("Map variable", HO.getReconstructionTypeMap() == NULL);
  }

  /* Test 35:*/
  template<>
  template<>
  void HighOrder2D_object::test<35>()
  {
    set_test_name("Check reconstruction type map, All main boundaries constrained");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    RunRegression = ON;

    MasterFile = "ReconstructionTypeMap_AllBoundariesConstrained.dat";
    CurrentFile = "Current_ReconstructionTypeMap_AllBoundariesConstrained.dat";

    HighOrder2D<double> HO;
    int RecOrder(3);
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;
    CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL = ON;

    // Generate a geometry
    Grid2D_Quad_Block_HO Grid;

    // Read the geometry from input file
    Open_Input_File("CartesianMesh.dat");
    in() >> Grid;

    // Initialize high-order variable
    HO.InitializeVariable(RecOrder,Grid,true);

    // Change spline type (i.e. create scenario)
    Grid.BndEastSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.BndEastSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.ExtendSouth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendWest_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);

    Grid.BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendNorth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendWest_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem); 
    
    // Update
    HO.AssociateGeometry(Grid);

    if (RunRegression){
      Open_Output_File(CurrentFile);      
      out() << "Reconstruction type map all boundaries constrained\n";
      HO.outputReconstructionTypeMap(out());

      // == Check all setups
      RunRegressionTest("Reconstruction type map", CurrentFile, MasterFile, 1.0e-12);
    } else {
      // Generate the master file
      Open_Output_File(MasterFile);
      out() << "Reconstruction type map all boundaries constrained\n";
      HO.outputReconstructionTypeMap(out());
    }
  }

  /* Test 36:*/
  template<>
  template<>
  void HighOrder2D_object::test<36>()
  {
    set_test_name("Check reconstruction map with complicated setup");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    RunRegression = ON;

    MasterFile = "ReconstructionTypeMap_1.dat";
    CurrentFile = "Current_ReconstructionTypeMap_1.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);     
    } else {
      // Open the master file
      Open_Output_File(MasterFile);
    }

    HighOrder2D<double> HO;
    int RecOrder(3);
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;
    CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL = ON;

    // Generate a geometry
    Grid2D_Quad_Block_HO Grid;

    // Read the geometry from input file
    Open_Input_File("CartesianMesh.dat");
    in() >> Grid;

    // Initialize high-order variable
    HO.InitializeVariable(RecOrder,Grid,true);

    // ============= First Setup =================

    // Change spline type (i.e. create scenario)
    Grid.BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendEast_BndNorthSpline = Grid.BndNorthSpline;
    Grid.ExtendEast_BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendWest_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);

    Grid.BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndSouthSpline = Grid.BndSouthSpline;
    Grid.ExtendEast_BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux); 

    Grid.BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndWestSpline = Grid.BndWestSpline;
    Grid.ExtendSouth_BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux); 
    Grid.ExtendNorth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem); 
    
    // Update
    HO.AssociateGeometry(Grid);

    if (RunRegression){
      out() << "First setup\n";
      HO.outputReconstructionTypeMap(out());
    } else {
      // Generate the master file
      out() << "First setup\n";
      HO.outputReconstructionTypeMap(out());
    }


    // ============= Second Setup (Rotate first setup counterclockwise) =================
    // Change spline type (i.e. create scenario)
    Grid.BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendNorth_BndWestSpline = Grid.BndWestSpline;
    Grid.ExtendNorth_BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendSouth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);

    Grid.BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndEastSpline = Grid.BndEastSpline;
    Grid.ExtendNorth_BndEastSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.ExtendNorth_BndEastSpline.setFluxCalcMethod(ReconstructionBasedFlux); 

    Grid.BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndSouthSpline = Grid.BndSouthSpline;
    Grid.ExtendEast_BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux); 
    Grid.ExtendWest_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem); 

    // Update
    HO.AssociateGeometry(Grid);

    if (RunRegression){
      out() << "Second setup (Rotate first setup counterclockwise)\n";
      HO.outputReconstructionTypeMap(out());
    } else {
      // Generate the master file
      out() << "Second setup (Rotate first setup counterclockwise)\n";
      HO.outputReconstructionTypeMap(out());
    }


    // ============= Third Setup (Rotate second setup counterclockwise) =================

    // Change spline type (i.e. create scenario)
    Grid.BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendWest_BndSouthSpline = Grid.BndWestSpline;
    Grid.ExtendWest_BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendEast_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);

    Grid.BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndNorthSpline = Grid.BndNorthSpline;
    Grid.ExtendWest_BndNorthSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.ExtendWest_BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux); 

    Grid.BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndEastSpline = Grid.BndEastSpline;
    Grid.ExtendNorth_BndEastSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.ExtendNorth_BndEastSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendSouth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem); 

    // Update
    HO.AssociateGeometry(Grid);

    if (RunRegression){
      out() << "Third setup (Rotate second setup counterclockwise)\n";
      HO.outputReconstructionTypeMap(out());
    } else {
      // Generate the master file
      out() << "Third setup (Rotate second setup counterclockwise)\n";
      HO.outputReconstructionTypeMap(out());
    }


    // ============= Fourth Setup (Rotate third setup counterclockwise) =================

    // Change spline type (i.e. create scenario)
    Grid.BndEastSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.BndEastSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.ExtendSouth_BndEastSpline = Grid.BndEastSpline;
    Grid.ExtendSouth_BndEastSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendNorth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);

    Grid.BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndWestSpline = Grid.BndWestSpline;
    Grid.ExtendSouth_BndWestSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.ExtendSouth_BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux); 

    Grid.BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndNorthSpline = Grid.BndNorthSpline;
    Grid.ExtendWest_BndNorthSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.ExtendWest_BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendEast_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem); 

    // Update
    HO.AssociateGeometry(Grid);

    if (RunRegression){
      out() << "Fourth setup (Rotate third setup counterclockwise)\n";
      HO.outputReconstructionTypeMap(out());

      // == Check all setups
      RunRegressionTest("Reconstruction type map", CurrentFile, MasterFile, 1.0e-12);
    } else {
      // Generate the master file
      out() << "Fourth setup (Rotate third setup counterclockwise)\n";
      HO.outputReconstructionTypeMap(out());
    }
  }

  /* Test 37:*/
  template<>
  template<>
  void HighOrder2D_object::test<37>()
  {
    set_test_name("Check reconstruction map, One main spline and extensions constrained at a time");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    RunRegression = ON;

    MasterFile = "ReconstructionTypeMap_2.dat";
    CurrentFile = "Current_ReconstructionTypeMap_2.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);     
    } else {
      // Open the master file
      Open_Output_File(MasterFile);
    }

    HighOrder2D<double> HO;
    int RecOrder(3);
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;
    CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL = ON;

    // Generate a geometry
    Grid2D_Quad_Block_HO Grid;

    // Read the geometry from input file
    Open_Input_File("CartesianMesh.dat");
    in() >> Grid;

    // Initialize high-order variable
    HO.InitializeVariable(RecOrder,Grid,true);

    // ====================== EAST side constrained ===========================

    // Change spline type (i.e. create scenario)
    Grid.BndEastSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.BndEastSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.ExtendSouth_BndEastSpline = Grid.BndEastSpline;
    Grid.ExtendNorth_BndEastSpline = Grid.BndEastSpline;

    Grid.BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);

    Grid.BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem); 
    
    // Update
    HO.AssociateGeometry(Grid);

    if (RunRegression){
      out() << "East side constrained\n";
      HO.outputReconstructionTypeMap(out());
    } else {
      // Generate the master file
      out() << "East side constrained\n";
      HO.outputReconstructionTypeMap(out());
    }

    // ====================== SOUTH side constrained ===========================

    // Change spline type (i.e. create scenario)
    Grid.BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndEastSpline = Grid.BndEastSpline;
    Grid.ExtendNorth_BndEastSpline = Grid.BndEastSpline;

    Grid.BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendWest_BndSouthSpline = Grid.BndSouthSpline;
    Grid.ExtendEast_BndSouthSpline = Grid.BndSouthSpline;

    Grid.BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem); 
    
    // Update
    HO.AssociateGeometry(Grid);

    if (RunRegression){
      out() << "South side constrained\n";
      HO.outputReconstructionTypeMap(out());
    } else {
      // Generate the master file
      out() << "South side constrained\n";
      HO.outputReconstructionTypeMap(out());
    }

    // ====================== WEST side constrained ===========================

    // Change spline type (i.e. create scenario)
    Grid.BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndEastSpline = Grid.BndEastSpline;
    Grid.ExtendNorth_BndEastSpline = Grid.BndEastSpline;

    Grid.BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndSouthSpline = Grid.BndSouthSpline;
    Grid.ExtendEast_BndSouthSpline = Grid.BndSouthSpline;

    Grid.BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendNorth_BndWestSpline = Grid.BndWestSpline;
    Grid.ExtendSouth_BndWestSpline = Grid.BndWestSpline;

    Grid.BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem); 
    
    // Update
    HO.AssociateGeometry(Grid);

    if (RunRegression){
      out() << "West side constrained\n";
      HO.outputReconstructionTypeMap(out());
    } else {
      // Generate the master file
      out() << "West side constrained\n";
      HO.outputReconstructionTypeMap(out());
    }

    // ====================== NORTH side constrained ===========================

    // Change spline type (i.e. create scenario)
    Grid.BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndEastSpline = Grid.BndEastSpline;
    Grid.ExtendNorth_BndEastSpline = Grid.BndEastSpline;

    Grid.BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndSouthSpline = Grid.BndSouthSpline;
    Grid.ExtendEast_BndSouthSpline = Grid.BndSouthSpline;

    Grid.BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndWestSpline = Grid.BndWestSpline;
    Grid.ExtendSouth_BndWestSpline = Grid.BndWestSpline;

    Grid.BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendWest_BndNorthSpline = Grid.BndNorthSpline;
    Grid.ExtendEast_BndNorthSpline = Grid.BndNorthSpline;
    
    // Update
    HO.AssociateGeometry(Grid);

    if (RunRegression){
      out() << "North side constrained\n";
      HO.outputReconstructionTypeMap(out());

      // == Check all setups
      RunRegressionTest("Reconstruction type map", CurrentFile, MasterFile, 1.0e-12);
    } else {
      // Generate the master file
      out() << "North side constrained\n";
      HO.outputReconstructionTypeMap(out());
    }
  }

  /* Test 38:*/
  template<>
  template<>
  void HighOrder2D_object::test<38>()
  {
    set_test_name("Check reconstruction map, Two main splines that meet in a corner constrained at a time");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    RunRegression = ON;

    MasterFile = "ReconstructionTypeMap_3.dat";
    CurrentFile = "Current_ReconstructionTypeMap_3.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);     
    } else {
      // Open the master file
      Open_Output_File(MasterFile);
    }

    HighOrder2D<double> HO;
    int RecOrder(3);
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;
    CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL = ON;

    // Generate a geometry
    Grid2D_Quad_Block_HO Grid;

    // Read the geometry from input file
    Open_Input_File("CartesianMesh.dat");
    in() >> Grid;

    // Initialize high-order variable
    HO.InitializeVariable(RecOrder,Grid,true);

    // ====================== SOUTH-EAST corner constrained ===========================

    // Change spline type (i.e. create scenario)
    Grid.BndEastSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.BndEastSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.ExtendSouth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendWest_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);

    Grid.BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem); 
    
    // Update
    HO.AssociateGeometry(Grid);

    if (RunRegression){
      out() << "South-East corner constrained\n";
      HO.outputReconstructionTypeMap(out());
    } else {
      // Generate the master file
      out() << "South-East corner constrained\n";
      HO.outputReconstructionTypeMap(out());
    }

    // ====================== SOUTH-WEST corner constrained ===========================

    // Change spline type (i.e. create scenario)
    Grid.BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.BndEastSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.ExtendSouth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendWest_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);

    Grid.BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendNorth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem); 
    
    // Update
    HO.AssociateGeometry(Grid);

    if (RunRegression){
      out() << "South-West corner constrained\n";
      HO.outputReconstructionTypeMap(out());
    } else {
      // Generate the master file
      out() << "South-West corner constrained\n";
      HO.outputReconstructionTypeMap(out());
    }

    // ====================== NORTH-WEST corner constrained ===========================

    // Change spline type (i.e. create scenario)
    Grid.BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.BndEastSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.ExtendSouth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);

    Grid.BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendNorth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendWest_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem); 
    
    // Update
    HO.AssociateGeometry(Grid);

    if (RunRegression){
      out() << "North-West corner constrained\n";
      HO.outputReconstructionTypeMap(out());
    } else {
      // Generate the master file
      out() << "North-West corner constrained\n";
      HO.outputReconstructionTypeMap(out());
    }

    // ====================== NORTH-EAST corner constrained ===========================

    // Change spline type (i.e. create scenario)
    Grid.BndEastSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.BndEastSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.ExtendSouth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);

    Grid.BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendWest_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem); 
    
    // Update
    HO.AssociateGeometry(Grid);

    if (RunRegression){
      out() << "North-East corner constrained\n";
      HO.outputReconstructionTypeMap(out());

      // == Check all setups
      RunRegressionTest("Reconstruction type map", CurrentFile, MasterFile, 1.0e-12);
    } else {
      // Generate the master file
      out() << "North-East corner constrained\n";
      HO.outputReconstructionTypeMap(out());
    }
  }

  /* Test 39:*/
  template<>
  template<>
  void HighOrder2D_object::test<39>()
  {
    set_test_name("Check stencil setting for complex opaque boundaries, Second setup");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    RunRegression = ON;

    HighOrder2D<double> HO;
    int RecOrder(3);
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;
    CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL = ON;

    // Generate a geometry
    Grid2D_Quad_Block_HO Grid;

    // Read the geometry from input file
    Open_Input_File("CartesianMesh.dat");
    in() >> Grid;

    // Initialize high-order variable
    HO.InitializeVariable(RecOrder,Grid,true);

    // Change spline type (i.e. create scenario)
    Grid.BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendNorth_BndWestSpline = Grid.BndWestSpline;
    Grid.ExtendNorth_BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendSouth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);

    Grid.BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndEastSpline = Grid.BndEastSpline;
    Grid.ExtendNorth_BndEastSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.ExtendNorth_BndEastSpline.setFluxCalcMethod(ReconstructionBasedFlux); 

    Grid.BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndSouthSpline = Grid.BndSouthSpline;
    Grid.ExtendEast_BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux); 
    Grid.ExtendWest_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem); 
    
    // Update
    HO.AssociateGeometry(Grid);

    int iCell, jCell, i, j;
    IndexType i_index, j_index;
    rings = HO.Rings();

    MasterFile = "ReconstructionStencilSetupForComplexOpaqueBoundaryConfiguration_II.dat";
    CurrentFile = "Current_ReconstructionStencilSetupForComplexOpaqueBoundaryConfiguration_II.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);      
      out() << "Reconstruction type map for complex configuration\n";
      HO.outputReconstructionTypeMap(out());

      out() << endl 
	    << "Stencils for each cell between ICl-NghostHO, ICu+NghostHO, JCl-NghostHO, JCu+NghostHO\n"
	    << "Warning! In the cells with 'n' reconstruction the stencils are probably wrong! That's okay!\n";

      for (iCell = Grid.ICl-HO.NghostHO(); iCell<=Grid.ICu+HO.NghostHO(); ++iCell){
	for (jCell = Grid.JCl-HO.NghostHO(); jCell<=Grid.JCu+HO.NghostHO(); ++jCell){
	  HO.displayDeviatedReconstructionStencil(out(), iCell, jCell, rings);
	}
      }

      // == Compare current file against master
      RunRegressionTest("Reconstruction Map + Cell Stencils", CurrentFile, MasterFile, 1.0e-12);

    } else {
      // Generate the master file
      Open_Output_File(MasterFile);

      out() << "Reconstruction type map for complex configuration\n";
      HO.outputReconstructionTypeMap(out());

      out() << endl 
	    << "Stencils for each cell between ICl-NghostHO, ICu+NghostHO, JCl-NghostHO, JCu+NghostHO\n"
	    << "Warning! In the cells with 'n' reconstruction the stencils are probably wrong! That's okay!\n";

      for (iCell = Grid.ICl-HO.NghostHO(); iCell<=Grid.ICu+HO.NghostHO(); ++iCell){
	for (jCell = Grid.JCl-HO.NghostHO(); jCell<=Grid.JCu+HO.NghostHO(); ++jCell){
	  HO.displayDeviatedReconstructionStencil(out(), iCell, jCell, rings);
	}
      }
    }
  }

  /* Test 40:*/
  template<>
  template<>
  void HighOrder2D_object::test<40>()
  {
    set_test_name("Check stencil setting for complex opaque boundaries, Third setup");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    RunRegression = ON;

    HighOrder2D<double> HO;
    int RecOrder(3);
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;
    CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL = ON;

    // Generate a geometry
    Grid2D_Quad_Block_HO Grid;

    // Read the geometry from input file
    Open_Input_File("CartesianMesh.dat");
    in() >> Grid;

    // Initialize high-order variable
    HO.InitializeVariable(RecOrder,Grid,true);

    // Change spline type (i.e. create scenario)
    Grid.BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendWest_BndSouthSpline = Grid.BndWestSpline;
    Grid.ExtendWest_BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendEast_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);

    Grid.BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndNorthSpline = Grid.BndNorthSpline;
    Grid.ExtendWest_BndNorthSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.ExtendWest_BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux); 

    Grid.BndEastSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndEastSpline = Grid.BndEastSpline;
    Grid.ExtendNorth_BndEastSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.ExtendNorth_BndEastSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendSouth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem); 
    
    // Update
    HO.AssociateGeometry(Grid);

    int iCell, jCell, i, j;
    IndexType i_index, j_index;
    rings = HO.Rings();

    MasterFile = "ReconstructionStencilSetupForComplexOpaqueBoundaryConfiguration_III.dat";
    CurrentFile = "Current_ReconstructionStencilSetupForComplexOpaqueBoundaryConfiguration_III.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);      
      out() << "Reconstruction type map for complex configuration\n";
      HO.outputReconstructionTypeMap(out());

      out() << endl 
	    << "Stencils for each cell between ICl-NghostHO, ICu+NghostHO, JCl-NghostHO, JCu+NghostHO\n"
	    << "Warning! In the cells with 'n' reconstruction the stencils are probably wrong! That's okay!\n";

      for (iCell = Grid.ICl-HO.NghostHO(); iCell<=Grid.ICu+HO.NghostHO(); ++iCell){
	for (jCell = Grid.JCl-HO.NghostHO(); jCell<=Grid.JCu+HO.NghostHO(); ++jCell){
	  HO.displayDeviatedReconstructionStencil(out(), iCell, jCell, rings);
	}
      }

      // == Compare current file against master
      RunRegressionTest("Reconstruction Map + Cell Stencils", CurrentFile, MasterFile, 1.0e-12);

    } else {
      // Generate the master file
      Open_Output_File(MasterFile);

      out() << "Reconstruction type map for complex configuration\n";
      HO.outputReconstructionTypeMap(out());

      out() << endl 
	    << "Stencils for each cell between ICl-NghostHO, ICu+NghostHO, JCl-NghostHO, JCu+NghostHO\n"
	    << "Warning! In the cells with 'n' reconstruction the stencils are probably wrong! That's okay!\n";

      for (iCell = Grid.ICl-HO.NghostHO(); iCell<=Grid.ICu+HO.NghostHO(); ++iCell){
	for (jCell = Grid.JCl-HO.NghostHO(); jCell<=Grid.JCu+HO.NghostHO(); ++jCell){
	  HO.displayDeviatedReconstructionStencil(out(), iCell, jCell, rings);
	}
      }
    }
  }

  /* Test 41:*/
  template<>
  template<>
  void HighOrder2D_object::test<41>()
  {
    set_test_name("Check stencil setting for complex opaque boundaries, Fourth setup");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    RunRegression = ON;

    HighOrder2D<double> HO;
    int RecOrder(3);
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;
    CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL = ON;

    // Generate a geometry
    Grid2D_Quad_Block_HO Grid;

    // Read the geometry from input file
    Open_Input_File("CartesianMesh.dat");
    in() >> Grid;

    // Initialize high-order variable
    HO.InitializeVariable(RecOrder,Grid,true);

    // Change spline type (i.e. create scenario)
    Grid.BndEastSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.BndEastSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.ExtendSouth_BndEastSpline = Grid.BndEastSpline;
    Grid.ExtendSouth_BndEastSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendNorth_BndEastSpline.setFluxCalcMethod(SolveRiemannProblem); 

    Grid.BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendEast_BndSouthSpline.setFluxCalcMethod(SolveRiemannProblem);

    Grid.BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendNorth_BndWestSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendSouth_BndWestSpline = Grid.BndWestSpline;
    Grid.ExtendSouth_BndWestSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.ExtendSouth_BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux); 

    Grid.BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem);
    Grid.ExtendWest_BndNorthSpline = Grid.BndNorthSpline;
    Grid.ExtendWest_BndNorthSpline.setBCtype(BC_SYMMETRY_PLANE);
    Grid.ExtendWest_BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Grid.ExtendEast_BndNorthSpline.setFluxCalcMethod(SolveRiemannProblem); 
    
    // Update
    HO.AssociateGeometry(Grid);

    int iCell, jCell, i, j;
    IndexType i_index, j_index;
    rings = HO.Rings();

    MasterFile = "ReconstructionStencilSetupForComplexOpaqueBoundaryConfiguration_IV.dat";
    CurrentFile = "Current_ReconstructionStencilSetupForComplexOpaqueBoundaryConfiguration_IV.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);      
      out() << "Reconstruction type map for complex configuration\n";
      HO.outputReconstructionTypeMap(out());

      out() << endl 
	    << "Stencils for each cell between ICl-NghostHO, ICu+NghostHO, JCl-NghostHO, JCu+NghostHO\n"
	    << "Warning! In the cells with 'n' reconstruction the stencils are probably wrong! That's okay!\n";

      for (iCell = Grid.ICl-HO.NghostHO(); iCell<=Grid.ICu+HO.NghostHO(); ++iCell){
	for (jCell = Grid.JCl-HO.NghostHO(); jCell<=Grid.JCu+HO.NghostHO(); ++jCell){
	  HO.displayDeviatedReconstructionStencil(out(), iCell, jCell, rings);
	}
      }

      // == Compare current file against master
      RunRegressionTest("Reconstruction Map + Cell Stencils", CurrentFile, MasterFile, 1.0e-12);

    } else {
      // Generate the master file
      Open_Output_File(MasterFile);

      out() << "Reconstruction type map for complex configuration\n";
      HO.outputReconstructionTypeMap(out());

      out() << endl 
	    << "Stencils for each cell between ICl-NghostHO, ICu+NghostHO, JCl-NghostHO, JCu+NghostHO\n"
	    << "Warning! In the cells with 'n' reconstruction the stencils are probably wrong! That's okay!\n";

      for (iCell = Grid.ICl-HO.NghostHO(); iCell<=Grid.ICu+HO.NghostHO(); ++iCell){
	for (jCell = Grid.JCl-HO.NghostHO(); jCell<=Grid.JCu+HO.NghostHO(); ++jCell){
	  HO.displayDeviatedReconstructionStencil(out(), iCell, jCell, rings);
	}
      }
    }
  }

  /* Test 42:*/
  template<>
  template<>
  void HighOrder2D_object::test<42>()
  {
    set_test_name("Check stencil setup for multi-block mesh with complex opaque boundaries");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    // Local variables
    HighOrder2D<double> HO[9];
    int RecOrder(4);
    int iCell, jCell, i, j, index;
    IndexType i_index, j_index;
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;
    CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL = ON;
    
    // Read geometry
    Grid2D_Quad_MultiBlock_HO Mesh;
    Grid2DTesting_Input_Parameters IP;
    strcpy(IP.Grid_File_Name, "HighOrderReconstruction/UnitTests/HighOrder2D_Data/CartesianMesh_3x3.grid");
    strcpy(IP.Output_File_Name, "StencilSettingMesh");
    Mesh.Read_Multi_Block_Grid_Using_IP(IP);

    // Set splines properties to create the desired scenario in Block (1,1)
    Mesh(1,1).BndNorthSpline = Mesh(1,2).BndNorthSpline;
    Mesh(1,1).BndNorthSpline.setBCtype(BC_SYMMETRY_PLANE);
    Mesh(1,1).BndNorthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Mesh(1,2).BndSouthSpline = Mesh(1,1).BndNorthSpline;    
    Mesh(0,2).ExtendEast_BndSouthSpline = Mesh(1,1).BndNorthSpline;
    Mesh(0,1).ExtendEast_BndNorthSpline = Mesh(1,1).BndNorthSpline;

    Mesh(1,1).ExtendEast_BndNorthSpline = Mesh(1,1).BndNorthSpline;
    Mesh(1,1).ExtendEast_BndNorthSpline.setBCtype(BC_SYMMETRY_PLANE);
    Mesh(2,2).BndSouthSpline = Mesh(1,1).ExtendEast_BndNorthSpline;
    Mesh(2,2).ExtendWest_BndSouthSpline = Mesh(1,1).BndNorthSpline;
    Mesh(2,1).BndNorthSpline = Mesh(1,1).ExtendEast_BndNorthSpline;
    Mesh(2,1).ExtendWest_BndNorthSpline = Mesh(1,1).BndNorthSpline;
    Mesh(1,2).ExtendEast_BndSouthSpline = Mesh(1,1).ExtendEast_BndNorthSpline;

    Mesh(1,1).ExtendSouth_BndWestSpline = Mesh(0,0).BndWestSpline;
    Mesh(1,1).ExtendSouth_BndWestSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Mesh(1,1).ExtendSouth_BndWestSpline.setBCtype(BC_SYMMETRY_PLANE);
    Mesh(0,0).BndEastSpline = Mesh(1,1).ExtendSouth_BndWestSpline;
    Mesh(1,0).BndWestSpline = Mesh(1,1).ExtendSouth_BndWestSpline;
    Mesh(0,1).ExtendSouth_BndEastSpline = Mesh(1,1).ExtendSouth_BndWestSpline;

    Mesh(1,1).ExtendEast_BndSouthSpline = Mesh(2,0).BndSouthSpline;
    Mesh(1,1).ExtendEast_BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux);
    Mesh(1,1).ExtendEast_BndSouthSpline.setBCtype(BC_SYMMETRY_PLANE);
    Mesh(2,1).BndSouthSpline = Mesh(1,1).ExtendEast_BndSouthSpline;
    Mesh(2,0).BndNorthSpline = Mesh(1,1).ExtendEast_BndSouthSpline;
    Mesh(1,0).ExtendEast_BndNorthSpline = Mesh(1,1).ExtendEast_BndSouthSpline;

    // Initialize high-order variables
    for (j=0, index=0; j<3; ++j){
      for (i = 0; i<3; ++i, ++index){
	HO[index].InitializeVariable(RecOrder,Mesh(i,j),true);
      }
    }

    // == Check stencils
    int ICl, ICu, JCl, JCu;
    ICl = 5;
    ICu = 14;
    JCl = 5;
    JCu = 24;

    // Set stencil rings
    rings = HO[0].Rings();

    CheckHighOrderReconstructionStencilConsistency(HO[4],
						   ICu+1, ICu+3,
						   JCl  , JCu,
						   HO[5],
						   ICl, JCl,
						   "HO[4] East against HO[5]");
    CheckHighOrderReconstructionStencilConsistency(HO[5],
						   ICl-1, ICl-3,
						   JCl  , JCu,
						   HO[4],
						   ICu, JCl,
						   "HO[5] West against HO[4]");


    CheckHighOrderReconstructionStencilConsistency(HO[4],
						   ICl  , ICl+2,
						   JCl-1, JCl-3,
						   HO[1],
						   ICl, JCu,
						   "HO[4] South against HO[1]");
    CheckHighOrderReconstructionStencilConsistency(HO[1],
						   ICl  , ICl+2,
						   JCu+1, JCu+3,
						   HO[4],
						   ICl, JCl,
						   "HO[1] North against HO[4]");


    CheckHighOrderReconstructionStencilConsistency(HO[0],
						   ICu  , ICu-2,
						   JCu+1, JCu+3,
						   HO[3],
						   ICu, JCl,
						   "HO[0] North against HO[3]");
    CheckHighOrderReconstructionStencilConsistency(HO[3],
						   ICu  , ICu-2,
						   JCl-1, JCl-3,
						   HO[0],
						   ICu, JCu,
						   "HO[3] South against HO[0]");


    CheckHighOrderReconstructionStencilConsistency(HO[3],
						   ICu+1, ICu+3,
						   JCu  , JCu-2,
						   HO[4],
						   ICl, JCu,
						   "HO[3] East against HO[4]");
    CheckHighOrderReconstructionStencilConsistency(HO[4],
						   ICl-1, ICl-3,
						   JCu, JCu-2,
						   HO[3],
						   ICu, JCu,
						   "HO[4] West against HO[3]");


    CheckHighOrderReconstructionStencilConsistency(HO[6],
						   ICu+1, ICu+3,
						   JCl  , JCl+2,
						   HO[7],
						   ICl, JCl,
						   "HO[6] East against HO[7]");
    CheckHighOrderReconstructionStencilConsistency(HO[7],
						   ICl-1, ICl-3,
						   JCl, JCl+2,
						   HO[6],
						   ICu, JCl,
						   "HO[7] West against HO[6]");


    CheckHighOrderReconstructionStencilConsistency(HO[7],
						   ICu+1, ICu+3,
						   JCl  , JCl+2,
						   HO[8],
						   ICl, JCl,
						   "HO[7] East against HO[8]");
    CheckHighOrderReconstructionStencilConsistency(HO[8],
						   ICl-1, ICl-3,
						   JCl, JCl+2,
						   HO[7],
						   ICu, JCl,
						   "HO[8] West against HO[7]");


    CheckHighOrderReconstructionStencilConsistency(HO[1],
						   ICu+1, ICu+3,
						   JCu  , JCu-2,
						   HO[2],
						   ICl, JCu,
						   "HO[1] East against HO[2]");
    CheckHighOrderReconstructionStencilConsistency(HO[2],
						   ICl-1, ICl-3,
						   JCu, JCu-2,
						   HO[1],
						   ICu, JCu,
						   "HO[2] West against HO[1]");


    // Set stencil rings
    rings = 1;

    CheckHighOrderReconstructionStencilConsistency(HO[4],
						   ICu+1, ICu+2,
						   JCl  , JCu,
						   HO[5],
						   ICl, JCl,
						   "HO[4] East against HO[5], 1 ring");
    CheckHighOrderReconstructionStencilConsistency(HO[5],
						   ICl-1, ICl-2,
						   JCl  , JCu,
						   HO[4],
						   ICu, JCl,
						   "HO[5] West against HO[4], 1 ring");


    CheckHighOrderReconstructionStencilConsistency(HO[4],
						   ICl  , ICl+1,
						   JCl-1, JCl-2,
						   HO[1],
						   ICl, JCu,
						   "HO[4] South against HO[1], 1 ring");
    CheckHighOrderReconstructionStencilConsistency(HO[1],
						   ICl  , ICl+1,
						   JCu+1, JCu+2,
						   HO[4],
						   ICl, JCl,
						   "HO[1] North against HO[4], 1 ring");


    CheckHighOrderReconstructionStencilConsistency(HO[0],
						   ICu  , ICu-1,
						   JCu+1, JCu+2,
						   HO[3],
						   ICu, JCl,
						   "HO[0] North against HO[3], 1 ring");
    CheckHighOrderReconstructionStencilConsistency(HO[3],
						   ICu  , ICu-1,
						   JCl-1, JCl-2,
						   HO[0],
						   ICu, JCu,
						   "HO[3] South against HO[0], 1 ring");


    CheckHighOrderReconstructionStencilConsistency(HO[3],
						   ICu+1, ICu+2,
						   JCu  , JCu-1,
						   HO[4],
						   ICl, JCu,
						   "HO[3] East against HO[4], 1 ring");
    CheckHighOrderReconstructionStencilConsistency(HO[4],
						   ICl-1, ICl-2,
						   JCu, JCu-1,
						   HO[3],
						   ICu, JCu,
						   "HO[4] West against HO[3], 1 ring");


    CheckHighOrderReconstructionStencilConsistency(HO[6],
						   ICu+1, ICu+2,
						   JCl  , JCl+1,
						   HO[7],
						   ICl, JCl,
						   "HO[6] East against HO[7], 1 ring");
    CheckHighOrderReconstructionStencilConsistency(HO[7],
						   ICl-1, ICl-2,
						   JCl, JCl+1,
						   HO[6],
						   ICu, JCl,
						   "HO[7] West against HO[6], 1 ring");


    CheckHighOrderReconstructionStencilConsistency(HO[7],
						   ICu+1, ICu+2,
						   JCl  , JCl+1,
						   HO[8],
						   ICl, JCl,
						   "HO[7] East against HO[8], 1 ring");
    CheckHighOrderReconstructionStencilConsistency(HO[8],
						   ICl-1, ICl-2,
						   JCl, JCl+1,
						   HO[7],
						   ICu, JCl,
						   "HO[8] West against HO[7], 1 ring");


    CheckHighOrderReconstructionStencilConsistency(HO[1],
						   ICu+1, ICu+2,
						   JCu  , JCu-1,
						   HO[2],
						   ICl, JCu,
						   "HO[1] East against HO[2], 1 ring");
    CheckHighOrderReconstructionStencilConsistency(HO[2],
						   ICl-1, ICl-2,
						   JCu, JCu-1,
						   HO[1],
						   ICu, JCu,
						   "HO[2] West against HO[1], 1 ring");
  }


  /* Test 43:*/
  template<>
  template<>
  void HighOrder2D_object::test<43>()
  {
    set_test_name("Calculate normal gradient");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    HighOrder2D<double> HO;
    int RecOrder(4);
    int p1,p2;
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = ON;

    // Generate a geometry
    Grid2D_Quad_Block_HO Grid;

    // Read the geometry from input file
    Open_Input_File("CartesianMesh.dat");
    in() >> Grid;

    // Initialize high-order variable
    HO.InitializeVariable(RecOrder,Grid,false);

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

    // == check solution
    double Result = 0.59509999999999985;
    Vector2D Point(-1.2, -1.8);	// This corresponds to DeltaX = 0.1 and DeltaY = -0.5

    ensure_distance("Point value with cell (3,3)",
		    HO.SolutionAtCoordinates(3,3,Point.x, Point.y), Result, AcceptedError(Result));

    ensure_distance("Point value with cell (3,3)",
		    HO.SolutionAtCoordinates(3,3,Point.x, Point.y, 1), Result, AcceptedError(Result));

    

    // === check gradients ===
    double dRdx(0), dRdy(0);
    double DeltaX(0.1), DeltaY(-0.5);

    // Calculate analytical results
    for ( p1 = 0; p1 <= RecOrder; ++p1){
      for ( p2 = 0; p2 <= RecOrder; ++p2){
	if (p1 + p2 <= RecOrder){
	  // == compute the x-derivative
	  dRdx += p1*pow(DeltaX,p1-1)*pow(DeltaY,p2)*HO.CellTaylorDerivState(3,3,p1,p2);

	  // == compute the y-derivative
	  dRdy += p2*pow(DeltaX,p1)*pow(DeltaY,p2-1)*HO.CellTaylorDerivState(3,3,p1,p2);
	}
      }
    }

    ensure_distance("X-Gradient value with cell (3,3)",
		    HO.XGradientStateAtLocation(3,3,Point), dRdx, AcceptedError(dRdx));

    ensure_distance("Y-Gradient value with cell (3,3)",
		    HO.YGradientStateAtLocation(3,3,Point), dRdy, AcceptedError(dRdy,1.0e-13));


    // === check gradient in the normal direction
    Vector2D normal_dir(-0.5, - 0.5*sqrt(3.0));
    int parameter(1);
    double dRdn;

    dRdn = HO.XGradientStateAtLocation(3,3,Point) * normal_dir.x +  HO.YGradientStateAtLocation(3,3,Point) * normal_dir.y;

    ensure_distance("Normal Gradient value with cell (3,3)",
		    HO.NormalGradientStateAtCoordinates(3,3,Point.x,Point.y,normal_dir),
		    dRdn, AcceptedError(dRdn));
    ensure_distance("Normal Gradient value with cell (3,3)",
		    HO.NormalGradientAtCoordinates(3,3,Point.x,Point.y,normal_dir, parameter),
		    dRdn, AcceptedError(dRdn));
    ensure_distance("Normal Gradient value with cell (3,3)",
		    HO.NormalGradientStateAtLocation(3,3,Point,normal_dir),
		    dRdn, AcceptedError(dRdn));
    ensure_distance("Normal Gradient value with cell (3,3)",
		    HO.NormalGradientAtLocation(3,3,Point,normal_dir, parameter),
		    dRdn, AcceptedError(dRdn));
  }


  /* Test 44:*/
  template<>
  template<>
  void HighOrder2D_object::test<44>()
  {
    set_test_name("Calculate solution prolongation");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    HighOrder2D<double> HO;
    int RecOrder(4);
    int p1,p2;
    double Ubar, U_SW, U_SE, U_NE, U_NW;
    double Ufine_SW, Ufine_SE, Ufine_NE, Ufine_NW;
    double areaSW, areaSE, areaNE, areaNW;
    Node2D_HO nodeSW,nodeSE,nodeNW,nodeNE,MidN,MidS,MidE,MidW,CC;
    int iCell,jCell;
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = ON;

    // Generate a geometry
    Grid2D_Quad_Block_HO Grid;

    // Read the geometry from input file
    Open_Input_File("CartesianMesh.dat");
    in() >> Grid;

    // Initialize high-order variable
    HO.InitializeVariable(RecOrder,Grid,false);

    // Set cell indexes and node coordinates
    iCell = 3;
    jCell = 3;
    nodeSW = Grid.nodeSW(iCell,jCell);
    nodeSE = Grid.nodeSE(iCell,jCell);
    nodeNW = Grid.nodeNW(iCell,jCell);
    nodeNE = Grid.nodeNE(iCell,jCell);
    MidN.setloc(Grid.xfaceN(iCell,jCell));
    MidS.setloc(Grid.xfaceS(iCell,jCell));
    MidE.setloc(Grid.xfaceE(iCell,jCell));
    MidW.setloc(Grid.xfaceW(iCell,jCell));
    CC.setloc(Grid.getNodeAverage(iCell,jCell));

    // Set the Taylor derivatives
    HO.CellTaylorDerivState(iCell,jCell,0,0) = 1.0;
    HO.CellTaylorDerivState(iCell,jCell,0,1) = 2.0;
    HO.CellTaylorDerivState(iCell,jCell,0,2) = 3.0;
    HO.CellTaylorDerivState(iCell,jCell,0,3) = 4.0;
    HO.CellTaylorDerivState(iCell,jCell,0,4) = 5.0;
    HO.CellTaylorDerivState(iCell,jCell,1,0) = 1.0;
    HO.CellTaylorDerivState(iCell,jCell,1,1) = 2.0;
    HO.CellTaylorDerivState(iCell,jCell,1,2) = 3.0;
    HO.CellTaylorDerivState(iCell,jCell,1,3) = 4.0;
    HO.CellTaylorDerivState(iCell,jCell,2,0) = 1.0;
    HO.CellTaylorDerivState(iCell,jCell,2,1) = 2.0;
    HO.CellTaylorDerivState(iCell,jCell,2,2) = 3.0;
    HO.CellTaylorDerivState(iCell,jCell,3,0) = 1.0;
    HO.CellTaylorDerivState(iCell,jCell,3,1) = 2.0;
    HO.CellTaylorDerivState(iCell,jCell,4,0) = 1.0;

    // Calculate average solutions for each sector and (iCell,jCell) cell
    Ubar = ( HO.IntegrateCellReconstructionOverQuadDomain(iCell,jCell,
							  nodeSW, nodeNW, 
							  nodeNE, nodeSE)/ Grid.CellArea(iCell,jCell) );

    areaSW = Grid.area(nodeSW,MidW,CC,MidS);
    U_SW = ( HO.IntegrateCellReconstructionOverQuadDomain(iCell,jCell,
							  nodeSW,MidW,CC,MidS)/ areaSW );

    areaSE = Grid.area(MidS,CC,MidE,nodeSE);
    U_SE = ( HO.IntegrateCellReconstructionOverQuadDomain(iCell,jCell,
							  MidS,CC,MidE,nodeSE)/ areaSE );

    areaNW = Grid.area(MidW,nodeNW,MidN,CC);
    U_NW = ( HO.IntegrateCellReconstructionOverQuadDomain(iCell,jCell,
							  MidW,nodeNW,MidN,CC)/ areaNW );

    areaNE = Grid.area(CC,MidN,nodeNE,MidE);
    U_NE = ( HO.IntegrateCellReconstructionOverQuadDomain(iCell,jCell,
							  CC,MidN,nodeNE,MidE)/ areaNE );

    // Calculate the average values for each sector provided by the high-order prolongation
    HO.ComputeHighOrderSolutionProlongation(iCell,jCell, Ubar,
					    Ufine_SW, Ufine_NW, 
					    Ufine_SE, Ufine_NE);

    // === Check results ===
    ensure_distance("Prolonged average solution to the SW sector", Ufine_SW, U_SW, AcceptedError(U_SW));
    ensure_distance("Prolonged average solution to the SE sector", Ufine_SE, U_SE, AcceptedError(U_SE));
    ensure_distance("Prolonged average solution to the NW sector", Ufine_SW, U_SW, AcceptedError(U_NW));
    ensure_distance("Prolonged average solution to the NE sector", Ufine_SW, U_SW, AcceptedError(U_NE));
  }

  /* Test 45:*/
  template<>
  template<>
  void HighOrder2D_object::test<45>()
  {
    set_test_name("Calculate solution prolongation with disturbed mesh");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    HighOrder2D<double> HO;
    int RecOrder(4);
    int p1,p2;
    double Ubar, U_SW, U_SE, U_NE, U_NW;
    double Ufine_SW, Ufine_SE, Ufine_NE, Ufine_NW;
    double areaSW, areaSE, areaNE, areaNW;
    Node2D_HO nodeSW,nodeSE,nodeNW,nodeNE,MidN,MidS,MidE,MidW,CC;
    int iCell,jCell;
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = ON;

    // Generate a geometry
    Grid2D_Quad_Block_HO Grid;

    // Read the geometry from input file
    Open_Input_File("CartesianMesh.dat");
    in() >> Grid;
    
    // Disturb mesh
    Grid.Disturb_Interior_Nodes(100);

    // Initialize high-order variable
    HO.InitializeVariable(RecOrder,Grid,false);

    for (iCell = Grid.ICl; iCell <= Grid.ICu; ++iCell){
      for (jCell = Grid.JCl; jCell <= Grid.JCu; ++jCell){

	// Set cell indexes and node coordinates
	nodeSW = Grid.nodeSW(iCell,jCell);
	nodeSE = Grid.nodeSE(iCell,jCell);
	nodeNW = Grid.nodeNW(iCell,jCell);
	nodeNE = Grid.nodeNE(iCell,jCell);
	MidN.setloc(Grid.xfaceN(iCell,jCell));
	MidS.setloc(Grid.xfaceS(iCell,jCell));
	MidE.setloc(Grid.xfaceE(iCell,jCell));
	MidW.setloc(Grid.xfaceW(iCell,jCell));
	CC.setloc(Grid.getNodeAverage(iCell,jCell));

	// Set the Taylor derivatives
	HO.CellTaylorDerivState(iCell,jCell,0,0) = 1.0;
	HO.CellTaylorDerivState(iCell,jCell,0,1) = 2.0;
	HO.CellTaylorDerivState(iCell,jCell,0,2) = 3.0;
	HO.CellTaylorDerivState(iCell,jCell,0,3) = 4.0;
	HO.CellTaylorDerivState(iCell,jCell,0,4) = 5.0;
	HO.CellTaylorDerivState(iCell,jCell,1,0) = 1.0;
	HO.CellTaylorDerivState(iCell,jCell,1,1) = 2.0;
	HO.CellTaylorDerivState(iCell,jCell,1,2) = 3.0;
	HO.CellTaylorDerivState(iCell,jCell,1,3) = 4.0;
	HO.CellTaylorDerivState(iCell,jCell,2,0) = 1.0;
	HO.CellTaylorDerivState(iCell,jCell,2,1) = 2.0;
	HO.CellTaylorDerivState(iCell,jCell,2,2) = 3.0;
	HO.CellTaylorDerivState(iCell,jCell,3,0) = 1.0;
	HO.CellTaylorDerivState(iCell,jCell,3,1) = 2.0;
	HO.CellTaylorDerivState(iCell,jCell,4,0) = 1.0;

	// Calculate average solutions for each sector and (iCell,jCell) cell
	Ubar = ( HO.IntegrateCellReconstructionOverQuadDomain(iCell,jCell,
							      nodeSW, nodeNW, 
							      nodeNE, nodeSE)/ Grid.CellArea(iCell,jCell) );

	areaSW = Grid.area(nodeSW,MidW,CC,MidS);
	U_SW = ( HO.IntegrateCellReconstructionOverQuadDomain(iCell,jCell,
							      nodeSW,MidW,CC,MidS)/ areaSW );

	areaSE = Grid.area(MidS,CC,MidE,nodeSE);
	U_SE = ( HO.IntegrateCellReconstructionOverQuadDomain(iCell,jCell,
							      MidS,CC,MidE,nodeSE)/ areaSE );

	areaNW = Grid.area(MidW,nodeNW,MidN,CC);
	U_NW = ( HO.IntegrateCellReconstructionOverQuadDomain(iCell,jCell,
							      MidW,nodeNW,MidN,CC)/ areaNW );

	areaNE = Grid.area(CC,MidN,nodeNE,MidE);
	U_NE = ( HO.IntegrateCellReconstructionOverQuadDomain(iCell,jCell,
							      CC,MidN,nodeNE,MidE)/ areaNE );

	// Calculate the average values for each sector provided by the high-order prolongation
	HO.ComputeHighOrderSolutionProlongation(iCell,jCell, Ubar,
						Ufine_SW, Ufine_NW, 
						Ufine_SE, Ufine_NE);
    
	// === Check results ===
	ostm() << "Checked Cell (" << iCell << "," << jCell << "), " 
	       << "prolonged solution to the SW sector\n";
	ensure_distance(ostm().str(), Ufine_SW, U_SW, AcceptedError(U_SW));
	ostmClear();

	ostm() << "Checked Cell (" << iCell << "," << jCell << "), " 
	       << "prolonged solution to the SE sector\n";
	ensure_distance(ostm().str(), Ufine_SE, U_SE, AcceptedError(U_SE));
	ostmClear();

	ostm() << "Checked Cell (" << iCell << "," << jCell << "), " 
	       << "prolonged solution to the NW sector\n";
	ensure_distance(ostm().str(), Ufine_SW, U_SW, AcceptedError(U_NW));
	ostmClear();

	ostm() << "Checked Cell (" << iCell << "," << jCell << "), " 
	       << "prolonged solution to the NE sector\n";
	ensure_distance(ostm().str(), Ufine_SW, U_SW, AcceptedError(U_NE));
	ostmClear();
      }
    }

  }

  /* Test 46:*/
  template<>
  template<>
  void HighOrder2D_object::test<46>()
  {
    set_test_name("Calculate solution prolongation with high-order circular mesh");
    set_local_input_path("HighOrder2D_Data");
    set_local_output_path("HighOrder2D_Data");

    HighOrder2D<double> HO;
    int RecOrder(4);
    int p1,p2;
    double Ubar, U_SW, U_SE, U_NE, U_NW;
    double Ufine_SW, Ufine_SE, Ufine_NE, Ufine_NW;
    double areaSW, areaSE, areaNE, areaNW;
    Node2D_HO nodeSW,nodeSE,nodeNW,nodeNE,MidN,MidS,MidE,MidW,CC;
    int iCell,jCell;
    
    // Set execution mode
    CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = ON;

    // Generate a geometry
    Grid2D_Quad_MultiBlock_HO MeshBlk;

    // Add test particular input parameters
    IP.i_Grid = GRID_CIRCULAR_CYLINDER;
    IP.Cylinder_Radius = 1;
    IP.Cylinder_Radius2 = 6;
    IP.Number_of_Blocks_Jdir = 1;
    IP.Number_of_Blocks_Idir = 2;
    IP.Number_of_Cells_Idir = 20;
    IP.Number_of_Cells_Jdir = 10;
    IP.Number_of_Ghost_Cells = 5;
    IP.Space_Accuracy = 4;
    IP.IncludeHighOrderBoundariesRepresentation = ON;
    IP.i_Smooth_Quad_Block = OFF;
    Grid2D_Quad_Block_HO::setContourIntegrationBasedOnLinearSegments();
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    IP.i_Mesh_Stretching = ON;
    IP.Mesh_Stretching_Type_Idir = STRETCHING_FCN_MINMAX_CLUSTERING;
    IP.Mesh_Stretching_Type_Jdir = STRETCHING_FCN_MIN_CLUSTERING;
    IP.Mesh_Stretching_Factor_Idir = 1.025;
    IP.Mesh_Stretching_Factor_Jdir = 1.001;
    IP.i_Reconstruction = RECONSTRUCTION_HIGH_ORDER;
    IP.X_Scale = 1.0;
    IP.X_Rotate = 0.0;
    IP.X_Shift = Vector2D(0);

    // Build the mesh
    CreateMesh(MeshBlk,IP);
    MeshBlk(0,0).BndSouthSpline.setFluxCalcMethod(ReconstructionBasedFlux);

    // Initialize high-order variable
    HO.InitializeVariable(RecOrder,MeshBlk(0,0),false);

    for (iCell = MeshBlk(0,0).ICl; iCell <= MeshBlk(0,0).ICu; ++iCell){
      for (jCell = MeshBlk(0,0).JCl; jCell <= MeshBlk(0,0).JCu; ++jCell){

	// Set cell indexes and node coordinates
	nodeSW = MeshBlk(0,0).nodeSW(iCell,jCell);
	nodeSE = MeshBlk(0,0).nodeSE(iCell,jCell);
	nodeNW = MeshBlk(0,0).nodeNW(iCell,jCell);
	nodeNE = MeshBlk(0,0).nodeNE(iCell,jCell);
	MidN.setloc(MeshBlk(0,0).xfaceN(iCell,jCell));
	MidS.setloc(MeshBlk(0,0).xfaceS(iCell,jCell));
	MidE.setloc(MeshBlk(0,0).xfaceE(iCell,jCell));
	MidW.setloc(MeshBlk(0,0).xfaceW(iCell,jCell));
	CC.setloc(MeshBlk(0,0).getNodeAverage(iCell,jCell));

	// Set the Taylor derivatives
	HO.CellTaylorDerivState(iCell,jCell,0,0) = 1.0;
	HO.CellTaylorDerivState(iCell,jCell,0,1) = 2.0;
	HO.CellTaylorDerivState(iCell,jCell,0,2) = 3.0;
	HO.CellTaylorDerivState(iCell,jCell,0,3) = 4.0;
	HO.CellTaylorDerivState(iCell,jCell,0,4) = 5.0;
	HO.CellTaylorDerivState(iCell,jCell,1,0) = 1.0;
	HO.CellTaylorDerivState(iCell,jCell,1,1) = 2.0;
	HO.CellTaylorDerivState(iCell,jCell,1,2) = 3.0;
	HO.CellTaylorDerivState(iCell,jCell,1,3) = 4.0;
	HO.CellTaylorDerivState(iCell,jCell,2,0) = 1.0;
	HO.CellTaylorDerivState(iCell,jCell,2,1) = 2.0;
	HO.CellTaylorDerivState(iCell,jCell,2,2) = 3.0;
	HO.CellTaylorDerivState(iCell,jCell,3,0) = 1.0;
	HO.CellTaylorDerivState(iCell,jCell,3,1) = 2.0;
	HO.CellTaylorDerivState(iCell,jCell,4,0) = 1.0;


	// Calculate average solutions for each sector and (iCell,jCell) cell
	// 	Ubar = ( HO.IntegrateCellReconstructionOverAnotherCellDomain(iCell,jCell,
	// 								     HO,
	// 								     iCell,jCell)/ MeshBlk(0,0).CellArea(iCell,jCell) );
	Ubar = 10.0;

	// Calculate the average values for each sector provided by the high-order prolongation
	HO.ComputeHighOrderSolutionProlongation(iCell,jCell, Ubar,
						Ufine_SW, Ufine_NW, 
						Ufine_SE, Ufine_NE);
      }
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

