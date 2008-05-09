/*!\file test_AdvectDiffuse2D_Simulations.cc
  \brief Regression tests for 2D advection-diffusion solver. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../AdvectDiffuse2DQuad.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_AdvectDiffuse2DSolver : public TestData {

    // Local variables
  public:
    int error_flag;

    // Member functions

    // Constructor
    Data_AdvectDiffuse2DSolver(){
      /* Set the global path for this test suite. 
	 It's a good practice to put it in the constructor of the data class in order to be set
	 automatically for each individual test. Declare it relative to the /src_2D directory,
	 otherwise the framework might not find the input and output files. */
      
      set_test_suite_path("AdvectDiffuse2D/UnitTests/");
    }

    ~Data_AdvectDiffuse2DSolver(){};

    // Solve Problem --> run AdvectDiffuse2DQuadSolver() subroutine with the proper input_file_name
    //                   and analyze the returned error_flag
    void Solve_Problem(int batch_flag = 1);

  private:
    
  };

  void Data_AdvectDiffuse2DSolver::Solve_Problem(int batch_flag){
    // Call solver
    error_flag = AdvectDiffuse2DQuadSolver(input_file_name,
					   batch_flag);
      
    // Check error_flag
    if (error_flag) {
      throw runtime_error("AdvectDiffuse2DQuadSolver() ERROR: Runtime error! For details, run verbose the regression test.");
    } /* endif */
  }


  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. 'NAME'_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_AdvectDiffuse2DSolver> AdvectDiffuse2DSolver_TestSuite;
  typedef AdvectDiffuse2DSolver_TestSuite::object AdvectDiffuse2DSolver_object;


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
  void AdvectDiffuse2DSolver_object::test<1>()
  {

    set_test_name("2nd-order advection-diffusion in rectangular channel");
    set_local_output_path("SimulationsData/AdvectionDiffusionInRectangularChannel/");
    set_local_input_path("SimulationsData/AdvectionDiffusionInRectangularChannel/");

    RunRegression = ON;
    
    // Set input file name
    Open_Input_File("AdvectionDiffusionInRectangularChannel.in");

    // call solver
    Solve_Problem();

    if (RunRegression){

      //===== Check nodal solution
      MasterFile  = "AdvectionDiffusionInRectangularChannel_cpu000000.dat";
      CurrentFile = "Current_AdvectionDiffusionInRectangularChannel_cpu000000.dat";
      // check
      RunRegressionTest("Nodal Solution", CurrentFile, MasterFile, 5.0e-9, 5.0e-9);

      //===== Check cell solution
      MasterFile  = "AdvectionDiffusionInRectangularChannel_cells_cpu000000.dat";
      CurrentFile = "Current_AdvectionDiffusionInRectangularChannel_cells_cpu000000.dat";
      // check
      RunRegressionTest("Cell Solution", CurrentFile, MasterFile, 5.0e-9, 5.0e-9);

      //===== Check solution error norms
      MasterFile  = "AdvectionDiffusionInRectangularChannel_ErrorNorms.dat";
      CurrentFile = "Current_AdvectionDiffusionInRectangularChannel_ErrorNorms.dat";
      // check
      RunRegressionTest("Solution Errors", CurrentFile, MasterFile, 5.0e-9, 5.0e-9);
    }
  }

  /* Test 2:*/
  template<>
  template<>
  void AdvectDiffuse2DSolver_object::test<2>()
  {

    set_test_name("2nd-order heat transfer in rectangular box");
    set_local_output_path("SimulationsData/HeatTransferInRectangularBox/");
    set_local_input_path("SimulationsData/HeatTransferInRectangularBox/");

    RunRegression = ON;
    
    // Set input file name
    Open_Input_File("StationaryHeatTransferWithLinearSource.in");

    // call solver
    Solve_Problem();

    if (RunRegression){

      //===== Check nodal solution 00
      MasterFile  = "StationaryHeatTransferWithLinearSource_00_cpu000000.dat";
      CurrentFile = "Current_StationaryHeatTransferWithLinearSource_00_cpu000000.dat";
      RunRegressionTest("Nodal solution 00", CurrentFile, MasterFile, 5.0e-9, 5.0e-9);

      //===== Check nodal solution 01
      MasterFile  = "StationaryHeatTransferWithLinearSource_01_cpu000000.dat";
      CurrentFile = "Current_StationaryHeatTransferWithLinearSource_01_cpu000000.dat";
      RunRegressionTest("Nodal solution 01", CurrentFile, MasterFile, 5.0e-9, 5.0e-9);

      //===== Check nodal solution 02
      MasterFile  = "StationaryHeatTransferWithLinearSource_02_cpu000000.dat";
      CurrentFile = "Current_StationaryHeatTransferWithLinearSource_02_cpu000000.dat";
      RunRegressionTest("Nodal solution 02", CurrentFile, MasterFile, 5.0e-9, 5.0e-9);

      //===== Check nodal solution 03
      MasterFile  = "StationaryHeatTransferWithLinearSource_03_cpu000000.dat";
      CurrentFile = "Current_StationaryHeatTransferWithLinearSource_03_cpu000000.dat";
      RunRegressionTest("Nodal solution 03", CurrentFile, MasterFile, 5.0e-9, 5.0e-9);

      //===== Check cell solution 00
      MasterFile  = "StationaryHeatTransferWithLinearSource_cells_00_cpu000000.dat";
      CurrentFile = "Current_StationaryHeatTransferWithLinearSource_cells_00_cpu000000.dat";
      RunRegressionTest("Cell solution 00", CurrentFile, MasterFile, 5.0e-9, 5.0e-9);

      //===== Check cell solution 01
      MasterFile  = "StationaryHeatTransferWithLinearSource_cells_01_cpu000000.dat";
      CurrentFile = "Current_StationaryHeatTransferWithLinearSource_cells_01_cpu000000.dat";
      RunRegressionTest("Cell solution 01", CurrentFile, MasterFile, 5.0e-9, 5.0e-9);

      //===== Check cell solution 02
      MasterFile  = "StationaryHeatTransferWithLinearSource_cells_02_cpu000000.dat";
      CurrentFile = "Current_StationaryHeatTransferWithLinearSource_cells_02_cpu000000.dat";
      RunRegressionTest("Cell solution 02", CurrentFile, MasterFile, 5.0e-9, 5.0e-9);

      //===== Check cell solution 03
      MasterFile  = "StationaryHeatTransferWithLinearSource_cells_03_cpu000000.dat";
      CurrentFile = "Current_StationaryHeatTransferWithLinearSource_cells_03_cpu000000.dat";
      RunRegressionTest("Cell solution 03", CurrentFile, MasterFile, 5.0e-9, 5.0e-9);

      //===== Check solution error norms
      MasterFile  = "StationaryHeatTransferWithLinearSource_ErrorNorms.dat";
      CurrentFile = "Current_StationaryHeatTransferWithLinearSource_ErrorNorms.dat";
      RunRegressionTest("Solution Errors", CurrentFile, MasterFile, 5.0e-9, 5.0e-9);
    }
  }

  /* Test 3:*/
  template<>
  template<>
  void AdvectDiffuse2DSolver_object::test<3>()
  {

    set_test_name("2nd-order periodic advection in rectangular channel");
    set_local_output_path("SimulationsData/PeriodicAdvectionInRectangularBox/");
    set_local_input_path("SimulationsData/PeriodicAdvectionInRectangularBox/");

    RunRegression = ON;
    
    // Set input file name
    Open_Input_File("PeriodicWaveAdvection.in");

    // call solver
    Solve_Problem();

    if (RunRegression){

      //===== Check nodal solution
      MasterFile  = "PeriodicWaveAdvection_cpu000000.dat";
      CurrentFile = "Current_PeriodicWaveAdvection_cpu000000.dat";
      // check
      RunRegressionTest("Nodal solution", CurrentFile, MasterFile, 5.0e-9, 5.0e-9);

      //===== Check cell solution
      MasterFile  = "PeriodicWaveAdvection_cells_cpu000000.dat";
      CurrentFile = "Current_PeriodicWaveAdvection_cells_cpu000000.dat";
      // check
      RunRegressionTest("Cell solution", CurrentFile, MasterFile, 5.0e-9, 5.0e-9);

      //===== Check solution erros
      MasterFile  = "PeriodicWaveAdvection_ErrorNorms.dat";
      CurrentFile = "Current_PeriodicWaveAdvection_ErrorNorms.dat";
      // check
      RunRegressionTest("Solution Errors", CurrentFile, MasterFile, 5.0e-9, 5.0e-9);
    }
  }

  /* Test 4:*/
  template<>
  template<>
  void AdvectDiffuse2DSolver_object::test<4>()
  {

    set_test_name("2nd-order circular advection in rectangular box");
    set_local_output_path("SimulationsData/CircularAdvectionInRectangularBox/");
    set_local_input_path("SimulationsData/CircularAdvectionInRectangularBox/");

    RunRegression = ON;
    
    // Set input file name
    Open_Input_File("CircularAdvectionInRectangularBox.in");

    // call solver
    Solve_Problem();

    if (RunRegression){

      //===== Check nodal solution
      MasterFile  = "CircularAdvectionInRectangularBox_cpu000000.dat";
      CurrentFile = "Current_CircularAdvectionInRectangularBox_cpu000000.dat";
      // check
      RunRegressionTest("Nodal solution", CurrentFile, MasterFile, 5.0e-3, 5.0e-9);

      //===== Check cell solution
      MasterFile  = "CircularAdvectionInRectangularBox_cells_cpu000000.dat";
      CurrentFile = "Current_CircularAdvectionInRectangularBox_cells_cpu000000.dat";
      // check
      RunRegressionTest("Cell solution", CurrentFile, MasterFile, 5.0e-9, 5.0e-9);

      //===== Check solution errors
      MasterFile  = "CircularAdvectionInRectangularBox_ErrorNorms.dat";
      CurrentFile = "Current_CircularAdvectionInRectangularBox_ErrorNorms.dat";
      // check
      RunRegressionTest("Solution Errors", CurrentFile, MasterFile, 5.0e-9, 5.0e-9);
    }
  }

}



// Test suite constructor
tut::AdvectDiffuse2DSolver_TestSuite AdvectDiffuse2DSolverTestSuite("Solver: AdvectDiffuse2D");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

