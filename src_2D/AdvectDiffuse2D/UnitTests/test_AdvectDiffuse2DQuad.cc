/*!\file test_AdvectDiffuse2DQuad.cc
  \brief Regression tests for class AdvectDiffuse2D_Quad_Block. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../AdvectDiffuse2DQuad.h"
#include "../Grid/HO_Grid2DQuadMultiBlock.h" /* Include 2D quadrilateral multiblock grid header file */
#include "../../HighOrderReconstruction/AccuracyAssessment2DMultiBlock.h"
#include "../../HighOrderReconstruction/HighOrder2D_MultiBlock.h" /* Include 2D high-order header file for multi-block level. */

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_AdvectDiffuse2D_Quad_Block : public TestData {

    // Local variables
  public:

    // ==== Member data =====

    // AdvectDiffuse2D input variables and parameters:
    AdvectDiffuse2D_Input_Parameters IP;
    
    /* Multi-block solution-adaptive quadrilateral mesh 
       solution variables. */
    Grid2D_Quad_MultiBlock_HO     MeshBlk;
    QuadTreeBlock_DataStructure  QuadTree;
    AdaptiveBlockResourceList    GlobalList_Soln_Blocks;
    AdaptiveBlock2D_List         LocalList_Soln_Blocks;
    AdvectDiffuse2D_Quad_Block   *SolnBlk;
    AdvectDiffuse2D_Quad_Block    Soln;
    CPUTime processor_cpu_time;

    int error_flag;
    int Status;                        // shows if the computational domain has been initialized

    //===== Member functions =====

    // Default Constructor
    Data_AdvectDiffuse2D_Quad_Block();

    // Destructor
    ~Data_AdvectDiffuse2D_Quad_Block(void);

    // Initialization of the computational domain
    void InitializeComputationalDomain(Grid2D_Quad_MultiBlock_HO & _MeshBlk_,
				       QuadTreeBlock_DataStructure & _QuadTree_,
				       AdaptiveBlockResourceList & _GlobalList_Soln_Blocks_,
				       AdaptiveBlock2D_List & _LocalList_Soln_Blocks_,
				       AdvectDiffuse2D_Quad_Block *& _SolnBlk_,
				       AdvectDiffuse2D_Input_Parameters & _IP_) throw(std::runtime_error);

    // Output_Block()
    void Output_Block(AdvectDiffuse2D_Quad_Block & SolnBlock,
		      AdvectDiffuse2D_Input_Parameters & _IP_){};

  private:
    
  };


  Data_AdvectDiffuse2D_Quad_Block::Data_AdvectDiffuse2D_Quad_Block(){
    /* Set the global path for this test suite. 
       It's a good practice to put it in the constructor of the data class in order to be set
       automatically for each individual test. Declare it relative to the /src_2D directory,
       otherwise the framework might not find the input and output files. */
    
    set_test_suite_path("AdvectDiffuse2D/UnitTests/");  

    // Initialize IP to default values
    Set_Default_Input_Parameters(IP);

    Status = OFF;
  }

  Data_AdvectDiffuse2D_Quad_Block::~Data_AdvectDiffuse2D_Quad_Block(void){
    if (Status == ON){
      SolnBlk = Deallocate(SolnBlk, IP);
      LocalList_Soln_Blocks.deallocate();
      GlobalList_Soln_Blocks.deallocate();
      QuadTree.deallocate();
    }
  }

  // === InitializeComputationalDomain()
  void Data_AdvectDiffuse2D_Quad_Block::InitializeComputationalDomain(Grid2D_Quad_MultiBlock_HO & _MeshBlk_,
								      QuadTreeBlock_DataStructure & _QuadTree_,
								      AdaptiveBlockResourceList & _GlobalList_Soln_Blocks_,
								      AdaptiveBlock2D_List & _LocalList_Soln_Blocks_,
								      AdvectDiffuse2D_Quad_Block *& _SolnBlk_,
								      AdvectDiffuse2D_Input_Parameters & _IP_)
    throw(std::runtime_error){

    error_flag = _MeshBlk_.Multi_Block_Grid(_IP_);
      
    if (error_flag) {
      throw runtime_error("CreateMesh() ERROR: Unable to create valid Euler2D multi-block mesh.");
    }

    _SolnBlk_ = Create_Initial_Solution_Blocks(_MeshBlk_.Grid_ptr,
					       _SolnBlk_,
					       _IP_,
					       _QuadTree_,
					       _GlobalList_Soln_Blocks_,
					       _LocalList_Soln_Blocks_);

    HighOrder2D_MultiBlock::Create_Initial_HighOrder_Variables(_SolnBlk_,
							       _LocalList_Soln_Blocks_);

    if (_SolnBlk_ == NULL) {
      throw runtime_error("Create_Initial_Solution_Blocks() ERROR: Unable to create initial Euler2D solution blocks.");
    }

    Status = ON;
  }

  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. 'NAME'_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_AdvectDiffuse2D_Quad_Block> AdvectDiffuse2D_Quad_Block_TestSuite;
  typedef AdvectDiffuse2D_Quad_Block_TestSuite::object AdvectDiffuse2D_Quad_Block_object;


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
  void AdvectDiffuse2D_Quad_Block_object::test<1>()
  {

    set_test_name("Generate advection-diffusion solution block");
    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    RunRegression = ON;
 
    // Set input file name
    Open_Input_File("CircularAdvectionDiffusion.in");

    // Parse the input file
    IP.Verbose() = false;
    IP.Parse_Input_File(input_file_name);

    // Create computational domain
    InitializeComputationalDomain(MeshBlk,QuadTree,
				  GlobalList_Soln_Blocks, LocalList_Soln_Blocks, 
				  SolnBlk, IP);

    // Apply initial condition
    ICs(SolnBlk,LocalList_Soln_Blocks,IP);

    // Output Tecplot
    Output_Tecplot(SolnBlk,LocalList_Soln_Blocks,IP,0,0);

    if (RunRegression){
      //===== Check solution
      MasterFile  = "Circular_Advection_In_Rectangular_Box_cpu000000.dat";
      CurrentFile = "Current_Circular_Advection_In_Rectangular_Box_cpu000000.dat";
      
      // check
      RunRegressionTest("Multiblock Tecplot Output", CurrentFile, MasterFile, 5.0e-7, 5.0e-7);
    }
  }

  /* Test 2:*/
  template<>
  template<>
  void AdvectDiffuse2D_Quad_Block_object::test<2>()
  {

    set_test_name("Input-output operators (good for restart!)");
    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    RunRegression = ON;

    // Set input file name
    Open_Input_File("CircularAdvectionDiffusion.in");

    // Parse the input file
    IP.Verbose() = false;
    IP.Parse_Input_File(input_file_name);

    // Create computational domain
    InitializeComputationalDomain(MeshBlk,QuadTree,
				  GlobalList_Soln_Blocks, LocalList_Soln_Blocks, 
				  SolnBlk, IP);

    // Apply initial condition
    ICs(SolnBlk,LocalList_Soln_Blocks,IP);

    // === check operator << 
    MasterFile  = "SolnBlk_blk000003.dat";
    CurrentFile = "Current_SolnBlk_blk000003.dat";

    if (RunRegression){     

      // Open file
      Open_Output_File(CurrentFile);

      // Print checked data to file
      Print_File(SolnBlk[3],out());

      // === check data
      RunRegressionTest("operator <<", CurrentFile, MasterFile, 5.0e-7, 5.0e-7);


      // === check input-output operator
      Check_Input_Output_Operator("Solution Block 4", SolnBlk[3]);

    } else {
      // produce master

      // Open file
      Open_Output_File(MasterFile);

      // Print block 4 to file
      Print_File(SolnBlk[3],out());
    }
  }

  /* Test 3:*/
  template<>
  template<>
  void AdvectDiffuse2D_Quad_Block_object::test<3>()
  {

    set_test_name("Assess accuracy");
    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    // Set input file name
    Open_Input_File("CircularAdvectionDiffusion.in");

    // Parse the input file
    IP.Verbose() = false;
    IP.Parse_Input_File(input_file_name);

    // Create computational domain
    InitializeComputationalDomain(MeshBlk,QuadTree,
				  GlobalList_Soln_Blocks, LocalList_Soln_Blocks, 
				  SolnBlk, IP);

    // Apply initial condition
    ICs(SolnBlk,LocalList_Soln_Blocks,IP);

    // Assess solution accuracy
    AccuracyAssessment2D_MultiBlock::PrintErrorNorms(SolnBlk,
 						     LocalList_Soln_Blocks,
 						     IP, cout);

    // == check errors
    ensure_equals("Number of cells", AccuracyAssessment2D_MultiBlock::TotalNumberOfCells(), 144);
    ensure_distance("L1_Norm", AccuracyAssessment2D_MultiBlock::L1(), 0.06550436942, AcceptedError(0.06550436942,1.0e-7) );
    ensure_distance("L2_Norm", AccuracyAssessment2D_MultiBlock::L2(), 0.1153247286, AcceptedError(0.1153247286,1.0e-7) );
    ensure_distance("LMax_Norm", AccuracyAssessment2D_MultiBlock::LMax(), 0.1961452823, AcceptedError(0.1961452823,1.0e-7) );
  }

  /* Test 4:*/
  template<>
  template<>
  void AdvectDiffuse2D_Quad_Block_object::test<4>()
  {

    set_test_name("Test InterfaceSolutionGradient() with diamond path reconstruction");
    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    // Set input file name
    Open_Input_File("LaplaceInRectangularBox.in");

    // Parse the input file
    IP.Verbose() = false;
    IP.Parse_Input_File(input_file_name);

    // Create computational domain
    InitializeComputationalDomain(MeshBlk,QuadTree,
				  GlobalList_Soln_Blocks, LocalList_Soln_Blocks, 
				  SolnBlk, IP);

    // Apply initial condition
    ICs(SolnBlk,LocalList_Soln_Blocks,IP);

    // Compute nodal values
    SolnBlk[0].Calculate_Nodal_Solutions();

    // Analytic result for the gradient at the interface between cell (7,6) and (8,6) (i-direction interface)
    Vector2D Location = SolnBlk[0].Grid.xfaceE(7,6);
    Vector2D Result = SolnBlk[0].ExactSoln->Gradient(Location.x,Location.y);

    ensure_distance("Gradient betwen cells (7,6) and (8,6)",
    		    SolnBlk[0].InterfaceSolutionGradient(7,6,8,6,
							 VISCOUS_RECONSTRUCTION_DIAMOND_PATH,
							 DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION),
    		    Result,
    		    AcceptedError(Result,1.0e-12));

    // Analytic result for the gradient at the interface between cell (7,6) and (7,7) (j-direction interface)
    Location = SolnBlk[0].Grid.xfaceN(7,6);
    Result = SolnBlk[0].ExactSoln->Gradient(Location.x,Location.y);

    ensure_distance("Gradient betwen cells (7,6) and (7,7)",
		    SolnBlk[0].InterfaceSolutionGradient(7,6,7,7,
							 VISCOUS_RECONSTRUCTION_DIAMOND_PATH,
							 DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION),
		    Result,
		    AcceptedError(Result,1.0e-12));
  }

  /* Test 5:*/
  template<>
  template<>
  void AdvectDiffuse2D_Quad_Block_object::test<5>()
  {

    set_test_name("EllipticFluxStateAtInteriorInterface()");
    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    // Set input file name
    Open_Input_File("LaplaceInRectangularBox.in");

    // Parse the input file
    IP.Verbose() = false;
    IP.Parse_Input_File(input_file_name);

    // Create computational domain
    InitializeComputationalDomain(MeshBlk,QuadTree,
				  GlobalList_Soln_Blocks, LocalList_Soln_Blocks, 
				  SolnBlk, IP);

    // Apply initial condition
    ICs(SolnBlk,LocalList_Soln_Blocks,IP);

    // Compute nodal values
    SolnBlk[0].Calculate_Nodal_Solutions();

    // Analytic result for the solution at the interface between cell (7,6) and (8,6) (i-direction interface)
    Vector2D Location = SolnBlk[0].Grid.xfaceE(7,6);
    double Result = SolnBlk[0].ExactSoln->Solution(Location.x,Location.y);

    // Numeric result
    AdvectDiffuse2D_State U_face;
    SolnBlk[0].EllipticFluxStateAtInteriorInterface(7,6,8,6,Location,U_face);

    ensure_distance("Solution betwen cells (7,6) and (8,6)", U_face.u, Result, AcceptedError(Result,1.0e-8));

    // Analytic result for the solution at the interface between cell (7,6) and (7,7) (j-direction interface)
    Location = SolnBlk[0].Grid.xfaceN(7,6);
    Result = SolnBlk[0].ExactSoln->Solution(Location.x,Location.y);

    // Numeric result
    SolnBlk[0].EllipticFluxStateAtInteriorInterface(7,6,7,7,Location,U_face);

    ensure_distance("Gradient betwen cells (7,6) and (7,7)", U_face.u, Result, AcceptedError(Result,1.0e-8));
  }

  /* Test 6:*/
  template<>
  template<>
  void AdvectDiffuse2D_Quad_Block_object::test<6>()
  {

    set_test_name("Generate advection-diffusion solution block with high-order variables");
    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    // Set input file name
    Open_Input_File("CircularAdvectionDiffusion_HighOrder.in");

    // Parse the input file
    IP.Verbose() = false;
    IP.Parse_Input_File(input_file_name);
    HighOrder2D_Input::Set_Final_Parameters(IP);

    // Create computational domain
    InitializeComputationalDomain(MeshBlk,QuadTree,
				  GlobalList_Soln_Blocks, LocalList_Soln_Blocks, 
				  SolnBlk, IP);

    // == check correct initialization
    ensure("High-order variables", SolnBlk[0].HighOrderVariables() != NULL);
    ensure_equals("Main High-order ", SolnBlk[0].HighOrderVariable(0).RecOrder(), 3);
    ensure_equals("Second High-order ", SolnBlk[0].HighOrderVariable(1).RecOrder(), 2);
    ensure_equals("Third High-order ", SolnBlk[0].HighOrderVariable(2).RecOrder(), 4);

    // Change the reconstruction orders of the objects
    vector<short int> RecOrders(3);
    RecOrders[0] = 0; RecOrders[1] = 4; RecOrders[2] = 3; 

    SolnBlk[0].allocate_HighOrder(3,RecOrders);

    // == check correct re-initialization of the high-order objects
    ensure_equals("Main High-order ", SolnBlk[0].HighOrderVariable(0).RecOrder(), RecOrders[0]);
    ensure_equals("Second High-order ", SolnBlk[0].HighOrderVariable(1).RecOrder(), RecOrders[1]);
    ensure_equals("Third High-order ", SolnBlk[0].HighOrderVariable(2).RecOrder(), RecOrders[2]);

    // == check calculation of the pseudo-inverse
    

    // associate a different grid to the third high-order object
    SolnBlk[0].HighOrderVariable(2).AssociateGeometry(SolnBlk[1].Grid);

    // == check calculation of the pseudo-inverse with the new grid
    
  }

  /* Test 7:*/
  template<>
  template<>
  void AdvectDiffuse2D_Quad_Block_object::test<7>()
  {

    set_test_name("Copy advection-diffusion solution block with high-order variables");
    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    // Set input file name
    Open_Input_File("CircularAdvectionDiffusion_HighOrder.in");

    // Parse the input file
    IP.Verbose() = false;
    IP.Parse_Input_File(input_file_name);
    HighOrder2D_Input::Set_Final_Parameters(IP);

    // Create computational domain
    InitializeComputationalDomain(MeshBlk,QuadTree,
				  GlobalList_Soln_Blocks, LocalList_Soln_Blocks, 
				  SolnBlk, IP);

    // == check correct initialization
    ensure("High-order variables", SolnBlk[0].HighOrderVariables() != NULL);
    ensure_equals("Main High-order ", SolnBlk[0].HighOrderVariable(0).RecOrder(), 3);
    ensure_equals("Second High-order ", SolnBlk[0].HighOrderVariable(1).RecOrder(), 2);
    ensure_equals("Third High-order ", SolnBlk[0].HighOrderVariable(2).RecOrder(), 4);

    // copy quad block
    Copy_Solution_Block(Soln,SolnBlk[0]);

    // == check copy
    ensure("High-order variables", Soln.HighOrderVariables() != NULL);
    ensure_equals("Main High-order ", Soln.HighOrderVariable(0).RecOrder(), 3);
    ensure_equals("Second High-order ", Soln.HighOrderVariable(1).RecOrder(), 2);
    ensure_equals("Third High-order ", Soln.HighOrderVariable(2).RecOrder(), 4);
  }

  

}



// Test suite constructor
tut::AdvectDiffuse2D_Quad_Block_TestSuite AdvectDiffuse2D_Quad_BlockTestSuite("Class:AdvectDiffuse2D_Quad_Block");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

