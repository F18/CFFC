/*!\file test_NavierStokes2DQuadBlock.cc
  \brief Regression tests for class NavierStokes2D_Quad_Block. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../NavierStokes2DQuad.h"
#include "../Grid/HO_Grid2DQuadMultiBlock.h" /* Include 2D quadrilateral multiblock grid header file */
#include "../../HighOrderReconstruction/AccuracyAssessment2DMultiBlock.h"
#include "../../HighOrderReconstruction/HighOrder2D_MultiBlock.h" /* Include 2D high-order header file for multi-block level. */

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_NavierStokes2D_Quad_Block : public TestData {

    // Local variables
  public:

    // ==== Member data =====

    // NavierStokes2D input variables and parameters:
    NavierStokes2D_Input_Parameters IP;
    
    /* Multi-block solution-adaptive quadrilateral mesh 
       solution variables. */
    Grid2D_Quad_MultiBlock_HO     MeshBlk;
    QuadTreeBlock_DataStructure  QuadTree;
    AdaptiveBlockResourceList    GlobalList_Soln_Blocks;
    AdaptiveBlock2D_List         LocalList_Soln_Blocks;
    NavierStokes2D_Quad_Block   *SolnBlk;
    NavierStokes2D_Quad_Block    Soln;
    CPUTime processor_cpu_time;

    int error_flag;
    int Status;                        // shows if the computational domain has been initialized

    //===== Member functions =====

    // Default Constructor
    Data_NavierStokes2D_Quad_Block(void);

    // Destructor
    ~Data_NavierStokes2D_Quad_Block(void);

    // Initialization of the computational domain
    void InitializeComputationalDomain(Grid2D_Quad_MultiBlock_HO & _MeshBlk_,
				       QuadTreeBlock_DataStructure & _QuadTree_,
				       AdaptiveBlockResourceList & _GlobalList_Soln_Blocks_,
				       AdaptiveBlock2D_List & _LocalList_Soln_Blocks_,
				       NavierStokes2D_Quad_Block *& _SolnBlk_,
				       NavierStokes2D_Input_Parameters & _IP_) throw(std::runtime_error);

    // Output_Block()
    void Output_Block(NavierStokes2D_Quad_Block & SolnBlock,
		      NavierStokes2D_Input_Parameters & _IP_){};

    // Set the local time step to value for all solution blocks.
    void SetLocalTimeStepToValue(NavierStokes2D_Quad_Block *& _SolnBlk_,
				 AdaptiveBlock2D_List & _LocalList_Soln_Blocks_,
				 const double & ValueToBeSet);    

  private:
    
  };


  Data_NavierStokes2D_Quad_Block::Data_NavierStokes2D_Quad_Block(){
    /* Set the global path for this test suite. 
       It's a good practice to put it in the constructor of the data class in order to be set
       automatically for each individual test. Declare it relative to the /src_2D directory,
       otherwise the framework might not find the input and output files. */
    
    set_test_suite_path("NavierStokes2D/UnitTests");

    // Initialize IP to default values
    Set_Default_Input_Parameters(IP);

    Status = OFF;
  }

  Data_NavierStokes2D_Quad_Block::~Data_NavierStokes2D_Quad_Block(void){
    if (Status == ON){
      SolnBlk = Deallocate(SolnBlk, IP);
      LocalList_Soln_Blocks.deallocate();
      GlobalList_Soln_Blocks.deallocate();
      QuadTree.deallocate();
    }
    HO_Grid2D_Execution_Mode::SetDefaults();
  }

  // === InitializeComputationalDomain()
  void Data_NavierStokes2D_Quad_Block::InitializeComputationalDomain(Grid2D_Quad_MultiBlock_HO & _MeshBlk_,
								     QuadTreeBlock_DataStructure & _QuadTree_,
								     AdaptiveBlockResourceList & _GlobalList_Soln_Blocks_,
								     AdaptiveBlock2D_List & _LocalList_Soln_Blocks_,
								     NavierStokes2D_Quad_Block *& _SolnBlk_,
								     NavierStokes2D_Input_Parameters & _IP_)
    throw(std::runtime_error){

    error_flag = _MeshBlk_.Multi_Block_Grid(_IP_);
      
    if (error_flag) {
      throw runtime_error("CreateMesh() ERROR: Unable to create valid NavierStokes2D multi-block mesh.");
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
      throw runtime_error("Create_Initial_Solution_Blocks() ERROR: Unable to create initial NavierStokes2D solution blocks.");
    }

    for (int n = 0 ; n <= _LocalList_Soln_Blocks_.Nblk-1 ; ++n ) {
      if (_LocalList_Soln_Blocks_.Block[n].used == ADAPTIVEBLOCK2D_USED) {
	if (_SolnBlk_[n].Grid.Check_Quad_Block_Completely() ){
	  throw runtime_error("Create_Initial_Solution_Blocks() ERROR: Invalid quadrilateral block detected!");
	}
      } /* endif */
    }  /* endfor */

    Status = ON;
  }


  // === SetLocalTimeStepToValue()
  void Data_NavierStokes2D_Quad_Block::SetLocalTimeStepToValue(NavierStokes2D_Quad_Block *& _SolnBlk_,
							AdaptiveBlock2D_List & _LocalList_Soln_Blocks_,
							const double & ValueToBeSet){

    int n, i, j;

    for ( n = 0 ; n <= _LocalList_Soln_Blocks_.Nblk-1 ; ++n ) {
      if (_LocalList_Soln_Blocks_.Block[n].used == ADAPTIVEBLOCK2D_USED) {
	for (j = _SolnBlk_[n].JCl; j <= _SolnBlk_[n].JCu; ++j ){
	  for (i = _SolnBlk_[n].ICl; i <= _SolnBlk_[n].ICu; ++i ){
	    _SolnBlk_[n].dt[i][j] = ValueToBeSet;
	  } /* endfor */
	} /* endfor */	
      } /* endif */

    }  /* endfor */
  }

  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. NavierStokes2D_Quad_Block_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_NavierStokes2D_Quad_Block> NavierStokes2D_Quad_Block_TestSuite;
  typedef NavierStokes2D_Quad_Block_TestSuite::object NavierStokes2D_Quad_Block_object;


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
  void NavierStokes2D_Quad_Block_object::test<1>()
  {

    set_test_name("Flat plate");

    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    RunRegression = OFF;
 
    // Set input file name
    Open_Input_File("flatplate.in");

    // Parse the input file
    IP.Verbose() = false;
    IP.Parse_Input_File(input_file_name);
    HighOrder2D_Input::Set_Final_Parameters(IP);
    CENO_Execution_Mode::CENO_SPEED_EFFICIENT = OFF;

    // Create computational domain
    InitializeComputationalDomain(MeshBlk,QuadTree,
				  GlobalList_Soln_Blocks, LocalList_Soln_Blocks, 
				  SolnBlk, IP);

    // Apply initial condition
    ICs(SolnBlk,LocalList_Soln_Blocks,IP);


    // Send solution information between neighbouring blocks to complete
    // prescription of initial data.
    error_flag = Send_All_Messages(SolnBlk,
				   LocalList_Soln_Blocks,
				   NUM_COMP_VECTOR2D,
				   ON);
    if (!error_flag) error_flag = Send_All_Messages(SolnBlk,
						    LocalList_Soln_Blocks,
						    NUM_VAR_NAVIERSTOKES2D,
						    OFF);
    
    // Prescribe boundary data consistent with initial data.
    BCs(SolnBlk,LocalList_Soln_Blocks,IP);
    
    error_flag = Boundary_Adaptive_Mesh_Refinement(SolnBlk,
						   IP,
						   QuadTree,
						   GlobalList_Soln_Blocks,
						   LocalList_Soln_Blocks);
    
  }

}



// Test suite constructor
tut::NavierStokes2D_Quad_Block_TestSuite NavierStokes2D_Quad_BlockTestSuite("Class:NavierStokes2D_Quad_Block");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

