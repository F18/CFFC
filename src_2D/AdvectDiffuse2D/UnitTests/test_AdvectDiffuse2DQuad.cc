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

    // Set the local time step to value for all solution blocks.
    void SetLocalTimeStepToValue(AdvectDiffuse2D_Quad_Block *& _SolnBlk_,
				 AdaptiveBlock2D_List & _LocalList_Soln_Blocks_,
				 const double & ValueToBeSet);

    // Compute right-hand-side term of the equation
    void ComputeEquationRightHandSideTerm(AdvectDiffuse2D_Quad_Block & _SolnBlk_,
					  AdvectDiffuse2D_Input_Parameters & _IP_,
					  const int & k_residual);

    // Compute residual errors
    void ComputeResidualErrors(AdvectDiffuse2D_Quad_Block & _SolnBlk_,
			       const int & k_residual_I, const int & k_residual_II,
			       double & L1, double & L2, double & LMax);
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
  void Data_AdvectDiffuse2D_Quad_Block::SetLocalTimeStepToValue(AdvectDiffuse2D_Quad_Block *& _SolnBlk_,
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

  // === Compute the integral of the right-hand-side term of the equation based on what the exact solution has set
  void Data_AdvectDiffuse2D_Quad_Block::ComputeEquationRightHandSideTerm(AdvectDiffuse2D_Quad_Block & _SolnBlk_,
									 AdvectDiffuse2D_Input_Parameters & _IP_,
									 const int & k_residual){

    int i,j;
    double IntResult;

    // Calculate the integral of the right hand side term over the domain of each interior cell divided by the local area
    // Use the ExactSoln pointer to access the exact solution
    if (_SolnBlk_.ExactSoln->IsExactSolutionSet()) {
      for ( j  = _SolnBlk_.JCl ; j <= _SolnBlk_.JCu ; ++j ) {
 	for ( i = _SolnBlk_.ICl ; i <= _SolnBlk_.ICu ; ++i ) {
	  IntResult = 
	    _SolnBlk_.Grid.Integration.IntegrateFunctionOverCell(i,j,
								 wrapped_member_function(_SolnBlk_.ExactSoln,
											 &AdvectDiffuse2D_ExactSolutions::
											 PDE_RightHandSide,
											 IntResult),
								 wrapped_member_function(_SolnBlk_.ExactSoln,
											 &AdvectDiffuse2D_ExactSolutions::
											 XDependencyIntegrated_PDE_RightHandSide,
											 IntResult),
								 IP.Exact_Integration_Digits,
								 IntResult)/_SolnBlk_.Grid.Cell[i][j].A;
	  _SolnBlk_.dUdt[i][j][k_residual] = IntResult;
	} /* endfor */
      } /* endfor */
    } else {
      // There is no exact solution set for this problem
      throw runtime_error("ICs() ERROR! No exact solution has been set!");
    }
  }


  // === Compute the errors between the left and right hand sides of the equation.
  //     The k_residual_I and k_residual_II indicate where the two residuals are stored.
  //     The residual of the first index is going to be overwritten for error plotting.
  //     Return error statistics in L1, L2 and LMax.
  void Data_AdvectDiffuse2D_Quad_Block::ComputeResidualErrors(AdvectDiffuse2D_Quad_Block & _SolnBlk_,
							      const int & k_residual_I, const int & k_residual_II,
							      double & L1, double & L2, double & LMax){

    int i,j, counter(0);

    // Initialize error norms
    L1 = L2 = LMax = 0.0;

    for ( j  = _SolnBlk_.JCl ; j <= _SolnBlk_.JCu ; ++j ) {
      for ( i = _SolnBlk_.ICl ; i <= _SolnBlk_.ICu ; ++i, ++counter ) {

	// Compute error
	_SolnBlk_.dUdt[i][j][k_residual_I] = fabs(_SolnBlk_.dUdt[i][j][k_residual_I] - 
						  _SolnBlk_.dUdt[i][j][k_residual_II]);

	// Compute error norms
	L1 += _SolnBlk_.dUdt[i][j][k_residual_I].u;
	L2 += sqr(_SolnBlk_.dUdt[i][j][k_residual_I].u);
	LMax = max(LMax,_SolnBlk_.dUdt[i][j][k_residual_I].u);
	
      } /* endfor */
    } /* endfor */

    // Compute final errors
    L1 /= counter;
    L2 /= counter; L2 = sqrt(L2);

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
    CENO_Execution_Mode::CENO_SPEED_EFFICIENT = ON;

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
    vector<int> RecOrders(3);
    RecOrders[0] = 0; RecOrders[1] = 4; RecOrders[2] = 3; 

    SolnBlk[0].allocate_HighOrder(3,RecOrders);

    // == check correct re-initialization of the high-order objects
    ensure_equals("Main High-order ", SolnBlk[0].HighOrderVariable(0).RecOrder(), RecOrders[0]);
    ensure_equals("Second High-order ", SolnBlk[0].HighOrderVariable(1).RecOrder(), RecOrders[1]);
    ensure_equals("Third High-order ", SolnBlk[0].HighOrderVariable(2).RecOrder(), RecOrders[2]);

    // == check calculation of the pseudo-inverse
    ensure_equals("Main High-order pseudo-inv", SolnBlk[0].HighOrderVariable(0).IsPseudoInversePreComputed(), true);
    ensure_equals("Second High-order pseudo-inv", SolnBlk[0].HighOrderVariable(1).IsPseudoInversePreComputed(), true);
    ensure_equals("Third High-order pseudo-inv", SolnBlk[0].HighOrderVariable(2).IsPseudoInversePreComputed(), true);

    // associate a different grid to the third high-order object
    SolnBlk[0].HighOrderVariable(2).AssociateGeometry(SolnBlk[1].Grid);

    // == check calculation of the pseudo-inverse with the new grid
    ensure_equals("Main High-order pseudo-inv", SolnBlk[0].HighOrderVariable(0).IsPseudoInversePreComputed(), true);
    ensure_equals("Second High-order pseudo-inv", SolnBlk[0].HighOrderVariable(1).IsPseudoInversePreComputed(), true);
    ensure_equals("Third High-order pseudo-inv", SolnBlk[0].HighOrderVariable(2).IsPseudoInversePreComputed(), true);    
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

  /* Test 8:*/
  template<>
  template<>
  void AdvectDiffuse2D_Quad_Block_object::test<8>()
  {

    set_test_name("Output advection-diffusion solution block with high-order variables");
    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    RunRegression = ON ;

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

    MasterFile = "Output_AdvectionDiffusionBlock.dat";
    CurrentFile = "Current_Output_AdvectionDiffusionBlock.dat";

    if (RunRegression){
      
      // open file for output
      Open_Output_File(CurrentFile);

      // generate current file
      out() << SolnBlk[0];

      // check result
      RunRegressionTest("operator <<", CurrentFile, MasterFile, 5.0e-12, 5.0e-12);

    } else {

      // open file for output
      Open_Output_File(MasterFile);
      
      // generate master
      out() << SolnBlk[0];
    }
    
  }

  /* Test 9:*/
  template<>
  template<>
  void AdvectDiffuse2D_Quad_Block_object::test<9>()
  {

    set_test_name("Check input-output operators");
    set_local_input_path("QuadBlockData");

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

    // == check 
    Check_Input_Output_Operator("SolnBlk[0] variable", SolnBlk[0]);
  }
  
  /* Test 10:*/
  template<>
  template<>
  void AdvectDiffuse2D_Quad_Block_object::test<10>()
  {

    set_test_name("Compute high-order reconstruction with pseudo-inverse");
    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    // Set input file name
    Open_Input_File("CircularAdvectionDiffusion_HighOrder_ErrorStudy.in");

    // Parse the input file
    IP.Verbose() = false;
    IP.Parse_Input_File(input_file_name);
    HighOrder2D_Input::Set_Final_Parameters(IP);
    CENO_Execution_Mode::CENO_SPEED_EFFICIENT = OFF;

    // Create computational domain
    InitializeComputationalDomain(MeshBlk,QuadTree,
				  GlobalList_Soln_Blocks, LocalList_Soln_Blocks, 
				  SolnBlk, IP);

    // == check correct initialization
    ensure("High-order variables", SolnBlk[0].HighOrderVariables() != NULL);
    ensure_equals("Main High-order ", SolnBlk[0].HighOrderVariable(0).RecOrder(), 3);
    ensure_equals("Second High-order ", SolnBlk[0].HighOrderVariable(1).RecOrder(), 2);
    ensure_equals("Third High-order ", SolnBlk[0].HighOrderVariable(2).RecOrder(), 4);
    
    // Apply initial condition
    ICs(SolnBlk,LocalList_Soln_Blocks,IP);

    // Reconstruct solution
    SolnBlk[0].HighOrderVariable(0).ComputeUnlimitedSolutionReconstruction(SolnBlk[0]);
    SolnBlk[0].HighOrderVariable(1).ComputeUnlimitedSolutionReconstruction(SolnBlk[0]);
    SolnBlk[0].HighOrderVariable(2).ComputeUnlimitedSolutionReconstruction(SolnBlk[0]);

    // Compute solution error
    double Error;

    SolnBlk[0].HighOrderVariable(0).ComputeSolutionErrors(wrapped_member_function(SolnBlk[0].ExactSolution(),
										  &AdvectDiffuse2D_Quad_Block::
										  Exact_Solution_Type::Solution,
										  Error),
							  1, 12);
    SolnBlk[0].HighOrderVariable(1).ComputeSolutionErrors(wrapped_member_function(SolnBlk[0].ExactSolution(),
										  &AdvectDiffuse2D_Quad_Block::
										  Exact_Solution_Type::Solution,
										  Error),
							  1, 12);
    SolnBlk[0].HighOrderVariable(2).ComputeSolutionErrors(wrapped_member_function(SolnBlk[0].ExactSolution(),
										  &AdvectDiffuse2D_Quad_Block::
										  Exact_Solution_Type::Solution,
										  Error),
							  1, 12);

    // Check errors against values determined in a grid convergence study which reproduced the expected order of accuracy
    ensure_distance("HO_0, L1", SolnBlk[0].HighOrderVariable(0).L1(),0.001801828715892732 , AcceptedError(0.001801828715892732));
    ensure_distance("HO_0, L2", SolnBlk[0].HighOrderVariable(0).L2(),2.22858063792015e-05 , AcceptedError(2.22858063792015e-05));
    ensure_distance("HO_0, LMax", SolnBlk[0].HighOrderVariable(0).LMax(),0.01647612251871387, AcceptedError(0.01647612251871387));
    ensure_distance("HO_0, BlockArea", SolnBlk[0].HighOrderVariable(0).BlockArea(), 1.0, AcceptedError(1.0));
    ensure_distance("HO_0, Block L1", SolnBlk[0].HighOrderVariable(0).BlockL1Norm(), 
		    0.001801828715892732, AcceptedError(0.001801828715892732));
    ensure_distance("HO_0, Block L2", SolnBlk[0].HighOrderVariable(0).BlockL2Norm(), 
		    0.00472078450887154, AcceptedError(0.00472078450887154));
    ensure_distance("HO_0, Block LMax", SolnBlk[0].HighOrderVariable(0).BlockLMaxNorm(), 
		    0.01647612251871387, AcceptedError(0.01647612251871387));
    ensure_equals("HO_1, Cells used", SolnBlk[0].HighOrderVariable(0).UsedCells(), 64);


    ensure_distance("HO_1, L1", SolnBlk[0].HighOrderVariable(1).L1(),0.008666756828796066 , AcceptedError(0.008666756828796066));
    ensure_distance("HO_1, L2", SolnBlk[0].HighOrderVariable(1).L2(),0.000483179638195169 , AcceptedError(0.000483179638195169));
    ensure_distance("HO_1, LMax", SolnBlk[0].HighOrderVariable(1).LMax(),0.07216104037829933, AcceptedError(0.07216104037829933));
    ensure_distance("HO_1, BlockArea", SolnBlk[0].HighOrderVariable(1).BlockArea(), 1.0, AcceptedError(1.0));
    ensure_distance("HO_1, Block L1", SolnBlk[0].HighOrderVariable(1).BlockL1Norm(), 
		    0.008666756828796066, AcceptedError(0.008666756828796066));
    ensure_distance("HO_1, Block L2", SolnBlk[0].HighOrderVariable(1).BlockL2Norm(), 
		    0.02198134750635568, AcceptedError(0.02198134750635568));
    ensure_distance("HO_1, Block LMax", SolnBlk[0].HighOrderVariable(1).BlockLMaxNorm(), 
		    0.07216104037829933, AcceptedError(0.07216104037829933));
    ensure_equals("HO_1, Cells used", SolnBlk[0].HighOrderVariable(1).UsedCells(), 64);


    ensure_distance("HO_2, L1", SolnBlk[0].HighOrderVariable(2).L1(),0.001390104987434202 , AcceptedError(0.001390104987434202));
    ensure_distance("HO_2, L2", SolnBlk[0].HighOrderVariable(2).L2(),1.383174295956005e-05, AcceptedError(1.383174295956005e-05));
    ensure_distance("HO_2, LMax", SolnBlk[0].HighOrderVariable(2).LMax(),0.01329552268181446, AcceptedError(0.01329552268181446));
    ensure_distance("HO_2, BlockArea", SolnBlk[0].HighOrderVariable(2).BlockArea(), 1.0, AcceptedError(1.0));
    ensure_distance("HO_2, Block L1", SolnBlk[0].HighOrderVariable(2).BlockL1Norm(), 
		    0.001390104987434202, AcceptedError(0.001390104987434202));
    ensure_distance("HO_2, Block L2", SolnBlk[0].HighOrderVariable(2).BlockL2Norm(), 
		    0.003719105128866358, AcceptedError(0.003719105128866358));
    ensure_distance("HO_2, Block LMax", SolnBlk[0].HighOrderVariable(2).BlockLMaxNorm(), 
		    0.01329552268181446, AcceptedError(0.01329552268181446));
    ensure_equals("HO_2, Cells used", SolnBlk[0].HighOrderVariable(2).UsedCells(), 64);

  }

  /* Test 11:*/
  template<>
  template<>
  void AdvectDiffuse2D_Quad_Block_object::test<11>()
  {

    set_test_name("Check high-order output functions");
    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    RunRegression = ON;

    // Set input file name
    Open_Input_File("CircularAdvectionDiffusion_HighOrder_ErrorStudy.in");

    // Parse the input file
    IP.Verbose() = false;
    IP.Parse_Input_File(input_file_name);
    HighOrder2D_Input::Set_Final_Parameters(IP);
    CENO_Execution_Mode::CENO_SPEED_EFFICIENT = OFF;

    // Create computational domain
    InitializeComputationalDomain(MeshBlk,QuadTree,
				  GlobalList_Soln_Blocks, LocalList_Soln_Blocks, 
				  SolnBlk, IP);

    // == check correct initialization
    ensure("High-order variables", SolnBlk[0].HighOrderVariables() != NULL);
    ensure_equals("Main High-order ", SolnBlk[0].HighOrderVariable(0).RecOrder(), 3);
    ensure_equals("Second High-order ", SolnBlk[0].HighOrderVariable(1).RecOrder(), 2);
    ensure_equals("Third High-order ", SolnBlk[0].HighOrderVariable(2).RecOrder(), 4);
    
    // Apply initial condition
    ICs(SolnBlk,LocalList_Soln_Blocks,IP);

    // Reconstruct the solution
    SolnBlk[0].HighOrderVariable(2).ComputeUnlimitedSolutionReconstruction(SolnBlk[0]);

    if (RunRegression){
      //===== Check solution interior nodal solution
      MasterFile  = "HighOrder_Circular_Advection_In_Rectangular_Box_cpu000000_block000.dat";
      CurrentFile = "Current_HighOrder_Circular_Advection_In_Rectangular_Box_cpu000000_block000.dat";

      // Output 
      SolnBlk[0].Output_Tecplot_HighOrder_Debug_Mode(LocalList_Soln_Blocks, IP, 0, 2);
      
      // check
      RunRegressionTest("HighOrder Interior Nodal Solution Tecplot", CurrentFile, MasterFile, 5.0e-7, 5.0e-7);

      //===== Check solution interior nodal solution
      MasterFile  = "HighOrder_Circular_Advection_In_Rectangular_Box_nodes_cpu000000_block000.dat";
      CurrentFile = "Current_HighOrder_Circular_Advection_In_Rectangular_Box_nodes_cpu000000_block000.dat";
      
      // Output 
      SolnBlk[0].Output_Nodes_Tecplot_HighOrder_Debug_Mode(LocalList_Soln_Blocks, IP, 0, 2);

      // check
      RunRegressionTest("HighOrder Interior+Ghost Nodal Solution Tecplot", CurrentFile, MasterFile, 5.0e-7, 5.0e-7);

      //===== Check solution centroid solution
      MasterFile  = "HighOrder_Circular_Advection_In_Rectangular_Box_cells_cpu000000_block000.dat";
      CurrentFile = "Current_HighOrder_Circular_Advection_In_Rectangular_Box_cells_cpu000000_block000.dat";

      // Output 
      SolnBlk[0].Output_Cells_Tecplot_HighOrder_Debug_Mode(LocalList_Soln_Blocks, IP, 0, 2);
      
      // check
      RunRegressionTest("HighOrder Centroid Solution Tecplot", CurrentFile, MasterFile, 5.0e-7, 5.0e-7);

    } else {
      // generate the master files

      // Output 
      SolnBlk[0].Output_Tecplot_HighOrder_Debug_Mode(LocalList_Soln_Blocks, IP, 0, 2);
      SolnBlk[0].Output_Nodes_Tecplot_HighOrder_Debug_Mode(LocalList_Soln_Blocks, IP, 0, 2);
      SolnBlk[0].Output_Cells_Tecplot_HighOrder_Debug_Mode(LocalList_Soln_Blocks, IP, 0, 2);
    }
  }

  /* Test 12:*/
  template<>
  template<>
  void AdvectDiffuse2D_Quad_Block_object::test<12>()
  {

    set_test_name("Piecewise Linear Reconstruction");
    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    RunRegression = ON;

    // Indexes of the tested cell
    int iCell,jCell; 
    iCell = 7;
    jCell = 7;

    // Set input file name
    Open_Input_File("CircularAdvectionDiffusion_HighOrder_With_DropOrder.in");

    // Parse the input file
    IP.Verbose() = false;
    IP.Parse_Input_File(input_file_name);
    HighOrder2D_Input::Set_Final_Parameters(IP);
    CENO_Execution_Mode::CENO_SPEED_EFFICIENT = OFF;

    // Create computational domain
    InitializeComputationalDomain(MeshBlk,QuadTree,
				  GlobalList_Soln_Blocks, LocalList_Soln_Blocks, 
				  SolnBlk, IP);

    // == check correct initialization
    ensure("High-order variables", SolnBlk[0].HighOrderVariables() != NULL);
    ensure_equals("Main High-order ", SolnBlk[0].HighOrderVariable(0).RecOrder(), 4);
    ensure_equals("Main High-order ", SolnBlk[0].HighOrderVariable(1).RecOrder(), 1);
    
    // Apply initial condition
    ICs(SolnBlk,LocalList_Soln_Blocks,IP);

    // Reconstruct the solution without geometric weighting
    CENO_Execution_Mode::CENO_APPLY_GEOMETRIC_WEIGHTING = OFF;
    SolnBlk[0].HighOrderVariable(0).ComputeUnlimitedSolutionReconstruction(SolnBlk[0]);
    SolnBlk[0].HighOrderVariable(1).ComputeUnlimitedSolutionReconstruction(SolnBlk[0]);
    // Reconstruct with PWL 
    SolnBlk[0].HighOrderVariable(0).CellInadequateFitValue(iCell,jCell,1) = ON;
    SolnBlk[0].HighOrderVariable(0).ComputeLimitedPiecewiseLinearSolutionReconstruction(SolnBlk[0],
											iCell,jCell,
											IP.Limiter());
    SolnBlk[0].HighOrderVariable(1).CellInadequateFitValue(iCell,jCell,1) = ON;
    SolnBlk[0].HighOrderVariable(1).ComputeLimitedPiecewiseLinearSolutionReconstruction(SolnBlk[0],
											iCell,jCell,
											IP.Limiter());

    Linear_Reconstruction_LeastSquares(SolnBlk[0],IP.Limiter());

    // == check HighOrderVariable(0)
    ensure_distance("D00, I",
		    SolnBlk[0].HighOrderVariable(0).CellTaylorDerivState(iCell,jCell,0),
		    SolnBlk[0].U[iCell][jCell],
		    AcceptedError(SolnBlk[0].U[iCell][jCell]));
    ensure_distance("D01, Drop Order, I",
		    SolnBlk[0].HighOrderVariable(0).get_dUdy(),
		    SolnBlk[0].dUdy[iCell][jCell],
		    AcceptedError(SolnBlk[0].dUdy[iCell][jCell]));
    ensure_distance("D10, Drop Order, I",
		    SolnBlk[0].HighOrderVariable(0).get_dUdx(),
		    SolnBlk[0].dUdx[iCell][jCell],
		    AcceptedError(SolnBlk[0].dUdx[iCell][jCell]));
    ensure_distance("D01, I",
		    SolnBlk[0].HighOrderVariable(0).CellTaylorDerivState(iCell,jCell,0,1),
		    SolnBlk[0].dUdy[iCell][jCell],
		    AcceptedError(SolnBlk[0].dUdy[iCell][jCell]));
    ensure_distance("D10, I",
		    SolnBlk[0].HighOrderVariable(0).CellTaylorDerivState(iCell,jCell,1,0),
		    SolnBlk[0].dUdx[iCell][jCell],
		    AcceptedError(SolnBlk[0].dUdx[iCell][jCell]));
    ensure_distance("phi, Drop Order, I",
		    SolnBlk[0].HighOrderVariable(0).CellTaylorDeriv(iCell,jCell).Limiter(1),
		    9.1528284874531129e-02,
		    AcceptedError(9.1528284874531129e-02));
    ensure_distance("D02, I", SolnBlk[0].HighOrderVariable(0).CellTaylorDerivState(iCell,jCell,0,2)[1], 0.0, AcceptedError(0.0));
    ensure_distance("D03, I", SolnBlk[0].HighOrderVariable(0).CellTaylorDerivState(iCell,jCell,0,3)[1], 0.0, AcceptedError(0.0));
    ensure_distance("D04, I", SolnBlk[0].HighOrderVariable(0).CellTaylorDerivState(iCell,jCell,0,4)[1], 0.0, AcceptedError(0.0));
    ensure_distance("D11, I", SolnBlk[0].HighOrderVariable(0).CellTaylorDerivState(iCell,jCell,1,1)[1], 0.0, AcceptedError(0.0));
    ensure_distance("D12, I", SolnBlk[0].HighOrderVariable(0).CellTaylorDerivState(iCell,jCell,1,2)[1], 0.0, AcceptedError(0.0));
    ensure_distance("D13, I", SolnBlk[0].HighOrderVariable(0).CellTaylorDerivState(iCell,jCell,1,3)[1], 0.0, AcceptedError(0.0));
    ensure_distance("D20, I", SolnBlk[0].HighOrderVariable(0).CellTaylorDerivState(iCell,jCell,2,0)[1], 0.0, AcceptedError(0.0));
    ensure_distance("D21, I", SolnBlk[0].HighOrderVariable(0).CellTaylorDerivState(iCell,jCell,2,1)[1], 0.0, AcceptedError(0.0));
    ensure_distance("D22, I", SolnBlk[0].HighOrderVariable(0).CellTaylorDerivState(iCell,jCell,2,2)[1], 0.0, AcceptedError(0.0));
    ensure_distance("D31, I", SolnBlk[0].HighOrderVariable(0).CellTaylorDerivState(iCell,jCell,3,1)[1], 0.0, AcceptedError(0.0));
    ensure_distance("D40, I", SolnBlk[0].HighOrderVariable(0).CellTaylorDerivState(iCell,jCell,4,0)[1], 0.0, AcceptedError(0.0));

    // == check HighOrderVariable(1)
    ensure_distance("D00, II",
		    SolnBlk[0].HighOrderVariable(1).CellTaylorDerivState(iCell,jCell,0),
		    SolnBlk[0].U[iCell][jCell],
		    AcceptedError(SolnBlk[0].U[iCell][jCell]));
    ensure_distance("D01, II",
		    SolnBlk[0].HighOrderVariable(1).CellTaylorDerivState(iCell,jCell,1),
		    SolnBlk[0].dUdy[iCell][jCell],
		    AcceptedError(SolnBlk[0].dUdy[iCell][jCell]));
    ensure_distance("D10, II",
		    SolnBlk[0].HighOrderVariable(1).CellTaylorDerivState(iCell,jCell,2),
		    SolnBlk[0].dUdx[iCell][jCell],
		    AcceptedError(SolnBlk[0].dUdx[iCell][jCell]));
    ensure_distance("phi, Drop Order, II",
		    SolnBlk[0].HighOrderVariable(1).CellTaylorDeriv(iCell,jCell).Limiter(1),
		    9.1528284874531129e-02,
		    AcceptedError(9.1528284874531129e-02));
  }

  /* Test 13:*/
  template<>
  template<>
  void AdvectDiffuse2D_Quad_Block_object::test<13>()
  {

    set_test_name("Check smoothness indicator");
    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    RunRegression = ON;

    // Set input file name
    Open_Input_File("CircularAdvectionDiffusion_HighOrder_SmoothnessIndicatorStudy.in");

    // Parse the input file
    IP.Verbose() = false;
    IP.Parse_Input_File(input_file_name);

    // Create computational domain
    InitializeComputationalDomain(MeshBlk,QuadTree,
				  GlobalList_Soln_Blocks, LocalList_Soln_Blocks, 
				  SolnBlk, IP);

    // == check correct initialization
    ensure("High-order variables", SolnBlk[0].HighOrderVariables() != NULL);
    ensure_equals("Main High-order ", SolnBlk[0].HighOrderVariable(0).RecOrder(), 3);
    ensure_equals("Second High-order ", SolnBlk[0].HighOrderVariable(1).RecOrder(), 2);
    ensure_equals("Third High-order ", SolnBlk[0].HighOrderVariable(2).RecOrder(), 4);
    
    // Apply initial condition
    ICs(SolnBlk,LocalList_Soln_Blocks,IP);

    // Reconstruct the solution
    SolnBlk[0].HighOrderVariable(0).ComputeUnlimitedSolutionReconstruction(SolnBlk[0]);
    SolnBlk[0].HighOrderVariable(1).ComputeUnlimitedSolutionReconstruction(SolnBlk[0]);
    SolnBlk[0].HighOrderVariable(2).ComputeUnlimitedSolutionReconstruction(SolnBlk[0]);

    // Calculate the smoothness indicator for the reconstructions
    SolnBlk[0].HighOrderVariable(0).ComputeSmoothnessIndicator(SolnBlk[0]);
    SolnBlk[0].HighOrderVariable(1).ComputeSmoothnessIndicator(SolnBlk[0]);
    SolnBlk[0].HighOrderVariable(2).ComputeSmoothnessIndicator(SolnBlk[0]);

    if (RunRegression) {

      //===== Check interior nodal values
      // Output Order 3
      MasterFile = "HighOrder_SmoothnessIndicator_Order3.dat";
      CurrentFile = "Current_HighOrder_SmoothnessIndicator_Order3.dat";
      Open_Output_File(CurrentFile);

      SolnBlk[0].Output_Nodes_Tecplot_HighOrder(0,0,0, 1, out(), 0);

      // check
      RunRegressionTest("Order 3", CurrentFile, MasterFile, 5.0e-10, 5.0e-10);

      // Output Order 2
      MasterFile = "HighOrder_SmoothnessIndicator_Order2.dat";
      CurrentFile = "Current_HighOrder_SmoothnessIndicator_Order2.dat";
      Open_Output_File(CurrentFile);

      SolnBlk[0].Output_Nodes_Tecplot_HighOrder(0,0,0, 1, out(), 1);

      // check
      RunRegressionTest("Order 2", CurrentFile, MasterFile, 5.0e-10, 5.0e-10);

      // Output Order 4
      MasterFile = "HighOrder_SmoothnessIndicator_Order4.dat";
      CurrentFile = "Current_HighOrder_SmoothnessIndicator_Order4.dat";
      Open_Output_File(CurrentFile);

      SolnBlk[0].Output_Nodes_Tecplot_HighOrder(0,0,0, 1, out(), 2);

      // check
      RunRegressionTest("Order 4", CurrentFile, MasterFile, 5.0e-10, 5.0e-10);

    } else {
      // Generate master files

      // Output Order 3
      MasterFile = "HighOrder_SmoothnessIndicator_Order3.dat";
      Open_Output_File(MasterFile);

      SolnBlk[0].Output_Nodes_Tecplot_HighOrder(0,0,0, 1, out(), 0);

      // Output Order 2
      MasterFile = "HighOrder_SmoothnessIndicator_Order2.dat";
      Open_Output_File(MasterFile);

      SolnBlk[0].Output_Nodes_Tecplot_HighOrder(0,0,0, 1, out(), 1);

      // Output Order 4
      MasterFile = "HighOrder_SmoothnessIndicator_Order4.dat";
      Open_Output_File(MasterFile);

      SolnBlk[0].Output_Nodes_Tecplot_HighOrder(0,0,0, 1, out(), 2);
    }

  }

  /* Test 14:*/
  template<>
  template<>
  void AdvectDiffuse2D_Quad_Block_object::test<14>()
  {

    set_test_name("Check ComputeHighOrderSolutionReconstruction()");
    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    RunRegression = ON;

    // Set input file name
    Open_Input_File("CircularAdvectionDiffusion_HighOrder_SmoothnessIndicatorStudy.in");

    // Parse the input file
    IP.Verbose() = false;
    IP.Parse_Input_File(input_file_name);

    // Create computational domain
    InitializeComputationalDomain(MeshBlk,QuadTree,
				  GlobalList_Soln_Blocks, LocalList_Soln_Blocks, 
				  SolnBlk, IP);

    // == check correct initialization
    ensure("High-order variables", SolnBlk[0].HighOrderVariables() != NULL);
    ensure_equals("Main High-order ", SolnBlk[0].HighOrderVariable(0).RecOrder(), 3);
    ensure_equals("Second High-order ", SolnBlk[0].HighOrderVariable(1).RecOrder(), 2);
    ensure_equals("Third High-order ", SolnBlk[0].HighOrderVariable(2).RecOrder(), 4);
    
    // Apply initial condition
    ICs(SolnBlk,LocalList_Soln_Blocks,IP);

    // Reconstruct the solution
    SolnBlk[0].HighOrderVariable(0).ComputeHighOrderSolutionReconstruction(SolnBlk[0],
									   IP.Limiter());
    SolnBlk[0].HighOrderVariable(1).ComputeHighOrderSolutionReconstruction(SolnBlk[0],
									   IP.Limiter());
    SolnBlk[0].HighOrderVariable(2).ComputeHighOrderSolutionReconstruction(SolnBlk[0],
									   IP.Limiter());

    if (RunRegression) {

      //===== Check interior nodal values
      // Output Order 3
      MasterFile = "HighOrder_SolutionReconstruction_Order3.dat";
      CurrentFile = "Current_HighOrder_SolutionReconstruction_Order3.dat";
      Open_Output_File(CurrentFile);

      SolnBlk[0].Output_Nodes_Tecplot_HighOrder(0,0,0, 1, out(), 0);

      // check
      RunRegressionTest("Order 3", CurrentFile, MasterFile, 5.0e-10, 5.0e-10);

      // Output Order 2
      MasterFile = "HighOrder_SolutionReconstruction_Order2.dat";
      CurrentFile = "Current_HighOrder_SolutionReconstruction_Order2.dat";
      Open_Output_File(CurrentFile);

      SolnBlk[0].Output_Nodes_Tecplot_HighOrder(0,0,0, 1, out(), 1);

      // check
      RunRegressionTest("Order 2", CurrentFile, MasterFile, 5.0e-10, 5.0e-10);

      // Output Order 4
      MasterFile = "HighOrder_SolutionReconstruction_Order4.dat";
      CurrentFile = "Current_HighOrder_SolutionReconstruction_Order4.dat";
      Open_Output_File(CurrentFile);

      SolnBlk[0].Output_Nodes_Tecplot_HighOrder(0,0,0, 1, out(), 2);

      // check
      RunRegressionTest("Order 4", CurrentFile, MasterFile, 5.0e-10, 5.0e-10);

    } else {
      // Generate master files

      // Output Order 3
      MasterFile = "HighOrder_SolutionReconstruction_Order3.dat";
      Open_Output_File(MasterFile);

      SolnBlk[0].Output_Nodes_Tecplot_HighOrder(0,0,0, 1, out(), 0);

      // Output Order 2
      MasterFile = "HighOrder_SolutionReconstruction_Order2.dat";
      Open_Output_File(MasterFile);

      SolnBlk[0].Output_Nodes_Tecplot_HighOrder(0,0,0, 1, out(), 1);

      // Output Order 4
      MasterFile = "HighOrder_SolutionReconstruction_Order4.dat";
      Open_Output_File(MasterFile);

      SolnBlk[0].Output_Nodes_Tecplot_HighOrder(0,0,0, 1, out(), 2);
    }

  }

  /* Test 15:*/
  template<>
  template<>
  void AdvectDiffuse2D_Quad_Block_object::test<15>()
  {

    set_test_name("Check dUdt_Multistage_Explicit_HighOrder()");
    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    RunRegression = ON;

    // Error norms
    double L1, L2, LMax;
    double L1_M, L2_M, LMax_M;	// master errors

    // Set input file name
    Open_Input_File("HighOrder_Residual_Study.in");

    // Parse the input file
    IP.Verbose() = false;
    IP.Parse_Input_File(input_file_name);

    // Create computational domain
    InitializeComputationalDomain(MeshBlk,QuadTree,
				  GlobalList_Soln_Blocks, LocalList_Soln_Blocks, 
				  SolnBlk, IP);

    // == check correct initialization
    ensure("High-order variables", SolnBlk[0].HighOrderVariables() != NULL);
    ensure_equals("Main High-order ", SolnBlk[0].HighOrderVariable(0).RecOrder(), 4);
    ensure_equals("2nd High-order " , SolnBlk[0].HighOrderVariable(1).RecOrder(), 1);
    ensure_equals("3rd High-order " , SolnBlk[0].HighOrderVariable(2).RecOrder(), 2);
    ensure_equals("4th High-order " , SolnBlk[0].HighOrderVariable(3).RecOrder(), 3);
    
    // Apply initial condition
    ICs(SolnBlk,LocalList_Soln_Blocks,IP);

    // Set local time step
    SetLocalTimeStepToValue(SolnBlk,
			    LocalList_Soln_Blocks,
			    1.0);

    // Compute integral of the RHS term and write it to the k_residual = 2
    ComputeEquationRightHandSideTerm(SolnBlk[0], IP, 2);

    // ========= Compute with HighOrderVariable(0) ========

    // Compute residuals for stage 1
    SolnBlk[0].dUdt_Multistage_Explicit_HighOrder(1,IP);

    // Compute residual errors
    ComputeResidualErrors(SolnBlk[0], 0, 2, L1, L2, LMax);

    // === check errors
    L1_M = 1.893420021730188e-06; L2_M = 3.843827746933967e-06; LMax_M = 1.646626337714271e-05;
    ensure_distance("L1, k=4"  , L1, L1_M, AcceptedError(L1_M, 1.0e-7) );
    ensure_distance("L2, k=4"  , L2, L2_M, AcceptedError(L2_M, 1.0e-7) );
    ensure_distance("LMax, k=4", LMax, LMax_M, AcceptedError(LMax_M, 1.0e-7) );


    // ========= Compute with HighOrderVariable(1) ========

    // Compute residuals for stage 1
    SolnBlk[0].dUdt_Multistage_Explicit_HighOrder(1,IP,1);

    // Compute residual errors
    ComputeResidualErrors(SolnBlk[0], 0, 2, L1, L2, LMax);

    // === check errors
    L1_M = 0.0001285957220702368; L2_M = 0.0001736157714376252; LMax_M = 0.0005888891194135541;
    ensure_distance("L1, k=1"  , L1, L1_M, AcceptedError(L1_M, 1.0e-7) );
    ensure_distance("L2, k=1"  , L2, L2_M, AcceptedError(L2_M, 1.0e-7) );
    ensure_distance("LMax, k=1", LMax, LMax_M, AcceptedError(LMax_M, 1.0e-7) );

    // ========= Compute with HighOrderVariable(2) ========

    // Compute residuals for stage 1
    SolnBlk[0].dUdt_Multistage_Explicit_HighOrder(1,IP,2);

    // Compute residual errors
    ComputeResidualErrors(SolnBlk[0], 0, 2, L1, L2, LMax);

    // === check errors
    L1_M = 0.0002726793258866349; L2_M = 0.0004199290419628109; LMax_M = 0.001404192919094511;
    ensure_distance("L1, k=2"  , L1, L1_M, AcceptedError(L1_M, 1.0e-7) );
    ensure_distance("L2, k=2"  , L2, L2_M, AcceptedError(L2_M, 1.0e-7) );
    ensure_distance("LMax, k=2", LMax, LMax_M, AcceptedError(LMax_M, 1.0e-7) );

    // ========= Compute with HighOrderVariable(3) ========

    // Compute residuals for stage 1
    SolnBlk[0].dUdt_Multistage_Explicit_HighOrder(1,IP,3);

    // Compute residual errors
    ComputeResidualErrors(SolnBlk[0], 0, 2, L1, L2, LMax);

    // === check errors
    L1_M = 4.32976294122695e-06; L2_M = 7.498140474432567e-06; LMax_M = 3.275058817210977e-05;
    ensure_distance("L1, k=3"  , L1, L1_M, AcceptedError(L1_M, 1.0e-7) );
    ensure_distance("L2, k=3"  , L2, L2_M, AcceptedError(L2_M, 1.0e-7) );
    ensure_distance("LMax, k=3", LMax, LMax_M, AcceptedError(LMax_M, 1.0e-7) );

    if (RunRegression == OFF){ 
      // Print errors
      cout << endl
	   << SolnBlk[0].ICu - SolnBlk[0].ICl + 1 << "x" <<  SolnBlk[0].JCu - SolnBlk[0].JCl + 1 << endl
	   << "L1_Norm = " << setprecision(16) << L1 << endl
	   << "L2_Norm = " << setprecision(16) << L2 << endl
	   << "Max_Norm = " << setprecision(16) << LMax << endl;

      // Output solution to check residual errors
      CurrentFile = "Current_HighOrder_Residual_Study.dat";
      Open_Output_File(CurrentFile);
      
      SolnBlk[0].Output_Nodes_Tecplot_HighOrder(0,0,0, 1, out(), 0);
    }

  }

  /* Test 16:*/
  template<>
  template<>
  void AdvectDiffuse2D_Quad_Block_object::test<16>()
  {

    set_test_name("Check dUdt_Residual_Evaluation_HighOrder()");
    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    RunRegression = ON;

    // Error norms
    double L1, L2, LMax;
    double L1_M, L2_M, LMax_M;	// master errors

    // Set input file name
    Open_Input_File("HighOrder_Residual_Study.in");

    // Parse the input file
    IP.Verbose() = false;
    IP.Parse_Input_File(input_file_name);

    // Create computational domain
    InitializeComputationalDomain(MeshBlk,QuadTree,
				  GlobalList_Soln_Blocks, LocalList_Soln_Blocks, 
				  SolnBlk, IP);

    // == check correct initialization
    ensure("High-order variables", SolnBlk[0].HighOrderVariables() != NULL);
    ensure_equals("Main High-order ", SolnBlk[0].HighOrderVariable(0).RecOrder(), 4);
    ensure_equals("2nd High-order " , SolnBlk[0].HighOrderVariable(1).RecOrder(), 1);
    ensure_equals("3rd High-order " , SolnBlk[0].HighOrderVariable(2).RecOrder(), 2);
    ensure_equals("4th High-order " , SolnBlk[0].HighOrderVariable(3).RecOrder(), 3);
    
    // Apply initial condition
    ICs(SolnBlk,LocalList_Soln_Blocks,IP);

    // Set local time step
    SetLocalTimeStepToValue(SolnBlk,
			    LocalList_Soln_Blocks,
			    1.0);

    // Compute integral of the RHS term and write it to the k_residual = 2
    ComputeEquationRightHandSideTerm(SolnBlk[0], IP, 2);

    // ========= Compute with HighOrderVariable(0) ========

    // Compute residuals for stage 1
    SolnBlk[0].dUdt_Residual_Evaluation_HighOrder(IP);

    // Compute residual errors
    ComputeResidualErrors(SolnBlk[0], 0, 2, L1, L2, LMax);

    // === check errors
    L1_M = 1.893420021730188e-06; L2_M = 3.843827746933967e-06; LMax_M = 1.646626337714271e-05;
    ensure_distance("L1, k=4"  , L1, L1_M, AcceptedError(L1_M, 1.0e-7) );
    ensure_distance("L2, k=4"  , L2, L2_M, AcceptedError(L2_M, 1.0e-7) );
    ensure_distance("LMax, k=4", LMax, LMax_M, AcceptedError(LMax_M, 1.0e-7) );


    // ========= Compute with HighOrderVariable(1) ========

    // Compute residuals for stage 1
    SolnBlk[0].dUdt_Residual_Evaluation_HighOrder(IP,1);

    // Compute residual errors
    ComputeResidualErrors(SolnBlk[0], 0, 2, L1, L2, LMax);

    // === check errors
    L1_M = 0.0001285957220702368; L2_M = 0.0001736157714376252; LMax_M = 0.0005888891194135541;
    ensure_distance("L1, k=1"  , L1, L1_M, AcceptedError(L1_M, 1.0e-7) );
    ensure_distance("L2, k=1"  , L2, L2_M, AcceptedError(L2_M, 1.0e-7) );
    ensure_distance("LMax, k=1", LMax, LMax_M, AcceptedError(LMax_M, 1.0e-7) );

    // ========= Compute with HighOrderVariable(2) ========

    // Compute residuals for stage 1
    SolnBlk[0].dUdt_Residual_Evaluation_HighOrder(IP,2);

    // Compute residual errors
    ComputeResidualErrors(SolnBlk[0], 0, 2, L1, L2, LMax);

    // === check errors
    L1_M = 0.0002726793258866349; L2_M = 0.0004199290419628109; LMax_M = 0.001404192919094511;
    ensure_distance("L1, k=2"  , L1, L1_M, AcceptedError(L1_M, 1.0e-7) );
    ensure_distance("L2, k=2"  , L2, L2_M, AcceptedError(L2_M, 1.0e-7) );
    ensure_distance("LMax, k=2", LMax, LMax_M, AcceptedError(LMax_M, 1.0e-7) );

    // ========= Compute with HighOrderVariable(3) ========

    // Compute residuals for stage 1
    SolnBlk[0].dUdt_Residual_Evaluation_HighOrder(IP,3);

    // Compute residual errors
    ComputeResidualErrors(SolnBlk[0], 0, 2, L1, L2, LMax);

    // === check errors
    L1_M = 4.32976294122695e-06; L2_M = 7.498140474432567e-06; LMax_M = 3.275058817210977e-05;
    ensure_distance("L1, k=3"  , L1, L1_M, AcceptedError(L1_M, 1.0e-7) );
    ensure_distance("L2, k=3"  , L2, L2_M, AcceptedError(L2_M, 1.0e-7) );
    ensure_distance("LMax, k=3", LMax, LMax_M, AcceptedError(LMax_M, 1.0e-7) );

    if (RunRegression == OFF){ 
      // Print errors
      cout << endl
	   << SolnBlk[0].ICu - SolnBlk[0].ICl + 1 << "x" <<  SolnBlk[0].JCu - SolnBlk[0].JCl + 1 << endl
	   << "L1_Norm = " << setprecision(16) << L1 << endl
	   << "L2_Norm = " << setprecision(16) << L2 << endl
	   << "Max_Norm = " << setprecision(16) << LMax << endl;

      // Output solution to check residual errors
      CurrentFile = "Current_HighOrder_Residual_Study.dat";
      Open_Output_File(CurrentFile);
      
      SolnBlk[0].Output_Nodes_Tecplot_HighOrder(0,0,0, 1, out(), 0);
    }

  }


  /* Test 17:*/
  template<>
  template<>
  void AdvectDiffuse2D_Quad_Block_object::test<17>()
  {

    set_test_name("Check AnalyseCellFaces");
    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    RunRegression = ON;

    // Set input file name
    Open_Input_File("HighOrder_CurvedBoundaries_Residual_Study.in");

    // Parse the input file
    IP.Verbose() = false;
    IP.Parse_Input_File(input_file_name);

    // Create computational domain
    InitializeComputationalDomain(MeshBlk,QuadTree,
				  GlobalList_Soln_Blocks, LocalList_Soln_Blocks, 
				  SolnBlk, IP);

    // Study interior cell
    SolnBlk[0].Grid.Integration.AnalyseCellFaces(8,10);

    ensure_equals("W face, Interior", SolnBlk[0].Grid.Integration.getWestFaceInfo(), false);
    ensure_equals("S face, Interior", SolnBlk[0].Grid.Integration.getSouthFaceInfo(), false);
    ensure_equals("E face, Interior", SolnBlk[0].Grid.Integration.getEastFaceInfo(), false);
    ensure_equals("N face, Interior", SolnBlk[0].Grid.Integration.getNorthFaceInfo(), false);
    ensure_equals("Face-block, Interior", SolnBlk[0].Grid.Integration.getFaceBlockEdgeCorrelation(), true);

    // Study SW corner
    SolnBlk[0].Grid.Integration.AnalyseCellFaces(4,4);

    ensure_equals("W face, I", SolnBlk[0].Grid.Integration.getWestFaceInfo(), true);
    ensure_equals("S face, I", SolnBlk[0].Grid.Integration.getSouthFaceInfo(), true);
    ensure_equals("E face, I", SolnBlk[0].Grid.Integration.getEastFaceInfo(), false);
    ensure_equals("N face, I", SolnBlk[0].Grid.Integration.getNorthFaceInfo(), false);
    ensure_equals("Face-block, I", SolnBlk[0].Grid.Integration.getFaceBlockEdgeCorrelation(), true);

    // Study ghost cell near West block edge
    SolnBlk[0].Grid.Integration.AnalyseCellFaces(3,4);

    ensure_equals("W face, II", SolnBlk[0].Grid.Integration.getWestFaceInfo(), false);
    ensure_equals("S face, II", SolnBlk[0].Grid.Integration.getSouthFaceInfo(), false);
    ensure_equals("E face, II", SolnBlk[0].Grid.Integration.getEastFaceInfo(), true);
    ensure_equals("N face, II", SolnBlk[0].Grid.Integration.getNorthFaceInfo(), false);
    ensure_equals("Face-block, II", SolnBlk[0].Grid.Integration.getFaceBlockEdgeCorrelation(), false);

    // Study SE corner
    SolnBlk[0].Grid.Integration.AnalyseCellFaces(19,4);

    ensure_equals("W face, III", SolnBlk[0].Grid.Integration.getWestFaceInfo(), false);
    ensure_equals("S face, III", SolnBlk[0].Grid.Integration.getSouthFaceInfo(), true);
    ensure_equals("E face, III", SolnBlk[0].Grid.Integration.getEastFaceInfo(), true);
    ensure_equals("N face, III", SolnBlk[0].Grid.Integration.getNorthFaceInfo(), false);
    ensure_equals("Face-block, III", SolnBlk[0].Grid.Integration.getFaceBlockEdgeCorrelation(), true);

    // Study ghost cell near East block edge
    SolnBlk[0].Grid.Integration.AnalyseCellFaces(20,4);

    ensure_equals("W face, IV", SolnBlk[0].Grid.Integration.getWestFaceInfo(), true);
    ensure_equals("S face, IV", SolnBlk[0].Grid.Integration.getSouthFaceInfo(), false);
    ensure_equals("E face, IV", SolnBlk[0].Grid.Integration.getEastFaceInfo(), false);
    ensure_equals("N face, IV", SolnBlk[0].Grid.Integration.getNorthFaceInfo(), false);
    ensure_equals("Face-block, IV", SolnBlk[0].Grid.Integration.getFaceBlockEdgeCorrelation(), false);

    // Study NE corner
    SolnBlk[0].Grid.Integration.AnalyseCellFaces(19,19);

    ensure_equals("W face, V", SolnBlk[0].Grid.Integration.getWestFaceInfo(), false);
    ensure_equals("S face, V", SolnBlk[0].Grid.Integration.getSouthFaceInfo(), false);
    ensure_equals("E face, V", SolnBlk[0].Grid.Integration.getEastFaceInfo(), true);
    ensure_equals("N face, V", SolnBlk[0].Grid.Integration.getNorthFaceInfo(), true);
    ensure_equals("Face-block, V", SolnBlk[0].Grid.Integration.getFaceBlockEdgeCorrelation(), true);

    // Study ghost cell near North block edge
    SolnBlk[0].Grid.Integration.AnalyseCellFaces(9,20);

    ensure_equals("W face, VI", SolnBlk[0].Grid.Integration.getWestFaceInfo(), false);
    ensure_equals("S face, VI", SolnBlk[0].Grid.Integration.getSouthFaceInfo(), true);
    ensure_equals("E face, VI", SolnBlk[0].Grid.Integration.getEastFaceInfo(), false);
    ensure_equals("N face, VI", SolnBlk[0].Grid.Integration.getNorthFaceInfo(), false);
    ensure_equals("Face-block, VI", SolnBlk[0].Grid.Integration.getFaceBlockEdgeCorrelation(), false);

    // Study NW corner
    SolnBlk[0].Grid.Integration.AnalyseCellFaces(4,19);

    ensure_equals("W face, VII", SolnBlk[0].Grid.Integration.getWestFaceInfo(), true);
    ensure_equals("S face, VII", SolnBlk[0].Grid.Integration.getSouthFaceInfo(), false);
    ensure_equals("E face, VII", SolnBlk[0].Grid.Integration.getEastFaceInfo(), false);
    ensure_equals("N face, VII", SolnBlk[0].Grid.Integration.getNorthFaceInfo(), true);
    ensure_equals("Face-block, VII", SolnBlk[0].Grid.Integration.getFaceBlockEdgeCorrelation(), true);

    // Study ghost cell near South block edge
    SolnBlk[0].Grid.Integration.AnalyseCellFaces(9,3);

    ensure_equals("W face, VIII", SolnBlk[0].Grid.Integration.getWestFaceInfo(), false);
    ensure_equals("S face, VIII", SolnBlk[0].Grid.Integration.getSouthFaceInfo(), false);
    ensure_equals("E face, VIII", SolnBlk[0].Grid.Integration.getEastFaceInfo(), false);
    ensure_equals("N face, VIII", SolnBlk[0].Grid.Integration.getNorthFaceInfo(), true);
    ensure_equals("Face-block, VIII", SolnBlk[0].Grid.Integration.getFaceBlockEdgeCorrelation(), false);

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

