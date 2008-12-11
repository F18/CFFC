/*!\file test_Euler2DQuadBlock.cc
  \brief Regression tests for class Euler2D_Quad_Block. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../Euler2DQuad.h"
#include "../Grid/HO_Grid2DQuadMultiBlock.h" /* Include 2D quadrilateral multiblock grid header file */
#include "../../HighOrderReconstruction/AccuracyAssessment2DMultiBlock.h"
#include "../../HighOrderReconstruction/HighOrder2D_MultiBlock.h" /* Include 2D high-order header file for multi-block level. */

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_Euler2D_Quad_Block : public TestData {

    // Local variables
  public:

    // ==== Member data =====

    // Euler2D input variables and parameters:
    Euler2D_Input_Parameters IP;
    
    /* Multi-block solution-adaptive quadrilateral mesh 
       solution variables. */
    Grid2D_Quad_MultiBlock_HO     MeshBlk;
    QuadTreeBlock_DataStructure  QuadTree;
    AdaptiveBlockResourceList    GlobalList_Soln_Blocks;
    AdaptiveBlock2D_List         LocalList_Soln_Blocks;
    Euler2D_Quad_Block   *SolnBlk;
    Euler2D_Quad_Block    Soln;
    CPUTime processor_cpu_time;

    int error_flag;
    int Status;                        // shows if the computational domain has been initialized

    //===== Member functions =====

    // Default Constructor
    Data_Euler2D_Quad_Block(void);

    // Destructor
    ~Data_Euler2D_Quad_Block(void);

    // Initialization of the computational domain
    void InitializeComputationalDomain(Grid2D_Quad_MultiBlock_HO & _MeshBlk_,
				       QuadTreeBlock_DataStructure & _QuadTree_,
				       AdaptiveBlockResourceList & _GlobalList_Soln_Blocks_,
				       AdaptiveBlock2D_List & _LocalList_Soln_Blocks_,
				       Euler2D_Quad_Block *& _SolnBlk_,
				       Euler2D_Input_Parameters & _IP_) throw(std::runtime_error);

    // Output_Block()
    void Output_Block(Euler2D_Quad_Block & SolnBlock,
		      Euler2D_Input_Parameters & _IP_){};

    // Set the local time step to value for all solution blocks.
    void SetLocalTimeStepToValue(Euler2D_Quad_Block *& _SolnBlk_,
				 AdaptiveBlock2D_List & _LocalList_Soln_Blocks_,
				 const double & ValueToBeSet);    

    // Check consistency of geometric properties for a multi-block mesh
    void CheckGeomPropertiesConsistencyOfBlocks(const Euler2D_Quad_Block &CheckedBlock,
						const int &i_Start, const int &i_End,
						const int &j_Start, const int &j_End,
						const Euler2D_Quad_Block &MasterBlock,
						const int &i_Master, const int &j_Master,
						const std::string & BaseMsg = "");

  private:
    
  };


  Data_Euler2D_Quad_Block::Data_Euler2D_Quad_Block(){
    /* Set the global path for this test suite. 
       It's a good practice to put it in the constructor of the data class in order to be set
       automatically for each individual test. Declare it relative to the /src_2D directory,
       otherwise the framework might not find the input and output files. */
    
    set_test_suite_path("Euler2D/UnitTests");

    // Initialize IP to default values
    Set_Default_Input_Parameters(IP);

    Status = OFF;
  }

  Data_Euler2D_Quad_Block::~Data_Euler2D_Quad_Block(void){
    if (Status == ON){
      SolnBlk = Deallocate(SolnBlk, IP);
      LocalList_Soln_Blocks.deallocate();
      GlobalList_Soln_Blocks.deallocate();
      QuadTree.deallocate();
    }
    HO_Grid2D_Execution_Mode::SetDefaults();
  }

  // === InitializeComputationalDomain()
  void Data_Euler2D_Quad_Block::InitializeComputationalDomain(Grid2D_Quad_MultiBlock_HO & _MeshBlk_,
							      QuadTreeBlock_DataStructure & _QuadTree_,
							      AdaptiveBlockResourceList & _GlobalList_Soln_Blocks_,
							      AdaptiveBlock2D_List & _LocalList_Soln_Blocks_,
							      Euler2D_Quad_Block *& _SolnBlk_,
							      Euler2D_Input_Parameters & _IP_)
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
  void Data_Euler2D_Quad_Block::SetLocalTimeStepToValue(Euler2D_Quad_Block *& _SolnBlk_,
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

  // === CheckGeomPropertiesConsistencyOfBlocks()
  void Data_Euler2D_Quad_Block::CheckGeomPropertiesConsistencyOfBlocks(const Euler2D_Quad_Block &CheckedBlock,
								       const int &i_Start, const int &i_End,
								       const int &j_Start, const int &j_End,
								       const Euler2D_Quad_Block &MasterBlock,
								       const int &i_Master, const int &j_Master,
								       const std::string & BaseMsg){

    int iCell,jCell;		// cell indexes for the checked block
    int iMast,jMast;		// cell indexes for the master block that corresponds to the CheckedBlock cell
    int iShift, jShift;		// i- and j-shift between the indexes of the two blocks
    bool ICond, JCond;		// indicators for how to loop over indexes
    int GMom;
    double Tolerance(1.0E-12);

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

	// === Check nodal values
	ostm() << "Node NW, cell (" << iCell << "," << jCell << "), " << BaseMsg << "\n"; 
	ensure_distance(ostm().str(), 
			CheckedBlock.Grid.nodeNW(iCell,jCell).X,
			MasterBlock.Grid.nodeNW(iMast,jMast).X,
			AcceptedError(MasterBlock.Grid.nodeNW(iMast,jMast).X,Tolerance));
	ostmClear();

	ostm() << "Node NE, cell (" << iCell << "," << jCell << "), " << BaseMsg << "\n"; 
	ensure_distance(ostm().str(), 
			CheckedBlock.Grid.nodeNE(iCell,jCell).X,
			MasterBlock.Grid.nodeNE(iMast,jMast).X,
			AcceptedError(MasterBlock.Grid.nodeNE(iMast,jMast).X,Tolerance));
	ostmClear();

	ostm() << "Node SE, cell (" << iCell << "," << jCell << "), " << BaseMsg << "\n"; 
	ensure_distance(ostm().str(), 
			CheckedBlock.Grid.nodeSE(iCell,jCell).X,
			MasterBlock.Grid.nodeSE(iMast,jMast).X,
			AcceptedError(MasterBlock.Grid.nodeSE(iMast,jMast).X,Tolerance));
	ostmClear();

	ostm() << "Node SW, cell (" << iCell << "," << jCell << "), " << BaseMsg << "\n"; 
	ensure_distance(ostm().str(), 
			CheckedBlock.Grid.nodeSW(iCell,jCell).X,
			MasterBlock.Grid.nodeSW(iMast,jMast).X,
			AcceptedError(MasterBlock.Grid.nodeSW(iMast,jMast).X,Tolerance));
	ostmClear();


	// === Check centroids
	ostm() << "Centroid, cell (" << iCell << "," << jCell << "), " << BaseMsg << "\n"; 
	ensure_distance(ostm().str(), 
			CheckedBlock.Grid.CellCentroid(iCell,jCell),
			MasterBlock.Grid.CellCentroid(iMast,jMast),
			AcceptedError(MasterBlock.Grid.CellCentroid(iMast,jMast),Tolerance));
	ostmClear();

	// === Check centroids
	ostm() << "Area, cell (" << iCell << "," << jCell << "), " << BaseMsg << "\n"; 
	ensure_distance(ostm().str(), 
			CheckedBlock.Grid.CellArea(iCell,jCell),
			MasterBlock.Grid.CellArea(iMast,jMast),
			AcceptedError(MasterBlock.Grid.CellArea(iMast,jMast),Tolerance));
	ostmClear();

	// === Check geometric moments
	for (GMom = 0; GMom<=CheckedBlock.Grid.CellGeomCoeff(iCell,jCell).LastElem(); ++GMom){
	  ostm() << "Geometric moment " << GMom << ", cell (" << iCell << "," << jCell << "), " << BaseMsg << "\n"; 
	  ensure_distance(ostm().str(), 
			  CheckedBlock.Grid.CellGeomCoeffValue(iCell,jCell,GMom),
			  MasterBlock.Grid.CellGeomCoeffValue(iMast,jMast,GMom),
			  AcceptedError(MasterBlock.Grid.CellGeomCoeffValue(iMast,jMast,GMom),Tolerance));
	  ostmClear();
	}
	
      }
    }
  }


  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. Euler2D_Quad_Block_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_Euler2D_Quad_Block> Euler2D_Quad_Block_TestSuite;
  typedef Euler2D_Quad_Block_TestSuite::object Euler2D_Quad_Block_object;


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
  void Euler2D_Quad_Block_object::test<1>()
  {

    set_test_name("High-order reconstruction");

    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    RunRegression = OFF;
 
    // Set input file name
    Open_Input_File("HO_Reconstruction.in");

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
    
    // Reconstruct solution
    SolnBlk[0].HighOrderVariable(0).ComputeHighOrderSolutionReconstruction(SolnBlk[0],
									   IP.Limiter());
   
    // Output Tecplot
    Output_Tecplot(SolnBlk,LocalList_Soln_Blocks,IP,0,0);

    Output_Cells_Tecplot(SolnBlk,LocalList_Soln_Blocks,IP,0,0);

  }


  /* Test 2:*/
  template<>
  template<>
  void Euler2D_Quad_Block_object::test<2>()
  {

    set_test_name("Consistency of geometric properties in ghost cells");

    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    int i0, j0, i1, j1, ng_I, ng_J, GMom;
    RunRegression = ON;
 
    // Set input file name
    Open_Input_File("CircularCylinderGrid.in");

    // Parse the input file
    IP.Verbose() = false;
    IP.Parse_Input_File(input_file_name);
    HighOrder2D_Input::Set_Final_Parameters(IP);
    CENO_Execution_Mode::CENO_SPEED_EFFICIENT = OFF;
    Grid2D_Quad_Block_HO::setLowOrderBoundaryRepresentation();    

    // Create computational domain
    InitializeComputationalDomain(MeshBlk,QuadTree,
				  GlobalList_Soln_Blocks, LocalList_Soln_Blocks, 
				  SolnBlk, IP);

    // Apply initial condition
    ICs(SolnBlk,LocalList_Soln_Blocks,IP);


    // Send messages between blocks
    error_flag = Send_All_Messages(SolnBlk, 
				   LocalList_Soln_Blocks,
				   NUM_COMP_VECTOR2D,
				   ON);

    if (RunRegression){

      // === Check corner geometric properties (not the last row of cells in i-direction) ===
      for (ng_I = 1; ng_I < SolnBlk[0].Nghost; ++ng_I){
	for (ng_J = 1; ng_J <= SolnBlk[0].Nghost; ++ng_J){
	  // Block 0 indexes
	  i0 = SolnBlk[0].ICu + ng_I;
	  j0 = SolnBlk[0].JCl - ng_J;
	  // Block 1 indexes
	  i1 = SolnBlk[1].ICl + ng_I - 1;
	  j1 = SolnBlk[1].JCl - ng_J;

	  // === Test properties
	  ensure_distance("Centroid",
			  SolnBlk[0].Grid.CellCentroid(i0,j0),
			  SolnBlk[1].Grid.CellCentroid(i1,j1),
			  AcceptedError(SolnBlk[0].Grid.CellCentroid(i0,j0)));
	  ensure_distance("Area",
			  SolnBlk[0].Grid.CellArea(i0,j0),
			  SolnBlk[1].Grid.CellArea(i1,j1),
			  AcceptedError(SolnBlk[0].Grid.CellArea(i0,j0)));

	  for (GMom = 0; GMom<=SolnBlk[0].Grid.CellGeomCoeff(i0,j0).LastElem(); ++GMom){
	    ensure_distance("Geometric moment",
			    SolnBlk[0].Grid.CellGeomCoeffValue(i0,j0,GMom),
			    SolnBlk[1].Grid.CellGeomCoeffValue(i1,j1,GMom),
			    AcceptedError(SolnBlk[0].Grid.CellGeomCoeffValue(i0,j0,GMom)));
	  }

	}
      }

    } // endif (RunRegression)

  }

  /* Test 3:*/
  template<>
  template<>
  void Euler2D_Quad_Block_object::test<3>()
  {

    set_test_name("Consistent geometric properties check with curved boundaries");

    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    int i0, j0, i1, j1, ng_I, ng_J, GMom;
    string ErrorMsg, BaseMsg;
    RunRegression = ON;
 
    // Set input file name
    Open_Input_File("CircularCylinderGrid.in");

    // Parse the input file
    IP.Verbose() = false;
    IP.Parse_Input_File(input_file_name);
    HighOrder2D_Input::Set_Final_Parameters(IP);
    CENO_Execution_Mode::CENO_SPEED_EFFICIENT = OFF;
    Grid2D_Quad_Block_HO::setHighOrderBoundaryRepresentation();

    // Create computational domain
    InitializeComputationalDomain(MeshBlk,QuadTree,
				  GlobalList_Soln_Blocks, LocalList_Soln_Blocks, 
				  SolnBlk, IP);

    // Apply initial condition
    ICs(SolnBlk,LocalList_Soln_Blocks,IP);


    // Send messages between blocks
    error_flag = Send_All_Messages(SolnBlk, 
				   LocalList_Soln_Blocks,
				   NUM_COMP_VECTOR2D,
				   ON);   

    if (RunRegression){

      // ==  Check nodes SE corner of Block 0
      for (ng_I = 0; ng_I < SolnBlk[0].Nghost; ++ng_I){
	for (ng_J = 0; ng_J <= SolnBlk[0].Nghost; ++ng_J){
	  // Block 0 indexes
	  i0 = SolnBlk[0].Grid.INu + ng_I;
	  j0 = SolnBlk[0].Grid.JNl - ng_J;
	  // Block 1 indexes
	  i1 = SolnBlk[1].Grid.INl + ng_I;
	  j1 = SolnBlk[1].Grid.JNl - ng_J;

	  ensure_distance("Nodes",
			  SolnBlk[0].Grid.Node[i0][j0].X,
			  SolnBlk[1].Grid.Node[i1][j1].X,
			  AcceptedError(SolnBlk[1].Grid.Node[i1][j1].X) );
			  // 	  Print_(SolnBlk[0].Grid.Node[SolnBlk[0].Grid.INu + ng_I][SolnBlk[0].Grid.JNl - ng_J]);
			  // 	  Print_(SolnBlk[1].Grid.Node[SolnBlk[1].Grid.INl + ng_I][SolnBlk[1].Grid.JNl - ng_J]);
	}
      }

      // === Check corner geometric properties (not the last row of cells in i-direction) ===
      for (ng_I = 1; ng_I < SolnBlk[0].Nghost; ++ng_I){
	for (ng_J = 1; ng_J <= SolnBlk[0].Nghost; ++ng_J){
	  // Block 0 indexes
	  i0 = SolnBlk[0].ICu + ng_I;
	  j0 = SolnBlk[0].JCl - ng_J;
	  // Block 1 indexes
	  i1 = SolnBlk[1].ICl + ng_I - 1;
	  j1 = SolnBlk[1].JCl - ng_J;

	  // === Test properties
	  ensure_distance("Centroid",
			  SolnBlk[0].Grid.CellCentroid(i0,j0),
			  SolnBlk[1].Grid.CellCentroid(i1,j1),
			  AcceptedError(SolnBlk[0].Grid.CellCentroid(i0,j0)));
	  ensure_distance("Area",
			  SolnBlk[0].Grid.CellArea(i0,j0),
			  SolnBlk[1].Grid.CellArea(i1,j1),
			  AcceptedError(SolnBlk[0].Grid.CellArea(i0,j0)));

	  for (GMom = 0; GMom<=SolnBlk[0].Grid.CellGeomCoeff(i0,j0).LastElem(); ++GMom){
	    ensure_distance("Geometric moment",
			    SolnBlk[0].Grid.CellGeomCoeffValue(i0,j0,GMom),
			    SolnBlk[1].Grid.CellGeomCoeffValue(i1,j1,GMom),
			    AcceptedError(SolnBlk[0].Grid.CellGeomCoeffValue(i0,j0,GMom)));
	  }

	}
      }

      // Check other cells on East and West sides
      for (ng_I = 1; ng_I < 4; ++ng_I){

	switch(ng_I){
	case 1:
	  // Block0 ghost cell near curved bnd with Block1 ghost cell
	  i0 = SolnBlk[0].ICu + 1;
	  j0 = SolnBlk[0].JCl - 1;
	  i1 = SolnBlk[1].ICl;
	  j1 = SolnBlk[1].JCl - 1;
	  BaseMsg = "Block0 ghost cell near curved bnd with Block1 ghost cell (1)";
	  break;

	case 2:
	  // Block0 ghost cell near curved bnd with Block1 interior cell
	  i0 = SolnBlk[0].ICu + 1;
	  j0 = SolnBlk[0].JCl;
	  i1 = SolnBlk[1].ICl;
	  j1 = SolnBlk[1].JCl;
	  BaseMsg = "Block0 ghost cell near curved bnd with Block1 interior cell (2)";
	  break;

	case 3:
	  // Block0 ghost cell near straight bnd with Block1 interior cell
	  i0 = SolnBlk[0].ICu + 1;
	  j0 = SolnBlk[0].JCl + 1;
	  i1 = SolnBlk[1].ICl;
	  j1 = SolnBlk[1].JCl + 1;
	  BaseMsg = "Block0 ghost cell near straight bnd with Block1 interior cell";
	  break;

	case 4:
	  // Block0 ghost cell near curved bnd with Block1 interior cell
	  i0 = SolnBlk[0].ICu + 2;
	  j0 = SolnBlk[0].JCl;
	  i1 = SolnBlk[1].ICl + 1;
	  j1 = SolnBlk[1].JCl;
	  BaseMsg = "Block0 ghost cell near curved bnd with Block1 interior cell (3)"; 
	  break;
	}

	// === Test properties
	ErrorMsg = "Centroid of " + BaseMsg;
	ensure_distance(ErrorMsg.c_str(),
			SolnBlk[0].Grid.CellCentroid(i0,j0),
			SolnBlk[1].Grid.CellCentroid(i1,j1),
			AcceptedError(SolnBlk[0].Grid.CellCentroid(i0,j0)));

	ErrorMsg = "Area of " + BaseMsg;
	ensure_distance(ErrorMsg.c_str(),
			SolnBlk[0].Grid.CellArea(i0,j0),
			SolnBlk[1].Grid.CellArea(i1,j1),
			AcceptedError(SolnBlk[0].Grid.CellArea(i0,j0)));
	
	ErrorMsg = "Geometric moment of " + BaseMsg;
	for (GMom = 0; GMom<=SolnBlk[0].Grid.CellGeomCoeff(i0,j0).LastElem(); ++GMom){
	  ensure_distance(ErrorMsg.c_str(),
			  SolnBlk[0].Grid.CellGeomCoeffValue(i0,j0,GMom),
			  SolnBlk[1].Grid.CellGeomCoeffValue(i1,j1,GMom),
			  AcceptedError(SolnBlk[0].Grid.CellGeomCoeffValue(i0,j0,GMom)));
	}

      }	// endfor (ng_I)

    } else {

      // Output Tecplot
      Output_Tecplot(SolnBlk,LocalList_Soln_Blocks,IP,0,0);
      Output_Cells_Tecplot(SolnBlk,LocalList_Soln_Blocks,IP,0,0);
      Output_Nodes_Tecplot(SolnBlk,LocalList_Soln_Blocks,IP,0,0);
    } // endif (RunRegression)

  }


  /* Test 4:*/
  template<>
  template<>
  void Euler2D_Quad_Block_object::test<4>()
  {

    set_test_name("Consistent geometric properties check for single block");

    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    int i0, j0, i1, j1, ng_I, ng_J, GMom;
    string ErrorMsg, BaseMsg;
    CPUTime processor_cpu_time;
    RunRegression = ON;
 
    // Set input file name
    Open_Input_File("RinglebGrid.in");

    // Parse the input file
    IP.Verbose() = false;
    IP.Parse_Input_File(input_file_name);
    HighOrder2D_Input::Set_Final_Parameters(IP);
    CENO_Execution_Mode::CENO_SPEED_EFFICIENT = OFF;
    Grid2D_Quad_Block_HO::setHighOrderBoundaryRepresentation();

    // Create computational domain
    InitializeComputationalDomain(MeshBlk,QuadTree,
				  GlobalList_Soln_Blocks, LocalList_Soln_Blocks, 
				  SolnBlk, IP);

    // Apply initial condition
    ICs(SolnBlk,LocalList_Soln_Blocks,IP);


    // Send messages between blocks
    error_flag = Send_All_Messages(SolnBlk, 
				   LocalList_Soln_Blocks,
				   NUM_COMP_VECTOR2D,
				   ON);   

    if (RunRegression){
      Write_Restart_Solution(SolnBlk, 
			     LocalList_Soln_Blocks, 
			     IP,
			     0, 0, processor_cpu_time);

      // == check
      CurrentFile = "Ringleb_Grid_blk000000.soln";
      MasterFile  = "Master_Ringleb_Grid_blk000000.soln";
      RunRegressionTest("Ringleb grid with curved bnds", CurrentFile, MasterFile, 5.0e-12, 5.0e-12);

    } else {
      // Output Tecplot
      Output_Tecplot(SolnBlk,LocalList_Soln_Blocks,IP,0,0);
      Output_Cells_Tecplot(SolnBlk,LocalList_Soln_Blocks,IP,0,0);
      Output_Nodes_Tecplot(SolnBlk,LocalList_Soln_Blocks,IP,0,0);
    }

  }


  /* Test 5:*/
  template<>
  template<>
  void Euler2D_Quad_Block_object::test<5>()
  {

    set_test_name("Consistent geometric properties check for mesh refinement");

    set_local_input_path("QuadBlockData");
    set_local_output_path("QuadBlockData");

    int i0, j0, i1, j1, ng_I, ng_J, GMom;
    int iCell, jCell;
    string ErrorMsg, BaseMsg;
    CPUTime processor_cpu_time;
    RunRegression = ON;
 
    // Set input file name
    Open_Input_File("RinglebGrid.in");

    // Parse the input file
    IP.Verbose() = false;
    IP.Parse_Input_File(input_file_name);
    HighOrder2D_Input::Set_Final_Parameters(IP);
    CENO_Execution_Mode::CENO_SPEED_EFFICIENT = OFF;
    Grid2D_Quad_Block_HO::setHighOrderBoundaryRepresentation();

    // Create computational domain
    InitializeComputationalDomain(MeshBlk,QuadTree,
				  GlobalList_Soln_Blocks, LocalList_Soln_Blocks, 
				  SolnBlk, IP);

    // Apply initial condition
    ICs(SolnBlk,LocalList_Soln_Blocks,IP);


    // Send messages between blocks
    error_flag = Send_All_Messages(SolnBlk, 
				   LocalList_Soln_Blocks,
				   NUM_COMP_VECTOR2D,
				   ON);   

    if (!error_flag) error_flag = Send_All_Messages(SolnBlk, 
						    LocalList_Soln_Blocks,
						    NUM_VAR_EULER2D,
						    OFF);

    /* Prescribe boundary data consistent with initial data. */
    BCs(SolnBlk, LocalList_Soln_Blocks, IP);

    // Perform uniform AMR
    error_flag = Uniform_AMR(SolnBlk,
			     IP,
			     QuadTree,
			     GlobalList_Soln_Blocks,
			     LocalList_Soln_Blocks);

    MasterFile = "RefinedRinglebMesh_GeometricProperties.dat";
    CurrentFile = "Current_RefinedRinglebMesh_GeometricProperties.dat";

    if (RunRegression){

      Open_Output_File(CurrentFile);

      // Output NW corner of Block 2
      out() << "NW corner of Block 2\n";
      for (iCell = SolnBlk[2].ICl - 2; iCell<= SolnBlk[2].ICl + 1; ++iCell){
	for (jCell = SolnBlk[2].JCu - 1; jCell<= SolnBlk[2].JCu + 2; ++jCell){
	  out() << SolnBlk[2].Grid.Cell[iCell][jCell];
	}
      }
    
      // Output NE corner of Block 3
      out() << "NE corner of Block 3\n";
      for (iCell = SolnBlk[3].ICu - 1; iCell<= SolnBlk[3].ICu + 2; ++iCell){
	for (jCell = SolnBlk[3].JCu - 1; jCell<= SolnBlk[3].JCu + 2; ++jCell){
	  out() << SolnBlk[3].Grid.Cell[iCell][jCell];
	}
      }

      // Output SE corner of Block 1
      out() << "SE corner of Block 1\n";
      for (iCell = SolnBlk[1].ICu - 1; iCell<= SolnBlk[1].ICu + 2; ++iCell){
	for (jCell = SolnBlk[1].JCl - 2; jCell<= SolnBlk[1].JCl + 1; ++jCell){
	  out() << SolnBlk[1].Grid.Cell[iCell][jCell];
	}
      }

      // Output SW corner of Block 0
      out() << "SW corner of Block 0\n";
      for (iCell = SolnBlk[0].ICl - 2; iCell<= SolnBlk[0].ICl + 1; ++iCell){
	for (jCell = SolnBlk[0].JCl - 2; jCell<= SolnBlk[0].JCl + 1; ++jCell){
	  out() << SolnBlk[0].Grid.Cell[iCell][jCell];
	}
      }

      RunRegressionTest("Check geometric properties of cells in domain corners",
			CurrentFile,
			MasterFile,
			5.0e-12, 5.0e-12);
      
      //==== Check correlations between the cells of adjacent blocks ====

      // Block 0
      CheckGeomPropertiesConsistencyOfBlocks(SolnBlk[0],
					     SolnBlk[0].ICl-1, SolnBlk[0].ICl-SolnBlk[0].Nghost,
					     SolnBlk[0].JCu+1, SolnBlk[0].JCu+SolnBlk[0].Nghost,
					     SolnBlk[2],
					     SolnBlk[2].ICl-1, SolnBlk[2].JCl,
					     "Check Block 0 with West ghost cells of Block 2");

      CheckGeomPropertiesConsistencyOfBlocks(SolnBlk[0],
					     SolnBlk[0].ICl, SolnBlk[0].ICu,
					     SolnBlk[0].JCu+1, SolnBlk[0].JCu+SolnBlk[0].Nghost,
					     SolnBlk[2],
					     SolnBlk[2].ICl, SolnBlk[2].JCl,
					     "Check Block 0 with South interior cells of Block 2");

      CheckGeomPropertiesConsistencyOfBlocks(SolnBlk[0],
					     SolnBlk[0].ICu+1, SolnBlk[0].ICu+SolnBlk[0].Nghost,
					     SolnBlk[0].JCu+1, SolnBlk[0].JCu+SolnBlk[0].Nghost,
					     SolnBlk[3],
					     SolnBlk[3].ICl  , SolnBlk[3].JCl,
					     "Check Block0 with interior SW corner of Block 3");

      CheckGeomPropertiesConsistencyOfBlocks(SolnBlk[0],
					     SolnBlk[0].ICu+1, SolnBlk[0].ICu+SolnBlk[0].Nghost,
					     SolnBlk[0].JCl, SolnBlk[0].JCu,
					     SolnBlk[1],
					     SolnBlk[1].ICl, SolnBlk[1].JCl,
					     "Check Block0 with West interior cells of Block 1");

      CheckGeomPropertiesConsistencyOfBlocks(SolnBlk[0],
					     SolnBlk[0].ICu+1, SolnBlk[0].ICu+SolnBlk[0].Nghost,
					     SolnBlk[0].JCl-1, SolnBlk[0].JCl-SolnBlk[0].Nghost,
					     SolnBlk[1],
					     SolnBlk[1].ICl, SolnBlk[1].JCl-1,
					     "Check Block0 with South ghost cells of Block 1");


      // Block 1
      CheckGeomPropertiesConsistencyOfBlocks(SolnBlk[1],
					     SolnBlk[1].ICl-1, SolnBlk[1].ICl-SolnBlk[1].Nghost,
					     SolnBlk[1].JCl-1, SolnBlk[1].JCl-SolnBlk[1].Nghost,
					     SolnBlk[0],
					     SolnBlk[0].ICu, SolnBlk[0].JCl-1,
					     "Check Block 1 with South ghost cells of Block 0");

      CheckGeomPropertiesConsistencyOfBlocks(SolnBlk[1],
					     SolnBlk[1].ICl-1, SolnBlk[1].ICl-SolnBlk[1].Nghost,
					     SolnBlk[1].JCl, SolnBlk[1].JCu,
					     SolnBlk[0],
					     SolnBlk[0].ICu, SolnBlk[0].JCl,
					     "Check Block 1 with interior East cells of Block 0");

      CheckGeomPropertiesConsistencyOfBlocks(SolnBlk[1],
					     SolnBlk[1].ICl-1, SolnBlk[1].ICl-SolnBlk[1].Nghost,
					     SolnBlk[1].JCu+1, SolnBlk[1].JCu+SolnBlk[1].Nghost,
					     SolnBlk[2],
					     SolnBlk[2].ICu, SolnBlk[2].JCl,
					     "Check Block 1 with interior SE corner of Block 2");

      CheckGeomPropertiesConsistencyOfBlocks(SolnBlk[1],
					     SolnBlk[1].ICl, SolnBlk[1].ICu,
					     SolnBlk[1].JCu+1, SolnBlk[1].JCu+SolnBlk[1].Nghost,
					     SolnBlk[3],
					     SolnBlk[3].ICl, SolnBlk[3].JCl,
					     "Check Block 1 with interior South cells of Block 3");

      CheckGeomPropertiesConsistencyOfBlocks(SolnBlk[1],
					     SolnBlk[1].ICu+1, SolnBlk[1].ICu+SolnBlk[1].Nghost,
					     SolnBlk[1].JCu+1, SolnBlk[1].JCu+SolnBlk[1].Nghost,
					     SolnBlk[3],
					     SolnBlk[3].ICu+1, SolnBlk[3].JCl,
					     "Check Block 1 with East ghost cells of Block 3");

      // Block 2
      CheckGeomPropertiesConsistencyOfBlocks(SolnBlk[2],
					     SolnBlk[2].ICu+1, SolnBlk[2].ICu+SolnBlk[2].Nghost,
					     SolnBlk[2].JCu+1, SolnBlk[2].JCu+SolnBlk[2].Nghost,
					     SolnBlk[3],
					     SolnBlk[3].ICl, SolnBlk[3].JCu+1,
					     "Check Block 2 with North ghost cells of Block 3");
    
      CheckGeomPropertiesConsistencyOfBlocks(SolnBlk[2],
					     SolnBlk[2].ICu+1, SolnBlk[2].ICu+SolnBlk[2].Nghost,
					     SolnBlk[2].JCl, SolnBlk[2].JCu,
					     SolnBlk[3],
					     SolnBlk[3].ICl, SolnBlk[3].JCl,
					     "Check Block 2 with interior West cells of Block 3");

      CheckGeomPropertiesConsistencyOfBlocks(SolnBlk[2],
					     SolnBlk[2].ICu+1, SolnBlk[2].ICu+SolnBlk[2].Nghost,
					     SolnBlk[2].JCl-1, SolnBlk[2].JCl-SolnBlk[2].Nghost,
					     SolnBlk[1],
					     SolnBlk[1].ICl, SolnBlk[1].JCu,
					     "Check Block 2 with interior NW corner of Block 1");

      CheckGeomPropertiesConsistencyOfBlocks(SolnBlk[2],
					     SolnBlk[2].ICl, SolnBlk[2].ICu,
					     SolnBlk[2].JCl-1, SolnBlk[2].JCl-SolnBlk[2].Nghost,
					     SolnBlk[0],
					     SolnBlk[0].ICl, SolnBlk[0].JCu,
					     "Check Block 2 with interior North cells of Block 0");

      CheckGeomPropertiesConsistencyOfBlocks(SolnBlk[2],
					     SolnBlk[2].ICl-1, SolnBlk[2].ICl-SolnBlk[2].Nghost,
					     SolnBlk[2].JCl-1, SolnBlk[2].JCl-SolnBlk[2].Nghost,
					     SolnBlk[0],
					     SolnBlk[0].ICl-1, SolnBlk[0].JCu,
					     "Check Block 2 with West ghost cells of Block 0");

      // Block 3
      CheckGeomPropertiesConsistencyOfBlocks(SolnBlk[3],
					     SolnBlk[3].ICl-1, SolnBlk[3].ICl-SolnBlk[3].Nghost,
					     SolnBlk[3].JCu+1, SolnBlk[3].JCu+SolnBlk[3].Nghost,
					     SolnBlk[2],
					     SolnBlk[2].ICu, SolnBlk[2].JCu+1,
					     "Check Block 3 with North ghost cells of Block 2");
    
      CheckGeomPropertiesConsistencyOfBlocks(SolnBlk[3],
					     SolnBlk[3].ICl-1, SolnBlk[3].ICl-SolnBlk[3].Nghost,
					     SolnBlk[3].JCl, SolnBlk[3].JCu,
					     SolnBlk[2],
					     SolnBlk[2].ICu, SolnBlk[2].JCl,
					     "Check Block 3 with interior East cells of Block 2");

      CheckGeomPropertiesConsistencyOfBlocks(SolnBlk[3],
					     SolnBlk[3].ICl-1, SolnBlk[3].ICl-SolnBlk[3].Nghost,
					     SolnBlk[3].JCl-1, SolnBlk[3].JCl-SolnBlk[3].Nghost,
					     SolnBlk[0],
					     SolnBlk[0].ICu, SolnBlk[0].JCu,
					     "Check Block 3 with interior NE corner of Block 0");

      CheckGeomPropertiesConsistencyOfBlocks(SolnBlk[3],
					     SolnBlk[3].ICl, SolnBlk[3].ICu,
					     SolnBlk[3].JCl-1, SolnBlk[3].JCl-SolnBlk[3].Nghost,
					     SolnBlk[1],
					     SolnBlk[1].ICl, SolnBlk[1].JCu,
					     "Check Block 3 with interior North cells of Block 1");

      CheckGeomPropertiesConsistencyOfBlocks(SolnBlk[3],
					     SolnBlk[3].ICu+1, SolnBlk[3].ICu+SolnBlk[3].Nghost,
					     SolnBlk[3].JCl-1, SolnBlk[3].JCl-SolnBlk[3].Nghost,
					     SolnBlk[1],
					     SolnBlk[1].ICu+1, SolnBlk[1].JCu,
					     "Check Block 3 with East ghost cells of Block 1");

    } else {
      
      // Output Tecplot
      Output_Tecplot(SolnBlk,LocalList_Soln_Blocks,IP,0,0);
      Output_Cells_Tecplot(SolnBlk,LocalList_Soln_Blocks,IP,0,0);
      Output_Nodes_Tecplot(SolnBlk,LocalList_Soln_Blocks,IP,0,0);
      
      Open_Output_File(MasterFile);

      // Output NW corner of Block 2
      out() << "NW corner of Block 2\n";
      for (iCell = SolnBlk[2].ICl - 2; iCell<= SolnBlk[2].ICl + 1; ++iCell){
	for (jCell = SolnBlk[2].JCu - 1; jCell<= SolnBlk[2].JCu + 2; ++jCell){
	  out() << SolnBlk[2].Grid.Cell[iCell][jCell];
	}
      }
    
      // Output NE corner of Block 3
      out() << "NE corner of Block 3\n";
      for (iCell = SolnBlk[3].ICu - 1; iCell<= SolnBlk[3].ICu + 2; ++iCell){
	for (jCell = SolnBlk[3].JCu - 1; jCell<= SolnBlk[3].JCu + 2; ++jCell){
	  out() << SolnBlk[3].Grid.Cell[iCell][jCell];
	}
      }

      // Output SE corner of Block 1
      out() << "SE corner of Block 1\n";
      for (iCell = SolnBlk[1].ICu - 1; iCell<= SolnBlk[1].ICu + 2; ++iCell){
	for (jCell = SolnBlk[1].JCl - 2; jCell<= SolnBlk[1].JCl + 1; ++jCell){
	  out() << SolnBlk[1].Grid.Cell[iCell][jCell];
	}
      }

      // Output SW corner of Block 0
      out() << "SW corner of Block 0\n";
      for (iCell = SolnBlk[0].ICl - 2; iCell<= SolnBlk[0].ICl + 1; ++iCell){
	for (jCell = SolnBlk[0].JCl - 2; jCell<= SolnBlk[0].JCl + 1; ++jCell){
	  out() << SolnBlk[0].Grid.Cell[iCell][jCell];
	}
      }

    } // endif (RunRegression)

  }

}



// Test suite constructor
tut::Euler2D_Quad_Block_TestSuite Euler2D_Quad_BlockTestSuite("Class:Euler2D_Quad_Block");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/


