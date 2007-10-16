/*!\file test_HighOrder1D_Euler1DUniformMesh.cc
  \brief Integration tests for template class HighOrder1D with Euler1D_UniformMesh. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../../HighOrderReconstruction/HighOrder1D.h"
#include "../Euler1D.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_HighOrder1D_EulerUniformMesh : public TestData {

    // Local variables
  public:
    Euler1D_UniformMesh *SolnBlk;
    CFD1D_Input_Parameters IP;

    // Member functions
    Data_HighOrder1D_EulerUniformMesh();
    ~Data_HighOrder1D_EulerUniformMesh();

    void SetUpDomain(char *Input_File_Name);

  private:
    
  };


  // Constructor
  Data_HighOrder1D_EulerUniformMesh::Data_HighOrder1D_EulerUniformMesh(): SolnBlk(NULL){

    set_test_suite_path("Euler1D/UnitTests/");
    set_local_input_path("test_HighOrder1D");
    set_local_output_path("test_HighOrder1D");

    Set_Default_Input_Parameters(IP);
    IP.Verbose() = OFF;
  }

  // Destructor
  Data_HighOrder1D_EulerUniformMesh::~Data_HighOrder1D_EulerUniformMesh(){
    Deallocate(SolnBlk);
  }


  void Data_HighOrder1D_EulerUniformMesh::SetUpDomain(char *Input_File_Name){
    // Set input file name
    Open_Input_File(Input_File_Name);

    // Parse input file
    IP.Parse_Input_File(input_file_name);
    
    /* Allocate memory for 1D Euler equation solution on
       uniform mesh. */
    if (IP.Verbose()) cout << "\n Creating memory for Euler1D solution variables.";

    if (SolnBlk != NULL){	// deallocate memory
      if (IP.Verbose()) cout << "\n Already allcoated memory.\n Deallocating Euler1D solution variables.";
      SolnBlk=Deallocate(SolnBlk);
    }
    SolnBlk=Allocate(SolnBlk,IP);
    
    if (SolnBlk == NULL){
      Message = "\n Euler1DSolvers::Allocate() Error! Probably not enough memory!";
      if (IP.Verbose()){
	cout << Message << endl;
      }
      throw runtime_error(Message.c_str());
    }

    /* Create uniform mesh. */
    if (IP.Verbose()){
      cout << "\n Creating uniform mesh.";
    }

    Grid(SolnBlk, IP.X_Min, IP.X_Max, IP.Number_of_Cells);

    /********************************************************  
     * Initialize Euler1D solution variables.               *
     ********************************************************/
   
    /* Initialize the conserved and primitive state
       solution variables. */
  
    if (IP.Verbose()){
      cout << "\n Prescribing Euler1D initial data.";
    }
    ICs(SolnBlk, 
	"AIR", 
	IP.i_ICs, 
	IP.Number_of_Cells,
	IP);
  }


  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. 'NAME'_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_HighOrder1D_EulerUniformMesh> HighOrder1D_EulerUniformMesh_TestSuite;
  typedef HighOrder1D_EulerUniformMesh_TestSuite::object HighOrder1D_EulerUniformMesh_object;


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
  void HighOrder1D_EulerUniformMesh_object::test<1>()
  {
    set_test_name("Test HighOrder1D functionality");
    set_local_input_path("SodProblem");

    IP.Verbose() = OFF;

    // Set up domain
    SetUpDomain("sod_HighOrder.in");

    // Reconstruct the solution in one of the cells
    SolnBlk[5].CellHighOrder().ComputeReconstructionPseudoInverse(SolnBlk,5);
    SolnBlk[5].CellHighOrder().ComputeUnlimitedSolutionReconstruction(SolnBlk,5,RECONSTRUCTION_CENO);
    //     Print_(SolnBlk[5].CellHighOrder().LHS());
    //     Print_(SolnBlk[5].CellHighOrder().GeomWeights());
    //     Print_(SolnBlk[5].CellHighOrder().CellDeriv());

  }

}


// Test suite constructor
tut::HighOrder1D_EulerUniformMesh_TestSuite 
HighOrder1D_EulerUniformMesh_TestSuite("Template Class:HighOrder1D && Euler1D_UniformMesh");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

