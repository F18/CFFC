/*!\file test_Cauchy_BoundaryConditions_Double.cc
  \brief Regression tests for template class Cauchy_BCs with double datatype. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../Cauchy_BoundaryConditions.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_CauchyBCs : public TestData {

    // Local variables
  public:

    // Constructor
    Data_CauchyBCs();

    // Fill data into container
    void Initialize_Container(Cauchy_BCs<double> & BC, const int & NumOfObjects);

  private:
    
  };

  // Constructor
  Data_CauchyBCs::Data_CauchyBCs(){

    set_test_suite_path("HighOrderReconstruction/UnitTests/");
    set_local_input_path("CauchyBCs_Data");
    set_local_output_path("CauchyBCs_Data");
    
  }

    // Fill data into container
  void Data_CauchyBCs::Initialize_Container(Cauchy_BCs<double> & BC, const int & NumOfObjects){
    
    BC.allocate(NumOfObjects);

    for (int n = 0; n < NumOfObjects; ++n){
      BC.DirichletBC(n+1) = 1.0*(n+1);
      BC.NeumannBC(n+1) = 2.0*(n+1);
      BC.a(n+1) = 3.0*(n+1);
      BC.b(n+1) = 4.0*(n+1);
    }
  }

  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. 'NAME'_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_CauchyBCs> CauchyBCs_TestSuite;
  typedef CauchyBCs_TestSuite::object CauchyBCs_object;


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
  void CauchyBCs_object::test<1>()
  {
    set_test_name("Defalt Constructor");

    Cauchy_BCs<double> BCs;

    ensure("Dirichlet", BCs.DirichletBC() == NULL);
    ensure("Neumann", BCs.NeumannBC() == NULL);
    ensure("a", BCs.a() == NULL);
    ensure("b", BCs.b() == NULL);
    ensure_equals("NumLoc", BCs.NumOfPoints(), 0);
  }

  /* Test 2:*/
  template<>
  template<>
  void CauchyBCs_object::test<2>()
  {
    set_test_name("allocate()");

    Cauchy_BCs<double> BCs;

    // call allocate()
    BCs.allocate(5);

    ensure("Dirichlet", BCs.DirichletBC() != NULL);
    ensure("Neumann", BCs.NeumannBC() != NULL);
    ensure("a", BCs.a() != NULL);
    ensure("b", BCs.b() != NULL);
    ensure_equals("NumLoc", BCs.NumOfPoints(), 5);
  }

  /* Test 3:*/
  template<>
  template<>
  void CauchyBCs_object::test<3>()
  {
    set_test_name("Check memory leaks");

    Cauchy_BCs<double> BCs;

    // call allocate()
    BCs.allocate(5);

    ensure("Dirichlet", BCs.DirichletBC() != NULL);
    ensure("Neumann", BCs.NeumannBC() != NULL);
    ensure("a", BCs.a() != NULL);
    ensure("b", BCs.b() != NULL);
    ensure_equals("NumLoc", BCs.NumOfPoints(), 5);

    // call allocate()
    BCs.allocate(15);

    ensure_equals("NumLoc", BCs.NumOfPoints(), 15);

    // call deallocate()
    BCs.deallocate();

    ensure("Dirichlet", BCs.DirichletBC() == NULL);
    ensure("Neumann", BCs.NeumannBC() == NULL);
    ensure("a", BCs.a() == NULL);
    ensure("b", BCs.b() == NULL);
    ensure_equals("NumLoc", BCs.NumOfPoints(), 0);
  }

  /* Test 4:*/
  template<>
  template<>
  void CauchyBCs_object::test<4>()
  {
    set_test_name("Output operator");

    RunRegression = ON;

    Cauchy_BCs<double> BCs;
    Initialize_Container(BCs,4);

    MasterFile = "Output_Operator.dat";
    CurrentFile = "Current_Output_Operator.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);
      
      out() << BCs;

      RunRegressionTest("operator <<", CurrentFile, MasterFile, 1.0e-12);

    } else {
      // Generate the master file
      Open_Output_File(MasterFile);

      out() << BCs;
    }
  }

  /* Test 5:*/
  template<>
  template<>
  void CauchyBCs_object::test<5>()
  {
    set_test_name("Input-output operators");

    Cauchy_BCs<double> BCs;

    // initialize variable
    Initialize_Container(BCs,4);

    // == check
    Check_Input_Output_Operator(BCs);
  }

  /* Test 6:*/
  template<>
  template<>
  void CauchyBCs_object::test<6>()
  {
    set_test_name("Copy constructor");

    Cauchy_BCs<double> BCs;
    Initialize_Container(BCs,5);

    Cauchy_BCs<double> BCs_Copy(BCs);

    // == check
    for (int n = 1; n <= BCs.NumOfPoints(); ++n){
      ensure_equals("Dirichlet", BCs_Copy.DirichletBC(n), BCs.DirichletBC(n) );
      ensure_equals("Neumann", BCs_Copy.NeumannBC(n), BCs.NeumannBC(n) );
      ensure_equals("a", BCs_Copy.a(n), BCs.a(n) );
      ensure_equals("b", BCs_Copy.b(n), BCs.b(n) );
    }

  }

  /* Test 7:*/
  template<>
  template<>
  void CauchyBCs_object::test<7>()
  {
    set_test_name("Assignment operator");

    Cauchy_BCs<double> BCs, BCs_Copy;
    Initialize_Container(BCs,5);
    Initialize_Container(BCs_Copy,15);

    // call 2nd time
    Initialize_Container(BCs_Copy,15);

    // Copy
    BCs_Copy = BCs;

    // == check
    for (int n = 1; n <= BCs.NumOfPoints(); ++n){
      ensure_equals("Dirichlet", BCs_Copy.DirichletBC(n), BCs.DirichletBC(n) );
      ensure_equals("Neumann", BCs_Copy.NeumannBC(n), BCs.NeumannBC(n) );
      ensure_equals("a", BCs_Copy.a(n), BCs.a(n) );
      ensure_equals("b", BCs_Copy.b(n), BCs.b(n) );
    }

  }


}



// Test suite constructor
tut::CauchyBCs_TestSuite CauchyBCsTestSuite("Template Class:Cauchy_BCs && double");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

