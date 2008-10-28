/*!\file test_MathContainers.cc
  \brief Regression tests for math container classes (i.e. DenseMatrix). */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../Matrix.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_DenseMatrix : public TestData {

    // Local variables
  public:

    // Constructor
    Data_DenseMatrix();

  private:
    
  };

  // Constructor
  Data_DenseMatrix::Data_DenseMatrix(){
    /* Set the global path for this test suite. 
       It's a good practice to put it in the constructor of the data class in order to be set
       automatically for each individual test. Declare it relative to the /src_2D directory,
       otherwise the framework might not find the input and output files. */
    
    set_test_suite_path("/Math/UnitTests");
  }


  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. 'NAME'_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_DenseMatrix> DenseMatrix_TestSuite;
  typedef DenseMatrix_TestSuite::object DenseMatrix_object;


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
  void DenseMatrix_object::test<1>()
  {
    set_test_name("Incorporate matrix member function");
   
    // Data
    DenseMatrix A(6,5), B(4,3), Result;

    // Initialize test matrices
    A(0,0)= -3.300330e-01; A(0,1)=  5.000000e-01; A(0,2)= -5.190278e-01;
    A(0,3)=  4.037896e-01; A(0,4)= -2.534612e-01;
    A(1,0)= -4.950495e-01; A(1,1)=  5.000000e-01; A(1,2)= -3.577083e-01;
    A(1,3)=  1.912687e-01; A(1,4)= -8.326478e-02;
    A(2,0)= -9.900990e-01; A(2,1)=  5.000000e-01; A(2,2)= -2.104167e-01;
    A(2,3)=  6.375625e-02; A(2,4)= -1.627732e-02;
    A(3,0)=  9.900990e-01; A(3,1)=  5.000000e-01; A(3,2)=  2.104167e-01;
    A(3,3)=  6.375625e-02; A(3,4)=  1.627732e-02;
    A(4,0)=  4.950495e-01; A(4,1)=  5.000000e-01; A(4,2)=  3.577083e-01;
    A(4,3)=  1.912687e-01; A(4,4)=  8.326478e-02;
    A(5,0)=  2.501813e-01; A(5,1)=  3.790247e-01; A(5,2)=  3.934487e-01;
    A(5,3)=  3.060924e-01; A(5,4)= 1.921361e-01;

    B(0,0)= -2.093656e-01; B(0,1)= 2.114593e-01; B(0,2)= -1.512815e-01;
    B(1,0)= -9.900990e-01; B(1,1)= 5.000000e-01; B(1,2)= -2.104167e-01;
    B(2,0)=  9.900990e-01; B(2,1)= 5.000000e-01; B(2,2)=  2.104167e-01;
    B(3,0)=  4.950495e-01; B(3,1)= 5.000000e-01; B(3,2)= 3.577083e-01;

    // Generate Result
    Result = A;

    Result(2,1)= -2.093656e-01; Result(2,2)= 2.114593e-01; Result(2,3)= -1.512815e-01;
    Result(3,1)= -9.900990e-01; Result(3,2)= 5.000000e-01; Result(3,3)= -2.104167e-01;
    Result(4,1)=  9.900990e-01; Result(4,2)= 5.000000e-01; Result(4,3)=  2.104167e-01;
    Result(5,1)=  4.950495e-01; Result(5,2)= 5.000000e-01; Result(5,3)= 3.577083e-01;

    // Incorporate matrix B into A
    A.incorporate_matrix(2,1,B);

    // === Check 
    ensure("Compare to solution", A == Result);

  }

  /* Test 2:*/
  template<>
  template<>
  void DenseMatrix_object::test<2>()
  {
    set_test_name("Incorporate matrix member function");
   
    // Data
    DenseMatrix A(6,5), B(4,3), Result;

    // Initialize test matrices
    A(0,0)= -3.300330e-01; A(0,1)=  5.000000e-01; A(0,2)= -5.190278e-01;
    A(0,3)=  4.037896e-01; A(0,4)= -2.534612e-01;
    A(1,0)= -4.950495e-01; A(1,1)=  5.000000e-01; A(1,2)= -3.577083e-01;
    A(1,3)=  1.912687e-01; A(1,4)= -8.326478e-02;
    A(2,0)= -9.900990e-01; A(2,1)=  5.000000e-01; A(2,2)= -2.104167e-01;
    A(2,3)=  6.375625e-02; A(2,4)= -1.627732e-02;
    A(3,0)=  9.900990e-01; A(3,1)=  5.000000e-01; A(3,2)=  2.104167e-01;
    A(3,3)=  6.375625e-02; A(3,4)=  1.627732e-02;
    A(4,0)=  4.950495e-01; A(4,1)=  5.000000e-01; A(4,2)=  3.577083e-01;
    A(4,3)=  1.912687e-01; A(4,4)=  8.326478e-02;
    A(5,0)=  2.501813e-01; A(5,1)=  3.790247e-01; A(5,2)=  3.934487e-01;
    A(5,3)=  3.060924e-01; A(5,4)= 1.921361e-01;

    B(0,0)= -2.093656e-01; B(0,1)= 2.114593e-01; B(0,2)= -1.512815e-01;
    B(1,0)= -9.900990e-01; B(1,1)= 5.000000e-01; B(1,2)= -2.104167e-01;
    B(2,0)=  9.900990e-01; B(2,1)= 5.000000e-01; B(2,2)=  2.104167e-01;
    B(3,0)=  4.950495e-01; B(3,1)= 5.000000e-01; B(3,2)= 3.577083e-01;

    // Generate Result
    Result = A;

    Result(0,0)= -2.093656e-01; Result(0,1)= 2.114593e-01; Result(0,2)= -1.512815e-01;
    Result(1,0)= -9.900990e-01; Result(1,1)= 5.000000e-01; Result(1,2)= -2.104167e-01;
    Result(2,0)=  9.900990e-01; Result(2,1)= 5.000000e-01; Result(2,2)=  2.104167e-01;
    Result(3,0)=  4.950495e-01; Result(3,1)= 5.000000e-01; Result(3,2)= 3.577083e-01;

    // Incorporate matrix B into A
    A.incorporate_matrix(0,0,B);

    // === Check 
    ensure("Compare to solution", A == Result);

  }


}



// Test suite constructor
tut::DenseMatrix_TestSuite DenseMatrixTestSuite("Class:DenseMatrix");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

