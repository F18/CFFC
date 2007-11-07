/*!\file test_Levermore1DVector.cc
  \brief Regression tests for template class Levermore1D_Vector datatype. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../Levermore1DVector.h"

/* define useful constants for tests */
#define LEVERMORE1D_VECTOR_LENGTH    50

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_Levermore1DVector : public TestData {

    // Local variables
  public:

    // Constructor
    Data_Levermore1DVector(){ }

  private:
    
  };

  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. 'NAME'_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_Levermore1DVector> Levermore1DVector_TestSuite;
  typedef Levermore1DVector_TestSuite::object Levermore1DVector_object;


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
  void Levermore1DVector_object::test<1>()
  {
    set_test_name("Constructors");

    int i(0);
    Levermore1D_Vector<LEVERMORE1D_VECTOR_LENGTH> V1;
    
    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      ensure_distance("Default Constructor, set it to zero.",V1[i],0.0,tol);
      //set to a value for copy constructor test
      V1[i] = pow((double)i,1.23456) / 98.765;
    }

    Levermore1D_Vector<LEVERMORE1D_VECTOR_LENGTH> V2(V1);

    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      ensure_distance("V2=V1",V2[i],V1[i],fabs(V1[i])*tol);
    }

  }

  /* Test 2:*/
  template<>
  template<>
  void Levermore1DVector_object::test<2>()
  {
    set_test_name("Assignment Operator");

    int i(0);
    Levermore1D_Vector<LEVERMORE1D_VECTOR_LENGTH> V1, V2, V3;
    
    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      V1[i] = pow((double)i,1.23456) / 98.765;
    }

    V3 = V2 = V1;

    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      ensure_distance("V2=V1",V2[i],V1[i],fabs(V1[i])*tol);
      ensure_distance("V3=V1",V3[i],V1[i],fabs(V1[i])*tol);
    }
  }

  /* Test 3:*/
  template<>
  template<>
  void Levermore1DVector_object::test<3>()
  {
    set_test_name("Addition Operator");

    int i(0);
    double a(0.0), b(0.0);
    Levermore1D_Vector<LEVERMORE1D_VECTOR_LENGTH> V1, V2, V3;
    
    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      V1[i] = a;
      V2[i] = b;
    }

    V3 = V1 + V2;

    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      ensure_distance("V3 = V1 + V2", V3[i], a+b, fabs(a+b)*tol);
    }
  }

  /* Test 4:*/
  template<>
  template<>
  void Levermore1DVector_object::test<4>()
  {
    set_test_name("Addition and Assignment Operator");

    int i(0);
    double a(0.0), b(0.0);
    Levermore1D_Vector<LEVERMORE1D_VECTOR_LENGTH> V1, V2, V3;
    
    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      V1[i] = a;
      V2[i] = b;
    }

    V3 = (V2 += V1);

    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);

      ensure_distance("V2 += V1",V2[i],a+b,fabs(a+b)*tol);
      ensure_distance("V3 = (V2+=V1)",V3[i],a+b,fabs(a+b)*tol);
    }
  }

  /* Test 5:*/
  template<>
  template<>
  void Levermore1DVector_object::test<5>()
  {
    set_test_name("Subtraction Operator");

    int i(0);
    double a(0.0), b(0.0);
    Levermore1D_Vector<LEVERMORE1D_VECTOR_LENGTH> V1, V2, V3;
    
    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      V1[i] = a;
      V2[i] = b;
    }

    V3 = V1 - V2;

    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      ensure_distance("V3 = V1 - V2", V3[i], a-b, fabs(a-b)*tol);
    }
  }

  /* Test 6:*/
  template<>
  template<>
  void Levermore1DVector_object::test<6>()
  {
    set_test_name("Subtraction and Assignment Operator");

    int i(0);
    double a(0.0), b(0.0);
    Levermore1D_Vector<LEVERMORE1D_VECTOR_LENGTH> V1, V2, V3;
    
    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      V1[i] = a;
      V2[i] = b;
    }

    V3 = (V2 -= V1);

    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);

      ensure_distance("V2 -= V1",V2[i],b-a,fabs(b-a)*tol);
      ensure_distance("V3 = (V2-=V1)",V3[i],b-a,fabs(b-a)*tol);
    }
  }

  /* Test 7:*/
  template<>
  template<>
  void Levermore1DVector_object::test<7>()
  {
    set_test_name("Dot Product Operator");

    int i(0);
    double a(0.0), b(0.0), dot(0.0), dot_levermore(0.0);
    Levermore1D_Vector<LEVERMORE1D_VECTOR_LENGTH> V1, V2;
    
    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      dot += (a*b);
      V1[i] = a;
      V2[i] = b;
    }

    dot_levermore = V1 * V2;
    ensure_distance("dot_prod = V1*V2)", dot_levermore, dot, fabs(dot)*tol);
    dot_levermore = V2 * V1;
    ensure_distance("dot_prod = V2*V1)", dot_levermore, dot, fabs(dot)*tol);

  }

  /* Test 8:*/
  template<>
  template<>
  void Levermore1DVector_object::test<8>()
  {
    set_test_name("Term-wise Product Operator");

    int i(0);
    double a(0.0), b(0.0);
    Levermore1D_Vector<LEVERMORE1D_VECTOR_LENGTH> V1, V2, V3;
    
    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      V1[i] = a;
      V2[i] = b;
    }

    V3 = V1 ^ V2;
    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      ensure_distance("V3 = V1^V2)", V3[i], a*b, fabs(a*b)*tol);
    }

    V3 = V2 ^ V1;
    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      ensure_distance("V3 = V1^V2)", V3[i], a*b, fabs(a*b)*tol);
    }
  }

  /* Test 9:*/
  template<>
  template<>
  void Levermore1DVector_object::test<9>()
  {
    set_test_name("copy_form Function");
    int i(0);
    Levermore1D_Vector<LEVERMORE1D_VECTOR_LENGTH> V1, V2;
    
    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      V1[i] = pow((double)i,1.23456) / 98.765;
    }

    V2.copy_from(V1);

    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      ensure_distance("V2=V1",V2[i],V1[i],fabs(V1[i])*tol);
    }
  }

  /* Test 10:*/
  template<>
  template<>
  void Levermore1DVector_object::test<10>()
  {
    set_test_name("zero, one, and set_all Function");
    int i(0);
    Levermore1D_Vector<LEVERMORE1D_VECTOR_LENGTH> V1, V2, V3;

    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      V1[i] = pow((double)i,1.23456) / 98.765;
      V2[i] = pow((double)i,1.23456) / 98.765;
      V3[i] = pow((double)i,1.23456) / 98.765;
      ensure("V1!=0.0", V1[i] != 0.0);
      ensure("V2!=1.0", V1[i] != 1.0);
      ensure("V3!=2.0", V1[i] != 2.0);
    }

    V1.zero();
    V2.one();
    V3.set_all(2.0);

    for(i=1;i<=LEVERMORE1D_VECTOR_LENGTH;++i) {
      ensure_distance("V1 = 0.0",V1[i],0.0,tol);
      ensure_distance("V2 = 1.0",V2[i],1.0,tol);
      ensure_distance("V3 = 2.0",V3[i],2.0,tol);
    }
  }

  //end tests
}



// Test suite constructor
tut::Levermore1DVector_TestSuite Levermore1DVectorTestSuite("Template Class:Levermore1D_Vector");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

