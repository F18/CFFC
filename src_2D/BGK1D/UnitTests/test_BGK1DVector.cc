/*!\file test_BGK1DVector.cc
  \brief Regression tests for template class BGK1D_Vector datatype. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../BGK1DVector.h"
#include "test_BGK1D_defines.h"

namespace tut
{

  /* Define the test-specific data class and add data members
     when tests have complex or repeating creation phase. */
  class Data_BGK1DVector : public TestData {

    // Local variables
  public:

    // Constructor
    Data_BGK1DVector(){
      if(!BGK1D_Vector::length_is_set())
	BGK1D_Vector::set_length(BGK1D_TEST_VECTOR_LENGTH);
    }

  private:

  };

  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. 'NAME'_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_BGK1DVector> BGK1DVector_TestSuite;
  typedef BGK1DVector_TestSuite::object BGK1DVector_object;


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
  void BGK1DVector_object::test<1>()
  {
    set_test_name("Copy Constructor");

    int i(0);
    BGK1D_Vector V1;

    for(i=0;i<BGK1D_Vector::get_length();++i) {
      //set to a value for copy constructor test
      V1[i] = pow((double)i,1.23456) / 98.765;
    }

    BGK1D_Vector V2(V1);

    for(i=0;i<BGK1D_Vector::get_length();++i) {
      ensure_distance("V2=V1",V2[i],V1[i],fabs(V1[i])*tol);
    }

  }

  /* Test 2:*/
  template<>
  template<>
  void BGK1DVector_object::test<2>()
  {
    set_test_name("Assignment Operator");

    int i(0);
    BGK1D_Vector V1, V2, V3;

    for(i=0;i<BGK1D_Vector::get_length();++i) {
      V1[i] = pow((double)i,1.23456) / 98.765;
    }

    V3 = V2 = V1;

    for(i=0;i<BGK1D_Vector::get_length();++i) {
      ensure_distance("V2=V1",V2[i],V1[i],fabs(V1[i])*tol);
      ensure_distance("V3=V1",V3[i],V1[i],fabs(V1[i])*tol);
    }
  }

  /* Test 3:*/
  template<>
  template<>
  void BGK1DVector_object::test<3>()
  {
    set_test_name("Addition Operator");

    int i(0);
    double a(0.0), b(0.0);
    BGK1D_Vector V1, V2, V3;

    for(i=0;i<BGK1D_Vector::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      V1[i] = a;
      V2[i] = b;
    }

    V3 = V1 + V2;

    for(i=0;i<BGK1D_Vector::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      ensure_distance("V3 = V1 + V2", V3[i], a+b, fabs(a+b)*tol);
    }
  }

  /* Test 4:*/
  template<>
  template<>
  void BGK1DVector_object::test<4>()
  {
    set_test_name("Addition and Assignment Operator");

    int i(0);
    double a(0.0), b(0.0);
    BGK1D_Vector V1, V2, V3;

    for(i=0;i<BGK1D_Vector::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      V1[i] = a;
      V2[i] = b;
    }

    V3 = (V2 += V1);

    for(i=0;i<BGK1D_Vector::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);

      ensure_distance("V2 += V1",V2[i],a+b,fabs(a+b)*tol);
      ensure_distance("V3 = (V2+=V1)",V3[i],a+b,fabs(a+b)*tol);
    }
  }

  /* Test 5:*/
  template<>
  template<>
  void BGK1DVector_object::test<5>()
  {
    set_test_name("Subtraction Operator");

    int i(0);
    double a(0.0), b(0.0);
    BGK1D_Vector V1, V2, V3;

    for(i=0;i<BGK1D_Vector::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      V1[i] = a;
      V2[i] = b;
    }

    V3 = V1 - V2;

    for(i=0;i<BGK1D_Vector::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      ensure_distance("V3 = V1 - V2", V3[i], a-b, fabs(a-b)*tol);
    }
  }

  /* Test 6:*/
  template<>
  template<>
  void BGK1DVector_object::test<6>()
  {
    set_test_name("Subtraction and Assignment Operator");

    int i(0);
    double a(0.0), b(0.0);
    BGK1D_Vector V1, V2, V3;

    for(i=0;i<BGK1D_Vector::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      V1[i] = a;
      V2[i] = b;
    }

    V3 = (V2 -= V1);

    for(i=0;i<BGK1D_Vector::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);

      ensure_distance("V2 -= V1",V2[i],b-a,fabs(b-a)*tol);
      ensure_distance("V3 = (V2-=V1)",V3[i],b-a,fabs(b-a)*tol);
    }
  }

  /* Test 7:*/
  template<>
  template<>
  void BGK1DVector_object::test<7>()
  {
    set_test_name("Dot Product Operator");

    int i(0);
    double a(0.0), b(0.0), dot(0.0), dot_levermore(0.0);
    BGK1D_Vector V1, V2;

    for(i=0;i<BGK1D_Vector::get_length();++i) {
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
  void BGK1DVector_object::test<8>()
  {
    set_test_name("zero and one functions");
    int i(0);
    BGK1D_Vector V1, V2;

    for(i=0;i<BGK1D_Vector::get_length();++i) {
      V1[i] = pow((double)(i+1),1.23456) / 98.765;
      V2[i] = pow((double)(i+1),1.23456) / 98.765;
      ensure("V1!=0.0", V1[i] != 0.0);
      ensure("V2!=1.0", V2[i] != 1.0);
    }

    V1.zero();
    V2.one();

    for(i=0;i<BGK1D_Vector::get_length();++i) {
      ensure_distance("V1 == 0.0",V1[i],0.0,tol);
      ensure_distance("V2 == 1.0",V2[i],1.0,tol);
    }
  }


  //end tests
}



// Test suite constructor
tut::BGK1DVector_TestSuite BGK1DVectorTestSuite("Class:BGK1D_Vector");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

