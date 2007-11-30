/*!\file test_Levermore1D_cState.cc
  \brief Regression tests for template class Levermore1D_cState datatype. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../Levermore1DState.h"
#include "test_Levermore1D_defines.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_Levermore1D_cState : public TestData {

    // Local variables
  public:

    // Constructor
    Data_Levermore1D_cState(){
      if(!Levermore1D_Vector::length_is_set())
	Levermore1D_Vector::set_length(LEVERMORE1D_VECTOR_LENGTH);
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
  typedef test_group<Data_Levermore1D_cState> Levermore1D_cState_TestSuite;
  typedef Levermore1D_cState_TestSuite::object Levermore1D_cState_object;


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
  void Levermore1D_cState_object::test<1>()
  {
    set_test_name("Copy Constructors");

    int i(0);
    Levermore1D_cState U1;
    
    for(i=1;i<=Levermore1D_Vector::get_length();++i) {
      //set to a value for copy constructor test
      U1[i] = pow((double)i,4.54321) / 98.765;
    }

    Levermore1D_cState U2(U1);

    for(i=1;i<=Levermore1D_Vector::get_length();++i) {
      ensure_distance("U2=U1",U2[i],U1[i],fabs(U1[i])*tol);
    }

  }

  /* Test 2:*/
  template<>
  template<>
  void Levermore1D_cState_object::test<2>()
  {
    set_test_name("Constructor from pState");

    int i(0);
    Levermore1D_pState W;

    double rho(1.005), u(-210.1), p(75663.0), q(-765.4e4), r(1.2345e8),
      rm5(3.222e9), rm6(2.10e12);

    double u1, u2, u3, u4, u5, u6, u7;

    u1 = rho;
    u2 = rho*u;
    u3 = rho*u*u + p;
    u4 = rho*u*u*u + 3.0*u*p + q;
    u5 = rho*u*u*u*u + 6.0*u*u*p + 4.0*u*q + r;
    u6 = rho*u*u*u*u*u + 10.0*u*u*u*p + 10.0*u*u*q + 5.0*u*r + rm5;
    u7 = rho*u*u*u*u*u*u + 15.0*u*u*u*u*p + 20.0*u*u*u*q + 15.0*u*u*r + 6.0*u*rm5 + rm6;

    W[1] = rho;
    W[2] = u;
    W[3] = p;
    if(Levermore1D_Vector::get_length() > 3) {
    W[4] = q;
    W[5] = r;
    }
    if(Levermore1D_Vector::get_length() > 5) {
    W[6] = rm5;
    W[7] = rm6;
    }

    Levermore1D_cState U(W);

    ensure_distance("u1=u1",u1,U[1],fabs(u1)*tol);
    ensure_distance("u2=u2",u2,U[2],fabs(u2)*tol);
    ensure_distance("u3=u3",u3,U[3],fabs(u3)*tol);
    if(Levermore1D_Vector::get_length() > 3) {
      ensure_distance("u4=u4",u4,U[4],fabs(u4)*tol);
      ensure_distance("u5=u5",u5,U[5],fabs(u5)*tol);
    }
    if(Levermore1D_Vector::get_length() > 5) {
      ensure_distance("u6=u6",u6,U[6],fabs(u6)*tol);
      ensure_distance("u7=u7",u7,U[7],fabs(u7)*tol);
    }
  }


  //end tests
}



// Test suite constructor
tut::Levermore1D_cState_TestSuite Levermore1D_cStateTestSuite("Class:Levermore1D_cState");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

