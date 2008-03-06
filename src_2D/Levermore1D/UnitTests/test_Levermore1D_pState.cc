/*!\file test_Levermore1D_pState.cc
  \brief Regression tests for template class Levermore1D_pState datatype. */

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
  class Data_Levermore1D_pState : public TestData {

    // Local variables
  public:

    // Constructor
    Data_Levermore1D_pState(){
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
  typedef test_group<Data_Levermore1D_pState> Levermore1D_pState_TestSuite;
  typedef Levermore1D_pState_TestSuite::object Levermore1D_pState_object;


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
  void Levermore1D_pState_object::test<1>()
  {
    set_test_name("Copy Constructors");

    int i(0);
    Levermore1D_pState P1;

    for(i=1;i<=Levermore1D_Vector::get_length();++i) {
      //set to a value for copy constructor test
      P1[i] = pow((double)i,6.54321) / 98.765;
    }

    Levermore1D_pState P2(P1);

    for(i=1;i<=Levermore1D_Vector::get_length();++i) {
      ensure_distance("P2=P1",P2[i],P1[i],fabs(P1[i])*tol);
    }

  }

  /* Test 2:*/
  template<>
  template<>
  void Levermore1D_pState_object::test<2>()
  {
    set_test_name("Constructor from cState");

    int i(0);
    Levermore1D_cState U;

    double rho(1.225), u(20.1), p(101325.0), q(123456.2), r(5.4321e8),
      rm5(2.333e9), rm6(7.89e12);

    U[1] = rho;
    U[2] = rho*u;
    U[3] = rho*u*u + p;
    if(Levermore1D_Vector::get_length() > 3) {
      U[4] = rho*u*u*u + 3.0*u*p + q;
      U[5] = rho*u*u*u*u + 6.0*u*u*p + 4.0*u*q + r;
    }
    if(Levermore1D_Vector::get_length() > 5) {
      U[6] = rho*u*u*u*u*u + 10.0*u*u*u*p + 10.0*u*u*q + 5.0*u*r + rm5;
      U[7] = rho*u*u*u*u*u*u + 15.0*u*u*u*u*p + 20.0*u*u*u*q + 15.0*u*u*r + 6.0*u*rm5 + rm6;
    }

    Levermore1D_pState W(U);

    ensure_distance("rho=rho",rho,W[1],fabs(rho)*tol);
    ensure_distance("u=u",u,W[2],fabs(u)*tol);
    ensure_distance("p=p",p,W[3],fabs(p)*tol);
    if(Levermore1D_Vector::get_length() > 3) {
      ensure_distance("q=q",q,W[4],fabs(q)*tol);
      ensure_distance("r=r",r,W[5],fabs(r)*tol);
    }
    if(Levermore1D_Vector::get_length() > 5) {
      ensure_distance("rm5=rm5",rm5,W[6],fabs(rm5)*tol);
      ensure_distance("rm6=rm6",rm6,W[7],fabs(rm6)*tol);
    }
  }

  /* Test 3:*/
  template<>
  template<>
  void Levermore1D_pState_object::test<3>()
  {
    set_test_name("Constructor from weights");

    double rho(1.554);
    double u(-223.1);
    double p(97322.1);

    double B(rho/(2.0*p));

    Levermore1D_weights A;
    A.zero();

    A[1] = -B*u*u+log(rho*sqrt(B/PI));
    A[2] = 2.0*B*u;
    A[3] = -B;

    Levermore1D_pState W(A,u);

    double rho2(A.integrate_conserved_moment(0,u));
    double u2(A.integrate_conserved_moment(1,u)/rho2);
    double p2(A.integrate_random_moment(2,u2,u));

    ensure_distance("density is equal", rho2, W[1], fabs(rho2)*1e-10+1e-10);
    ensure_distance("velocity is equal", u2, W[2], fabs(u2)*1e-10+1e-10);
    ensure_distance("pressure is equal", p2, W[3], fabs(p2)*1e-10+1e-10);

  }

  /* Test 4:*/
  template<>
  template<>
  void Levermore1D_pState_object::test<4>()
  {
    set_test_name("Assignment Operator");

    int i(0);
    Levermore1D_pState W1, W2, W3;

    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      W1[i] = pow((double)i,1.23456) / 98.765;
    }

    W3 = W2 = W1;

    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      ensure_distance("W2=W1",W2[i],W1[i],fabs(W1[i])*tol);
      ensure_distance("W3=W1",W3[i],W1[i],fabs(W1[i])*tol);
    }
  }

  /* Test 5:*/
  template<>
  template<>
  void Levermore1D_pState_object::test<5>()
  {
    set_test_name("Addition Operator");

    int i(0);
    double a(0.0), b(0.0);
    Levermore1D_pState W1, W2, W3;

    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      W1[i] = a;
      W2[i] = b;
    }

    W3 = W1 + W2;

    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      ensure_distance("W3 = W1 + W2", W3[i], a+b, fabs(a+b)*tol);
    }
  }

  /* Test 6:*/
  template<>
  template<>
  void Levermore1D_pState_object::test<6>()
  {
    set_test_name("Addition and Assignment Operator");

    int i(0);
    double a(0.0), b(0.0);
    Levermore1D_pState W1, W2, W3;

    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      W1[i] = a;
      W2[i] = b;
    }

    W3 = (W2 += W1);

    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);

      ensure_distance("W2 += W1",W2[i],a+b,fabs(a+b)*tol);
      ensure_distance("W3 = (W2+=W1)",W3[i],a+b,fabs(a+b)*tol);
    }
  }

  /* Test 7:*/
  template<>
  template<>
  void Levermore1D_pState_object::test<7>()
  {
    set_test_name("Subtraction Operator");

    int i(0);
    double a(0.0), b(0.0);
    Levermore1D_pState W1, W2, W3;

    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      W1[i] = a;
      W2[i] = b;
    }

    W3 = W1 - W2;

    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      ensure_distance("W3 = W1 - W2", W3[i], a-b, fabs(a-b)*tol);
    }
  }

  /* Test 8:*/
  template<>
  template<>
  void Levermore1D_pState_object::test<8>()
  {
    set_test_name("Subtraction and Assignment Operator");

    int i(0);
    double a(0.0), b(0.0);
    Levermore1D_pState W1, W2, W3;

    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      W1[i] = a;
      W2[i] = b;
    }

    W3 = (W2 -= W1);

    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);

      ensure_distance("W2 -= W1",W2[i],b-a,fabs(b-a)*tol);
      ensure_distance("W3 = (W2-=W1)",W3[i],b-a,fabs(b-a)*tol);
    }
  }

  /* Test 9:*/
  template<>
  template<>
  void Levermore1D_pState_object::test<10>()
  {
    set_test_name("Dot Product Operator");

    int i(0);
    double a(0.0), b(0.0), dot(0.0), dot_levermore(0.0);
    Levermore1D_pState W1, W2;

    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      dot += (a*b);
      W1[i] = a;
      W2[i] = b;
    }

    dot_levermore = W1 * W2;
    ensure_distance("dot_prod = W1*W2)", dot_levermore, dot, fabs(dot)*tol);
    dot_levermore = W2 * W1;
    ensure_distance("dot_prod = W2*W1)", dot_levermore, dot, fabs(dot)*tol);

  }

  /* Test 11:*/
  template<>
  template<>
  void Levermore1D_pState_object::test<11>()
  {
    set_test_name("Term-wise Product Operator");

    int i(0);
    double a(0.0), b(0.0);
    Levermore1D_pState W1, W2, W3;

    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      W1[i] = a;
      W2[i] = b;
    }

    W3 = W1 ^ W2;
    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      ensure_distance("W3 = W1^W2)", W3[i], a*b, fabs(a*b)*tol);
    }

    W3 = W2 ^ W1;
    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      ensure_distance("W3 = W1^W2)", W3[i], a*b, fabs(a*b)*tol);
    }
  }

  /* Test 12:*/
  template<>
  template<>
  void Levermore1D_pState_object::test<12>()
  {
    set_test_name("copy_form Function");
    int i(0);
    Levermore1D_pState W1, W2;

    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      W1[i] = pow((double)i,1.23456) / 98.765;
    }

    W2.copy_from(W1);

    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      ensure_distance("W2=W1",W2[i],W1[i],fabs(W1[i])*tol);
    }
  }

  /* Test 13:*/
  template<>
  template<>
  void Levermore1D_pState_object::test<13>()
  {
    set_test_name("vacuum, zero, one, and set_all Function");
    int i(0);
    Levermore1D_pState W1, W2, W3, W4;

    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      W1[i] = pow((double)i,1.23456) / 98.765;
      W2[i] = pow((double)i,1.23456) / 98.765;
      W3[i] = pow((double)i,1.23456) / 98.765;
      W4[i] = pow((double)i,1.23456) / 98.765;
      ensure("W1!=0.0", W1[i] != 0.0);
      ensure("W2!=1.0", W2[i] != 1.0);
      ensure("W3!=2.0", W3[i] != 2.0);
      ensure("W4!=2.0", W4[i] != 0.0);
    }

    W1.zero();
    W2.one();
    W3.set_all(2.0);
    W4.Vacuum();

    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      ensure_distance("W1 = 0.0",W1[i],0.0,tol);
      ensure_distance("W2 = 1.0",W2[i],1.0,tol);
      ensure_distance("W3 = 2.0",W3[i],2.0,tol);
      ensure_distance("W4 = 0.0",W4[i],0.0,tol);
    }
  }

  //end tests
}



// Test suite constructor
tut::Levermore1D_pState_TestSuite Levermore1D_pStateTestSuite("Class:Levermore1D_pState");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

