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

  /* Test 3:*/
  template<>
  template<>
  void Levermore1D_cState_object::test<3>()
  {
    set_test_name("Constructor from weights");

    double rho(1.225);
    double u(132.22);
    double p(101325.0);

    double momentum(rho*u);
    double e(p+rho*u*u);

    double m = Levermore1D_weights::m();
    double n(rho/m); //number density

    double B(rho/(2.0*p));

    Levermore1D_weights A;
    A.zero();

    A[1] = -B*u*u+log(n*sqrt(B/PI));
    A[2] = 2.0*B*u;
    A[3] = -B;

    Levermore1D_cState U(A);

    ensure_distance("density is equal", rho, U[1], fabs(rho)*1e-12+1e-12);
    ensure_distance("momentum is equal", momentum, U[2], fabs(momentum)*1e-12+1e-12);
    ensure_distance("energy is equal", e, U[3], fabs(e)*1e-12+1e-12);

  }

  /* Test 4:*/
  template<>
  template<>
  void Levermore1D_cState_object::test<4>()
  {
    set_test_name("Assignment Operator");

    int i(0);
    Levermore1D_cState U1, U2, U3;

    for(i=1;i<=Levermore1D_cState::get_length();++i) {
      U1[i] = pow((double)i,1.23456) / 98.765;
    }

    U3 = U2 = U1;

    for(i=1;i<=Levermore1D_cState::get_length();++i) {
      ensure_distance("U2=U1",U2[i],U1[i],fabs(U1[i])*tol);
      ensure_distance("U3=U1",U3[i],U1[i],fabs(U1[i])*tol);
    }
  }

  /* Test 5:*/
  template<>
  template<>
  void Levermore1D_cState_object::test<5>()
  {
    set_test_name("Addition Operator");

    int i(0);
    double a(0.0), b(0.0);
    Levermore1D_cState U1, U2, U3;

    for(i=1;i<=Levermore1D_cState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      U1[i] = a;
      U2[i] = b;
    }

    U3 = U1 + U2;

    for(i=1;i<=Levermore1D_cState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      ensure_distance("U3 = U1 + U2", U3[i], a+b, fabs(a+b)*tol);
    }
  }

  /* Test 6:*/
  template<>
  template<>
  void Levermore1D_cState_object::test<6>()
  {
    set_test_name("Addition and Assignment Operator");

    int i(0);
    double a(0.0), b(0.0);
    Levermore1D_cState U1, U2, U3;

    for(i=1;i<=Levermore1D_cState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      U1[i] = a;
      U2[i] = b;
    }

    U3 = (U2 += U1);

    for(i=1;i<=Levermore1D_cState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);

      ensure_distance("U2 += U1",U2[i],a+b,fabs(a+b)*tol);
      ensure_distance("U3 = (U2+=U1)",U3[i],a+b,fabs(a+b)*tol);
    }
  }

  /* Test 7:*/
  template<>
  template<>
  void Levermore1D_cState_object::test<7>()
  {
    set_test_name("Subtraction Operator");

    int i(0);
    double a(0.0), b(0.0);
    Levermore1D_cState U1, U2, U3;

    for(i=1;i<=Levermore1D_cState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      U1[i] = a;
      U2[i] = b;
    }

    U3 = U1 - U2;

    for(i=1;i<=Levermore1D_cState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      ensure_distance("U3 = U1 - U2", U3[i], a-b, fabs(a-b)*tol);
    }
  }

  /* Test 8:*/
  template<>
  template<>
  void Levermore1D_cState_object::test<8>()
  {
    set_test_name("Subtraction and Assignment Operator");

    int i(0);
    double a(0.0), b(0.0);
    Levermore1D_cState U1, U2, U3;

    for(i=1;i<=Levermore1D_cState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      U1[i] = a;
      U2[i] = b;
    }

    U3 = (U2 -= U1);

    for(i=1;i<=Levermore1D_cState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);

      ensure_distance("U2 -= U1",U2[i],b-a,fabs(b-a)*tol);
      ensure_distance("U3 = (U2-=U1)",U3[i],b-a,fabs(b-a)*tol);
    }
  }

  /* Test 9:*/
  template<>
  template<>
  void Levermore1D_cState_object::test<9>()
  {
    set_test_name("Dot Product Operator");

    int i(0);
    double a(0.0), b(0.0), dot(0.0), dot_levermore(0.0);
    Levermore1D_cState U1, U2;

    for(i=1;i<=Levermore1D_cState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      dot += (a*b);
      U1[i] = a;
      U2[i] = b;
    }

    dot_levermore = U1 * U2;
    ensure_distance("dot_prod = U1*U2)", dot_levermore, dot, fabs(dot)*tol);
    dot_levermore = U2 * U1;
    ensure_distance("dot_prod = U2*U1)", dot_levermore, dot, fabs(dot)*tol);

  }

  /* Test 10:*/
  template<>
  template<>
  void Levermore1D_cState_object::test<10>()
  {
    set_test_name("Term-wise Product Operator");

    int i(0);
    double a(0.0), b(0.0);
    Levermore1D_cState U1, U2, U3;

    for(i=1;i<=Levermore1D_cState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      U1[i] = a;
      U2[i] = b;
    }

    U3 = U1 ^ U2;
    for(i=1;i<=Levermore1D_cState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      ensure_distance("U3 = U1^U2)", U3[i], a*b, fabs(a*b)*tol);
    }

    U3 = U2 ^ U1;
    for(i=1;i<=Levermore1D_cState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp (i);
      ensure_distance("U3 = U1^U2)", U3[i], a*b, fabs(a*b)*tol);
    }
  }

  /* Test 11:*/
  template<>
  template<>
  void Levermore1D_cState_object::test<11>()
  {
    set_test_name("copy_form Function");
    int i(0);
    Levermore1D_cState U1, U2;

    for(i=1;i<=Levermore1D_cState::get_length();++i) {
      U1[i] = pow((double)i,1.23456) / 98.765;
    }

    U2.copy_from(U1);

    for(i=1;i<=Levermore1D_cState::get_length();++i) {
      ensure_distance("U2=U1",U2[i],U1[i],fabs(U1[i])*tol);
    }
  }

  /* Test 12:*/
  template<>
  template<>
  void Levermore1D_cState_object::test<12>()
  {
    set_test_name("Vacuum, zero, one, and set_all Function");
    int i(0);
    Levermore1D_cState U1, U2, U3, U4;

    for(i=1;i<=Levermore1D_cState::get_length();++i) {
      U1[i] = pow((double)i,1.23456) / 98.765;
      U2[i] = pow((double)i,1.23456) / 98.765;
      U3[i] = pow((double)i,1.23456) / 98.765;
      U4[i] = pow((double)i,1.23456) / 98.765;
      ensure("U1!=0.0", U1[i] != 0.0);
      ensure("U2!=1.0", U2[i] != 1.0);
      ensure("U3!=2.0", U3[i] != 2.0);
      ensure("U4!=0.0", U4[i] != 0.0);
    }

    U1.zero();
    U2.one();
    U3.set_all(2.0);
    U4.Vacuum();

    for(i=1;i<=Levermore1D_cState::get_length();++i) {
      ensure_distance("U1 = 0.0",U1[i],0.0,tol);
      ensure_distance("U2 = 1.0",U2[i],1.0,tol);
      ensure_distance("U3 = 2.0",U3[i],2.0,tol);
      ensure_distance("U4 = 0.0",U4[i],0.0,tol);
    }
  }

  /* Test 13:*/
  template<>
  template<>
  void Levermore1D_cState_object::test<13>()
  {
    set_test_name("calculate moments");
    double rho(1.225);
    double u(132.22);
    double p(101325.0);

    double momentum(rho*u);
    double e(p+rho*u*u);

    double m = Levermore1D_weights::m();
    double n(rho/m); //number density

    double B(rho/(2.0*p));

    Levermore1D_weights A;
    A.zero();

    A[1] = -B*u*u+log(n*sqrt(B/PI));
    A[2] = 2.0*B*u;
    A[3] = -B;

    Levermore1D_cState U(A);
    double momentU, momentA;
    char testname[256];

    //check series of moments
    for(int i=0;i<15;++i) {
      sprintf(testname, "Integrate moment %d.", i);
      momentU = U.moment(i,A);
      momentA = A.integrate_conserved_moment(i);
      ensure_distance(testname, momentU, momentA, fabs(momentU)*1e-10+1e-10);
    }

  }

  /* Test 14:*/
  template<>
  template<>
  void Levermore1D_cState_object::test<14>()
  {
    set_test_name("calculate moments (part 2)");
    //same as before, but with zero velocity.
    double rho(1.225);
    double u(0.00);
    double p(101325.0);

    double momentum(rho*u);
    double e(p+rho*u*u);

    double m = Levermore1D_weights::m();
    double n(rho/m); //number density

    double B(rho/(2.0*p));

    Levermore1D_weights A;
    A.zero();

    A[1] = -B*u*u+log(n*sqrt(B/PI));
    A[2] = 2.0*B*u;
    A[3] = -B;

    Levermore1D_cState U(A);
    double momentU, momentA;
    char testname[256];

    //check series of moments
    for(int i=0;i<15;++i) {
      sprintf(testname, "Integrate moment %d.", i);
      momentU = U.moment(i,A);
      momentA = A.integrate_conserved_moment(i);
      ensure_distance(testname, momentU, momentA, fabs(momentU)*1e-10+1e-10);
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

