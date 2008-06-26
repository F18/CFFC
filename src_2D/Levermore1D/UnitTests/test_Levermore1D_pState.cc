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
      b = sqrt((double)i) * exp((double)i);
      W1[i] = a;
      W2[i] = b;
    }

    W3 = W1 + W2;

    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp((double)i);
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
      b = sqrt((double)i) * exp((double)i);
      W1[i] = a;
      W2[i] = b;
    }

    W3 = (W2 += W1);

    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp((double)i);

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
      b = sqrt((double)i) * exp((double)i);
      W1[i] = a;
      W2[i] = b;
    }

    W3 = W1 - W2;

    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp((double)i);
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
      b = sqrt((double)i) * exp((double)i);
      W1[i] = a;
      W2[i] = b;
    }

    W3 = (W2 -= W1);

    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp((double)i);

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
      b = sqrt((double)i) * exp((double)i);
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
      b = sqrt((double)i) * exp((double)i);
      W1[i] = a;
      W2[i] = b;
    }

    W3 = W1 ^ W2;
    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp((double)i);
      ensure_distance("W3 = W1^W2)", W3[i], a*b, fabs(a*b)*tol);
    }

    W3 = W2 ^ W1;
    for(i=1;i<=Levermore1D_pState::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp((double)i);
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

  /* Test 14:*/
  template<>
  template<>
  void Levermore1D_pState_object::test<14>()
  {
    set_test_name("Maxwell-Boltzmann");

    double rho(2.1), u(-120.2), p(113998.0), expected(0.0);

    Levermore1D_pState W(rho,u,p);

    ensure_distance("rho==rho",W[1],rho,fabs(rho)*tol);
    ensure_distance("u==u",W[2],u,fabs(u)*tol);
    ensure_distance("p==p",W[3],p,fabs(p)*tol);
    for(int i=4; i <= Levermore1D_Vector::get_length(); i=i+2) {
      ensure_distance("Odd moment == zero",W[i],0.0,tol);
      expected = ((double)Double_Factorial(i-1)) * pow(p,(i+1)/2 ) / pow(rho,(i+1)/2-1);
      ensure_distance("Even moment == expected",W[i+1],expected,fabs(expected)*tol);
    }

  }

  /* Test 15:*/
  template<>
  template<>
  void Levermore1D_pState_object::test<15>()
  {
    set_test_name("dUdW (cState version as well)");

    double rho(1.005), u(-210.1), p(75663.0), q(-765.4e4), r(1.2345e8);

    Levermore1D_pState W(rho,u,p);
    if(Levermore1D_Vector::get_length() > 3) {
      W[4] = q;
      W[5] = r;
    }

    Levermore1D_cState U(W);

    DenseMatrix dUdW1, dUdW1a, dUdW2(Levermore1D_Vector::get_length(),
				     Levermore1D_Vector::get_length());

    dUdW1  = W.dUdW();
    dUdW1a = U.dUdW();

    dUdW2.zero();

    dUdW2(0,0) = 1.0;
    dUdW2(1,0) = u;
    dUdW2(2,0) = u*u;

    dUdW2(1,1) = rho;
    dUdW2(2,1) = 2.0*rho*u;

    dUdW2(2,2) = 1.0;

    if(Levermore1D_Vector::get_length() > 3) {
      dUdW2(3,0) = u*u*u;
      dUdW2(4,0) = u*u*u*u;

      dUdW2(3,1) = 3.0*(rho*u*u + p);
      dUdW2(4,1) = 4.0*rho*u*u*u + 12.0*u*p + 4.0*q;

      dUdW2(3,2) = 3.0*u;
      dUdW2(4,2) = 6.0*u*u;

      dUdW2(3,3) = 1.0;
      dUdW2(4,3) = 4.0*u;

      dUdW2(4,4) = 1.0;
    }

    if(Levermore1D_Vector::get_length() > 5) {
      dUdW2(5,0) = u*u*u*u*u;
      dUdW2(6,0) = u*u*u*u*u*u;

      dUdW2(5,1) = 5.0*rho*u*u*u*u + 30.0*u*u*p + 20.0*u*q + 5.0*r;
      dUdW2(6,1) = 6.0*rho*u*u*u*u*u + 60.0*u*u*u*p + 60.0*u*u*q + 30.0*u*r + 6.0*W[6];

      dUdW2(5,2) = 10.0*u*u*u;
      dUdW2(6,2) = 15.0*u*u*u*u;

      dUdW2(5,3) = 10.0*u*u;
      dUdW2(6,3) = 20.0*u*u*u;

      dUdW2(5,4) = 5.0*u;
      dUdW2(6,4) = 15.0*u*u;

      dUdW2(5,5) = 1.0;
      dUdW2(6,5) = 6.0*u;

      dUdW2(6,6) = 1.0;
    }

    if(Levermore1D_Vector::get_length() > 7) {
      dUdW2(7,0) = u*u*u*u*u*u*u;
      dUdW2(8,0) = u*u*u*u*u*u*u*u;

      dUdW2(7,1) = 7.0*rho*u*u*u*u*u*u + 105.0*u*u*u*u*p + 140.0*u*u*u*q +
                   105*u*u*r + 42.0*u*W[6] + 7.0*W[7];
      dUdW2(8,1) = 8.0*rho*u*u*u*u*u*u*u + 168.0*u*u*u*u*u*p + 280.0*u*u*u*u*q +
	           280.0*u*u*u*r + 168.0*u*u*W[6] + 56.0*u*W[7] + 8.0*W[8];

      dUdW2(7,2) = 21.0*u*u*u*u*u;
      dUdW2(8,2) = 28.0*u*u*u*u*u*u;

      dUdW2(7,3) = 35.0*u*u*u*u;
      dUdW2(8,3) = 56.0*u*u*u*u*u;

      dUdW2(7,4) = 35.0*u*u*u;
      dUdW2(8,4) = 70.0*u*u*u*u;

      dUdW2(7,5) = 21.0*u*u;
      dUdW2(8,5) = 56.0*u*u*u;

      dUdW2(7,6) = 7.0*u;
      dUdW2(8,6) = 28.0*u*u;

      dUdW2(7,7) = 1.0;
      dUdW2(8,7) = 8.0*u;

      dUdW2(8,8) = 1.0;
    }

    for(int i=0; i < min(9,Levermore1D_Vector::get_length()); ++i) {
      for(int j=0; j < min(9,Levermore1D_Vector::get_length()); ++j) {
	ensure_distance("dUdW1==dUdW2",dUdW1(i,j),dUdW2(i,j),fabs(dUdW1(i,j)*tol+tol));
	ensure_distance("dUdW1a==dUdW2",dUdW1a(i,j),dUdW2(i,j),fabs(dUdW1a(i,j)*tol+tol));
      }
    }
  }

  /* Test 16:*/
  template<>
  template<>
  void Levermore1D_pState_object::test<16>()
  {
    set_test_name("dU_MBdW (cState version as well)");

    double rho(1.1), u(-102.1), p(85222.0), q(-589.4e4), r(1.5432e8);

    Levermore1D_pState W(rho,u,p);
    if(Levermore1D_Vector::get_length() > 3) {
      W[4] = q;
      W[5] = r;
    }

    Levermore1D_cState U(W);

    DenseMatrix dU_MBdW1, dU_MBdW1a, dU_MBdW2(Levermore1D_Vector::get_length(),
					      Levermore1D_Vector::get_length());

    dU_MBdW1  = W.dU_MBdW();
    dU_MBdW1a = U.dU_MBdW();

    dU_MBdW2.zero();

    dU_MBdW2(0,0) = 1.0;
    dU_MBdW2(1,0) = u;
    dU_MBdW2(2,0) = u*u;

    dU_MBdW2(1,1) = rho;
    dU_MBdW2(2,1) = 2.0*rho*u;

    dU_MBdW2(2,2) = 1.0;

    if(Levermore1D_Vector::get_length() > 3) {
      dU_MBdW2(3,0) = u*u*u;
      dU_MBdW2(4,0) = u*u*u*u-3*p*p/rho/rho;

      dU_MBdW2(3,1) = 3.0*(rho*u*u + p);
      dU_MBdW2(4,1) = 4.0*rho*u*u*u + 12.0*u*p;

      dU_MBdW2(3,2) = 3.0*u;
      dU_MBdW2(4,2) = 6.0*u*u+6.0*p/rho;
    }

    if(Levermore1D_Vector::get_length() > 5) {
      dU_MBdW2(5,0) = u*u*u*u*u - 15.0*u*p*p/rho/rho;
      dU_MBdW2(6,0) = u*u*u*u*u*u - 30.0*p*p*p/rho/rho/rho - 45.0*u*u*p*p/rho/rho;

      dU_MBdW2(5,1) = 5.0*rho*u*u*u*u + 30.0*u*u*p + 15.0*p*p/rho;
      dU_MBdW2(6,1) = 6.0*rho*u*u*u*u*u + 60.0*u*u*u*p + 90.0*u*p*p/rho;

      dU_MBdW2(5,2) = 10.0*u*u*u + 30.0*u*p/rho;
      dU_MBdW2(6,2) = 15.0*u*u*u*u + 90.0*u*u*p/rho + 45.0*p*p/rho/rho;
    }

    if(Levermore1D_Vector::get_length() > 7) {
      dU_MBdW2(7,0) = u*u*u*u*u*u*u - 105.0*u*u*u*p*p/rho/rho - 210.0*u*p*p*p/rho/rho/rho;
      dU_MBdW2(8,0) = u*u*u*u*u*u*u*u - 210.0*u*u*u*u*p*p/rho/rho
	              - 840.0*u*u*p*p*p/rho/rho/rho - 315.0*p*p*p*p/rho/rho/rho/rho;

      dU_MBdW2(7,1) = 7.0*rho*u*u*u*u*u*u + 105*u*u*u*u*p + 315*u*u*p*p/rho + 105.0*p*p*p/rho/rho;
      dU_MBdW2(8,1) = 8.0*rho*u*u*u*u*u*u*u + 168.0*u*u*u*u*u*p + 840.0*u*u*u*p*p/rho + 840.0*u*p*p*p/rho/rho;

      dU_MBdW2(7,2) = 21.0*u*u*u*u*u + 210.0*u*u*u*p/rho + 315.0*u*p*p/rho/rho;
      dU_MBdW2(8,2) = 28.0*u*u*u*u*u*u + 420.0*u*u*u*u*p/rho + 1260.0*u*u*p*p/rho/rho + 420.0*p*p*p/rho/rho/rho;
    }

    for(int i=0; i < min(9,Levermore1D_Vector::get_length()); ++i) {
      for(int j=0; j < min(9,Levermore1D_Vector::get_length()); ++j) {
	ensure_distance("dU_MBdW1==dU_MBdW2",dU_MBdW1(i,j),dU_MBdW2(i,j),fabs(dU_MBdW1(i,j)*tol+tol));
	ensure_distance("dU_MBdW1a==dU_MBdW2",dU_MBdW1a(i,j),dU_MBdW2(i,j),fabs(dU_MBdW1a(i,j)*tol+tol));
      }
    }
  }

  /* Test 17:*/
  template<>
  template<>
  void Levermore1D_pState_object::test<17>()
  {
    set_test_name("dU_MBdW (cState version as well)");

    double rho(1.1), u(0.0), p(85222.0), q(0.0), r(1.5432e8);

    Levermore1D_pState W(rho,u,p);
    if(Levermore1D_Vector::get_length() > 3) {
      W[4] = q;
      W[5] = r;
    }

    Levermore1D_cState U(W);

    DenseMatrix dU_MBdW1, dU_MBdW1a, dU_MBdW2(Levermore1D_Vector::get_length(),
					      Levermore1D_Vector::get_length());

    dU_MBdW1  = W.dU_MBdW();
    dU_MBdW1a = U.dU_MBdW();

    dU_MBdW2.zero();

    dU_MBdW2(0,0) = 1.0;
    dU_MBdW2(1,0) = u;
    dU_MBdW2(2,0) = u*u;

    dU_MBdW2(1,1) = rho;
    dU_MBdW2(2,1) = 2.0*rho*u;

    dU_MBdW2(2,2) = 1.0;

    if(Levermore1D_Vector::get_length() > 3) {
      dU_MBdW2(3,0) = u*u*u;
      dU_MBdW2(4,0) = u*u*u*u-3*p*p/rho/rho;

      dU_MBdW2(3,1) = 3.0*(rho*u*u + p);
      dU_MBdW2(4,1) = 4.0*rho*u*u*u + 12.0*u*p;

      dU_MBdW2(3,2) = 3.0*u;
      dU_MBdW2(4,2) = 6.0*u*u+6.0*p/rho;
    }

    if(Levermore1D_Vector::get_length() > 5) {
      dU_MBdW2(5,0) = u*u*u*u*u - 15.0*u*p*p/rho/rho;
      dU_MBdW2(6,0) = u*u*u*u*u*u - 30.0*p*p*p/rho/rho/rho - 45.0*u*u*p*p/rho/rho;

      dU_MBdW2(5,1) = 5.0*rho*u*u*u*u + 30.0*u*u*p + 15.0*p*p/rho;
      dU_MBdW2(6,1) = 6.0*rho*u*u*u*u*u + 60.0*u*u*u*p + 90.0*u*p*p/rho;

      dU_MBdW2(5,2) = 10.0*u*u*u + 30.0*u*p/rho;
      dU_MBdW2(6,2) = 15.0*u*u*u*u + 90.0*u*u*p/rho + 45.0*p*p/rho/rho;
    }

    if(Levermore1D_Vector::get_length() > 7) {
      dU_MBdW2(7,0) = u*u*u*u*u*u*u - 105.0*u*u*u*p*p/rho/rho - 210.0*u*p*p*p/rho/rho/rho;
      dU_MBdW2(8,0) = u*u*u*u*u*u*u*u - 210.0*u*u*u*u*p*p/rho/rho
	              - 840.0*u*u*p*p*p/rho/rho/rho - 315.0*p*p*p*p/rho/rho/rho/rho;

      dU_MBdW2(7,1) = 7.0*rho*u*u*u*u*u*u + 105*u*u*u*u*p + 315*u*u*p*p/rho + 105.0*p*p*p/rho/rho;
      dU_MBdW2(8,1) = 8.0*rho*u*u*u*u*u*u*u + 168.0*u*u*u*u*u*p + 840.0*u*u*u*p*p/rho + 840.0*u*p*p*p/rho/rho;

      dU_MBdW2(7,2) = 21.0*u*u*u*u*u + 210.0*u*u*u*p/rho + 315.0*u*p*p/rho/rho;
      dU_MBdW2(8,2) = 28.0*u*u*u*u*u*u + 420.0*u*u*u*u*p/rho + 1260.0*u*u*p*p/rho/rho + 420.0*p*p*p/rho/rho/rho;
    }

    for(int i=0; i < min(9,Levermore1D_Vector::get_length()); ++i) {
      for(int j=0; j < min(9,Levermore1D_Vector::get_length()); ++j) {
	ensure_distance("dU_MBdW1==dU_MBdW2",dU_MBdW1(i,j),dU_MBdW2(i,j),fabs(dU_MBdW1(i,j)*tol+tol));
	ensure_distance("dU_MBdW1a==dU_MBdW2",dU_MBdW1a(i,j),dU_MBdW2(i,j),fabs(dU_MBdW1a(i,j)*tol+tol));
      }
    }
  }

  /* Test 18:*/
  template<>
  template<>
  void Levermore1D_pState_object::test<18>()
  {
    set_test_name("dSdW");

    double rho(1.1), u(-102.1), p(85222.0), q(-589.4e4), r(1.5432e8), tau,
      M5(2.2113e10), M6(1.112e12), M7(9.876e13), M8(5.434e16);

    Levermore1D_pState W(rho,u,p);

    if(Levermore1D_Vector::get_length() > 3) {
      W[4] = q;
      W[5] = r;
    }

    if(Levermore1D_Vector::get_length() > 5) {
      W[6] = M5;
      W[7] = M6;
    }

    if(Levermore1D_Vector::get_length() > 7) {
      W[8] = M7;
      W[9] = M8;
    }

    tau = W.relaxation_time();

    DenseMatrix dSdW1, dSdW2(max(9,Levermore1D_Vector::get_length()),
			     max(9,Levermore1D_Vector::get_length()));

    dSdW1  = W.dSdW();

    dSdW2.zero();

    // from Maple:
    dSdW2(0,0) = 0.0;
    dSdW2(0,1) = 0.0;
    dSdW2(0,2) = 0.0;
    dSdW2(0,3) = 0.0;
    dSdW2(0,4) = 0.0;
    dSdW2(0,5) = 0.0;
    dSdW2(0,6) = 0.0;
    dSdW2(0,7) = 0.0;
    dSdW2(0,8) = 0.0;
    dSdW2(1,0) = 0.0;
    dSdW2(1,1) = 0.0;
    dSdW2(1,2) = 0.0;
    dSdW2(1,3) = 0.0;
    dSdW2(1,4) = 0.0;
    dSdW2(1,5) = 0.0;
    dSdW2(1,6) = 0.0;
    dSdW2(1,7) = 0.0;
    dSdW2(1,8) = 0.0;
    dSdW2(2,0) = 0.0;
    dSdW2(2,1) = 0.0;
    dSdW2(2,2) = 0.0;
    dSdW2(2,3) = 0.0;
    dSdW2(2,4) = 0.0;
    dSdW2(2,5) = 0.0;
    dSdW2(2,6) = 0.0;
    dSdW2(2,7) = 0.0;
    dSdW2(2,8) = 0.0;
    dSdW2(3,0) = 0.0;
    dSdW2(3,1) = 0.0;
    dSdW2(3,2) = 0.0;
    dSdW2(3,3) = -1/tau;
    dSdW2(3,4) = 0.0;
    dSdW2(3,5) = 0.0;
    dSdW2(3,6) = 0.0;
    dSdW2(3,7) = 0.0;
    dSdW2(3,8) = 0.0;
    dSdW2(4,0) = -3.0/tau*p*p/(rho*rho);
    dSdW2(4,1) = -4.0/tau*q;
    dSdW2(4,2) = 6.0/tau*p/rho;
    dSdW2(4,3) = -4.0/tau*u;
    dSdW2(4,4) = -1/tau;
    dSdW2(4,5) = 0.0;
    dSdW2(4,6) = 0.0;
    dSdW2(4,7) = 0.0;
    dSdW2(4,8) = 0.0;
    dSdW2(5,0) = -15.0/tau*u*p*p/(rho*rho);
    dSdW2(5,1) = 1/tau*(15.0*p*p/rho-20.0*u*q-5.0*r);
    dSdW2(5,2) = 30.0/tau*u*p/rho;
    dSdW2(5,3) = -10.0/tau*u*u;
    dSdW2(5,4) = -5.0/tau*u;
    dSdW2(5,5) = -1/tau;
    dSdW2(5,6) = 0.0;
    dSdW2(5,7) = 0.0;
    dSdW2(5,8) = 0.0;
    dSdW2(6,0) = 1/tau*(-30.0*p*p*p/(rho*rho*rho)-45.0*u*u*p*p/(rho*rho));
    dSdW2(6,1) = 1/tau*(90.0*u*p*p/rho-60.0*u*u*q-30.0*u*r-6.0*M5);
    dSdW2(6,2) = 1/tau*(45.0*p*p/(rho*rho)+90.0*u*u*p/rho);
    dSdW2(6,3) = -20.0/tau*u*u*u;
    dSdW2(6,4) = -15.0/tau*u*u;
    dSdW2(6,5) = -6.0/tau*u;
    dSdW2(6,6) = -1/tau;
    dSdW2(6,7) = 0.0;
    dSdW2(6,8) = 0.0;
    dSdW2(7,0) = 1/tau*(-210.0*u*p*p*p/(rho*rho*rho)-105.0*u*u*u*p*p/(rho*rho));
    dSdW2(7,1) = 1/tau*(105.0*p*p*p/(rho*rho)+315.0*u*u*p*p/rho-140.0*u*u*u*q-105.0*u*u*r-42.0*u*M5-7.0*M6);
    dSdW2(7,2) = 1/tau*(315.0*u*p*p/(rho*rho)+210.0*u*u*u*p/rho);
    dSdW2(7,3) = -35.0/tau*u*u*u*u;
    dSdW2(7,4) = -35.0/tau*u*u*u;
    dSdW2(7,5) = -21.0/tau*u*u;
    dSdW2(7,6) = -7.0/tau*u;
    dSdW2(7,7) = -1/tau;
    dSdW2(7,8) = 0.0;
    dSdW2(8,0) = 1/tau*(-840.0*u*u*p*p*p/(rho*rho*rho)-210.0*u*u*u*u*p*p/(rho*rho)-315.0*p*p*p*p/(rho*rho*rho*rho));
    dSdW2(8,1) = 1/tau*(840.0*u*p*p*p/(rho*rho)+840.0*u*u*u*p*p/rho-280.0*u*u*u*u*q-280.0*u*u*u*r-168.0*u*u*M5-56.0*u*M6-8.0*M7);
    dSdW2(8,2) = 1/tau*(1260.0*u*u*p*p/(rho*rho)+420.0*u*u*u*u*p/rho+420.0*p*p*p/(rho*rho*rho));
    dSdW2(8,3) = -56.0/tau*u*u*u*u*u;
    dSdW2(8,4) = -70.0/tau*u*u*u*u;
    dSdW2(8,5) = -56.0/tau*u*u*u;
    dSdW2(8,6) = -28.0/tau*u*u;
    dSdW2(8,7) = -8.0/tau*u;
    dSdW2(8,8) = -1/tau;

    for(int i=0; i < min(9,Levermore1D_Vector::get_length()); ++i) {
      for(int j=0; j < min(9,Levermore1D_Vector::get_length()); ++j) {
	ensure_distance("dSdW1==dSdW2",dSdW1(i,j),dSdW2(i,j),fabs(dSdW1(i,j)*tol+tol));
      }
    }
  }

  /* Test 19:*/
  template<>
  template<>
  void Levermore1D_pState_object::test<19>()
  {
    set_test_name("dSdU");

    double rho(2.1), u(-987.1), p(101222.0), q(-589.4e4), r(1.432e8), tau,
      M5(8.5113e10), M6(1.882e12), M7(-9.7611e13), M8(8.434e16);

    Levermore1D_pState W(rho,u,p);

    if(Levermore1D_Vector::get_length() > 3) {
      W[4] = q;
      W[5] = r;
    }

    if(Levermore1D_Vector::get_length() > 5) {
      W[6] = M5;
      W[7] = M6;
    }

    if(Levermore1D_Vector::get_length() > 7) {
      W[8] = M7;
      W[9] = M8;
    }

    tau = W.relaxation_time();

    DenseMatrix dSdU1, dSdU2(max(9,Levermore1D_Vector::get_length()),
			     max(9,Levermore1D_Vector::get_length()));

    dSdU1  = W.dSdU();

    dSdU2.zero();

    // from Maple:
    dSdU2(0,0) = 0.0;
    dSdU2(0,1) = 0.0;
    dSdU2(0,2) = 0.0;
    dSdU2(0,3) = 0.0;
    dSdU2(0,4) = 0.0;
    dSdU2(0,5) = 0.0;
    dSdU2(0,6) = 0.0;
    dSdU2(0,7) = 0.0;
    dSdU2(0,8) = 0.0;
    dSdU2(1,0) = 0.0;
    dSdU2(1,1) = 0.0;
    dSdU2(1,2) = 0.0;
    dSdU2(1,3) = 0.0;
    dSdU2(1,4) = 0.0;
    dSdU2(1,5) = 0.0;
    dSdU2(1,6) = 0.0;
    dSdU2(1,7) = 0.0;
    dSdU2(1,8) = 0.0;
    dSdU2(2,0) = 0.0;
    dSdU2(2,1) = 0.0;
    dSdU2(2,2) = 0.0;
    dSdU2(2,3) = 0.0;
    dSdU2(2,4) = 0.0;
    dSdU2(2,5) = 0.0;
    dSdU2(2,6) = 0.0;
    dSdU2(2,7) = 0.0;
    dSdU2(2,8) = 0.0;
    dSdU2(3,0) = 1/tau*u*(rho*u*u-3.0*p)/rho;
    dSdU2(3,1) = -3.0/tau*(rho*u*u-p)/rho;
    dSdU2(3,2) = 3.0/tau*u;
    dSdU2(3,3) = -1/tau;
    dSdU2(3,4) = 0.0;
    dSdU2(3,5) = 0.0;
    dSdU2(3,6) = 0.0;
    dSdU2(3,7) = 0.0;
    dSdU2(3,8) = 0.0;
    dSdU2(4,0) = -3.0/tau*p*p/(rho*rho)+4.0/tau*q*u/rho+6.0/tau*p/rho*u*u+
      4.0/tau*u*u*(rho*u*u-3.0*p)/rho-1/tau*u*(rho*u*u*u+4.0*q)/rho;
    dSdU2(4,1) = -4.0/tau*q/rho-12.0/tau*u*p/rho-12.0/tau*u*(rho*u*u-p)/rho
      +4.0/tau*(rho*u*u*u+q)/rho;
    dSdU2(4,2) = 6.0/tau*p/rho+6.0/tau*u*u;
    dSdU2(4,3) = 0.0;
    dSdU2(4,4) = -1/tau;
    dSdU2(4,5) = 0.0;
    dSdU2(4,6) = 0.0;
    dSdU2(4,7) = 0.0;
    dSdU2(4,8) = 0.0;
    dSdU2(5,0) = -15.0/tau*u*p*p/(rho*rho)-1/tau*
      (15.0*p*p/rho-20.0*u*q-5.0*r)*u/rho+30.0/tau*u*u*u*p/rho+10.0/tau*u*u*u*
      (rho*u*u-3.0*p)/rho-5.0/tau*u*u*(rho*u*u*u+4.0*q)/rho-1/tau*u*(-rho*u*u*u*u+5.0*r)/rho;
    dSdU2(5,1) = 1/tau*(15.0*p*p/rho-20.0*u*q-5.0*r)/rho-60.0/tau*p/rho*u*u
      -30.0/tau*u*u*(rho*u*u-p)/rho+20.0/tau*u*(rho*u*u*u+q)/rho+5.0/tau*(-rho*u*u*u*u+r)/rho;
    dSdU2(5,2) = 30.0/tau*u*p/rho+10.0/tau*u*u*u;
    dSdU2(5,3) = 0.0;
    dSdU2(5,4) = 0.0;
    dSdU2(5,5) = -1/tau;
    dSdU2(5,6) = 0.0;
    dSdU2(5,7) = 0.0;
    dSdU2(5,8) = 0.0;
    dSdU2(6,0) = 1/tau*(-30.0*p*p*p/(rho*rho*rho)-45.0*u*u*p*p/(rho*rho))-
      1/tau*(90.0*u*p*p/rho-60.0*u*u*q-30.0*u*r-6.0*M5)*u/rho+1/tau*
      (45.0*p*p/(rho*rho)+90.0*u*u*p/rho)*u*u+20.0/tau*u*u*u*u*(rho*u*u-3.0*p)/rho-15.0/tau*u*u*u*
      (rho*u*u*u+4.0*q)/rho-6.0/tau*u*u*(-rho*u*u*u*u+5.0*r)/rho-1/tau*u*(rho*u*u*u*u*u+6.0*M5)/rho;
    dSdU2(6,1) = 1/tau*(90.0*u*p*p/rho-60.0*u*u*q-30.0*u*r-6.0*M5)/rho-2.0/
      tau*(45.0*p*p/(rho*rho)+90.0*u*u*p/rho)*u-60.0/tau*u*u*u*(rho*u*u-p)/rho+60.0/
      tau*u*u*(rho*u*u*u+q)/rho+30.0/tau*u*(-rho*u*u*u*u+r)/rho+6.0/tau*(rho*u*u*u*u*
									 u+M5)/rho;
    dSdU2(6,2) = 1/tau*(45.0*p*p/(rho*rho)+90.0*u*u*p/rho)+15.0/tau*u*u*u*u;
    dSdU2(6,3) = 0.0;
    dSdU2(6,4) = 0.0;
    dSdU2(6,5) = 0.0;
    dSdU2(6,6) = -1/tau;
    dSdU2(6,7) = 0.0;
    dSdU2(6,8) = 0.0;
    dSdU2(7,0) = 1/tau*(-210.0*u*p*p*p/(rho*rho*rho)-105.0*u*u*u*p*p/
			(rho*rho))-1/tau*(105.0*p*p*p/(rho*rho)+315.0*u*u*p*p/rho-
					  140.0*u*u*u*q-105.0*u*u*r-42.0*u*M5-7.0*M6)*u/rho+
      1/tau*(315.0*u*p*p/(rho*rho)+210.0*u*u*u*p/rho)*u*u+
      35.0/tau*u*u*u*u*u*(rho*u*u-3.0*p)/rho-35.0/tau*u*u*u*u*(rho*u*u*u+4.0*q)/rho
      -21.0/tau*u*u*u*(-rho*u*u*u*u+5.0*r)/rho-7.0/tau*u*u*(rho*u*u*u*u*u+6.0*M5)/rho
      -1/tau*u*(-rho*u*u*u*u*u*u+7.0*M6)/rho;
    dSdU2(7,1) = 1/tau*(105.0*p*p*p/(rho*rho)+315.0*u*u*p*p/rho-140.0*u*u*u*q-
			105.0*u*u*r-42.0*u*M5-7.0*M6)/rho-2.0/tau*
      (315.0*u*p*p/(rho*rho)+210.0*u*u*u*p/rho)*u-105.0/tau*u*u*u*u*(rho*u*u-p)/rho+140.0/tau*u*u*u*(rho*u*u*u+q)/rho+
      105.0/tau*u*u*(-rho*u*u*u*u+r)/rho+42.0/tau*u*(rho*u*u*u*u*u+M5)/rho+7.0/tau*(-rho*u*u*u*u*u*u+M6)/rho;
    dSdU2(7,2) = 1/tau*(315.0*u*p*p/(rho*rho)+210.0*u*u*u*p/rho)+21.0/tau*u*u*u*u*u;
    dSdU2(7,3) = 0.0;
    dSdU2(7,4) = 0.0;
    dSdU2(7,5) = 0.0;
    dSdU2(7,6) = 0.0;
    dSdU2(7,7) = -1/tau;
    dSdU2(7,8) = 0.0;
    dSdU2(8,0) = 1/tau*(-840.0*u*u*p*p*p/(rho*rho*rho)-210.0*u*u*u*u*p*p/
			(rho*rho)-315.0*p*p*p*p/(rho*rho*rho*rho))-
      1/tau*(840.0*u*p*p*p/(rho*rho)+840.0*u*u*u*p*p/rho-280.0*u*u*u*u*q-280.0*u*u*u*r-
	     168.0*u*u*M5-56.0*u*M6-8.0*M7)*u/rho+1/tau*
      (1260.0*u*u*p*p/(rho*rho)+420.0*u*u*u*u*p/rho+420.0*p*p*p/(rho*rho*rho))*u*u+
      56.0/tau*u*u*u*u*u*u*(rho*u*u-3.0*p)/rho-70.0/tau*u*u*u*u*u*(rho*u*u*u+4.0*q)/rho
      -56.0/tau*u*u*u*u*(-rho*u*u*u*u+5.0*r)/rho-28.0/tau*u*u*u*(rho*u*u*u*u*u+6.0*M5)/rho
      -8.0/tau*u*u*(-rho*u*u*u*u*u*u+7.0*M6)/rho-1/tau*u*(rho*u*u*u*u*u*u*u+8.0*M7)/rho;
    dSdU2(8,1) = 1/tau*(840.0*u*p*p*p/(rho*rho)+840.0*u*u*u*p*p/rho-280.0*u*u*u*u*q
			-280.0*u*u*u*r-168.0*u*u*M5-56.0*u*M6-8.0*M7)/rho
      -2.0/tau*(1260.0*u*u*p*p/(rho*rho)+420.0*u*u*u*u*p/rho+420.0*p*p*p/(rho*rho*rho))*u
      -168.0/tau*u*u*u*u*u*(rho*u*u-p)/rho+280.0/tau*u*u*u*u*(rho*u*u*u+q)/rho
      +280.0/tau*u*u*u*(-rho*u*u*u*u+r)/rho+168.0/tau*u*u*(rho*u*u*u*u*u+M5)/rho
      +56.0/tau*u*(-rho*u*u*u*u*u*u+M6)/rho+8.0/tau*(rho*u*u*u*u*u*u*u+M7)/rho;
    dSdU2(8,2) = 1/tau*(1260.0*u*u*p*p/(rho*rho)+420.0*u*u*u*u*p/rho+420.0*
			p*p*p/(rho*rho*rho))+28.0/tau*u*u*u*u*u*u;
    dSdU2(8,3) = 0.0;
    dSdU2(8,4) = 0.0;
    dSdU2(8,5) = 0.0;
    dSdU2(8,6) = 0.0;
    dSdU2(8,7) = 0.0;
    dSdU2(8,8) = -1/tau;

    for(int i=0; i < min(9,Levermore1D_Vector::get_length()); ++i) {
      for(int j=0; j < min(9,Levermore1D_Vector::get_length()); ++j) {
	ensure_distance("dSdU1==dSdU2",dSdU1(i,j),dSdU2(i,j),fabs(dSdU1(i,0)*tol*HUNDRED+tol*HUNDRED));
      }
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

