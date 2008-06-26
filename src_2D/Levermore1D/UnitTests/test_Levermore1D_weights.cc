/*!\file test_Levermore1D_weights.cc
  \brief Regression tests for template class Levermore1D_weights datatype. */

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
  class Data_Levermore1D_weights : public TestData {

    // Local variables
  public:

    // Constructor
    Data_Levermore1D_weights(){
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
  typedef test_group<Data_Levermore1D_weights> Levermore1D_weights_TestSuite;
  typedef Levermore1D_weights_TestSuite::object Levermore1D_weights_object;


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
  void Levermore1D_weights_object::test<1>()
  {
    set_test_name("Copy Constructors");

    int i(0);
    Levermore1D_weights alpha1;

    for(i=1;i<=Levermore1D_Vector::get_length();++i) {
      //set to a value for copy constructor test
      alpha1[i] = pow((double)i,3.54321) / 18.765;
    }

    Levermore1D_weights alpha2(alpha1);

    for(i=1;i<=Levermore1D_Vector::get_length();++i) {
      ensure_distance("alpha2=alpha1",alpha2[i],alpha1[i],fabs(alpha1[i])*tol);
    }

  }

  /* Test 2:*/
  template<>
  template<>
  void Levermore1D_weights_object::test<2>()
  {
    set_test_name("Constructor from pState");

    double rho(1.333);
    double u(-92.22);
    double p(90325.0);
    double m = Levermore1D_weights::m();
    double n(rho/m); //number density
    double B(rho/(2.0*p));

    Levermore1D_pState W1(rho,u,p);

    Levermore1D_weights A(W1);

    Levermore1D_pState W2(A,u);

    for(int i=1;i<=Levermore1D_Vector::get_length();i=i+2) {
      double toler = fabs(W1[i])*1e-5+1e-12;
      ensure_distance("W2=W1",W2[i],W1[i],toler);
      if(i!=Levermore1D_Vector::get_length()) {
	ensure_distance("W2=W1",W2[i+1],W1[i+1],toler);
      }
    }

  }

  /* Test 3:*/
  template<>
  template<>
  void Levermore1D_weights_object::test<3>()
  {
    set_test_name("Constructor from cState");

    double rho(1.225);
    double u(132.22);
    double p(101325.0);
    double m = Levermore1D_weights::m();
    double n(rho/m); //number density
    double B(rho/(2.0*p));

    Levermore1D_pState W1(rho,u,p);
    Levermore1D_cState U1(W1);

    Levermore1D_weights A(U1);
    Levermore1D_cState U2(A,u);

    for(int i=1;i<=Levermore1D_Vector::get_length();i=i+2) {
      double toler = fabs(U1[i])*1e-5+1e-12;
      ensure_distance("U2=U1",U2[i],U1[i],toler);
      if(i!=Levermore1D_Vector::get_length()) {
	ensure_distance("U2=U1",U2[i+1],U1[i+1],toler);
      }
    }

  }

  /* Test 4:*/
  template<>
  template<>
  void Levermore1D_weights_object::test<4>()
  {
    set_test_name("Assignment Operator");

    int i(0);
    Levermore1D_weights A1, A2, A3;

    for(i=1;i<=Levermore1D_weights::get_length();++i) {
      A1[i] = pow((double)i,1.23456) / 98.765;
    }

    A3 = A2 = A1;

    for(i=1;i<=Levermore1D_weights::get_length();++i) {
      ensure_distance("A2=A1",A2[i],A1[i],fabs(A1[i])*tol);
      ensure_distance("A3=A1",A3[i],A1[i],fabs(A1[i])*tol);
    }
  }

  /* Test 5:*/
  template<>
  template<>
  void Levermore1D_weights_object::test<5>()
  {
    set_test_name("Addition Operator");

    int i(0);
    double a(0.0), b(0.0);
    Levermore1D_weights A1, A2, A3;

    for(i=1;i<=Levermore1D_weights::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp((double)i);
      A1[i] = a;
      A2[i] = b;
    }

    A3 = A1 + A2;

    for(i=1;i<=Levermore1D_weights::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp((double)i);
      ensure_distance("A3 = A1 + A2", A3[i], a+b, fabs(a+b)*tol);
    }
  }

  /* Test 6:*/
  template<>
  template<>
  void Levermore1D_weights_object::test<6>()
  {
    set_test_name("Addition and Assignment Operator");

    int i(0);
    double a(0.0), b(0.0);
    Levermore1D_weights A1, A2, A3;

    for(i=1;i<=Levermore1D_weights::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp((double)i);
      A1[i] = a;
      A2[i] = b;
    }

    A3 = (A2 += A1);

    for(i=1;i<=Levermore1D_weights::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp((double)i);

      ensure_distance("A2 += A1",A2[i],a+b,fabs(a+b)*tol);
      ensure_distance("A3 = (A2+=A1)",A3[i],a+b,fabs(a+b)*tol);
    }
  }

  /* Test 7:*/
  template<>
  template<>
  void Levermore1D_weights_object::test<7>()
  {
    set_test_name("Subtraction Operator");

    int i(0);
    double a(0.0), b(0.0);
    Levermore1D_weights A1, A2, A3;

    for(i=1;i<=Levermore1D_weights::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp((double)i);
      A1[i] = a;
      A2[i] = b;
    }

    A3 = A1 - A2;

    for(i=1;i<=Levermore1D_weights::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp((double)i);
      ensure_distance("A3 = A1 - A2", A3[i], a-b, fabs(a-b)*tol);
    }
  }

  /* Test 8:*/
  template<>
  template<>
  void Levermore1D_weights_object::test<8>()
  {
    set_test_name("Subtraction and Assignment Operator");

    int i(0);
    double a(0.0), b(0.0);
    Levermore1D_weights A1, A2, A3;

    for(i=1;i<=Levermore1D_weights::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp((double)i);
      A1[i] = a;
      A2[i] = b;
    }

    A3 = (A2 -= A1);

    for(i=1;i<=Levermore1D_weights::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp((double)i);

      ensure_distance("A2 -= A1",A2[i],b-a,fabs(b-a)*tol);
      ensure_distance("A3 = (A2-=A1)",A3[i],b-a,fabs(b-a)*tol);
    }
  }

  /* Test 9:*/
  template<>
  template<>
  void Levermore1D_weights_object::test<9>()
  {
    set_test_name("Dot Product Operator");

    int i(0);
    double a(0.0), b(0.0), dot(0.0), dot_levermore(0.0);
    Levermore1D_weights A1, A2;

    for(i=1;i<=Levermore1D_weights::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp((double)i);
      dot += (a*b);
      A1[i] = a;
      A2[i] = b;
    }

    dot_levermore = A1 * A2;
    ensure_distance("dot_prod = A1*A2)", dot_levermore, dot, fabs(dot)*tol);
    dot_levermore = A2 * A1;
    ensure_distance("dot_prod = A2*A1)", dot_levermore, dot, fabs(dot)*tol);

  }

  /* Test 10:*/
  template<>
  template<>
  void Levermore1D_weights_object::test<10>()
  {
    set_test_name("Term-wise Product Operator");

    int i(0);
    double a(0.0), b(0.0);
    Levermore1D_weights A1, A2, A3;

    for(i=1;i<=Levermore1D_weights::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp((double)i);
      A1[i] = a;
      A2[i] = b;
    }

    A3 = A1 ^ A2;
    for(i=1;i<=Levermore1D_weights::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp((double)i);
      ensure_distance("A3 = A1^A2)", A3[i], a*b, fabs(a*b)*tol);
    }

    A3 = A2 ^ A1;
    for(i=1;i<=Levermore1D_weights::get_length();++i) {
      a = pow((double)i,1.23456) / 98.765;
      b = sqrt((double)i) * exp((double)i);
      ensure_distance("A3 = A1^A2)", A3[i], a*b, fabs(a*b)*tol);
    }
  }

  /* Test 11:*/
  template<>
  template<>
  void Levermore1D_weights_object::test<11>()
  {
    set_test_name("copy_form Function");
    int i(0);
    Levermore1D_weights A1, A2;

    for(i=1;i<=Levermore1D_weights::get_length();++i) {
      A1[i] = pow((double)i,1.23456) / 98.765;
    }

    A2.copy_from(A1);

    for(i=1;i<=Levermore1D_weights::get_length();++i) {
      ensure_distance("A2=A1",A2[i],A1[i],fabs(A1[i])*tol);
    }
  }

  /* Test 12:*/
  template<>
  template<>
  void Levermore1D_weights_object::test<12>()
  {
    set_test_name("zero, one, and set_all Function");
    int i(0);
    Levermore1D_weights A1, A2, A3;

    for(i=1;i<=Levermore1D_weights::get_length();++i) {
      A1[i] = pow((double)i,1.23456) / 98.765;
      A2[i] = pow((double)i,1.23456) / 98.765;
      A3[i] = pow((double)i,1.23456) / 98.765;
      ensure("A1!=0.0", A1[i] != 0.0);
      ensure("A2!=1.0", A2[i] != 1.0);
      ensure("A3!=2.0", A3[i] != 2.0);
    }

    A1.zero();
    A2.one();
    A3.set_all(2.0);

    for(i=1;i<=Levermore1D_weights::get_length();++i) {
      ensure_distance("A1 = 0.0",A1[i],0.0,tol);
      ensure_distance("A2 = 1.0",A2[i],1.0,tol);
      ensure_distance("A3 = 2.0",A3[i],2.0,tol);
    }
  }

  /* Test 13:*/
  template<>
  template<>
  void Levermore1D_weights_object::test<13>()
  {
    set_test_name("value_at");

    Levermore1D_weights alpha;
    double exponent, expected;
    double a1(1.223), a2(0.552), a3(0.456), a4(-0.223), a5(0.003);
    double v(1.3), us(50.0); //anything for us

    alpha.zero();

    alpha[1] = a1;
    alpha[2] = a2;
    alpha[3] = a3;
    if(Levermore1D_Vector::get_length() > 3) {
      alpha[4] = a4;
      alpha[5] = a5;
    }

    exponent = a1 + a2*v + a3*v*v - STABILIZATION*pow(v-us,Levermore1D_Vector::get_length()+1);
    if(Levermore1D_Vector::get_length() > 3) {
      exponent += a4*v*v*v + a5*v*v*v*v;
    }

    expected = exp(exponent);

    ensure_distance("value_at",expected,alpha.value_at(v,us),fabs(expected)*tol);

  }

  /* Test 14:*/
  template<>
  template<>
  void Levermore1D_weights_object::test<14>()
  {
    set_test_name("Generate domain of realizability for 5-moment system");

//    if(Levermore1D_Vector::get_length() == 5) {
//
//      ofstream fout;
//      fout.open("L1D_5M_realizability.dat");
//      fout.precision(16);
//
//      char output[4];
//      output[0] = '|';
//      output[1] = '/';
//      output[2] = '-';
//      output[3] = '\\';
//      int out_count(0);
//      cout << output[out_count%4]; cout.flush();
//      out_count++;
//
//      double rho(1.0), p(1.0), q, r, r_min;
//      double step(0.005), factor(1.05);
//
//      Levermore1D_cState U;
//      Levermore1D_weights A;
//      U.zero();
//      U[1] = rho;
//      U[3] = p;
//
//      for(q = -2.0; q<=2.0; q+=step) {
//	r_min = factor*(p*p/rho+q*q/p);
//	r = r_min;
//	U[4] = q;
//	U[5] = r;
//	A.MaxBoltz(U);
//	while( ! A.set_from_U(U) ) {
//	  r += step;
//	  U[5] = r;
//	  A.MaxBoltz(U);
//	}
//	if(r!=r_min) r -= step;
//	fout << q << "   " << r_min << "   " << r << endl;
//	cout << "\b" << output[out_count%4]; cout.flush();
//	out_count++;
//      }
//      cout << "\b";
//    }
  }

  /* Test 15:*/
  template<>
  template<>
  void Levermore1D_weights_object::test<15>()
  {
    set_test_name("Set from U. I was debugging the length=5 case.");

    if(Levermore1D_Vector::get_length() == 5) {

      Levermore1D_cState U, U2;
      U.zero();
      U[1] = 1.0;
      U[3] = 1.0e5;
      U[5] = 3.5e10;

      //U = Levermore1D_cState(Aa);

      Levermore1D_weights A(U);
      U2.set_from_A(A,U[2]/U[1]);

      for(int i=1;i<=Levermore1D_Vector::get_length();++i) {
	ensure_distance("U2=U",U2[i],U[i],fabs(U[i])*1e-6+1e-10);
      }

    }
  }

  /* Test 16:*/
  template<>
  template<>
  void Levermore1D_weights_object::test<16>()
  {
    set_test_name("Generate output (temporary test).");
//    int i(0), j(0);
//    double Func(0.0), diff(0.0);
//    Levermore1D_weights A;
//    Levermore1D_cState U, Ud; //Ud = U desired
//    ofstream fout;
//    fout.open("Levermore.out");
//    fout.precision(16);
//
//    Ud.zero();
//    Ud[1] = 1.0;
//    Ud[3] = 1.0;
//    Ud[5] = 5.5;
//
//    for(double C1=-0.01; fabs(C1)<1.0e0; C1-=0.01) {
//      j++;
//    }
//    for(double C2=-0.01; fabs(C2)<1.0e1; C2-=0.01) {
//      i++;
//    }
//
//    fout << "TITLE = \"" << CFFC_Name() << ": 1D Levermore Solution "
//	 << "VARIABLES = \"C1\" \\ \n"
//	 << "\"C2\" \\ \n"
//	 << "\"Func\" \\ \n"
//	 << "\"Diff\" \\ \n"
//	 << "\"rho\" \\ \n"
//	 << "\"p\" \\ \n"
//	 << "\"r\" \n"
//	 << "ZONE T =  Zone1 \\ \n"
//	 << "I = " << i << " \\ \n"
//	 << "J = " << j << " \\ \n"
//	 << "F = POINT \n";
//
//    for(double C1=-0.01; fabs(C1)<1.0e0; C1-=0.01) {
//      for(double C2=-0.01; fabs(C2)<1.0e1; C2-=0.01) {
//	A.zero();
//	A[3] = C2;
//	A[5] = C1;
//
//	U = Levermore1D_cState(A);
//	A[1] = -log(U[1]);
//	U = Levermore1D_cState(A);
//
//	Func = U[1]-Ud[1]*A[1]-Ud[3]*A[3]-Ud[5]*A[5];
//	diff = fabs( (U[1]-Ud[1])/Ud[1]) +
//	  fabs( (U[3]-Ud[3])/Ud[3]) +
//	  fabs( (U[5]-Ud[5])/Ud[5]);
//	fout << C1 << "  " << C2 << "  " << Func << "  " << diff << " " << U[1] << "  " << U[3] << "  " << U[5] << endl;
//
//      }
//    }
//
//    fout << endl << "ZONE T =  Zone2 \\ \n"
//	 << "I = " << i << " \\ \n"
//	 << "J = " << j << " \\ \n"
//	 << "F = POINT \n";
//
//    for(double C1=-0.01; fabs(C1)<1.0e0; C1-=0.01) {
//      for(double C2=0.01; fabs(C2)<1.0e1; C2+=0.01) {
//	A.zero();
//	A[3] = C2;
//	A[5] = C1;
//
//	U = Levermore1D_cState(A);
//	A[1] = -log(U[1]);
//	U = Levermore1D_cState(A);
//
//	Func = U[1]-Ud[1]*A[1]-Ud[3]*A[3]-Ud[5]*A[5];
//	diff = fabs( (U[1]-Ud[1])/Ud[1]) +
//	  fabs( (U[3]-Ud[3])/Ud[3]) +
//	  fabs( (U[5]-Ud[5])/Ud[5]);
//	fout << C1 << "  " << C2 << "  " << Func << "  " << diff << " " << U[1] << "  " << U[3] << "  " << U[5] << endl;
//
//      }
//    }
//
//
  }


  /* Test 17: Galilean invariant stabilization?*/
  template<>
  template<>
  void Levermore1D_weights_object::test<17>()
  {
    set_test_name("Check if stabilization is applied in a Galilean invariant way.");

    double rho, u, p;
    char testname[256];
    double moment1, moment2;

    rho = 1.225;
    u = 555.221;
    p = 101325.0;

    Levermore1D_pState W1(rho,ZERO,p), W2(rho,u,p);

    Levermore1D_weights A1(W1), A2(W2);

    for(int i=Levermore1D_Vector::get_length()+1; i<15; i=i+2) {
      sprintf(testname, "Integrate random moment %d.", i);
      moment1 = A1.integrate_random_moment(i,0.0,0.0);
      moment2 = A2.integrate_random_moment(i,u,u);
      ensure_distance(testname, moment1, moment2, fabs(moment2)*1e-6+1e-6);
    }


  }

  /* Test 18:*/
  template<>
  template<>
  void Levermore1D_weights_object::test<18>()
  {
    set_test_name("Set from U. I was debugging the length=5 case.");

    if(Levermore1D_Vector::get_length() == 5) {

      Levermore1D_cState U, U2;
      U.zero();
      U[1] = 0.1967128225839452;
      U[2] = 12.73347474627599;
      U[3] = 15157.91196934675;
      U[4] = 3385491.534787332;
      U[5] = 3676070793.14471;

      Levermore1D_weights A(U);

      U2.set_from_A(A,U[2]/U[1]);

      Levermore1D_pState W(A,U[2]/U[1]);

      for(int i=1;i<=Levermore1D_Vector::get_length();++i) {
	ensure_distance("U2=U",U2[i],U[i],fabs(U[i])*1e-6+1e-10);
      }

    }
  }

  /* Test 19:*/
  template<>
  template<>
  void Levermore1D_weights_object::test<19>()
  {
    set_test_name("Set from U. I was debugging the length=5 case.");

    if(Levermore1D_Vector::get_length() == 5) {

      Levermore1D_cState U, U2;
      U.zero();
      U[1] = 0.6956844155639194;
      U[2] = 114.6564144829914;
      U[3] = 57340.99014363364;
      U[4] = 21254866.53602972;
      U[5] = 14036286526.45919;

      Levermore1D_weights A;

      A[1] = -6.781748028337478;
      A[2] = 0.002939120747963637;
      A[3] = -1.058450322321528e-05;
      A[4] = -1.967976145090454e-09;
      A[5] = 4.743970027691355e-12;

      A.set_from_U(U);

      U2.set_from_A(A,U[2]/U[1]);

      Levermore1D_pState W(A,U[2]/U[1]);

      for(int i=1;i<=Levermore1D_Vector::get_length();++i) {
	ensure_distance("U2=U",U2[i],U[i],fabs(U[i])*1e-6+1e-10);
      }

    }
  }

  /* Test 20:*/
  template<>
  template<>
  void Levermore1D_weights_object::test<20>()
  {
    set_test_name("Set from U. I was debugging the length=5 case.");

    if(Levermore1D_Vector::get_length() == 5) {

      Levermore1D_cState U, U2;
      U.zero();
      U[1] = 0.6935435631827335;
      U[2] = 117.203654476427;
      U[3] = 56976.78791740329;
      U[4] = 21171789.15245417;
      U[5] = 14229288434.8971;

      Levermore1D_weights A;

      A[1] = -6.803976689981532;
      A[2] = 0.003215016590541076;
      A[3] = -1.079347487660824e-05;
      A[4] = -2.81495124853346e-09;
      A[5] = 5.408636132588052e-12;

//      A[1] = -6.781748028337478;
//      A[2] = 0.002939120747963637;
//      A[3] = -1.058450322321528e-05;
//      A[4] = -1.967976145090454e-09;
//      A[5] = 4.743970027691355e-12;

      A.set_from_U(U);

      U2.set_from_A(A,U[2]/U[1]);

      Levermore1D_pState W(A,U[2]/U[1]);

      for(int i=1;i<=Levermore1D_Vector::get_length();++i) {
	ensure_distance("U2=U",U2[i],U[i],fabs(U[i])*1e-6+1e-10);
      }

    }
  }

  //end tests
}



// Test suite constructor
tut::Levermore1D_weights_TestSuite Levermore1D_weightsTestSuite("Class:Levermore1D_weights");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

