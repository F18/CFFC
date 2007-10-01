/*!\file test_Euler3DPolytropicState.cc
  \brief Regression tests for class 'CLASS_NAME' . */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "Euler3DPolytropicState.h"
#include "../Utilities/Utilities.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_Euler3DPolytropicState : public TestData {

    // Local variables
  public:

    // Constructor
    Data_Euler3DPolytropicState(){
      /* Set the global path for this test suite. 
	 It's a good practice to put it in the constructor of the data class in order to be set
	 automatically for each individual test. Declare it relative to the /src_2D directory,
	 otherwise the framework might not find the input and output files. */
      
      set_test_suit_path("Euler");
    }

  private:
    
  };

void ensure_distance_cState(string assertion, const Euler3D_Polytropic_cState &U1, const Euler3D_Polytropic_cState &U2, double tol) {
	ensure_distance (assertion, U1*U1 , U2*U2 , tol );
}

void ensure_distance_pState(string assertion, const Euler3D_Polytropic_pState &W1, const Euler3D_Polytropic_pState &W2, double tol) {
	ensure_distance (assertion, W1*W1 , W2*W2 , tol );
}

void ensure_distance_DenseMatrix(string assertion, DenseMatrix &M1, DenseMatrix M2, double tol) {
	ensure("dimensions", M1.get_n() == M2.get_n());
	ensure("dimensions", M1.get_m() == M2.get_m());
	int n = M1.get_m(), m = M1.get_n();
	for (int i=0; i<m; i++) {
		for (int j=0; j<n; j++) {
			ensure_distance(assertion, M1(i,j) , M2(i,j) , tol);
		}
	}
		
			
}


  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. 'NAME'_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_Euler3DPolytropicState> Euler3DPolytropicState_TestSuite;
  typedef Euler3DPolytropicState_TestSuite::object Euler3DPolytropicState_object;


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
  void Euler3DPolytropicState_object::test<1>()
  {
	
    set_test_name("Check constructors");
	
	double tol = 1e-10;
	
	/* pState */
	Euler3D_Polytropic_pState Wcr;
	ensure ("creation W.rho", Wcr.rho == DENSITY_STDATM);
	ensure ("creation W.v" , Wcr.v == Vector3D(ZERO,ZERO,ZERO));
	ensure ("creation W.p" , Wcr.p == PRESSURE_STDATM); 
		  
	double rho = 2.0;
	Vector3D v(1.0,2.0,3.0);
	double p = 4.0;

	Euler3D_Polytropic_pState W(rho, v.x, v.y, v.z, p);
	ensure ("assignment W.rho", W.rho == rho);
	ensure ("assignment W.v.x", W.v.x == v.x);
	ensure ("assignment W.v.y", W.v.y == v.y);
	ensure ("assignment W.v.z", W.v.z == v.z);
	ensure ("assignment W.p" , W.p == p);
	
	Euler3D_Polytropic_pState Wco(W);
	ensure ("copy W", Wco == W);
	
	Euler3D_Polytropic_pState Was(rho, Vector3D(v.x,v.y,v.z), p);
	ensure ("assignment W", Was == W);
	
	/* cState */
	Euler3D_Polytropic_cState Ucr;
	ensure ("creation U.rho", Ucr.rho == DENSITY_STDATM);
	ensure ("creation U.rhov" , Ucr.rhov == Vector3D(ZERO,ZERO,ZERO));
	ensure ("creation U.E" , Ucr.E == PRESSURE_STDATM/(GAMMA_AIR-ONE)); 
	
	Vector3D rhov = rho*v;
	double E = p/(GAMMA_AIR-ONE) + HALF*rho*v.sqr();
	Euler3D_Polytropic_cState U(rho, rhov.x, rhov.y, rhov.z, E);
	ensure ("assignment U.rho", U.rho == rho);
	ensure ("assignment U.rhov.x", U.rhov.x == rhov.x);
	ensure ("assignment U.rhov.y", U.rhov.y == rhov.y);
	ensure ("assignment U.rhov.z", U.rhov.z == rhov.z);
	ensure ("assignment U.E" , U.E == E);
	
	Euler3D_Polytropic_cState Uas(rho, Vector3D(rhov.x,rhov.y,rhov.z), E);
	ensure ("assignment U", Uas == U);
	
	Euler3D_Polytropic_cState Uco(U);
	ensure ("copy U", Uco == U);
	
	/* check if W == U.W() and U == W.U() */
	ensure_distance_pState ("assignment W -> U", W , U.W() , tol );
	ensure_distance_cState ("assignment U -> W", U , W.U() , tol );

  }








 /* Test 2:*/
  template<>
  template<>
  void Euler3DPolytropicState_object::test<2>()
  {
	
    set_test_name("Check static variables");
	
	double tol = 1e-10;
	
	double rho = 2.0;
	Vector3D v(1.0,2.0,3.0);
	double p = 4.0;

	Euler3D_Polytropic_pState W(rho, v.x, v.y, v.z, p);
	Euler3D_Polytropic_cState U = W.U();

	ensure_distance("p", W.p , U.p() , tol);
	ensure_distance("h", W.h() , U.h() , tol);
	ensure_distance("H", W.H() , U.H() , tol);
	ensure_distance("a", W.a() , U.a() , tol);
	ensure_distance("a2", W.a2() , U.a2() , tol);
	ensure("v", W.v == U.v());
	ensure("rhov", W.rhov() == U.rhov);
	ensure_distance("e", W.e() , U.e() , tol);
	ensure_distance("E", W.E() , U.E , tol);
	ensure_distance("M", W.M() , U.M() , tol);
	ensure_distance("s", W.s() , U.s() , tol);
	
	ensure("enthalpy W", W.H()/W.rho == W.h());
	ensure("enthalpy U", U.H()/U.rho == U.h());
	
	ensure("enthalpy WH", W.H() == W.E() + W.p );

  }







 /* Test 3:*/
  template<>
  template<>
  void Euler3DPolytropicState_object::test<3>()
  {
	
    set_test_name("Check operators");
	
	double tol = 1e-10;
	
	double rho = 2.0;
	Vector3D v(1.0,2.0,3.0);
	double p = 4.0;
	Vector3D rhov = rho*v;
	double E =  p/(GAMMA_AIR-ONE) + HALF*rho*v.sqr();
	double a=3.0;

	Euler3D_Polytropic_pState W(rho, v.x, v.y, v.z, p);
	Euler3D_Polytropic_cState U = W.U();

	ensure_distance("W[1]", W.rho , W[1] , tol);
	ensure_distance("W[2]", W.v.x , W[2] , tol);
	ensure_distance("W[3]", W.v.y , W[3] , tol);
	ensure_distance("W[4]", W.v.z , W[4] , tol);
	ensure_distance("W[5]", W.p   , W[5] , tol);
	
	ensure_distance("U[1]", U.rho    , U[1] , tol);
	ensure_distance("U[2]", U.rhov.x , U[2] , tol);
	ensure_distance("U[3]", U.rhov.y , U[3] , tol);
	ensure_distance("U[4]", U.rhov.z , U[4] , tol);
	ensure_distance("U[5]", U.E      , U[5] , tol);
	
	ensure_distance_pState("W+W" , W+W , Euler3D_Polytropic_pState(rho+rho,v+v,p+p) , tol);
	ensure_distance_pState("W-W" , W-W , Euler3D_Polytropic_pState(rho-rho,v-v,p-p) , tol);
	ensure_distance		  ("W*W" , W*W , rho*rho+v*v+p*p , tol);
	ensure_distance_pState("W^W" , W^W , Euler3D_Polytropic_pState(rho*rho,v.x*v.x,v.y*v.y,v.z*v.z,p*p) , tol);
	ensure_distance_pState("W*a" , W*a , Euler3D_Polytropic_pState(rho*a,v*a,p*a) , tol);
	ensure_distance_pState("a*W" , a*W , Euler3D_Polytropic_pState(a*rho,a*v,a*p) , tol);

	ensure_distance_cState("U+U" , U+U , Euler3D_Polytropic_cState(rho+rho,rhov.x+rhov.x,rhov.y+rhov.y,rhov.z+rhov.z,E+E) , tol);
	ensure_distance_cState("U-U" , U-U , Euler3D_Polytropic_cState(rho-rho,rhov-rhov,E-E) , tol);
	ensure_distance		  ("U*U" , U*U , rho*rho+rhov*rhov+E*E , tol);
	ensure_distance_cState("U^U" , U^U , Euler3D_Polytropic_cState(rho*rho,rhov.x*rhov.x,rhov.y*rhov.y,rhov.z*rhov.z,E*E) , tol);
	ensure_distance_cState("U*a" , U*a , Euler3D_Polytropic_cState(rho*a,rhov*a,E*a) , tol);
	ensure_distance_cState("a*U" , a*U , Euler3D_Polytropic_cState(a*rho,a*rhov,a*E) , tol);
	
  }
  
  
   /* Test 4:*/
  template<>
  template<>
  void Euler3DPolytropicState_object::test<4>()
  {
	
    set_test_name("Check Fluxes");
	
	double tol = 1e-10;
	
	double rho = 2.0;
	Vector3D v(1.0,2.0,3.0);
	double p = 4.0;

	Euler3D_Polytropic_pState W(rho, v.x, v.y, v.z, p);
	Euler3D_Polytropic_cState U = W.U();
	Euler3D_Polytropic_cState Fx = W.Fx();	
	Euler3D_Polytropic_cState Fy = W.Fy();
	Euler3D_Polytropic_cState Fz = W.Fz();
	DenseMatrix dFxdU(5,5), dFxdU_num(5,5);
	DenseMatrix dFydU(5,5), dFydU_num(5,5);
	DenseMatrix dFzdU(5,5), dFzdU_num(5,5);
	DenseMatrix dUdW(5,5), dUdW_num(5,5);
	DenseMatrix dWdU(5,5), dWdU_num(5,5);
	dFxdU.zero();	dFxdU_num.zero();
	dFydU.zero();	dFydU_num.zero();
	dFzdU.zero();	dFzdU_num.zero();
	dUdW.zero();	dUdW_num.zero();
	dWdU.zero();	dWdU_num.zero();
	W.dFxdU(dFxdU);
	W.dFydU(dFydU);
	W.dFzdU(dFzdU);
	W.dUdW(dUdW);
	W.dWdU(dWdU);
		
	ensure("Fx[1]", Fx[1] == 2.0);
	ensure("Fx[2]", Fx[2] == 6.0);
	ensure("Fx[3]", Fx[3] == 4.0);
	ensure("Fx[4]", Fx[4] == 6.0);
	ensure("Fx[5]", Fx[5] == 28.0);
	
	ensure("Fy[1]", Fy[1] == 4.0);
	ensure("Fy[2]", Fy[2] == 4.0);
	ensure("Fy[3]", Fy[3] == 12.0);
	ensure("Fy[4]", Fy[4] == 12.0);
	ensure("Fy[5]", Fy[5] == 56.0);
	
	ensure("Fz[1]", Fz[1] == 6.0);
	ensure("Fz[2]", Fz[2] == 6.0);
	ensure("Fz[3]", Fz[3] == 12.0);
	ensure("Fz[4]", Fz[4] == 22.0);
	ensure("Fz[5]", Fz[5] == 84.0);
	
	ensure_distance_cState("U.Fx()=W.Fx()" , U.Fx() , Fx , tol);
	ensure_distance_cState("U.Fy()=W.Fy()" , U.Fy() , Fy , tol);
	ensure_distance_cState("U.Fz()=W.Fz()" , U.Fz() , Fz , tol);
	
	/* Mathematica results for dFxdU */
	dFxdU_num(0,1) = 1.0;
	dFxdU_num(1,0) = 1.8;
	dFxdU_num(1,1) = 1.6;
	dFxdU_num(1,2) = -0.8;
	dFxdU_num(1,3) = -1.2;
	dFxdU_num(1,4) = 0.4;
	dFxdU_num(2,0) = -2.0;
	dFxdU_num(2,1) = 2.0;
	dFxdU_num(2,2) = 1.0;
	dFxdU_num(3,0) = -3.0;
	dFxdU_num(3,1) = 3.0;
	dFxdU_num(3,3) = 1.0;
	dFxdU_num(4,0) = -11.2;
	dFxdU_num(4,1) = 13.6;
	dFxdU_num(4,2) = -0.8;
	dFxdU_num(4,3) = -1.2;
	dFxdU_num(4,4) = 1.4;

	ensure_distance_DenseMatrix("dFxdU", dFxdU, dFxdU_num, tol);
	
	/* Mathematica results for dFydU */
	dFydU_num(0,2) = 1.0;
	dFydU_num(1,0) = -2.0;
	dFydU_num(1,1) = 2.0;
	dFydU_num(1,2) = 1.0;
	dFydU_num(2,0) = -1.2;
	dFydU_num(2,1) = -0.4;
	dFydU_num(2,2) = 3.2;
	dFydU_num(2,3) = -1.2;
	dFydU_num(2,4) = 0.4;
	dFydU_num(3,0) = -6.0;
	dFydU_num(3,2) = 3.0;
	dFydU_num(3,3) = 2.0;
	dFydU_num(4,0) = -22.4;
	dFydU_num(4,1) = -0.8;
	dFydU_num(4,2) = 12.4;
	dFydU_num(4,3) = -2.4;
	dFydU_num(4,4) = 2.8;
	
	ensure_distance_DenseMatrix("dFydU", dFydU, dFydU_num, tol);

	/* Mathematica results for dFzdU */
	dFzdU_num(0,3) = 1.0;
	dFzdU_num(1,0) = -3.0;
	dFzdU_num(1,1) = 3.0;
	dFzdU_num(1,3) = 1.0;
	dFzdU_num(2,0) = -6.0;
	dFzdU_num(2,2) = 3.0;
	dFzdU_num(2,3) = 2.0;
	dFzdU_num(3,0) = -6.2;
	dFzdU_num(3,1) = -0.4;
	dFzdU_num(3,2) = -0.8;
	dFzdU_num(3,3) = 4.8;
	dFzdU_num(3,4) = 0.4;
	dFzdU_num(4,0) = -33.6;
	dFzdU_num(4,1) = -1.2;
	dFzdU_num(4,2) = -2.4;
	dFzdU_num(4,3) = 10.4;
	dFzdU_num(4,4) = 4.2;
	
	ensure_distance_DenseMatrix("dFzdU", dFzdU, dFzdU_num, tol);
	
	/* Mathematica results for dUdW*/
	dUdW_num(0,0) = 1.0;
	dUdW_num(1,0) = 1.0;
	dUdW_num(1,1) = 2.0;
	dUdW_num(2,0) = 2.0;
	dUdW_num(2,2) = 2.0;
	dUdW_num(3,0) = 3.0;
	dUdW_num(3,3) = 2.0;
	dUdW_num(4,0) = 7.0;
	dUdW_num(4,1) = 2.0;
	dUdW_num(4,2) = 4.0;
	dUdW_num(4,3) = 6.0;
	dUdW_num(4,4) = 2.5;
	
	ensure_distance_DenseMatrix("dUdW", dUdW, dUdW_num, tol);
	
	/* Mathematica results for dUdW*/
	dWdU_num(0,0) = 1.0;
	dWdU_num(1,0) = -0.5;
	dWdU_num(1,1) = 0.5;
	dWdU_num(2,0) = -1.0;
	dWdU_num(2,2) = 0.5;
	dWdU_num(3,0) = -1.5;
	dWdU_num(3,3) = 0.5;
	dWdU_num(4,0) = 2.8;
	dWdU_num(4,1) = -0.4;
	dWdU_num(4,2) = -0.8;
	dWdU_num(4,3) = -1.2;
	dWdU_num(4,4) = 0.4;
	
	ensure_distance_DenseMatrix("dWdU", dWdU, dWdU_num, tol);
	
	
	/* Mathematica results for rc_x*/
	Euler3D_Polytropic_cState rc_x1_num(1.0,-0.6733200530681511,2.0,3.0,12.32667994693185);
	Euler3D_Polytropic_cState rc_x2_num(1.0,1.0,2.0,3.0,7.0);
	Euler3D_Polytropic_cState rc_x3_num(0.0,0.0,2.0,0.0,4.0);
	Euler3D_Polytropic_cState rc_x4_num(0.0,0.0,0.0,2.0,6.0);
	Euler3D_Polytropic_cState rc_x5_num(1.0,2.673320053068151,2.0,3.0,15.673320053068153);

	ensure_distance_cState("rc_x1",rc_x1_num,W.rc_x(1),tol);
	ensure_distance_cState("rc_x2",rc_x2_num,W.rc_x(2),tol);
	ensure_distance_cState("rc_x3",rc_x3_num,W.rc_x(3),tol);
	ensure_distance_cState("rc_x4",rc_x4_num,W.rc_x(4),tol);
	ensure_distance_cState("rc_x5",rc_x5_num,W.rc_x(5),tol);
	
	/* Mathematica results for lp_x*/
	Euler3D_Polytropic_pState lp_x1_num(0.0,-0.5976143046671968,0.0,0.0,0.17857142857142855);
	Euler3D_Polytropic_pState lp_x2_num(1.0,0.0,0.0,0.0,-0.3571428571428571);
	Euler3D_Polytropic_pState lp_x3_num(0.0,0.0,1.0,0.0,0.0);
	Euler3D_Polytropic_pState lp_x4_num(0.0,0.0,0.0,1.0,0.0);
	Euler3D_Polytropic_pState lp_x5_num(0.0,0.5976143046671968,0.0,0.0,0.17857142857142855);
	
	ensure_distance_pState("lp_x1",lp_x1_num,W.lp_x(1),tol);
	ensure_distance_pState("lp_x2",lp_x2_num,W.lp_x(2),tol);
	ensure_distance_pState("lp_x3",lp_x3_num,W.lp_x(3),tol);
	ensure_distance_pState("lp_x4",lp_x4_num,W.lp_x(4),tol);
	ensure_distance_pState("lp_x5",lp_x5_num,W.lp_x(5),tol);
	
	
	/* Mathematica results for rc_y*/
	Euler3D_Polytropic_cState rc_y1_num(1.0,1.0,0.3266799469318489,3.0,10.6533598938637);
	Euler3D_Polytropic_cState rc_y2_num(1.0,1.0,2.0,3.0,7.0);
	Euler3D_Polytropic_cState rc_y3_num(0.0,2.0,0.0,0.0,2.0);
	Euler3D_Polytropic_cState rc_y4_num(0.0,0.0,0.0,2.0,6.0);
	Euler3D_Polytropic_cState rc_y5_num(1.0,1.0,3.673320053068151,3.0,17.346640106136306);
	
	ensure_distance_cState("rc_y1",rc_y1_num,W.rc_y(1),tol);
	ensure_distance_cState("rc_y2",rc_y2_num,W.rc_y(2),tol);
	ensure_distance_cState("rc_y3",rc_y3_num,W.rc_y(3),tol);
	ensure_distance_cState("rc_y4",rc_y4_num,W.rc_y(4),tol);
	ensure_distance_cState("rc_y5",rc_y5_num,W.rc_y(5),tol);
	
	/* Mathematica results for lp_y*/
	Euler3D_Polytropic_pState lp_y1_num(0.0,0.0,-0.5976143046671968,0.0,0.17857142857142855);
	Euler3D_Polytropic_pState lp_y2_num(1.0,0.0,0.0,0.0,-0.3571428571428571);
	Euler3D_Polytropic_pState lp_y3_num(0.0,1.0,0.0,0.0,0.0);
	Euler3D_Polytropic_pState lp_y4_num(0.0,0.0,0.0,1.0,0.0);
	Euler3D_Polytropic_pState lp_y5_num(0.0,0.0,0.5976143046671968,0.0,0.17857142857142855);
	
	ensure_distance_pState("lp_y1",lp_y1_num,W.lp_y(1),tol);
	ensure_distance_pState("lp_y2",lp_y2_num,W.lp_y(2),tol);
	ensure_distance_pState("lp_y3",lp_y3_num,W.lp_y(3),tol);
	ensure_distance_pState("lp_y4",lp_y4_num,W.lp_y(4),tol);
	ensure_distance_pState("lp_y5",lp_y5_num,W.lp_y(5),tol);
	
	
	
	/* Mathematica results for rc_z*/
	Euler3D_Polytropic_cState rc_z1_num(1.0,1.0,2.0,1.3266799469318489,8.980039840795548);
	Euler3D_Polytropic_cState rc_z2_num(1.0,1.0,2.0,3.0,7.0);
	Euler3D_Polytropic_cState rc_z3_num(0.0,2.0,0.0,0.0,2.0);
	Euler3D_Polytropic_cState rc_z4_num(0.0,0.0,2.0,0.0,4.0);
	Euler3D_Polytropic_cState rc_z5_num(1.0,1.0,2.0,4.673320053068151,19.019960159204455);
	
	ensure_distance_cState("rc_z1",rc_z1_num,W.rc_z(1),tol);
	ensure_distance_cState("rc_z2",rc_z2_num,W.rc_z(2),tol);
	ensure_distance_cState("rc_z3",rc_z3_num,W.rc_z(3),tol);
	ensure_distance_cState("rc_z4",rc_z4_num,W.rc_z(4),tol);
	ensure_distance_cState("rc_z5",rc_z5_num,W.rc_z(5),tol);
	
	/* Mathematica results for lp_z*/
	Euler3D_Polytropic_pState lp_z1_num(0.0,0.0,0.0,-0.5976143046671968,0.17857142857142855);
	Euler3D_Polytropic_pState lp_z2_num(1.0,0.0,0.0,0.0,-0.3571428571428571);
	Euler3D_Polytropic_pState lp_z3_num(0.0,1.0,0.0,0.0,0.0);
	Euler3D_Polytropic_pState lp_z4_num(0.0,0.0,1.0,0.0,0.0);
	Euler3D_Polytropic_pState lp_z5_num(0.0,0.0,0.0,0.5976143046671968,0.17857142857142855);
		
	ensure_distance_pState("lp_z1",lp_z1_num,W.lp_z(1),tol);
	ensure_distance_pState("lp_z2",lp_z2_num,W.lp_z(2),tol);
	ensure_distance_pState("lp_z3",lp_z3_num,W.lp_z(3),tol);
	ensure_distance_pState("lp_z4",lp_z4_num,W.lp_z(4),tol);
	ensure_distance_pState("lp_z5",lp_z5_num,W.lp_z(5),tol);
	
	
	
	
  }


	/* Test 5:*/
	template<>
	template<>
	void Euler3DPolytropicState_object::test<5>()
	{
		
		set_test_name("Check Flux functions");
		
		double tol = 1e-10;
		
		double rho = 2.0;
		Vector3D v(1.0,2.0,3.0);
		double p = 4.0;
		
		Euler3D_Polytropic_pState W1(rho, v.x, v.y, v.z, p);
		Euler3D_Polytropic_pState W2, Roe_x, Roe_y, Roe_z, Roe_n;
		Euler3D_Polytropic_pState HLLE_x, HLLE_y, HLLE_z, HLLE_n;
		W2 = 2.0*W1;

		Roe_x = Euler3D_Polytropic_pState::FluxRoe_x(W1,W2);
		Roe_y = Euler3D_Polytropic_pState::FluxRoe_y(W1,W2);
		Roe_z = Euler3D_Polytropic_pState::FluxRoe_z(W1,W2);
		Roe_n = Euler3D_Polytropic_pState::FluxRoe_n(W1,W2,Vector3D(0,1,0));
		
		HLLE_x = Euler3D_Polytropic_pState::FluxHLLE_x(W1,W2);
		HLLE_y = Euler3D_Polytropic_pState::FluxHLLE_y(W1,W2);
		HLLE_z = Euler3D_Polytropic_pState::FluxHLLE_z(W1,W2);
		HLLE_n = Euler3D_Polytropic_pState::FluxHLLE_n(W1,W2,Vector3D(0,1,0));

		ensure_distance_cState("FluxRoe_x", Roe_x , Euler3D_Polytropic_pState::FluxRoe_n(W1,W2,Vector3D(1,0,0)) , tol);
		ensure_distance_cState("FluxRoe_y", Roe_y , Euler3D_Polytropic_pState::FluxRoe_n(W1,W2,Vector3D(0,1,0)) , tol);
		ensure_distance_cState("FluxRoe_z", Roe_z , Euler3D_Polytropic_pState::FluxRoe_n(W1,W2,Vector3D(0,0,1)) , tol);

		ensure_distance_cState("FluxHLLE_x", HLLE_x , Euler3D_Polytropic_pState::FluxHLLE_n(W1,W2,Vector3D(1,0,0)) , tol);
		ensure_distance_cState("FluxHLLE_y", HLLE_y , Euler3D_Polytropic_pState::FluxHLLE_n(W1,W2,Vector3D(0,1,0)) , tol);
		ensure_distance_cState("FluxHLLE_z", HLLE_z , Euler3D_Polytropic_pState::FluxHLLE_n(W1,W2,Vector3D(0,0,1)) , tol);

	}

}



// Test suite constructor
tut::Euler3DPolytropicState_TestSuite Euler3DPolytropicStateTestSuite("Class:Euler3DPolytropicState");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

