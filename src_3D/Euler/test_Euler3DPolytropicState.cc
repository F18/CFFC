/*!\file test_Euler3DPolytropicState.cc
  \brief Regression tests for class 'CLASS_NAME' . */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "Euler3DPolytropicState.h"

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
	Euler3D_Polytropic_cState F = W.F();
		
	ensure("F[1]", F[1] == 2.0);
	ensure("F[2]", F[2] == 6.0);
	ensure("F[3]", F[3] == 4.0);
	ensure("F[4]", F[4] == 6.0);
	ensure("F[5]", F[5] == 28.0);
	
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

