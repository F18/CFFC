/*!\file test_MHD3DState.cc
  \brief Regression tests for MHD3D_pState and MHD3D_cState classes. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../MHD3DState.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_MHD3D : public TestData {

    // Local variables
  public:

    // Constructor
    Data_MHD3D(){
      /* Set the global path for this test suite. 
	 It's a good practice to put it in the constructor of the data class in order to be set
	 automatically for each individual test. Declare it relative to the /src_2D directory,
	 otherwise the framework might not find the input and output files. */
      
      set_test_suite_path("IdealMHD2D/UnitTests");
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
  typedef test_group<Data_MHD3D> MHD3D_TestSuite;
  typedef MHD3D_TestSuite::object MHD3D_object;


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
  void MHD3D_object::test<1>()
  {

    set_test_name("State constructors for planar pState");

    double d, vx, vy, p, B0x, B0y, B1x, B1y;

    // Initialize variables
    d = 1.224;
    vx = -1.23;
    vy = 101.3434;
    p = -1.2323e5;
    B0x = 10.23;
    B0y = 1.2323;
    B1x = 56454.3434;
    B1y = 0.00002345;

    Vector3D V(vx,vy), B0(B0x,B0y), B1(B1x,B1y);

    // Create states
    MHD3D_pState W1,
      W2(d),
      W3(d, Vector2D(vx, vy), Vector2D(B1x,B1y), p),
      W4(d, Vector2D(vx, vy), Vector2D(B1x,B1y), Vector2D(B0x,B0y), p),
      W5(d, vx, vy, B1x, B1y, p),
      W6(d, vx, vy, B1x, B1y, B0x, B0y, p, "2D");

    // == Check correct setup

    // = Default constructor =
    ensure_equals("Default constructor density", W1.d(), 1.0);
    ensure_equals("Default constructor velocity", W1.v(), Vector3D(0));
    ensure_equals("Default constructor velocity x", W1.vx(), 0.0);
    ensure_equals("Default constructor velocity y", W1.vy(), 0.0);
    ensure_equals("Default constructor pressure", W1.p(), 1.0);
    ensure_equals("Default constructor perturbative field", W1.B1(), Vector3D(0.0));
    ensure_equals("Default constructor perturbative field x", W1.B1x(), 0.0);
    ensure_equals("Default constructor perturbative field y", W1.B1y(), 0.0);
    ensure_equals("Default constructor intrinsic field", W1.B0(), Vector3D(0.0));
    ensure_equals("Default constructor intrinsic field x", W1.B0x(), 0.0);
    ensure_equals("Default constructor intrinsic field y", W1.B0y(), 0.0);

    // = Value constructor =
    ensure_equals("Value constructor density", W2.d(), d);
    ensure_equals("Value constructor velocity", W2.v(), Vector3D(d));
    ensure_equals("Value constructor velocity x", W2.vx(), d);
    ensure_equals("Value constructor velocity y", W2.vy(), d);
    ensure_equals("Value constructor pressure", W2.p(), d);
    ensure_equals("Value constructor perturbative field", W2.B1(), Vector3D(d));
    ensure_equals("Value constructor perturbative field x", W2.B1x(), d);
    ensure_equals("Value constructor perturbative field y", W2.B1y(), d);
    ensure_equals("Value constructor intrinsic field", W2.B0(), Vector3D(d));
    ensure_equals("Value constructor intrinsic field x", W2.B0x(), d);
    ensure_equals("Value constructor intrinsic field y", W2.B0y(), d);

    // = Constructor with density, velocity, pressure and perturbative magnetic field set =
    ensure_equals("W3 set constructor density", W3.d(), d);
    ensure_equals("W3 set constructor velocity", W3.v(), V);
    ensure_equals("W3 set constructor velocity x", W3.vx(), vx);
    ensure_equals("W3 set constructor velocity y", W3.vy(), vy);
    ensure_equals("W3 set constructor pressure", W3.p(), p);
    ensure_equals("W3 set constructor perturbative field", W3.B1(), B1);
    ensure_equals("W3 set constructor perturbative field x", W3.B1x(), B1x);
    ensure_equals("W3 set constructor perturbative field y", W3.B1y(), B1y);
    ensure_equals("W3 set constructor intrinsic field", W3.B0(), Vector3D(0));
    ensure_equals("W3 set constructor intrinsic field x", W3.B0x(), 0.0);
    ensure_equals("W3 set constructor intrinsic field y", W3.B0y(), 0.0);

    // = Constructor with density, velocity, pressure, perturbative and intrinsic magnetic field set =
    ensure_equals("W4 set constructor density", W4.d(), d);
    ensure_equals("W4 set constructor velocity", W4.v(), V);
    ensure_equals("W4 set constructor velocity x", W4.vx(), vx);
    ensure_equals("W4 set constructor velocity y", W4.vy(), vy);
    ensure_equals("W4 set constructor pressure", W4.p(), p);
    ensure_equals("W4 set constructor perturbative field", W4.B1(), B1);
    ensure_equals("W4 set constructor perturbative field x", W4.B1x(), B1x);
    ensure_equals("W4 set constructor perturbative field y", W4.B1y(), B1y);
    ensure_equals("W4 set constructor intrinsic field", W4.B0(), B0);
    ensure_equals("W4 set constructor intrinsic field x", W4.B0x(), B0x);
    ensure_equals("W4 set constructor intrinsic field y", W4.B0y(), B0y);

    // = Constructor with all components set and zero intrinsic magnetic field
    ensure_equals("W5 set constructor density", W5.d(), d);
    ensure_equals("W5 set constructor velocity", W5.v(), V);
    ensure_equals("W5 set constructor velocity x", W5.vx(), vx);
    ensure_equals("W5 set constructor velocity y", W5.vy(), vy);
    ensure_equals("W5 set constructor pressure", W5.p(), p);
    ensure_equals("W5 set constructor perturbative field", W5.B1(), B1);
    ensure_equals("W5 set constructor perturbative field x", W5.B1x(), B1x);
    ensure_equals("W5 set constructor perturbative field y", W5.B1y(), B1y);
    ensure_equals("W5 set constructor intrinsic field", W5.B0(), Vector3D(0));
    ensure_equals("W5 set constructor intrinsic field x", W5.B0x(), 0.0);
    ensure_equals("W5 set constructor intrinsic field y", W5.B0y(), 0.0);

    // = Constructor with all components set =
    ensure_equals("All components set constructor density", W6.d(), d);
    ensure_equals("All components set constructor velocity", W6.v(), V);
    ensure_equals("All components set constructor velocity x", W6.vx(), vx);
    ensure_equals("All components set constructor velocity y", W6.vy(), vy);
    ensure_equals("All components set constructor pressure", W6.p(), p);
    ensure_equals("All components set constructor perturbative field", W6.B1(), B1);
    ensure_equals("All components set constructor perturbative field x", W6.B1x(), B1x);
    ensure_equals("All components set constructor perturbative field y", W6.B1y(), B1y);
    ensure_equals("All components set constructor intrinsic field", W6.B0(), B0);
    ensure_equals("All components set constructor intrinsic field x", W6.B0x(), B0x);
    ensure_equals("All components set constructor intrinsic field y", W6.B0y(), B0y);

  }


  /* Test 2:*/
  template<>
  template<>
  void MHD3D_object::test<2>()
  {

    set_test_name("Roe and HLLE solver for Dai-Woodward IVP");

    MHD3D_pState Wl, Wr;	// Left and right interface solutions
    MHD3D_pState W_exact;	/* Exact solution at x/t = 0, determined by an exact
				   Riemann solver (Ryu & Jones, 1995; M. Torrilhon, 2003) */
    MHD3D_cState FExact, FRoe, FHLLE;

    // Set Left state
    Wl = MHD3D_pState(1.08,
		      1.2, 0.01, 0.5,
		      0.5641895835477563, 1.0155412503859613, 0.5641895835477563,
		      ZERO, ZERO, ZERO,
		      0.95);

    // Set Right state
    Wr = MHD3D_pState(1,
		      0, 0, 0,
		      0.5641895835477563, 1.1283791670955126, 0.5641895835477563,
		      ZERO, ZERO, ZERO,
		      1);

    // Set exact solution
    W_exact = MHD3D_pState(1.4903369980,
			   0.6058783288, 0.1123507383, 0.5568615213,
			   0.5641895835, 1.4383158354, 0.7990643530,
			   ZERO, ZERO, ZERO,
			   1.6557697598);

    // Set exact flux
    FExact = MHD3D_cState(9.0296288969704908e-01,
			  3.3973286046119693e+00, -7.1003426480081577e-01, 5.2001503925375836e-02,
			  0.0000000000000000e+00, 8.0805727839133357e-01, 1.6996030502986992e-01,
			  0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
			  4.1174729436401369e+00);

    // Set Roe flux
    FRoe = MHD3D_cState(8.8840343980143133e-01,
			3.5083917212375098e+00, -7.0701164371581227e-01, 1.4962535682394074e-01,
			0.0000000000000000e+00, 7.4862081866183550e-01, 1.5076583453298822e-01,
			0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
			4.0604269556077419e+00);

    // Set HLLE flux
    FHLLE = MHD3D_cState(9.3125015277644130e-01,
			 3.6002140639449891e+00, -5.7670381847206553e-01, 5.4751950995068110e-01,
			 0.0000000000000000e+00, 7.2136671694854060e-01, 2.6426893592634448e-01,
			 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
			 4.1869629297799671e+00);

    // === Check fluxes ===
    ensure_distance("Exact flux", W_exact.F(), FExact, AcceptedError(FExact));
    ensure_distance("Roe flux", FluxRoe(Wl,Wr), FRoe, AcceptedError(FRoe));
    ensure_distance("HLLE flux", FluxHLLE(Wl,Wr), FHLLE, AcceptedError(FHLLE));

  }


}



// Test suite constructor
tut::MHD3D_TestSuite MHD3DTestSuite("Class:MHD3D_States");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

