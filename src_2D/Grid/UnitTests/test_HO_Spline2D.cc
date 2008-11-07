/*!\file test_HO_Spline2D.cc
  \brief Regression tests for class Spline2D_HO. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../HO_Spline2D.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_Spline2D_HO : public TestData {

    // Local variables
  public:

    // === Member variables
    Spline2D_HO S, S1;
    Spline2D_HO Curve;
    Spline2D_HO Circle1, Circle2;
    Spline2D_HO RectangleCurve;

    Vector2D NumericNormal, NumericTangent;
    Vector2D P;
    Vector2D V1, V2;
    Vector2D AnalyticTangent, AnalyticNormal;

    double Radius, Angle1, Angle2;
    double theta;
    double degree;		// to transform an angle given in degrees to rad, just multiply with 'degree'
    int NumSplinePoints;
    int TotalPoints;

    double Length, Width;

    LinkedList<Vector2D> Vlist;


    // === Member functions

    // Constructor
    Data_Spline2D_HO(void);

    Vector2D SplinePoint(const int &i);
    void InitializeSpline(Spline2D_HO &S);
    void InitializeSpline(Spline2D_HO &S, const int &N);
    void PrintPathlength(Spline2D_HO &S);    

  private:
    
  };
  
  Data_Spline2D_HO::Data_Spline2D_HO(void){
    /* Set the global path for this test suite. 
       It's a good practice to put it in the constructor of the data class in order to be set
       automatically for each individual test. Declare it relative to the /src_2D directory,
       otherwise the framework might not find the input and output files. */
    
    set_test_suite_path("Grid/UnitTests");
    set_local_input_path("HO_Spline2D");
    set_local_output_path("HO_Spline2D");

    Radius = 1.0232323e-1;
    degree = TWO*PI/360.0;
    Angle1 = 0.0;
    Angle2 = 90.0;
    NumSplinePoints = 10;
    Length = 3.0;
    Width  = 2.4;
    TotalPoints = 0;

    // Reset solid body counter
    Spline2D_HO::ResetCounter();
  }

  Vector2D Data_Spline2D_HO::SplinePoint(const int &i){
    return Vector2D(i*0.34345, (i+2.2323)*0.3434);
  }

  void Data_Spline2D_HO::InitializeSpline(Spline2D_HO & S){
    S.settype(SPLINE2D_QUINTIC);
    S.setFluxCalcMethod(ReconstructionBasedFlux); 
    
    for (int i=0; i<=S.np-1; ++i){
      S.Xp[i] = SplinePoint(i);
      S.tp[i] = SPLINE2D_POINT_NORMAL;
    }
    
    // adjust the ends
    S.tp[0] = S.tp[S.np-1] = SPLINE2D_POINT_SHARP_CORNER;
    S.pathlength();
    // set BC type
    S.setBCtype(BC_CHARACTERISTIC);
  }

  void Data_Spline2D_HO::InitializeSpline(Spline2D_HO &S, const int &N){
    S.allocate(N);
    InitializeSpline(S);
  }

  void Data_Spline2D_HO::PrintPathlength(Spline2D_HO &S){
    for (int i=0; i<=S.np-1; ++i){
      cout << "sp[" << i << "] = " << S.sp[i] << endl;
    }
  }
  

  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. 'NAME'_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_Spline2D_HO,100> Spline2D_HO_TestSuite;
  typedef Spline2D_HO_TestSuite::object Spline2D_HO_object;


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

  /* Test 1: */
  template<>
  template<>
  void Spline2D_HO_object::test<1>()
  {
    set_test_name("Constructor");

    Vector2D *Xp(NULL); double *sp(NULL);
    int *tp(NULL), *bc(NULL);

    // check
    ensure_equals("type",S.type,SPLINE2D_CONSTANT);
    ensure_equals("np",S.np,0);
    ensure_equals("getFluxCalcMethod()",S.getFluxCalcMethod(),SolveRiemannProblem);
    ensure_equals("Xp",S.Xp,Xp);
    ensure_equals("tp",S.tp,tp);
    ensure_equals("bc",S.bc,bc);
    ensure_equals("sp",S.sp,sp);
    ensure_equals("IsSolidBoundary()",S.IsSolidBoundary(),false);
    ensure_equals("getBodyID()",S.getBodyID(),0);
    ensure_equals("NumberOfSolidBodies()",S.NumberOfSolidBodies(),0);
  }

  /* Test 2: */
  template<>
  template<>
  void Spline2D_HO_object::test<2>()
  {
    set_test_name("Allocate");

    S.allocate(3,SPLINE2D_POINT_SHARP_CORNER);

    // try to give values
    S.Xp[2] = 2.0;
    S.sp[1] = 1.0;
    S.bc[0] = 3;

    // check
    ensure_equals("type",S.type,SPLINE2D_CONSTANT);
    ensure_equals("np",S.np,3);
    ensure_equals("getFluxCalcMethod()",S.getFluxCalcMethod(),SolveRiemannProblem);
    ensure_equals("Xp",S.Xp[2],2.0);
    ensure_equals("tp",S.tp[2],1);
    ensure_equals("bc",S.bc[0],3);
    ensure_equals("sp",S.sp[1],1.0);
    ensure_equals("IsSolidBoundary()",S.IsSolidBoundary(),false);
    ensure_equals("getBodyID()",S.getBodyID(),0);
    ensure_equals("NumberOfSolidBodies()",S.NumberOfSolidBodies(),0);
  }

  /* Test 3: */
  template<>
  template<>
  void Spline2D_HO_object::test<3>()
  {
    set_test_name("Destructor");

    Vector2D *Xp(NULL); double *sp(NULL);
    int *tp(NULL), *bc(NULL);

    S.allocate(3,SPLINE2D_POINT_SHARP_CORNER);

    // try to give values
    S.Xp[2] = 2.0;
    S.sp[1] = 1.0;
    S.bc[0] = 3;

    S.deallocate();

    // check
    ensure_equals("type",S.type,SPLINE2D_CONSTANT);
    ensure_equals("np",S.np,0);
    ensure_equals("getFluxCalcMethod()",S.getFluxCalcMethod(),SolveRiemannProblem);
    ensure_equals("Xp",S.Xp,Xp);
    ensure_equals("tp",S.tp,tp);
    ensure_equals("bc",S.bc,bc);
    ensure_equals("sp",S.sp,sp);
    ensure_equals("IsSolidBoundary()",S.IsSolidBoundary(),false);
    ensure_equals("getBodyID()",S.getBodyID(),0);
    ensure_equals("NumberOfSolidBodies()",S.NumberOfSolidBodies(),0);
  }

  /* Test 4: */
  template<>
  template<>
  void Spline2D_HO_object::test<4>()
  {
    set_test_name("settype()");
    S.settype(SPLINE2D_QUINTIC);
    // check
    ensure_equals("type",S.type,SPLINE2D_QUINTIC);
  }

  /* Test 5: */
  template<>
  template<>
  void Spline2D_HO_object::test<5>()
  {
    set_test_name("setFluxCalcMethod()");

    // check default value
    ensure_equals("getFluxCalcMethod()",S.getFluxCalcMethod(),SolveRiemannProblem);

    // assign new value
    S.setFluxCalcMethod(ReconstructionBasedFlux);

    // check
    ensure_equals("getFluxCalcMethod()",S.getFluxCalcMethod(),ReconstructionBasedFlux);
  }

  /* Test 6: */
  template<>
  template<>
  void Spline2D_HO_object::test<6>()
  {
    set_test_name("setBCtype()");
    S.allocate(3);
    S.setBCtype(BC_CHARACTERISTIC);
    ensure_equals("BC",S.bc[2],BC_CHARACTERISTIC);
  }

  /* Test 7: */
  template<>
  template<>
  void Spline2D_HO_object::test<7>()
  {
    set_test_name("pathlength()");
    S.allocate(3);
    S.Xp[0] = Vector2D(0.0,0.0); S.Xp[1] = Vector2D(0.0,2.33); S.Xp[2] = Vector2D(5.666,2.33);
    S.pathlength();

    ensure_equals("sp[0]",S.sp[0],0);
    ensure_equals("sp[1]",S.sp[1],2.33);
    ensure_equals("sp[2]",S.sp[2],2.33+5.666);
  }

  /* Test 8: */
  template<>
  template<>
  void Spline2D_HO_object::test<8>()
  {
    set_test_name("Copy Constructor");

    // Initialize S
    InitializeSpline(S,10);
    // Make spline a solid body boundary
    S.makeSplineSolidBoundary();

    // Create Copy
    Spline2D_HO Copy(S);

    ensure_equals("type",Copy.type,S.type);
    ensure_equals("np",Copy.np,S.np);
    ensure_equals("getFluxCalcMethod()",Copy.getFluxCalcMethod(),S.getFluxCalcMethod());
    ensure_equals("getBodyID()",S.getBodyID(), 1);
    ensure_equals("getBodyID()",Copy.getBodyID(),S.getBodyID());
    ensure_equals("NumberOfSolidBodies()", Spline2D_HO::NumberOfSolidBodies(), 1);

    for(int i=0; i<=Copy.np-1; ++i){
      ensure_equals("Xp",Copy.Xp[i],SplinePoint(i));
      ensure_equals("tp",Copy.tp[i],S.tp[i]);
      ensure_equals("bc",Copy.bc[i],S.bc[i]);
      ensure_equals("sp",Copy.sp[i],S.sp[i]);
    }
  }

  /* Test 9: */
  template<>
  template<>
  void Spline2D_HO_object::test<9>()
  {
    set_test_name("Assignment operator");

    InitializeSpline(S,10);

    // Make spline a solid body boundary
    S.makeSplineSolidBoundary();

    Spline2D_HO Copy(S);

    ensure_equals("type",Copy.type,S.type);
    ensure_equals("np",Copy.np,S.np);
    ensure_equals("getFluxCalcMethod()",Copy.getFluxCalcMethod(),S.getFluxCalcMethod());
    ensure_equals("getBodyID()",S.getBodyID(), 1);
    ensure_equals("getBodyID()",Copy.getBodyID(),S.getBodyID());
    ensure_equals("NumberOfSolidBodies()", Spline2D_HO::NumberOfSolidBodies(), 1);

    for(int i=0; i<=Copy.np-1; ++i){
      ensure_equals("Xp",Copy.Xp[i],SplinePoint(i));
      ensure_equals("tp",Copy.tp[i],S.tp[i]);
      ensure_equals("bc",Copy.bc[i],S.bc[i]);
      ensure_equals("sp",Copy.sp[i],S.sp[i]);
    }
  }

  /* Test 10: */
  template<>
  template<>
  void Spline2D_HO_object::test<10>()
  {
    set_test_name("Assignment operator");

    InitializeSpline(S,10);
    // Make spline a solid body boundary
    S.makeSplineSolidBoundary();

    Spline2D_HO Copy;
    InitializeSpline(Copy,5);

    // Assignment
    Copy = S;

    ensure_equals("type",Copy.type,S.type);
    ensure_equals("np",Copy.np,S.np);
    ensure_equals("getFluxCalcMethod()",Copy.getFluxCalcMethod(),S.getFluxCalcMethod());
    ensure_equals("getBodyID()",S.getBodyID(), 1);
    ensure_equals("getBodyID()",Copy.getBodyID(),S.getBodyID());
    ensure_equals("NumberOfSolidBodies()", Spline2D_HO::NumberOfSolidBodies(), 1);

    for(int i=0; i<=Copy.np-1; ++i){
      ensure_equals("Xp",Copy.Xp[i],SplinePoint(i));
      ensure_equals("tp",Copy.tp[i],S.tp[i]);
      ensure_equals("bc",Copy.bc[i],S.bc[i]);
      ensure_equals("sp",Copy.sp[i],S.sp[i]);
    }
  }


  /* Test 11: */
  template<>
  template<>
  void Spline2D_HO_object::test<11>()
  {
    set_test_name("Friend pathlength()");

    InitializeSpline(S,10);

    Spline2D_HO Copy(S);

    // call the friend function
    pathlength(Copy);

    ensure_equals("type",Copy.type,S.type);
    ensure_equals("np",Copy.np,S.np);
    ensure_equals("getFluxCalcMethod()",Copy.getFluxCalcMethod(),S.getFluxCalcMethod());

    for(int i=0; i<=Copy.np-1; ++i){
      ensure_equals("Xp",Copy.Xp[i],SplinePoint(i));
      ensure_equals("tp",Copy.tp[i],S.tp[i]);
      ensure_equals("bc",Copy.bc[i],S.bc[i]);
      ensure_equals("sp",Copy.sp[i],S.sp[i]);
    }
  }

  /* Test 12: */
  template<>
  template<>
  void Spline2D_HO_object::test<12>()
  {
    set_test_name("== operator");

    InitializeSpline(S,10);
    Spline2D_HO Copy(S);

    ensure_equals("From Copy Constructor", Copy == S, true);
  }

  /* Test 13: */
  template<>
  template<>
  void Spline2D_HO_object::test<13>()
  {
    set_test_name("== operator");

    InitializeSpline(S,10);
    Spline2D_HO Copy(S);

    ensure_equals("From operator =", Copy == S, true);
  }

  /* Test 14: */
  template<>
  template<>
  void Spline2D_HO_object::test<14>()
  {
    set_test_name("== operator");

    InitializeSpline(S,10);
    Spline2D_HO Copy(S);

    // modify Copy
    Copy.Xp[9] = 0.000123;
    ensure_equals("1) Xp", Copy == S, false);

    Copy = S;
    // modify Copy
    Copy.sp[5] = 0.000123;
    ensure_equals("2) sp", Copy == S, false);

    Copy = S;
    // modify Copy
    Copy.tp[3] = SPLINE2D_POINT_SHARP_CORNER;
    ensure_equals("3) tp", Copy == S, false);

    Copy = S;
    // modify Copy
    Copy.bc[8] = BC_NONE;
    ensure_equals("4) bc", Copy == S, false);

    Copy = S;
    // modify Copy
    Copy.setFluxCalcMethod(SolveRiemannProblem);
    ensure_equals("5) FluxMethod", Copy == S, false);

    Copy = S;
    // modify Copy
    InitializeSpline(Copy,5);
    ensure_equals("6) All", Copy == S, false);
  }

  /* Test 15: */
  template<>
  template<>
  void Spline2D_HO_object::test<15>()
  {
    set_test_name("== operator");

    InitializeSpline(S,10);
    Spline2D_HO Copy(S);

    // modify
    Copy.Xp[1] = 1.23e5;
    Copy.Xp[1] = S.Xp[1];
    // check
    ensure_equals("Direct copy", Copy == S, true);
  }

  /* Test 16: */
  template<>
  template<>
  void Spline2D_HO_object::test<16>()
  {
    set_test_name("!= operator");

    InitializeSpline(S,10);
    Spline2D_HO Copy(S);

    ensure_equals("From operator =", Copy != S, false);
  }

  /* Test 17: */
  template<>
  template<>
  void Spline2D_HO_object::test<17>()
  {
    set_test_name("!= operator");

    InitializeSpline(S,10);
    Spline2D_HO Copy(S);

    // modify Copy
    Copy.Xp[9] = 0.000123;
    ensure_equals("1) Xp", Copy != S, true);

    Copy = S;
    // modify Copy
    Copy.sp[5] = 0.000123;
    ensure_equals("2) sp", Copy != S, true);

    Copy = S;
    // modify Copy
    Copy.tp[3] = SPLINE2D_POINT_SHARP_CORNER;
    ensure_equals("3) tp", Copy != S, true);

    Copy = S;
    // modify Copy
    Copy.bc[8] = BC_NONE;
    ensure_equals("4) bc", Copy != S, true);

    Copy = S;
    // modify Copy
    Copy.setFluxCalcMethod(SolveRiemannProblem);
    ensure_equals("5) FluxMethod", Copy != S, true);

    Copy = S;
    // modify Copy
    InitializeSpline(Copy,5);
    ensure_equals("6) All", Copy != S, true);
  }


  /* Test 18: */
  template<>
  template<>
  void Spline2D_HO_object::test<18>()
  {
    set_test_name("Concatenation operator");
    InitializeSpline(S,10);
    InitializeSpline(S1,5);

    Spline2D_HO Copy(S), Copy1(S1);

    Curve = S + S1;
    const Vector2D Shift = S.Xp[9] - S1.Xp[0];

    // ensure that the initial splines weren't modified
    ensure_equals("S", S, Copy);
    ensure_equals("S1", S1, Copy1);

    // ensure that the concatenated spline is correct
    ensure_equals("np",Curve.np,14);
    ensure_equals("type",Curve.type,S.type);
    ensure_equals("FluxMethod",Curve.getFluxCalcMethod(),S.getFluxCalcMethod());
    for (int i=0; i<=S.np-1; ++i){
      ensure_equals("1) Xp",Curve.Xp[i],S.Xp[i]);
      ensure_equals("1) tp",Curve.tp[i],S.tp[i]);
      ensure_equals("1) bc",Curve.bc[i],S.bc[i]);
    }

    for (int i=S.np; i<=Curve.np-1; ++i){
      ensure_equals("2) Xp",Curve.Xp[i],S1.Xp[i-S.np+1] + Shift);
      ensure_equals("2) tp",Curve.tp[i],S1.tp[i-S.np+1]);
      ensure_equals("2) bc",Curve.bc[i],S1.bc[i-S.np+1]);
    }
  }

  /* Test 19: */
  template<>
  template<>
  void Spline2D_HO_object::test<19>()
  {
    set_test_name("Translate_Spline()");
    InitializeSpline(S,10);

    Spline2D_HO Copy(S);
    const Vector2D Shift(2.334, 12.234);

    // translate
    S.Translate_Spline(Shift);

    ensure_equals("np",S.np,Copy.np);
    ensure_equals("type",S.type,Copy.type);
    ensure_equals("FluxMethod",S.getFluxCalcMethod(),Copy.getFluxCalcMethod());
    for (int i=0; i<=S.np-1; ++i){
      ensure_equals("Xp",S.Xp[i],Copy.Xp[i] + Shift);
      ensure_equals("tp",S.tp[i],Copy.tp[i]);
      ensure_equals("bc",S.bc[i],Copy.bc[i]);
    }
  }

  /* Test 20: */
  template<>
  template<>
  void Spline2D_HO_object::test<20>()
  {
    set_test_name("Translate operator+");
    InitializeSpline(S,10);

    Spline2D_HO Copy(S);
    const Vector2D Shift(2.334, 12.234);

    // translate
    Curve = S + Shift;

    // Check that S is not modified
    ensure_equals("Same", S, Copy);

    // check Curve
    ensure_equals("np",Curve.np,S.np);
    ensure_equals("type",Curve.type,S.type);
    ensure_equals("FluxMethod",Curve.getFluxCalcMethod(),S.getFluxCalcMethod());
    for (int i=0; i<=S.np-1; ++i){
      ensure_equals("Xp",Curve.Xp[i],S.Xp[i] + Shift);
      ensure_equals("tp",Curve.tp[i],S.tp[i]);
      ensure_equals("bc",Curve.bc[i],S.bc[i]);
    }
  }


  /* Test 21: */
  template<>
  template<>
  void Spline2D_HO_object::test<21>()
  {
    set_test_name("Friend--> Shift operator+");
    InitializeSpline(S,10);

    Spline2D_HO Copy(S);
    const Vector2D Shift(2.334, 12.234);

    Curve = Shift + S;

    // Check that S is not modified
    ensure_equals("Same", S, Copy);

    // check Curve
    ensure_equals("np",Curve.np,10);
    ensure_equals("type",Curve.type,S.type);
    ensure_equals("FluxMethod",Curve.getFluxCalcMethod(),S.getFluxCalcMethod());
    for (int i=0; i<=S.np-1; ++i){
      ensure_equals("Xp",Curve.Xp[i],S.Xp[i] + Shift);
      ensure_equals("tp",Curve.tp[i],S.tp[i]);
      ensure_equals("bc",Curve.bc[i],S.bc[i]);
    }
  }

  /* Test 22: */
  template<>
  template<>
  void Spline2D_HO_object::test<22>()
  {
    set_test_name("Shift operator-");
    InitializeSpline(S,10);

    Spline2D_HO Copy(S);
    const Vector2D Shift(2.334, 12.234);

    Curve = S - Shift;

    // Check that S is not modified
    ensure_equals("Same", S, Copy);

    // check Curve
    ensure_equals("np",Curve.np,10);
    ensure_equals("type",Curve.type,S.type);
    ensure_equals("FluxMethod",Curve.getFluxCalcMethod(),S.getFluxCalcMethod());
    for (int i=0; i<=S.np-1; ++i){
      ensure_equals("Xp",Curve.Xp[i],S.Xp[i] - Shift);
      ensure_equals("tp",Curve.tp[i],S.tp[i]);
      ensure_equals("bc",Curve.bc[i],S.bc[i]);
    }
  }


  /* Test 23: */
  template<>
  template<>
  void Spline2D_HO_object::test<23>()
  {
    set_test_name("Shift operator+");
    InitializeSpline(S,10);

    Spline2D_HO Copy(S);
    const Vector2D Shift(2.334, 12.234);

    // operator +
    Curve = Shift + S - Shift;

    // Check that S is not modified
    ensure_equals("Same", S, Copy);

    ensure_equals("np",Curve.np,10);
    ensure_equals("type",Curve.type,Copy.type);
    ensure_equals("FluxMethod",Curve.getFluxCalcMethod(),Copy.getFluxCalcMethod());
    for (int i=0; i<=Curve.np-1; ++i){
      ensure_distance("Xp.x",Curve.Xp[i].x,Copy.Xp[i].x, tol*Copy.Xp[i].x + tol);
      ensure_distance("Xp.y",Curve.Xp[i].y,Copy.Xp[i].y, tol*Copy.Xp[i].y + tol);
      ensure_equals("tp",Curve.tp[i],Copy.tp[i]);
      ensure_equals("bc",Curve.bc[i],Copy.bc[i]);
    }
  }

  /* Test 24: */
  template<>
  template<>
  void Spline2D_HO_object::test<24>()
  {
    set_test_name("operator *");
    InitializeSpline(S,10);
    double a(2);

    Spline2D_HO Copy(S);

    // operator *
    Curve = S*a;

    // Check that S is not modified
    ensure_equals("Same", S, Copy);

    ensure_equals("np",Curve.np,10);
    ensure_equals("type",Curve.type,S.type);
    ensure_equals("FluxMethod",Curve.getFluxCalcMethod(),S.getFluxCalcMethod());
    for (int i=0; i<=S.np-1; ++i){
      ensure_equals("Xp",Curve.Xp[i],S.Xp[i]*a);
      ensure_equals("tp",Curve.tp[i],S.tp[i]);
      ensure_equals("bc",Curve.bc[i],S.bc[i]);
      ensure_equals("sp",Curve.sp[i],S.sp[i]*a);
    }
  }


  /* Test 25: */
  template<>
  template<>
  void Spline2D_HO_object::test<25>()
  {
    set_test_name("operator /");
    InitializeSpline(S,10);
    double a(2.3);

    Spline2D_HO Copy(S);

    // operator *
    Curve = S/a;

    // Check that S is not modified
    ensure_equals("Same", S, Copy);

    ensure_equals("np",Curve.np,10);
    ensure_equals("type",Curve.type,S.type);
    ensure_equals("FluxMethod",Curve.getFluxCalcMethod(),S.getFluxCalcMethod());
    for (int i=0; i<=S.np-1; ++i){
      ensure_distance("Xp",Curve.Xp[i],S.Xp[i]/a, tol*S.Xp[i] + Vector2D(tol));
      ensure_equals("tp",Curve.tp[i],S.tp[i]);
      ensure_equals("bc",Curve.bc[i],S.bc[i]);
      ensure_distance("sp",Curve.sp[i],S.sp[i]/a, tol*S.sp[i] + tol);
    }
  }

  /* Test 26: */
  template<>
  template<>
  void Spline2D_HO_object::test<26>()
  {
    set_test_name("operator * /");
    InitializeSpline(S,10);
    double a(2);

    Spline2D_HO Copy(S);

    // operator *
    Curve = (S*a)/a;

    // Check that S is not modified
    ensure_equals("Same", S, Copy);

    ensure_equals("np",Curve.np,10);
    ensure_equals("type",Curve.type,S.type);
    ensure_equals("FluxMethod",Curve.getFluxCalcMethod(),S.getFluxCalcMethod());
    for (int i=0; i<=S.np-1; ++i){
      ensure_equals("Xp",Curve.Xp[i],S.Xp[i]);
      ensure_equals("tp",Curve.tp[i],S.tp[i]);
      ensure_equals("bc",Curve.bc[i],S.bc[i]);
      ensure_equals("sp",Curve.sp[i],S.sp[i]);
    }
  }

  /* Test 27: */
  template<>
  template<>
  void Spline2D_HO_object::test<27>()
  {
    set_test_name("Scale_Spline()");
    InitializeSpline(S,10);
    double a(2);
    Spline2D_HO Copy(S);

    // Scale_Spline()
    S.Scale_Spline(a);

    ensure_equals("np",S.np,10);
    ensure_equals("type",S.type,Copy.type);
    ensure_equals("FluxMethod",S.getFluxCalcMethod(),Copy.getFluxCalcMethod());
    for (int i=0; i<=S.np-1; ++i){
      ensure_equals("Xp",S.Xp[i],Copy.Xp[i]*a);
      ensure_equals("tp",S.tp[i],Copy.tp[i]);
      ensure_equals("bc",S.bc[i],Copy.bc[i]);
      ensure_equals("sp",S.sp[i],Copy.sp[i]*a);
    }
  }

  /* Test 28: */
  template<>
  template<>
  void Spline2D_HO_object::test<28>()
  {
    set_test_name("rotate()");
    InitializeSpline(S,10);

    double a(90*PI/180);
    Spline2D_HO Copy(S);

    // rotate
    S.Rotate_Spline(a);

    ensure_equals("np",S.np,10);
    ensure_equals("type",S.type,Copy.type);
    ensure_equals("FluxMethod",S.getFluxCalcMethod(),Copy.getFluxCalcMethod());
    for (int i=0; i<=S.np-1; ++i){
      ensure_distance("Xp.x",S.Xp[i].x,-Copy.Xp[i].y,tol*(1.0 + fabs(Copy.Xp[i].y)));
      ensure_distance("Xp.y",S.Xp[i].y, Copy.Xp[i].x,tol*(1.0 + fabs(Copy.Xp[i].x)));
      ensure_equals("tp",S.tp[i],Copy.tp[i]);
      ensure_equals("bc",S.bc[i],Copy.bc[i]);
      ensure_equals("sp",S.sp[i],Copy.sp[i]);
    }
  }

  /* Test 29: */
  template<>
  template<>
  void Spline2D_HO_object::test<29>()
  {
    set_test_name("operator ^");
    InitializeSpline(S,10);

    double a(90*PI/180);
    Spline2D_HO Copy(S);

    // operator ^
    Curve = S ^ a;

    // Check that S is not modified
    ensure_equals("Same", S, Copy);

    ensure_equals("np",Curve.np,10);
    ensure_equals("type",Curve.type,Copy.type);
    ensure_equals("FluxMethod",Curve.getFluxCalcMethod(),Copy.getFluxCalcMethod());
    for (int i=0; i<=S.np-1; ++i){
      ensure_distance("Xp.x",Curve.Xp[i].x,-Copy.Xp[i].y,tol*(1.0 + fabs(Copy.Xp[i].y)));
      ensure_distance("Xp.y",Curve.Xp[i].y, Copy.Xp[i].x,tol*(1.0 + fabs(Copy.Xp[i].x)));
      ensure_equals("tp",Curve.tp[i],Copy.tp[i]);
      ensure_equals("bc",Curve.bc[i],Copy.bc[i]);
      ensure_equals("sp",Curve.sp[i],Copy.sp[i]);
    }
  }

  /* Test 30: */
  template<>
  template<>
  void Spline2D_HO_object::test<30>()
  {
    set_test_name("Copy_Spline()");

    // Initialize S
    InitializeSpline(S,10);

    //     // Create Copy
    //     Copy_Spline(Curve,S);
    
    //     ensure_equals("Copy",Curve,S);
  }

  /* Test 31: */
  template<>
  template<>
  void Spline2D_HO_object::test<31>()
  {
    set_test_name("Translate_Spline()");
    InitializeSpline(S,10);

    Spline2D_HO Copy(S);
    const Vector2D Shift(2.334, 12.234);

    // Translate
    Curve = Copy + Shift;
    S.Translate_Spline(Shift);
    
    // check
    ensure_equals("Translate", S, Curve);
  }
  
  /* Test 32: */
  template<>
  template<>
  void Spline2D_HO_object::test<32>()
  {
    set_test_name("Scale_Spline()");
    InitializeSpline(S,10);

    Spline2D_HO Copy(S);
    double a(3.5);

    // Scale
    Curve = Copy*a;
    S.Scale_Spline(a);
    
    // check
    ensure_equals("Scale", S, Curve);
  }

  /* Test 33: */
  template<>
  template<>
  void Spline2D_HO_object::test<33>()
  {
    set_test_name("Rotate_Spline()");
    InitializeSpline(S,10);

    double a(90*PI/180);
    Spline2D_HO Copy(S);

    // rotate
    Curve = S ^ a;
    S.Rotate_Spline(a);
    
    // check
    ensure_equals("Rotate", S, Curve);
  }

  /* Test 34: */
  template<>
  template<>
  void Spline2D_HO_object::test<34>()
  {
    set_test_name("operator << >>");
    InitializeSpline(S,10);
    S.makeSplineSolidBoundary();

    // check
    Check_Input_Output_Operator(S);
  }


  /* Test 35: */
  template<>
  template<>
  void Spline2D_HO_object::test<35>()
  {
    set_test_name("Reflect_Spline()");
    RunRegression = ON;
    InitializeSpline(S,10);

    MasterFile = "Reflect_Spline.dat";
    CurrentFile ="Current_Reflect_Spline.dat";

    // reflect
    S.Reflect_Spline();

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // write new data
      S.OutputTecplot(out());

      // check
      RunRegressionTest("Reflect_Spline()", CurrentFile, MasterFile, 5.0e-9, 5.0e-9);

    } else {
      Open_Output_File(MasterFile);

      // write data
      S.OutputTecplot(out());
    }
  }

  /* Test 36: */
  template<>
  template<>
  void Spline2D_HO_object::test<36>()
  {
    set_test_name("Reverse_Spline()");
    RunRegression = ON;
    InitializeSpline(S,10);

    MasterFile = "Reverse_Spline.dat";
    CurrentFile ="Current_Reverse_Spline.dat";

    // reverse
    S.Reverse_Spline();

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // write new data
      S.OutputTecplot(out());

      // check
      RunRegressionTest("Reverse_Spline()", CurrentFile, MasterFile, 5.0e-9, 5.0e-9);

    } else {
      Open_Output_File(MasterFile);

      // write data
      S.OutputTecplot(out());
    }
  }

  /* Test 37: */
  template<>
  template<>
  void Spline2D_HO_object::test<37>()
  {
    set_test_name("Create_Spline_Line()");
    RunRegression = ON;

    // Accepted tolerance
    tol = 5.0e-9;

    MasterFile = "Create_Spline_Line.dat";
    CurrentFile ="Current_Create_Spline_Line.dat";

    // Create line
    S.Create_Spline_Line(Vector2D(2.3,4.0),Vector2D(6.0,-5.0),3);

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // write new data
      S.OutputTecplot(out());

      // check
      RunRegressionTest("Create_Spline_Line()", CurrentFile, MasterFile, tol);

    } else {
      Open_Output_File(MasterFile);

      // write data
      S.OutputTecplot(out());
    }
  }

  /* Test 38: */
  template<>
  template<>
  void Spline2D_HO_object::test<38>()
  {
    set_test_name("Concatenate_Splines()");
    InitializeSpline(S,10);
    InitializeSpline(S1,5);

    Spline2D_HO Curve2, Copy(S), Copy1(S1);

    Curve = S + S1;
    Curve2 = Concatenate_Splines(S,S1);

    // ensure that the initial splines weren't modified
    ensure_equals("S", S, Copy);
    ensure_equals("S1", S1, Copy1);

    // ensure that the concatenated spline is correct
    ensure_equals("Concatenated",Curve2,Curve);
  }

  /* Test 39: */
  template<>
  template<>
  void Spline2D_HO_object::test<39>()
  {
    set_test_name("Create_Spline_Circular_Arc()");
    RunRegression = ON;

    MasterFile = "Create_Spline_Circular_Arc.dat";
    CurrentFile ="Current_Create_Spline_Circular_Arc.dat";

    // Create arc
    S.Create_Spline_Circular_Arc(Vector2D(2.3,4.0),1,5,120,100);

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // write new data
      S.OutputTecplot(out());

      // check
      RunRegressionTest("Create_Spline_Circular_Arc()", CurrentFile, MasterFile, 5.0e-9);

    } else {
      Open_Output_File(MasterFile);

      // write data
      S.OutputTecplot(out());
    }
  }

  /* Test 40: */
  template<>
  template<>
  void Spline2D_HO_object::test<40>()
  {
    set_test_name("Create_Spline_Ellipsoidal_Arc()");
    RunRegression = ON;

    MasterFile = "Create_Spline_Ellipsoidal_Arc.dat";
    CurrentFile ="Current_Create_Spline_Ellipsoidal_Arc.dat";

    // Create arc
    S.Create_Spline_Ellipsoidal_Arc(Vector2D(2.3,4.0),1,2,5,120,100);

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // write new data
      S.OutputTecplot(out());

      // check
      RunRegressionTest("Create_Spline_Ellipsoidal_Arc()", CurrentFile, MasterFile, 5.0e-9);

    } else {
      Open_Output_File(MasterFile);

      // write data
      S.OutputTecplot(out());
    }
  }

  /* Test 41: */
  template<>
  template<>
  void Spline2D_HO_object::test<41>()
  {
    set_test_name("Create_Spline_Sinusoidal_Line()");
    RunRegression = ON;

    MasterFile = "Create_Spline_Sinusoidal_Line.dat";
    CurrentFile ="Current_Create_Spline_Sinusoidal_Line.dat";

    // Create arc
    S.Create_Spline_Sinusoidal_Line(Vector2D(2.3,4.0),1,2,0,270,100);

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // write new data
      S.OutputTecplot(out());

      // check
      RunRegressionTest("Create_Spline_Sinusoidal_Line()", CurrentFile, MasterFile, 5.0e-9);

    } else {
      Open_Output_File(MasterFile);

      // write data
      S.OutputTecplot(out());
    }
  }

  /* Test 42: */
  template<>
  template<>
  void Spline2D_HO_object::test<42>()
  {
    set_test_name("Create_Spline_NACA_Aerofoil()");
    RunRegression = ON;

    MasterFile = "Create_Spline_NACA_Aerofoil.dat";
    CurrentFile ="Current_Create_Spline_NACA_Aerofoil.dat";

    // Create arc
    S.Create_Spline_NACA_Aerofoil("23015",10,0,100);

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // write new data
      S.OutputTecplot(out());

      // check
      RunRegressionTest("Create_Spline_NACA_Aerofoil()", CurrentFile, MasterFile, 5.0e-9);

    } else {
      Open_Output_File(MasterFile);

      // write data
      S.OutputTecplot(out());
    }
  }

  /* Test 43: */
  template<>
  template<>
  void Spline2D_HO_object::test<43>()
  {
    set_test_name("Create_Spline_Rectangle()");
    RunRegression = ON;

    MasterFile = "Create_Spline_Rectangle.dat";
    CurrentFile ="Current_Create_Spline_Rectangle.dat";

    // Create arc
    S.Create_Spline_Rectangle(Vector2D(0.0,0.0),3.4,5.3);

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // write new data
      S.OutputTecplot(out());

      // check
      RunRegressionTest("Create_Spline_Rectangle()", CurrentFile, MasterFile, 5.0e-9);

    } else {
      Open_Output_File(MasterFile);

      // write data
      S.OutputTecplot(out());
    }
  }

  /* Test 44: */
  template<>
  template<>
  void Spline2D_HO_object::test<44>()
  {
    set_test_name("Create_Spline_Ringleb_Flow()");
    RunRegression = ON;

    MasterFile = "Create_Spline_Ringleb_Flow_Tecplot.dat";
    CurrentFile ="Current_Create_Spline_Ringleb_Flow_Tecplot.dat";

    // Create spline
    S.Create_Spline_Ringleb_Flow();

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // write new data
      S.OutputTecplot(out());

      // check
      RunRegressionTest("Create_Spline_Ringleb_Flow()", CurrentFile, MasterFile, 5.0e-9);

    } else {
      Open_Output_File(MasterFile);

      // write data
      S.OutputTecplot(out());
    }
  }

  /* Test 45: */
  template<>
  template<>
  void Spline2D_HO_object::test<45>()
  {
    set_test_name("Create_Spline_Ringleb_Flow()");
    RunRegression = ON;

    MasterFile = "Create_Spline_Ringleb_Flow_Output.dat";
    CurrentFile ="Current_Create_Spline_Ringleb_Flow_Output.dat";

    // Create spline
    S.Create_Spline_Ringleb_Flow();

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // write new data
      out() << S;

      // check
      RunRegressionTest("Create_Spline_Ringleb_Flow()", CurrentFile, MasterFile, 5.0e-9);

    } else {
      Open_Output_File(MasterFile);

      // write data
      out() << S;
    }
  }


  /* Test 46: */
  template<>
  template<>
  void Spline2D_HO_object::test<46>()
  {
    set_test_name("find_subinterval()");
    int il,ir;

    CurrentFile ="NACA_Aerofoil.dat";

    // Create arc
    S.Create_Spline_NACA_Aerofoil("23015",10,0,10);

    // check s=0
    S.find_subinterval(S.sp[0],il,ir);
    ensure_equals("il: s=0",il,0);
    ensure_equals("ir: s=0",ir,1);

    // check s=latest
    S.find_subinterval(S.sp[9],il,ir);
    ensure_equals("il: s=latest",il,8);
    ensure_equals("ir: s=latest",ir,9);

    // check s=middle
    S.find_subinterval(S.sp[5],il,ir);
    ensure_equals("il: s=middle",il,4);
    ensure_equals("ir: s=middle",ir,5);

    // check s=middle
    S.find_subinterval(9.9,il,ir);
    ensure_equals("il: s=9.9",il,4);
    ensure_equals("ir: s=9.9",ir,5);

    // check s=-0.0001
    S.find_subinterval(-0.0001,il,ir);
    ensure_equals("il: s=-0.0001",il,0);
    ensure_equals("ir: s=-0.0001",ir,0);

    // check s=20.27546326378272
    S.find_subinterval(20.27546326378272,il,ir);
    ensure_equals("il: s=20.27546326378272",il,0);
    ensure_equals("ir: s=20.27546326378272",ir,0);

  }

  /* Test 47: */
  template<>
  template<>
  void Spline2D_HO_object::test<47>()
  {
    set_test_name("Vector2D_LL constructor");

    // check
    ensure_equals("next", Vlist.next==NULL, true);
    ensure_equals("prior", Vlist.prior==NULL, true);
    ensure_equals("np", Vlist.np, 0);
  }

  /* Test 48: */
  template<>
  template<>
  void Spline2D_HO_object::test<48>()
  {
    set_test_name("LinkedList<Vector2D>");

    Vector2D V1(0.0,1.0), V2(2.3,4.5), V3(4.3,5.6);

    Vlist.add(V1);
    Vlist.add(V2);
    Vlist.add(V3);

    // check
    ensure_equals("np", Vlist.np, 3);
    ensure_equals("V1", Vlist[0], V1);
    ensure_equals("V2", Vlist[1], V2);
    ensure_equals("V3", Vlist[2], V3);
    ensure_equals("V-1", Vlist[-1], V1);

  }

  /* Test 49: */
  template<>
  template<>
  void Spline2D_HO_object::test<49>()
  {
    set_test_name("getS()");
    Msg = new char [100];
    string msg;

    // Create spline
    S.Create_Spline_Ringleb_Flow();

    for (int i=0; i<=S.np-1; ++i){
      ostm() << i;
      msg = "sp[" + ostm().str() + "]"; 
      strcpy(Msg,msg.c_str());
      ensure_distance(Msg,getS(S.Xp[i],S),S.sp[i],tol*S.sp[i]+tol);
      ostmClear();		// clear the string stream
    }
  }

  /* Test 51: */
  template<>
  template<>
  void Spline2D_HO_object::test<51>()
  {
    set_test_name("1) getS()");
 
    Radius = 100.0;
    double PosTol = 1.0e-15;	// position deviation
    double Y = 1.2247e-14;
    double X = -Radius;
    double path;
    tol = 1.0e-13;
    Vector2D ReturnedPoint;

    // tested points
    Vector2D Point[] = { Vector2D(X       ,Y       ),
			 Vector2D(X+PosTol,Y       ),
			 Vector2D(X-PosTol,Y       ),
			 Vector2D(X       ,Y+PosTol),
			 Vector2D(X       ,Y-PosTol),
			 Vector2D(X+PosTol,Y+PosTol),
			 Vector2D(X+PosTol,Y-PosTol),
			 Vector2D(X-PosTol,Y+PosTol),
			 Vector2D(X-PosTol,Y-PosTol) };

    char *Msg[] = {"X = -Radius, Y = 0.0",
		   "+ X",
		   "- X",
		   "+ Y",
		   "- Y",
		   "+ X, +Y",
		   "+ X, -Y",
		   "- X, +Y",
		   "- X, -Y" };

    TotalPoints = 3;	// there are limitations for the rest of the tests

    // Create splines
    Circle1.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,180.0,250);

    for (int i = 0; i<TotalPoints; ++i ){
      path = getS(Point[i],Circle1);
      ReturnedPoint = Circle1.Spline(path);
      ensure_distance(Msg[i],ReturnedPoint.x,Point[i].x,tol*fabs(Point[i].x) + tol);
      ensure_distance(Msg[i],ReturnedPoint.y,Point[i].y,tol*fabs(Point[i].y) + tol);
    }

  }

  /* Test 52: */
  template<>
  template<>
  void Spline2D_HO_object::test<52>()
  {
    set_test_name("2) getS()");
 
    Radius = 100.0;
    double PosTol = 1.0e-15;	// position deviation
    double Y = Radius;
    double X = 0.0;
    double path;
    tol = 1.0e-10;
    Vector2D ReturnedPoint;

    // tested points
    Vector2D Point[] = { Vector2D(X       ,Y       ),
			 Vector2D(X+PosTol,Y       ),
			 Vector2D(X-PosTol,Y       ),
			 Vector2D(X       ,Y+PosTol),
			 Vector2D(X       ,Y-PosTol),
			 Vector2D(X+PosTol,Y+PosTol),
			 Vector2D(X+PosTol,Y-PosTol),
			 Vector2D(X-PosTol,Y+PosTol),
			 Vector2D(X-PosTol,Y-PosTol) };

    char *Msg[] = {"X = 0.0, Y= Radius",
		   "+ X",
		   "- X",
		   "+ Y",
		   "- Y",
		   "+ X, +Y",
		   "+ X, -Y",
		   "- X, +Y",
		   "- X, -Y" };

    TotalPoints = 9;

    // Create splines
    Circle1.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,180.0,400);

    for (int i = 0; i<TotalPoints; ++i ){
      path = Circle1.getS(Point[i]);
      ReturnedPoint = Circle1.Spline(path);
      ensure_distance(Msg[i],ReturnedPoint,Point[i],tol*Point[i] + Vector2D(tol));
    }

  }

  /* Test 53: */
  template<>
  template<>
  void Spline2D_HO_object::test<53>()
  {
    set_test_name("3) getS()");

    Radius = 1.0;
    theta = 30*degree;
    double PosTol = 1.0e-15;	// position deviation
    double Y = Radius*sin(theta);
    double X = Radius*cos(theta);
    double path;
    tol = 1.0e-10;
    Vector2D ReturnedPoint;

    // tested points
    Vector2D Point[] = { Vector2D(X       ,Y       ),
			 Vector2D(X+PosTol,Y       ),
			 Vector2D(X-PosTol,Y       ),
			 Vector2D(X       ,Y+PosTol),
			 Vector2D(X       ,Y-PosTol),
			 Vector2D(X+PosTol,Y+PosTol),
			 Vector2D(X+PosTol,Y-PosTol),
			 Vector2D(X-PosTol,Y+PosTol),
			 Vector2D(X-PosTol,Y-PosTol) };

    char *Msg[] = {"X = R*cos(30), Y= R*sin(30)",
		   "+ X",
		   "- X",
		   "+ Y",
		   "- Y",
		   "+ X, +Y",
		   "+ X, -Y",
		   "- X, +Y",
		   "- X, -Y" };

    TotalPoints = 9;

    // Create splines
    Circle1.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,180.0,300);

    for (int i = 0; i<TotalPoints; ++i ){
      path = Circle1.getS(Point[i]);
      ReturnedPoint = Circle1.Spline(path);
      ensure_distance(Msg[i],ReturnedPoint,Point[i],tol*Point[i] + Vector2D(tol));
    }

  }

  /* Test 54: */
  template<>
  template<>
  void Spline2D_HO_object::test<54>()
  {
    set_test_name("4) getS()");

    Radius = 100.0;
    theta = 60*degree;
    double PosTol = 1.0e-15;	// position deviation
    double Y = Radius*sin(theta);
    double X = Radius*cos(theta);
    double path;
    tol = 1.0e-13;
    Vector2D ReturnedPoint;

    // tested points
    Vector2D Point[] = { Vector2D(X       ,Y       ),
			 Vector2D(X+PosTol,Y       ),
			 Vector2D(X-PosTol,Y       ),
			 Vector2D(X       ,Y+PosTol),
			 Vector2D(X       ,Y-PosTol),
			 Vector2D(X+PosTol,Y+PosTol),
			 Vector2D(X+PosTol,Y-PosTol),
			 Vector2D(X-PosTol,Y+PosTol),
			 Vector2D(X-PosTol,Y-PosTol) };

    char *Msg[] = {"X = R*cos(60), Y= R*sin(60)",
		   "+ X",
		   "- X",
		   "+ Y",
		   "- Y",
		   "+ X, +Y",
		   "+ X, -Y",
		   "- X, +Y",
		   "- X, -Y" };

    TotalPoints = 9;

    // Create splines
    Circle1.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,180.0,250);

    for (int i = 0; i<TotalPoints; ++i ){
      path = Circle1.getS(Point[i]);
      ReturnedPoint = Circle1.Spline(path);
      ensure_distance(Msg[i],ReturnedPoint,Point[i],tol*Point[i] + Vector2D(tol));
    }

  }

  /* Test 55: */
  template<>
  template<>
  void Spline2D_HO_object::test<55>()
  {
    set_test_name("5) getS()");

    Radius = 100.0;
    theta = 45*degree;
    double PosTol = 1.0e-15;	// position deviation
    double Y = Radius*sin(theta);
    double X = Radius*cos(theta);
    double path;
    tol = 1.0e-10;
    Vector2D ReturnedPoint;

    // tested points
    Vector2D Point[] = { Vector2D(X       ,Y       ),
			 Vector2D(X+PosTol,Y       ),
			 Vector2D(X-PosTol,Y       ),
			 Vector2D(X       ,Y+PosTol),
			 Vector2D(X       ,Y-PosTol),
			 Vector2D(X+PosTol,Y+PosTol),
			 Vector2D(X+PosTol,Y-PosTol),
			 Vector2D(X-PosTol,Y+PosTol),
			 Vector2D(X-PosTol,Y-PosTol) };

    char *Msg[] = {"X = R*cos(45), Y= R*sin(45)",
		   "+ X",
		   "- X",
		   "+ Y",
		   "- Y",
		   "+ X, +Y",
		   "+ X, -Y",
		   "- X, +Y",
		   "- X, -Y" };

    TotalPoints = 9;

    // Create splines
    Circle1.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,180.0,400);

    for (int i = 0; i<TotalPoints; ++i ){
      path = Circle1.getS(Point[i]);
      ReturnedPoint = Circle1.Spline(path);
      ensure_distance(Msg[i],ReturnedPoint,Point[i],tol*Point[i] + Vector2D(tol));
    }

  }

  /* Test 56: */
  template<>
  template<>
  void Spline2D_HO_object::test<56>()
  {
    set_test_name("6) getS()");
 
    Radius = 100.0;
    theta = 14.56*degree;
    double PosTol = 1.0e-14;	// position deviation
    double Y = Radius*sin(theta);
    double X = Radius*cos(theta);
    double path;
    Vector2D ReturnedPoint;

    tol = 1.0e-10;

    // tested points
    Vector2D Point[] = { Vector2D(X       ,Y       ),
			 Vector2D(X+PosTol,Y       ),
			 Vector2D(X-PosTol,Y       ),
			 Vector2D(X       ,Y+PosTol),
			 Vector2D(X       ,Y-PosTol),
			 Vector2D(X+PosTol,Y+PosTol),
			 Vector2D(X+PosTol,Y-PosTol),
			 Vector2D(X-PosTol,Y+PosTol),
			 Vector2D(X-PosTol,Y-PosTol) };

    char *Msg[] = {"X = R*cos(14.56), Y=R*sin(14.56) ",
		   "+ X",
		   "- X",
		   "+ Y",
		   "- Y",
		   "+ X, +Y",
		   "+ X, -Y",
		   "- X, +Y",
		   "- X, -Y" };

    TotalPoints = 9;

    // Create splines
    Circle1.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,180.0,250);

    for (int i = 0; i<TotalPoints; ++i ){
      path = Circle1.getS(Point[i]);
      ReturnedPoint = Circle1.Spline(path);
      ensure_distance(Msg[i],ReturnedPoint,Point[i],tol*Point[i] + Vector2D(tol));
    }
  }

  /* Test 57: */
  template<>
  template<>
  void Spline2D_HO_object::test<57>()
  {
    set_test_name("7) getS()");
 
    Radius = 100.0;
    theta = 140.56*degree;
    double PosTol = 1.0e-14;	// position deviation
    double Y = Radius*sin(theta);
    double X = Radius*cos(theta);
    double path;
    Vector2D ReturnedPoint;

    tol = 1.0e-10;

    // tested points
    Vector2D Point[] = { Vector2D(X       ,Y       ),
			 Vector2D(X+PosTol,Y       ),
			 Vector2D(X-PosTol,Y       ),
			 Vector2D(X       ,Y+PosTol),
			 Vector2D(X       ,Y-PosTol),
			 Vector2D(X+PosTol,Y+PosTol),
			 Vector2D(X+PosTol,Y-PosTol),
			 Vector2D(X-PosTol,Y+PosTol),
			 Vector2D(X-PosTol,Y-PosTol) };

    char *Msg[] = {"X = R*cos(140.56), Y=R*sin(140.56) ",
		   "+ X",
		   "- X",
		   "+ Y",
		   "- Y",
		   "+ X, +Y",
		   "+ X, -Y",
		   "- X, +Y",
		   "- X, -Y" };

    TotalPoints = 9;

    // Create splines
    Circle1.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,180.0,430);

    for (int i = 0; i<TotalPoints; ++i ){
      path = Circle1.getS(Point[i]);
      ReturnedPoint = Circle1.Spline(path);
      ensure_distance(Msg[i],ReturnedPoint,Point[i],AcceptedError(Point[i],tol) );
    }
  }


  /* Test 58: */
  template<>
  template<>
  void Spline2D_HO_object::test<58>()
  {
    set_test_name("nSpline(), Circular_Arc, 45 degrees");

    // Accepted tolerance
    tol = 1.0e-10;
    // Create the circular arc
    NumSplinePoints = 200;
    V1 = Vector2D(3.0,2.0);	// centroid
    S.Create_Spline_Circular_Arc(V1,Radius,Angle1,Angle2,NumSplinePoints);

    // Analytic results
    AnalyticNormal = Vector2D(0.5*sqrt(2.0),0.5*sqrt(2.0));

    // The point on the spline where the normal is computed
    Vector2D Point = Vector2D(V1.x + Radius*0.5*sqrt(2.0),V2.y + Radius*0.5*sqrt(2.0));

    NumericNormal = S.nSpline(Point);
    ensure_distance("Normal at 45 degrees",NumericNormal,AnalyticNormal,AcceptedError(AnalyticNormal,tol));
  }

  /* Test 59: */
  template<>
  template<>
  void Spline2D_HO_object::test<59>()
  {
    set_test_name("nSpline(), Circular_Arc, 30 degrees");

    // Accepted tolerance
    tol = 1.0e-7;
    // Create the circular arc
    NumSplinePoints = 200;
    V1 = Vector2D(3.0,2.0);	// centroid
    S.Create_Spline_Circular_Arc(V1,Radius,Angle1,Angle2,NumSplinePoints);

    // Analytic results
    AnalyticNormal = Vector2D(0.5*sqrt(3.0),0.5);

    // The point on the spline where the normal is computed
    Vector2D Point = Vector2D(V1.x + Radius*0.5*sqrt(3.0),V2.y + Radius*0.5*sqrt(3.0));

    // check
    NumericNormal = S.nSpline(Point);
    ensure_distance("Normal at 30 degrees",NumericNormal,AnalyticNormal,AcceptedError(AnalyticNormal,tol));
  }

  /* Test 60: */
  template<>
  template<>
  void Spline2D_HO_object::test<60>()
  {
    set_test_name("nSpline(), Circular_Arc, 60 degrees");

    // Accepted tolerance
    tol = 1.0e-7;
    // Create the circular arc
    NumSplinePoints = 200;
    V1 = Vector2D(3.0,2.0);	// centroid
    S.Create_Spline_Circular_Arc(V1,Radius,Angle1,Angle2,NumSplinePoints);

    // Analytic results
    AnalyticNormal = Vector2D(0.5,0.5*sqrt(3.0));

    // The point on the spline where the normal is computed
    Vector2D Point = Vector2D(V1.x + Radius*0.5 ,V2.y + Radius*0.5*sqrt(3.0));

    // check
    NumericNormal = S.nSpline(Point);
    ensure_distance("Normal at 60 degrees",NumericNormal,AnalyticNormal,AcceptedError(AnalyticNormal,tol));
  }

  /* Test 61: */
  template<>
  template<>
  void Spline2D_HO_object::test<61>()
  {
    set_test_name("nSpline(), Circular_Arc, 90 degrees");

    // Accepted tolerance
    tol = 1.0e-6;

    // Create the circular arc
    NumSplinePoints = 200;
    V1 = Vector2D(3.0,2.0);	// centroid
    S.Create_Spline_Circular_Arc(V1,Radius,Angle1,Angle2,NumSplinePoints);

    // Analytic results
    AnalyticNormal = Vector2D(0.0,1.0);

    // The point on the spline where the normal is computed
    Vector2D Point = Vector2D(V1.x + Radius*0.0 ,V2.y + Radius*1.0);

    // check
    NumericNormal = S.nSpline(Point);
    ensure_distance("Normal at 90 degrees",NumericNormal,AnalyticNormal,AcceptedError(AnalyticNormal,tol));
  }

  /* Test 62: */
  template<>
  template<>
  void Spline2D_HO_object::test<62>()
  {
    set_test_name("nSpline(), Circular_Arc, 0 degrees");

    // Accepted tolerance
    tol = 1.0e-6;

    // Create the circular arc
    NumSplinePoints = 200;
    V1 = Vector2D(3.0,2.0);	// centroid
    S.Create_Spline_Circular_Arc(V1,Radius,Angle1,Angle2,NumSplinePoints);

    // Analytic results
    AnalyticNormal = Vector2D(1.0,0.0);

    // The point on the spline where the normal is computed
    Vector2D Point = Vector2D(V1.x + Radius*1.0 ,V2.y + Radius*0.0);

    // check
    NumericNormal = S.nSpline(Point);
    ensure_distance("Normal at 0 degrees",NumericNormal,AnalyticNormal,AcceptedError(AnalyticNormal,tol));
  }

  /* Test 63: */
  template<>
  template<>
  void Spline2D_HO_object::test<63>()
  {
    set_test_name("nSpline(), Rectangle 1");

    // Accepted tolerance
    tol = 1.0e-10;

    // Create the rectangle
    V1 = Vector2D(0.5,1.3);	// centroid
    S.Create_Spline_Rectangle(V1,Length,Width);

    // Analytic results
    AnalyticNormal = Vector2D(0.0,-1.0);

    // The point on the spline where the normal is computed
    Vector2D Point = Vector2D(0.5,-0.2);

    // check
    NumericNormal = S.nSpline(Point);
    ensure_distance("Rectangle 1.",NumericNormal,AnalyticNormal,AcceptedError(AnalyticNormal,tol));
  }

  /* Test 64: */
  template<>
  template<>
  void Spline2D_HO_object::test<64>()
  {
    set_test_name("nSpline(), Rectangle 2");

    // Accepted tolerance
    tol = 1.0e-10;

    // Create the rectangle
    V1 = Vector2D(0.5,1.3);	// centroid
    S.Create_Spline_Rectangle(V1,Length,Width);

    // Analytic results
    AnalyticNormal = Vector2D(0.0,-1.0);

    // The point on the spline where the normal is computed
    Vector2D Point = Vector2D(-0.7,-0.2);

    // check
    NumericNormal = S.nSpline(Point);
    ensure_distance("Rectangle 2.",NumericNormal,AnalyticNormal,AcceptedError(AnalyticNormal,tol));
  }

  /* Test 65: */
  template<>
  template<>
  void Spline2D_HO_object::test<65>()
  {
    set_test_name("nSpline(), Rectangle 3");

    // Accepted tolerance
    tol = 1.0e-10;

    // Create the rectangle
    V1 = Vector2D(0.5,1.3);	// centroid
    S.Create_Spline_Rectangle(V1,Length,Width);

    // Analytic results
    AnalyticNormal = Vector2D(-1.0,0.0);

    // The point on the spline where the normal is computed
    Vector2D Point = Vector2D(-0.7,0.0);

    // check
    NumericNormal = S.nSpline(Point);
    ensure_distance("Rectangle 3.",NumericNormal,AnalyticNormal,AcceptedError(AnalyticNormal,tol));
  }

   /* Test 66: */
  template<>
  template<>
  void Spline2D_HO_object::test<66>()
  {
    set_test_name("nSpline(), Rectangle 4");

    // Accepted tolerance
    tol = 1.0e-10;

    // Create the rectangle
    V1 = Vector2D(0.5,1.3);	// centroid
    S.Create_Spline_Rectangle(V1,Length,Width);

    // Analytic results
    AnalyticNormal = Vector2D(0.0,1.0);

    // The point on the spline where the normal is computed
    Vector2D Point = Vector2D(0.5,2.8);

    // check
    NumericNormal = S.nSpline(Point);
    ensure_distance("Rectangle 4.",NumericNormal,AnalyticNormal,AcceptedError(AnalyticNormal,tol));
  }

   /* Test 67: */
  template<>
  template<>
  void Spline2D_HO_object::test<67>()
  {
    set_test_name("nSpline(), Rectangle 5");

    // Accepted tolerance
    tol = 1.0e-10;

    // Create the rectangle
    V1 = Vector2D(0.5,1.3);	// centroid
    S.Create_Spline_Rectangle(V1,Length,Width);

    // Analytic results
    AnalyticNormal = Vector2D(1.0,0.0);

    // The point on the spline where the normal is computed
    Vector2D Point = Vector2D(1.7,2.0);

    // check
    NumericNormal = S.nSpline(Point);
    ensure_distance("Rectangle 5.",NumericNormal,AnalyticNormal,AcceptedError(AnalyticNormal,tol));
  }

   /* Test 68: */
  template<>
  template<>
  void Spline2D_HO_object::test<68>()
  {
    set_test_name("nSpline(), Rectangle 6");

    // Accepted tolerance
    tol = 1.0e-10;

    // Create the rectangle
    V1 = Vector2D(0.5,1.3);	// centroid
    S.Create_Spline_Rectangle(V1,Length,Width);

    // Analytic results
    AnalyticNormal = Vector2D(0.0,-1.0);

    // The point on the spline where the normal is computed
    Vector2D Point = Vector2D(1.7,-0.2);

    // check
    NumericNormal = S.nSpline(Point);
    ensure_distance("Rectangle 6.",NumericNormal,AnalyticNormal,AcceptedError(AnalyticNormal,tol));
  }

   /* Test 69: */
  template<>
  template<>
  void Spline2D_HO_object::test<69>()
  {
    set_test_name("nSpline(), Rectangle 7");

    // Accepted tolerance
    tol = 1.0e-10;

    // Create the rectangle
    V1 = Vector2D(0.5,1.3);	// centroid
    S.Create_Spline_Rectangle(V1,Length,Width);

    // Analytic results
    AnalyticNormal = Vector2D(-1.0,0.0);

    // The point on the spline where the normal is computed
    Vector2D Point = Vector2D(-0.7,-0.1999);

    // check
    NumericNormal = S.nSpline(Point);
    ensure_distance("Rectangle 7.",NumericNormal,AnalyticNormal,AcceptedError(AnalyticNormal,tol));
  }


   /* Test 70: */
  template<>
  template<>
  void Spline2D_HO_object::test<70>()
  {
    set_test_name("tSpline(), Rectangle 1");

    // Create the rectangle
    V1 = Vector2D(0.5,1.3);	// centroid
    S.Create_Spline_Rectangle(V1,Length,Width);

    // Analytic results
    AnalyticTangent = Vector2D(1.0,0.0);

    // Compute the numerical tangent
    double s(1.2);
    NumericTangent = S.tSpline(S.Spline(s), LINEAR);

    ensure_distance("Tangent", NumericTangent, AnalyticTangent, AcceptedError(AnalyticTangent,tol));
  }

   /* Test 71: */
  template<>
  template<>
  void Spline2D_HO_object::test<71>()
  {
    set_test_name("tSpline(), Rectangle 2, Corner");

    // Create the rectangle
    V1 = Vector2D(0.5,1.3);	// centroid
    S.Create_Spline_Rectangle(V1,Length,Width);

    // Analytic results
    AnalyticTangent = Vector2D(1.0,0.0);

    // Compute the numerical tangent
    NumericTangent = S.tSpline(S.Xp[0],CUBIC);

    ensure_distance("Tangent", NumericTangent, AnalyticTangent, AcceptedError(AnalyticTangent,tol));
  }

   /* Test 72: */
  template<>
  template<>
  void Spline2D_HO_object::test<72>()
  {
    set_test_name("tSpline(), Circle , 0 degrees, 3rd-Order polynom");

    // Accepted tolerance
    tol = 2.0e-7;

    // Create the circular arc
    NumSplinePoints = 200;
    V1 = Vector2D(3.0,2.0);	// centroid
    S.Create_Spline_Circular_Arc(V1,Radius,Angle1,Angle2,NumSplinePoints);

    // Analytic results
    AnalyticTangent = Vector2D(0.0,1.0);

    // The point on the spline where the normal is computed
    Vector2D Point = Vector2D(V1.x + Radius*1.0 ,V2.y + Radius*0.0);

    // Compute the numerical tangent
    NumericTangent = S.tSpline(Point,CUBIC);

    ensure_distance("Tangent", NumericTangent, AnalyticTangent, AcceptedError(AnalyticTangent,tol));
  }

   /* Test 73: */
  template<>
  template<>
  void Spline2D_HO_object::test<73>()
  {
    set_test_name("tSpline(), Circle, 60 degrees, 2nd-Order polynom");

    // Accepted tolerance
    tol = 1.0e-7;

    // Create the circular arc
    NumSplinePoints = 200;
    V1 = Vector2D(3.0,2.0);	// centroid
    S.Create_Spline_Circular_Arc(V1,Radius,Angle1,Angle2,NumSplinePoints);

    // Analytic results
    AnalyticTangent = Vector2D(-0.5*sqrt(3.0),0.5);

    // The point on the spline where the normal is computed
    Vector2D Point = Vector2D(V1.x + Radius*0.5 ,V2.y + Radius*0.5*sqrt(3.0));

    // Compute the numerical tangent
    NumericTangent = S.tSpline(Point,QUADRATIC);

    ensure_distance("Tangent", NumericTangent, AnalyticTangent, AcceptedError(AnalyticTangent,tol));
  }

   /* Test 74: */
  template<>
  template<>
  void Spline2D_HO_object::test<74>()
  {
    set_test_name("tSpline(), Circle, 60 degrees, 3rd-Order polynom");

    // Accepted tolerance
    tol = 1.0e-7;

    // Create the circular arc
    NumSplinePoints = 200;
    V1 = Vector2D(3.0,2.0);	// centroid
    S.Create_Spline_Circular_Arc(V1,Radius,Angle1,Angle2,NumSplinePoints);

    // Analytic results
    AnalyticTangent = Vector2D(-0.5*sqrt(3.0),0.5);

    // The point on the spline where the normal is computed
    Vector2D Point = Vector2D(V1.x + Radius*0.5 ,V2.y + Radius*0.5*sqrt(3.0));

    // Compute the numerical tangent
    NumericTangent = S.tSpline(Point,CUBIC);

    ensure_distance("Tangent", NumericTangent, AnalyticTangent, AcceptedError(AnalyticTangent,tol));
  }

  /* Test 75: */
  template<>
  template<>
  void Spline2D_HO_object::test<75>()
  {
    set_test_name("tSpline(), Circular_Arc, 90 degrees");

    // Accepted tolerance
    tol = 1.0e-6;

    // Create the circular arc
    NumSplinePoints = 200;
    V1 = Vector2D(3.0,2.0);	// centroid
    S.Create_Spline_Circular_Arc(V1,Radius,Angle1,Angle2,NumSplinePoints);

    // Analytic results
    AnalyticTangent = Vector2D(-1.0,0.0);

    // The point on the spline where the normal is computed
    Vector2D Point = Vector2D(V1.x + Radius*0.0 ,V2.y + Radius*1.0);

    // check
    NumericTangent = S.tSpline(Point,CUBIC);
    ensure_distance("Tangent at 90 degrees",NumericTangent,AnalyticTangent,AcceptedError(AnalyticTangent,tol));
  }

  /* Test 76: */
  template<>
  template<>
  void Spline2D_HO_object::test<76>()
  {
    set_test_name("nSpline(), Circular_Arc, 120 degrees");

    // Accepted tolerance
    tol = 1.0e-7;
    // Create the circular arc
    NumSplinePoints = 200;
    V1 = Vector2D(3.0,2.0);	// centroid
    Angle1 = 90;
    Angle2 = 180;
    S.Create_Spline_Circular_Arc(V1,Radius,Angle1,Angle2,NumSplinePoints);

    // Analytic results
    AnalyticNormal = Vector2D(-0.5,0.5*sqrt(3.0));

    // The point on the spline where the normal is computed
    Vector2D Point = Vector2D(V1.x - Radius*0.5 ,V2.y + Radius*0.5*sqrt(3.0));

    // check
    NumericNormal = S.nSpline(Point);
    ensure_distance("Normal at 120 degrees",NumericNormal,AnalyticNormal,AcceptedError(AnalyticNormal,tol));

  }

  /* Test 77: */
  template<>
  template<>
  void Spline2D_HO_object::test<77>()
  {
    set_test_name("Two solid bodies");

    // Create the circular arc
    NumSplinePoints = 200;
    V1 = Vector2D(3.0,2.0);	// centroid
    Angle1 = 90;
    Angle2 = 180;
    S.Create_Spline_Circular_Arc(V1,Radius,Angle1,Angle2,NumSplinePoints);
    // Make S solid boundary of a new body
    S.makeSplineSolidBoundary();

    // Create Copy
    Spline2D_HO Copy(S);
    // Make Copy solid boundary of a new body
    Copy.makeSplineSolidBoundary();
    
    // == check
    ensure_equals("total number of solid bodies", Spline2D_HO::NumberOfSolidBodies(), 2);
    ensure_equals("bodyID for S", S.getBodyID(), 1);
    ensure_equals("bodyID for Copy", Copy.getBodyID(), 2);

    // Change body for spline 1
    S.makeSplineSolidBoundary(2);

    ensure_equals("bodyID for S", S.getBodyID(), 2);

  }

  /* Test 78: */
  template<>
  template<>
  void Spline2D_HO_object::test<78>()
  {
    set_test_name("Concatenate two solid bodies");

    // Create the circular arc
    NumSplinePoints = 200;
    V1 = Vector2D(3.0,2.0);	// centroid
    Angle1 = 90;
    Angle2 = 180;
    S.Create_Spline_Circular_Arc(V1,Radius,Angle1,Angle2,NumSplinePoints);
    // Make S solid boundary of a new body
    S.makeSplineSolidBoundary();

    // Create Copy
    Spline2D_HO Copy(S);
    // Make Copy solid boundary of a new body
    Copy.makeSplineSolidBoundary();
    // Translate Copy
    Copy.Translate_Spline(Vector2D(3.5,4.5));
    
    // == check
    ensure_equals("total number of solid bodies", Spline2D_HO::NumberOfSolidBodies(), 2);
    ensure_equals("bodyID for S", S.getBodyID(), 1);
    ensure_equals("bodyID for Copy", Copy.getBodyID(), 2);

    // Concatenate splines
    Spline2D_HO MergedCurve;
    MergedCurve = S + Copy;

    // == check
    ensure_equals("total number of solid bodies", Spline2D_HO::NumberOfSolidBodies(), 2);
    ensure_equals("bodyID for MergedCurve", MergedCurve.getBodyID(), 1);
  }

}



// Test suite constructor
tut::Spline2D_HO_TestSuite Spline2D_HOTestSuite("Class:Spline2D_HO");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

