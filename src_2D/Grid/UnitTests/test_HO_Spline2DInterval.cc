/*!\file test_HO_Spline2DInterval.cc
  \brief Regression tests for class Spline2DInterval_HO. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../HO_Spline2DInterval.h"
#include "../../HighOrderReconstruction/TaylorDerivatives.h"

namespace tut
{

  // Define a test function
  inline double TestFunction(const double & x, const double & y){ return exp(y); }

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_Spline2DInterval_HO : public TestData {

    // Local variables
  public:
    typedef TaylorDerivativesContainer<TwoD,double> GeometricMoments;

    Spline2DInterval_HO SInfo, SInfo1;
    int NumSubIntervals, N_GQP;

    Spline2D_HO S, S1;
    Spline2D_HO Curve;
    Spline2D_HO Circle1, Circle2;
    Spline2D_HO RectangleCurve;

    Vector2D AnalyticNormal, NumericNormal;
    Vector2D Normal, P, P1, P2;
    Vector2D V1, V2, Centroid;

    double Radius, Angle1, Angle2;
    double theta;
    double degree;		// to transform an angle given in degrees to rad, just multiply with 'degree'
    int NumSplinePoints;

    double Length, Width;
    double AnalyticResult;

    // Null pointers
    Vector2D *NullVector;
    Spline2D_HO *NullSplinePtr;
    double *NullDouble;

    // Constructor
    Data_Spline2DInterval_HO(void);

    Vector2D SplinePoint(const int &i);
    void InitializeSpline(Spline2D_HO &S);
    void InitializeSpline(Spline2D_HO &S, const int &N);
    void PrintPathlength(Spline2D_HO &S);

  private:
    
  };

  // Constructor
  Data_Spline2DInterval_HO::Data_Spline2DInterval_HO(void){
    /* Set the global path for this test suite. 
       It's a good practice to put it in the constructor of the data class in order to be set
       automatically for each individual test. Declare it relative to the /src_2D directory,
       otherwise the framework might not find the input and output files. */
    set_test_suite_path("Grid/UnitTests");

    // set default class parameters for each unit test
    Spline2DInterval_HO::Set_Default_Parameters();

    Radius = 1.0232323e-1;
    degree = TWO*PI/360.0;
    Angle1 = 0.0;
    Angle2 = 90.0;
    NumSplinePoints = 10;
    Length = 3.0;
    Width  = 2.4;

    // initialize the Null pointers
    NullVector = NULL;
    NullSplinePtr = NULL;
    NullDouble = NULL;
  }

  // Generate spline points
  Vector2D Data_Spline2DInterval_HO::SplinePoint(const int &i){
    return Vector2D(i*0.34345, (i+2.2323)*0.3434);
  }

  // Initialize spline
  void Data_Spline2DInterval_HO::InitializeSpline(Spline2D_HO & S){
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

  // Initialize a spline with a specified number of control points
  void Data_Spline2DInterval_HO::InitializeSpline(Spline2D_HO &S, const int &N){
    S.allocate(N);
    InitializeSpline(S);
  }

  // Output the pathlength data to standard output
  void Data_Spline2DInterval_HO::PrintPathlength(Spline2D_HO &S){
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
  typedef test_group<Data_Spline2DInterval_HO> Spline2DInterval_HO_TestSuite;
  typedef Spline2DInterval_HO_TestSuite::object Spline2DInterval_HO_object;


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
  void Spline2DInterval_HO_object::test<1>()
  {

    set_test_name("Constructor");

    // check
    ensure_equals("N_SubIntervals",SInfo.NumOfSubIntervals(),0);
    ensure_equals("N_GQP",SInfo.GQPointsPerSubInterval(),0);
    ensure_equals("Total N_GQP",SInfo.NumGQPoints(),0);
    ensure_equals("IntervalLength",SInfo.IntLength(), NullDouble);
    ensure_equals("GQP",SInfo.GQPoint(),NullVector);
    ensure_equals("NormalsGQP",SInfo.NormalGQPoint(), NullVector);
    ensure_equals("Total GQP_ContourInt",SInfo.NumGQPoints_ContourIntegral(),0);
    ensure_equals("GQP_ContourInt",SInfo.GQPointContourIntegral(),NullVector);
    ensure_equals("dYdS",SInfo.dYdS(), NullDouble);
  }

   /* Test 2: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<2>()
  {
    set_test_name("allocate()");
    NumSubIntervals = 4;
    N_GQP = 2;

    // allocate
    SInfo.allocate(N_GQP,NumSubIntervals);

    // check
    ensure_equals("Number of subintervals", SInfo.NumOfSubIntervals(),NumSubIntervals);
    ensure_equals("Gauss Points per subinterval", SInfo.GQPointsPerSubInterval(), N_GQP);
    ensure_equals("Last Gauss Point", SInfo.GQPoint(NumSubIntervals*N_GQP),Vector2D(0.0));
    ensure_equals("Last Normal at GQP", SInfo.NormalGQPoint(NumSubIntervals*N_GQP),Vector2D(0.0));
    ensure_equals("Total GQP_ContourInt",SInfo.NumGQPoints_ContourIntegral(),12);
    ensure_equals("Last Gauss Point Contour Int", SInfo.GQPointContourIntegral(NumSubIntervals*3),Vector2D(0.0));
    ensure_equals("Last dYdS", SInfo.dYdS(NumSubIntervals), 0.0);
    ensure_equals("Last SubInterval length", SInfo.IntLength(NumSubIntervals), 0.0);

  }

  /* Test 3: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<3>()
  {
    set_test_name("deallocate()");

    NumSubIntervals = 4;
    N_GQP = 2;

    // allocate
    SInfo.allocate(N_GQP,NumSubIntervals);
    SInfo.deallocate();

    // check
    ensure_equals("IntervalLength",SInfo.IntLength(), NullDouble);
    ensure_equals("GQP",SInfo.GQPoint(),NullVector);
    ensure_equals("NormalsGQP",SInfo.NormalGQPoint(), NullVector);
    ensure_equals("Total GQP_ContourInt",SInfo.NumGQPoints_ContourIntegral(),0);
    ensure_equals("MomentGQP",SInfo.GQPointContourIntegral(),NullVector);
    ensure_equals("dYdS",SInfo.dYdS(), NullDouble);
    ensure_equals("Number of subintervals", SInfo.NumOfSubIntervals(),0);
    ensure_equals("Gauss Points per subinterval", SInfo.GQPointsPerSubInterval(), 0);
  }


  /* Test 4: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<4>()
  {
    set_test_name("2) allocate()");
    NumSubIntervals = 4;
    N_GQP = 2;

    // allocate
    SInfo.allocate(N_GQP,NumSubIntervals);

    NumSubIntervals = 10;
    N_GQP = 1;

    // allocate again
    SInfo.allocate(N_GQP,NumSubIntervals);

    // check
    ensure_equals("Number of subintervals", SInfo.NumOfSubIntervals(),NumSubIntervals);
    ensure_equals("Gauss Points per subinterval", SInfo.GQPointsPerSubInterval(), N_GQP);
    ensure_equals("Total GQP_ContourInt",SInfo.NumGQPoints_ContourIntegral(),3*NumSubIntervals);
    ensure_equals("Last Gauss Point", SInfo.GQPoint(NumSubIntervals*N_GQP),Vector2D(0.0));
    ensure_equals("Last Normal at GQP", SInfo.NormalGQPoint(NumSubIntervals*N_GQP),Vector2D(0.0));
    ensure_equals("Last Gauss Point Moment", SInfo.GQPointContourIntegral(NumSubIntervals*3),Vector2D(0.0));
    ensure_equals("Last dYdS", SInfo.dYdS(NumSubIntervals*N_GQP), 0.0);
    ensure_equals("Last SubInterval length", SInfo.IntLength(NumSubIntervals), 0.0);
  }


  /* Test 5: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<5>()
  {
    set_test_name("1) InitializeInterval() for a line");

    V1 = Vector2D(0.0);
    V2 = Vector2D(4.0,3.0);

    // Create Spline
    S.Create_Spline_Line(V1,V2,10);

    // Initialize interval
    SInfo.InitializeInterval(S,S.Xp[0],S.Xp[9],2);

    // check
    ensure_equals("Number of subintervals", SInfo.NumOfSubIntervals(),1);
    ensure_equals("Gauss Points per subinterval", SInfo.GQPointsPerSubInterval(), 2);

    ensure_distance("SubInterval length", SInfo.IntLength(1), 5.0, tol);

    P = Vector2D(8.452994616207484709817024e-1, 6.339745962155613532362768e-1);
    ensure_distance("First Gauss Point", SInfo.GQPoint(1), P, AcceptedError(P,tol));
    P = Vector2D(3.1547005383792515290182976e+0,2.3660254037844386467637232e+0);
    ensure_distance("Second Gauss Point", SInfo.GQPoint(2), P, AcceptedError(P,tol));

    AnalyticNormal = Vector2D(0.6,-0.8);
    ensure_distance("First Normal at GQP" , SInfo.NormalGQPoint(1), AnalyticNormal, AcceptedError(AnalyticNormal,tol));
    ensure_distance("Second Normal at GQP", SInfo.NormalGQPoint(2), AnalyticNormal, AcceptedError(AnalyticNormal,tol));

    ensure_distance("First dYdS at GQP_ContourInt" , SInfo.dYdS(1), AnalyticNormal.x, AcceptedError(AnalyticNormal.x,tol));
    ensure_distance("Second dYdS at GQP_ContourInt", SInfo.dYdS(2), AnalyticNormal.x, AcceptedError(AnalyticNormal.x,tol));
    ensure_distance("Third dYdS at GQP_ContourInt" , SInfo.dYdS(3), AnalyticNormal.x, AcceptedError(AnalyticNormal.x,tol));
  }

  /* Test 6: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<6>()
  {
    set_test_name("2) InitializeInterval() for a quater of a circle");

    double path_coord;

    tol = 1.0e-10;

    // Create Spline
    S.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,90.0,200);

    // Initialize interval
    SInfo.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],2);

    // check
    ensure_equals("Number of subintervals", SInfo.NumOfSubIntervals(),1);
    ensure_equals("Gauss Points per subinterval", SInfo.GQPointsPerSubInterval(), 2);

    double Answer =  0.5*PI*Radius;
    ensure_distance("SubInterval length", SInfo.IntLength(1), Answer, AcceptedError(Answer,tol));

    // First Gauss Point
    path_coord = 0.5*PI*Radius*0.2113248654051871177454256;
    P1 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("First Gauss Point", SInfo.GQPoint(1), P1, AcceptedError(P1,tol));

    // Second Gauss Point
    path_coord = 0.5*PI*Radius*0.7886751345948128822545744;
    P2 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("Second Gauss Point", SInfo.GQPoint(2), P2, AcceptedError(P2,tol));


    // Check Normals
    tol = 1.0e-7;		// lower tolerance is necessary
    AnalyticNormal = P1/abs(P1);
    ensure_distance("First Normal at GQP", SInfo.NormalGQPoint(1), AnalyticNormal, AcceptedError(AnalyticNormal,tol));

    AnalyticNormal = P2/abs(P2);
    ensure_distance("Second Normal at GQP", SInfo.NormalGQPoint(2), AnalyticNormal,AcceptedError(AnalyticNormal,tol));

    // First Moment Gauss Point
    path_coord = 0.5*PI*Radius*0.1127016653792583114820736;
    P1 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("First Moment Gauss Point", SInfo.GQPointContourIntegral(1), P1, AcceptedError(P1,tol));

    // Second Moment Gauss Point
    path_coord = 0.5*PI*Radius*0.5;
    P1 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("Second Moment Gauss Point", SInfo.GQPointContourIntegral(2), P1, AcceptedError(P1,tol));

    // Third Moment Gauss Point
    path_coord = 0.5*PI*Radius*0.8872983346207416885179264;
    P1 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("Third Moment Gauss Point", SInfo.GQPointContourIntegral(3), P1, AcceptedError(P1,tol));
  }

  /* Test 7: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<7>()
  {
    set_test_name("UpdateInterval()");

    double path_coord;

    tol = 1.0e-10;

    V1 = Vector2D(0.0);
    V2 = Vector2D(4.0,3.0);

    // Create Spline 1
    S.Create_Spline_Line(V1,V2,10);

    // Initialize interval
    SInfo.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],1);

    // Create Spline 2
    S.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,90.0,200);

    // Update interval (different number of GQP!!)
    SInfo.UpdateInterval(S,S.Xp[0],S.Xp[S.np-1],2);

    // check
    ensure_equals("Number of subintervals", SInfo.NumOfSubIntervals(),1);
    ensure_equals("Gauss Points per subinterval", SInfo.GQPointsPerSubInterval(), 2);

    double Answer =  0.5*PI*Radius;
    ensure_distance("SubInterval length", SInfo.IntLength(1), Answer, AcceptedError(Answer,tol));

    // First Gauss Point
    path_coord = 0.5*PI*Radius*0.2113248654051871177454256;
    P1 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("First Gauss Point", SInfo.GQPoint(1), P1, AcceptedError(P1,tol));

    // Second Gauss Point
    path_coord = 0.5*PI*Radius*0.7886751345948128822545744;
    P2 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("Second Gauss Point", SInfo.GQPoint(2), P2, AcceptedError(P2,tol));

    // Check Normals
    tol = 1.0e-7;		// lower tolerance is necessary
    AnalyticNormal = P1/abs(P1);
    ensure_distance("First Normal at GQP", SInfo.NormalGQPoint(1), AnalyticNormal, AcceptedError(AnalyticNormal,tol));

    AnalyticNormal = P2/abs(P2);
    ensure_distance("Second Normal at GQP", SInfo.NormalGQPoint(2), AnalyticNormal,AcceptedError(AnalyticNormal,tol));

    // First Moment Gauss Point
    path_coord = 0.5*PI*Radius*0.1127016653792583114820736;
    P1 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("First Moment Gauss Point", SInfo.GQPointContourIntegral(1), P1, AcceptedError(P1,tol));

    // Second Moment Gauss Point
    path_coord = 0.5*PI*Radius*0.5;
    P1 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("Second Moment Gauss Point", SInfo.GQPointContourIntegral(2), P1, AcceptedError(P1,tol));

    // Third Moment Gauss Point
    path_coord = 0.5*PI*Radius*0.8872983346207416885179264;
    P1 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("Third Moment Gauss Point", SInfo.GQPointContourIntegral(3), P1, AcceptedError(P1,tol));

  }

  /* Test 8: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<8>()
  {
    set_test_name("Check << & >> operators");

    // Create Spline
    S.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,90.0,200);

    // Update interval
    SInfo.UpdateInterval(S,S.Xp[0],S.Xp[S.np-1],2);

    Check_Input_Output_Operator(SInfo);

  }

  /* Test 9: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<9>()
  {
    set_test_name("Copy constructor");

    // Create Spline
    S.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,90.0,200);

    // Update interval
    SInfo.UpdateInterval(S,S.Xp[0],S.Xp[S.np-1],2);

    Spline2DInterval_HO SInfo_Copy(SInfo);

    // check
    ensure_equals("N_SubIntervals",SInfo_Copy.NumOfSubIntervals(),SInfo.NumOfSubIntervals());
    ensure_equals("N_GQP",SInfo_Copy.GQPointsPerSubInterval(),SInfo.GQPointsPerSubInterval());
    ensure_equals("IntervalLength",SInfo_Copy.IntLength(1),SInfo.IntLength(1));
    ensure_equals("GQP1",SInfo_Copy.GQPoint(1),SInfo.GQPoint(1));
    ensure_equals("GQP2",SInfo_Copy.GQPoint(2),SInfo.GQPoint(2));
    ensure_equals("NormalsGQP1",SInfo_Copy.NormalGQPoint(1),SInfo.NormalGQPoint(1));
    ensure_equals("NormalsGQP2",SInfo_Copy.NormalGQPoint(2),SInfo.NormalGQPoint(2));
    ensure_equals("MomentGQP1",SInfo_Copy.GQPointContourIntegral(1),SInfo.GQPointContourIntegral(1));
    ensure_equals("MomentGQP2",SInfo_Copy.GQPointContourIntegral(2),SInfo.GQPointContourIntegral(2));
    ensure_equals("MomentGQP3",SInfo_Copy.GQPointContourIntegral(3),SInfo.GQPointContourIntegral(3));
    ensure_equals("dYdS1",SInfo_Copy.dYdS(1),SInfo.dYdS(1));
    ensure_equals("dYdS2",SInfo_Copy.dYdS(2),SInfo.dYdS(2));
    ensure_equals("dYdS3",SInfo_Copy.dYdS(3),SInfo.dYdS(3));

  }

  /* Test 10: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<10>()
  {
    set_test_name("Assignment operator");

    // Create Spline
    S.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,90.0,200);

    // Update interval
    SInfo.UpdateInterval(S,S.Xp[0],S.Xp[S.np-1],2);

    Spline2DInterval_HO SInfo_Copy;

    // Assign spline interval
    SInfo_Copy = SInfo;

    // check
    ensure_equals("N_SubIntervals",SInfo_Copy.NumOfSubIntervals(),SInfo.NumOfSubIntervals());
    ensure_equals("N_GQP",SInfo_Copy.GQPointsPerSubInterval(),SInfo.GQPointsPerSubInterval());
    ensure_equals("IntervalLength",SInfo_Copy.IntLength(1),SInfo.IntLength(1));
    ensure_equals("GQP1",SInfo_Copy.GQPoint(1),SInfo.GQPoint(1));
    ensure_equals("GQP2",SInfo_Copy.GQPoint(2),SInfo.GQPoint(2));
    ensure_equals("NormalsGQP1",SInfo_Copy.NormalGQPoint(1),SInfo.NormalGQPoint(1));
    ensure_equals("NormalsGQP2",SInfo_Copy.NormalGQPoint(2),SInfo.NormalGQPoint(2));
    ensure_equals("MomentGQP1",SInfo_Copy.GQPointContourIntegral(1),SInfo.GQPointContourIntegral(1));
    ensure_equals("MomentGQP2",SInfo_Copy.GQPointContourIntegral(2),SInfo.GQPointContourIntegral(2));
    ensure_equals("MomentGQP3",SInfo_Copy.GQPointContourIntegral(3),SInfo.GQPointContourIntegral(3));
    ensure_equals("dYdS1",SInfo_Copy.dYdS(1),SInfo.dYdS(1));
    ensure_equals("dYdS2",SInfo_Copy.dYdS(2),SInfo.dYdS(2));
    ensure_equals("dYdS3",SInfo_Copy.dYdS(3),SInfo.dYdS(3));
  }

  /* Test 11: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<11>()
  {

    set_test_name("IntegratePolynomialTerm(), Orders up to 3 inclusiv");

    // Accepted tolerance
    tol = 1.0e-13;

    // Initialize variables
    GeometricMoments GeomCoeff(3); //

    P1.x = 1.1000;  P1.y = 2.2300;
    P2.x = 10.2390; P2.y = 5.32434;
    Centroid.x = -3.0023434; Centroid.y = 1.002323;

    // Create spline
    S.Create_Spline_Line(P1,P2,10);
    // Initialize interval
    SInfo.InitializeInterval(S,P1,P2,1);

    // Maple solution
    double MapleSolution[10];
    MapleSolution[0] = 26.833631906356000000;  MapleSolution[1] = 81.751337743235260865;
    MapleSolution[2] = 268.49282505383518942;  MapleSolution[3] = 930.47115942699369454;
    MapleSolution[4] = 127.11700039021218166;  MapleSolution[5] = 415.96630433802121991;
    MapleSolution[6] = 1438.0152612494315498;  MapleSolution[7] = 859.40250191036158901;
    MapleSolution[8] = 2963.5361653827905064;  MapleSolution[9] = 6871.6055532023448135;

    int i;
    for (i=0; i<=9 ; ++i){
      GeomCoeff(i).D() = SInfo.IntegratePolynomialTermAndDivide(Centroid,GeomCoeff(i).P1(), GeomCoeff(i).P2());
      ensure_distance("Check integration", GeomCoeff(i).D(), MapleSolution[i], AcceptedError(MapleSolution[i],tol));
    }
  }

  /* Test 12: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<12>()
  {
    set_test_name("IntegratePolynomialTerm(), Orders up to 3 inclusiv, 2 segments");

    // Accepted tolerance
    tol = 1.0e-13;

    // Initialize variables
    GeometricMoments LineInt1(3), LineInt2(3), LineIntTotal(3);
    Spline2DInterval_HO SI1, SI2;

    double MapleSolution[10];
    MapleSolution[0] = 26.833631906356000000;  MapleSolution[1] = 81.751337743235260865;
    MapleSolution[2] = 268.49282505383518942;  MapleSolution[3] = 930.47115942699369454;
    MapleSolution[4] = 127.11700039021218166;  MapleSolution[5] = 415.96630433802121991;
    MapleSolution[6] = 1438.0152612494315498;  MapleSolution[7] = 859.40250191036158901;
    MapleSolution[8] = 2963.5361653827905064;  MapleSolution[9] = 6871.6055532023448135;

    P1.x = 1.1000; P1.y = 2.2300;
    P2.x = 10.2390; P2.y = 5.32434;
    P = 0.5*(P1 + P2);		// middle point
    Centroid.x = -3.0023434; Centroid.y = 1.002323;

    // Create spline 1
    S.Create_Spline_Line(P1,P,10);
    // Initialize interval 1
    SI1.InitializeInterval(S,P1,P,2);

    // Create spline 2
    S.Create_Spline_Line(P,P2,10);
    // Initialize interval 2
    SI2.InitializeInterval(S,P,P2,2);

    for (int i=0; i<=9 ; ++i){
      LineInt1(i).D() = SI1.IntegratePolynomialTermAndDivide(Centroid,LineInt1(i).P1(),LineInt1(i).P2());
      LineInt2(i).D() = SI2.IntegratePolynomialTermAndDivide(Centroid,LineInt2(i).P1(),LineInt2(i).P2());
    }
    
    LineIntTotal = LineInt1 + LineInt2;

    // Check results
    for (int i=0; i<=9 ; ++i){
      ensure_distance("Check integration", LineIntTotal(i).D(), MapleSolution[i], AcceptedError(MapleSolution[i],tol));
    }
  }

  /* Test 13: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<13>()
  {
    set_test_name("Rectangle area");

    double Length = 3.33,
           Width  = 2.46;
    double RectangleArea = 0.0;
    Spline2DInterval_HO SI1, SI2, SI3, SI4;

    // Analytic Result
    AnalyticResult = Length * Width;

    // Create the first line
    S.Create_Spline_Line(Vector2D(0.5*Length,-0.5*Width),Vector2D(-0.5*Length,-0.5*Width),5);
    // Initialize interval
    SI1.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],2);

    // Create the second line
    S.Create_Spline_Line(Vector2D(-0.5*Length,-0.5*Width),Vector2D(-0.5*Length,0.5*Width),5);
    SI2.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],2);

    // Create the third line
    S.Create_Spline_Line(Vector2D(-0.5*Length,0.5*Width),Vector2D(0.5*Length,0.5*Width),5);
    SI3.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],2);

    // Create the fourth line
    S.Create_Spline_Line(Vector2D(0.5*Length,0.5*Width),Vector2D(0.5*Length,-0.5*Width),5);
    SI4.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],2);

    // The numerical area is negative because of the considered direction
    RectangleArea = SI1.AreaContribution() + SI2.AreaContribution() + SI3.AreaContribution() + SI4.AreaContribution();

    // Check results
    ensure_distance("RectangleArea",-RectangleArea,AnalyticResult,AcceptedError(AnalyticResult,tol));
  }

  /* Test 14: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<14>()
  {

    set_test_name("IntegratePolynomialTerm(), Orders up to 3 inclusiv, circular arc");

    // Accepted tolerance
    tol = 1.0e-6;

    // Initialize variables
    GeometricMoments GeomCoeff(3);

    Angle1 = 0.0; Angle2 = 18.0;
    Radius = 1.0;
    Centroid.x = 0.0; Centroid.y = 0.0;

    // Create Spline
    S.Create_Spline_Circular_Arc(Centroid,Radius,Angle1,Angle2,65);

    // Initialize interval
    SInfo.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],1);

    // Maple solution
    double MapleSolution[10];
    MapleSolution[0] = 0.30402594575260794420;    MapleSolution[1] = 0.04658709990183884620;
    MapleSolution[2] = 0.0095493920356488663538;  MapleSolution[3] = 0.0022057427070126132888;
    MapleSolution[4] = 0.14959041432289476003;    MapleSolution[5] = 0.022733047314457232490;
    MapleSolution[6] = 0.0046363017903618450050;  MapleSolution[7] = 0.098158851238986359285;
    MapleSolution[8] = 0.014793785731608744303;   MapleSolution[9] = 0.072477056266266457516;

    int i;
    for (i=0; i<=9 ; ++i){
      GeomCoeff(i).D() = SInfo.IntegratePolynomialTermAndDivide(Centroid,GeomCoeff(i).P1(), GeomCoeff(i).P2());
    }

    ensure_distance("Check Coeff 0", GeomCoeff(0).D(), MapleSolution[0], AcceptedError(MapleSolution[0],tol));
    ensure_distance("Check Coeff 1", GeomCoeff(1).D(), MapleSolution[1], AcceptedError(MapleSolution[1],tol));
    ensure_distance("Check Coeff 2", GeomCoeff(2).D(), MapleSolution[2], AcceptedError(MapleSolution[2],tol));
    ensure_distance("Check Coeff 3", GeomCoeff(3).D(), MapleSolution[3], AcceptedError(MapleSolution[3],tol));
    ensure_distance("Check Coeff 4", GeomCoeff(4).D(), MapleSolution[4], AcceptedError(MapleSolution[4],tol));
    ensure_distance("Check Coeff 5", GeomCoeff(5).D(), MapleSolution[5], AcceptedError(MapleSolution[5],tol));
    ensure_distance("Check Coeff 6", GeomCoeff(6).D(), MapleSolution[6], AcceptedError(MapleSolution[6],tol));
    ensure_distance("Check Coeff 7", GeomCoeff(7).D(), MapleSolution[7], AcceptedError(MapleSolution[7],tol));
    ensure_distance("Check Coeff 8", GeomCoeff(8).D(), MapleSolution[8], AcceptedError(MapleSolution[8],tol));
    ensure_distance("Check Coeff 9", GeomCoeff(9).D(), MapleSolution[9], AcceptedError(MapleSolution[9],tol));
  }

  /* Test 15: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<15>()
  {

    set_test_name("IntegratePolynomialTerm(), Orders up to 3 inclusiv, circular arc 9 degrees");

    // Accepted tolerance
    tol = 1.0e-6;

    // Initialize variables
    GeometricMoments GeomCoeff(3);

    Angle1 = 0.0; Angle2 = 9.0;
    Radius = 1.0;
    Centroid.x = 0.0; Centroid.y = 0.0;

    // Create Spline
    S.Create_Spline_Circular_Arc(Centroid,Radius,Angle1,Angle2,15);

    // Initialize interval
    SInfo.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],1);

    // Maple solution
    double MapleSolution[10];
    MapleSolution[0] = 0.15579406493348168699;    MapleSolution[1] = 0.01216070450218491326;
    MapleSolution[2] = 0.0012666649507964224530;  MapleSolution[3] = 0.00014848960318196448130;
    MapleSolution[4] = 0.077579195212567692199;   MapleSolution[5] = 0.0060430771944443482419;
    MapleSolution[6] = 0.00062866897697822820700; MapleSolution[7] = 0.051509133327561754847;
    MapleSolution[8] = 0.004004071633000982925;   MapleSolution[9] = 0.038475263117794731996;

    int i;
    for (i=0; i<=9 ; ++i){
      GeomCoeff(i).D() = SInfo.IntegratePolynomialTermAndDivide(Centroid,GeomCoeff(i).P1(), GeomCoeff(i).P2());
    }

    ensure_distance("Check Coeff 0", GeomCoeff(0).D(), MapleSolution[0], AcceptedError(MapleSolution[0],tol));
    ensure_distance("Check Coeff 1", GeomCoeff(1).D(), MapleSolution[1], AcceptedError(MapleSolution[1],tol));
    ensure_distance("Check Coeff 2", GeomCoeff(2).D(), MapleSolution[2], AcceptedError(MapleSolution[2],tol));
    ensure_distance("Check Coeff 3", GeomCoeff(3).D(), MapleSolution[3], AcceptedError(MapleSolution[3],tol));
    ensure_distance("Check Coeff 4", GeomCoeff(4).D(), MapleSolution[4], AcceptedError(MapleSolution[4],tol));
    ensure_distance("Check Coeff 5", GeomCoeff(5).D(), MapleSolution[5], AcceptedError(MapleSolution[5],tol));
    ensure_distance("Check Coeff 6", GeomCoeff(6).D(), MapleSolution[6], AcceptedError(MapleSolution[6],tol));
    ensure_distance("Check Coeff 7", GeomCoeff(7).D(), MapleSolution[7], AcceptedError(MapleSolution[7],tol));
    ensure_distance("Check Coeff 8", GeomCoeff(8).D(), MapleSolution[8], AcceptedError(MapleSolution[8],tol));
    ensure_distance("Check Coeff 9", GeomCoeff(9).D(), MapleSolution[9], AcceptedError(MapleSolution[9],tol));
  }


  /* Test 16: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<16>()
  {

    set_test_name("Circle area");

    // Accepted tolerance
    tol = 5.0e-2;

    Radius = 100.0;
    double CircleArea = 0.0;

    // Analytic result
    AnalyticResult = PI * Radius * Radius;

    // Create circle
    S.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,180.0,400);
    S1.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,180.0,360.0,400);

    // Initialize intervals
    SInfo.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],1);
    SInfo1.InitializeInterval(S1,S1.Xp[0],S1.Xp[S1.np-1],1);

    CircleArea = SInfo.AreaContribution() + SInfo1.AreaContribution();

    // Analytic location for the Gauss quadrature points
    Vector2D GQP1, GQP2, GQP3;
    double l, L(PI*Radius), alpha, dYdS;
    
    // GQP1
    l = L*GaussQuadratureData::GQ3_Abscissa[0];
    alpha = l/Radius;

    GQP1.x = Radius * cos(alpha);
    GQP1.y = Radius * sin(alpha);
    dYdS = cos(alpha);

    ensure_distance("GQP1", SInfo.GQPointContourIntegral(1), GQP1, AcceptedError(GQP1,1.0e-12));
    ensure_distance("dYdS 1", SInfo.dYdS(1), dYdS, AcceptedError(dYdS,1.0e-5));

    // GQP2
    l = L*GaussQuadratureData::GQ3_Abscissa[1];
    alpha = l/Radius;

    GQP2.x = Radius * cos(alpha);
    GQP2.y = Radius * sin(alpha);
    dYdS = cos(alpha);

    ensure_distance("GQP2", SInfo.GQPointContourIntegral(2), GQP2, AcceptedError(GQP2,1.0e-10));
    ensure_distance("dYdS 2", SInfo.dYdS(2), dYdS, AcceptedError(dYdS,1.0e-12));

    // GQP2
    l = L*GaussQuadratureData::GQ3_Abscissa[2];
    alpha = l/Radius;

    GQP3.x = Radius * cos(alpha);
    GQP3.y = Radius * sin(alpha);
    dYdS = cos(alpha);

    ensure_distance("GQP3", SInfo.GQPointContourIntegral(3), GQP3, AcceptedError(GQP3,1.0e-12));
    ensure_distance("dYdS 3", SInfo.dYdS(3), dYdS, AcceptedError(dYdS,1.0e-5));

    // Check results
    ensure_distance("CircleArea",CircleArea,AnalyticResult,AcceptedError(AnalyticResult,tol));
  }

  /* Test 17: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<17>()
  {

    set_test_name("Circle Centroid");

    // Accepted tolerance
    tol = 1.0e-12;

    Radius = 100.0;
    double CircleArea = 0.0;
    Vector2D AnalyticCentroid(-2.000002323,100.3434);

    // Create splines
    S.Create_Spline_Circular_Arc(AnalyticCentroid,Radius,0.0,180.0,200);
    S1.Create_Spline_Circular_Arc(AnalyticCentroid,Radius,180.0,360.0,200);

    // Nodes
    Node2D_HO P1(S.Xp[0]), P2(S.Xp[S.np-1]);
    Node2D_HO P3(S1.Xp[0]), P4(S1.Xp[S1.np-1]);

    // Initialize intervals
    SInfo.InitializeInterval(S,P1,P2,1);
    SInfo1.InitializeInterval(S1,P3,P4,1);

    // Calculate the circle area
    CircleArea = SInfo.AreaContribution() + SInfo1.AreaContribution();

    // Calculation of the centroid coordinates
    Centroid.x = ( SInfo.XCentroidContribution() + SInfo1.XCentroidContribution() )*0.5/CircleArea;

    Centroid.y = ( SInfo.YCentroidContribution() + SInfo1.YCentroidContribution() ) /CircleArea;

    // Check results
    ensure_distance("Centroid", Centroid, AnalyticCentroid, AcceptedError(AnalyticCentroid,tol));

  }

  /* Test 18: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<18>()
  {

    set_test_name("Circle area for splines with sharp corners");

    // Accepted tolerance
    tol = 5.0e-7;

    Radius = 100.0;
    double CircleArea = 0.0;

    // Analytic result
    AnalyticResult = PI * Radius * Radius;

    // Create circle
    S.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,180.0,400);
    S1.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,180.0,360.0,400);

    // Switch some of the smooth control points to sharp corners
    S.tp[100] = SPLINE2D_POINT_SHARP_CORNER;
    S.tp[200] = SPLINE2D_POINT_SHARP_CORNER;
    S.tp[300] = SPLINE2D_POINT_SHARP_CORNER;

    S1.tp[100] = SPLINE2D_POINT_SHARP_CORNER;
    S1.tp[200] = SPLINE2D_POINT_SHARP_CORNER;
    S1.tp[300] = SPLINE2D_POINT_SHARP_CORNER;

    // Initialize intervals
    SInfo.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],1);
    SInfo1.InitializeInterval(S1,S1.Xp[0],S1.Xp[S1.np-1],1);

    CircleArea = SInfo.AreaContribution() + SInfo1.AreaContribution();

    // Check results
    ensure_distance("CircleArea",CircleArea,AnalyticResult,AcceptedError(AnalyticResult,tol));
  }

  /* Test 19: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<19>()
  {

    set_test_name("IntegratePolynomialTerm(), Orders up to 3 inclusiv, 5-point Gauss quad");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    // Accepted tolerance
    tol = 1.0e-13;

    // Initialize variables
    GeometricMoments GeomCoeff(3); //

    P1.x = 1.1000;  P1.y = 2.2300;
    P2.x = 10.2390; P2.y = 5.32434;
    Centroid.x = -3.0023434; Centroid.y = 1.002323;

    // Create spline
    S.Create_Spline_Line(P1,P2,10);
    // Initialize interval
    SInfo.InitializeInterval(S,P1,P2,1);

    // Maple solution
    double MapleSolution[10];
    MapleSolution[0] = 26.833631906356000000;  MapleSolution[1] = 81.751337743235260865;
    MapleSolution[2] = 268.49282505383518942;  MapleSolution[3] = 930.47115942699369454;
    MapleSolution[4] = 127.11700039021218166;  MapleSolution[5] = 415.96630433802121991;
    MapleSolution[6] = 1438.0152612494315498;  MapleSolution[7] = 859.40250191036158901;
    MapleSolution[8] = 2963.5361653827905064;  MapleSolution[9] = 6871.6055532023448135;

    int i;
    for (i=0; i<=9 ; ++i){
      GeomCoeff(i).D() = SInfo.IntegratePolynomialTermAndDivide(Centroid,GeomCoeff(i).P1(), GeomCoeff(i).P2());
      ensure_distance("Check integration", GeomCoeff(i).D(), MapleSolution[i], AcceptedError(MapleSolution[i],tol));
    }
  }

  /* Test 20: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<20>()
  {
    set_test_name("IntegratePolynomialTerm(), Orders up to 3 inclusiv, 2 segments, 5-point Gauss quad");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    // Accepted tolerance
    tol = 1.0e-13;

    // Initialize variables
    GeometricMoments LineInt1(3), LineInt2(3), LineIntTotal(3);
    Spline2DInterval_HO SI1, SI2;

    double MapleSolution[10];
    MapleSolution[0] = 26.833631906356000000;  MapleSolution[1] = 81.751337743235260865;
    MapleSolution[2] = 268.49282505383518942;  MapleSolution[3] = 930.47115942699369454;
    MapleSolution[4] = 127.11700039021218166;  MapleSolution[5] = 415.96630433802121991;
    MapleSolution[6] = 1438.0152612494315498;  MapleSolution[7] = 859.40250191036158901;
    MapleSolution[8] = 2963.5361653827905064;  MapleSolution[9] = 6871.6055532023448135;

    P1.x = 1.1000; P1.y = 2.2300;
    P2.x = 10.2390; P2.y = 5.32434;
    P = 0.5*(P1 + P2);		// middle point
    Centroid.x = -3.0023434; Centroid.y = 1.002323;

    // Create spline 1
    S.Create_Spline_Line(P1,P,10);
    // Initialize interval 1
    SI1.InitializeInterval(S,P1,P,2);

    // Create spline 2
    S.Create_Spline_Line(P,P2,10);
    // Initialize interval 2
    SI2.InitializeInterval(S,P,P2,2);

    for (int i=0; i<=9 ; ++i){
      LineInt1(i).D() = SI1.IntegratePolynomialTermAndDivide(Centroid,LineInt1(i).P1(),LineInt1(i).P2());
      LineInt2(i).D() = SI2.IntegratePolynomialTermAndDivide(Centroid,LineInt2(i).P1(),LineInt2(i).P2());
    }
    
    LineIntTotal = LineInt1 + LineInt2;

    // Check results
    for (int i=0; i<=9 ; ++i){
      ensure_distance("Check integration", LineIntTotal(i).D(), MapleSolution[i], AcceptedError(MapleSolution[i],tol));
    }
  }

  /* Test 21: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<21>()
  {
    set_test_name("Rectangle area, 5-point Gauss quad");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    double Length = 3.33,
           Width  = 2.46;
    double RectangleArea = 0.0;
    Spline2DInterval_HO SI1, SI2, SI3, SI4;

    // Analytic Result
    AnalyticResult = Length * Width;

    // Create the first line
    S.Create_Spline_Line(Vector2D(0.5*Length,-0.5*Width),Vector2D(-0.5*Length,-0.5*Width),5);
    // Initialize interval
    SI1.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],2);

    // Create the second line
    S.Create_Spline_Line(Vector2D(-0.5*Length,-0.5*Width),Vector2D(-0.5*Length,0.5*Width),5);
    SI2.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],2);

    // Create the third line
    S.Create_Spline_Line(Vector2D(-0.5*Length,0.5*Width),Vector2D(0.5*Length,0.5*Width),5);
    SI3.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],2);

    // Create the fourth line
    S.Create_Spline_Line(Vector2D(0.5*Length,0.5*Width),Vector2D(0.5*Length,-0.5*Width),5);
    SI4.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],2);

    // The numerical area is negative because of the considered direction
    RectangleArea = SI1.AreaContribution() + SI2.AreaContribution() + SI3.AreaContribution() + SI4.AreaContribution();

    // Check results
    ensure_distance("RectangleArea",-RectangleArea,AnalyticResult,AcceptedError(AnalyticResult,tol));
  }

  /* Test 22: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<22>()
  {

    set_test_name("IntegratePolynomialTerm(), Orders up to 3 inclusiv, circular arc, 5-point Gauss quad");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    // Accepted tolerance
    tol = 1.0e-6;

    // Initialize variables
    GeometricMoments GeomCoeff(3);

    Angle1 = 0.0; Angle2 = 18.0;
    Radius = 1.0;
    Centroid.x = 0.0; Centroid.y = 0.0;

    // Create Spline
    S.Create_Spline_Circular_Arc(Centroid,Radius,Angle1,Angle2,65);

    // Initialize interval
    SInfo.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],1);

    // Maple solution
    double MapleSolution[10];
    MapleSolution[0] = 0.30402594575260794420;    MapleSolution[1] = 0.04658709990183884620;
    MapleSolution[2] = 0.0095493920356488663538;  MapleSolution[3] = 0.0022057427070126132888;
    MapleSolution[4] = 0.14959041432289476003;    MapleSolution[5] = 0.022733047314457232490;
    MapleSolution[6] = 0.0046363017903618450050;  MapleSolution[7] = 0.098158851238986359285;
    MapleSolution[8] = 0.014793785731608744303;   MapleSolution[9] = 0.072477056266266457516;

    int i;
    for (i=0; i<=9 ; ++i){
      GeomCoeff(i).D() = SInfo.IntegratePolynomialTermAndDivide(Centroid,GeomCoeff(i).P1(), GeomCoeff(i).P2());
    }

    ensure_distance("Check Coeff 0", GeomCoeff(0).D(), MapleSolution[0], AcceptedError(MapleSolution[0],tol));
    ensure_distance("Check Coeff 1", GeomCoeff(1).D(), MapleSolution[1], AcceptedError(MapleSolution[1],tol));
    ensure_distance("Check Coeff 2", GeomCoeff(2).D(), MapleSolution[2], AcceptedError(MapleSolution[2],tol));
    ensure_distance("Check Coeff 3", GeomCoeff(3).D(), MapleSolution[3], AcceptedError(MapleSolution[3],tol));
    ensure_distance("Check Coeff 4", GeomCoeff(4).D(), MapleSolution[4], AcceptedError(MapleSolution[4],tol));
    ensure_distance("Check Coeff 5", GeomCoeff(5).D(), MapleSolution[5], AcceptedError(MapleSolution[5],tol));
    ensure_distance("Check Coeff 6", GeomCoeff(6).D(), MapleSolution[6], AcceptedError(MapleSolution[6],tol));
    ensure_distance("Check Coeff 7", GeomCoeff(7).D(), MapleSolution[7], AcceptedError(MapleSolution[7],tol));
    ensure_distance("Check Coeff 8", GeomCoeff(8).D(), MapleSolution[8], AcceptedError(MapleSolution[8],tol));
    ensure_distance("Check Coeff 9", GeomCoeff(9).D(), MapleSolution[9], AcceptedError(MapleSolution[9],tol));
  }

  /* Test 23: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<23>()
  {

    set_test_name("IntegratePolynomialTerm(), Orders up to 3 inclusiv, circular arc 9 degrees, 5-point Gauss quad");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    // Accepted tolerance
    tol = 1.0e-6;

    // Initialize variables
    GeometricMoments GeomCoeff(3);

    Angle1 = 0.0; Angle2 = 9.0;
    Radius = 1.0;
    Centroid.x = 0.0; Centroid.y = 0.0;

    // Create Spline
    S.Create_Spline_Circular_Arc(Centroid,Radius,Angle1,Angle2,15);

    // Initialize interval
    SInfo.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],1);

    // Maple solution
    double MapleSolution[10];
    MapleSolution[0] = 0.15579406493348168699;    MapleSolution[1] = 0.01216070450218491326;
    MapleSolution[2] = 0.0012666649507964224530;  MapleSolution[3] = 0.00014848960318196448130;
    MapleSolution[4] = 0.077579195212567692199;   MapleSolution[5] = 0.0060430771944443482419;
    MapleSolution[6] = 0.00062866897697822820700; MapleSolution[7] = 0.051509133327561754847;
    MapleSolution[8] = 0.004004071633000982925;   MapleSolution[9] = 0.038475263117794731996;

    int i;
    for (i=0; i<=9 ; ++i){
      GeomCoeff(i).D() = SInfo.IntegratePolynomialTermAndDivide(Centroid,GeomCoeff(i).P1(), GeomCoeff(i).P2());
    }

    ensure_distance("Check Coeff 0", GeomCoeff(0).D(), MapleSolution[0], AcceptedError(MapleSolution[0],tol));
    ensure_distance("Check Coeff 1", GeomCoeff(1).D(), MapleSolution[1], AcceptedError(MapleSolution[1],tol));
    ensure_distance("Check Coeff 2", GeomCoeff(2).D(), MapleSolution[2], AcceptedError(MapleSolution[2],tol));
    ensure_distance("Check Coeff 3", GeomCoeff(3).D(), MapleSolution[3], AcceptedError(MapleSolution[3],tol));
    ensure_distance("Check Coeff 4", GeomCoeff(4).D(), MapleSolution[4], AcceptedError(MapleSolution[4],tol));
    ensure_distance("Check Coeff 5", GeomCoeff(5).D(), MapleSolution[5], AcceptedError(MapleSolution[5],tol));
    ensure_distance("Check Coeff 6", GeomCoeff(6).D(), MapleSolution[6], AcceptedError(MapleSolution[6],tol));
    ensure_distance("Check Coeff 7", GeomCoeff(7).D(), MapleSolution[7], AcceptedError(MapleSolution[7],tol));
    ensure_distance("Check Coeff 8", GeomCoeff(8).D(), MapleSolution[8], AcceptedError(MapleSolution[8],tol));
    ensure_distance("Check Coeff 9", GeomCoeff(9).D(), MapleSolution[9], AcceptedError(MapleSolution[9],tol));
  }


  /* Test 24: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<24>()
  {

    set_test_name("Circle area, 5-point Gauss quad");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    // Accepted tolerance
    tol = 5.0e-4;

    Radius = 100.0;
    double CircleArea = 0.0;

    // Analytic result
    AnalyticResult = PI * Radius * Radius;

    // Create circle
    S.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,180.0,400);
    S1.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,180.0,360.0,400);

    // Initialize intervals
    SInfo.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],1);
    SInfo1.InitializeInterval(S1,S1.Xp[0],S1.Xp[S1.np-1],1);

    CircleArea = SInfo.AreaContribution() + SInfo1.AreaContribution();

    // Analytic location for the Gauss quadrature points
    Vector2D GQP1, GQP2, GQP3, GQP4, GQP5;
    double l, L(PI*Radius), alpha, dYdS;
    
    // GQP1
    l = L*GaussQuadratureData::GQ5_Abscissa[0];
    alpha = l/Radius;

    GQP1.x = Radius * cos(alpha);
    GQP1.y = Radius * sin(alpha);
    dYdS = cos(alpha);

    ensure_distance("GQP1", SInfo.GQPointContourIntegral(1), GQP1, AcceptedError(GQP1,1.0e-10));
    ensure_distance("dYdS 1", SInfo.dYdS(1), dYdS, AcceptedError(dYdS,1.0e-5));

    // GQP2
    l = L*GaussQuadratureData::GQ5_Abscissa[1];
    alpha = l/Radius;

    GQP2.x = Radius * cos(alpha);
    GQP2.y = Radius * sin(alpha);
    dYdS = cos(alpha);

    ensure_distance("GQP2", SInfo.GQPointContourIntegral(2), GQP2, AcceptedError(GQP2,1.0e-10));
    ensure_distance("dYdS 2", SInfo.dYdS(2), dYdS, AcceptedError(dYdS,1.0e-5));

    // GQP3
    l = L*GaussQuadratureData::GQ5_Abscissa[2];
    alpha = l/Radius;

    GQP3.x = Radius * cos(alpha);
    GQP3.y = Radius * sin(alpha);
    dYdS = cos(alpha);

    ensure_distance("GQP3", SInfo.GQPointContourIntegral(3), GQP3, AcceptedError(GQP3,1.0e-10));
    ensure_distance("dYdS 3", SInfo.dYdS(3), dYdS, AcceptedError(dYdS,1.0e-5));

    // GQP4
    l = L*GaussQuadratureData::GQ5_Abscissa[3];
    alpha = l/Radius;

    GQP4.x = Radius * cos(alpha);
    GQP4.y = Radius * sin(alpha);
    dYdS = cos(alpha);

    ensure_distance("GQP4", SInfo.GQPointContourIntegral(4), GQP4, AcceptedError(GQP4,1.0e-10));
    ensure_distance("dYdS 4", SInfo.dYdS(4), dYdS, AcceptedError(dYdS,1.0e-5));

    // GQP5
    l = L*GaussQuadratureData::GQ5_Abscissa[4];
    alpha = l/Radius;

    GQP5.x = Radius * cos(alpha);
    GQP5.y = Radius * sin(alpha);
    dYdS = cos(alpha);

    ensure_distance("GQP5", SInfo.GQPointContourIntegral(5), GQP5, AcceptedError(GQP5,1.0e-10));
    ensure_distance("dYdS 5", SInfo.dYdS(5), dYdS, AcceptedError(dYdS,1.0e-5));

    // Check results
    ensure_distance("CircleArea",CircleArea,AnalyticResult,AcceptedError(AnalyticResult,tol));
  }

  /* Test 25: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<25>()
  {

    set_test_name("Circle Centroid, 5-point Gauss quad");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    // Accepted tolerance
    tol = 1.0e-12;

    Radius = 100.0;
    double CircleArea = 0.0;
    Vector2D AnalyticCentroid(-2.000002323,100.3434);

    // Create splines
    S.Create_Spline_Circular_Arc(AnalyticCentroid,Radius,0.0,180.0,200);
    S1.Create_Spline_Circular_Arc(AnalyticCentroid,Radius,180.0,360.0,200);

    // Nodes
    Node2D_HO P1(S.Xp[0]), P2(S.Xp[S.np-1]);
    Node2D_HO P3(S1.Xp[0]), P4(S1.Xp[S1.np-1]);

    // Initialize intervals
    SInfo.InitializeInterval(S,P1,P2,1);
    SInfo1.InitializeInterval(S1,P3,P4,1);

    // Calculate the circle area
    CircleArea = SInfo.AreaContribution() + SInfo1.AreaContribution();

    // Calculation of the centroid coordinates
    Centroid.x = ( SInfo.XCentroidContribution() + SInfo1.XCentroidContribution() )*0.5/CircleArea;

    Centroid.y = ( SInfo.YCentroidContribution() + SInfo1.YCentroidContribution() ) /CircleArea;

    // Check results
    ensure_distance("Centroid", Centroid, AnalyticCentroid, AcceptedError(AnalyticCentroid,tol));

  }

  /* Test 26: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<26>()
  {

    set_test_name("Circle area for splines with sharp corners, 5-point Gauss quad");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    // Accepted tolerance
    tol = 5.0e-6;

    Radius = 100.0;
    double CircleArea = 0.0;

    // Analytic result
    AnalyticResult = PI * Radius * Radius;

    // Create circle
    S.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,180.0,400);
    S1.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,180.0,360.0,400);

    // Switch some of the smooth control points to sharp corners
    S.tp[100] = SPLINE2D_POINT_SHARP_CORNER;
    S.tp[200] = SPLINE2D_POINT_SHARP_CORNER;
    S.tp[300] = SPLINE2D_POINT_SHARP_CORNER;

    S1.tp[100] = SPLINE2D_POINT_SHARP_CORNER;
    S1.tp[200] = SPLINE2D_POINT_SHARP_CORNER;
    S1.tp[300] = SPLINE2D_POINT_SHARP_CORNER;

    // Initialize intervals
    SInfo.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],1);
    SInfo1.InitializeInterval(S1,S1.Xp[0],S1.Xp[S1.np-1],1);

    CircleArea = SInfo.AreaContribution() + SInfo1.AreaContribution();

    // Check results
    ensure_distance("CircleArea",CircleArea,AnalyticResult,AcceptedError(AnalyticResult,tol));
  }

   /* Test 27: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<27>()
  {
    set_test_name("allocate(), 5-point Gauss quad");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    NumSubIntervals = 4;
    N_GQP = 2;

    // allocate
    SInfo.allocate(N_GQP,NumSubIntervals);

    // check
    ensure_equals("Number of subintervals", SInfo.NumOfSubIntervals(),NumSubIntervals);
    ensure_equals("Gauss Points per subinterval", SInfo.GQPointsPerSubInterval(), N_GQP);
    ensure_equals("Last Gauss Point", SInfo.GQPoint(NumSubIntervals*N_GQP),Vector2D(0.0));
    ensure_equals("Last Normal at GQP", SInfo.NormalGQPoint(NumSubIntervals*N_GQP),Vector2D(0.0));
    ensure_equals("Total GQP_ContourInt",SInfo.NumGQPoints_ContourIntegral(),20);
    ensure_equals("Last Gauss Point Contour Int", SInfo.GQPointContourIntegral(NumSubIntervals*5),Vector2D(0.0));
    ensure_equals("Last dYdS", SInfo.dYdS(NumSubIntervals), 0.0);
    ensure_equals("Last SubInterval length", SInfo.IntLength(NumSubIntervals), 0.0);

  }

  /* Test 28: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<28>()
  {
    set_test_name("deallocate(), 5-point Gauss quad");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    NumSubIntervals = 4;
    N_GQP = 2;

    // allocate
    SInfo.allocate(N_GQP,NumSubIntervals);
    SInfo.deallocate();

    // check
    ensure_equals("IntervalLength",SInfo.IntLength(), NullDouble);
    ensure_equals("GQP",SInfo.GQPoint(),NullVector);
    ensure_equals("NormalsGQP",SInfo.NormalGQPoint(), NullVector);
    ensure_equals("Total GQP_ContourInt",SInfo.NumGQPoints_ContourIntegral(),0);
    ensure_equals("MomentGQP",SInfo.GQPointContourIntegral(),NullVector);
    ensure_equals("dYdS",SInfo.dYdS(), NullDouble);
    ensure_equals("Number of subintervals", SInfo.NumOfSubIntervals(),0);
    ensure_equals("Gauss Points per subinterval", SInfo.GQPointsPerSubInterval(), 0);
  }


  /* Test 29: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<29>()
  {
    set_test_name("2) allocate(), 5-point Gauss quad");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    NumSubIntervals = 4;
    N_GQP = 2;

    // allocate
    SInfo.allocate(N_GQP,NumSubIntervals);

    NumSubIntervals = 10;
    N_GQP = 1;

    // allocate again
    SInfo.allocate(N_GQP,NumSubIntervals);

    // check
    ensure_equals("Number of subintervals", SInfo.NumOfSubIntervals(),NumSubIntervals);
    ensure_equals("Gauss Points per subinterval", SInfo.GQPointsPerSubInterval(), N_GQP);
    ensure_equals("Total GQP_ContourInt",SInfo.NumGQPoints_ContourIntegral(),5*NumSubIntervals);
    ensure_equals("Last Gauss Point", SInfo.GQPoint(NumSubIntervals*N_GQP),Vector2D(0.0));
    ensure_equals("Last Normal at GQP", SInfo.NormalGQPoint(NumSubIntervals*N_GQP),Vector2D(0.0));
    ensure_equals("Last Gauss Point Moment", SInfo.GQPointContourIntegral(NumSubIntervals*5),Vector2D(0.0));
    ensure_equals("Last dYdS", SInfo.dYdS(NumSubIntervals*N_GQP), 0.0);
    ensure_equals("Last SubInterval length", SInfo.IntLength(NumSubIntervals), 0.0);
  }


  /* Test 30: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<30>()
  {
    set_test_name("1) InitializeInterval() for a line, 5-point Gauss quad");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    V1 = Vector2D(0.0);
    V2 = Vector2D(4.0,3.0);

    // Create Spline
    S.Create_Spline_Line(V1,V2,10);

    // Initialize interval
    SInfo.InitializeInterval(S,S.Xp[0],S.Xp[9],2);

    // check
    ensure_equals("Number of subintervals", SInfo.NumOfSubIntervals(),1);
    ensure_equals("Gauss Points per subinterval", SInfo.GQPointsPerSubInterval(), 2);

    ensure_distance("SubInterval length", SInfo.IntLength(1), 5.0, tol);

    P = Vector2D(8.452994616207484709817024e-1, 6.339745962155613532362768e-1);
    ensure_distance("First Gauss Point", SInfo.GQPoint(1), P, AcceptedError(P,tol));
    P = Vector2D(3.1547005383792515290182976e+0,2.3660254037844386467637232e+0);
    ensure_distance("Second Gauss Point", SInfo.GQPoint(2), P, AcceptedError(P,tol));

    AnalyticNormal = Vector2D(0.6,-0.8);
    ensure_distance("First Normal at GQP" , SInfo.NormalGQPoint(1), AnalyticNormal, AcceptedError(AnalyticNormal,tol));
    ensure_distance("Second Normal at GQP", SInfo.NormalGQPoint(2), AnalyticNormal, AcceptedError(AnalyticNormal,tol));

    ensure_distance("First dYdS at GQP_ContourInt" , SInfo.dYdS(1), AnalyticNormal.x, AcceptedError(AnalyticNormal.x,tol));
    ensure_distance("Second dYdS at GQP_ContourInt", SInfo.dYdS(2), AnalyticNormal.x, AcceptedError(AnalyticNormal.x,tol));
    ensure_distance("Third dYdS at GQP_ContourInt" , SInfo.dYdS(3), AnalyticNormal.x, AcceptedError(AnalyticNormal.x,tol));
    ensure_distance("Second dYdS at GQP_ContourInt", SInfo.dYdS(4), AnalyticNormal.x, AcceptedError(AnalyticNormal.x,tol));
    ensure_distance("Third dYdS at GQP_ContourInt" , SInfo.dYdS(5), AnalyticNormal.x, AcceptedError(AnalyticNormal.x,tol));
  }

  /* Test 31: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<31>()
  {
    set_test_name("2) InitializeInterval() for a quater of a circle, 5-point Gauss quad");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    double path_coord;

    tol = 1.0e-10;

    // Create Spline
    S.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,90.0,200);

    // Initialize interval
    SInfo.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],2);

    // check
    ensure_equals("Number of subintervals", SInfo.NumOfSubIntervals(),1);
    ensure_equals("Gauss Points per subinterval", SInfo.GQPointsPerSubInterval(), 2);

    double Answer =  0.5*PI*Radius;
    ensure_distance("SubInterval length", SInfo.IntLength(1), Answer, AcceptedError(Answer,tol));

    // First Gauss Point
    path_coord = 0.5*PI*Radius*0.2113248654051871177454256;
    P1 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("First Gauss Point", SInfo.GQPoint(1), P1, AcceptedError(P1,tol));

    // Second Gauss Point
    path_coord = 0.5*PI*Radius*0.7886751345948128822545744;
    P2 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("Second Gauss Point", SInfo.GQPoint(2), P2, AcceptedError(P2,tol));

    // Check Normals
    tol = 1.0e-7;		// lower tolerance is necessary
    AnalyticNormal = P1/abs(P1);
    ensure_distance("First Normal at GQP", SInfo.NormalGQPoint(1), AnalyticNormal, AcceptedError(AnalyticNormal,tol));

    AnalyticNormal = P2/abs(P2);
    ensure_distance("Second Normal at GQP", SInfo.NormalGQPoint(2), AnalyticNormal,AcceptedError(AnalyticNormal,tol));

    // First Gauss Contour Integration Point
    path_coord = 0.5*PI*Radius*4.6910077030668e-2;
    P1 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("First Moment Gauss Point", SInfo.GQPointContourIntegral(1), P1, AcceptedError(P1,tol));

    // Second Gauss Contour Integration Point
    path_coord = 0.5*PI*Radius*2.30765344947158e-1;
    P1 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("Second Moment Gauss Point", SInfo.GQPointContourIntegral(2), P1, AcceptedError(P1,tol));

    // Third Gauss Contour Integration Point
    path_coord = 0.5*PI*Radius*0.5;
    P1 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("Third Moment Gauss Point", SInfo.GQPointContourIntegral(3), P1, AcceptedError(P1,tol));

    // Fourth Gauss Contour Integration Point
    path_coord = 0.5*PI*Radius*7.69234655052842e-1;
    P1 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("Fourth Moment Gauss Point", SInfo.GQPointContourIntegral(4), P1, AcceptedError(P1,tol));

    // Fifth Gauss Contour Integration Point
    path_coord = 0.5*PI*Radius*9.53089922969332e-1;
    P1 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("Fifth Moment Gauss Point", SInfo.GQPointContourIntegral(5), P1, AcceptedError(P1,tol));
  }

  /* Test 32: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<32>()
  {
    set_test_name("UpdateInterval(), 5-point Gauss quad");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    double path_coord;

    tol = 1.0e-10;

    V1 = Vector2D(0.0);
    V2 = Vector2D(4.0,3.0);

    // Create Spline 1
    S.Create_Spline_Line(V1,V2,10);

    // Initialize interval
    SInfo.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],1);

    // Create Spline 2
    S.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,90.0,200);

    // Update interval (different number of GQP!!)
    SInfo.UpdateInterval(S,S.Xp[0],S.Xp[S.np-1],2);

    // check
    ensure_equals("Number of subintervals", SInfo.NumOfSubIntervals(),1);
    ensure_equals("Gauss Points per subinterval", SInfo.GQPointsPerSubInterval(), 2);

    double Answer =  0.5*PI*Radius;
    ensure_distance("SubInterval length", SInfo.IntLength(1), Answer, AcceptedError(Answer,tol));

    // First Gauss Point
    path_coord = 0.5*PI*Radius*0.2113248654051871177454256;
    P1 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("First Gauss Point", SInfo.GQPoint(1), P1, AcceptedError(P1,tol));

    // Second Gauss Point
    path_coord = 0.5*PI*Radius*0.7886751345948128822545744;
    P2 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("Second Gauss Point", SInfo.GQPoint(2), P2, AcceptedError(P2,tol));

    // Check Normals
    tol = 1.0e-7;		// lower tolerance is necessary
    AnalyticNormal = P1/abs(P1);
    ensure_distance("First Normal at GQP", SInfo.NormalGQPoint(1), AnalyticNormal, AcceptedError(AnalyticNormal,tol));

    AnalyticNormal = P2/abs(P2);
    ensure_distance("Second Normal at GQP", SInfo.NormalGQPoint(2), AnalyticNormal,AcceptedError(AnalyticNormal,tol));

    // First Gauss Contour Integration Point
    path_coord = 0.5*PI*Radius*4.6910077030668e-2;
    P1 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("First Moment Gauss Point", SInfo.GQPointContourIntegral(1), P1, AcceptedError(P1,tol));

    // Second Gauss Contour Integration Point
    path_coord = 0.5*PI*Radius*2.30765344947158e-1;
    P1 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("Second Moment Gauss Point", SInfo.GQPointContourIntegral(2), P1, AcceptedError(P1,tol));

    // Third Gauss Contour Integration Point
    path_coord = 0.5*PI*Radius*0.5;
    P1 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("Third Moment Gauss Point", SInfo.GQPointContourIntegral(3), P1, AcceptedError(P1,tol));

    // Fourth Gauss Contour Integration Point
    path_coord = 0.5*PI*Radius*7.69234655052842e-1;
    P1 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("Fourth Moment Gauss Point", SInfo.GQPointContourIntegral(4), P1, AcceptedError(P1,tol));

    // Fifth Gauss Contour Integration Point
    path_coord = 0.5*PI*Radius*9.53089922969332e-1;
    P1 = Radius * Vector2D(cos(path_coord/Radius), sin(path_coord/Radius));
    ensure_distance("Fifth Moment Gauss Point", SInfo.GQPointContourIntegral(5), P1, AcceptedError(P1,tol));
  }

  /* Test 33: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<33>()
  {
    set_test_name("Check << & >> operators, 5-point Gauss quad");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    // Create Spline
    S.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,90.0,200);

    // Update interval
    SInfo.UpdateInterval(S,S.Xp[0],S.Xp[S.np-1],2);

    Check_Input_Output_Operator(SInfo);
  }

  /* Test 34: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<34>()
  {
    set_test_name("Copy constructor, 5-point Gauss quad");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    // Create Spline
    S.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,90.0,200);

    // Update interval
    SInfo.UpdateInterval(S,S.Xp[0],S.Xp[S.np-1],2);

    Spline2DInterval_HO SInfo_Copy(SInfo);

    // check
    ensure_equals("N_SubIntervals",SInfo_Copy.NumOfSubIntervals(),SInfo.NumOfSubIntervals());
    ensure_equals("N_GQP",SInfo_Copy.GQPointsPerSubInterval(),SInfo.GQPointsPerSubInterval());
    ensure_equals("IntervalLength",SInfo_Copy.IntLength(1),SInfo.IntLength(1));
    ensure_equals("GQP1",SInfo_Copy.GQPoint(1),SInfo.GQPoint(1));
    ensure_equals("GQP2",SInfo_Copy.GQPoint(2),SInfo.GQPoint(2));
    ensure_equals("NormalsGQP1",SInfo_Copy.NormalGQPoint(1),SInfo.NormalGQPoint(1));
    ensure_equals("NormalsGQP2",SInfo_Copy.NormalGQPoint(2),SInfo.NormalGQPoint(2));
    ensure_equals("MomentGQP1",SInfo_Copy.GQPointContourIntegral(1),SInfo.GQPointContourIntegral(1));
    ensure_equals("MomentGQP2",SInfo_Copy.GQPointContourIntegral(2),SInfo.GQPointContourIntegral(2));
    ensure_equals("MomentGQP3",SInfo_Copy.GQPointContourIntegral(3),SInfo.GQPointContourIntegral(3));
    ensure_equals("MomentGQP4",SInfo_Copy.GQPointContourIntegral(4),SInfo.GQPointContourIntegral(4));
    ensure_equals("MomentGQP5",SInfo_Copy.GQPointContourIntegral(5),SInfo.GQPointContourIntegral(5));
    ensure_equals("dYdS1",SInfo_Copy.dYdS(1),SInfo.dYdS(1));
    ensure_equals("dYdS2",SInfo_Copy.dYdS(2),SInfo.dYdS(2));
    ensure_equals("dYdS3",SInfo_Copy.dYdS(3),SInfo.dYdS(3));
    ensure_equals("dYdS4",SInfo_Copy.dYdS(4),SInfo.dYdS(4));
    ensure_equals("dYdS5",SInfo_Copy.dYdS(5),SInfo.dYdS(5));
  }

  /* Test 35: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<35>()
  {
    set_test_name("Assignment operator, 5-point Gauss quad");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    // Create Spline
    S.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,90.0,200);

    // Update interval
    SInfo.UpdateInterval(S,S.Xp[0],S.Xp[S.np-1],2);

    Spline2DInterval_HO SInfo_Copy;

    // Assign spline interval
    SInfo_Copy = SInfo;

    // check
    ensure_equals("N_SubIntervals",SInfo_Copy.NumOfSubIntervals(),SInfo.NumOfSubIntervals());
    ensure_equals("N_GQP",SInfo_Copy.GQPointsPerSubInterval(),SInfo.GQPointsPerSubInterval());
    ensure_equals("IntervalLength",SInfo_Copy.IntLength(1),SInfo.IntLength(1));
    ensure_equals("GQP1",SInfo_Copy.GQPoint(1),SInfo.GQPoint(1));
    ensure_equals("GQP2",SInfo_Copy.GQPoint(2),SInfo.GQPoint(2));
    ensure_equals("NormalsGQP1",SInfo_Copy.NormalGQPoint(1),SInfo.NormalGQPoint(1));
    ensure_equals("NormalsGQP2",SInfo_Copy.NormalGQPoint(2),SInfo.NormalGQPoint(2));
    ensure_equals("MomentGQP1",SInfo_Copy.GQPointContourIntegral(1),SInfo.GQPointContourIntegral(1));
    ensure_equals("MomentGQP2",SInfo_Copy.GQPointContourIntegral(2),SInfo.GQPointContourIntegral(2));
    ensure_equals("MomentGQP3",SInfo_Copy.GQPointContourIntegral(3),SInfo.GQPointContourIntegral(3));
    ensure_equals("MomentGQP4",SInfo_Copy.GQPointContourIntegral(4),SInfo.GQPointContourIntegral(4));
    ensure_equals("MomentGQP5",SInfo_Copy.GQPointContourIntegral(5),SInfo.GQPointContourIntegral(5));
    ensure_equals("dYdS1",SInfo_Copy.dYdS(1),SInfo.dYdS(1));
    ensure_equals("dYdS2",SInfo_Copy.dYdS(2),SInfo.dYdS(2));
    ensure_equals("dYdS3",SInfo_Copy.dYdS(3),SInfo.dYdS(3));
    ensure_equals("dYdS2",SInfo_Copy.dYdS(4),SInfo.dYdS(4));
    ensure_equals("dYdS3",SInfo_Copy.dYdS(5),SInfo.dYdS(5));
  }
 
  /* Test 36: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<36>()
  {

    set_test_name("Invariant to shifted coordinate system");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    // Accepted tolerance
    tol = 1.0e-6;

    // Initialize variables
    GeometricMoments GeomCoeff(3);
    Vector2D Shift = Vector2D(-1.0e2, 2.0e5);

    Angle1 = 0.0; Angle2 = 9.0;
    Radius = 1.0;
    Centroid.x = 0.0; Centroid.y = 0.0;

    // Create the original Spline and the Translated Spline
    S.Create_Spline_Circular_Arc(Centroid,Radius,Angle1,Angle2,15);
    S1 = S + Shift;

    // Initialize interval
    SInfo.InitializeInterval(S,S.Xp[0],S.Xp[S.np-1],1);
    SInfo1.InitializeInterval(S1,S1.Xp[0],S1.Xp[S1.np-1],1);

    // Maple solution for the spline that is NOT translated
    double MapleSolution[10];
    MapleSolution[0] = 0.15579406493348168699;    MapleSolution[1] = 0.01216070450218491326;
    MapleSolution[2] = 0.0012666649507964224530;  MapleSolution[3] = 0.00014848960318196448130;
    MapleSolution[4] = 0.077579195212567692199;   MapleSolution[5] = 0.0060430771944443482419;
    MapleSolution[6] = 0.00062866897697822820700; MapleSolution[7] = 0.051509133327561754847;
    MapleSolution[8] = 0.004004071633000982925;   MapleSolution[9] = 0.038475263117794731996;

    int i;

    // check area contribution
    ensure_distance("Check Area Contribution",
		    SInfo1.AreaContribution(Shift),
		    SInfo.AreaContribution(),
		    AcceptedError(SInfo.AreaContribution(),tol));

    ensure_distance("Check X Centroid Contribution",
		    SInfo1.XCentroidContribution(Shift),
		    SInfo.XCentroidContribution(),
		    AcceptedError(SInfo.XCentroidContribution(),tol));

    ensure_distance("Check Y Centroid Contribution",
		    SInfo1.YCentroidContribution(Shift),
		    SInfo.YCentroidContribution(),
		    AcceptedError(SInfo.YCentroidContribution(),tol));

    for (i=0; i<=9 ; ++i){
      GeomCoeff(i).D() = SInfo1.IntegratePolynomialTermAndDivide(Centroid+Shift,GeomCoeff(i).P1(), GeomCoeff(i).P2());
    }

    ensure_distance("Check Coeff 0", GeomCoeff(0).D(), MapleSolution[0], AcceptedError(MapleSolution[0],tol));
    ensure_distance("Check Coeff 1", GeomCoeff(1).D(), MapleSolution[1], AcceptedError(MapleSolution[1],tol));
    ensure_distance("Check Coeff 2", GeomCoeff(2).D(), MapleSolution[2], AcceptedError(MapleSolution[2],tol));
    ensure_distance("Check Coeff 3", GeomCoeff(3).D(), MapleSolution[3], AcceptedError(MapleSolution[3],tol));
    ensure_distance("Check Coeff 4", GeomCoeff(4).D(), MapleSolution[4], AcceptedError(MapleSolution[4],tol));
    ensure_distance("Check Coeff 5", GeomCoeff(5).D(), MapleSolution[5], AcceptedError(MapleSolution[5],tol));
    ensure_distance("Check Coeff 6", GeomCoeff(6).D(), MapleSolution[6], AcceptedError(MapleSolution[6],tol));
    ensure_distance("Check Coeff 7", GeomCoeff(7).D(), MapleSolution[7], AcceptedError(MapleSolution[7],tol));
    ensure_distance("Check Coeff 8", GeomCoeff(8).D(), MapleSolution[8], AcceptedError(MapleSolution[8],tol));
    ensure_distance("Check Coeff 9", GeomCoeff(9).D(), MapleSolution[9], AcceptedError(MapleSolution[9],tol));
  }

  /* Test 37: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<37>()
  {
    set_test_name("IntegrateFunctionWithRespectToY() along a line");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    V1 = Vector2D(0.0);
    V2 = Vector2D(4.0,3.0);

    // Create Spline
    S.Create_Spline_Line(V1,V2,10);

    // Initialize interval
    SInfo.InitializeInterval(S,S.Xp[0],S.Xp[9],2);

    double Int;

    // Analytic result
    AnalyticResult = exp(V1.y)*(exp(V2.y - V1.y) - 1.0);

    // Numerical result
    Int = SInfo.IntegrateFunctionWithRespectToY(TestFunction,Int);

    // == check
    ensure_distance("Integral value", Int, AnalyticResult, AcceptedError(AnalyticResult, 1.0e-6));

  }

  /* Test 38: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<38>()
  {
    set_test_name("IntegrateFunctionProjectionOnNormalDirections() along a line");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    V1 = Vector2D(0.0);
    V2 = Vector2D(4.0,3.0);

    // Create Spline
    S.Create_Spline_Line(V1,V2,10);

    // Initialize interval
    SInfo.InitializeInterval(S,S.Xp[0],S.Xp[9],2);

    double IntX(3.0), IntY(5.0);

    // Analytic result
    AnalyticResult = 31.8092282053127796;

    // Numerical result
    SInfo.IntegrateFunctionProjectionOnNormalDirections(S, TestFunction, IntX, IntY);

    // == check
    ensure_distance("X Component Integral", IntX, 3.0 + AnalyticResult*0.6, AcceptedError(3.0 + AnalyticResult*0.6, 1.0e-6));
    ensure_distance("Y Component Integral", IntY, 5.0 + AnalyticResult*(-0.8), AcceptedError(5.0 + AnalyticResult*(-0.8), 1.0e-6));
  }

  /* Test 39: */
  template<>
  template<>
  void Spline2DInterval_HO_object::test<39>()
  {
    set_test_name("IntegrateFunctionAlongInterval() along a line");

    // Set 5-point Gauss integration
    Spline2DInterval_HO::setFivePointGaussQuadContourIntegration();

    V1 = Vector2D(0.0);
    V2 = Vector2D(4.0,3.0);

    // Create Spline
    S.Create_Spline_Line(V1,V2,10);

    // Initialize interval
    SInfo.InitializeInterval(S,S.Xp[0],S.Xp[9],2);

    double Int;

    // Analytic result
    AnalyticResult = 31.8092282053127796;

    // Numerical result
    Int = SInfo.IntegrateFunctionAlongInterval(TestFunction,
					       Int);

    // == check
    ensure_distance("Integral value", Int, AnalyticResult, AcceptedError(AnalyticResult, 1.0e-6));

  }

}



// Test suite constructor
tut::Spline2DInterval_HO_TestSuite Spline2DInterval_HOTestSuite("Class:Spline2DInterval_HO");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

