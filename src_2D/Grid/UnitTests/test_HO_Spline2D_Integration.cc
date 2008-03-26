/*!\file test_HO_Spline2D_Integration.cc
  \brief Regression tests for integration functions in class Spline2D_HO. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../HO_Spline2D.h"
#include "../../HighOrderReconstruction/TaylorDerivatives.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_ContourIntegration_Along_Spline2D_HO : public TestData {

    // Local variables
  public:

    // === Member variables
    typedef TaylorDerivativesContainer<TwoD,double> GeometricMoments;

    double TheoreticResult;	// Result obtained with the exact integration
    double NumericResult;
    double StartPointX, StartPointY;		
    double EndPointX, EndPointY;
    double CentroidX, CentroidY;
    double AnalyticResult;

    Spline2D_HO Curve1, Curve2, Curve, CurveSegment;
    Spline2D_HO Circle;

    // === Member functions

    // Constructor
    Data_ContourIntegration_Along_Spline2D_HO(void);

  private:
    
  };
  
  Data_ContourIntegration_Along_Spline2D_HO::Data_ContourIntegration_Along_Spline2D_HO(void){
    /* Set the global path for this test suite. 
       It's a good practice to put it in the constructor of the data class in order to be set
       automatically for each individual test. Declare it relative to the /src_2D directory,
       otherwise the framework might not find the input and output files. */
    
    set_test_suite_path("Grid/UnitTests");
  }

  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. 'NAME'_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_ContourIntegration_Along_Spline2D_HO> ContourInt_Spline2D_HO_TestSuite;
  typedef ContourInt_Spline2D_HO_TestSuite::object ContourInt_Spline2D_HO_object;


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


 /* Test 1: PolynomLineIntegration along 1 segment line */
  template<>
  template<>
  void ContourInt_Spline2D_HO_object::test<1>()
  {

    set_test_name("PolynomLineIntegration(), Orders up to 3 inclusiv");

    // Initialize variables
    GeometricMoments GeomCoeff(3); //

    StartPointX = 1.1000;     StartPointY = 2.2300;
    EndPointX   = 10.2390;    EndPointY   = 5.32434;
    CentroidX   = -3.0023434; CentroidY   = 1.002323;

    double MapleSolution[10];
    MapleSolution[0] = 26.833631906356;  MapleSolution[1] = 81.7513377432353;
    MapleSolution[2] = 268.492825053835; MapleSolution[3] = 930.471159426994;
    MapleSolution[4] = 254.234000780424; MapleSolution[5] = 831.932608676042;
    MapleSolution[6] = 2876.03052249886; MapleSolution[7] = 2578.20750573109;
    MapleSolution[8] = 8890.60849614838; MapleSolution[9] = 27486.4222128094;

    for (int i=0; i<=9 ; ++i){
      GeomCoeff(i).D() = PolynomLineIntegration(StartPointX,StartPointY,EndPointX,EndPointY,CentroidX,CentroidY,
						GeomCoeff(i).P1(), GeomCoeff(i).P2());
      
      ensure_distance("Check integration", GeomCoeff(i).D(), MapleSolution[i], tol*(1.0 + fabs(MapleSolution[i])));
    }
  }

  /* Test 2: PolynomLineIntegration along 2 segments */
  template<>
  template<>
  void ContourInt_Spline2D_HO_object::test<2>()
  {
    set_test_name("PolynomLineIntegration(), Orders up to 3 inclusiv, 2 segments");

    // Initialize variables
    double MiddlePointX, MiddlePointY;
    GeometricMoments LineInt1(3), LineInt2(3), LineIntTotal(3);

    double MapleSolution[10];
    MapleSolution[0] = 26.833631906356;  MapleSolution[1] = 81.7513377432353;
    MapleSolution[2] = 268.492825053835; MapleSolution[3] = 930.471159426994;
    MapleSolution[4] = 254.234000780424; MapleSolution[5] = 831.932608676042;
    MapleSolution[6] = 2876.03052249886; MapleSolution[7] = 2578.20750573109;
    MapleSolution[8] = 8890.60849614838; MapleSolution[9] = 27486.4222128094;

    StartPointX = 1.1000; StartPointY = 2.2300;
    EndPointX = 10.2390; EndPointY = 5.32434;
    MiddlePointX =0.5*(EndPointX+StartPointX); MiddlePointY = 0.5*(EndPointY+StartPointY);
    CentroidX = -3.0023434; CentroidY = 1.002323;

    for (int i=0; i<=9 ; ++i){
      LineInt1(i).D() = PolynomLineIntegration(StartPointX,StartPointY,MiddlePointX,MiddlePointY,
					       CentroidX,CentroidY,
					       LineInt1(i).P1(), LineInt1(i).P2());

      LineInt2(i).D() = PolynomLineIntegration(MiddlePointX,MiddlePointY,EndPointX,EndPointY,
					       CentroidX,CentroidY,
					       LineInt2(i).P1(), LineInt2(i).P2());
    }
    
    LineIntTotal = LineInt1 + LineInt2;

    // Check results
    for (int i=0; i<=9 ; ++i){
      ensure_distance("Check integration", LineIntTotal(i).D(), MapleSolution[i], tol*(1.0 + fabs(MapleSolution[i])));
    }

  }


  /* Test 3: PolynomLineIntegration along a line spline */
  template<>
  template<>
  void ContourInt_Spline2D_HO_object::test<3>()
  {

    set_test_name("Geometric moments along a line");

    Vector2D V1(1.1,2.23), V2(10.2390,5.32434);
    Spline2D_HO Line;
    CentroidX = -3.0023434; CentroidY = 1.002323;
    Vector2D Centroid(CentroidX,CentroidY);

    Line.Create_Spline_Line(V1,V2,5);
    Line.setBCtype(BC_NULL);

    double MapleSolution[10];
    MapleSolution[0] = 26.833631906356;  MapleSolution[1] = 81.7513377432353;
    MapleSolution[2] = 268.492825053835; MapleSolution[3] = 930.471159426994;
    MapleSolution[4] = 254.234000780424; MapleSolution[5] = 831.932608676042;
    MapleSolution[6] = 2876.03052249886; MapleSolution[7] = 2578.20750573109;
    MapleSolution[8] = 8890.60849614838; MapleSolution[9] = 27486.4222128094;

    GeometricMoments GeomCoeff(3);

    for (int i=0; i<=9 ; ++i){
      GeomCoeff(i).D() = Line.PolynomOrderIntegration(V1,V2,Centroid,15,
						      GeomCoeff(i).P1(),
						      GeomCoeff(i).P2());
      ensure_distance("Check integration", GeomCoeff(i).D(), MapleSolution[i], tol*(1.0 + fabs(MapleSolution[i])));
    }

  }

  /* Test 4: PolynomLineIntegration along a circular spline */
  template<>
  template<>
  void ContourInt_Spline2D_HO_object::test<4>()
  {

    set_test_name("Centroid contributions");

    Vector2D V1(3.0,2.0);
    Curve.Create_Spline_Circular_Arc(V1,1.0,10.0,90.0,7);
    Curve.setBCtype(BC_NULL);

    GeometricMoments GeomCoeff(3);

    for (int i=0; i<=9 ; ++i){
      GeomCoeff(i).D() = Curve.PolynomOrderIntegration(Curve.Xp[0],Curve.Xp[Curve.np-1],V1,15,
						       GeomCoeff(i).P1(),
						       GeomCoeff(i).P2());

    }

    double CentroidIntX, CentroidIntY;

    CentroidIntX = Curve.PolynomOrderIntegration(Curve.Xp[0],Curve.Xp[Curve.np-1],Vector2D(0.0,0.0),15,
						 1,0);

    CentroidIntY = Curve.PolynomOrderIntegration(Curve.Xp[0],Curve.Xp[Curve.np-1],Vector2D(0.0,0.0),15,
						 0,1);

    ensure_distance("CentroidX", CentroidIntX, 11.60730746954379, tol*(11.60730746954379 + 1.0));
    ensure_distance("CentroidY", CentroidIntY, 7.956360557093456, tol*(7.956360557093456 + 1.0));
  }

  /* Test 5: Computation of a rectangle area using Green-Gauss transformation */
  template<>
  template<>
  void ContourInt_Spline2D_HO_object::test<5>()
  {
    set_test_name("Rectangle area");

    // There are 2 curves used for defining the boundary of the rectangle
    double Length = 3.33,
           Width  = 2.46;
    double RectangleArea = 0.0;
    AnalyticResult = Length * Width;

    // Create the first curve
    Curve1.Create_Spline_Line(Vector2D(0.0,-0.5*Width),Vector2D(-0.5*Length,-0.5*Width),5);
    Curve.Create_Spline_Line(Vector2D(-0.5*Length,-0.5*Width),Vector2D(-0.5*Length,0.5*Width),5);
    Curve1 = Curve1 + Curve;

    Curve.Create_Spline_Line(Vector2D(-0.5*Length,0.5*Width),Vector2D(0.0,0.5*Width),5);
    Curve1 = Curve1 + Curve;

    // Create the second curve
    Curve2.Create_Spline_Line(Vector2D(0.0,-0.5*Width),Vector2D(0.5*Length,-0.5*Width),5);
    Curve.Create_Spline_Line(Vector2D(0.5*Length,-0.5*Width),Vector2D(0.5*Length,0.5*Width),5);
    Curve2 = Curve2 + Curve;

    Curve.Create_Spline_Line(Vector2D(0.5*Length,0.5*Width),Vector2D(0.0,0.5*Width),5);
    Curve2 = Curve2 + Curve;

    Vector2D P1(0.0,-Width*0.5), P2(0.0,Width*0.5);
    RectangleArea = Curve1.ZeroOrderIntegration(P2,P1,15) + Curve2.ZeroOrderIntegration(P1,P2,15);

    // Check results
    ensure_distance("RectangleArea",RectangleArea,AnalyticResult,tol*(1.0 + fabs(AnalyticResult)));
  }

  /* Test 6: Computation of a circle area using Green-Gauss transformation */
  template<>
  template<>
  void ContourInt_Spline2D_HO_object::test<6>()
  {
    set_test_name("Circle area");

    // Accepted tolerance
    tol = 1.0e-8;

    Spline2D_HO Circle1, Circle2;
    double Radius = 100.0;
    double CircleArea = 0.0;
    AnalyticResult = PI * Radius * Radius;

    Circle1.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,0.0,180.0,200);
    Circle2.Create_Spline_Circular_Arc(Vector2D(0.0,0.0),Radius,180.0,360.0,200);

    Node2D_HO P1(Circle1.Xp[0]), P2(Circle1.Xp[Circle1.np-1]);
    Node2D_HO P3(Circle2.Xp[0]), P4(Circle2.Xp[Circle2.np-1]);

    CircleArea = Circle1.ZeroOrderIntegration(P1,P2,14) + Circle2.ZeroOrderIntegration(P3,P4,14);

    // Check results
    ensure_distance("CircleArea",CircleArea,AnalyticResult,tol*(1.0 + fabs(AnalyticResult)));
  }


  /* Test 7: Computation of a complex figure area using Green-Gauss transformation */
  template<>
  template<>
  void ContourInt_Spline2D_HO_object::test<7>()
  {
    set_test_name("Complex figure area");

    // Accepted tolerance
    tol = 1.0e-8;

    double FigureArea(0.0), AnalyticArea;
    
    AnalyticArea = 8.59*5.71576 + PI*2.22625;

    // ************* Create The First Spline **************************
    // create the sinusoidal
    Curve1.Create_Spline_Sinusoidal_Line(Vector2D(0.0,0.0),2.33,6.1,0.0,360.0,50);
    // create the horizontal line
    CurveSegment.Create_Spline_Line(Vector2D(6.1,0.0),Vector2D(11.28,0.0),2);
    Curve1 = Curve1 + CurveSegment;
    // create the first oblique line
    CurveSegment.Create_Spline_Line(Vector2D(11.28,0.0),Vector2D(9.2,-5.71576),2);
    Curve1 = Curve1 + CurveSegment;
    // create the first circle
    CurveSegment.Create_Spline_Circular_Arc(Vector2D(7.5,-5.71576),1.7,180.0,360.0,200);
    CurveSegment.Reverse_Spline();
    Curve1 = Curve1 + CurveSegment;

    // ************* Create The Second Spline **************************
    // create the second circle
    Curve2.Create_Spline_Circular_Arc(Vector2D(4.55,-5.71576),1.25,180.0,360.0,200);
    Curve2.Reverse_Spline();
    // create the second oblique line
    CurveSegment.Create_Spline_Line(Vector2D(3.3,-5.71576),Vector2D(0.0,0.0),2);
    Curve2 = Curve2 + CurveSegment;

    // using the right-hand rule in order to obtain positive quantity
    FigureArea = ( Curve1.ZeroOrderIntegration(Curve1.Xp[Curve1.np-1], Curve1.Xp[0],15) + 
		   Curve2.ZeroOrderIntegration(Curve2.Xp[Curve2.np-1], Curve2.Xp[0],15) );

    // Check results
    ensure_distance("FigureArea",FigureArea,AnalyticArea,tol*(1.0 + fabs(AnalyticArea)));

  }

  /* Test 8: Computation of the centroid of a circle */
  template<>
  template<>
  void ContourInt_Spline2D_HO_object::test<8>()
  {
    set_test_name("Circle Centroid");

    Spline2D_HO SemiCircle1, SemiCircle2;
    double Radius = 100.0;
    double CircleArea = 0.0;
    Vector2D Centroid, AnalyticCentroid(-2.000002323,100.3434);

    SemiCircle1.Create_Spline_Circular_Arc(AnalyticCentroid,Radius,0.0,180.0,200);
    SemiCircle2.Create_Spline_Circular_Arc(AnalyticCentroid,Radius,180.0,360.0,200);

    Node2D_HO P1(SemiCircle1.Xp[0]), P2(SemiCircle1.Xp[SemiCircle1.np-1]);
    Node2D_HO P3(SemiCircle2.Xp[0]), P4(SemiCircle2.Xp[SemiCircle2.np-1]);

    // Calculation of the circle area
    CircleArea = SemiCircle1.ZeroOrderIntegration(P1,P2,14) + SemiCircle2.ZeroOrderIntegration(P3,P4,14);

    // Calculation of the centroid coordinates
    Centroid.x = (SemiCircle1.PolynomOrderIntegration(P1,P2,Vector2D(0.0,0.0),15,1,0) + 
		  SemiCircle2.PolynomOrderIntegration(P3,P4,Vector2D(0.0,0.0),15,1,0)) /2.0/CircleArea;

    Centroid.y = (SemiCircle1.PolynomOrderIntegration(P1,P2,Vector2D(0.0,0.0),15,0,1) + 
		  SemiCircle2.PolynomOrderIntegration(P3,P4,Vector2D(0.0,0.0),15,0,1)) / CircleArea;

    // Check results
    ensure_distance("Centroid", Centroid, AnalyticCentroid, AcceptedError(AnalyticCentroid,tol));

  }

  /* Test 9: */
  template<>
  template<>
  void ContourInt_Spline2D_HO_object::test<9>()
  {
    set_test_name("Geometric moments along a circular arc");

    // Accepted tolerance
    tol = 1.0e-9;

    // Initialize variables
    GeometricMoments GeomCoeff(3); //

    double Angle1, Angle2, Radius;
    Vector2D Centroid;

    Angle1 = 0.0; Angle2 = 18.0;
    Radius = 1.0;
    Centroid.x = 0.0; Centroid.y = 0.0;

    // Create Spline
    Curve.Create_Spline_Circular_Arc(Centroid,Radius,Angle1,Angle2,65);

    // Maple solution
    double MapleSolution[10];
    MapleSolution[0] = 0.30402594575260794420;    MapleSolution[1] = 0.04658709990183884620;
    MapleSolution[2] = 0.0095493920356488663538;  MapleSolution[3] = 0.0022057427070126132888;
    MapleSolution[4] = 0.14959041432289476003;    MapleSolution[5] = 0.022733047314457232490;
    MapleSolution[6] = 0.0046363017903618450050;  MapleSolution[7] = 0.098158851238986359285;
    MapleSolution[8] = 0.014793785731608744303;   MapleSolution[9] = 0.072477056266266457516;

    // Get the numerical solution
    for (int i=0; i<=9 ; ++i){
      GeomCoeff(i).D() = Curve.PolynomOrderIntegration(Curve.Xp[0],Curve.Xp[Curve.np-1],Centroid,15,
						       GeomCoeff(i).P1(),
						       GeomCoeff(i).P2()) / (GeomCoeff(i).P1() + 1.0);
    }

    ensure_distance("Check Coeff 0", GeomCoeff(0).D(), MapleSolution[0], tol*(1.0 + fabs(MapleSolution[0])));
    ensure_distance("Check Coeff 1", GeomCoeff(1).D(), MapleSolution[1], tol*(1.0 + fabs(MapleSolution[1])));
    ensure_distance("Check Coeff 2", GeomCoeff(2).D(), MapleSolution[2], tol*(1.0 + fabs(MapleSolution[2])));
    ensure_distance("Check Coeff 3", GeomCoeff(3).D(), MapleSolution[3], tol*(1.0 + fabs(MapleSolution[3])));
    ensure_distance("Check Coeff 4", GeomCoeff(4).D(), MapleSolution[4], tol*(1.0 + fabs(MapleSolution[4])));
    ensure_distance("Check Coeff 5", GeomCoeff(5).D(), MapleSolution[5], tol*(1.0 + fabs(MapleSolution[5])));
    ensure_distance("Check Coeff 6", GeomCoeff(6).D(), MapleSolution[6], tol*(1.0 + fabs(MapleSolution[6])));
    ensure_distance("Check Coeff 7", GeomCoeff(7).D(), MapleSolution[7], tol*(1.0 + fabs(MapleSolution[7])));
    ensure_distance("Check Coeff 8", GeomCoeff(8).D(), MapleSolution[8], tol*(1.0 + fabs(MapleSolution[8])));
    ensure_distance("Check Coeff 9", GeomCoeff(9).D(), MapleSolution[9], tol*(1.0 + fabs(MapleSolution[9])));

  }
}



// Test suite constructor
tut::ContourInt_Spline2D_HO_TestSuite ContourInt_Spline2D_HOTestSuite("Class:Spline2D_HO (integration)");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

