/*!\file test_NumericalLibrary.cc
   \brief Regression tests for functions prototyped in NumericalLibrary.h */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */

/* Include CFFC header files */
#include "TestData.h"
#include "../NumericalLibrary.h"
#include "TestFunctions/TestFunctions_1D.h"
#include "TestFunctions/TestFunctions_2D.h"
#include "TestFunctions/TestFunctions_3D.h"


namespace tut

{
  // Data used for testing
  class Data_GaussLobattoQuadrature: public TestData {
  public:
    FunctionType1D func;        // Use 1D function for integration

    double TheoreticResult;    // Result obtained with the exact integration
    double NumericResultA, NumericResultB;
    double StartPoint;		// One end of the interval
    double EndPoint;		// The other end of the interval
    double MatlabResult, MapleResult;	// The result obtained using Matlab, Maple
    double DummyParam;
  };

  typedef test_group<Data_GaussLobattoQuadrature> GLQ_TestSuite;
  typedef GLQ_TestSuite::object GLQ_Object;

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

  /* Test 1:*/
  template<>
  template<>
  void GLQ_Object::test<1>()
  {
    // Initialize variables
    func = Test_Default1D;
    MatlabResult = 6.178088507844034e+13;

    // Negative interval
    StartPoint = -125.00;
    EndPoint = -89.0;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Default1D_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, fabs(tol*TheoreticResult));
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, fabs(tol*NumericResultA));
    ensure_distance("GLQ: Matlab Value", NumericResultA, MatlabResult, fabs(tol*MatlabResult));
  }


  /* Test 2:*/
  template<>
  template<>
  void GLQ_Object::test<2>()
  {
    // Initialize variables
    func = Test_Default1D;

    // Negative interval and zero
    StartPoint = -2500.00;
    EndPoint = 0.0;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Default1D_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);


    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);

  }

  /* Test 3:*/
  template<>
  template<>
  void GLQ_Object::test<3>()
  {
    // Initialize variables
    func = Test_Default1D;

    // Values close to zero
    StartPoint = -0.00000016;
    EndPoint = 1.0E-14;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Default1D_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);


    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, fabs(tol*TheoreticResult));
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, fabs(tol*NumericResultA));
  }

  /* Test 4:*/
  template<>
  template<>
  void GLQ_Object::test<4>()
  {
    // Initialize variables
    func = Test_Default1D;

    // Values close to zero
    StartPoint = -1.15E-14;
    EndPoint = 0.2345;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Default1D_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);


    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);

  }

  /* Test 5:*/
  template<>
  template<>
  void GLQ_Object::test<5>()
  {
    // Initialize variables
    func = Test_Default1D;

    // Pozitive interval
    StartPoint = 125.00;
    EndPoint = 89.0;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Default1D_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);


    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, fabs(tol*TheoreticResult));
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, fabs(tol*NumericResultA));
  }

  /* Test 6:*/
  template<>
  template<>
  void GLQ_Object::test<6>()
  {
    // Initialize variables
    func = Test_Example1;

    // Negative interval
    StartPoint = -100.2334;
    EndPoint = -103.113;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Example1_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);


    // Check
    tol = 1.0e-13;
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, fabs(tol*TheoreticResult));
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, fabs(tol*NumericResultA));
  }

  /* Test 7:*/
  template<>
  template<>
  void GLQ_Object::test<7>()
  {
    // Initialize variables
    func = Test_Example1;

    // Negative interval and Zero
    StartPoint = -0.00004540;
    EndPoint = 0.0;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Example1_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);


    // Check
    tol = 1.0e-12;
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*fabs(TheoreticResult));
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, tol*fabs(NumericResultA));
  }

  /* Test 8:*/
  template<>
  template<>
  void GLQ_Object::test<8>()
  {
    // Initialize variables
    func = Test_Example1;

    // Interval close to zero
    StartPoint = -0.0016;
    EndPoint = 1.0E-10;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Example1_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);


    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);

  }

  /* Test 9:*/
  template<>
  template<>
  void GLQ_Object::test<9>()
  {
    // Initialize variables
    func = Test_Example1;

    // Interval close to zero
    StartPoint = -1.15E-14;
    EndPoint = 0.2345;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Example1_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);


    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);

  }

  /* Test 10:*/
  template<>
  template<>
  void GLQ_Object::test<10>()
  {
    // Initialize variables
    func = Test_Example1;

    // Pozitive interval
    StartPoint = 105.00;
    EndPoint = 102.0;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Example1_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*fabs(TheoreticResult));
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, tol*fabs(NumericResultA));

  }

  /* Test 11:*/
  template<>
  template<>
  void GLQ_Object::test<11>()
  {
    // Initialize variables
    func = Test_Example3;

    // Negative interval
    StartPoint = -125.00;
    EndPoint = -89.0;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Example3_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);


    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);

  }

  /* Test 12:*/
  template<>
  template<>
  void GLQ_Object::test<12>()
  {
    // Initialize variables
    func = Test_Example3;

    // Negative interval and Zero
    StartPoint = -2500.00;
    EndPoint = 0.0;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Example3_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);


    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);

  }

  /* Test 13:*/
  template<>
  template<>
  void GLQ_Object::test<13>()
  {
    // Initialize variables
    func = Test_Example3;

    // Interval close to zero
    StartPoint = -0.00000016;
    EndPoint = 1.0E-14;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Example3_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);


    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);

  }

  /* Test 14:*/
  template<>
  template<>
  void GLQ_Object::test<14>()
  {
    // Initialize variables
    func = Test_Example3;

    // Interval close to zero
    StartPoint = -1.15E-14;
    EndPoint = 0.2345;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Example3_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);


    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);

  }

  /* Test 15:*/
  template<>
  template<>
  void GLQ_Object::test<15>()
  {
    // Initialize variables
    func = Test_Example3;

    // Pozitive interval
    StartPoint = 125.00;
    EndPoint = 89.0;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Example3_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);


    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*fabs(TheoreticResult));
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, tol*fabs(NumericResultA));
  }


  /* Test 16:*/
  template<>
  template<>
  void GLQ_Object::test<16>()
  {
    // Initialize variables
    func = Test_Example2;
    tol = 1.0e-11;

    // Negative interval
    StartPoint = -124.00;
    EndPoint = -122.5;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Example2_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);


    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*fabs(TheoreticResult));
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, tol*fabs(NumericResultA));
  }

  /* Test 17:*/
  template<>
  template<>
  void GLQ_Object::test<17>()
  {
    // Initialize variables
    func = Test_Example2;
    tol = 1.0e-13;

    // Negative interval and Zero
    StartPoint = -1.232;
    EndPoint = 0.0;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Example2_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);


    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);

  }

  /* Test 18:*/
  template<>
  template<>
  void GLQ_Object::test<18>()
  {
    // Initialize variables
    func = Test_Example2;

    // Interval close to zero
    StartPoint = -0.016;
    EndPoint = 1.0E-14;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Example2_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);


    // Check
    tol = 1.0e-13;
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*fabs(TheoreticResult));
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, tol*fabs(NumericResultA));
  }

  /* Test 19:*/
  template<>
  template<>
  void GLQ_Object::test<19>()
  {
    // Initialize variables
    func = Test_Example2;

    // Interval close to zero
    StartPoint = -1.15E-14;
    EndPoint = 0.2345;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Example2_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);

  }

  /* Test 20:*/
  template<>
  template<>
  void GLQ_Object::test<20>()
  {
    // Initialize variables
    func = Test_Example2;
    tol = 1.0e-10;

    // Pozitive interval
    StartPoint = 122.33;
    EndPoint = 121.450;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Example2_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);


    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*fabs(TheoreticResult));
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, tol*fabs(NumericResultA));
  }

  /* Test 21:*/
  template<>
  template<>
  void GLQ_Object::test<21>()
  {
    // Initialize variables
    func = Test_Example4;

    // Pozitive interval
    StartPoint = 125.00;
    EndPoint = 89.0;

    // GaussLobattoQuadrature
    TheoreticResult = Test_Example4_Integral(StartPoint,EndPoint);
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,DummyParam);


    // Check
    tol = 1.0e-13;
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*fabs(TheoreticResult));
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, tol*fabs(NumericResultA));
  }

  /* Test 22:*/

  /* Test 23:*/

  /* Test 24:*/
  template<>
  template<>
  void GLQ_Object::test<24>()
  {
    // Initialize variables
    func = Test_Example10;
    StartPoint = 0.0;
    EndPoint = 1.0;

    TheoreticResult = Test_Example10_Integral(StartPoint,EndPoint);
    // GaussLobattoQuadrature
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,NumericResultA);

    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
  }

  /* Test 25:*/
  template<>
  template<>
  void GLQ_Object::test<25>()
  {
    // Initialize variables
    func = Test_Example11;
    StartPoint = 0.5;
    EndPoint = 1.0;

    TheoreticResult = Test_Example11_Integral(StartPoint,EndPoint);
    // GaussLobattoQuadrature
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,NumericResultA);
 
    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
  }

  /* Test 26:*/
  template<>
  template<>
  void GLQ_Object::test<26>()
  {
    // Initialize variables
    func = Test_Example12;
    tol = 1.0e-10;
    StartPoint = 0.3;
    EndPoint = 1.0;

    TheoreticResult = Test_Example12_Integral(StartPoint,EndPoint);
    // GaussLobattoQuadrature
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,NumericResultA);

    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*fabs(TheoreticResult));

  }


  /* Test 27:*/
  template<>
  template<>
  void GLQ_Object::test<27>()
  {
    // Initialize variables
    func = Test_Example13;
    StartPoint = 0.0;
    EndPoint = 1.0;

    TheoreticResult = Test_Example13_Integral(StartPoint,EndPoint);
    // GaussLobattoQuadrature
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,NumericResultA);

    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
  }

  /* Test 28:*/
  template<>
  template<>
  void GLQ_Object::test<28>()
  {
    // Initialize variables
    func = Test_Example14;
    StartPoint = 0.4;
    EndPoint = 1.0;

    TheoreticResult = Test_Example14_Integral(StartPoint,EndPoint);
    // GaussLobattoQuadrature
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,NumericResultA);

    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
  }

  /* Test 29:*/
  template<>
  template<>
  void GLQ_Object::test<29>()
  {
    // Initialize variables
    func = Test_Example15;

    StartPoint = 0.0;
    EndPoint = 1.0;

    TheoreticResult = Test_Example15_Integral(StartPoint,EndPoint);
    // GaussLobattoQuadrature
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,NumericResultA);

    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
  }

  /* Test 30:*/
  template<>
  template<>
  void GLQ_Object::test<30>()
  {
    // Initialize variables
    func = Test_Example16;
    StartPoint = -100.0;
    EndPoint = 100.0;

    TheoreticResult = Test_Example16_Integral(StartPoint,EndPoint);
    // GaussLobattoQuadrature
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,NumericResultA);

    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
  }

  /* Test 31:*/
  template<>
  template<>
  void GLQ_Object::test<31>()
  {
    // Initialize variables
    func = Test_Example17;
    tol = 1.0e-10;
    StartPoint = 0.0;
    EndPoint = 1.0;

    TheoreticResult = Test_Example17_Integral(StartPoint,EndPoint);
    // GaussLobattoQuadrature
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,NumericResultA);

    // Check
    ensure_distance("GLQ: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
  }


  /* Test 32:*/
  template<>
  template<>
  void GLQ_Object::test<32>()
  {
    // Initialize variables
    func = Test_Example18;
    MapleResult = .12989473684210526316e38;
    StartPoint = -100.0;
    EndPoint = 100.0;

    // GaussLobattoQuadrature
    NumericResultA = GaussLobattoQuadrature(func,StartPoint,EndPoint,NumericResultA);
    NumericResultB = GaussLobattoQuadrature(func,EndPoint,StartPoint,NumericResultB);
 
    // Check
    ensure_distance("GLQ: Maple Value", NumericResultA, MapleResult, tol*MapleResult);
    ensure_distance("GLQ: Change limits", NumericResultA, -NumericResultB, tol*fabs(NumericResultA));
  }

}


namespace tut

{
  // Data used for testing
  class Data_AdaptiveGaussianQuadrature: public TestData {
  public:
    FunctionType1D func;        // Use 1D function for integration

    double TheoreticResult;    // Result obtained with the exact integration
    double NumericResultA, NumericResultB;
    double StartPoint;		// One end of the interval
    double EndPoint;		// The other end of the interval
    double MatlabResult, MapleResult;	// The result obtained using Matlab, Maple
    double DummyParam;
  };

  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (LinearSystems_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_AdaptiveGaussianQuadrature> AGQ_TestSuite;
  typedef AGQ_TestSuite::object AGQ_Object;

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

  /* Test 1:*/
  template<>
  template<>
  void AGQ_Object::test<1>()
  {
    // Initialize variables
    func = Test_Default1D;
    MatlabResult = 6.178088507844034e+13;

    // Negative interval
    StartPoint = -125.00;
    EndPoint = -89.0;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Default1D_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, fabs(tol*TheoreticResult));
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, fabs(tol*NumericResultA));
    ensure_distance("ADG: Matlab Value", NumericResultA, MatlabResult, fabs(tol*MatlabResult));
  }


  /* Test 2:*/
  template<>
  template<>
  void AGQ_Object::test<2>()
  {
    // Initialize variables
    func = Test_Default1D;

    // Negative interval and zero
    StartPoint = -2500.00;
    EndPoint = 0.0;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Default1D_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);
  }

  /* Test 3:*/
  template<>
  template<>
  void AGQ_Object::test<3>()
  {
    // Initialize variables
    func = Test_Default1D;

    // Values close to zero
    StartPoint = -0.00000016;
    EndPoint = 1.0E-14;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Default1D_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, fabs(tol*TheoreticResult));
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, fabs(tol*NumericResultA));
  }

  /* Test 4:*/
  template<>
  template<>
  void AGQ_Object::test<4>()
  {
    // Initialize variables
    func = Test_Default1D;

    // Values close to zero
    StartPoint = -1.15E-14;
    EndPoint = 0.2345;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Default1D_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);
  }

  /* Test 5:*/
  template<>
  template<>
  void AGQ_Object::test<5>()
  {
    // Initialize variables
    func = Test_Default1D;

    // Pozitive interval
    StartPoint = 125.00;
    EndPoint = 89.0;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Default1D_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, fabs(tol*TheoreticResult));
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, fabs(tol*NumericResultA));
  }

  /* Test 6:*/
  template<>
  template<>
  void AGQ_Object::test<6>()
  {
    // Initialize variables
    func = Test_Example1;

    // Negative interval
    StartPoint = -122.45;
    EndPoint = -113.3334;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Example1_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    tol = 1.0e-13;
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, fabs(tol*TheoreticResult));
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, fabs(tol*NumericResultA));
  }

  /* Test 7:*/
  template<>
  template<>
  void AGQ_Object::test<7>()
  {
    // Initialize variables
    func = Test_Example1;

    // Negative interval and Zero
    StartPoint = -123.453;
    EndPoint = 0.0;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Example1_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*fabs(TheoreticResult));
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, tol*fabs(NumericResultA));
  }

  /* Test 8:*/
  template<>
  template<>
  void AGQ_Object::test<8>()
  {
    // Initialize variables
    func = Test_Example1;

    // Interval close to zero
    StartPoint = -0.016;
    EndPoint = 1.0E-14;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Example1_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);
  }

  /* Test 9:*/
  template<>
  template<>
  void AGQ_Object::test<9>()
  {
    // Initialize variables
    func = Test_Example1;

    // Interval close to zero
    StartPoint = -1.15E-14;
    EndPoint = 0.2345;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Example1_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);
  }

  /* Test 10:*/
  template<>
  template<>
  void AGQ_Object::test<10>()
  {
    // Initialize variables
    func = Test_Example1;

    // Pozitive interval
    StartPoint = 125.00;
    EndPoint = 89.0;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Example1_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*fabs(TheoreticResult));
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, tol*fabs(NumericResultA));
  }

  /* Test 11:*/
  template<>
  template<>
  void AGQ_Object::test<11>()
  {
    // Initialize variables
    func = Test_Example3;

    // Negative interval
    StartPoint = -125.00;
    EndPoint = -89.0;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Example3_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);
  }

  /* Test 12:*/
  template<>
  template<>
  void AGQ_Object::test<12>()
  {
    // Initialize variables
    func = Test_Example3;

    // Negative interval and Zero
    StartPoint = -2500.00;
    EndPoint = 0.0;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Example3_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);
  }

  /* Test 13:*/
  template<>
  template<>
  void AGQ_Object::test<13>()
  {
    // Initialize variables
    func = Test_Example3;

    // Interval close to zero
    StartPoint = -0.00000016;
    EndPoint = 1.0E-14;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Example3_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);
  }

  /* Test 14:*/
  template<>
  template<>
  void AGQ_Object::test<14>()
  {
    // Initialize variables
    func = Test_Example3;

    // Interval close to zero
    StartPoint = -1.15E-14;
    EndPoint = 0.2345;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Example3_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);
  }

  /* Test 15:*/
  template<>
  template<>
  void AGQ_Object::test<15>()
  {
    // Initialize variables
    func = Test_Example3;

    // Pozitive interval
    StartPoint = 125.00;
    EndPoint = 89.0;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Example3_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*fabs(TheoreticResult));
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, tol*fabs(NumericResultA));
  }


  /* Test 16:*/
  template<>
  template<>
  void AGQ_Object::test<16>()
  {
    // Initialize variables
    func = Test_Example2;

    // Negative interval
    StartPoint = -125.00;
    EndPoint = -121.0;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Example2_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*fabs(TheoreticResult));
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, tol*fabs(NumericResultA));
  }

  /* Test 17:*/
  template<>
  template<>
  void AGQ_Object::test<17>()
  {
    // Initialize variables
    func = Test_Example2;

    // Negative interval and Zero
    StartPoint = -30.00;
    EndPoint = 0.0;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Example2_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);
  }

  /* Test 18:*/
  template<>
  template<>
  void AGQ_Object::test<18>()
  {
    // Initialize variables
    func = Test_Example2;

    // Interval close to zero
    StartPoint = -0.016;
    EndPoint = 1.0E-14;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Example2_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    tol = 1.0e-13;
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*fabs(TheoreticResult));
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, tol*fabs(NumericResultA));
  }

  /* Test 19:*/
  template<>
  template<>
  void AGQ_Object::test<19>()
  {
    // Initialize variables
    func = Test_Example2;

    // Interval close to zero
    StartPoint = -1.15E-14;
    EndPoint = 0.2345;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Example2_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);
  }

  /* Test 20:*/
  template<>
  template<>
  void AGQ_Object::test<20>()
  {
    // Initialize variables
    func = Test_Example2;

    // Pozitive interval
    StartPoint = 150.00;
    EndPoint = 1.0;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Example2_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, tol*NumericResultA);
  }

  /* Test 21:*/
  template<>
  template<>
  void AGQ_Object::test<21>()
  {
    // Initialize variables
    func = Test_Example4;

    // Pozitive interval
    StartPoint = 125.00;
    EndPoint = 89.0;

    // AdaptiveGaussianQuadrature
    TheoreticResult = Test_Example4_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,DummyParam);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*fabs(TheoreticResult));
    ensure_distance("ADG: Change limits", NumericResultA, -NumericResultB, tol*fabs(NumericResultA));
  }

  /* Test 22:*/
  template<>
  template<>
  void AGQ_Object::test<22>()
  {
    // Initialize variables
    func = Test_Default1D;
    SuperFunctionType1D super_func = Error1D;
    FunctionType1D func2 = Test_Example6;
    MatlabResult = 9.639642857142857;

    // Pozitive interval
    StartPoint = -2.0;
    EndPoint = -1.0;

    NumericResultA = AdaptiveGaussianQuadrature(super_func,func,func2,StartPoint,EndPoint,digits);

    // Check
    ensure_distance("ADG: Matlab Value", NumericResultA, MatlabResult, tol*MatlabResult);
  }

  /* Test 23:*/
  template<>
  template<>
  void AGQ_Object::test<23>()
  {
    // Initialize variables
    func = Test_Example9;
    MapleResult = 0.81785715647;
    tol = 1.0e-12;

    // Pozitive interval
    StartPoint = 1.0;
    EndPoint = 1000.0;

    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,NumericResultA,digits);

    // Check
    ensure_distance("ADG: Matlab Value", NumericResultA, MapleResult, tol*MapleResult);
  }

  /* Test 24:*/
  template<>
  template<>
  void AGQ_Object::test<24>()
  {
    // Initialize variables
    func = Test_Example10;
    StartPoint = 0.0;
    EndPoint = 1.0;

    TheoreticResult = Test_Example10_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,NumericResultA);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
  }

  /* Test 25:*/
  template<>
  template<>
  void AGQ_Object::test<25>()
  {
    // Initialize variables
    func = Test_Example11;
    StartPoint = 0.0;
    EndPoint = 1.0;

    TheoreticResult = Test_Example11_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,NumericResultA);
 
    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
  }

  /* Test 26:*/
  template<>
  template<>
  void AGQ_Object::test<26>()
  {
    // Initialize variables
    func = Test_Example12;
    StartPoint = 0.0;
    EndPoint = 1.0;
    tol = 5.0e-10;

    TheoreticResult = Test_Example12_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,NumericResultA,digits);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
  }


  /* Test 27:*/
  template<>
  template<>
  void AGQ_Object::test<27>()
  {
    // Initialize variables
    func = Test_Example13;
    StartPoint = 0.0;
    EndPoint = 1.0;

    TheoreticResult = Test_Example13_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,NumericResultA);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
  }

  /* Test 28:*/
  template<>
  template<>
  void AGQ_Object::test<28>()
  {
    // Initialize variables
    func = Test_Example14;
    StartPoint = 0.0;
    EndPoint = 1.0;

    TheoreticResult = Test_Example14_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,NumericResultA);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
  }

  /* Test 29:*/
  template<>
  template<>
  void AGQ_Object::test<29>()
  {
    // Initialize variables
    func = Test_Example15;

    StartPoint = 0.0;
    EndPoint = 1.0;

    TheoreticResult = Test_Example15_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,NumericResultA);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
  }

  /* Test 30:*/
  template<>
  template<>
  void AGQ_Object::test<30>()
  {
    // Initialize variables
    func = Test_Example16;
    StartPoint = -100.0;
    EndPoint = 100.0;

    TheoreticResult = Test_Example16_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,NumericResultA);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
  }

  /* Test 31:*/
  template<>
  template<>
  void AGQ_Object::test<31>()
  {
    // Initialize variables
    func = Test_Example17;
    StartPoint = 0.0;
    EndPoint = 1.0;

    TheoreticResult = Test_Example17_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,NumericResultA);

    // Check
    ensure_distance("ADG: Integral Value", NumericResultA, TheoreticResult, tol*TheoreticResult);
  }


  /* Test 32:*/
  template<>
  template<>
  void AGQ_Object::test<32>()
  {
    // Initialize variables
    func = Test_Example18;
    MapleResult = .12989473684210526316e38;
    StartPoint = -100.0;
    EndPoint = 100.0;

    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,NumericResultA);
 
    // Check
    ensure_distance("ADG: Maple Value", NumericResultA, MapleResult, tol*MapleResult);
  }

}

namespace tut

{
  // Data used for testing
  class Data_AdaptiveGaussianQuadrature2D: public TestData{
  public:
    FunctionType2D func;        // Use 1D function for integration
    double TheoreticResult;	// Result obtained with the exact integration
    double NumericResultA, NumericResultB, NumericResult;
    double StartPointX;		// One end of the interval (X dir)
    double EndPointX;		// The other end of the interval (X dir)
    double StartPointY;		// One end of the interval (Y dir)
    double EndPointY;		// The other end of the interval (Y dir)
    double AcceptError;	        // Acceptable error
    double MapleResult;		// Result obtained with Maple

    void IntegrateFunctionNumerically(void);
    void CheckSolution(double & Result);
  };

  void Data_AdaptiveGaussianQuadrature2D::IntegrateFunctionNumerically(void){
    NumericResult = AdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,
					       digits,NumericResult);
    // Switch X limits
    NumericResultA = AdaptiveGaussianQuadrature(func,EndPointX,StartPointX,StartPointY,EndPointY,
						digits,NumericResultA);
    // Switch Y limits
    NumericResultB = AdaptiveGaussianQuadrature(func,StartPointX,EndPointX,EndPointY,StartPointY,
						digits,NumericResultB);

    // check independency of limits
    AcceptError = tol*fabs(NumericResult);
    ensure_distance("AGQ2D: Change limits X", NumericResult, -NumericResultA, AcceptError);
    ensure_distance("AGQ2D: Change limits Y", NumericResult, -NumericResultB, AcceptError);
  }

  void Data_AdaptiveGaussianQuadrature2D::CheckSolution(double & Result){
    AcceptError = tol*fabs(Result);
    ensure_distance("AGQ2D: Integral Value", NumericResult, Result, AcceptError);
  }

  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (LinearSystems_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_AdaptiveGaussianQuadrature2D> AGQ2D_TestSuite;
  typedef AGQ2D_TestSuite::object  AGQ2D_Object;

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

  /* Test 1:*/
  template<>
  template<>
  void AGQ2D_Object::test<1>()
  {
    // Initialize variables
    func = Test_Default2D;

    // Negative interval
    StartPointX = -125.00;
    EndPointX = -89.0;
    StartPointY = -125.00;
    EndPointY = -89.0;

    TheoreticResult = Test_Default2D_Integral(StartPointX,EndPointX,StartPointY,EndPointY);
    // Numerical Integration 
    IntegrateFunctionNumerically();

    // Check solution against theoretical result
    CheckSolution(TheoreticResult);
  }

  /* Test 2:*/
  template<>
  template<>
  void AGQ2D_Object::test<2>()
  {
    // Initialize variables
    func = Test_Default2D;

    // Negative interval and Zero
    StartPointX = -30.00;
    EndPointX = 0.0;
    StartPointY = -125.00;
    EndPointY = -89.0;

    TheoreticResult = Test_Default2D_Integral(StartPointX,EndPointX,StartPointY,EndPointY);
    // Numerical Integration 
    IntegrateFunctionNumerically();

    // Check solution against theoretical result
    CheckSolution(TheoreticResult);

  }

  /* Test 3:*/
  template<>
  template<>
  void AGQ2D_Object::test<3>()
  {
    // Initialize variables
    func = Test_Default2D;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 1.0E-14;
    StartPointY = -25.00;
    EndPointY = 0.0;

    TheoreticResult = Test_Default2D_Integral(StartPointX,EndPointX,StartPointY,EndPointY);
    // Numerical Integration 
    IntegrateFunctionNumerically();

    // Check solution against theoretical result
    CheckSolution(TheoreticResult);

  }

  /* Test 4:*/
  template<>
  template<>
  void AGQ2D_Object::test<4>()
  {
    // Initialize variables
    func = Test_Default2D;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 1.0E-14;
    // Interval close to zero
    StartPointY = -1.15E-14;
    EndPointY = 0.2345;

    TheoreticResult = Test_Default2D_Integral(StartPointX,EndPointX,StartPointY,EndPointY);
    // Numerical Integration 
    IntegrateFunctionNumerically();

    // Check solution against theoretical result
    CheckSolution(TheoreticResult);

  }

  /* Test 5:*/
  template<>
  template<>
  void AGQ2D_Object::test<5>()
  {
    // Initialize variables
    func = Test_Default2D;

    // Positive interval
    StartPointX = 89.00;
    EndPointX = 125.0;
    StartPointY = 1250.00;
    EndPointY = 8900.0;

    TheoreticResult = Test_Default2D_Integral(StartPointX,EndPointX,StartPointY,EndPointY);
    // Numerical Integration 
    IntegrateFunctionNumerically();

    // Check solution against theoretical result
    CheckSolution(TheoreticResult);

  }

  /* Test 6:*/
  template<>
  template<>
  void AGQ2D_Object::test<6>()
  {
    // Initialize variables
    func = Test_Example2;
    tol = 1.0e-10;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 1.0E-14;
    // Interval close to zero
    StartPointY = -1.15E-14;
    EndPointY = 0.24;

    TheoreticResult = Test_Example2_Integral(StartPointX,EndPointX,StartPointY,EndPointY);
    // Numerical Integration 
    IntegrateFunctionNumerically();

    // Check solution against theoretical result
    CheckSolution(TheoreticResult);

  }

  /* Test 7:*/
  template<>
  template<>
  void AGQ2D_Object::test<7>()
  {
    // Initialize variables
    func = Test_Example2;
    tol = 1.0e-11;

    // Interval close to zero
    StartPointX = -16;
    EndPointX = 2.0;
    // Interval close to zero
    StartPointY = -1.3;
    EndPointY = 20;

    TheoreticResult = Test_Example2_Integral(StartPointX,EndPointX,StartPointY,EndPointY);
    // Numerical Integration 
    IntegrateFunctionNumerically();

    // Check solution against theoretical result
    CheckSolution(TheoreticResult);

  }

  /* Test 8:*/
  template<>
  template<>
  void AGQ2D_Object::test<8>()
  {
    // Initialize variables
    func = Test_Example3;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 1.0E-14;
    // Interval close to zero
    StartPointY = -1.15E-14;
    EndPointY = 0.24;

    TheoreticResult = Test_Example3_Integral(StartPointX,EndPointX,StartPointY,EndPointY);
    // Numerical Integration 
    IntegrateFunctionNumerically();

    // Check solution against theoretical result
    CheckSolution(TheoreticResult);

  }

  /* Test 9:*/
  template<>
  template<>
  void AGQ2D_Object::test<9>()
  {
    // Initialize variables
    func = Test_Example6;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 2.0;
    // Interval close to zero
    StartPointY = -1.15E-14;
    EndPointY = 0.24;

    TheoreticResult = Test_Example6_Integral(StartPointX,EndPointX,StartPointY,EndPointY);
    // Numerical Integration 
    IntegrateFunctionNumerically();

    // Check solution against theoretical result
    CheckSolution(TheoreticResult);

  }

  /* Test 10:*/
  template<>
  template<>
  void AGQ2D_Object::test<10>()
  {
    // Initialize variables
    func = Test_Example7;
    tol = 1.0e-8;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 2.0;
    // Interval close to zero
    StartPointY = -1.15E-14;
    EndPointY = 0.24;

    TheoreticResult = Test_Example7_Integral(StartPointX,EndPointX,StartPointY,EndPointY);
    // Numerical Integration 
    IntegrateFunctionNumerically();

    // Check solution against theoretical result
    CheckSolution(TheoreticResult);
  }

  /* Test 11:*/
  template<>
  template<>
  void AGQ2D_Object::test<11>()
  {
    // Initialize variables
    func = Test_Example11;

    // Interval close to zero
    StartPointX = -1.6;
    EndPointX = 2.4;
    // Interval close to zero
    StartPointY = -1.15E-14;
    EndPointY = 0.24;

    MapleResult = 0.0181317415192538;

    // Numerical Integration 
    IntegrateFunctionNumerically();

    // Check solution against Maple result
    CheckSolution(MapleResult);

  }

}

namespace tut

{
  // Data used for testing
  class Data_quad2dAdaptiveGaussianQuadrature: public TestData{
  public:
    FunctionType2D func;        // Use 1D function for integration
    double TheoreticResult;	// Result obtained with the exact integration
    double NumericResultA, NumericResultB, NumericResult;
    double StartPointX;		// One end of the interval (X dir)
    double EndPointX;		// The other end of the interval (X dir)
    double StartPointY;		// One end of the interval (Y dir)
    double EndPointY;		// The other end of the interval (Y dir)
    double AcceptError;	        // Acceptable error
    double MapleResult;		// Result obtained with Maple

    void IntegrateFunctionNumerically(void);
    void CheckSolution(double & Result);
  };

  void Data_quad2dAdaptiveGaussianQuadrature::IntegrateFunctionNumerically(void){
    NumericResult = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,digits);
    // Switch X limits
    NumericResultA = quad2dAdaptiveGaussianQuadrature(func,EndPointX,StartPointX,StartPointY,EndPointY,digits);
    // Switch Y limits
    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,EndPointY,StartPointY,digits);

    // check independency of limits
    AcceptError = tol*fabs(NumericResult);
    ensure_distance("AGQ2D: Change limits X", NumericResult, -NumericResultA, AcceptError);
    ensure_distance("AGQ2D: Change limits Y", NumericResult, -NumericResultB, AcceptError);
  }

  void Data_quad2dAdaptiveGaussianQuadrature::CheckSolution(double & Result){
    AcceptError = tol*fabs(Result);
    ensure_distance("AGQ2D: Integral Value", NumericResult, Result, AcceptError);
  }

  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (LinearSystems_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_quad2dAdaptiveGaussianQuadrature> quad2dAGQ_TestSuite;
  typedef quad2dAGQ_TestSuite::object  quad2dAGQ_Object;

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

  /* Test 1:*/
  template<>
  template<>
  void quad2dAGQ_Object::test<1>()
  {
    // Initialize variables
    func = Test_Default2D;

    // Negative interval
    StartPointX = -125.00;
    EndPointX = -89.0;
    StartPointY = -125.00;
    EndPointY = -89.0;

    // Analytic Integration
    TheoreticResult = Test_Default2D_Integral(StartPointX,EndPointX,StartPointY,EndPointY);

    // Numeric Integration
    IntegrateFunctionNumerically();

    // Check against anlytic result
    CheckSolution(TheoreticResult);
  }

  /* Test 2:*/
  template<>
  template<>
  void quad2dAGQ_Object::test<2>()
  {
    // Initialize variables
    func = Test_Default2D;

    // Negative interval and Zero
    StartPointX = -30.00;
    EndPointX = 0.0;
    StartPointY = -125.00;
    EndPointY = -89.0;

    TheoreticResult = Test_Default2D_Integral(StartPointX,EndPointX,StartPointY,EndPointY);

    // Numeric Integration
    IntegrateFunctionNumerically();

    // Check against anlytic result
    CheckSolution(TheoreticResult);

  }

  /* Test 3:*/
  template<>
  template<>
  void quad2dAGQ_Object::test<3>()
  {
    // Initialize variables
    func = Test_Default2D;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 1.0E-14;
    StartPointY = -25.00;
    EndPointY = 0.0;

    TheoreticResult = Test_Default2D_Integral(StartPointX,EndPointX,StartPointY,EndPointY);

    // Numeric Integration
    IntegrateFunctionNumerically();

    // Check against anlytic result
    CheckSolution(TheoreticResult);

  }

  /* Test 4:*/
  template<>
  template<>
  void quad2dAGQ_Object::test<4>()
  {
    // Initialize variables
    func = Test_Default2D;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 1.0E-14;
    // Interval close to zero
    StartPointY = -1.15E-14;
    EndPointY = 0.2345;

    TheoreticResult = Test_Default2D_Integral(StartPointX,EndPointX,StartPointY,EndPointY);

    // Numeric Integration
    IntegrateFunctionNumerically();

    // Check against anlytic result
    CheckSolution(TheoreticResult);

  }

  /* Test 5:*/
  template<>
  template<>
  void quad2dAGQ_Object::test<5>()
  {
    // Initialize variables
    func = Test_Default2D;

    // Positive interval
    StartPointX = 89.00;
    EndPointX = 125.0;
    StartPointY = 1250.00;
    EndPointY = 8900.0;

    TheoreticResult = Test_Default2D_Integral(StartPointX,EndPointX,StartPointY,EndPointY);

    // Numeric Integration
    IntegrateFunctionNumerically();

    // Check against anlytic result
    CheckSolution(TheoreticResult);

  }

  /* Test 6:*/
  template<>
  template<>
  void quad2dAGQ_Object::test<6>()
  {
    // Initialize variables
    func = Test_Example2;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 1.0E-14;
    // Interval close to zero
    StartPointY = -1.15E-14;
    EndPointY = 0.24;

    TheoreticResult = Test_Example2_Integral(StartPointX,EndPointX,StartPointY,EndPointY);

    // Numeric Integration
    IntegrateFunctionNumerically();

    // Check against anlytic result
    tol = 1.0e-10;
    CheckSolution(TheoreticResult);

  }

  /* Test 7:*/
  template<>
  template<>
  void quad2dAGQ_Object::test<7>()
  {
    // Initialize variables
    func = Test_Example2;

    // Interval close to zero
    StartPointX = -16;
    EndPointX = 2.0;
    // Interval close to zero
    StartPointY = -1.3;
    EndPointY = 20;

    TheoreticResult = Test_Example2_Integral(StartPointX,EndPointX,StartPointY,EndPointY);

    // Numeric Integration
    IntegrateFunctionNumerically();

    // Check against anlytic result
    tol = 1.0e-10;
    CheckSolution(TheoreticResult);

  }

  /* Test 8:*/
  template<>
  template<>
  void quad2dAGQ_Object::test<8>()
  {
    // Initialize variables
    func = Test_Example3;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 1.0E-14;
    // Interval close to zero
    StartPointY = -1.15E-14;
    EndPointY = 0.24;

    TheoreticResult = Test_Example3_Integral(StartPointX,EndPointX,StartPointY,EndPointY);

    // Numeric Integration
    IntegrateFunctionNumerically();

    // Check against anlytic result
    CheckSolution(TheoreticResult);

  }

  /* Test 9:*/
  template<>
  template<>
  void quad2dAGQ_Object::test<9>()
  {
    // Initialize variables
    func = Test_Example6;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 2.0;
    // Interval close to zero
    StartPointY = -1.15E-14;
    EndPointY = 0.24;

    TheoreticResult = Test_Example6_Integral(StartPointX,EndPointX,StartPointY,EndPointY);

    // Numeric Integration
    IntegrateFunctionNumerically();

    // Check against anlytic result
    CheckSolution(TheoreticResult);

  }

  /* Test 10:*/
  template<>
  template<>
  void quad2dAGQ_Object::test<10>()
  {
    // Initialize variables
    func = Test_Example7;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 2.0;
    // Interval close to zero
    StartPointY = -1.15E-14;
    EndPointY = 0.24;

    TheoreticResult = Test_Example7_Integral(StartPointX,EndPointX,StartPointY,EndPointY);

    // Numeric Integration
    IntegrateFunctionNumerically();

    // Check against anlytic result
    tol = 1.0e-8;
    CheckSolution(TheoreticResult);

  }


  /* Test 11:*/
  template<>
  template<>
  void quad2dAGQ_Object::test<11>()
  {
    // Initialize variables
    func = Test_Example11;

    // Interval close to zero
    StartPointX = -1.6;
    EndPointX = 2.4;
    // Interval close to zero
    StartPointY = -1.15E-14;
    EndPointY = 0.24;

    MapleResult = 0.0181317415192538;

    // Numeric Integration
    IntegrateFunctionNumerically();

    // Check against anlytic result
    CheckSolution(MapleResult);
  }

}

namespace tut
{
  // Data used for testing
  class Data_AdaptiveGaussianQuadrature3D: public TestData{
  public:
    FunctionType3D func;        // Use 3D function for integration
    double TheoreticResult;	// Result obtained with the exact integration
    double NumericResultA,      // Result for limits changed
           NumericResultB,      // Result for limits changed
           NumericResultC,      // Result for limits changed
           NumericResult;	// Main result
    double StartPointX;		// One end of the interval (X dir)
    double EndPointX;		// The other end of the interval (X dir)
    double StartPointY;		// One end of the interval (Y dir)
    double EndPointY;		// The other end of the interval (Y dir)
    double StartPointZ;		// One end of the interval (Z dir)
    double EndPointZ;		// The other end of the interval (Z dir)
    double AcceptError;	        // Acceptable error
    double MapleResult;		// Result obtained with Maple


    void IntegrateFunctionNumerically(void);
    void CheckSolution(double & Result);

  };

  void Data_AdaptiveGaussianQuadrature3D::IntegrateFunctionNumerically(void){
    NumericResult = AdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,
					       StartPointZ,EndPointZ, digits,NumericResult);

    // Switch X limits
    NumericResultA = AdaptiveGaussianQuadrature(func,EndPointX,StartPointX,StartPointY,EndPointY,
						StartPointZ,EndPointZ, digits,NumericResult);
    
    // Switch Y limits
    NumericResultB = AdaptiveGaussianQuadrature(func,StartPointX,EndPointX,EndPointY,StartPointY,
						StartPointZ,EndPointZ, digits,NumericResult);
    
    // Switch Z limits
    NumericResultC = AdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,
						EndPointZ,StartPointZ,digits,NumericResult);

    // check independency of limits
    AcceptError = tol*fabs(NumericResult);
    ensure_distance("AGQ3D: Change limits X", NumericResult, -NumericResultA, AcceptError);
    ensure_distance("AGQ3D: Change limits Y", NumericResult, -NumericResultB, AcceptError);
    ensure_distance("AGQ3D: Change limits Z", NumericResult, -NumericResultC, AcceptError);
  }

  void Data_AdaptiveGaussianQuadrature3D::CheckSolution(double & Result){
    AcceptError = tol*fabs(Result);
    ensure_distance("AGQ3D: Integral Value", NumericResult, Result, AcceptError);
  }

  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   */
  typedef test_group<Data_AdaptiveGaussianQuadrature3D> AGQ3D_TestSuite;
  typedef AGQ3D_TestSuite::object  AGQ3D_Object;

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

  /* Test 1:*/
  template<>
  template<>
  void AGQ3D_Object::test<1>()
  {
    // Initialize variables
    func = Test_Default3D;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 2.0;
    // Interval close to zero
    StartPointY = -1.0;
    EndPointY = 0.24;

    StartPointZ = -1.34;
    EndPointZ = -0.04;

    // Analytical Integration
    TheoreticResult = Test_Default3D_Integral(StartPointX,EndPointX,StartPointY,EndPointY,
					      StartPointZ,EndPointZ);

    // Numeric Integration
    IntegrateFunctionNumerically();

    // Check
    CheckSolution(TheoreticResult);
  }

}

namespace tut
{
  // 1D
  GLQ_TestSuite GLQ1D_Test("Integration:GaussLobattoQuadrature() in 1D");
  AGQ_TestSuite AGQ1D_Test("Integration:AdaptiveGaussianQuadrature() in 1D");

  // 2D
  AGQ2D_TestSuite AGQ2D_Test("Integration:AdaptiveGaussianQuadrature() in 2D");
  quad2dAGQ_TestSuite quad2dAGQ_Test("Integration:quad2DAdaptiveGaussianQuadrature() in 2D");

  // 3D
  AGQ3D_TestSuite AGQ3D_Test("Integration:AdaptiveGaussianQuadrature() in 3D");

}
