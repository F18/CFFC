/*!\file test_NumericalLibrary.cc
   \brief Regression tests for functions prototyped in NumericalLibrary.h */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */

/* Include CFFC header files */
#include "TestData.h"
#include "../../../../src_2D/Math//NumericalLibrary.h"
#include "TestFunctions/TestFunctions_1D.h"
#include "TestFunctions/TestFunctions_2D.h"
#include "TestFunctions/TestFunctions_3D.h"


namespace tut
{

  // **********************************************
  //                 TEST SUITE: GLQ_TestSuite
  // **********************************************
  // Data used for testing
  // **********************
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
    func = Test_Example20;
    StartPoint = 0.0;
    EndPoint = 1.0;

    TheoreticResult = Test_Example20_Integral(StartPoint,EndPoint);
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


  // **********************************************
  //                 TEST SUITE: AGQ_TestSuite
  // **********************************************
  // Data used for testing
  // **********************
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
    FunctionType1D func2 = Test_Example6;
    ErrorFunc<double,FunctionType1D> Error1D = error_function(func,func2,NumericResultA);
    MatlabResult = 9.639642857142857;

    // Pozitive interval
    StartPoint = -2.0;
    EndPoint = -1.0;

    NumericResultA = AdaptiveGaussianQuadrature(Error1D,StartPoint,EndPoint,DummyParam,digits);

    // Check
    ensure_distance("ADG: Matlab Value", NumericResultA, MatlabResult, tol*MatlabResult);
  }

  /* Test 23:*/
  template<>
  template<>
  void AGQ_Object::test<23>()
  {
    // Initialize variables
    func = Test_Example19;
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
    func = Test_Example20;
    StartPoint = 0.0;
    EndPoint = 1.0;

    TheoreticResult = Test_Example20_Integral(StartPoint,EndPoint);
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




  // **********************************************
  //                 TEST SUITE: AGQ2_TestSuite
  // **********************************************
  // Data used for testing
  // **********************
  class Data_AdaptiveGaussianQuadrature2D: public TestData{
  public:
    FunctionType2D func;        // Use 2D function for integration
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

  /* Test 12:*/
  template<>
  template<>
  void AGQ2D_Object::test<12>()
  {
    // Initialize variables
    func = Test_Example11;
    FunctionType2D func2 = Test_Example6;
    ErrorFunc<double,FunctionType2D> Error2D(func,func2);

    // Interval close to zero
    StartPointX = -1.6;
    EndPointX = 2.4;
    // Interval close to zero
    StartPointY = -1.15E-14;
    EndPointY = 0.24;

    MapleResult = 0.963948258480793;

    // Numerical Integration 
    NumericResult = AdaptiveGaussianQuadrature(Error2D,StartPointX,EndPointX,StartPointY,EndPointY,
					       digits,NumericResult);

    // Check solution against Maple result
    CheckSolution(MapleResult);

  }



  // *************************************************
  //                 TEST SUITE: quad2dAGQ_TestSuite
  // *************************************************
  // Data used for testing
  // **********************
  class Data_quad2dAdaptiveGaussianQuadrature: public TestData{
  public:
    FunctionType2D func;        // Use 2D function for integration
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




  // **********************************************
  //                 TEST SUITE: AGQ3D_TestSuite
  // **********************************************
  // Data used for testing
  // **********************
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



  // ******************************************************
  //                 TEST SUITE: ErrorSubroutines_TestSuite
  // ******************************************************
  // Data used for testing
  // **********************
  class Data_ErrorSubroutines: public TestData{
  public:
    double TheoreticResult;	// Result obtained with the exact integration
    double NumericResultA, NumericResultB, NumericResult;
    double StartPointX;		// One end of the interval (X dir)
    double EndPointX;		// The other end of the interval (X dir)
    double StartPointY;		// One end of the interval (Y dir)
    double EndPointY;		// The other end of the interval (Y dir)
    double AcceptError;	        // Acceptable error
    double MapleResult;		// Result obtained with Maple
    double DummyParam;


    class TestClass{
    private:
      TestClass();
      double val;

    public:
      typedef double (TestClass::*MemberFunctionType1D)(const double &);
      typedef double (TestClass::*MemberFunctionType2D)(const double &, const double &);
      typedef double (TestClass::*MemberFunctionType3D)(const double &, const double &, const double &);

      TestClass(double _val_): val(_val_){}

      // member functions
      double SolutionAtCoordinate(const double &coord){return val*coord;}
      double SquareSolutionAtCoordinate(const double &coord){return val*val*coord*coord;}

      double SolutionAtCoordinate(const double &coord, const double & coord2){return val*coord*coord2;}
      double SquareSolutionAtCoordinate(const double &coord, const double &coord2){return val*val*coord2*coord2;}
      
      double SolutionAtCoordinate(const double &coord, const double &coord2, const double &coord3){
	return val*coord*coord2*coord3;
      }
      double SquareSolutionAtCoordinate(const double &coord, const double &coord2, const double &coord3){
	return val*val*coord2*coord2*coord3*coord3;
      }
    };

  };

  typedef test_group<Data_ErrorSubroutines> ErrorSubroutines_TestSuite;
  typedef ErrorSubroutines_TestSuite::object  ErrorSub_Object;

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
  void ErrorSub_Object::test<1>()
  {
    set_test_name("check _Mapped_Function_Wrapper_");

    FunctionType1D func = Test_Example1;
    _Mapped_Function_Wrapper_<FunctionType1D,double> MappedFunction(func,-100.0, -50.0, -1.0, 1.0);

    ensure_equals("Evaluate funct at the upper bound", MappedFunction(-50), func(1));
    ensure_equals("Evaluate funct at the lower bound", MappedFunction(-100), func(-1));
    ensure_equals("Evaluate funct in the middle", MappedFunction(-75), func(0));
  }

  /* Test 2:*/
  template<>
  template<>
  void ErrorSub_Object::test<2>()
  {
    set_test_name("check _Member_Function_Wrapper_ in 1D");

    TestClass Obj(3.5);		// create object

    // First wrapper
    _Member_Function_Wrapper_<TestClass,
      TestClass::MemberFunctionType1D,double> WrappedMemberFunction(&Obj,
								    &TestClass::SolutionAtCoordinate);

    // Second wrapper
    _Member_Function_Wrapper_<TestClass,
      TestClass::MemberFunctionType1D,double> WrappedMemberFunctionSquare(&Obj,
									  &TestClass::SquareSolutionAtCoordinate);

    ensure_equals("Evaluate SolutionAtCoordinate() at 6.0", WrappedMemberFunction(6.0), Obj.SolutionAtCoordinate(6.0));
    ensure_equals("Evaluate SquareSolutionAtCoordinate() at 6.0",
		  WrappedMemberFunctionSquare(6.0),
		  Obj.SquareSolutionAtCoordinate(6.0));
  }

  /* Test 3:*/
  template<>
  template<>
  void ErrorSub_Object::test<3>()
  {
    set_test_name("check _Member_Function_Wrapper_ in 2D");

    TestClass Obj(5.55);		// create object

    // First wrapper
    _Member_Function_Wrapper_<TestClass,
      TestClass::MemberFunctionType2D,double> WrappedMemberFunction(&Obj,
								    &TestClass::SolutionAtCoordinate);

    // Second wrapper
    _Member_Function_Wrapper_<TestClass,
      TestClass::MemberFunctionType2D,double> WrappedMemberFunctionSquare(&Obj,
									  &TestClass::SquareSolutionAtCoordinate);

    ensure_equals("Evaluate SolutionAtCoordinate() at (6.0,1.23)",
		  WrappedMemberFunction(6.0,1.23),
		  Obj.SolutionAtCoordinate(6.0,1.23));
    ensure_equals("Evaluate SquareSolutionAtCoordinate() at (6.0,1.23)",
		  WrappedMemberFunctionSquare(6.0,1.23),
		  Obj.SquareSolutionAtCoordinate(6.0,1.23));
  }

  /* Test 4:*/
  template<>
  template<>
  void ErrorSub_Object::test<4>()
  {
    set_test_name("check _Member_Function_Wrapper_ in 3D");

    TestClass Obj(-83.5);		// create object

    // First wrapper
    _Member_Function_Wrapper_<TestClass,
      TestClass::MemberFunctionType3D,double> WrappedMemberFunction(&Obj,
								    &TestClass::SolutionAtCoordinate);

    // Second wrapper
    _Member_Function_Wrapper_<TestClass,
      TestClass::MemberFunctionType3D,double> WrappedMemberFunctionSquare(&Obj,
									  &TestClass::SquareSolutionAtCoordinate);

    ensure_equals("Evaluate SolutionAtCoordinate() at (1.1,2.2,-3.3)",
		  WrappedMemberFunction(1.1,2.2,-3.3),
		  Obj.SolutionAtCoordinate(1.1,2.2,-3.3));
    ensure_equals("Evaluate SquareSolutionAtCoordinate() at (1.1,2.2,-3.3)",
		  WrappedMemberFunctionSquare(1.1,2.2,-3.3),
		  Obj.SquareSolutionAtCoordinate(1.1,2.2,-3.3));
  }

  /* Test 5:*/
  template<>
  template<>
  void ErrorSub_Object::test<5>()
  {
    set_test_name("integrate with _Member_Function_Wrapper_ in 1D");

    TestClass Obj(-3.5);		// create object

    // First wrapper
    _Member_Function_Wrapper_<TestClass,
      TestClass::MemberFunctionType1D,double> WrappedMemberFunction(&Obj,
								    &TestClass::SolutionAtCoordinate);

    // Second wrapper
    _Member_Function_Wrapper_<TestClass,
      TestClass::MemberFunctionType1D,double> WrappedMemberFunctionSquare(&Obj,
									  &TestClass::SquareSolutionAtCoordinate);

    // Integrate wrapper

    // Define interval
    StartPointX = -5.00;
    EndPointX = 0.0;

    // AdaptiveGaussianQuadrature
    NumericResultA = AdaptiveGaussianQuadrature(WrappedMemberFunction,StartPointX,EndPointX,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(WrappedMemberFunctionSquare,StartPointX,EndPointX,DummyParam);

    ensure_equals("NumericResultA", NumericResultA, 43.75);
    ensure_distance("NumericResultB", NumericResultB, 510.4166666666667, 1.0e-13);
  }

  /* Test 6:*/
  template<>
  template<>
  void ErrorSub_Object::test<6>()
  {
    set_test_name("integrate error with _Member_Function_Wrapper_ in 1D");

    TestClass Obj(-3.5);		// create object

    // First wrapper
    _Member_Function_Wrapper_<TestClass,
      TestClass::MemberFunctionType1D,double> WrappedMemberFunction(&Obj,
								    &TestClass::SolutionAtCoordinate);

    // Second wrapper
    _Member_Function_Wrapper_<TestClass,
      TestClass::MemberFunctionType1D,double> WrappedMemberFunctionSquare(&Obj,
									  &TestClass::SquareSolutionAtCoordinate);

    ErrorFunc<double,FunctionType1D,_Member_Function_Wrapper_<TestClass,TestClass::MemberFunctionType1D,double> >
      error(Test_Example1,
	    WrappedMemberFunction);

    // Integrate wrapper

    // Define interval
    StartPointX = -5.00;
    EndPointX = 0.0;

    // AdaptiveGaussianQuadrature

    // Use error_function
    NumericResultA = AdaptiveGaussianQuadrature(error_function(WrappedMemberFunction,
							       WrappedMemberFunctionSquare,
							       NumericResultA),
					       StartPointX,EndPointX,DummyParam);
    
    // Use error object
    NumericResultB = AdaptiveGaussianQuadrature(error,StartPointX,EndPointX,DummyParam);

    ensure_distance("NumericResultA", NumericResultA, 466.7619047619047619, 1.0e-13);
    ensure_distance("NumericResultB", NumericResultB, 47.2950836148489, 1.0e-13);
  }


  /* Test 7:*/
  template<>
  template<>
  void ErrorSub_Object::test<7>()
  {
    set_test_name("error_function with wrapped_member_function and mapped_function");

    TestClass Obj(-3.5);		// create object

    FunctionType1D func = Test_Example1;
    TestClass::MemberFunctionType1D mem_func = &TestClass::SquareSolutionAtCoordinate;
    double Coord = -2.5;
    TheoreticResult = 74.5625;

    ensure_distance("error_function",
		    error_function(mapped_function(func,NumericResult,-5.0,0.0,-1.0,1.0), // first function
				   wrapped_member_function(&Obj,                          // second function
							   mem_func,
							   NumericResult), 
				   NumericResult)(Coord),                                 // solution type and the X coordinate 
		    TheoreticResult,
		    1.0e-13);
  }

  /* Test 8:*/
  template<>
  template<>
  void ErrorSub_Object::test<8>()
  {
    set_test_name("integrate error with _Member_Function_Wrapper_ and _Mapped_Function_Wrapper_ in 1D");

    TestClass Obj(-3.5);		// create object

    // Wrap member function
    _Member_Function_Wrapper_<TestClass,
      TestClass::MemberFunctionType1D,double> WrappedMemberFunctionSquare(&Obj,
									  &TestClass::SquareSolutionAtCoordinate);

    FunctionType1D func = Test_Example1;

    // Integrate wrapper

    // Define interval
    StartPointX = -5.00;
    EndPointX = 0.0;

    NumericResultA = AdaptiveGaussianQuadrature(error_function(mapped_function(func,NumericResult,-5.0,0.0,-1.0,1.0),
							       WrappedMemberFunctionSquare,
							       NumericResult),
						StartPointX,EndPointX,DummyParam);
    
    ensure_distance("Result", NumericResultA, 503.2298836292962, 1.0e-12);
  }

  /* Test 9:*/
  template<>
  template<>
  void ErrorSub_Object::test<9>()
  {
    set_test_name("integrate with wrapped_member_function");

    TestClass Obj(-3.5);		// create object

    // Integrate member function

    // Define interval
    StartPointX = -5.00;
    EndPointX = 0.0;

    TestClass::MemberFunctionType1D func = &TestClass::SolutionAtCoordinate;

    // Use wrapped_member_function
    NumericResultA = AdaptiveGaussianQuadrature(wrapped_member_function(&Obj,
									func,
									NumericResult),
						StartPointX,EndPointX,DummyParam);
    
    // Test result
    ensure_distance("NumericResultA", NumericResultA, 43.75, 1.0e-13);
  }

  /* Test 10:*/
  template<>
  template<>
  void ErrorSub_Object::test<10>()
  {
    set_test_name("integrate square_error_function with wrapped_memeber_function in 1D");

    TestClass Obj(-3.5);		// create object

    TestClass::MemberFunctionType1D func = &TestClass::SolutionAtCoordinate;
    FunctionType1D func1D = Test_Example1;

    // Define interval
    StartPointX = -5.00;
    EndPointX = 0.0;

    // Use square_error_function
    NumericResultA = AdaptiveGaussianQuadrature(square_error_function(func1D,
								      wrapped_member_function(&Obj,func,NumericResultA),
								      NumericResultA),
						StartPointX,EndPointX,DummyParam);

    // Test result
    ensure_distance("Result", NumericResultA, 610.53786993983962760, 1.0e-13);
  }

}


/************************
 *  CREATE TEST SUITES  *
 ************************/

// 1D
tut::GLQ_TestSuite GLQ1D_Test("Integration:GaussLobattoQuadrature() in 1D");
tut::AGQ_TestSuite AGQ1D_Test("Integration:AdaptiveGaussianQuadrature() in 1D");

// 2D
tut::AGQ2D_TestSuite AGQ2D_Test("Integration:AdaptiveGaussianQuadrature() in 2D");
tut::quad2dAGQ_TestSuite quad2dAGQ_Test("Integration:quad2DAdaptiveGaussianQuadrature() in 2D");

// 3D
tut::AGQ3D_TestSuite AGQ3D_Test("Integration:AdaptiveGaussianQuadrature() in 3D");

// ErrorSubroutines
tut::ErrorSubroutines_TestSuite ErrorSubroutines_Test("Numerical Library:Error functors & wrappers");
