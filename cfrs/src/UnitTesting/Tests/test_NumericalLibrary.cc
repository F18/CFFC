// Test file for functions used in Numerical Library
#include <tut.h>
#include "Math/NumericalLibrary.h"
#include "TestFunctions/TestFunctions.h"
#include "TestFunctions/TestFunctions_2D.h"
#include "TestFunctions/TestFunctions_3D.h"
#include <cmath>
#include <limits>

namespace tut

{
  // Data used for testing
  struct Data_AdaptiveGaussianQuadrature{
    int digits;			// Number of exact imposed digits --> gives the precision
    long double tol;	                // Numerical tolerance for comparing the solutions
    FunctionType1D func;        // Use 1D function for integration
    long double TheoreticResult;    // Result obtained with the exact integration
    long double NumericResultA, NumericResultB;
    double StartPoint;		// One end of the interval
    double EndPoint;		// The other end of the interval
    long double RelErr;		// Relative error
    long double EpsMachine;		// Machine accuracy
    double MatlabResult;	// The result obtained using Matlab
    long double DummyParam;
 };

   /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (LinearSystems_TestGroup) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_AdaptiveGaussianQuadrature> AdaptiveGaussianQuadrature_TestGroup;
  typedef AdaptiveGaussianQuadrature_TestGroup::object
                                             AdaptiveGaussianQuadrature_TestObject;
}

/**********************************************
 *                THE TESTS                   *
 *********************************************/

namespace tut
{

  /* Test 1:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<1>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Default1D;
    MatlabResult = 6.178088507844034e+13;

    // Negative interval
    StartPoint = -125.00;
    EndPoint = -89.0;

    TheoreticResult = Test_Default1D_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);

    if ( fabs(TheoreticResult) > EpsMachine){
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    }
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

     // Check results
    ensure("Check integration", RelErr <= tol);
    if ( fabs(NumericResultA) > EpsMachine){
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    }
    else
      RelErr = fabs(NumericResultA + NumericResultB);
    ensure("Check independency of limits", RelErr <= tol);

    if ( fabs(MatlabResult) > EpsMachine)
      RelErr = fabs((MatlabResult - NumericResultA)/MatlabResult);
    else
      RelErr = fabs(MatlabResult - NumericResultA);
    ensure("Check Matlab Result", RelErr <= tol);
  }

  /* Test 2:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<2>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    EpsMachine = numeric_limits<long double>::epsilon( );
    tol = 0.5*pow(10.0,1.0-digits);
    func = Test_Default1D;

    // Negative interval and zero
    StartPoint = -2500.00;
    EndPoint = 0.0;

    TheoreticResult = Test_Default1D_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);
    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    if ( fabs(NumericResultA) > EpsMachine)
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    else
      RelErr = fabs(NumericResultA + NumericResultB);
    ensure("Check independency of limits", RelErr <= tol);
  }

  /* Test 3:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<3>()
  {
    // Initialize variables
    digits = 15;
    EpsMachine = 1.0E-15;
    tol = 0.5*pow(10.0,1.0-digits);
    func = Test_Default1D;

    // Values close to zero
    StartPoint = -0.00000016;
    EndPoint = 1.0E-14;

    TheoreticResult = Test_Default1D_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    if ( fabs(NumericResultA) > EpsMachine)
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    else
      RelErr = fabs(NumericResultA + NumericResultB);
    ensure("Check independency of limits", RelErr <= tol);
  }

  /* Test 4:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<4>()
  {
    // Initialize variables
    digits = 15;
    EpsMachine = 1.0E-15;
    tol = 0.5*pow(10.0,1.0-digits);
    func = Test_Default1D;

    // Values close to zero
    StartPoint = -1.15E-14;
    EndPoint = 0.2345;

    TheoreticResult = Test_Default1D_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);
    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    if ( fabs(NumericResultA) > EpsMachine)
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    else
      RelErr = fabs(NumericResultA + NumericResultB);
    ensure("Check independency of limits", RelErr <= tol);
  }

  /* Test 5:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<5>()
  {
    // Initialize variables
    digits = 15;
    EpsMachine = 1.0E-15;
    tol = 0.5*pow(10.0,1.0-digits);
    func = Test_Default1D;

    // Pozitive interval
    StartPoint = 125.00;
    EndPoint = 89.0;

    TheoreticResult = Test_Default1D_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);
    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    if ( fabs(NumericResultA) > EpsMachine)
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    else
      RelErr = fabs(NumericResultA + NumericResultB);
    ensure("Check independency of limits", RelErr <= tol);
  }

  /* Test 6:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<6>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example1;

    // Negative interval
    StartPoint = -125.00;
    EndPoint = -89.0;

    TheoreticResult = Test_Example1_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

   if ( fabs(NumericResultA) > EpsMachine)
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    else
      RelErr = fabs(NumericResultA + NumericResultB);
    ensure("Check independency of limits", RelErr <= tol);
  }

  /* Test 7:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<7>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example1;

    // Negative interval and Zero
    StartPoint = -2500.00;
    EndPoint = 0.0;

    TheoreticResult = Test_Example1_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);
    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    if ( fabs(NumericResultA) > EpsMachine){
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    }
    else
      RelErr = fabs(NumericResultA + NumericResultB);

    ensure("Check independency of limits", RelErr <= tol);
  }

  /* Test 8:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<8>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10-1;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example1;

    // Interval close to zero
    StartPoint = -0.00000016;
    EndPoint = 1.0E-14;

    TheoreticResult = Test_Example1_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

   if ( fabs(NumericResultA) > EpsMachine)
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    else
      RelErr = fabs(NumericResultA + NumericResultB);

    ensure("Check independency of limits", RelErr <= tol);

  }

  /* Test 9:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<9>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example1;

    // Interval close to zero
    StartPoint = -1.15E-14;
    EndPoint = 0.2345;

    TheoreticResult = Test_Example1_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);
    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    if ( fabs(NumericResultA) > EpsMachine)
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    else
      RelErr = fabs(NumericResultA + NumericResultB);

    ensure("Check independency of limits", RelErr <= tol);
  }

  /* Test 10:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<10>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example1;

    // Pozitive interval
    StartPoint = 125.00;
    EndPoint = 89.0;

    TheoreticResult = Test_Example1_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);
    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    if ( fabs(NumericResultA) > EpsMachine)
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    else
      RelErr = fabs(NumericResultA + NumericResultB);
    ensure("Check independency of limits", RelErr <= tol);
  }

 /* Test 11:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<11>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example3;

    // Negative interval
    StartPoint = -125.00;
    EndPoint = -89.0;

    TheoreticResult = Test_Example3_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

   if ( fabs(NumericResultA) > EpsMachine)
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    else
      RelErr = fabs(NumericResultA + NumericResultB);
    ensure("Check independency of limits", RelErr <= tol);
  }

  /* Test 12:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<12>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example3;

    // Negative interval and Zero
    StartPoint = -2500.00;
    EndPoint = 0.0;

    TheoreticResult = Test_Example3_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);
    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

   if ( fabs(NumericResultA) > EpsMachine)
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    else
      RelErr = fabs(NumericResultA + NumericResultB);

    ensure("Check independency of limits", RelErr <= tol);
  }

  /* Test 13:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<13>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example3;

    // Interval close to zero
    StartPoint = -0.00000016;
    EndPoint = 1.0E-14;

    TheoreticResult = Test_Example3_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

   if ( fabs(NumericResultA) > EpsMachine)
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    else
      RelErr = fabs(NumericResultA + NumericResultB);

    ensure("Check independency of limits", RelErr <= tol);

  }

  /* Test 14:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<14>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example3;

    // Interval close to zero
    StartPoint = -1.15E-14;
    EndPoint = 0.2345;

    TheoreticResult = Test_Example3_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);
    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    if ( fabs(NumericResultA) > EpsMachine)
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    else
      RelErr = fabs(NumericResultA + NumericResultB);

    ensure("Check independency of limits", RelErr <= tol);
  }

  /* Test 15:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<15>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example3;

    // Pozitive interval
    StartPoint = 125.00;
    EndPoint = 89.0;

    TheoreticResult = Test_Example3_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);
    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    if ( fabs(NumericResultA) > EpsMachine)
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    else
      RelErr = fabs(NumericResultA + NumericResultB);
    ensure("Check independency of limits", RelErr <= tol);
  }


 /* Test 16:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<16>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10-1;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example2;

    // Negative interval
    StartPoint = -125.00;
    EndPoint = -89.0;

    TheoreticResult = Test_Example2_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

   if ( fabs(NumericResultA) > EpsMachine)
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    else
      RelErr = fabs(NumericResultA + NumericResultB);
    ensure("Check independency of limits", RelErr <= tol);
  }

  /* Test 17:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<17>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example2;

    // Negative interval and Zero
    StartPoint = -30.00;
    EndPoint = 0.0;

    TheoreticResult = Test_Example2_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);
    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

   if ( fabs(NumericResultA) > EpsMachine)
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    else
      RelErr = fabs(NumericResultA + NumericResultB);

    ensure("Check independency of limits", RelErr <= tol);
  }

  /* Test 18:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<18>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = 7.0E-14;
    func = Test_Example2;

    // Interval close to zero
    StartPoint = -0.00000016;
    EndPoint = 1.0E-14;

    TheoreticResult = Test_Example2_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

   if ( fabs(NumericResultA) > EpsMachine)
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    else
      RelErr = fabs(NumericResultA + NumericResultB);

    ensure("Check independency of limits", RelErr <= tol);

  }

  /* Test 19:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<19>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example2;

    // Interval close to zero
    StartPoint = -1.15E-14;
    EndPoint = 0.2345;

    TheoreticResult = Test_Example2_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);
    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    if ( fabs(NumericResultA) > EpsMachine)
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    else
      RelErr = fabs(NumericResultA + NumericResultB);

    ensure("Check independency of limits", RelErr <= tol);
  }

  /* Test 20:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<20>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example2;

    // Pozitive interval
    StartPoint = 125.00;
    EndPoint = 89.0;

    TheoreticResult = Test_Example2_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);
    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    if ( fabs(NumericResultA) > EpsMachine)
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    else
      RelErr = fabs(NumericResultA + NumericResultB);
    ensure("Check independency of limits", RelErr <= tol);
  }

  /* Test 21:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<21>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example4;

    // Pozitive interval
    StartPoint = 125.00;
    EndPoint = 89.0;

    TheoreticResult = Test_Example4_Integral(StartPoint,EndPoint);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPoint,EndPoint,digits,DummyParam);
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPoint,StartPoint,digits,DummyParam);
    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    if ( fabs(NumericResultA) > EpsMachine)
      RelErr = fabs((NumericResultA + NumericResultB)/ NumericResultA);
    else
      RelErr = fabs(NumericResultA + NumericResultB);
    ensure("Check independency of limits", RelErr <= tol);
  }

  /* Test 22:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature_TestObject::test<22>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Default1D;
    SuperFunctionType1D super_func = Error1D;
    FunctionType1D func2 = Test_Example6;
    double MatlabResult = 9.639642857142857;
    // Pozitive interval
    StartPoint = -2.0;
    EndPoint = -1.0;

    TheoreticResult = MatlabResult;
    NumericResultA = AdaptiveGaussianQuadrature(super_func,func,func2,StartPoint,EndPoint,digits);
 
    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);
  }

//   /* Test 23:*/
//   template<>
//   template<>
//   void AdaptiveGaussianQuadrature_TestObject::test<23>()
//   {
//     // Initialize variables
//     typedef Gaussian2D_pState (*GaussianType) (const double);

//     digits = numeric_limits<double>::digits10;
//     tol = 0.5*pow(10.0,1.0-digits);
//     EpsMachine = numeric_limits<long double>::epsilon( );
//     GaussianType GaussianFunction = GaussianTestFunction;
//     Gaussian2D_pState Result, RelError;

//     // Pozitive interval
//     StartPoint = -2.0;
//     EndPoint = -1.0;

//     Gaussian2D_pState MatlabResult = Gaussian2D_pState(0.5*(pow(EndPoint,2) - pow(StartPoint,2)));

//     Result = AdaptiveGaussianQuadrature(GaussianFunction,StartPoint,EndPoint,digits,Result);
 
// //     if ( fabs(Result) > (Gaussian2D_pState)EpsMachine)
// //       RelError = fabs((MatlabResult - Result)/MatlabResult);
// //     else
//     RelError = fabs(MatlabResult - Result);

//     std::cout << "Integration Result: " << Result << std::endl; 
//     std::cout << "MatlabResult: " << MatlabResult << std::endl; 
//     std::cout << "RelError: " << RelError << std::endl; 
//     // Check results
//     //    ensure("Check integration", RelError < (Gaussian2D_pState)tol);
//   }


};


namespace tut

{
  // Data used for testing
  struct Data_AdaptiveGaussianQuadrature2D{
    int digits;			// Number of exact imposed digits --> gives the precision
    double tol;	                // Numerical tolerance for comparing the solutions
    FunctionType2D func;        // Use 1D function for integration
    double TheoreticResult;	// Result obtained with the exact integration
    double NumericResultA, NumericResultB, NumericResult;
    double StartPointX;		// One end of the interval (X dir)
    double EndPointX;		// The other end of the interval (X dir)
    double StartPointY;		// One end of the interval (Y dir)
    double EndPointY;		// The other end of the interval (Y dir)
    double RelErr;		// Relative error
    double EpsMachine;		// Machine accuracy
 };

   /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (LinearSystems_TestGroup) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_AdaptiveGaussianQuadrature2D> AdaptiveGaussianQuadrature2D_TestGroup;
  typedef AdaptiveGaussianQuadrature2D_TestGroup::object
                                             AdaptiveGaussianQuadrature2D_TestObject;
}

/**********************************************
 *                THE TESTS                   *
 *********************************************/

namespace tut
{

  /* Test 1:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature2D_TestObject::test<1>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Default2D;

    // Negative interval
    StartPointX = -125.00;
    EndPointX = -89.0;
    StartPointY = -125.00;
    EndPointY = -89.0;

    TheoreticResult = Test_Default2D_Integral(StartPointX,EndPointX,StartPointY,EndPointY);
    NumericResultA = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,digits);
    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,EndPointY,StartPointY,digits);
    NumericResult = AdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,
					       digits,NumericResult);
    ensure("Different Methods", fabs(NumericResult-NumericResultA) <= tol);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits Y", RelErr <= tol);

    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,EndPointX,StartPointX,StartPointY,EndPointY,digits);
    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits X", RelErr <= tol);

  }

  /* Test 2:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature2D_TestObject::test<2>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Default2D;

    // Negative interval and Zero
    StartPointX = -30.00;
    EndPointX = 0.0;
    StartPointY = -125.00;
    EndPointY = -89.0;

    TheoreticResult = Test_Default2D_Integral(StartPointX,EndPointX,StartPointY,EndPointY);
    NumericResultA = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,digits);
    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,EndPointY,StartPointY,digits);
    NumericResult = AdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,
					       digits,NumericResult);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResult)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResult);

    // Check results
    ensure("Check integration", RelErr <= tol);


    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits Y", RelErr <= tol);

    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,EndPointX,StartPointX,StartPointY,EndPointY,digits);
    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits X", RelErr <= tol);
  }

  /* Test 3:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature2D_TestObject::test<3>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Default2D;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 1.0E-14;
    StartPointY = -25.00;
    EndPointY = 0.0;

    TheoreticResult = Test_Default2D_Integral(StartPointX,EndPointX,StartPointY,EndPointY);
    NumericResultA = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,digits);
    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,EndPointY,StartPointY,digits);
    NumericResult = AdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,
					       digits,NumericResult);
    ensure("Different Methods", fabs(NumericResult-NumericResultA) <= tol);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits Y", RelErr <= tol);

    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,EndPointX,StartPointX,StartPointY,EndPointY,digits);
    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits X", RelErr <= tol);
  }

  /* Test 4:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature2D_TestObject::test<4>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Default2D;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 1.0E-14;
    // Interval close to zero
    StartPointY = -1.15E-14;
    EndPointY = 0.2345;

    TheoreticResult = Test_Default2D_Integral(StartPointX,EndPointX,StartPointY,EndPointY);
    NumericResultA = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,digits);
    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,EndPointY,StartPointY,digits);
    NumericResult = AdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,
					       digits,NumericResult);
    ensure("Different Methods", fabs(NumericResult-NumericResultA) <= tol);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits Y", RelErr <= tol);

    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,EndPointX,StartPointX,StartPointY,EndPointY,digits);
    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits X", RelErr <= tol);
  }

  /* Test 5:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature2D_TestObject::test<5>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Default2D;

    // Positive interval
    StartPointX = 89.00;
    EndPointX = 125.0;
    StartPointY = 1250.00;
    EndPointY = 8900.0;

    TheoreticResult = Test_Default2D_Integral(StartPointX,EndPointX,StartPointY,EndPointY);
    NumericResultA = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,digits);
    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,EndPointY,StartPointY,digits);
    NumericResult = AdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,
					       digits,NumericResult);
    ensure("Different Methods", fabs(NumericResult-NumericResultA) <= tol);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits Y", RelErr <= tol);

    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,EndPointX,StartPointX,StartPointY,EndPointY,digits);
    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits X", RelErr <= tol);
  }

  /* Test 6:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature2D_TestObject::test<6>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = 1.0E-10;
    func = Test_Example2;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 1.0E-14;
    // Interval close to zero
    StartPointY = -1.15E-14;
    EndPointY = 0.24;

    TheoreticResult = Test_Example2_Integral(StartPointX,EndPointX,StartPointY,EndPointY);
    NumericResultA = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,digits);
    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,EndPointY,StartPointY,digits);
    NumericResult = AdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,
					       digits,NumericResult);
    ensure("Different Methods", fabs(NumericResult-NumericResultA) <= tol);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits Y", RelErr <= tol);

    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,EndPointX,StartPointX,StartPointY,EndPointY,digits);
    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits X", RelErr <= tol);
  }

  /* Test 7:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature2D_TestObject::test<7>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10-4;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example2;

    // Interval close to zero
    StartPointX = -16;
    EndPointX = 2.0;
    // Interval close to zero
    StartPointY = -1.3;
    EndPointY = 20;

    TheoreticResult = Test_Example2_Integral(StartPointX,EndPointX,StartPointY,EndPointY);
    NumericResultA = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,digits);
    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,EndPointY,StartPointY,digits);
    NumericResult = AdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,
					       digits,NumericResult);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResult)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResult);

    // Check results
    ensure("Check integration", RelErr <= tol);


    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits Y", RelErr <= tol);

    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,EndPointX,StartPointX,StartPointY,EndPointY,digits);
    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits X", RelErr <= tol);
  }

  /* Test 8:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature2D_TestObject::test<8>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example3;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 1.0E-14;
    // Interval close to zero
    StartPointY = -1.15E-14;
    EndPointY = 0.24;

    TheoreticResult = Test_Example3_Integral(StartPointX,EndPointX,StartPointY,EndPointY);
    NumericResultA = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,digits);
    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,EndPointY,StartPointY,digits);
    NumericResult = AdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,
					       digits,NumericResult);
    ensure("Different Methods", fabs(NumericResult-NumericResultA) <= tol);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits Y", RelErr <= tol);

    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,EndPointX,StartPointX,StartPointY,EndPointY,digits);
    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits X", RelErr <= tol);
  }

  /* Test 9:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature2D_TestObject::test<9>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10-8; // this function has singularities --> lower accuracy
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example6;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 2.0;
    // Interval close to zero
    StartPointY = -1.15E-14;
    EndPointY = 0.24;

    TheoreticResult = Test_Example6_Integral(StartPointX,EndPointX,StartPointY,EndPointY);
    NumericResultA = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,digits);
    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,EndPointY,StartPointY,digits);
    NumericResult = AdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,
					       digits,NumericResult);
    ensure("Different Methods", fabs(NumericResult-NumericResultA) <= tol);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits Y", RelErr <= tol);

    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,EndPointX,StartPointX,StartPointY,EndPointY,digits);
    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits X", RelErr <= tol);
  }

 /* Test 10:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature2D_TestObject::test<10>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10-7;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example7;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 2.0;
    // Interval close to zero
    //    StartPointY = -1.15E-14;
    StartPointY = -1.0;
    EndPointY = 0.24;

    TheoreticResult = Test_Example7_Integral(StartPointX,EndPointX,StartPointY,EndPointY);
    NumericResultA = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,digits);
    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,StartPointX,EndPointX,EndPointY,StartPointY,digits);
    NumericResult = AdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,
					       digits,NumericResult);
    ensure("Different Methods", fabs(NumericResult-NumericResultA) <= tol);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits Y", RelErr <= tol);

    NumericResultB = quad2dAdaptiveGaussianQuadrature(func,EndPointX,StartPointX,StartPointY,EndPointY,digits);
    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits X", RelErr <= tol);
  }
};

namespace tut
{
  // Data used for testing
  struct Data_AdaptiveGaussianQuadrature3D{
    int digits;			// Number of exact imposed digits --> gives the precision
    double tol;	                // Numerical tolerance for comparing the solutions
    FunctionType3D func;        // Use 3D function for integration
    double TheoreticResult;	// Result obtained with the exact integration
    double NumericResultA, NumericResultB, NumericResult;
    double StartPointX;		// One end of the interval (X dir)
    double EndPointX;		// The other end of the interval (X dir)
    double StartPointY;		// One end of the interval (Y dir)
    double EndPointY;		// The other end of the interval (Y dir)
    double StartPointZ;		// One end of the interval (Z dir)
    double EndPointZ;		// The other end of the interval (Z dir)
    double RelErr;		// Relative error
    double EpsMachine;		// Machine accuracy
 };

   /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   */
  typedef test_group<Data_AdaptiveGaussianQuadrature3D> AdaptiveGaussianQuadrature3D_TestGroup;
  typedef AdaptiveGaussianQuadrature3D_TestGroup::object
                                             AdaptiveGaussianQuadrature3D_TestObject;
}

/**********************************************
 *                THE TESTS                   *
 *********************************************/

namespace tut
{


 /* Test 1:*/
  template<>
  template<>
  void AdaptiveGaussianQuadrature3D_TestObject::test<1>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10-7;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Default3D;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 2.0;
    // Interval close to zero
    StartPointY = -1.0;
    EndPointY = 0.24;

    StartPointZ = -1.34;
    EndPointZ = -0.04;

    // Compare with the theoretical value
    TheoreticResult = Test_Default3D_Integral(StartPointX,EndPointX,StartPointY,EndPointY,
					      StartPointZ,EndPointZ);
    NumericResultA = AdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,
						StartPointZ,EndPointZ, digits,NumericResult);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    // Compare with the value obtained by swhitching the X limits
    NumericResultB = AdaptiveGaussianQuadrature(func,EndPointX,StartPointX,StartPointY,EndPointY,
						StartPointZ,EndPointZ, digits,NumericResult);

    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits X", RelErr <= tol);

    // Compare with the value obtained by swhitching the Y limits
    NumericResultB = AdaptiveGaussianQuadrature(func,StartPointX,EndPointX,EndPointY,StartPointY,
						StartPointZ,EndPointZ, digits,NumericResult);

    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits Y", RelErr <= tol);

    // Compare with the value obtained by swhitching the Z limits
    NumericResultB = AdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,
						EndPointZ,StartPointZ,digits,NumericResult);

    RelErr = fabs((NumericResultA + NumericResultB));
    ensure("Check independency of limits Z", RelErr <= tol);
  }

}

/************************************************************************************
 Tests for  QuadrilateralQuadrature function
 ***********************************************************************************/

namespace tut
{
  // Data used for testing
  struct Data_QuadrilateralQuadrature{
    int digits;			// Number of exact imposed digits --> gives the precision
    double tol;	                // Numerical tolerance for comparing the solutions
    FunctionType2D func;        // Use 3D function for integration
    double TheoreticResult;	// Result obtained with the exact integration
    double NumericResultA, NumericResult;
    double StartPointX;		// One end of the interval (X dir)
    double EndPointX;		// The other end of the interval (X dir)
    double StartPointY;		// One end of the interval (Y dir)
    double EndPointY;		// The other end of the interval (Y dir)
    double RelErr;		// Relative error
    double EpsMachine;		// Machine accuracy
 };

   /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   */
  typedef test_group<Data_QuadrilateralQuadrature> QuadrilateralQuadrature_TestGroup;
  typedef QuadrilateralQuadrature_TestGroup::object
                                             QuadrilateralQuadrature_TestObject;
}

/**********************************************
 *                THE TESTS                   *
 *********************************************/

namespace tut
{

 /* Test 1: QuadrilateralQuadrature() over a rectangle*/
  template<>
  template<>
  void QuadrilateralQuadrature_TestObject::test<1>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10-8;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example6;

    // Interval close to zero
    StartPointX = -0.00000016;
    EndPointX = 2.0;
    // Interval close to zero
    StartPointY = -1.15E-14;
    EndPointY = 0.24;

    // Generate the nodes
    Vector2D Node1, Node2, Node3, Node4;

    Node1.x = StartPointX; Node1.y = StartPointY;
    Node2.x = StartPointX; Node2.y = EndPointY;
    Node3.x = EndPointX; Node3.y = EndPointY;
    Node4.x = EndPointX; Node4.y = StartPointY;


    // Compare with the theoretical value
    TheoreticResult = Test_Example6_Integral(StartPointX,EndPointX,StartPointY,EndPointY);

    NumericResultA = AdaptiveGaussianQuadrature(func,StartPointX,EndPointX,StartPointY,EndPointY,
						digits,NumericResultA);

    NumericResult = QuadrilateralQuadrature(func,Node1,Node2,Node3,Node4,digits,NumericResult);
    // Obs: The order of the nodes must be clockwise in order to obtain the correct value of the integral 

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResultA)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResultA);

    // Check results
    ensure("Check integration", RelErr <= tol);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResult)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResult);

    // Check results
    ensure("Check integration", RelErr <= tol);

  }


 /* Test 2: The value returned by BilinearTransformFunctionInPlan() on a unit square*/
  template<>
  template<>
  void QuadrilateralQuadrature_TestObject::test<2>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10-8;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example6;

    StartPointX = 0.0;
    EndPointX = 1.0;
    // Interval close to zero
    StartPointY = 0.0;
    EndPointY = 1.0;

    // Generate the nodes
    Vector2D Node1, Node2, Node3, Node4;

    Node1.x = StartPointX; Node1.y = StartPointY;
    Node2.x = StartPointX; Node2.y = EndPointY;
    Node3.x = EndPointX; Node3.y = EndPointY;
    Node4.x = EndPointX; Node4.y = StartPointY;

    double x,y;
    x = 0.5*(StartPointX+EndPointX);
    y = 0.5*(StartPointY+EndPointY);

    // Compare with the theoretical value
    TheoreticResult = func(x,y);

    BilinearTransformFunctionInPlan<FunctionType2D,Vector2D,double> MyFunction(func,Node1,Node2,Node3,Node4);

    NumericResult = MyFunction(x,y);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResult)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResult);

    // Check results
    ensure("Check Function Value", RelErr <= tol);
  }

 /* Test 3: QuadrilateralQuadrature() over a trapezoid*/
  template<>
  template<>
  void QuadrilateralQuadrature_TestObject::test<3>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example10;

    // Generate the nodes
    Vector2D Node1, Node2, Node3, Node4;

    Node1.x =  0.0; Node1.y = -2.0;
    Node2.x =  0.0; Node2.y = -1.0;
    Node3.x =  1.0; Node3.y = 0.0;
    Node4.x =  2.0; Node4.y = 0.0;

    // Compare with the theoretical value
    TheoreticResult = 3.0/4.0*(exp(1.0) - exp(-1.0));

    NumericResult = QuadrilateralQuadrature(func,Node1,Node2,Node3,Node4,digits,NumericResult);
    // Obs: The order of the nodes must be clockwise in order to obtain the correct value of the integral 

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResult)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResult);

    // Check results
    ensure("Check integration", RelErr <= tol);

  }

 /* Test 4: QuadrilateralQuadrature() over a triangle*/
  template<>
  template<>
  void QuadrilateralQuadrature_TestObject::test<4>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example11;

    // Generate the nodes
    Vector2D Node1, Node2, Node3, Node4;

    Node1.x =  0.0; Node1.y = 0.0;
    Node2.x =  0.0; Node2.y = 20.0;
    Node3.x =  13.0; Node3.y = 7.0;
    Node4.x =  20.0; Node4.y = 0.0;

    // Compare with the theoretical value
    TheoreticResult = 1 + exp(-4.0) - 2*exp(-2.0);

    NumericResult = QuadrilateralQuadrature(func,Node1,Node2,Node3,Node4,digits,NumericResult);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResult)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResult);

    // Check results
    ensure("Check integration", RelErr <= tol);

  }


 /* Test 5: QuadrilateralQuadrature() over a quadrilateral cell
    -> Computation of the cell area*/
  template<>
  template<>
  void QuadrilateralQuadrature_TestObject::test<5>()
  {
    // Initialize variables
    digits = numeric_limits<double>::digits10;
    tol = 0.5*pow(10.0,1.0-digits);
    EpsMachine = numeric_limits<long double>::epsilon( );
    func = Test_Example12;

    // Generate the nodes
    Vector2D Node1, Node2, Node3, Node4;
    Vector2D SE, SW, NW, NE;

    Node1.x = -1.2344-1000; Node1.y = -2.7635+1000;
    Node2.x = -5.121-1000; Node2.y = 3.1234+1000;
    Node3.x =  2.56-1000; Node3.y = 1.2334+1000;
    Node4.x =  10.2323-1000; Node4.y = 2.89762+1010;

    SW = Node1; NW = Node2; NE = Node3; SE = Node4;

    // Compare with the theoretical value
    TheoreticResult = (0.5*(((SE-SW)^(NW-SW))+((NE-NW)^(NE-SE))));;

    NumericResult = QuadrilateralQuadrature(func,Node1,Node2,Node3,Node4,digits,NumericResult);

    if ( fabs(TheoreticResult) > EpsMachine)
      RelErr = fabs((TheoreticResult - NumericResult)/TheoreticResult);
    else
      RelErr = fabs(TheoreticResult - NumericResult);

    // Check results
    ensure("Check integration", RelErr <= tol);

  }

}

namespace tut
{
  AdaptiveGaussianQuadrature_TestGroup AdaptiveGaussianQuadrature1D_Test("AdaptiveGaussianQuadrature1D_Test");
  AdaptiveGaussianQuadrature2D_TestGroup AdaptiveGaussianQuadrature2D_Test("AdaptiveGaussianQuadrature2D_Test");
  AdaptiveGaussianQuadrature3D_TestGroup AdaptiveGaussianQuadrature3D_Test("AdaptiveGaussianQuadrature3D_Test");
  QuadrilateralQuadrature_TestGroup QuadrilateralQuadrature_Test("QuadrilateralQuadrature_Test");
}
