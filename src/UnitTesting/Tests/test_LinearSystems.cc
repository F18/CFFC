// Test file for class Cell1D_Uniform
#include <tut.h>
#include "Math/LinearSystems.h"
#include <cmath>
#include "include/require.h"

namespace tut

{
  // Data used for testing
  struct Data_LinearSystems{
    int krank;			// The rank of the matrix
    double RezNorm;		// The residual norm of the Least Square problem
    double RezNormMatlab;	// The residual norm based on Matlab solution
    double tau;           
    double tol;	                // Numerical tolerance for comparing the solutions
  };

   /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (LinearSystems_TestGroup) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_LinearSystems> LinearSystems_TestGroup;
  typedef LinearSystems_TestGroup::object LinearSystems_TestObject;
}

/**********************************************
 *                THE TESTS                   *
 *********************************************/

// ************* Matrix Model ***************
//     A(0,0)=  ; A(0,1)= ; A(0,2)= ;
//     A(1,0)=  ; A(1,1)= ; A(1,2)= ;
//     A(2,0)=  ; A(2,1)= ; A(2,2)= ;
//     A(3,0)=  ; A(3,1)= ; A(3,2)= ;


namespace tut
{

  /* Test 1:*/
  template<>
  template<>
  void LinearSystems_TestObject::test<1>()
  {
    DenseMatrix A(4,3); 	// Test matrix
    ColumnVector B(4);		// Test free term
    ColumnVector X(3);		// Solution vector
    ColumnVector XMatlab(3);	// Matlab solution vector
    ColumnVector SolutionDiff(3);	// Difference between the two solutions
    // Initialize variables
    tau = 10.0;
    tol = 1.0e-14;

    // Initialize test matrix & free term
    A(0,0)= -2.093656e-01; A(0,1)= 2.114593e-01; A(0,2)= -1.512815e-01;
    A(1,0)= -9.900990e-01; A(1,1)= 5.000000e-01; A(1,2)= -2.104167e-01;
    A(2,0)=  9.900990e-01; A(2,1)= 5.000000e-01; A(2,2)=  2.104167e-01;
    A(3,0)=  4.950495e-01; A(3,1)= 5.000000e-01; A(3,2)= 3.577083e-01;

    B(0) = 4.364948e-01; B(1) = 4.353389e-16; B(2) = 0.000000e+00; B(3) = 0.000000e+00;

    // Initialize Matlab solution
    XMatlab(0) = 0.20223247422880; XMatlab(1) = 0.24027609285144;
    XMatlab(2) =-0.95158872955264;
    RezNormMatlab = 0.35212869025901;

    Solve_LS_Householder(A,B,X,krank,RezNorm);   

    // check condition
    for (int i=0; i<A.size(1); i++){
      SolutionDiff(i) = fabs(X(i)-XMatlab(i));
      ensure ("System 4x3",SolutionDiff(i)<=tol);
    }
    ensure ("System 4x3",fabs(RezNormMatlab-RezNorm)<=tol);
  }

//   /*Test 2:*/
//   template<>
//   template<>
//   void LinearSystems_TestObject::test<2>()
//   {
    
//     DenseMatrix A(4,5); 	// Test matrix
//     ColumnVector B(4);		// Test free term
//     ColumnVector X(5);		// Solution vector
//     ColumnVector XMatlab(5);	// Matlab solution vector
//     ColumnVector SolutionDiff(5);	// Difference between the two solutions
//     // Initialize variables
//     tau = 10.0;
//     tol = 1.0e-14;
    
//     // Initialize test matrix & free term
//     A(0,0)= 1.000000e-02; A(0,1)= 1.000000e-01; A(0,2)= 1.000000e+00;
//     A(0,3)= 1.000000e+01; A(0,4)= 1.000000e+02;
//     A(1,0)= 2.500000e-02; A(1,1)= 2.500000e-01; A(1,2)= 2.500000e+00;
//     A(1,3)= 2.500000e+01; A(1,4)= 2.500000e+02;
//     A(2,0)= 3.500000e+01; A(2,1)= 3.500000e+00; A(2,2)= 3.500000e-01;
//     A(2,3)= 3.500000e-02; A(2,4)= 3.500000e-03;
//     A(3,0)= 4.350000e+02; A(3,1)= 4.350000e-01; A(3,2)= 4.350000e-02;
//     A(3,3)= 4.350000e-03; A(3,4)= 4.350000e-06;
    
//     B(0) = 1.0/3.0*0.000000e+00; B(1) = 1.0/3.0*2.400000e+01;
//     B(2) = 1.0/3.0*3.200000e+01; B(3) = 1.0/3.0*-4.300000e-16;
    
//     // Initialize Matlab solution
//     XMatlab(0) = -0.00307838478822; XMatlab(1) = 3.04757383890261;
//     XMatlab(2) = 0.30490798411287; XMatlab(3)=0.03199680063742;
//     XMatlab(4) = 0.01829018099126;
//     RezNormMatlab = 2.97112541083283;
    
//     Solve_LS_Householder(A,B,X,krank,RezNorm);   
    
//     // check condition
//     for (int i=0; i<A.size(1); i++){
//       SolutionDiff(i) = fabs(X(i)-XMatlab(i));
//       ensure ("Solution",SolutionDiff(i)<=tol);
//     }
//     ensure ("Residual",fabs(RezNormMatlab-RezNorm)<=tol);
//   }
  
  /* Test 3 */
  template<>
  template<>
  void LinearSystems_TestObject::test<3>()
  {
    DenseMatrix A(6,5); 	// Test matrix
    ColumnVector B(6);		// Test free term
    ColumnVector X(5);		// Solution vector
    ColumnVector XMatlab(5);	// Matlab solution vector
    ColumnVector SolutionDiff(5);	// Difference between the two solutions
    // Initialize variables
    tau = 10.0;
    tol = 1.0e-14;
  
    // Initialize test matrix & free term
    A(0,0)= -3.300330e-01; A(0,1)=  5.000000e-01; A(0,2)= -5.190278e-01;
    A(0,3)=  4.037896e-01; A(0,4)= -2.534612e-01;
    A(1,0)= -4.950495e-01; A(1,1)=  5.000000e-01; A(1,2)= -3.577083e-01;
    A(1,3)=  1.912687e-01; A(1,4)= -8.326478e-02;
    A(2,0)= -9.900990e-01; A(2,1)=  5.000000e-01; A(2,2)= -2.104167e-01;
    A(2,3)=  6.375625e-02; A(2,4)= -1.627732e-02;
    A(3,0)=  9.900990e-01; A(3,1)=  5.000000e-01; A(3,2)=  2.104167e-01;
    A(3,3)=  6.375625e-02; A(3,4)=  1.627732e-02;
    A(4,0)=  4.950495e-01; A(4,1)=  5.000000e-01; A(4,2)=  3.577083e-01;
    A(4,3)=  1.912687e-01; A(4,4)=  8.326478e-02;
    A(5,0)=  2.501813e-01; A(5,1)=  3.790247e-01; A(5,2)=  3.934487e-01;
    A(5,3)=  3.060924e-01; A(5,4)= 1.921361e-01;

    B(0) = 1.934840e-16; B(1) = 4.353389e-16;
    B(2) = 0.000000e+00; B(3) = 1.741356e-15;
    B(4) = 4.353389e-16; B(5) = -3.995151e-01;

    // Initialize Matlab solution
    XMatlab(0) = -0.10416722683563; XMatlab(1) = 0.17943718969498;
    XMatlab(2) = 0.66236548323466; XMatlab(3)= -0.79538109904828;
    XMatlab(4) = -2.22621976242727;
    RezNormMatlab = 0.11106331948011;

    Solve_LS_Householder(A,B,X,krank,RezNorm);   

    // check condition
    for (int i=0; i<A.size(1); i++){
      SolutionDiff(i) = fabs(X(i)-XMatlab(i));
      ensure ("System 6x5",SolutionDiff(i)<=tol);
    }
    ensure ("System 6x5",fabs(RezNormMatlab-RezNorm)<=tol);
  }

  /* Test 4 */
  template<>
  template<>
  void LinearSystems_TestObject::test<4>()
  {
    DenseMatrix A(4,3); 	// Test matrix
    ColumnVector B(4);		// Test free term
    ColumnVector X(3);		// Solution vector
    ColumnVector XMatlab(3);	// Matlab solution vector
    ColumnVector SolutionDiff(3);	// Difference between the two solutions
    // Initialize variables
    tau = 10.0;
    tol = 1.0e-14;

    // Initialize test matrix & free term
    A(0,0)=  9.480438e-02; A(0,1)= 9.575242e-02; A(0,2)= 6.850288e-02;
    A(1,0)= -0.000000e+00; A(1,1)= 0.000000e+00; A(1,2)=-0.000000e+00;
    A(2,0)=  0.000000e+00; A(2,1)= 0.000000e+00; A(2,2)= 0.000000e+00;
    A(3,0)= -0.000000e+00; A(3,1)= 0.000000e+00; A(3,2)= -0.000000e+00;

    B(0) = 5.121065e-02; B(1) = -0.000000e+00;
    B(2) = 0.000000e+00; B(3) =  0.000000e+00;

    // Initialize Matlab solution
    XMatlab(0) = 0.21248129970697; XMatlab(1) = 0.21460610418725;
    XMatlab(2) = 0.15353279011023;
    RezNormMatlab = 1.387778780781446e-17;

    Solve_LS_Householder(A,B,X,krank,RezNorm);   

    // check condition
    for (int i=0; i<A.size(1); i++){
      SolutionDiff(i) = fabs(X(i)-XMatlab(i));
      ensure ("System 4x3",SolutionDiff(i)<=tol);
    }
    ensure ("System 4x3",fabs(RezNormMatlab-RezNorm)<=tol);
  }

  /* Test 5 */
  template<>
  template<>
  void LinearSystems_TestObject::test<5>()
  {
    DenseMatrix A(4,3); 	// Test matrix
    ColumnVector B(4);		// Test free term
    ColumnVector X(3);		// Solution vector
    ColumnVector XMatlab(3);	// Matlab solution vector
    ColumnVector SolutionDiff(3);	// Difference between the two solutions
    // Initialize variables
    tau = 10.0;
    tol = 1.0e-14;

    // Initialize test matrix & free term
    A(0,0)= -4.942834e-01; A(0,1)= 4.992262e-01; A(0,2)= -3.571548e-01;
    A(1,0)= -9.900990e-01; A(1,1)= 5.000000e-01; A(1,2)= -2.104167e-01;
    A(2,0)=  9.900990e-01; A(2,1)= 5.000000e-01; A(2,2)=  2.104167e-01;
    A(3,0)=  4.950495e-01; A(3,1)= 5.000000e-01; A(3,2)=  3.577083e-01;

    B(0) = -2.669979e-01; B(1) = -7.080900e-03;
    B(2) = -2.464985e-04; B(3) = -6.289122e-05;

    // Initialize Matlab solution
    XMatlab(0) = -0.10753505654966; XMatlab(1) = -0.13730018989479;
    XMatlab(2) = 0.52223777268398;
    RezNormMatlab = 0.13002318505849;

    Solve_LS_Householder(A,B,X,krank,RezNorm);   

    // check condition
    for (int i=0; i<A.size(1); i++){
      SolutionDiff(i) = fabs(X(i)-XMatlab(i));
      ensure ("System 4x3",SolutionDiff(i)<=tol);
    }
    ensure ("System 4x3",fabs(RezNormMatlab-RezNorm)<=tol);
  }

  /* Test 6 */
  template<>
  template<>
  void LinearSystems_TestObject::test<6>()
  {
    DenseMatrix A(6,4); 	// Test matrix
    ColumnVector B(6);		// Test free term
    ColumnVector X(4);		// Solution vector
    ColumnVector XMatlab(4);	// Matlab solution vector
    ColumnVector SolutionDiff(4);	// Difference between the two solutions
    // Initialize variables
    tau = 10.0;
    tol = 1.0e-14;

    // Initialize test matrix & free term
    A(0,0)= -0.000000e+00; A(0,1)= 0.000000e+00; A(0,2)= -0.000000e+00;
    A(0,3)=  0.000000e+00;
    A(1,0)= -0.000000e+00; A(1,1)= 0.000000e+00; A(1,2)= -0.000000e+00;
    A(1,3)=  0.000000e+00;
    A(2,0)= -0.000000e+00; A(2,1)= 0.000000e+00; A(2,2)= -0.000000e+00;
    A(2,3)=  0.000000e+00;
    A(3,0)=  0.000000e+00; A(3,1)= 0.000000e+00; A(3,2)=  0.000000e+00;
    A(3,3)=  0.000000e+00;
    A(4,0)=  0.000000e+00; A(4,1)= 0.000000e+00; A(4,2)=  0.000000e+00;
    A(4,3)=  0.000000e+00;
    A(5,0)=  4.212746e-02; A(5,1)= 6.382311e-02; A(5,2)=  6.625193e-02;
    A(5,3)=  5.154221e-02;

    B(0) =  0.000000e+00; B(1) =  0.000000e+00;
    B(2) = -0.000000e+00; B(3) =  0.000000e+00;
    B(4) =  0.000000e+00; B(5) =  1.516721e-02;

    // Initialize Matlab solution
    XMatlab(0) = 0.04955440911134; XMatlab(1) = 0.07507493933168;
    XMatlab(2) = 0.07793195325889; XMatlab(3) = 0.06062895225211;
    RezNormMatlab = 3.469446951953614e-18;
    
    Solve_LS_Householder(A,B,X,krank,RezNorm);   

    // check condition
    for (int i=0; i<A.size(1); i++){
      SolutionDiff(i) = fabs(X(i)-XMatlab(i));
      ensure ("System 6x4",SolutionDiff(i)<=tol);
    }
    ensure ("System 6x4",fabs(RezNormMatlab-RezNorm)<=tol);
  }

  /* Test 7 */
  template<>
  template<>
  void LinearSystems_TestObject::test<7>()
  {
    DenseMatrix A(2,1); 	// Test matrix
    ColumnVector B(2);		// Test free term
    ColumnVector X(1);		// Solution vector
    ColumnVector XMatlab(1);	// Matlab solution vector
    ColumnVector SolutionDiff(1);	// Difference between the two solutions
    // Initialize variables
    tau = 10.0;
    tol = 1.0e-14;

    // Initialize test matrix & free term
    A(0,0)= -3.000000e+00;
    A(1,0)= 3.000000e+00;

    B(0) =  0.000000e+00; B(1) = -4.285714e-01;

    // Initialize Matlab solution
    XMatlab(0) = -0.07142856666667;
    RezNormMatlab = 0.30304574316261;
    
    Solve_LS_Householder(A,B,X,krank,RezNorm);   

    // check condition
    for (int i=0; i<A.size(1); i++){
      SolutionDiff(i) = fabs(X(i)-XMatlab(i));
      ensure ("System 2x1",SolutionDiff(i)<=tol);
    }
    ensure ("System 2x1",fabs(RezNormMatlab-RezNorm)<=tol);
  }

  /* Test 8 */
  template<>
  template<>
  void LinearSystems_TestObject::test<8>()
  {
    DenseMatrix A(2,1); 	// Test matrix
    ColumnVector B(2);		// Test free term
    ColumnVector X(1);		// Solution vector
    ColumnVector XMatlab(1);	// Matlab solution vector
    ColumnVector SolutionDiff(1);	// Difference between the two solutions
    // Initialize variables
    tau = 10.0;
    tol = 1.0e-14;

    // Initialize test matrix & free term
    A(0,0)= -3.000000e+00;
    A(1,0)= 3.000000e+00;
    B(0) =  0.000000e+00; B(1) = -4.714286e+00;

    // Initialize Matlab solution

    XMatlab(0) = -0.78571433333333;
    RezNormMatlab = 3.33350359905281;
    
    Solve_LS_Householder(A,B,X,krank,RezNorm);   

    // check condition
    for (int i=0; i<A.size(1); i++){
      SolutionDiff(i) = fabs(X(i)-XMatlab(i));
      ensure ("System 2x1",SolutionDiff(i)<=tol);
    }
    ensure ("System 2x1",fabs(RezNormMatlab-RezNorm)<=tol);
  }

  /* Test 9 */
  template<>
  template<>
  void LinearSystems_TestObject::test<9>()
  {
    DenseMatrix A(4,2); 	// Test matrix
    ColumnVector B(4);		// Test free term
    ColumnVector X(2);		// Solution vector
    ColumnVector XMatlab(2);	// Matlab solution vector
    ColumnVector SolutionDiff(2);	// Difference between the two solutions
    // Initialize variables
    tau = 10.0;
    tol = 1.0e-14;

    // Initialize test matrix & free term
    A(0,0)= -1.500000e+00; A(0,1)= 5.000000e-01;
    A(1,0)= -3.000000e+00; A(1,1)= 5.000000e-01;
    A(2,0)=  3.000000e+00; A(2,1)= 5.000000e-01;
    A(3,0)=  1.500000e+00; A(3,1)= 5.000000e-01;

    B(0) = 0.000000e+00; B(1) = 0.000000e+00;
    B(2) =-8.914286e+00; B(3) =-2.228571e+00;

    // Initialize Matlab solution
    
    XMatlab(0) = -1.33714286666667; XMatlab(1) = -5.57142850000000;
    RezNormMatlab = 3.62785191662778;
    
    Solve_LS_Householder(A,B,X,krank,RezNorm);   

    // check condition
    for (int i=0; i<A.size(1); i++){
      SolutionDiff(i) = fabs(X(i)-XMatlab(i));
      ensure ("System 4x2",SolutionDiff(i)<=tol);
    }
    ensure ("System 4x2",fabs(RezNormMatlab-RezNorm)<=tol);
  }

  /* Test 10 */
  template<>
  template<>
  void LinearSystems_TestObject::test<10>()
  {
    DenseMatrix A(4,2); 	// Test matrix
    ColumnVector B(4);		// Test free term
    ColumnVector X(2);		// Solution vector
    ColumnVector XMatlab(2);	// Matlab solution vector
    ColumnVector SolutionDiff(2);	// Difference between the two solutions
    // Initialize variables
    tau = 10.0;
    tol = 1.0e-14;

    // Initialize test matrix & free term
    A(0,0)= -1.500000e+00; A(0,1)= 5.000000e-01;
    A(1,0)= -3.000000e+00; A(1,1)= 5.000000e-01;
    A(2,0)=  3.000000e+00; A(2,1)= 5.000000e-01;
    A(3,0)=  1.500000e+00; A(3,1)= 5.000000e-01;

    B(0) = 9.516197e-16; B(1) = 0.000000e+00;
    B(2) = 0.000000e+00; B(3) =-1.071429e-01;

    // Initialize Matlab solution
    XMatlab(0) = -0.00714286000000; XMatlab(1) = -0.05357145000000;
    RezNormMatlab = 0.08638136757002;
    
    Solve_LS_Householder(A,B,X,krank,RezNorm);   

    // check condition
    for (int i=0; i<A.size(1); i++){
      SolutionDiff(i) = fabs(X(i)-XMatlab(i));
      ensure ("System 4x2",SolutionDiff(i)<=tol);
    }
    ensure ("System 4x2",fabs(RezNormMatlab-RezNorm)<=tol);
  }

  /* Test 11 */
  template<>
  template<>
  void LinearSystems_TestObject::test<11>()
  {
    DenseMatrix A(6,4); 	// Test matrix
    ColumnVector B(6);		// Test free term
    ColumnVector X(4);		// Solution vector
    ColumnVector XMatlab(4);	// Matlab solution vector
    ColumnVector SolutionDiff(4);	// Difference between the two solutions
    // Initialize variables
    tau = 10.0;
    tol = 1.0e-14;

    // Initialize test matrix & free term
    A(0,0)= -1.000000e-00; A(0,1)= 5.000000e-01; A(0,2)=-1.712963e-01; A(0,3)= 4.398148e-02;
    A(1,0)= -1.500000e+00; A(1,1)= 5.000000e-01; A(1,2)=-1.180556e-01; A(1,3)= 2.083333e-02;
    A(2,0)= -3.000000e+00; A(2,1)= 5.000000e-01; A(2,2)=-6.944444e-02; A(2,3)= 6.944444e-03;
    A(3,0)=  3.000000e+00; A(3,1)= 5.000000e-01; A(3,2)= 6.944444e-02; A(3,3)= 6.944444e-03;
    A(4,0)=  1.500000e+00; A(4,1)= 5.000000e-01; A(4,2)= 1.180556e-01; A(4,3)= 2.083333e-02;
    A(5,0)=  1.000000e-00; A(5,1)= 5.000000e-01; A(5,2)= 1.712963e-01; A(5,3)= 4.398148e-02;


    B(0) = 4.229421e-16; B(1) = 0.000000e+00;
    B(2) = 0.000000e+00; B(3) =-4.285714e-01;
    B(4) =-1.071429e-01; B(5) =-4.761905e-02;

    // Initialize Matlab solution
    XMatlab(0) = -0.07825242962236; XMatlab(1) = -0.42164719943979;
    XMatlab(2) =  0.38003146194418; XMatlab(3) =  4.74927072346217; 
    RezNormMatlab = 0.10689757419026;
    
    Solve_LS_Householder(A,B,X,krank,RezNorm);   

    // check condition
    for (int i=0; i<A.size(1); i++){
      SolutionDiff(i) = fabs(X(i)-XMatlab(i));
      ensure ("System 6x4",SolutionDiff(i)<=tol);
    }
    ensure ("System 6x4",fabs(RezNormMatlab-RezNorm)<=tol);
  }

  /* Test 12 */
  template<>
  template<>
  void LinearSystems_TestObject::test<12>()
  {
    DenseMatrix A(8,6); 	// Test matrix
    ColumnVector B(8);		// Test free term
    ColumnVector X(6);		// Solution vector
    ColumnVector XMatlab(6);	// Matlab solution vector
    ColumnVector SolutionDiff(6);	// Difference between the two solutions
    // Initialize variables
    tau = 1.0;
    tol = 1.0e-10;   		// The condition number of this matrix is k = 5.5978e+04

    // Initialize test matrix & free term
    A(0,0)= -1.250000e+00; A(0,1)= 5.000000e-01; A(0,2)= -1.354167e-01; A(0,3)=  2.750000e-02;
    A(0,4)= -4.489931e-03; A(0,5)= 6.137500e-04 ;
    A(1,0)= -1.666667e+00; A(1,1)= 5.000000e-01; A(1,2)= -1.027778e-01; A(1,3)=  1.583333e-02;
    A(1,4)= -1.968056e-03; A(1,5)= 2.054167e-04;
    A(2,0)= -2.500000e+00; A(2,1)= 5.000000e-01; A(2,2)= -7.083333e-02; A(2,3)=  7.500000e-03;
    A(2,4)= -6.465278e-04; A(2,5)= 4.708333e-05;
    A(3,0)= -5.000000e+00; A(3,1)= 5.000000e-01; A(3,2)= -4.166667e-02; A(3,3)=  2.500000e-03;
    A(3,4)= -1.263889e-04; A(3,5)= 5.416667e-06;
    A(4,0)=  5.000000e+00; A(4,1)= 5.000000e-01; A(4,2)=  4.166667e-02; A(4,3)=  2.500000e-03;
    A(4,4)=  1.263889e-04; A(4,5)= 5.416667e-06;
    A(5,0)=  2.500000e+00; A(5,1)= 5.000000e-01; A(5,2)=  7.083333e-02; A(5,3)=  7.500000e-03;
    A(5,4)=  6.465278e-04; A(5,5)= 4.708333e-05;
    A(6,0)=  1.666667e+00; A(6,1)= 5.000000e-01; A(6,2)=  1.027778e-01; A(6,3)=  1.583333e-02;
    A(6,4)=  1.968056e-03; A(6,5)= 2.054167e-04;
    A(7,0)=  1.250000e+00; A(7,1)= 5.000000e-01; A(7,2)=  1.354167e-01; A(7,3)=  2.750000e-02;
    A(7,4)=  4.489931e-03; A(7,5)= 6.137500e-04;

    B(0) = (1/16.0630)* 4.166255e+01; B(1) =(1/16.0630)*1.967277e+01;
    B(2) = (1/16.0630)*(-4.440800e+00); B(3) =(1/16.0630)*(-4.781350e+01);
    B(4) = (1/16.0630)*7.203730e+01; B(5) =(1/16.0630)* 3.482060e+01;
    B(6) = (1/16.0630)*2.109503e+01; B(7) =(1/16.0630)*1.391725e+01;

    // Initialize Matlab solution
    XMatlab(0) =  0.84881784906768; XMatlab(1) = 1.31781097419130;
    XMatlab(2) = -12.13223678596093; XMatlab(3) = 37.95059077817579;
    XMatlab(4) = -62.75261698835658; XMatlab(5) = 44.82268878668197;
    RezNormMatlab = 2.701045096544557e-07;
    
    Solve_LS_Householder(A,B,X,krank,RezNorm);   

    // check condition
    for (int i=0; i<A.size(1); i++){
      SolutionDiff(i) = fabs(X(i)-XMatlab(i));
      ensure ("Solution",SolutionDiff(i)<=tol);
    }
    ensure ("Residual",fabs(RezNormMatlab-RezNorm)<=tol);
  }

  /* Test 13 */
  template<>
  template<>
  void LinearSystems_TestObject::test<13>()
  {
    DenseMatrix A(8,5); 	// Test matrix
    ColumnVector B(8);		// Test free term
    ColumnVector X(5);		// Solution vector
    ColumnVector XMatlab(5);	// Matlab solution vector
    ColumnVector SolutionDiff(5);	// Difference between the two solutions
    // Initialize variables
    tau = 10.0;
    tol = 1.0e-14;

    // Initialize test matrix & free term
    A(0,0)= -2.911958e-01; A(0,1)=  9.372864e-02; A(0,2)= -3.204096e-01;
    A(0,3)=  2.062637e-01; A(0,4)=  1.134784e-01;
    A(1,0)=  0.000000e+00; A(1,1)=  0.000000e+00; A(1,2)= -5.850551e-01;
    A(1,3)= -0.000000e+00; A(1,4)=  2.072070e-01;
    A(2,0)=  2.911958e-01; A(2,1)=  9.372864e-02; A(2,2)= -3.204096e-01;
    A(2,3)= -2.062637e-01; A(2,4)=  1.134784e-01;
    A(3,0)= -6.437500e-01; A(3,1)=  2.072070e-01; A(3,2)=  0.000000e+00;
    A(3,3)= -0.000000e+00; A(3,4)=  0.000000e+00;
    A(4,0)= -2.911958e-01; A(4,1)=  9.372864e-02; A(4,2)=  3.204096e-01;
    A(4,3)= -2.062637e-01; A(4,4)=  1.134784e-01; 
    A(5,0)=  0.000000e+00; A(5,1)=  0.000000e+00; A(5,2)=  5.850551e-01; 
    A(5,3)=  0.000000e+00; A(5,4)=  2.072070e-01; 
    A(6,0)=  2.911958e-01; A(6,1)=  9.372864e-02; A(6,2)=  3.204096e-01; 
    A(6,3)=  2.062637e-01; A(6,4)=  1.134784e-01; 
    A(7,0)=  6.437500e-01; A(7,1)=  2.072070e-01; A(7,2)=  0.000000e+00; 
    A(7,3)=  0.000000e+00; A(7,4)=  0.000000e+00; 

    B(0) = 2.181889e-02; B(1) = 3.984043e-02;
    B(2) = 2.181889e-02; B(3) = 0.000000e+00;
    B(4) =-7.982522e-04; B(5) =-1.457577e-03;
    B(6) =-7.982522e-04; B(7) = -2.220446e-16;

    // Initialize Matlab solution
    XMatlab(0) = -0.00000000000000; XMatlab(1) = -0.00000000759984;
    XMatlab(2) = -0.03529411408438; XMatlab(3) = 0.00000000000000;
    XMatlab(4) =  0.09261957652707; 
    RezNormMatlab = 6.846461854915191e-09;
    
    Solve_LS_Householder(A,B,X,krank,RezNorm);   

    // check condition
    for (int i=0; i<A.size(1); i++){
      SolutionDiff(i) = fabs(X(i)-XMatlab(i));
      ensure ("Solution",SolutionDiff(i)<=tol);
    }
    ensure ("Residual",fabs(RezNormMatlab-RezNorm)<=tol);
  }

  /* Test 14 */
  template<>
  template<>
  void LinearSystems_TestObject::test<14>()
  {
    DenseMatrix A(2,1); 	// Test matrix
    ColumnVector B(2);		// Test free term
    ColumnVector X(1);		// Solution vector
    ColumnVector XMatlab(1);	// Matlab solution vector
    ColumnVector SolutionDiff(1);	// Difference between the two solutions
    // Initialize variables
    tau = 10.0;
    tol = 1.0e-12;

    // Initialize test matrix & free term
    A(0,0)= -1.16226507551140710e+01;
    A(1,0)=  1.06140413102548088e+01;

    B(0) = -3.29961094096233819e+03;
    B(1) =  3.59331102526692803e+03;

    // Initialize Matlab solution
    XMatlab(0) = 3.087453711199031e+02;
    RezNormMatlab = 4.283127593574232e+02;
    
    Solve_LS_Householder(A,B,X,krank,RezNorm);   

    // check condition
    for (int i=0; i<A.size(1); i++){
      SolutionDiff(i) = fabs(X(i)-XMatlab(i));
      ensure ("Solution",SolutionDiff(i)<=tol);
    }
    ensure ("Residual",fabs(RezNormMatlab-RezNorm)<=tol);
  }

  /* Test 15 */
  template<>
  template<>
  void LinearSystems_TestObject::test<15>()
  {
    DenseMatrix A(4,3); 	// Test matrix
    ColumnVector B(4);		// Test free term
    ColumnVector X(3);		// Solution vector
    ColumnVector XMatlab(3);	// Matlab solution vector
    ColumnVector SolutionDiff(3);	// Difference between the two solutions
    // Initialize variables
    tau = 1.0;
    tol = 1.0e-14;             // The condition number of this matrix is k = 2.2682e+26

    // Initialize test matrix & free term

    A(0,0)= -2.758547e+25; A(0,1)= 2.251875e+24; A(0,2)= -1.953157e+23;
    A(1,0)= -8.827352e+26; A(1,1)= 3.603001e+25; A(1,2)= -1.838266e+24;
    A(2,0)=  3.919987e+02; A(2,1)= 1.599995e+01; A(2,2)=  8.163239e-01;
    A(3,0)=  1.224996e+01; A(3,1)= 9.999967e-01; A(3,2)=  8.673441e-02;

    B(0) = 0.000000e+00;
    B(1) = 0.000000e+00;
    B(2) = -5.042079e+09;
    B(3) = -1.575650e+08;

    // Initialize Matlab solution
  
    XMatlab(0) = -0.00509620036129*1.0e-36; XMatlab(1) = -0.12389727709673*1.0e-36;
    XMatlab(2) = 0.01518595433954*1.0e-36;
    RezNormMatlab = 5.044540352843458e+09;
    
    Solve_LS_Householder(A,B,X,krank,RezNorm);   

    // check condition
    for (int i=0; i<A.size(1); i++){
      SolutionDiff(i) = fabs(X(i)-XMatlab(i));
      ensure ("Solution",SolutionDiff(i)<=tol);
    }
    ensure ("Residual",fabs(RezNormMatlab-RezNorm)/RezNorm<=tol);
  }

  /* Test 16 */
  template<>
  template<>
  void LinearSystems_TestObject::test<16>()
  {
    DenseMatrix A(4,3); 	// Test matrix
    ColumnVector B(4);		// Test free term
    ColumnVector X(3);		// Solution vector
    ColumnVector XMatlab(3);	// Matlab solution vector
    ColumnVector SolutionDiff(3);	// Difference between the two solutions
    ColumnVector Constant(4);
    // Initialize variables
    tau = 10.0;
    tol = 1.0e-14;
    // Initialize test matrix & free term

    A(0,0)= -0.2759; A(0,1)= 0.0225; A(0,2)= -0.0020;
    A(1,0)= -8.8274; A(1,1)= 0.3603; A(1,2)= -0.0184;
    A(2,0)= 0.0; A(2,1)= 0.0; A(2,2)= 0.0;
    A(3,0)= 0.0; A(3,1)= 0.0; A(3,2)= 0.0;

    B(0) = 0.0;
    B(1) = 0.0;
    B(2) = -5.0421;
    B(3) = -1.5756;

    // Initialize Matlab solution
  
    XMatlab(0) = 0.0; XMatlab(1) = 0.0; XMatlab(2) = 0.0;
    RezNormMatlab = 5.04456243593039;
    
    Solve_LS_Householder(A,B,X,krank,RezNorm);   

    // check condition
    for (int i=0; i<A.size(1); i++){
      SolutionDiff(i) = fabs(X(i)-XMatlab(i));
      ensure ("Solution",SolutionDiff(i)<=tol);
    }

    ensure ("Residual",fabs(RezNormMatlab-RezNorm)/RezNorm<=0.1);
  }

  /* Test 17 */
  template<>
  template<>
  void LinearSystems_TestObject::test<17>()
  {

    DenseMatrix A(4,3); 	// Test matrix
    ColumnVector B(4);		// Test free term
    ColumnVector X(3);		// Solution vector
    ColumnVector XMatlab(3);	// Matlab solution vector
    ColumnVector SolutionDiff(3);	// Difference between the two solutions
    ColumnVector Constant(4);
    // Initialize variables
    tau = 10.0;
    tol = 1.0e-14;
    // Initialize test matrix & free term

    A(0,0)= -0.2759; A(0,1)= 0.0225; A(0,2)= -0.0020;
    A(1,0)= -8.8274; A(1,1)= 0.3603; A(1,2)= -0.0184;
    A(2,0)= 0.0; A(2,1)= 0.0; A(2,2)= 0.0;
    A(3,0)= 0.0; A(3,1)= 0.0; A(3,2)= 0.0;

    B(0) = 0.0;
    B(1) = 0.0;
    B(2) = -5.0421;
    B(3) = -1.5756;


    //     DenseMatrix B_copy(&B(0),B.size(),1,MV_Matrix_::ref);
    //     DenseMatrix X_copy(&X(0),X.size(),1,MV_Matrix_::ref);
    //     ColumnVector Rnorm(1);

    //     B_copy(3,0) = 20.2020;


    //     DenseMatrix B_transp(&B_copy(0,0),B_copy.size(1),B_copy.size(0),MV_Matrix_::ref);

    //     B_transp(0,2) = 1.11111;
    //     Print(B_transp);

    //     B_transp.permute_col(1,2);
    //     Print(B_transp);
    //     Print(B_copy);
    //     A.permute_col(0,2);

    //     A.permute_row(0,3);
    //     Print(A);

    // Initialize Matlab solution
  
    XMatlab(0) = 0.0; XMatlab(1) = 0.0; XMatlab(2) = 0.0;
    RezNormMatlab = 5.04456243593039;
    
    Solve_LS_Householder(A,B,X,krank,RezNorm);   
    //    Solve_LS_Householder(A,B_copy,X_copy,krank,Rnorm);   

    // check condition
    for (int i=0; i<A.size(1); i++){
      SolutionDiff(i) = fabs(X(i)-XMatlab(i));
      ensure ("Solution",SolutionDiff(i)<=tol);
    }

    //    ensure ("Residual",fabs(RezNormMatlab-RezNorm)/RezNorm<=tol);
  }

  /* Test 18:*/
  template<>
  template<>
  void LinearSystems_TestObject::test<18>()
  {
    integer NROW, NCOL;
    NROW = 4;
    NCOL = 3;

    DenseMatrix A(NROW,NCOL);	// Test matrix
    ColumnVector B(NROW);	// Test free term
    ColumnVector X(NCOL);		// Solution vector
    ColumnVector XMatlab(NCOL);	// Matlab solution vector
    ColumnVector SolutionDiff(NCOL);	// Difference between the two solutions
    // Initialize variables
    tol = 1.0e-14;

    // Initialize test matrix & free term
    A(0,0)= -2.093656e-01; A(0,1)= 2.114593e-01; A(0,2)= -1.512815e-01;
    A(1,0)= -9.900990e-01; A(1,1)= 5.000000e-01; A(1,2)= -2.104167e-01;
    A(2,0)=  9.900990e-01; A(2,1)= 5.000000e-01; A(2,2)=  2.104167e-01;
    A(3,0)=  4.950495e-01; A(3,1)= 5.000000e-01; A(3,2)= 3.577083e-01;

    B(0) = 4.364948e-01; B(1) = 4.353389e-16; B(2) = 0.000000e+00; B(3) = 0.000000e+00;

    // Initialize Matlab solution
    XMatlab(0) = 0.20223247422880; XMatlab(1) = 0.24027609285144;
    XMatlab(2) =-0.95158872955264;
    RezNormMatlab = 0.35212869025901;

    // solve the least-squares system with local subroutine
    //##    Solve_LS_Householder(A,B,X,krank,RezNorm);

    // solve the least-squares system with Lapack subroutine
    Solve_LS_Householder_F77(A,B,krank,NROW,NCOL);

    // check condition
    for (int i=0; i<(int)A.size(1); ++i){
      SolutionDiff(i) = fabs(B(i)-XMatlab(i))/(1.0 + fabs(XMatlab(i)));
      ensure ("System 4x3",SolutionDiff(i)<=tol);
    }

  }

};

namespace tut
{
  LinearSystems_TestGroup Solve_LS_Householder_Test("Solve_LS_Householder");
}
