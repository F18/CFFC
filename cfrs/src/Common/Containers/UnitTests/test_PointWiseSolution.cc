/*!\file test_PointWiseSolution.cc
   \brief Regression tests for PointWiseSolution class.*/

/* Include required C++ libraries. */
// None

/* Using std namespace functions */

/* Include CFFC header files */
#include "TestData.h"
#include "../PointWiseSolution.h"
#include "Grid/Grid2D/Cell2D.h"
#include "CFD/Gaussian2DState.h"

using namespace std;

namespace tut
{

  // **********************************************
  //     TEST SUITE: PointWiseSolution_TestSuite
  // **********************************************
  // Data used for testing
  // **********************
  class Data_PointWiseSolution: public TestData {
  public:
    PointWiseSolution<Node2D,double> SolutionPoint;
    PointWiseSolution<Node2D,Gaussian2D_pState> G_Solution;

    Data_PointWiseSolution(void){
      set_test_suite_path("Common/Containers/UnitTests");
    }

    ~Data_PointWiseSolution(){}

  };
  
  typedef test_group<Data_PointWiseSolution> PointWiseSolution_TestSuite;
  typedef PointWiseSolution_TestSuite::object  PWS_Object;
  
  
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
  void PWS_Object::test<1>()
  {
    set_test_name("Assignment operator");

    Node2D Node(1,1);
    PointWiseSolution<Node2D> A(Node, 10);
    PointWiseSolution<Node2D> B;
    B = A;
    double ValB = B.GetValue();
    double ValA = A.GetValue();
    ensure("Assignment & equal operators", A == B);
    ensure("Assignment & equal operators", ValA == ValB);
  }

  /* Test 2:*/
  template<>
  template<>
  void PWS_Object::test<2>()
  {

    set_test_name("SetValue() member function");

    Node2D Node(1,1);
    double value = 10.0;
    PointWiseSolution<Node2D> A(Node, value);
    PointWiseSolution<Node2D> B;

    // Set values using member functions
    B.SetNode(Node);
    B.SetValue(value);

    ensure("SetValue() & SetNode()", A == B);
  }

  /* Test 3:*/
  template<>
  template<>
  void PWS_Object::test<3>()
  {

    set_test_name("operator << ");
    RunRegression = ON;

    // Generate data
    Node2D Node(1,1);
    double value = 10.0;
    PointWiseSolution<Node2D> A(Node, value);

    MasterFile = "PointWiseSolution_Master.dat";
    CurrentFile = "Current.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // output object
      out() << A;

      // == check against Master
      RunRegressionTest("Output Operator", CurrentFile, MasterFile, 1.0e-12);

    } else {
      Open_Output_File(MasterFile);

      // output object
      out() << A;
    }
  }

  /* Test 4: */
  template<>
  template<>
  void PWS_Object::test<4>()
  {

    set_test_name("Copy Constructor");

    // Generate data
    Node2D Node(1,1);
    Gaussian2D_pState GS;
    GS[1] = 1;
    GS[2] = 2;
    GS[3] = 3;
    G_Solution.SetNode(Node);
    G_Solution.SetValue(GS);

    // Use copy constructor
    PointWiseSolution<Node2D, Gaussian2D_pState> A = G_Solution;

    // == check
    ensure("Copy Constructor", A == G_Solution);
  }

  /* Test 5:*/
  template<>
  template<>
  void PWS_Object::test<5>()
  {
    set_test_name(" Assignment Operator");

    // Generate data
    Node2D Node(1,1);
    Gaussian2D_pState GS;
    GS[1] = 1;
    GS[2] = 2;
    GS[3] = 3;
    G_Solution.SetNode(Node);
    G_Solution.SetValue(GS);

    // Use assigment operator 
    PointWiseSolution<Node2D, Gaussian2D_pState> A;
    A = G_Solution;

    // == check
    ensure("Assignment Operator: Test1", A == G_Solution);

    // == check
    GS[1] = 2;
    A.SetValue(GS);
    ensure("Assignment Operator: Test2", A != G_Solution);
  }

  /* Test 6:*/
  template<>
  template<>
  void PWS_Object::test<6>()
  {

    set_test_name(" Assignment Operator");

    // Generate data
    PointWiseSolution<Node3D,Gaussian2D_pState> G_Solution3D;
    Node3D Node(1,1,0.5);
    Gaussian2D_pState GS;
    GS[1] = 1;
    GS[2] = 2;
    GS[3] = 3;
    G_Solution3D.SetNode(Node);
    G_Solution3D.SetValue(GS);

    // use operator 
    PointWiseSolution<Node3D, Gaussian2D_pState> A;
    A = G_Solution3D;

    // == check
    ensure("Copy Constructor", A == G_Solution3D);
  }

}

/************************
 *  CREATE TEST SUITES  *
 ************************/

tut::PointWiseSolution_TestSuite PointWiseSolution_Test("Template Class:PointWiseSolution");

