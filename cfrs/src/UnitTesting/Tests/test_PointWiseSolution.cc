#include <tut.h>
#include "Common/Containers/PointWiseSolution.h"
#include "Grid/Grid2D/Cell2D.h"
#include "CFD/Gaussian2DState.h"


using namespace std;

namespace tut

{
  // Data used for testing
  struct Data_PointWiseSolution{
    PointWiseSolution<Node2D,double> SolutionPoint;
    PointWiseSolution<Node2D,Gaussian2D_pState> G_Solution;
  };

   /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   */
  typedef test_group<Data_PointWiseSolution> PointWiseSolution_TestGroup;
  typedef PointWiseSolution_TestGroup::object  PointWiseSolution_TestObject;
}

/**********************************************
 *                THE TESTS                   *
 *********************************************/

namespace tut
{

  /* Test 1: Assignment operator*/
  template<>
  template<>
  void PointWiseSolution_TestObject::test<1>()
  {
    Node2D Node(1,1);
    PointWiseSolution<Node2D> A(Node, 10);
    PointWiseSolution<Node2D> B;
    B = A;
    double ValB = B.GetValue();
    double ValA = A.GetValue();
    ensure("Assignment & equal operators", A == B);
    ensure("Assignment & equal operators", ValA == ValB);
  }

  /* Test 2: SetValue() member function*/
  template<>
  template<>
  void PointWiseSolution_TestObject::test<2>()
  {
    Node2D Node(1,1);
    double value = 10.0;
    PointWiseSolution<Node2D> A(Node, value);
    PointWiseSolution<Node2D> B;

    // Set values using member functions
    B.SetNode(Node);
    B.SetValue(value);

    ensure("SetValue() & SetNode()", A == B);
  }

  /* Test 3: Print on the screen*/
  template<>
  template<>
  void PointWiseSolution_TestObject::test<3>()
  {
    Node2D Node(1,1);
    double value = 10.0;
    PointWiseSolution<Node2D> A(Node, value);
    fail("cout test -> This test is not active");
    //    cout << A ;
  }

  /* Test 4: Copy Constructor*/
  template<>
  template<>
  void PointWiseSolution_TestObject::test<4>()
  {
    Node2D Node(1,1);
    Gaussian2D_pState GS;
    GS[1] = 1;
    GS[2] = 2;
    GS[3] = 3;
    G_Solution.SetNode(Node);
    G_Solution.SetValue(GS);

    PointWiseSolution<Node2D, Gaussian2D_pState> A = G_Solution;

    ensure("Copy Constructor", A == G_Solution);
  }

  /* Test 5: Assignment Operator*/
  template<>
  template<>
  void PointWiseSolution_TestObject::test<5>()
  {
    Node2D Node(1,1);
    Gaussian2D_pState GS;
    GS[1] = 1;
    GS[2] = 2;
    GS[3] = 3;
    G_Solution.SetNode(Node);
    G_Solution.SetValue(GS);

    PointWiseSolution<Node2D, Gaussian2D_pState> A;

    A = G_Solution;

    ensure("Assignment Operator: Test1", A == G_Solution);

    GS[1] = 2;
    A.SetValue(GS);
    ensure("Assignment Operator: Test2", A != G_Solution);
  }

  /* Test 6: Assignment Operator*/
  template<>
  template<>
  void PointWiseSolution_TestObject::test<6>()
  {
    PointWiseSolution<Node3D,Gaussian2D_pState> G_Solution3D;
    Node3D Node(1,1,0.5);
    Gaussian2D_pState GS;
    GS[1] = 1;
    GS[2] = 2;
    GS[3] = 3;
    G_Solution3D.SetNode(Node);
    G_Solution3D.SetValue(GS);

    PointWiseSolution<Node3D, Gaussian2D_pState> A;

    A = G_Solution3D;

    ensure("Copy Constructor", A == G_Solution3D);
  }


}

namespace tut
{
  PointWiseSolution_TestGroup PointWiseSolution_Test("PointWiseSolution_Test");
}
