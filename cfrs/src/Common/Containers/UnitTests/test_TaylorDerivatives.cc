/*!\file test_TaylorDerivatives.cc
   \brief Regression tests for TaylorDerivatives class.*/

/* Include required C++ libraries. */
// None

/* Using std namespace functions */

/* Include CFFC header files */
#include "TestData.h"
#include "include/TaylorDerivatives.h"
#include "CFD/Gaussian2DState.h"

using namespace std;

namespace tut
{

  // **********************************************
  //     TEST SUITE: TaylorDerivatives_TestSuite
  // **********************************************
  // Data used for testing
  // **********************
  class Data_TaylorDerivatives: public TestData{
  public:
    TaylorDerivativesContainer<OneD,double> DerivContainer1D;
    TaylorDerivativesContainer<TwoD,double> DerivContainer2D;
    TaylorDerivativesContainer<ThreeD,double> DerivContainer3D;
    DerivativeObj<OneD,double> Obj1D;
    DerivativeObj<TwoD,double> Obj2D;
    DerivativeObj<ThreeD,double> Obj3D;

    TaylorDerivativesContainer<TwoD,Gaussian2D_pState> GaussianHOD;
    TaylorDerivativesContainer<OneD,double> Container1D;

    Data_TaylorDerivatives(){
      set_test_suite_path("Common/Containers/UnitTests");
    }

    ~Data_TaylorDerivatives(){}
  };

  typedef test_group<Data_TaylorDerivatives> TaylorDerivatives_TestGroup;
  typedef TaylorDerivatives_TestGroup::object TaylorDerivatives_TestObject;


  /*Test1: */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<1>()
  {
    set_test_name("operator ==");

    TaylorDerivativesContainer<TwoD,double> TestDeriv;
    TaylorDerivativesContainer<TwoD,double> DerivContainer2D_New(2);
    DerivContainer2D.GenerateContainer(2);

    ensure("Equal containers", DerivContainer2D_New == DerivContainer2D);

    TaylorDerivativesContainer<ThreeD,double> DerivContainer3D_New(2);
    DerivContainer3D.GenerateContainer(2);
    ensure("Equal containers", DerivContainer3D_New == DerivContainer3D);
  }

  /*Test2: */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<2>()
  {
    set_test_name("operator =");

    DerivContainer2D.GenerateContainer(2);
    TaylorDerivativesContainer<TwoD,double> DerivContainer2D_New;

    DerivContainer2D_New = DerivContainer2D;
    ensure("Assignment operator", DerivContainer2D_New == DerivContainer2D);
  }

  /*Test3: */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<3>()
  {
    set_test_name("copy constructor");

    DerivContainer3D.GenerateContainer(3);
    TaylorDerivativesContainer<ThreeD,double> DerivContainer3D_New = DerivContainer3D;

    ensure("Copy Constructor()", DerivContainer3D_New == DerivContainer3D);
  }

  /*Test4: */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<4>()
  {
    set_test_name("copy constructor");

    DerivativeObj<TwoD,double> NewObj = Obj2D;
    ensure("Copy Constructor", !(NewObj != Obj2D));
  }

  /*Test5: */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<5>()
  {

    set_test_name("Assignment operator");

    DerivativeObj<TwoD,double> NewObj;
    NewObj = Obj2D;
    ensure("Assignment operator", NewObj == Obj2D);
  }

  /*Test6: */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<6>()
  {
    set_test_name("operator <<");
    RunRegression = ON;

    // Generate data
    DerivativeObj<TwoD,double> NewObj;
    NewObj = Obj2D;

    MasterFile = "DerivativeObj_Master.dat";
    CurrentFile = "Current.dat";

    if (RunRegression){
      Open_Output_File(CurrentFile);

      // output object
      out() << NewObj;

      // == check against Master
      RunRegressionTest("Output Operator", CurrentFile, MasterFile, 1.0e-12);

    } else {
      Open_Output_File(MasterFile);

      // output object
      out() << NewObj;
    }
  }

  /*Test7: SetPowers() functions */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<7>()
  {
    int p1, p2, p3;
    p1 = 1; p2 = 2; p3 = 4;
    std::vector<int> XYZ, XYZ3;
    XYZ.push_back(p1); XYZ.push_back(p2);
    DerivativeObj<TwoD,double> NewObj(XYZ,0);
    Obj2D.SetPowers(p1,p2);
    ensure("SetPowers()", Obj2D == NewObj);

    XYZ3 = XYZ; XYZ3.push_back(p3);
    DerivativeObj<ThreeD,double> NewObj3D(XYZ3,0);
    Obj3D.SetPowers(p1,p2,p3);
    ensure("SetPowers()", Obj3D == NewObj3D);
  }

  /*Test8: SetPowers() functions */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<8>()
  {
    int p1, p2, p3;
    p1 = 1; p2 = 2; p3 = 4;
    std::vector<int> XYZ, XYZ3;
    XYZ.push_back(p1); XYZ.push_back(p2);
    DerivativeObj<TwoD,double> NewObj(XYZ,0);
    Obj2D.SetPowers(p1,p2);
    ensure("SetPowers()", Obj2D == NewObj);

    XYZ3 = XYZ; XYZ3.push_back(p3);
    DerivativeObj<ThreeD,double> NewObj3D(XYZ3,0);
    Obj3D.SetPowers(p1,p2,p3);
    ensure("SetPowers()", Obj3D == NewObj3D);
  }

  /*Test9: IsPowerEqualTo() functions */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<9>()
  {
    int p1, p2, p3;
    p1 = 1; p2 = 2; p3 = 4;
    double val = 10.1212;

    Obj2D.SetPowers(p1,p2);
    Obj2D.SetValue(val);

    ensure("IsPowerEqualTo()", Obj2D.IsPowerEqualTo(p1,p2));

    Obj3D.SetPowers(p1,p2,p3);
    ensure("IsPowerEqualTo()", Obj3D.IsPowerEqualTo(p1,p2,p3));
  }

  /*Test10: D() functions */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<10>()
  {
    int p1, p2, p3;
    p1 = 1; p2 = 2; p3 = 4;
    double val = 10.1212;

    Obj2D.SetPowers(p1,p2);
    Obj2D.SetValue(val);

    ensure("D()", Obj2D.D() == val);

    Obj3D.SetPowers(p1,p2,p3);
    Obj3D.SetValue(val*10);
    ensure("D()", Obj3D.D() == val*10);
  }


  /*Test11: P1(), P2(), P3() functions */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<11>()
  {
    int p1, p2, p3;
    p1 = 1; p2 = 2; p3 = 4;
    double val = 10.1212;

    Obj1D.SetPowers(p1);
    ensure("P1()", Obj1D.P1() == p1);

    Obj2D.SetPowers(p1,p2);
    Obj2D.SetValue(val);

    ensure("P1()", Obj2D.P1() == p1);
    ensure("P2()", Obj2D.P2() == p2);
    
    Obj3D.SetPowers(p1,p2,p3);
    Obj3D.SetValue(val*10);
    ensure("P1()", Obj3D.P1() == p1);
    ensure("P2()", Obj3D.P2() == p2);
    ensure("P3()", Obj3D.P3() == p3);

  }

  /*Test12: TaylorDerivativesContainer:: << operator */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<12>()
  {
    DerivContainer3D.GenerateContainer(3);
    //    cout << DerivContainer3D;
    DerivContainer2D.GenerateContainer(3);
    //    cout << DerivContainer2D;
    DerivContainer1D.GenerateContainer(3);
    //    cout << DerivContainer1D;

  }

  /*Test13: TaylorDerivativesContainer::operator(int) */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<13>()
  {
    int p1, p2, p3;
    p1 = 3; p2 = 2; p3 = 4;
    double val = 10;
    DerivContainer1D.GenerateContainer(3);
    DerivContainer1D(p1) = val;
    ensure("operator(int)", DerivContainer1D(p1) == val); 

    DerivContainer2D.GenerateContainer(5);
    DerivContainer2D(p1,p2) = val;
    ensure("operator(int,int)", DerivContainer2D(p1,p2) == val); 

    DerivContainer3D.GenerateContainer(9);
    DerivContainer3D(p1,p2,p3) = val;
    DerivContainer3D(0,0,1) = val;
    DerivContainer3D(1,1,7) = 12.2323;
    DerivContainer3D(1,0,8) = DerivContainer3D(9,0,0) = 56.454;

    ensure("operator(int,int,int)", DerivContainer3D(p1,p2,p3) ==  val);
    ensure("operator(int,int,int)", DerivContainer3D(0,0,1) ==  val);
    ensure("operator(int,int,int)", DerivContainer3D(1,1,7) ==  12.2323);
    ensure("operator(int,int,int)", DerivContainer3D(1,0,8) ==  56.454);
  }

  /*Test14: TaylorDerivativesContainer:: GaussianHOD */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<14>()
  {
    double val = 0.232345;
    GaussianHOD.GenerateContainer(3);
    GaussianHOD(1,2)[1] = val;
    GaussianHOD(3,0)[NUM_VAR_GAUSSIAN2D] = val;

    ensure("GaussianHOD", GaussianHOD(1,2)[1] == val );
    ensure("GaussianHOD", GaussianHOD(3,0)[NUM_VAR_GAUSSIAN2D] == val);
  }

  /*Test15: */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<15>()
  {

    set_test_name("ComputeSolutionFor(DeltaX)");

    // Set test data
    double DeltaX = 0.1;
    int RO;
    RO = 3;
    DerivContainer1D.GenerateContainer(RO);
    for (int i=0; i<=RO; i++)
      DerivContainer1D(i) = i+1;

    // == check
    ensure_equals("ComputeSolutionFor(double)", DerivContainer1D.ComputeSolutionFor(DeltaX), 1.234);
  }

  /*Test16: */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<16>()
  {

    set_test_name("ComputeSolutionFor(DeltaX, DeltaY)");

    double DeltaX = 1.2333333; double DeltaY = 0.1235;
    int RO;
    RO = 4;
    double result(0.0);
    DerivContainer2D.GenerateContainer(RO);
    DerivContainer2D.Limiter() = 0.343333427;

    for (int i=0; i<=RO; i++)
      for (int j=0; j<=RO-i; j++){
	DerivContainer2D(i,j) = i+j+1;
	if( (i==0)&&(j==0)){
	  result += (i+j+1) * pow(DeltaX,i) * pow(DeltaY,j);
	} else{
	  result += DerivContainer2D.Limiter() * (i+j+1) * pow(DeltaX,i) * pow(DeltaY,j);
	}
      }

    ensure_distance("ComputeSolutionFor(double)", DerivContainer2D.ComputeSolutionFor(DeltaX,DeltaY), result, 5.0e-14);
  }

  /*Test17: */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<17>()
  {

    set_test_name("ComputeSolutionFor(DeltaX, DeltaY, DeltaZ)");

    double DeltaX = 0.1; double DeltaY = 0.2; double DeltaZ = 0.3;
    int RO;
    RO = 4;
    double result = 0.0;
    DerivContainer3D.GenerateContainer(RO);
    for (int i=0; i<=RO; i++)
      for (int j=0; j<=RO-i; j++)
	for (int k=0; k<=RO-i-j; k++){
	  DerivContainer3D(i,j,k) = (i+j+k+1)*0.2;
	  result += (i+j+k+1)*0.2 * pow(DeltaX,i) * pow(DeltaY,j) * pow(DeltaZ,k);
      }

    ensure_distance("ComputeSolutionFor(double)", DerivContainer3D.ComputeSolutionFor(DeltaX,DeltaY,DeltaZ), result, 5.0e-14);
  }

  /*Test18: ComputeSolutionFor(DeltaX, DeltaY, DeltaZ) */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<18>()
  {
    set_test_name("ComputeSolutionFor(DeltaX, DeltaY, DeltaZ)");

    double DeltaX = 0.1; double DeltaY = 0.2; double DeltaZ = 0.3;
    int RO;
    RO = 4;
    double result = 0.0;
    TaylorDerivativesContainer<ThreeD, Gaussian2D_pState> GaussianHOD_;
    GaussianHOD_.GenerateContainer(RO);
    for (int i=0; i<=RO; i++)
      for (int j=0; j<=RO-i; j++)
	for (int k=0; k<=RO-i-j; k++){
	  GaussianHOD_(i,j,k)[3] = (i+j+k+1)*0.2;
	  result += (i+j+k+1)*0.2 * pow(DeltaX,i) * pow(DeltaY,j) * pow(DeltaZ,k);
      }


    ensure_distance("ComputeSolutionFor(double)", GaussianHOD_.ComputeSolutionFor(DeltaX,DeltaY,DeltaZ)[3], result, 5.0e-14);
  }

  /****************************************************************
   *     Tests for the TaylorDerivativesContainer::iterator       *
   ***************************************************************/

  /*Test19: Constructor */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<19>()
  {
#ifdef __Use_Iterator__
    TaylorDerivativesContainer<OneD,double>::iterator Iter = DerivContainer1D;
    TaylorDerivativesContainer<OneD,double> DC1D;
    TaylorDerivativesContainer<OneD,double>::iterator Iter2 = DC1D;

    int RO,i;
    RO = 4;
    DerivContainer1D.GenerateContainer(RO);
    for (i=0; i<=RO; i++)
      DerivContainer1D(i) = i+1.0;

    // check the values by moving through the container using the interator
    for (i=0,Iter = DerivContainer1D.begin(); Iter!=DerivContainer1D.end(); Iter++, i++){
      ensure("TaylorDerivativesContainer::iterator", Iter->D() == (i+1.0));
    }

    DC1D.GenerateContainer(RO-1);
    for (int i=0; i<=RO-1; i++)
      DC1D(i) = i*3.14;

    Iter = Iter2;
    for (i=0, Iter = DC1D.begin(); Iter!=DC1D.end(); Iter++ , i++){
      ensure("TaylorDerivativesContainer::iterator", Iter->D() == (i*3.14));
    }
#endif
  }

  /*Test20: Jump operators */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<20>()
  {
#ifdef __Use_Iterator__
    TaylorDerivativesContainer<OneD,double>::iterator Iter = DerivContainer1D;
    TaylorDerivativesContainer<OneD,double>::iterator Iter2 = DerivContainer1D;
    int RO,i;
    RO = 4;
    DerivContainer1D.GenerateContainer(RO);
    for (i=0; i<=RO; i++)
      DerivContainer1D(i) = i+1.0;

    Iter = DerivContainer1D.back();
    ensure("Iter.back()", Iter->P1() == RO);

    Iter--;
    ensure("Iter--", Iter->P1() == RO-1 );

    Iter-=1;
    ensure("Iter-=1", Iter->P1() == RO-2 );

    Iter2++;
    ensure("Compare iterators", Iter > Iter2);
    ensure("Compare iterators", Iter >= Iter2);

    Iter2++;
    ensure("Compare iterators", Iter == Iter2);

    Iter2++;
    ensure("Compare iterators", Iter < Iter2);
    ensure("Compare iterators", Iter <= Iter2);
    ensure("Compare iterators", Iter != Iter2);
#endif
  }

  /*Test21: TaylorDerivatives::ComputeSolutionFor(DeltaX) */
  template<>
  template<>
  void TaylorDerivatives_TestObject::test<21>()
  {
    double DeltaX = 0.1;
    int RO;
    RO = 3;
    Container1D.GenerateContainer(RO);
    for (int i=0; i<=RO; ++i)
      Container1D(i) = i+1;

    Container1D.Limiter() = 0.5;

    //    ensure("ComputeSolutionFor(double)", Container1D.ComputeSolutionFor(DeltaX) == 1.234 );
    ensure_equals("ComputeSolutionFor(double)", Container1D.ComputeSolutionFor(DeltaX), 1.117 );

    TaylorDerivativesContainer<OneD,double> SecondContainer1D = Container1D;
    TaylorDerivativesContainer<OneD,double> ThirdContainer1D;

    ThirdContainer1D = Container1D;

    ensure("Copy Constructor", SecondContainer1D == Container1D);
    ensure("Assignment Operator", ThirdContainer1D == SecondContainer1D);

    ThirdContainer1D = SecondContainer1D;
    ensure("Assignment Operator", ThirdContainer1D == Container1D);

  }

}

tut::TaylorDerivatives_TestGroup TaylorDerivatives_Test("Template Class:TaylorDerivatives");
