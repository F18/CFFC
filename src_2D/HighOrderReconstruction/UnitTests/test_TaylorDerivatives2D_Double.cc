/*!\file test_TaylorDerivatives2D.cc
  \brief Regression tests for the template (container) class TaylorDerivatives for 2D. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../TaylorDerivatives.h"


namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  struct Data_TaylorDerivatives2D: public TestData {

    // Public data and member functions
  public:
    
    typedef TaylorDerivativesContainer<TwoD,double> TD_Type;
    typedef DerivativeObj<TwoD,double> Deriv_Type;
    
    // DerivativeObj
    Deriv_Type DObj;
    
    // TaylorDerivativesContainer
    TD_Type TD, TD_1, TD0, TD1, TD2;

    // Default Constructor
    Data_TaylorDerivatives2D(void);

    // Destructor
    ~Data_TaylorDerivatives2D(void){};

    // Initialize container
    void InitContainer(TD_Type & Obj);

    // Initalization
    double Initialization(const int & i){ return (i+1)*3.03; }

    // EnsureValues --> Values of left are equal to values of right
    void EnsureValues(TD_Type & left, TD_Type & right);
    
  };
  

  Data_TaylorDerivatives2D::Data_TaylorDerivatives2D(void){
    // Initialize DObj
    DObj.SetPowers(1,3);
    DObj.D() = 10.343434;
    
    // Initialize Container2D
    // GenerateContainer
    TD0.GenerateContainer(0);
    TD1.GenerateContainer(1);
    TD2.GenerateContainer(2);

    // Set paths
    set_test_suite_path("HighOrderReconstruction/UnitTests");
    set_local_output_path("TaylorDerivatives2D_Data");

  }

  void Data_TaylorDerivatives2D::InitContainer(TD_Type & Obj){
    for (int i = 0; i<=Obj.LastElem(); ++i){
      Obj(i).D() = Initialization(i);
    }
  }

  void Data_TaylorDerivatives2D::EnsureValues(TD_Type & left, TD_Type & right){
    char msg[20], basic[10], number[3];
    strcpy(basic,"Entry #");

    ensure_equals("EnsureValues()--> size()", left.size(), right.size());

    for (int i=left.FirstElem(); i<=left.LastElem(); ++i){
      strcpy(msg,basic);
      sprintf (number, "%d", i);
      strcat(msg,number);
      ensure_equals(msg,left(i).D(), right(i).D());
    }
  }

  typedef test_group<Data_TaylorDerivatives2D> TaylorDerivatives2D_TestGroup;
  typedef TaylorDerivatives2D_TestGroup::object object;


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


  /*******************************************************
   * Tests for CLASS Template: DerivativeObj             *
   ******************************************************/

  /*Test1: */
  template<>
  template<>
  void object::test<1>()
  {
    set_test_name("DerivativeObj-> Constructor");
    ensure_equals("Power1", DObj.P1(), 1);
    ensure_equals("Power2", DObj.P2(), 3);
    ensure_equals("Value", DObj.D(),10.343434);

  }

  /*Test2: */
  template<>
  template<>
  void object::test<2>()
  {
    set_test_name("DerivativeObj-> Copy Constructor");

    Deriv_Type DObj2(DObj);

    ensure_equals("Power1", DObj2.P1(), DObj.P1());
    ensure_equals("Power2", DObj2.P2(), DObj.P2());
    ensure_equals("Value", DObj2.D(), DObj.D());
  }

  /*Test3: */
  template<>
  template<>
  void object::test<3>()
  {
    set_test_name("DerivativeObj-> Assignment operator");

    Deriv_Type DObj2;
    DObj2 = DObj;

    ensure_equals("Power1", DObj2.P1(), DObj.P1());
    ensure_equals("Power2", DObj2.P2(), DObj.P2());
    ensure_equals("Value", DObj2.D(), DObj.D());
  }

  /*Test4: */
  template<>
  template<>
  void object::test<4>()
  {
    set_test_name("DerivativeObj-> operator ==");

    Deriv_Type DObj2;
    DObj2 = DObj;

    ensure_equals("==", DObj2 == DObj, true);
  }

  /*Test5: */
  template<>
  template<>
  void object::test<5>()
  {
    set_test_name("DerivativeObj-> operator !=");

    Deriv_Type DObj2(DObj);
    DObj2.P1() = 2;
    ensure_equals("1: !=", DObj2 != DObj, true);

    DObj2.P2() = 2;
    ensure_equals("2: !=", DObj2 != DObj, true);

    DObj2.D() = 1002E-2;
    ensure_equals("3: !=", DObj2 != DObj, true);
  }

  /*Test6: */
  template<>
  template<>
  void object::test<6>()
  {
    set_test_name("DerivativeObj-> IsPowerEqualTo()");

    Deriv_Type DObj2(DObj);
    ensure_equals("1: IsPowerEqualTo()", DObj2.IsPowerEqualTo(1,3), true);

    DObj2.P1() = 2;
    ensure_equals("2: IsPowerEqualTo()", DObj2.IsPowerEqualTo(2,3), true);

    ensure_equals("3: IsPowerEqualTo()", DObj2.IsPowerEqualTo(0,0), false);
  }

  /*Test7: */
  template<>
  template<>
  void object::test<7>()
  {
    set_test_name("DerivativeObj-> D(int)");

    ensure_equals("1: Val access", DObj.D(0), 10.343434);
    ensure_equals("2: Val access", DObj.D(10), 10.343434);
    ensure_equals("3: Val access", DObj.D(-10), 10.343434);

  }

  /*Test8: */
  template<>
  template<>
  void object::test<8>()
  {
    set_test_name("DerivativeObj-> operator << & >>");
    Deriv_Type DObj2(1,1,1.0); // this object is different than DObj!

    Check_Input_Output_Operator(DObj2);
  }

  /*Test9: */
  template<>
  template<>
  void object::test<9>()
  {
    set_test_name("DerivativeObj-> Read_Derivative");
    Deriv_Type DObj2;

    // open the output file
    Open_Output_File("DObj1");
    // output data
    output_file << DObj << endl;

    // open the output file for reading
    Open_Output_File_For_Reading();

    // read data
    DObj2.Read_Derivative(in_file);

    ensure_equals("Read_Derivative", DObj2, DObj);

    Remove_Output_File();
  }




  /**************************************************************************
   * Tests for CLASS Templete: TaylorDerivatives 2D Container Declaration   *
   *************************************************************************/

  /*Test10: */
  template<>
  template<>
  void object::test<10>()
  {
    set_test_name("TaylorDeriv-> Default Constructor");
    ensure_equals("size()", TD.size(), 0);
    ensure_equals("Limiter()", TD.Limiter(), 1.0);
    ensure_equals("RecOrder()", TD.RecOrder(), -1.0);
    ensure_equals("FirstElem()", TD.FirstElem(), 0);
    ensure_equals("LastElem()", TD.LastElem(), -1);
  }

  /*Test11: */
  template<>
  template<>
  void object::test<11>()
  {
    set_test_name("TaylorDeriv-> Constructor(int)");

    TD_Type TD3(3);
    ensure_equals("size()", TD3.size(), 10);
    ensure_equals("Limiter()", TD3.Limiter(), 1.0);
    ensure_equals("RecOrder()", TD3.RecOrder(), 3);
  }

  /*Test12: */
  template<>
  template<>
  void object::test<12>()
  {
    set_test_name("TaylorDeriv-> GenerateContainer ");

    ensure_equals("2: size()", TD0.size(), 1);
    ensure_equals("2: Limiter()", TD0.Limiter(), 1.0);
    ensure_equals("2: RecOrder()", TD0.RecOrder(), 0);

    ensure_equals("2: size()", TD1.size(),3);
    ensure_equals("2: Limiter()", TD1.Limiter(), 1.0);
    ensure_equals("2: RecOrder()", TD1.RecOrder(),1);

    ensure_equals("3: size()", TD2.size(), 6);
    ensure_equals("3: Limiter()", TD2.Limiter(), 1.0);
    ensure_equals("3: RecOrder()", TD2.RecOrder(), 2);
  }

  /*Test13: */
  template<>
  template<>
  void object::test<13>()
  {
    set_test_name("TaylorDeriv-> Copy constructor");

    InitContainer(TD2);
    TD2.Limiter() = 0.232345;
    TD_Type TD_Copy(TD2);

    ensure_equals("size()", TD_Copy.size(), TD2.size());
    ensure_equals("Limiter()", TD_Copy.Limiter(), TD2.Limiter());
    ensure_equals("RecOrder()", TD_Copy.RecOrder(), TD2.RecOrder());
    ensure_equals("FirstElem()", TD_Copy.FirstElem(), TD2.FirstElem());
    ensure_equals("LastElem()", TD_Copy.LastElem(), TD2.LastElem());

    // check the values too
    EnsureValues(TD_Copy,TD2);
  }

  /*Test14: */
  template<>
  template<>
  void object::test<14>()
  {
    set_test_name("TaylorDeriv-> Assignment operator");

    InitContainer(TD2);
    TD2.Limiter() = 0.232345;
    TD_Type TD_Copy;

    // operator
    TD_Copy = TD2;

    // check 
    ensure_equals("size()", TD_Copy.size(), TD2.size());
    ensure_equals("Limiter()", TD_Copy.Limiter(), TD2.Limiter());
    ensure_equals("RecOrder()", TD_Copy.RecOrder(), TD2.RecOrder());
    ensure_equals("FirstElem()", TD_Copy.FirstElem(), TD2.FirstElem());
    ensure_equals("LastElem()", TD_Copy.LastElem(), TD2.LastElem());

    // check the values too
    EnsureValues(TD_Copy,TD2);
  }

  /*Test15: */
  template<>
  template<>
  void object::test<15>()
  {
    set_test_name("TaylorDeriv-> IndexOrder() ");
    TD_Type TD3(3);

    ensure_equals("pair (0,0)", TD3.IndexOrder(0,0), 0);
    ensure_equals("pair (0,1)", TD3.IndexOrder(0,1), 1);
    ensure_equals("pair (0,2)", TD3.IndexOrder(0,2), 2);
    ensure_equals("pair (0,3)", TD3.IndexOrder(0,3), 3);

    ensure_equals("pair (1,0)", TD3.IndexOrder(1,0), 4);
    ensure_equals("pair (1,1)", TD3.IndexOrder(1,1), 5);
    ensure_equals("pair (1,2)", TD3.IndexOrder(1,2), 6);

    ensure_equals("pair (2,0)", TD3.IndexOrder(2,0), 7);
    ensure_equals("pair (2,1)", TD3.IndexOrder(2,1), 8);
    ensure_equals("pair (3,0)", TD3.IndexOrder(3,0), 9);

    // passing powers out of range is not detected (the code is faster this way)
  }


  /*Test16: */
  template<>
  template<>
  void object::test<16>()
  {
    set_test_name("TaylorDeriv-> ResetLimiter() ");
    TD.Limiter() = 0.5;
    ensure_equals("Val", TD.Limiter(), 0.5);
    TD.ResetLimiter();
    ensure_equals("Limiter=1.0", TD.Limiter(3), 1.0);
  }

  /*Test17: */
  template<>
  template<>
  void object::test<17>()
  {
    set_test_name("TaylorDeriv-> operator == ");
    TD_Type TD3(3);

    // Initialize Container
    InitContainer(TD3);

    // Assignment
    TD2 = TD3;

    // Check
    ensure_equals("==", TD2 == TD3, true);
  }

  /*Test18: */
  template<>
  template<>
  void object::test<18>()
  {
    set_test_name("TaylorDeriv-> operator != ");
    TD_Type TD3(3);

    // Assignment
    TD2 = TD3;
    TD3.GenerateContainer(2);
    // Check
    ensure_equals("1: !=", TD2 != TD3, true);

    // Initialize Container
    InitContainer(TD3);
    TD2 = TD3;
    TD3(2,0) = 10.232323;

    // Check
    ensure_equals("2: !=", TD2 != TD3, true);

  }

  /*Test19: */
  template<>
  template<>
  void object::test<19>()
  {
    set_test_name("TaylorDeriv-> : Out->In->Compare ");

    TD_Type TD3(3);
    // Initialize Container
    InitContainer(TD3);
    Check_Input_Output_Operator(TD3);
  }


  /*Test20: */
  template<>
  template<>
  void object::test<20>()
  {
    set_test_name("TaylorDeriv-> : ComputeSolutionFor() ");

    double DeltaX, DeltaY, Result;
    int RO;			// Reconstruction Order
    char number[4], msg[20];

    // Input Data
    DeltaX = 1.2333333;
    DeltaY = 0.1235;

    for (RO = 4; RO>=1; --RO){
      TD.GenerateContainer(RO);

      // Initialize container
      InitContainer(TD);
      
      // Get analytic result
      Result = 0.0;
      for (int i = 0; i<=TD.LastElem(); ++i){
	Result += Initialization(i) * pow(DeltaX,TD(i).P1()) * pow(DeltaY,TD(i).P2());
      }

      // check
      sprintf(number,"%d",RO);
      strcpy(msg, "Order ");
      strcat(msg, number);
      ensure_distance(msg, TD.ComputeSolutionFor(DeltaX, DeltaY), Result, tol*(1.0 + Result));
    }

  }


  /*Test21: */
  template<>
  template<>
  void object::test<21>()
  {
    set_test_name("TaylorDeriv-> : ComputeSolutionFor() ");

    double DeltaX, DeltaY, Result, Limiter;
    int RO;			// Reconstruction Order
    char number[4], msg[20];

    // Input Data
    DeltaX = 1.2333333;
    DeltaY = 0.1235;
    Limiter = 0.3333;

    for (RO = 4; RO>=1; --RO){
      TD.GenerateContainer(RO);
      TD.Limiter() = Limiter;

      // Initialize container
      InitContainer(TD);
      
      // Get analytic result
      Result = 0.0;
      for (int i = TD.FirstElem()+1; i<=TD.LastElem(); ++i){
	Result += Initialization(i) * pow(DeltaX,TD(i).P1()) * pow(DeltaY,TD(i).P2());
      }
      Result = Initialization(0) + Limiter*Result;

      // check
      sprintf(number,"%d",RO);
      strcpy(msg, "Order ");
      strcat(msg, number);
      ensure_distance(msg, TD.ComputeSolutionFor(DeltaX, DeltaY), Result, tol*Result);
    }

  }

  /*Test22: */
  template<>
  template<>
  void object::test<22>()
  {
    set_test_name("TaylorDeriv-> : ComputeXGradientFor() ");

    double DeltaX, DeltaY, Result;
    int RO;			// Reconstruction Order
    char number[4], msg[20];

    // Input Data
    DeltaX = 1.2333333;
    DeltaY = 0.1235;

    for (RO = 4; RO>=1; --RO){
      TD.GenerateContainer(RO);

      // Initialize container
      InitContainer(TD);
      
      // Get analytic result
      Result = 0.0;
      for (int i = 0; i<=TD.LastElem(); ++i){
	Result += ( Initialization(i) * TD(i).P1()*pow(DeltaX,TD(i).P1()-1.0) *
		    pow(DeltaY,TD(i).P2()) );
      }

      // check
      sprintf(number,"%d",RO);
      strcpy(msg, "Order ");
      strcat(msg, number);
      ensure_distance(msg, TD.ComputeXGradientFor(DeltaX, DeltaY), Result, tol*(1.0 + Result));
    }

  }

  /*Test23: */
  template<>
  template<>
  void object::test<23>()
  {
    set_test_name("TaylorDeriv-> : ComputeYGradientFor() ");

    double DeltaX, DeltaY, Result;
    int RO;			// Reconstruction Order
    char number[4], msg[20];

    // Input Data
    DeltaX = 1.2333333;
    DeltaY = 0.1235;

    for (RO = 4; RO>=1; --RO){
      TD.GenerateContainer(RO);

      // Initialize container
      InitContainer(TD);
      
      // Get analytic result
      Result = 0.0;
      for (int i = 0; i<=TD.LastElem(); ++i){
	Result += ( Initialization(i) * TD(i).P2()*pow(DeltaX,TD(i).P1()) *
		    pow(DeltaY,TD(i).P2() - 1.0) );
      }

      // check
      sprintf(number,"%d",RO);
      strcpy(msg, "Order ");
      strcat(msg, number);
      ensure_distance(msg, TD.ComputeYGradientFor(DeltaX, DeltaY), Result, tol*(1.0 + Result));
    }

  }

}



// Test suite constructor
tut::TaylorDerivatives2D_TestGroup TaylorDerivatives2D_Test("Template Class:TaylorDerivatives 2D && double");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

