/*!\file test_AdvectDiffuse2DState.cc
  \brief Regression tests for AdvectDiffuse2D_State class. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../New_AdvectDiffuse2DState.h"

namespace tut
{

  /* Define the test-specific data class and add data members 
     when tests have complex or repeating creation phase. */
  class Data_AdvectDiffuse2DState : public TestData {

    // Local variables
  public:

    // Constructor
    Data_AdvectDiffuse2DState(){
      /* Set the global path for this test suite. 
	 It's a good practice to put it in the constructor of the data class in order to be set
	 automatically for each individual test. Declare it relative to the /src_2D directory,
	 otherwise the framework might not find the input and output files. */
      
      set_test_suite_path("AdvectDiffuse2D/UnitTests/");
    }

  private:
    
  };

  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. 'NAME'_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_AdvectDiffuse2DState> AdvectDiffuse2DState_TestSuite;
  typedef AdvectDiffuse2DState_TestSuite::object AdvectDiffuse2DState_object;


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
  void AdvectDiffuse2DState_object::test<1>()
  {

    set_test_name("Uniform velocity field, 0 degrees");

    VelocityFields::Set_UniformFlow_MagnitudeAngleVelocity(10.0,0);

    // === check
    ensure_equals("XAxis", VelocityFields::Uniform_Flow_XAxis(1.0,2.0), Vector2D(10.0,0.0));
    ensure_equals("YAxis", VelocityFields::Uniform_Flow_YAxis(1.0,2.0), Vector2D(0.0,0.0));
    ensure_equals("Uniform flow", VelocityFields::Uniform_Flow(1.0,2.0), Vector2D(10.0,0.0));
  }

  /* Test 2:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<2>()
  {

    set_test_name("Uniform velocity field, 30 degrees");

    VelocityFields::Set_UniformFlow_MagnitudeAngleVelocity(10.0,30);

    Vector2D tol(1.0e-13);

    // === check
    ensure_distance("XAxis", VelocityFields::Uniform_Flow_XAxis(1.0,2.0), Vector2D(8.660254037844386467637e+0,0.0), tol);
    ensure_distance("YAxis", VelocityFields::Uniform_Flow_YAxis(1.0,2.0), Vector2D(0.0,5.0), tol);
    ensure_distance("Uniform flow", VelocityFields::Uniform_Flow(1.0,2.0), Vector2D(8.660254037844386467637e+0,5.0), tol);
  }

  /* Test 3:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<3>()
  {

    set_test_name("Uniform velocity field, 330 degrees");

    VelocityFields::Set_UniformFlow_MagnitudeAngleVelocity(10.0,330);

    Vector2D tol(1.0e-13);

    // === check
    ensure_distance("XAxis", VelocityFields::Uniform_Flow_XAxis(1.0,2.0), Vector2D(8.660254037844386467637e+0,0.0), tol);
    ensure_distance("YAxis", VelocityFields::Uniform_Flow_YAxis(1.0,2.0), Vector2D(0.0,-5.0), tol);
    ensure_distance("Uniform flow", VelocityFields::Uniform_Flow(1.0,2.0), Vector2D(8.660254037844386467637e+0,-5.0), tol);
  }

  /* Test 4:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<4>()
  {

    set_test_name("Uniform velocity field, 210 degrees");

    VelocityFields::Set_UniformFlow_MagnitudeAngleVelocity(10.0,210);

    Vector2D tol(1.0e-13);

    // === check
    ensure_distance("XAxis", VelocityFields::Uniform_Flow_XAxis(1.0,2.0), Vector2D(-8.660254037844386467637e+0,0.0), tol);
    ensure_distance("YAxis", VelocityFields::Uniform_Flow_YAxis(1.0,2.0), Vector2D(0.0,-5.0), tol);
    ensure_distance("Uniform flow", VelocityFields::Uniform_Flow(1.0,2.0), Vector2D(-8.660254037844386467637e+0,-5.0), tol);
  }


  /* Test 5:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<5>()
  {

    set_test_name("Rotational velocity field about origin");

    VelocityFields::Set_UniformRotationalFlow_Parameters(10.0);
    
    // === check
    ensure_equals("Velocity", VelocityFields::Rotational_Flow_About_Origin(3.0,4.0), Vector2D(-40,30));
  }

  /* Test 6:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<6>()
  {

    set_test_name("Rotational velocity field about arbitrary point");

    VelocityFields::Set_UniformRotationalFlow_Parameters(10.0, Vector2D(2.5,4.5));

    // === check
    ensure_equals("Velocity", VelocityFields::Rotational_Flow_About_Arbitrary_Point(5.5,8.5), Vector2D(-40,30));
  }

  /* Test 7:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<7>()
  {

    set_test_name("Rotational velocity field about origin, inverse variation angular velocity");

    VelocityFields::Set_RotationalFlow_Parameters(10.0, "Inverse_Proportional_Distance");

    // === check
    ensure_equals("Velocity First Point", VelocityFields::Rotational_Flow_About_Arbitrary_Point(3,4), Vector2D(-8,6));
    ensure_equals("Velocity Second Point", VelocityFields::Rotational_Flow_About_Arbitrary_Point(6,8), Vector2D(-8,6));
  }

  /* Test 8:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<8>()
  {

    set_test_name("Rotational velocity field about arbitrary point, inverse variation angular velocity");

    VelocityFields::Set_RotationalFlow_Parameters(10.0, "Inverse_Proportional_Distance", Vector2D(2.5,4.5));

    // === check
    ensure_equals("Velocity First Point", VelocityFields::Rotational_Flow_About_Arbitrary_Point(5.5,8.5), Vector2D(-8,6));
    ensure_equals("Velocity Second Point", VelocityFields::Rotational_Flow_About_Arbitrary_Point(8.5,12.5), Vector2D(-8,6));
  }

  /* Test 9:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<9>()
  {

    set_test_name("Connect_Pointer_To_Flow_Field()");

    AdvectDiffuse2D_State_New State;

    VelocityFields::Set_RotationalFlow_Parameters(10.0, "Inverse_Proportional_Distance", Vector2D(2.5,4.5));

    VelocityFields::Connect_Pointer_To_Flow_Field(State.V);

    // === check
    ensure_equals("Velocity First Point", State.V(5.5,8.5), Vector2D(-8,6));
    ensure_equals("Velocity Second Point", State.V(8.5,12.5), Vector2D(-8,6));
  }

  /* Test 10:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<10>()
  {

    set_test_name("Connect_Pointer_To_Diffusion_Field()");

    AdvectDiffuse2D_State_New State;

    DiffusionFields::Set_ConstantDiffusionField(10.0);

    DiffusionFields::Connect_Pointer_To_Diffusion_Field(State.k);

    // === check
    ensure_equals("Diffusion coeff. @ location", State.k(5.5,8.5,0.0), 10.0 );
  }

  /* Test 11:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<11>()
  {

    set_test_name("operator ==");

    AdvectDiffuse2D_State_New State(12.2323), State_Copy(12.2323);

    // === check
    ensure_equals("equal states", State == State_Copy, true );
  }
  
  /* Test 12:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<12>()
  {

    set_test_name("operator ==");

    AdvectDiffuse2D_State_New State(12.2323), State_Copy(12.232);

    // === check
    ensure_equals("different states", State == State_Copy, false );
  }

  /* Test 13:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<13>()
  {

    set_test_name("operator !=");

    AdvectDiffuse2D_State_New State(12.2323), State_Copy(12.232);

    // === check
    ensure_equals("different states", State != State_Copy, true );
  }
  
  /* Test 14:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<14>()
  {

    set_test_name("operator !=");

    AdvectDiffuse2D_State_New State(12.2323), State_Copy(12.2323);

    // === check
    ensure_equals("equal states", State != State_Copy, false);
  }

  /* Test 15:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<15>()
  {

    set_test_name("Assignment operator");

    AdvectDiffuse2D_State_New State(12.2323), State_Copy;

    // === check
    ensure_equals("different states", State != State_Copy, true);

    // use operator =
    State_Copy = State;

    // === check operator
    ensure_equals("equal states", State == State_Copy, true);
  }

  /* Test 16:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<16>()
  {

    set_test_name("Copy constructor");

    AdvectDiffuse2D_State_New State(12.2323), State_Copy(State);

    // === check
    ensure_equals("equal states", State == State_Copy, true);
  }

  /* Test 17:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<17>()
  {

    set_test_name("operator +");

    AdvectDiffuse2D_State_New A(12.2323), B(A), C, Result(12.2323*2);

    // summation
    C = A + B;

    // === check
    ensure_equals("unmodified A", A, AdvectDiffuse2D_State_New(12.2323) );
    ensure_equals("unmodified B", B, AdvectDiffuse2D_State_New(12.2323) );
    ensure_equals("summation", C , Result);
  }

  /* Test 18:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<18>()
  {

    set_test_name("operator -");

    AdvectDiffuse2D_State_New A(12.2323), B(12.2325), C, Result(-0.0002), tol(1.0e-14);

    // subtraction
    C = A - B;

    // === check
    ensure_equals("unmodified A", A, AdvectDiffuse2D_State_New(12.2323) );
    ensure_equals("unmodified B", B, AdvectDiffuse2D_State_New(12.2325) );
    ensure_distance("subtraction", C , Result, tol);
  }

  /* Test 19:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<19>()
  {

    set_test_name("operator *");

    AdvectDiffuse2D_State_New A(12.2323), B(12.2325), C, Result(12.2323*12.2325);

    // multiplication
    C = A * B;

    // === check
    ensure_equals("unmodified A", A, AdvectDiffuse2D_State_New(12.2323) );
    ensure_equals("unmodified B", B, AdvectDiffuse2D_State_New(12.2325) );
    ensure_equals("multiplication", C , Result);
  }

  /* Test 20:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<20>()
  {

    set_test_name("operator * ");

    double a(-0.2323);
    AdvectDiffuse2D_State_New A(12.2323), B, C, Result(12.2323*a);

    // multiplication
    B = A * a;
    C = a * A;

    // === check
    ensure_equals("unmodified A", A, AdvectDiffuse2D_State_New(12.2323) );
    ensure_equals("scalar multiplication I", B , Result);
    ensure_equals("scalar multiplication II", C , Result);
  }

  /* Test 21:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<21>()
  {

    set_test_name("operator / ");

    double a(-0.2323);
    AdvectDiffuse2D_State_New A(12.2323), B, C, Result(12.2323/a);

    // division
    B = A / a;

    // === check
    ensure_equals("unmodified A", A, AdvectDiffuse2D_State_New(12.2323) );
    ensure_equals("division with scalar", B , Result);
  }

  /* Test 22:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<22>()
  {

    set_test_name("unary operator + ");

    AdvectDiffuse2D_State_New A(12.2323), B, Result(12.2323);

    // unary operator
    B = + A;

    // === check
    ensure_equals("unmodified A", A, AdvectDiffuse2D_State_New(12.2323) );
    ensure_equals("result", B , Result);
  }

  /* Test 23:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<23>()
  {

    set_test_name("unary operator - ");

    AdvectDiffuse2D_State_New A(12.2323), B, Result(-12.2323);

    // unary operator
    B = - A;

    // === check
    ensure_equals("unmodified A", A, AdvectDiffuse2D_State_New(12.2323) );
    ensure_equals("result", B , Result);
  }

  /* Test 24:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<24>()
  {

    set_test_name("operator +=");

    AdvectDiffuse2D_State_New A(12.2323), B(A), Result(12.2323*2);

    // summation
    A += B;

    // === check
    ensure_equals("unmodified B", B, AdvectDiffuse2D_State_New(12.2323) );
    ensure_equals("summation", A , Result);
  }

  /* Test 25:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<25>()
  {

    set_test_name("operator -=");

    AdvectDiffuse2D_State_New A(12.2323), B(12.2325),  Result(-0.0002), tol(1.0e-13);

    // subtraction
    A -= B;

    // === check
    ensure_equals("unmodified B", B, AdvectDiffuse2D_State_New(12.2325) );
    ensure_distance("subtraction", A , Result, tol);
  }

  /* Test 26:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<26>()
  {

    set_test_name("operator *=");

    AdvectDiffuse2D_State_New A(12.2323), B(12.2325), Result(12.2323*12.2325);

    // multiplication
    A *= B;

    // === check
    ensure_equals("unmodified B", B, AdvectDiffuse2D_State_New(12.2325) );
    ensure_equals("multiplication", A , Result);
  }

  /* Test 27:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<27>()
  {

    set_test_name("operator * ");

    double a(-0.2323);
    AdvectDiffuse2D_State_New A(12.2323), Result(12.2323*a);

    // multiplication
    A *= a;

    // === check
    ensure_equals("scalar multiplication", A , Result);
  }

  /* Test 28:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<28>()
  {

    set_test_name("operator / ");

    double a(-0.2323);
    AdvectDiffuse2D_State_New A(12.2323), Result(12.2323/a);

    // division
    A /= a;

    // === check
    ensure_equals("division with scalar", A , Result);
  }

  /* Test 29:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<29>()
  {

    set_test_name("combination of operators ");

    double a(-0.2323);
    AdvectDiffuse2D_State_New A(12.2323), B(-0.2324), C(-34.24), Result(12.2323 - (12.2323 + (-0.2324))*(-34.24)/a );

    // division
    A += -(A + B)*C/a;

    // === check
    ensure_equals("division with scalar", A , Result);
  }

  /* Test 30:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<30>()
  {

    set_test_name("input-output operators ");

    AdvectDiffuse2D_State_New A(12.2323);

    // === check
    Check_Input_Output_Operator(A);
  }

  /* Test 31:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<31>()
  {

    set_test_name("Advective flux");

    VelocityFields::Set_RotationalFlow_Parameters(10.0, "Inverse_Proportional_Distance", Vector2D(2.5,4.5));
    VelocityFields::Connect_Pointer_To_Flow_Field(AdvectDiffuse2D_State_New::V);

    AdvectDiffuse2D_State_New A(12.2323);
    Vector2D Result = Vector2D(-8,6)*12.2323, tol(1.0e-13);

    // === check advective flux
    ensure_distance("Fa at Point", A.Fa(8.5,12.5), Result, tol);
    ensure_distance("Second Fa call", Fa(A,8.5,12.5), Result, tol);
  }

  /* Test 32:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<32>()
  {

    set_test_name("Diffusive flux");

    // Set velocity field
    VelocityFields::Set_RotationalFlow_Parameters(10.0, "Inverse_Proportional_Distance", Vector2D(2.5,4.5));
    VelocityFields::Connect_Pointer_To_Flow_Field(AdvectDiffuse2D_State_New::V);

    // Set diffusion field
    DiffusionFields::Set_ConstantDiffusionField(10.0);
    DiffusionFields::Connect_Pointer_To_Diffusion_Field(AdvectDiffuse2D_State_New::k);

    // Set initial data
    AdvectDiffuse2D_State_New A(12.2323);
    double dUdx(1.5), dUdy(3.4);
    Vector2D tol(1.0e-13), GradU(dUdx,dUdy), Result(-10.0*dUdx,-10*dUdy);

    // === check diffusive flux
    ensure_distance("Fd at Point", A.Fd(dUdx,dUdy,8.5,12.5), Result, tol);
    ensure_distance("Fd call 2"  , A.Fd(GradU,8.5,12.5), Result, tol);
    ensure_distance("Fd call 3"  , Fd(A,dUdx,dUdy,8.5,12.5), Result, tol);
    ensure_distance("Fd call 4"  , Fd(A,GradU,8.5,12.5), Result, tol);
  }

  /* Test 33:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<33>()
  {

    set_test_name("Upwind advective flux");

    // Set velocity field
    VelocityFields::Set_RotationalFlow_Parameters(10.0, "Inverse_Proportional_Distance", Vector2D(2.5,4.5));
    VelocityFields::Connect_Pointer_To_Flow_Field(AdvectDiffuse2D_State_New::V);

    // Set initial data
    AdvectDiffuse2D_State_New A(12.0), B(5.0);
    Vector2D Normal(1,0);
    
    // === check advective flux
    ensure_distance("Fa I", Fa(A,B,8.5,12.5,Normal), -40.0, tol);
  }

  /* Test 34:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<34>()
  {

    set_test_name("Upwind advective flux");

    // Set velocity field
    VelocityFields::Set_RotationalFlow_Parameters(10.0, "Inverse_Proportional_Distance", Vector2D(2.5,4.5));
    VelocityFields::Connect_Pointer_To_Flow_Field(AdvectDiffuse2D_State_New::V);

    // Set initial data
    AdvectDiffuse2D_State_New A(12.0), B(5.0);
    Vector2D Normal(0,1);
    
    // === check advective flux
    ensure_distance("Fa I", Fa(A,B,8.5,12.5,Normal), 72.0, tol);
  }

  /* Test 35:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<35>()
  {

    set_test_name("Diffusive flux");

    // Set velocity field
    VelocityFields::Set_RotationalFlow_Parameters(10.0, "Inverse_Proportional_Distance", Vector2D(2.5,4.5));
    VelocityFields::Connect_Pointer_To_Flow_Field(AdvectDiffuse2D_State_New::V);

    // Set diffusion field
    DiffusionFields::Set_ConstantDiffusionField(10.0);
    DiffusionFields::Connect_Pointer_To_Diffusion_Field(AdvectDiffuse2D_State_New::k);

    // Set initial data
    AdvectDiffuse2D_State_New Ul(12.2323), Ur(5.34);
    double dUdx_L(1.5), dUdy_L(3.4), dUdx_R(1.5), dUdy_R(3.4);
    Vector2D Normal(1,0);
    double Result(-10*1.5);

    // === check diffusive flux
    ensure_distance("Fd at Point", Fd(Ul,dUdx_L,dUdy_L,Ur,dUdx_R,dUdy_R,8.5,12.5,Normal), Result, tol);
  }

  /* Test 36:*/
  template<>
  template<>
  void AdvectDiffuse2DState_object::test<36>()
  {

    set_test_name("Diffusive flux");

    // Set velocity field
    VelocityFields::Set_RotationalFlow_Parameters(10.0, "Inverse_Proportional_Distance", Vector2D(2.5,4.5));
    VelocityFields::Connect_Pointer_To_Flow_Field(AdvectDiffuse2D_State_New::V);

    // Set diffusion field
    DiffusionFields::Set_ConstantDiffusionField(10.0);
    DiffusionFields::Connect_Pointer_To_Diffusion_Field(AdvectDiffuse2D_State_New::k);

    // Set initial data
    AdvectDiffuse2D_State_New Ul(12.2323), Ur(5.34);
    double dUdx_L(1.5), dUdy_L(6.4), dUdx_R(3.5), dUdy_R(3.4);
    Vector2D Normal(0,1);
    double Result(-10*0.5*(6.4+3.4));

    // === check diffusive flux
    ensure_distance("Fd at Point", Fd(Ul,dUdx_L,dUdy_L,Ur,dUdx_R,dUdy_R,8.5,12.5,Normal), Result, tol);
  }

}



// Test suite constructor
tut::AdvectDiffuse2DState_TestSuite AdvectDiffuse2DStateTestSuite("Class:AdvectDiffuse2D_State");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

