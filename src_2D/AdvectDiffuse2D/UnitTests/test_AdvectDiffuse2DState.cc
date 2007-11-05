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

