/*!\file test_Perfect_Gas_Shocks.cc
  \brief Regression tests for perfect-gas normal-shock relation tests. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../Perfect_Gas_Shocks.h"

namespace tut
{

  /* Define the test-specific data class and add data members
     when tests have complex or repeating creation phase. */
  class Data_Perfect_Gas_Shocks : public TestData {

    // Local variables
  public:

    // Constructor
    Data_Perfect_Gas_Shocks(){ }

  private:

  };

  /**
   * This group of declarations is just to register
   * test group in test-application-wide singleton.
   * Name of test group object (e.g. 'NAME'_TestSuite) shall
   * be unique in tut:: namespace. Alternatively, you
   * you may put it into anonymous namespace.
   */
  typedef test_group<Data_Perfect_Gas_Shocks> Perfect_Gas_Shocks_TestSuite;
  typedef Perfect_Gas_Shocks_TestSuite::object Perfect_Gas_Shocks_object;


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
  void Perfect_Gas_Shocks_object::test<1>()
  {
    set_test_name("Density-Ratio Calculation");

    double Mach(1.2345), Gamma(1.4), ratio, expected;

    ratio = normal_shock_density_ratio(Mach, Gamma);
    expected = 1.40158724179577;

    ensure_distance("ratio == expected", ratio, expected, expected*tol);
  }

  /* Test 2:*/
  template<>
  template<>
  void Perfect_Gas_Shocks_object::test<2>()
  {
    set_test_name("Veclocity-Ratio Calculation");

    double Mach(2.22345), Gamma(1.666666666), ratio, expected;

    ratio = normal_shock_velocity_ratio(Mach, Gamma);
    expected = 0.4017073169554191;

    ensure_distance("ratio == expected", ratio, expected, expected*tol);
  }

  /* Test 3:*/
  template<>
  template<>
  void Perfect_Gas_Shocks_object::test<3>()
  {
    set_test_name("Pressure-Ratio Calculation");

    double Mach(1.5563), Gamma(3.0), ratio, expected;

    ratio = normal_shock_pressure_ratio(Mach, Gamma);
    expected = 3.133104535;

    ensure_distance("ratio == expected", ratio, expected, expected*tol);
  }

  /* Test 4:*/
  template<>
  template<>
  void Perfect_Gas_Shocks_object::test<4>()
  {
    set_test_name("Temperature-Ratio Calculation");

    double Mach(3.2256), Gamma(1.4), ratio, expected;

    ratio = normal_shock_temperature_ratio(Mach, Gamma);
    expected = 2.954191833613753;

    ensure_distance("ratio == expected", ratio, expected, expected*tol);
  }

  //end tests
}



// Test suite constructor
tut::Perfect_Gas_Shocks_TestSuite Perfect_Gas_Shocks_TestSuite("Functions:Perfect_Gas_Shocks");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

