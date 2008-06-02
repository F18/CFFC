/*!\file test_Levermore1D_UniformMesh.cc
  \brief Regression tests for template class Levermore1D_UniformMesh datatype. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "TestData.h"
#include "../Levermore1D.h"
#include "test_Levermore1D_defines.h"

namespace tut
{

  /* Define the test-specific data class and add data members
     when tests have complex or repeating creation phase. */
  class Data_Levermore1D_UniformMesh : public TestData {

    // Local variables
  public:

    // Constructor
    Data_Levermore1D_UniformMesh(){
      if(!Levermore1D_Vector::length_is_set())
	Levermore1D_Vector::set_length(LEVERMORE1D_VECTOR_LENGTH);
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
  typedef test_group<Data_Levermore1D_UniformMesh> Levermore1D_UniformMesh_TestSuite;
  typedef Levermore1D_UniformMesh_TestSuite::object Levermore1D_UniformMesh_object;


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
  void Levermore1D_UniformMesh_object::test<1>()
  {
    set_test_name("dUdA_inverse construction");

    double 
      rho(1.225),
      u(123.0),
      p(101325.0);

    Levermore1D_pState W1(rho,u,p);
    Levermore1D_cState U1(W1);
    Levermore1D_weights A1(U1);

    Cell1D_Uniform dummy;

    Levermore1D_UniformMesh L1(W1,dummy);

    DenseMatrix identity;

    L1.calculate_Hessians();

    identity = L1.dUdA * L1.dUdA_inv;

    for(int i=0; i < Levermore1D_Vector::get_length(); ++i) {
      for(int j=0; j < Levermore1D_Vector::get_length(); ++j) {
	if(i==j) {
	  ensure_distance("Diagonal == 1",identity(i,j),1.0, tol*U1.moment( 2*((i+j)/2), A1, u)*100.0 );
	} else {
	  ensure_distance("Off-diagonal == 0",identity(i,j),0.0, tol*U1.moment( 2*((i+j)/2), A1, u)*100.0 );
	}
      }
    }
  } //end test 1

  /* Test 2:*/
  template<>
  template<>
  void Levermore1D_UniformMesh_object::test<2>()
  {
    set_test_name("Checking dUdA_inv * delta U");

//    double 
//      rho(0.556),
//      u(-101.0),
//      p(40123.0);
//
//    Levermore1D_pState W1(rho,u,p);
//    Cell1D_Uniform dummy;
//
//    Levermore1D_UniformMesh L1(W1,dummy);
//    L1.calculate_Hessians();
//
//    Levermore1D_Vector dU;
//
//    for(int i=1; i <= Levermore1D_Vector::get_length(); ++i) {
//      dU[i] = L1.U[i] * ((double)i)*0.01;
//    }
//
//    Levermore1D_UniformMesh  L2(L1);
//
//    L2.U += dU;
//    L2.A += L1.dUdA_inv*dU;
//
//    cout << endl << L1.U << endl << L2.U << endl << L1.A << endl << L2.A << endl << endl;
//
//    cout << L1.U.detector_value(L1.A) << endl;
//    cout << L2.U.detector_value(L2.A) << endl << endl;
//
//    Levermore1D_UniformMesh L3(L2);
//    L3.A.set_from_U(L3.U);
//
//    cout << L2.A << endl;
//    cout << L3.A << endl << endl;
//
//    double us(L2.U[2]/L2.U[1]);
//
//    cout << L2.U.moment(Levermore1D_Vector::get_length(),L2.A,us) << endl;
//    cout << L3.U.moment(Levermore1D_Vector::get_length(),L3.A,us) << endl << endl;


  }//end test 2

  //end tests
}



// Test suite constructor
tut::Levermore1D_UniformMesh_TestSuite Levermore1D_UniformMeshTestSuite("Class:Levermore1D_UniformMesh");

/*************************************************************
 Guidelines for naming "Test Suite Name".
 Write the name as "Category:Representative Name"

  e.g. "Class:NameOfTheClass"
       "Template Class:NameOfTheTemplateClass [&& type used for testing]"
       "Integration:NameOfTheMethod"
       "Linear Systems Solvers:NameOfTheMethod"

*/

