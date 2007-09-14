/* test_BasicTestData.cc:  Regression tests for the BasicTestData class. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFDkit+caboodle header files */
#include "BasicTestData.h"


namespace tut
{
  typedef test_group<TestData> TestGroup;
  typedef TestGroup::object object;
}


namespace tut
{

  /*******************************************************
   * Tests for CLASS : BasicTestData                     *
   ******************************************************/

  /*Test1: */
  template<>
  template<>
  void object::test<1>()
  {
    set_test_name("Constructor");
    TestData VoidConstruct;

    #ifdef _GNU_GCC_V3
    ensure_equals("digits", digits, 15); // not quite sure
    #else
    ensure_equals("digits", digits, 14);
    #endif

    ensure_equals("output_file_name",strcmp(output_file_name, Global_Output_Path) == 0, true);
    ensure_equals("input_file_name",strcmp(input_file_name, Global_Input_Path) == 0, true);
  }

  /*Test2: */
  template<>
  template<>
  void object::test<2>()
  {
    set_test_name("open_output_stream()");

    try {
      ofstream output_stream;
      char * output_stream_file_name = "UnitTesting/Results_Tests/test_BasicTestData/BasicTestData_Test";
      open_output_stream(output_stream,output_stream_file_name);

      output_stream.close();
    } catch(runtime_error &){
      // Failure
      throw tut::failure("open_output_stream()");
    }
  }

  /*Test3: */
  template<>
  template<>
  void object::test<3>()
  {
    set_test_name("open_output_stream()");

    try {
      ofstream output_stream;
      char * output_stream_file_name = "UnitTesting/Results_Tests/BasicTestData_Test/FailureTest";
      open_output_stream(output_stream,output_stream_file_name);

      output_stream.close();
    } catch(runtime_error &){
      // OK
    }
  }

  /*Test4: */
  template<>
  template<>
  void object::test<4>()
  {
    set_test_name("open_input_stream()");

    try {
      ifstream input_stream;
      char * input_stream_file_name = "UnitTesting/Results_Tests/test_BasicTestData/BasicTestData_Test";
      open_input_stream(input_stream,input_stream_file_name);

      input_stream.close();
    } catch(runtime_error &){
      // Failure
      throw tut::failure("open_input_stream()");
    }
  }

  /*Test5: */
  template<>
  template<>
  void object::test<5>()
  {
    set_test_name("open_input_stream()");

    try {
      ifstream input_stream;
      char * input_stream_file_name = "UnitTesting/Results_Tests/BasicTestData_Test/FailureTest";
      open_input_stream(input_stream,input_stream_file_name);

      input_stream.close();
    } catch(runtime_error &){
      // OK
    }
  }

  /*Test6: */
  template<>
  template<>
  void object::test<6>()
  {
    set_test_name("stream.is_open()");

    try {
      ifstream input_stream;

      if(input_stream.is_open()){
	throw tut::failure("input_stream");
      } else {
	throw runtime_error("expected");
      }

    } catch(runtime_error &){
      // OK
    }
  }

  /*Test7: */
  template<>
  template<>
  void object::test<7>()
  {
    set_test_name("stream.is_open()");

    try {
      ofstream output_stream;

      if(!output_stream.is_open()){
	throw runtime_error("expected");
      } else {
	throw tut::failure("output_stream");
      }

    } catch(runtime_error &){
      // OK
    }
  }


  /*Test8: */
  template<>
  template<>
  void object::test<8>()
  {
    set_test_name("set_local_output_path()");

    set_local_output_path("test_BasicTestData");
    ensure_equals("Local_Output_Path",
		  strcmp(Local_Output_Path,"UnitTesting/Results_Tests/test_BasicTestData/")==0,
		  true);
  }

  /*Test9: */
  template<>
  template<>
  void object::test<9>()
  {
    set_test_name("set_local_input_path()");

    try{
      set_local_output_path("__test_BasicTestData__/Test/");
      throw tut::failure("Local_Input_Path");

    } catch(tut::bad_ctor &){
      // OK
    }
  }

  /*Test10: */
  template<>
  template<>
  void object::test<10>()
  {
    set_test_name("in()/out()");
    set_local_output_path("test_BasicTestData");

    const char *File = "OutputInput_Test";

    Open_Output_File(File);
    Open_Output_File_For_Reading();

    double a(1231313.131323434), a_read(0.0);

    out() << setprecision(16) << a << endl;
    in()  >> a_read;

    ensure_equals("Same double", a_read, a);
  }

  /*Test11: */
  template<>
  template<>
  void object::test<11>()
  {
    set_test_name("Check_Input_Output_Operator()");
    __LocalTest__ MyData;

    MyData.x = 10; MyData.y = 1.3433e-29; MyData.z = 32.03232e5;
    Check_Input_Output_Operator("first check", MyData);

    MyData.x = 10.12312312312313; MyData.y = 1.3433e-29; MyData.z = 32.03232e5;
    Check_Input_Output_Operator("second check", MyData);

  }

  /*Test12: */
  template<>
  template<>
  void object::test<12>()
  {
    CompareFilesWithNumericalTolerance("Test",
				       "UnitTesting/Results_Tests/test_BasicTestData/Compare.dat",
				       "UnitTesting/Results_Tests/test_BasicTestData/Master.dat",
				       1.0e-15, 1.0e-15);
  }
}

tut::TestGroup BasicTestData("selftest");
