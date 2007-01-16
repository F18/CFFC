// Test_Run function: Runs the tests and prints the results
#include <iostream>
#include "tut.h"
#include "tut_reporter.h"

using namespace tut;
namespace tut
{
  test_runner_singleton runner;
}

void Test_Run( std::string TestSuit, int TestNumber) {

  reporter visi;
  runner.get().set_callback(&visi);

  try{
    if (TestSuit.empty())
    // run all tests in all groups
      runner.get().run_tests();
    else if (TestNumber>0)
      runner.get().run_test(TestSuit, TestNumber); 
    else
      runner.get().run_tests(TestSuit); 
  }
  catch( const std::exception& ex){
    std::cerr << "tut raised ex: " << ex.what() << std::endl;
  }
};
