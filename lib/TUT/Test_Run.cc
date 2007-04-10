// Test_Run function: Runs the tests and prints the results
#include "Test_Run.h"
#include "tut_restartable.h"

using namespace tut;
using std::cout;
using std::endl;
using std::cerr;

namespace tut
{
  test_runner_singleton runner;
}

void Test_Run( std::string TestSuit, int TestNumber) {

  try{
    reporter visi;
    restartable_wrapper restartable;

    if (TestSuit.empty()){
      cout << "\n Regression Test: run all registered tests.\n"
	   << "*********************************************\n";

      // run all tests in all groups
      runner.get().set_callback(&visi);
      runner.get().run_tests();

    } else if (TestSuit == "regression") {
      cout << "\n Code Regression Checking: run all registered tests.\n"
	   << " Test results are printed only after all the tests finished.\n"
	   << " NB: If the application is terminated by OS, just restart it!\n"
	   << "*************************************************************\n";

      restartable.set_callback(&visi);
      restartable.run_tests();

    } else if (TestSuit == "list"){
      runner.get().set_callback(&visi);
      cout << "\n Registered test groups:" << endl;
      groupnames gl = tut::runner.get().list_groups();
      groupnames::const_iterator i = gl.begin();
      groupnames::const_iterator e = gl.end();
      while(i != e)
	{
	  cout << "  " << *i << endl;
	  ++i;
	}
      cout << endl;

    }else if (TestNumber>0){
      cout << "\n Run test: " << TestSuit << ", #" << TestNumber << endl
	   << "*****************************************************\n";
      runner.get().set_callback(&visi);
      runner.get().run_test(TestSuit, TestNumber); 

    } else {
      cout << "\n Run test suit: " << TestSuit << endl
	   << "*****************************************************\n";
      runner.get().set_callback(&visi);
      runner.get().run_tests(TestSuit); 
    } // endif

  }
  catch( const std::exception& ex){
    cerr << "tut raised ex: " << ex.what() << endl;
  }

};
