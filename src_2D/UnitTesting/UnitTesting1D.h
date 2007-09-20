/*!\file UnitTesting1D.h
  \brief Header file defining unit testing framework for 1D solvers. */

#ifndef _UNITTESTING_1D_INCLUDED
#define _UNITTESTING_1D_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>

using namespace std;

/* Include CFFC CFD header file. */

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

/* Include TestRun header file for running available tests using TUT. */

#ifndef _NO_TUT_TESTING
#ifndef _TEST_RUN_INCLUDED
#include "TestRun.h"
#endif // _TEST_RUN_INCLUDED
#endif

/********************************************************
 * Subroutine for performing TUT unit testing.          *
 ********************************************************/
inline int Perform_UnitTesting(string TestSuite, int TestNumber) {

  int error_flag = 0;

#ifndef _NO_TUT_TESTING
  Test_Run(TestSuite, TestNumber);
  error_flag = 0;
#else
  cout << "\n Current program has not been compiled with TUT unit testing on."
       << "\n To include the TUT unit testing framework, re-compile code with option: "
       << "\n TUT_TESTING = ON.\n";
  cout.flush();
  error_flag = 1;   
#endif

  return error_flag;

}

#endif /* _UNITTESTING_INCLUDED */
