/*!\file UnitTesting.h
  \brief Header file defining unit testing framework. */

#ifndef _UNITTESTING_INCLUDED
#define _UNITTESTING_INCLUDED

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

/* Include MPI header file. */

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

/* Include TestRun header file for running available tests using TUT. */

#ifndef _NO_TUT_TESTING
#ifndef _TEST_RUN_INCLUDED
#include "TestRun.h"
#endif // _TEST_RUN_INCLUDED

#include "TestData.h"

#endif

/********************************************************
 * Subroutine for performing TUT unit testing.          *
 ********************************************************/
inline int Perform_UnitTesting(string TestSuite, int TestNumber, string Path_to_Root_Location) {

  int error_flag = 0;

#ifndef _NO_TUT_TESTING
  if (CFFC_Primary_MPI_Processor() ) {

    int Result;
    // change directory
    Result = chdir(Path_to_Root_Location.c_str());

    Print_(Result);
    Print_(getenv("PWD"))

    // set TestData::Root_Path
    if (Path_to_Root_Location.empty()){
      // check if the environment variable CFFC_UnitTest_Path is set
      

    } else {
      TestData::Root_Path = Path_to_Root_Location.c_str();
    }

    Test_Run(TestSuite, TestNumber);
    error_flag = 0;
  } /* endif */
  CFFC_Broadcast_MPI(&error_flag, 1);
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
