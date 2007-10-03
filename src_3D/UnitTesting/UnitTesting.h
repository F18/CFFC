/* UnitTesting.h:  Header file defining unit testing framework. */

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
#include "../System/System_Linux.h"
#endif

/********************************************************
 * Subroutine for performing TUT unit testing.          *
 ********************************************************/
inline int Perform_UnitTesting(string TestSuite, int TestNumber, string Path_to_Root_Location) {

  int error_flag = 0;

#ifndef _NO_TUT_TESTING
  if (CFFC_Primary_MPI_Processor() ) {

    char * CFFC_UnitTest_Path(NULL);
    char * CurrentDirectory = new char [_MAX_PATH_];
    char * NewCurrentDirectory = new char [_MAX_PATH_];
    int RunTests = ON;
    

    try {
    
      // Get the current directory
      System::Get_Current_Path(CurrentDirectory);

      /* Change the current working directory to 'Path_to_Root_Location' or 'CFFC_UnitTest_Path_3D'.
	 The UnitTesting framework can only find the input and
	 output files if it is always launched from "src_3D" directory.
	 Therefore, if the executable is launched from a different directory, 
	 the 'Path_to_Root_Location' or 'CFFC_UnitTest_Path_3D' variables are 
	 used to determined the location of the "src_3D" directory.
	 Consequently, the current working directory is changed to 'Path_to_Root_Location' 
	 (path provided on the command line) if that exists, or to 'CFFC_UnitTest_Path_3D'
	 (path provided as an environment variable) if that is set up.
	 If none of these variables are set, it is assumed that the current working directory is "src_3D".
      */

      // check if the environment variable CFFC_UnitTest_Path_3D is set
      CFFC_UnitTest_Path = getenv("CFFC_UnitTest_Path_3D");

      if (!Path_to_Root_Location.empty()){
	// change the current working directory
	if ( System::Change_Current_Working_Directory(Path_to_Root_Location) == 0){

	  // determine the new current working directory
	  System::Get_Current_Path(NewCurrentDirectory);
	
	  if (strcmp(CurrentDirectory, NewCurrentDirectory) != 0){
	    // warn the user about the directory change
	    cout << "\nUsing directory " << Path_to_Root_Location << " to run the UnitTesting framework.\n";
	  }
	} else {
	  cerr << "\nUnitTesting ERROR: The path provided on the command line with option '-path' "
	       << "\ndoes not represent a valid path to 'src_3D' directory!"
	       << "\nPlease provide a valid path to run the UnitTesting framework.\n";
	
	  error_flag = 1;
	  RunTests = OFF;
	}

      } else if (CFFC_UnitTest_Path != NULL){
	// if CFFC_UnitTest_Path is set, changed working directory to that path

	// change the current working directory
	if ( System::Change_Current_Working_Directory(string(CFFC_UnitTest_Path)) == 0){
	  
	  // determine the new current working directory
	  System::Get_Current_Path(NewCurrentDirectory);

	  if (strcmp(CurrentDirectory, NewCurrentDirectory) != 0){
	    // warn the user about changing directory
	    cout << "\nUsing directory " << CFFC_UnitTest_Path << " to run the UnitTesting framework.\n";
	  }
	} else {
	  cerr << "\nUnitTesting ERROR: The environment variable 'CFFC_UnitTest_Path_3D' "
	       << "\ndoes not represent a valid path to 'src_3D' directory!"
	       << "\nPlease provide a valid path to run the UnitTesting framework.\n";
	
	  error_flag = 1;
	  RunTests = OFF;
	}

      } else {
	// check if the current directory ends in '/src_3D' (SEVEN characters before the end of string)
	char * buffer(NULL);
	if (strlen(CurrentDirectory) > 7){
	  buffer = & CurrentDirectory[strlen(CurrentDirectory) - 7];
	} 
	
	if (strcmp(buffer,"/src_3D") != 0){
	  cerr << "\nUnitTesting ERROR: The current directory is not '/src_3D'."
	       << "\nPlease run the executable from '/src_3D' directory or provide a valid"
	       << "\npath to this directory in order to run the UnitTesting framework."
	       << "\nUse: -path [path_to_src_3D] or set 'CFFC_UnitTest_Path_3D' environment variable.\n";
	  
	  error_flag = 1;
	  RunTests = OFF;
	}
      }

      // Deallocate memory
      delete [] CurrentDirectory; CurrentDirectory = NULL;
      delete [] NewCurrentDirectory; NewCurrentDirectory = NULL;

      // Run tests
      if (RunTests){
	Test_Run(TestSuite, TestNumber);
	error_flag = 0;
      }
    }
    catch (std::runtime_error &){
      // Deallocate memory
      delete [] CurrentDirectory; CurrentDirectory = NULL;
      delete [] NewCurrentDirectory; NewCurrentDirectory = NULL;

      // MPI calls
      error_flag = 1;
      CFFC_Broadcast_MPI(&error_flag, 1);
      CFFC_Finalize_MPI();

      // Re-throw the exception
      throw;
    }

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
