/*******************************************************************
 *******************************************************************
 **********                                               **********
 ******                DUSTY2D Version 1.00                   ******
 ****                                                           ****
 ****                                                           ****
 ******       UTIAS, CFD & Propulsion Group, 1999-2007        ******
 **********                                               **********
 *******************************************************************
 *******************************************************************/


// Include the header files defining various variable types and data
// structures, classes, functions, operators, global variables,
// etc..

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

// Include CFFC header files.

#include "Dusty2DQuad.h"
#include "../MPI/MPI.h"
#include "../ICEM/ICEMCFD.h"
#include "../UnitTesting/UnitTesting.h"

/* Begin Dusty2D program. */

int main(int num_arg, char *arg_ptr[]) {

  /********************************************************************
   * VARIABLE DECLARATIONS                                            *
   ********************************************************************/

  // Command line flags:  
  int version_flag,
      help_flag,
      batch_flag,
      pde_flag,
      file_flag,
      error_flag,
      mpirun_flag,
      test_flag;

  // Name of program entered on the command line:
  char *command_name_ptr;

  // Title of code:
  char *program_title_ptr = 
     "Dusty2D: Numerical Solvers for the Dusty2D PDEs.";

  // Version of code:
  char *program_version_ptr = 
     "Version 1.00, UTIAS CFD & Propulsion Group, 1999-2007.";
  
  // Input file name:
  char *Input_File_Name_ptr = "dusty2D.in";

  // Equation type indicator:
  char *Equation_Type_ptr = "Dusty2D";
  char Equation_Type[256];

  // Input file stream:
  ifstream Input_File;

  // Unit testing flags:
  string TestSuite, TestRootPath;
  int TestNumber=0;

  /********************************************************************  
   * LOCAL VARIABLE INITIALIZATIONS                                   *
   ********************************************************************/

  // Set default equation type.
  strcpy(Equation_Type,Equation_Type_ptr);

  /********************************************************************
   * PARSE COMMAND LINE ARGUMENTS                                     *
   ********************************************************************/

  // Initialize command line flags.
  version_flag = 0;
  help_flag = 0;
  batch_flag = 0;
  pde_flag = 0;
  file_flag = 0;
  error_flag = 0;
  mpirun_flag = 0;
  test_flag = 0;

  // Save the command line name of the program for future use.
  command_name_ptr = arg_ptr[0];

  // Parse and interpret command line arguments.  Note that there
  // are several different possible arguments which are:
  // 1) -v  lists program version information to standard output,
  // 2) -h  lists all possible optional arguments for program,
  // 3) -i  execute program in interactive mode (default mode),
  // 4) -b  execute program in batch mode,
  // 5) -pde type   sets type of partial differential equation
  //                to be solve to "type" (default is "Dusty2D"),
  // 6) -f name  uses "name" as the input data file rather than
  //             the standard input data file "dusty2D.in".

  if (num_arg >= 2) {
    for (int i = 1; i < num_arg; i++) {
      if (strcmp(arg_ptr[i],"-v") == 0 ||
          strcmp(arg_ptr[i],"--version") == 0) {
        version_flag = 1;
      } else if (strcmp(arg_ptr[i],"-h") == 0 ||
                 strcmp(arg_ptr[i],"--help") == 0) {
        help_flag = 1;
      } else if (strcmp(arg_ptr[i],"-i") == 0 ||
                 strcmp(arg_ptr[i],"--inter") == 0) {
        batch_flag = 0;
      } else if (strcmp(arg_ptr[i],"-b") == 0 ||
                 strcmp(arg_ptr[i],"--batch") == 0) {
        batch_flag = 1;
      } else if (strcmp(arg_ptr[i],"-pde") == 0) {
        pde_flag = 1;
	//} else if (strcmp(arg_ptr[i-1],"-pde") == 0) {
        //Equation_Type_ptr = arg_ptr[i];
        //strcpy(Equation_Type, Equation_Type_ptr);
      } else if (strcmp(arg_ptr[i],"-f") == 0) {
        file_flag = 1;
      } else if (strcmp(arg_ptr[i-1],"-f") == 0) {
        Input_File_Name_ptr = arg_ptr[i];
      } else if (strcmp(arg_ptr[i],"-p4pg") == 0) {
        mpirun_flag = 1;
      } else if (strcmp(arg_ptr[i-1],"-p4pg") == 0) {
        //mpirun_flag = 1;
      } else if (strcmp(arg_ptr[i],"-p4wd") == 0) {
        mpirun_flag = 1;
      } else if (strcmp(arg_ptr[i-1],"-p4wd") == 0) {
        //mpirun_flag = 1;
      } else if (strcmp(arg_ptr[i], "-t") == 0) {
	test_flag=1;
	if (num_arg-1>i && strcmp(arg_ptr[i+1], "-path") != 0 ){
	  // Read TestSuite name
	  TestSuite= arg_ptr[++i];
	} /* endif */
	if (num_arg-1>i && strcmp(arg_ptr[i+1], "-path") != 0){
	  // Read TestNumber
	  TestNumber = atoi(arg_ptr[++i]);
	} /* endif */
      } else if (strcmp(arg_ptr[i], "-path") == 0) {
	if (num_arg-1>i){
	  TestRootPath = arg_ptr[++i];
	}  /* endif */
      } else {
        error_flag = 1;
      }
    }
  }

  // Display command line argument error message and terminate the
  // program as required.
  if (error_flag) {
    cout << "\n Dusty2D ERROR: Invalid command line argument.\n";
    return error_flag;
  }

  /********************************************************************
   * INITIALIZE MPI                                                   *
   ********************************************************************/
  CFFC_Initialize_MPI(num_arg,arg_ptr);
  if (!CFFC_Primary_MPI_Processor()) batch_flag = 1;

  /********************************************************************
   * DISPLAY THE PROGRAM TITLE AND VERSION INFORMATION AS REQUIRED.   *
   ********************************************************************/
  if (CFFC_Primary_MPI_Processor() && (version_flag || help_flag || !batch_flag)) {
     cout << endl << program_title_ptr << endl;
     cout << program_version_ptr << endl;
     cout << "Built using " << CFFC_Version() << endl;
     cout << CFFC_Version_MPI() << endl;
     cout << ICEMCFD_Version() << "\n";
     cout << "Built using MV++, SparseLib++, IML++, and BPKIT Libraries\n";
  }
  if (version_flag) {
     CFFC_Finalize_MPI();
     return (0);
  } /* endif */

  /********************************************************************
   * DISPLAY THE PROGRAM HELP INFORMATION AS REQUIRED.                *
   ********************************************************************/

  if (CFFC_Primary_MPI_Processor() && help_flag) {
     cout << "Usage:\n";
     cout << "dusty2D [-v] [-h] [-i] [-b] [-f name] [-t name number]\n";
     cout << " -v (--version)  display version information\n";
     cout << " -h (--help)  show this help\n";
     cout << " -i (--inter)  execute in interactive mode (default)\n";
     cout << " -b (--batch)  execute in batch mode\n";
     cout << " -f name  use `name' input data file (`dusty2D.in' is default)\n";
     cout << " -t              run all available test suites\n" 
	  << " -t list         list available test suites\n"
	  << " -t regression   run all tests in batch mode with fault recovery\n"
	  << " -t test-suite-name [test-number] \n"
	  << "                 run all tests in test-suite-name or \n"
	  << "                 a particular test-number \n"
	  << "                 example: dusty2D -t MyTestSuit 3\n"
	  << " -path name      use 'name' as path to '/src_2D' directory\n"
	  << "                 use this option if you run UnitTesting framework\n"
	  << "                 from a directory different than '/src_2D'\n";  
     cout.flush();
  }
  if (help_flag) {
     CFFC_Finalize_MPI();
     return (0);
  } /* endif */

  /*********************************************************************
   * RUN UNIT TESTS AS REQUIRED AND STOP (tests run using only 1 CPU). *
   *********************************************************************/

#ifndef _NO_TUT_TESTING
  if (test_flag) {
     error_flag = Perform_UnitTesting(TestSuite, TestNumber, TestRootPath);
     CFFC_Finalize_MPI();
     return (error_flag);
  } /* endif */
#endif

  /********************************************************************
   * PERFORM REQUIRED CALCULATIONS.                                   *
   ********************************************************************/

  error_flag = Dusty2DQuadSolver(Input_File_Name_ptr,batch_flag);

  if (error_flag) {
    CFFC_Finalize_MPI();
    return error_flag;
  }

  /********************************************************************
   * FINALIZE MPI                                                     *
   ********************************************************************/
  CFFC_Finalize_MPI();

  /********************************************************************
   * TERMINATE PROGRAM EXECUTION                                      *
   ********************************************************************/
  if (CFFC_Primary_MPI_Processor() && !batch_flag) 
    cout << "\n\nDusty2D: Execution complete.\n";

  // Dusty2D program finished.
  return 0;

}
