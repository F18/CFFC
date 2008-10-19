
/*******************************************************************
 *******************************************************************
 **********          Solvers1D - Built using CFFC            *********
 ******                                                       ******   
 ****                                                           ****
 ******       UTIAS, CFD & Propulsion Group, 1999-2007        ******
 **********                                                *********
 *******************************************************************
 ******************************************************************* 

         This computer program can be used to solve a selected
 set of linear and nonlinear partial differential equations (PDEs)
 on  one-dimensional spatial domains using a variety of solution 
 techniques appropriate for the problem.
 The PDEs that can be treated include:

 -- the scalar advection equation,
 -- the inviscid Burger's equation,
 -- the scalar heat equation,
 -- the hyperbolic heat equations (Maxwell-Cattaneo equations),
 -- the Euler equations of continuum fluid dynamics,
 -- the ideal MHD equations,

 The program is built using the CFFC library.

 *******************************************************************/

/* Include the header files defining various variable types and data
   structures, classes, functions, operators, global variables,
   etc.. */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>

using namespace std;

// Include CFFC header files.
#include "Scalar1D.h"
#include "Heat1D.h"
#include "HyperHeat1D.h"
#include "MHD1D.h"
#include "../Euler1D/Euler1D.h"
#include "../UnitTesting/UnitTesting1D.h"

/* Begin Solvers1D program. */

int main(int num_arg, char *arg_ptr[]) {

  /********************************************************  
   * VARIABLE DECLARATIONS                                *
   ********************************************************/

  // Command line flags:  
  int version_flag, 
      help_flag, 
      batch_flag, 
      pde_flag,
      file_flag, 
      error_flag,
      mpirun_flag,
      test_flag;

  // Other local integer variables:
  int i;

  // Name of program entered on the command line:
  char *command_name_ptr;

  // Title of code:
  char *program_title_ptr = 
     "Solvers1D: Stand alone solver for all implemented CFFC equations in 1 space dimension.";

  // Version of code:
  char *program_version_ptr = 
     "Version 1.00, UTIAS CFD & Propulsion Group, 1999-2007.";
  
  // Input file name:
  char *Input_File_Name_ptr = "euler1D.in";

  // Equation type indicator:
  char *Equation_Type_ptr = "Euler1D";
  char Equation_Type[256];

  // Input file stream:
  ifstream Input_File;

  // Unit testing flags:
  string TestSuite, TestRootPath;
  int TestNumber=0;

  /********************************************************  
   * LOCAL VARIABLE INITIALIZATIONS                       *
   ********************************************************/

  /* Set default equation type. */

  strcpy(Equation_Type, Equation_Type_ptr);

  /********************************************************  
   * PARSE COMMAND LINE ARGUMENTS                         *
   ********************************************************/

  /* Initialize command line flags. */

  version_flag = 0;
  help_flag = 0;
  batch_flag = 0;
  pde_flag = 0;
  file_flag = 0;
  error_flag = 0;
  mpirun_flag = 0;
  test_flag = 0;

  /* Save the command line name of the program for future use. */

  command_name_ptr = arg_ptr[0];
  
  /* Parse and interpret command line arguments.  Note that there
     are several different possible arguments which are:
     1) -v  lists program version information to standard output,
     2) -h  lists all possible optional arguments for program,
     3) -i  execute program in interactive mode (default mode),
     4) -b  execute program in batch mode,
     5) -f name  uses "name" as the input data file rather than
                 the standard input data file "euler1D.in". */

  if (num_arg >= 2) {
    for (i = 1; i <= num_arg - 1; ++i) {
      if (strcmp(arg_ptr[i], "-v") == 0 ||
          strcmp(arg_ptr[i], "--version") == 0) {
        version_flag = 1;
      } else if (strcmp(arg_ptr[i], "-h") == 0 ||
                 strcmp(arg_ptr[i], "--help") == 0) {
        help_flag = 1;
      } else if (strcmp(arg_ptr[i], "-i") == 0 ||
                 strcmp(arg_ptr[i], "--inter") == 0) {
        batch_flag = 0;
      } else if (strcmp(arg_ptr[i], "-b") == 0 ||
                 strcmp(arg_ptr[i], "--batch") == 0) {
        batch_flag = 1;
      } else if (strcmp(arg_ptr[i], "-pde") == 0) {
        pde_flag = 1;
//       } else if (strcmp(arg_ptr[i-1], "-pde") == 0) {
//         Equation_Type_ptr = arg_ptr[i];
//         strcpy(Equation_Type, Equation_Type_ptr);
      } else if (strcmp(arg_ptr[i], "-f") == 0) {
        file_flag = 1;
      } else if (strcmp(arg_ptr[i-1], "-f") == 0) {
        Input_File_Name_ptr = arg_ptr[i];
      } else if (strcmp(arg_ptr[i], "-p4pg") == 0) {
        mpirun_flag = 1;
      } else if (strcmp(arg_ptr[i-1], "-p4pg") == 0) {
        //mpirun_flag = 1;
      } else if (strcmp(arg_ptr[i], "-p4wd") == 0) {
        mpirun_flag = 1;
      } else if (strcmp(arg_ptr[i-1], "-p4wd") == 0) {
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
      } /* endif */
    } /* endfor */
  } /* endif */

  /* Display command line argument error message and
     terminate the program as required. */

  if (error_flag) {
    cout << "\n Solvers1D ERROR: Invalid command line argument.\n";
    return (error_flag);
  } /* endif */


  /******************************************************************
   * DISPLAY THE PROGRAM TITLE AND VERSION INFORMATION AS REGUIRED. *
   ******************************************************************/

  if ((version_flag || help_flag || !batch_flag)) {
     cout << '\n' << program_title_ptr << '\n';
     cout << program_version_ptr << '\n';
     cout << "Built using " << CFFC_Version() << "\n";
     cout << "Built using MV++, SparseLib++, IML++, and BPKIT Libraries\n";
     cout.flush();
  } /* endif */
  if (version_flag) {
     return (0);
  } /* endif */

  /******************************************************************
   * DISPLAY THE PROGRAM HELP INFORMATION AS REQUIRED.              *
   ******************************************************************/

  if ( help_flag) {
     cout << "Usage:\n";
     cout << "solvers1D [-v] [-h] [-i] [-b] [-f name] [-t name number]\n";
     cout << " -v (--version)  display version information\n";
     cout << " -h (--help)     show this help\n";
     cout << " -i (--inter)    execute in interactive mode (default)\n";
     cout << " -b (--batch)    execute in batch mode\n";
     cout << " -pde type       solve `type' PDEs (`Euler1D' is default)\n";
     cout << " -f name         use `name' input data file (`euler1D.in' is default)\n";
     cout << " -t              run all available test suites\n" 
	  << " -t list         list available test suites\n"
	  << " -t regression   run all tests in batch mode with fault recovery\n"
	  << " -t test-suite-name [test-number] \n"
	  << "                 run all tests in test-suite-name or \n"
	  << "                 a particular test-number \n"
	  << "                 example: solvers1D -t MyTestSuit 3\n"
	  << " -path name      use 'name' as path to '/src_2D' directory\n"
	  << "                 use this option if you run UnitTesting framework\n"
	  << "                 from a directory different than '/src_2D'\n";     
     cout.flush();
  } /* endif */
  if (help_flag) {
     return (0);
  } /* endif */

  /*********************************************************************
   * RUN UNIT TESTS AS REQUIRED AND STOP (tests run using only 1 CPU). *
   *********************************************************************/

#ifndef _NO_TUT_TESTING
  if (test_flag) {
     error_flag = Perform_UnitTesting(TestSuite, TestNumber, TestRootPath);
     return (error_flag);
  } /* endif */
#endif

  /***********************************************************  
   * PERFORM REQUIRED CALCULATIONS.                          *
   ***********************************************************/

  /***********************************************************  
   * PERFORM REQUIRED CALCULATIONS.                          *
   ***********************************************************/

  /* For the PDE of interest, solve the corresponding initial-boundary
     value problem(s)/boundary value problem(s) (IBVP/BVP) using the 
     appropriate equation solver. */

  /* Scalar1D */
  if (strcmp(Equation_Type, "Scalar1D") == 0) {
    error_flag = Scalar1DSolver(Input_File_Name_ptr,
				batch_flag);

  /* Heat1D */
  } else if (strcmp(Equation_Type, "Heat1D") == 0) {
    error_flag = Heat1DSolver(Input_File_Name_ptr,
			      batch_flag);

  /* HyperHeat1D */
  } else if (strcmp(Equation_Type, "HyperHeat1D") == 0) {
    error_flag = HyperHeat1DSolver(Input_File_Name_ptr,
				   batch_flag);

  /* Euler1D */
  } else if (strcmp(Equation_Type, "Euler1D") == 0) {
    error_flag = Euler1DSolver(Input_File_Name_ptr,
			       batch_flag);
  
  /* MHD1D */
  } else if (strcmp(Equation_Type, "MHD1D") == 0) {
    error_flag = MHD1DSolver(Input_File_Name_ptr,
			     batch_flag);
  } /* endif */

  
  if (error_flag) {
     return (error_flag);
  } /* endif */


  /********************************************************  
   * TERMINATE PROGRAM EXECUTION                          *
   ********************************************************/

  if (!batch_flag) 
     cout << "\n\nSolvers1D: Execution complete.\n";

  //Ending properly
  return (0);

/* End Solvers1D program. */

}
