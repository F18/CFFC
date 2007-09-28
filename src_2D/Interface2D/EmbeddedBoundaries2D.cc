/*******************************************************************
 *******************************************************************
 **********                                               **********
 ******            EMBEDDEDBOUNDARIES Version 1.00            ******
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

// Include required CFFC header files.

#include "EmbeddedBoundaries2D.h"
#include "EmbeddedBoundaries2D_Solver.h"
#include "EmbeddedBoundaries2D_FASMultigrid.h"
#include "EmbeddedBoundaries2D_Euler.h"
#include "EmbeddedBoundaries2D_NavierStokes.h"
#include "EmbeddedBoundaries2D_Dusty.h"
// #include "EmbeddedBoundaries2D_LESPremixed.h"
#include "../MPI/MPI.h"
#include "../ICEM/ICEMCFD.h"
#include "../UnitTesting/UnitTesting.h"

// Begin embeddedboundaries2D program.

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
      automate_flag,
      test_flag;

  // Name of program entered on the command line:
  char *command_name_ptr;

  // Title of code:
  char *program_title_ptr =
     "EmbeddedBoundaries2D: Numerical Solvers for various hyperbolic/elliptic systems of PDEs with embedded boundaries.";

  // Version of code:
  char *program_version_ptr = 
     "Version 1.00, UTIAS CFD & Propulsion Group, 1999-2007.";
  
  // Input file name:
  char *Input_File_Name_ptr = "embeddedboundaries2D.in";

  // Equation type indicator:
  char *Equation_Type_ptr = "Euler2D";
  char Equation_Type[256];

  // Input file stream:
  ifstream Input_File;

  // Unit testing flags:
  string TestSuite;
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
  automate_flag = 0;

  // Save the command line name of the program for future use.
  command_name_ptr = arg_ptr[0];

  // Parse and interpret command line arguments.  Note that there
  // are several different possible arguments which are:
  // 1) -v  lists program version information to standard output,
  // 2) -h  lists all possible optional arguments for program,
  // 3) -i  execute program in interactive mode (default mode),
  // 4) -b  execute program in batch mode,
  // 5) -pde type  sets type of partial differential equation to
  //               be solve to "type" (default is "Euler2D"),
  // 6) -f name  uses "name" as the input data file rather than
  //             the standard input data file "embeddedboundaries2D.in".
  // 7) -a  automatically sets up next case for unsteady flows.

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
      } else if (strcmp(arg_ptr[i-1],"-pde") == 0) {
        Equation_Type_ptr = arg_ptr[i];
        strcpy(Equation_Type, Equation_Type_ptr);
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
      } else if (strcmp(arg_ptr[i],"-a") == 0) {
	automate_flag = 1;
      } else if (strcmp(arg_ptr[i], "-t") == 0) {
	test_flag=1;
	if (num_arg-1>i){
	  TestSuite= arg_ptr[i+1];
	} /* endif */
	if (num_arg-1>i+1){
	  TestNumber = atoi(arg_ptr[i+2]);
	} /* endif */
	break;
      } else {
        error_flag = 1;
      }
    }
  }

  // Display command line argument error message and terminate the
  // program as required.
  if (error_flag) {
    cout << "\n EmbeddedBoundaries2D ERROR: Invalid command line argument.\n";
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
     cout << "Built using MV++, SparseLib++, IML++, BPKIT, and FFTW Libraries\n";
     cout.flush();
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
     cout << "embeddedboundaries2D [-v] [-h] [-i] [-b] [-pde type] [-f name] [-t name number]\n";
     cout << " -v (--version)  display version information\n";
     cout << " -h (--help)  show this help\n";
     cout << " -i (--inter)  execute in interactive mode (default)\n";
     cout << " -b (--batch)  execute in batch mode\n";
     cout << " -pde type  solve `type' PDEs (`Euler2D' is default)\n";
     cout << " -f name  use `name' input data file (`embeddedboundaries2D.in' is default)\n";
     cout << " -t              run all available test suites\n" 
	  << " -t list         list available test suites\n"
	  << " -t regression   run all tests in batch mode with fault recovery\n"
	  << " -t test-suite-name [test-number] \n"
	  << "                 run all tests in test-suite-name or \n"
	  << "                 a particular test-number \n"
	  << "                 example: embeddedboundaries2D -t MyTestSuit 3\n";
     cout.flush();
  }
  if (help_flag) {
     CFFC_Finalize_MPI();
     return (0);
  } /* endif */

  /*********************************************************************
   * RUN UNIT TESTS AS REQUIRED AND STOP (tests run using only 1 CPU). *
   *********************************************************************/

  if (test_flag) {
     error_flag = Perform_UnitTesting(TestSuite, TestNumber);
     CFFC_Finalize_MPI();
     return (error_flag);
  } /* endif */

  /********************************************************************
   * PERFORM REQUIRED CALCULATIONS.                                   *
   ********************************************************************/

  // Used for burning crv7.
  if (CFFC_Primary_MPI_Processor() && automate_flag) {
    cout << "\n Generating input file(s) and script for the next frame.";
    std::system("/nfs/fe01/d1/cfd/jai/bin/setup.pl > ./setup.log");
    cout << endl;
  }
  CFFC_Barrier_MPI();

  // For the PDE of interest, solve the corresponding initial-boundary
  // value problem(s)/boundary value problem(s) (IBVP/BVP) using the 
  // appropriate equation solver.
  if (strcmp(Equation_Type,"Euler2D") == 0) {
    error_flag = EmbeddedBoundaries2D_Solver<Euler2D_cState,
                                             Euler2D_pState,
                                             Euler2D_Quad_Block,
                                             Euler2D_Input_Parameters>(Input_File_Name_ptr,
								       batch_flag);
  } else if (strcmp(Equation_Type,"NavierStokes2D") == 0) {
    error_flag = EmbeddedBoundaries2D_Solver<NavierStokes2D_cState,
                                             NavierStokes2D_pState,
                                             NavierStokes2D_Quad_Block, 
                                             NavierStokes2D_Input_Parameters>(Input_File_Name_ptr,
									      batch_flag);
  } else if (strcmp(Equation_Type,"Dusty2D") == 0) {
    error_flag = EmbeddedBoundaries2D_Solver<Dusty2D_cState,
                                             Dusty2D_pState,
                                             Dusty2D_Quad_Block, 
                                             Dusty2D_Input_Parameters>(Input_File_Name_ptr,
 								       batch_flag);
   } else if (strcmp(Equation_Type,"Gaussian2D") == 0) {
     error_flag = EmbeddedBoundaries2D_Solver<Gaussian2D_cState,
                                              Gaussian2D_pState,
                                              Gaussian2D_Quad_Block, 
                                              Gaussian2D_Input_Parameters>(Input_File_Name_ptr,
									   batch_flag);
//    } else if (strcmp(Equation_Type,"LESPremix2D") == 0) {
//      error_flag = EmbeddedBoundaries2D_Solver<LESPremix2D_cState,
//                                               LESPremix2D_pState,
//                                               LESPremix2D_Quad_Block, 
//                                               LESPremix2D_Input_Parameters>(Input_File_Name_ptr,
// 									    batch_flag);
  } else {
    if (CFFC_Primary_MPI_Processor() && !batch_flag)
      cout << "\n\n EmbeddedBoundaries2D ERROR: Specified equation set is not supported.\n";
    error_flag = 1;
  }

  if (error_flag) {
    CFFC_Finalize_MPI();
    return error_flag;
  }

  // Used for burning crv7.
  if (CFFC_Primary_MPI_Processor() && automate_flag) {
    cout << "\n\n Copying restart files for the next frame.";
    std::system("/nfs/fe01/d1/cfd/jai/bin/copy.pl >> ./setup.log");
    //std::system("/nfs/fe01/d1/cfd/jai/bin/setup_next.pl > ./setup.log");
  }
  CFFC_Barrier_MPI();

  /********************************************************************
   * FINALIZE MPI                                                     *
   ********************************************************************/
  CFFC_Finalize_MPI();

  /********************************************************************
   * TERMINATE PROGRAM EXECUTION                                      *
   ********************************************************************/
  if (CFFC_Primary_MPI_Processor() && !batch_flag) 
    cout << "\n\nEmbeddedBoundaries2D: Execution complete.\n";

  // Embeddedboundaries2D program finished.
  return 0;

}
