/*! \mainpage

\verbatim
 *******************************************************************
 *******************************************************************
 **********                                               **********
 ******                 CFFC2D Version 0.01                  *******
 ****                       (07/19/07)                          ****
 ****                                                           ****
 ****  Computational Framework for Fluids and Combustion (CFFC) ****
 ****              in One and Two Space Dimensions              ****
 ****                                                           ****
 ******       UTIAS, CFD & Propulsion Group, 1999-2007        ******
 **********                                               **********
 *******************************************************************
 ******************************************************************* 
\endverbatim
*/

/*! \file cffc2D.cc

    \brief Main executable for CFFC2D framework.

 This computer program can be used to solve a selected
 set of linear and non-linear partial differential equations (PDEs)
 on  one- and two-dimensional spatial domains using a
 variety of solution techniques appropriate for the problem.  The
 PDEs that can be treated include:

 - the scalar advection equation,
 - the inviscid Burger's equation,
 - the Poisson equation,
 - the scalar heat equation,
 - the hyperbolic heat equations (Maxwell-Cattaneo equations),
 - the advection diffusion equation,
 - the Euler equations of continuum fluid dynamics,
 - the Euler/dusty-gas equations,
 - the Euler multispecies with finite rate chemistry equation system,
 - the Navier-Stokes equations of continuum fluid dynamics,
 - the ideal MHD equations,
 - the 10-moment equations of the Gaussian closure,
 - the generalized 10-moment equations with heat conduction,
 - the 35-moment equations of the Gaussian-based 35-moment closure.

 The program is built using the CFFC library.

*/

/* Include the header files defining various variable types and data
   structures, classes, functions, operators, global variables,
   etc.. */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>

using namespace std;

#include "Solvers1D/Scalar1D.h"
#include "Solvers1D/Heat1D.h"
#include "Solvers1D/HyperHeat1D.h"
#include "Euler1D/Euler1D.h"
#include "Solvers1D/MHD1D.h"
#include "AdvectDiffuse2D/AdvectDiffuse2DQuad.h"
#include "AdvectDiffuse2D/AdvectDiffuse2DQuad_NKS.h"
#include "Euler2D/Euler2DQuad.h"
#include "NavierStokes2D/NavierStokes2DQuad.h"
#include "Chem2D/Chem2DQuad.h"
#include "HighTemp2D/HighTemp2DQuad.h"
#include "Dusty2D/Dusty2DQuad.h"
#include "LESPremixed2D/LESPremixed2DQuad.h"
#include "Gaussian2D/Gaussian2DCartesian.h"
#include "Gaussian2D/Gaussian2DQuad.h"
#include "Ion5Moment2D/Ion5Moment2DQuad.h"
#include "Electrostatic2D/Electrostatic2DQuad.h"
#include "Rte2D/Rte2DQuad.h"
#include "LevelSet2D/LevelSet2DQuad.h"
#include "MPI/MPI.h"
#include "ICEM/ICEMCFD.h"
#include "UnitTesting/UnitTesting.h"

/* Begin CFFC2D program. */

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
     "CFFC2D: Computational Framework for Fluids and Combustion in 1 and 2 Space Dimensions.";

  // Version of code:
  char *program_version_ptr = 
     "Version 1.00, UTIAS CFD & Propulsion Group, 1999-2007.";
  
  // Input file name:
  char *Input_File_Name_ptr = "cffc2d.in";

  // Equation type indicator:
  char *Equation_Type_ptr = "Euler2D";
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
     5) -pde type   sets type of partial differential equation
                    to be solve to "type" (default is "Euler2D"),
     6) -f name  uses "name" as the input data file rather than
                 the standard input data file "cffc2d.in". */

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
      } else if (strcmp(arg_ptr[i-1], "-pde") == 0) {
        Equation_Type_ptr = arg_ptr[i];
        strcpy(Equation_Type, Equation_Type_ptr);
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
    cout << "\n CFFC2D ERROR: Invalid command line argument.\n";
    return (error_flag);
  } /* endif */

  /********************************************************  
   * INITIALIZE MPI                                       *
   ********************************************************/

  CFFC_Initialize_MPI(num_arg, arg_ptr);
  if (!CFFC_Primary_MPI_Processor()) batch_flag = 1;

  /******************************************************************
   * DISPLAY THE PROGRAM TITLE AND VERSION INFORMATION AS REQUIRED. *
   ******************************************************************/

  if (CFFC_Primary_MPI_Processor() && (version_flag || help_flag || !batch_flag)) {
     cout << '\n' << program_title_ptr << '\n'
	  << program_version_ptr << '\n'
	  << "Built using " << CFFC_Version() << "\n"
	  << CFFC_Version_MPI() << "\n"
	  << Cantera_Version() << "\n"
	  << ICEMCFD_Version() << "\n"
	  << "Built using MV++, SparseLib++, IML++, BPKIT, and FFTW Libraries\n"
	  << "Built using CEA Thermodynamic and Transport Data, NASA Glenn Research Center\n";
     cout.flush();
  } /* endif */

  if (version_flag) {
     CFFC_Finalize_MPI();
     return (0);
  } /* endif */

  /******************************************************************
   * DISPLAY THE PROGRAM HELP INFORMATION AS REQUIRED.              *
   ******************************************************************/

  if (CFFC_Primary_MPI_Processor() && help_flag) {
     cout << "Usage:\n";
     cout << "cffc2d [-v] [-h] [-i] [-b] [-pde type] [-f name] [-t name number]\n";
     cout << " -v (--version)  display version information\n";
     cout << " -h (--help)  show this help\n";
     cout << " -i (--inter)  execute in interactive mode (default)\n";
     cout << " -b (--batch)  execute in batch mode\n";
     cout << " -pde type  solve `type' PDEs (`Euler2D' is default)\n";
     cout << " -f name  use `name' input data file (`cffc2d.in' is default)\n";
     cout << " -t              run all available test suites\n" 
	  << " -t list         list available test suites\n"
	  << " -t regression   run all tests in batch mode with fault recovery\n"
	  << " -t test-suite-name [test-number] \n"
	  << "                 run all tests in test-suite-name or \n"
	  << "                 a particular test-number \n"
 	  << "                 example: cffc2D -t MyTestSuit 3\n"
	  << " -path name      use 'name' as path to '/src_2D' directory\n"
	  << "                 use this option if you run UnitTesting framework\n"
	  << "                 from a directory different than '/src_2D'\n";  
     cout.flush();
  } /* endif */
  if (help_flag) {
     CFFC_Finalize_MPI();
     return (0);
  } /* endif */

  /*********************************************************************
   * RUN UNIT TESTS AS REQUIRED AND STOP (tests run using only 1 CPU). *
   *********************************************************************/

  if (test_flag) {
     error_flag = Perform_UnitTesting(TestSuite, TestNumber, TestRootPath);
     CFFC_Finalize_MPI();
     return (error_flag);
  } /* endif */

  /***********************************************************  
   * PERFORM REQUIRED CALCULATIONS.                          *
   ***********************************************************/

  /* For the PDE of interest, solve the corresponding initial-boundary
     value problem(s)/boundary value problem(s) (IBVP/BVP) using the 
     appropriate equation solver. */

  /* Scalar1D */
  if (strcmp(Equation_Type, "Scalar1D") == 0) {
      if (CFFC_Primary_MPI_Processor()) {
         error_flag = Scalar1DSolver(Input_File_Name_ptr,
                                     batch_flag);
      } /* endif */
      CFFC_Broadcast_MPI(&error_flag, 1);

  /* Heat1D */
  } else if (strcmp(Equation_Type, "Heat1D") == 0) {
      if (CFFC_Primary_MPI_Processor()) {
         error_flag = Heat1DSolver(Input_File_Name_ptr,
                                   batch_flag);
      } /* endif */
      CFFC_Broadcast_MPI(&error_flag, 1);
  
  /* HyperHeat1D */
  } else if (strcmp(Equation_Type, "HyperHeat1D") == 0) {
      if (CFFC_Primary_MPI_Processor()) {
         error_flag = HyperHeat1DSolver(Input_File_Name_ptr,
                                        batch_flag);
      } /* endif */
      CFFC_Broadcast_MPI(&error_flag, 1);
  
  /* Euler1D */
  } else if (strcmp(Equation_Type, "Euler1D") == 0) {
      if (CFFC_Primary_MPI_Processor()) {
         error_flag = Euler1DSolver(Input_File_Name_ptr,
                                    batch_flag);
      } /* endif */
      CFFC_Broadcast_MPI(&error_flag, 1);
  
  /* MHD1D */
  } else if (strcmp(Equation_Type, "MHD1D") == 0) {
      if (CFFC_Primary_MPI_Processor()) {
         error_flag = MHD1DSolver(Input_File_Name_ptr,
                                  batch_flag);
      } /* endif */
      CFFC_Broadcast_MPI(&error_flag, 1);

  /* AdvectDiffuse2D */
  } else if (strcmp(Equation_Type, "AdvectDiffuse2D") == 0) {
    error_flag = AdvectDiffuse2DQuadSolver(Input_File_Name_ptr,
					   batch_flag);

  /* Euler2D */
  } else if (strcmp(Equation_Type, "Euler2D") == 0) {
    error_flag = Euler2DQuadSolver(Input_File_Name_ptr,
                                     batch_flag);

  /* Dusty2D */
  } else if (strcmp(Equation_Type, "Dusty2D") == 0) {
    error_flag = Dusty2DQuadSolver(Input_File_Name_ptr,
                                   batch_flag);
  
  /* Ion5Moment2D */
  } else if (strcmp(Equation_Type, "Ion5Moment2D") == 0) {
      error_flag = Ion5Moment2DQuadSolver(Input_File_Name_ptr,
                                          batch_flag);
  
  /* ElectroStatic2D */
  } else if (strcmp(Equation_Type, "ElectroStatic2D") == 0) {
      error_flag = Electrostatic2DQuadSolver(Input_File_Name_ptr,
                                          batch_flag);

  /* LESPremixed2D */
  } else if (strcmp(Equation_Type, "LESPremixed2D") == 0) {
      error_flag = LESPremixed2DQuadSolver(Input_File_Name_ptr,
					   batch_flag);

  /* LevelSet2D */
  } else if (strcmp(Equation_Type, "LevelSet2D") == 0) {
      error_flag = LevelSet2DQuadSolver(Input_File_Name_ptr,
					batch_flag);
  
  /* Chem2D */
  } else if (strcmp(Equation_Type, "Chem2D") == 0) {
      error_flag = Chem2DQuadSolver(Input_File_Name_ptr,
			  	    batch_flag);

  /* NavierStokes2D */
  } else if (strcmp(Equation_Type, "NavierStokes2D") == 0) {
      error_flag = NavierStokes2DQuadSolver(Input_File_Name_ptr,
			  	            batch_flag);

  /* Gaussian2D_Cartesian */
  } else if (strcmp(Equation_Type, "Gaussian2D_Cartesian") == 0) {
    if(CFFC_Primary_MPI_Processor()) {
      error_flag = Gaussian2DCartesianSolver(Input_File_Name_ptr,
					     batch_flag);
    } /* endif */
    CFFC_Broadcast_MPI(&error_flag, 1);

  /* Gaussian2D */
  } else if(strcmp(Equation_Type, "Gaussian2D") == 0) {
      error_flag = Gaussian2DQuadSolver(Input_File_Name_ptr,
	  			        batch_flag);

  /* HighTemp2D */
  } else if (strcmp(Equation_Type, "HighTemp2D") == 0) {
      error_flag = HighTemp2DQuadSolver(Input_File_Name_ptr,
			  	        batch_flag);

  /* Rte2D */
  } else if (strcmp(Equation_Type, "Rte2D") == 0) {
      error_flag = Rte2DQuadSolver(Input_File_Name_ptr,
                                   batch_flag);

  } /* endif */

  if (error_flag) {
     CFFC_Finalize_MPI();
     return (error_flag);
  } /* endif */

  /********************************************************  
   * FINALIZE MPI                                         *
   ********************************************************/

  CFFC_Finalize_MPI();

  /********************************************************  
   * TERMINATE PROGRAM EXECUTION                          *
   ********************************************************/

  if (CFFC_Primary_MPI_Processor() && !batch_flag) 
     cout << "\n\nCFFC2D: Execution complete.\n";

  //Ending properly
  return (0);

/* End CFFC program. */

}
