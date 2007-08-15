/*******************************************************************
 *******************************************************************
 **********         Gaussian2D - Built using CFFC          *********
 ******                                                       ******   
 ****                                                           ****
 ******       UTIAS, CFD & Propulsion Group, 1999-2004        ******
 **********                                                *********
 *******************************************************************
 ******************************************************************* 

 This program is a standalone version of of the Gaussian2D Solver that 
 is part of the PDES++ library of solvers.  

 The program is built using the CFFC library.

 *******************************************************************/

/* Include the header files defining various variable types and data
   structures, classes, functions, operators, global variables,
   etc.. */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#include "../MPI/MPI.h"
#include "../ICEM/ICEMCFD.h"

#include "Gaussian2DCartesian.h"
#include "Gaussian2DQuad.h"


/* Begin PDES++ program. */
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
      mpirun_flag;

  // Other local integer variables:
  int i;

  // Name of program entered on the command line:
  char *command_name_ptr;

  // Title of code:
  char *program_title_ptr = 
     "Gaussian2D: Stand alone Solver.";

  // Version of code:
  char *program_version_ptr = 
     "Beta Version, UTIAS CFD & Propulsion Group, 1999-2003.";
  
  // Input file name:
  char *Input_File_Name_ptr = "Gaussian2D.in";

  // Equation type indicator:
  char *Equation_Type_ptr = "Gaussian2D";
  char Equation_Type[256];

  // Input file stream:
  ifstream Input_File;

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

  /* Save the command line name of the program for future use. */

  command_name_ptr = arg_ptr[0];
  
  /* Parse and interpret command line arguments.  Note that there
     are several different possible arguments which are:

     1) -v  lists PDES++ version information to standard output,
     2) -h  lists all possible optional arguments for PDES++,
     3) -i  execute PDES++ in interactive mode (default mode),
     4) -b  execute PDES++ in batch mode,
     5) -pde type   sets type of partial differential equation
                    to be solve to "type" (default is "Gaussian2D"),
     6) -f name  uses "name" as the input data file rather than
                 the standard input data file "pdes++.in". */

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
      } else {
        error_flag = 1;
      } /* endif */
    } /* endfor */
  } /* endif */

  /* Display command line argument error message and
     terminate the program as required. */

  if (error_flag) {
    cout << "\n ERROR: Invalid command line argument.\n";
    return (error_flag);
  } /* endif */

  /********************************************************  
   * INITIALIZE MPI                                       *
   ********************************************************/

  CFFC_Initialize_MPI(num_arg, arg_ptr);
  if (!CFFC_Primary_MPI_Processor()) batch_flag = 1;

  /******************************************************************
   * DISPLAY THE PROGRAM TITLE AND VERSION INFORMATION AS REGUIRED. *
   ******************************************************************/

  if (CFFC_Primary_MPI_Processor() && (version_flag || help_flag || !batch_flag)) {
     cout << '\n' << program_title_ptr << '\n';
     cout << program_version_ptr << '\n';
     cout << "Built using " << CFFC_Version() << "\n";
     cout << CFFC_Version_MPI() << "\n";
     cout << ICEMCFD_Version() << "\n";
     cout << "Built using MV++, SparseLib++, IML++, and BPKIT Libraries.\n";
     cout.flush();
     if (version_flag) return (0);
  } /* endif */

  /******************************************************************
   * DISPLAY THE PROGRAM HELP INFORMATION AS REQUIRED.              *
   ******************************************************************/

  if (CFFC_Primary_MPI_Processor() && help_flag) {
     cout << "Usage:\n";
     cout << "pdes++ [-v] [-h] [-i] [-b] [-pde type] [-f name]\n";
     cout << " -v (--version)  display version information\n";
     cout << " -h (--help)  show this help\n";
     cout << " -i (--inter)  execute in interactive mode (default)\n";
     cout << " -b (--batch)  execute in batch mode\n";
     cout << " -pde type  solve `type' PDEs (`Gaussian2D' is default)\n";
     cout << " -f name  use `name' input data file (`pdes++.in' is default)\n";
     cout.flush();
     return (0);
  } /* endif */

  /***********************************************************  
   * PERFORM REQUIRED CALCULATIONS.                          *
   ***********************************************************/


  if (strcmp(Equation_Type, "Gaussian2D_Cartesian") == 0){
      if(CFFC_Primary_MPI_Processor()) {
        error_flag = Gaussian2DCartesianSolver(Input_File_Name_ptr,
	  				       batch_flag);
      } /* endif */
      CFFC_Broadcast_MPI(&error_flag, 1);

  /* Gaussian2D */
  } else if(strcmp(Equation_Type, "Gaussian2D") == 0){
      error_flag = Gaussian2DQuadSolver(Input_File_Name_ptr,
	  				       batch_flag);
  }/* endif */

 
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
     cout << "\n\nGaussian2D: Execution complete.\n";
  return (0);

/* End PDES++ program. */

}
