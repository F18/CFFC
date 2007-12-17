/*******************************************************************
 *******************************************************************
 **********                                               **********
 ******         FANS3DThermallyPerfect Version 1.00          *******
 ****                       (07/19/07)                          ****
 ****                                                           ****
 ****        Solver for FANS Equations Governing Flow of        ****
 ****     Thermally Perfect Gases in Three Space Dimensions     ****
 ****                                                           ****
 ******       UTIAS, CFD & Propulsion Group, 1999-2007        ******
 **********                                               **********
 *******************************************************************
 ******************************************************************* 
 
 Include the header files defining various variable types and data
 structures, classes, functions, operators, global variables,
 etc.. */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>

using namespace std;

// Include CFFC header files.

#include "FANS3DThermallyPerfect.h"
#include "../UnitTesting/UnitTesting.h"

/* Begin FANS3DThermallyPerfect program. */

int main(int num_arg, char *arg_ptr[]) {
   
   /********************************************************  
    * VARIABLE DECLARATIONS                                *
    ********************************************************/

   // Command line flags:  
   int version_flag, 
       help_flag, 
       batch_flag, 
       file_flag, 
       error_flag,
       mpirun_flag,
       test_flag;

   // Other local integer variables:
   int i;

   // Name of program entered on the command line:
   char *command_name_ptr;

   char *program_title_ptr = 
     "FANS3DThermallyPerfect: Solver for the FANS Equations Governing Thermally Perfect Gases in 3 Space Dimensions.";

   // Version of code:
   char *program_version_ptr = 
      "Version 1.00, UTIAS CFD & Propulsion Group, 1999-2007.";
  
   // Input file name:
   char *Input_File_Name_ptr = "fans3d.in";

   // Input file stream:
   ifstream Input_File;

   // Testing commands
   string TestSuite, TestRootPath;
   int TestNumber=0;

   /********************************************************  
    * PARSE COMMAND LINE ARGUMENTS                         *
    ********************************************************/

   /* Initialize command line flags. */

   version_flag = 0;
   help_flag = 0;
   batch_flag = 0;
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
      the standard input data file "fans3D.in". */

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

   if (error_flag) {
     cout << "\n FANS3DThermallyPerfect ERROR: Invalid command line argument.\n";
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
      cout << '\n' << program_title_ptr << '\n';
      cout << program_version_ptr << '\n';
      cout << "Built using " << CFFC_Version() << "\n";
      cout << CFFC_Version_MPI() << "\n";
      cout << Cantera_Version() << "\n";
      cout << ICEMCFD_Version() << "\n";
      cout << "Built using MV++, SparseLib++, IML++, BPKIT, and FFTW Libraries\n";
      cout << "Built using CEA Thermodynamic and Transport Data, NASA Glenn Research Center\n";
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
      cout << "fans3Dthermallyperfect [-v] [-h] [-i] [-b] [-f name] [-t name number]\n";
      cout << " -v (--version)  display version information\n";
      cout << " -h (--help)  show this help\n";
      cout << " -i (--inter)  execute in interactive mode (default)\n";
      cout << " -b (--batch)  execute in batch mode\n";
      cout << " -f name  use `name' input data file (`fans3D.in' is default)\n";
      cout << " -t              run all available test suites\n" 
	   << " -t list         list available test suites\n"
	   << " -t regression   run all tests in batch mode with fault recovery\n"
	   << " -t test-suite-name [test-number] \n"
	   << "                 run all tests in test-suite-name or \n"
	   << "                 a particular test-number \n"
 	   << "                 example: fans3Dthermallyperfect -t MyTestSuit 3\n"
	   << " -path name      use 'name' as path to '/src_3D' directory\n"
	   << "                 use this option if you run UnitTesting framework\n"
	   << "                 from a directory different than '/src_3D'\n";    
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

   /* Solve the corresponding initial-boundary value problem(s)/boundary 
      value problem(s) (IBVP/BVP). */

   error_flag = HexaSolver<FANS3D_ThermallyPerfect_KOmega_pState, 
                           FANS3D_ThermallyPerfect_KOmega_cState>
               (Input_File_Name_ptr, batch_flag);
 
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
      cout << "\n\nFANS3DThermallyPerfect: Execution complete.\n\n";

   //Ending properly
   return (0);
   
/* End FANS3DThermallyPerfect program. */

}
