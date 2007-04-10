/* Include the header files defining various variable types and data
   structures, classes, functions, operators, global variables,
   etc.. */

#include "Reconstruction/Reconstruction.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <Test_Run.h>
#include "Common/SourceRevisionData.h"

using namespace std;

/* Begin Reconstruction main program. */

int main(int num_arg, char *arg_ptr[]){

  /********************************************************  
   * VARIABLE DECLARATIONS                                *
   ********************************************************/
  // Command line flags:  
  int version_flag, 
    help_flag, 
    file_flag, 
    error_flag,
    test_flag;

  vector<string> arguments(num_arg);
  // Other local integer variables:
  int i;

  // Name of program entered on the command line:
  char *command_name_ptr;

  // Title of code:
  char *program_title_ptr = 
    "'Reconstruction++' is a program for reconstructing functions \nin one-, two- and three-dimensions based on cell averages.";

  // Version of code:
  char *program_version_ptr = 
    "UTIAS CFD & Propulsion Group, 1999 -";
  
  // Input file name:
  char *Input_File_Name_ptr = "reconstruct_DEFAULT.in";

  // Equation type indicator:
  char *Reconstruct_Type_ptr = "1D";
  char Reconstruct_Type[256];

  // Input file stream:
  ifstream Input_File;

  // Testing commands
  string TestSuit;
  int TestNumber=0;

  /********************************************************  
   * LOCAL VARIABLE INITIALIZATIONS                       *
   ********************************************************/

  /* Set default equation type. */

  strcpy(Reconstruct_Type, Reconstruct_Type_ptr);

  /********************************************************  
   * PARSE COMMAND LINE ARGUMENTS                         *
   ********************************************************/

  /* Initialize command line flags. */

  version_flag = 0;
  help_flag = 0;
  file_flag = 0;
  error_flag = 0;
  test_flag = 0;

  /* Save the command line name of the program for future use. */

  command_name_ptr = arg_ptr[0];
  
  /* Parse and interpret command line arguments.  Note that there
     are several different possible arguments which are:

     1) -v  lists RECONSTRUCT version information to standard output,
     2) -h  lists all possible optional arguments for RECONSTRUCT,
     5) -rec type   sets type of reconstruction 
     to be solve to "type" (default is "Reconstruct1D"),
     6) -f name  uses "name" as the input data file rather than
     the standard input data file "reconstruct.in". */

  if (num_arg >= 2) {
    for (i = 1; i <= num_arg - 1; ++i) {
      if (strcmp(arg_ptr[i], "-v") == 0 ||
          strcmp(arg_ptr[i], "--version") == 0) {
        version_flag = 1;
      } else if (strcmp(arg_ptr[i], "-h") == 0 ||
                 strcmp(arg_ptr[i], "--help") == 0) {
        help_flag = 1;
      } else if (strcmp(arg_ptr[i], "-rec") == 0) {
	// do nothing (only to show that this is a valid argument)
      } else if (strcmp(arg_ptr[i-1], "-rec") == 0) {
        Reconstruct_Type_ptr = arg_ptr[i];
        strcpy(Reconstruct_Type, Reconstruct_Type_ptr);
      } else if (strcmp(arg_ptr[i], "-f") == 0) {
        file_flag = 1;
      } else if (strcmp(arg_ptr[i-1], "-f") == 0) {
        Input_File_Name_ptr = arg_ptr[i];
      } else if (strcmp(arg_ptr[i], "-t") == 0) {
	test_flag=1;
	if (num_arg-1>i){
	  TestSuit= arg_ptr[i+1];
	}
	if (num_arg-1>i+1){
	  TestNumber = atoi(arg_ptr[i+2]);
	}
	break;
      } else {
        error_flag = 1;
      } /* endif */
    } /* endfor */
  } /* endif */

  /* Display command line argument error message and
     terminate the program as required. */

  if (error_flag) {
    cout << "\nRECONSTRUCT ERROR: Invalid command line argument.\n";
    return (error_flag);
  } /* endif */

  /******************************************************************
   * DISPLAY THE PROGRAM TITLE AND VERSION INFORMATION AS REGUIRED. *
   ******************************************************************/

  if (version_flag || help_flag) {
    cout << '\n' << program_title_ptr << '\n'
	 << program_version_ptr << SourceCode::LastCommitted_Date() << '\n'
	 << " --> repository version: rev. " << SourceCode::LastCommitted_Revision() << '\n'
	 << " --> local version: rev. " << SourceCode::RevisionAtCompileTime() << '\n'
	 << " --> compiled on: " << SourceCode::TimeAtCompilation() << '\n';
    cout.flush();
    if (version_flag) return (0);
  } /* endif */

  /******************************************************************
   * DISPLAY THE PROGRAM HELP INFORMATION AS REQUIRED.              *
   ******************************************************************/

  if (help_flag) {
    cout << "Usage:\n";
    cout << "RECONSTRUCT [-v] [-h] [-i] [-b] [-rec type] [-f name]\n";
    cout << " -v (--version)  displays version information\n";
    cout << " -h (--help)     shows this help\n";
    cout << " -rec type       reconstructs `type' reconstruction\n" 
	 << "                 Possible types: '1D', '2D' or '3D'"
	 << " \n                (`1D' is default)\n";
    cout << " -f name         uses `name' input data file "
	 << "(`reconstruct_DEFAULT.in' is default)\n";
    cout.flush();
    return (0);
  } /* endif */

  /**********************************************************
   * UNIT TESTING                                           *
   *********************************************************/
  if (test_flag) {
    Test_Run(TestSuit,TestNumber);
    exit(0);
  }

  /***********************************************************  
   * PERFORM REQUIRED CALCULATIONS.                          *
   ***********************************************************/
  /* For the function of interest, solve the reconstruction 1D, 2D or 3D. */
  if (strcmp(Reconstruct_Type, "3D") == 0) {
    /* Reconstruction3D */
    error_flag = Reconstruction3DSolver(Input_File_Name_ptr);

  }else if (strcmp(Reconstruct_Type, "2D") == 0) {
    /* Reconstruction2D */
    error_flag = Reconstruction2DSolver(Input_File_Name_ptr);

  }else if (strcmp(Reconstruct_Type, "1D") == 0) {
    /* Reconstruction1D */
    error_flag = Reconstruction1DSolver(Input_File_Name_ptr);

  } else {
    /* try Reconstruction1D */
    std::cout << "Running Reconstruction1D ... " << std::endl;
    error_flag = Reconstruction1DSolver(Input_File_Name_ptr);
  }/* endif */

  if (error_flag) {
    return (error_flag);
  } /* endif */

  /********************************************************  
   * TERMINATE PROGRAM EXECUTION                          *
   ********************************************************/

  cout << "\n\nRECONSTRUCT: Execution complete.\n";
  return (0);

/* End RECONSTRUCT program. */

}
