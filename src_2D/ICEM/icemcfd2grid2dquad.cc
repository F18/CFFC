/*******************************************************************
 *******************************************************************
 **********                                               **********
 ******                   icemcfd2grid2dquad                 *******
 ****                                                           ****
 ****               Computer program for converting             ****
 ****    ICEM CFD mesh to CFDkit+caboodle Grid2D_Quad_Block     ****
 ****        format for subsequent reading from input file      ****
 ****                                                           ****
 ****                  Written by Ming Yao Ding                 ****
 ****              email: mingyao.ding@utoronto.ca              ****
 ****                  Modified by Stefan Neata                 ****
 ****              email: stefan.neata@utoronto.ca              ****
 ******       UTIAS, CFD & Propulsion Group, 1999-2007        ******
 **********                                               **********
 *******************************************************************
 *******************************************************************/

#include "ICEMCFD.h"

main(int num_arg, char *arg_ptr[]){

  /********************************************************  
   * VARIABLE DECLARATIONS                                *
   ********************************************************/

  // Name of program entered on the command line:
  char *command_name_ptr;

  // Title of code:
  char *program_title_ptr = 
     "icemcfd2grid2dquad: Converts and writes ICEM CFD Mesh to Grid2D_Quad_Block Format.";

  // Version of code:
  char *program_version_ptr = 
     "Version 1.30, UTIAS CFD & Propulsion Group, 1999-2007.";
  
  // Output file name:
  char *Output_File_Name_ptr = "icemcfd_mesh.grid";

  // Grid:
  Grid2D_Quad_Block          **Grid_ptr_Normal;

  // Grid type indicator:
  char *Grid_Type_ptr = "Normal";
  char Grid_Type[256];

  // Output file:
  ofstream output_file;

  // Command line flags:  
  int version_flag, 
      help_flag, 
      file_flag, 
      error_flag;

  // Other local variables:
  int i, nblocks_idir, nblocks_jdir;
  char **icemcfd_file_names;
  
  /********************************************************  
   * PARSE COMMAND LINE ARGUMENTS                         *
   ********************************************************/

  /* Initialize command line flags. */

  version_flag = 0;
  help_flag = 0;
  file_flag = 0;
  error_flag = 0;

  /* Save the command line name of the program for future use. */

  command_name_ptr = arg_ptr[0];
  
  /* Parse and interpret command line arguments.  Note that there
     are several different possible arguments which are:

     1) -v  lists program version information to standard output,
     2) -h  lists all possible optional arguments for program,
     3) -f name  uses "name" as the output data file rather than
                 the standard output data file "icemcfd_mesh.grid". */

  if (num_arg >= 2) {
    for (i = 1; i <= num_arg - 1; ++i) {
      if (strcmp(arg_ptr[i], "-v") == 0||
          strcmp(arg_ptr[i], "--version") == 0) {
        version_flag = 1;
      } else if (strcmp(arg_ptr[i], "-h") == 0 ||
                 strcmp(arg_ptr[i], "--help") == 0) {
        help_flag = 1;
      } else if (strcmp(arg_ptr[i], "-f") == 0) {
        file_flag = 1;
      } else if (strcmp(arg_ptr[i-1], "-f") == 0) {
        Output_File_Name_ptr = arg_ptr[i];
      } else {
        error_flag = 1;
      } /* endif */
    } /* endfor */
  } /* endif */
  strcpy(Grid_Type, Grid_Type_ptr);

  /* Display command line argument error message and
     terminate the program as required. */

  if (error_flag) {
    cout << "\nicemcfd2grid2dquad ERROR: Invalid command line argument.\n";
    return (error_flag);
  } /* endif */

  /******************************************************************
   * DISPLAY THE PROGRAM TITLE AND VERSION INFORMATION AS REGUIRED. *
   ******************************************************************/

  cout << '\n' << program_title_ptr << '\n';
  cout << program_version_ptr << '\n';
  cout << "Built using " << CFDkit_Version() << "\n";
  cout << ICEMCFD_Version() << "\n";
  cout.flush();
  if (version_flag) return (0);

  /******************************************************************
   * DISPLAY THE PROGRAM HELP INFORMATION AS REQUIRED.              *
   ******************************************************************/

  if (help_flag) {
     cout << "Usage:\n";
     cout << "icemcfd2grid2dquad [-v] [-h] [-f name]\n";
     cout << " -v (--version)  display version information\n";
     cout << " -h (--help)  show this help\n";
     cout << " -grid  sets the grid type (Normal is default)\n";
     cout << " -f name  use `name' output data file (`icemcfd_mesh.grid' is default)\n";
     cout.flush();
     return (0);
  } /* endif */

  /***********************************************************  
   * PROCESS ICEM CFD FILE.                                  *
   ***********************************************************/

  // Set default ICEM CFD topology, and family boco and topo files.
  icemcfd_file_names = ICEMCFD_get_filenames();
  
  // Read the ICEM CFD mesh and save it as a multiblock quadrilateral 2D mesh. 
  Grid_ptr_Normal = ICEMCFD_Read(icemcfd_file_names,
				 Grid_ptr_Normal,
				 &nblocks_idir,
				 &nblocks_jdir);
  if (Grid_ptr_Normal == NULL) {
    cout << "\n ICEMCFD ERROR: Unable to read ICEM CFD domain, boco, and/or topo file(s).\n";
    cout.flush();
    return (0);
  } /* endif */

  // Output the mesh to a data file in a format suitable for reading by
  // CFDkit+caboodle software.
  output_file.open(Output_File_Name_ptr,ios::out);

  output_file << "ICEMCFD" << "\n" 
              << GRID_ICEMCFD << "\n";
  Write_Multi_Block_Grid(Grid_ptr_Normal,
			 nblocks_idir,
			 nblocks_jdir,
			 output_file);

  output_file.close();
  
  /********************************************************  
   * TERMINATE PROGRAM EXECUTION                          *
   ********************************************************/

  cout << "\nicemcfd2grid2dquad: Execution complete.\n";
  return (0);

/* End icemcfd2grid2dquad program. */

}
