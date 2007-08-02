/*******************************************************************
 *******************************************************************
 **********                                               **********
 ******                   icemcfd2tecplot                    *******
 ****                                                           ****
 ****               Computer program for converting             ****
 ****    ICEM CFD mesh to TECPLOT output format for plotting    ****
 ****                                                           ****
 ****                  Written by Ming Yao Ding                 ****
 ****              email: mingyao.ding@utoronto.ca              ****
 ****             Modified & Written by Stefan Neata            ****
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
     "icemcfd2tecplot: Converts and writes ICEM CFD Mesh to Tecplot Format.";

  // Version of code:
  char *program_version_ptr = 
     "Version 1.30, UTIAS CFD & Propulsion Group, 1999-2007.";
  
  // Output file name:
  char *Output_File_Name_ptr = "icemcfd_mesh.dat";

  // 2D and 3D Grid data structures:
  Grid2D_Quad_Block **Grid2D_ptr;
  Grid3D_Hexa_Block ***Grid3D_ptr;

  // Output file:
  ofstream output_file;

  // Command line flags:  
  int version_flag, 
      help_flag, 
      file_flag, 
      error_flag,
      node_flag,
      three_D_flag;

  // Other local variables:
  int i, nblocks_idir, nblocks_jdir, nblocks_kdir;
  char **icemcfd_file_names;

  /********************************************************  
   * PARSE COMMAND LINE ARGUMENTS                         *
   ********************************************************/

  /* Initialize command line flags. */

  version_flag = 0;
  help_flag = 0;
  file_flag = 0;
  error_flag = 0;
  node_flag = 1;
  three_D_flag = 0;

  /* Save the command line name of the program for future use. */

  command_name_ptr = arg_ptr[0];
  
  /* Parse and interpret command line arguments.  Note that there
     are several different possible arguments which are:

     1) -v  lists program version information to standard output,
     2) -h  lists all possible optional arguments for program,
     3) -f name  uses "name" as the output data file rather than
                 the standard output data file "icemcfd_mesh.dat", 
     4) -c outputs cell information to TECPLOT data file,
     5) -n outputs node information to TECPLOT data file,
     6) -2D instructs program that grid is a 2D mesh,
     7) -3D instructs program that grid is a 3D mesh. */

  if (num_arg >= 2) {
    for (i = 1; i <= num_arg - 1; ++i) {
      if (strcmp(arg_ptr[i], "-v") == 0 ||
          strcmp(arg_ptr[i], "--version") == 0) {
        version_flag = 1;
      } else if (strcmp(arg_ptr[i], "-h") == 0 ||
                 strcmp(arg_ptr[i], "--help") == 0) {
        help_flag = 1;
      } else if (strcmp(arg_ptr[i], "-f") == 0) {
        file_flag = 1;
      } else if (strcmp(arg_ptr[i-1], "-f") == 0) {
        Output_File_Name_ptr = arg_ptr[i];
      } else if (strcmp(arg_ptr[i], "-n") == 0 ||
                 strcmp(arg_ptr[i], "--nodes") == 0) {
        node_flag = 1;
      } else if (strcmp(arg_ptr[i], "-c") == 0 ||
                 strcmp(arg_ptr[i], "--cells") == 0) {
        node_flag = 0;
     } else if (strcmp(arg_ptr[i], "-2D") == 0 ||
                 strcmp(arg_ptr[i], "--2D") == 0) {
        three_D_flag = 0;
      } else if (strcmp(arg_ptr[i], "-3D") == 0 ||
                 strcmp(arg_ptr[i], "--3D") == 0) {
        three_D_flag = 1;
      } else {
        error_flag = 1;
      } /* endif */
    } /* endfor */
  } /* endif */

  /* Display command line argument error message and
     terminate the program as required. */

  if (error_flag) {
    cout << "\nicemcfd2tecplot ERROR: Invalid command line argument.\n";
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
     cout << "icemcfd2tecplot [-v] [-h] [-f name] [-c -n] [-2D -3D]\n";
     cout << " -v (--version)  display version information\n";
     cout << " -h (--help)  show this help\n";
     cout << " -f name  use `name' output data file (`icemcfd_mesh.dat' is default)\n";
     cout << " -n (--nodes)  output nodes of multiblock mesh (default)\n";
     cout << " -c (--cells)  output cells of multiblock mesh\n";  
     cout << " -2D (--2D) grid is a 2D multiblock mesh (default)\n";
     cout << " -3D (--3D) grid is a 3D multiblock mesh\n";  
     cout.flush();
     return (0);
  } /* endif */

  /***********************************************************  
   * PROCESS ICEM CFD FILE.                                  *
   ***********************************************************/

  // Set default ICEM CFD topology, and family boco and topo files.  
  icemcfd_file_names = ICEMCFD_get_filenames();

  if (three_D_flag) {
     // Read the ICEM CFD mesh and save it as a multiblock hexahedral 3D mesh.
     Grid3D_ptr = ICEMCFD_Read(icemcfd_file_names,
                               Grid3D_ptr,
                               &nblocks_idir,
                               &nblocks_jdir,
                               &nblocks_kdir);
     if (Grid3D_ptr == NULL) {
        cout << "\n ICEMCFD ERROR: Unable to read ICEM CFD domain, boco, and/or topo file(s).\n";
        cout.flush();
        return (0);
  } /* endif */
  } else {
     // Read the ICEM CFD mesh and save it as a multiblock quadrilateral 2D mesh.
     Grid2D_ptr = ICEMCFD_Read(icemcfd_file_names,
                               Grid2D_ptr,
                               &nblocks_idir,
                               &nblocks_jdir);
     if (Grid2D_ptr == NULL) {
        cout << "\n ICEMCFD ERROR: Unable to read ICEM CFD domain, boco, and/or topo file(s).\n";
        cout.flush();
        return (0);
     } /* endif */
  } /* endif */

  // Output the mesh to a data file in a format suitable for plotting with Tecplot.
  output_file.open(Output_File_Name_ptr,ios::out);

  if (three_D_flag) {
     if (node_flag) {
        Output_Tecplot(Grid3D_ptr,
                       nblocks_idir,
                       nblocks_jdir,
		       nblocks_kdir,
                       output_file);
  } else {
        Output_Cells_Tecplot(Grid3D_ptr,
                             nblocks_idir,
                             nblocks_jdir,
                             nblocks_kdir,
                             output_file);
  } /* endif */
  } else {
     if (node_flag) {
        Output_Tecplot(Grid2D_ptr,
                       nblocks_idir,
                       nblocks_jdir,
                       output_file);
     } else {
        Output_Cells_Tecplot(Grid2D_ptr,
                             nblocks_idir,
                             nblocks_jdir,
                             output_file);
     } /* endif */
  } /* endif */

  output_file.close();
  
  /********************************************************  
   * TERMINATE PROGRAM EXECUTION                          *
   ********************************************************/

  cout << "\nicemcfd2tecplot: Execution complete.\n";
  return (0);

/* End icemcfd2tecplot program. */

}
