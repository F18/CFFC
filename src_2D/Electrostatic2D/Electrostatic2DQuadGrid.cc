/**********************************************************************
 * Electrostatic2DQuadGrid.cc: Grid manipulation subroutines for 2D   *
 *                             Electrostatic multi-block              *
 *                             quadrilateral mesh solution classes.   *
 **********************************************************************/

// Include 2D Electrostatic quadrilateral mesh solution header file.

#ifndef _ELECTROSTATIC2D_QUAD_INCLUDED
#include "Electrostatic2DQuad.h"
#endif // _ELECTROSTATIC2D_QUAD_INCLUDED

/**********************************************************************
 * Electrostatic2D_Quad_Block -- Multiple Block External Subroutines  *
 *                               for Quadrilateral Mesh.              *
 **********************************************************************/

/**********************************************************************
 * Routine: Multi_Block_Grid                                          *
 *                                                                    *
 * Generates multi-block quadilateral mesh.                           *
 *                                                                    *
 **********************************************************************/
Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
				     Electrostatic2D_Input_Parameters &IP) {

  // Generate appropriate mesh.
  switch(IP.i_Grid) {
  case GRID_READ_FROM_DEFINITION_FILE :
    Grid_ptr = Read_Multi_Block_Grid_Definition(Grid_ptr,IP);
    if (Grid_ptr == NULL) {
      cout << "\n " << CFDkit_Name() 
	   << " Electrostatic2D ERROR: Unable to open multi-block mesh definition file "
	   << IP.Grid_Definition_File_Name << ".\n";
    }
    break;
  case GRID_READ_FROM_GRID_DATA_FILE :
    Grid_ptr = Read_Multi_Block_Grid(Grid_ptr,IP);
    if (Grid_ptr == NULL) {
      cout << "\n " << CFDkit_Name() 
	   << " Electrostatic2D ERROR: Unable to open multi-block mesh data file "
	   << IP.Grid_File_Name << ".\n";
    }
    break;
  case GRID_SQUARE :
    Grid_ptr = Grid_Rectangular_Box(Grid_ptr,
				    IP.Number_of_Blocks_Idir,
				    IP.Number_of_Blocks_Jdir,
				    IP.Box_Width,
  				    IP.Box_Width,
				    IP.Number_of_Cells_Idir,
				    IP.Number_of_Cells_Jdir,
				    IP.Number_of_Ghost_Cells);
    break;
  case GRID_RECTANGULAR_BOX :
    Grid_ptr = Grid_Rectangular_Box(Grid_ptr,
				    IP.Number_of_Blocks_Idir,
				    IP.Number_of_Blocks_Jdir,
				    IP.Box_Width,
				    IP.Box_Height,
				    IP.Number_of_Cells_Idir,
				    IP.Number_of_Cells_Jdir,
				    IP.Number_of_Ghost_Cells);
    break;
  case GRID_CIRCULAR_CYLINDER :
    Grid_ptr = Grid_Circular_Cylinder(Grid_ptr,
				      IP.Number_of_Blocks_Idir,
				      IP.Number_of_Blocks_Jdir,
				      IP.Cylinder_Radius,
				      IP.Mesh_Stretching_Type_Idir,
				      IP.Mesh_Stretching_Type_Jdir,
				      IP.Mesh_Stretching_Factor_Idir,
				      IP.Mesh_Stretching_Factor_Jdir,
				      IP.Number_of_Cells_Idir,
				      IP.Number_of_Cells_Jdir,
				      IP.Number_of_Ghost_Cells);
    break;
  case GRID_FREE_JET :
    Grid_ptr = Grid_Free_Jet(Grid_ptr,
			     IP.Number_of_Blocks_Idir,
			     IP.Number_of_Blocks_Jdir,
			     IP.Orifice_Radius,
			     IP.Number_of_Cells_Idir,
			     IP.Number_of_Cells_Jdir,
			     IP.Number_of_Ghost_Cells);
    break;
  case GRID_DESOLVATION_CHAMBER :
    Grid_ptr = Grid_Desolvation_Chamber(Grid_ptr,
					BC_FIXED,
					IP.Number_of_Blocks_Idir,
					IP.Number_of_Blocks_Jdir,
					IP.Number_of_Cells_Idir,
					IP.Number_of_Cells_Jdir,
					IP.Number_of_Ghost_Cells);
    for (int jBlk = 0; jBlk < IP.Number_of_Blocks_Jdir; jBlk++) {
      for (int iBlk = 0; iBlk < IP.Number_of_Blocks_Idir; iBlk++) {
	if (iBlk == 0 && jBlk == 0) {
	  Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(BC_FIXED);
	} else if (iBlk == 4 && jBlk == 0) {
	  Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(BC_FIXED);
	} else if (iBlk == 0 && jBlk == 1) {
	  Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(BC_FIXED);
	} else if (iBlk == 1 && jBlk == 1) {
	  Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(BC_FIXED);
	} else if (iBlk == 2 && jBlk == 1) {
	  Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(BC_FIXED);
	} else if (iBlk == 3 && jBlk == 1) {
	  Grid_ptr[iBlk][jBlk].BndEastSpline.setBCtype(BC_FIXED);
	} else if (iBlk == 5 && jBlk == 1) {
	  Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(BC_FIXED);
	} else if (iBlk == 0 && jBlk == 2) {
	  Grid_ptr[iBlk][jBlk].BndEastSpline.setBCtype(BC_FIXED);
	  Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(BC_FIXED);
	} else if (iBlk == 3 && jBlk == 2) {
	  Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(BC_FIXED);
	  Grid_ptr[iBlk][jBlk].BndEastSpline.setBCtype(BC_FIXED);
	  Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(BC_FIXED);
	} else if (iBlk == 5 && jBlk == 2) {
	  Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(BC_FIXED);
	} else if (iBlk == 0 && jBlk == 3) {
	  Grid_ptr[iBlk][jBlk].BndEastSpline.setBCtype(BC_FIXED);
	  Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(BC_FIXED);
	} else if (iBlk == 5 && jBlk == 3) {
	  Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(BC_FIXED);
	} else if (iBlk == 0 && jBlk == 4) {
	  Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(BC_FIXED);
	  Grid_ptr[iBlk][jBlk].BndEastSpline.setBCtype(BC_FIXED);
	  Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(BC_FIXED);
	} else if (iBlk == 5 && jBlk == 4) {
	  Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(BC_FIXED);
	} else if (iBlk == 5 && jBlk == 5) {
	  Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(BC_FIXED);
	  Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(BC_FIXED);
	} else if (iBlk == 6 && jBlk == 5) {
	  Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(BC_FIXED);
	}
	if (Grid_ptr[iBlk][jBlk].Node != NULL) {
	  Set_BCs(Grid_ptr[iBlk][jBlk]);
	  Update_Exterior_Nodes(Grid_ptr[iBlk][jBlk]);
	  Update_Cells(Grid_ptr[iBlk][jBlk]);
	}
      }
    }
    break;
  case GRID_ICEMCFD :
    Grid_ptr = ICEMCFD_Read(IP.ICEMCFD_FileNames,
    			    Grid_ptr,
    			    &IP.Number_of_Blocks_Idir,
    			    &IP.Number_of_Blocks_Jdir);
    break;
  default:
    Grid_ptr = Grid_Rectangular_Box(Grid_ptr,
				    IP.Number_of_Blocks_Idir,
				    IP.Number_of_Blocks_Jdir,
				    IP.Box_Width,
				    IP.Box_Height,
				    IP.Number_of_Cells_Idir,
				    IP.Number_of_Cells_Jdir,
				    IP.Number_of_Ghost_Cells);
    break;
  }

  // Reset boundary conditions if required.
  if (IP.BCs_Specified) {
    for (int jBlk = 0; jBlk < IP.Number_of_Blocks_Jdir; jBlk++) {
      for (int iBlk = 0; iBlk < IP.Number_of_Blocks_Idir; iBlk++) {
	if (jBlk == IP.Number_of_Blocks_Jdir-1)
	  Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(IP.BC_North);
	if (jBlk == 0)
	  Grid_ptr[iBlk][jBlk].BndSouthSpline.setBCtype(IP.BC_South);
	if (iBlk == IP.Number_of_Blocks_Idir-1)
	  Grid_ptr[iBlk][jBlk].BndEastSpline.setBCtype(IP.BC_East);
	if (iBlk == 0)
	  Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(IP.BC_West);
	Set_BCs(Grid_ptr[iBlk][jBlk]);
	Update_Exterior_Nodes(Grid_ptr[iBlk][jBlk]);
	Update_Cells(Grid_ptr[iBlk][jBlk]);
      }
    }
  }

  // First translate quadrilateral mesh as specified by input parameters.
  if (abs(IP.X_Shift) > TOLER)
    Translate_Multi_Block_Grid(Grid_ptr,
			       IP.Number_of_Blocks_Idir,
			       IP.Number_of_Blocks_Jdir,
			       IP.X_Shift);

  // Next scale quadrilateral mesh as specified by input parameters.
  if (fabs(IP.X_Scale-ONE) > TOLER)
    Scale_Multi_Block_Grid(Grid_ptr,
			   IP.Number_of_Blocks_Idir,
			   IP.Number_of_Blocks_Jdir,
			   IP.X_Scale);

  // Finally rotate quadrilateral mesh as specified by input parameters.
  if (fabs(IP.X_Rotate) > TOLER)
    Rotate_Multi_Block_Grid(Grid_ptr,
			    IP.Number_of_Blocks_Idir,
			    IP.Number_of_Blocks_Jdir,
			    TWO*PI*IP.X_Rotate/360.00);

  // Return the multi-block quadilateral mesh.
  return Grid_ptr;

}

/**********************************************************************
 * Routine: Broadcast_Multi_Block_Grid                                *
 *                                                                    *
 * Broadcast multi-block quadilateral mesh to all processors involved *
 * in the calculation from the primary processor using the MPI        *
 * broadcast routine.                                                 *
 *                                                                    *
 **********************************************************************/
Grid2D_Quad_Block** Broadcast_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
					       Electrostatic2D_Input_Parameters &IP) {

  // Broadcast the multi-block quadilateral mesh.
#ifdef _MPI_VERSION
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Idir),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(IP.Number_of_Blocks_Jdir),
			1,
			MPI::INT,0);
  Grid_ptr = Broadcast_Multi_Block_Grid(Grid_ptr,
					IP.Number_of_Blocks_Idir,
					IP.Number_of_Blocks_Jdir);
#endif

  // Return the multi-block quadilateral mesh.
  return Grid_ptr;

}

/**********************************************************************
 * Routine: Write_Multi_Block_Grid_Definition                         *
 *                                                                    *
 * Writes a grid definition file for a multi-block quadilateral mesh  *
 * in a format suitable for retrieval and re-use purposes.  Returns a *
 * non-zero value if unable to write the grid definition file.        *
 *                                                                    *
 **********************************************************************/
int Write_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                      Electrostatic2D_Input_Parameters &IP) {

  char *mesh_definition_file_name_ptr;
  ofstream mesh_definition_file;

  // Open the grid definition file.
  mesh_definition_file_name_ptr = IP.Grid_Definition_File_Name;
  mesh_definition_file.open(mesh_definition_file_name_ptr,ios::out);
  if (mesh_definition_file.bad()) return 1;

  // Write grid type information.
  mesh_definition_file << IP.Grid_Type << endl
		       << IP.i_Grid << endl;

  // Write the grid definition information for each 
  // quadrilateral grid block.
  Write_Multi_Block_Grid_Definition(Grid_ptr,
				    IP.Number_of_Blocks_Idir,
				    IP.Number_of_Blocks_Jdir,
				    mesh_definition_file);

  // Close the grid definition file.
  mesh_definition_file.close();

  // Writing of grid definition file complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Read_Multi_Block_Grid_Definition                          *
 *                                                                    *
 * Reads a grid definition file for a multi-block quadilateral mesh.  *
 * Returns a pointer to the mesh.                                     *
 *                                                                    *
 **********************************************************************/
Grid2D_Quad_Block** Read_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
						     Electrostatic2D_Input_Parameters &IP) {

  char buffer[256];
  char *mesh_definition_file_name_ptr;
  ifstream mesh_definition_file;

  // Open the grid definition file.
  mesh_definition_file_name_ptr = IP.Grid_Definition_File_Name;
  mesh_definition_file.open(mesh_definition_file_name_ptr,ios::in);
  if (mesh_definition_file.bad()) return NULL;

  // Read grid type information.
  mesh_definition_file.getline(buffer,sizeof(buffer));
  strcpy(IP.Grid_Type,buffer);
  mesh_definition_file >> IP.i_Grid;

  // Read the grid definition information for each
  // quadrilateral grid block.
  Grid_ptr = Read_Multi_Block_Grid_Definition(Grid_ptr,
					      IP.Number_of_Blocks_Idir,
					      IP.Number_of_Blocks_Jdir,
					      mesh_definition_file);

  // Close the grid definition file.
  mesh_definition_file.close();

  // Reading of grid definition file complete.  Return the mesh.
  return Grid_ptr;

}

/**********************************************************************
 * Routine: Write_Multi_Block_Grid                                    *
 *                                                                    *
 * Writes multi-block quadilateral mesh to a grid data file in a      *
 * format suitable for retrieval and re-use purposes.  Returns a      *
 * non-zero value if unable to write the grid data file.              *
 *                                                                    *
 **********************************************************************/
int Write_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                           Electrostatic2D_Input_Parameters &IP) {

  char *mesh_file_name_ptr;
  ofstream mesh_file;

  // Open the grid data output file.
  mesh_file_name_ptr = IP.Grid_File_Name;
  mesh_file.open(mesh_file_name_ptr,ios::out);
  if (mesh_file.bad()) return 1;

  // Write grid type information.
  mesh_file << IP.Grid_Type << "\n" << IP.i_Grid << "\n";

  // Write the grid data for each quadrilateral grid block.
  Write_Multi_Block_Grid(Grid_ptr,
			 IP.Number_of_Blocks_Idir,
			 IP.Number_of_Blocks_Jdir,
			 mesh_file);

  // Close the grid data output file.
  mesh_file.close();

  // Writing of grid data output file complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Read_Multi_Block_Grid                                     *
 *                                                                    *
 * Reads multi-block quadilateral mesh from a grid data file.         *
 * Returns a pointer to the mesh.                                     *
 *                                                                    *
 **********************************************************************/
Grid2D_Quad_Block** Read_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
					  Electrostatic2D_Input_Parameters &IP) {

  char buffer[256];
  char *mesh_file_name_ptr;
  ifstream mesh_file;

  // Open the grid data input file.
  mesh_file_name_ptr = IP.Grid_File_Name;
  mesh_file.open(mesh_file_name_ptr,ios::in);
  if (mesh_file.bad()) return NULL;

  // Read grid type information.
  mesh_file.getline(buffer,sizeof(buffer));
  strcpy(IP.Grid_Type,buffer);
  mesh_file >> IP.i_Grid;

  // Read the grid data for each quadrilateral grid block.
  Grid_ptr = Read_Multi_Block_Grid(Grid_ptr,
				   IP.Number_of_Blocks_Idir,
				   IP.Number_of_Blocks_Jdir,
				   mesh_file);

  // Close the grid data input file.
  mesh_file.close();

  // Reading of grid data input file complete.  Return the mesh.
  return Grid_ptr;

}

/**********************************************************************
 * Routine: Output_Tecplot                                            *
 *                                                                    *
 * Writes the nodes of a multi-block quadilateral mesh in a format    *
 * suitable for plotting with TECPLOT.   Returns a non-zero value if  *
 * unable to write the TECPLOT file.                                  *
 *                                                                    *
 **********************************************************************/
int Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                   Electrostatic2D_Input_Parameters &IP) {

  int i;
  char prefix[256], extension[256], mesh_file_name[256];
  char *mesh_file_name_ptr;
  ofstream mesh_file;

  // Determine prefix of grid data output file name.
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i]=IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_mesh");

  // Determine grid data output file name.
  strcpy(extension,".dat");
  strcpy(mesh_file_name,prefix);
  strcat(mesh_file_name,extension);
  mesh_file_name_ptr = mesh_file_name;

  // Open the grid data output file.
  mesh_file.open(mesh_file_name_ptr,ios::out);
  if (mesh_file.bad()) return 1;

  // Write the node locations for each quadrilateral grid block.
  Output_Tecplot(Grid_ptr,
		 IP.Number_of_Blocks_Idir,
		 IP.Number_of_Blocks_Jdir,
		 mesh_file);

  // Close the grid data output file.
  mesh_file.close();

  // Writing of grid data output file complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Output_Nodes_Tecplot                                      *
 *                                                                    *
 * Writes the nodes of a multi-block quadilateral mesh in a format    *
 * suitable for plotting with TECPLOT.  Includes boundary nodes.      *
 * Returns a non-zero value if unable to write the TECPLOT file.      *
 *                                                                    *
 **********************************************************************/
int Output_Nodes_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                         Electrostatic2D_Input_Parameters &IP) {

  int i;
  char prefix[256], extension[256], mesh_file_name[256];
  char *mesh_file_name_ptr;
  ofstream mesh_file;

  // Determine prefix of grid data output file name.
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i]=IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_mesh_nodes");

  // Determine grid data output file name.
  strcpy(extension,".dat");
  strcpy(mesh_file_name,prefix);
  strcat(mesh_file_name,extension);
  mesh_file_name_ptr = mesh_file_name;

  // Open the grid data output file.
  mesh_file.open(mesh_file_name_ptr,ios::out);
  if (mesh_file.bad()) return 1;

  // Write the node locations for each quadrilateral grid block.
  Output_Nodes_Tecplot(Grid_ptr,
		       IP.Number_of_Blocks_Idir,
		       IP.Number_of_Blocks_Jdir,
		       mesh_file);

  // Close the grid data output file.
  mesh_file.close();

  // Writing of grid data output file complete.  Return zero value.
  return 0;

}

/**********************************************************************
 * Routine: Output_Cells_Tecplot                                      *
 *                                                                    *
 * Writes the cells of a multi-block quadilateral mesh in a format    *
 * suitable for plotting with TECPLOT.  Returns a non-zero value if   *
 * unable to write the TECPLOT file.                                  *
 *                                                                    *
 **********************************************************************/
int Output_Cells_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                         Electrostatic2D_Input_Parameters &IP) {

  int i;
  char prefix[256], extension[256], mesh_file_name[256];
  char *mesh_file_name_ptr;
  ofstream mesh_file;

  // Determine prefix of grid data output file name.
  i = 0;
  while (1) {
    if (IP.Output_File_Name[i] == ' ' ||
	IP.Output_File_Name[i] == '.') break;
    prefix[i]=IP.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(IP.Output_File_Name)) break;
  }
  prefix[i] = '\0';
  strcat(prefix,"_mesh_cells");

  // Determine grid data output file name.
  strcpy(extension,".dat");
  strcpy(mesh_file_name,prefix);
  strcat(mesh_file_name,extension);
  mesh_file_name_ptr = mesh_file_name;

  // Open the grid data output file.
  mesh_file.open(mesh_file_name_ptr,ios::out);
  if (mesh_file.bad()) return 1;

  // Write the node locations for each quadrilateral grid block.
  Output_Cells_Tecplot(Grid_ptr,
		       IP.Number_of_Blocks_Idir,
		       IP.Number_of_Blocks_Jdir,
		       mesh_file);

  // Close the grid data output file.
  mesh_file.close();

  // Writing of grid data output file complete.  Return zero value.
  return 0;

}
