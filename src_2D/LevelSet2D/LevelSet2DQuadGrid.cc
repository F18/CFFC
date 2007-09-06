/**********************************************************************
 * LevelSet2DQuadGrid.cc                                              *
 *                                                                    *
 * Grid manipulation subroutines for 2D Level Set multi-block         *
 * quadrilateral mesh solution classes.                               *
 *                                                                    *
 **********************************************************************/

// Include 2D LevelSet quadrilateral mesh solution header file.

#ifndef _LEVELSET2D_QUAD_INCLUDED
#include "LevelSet2DQuad.h"
#endif // _LEVELSET2D_QUAD_INCLUDED

/**********************************************************************
 * LevelSet2D_Quad_Block -- Multiple Block External Subroutines for   *
 *                          Mesh.                                     *
 **********************************************************************/

/**********************************************************************
 * Routine: Multi_Block_Grid                                          *
 *                                                                    *
 * Generates multi-block quadilateral mesh.                           *
 *                                                                    *
 **********************************************************************/
Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
				     LevelSet2D_Input_Parameters &Input_Parameters) {

  
  // Generate appropriate mesh.
  switch(Input_Parameters.i_Grid) {
  case GRID_READ_FROM_DEFINITION_FILE :
    //Grid_ptr = Read_Multi_Block_Grid_Definition(Grid_ptr,
    //					Input_Parameters);
    if (Grid_ptr == NULL)
      cout << "\n " << CFFC_Name() 
	   << " LevelSet2D ERROR: Unable to open multi-block mesh definition file "
	   << Input_Parameters.Grid_Definition_File_Name << ".\n";
    break;

  case GRID_READ_FROM_GRID_DATA_FILE :
    //Grid_ptr = Read_Multi_Block_Grid(Grid_ptr,
    //			     Input_Parameters);
    if (Grid_ptr == NULL)
      cout << "\n " << CFFC_Name() 
	   << " LevelSet2D ERROR: Unable to open multi-block mesh data file "
	   << Input_Parameters.Grid_File_Name << ".\n";
    break;

  case GRID_SQUARE :
    Grid_ptr = Grid_Rectangular_Box(Grid_ptr,
				    Input_Parameters.Number_of_Blocks_Idir,
				    Input_Parameters.Number_of_Blocks_Jdir,
				    Input_Parameters.Box_Width,
				    Input_Parameters.Box_Width,
				    Input_Parameters.Number_of_Cells_Idir,
				    Input_Parameters.Number_of_Cells_Jdir,
				    Input_Parameters.Number_of_Ghost_Cells);
    break;

  case GRID_RECTANGULAR_BOX :
    Grid_ptr = Grid_Rectangular_Box(Grid_ptr,
				    Input_Parameters.Number_of_Blocks_Idir,
				    Input_Parameters.Number_of_Blocks_Jdir,
				    Input_Parameters.Box_Width,
				    Input_Parameters.Box_Height,
				    Input_Parameters.Number_of_Cells_Idir,
				    Input_Parameters.Number_of_Cells_Jdir,
				    Input_Parameters.Number_of_Ghost_Cells);
    break;

  case GRID_ROCKET_MOTOR :
    Input_Parameters.Box_Width = Input_Parameters.Grain_Length + Input_Parameters.Nozzle_Length;
    Input_Parameters.Box_Height = Input_Parameters.Grain_Radius;
    Grid_ptr = Grid_Rectangular_Box(Grid_ptr,
				    Input_Parameters.Number_of_Blocks_Idir,
				    Input_Parameters.Number_of_Blocks_Jdir,
				    Input_Parameters.Box_Width,
				    Input_Parameters.Box_Height,
				    Input_Parameters.Number_of_Cells_Idir,
				    Input_Parameters.Number_of_Cells_Jdir,
				    Input_Parameters.Number_of_Ghost_Cells);
    Input_Parameters.X_Shift = Vector2D(Input_Parameters.Nozzle_Length - HALF*(Input_Parameters.Box_Width),
					HALF*Input_Parameters.Box_Height);
    break;

  default:
    Grid_ptr = Grid_Rectangular_Box(Grid_ptr,
				    Input_Parameters.Number_of_Blocks_Idir,
				    Input_Parameters.Number_of_Blocks_Jdir,
				    Input_Parameters.Box_Width,
				    Input_Parameters.Box_Height,
				    Input_Parameters.Number_of_Cells_Idir,
				    Input_Parameters.Number_of_Cells_Jdir,
				    Input_Parameters.Number_of_Ghost_Cells);
    break;

  }

  // Set all non-BC_NONE boundaries to the specified BC type 
  // (BC_CONSTANT_EXTRAPOLATION or BC_LINEAR_EXTRAPOLATION).
  for (int jBlk = 0; jBlk < Input_Parameters.Number_of_Blocks_Jdir; jBlk++) {
    for (int iBlk = 0; iBlk < Input_Parameters.Number_of_Blocks_Idir; iBlk++) {
      if (jBlk == Input_Parameters.Number_of_Blocks_Jdir-1)
	Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(Input_Parameters.i_BC_Type);
      if (jBlk == 0)
	Grid_ptr[iBlk][jBlk].BndSouthSpline.setBCtype(Input_Parameters.i_BC_Type);
      if (iBlk == Input_Parameters.Number_of_Blocks_Idir-1)
	Grid_ptr[iBlk][jBlk].BndEastSpline.setBCtype(Input_Parameters.i_BC_Type);
      if (iBlk == 0)
	Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(Input_Parameters.i_BC_Type);
      Set_BCs(Grid_ptr[iBlk][jBlk]);
      Update_Exterior_Nodes(Grid_ptr[iBlk][jBlk]);
      Update_Cells(Grid_ptr[iBlk][jBlk]);
    }
  }

  // First translate quadrilateral mesh as specified by input parameters.
  if (abs(Input_Parameters.X_Shift) > TOLER)
    Translate_Multi_Block_Grid(Grid_ptr,
			       Input_Parameters.Number_of_Blocks_Idir,
			       Input_Parameters.Number_of_Blocks_Jdir,
			       Input_Parameters.X_Shift);
  
  // Next scale quadrilateral mesh as specified by input parameters.
  if (fabs(Input_Parameters.X_Scale-ONE) > TOLER)
    Scale_Multi_Block_Grid(Grid_ptr,
			   Input_Parameters.Number_of_Blocks_Idir,
			   Input_Parameters.Number_of_Blocks_Jdir,
			   Input_Parameters.X_Scale);
  
  // Finally rotate quadrilateral mesh as specified by input parameters.
  if (fabs(Input_Parameters.X_Rotate) > TOLER)
    Rotate_Multi_Block_Grid(Grid_ptr,
			    Input_Parameters.Number_of_Blocks_Idir,
			    Input_Parameters.Number_of_Blocks_Jdir,
			    TWO*PI*Input_Parameters.X_Rotate/360.00);
  
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
                                               LevelSet2D_Input_Parameters &Input_Parameters) {

  // Broadcast the multi-block quadilateral mesh.

#ifdef _MPI_VERSION
  MPI::COMM_WORLD.Bcast(&(Input_Parameters.Number_of_Blocks_Idir),
			1,
			MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&(Input_Parameters.Number_of_Blocks_Jdir),
			1,
			MPI::INT,0);
  Grid_ptr = Broadcast_Multi_Block_Grid(Grid_ptr,
					Input_Parameters.Number_of_Blocks_Idir,
					Input_Parameters.Number_of_Blocks_Jdir);
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
                                      LevelSet2D_Input_Parameters &Input_Parameters) {

  char *mesh_definition_file_name_ptr;
  ofstream mesh_definition_file;
  
  // Open the grid definition file.
  mesh_definition_file_name_ptr = Input_Parameters.Grid_Definition_File_Name;
  mesh_definition_file.open(mesh_definition_file_name_ptr,ios::out);
  if (mesh_definition_file.bad()) return 1;
  
  // Write grid type information.
  mesh_definition_file << Input_Parameters.Grid_Type << "\n" 
		       << Input_Parameters.i_Grid << "\n";
  
  // Write the grid definition information for each quadrilateral grid
  // block.
  Write_Multi_Block_Grid_Definition(Grid_ptr,
				    Input_Parameters.Number_of_Blocks_Idir,
				    Input_Parameters.Number_of_Blocks_Jdir,
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
                                                     LevelSet2D_Input_Parameters &Input_Parameters) {

  char buffer[256];
  char *mesh_definition_file_name_ptr;
  ifstream mesh_definition_file;
  
  // Open the grid definition file.
  mesh_definition_file_name_ptr = Input_Parameters.Grid_Definition_File_Name;
  mesh_definition_file.open(mesh_definition_file_name_ptr,ios::in);
  if (mesh_definition_file.bad()) return NULL;
  
  // Read grid type information.
  mesh_definition_file.getline(buffer,sizeof(buffer));
  strcpy(Input_Parameters.Grid_Type,buffer);
  mesh_definition_file >> Input_Parameters.i_Grid;

  // Read the grid definition information for each quadrilateral grid 
  // block.
  Grid_ptr = Read_Multi_Block_Grid_Definition(Grid_ptr,
					      Input_Parameters.Number_of_Blocks_Idir,
					      Input_Parameters.Number_of_Blocks_Jdir,
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
                           LevelSet2D_Input_Parameters &Input_Parameters) {

  char *mesh_file_name_ptr;
  ofstream mesh_file;
  
  // Open the grid data output file.
  mesh_file_name_ptr = Input_Parameters.Grid_File_Name;
  mesh_file.open(mesh_file_name_ptr,ios::out);
  if (mesh_file.bad()) return 1;
  
  // Write grid type information.
  mesh_file << Input_Parameters.Grid_Type << "\n" 
	    << Input_Parameters.i_Grid << "\n";
  
  // Write the grid data for each quadrilateral grid block.
  Write_Multi_Block_Grid(Grid_ptr,
			 Input_Parameters.Number_of_Blocks_Idir,
			 Input_Parameters.Number_of_Blocks_Jdir,
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
						   LevelSet2D_Input_Parameters &Input_Parameters) {

  char buffer[256];
  char *mesh_file_name_ptr;
  ifstream mesh_file;
  
  // Open the grid data input file.
  mesh_file_name_ptr = Input_Parameters.Grid_File_Name;
  mesh_file.open(mesh_file_name_ptr,ios::in);
  if (mesh_file.bad()) return NULL;
  
  // Read grid type information.
  mesh_file.getline(buffer,sizeof(buffer));
  strcpy(Input_Parameters.Grid_Type,buffer);
  mesh_file >> Input_Parameters.i_Grid;

  // Read the grid data for each quadrilateral grid block.
  Grid_ptr = Read_Multi_Block_Grid(Grid_ptr,
				   Input_Parameters.Number_of_Blocks_Idir,
				   Input_Parameters.Number_of_Blocks_Jdir,
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
 * suitable for plotting with TECPLOT.  Returns a non-zero value if   *
 * unable to write the TECPLOT file.                                  *
 *                                                                    *
 **********************************************************************/
int Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                   LevelSet2D_Input_Parameters &Input_Parameters) {

  int i, j;
  char prefix[256], extension[256], mesh_file_name[256];
  char *mesh_file_name_ptr;
  ofstream mesh_file;  
  
  // Determine prefix of grid data output file name.
  i = 0;
  while (1) {
    if (Input_Parameters.Output_File_Name[i] == ' ' ||
	Input_Parameters.Output_File_Name[i] == '.') break;
    prefix[i]=Input_Parameters.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(Input_Parameters.Output_File_Name)) break;
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
		 Input_Parameters.Number_of_Blocks_Idir,
		 Input_Parameters.Number_of_Blocks_Jdir,
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
                         LevelSet2D_Input_Parameters &Input_Parameters) {

  int i, j;
  char prefix[256], extension[256], mesh_file_name[256];
  char *mesh_file_name_ptr;
  ofstream mesh_file;  
  
  // Determine prefix of grid data output file name.
  i = 0;
  while (1) {
    if (Input_Parameters.Output_File_Name[i] == ' ' ||
	Input_Parameters.Output_File_Name[i] == '.') break;
    prefix[i]=Input_Parameters.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(Input_Parameters.Output_File_Name)) break;
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
		       Input_Parameters.Number_of_Blocks_Idir,
		       Input_Parameters.Number_of_Blocks_Jdir,
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
                         LevelSet2D_Input_Parameters &Input_Parameters) {

  int i, j;
  char prefix[256], extension[256], mesh_file_name[256];
  char *mesh_file_name_ptr;
  ofstream mesh_file;  

  // Determine prefix of grid data output file name.
  i = 0;
  while (1) {
    if (Input_Parameters.Output_File_Name[i] == ' ' ||
	Input_Parameters.Output_File_Name[i] == '.') break;
    prefix[i]=Input_Parameters.Output_File_Name[i];
    i = i + 1;
    if (i > strlen(Input_Parameters.Output_File_Name)) break;
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
		       Input_Parameters.Number_of_Blocks_Idir,
		       Input_Parameters.Number_of_Blocks_Jdir,
		       mesh_file);
  
  // Close the grid data output file.
  mesh_file.close();
  
  // Writing of grid data output file complete.  Return zero value.
  return 0;
    
}







