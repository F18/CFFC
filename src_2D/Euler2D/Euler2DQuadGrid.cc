/* Euler2DQuadGrid.cc:  Grid Manipulation Subroutines for 2D Euler 
                        Multi-Block Quadrilateral Mesh 
                        Solution Classes. */

/* Include 2D Euler quadrilateral mesh solution header file. */

#ifndef _EULER2D_QUAD_INCLUDED
#include "Euler2DQuad.h"
#endif // _EULER2D_QUAD_INCLUDED

/**************************************************************************
 * Euler2D_Quad_Block -- Multiple Block External Subroutines for Mesh.    *
 **************************************************************************/

/********************************************************
 * Routine: Multi_Block_Grid                            *
 *                                                      *
 * Generates multi-block quadilateral mesh.             *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                     Euler2D_Input_Parameters &Input_Parameters) {

    int iBlk, jBlk;

    /* Generate appropriate mesh. */

    switch(Input_Parameters.i_Grid) {
      case GRID_READ_FROM_DEFINITION_FILE :
        Grid_ptr = Read_Multi_Block_Grid_Definition(Grid_ptr,
                                                    Input_Parameters);
        if (Grid_ptr == NULL) {
	   cout << "\n " << CFFC_Name() 
                << " Euler2D ERROR: Unable to open multi-block mesh definition file "
                << Input_Parameters.Grid_Definition_File_Name << ".\n";
        } /* endif */
        if (Input_Parameters.i_Grid == GRID_NASA_ROTOR_37) {
           Input_Parameters.NASA_Rotor37.init(Input_Parameters.Rotor_Flow_Type,
                                              Input_Parameters.NASA_Rotor37_Data_Directory);
           Input_Parameters.Wo.setgas("AIR");
 	   Input_Parameters.NASA_Rotor37.getPstateREL_up(Input_Parameters.Wo,Input_Parameters.Rotor_Percent_Span);
           Input_Parameters.Pressure = Input_Parameters.Wo.p;
           Input_Parameters.Temperature = Input_Parameters.Wo.T();
           Input_Parameters.Mach_Number = 
              Input_Parameters.NASA_Rotor37.getMachREL_up(Input_Parameters.Rotor_Percent_Span);
           Input_Parameters.Flow_Angle = atan2(Input_Parameters.Wo.v.y, Input_Parameters.Wo.v.x); 
           if (Input_Parameters.Flow_Angle < ZERO) 
              Input_Parameters.Flow_Angle = TWO*PI + Input_Parameters.Flow_Angle;
           Input_Parameters.Flow_Angle = 180.00*Input_Parameters.Flow_Angle/PI;
        } else if (Input_Parameters.i_Grid == GRID_NASA_ROTOR_67) {
           Input_Parameters.NASA_Rotor67.init(Input_Parameters.Rotor_Flow_Type,
                                              Input_Parameters.NASA_Rotor67_Data_Directory);
           Input_Parameters.Wo.setgas("AIR");
	   Input_Parameters.NASA_Rotor67.getPstateREL_up(Input_Parameters.Wo,Input_Parameters.Rotor_Percent_Span);
           Input_Parameters.Pressure = Input_Parameters.Wo.p;
           Input_Parameters.Temperature = Input_Parameters.Wo.T();
           Input_Parameters.Mach_Number = 
              Input_Parameters.NASA_Rotor67.getMachREL_up(Input_Parameters.Rotor_Percent_Span);
           Input_Parameters.Flow_Angle = atan2(Input_Parameters.Wo.v.y, Input_Parameters.Wo.v.x); 
           if (Input_Parameters.Flow_Angle < ZERO) 
              Input_Parameters.Flow_Angle = TWO*PI + Input_Parameters.Flow_Angle;
           Input_Parameters.Flow_Angle = 180.00*Input_Parameters.Flow_Angle/PI;
        } /* endif */
        break;
      case GRID_READ_FROM_GRID_DATA_FILE :
        Grid_ptr = Read_Multi_Block_Grid(Grid_ptr,
                                         Input_Parameters);
        if (Grid_ptr == NULL) {
	   cout << "\n " << CFFC_Name() 
                << " Euler2D ERROR: Unable to open multi-block mesh data file "
                << Input_Parameters.Grid_File_Name << ".\n";
        } /* endif */
        if (Input_Parameters.i_Grid == GRID_NASA_ROTOR_37) {
           Input_Parameters.NASA_Rotor37.init(Input_Parameters.Rotor_Flow_Type,
                                              Input_Parameters.NASA_Rotor37_Data_Directory);
           Input_Parameters.Wo.setgas("AIR");
	   Input_Parameters.NASA_Rotor37.getPstateREL_up(Input_Parameters.Wo,Input_Parameters.Rotor_Percent_Span);
        } else if (Input_Parameters.i_Grid == GRID_NASA_ROTOR_67) {
           Input_Parameters.NASA_Rotor67.init(Input_Parameters.Rotor_Flow_Type,
                                              Input_Parameters.NASA_Rotor67_Data_Directory);
           Input_Parameters.Wo.setgas("AIR");
	   Input_Parameters.NASA_Rotor67.getPstateREL_up(Input_Parameters.Wo,Input_Parameters.Rotor_Percent_Span);
        } /* endif */
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
	if (!Input_Parameters.i_Mesh_Stretching) {
	  Grid_ptr = Grid_Rectangular_Box(Grid_ptr,
					  Input_Parameters.Number_of_Blocks_Idir,
					  Input_Parameters.Number_of_Blocks_Jdir,
					  Input_Parameters.Box_Width,
					  Input_Parameters.Box_Height,
					  Input_Parameters.Number_of_Cells_Idir,
					  Input_Parameters.Number_of_Cells_Jdir,
					  Input_Parameters.Number_of_Ghost_Cells);
	} else {
	  Grid_ptr = Grid_Rectangular_Box(Grid_ptr,
					  Input_Parameters.Number_of_Blocks_Idir,
					  Input_Parameters.Number_of_Blocks_Jdir,
					  Input_Parameters.Box_Width,
					  Input_Parameters.Box_Height,
					  Input_Parameters.i_Mesh_Stretching,
					  Input_Parameters.Mesh_Stretching_Type_Idir,
					  Input_Parameters.Mesh_Stretching_Type_Jdir,
					  Input_Parameters.Mesh_Stretching_Factor_Idir,
					  Input_Parameters.Mesh_Stretching_Factor_Jdir,
					  Input_Parameters.Number_of_Cells_Idir,
					  Input_Parameters.Number_of_Cells_Jdir,
					  Input_Parameters.Number_of_Ghost_Cells);
	}
        break;
      case GRID_FLAT_PLATE :
	Grid_ptr = Grid_Flat_Plate(Grid_ptr,
				   Input_Parameters.Number_of_Blocks_Idir,
 		                   Input_Parameters.Number_of_Blocks_Jdir,
				   Input_Parameters.Plate_Length,
				   Input_Parameters.BC_South,
				   Input_Parameters.i_Mesh_Stretching,
				   Input_Parameters.Mesh_Stretching_Factor_Idir,
				   Input_Parameters.Mesh_Stretching_Factor_Jdir,
  		                   Input_Parameters.Number_of_Cells_Idir,
 		                   Input_Parameters.Number_of_Cells_Jdir,
 				   Input_Parameters.Number_of_Ghost_Cells);
        break;
      case GRID_PIPE :
        Grid_ptr = Grid_Pipe(Grid_ptr,
                             Input_Parameters.Number_of_Blocks_Idir,
		             Input_Parameters.Number_of_Blocks_Jdir,
                             Input_Parameters.Pipe_Length,
                             Input_Parameters.Pipe_Radius,0,ZERO,
 		             Input_Parameters.Number_of_Cells_Idir,
		             Input_Parameters.Number_of_Cells_Jdir,
			     Input_Parameters.Number_of_Ghost_Cells);
        break;
      case GRID_BLUNT_BODY :
        Grid_ptr = Grid_Blunt_Body(Grid_ptr,
                                   Input_Parameters.Number_of_Blocks_Idir,
                                   Input_Parameters.Number_of_Blocks_Jdir,
                                   Input_Parameters.Blunt_Body_Radius,
                                   Input_Parameters.Blunt_Body_Mach_Number,
    	                           Input_Parameters.Number_of_Cells_Idir,
    	                           Input_Parameters.Number_of_Cells_Jdir,
				   Input_Parameters.Number_of_Ghost_Cells);
        break;
      case GRID_ROCKET_MOTOR :
        Grid_ptr = Grid_Rocket_Motor(Grid_ptr,
                                     Input_Parameters.Number_of_Blocks_Idir,
                                     Input_Parameters.Number_of_Blocks_Jdir,
				     Input_Parameters.Chamber_Length,
				     Input_Parameters.Chamber_Radius,
				     Input_Parameters.Chamber_To_Throat_Length,
				     Input_Parameters.Nozzle_Length,
				     Input_Parameters.Nozzle_Radius_Exit,
				     Input_Parameters.Nozzle_Radius_Throat,
				     Input_Parameters.Grain_Radius,
				     Input_Parameters.Nozzle_Type,
				     Input_Parameters.BC_North,
				     Input_Parameters.i_Mesh_Stretching,
				     Input_Parameters.Mesh_Stretching_Type_Jdir,
				     Input_Parameters.Mesh_Stretching_Factor_Idir,
				     Input_Parameters.Mesh_Stretching_Factor_Jdir,
   	                             Input_Parameters.Number_of_Cells_Idir,
    	                             Input_Parameters.Number_of_Cells_Jdir,
				     Input_Parameters.Number_of_Ghost_Cells);
	if (Input_Parameters.Interface_IP.Component_List.Ni) {
	  for (int jBlk = 0; jBlk < Input_Parameters.Number_of_Blocks_Jdir; jBlk++) {
	    for (int iBlk = 0; iBlk < Input_Parameters.Number_of_Blocks_Idir; iBlk++) {
	      if (Grid_ptr[iBlk][jBlk].Node != NULL) {
		if (jBlk == Input_Parameters.Number_of_Blocks_Jdir-1)
		  Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(Input_Parameters.BC_North);
		Set_BCs(Grid_ptr[iBlk][jBlk]);
		Update_Exterior_Nodes(Grid_ptr[iBlk][jBlk]);
		Update_Cells(Grid_ptr[iBlk][jBlk]);
	      }
	    }
	  }
	}
        break;
    case GRID_NOZZLELESS_ROCKET_MOTOR :
      Grid_ptr = Grid_Nozzleless_Rocket_Motor(Grid_ptr,
					      Input_Parameters.Number_of_Blocks_Idir,
					      Input_Parameters.Number_of_Blocks_Jdir,
					      Input_Parameters.Chamber_Length,
					      Input_Parameters.Chamber_Radius,
					      Input_Parameters.Nozzle_Length,
					      Input_Parameters.Nozzle_Radius_Exit,
					      BC_REFLECTION,
					      Input_Parameters.i_Mesh_Stretching,
					      Input_Parameters.Mesh_Stretching_Type_Jdir,
					      Input_Parameters.Mesh_Stretching_Factor_Idir,
					      Input_Parameters.Mesh_Stretching_Factor_Jdir,
					      Input_Parameters.Number_of_Cells_Idir,
					      Input_Parameters.Number_of_Cells_Jdir,
					      Input_Parameters.Number_of_Ghost_Cells);
        break;
      case GRID_CIRCULAR_CYLINDER :
        Grid_ptr = Grid_Circular_Cylinder(Grid_ptr,
                                          Input_Parameters.Number_of_Blocks_Idir,
		                          Input_Parameters.Number_of_Blocks_Jdir,
                                          Input_Parameters.Cylinder_Radius,
                                          Input_Parameters.Cylinder_Radius2,
					  Input_Parameters.Mesh_Stretching_Type_Idir,
					  Input_Parameters.Mesh_Stretching_Type_Jdir,
					  Input_Parameters.Mesh_Stretching_Factor_Idir,
					  Input_Parameters.Mesh_Stretching_Factor_Jdir,
 		                          Input_Parameters.Number_of_Cells_Idir,
		                          Input_Parameters.Number_of_Cells_Jdir,
					  Input_Parameters.Number_of_Ghost_Cells);
        break;
      case GRID_ELLIPSE :
        Grid_ptr = Grid_Ellipse(Grid_ptr,
                                Input_Parameters.Number_of_Blocks_Idir,
		                Input_Parameters.Number_of_Blocks_Jdir,
                                Input_Parameters.Ellipse_Length_X_Axis,
                                Input_Parameters.Ellipse_Length_Y_Axis,
 		                Input_Parameters.Number_of_Cells_Idir,
		                Input_Parameters.Number_of_Cells_Jdir,
				Input_Parameters.Number_of_Ghost_Cells);
        break;
      case GRID_NACA_AEROFOIL :
        Grid_ptr = Grid_NACA_Aerofoil(Grid_ptr,
                                      Input_Parameters.Number_of_Blocks_Idir,
		                      Input_Parameters.Number_of_Blocks_Jdir,
                                      Input_Parameters.NACA_Aerofoil_Type,
                                      Input_Parameters.Chord_Length,
 		                      Input_Parameters.Number_of_Cells_Idir,
		                      Input_Parameters.Number_of_Cells_Jdir,
				      Input_Parameters.Number_of_Ghost_Cells);
        break;
      case GRID_FREE_JET :
        Grid_ptr = Grid_Free_Jet(Grid_ptr,
                                 Input_Parameters.Number_of_Blocks_Idir,
		                 Input_Parameters.Number_of_Blocks_Jdir,
                                 Input_Parameters.Orifice_Radius,
 		                 Input_Parameters.Number_of_Cells_Idir,
		                 Input_Parameters.Number_of_Cells_Jdir,
				 Input_Parameters.Number_of_Ghost_Cells);
        break;
      case GRID_WEDGE :
	Grid_ptr = Grid_Wedge(Grid_ptr,
                              Input_Parameters.Number_of_Blocks_Idir,
		              Input_Parameters.Number_of_Blocks_Jdir,
			      Input_Parameters.Wedge_Angle,
			      Input_Parameters.Wedge_Length,
			      Input_Parameters.BC_South,
			      Input_Parameters.i_Mesh_Stretching,
			      Input_Parameters.Mesh_Stretching_Factor_Idir,
			      Input_Parameters.Mesh_Stretching_Factor_Jdir,
			      Input_Parameters.Number_of_Cells_Idir,
			      Input_Parameters.Number_of_Cells_Jdir,
			      Input_Parameters.Number_of_Ghost_Cells);
	break;
      case GRID_UNSTEADY_BLUNT_BODY :
	Grid_ptr = Grid_Unsteady_Blunt_Body(Grid_ptr,
                                            Input_Parameters.Number_of_Blocks_Idir,
		                            Input_Parameters.Number_of_Blocks_Jdir,
					    Input_Parameters.Blunt_Body_Radius,
					    Input_Parameters.Blunt_Body_Mach_Number,
					    Input_Parameters.Number_of_Cells_Idir,
					    Input_Parameters.Number_of_Cells_Jdir,
					    Input_Parameters.Number_of_Ghost_Cells);
	break;
      case GRID_RINGLEB_FLOW :
        Grid_ptr = Grid_Ringleb_Flow(Grid_ptr,
                                     Input_Parameters.Number_of_Blocks_Idir,
	    	                     Input_Parameters.Number_of_Blocks_Jdir,
				     Input_Parameters.Inner_Streamline_Number,
				     Input_Parameters.Outer_Streamline_Number,
				     Input_Parameters.Isotach_Line,
                                     Input_Parameters.Number_of_Cells_Idir,
                                     Input_Parameters.Number_of_Cells_Jdir,
				     Input_Parameters.Number_of_Ghost_Cells);
        break;
      case GRID_BUMP_CHANNEL_FLOW :
        Grid_ptr = Grid_Bump_Channel_Flow(Grid_ptr,
                                          Input_Parameters.Number_of_Blocks_Idir,
		                          Input_Parameters.Number_of_Blocks_Jdir,
					  Input_Parameters.Smooth_Bump,
 		                          Input_Parameters.Number_of_Cells_Idir,
		                          Input_Parameters.Number_of_Cells_Jdir,
					  Input_Parameters.Number_of_Ghost_Cells);
        break;
      case GRID_NASA_ROTOR_37 :
        Input_Parameters.Number_of_Blocks_Idir = 3;
        Input_Parameters.Number_of_Blocks_Jdir = 2;
	if (!Input_Parameters.Interface_IP.Component_List.Ni) {
	  Grid_ptr = Input_Parameters.NASA_Rotor37.genMeshH_3x2_AUTO(Grid_ptr,
								     Input_Parameters.Rotor_Percent_Span,
								     Input_Parameters.NASA_Rotor37.z_up,
								     Input_Parameters.NASA_Rotor37.z_d,
								     Input_Parameters.Number_of_Cells_Idir,
								     Input_Parameters.Number_of_Cells_Jdir,
								     Input_Parameters.Number_of_Ghost_Cells,
								     min(250, 8*max(Input_Parameters.Number_of_Cells_Idir,
										    Input_Parameters.Number_of_Cells_Jdir)));
	} else {
	  Grid_ptr = Grid_NASA_Rotor_37(Grid_ptr,
					Input_Parameters.Number_of_Blocks_Idir,
					Input_Parameters.Number_of_Blocks_Jdir,
					Input_Parameters.Rotor_Percent_Span,
					Input_Parameters.Number_of_Cells_Idir, 
					Input_Parameters.Number_of_Cells_Jdir,
					Input_Parameters.Number_of_Ghost_Cells);
	}
        break;
      case GRID_NASA_ROTOR_67 :
        Input_Parameters.Number_of_Blocks_Idir = 3;
        Input_Parameters.Number_of_Blocks_Jdir = 2;
	if (!Input_Parameters.Interface_IP.Component_List.Ni) {
	  Grid_ptr = Input_Parameters.NASA_Rotor67.genMeshH_3x2_AUTO(Grid_ptr,
								     Input_Parameters.Rotor_Percent_Span,
								     Input_Parameters.NASA_Rotor67.z_up,
								     Input_Parameters.NASA_Rotor67.z_d,
								     Input_Parameters.Number_of_Cells_Idir,
								     Input_Parameters.Number_of_Cells_Jdir,
								     Input_Parameters.Number_of_Ghost_Cells,
								     min(250, 8*max(Input_Parameters.Number_of_Cells_Idir,
										    Input_Parameters.Number_of_Cells_Jdir)));
	} else {
	  Grid_ptr = Grid_NASA_Rotor_67(Grid_ptr,
					Input_Parameters.Number_of_Blocks_Idir,
					Input_Parameters.Number_of_Blocks_Jdir,
					Input_Parameters.Rotor_Percent_Span,
					Input_Parameters.Number_of_Cells_Idir,
					Input_Parameters.Number_of_Cells_Jdir,
					Input_Parameters.Number_of_Ghost_Cells);
	}
        break;
      case GRID_ICEMCFD :
	Grid_ptr = ICEMCFD_Read(Input_Parameters.ICEMCFD_FileNames,
				Grid_ptr,
				&Input_Parameters.Number_of_Blocks_Idir,
				&Input_Parameters.Number_of_Blocks_Jdir);
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
    } /* endswitch */

    /* Reset boundary conditions if required. */

    if (Input_Parameters.BCs_Specified) {
      for (int jBlk = 0; jBlk < Input_Parameters.Number_of_Blocks_Jdir; jBlk++) {
	for (int iBlk = 0; iBlk < Input_Parameters.Number_of_Blocks_Idir; iBlk++) {
	  if (jBlk == Input_Parameters.Number_of_Blocks_Jdir-1)
	    Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(Input_Parameters.BC_North);
	  if (jBlk == 0)
	    Grid_ptr[iBlk][jBlk].BndSouthSpline.setBCtype(Input_Parameters.BC_South);
	  if (iBlk == Input_Parameters.Number_of_Blocks_Idir-1)
	    Grid_ptr[iBlk][jBlk].BndEastSpline.setBCtype(Input_Parameters.BC_East);
	  if (iBlk == 0)
	    Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(Input_Parameters.BC_West);
	  Set_BCs(Grid_ptr[iBlk][jBlk]);
	  Update_Exterior_Nodes(Grid_ptr[iBlk][jBlk]);
	  Update_Cells(Grid_ptr[iBlk][jBlk]);
	}
      }
    }

    /* First translate quadrilateral mesh as specified by input parameters. */

    if (abs(Input_Parameters.X_Shift) > TOLER) {
       Translate_Multi_Block_Grid(Grid_ptr,
                                  Input_Parameters.Number_of_Blocks_Idir,
                                  Input_Parameters.Number_of_Blocks_Jdir,
                                  Input_Parameters.X_Shift);
    } /* endif */

    /* Next scale quadrilateral mesh as specified by input parameters. */

    if (fabs(Input_Parameters.X_Scale-ONE) > TOLER) {
       Scale_Multi_Block_Grid(Grid_ptr,
                              Input_Parameters.Number_of_Blocks_Idir,
                              Input_Parameters.Number_of_Blocks_Jdir,
                              Input_Parameters.X_Scale);
    } /* endif */

    /* Finally rotate quadrilateral mesh as specified by input parameters. */

    if (fabs(Input_Parameters.X_Rotate) > TOLER) {
       Rotate_Multi_Block_Grid(Grid_ptr,
                               Input_Parameters.Number_of_Blocks_Idir,
                               Input_Parameters.Number_of_Blocks_Jdir,
                               TWO*PI*Input_Parameters.X_Rotate/360.00);
    } /* endif */

    /* Return the multi-block quadilateral mesh. */

    return(Grid_ptr);

}

/********************************************************
 * Routine: Broadcast_Multi_Block_Grid                  *
 *                                                      *
 * Broadcast multi-block quadilateral mesh to all       *
 * processors involved in the calculation from the      *
 * primary processor using the MPI broadcast routine.   *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Broadcast_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                               Euler2D_Input_Parameters &Input_Parameters) {

    /* Broadcast the multi-block quadilateral mesh. */

#ifdef _MPI_VERSION
    MPI::COMM_WORLD.Bcast(&(Input_Parameters.Number_of_Blocks_Idir), 
                          1, 
                          MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Input_Parameters.Number_of_Blocks_Jdir), 
                          1, 
                          MPI::INT, 0);
    Grid_ptr = Broadcast_Multi_Block_Grid(Grid_ptr,
			                  Input_Parameters.Number_of_Blocks_Idir,
		                          Input_Parameters.Number_of_Blocks_Jdir);
#endif

    /* Return the multi-block quadilateral mesh. */

    return(Grid_ptr);

}

/********************************************************
 * Routine: Write_Multi_Block_Grid_Definition           *
 *                                                      *
 * Writes a grid definition file for a multi-block      *
 * quadilateral mesh in a format suitable for retrieval *
 * and re-use purposes.  Returns a non-zero value if    *
 * unable to write the grid definition file.            *
 *                                                      *
 ********************************************************/
int Write_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                      Euler2D_Input_Parameters &Input_Parameters) {

    char *mesh_definition_file_name_ptr;
    ofstream mesh_definition_file;

    /* Open the grid definition file. */

    mesh_definition_file_name_ptr = Input_Parameters.Grid_Definition_File_Name;
    mesh_definition_file.open(mesh_definition_file_name_ptr, ios::out);
    if (mesh_definition_file.fail()) return (1);

    /* Write grid type information. */

    mesh_definition_file << Input_Parameters.Grid_Type << "\n" 
                         << Input_Parameters.i_Grid << "\n";

    /* Write the grid definition information for each 
       quadrilateral grid block. */

    Write_Multi_Block_Grid_Definition(Grid_ptr,
			              Input_Parameters.Number_of_Blocks_Idir,
		                      Input_Parameters.Number_of_Blocks_Jdir,
                                      mesh_definition_file);

    /* Close the grid definition file. */

    mesh_definition_file.close();

    /* Writing of grid definition file complete.  Return zero value. */

    return(0);

}

/********************************************************
 * Routine: Read_Multi_Block_Grid_Definition            *
 *                                                      *
 * Reads a grid definition file for a multi-block       *
 * quadilateral mesh.  Returns a pointer to the mesh.   *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Read_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                                     Euler2D_Input_Parameters &Input_Parameters) {

    char buffer[256];
    char *mesh_definition_file_name_ptr;
    ifstream mesh_definition_file;

    /* Open the grid definition file. */

    mesh_definition_file_name_ptr = Input_Parameters.Grid_Definition_File_Name;
    mesh_definition_file.open(mesh_definition_file_name_ptr, ios::in);
    if (mesh_definition_file.fail()) return (NULL);

    /* Read grid type information. */

    mesh_definition_file.getline(buffer, sizeof(buffer));
    strcpy(Input_Parameters.Grid_Type, buffer);
    mesh_definition_file >> Input_Parameters.i_Grid;

    /* Read the grid definition information for each 
       quadrilateral grid block. */

    Grid_ptr = Read_Multi_Block_Grid_Definition(Grid_ptr,
		              	                Input_Parameters.Number_of_Blocks_Idir,
		                                Input_Parameters.Number_of_Blocks_Jdir,
                                                mesh_definition_file);

    /* Close the grid definition file. */

    mesh_definition_file.close();

    /* Reading of grid definition file complete.  Return the mesh. */

    return(Grid_ptr);

}

/********************************************************
 * Routine: Write_Multi_Block_Grid                      *
 *                                                      *
 * Writes multi-block quadilateral mesh to a grid data  *
 * file in a format suitable for retrieval and re-use   *
 * purposes.  Returns a non-zero value if unable to     *
 * write the grid data file.                            *
 *                                                      *
 ********************************************************/
int Write_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                           Euler2D_Input_Parameters &Input_Parameters) {

    char *mesh_file_name_ptr;
    ofstream mesh_file;

    /* Open the grid data output file. */

    mesh_file_name_ptr = Input_Parameters.Grid_File_Name;
    mesh_file.open(mesh_file_name_ptr, ios::out);
    if (mesh_file.fail()) return (1);

    /* Write grid type information. */

    mesh_file << Input_Parameters.Grid_Type << "\n" 
              << Input_Parameters.i_Grid << "\n";

    /* Write the grid data for each quadrilateral grid block. */

    Write_Multi_Block_Grid(Grid_ptr,
			   Input_Parameters.Number_of_Blocks_Idir,
		           Input_Parameters.Number_of_Blocks_Jdir,
                           mesh_file);

    /* Close the grid data output file. */

    mesh_file.close();

    /* Writing of grid data output file complete.  Return zero value. */

    return(0);

}

/********************************************************
 * Routine: Read_Multi_Block_Grid                       *
 *                                                      *
 * Reads multi-block quadilateral mesh from a grid data *
 * file.  Returns a pointer to the mesh.                *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Read_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                          Euler2D_Input_Parameters &Input_Parameters) {

    char buffer[256];
    char *mesh_file_name_ptr;
    ifstream mesh_file;

    /* Open the grid data input file. */

    mesh_file_name_ptr = Input_Parameters.Grid_File_Name;
    mesh_file.open(mesh_file_name_ptr, ios::in);
    if (mesh_file.fail()) return (NULL);

    /* Read grid type information. */

    mesh_file.getline(buffer, sizeof(buffer));
    strcpy(Input_Parameters.Grid_Type, buffer);
    mesh_file >> Input_Parameters.i_Grid;

    /* Read the grid data for each quadrilateral grid block. */

    Grid_ptr = Read_Multi_Block_Grid(Grid_ptr,
                                     Input_Parameters.Number_of_Blocks_Idir,
		                     Input_Parameters.Number_of_Blocks_Jdir,
                                     mesh_file);

    /* Close the grid data input file. */

    mesh_file.close();

    /* Reading of grid data input file complete.  Return the mesh. */

    return(Grid_ptr);

}

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the nodes of a multi-block quadilateral mesh  *
 * in a format suitable for plotting with TECPLOT.      *
 * Returns a non-zero value if unable to write the      *
 * TECPLOT file.                                        *
 *                                                      *
 ********************************************************/
int Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                   Euler2D_Input_Parameters &Input_Parameters) {

    int i, j;
    char prefix[256], extension[256], mesh_file_name[256];
    char *mesh_file_name_ptr;
    ofstream mesh_file;  

    /* Determine prefix of grid data output file name. */

    i = 0;
    while (1) {
       if (Input_Parameters.Output_File_Name[i] == ' ' ||
           Input_Parameters.Output_File_Name[i] == '.') break;
       prefix[i]=Input_Parameters.Output_File_Name[i];
       i = i + 1;
       if (i > strlen(Input_Parameters.Output_File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';
    strcat(prefix, "_mesh");

    /* Determine grid data output file name. */

    strcpy(extension, ".dat");
    strcpy(mesh_file_name, prefix);
    strcat(mesh_file_name, extension);
    mesh_file_name_ptr = mesh_file_name;

    /* Open the grid data output file. */

    mesh_file.open(mesh_file_name_ptr, ios::out);
    if (mesh_file.fail()) return (1);

    /* Write the node locations for each quadrilateral grid block. */

    Output_Tecplot(Grid_ptr,
		   Input_Parameters.Number_of_Blocks_Idir,
		   Input_Parameters.Number_of_Blocks_Jdir,
	           mesh_file);

    /* Close the grid data output file. */

    mesh_file.close();

    /* Writing of grid data output file complete.  Return zero value. */

    return(0);

}

/********************************************************
 * Routine: Output_Nodes_Tecplot                        *
 *                                                      *
 * Writes the nodes of a multi-block quadilateral mesh  *
 * in a format suitable for plotting with TECPLOT.      *
 * Includes boundary nodes.  Returns a non-zero value   *
 * if unable to write the TECPLOT file.                 *
 *                                                      *
 ********************************************************/
int Output_Nodes_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                         Euler2D_Input_Parameters &Input_Parameters) {

    int i, j;
    char prefix[256], extension[256], mesh_file_name[256];
    char *mesh_file_name_ptr;
    ofstream mesh_file;  

    /* Determine prefix of grid data output file name. */

    i = 0;
    while (1) {
       if (Input_Parameters.Output_File_Name[i] == ' ' ||
           Input_Parameters.Output_File_Name[i] == '.') break;
       prefix[i]=Input_Parameters.Output_File_Name[i];
       i = i + 1;
       if (i > strlen(Input_Parameters.Output_File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';
    strcat(prefix, "_mesh_nodes");

    /* Determine grid data output file name. */

    strcpy(extension, ".dat");
    strcpy(mesh_file_name, prefix);
    strcat(mesh_file_name, extension);
    mesh_file_name_ptr = mesh_file_name;

    /* Open the grid data output file. */

    mesh_file.open(mesh_file_name_ptr, ios::out);
    if (mesh_file.fail()) return (1);

    /* Write the node locations for each quadrilateral grid block. */

    Output_Nodes_Tecplot(Grid_ptr,
		         Input_Parameters.Number_of_Blocks_Idir,
		         Input_Parameters.Number_of_Blocks_Jdir,
	                 mesh_file);

    /* Close the grid data output file. */

    mesh_file.close();

    /* Writing of grid data output file complete.  Return zero value. */

    return(0);

}

/********************************************************
 * Routine: Output_Cells_Tecplot                        *
 *                                                      *
 * Writes the cells of a multi-block quadilateral mesh  *
 * in a format suitable for plotting with TECPLOT.      *
 * Returns a non-zero value if unable to write the      *
 * TECPLOT file.                                        *
 *                                                      *
 ********************************************************/
int Output_Cells_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                         Euler2D_Input_Parameters &Input_Parameters) {

    int i, j;
    char prefix[256], extension[256], mesh_file_name[256];
    char *mesh_file_name_ptr;
    ofstream mesh_file;  

    /* Determine prefix of grid data output file name. */

    i = 0;
    while (1) {
       if (Input_Parameters.Output_File_Name[i] == ' ' ||
           Input_Parameters.Output_File_Name[i] == '.') break;
       prefix[i]=Input_Parameters.Output_File_Name[i];
       i = i + 1;
       if (i > strlen(Input_Parameters.Output_File_Name) ) break;
    } /* endwhile */
    prefix[i] = '\0';
    strcat(prefix, "_mesh_cells");

    /* Determine grid data output file name. */

    strcpy(extension, ".dat");
    strcpy(mesh_file_name, prefix);
    strcat(mesh_file_name, extension);
    mesh_file_name_ptr = mesh_file_name;

    /* Open the grid data output file. */

    mesh_file.open(mesh_file_name_ptr, ios::out);
    if (mesh_file.fail()) return (1);

    /* Write the node locations for each quadrilateral grid block. */

    Output_Cells_Tecplot(Grid_ptr,
		         Input_Parameters.Number_of_Blocks_Idir,
		         Input_Parameters.Number_of_Blocks_Jdir,
	                 mesh_file);

    /* Close the grid data output file. */

    mesh_file.close();

    /* Writing of grid data output file complete.  Return zero value. */

    return(0);

}
