/* Flame2DQuadGrid.cc:  Grid Manipulation Subroutines for 2D Euler 
                        Multi-Block Quadrilateral Mesh 
                        Solution Classes. */

/* Include 2D Euler quadrilateral mesh solution header file. */

#ifndef _FLAME2D_QUAD_INCLUDED
#include "Flame2DQuad.h"
#endif // _FLAME2D_QUAD_INCLUDED

/**************************************************************************
 * Flame2D_Quad_Block -- Multiple Block External Subroutines for Mesh.    *
 **************************************************************************/

/********************************************************
 * Routine: Multi_Block_Grid                            *
 *                                                      *
 * Generates multi-block quadilateral mesh.             *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                     Flame2D_Input_Parameters &Input_Parameters) {

    int iBlk, jBlk;

    /* Generate appropriate mesh. */

    switch(Input_Parameters.i_Grid) {
      case GRID_READ_FROM_DEFINITION_FILE :
        Grid_ptr = Read_Multi_Block_Grid_Definition(Grid_ptr,
                                                    Input_Parameters);
        if (Grid_ptr == NULL) {
	   cout << "\n " << CFFC_Name() 
                << " Flame2D ERROR: Unable to open multi-block mesh definition file "
                << Input_Parameters.Grid_Definition_File_Name << ".\n";
        } /* endif */
        break;
      case GRID_READ_FROM_GRID_DATA_FILE :
        Grid_ptr = Read_Multi_Block_Grid(Grid_ptr,
                                         Input_Parameters);
        if (Grid_ptr == NULL) {
	   cout << "\n " << CFFC_Name() 
                << " Flame2D ERROR: Unable to open multi-block mesh data file "
                << Input_Parameters.Grid_File_Name << ".\n";
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
        for ( jBlk = 0; jBlk <= Input_Parameters.Number_of_Blocks_Jdir-1; ++jBlk ) {
	  for ( iBlk = 0; iBlk <= Input_Parameters.Number_of_Blocks_Idir-1; ++iBlk ) {
	    if (jBlk == Input_Parameters.Number_of_Blocks_Jdir-1) {
	      Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(BC_CONSTANT_EXTRAPOLATION);
	    } else {
	      Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(BC_NONE);
	    } /* endif */
	    if (jBlk == 0) {
	      Grid_ptr[iBlk][jBlk].BndSouthSpline.setBCtype(BC_CONSTANT_EXTRAPOLATION);  //BC_FIXED
	    } else {
	      Grid_ptr[iBlk][jBlk].BndSouthSpline.setBCtype(BC_NONE);
	    } /* endif */
	    if (iBlk == Input_Parameters.Number_of_Blocks_Idir-1) {
	      Grid_ptr[iBlk][jBlk].BndEastSpline.setBCtype(BC_CONSTANT_EXTRAPOLATION);
	    } else {
	      Grid_ptr[iBlk][jBlk].BndEastSpline.setBCtype(BC_NONE);
	    } /* endif */
	    if (iBlk == 0) {
		 Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(BC_CONSTANT_EXTRAPOLATION);
	    } else {
	      Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(BC_NONE);
	    } /* endif */
	    Set_BCs(Grid_ptr[iBlk][jBlk]);
	    Update_Exterior_Nodes(Grid_ptr[iBlk][jBlk]);
	    Update_Cells(Grid_ptr[iBlk][jBlk]);
	  } 
        } 
	
        break;
    
	//***************** FOR VISCOUS FLOWS *************************//
    case GRID_COUETTE :
        Grid_ptr = Grid_Rectangular_Box(Grid_ptr,
                                        Input_Parameters.Number_of_Blocks_Idir,
		                        Input_Parameters.Number_of_Blocks_Jdir,
                                        Input_Parameters.Box_Width,
                                        Input_Parameters.Box_Height,
 		                        Input_Parameters.Number_of_Cells_Idir,
		                        Input_Parameters.Number_of_Cells_Jdir,
					Input_Parameters.Number_of_Ghost_Cells);
        for ( jBlk = 0; jBlk <= Input_Parameters.Number_of_Blocks_Jdir-1; ++jBlk ) {
	  for ( iBlk = 0; iBlk <= Input_Parameters.Number_of_Blocks_Idir-1; ++iBlk ) {
	    // North (Top)
	    if (jBlk == Input_Parameters.Number_of_Blocks_Jdir-1) {
	      Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(BC_MOVING_WALL_HEATFLUX);
	    } else {
                 Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(BC_NONE);
	    } 
	    // South (Bottom)
	    if (jBlk == 0) {
	      Grid_ptr[iBlk][jBlk].BndSouthSpline.setBCtype(BC_WALL_VISCOUS_HEATFLUX);
	    } else {
	      Grid_ptr[iBlk][jBlk].BndSouthSpline.setBCtype(BC_NONE);
	    } 
	    // East (Outflow)
	    if (iBlk == Input_Parameters.Number_of_Blocks_Idir-1) {
	      Grid_ptr[iBlk][jBlk].BndEastSpline.setBCtype(BC_OUTFLOW_SUBSONIC);
	    } else {
	      Grid_ptr[iBlk][jBlk].BndEastSpline.setBCtype(BC_NONE);
	    }
	    // West (Inflow)
	    if (iBlk == 0) {  
	      Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(BC_INFLOW_SUBSONIC);
	    } else {
	      Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(BC_NONE);
	    } 
	    
	    Set_BCs(Grid_ptr[iBlk][jBlk]);
	    Update_Exterior_Nodes(Grid_ptr[iBlk][jBlk]);
	    Update_Cells(Grid_ptr[iBlk][jBlk]);
	  } 
        } 
        break;

	/*************** TEST MULTIBLOCK GRID *******************************/
	// Primarily used to test BC's convergence issues with NKS
    case GRID_TEST :
      Grid_ptr = Grid_Rectangular_Box(Grid_ptr,
				      Input_Parameters.Number_of_Blocks_Idir,
				      Input_Parameters.Number_of_Blocks_Jdir,
				      Input_Parameters.Box_Width,
				      Input_Parameters.Box_Height,
				      Input_Parameters.Number_of_Cells_Idir,
				      Input_Parameters.Number_of_Cells_Jdir,
				      Input_Parameters.Number_of_Ghost_Cells);
      for ( jBlk = 0; jBlk <= Input_Parameters.Number_of_Blocks_Jdir-1; ++jBlk ) {
	for ( iBlk = 0; iBlk <= Input_Parameters.Number_of_Blocks_Idir-1; ++iBlk ) {
	  // North (Top)
	  if (jBlk == Input_Parameters.Number_of_Blocks_Jdir-1) {
	    Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(BC_2DFLAME_OUTFLOW);
	  } else {
	    Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(BC_NONE);
	  } 
	  // South (Bottom)
	  if (iBlk == 0) {
	    Grid_ptr[iBlk][jBlk].BndSouthSpline.setBCtype(BC_FIXED); //BC_2DFLAME_INFLOW); //
	  } else if ( iBlk == 1) {
	    Grid_ptr[iBlk][jBlk].BndSouthSpline.setBCtype(BC_2DFLAME_INFLOW); //BC_WALL_VISCOUS_HEATFLUX);
	  } else {
	    Grid_ptr[iBlk][jBlk].BndSouthSpline.setBCtype(BC_FIXED); 
	  } 

	  // East 
	  if (iBlk == Input_Parameters.Number_of_Blocks_Idir-1) {
	    Grid_ptr[iBlk][jBlk].BndEastSpline.setBCtype(BC_REFLECTION);
	  } else {
	    Grid_ptr[iBlk][jBlk].BndEastSpline.setBCtype(BC_NONE); 
	  }
	  // West
	  if (iBlk == 0) {  
	    Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(BC_REFLECTION);
	  } else {
	    Grid_ptr[iBlk][jBlk].BndWestSpline.setBCtype(BC_NONE);
	  } 
	    
	  Set_BCs(Grid_ptr[iBlk][jBlk]);
	  Update_Exterior_Nodes(Grid_ptr[iBlk][jBlk]);
	  Update_Cells(Grid_ptr[iBlk][jBlk]);
	} 
      } 
      break;
      
    case GRID_DRIVEN_CAVITY_FLOW :
      Grid_ptr = Grid_Driven_Cavity_Flow(Grid_ptr,
					 Input_Parameters.Number_of_Blocks_Idir,
					 Input_Parameters.Number_of_Blocks_Jdir,
					 Input_Parameters.Box_Width,
					 Input_Parameters.Box_Height,
					 Input_Parameters.Mesh_Stretching_Type_Idir,
					 Input_Parameters.Mesh_Stretching_Type_Jdir,
					 Input_Parameters.Mesh_Stretching_Factor_Idir,
					 Input_Parameters.Mesh_Stretching_Factor_Jdir,		
					 Input_Parameters.Number_of_Cells_Idir,
					 Input_Parameters.Number_of_Cells_Jdir,
					 Input_Parameters.Number_of_Ghost_Cells);

       break;

	//***************** FOR 1D Premixed FLAMES *************************//
    case GRID_1DFLAME :
      Grid_ptr = Grid_1D_Flame(Grid_ptr,
			       Input_Parameters.Number_of_Blocks_Idir,
			       Input_Parameters.Number_of_Blocks_Jdir,
			       Input_Parameters.Box_Width,
			       Input_Parameters.Box_Height,
			       Input_Parameters.Number_of_Cells_Idir,
			       Input_Parameters.Number_of_Cells_Jdir,
			       Input_Parameters.Number_of_Ghost_Cells);
      break;
      
      //************ FOR 2D Axisymmetric Diffusions Flames **************//
    case GRID_LAMINAR_FLAME :
      Grid_ptr = Grid_2D_Laminar_Flame(Grid_ptr,
				       Input_Parameters.Number_of_Blocks_Idir,
				       Input_Parameters.Number_of_Blocks_Jdir,
				       Input_Parameters.Pipe_Radius,
				       Input_Parameters.Pipe_Length,
				       Input_Parameters.Number_of_Cells_Idir,
				       Input_Parameters.Number_of_Cells_Jdir,
				       Input_Parameters.Number_of_Ghost_Cells,
				       Input_Parameters.i_ICs);
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
      case GRID_FLAT_PLATE :
        Grid_ptr = Grid_Flat_Plate_NK(Grid_ptr,
				      Input_Parameters.Number_of_Blocks_Idir,
				      Input_Parameters.Number_of_Blocks_Jdir,
				      Input_Parameters.Plate_Length, 
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
                             Input_Parameters.Pipe_Radius,
			     Input_Parameters.Axisymmetric,
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
	if (Input_Parameters.FlowType == FLOWTYPE_INVISCID) Input_Parameters.BC_North = BC_REFLECTION;
	else Input_Parameters.BC_North = BC_WALL_VISCOUS_ISOTHERMAL;
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
        break;
      case GRID_CIRCULAR_CYLINDER :
        Grid_ptr = Grid_Circular_Cylinder(Grid_ptr,
                                          Input_Parameters.Number_of_Blocks_Idir,
		                          Input_Parameters.Number_of_Blocks_Jdir,
                                          Input_Parameters.Cylinder_Radius,
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
			      Input_Parameters.Mesh_Stretching_Type_Idir,
			      Input_Parameters.Mesh_Stretching_Type_Jdir,
			      Input_Parameters.Mesh_Stretching_Factor_Idir,
			      Input_Parameters.Mesh_Stretching_Factor_Jdir,
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

    // Reset boundary conditions if required.
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
                                               Flame2D_Input_Parameters &Input_Parameters) {

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
                                      Flame2D_Input_Parameters &Input_Parameters) {

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
                                                     Flame2D_Input_Parameters &Input_Parameters) {

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
                           Flame2D_Input_Parameters &Input_Parameters) {

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
                                          Flame2D_Input_Parameters &Input_Parameters) {

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
                   Flame2D_Input_Parameters &Input_Parameters) {

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
                         Flame2D_Input_Parameters &Input_Parameters) {

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
                         Flame2D_Input_Parameters &Input_Parameters) {

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
