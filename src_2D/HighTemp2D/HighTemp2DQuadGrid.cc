/**********************************************************************
 * HighTemp2DQuadGrid.cc: Grid manipulation subroutines for 2D HT     *
 *                            N-S multi-block quadrilateral           *
 *                            mesh solution classes.                  *
 **********************************************************************/

#include "HighTemp2DQuad.h"

/**********************************************************************
 * HighTemp2D_Quad_Block -- Multiple Block External Subroutines   *
 *                              for Quadrilateral Mesh.               *
 **********************************************************************/

/**********************************************************************
 * Routine: Multi_Block_Grid                                          *
 *                                                                    *
 * Generates multi-block quadilateral mesh.                           *
 *                                                                    *
 **********************************************************************/
Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
				     HighTemp2D_Input_Parameters &IP) {

  //char NASA_Rotor_Data_Directory[128];
  //int Rotor_Flow_Type;

  // Generate appropriate mesh.
  switch(IP.i_Grid) {
  case GRID_READ_FROM_DEFINITION_FILE :
    Grid_ptr = Read_Multi_Block_Grid_Definition(Grid_ptr,IP);
    if (Grid_ptr == NULL) {
      cout << "\n " << CFFC_Name() 
	   << " HighTemp2D ERROR: Unable to open multi-block mesh definition file "
	   << IP.Grid_Definition_File_Name << ".\n";
    }
    if (IP.i_Grid == GRID_NASA_ROTOR_37) {
      IP.NASA_Rotor37.init(IP.Rotor_Flow_Type,IP.NASA_Rotor37_Data_Directory);
      IP.Wo.set_gas("AIR");
      IP.NASA_Rotor37.getPstateREL_up(IP.Wo,IP.Rotor_Percent_Span);
      IP.Pressure = IP.Wo.p;
      IP.Temperature = IP.Wo.T();
      IP.Mach_Number = IP.NASA_Rotor37.getMachREL_up(IP.Rotor_Percent_Span);
      IP.Flow_Angle = atan2(IP.Wo.v.y, IP.Wo.v.x); 
      if (IP.Flow_Angle < ZERO) IP.Flow_Angle = TWO*PI + IP.Flow_Angle;
      IP.Flow_Angle = 180.00*IP.Flow_Angle/PI;
    } else if (IP.i_Grid == GRID_NASA_ROTOR_67) {
      IP.NASA_Rotor67.init(IP.Rotor_Flow_Type,IP.NASA_Rotor67_Data_Directory);
      IP.Wo.set_gas("AIR");
      IP.NASA_Rotor67.getPstateREL_up(IP.Wo,IP.Rotor_Percent_Span);
      IP.Pressure = IP.Wo.p;
      IP.Temperature = IP.Wo.T();
      IP.Mach_Number = IP.NASA_Rotor67.getMachREL_up(IP.Rotor_Percent_Span);
      IP.Flow_Angle = atan2(IP.Wo.v.y, IP.Wo.v.x); 
      if (IP.Flow_Angle < ZERO) IP.Flow_Angle = TWO*PI + IP.Flow_Angle;
      IP.Flow_Angle = 180.00*IP.Flow_Angle/PI;
    }
    break;
  case GRID_READ_FROM_GRID_DATA_FILE :
    Grid_ptr = Read_Multi_Block_Grid(Grid_ptr,IP);
    if (Grid_ptr == NULL) {
      cout << "\n " << CFFC_Name() 
	   << " HighTemp2D ERROR: Unable to open multi-block mesh data file "
	   << IP.Grid_File_Name << ".\n";
    }
    if (IP.i_Grid == GRID_NASA_ROTOR_37) {
      IP.NASA_Rotor37.init(IP.Rotor_Flow_Type,IP.NASA_Rotor37_Data_Directory);
      IP.Wo.set_gas("AIR");
      IP.NASA_Rotor37.getPstateREL_up(IP.Wo,IP.Rotor_Percent_Span);
    } else if (IP.i_Grid == GRID_NASA_ROTOR_67) {
      IP.NASA_Rotor67.init(IP.Rotor_Flow_Type,IP.NASA_Rotor67_Data_Directory);
      IP.Wo.set_gas("AIR");
      IP.NASA_Rotor67.getPstateREL_up(IP.Wo,IP.Rotor_Percent_Span);
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
    if (!IP.i_Mesh_Stretching) {
      if (IP.i_ICs == IC_VISCOUS_CHANNEL_FLOW) {
	IP.Box_Width += IP.Box_Width/double(IP.Number_of_Cells_Idir-1);
	Grid_ptr = Grid_Rectangular_Box(Grid_ptr,
					IP.Number_of_Blocks_Idir,
					IP.Number_of_Blocks_Jdir,
					IP.Box_Width,
					IP.Box_Height,
					IP.Number_of_Cells_Idir,
					IP.Number_of_Cells_Jdir,
					IP.Number_of_Ghost_Cells);
	IP.Box_Width /= (ONE + ONE/double(IP.Number_of_Cells_Idir-1));
      } else {
	Grid_ptr = Grid_Rectangular_Box(Grid_ptr,
					IP.Number_of_Blocks_Idir,
					IP.Number_of_Blocks_Jdir,
					IP.Box_Width,
					IP.Box_Height,
					IP.Number_of_Cells_Idir,
					IP.Number_of_Cells_Jdir,
					IP.Number_of_Ghost_Cells);
      }
    } else {
      Grid_ptr = Grid_Rectangular_Box(Grid_ptr,
				      IP.Number_of_Blocks_Idir,
				      IP.Number_of_Blocks_Jdir,
				      IP.Box_Width,
				      IP.Box_Height,
				      IP.i_Mesh_Stretching,
				      IP.Mesh_Stretching_Type_Idir,
				      IP.Mesh_Stretching_Type_Jdir,
				      IP.Mesh_Stretching_Factor_Idir,
				      IP.Mesh_Stretching_Factor_Jdir,
				      IP.Number_of_Cells_Idir,
				      IP.Number_of_Cells_Jdir,
				      IP.Number_of_Ghost_Cells);
    }
    break;
  case GRID_FLAT_PLATE :
    Grid_ptr = Grid_Flat_Plate(Grid_ptr,
				IP.Number_of_Blocks_Idir,
				IP.Number_of_Blocks_Jdir,
				IP.Plate_Length,
				IP.BC_South, 
				IP.i_Mesh_Stretching,
				IP.Mesh_Stretching_Factor_Idir,
				IP.Mesh_Stretching_Factor_Jdir,
				IP.Number_of_Cells_Idir,
				IP.Number_of_Cells_Jdir,
				IP.Number_of_Ghost_Cells);
    break;
  case GRID_FLAT_PLATE_NK : 
        Grid_ptr = Grid_Flat_Plate_NK(Grid_ptr,
				      IP.Number_of_Blocks_Idir,
				      IP.Number_of_Blocks_Jdir,
				      IP.Plate_Length, 
				      IP.i_Mesh_Stretching,
				      IP.Mesh_Stretching_Factor_Idir,
				      IP.Mesh_Stretching_Factor_Jdir,
				      IP.Number_of_Cells_Idir,
				      IP.Number_of_Cells_Jdir,
				      IP.Number_of_Ghost_Cells);
    break;
  case GRID_PIPE :
    if (IP.FlowType == FLOWTYPE_LAMINAR) {
      IP.Pipe_Length += IP.Pipe_Length/double(IP.Number_of_Cells_Idir-1);
      Grid_ptr = Grid_Pipe(Grid_ptr,
			   IP.Number_of_Blocks_Idir,
			   IP.Number_of_Blocks_Jdir,
			   IP.Pipe_Length,
			   IP.Pipe_Radius,
			   IP.i_Mesh_Stretching,
			   IP.Mesh_Stretching_Factor_Jdir,
			   IP.Number_of_Cells_Idir,
			   IP.Number_of_Cells_Jdir,
			   IP.Number_of_Ghost_Cells);
      IP.Pipe_Length /= (ONE + ONE/double(IP.Number_of_Cells_Idir-1));
      IP.X_Shift.x = -HALF*IP.Pipe_Length/double(IP.Number_of_Cells_Idir-1);
    } else {
      Grid_ptr = Grid_Pipe(Grid_ptr,
			   IP.Number_of_Blocks_Idir,
			   IP.Number_of_Blocks_Jdir,
			   IP.Pipe_Length,
			   IP.Pipe_Radius,
			   IP.i_Mesh_Stretching,
			   IP.Mesh_Stretching_Factor_Jdir,
			   IP.Number_of_Cells_Idir,
			   IP.Number_of_Cells_Jdir,
			   IP.Number_of_Ghost_Cells);
    }
    break;
  case GRID_BLUNT_BODY :
    Grid_ptr = Grid_Blunt_Body(Grid_ptr,
			       IP.Number_of_Blocks_Idir,
			       IP.Number_of_Blocks_Jdir,
			       IP.Blunt_Body_Radius,
			       IP.Blunt_Body_Mach_Number,
			       IP.Number_of_Cells_Idir,
			       IP.Number_of_Cells_Jdir,
			       IP.Number_of_Ghost_Cells);
    break;
  case GRID_ROCKET_MOTOR :
    if (IP.FlowType == FLOWTYPE_INVISCID) IP.BC_North = BC_REFLECTION;
    else IP.BC_North = BC_WALL_VISCOUS_ISOTHERMAL;
    Grid_ptr = Grid_Rocket_Motor(Grid_ptr,
				 IP.Number_of_Blocks_Idir,
				 IP.Number_of_Blocks_Jdir,
				 IP.Chamber_Length,
				 IP.Chamber_Radius,
				 IP.Chamber_To_Throat_Length,
				 IP.Nozzle_Length,
				 IP.Nozzle_Radius_Exit,
				 IP.Nozzle_Radius_Throat,
				 IP.Grain_Radius,
				 IP.Nozzle_Type,
				 IP.BC_North,
				 IP.i_Mesh_Stretching,
				 IP.Mesh_Stretching_Type_Jdir,
				 IP.Mesh_Stretching_Factor_Idir,
				 IP.Mesh_Stretching_Factor_Jdir,
				 IP.Number_of_Cells_Idir,
				 IP.Number_of_Cells_Jdir,
				 IP.Number_of_Ghost_Cells);
    if (IP.Interface_IP.Component_List.Ni) {
      for (int jBlk = 0; jBlk < IP.Number_of_Blocks_Jdir; jBlk++) {
	for (int iBlk = 0; iBlk < IP.Number_of_Blocks_Idir; iBlk++) {
	  if (Grid_ptr[iBlk][jBlk].Node != NULL) {
	    if (jBlk == IP.Number_of_Blocks_Jdir-1)
	      Grid_ptr[iBlk][jBlk].BndNorthSpline.setBCtype(IP.BC_North);
	    Set_BCs(Grid_ptr[iBlk][jBlk]);
	    Update_Exterior_Nodes(Grid_ptr[iBlk][jBlk]);
	    Update_Cells(Grid_ptr[iBlk][jBlk]);
	  }
	}
      }
    }
    break;
  case GRID_NOZZLE :
    Grid_ptr = Grid_Nozzle(Grid_ptr,
			   IP.Number_of_Blocks_Idir,
			   IP.Number_of_Blocks_Jdir,
			   IP.Nozzle_Length,
			   IP.Chamber_Radius,
			   IP.Nozzle_Radius_Exit,
			   IP.Nozzle_Radius_Throat,
			   IP.Nozzle_Type,
			   IP.i_Mesh_Stretching,
			   IP.Mesh_Stretching_Type_Idir,
			   IP.Mesh_Stretching_Type_Jdir,
			   IP.Mesh_Stretching_Factor_Idir,
			   IP.Mesh_Stretching_Factor_Jdir,
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
  case GRID_ELLIPSE :
    Grid_ptr = Grid_Ellipse(Grid_ptr,
			    IP.Number_of_Blocks_Idir,
			    IP.Number_of_Blocks_Jdir,
			    IP.Ellipse_Length_X_Axis,
			    IP.Ellipse_Length_Y_Axis,
			    IP.Number_of_Cells_Idir,
			    IP.Number_of_Cells_Jdir,
			    IP.Number_of_Ghost_Cells);
    break;
  case GRID_NACA_AEROFOIL :
    Grid_ptr = Grid_NACA_Aerofoil(Grid_ptr,
				  IP.Number_of_Blocks_Idir,
				  IP.Number_of_Blocks_Jdir,
				  IP.NACA_Aerofoil_Type,
				  IP.Chord_Length,
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
  case GRID_WEDGE :
    if (IP.FlowType == FLOWTYPE_INVISCID) IP.BC_South = BC_REFLECTION;
    Grid_ptr = Grid_Wedge(Grid_ptr,
			  IP.Number_of_Blocks_Idir,
			  IP.Number_of_Blocks_Jdir,
			  IP.Wedge_Angle,
			  IP.Wedge_Length,
			  IP.BC_South,
			  IP.i_Mesh_Stretching,
			  IP.Mesh_Stretching_Factor_Idir,
			  IP.Mesh_Stretching_Factor_Jdir,
			  IP.Number_of_Cells_Idir,
			  IP.Number_of_Cells_Jdir,
			  IP.Number_of_Ghost_Cells);
    break;
  case GRID_UNSTEADY_BLUNT_BODY :
    Grid_ptr = Grid_Unsteady_Blunt_Body(Grid_ptr,
					IP.Number_of_Blocks_Idir,
					IP.Number_of_Blocks_Jdir,
					IP.Blunt_Body_Radius,
					IP.Blunt_Body_Mach_Number,
					IP.Number_of_Cells_Idir,
					IP.Number_of_Cells_Jdir,
					IP.Number_of_Ghost_Cells);
    break;
  case GRID_RINGLEB_FLOW :
    Grid_ptr = Grid_Ringleb_Flow(Grid_ptr,
				 IP.Number_of_Blocks_Idir,
				 IP.Number_of_Blocks_Jdir,
				 IP.Inner_Streamline_Number,
				 IP.Outer_Streamline_Number,
				 IP.Isotach_Line,
				 IP.Number_of_Cells_Idir,
				 IP.Number_of_Cells_Jdir,
				 IP.Number_of_Ghost_Cells);
    break;
  case GRID_BUMP_CHANNEL_FLOW :
    Grid_ptr = Grid_Bump_Channel_Flow(Grid_ptr,
				      IP.Number_of_Blocks_Idir,
				      IP.Number_of_Blocks_Jdir,
				      IP.Smooth_Bump,
				      IP.Number_of_Cells_Idir,
				      IP.Number_of_Cells_Jdir,
				      IP.Number_of_Ghost_Cells);
    break;
  case GRID_DRIVEN_CAVITY_FLOW :
    Grid_ptr = Grid_Driven_Cavity_Flow(Grid_ptr,
				       IP.Number_of_Blocks_Idir,
				       IP.Number_of_Blocks_Jdir,
				       IP.Box_Width,
				       IP.Box_Height,
				       IP.Mesh_Stretching_Type_Idir,
				       IP.Mesh_Stretching_Type_Jdir,
				       IP.Mesh_Stretching_Factor_Idir,
				       IP.Mesh_Stretching_Factor_Jdir,
				       IP.Number_of_Cells_Idir,
				       IP.Number_of_Cells_Jdir,
				       IP.Number_of_Ghost_Cells);
    break;
  case GRID_BACKWARD_FACING_STEP :
    Grid_ptr = Grid_Backward_Facing_Step(Grid_ptr,
					 IP.Number_of_Blocks_Idir,
					 IP.Number_of_Blocks_Jdir,
					 IP.Step_Height,
					 IP.Top_Wall_Deflection,
					 IP.Mesh_Stretching_Factor_Idir,
					 IP.Mesh_Stretching_Factor_Jdir,
					 IP.Number_of_Cells_Idir,
					 IP.Number_of_Cells_Jdir,
					 IP.Number_of_Ghost_Cells);
    break;
  case GRID_FORWARD_FACING_STEP :
    Grid_ptr = Grid_Forward_Facing_Step(Grid_ptr,
					 IP.Number_of_Blocks_Idir,
					 IP.Number_of_Blocks_Jdir,
					 IP.Step_Height,
					 IP.Channel_Gap,
					 IP.Mesh_Stretching_Factor_Idir,
					 IP.Mesh_Stretching_Factor_Jdir,
					 IP.Number_of_Cells_Idir,
					 IP.Number_of_Cells_Jdir,
					 IP.Number_of_Ghost_Cells);
    break;
 /* case GRID_MIXING_LAYER:
    Grid_ptr = Grid_Mixing_Layer(Grid_ptr,
				 IP.Number_of_Blocks_Idir,
				 IP.Number_of_Blocks_Jdir,
				 IP.Box_Width,
				 IP.Mesh_Stretching_Factor_Idir,
				 IP.Mesh_Stretching_Factor_Jdir,
				 IP.Number_of_Cells_Idir,
				 IP.Number_of_Cells_Jdir,
				 IP.Number_of_Ghost_Cells);
    break ;
    */    
  case GRID_NASA_ROTOR_37 :
    IP.Number_of_Blocks_Idir = 3;
    IP.Number_of_Blocks_Jdir = 2;
    if (!IP.Interface_IP.Component_List.Ni) {
      Grid_ptr = IP.NASA_Rotor37.genMeshH_3x2_AUTO(Grid_ptr,
						   IP.Rotor_Percent_Span, 
						   IP.NASA_Rotor37.z_up,
						   IP.NASA_Rotor37.z_d,
						   IP.Number_of_Cells_Idir, 
						   IP.Number_of_Cells_Jdir,
						   IP.Number_of_Ghost_Cells,
						   min(250,8*max(IP.Number_of_Cells_Idir,
								 IP.Number_of_Cells_Jdir)));
    } else {
      Grid_ptr = Grid_NASA_Rotor_37(Grid_ptr,
				    IP.Number_of_Blocks_Idir,
				    IP.Number_of_Blocks_Jdir,
				    IP.Rotor_Percent_Span,
				    IP.Number_of_Cells_Idir, 
				    IP.Number_of_Cells_Jdir,
				    IP.Number_of_Ghost_Cells);
    }
    break;
  case GRID_NASA_ROTOR_67 :
    IP.Number_of_Blocks_Idir = 3;
    IP.Number_of_Blocks_Jdir = 2;
    if (!IP.Interface_IP.Component_List.Ni) {
      Grid_ptr = IP.NASA_Rotor67.genMeshH_3x2_AUTO(Grid_ptr,
						   IP.Rotor_Percent_Span, 
						   IP.NASA_Rotor67.z_up,
						   IP.NASA_Rotor67.z_d,
						   IP.Number_of_Cells_Idir, 
						   IP.Number_of_Cells_Jdir,
						   IP.Number_of_Ghost_Cells,
						   min(250,8*max(IP.Number_of_Cells_Idir,
								 IP.Number_of_Cells_Jdir)));
    } else {
      Grid_ptr = Grid_NASA_Rotor_67(Grid_ptr,
				    IP.Number_of_Blocks_Idir,
				    IP.Number_of_Blocks_Jdir,
				    IP.Rotor_Percent_Span,
				    IP.Number_of_Cells_Idir, 
				    IP.Number_of_Cells_Jdir,
				    IP.Number_of_Ghost_Cells);
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

  if (IP.Interface_IP.Component_List.Ni > 0) {
    if (IP.Interface_IP.Component_List[1].Type == INTERFACE_FLAT_PLATE) {
      Translate_Multi_Block_Grid(Grid_ptr,
				 IP.Number_of_Blocks_Idir,
				 IP.Number_of_Blocks_Jdir,
				 Vector2D(HALF*IP.Interface_IP.Component_List[1].Length1,ZERO));
    }
  }

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
							HighTemp2D_Input_Parameters &IP) {

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
                                      HighTemp2D_Input_Parameters &IP) {

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
							      HighTemp2D_Input_Parameters &IP) {

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
                           HighTemp2D_Input_Parameters &IP) {

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
						   HighTemp2D_Input_Parameters &IP) {

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
                   HighTemp2D_Input_Parameters &IP) {

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
                         HighTemp2D_Input_Parameters &IP) {

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
                         HighTemp2D_Input_Parameters &IP) {

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
