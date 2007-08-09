/* Grid3DHexaMultiBlock.cc:  Member functions for 
                             3D hexahedral multiblock grid class. */

/* Include 3D hexahedral block grid header file. */

#ifndef GRID3D_HEXA_MULTIBLOCK_INCLUDED
#include "Grid3DHexaMultiBlock.h"
#endif //_GRID3D_HEXA_MULTIBLOCK_INCLUDED

/********************************************************
 * Routine: Allocate_Grid_Blocks                        *
 *                                                      *
 * Allocate memory for a 3D array of 3D hexahedral      *
 * multi-block grids.                                   *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Allocate(const int Ni, 
                                       const int Nj, 
                                       const int Nk) {

   assert( Ni >= 1 && Nj >= 1 && Nk >= 1 && !Allocated);

   NBlk_Idir = Ni; NBlk_Jdir = Nj; NBlk_Kdir = Nk; Allocated = 1;

   Grid_Blks = new Grid3D_Hexa_Block**[NBlk_Idir];
   for (int i = 0 ; i < NBlk_Idir ; ++i ) {
      Grid_Blks[i] = new Grid3D_Hexa_Block*[NBlk_Jdir];
      for (int j = 0 ; j <NBlk_Jdir ; ++j ) {
         Grid_Blks[i][j] = new Grid3D_Hexa_Block[NBlk_Kdir];
      }  /* endfor */
   }/* endfor */ 

   Allocated = 1;

}

/********************************************************
 * Routine: Deallocate_Multi_Block_Grid                 *
 *                                                      *
 * Deallocate memory for a 3D array of 3D hexahedral    *
 * multi-block grids.                                   *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Deallocate(void) {

   assert( NBlk_Idir >= 1 && NBlk_Jdir >= 1 && NBlk_Kdir >= 1 && Allocated); 

   for (int i = 0 ; i <= NBlk_Idir-1 ; ++i ) {
      for ( int j=0 ; j <= NBlk_Jdir-1 ; ++j ) {
	 for ( int k = NBlk_Kdir-1 ; k >= 0; --k ) {
            if (Grid_Blks[i][j][k].Used) Grid_Blks[i][j][k].deallocate();
	 }/* end for */
	 delete []Grid_Blks[i][j];
	 Grid_Blks[i][j]=NULL;
       }  /* endfor */
       
       delete []Grid_Blks[i];
       Grid_Blks[i] = NULL;
    }  /* endfor */

    delete []Grid_Blks;
    Grid_Blks = NULL;

    Allocated = 0;

}

/********************************************************
 * Routine: Copy                                        *
 *                                                      *
 * Make a copy of multiblock hexahedral grid Grid2.     *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Copy(Grid3D_Hexa_Multi_Block &Grid2) {

  if (Grid2.Allocated) {

    /* Ensure multiblock grid arrays have same dimensions. */

    if (Allocated && (NBlk_Idir != Grid2.NBlk_Idir ||
                      NBlk_Jdir != Grid2.NBlk_Jdir ||
                      NBlk_Kdir != Grid2.NBlk_Kdir) ) {
      Deallocate();
      Allocate(Grid2.NBlk_Idir, Grid2.NBlk_Jdir, Grid2.NBlk_Kdir);
    } else if (!Allocated) {
      Allocate(Grid2.NBlk_Idir, Grid2.NBlk_Jdir, Grid2.NBlk_Kdir);
    } /* endif */

    /* Copy each grid block. */

    for (int  k = 0 ; k < NBlk_Kdir ; ++k ) {
       for (int  j = 0 ; j < NBlk_Jdir ; ++j ) {
          for (int  i = 0 ; i < NBlk_Idir ; ++i ) {
	     if (Grid2.Grid_Blks[i][j][k].Used) 
                Grid_Blks[i][j][k].Copy(Grid2.Grid_Blks[i][j][k]);
          } /* endfor */
       } /* endfor */
    } /* endfor */

  } /* endif */

}

/********************************************************
 * Routine: Output                                      *
 *                                                      *
 * Writes a 3D array of 3D hexahedral multi-block       *
 * grids to the specified output stream for retrieval   *
 * and re-use purposes.                                 *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Output(ostream &Out_File) {
   
  Out_File << NBlk_Idir << " " 
           << NBlk_Jdir << " "
           << NBlk_Kdir << "\n";
  
  for (int  k = 0 ; k < NBlk_Kdir ; ++k ) {
     for (int  j = 0 ; j < NBlk_Jdir ; ++j ) {
       for (int  i = 0 ; i < NBlk_Idir ; ++i ) {
          Out_File << setprecision(14) << Grid_Blks[i][j][k] << setprecision(6);
       } /* endfor */
     } /* endfor */
  } /* endfor */

}

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the nodes of a 3D array of 3D hexahedral      *
 * multi-block grids to the specified output stream in  *
 * a format suitable for plotting the grid with         *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Output_Tecplot(ostream &Out_File) {

    int block_number, i_output_title;
    block_number = 0;
    i_output_title = 1;
    
    for (int  k = 0 ; k < NBlk_Kdir ; ++k ) {
       for (int  j = 0; j < NBlk_Jdir ; ++j ) {
	  for (int i = 0 ; i < NBlk_Idir ; ++i ) {
             if (Grid_Blks[i][j][k].Used) {
                Grid_Blks[i][j][k].Output_Tecplot(block_number,
                                                  i_output_title,
                                                  Out_File);
	        block_number = block_number + 1;
                if (i_output_title) i_output_title = 0;
             } /* endif */
	  } /* endfor */
       } /* endfor */
    } /* endfor */

}

/********************************************************
 * Routine: Output_Nodes_Tecplot                        *
 *                                                      *
 * Writes the nodes of a 3D array of 3D hexahedral      *
 * multi-block grids to the specified output stream in  *
 * a format suitable for plotting the grid with         *
 * TECPLOT.  Include boundary nodes.                    *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Output_Nodes_Tecplot(ostream &Out_File) {

    int block_number, i_output_title;
    block_number = 0;
    i_output_title = 1;
  
    for (int  k = 0 ; k < NBlk_Kdir ; ++k ) {
       for (int  j = 0; j < NBlk_Jdir ; ++j ) {
          for (int i = 0 ; i < NBlk_Idir ; ++i ) {
	     if (Grid_Blks[i][j][k].Used) {
                Grid_Blks[i][j][k].Output_Nodes_Tecplot(block_number,
                                                        i_output_title,
                                                        Out_File);
                block_number = block_number + 1;
                if (i_output_title) i_output_title = 0;
	     } /* endif */
          } /* endfor */
       } /* endfor */
    } /* endfor */ 

}

/********************************************************
 * Routine: Output_Cells_Tecplot                        *
 *                                                      *
 * Writes the cells of a 3D array of 3D hexahedral      *
 * multi-block grids to the specified output stream in  *
 * a format suitable for plotting the grid with         *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Output_Cells_Tecplot(ostream &Out_File) {

    int block_number, i_output_title;
    block_number = 0;
    i_output_title = 1;

    for (int  k = 0 ; k <NBlk_Kdir ; ++k ) {
       for (int  j = 0; j <NBlk_Jdir ; ++j ) {
          for (int i = 0 ; i <NBlk_Idir ; ++i ) {
	     if (Grid_Blks[i][j][k].Used) {
                Grid_Blks[i][j][k].Output_Cells_Tecplot(block_number,
                                                        i_output_title,
                                                        Out_File);
                block_number = block_number + 1;
                if (i_output_title) i_output_title = 0;
	     } /* endif */
          } /* endfor */
       } /* endfor */
    }/* endfor */ 

}

/********************************************************
 * Routine: Output_Gnuplot                              *
 *                                                      *
 * Writes the nodes of a 3D array of 3D hexahedral      *
 * multi-block grids to the specified output stream in  *
 * a format suitable for plotting the grid with         *
 * GNUPLOT.                                             *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Output_Gnuplot(ostream &Out_File) {
  
    int  block_number, i_output_title;
    block_number = 0;
    i_output_title = 1;
    
    for (int k = 0 ; k < NBlk_Kdir ; ++k ) {
        for (int j = 0; j < NBlk_Jdir ; ++j ) {
	   for (int i =0 ; i < NBlk_Idir ; ++i ) {
	      if (Grid_Blks[i][j][k].Used) {
	          Grid_Blks[i][j][k].Output_Gnuplot(block_number,
			                            i_output_title,
			                            Out_File);
                 block_number = block_number + 1;
	         if (i_output_title) i_output_title = 0;
	      } /* endif */
	   } /* endfor */
        } /* endfor */
    } /* endfor */
    
}

/********************************************************
 * Routine: Create_Grid                                 *
 *                                                      *
 * Generates a 3D multiblock mesh depending on          *
 * input parameters.                                    *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Create_Grid(Grid3D_Input_Parameters &Input){

    switch(Input.i_Grid) {
      case GRID_CUBE :
        Create_Grid_Cube(Input);
        break;

      case GRID_CHANNEL_XDIR :
      case GRID_CHANNEL_YDIR:
      case GRID_CHANNEL_ZDIR:
        Create_Grid_Channel(Input);
        break;

      case GRID_COUETTE_XDIR :
      case GRID_COUETTE_YDIR:
      case GRID_COUETTE_ZDIR:
        Create_Grid_Couette(Input);
        break;

      case GRID_PIPE :
        Create_Grid_Pipe(Input);
        break;

      case GRID_BLUFF_BODY_BURNER :
        Create_Grid_Bluff_Body_Burner(Input);
        break;

      default:
        Create_Grid_Cube(Input);
        break;
    } /* endswitch */

}

/********************************************************
 * Routine: Create_Grid_Cube                            *
 *                                                      *
 * Generates a 3D Cartesian multiblock mesh for a cube. *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Create_Grid_Cube(Grid3D_Input_Parameters &Input){

    int BC_top, BC_bottom;
    Grid2D_Quad_Block **Grid2D_Box_XYplane;

    /* Allocate required memory. */

    Allocate(Input.NBlk_Idir, Input.NBlk_Jdir, Input.NBlk_Kdir);

    /* Creat 2D cross-section grids from which the 3D grid
       will be extruded. */
    
    Grid2D_Box_XYplane = Grid_Rectangular_Box(Grid2D_Box_XYplane,
                                              Input.NBlk_Idir, 
                                              Input.NBlk_Jdir,
                                              Input.Box_Width,
                                              Input.Box_Height,
					      ON,
					      Input.Stretching_Type_Idir,
					      Input.Stretching_Type_Jdir,
					      Input.Stretching_Factor_Idir,
					      Input.Stretching_Factor_Jdir,
                                              Input.NCells_Idir,
                                              Input.NCells_Jdir,
					      Input.Nghost);


    /* Create the mesh for each block representing
       the complete grid. */

    for (int kBlk = 0; kBlk <= Input.NBlk_Kdir-1; ++kBlk) {
       for (int jBlk = 0; jBlk <= Input.NBlk_Jdir-1; ++jBlk) {
          for (int iBlk = 0; iBlk <= Input.NBlk_Idir-1; ++iBlk) {

             /* Extrude each of the grid blocks from the
                appropriate 2D grid in XY-plane. */

             Grid_Blks[iBlk][jBlk][kBlk].Extrude(Grid2D_Box_XYplane[iBlk][jBlk],
                                                 Input.NCells_Kdir,
			                         Input.Stretching_Type_Kdir,
						 Input.Stretching_Factor_Kdir,
                                                 -HALF*Input.Box_Length+
                                                 (double(kBlk)/double(Input.NBlk_Kdir))*Input.Box_Length,
                                                 -HALF*Input.Box_Length+
                                                 (double(kBlk+1)/double(Input.NBlk_Kdir))*Input.Box_Length);

             /* Assign top and bottom boundary conditions. */

	     if (kBlk == Input.NBlk_Kdir-1) {
                BC_top = BC_REFLECTION;
             } else {
                BC_top = BC_NONE;
             } /* endif */
             if (kBlk == 0) {
                BC_bottom = BC_REFLECTION;
             } else {
                BC_bottom = BC_NONE;
             } /* endif */

             Grid_Blks[iBlk][jBlk][kBlk].Set_BCs_Zdir(BC_top, BC_bottom);

	  } /* endfor */
       } /* endfor */
    } /* endfor */

    /* Deallocate 2D grid. */

    Grid2D_Box_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Box_XYplane,
                                                     Input.NBlk_Idir, 
                                                     Input.NBlk_Jdir);

}

/********************************************************
 * Routine: Create_Grid_Channel                         *
 *                                                      *
 * Generates a 3D Cartesian multiblock mesh for a       *
 * channel flow.                                        *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Create_Grid_Channel(Grid3D_Input_Parameters &Input){

    int BC_east, BC_west, BC_north, BC_south, BC_top, BC_bottom;
    Grid2D_Quad_Block **Grid2D_Box_XYplane;

    /* Allocate required memory. */

    Allocate(Input.NBlk_Idir, Input.NBlk_Jdir, Input.NBlk_Kdir);

    /* Creat 2D cross-section grids from which the 3D grid
       will be extruded. */
    
    Grid2D_Box_XYplane = Grid_Rectangular_Box(Grid2D_Box_XYplane,
                                              Input.NBlk_Idir, 
                                              Input.NBlk_Jdir,
                                              Input.Box_Width,
                                              Input.Box_Height,
					      ON,
					      Input.Stretching_Type_Idir,
					      Input.Stretching_Type_Jdir,
					      Input.Stretching_Factor_Idir,
					      Input.Stretching_Factor_Jdir,
                                              Input.NCells_Idir,
                                              Input.NCells_Jdir,
					      Input.Nghost);

    /* Create the mesh for each block representing
       the complete grid. */

    for (int kBlk = 0; kBlk <= Input.NBlk_Kdir-1; ++kBlk) {
       for (int jBlk = 0; jBlk <= Input.NBlk_Jdir-1; ++jBlk) {
          for (int iBlk = 0; iBlk <= Input.NBlk_Idir-1; ++iBlk) {

             /* Extrude each of the grid blocks from the
                appropriate 2D grid in XY-plane. */

             Grid_Blks[iBlk][jBlk][kBlk].Extrude(Grid2D_Box_XYplane[iBlk][jBlk],
                                                 Input.NCells_Kdir,
			                         Input.Stretching_Type_Kdir,
						 Input.Stretching_Factor_Kdir,
                                                 (double(kBlk)/double(Input.NBlk_Kdir))*Input.Box_Length,
                                                 (double(kBlk+1)/double(Input.NBlk_Kdir))*Input.Box_Length);

             /* Assign top and bottom boundary conditions. */

             if (Input.i_Grid == GRID_CHANNEL_ZDIR) {

   	        if (iBlk == Input.NBlk_Idir-1) {
                   BC_east = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_east = BC_NONE;
                } /* endif */
                if (iBlk == 0) {
                   BC_west = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_west = BC_NONE;
                } /* endif */

	        if (jBlk == Input.NBlk_Jdir-1) {
                   BC_north = BC_WALL_VISCOUS;
                } else {
                   BC_north = BC_NONE;
                } /* endif */
                if (jBlk == 0) {
                   BC_south = BC_WALL_VISCOUS;
                } else {
                   BC_south = BC_NONE;
                } /* endif */

	        if (kBlk == Input.NBlk_Kdir-1) {
                   BC_top = BC_CHANNEL_OUTFLOW;
                } else {
                   BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                   BC_bottom = BC_CHANNEL_INFLOW;
                } else {
                   BC_bottom = BC_NONE;
                } /* endif */

             } else if (Input.i_Grid == GRID_CHANNEL_XDIR) {

    	        if (iBlk == Input.NBlk_Idir-1) {
                   BC_east = BC_CHANNEL_OUTFLOW;
                } else {
                   BC_east = BC_NONE;
                } /* endif */
                if (iBlk == 0) {
                   BC_west = BC_CHANNEL_INFLOW;
                } else {
                   BC_west = BC_NONE;
                } /* endif */

	        if (jBlk == Input.NBlk_Jdir-1) {
                   BC_north = BC_WALL_VISCOUS;
                } else {
                   BC_north = BC_NONE;
                } /* endif */
                if (jBlk == 0) {
                   BC_south = BC_WALL_VISCOUS;
                } else {
                   BC_south = BC_NONE;
                } /* endif */

	        if (kBlk == Input.NBlk_Kdir-1) {
                   BC_top = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                   BC_bottom = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_bottom = BC_NONE;
                } /* endif */

             } else if (Input.i_Grid == GRID_CHANNEL_YDIR) {

   	        if (iBlk == Input.NBlk_Idir-1) {
                   BC_east = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_east = BC_NONE;
                } /* endif */
                if (iBlk == 0) {
                   BC_west = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_west = BC_NONE;
                } /* endif */

	        if (jBlk == Input.NBlk_Jdir-1) {
                   BC_north = BC_CHANNEL_OUTFLOW;
                } else {
                   BC_north = BC_NONE;
                } /* endif */
                if (jBlk == 0) {
                   BC_south = BC_CHANNEL_INFLOW;
                } else {
                   BC_south = BC_NONE;
                } /* endif */

	        if (kBlk == Input.NBlk_Kdir-1) {
                   BC_top = BC_WALL_VISCOUS;
                } else {
                   BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                   BC_bottom = BC_WALL_VISCOUS;
                } else {
                   BC_bottom = BC_NONE;
                } /* endif */
	     } /* endif */

             Grid_Blks[iBlk][jBlk][kBlk].Set_BCs(BC_east, 
                                                 BC_west, 
                                                 BC_north, 
                                                 BC_south, 
                                                 BC_top, 
                                                 BC_bottom);
 
	  } /* endfor */
       } /* endfor */
    } /* endfor */

    /* Deallocate 2D grid. */

    Grid2D_Box_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Box_XYplane,
                                                     Input.NBlk_Idir, 
                                                     Input.NBlk_Jdir);

}

/********************************************************
 * Routine: Create_Grid_Couette                         *
 *                                                      *
 * Generates a 3D Cartesian multiblock mesh for a       *
 * couette flow.                                        *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Create_Grid_Couette(Grid3D_Input_Parameters &Input){
      
    int BC_east, BC_west, BC_north, BC_south, BC_top, BC_bottom;
    Grid2D_Quad_Block **Grid2D_Box_XYplane;

    /* Allocate required memory. */

    Allocate(Input.NBlk_Idir, Input.NBlk_Jdir, Input.NBlk_Kdir);

    /* Creat 2D cross-section grids from which the 3D grid
       will be extruded. */
    
    Grid2D_Box_XYplane = Grid_Rectangular_Box(Grid2D_Box_XYplane,
                                              Input.NBlk_Idir, 
                                              Input.NBlk_Jdir,
                                              Input.Box_Width,
                                              Input.Box_Height,
					      ON,
					      Input.Stretching_Type_Idir,
					      Input.Stretching_Type_Jdir,
					      Input.Stretching_Factor_Idir,
					      Input.Stretching_Factor_Jdir,
                                              Input.NCells_Idir,
                                              Input.NCells_Jdir,
					      Input.Nghost);

    /* Create the mesh for each block representing
       the complete grid. */

    for (int kBlk = 0; kBlk <= Input.NBlk_Kdir-1; ++kBlk) {
       for (int jBlk = 0; jBlk <= Input.NBlk_Jdir-1; ++jBlk) {
          for (int iBlk = 0; iBlk <= Input.NBlk_Idir-1; ++iBlk) {

             /* Extrude each of the grid blocks from the
                appropriate 2D grid in XY-plane. */

             Grid_Blks[iBlk][jBlk][kBlk].Extrude(Grid2D_Box_XYplane[iBlk][jBlk],
                                                 Input.NCells_Kdir,
			                         Input.Stretching_Type_Kdir,
						 Input.Stretching_Factor_Kdir,
                                                 (double(kBlk)/double(Input.NBlk_Kdir))*Input.Box_Length,
                                                 (double(kBlk+1)/double(Input.NBlk_Kdir))*Input.Box_Length);

             /* Assign top and bottom boundary conditions. */

             if (Input.i_Grid == GRID_COUETTE_ZDIR) {

   	        if (iBlk == Input.NBlk_Idir-1) {
                   BC_east = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_east = BC_NONE;
                } /* endif */
                if (iBlk == 0) {
                   BC_west = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_west = BC_NONE;
                } /* endif */

	        if (jBlk == Input.NBlk_Jdir-1) {
		   BC_north = BC_MOVING_WALL;
                } else {
                   BC_north = BC_NONE;
                } /* endif */
                if (jBlk == 0) {
                   BC_south = BC_NO_SLIP;
                } else {
                   BC_south = BC_NONE;
                } /* endif */

	        if (kBlk == Input.NBlk_Kdir-1) {
                   BC_top = BC_FIXED_PRESSURE;
                } else {
                   BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                   BC_bottom = BC_FIXED_PRESSURE;
                } else {
                   BC_bottom = BC_NONE;
                } /* endif */

             } else if (Input.i_Grid == GRID_COUETTE_XDIR) {

   	        if (iBlk == Input.NBlk_Idir-1) {
                   BC_east = BC_FIXED_PRESSURE;
                } else {
                   BC_east = BC_NONE;
                } /* endif */
                if (iBlk == 0) {
                   BC_west = BC_FIXED_PRESSURE;
                } else {
                   BC_west = BC_NONE;
                } /* endif */

	        if (jBlk == Input.NBlk_Jdir-1) {
                   BC_north = BC_MOVING_WALL;
                } else {
                   BC_north = BC_NONE;
                } /* endif */
                if (jBlk == 0) {
                   BC_south = BC_NO_SLIP;
                } else {
                   BC_south = BC_NONE;
                } /* endif */

	        if (kBlk == Input.NBlk_Kdir-1) {
                   BC_top = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                   BC_bottom = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_bottom = BC_NONE;
                } /* endif */

             } else if (Input.i_Grid == GRID_COUETTE_YDIR) {

   	        if (iBlk == Input.NBlk_Idir-1) {
                   BC_east = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_east = BC_NONE;
                } /* endif */
                if (iBlk == 0) {
                   BC_west = BC_CONSTANT_EXTRAPOLATION;
                } else {
                   BC_west = BC_NONE;
                } /* endif */

	        if (jBlk == Input.NBlk_Jdir-1) {
                   BC_north = BC_FIXED_PRESSURE;
                } else {
                   BC_north = BC_NONE;
                } /* endif */
                if (jBlk == 0) {
                   BC_south = BC_FIXED_PRESSURE;
                } else {
                   BC_south = BC_NONE;
                } /* endif */

	        if (kBlk == Input.NBlk_Kdir-1) {
                   BC_top = BC_MOVING_WALL;
                } else {
                   BC_top = BC_NONE;
                } /* endif */
                if (kBlk == 0) {
                   BC_bottom = BC_NO_SLIP;
                } else {
                   BC_bottom = BC_NONE;
                } /* endif */
	     } /* endif */

             Grid_Blks[iBlk][jBlk][kBlk].Set_BCs(BC_east, 
                                                 BC_west, 
                                                 BC_north, 
                                                 BC_south, 
                                                 BC_top, 
                                                 BC_bottom);
 
	  } /* endfor */
       } /* endfor */
    } /* endfor */

    /* Deallocate 2D grid. */

    Grid2D_Box_XYplane = Deallocate_Multi_Block_Grid(Grid2D_Box_XYplane,
                                                     Input.NBlk_Idir, 
                                                     Input.NBlk_Jdir);

}

/********************************************************
 * Routine: Create_Grid_Pipe                            *
 *                                                      *
 * Generates a 3D multiblock mesh for pipe flows.       *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Create_Grid_Pipe(Grid3D_Input_Parameters &Input){

    int numblk_idir_pipe, numblk_jdir_pipe;
    Grid2D_Quad_Block **Grid2D_Pipe_XYplane;

    /* Allocate required memory. */

    Input.NBlk_Idir = 5;
    Input.NBlk_Jdir = 1;
    Input.NBlk_Kdir = 1;

    Allocate(Input.NBlk_Idir, Input.NBlk_Jdir, Input.NBlk_Kdir);

    /* Creat 2D cross-section grids from which the 3D grid
       will be extruded. */

    Grid2D_Pipe_XYplane = Grid_Tube_2D(
                                   Grid2D_Pipe_XYplane,
                                   numblk_idir_pipe,
		                   numblk_jdir_pipe,
                                   Input.Pipe_Radius,
                                   Input.NCells_Idir,
                                   Input.NCells_Jdir,
				   Input.Nghost,
                                   STRETCHING_FCN_MAX_CLUSTERING,
                                   1.25);

    /* Create the mesh for each block representing
       the complete grid. */

    for (int iBlk = 0; iBlk <= Input.NBlk_Idir-1; ++iBlk) {

        /* Extrude each of the grid blocks from the
           appropriate 2D grid in XY-plane. */

        Grid_Blks[iBlk][0][0].Extrude(Grid2D_Pipe_XYplane[iBlk][0],
                                      Input.NCells_Kdir,
       	  	                      Input.Stretching_Type_Kdir,
				      Input.Stretching_Factor_Kdir,
                                      ZERO,
                                      Input.Pipe_Length);
        if (iBlk == 0) {
           Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
         	                         BC_NONE,
                                         BC_NONE,
       	                                 BC_NONE,
                                         BC_NONE,
                                         BC_DIRICHLET);
        } else {
           Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
       	                                 BC_NONE,
					 BC_WALL_VISCOUS,
       	                                 BC_NONE,
                                         BC_NONE,
                                         BC_DIRICHLET);
        } /* endif */

    } /* endfor */

    /* Deallocate 2D grid. */

    Grid2D_Pipe_XYplane = Deallocate_Multi_Block_Grid(
                                   Grid2D_Pipe_XYplane,
                                   numblk_idir_pipe,
		                   numblk_jdir_pipe);

}

/********************************************************
 * Routine: Create_Grid_Bluff_Body_Burner               *
 *                                                      *
 * Generates a 3D multiblock mesh for TNF bluff body    *
 * burner.                                              *
 *                                                      *
 ********************************************************/
void Grid3D_Hexa_Multi_Block::Create_Grid_Bluff_Body_Burner(Grid3D_Input_Parameters &Input){

    int numblk_idir_fuel, numblk_jdir_fuel,
        numblk_idir_bluffbody, numblk_jdir_bluffbody,
        numblk_idir_coflow, numblk_jdir_coflow;
    Grid2D_Quad_Block **Grid2D_Fuel_Line_XYplane,
                      **Grid2D_Bluff_Body_Inner_XYplane,
                      **Grid2D_Bluff_Body_Outer_XYplane,
                      **Grid2D_Coflow_XYplane;

    /* Allocate required memory. */

    Input.NBlk_Idir = 42;
    Input.NBlk_Jdir = 1;
    Input.NBlk_Kdir = 1;

    Allocate(Input.NBlk_Idir, Input.NBlk_Jdir, Input.NBlk_Kdir);

    /* Creat 2D cross-section grids from which the 3D grid
       will be extruded. */

    Grid2D_Fuel_Line_XYplane = Grid_Tube_2D(
                                   Grid2D_Fuel_Line_XYplane,
                                   numblk_idir_fuel,
		                   numblk_jdir_fuel,
                                   Input.Radius_Fuel_Line,
                                   Input.NCells_Idir,
                                   Input.NCells_Jdir,
				   Input.Nghost,
                                   STRETCHING_FCN_MAX_CLUSTERING,
                                   1.25);

    Grid2D_Bluff_Body_Inner_XYplane = Grid_Annulus_2D(
                                   Grid2D_Bluff_Body_Inner_XYplane,
                                   numblk_idir_bluffbody,
		                   numblk_jdir_bluffbody,
                                   Input.Radius_Fuel_Line,
 			           HALF*Input.Radius_Bluff_Body,
                                   Input.NCells_Idir,
                                   Input.NCells_Jdir,
				   Input.Nghost,
                                   STRETCHING_FCN_MIN_CLUSTERING,
                                   1.10);

    Grid2D_Bluff_Body_Outer_XYplane = Grid_Annulus_2D(
                                   Grid2D_Bluff_Body_Outer_XYplane,
                                   numblk_idir_bluffbody,
		                   numblk_jdir_bluffbody,
                                   HALF*Input.Radius_Bluff_Body,
 			           Input.Radius_Bluff_Body,
                                   Input.NCells_Idir,
                                   Input.NCells_Jdir,
				   Input.Nghost,
                                   STRETCHING_FCN_MAX_CLUSTERING,
                                   1.10);

    Grid2D_Coflow_XYplane = Grid_Annulus_2D(
                                   Grid2D_Coflow_XYplane,
                                   numblk_idir_coflow,
		                   numblk_jdir_coflow,
 			           Input.Radius_Bluff_Body,
                                   Input.Radius_Coflow_Inlet_Pipe,
                                   Input.NCells_Idir,
                                   Input.NCells_Jdir,
				   Input.Nghost,
                                   STRETCHING_FCN_MIN_CLUSTERING,
                                   1.10);

    /* Create the mesh for each block representing
       the complete grid. */

    for (int iBlk = 0; iBlk <= Input.NBlk_Idir-1; ++iBlk) {

        /* Extrude each of the grid blocks from the
           appropriate 2D grid in XY-plane. */

        if (iBlk <= 4) {
           Grid_Blks[iBlk][0][0].Extrude(
                              Grid2D_Fuel_Line_XYplane[iBlk][0],
                              Input.NCells_Kdir,
			      STRETCHING_FCN_MIN_CLUSTERING,
                              1.25,
                              ZERO,
                              0.25*Input.Length_Combustor_Tube);
           if (iBlk == 0) {
	      Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                            BC_NONE,
                                            BC_NONE,
		                            BC_NONE,
                                            BC_NONE,
                                            BC_DIRICHLET);
	   } else {
              Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                            BC_NONE,
                                            BC_NONE,
		                            BC_NONE,
                                            BC_NONE,
                                            BC_DIRICHLET);
           }  /* endif */

        } else if (iBlk >= 5 && iBlk <= 9) {
           Grid_Blks[iBlk][0][0].Extrude(
                              Grid2D_Fuel_Line_XYplane[iBlk-5][0],
                              Input.NCells_Kdir,
			      STRETCHING_FCN_LINEAR,
                              ZERO,
                              0.25*Input.Length_Combustor_Tube,
                              Input.Length_Combustor_Tube);
           if (iBlk-5 == 0) {
	      Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                            BC_NONE,
                                            BC_NONE,
		                            BC_NONE,
                                            BC_FIXED_PRESSURE,
                                            BC_NONE);
	   } else {
              Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                            BC_NONE,
                                            BC_NONE,
		                            BC_NONE,
                                            BC_FIXED_PRESSURE,
                                            BC_NONE);
           } /* endif */

        } else if (iBlk >= 10 && iBlk <= 13) {
           Grid_Blks[iBlk][0][0].Extrude(
                              Grid2D_Bluff_Body_Inner_XYplane[iBlk-10][0],
                              Input.NCells_Kdir,
			      STRETCHING_FCN_MIN_CLUSTERING,
                              1.25,
                              ZERO,
                              0.25*Input.Length_Combustor_Tube);
	   Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                         BC_NONE,
                                         BC_NONE,
		                         BC_NONE,
                                         BC_NONE,
                                         BC_WALL_VISCOUS);

        } else if (iBlk >= 14 && iBlk <= 17) {
           Grid_Blks[iBlk][0][0].Extrude(
                              Grid2D_Bluff_Body_Inner_XYplane[iBlk-14][0],
                              Input.NCells_Kdir,
			      STRETCHING_FCN_LINEAR,
                              ZERO,
                              0.25*Input.Length_Combustor_Tube,
                              Input.Length_Combustor_Tube);
	   Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                         BC_NONE,
                                         BC_NONE,
		                         BC_NONE,
                                         BC_FIXED_PRESSURE,
                                         BC_NONE);

        } else if (iBlk >= 18 && iBlk <= 21) {
           Grid_Blks[iBlk][0][0].Extrude(
                              Grid2D_Bluff_Body_Outer_XYplane[iBlk-18][0],
                              Input.NCells_Kdir,
			      STRETCHING_FCN_MIN_CLUSTERING,
                              1.25,
                              ZERO,
                              0.25*Input.Length_Combustor_Tube);
	   Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                         BC_NONE,
                                         BC_NONE,
		                         BC_NONE,
                                         BC_NONE,
                                         BC_WALL_VISCOUS);

        } else if (iBlk >= 22 && iBlk <= 25) {
           Grid_Blks[iBlk][0][0].Extrude(
                              Grid2D_Bluff_Body_Outer_XYplane[iBlk-22][0],
                              Input.NCells_Kdir,
			      STRETCHING_FCN_LINEAR,
                              ONE,
                              0.25*Input.Length_Combustor_Tube,
                              Input.Length_Combustor_Tube);
	   Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                         BC_NONE,
                                         BC_NONE,
		                         BC_NONE,
                                         BC_FIXED_PRESSURE,
                                         BC_NONE);

        } else if (iBlk >= 26 && iBlk <= 29) {
           Grid_Blks[iBlk][0][0].Extrude(
                              Grid2D_Coflow_XYplane[iBlk-26][0],
                              Input.NCells_Kdir,
			      STRETCHING_FCN_MIN_CLUSTERING,
                              1.25,
                              ZERO,
                              0.25*Input.Length_Combustor_Tube);
	   Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                         BC_NONE,
                                         BC_REFLECTION,
		                         BC_NONE,
                                         BC_NONE,
                                         BC_NONE);

        } else if (iBlk >= 30 && iBlk <= 33) {
           Grid_Blks[iBlk][0][0].Extrude(
                              Grid2D_Coflow_XYplane[iBlk-30][0],
                              Input.NCells_Kdir,
			      STRETCHING_FCN_LINEAR,
                              ZERO,
                              0.25*Input.Length_Combustor_Tube,
                              Input.Length_Combustor_Tube);
	   Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                         BC_NONE,
                                         BC_REFLECTION,
		                         BC_NONE,
                                         BC_FIXED_PRESSURE,
                                         BC_NONE);

        } else if (iBlk >= 34 && iBlk <= 37) {
           Grid_Blks[iBlk][0][0].Extrude(
                              Grid2D_Coflow_XYplane[iBlk-34][0],
                              Input.NCells_Kdir,
			      STRETCHING_FCN_MAX_CLUSTERING,
                              1.25,
                              -0.25*Input.Length_Coflow_Inlet_Pipe,
                              ZERO);
	   Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                         BC_NONE,
                                         BC_REFLECTION,
		                         BC_WALL_VISCOUS,
                                         BC_NONE,
                                         BC_NONE);

        } else if (iBlk >= 38 && iBlk <= 41) {
           Grid_Blks[iBlk][0][0].Extrude(
                              Grid2D_Coflow_XYplane[iBlk-38][0],
                              Input.NCells_Kdir,
			      STRETCHING_FCN_LINEAR,
                              ONE,
                              -Input.Length_Coflow_Inlet_Pipe,
                              -0.25*Input.Length_Coflow_Inlet_Pipe);
	   Grid_Blks[iBlk][0][0].Set_BCs(BC_NONE,
		                         BC_NONE,
                                         BC_REFLECTION,
		                         BC_WALL_VISCOUS,
                                         BC_DIRICHLET,
                                         BC_NONE);

        } /* endif */

    } /* endfor */

    /* Deallocate 2D grids. */

    Grid2D_Fuel_Line_XYplane = Deallocate_Multi_Block_Grid(
                                   Grid2D_Fuel_Line_XYplane,
                                   numblk_idir_fuel,
		                   numblk_jdir_fuel);

    Grid2D_Bluff_Body_Inner_XYplane = Deallocate_Multi_Block_Grid(
                                   Grid2D_Bluff_Body_Inner_XYplane,
                                   numblk_idir_bluffbody,
		                   numblk_jdir_bluffbody);

    Grid2D_Bluff_Body_Outer_XYplane = Deallocate_Multi_Block_Grid(
                                   Grid2D_Bluff_Body_Outer_XYplane,
                                   numblk_idir_bluffbody,
		                   numblk_jdir_bluffbody);

    Grid2D_Coflow_XYplane = Deallocate_Multi_Block_Grid(
                                   Grid2D_Coflow_XYplane,
                                   numblk_idir_coflow,
		                   numblk_jdir_coflow);

}
