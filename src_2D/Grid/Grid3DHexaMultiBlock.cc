/* Grid3DHexaMultiBlock.cc:  Multi-block subroutines for 
                             3D hexahedral block grid class. */

/* Include 3D quadrilateral block grid type header file. */

#ifndef _GRID3D_HEXA_BLOCK_INCLUDED
#include "Grid3DHexa.h"
#endif // _GRID3D_HEXA_BLOCK_INCLUDED

/*************************************************************************
 * Grid3D_Hexa_Block -- External subroutines for 3D array of grid blocks.*
 *************************************************************************/

/********************************************************
 * Routine: Allocate_Multi_Block_Grid                   *
 *                                                      *
 * Allocate memory for a 3D array of 3D hexahedral      *
 * multi-block grids.                                   *
 *                                                      *
 ********************************************************/
Grid3D_Hexa_Block***  Allocate_Multi_Block_Grid(Grid3D_Hexa_Block ***Grid_ptr,
						const int Number_of_Blocks_Idir,
						const int Number_of_Blocks_Jdir,
						const int Number_of_Blocks_Kdir) {
    
    /* Allocate memory. */
    
    Grid_ptr = new Grid3D_Hexa_Block**[Number_of_Blocks_Idir];
    for (int i = 0 ; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      Grid_ptr[i] = new Grid3D_Hexa_Block*[Number_of_Blocks_Jdir];
      for (int j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
	Grid_ptr[i][j] = new Grid3D_Hexa_Block[Number_of_Blocks_Kdir];
      }  /* endfor */
    }/* endfor */
    
    /* Return memory location. */
    return(Grid_ptr);
    
}

/********************************************************
 * Routine: Deallocate_Multi_Block_Grid                 *
 *                                                      *
 * Deallocate memory for a 3D array of 3D hexahedral    *
 * multi-block grids.                                   *
 *                                                      *
 ********************************************************/
Grid3D_Hexa_Block*** Deallocate_Multi_Block_Grid(Grid3D_Hexa_Block ***Grid_ptr,
						 const int Number_of_Blocks_Idir,
						 const int Number_of_Blocks_Jdir,
						 const int Number_of_Blocks_Kdir) {
 
    /* Deallocate memory. */

    for (int i = 0 ; i <= Number_of_Blocks_Idir-1 ; ++i ) {
       for ( int j=0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
	 for ( int k = Number_of_Blocks_Kdir-1 ; k >= 0; --k ) {
          if (Grid_ptr[i][j][k].Node != NULL && Grid_ptr[i][j][k].Cell != NULL) { 
             Grid_ptr[i][j][k].deallocate();
          } else if (Grid_ptr[i][j][k].Node != NULL) {
             Grid_ptr[i][j][k].deallocateNodes();
          } else if (Grid_ptr[i][j][k].Cell != NULL) {
             Grid_ptr[i][j][k].deallocateCells();
          } /* endif */
	 }/* end for */
	 delete []Grid_ptr[i][j];
	 Grid_ptr[i][j]=NULL;
       }  /* endfor */
       
       delete []Grid_ptr[i];
       Grid_ptr[i] = NULL;
    }  /* endfor */
    delete []Grid_ptr;
    Grid_ptr = NULL;

    /* Return memory location. */

    return(Grid_ptr);

}

/********************************************************
 * Routine: Write_Multi_Block_Grid                      *
 *                                                      *
 * Writes a 3D array of 3D hexahedral multi-block       *
 * grids to the specified output stream for retrieval   *
 * and re-use purposes.                                 *
 *                                                      *
 ********************************************************/
void Write_Multi_Block_Grid(Grid3D_Hexa_Block ***Grid_ptr,
			    const int Number_of_Blocks_Idir,
		            const int Number_of_Blocks_Jdir,
			    const int Number_of_Blocks_Kdir,
                            ostream &Out_File) {
 
  Out_File << Number_of_Blocks_Idir << " " 
	   << Number_of_Blocks_Jdir << " "
	   << Number_of_Blocks_Kdir << "\n";
  
  for (int  k = 0 ; k <= Number_of_Blocks_Kdir-1 ; ++k ) {
    for (int  j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
      for (int  i = 0 ; i <= Number_of_Blocks_Idir-1 ; ++i ) {
	Write_Hexa_Block(Grid_ptr[i][j][k], Out_File);
      }  /* endfor */
    }  /* endfor */
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
void Output_Tecplot(Grid3D_Hexa_Block ***Grid_ptr,
		    const int Number_of_Blocks_Idir,
		    const int Number_of_Blocks_Jdir,
		    const int Number_of_Blocks_Kdir,
	            ostream &Out_File) {

    int block_number, i_output_title;

    block_number = 0;
    i_output_title = 1;
    
    for (int  k = 0 ; k <= Number_of_Blocks_Kdir-1 ; ++k ) {
      for (int  j = 0; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
	for (int i = 0 ; i <= Number_of_Blocks_Idir-1 ; ++i ) {
          if (Grid_ptr[i][j][k].Node != NULL) {
	    Output_Tecplot(Grid_ptr[i][j][k],
			   block_number,
			   i_output_title,
			   Out_File);
	    block_number = block_number + 1;
             if (i_output_title) i_output_title = 0;
          } /* endif */
	}  /* endfor */
      }  /* endfor */
    }/* endfor */

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
void Output_Nodes_Tecplot(Grid3D_Hexa_Block ***Grid_ptr,
                          const int Number_of_Blocks_Idir,
                          const int Number_of_Blocks_Jdir,
			  const int Number_of_Blocks_Kdir,
	                  ostream &Out_File) {
  int  block_number = 0;
  int  i_output_title = 1;
  
  for (int  k = 0 ; k <= Number_of_Blocks_Kdir-1 ; ++k ) {
    for (int  j = 0; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
      for (int i = 0 ; i <= Number_of_Blocks_Idir-1 ; ++i ) {
	if (Grid_ptr[i][j][k].Node != NULL) {
	  Output_Nodes_Tecplot(Grid_ptr[i][j][k],
			       block_number,
			       i_output_title,
			       Out_File);
	  block_number = block_number + 1;
	  if (i_output_title) i_output_title = 0;
	} /* endif */
      }  /* endfor */
    }  /* endfor */
  }/* endfor */ 
}

/********************************************************
 * Routine: Output_Cells_Tecplot                        *
 *                                                      *
 * Writes the cells of a 2D array of 2D hexahedral      *
 * multi-block grids to the specified output stream in  *
 * a format suitable for plotting the grid with         *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
void Output_Cells_Tecplot(Grid3D_Hexa_Block ***Grid_ptr,
                          const int Number_of_Blocks_Idir,
                          const int Number_of_Blocks_Jdir,
			  const int Number_of_Blocks_Kdir,
 	                  ostream &Out_File) {

    int block_number, i_output_title;

    block_number = 0;
    i_output_title = 1;
  for (int  k = 0 ; k <= Number_of_Blocks_Kdir-1 ; ++k ) {
    for (int  j = 0; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
      for (int i = 0 ; i <= Number_of_Blocks_Idir-1 ; ++i ) {
	if (Grid_ptr[i][j][k].Node != NULL) {
	  Output_Cells_Tecplot(Grid_ptr[i][j][k],
			       block_number,
			       i_output_title,
			       Out_File);
	  block_number = block_number + 1;
	  if (i_output_title) i_output_title = 0;
	} /* endif */
      }  /* endfor */
    }  /* endfor */
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
void Output_Gnuplot(Grid3D_Hexa_Block ***Grid_ptr,
		    const int Number_of_Blocks_Idir,
		    const int Number_of_Blocks_Jdir,
		    const int Number_of_Blocks_Kdir,
	            ostream &Out_File) {
  
  int  block_number, i_output_title;

    block_number = 0;
    i_output_title = 1;
    
    for (int k = 0 ; k <= Number_of_Blocks_Kdir-1 ; ++k ) {
      for (int j = 0; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
	for(int i =0 ; i <= Number_of_Blocks_Kdir-1 ; ++i ) {
	  if (Grid_ptr[i][j][k].Node != NULL) {
	    Output_Gnuplot(Grid_ptr[i][j][k],
			  block_number,
			  i_output_title,
			  Out_File);
	   block_number = block_number + 1;
	   if (i_output_title) i_output_title = 0;
	  } /* endif */
	}  /* endfor */
      }  /* endfor */
    } /* endfor */
    
}

/********************************************************
 * Routine: Grid_Bluff_Body_Burner_3D                   *
 ********************************************************/
Grid3D_Hexa_Block*** Grid_Bluff_Body_Burner_3D(Grid3D_Hexa_Block ***Grid_ptr,
				               int &Number_of_Blocks_Idir,
				               int &Number_of_Blocks_Jdir,
                                               int &Number_of_Blocks_Kdir,
    				               const double &Radius_Fuel_Line,
                                               const double &Radius_Bluff_Body,
                                               const double &Radius_Coflow,
                                               const double &Length_Upstream,
                                               const double &Length_Downstream,
                                               const int Number_of_Cells_Idir,
                                               const int Number_of_Cells_Jdir,
                                               const int Number_of_Cells_Kdir,
                                               const int Number_of_Ghost_Cells) {

    int iBlk;
    int numblk_idir_fuel, numblk_jdir_fuel,
        numblk_idir_bluffbody, numblk_jdir_bluffbody,
        numblk_idir_coflow, numblk_jdir_coflow;
    Grid2D_Quad_Block **Grid2D_Fuel_Line_XYplane,
                      **Grid2D_Bluff_Body_Inner_XYplane,
                      **Grid2D_Bluff_Body_Outer_XYplane,
                      **Grid2D_Coflow_XYplane;

    /* Allocate required memory. */

    Number_of_Blocks_Idir = 42;
    Number_of_Blocks_Jdir = 1;
    Number_of_Blocks_Kdir = 1;

    Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
                                         Number_of_Blocks_Idir, 
                                         Number_of_Blocks_Jdir,
                                         Number_of_Blocks_Kdir);

    /* Creat 2D cross-section grids from which the 3D grid
       will be extruded. */

    Grid2D_Fuel_Line_XYplane = Grid_Tube_2D(
                                   Grid2D_Fuel_Line_XYplane,
                                   numblk_idir_fuel,
		                   numblk_jdir_fuel,
                                   Radius_Fuel_Line,
 		                   Number_of_Cells_Idir,
		                   Number_of_Cells_Jdir,
		                   Number_of_Ghost_Cells,
                                   STRETCHING_FCN_MAX_CLUSTERING,
                                   1.25);

    Grid2D_Bluff_Body_Inner_XYplane = Grid_Annulus_2D(
                                   Grid2D_Bluff_Body_Inner_XYplane,
                                   numblk_idir_bluffbody,
		                   numblk_jdir_bluffbody,
                                   Radius_Fuel_Line,
 			           HALF*Radius_Bluff_Body,
 		                   Number_of_Cells_Idir,
		                   Number_of_Cells_Jdir,
		                   Number_of_Ghost_Cells,
                                   STRETCHING_FCN_MIN_CLUSTERING,
                                   1.10);

    Grid2D_Bluff_Body_Outer_XYplane = Grid_Annulus_2D(
                                   Grid2D_Bluff_Body_Outer_XYplane,
                                   numblk_idir_bluffbody,
		                   numblk_jdir_bluffbody,
                                   HALF*Radius_Bluff_Body,
 			           Radius_Bluff_Body,
 		                   Number_of_Cells_Idir,
		                   Number_of_Cells_Jdir,
		                   Number_of_Ghost_Cells,
                                   STRETCHING_FCN_MAX_CLUSTERING,
                                   1.10);

    Grid2D_Coflow_XYplane = Grid_Annulus_2D(
                                   Grid2D_Coflow_XYplane,
                                   numblk_idir_coflow,
		                   numblk_jdir_coflow,
 			           Radius_Bluff_Body,
                                   Radius_Coflow,
 		                   Number_of_Cells_Idir,
		                   Number_of_Cells_Jdir,
		                   Number_of_Ghost_Cells,
                                   STRETCHING_FCN_MINMAX_CLUSTERING,
                                   1.10);

    /* Create the mesh for each block representing
       the complete grid. */

    for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {
        /* Create the splines defining the north, south,
           east, and west boundaries of the grid. */

        if (iBlk <= 4) {
           Extrude_Hexa_Block(Grid_ptr[iBlk][0][0],
                              Grid2D_Fuel_Line_XYplane[iBlk][0],
                              Number_of_Cells_Kdir,
			      STRETCHING_FCN_MIN_CLUSTERING,
                              1.25,
                              ZERO,
                              0.25*Length_Downstream);

        } else if (iBlk >= 5 && iBlk <= 9) {
           Extrude_Hexa_Block(Grid_ptr[iBlk][0][0],
                              Grid2D_Fuel_Line_XYplane[iBlk-5][0],
                              Number_of_Cells_Kdir,
			      STRETCHING_FCN_LINEAR,
                              ONE,
                              0.25*Length_Downstream,
                              Length_Downstream);

        } else if (iBlk >= 10 && iBlk <= 13) {
           Extrude_Hexa_Block(Grid_ptr[iBlk][0][0],
                              Grid2D_Bluff_Body_Inner_XYplane[iBlk-10][0],
                              Number_of_Cells_Kdir,
			      STRETCHING_FCN_MIN_CLUSTERING,
                              1.25,
                              ZERO,
                              0.25*Length_Downstream);

        } else if (iBlk >= 14 && iBlk <= 17) {
           Extrude_Hexa_Block(Grid_ptr[iBlk][0][0],
                              Grid2D_Bluff_Body_Inner_XYplane[iBlk-14][0],
                              Number_of_Cells_Kdir,
			      STRETCHING_FCN_LINEAR,
                              ONE,
                              0.25*Length_Downstream,
                              Length_Downstream);

        } else if (iBlk >= 18 && iBlk <= 21) {
           Extrude_Hexa_Block(Grid_ptr[iBlk][0][0],
                              Grid2D_Bluff_Body_Outer_XYplane[iBlk-18][0],
                              Number_of_Cells_Kdir,
			      STRETCHING_FCN_MIN_CLUSTERING,
                              1.25,
                              ZERO,
                              0.25*Length_Downstream);

        } else if (iBlk >= 22 && iBlk <= 25) {
           Extrude_Hexa_Block(Grid_ptr[iBlk][0][0],
                              Grid2D_Bluff_Body_Outer_XYplane[iBlk-22][0],
                              Number_of_Cells_Kdir,
			      STRETCHING_FCN_LINEAR,
                              ONE,
                              0.25*Length_Downstream,
                              Length_Downstream);

        } else if (iBlk >= 26 && iBlk <= 29) {
           Extrude_Hexa_Block(Grid_ptr[iBlk][0][0],
                              Grid2D_Coflow_XYplane[iBlk-26][0],
                              Number_of_Cells_Kdir,
			      STRETCHING_FCN_MIN_CLUSTERING,
                              1.25,
                              ZERO,
                              0.25*Length_Downstream);

        } else if (iBlk >= 30 && iBlk <= 33) {
           Extrude_Hexa_Block(Grid_ptr[iBlk][0][0],
                              Grid2D_Coflow_XYplane[iBlk-30][0],
                              Number_of_Cells_Kdir,
			      STRETCHING_FCN_LINEAR,
                              ONE,
                              0.25*Length_Downstream,
                              Length_Downstream);

        } else if (iBlk >= 34 && iBlk <= 37) {
           Extrude_Hexa_Block(Grid_ptr[iBlk][0][0],
                              Grid2D_Coflow_XYplane[iBlk-34][0],
                              Number_of_Cells_Kdir,
			      STRETCHING_FCN_MAX_CLUSTERING,
                              1.25,
                              -0.25*Length_Downstream,
                              ZERO);

        } else if (iBlk >= 38 && iBlk <= 41) {
           Extrude_Hexa_Block(Grid_ptr[iBlk][0][0],
                              Grid2D_Coflow_XYplane[iBlk-38][0],
                              Number_of_Cells_Kdir,
			      STRETCHING_FCN_LINEAR,
                              ONE,
                              -Length_Downstream,
                              -0.25*Length_Downstream);

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

    /* Return the grid. */

    return(Grid_ptr);

}
