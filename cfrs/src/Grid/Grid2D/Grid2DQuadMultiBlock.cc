/* Grid2DQuadMultiBlock.cc:  Multi-block subroutines for 
                             2D quadrilateral block grid class. */

/* Include 2D quadrilateral block grid type header file. */

#ifndef _GRID2D_QUAD_BLOCK_INCLUDED
#include "Grid2DQuad.h"
#endif // _GRID2D_QUAD_BLOCK_INCLUDED

/*************************************************************************
 * Grid2D_Quad_Block -- External subroutines for 2D array of grid blocks.*
 *************************************************************************/

/********************************************************
 * Routine: Allocate_Multi_Block_Grid                   *
 *                                                      *
 * Allocate memory for a 2D array of 2D quadrilateral   *
 * multi-block grids.                                   *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Allocate_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
					      const int Number_of_Blocks_Idir,
		                              const int Number_of_Blocks_Jdir) {

    int i;
 
    /* Allocate memory. */

    Grid_ptr = new Grid2D_Quad_Block*[Number_of_Blocks_Idir];
    for ( i = 0 ; i <= Number_of_Blocks_Idir-1 ; ++i ) {
      Grid_ptr[i] = new Grid2D_Quad_Block[Number_of_Blocks_Jdir];
    }  /* endfor */

    /* Return memory location. */

    return(Grid_ptr);

}

/********************************************************
 * Routine: Deallocate_Multi_Block_Grid                 *
 *                                                      *
 * Deallocate memory for a 2D array of 2D quadrilateral *
 * multi-block grids.                                   *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Deallocate_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
					        const int Number_of_Blocks_Idir,
		                                const int Number_of_Blocks_Jdir) {

    int i, j;
 
    /* Deallocate memory if memory was allocated. */

    if (Grid_ptr != NULL){
      for ( i = 0 ; i <= Number_of_Blocks_Idir-1 ; ++i ) {
	for ( j = Number_of_Blocks_Jdir-1 ; j >= 0 ; --j ) {
          if (Grid_ptr[i][j].Node != NULL && Grid_ptr[i][j].Cell != NULL) { 
	    Grid_ptr[i][j].deallocate();
          } else if (Grid_ptr[i][j].Node != NULL) {
	    Grid_ptr[i][j].deallocateNodes();
          } else if (Grid_ptr[i][j].Cell != NULL) {
	    Grid_ptr[i][j].deallocateCells();
          } /* endif */
	}  /* endfor */

	delete []Grid_ptr[i];
	Grid_ptr[i] = NULL;
      }  /* endfor */
      delete []Grid_ptr;
      Grid_ptr = NULL;
    }
    /* Return memory location. */
    return(Grid_ptr);
}

/********************************************************
 * Routine: Broadcast_Multi_Block_Grid                  *
 *                                                      *
 * Broadcasts a 2D array of 2D quadrilateral            *
 * multi-block grids to all processors involved in the  *
 * calculation from the primary processor using the MPI *
 * broadcast routine.                                   *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Broadcast_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
			                       int &Number_of_Blocks_Idir,
		                               int &Number_of_Blocks_Jdir) {
    
#ifdef _MPI_VERSION
    int i, j;

    if (!CFDkit_Primary_MPI_Processor()) {
       if (Grid_ptr != NULL) Grid_ptr = Deallocate_Multi_Block_Grid(Grid_ptr,
                                                                    Number_of_Blocks_Idir,
                                                                    Number_of_Blocks_Jdir);
       Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr,
                                            Number_of_Blocks_Idir,
                                            Number_of_Blocks_Jdir);
    } /* endif */

    for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
       for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
 	  Broadcast_Quad_Block(Grid_ptr[i][j]);
       }  /* endfor */
    }  /* endfor */
#endif

    return(Grid_ptr);

}

/********************************************************
 * Routine: Write_Multi_Block_Grid_Definition           *
 *                                                      *
 * Writes definition file information for a 2D array of *
 * 2D quadrilateral multi-block grids to the specified  *
 * output stream for retrieval and re-use purposes.     *
 *                                                      *
 ********************************************************/
void Write_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
			               const int Number_of_Blocks_Idir,
		                       const int Number_of_Blocks_Jdir,
                                       ostream &Out_File) {

    int i, j;
 
    Out_File << Number_of_Blocks_Idir << " " 
             << Number_of_Blocks_Jdir << "\n";

    for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
       for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
 	  Write_Quad_Block_Definition(Grid_ptr[i][j], Out_File);
       }  /* endfor */
    }  /* endfor */

}

/********************************************************
 * Routine: Read_Multi_Block_Grid_Definition            *
 *                                                      *
 * Reads definition file information for a 2D array of  *
 * 2D quadrilateral multi-block grids from the          *
 * specified input stream.                              *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Read_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
			                             int &Number_of_Blocks_Idir,
		                                     int &Number_of_Blocks_Jdir,
                                                     istream &In_File) {
    
    int i, j;

    In_File.setf(ios::skipws);
    In_File >> i >> j;
    In_File.unsetf(ios::skipws);
    
    if (i != Number_of_Blocks_Idir ||
        j != Number_of_Blocks_Jdir ||
        Grid_ptr == NULL) {
       if (Grid_ptr != NULL) Grid_ptr = Deallocate_Multi_Block_Grid(Grid_ptr,
                                                                    Number_of_Blocks_Idir,
                                                                    Number_of_Blocks_Jdir);
       Number_of_Blocks_Idir = i;
       Number_of_Blocks_Jdir = j;
       Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr,
                                            Number_of_Blocks_Idir,
                                            Number_of_Blocks_Jdir);
    } /* endif */

    for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
       for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
 	  Read_Quad_Block_Definition(Grid_ptr[i][j], In_File);

	  if (Grid_ptr[i][j].Node != NULL) {
	    Smooth_Quad_Block(Grid_ptr[i][j],
			      2*max(Grid_ptr[i][j].NCi, 
				    Grid_ptr[i][j].NCj));
	  } /* endif */

       }  /* endfor */
    }  /* endfor */

    return(Grid_ptr);

}

/********************************************************
 * Routine: Write_Multi_Block_Grid                      *
 *                                                      *
 * Writes a 2D array of 2D quadrilateral multi-block    *
 * grids to the specified output stream for retrieval   *
 * and re-use purposes.                                 *
 *                                                      *
 ********************************************************/
void Write_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
			    const int Number_of_Blocks_Idir,
		            const int Number_of_Blocks_Jdir,
                            ostream &Out_File) {

    int i, j;
 
    Out_File << Number_of_Blocks_Idir << " " 
             << Number_of_Blocks_Jdir << "\n";

    for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
       for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
 	  Write_Quad_Block(Grid_ptr[i][j], Out_File);
       }  /* endfor */
    }  /* endfor */

}

/********************************************************
 * Routine: Read_Multi_Block_Grid                       *
 *                                                      *
 * Writes a 2D array of 2D quadrilateral multi-block    *
 * grids from the specified input stream.               *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Read_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
			                  int &Number_of_Blocks_Idir,
		                          int &Number_of_Blocks_Jdir,
                                          istream &In_File) {
    
    int i, j;

    In_File.setf(ios::skipws);
    In_File >> i >> j;
    In_File.unsetf(ios::skipws);

    if (i != Number_of_Blocks_Idir ||
        j != Number_of_Blocks_Jdir ||
        Grid_ptr == NULL) {
       if (Grid_ptr != NULL) Grid_ptr = Deallocate_Multi_Block_Grid(Grid_ptr,
                                                                    Number_of_Blocks_Idir,
                                                                    Number_of_Blocks_Jdir);
       Number_of_Blocks_Idir = i;
       Number_of_Blocks_Jdir = j;
       Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr,
                                            Number_of_Blocks_Idir,
                                            Number_of_Blocks_Jdir);
    } /* endif */

    for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
       for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
 	  Read_Quad_Block(Grid_ptr[i][j], In_File);
       }  /* endfor */
    }  /* endfor */

    return(Grid_ptr);

}

/********************************************************
 * Routine: Translate_Multi_Block_Grid                  *
 *                                                      *
 * Translate the nodes of a 2D array of 2D              *
 * quadrilateral multi-block grids.                     *
 *                                                      *
 ********************************************************/
void Translate_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
		                const int Number_of_Blocks_Idir,
		                const int Number_of_Blocks_Jdir,
	      	                const Vector2D &V) {

    int i, j;

    for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
       for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
          if (Grid_ptr[i][j].Node != NULL) {
             Translate_Quad_Block(Grid_ptr[i][j],
	                          V);
          } /* endif */
       }  /* endfor */
    }  /* endfor */

}

/********************************************************
 * Routine: Scale_Multi_Block_Grid                      *
 *                                                      *
 * Scales 2D array of 2D quadrilateral multi-block      *
 * grids.                                               *
 *                                                      *
 ********************************************************/
void Scale_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
		            const int Number_of_Blocks_Idir,
		            const int Number_of_Blocks_Jdir,
	      	            const double &Scaling_Factor) {

    int i, j;

    for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
       for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
          if (Grid_ptr[i][j].Node != NULL) {
             Scale_Quad_Block(Grid_ptr[i][j],
	                      Scaling_Factor);
          } /* endif */
       }  /* endfor */
    }  /* endfor */

}

/********************************************************
 * Routine: Rotate_Multi_Block_Grid                     *
 *                                                      *
 * Rotates 2D array of 2D quadrilateral multi-block     *
 * grids.                                               *
 *                                                      *
 ********************************************************/
void Rotate_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
		             const int Number_of_Blocks_Idir,
		             const int Number_of_Blocks_Jdir,
	      	             const double &Angle) {

    int i, j;

    for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
       for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
          if (Grid_ptr[i][j].Node != NULL) {
             Rotate_Quad_Block(Grid_ptr[i][j],
	                       Angle);
          } /* endif */
       }  /* endfor */
    }  /* endfor */

}

/********************************************************
 * Routine: Reflect_Multi_Block_Grid                    *
 *                                                      *
 * Reflects 2D array of 2D quadrilateral multi-block    *
 * grids about y=0 axis.                                *
 *                                                      *
 ********************************************************/
void Reflect_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
		              const int Number_of_Blocks_Idir,
		              const int Number_of_Blocks_Jdir) {

    int i, j;

    for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
       for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
          if (Grid_ptr[i][j].Node != NULL) {
	    Reflect_Quad_Block(Grid_ptr[i][j]);
          } /* endif */
       }  /* endfor */
    }  /* endfor */

}

/********************************************************
 * Routine: Check_Multi_Block_Grid                      *
 *                                                      *
 * Check the validity of 2D array of 2D quadrilateral   *
 * multi-block grids.  Returns a non-zero result if     *
 * mesh is not valid.                                   *
 *                                                      *
 ********************************************************/
int Check_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
		            const int Number_of_Blocks_Idir,
		            const int Number_of_Blocks_Jdir) {

    int i, j;

    for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
       for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
          if (Grid_ptr[i][j].Node != NULL) {
	     if (Check_Quad_Block(Grid_ptr[i][j])) return(1);
          } /* endif */
       }  /* endfor */
    }  /* endfor */

    return(0);

}

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the nodes of a 2D array of 2D quadrilateral   *
 * multi-block grids to the specified output stream in  *
 * a format suitable for plotting the grid with         *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
void Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
		    const int Number_of_Blocks_Idir,
		    const int Number_of_Blocks_Jdir,
	            ostream &Out_File) {

    int i, j, block_number, i_output_title;

    block_number = 0;
    i_output_title = 1;

    for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
       for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
          if (Grid_ptr[i][j].Node != NULL) {
             Output_Tecplot(Grid_ptr[i][j],
                            block_number,
                            i_output_title,
	                    Out_File);
             block_number = block_number + 1;
             if (i_output_title) i_output_title = 0;
          } /* endif */
       }  /* endfor */
    }  /* endfor */

}

/********************************************************
 * Routine: Output_Nodes_Tecplot                        *
 *                                                      *
 * Writes the nodes of a 2D array of 2D quadrilateral   *
 * multi-block grids to the specified output stream in  *
 * a format suitable for plotting the grid with         *
 * TECPLOT.  Include boundary nodes.                    *
 *                                                      *
 ********************************************************/
void Output_Nodes_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                          const int Number_of_Blocks_Idir,
                          const int Number_of_Blocks_Jdir,
	                  ostream &Out_File) {

    int i, j, block_number, i_output_title;

    block_number = 0;
    i_output_title = 1;

    for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
       for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
          if (Grid_ptr[i][j].Node != NULL) {
             Output_Nodes_Tecplot(Grid_ptr[i][j],
                                  block_number,
                                  i_output_title,
	                          Out_File);
             block_number = block_number + 1;
             if (i_output_title) i_output_title = 0;
          } /* endif */
       }  /* endfor */
    }  /* endfor */

}

/********************************************************
 * Routine: Output_Cells_Tecplot                        *
 *                                                      *
 * Writes the cells of a 2D array of 2D quadrilateral   *
 * multi-block grids to the specified output stream in  *
 * a format suitable for plotting the grid with         *
 * TECPLOT.                                             *
 *                                                      *
 ********************************************************/
void Output_Cells_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                          const int Number_of_Blocks_Idir,
                          const int Number_of_Blocks_Jdir,
 	                  ostream &Out_File) {

    int i, j, block_number, i_output_title;

    block_number = 0;
    i_output_title = 1;

    for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
       for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
          if (Grid_ptr[i][j].Node != NULL) {
             Output_Cells_Tecplot(Grid_ptr[i][j],
                                  block_number,
                                  i_output_title,
	                          Out_File);
             block_number = block_number + 1;
             if (i_output_title) i_output_title = 0;
          } /* endif */
       }  /* endfor */
    }  /* endfor */

}

/********************************************************
 * Routine: Output_Gnuplot                              *
 *                                                      *
 * Writes the nodes of a 2D array of 2D quadrilateral   *
 * multi-block grids to the specified output stream in  *
 * a format suitable for plotting the grid with         *
 * GNUPLOT.                                             *
 *                                                      *
 ********************************************************/
void Output_Gnuplot(Grid2D_Quad_Block **Grid_ptr,
		    const int Number_of_Blocks_Idir,
		    const int Number_of_Blocks_Jdir,
	            ostream &Out_File) {

    int i, j, block_number, i_output_title;

    block_number = 0;
    i_output_title = 1;

    for ( j = 0 ; j <= Number_of_Blocks_Jdir-1 ; ++j ) {
       for ( i = 0; i <= Number_of_Blocks_Idir-1 ; ++i ) {
          if (Grid_ptr[i][j].Node != NULL) {
             Output_Gnuplot(Grid_ptr[i][j],
                            block_number,
                            i_output_title,
	                    Out_File);
             block_number = block_number + 1;
             if (i_output_title) i_output_title = 0;
          } /* endif */
       }  /* endfor */
    }  /* endfor */

}

/********************************************************
 * Routine: Grid_Rectangular_Box                        *
 *                                                      *
 * Generates a uniform 2D Cartesian mesh for a          *
 * rectangular box shaped domain.                       *
 *                                                      *
 * Usage: Grid_ptr = Grid_Rectangular_Box(Grid_ptr,     *
 *                                        nblk_i,       *
 *                                        nblk_j,       *
 *                                        TEN,          *
 *                                        FIVE,         *
 *   	                                  100,          *
 *  	                                  50,           *
 *                                        2);           *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Grid_Rectangular_Box(Grid2D_Quad_Block **Grid_ptr,
                                         int &Number_of_Blocks_Idir,
                                         int &Number_of_Blocks_Jdir,
                                         const double &Width,
                                         const double &Height,
 	                                 const int Number_of_Cells_Idir,
	                                 const int Number_of_Cells_Jdir,
					 const int Number_of_Ghost_Cells,
					 const int NumOfIter_UnsmoothMesh,
					 const int Stretch_I, const int Stretch_J,
					 const double Beta_I, const double Beta_J,
					 const double Tau_I, const double Tau_J) {

    int iBlk, jBlk, n_cells_i, n_cells_j, 
        Orthogonal_North, Orthogonal_South,
        Orthogonal_East, Orthogonal_West;
    Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
    Spline2D Bnd_Spline_North, Bnd_Spline_South,
             Bnd_Spline_East, Bnd_Spline_West;

    /* Allocate memory for grid block. */

    if (Number_of_Blocks_Idir < 0) Number_of_Blocks_Idir = 1;
    if (Number_of_Blocks_Jdir < 0) Number_of_Blocks_Jdir = 1;
    Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
                                         Number_of_Blocks_Idir, 
                                         Number_of_Blocks_Jdir);

    /* Create the mesh for each block representing
       the complete grid. */

    for ( jBlk = 0; jBlk <= Number_of_Blocks_Jdir-1; ++jBlk ) {
       for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {

          /* Assign values to the locations of the corners
             of the rectangular box shaped domain. */

          xc_NW = Vector2D(-HALF*Width+double(iBlk)*Width/double(Number_of_Blocks_Idir), 
                           -HALF*Height+double(jBlk+1)*Height/double(Number_of_Blocks_Jdir));
          xc_NE = Vector2D(-HALF*Width+double(iBlk+1)*Width/double(Number_of_Blocks_Idir), 
                           -HALF*Height+double(jBlk+1)*Height/double(Number_of_Blocks_Jdir));
          xc_SE = Vector2D(-HALF*Width+double(iBlk+1)*Width/double(Number_of_Blocks_Idir), 
                           -HALF*Height+double(jBlk)*Height/double(Number_of_Blocks_Jdir));
          xc_SW = Vector2D(-HALF*Width+double(iBlk)*Width/double(Number_of_Blocks_Idir), 
                           -HALF*Height+double(jBlk)*Height/double(Number_of_Blocks_Jdir));

          /* Create the splines defining the north, south,
             east, and west boundaries of the rectangular box. */

          Create_Spline_Line(Bnd_Spline_North, xc_NW, xc_NE, 2);
          Create_Spline_Line(Bnd_Spline_South, xc_SW, xc_SE, 2);
          Create_Spline_Line(Bnd_Spline_East, xc_SE, xc_NE, 2);
          Create_Spline_Line(Bnd_Spline_West, xc_SW, xc_NW, 2);

          /* Set the boundary condition types for each of the
             boundary splines. */

          if (jBlk == Number_of_Blocks_Jdir-1) {
	     Bnd_Spline_North.setBCtype(BC_REFLECTION);
          } else {
             Bnd_Spline_North.setBCtype(BC_NONE);
          } /* endif */
          if (jBlk == 0) {
	     Bnd_Spline_South.setBCtype(BC_REFLECTION);
          } else {
             Bnd_Spline_South.setBCtype(BC_NONE);
          } /* endif */
          if (iBlk == Number_of_Blocks_Idir-1) {
	     Bnd_Spline_East.setBCtype(BC_REFLECTION);
          } else {
             Bnd_Spline_East.setBCtype(BC_NONE);
          } /* endif */
          if (iBlk == 0) {
	     Bnd_Spline_West.setBCtype(BC_REFLECTION);
          } else {
             Bnd_Spline_West.setBCtype(BC_NONE);
          } /* endif */

          /* Assign values to the boundary grid line orthogonality parameters. */

          Orthogonal_North = 0;
          Orthogonal_South = 0;
          Orthogonal_East = 0;
          Orthogonal_West = 0;

          /* Create the 2D quadrilateral grid block representing
             the mesh. */

          Create_Quad_Block(Grid_ptr[iBlk][jBlk],
                            Bnd_Spline_North,
                            Bnd_Spline_South,
                            Bnd_Spline_East,
                            Bnd_Spline_West,
                            Number_of_Cells_Idir/Number_of_Blocks_Idir,
  	                    Number_of_Cells_Jdir/Number_of_Blocks_Jdir,
			    Number_of_Ghost_Cells,
                            GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
                            Stretch_I,
                            Beta_I, 
                            Tau_I,
                            Stretch_J,
                            Beta_J,
                            Tau_J,
                            Orthogonal_North,
	      	            Orthogonal_South,
     		            Orthogonal_East,
                            Orthogonal_West);

	  /* Smooth the 2D quadrilateral grid block. */
	  Smooth_Quad_Block(Grid_ptr[iBlk][jBlk],
			    min(250, 2*max(Number_of_Cells_Idir,Number_of_Cells_Jdir)));

	  /* Unsmooth the 2D quadrilateral grid block. */
	  Unsmooth_Interior_Nodes(Grid_ptr[iBlk][jBlk],NumOfIter_UnsmoothMesh);

          /* Deallocate the memory for the boundary splines. */

          Bnd_Spline_North.deallocate();
          Bnd_Spline_South.deallocate();
          Bnd_Spline_East.deallocate();
          Bnd_Spline_West.deallocate();

       } /* endfor */
    } /* endfor */

    /* Return the grid. */
    return(Grid_ptr);

}

/********************************************************
 * Routine: Grid_Flat_Plate                             *
 *                                                      *
 * Generates a quadilateral mesh with clustering        *
 * consisting of two grid blocks for predicting viscous *
 * flow and boundary layer development over a flat      *
 * plate.                                               *
 *                                                      *
 * Usage: Grid_ptr = Grid_Flat_Plate(Grid_ptr,          *
 *                                   nblk_i,            *
 *                                   nblk_j,            *
 *                                   TWO,               *
 *         		             100,               *
 *         		             100);              *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Grid_Flat_Plate(Grid2D_Quad_Block **Grid_ptr,
                                    int &Number_of_Blocks_Idir,
                                    int &Number_of_Blocks_Jdir,
                                    const double &Length,
 		                    const int Number_of_Cells_Idir,
		                    const int Number_of_Cells_Jdir,
				    const int Number_of_Ghost_Cells) {

    int iBlk, n_cells_i, n_cells_j, 
        Stretch_I, Stretch_J,
        Orthogonal_North, Orthogonal_South,
        Orthogonal_East, Orthogonal_West;
    double Beta_I, Tau_I, Beta_J, Tau_J;
    Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
    Spline2D Bnd_Spline_North, Bnd_Spline_South,
             Bnd_Spline_East, Bnd_Spline_West;

    /* Allocate memory for grid blocks.  There are two grid
       blocks for this mesh. */

    Number_of_Blocks_Idir = 2;
    Number_of_Blocks_Jdir = 1;
    Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
                                         Number_of_Blocks_Idir, 
                                         Number_of_Blocks_Jdir);

    /* Create the mesh for each block representing
       the complete grid. */

    for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {

        /* Assign values to the locations of the corners
           of the rectangular box shaped domain representing
           each of the blocks in the grid. */

        if (iBlk == 0) {
           xc_NW = Vector2D(-Length, TWO*Length);
           xc_NE = Vector2D(ZERO   , TWO*Length);
           xc_SE = Vector2D(ZERO   , ZERO);
           xc_SW = Vector2D(-Length, ZERO);
        } else {
           xc_NW = Vector2D(ZERO  , TWO*Length);
           xc_NE = Vector2D(Length, TWO*Length);
           xc_SE = Vector2D(Length, ZERO);
           xc_SW = Vector2D(ZERO  , ZERO);
        } /* endif */
   
        /* Create the splines defining the north, south,
           east, and west boundaries of the rectangular boxes. */

        Create_Spline_Line(Bnd_Spline_North, xc_NW, xc_NE, 2);
        Create_Spline_Line(Bnd_Spline_South, xc_SW, xc_SE, 2);
        Create_Spline_Line(Bnd_Spline_East, xc_SE, xc_NE, 2);
        Create_Spline_Line(Bnd_Spline_West, xc_SW, xc_NW, 2);

        /* Set the boundary condition types for each of the
           boundary splines. */

        if (iBlk == 0) {
           Bnd_Spline_North.setBCtype(BC_LINEAR_EXTRAPOLATION);
           Bnd_Spline_South.setBCtype(BC_REFLECTION);
           Bnd_Spline_East.setBCtype(BC_NONE);
           Bnd_Spline_West.setBCtype(BC_FIXED);
        } else {
           Bnd_Spline_North.setBCtype(BC_LINEAR_EXTRAPOLATION);
           Bnd_Spline_South.setBCtype(BC_FIXED_TEMP_WALL);
           Bnd_Spline_East.setBCtype(BC_LINEAR_EXTRAPOLATION);
           Bnd_Spline_West.setBCtype(BC_NONE);
        } /* endif */

        /* Determine the number of cells for this block. */

        n_cells_i = Number_of_Cells_Idir/2;
	n_cells_j = Number_of_Cells_Jdir;

        /* Assign values to the stretching function parameters
           and boundary grid line orthogonality parameters. */

        if (iBlk == 0) {
           Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
           Beta_I = 1.075; 
           Tau_I = ZERO;
           Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
           Beta_J = 1.01;
           Tau_J = ZERO;
           Orthogonal_North = 0;
           Orthogonal_South = 0;
           Orthogonal_East = 0;
           Orthogonal_West = 0;
        } else {
           Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
           Beta_I = 1.075; 
           Tau_I = ZERO;
           Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
           Beta_J = 1.01;
           Tau_J = ZERO;
           Orthogonal_North = 0;
           Orthogonal_South = 0;
           Orthogonal_East = 0;
           Orthogonal_West = 0;
        } /* endif */

        /* Create the 2D quadrilateral grid block. */

        Create_Quad_Block(Grid_ptr[iBlk][0],
                          Bnd_Spline_North,
                          Bnd_Spline_South,
                          Bnd_Spline_East,
                          Bnd_Spline_West,
                          n_cells_i,
  	                  n_cells_j,
			  Number_of_Ghost_Cells,
                          GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
                          Stretch_I,
                          Beta_I, 
                          Tau_I,
                          Stretch_J,
                          Beta_J,
                          Tau_J,
                          Orthogonal_North,
		          Orthogonal_South,
     		          Orthogonal_East,
                          Orthogonal_West);

        /* Deallocate the memory for the boundary splines. */

        Bnd_Spline_North.deallocate();
        Bnd_Spline_South.deallocate();
        Bnd_Spline_East.deallocate();
        Bnd_Spline_West.deallocate();

    } /* endfor */    

    /* Return the grid. */

    return(Grid_ptr);

}

/********************************************************
 * Routine: Grid_1D_Flame                               *
 *                                                      *
 * Generates a quadilateral mesh with clustering        *
 * consisting of one block clustered in the middle      *
 * for predicting 1D flame speeds of a particular       *
 * Chemical Mechanism.                                  *
 *                                                      *
 * Usage: Grid_ptr = Grid_1D_Flame (Grid_ptr,           *
 *                                   nblk_i,            *
 *                                   nblk_j,            *
 *                                   TWO,               *
 *                                   0.2,               *
 *         		             100,               *
 *         		             10,                *
 *                                   2);                *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Grid_1D_Flame(Grid2D_Quad_Block **Grid_ptr,
				  int &Number_of_Blocks_Idir,
				  int &Number_of_Blocks_Jdir,
				  const double &Length,
				  const double &Heigth,
				  const int Number_of_Cells_Idir,
				  const int Number_of_Cells_Jdir,
				  const int Number_of_Ghost_Cells) {

  int iBlk, n_cells_i, n_cells_j, 
    Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D Bnd_Spline_North, Bnd_Spline_South,
    Bnd_Spline_East, Bnd_Spline_West;
  
  /* Allocate memory for grid blocks.  There are two grid
     blocks for this mesh. */
  
  Number_of_Blocks_Idir = 1;
  Number_of_Blocks_Jdir = 1;
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
				       Number_of_Blocks_Idir, 
				       Number_of_Blocks_Jdir);
  
    /* Create the mesh for each block representing
       the complete grid. */

  for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {
      
    /* Assign values to the locations of the corners
       of the rectangular box shaped domain representing
       each of the blocks in the grid. */
    
    xc_NW = Vector2D(ZERO  , Heigth);
    xc_NE = Vector2D(Length/Number_of_Blocks_Idir, Heigth);
    xc_SE = Vector2D(Length/Number_of_Blocks_Idir, ZERO);
    xc_SW = Vector2D(ZERO  , ZERO);
    
    /* Create the splines defining the north, south,
       east, and west boundaries of the rectangular boxes. */
    
    Create_Spline_Line(Bnd_Spline_North, xc_NW, xc_NE, 2);
    Create_Spline_Line(Bnd_Spline_South, xc_SW, xc_SE, 2);
    Create_Spline_Line(Bnd_Spline_East, xc_SE, xc_NE, 2);
    Create_Spline_Line(Bnd_Spline_West, xc_SW, xc_NW, 2);
    
    /* Set the boundary condition types for each of the
       boundary splines. */

    Bnd_Spline_West.setBCtype(BC_FLAME_INFLOW);
    Bnd_Spline_North.setBCtype(BC_CONSTANT_EXTRAPOLATION);
    Bnd_Spline_South.setBCtype(BC_CONSTANT_EXTRAPOLATION);
    Bnd_Spline_East.setBCtype(BC_FLAME_OUTFLOW);

    /* Determine the number of cells for this block. */
    
    n_cells_i = Number_of_Cells_Idir; 
    n_cells_j = Number_of_Cells_Jdir;
    
    /* Assign values to the stretching function parameters
           and boundary grid line orthogonality parameters. */
    
    Stretch_J = STRETCHING_FCN_LINEAR;
    Beta_J = ZERO;
    Tau_J = ZERO;
    Orthogonal_North = 0;
    Orthogonal_South = 0;
    Orthogonal_East = 0;
    Orthogonal_West = 0;
    Stretch_I = STRETCHING_FCN_MIDPT_CLUSTERING;
    Beta_I=0.5;
    Tau_I= 8.0;
    //Stretch_I = STRETCHING_FCN_LINEAR;
    
    /* Create the 2D quadrilateral grid block. */
    
    Create_Quad_Block(Grid_ptr[iBlk][0],
		      Bnd_Spline_North,
		      Bnd_Spline_South,
		      Bnd_Spline_East,
		      Bnd_Spline_West,
		      n_cells_i,
		      n_cells_j,
		      Number_of_Ghost_Cells,
		      GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
		      Stretch_I,
		      Beta_I, 
		      Tau_I,
		      Stretch_J,
		      Beta_J,
		      Tau_J,
		      Orthogonal_North,
		      Orthogonal_South,
		      Orthogonal_East,
		      Orthogonal_West);
    
    /* Deallocate the memory for the boundary splines. */
    Bnd_Spline_North.deallocate();
    Bnd_Spline_South.deallocate();
    Bnd_Spline_East.deallocate();
    Bnd_Spline_West.deallocate();
    
  }

  /* Return the grid. */
  return(Grid_ptr);

}

/********************************************************
 * Routine: Grid_2D_Lamiar_Flame (Vertical)             *
 *                                                      *
 * Generates a quadilateral mesh with clustering        *
 * along the centerline (West) and entry (south)        *
 * for the predicition of 2D laminar diffusion flames   *                      
 *                                                      *
 * Usage: Grid_ptr = Grid_2D_Laminar_Flame (Grid_ptr,   *
 *                                   nblk_i,            *
 *                                   nblk_j,            *
 *                                   TWO,               *
 *                                   0.2,               *
 *         		             100,               *
 *         		             10,                *
 *                                   2);                *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Grid_2D_Laminar_Flame(Grid2D_Quad_Block **Grid_ptr,
					  int &Number_of_Blocks_Idir,
					  int &Number_of_Blocks_Jdir,
					  const double &Length,
					  const double &Heigth,
					  const int Number_of_Cells_Idir,
					  const int Number_of_Cells_Jdir,
					  const int Number_of_Ghost_Cells) {


    int  n_cells_i, n_cells_j, Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
    double Beta_I, Tau_I, Beta_J, Tau_J, Top, Bot, East,West;
    Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
    Spline2D Bnd_Spline_North, Bnd_Spline_South,
      Bnd_Spline_East, Bnd_Spline_West;

    //based on Heigth = 10cm and Length =5cm
    double fuel_spacing = 0.002;                   //m 
    double tube_spacing = fuel_spacing + 0.00038;  //m 
    double air_spacing = 0.025 - tube_spacing;     //m 

    //I-direction (inlet) blocks 
    int Number_of_Blocks_Fuel = 4;
    int Number_of_Blocks_Gap  = 1;
    int Number_of_Blocks_Air = 6;
    int Number_of_Blocks_Free = 1;

    if( Number_of_Blocks_Idir != Number_of_Blocks_Fuel + Number_of_Blocks_Gap + Number_of_Blocks_Air + Number_of_Blocks_Free ){
      cout<<"\n WARNING: Grid_2D_Laminar_Flame has a fixed initial number of Blocks in the x-direction to insure proper BC's, ";
      cout<<" currently it is set to "<< Number_of_Blocks_Fuel + Number_of_Blocks_Gap + Number_of_Blocks_Air + Number_of_Blocks_Free; 
    }

    Number_of_Blocks_Idir = Number_of_Blocks_Fuel + Number_of_Blocks_Gap + Number_of_Blocks_Air + Number_of_Blocks_Free;

    /* Allocate memory for grid blocks. */
    Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
                                         Number_of_Blocks_Idir, 
                                         Number_of_Blocks_Jdir);


    /* Create the mesh for each block representing the complete grid. */

    for ( int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++ ) {
      for ( int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++ ) {
	
	/* Assign values to the locations of the corners  of the rectangular box shaped domain representing
           each of the blocks in the grid. */
	       
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.05; 
	Tau_J = ZERO;
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.05; 
	Tau_I = ZERO;

	//Stretching for J Blocks
	Top = StretchingFcn(double(jBlk + 1)/double(Number_of_Blocks_Jdir), Beta_J, Tau_J, Stretch_J);
	Bot = StretchingFcn(double(jBlk)/double(Number_of_Blocks_Jdir), Beta_J, Tau_J, Stretch_J);  

	/**************** INLET **********************/
	// Fuel Inlet
	if(iBlk < Number_of_Blocks_Fuel ) {
	  xc_NW = Vector2D( double(iBlk) * fuel_spacing/double(Number_of_Blocks_Fuel), Top*Heigth);
	  xc_NE = Vector2D( double(iBlk + 1) * fuel_spacing/double(Number_of_Blocks_Fuel), Top*Heigth);			    
	  xc_SE = Vector2D( double(iBlk + 1) * fuel_spacing/double(Number_of_Blocks_Fuel), Bot*Heigth);			    
	  xc_SW = Vector2D( double(iBlk) * fuel_spacing/double(Number_of_Blocks_Fuel), Bot*Heigth);
			    
	  //Dead space
	} else if( iBlk == Number_of_Blocks_Gap + Number_of_Blocks_Fuel - 1 ) {     
	  xc_NW = Vector2D( fuel_spacing, Top*Heigth);
	  xc_NE = Vector2D( tube_spacing,Top*Heigth); 
	  xc_SE = Vector2D( tube_spacing, Bot*Heigth); 
	  xc_SW = Vector2D( fuel_spacing, Bot*Heigth); 

	  //Air Inlet
	} else if( iBlk < Number_of_Blocks_Gap + Number_of_Blocks_Fuel + Number_of_Blocks_Air) {

	  //Air Coflow Block Stretching I-direction
	  West = StretchingFcn(double(iBlk - (Number_of_Blocks_Fuel + Number_of_Blocks_Gap))/double(Number_of_Blocks_Air), Beta_I, Tau_I, Stretch_I);
	  East = StretchingFcn(double(iBlk- (Number_of_Blocks_Fuel + Number_of_Blocks_Gap) + 1)/double(Number_of_Blocks_Air), Beta_I, Tau_I, Stretch_I);  

	  xc_NW = Vector2D( tube_spacing + West*air_spacing, Top*Heigth);			    
	  xc_NE = Vector2D( tube_spacing + East*air_spacing, Top*Heigth);					 			   
	  xc_SE = Vector2D( tube_spacing + East*air_spacing, Bot*Heigth);					   
	  xc_SW = Vector2D( tube_spacing + West*air_spacing, Bot*Heigth);
			   	  
	  //Quiesent Air
	} else {
	  xc_NW = Vector2D( tube_spacing + air_spacing + double(iBlk - (Number_of_Blocks_Fuel + Number_of_Blocks_Gap + Number_of_Blocks_Air)) 
			    * (Length - air_spacing - tube_spacing)/double(Number_of_Blocks_Free), Top*Heigth);
	  xc_NE = Vector2D( tube_spacing + air_spacing + double(iBlk - (Number_of_Blocks_Fuel + Number_of_Blocks_Gap + Number_of_Blocks_Air) + 1) 
			    * (Length - air_spacing - tube_spacing)/double(Number_of_Blocks_Free), Top*Heigth);
	  xc_SE = Vector2D( tube_spacing + air_spacing + double(iBlk - (Number_of_Blocks_Fuel + Number_of_Blocks_Gap + Number_of_Blocks_Air) + 1) 
			    * (Length - air_spacing - tube_spacing)/double(Number_of_Blocks_Free), Bot*Heigth);	
	  xc_SW = Vector2D( tube_spacing + air_spacing + double(iBlk - (Number_of_Blocks_Fuel + Number_of_Blocks_Gap + Number_of_Blocks_Air ))
			    * (Length - air_spacing - tube_spacing)/double(Number_of_Blocks_Free), Bot*Heigth);	
	}
	
	Create_Spline_Line(Bnd_Spline_North, xc_NW, xc_NE, 2);
	Create_Spline_Line(Bnd_Spline_South, xc_SW, xc_SE, 2);
	Create_Spline_Line(Bnd_Spline_East, xc_SE, xc_NE, 2);
	Create_Spline_Line(Bnd_Spline_West, xc_SW, xc_NW, 2);
	
	/* Set the boundary condition types for each of the  boundary splines. */

	if (iBlk == 0 ) {
	  Bnd_Spline_West.setBCtype(BC_REFLECTION);    //centerline
	   Bnd_Spline_East.setBCtype(BC_NONE);
	} else if (iBlk == Number_of_Blocks_Idir-1 )  {
	  Bnd_Spline_West.setBCtype(BC_NONE);      
	  Bnd_Spline_East.setBCtype(BC_FREE_SLIP);    //farfield right  
	} else { 
	  Bnd_Spline_West.setBCtype(BC_NONE);      
	  Bnd_Spline_East.setBCtype(BC_NONE);
	}

	if (jBlk == 0 && iBlk == Number_of_Blocks_Gap + Number_of_Blocks_Fuel - 1 ) {
	  Bnd_Spline_South.setBCtype(BC_NO_SLIP);      //Gap wall
	  Bnd_Spline_North.setBCtype(BC_NONE);
	} else if (jBlk == 0 && iBlk < Number_of_Blocks_Gap + Number_of_Blocks_Fuel + Number_of_Blocks_Air ) {
	  Bnd_Spline_South.setBCtype(BC_FIXED);        //Bottom Inflow Left
	  Bnd_Spline_North.setBCtype(BC_NONE);
	} else if(jBlk == 0) {
	  Bnd_Spline_South.setBCtype(BC_NO_SLIP);      //Bottom Right        
	  Bnd_Spline_North.setBCtype(BC_NONE);
	} else if(jBlk == Number_of_Blocks_Jdir-1 ){  
	  Bnd_Spline_North.setBCtype(BC_CHARACTERISTIC); //Top outflow 
	  Bnd_Spline_South.setBCtype(BC_NONE);
	} else {
	  Bnd_Spline_South.setBCtype(BC_NONE);
	  Bnd_Spline_North.setBCtype(BC_NONE);
	}


      /* Assign values to the stretching function parameters
 	 and boundary grid line orthogonality parameters. */
	Stretch_I = STRETCHING_FCN_LINEAR; 
	Beta_I = ZERO; 
	Tau_I = ZERO;
	Orthogonal_East = 0;
	Orthogonal_West = 0;
	Stretch_J = STRETCHING_FCN_LINEAR;
	Beta_J = ZERO;
	Tau_J = ZERO; 
	Orthogonal_North = 0;
	Orthogonal_South = 0;
	
	/* Determine the number of cells for this block. */
	n_cells_i = Number_of_Cells_Idir/Number_of_Blocks_Idir;
	n_cells_j = Number_of_Cells_Jdir/Number_of_Blocks_Jdir;
	
        /* Create the 2D quadrilateral grid block. */
        Create_Quad_Block(Grid_ptr[iBlk][jBlk],
                          Bnd_Spline_North,
                          Bnd_Spline_South,
                          Bnd_Spline_East,
                          Bnd_Spline_West,
                          n_cells_i,
  	                  n_cells_j,
			  Number_of_Ghost_Cells,
                          GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
                          Stretch_I,
                          Beta_I, 
                          Tau_I,
                          Stretch_J,
                          Beta_J,
                          Tau_J,
                          Orthogonal_North,
		          Orthogonal_South,
     		          Orthogonal_East,
                          Orthogonal_West);

        /* Deallocate the memory for the boundary splines. */
        Bnd_Spline_North.deallocate();
        Bnd_Spline_South.deallocate();
        Bnd_Spline_East.deallocate();
        Bnd_Spline_West.deallocate();
      } 
    }
    /* Return the grid. */
    return(Grid_ptr);
}

// // Scott's Original Vertical Pipe
// /********************************************************
//  * Routine: Grid_2D_Lamiar_Flame                        *
//  *                                                      *
//  * Generates a quadilateral mesh with clustering        *
//  * along the centerline (West) and entry (south)        *
//  * for the predicition of 2D laminar diffusion flames   *                      
//  *                                                      *
//  * Usage: Grid_ptr = Grid_2D_Laminar_Flame (Grid_ptr,   *
//  *                                   nblk_i,            *
//  *                                   nblk_j,            *
//  *                                   TWO,               *
//  *                                   0.2,               *
//  *         		             100,               *
//  *         		             10);               *
//  *                                                      *
//  ********************************************************/

// Grid2D_Quad_Block** Grid_2D_Laminar_Flame(Grid2D_Quad_Block **Grid_ptr,
// 					  int &Number_of_Blocks_Idir,
// 					  int &Number_of_Blocks_Jdir,
// 					  const double &Length,
// 					  const double &Heigth,
// 					  const int Number_of_Cells_Idir,
// 					  const int Number_of_Cells_Jdir,
//                                        const int Number_of_Ghost_Cells) {


//   if( Number_of_Blocks_Idir !=2 ){
//     cout<<"\n WARNING: Grid_2D_Laminar_Flame must have only 2 I direction Blocks for proper BCs."; 
//   }

//     int  n_cells_i, n_cells_j, 
//         Stretch_I, Stretch_J,
//         Orthogonal_North, Orthogonal_South,
//         Orthogonal_East, Orthogonal_West;
//     double Beta_I, Tau_I, Beta_J, Tau_J;
//     Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
//     Spline2D Bnd_Spline_North, Bnd_Spline_South,
//              Bnd_Spline_East, Bnd_Spline_West;

//     /* Allocate memory for grid blocks. */
//     Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
//                                          Number_of_Blocks_Idir, 
//                                          Number_of_Blocks_Jdir);

//     /* Create the mesh for each block representing the complete grid. */

//     for ( int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++ ) {
//       for ( int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++ ) {

//         /* Assign values to the locations of the corners
//            of the rectangular box shaped domain representing
//            each of the blocks in the grid. */

// 	xc_NW = Vector2D( double(iBlk) * Length/double(Number_of_Blocks_Idir),
// 			  double(jBlk + 1) * Heigth/double(Number_of_Blocks_Jdir));
// 	xc_NE = Vector2D( double(iBlk + 1) * Length/double(Number_of_Blocks_Idir), 
// 			  double(jBlk + 1) * Heigth/double(Number_of_Blocks_Jdir));
// 	xc_SE = Vector2D( double(iBlk + 1) * Length/double(Number_of_Blocks_Idir), 
// 			  double(jBlk) * Heigth/double(Number_of_Blocks_Jdir));
// 	xc_SW = Vector2D( double(iBlk) * Length/double(Number_of_Blocks_Idir), 
// 			  double(jBlk) * Heigth/double(Number_of_Blocks_Jdir));
     
// 	/* Create the splines defining the north, south, 
// 	   east, and west boundaries of the rectangular boxes. */
	
// 	Create_Spline_Line(Bnd_Spline_North, xc_NW, xc_NE, 2);
// 	Create_Spline_Line(Bnd_Spline_South, xc_SW, xc_SE, 2);
// 	Create_Spline_Line(Bnd_Spline_East, xc_SE, xc_NE, 2);
// 	Create_Spline_Line(Bnd_Spline_West, xc_SW, xc_NW, 2);
	
// 	/* Set the boundary condition types for each of the
// 	   boundary splines. */

// 	if (iBlk == 0 ) {
// 	  Bnd_Spline_West.setBCtype(BC_REFLECTION);    //centerline
// 	   Bnd_Spline_East.setBCtype(BC_NONE);
// 	} else if (iBlk == Number_of_Blocks_Idir-1 )  {
// 	  Bnd_Spline_West.setBCtype(BC_NONE);      
// 	  Bnd_Spline_East.setBCtype(BC_FREE_SLIP);    //farfield right  
// 	} else { 
// 	  Bnd_Spline_West.setBCtype(BC_NONE);      
// 	  Bnd_Spline_East.setBCtype(BC_NONE);
// 	}

// 	if (jBlk == 0 && iBlk == 0) {
// 	  Bnd_Spline_South.setBCtype(BC_FIXED);        //Bottom Inflow Left
// 	  Bnd_Spline_North.setBCtype(BC_NONE);
// 	} else if(jBlk == 0 && iBlk != 0) {
// 	  Bnd_Spline_South.setBCtype(BC_NO_SLIP);      //Bottom Right        
// 	  Bnd_Spline_North.setBCtype(BC_NONE);
// 	} else if( jBlk == Number_of_Blocks_Jdir-1 ){  
// 	  Bnd_Spline_North.setBCtype(BC_CHARACTERISTIC); //Top outflow 
// 	  Bnd_Spline_South.setBCtype(BC_NONE);
// 	} else {
// 	  Bnd_Spline_South.setBCtype(BC_NONE);
// 	  Bnd_Spline_North.setBCtype(BC_NONE);
// 	}


//       /* Assign values to the stretching function parameters
//  	 and boundary grid line orthogonality parameters. */
// 	if (iBlk == 0) {
// 	  Stretch_I = STRETCHING_FCN_MIN_CLUSTERING; 
// 	  Beta_I = 1.05;  //1.075; 
// 	  Tau_I = ZERO;
// 	  Orthogonal_East = 0;
// 	  Orthogonal_West = 0;
//         } else {
// 	  Stretch_I = STRETCHING_FCN_LINEAR; 
// 	  Beta_I = ZERO; 
// 	  Tau_I = ZERO;
// 	  Orthogonal_East = 0;
// 	  Orthogonal_West = 0;
//  	}

// 	if (jBlk == 0) {
// 	  Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
// 	  Beta_J = 1.05; //1.075;
// 	  Tau_J = ZERO;
// 	  Orthogonal_North = 0;
// 	  Orthogonal_South = 0;
//         } else {
// 	  Stretch_J = STRETCHING_FCN_LINEAR;
// 	  Beta_J = ZERO;
// 	  Tau_J = ZERO; 
// 	  Orthogonal_North = 0;
// 	  Orthogonal_South = 0;
// 	}
     
	
// 	/* Determine the number of cells for this block. */
// 	n_cells_i = Number_of_Cells_Idir/Number_of_Blocks_Idir;
// 	n_cells_j = Number_of_Cells_Jdir/Number_of_Blocks_Jdir;
	
//         /* Create the 2D quadrilateral grid block. */
//         Create_Quad_Block(Grid_ptr[iBlk][jBlk],
//                           Bnd_Spline_North,
//                           Bnd_Spline_South,
//                           Bnd_Spline_East,
//                           Bnd_Spline_West,
//                           n_cells_i,
//   	                     n_cells_j,
//                           Number_of_Ghost_Cells,
//                           GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
//                           Stretch_I,
//                           Beta_I, 
//                           Tau_I,
//                           Stretch_J,
//                           Beta_J,
//                           Tau_J,
//                           Orthogonal_North,
// 		          Orthogonal_South,
//      		          Orthogonal_East,
//                           Orthogonal_West);

//         /* Deallocate the memory for the boundary splines. */
//         Bnd_Spline_North.deallocate();
//         Bnd_Spline_South.deallocate();
//         Bnd_Spline_East.deallocate();
//         Bnd_Spline_West.deallocate();
//       } 
//     }
//     /* Return the grid. */
//     return(Grid_ptr);
// }

/********************************************************
 * Routine: Grid_Pipe                                   *
 *                                                      *
 * Generates a single block quadilateral mesh with      *
 * clustering for predicting viscous flow and boundary  *
 * layer development in a cylindrical duct or pipe.     *
 *                                                      *
 * Usage: Grid_ptr = Grid_Pipe(Grid_ptr,                *
 *                             nblk_i,                  *
 *                             nblk_j,                  *
 *                             TEN,                     *
 *                             HALF,                    *
 *   	                       100,                     *
 *  	                       50,                      *
 *                             2);                      *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Grid_Pipe(Grid2D_Quad_Block **Grid_ptr,
                              int &Number_of_Blocks_Idir,
                              int &Number_of_Blocks_Jdir,
                              const double &Length,
                              const double &Radius,
			      const int Number_of_Cells_Idir,
	                      const int Number_of_Cells_Jdir,
			      const int Number_of_Ghost_Cells) {

  int Stretch_I, Stretch_J,
    Orthogonal_North, Orthogonal_South,
    Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D Bnd_Spline_North, Bnd_Spline_South,
    Bnd_Spline_East, Bnd_Spline_West;

  /* Allocate memory for grid block. */
  
  Number_of_Blocks_Idir = 1;
  Number_of_Blocks_Jdir = 1;
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
				       Number_of_Blocks_Idir, 
				       Number_of_Blocks_Jdir);
  
  /* Assign values to the locations of the corners
     of the rectangular box defining the pipe geometry. */
  xc_NW = Vector2D(ZERO , Radius);
  xc_NE = Vector2D(Length, Radius);
  xc_SE = Vector2D(Length, ZERO);
  xc_SW = Vector2D(ZERO , ZERO);
 
  /* Create the splines defining the north, south,
     east, and west boundaries of the grid. */
  
  Create_Spline_Line(Bnd_Spline_North, xc_NW, xc_NE, 2);
  Create_Spline_Line(Bnd_Spline_South, xc_SW, xc_SE, 2);
  Create_Spline_Line(Bnd_Spline_East, xc_SE, xc_NE, 2);
  Create_Spline_Line(Bnd_Spline_West, xc_SW, xc_NW, 2);
  
  /* Set the boundary condition types for each of the
     boundary splines. */
   
  Bnd_Spline_North.setBCtype(BC_NO_SLIP); //(BC_FIXED_TEMP_WALL);
  Bnd_Spline_South.setBCtype(BC_REFLECTION);
  Bnd_Spline_East.setBCtype(BC_FIXED);//(BC_LINEAR_EXTRAPOLATION); **Note LSR_2 xinfeng
  Bnd_Spline_West.setBCtype(BC_FIXED);

  /* Assign values to the stretching function parameters
     and boundary grid line orthogonality parameters. */
  Stretch_I = STRETCHING_FCN_LINEAR;
  Beta_I = ZERO; 
  Tau_I = ZERO;
  Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
  Beta_J = 1.0001;
  Tau_J = ZERO; 
    
  Orthogonal_North = 0;
  Orthogonal_South = 0;
  Orthogonal_East = 0;
  Orthogonal_West = 0;

  /* Create the 2D quadrilateral grid block representing
       the mesh. */
  Create_Quad_Block(Grid_ptr[0][0],
		    Bnd_Spline_North,
		    Bnd_Spline_South,
		    Bnd_Spline_East,
		    Bnd_Spline_West,
		    Number_of_Cells_Idir,
		    Number_of_Cells_Jdir,
		    Number_of_Ghost_Cells,
		    GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
		    Stretch_I,
		    Beta_I, 
		    Tau_I,
		    Stretch_J,
		    Beta_J,
		    Tau_J,
		    Orthogonal_North,
		    Orthogonal_South,
		    Orthogonal_East,
		    Orthogonal_West);
  
  /* Deallocate the memory for the boundary splines. */
  
  Bnd_Spline_North.deallocate();
  Bnd_Spline_South.deallocate();
  Bnd_Spline_East.deallocate();
  Bnd_Spline_West.deallocate();
  
  /* Return the grid. */
  
  return(Grid_ptr);
  
}

/********************************************************
 * Routine: Grid_Pipe                                   *
 *                                                      *
 * Generates a single block quadilateral mesh with      *
 * clustering for predicting viscous flow and boundary  *
 * layer development in a cylindrical duct or pipe.     *
 *                                                      *
 * Usage: Grid_ptr = Grid_Pipe(Grid_ptr,                *
 *                             nblk_i,                  *
 *                             nblk_j,                  *
 *                             TEN,                     *
 *                             HALF,                    *
 *   	                       100,                     *
 *  	                       50,                      *
 *                             2);                      *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Grid_Pipe(Grid2D_Quad_Block **Grid_ptr,
                              int &Number_of_Blocks_Idir,
                              int &Number_of_Blocks_Jdir,
                              const double &Length,
                              const double &Radius,
			      const int &Axisymmetric,
 	                      const int Number_of_Cells_Idir,
	                      const int Number_of_Cells_Jdir,
			      const int Number_of_Ghost_Cells) {

    int Stretch_I, Stretch_J,
        Orthogonal_North, Orthogonal_South,
        Orthogonal_East, Orthogonal_West;
    double Beta_I, Tau_I, Beta_J, Tau_J;
    Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
    Spline2D Bnd_Spline_North, Bnd_Spline_South,
             Bnd_Spline_East, Bnd_Spline_West;

    /* Allocate memory for grid block. */

    Number_of_Blocks_Idir = 1;
    Number_of_Blocks_Jdir = 1;
    Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
                                         Number_of_Blocks_Idir, 
                                         Number_of_Blocks_Jdir);

    /* Assign values to the locations of the corners
       of the rectangular box defining the pipe geometry. */

    if(Axisymmetric ==1){
      xc_NW = Vector2D(ZERO , Radius);
      xc_NE = Vector2D(Length, Radius);
      xc_SE = Vector2D(Length, ZERO);
      xc_SW = Vector2D(ZERO , ZERO);
    }
    if(Axisymmetric ==2){
      xc_NW = Vector2D(ZERO , Length);
      xc_NE = Vector2D(Radius,Length);
      xc_SE = Vector2D(Radius, ZERO);
      xc_SW = Vector2D(ZERO , ZERO); //axisymmetric 2
    }
    /* Create the splines defining the north, south,
       east, and west boundaries of the grid. */

    Create_Spline_Line(Bnd_Spline_North, xc_NW, xc_NE, 2);
    Create_Spline_Line(Bnd_Spline_South, xc_SW, xc_SE, 2);
    Create_Spline_Line(Bnd_Spline_East, xc_SE, xc_NE, 2);
    Create_Spline_Line(Bnd_Spline_West, xc_SW, xc_NW, 2);

    /* Set the boundary condition types for each of the
       boundary splines. */
    if(Axisymmetric ==1){
      Bnd_Spline_North.setBCtype(BC_NO_SLIP); //(BC_FIXED_TEMP_WALL);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_FIXED);//(BC_LINEAR_EXTRAPOLATION); **Note LSR_2 xinfeng
      Bnd_Spline_West.setBCtype(BC_FIXED);

      /* Assign values to the stretching function parameters
	 and boundary grid line orthogonality parameters. */
      
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ZERO; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
      Beta_J = 1.0001;
      Tau_J = ZERO;  //axisymmetric #1
      
    }
    if(Axisymmetric ==2){
      Bnd_Spline_East.setBCtype(BC_NO_SLIP); //(BC_FIXED_TEMP_WALL);
      Bnd_Spline_West.setBCtype(BC_REFLECTION);
      Bnd_Spline_North.setBCtype(BC_FIXED);//(BC_CONSTANT_EXTRAPOLATION); //(BC_LINEAR_EXTRAPOLATION); **Note LSR_2 xinfeng // axisymmetric 2
      Bnd_Spline_South.setBCtype(BC_FIXED);
        /* Assign values to the stretching function parameters
	 and boundary grid line orthogonality parameters. */
      Stretch_I =  STRETCHING_FCN_MAX_CLUSTERING;
      Beta_I = 1.0001;
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;
      Beta_J =  ZERO; 
      Tau_J = ZERO;
  
    }

    Orthogonal_North = 0;
    Orthogonal_South = 0;
    Orthogonal_East = 0;
    Orthogonal_West = 0;

    /* Create the 2D quadrilateral grid block representing
       the mesh. */

    Create_Quad_Block(Grid_ptr[0][0],
                      Bnd_Spline_North,
                      Bnd_Spline_South,
                      Bnd_Spline_East,
                      Bnd_Spline_West,
                      Number_of_Cells_Idir,
  	              Number_of_Cells_Jdir,
		      Number_of_Ghost_Cells,
                      GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
                      Stretch_I,
                      Beta_I, 
                      Tau_I,
                      Stretch_J,
                      Beta_J,
                      Tau_J,
                      Orthogonal_North,
		      Orthogonal_South,
     		      Orthogonal_East,
                      Orthogonal_West);

    /* Deallocate the memory for the boundary splines. */

    Bnd_Spline_North.deallocate();
    Bnd_Spline_South.deallocate();
    Bnd_Spline_East.deallocate();
    Bnd_Spline_West.deallocate();

    /* Return the grid. */

    return(Grid_ptr);

}

/********************************************************
 * Routine: Grid_Blunt_Body                             *
 *                                                      *
 * Generates a single block quadilateral mesh with      *
 * clustering for predicting supersonic flow around     *
 * a cirucular cylinder blunt body.                     *
 *                                                      *
 * Usage: Grid_ptr = Grid_Blunt_Body(Grid_ptr,          *
 *                                   nblk_i,            *
 *                                   nblk_j,            *
 *                                   ONE,               *
 *                                   FOUR,              *
 *   	                             150,               *
 *  	                             50,                *
 *                                   2);                *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Grid_Blunt_Body(Grid2D_Quad_Block **Grid_ptr,
                                    int &Number_of_Blocks_Idir,
                                    int &Number_of_Blocks_Jdir,
                                    const double &Radius,
                                    const double &Mach_Number,
 	                            const int Number_of_Cells_Idir,
	                            const int Number_of_Cells_Jdir,
				    const int Number_of_Ghost_Cells) {

    int Stretch_I, Stretch_J,
        Orthogonal_North, Orthogonal_South,
        Orthogonal_East, Orthogonal_West;
    double Beta_I, Tau_I, Beta_J, Tau_J;
    Vector2D x1, x2;
    Spline2D Bnd_Spline_North, Bnd_Spline_South,
             Bnd_Spline_East, Bnd_Spline_West;

    /* Allocate memory for grid block. */

    Number_of_Blocks_Idir = 1;
    Number_of_Blocks_Jdir = 1;
    Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
                                         Number_of_Blocks_Idir, 
                                         Number_of_Blocks_Jdir);

    /* Create the splines defining the north, south,
       east, and west boundaries of the grid. */

    Create_Spline_Bow_Shock(Bnd_Spline_North,
                            Radius,
			    Mach_Number,
			    1,
			    181);
    x1 = Vector2D(ZERO , ZERO);
    Create_Spline_Circular_Arc(Bnd_Spline_South,
	      		       x1,
			       Radius,
                               180.00,
			       90.00,
  	                       181);
    x1 = Bnd_Spline_South.Xp[0];
    x2 = Bnd_Spline_North.Xp[0];
    Create_Spline_Line(Bnd_Spline_West, x1, x2, 2);
    x1 = Bnd_Spline_South.Xp[Bnd_Spline_South.np-1];
    x2 = Bnd_Spline_North.Xp[Bnd_Spline_North.np-1];
    Create_Spline_Line(Bnd_Spline_East, x1, x2, 2);

    /* Set the boundary condition types for each of the
       boundary splines. */

    Bnd_Spline_North.setBCtype(BC_FIXED);
    Bnd_Spline_South.setBCtype(BC_REFLECTION);
    Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
    Bnd_Spline_West.setBCtype(BC_REFLECTION);

    /* Assign values to the stretching function parameters
       and boundary grid line orthogonality parameters. */

    Stretch_I = STRETCHING_FCN_LINEAR;
    Beta_I = ZERO; 
    Tau_I = ZERO;
    Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
    Beta_J = 1.05;
    Tau_J = ZERO;
    Orthogonal_North = 0;
    Orthogonal_South = 1;
    Orthogonal_East = 1;
    Orthogonal_West = 1;

    /* Create the 2D quadrilateral grid block representing
       the mesh. */

    Create_Quad_Block(Grid_ptr[0][0],
                      Bnd_Spline_North,
                      Bnd_Spline_South,
                      Bnd_Spline_East,
                      Bnd_Spline_West,
                      Number_of_Cells_Idir,
  	              Number_of_Cells_Jdir,
		      Number_of_Ghost_Cells,
                      GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
                      Stretch_I,
                      Beta_I, 
                      Tau_I,
                      Stretch_J,
                      Beta_J,
                      Tau_J,
                      Orthogonal_North,
		      Orthogonal_South,
     		      Orthogonal_East,
                      Orthogonal_West);

    /* Smooth the 2D quadrilateral grid block. */
    Smooth_Quad_Block(Grid_ptr[0][0],
                      min(250, 2*max(Number_of_Cells_Idir,Number_of_Cells_Jdir)));

    /* Deallocate the memory for the boundary splines. */

    Bnd_Spline_North.deallocate();
    Bnd_Spline_South.deallocate();
    Bnd_Spline_East.deallocate();
    Bnd_Spline_West.deallocate();

    /* Return the grid. */

    return(Grid_ptr);

}

/********************************************************
 * Routine: Grid_Rocket_Motor                           *
 *                                                      *
 * Generates a quadilateral mesh with clustering        *
 * consisting of two grid blocks for predicting the     *
 * axisymmetric core flow in a cylindrical grain solid  *
 * propellant rocket motor.                             *
 *                                                      *
 * Usage: Grid_ptr = Grid_Rocket_Motor(Grid_ptr,        *
 *                                     nblk_i,          *
 *                                     nblk_j,          *
 *                                     0.8350,          *
 *                                     0.0200,          *
 *                                     0.0500,          *
 *                                     0.1500,          *
 *                                     0.0300,          *
 *                                     0.0100,          *
 *         		               100,             *
 *         		               100,             *
 *                                     2);              *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Grid_Rocket_Motor(Grid2D_Quad_Block **Grid_ptr,
                                      int &Number_of_Blocks_Idir,
                                      int &Number_of_Blocks_Jdir,
                                      const double &Length_Grain,
                                      const double &Radius_Grain,
                                      const double &Length_Grain_To_Throat,
                                      const double &Length_Nozzle,
                                      const double &Radius_Nozzle_Exit,
                                      const double &Radius_Nozzle_Throat,
 		                      const int Number_of_Cells_Idir,
		                      const int Number_of_Cells_Jdir,
				      const int Number_of_Ghost_Cells) {

  int iBlk, n_cells_i, n_cells_j, 
      Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D Bnd_Spline_North, Bnd_Spline_South,
           Bnd_Spline_East, Bnd_Spline_West;
  double num1 = ZERO;
  double num2 = ZERO;

  /* Allocate memory for grid blocks.  There are two grid
     blocks for this mesh. */

  Number_of_Blocks_Idir = max(2, Number_of_Blocks_Idir);
  Number_of_Blocks_Jdir = 1;
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
				       Number_of_Blocks_Idir, 
				       Number_of_Blocks_Jdir);
  
  /* Create the mesh for each block representing
     the complete grid. */
  
  for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {
    
    num1 = Length_Grain*double(iBlk)/double(Number_of_Blocks_Idir-1);
    num2 = Length_Grain*double(iBlk+1)/double(Number_of_Blocks_Idir-1);

    /* Create the splines defining the north, south,
       east, and west boundaries of the grid block. */
    
    if (iBlk < Number_of_Blocks_Idir-1) {
      xc_NW = Vector2D(-Length_Grain + num1, Radius_Grain);
      xc_NE = Vector2D(-Length_Grain + num2, Radius_Grain);
      xc_SE = Vector2D(-Length_Grain + num2, ZERO);
      xc_SW = Vector2D(-Length_Grain + num1, ZERO);
      Create_Spline_Line(Bnd_Spline_North, xc_NW, xc_NE, 2);
      Create_Spline_Line(Bnd_Spline_South, xc_SW, xc_SE, 2);
      Create_Spline_Line(Bnd_Spline_East, xc_SE, xc_NE, 2);
      Create_Spline_Line(Bnd_Spline_West, xc_SW, xc_NW, 2);
    } else {
      xc_NW = Vector2D(ZERO         , Radius_Grain);
      xc_NE = Vector2D(Length_Nozzle, Radius_Nozzle_Exit);
      xc_SE = Vector2D(Length_Nozzle, ZERO);
      xc_SW = Vector2D(ZERO         , ZERO);
      Create_Spline_Line(Bnd_Spline_South, xc_SW, xc_SE, 2);
      Create_Spline_Line(Bnd_Spline_East, xc_SE, xc_NE, 2);
      Create_Spline_Line(Bnd_Spline_West, xc_SW, xc_NW, 2);
      Create_Spline_Area_Variation(Bnd_Spline_North,
				   ZERO,
				   Length_Grain_To_Throat,
				   Length_Nozzle,
				   Radius_Grain,
				   Radius_Nozzle_Throat,
				   Radius_Nozzle_Exit,
				   501);
    } /* endif */
    
      /* Set the boundary condition types for each of the
	 boundary splines. */
      
      if (iBlk < Number_of_Blocks_Idir-1) {
	if (iBlk == 0) {
	  Bnd_Spline_North.setBCtype(BC_BURNING_SURFACE);
	  Bnd_Spline_South.setBCtype(BC_REFLECTION);
	  Bnd_Spline_East.setBCtype(BC_NONE);
	  Bnd_Spline_West.setBCtype(BC_REFLECTION);
	} else {
	  Bnd_Spline_North.setBCtype(BC_BURNING_SURFACE);
	  Bnd_Spline_South.setBCtype(BC_REFLECTION);
	  Bnd_Spline_East.setBCtype(BC_NONE);
	  Bnd_Spline_West.setBCtype(BC_NONE);
	}
      } else {
	Bnd_Spline_North.setBCtype(BC_REFLECTION);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	Bnd_Spline_East.setBCtype(BC_LINEAR_EXTRAPOLATION);
	Bnd_Spline_West.setBCtype(BC_NONE);
      } /* endif */
      
      /* Determine the number of cells for this block. */
      
      if (iBlk < Number_of_Blocks_Idir-1) {
	n_cells_i = int(double(Number_of_Cells_Idir)*
			(Length_Grain/(Length_Grain+
				       Length_Nozzle)));
      } else {
	n_cells_i = Number_of_Cells_Idir -
	  int(double(Number_of_Cells_Idir)*
	      (Length_Grain/(Length_Grain+
			     Length_Nozzle)));
      } /* endif */
      n_cells_i = Number_of_Cells_Idir;
      n_cells_j = Number_of_Cells_Jdir;
      
      /* Assign values to the stretching function parameters
	 and boundary grid line orthogonality parameters. */
      
      if (iBlk < Number_of_Blocks_Idir-1) {
	Stretch_I = STRETCHING_FCN_LINEAR;
	Beta_I = ONE; 
	Tau_I = ZERO;
	Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_J = 1.05;
	Tau_J = ZERO;
	Orthogonal_North = 0;
	Orthogonal_South = 0;
	Orthogonal_East = 0;
	Orthogonal_West = 0;
      } else {
	Stretch_I = STRETCHING_FCN_MIDPT_CLUSTERING;
	Beta_I = Length_Grain_To_Throat/Length_Nozzle; 
	Tau_I = 3.00;
	Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_J = 1.05;
	Tau_J = ZERO;
	Orthogonal_North = 1;
	Orthogonal_South = 1;
	Orthogonal_East = 1;
	Orthogonal_West = 1;
      } /* endif */
      
      /* Create the 2D quadrilateral grid block. */
      
      Create_Quad_Block(Grid_ptr[iBlk][0],
			Bnd_Spline_North,
			Bnd_Spline_South,
			Bnd_Spline_East,
			Bnd_Spline_West,
			n_cells_i,
			n_cells_j,
			Number_of_Ghost_Cells,
			GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
			Stretch_I,
			Beta_I, 
			Tau_I,
			Stretch_J,
			Beta_J,
			Tau_J,
			Orthogonal_North,
			Orthogonal_South,
			Orthogonal_East,
			Orthogonal_West);
      
      /* Smooth the 2D quadrilateral grid block. */
      
      if (iBlk == Number_of_Blocks_Idir-1) 
	Smooth_Quad_Block(Grid_ptr[iBlk][0], min(250, 2*max(n_cells_i,n_cells_j)));
      
      /* Deallocate the memory for the boundary splines. */
      
      Bnd_Spline_North.deallocate();
      Bnd_Spline_South.deallocate();
      Bnd_Spline_East.deallocate();
      Bnd_Spline_West.deallocate();
      
    } /* endfor */

  /* Return the grid. */
  
  return(Grid_ptr);
  
} 

/**********************************************************************
 * Routine: Grid_Rocket_Motor_Cold_Flow                               *
 *                                                                    *
 * Generates a quadilateral mesh with clustering consisting of two    *
 * grid blocks for predicting the axisymmetric core flow geometry     *
 * given by Prof. Greatrix's cold flow aparatus for studing shock     *
 * reflections at the converging section of the nozzle.               *
 *                                                                    *
 * Usage: Grid_ptr = Grid_Rocket_Motor_Cold_Flow(Grid_ptr,            *
 *                                     ncells_i,                      *
 *                                     ncells_j,                      *
 *                                     Nghost_Cells);                 *
 *                                                                    *
 **********************************************************************/
Grid2D_Quad_Block** Grid_Rocket_Motor_Cold_Flow(Grid2D_Quad_Block **Grid_ptr,
						const int Number_of_Cells_Idir,
						const int Number_of_Cells_Jdir,
						const int Number_of_Ghost_Cells) {

  int n_cells_i, n_cells_j, 
      Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D Bnd_Spline_North, Bnd_Spline_South,
           Bnd_Spline_East, Bnd_Spline_West;
  double num1 = ZERO;
  double num2 = ZERO;
  int Number_of_Blocks_Idir = 14;
  int Number_of_Blocks_Jdir = 1;
  double Radius_Chamber = 0.0254;
  double Radius_Nozzle_Exit = 0.0105;
  double Radius_Nozzle_Throat = 0.00395;
  double Length_Main_Chamber = 2.0;
  double Length_High_Pressure_Chamber = 1.5;
  double Length_Converging_Section = 0.0215;
  double Length_Diverging_Section = 0.023;
  double Length_Nozzle_Throat = 0.001;

  // Allocate memory for grid blocks.  There are two grid blocks for 
  // this mesh.
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
				       Number_of_Blocks_Idir, 
				       Number_of_Blocks_Jdir);

  // Create the mesh for each block representing the complete grid.
  for (int iBlk = 0; iBlk <  Number_of_Blocks_Idir; iBlk++) {

    if (iBlk == 0 || iBlk == 1 || iBlk == 2) {
      // Define corner points.
      xc_NW = Vector2D(-3.5+double(iBlk)*HALF,Radius_Chamber);
      xc_NE = Vector2D(-3.0+double(iBlk)*HALF,Radius_Chamber);
      xc_SE = Vector2D(-3.0+double(iBlk)*HALF,ZERO);
      xc_SW = Vector2D(-3.5+double(iBlk)*HALF,ZERO);
      // Create the splines defining the north, south, east, and west 
      // boundaries of the grid block.
      Create_Spline_Line(Bnd_Spline_North, xc_NW, xc_NE, 2);
      Create_Spline_Line(Bnd_Spline_South, xc_SW, xc_SE, 2);
      Create_Spline_Line(Bnd_Spline_East,  xc_SE, xc_NE, 2);
      Create_Spline_Line(Bnd_Spline_West,  xc_SW, xc_NW, 2);
      // Set the boundary condition types for each of the boundary splines.
      Bnd_Spline_North.setBCtype(BC_REFLECTION);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      if (iBlk == 2) Bnd_Spline_East.setBCtype(BC_FIXED);
      else Bnd_Spline_East.setBCtype(BC_NONE);
      if (iBlk == 0) Bnd_Spline_West.setBCtype(BC_REFLECTION);
      else Bnd_Spline_West.setBCtype(BC_NONE);
      // Assign values to the stretching function parameters and boundary 
      // grid line orthogonality parameters.
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ONE; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
      Beta_J = 1.05;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;
      // Create the 2D quadrilateral grid block.
      Create_Quad_Block(Grid_ptr[iBlk][0],
			Bnd_Spline_North,
			Bnd_Spline_South,
			Bnd_Spline_East,
			Bnd_Spline_West,
			Number_of_Cells_Idir,
			Number_of_Cells_Jdir,
			Number_of_Ghost_Cells,
			GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
			Stretch_I,
			Beta_I, 
			Tau_I,
			Stretch_J,
			Beta_J,
			Tau_J,
			Orthogonal_North,
			Orthogonal_South,
			Orthogonal_East,
			Orthogonal_West);
      // Deallocate the memory for the boundary splines.
      Bnd_Spline_North.deallocate();
      Bnd_Spline_South.deallocate();
      Bnd_Spline_East.deallocate();
      Bnd_Spline_West.deallocate();

    } else if (iBlk == 3 || iBlk == 4 || iBlk == 5 || iBlk ==  6 ||
	       iBlk == 7 || iBlk == 8 || iBlk == 9 || iBlk == 10) {
      // Define corner points.
      xc_NW = Vector2D(-2.00+double(iBlk-3)*0.25,Radius_Chamber);
      xc_NE = Vector2D(-1.75+double(iBlk-3)*0.25,Radius_Chamber);
      xc_SE = Vector2D(-1.75+double(iBlk-3)*0.25,ZERO);
      xc_SW = Vector2D(-2.00+double(iBlk-3)*0.25,ZERO);
      // Create the splines defining the north, south, east, and west 
      // boundaries of the grid block.
      Create_Spline_Line(Bnd_Spline_North, xc_NW, xc_NE, 2);
      Create_Spline_Line(Bnd_Spline_South, xc_SW, xc_SE, 2);
      Create_Spline_Line(Bnd_Spline_East,  xc_SE, xc_NE, 2);
      Create_Spline_Line(Bnd_Spline_West,  xc_SW, xc_NW, 2);
      // Set the boundary condition types for each of the boundary splines.
      Bnd_Spline_North.setBCtype(BC_REFLECTION);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_NONE);
      if (iBlk == 3) Bnd_Spline_West.setBCtype(BC_FIXED);
      else Bnd_Spline_West.setBCtype(BC_NONE);
      // Assign values to the stretching function parameters and boundary 
      // grid line orthogonality parameters.
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ONE; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
      Beta_J = 1.05;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;
      // Create the 2D quadrilateral grid block.
      Create_Quad_Block(Grid_ptr[iBlk][0],
			Bnd_Spline_North,
			Bnd_Spline_South,
			Bnd_Spline_East,
			Bnd_Spline_West,
			Number_of_Cells_Idir,
			Number_of_Cells_Jdir,
			Number_of_Ghost_Cells,
			GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
			Stretch_I,
			Beta_I, 
			Tau_I,
			Stretch_J,
			Beta_J,
			Tau_J,
			Orthogonal_North,
			Orthogonal_South,
			Orthogonal_East,
			Orthogonal_West);
      // Deallocate the memory for the boundary splines.
      Bnd_Spline_North.deallocate();
      Bnd_Spline_South.deallocate();
      Bnd_Spline_East.deallocate();
      Bnd_Spline_West.deallocate();

    } else if (iBlk == 11) {
      // Define corner points.
      xc_NW = Vector2D(ZERO,Radius_Chamber);
      xc_NE = Vector2D(Length_Converging_Section,Radius_Nozzle_Throat);
      xc_SE = Vector2D(Length_Converging_Section,ZERO);
      xc_SW = Vector2D(ZERO,ZERO);
      // Create the splines defining the north, south, east, and west 
      // boundaries of the grid block.
      Create_Spline_Line(Bnd_Spline_North, xc_NW, xc_NE, 2);
      Create_Spline_Line(Bnd_Spline_South, xc_SW, xc_SE, 2);
      Create_Spline_Line(Bnd_Spline_East,  xc_SE, xc_NE, 2);
      Create_Spline_Line(Bnd_Spline_West,  xc_SW, xc_NW, 2);
      // Set the boundary condition types for each of the boundary splines.
      Bnd_Spline_North.setBCtype(BC_REFLECTION);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_NONE);
      // Assign values to the stretching function parameters and boundary 
      // grid line orthogonality parameters.
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ONE; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
      Beta_J = 1.05;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;
      // Create the 2D quadrilateral grid block.
      Create_Quad_Block(Grid_ptr[iBlk][0],
			Bnd_Spline_North,
			Bnd_Spline_South,
			Bnd_Spline_East,
			Bnd_Spline_West,
			Number_of_Cells_Idir,
			Number_of_Cells_Jdir,
			Number_of_Ghost_Cells,
			GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
			Stretch_I,
			Beta_I, 
			Tau_I,
			Stretch_J,
			Beta_J,
			Tau_J,
			Orthogonal_North,
			Orthogonal_South,
			Orthogonal_East,
			Orthogonal_West);
      // Deallocate the memory for the boundary splines.
      Bnd_Spline_North.deallocate();
      Bnd_Spline_South.deallocate();
      Bnd_Spline_East.deallocate();
      Bnd_Spline_West.deallocate();

    } else if (iBlk == 12) {
      // Define corner points.
      xc_NW = Vector2D(Length_Converging_Section,Radius_Nozzle_Throat);
      xc_NE = Vector2D(Length_Converging_Section+Length_Nozzle_Throat,Radius_Nozzle_Throat);
      xc_SE = Vector2D(Length_Converging_Section+Length_Nozzle_Throat,ZERO);
      xc_SW = Vector2D(Length_Converging_Section,ZERO);
      // Create the splines defining the north, south, east, and west 
      // boundaries of the grid block.
      Create_Spline_Line(Bnd_Spline_North, xc_NW, xc_NE, 2);
      Create_Spline_Line(Bnd_Spline_South, xc_SW, xc_SE, 2);
      Create_Spline_Line(Bnd_Spline_East,  xc_SE, xc_NE, 2);
      Create_Spline_Line(Bnd_Spline_West,  xc_SW, xc_NW, 2);
      // Set the boundary condition types for each of the boundary splines.
      Bnd_Spline_North.setBCtype(BC_REFLECTION);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_NONE);
      // Assign values to the stretching function parameters and boundary 
      // grid line orthogonality parameters.
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ONE; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
      Beta_J = 1.05;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;
      // Create the 2D quadrilateral grid block.
      Create_Quad_Block(Grid_ptr[iBlk][0],
			Bnd_Spline_North,
			Bnd_Spline_South,
			Bnd_Spline_East,
			Bnd_Spline_West,
			Number_of_Cells_Idir,
			Number_of_Cells_Jdir,
			Number_of_Ghost_Cells,
			GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
			Stretch_I,
			Beta_I, 
			Tau_I,
			Stretch_J,
			Beta_J,
			Tau_J,
			Orthogonal_North,
			Orthogonal_South,
			Orthogonal_East,
			Orthogonal_West);
      // Deallocate the memory for the boundary splines.
      Bnd_Spline_North.deallocate();
      Bnd_Spline_South.deallocate();
      Bnd_Spline_East.deallocate();
      Bnd_Spline_West.deallocate();

    } else if (iBlk == 13) {
      // Define corner points.
      xc_NW = Vector2D(Length_Converging_Section+Length_Nozzle_Throat,Radius_Nozzle_Throat);
      xc_NE = Vector2D(Length_Converging_Section+Length_Nozzle_Throat+Length_Diverging_Section,Radius_Nozzle_Exit);
      xc_SE = Vector2D(Length_Converging_Section+Length_Nozzle_Throat+Length_Diverging_Section,ZERO);
      xc_SW = Vector2D(Length_Converging_Section+Length_Nozzle_Throat,ZERO);
      // Create the splines defining the north, south, east, and west 
      // boundaries of the grid block.
      Create_Spline_Line(Bnd_Spline_North, xc_NW, xc_NE, 2);
      Create_Spline_Line(Bnd_Spline_South, xc_SW, xc_SE, 2);
      Create_Spline_Line(Bnd_Spline_East,  xc_SE, xc_NE, 2);
      Create_Spline_Line(Bnd_Spline_West,  xc_SW, xc_NW, 2);
      // Set the boundary condition types for each of the boundary splines.
      Bnd_Spline_North.setBCtype(BC_REFLECTION);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_West.setBCtype(BC_NONE);
      // Assign values to the stretching function parameters and boundary 
      // grid line orthogonality parameters.
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ONE; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
      Beta_J = 1.05;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;
      // Create the 2D quadrilateral grid block.
      Create_Quad_Block(Grid_ptr[iBlk][0],
			Bnd_Spline_North,
			Bnd_Spline_South,
			Bnd_Spline_East,
			Bnd_Spline_West,
			Number_of_Cells_Idir,
			Number_of_Cells_Jdir,
			Number_of_Ghost_Cells,
			GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
			Stretch_I,
			Beta_I, 
			Tau_I,
			Stretch_J,
			Beta_J,
			Tau_J,
			Orthogonal_North,
			Orthogonal_South,
			Orthogonal_East,
			Orthogonal_West);
      // Deallocate the memory for the boundary splines.
      Bnd_Spline_North.deallocate();
      Bnd_Spline_South.deallocate();
      Bnd_Spline_East.deallocate();
      Bnd_Spline_West.deallocate();

    }

  }
  //if (iBlk == Number_of_Blocks_Idir-1) 
  //  Smooth_Quad_Block(Grid_ptr[iBlk][0], min(250, 2*max(n_cells_i,n_cells_j)));

  // Return the grid.
  return Grid_ptr;
  
} 

/********************************************************
 * Routine: Grid_Circular_Cylinder                      *
 *                                                      *
 * Generates a double-block O-type grid for predicting  *
 * flow past a circular cylinder.                       *
 *                                                      *
 * Usage: Grid_ptr = Grid_Circular_Cylinder(Grid_ptr,   *
 *                                          nblk_i,     *
 *                                          nblk_j,     *
 *                                          THREE,      *
 *   		                            100,        *
 *  		                            50,         *
 *                                          2);         *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Grid_Circular_Cylinder(Grid2D_Quad_Block **Grid_ptr,
                                           int &Number_of_Blocks_Idir,
                                           int &Number_of_Blocks_Jdir,
                                           const double &Radius,
 		                           const int Number_of_Cells_Idir,
		                           const int Number_of_Cells_Jdir,
					   const int Number_of_Ghost_Cells) {

    int iBlk, n_cells_i, n_cells_j, Stretch_I, Stretch_J,
        Orthogonal_North, Orthogonal_South,
        Orthogonal_East, Orthogonal_West;
    double Beta_I, Tau_I, Beta_J, Tau_J;
    Vector2D x1, x2;
    Spline2D Bnd_Spline_North, Bnd_Spline_South,
             Bnd_Spline_East, Bnd_Spline_West;

    /* Allocate memory for grid blocks.  There are two grid
       blocks for this mesh. */

    Number_of_Blocks_Idir = 2;
    Number_of_Blocks_Jdir = 1;
    Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
                                         Number_of_Blocks_Idir, 
                                         Number_of_Blocks_Jdir);

    /* Create the mesh for each block representing
       the complete grid. */

    for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {

        /* Create the splines defining the north, south,
           east, and west boundaries of the grid. */

        if (iBlk == 0) {
           x1 = Vector2D(ZERO , ZERO);
           Create_Spline_Circular_Arc(Bnd_Spline_North,
	      		              x1,
			              32.00*Radius,
                                      360.00,
			              180.00,
  	                              361);
           Create_Spline_Circular_Arc(Bnd_Spline_South,
			              x1,
			              Radius,
                                      360.00,
			              180.00,
  	                              361);
           x1 = Vector2D(Radius, ZERO);
           x2 = Vector2D(32.00*Radius, ZERO);
           Create_Spline_Line(Bnd_Spline_West, x1, x2, 2);
           x1 = Vector2D(-Radius, ZERO);
           x2 = Vector2D(-32.00*Radius, ZERO);
           Create_Spline_Line(Bnd_Spline_East, x1, x2, 2);
        } else {
           x1 = Vector2D(ZERO , ZERO);
           Create_Spline_Circular_Arc(Bnd_Spline_North,
	      		              x1,
			              32.00*Radius,
                                      180.00,
			              ZERO,
  	                              361);
           Create_Spline_Circular_Arc(Bnd_Spline_South,
			              x1,
			              Radius,
                                      180.00,
			              ZERO,
  	                              361);
           x1 = Vector2D(-Radius, ZERO);
           x2 = Vector2D(-32.00*Radius, ZERO);
           Create_Spline_Line(Bnd_Spline_West, x1, x2, 2);
           x1 = Vector2D(Radius, ZERO);
           x2 = Vector2D(32.00*Radius, ZERO);
           Create_Spline_Line(Bnd_Spline_East, x1, x2, 2);
        } /* endif */

        /* Set the boundary condition types for each of the
           boundary splines. */

        if (iBlk == 0) {
           Bnd_Spline_North.setBCtype(BC_FIXED);
           Bnd_Spline_South.setBCtype(BC_REFLECTION);
           Bnd_Spline_East.setBCtype(BC_NONE);
           Bnd_Spline_West.setBCtype(BC_NONE);
        } else {
           Bnd_Spline_North.setBCtype(BC_FIXED);
           Bnd_Spline_South.setBCtype(BC_REFLECTION);
           Bnd_Spline_East.setBCtype(BC_NONE);
           Bnd_Spline_West.setBCtype(BC_NONE);
        } /* endif */

        /* Determine the number of cells for this block. */

        n_cells_i = Number_of_Cells_Idir/2;
	n_cells_j = Number_of_Cells_Jdir;

        /* Assign values to the stretching function parameters
           and boundary grid line orthogonality parameters. */

        if (iBlk == 0) {
           Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
           Beta_I = 1.025;
           Tau_I = ZERO;
           Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
           Beta_J = 1.001;
           Tau_J = ZERO;
           Orthogonal_North = 0;
           Orthogonal_South = 0;
           Orthogonal_East = 0;
           Orthogonal_West = 0;
        } else {
           Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
           Beta_I = 1.025; 
           Tau_I = ZERO;
           Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
           Beta_J = 1.001;
           Tau_J = ZERO;
           Orthogonal_North = 0;
           Orthogonal_South = 0;
           Orthogonal_East = 0;
           Orthogonal_West = 0;
        } /* endif */

        /* Create the 2D quadrilateral grid block. */

        Create_Quad_Block(Grid_ptr[iBlk][0],
                          Bnd_Spline_North,
                          Bnd_Spline_South,
                          Bnd_Spline_East,
                          Bnd_Spline_West,
                          n_cells_i,
  	                  n_cells_j,
			  Number_of_Ghost_Cells,
                          GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
                          Stretch_I,
                          Beta_I, 
                          Tau_I,
                          Stretch_J,
                          Beta_J,
                          Tau_J,
                          Orthogonal_North,
		          Orthogonal_South,
     		          Orthogonal_East,
                          Orthogonal_West);

        /* Deallocate the memory for the boundary splines. */

        Bnd_Spline_North.deallocate();
        Bnd_Spline_South.deallocate();
        Bnd_Spline_East.deallocate();
        Bnd_Spline_West.deallocate();

    } /* endfor */

    /* Return the grid. */

    return(Grid_ptr);

}

/********************************************************
 * Routine: Grid_Ellipse                                *
 *                                                      *
 * Generates a double-block O-type grid for predicting  *
 * flow past an ellipse.                                *
 *                                                      *
 * Usage: Grid_ptr = Grid_Ellipse(Grid_ptr,             *
 *                                nblk_i,               *
 *                                nblk_j,               *
 *                                FOUR,                 *
 *                                ONE,                  *
 *   		                  100,                  *
 *  		                  50);                  *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Grid_Ellipse(Grid2D_Quad_Block **Grid_ptr,
                                 int &Number_of_Blocks_Idir,
                                 int &Number_of_Blocks_Jdir,
                                 const double &A,
                                 const double &B,
 		                 const int Number_of_Cells_Idir,
		                 const int Number_of_Cells_Jdir,
				 const int Number_of_Ghost_Cells) {

    int iBlk, n_cells_i, n_cells_j, Stretch_I, Stretch_J,
        Orthogonal_North, Orthogonal_South,
        Orthogonal_East, Orthogonal_West;
    double Beta_I, Tau_I, Beta_J, Tau_J;
    Vector2D x1, x2;
    Spline2D Bnd_Spline_North, Bnd_Spline_South,
             Bnd_Spline_East, Bnd_Spline_West;

    /* Allocate memory for grid blocks.  There are two grid
       blocks for this mesh. */

    Number_of_Blocks_Idir = 2;
    Number_of_Blocks_Jdir = 1;
    Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
                                         Number_of_Blocks_Idir, 
                                         Number_of_Blocks_Jdir);

    /* Create the mesh for each block representing
       the complete grid. */

    for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {

        /* Create the splines defining the north, south,
           east, and west boundaries of the grid. */

        if (iBlk == 0) {
           x1 = Vector2D(ZERO , ZERO);
           Create_Spline_Circular_Arc(Bnd_Spline_North,
	      		              x1,
			              32.00*A,
                                      360.00,
			              180.00,
  	                              361);
           Create_Spline_Ellipsoidal_Arc(Bnd_Spline_South,
			                 x1,
			                 A,
                                         B,
                                         360.00,
			                 180.00,
  	                                 361);
           x1 = Vector2D(A, ZERO);
           x2 = Vector2D(32.00*A, ZERO);
           Create_Spline_Line(Bnd_Spline_West, x1, x2, 2);
           x1 = Vector2D(-A, ZERO);
           x2 = Vector2D(-32.00*A, ZERO);
           Create_Spline_Line(Bnd_Spline_East, x1, x2, 2);
        } else {
           x1 = Vector2D(ZERO , ZERO);
           Create_Spline_Circular_Arc(Bnd_Spline_North,
	      		              x1,
			              32.00*A,
                                      180.00,
			              ZERO,
  	                              361);
           Create_Spline_Ellipsoidal_Arc(Bnd_Spline_South,
			                 x1,
			                 A,
                                         B,
                                         180.00,
			                 ZERO,
  	                                 361);
           x1 = Vector2D(-A, ZERO);
           x2 = Vector2D(-32.00*A, ZERO);
           Create_Spline_Line(Bnd_Spline_West, x1, x2, 2);
           x1 = Vector2D(A, ZERO);
           x2 = Vector2D(32.00*A, ZERO);
           Create_Spline_Line(Bnd_Spline_East, x1, x2, 2);
        } /* endif */

        /* Set the boundary condition types for each of the
           boundary splines. */

        if (iBlk == 0) {
           Bnd_Spline_North.setBCtype(BC_FIXED);
           Bnd_Spline_South.setBCtype(BC_REFLECTION);
           Bnd_Spline_East.setBCtype(BC_NONE);
           Bnd_Spline_West.setBCtype(BC_NONE);
        } else {
           Bnd_Spline_North.setBCtype(BC_FIXED);
           Bnd_Spline_South.setBCtype(BC_REFLECTION);
           Bnd_Spline_East.setBCtype(BC_NONE);
           Bnd_Spline_West.setBCtype(BC_NONE);
        } /* endif */

        /* Determine the number of cells for this block. */

        n_cells_i = Number_of_Cells_Idir/2;
	n_cells_j = Number_of_Cells_Jdir;

        /* Assign values to the stretching function parameters
           and boundary grid line orthogonality parameters. */

        if (iBlk == 0) {
           Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
           Beta_I = 1.025; 
           Tau_I = ZERO;
           Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
           Beta_J = 1.001;
           Tau_J = ZERO;
           Orthogonal_North = 1;
           Orthogonal_South = 1;
           Orthogonal_East = 1;
           Orthogonal_West = 1;
        } else {
           Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
           Beta_I = 1.025; 
           Tau_I = ZERO;
           Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
           Beta_J = 1.001;
           Tau_J = ZERO;
           Orthogonal_North = 1;
           Orthogonal_South = 1;
           Orthogonal_East = 1;
           Orthogonal_West = 1;
        } /* endif */

        /* Create the 2D quadrilateral grid block. */

        Create_Quad_Block(Grid_ptr[iBlk][0],
                          Bnd_Spline_North,
                          Bnd_Spline_South,
                          Bnd_Spline_East,
                          Bnd_Spline_West,
                          n_cells_i,
  	                  n_cells_j,
			  Number_of_Ghost_Cells,
                          GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
                          Stretch_I,
                          Beta_I, 
                          Tau_I,
                          Stretch_J,
                          Beta_J,
                          Tau_J,
                          Orthogonal_North,
		          Orthogonal_South,
     		          Orthogonal_East,
                          Orthogonal_West);

        /* Smooth the 2D quadrilateral grid block. */

        Smooth_Quad_Block(Grid_ptr[iBlk][0],
                          min(250, 2*max(n_cells_i,n_cells_j)));

        /* Deallocate the memory for the boundary splines. */

        Bnd_Spline_North.deallocate();
        Bnd_Spline_South.deallocate();
        Bnd_Spline_East.deallocate();
        Bnd_Spline_West.deallocate();

    } /* endfor */

    /* Return the grid. */

    return(Grid_ptr);

}

/********************************************************
 * Routine: Grid_NACA_Aerofoil                          *
 *                                                      *
 * Generates a C-type grid consisting of four           *
 * quadrilateral grid blocks for predicting flow past   *
 * NACA 4-digit and 5-digit aerofoils.                  *
 *                                                      *
 * Usage:  Grid_ptr = Grid_NACA_Aerofoil(Grid_ptr,      *
 *                                       nblk_i,        *
 *                                       nblk_j,        *
 *                                       "4412",        *
 *                                       ONE,           *
 *		                         120,           *
 *		                         50,            *
 *                                       2);            *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Grid_NACA_Aerofoil(Grid2D_Quad_Block **Grid_ptr,
                                       int &Number_of_Blocks_Idir,
                                       int &Number_of_Blocks_Jdir,
                                       char *NACA_Aerofoil_Type_ptr,
                                       const double &Chord_Length,
 		                       const int Number_of_Cells_Idir,
		                       const int Number_of_Cells_Jdir,
				       const int Number_of_Ghost_Cells) {

    int iBlk, n_cells_i, n_cells_j,
        Stretch_I, Stretch_J,
        Orthogonal_North, Orthogonal_South,
        Orthogonal_East, Orthogonal_West;
    double Beta_I, Tau_I, Beta_J, Tau_J;
    Vector2D x1, x2;
    Spline2D Bnd_Spline_North, Bnd_Spline_South,
             Bnd_Spline_East, Bnd_Spline_West,
             s_tmp1, s_tmp2;

    /* Allocate memory for grid blocks.  There are four grid
       blocks for this mesh. */

    Number_of_Blocks_Idir = 4;
    Number_of_Blocks_Jdir = 1;
    Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
                                         Number_of_Blocks_Idir, 
                                         Number_of_Blocks_Jdir);

    /* Create the mesh for each block representing
       the complete grid. */

    for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {

        /* Create the splines defining the north, south,
           east, and west boundaries of the grid. */

        if (iBlk == 0) {
           x1 = Vector2D(32.00*Chord_Length, -32.00*Chord_Length);
           x2 = Vector2D(Chord_Length, -32.00*Chord_Length);
           Create_Spline_Line(Bnd_Spline_North, x1, x2, 2);
           x1 = Vector2D(32.00*Chord_Length, ZERO);
           x2 = Vector2D(Chord_Length, ZERO);
           Create_Spline_Line(Bnd_Spline_South, x1, x2, 2);
           x1 = Vector2D(32.00*Chord_Length, ZERO);
           x2 = Vector2D(32.00*Chord_Length, -32.00*Chord_Length);
           Create_Spline_Line(Bnd_Spline_West, x1, x2, 2);
           x1 = Vector2D(Chord_Length, ZERO);
           x2 = Vector2D(Chord_Length, -32.00*Chord_Length);
           Create_Spline_Line(Bnd_Spline_East, x1, x2, 2);
        } else if (iBlk == 1) {
           x1 = Vector2D(Chord_Length, -32.00*Chord_Length);
           x2 = Vector2D(ZERO, -32.00*Chord_Length);
           Create_Spline_Line(s_tmp1, x1, x2, 2);
           x1 = Vector2D(ZERO , ZERO);
           Create_Spline_Circular_Arc(s_tmp2,
	      		              x1,
			              32.00*Chord_Length,
                                      270.00,
			              180.00,
  	                              181);
           Bnd_Spline_North = Concatenate_Splines(s_tmp1, s_tmp2);
           s_tmp1.deallocate();
           s_tmp2.deallocate();
           Create_Spline_NACA_Aerofoil(Bnd_Spline_South,
                                       NACA_Aerofoil_Type_ptr,
                                       Chord_Length,
                                       -1,
	  		               501);
           x1 = Vector2D(Chord_Length, ZERO);
           x2 = Vector2D(Chord_Length, -32.00*Chord_Length);
           Create_Spline_Line(Bnd_Spline_West, x1, x2, 2);
           x1 = Vector2D(ZERO, ZERO);
           x2 = Vector2D(-32.00*Chord_Length, ZERO);
           Create_Spline_Line(Bnd_Spline_East, x1, x2, 2);
        } else if (iBlk == 2) {
           x1 = Vector2D(ZERO , ZERO);
           Create_Spline_Circular_Arc(s_tmp1,
	      		              x1,
			              32.00*Chord_Length,
                                      180.00,
			              90.00,
  	                              181);
           x1 = Vector2D(ZERO, 32.00*Chord_Length);
           x2 = Vector2D(Chord_Length, 32.00*Chord_Length);
           Create_Spline_Line(s_tmp2, x1, x2, 2);
           Bnd_Spline_North = Concatenate_Splines(s_tmp1, s_tmp2);
           s_tmp1.deallocate();
           s_tmp2.deallocate();
           Create_Spline_NACA_Aerofoil(Bnd_Spline_South,
                                       NACA_Aerofoil_Type_ptr,
                                       Chord_Length,
                                       1,
	  		               501);
           x1 = Vector2D(ZERO, ZERO);
           x2 = Vector2D(-32.00*Chord_Length, ZERO);
           Create_Spline_Line(Bnd_Spline_West, x1, x2, 2);
           x1 = Vector2D(Chord_Length, ZERO);
           x2 = Vector2D(Chord_Length, 32.00*Chord_Length);
           Create_Spline_Line(Bnd_Spline_East, x1, x2, 2);
        } else {
           x1 = Vector2D(Chord_Length, 32.00*Chord_Length);
           x2 = Vector2D(32.00*Chord_Length, 32.00*Chord_Length);
           Create_Spline_Line(Bnd_Spline_North, x1, x2, 2);
           x1 = Vector2D(Chord_Length, ZERO);
           x2 = Vector2D(32.00*Chord_Length, ZERO);
           Create_Spline_Line(Bnd_Spline_South, x1, x2, 2);
           x1 = Vector2D(Chord_Length, ZERO);
           x2 = Vector2D(Chord_Length, 32.00*Chord_Length);
           Create_Spline_Line(Bnd_Spline_West, x1, x2, 2);
           x1 = Vector2D(32.00*Chord_Length, ZERO);
           x2 = Vector2D(32.00*Chord_Length, 32.00*Chord_Length);
           Create_Spline_Line(Bnd_Spline_East, x1, x2, 2);
        } /* endif */

        /* Set the boundary condition types for each of the
           boundary splines. */

        if (iBlk == 0) {
           Bnd_Spline_North.setBCtype(BC_FIXED);
           Bnd_Spline_South.setBCtype(BC_NONE);
           Bnd_Spline_East.setBCtype(BC_NONE);
           Bnd_Spline_West.setBCtype(BC_FIXED);
        } else if (iBlk == 1) {
           Bnd_Spline_North.setBCtype(BC_FIXED);
           Bnd_Spline_South.setBCtype(BC_REFLECTION);
           Bnd_Spline_East.setBCtype(BC_NONE);
           Bnd_Spline_West.setBCtype(BC_NONE);
        } else if (iBlk == 2) {
           Bnd_Spline_North.setBCtype(BC_FIXED);
           Bnd_Spline_South.setBCtype(BC_REFLECTION);
           Bnd_Spline_East.setBCtype(BC_NONE);
           Bnd_Spline_West.setBCtype(BC_NONE);
        } else {
           Bnd_Spline_North.setBCtype(BC_FIXED);
           Bnd_Spline_South.setBCtype(BC_NONE);
           Bnd_Spline_East.setBCtype(BC_FIXED);
           Bnd_Spline_West.setBCtype(BC_NONE);
        } /* endif */

        /* Determine the number of cells for this block. */

        if (iBlk == 1 || iBlk == 2) {
           n_cells_i = 2*(((Number_of_Cells_Idir/3)+1)/2);
	   n_cells_j = Number_of_Cells_Jdir;
        } else {
           n_cells_i = 2*(((Number_of_Cells_Idir/6)+1)/2);
	   n_cells_j = Number_of_Cells_Jdir;
        } /* endif */
        n_cells_i = Number_of_Cells_Idir;
        n_cells_j = Number_of_Cells_Jdir;

        /* Assign values to the stretching function parameters
           and boundary grid line orthogonality parameters. */

        if (iBlk == 0) {
           Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
           Beta_I = 1.000075; 
           Tau_I = ZERO;
           Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
           Beta_J = 1.00015;
           Tau_J = ZERO;
           Orthogonal_North = 0;
           Orthogonal_South = 0;
           Orthogonal_East = 0;
           Orthogonal_West = 0;
        } else if (iBlk == 1) {
           Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
           Beta_I = 1.01; 
           Tau_I = ZERO;
           Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
           Beta_J = 1.00015;
           Tau_J = ZERO;
           Orthogonal_North = 1;
           Orthogonal_South = 1;
           Orthogonal_East = 1;
           Orthogonal_West = 1;
        } else if (iBlk == 2) {
           Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
           Beta_I = 1.01; 
           Tau_I = ZERO;
           Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
           Beta_J = 1.00015;
           Tau_J = ZERO;
           Orthogonal_North = 1;
           Orthogonal_South = 1;
           Orthogonal_East = 1;
           Orthogonal_West = 1;
        } else {
           Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
           Beta_I = 1.000075; 
           Tau_I = ZERO;
           Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
           Beta_J = 1.00015;
           Tau_J = ZERO;
           Orthogonal_North = 0;
           Orthogonal_South = 0;
           Orthogonal_East = 0;
           Orthogonal_West = 0;
        } /* endif */

        /* Create the 2D quadrilateral grid block. */

        Create_Quad_Block(Grid_ptr[iBlk][0],
                          Bnd_Spline_North,
                          Bnd_Spline_South,
                          Bnd_Spline_East,
                          Bnd_Spline_West,
                          n_cells_i,
  	                  n_cells_j,
			  Number_of_Ghost_Cells,
                          GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
                          Stretch_I,
                          Beta_I, 
                          Tau_I,
                          Stretch_J,
                          Beta_J,
                          Tau_J,
                          Orthogonal_North,
		          Orthogonal_South,
     		          Orthogonal_East,
                          Orthogonal_West);

        /* Smooth the 2D quadrilateral grid block. */

        if (iBlk == 1 || iBlk == 2) 
           Smooth_Quad_Block(Grid_ptr[iBlk][0], min(250, 2*max(n_cells_i,n_cells_j)));

        /* Deallocate the memory for the boundary splines. */

        Bnd_Spline_North.deallocate();
        Bnd_Spline_South.deallocate();
        Bnd_Spline_East.deallocate();
        Bnd_Spline_West.deallocate();

        /* Force the mesh to be symmetric. */

//          if (iBlk == 2) {
//             Copy_Quad_Block(Grid_ptr[2][0],
//                             Grid_ptr[1][0]);
//             Copy_Spline(Grid_ptr[2][0].BndNorthSpline,
//                         Grid_ptr[1][0].BndNorthSpline);
//             Copy_Spline(Grid_ptr[2][0].BndSouthSpline,
//                         Grid_ptr[1][0].BndSouthSpline);
//             Copy_Spline(Grid_ptr[2][0].BndEastSpline,
//                         Grid_ptr[1][0].BndEastSpline);
//             Copy_Spline(Grid_ptr[2][0].BndWestSpline,
//                         Grid_ptr[1][0].BndWestSpline);
//             Reflect_Quad_Block(Grid_ptr[2][0]);
//          } else if (iBlk == 3) {
//             Copy_Quad_Block(Grid_ptr[3][0],
//                             Grid_ptr[0][0]);
//             Copy_Spline(Grid_ptr[3][0].BndNorthSpline,
//                         Grid_ptr[0][0].BndNorthSpline);
//             Copy_Spline(Grid_ptr[3][0].BndSouthSpline,
//                         Grid_ptr[0][0].BndSouthSpline);
//             Copy_Spline(Grid_ptr[3][0].BndEastSpline,
//                         Grid_ptr[0][0].BndEastSpline);
//             Copy_Spline(Grid_ptr[3][0].BndWestSpline,
//                         Grid_ptr[0][0].BndWestSpline);
//             Reflect_Quad_Block(Grid_ptr[3][0]);
//          } /* endif */

    } /* endfor */

    /* Return the grid. */

    return(Grid_ptr);

}

/********************************************************
 * Routine: Grid_Free_Jet                               *
 *                                                      *
 * Generates a multi-block grid for predicting free-jet *
 * flows associated with expansion through an orifice.  *
 *                                                      *
 * Usage: Grid_ptr = Grid_Free_Jet(Grid_ptr,            *
 *                                 nblk_i,              *
 *                                 nblk_j,              *
 *                                 THREE,               *
 *   		                   100,                 *
 *  		                   50,                  *
 *                                 2);                  *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Grid_Free_Jet(Grid2D_Quad_Block **Grid_ptr,
                                  int &Number_of_Blocks_Idir,
                                  int &Number_of_Blocks_Jdir,
                                  const double &Radius,
 		                  const int Number_of_Cells_Idir,
		                  const int Number_of_Cells_Jdir,
				  const int Number_of_Ghost_Cells) {

    int iBlk, jBlk, n_cells_i, n_cells_j, Stretch_I, Stretch_J,
        Orthogonal_North, Orthogonal_South,
        Orthogonal_East, Orthogonal_West;
    double Beta_I, Tau_I, Beta_J, Tau_J;
    Vector2D x1, x2;
    Spline2D Bnd_Spline_North, Bnd_Spline_South,
             Bnd_Spline_East, Bnd_Spline_West;

    double x_orifice_up, x_orifice_down, x_up, x_down, r_max;

    /* Allocate memory for grid blocks.  There are two grid
       blocks for this mesh. */

    Number_of_Blocks_Idir = 3;
    Number_of_Blocks_Jdir = 2;
    Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
                                         Number_of_Blocks_Idir, 
                                         Number_of_Blocks_Jdir);

    /* Create the mesh for each block representing
       the complete grid. */

    for ( jBlk = 0; jBlk <= Number_of_Blocks_Jdir-1; ++jBlk ) {
       for ( iBlk = 0; iBlk <= Number_of_Blocks_Idir-1; ++iBlk ) {

           x_orifice_up = -EIGHT*Radius;
           x_orifice_down = ZERO;
           x_up = -24*Radius;
           x_down = 60*Radius;
           r_max = 40*Radius;

           if (iBlk == 0 && jBlk == 0) {
              x1 = Vector2D(x_up, Radius);
              x2 = Vector2D(x_orifice_up, Radius);
              Create_Spline_Line(Bnd_Spline_North, x1, x2, 2);
              x1 = Vector2D(x_up, ZERO);
              x2 = Vector2D(x_orifice_up, ZERO);
              Create_Spline_Line(Bnd_Spline_South, x1, x2, 2);
              x1 = Vector2D(x_up, ZERO);
              x2 = Vector2D(x_up, Radius);
              Create_Spline_Line(Bnd_Spline_West, x1, x2, 2);
              x1 = Vector2D(x_orifice_up, ZERO);
              x2 = Vector2D(x_orifice_up, Radius);
              Create_Spline_Line(Bnd_Spline_East, x1, x2, 2);

              Bnd_Spline_North.setBCtype(BC_NONE);
              Bnd_Spline_South.setBCtype(BC_REFLECTION);
              Bnd_Spline_East.setBCtype(BC_NONE);
              Bnd_Spline_West.setBCtype(BC_FIXED);

              n_cells_i = Number_of_Cells_Idir/3;
      	      n_cells_j = Number_of_Cells_Jdir/2;

              Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
              Beta_I = 1.05; 
              Tau_I = ZERO;
              Stretch_J = STRETCHING_FCN_MINMAX_CLUSTERING;
              Beta_J = 1.5;
              Tau_J = ZERO;
              Orthogonal_North = 0;
              Orthogonal_South = 0;
              Orthogonal_East = 0;
              Orthogonal_West = 0;

              Create_Quad_Block(Grid_ptr[iBlk][jBlk],
                                Bnd_Spline_North,
                                Bnd_Spline_South,
                                Bnd_Spline_East,
                                Bnd_Spline_West,
                                n_cells_i,
  	                        n_cells_j,
				Number_of_Ghost_Cells,
                                GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
                                Stretch_I,
                                Beta_I, 
                                Tau_I,
                                Stretch_J,
                                Beta_J,
                                Tau_J,
                                Orthogonal_North,
		                Orthogonal_South,
     		                Orthogonal_East,
                                Orthogonal_West);

              Bnd_Spline_North.deallocate();
              Bnd_Spline_South.deallocate();
              Bnd_Spline_East.deallocate();
              Bnd_Spline_West.deallocate();

           } else if (iBlk == 1 && jBlk == 0) {
              x1 = Vector2D(x_orifice_up, Radius);
              x2 = Vector2D(x_orifice_down, Radius);
              Create_Spline_Line(Bnd_Spline_North, x1, x2, 2);
              x1 = Vector2D(x_orifice_up, ZERO);
              x2 = Vector2D(x_orifice_down, ZERO);
              Create_Spline_Line(Bnd_Spline_South, x1, x2, 2);
              x1 = Vector2D(x_orifice_up, ZERO);
              x2 = Vector2D(x_orifice_up, Radius);
              Create_Spline_Line(Bnd_Spline_West, x1, x2, 2);
              x1 = Vector2D(x_orifice_down, ZERO);
              x2 = Vector2D(x_orifice_down, Radius);
              Create_Spline_Line(Bnd_Spline_East, x1, x2, 2);

              Bnd_Spline_North.setBCtype(BC_REFLECTION);
              Bnd_Spline_South.setBCtype(BC_REFLECTION);
              Bnd_Spline_East.setBCtype(BC_NONE);
              Bnd_Spline_West.setBCtype(BC_NONE);

              n_cells_i = Number_of_Cells_Idir/3;
      	      n_cells_j = Number_of_Cells_Jdir/2;

              Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
              Beta_I = 1.1; 
              Tau_I = ZERO;
              Stretch_J = STRETCHING_FCN_MINMAX_CLUSTERING;
              Beta_J = 1.5;
              Tau_J = ZERO;
              Orthogonal_North = 0;
              Orthogonal_South = 0;
              Orthogonal_East = 0;
              Orthogonal_West = 0;

              Create_Quad_Block(Grid_ptr[iBlk][jBlk],
                                Bnd_Spline_North,
                                Bnd_Spline_South,
                                Bnd_Spline_East,
                                Bnd_Spline_West,
                                n_cells_i,
  	                        n_cells_j,
				Number_of_Ghost_Cells,
                                GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
                                Stretch_I,
                                Beta_I, 
                                Tau_I,
                                Stretch_J,
                                Beta_J,
                                Tau_J,
                                Orthogonal_North,
		                Orthogonal_South,
     		                Orthogonal_East,
                                Orthogonal_West);

              Bnd_Spline_North.deallocate();
              Bnd_Spline_South.deallocate();
              Bnd_Spline_East.deallocate();
              Bnd_Spline_West.deallocate();

           } else if (iBlk == 2 && jBlk == 0) {
              x1 = Vector2D(x_orifice_down, Radius);
              x2 = Vector2D(x_down, Radius);
              Create_Spline_Line(Bnd_Spline_North, x1, x2, 2);
              x1 = Vector2D(x_orifice_down, ZERO);
              x2 = Vector2D(x_down, ZERO);
              Create_Spline_Line(Bnd_Spline_South, x1, x2, 2);
              x1 = Vector2D(x_orifice_down, ZERO);
              x2 = Vector2D(x_orifice_down, Radius);
              Create_Spline_Line(Bnd_Spline_West, x1, x2, 2);
              x1 = Vector2D(x_down, ZERO);
              x2 = Vector2D(x_down, Radius);
              Create_Spline_Line(Bnd_Spline_East, x1, x2, 2);

              Bnd_Spline_North.setBCtype(BC_NONE);
              Bnd_Spline_South.setBCtype(BC_REFLECTION);
              Bnd_Spline_East.setBCtype(BC_FIXED);
              Bnd_Spline_West.setBCtype(BC_NONE);

              n_cells_i = Number_of_Cells_Idir/3;
      	      n_cells_j = Number_of_Cells_Jdir/2;

              Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
              Beta_I = 1.01; 
              Tau_I = ZERO;
              Stretch_J = STRETCHING_FCN_MINMAX_CLUSTERING;
              Beta_J = 1.5;
              Tau_J = ZERO;
              Orthogonal_North = 0;
              Orthogonal_South = 0;
              Orthogonal_East = 0;
              Orthogonal_West = 0;

              Create_Quad_Block(Grid_ptr[iBlk][jBlk],
                                Bnd_Spline_North,
                                Bnd_Spline_South,
                                Bnd_Spline_East,
                                Bnd_Spline_West,
                                n_cells_i,
  	                        n_cells_j,
				Number_of_Ghost_Cells,
                                GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
                                Stretch_I,
                                Beta_I, 
                                Tau_I,
                                Stretch_J,
                                Beta_J,
                                Tau_J,
                                Orthogonal_North,
		                Orthogonal_South,
     		                Orthogonal_East,
                                Orthogonal_West);

              Bnd_Spline_North.deallocate();
              Bnd_Spline_South.deallocate();
              Bnd_Spline_East.deallocate();
              Bnd_Spline_West.deallocate();

           } else if (iBlk == 0 && jBlk == 1) {
              x1 = Vector2D(x_up, r_max);
              x2 = Vector2D(x_orifice_up, r_max);
              Create_Spline_Line(Bnd_Spline_North, x1, x2, 2);
              x1 = Vector2D(x_up, Radius);
              x2 = Vector2D(x_orifice_up, Radius);
              Create_Spline_Line(Bnd_Spline_South, x1, x2, 2);
              x1 = Vector2D(x_up, Radius);
              x2 = Vector2D(x_up, r_max);
              Create_Spline_Line(Bnd_Spline_West, x1, x2, 2);
              x1 = Vector2D(x_orifice_up, Radius);
              x2 = Vector2D(x_orifice_up, r_max);
              Create_Spline_Line(Bnd_Spline_East, x1, x2, 2);

              Bnd_Spline_North.setBCtype(BC_FIXED);
              Bnd_Spline_South.setBCtype(BC_NONE);
              Bnd_Spline_East.setBCtype(BC_REFLECTION);
              Bnd_Spline_West.setBCtype(BC_FIXED);

              n_cells_i = Number_of_Cells_Idir/3;
      	      n_cells_j = Number_of_Cells_Jdir/2;

              Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
              Beta_I = 1.05; 
              Tau_I = ZERO;
              Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
              Beta_J = 1.004;
              Tau_J = ZERO;
              Orthogonal_North = 0;
              Orthogonal_South = 0;
              Orthogonal_East = 0;
              Orthogonal_West = 0;

              Create_Quad_Block(Grid_ptr[iBlk][jBlk],
                                Bnd_Spline_North,
                                Bnd_Spline_South,
                                Bnd_Spline_East,
                                Bnd_Spline_West,
                                n_cells_i,
  	                        n_cells_j,
				Number_of_Ghost_Cells,
                                GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
                                Stretch_I,
                                Beta_I, 
                                Tau_I,
                                Stretch_J,
                                Beta_J,
                                Tau_J,
                                Orthogonal_North,
		                Orthogonal_South,
     		                Orthogonal_East,
                                Orthogonal_West);

              Bnd_Spline_North.deallocate();
              Bnd_Spline_South.deallocate();
              Bnd_Spline_East.deallocate();
              Bnd_Spline_West.deallocate();

           } else if (iBlk == 1 && jBlk == 1) {

           } else if (iBlk == 2 && jBlk == 1) {
              x1 = Vector2D(x_orifice_down, r_max);
              x2 = Vector2D(x_down, r_max);
              Create_Spline_Line(Bnd_Spline_North, x1, x2, 2);
              x1 = Vector2D(x_orifice_down, Radius);
              x2 = Vector2D(x_down, Radius);
              Create_Spline_Line(Bnd_Spline_South, x1, x2, 2);
              x1 = Vector2D(x_orifice_down, Radius);
              x2 = Vector2D(x_orifice_down, r_max);
              Create_Spline_Line(Bnd_Spline_West, x1, x2, 2);
              x1 = Vector2D(x_down, Radius);
              x2 = Vector2D(x_down, r_max);
              Create_Spline_Line(Bnd_Spline_East, x1, x2, 2);

              Bnd_Spline_North.setBCtype(BC_FIXED);
              Bnd_Spline_South.setBCtype(BC_NONE);
              Bnd_Spline_East.setBCtype(BC_FIXED);
              Bnd_Spline_West.setBCtype(BC_REFLECTION);

              n_cells_i = Number_of_Cells_Idir/3;
      	      n_cells_j = Number_of_Cells_Jdir/2;

              Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
              Beta_I = 1.01; 
              Tau_I = ZERO;
              Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
              Beta_J = 1.004;
              Tau_J = ZERO;
              Orthogonal_North = 0;
              Orthogonal_South = 0;
              Orthogonal_East = 0;
              Orthogonal_West = 0;

              Create_Quad_Block(Grid_ptr[iBlk][jBlk],
                                Bnd_Spline_North,
                                Bnd_Spline_South,
                                Bnd_Spline_East,
                                Bnd_Spline_West,
                                n_cells_i,
  	                        n_cells_j,
				Number_of_Ghost_Cells,
                                GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
                                Stretch_I,
                                Beta_I, 
                                Tau_I,
                                Stretch_J,
                                Beta_J,
                                Tau_J,
                                Orthogonal_North,
		                Orthogonal_South,
     		                Orthogonal_East,
                                Orthogonal_West);

              Bnd_Spline_North.deallocate();
              Bnd_Spline_South.deallocate();
              Bnd_Spline_East.deallocate();
              Bnd_Spline_West.deallocate();
           } /* endif */

       } /* endfor */
    } /* endfor */

    /* Return the grid. */

    return(Grid_ptr);

}

/**********************************************************************
 * Routine: Grid_Wedge                                                *
 *                                                                    *
 * Generates a quadilateral mesh with clustering consisting of two    *
 * grid blocks for predicting a wedge flow for looking at shock wave  *
 * interactions at a wedge.                                           *
 *                                                                    *
 * Usage: Grid_ptr = Grid_Wedge(Grid_ptr,                             *
 *                              nblk_i,                               *
 *                              nblk_j,                               *
 *                              45.0,                                 *
 *                              0.50,                                 *
 *         		        100,                                  *
 *         		        100,                                  *
 *                              2);                                   *
 *                                                                    *
 **********************************************************************/
Grid2D_Quad_Block** Grid_Wedge(Grid2D_Quad_Block **Grid_ptr,
                               int &Number_of_Blocks_Idir,
                               int &Number_of_Blocks_Jdir,
			       const double &Wedge_Angle,
			       const double &Wedge_Length,
			       const int Number_of_Cells_Idir,
			       const int Number_of_Cells_Jdir,
			       const int Number_of_Ghost_Cells) {

  int Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D Bnd_Spline_North, Bnd_Spline_South,
           Bnd_Spline_East, Bnd_Spline_West;

  // Allocate memory for grid blocks.  There are two grid blocks for 
  // this mesh.
  Number_of_Blocks_Idir = 2;
  Number_of_Blocks_Jdir = 1;
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
				       Number_of_Blocks_Idir, 
				       Number_of_Blocks_Jdir);
  
  // Create the mesh for each block representing the complete grid.
  for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

    if (iBlk == 0) {
      xc_NW = Vector2D(-Wedge_Length,Wedge_Length);
      xc_NE = Vector2D(ZERO         ,Wedge_Length);
      xc_SE = Vector2D(ZERO         ,ZERO);
      xc_SW = Vector2D(-Wedge_Length,ZERO);
    } else if (iBlk == 1) {
      xc_NW = Vector2D(ZERO                                  ,Wedge_Length);
      xc_NE = Vector2D(Wedge_Length*cos(PI*Wedge_Angle/180.0),Wedge_Length);
      xc_SE = Vector2D(Wedge_Length*cos(PI*Wedge_Angle/180.0),Wedge_Length*sin(PI*Wedge_Angle/180.0));
      xc_SW = Vector2D(ZERO                                  ,ZERO);
    }
    // Create the splines defining the north, south, east, and west 
    // boundaries of the rectangular box.
    Create_Spline_Line(Bnd_Spline_North, xc_NW, xc_NE, 2);
    Create_Spline_Line(Bnd_Spline_South, xc_SW, xc_SE, 2);
    Create_Spline_Line(Bnd_Spline_East, xc_SE, xc_NE, 2);
    Create_Spline_Line(Bnd_Spline_West, xc_SW, xc_NW, 2);
    // Set the boundary condition types for each of the boundary splines.
    if (iBlk == 0) {
      Bnd_Spline_North.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_FIXED);
    } else if (iBlk == 1) {
      Bnd_Spline_North.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_West.setBCtype(BC_NONE);
    }
    // Assign values to the stretching function parameters and
    // boundary grid line orthogonality parameters.
    Stretch_I = STRETCHING_FCN_LINEAR;
    Beta_I = ONE; 
    Tau_I = ZERO;
    Stretch_J = STRETCHING_FCN_LINEAR;
    Beta_J = ONE;
    Tau_J = ZERO;
    Orthogonal_North = 0;
    Orthogonal_South = 0;
    Orthogonal_East = 0;
    Orthogonal_West = 0;
    // Create the 2D quadrilateral grid block.
    Create_Quad_Block(Grid_ptr[iBlk][0],
		      Bnd_Spline_North,
		      Bnd_Spline_South,
		      Bnd_Spline_East,
		      Bnd_Spline_West,
		      Number_of_Cells_Idir,
		      Number_of_Cells_Jdir,
		      Number_of_Ghost_Cells,
		      GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
		      Stretch_I,
		      Beta_I, 
		      Tau_I,
		      Stretch_J,
		      Beta_J,
		      Tau_J,
		      Orthogonal_North,
		      Orthogonal_South,
		      Orthogonal_East,
		      Orthogonal_West);
    Smooth_Quad_Block(Grid_ptr[iBlk][0], min(250, 2*max(Number_of_Cells_Idir,Number_of_Cells_Jdir)));
    // Deallocate the memory for the boundary splines.
    Bnd_Spline_North.deallocate();
    Bnd_Spline_South.deallocate();
    Bnd_Spline_East.deallocate();
    Bnd_Spline_West.deallocate();
  }

  // Return the grid.
  return Grid_ptr;

}

/**********************************************************************
 * Routine: Grid_Unsteady_Blunt_Body                                  *
 *                                                                    *
 * Generates a single block quadilateral mesh with clustering for     *
 * predicting supersonic flow around a cirucular cylinder blunt body. *
 *                                                                    *
 * Usage: Grid_ptr = Grid_Unsteady_Blunt_Body(Grid_ptr,               *
 *                                            nblk_i,                 *
 *                                            nblk_j,                 *
 *                                            ONE,                    *
 *                                            FOUR,                   *
 *   	                                      150,                    *
 *  	                                      50,                     *
 *                                            2);                     *
 *                                                                    *
 **********************************************************************/
Grid2D_Quad_Block** Grid_Unsteady_Blunt_Body(Grid2D_Quad_Block **Grid_ptr,
                                             int &Number_of_Blocks_Idir,
                                             int &Number_of_Blocks_Jdir,
					     const double &Radius,
					     const double &Mach_Number,
					     const int Number_of_Cells_Idir,
					     const int Number_of_Cells_Jdir,
					     const int Number_of_Ghost_Cells) {

  int Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D Bnd_Spline_North, Bnd_Spline_South,
           Bnd_Spline_East,  Bnd_Spline_West, Bow_Spline;

  // Allocate memory for grid block.
  Number_of_Blocks_Idir = 1;
  Number_of_Blocks_Jdir = 3;
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
				       Number_of_Blocks_Idir, 
				       Number_of_Blocks_Jdir);

  // Create the mesh for each block representing the complete grid.
  for (int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++) {
    // Create the splines defining the north, south, east, and west 
    // boundaries of the rectangular box.
    if (jBlk == 2) {
      Create_Spline_Bow_Shock(Bow_Spline,
			      Radius,
			      Mach_Number,
			      1,
			      181);
      xc_SW = Bow_Spline.Xp[0] - Vector2D(TWO*Radius,ZERO);
      xc_SE = Vector2D(xc_SW.x,Bow_Spline.Xp[180].y);
      xc_NW = Bow_Spline.Xp[0] - Vector2D(FOUR*Radius,ZERO);
      xc_NE = Vector2D(xc_NW.x,Bow_Spline.Xp[180].y);
      Create_Spline_Line(Bnd_Spline_North,xc_NW,xc_NE,2);
      Create_Spline_Line(Bnd_Spline_South,xc_SW,xc_SE,2);
      Create_Spline_Line(Bnd_Spline_East, xc_SE,xc_NE,2);
      Create_Spline_Line(Bnd_Spline_West, xc_SW,xc_NW,2);
    } else if (jBlk == 1) {
      Create_Spline_Bow_Shock(Bnd_Spline_South,
			      Radius,
			      Mach_Number,
			      1,
			      181);
      xc_SW = Bnd_Spline_South.Xp[0];
      xc_SE = Bnd_Spline_South.Xp[180];
      xc_NW = Bnd_Spline_South.Xp[0] - Vector2D(TWO*Radius,ZERO);
      xc_NE = Vector2D(xc_NW.x,Bnd_Spline_South.Xp[180].y);
      Create_Spline_Line(Bnd_Spline_North,xc_NW,xc_NE,2);
      Create_Spline_Line(Bnd_Spline_East,xc_SE,xc_NE,2);
      Create_Spline_Line(Bnd_Spline_West,xc_SW,xc_NW,2);
    } else if (jBlk == 0) {
      Create_Spline_Bow_Shock(Bnd_Spline_North,
			      Radius,
			      Mach_Number,
			      1,
			      181);
      Create_Spline_Circular_Arc(Bnd_Spline_South,
				 Vector2D(ZERO,ZERO),
				 Radius,
				 180.00,
				 90.00,
				 181);
      xc_SW = Bnd_Spline_South.Xp[0];
      xc_NW = Bnd_Spline_North.Xp[0];
      Create_Spline_Line(Bnd_Spline_West, xc_SW, xc_NW, 2);
      xc_SE = Bnd_Spline_South.Xp[Bnd_Spline_South.np-1];
      xc_NE = Bnd_Spline_North.Xp[Bnd_Spline_North.np-1];
      Create_Spline_Line(Bnd_Spline_East, xc_SE, xc_NE, 2);
    }
    // Set the boundary condition types for each of the boundary splines.
    if (jBlk == 2) {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(BC_NONE);
      Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_West.setBCtype(BC_REFLECTION);
    } else if (jBlk == 1) {
      Bnd_Spline_North.setBCtype(BC_NONE);
      Bnd_Spline_South.setBCtype(BC_NONE);
      Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_West.setBCtype(BC_REFLECTION);
    } else if (jBlk == 0) {
      Bnd_Spline_North.setBCtype(BC_NONE);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_West.setBCtype(BC_REFLECTION);
    }
    // Assign values to the stretching function parameters and
    // boundary grid line orthogonality parameters.
    if (jBlk == 2 || jBlk == 1) {
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ONE; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;
      Beta_J = ONE;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;
    } else if (jBlk == 0) {
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ZERO; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = 1.05;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 1;
      Orthogonal_East = 1;
      Orthogonal_West = 1;
    }
    // Create the 2D quadrilateral grid block.
    Create_Quad_Block(Grid_ptr[0][jBlk],
		      Bnd_Spline_North,
		      Bnd_Spline_South,
		      Bnd_Spline_East,
		      Bnd_Spline_West,
		      Number_of_Cells_Idir,
		      Number_of_Cells_Jdir,
		      Number_of_Ghost_Cells,
		      GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
		      Stretch_I,
		      Beta_I, 
		      Tau_I,
		      Stretch_J,
		      Beta_J,
		      Tau_J,
		      Orthogonal_North,
		      Orthogonal_South,
		      Orthogonal_East,
		      Orthogonal_West);
    // Deallocate the memory for the boundary splines.
    Bnd_Spline_North.deallocate();
    Bnd_Spline_South.deallocate();
    Bnd_Spline_East.deallocate();
    Bnd_Spline_West.deallocate();
    // Smooth the 2D quadrilateral grid block.
    Smooth_Quad_Block(Grid_ptr[0][jBlk],
		      min(250, 2*max(Number_of_Cells_Idir,Number_of_Cells_Jdir)));
  }

  if (Bow_Spline.np > 0) Bow_Spline.deallocate();

  // Return the grid.
  return Grid_ptr;

}

/**********************************************************************
 * Routine: Grid_Ringleb_Flow                                         *
 *                                                                    *
 * Generates a uniform 2D mesh for Ringleb's flow.                    *
 *                                                                    *
 * Usage: Grid_ptr = Grid_Ringleb_Flow(Grid_ptr,                      *
 *                                     nblk_i,                        *
 *                                     nblk_j,                        *
 *   	                               ncells_i,                      *
 *  	                               ncells_j,                      *
 *                                     Nghost_Cells);                 *
 *                                                                    *
 **********************************************************************/
Grid2D_Quad_Block** Grid_Ringleb_Flow(Grid2D_Quad_Block **Grid_ptr,
                                      int &Number_of_Blocks_Idir,
                                      int &Number_of_Blocks_Jdir,
				      const int Number_of_Cells_Idir,
				      const int Number_of_Cells_Jdir,
				      const int Number_of_Ghost_Cells,
				      const int NumOfIter_UnsmoothMesh,
				      const int Stretch_I, const int Stretch_J,
				      const double Beta_I, const double Beta_J,
				      const double Tau_I, const double Tau_J) {

  int nk, nq,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D Bnd_Spline_North, Bnd_Spline_South,
           Bnd_Spline_East, Bnd_Spline_West;
  double **rho;
  double delta_q;
  double delta_k;
  double k_init, q_init, q_final;
  double **q, **k, qo, ko, c, J;
  double g = 1.40;
  Vector2D norm_dir, X_norm, X_tan;

  // Allocate memory for grid block.
  Number_of_Blocks_Idir = 1;
  Number_of_Blocks_Jdir = 1;
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
				       Number_of_Blocks_Idir, 
				       Number_of_Blocks_Jdir);

  nk = 32;
  nq = 50;

  // Create the mesh for each block representing the complete grid.
  k = new double*[nk];
  for (int i = 0; i < nk; ++i) k[i] = new double[nq];
  q = new double*[nk];
  for (int i = 0; i < nk; ++i) q[i] = new double[nq];
  rho = new double*[nk];
  for (int i = 0; i < nk; ++i) rho[i] = new double[nq];

  delta_k = (1.50-0.75)/double(nk-1);
  k_init = 0.75;
  q_init = 0.50;
  for (int i = 0; i < nk; ++i){
    ko = k_init + double(i)*delta_k;
    q_final = ko; // condition y = 0
    delta_q = (q_final - q_init)/double(nq-1);
    for (int j = 0; j < nq; ++j) {
      qo = q_init + double(j)*delta_q;
      k[i][j] = ko;
      q[i][j] = qo;
    }
  }

  Bnd_Spline_North.allocate(nk);  Bnd_Spline_North.settype(SPLINE2D_QUINTIC);
  Bnd_Spline_South.allocate(nk);  Bnd_Spline_South.settype(SPLINE2D_QUINTIC);
  Bnd_Spline_East.allocate(nq);	  Bnd_Spline_East.settype(SPLINE2D_QUINTIC);
  Bnd_Spline_West.allocate(nq);	  Bnd_Spline_West.settype(SPLINE2D_QUINTIC);

  for (int i = 0; i < nk; ++i) {
    for (int j = 0; j < nq; ++j){

      c = sqrt(ONE - ((g-ONE)/TWO)*q[i][j]*q[i][j]);
      rho[i][j] = pow(c,TWO/(g-ONE));
      J = ONE/c + ONE/(THREE*pow(c,THREE)) + ONE/(FIVE*pow(c,FIVE)) - HALF*log((ONE+c)/(ONE-c));
      // NORTH spline.
      if (j == 0) {
	Bnd_Spline_North.Xp[nk-1-i].x = (HALF/rho[i][j])*(TWO/(k[i][j]*k[i][j]) - ONE/(q[i][j]*q[i][j])) - HALF*J;
	Bnd_Spline_North.Xp[nk-1-i].y = (ONE/(k[i][j]*rho[i][j]*q[i][j]))*sqrt(ONE - (q[i][j]*q[i][j])/(k[i][j]*k[i][j]));
	Bnd_Spline_North.bc[nk-1-i] = BC_RINGLEB_FLOW;
	if (i == 0 || i == nk-1) Bnd_Spline_North.tp[nk-1-i] = SPLINE2D_POINT_SHARP_CORNER;
	else Bnd_Spline_North.tp[nk-1-i] = SPLINE2D_POINT_NORMAL;
      }
      // SOUTH spline.
      if (j == nq-1) {
	Bnd_Spline_South.Xp[nk-1-i].x = (HALF/rho[i][j])*(TWO/(k[i][j]*k[i][j]) - ONE/(q[i][j]*q[i][j])) - HALF*J;
	Bnd_Spline_South.Xp[nk-1-i].y = ZERO;
	Bnd_Spline_South.bc[nk-1-i] = BC_RINGLEB_FLOW;
	if (i == 0 || i == nk-1) Bnd_Spline_South.tp[nk-1-i] = SPLINE2D_POINT_SHARP_CORNER;
	else Bnd_Spline_South.tp[nk-1-i] = SPLINE2D_POINT_NORMAL;
      }
      // EAST spline.
      if (i == 0) {
	Bnd_Spline_East.Xp[nq-1-j].x = (HALF/rho[i][j])*(TWO/(k[i][j]*k[i][j]) - ONE/(q[i][j]*q[i][j])) - HALF*J;
	Bnd_Spline_East.Xp[nq-1-j].y = (ONE/(k[i][j]*rho[i][j]*q[i][j]))*sqrt(ONE - (q[i][j]*q[i][j])/(k[i][j]*k[i][j]));
	Bnd_Spline_East.bc[nq-1-j] = BC_REFLECTION;
	if (j == 0 || j == nq-1) Bnd_Spline_East.tp[nq-1-j] = SPLINE2D_POINT_SHARP_CORNER;
	else Bnd_Spline_East.tp[nq-1-j] = SPLINE2D_POINT_NORMAL;
      }
      // WEST spline.
      if (i == nk-1) {
	Bnd_Spline_West.Xp[nq-1-j].x = (HALF/rho[i][j])*(TWO/(k[i][j]*k[i][j]) - ONE/(q[i][j]*q[i][j])) - HALF*J;
	Bnd_Spline_West.Xp[nq-1-j].y = (ONE/(k[i][j]*rho[i][j]*q[i][j]))*sqrt(ONE - (q[i][j]*q[i][j])/(k[i][j]*k[i][j]));
	Bnd_Spline_West.bc[nq-1-j] = BC_REFLECTION;
	if (j == 0 || j == nq-1) Bnd_Spline_West.tp[nq-1-j] = SPLINE2D_POINT_SHARP_CORNER;
	else Bnd_Spline_West.tp[nq-1-j] = SPLINE2D_POINT_NORMAL;
      }

    }
  }

  Bnd_Spline_North.pathlength();
  Bnd_Spline_South.pathlength();
  Bnd_Spline_East.pathlength();
  Bnd_Spline_West.pathlength();

  // Set the boundary condition types for each of the boundary 
  // splines.
  Bnd_Spline_North.setBCtype(BC_RINGLEB_FLOW);
  Bnd_Spline_South.setBCtype(BC_RINGLEB_FLOW);
  //   Bnd_Spline_East.setBCtype(BC_REFLECTION);
  //   Bnd_Spline_West.setBCtype(BC_REFLECTION);
  Bnd_Spline_East.setBCtype(BC_RINGLEB_FLOW);
  Bnd_Spline_West.setBCtype(BC_RINGLEB_FLOW);

  // Assign values to the stretching function parameters and
  // boundary grid line orthogonality parameters.

  Orthogonal_North = 1;
  Orthogonal_South = 1;
  Orthogonal_East = 1;
  Orthogonal_West = 1;

  // Create the 2D quadrilateral grid block representing the mesh.
  Create_Quad_Block(Grid_ptr[0][0],
		    Bnd_Spline_North,
		    Bnd_Spline_South,
		    Bnd_Spline_East,
		    Bnd_Spline_West,
		    Number_of_Cells_Idir,
		    Number_of_Cells_Jdir,
		    Number_of_Ghost_Cells,
		    GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST,
		    Stretch_I,
		    Beta_I,
		    Tau_I,
		    Stretch_J,
		    Beta_J,
		    Tau_J,
		    Orthogonal_North,
		    Orthogonal_South,
		    Orthogonal_East,
		    Orthogonal_West);

  /* Smooth the 2D quadrilateral grid block. */
  Smooth_Quad_Block(Grid_ptr[0][0],min(250,2*max(Number_of_Cells_Idir,Number_of_Cells_Jdir)));

  /* Unsmooth the 2D quadrilateral grid block. */
  Unsmooth_Interior_Nodes(Grid_ptr[0][0],NumOfIter_UnsmoothMesh);

  // Deallocate the memory for the boundary splines.
  Bnd_Spline_North.deallocate();
  Bnd_Spline_South.deallocate();
  Bnd_Spline_East.deallocate();
  Bnd_Spline_West.deallocate();

  // Deallocate memory for point, kq, and rho.
  for (int i = 0; i < nk; i++) {
    delete []k[i];   k[i]   = NULL;
    delete []q[i];   q[i]   = NULL;
    delete []rho[i]; rho[i] = NULL;
  }
  delete []k;   k   = NULL;
  delete []q;   q   = NULL;
  delete []rho; rho = NULL;

  // Return the grid.
  return Grid_ptr;

}

/**********************************************************************
 * Routine: Grid_Bump_Channel_Flow                                    *
 *                                                                    *
 * Generates a single block quadilateral mesh with clustering for     *
 * predicting supersonic flow around a cirucular cylinder blunt body. *
 *                                                                    *
 * Usage: Grid_ptr = Grid_ump_Channel_Flow(Grid_ptr,                  *
 *                                         nblk_i,                    *
 *                                         nblk_j,                    *
 *                                         ONE,                       *
 *   	                                   150,                       *
 *  	                                   50,                        *
 *                                         2);                        *
 *                                                                    *
 **********************************************************************/
Grid2D_Quad_Block** Grid_Bump_Channel_Flow(Grid2D_Quad_Block **Grid_ptr,
					   int &Number_of_Blocks_Idir,
					   int &Number_of_Blocks_Jdir,
					   const double &Radius,
					   const int Number_of_Cells_Idir,
					   const int Number_of_Cells_Jdir,
					   const int Number_of_Ghost_Cells) {

  int Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D x1, x2;
  Spline2D Bnd_Spline_North, Bnd_Spline_South,
           Bnd_Spline_East, Bnd_Spline_West;
  double R, Theta;

  // Allocate memory for grid block.
  Number_of_Blocks_Idir = 4;
  Number_of_Blocks_Jdir = 2;
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
				       Number_of_Blocks_Idir, 
				       Number_of_Blocks_Jdir);

  // Create the mesh for each block representing the complete grid.
  for (int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++) {
    for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

      if (iBlk == 0 && jBlk == 0) {
	// Create the splines defining the north, south, east, and west 
	// boundaries of the grid.
	x1 = Vector2D(-1.0,1.0);
	x2 = Vector2D( 0.0,1.0);
	Create_Spline_Line(Bnd_Spline_North,x1,x2,2);
	x1 = Vector2D(-1.0,0.0);
	x2 = Vector2D( 0.0,0.0);
	Create_Spline_Line(Bnd_Spline_South,x1,x2,2);
	x1 = Vector2D( 0.0,0.0);
	x2 = Vector2D( 0.0,1.0);
	Create_Spline_Line(Bnd_Spline_East,x1,x2,2);
	x1 = Vector2D(-1.0,0.0);
	x2 = Vector2D(-1.0,1.0);
	Create_Spline_Line(Bnd_Spline_West,x1,x2,2);
	// Set the boundary condition types for each of the boundary splines.
	Bnd_Spline_North.setBCtype(BC_NONE);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	Bnd_Spline_East.setBCtype(BC_NONE);
	Bnd_Spline_West.setBCtype(BC_FIXED);//BC_CHARACTERISTIC);//BC_CONSTANT_EXTRAPOLATION);
      } else if (iBlk == 1 && jBlk == 0) {
	// Create the splines defining the north, south, east, and west 
	// boundaries of the grid.
	x1 = Vector2D( 0.0,1.0);
	x2 = Vector2D( 1.0,1.0);
	Create_Spline_Line(Bnd_Spline_North,x1,x2,2);
	R = (0.25 + 0.042*0.042)/(TWO*0.042);
	Theta = acos((R - 0.042)/R);
	Create_Spline_Circular_Arc(Bnd_Spline_South,
				   Vector2D_ZERO,
				   R,
				   ZERO,
				   -TWO*Theta*180.0/PI,
				   31);
	Rotate_Spline(Bnd_Spline_South,HALF*PI+Theta);
	Translate_Spline(Bnd_Spline_South,Vector2D(HALF,0.042-R));
	x1 = Vector2D( 1.0,0.0);
	x2 = Vector2D( 1.0,1.0);
	Create_Spline_Line(Bnd_Spline_East,x1,x2,2);
	x1 = Vector2D( 0.0,0.0);
	x2 = Vector2D( 0.0,1.0);
	Create_Spline_Line(Bnd_Spline_West,x1,x2,2);
	// Set the boundary condition types for each of the boundary splines.
	Bnd_Spline_North.setBCtype(BC_NONE);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	Bnd_Spline_East.setBCtype(BC_NONE);
	Bnd_Spline_West.setBCtype(BC_NONE);
      } else if (iBlk == 2 && jBlk == 0) {
	// Create the splines defining the north, south, east, and west 
	// boundaries of the grid.
	x1 = Vector2D( 1.00,1.00);
	x2 = Vector2D( 2.75,1.00);
	Create_Spline_Line(Bnd_Spline_North,x1,x2,2);
	x1 = Vector2D( 1.00,0.00);
	x2 = Vector2D( 2.75,0.00);
	Create_Spline_Line(Bnd_Spline_South,x1,x2,2);
	x1 = Vector2D( 2.75,0.00);
	x2 = Vector2D( 2.75,1.00);
	Create_Spline_Line(Bnd_Spline_East,x1,x2,2);
	x1 = Vector2D( 1.00,0.00);
	x2 = Vector2D( 1.00,1.00);
	Create_Spline_Line(Bnd_Spline_West,x1,x2,2);
	// Set the boundary condition types for each of the boundary splines.
	Bnd_Spline_North.setBCtype(BC_NONE);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	Bnd_Spline_East.setBCtype(BC_NONE);
	Bnd_Spline_West.setBCtype(BC_NONE);
      } else if (iBlk == 3 && jBlk == 0) {
	// Create the splines defining the north, south, east, and west 
	// boundaries of the grid.
	x1 = Vector2D( 2.75,1.00);
	x2 = Vector2D( 4.50,1.00);
	Create_Spline_Line(Bnd_Spline_North,x1,x2,2);
	x1 = Vector2D( 2.75,0.00);
	x2 = Vector2D( 4.50,0.00);
	Create_Spline_Line(Bnd_Spline_South,x1,x2,2);
	x1 = Vector2D( 4.50,0.00);
	x2 = Vector2D( 4.50,1.00);
	Create_Spline_Line(Bnd_Spline_East,x1,x2,2);
	x1 = Vector2D( 2.75,0.00);
	x2 = Vector2D( 2.75,1.00);
	Create_Spline_Line(Bnd_Spline_West,x1,x2,2);
	// Set the boundary condition types for each of the boundary splines.
	Bnd_Spline_North.setBCtype(BC_NONE);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	Bnd_Spline_East.setBCtype(BC_CHARACTERISTIC);//BC_CONSTANT_EXTRAPOLATION);
	Bnd_Spline_West.setBCtype(BC_NONE);
      } else if (iBlk == 0 && jBlk == 1) {
	// Create the splines defining the north, south, east, and west 
	// boundaries of the grid.
	x1 = Vector2D(-1.00,2.00);
	x2 = Vector2D( 0.00,2.00);
	Create_Spline_Line(Bnd_Spline_North,x1,x2,2);
	x1 = Vector2D(-1.00,1.00);
	x2 = Vector2D( 0.00,1.00);
	Create_Spline_Line(Bnd_Spline_South,x1,x2,2);
	x1 = Vector2D( 0.00,1.00);
	x2 = Vector2D( 0.00,2.00);
	Create_Spline_Line(Bnd_Spline_East,x1,x2,2);
	x1 = Vector2D(-1.00,1.00);
	x2 = Vector2D(-1.00,2.00);
	Create_Spline_Line(Bnd_Spline_West,x1,x2,2);
	// Set the boundary condition types for each of the boundary splines.
	Bnd_Spline_North.setBCtype(BC_REFLECTION);
	Bnd_Spline_South.setBCtype(BC_NONE);
	Bnd_Spline_East.setBCtype(BC_NONE);
	Bnd_Spline_West.setBCtype(BC_FIXED);//BC_CHARACTERISTIC);//BC_CONSTANT_EXTRAPOLATION);
      } else if (iBlk == 1 && jBlk == 1) {
	// Create the splines defining the north, south, east, and west 
	// boundaries of the grid.
	x1 = Vector2D( 0.00,2.00);
	x2 = Vector2D( 1.00,2.00);
	Create_Spline_Line(Bnd_Spline_North,x1,x2,2);
	x1 = Vector2D( 0.00,1.00);
	x2 = Vector2D( 1.00,1.00);
	Create_Spline_Line(Bnd_Spline_South,x1,x2,2);
	x1 = Vector2D( 1.00,1.00);
	x2 = Vector2D( 1.00,2.00);
	Create_Spline_Line(Bnd_Spline_East,x1,x2,2);
	x1 = Vector2D( 0.00,1.00);
	x2 = Vector2D( 0.00,2.00);
	Create_Spline_Line(Bnd_Spline_West,x1,x2,2);
	// Set the boundary condition types for each of the boundary splines.
	Bnd_Spline_North.setBCtype(BC_REFLECTION);
	Bnd_Spline_South.setBCtype(BC_NONE);
	Bnd_Spline_East.setBCtype(BC_NONE);
	Bnd_Spline_West.setBCtype(BC_NONE);
      } else if (iBlk == 2 && jBlk == 1) {
	// Create the splines defining the north, south, east, and west 
	// boundaries of the grid.
	x1 = Vector2D( 1.00,2.00);
	x2 = Vector2D( 2.75,2.00);
	Create_Spline_Line(Bnd_Spline_North,x1,x2,2);
	x1 = Vector2D( 1.00,1.00);
	x2 = Vector2D( 2.75,1.00);
	Create_Spline_Line(Bnd_Spline_South,x1,x2,2);
	x1 = Vector2D( 2.75,1.00);
	x2 = Vector2D( 2.75,2.00);
	Create_Spline_Line(Bnd_Spline_East,x1,x2,2);
	x1 = Vector2D( 1.00,1.00);
	x2 = Vector2D( 1.00,2.00);
	Create_Spline_Line(Bnd_Spline_West,x1,x2,2);
	// Set the boundary condition types for each of the boundary splines.
	Bnd_Spline_North.setBCtype(BC_REFLECTION);
	Bnd_Spline_South.setBCtype(BC_NONE);
	Bnd_Spline_East.setBCtype(BC_NONE);
	Bnd_Spline_West.setBCtype(BC_NONE);
      } else if (iBlk == 3 && jBlk == 1) {
	// Create the splines defining the north, south, east, and west 
	// boundaries of the grid.
	x1 = Vector2D( 2.75,2.00);
	x2 = Vector2D( 4.50,2.00);
	Create_Spline_Line(Bnd_Spline_North,x1,x2,2);
	x1 = Vector2D( 2.75,1.00);
	x2 = Vector2D( 4.50,1.00);
	Create_Spline_Line(Bnd_Spline_South,x1,x2,2);
	x1 = Vector2D( 4.50,1.00);
	x2 = Vector2D( 4.50,2.00);
	Create_Spline_Line(Bnd_Spline_East,x1,x2,2);
	x1 = Vector2D( 2.75,1.00);
	x2 = Vector2D( 2.75,2.00);
	Create_Spline_Line(Bnd_Spline_West,x1,x2,2);
	// Set the boundary condition types for each of the boundary splines.
	Bnd_Spline_North.setBCtype(BC_REFLECTION);
	Bnd_Spline_South.setBCtype(BC_NONE);
	Bnd_Spline_East.setBCtype(BC_CHARACTERISTIC);//BC_CONSTANT_EXTRAPOLATION);
	Bnd_Spline_West.setBCtype(BC_NONE);
      }

      // Assign values to the stretching function parameters and boundary 
      // grid line orthogonality parameters.
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ZERO; 
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;
      Beta_J = ZERO;
      Tau_J = ZERO;
      Orthogonal_North = 1;
      Orthogonal_South = 1;
      Orthogonal_East = 1;
      Orthogonal_West = 1;

      // Create the 2D quadrilateral grid block representing the mesh.
      Create_Quad_Block(Grid_ptr[iBlk][jBlk],
			Bnd_Spline_North,
			Bnd_Spline_South,
			Bnd_Spline_East,
			Bnd_Spline_West,
			Number_of_Cells_Idir,
			Number_of_Cells_Jdir,
			Number_of_Ghost_Cells,
			GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH,
			Stretch_I,
			Beta_I, 
			Tau_I,
			Stretch_J,
			Beta_J,
			Tau_J,
			Orthogonal_North,
			Orthogonal_South,
			Orthogonal_East,
			Orthogonal_West);

      // Smooth the 2D quadrilateral grid block.
      Smooth_Quad_Block(Grid_ptr[iBlk][jBlk],min(250,2*max(Number_of_Cells_Idir,Number_of_Cells_Jdir)));

      // Deallocate the memory for the boundary splines.
      Bnd_Spline_North.deallocate();
      Bnd_Spline_South.deallocate();
      Bnd_Spline_East.deallocate();
      Bnd_Spline_West.deallocate();

    }
  }

  // Return the grid.
  return Grid_ptr;

}
