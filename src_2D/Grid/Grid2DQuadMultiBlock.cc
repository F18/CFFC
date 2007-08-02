/* Grid2DQuadMultiBlock.cc:  Multi-block subroutines for 
                             2D quadrilateral block grid class. */

/* Include 2D quadrilateral block grid type header file. */

#ifndef _GRID2D_QUAD_BLOCK_INCLUDED
#include "Grid2DQuad.h"
#endif // _GRID2D_QUAD_BLOCK_INCLUDED

// Include NASA rotor 37 and 67 header files.

#ifndef _NASA_ROTOR37_INCLUDED
#include "../Grid/NASARotor37.h"
#endif // _NASA_ROTOR37_INCLUDED

#ifndef _NASA_ROTOR67_INCLUDED
#include "../Grid/NASARotor67.h"
#endif // _NASA_ROTOR67_INCLUDED

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
 
    /* Deallocate memory. */

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

//            if (Grid_ptr[i][j].Node != NULL) {
//               Smooth_Quad_Block(Grid_ptr[i][j],
//                                 2*max(Grid_ptr[i][j].NCi, 
//                                       Grid_ptr[i][j].NCj));
//            } /* endif */

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
					 const int Number_of_Ghost_Cells) {

    int iBlk, jBlk, n_cells_i, n_cells_j, 
        Stretch_I, Stretch_J,
        Orthogonal_North, Orthogonal_South,
        Orthogonal_East, Orthogonal_West;
    double Beta_I, Tau_I, Beta_J, Tau_J;
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

          /* Assign values to the stretching function parameters
             and boundary grid line orthogonality parameters. */

          Stretch_I = STRETCHING_FCN_LINEAR;
          Beta_I = ZERO; 
          Tau_I = ZERO;
          Stretch_J = STRETCHING_FCN_LINEAR;
          Beta_J = ZERO;
          Tau_J = ZERO;
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

Grid2D_Quad_Block** Grid_Rectangular_Box(Grid2D_Quad_Block **Grid_ptr,
                                         int &Number_of_Blocks_Idir,
                                         int &Number_of_Blocks_Jdir,
                                         const double &Width,
                                         const double &Height,
					 const int &Stretching_Flag,
					 const int &Stretching_Type_Idir,
					 const int &Stretching_Type_Jdir,
					 const double &Stretching_Factor_Idir,
					 const double &Stretching_Factor_Jdir,
 	                                 const int Number_of_Cells_Idir,
	                                 const int Number_of_Cells_Jdir,
					 const int Number_of_Ghost_Cells) {

  int iBlk, jBlk, n_cells_i, n_cells_j,
      Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
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

      /* Assign values to the stretching function parameters
	 and boundary grid line orthogonality parameters. */

      if (Stretching_Flag) {
	Stretch_I = Stretching_Type_Idir;
	Stretch_J = Stretching_Type_Jdir;
	Beta_I = Stretching_Factor_Idir;
	Beta_J = Stretching_Factor_Jdir;
      } else {
	Stretch_I = STRETCHING_FCN_LINEAR;
	Stretch_J = STRETCHING_FCN_LINEAR;
	Beta_I = ZERO; 
	Beta_J = ZERO;
      }
      Tau_I = ZERO;
      Tau_J = ZERO;
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
 *         		             100,               *
 *         		             2);                *
 *                                                      *
 ********************************************************/
Grid2D_Quad_Block** Grid_Flat_Plate4(Grid2D_Quad_Block **Grid_ptr,
				     int &Number_of_Blocks_Idir,
				     int &Number_of_Blocks_Jdir,
				     const double &Length,
				     const int &Stretching_Flag,
				     const double &Stretching_Factor_Idir,
				     const double &Stretching_Factor_Jdir,
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

  // Allocate memory for grid blocks.  There are two grid blocks for 
  // this mesh.
  Number_of_Blocks_Idir = 3;
  Number_of_Blocks_Jdir = 1;
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr,
				       Number_of_Blocks_Idir,
				       Number_of_Blocks_Jdir);

  // Create the mesh for each block representing the complete grid.
  for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

    // Assign values to the locations of the corners of the 
    // rectangular box shaped domain representing each of the 
    // blocks in the grid.
    if (iBlk == 0) {
      xc_NW = Vector2D(-4.0*Length,4.0*Length);
      xc_NE = Vector2D(ZERO,4.0*Length);
      xc_SE = Vector2D(ZERO,ZERO);
      xc_SW = Vector2D(-4.0*Length,ZERO);
    } else if (iBlk == 1) {
      xc_NW = Vector2D(ZERO,4.0*Length);
      xc_NE = Vector2D(Length,4.0*Length);
      xc_SE = Vector2D(Length,ZERO);
      xc_SW = Vector2D(ZERO,ZERO);
    } else {
      xc_NW = Vector2D(Length,4.0*Length);
      xc_NE = Vector2D(5.0*Length,4.0*Length);
      xc_SE = Vector2D(5.0*Length,ZERO);
      xc_SW = Vector2D(Length,ZERO);
    }

    // Create the splines defining the north, south, east, and west 
    // boundaries of the rectangular boxes.
    Create_Spline_Line(Bnd_Spline_North,xc_NW,xc_NE,2);
    Create_Spline_Line(Bnd_Spline_South,xc_SW,xc_SE,2);
    Create_Spline_Line(Bnd_Spline_East,xc_SE,xc_NE,2);
    Create_Spline_Line(Bnd_Spline_West,xc_SW,xc_NW,2);

    // Set the boundary condition types for each of the boundary splines.
    if (iBlk == 0) {
      Bnd_Spline_North.setBCtype(BC_FIXED);//BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_FIXED);
    } else if (iBlk == 1) {
      Bnd_Spline_North.setBCtype(BC_FIXED);//BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_South.setBCtype(BC_WALL_VISCOUS_HEATFLUX);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_NONE);
    } else {
      Bnd_Spline_North.setBCtype(BC_FIXED);//BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_FIXED);//BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_West.setBCtype(BC_NONE);
    }

    // Determine the number of cells for this block.
    n_cells_i = Number_of_Cells_Idir;
    n_cells_j = Number_of_Cells_Jdir;

    // Assign values to the stretching function parameters and boundary
    // grid line orthogonality parameters.
    if (iBlk == 0) {
      Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
      Beta_I = ONE + (Stretching_Factor_Idir-ONE)/6.50;
    } else if (iBlk == 1) {
      Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
      Beta_I = Stretching_Factor_Idir;
    } else {
      Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_I = ONE + (Stretching_Factor_Idir-ONE)/6.50;
    }
    Tau_I = ZERO;
    Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
    Beta_J = Stretching_Factor_Jdir;
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

    // Deallocate the memory for the boundary splines.
    Bnd_Spline_North.deallocate();
    Bnd_Spline_South.deallocate();
    Bnd_Spline_East.deallocate();
    Bnd_Spline_West.deallocate();
    
  }

  // Return the grid.
  return Grid_ptr;

}
Grid2D_Quad_Block** Grid_Flat_Plate3(Grid2D_Quad_Block **Grid_ptr,
				    int &Number_of_Blocks_Idir,
				    int &Number_of_Blocks_Jdir,
				    const double &Length,
				    const int &Stretching_Flag,
				    const double &Stretching_Factor_Idir,
				    const double &Stretching_Factor_Jdir,
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

  // Allocate memory for grid blocks.  There are two grid blocks for 
  // this mesh.
  Number_of_Blocks_Idir = 2;
  Number_of_Blocks_Jdir = 1;
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr,
				       Number_of_Blocks_Idir,
				       Number_of_Blocks_Jdir);

  // Create the mesh for each block representing the complete grid.
  for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

    // Assign values to the locations of the corners of the 
    // rectangular box shaped domain representing each of the 
    // blocks in the grid.
    if (iBlk == 0) {
      xc_NW = FOUR*Vector2D(-Length,Length);
      xc_NE = FOUR*Vector2D(ZERO,Length);
      xc_SE = FOUR*Vector2D(ZERO,ZERO);
      xc_SW = FOUR*Vector2D(-Length,ZERO);
    } else {
      xc_NW = Vector2D(ZERO,FOUR*Length);
      xc_NE = Vector2D(Length,FOUR*Length);
      xc_SE = Vector2D(Length,ZERO);
      xc_SW = Vector2D(ZERO,  ZERO);
    }

    // Create the splines defining the north, south, east, and west 
    // boundaries of the rectangular boxes.
    Create_Spline_Line(Bnd_Spline_North,xc_NW,xc_NE,2);
    Create_Spline_Line(Bnd_Spline_South,xc_SW,xc_SE,2);
    Create_Spline_Line(Bnd_Spline_East,xc_SE,xc_NE,2);
    Create_Spline_Line(Bnd_Spline_West,xc_SW,xc_NW,2);

    // Set the boundary condition types for each of the boundary splines.
    if (iBlk == 0) {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_FIXED);
    } else {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(BC_WALL_VISCOUS_HEATFLUX);
      Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_West.setBCtype(BC_NONE);
    }

    // Determine the number of cells for this block.
    n_cells_i = Number_of_Cells_Idir;///2;
    n_cells_j = Number_of_Cells_Jdir;

    // Assign values to the stretching function parameters and boundary
    // grid line orthogonality parameters.
    Stretch_I = STRETCHING_FCN_LINEAR;
    Beta_I = ZERO;
    Tau_I = ZERO;
    Stretch_J = STRETCHING_FCN_LINEAR;
    Beta_J = ZERO;
    Tau_J = ZERO;
    if (Stretching_Flag) {
      if (iBlk == 0) {
	Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_I = ONE + (Stretching_Factor_Idir-ONE)/6.50;
      } else {
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = Stretching_Factor_Idir;
      }
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = Stretching_Factor_Jdir;
    }
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

    // Deallocate the memory for the boundary splines.
    Bnd_Spline_North.deallocate();
    Bnd_Spline_South.deallocate();
    Bnd_Spline_East.deallocate();
    Bnd_Spline_West.deallocate();
    
  }

  // Return the grid.
  return Grid_ptr;

}
Grid2D_Quad_Block** Grid_Flat_Plate9(Grid2D_Quad_Block **Grid_ptr,
				    int &Number_of_Blocks_Idir,
				    int &Number_of_Blocks_Jdir,
				    const double &Length,
				    const int &Flat_Plate_BC_Type,
				    const int &Stretching_Flag,
				    const double &Stretching_Factor_Idir,
				    const double &Stretching_Factor_Jdir,
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

  // Allocate memory for grid blocks.  There are two grid blocks for 
  // this mesh.
  Number_of_Blocks_Idir = 2;
  Number_of_Blocks_Jdir = 1;
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr,
				       Number_of_Blocks_Idir,
				       Number_of_Blocks_Jdir);

  // Create the mesh for each block representing the complete grid.
  for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

    // Assign values to the locations of the corners of the 
    // rectangular box shaped domain representing each of the 
    // blocks in the grid.
    if (iBlk == 0) {
      xc_NW = TWO*Vector2D(-Length,Length);
      xc_NE = TWO*Vector2D(ZERO,Length);
      xc_SE = TWO*Vector2D(ZERO,ZERO);
      xc_SW = TWO*Vector2D(-Length,ZERO);
    } else {
      xc_NW = TWO*Vector2D(ZERO,Length);
      xc_NE = TWO*Vector2D(Length,Length);
      xc_SE = TWO*Vector2D(Length,ZERO);
      xc_SW = TWO*Vector2D(ZERO,  ZERO);
    }

    // Create the splines defining the north, south, east, and west 
    // boundaries of the rectangular boxes.
    Create_Spline_Line(Bnd_Spline_North,xc_NW,xc_NE,2);
    Create_Spline_Line(Bnd_Spline_South,xc_SW,xc_SE,2);
    Create_Spline_Line(Bnd_Spline_East,xc_SE,xc_NE,2);
    Create_Spline_Line(Bnd_Spline_West,xc_SW,xc_NW,2);

    // Set the boundary condition types for each of the boundary splines.
    if (iBlk == 0) {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_FIXED);
    } else {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(Flat_Plate_BC_Type);
      Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_West.setBCtype(BC_NONE);
    }

    // Determine the number of cells for this block.
    n_cells_i = Number_of_Cells_Idir;///2;
    n_cells_j = Number_of_Cells_Jdir;

    // Assign values to the stretching function parameters and boundary
    // grid line orthogonality parameters.
    Stretch_I = STRETCHING_FCN_LINEAR;
    Beta_I = ZERO;
    Tau_I = ZERO;
    Stretch_J = STRETCHING_FCN_LINEAR;
    Beta_J = ZERO;
    Tau_J = ZERO;
//     if (Stretching_Flag) {
//       if (iBlk == 0) {
// 	Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
// 	Beta_I = ONE + (Stretching_Factor_Idir-ONE)/6.50;
//       } else {
// 	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
// 	Beta_I = Stretching_Factor_Idir;
//       }
//       Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
//       Beta_J = Stretching_Factor_Jdir;
//     }
    if (Stretching_Flag) {
      if (iBlk == 0) Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
      else Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_I = Stretching_Factor_Idir;
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = Stretching_Factor_Jdir;
    }
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

    // Deallocate the memory for the boundary splines.
    Bnd_Spline_North.deallocate();
    Bnd_Spline_South.deallocate();
    Bnd_Spline_East.deallocate();
    Bnd_Spline_West.deallocate();

  }

  // Return the grid.
  return Grid_ptr;

}
Grid2D_Quad_Block** Grid_Flat_Plate(Grid2D_Quad_Block **Grid_ptr,
				    int &Number_of_Blocks_Idir,
				    int &Number_of_Blocks_Jdir,
				    const double &Length,
				    const int &Flat_Plate_BC_Type,
				    const int &Stretching_Flag,
				    const double &Stretching_Factor_Idir,
				    const double &Stretching_Factor_Jdir,
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

  // Allocate memory for grid blocks.  There are two grid blocks for 
  // this mesh.
  if (Flat_Plate_BC_Type == BC_BURNING_SURFACE) Number_of_Blocks_Idir = 3;
  else Number_of_Blocks_Idir = 2;
  Number_of_Blocks_Jdir = 1;
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr,
				       Number_of_Blocks_Idir,
				       Number_of_Blocks_Jdir);

  // Create the mesh for each block representing the complete grid.
  for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

    // Assign values to the locations of the corners of the 
    // rectangular box shaped domain representing each of the 
    // blocks in the grid.
    if (iBlk == 0) {
      xc_NW = TWO*Vector2D(-Length,Length);
      xc_NE = TWO*Vector2D(ZERO,Length);
      xc_SE = TWO*Vector2D(ZERO,ZERO);
      xc_SW = TWO*Vector2D(-Length,ZERO);
    } else if (iBlk == 1) {
      xc_NW = TWO*Vector2D(ZERO,Length);
      xc_NE = TWO*Vector2D(Length,Length);
      xc_SE = TWO*Vector2D(Length,ZERO);
      xc_SW = TWO*Vector2D(ZERO,  ZERO);
    } else {
      xc_NW = TWO*Vector2D(Length,Length);
      xc_NE = TWO*Vector2D(TWO*Length,Length);
      xc_SE = TWO*Vector2D(TWO*Length,ZERO);
      xc_SW = TWO*Vector2D(Length,ZERO);
    }

    // Create the splines defining the north, south, east, and west 
    // boundaries of the rectangular boxes.
    Create_Spline_Line(Bnd_Spline_North,xc_NW,xc_NE,2);
    Create_Spline_Line(Bnd_Spline_South,xc_SW,xc_SE,2);
    Create_Spline_Line(Bnd_Spline_East,xc_SE,xc_NE,2);
    Create_Spline_Line(Bnd_Spline_West,xc_SW,xc_NW,2);

    // Set the boundary condition types for each of the boundary splines.
    if (iBlk == 0) {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_FIXED);
    } else if (iBlk == 1) {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      if (Flat_Plate_BC_Type != BC_BURNING_SURFACE) {
	Bnd_Spline_South.setBCtype(Flat_Plate_BC_Type);
	Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      } else {
	Bnd_Spline_South.setBCtype(BC_WALL_VISCOUS_ISOTHERMAL);
	Bnd_Spline_East.setBCtype(BC_NONE);
      }
      Bnd_Spline_West.setBCtype(BC_NONE);
    } else {
      Bnd_Spline_North.setBCtype(BC_FIXED);
      Bnd_Spline_South.setBCtype(Flat_Plate_BC_Type);
      Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_West.setBCtype(BC_NONE);
    }

    // Determine the number of cells for this block.
    n_cells_i = Number_of_Cells_Idir;///2;
    n_cells_j = Number_of_Cells_Jdir;

    // Assign values to the stretching function parameters and boundary
    // grid line orthogonality parameters.
    Stretch_I = STRETCHING_FCN_LINEAR;
    Beta_I = ZERO;
    Tau_I = ZERO;
    Stretch_J = STRETCHING_FCN_LINEAR;
    Beta_J = ZERO;
    Tau_J = ZERO;
    if (Stretching_Flag) {
      if (iBlk == 0) {
	Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
      } else if (iBlk == 1) {
	if (Flat_Plate_BC_Type != BC_BURNING_SURFACE) {
	  Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	} else {
	  Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
	}
      } else {
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
      }
      Beta_I = Stretching_Factor_Idir;
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = Stretching_Factor_Jdir;
    }
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

    // Deallocate the memory for the boundary splines.
    Bnd_Spline_North.deallocate();
    Bnd_Spline_South.deallocate();
    Bnd_Spline_East.deallocate();
    Bnd_Spline_West.deallocate();

  }

  // Return the grid.
  return Grid_ptr;

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
			      const int Stretching_Flag,
			      const double Stretching_Factor,
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
   
  Bnd_Spline_North.setBCtype(BC_WALL_VISCOUS_ISOTHERMAL);
  Bnd_Spline_South.setBCtype(BC_REFLECTION);
  Bnd_Spline_East.setBCtype(BC_FIXED);
  Bnd_Spline_West.setBCtype(BC_FIXED);

  /* Assign values to the stretching function parameters
     and boundary grid line orthogonality parameters. */
  Stretch_I = STRETCHING_FCN_LINEAR;
  Beta_I = ZERO; 
  Tau_I = ZERO;
  Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
  Beta_J = ZERO;
  Tau_J = ZERO; 
  if (Stretching_Flag) Beta_J = Stretching_Factor;
  else Stretch_J = STRETCHING_FCN_LINEAR;

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

    if (Axisymmetric == 1) {
      Bnd_Spline_North.setBCtype(BC_WALL_VISCOUS_ISOTHERMAL);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_FIXED);
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

    if (Axisymmetric == 2) {
      Bnd_Spline_East.setBCtype(BC_WALL_VISCOUS_ISOTHERMAL);
      Bnd_Spline_West.setBCtype(BC_REFLECTION);
      Bnd_Spline_North.setBCtype(BC_FIXED);
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
				      const double &Length_Chamber,
				      const double &Radius_Chamber,
				      const double &Length_Chamber_To_Throat,
				      const double &Length_Nozzle,
				      const double &Radius_Nozzle_Exit,
				      const double &Radius_Nozzle_Throat,
				      const double &Radius_Grain,
				      const int &Nozzle_Type,
				      const int &Chamber_BC_Type,
				      const int &Stretching_Flag,
				      const int Stretching_Type_Jdir,
				      const double &Stretching_Factor_Idir,
				      const double &Stretching_Factor_Jdir,
				      const int Number_of_Cells_Idir,
				      const int Number_of_Cells_Jdir,
				      const int Number_of_Ghost_Cells) {

  int error_flag, block_flag;
  int Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D Bnd_Spline_North, Bnd_Spline_South,
           Bnd_Spline_East, Bnd_Spline_West;
  double num1 = ZERO, num2 = ZERO;
  int iBlk_Chamber, iBlk_Nozzle;

  // Allocate memory for grid blocks.
  if (Radius_Grain > Radius_Chamber) return NULL;
  if (Radius_Grain > ZERO) {
    iBlk_Nozzle = 2;
    Number_of_Blocks_Jdir = 2;
  } else if (Radius_Grain < ZERO) {
    iBlk_Nozzle = 1;
    Number_of_Blocks_Jdir = 2;
  } else {
    iBlk_Nozzle = 1;
    Number_of_Blocks_Jdir = 1;
  }
  iBlk_Chamber = Number_of_Blocks_Idir-1;
  Number_of_Blocks_Idir = iBlk_Chamber + iBlk_Nozzle;
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr,
				       Number_of_Blocks_Idir,
				       Number_of_Blocks_Jdir);

  // Create the mesh for each block representing the complete grid.
  for (int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++) {
    for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

      // Default the current block to off (a hole).
      block_flag = OFF;

      // Default the stretching parameters.
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ONE;
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;
      Beta_J = ONE;
      Tau_J = ZERO;
      Orthogonal_North = 1;
      Orthogonal_South = 1;
      Orthogonal_East = 1;
      Orthogonal_West = 1;

//       if (iBlk < iBlk_Chamber) {
// 	cout << endl << StretchingFcn(double(iBlk*Number_of_Cells_Idir)/double(iBlk_Chamber*Number_of_Cells_Idir),Stretching_Factor_Idir,ZERO,STRETCHING_FCN_MIN_CLUSTERING);
// 	cout << " " << StretchingFcn(double((iBlk+1)*Number_of_Cells_Idir)/double(iBlk_Chamber*Number_of_Cells_Idir),Stretching_Factor_Idir,ZERO,STRETCHING_FCN_MIN_CLUSTERING);
//       }

      // Create the splines defining the north, south, east, and west 
      // boundaries of the grid block.
      if (jBlk == 0 && iBlk < iBlk_Chamber) {
	num1 = Length_Chamber*double(iBlk)/double(Number_of_Blocks_Idir-iBlk_Nozzle);
	num2 = Length_Chamber*double(iBlk+1)/double(Number_of_Blocks_Idir-iBlk_Nozzle);
//  	num1 = Length_Chamber*StretchingFcn(double(iBlk)/double(iBlk_Chamber),Stretching_Factor_Idir,ZERO,STRETCHING_FCN_MIN_CLUSTERING);
//  	num2 = Length_Chamber*StretchingFcn(double(iBlk+1)/double(iBlk_Chamber),Stretching_Factor_Idir,ZERO,STRETCHING_FCN_MIN_CLUSTERING);
	xc_SE = Vector2D(-Length_Chamber+num2,ZERO);
	xc_SW = Vector2D(-Length_Chamber+num1,ZERO);
	if (iBlk_Nozzle == 1) {
	  xc_NW = Vector2D(-Length_Chamber+num1,Radius_Chamber);
	  xc_NE = Vector2D(-Length_Chamber+num2,Radius_Chamber);
	} else {
	  xc_NW = Vector2D(-Length_Chamber+num1,Radius_Grain);
	  xc_NE = Vector2D(-Length_Chamber+num2,Radius_Grain);
	}
	Create_Spline_Line(Bnd_Spline_North,xc_NW,xc_NE,2);
	Create_Spline_Line(Bnd_Spline_South,xc_SW,xc_SE,2);
	Create_Spline_Line(Bnd_Spline_East,xc_SE,xc_NE,2);
	Create_Spline_Line(Bnd_Spline_West,xc_SW,xc_NW,2);
	if (iBlk_Nozzle == 1 && Radius_Grain < ZERO) Bnd_Spline_North.setBCtype(BC_NONE); 
	else Bnd_Spline_North.setBCtype(BC_BURNING_SURFACE);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	Bnd_Spline_East.setBCtype(BC_NONE);
	if (iBlk == 0) Bnd_Spline_West.setBCtype(Chamber_BC_Type);
	else Bnd_Spline_West.setBCtype(BC_NONE);
	if (Stretching_Flag) {
	  if (iBlk == 0) {
	    Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	    Beta_I = Stretching_Factor_Idir;
	    Tau_I = ZERO;
	  } else if (iBlk == iBlk_Chamber-1) {
	    Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	    Stretch_I = STRETCHING_FCN_LINEAR;
	    Beta_I = Stretching_Factor_Idir;
	    Tau_I = ZERO;
	  }
	}
	Stretch_J = Stretching_Type_Jdir;
	Beta_J = Stretching_Factor_Jdir;
	Tau_J = ZERO;
	block_flag = ON;

      } else if (jBlk == 1 && iBlk < iBlk_Chamber && Radius_Grain < ZERO) {
	num1 = Length_Chamber*double(iBlk)/double(Number_of_Blocks_Idir-iBlk_Nozzle);
	num2 = Length_Chamber*double(iBlk+1)/double(Number_of_Blocks_Idir-iBlk_Nozzle);
//  	num1 = Length_Chamber*StretchingFcn(double(iBlk)/double(iBlk_Chamber),Stretching_Factor_Idir,ZERO,STRETCHING_FCN_MIN_CLUSTERING);
//  	num2 = Length_Chamber*StretchingFcn(double(iBlk+1)/double(iBlk_Chamber),Stretching_Factor_Idir,ZERO,STRETCHING_FCN_MIN_CLUSTERING);
	xc_NW = Vector2D(-Length_Chamber+num1,Radius_Chamber-Radius_Grain);
	xc_NE = Vector2D(-Length_Chamber+num2,Radius_Chamber-Radius_Grain);
	xc_SE = Vector2D(-Length_Chamber+num2,Radius_Chamber);
	xc_SW = Vector2D(-Length_Chamber+num1,Radius_Chamber);
	Create_Spline_Line(Bnd_Spline_North,xc_NW,xc_NE,2);
	Create_Spline_Line(Bnd_Spline_South,xc_SW,xc_SE,2);
	Create_Spline_Line(Bnd_Spline_East,xc_SE,xc_NE,2);
	Create_Spline_Line(Bnd_Spline_West,xc_SW,xc_NW,2);
	Bnd_Spline_North.setBCtype(BC_BURNING_SURFACE);
	Bnd_Spline_South.setBCtype(BC_NONE);
	Bnd_Spline_East.setBCtype(BC_NONE);
	if (iBlk == iBlk_Chamber-1) Bnd_Spline_East.setBCtype(Chamber_BC_Type);
	else Bnd_Spline_East.setBCtype(BC_NONE);
	if (iBlk == 0) Bnd_Spline_West.setBCtype(Chamber_BC_Type);
	else Bnd_Spline_West.setBCtype(BC_NONE);
	if (Stretching_Flag) {
	  if (iBlk == 0) {
	    Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	    Beta_I = Stretching_Factor_Idir;
	    Tau_I = ZERO;
	  } else if (iBlk == iBlk_Chamber-1) {
	    Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	    Stretch_I = STRETCHING_FCN_LINEAR;
	    Beta_I = Stretching_Factor_Idir;
	    Tau_I = ZERO;
	  }
	}
	//Stretch_J = Stretching_Type_Jdir;
	//Beta_J = Stretching_Factor_Jdir;
	//Tau_J = ZERO;
	block_flag = ON;

      } else if (jBlk == 0 && iBlk >= iBlk_Chamber && iBlk_Nozzle == 1) {
	Create_Spline_Area_Variation(Bnd_Spline_North,
				     ZERO,
				     Length_Chamber_To_Throat,
				     Length_Nozzle,
				     Radius_Chamber,
				     Radius_Nozzle_Throat,
				     Radius_Nozzle_Exit,
				     Nozzle_Type,
				     501);
	xc_NW = Vector2D(ZERO,Radius_Chamber);
	xc_NE = Bnd_Spline_North.Xp[500];
	xc_SE = Vector2D(Length_Nozzle,ZERO);
	xc_SW = Vector2D(ZERO,ZERO);
	Create_Spline_Line(Bnd_Spline_East,xc_SE,xc_NE,2);
	Create_Spline_Line(Bnd_Spline_West,xc_SW,xc_NW,2);
	Create_Spline_Line(Bnd_Spline_South,xc_SW,xc_SE,2);
	Bnd_Spline_North.setBCtype(Chamber_BC_Type);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
	Bnd_Spline_East.setBCtype(BC_CHARACTERISTIC);
	Bnd_Spline_West.setBCtype(BC_NONE);
	Stretch_I = STRETCHING_FCN_MIDPT_CLUSTERING;
	Beta_I = Length_Chamber_To_Throat/Length_Nozzle;
	Tau_I = 3.00;
	Stretch_J = Stretching_Type_Jdir;
	Beta_J = Stretching_Factor_Jdir;
	Tau_J = ZERO;
	block_flag = ON;

      } else if (jBlk == 0 && iBlk - iBlk_Chamber == 0 && iBlk_Nozzle == 2) {
	xc_SW = Vector2D(ZERO,ZERO);
	xc_SE = Vector2D(Length_Chamber_To_Throat,ZERO);
	xc_NE = Vector2D(Length_Chamber_To_Throat,HALF*Radius_Nozzle_Throat);
	xc_NW = Vector2D(ZERO,Radius_Grain);
	Create_Spline_Line(Bnd_Spline_East,xc_SE,xc_NE,2);
	Create_Spline_Line(Bnd_Spline_West,xc_SW,xc_NW,2);
	Create_Spline_Line(Bnd_Spline_South,xc_SW,xc_SE,2);
	Create_Spline_Converging_Nozzle(Bnd_Spline_North,
					ZERO,
					Length_Chamber_To_Throat,
					Radius_Grain,
					HALF*Radius_Nozzle_Throat,
					251);
	Bnd_Spline_North.setBCtype(BC_NONE);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	Bnd_Spline_East.setBCtype(BC_NONE);
	Bnd_Spline_West.setBCtype(BC_NONE);
	if (Stretching_Flag) {
	  Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
	  Beta_I = Stretching_Factor_Idir;
	  Tau_I = ZERO;
	  Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	  Beta_J = Stretching_Factor_Jdir;
	  Tau_J = ZERO;
	}
	block_flag = ON;

      } else if (jBlk == 0 && iBlk - iBlk_Chamber == 1 && iBlk_Nozzle == 2) {
	xc_SW = Vector2D(Length_Chamber_To_Throat,ZERO);
	xc_SE = Vector2D(Length_Nozzle,ZERO);
	xc_NE = Vector2D(Length_Nozzle,HALF*Radius_Nozzle_Exit);
	xc_NW = Vector2D(Length_Chamber_To_Throat,HALF*Radius_Nozzle_Throat);
	Create_Spline_Line(Bnd_Spline_East,xc_SE,xc_NE,2);
	Create_Spline_Line(Bnd_Spline_West,xc_SW,xc_NW,2);
	Create_Spline_Line(Bnd_Spline_South,xc_SW,xc_SE,2);
	Create_Spline_Diverging_Nozzle(Bnd_Spline_North,
				       Length_Chamber_To_Throat,
				       Length_Nozzle,
				       HALF*Radius_Nozzle_Throat,
				       HALF*Radius_Nozzle_Exit,
				       251);
	Bnd_Spline_North.setBCtype(BC_NONE);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
	Bnd_Spline_East.setBCtype(BC_CHARACTERISTIC);
	Bnd_Spline_West.setBCtype(BC_NONE);
	if (Stretching_Flag) {
	  Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	  Beta_I = Stretching_Factor_Idir;
	  Tau_I = ZERO;
	  Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	  Beta_J = Stretching_Factor_Jdir;
	  Tau_J = ZERO;
	}
	block_flag = ON;

      } else if (jBlk == 1 && iBlk - iBlk_Chamber == 0 && iBlk_Nozzle == 2) {
	xc_SW = Vector2D(ZERO,Radius_Grain);
	xc_SE = Vector2D(Length_Chamber_To_Throat,HALF*Radius_Nozzle_Throat);
	xc_NE = Vector2D(Length_Chamber_To_Throat,Radius_Nozzle_Throat);
	xc_NW = Vector2D(ZERO,Radius_Chamber);
	Create_Spline_Line(Bnd_Spline_East,xc_SE,xc_NE,2);
	Create_Spline_Line(Bnd_Spline_West,xc_SW,xc_NW,2);
	Create_Spline_Converging_Nozzle(Bnd_Spline_South,
					ZERO,
					Length_Chamber_To_Throat,
					Radius_Grain,
					HALF*Radius_Nozzle_Throat,
					251);
	Create_Spline_Converging_Nozzle(Bnd_Spline_North,
					ZERO,
					Length_Chamber_To_Throat,
					Radius_Chamber,
					Radius_Nozzle_Throat,
					251);
	Bnd_Spline_North.setBCtype(Chamber_BC_Type);
	Bnd_Spline_South.setBCtype(BC_NONE);
	Bnd_Spline_East.setBCtype(BC_NONE);
	Bnd_Spline_West.setBCtype(BC_NONE);
	if (Stretching_Flag) {
	  Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
	  Beta_I = Stretching_Factor_Idir;
	  Tau_I = ZERO;
	  Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	  Beta_J = Stretching_Factor_Jdir;
	  Tau_J = ZERO;
	}
	block_flag = ON;

      } else if (jBlk == 1 && iBlk - iBlk_Chamber == 1 && iBlk_Nozzle == 2) {
	xc_SW = Vector2D(Length_Chamber_To_Throat,HALF*Radius_Nozzle_Throat);
	xc_SE = Vector2D(Length_Nozzle,HALF*Radius_Nozzle_Exit);
	xc_NE = Vector2D(Length_Nozzle,Radius_Nozzle_Exit);
	xc_NW = Vector2D(Length_Chamber_To_Throat,Radius_Nozzle_Throat);
	Create_Spline_Line(Bnd_Spline_East,xc_SE,xc_NE,2);
	Create_Spline_Line(Bnd_Spline_West,xc_SW,xc_NW,2);
	Create_Spline_Diverging_Nozzle(Bnd_Spline_South,
				       Length_Chamber_To_Throat,
				       Length_Nozzle,
				       HALF*Radius_Nozzle_Throat,
				       HALF*Radius_Nozzle_Exit,
				       251);
	Create_Spline_Diverging_Nozzle(Bnd_Spline_North,
				       Length_Chamber_To_Throat,
				       Length_Nozzle,
				       Radius_Nozzle_Throat,
				       Radius_Nozzle_Exit,
				       251);
	Bnd_Spline_North.setBCtype(Chamber_BC_Type);
	Bnd_Spline_South.setBCtype(BC_NONE);
	Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
	Bnd_Spline_East.setBCtype(BC_CHARACTERISTIC);
	Bnd_Spline_West.setBCtype(BC_NONE);
	if (Stretching_Flag) {
	  Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	  Beta_I = Stretching_Factor_Idir;
	  Tau_I = ZERO;
	  Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	  Beta_J = Stretching_Factor_Jdir;
	  Tau_J = ZERO;
	}
	block_flag = ON;

      } else {
	block_flag = OFF;

      }

      if (block_flag) {
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
// 	Smooth_Quad_Block(Grid_ptr[iBlk][jBlk],
// 			  min(250,2*max(Number_of_Cells_Idir,
// 					Number_of_Cells_Jdir)));
	Smooth_Rocket_Motor(Grid_ptr[iBlk][jBlk],
			    Length_Chamber,
			    Radius_Chamber,
			    Length_Chamber_To_Throat,
			    Length_Nozzle,
			    Radius_Nozzle_Exit,
			    Radius_Nozzle_Throat,
			    Radius_Grain,
			    Nozzle_Type,
			    Stretching_Factor_Idir,
			    Stretching_Factor_Jdir,
			    -1,
			    0,
			    0,0,
			    iBlk,jBlk,
			    Number_of_Blocks_Idir,
			    Number_of_Blocks_Jdir);

	// Deallocate the memory for the boundary splines.
	Bnd_Spline_North.deallocate();
	Bnd_Spline_South.deallocate();
	Bnd_Spline_East.deallocate();
	Bnd_Spline_West.deallocate();

      }

    }
  }

  // Return the grid.
  return Grid_ptr;

}

/**********************************************************************
 * Routine: Grid_Nozzleless_Rocket_Motor                              *
 *                                                                    *
 **********************************************************************/
Grid2D_Quad_Block** Grid_Nozzleless_Rocket_Motor(Grid2D_Quad_Block **Grid_ptr,
						 int &Number_of_Blocks_Idir,
						 int &Number_of_Blocks_Jdir,
						 const double &Length_Chamber,
						 const double &Radius_Chamber,
						 const double &Length_Nozzle,
						 const double &Radius_Nozzle_Exit,
						 const int &Chamber_BC_Type,
						 const int &Stretching_Flag,
						 const int Stretching_Type_Jdir,
						 const double &Stretching_Factor_Idir,
						 const double &Stretching_Factor_Jdir,
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

  // Allocate memory for grid blocks.
  Number_of_Blocks_Jdir = 1;
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr,
				       Number_of_Blocks_Idir,
				       Number_of_Blocks_Jdir);

  // Create the mesh for each block representing the complete grid.
  for (int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++) {
    for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

      // Default the stretching parameters.
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ONE;
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;
      Beta_J = ONE;
      Tau_J = ZERO;
      Orthogonal_North = 1;
      Orthogonal_South = 1;
      Orthogonal_East = 1;
      Orthogonal_West = 1;

      // Create the splines defining the north, south, east, and west 
      // boundaries of the grid block.
      if (iBlk < Number_of_Blocks_Idir-1) {//2) {
// 	xc_SW = Vector2D((Length_Chamber-0.002)*double(iBlk)/double(Number_of_Blocks_Idir-2),ZERO);
// 	xc_SE = Vector2D((Length_Chamber-0.002)*double(iBlk+1)/double(Number_of_Blocks_Idir-2),ZERO);
// 	xc_NW = Vector2D((Length_Chamber-0.002)*double(iBlk)/double(Number_of_Blocks_Idir-2),Radius_Chamber);
// 	xc_NE = Vector2D((Length_Chamber-0.002)*double(iBlk+1)/double(Number_of_Blocks_Idir-2),Radius_Chamber);
	xc_SW = Vector2D(Length_Chamber*double(iBlk)/double(Number_of_Blocks_Idir-1),ZERO);
	xc_SE = Vector2D(Length_Chamber*double(iBlk+1)/double(Number_of_Blocks_Idir-1),ZERO);
	xc_NW = Vector2D(Length_Chamber*double(iBlk)/double(Number_of_Blocks_Idir-1),Radius_Chamber);
	xc_NE = Vector2D(Length_Chamber*double(iBlk+1)/double(Number_of_Blocks_Idir-1),Radius_Chamber);
	Create_Spline_Line(Bnd_Spline_North,xc_NW,xc_NE,2);
	Create_Spline_Line(Bnd_Spline_South,xc_SW,xc_SE,2);
	Create_Spline_Line(Bnd_Spline_East,xc_SE,xc_NE,2);
	Create_Spline_Line(Bnd_Spline_West,xc_SW,xc_NW,2);
	Bnd_Spline_North.setBCtype(BC_MASS_INJECTION);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	Bnd_Spline_East.setBCtype(BC_NONE);
	if (iBlk == 0) Bnd_Spline_West.setBCtype(Chamber_BC_Type);
	else Bnd_Spline_West.setBCtype(BC_NONE);
//       } else if (iBlk == Number_of_Blocks_Idir-2) {
// 	xc_SW = Vector2D(Length_Chamber-0.002,ZERO);
// 	xc_SE = Vector2D(Length_Chamber,ZERO);
// 	xc_NW = Vector2D(Length_Chamber-0.002,Radius_Chamber);
// 	xc_NE = Vector2D(Length_Chamber,Radius_Chamber);
// 	Create_Spline_Line(Bnd_Spline_North,xc_NW,xc_NE,2);
// 	Create_Spline_Line(Bnd_Spline_South,xc_SW,xc_SE,2);
// 	Create_Spline_Line(Bnd_Spline_East,xc_SE,xc_NE,2);
// 	Create_Spline_Line(Bnd_Spline_West,xc_SW,xc_NW,2);
// 	Bnd_Spline_North.setBCtype(Chamber_BC_Type);
// 	Bnd_Spline_South.setBCtype(BC_REFLECTION);
// 	Bnd_Spline_East.setBCtype(BC_NONE);
// 	if (iBlk == 0) Bnd_Spline_West.setBCtype(Chamber_BC_Type);
// 	else Bnd_Spline_West.setBCtype(BC_NONE);
      } else {
	xc_SW = Vector2D(Length_Chamber,ZERO);
	xc_SE = Vector2D(Length_Chamber+Length_Nozzle,ZERO);
	xc_NW = Vector2D(Length_Chamber,Radius_Chamber);
	xc_NE = Vector2D(Length_Chamber+Length_Nozzle,Radius_Nozzle_Exit);
	Create_Spline_Line(Bnd_Spline_North,xc_NW,xc_NE,2);
	Create_Spline_Line(Bnd_Spline_South,xc_SW,xc_SE,2);
	Create_Spline_Line(Bnd_Spline_East,xc_SE,xc_NE,2);
	Create_Spline_Line(Bnd_Spline_West,xc_SW,xc_NW,2);
	Bnd_Spline_North.setBCtype(Chamber_BC_Type);
	Bnd_Spline_South.setBCtype(BC_REFLECTION);
	//Bnd_Spline_East.setBCtype(BC_CHARACTERISTIC);
	//Bnd_Spline_East.setBCtype(BC_FIXED_PRESSURE);
	Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
	Bnd_Spline_West.setBCtype(BC_NONE);
      }
      if (Stretching_Flag) {
	Stretch_J = Stretching_Type_Jdir;
	Beta_J = Stretching_Factor_Jdir;
	Tau_J = ZERO;
	if (iBlk == 0) {
	  Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	  Beta_I = Stretching_Factor_Idir;
	} else if (iBlk == Number_of_Blocks_Idir-2) {//3) {
	  Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	  Beta_I = Stretching_Factor_Idir;
	} else if (iBlk == Number_of_Blocks_Idir-1) {
	  Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	  Beta_I = Stretching_Factor_Idir;
	}
      }

      // Create the quadrilateral solution block.
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

/**********************************************************************
 * Routine: Grid_Nozzle                                               *
 *                                                                    *
 * Generates a two-block grid for a rocket motor nozzle.              *
 *                                                                    *
 **********************************************************************/
Grid2D_Quad_Block** Grid_Nozzle(Grid2D_Quad_Block **Grid_ptr,
				int &Number_of_Blocks_Idir,
				int &Number_of_Blocks_Jdir,
				const double &Length_Nozzle,
				const double &Radius_Chamber,
				const double &Radius_Nozzle_Exit,
				const double &Radius_Nozzle_Throat,
				const int &Nozzle_Type,
				const int &Stretching_Flag,
				const int &Stretching_Type_Idir,
				const int &Stretching_Type_Jdir,
				const double &Stretching_Factor_Idir,
				const double &Stretching_Factor_Jdir,
				const int Number_of_Cells_Idir,
				const int Number_of_Cells_Jdir,
				const int Number_of_Ghost_Cells) {

  int error_flag, block_flag;
  int Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xc_NW, xc_NE, xc_SE, xc_SW;
  Spline2D Bnd_Spline_North, Bnd_Spline_South,
           Bnd_Spline_East, Bnd_Spline_West;

  // Allocate memory for grid blocks.
  if (!Nozzle_Type) Number_of_Blocks_Idir = 1;
  else Number_of_Blocks_Idir = 2;
  Number_of_Blocks_Jdir = 1;
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr,
				       Number_of_Blocks_Idir,
				       Number_of_Blocks_Jdir);

  // Create the mesh for each block representing the complete grid.
  for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

    if (iBlk == 0) {
      // Converging section of the nozzle.
      xc_SW = Vector2D(ZERO,ZERO);
      xc_SE = Vector2D(HALF*Length_Nozzle,ZERO);
      xc_NE = Vector2D(HALF*Length_Nozzle,Radius_Nozzle_Throat);
      xc_NW = Vector2D(ZERO,Radius_Chamber);
      Create_Spline_Converging_Nozzle(Bnd_Spline_North,
				      ZERO,
				      HALF*Length_Nozzle,
				      Radius_Chamber,
				      Radius_Nozzle_Throat,
				      251);
      Create_Spline_Line(Bnd_Spline_South,xc_SW,xc_SE,2);
      Create_Spline_Line(Bnd_Spline_East,xc_SE,xc_NE,2);
      Create_Spline_Line(Bnd_Spline_West,xc_SW,xc_NW,2);
      Bnd_Spline_North.setBCtype(BC_REFLECTION);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_FIXED);
    } else {
      // Diverging section of the nozzle.
      xc_SW = Vector2D(HALF*Length_Nozzle,ZERO);
      xc_SE = Vector2D(Length_Nozzle,ZERO);
      xc_NE = Vector2D(Length_Nozzle,Radius_Nozzle_Exit);
      xc_NW = Vector2D(HALF*Length_Nozzle,Radius_Nozzle_Throat);
      Create_Spline_Diverging_Nozzle(Bnd_Spline_North,
				     HALF*Length_Nozzle,
				     Length_Nozzle,
				     Radius_Nozzle_Throat,
				     Radius_Nozzle_Exit,
				     251);
      Create_Spline_Line(Bnd_Spline_South,xc_SW,xc_SE,2);
      Create_Spline_Line(Bnd_Spline_East,xc_SE,xc_NE,2);
      Create_Spline_Line(Bnd_Spline_West,xc_SW,xc_NW,2);
      Bnd_Spline_North.setBCtype(BC_REFLECTION);
      Bnd_Spline_South.setBCtype(BC_REFLECTION);
      Bnd_Spline_East.setBCtype(BC_FIXED);
      Bnd_Spline_West.setBCtype(BC_NONE);
    }

    // Assign values to the stretching function parameters and 
    // boundary grid line orthogonality parameters.
    Stretch_I = STRETCHING_FCN_LINEAR; Beta_I = ONE; Tau_I = ZERO;
    Stretch_J = STRETCHING_FCN_LINEAR; Beta_J = ONE; Tau_J = ZERO;
    if (Stretching_Flag) {
      if (Stretching_Type_Idir != STRETCHING_FCN_LINEAR) {
	if (iBlk == 0) Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	else Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_I = Stretching_Factor_Idir;
      }
      Stretch_J = Stretching_Type_Jdir;
      Beta_J = Stretching_Factor_Jdir;
    }
    Orthogonal_North = 1;
    Orthogonal_South = 1;
    Orthogonal_East = 1;
    Orthogonal_West = 1;

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

    // Smooth the 2D quadrilateral grid block.
    Smooth_Quad_Block(Grid_ptr[iBlk][0],
		      min(250,2*max(Number_of_Cells_Idir,
				    Number_of_Cells_Jdir)));

    // Deallocate the memory for the boundary splines.
    Bnd_Spline_North.deallocate();
    Bnd_Spline_South.deallocate();
    Bnd_Spline_East.deallocate();
    Bnd_Spline_West.deallocate();

  }

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
					   const int &Stretching_Type_Idir,
					   const int &Stretching_Type_Jdir,
					   const double &Stretching_Factor_Idir,
					   const double &Stretching_Factor_Jdir,
 		                           const int Number_of_Cells_Idir,
		                           const int Number_of_Cells_Jdir,
					   const int Number_of_Ghost_Cells) {

  return Grid_Circular_Cylinder(Grid_ptr,
				Number_of_Blocks_Idir,
                                Number_of_Blocks_Jdir,
                                Radius,
				32.00*Radius,
				Stretching_Type_Idir,
				Stretching_Type_Jdir,
				Stretching_Factor_Idir,
				Stretching_Factor_Jdir,
 		                Number_of_Cells_Idir,
		                Number_of_Cells_Jdir,
				Number_of_Ghost_Cells);

}

Grid2D_Quad_Block** Grid_Circular_Cylinder(Grid2D_Quad_Block **Grid_ptr,
                                           int &Number_of_Blocks_Idir,
                                           int &Number_of_Blocks_Jdir,
                                           const double &Inner_Radius,
                                           const double &Outer_Radius,
					   const int &Stretching_Type_Idir,
					   const int &Stretching_Type_Jdir,
					   const double &Stretching_Factor_Idir,
					   const double &Stretching_Factor_Jdir,
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
           x1 = Vector2D(ZERO,ZERO);
           Create_Spline_Circular_Arc(Bnd_Spline_North,
	      		              x1,
			              Outer_Radius,
                                      360.00,
			              180.00,
  	                              361);
           Create_Spline_Circular_Arc(Bnd_Spline_South,
			              x1,
			              Inner_Radius,
                                      360.00,
			              180.00,
  	                              361);
           x1 = Vector2D(Inner_Radius, ZERO);
           x2 = Vector2D(Outer_Radius, ZERO);
           Create_Spline_Line(Bnd_Spline_West, x1, x2, 2);
           x1 = Vector2D(-Inner_Radius, ZERO);
           x2 = Vector2D(-Outer_Radius, ZERO);
           Create_Spline_Line(Bnd_Spline_East, x1, x2, 2);
        } else {
           x1 = Vector2D(ZERO,ZERO);
           Create_Spline_Circular_Arc(Bnd_Spline_North,
	      		              x1,
			              Outer_Radius,
                                      180.00,
			              ZERO,
  	                              361);
           Create_Spline_Circular_Arc(Bnd_Spline_South,
			              x1,
			              Inner_Radius,
                                      180.00,
			              ZERO,
  	                              361);
           x1 = Vector2D(-Inner_Radius, ZERO);
           x2 = Vector2D(-Outer_Radius, ZERO);
           Create_Spline_Line(Bnd_Spline_West, x1, x2, 2);
           x1 = Vector2D(Inner_Radius, ZERO);
           x2 = Vector2D(Outer_Radius, ZERO);
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
	   Stretch_I = Stretching_Type_Idir;
           Beta_I = Stretching_Factor_Idir;
           Tau_I = ZERO;
           Stretch_J = Stretching_Type_Jdir;
           Beta_J = Stretching_Factor_Jdir;
           Tau_J = ZERO;
        } else {
           Stretch_I = Stretching_Type_Idir;
           Beta_I = Stretching_Factor_Idir;
           Tau_I = ZERO;
           Stretch_J = Stretching_Type_Jdir;
           Beta_J = Stretching_Factor_Jdir;
           Tau_J = ZERO;
        } /* endif */

	Orthogonal_North = 0;
	Orthogonal_South = 0;
	Orthogonal_East = 0;
	Orthogonal_West = 0;

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
 *  		                  50,                   *
 *  		                  2);                   *
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
           Bnd_Spline_South.setBCtype(BC_WALL_VISCOUS_HEATFLUX);//BC_REFLECTION);
           Bnd_Spline_East.setBCtype(BC_NONE);
           Bnd_Spline_West.setBCtype(BC_NONE);
        } else if (iBlk == 2) {
           Bnd_Spline_North.setBCtype(BC_FIXED);
           Bnd_Spline_South.setBCtype(BC_WALL_VISCOUS_HEATFLUX);//BC_REFLECTION);
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
 * Routine: Grid_Mixing_Layer                                         *
 *                                                                    *
 * Generates a quadilateral mesh with for predicting mixing layer     *
 * flows.                                                             *
 *                                                                    *
 * Usage: Grid_ptr = Grid_Mixing_Layer(Grid_ptr,                      *
 *                                     nblk_i,                        *
 *                                     nblk_j,                        *
 *                                     length,                        *
 *				       1.01,                          *  
 *				       1.01,                          *  
 *				       32,                            *  
 *				       24,                            *  
 *                                     2)                             *
 *                                                                    *
 **********************************************************************/
Grid2D_Quad_Block** Grid_Mixing_Layer(Grid2D_Quad_Block **Grid_ptr,
				      int &Number_of_Blocks_Idir,
				      int &Number_of_Blocks_Jdir,
				      const double &Length,
				      const double &Stretching_Factor_Idir,
				      const double &Stretching_Factor_Jdir,
				      const int Number_of_Cells_Idir,
				      const int Number_of_Cells_Jdir,
				      const int Number_of_Ghost_Cells) {

  int block_flag,
      BCtypeN, BCtypeS, BCtypeE, BCtypeW,
      Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xNW, xNE, xSE, xSW;
  Spline2D Bnd_Spline_North, Bnd_Spline_South,
           Bnd_Spline_East, Bnd_Spline_West;

  // Allocate memory for grid blocks. There are two grid blocks for
  // this mesh.
  Number_of_Blocks_Idir = 1;
  Number_of_Blocks_Jdir = 2;
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr,
                                       Number_of_Blocks_Idir,
                                       Number_of_Blocks_Jdir);

  // Create the mesh for each block representing the complete grid.
  for (int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++) {
    
    // Preset all boundary conditions to none and reset when required:
    BCtypeN = BC_NONE;
    BCtypeS = BC_NONE;
    BCtypeE = BC_FIXED_PRESSURE;
    BCtypeW = BC_FIXED;
 
    // Preset all stretching parameters and reset when required:
    Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
    Beta_I = Stretching_Factor_Idir;
    Tau_I = ZERO;
    Stretch_J = STRETCHING_FCN_LINEAR;
    Beta_J = Stretching_Factor_Jdir;
    Tau_J = ZERO;
    Orthogonal_North = 0;
    Orthogonal_South = 0;
    Orthogonal_East = 0;
    Orthogonal_West = 0;

    if (jBlk == 0) {
      xSW = Vector2D(ZERO,-0.25*Length);
      xSE = Vector2D(Length,-0.25*Length);
      xNW = Vector2D(ZERO,ZERO);
      xNE = Vector2D(Length,ZERO);
      BCtypeS = BC_CONSTANT_EXTRAPOLATION;
      Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	
    } else{
      xSW = Vector2D(ZERO,ZERO);
      xSE = Vector2D(Length,ZERO);
      xNW = Vector2D(ZERO,0.25*Length);
      xNE = Vector2D(Length,0.25*Length);
      BCtypeN = BC_CONSTANT_EXTRAPOLATION;
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = Stretching_Factor_Jdir; 
    }

    // Create the splines defining the north, south, east, and west
    // boundaries of the rectangular boxes.
    Create_Spline_Line(Bnd_Spline_North,xNW,xNE,2);
    Create_Spline_Line(Bnd_Spline_South,xSW,xSE,2);
    Create_Spline_Line(Bnd_Spline_East,xSE,xNE,2);
    Create_Spline_Line(Bnd_Spline_West,xSW,xNW,2);

    // Set the boundary condition types for each of the boundary
    // splines:
    Bnd_Spline_North.setBCtype(BCtypeN);
    Bnd_Spline_South.setBCtype(BCtypeS);
    Bnd_Spline_East.setBCtype(BCtypeE);
    Bnd_Spline_West.setBCtype(BCtypeW);

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

  }

  // Return the grid.
  return Grid_ptr;

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
			       const int &Wedge_BC_Type,
			       const int &Stretching_Flag,
			       const double &Stretching_Factor_Idir,
			       const double &Stretching_Factor_Jdir,
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
      xc_NE = Vector2D(ZERO,Wedge_Length);
      xc_SE = Vector2D(ZERO,ZERO);
      xc_SW = Vector2D(-Wedge_Length,ZERO);
    } else if (iBlk == 1) {
      xc_NW = Vector2D(ZERO,Wedge_Length);
      xc_NE = Vector2D(Wedge_Length*cos(PI*Wedge_Angle/180.0),Wedge_Length);
      xc_SE = Vector2D(Wedge_Length*cos(PI*Wedge_Angle/180.0),Wedge_Length*sin(PI*Wedge_Angle/180.0));
      xc_SW = Vector2D(ZERO,ZERO);
    }
    // Create the splines defining the north, south, east, and west 
    // boundaries of the rectangular box.
    Create_Spline_Line(Bnd_Spline_North, xc_NW, xc_NE, 2);
    Create_Spline_Line(Bnd_Spline_South, xc_SW, xc_SE, 2);
    Create_Spline_Line(Bnd_Spline_East, xc_SE, xc_NE, 2);
    Create_Spline_Line(Bnd_Spline_West, xc_SW, xc_NW, 2);
    // Set the boundary condition types for each of the boundary splines.
    Bnd_Spline_North.setBCtype(BC_REFLECTION);
    Bnd_Spline_South.setBCtype(Wedge_BC_Type);
    if (iBlk == 0) {
      Bnd_Spline_East.setBCtype(BC_NONE);
      Bnd_Spline_West.setBCtype(BC_FIXED);
    } else if (iBlk == 1) {
      Bnd_Spline_East.setBCtype(BC_CONSTANT_EXTRAPOLATION);
      Bnd_Spline_West.setBCtype(BC_NONE);
    }
    // Assign values to the stretching function parameters and
    // boundary grid line orthogonality parameters.
    if (Stretching_Flag) {
      if (iBlk == 0) Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
      else Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_I = Stretching_Factor_Idir;
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
      Beta_J = Stretching_Factor_Jdir;
    } else {
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ONE;
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;
      Beta_J = ONE;
      Tau_J = ZERO;
    }
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
    //Smooth_Quad_Block(Grid_ptr[iBlk][0], min(250, 2*max(Number_of_Cells_Idir,Number_of_Cells_Jdir)));
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
				      const double &Inner_Streamline_Number,
				      const double &Outer_Streamline_Number,
				      const double &Isotach_Line,
				      const int Number_of_Cells_Idir,
				      const int Number_of_Cells_Jdir,
				      const int Number_of_Ghost_Cells) {

  assert(Inner_Streamline_Number > Outer_Streamline_Number);
  assert(Inner_Streamline_Number < 5.0/3.0);
  assert(Outer_Streamline_Number > Isotach_Line); 

  int nk, nq,
      Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
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


  // Create the mesh for each block representing the complete grid.
  nk = 32;
  nq = 50;
  k = new double*[nk];
  q = new double*[nk];
  rho = new double*[nk];
  for (int i = 0; i < nk; i++) {
    k[i] = new double[nq];
    q[i] = new double[nq];
    rho[i] = new double[nq];
    for (int j = 0; j < nq; j++) {
      k[i][j] = ZERO;
      q[i][j] = ZERO;
      rho[i][j] = ZERO;
    }
  }

  delta_k = (Inner_Streamline_Number-Outer_Streamline_Number)/double(nk-1);
  k_init = Outer_Streamline_Number;
  q_init = Isotach_Line;
  for (int i = 0; i < nk; i++){
    ko = k_init + double(i)*delta_k;
    q_final = ko; // condition y = 0
    delta_q = (q_final - q_init)/double(nq-1);
    for (int j = 0; j < nq; j++) {
      if (j == nq-1) qo = q_final;
      else qo = q_init + double(j)*delta_q;
      k[i][j] = ko;
      q[i][j] = qo;
    }
  }

  Bnd_Spline_North.allocate(nk); Bnd_Spline_North.settype(SPLINE2D_QUINTIC);
  Bnd_Spline_South.allocate(nk); Bnd_Spline_South.settype(SPLINE2D_QUINTIC);
  Bnd_Spline_East.allocate(nq);	 Bnd_Spline_East.settype(SPLINE2D_QUINTIC);
  Bnd_Spline_West.allocate(nq);	 Bnd_Spline_West.settype(SPLINE2D_QUINTIC);

  for (int i = 0; i < nk; i++) {
    for (int j = 0; j < nq; j++){

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
	Bnd_Spline_East.bc[nq-1-j] = BC_RINGLEB_FLOW;//BC_REFLECTION;
	if (j == 0 || j == nq-1) Bnd_Spline_East.tp[nq-1-j] = SPLINE2D_POINT_SHARP_CORNER;
	else Bnd_Spline_East.tp[nq-1-j] = SPLINE2D_POINT_NORMAL;
      }
      // WEST spline.
      if (i == nk-1) {
	Bnd_Spline_West.Xp[nq-1-j].x = (HALF/rho[i][j])*(TWO/(k[i][j]*k[i][j]) - ONE/(q[i][j]*q[i][j])) - HALF*J;
	Bnd_Spline_West.Xp[nq-1-j].y = (ONE/(k[i][j]*rho[i][j]*q[i][j]))*sqrt(ONE - (q[i][j]*q[i][j])/(k[i][j]*k[i][j]));
	Bnd_Spline_West.bc[nq-1-j] = BC_RINGLEB_FLOW;//BC_REFLECTION;
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
  Bnd_Spline_East.setBCtype(BC_REFLECTION);
  Bnd_Spline_West.setBCtype(BC_REFLECTION);
  Bnd_Spline_East.setBCtype(BC_RINGLEB_FLOW);
  Bnd_Spline_West.setBCtype(BC_RINGLEB_FLOW);

  // Assign values to the stretching function parameters and
  // boundary grid line orthogonality parameters.
  Stretch_I = STRETCHING_FCN_LINEAR;
  Beta_I = ZERO; 
  Tau_I = ZERO;
  Stretch_J = STRETCHING_FCN_LINEAR;
  Beta_J = ONE;
  Tau_J = THREE;
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

  Smooth_Quad_Block(Grid_ptr[0][0],min(250,2*max(Number_of_Cells_Idir,Number_of_Cells_Jdir)));

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
					   const int Smooth_Bump,
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
      } else if (iBlk == 1 && jBlk == 0) {
	// Create the splines defining the north, south, east, and west 
	// boundaries of the grid.
	x1 = Vector2D( 0.0,1.0);
	x2 = Vector2D( 1.0,1.0);
	Create_Spline_Line(Bnd_Spline_North,x1,x2,2);
	if (!Smooth_Bump) {
	  // Non-smooth circular bump.
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
	} else {
	  // Smooth sin^2 bump.
	  Bnd_Spline_South.allocate(31);
	  Bnd_Spline_South.settype(SPLINE2D_QUINTIC);
	  for (int i = 0; i < 31; i++) {
	    Bnd_Spline_South.Xp[i].x = double(i)/30.0;
	    Bnd_Spline_South.Xp[i].y = 0.042*sqr(sin(PI*double(i)/30.0));
	    Bnd_Spline_South.bc[i] = BC_REFLECTION;
	    Bnd_Spline_South.tp[i] = SPLINE2D_POINT_NORMAL;
	  }
	  Bnd_Spline_South.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
	  Bnd_Spline_South.tp[30] = SPLINE2D_POINT_SHARP_CORNER;
	  Bnd_Spline_South.pathlength();
	}
	x1 = Vector2D( 1.0,0.0);
	x2 = Vector2D( 1.0,1.0);
	Create_Spline_Line(Bnd_Spline_East,x1,x2,2);
	x1 = Vector2D( 0.0,0.0);
	x2 = Vector2D( 0.0,1.0);
	Create_Spline_Line(Bnd_Spline_West,x1,x2,2);
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
      }

      // Set the boundary condition types for each of the boundary splines.
      if (jBlk == 0) Bnd_Spline_South.setBCtype(BC_REFLECTION);
      else Bnd_Spline_South.setBCtype(BC_NONE);
      if (jBlk == 1) Bnd_Spline_North.setBCtype(BC_REFLECTION);
      else Bnd_Spline_North.setBCtype(BC_NONE);
      if (iBlk == 0) Bnd_Spline_West.setBCtype(BC_FIXED);
      else Bnd_Spline_West.setBCtype(BC_NONE);
      if (iBlk == 3) Bnd_Spline_East.setBCtype(BC_FIXED);//BC_CHARACTERISTIC);//BC_CONSTANT_EXTRAPOLATION);
      else Bnd_Spline_East.setBCtype(BC_NONE);

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

/**********************************************************************
 * Routine: Grid_Backward_Facing_Step                                 *
 *                                                                    *
 * Generates a quadilateral mesh with clustering consisting of five   *
 * grid blocks for predicting viscous flow over a backward facing     *
 * step.                                                              *
 *                                                                    *
 * Usage: Grid_ptr = Grid_Backward_Facing_Step(Grid_ptr,              *
 *                                             nblk_i,                *
 *                                             nblk_j,                *
 *                                             ncells_i,              *
 *                                             ncells_j,              *
 *                                             1 or 0);               *
 *                                                                    *
 **********************************************************************/
Grid2D_Quad_Block** Grid_Backward_Facing_Step(Grid2D_Quad_Block **Grid_ptr,
					      int &Number_of_Blocks_Idir,
					      int &Number_of_Blocks_Jdir,
					      const double &Step_Height,
					      const double &Top_Wall_Deflection,
					      const double Stretching_Factor_Idir,
					      const double Stretching_Factor_Jdir,
					      const int Number_of_Cells_Idir,
					      const int Number_of_Cells_Jdir,
					      const int Number_of_Ghost_Cells) {

  int block_flag,
      BCtypeN, BCtypeS, BCtypeE, BCtypeW,
      Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xNW, xNE, xSE, xSW;
  Spline2D Bnd_Spline_North, Bnd_Spline_South,
           Bnd_Spline_East, Bnd_Spline_West;
  double Inlet_Height, Inlet_Length, Outlet_Length;

  // Set geometric constants.
  Inlet_Height = 8.0*Step_Height;
  Inlet_Length = 20.0*Step_Height;
  Outlet_Length = 50.0*Step_Height;

  // Allocate memory for grid blocks.  There are two grid blocks for 
  // this mesh.
  Number_of_Blocks_Idir = 7;
  Number_of_Blocks_Jdir = 3;
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr,
				       Number_of_Blocks_Idir,
				       Number_of_Blocks_Jdir);

  // Create the mesh for each block representing the complete grid.
  for (int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++) {
    for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

      // Preset the block indicator and turn off when required:
      block_flag = ON;

      // Preset all boundary conditions to none and reset when required:
      BCtypeN = BC_NONE;
      BCtypeS = BC_NONE;
      BCtypeE = BC_NONE;
      BCtypeW = BC_NONE;

      // Preset all stretching parameters and reset when required:
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ZERO;
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;
      Beta_J = ZERO;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;

      // Determine the block indicator and all boundary splines and
      // reset all boundary conditions and stretching parameters for the
      // current block: 
      if (iBlk == 0 && jBlk == 0) {
	block_flag = OFF;
      } else if (iBlk == 1 && jBlk == 0) {
	block_flag = OFF;
      } else if (iBlk == 2 && jBlk == 0) {
	xSW = Vector2D(ZERO,ZERO);
	xSE = Vector2D(0.2*Outlet_Length,ZERO);
	xNW = Vector2D(ZERO,Step_Height);
	xNE = Vector2D(0.2*Outlet_Length,Step_Height);
	BCtypeS = BC_WALL_VISCOUS_HEATFLUX;
	BCtypeW = BC_WALL_VISCOUS_HEATFLUX;
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.10;
	Stretch_J = STRETCHING_FCN_MINMAX_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 3 && jBlk == 0) {
	xSW = Vector2D(0.2*Outlet_Length,ZERO);
	xSE = Vector2D(0.4*Outlet_Length,ZERO);
	xNW = Vector2D(0.2*Outlet_Length,Step_Height);
	xNE = Vector2D(0.4*Outlet_Length,Step_Height);
	BCtypeS = BC_WALL_VISCOUS_HEATFLUX;
	Stretch_J = STRETCHING_FCN_MINMAX_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 4 && jBlk == 0) {
	xSW = Vector2D(0.4*Outlet_Length,ZERO);
	xSE = Vector2D(0.6*Outlet_Length,ZERO);
	xNW = Vector2D(0.4*Outlet_Length,Step_Height);
	xNE = Vector2D(0.6*Outlet_Length,Step_Height);
	BCtypeS = BC_WALL_VISCOUS_HEATFLUX;
	Stretch_J = STRETCHING_FCN_MINMAX_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 5 && jBlk == 0) {
	xSW = Vector2D(0.6*Outlet_Length,ZERO);
	xSE = Vector2D(0.8*Outlet_Length,ZERO);
	xNW = Vector2D(0.6*Outlet_Length,Step_Height);
	xNE = Vector2D(0.8*Outlet_Length,Step_Height);
	BCtypeS = BC_WALL_VISCOUS_HEATFLUX;
	Stretch_J = STRETCHING_FCN_MINMAX_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 6 && jBlk == 0) {
	xSW = Vector2D(0.8*Outlet_Length,ZERO);
	xSE = Vector2D(Outlet_Length,ZERO);
	xNW = Vector2D(0.8*Outlet_Length,Step_Height);
	xNE = Vector2D(Outlet_Length,Step_Height);
	BCtypeS = BC_WALL_VISCOUS_HEATFLUX;
	BCtypeE = BC_OUTFLOW_SUBSONIC;
	Stretch_J = STRETCHING_FCN_MINMAX_CLUSTERING;
	Beta_J = 1.10;

      } else if (iBlk == 0 && jBlk == 1) {
	xSW = Vector2D(-Inlet_Length,Step_Height);
	xSE = Vector2D(-HALF*Inlet_Length,Step_Height);
	xNW = Vector2D(-Inlet_Length,Step_Height+HALF*Inlet_Height);
	xNE = Vector2D(-HALF*Inlet_Length,Step_Height+HALF*Inlet_Height);
	BCtypeS = BC_WALL_VISCOUS_HEATFLUX;
	BCtypeW = BC_INFLOW_SUBSONIC;
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 1 && jBlk == 1) {
	xSW = Vector2D(-HALF*Inlet_Length,Step_Height);
	xSE = Vector2D(ZERO,Step_Height);
	xNW = Vector2D(-HALF*Inlet_Length,Step_Height+HALF*Inlet_Height);
	xNE = Vector2D(ZERO,Step_Height+HALF*Inlet_Height);
	BCtypeS = BC_WALL_VISCOUS_HEATFLUX;
	Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_I = 1.10;
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 2 && jBlk == 1) {
	xSW = Vector2D(ZERO,Step_Height);
	xSE = Vector2D(0.2*Outlet_Length,Step_Height);
	xNW = Vector2D(ZERO,Step_Height+HALF*Inlet_Height);
	xNE = Vector2D(0.2*Outlet_Length,Step_Height+HALF*Inlet_Height);
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.10;
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 3 && jBlk == 1) {
	xSW = Vector2D(0.2*Outlet_Length,Step_Height);
	xSE = Vector2D(0.4*Outlet_Length,Step_Height);
	xNW = Vector2D(0.2*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xNE = Vector2D(0.4*Outlet_Length,Step_Height+HALF*Inlet_Height);
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 4 && jBlk == 1) {
	xSW = Vector2D(0.4*Outlet_Length,Step_Height);
	xSE = Vector2D(0.6*Outlet_Length,Step_Height);
	xNW = Vector2D(0.4*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xNE = Vector2D(0.6*Outlet_Length,Step_Height+HALF*Inlet_Height);
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 5 && jBlk == 1) {
	xSW = Vector2D(0.6*Outlet_Length,Step_Height);
	xSE = Vector2D(0.8*Outlet_Length,Step_Height);
	xNW = Vector2D(0.6*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xNE = Vector2D(0.8*Outlet_Length,Step_Height+HALF*Inlet_Height);
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 6 && jBlk == 1) {
	xSW = Vector2D(0.8*Outlet_Length,Step_Height);
	xSE = Vector2D(Outlet_Length,Step_Height);
	xNW = Vector2D(0.8*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xNE = Vector2D(Outlet_Length,Step_Height+HALF*Inlet_Height);
	BCtypeE = BC_OUTFLOW_SUBSONIC;
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.10;

      } else if (iBlk == 0 && jBlk == 2) {
	xSW = Vector2D(-Inlet_Length,Step_Height+HALF*Inlet_Height);
	xSE = Vector2D(-HALF*Inlet_Length,Step_Height+HALF*Inlet_Height);
	xNW = Vector2D(-Inlet_Length,Step_Height+Inlet_Height);
	xNE = Vector2D(-HALF*Inlet_Length,Step_Height+Inlet_Height);
	BCtypeN = BC_WALL_VISCOUS_HEATFLUX;
	BCtypeW = BC_INFLOW_SUBSONIC;
	Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 1 && jBlk == 2) {
	xSW = Vector2D(-HALF*Inlet_Length,Step_Height+HALF*Inlet_Height);
	xSE = Vector2D(ZERO,Step_Height+HALF*Inlet_Height);
	xNW = Vector2D(-HALF*Inlet_Length,Step_Height+Inlet_Height);
	xNE = Vector2D(ZERO,Step_Height+Inlet_Height);
	BCtypeN = BC_WALL_VISCOUS_HEATFLUX;
	Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_I = 1.10;
	Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 2 && jBlk == 2) {
	xSW = Vector2D(ZERO,Step_Height+HALF*Inlet_Height);
	xSE = Vector2D(0.2*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xNW = Vector2D(ZERO,Step_Height+Inlet_Height);
	xNE = Vector2D(0.2*Outlet_Length,Step_Height+Inlet_Height+0.20*Outlet_Length*tan(Top_Wall_Deflection*PI/180.0));
	BCtypeN = BC_WALL_VISCOUS_HEATFLUX;
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.10;
	Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 3 && jBlk == 2) {
	xSW = Vector2D(0.2*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xSE = Vector2D(0.4*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xNW = Vector2D(0.2*Outlet_Length,Step_Height+Inlet_Height+0.20*Outlet_Length*tan(Top_Wall_Deflection*PI/180.0));
	xNE = Vector2D(0.4*Outlet_Length,Step_Height+Inlet_Height+0.40*Outlet_Length*tan(Top_Wall_Deflection*PI/180.0));
	BCtypeN = BC_WALL_VISCOUS_HEATFLUX;
	Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 4 && jBlk == 2) {
	xSW = Vector2D(0.4*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xSE = Vector2D(0.6*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xNW = Vector2D(0.4*Outlet_Length,Step_Height+Inlet_Height+0.40*Outlet_Length*tan(Top_Wall_Deflection*PI/180.0));
	xNE = Vector2D(0.6*Outlet_Length,Step_Height+Inlet_Height+0.60*Outlet_Length*tan(Top_Wall_Deflection*PI/180.0));
	BCtypeN = BC_WALL_VISCOUS_HEATFLUX;
	Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 5 && jBlk == 2) {
	xSW = Vector2D(0.6*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xSE = Vector2D(0.8*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xNW = Vector2D(0.6*Outlet_Length,Step_Height+Inlet_Height+0.60*Outlet_Length*tan(Top_Wall_Deflection*PI/180.0));
	xNE = Vector2D(0.8*Outlet_Length,Step_Height+Inlet_Height+0.80*Outlet_Length*tan(Top_Wall_Deflection*PI/180.0));
	BCtypeN = BC_WALL_VISCOUS_HEATFLUX;
	Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_J = 1.10;
      } else if (iBlk == 6 && jBlk == 2) {
	xSW = Vector2D(0.8*Outlet_Length,Step_Height+HALF*Inlet_Height);
	xSE = Vector2D(Outlet_Length,Step_Height+HALF*Inlet_Height);
	xNW = Vector2D(0.8*Outlet_Length,Step_Height+Inlet_Height+0.80*Outlet_Length*tan(Top_Wall_Deflection*PI/180.0));
	xNE = Vector2D(Outlet_Length,Step_Height+Inlet_Height+Outlet_Length*tan(Top_Wall_Deflection*PI/180.0));
	BCtypeN = BC_WALL_VISCOUS_HEATFLUX;
	BCtypeE = BC_OUTFLOW_SUBSONIC;
	Stretch_J = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_J = 1.10;

      }

      if (block_flag) {

	// Create the splines defining the north, south, east, and west 
	// boundaries of the rectangular boxes.
	Create_Spline_Line(Bnd_Spline_North,xNW,xNE,2);
	Create_Spline_Line(Bnd_Spline_South,xSW,xSE,2);
	Create_Spline_Line(Bnd_Spline_East,xSE,xNE,2);
	Create_Spline_Line(Bnd_Spline_West,xSW,xNW,2);

	// Set the boundary condition types for each of the boundary
	// splines:
	Bnd_Spline_North.setBCtype(BCtypeN);
	Bnd_Spline_South.setBCtype(BCtypeS);
	Bnd_Spline_East.setBCtype(BCtypeE);
	Bnd_Spline_West.setBCtype(BCtypeW);
    
	// Create the 2D quadrilateral grid block.
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

	// Deallocate the memory for the boundary splines.
	Bnd_Spline_North.deallocate();
	Bnd_Spline_South.deallocate();
	Bnd_Spline_East.deallocate();
	Bnd_Spline_West.deallocate();

      }

    }    
  }

  // Return the grid.
  return Grid_ptr;

}

/**********************************************************************
 * Routine: Grid_Desolvation_Chamber                                  *
 *                                                                    *
 * This routine creates a mesh corresponding to a simplified          *
 * desolvation chamber designed by MDS-SCIEX.                         *
 *                                                                    *
 **********************************************************************/
Grid2D_Quad_Block** Grid_Desolvation_Chamber(Grid2D_Quad_Block **Grid_ptr,
					     const int &Chamber_BC_Type,
					     int &Number_of_Blocks_Idir,
					     int &Number_of_Blocks_Jdir,
					     const int Number_of_Cells_Idir,
					     const int Number_of_Cells_Jdir,
					     const int Number_of_Ghost_Cells) {

  int block_flag,
      BCtypeN, BCtypeS, BCtypeE, BCtypeW,
      Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xNW, xNE, xSE, xSW;
  Spline2D Bnd_Spline_North, Bnd_Spline_South,
           Bnd_Spline_East,  Bnd_Spline_West, Bow_Spline;

  // Allocate memory for grid block.
  Number_of_Blocks_Idir = 7;
  Number_of_Blocks_Jdir = 6;
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
				       Number_of_Blocks_Idir, 
				       Number_of_Blocks_Jdir);


  // Create the mesh for each block representing the complete grid.
  for (int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++) {
    for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

      // Preset the block indicator and turn off when required:
      block_flag = ON;

      // Preset all boundary conditions to none and reset when required:
      BCtypeN = BC_NONE;
      BCtypeS = BC_NONE;
      BCtypeE = BC_NONE;
      BCtypeW = BC_NONE;

      // Preset all stretching parameters and reset when required:
      Stretch_I = STRETCHING_FCN_LINEAR;
      Beta_I = ZERO;
      Tau_I = ZERO;
      Stretch_J = STRETCHING_FCN_LINEAR;
      Beta_J = ZERO;
      Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;

      // Determine the block indicator and all boundary splines and
      // reset all boundary conditions and stretching parameters for the
      // current block: 
      if (iBlk == 0 && jBlk == 0) {
	xSW = Vector2D(-0.005000,0.000000);
	xSE = Vector2D( 0.000000,0.000000);
	xNW = Vector2D(-0.005000,0.000250);
	xNE = Vector2D( 0.000000,0.000250);
	BCtypeS = Chamber_BC_Type;
	BCtypeW = Chamber_BC_Type;
      } else if (iBlk == 1 && jBlk == 0) {
	xSW = Vector2D( 0.000000,0.000000);
	xSE = Vector2D( 0.006500,0.000000);
	xNW = Vector2D( 0.000000,0.000250);
	xNE = Vector2D( 0.006500,0.000250);
	BCtypeS = Chamber_BC_Type;
      } else if (iBlk == 2 && jBlk == 0) {
	xSW = Vector2D( 0.006500,0.000000);
	xSE = Vector2D( 0.013000,0.000000);
	xNW = Vector2D( 0.006500,0.000250);
	xNE = Vector2D( 0.013000,0.000250);
	BCtypeS = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_I = 1.01;
	Tau_I = ZERO;
      } else if (iBlk == 3 && jBlk == 0) {
	xSW = Vector2D( 0.013000,0.000000);
	xSE = Vector2D( 0.014000,0.000000);
	xNW = Vector2D( 0.013000,0.000250);
	xNE = Vector2D( 0.014000,0.000250);
	BCtypeS = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
	Beta_I = 1.05;
	Tau_I = ZERO;
      } else if (iBlk == 4 && jBlk == 0) {
	xSW = Vector2D( 0.014000,0.000000);
	xSE = Vector2D( 0.014250,0.000000);
	xNW = Vector2D( 0.014000,0.000250);
	xNE = Vector2D( 0.014250,0.000250);
	BCtypeN = Chamber_BC_Type;
	BCtypeS = Chamber_BC_Type;
      } else if (iBlk == 5 && jBlk == 0) {
	xSW = Vector2D( 0.014250,0.000000);
	xSE = Vector2D( 0.040000,0.000000);
	xNW = Vector2D( 0.014250,0.000250);
	xNE = Vector2D( 0.040000,0.000250);
	BCtypeS = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.001;
	Tau_I = ZERO;
      } else if (iBlk == 6 && jBlk == 0) {
	xSW = Vector2D( 0.040000,0.000000);
	xSE = Vector2D( 0.100000,0.000000);
	xNW = Vector2D( 0.040000,0.000250);
	xNE = Vector2D( 0.100000,0.000250);
	BCtypeS = Chamber_BC_Type;
	BCtypeE = BC_CONSTANT_EXTRAPOLATION;

      } else if (iBlk == 0 && jBlk == 1) {
	xSW = Vector2D(-0.005000,0.000250);
	xSE = Vector2D( 0.000000,0.000250);
	xNW = Vector2D(-0.005000,0.000500);
	xNE = Vector2D( 0.000000,0.000500);
	BCtypeW = Chamber_BC_Type;
      } else if (iBlk == 1 && jBlk == 1) {
	xSW = Vector2D( 0.000000,0.000250);
	xSE = Vector2D( 0.006500,0.000250);
	xNW = Vector2D( 0.000000,0.000500);
	xNE = Vector2D( 0.006500,0.000500);
	BCtypeN = Chamber_BC_Type;
      } else if (iBlk == 2 && jBlk == 1) {
	xSW = Vector2D( 0.006500,0.000250);
	xSE = Vector2D( 0.013000,0.000250);
	xNW = Vector2D( 0.006500,0.000500);
	xNE = Vector2D( 0.013000,0.000500);
	BCtypeN = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MAX_CLUSTERING;
	Beta_I = 1.01;
	Tau_I = ZERO;
      } else if (iBlk == 3 && jBlk == 1) {
	xSW = Vector2D( 0.013000,0.000250);
	xSE = Vector2D( 0.014000,0.000250);
	xNW = Vector2D( 0.013000,0.000500);
	xNE = Vector2D( 0.014000,0.000500);
	BCtypeE = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
	Beta_I = 1.05;
	Tau_I = ZERO;
      } else if (iBlk == 4 && jBlk == 1) {
	block_flag = OFF;
      } else if (iBlk == 5 && jBlk == 1) {
	xSW = Vector2D( 0.014250,0.000250);
	xSE = Vector2D( 0.040000,0.000250);
	xNW = Vector2D( 0.014250,0.000500);
	xNE = Vector2D( 0.040000,0.000500);
	BCtypeW = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.001;
	Tau_I = ZERO;
      } else if (iBlk == 6 && jBlk == 1) {
	xSW = Vector2D( 0.040000,0.000250);
	xSE = Vector2D( 0.100000,0.000250);
	xNW = Vector2D( 0.040000,0.000500);
	xNE = Vector2D( 0.100000,0.000500);
	BCtypeE = BC_CONSTANT_EXTRAPOLATION;

      } else if (iBlk == 0 && jBlk == 2) {
	xSW = Vector2D(-0.005000,0.000500);
	xSE = Vector2D( 0.000000,0.000500);
	xNW = Vector2D(-0.005000,0.000925);
	xNE = Vector2D( 0.000000,0.000925);
	BCtypeE = Chamber_BC_Type;
	BCtypeW = Chamber_BC_Type;
      } else if (iBlk == 1 && jBlk == 2) {
	block_flag = OFF;
      } else if (iBlk == 2 && jBlk == 2) {
	block_flag = OFF;
      } else if (iBlk == 3 && jBlk == 2) {
	xSW = Vector2D( 0.01300,0.0005000);
	xSE = Vector2D( 0.01400,0.0005000);
	xNW = Vector2D( 0.01300,0.0009250);
	xNE = Vector2D( 0.01400,0.0009250);
	BCtypeN = Chamber_BC_Type;
	BCtypeE = Chamber_BC_Type;
	BCtypeW = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MINMAX_CLUSTERING;
	Beta_I = 1.05;
	Tau_I = ZERO;
      } else if (iBlk == 4 && jBlk == 2) {
	block_flag = OFF;
      } else if (iBlk == 5 && jBlk == 2) {
	xSW = Vector2D( 0.014250,0.000500);
	xSE = Vector2D( 0.040000,0.000500);
	xNW = Vector2D( 0.014250,0.000925);
	xNE = Vector2D( 0.040000,0.000925);
	BCtypeW = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.001;
	Tau_I = ZERO;
      } else if (iBlk == 6 && jBlk == 2) {
	xSW = Vector2D( 0.040000,0.000500);
	xSE = Vector2D( 0.100000,0.000500);
	xNW = Vector2D( 0.040000,0.000925);
	xNE = Vector2D( 0.100000,0.000925);
	BCtypeE = BC_CONSTANT_EXTRAPOLATION;

      } else if (iBlk == 0 && jBlk == 3) {
	xSW = Vector2D(-0.005000,0.000925);
	xSE = Vector2D( 0.000000,0.000925);
	xNW = Vector2D(-0.005000,0.004500);
	xNE = Vector2D( 0.000000,0.004500);
	BCtypeE = Chamber_BC_Type;
	BCtypeW = Chamber_BC_Type;
      } else if (iBlk == 1 && jBlk == 3) {
	block_flag = OFF;
      } else if (iBlk == 2 && jBlk == 3) {
	block_flag = OFF;
      } else if (iBlk == 3 && jBlk == 3) {
	block_flag = OFF;
      } else if (iBlk == 4 && jBlk == 3) {
	block_flag = OFF;
      } else if (iBlk == 5 && jBlk == 3) {
	xSW = Vector2D( 0.014250,0.000925);
	xSE = Vector2D( 0.040000,0.000925);
	xNW = Vector2D( 0.014250,0.004500);
	xNE = Vector2D( 0.040000,0.004500);
	BCtypeW = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.001;
	Tau_I = ZERO;
      } else if (iBlk == 6 && jBlk == 3) {
	xSW = Vector2D( 0.040000,0.000925);
	xSE = Vector2D( 0.100000,0.000925);
	xNW = Vector2D( 0.040000,0.004500);
	xNE = Vector2D( 0.100000,0.004500);
	BCtypeE = BC_CONSTANT_EXTRAPOLATION;

      } else if (iBlk == 0 && jBlk == 4) {
	xSW = Vector2D(-0.005000,0.004500);
	xSE = Vector2D( 0.000000,0.004500);
	xNW = Vector2D(-0.005000,0.009000);
	xNE = Vector2D( 0.000000,0.009000);
	BCtypeN = Chamber_BC_Type;
	BCtypeE = Chamber_BC_Type;
	BCtypeW = Chamber_BC_Type;
      } else if (iBlk == 1 && jBlk == 4) {
	block_flag = OFF;
      } else if (iBlk == 2 && jBlk == 4) {
	block_flag = OFF;
      } else if (iBlk == 3 && jBlk == 4) {
	block_flag = OFF;
      } else if (iBlk == 4 && jBlk == 4) {
	block_flag = OFF;
      } else if (iBlk == 5 && jBlk == 4) {
	xSW = Vector2D( 0.014250,0.004500);
	xSE = Vector2D( 0.040000,0.004500);
	xNW = Vector2D( 0.014250,0.009000);
	xNE = Vector2D( 0.040000,0.009000);
	BCtypeW = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.001;
	Tau_I = ZERO;
      } else if (iBlk == 6 && jBlk == 4) {
	xSW = Vector2D( 0.040000,0.004500);
	xSE = Vector2D( 0.100000,0.004500);
	xNW = Vector2D( 0.040000,0.009000);
	xNE = Vector2D( 0.100000,0.009000);
	BCtypeE = BC_CONSTANT_EXTRAPOLATION;

      } else if (iBlk == 0 && jBlk == 5) {
	block_flag = OFF;
      } else if (iBlk == 1 && jBlk == 5) {
	block_flag = OFF;
      } else if (iBlk == 2 && jBlk == 5) {
	block_flag = OFF;
      } else if (iBlk == 3 && jBlk == 5) {
	block_flag = OFF;
      } else if (iBlk == 4 && jBlk == 5) {
	block_flag = OFF;
      } else if (iBlk == 5 && jBlk == 5) {
	xSW = Vector2D( 0.014250,0.009000);
	xSE = Vector2D( 0.040000,0.009000);
	xNW = Vector2D( 0.014250,0.070000);
	xNE = Vector2D( 0.040000,0.070000);
	BCtypeN = Chamber_BC_Type;
	BCtypeW = Chamber_BC_Type;
	Stretch_I = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_I = 1.001;
	Tau_I = ZERO;
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.015;
	Tau_J = ZERO;
      } else if (iBlk == 6 && jBlk == 5) {
	xSW = Vector2D( 0.040000,0.009000);
	xSE = Vector2D( 0.100000,0.009000);
	xNW = Vector2D( 0.040000,0.070000);
	xNE = Vector2D( 0.100000,0.070000);
	BCtypeN = Chamber_BC_Type;
	BCtypeE = BC_CONSTANT_EXTRAPOLATION;
	Stretch_J = STRETCHING_FCN_MIN_CLUSTERING;
	Beta_J = 1.015;
	Tau_J = ZERO;

      }

      if (block_flag) {

	// Create the splines defining the north, south, east, and west 
	// boundaries of the current block:
	Create_Spline_Line(Bnd_Spline_North,xNW,xNE,2);
	Create_Spline_Line(Bnd_Spline_South,xSW,xSE,2);
	Create_Spline_Line(Bnd_Spline_East,xSE,xNE,2);
	Create_Spline_Line(Bnd_Spline_West,xSW,xNW,2);

	// Set the boundary condition types for each of the boundary
	// splines:
	Bnd_Spline_North.setBCtype(BCtypeN);
	Bnd_Spline_South.setBCtype(BCtypeS);
	Bnd_Spline_East.setBCtype(BCtypeE);
	Bnd_Spline_West.setBCtype(BCtypeW);

	// Create the 2D quadrilateral grid block.
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

	// Deallocate the memory for the boundary splines.
	Bnd_Spline_North.deallocate();
	Bnd_Spline_South.deallocate();
	Bnd_Spline_East.deallocate();
	Bnd_Spline_West.deallocate();

      }

    }
  }

  // Return the grid.
  return Grid_ptr;

}

/**********************************************************************
 * Routine: Grid_NASA_Rotor_37                                        *
 *                                                                    *
 * This routine creates a mesh to be used in conjunction with and     *
 * embedded interface representing a section of a NASA rotor 37 blade *
 * at a given percent span.                                           *
 *                                                                    *
 *  +---------------+-------------------+---------------+             *
 *  |               |                   |               |             *
 *  |               |                   |               |             *
 *  |  Block (0,1)  |    Block (1,1)    |  Block (2,1)  |             *
 *  |               |                   |               |             *
 *  |               |                   |               |             *
 *  +---------------<<<<<<<<BLADE>>>>>>>>---------------+             *
 *  |               |                   |               |             *
 *  |               |                   |               |             *
 *  |  Block (0,0)  |    Block (1,0)    |  Block (2,0)  |             *
 *  |               |                   |               |             *
 *  |               |                   |               |             *
 *  +---------------+-------------------+---------------+             *
 *                                                                    *
 **********************************************************************/
Grid2D_Quad_Block** Grid_NASA_Rotor_37(Grid2D_Quad_Block **Grid_ptr,
				       int &Number_of_Blocks_Idir,
				       int &Number_of_Blocks_Jdir,
				       const double &Rotor_Percent_Span,
				       const int Number_of_Cells_Idir,
				       const int Number_of_Cells_Jdir,
				       const int Number_of_Ghost_Cells) {

  int error_flag;
  Spline2D upperB, lowerB, upperMiddleB, lowerMiddleB, camberTrail, camberLead, camberBlade, Rotor_Spline;
  int iTrail, iLead, pos, mIndex;
  Vector2D swap, leadV, trailV, x_temp;
  double zlU, zlL, zrU, zrL, mLead, mTrail, mLeadc, mTrailc, 
         mTop, mBot, m1, m2, z, dz, A, m, zLead, zTrail, dm, swapm;
  double beta_i, beta_j;

  // NASA Rotor 37 variables.
  char NASA_Rotor_Data_Directory[128];
  int Rotor_Flow_Type;
  NASARotor37 NASA_Rotor_37;

  // Initialize NASA rotor 37 class.
  strcpy(NASA_Rotor_Data_Directory,"/nfs/fe01/d1/cfd/jai/CFDkit+caboodle/data/NASA_Rotors/R37/");
  Rotor_Flow_Type = 1;//PEAK_FLOW;
  NASA_Rotor_37.init(Rotor_Flow_Type,NASA_Rotor_Data_Directory);

  // Allocate memory for grid block.
  Number_of_Blocks_Idir = 3;
  Number_of_Blocks_Jdir = 2;
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
				       Number_of_Blocks_Idir, 
				       Number_of_Blocks_Jdir);

  // 'Hard-wire' stretching parameters.
  if (Number_of_Cells_Idir/2 < 10) {
    beta_i = 1.05;
  } else if (Number_of_Cells_Idir/2 < 25) {
    beta_i = 1.01;
  } else if (Number_of_Cells_Idir/2 < 50) {
    beta_i = 1.005;
  } else if (Number_of_Cells_Idir/2 < 100) {
    beta_i = 1.0025;
  } else {
    beta_i = 1.001;
  } 
  if (Number_of_Cells_Jdir < 10) {
    beta_j = 1.05;
  } else if (Number_of_Cells_Jdir < 25) {
    beta_j = 1.01;
  } else if (Number_of_Cells_Jdir < 50) {
    beta_j = 1.005;
  } else if (Number_of_Cells_Jdir < 100) {
    beta_j = 1.0025;
  } else {
    beta_j = 1.001;
  }

  // Find the index for the leading and trailing edge.
  NASA_Rotor_37.findLT(Rotor_Percent_Span,iTrail,iLead);

  // Get cross section spline.
  Rotor_Spline = NASA_Rotor_37.getBladeCS(Rotor_Percent_Span);

  // Re-organize rotor spline points such that they ascend in the 
  // clockwise direction
  for (int i = 0; i < 100; i++) {
    swap = Rotor_Spline.Xp[2*100-2-i];
    Rotor_Spline.Xp[2*100-2-i] = Rotor_Spline.Xp[i];
    Rotor_Spline.Xp[i] = swap;
  }

  // Exchange iTrail and iLead accordingly.
  iTrail = Rotor_Spline.np-1 - iTrail;
  iLead = Rotor_Spline.np-1 - iLead;

  // Save leading and trailing edge vectors.
  leadV = Rotor_Spline.Xp[iLead];
  trailV = Rotor_Spline.Xp[iTrail];

  // Get camber line splines - 40 pts.
  camberTrail = NASA_Rotor_37.getCamberTrail(Rotor_Percent_Span,40,NASA_Rotor_37.z_d); 
  camberLead = NASA_Rotor_37.getCamberLead(Rotor_Percent_Span,40,NASA_Rotor_37.z_up);
  camberBlade = NASA_Rotor_37.getCamberBlade(Rotor_Percent_Span,40);
  
  // Get upper, lower, upperMiddle, and lowerMiddle boundaries.
  upperB = NASA_Rotor_37.getCamberAndShift(Rotor_Percent_Span,40,NASA_Rotor_37.z_up,NASA_Rotor_37.z_d,PI/NASA_Rotor_37.num_blades);
  lowerB = NASA_Rotor_37.getCamberAndShift(Rotor_Percent_Span,40,NASA_Rotor_37.z_up,NASA_Rotor_37.z_d,-PI/NASA_Rotor_37.num_blades);
  upperMiddleB = NASA_Rotor_37.getCamberAndShift(Rotor_Percent_Span,40,NASA_Rotor_37.z_up,NASA_Rotor_37.z_d,PI/(2*NASA_Rotor_37.num_blades));
  lowerMiddleB = NASA_Rotor_37.getCamberAndShift(Rotor_Percent_Span,40,NASA_Rotor_37.z_up,NASA_Rotor_37.z_d,-PI/(2*NASA_Rotor_37.num_blades));

  // Locate boundary corners to ensure that boundaries meeting at the leading
  // and trailing edges are separated by approximately 45deg
   
  // Calculate slope of leading and trailing camber lines.
  mLeadc = (camberLead.Xp[camberLead.np-1].y - camberLead.Xp[0].y)/
           (camberLead.Xp[camberLead.np-1].x - camberLead.Xp[0].x);

  mTrailc = (camberTrail.Xp[camberTrail.np-1].y-camberTrail.Xp[0].y)/
            (camberTrail.Xp[camberTrail.np-1].x-camberTrail.Xp[0].x);

  // Calculate slope of tangent to cross-section at the leading and 
  // trailing edge using a central-difference
  z = Rotor_Spline.Xp[iLead].x;
  dz = (Rotor_Spline.Xp[iLead].x - Rotor_Spline.Xp[iLead+1].x)/TEN;
  mLead = (-   getY(z+2*dz,Rotor_Spline)[0].y 
	   + 8*getY(z+  dz,Rotor_Spline)[0].y
	   - 8*getY(z-  dz,Rotor_Spline)[0].y + 
               getY(z-2*dz,Rotor_Spline)[0].y)/(12*dz);

  z = Rotor_Spline.Xp[iTrail].x;
  dz = (Rotor_Spline.Xp[iTrail].x - Rotor_Spline.Xp[iTrail-1].x)/TEN;
  mTrail = (-   getY(z+2*dz,Rotor_Spline)[1].y
	    + 8*getY(z+  dz,Rotor_Spline)[1].y 
	    - 8*getY(z-  dz,Rotor_Spline)[1].y
	    +   getY(z-2*dz,Rotor_Spline)[1].y)/(12*dz);

  // Get slopes of lines approx 45deg to camber line at the LEADING EDGE.
  A = (mLead+mLeadc)/(1-mLead*mLeadc);
  mTop = (-1+sqrt(1+A*A))/A;
  mBot = (-1-sqrt(1+A*A))/A;

  // Find point on upperMiddle spline - leading edge.
  dm = 1e15;
  mIndex = -1;
  for (int i = 0; i < upperMiddleB.np; i++) {
    if (fabs(upperMiddleB.Xp[i].x-leadV.x) > TOLER) {
      m = (upperMiddleB.Xp[i].y-leadV.y)/(upperMiddleB.Xp[i].x-leadV.x);
      if (fabs(m-mTop) < dm) {
	dm = fabs(m-mTop);
	mIndex = i;
      }
    }
  }
  assert(mIndex != -1);
  zlU = upperMiddleB.Xp[mIndex].x;  

  // Find point on lowerMiddle spline - leading edge.
  dm = 1e15;
  mIndex = -1;
  for (int i = 0; i < lowerMiddleB.np; i++) {
    if (fabs(lowerMiddleB.Xp[i].x-leadV.x) > TOLER) {
      m = (lowerMiddleB.Xp[i].y-leadV.y)/(lowerMiddleB.Xp[i].x-leadV.x);
      if (fabs(m-mBot) < dm) {
	dm=fabs(m-mBot);
	mIndex=i;
      }
    }
  }
  assert(mIndex != -1);
  zlL = lowerMiddleB.Xp[mIndex].x;  

  // Get slopes of lines approx 45deg to camber line at the TRAILING EDGE.
  A = (mTrail+mTrailc)/(1-mTrail*mTrailc);
  mBot = (-1+sqrt(1+A*A))/A;
  mTop = (-1-sqrt(1+A*A))/A;
  
  // Find point on upperMiddle spline - trailing edge.
  dm = 1e15;
  mIndex = -1;
  for (int i=0; i < upperMiddleB.np; i++) {
    if (fabs(upperMiddleB.Xp[i].x-trailV.x) > TOLER) {
      m = (upperMiddleB.Xp[i].y-trailV.y)/(upperMiddleB.Xp[i].x-trailV.x);
      if (fabs(m-mTop) < dm) {
	dm=fabs(m-mTop);
	mIndex=i;
      }
    }
  }
  assert(mIndex != -1);
  zrU = upperMiddleB.Xp[mIndex].x;  

  dm = 1e15;
  mIndex = -1;
  for (int i = 0; i < upperMiddleB.np; i++) {
    if (fabs(upperMiddleB.Xp[i].x-trailV.x) > TOLER) {
      m = (upperMiddleB.Xp[i].y-trailV.y)/(upperMiddleB.Xp[i].x-trailV.x);
      if (fabs(m-mBot) < dm) {
	dm = fabs(m-mBot);
	mIndex = i;
      }
    }
  }
  assert(mIndex != -1);
  
  // Determine which to use.
  if (upperMiddleB.Xp[mIndex].x>zrU) {
    zrU = upperMiddleB.Xp[mIndex].x;  
    swapm = mBot;
    mBot = mTop;
    mTop = swapm;
  }

  // Find point on lowerMiddle spline - trailing edge.
  dm = 1e15;
  mIndex = -1;
  for (int i=0; i < lowerMiddleB.np; i++) {
    if (fabs(lowerMiddleB.Xp[i].x-trailV.x) > TOLER) {
      m = (lowerMiddleB.Xp[i].y-trailV.y)/(lowerMiddleB.Xp[i].x-trailV.x);
      if (fabs(m-mBot) < dm) {
	dm = fabs(m-mBot);
	mIndex = i;
      }
    }
  }
  assert(mIndex != -1);
  zrL = lowerMiddleB.Xp[mIndex].x;

  /////////////////
  // BLOCK (1,1) //
  /////////////////
  // North spline.
  Grid_ptr[1][1].BndNorthSpline.allocate(70);
  Grid_ptr[1][1].BndNorthSpline.settype(SPLINE2D_QUINTIC);
  zLead = (zlL+zlU)/2;
  zTrail = (zrL+zrU)/2;
  dz = (zTrail-zLead)/69;
  for (int i = 0; i < 70; i++) {
    Grid_ptr[1][1].BndNorthSpline.Xp[i] = getY(zLead+i*dz,upperB)[0];
    Grid_ptr[1][1].BndNorthSpline.bc[i] = BC_NONE;
    if (i == 0 || i == Grid_ptr[1][1].BndNorthSpline.np-1) {
      Grid_ptr[1][1].BndNorthSpline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      Grid_ptr[1][1].BndNorthSpline.tp[i] = SPLINE2D_POINT_NORMAL;
    }
  }
  Grid_ptr[1][1].BndNorthSpline.pathlength();
  // South spline.
  Grid_ptr[1][1].BndSouthSpline.allocate(iTrail-iLead+1);
  Grid_ptr[1][1].BndSouthSpline.settype(SPLINE2D_QUINTIC);
  for (int i = iLead; i <= iTrail; i++)
    Grid_ptr[1][1].BndSouthSpline.Xp[i-iLead] = Rotor_Spline.Xp[i];
  Grid_ptr[1][1].BndSouthSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][1].BndSouthSpline.tp[Grid_ptr[1][1].BndSouthSpline.np-1] = SPLINE2D_POINT_SHARP_CORNER;
  for (int i = 1; i < Grid_ptr[1][1].BndSouthSpline.np-1; i++)
    Grid_ptr[1][1].BndSouthSpline.tp[i]=SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][1].BndSouthSpline.pathlength();
  // East spline.
  Grid_ptr[1][1].BndEastSpline.allocate(5);
  Grid_ptr[1][1].BndEastSpline.settype(SPLINE2D_CUBIC);
  Grid_ptr[1][1].BndEastSpline.Xp[0] = Rotor_Spline.Xp[iTrail];
  x_temp = getY(zrU,upperMiddleB)[0];
  Grid_ptr[1][1].BndEastSpline.Xp[1] = Rotor_Spline.Xp[iTrail]+HALF*(x_temp-Rotor_Spline.Xp[iTrail]);
  Grid_ptr[1][1].BndEastSpline.Xp[2] = x_temp;
  Grid_ptr[1][1].BndEastSpline.Xp[3] = x_temp + HALF*(Grid_ptr[1][1].BndNorthSpline.Xp[Grid_ptr[1][1].BndNorthSpline.np-1] - x_temp);
  Grid_ptr[1][1].BndEastSpline.Xp[4] = Grid_ptr[1][1].BndNorthSpline.Xp[Grid_ptr[1][1].BndNorthSpline.np-1];
  Grid_ptr[1][1].BndEastSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][1].BndEastSpline.tp[1] = SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][1].BndEastSpline.tp[2] = SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][1].BndEastSpline.tp[3] = SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][1].BndEastSpline.tp[4] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][1].BndEastSpline.pathlength();
  // West spline.
  Grid_ptr[1][1].BndWestSpline.allocate(6);
  Grid_ptr[1][1].BndWestSpline.settype(SPLINE2D_CUBIC);
  Grid_ptr[1][1].BndWestSpline.Xp[0] = Rotor_Spline.Xp[iLead];
  x_temp = getY(zlU,upperMiddleB)[0];
  Grid_ptr[1][1].BndWestSpline.Xp[1] = Rotor_Spline.Xp[iLead]+0.3333*(x_temp-Rotor_Spline.Xp[iLead]);
  Grid_ptr[1][1].BndWestSpline.Xp[2] = HALF*(Rotor_Spline.Xp[iLead]+HALF*(x_temp-Rotor_Spline.Xp[iLead])+
			               x_temp+0.20*(Grid_ptr[1][1].BndNorthSpline.Xp[0]-x_temp));
  Grid_ptr[1][1].BndWestSpline.Xp[3] = HALF*(Rotor_Spline.Xp[iLead]+0.6666*(x_temp-Rotor_Spline.Xp[iLead])+
			               x_temp+0.3333*(Grid_ptr[1][1].BndNorthSpline.Xp[0]-x_temp));
  Grid_ptr[1][1].BndWestSpline.Xp[4] = x_temp+HALF*(Grid_ptr[1][1].BndNorthSpline.Xp[0]-x_temp);
  Grid_ptr[1][1].BndWestSpline.Xp[5] = Grid_ptr[1][1].BndNorthSpline.Xp[0];
  Grid_ptr[1][1].BndWestSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][1].BndWestSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][1].BndWestSpline.tp[2] = SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][1].BndWestSpline.tp[3] = SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][1].BndWestSpline.tp[4] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][1].BndWestSpline.tp[5] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][1].BndWestSpline.pathlength();
  // Create quad block (1,1).
  Create_Quad_Block(Grid_ptr[1][1],
		    Grid_ptr[1][1].BndNorthSpline,
		    Grid_ptr[1][1].BndSouthSpline,
		    Grid_ptr[1][1].BndEastSpline,
		    Grid_ptr[1][1].BndWestSpline,
		    Number_of_Cells_Idir,
		    Number_of_Cells_Jdir,
		    Number_of_Ghost_Cells,
		    GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_XY,
		    STRETCHING_FCN_MINMAX_CLUSTERING,
		    ONE + (beta_i-ONE)/ONE,
		    0,
		    STRETCHING_FCN_MIN_CLUSTERING,
		    beta_j,
		    0,
		    ORTHOGONAL,
		    ORTHOGONAL,
		    ORTHOGONAL,
		    ORTHOGONAL);
  // Smooth quad block (1,1).
  Smooth_Quad_Block(Grid_ptr[1][1],
		    min(250,2*max(Number_of_Cells_Idir,
				  Number_of_Cells_Jdir)));

  /////////////////
  // BLOCK (1,0) //
  /////////////////
  // North spline.
  Grid_ptr[1][0].BndNorthSpline.allocate(Rotor_Spline.np-iTrail+iLead);
  Grid_ptr[1][0].BndNorthSpline.settype(SPLINE2D_QUINTIC);
  for (int i = iLead; i >= 0; i--)
    Grid_ptr[1][0].BndNorthSpline.Xp[iLead-i] = Rotor_Spline.Xp[i];
  for (int i=Rotor_Spline.np-2; i >= iTrail; i--)
    Grid_ptr[1][0].BndNorthSpline.Xp[iLead+1+(Rotor_Spline.np-2)-i] = Rotor_Spline.Xp[i];
  for (int i = 0; i < Grid_ptr[1][0].BndNorthSpline.np; i++) {
    Grid_ptr[1][0].BndNorthSpline.bc[i] = BC_NONE;
    if (i == 0 || i == Grid_ptr[1][0].BndNorthSpline.np-1) {
      Grid_ptr[1][0].BndNorthSpline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      Grid_ptr[1][0].BndNorthSpline.tp[i] = SPLINE2D_POINT_NORMAL;
    }
  }
  Grid_ptr[1][0].BndNorthSpline.pathlength();  
  // South spline.
  Grid_ptr[1][0].BndSouthSpline.allocate(70);
  Grid_ptr[1][0].BndSouthSpline.settype(SPLINE2D_QUINTIC);
  zLead = (zlL+zlU)/2;;
  zTrail = (zrL+zrU)/2;;
  dz = (zTrail-zLead)/69;
  for (int i = 0; i < 70; i++) 
    Grid_ptr[1][0].BndSouthSpline.Xp[i] = getY(zLead+i*dz,lowerB)[0];
  for (int i = 0; i < Grid_ptr[1][0].BndSouthSpline.np; i++) {
    Grid_ptr[1][0].BndSouthSpline.bc[i] = BC_NONE;
    if (i == 0 || i == Grid_ptr[1][0].BndSouthSpline.np-1) {
      Grid_ptr[1][0].BndSouthSpline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      Grid_ptr[1][0].BndSouthSpline.tp[i] = SPLINE2D_POINT_NORMAL;
    }
  }
  Grid_ptr[1][0].BndSouthSpline.pathlength();
  // West spline.
  Grid_ptr[1][0].BndWestSpline.allocate(5);
  Grid_ptr[1][0].BndWestSpline.settype(SPLINE2D_CUBIC);
  Grid_ptr[1][0].BndWestSpline.Xp[0] = Grid_ptr[1][0].BndSouthSpline.Xp[0];
  x_temp = getY(zlL,lowerMiddleB)[0];
  Grid_ptr[1][0].BndWestSpline.Xp[1] = Grid_ptr[1][0].BndSouthSpline.Xp[0]+
                                       HALF*(x_temp-Grid_ptr[1][0].BndSouthSpline.Xp[0]);
  Grid_ptr[1][0].BndWestSpline.Xp[2] = x_temp;
  Grid_ptr[1][0].BndWestSpline.Xp[3] = x_temp+HALF*(Rotor_Spline.Xp[iLead]-x_temp);
  Grid_ptr[1][0].BndWestSpline.Xp[4] = Rotor_Spline.Xp[iLead];
  for (int i = 0; i < 5; i++) Grid_ptr[1][0].BndWestSpline.bc[i] = BC_NONE;
  Grid_ptr[1][0].BndWestSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][0].BndWestSpline.tp[1] = SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][0].BndWestSpline.tp[2] = SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][0].BndWestSpline.tp[3] = SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][0].BndWestSpline.tp[4] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][0].BndWestSpline.pathlength();
  // East spline.
  Grid_ptr[1][0].BndEastSpline.allocate(6);
  Grid_ptr[1][0].BndEastSpline.settype(SPLINE2D_CUBIC);
  Grid_ptr[1][0].BndEastSpline.Xp[0] = Grid_ptr[1][0].BndSouthSpline.Xp[Grid_ptr[1][0].BndSouthSpline.np-1];
  x_temp = getY(zrL,lowerMiddleB)[0];
  Grid_ptr[1][0].BndEastSpline.Xp[1] = x_temp + HALF*(Grid_ptr[1][0].BndSouthSpline.Xp[Grid_ptr[1][0].BndSouthSpline.np-1]-x_temp);
  Grid_ptr[1][0].BndEastSpline.Xp[2] = HALF*(Rotor_Spline.Xp[iTrail]+0.6666*(x_temp-Rotor_Spline.Xp[iTrail])+
                                       x_temp+0.3333*(Grid_ptr[1][0].BndSouthSpline.Xp[Grid_ptr[1][0].BndSouthSpline.np-1]-x_temp));
  Grid_ptr[1][0].BndEastSpline.Xp[3] = 0.95*(Rotor_Spline.Xp[iTrail]+HALF*(x_temp-Rotor_Spline.Xp[iTrail]))+
                                       0.05*(x_temp+0.80*(Grid_ptr[1][0].BndSouthSpline.Xp[Grid_ptr[1][0].BndSouthSpline.np-1]-x_temp));
  Grid_ptr[1][0].BndEastSpline.Xp[4] = Rotor_Spline.Xp[iTrail]+0.3333*(x_temp-Rotor_Spline.Xp[iTrail]);
  Grid_ptr[1][0].BndEastSpline.Xp[5] = Rotor_Spline.Xp[iTrail];
  for (int i = 0; i < 6; i++) Grid_ptr[1][0].BndEastSpline.bc[i] = BC_NONE;
  Grid_ptr[1][0].BndEastSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][0].BndEastSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][0].BndEastSpline.tp[2] = SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][0].BndEastSpline.tp[3] = SPLINE2D_POINT_NORMAL;
  Grid_ptr[1][0].BndEastSpline.tp[4] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][0].BndEastSpline.tp[5] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[1][0].BndEastSpline.pathlength();
  // Create quad block (1,0).
  Create_Quad_Block(Grid_ptr[1][0],
		    Grid_ptr[1][0].BndNorthSpline,
		    Grid_ptr[1][0].BndSouthSpline,
		    Grid_ptr[1][0].BndEastSpline,
		    Grid_ptr[1][0].BndWestSpline,
		    Number_of_Cells_Idir,
		    Number_of_Cells_Jdir,
		    Number_of_Ghost_Cells,
		    GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_XY,
		    STRETCHING_FCN_MINMAX_CLUSTERING,
		    ONE + (beta_i-ONE)/ONE,
		    0,
		    STRETCHING_FCN_MAX_CLUSTERING,
		    beta_j,
		    0,
		    ORTHOGONAL,
		    ORTHOGONAL,
		    ORTHOGONAL,
		    ORTHOGONAL);
  // Smooth quad block (1,0).
  Smooth_Quad_Block(Grid_ptr[1][0],
		    min(250,2*max(Number_of_Cells_Idir,
				  Number_of_Cells_Jdir)));

  /////////////////
  // BLOCK (0,0) //
  /////////////////
  // North spline.
  Grid_ptr[0][0].BndNorthSpline.allocate(2);
  Grid_ptr[0][0].BndNorthSpline.settype(SPLINE2D_LINEAR);
  Grid_ptr[0][0].BndNorthSpline.Xp[0] = camberLead.Xp[0];
  Grid_ptr[0][0].BndNorthSpline.Xp[1] = Rotor_Spline.Xp[iLead];
  Grid_ptr[0][0].BndNorthSpline.bc[0] = BC_NONE;
  Grid_ptr[0][0].BndNorthSpline.bc[1] = BC_NONE;
  Grid_ptr[0][0].BndNorthSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[0][0].BndNorthSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[0][0].BndNorthSpline.pathlength();
  // South spline.
  Grid_ptr[0][0].BndSouthSpline.allocate(30);
  Grid_ptr[0][0].BndSouthSpline.settype(SPLINE2D_QUINTIC);
  zLead = lowerB.Xp[0].x;
  zTrail = Grid_ptr[1][0].BndSouthSpline.Xp[0].x;
  dz = (zTrail-zLead)/29;
  for (int i = 0; i < 30; i++)
    Grid_ptr[0][0].BndSouthSpline.Xp[i] = getY(zLead+i*dz, lowerB)[0];
  Grid_ptr[0][0].BndSouthSpline.Xp[0] = lowerB.Xp[0];
  for (int i = 0; i < Grid_ptr[0][0].BndSouthSpline.np; i++) {
    Grid_ptr[0][0].BndSouthSpline.bc[i] = BC_NONE;
    if (i == 0 || i == Grid_ptr[0][0].BndSouthSpline.np-1) {
      Grid_ptr[0][0].BndSouthSpline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      Grid_ptr[0][0].BndSouthSpline.tp[i] = SPLINE2D_POINT_NORMAL;
    }
  }
  Grid_ptr[0][0].BndSouthSpline.pathlength();
  // East spline.
  Copy_Spline(Grid_ptr[0][0].BndEastSpline,Grid_ptr[1][0].BndWestSpline);
  // West spline.
  Grid_ptr[0][0].BndWestSpline.allocate(2);
  Grid_ptr[0][0].BndWestSpline.settype(SPLINE2D_LINEAR);
  Grid_ptr[0][0].BndWestSpline.Xp[0] = Grid_ptr[0][0].BndSouthSpline.Xp[0];
  Grid_ptr[0][0].BndWestSpline.Xp[1] = Grid_ptr[0][0].BndNorthSpline.Xp[0];
  Grid_ptr[0][0].BndWestSpline.bc[0] = BC_FIXED;
  Grid_ptr[0][0].BndWestSpline.bc[1] = BC_FIXED;
  Grid_ptr[0][0].BndWestSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[0][0].BndWestSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[0][0].BndWestSpline.pathlength();
  // Create quad block (0,0).
  Create_Quad_Block(Grid_ptr[0][0],
		    Grid_ptr[0][0].BndNorthSpline,
		    Grid_ptr[0][0].BndSouthSpline,
		    Grid_ptr[0][0].BndEastSpline,
		    Grid_ptr[0][0].BndWestSpline,
		    Number_of_Cells_Idir,
		    Number_of_Cells_Jdir,
		    Number_of_Ghost_Cells,
		    GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST,
		    STRETCHING_FCN_MAX_CLUSTERING,
		    beta_i,
		    0,
		    STRETCHING_FCN_MAX_CLUSTERING,
		    beta_j,
		    0,
		    ORTHOGONAL,
		    ORTHOGONAL,
		    ORTHOGONAL,
		    NOT_ORTHOGONAL);
  // Smooth quad block (0,0).
  Smooth_Quad_Block(Grid_ptr[0][0],
		    min(250,2*max(Number_of_Cells_Idir,
				  Number_of_Cells_Jdir)));

  /////////////////
  // BLOCK (2,0) //
  /////////////////
  // North spline.
  Grid_ptr[2][0].BndNorthSpline.allocate(2);
  Grid_ptr[2][0].BndNorthSpline.settype(SPLINE2D_LINEAR);
  Grid_ptr[2][0].BndNorthSpline.Xp[0] = Rotor_Spline.Xp[iTrail];
  Grid_ptr[2][0].BndNorthSpline.Xp[1] = camberTrail.Xp[camberTrail.np-1];
  Grid_ptr[2][0].BndNorthSpline.bc[0] = BC_NONE;
  Grid_ptr[2][0].BndNorthSpline.bc[1] = BC_NONE;
  Grid_ptr[2][0].BndNorthSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[2][0].BndNorthSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[2][0].BndNorthSpline.pathlength();
  // South spline.
  Grid_ptr[2][0].BndSouthSpline.allocate(30);
  Grid_ptr[2][0].BndSouthSpline.settype(SPLINE2D_QUINTIC);
  zLead = (zrL+zrU)/2;;
  zTrail = lowerB.Xp[lowerB.np-1].x;
  dz = (zTrail-zLead)/29;
  for (int i = 0; i < 30; i++) {
    Grid_ptr[2][0].BndSouthSpline.Xp[i].x = zLead+i*dz;
    if(Grid_ptr[2][0].BndSouthSpline.Xp[i].x>zTrail) Grid_ptr[2][0].BndSouthSpline.Xp[i].x = zTrail;
    Grid_ptr[2][0].BndSouthSpline.Xp[i].y = getY(Grid_ptr[2][0].BndSouthSpline.Xp[i].x,lowerB)[0].y;
  }
  Grid_ptr[2][0].BndSouthSpline.Xp[Grid_ptr[2][0].BndSouthSpline.np-1] = lowerB.Xp[lowerB.np-1];
  for (int i = 0; i < Grid_ptr[2][0].BndSouthSpline.np; i++) {
    Grid_ptr[2][0].BndSouthSpline.tp[i] = BC_NONE;
    if (i == 0 || i == Grid_ptr[2][0].BndSouthSpline.np-1) {
      Grid_ptr[2][0].BndSouthSpline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      Grid_ptr[2][0].BndSouthSpline.tp[i] = SPLINE2D_POINT_NORMAL;
    }
  }
  Grid_ptr[2][0].BndSouthSpline.pathlength();
  // East spline.
  Grid_ptr[2][0].BndEastSpline.allocate(2);
  Grid_ptr[2][0].BndEastSpline.settype(SPLINE2D_LINEAR);
  Grid_ptr[2][0].BndEastSpline.Xp[0] = Grid_ptr[2][0].BndSouthSpline.Xp[Grid_ptr[2][0].BndSouthSpline.np-1];
  Grid_ptr[2][0].BndEastSpline.Xp[1] = Grid_ptr[2][0].BndNorthSpline.Xp[Grid_ptr[2][0].BndNorthSpline.np-1];
  Grid_ptr[2][0].BndEastSpline.bc[0] = BC_CHARACTERISTIC;
  Grid_ptr[2][0].BndEastSpline.bc[1] = BC_CHARACTERISTIC;
  Grid_ptr[2][0].BndEastSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[2][0].BndEastSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[2][0].BndEastSpline.pathlength();
  // West spline.
  Copy_Spline(Grid_ptr[2][0].BndWestSpline,Grid_ptr[1][0].BndEastSpline); 
  // Create quad block (2,0).
  Create_Quad_Block(Grid_ptr[2][0],
		    Grid_ptr[2][0].BndNorthSpline,
		    Grid_ptr[2][0].BndSouthSpline,
		    Grid_ptr[2][0].BndEastSpline,
		    Grid_ptr[2][0].BndWestSpline,
		    Number_of_Cells_Idir,
		    Number_of_Cells_Jdir,
		    Number_of_Ghost_Cells,
		    GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST,
		    STRETCHING_FCN_MIN_CLUSTERING,
		    beta_i,
		    0,
		    STRETCHING_FCN_MAX_CLUSTERING,
		    beta_j,
		    0,
		    ORTHOGONAL,//NOT_
		    ORTHOGONAL,
		    NOT_ORTHOGONAL,
		    ORTHOGONAL);
  // Smooth quad block (2,0).
  Smooth_Quad_Block(Grid_ptr[2][0],
		    min(250,2*max(Number_of_Cells_Idir,
				  Number_of_Cells_Jdir)));

  /////////////////
  // BLOCK (0,1) //
  /////////////////
  // North spline.
  Grid_ptr[0][1].BndNorthSpline.allocate(30);
  Grid_ptr[0][1].BndNorthSpline.settype(SPLINE2D_QUINTIC);
  zLead = upperB.Xp[0].x;
  zTrail = (zlL+zlU)/2;
  dz = (zTrail-zLead)/29;
  for (int i = 0; i < 30; i++)
    Grid_ptr[0][1].BndNorthSpline.Xp[i] = getY(zLead+i*dz,upperB)[0];
  Grid_ptr[0][1].BndNorthSpline.Xp[0] = upperB.Xp[0];
  for (int i = 0; i < Grid_ptr[0][1].BndNorthSpline.np; i++) {
    Grid_ptr[0][1].BndNorthSpline.bc[i] = BC_NONE;
    if (i == 0 || i == Grid_ptr[0][1].BndNorthSpline.np-1) {
      Grid_ptr[0][1].BndNorthSpline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      Grid_ptr[0][1].BndNorthSpline.tp[i] = SPLINE2D_POINT_NORMAL;
    }
  }
  Grid_ptr[0][1].BndNorthSpline.pathlength();
  // South spline.
  Copy_Spline(Grid_ptr[0][1].BndSouthSpline,Grid_ptr[0][0].BndNorthSpline);
  // East spline.
  Copy_Spline(Grid_ptr[0][1].BndEastSpline,Grid_ptr[1][1].BndWestSpline);
  // West spline.
  Grid_ptr[0][1].BndWestSpline.allocate(2);
  Grid_ptr[0][1].BndWestSpline.settype(SPLINE2D_LINEAR);
  Grid_ptr[0][1].BndWestSpline.Xp[0] = Grid_ptr[0][1].BndSouthSpline.Xp[0];
  Grid_ptr[0][1].BndWestSpline.Xp[1] = Grid_ptr[0][1].BndNorthSpline.Xp[0];
  Grid_ptr[0][1].BndWestSpline.bc[0] = BC_FIXED;
  Grid_ptr[0][1].BndWestSpline.bc[1] = BC_FIXED;
  Grid_ptr[0][1].BndWestSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[0][1].BndWestSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[0][1].BndWestSpline.pathlength();
  // Create quad block (0,1).
  Create_Quad_Block(Grid_ptr[0][1],
		    Grid_ptr[0][1].BndNorthSpline,
		    Grid_ptr[0][1].BndSouthSpline,
		    Grid_ptr[0][1].BndEastSpline,
		    Grid_ptr[0][1].BndWestSpline,
		    Number_of_Cells_Idir,
		    Number_of_Cells_Jdir,
		    Number_of_Ghost_Cells,
		    GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST,
		    STRETCHING_FCN_MAX_CLUSTERING,
		    beta_i,
		    0,
		    STRETCHING_FCN_MIN_CLUSTERING,
		    beta_j,
		    0,
		    ORTHOGONAL,
		    ORTHOGONAL,
		    ORTHOGONAL,
		    NOT_ORTHOGONAL);
  // Smooth quad block (0,1).
  Smooth_Quad_Block(Grid_ptr[0][1],
		    min(250,2*max(Number_of_Cells_Idir,
				  Number_of_Cells_Jdir)));

  /////////////////
  // BLOCK (2,1) //
  /////////////////
  // North spline.
  Grid_ptr[2][1].BndNorthSpline.allocate(30);
  Grid_ptr[2][1].BndNorthSpline.settype(SPLINE2D_QUINTIC);
  zLead = (zrL+zrU)/2;
  zTrail = upperB.Xp[upperB.np-1].x;
  dz = (zTrail-zLead)/29;
  for (int i = 0; i < 30; i++) {
    Grid_ptr[2][1].BndNorthSpline.Xp[i].x = zLead+i*dz;
    if (Grid_ptr[2][1].BndNorthSpline.Xp[i].x > zTrail)
      Grid_ptr[2][1].BndNorthSpline.Xp[i].x = zTrail;
    Grid_ptr[2][1].BndNorthSpline.Xp[i].y = getY(Grid_ptr[2][1].BndNorthSpline.Xp[i].x,upperB)[0].y;
  }
  Grid_ptr[2][1].BndNorthSpline.Xp[Grid_ptr[2][1].BndNorthSpline.np-1] = upperB.Xp[upperB.np-1];
  for (int i = 0; i < Grid_ptr[2][1].BndNorthSpline.np; i++) {
    Grid_ptr[2][1].BndNorthSpline.bc[i] = BC_NONE;
    if (i == 0 || i == Grid_ptr[2][1].BndNorthSpline.np-1) {
      Grid_ptr[2][1].BndNorthSpline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      Grid_ptr[2][1].BndNorthSpline.tp[i] = SPLINE2D_POINT_NORMAL;
    }
  }
  Grid_ptr[2][1].BndNorthSpline.pathlength();
  // South spline.
  Copy_Spline(Grid_ptr[2][1].BndSouthSpline,Grid_ptr[2][0].BndNorthSpline);
  // East spline.
  Grid_ptr[2][1].BndEastSpline.allocate(2);
  Grid_ptr[2][1].BndEastSpline.settype(SPLINE2D_LINEAR);
  Grid_ptr[2][1].BndEastSpline.Xp[0] = Grid_ptr[2][1].BndSouthSpline.Xp[Grid_ptr[2][1].BndSouthSpline.np-1];
  Grid_ptr[2][1].BndEastSpline.Xp[1] = Grid_ptr[2][1].BndNorthSpline.Xp[Grid_ptr[2][1].BndNorthSpline.np-1];
  Grid_ptr[2][1].BndEastSpline.bc[0] = BC_CHARACTERISTIC;
  Grid_ptr[2][1].BndEastSpline.bc[1] = BC_CHARACTERISTIC;
  Grid_ptr[2][1].BndEastSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[2][1].BndEastSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  Grid_ptr[2][1].BndEastSpline.pathlength();
  // West spline.
  Copy_Spline(Grid_ptr[2][1].BndWestSpline,Grid_ptr[1][1].BndEastSpline);
  // Create quad block (2,1).
  Create_Quad_Block(Grid_ptr[2][1],
		    Grid_ptr[2][1].BndNorthSpline,
		    Grid_ptr[2][1].BndSouthSpline,
		    Grid_ptr[2][1].BndEastSpline,
		    Grid_ptr[2][1].BndWestSpline,
		    Number_of_Cells_Idir,
		    Number_of_Cells_Jdir,
		    Number_of_Ghost_Cells,
		    GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST,
		    STRETCHING_FCN_MIN_CLUSTERING,
		    beta_i,
		    0,
		    STRETCHING_FCN_MIN_CLUSTERING,
		    beta_j,
		    0,
		    ORTHOGONAL,
		    ORTHOGONAL,
		    NOT_ORTHOGONAL,
		    ORTHOGONAL);
  // Smooth quad block (2,1).
  Smooth_Quad_Block(Grid_ptr[2][1],
		    min(250,2*max(Number_of_Cells_Idir,
				  Number_of_Cells_Jdir)));

  // Deallocate splines.
  upperB.deallocate();
  lowerB.deallocate();
  upperMiddleB.deallocate();
  lowerMiddleB.deallocate();
  camberTrail.deallocate();
  camberLead.deallocate();
  camberBlade.deallocate();
  Rotor_Spline.deallocate();

  // Return the grid.
  return Grid_ptr;

}

/**********************************************************************
 * Routine: Grid_NASA_Rotor_67                                        *
 *                                                                    *
 * This routine creates a mesh to be used in conjunction with and     *
 * embedded interface representing a section of a NASA rotor 67 blade *
 * at a given percent span.                                           *
 *                                                                    *
 *  +---------------+-------------------+---------------+             *
 *  |               |                   |               |             *
 *  |               |                   |               |             *
 *  |  Block (0,1)  |    Block (1,1)    |  Block (2,1)  |             *
 *  |               |                   |               |             *
 *  |               |                   |               |             *
 *  +---------------<<<<<<<<BLADE>>>>>>>>---------------+             *
 *  |               |                   |               |             *
 *  |               |                   |               |             *
 *  |  Block (0,0)  |    Block (1,0)    |  Block (2,0)  |             *
 *  |               |                   |               |             *
 *  |               |                   |               |             *
 *  +---------------+-------------------+---------------+             *
 *                                                                    *
 **********************************************************************/
Grid2D_Quad_Block** Grid_NASA_Rotor_67(Grid2D_Quad_Block **Grid_ptr,
				       int &Number_of_Blocks_Idir,
				       int &Number_of_Blocks_Jdir,
				       const double &Rotor_Percent_Span,
				       const int Number_of_Cells_Idir,
				       const int Number_of_Cells_Jdir,
				       const int Number_of_Ghost_Cells) {

  int error_flag;
  Spline2D upperB, lowerB, upperMiddleB, lowerMiddleB, camberTrail, camberLead, camberBlade, Rotor_Spline;
  int iTrail, iLead, pos, mIndex;
  Vector2D swap, leadV, trailV, x_temp;
  double zlU, zlL, zrU, zrL, mLead, mTrail, mLeadc, mTrailc, 
         mTop, mBot, m1, m2, z, dz, A, m, zLead, zTrail, dm, swapm;
  Spline2D BndNorthSpline, tempBndNorthSpline, BndSouthSpline, tempBndSouthSpline, BndEastSpline, BndWestSpline;
  int Stretch_I, Stretch_J,
      Orthogonal_North, Orthogonal_South,
      Orthogonal_East, Orthogonal_West;
  double Beta_I, Tau_I, Beta_J, Tau_J;
  Vector2D xNW, xNE, xSE, xSW;
  Spline2D Bnd_Spline_North, Bnd_Spline_South,
           Bnd_Spline_East,  Bnd_Spline_West, Bow_Spline;

  // NASA Rotor 67 variables.
  char NASA_Rotor_Data_Directory[128];
  int Rotor_Flow_Type;
  NASARotor67 NASA_Rotor_67;

  // Initialize NASA rotor 67 class.
  strcpy(NASA_Rotor_Data_Directory,"/nfs/fe01/d1/cfd/jai/CFDkit+caboodle/data/NASA_Rotors/R67/");
  Rotor_Flow_Type = 1;//PEAK_FLOW;
  NASA_Rotor_67.init(Rotor_Flow_Type,NASA_Rotor_Data_Directory);

  // Allocate memory for grid block.
  Number_of_Blocks_Idir = 3;
  Number_of_Blocks_Jdir = 2;
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
				       Number_of_Blocks_Idir, 
				       Number_of_Blocks_Jdir);

  // Find the index for the leading and trailing edge.
  NASA_Rotor_67.findLT(Rotor_Percent_Span,iTrail,iLead);

  // Get cross section spline.
  Rotor_Spline = NASA_Rotor_67.getBladeCS(Rotor_Percent_Span);

  // Re-organize rotor spline points such that they ascend in the 
  // clockwise direction
  for (int i = 0; i < 70; i++ ){
    swap = Rotor_Spline.Xp[2*70-2-i];
    Rotor_Spline.Xp[2*70-2-i] = Rotor_Spline.Xp[i];
    Rotor_Spline.Xp[i] = swap;
  }

  // Exchange iTrail and iLead accordingly.
  iTrail = Rotor_Spline.np-1 - iTrail;
  iLead = Rotor_Spline.np-1 - iLead;

  // Save leading and trailing edge vectors.
  leadV = Rotor_Spline.Xp[iLead];
  trailV = Rotor_Spline.Xp[iTrail];

  // Get camber line splines - 40 pts.
  camberTrail = NASA_Rotor_67.getCamberTrail(Rotor_Percent_Span,40,NASA_Rotor_67.z_d); 
  camberLead = NASA_Rotor_67.getCamberLead(Rotor_Percent_Span,40,NASA_Rotor_67.z_up);
  camberBlade = NASA_Rotor_67.getCamberBlade(Rotor_Percent_Span,40);
  
  // Get upper, lower, upperMiddle, and lowerMiddle boundaries.
  upperB = NASA_Rotor_67.getCamberAndShift(Rotor_Percent_Span,40,NASA_Rotor_67.z_up,NASA_Rotor_67.z_d,PI/NASA_Rotor_67.num_blades);
  lowerB = NASA_Rotor_67.getCamberAndShift(Rotor_Percent_Span,40,NASA_Rotor_67.z_up,NASA_Rotor_67.z_d,-PI/NASA_Rotor_67.num_blades);
  upperMiddleB = NASA_Rotor_67.getCamberAndShift(Rotor_Percent_Span,40,NASA_Rotor_67.z_up,NASA_Rotor_67.z_d,PI/(2*NASA_Rotor_67.num_blades));
  lowerMiddleB = NASA_Rotor_67.getCamberAndShift(Rotor_Percent_Span,40,NASA_Rotor_67.z_up,NASA_Rotor_67.z_d,-PI/(2*NASA_Rotor_67.num_blades));

  // Locate boundary corners to ensure that boundaries meeting at the leading
  // and trailing edges are separated by approximately 45deg
   
  // Calculate slope of leading and trailing camber lines.
  mLeadc = (camberLead.Xp[camberLead.np-1].y - camberLead.Xp[0].y)/
           (camberLead.Xp[camberLead.np-1].x - camberLead.Xp[0].x);

  mTrailc = (camberTrail.Xp[camberTrail.np-1].y-camberTrail.Xp[0].y)/
            (camberTrail.Xp[camberTrail.np-1].x-camberTrail.Xp[0].x);

  // Calculate slope of tangent to cross-section at the leading and 
  // trailing edge using a central-difference
  z = Rotor_Spline.Xp[iLead].x;
  dz = (Rotor_Spline.Xp[iLead].x - Rotor_Spline.Xp[iLead+1].x)/TEN;
  mLead = (-   getY(z+2*dz,Rotor_Spline)[0].y 
	   + 8*getY(z+  dz,Rotor_Spline)[0].y
	   - 8*getY(z-  dz,Rotor_Spline)[0].y + 
               getY(z-2*dz,Rotor_Spline)[0].y)/(12*dz);

  z = Rotor_Spline.Xp[iTrail].x;
  dz = (Rotor_Spline.Xp[iTrail].x - Rotor_Spline.Xp[iTrail-1].x)/TEN;
  mTrail = (-   getY(z+2*dz,Rotor_Spline)[1].y
	    + 8*getY(z+  dz,Rotor_Spline)[1].y 
	    - 8*getY(z-  dz,Rotor_Spline)[1].y
	    +   getY(z-2*dz,Rotor_Spline)[1].y)/(12*dz);

  // Get slopes of lines approx 45deg to camber line at the LEADING EDGE.
  A = (mLead+mLeadc)/(1-mLead*mLeadc);
  mTop = (-1+sqrt(1+A*A))/A;
  mBot = (-1-sqrt(1+A*A))/A;

  // Find point on upperMiddle spline - leading edge.
  dm = 1e15;
  mIndex = -1;
  for (int i = 0; i < upperMiddleB.np; i++) {
    if (fabs(upperMiddleB.Xp[i].x-leadV.x) > TOLER) {
      m = (upperMiddleB.Xp[i].y-leadV.y)/(upperMiddleB.Xp[i].x-leadV.x);
      if (fabs(m-mTop) < dm) {
	dm = fabs(m-mTop);
	mIndex = i;
      }
    }
  }
  assert(mIndex != -1);
  zlU = upperMiddleB.Xp[mIndex].x;  

  // Find point on lowerMiddle spline - leading edge.
  dm = 1e15;
  mIndex = -1;
  for (int i = 0; i < lowerMiddleB.np; i++) {
    if (fabs(lowerMiddleB.Xp[i].x-leadV.x) > TOLER) {
      m = (lowerMiddleB.Xp[i].y-leadV.y)/(lowerMiddleB.Xp[i].x-leadV.x);
      if (fabs(m-mBot) < dm) {
	dm=fabs(m-mBot);
	mIndex=i;
      }
    }
  }
  assert(mIndex != -1);
  zlL = lowerMiddleB.Xp[mIndex].x;  

  // Get slopes of lines approx 45deg to camber line at the TRAILING EDGE.
  A = (mTrail+mTrailc)/(1-mTrail*mTrailc);
  mBot = (-1+sqrt(1+A*A))/A;
  mTop = (-1-sqrt(1+A*A))/A;
  
  // Find point on upperMiddle spline - trailing edge.
  dm = 1e15;
  mIndex = -1;
  for (int i=0; i < upperMiddleB.np; i++) {
    if (fabs(upperMiddleB.Xp[i].x-trailV.x) > TOLER) {
      m = (upperMiddleB.Xp[i].y-trailV.y)/(upperMiddleB.Xp[i].x-trailV.x);
      if (fabs(m-mTop) < dm) {
	dm=fabs(m-mTop);
	mIndex=i;
      }
    }
  }
  assert(mIndex != -1);
  zrU = upperMiddleB.Xp[mIndex].x;  

  dm = 1e15;
  mIndex = -1;
  for (int i = 0; i < upperMiddleB.np; i++) {
    if (fabs(upperMiddleB.Xp[i].x-trailV.x) > TOLER) {
      m = (upperMiddleB.Xp[i].y-trailV.y)/(upperMiddleB.Xp[i].x-trailV.x);
      if (fabs(m-mBot) < dm) {
	dm = fabs(m-mBot);
	mIndex = i;
      }
    }
  }
  assert(mIndex != -1);
  
  // Determine which to use.
  if (upperMiddleB.Xp[mIndex].x>zrU) {
    zrU = upperMiddleB.Xp[mIndex].x;  
    swapm = mBot;
    mBot = mTop;
    mTop = swapm;
  }

  // Find point on lowerMiddle spline - trailing edge.
  dm = 1e15;
  mIndex = -1;
  for (int i=0; i < lowerMiddleB.np; i++) {
    if (fabs(lowerMiddleB.Xp[i].x-trailV.x) > TOLER) {
      m = (lowerMiddleB.Xp[i].y-trailV.y)/(lowerMiddleB.Xp[i].x-trailV.x);
      if (fabs(m-mBot) < dm) {
	dm = fabs(m-mBot);
	mIndex = i;
      }
    }
  }
  assert(mIndex != -1);
  zrL = lowerMiddleB.Xp[mIndex].x;

  // NORTH BOUNDARY SPLINE.
  // First segment.
  BndNorthSpline.allocate(30);
  BndNorthSpline.settype(SPLINE2D_QUINTIC);
  zLead = upperB.Xp[0].x;
  zTrail = (zlL+zlU)/2;
  dz = (zTrail-zLead)/29;
  for (int i = 0; i < 30; i++)
    BndNorthSpline.Xp[i] = getY(zLead+i*dz,upperB)[0];
  BndNorthSpline.Xp[0] = upperB.Xp[0];
  BndNorthSpline.pathlength();
  // Second segment.
  tempBndNorthSpline.allocate(70);
  tempBndNorthSpline.settype(SPLINE2D_QUINTIC);
  zLead = (zlL+zlU)/2;
  zTrail = (zrL+zrU)/2;
  dz = (zTrail-zLead)/69;
  for (int i = 0; i < 70; i++)
    tempBndNorthSpline.Xp[i] = getY(zLead+i*dz,upperB)[0];
  tempBndNorthSpline.pathlength();
  // Concatenate first and second segments.
  BndNorthSpline = Concatenate_Splines(BndNorthSpline,tempBndNorthSpline);
  // Final segment.
  tempBndNorthSpline.deallocate();
  tempBndNorthSpline.allocate(30);
  tempBndNorthSpline.settype(SPLINE2D_QUINTIC);
  zLead = (zrL+zrU)/2;
  zTrail = upperB.Xp[upperB.np-1].x;
  dz = (zTrail-zLead)/29;
  for (int i = 0; i < 30; i++) {
    tempBndNorthSpline.Xp[i].x = zLead+i*dz;
    if (tempBndNorthSpline.Xp[i].x > zTrail)
      tempBndNorthSpline.Xp[i].x = zTrail;
    tempBndNorthSpline.Xp[i].y = getY(tempBndNorthSpline.Xp[i].x,upperB)[0].y;
  }
  tempBndNorthSpline.Xp[tempBndNorthSpline.np-1] = upperB.Xp[upperB.np-1];
  tempBndNorthSpline.pathlength();
  // Concatenate first and final segments.
  BndNorthSpline = Concatenate_Splines(BndNorthSpline,tempBndNorthSpline);
  tempBndNorthSpline.deallocate();
  for (int i = 0; i < BndNorthSpline.np; i++) {
    BndNorthSpline.bc[i] = BC_NONE;
    if (i == 0 || i == BndNorthSpline.np-1) {
      BndNorthSpline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      BndNorthSpline.tp[i] = SPLINE2D_POINT_NORMAL;
    }
  }
//    BndNorthSpline.allocate(2);
//    BndNorthSpline.settype(SPLINE2D_LINEAR);
//    BndNorthSpline.Xp[0] = upperB.Xp[0];
//    BndNorthSpline.Xp[1] = upperB.Xp[lowerB.np-1];
//    BndNorthSpline.bc[0] = BC_NONE;
//    BndNorthSpline.bc[1] = BC_NONE;
//    BndNorthSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
//    BndNorthSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
//    BndNorthSpline.pathlength();

  // SOUTH BOUNDARY SPLINE.
  // Second segment.
  tempBndSouthSpline.allocate(70);
  tempBndSouthSpline.settype(SPLINE2D_QUINTIC);
  zLead = (zlL+zlU)/2;;
  zTrail = (zrL+zrU)/2;;
  dz = (zTrail-zLead)/69;
  for (int i = 0; i < tempBndSouthSpline.np; i++)
    tempBndSouthSpline.Xp[i] = getY(zLead+i*dz,lowerB)[0];
  tempBndSouthSpline.pathlength();
  // First segment.
  BndSouthSpline.allocate(30);
  BndSouthSpline.settype(SPLINE2D_QUINTIC);
  zLead = lowerB.Xp[0].x;
  zTrail = tempBndSouthSpline.Xp[0].x;
  dz = (zTrail-zLead)/29;
  for (int i = 0; i < 30; i++)
    BndSouthSpline.Xp[i] = getY(zLead+i*dz,lowerB)[0];
  BndSouthSpline.Xp[0] = lowerB.Xp[0];
  BndSouthSpline.pathlength();
  // Concatenate first and second segments.
  BndSouthSpline = Concatenate_Splines(BndSouthSpline,tempBndSouthSpline);
  tempBndSouthSpline.deallocate();
  // Final segment.
  tempBndSouthSpline.allocate(30);
  tempBndSouthSpline.settype(SPLINE2D_QUINTIC);
  zLead = (zrL+zrU)/2;;
  zTrail = lowerB.Xp[lowerB.np-1].x;
  dz = (zTrail-zLead)/29;
  for (int i = 0; i < 30; i++) {
    tempBndSouthSpline.Xp[i].x = zLead+i*dz;
    if(tempBndSouthSpline.Xp[i].x>zTrail) tempBndSouthSpline.Xp[i].x = zTrail;
    tempBndSouthSpline.Xp[i].y = getY(tempBndSouthSpline.Xp[i].x,lowerB)[0].y;
  }
  tempBndSouthSpline.Xp[tempBndSouthSpline.np-1] = lowerB.Xp[lowerB.np-1];
  tempBndSouthSpline.pathlength();
  // Concatenate first and final segments.
  BndSouthSpline = Concatenate_Splines(BndSouthSpline,tempBndSouthSpline);
  tempBndSouthSpline.deallocate();
  for (int i = 0; i < BndSouthSpline.np; i++) {
    BndSouthSpline.bc[i] = BC_NONE;
    if (i == 0 || i == BndSouthSpline.np-1) {
      BndSouthSpline.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
    } else {
      BndSouthSpline.tp[i] = SPLINE2D_POINT_NORMAL;
    }
  }
//    BndSouthSpline.allocate(2);
//    BndSouthSpline.settype(SPLINE2D_LINEAR);
//    BndSouthSpline.Xp[0] = lowerB.Xp[0];
//    BndSouthSpline.Xp[1] = lowerB.Xp[lowerB.np-1];
//    BndSouthSpline.bc[0] = BC_NONE;
//    BndSouthSpline.bc[1] = BC_NONE;
//    BndSouthSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
//    BndSouthSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
//    BndSouthSpline.pathlength();

  // EAST BOUNDARY SPLINE.
  BndEastSpline.allocate(2);
  BndEastSpline.settype(SPLINE2D_LINEAR);
  BndEastSpline.Xp[0] = BndSouthSpline.Xp[BndSouthSpline.np-1];
  BndEastSpline.Xp[1] = BndNorthSpline.Xp[BndNorthSpline.np-1];
  BndEastSpline.bc[0] = BC_CHARACTERISTIC;
  BndEastSpline.bc[1] = BC_CHARACTERISTIC;
  BndEastSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  BndEastSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  BndEastSpline.pathlength();

  // WEST BOUNDARY SPLINE.
  BndWestSpline.allocate(2);
  BndWestSpline.settype(SPLINE2D_LINEAR);
  BndWestSpline.Xp[0] = BndSouthSpline.Xp[0];
  BndWestSpline.Xp[1] = BndNorthSpline.Xp[0];
  BndWestSpline.bc[0] = BC_FIXED;
  BndWestSpline.bc[1] = BC_FIXED;
  BndWestSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
  BndWestSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
  BndWestSpline.pathlength();

  // Assign values to the stretching function parameters and
  // boundary grid line orthogonality parameters.
  Stretch_I = STRETCHING_FCN_LINEAR;
  Beta_I = ZERO; 
  Tau_I = ZERO;
  Stretch_J = STRETCHING_FCN_LINEAR;
  Beta_J = ZERO;
  Tau_J = ZERO;
  Orthogonal_North = 0;
  Orthogonal_South = 0;
  Orthogonal_East = 0;
  Orthogonal_West = 0;

  // Create the mesh for each block representing the complete grid.
  for (int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++) {
    for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {
      // Create the splines defining the north, south, east, and west 
      // boundaries of the block.
      xSW.x = BndWestSpline.Xp[0].x + ((double(iBlk)/double(Number_of_Blocks_Idir))*
				       (BndEastSpline.Xp[0].x - BndWestSpline.Xp[0].x));
      xSW.y = getY(xSW.x,BndSouthSpline)[0].y + ((double(jBlk)/double(Number_of_Blocks_Jdir))*
						 (getY(xSW.x,BndNorthSpline)[0].y - getY(xSW.x,BndSouthSpline)[0].y));
      xSE.x = BndWestSpline.Xp[0].x + ((double(iBlk+1)/double(Number_of_Blocks_Idir))*
				       (BndEastSpline.Xp[0].x - BndWestSpline.Xp[0].x));
      xSE.y = getY(xSE.x,BndSouthSpline)[0].y + ((double(jBlk)/double(Number_of_Blocks_Jdir))*
						 (getY(xSE.x,BndNorthSpline)[0].y - getY(xSE.x,BndSouthSpline)[0].y));
      xNW.x = xSW.x;
      xNW.y = getY(xSW.x,BndSouthSpline)[0].y + ((double(jBlk+1)/double(Number_of_Blocks_Jdir))*
						 (getY(xNW.x,BndNorthSpline)[0].y - getY(xNW.x,BndSouthSpline)[0].y));
      xNE.x = BndWestSpline.Xp[0].x + ((double(iBlk+1)/double(Number_of_Blocks_Idir))*
				       (BndEastSpline.Xp[0].x - BndWestSpline.Xp[0].x));
      xNE.y = getY(xNE.x,BndSouthSpline)[0].y + ((double(jBlk+1)/double(Number_of_Blocks_Jdir))*
						 (getY(xNE.x,BndNorthSpline)[0].y - getY(xNE.x,BndSouthSpline)[0].y));
      // North Spline.
      Bnd_Spline_North.allocate(10);
      Bnd_Spline_North.settype(SPLINE2D_QUINTIC);
      for (int i = 0; i < Bnd_Spline_North.np; i++) {
	Bnd_Spline_North.Xp[i].x = xNW.x + (double(i)/double(Bnd_Spline_North.np-1))*(xNE.x - xNW.x);
	Bnd_Spline_North.Xp[i].y = getY(Bnd_Spline_North.Xp[i].x,BndSouthSpline)[0].y + 
                                   (double(jBlk+1)/double(Number_of_Blocks_Jdir))*
                                   (getY(Bnd_Spline_North.Xp[i].x,BndNorthSpline)[0].y - 
				    getY(Bnd_Spline_North.Xp[i].x,BndSouthSpline)[0].y);
  	Bnd_Spline_North.bc[i] = BC_NONE;
	if (i == 0 || i == Bnd_Spline_North.np-1) {
	  Bnd_Spline_North.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
	} else {
	  Bnd_Spline_North.tp[i] = SPLINE2D_POINT_NORMAL;
	}
      }
      Bnd_Spline_North.pathlength();
      // South Spline.
      Bnd_Spline_South.allocate(10);
      Bnd_Spline_South.settype(SPLINE2D_QUINTIC);
      for (int i = 0; i < Bnd_Spline_South.np; i++) {
	Bnd_Spline_South.Xp[i].x = xSW.x + (double(i)/double(Bnd_Spline_South.np-1))*(xSE.x - xSW.x);
	Bnd_Spline_South.Xp[i].y = getY(Bnd_Spline_South.Xp[i].x,BndSouthSpline)[0].y + 
                                   (double(jBlk)/double(Number_of_Blocks_Jdir))*
	                           (getY(Bnd_Spline_South.Xp[i].x,BndNorthSpline)[0].y - 
				    getY(Bnd_Spline_South.Xp[i].x,BndSouthSpline)[0].y);
	Bnd_Spline_South.bc[i] = BC_NONE;
	if (i == 0 || i == Bnd_Spline_South.np-1) {
	  Bnd_Spline_South.tp[i] = SPLINE2D_POINT_SHARP_CORNER;
	} else {
	  Bnd_Spline_South.tp[i] = SPLINE2D_POINT_NORMAL;
	}
      }
      Bnd_Spline_South.pathlength();
      // East Spline.
      Bnd_Spline_East.allocate(2);
      Bnd_Spline_East.settype(SPLINE2D_LINEAR);
      Bnd_Spline_East.Xp[0] = xSE;
      Bnd_Spline_East.Xp[1] = xNE;
      if (iBlk == Number_of_Blocks_Idir-1) {
	Bnd_Spline_East.bc[0] = BC_CHARACTERISTIC;
	Bnd_Spline_East.bc[1] = BC_CHARACTERISTIC;
      } else {
	Bnd_Spline_East.bc[0] = BC_NONE;
	Bnd_Spline_East.bc[1] = BC_NONE;
      }
      Bnd_Spline_East.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
      Bnd_Spline_East.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
      Bnd_Spline_East.pathlength();
      // West Spline.
      Bnd_Spline_West.allocate(2);
      Bnd_Spline_West.settype(SPLINE2D_LINEAR);
      Bnd_Spline_West.Xp[0] = xSW;
      Bnd_Spline_West.Xp[1] = xNW;
      if (iBlk == 0) {
	Bnd_Spline_West.bc[0] = BC_FIXED;
	Bnd_Spline_West.bc[1] = BC_FIXED;
      } else {
	Bnd_Spline_West.bc[0] = BC_NONE;
	Bnd_Spline_West.bc[1] = BC_NONE;
      }
      Bnd_Spline_West.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
      Bnd_Spline_West.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
      Bnd_Spline_West.pathlength();
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
      Smooth_Quad_Block(Grid_ptr[iBlk][jBlk],
			min(250,2*max(Number_of_Cells_Idir,
				      Number_of_Cells_Jdir)));
      // Deallocate the memory for the boundary splines.
      Bnd_Spline_North.deallocate();
      Bnd_Spline_South.deallocate();
      Bnd_Spline_East.deallocate();
      Bnd_Spline_West.deallocate();
    }
  }

  // Deallocate splines.
  upperB.deallocate();
  lowerB.deallocate();
  upperMiddleB.deallocate();
  lowerMiddleB.deallocate();
  camberTrail.deallocate();
  camberLead.deallocate();
  camberBlade.deallocate();
  Rotor_Spline.deallocate();
  BndNorthSpline.deallocate();
  BndSouthSpline.deallocate();
  BndEastSpline.deallocate();
  BndWestSpline.deallocate();

  // Return the grid.
  return Grid_ptr;

}

/**********************************************************************
 * Routine: Grid_Driven_Cavity_Flow                                   *
 *                                                                    *
 * Generates a mesh for the driven cavity flow.                       *
 *                                                                    *
 * Usage: Grid_ptr = Grid_Driven_Cavity(Grid_ptr,                     *
 *                                      nblk_i,                       *
 *                                      nblk_j,                       *
 *                                      TEN,                          *
 *                                      FIVE,                         *
 *   	                                100,                          *
 *  	                                50,                           *
 *                                      2);                           *
 *                                                                    *
 **********************************************************************/
Grid2D_Quad_Block** Grid_Driven_Cavity_Flow(Grid2D_Quad_Block **Grid_ptr,
					    int &Number_of_Blocks_Idir,
					    int &Number_of_Blocks_Jdir,
					    const double &Width,
					    const double &Height,
					    const int &Stretching_Type_Idir,
					    const int &Stretching_Type_Jdir,
					    const double &Stretching_Factor_Idir,
					    const double &Stretching_Factor_Jdir,
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

  // Allocate memory for grid block.

  Number_of_Blocks_Idir = 1;
  Number_of_Blocks_Jdir = 1;
  Grid_ptr = Allocate_Multi_Block_Grid(Grid_ptr, 
				       Number_of_Blocks_Idir, 
				       Number_of_Blocks_Jdir);

  // Create the mesh for each block representing the complete grid.

  for (int jBlk = 0; jBlk < Number_of_Blocks_Jdir; jBlk++) {
    for (int iBlk = 0; iBlk < Number_of_Blocks_Idir; iBlk++) {

      // Assign values to the locations of the corners of the
      // rectangular box shaped domain.

      xc_NW = Vector2D(-HALF*Width+double(iBlk)*Width/double(Number_of_Blocks_Idir), 
		       -HALF*Height+double(jBlk+1)*Height/double(Number_of_Blocks_Jdir));
      xc_NE = Vector2D(-HALF*Width+double(iBlk+1)*Width/double(Number_of_Blocks_Idir), 
		       -HALF*Height+double(jBlk+1)*Height/double(Number_of_Blocks_Jdir));
      xc_SE = Vector2D(-HALF*Width+double(iBlk+1)*Width/double(Number_of_Blocks_Idir), 
		       -HALF*Height+double(jBlk)*Height/double(Number_of_Blocks_Jdir));
      xc_SW = Vector2D(-HALF*Width+double(iBlk)*Width/double(Number_of_Blocks_Idir), 
		       -HALF*Height+double(jBlk)*Height/double(Number_of_Blocks_Jdir));

      // Create the splines defining the north, south, east, and west
      // boundaries of the rectangular box.
      Create_Spline_Line(Bnd_Spline_North, xc_NW, xc_NE, 2);
      Create_Spline_Line(Bnd_Spline_South, xc_SW, xc_SE, 2);
      Create_Spline_Line(Bnd_Spline_East, xc_SE, xc_NE, 2);
      Create_Spline_Line(Bnd_Spline_West, xc_SW, xc_NW, 2);

      // Set the boundary condition types for each of the boundary
      // splines.
      if (jBlk == Number_of_Blocks_Jdir-1) {
	Bnd_Spline_North.setBCtype(BC_MOVING_WALL_ISOTHERMAL);
      } else {
	Bnd_Spline_North.setBCtype(BC_NONE);
      }
      if (jBlk == 0) {
	Bnd_Spline_South.setBCtype(BC_WALL_VISCOUS_ISOTHERMAL);
      } else {
	Bnd_Spline_South.setBCtype(BC_NONE);
      }
      if (iBlk == Number_of_Blocks_Idir-1) {
	Bnd_Spline_East.setBCtype(BC_WALL_VISCOUS_ISOTHERMAL);
      } else {
	Bnd_Spline_East.setBCtype(BC_NONE);
      }
      if (iBlk == 0) {
	Bnd_Spline_West.setBCtype(BC_WALL_VISCOUS_ISOTHERMAL);
      } else {
	Bnd_Spline_West.setBCtype(BC_NONE);
      }

      // Assign values to the stretching function parameters and
      // boundary grid line orthogonality parameters.
      Stretch_I = Stretching_Type_Idir;
      Stretch_J = Stretching_Type_Jdir;
      Beta_I = Stretching_Factor_Idir; Tau_I = ZERO;
      Beta_J = Stretching_Factor_Jdir; Tau_J = ZERO;
      Orthogonal_North = 0;
      Orthogonal_South = 0;
      Orthogonal_East = 0;
      Orthogonal_West = 0;

      // Create the 2D quadrilateral grid block representing the mesh.
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
