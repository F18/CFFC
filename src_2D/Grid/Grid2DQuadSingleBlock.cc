/* Grid2DQuadSingleBlock.cc:  Single-block subroutines for 
                              2D quadrilateral block grid class. */

/* Include 2D quadrilateral block grid type header file. */

#ifndef _GRID2D_QUAD_BLOCK_INCLUDED
#include "Grid2DQuad.h"
#endif // _GRID2D_QUAD_BLOCK_INCLUDED

/*************************************************************************
 * Grid2D_Quad_Block -- External subroutines for single grid block.      *
 *************************************************************************/

/********************************************************
 * Routine: Create_Quad_Block                           *
 *                                                      *
 * Create quadrilateral grid block for a Cartesian      *
 * mesh defined on a square box with four corner        *
 * coordinates (X,Y): (-1,1), (1,1), (1,-1), (-1,1).    *
 *                                                      *
 ********************************************************/
void Create_Quad_Block(Grid2D_Quad_Block &Grid,
                       const int Number_of_Cells_Idir,
  	               const int Number_of_Cells_Jdir,
		       const int Number_of_Ghost_Cells) {

    int i, j;
    double S_i, S_j, 
           s_north, s_south, s_east, s_west,
           smax_north, smax_south, smax_east, smax_west, 
           s_i, s_j, smax_i, smax_j,
           w_north, w_south, w_east, w_west, w_total;
    Vector2D x_north, x_south, x_east, x_west;
 
    /* Allocate (re-allocate) memory for the cells and nodes 
       of the quadrilateral mesh block. */

    if (Grid.Node != NULL && Grid.Cell != NULL) { 
       Grid.deallocate();
    } else if (Grid.Node != NULL) {
       Grid.deallocateNodes();
    } else if (Grid.Cell != NULL) {
       Grid.deallocateCells();
    } /* endif */
    Grid.allocate(Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Ghost_Cells);

    /* Allocate (re-allocate) memory for the boundary splines 
       defining this quadrilateral mesh block. */

    if (Grid.BndNorthSpline.np != 0) Grid.BndNorthSpline.deallocate();
    Grid.BndNorthSpline.allocate(2);

    if (Grid.BndSouthSpline.np != 0) Grid.BndSouthSpline.deallocate();
    Grid.BndSouthSpline.allocate(2);

    if (Grid.BndEastSpline.np != 0) Grid.BndEastSpline.deallocate();
    Grid.BndEastSpline.allocate(2);

    if (Grid.BndWestSpline.np != 0) Grid.BndWestSpline.deallocate();
    Grid.BndWestSpline.allocate(2);

    /* Set the boundary spline types to linear, assign
       spline points, and calculate the spline pathlengths. */

    Grid.BndNorthSpline.settype(SPLINE2D_LINEAR);
    Grid.BndNorthSpline.Xp[0] = Vector2D(-ONE, ONE);
    Grid.BndNorthSpline.Xp[1] = Vector2D( ONE, ONE);
    Grid.BndNorthSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
    Grid.BndNorthSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
    Grid.BndNorthSpline.bc[0] = BC_REFLECTION;
    Grid.BndNorthSpline.bc[1] = BC_REFLECTION;
    Grid.BndNorthSpline.pathlength();
    Grid.SminN = Grid.BndNorthSpline.sp[0];
    Grid.SmaxN = Grid.BndNorthSpline.sp[Grid.BndNorthSpline.np-1];

    Grid.BndSouthSpline.settype(SPLINE2D_LINEAR);
    Grid.BndSouthSpline.Xp[0] = Vector2D(-ONE,-ONE);
    Grid.BndSouthSpline.Xp[1] = Vector2D( ONE,-ONE);
    Grid.BndSouthSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
    Grid.BndSouthSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
    Grid.BndSouthSpline.bc[0] = BC_REFLECTION;
    Grid.BndSouthSpline.bc[1] = BC_REFLECTION;
    Grid.BndSouthSpline.pathlength();
    Grid.SminS = Grid.BndSouthSpline.sp[0];
    Grid.SmaxS = Grid.BndSouthSpline.sp[Grid.BndSouthSpline.np-1];

    Grid.BndEastSpline.settype(SPLINE2D_LINEAR);
    Grid.BndEastSpline.Xp[0] = Vector2D( ONE,-ONE);
    Grid.BndEastSpline.Xp[1] = Vector2D( ONE, ONE);
    Grid.BndEastSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
    Grid.BndEastSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
    Grid.BndEastSpline.bc[0] = BC_REFLECTION;
    Grid.BndEastSpline.bc[1] = BC_REFLECTION;
    Grid.BndEastSpline.pathlength();
    Grid.SminE = Grid.BndEastSpline.sp[0];
    Grid.SmaxE = Grid.BndEastSpline.sp[Grid.BndEastSpline.np-1];

    Grid.BndWestSpline.settype(SPLINE2D_LINEAR);
    Grid.BndWestSpline.Xp[0] = Vector2D(-ONE,-ONE);
    Grid.BndWestSpline.Xp[1] = Vector2D(-ONE, ONE);
    Grid.BndWestSpline.tp[0] = SPLINE2D_POINT_SHARP_CORNER;
    Grid.BndWestSpline.tp[1] = SPLINE2D_POINT_SHARP_CORNER;
    Grid.BndWestSpline.bc[0] = BC_REFLECTION;
    Grid.BndWestSpline.bc[1] = BC_REFLECTION;
    Grid.BndWestSpline.pathlength();
    Grid.SminW = Grid.BndWestSpline.sp[0];
    Grid.SmaxW = Grid.BndWestSpline.sp[Grid.BndWestSpline.np-1];

    /* Compute the interior nodes for the quadrilateral mesh block. */

    for ( j = Grid.JNl ; j <= Grid.JNu ; ++j) {
        S_j = double(j-Grid.JNl)/double(Grid.JNu-Grid.JNl);

        smax_east = Grid.BndEastSpline.sp[Grid.BndEastSpline.np-1]-
                    Grid.BndEastSpline.sp[0];
        s_east = S_j * smax_east + Grid.BndEastSpline.sp[0];
        x_east = Spline(s_east, Grid.BndEastSpline);
	
        smax_west = Grid.BndWestSpline.sp[Grid.BndWestSpline.np-1]-
                    Grid.BndWestSpline.sp[0];
        s_west = S_j * smax_west + Grid.BndWestSpline.sp[0];
        x_west = Spline(s_west, Grid.BndWestSpline);

        for ( i = Grid.INl ; i <= Grid.INu ; ++i) {
            S_i = double(i-Grid.INl)/double(Grid.INu-Grid.INl);

            smax_north = Grid.BndNorthSpline.sp[Grid.BndNorthSpline.np-1]-
	                 Grid.BndNorthSpline.sp[0];
            s_north = S_i * smax_north + Grid.BndNorthSpline.sp[0];
            x_north = Spline(s_north, Grid.BndNorthSpline);

            smax_south = Grid.BndSouthSpline.sp[Grid.BndSouthSpline.np-1]-
                         Grid.BndSouthSpline.sp[0];
            s_south = S_i * smax_south + Grid.BndSouthSpline.sp[0];
            x_south = Spline(s_south, Grid.BndSouthSpline);

            if (j == Grid.JNl) {
	      Grid.Node[i][j].X = x_south;
            } else if (j == Grid.JNu) {
	      Grid.Node[i][j].X = x_north;
            } else if (i == Grid.INl) {
	      Grid.Node[i][j].X = x_west;
            } else if (i == Grid.INu) {
	      Grid.Node[i][j].X = x_east;
            } else {
              s_i = (ONE-S_j)*s_south + S_j*s_north;
              s_j = (ONE-S_i)*s_west + S_i*s_east;

              smax_i = (ONE-S_j)*smax_south + S_j*smax_north;
              smax_j = (ONE-S_i)*smax_west + S_i*smax_east;
              
              w_south = (ONE-s_j/smax_j);
              w_north = (s_j/smax_j);
              w_total = w_south + w_north;
              Grid.Node[i][j].X = (w_south*x_south + 
                                   w_north*x_north)/w_total;
            } /* endif */

        } /* endfor */
    } /* endfor */

    /* Set the boundary condition types at the grid block
       boundaries. */
 
    Set_BCs(Grid);

    /* Compute the exterior nodes for the quadrilateral mesh block. */

    Update_Exterior_Nodes(Grid);

    /* Compute the cells for the quadrilateral mesh block. */

    Update_Cells(Grid);

}

/********************************************************
 * Routine: Create_Quad_Block                           *
 *                                                      *
 * Create a 2D quadrilateral grid block with the four   *
 * boundaries of the block defined by blended splines   *
 * which are given in the input boundary spline file:   *
 *                                                      *
 * Bnd_Spline_File_Name_ptr.                            *
 *                                                      *
 ********************************************************/
void Create_Quad_Block(Grid2D_Quad_Block &Grid,
                       char *Bnd_Spline_File_Name_ptr,
                       const int Number_of_Cells_Idir,
  	               const int Number_of_Cells_Jdir,
		       const int Number_of_Ghost_Cells) {

    ifstream bnd_spline_file;
    int i, j, k, kx, ky, kx_max, ky_max, npts, spline_type;
    int node_init_procedure;
    double S_i, S_j, 
           s_north, s_south, s_east, s_west,
           smax_north, smax_south, smax_east, smax_west, 
           s_i, s_j, smax_i, smax_j,
           w_north, w_south, w_east, w_west, w_total;
    double smax_S1_NS, smax_S2_NS, smax_S1_EW, smax_S2_EW,
           S1_i, S2_i, S1_j, S2_j, dS_i, dS_j;
    Vector2D x_north, x_south, x_east, x_west;
    Spline2D S1_NS, S2_NS, S1_EW, S2_EW;
 
    /* Allocate (re-allocate) memory for the cells and nodes 
       of the quadrilateral mesh block. */

    if (Grid.Node != NULL && Grid.Cell != NULL) { 
       Grid.deallocate();
    } else if (Grid.Node != NULL) {
       Grid.deallocateNodes();
    } else if (Grid.Cell != NULL) {
       Grid.deallocateCells();
    } /* endif */
    Grid.allocate(Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Ghost_Cells);

    /* Open the data file containing the boundary splines. */

    bnd_spline_file.open(Bnd_Spline_File_Name_ptr,ios::in);

    /* For each of the north, south, east, and west boundaries
       of this mesh block, read in the number of spline points, 
       allocate memory for the boundary splines, read in the 
       spline points, and finally calculate the spline 
       pathlengths. */

    bnd_spline_file.setf(ios::skipws);
    bnd_spline_file >> npts;
    bnd_spline_file.unsetf(ios::skipws);
    if (Grid.BndNorthSpline.np != 0) Grid.BndNorthSpline.deallocate();
    Grid.BndNorthSpline.allocate(npts);
    bnd_spline_file.setf(ios::skipws);
    bnd_spline_file >> spline_type;
    bnd_spline_file.unsetf(ios::skipws);
    Grid.BndNorthSpline.settype(spline_type);
    bnd_spline_file >> Grid.BndNorthSpline;
    Grid.BndNorthSpline.pathlength();
    Grid.SminN = Grid.BndNorthSpline.sp[0];
    Grid.SmaxN = Grid.BndNorthSpline.sp[Grid.BndNorthSpline.np-1];

    bnd_spline_file.setf(ios::skipws);
    bnd_spline_file >> npts;
    bnd_spline_file.unsetf(ios::skipws);
    if (Grid.BndSouthSpline.np != 0) Grid.BndSouthSpline.deallocate();
    Grid.BndSouthSpline.allocate(npts);
    bnd_spline_file.setf(ios::skipws);
    bnd_spline_file >> spline_type;
    bnd_spline_file.unsetf(ios::skipws);
    Grid.BndSouthSpline.settype(spline_type);
    bnd_spline_file >> Grid.BndSouthSpline;
    Grid.BndSouthSpline.pathlength();
    Grid.SminS = Grid.BndSouthSpline.sp[0];
    Grid.SmaxS = Grid.BndSouthSpline.sp[Grid.BndSouthSpline.np-1];

    bnd_spline_file.setf(ios::skipws);
    bnd_spline_file >> npts;
    bnd_spline_file.unsetf(ios::skipws);
    if (Grid.BndEastSpline.np != 0) Grid.BndEastSpline.deallocate();
    Grid.BndEastSpline.allocate(npts);
    bnd_spline_file.setf(ios::skipws);
    bnd_spline_file >> spline_type;
    bnd_spline_file.unsetf(ios::skipws);
    Grid.BndEastSpline.settype(spline_type);
    bnd_spline_file >> Grid.BndEastSpline;
    Grid.BndEastSpline.pathlength();
    Grid.SminE = Grid.BndEastSpline.sp[0];
    Grid.SmaxE = Grid.BndEastSpline.sp[Grid.BndEastSpline.np-1];

    bnd_spline_file.setf(ios::skipws);
    bnd_spline_file >> npts;
    bnd_spline_file.unsetf(ios::skipws);
    if (Grid.BndWestSpline.np != 0) Grid.BndWestSpline.deallocate();
    Grid.BndWestSpline.allocate(npts);
    bnd_spline_file.setf(ios::skipws);
    bnd_spline_file >> spline_type;
    bnd_spline_file.unsetf(ios::skipws);
    Grid.BndWestSpline.settype(spline_type);
    bnd_spline_file >> Grid.BndWestSpline;
    Grid.BndWestSpline.pathlength();
    Grid.SminW = Grid.BndWestSpline.sp[0];
    Grid.SmaxW = Grid.BndWestSpline.sp[Grid.BndWestSpline.np-1];

    /* Read the node initialization procedure for this 
       quadrilateral grid block. */

    bnd_spline_file.setf(ios::skipws);
    bnd_spline_file >> node_init_procedure;
    bnd_spline_file.unsetf(ios::skipws);

    /* Input the type of node distribution stretching 
       functions to be used for each of the coordinate
       directions. Also read in the related stretching
       parameters. */

    bnd_spline_file.setf(ios::skipws);
    bnd_spline_file >> Grid.StretchI >> Grid.BetaI >> Grid.TauI;
    bnd_spline_file >> Grid.StretchJ >> Grid.BetaJ >> Grid.TauJ;
    bnd_spline_file.unsetf(ios::skipws);

    /* Input the grid orthogonality specifiers for
       the north, south, east, and west boundaries. */

    bnd_spline_file.setf(ios::skipws);
    bnd_spline_file >> Grid.OrthogonalN >> Grid.OrthogonalS
                    >> Grid.OrthogonalE >> Grid.OrthogonalW;
    bnd_spline_file.unsetf(ios::skipws);

    /* Close the data file containing the boundary splines. */

    bnd_spline_file.close();

    /* Compute the interior nodes for the quadrilateral mesh block. */

    smax_north = Grid.BndNorthSpline.sp[Grid.BndNorthSpline.np-1]-
                 Grid.BndNorthSpline.sp[0];
    smax_south = Grid.BndSouthSpline.sp[Grid.BndSouthSpline.np-1]-
                 Grid.BndSouthSpline.sp[0];

    smax_east = Grid.BndEastSpline.sp[Grid.BndEastSpline.np-1]-
                Grid.BndEastSpline.sp[0];
    smax_west = Grid.BndWestSpline.sp[Grid.BndWestSpline.np-1]-
                Grid.BndWestSpline.sp[0];

    dS_i = ONE/(TEN*double(Grid.INu-Grid.INl));
    dS_j = ONE/(TEN*double(Grid.JNu-Grid.JNl));

    for ( j = Grid.JNl ; j <= Grid.JNu ; ++j) {
        for ( i = Grid.INl ; i <= Grid.INu ; ++i) {
            S_j = double(j-Grid.JNl)/double(Grid.JNu-Grid.JNl);
            S_j = StretchingFcn(S_j, Grid.BetaJ, Grid.TauJ, Grid.StretchJ);
            s_east = S_j * smax_east + Grid.BndEastSpline.sp[0];
            x_east = Spline(s_east, Grid.BndEastSpline);
            s_west = S_j * smax_west + Grid.BndWestSpline.sp[0];
            x_west = Spline(s_west, Grid.BndWestSpline);

            S_i = double(i-Grid.INl)/double(Grid.INu-Grid.INl);
            S_i = StretchingFcn(S_i, Grid.BetaI, Grid.TauI, Grid.StretchI);
            s_north = S_i * smax_north + Grid.BndNorthSpline.sp[0];
            x_north = Spline(s_north, Grid.BndNorthSpline);
            s_south = S_i * smax_south + Grid.BndSouthSpline.sp[0];
            x_south = Spline(s_south, Grid.BndSouthSpline);

            if (j == Grid.JNl) {
	      Grid.Node[i][j].X = x_south;
            } else if (j == Grid.JNu) {
	      Grid.Node[i][j].X = x_north;
            } else if (i == Grid.INl) {
	      Grid.Node[i][j].X = x_west;
            } else if (i == Grid.INu) {
	      Grid.Node[i][j].X = x_east;
            } else {
              s_i = (ONE-S_j)*s_south + S_j*s_north;
              s_j = (ONE-S_i)*s_west + S_i*s_east;

              smax_i = (ONE-S_j)*smax_south + S_j*smax_north;
              smax_j = (ONE-S_i)*smax_west + S_i*smax_east;

              switch(node_init_procedure) {
                //============================================================
                case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH :
                  w_south = (ONE-s_j/smax_j);
                  w_north = (s_j/smax_j);
                  w_total = w_south + w_north;
                  Grid.Node[i][j].X = (w_south*x_south + 
                                       w_north*x_north)/w_total;
                  break;
                //============================================================
                case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST :
                  w_west =  (ONE-s_i/smax_i);
                  w_east =  (s_i/smax_i);
                  w_total = w_east + w_west;
                  Grid.Node[i][j].X = (w_west*x_west + 
                                       w_east*x_east)/w_total;
                  break;
                //============================================================
                case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_XY :
                  kx_max = 5;
                  dS_i = ONE/double(kx_max-1);
                  kx = int(floor(S_i/dS_i));
                  S1_i = max(ZERO, double(kx)*dS_i);
                  S2_i = min(ONE, S1_i + dS_i);

                  ky_max = 7;
                  dS_j = ONE/double(ky_max-1);
                  ky = int(floor(S_j/dS_j));
                  S1_j = max(ZERO, double(ky)*dS_j);
                  S2_j = min(ONE, S1_j + dS_j);

		  if (S1_i - TOLER > ZERO) {
                     S1_NS.allocate(2);
                     S1_NS.settype(SPLINE2D_LINEAR);

                     s_south = S1_i * smax_south + Grid.BndSouthSpline.sp[0];
                     S1_NS.Xp[0] = Spline(s_south, Grid.BndSouthSpline);
                     S1_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                     s_north = S1_i * smax_north + Grid.BndNorthSpline.sp[0],
                     S1_NS.Xp[1] = Spline(s_north, Grid.BndNorthSpline);
                     S1_NS.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

                     S1_NS.pathlength();
                     smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
                                  S1_NS.sp[0];
                  } else {
                     Copy_Spline(S1_NS, Grid.BndWestSpline);
                     smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
                                  S1_NS.sp[0];
                  } /* endif */

                  if (S2_i + TOLER < ONE) {
                     S2_NS.allocate(2);
                     S2_NS.settype(SPLINE2D_LINEAR);

                     s_south = S2_i * smax_south + Grid.BndSouthSpline.sp[0];
                     S2_NS.Xp[0] = Spline(s_south, Grid.BndSouthSpline);
                     S2_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                     s_north = S2_i * smax_north + Grid.BndNorthSpline.sp[0],
                     S2_NS.Xp[1] = Spline(s_north, Grid.BndNorthSpline);
                     S2_NS.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

                     S2_NS.pathlength();
                     smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
                                  S2_NS.sp[0];
                  } else {
                     Copy_Spline(S2_NS, Grid.BndEastSpline);
                     smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
                                  S2_NS.sp[0];
                  } /* endif */

		  if (S1_j - TOLER > ZERO) {
                     S1_EW.allocate(kx_max);
                     S1_EW.settype(SPLINE2D_LINEAR);

                     s_west = S1_j * smax_west + Grid.BndWestSpline.sp[0];
                     S1_EW.Xp[0] = Spline(s_west, Grid.BndWestSpline);
                     S1_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                     s_east = S1_j * smax_east + Grid.BndEastSpline.sp[0];
                     S1_EW.Xp[kx_max-1] = Spline(s_east, Grid.BndEastSpline);
                     S1_EW.tp[kx_max-1] = SPLINE2D_POINT_SHARP_CORNER;

                     for ( k = 1 ; k <= kx_max-2 ; ++k) {                     
                        s_north = double(k)*smax_north/double(kx_max-1) + Grid.BndNorthSpline.sp[0];
                        x_north = Spline(s_north, Grid.BndNorthSpline);

                        s_south = double(k)*smax_south/double(kx_max-1) + Grid.BndSouthSpline.sp[0];
                        x_south = Spline(s_south, Grid.BndSouthSpline);

                        w_south =  ONE-S1_j;
                        w_north =  S1_j;
                        w_total = w_north + w_south;
                        S1_EW.Xp[k] = (w_south*x_south + w_north*x_north)/w_total;
                        S1_EW.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
                     } /* endfor */

                     S1_EW.pathlength();
                     smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
                                  S1_EW.sp[0];
                  } else {
                     Copy_Spline(S1_EW, Grid.BndSouthSpline);
                     smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
                                  S1_EW.sp[0];
                  } /* endif */

                  if (S2_j + TOLER < ONE) {
                     S2_EW.allocate(kx_max);
                     S2_EW.settype(SPLINE2D_LINEAR);

                     s_west = S2_j * smax_west + Grid.BndWestSpline.sp[0];
                     S2_EW.Xp[0] = Spline(s_west, Grid.BndWestSpline);
                     S2_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                     s_east = S2_j * smax_east + Grid.BndEastSpline.sp[0];
                     S2_EW.Xp[kx_max-1] = Spline(s_east, Grid.BndEastSpline);
                     S2_EW.tp[kx_max-1] = SPLINE2D_POINT_SHARP_CORNER;

                     for ( k = 1 ; k <= kx_max-2 ; ++k) {                     
                        s_north = double(k)*smax_north/double(kx_max-1) + Grid.BndNorthSpline.sp[0];
                        x_north = Spline(s_north, Grid.BndNorthSpline);

                        s_south = double(k)*smax_south/double(kx_max-1) + Grid.BndSouthSpline.sp[0];
                        x_south = Spline(s_south, Grid.BndSouthSpline);

                        w_south =  ONE-S2_j;
                        w_north =  S2_j;
                        w_total = w_north + w_south;
                        S2_EW.Xp[k] = (w_south*x_south + w_north*x_north)/w_total;
                        S2_EW.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
                     } /* endfor */

                     S2_EW.pathlength();
                     smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
                                  S2_EW.sp[0];
                  } else {
                     Copy_Spline(S2_EW, Grid.BndNorthSpline);
                     smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
                                  S2_EW.sp[0];
                  } /* endif */

                  s_north = S_i * smax_S2_EW + S2_EW.sp[0];
                  x_north = Spline(s_north, S2_EW);
                  s_south = S_i * smax_S1_EW + S1_EW.sp[0];
                  x_south = Spline(s_south, S1_EW);

                  s_east = S_j * smax_S2_NS + S2_NS.sp[0];
                  x_east = Spline(s_east, S2_NS);
                  s_west = S_j * smax_S1_NS + S1_NS.sp[0];
                  x_west = Spline(s_west, S1_NS);
      
                  if (ky == 0 || ky == ky_max - 2) {
                     w_south =  ONE-(S_j-S1_j)/(S2_j-S1_j);
                     w_north =  (S_j-S1_j)/(S2_j-S1_j);
                     w_total = w_south + w_north;
                     Grid.Node[i][j].X = (w_south*x_south + 
                                       w_north*x_north)/w_total;
                  } else {
                     w_west =  ONE-(S_i-S1_i)/(S2_i-S1_i);
                     w_east =  (S_i-S1_i)/(S2_i-S1_i);
                     w_total = w_east + w_west;
                     Grid.Node[i][j].X = (w_west*x_west + 
                                          w_east*x_east)/w_total;
                  } /* endif */

                  S1_NS.deallocate(); S2_NS.deallocate();
                  S2_EW.deallocate(); S2_EW.deallocate();
                  break;
                //============================================================
                case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_YX :
                  kx_max = 7;
                  dS_i = ONE/double(kx_max-1);
                  kx = int(floor(S_i/dS_i));
                  S1_i = max(ZERO, double(kx)*dS_i);
                  S2_i = min(ONE, S1_i + dS_i);

                  ky_max = 5;
                  dS_j = ONE/double(ky_max-1);
                  ky = int(floor(S_j/dS_j));
                  S1_j = max(ZERO, double(ky)*dS_j);
                  S2_j = min(ONE, S1_j + dS_j);

		  if (S1_i - TOLER > ZERO) {
                     S1_NS.allocate(ky_max);
                     S1_NS.settype(SPLINE2D_LINEAR);

                     s_south = S1_i * smax_south + Grid.BndSouthSpline.sp[0];
                     S1_NS.Xp[0] = Spline(s_south, Grid.BndSouthSpline);
                     S1_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                     s_north = S1_i * smax_north + Grid.BndNorthSpline.sp[0],
                     S1_NS.Xp[ky_max-1] = Spline(s_north, Grid.BndNorthSpline);
                     S1_NS.tp[ky_max-1] = SPLINE2D_POINT_SHARP_CORNER;

                     for ( k = 1 ; k <= ky_max-2 ; ++k) {                     
                        s_east = double(k)*smax_east/double(ky_max-1) + Grid.BndEastSpline.sp[0];
                        x_east = Spline(s_east, Grid.BndEastSpline);

                        s_west = double(k)*smax_west/double(ky_max-1) + Grid.BndWestSpline.sp[0];
                        x_west = Spline(s_west, Grid.BndWestSpline);

                        w_west =  ONE-S1_i;
                        w_east =  S1_i;
                        w_total = w_east + w_west;
                        S1_NS.Xp[k] = (w_west*x_west + w_east*x_east)/w_total;
                        S1_NS.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
                     } /* endfor */

                     S1_NS.pathlength();
                     smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
                                  S1_NS.sp[0];
                  } else {
                     Copy_Spline(S1_NS, Grid.BndWestSpline);
                     smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
                                  S1_NS.sp[0];
                  } /* endif */

                  if (S2_i + TOLER < ONE) {
                     S2_NS.allocate(ky_max);
                     S2_NS.settype(SPLINE2D_LINEAR);

                     s_south = S2_i * smax_south + Grid.BndSouthSpline.sp[0];
                     S2_NS.Xp[0] = Spline(s_south, Grid.BndSouthSpline);
                     S2_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                     s_north = S2_i * smax_north + Grid.BndNorthSpline.sp[0],
                     S2_NS.Xp[ky_max-1] = Spline(s_north, Grid.BndNorthSpline);
                     S2_NS.tp[ky_max-1] = SPLINE2D_POINT_SHARP_CORNER;

                     for ( k = 1 ; k <= ky_max-2 ; ++k) {
                        s_east = double(k)*smax_east/double(ky_max-1) + Grid.BndEastSpline.sp[0];
                        x_east = Spline(s_east, Grid.BndEastSpline);

                        s_west = double(k)*smax_west/double(ky_max-1) + Grid.BndWestSpline.sp[0];
                        x_west = Spline(s_west, Grid.BndWestSpline);

                        w_west =  ONE-S2_i;
                        w_east =  S2_i;
                        w_total = w_east + w_west;
                        S2_NS.Xp[k] = (w_west*x_west + w_east*x_east)/w_total;
                        S2_NS.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
                     } /* endfor */

                     S2_NS.pathlength();
                     smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
                                  S2_NS.sp[0];
                  } else {
                     Copy_Spline(S2_NS, Grid.BndEastSpline);
                     smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
                                  S2_NS.sp[0];
                  } /* endif */

		  if (S1_j - TOLER > ZERO) {
                     S1_EW.allocate(2);
                     S1_EW.settype(SPLINE2D_LINEAR);

                     s_west = S1_j * smax_west + Grid.BndWestSpline.sp[0];
                     S1_EW.Xp[0] = Spline(s_west, Grid.BndWestSpline);
                     S1_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                     s_east = S1_j * smax_east + Grid.BndEastSpline.sp[0];
                     S1_EW.Xp[1] = Spline(s_east, Grid.BndEastSpline);
                     S1_EW.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

                     S1_EW.pathlength();
                     smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
                                  S1_EW.sp[0];
                  } else {
                     Copy_Spline(S1_EW, Grid.BndSouthSpline);
                     smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
                                  S1_EW.sp[0];
                  } /* endif */

                  if (S2_j + TOLER < ONE) {
                     S2_EW.allocate(2);
                     S2_EW.settype(SPLINE2D_LINEAR);

                     s_west = S2_j * smax_west + Grid.BndWestSpline.sp[0];
                     S2_EW.Xp[0] = Spline(s_west, Grid.BndWestSpline);
                     S2_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                     s_east = S2_j * smax_east + Grid.BndEastSpline.sp[0];
                     S2_EW.Xp[1] = Spline(s_east, Grid.BndEastSpline);
                     S2_EW.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

                     S2_EW.pathlength();
                     smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
                                  S2_EW.sp[0];
                  } else {
                     Copy_Spline(S2_EW, Grid.BndNorthSpline);
                     smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
                                  S2_EW.sp[0];
                  } /* endif */

                  s_north = S_i * smax_S2_EW + S2_EW.sp[0];
                  x_north = Spline(s_north, S2_EW);
                  s_south = S_i * smax_S1_EW + S1_EW.sp[0];
                  x_south = Spline(s_south, S1_EW);

                  s_east = S_j * smax_S2_NS + S2_NS.sp[0];
                  x_east = Spline(s_east, S2_NS);
                  s_west = S_j * smax_S1_NS + S1_NS.sp[0];
                  x_west = Spline(s_west, S1_NS);
      
                  if (kx == 0 || kx == kx_max - 2) {
                     w_west =  ONE-(S_i-S1_i)/(S2_i-S1_i);
                     w_east =  (S_i-S1_i)/(S2_i-S1_i);
                     w_total = w_east + w_west;
                     Grid.Node[i][j].X = (w_west*x_west + 
                                          w_east*x_east)/w_total;
                  } else {
                     w_south =  ONE-(S_j-S1_j)/(S2_j-S1_j);
                     w_north =  (S_j-S1_j)/(S2_j-S1_j);
                     w_total = w_south + w_north;
                     Grid.Node[i][j].X = (w_south*x_south + 
                                       w_north*x_north)/w_total;
                  } /* endif */

                  S1_NS.deallocate(); S2_NS.deallocate();
                  S2_EW.deallocate(); S2_EW.deallocate();
                //============================================================
                default:
                  w_south = (ONE-s_j/smax_j);
                  w_north = (s_j/smax_j);
                  w_total = w_south + w_north;
                  Grid.Node[i][j].X = (w_south*x_south + 
                                       w_north*x_north)/w_total;
                  break;
              } /* endswitch */

            } /* endif */

        } /* endfor */
    } /* endfor */

    /* Set the boundary condition types at the quadrilateral 
       grid block boundaries. */

    Set_BCs(Grid);

    /* Compute the exterior nodes for the quadrilateral mesh block. */

    Update_Exterior_Nodes(Grid);

    /* Compute the cells for the quadrilateral mesh block. */

    Update_Cells(Grid);

}

/********************************************************
 * Routine: Create_Quad_Block                           *
 *                                                      *
 * Create a 2D quadrilateral grid block with the four   *
 * boundaries of the block defined by blended splines   *
 * which are given as input arguments to the routine.   *
 *                                                      *
 ********************************************************/
void Create_Quad_Block(Grid2D_Quad_Block &Grid,
                       Spline2D &Bnd_Spline_North,
                       Spline2D &Bnd_Spline_South,
                       Spline2D &Bnd_Spline_East,
                       Spline2D &Bnd_Spline_West,
                       const int Number_of_Cells_Idir,
  	               const int Number_of_Cells_Jdir,
		       const int Number_of_Ghost_Cells,
                       const int Node_Init_Procedure,
                       const int Stretch_I,
                       const double &Beta_I, 
                       const double &Tau_I,
                       const int Stretch_J,
                       const double &Beta_J,
                       const double &Tau_J,
                       const int Orthogonal_North,
		       const int Orthogonal_South,
		       const int Orthogonal_East,
                       const int Orthogonal_West) {

    int i, j, k, kx, ky, kx_max, ky_max, npts, spline_type;
    double S_i, S_j, 
           s_north, s_south, s_east, s_west,
           smax_north, smax_south, smax_east, smax_west, 
           s_i, s_j, smax_i, smax_j,
           w_north, w_south, w_east, w_west, w_total;
    double smax_S1_NS, smax_S2_NS, smax_S1_EW, smax_S2_EW,
           S1_i, S2_i, S1_j, S2_j, dS_i, dS_j;
    Vector2D x_north, x_south, x_east, x_west;
    Spline2D S1_NS, S2_NS, S1_EW, S2_EW;

    /* Allocate (re-allocate) memory for the cells and nodes 
       of the quadrilateral mesh block. */

    if (Grid.Node != NULL && Grid.Cell != NULL) { 
       Grid.deallocate();
    } else if (Grid.Node != NULL) {
       Grid.deallocateNodes();
    } else if (Grid.Cell != NULL) {
       Grid.deallocateCells();
    } /* endif */
    Grid.allocate(Number_of_Cells_Idir, Number_of_Cells_Jdir, Number_of_Ghost_Cells);

    /* For each of the north, south, east, and west boundaries
       of this mesh block, assign the boundary splines specified
       by the subroutine input parameters. */

    Copy_Spline(Grid.BndNorthSpline, Bnd_Spline_North);
    Copy_Spline(Grid.BndSouthSpline, Bnd_Spline_South);
    Copy_Spline(Grid.BndEastSpline, Bnd_Spline_East);
    Copy_Spline(Grid.BndWestSpline, Bnd_Spline_West);

    Grid.SminN = Grid.BndNorthSpline.sp[0];
    Grid.SmaxN = Grid.BndNorthSpline.sp[Grid.BndNorthSpline.np-1];
    Grid.SminS = Grid.BndSouthSpline.sp[0];
    Grid.SmaxS = Grid.BndSouthSpline.sp[Grid.BndSouthSpline.np-1];
    Grid.SminE = Grid.BndEastSpline.sp[0];
    Grid.SmaxE = Grid.BndEastSpline.sp[Grid.BndEastSpline.np-1];
    Grid.SminW = Grid.BndWestSpline.sp[0];
    Grid.SmaxW = Grid.BndWestSpline.sp[Grid.BndWestSpline.np-1];

    /* Assign values for the type of node distribution 
       stretching functions to be used for each of the 
       coordinate directions. Also assign values to the
       related stretching parameters. */

    Grid.StretchI = Stretch_I;
    Grid.BetaI = Beta_I;
    Grid.TauI = Tau_I;
    Grid.StretchJ = Stretch_J;
    Grid.BetaJ = Beta_J;
    Grid.TauJ = Tau_J;

    /* Assign values to the grid orthogonality specifiers for
       the north, south, east, and west boundaries. */

    Grid.OrthogonalN = Orthogonal_North;
    Grid.OrthogonalS = Orthogonal_South;
    Grid.OrthogonalE = Orthogonal_East;
    Grid.OrthogonalW = Orthogonal_West;

    /* Compute the interior nodes for the quadrilateral mesh block. */

    smax_north = Grid.BndNorthSpline.sp[Grid.BndNorthSpline.np-1]-
                 Grid.BndNorthSpline.sp[0];
    smax_south = Grid.BndSouthSpline.sp[Grid.BndSouthSpline.np-1]-
                 Grid.BndSouthSpline.sp[0];

    smax_east = Grid.BndEastSpline.sp[Grid.BndEastSpline.np-1]-
                Grid.BndEastSpline.sp[0];
    smax_west = Grid.BndWestSpline.sp[Grid.BndWestSpline.np-1]-
                Grid.BndWestSpline.sp[0];

    for ( j = Grid.JNl ; j <= Grid.JNu ; ++j) {
        for ( i = Grid.INl ; i <= Grid.INu ; ++i) {
            S_j = double(j-Grid.JNl)/double(Grid.JNu-Grid.JNl);
            S_j = StretchingFcn(S_j, Grid.BetaJ, Grid.TauJ, Grid.StretchJ);
            s_east = S_j * smax_east + Grid.BndEastSpline.sp[0];
            x_east = Spline(s_east, Grid.BndEastSpline);
            s_west = S_j * smax_west + Grid.BndWestSpline.sp[0];
            x_west = Spline(s_west, Grid.BndWestSpline);

            S_i = double(i-Grid.INl)/double(Grid.INu-Grid.INl);
            S_i = StretchingFcn(S_i, Grid.BetaI, Grid.TauI, Grid.StretchI);
            s_north = S_i * smax_north + Grid.BndNorthSpline.sp[0];
            x_north = Spline(s_north, Grid.BndNorthSpline);
            s_south = S_i * smax_south + Grid.BndSouthSpline.sp[0];
            x_south = Spline(s_south, Grid.BndSouthSpline);

            if (j == Grid.JNl) {
	      Grid.Node[i][j].X = x_south;
            } else if (j == Grid.JNu) {
	      Grid.Node[i][j].X = x_north;
            } else if (i == Grid.INl) {
	      Grid.Node[i][j].X = x_west;
            } else if (i == Grid.INu) {
	      Grid.Node[i][j].X = x_east;
            } else {
              s_i = (ONE-S_j)*s_south + S_j*s_north;
              s_j = (ONE-S_i)*s_west + S_i*s_east;

              smax_i = (ONE-S_j)*smax_south + S_j*smax_north;
              smax_j = (ONE-S_i)*smax_west + S_i*smax_east;

              switch(Node_Init_Procedure) {
                //============================================================
                case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH :
                  w_south = (ONE-s_j/smax_j);
                  w_north = (s_j/smax_j);
                  w_total = w_south + w_north;
                  Grid.Node[i][j].X = (w_south*x_south + 
                                       w_north*x_north)/w_total;
                  break;
                //============================================================
                case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST :
                  w_west =  (ONE-s_i/smax_i);
                  w_east =  (s_i/smax_i);
                  w_total = w_east + w_west;
                  Grid.Node[i][j].X = (w_west*x_west + 
                                       w_east*x_east)/w_total;
                  break;
                //============================================================
                case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_XY :
                  kx_max = 5;
                  dS_i = ONE/double(kx_max-1);
                  kx = int(floor(S_i/dS_i));
                  S1_i = max(ZERO, double(kx)*dS_i);
                  S2_i = min(ONE, S1_i + dS_i);

                  ky_max = 7;
                  dS_j = ONE/double(ky_max-1);
                  ky = int(floor(S_j/dS_j));
                  S1_j = max(ZERO, double(ky)*dS_j);
                  S2_j = min(ONE, S1_j + dS_j);

		  if (S1_i - TOLER > ZERO) {
                     S1_NS.allocate(2);
                     S1_NS.settype(SPLINE2D_LINEAR);

                     s_south = S1_i * smax_south + Grid.BndSouthSpline.sp[0];
                     S1_NS.Xp[0] = Spline(s_south, Grid.BndSouthSpline);
                     S1_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                     s_north = S1_i * smax_north + Grid.BndNorthSpline.sp[0],
                     S1_NS.Xp[1] = Spline(s_north, Grid.BndNorthSpline);
                     S1_NS.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

                     S1_NS.pathlength();
                     smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
                                  S1_NS.sp[0];
                  } else {
                     Copy_Spline(S1_NS, Grid.BndWestSpline);
                     smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
                                  S1_NS.sp[0];
                  } /* endif */

                  if (S2_i + TOLER < ONE) {
                     S2_NS.allocate(2);
                     S2_NS.settype(SPLINE2D_LINEAR);

                     s_south = S2_i * smax_south + Grid.BndSouthSpline.sp[0];
                     S2_NS.Xp[0] = Spline(s_south, Grid.BndSouthSpline);
                     S2_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                     s_north = S2_i * smax_north + Grid.BndNorthSpline.sp[0],
                     S2_NS.Xp[1] = Spline(s_north, Grid.BndNorthSpline);
                     S2_NS.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

                     S2_NS.pathlength();
                     smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
                                  S2_NS.sp[0];
                  } else {
                     Copy_Spline(S2_NS, Grid.BndEastSpline);
                     smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
                                  S2_NS.sp[0];
                  } /* endif */

		  if (S1_j - TOLER > ZERO) {
                     S1_EW.allocate(kx_max);
                     S1_EW.settype(SPLINE2D_LINEAR);

                     s_west = S1_j * smax_west + Grid.BndWestSpline.sp[0];
                     S1_EW.Xp[0] = Spline(s_west, Grid.BndWestSpline);
                     S1_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                     s_east = S1_j * smax_east + Grid.BndEastSpline.sp[0];
                     S1_EW.Xp[kx_max-1] = Spline(s_east, Grid.BndEastSpline);
                     S1_EW.tp[kx_max-1] = SPLINE2D_POINT_SHARP_CORNER;

                     for ( k = 1 ; k <= kx_max-2 ; ++k) {                     
                        s_north = double(k)*smax_north/double(kx_max-1) + Grid.BndNorthSpline.sp[0];
                        x_north = Spline(s_north, Grid.BndNorthSpline);

                        s_south = double(k)*smax_south/double(kx_max-1) + Grid.BndSouthSpline.sp[0];
                        x_south = Spline(s_south, Grid.BndSouthSpline);

                        w_south =  ONE-S1_j;
                        w_north =  S1_j;
                        w_total = w_north + w_south;
                        S1_EW.Xp[k] = (w_south*x_south + w_north*x_north)/w_total;
                        S1_EW.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
                     } /* endfor */

                     S1_EW.pathlength();
                     smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
                                  S1_EW.sp[0];
                  } else {
                     Copy_Spline(S1_EW, Grid.BndSouthSpline);
                     smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
                                  S1_EW.sp[0];
                  } /* endif */

                  if (S2_j + TOLER < ONE) {
                     S2_EW.allocate(kx_max);
                     S2_EW.settype(SPLINE2D_LINEAR);

                     s_west = S2_j * smax_west + Grid.BndWestSpline.sp[0];
                     S2_EW.Xp[0] = Spline(s_west, Grid.BndWestSpline);
                     S2_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                     s_east = S2_j * smax_east + Grid.BndEastSpline.sp[0];
                     S2_EW.Xp[kx_max-1] = Spline(s_east, Grid.BndEastSpline);
                     S2_EW.tp[kx_max-1] = SPLINE2D_POINT_SHARP_CORNER;

                     for ( k = 1 ; k <= kx_max-2 ; ++k) {                     
                        s_north = double(k)*smax_north/double(kx_max-1) + Grid.BndNorthSpline.sp[0];
                        x_north = Spline(s_north, Grid.BndNorthSpline);

                        s_south = double(k)*smax_south/double(kx_max-1) + Grid.BndSouthSpline.sp[0];
                        x_south = Spline(s_south, Grid.BndSouthSpline);

                        w_south =  ONE-S2_j;
                        w_north =  S2_j;
                        w_total = w_north + w_south;
                        S2_EW.Xp[k] = (w_south*x_south + w_north*x_north)/w_total;
                        S2_EW.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
                     } /* endfor */

                     S2_EW.pathlength();
                     smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
                                  S2_EW.sp[0];
                  } else {
                     Copy_Spline(S2_EW, Grid.BndNorthSpline);
                     smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
                                  S2_EW.sp[0];
                  } /* endif */

                  s_north = S_i * smax_S2_EW + S2_EW.sp[0];
                  x_north = Spline(s_north, S2_EW);
                  s_south = S_i * smax_S1_EW + S1_EW.sp[0];
                  x_south = Spline(s_south, S1_EW);

                  s_east = S_j * smax_S2_NS + S2_NS.sp[0];
                  x_east = Spline(s_east, S2_NS);
                  s_west = S_j * smax_S1_NS + S1_NS.sp[0];
                  x_west = Spline(s_west, S1_NS);
      
                  if (ky == 0 || ky == ky_max - 2) {
                     w_south =  ONE-(S_j-S1_j)/(S2_j-S1_j);
                     w_north =  (S_j-S1_j)/(S2_j-S1_j);
                     w_total = w_south + w_north;
                     Grid.Node[i][j].X = (w_south*x_south + 
                                       w_north*x_north)/w_total;
                  } else {
                     w_west =  ONE-(S_i-S1_i)/(S2_i-S1_i);
                     w_east =  (S_i-S1_i)/(S2_i-S1_i);
                     w_total = w_east + w_west;
                     Grid.Node[i][j].X = (w_west*x_west + 
                                          w_east*x_east)/w_total;
                  } /* endif */

                  S1_NS.deallocate(); S2_NS.deallocate();
                  S2_EW.deallocate(); S2_EW.deallocate();
                  break;
                //============================================================
                case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_YX :
                  kx_max = 7;
                  dS_i = ONE/double(kx_max-1);
                  kx = int(floor(S_i/dS_i));
                  S1_i = max(ZERO, double(kx)*dS_i);
                  S2_i = min(ONE, S1_i + dS_i);

                  ky_max = 5;
                  dS_j = ONE/double(ky_max-1);
                  ky = int(floor(S_j/dS_j));
                  S1_j = max(ZERO, double(ky)*dS_j);
                  S2_j = min(ONE, S1_j + dS_j);

		  if (S1_i - TOLER > ZERO) {
                     S1_NS.allocate(ky_max);
                     S1_NS.settype(SPLINE2D_LINEAR);

                     s_south = S1_i * smax_south + Grid.BndSouthSpline.sp[0];
                     S1_NS.Xp[0] = Spline(s_south, Grid.BndSouthSpline);
                     S1_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                     s_north = S1_i * smax_north + Grid.BndNorthSpline.sp[0],
                     S1_NS.Xp[ky_max-1] = Spline(s_north, Grid.BndNorthSpline);
                     S1_NS.tp[ky_max-1] = SPLINE2D_POINT_SHARP_CORNER;

                     for ( k = 1 ; k <= ky_max-2 ; ++k) {                     
                        s_east = double(k)*smax_east/double(ky_max-1) + Grid.BndEastSpline.sp[0];
                        x_east = Spline(s_east, Grid.BndEastSpline);

                        s_west = double(k)*smax_west/double(ky_max-1) + Grid.BndWestSpline.sp[0];
                        x_west = Spline(s_west, Grid.BndWestSpline);

                        w_west =  ONE-S1_i;
                        w_east =  S1_i;
                        w_total = w_east + w_west;
                        S1_NS.Xp[k] = (w_west*x_west + w_east*x_east)/w_total;
                        S1_NS.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
                     } /* endfor */

                     S1_NS.pathlength();
                     smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
                                  S1_NS.sp[0];
                  } else {
                     Copy_Spline(S1_NS, Grid.BndWestSpline);
                     smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
                                  S1_NS.sp[0];
                  } /* endif */

                  if (S2_i + TOLER < ONE) {
                     S2_NS.allocate(ky_max);
                     S2_NS.settype(SPLINE2D_LINEAR);

                     s_south = S2_i * smax_south + Grid.BndSouthSpline.sp[0];
                     S2_NS.Xp[0] = Spline(s_south, Grid.BndSouthSpline);
                     S2_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                     s_north = S2_i * smax_north + Grid.BndNorthSpline.sp[0],
                     S2_NS.Xp[ky_max-1] = Spline(s_north, Grid.BndNorthSpline);
                     S2_NS.tp[ky_max-1] = SPLINE2D_POINT_SHARP_CORNER;

                     for ( k = 1 ; k <= ky_max-2 ; ++k) {
                        s_east = double(k)*smax_east/double(ky_max-1) + Grid.BndEastSpline.sp[0];
                        x_east = Spline(s_east, Grid.BndEastSpline);

                        s_west = double(k)*smax_west/double(ky_max-1) + Grid.BndWestSpline.sp[0];
                        x_west = Spline(s_west, Grid.BndWestSpline);

                        w_west =  ONE-S2_i;
                        w_east =  S2_i;
                        w_total = w_east + w_west;
                        S2_NS.Xp[k] = (w_west*x_west + w_east*x_east)/w_total;
                        S2_NS.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
                     } /* endfor */

                     S2_NS.pathlength();
                     smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
                                  S2_NS.sp[0];
                  } else {
                     Copy_Spline(S2_NS, Grid.BndEastSpline);
                     smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
                                  S2_NS.sp[0];
                  } /* endif */

		  if (S1_j - TOLER > ZERO) {
                     S1_EW.allocate(2);
                     S1_EW.settype(SPLINE2D_LINEAR);

                     s_west = S1_j * smax_west + Grid.BndWestSpline.sp[0];
                     S1_EW.Xp[0] = Spline(s_west, Grid.BndWestSpline);
                     S1_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                     s_east = S1_j * smax_east + Grid.BndEastSpline.sp[0];
                     S1_EW.Xp[1] = Spline(s_east, Grid.BndEastSpline);
                     S1_EW.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

                     S1_EW.pathlength();
                     smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
                                  S1_EW.sp[0];
                  } else {
                     Copy_Spline(S1_EW, Grid.BndSouthSpline);
                     smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
                                  S1_EW.sp[0];
                  } /* endif */

                  if (S2_j + TOLER < ONE) {
                     S2_EW.allocate(2);
                     S2_EW.settype(SPLINE2D_LINEAR);

                     s_west = S2_j * smax_west + Grid.BndWestSpline.sp[0];
                     S2_EW.Xp[0] = Spline(s_west, Grid.BndWestSpline);
                     S2_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                     s_east = S2_j * smax_east + Grid.BndEastSpline.sp[0];
                     S2_EW.Xp[1] = Spline(s_east, Grid.BndEastSpline);
                     S2_EW.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

                     S2_EW.pathlength();
                     smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
                                  S2_EW.sp[0];
                  } else {
                     Copy_Spline(S2_EW, Grid.BndNorthSpline);
                     smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
                                  S2_EW.sp[0];
                  } /* endif */

                  s_north = S_i * smax_S2_EW + S2_EW.sp[0];
                  x_north = Spline(s_north, S2_EW);
                  s_south = S_i * smax_S1_EW + S1_EW.sp[0];
                  x_south = Spline(s_south, S1_EW);

                  s_east = S_j * smax_S2_NS + S2_NS.sp[0];
                  x_east = Spline(s_east, S2_NS);
                  s_west = S_j * smax_S1_NS + S1_NS.sp[0];
                  x_west = Spline(s_west, S1_NS);
      
                  if (kx == 0 || kx == kx_max - 2) {
                     w_west =  ONE-(S_i-S1_i)/(S2_i-S1_i);
                     w_east =  (S_i-S1_i)/(S2_i-S1_i);
                     w_total = w_east + w_west;
                     Grid.Node[i][j].X = (w_west*x_west + 
                                          w_east*x_east)/w_total;
                  } else {
                     w_south =  ONE-(S_j-S1_j)/(S2_j-S1_j);
                     w_north =  (S_j-S1_j)/(S2_j-S1_j);
                     w_total = w_south + w_north;
                     Grid.Node[i][j].X = (w_south*x_south + 
                                       w_north*x_north)/w_total;
                  } /* endif */

                  S1_NS.deallocate(); S2_NS.deallocate();
                  S2_EW.deallocate(); S2_EW.deallocate();
                  break;
                //============================================================
                default:
                  w_south = (ONE-s_j/smax_j);
                  w_north = (s_j/smax_j);
                  w_total = w_south + w_north;
                  Grid.Node[i][j].X = (w_south*x_south + 
                                       w_north*x_north)/w_total;
                  break;
              } /* endswitch */

            } /* endif */

        } /* endfor */
    } /* endfor */

    /* Set the boundary condition types at the quadrilateral 
       grid block boundaries. */

    Set_BCs(Grid);

    /* Compute the exterior nodes for the quadrilateral mesh block. */

    Update_Exterior_Nodes(Grid);

    /* Compute the cells for the quadrilateral mesh block. */

    Update_Cells(Grid);

}

/********************************************************
 * Routine: Broadcast_Quad_Block                        *
 *                                                      *
 * Broadcast quadrilateral grid block to all            *
 * processors involved in the calculation from the      *
 * primary processor using the MPI broadcast routine.   *
 *                                                      *
 ********************************************************/
void Broadcast_Quad_Block(Grid2D_Quad_Block &Grid) {

#ifdef _MPI_VERSION

    int i, j, ni, nj, ng, mesh_allocated, buffer_size;
    double *buffer;
 
    /* Broadcast the number of cells in each direction. */

    if (CFFC_Primary_MPI_Processor()) {
      ni = Grid.NCi;
      nj = Grid.NCj;
      ng = Grid.Nghost;
      if (Grid.Node != NULL) {
         mesh_allocated = 1;
      } else {
	 mesh_allocated = 0;
      } /* endif */ 
    } /* endif */

    MPI::COMM_WORLD.Bcast(&ni, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&nj, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&ng, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&mesh_allocated, 1, MPI::INT, 0);

    /* On non-primary MPI processors, allocate (re-allocate) 
       memory for the cells and nodes of the quadrilateral 
       mesh block as necessary. */

    if (!CFFC_Primary_MPI_Processor()) {
       if (Grid.NCi != ni || Grid.NCj != nj || Grid.Nghost != ng) { 
          if (Grid.Node != NULL && Grid.Cell != NULL) { 
             Grid.deallocate();
          } else if (Grid.Node != NULL) {
             Grid.deallocateNodes();
          } else if (Grid.Cell != NULL) {
             Grid.deallocateCells();
          } /* endif */
          if (mesh_allocated) Grid.allocate(ni-2*ng, nj-2*ng, ng); 
       } /* endif */
    } /* endif */

    /* Broadcast the north, south, east, and west 
       boundary splines. */

    Broadcast_Spline(Grid.BndNorthSpline);
    Broadcast_Spline(Grid.BndSouthSpline);
    Broadcast_Spline(Grid.BndEastSpline);
    Broadcast_Spline(Grid.BndWestSpline);

   /* Broadcast min/max pathlengths for boundary splines. */

    MPI::COMM_WORLD.Bcast(&(Grid.SminN), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Grid.SmaxN), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Grid.SminS), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Grid.SmaxS), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Grid.SminE), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Grid.SmaxE), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Grid.SminW), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Grid.SmaxW), 1, MPI::DOUBLE, 0);

    /* Broadcast node control parameters. */

    MPI::COMM_WORLD.Bcast(&(Grid.StretchI), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Grid.BetaI), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Grid.TauI), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Grid.StretchJ), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Grid.BetaJ), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(Grid.TauJ), 1, MPI::DOUBLE, 0);

    MPI::COMM_WORLD.Bcast(&(Grid.OrthogonalN), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Grid.OrthogonalS), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Grid.OrthogonalE), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(&(Grid.OrthogonalW), 1, MPI::INT, 0);

    /* Broadcast the node locations for grid block. */

    if (mesh_allocated) {
       ni = (Grid.INu+Grid.Nghost) - (Grid.INl-Grid.Nghost) + 1;
       nj = (Grid.JNu+Grid.Nghost) - (Grid.JNl-Grid.Nghost) + 1;
       buffer = new double[2*ni*nj];

       if (CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (j  = Grid.JNl-Grid.Nghost; j <= Grid.JNu+Grid.Nghost; ++j ) {
              for ( i = Grid.INl-Grid.Nghost; i <= Grid.INu+Grid.Nghost; ++i ) {
 	          buffer[buffer_size] = Grid.Node[i][j].X.x;
 	          buffer[buffer_size+1] = Grid.Node[i][j].X.y;
                  buffer_size = buffer_size + 2;
              } /* endfor */
          } /* endfor */
       } /* endif */

       buffer_size = 2*ni*nj;
       MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

       if (!CFFC_Primary_MPI_Processor()) {
          buffer_size = 0;
          for (j  = Grid.JNl-Grid.Nghost; j <= Grid.JNu+Grid.Nghost; ++j ) {
              for ( i = Grid.INl-Grid.Nghost; i <= Grid.INu+Grid.Nghost; ++i ) {
 	          Grid.Node[i][j].X.x = buffer[buffer_size];
 	          Grid.Node[i][j].X.y = buffer[buffer_size+1];
                  buffer_size = buffer_size + 2;
              } /* endfor */
          } /* endfor */
       } /* endif */

       delete []buffer;
       buffer = NULL;

    } /* endif */

    /*  On non-primary MPI processors, set the boundary condition types
        compute the exterior nodes, and compute the cells for the 
        quadrilateral mesh block. */

    if (mesh_allocated && !CFFC_Primary_MPI_Processor()) {
       Set_BCs(Grid);
       Update_Exterior_Nodes(Grid);
       Update_Cells(Grid);
    } /* endif */
#endif

}

#ifdef _MPI_VERSION
/********************************************************
 * Routine: Broadcast_Quad_Block                        *
 *                                                      *
 * Broadcast quadrilateral grid block to all processors *
 * associated with the specified communicator from the  *
 * specified processor using the MPI broadcast routine. *
 *                                                      *
 ********************************************************/
void Broadcast_Quad_Block(Grid2D_Quad_Block &Grid,
                          MPI::Intracomm &Communicator, 
                          const int Source_CPU) {

    int Source_Rank = 0;
    int i, j, ni, nj, ng, mesh_allocated, buffer_size;
    double *buffer;
 
    /* Broadcast the number of cells in each direction. */

    if (CFFC_MPI::This_Processor_Number == Source_CPU) {
      ni = Grid.NCi;
      nj = Grid.NCj;
      ng = Grid.Nghost;
      if (Grid.Node != NULL) {
         mesh_allocated = 1;
      } else {
	 mesh_allocated = 0;
      } /* endif */ 
    } /* endif */

    Communicator.Bcast(&ni, 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&nj, 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&ng, 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&mesh_allocated, 1, MPI::INT, Source_Rank);

    /* On non-source MPI processors, allocate (re-allocate) 
       memory for the cells and nodes of the quadrilateral 
       mesh block as necessary. */

    if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
       if (Grid.NCi != ni || Grid.NCj != nj || Grid.Nghost != ng) { 
          if (Grid.Node != NULL && Grid.Cell != NULL) { 
             Grid.deallocate();
          } else if (Grid.Node != NULL) {
             Grid.deallocateNodes();
          } else if (Grid.Cell != NULL) {
             Grid.deallocateCells();
          } /* endif */
          if (mesh_allocated) Grid.allocate(ni-2*ng, nj-2*ng, ng); 
       } /* endif */
    } /* endif */

    /* Broadcast the north, south, east, and west 
       boundary splines. */

    Broadcast_Spline(Grid.BndNorthSpline, Communicator, Source_CPU);
    Broadcast_Spline(Grid.BndSouthSpline, Communicator, Source_CPU);
    Broadcast_Spline(Grid.BndEastSpline, Communicator, Source_CPU);
    Broadcast_Spline(Grid.BndWestSpline, Communicator, Source_CPU);

   /* Broadcast min/max pathlengths for boundary splines. */

    Communicator.Bcast(&(Grid.SminN), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(Grid.SmaxN), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(Grid.SminS), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(Grid.SmaxS), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(Grid.SminE), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(Grid.SmaxE), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(Grid.SminW), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(Grid.SmaxW), 1, MPI::DOUBLE, Source_Rank);

    /* Broadcast node control parameters. */

    Communicator.Bcast(&(Grid.StretchI), 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&(Grid.BetaI), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(Grid.TauI), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(Grid.StretchJ), 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&(Grid.BetaJ), 1, MPI::DOUBLE, Source_Rank);
    Communicator.Bcast(&(Grid.TauJ), 1, MPI::DOUBLE, Source_Rank);

    Communicator.Bcast(&(Grid.OrthogonalN), 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&(Grid.OrthogonalS), 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&(Grid.OrthogonalE), 1, MPI::INT, Source_Rank);
    Communicator.Bcast(&(Grid.OrthogonalW), 1, MPI::INT, Source_Rank);

    /* Broadcast the node locations for grid block. */

    if (mesh_allocated) {
       ni = (Grid.INu+Grid.Nghost) - (Grid.INl-Grid.Nghost) + 1;
       nj = (Grid.JNu+Grid.Nghost) - (Grid.JNl-Grid.Nghost) + 1;
       buffer = new double[2*ni*nj];

       if (CFFC_MPI::This_Processor_Number == Source_CPU) {
          buffer_size = 0;
          for (j  = Grid.JNl-Grid.Nghost ; j <= Grid.JNu+Grid.Nghost ; ++j ) {
              for ( i = Grid.INl-Grid.Nghost ; i <= Grid.INu+Grid.Nghost ; ++i ) {
 	          buffer[buffer_size] = Grid.Node[i][j].X.x;
 	          buffer[buffer_size+1] = Grid.Node[i][j].X.y;
                  buffer_size = buffer_size + 2;
              } /* endfor */
          } /* endfor */
       } /* endif */

       buffer_size = 2*ni*nj;
       Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

       if (!(CFFC_MPI::This_Processor_Number == Source_CPU)) {
          buffer_size = 0;
          for (j  = Grid.JNl-Grid.Nghost; j <= Grid.JNu+Grid.Nghost; ++j ) {
              for ( i = Grid.INl-Grid.Nghost; i <= Grid.INu+Grid.Nghost; ++i ) {
 	          Grid.Node[i][j].X.x = buffer[buffer_size];
 	          Grid.Node[i][j].X.y = buffer[buffer_size+1];
                  buffer_size = buffer_size + 2;
              } /* endfor */
          } /* endfor */
       } /* endif */

       delete []buffer; 
       buffer = NULL;

    } /* endif */

    /*  On non-source MPI processors, set the boundary condition types
        compute the exterior nodes, and compute the cells for the 
        quadrilateral mesh block. */

    if (mesh_allocated && 
        !(CFFC_MPI::This_Processor_Number == Source_CPU)) {
       Set_BCs(Grid);
       Update_Exterior_Nodes(Grid);
       Update_Cells(Grid);
    } /* endif */

}
#endif

/********************************************************
 * Routine: Smooth_Quad_Block                           *
 *                                                      *
 * Smooths quadrilateral grid block using an elliptic   *
 * grid generation method (transformed Poisson's        *
 * equation method) with a successive over-relaxation   *
 * (SOR) Gauss-Seidel solution procedure.  The routine  *
 * follows the technique outlined by Sorenson (1980)    *
 * and enforces fixed grid spacing and orthogonality at * 
 * the bottom (south) and top (north) boundaries of the *
 * quadrilateral block as required.  Modifications have *
 * been introduced to permit the enforcement of grid    *
 * orthogonality at the left (west) and right (east)    *
 * boundaries as well.                                  *
 *                                                      *
 ********************************************************/
void Smooth_Quad_Block(Grid2D_Quad_Block &Grid,
  	               const int Number_of_Iterations) {

 /* Local Variable description:

    xij, yij                  Two-dimensional arrays containing
                              the coordinates of the nodes for
                              the grid block in physical (x,y) space.

    deltas_b, deltas_t,       One-dimensional arrays containing the
    deltas_l, deltas_r        node spacing in the normal direction
                              between the first two nodes at the upper
                              and lower boundaries of the grid block.

    xz_b, yz_b, xe_b, ye_b,   One-dimensional arrays containing the
    xzz_b, yzz_b, xee_b,      partial derivatives of x and y
    yee_b, xze_b, yze_b       coordinates of the physical space
                              with respect to the zeta and eta
    xz_t, yz_t, xe_t, ye_t,   coordinates of the computational space,
    xzz_t, yzz_t, xee_t,      all evaluated at the bottom (south),
    yee_t, xze_t, yze_t       top (north), left (west), and right 
                              (east) boundaries of the grid block.

    xz_l, yz_l, xe_l, ye_l,
    xzz_l, yzz_l, xee_l,   
    yee_l, xze_l, yze_l

    xz_r, yz_r, xe_r, ye_r,
    xzz_r, yzz_r, xee_r,   
    yee_r, xze_r, yze_r

    jacob_b, jacob_t,         One-dimensional arrays containing the
    jacob_l, jacob_r          Jacobian of the coordinate transformation
                              evaluated at the bottom (south), top 
                              (north), left (west), and right (east)
                              boundaries of the grid block.

    pb, pt, pl, pr,           One-dimensional arrays containing the
    qb, qt, ql, qr            values of the source terms of the
                              transformed Poisson's equations at the
                              bottom (south), top (north), left (west), 
                              and right (east) boundaries of the grid 
                              block.

    norm_b, norm_t,           Integer parameters indicating whether 
    norm_l, norm_r            orthogonality at the bottom (south), top
                              (north), left (west), and right (east) 
                              boundaries is to be enforced.  If value 
                              is zero then orthogonality is not enforced 
                              and any other value causes the orthogonality 
                              condition to be applied.

    wsor, wsur                SOR and successive under-relaxation (SUR)
                              iteration parameters.  wsor is used to
                              over-relax the solution for x and y and
                              wsur is used to under-relax the solution
                              for p, q, r, and s.  Typically, 
                              wsor=1.5-1.8 and wsur=0.30-0.70.

    a, b,                     Exponential decay parameters defining the
    c, d,                     effects of the source terms away from the
    e                         bottom (south), top (north), left (west), 
                              and right (east) boundaries, as well as the 
                              corners.  If a (or b, c, d, e) is small 
                              (i.e., a = 0.2) source terms have effects far 
                              from boundary.  If a (or b, c, d, e) is large
                              (i.e., a = 0.7) source terms have less 
                              effects far from boundary.

    n_zeta, n_eta             Number of nodes in the zeta and eta
                              directions of the computational domain for 
                              grid block.

    alpha, beta, gamma        Parameters used in the transformed Poisson's
                              equations.

    r1, r2, r3, r4            Parameters used to evaluate the source terms
                              of the transformed Poisson's equations at
                              the bottom (south), top (north), left (west), 
                              and right (east) boundaries of the grid block.

    dpb, dpt, dpl, dpr,       Changes in pb, pt, pl, pr, qb, qt, ql, and qr
    dqb, dqt, dql, dqr        before the under-relaxation and limiting
                              algorithm is applied.

    dpmax, dqmax              Maximum permitted changes in p and q at the
                              bottom, top, left, and right boundaries.
 
    xzij, yzij, xeij, yeij,   Local values of the transformation metrics
    jacobij                   and Jacobian.

    pij, qij                  Local values of the source terms.

    dxij, dyij                Changes in xij and yij for on step in the
                              iteration technique before the over-
                              relaxation algorithm is applied.

    fb, ft, fl, fr,           Exponentials use to evaluate the local
    fbl, fbr, ftl, ftr        values of the source terms. */

    double **xij, **yij,
           *deltas_b, *deltas_t, *deltas_l, *deltas_r,
           *xz_b, *yz_b, *xe_b, *ye_b, 
           *xzz_b, *yzz_b, *xee_b, *yee_b, *xze_b, *yze_b,
           *xz_t, *yz_t, *xe_t, *ye_t,
           *xzz_t, *yzz_t, *xee_t, *yee_t, *xze_t, *yze_t,
           *xz_l, *yz_l, *xe_l, *ye_l,
           *xzz_l, *yzz_l, *xee_l, *yee_l, *xze_l, *yze_l,
           *xz_r, *yz_r, *xe_r, *ye_r,
           *xzz_r, *yzz_r, *xee_r, *yee_r, *xze_r, *yze_r,
           *pb, *pt, *pl, *pr,
           *qb, *qt, *ql, *qr,
           *jacob_b, *jacob_t, *jacob_l, *jacob_r;

    int num_iter, i_poisson, 
        norm_b, norm_t, norm_l, norm_r,
        n_zeta, n_eta;
    double wsor, wsur, a, b, c, d, e, l_damping;

    int i, j, i_step, ii, jj;
    double alpha, beta, gamma;
    double r1, r2, r3, r4;
    double dpb, dpt, dpl, dpr, dqb, dqt, dql, dqr;
    double dpmax, dqmax;
    double xzij, yzij, xeij, yeij, jacobij, pij, qij, dxij, dyij;
    double fb, ft, fl, fr, fbl, fbr, ftl, ftr;

    /* Determine the size of and create the Poisson equation solution 
       variables for the elliptic smoothing algorithm. */

    n_zeta = Grid.INu - Grid.INl + 1;
    n_eta = Grid.JNu - Grid.JNl + 1;

    xij = new double*[n_zeta];
    yij = new double*[n_zeta];
    for ( i = 0; i <= n_zeta-1 ; ++i ) {
        xij[i] = new double[n_eta];
        yij[i] = new double[n_eta];
    } /* endfor */

    deltas_b = new double[n_zeta];
    deltas_t = new double[n_zeta];
    deltas_l = new double[n_eta];
    deltas_r = new double[n_eta];

    xz_b = new double[n_zeta]; 
    yz_b = new double[n_zeta];
    xe_b = new double[n_zeta];
    ye_b = new double[n_zeta];
    xzz_b = new double[n_zeta];
    yzz_b = new double[n_zeta];
    xee_b = new double[n_zeta];
    yee_b = new double[n_zeta];
    xze_b = new double[n_zeta];
    yze_b = new double[n_zeta];

    xz_t = new double[n_zeta];
    yz_t = new double[n_zeta];
    xe_t = new double[n_zeta]; 
    ye_t = new double[n_zeta];
    xzz_t = new double[n_zeta];
    yzz_t = new double[n_zeta];
    xee_t = new double[n_zeta];
    yee_t = new double[n_zeta];
    xze_t = new double[n_zeta];
    yze_t = new double[n_zeta];

    xz_l = new double[n_eta];
    yz_l = new double[n_eta];
    xe_l = new double[n_eta];
    ye_l = new double[n_eta];
    xzz_l = new double[n_eta];
    yzz_l = new double[n_eta];
    xee_l = new double[n_eta];
    yee_l = new double[n_eta];
    xze_l = new double[n_eta];
    yze_l = new double[n_eta];

    xz_r = new double[n_eta]; 
    yz_r = new double[n_eta];
    xe_r = new double[n_eta];
    ye_r = new double[n_eta];
    xzz_r = new double[n_eta];
    yzz_r = new double[n_eta];
    xee_r = new double[n_eta];
    yee_r = new double[n_eta];
    xze_r = new double[n_eta];
    yze_r = new double[n_eta];

    pb = new double[n_zeta];
    pt = new double[n_zeta];
    pl = new double[n_eta];
    pr = new double[n_eta];

    qb = new double[n_zeta]; 
    qt = new double[n_zeta]; 
    ql = new double[n_eta];
    qr = new double[n_eta];

    jacob_b = new double[n_zeta];
    jacob_t = new double[n_zeta];
    jacob_l = new double[n_eta];
    jacob_r = new double[n_eta];

    /* Assign initial values to the Poisson equation solution
       variables using the current node locations for the nodes
       in the quadrilateral grid block. */

    for (j  = Grid.JNl ; j <= Grid.JNu ; ++j ) {
       for ( i = Grid.INl ; i <= Grid.INu ; ++i ) {
	   ii = i - Grid.INl;
	   jj = j - Grid.JNl;
           xij[ii][jj] = Grid.Node[i][j].X.x;
           yij[ii][jj] = Grid.Node[i][j].X.y;
       } /* endfor */
    } /* endfor */

    /* Set the Poisson equation solution parameters. */

    norm_b = Grid.OrthogonalS;
    norm_t = Grid.OrthogonalN;
    norm_l = Grid.OrthogonalW;
    norm_r = Grid.OrthogonalE;

    wsor = 1.25;
    wsur = 0.10;
    l_damping = 1.00;

    if (norm_b == 0 && norm_t == 0 &&
        norm_l == 0 && norm_r == 0) {
      i_poisson = 0;
      a = ONE;
      b = ONE;
      c = ONE;
      d = ONE;
      e = HALF;
    } else if (norm_b != 0 && norm_t == 0 &&
               norm_l == 0 && norm_r == 0) {
      i_poisson = 1;
      a = l_damping;
      b = ONE;
      c = ONE;
      d = ONE;
      e = HALF;
    } else if (norm_b == 0 && norm_t != 0 &&
               norm_l == 0 && norm_r == 0) {
      i_poisson = 1;
      a = ONE;
      b = l_damping;
      c = ONE;
      d = ONE;
      e = HALF;
    } else if (norm_b == 0 && norm_t == 0 &&
               norm_l != 0 && norm_r == 0) {
      i_poisson = 1;
      a = ONE;
      b = ONE;
      c = l_damping;
      d = ONE;
      e = HALF;
    } else if (norm_b == 0 && norm_t == 0 &&
               norm_l == 0 && norm_r != 0) {
      i_poisson = 1;
      a = ONE;
      b = ONE;
      c = ONE;
      d = l_damping;
      e = HALF;
    } else if (norm_b != 0 && norm_t != 0 &&
               norm_l == 0 && norm_r == 0) {
      i_poisson = 1;
      a = l_damping;
      b = l_damping;
      c = ONE;
      d = ONE;
      e = HALF;
    } else if (norm_b == 0 && norm_t == 0 &&
               norm_l != 0 && norm_r != 0) {
      i_poisson = 1;
      a = ONE;
      b = ONE;
      c = l_damping;
      d = l_damping;
      e = HALF;
    } else if (norm_b != 0 && norm_t == 0 &&
               norm_l != 0 && norm_r == 0) {
      i_poisson = 1;
      a = l_damping;
      b = ONE;
      c = l_damping;
      d = ONE;
      e = HALF;
    } else if (norm_b == 0 && norm_t != 0 &&
               norm_l == 0 && norm_r != 0) {
      i_poisson = 1;
      a = ONE;
      b = l_damping;
      c = ONE;
      d = l_damping;
      e = HALF;
    } else if (norm_b == 0 && norm_t != 0 &&
               norm_l != 0 && norm_r != 0) {
      i_poisson = 1;
      a = ONE;
      b = l_damping;
      c = l_damping;
      d = l_damping;
      e = HALF;
    } else if (norm_b != 0 && norm_t == 0 &&
               norm_l != 0 && norm_r != 0) {
      i_poisson = 1;
      a = l_damping;
      b = ONE;
      c = l_damping;
      d = l_damping;
      e = HALF;
    } else if (norm_b != 0 && norm_t != 0 &&
               norm_l == 0 && norm_r != 0) {
      i_poisson = 1;
      a = l_damping;
      b = l_damping;
      c = ONE;
      d = l_damping;
      e = HALF;
    } else if (norm_b != 0 && norm_t != 0 &&
               norm_l != 0 && norm_r == 0) {
      i_poisson = 1;
      a = l_damping;
      b = l_damping;
      c = l_damping;
      d = ONE;
      e = HALF;
    } else {
      i_poisson = 1;
      a = l_damping;
      b = l_damping;
      c = l_damping;
      d = l_damping;
      e = HALF;
    } /* endif */

    /* Initialize the partial derivatives and source terms at the 
       bottom (south), top (north), left (west), and right (east) 
       boundaries.  The procedure for doing this depends on whether 
       or not orthogonality is enforced at the bottom (south), 
       top (north), left (west), and right (east) boundaries. */

    /* BOTTOM AND TOP */

    xz_b[0] = HALF * (xij[1][0] - xij[0][0]);
    xz_t[0] = HALF * (xij[1][n_eta-1] - xij[0][n_eta-1]);

    yz_b[0] = HALF * (yij[1][0] - yij[0][0]);
    yz_t[0] = HALF * (yij[1][n_eta-1] - yij[0][n_eta-1]);

    xe_b[0] = ZERO;
    xe_t[0] = ZERO;

    ye_b[0] = ZERO;
    ye_t[0] = ZERO;

    for ( i = 1; i < n_zeta - 1; ++i ) {
        xz_b[i] = HALF * (xij[i+1][0] - xij[i-1][0]);
        xz_t[i] = HALF * (xij[i+1][n_eta-1] - xij[i-1][n_eta-1]);
 
        yz_b[i] = HALF * (yij[i+1][0] - yij[i-1][0]);
        yz_t[i] = HALF * (yij[i+1][n_eta-1] - yij[i-1][n_eta-1]);

        xzz_b[i] = xij[i+1][0] - TWO * xij[i][0] + xij[i-1][0];
        xzz_t[i] = xij[i+1][n_eta-1] - TWO * xij[i][n_eta-1] +
                   xij[i-1][n_eta-1];

        yzz_b[i] = yij[i+1][0] - TWO * yij[i][0] + yij[i-1][0];
        yzz_t[i] = yij[i+1][n_eta-1] - TWO * yij[i][n_eta-1] +
                   yij[i-1][n_eta-1];
 
        if (norm_b == 0) {
           xe_b[i] = xij[i][1] - xij[i][0];
           ye_b[i] = yij[i][1] - yij[i][0];
        } else {
           deltas_b[i] = hypot( (xij[i][1]-xij[i][0]),
                                (yij[i][1]-yij[i][0]) );
           xe_b[i] = - deltas_b[i] * yz_b[i] /
                       hypot( xz_b[i], yz_b[i] );
           ye_b[i] = deltas_b[i] * xz_b[i] /
                     hypot( xz_b[i], yz_b[i] );
        } /* endif */

        if (norm_t == 0) {
          xe_t[i] = xij[i][n_eta-1] - xij[i][n_eta-2];
          ye_t[i] = yij[i][n_eta-1] - yij[i][n_eta-2];
        } else {
          deltas_t[i] = hypot( (xij[i][n_eta-1]-xij[i][n_eta-2]),
                               (yij[i][n_eta-1]-yij[i][n_eta-2]) );
          xe_t[i] = - deltas_t[i] * yz_t[i] /
                      hypot( xz_t[i], yz_t[i] );
          ye_t[i] = deltas_t[i] * xz_t[i] /
                    hypot( xz_t[i], yz_t[i] );
        } /* endif */      

        jacob_b[i] = xz_b[i] * ye_b[i] - xe_b[i] * yz_b[i];
        jacob_t[i] = xz_t[i] * ye_t[i] - xe_t[i] * yz_t[i];

        pb[i] = ZERO;
        pt[i] = ZERO;

        qb[i] = ZERO;
        qt[i] = ZERO;
    } /* endfor */

    xz_b[n_zeta-1] = HALF * (xij[n_zeta-1][0] - xij[n_zeta-2][0]);
    xz_t[n_zeta-1] = HALF * (xij[n_zeta-1][n_eta-1] -
                             xij[n_zeta-2][n_eta-1]);

    yz_b[n_zeta-1] = HALF * (yij[n_zeta-1][0] - yij[n_zeta-2][0]);
    yz_t[n_zeta-1] = HALF * (yij[n_zeta-1][n_eta-1] -
                             yij[n_zeta-2][n_eta-1]);

    xe_b[n_zeta-1] = ZERO;
    xe_t[n_zeta-1] = ZERO;

    ye_b[n_zeta-1] = ZERO;
    ye_t[n_zeta-1] = ZERO;

    for ( i = 1; i < n_zeta - 1; ++i ) {
        xze_b[i] = HALF * (xe_b[i+1] - xe_b[i-1]);
        xze_t[i] = HALF * (xe_t[i+1] - xe_t[i-1]);
 
        yze_b[i] = HALF * (ye_b[i+1] - ye_b[i-1]);
        yze_t[i] = HALF * (ye_t[i+1] - ye_t[i-1]);
    } /* endfor */

    /* LEFT AND RIGHT */

    xe_l[0] = HALF * (xij[0][1] - xij[0][0]);
    xe_r[0] = HALF * (xij[n_zeta-1][1] - xij[n_zeta-1][0]);

    ye_l[0] = HALF * (yij[0][1] - yij[0][0]);
    ye_r[0] = HALF * (yij[n_zeta-1][1] - yij[n_zeta-1][0]);

    xz_l[0] = ZERO;
    xz_r[0] = ZERO;

    yz_l[0] = ZERO;
    yz_r[0] = ZERO;

    for ( j = 1; j < n_eta - 1; ++j ) {
        xe_l[j] = HALF * (xij[0][j+1] - xij[0][j-1]);
        xe_r[j] = HALF * (xij[n_zeta-1][j+1] - xij[n_zeta-1][j-1]);

        ye_l[j] = HALF * (yij[0][j+1] - yij[0][j-1]);
        ye_r[j] = HALF * (yij[n_zeta-1][j+1] - yij[n_zeta-1][j-1]);

        xee_l[j] = xij[0][j+1] - TWO * xij[0][j] + xij[0][j-1];
        xee_r[j] = xij[n_zeta-1][j+1] - TWO * xij[n_zeta-1][j] +
                   xij[n_zeta-1][j-1];

        yee_l[j] = yij[0][j+1] - TWO * yij[0][j] + yij[0][j-1];
        yee_r[j] = yij[n_zeta-1][j+1] - TWO * yij[n_zeta-1][j] +
                   yij[n_zeta-1][j-1];

        if (norm_l == 0) {
           xz_l[j] = xij[1][j] - xij[0][j];
           yz_l[j] = yij[1][j] - yij[0][j];
        } else {
           deltas_l[j] = hypot( (xij[1][j]-xij[0][j]),
                                (yij[1][j]-yij[0][j]) );
           xz_l[j] = deltas_l[j] * ye_l[j] /
                     hypot( xe_l[j], ye_l[j] );
           yz_l[j] = - deltas_l[j] * xe_l[j] /
                       hypot( xe_l[j], ye_l[j] );
        } /* endif */

        if (norm_r == 0) {
           xz_r[j] = xij[n_zeta-1][j] - xij[n_zeta-2][j];
           yz_r[j] = yij[n_zeta-1][j] - yij[n_zeta-2][j];
        } else {
           deltas_r[j] = hypot( (xij[n_zeta-1][j]-xij[n_zeta-2][j]),
                                (yij[n_zeta-1][j]-yij[n_zeta-2][j]) );
           xz_r[j] = deltas_r[j] * ye_r[j] /
                     hypot( xe_r[j], ye_r[j] );
           yz_r[j] = -deltas_r[j] * xe_r[j] /
                      hypot( xe_r[j], ye_r[j] );
        } /* endif */

        jacob_l[j] = xz_l[j] * ye_l[j] - xe_l[j] * yz_l[j];
        jacob_r[j] = xz_r[j] * ye_r[j] - xe_r[j] * yz_r[j];

        pl[j] = ZERO;
        pr[j] = ZERO;

        ql[j] = ZERO;
        qr[j] = ZERO;
    } /* endfor */

    xe_l[n_eta-1] = HALF * (xij[0][n_eta-1] - xij[0][n_eta-2]);
    xe_r[n_eta-1] = HALF * (xij[n_zeta-1][n_eta-1] -
                            xij[n_zeta-1][n_eta-2]);

    ye_l[n_eta-1] = HALF * (yij[0][n_eta-1] - yij[0][n_eta-2]);
    ye_r[n_eta-1] = HALF * (yij[n_zeta-1][n_eta-1] -
                            yij[n_zeta-1][n_eta-2]);

    xz_l[n_eta-1] = ZERO;
    xz_r[n_eta-1] = ZERO;

    yz_l[n_eta-1] = ZERO;
    yz_r[n_eta-1] = ZERO;

    for ( j = 1; j < n_eta - 1; ++j ) {
        xze_l[j] = HALF * (xz_l[j+1] - xz_l[j-1]);
        xze_r[j] = HALF * (xz_r[j+1] - xz_r[j-1]);

        yze_l[j] = HALF * (yz_l[j+1] - yz_l[j-1]);
        yze_r[j] = HALF * (yz_r[j+1] - yz_r[j-1]);
    } /* endfor */

    /* Begin a new iteration cycle and use the current values of xij
       and yij to determine the derivatives xee and yee at the bottom
       (south) and top (north) boundaries as well as the derivatives
       xzz and yzz at the left (west) and right (east) boundaries, 
       repectively.  These values may then be used to determine pb, 
       pt, pl, pr, qb, qt, ql, and qr.  */

    num_iter = 0;

    next_iteration: ;

    num_iter = num_iter + 1;

    if (i_poisson == 0) goto no_source_term_boundary_evaluation;

    /* BOTTOM AND TOP */

    for ( i = 1; i < n_zeta - 1; ++i ) {
        xee_b[i] = HALF * (-SEVEN*xij[i][0] + EIGHT*xij[i][1] -
                            xij[i][2]) - THREE * xe_b[i];
        xee_t[i] = HALF * (-SEVEN*xij[i][n_eta-1] + 
                            EIGHT*xij[i][n_eta-2] -
                            xij[i][n_eta-3]) + THREE * xe_t[i];

        yee_b[i] = HALF * (-SEVEN*yij[i][0] + EIGHT*yij[i][1] -
                            yij[i][2]) - THREE * ye_b[i];
        yee_t[i] = HALF * (-SEVEN*yij[i][n_eta-1] + 
                            EIGHT*yij[i][n_eta-2] -
                            yij[i][n_eta-3]) + THREE * ye_t[i];

        alpha = xe_b[i] * xe_b[i] + ye_b[i] * ye_b[i];
        beta  = xz_b[i] * xe_b[i] + yz_b[i] * ye_b[i];
        gamma = xz_b[i] * xz_b[i] + yz_b[i] * yz_b[i];
        r1 = - ONE * (alpha * xzz_b[i] - TWO * beta * xze_b[i] +
             gamma * xee_b[i]) / (jacob_b[i] * jacob_b[i]);
        r2 = - ONE * (alpha * yzz_b[i] - TWO * beta * yze_b[i] +
             gamma * yee_b[i]) / (jacob_b[i] * jacob_b[i]);

        dpb = (ye_b[i] * r1 - xe_b[i] * r2) / jacob_b[i] - pb[i];
        if (fabs(pb[i]) <= ONE ) {
          dpmax = HALF;
        } else {
          dpmax = HALF * fabs(pb[i]);
        } /* endif */
        if (wsur * fabs(dpb) <= dpmax) {
          pb[i] = pb[i] + wsur * dpb;
        } else {
	    if (dpb >= ZERO) {
               pb[i] = pb[i] + dpmax;
            } else {
               pb[i] = pb[i] - dpmax;
            } /* endif */
        } /* endif */

        dqb = (xz_b[i] * r2 - yz_b[i] * r1) / jacob_b[i] - qb[i];
        if (fabs(qb[i]) <= ONE ) {
          dqmax = HALF;
        } else {
          dqmax = HALF * fabs(qb[i]);
        } /* endif */
        if (wsur * fabs(dqb) <= dqmax) {
          qb[i] = qb[i] + wsur * dqb;
        } else {
	    if (dqb >= ZERO) {
               qb[i] = qb[i] + dqmax;
            } else {
               qb[i] = qb[i] - dqmax;
            } /* endif */
        } /* endif */

        alpha = xe_t[i] * xe_t[i] + ye_t[i] * ye_t[i];
        beta  = xz_t[i] * xe_t[i] + yz_t[i] * ye_t[i];
        gamma = xz_t[i] * xz_t[i] + yz_t[i] * yz_t[i];
        r3 = - ONE * (alpha * xzz_t[i] - TWO * beta * xze_t[i] +
             gamma * xee_t[i]) / (jacob_t[i] * jacob_t[i]);
        r4 = - ONE * (alpha * yzz_t[i] - TWO * beta * yze_t[i] +
             gamma * yee_t[i]) / (jacob_t[i] * jacob_t[i]);

        dpt = (ye_t[i] * r3 - xe_t[i] * r4) / jacob_t[i] - pt[i];
        if (fabs(pt[i]) <= ONE ) {
          dpmax = HALF;
        } else {
          dpmax = HALF * fabs(pt[i]);
        } /* endif */
        if (wsur * fabs(dpt) <= dpmax) {
          pt[i] = pt[i] + wsur * dpt;
        } else {
	    if (dpt >= ZERO) {
               pt[i] = pt[i] + dpmax;
            } else {
               pt[i] = pt[i] - dpmax;
            } /* endif */
        } /* endif */

        dqt = (xz_t[i] * r4 - yz_t[i] * r3) / jacob_t[i] - qt[i];
        if (fabs(qt[i]) <= ONE ) {
          dqmax = HALF;
        } else {
          dqmax = HALF * fabs(qt[i]);
        } /* endif */
        if (wsur * fabs(dqt) <= dqmax) {
          qt[i] = qt[i] + wsur * dqt;
        } else {
	    if (dqt >= ZERO) {
               qt[i] = qt[i] + dqmax;
            } else {
               qt[i] = qt[i] - dqmax;
            } /* endif */
        } /* endif */
    } /* endfor */

    /* LEFT AND RIGHT */

    for ( j = 1; j < n_eta - 1; ++j ) {
        xzz_l[j] = HALF * (-SEVEN*xij[0][j] + EIGHT*xij[1][j] -
                            xij[2][j]) - THREE * xz_l[j];
        xzz_r[j] = HALF * (-SEVEN*xij[n_zeta-1][j] + 
                            EIGHT*xij[n_zeta-2][j] -
                            xij[n_zeta-3][j]) + THREE * xz_r[j];

        yzz_l[j] = HALF * (-SEVEN*yij[0][j] + EIGHT*yij[1][j] -
                            yij[2][j]) - THREE * yz_l[j];
        yzz_r[j] = HALF * (-SEVEN*yij[n_zeta-1][j] + 
                            EIGHT*yij[n_zeta-2][j] -
                            yij[n_zeta-3][j]) + THREE * yz_r[j];

        alpha = xe_l[j] * xe_l[j] + ye_l[j] * ye_l[j];
        beta  = xz_l[j] * xe_l[j] + yz_l[j] * ye_l[j];
        gamma = xz_l[j] * xz_l[j] + yz_l[j] * yz_l[j];
        r1 = - ONE * (alpha * xzz_l[j] - TWO * beta * xze_l[j] +
             gamma * xee_l[j]) / (jacob_l[j] * jacob_l[j]);
        r2 = - ONE * (alpha * yzz_l[j] - TWO * beta * yze_l[j] +
             gamma * yee_l[j]) / (jacob_l[j] * jacob_l[j]);

        dpl = (ye_l[j] * r1 - xe_l[j] * r2) / jacob_l[j] - pl[j];
        if (fabs(pl[j]) <= ONE ) {
          dpmax = HALF;
        } else {
          dpmax = HALF * fabs(pl[j]);
        } /* endif */
        if (wsur * fabs(dpl) <= dpmax) {
          pl[j] = pl[j] + wsur * dpl;
        } else {
	    if (dpl >= ZERO) {
               pl[j] = pl[j] + dpmax;
            } else {
               pl[j] = pl[j] - dpmax;
            } /* endif */
        } /* endif */

        dql = (xz_l[j] * r2 - yz_l[j] * r1) / jacob_l[j] - ql[j];
        if (fabs(ql[j]) <= ONE ) {
          dqmax = HALF;
        } else {
          dqmax = HALF * fabs(ql[j]);
        } /* endif */
        if (wsur * fabs(dql) <= dqmax) {
          ql[j] = ql[j] + wsur * dql;
        } else {
	    if (dql >= ZERO) {
               ql[j] = ql[j] + dqmax;
            } else {
               ql[j] = ql[j] - dqmax;
            } /* endif */
        } /* endif */

        alpha = xe_r[j] * xe_r[j] + ye_r[j] * ye_r[j];
        beta  = xz_r[j] * xe_r[j] + yz_r[j] * ye_r[j];
        gamma = xz_r[j] * xz_r[j] + yz_r[j] * yz_r[j];
        r3 = - ONE * (alpha * xzz_r[j] - TWO * beta * xze_r[j] +
             gamma * xee_r[j]) / (jacob_r[j] * jacob_r[j]);
        r4 = - ONE * (alpha * yzz_r[j] - TWO * beta * yze_r[j] +
             gamma * yee_r[j]) / (jacob_r[j] * jacob_r[j]);

        dpr = (ye_r[j] * r3 - xe_r[j] * r4) / jacob_r[j] - pr[j];
        if (fabs(pr[j]) <= ONE ) {
          dpmax = HALF;
        } else {
          dpmax = HALF * fabs(pr[j]);
        } /* endif */
        if (wsur * fabs(dpr) <= dpmax) {
          pr[j] = pr[j] + wsur * dpr;
        } else {
	    if (dpr >= ZERO) {
               pr[j] = pr[j] + dpmax;
            } else {
               pr[j] = pr[j] - dpmax;
            } /* endif */
	} /* endif */

        dqr = (xz_r[j] * r4 - yz_r[j] * r3) / jacob_r[j] - qr[j];
        if (fabs(qr[j]) <= ONE ) {
          dqmax = HALF;
        } else {
          dqmax = HALF * fabs(qr[j]);
        } /* endif */
        if (wsur * fabs(dqr) <= dqmax) {
          qr[j] = qr[j] + wsur * dqr;
        } else {
	    if (dqr >= ZERO) {
               qr[j] = qr[j] + dqmax;
            } else {
               qr[j] = qr[j] - dqmax;
            } /* endif */
        } /* endif */
    } /* endfor */

    /* Perform one step in the SOR Gauss-Seidel centred 
       finite-difference solution technique. */

    no_source_term_boundary_evaluation: ;

    for ( j = 1; j < n_eta - 1; ++j ) {
        if (i_poisson != 0) {
            fb = exp ( - a * double(j) );
            ft = exp ( b * (double(j) - double(n_eta) + ONE) );
        } /* endfor */

        for ( i = 1; i < n_zeta - 1; ++i ) { 
            if (i_poisson != 0) {
                fl = exp ( - c * double(i) );
                fr = exp ( d * (double(i) - double(n_zeta) + ONE) );
 
                fbl = ONE -
                      exp ( -e * sqrt( double(i)*double(i) +
                      double(j)*double(j) ) );
                fbr = ONE -
                      exp ( -e * sqrt( (double(i) - double(n_zeta) + ONE) *
                      (double(i) - double(n_zeta) + ONE) +
                      double(j)*double(j) ) );
                ftl = ONE -
                      exp ( -e * sqrt( double(i)*double(i) +
                      (double(j) - double(n_eta) + ONE) *
                      (double(j) - double(n_eta) + ONE) ) );
                ftr = ONE -
                      exp ( -e * sqrt( (double(i) - double(n_zeta) + ONE) *
                      (double(i) - double(n_zeta) + ONE) +
                      (double(j) - double(n_eta) + ONE) *
                      (double(j) - double(n_eta) + ONE) ) );
                pij = (pb[i] * fb + pt[i] * ft + pl[j] * fl + pr[j] * fr) *
                      (fbl * fbr * ftl * ftr);
                qij = (qb[i] * fb + qt[i] * ft + ql[j] * fl + qr[j] * fr) *
                      (fbl * fbr * ftl * ftr);
            } else {
                pij = ZERO;
                qij = ZERO;
            } /* endif */

            xzij = HALF * (xij[i+1][j] - xij[i-1][j]);
            yzij = HALF * (yij[i+1][j] - yij[i-1][j]);
            xeij = HALF * (xij[i][j+1] - xij[i][j-1]);
            yeij = HALF * (yij[i][j+1] - yij[i][j-1]);

            jacobij = xzij * yeij - xeij * yzij;

            alpha = xeij * xeij + yeij * yeij;
            beta  = xzij * xeij + yzij * yeij;
            gamma = xzij * xzij + yzij * yzij;

            dxij = HALF * (jacobij * jacobij * (xzij * pij +
                   xeij * qij) + alpha * (xij[i+1][j] + xij[i-1][j]) +
                   gamma * (xij[i][j+1] + xij[i][j-1]) - HALF * beta *
                   (xij[i+1][j+1] - xij[i+1][j-1] - xij[i-1][j+1] +
                   xij[i-1][j-1])) / (alpha + gamma) - xij[i][j];
            xij[i][j] = xij[i][j] + wsor * dxij;
          
            dyij = HALF * (jacobij * jacobij * (yzij * pij +
                   yeij * qij) + alpha * (yij[i+1][j] + yij[i-1][j]) +
                   gamma * (yij[i][j+1] + yij[i][j-1]) - HALF * beta *
                   (yij[i+1][j+1] - yij[i+1][j-1] - yij[i-1][j+1] +
                   yij[i-1][j-1])) / (alpha + gamma) - yij[i][j];
            yij[i][j] = yij[i][j] + wsor * dyij;
        } /* endfor */
    } /* endfor */

    /* Check to see if the grid block smoothing is complete.
       If not go to next_iteration. */

    if (num_iter < Number_of_Iterations) goto next_iteration;

    /* Save the newly computed interior node locations for 
       the quadrilateral grid block. */

    for (j  = Grid.JNl ; j <= Grid.JNu ; ++j ) {
       for ( i = Grid.INl ; i <= Grid.INu ; ++i ) {
	   ii = i - Grid.INl;
	   jj = j - Grid.JNl;
           Grid.Node[i][j].X.x = xij[ii][jj];
           Grid.Node[i][j].X.y = yij[ii][jj];
       } /* endfor */
    } /* endfor */

    /* Re-compute the exterior nodes for the quadrilateral mesh block. */

    Update_Exterior_Nodes(Grid);

    /* Re-compute the cell values for the quadrilateral mesh block. */

    Update_Cells(Grid);

    /* Delete (deallocate) the Poisson equation solution variables. */

    for ( i = 0; i <= n_zeta-1 ; ++i ) {
        delete []xij[i]; xij[i] = NULL;
        delete []yij[i]; yij[i] = NULL;
    } /* endfor */
    delete []xij; xij = NULL;
    delete []yij; yij = NULL;

    delete []deltas_b; deltas_b = NULL;
    delete []deltas_t; deltas_t = NULL;
    delete []deltas_l; deltas_l = NULL; 
    delete []deltas_r; deltas_r = NULL;

    delete []xz_b; xz_b = NULL; 
    delete []yz_b; yz_b = NULL;
    delete []xe_b; xe_b = NULL;
    delete []ye_b; ye_b = NULL;
    delete []xzz_b; xzz_b = NULL;
    delete []yzz_b; yzz_b = NULL;
    delete []xee_b; xee_b = NULL;
    delete []yee_b; yee_b = NULL;
    delete []xze_b; xze_b = NULL;
    delete []yze_b; yze_b = NULL;

    delete []xz_t; xz_t = NULL;
    delete []yz_t; yz_t = NULL;
    delete []xe_t; xe_t = NULL; 
    delete []ye_t; ye_t = NULL;
    delete []xzz_t; xzz_t = NULL;
    delete []yzz_t; yzz_t = NULL;
    delete []xee_t; xee_t = NULL;
    delete []yee_t; yee_t = NULL;
    delete []xze_t; xze_t = NULL;
    delete []yze_t; yze_t = NULL;

    delete []xz_l; xz_l = NULL;
    delete []yz_l; yz_l = NULL;
    delete []xe_l; xe_l = NULL;
    delete []ye_l; ye_l = NULL;
    delete []xzz_l; xzz_l = NULL;
    delete []yzz_l; yzz_l = NULL;
    delete []xee_l; xee_l = NULL;
    delete []yee_l; yee_l = NULL;
    delete []xze_l; xze_l = NULL;
    delete []yze_l; yze_l = NULL;

    delete []xz_r; xz_r = NULL;
    delete []yz_r; yz_r = NULL;
    delete []xe_r; xe_r = NULL;
    delete []ye_r; ye_r = NULL;
    delete []xzz_r; xzz_r = NULL;
    delete []yzz_r; yzz_r = NULL;
    delete []xee_r; xee_r = NULL;
    delete []yee_r; yee_r = NULL;
    delete []xze_r; xze_r = NULL;
    delete []yze_r; yze_r = NULL;

    delete []pb; pb = NULL;
    delete []pt; pt = NULL;
    delete []pl; pl = NULL;
    delete []pr; pr = NULL;

    delete []qb; qb = NULL;
    delete []qt; qt = NULL;
    delete []ql; ql = NULL;
    delete []qr; qr = NULL;

    delete []jacob_b; jacob_b = NULL;
    delete []jacob_t; jacob_t = NULL;
    delete []jacob_l; jacob_l = NULL;
    delete []jacob_r; jacob_r = NULL;

}

/**********************************************************************
 * Routine: Smooth_Rocket_Motor                                       *
 **********************************************************************/
void Smooth_Rocket_Motor(Grid2D_Quad_Block &Grid,
			 const double &Length_Chamber,
			 const double &Radius_Chamber,
			 const double &Length_Chamber_To_Throat,
			 const double &Length_Nozzle,
			 const double &Radius_Nozzle_Exit,
			 const double &Radius_Nozzle_Throat,
			 const double &Radius_Grain,
			 const int &Nozzle_Type,
			 const double &Stretching_Factor_Idir,
			 const double &Stretching_Factor_Jdir,
			 const int &sector,
			 const int &level,
			 const int &di,const int &dj,
			 const int &ri,const int &rj,
			 const int &NRi,const int &NRj) {

  // Exit immediately if the current block is in the chamber.
  if (Grid.Node[Grid.INl+5][Grid.JNl+5].X.x < ZERO) return ;

  Smooth_Quad_Block(Grid,min(250,4*max(Grid.NCi-4,Grid.NCi-4)));

  return ;

}

/********************************************************
 * Routine: Set_BCs                                     *
 *                                                      *
 * Set the boundary condition type for the north,       *
 * south, east, and west boundaries of the              *
 * quadrilateral mesh block.                            *
 *                                                      *
 ********************************************************/
void Set_BCs(Grid2D_Quad_Block &Grid) {

    int i, j, bc_type_left, bc_type_right;
    double S_i, S_j, 
           s_north, s_south, s_east, s_west,
           smax_north, smax_south, smax_east, smax_west;

    if (Grid.BndNorthSpline.np == 0) {
       for ( i = Grid.ICl-Grid.Nghost; i <= Grid.ICu+Grid.Nghost; ++i) {
	   Grid.BCtypeN[i] = BC_NONE;
       } /* endfor */
    } else {
       Grid.BCtypeN[Grid.ICl-1] = BCtype(Grid.SminN,
                                         Grid.BndNorthSpline);
       Grid.BCtypeN[Grid.ICu+1] = BCtype(Grid.SmaxN,
                                         Grid.BndNorthSpline);
       for (int GCell=2; GCell<= Grid.Nghost; ++GCell){
	 Grid.BCtypeN[Grid.ICl-GCell] = Grid.BCtypeN[Grid.ICl-GCell+1];
	 Grid.BCtypeN[Grid.ICu+GCell] = Grid.BCtypeN[Grid.ICu+GCell-1];
       }
       for ( i = Grid.INl ; i <= Grid.INu-1 ; ++i) {
	 s_north = getS(Grid.Node[i][Grid.JNu].X, Grid.BndNorthSpline);
	 bc_type_left = BCtype(s_north, Grid.BndNorthSpline);
	 s_north = getS(Grid.Node[i+1][Grid.JNu].X, Grid.BndNorthSpline);
	 bc_type_right = BCtype(s_north, Grid.BndNorthSpline);

	 if (bc_type_left == bc_type_right) {
	   Grid.BCtypeN[i] = bc_type_left;
	 } else {
	 } /* endif */
       } /* endfor */
    } /* endif */

    if (Grid.BndSouthSpline.np == 0) {
       for ( i = Grid.ICl-Grid.Nghost; i <= Grid.ICu+Grid.Nghost ; ++i) {
	   Grid.BCtypeS[i] = BC_NONE;
       } /* endfor */
    } else {
       Grid.BCtypeS[Grid.ICl-1] = BCtype(Grid.SminS, 
                                         Grid.BndSouthSpline);
       Grid.BCtypeS[Grid.ICu+1] = BCtype(Grid.SmaxS, 
                                         Grid.BndSouthSpline);
       for (int GCell=2; GCell<= Grid.Nghost; ++GCell){
	 Grid.BCtypeS[Grid.ICl-GCell] = Grid.BCtypeS[Grid.ICl-GCell+1];
	 Grid.BCtypeS[Grid.ICu+GCell] = Grid.BCtypeS[Grid.ICu+GCell-1];
       }
       for ( i = Grid.INl ; i <= Grid.INu-1 ; ++i) {
	 s_south = getS(Grid.Node[i][Grid.JNl].X, Grid.BndSouthSpline);
	 bc_type_left = BCtype(s_south, Grid.BndSouthSpline);
	 s_south = getS(Grid.Node[i+1][Grid.JNl].X, Grid.BndSouthSpline);
	 bc_type_right = BCtype(s_south, Grid.BndSouthSpline);

	 if (bc_type_left == bc_type_right) {
	   Grid.BCtypeS[i] = bc_type_left;
	 } else {
	 } /* endif */
       } /* endfor */
    } /* endif */

    if (Grid.BndEastSpline.np == 0) {
       for ( j = Grid.JCl-Grid.Nghost; j <= Grid.JCu+Grid.Nghost ; ++j) {
	   Grid.BCtypeE[j] = BC_NONE;
       } /* endfor */
    } else {
       Grid.BCtypeE[Grid.JCl-1] = BCtype(Grid.SminE,
                                         Grid.BndEastSpline);
       Grid.BCtypeE[Grid.JCu+1] = BCtype(Grid.SmaxE, 
                                         Grid.BndEastSpline);
       for (int GCell=2; GCell<= Grid.Nghost; ++GCell){
	 Grid.BCtypeE[Grid.JCl-GCell] = Grid.BCtypeE[Grid.JCl-GCell+1];
	 Grid.BCtypeE[Grid.JCu+GCell] = Grid.BCtypeE[Grid.JCu+GCell-1];
       }
       for ( j = Grid.JNl ; j <= Grid.JNu-1 ; ++j) {
	 s_east = getS(Grid.Node[Grid.INu][j].X, Grid.BndEastSpline);
	 bc_type_left = BCtype(s_east, Grid.BndEastSpline);
	 s_east = getS(Grid.Node[Grid.INu][j+1].X, Grid.BndEastSpline);
	 bc_type_right = BCtype(s_east, Grid.BndEastSpline);

	 if (bc_type_left == bc_type_right) {
	   Grid.BCtypeE[j] = bc_type_left;
	 } else {
	 } /* endif */
       } /* endfor */
    } /* endif */

    if (Grid.BndWestSpline.np == 0) {
       for ( j = Grid.JCl-Grid.Nghost ; j <= Grid.JCu+Grid.Nghost ; ++j) {
	   Grid.BCtypeW[j] = BC_NONE;
       } /* endfor */
    } else {
       Grid.BCtypeW[Grid.JCl-1] = BCtype(Grid.SminW, 
                                         Grid.BndWestSpline);
       Grid.BCtypeW[Grid.JCu+1] = BCtype(Grid.SmaxW, 
                                         Grid.BndWestSpline);
       for (int GCell=2; GCell<= Grid.Nghost; ++GCell){
	 Grid.BCtypeW[Grid.JCl-GCell] = Grid.BCtypeW[Grid.JCl-GCell+1];
	 Grid.BCtypeW[Grid.JCu+GCell] = Grid.BCtypeW[Grid.JCu+GCell-1];
       }
       for ( j = Grid.JNl ; j <= Grid.JNu-1 ; ++j) {
	 s_west = getS(Grid.Node[Grid.INl][j].X, Grid.BndWestSpline);
	 bc_type_left = BCtype(s_west, Grid.BndWestSpline);
	 s_west = getS(Grid.Node[Grid.INl][j+1].X, Grid.BndWestSpline);
	 bc_type_right = BCtype(s_west, Grid.BndWestSpline);

	 if (bc_type_left == bc_type_right) {
	   Grid.BCtypeW[j] = bc_type_left;
	 } else {
	 } /* endif */
       } /* endfor */
    } /* endif */

}

/********************************************************
 * Routine: Update_Exterior_Nodes                       *
 *                                                      *
 * Updates the exterior nodes for the quadrilateral     *
 * mesh block.                                          *
 *                                                      *
 ********************************************************/
void Update_Exterior_Nodes(Grid2D_Quad_Block &Grid) {

    int i, j;
    Vector2D norm_dir, X_norm, X_tan;

    // Update West and East boundary nodes.
    for ( j = Grid.JNl ; j <= Grid.JNu ; ++j) {
//         if (Grid.BCtypeW[j] == BC_NONE) {
//         } else 
	  if (Grid.BCtypeW[j] == Grid.BCtypeW[j-1] &&
		   Grid.BCtypeW[j] != BC_REFLECTION &&
		   Grid.BCtypeW[j] != BC_PERIODIC &&
		   Grid.BCtypeW[j] != BC_NO_SLIP &&
		   Grid.BCtypeW[j] != BC_WALL_VISCOUS_ISOTHERMAL &&
		   Grid.BCtypeW[j] != BC_WALL_VISCOUS_HEATFLUX &&
		   Grid.BCtypeW[j] != BC_MOVING_WALL &&
		   Grid.BCtypeW[j] != BC_MOVING_WALL_ISOTHERMAL &&
		   Grid.BCtypeW[j] != BC_MOVING_WALL_HEATFLUX &&
		   Grid.BCtypeW[j] != BC_BURNING_SURFACE &&
		   Grid.BCtypeW[j] != BC_MASS_INJECTION) {
	  for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
	    Grid.Node[Grid.INl-GCell][j].X = Grid.Node[Grid.INl][j].X -
	                                    (Grid.Node[Grid.INl+GCell][j].X - 
                                             Grid.Node[Grid.INl][j].X);
	  }
        } else if (Grid.BCtypeW[j] == Grid.BCtypeW[j-1] &&
                   (Grid.BCtypeW[j] == BC_REFLECTION ||
		    Grid.BCtypeW[j] == BC_NO_SLIP ||
		    Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		    Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
		    Grid.BCtypeW[j] == BC_MOVING_WALL ||
		    Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
		    Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
		    Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
		    Grid.BCtypeW[j] == BC_MASS_INJECTION)) {
	   if (j > Grid.JNl && j < Grid.JNu) {
 	      norm_dir = - HALF*(Grid.nfaceW(Grid.ICl, j) + 
                                 Grid.nfaceW(Grid.ICl, j-1));
	      for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
		X_norm = ((Grid.Node[Grid.INl+GCell][j].X - 
			   Grid.Node[Grid.INl][j].X) * norm_dir) * norm_dir;
		X_tan = (Grid.Node[Grid.INl+GCell][j].X - 
			 Grid.Node[Grid.INl][j].X) - X_norm;
		Grid.Node[Grid.INl-GCell][j].X = Grid.Node[Grid.INl][j].X -
		                                 X_norm + X_tan;
	      }
           } else if ( j == Grid.JNl ) {
	      //norm_dir = - Grid.nfaceW(Grid.ICl, j);
	      for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
		Grid.Node[Grid.INl-GCell][j].X = Grid.Node[Grid.INl][j].X -
		                                 (Grid.Node[Grid.INl+GCell][j].X - 
						  Grid.Node[Grid.INl][j].X);
	      }
           } else {
 	      //norm_dir = - Grid.nfaceW(Grid.ICl, j-1);
	     for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
	       Grid.Node[Grid.INl-GCell][j].X = Grid.Node[Grid.INl][j].X -
                                               (Grid.Node[Grid.INl+GCell][j].X - 
						Grid.Node[Grid.INl][j].X);
	     }
           } /* endif */
        } else if (Grid.BCtypeW[j] == Grid.BCtypeW[j-1] &&
                   Grid.BCtypeW[j] == BC_PERIODIC) {
	     for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
	       Grid.Node[Grid.INl-GCell][j].X = Grid.Node[Grid.INl][j].X -
                                               (Grid.Node[Grid.INu][j].X - 
						Grid.Node[Grid.INu-GCell][j].X);
	     }
        } /* endif */

//         if (Grid.BCtypeE[j] == BC_NONE) {
//         } else 
	  if (Grid.BCtypeE[j] == Grid.BCtypeE[j-1] &&
		   Grid.BCtypeE[j] != BC_REFLECTION &&
		   Grid.BCtypeE[j] != BC_PERIODIC &&
		   Grid.BCtypeE[j] != BC_NO_SLIP &&
		   Grid.BCtypeE[j] != BC_WALL_VISCOUS_ISOTHERMAL &&
		   Grid.BCtypeE[j] != BC_WALL_VISCOUS_HEATFLUX &&
		   Grid.BCtypeE[j] != BC_MOVING_WALL &&
		   Grid.BCtypeE[j] != BC_MOVING_WALL_ISOTHERMAL &&
		   Grid.BCtypeE[j] != BC_MOVING_WALL_HEATFLUX &&
		   Grid.BCtypeE[j] != BC_BURNING_SURFACE &&
		   Grid.BCtypeE[j] != BC_MASS_INJECTION) {
	  for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
           Grid.Node[Grid.INu+GCell][j].X = Grid.Node[Grid.INu][j].X +
                                           (Grid.Node[Grid.INu][j].X - 
					    Grid.Node[Grid.INu-GCell][j].X);
	  }
        } else if (Grid.BCtypeE[j] == Grid.BCtypeE[j-1] &&
                   (Grid.BCtypeE[j] == BC_REFLECTION ||	  
		    Grid.BCtypeE[j] == BC_NO_SLIP ||
		    Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		    Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		    Grid.BCtypeE[j] == BC_MOVING_WALL ||
		    Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
		    Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
		    Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
		    Grid.BCtypeE[j] == BC_MASS_INJECTION)) {
	   if (j > Grid.JNl && j < Grid.JNu) {
 	      norm_dir = HALF*(Grid.nfaceE(Grid.ICu, j) + 
                               Grid.nfaceE(Grid.ICu, j-1));
	      for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
		X_norm = ((Grid.Node[Grid.INu][j].X - 
			   Grid.Node[Grid.INu-GCell][j].X) * norm_dir) * norm_dir;
		X_tan = (Grid.Node[Grid.INu][j].X - 
			 Grid.Node[Grid.INu-GCell][j].X) - X_norm;
		Grid.Node[Grid.INu+GCell][j].X = Grid.Node[Grid.INu][j].X +
		                                 X_norm - X_tan;
	      }
           } else if ( j == Grid.JNl ) {
	      //norm_dir = Grid.nfaceE(Grid.ICu, j);
	      for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
		Grid.Node[Grid.INu+GCell][j].X = Grid.Node[Grid.INu][j].X +
		                                 (Grid.Node[Grid.INu][j].X - 
						  Grid.Node[Grid.INu-GCell][j].X);
	      }
           } else {
	      //norm_dir = Grid.nfaceE(Grid.ICu, j-1);
	      for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
		Grid.Node[Grid.INu+GCell][j].X = Grid.Node[Grid.INu][j].X +
                                                 (Grid.Node[Grid.INu][j].X - 
						  Grid.Node[Grid.INu-GCell][j].X);
	      }
           } /* endif */
        } else if (Grid.BCtypeE[j] == Grid.BCtypeE[j-1] &&
                   Grid.BCtypeE[j] == BC_PERIODIC) {
	      for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
		Grid.Node[Grid.INu+GCell][j].X = Grid.Node[Grid.INu][j].X +
                                                 (Grid.Node[Grid.INl+GCell][j].X - 
						  Grid.Node[Grid.INl][j].X);
	      }
        } /* endif */
    } /* endfor */

    // Update South and North boundary nodes.
    for ( i = Grid.INl ; i <= Grid.INu ; ++i) {
//         if (Grid.BCtypeS[j] == BC_NONE) {
//         } else 
	  if (Grid.BCtypeS[i] == Grid.BCtypeS[i-1] &&
		   Grid.BCtypeS[i] != BC_REFLECTION &&
		   Grid.BCtypeS[i] != BC_PERIODIC &&
		   Grid.BCtypeS[i] != BC_NO_SLIP &&
		   Grid.BCtypeS[i] != BC_WALL_VISCOUS_ISOTHERMAL &&
		   Grid.BCtypeS[i] != BC_WALL_VISCOUS_HEATFLUX &&
		   Grid.BCtypeS[i] != BC_MOVING_WALL &&
		   Grid.BCtypeS[i] != BC_MOVING_WALL_ISOTHERMAL &&
		   Grid.BCtypeS[i] != BC_MOVING_WALL_HEATFLUX &&
		   Grid.BCtypeS[i] != BC_BURNING_SURFACE &&
		   Grid.BCtypeS[i] != BC_MASS_INJECTION &&
		   Grid.BCtypeS[i] != BC_RINGLEB_FLOW) {
	  for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
           Grid.Node[i][Grid.JNl-GCell].X = Grid.Node[i][Grid.JNl].X -
                                            (Grid.Node[i][Grid.JNl+GCell].X - 
					     Grid.Node[i][Grid.JNl].X);
	  }
        } else if (Grid.BCtypeS[i] == Grid.BCtypeS[i-1] &&
                   (Grid.BCtypeS[i] == BC_REFLECTION ||  
		    Grid.BCtypeS[i] == BC_NO_SLIP ||
		    Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		    Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
		    Grid.BCtypeS[i] == BC_MOVING_WALL ||
		    Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
		    Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX ||
		    Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
		    Grid.BCtypeS[i] == BC_MASS_INJECTION ||
		    Grid.BCtypeS[i] == BC_RINGLEB_FLOW)) {
	   if (i > Grid.INl && i < Grid.INu) {
//  	      norm_dir = - HALF*(Grid.nfaceS(i, Grid.JCl) + 
//                                  Grid.nfaceS(i-1, Grid.JCl));
// 	      for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
// 		X_norm = ((Grid.Node[i][Grid.JNl+GCell].X - 
// 			   Grid.Node[i][Grid.JNl].X) * norm_dir) * norm_dir;
// 		X_tan = (Grid.Node[i][Grid.JNl+GCell].X - 
// 			 Grid.Node[i][Grid.JNl].X) - X_norm;
// 		Grid.Node[i][Grid.JNl-GCell].X = Grid.Node[i][Grid.JNl].X -
// 	                                         X_norm + X_tan;
// 	      }
	     if (Grid.lfaceS(i,Grid.JCl) > NANO &&
		 Grid.lfaceS(i-1,Grid.JCl) > NANO) {
	       norm_dir = - HALF*(Grid.nfaceS(i, Grid.JCl) + 
                                  Grid.nfaceS(i-1, Grid.JCl));
	     } else if (Grid.lfaceS(i,Grid.JCl) > NANO) {
	       norm_dir = - Grid.nfaceS(i, Grid.JCl);
	     } else if (Grid.lfaceS(i,Grid.JCl) > NANO) {
	       norm_dir = - Grid.nfaceS(i-1, Grid.JCl);
	     } else {
	     }
	     for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
	       X_norm = ((Grid.Node[i][Grid.JNl+GCell].X - 
			  Grid.Node[i][Grid.JNl].X) * norm_dir) * norm_dir;
	       X_tan = (Grid.Node[i][Grid.JNl+GCell].X - 
			Grid.Node[i][Grid.JNl].X) - X_norm;
	       Grid.Node[i][Grid.JNl-GCell].X = Grid.Node[i][Grid.JNl].X -
		                                X_norm + X_tan;
	     }
           } else if ( i == Grid.INl ) {
	      //norm_dir = - Grid.nfaceS(i, Grid.JCl);
	      for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
		Grid.Node[i][Grid.JNl-GCell].X = Grid.Node[i][Grid.JNl].X -
		                                 (Grid.Node[i][Grid.JNl+GCell].X - 
						  Grid.Node[i][Grid.JNl].X);
	      }
           } else {
	      //norm_dir = - Grid.nfaceS(i-1, Grid.JCl);
	      for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
		Grid.Node[i][Grid.JNl-GCell].X = Grid.Node[i][Grid.JNl].X -
              		                         (Grid.Node[i][Grid.JNl+GCell].X - 
						  Grid.Node[i][Grid.JNl].X);
	      }
           } /* endif */
        } else if (Grid.BCtypeS[i] == Grid.BCtypeS[i-1] &&
                   Grid.BCtypeS[i] == BC_PERIODIC) {
	  for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
           Grid.Node[i][Grid.JNl-GCell].X = Grid.Node[i][Grid.JNl].X -
                                            (Grid.Node[i][Grid.JNu].X - 
					     Grid.Node[i][Grid.JNu-GCell].X);
	  }
        } /* endif */

//         if (Grid.BCtypeN[j] == BC_NONE) {
//         } else 
	  if (Grid.BCtypeN[i] == Grid.BCtypeN[i-1] &&
		   Grid.BCtypeN[i] != BC_REFLECTION &&
		   Grid.BCtypeN[i] != BC_PERIODIC &&
		   Grid.BCtypeN[i] != BC_NO_SLIP &&
		   Grid.BCtypeN[i] != BC_WALL_VISCOUS_ISOTHERMAL &&
		   Grid.BCtypeN[i] != BC_WALL_VISCOUS_HEATFLUX &&
		   Grid.BCtypeN[i] != BC_MOVING_WALL &&
		   Grid.BCtypeN[i] != BC_MOVING_WALL_ISOTHERMAL &&
		   Grid.BCtypeN[i] != BC_MOVING_WALL_HEATFLUX  &&
		   Grid.BCtypeN[i] != BC_BURNING_SURFACE &&
		   Grid.BCtypeN[i] != BC_MASS_INJECTION) {
	  for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
	    Grid.Node[i][Grid.JNu+GCell].X = Grid.Node[i][Grid.JNu].X +
                                             (Grid.Node[i][Grid.JNu].X - 
					      Grid.Node[i][Grid.JNu-GCell].X);
	  }
        } else if (Grid.BCtypeN[i] == Grid.BCtypeN[i-1] &&
                   (Grid.BCtypeN[i] == BC_REFLECTION || 
		    Grid.BCtypeN[i] == BC_NO_SLIP ||
		    Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		    Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
		    Grid.BCtypeN[i] == BC_MOVING_WALL ||
		    Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
		    Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
		    Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
		    Grid.BCtypeN[i] == BC_MASS_INJECTION)) {
	   if (i > Grid.INl && i < Grid.INu) {
 	      norm_dir = HALF*(Grid.nfaceN(i, Grid.JCu) + 
                               Grid.nfaceN(i-1, Grid.JCu));
	      for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
		X_norm = ((Grid.Node[i][Grid.JNu].X - 
			   Grid.Node[i][Grid.JNu-GCell].X) * norm_dir) * norm_dir;
		X_tan = (Grid.Node[i][Grid.JNu].X - 
			 Grid.Node[i][Grid.JNu-GCell].X) - X_norm;
		Grid.Node[i][Grid.JNu+GCell].X = Grid.Node[i][Grid.JNu].X +
		                                 X_norm - X_tan;
	      }
           } else if ( i == Grid.INl ) {
	      //norm_dir = Grid.nfaceN(i, Grid.JCu);
	      for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
		Grid.Node[i][Grid.JNu+GCell].X = Grid.Node[i][Grid.JNu].X +
                                                 (Grid.Node[i][Grid.JNu].X - 
						  Grid.Node[i][Grid.JNu-GCell].X);
	      }
           } else {
	      //norm_dir = Grid.nfaceN(i-1, Grid.JCu);
	     for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
	       Grid.Node[i][Grid.JNu+GCell].X = Grid.Node[i][Grid.JNu].X +
                                                (Grid.Node[i][Grid.JNu].X - 
						 Grid.Node[i][Grid.JNu-GCell].X);
	     }
           } /* endif */
        } else if (Grid.BCtypeN[i] == Grid.BCtypeN[i-1] &&
                   Grid.BCtypeN[i] == BC_PERIODIC) {
	  for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
           Grid.Node[i][Grid.JNu+GCell].X = Grid.Node[i][Grid.JNu].X +
                                            (Grid.Node[i][Grid.JNl+GCell].X - 
					     Grid.Node[i][Grid.JNl].X);
	  }
        } /* endif */
    } /* endfor */

    // Update SW and SE corner nodes.
    for (int j = Grid.JNl-Grid.Nghost ; j <= Grid.JNl ; ++j) {
//         if (Grid.BCtypeW[j] == BC_NONE) {
//         } else 
      if (Grid.BCtypeW[j] != BC_REFLECTION &&
		   Grid.BCtypeW[j] != BC_PERIODIC &&  
		   Grid.BCtypeW[j] != BC_NO_SLIP &&
		   Grid.BCtypeW[j] != BC_WALL_VISCOUS_ISOTHERMAL &&
		   Grid.BCtypeW[j] != BC_WALL_VISCOUS_HEATFLUX &&
		   Grid.BCtypeW[j] != BC_MOVING_WALL &&
		   Grid.BCtypeW[j] != BC_MOVING_WALL_ISOTHERMAL &&
		   Grid.BCtypeW[j] != BC_MOVING_WALL_HEATFLUX &&
		   Grid.BCtypeW[j] != BC_BURNING_SURFACE &&
		   Grid.BCtypeW[j] != BC_MASS_INJECTION) {
	  for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
           Grid.Node[Grid.INl-GCell][j].X = Grid.Node[Grid.INl][j].X -
                                            (Grid.Node[Grid.INl+GCell][j].X - 
					     Grid.Node[Grid.INl][j].X);
	  }
        } else if (Grid.BCtypeW[j] == BC_REFLECTION || 
		   Grid.BCtypeW[j] == BC_NO_SLIP ||
		   Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		   Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
		   Grid.BCtypeW[j] == BC_MOVING_WALL ||
		   Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
		   Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
		   Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
		   Grid.BCtypeW[j] == BC_MASS_INJECTION) {
 	   if (j != Grid.JNl-Grid.Nghost) {
              norm_dir = - HALF*(Grid.nfaceW(Grid.ICl, j) + 
                                 Grid.nfaceW(Grid.ICl, j-1));
           } else {
	      norm_dir = - Grid.nfaceW(Grid.ICl, j); 
           } /* endif */
	   for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
	     X_norm = ((Grid.Node[Grid.INl+GCell][j].X - 
			Grid.Node[Grid.INl][j].X) * norm_dir) * norm_dir;
	     X_tan = (Grid.Node[Grid.INl+GCell][j].X - 
		      Grid.Node[Grid.INl][j].X) - X_norm;
	     Grid.Node[Grid.INl-GCell][j].X = Grid.Node[Grid.INl][j].X -
	                                  X_norm + X_tan;
	   }
        } else if (Grid.BCtypeW[j] == BC_PERIODIC) {
	  for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
	    Grid.Node[Grid.INl-GCell][j].X = Grid.Node[Grid.INl][j].X -
                                             (Grid.Node[Grid.INu][j].X - 
					      Grid.Node[Grid.INu-GCell][j].X);
	  }
        } /* endif */

//         if (Grid.BCtypeE[j] == BC_NONE) {
//         } else 
      if (Grid.BCtypeE[j] != BC_REFLECTION &&
		   Grid.BCtypeE[j] != BC_PERIODIC &&   
		   Grid.BCtypeE[j] != BC_NO_SLIP &&
		   Grid.BCtypeE[j] != BC_WALL_VISCOUS_ISOTHERMAL &&
		   Grid.BCtypeE[j] != BC_WALL_VISCOUS_HEATFLUX &&
		   Grid.BCtypeE[j] != BC_MOVING_WALL &&
		   Grid.BCtypeE[j] != BC_MOVING_WALL_ISOTHERMAL &&
		   Grid.BCtypeE[j] != BC_MOVING_WALL_HEATFLUX &&
		   Grid.BCtypeE[j] != BC_BURNING_SURFACE &&
		   Grid.BCtypeE[j] != BC_MASS_INJECTION) {
	  for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
	    Grid.Node[Grid.INu+GCell][j].X = Grid.Node[Grid.INu][j].X +
                                             (Grid.Node[Grid.INu][j].X - 
					      Grid.Node[Grid.INu-GCell][j].X);
	  }
        } else if (Grid.BCtypeE[j] == BC_REFLECTION ||
		   Grid.BCtypeE[j] == BC_NO_SLIP ||
		   Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		   Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		   Grid.BCtypeE[j] == BC_MOVING_WALL ||
		   Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
		   Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX  ||
		   Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
		   Grid.BCtypeE[j] == BC_MASS_INJECTION) {
 	   if (j != Grid.JNl-Grid.Nghost) {
              norm_dir = HALF*(Grid.nfaceE(Grid.ICu, j) + 
                               Grid.nfaceE(Grid.ICu, j-1));
           } else {
	      norm_dir = Grid.nfaceE(Grid.ICu, j); 
           } /* endif */
	   for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
	     X_norm = ((Grid.Node[Grid.INu][j].X - 
			Grid.Node[Grid.INu-GCell][j].X) * norm_dir) * norm_dir;
	     X_tan = (Grid.Node[Grid.INu][j].X - 
		      Grid.Node[Grid.INu-GCell][j].X) - X_norm;
	     Grid.Node[Grid.INu+GCell][j].X = Grid.Node[Grid.INu][j].X +
	                                      X_norm - X_tan;
	   }
        } else if (Grid.BCtypeE[j] == BC_PERIODIC) {
	   for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
	     Grid.Node[Grid.INu+GCell][j].X = Grid.Node[Grid.INu][j].X +
                                              (Grid.Node[Grid.INl+GCell][j].X - 
					       Grid.Node[Grid.INl][j].X);
	   }
        } /* endif */
    } /* endfor */

    // Update NW and NE corner nodes.
    for (int j = Grid.JNu ; j <= Grid.JNu+Grid.Nghost ; ++j) {
//         if (Grid.BCtypeW[j] == BC_NONE) {
//         } else 
      if (Grid.BCtypeW[j-1] != BC_REFLECTION &&
		   Grid.BCtypeW[j-1] != BC_PERIODIC &&  
		   Grid.BCtypeW[j-1] != BC_NO_SLIP &&
		   Grid.BCtypeW[j-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
		   Grid.BCtypeW[j-1] != BC_WALL_VISCOUS_HEATFLUX &&
		   Grid.BCtypeW[j-1] != BC_MOVING_WALL &&
		   Grid.BCtypeW[j-1] != BC_MOVING_WALL_ISOTHERMAL &&
		   Grid.BCtypeW[j-1] != BC_MOVING_WALL_HEATFLUX &&
		   Grid.BCtypeW[j-1] != BC_BURNING_SURFACE &&
		   Grid.BCtypeW[j-1] != BC_MASS_INJECTION) {
	  for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
	    Grid.Node[Grid.INl-GCell][j].X = Grid.Node[Grid.INl][j].X -
                                            (Grid.Node[Grid.INl+GCell][j].X - 
					     Grid.Node[Grid.INl][j].X);
	  }
        } else if (Grid.BCtypeW[j-1] == BC_REFLECTION || 
		   Grid.BCtypeW[j-1] == BC_NO_SLIP ||
		   Grid.BCtypeW[j-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
		   Grid.BCtypeW[j-1] == BC_WALL_VISCOUS_HEATFLUX ||
		   Grid.BCtypeW[j-1] == BC_MOVING_WALL ||
		   Grid.BCtypeW[j-1] == BC_MOVING_WALL_ISOTHERMAL ||
		   Grid.BCtypeW[j-1] == BC_MOVING_WALL_HEATFLUX ||
		   Grid.BCtypeW[j-1] == BC_BURNING_SURFACE ||
		   Grid.BCtypeW[j-1] == BC_MASS_INJECTION) {
 	   if (j != Grid.JNu+Grid.Nghost) {
              norm_dir = - HALF*(Grid.nfaceW(Grid.ICl, j) + 
                                 Grid.nfaceW(Grid.ICl, j-1));
           } else {
	      norm_dir = - Grid.nfaceW(Grid.ICl, j-1); 
           } /* endif */
	   for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
	     X_norm = ((Grid.Node[Grid.INl+GCell][j].X - 
			Grid.Node[Grid.INl][j].X) * norm_dir) * norm_dir;
	     X_tan = (Grid.Node[Grid.INl+GCell][j].X - 
		      Grid.Node[Grid.INl][j].X) - X_norm;
	     Grid.Node[Grid.INl-GCell][j].X = Grid.Node[Grid.INl][j].X -
	                                      X_norm + X_tan;
	   }
        } else if (Grid.BCtypeW[j-1] == BC_PERIODIC) {
	  for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
	    Grid.Node[Grid.INl-GCell][j].X = Grid.Node[Grid.INl][j].X -
                                             (Grid.Node[Grid.INu][j].X - 
					      Grid.Node[Grid.INu-GCell][j].X);
	  }
        } /* endif */

//         if (Grid.BCtypeE[j] == BC_NONE) {
//         } else 
      if (Grid.BCtypeE[j-1] != BC_REFLECTION &&
		   Grid.BCtypeE[j-1] != BC_PERIODIC && 
		   Grid.BCtypeE[j-1] != BC_NO_SLIP &&
		   Grid.BCtypeE[j-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
		   Grid.BCtypeE[j-1] != BC_WALL_VISCOUS_HEATFLUX &&
		   Grid.BCtypeE[j-1] != BC_MOVING_WALL &&
		   Grid.BCtypeE[j-1] != BC_MOVING_WALL_ISOTHERMAL &&
		   Grid.BCtypeE[j-1] != BC_MOVING_WALL_HEATFLUX &&
		   Grid.BCtypeE[j-1] != BC_BURNING_SURFACE &&
		   Grid.BCtypeE[j-1] != BC_MASS_INJECTION) {
	  for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
	    Grid.Node[Grid.INu+GCell][j].X = Grid.Node[Grid.INu][j].X +
                                             (Grid.Node[Grid.INu][j].X - 
					      Grid.Node[Grid.INu-GCell][j].X);
	  }
        } else if (Grid.BCtypeE[j-1] == BC_REFLECTION || 
		   Grid.BCtypeE[j-1] == BC_NO_SLIP ||
		   Grid.BCtypeE[j-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
		   Grid.BCtypeE[j-1] == BC_WALL_VISCOUS_HEATFLUX ||
		   Grid.BCtypeE[j-1] == BC_MOVING_WALL ||
		   Grid.BCtypeE[j-1] == BC_MOVING_WALL_ISOTHERMAL ||
		   Grid.BCtypeE[j-1] == BC_MOVING_WALL_HEATFLUX ||
		   Grid.BCtypeE[j-1] == BC_BURNING_SURFACE ||
		   Grid.BCtypeE[j-1] == BC_MASS_INJECTION) {
 	   if (j != Grid.JNu+Grid.Nghost) {
              norm_dir = HALF*(Grid.nfaceE(Grid.ICu, j) + 
                               Grid.nfaceE(Grid.ICu, j-1));
           } else {
	      norm_dir = Grid.nfaceE(Grid.ICu, j-1); 
           } /* endif */
	   for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
	     X_norm = ((Grid.Node[Grid.INu][j].X - 
			Grid.Node[Grid.INu-GCell][j].X) * norm_dir) * norm_dir;
	     X_tan = (Grid.Node[Grid.INu][j].X - 
		      Grid.Node[Grid.INu-GCell][j].X) - X_norm;
	     Grid.Node[Grid.INu+GCell][j].X = Grid.Node[Grid.INu][j].X +
	                                      X_norm - X_tan;
	   }
        } else if (Grid.BCtypeE[j-1] == BC_PERIODIC) {
	  for(int GCell=1; GCell<=Grid.Nghost; ++GCell){
           Grid.Node[Grid.INu+GCell][j].X = Grid.Node[Grid.INu][j].X +
                                            (Grid.Node[Grid.INl+GCell][j].X - 
					     Grid.Node[Grid.INl][j].X);
	  }
        } /* endif */
    } /* endfor */

    Update_Corner_Ghost_Nodes(Grid);

}

/********************************************************
 * Routine: Update_Corner_Ghost_Nodes                   *
 *                                                      *
 * Updates the corner ghost nodes for the quadrilateral *
 * mesh block.                                          *
 *                                                      *
 ********************************************************/
void Update_Corner_Ghost_Nodes(Grid2D_Quad_Block &Grid) {

  Vector2D norm_dir, X_norm, X_tan;

  // SOUTH-WEST corner:
  if (Grid.BCtypeS[Grid.INl-1] == BC_NONE &&
      Grid.BCtypeW[Grid.JNl-1] == BC_NONE) {
    // Do nothing.

  } else if ((Grid.BCtypeS[Grid.INl-1] != BC_REFLECTION &&
	      Grid.BCtypeS[Grid.INl-1] != BC_PERIODIC &&
	      Grid.BCtypeS[Grid.INl-1] != BC_NO_SLIP &&
	      Grid.BCtypeS[Grid.INl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      Grid.BCtypeS[Grid.INl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	      Grid.BCtypeS[Grid.INl-1] != BC_MOVING_WALL &&
	      Grid.BCtypeS[Grid.INl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	      Grid.BCtypeS[Grid.INl-1] != BC_MOVING_WALL_HEATFLUX &&
	      Grid.BCtypeS[Grid.INl-1] != BC_BURNING_SURFACE &&
	      Grid.BCtypeS[Grid.INl-1] != BC_MASS_INJECTION &&
	      Grid.BCtypeS[Grid.INl-1] != BC_RINGLEB_FLOW) &&
	     Grid.BCtypeW[Grid.JNl-1] == BC_NONE) {
    // Extrapolate cells south.
    for (int ng = 1; ng <= Grid.Nghost; ng++) {
      for (int i = Grid.INl-Grid.Nghost; i <= Grid.INl; i++) {
	Grid.Node[i][Grid.JNl-ng].X = Grid.Node[i][Grid.JNl].X +
                                      (Grid.Node[i][Grid.JNl].X - 
				       Grid.Node[i][Grid.JNl+ng].X);
      }
    }

  } else if ((Grid.BCtypeS[Grid.INl-1] == BC_NONE ||
	      (Grid.BCtypeS[Grid.INl-1] != BC_REFLECTION &&
	       Grid.BCtypeS[Grid.INl-1] != BC_PERIODIC &&
	       Grid.BCtypeS[Grid.INl-1] != BC_NO_SLIP &&
	       Grid.BCtypeS[Grid.INl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       Grid.BCtypeS[Grid.INl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	       Grid.BCtypeS[Grid.INl-1] != BC_MOVING_WALL &&
	       Grid.BCtypeS[Grid.INl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	       Grid.BCtypeS[Grid.INl-1] != BC_MOVING_WALL_HEATFLUX &&
	       Grid.BCtypeS[Grid.INl-1] != BC_BURNING_SURFACE &&
	       Grid.BCtypeS[Grid.INl-1] != BC_MASS_INJECTION &&
	       Grid.BCtypeS[Grid.INl-1] != BC_RINGLEB_FLOW)) &&
	     (Grid.BCtypeW[Grid.JNl-1] != BC_REFLECTION &&
	      Grid.BCtypeW[Grid.JNl-1] != BC_PERIODIC &&  
	      Grid.BCtypeW[Grid.JNl-1] != BC_NO_SLIP &&
	      Grid.BCtypeW[Grid.JNl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      Grid.BCtypeW[Grid.JNl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	      Grid.BCtypeW[Grid.JNl-1] != BC_MOVING_WALL &&
	      Grid.BCtypeW[Grid.JNl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	      Grid.BCtypeW[Grid.JNl-1] != BC_MOVING_WALL_HEATFLUX &&
	      Grid.BCtypeW[Grid.JNl-1] != BC_BURNING_SURFACE &&
	      Grid.BCtypeW[Grid.JNl-1] != BC_MASS_INJECTION)) {
    // Extrapolate cells west.
    for (int ng = 1; ng <= Grid.Nghost; ng++) {
      for (int j = Grid.JNl-Grid.Nghost; j <= Grid.JNl; j++) {
 	Grid.Node[Grid.INl-ng][j].X = Grid.Node[Grid.INl][j].X -
 	                              (Grid.Node[Grid.INl+ng][j].X - 
 				       Grid.Node[Grid.INl][j].X);
      }
    }

  } else if ((Grid.BCtypeS[Grid.INl-1] == BC_NONE ||
	      (Grid.BCtypeS[Grid.INl-1] != BC_REFLECTION &&
	       Grid.BCtypeS[Grid.INl-1] != BC_PERIODIC &&
	       Grid.BCtypeS[Grid.INl-1] != BC_NO_SLIP &&
	       Grid.BCtypeS[Grid.INl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       Grid.BCtypeS[Grid.INl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	       Grid.BCtypeS[Grid.INl-1] != BC_MOVING_WALL &&
	       Grid.BCtypeS[Grid.INl-1] != BC_MOVING_WALL_HEATFLUX &&
	       Grid.BCtypeS[Grid.INl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	       Grid.BCtypeS[Grid.INl-1] != BC_BURNING_SURFACE &&
	       Grid.BCtypeS[Grid.INl-1] != BC_MASS_INJECTION &&
	       Grid.BCtypeS[Grid.INl-1] != BC_RINGLEB_FLOW)) &&
	     (Grid.BCtypeW[Grid.JNl-1] == BC_REFLECTION ||
	      Grid.BCtypeW[Grid.JNl-1] == BC_PERIODIC ||
	      Grid.BCtypeW[Grid.JNl-1] == BC_NO_SLIP ||
	      Grid.BCtypeW[Grid.JNl-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	      Grid.BCtypeW[Grid.JNl-1] == BC_WALL_VISCOUS_HEATFLUX ||
	      Grid.BCtypeW[Grid.JNl-1] == BC_MOVING_WALL ||
	      Grid.BCtypeW[Grid.JNl-1] == BC_MOVING_WALL_ISOTHERMAL ||
	      Grid.BCtypeW[Grid.JNl-1] == BC_MOVING_WALL_HEATFLUX ||
	      Grid.BCtypeW[Grid.JNl-1] == BC_BURNING_SURFACE ||
	      Grid.BCtypeW[Grid.JNl-1] == BC_MASS_INJECTION)) {
    // Reflect cells west.
    for (int ng = 1; ng <= Grid.Nghost; ng++) {
      for (int j = Grid.JNl-Grid.Nghost; j <= Grid.JNl; j++) {
	if (j == Grid.JNl-Grid.Nghost) {
	  norm_dir = - Grid.nfaceW(Grid.ICl,Grid.JCl-Grid.Nghost);
	} else {
	  norm_dir = - HALF*(Grid.nfaceW(Grid.ICl,j) + 
			     Grid.nfaceW(Grid.ICl,j-1));
	}
	X_norm = ((Grid.Node[Grid.INl+ng][j].X - 
		   Grid.Node[Grid.INl][j].X)*norm_dir)*norm_dir;
	X_tan = (Grid.Node[Grid.INl+ng][j].X -
		 Grid.Node[Grid.INl][j].X) - X_norm;
	Grid.Node[Grid.INl-ng][j].X = Grid.Node[Grid.INl][j].X -
  	                              X_norm + X_tan;
      }
    }

  } else if (Grid.BCtypeS[Grid.INl-1] == BC_REFLECTION ||
	     Grid.BCtypeS[Grid.INl-1] == BC_PERIODIC ||
	     Grid.BCtypeS[Grid.INl-1] == BC_NO_SLIP ||
	     Grid.BCtypeS[Grid.INl-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     Grid.BCtypeS[Grid.INl-1] == BC_WALL_VISCOUS_HEATFLUX ||
	     Grid.BCtypeS[Grid.INl-1] == BC_MOVING_WALL ||
	     Grid.BCtypeS[Grid.INl-1] == BC_MOVING_WALL_ISOTHERMAL ||
	     Grid.BCtypeS[Grid.INl-1] == BC_MOVING_WALL_HEATFLUX ||
	     Grid.BCtypeS[Grid.INl-1] == BC_BURNING_SURFACE ||
	     Grid.BCtypeS[Grid.INl-1] == BC_MASS_INJECTION ||
	     Grid.BCtypeS[Grid.INl-1] == BC_RINGLEB_FLOW) {
    // Reflect cells south.
    for (int ng = 1; ng <= Grid.Nghost; ng++) {
      //for (int i = Grid.INl-Grid.Nghost; i <= Grid.INl; i++) {
      for (int i = Grid.INl; i >= Grid.INl-Grid.Nghost; i--) {
	if (i == Grid.INl-Grid.Nghost) {
	  norm_dir = - Grid.nfaceS(Grid.ICl-Grid.Nghost,Grid.JCl);
	} else {
	  norm_dir = - HALF*(Grid.nfaceS(i,Grid.JCl) +
			     Grid.nfaceS(i-1,Grid.JCl));
	}
	X_norm = ((Grid.Node[i][Grid.JNl+ng].X - 
		   Grid.Node[i][Grid.JNl].X)*norm_dir)*norm_dir;
	X_tan = (Grid.Node[i][Grid.JNl+ng].X - 
		 Grid.Node[i][Grid.JNl].X) - X_norm;
	Grid.Node[i][Grid.JNl-ng].X = Grid.Node[i][Grid.JNl].X -
	                              X_norm + X_tan;
      }
    }
  }

  // SOUTH-EAST corner:
  if (Grid.BCtypeS[Grid.INu+1] == BC_NONE &&
      Grid.BCtypeE[Grid.JNl-1] == BC_NONE) {
    // Do nothing.

  } else if ((Grid.BCtypeS[Grid.INu+1] != BC_REFLECTION &&
	      Grid.BCtypeS[Grid.INu+1] != BC_PERIODIC &&
	      Grid.BCtypeS[Grid.INu+1] != BC_NO_SLIP &&
	      Grid.BCtypeS[Grid.INu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      Grid.BCtypeS[Grid.INu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	      Grid.BCtypeS[Grid.INu+1] != BC_MOVING_WALL &&
	      Grid.BCtypeS[Grid.INu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	      Grid.BCtypeS[Grid.INu+1] != BC_MOVING_WALL_HEATFLUX &&
	      Grid.BCtypeS[Grid.INu+1] != BC_BURNING_SURFACE &&
	      Grid.BCtypeS[Grid.INu+1] != BC_MASS_INJECTION &&
	      Grid.BCtypeS[Grid.INu+1] != BC_RINGLEB_FLOW) &&
	     Grid.BCtypeE[Grid.JNl-1] == BC_NONE) {
    // Extrapolate cells south.
    for (int ng = 1; ng <= Grid.Nghost; ng++) {
      for (int i = Grid.INu; i <= Grid.INu+Grid.Nghost; i++) {
	Grid.Node[i][Grid.JNl-ng].X = Grid.Node[i][Grid.JNl].X +
                                      (Grid.Node[i][Grid.JNl].X - 
				       Grid.Node[i][Grid.JNl+ng].X);
      }
    }

  } else if ((Grid.BCtypeS[Grid.INu-1] == BC_NONE ||
	      (Grid.BCtypeS[Grid.INu+1] != BC_REFLECTION &&
	       Grid.BCtypeS[Grid.INu+1] != BC_PERIODIC &&
	       Grid.BCtypeS[Grid.INu+1] != BC_NO_SLIP &&
	       Grid.BCtypeS[Grid.INu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       Grid.BCtypeS[Grid.INu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	       Grid.BCtypeS[Grid.INu+1] != BC_MOVING_WALL &&
	       Grid.BCtypeS[Grid.INu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	       Grid.BCtypeS[Grid.INu+1] != BC_MOVING_WALL_HEATFLUX &&
	       Grid.BCtypeS[Grid.INu+1] != BC_BURNING_SURFACE &&
	       Grid.BCtypeS[Grid.INu+1] != BC_MASS_INJECTION &&
	       Grid.BCtypeS[Grid.INu+1] != BC_RINGLEB_FLOW)) &&
	     (Grid.BCtypeE[Grid.JNl-1] != BC_REFLECTION &&
	      Grid.BCtypeE[Grid.JNl-1] != BC_PERIODIC &&  
	      Grid.BCtypeE[Grid.JNl-1] != BC_NO_SLIP &&
	      Grid.BCtypeE[Grid.JNl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      Grid.BCtypeE[Grid.JNl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	      Grid.BCtypeE[Grid.JNl-1] != BC_MOVING_WALL &&
	      Grid.BCtypeE[Grid.JNl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	      Grid.BCtypeE[Grid.JNl-1] != BC_MOVING_WALL_HEATFLUX &&
	      Grid.BCtypeE[Grid.JNl-1] != BC_BURNING_SURFACE &&
	      Grid.BCtypeE[Grid.JNl-1] != BC_MASS_INJECTION)) {
    // Extrapolate cells east.
    for (int ng = 1; ng <= Grid.Nghost; ng++) {
      for (int j = Grid.JNl-Grid.Nghost; j <= Grid.JNl; j++) {
	Grid.Node[Grid.INu+ng][j].X = Grid.Node[Grid.INu][j].X +
	                              (Grid.Node[Grid.INu][j].X - 
				       Grid.Node[Grid.INu-ng][j].X);
      }
    }

  } else if ((Grid.BCtypeS[Grid.INu+1] == BC_NONE ||
	      (Grid.BCtypeS[Grid.INu+1] != BC_REFLECTION &&
	       Grid.BCtypeS[Grid.INu+1] != BC_PERIODIC &&
	       Grid.BCtypeS[Grid.INu+1] != BC_NO_SLIP &&
	       Grid.BCtypeS[Grid.INu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       Grid.BCtypeS[Grid.INu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	       Grid.BCtypeS[Grid.INu+1] != BC_MOVING_WALL &&
	       Grid.BCtypeS[Grid.INu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	       Grid.BCtypeS[Grid.INu+1] != BC_MOVING_WALL_HEATFLUX &&
	       Grid.BCtypeS[Grid.INu+1] != BC_BURNING_SURFACE &&
	       Grid.BCtypeS[Grid.INu+1] != BC_MASS_INJECTION &&
	       Grid.BCtypeS[Grid.INu+1] != BC_RINGLEB_FLOW)) &&
	     (Grid.BCtypeE[Grid.JNl-1] == BC_REFLECTION ||
	      Grid.BCtypeE[Grid.JNl-1] == BC_PERIODIC ||
	      Grid.BCtypeE[Grid.JNl-1] == BC_NO_SLIP ||
	      Grid.BCtypeE[Grid.JNl-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	      Grid.BCtypeE[Grid.JNl-1] == BC_WALL_VISCOUS_HEATFLUX ||
	      Grid.BCtypeE[Grid.JNl-1] == BC_MOVING_WALL ||
	      Grid.BCtypeE[Grid.JNl-1] == BC_MOVING_WALL_ISOTHERMAL ||
	      Grid.BCtypeE[Grid.JNl-1] == BC_MOVING_WALL_HEATFLUX ||
	      Grid.BCtypeE[Grid.JNl-1] == BC_BURNING_SURFACE ||
	      Grid.BCtypeE[Grid.JNl-1] == BC_MASS_INJECTION)) {
    // Reflect cells east.
    for (int ng = 1; ng <= Grid.Nghost; ng++) {
      for (int j = Grid.JNl; j >= Grid.JNl-Grid.Nghost; j--) {
	if (j == Grid.JNl-Grid.Nghost) {
	  norm_dir = Grid.nfaceE(Grid.ICu,Grid.JCl-Grid.Nghost);
	} else {
	  norm_dir = HALF*(Grid.nfaceE(Grid.ICu,j) + 
			   Grid.nfaceE(Grid.ICu,j+1));
	}
	norm_dir = Grid.nfaceE(Grid.ICu,j);
	X_norm = ((Grid.Node[Grid.INu][j].X - 
		   Grid.Node[Grid.INu-ng][j].X)*norm_dir)*norm_dir;
	X_tan = (Grid.Node[Grid.INu][j].X - 
		 Grid.Node[Grid.INu-ng][j].X) - X_norm;
	Grid.Node[Grid.INu+ng][j].X = Grid.Node[Grid.INu][j].X +
  	                              X_norm - X_tan;
      }
    }

  } else if (Grid.BCtypeS[Grid.INu+1] == BC_REFLECTION ||
	     Grid.BCtypeS[Grid.INu+1] == BC_PERIODIC ||
	     Grid.BCtypeS[Grid.INu+1] == BC_NO_SLIP ||
	     Grid.BCtypeS[Grid.INu+1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     Grid.BCtypeS[Grid.INu+1] == BC_WALL_VISCOUS_HEATFLUX ||
	     Grid.BCtypeS[Grid.INu+1] == BC_MOVING_WALL ||
	     Grid.BCtypeS[Grid.INu+1] == BC_MOVING_WALL_ISOTHERMAL ||
	     Grid.BCtypeS[Grid.INu+1] == BC_MOVING_WALL_HEATFLUX ||
	     Grid.BCtypeS[Grid.INu+1] == BC_BURNING_SURFACE ||
	     Grid.BCtypeS[Grid.INu+1] == BC_MASS_INJECTION ||
	     Grid.BCtypeS[Grid.INu+1] == BC_RINGLEB_FLOW) {
    // Reflect cells south.
    for (int ng = 1; ng <= Grid.Nghost; ng++) {
      for (int i = Grid.INu; i <= Grid.INu+Grid.Nghost; i++) {
	if (i == Grid.INu+Grid.Nghost) {
	  norm_dir = -Grid.nfaceS(Grid.ICu+Grid.Nghost,Grid.JCl);
	} else {
	  norm_dir = -HALF*(Grid.nfaceS(i,Grid.JCl) +
			    Grid.nfaceS(i-1,Grid.JCl));
	}
 	X_norm = ((Grid.Node[i][Grid.JNl+ng].X - 
 		   Grid.Node[i][Grid.JNl].X)*norm_dir)*norm_dir;
 	X_tan = (Grid.Node[i][Grid.JNl+ng].X - 
 		 Grid.Node[i][Grid.JNl].X) - X_norm;
  	Grid.Node[i][Grid.JNl-ng].X = Grid.Node[i][Grid.JNl].X -
  	                              X_norm + X_tan;
      }
    }
  }

  // NORTH-WEST corner:
  if (Grid.BCtypeN[Grid.INl-1] == BC_NONE &&
      Grid.BCtypeW[Grid.JNu+1] == BC_NONE) {
    // Do nothing.

  } else if ((Grid.BCtypeN[Grid.INl-1] != BC_REFLECTION &&
	      Grid.BCtypeN[Grid.INl-1] != BC_PERIODIC &&
	      Grid.BCtypeN[Grid.INl-1] != BC_NO_SLIP &&
	      Grid.BCtypeN[Grid.INl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      Grid.BCtypeN[Grid.INl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	      Grid.BCtypeN[Grid.INl-1] != BC_MOVING_WALL &&
	      Grid.BCtypeN[Grid.INl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	      Grid.BCtypeN[Grid.INl-1] != BC_MOVING_WALL_HEATFLUX &&
	      Grid.BCtypeN[Grid.INl-1] != BC_BURNING_SURFACE &&
	      Grid.BCtypeN[Grid.INl-1] != BC_MASS_INJECTION) &&
	     Grid.BCtypeW[Grid.JNu+1] == BC_NONE) {
    // Extrapolate cells north.
    for (int ng = 1; ng <= Grid.Nghost; ng++) {
      for (int i = Grid.INl-Grid.Nghost; i <= Grid.INl; i++) {
	Grid.Node[i][Grid.JNu+ng].X = Grid.Node[i][Grid.JNu].X +
                                      (Grid.Node[i][Grid.JNu].X - 
				       Grid.Node[i][Grid.JNu-ng].X);
      }
    }

  } else if ((Grid.BCtypeN[Grid.INl-1] == BC_NONE ||
	      (Grid.BCtypeN[Grid.INl-1] != BC_REFLECTION &&
	       Grid.BCtypeN[Grid.INl-1] != BC_PERIODIC &&
	       Grid.BCtypeN[Grid.INl-1] != BC_NO_SLIP &&
	       Grid.BCtypeN[Grid.INl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       Grid.BCtypeN[Grid.INl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	       Grid.BCtypeN[Grid.INl-1] != BC_MOVING_WALL &&
	       Grid.BCtypeN[Grid.INl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	       Grid.BCtypeN[Grid.INl-1] != BC_MOVING_WALL_HEATFLUX &&
	       Grid.BCtypeN[Grid.INl-1] != BC_BURNING_SURFACE &&
	       Grid.BCtypeN[Grid.INl-1] != BC_MASS_INJECTION)) &&
	     (Grid.BCtypeW[Grid.JNu+1] != BC_REFLECTION &&
	      Grid.BCtypeW[Grid.JNu+1] != BC_PERIODIC &&  
	      Grid.BCtypeW[Grid.JNu+1] != BC_NO_SLIP &&
	      Grid.BCtypeW[Grid.JNu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      Grid.BCtypeW[Grid.JNu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	      Grid.BCtypeW[Grid.JNu+1] != BC_MOVING_WALL &&
	      Grid.BCtypeW[Grid.JNu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	      Grid.BCtypeW[Grid.JNu+1] != BC_MOVING_WALL_HEATFLUX &&
	      Grid.BCtypeW[Grid.JNu+1] != BC_BURNING_SURFACE &&
	      Grid.BCtypeW[Grid.JNu+1] != BC_MASS_INJECTION)) {
    // Extrapolate cells west.
    for (int ng = 1; ng <= Grid.Nghost; ng++) {
      for (int j = Grid.JNu; j <= Grid.JNu+Grid.Nghost; j++) {
	Grid.Node[Grid.INl-ng][j].X = Grid.Node[Grid.INl][j].X -
	                              (Grid.Node[Grid.INl+ng][j].X - 
				       Grid.Node[Grid.INl][j].X);
      }
    }

  } else if ((Grid.BCtypeN[Grid.INl-1] == BC_NONE ||
	      (Grid.BCtypeN[Grid.INl-1] != BC_REFLECTION &&
	       Grid.BCtypeN[Grid.INl-1] != BC_PERIODIC &&
	       Grid.BCtypeN[Grid.INl-1] != BC_NO_SLIP &&
	       Grid.BCtypeN[Grid.INl-1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       Grid.BCtypeN[Grid.INl-1] != BC_WALL_VISCOUS_HEATFLUX &&
	       Grid.BCtypeN[Grid.INl-1] != BC_MOVING_WALL &&
	       Grid.BCtypeN[Grid.INl-1] != BC_MOVING_WALL_ISOTHERMAL &&
	       Grid.BCtypeN[Grid.INl-1] != BC_MOVING_WALL_HEATFLUX &&
	       Grid.BCtypeN[Grid.INl-1] != BC_BURNING_SURFACE &&
	       Grid.BCtypeN[Grid.INl-1] != BC_MASS_INJECTION)) &&
	     (Grid.BCtypeW[Grid.JNu+1] == BC_REFLECTION ||
	      Grid.BCtypeW[Grid.JNu+1] == BC_PERIODIC ||
	      Grid.BCtypeW[Grid.JNu+1] == BC_NO_SLIP ||
	      Grid.BCtypeW[Grid.JNu+1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	      Grid.BCtypeW[Grid.JNu+1] == BC_WALL_VISCOUS_HEATFLUX ||
	      Grid.BCtypeW[Grid.JNu+1] == BC_MOVING_WALL ||
	      Grid.BCtypeW[Grid.JNu+1] == BC_MOVING_WALL_ISOTHERMAL ||
	      Grid.BCtypeW[Grid.JNu+1] == BC_MOVING_WALL_HEATFLUX ||
	      Grid.BCtypeW[Grid.JNu+1] == BC_BURNING_SURFACE ||
	      Grid.BCtypeW[Grid.JNu+1] == BC_MASS_INJECTION)) {
    // Reflect cells west.
    for (int ng = 1; ng <= Grid.Nghost; ng++) {
      for (int j = Grid.JNu; j <= Grid.JNu+Grid.Nghost; j++) {
	if (j == Grid.JNu+Grid.Nghost) {
	  norm_dir = - Grid.nfaceW(Grid.ICl,Grid.JCu+Grid.Nghost);
	} else {
	  norm_dir = - HALF*(Grid.nfaceW(Grid.ICl,j) + 
			     Grid.nfaceW(Grid.ICl,j-1));
	}
	X_norm = ((Grid.Node[Grid.INl+ng][j].X -
		   Grid.Node[Grid.INl][j].X)*norm_dir)*norm_dir;
	X_tan = (Grid.Node[Grid.INl+ng][j].X - 
		 Grid.Node[Grid.INl][j].X) - X_norm;
	Grid.Node[Grid.INl-ng][j].X = Grid.Node[Grid.INl][j].X -
  	                              X_norm + X_tan;
      }
    }

  } else if (Grid.BCtypeN[Grid.INl-1] == BC_REFLECTION ||
	     Grid.BCtypeN[Grid.INl-1] == BC_PERIODIC ||
	     Grid.BCtypeN[Grid.INl-1] == BC_NO_SLIP ||
	     Grid.BCtypeN[Grid.INl-1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     Grid.BCtypeN[Grid.INl-1] == BC_WALL_VISCOUS_HEATFLUX ||
	     Grid.BCtypeN[Grid.INl-1] == BC_MOVING_WALL ||
	     Grid.BCtypeN[Grid.INl-1] == BC_MOVING_WALL_ISOTHERMAL ||
	     Grid.BCtypeN[Grid.INl-1] == BC_MOVING_WALL_HEATFLUX ||
	     Grid.BCtypeN[Grid.INl-1] == BC_BURNING_SURFACE ||
	     Grid.BCtypeN[Grid.INl-1] == BC_MASS_INJECTION) {
    // Reflect cells north.
    for (int ng = 1; ng <= Grid.Nghost; ng++) {
      for (int i = Grid.INl-Grid.Nghost; i <= Grid.INl; i++) {
	if (i == Grid.INl-Grid.Nghost) {
	  norm_dir = -Grid.nfaceN(Grid.ICl-Grid.Nghost,Grid.JCu);
	} else {
	  norm_dir = -HALF*(Grid.nfaceN(i,Grid.JCu) + 
			    Grid.nfaceN(i-1,Grid.JCu));
	}
	X_norm = ((Grid.Node[i][Grid.JNu].X - 
		   Grid.Node[i][Grid.JNu-ng].X) * norm_dir) * norm_dir;
	X_tan = (Grid.Node[i][Grid.JNu].X - 
		 Grid.Node[i][Grid.JNu-ng].X) - X_norm;
	Grid.Node[i][Grid.JNu+ng].X = Grid.Node[i][Grid.JNu].X +
		                      X_norm - X_tan;
      }
    }
  }

  // NORTH-EAST corner:
  if (Grid.BCtypeN[Grid.INu+1] == BC_NONE &&
      Grid.BCtypeE[Grid.JNu+1] == BC_NONE) {
    // Do nothing.

  } else if ((Grid.BCtypeN[Grid.INu+1] != BC_REFLECTION &&
	      Grid.BCtypeN[Grid.INu+1] != BC_PERIODIC &&
	      Grid.BCtypeN[Grid.INu+1] != BC_NO_SLIP &&
	      Grid.BCtypeN[Grid.INu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      Grid.BCtypeN[Grid.INu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	      Grid.BCtypeN[Grid.INu+1] != BC_MOVING_WALL &&
	      Grid.BCtypeN[Grid.INu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	      Grid.BCtypeN[Grid.INu+1] != BC_MOVING_WALL_HEATFLUX &&
	      Grid.BCtypeN[Grid.INu+1] != BC_BURNING_SURFACE &&
	      Grid.BCtypeN[Grid.INu+1] != BC_MASS_INJECTION) &&
	     Grid.BCtypeE[Grid.JNu+1] == BC_NONE) {
    // Extrapolate cells north.
    for (int ng = 1; ng <= Grid.Nghost; ng++) {
      for (int i = Grid.INu; i <= Grid.INu+Grid.Nghost; i++) {
	Grid.Node[i][Grid.JNu+ng].X = Grid.Node[i][Grid.JNu].X +
                                      (Grid.Node[i][Grid.JNu].X - 
				       Grid.Node[i][Grid.JNu-ng].X);
      }
    }

  } else if ((Grid.BCtypeN[Grid.INu+1] == BC_NONE ||
	      (Grid.BCtypeN[Grid.INu+1] != BC_REFLECTION &&
	       Grid.BCtypeN[Grid.INu+1] != BC_PERIODIC &&
	       Grid.BCtypeN[Grid.INu+1] != BC_NO_SLIP &&
	       Grid.BCtypeN[Grid.INu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       Grid.BCtypeN[Grid.INu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	       Grid.BCtypeN[Grid.INu+1] != BC_MOVING_WALL &&
	       Grid.BCtypeN[Grid.INu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	       Grid.BCtypeN[Grid.INu+1] != BC_MOVING_WALL_HEATFLUX &&
	       Grid.BCtypeN[Grid.INu+1] != BC_BURNING_SURFACE &&
	       Grid.BCtypeN[Grid.INu+1] != BC_MASS_INJECTION)) &&
	     (Grid.BCtypeE[Grid.JNu+1] != BC_REFLECTION &&
	      Grid.BCtypeE[Grid.JNu+1] != BC_PERIODIC &&  
	      Grid.BCtypeE[Grid.JNu+1] != BC_NO_SLIP &&
	      Grid.BCtypeE[Grid.JNu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	      Grid.BCtypeE[Grid.JNu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	      Grid.BCtypeE[Grid.JNu+1] != BC_MOVING_WALL &&
	      Grid.BCtypeE[Grid.JNu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	      Grid.BCtypeE[Grid.JNu+1] != BC_MOVING_WALL_HEATFLUX &&
	      Grid.BCtypeE[Grid.JNu+1] != BC_BURNING_SURFACE &&
	      Grid.BCtypeE[Grid.JNu+1] != BC_MASS_INJECTION)) {
    // Extrapolate cells east.
    for (int ng = 1; ng <= Grid.Nghost; ng++) {
      for (int j = Grid.JNu; j <= Grid.JNu+Grid.Nghost; j++) {
	Grid.Node[Grid.INu+ng][j].X = Grid.Node[Grid.INu][j].X +
	                              (Grid.Node[Grid.INu][j].X - 
				       Grid.Node[Grid.INu-ng][j].X);
      }
    }

  } else if ((Grid.BCtypeN[Grid.INu+1] == BC_NONE ||
	      (Grid.BCtypeN[Grid.INu+1] != BC_REFLECTION &&
	       Grid.BCtypeN[Grid.INu+1] != BC_PERIODIC &&
	       Grid.BCtypeN[Grid.INu+1] != BC_NO_SLIP &&
	       Grid.BCtypeN[Grid.INu+1] != BC_WALL_VISCOUS_ISOTHERMAL &&
	       Grid.BCtypeN[Grid.INu+1] != BC_WALL_VISCOUS_HEATFLUX &&
	       Grid.BCtypeN[Grid.INu+1] != BC_MOVING_WALL &&
	       Grid.BCtypeN[Grid.INu+1] != BC_MOVING_WALL_ISOTHERMAL &&
	       Grid.BCtypeN[Grid.INu+1] != BC_MOVING_WALL_HEATFLUX &&
	       Grid.BCtypeN[Grid.INu+1] != BC_BURNING_SURFACE &&
	       Grid.BCtypeN[Grid.INu+1] != BC_MASS_INJECTION)) &&
	     (Grid.BCtypeE[Grid.JNu+1] == BC_REFLECTION ||
	      Grid.BCtypeE[Grid.JNu+1] == BC_PERIODIC ||
	      Grid.BCtypeE[Grid.JNu+1] == BC_NO_SLIP ||
	      Grid.BCtypeE[Grid.JNu+1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	      Grid.BCtypeE[Grid.JNu+1] == BC_WALL_VISCOUS_HEATFLUX ||
	      Grid.BCtypeE[Grid.JNu+1] == BC_MOVING_WALL ||
	      Grid.BCtypeE[Grid.JNu+1] == BC_MOVING_WALL_ISOTHERMAL ||
	      Grid.BCtypeE[Grid.JNu+1] == BC_MOVING_WALL_HEATFLUX ||
	      Grid.BCtypeE[Grid.JNu+1] == BC_BURNING_SURFACE ||
	      Grid.BCtypeE[Grid.JNu+1] == BC_MASS_INJECTION)) {
    // Reflect cells east.
    for (int ng = 1; ng <= Grid.Nghost; ng++) {
      for (int j = Grid.JNu; j <= Grid.JNu+Grid.Nghost; j++) {
	if (j == Grid.JNu+Grid.Nghost) {
	  norm_dir = Grid.nfaceE(Grid.ICu,Grid.JCu+Grid.Nghost);
	} else {
	  norm_dir = HALF*(Grid.nfaceE(Grid.ICu,j) + 
			   Grid.nfaceE(Grid.ICu,j-1));
	}
	X_norm = ((Grid.Node[Grid.INu][j].X - 
		   Grid.Node[Grid.INu-ng][j].X)*norm_dir)*norm_dir;
	X_tan = (Grid.Node[Grid.INu][j].X - 
		 Grid.Node[Grid.INu-ng][j].X) - X_norm;
	Grid.Node[Grid.INu+ng][j].X = Grid.Node[Grid.INu][j].X +
  	                              X_norm - X_tan;
      }
    }

  } else if (Grid.BCtypeN[Grid.INu+1] == BC_REFLECTION ||
	     Grid.BCtypeN[Grid.INu+1] == BC_PERIODIC ||
	     Grid.BCtypeN[Grid.INu+1] == BC_NO_SLIP ||
	     Grid.BCtypeN[Grid.INu+1] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     Grid.BCtypeN[Grid.INu+1] == BC_WALL_VISCOUS_HEATFLUX ||
	     Grid.BCtypeN[Grid.INu+1] == BC_MOVING_WALL ||
	     Grid.BCtypeN[Grid.INu+1] == BC_MOVING_WALL_ISOTHERMAL ||
	     Grid.BCtypeN[Grid.INu+1] == BC_MOVING_WALL_HEATFLUX ||
	     Grid.BCtypeN[Grid.INu+1] == BC_BURNING_SURFACE ||
	     Grid.BCtypeN[Grid.INu+1] == BC_MASS_INJECTION) {
    // Reflect cells north.
    for (int ng = 1; ng <= Grid.Nghost; ng++) {
      for (int i = Grid.INu; i <= Grid.INu+Grid.Nghost; i++) {
	if (i == Grid.INu+Grid.Nghost) {
	  norm_dir = -Grid.nfaceN(Grid.ICu+Grid.Nghost,Grid.JCu);
	} else {
	  norm_dir = -HALF*(Grid.nfaceN(i,Grid.JCu) + 
			    Grid.nfaceN(i-1,Grid.JCu));
	}
	X_norm = ((Grid.Node[i][Grid.JNu].X - 
		   Grid.Node[i][Grid.JNu-ng].X) * norm_dir) * norm_dir;
	X_tan = (Grid.Node[i][Grid.JNu].X - 
		 Grid.Node[i][Grid.JNu-ng].X) - X_norm;
	Grid.Node[i][Grid.JNu+ng].X = Grid.Node[i][Grid.JNu].X +
		                      X_norm - X_tan;
      }
    }
  }

  // Update cell info.
  Update_Cells(Grid);

}

/********************************************************
 * Routine: Update_Cells                                *
 *                                                      *
 * Updates the cell information for the quadrilateral   *
 * mesh block.                                          *
 *                                                      *
 ********************************************************/
void Update_Cells(Grid2D_Quad_Block &Grid) {

    int i, j;

    for ( j = Grid.JCl-Grid.Nghost ; j <= Grid.JCu+Grid.Nghost ; ++j) {
        for ( i = Grid.ICl-Grid.Nghost ; i <= Grid.ICu+Grid.Nghost ; ++i) {
	  Grid.Cell[i][j].I = i;
	  Grid.Cell[i][j].J = j;
	  Grid.Cell[i][j].Xc = Grid.centroid(i, j);
	  Grid.Cell[i][j].A = Grid.area(i, j);
        } /* endfor */
    } /* endfor */

//   Vector2D X1, X2, X3, X4, Xc1, Xc2, X;
//   double A1, A2;
//   Polygon P;
//   cout << setprecision(14);
//     for ( j = Grid.JCl-Grid.Nghost ; j <= Grid.JCu+Grid.Nghost ; ++j) {
//         for ( i = Grid.ICl-Grid.Nghost ; i <= Grid.ICu+Grid.Nghost ; ++i) {
//   // Cell nodes in counter-clockwise order.
//   X1 = Grid.Node[i  ][j  ].X;
//   X2 = Grid.Node[i+1][j  ].X;
//   X3 = Grid.Node[i+1][j+1].X;
//   X4 = Grid.Node[i  ][j+1].X;
//   P.convert(X1,X2,X3,X4);
//   // Determine the centroid and area of the sub-triangles.
//   Xc1 = (X1+X2+X3)/3.0;
//   Xc2 = (X1+X3+X4)/3.0;
// //   A1 = HALF*((X1^X2) + (X2^X3) + (X3^X1));
// //   A2 = HALF*((X1^X3) + (X3^X4) + (X4^X1));
//   A1 = HALF*((X2-X1)^(X3-X1));
//   A2 = HALF*((X3-X4)^(X3-X2));
//   X = (A1*Xc1 + A2*Xc2)/(A1+A2);
// //   cout << endl << X1 << X2 << X3 << X4 << X << (A1*Xc1 + A2*Xc2) << " " << A1 << " " << A2 << 0.25*(X1+X2+X3+X4);
// //   if (!P.point_in_polygon(X)) {
// //     //cout << endl << X1 << X2 << X3 << X4 << Xc1 << Xc2 << X << 0.25*(X1+X2+X3+X4);
// //     cout << " FUCK";
// //   }
//   if (A1 < ZERO || A2 < ZERO) {
//     cout << endl << A1 << " " << A2 << X1 << X2 << X3 << X4 << X << 0.25*(X1+X2+X3+X4);
//   }
//   P.deallocate();
//         }
//     }

}

/********************************************************
 * Routine: Check_Quad_Block                            *
 *                                                      *
 * Check the validity of the quadrilateral mesh block.  *
 * Returns a non-zero result if mesh is not valid.      *
 *                                                      *
 ********************************************************/
int Check_Quad_Block(Grid2D_Quad_Block &Grid) {

    int i, j;

    for ( j = Grid.JCl ; j <= Grid.JCu ; ++j) {
        for ( i = Grid.ICl ; i <= Grid.ICu ; ++i) {
	  if (Grid.Cell[i][j].A <= ZERO) {
	    cout << endl << i << " " << j;
	    cout << endl << Grid.ICl << " " << Grid.ICu;
	    cout << endl << Grid.JCl << " " << Grid.JCu;
	    cout << endl << Grid.Cell[i][j].A;
	    cout << endl << Grid.Cell[i][j].Xc;
	    cout << endl << Grid.Node[i][j].X;
	    cout << endl << Grid.Node[i+1][j].X;
	    cout << endl << Grid.Node[i][j+1].X;
	    cout << endl << Grid.Node[i+1][j+1].X;
	    cout.flush();
   	    return(1);
          } /* endif */
        } /* endfor */
    } /* endfor */

    return(0);

}

/********************************************************
 * Routine: Write_Quad_Block_Definition                 *
 *                                                      *
 * Writes definition file information for 2D            *
 * quadrilateral grid block to the specified output     *
 * stream for retrieval and re-use purposes.            *
 *                                                      *
 ********************************************************/
void Write_Quad_Block_Definition(Grid2D_Quad_Block &Grid,
                                 ostream &Out_File) {

    Out_File << setprecision(14);
    if (Grid.NNi == 0 || Grid.NNj == 0) {
       Out_File << Grid.NCi << " " 
	        << Grid.NCj << " "
		<< Grid.Nghost << "\n";
    } else {
       Out_File << Grid.NCi << " " 
	        << Grid.NCj << " "
		<< Grid.Nghost << "\n" 
                << Grid.BndNorthSpline.np << "\n"
                << Grid.BndNorthSpline.type << "\n"
                << Grid.BndNorthSpline 
                << Grid.BndSouthSpline.np << "\n"
                << Grid.BndSouthSpline.type << "\n"
                << Grid.BndSouthSpline
                << Grid.BndEastSpline.np << "\n"
                << Grid.BndEastSpline.type << "\n"
                << Grid.BndEastSpline 
                << Grid.BndWestSpline.np << "\n"
                << Grid.BndWestSpline.type << "\n"
                << Grid.BndWestSpline 
                << GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH << "\n";
       Out_File.setf(ios::scientific);
       Out_File << Grid.StretchI << " " 
                << Grid.BetaI << " " 
                << Grid.TauI << "\n"
                << Grid.StretchJ << " " 
                << Grid.BetaJ << " " 
                << Grid.TauJ << "\n";
       Out_File.unsetf(ios::scientific);
       Out_File << Grid.OrthogonalN << " "
                << Grid.OrthogonalS << " "
                << Grid.OrthogonalE << " "
                << Grid.OrthogonalW << "\n";
    } /* endif */
    Out_File << setprecision(6);

}

/********************************************************
 * Routine: Read_Quad_Block_Definition                  *
 *                                                      *
 * Reads definition file for 2D quadrilateral grid      *
 * block from the specified output stream.              *
 *                                                      *
 ********************************************************/
void Read_Quad_Block_Definition(Grid2D_Quad_Block &Grid,
                                istream &In_File) {

    int i, j, ng, k, kx, ky, kx_max, ky_max, npts, spline_type;
    int node_init_procedure;
    double S_i, S_j, 
           s_north, s_south, s_east, s_west,
           smax_north, smax_south, smax_east, smax_west, 
           s_i, s_j, smax_i, smax_j,
           w_north, w_south, w_east, w_west, w_total;
    double smax_S1_NS, smax_S2_NS, smax_S1_EW, smax_S2_EW,
           S1_i, S2_i, S1_j, S2_j, dS_i, dS_j;
    Vector2D x_north, x_south, x_east, x_west;
    Spline2D S1_NS, S2_NS, S1_EW, S2_EW;

    In_File.setf(ios::skipws);
    In_File >> i >> j >> ng;
    In_File.unsetf(ios::skipws);

    if (i != 0 && j != 0 && ng != 0) {
       /* Allocate (re-allocate) memory for the cells and nodes 
          of the quadrilateral mesh block as required. */

       if (Grid.NCi != i || Grid.NCj != j || Grid.Nghost != ng) { 
          Grid.deallocate();
          Grid.allocate(i, j, ng);
       } /* endif */

       /* For each of the north, south, east, and west boundaries
          of this mesh block, read in the number of spline points, 
          allocate memory for the boundary splines, read in the 
          spline points, and finally calculate the spline 
          pathlengths. */

       In_File.setf(ios::skipws);
       In_File >> npts;
       In_File.unsetf(ios::skipws);
       if (Grid.BndNorthSpline.np != 0) Grid.BndNorthSpline.deallocate();
       Grid.BndNorthSpline.allocate(npts);
       In_File.setf(ios::skipws);
       In_File >> spline_type;
       In_File.unsetf(ios::skipws);
       Grid.BndNorthSpline.settype(spline_type);
       In_File >> Grid.BndNorthSpline;
       Grid.BndNorthSpline.pathlength();
       Grid.SminN = Grid.BndNorthSpline.sp[0];
       Grid.SmaxN = Grid.BndNorthSpline.sp[Grid.BndNorthSpline.np-1];

       In_File.setf(ios::skipws);
       In_File >> npts;
       In_File.unsetf(ios::skipws);
       if (Grid.BndSouthSpline.np != 0) Grid.BndSouthSpline.deallocate();
       Grid.BndSouthSpline.allocate(npts);
       In_File.setf(ios::skipws);
       In_File >> spline_type;
       In_File.unsetf(ios::skipws);
       Grid.BndSouthSpline.settype(spline_type);
       In_File >> Grid.BndSouthSpline;
       Grid.BndSouthSpline.pathlength();
       Grid.SminS = Grid.BndSouthSpline.sp[0];
       Grid.SmaxS = Grid.BndSouthSpline.sp[Grid.BndSouthSpline.np-1];

       In_File.setf(ios::skipws);
       In_File >> npts;
       In_File.unsetf(ios::skipws);
       if (Grid.BndEastSpline.np != 0) Grid.BndEastSpline.deallocate();
       Grid.BndEastSpline.allocate(npts);
       In_File.setf(ios::skipws);
       In_File >> spline_type;
       In_File.unsetf(ios::skipws);
       Grid.BndEastSpline.settype(spline_type);
       In_File >> Grid.BndEastSpline;
       Grid.BndEastSpline.pathlength();
       Grid.SminE = Grid.BndEastSpline.sp[0];
       Grid.SmaxE = Grid.BndEastSpline.sp[Grid.BndEastSpline.np-1];

       In_File.setf(ios::skipws);
       In_File >> npts;
       In_File.unsetf(ios::skipws);
       if (Grid.BndWestSpline.np != 0) Grid.BndWestSpline.deallocate();
       Grid.BndWestSpline.allocate(npts);
       In_File.setf(ios::skipws);
       In_File >> spline_type;
       In_File.unsetf(ios::skipws);
       Grid.BndWestSpline.settype(spline_type);
       In_File >> Grid.BndWestSpline;
       Grid.BndWestSpline.pathlength();
       Grid.SminW = Grid.BndWestSpline.sp[0];
       Grid.SmaxW = Grid.BndWestSpline.sp[Grid.BndWestSpline.np-1];

       /* Read the node initialization procedure for this 
          quadrilateral grid block. */

       In_File.setf(ios::skipws);
       In_File >> node_init_procedure;
       In_File.unsetf(ios::skipws);
      
       /* Input the type of node distribution stretching 
          functions to be used for each of the coordinate
          directions. Also read in the related stretching
          parameters. */

       In_File.setf(ios::skipws);
       In_File >> Grid.StretchI >> Grid.BetaI >> Grid.TauI;
       In_File >> Grid.StretchJ >> Grid.BetaJ >> Grid.TauJ;
       In_File.unsetf(ios::skipws);

       /* Input the grid orthogonality specifiers for
          the north, south, east, and west boundaries. */

       In_File.setf(ios::skipws);
       In_File >> Grid.OrthogonalN >> Grid.OrthogonalS
               >> Grid.OrthogonalE >> Grid.OrthogonalW;
       In_File.unsetf(ios::skipws);

       /* Compute the interior nodes for the quadrilateral mesh block. */

       smax_north = Grid.BndNorthSpline.sp[Grid.BndNorthSpline.np-1]-
                    Grid.BndNorthSpline.sp[0];
       smax_south = Grid.BndSouthSpline.sp[Grid.BndSouthSpline.np-1]-
                    Grid.BndSouthSpline.sp[0];

       smax_east = Grid.BndEastSpline.sp[Grid.BndEastSpline.np-1]-
                   Grid.BndEastSpline.sp[0];
       smax_west = Grid.BndWestSpline.sp[Grid.BndWestSpline.np-1]-
                   Grid.BndWestSpline.sp[0];

       dS_i = ONE/(TEN*double(Grid.INu-Grid.INl));
       dS_j = ONE/(TEN*double(Grid.JNu-Grid.JNl));

       for ( j = Grid.JNl ; j <= Grid.JNu ; ++j) {
           for ( i = Grid.INl ; i <= Grid.INu ; ++i) {
               S_j = double(j-Grid.JNl)/double(Grid.JNu-Grid.JNl);
               S_j = StretchingFcn(S_j, Grid.BetaJ, Grid.TauJ, Grid.StretchJ);
               s_east = S_j * smax_east + Grid.BndEastSpline.sp[0];
               x_east = Spline(s_east, Grid.BndEastSpline);
               s_west = S_j * smax_west + Grid.BndWestSpline.sp[0];
               x_west = Spline(s_west, Grid.BndWestSpline);

               S_i = double(i-Grid.INl)/double(Grid.INu-Grid.INl);
               S_i = StretchingFcn(S_i, Grid.BetaI, Grid.TauI, Grid.StretchI);
               s_north = S_i * smax_north + Grid.BndNorthSpline.sp[0];
               x_north = Spline(s_north, Grid.BndNorthSpline);
               s_south = S_i * smax_south + Grid.BndSouthSpline.sp[0];
               x_south = Spline(s_south, Grid.BndSouthSpline);

               if (j == Grid.JNl) {
	         Grid.Node[i][j].X = x_south;
               } else if (j == Grid.JNu) {
	         Grid.Node[i][j].X = x_north;
               } else if (i == Grid.INl) {
	         Grid.Node[i][j].X = x_west;
               } else if (i == Grid.INu) {
	         Grid.Node[i][j].X = x_east;
               } else {
                 s_i = (ONE-S_j)*s_south + S_j*s_north;
                 s_j = (ONE-S_i)*s_west + S_i*s_east;

                 smax_i = (ONE-S_j)*smax_south + S_j*smax_north;
                 smax_j = (ONE-S_i)*smax_west + S_i*smax_east;

                 switch(node_init_procedure) {
                   //============================================================
                   case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_NORTH_SOUTH :
                     w_south = (ONE-s_j/smax_j);
                     w_north = (s_j/smax_j);
                     w_total = w_south + w_north;
                     Grid.Node[i][j].X = (w_south*x_south + 
                                          w_north*x_north)/w_total;
                     break;
                   //============================================================
                   case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_EAST_WEST :
                     w_west =  (ONE-s_i/smax_i);
                     w_east =  (s_i/smax_i);
                     w_total = w_east + w_west;
                     Grid.Node[i][j].X = (w_west*x_west + 
                                          w_east*x_east)/w_total;
                     break;
                   //============================================================
                   case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_XY :
                     kx_max = 5;
                     dS_i = ONE/double(kx_max-1);
                     kx = int(floor(S_i/dS_i));
                     S1_i = max(ZERO, double(kx)*dS_i);
                     S2_i = min(ONE, S1_i + dS_i);

                     ky_max = 7;
                     dS_j = ONE/double(ky_max-1);
                     ky = int(floor(S_j/dS_j));
                     S1_j = max(ZERO, double(ky)*dS_j);
                     S2_j = min(ONE, S1_j + dS_j);

	      	     if (S1_i - TOLER > ZERO) {
                        S1_NS.allocate(2);
                        S1_NS.settype(SPLINE2D_LINEAR);

                        s_south = S1_i * smax_south + Grid.BndSouthSpline.sp[0];
                        S1_NS.Xp[0] = Spline(s_south, Grid.BndSouthSpline);
                        S1_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                        s_north = S1_i * smax_north + Grid.BndNorthSpline.sp[0],
                        S1_NS.Xp[1] = Spline(s_north, Grid.BndNorthSpline);
                        S1_NS.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

                        S1_NS.pathlength();
                        smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
                                     S1_NS.sp[0];
                     } else {
                        Copy_Spline(S1_NS, Grid.BndWestSpline);
                        smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
                                     S1_NS.sp[0];
                     } /* endif */

                     if (S2_i + TOLER < ONE) {
                        S2_NS.allocate(2);
                        S2_NS.settype(SPLINE2D_LINEAR);

                        s_south = S2_i * smax_south + Grid.BndSouthSpline.sp[0];
                        S2_NS.Xp[0] = Spline(s_south, Grid.BndSouthSpline);
                        S2_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                        s_north = S2_i * smax_north + Grid.BndNorthSpline.sp[0],
                        S2_NS.Xp[1] = Spline(s_north, Grid.BndNorthSpline);
                        S2_NS.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

                        S2_NS.pathlength();
                        smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
                                     S2_NS.sp[0];
                     } else {
                        Copy_Spline(S2_NS, Grid.BndEastSpline);
                        smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
                                     S2_NS.sp[0];
                     } /* endif */

		     if (S1_j - TOLER > ZERO) {
                        S1_EW.allocate(kx_max);
                        S1_EW.settype(SPLINE2D_LINEAR);

                        s_west = S1_j * smax_west + Grid.BndWestSpline.sp[0];
                        S1_EW.Xp[0] = Spline(s_west, Grid.BndWestSpline);
                        S1_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                        s_east = S1_j * smax_east + Grid.BndEastSpline.sp[0];
                        S1_EW.Xp[kx_max-1] = Spline(s_east, Grid.BndEastSpline);
                        S1_EW.tp[kx_max-1] = SPLINE2D_POINT_SHARP_CORNER;

                        for ( k = 1 ; k <= kx_max-2 ; ++k) {                     
                           s_north = double(k)*smax_north/double(kx_max-1) + Grid.BndNorthSpline.sp[0];
                           x_north = Spline(s_north, Grid.BndNorthSpline);

                           s_south = double(k)*smax_south/double(kx_max-1) + Grid.BndSouthSpline.sp[0];
                           x_south = Spline(s_south, Grid.BndSouthSpline);

                           w_south =  ONE-S1_j;
                           w_north =  S1_j;
                           w_total = w_north + w_south;
                           S1_EW.Xp[k] = (w_south*x_south + w_north*x_north)/w_total;
                           S1_EW.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
                        } /* endfor */

                        S1_EW.pathlength();
                        smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
                                     S1_EW.sp[0];
                     } else {
                        Copy_Spline(S1_EW, Grid.BndSouthSpline);
                        smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
                                     S1_EW.sp[0];
                     } /* endif */

                     if (S2_j + TOLER < ONE) {
                        S2_EW.allocate(kx_max);
                        S2_EW.settype(SPLINE2D_LINEAR);

                        s_west = S2_j * smax_west + Grid.BndWestSpline.sp[0];
                        S2_EW.Xp[0] = Spline(s_west, Grid.BndWestSpline);
                        S2_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                        s_east = S2_j * smax_east + Grid.BndEastSpline.sp[0];
                        S2_EW.Xp[kx_max-1] = Spline(s_east, Grid.BndEastSpline);
                        S2_EW.tp[kx_max-1] = SPLINE2D_POINT_SHARP_CORNER;

                        for ( k = 1 ; k <= kx_max-2 ; ++k) {                     
                           s_north = double(k)*smax_north/double(kx_max-1) + Grid.BndNorthSpline.sp[0];
                           x_north = Spline(s_north, Grid.BndNorthSpline);

                           s_south = double(k)*smax_south/double(kx_max-1) + Grid.BndSouthSpline.sp[0];
                           x_south = Spline(s_south, Grid.BndSouthSpline);

                           w_south =  ONE-S2_j;
                           w_north =  S2_j;
                           w_total = w_north + w_south;
                           S2_EW.Xp[k] = (w_south*x_south + w_north*x_north)/w_total;
                           S2_EW.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
                        } /* endfor */

                        S2_EW.pathlength();
                        smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
                                     S2_EW.sp[0];
                     } else {
                        Copy_Spline(S2_EW, Grid.BndNorthSpline);
                        smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
                                     S2_EW.sp[0];
                     } /* endif */

                     s_north = S_i * smax_S2_EW + S2_EW.sp[0];
                     x_north = Spline(s_north, S2_EW);
                     s_south = S_i * smax_S1_EW + S1_EW.sp[0];
                     x_south = Spline(s_south, S1_EW);

                     s_east = S_j * smax_S2_NS + S2_NS.sp[0];
                     x_east = Spline(s_east, S2_NS);
                     s_west = S_j * smax_S1_NS + S1_NS.sp[0];
                     x_west = Spline(s_west, S1_NS);
      
                     if (ky == 0 || ky == ky_max - 2) {
                        w_south =  ONE-(S_j-S1_j)/(S2_j-S1_j);
                        w_north =  (S_j-S1_j)/(S2_j-S1_j);
                        w_total = w_south + w_north;
                        Grid.Node[i][j].X = (w_south*x_south + 
                                             w_north*x_north)/w_total;
                     } else {
                        w_west =  ONE-(S_i-S1_i)/(S2_i-S1_i);
                        w_east =  (S_i-S1_i)/(S2_i-S1_i);
                        w_total = w_east + w_west;
                        Grid.Node[i][j].X = (w_west*x_west + 
                                             w_east*x_east)/w_total;
                     } /* endif */

                     S1_NS.deallocate(); S2_NS.deallocate();
                     S2_EW.deallocate(); S2_EW.deallocate();
                     break;
                   //============================================================
                   case GRID2D_QUAD_BLOCK_INIT_PROCEDURE_TRANS_FINITE_YX :
                     kx_max = 7;
                     dS_i = ONE/double(kx_max-1);
                     kx = int(floor(S_i/dS_i));
                     S1_i = max(ZERO, double(kx)*dS_i);
                     S2_i = min(ONE, S1_i + dS_i);

                     ky_max = 5;
                     dS_j = ONE/double(ky_max-1);
                     ky = int(floor(S_j/dS_j));
                     S1_j = max(ZERO, double(ky)*dS_j);
                     S2_j = min(ONE, S1_j + dS_j);

		     if (S1_i - TOLER > ZERO) {
                        S1_NS.allocate(ky_max);
                        S1_NS.settype(SPLINE2D_LINEAR);

                        s_south = S1_i * smax_south + Grid.BndSouthSpline.sp[0];
                        S1_NS.Xp[0] = Spline(s_south, Grid.BndSouthSpline);
                        S1_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                        s_north = S1_i * smax_north + Grid.BndNorthSpline.sp[0],
                        S1_NS.Xp[ky_max-1] = Spline(s_north, Grid.BndNorthSpline);
                        S1_NS.tp[ky_max-1] = SPLINE2D_POINT_SHARP_CORNER;

                        for ( k = 1 ; k <= ky_max-2 ; ++k) {                     
                           s_east = double(k)*smax_east/double(ky_max-1) + Grid.BndEastSpline.sp[0];
                           x_east = Spline(s_east, Grid.BndEastSpline);

                           s_west = double(k)*smax_west/double(ky_max-1) + Grid.BndWestSpline.sp[0];
                           x_west = Spline(s_west, Grid.BndWestSpline);

                           w_west =  ONE-S1_i;
                           w_east =  S1_i;
                           w_total = w_east + w_west;
                           S1_NS.Xp[k] = (w_west*x_west + w_east*x_east)/w_total;
                           S1_NS.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
                        } /* endfor */

                        S1_NS.pathlength();
                        smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
                                     S1_NS.sp[0];
                     } else {
                        Copy_Spline(S1_NS, Grid.BndWestSpline);
                        smax_S1_NS = S1_NS.sp[S1_NS.np-1]-
                                     S1_NS.sp[0];
                     } /* endif */

                     if (S2_i + TOLER < ONE) {
                        S2_NS.allocate(ky_max);
                        S2_NS.settype(SPLINE2D_LINEAR);

                        s_south = S2_i * smax_south + Grid.BndSouthSpline.sp[0];
                        S2_NS.Xp[0] = Spline(s_south, Grid.BndSouthSpline);
                        S2_NS.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                        s_north = S2_i * smax_north + Grid.BndNorthSpline.sp[0],
                        S2_NS.Xp[ky_max-1] = Spline(s_north, Grid.BndNorthSpline);
                        S2_NS.tp[ky_max-1] = SPLINE2D_POINT_SHARP_CORNER;

                        for ( k = 1 ; k <= ky_max-2 ; ++k) {
                           s_east = double(k)*smax_east/double(ky_max-1) + Grid.BndEastSpline.sp[0];
                           x_east = Spline(s_east, Grid.BndEastSpline);

                           s_west = double(k)*smax_west/double(ky_max-1) + Grid.BndWestSpline.sp[0];
                           x_west = Spline(s_west, Grid.BndWestSpline);

                           w_west =  ONE-S2_i;
                           w_east =  S2_i;
                           w_total = w_east + w_west;
                           S2_NS.Xp[k] = (w_west*x_west + w_east*x_east)/w_total;
                           S2_NS.tp[k] = SPLINE2D_POINT_SHARP_CORNER;
                        } /* endfor */

                        S2_NS.pathlength();
                        smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
                                     S2_NS.sp[0];
                     } else {
                        Copy_Spline(S2_NS, Grid.BndEastSpline);
                        smax_S2_NS = S2_NS.sp[S2_NS.np-1]-
                                     S2_NS.sp[0];
                     } /* endif */

		     if (S1_j - TOLER > ZERO) {
                        S1_EW.allocate(2);
                        S1_EW.settype(SPLINE2D_LINEAR);

                        s_west = S1_j * smax_west + Grid.BndWestSpline.sp[0];
                        S1_EW.Xp[0] = Spline(s_west, Grid.BndWestSpline);
                        S1_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                        s_east = S1_j * smax_east + Grid.BndEastSpline.sp[0];
                        S1_EW.Xp[1] = Spline(s_east, Grid.BndEastSpline);
                        S1_EW.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

                        S1_EW.pathlength();
                        smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
                                     S1_EW.sp[0];
                     } else {
                        Copy_Spline(S1_EW, Grid.BndSouthSpline);
                        smax_S1_EW = S1_EW.sp[S1_EW.np-1]-
                                     S1_EW.sp[0];
                     } /* endif */

                     if (S2_j + TOLER < ONE) {
                        S2_EW.allocate(2);
                        S2_EW.settype(SPLINE2D_LINEAR);
   
                        s_west = S2_j * smax_west + Grid.BndWestSpline.sp[0];
                        S2_EW.Xp[0] = Spline(s_west, Grid.BndWestSpline);
                        S2_EW.tp[0] = SPLINE2D_POINT_SHARP_CORNER;

                        s_east = S2_j * smax_east + Grid.BndEastSpline.sp[0];
                        S2_EW.Xp[1] = Spline(s_east, Grid.BndEastSpline);
                        S2_EW.tp[1] = SPLINE2D_POINT_SHARP_CORNER;

                        S2_EW.pathlength();
                        smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
                                     S2_EW.sp[0];
                     } else {
                        Copy_Spline(S2_EW, Grid.BndNorthSpline);
                        smax_S2_EW = S2_EW.sp[S2_EW.np-1]-
                                     S2_EW.sp[0];
                     } /* endif */

                     s_north = S_i * smax_S2_EW + S2_EW.sp[0];
                     x_north = Spline(s_north, S2_EW);
                     s_south = S_i * smax_S1_EW + S1_EW.sp[0];
                     x_south = Spline(s_south, S1_EW);

                     s_east = S_j * smax_S2_NS + S2_NS.sp[0];
                     x_east = Spline(s_east, S2_NS);
                     s_west = S_j * smax_S1_NS + S1_NS.sp[0];
                     x_west = Spline(s_west, S1_NS);
      
                     if (kx == 0 || kx == kx_max - 2) {
                        w_west =  ONE-(S_i-S1_i)/(S2_i-S1_i);
                        w_east =  (S_i-S1_i)/(S2_i-S1_i);
                        w_total = w_east + w_west;
                        Grid.Node[i][j].X = (w_west*x_west + 
                                             w_east*x_east)/w_total;
                     } else {
                        w_south =  ONE-(S_j-S1_j)/(S2_j-S1_j);
                        w_north =  (S_j-S1_j)/(S2_j-S1_j);
                        w_total = w_south + w_north;
                        Grid.Node[i][j].X = (w_south*x_south + 
                                             w_north*x_north)/w_total;
                     } /* endif */

                     S1_NS.deallocate(); S2_NS.deallocate();
                     S2_EW.deallocate(); S2_EW.deallocate();
                     break;
                   //============================================================
                   default:
                     w_south = (ONE-s_j/smax_j);
                     w_north = (s_j/smax_j);
                     w_total = w_south + w_north;
                     Grid.Node[i][j].X = (w_south*x_south + 
                                          w_north*x_north)/w_total;
                     break;
                 } /* endswitch */

               } /* endif */

           } /* endfor */
       } /* endfor */

       /* Set the boundary condition types at the quadrilateral 
          grid block boundaries. */

       Set_BCs(Grid);

       /* Compute the exterior nodes for the quadrilateral mesh block. */

       Update_Exterior_Nodes(Grid);

       /* Compute the cells for the quadrilateral mesh block. */

       Update_Cells(Grid);
    } /* endif */

}

/********************************************************
 * Routine: Write_Quad_Block                            *
 *                                                      *
 * Writes 2D quadrilateral grid block to the            *
 * specified output stream for retrieval and re-use     *
 * purposes.                                            *
 *                                                      *
 ********************************************************/
void Write_Quad_Block(Grid2D_Quad_Block &Grid,
	              ostream &Out_File) {

    Out_File << setprecision(14) << Grid << setprecision(6);

}

/********************************************************
 * Routine: Read_Quad_Block                             *
 *                                                      *
 * Reads 2D quadrilateral grid block from the           *
 * specified input stream.                              *
 *                                                      *
 ********************************************************/
void Read_Quad_Block(Grid2D_Quad_Block &Grid,
	             istream &In_File) {

    In_File >> Grid;

}

/********************************************************
 * Routine: Copy_Quad_Block                             *
 *                                                      *
 * Copies the node and cell locations of quadrilateral  *
 * mesh block Grid2 to Grid1.                           *
 *                                                      *
 ********************************************************/
void Copy_Quad_Block(Grid2D_Quad_Block &Grid1,
                     Grid2D_Quad_Block &Grid2) {

    int i, j, ni, nj;
 
    /* Allocate (re-allocate) memory for the cells and nodes 
       of the quadrilateral mesh block Grid1 as necessary. */

    if (Grid1.NNi != Grid2.NNi || Grid1.NNj != Grid2.NNj ||
        Grid1.NCi != Grid2.NCi || Grid1.NCj != Grid2.NCj ) { 
       if (Grid1.Node != NULL && Grid1.Cell != NULL) { 
          Grid1.deallocate();
       } else if (Grid1.Node != NULL) {
          Grid1.deallocateNodes();
       } else if (Grid1.Cell != NULL) {
          Grid1.deallocateCells();
       } /* endif */
       if (Grid2.Node != NULL) Grid1.allocate(Grid2.NCi-2*Grid2.Nghost, 
                                              Grid2.NCj-2*Grid2.Nghost,
					      Grid2.Nghost);
    } /* endif */

    /* Copy the node locations of grid block Grid2 to Grid1. */

    if (Grid2.Node != NULL) {
       for (j  = Grid2.JNl-Grid2.Nghost; j <= Grid2.JNu+Grid2.Nghost ; ++j ) {
           for ( i = Grid2.INl-Grid2.Nghost ; i <= Grid2.INu+Grid2.Nghost ; ++i ) {
 	       Grid1.Node[i][j].X = Grid2.Node[i][j].X;
           } /* endfor */
       } /* endfor */
    } /* endif */

    /* Copy the cell values of grid block Grid2 to Grid1. */

    if (Grid2.Node != NULL) {
        for ( j = Grid2.JCl-Grid2.Nghost; j <= Grid2.JCu+Grid2.Nghost ; ++j) {
            for ( i = Grid2.ICl-Grid2.Nghost ; i <= Grid2.ICu+Grid2.Nghost ; ++i) {
	        Grid1.Cell[i][j].I  = Grid2.Cell[i][j].I;
  	        Grid1.Cell[i][j].J  = Grid2.Cell[i][j].J;
	        Grid1.Cell[i][j].Xc = Grid2.Cell[i][j].Xc;
		Grid1.Cell[i][j].A  = Grid2.Cell[i][j].A;
            } /* endfor */
        } /* endfor */
    } /* endif */

    /* Copy the boundary condition type info of grid block 
       Grid2 to Grid1. */

    if (Grid2.Node != NULL) {
        for ( i = Grid2.ICl-Grid2.Nghost ; i <= Grid2.ICu+Grid2.Nghost ; ++i) {
	    Grid1.BCtypeN[i] = Grid2.BCtypeN[i];
	    Grid1.BCtypeS[i] = Grid2.BCtypeS[i];
        } /* endfor */
        for ( j = Grid2.JCl-Grid2.Nghost ; j <= Grid2.JCu+Grid2.Nghost ; ++j) {
	    Grid1.BCtypeE[j] = Grid2.BCtypeE[j];
	    Grid1.BCtypeW[j] = Grid2.BCtypeW[j];
        } /* endfor */
    } /* endif */

    /* Copy boundary spline info of grid block 
       Grid2 to Grid1. */

    if (Grid2.BndNorthSpline.np != 0) {
       Copy_Spline(Grid1.BndNorthSpline, Grid2.BndNorthSpline);
    } else if (Grid1.BndNorthSpline.np != 0) {
       Grid1.BndNorthSpline.deallocate();
    } /* endif */

    if (Grid2.BndSouthSpline.np != 0) {
       Copy_Spline(Grid1.BndSouthSpline, Grid2.BndSouthSpline);
    } else if (Grid1.BndSouthSpline.np != 0) {
       Grid1.BndSouthSpline.deallocate();
    } /* endif */

    if (Grid2.BndEastSpline.np != 0) {
       Copy_Spline(Grid1.BndEastSpline, Grid2.BndEastSpline);
    } else if (Grid1.BndEastSpline.np != 0) {
       Grid1.BndEastSpline.deallocate();
    } /* endif */

    if (Grid2.BndWestSpline.np != 0) {
       Copy_Spline(Grid1.BndWestSpline, Grid2.BndWestSpline);
    } else if (Grid1.BndWestSpline.np != 0) {
       Grid1.BndWestSpline.deallocate();
    } /* endif */

    /* Copy boundary spline pathlength info of grid block 
       Grid2 to Grid1. */

    Grid1.SminN = Grid2.SminN;
    Grid1.SmaxN = Grid2.SmaxN;
    Grid1.SminS = Grid2.SminS;
    Grid1.SmaxS = Grid2.SmaxS;
    Grid1.SminE = Grid2.SminE;
    Grid1.SmaxE = Grid2.SmaxE;
    Grid1.SminW = Grid2.SminW;
    Grid1.SmaxW = Grid2.SmaxW;

    /* Copy node stretching info of grid block 
       Grid2 to Grid1. */

    Grid1.StretchI = Grid2.StretchI;
    Grid1.BetaI = Grid2.BetaI;
    Grid1.TauI = Grid2.TauI;
    Grid1.StretchJ = Grid2.StretchJ;
    Grid1.BetaJ = Grid2.BetaJ;
    Grid1.TauJ = Grid2.TauJ;
    Grid1.OrthogonalN = Grid2.OrthogonalN;
    Grid1.OrthogonalS = Grid2.OrthogonalS;
    Grid1.OrthogonalE = Grid2.OrthogonalE;
    Grid1.OrthogonalW = Grid2.OrthogonalW;

}

/********************************************************
 * Routine: Translate_Quad_Block                        *
 *                                                      *
 * Translates or shifts the positions of the nodes of a *
 * quadrilateral grid block.                            *
 *                                                      *
 ********************************************************/
void Translate_Quad_Block(Grid2D_Quad_Block &Grid,
                          const Vector2D &V) {

    int i, j;

    for ( j = Grid.JNl-Grid.Nghost ; j <= Grid.JNu+Grid.Nghost; ++j ) {
       for ( i = Grid.INl-Grid.Nghost ; i <= Grid.INu+Grid.Nghost; ++i ) {
           Grid.Node[i][j].X += V;
       } /* endfor */
    } /* endfor */

    for ( j = Grid.JCl-Grid.Nghost ; j <= Grid.JCu+Grid.Nghost ; ++j) {
        for ( i = Grid.ICl-Grid.Nghost ; i <= Grid.ICu+Grid.Nghost ; ++i) {
  	    Grid.Cell[i][j].Xc = Grid.centroid(i, j);
            Grid.Cell[i][j].A = Grid.area(i, j);
        } /* endfor */
    } /* endfor */

    if (Grid.BndNorthSpline.np != 0 ) 
       Translate_Spline(Grid.BndNorthSpline, V);
    if (Grid.BndSouthSpline.np != 0 ) 
       Translate_Spline(Grid.BndSouthSpline, V);
    if (Grid.BndEastSpline.np != 0 ) 
       Translate_Spline(Grid.BndEastSpline, V);
    if (Grid.BndWestSpline.np != 0 )
       Translate_Spline(Grid.BndWestSpline, V);
 
}

/********************************************************
 * Routine: Scale_Quad_Block                            *
 *                                                      *
 * Scales the quadrilateral grid block.                 *
 *                                                      *
 ********************************************************/
void Scale_Quad_Block(Grid2D_Quad_Block &Grid,
                      const double &Scaling_Factor) {

    int i, j;

    for ( j = Grid.JNl-Grid.Nghost ; j <= Grid.JNu+Grid.Nghost; ++j ) {
       for ( i = Grid.INl-Grid.Nghost ; i <= Grid.INu+Grid.Nghost; ++i ) {
           Grid.Node[i][j].X = Grid.Node[i][j].X*Scaling_Factor;
       } /* endfor */
    } /* endfor */

    for ( j = Grid.JCl-Grid.Nghost ; j <= Grid.JCu+Grid.Nghost ; ++j) {
        for ( i = Grid.ICl-Grid.Nghost ; i <= Grid.ICu+Grid.Nghost ; ++i) {
  	    Grid.Cell[i][j].Xc = Grid.centroid(i, j);
            Grid.Cell[i][j].A = Grid.area(i, j);
        } /* endfor */
    } /* endfor */

    if (Grid.BndNorthSpline.np != 0 ) 
       Scale_Spline(Grid.BndNorthSpline, Scaling_Factor);
    if (Grid.BndSouthSpline.np != 0 ) 
       Scale_Spline(Grid.BndSouthSpline, Scaling_Factor);
    if (Grid.BndEastSpline.np != 0 ) 
       Scale_Spline(Grid.BndEastSpline, Scaling_Factor);
    if (Grid.BndWestSpline.np != 0 )
       Scale_Spline(Grid.BndWestSpline, Scaling_Factor);

    Grid.SminN = Grid.SminN*Scaling_Factor;
    Grid.SmaxN = Grid.SmaxN*Scaling_Factor;
    Grid.SminS = Grid.SminS*Scaling_Factor;
    Grid.SmaxS = Grid.SmaxS*Scaling_Factor;
    Grid.SminE = Grid.SminE*Scaling_Factor;
    Grid.SmaxE = Grid.SmaxE*Scaling_Factor;
    Grid.SminW = Grid.SminW*Scaling_Factor;
    Grid.SmaxW = Grid.SmaxW*Scaling_Factor;
 
}

/********************************************************
 * Routine: Rotate_Quad_Block                           *
 *                                                      *
 * Rotates the quadrilateral grid block.                *
 *                                                      *
 ********************************************************/
void Rotate_Quad_Block(Grid2D_Quad_Block &Grid,
                       const double &Angle) {

    int i, j;
    double cos_angle, sin_angle;
    Vector2D X;

    cos_angle = cos(-Angle); 
    sin_angle = sin(-Angle);

    for ( j = Grid.JNl-Grid.Nghost ; j <= Grid.JNu+Grid.Nghost; ++j ) {
       for ( i = Grid.INl-Grid.Nghost ; i <= Grid.INu+Grid.Nghost; ++i ) {
           X.x = Grid.Node[i][j].X.x*cos_angle +
                 Grid.Node[i][j].X.y*sin_angle;
           X.y = - Grid.Node[i][j].X.x*sin_angle +
                   Grid.Node[i][j].X.y*cos_angle;
           Grid.Node[i][j].X = X;
       } /* endfor */
    } /* endfor */

    for ( j = Grid.JCl-Grid.Nghost ; j <= Grid.JCu+Grid.Nghost ; ++j) {
        for ( i = Grid.ICl-Grid.Nghost ; i <= Grid.ICu+Grid.Nghost ; ++i) {
  	    Grid.Cell[i][j].Xc = Grid.centroid(i, j);
            Grid.Cell[i][j].A = Grid.area(i, j);
        } /* endfor */
    } /* endfor */

    if (Grid.BndNorthSpline.np != 0 ) 
       Rotate_Spline(Grid.BndNorthSpline, Angle);
    if (Grid.BndSouthSpline.np != 0 ) 
       Rotate_Spline(Grid.BndSouthSpline, Angle);
    if (Grid.BndEastSpline.np != 0 ) 
       Rotate_Spline(Grid.BndEastSpline, Angle);
    if (Grid.BndWestSpline.np != 0 )
       Rotate_Spline(Grid.BndWestSpline, Angle);
 
}

/********************************************************
 * Routine: Reflect_Quad_Block                          *
 *                                                      *
 * Re-computes the locations of the nodes and cells of  *
 * the quadrilateral grid block based on a mirror       *
 * reflection about the y=0 axis.  The cells and nodes  *
 * are also re-ordered in the i-direction.              *
 *                                                      *
 ********************************************************/
void Reflect_Quad_Block(Grid2D_Quad_Block &Grid) {

    int i, j;
    Vector2D *X;
    Spline2D S;

    X = new Vector2D[Grid.NNi];

    for ( j = Grid.JNl-Grid.Nghost ; j <= Grid.JNu+Grid.Nghost; ++j ) {
       for ( i = Grid.INl-Grid.Nghost ; i <= Grid.INu+Grid.Nghost; ++i ) {
	   X[Grid.NNi-1-i].x = Grid.Node[i][j].X.x;
           X[Grid.NNi-1-i].y = - Grid.Node[i][j].X.y;
       } /* endfor */
       for ( i = Grid.INl-Grid.Nghost ; i <= Grid.INu+Grid.Nghost; ++i ) {
	   Grid.Node[i][j].X = X[i];
       } /* endfor */
    } /* endfor */

    delete []X;
    X = NULL;

    for ( j = Grid.JCl-Grid.Nghost ; j <= Grid.JCu+Grid.Nghost ; ++j) {
        for ( i = Grid.ICl-Grid.Nghost ; i <= Grid.ICu+Grid.Nghost ; ++i) {
  	    Grid.Cell[i][j].Xc = Grid.centroid(i, j);
            Grid.Cell[i][j].A = Grid.area(i, j);
        } /* endfor */
    } /* endfor */

    if (Grid.BndNorthSpline.np != 0 ) 
       Reflect_Spline(Grid.BndNorthSpline);
    if (Grid.BndSouthSpline.np != 0 ) 
       Reflect_Spline(Grid.BndSouthSpline);
    if (Grid.BndEastSpline.np != 0 ) 
       Reflect_Spline(Grid.BndEastSpline);
    if (Grid.BndWestSpline.np != 0 ) {
       Reflect_Spline(Grid.BndWestSpline);
       Copy_Spline(S,Grid.BndEastSpline);
       Copy_Spline(Grid.BndEastSpline,Grid.BndWestSpline);
       Copy_Spline(Grid.BndWestSpline,S);
       if (S.np != 0) S.deallocate();
    } /* endif */

    /* Reset the boundary condition types at the quadrilateral 
       grid block boundaries. */

    Set_BCs(Grid);
 
}

/********************************************************
 * Routine: Output_Tecplot                              *
 *                                                      *
 * Writes the nodes of the quadrilateral mesh to the    *
 * specified output stream in a format suitable for     *
 * plotting the grid with TECPLOT.                      *
 *                                                      *
 ********************************************************/
void Output_Tecplot(Grid2D_Quad_Block &Grid,
                    const int Block_Number,
                    const int Output_Title,
	            ostream &Out_File) {

    int i, j;

    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFFC_Name()
                << ": 2D Structured Curvilinear Grid Block (Node Locations)"
                << "\"" << "\n"
	        << "VARIABLES = \"x\" \\ \n"
                << "\"y\" \n"
                << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << Grid.INu - Grid.INl + 1 << " \\ \n"
                << "J = " << Grid.JNu - Grid.JNl + 1 << " \\ \n"
                << "F = POINT \n";

    } else {
       Out_File << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << Grid.INu - Grid.INl + 1 << " \\ \n"
                << "J = " << Grid.JNu - Grid.JNl + 1 << " \\ \n"
                << "F = POINT \n";
    } /* endif */

    for (j  = Grid.JNl ; j <= Grid.JNu ; ++j ) {
       for ( i = Grid.INl ; i <= Grid.INu ; ++i ) {
	   Out_File << " " << Grid.Node[i][j].X << "\n";
       } /* endfor */
    } /* endfor */
    Out_File << setprecision(6);

}

/********************************************************
 * Routine: Output_Nodes_Tecplot                        *
 *                                                      *
 * Writes the nodes of the quadrilateral mesh to the    *
 * specified output stream in a format suitable for     *
 * plotting the grid with TECPLOT.  Includes boundary   *
 * nodes.                                               *
 *                                                      *
 ********************************************************/
void Output_Nodes_Tecplot(Grid2D_Quad_Block &Grid,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

    int i, j;

    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFFC_Name()
                << ": 2D Structured Curvilinear Grid Block (Node Locations)"
                << "\"" << "\n"
	        << "VARIABLES = \"x\" \\ \n"
                << "\"y\" \n"
                << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << Grid.INu - Grid.INl + 1 + 2*Grid.Nghost << " \\ \n"
                << "J = " << Grid.JNu - Grid.JNl + 1 + 2*Grid.Nghost << " \\ \n"
                << "F = POINT \n";

    } else {
       Out_File << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << Grid.INu - Grid.INl + 1 + 2*Grid.Nghost << " \\ \n"
                << "J = " << Grid.JNu - Grid.JNl + 1 + 2*Grid.Nghost << " \\ \n"
                << "F = POINT \n";
    } /* endif */

    for (j  = Grid.JNl-Grid.Nghost ; j <= Grid.JNu+Grid.Nghost ; ++j ) {
       for ( i = Grid.INl-Grid.Nghost ; i <= Grid.INu+Grid.Nghost ; ++i ) {
	   Out_File << " " << Grid.Node[i][j].X << "\n";
       } /* endfor */
    } /* endfor */
    Out_File << setprecision(6);

}

/********************************************************
 * Routine: Output_Cells_Tecplot                        *
 *                                                      *
 * Writes the cells of the quadrilateral mesh to the    *
 * specified output stream in a format suitable for     *
 * plotting the grid with TECPLOT.                      *
 *                                                      *
 ********************************************************/
void Output_Cells_Tecplot(Grid2D_Quad_Block &Grid,
                          const int Block_Number,
                          const int Output_Title,
	                  ostream &Out_File) {

    int i, j;

    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "TITLE = \"" << CFFC_Name()
                << ": 2D Structured Curvilinear Grid Block (Cell Locations)"
                << "\"" << "\n"
	        << "VARIABLES = \"x\" \\ \n"
                << "\"y\" \n"
                << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << Grid.ICu - Grid.ICl + 1 << " \\ \n"
                << "J = " << Grid.JCu - Grid.JCl + 1 << " \\ \n"
                << "F = POINT \n";

    } else {
       Out_File << "ZONE T =  \"Block Number = " << Block_Number
                << "\" \\ \n"
                << "I = " << Grid.ICu - Grid.ICl + 1 << " \\ \n"
                << "J = " << Grid.JCu - Grid.JCl + 1 << " \\ \n"
                << "F = POINT \n";
    } /* endif */

    for (j  = Grid.JCl ; j <= Grid.JCu ; ++j ) {
       for ( i = Grid.ICl ; i <= Grid.ICu ; ++i ) {
	   Out_File << " " << Grid.Cell[i][j].Xc << "\n";
       } /* endfor */
    } /* endfor */
    Out_File << setprecision(6);

}

/********************************************************
 * Routine: Output_Gnuplot                              *
 *                                                      *
 * Writes the nodes of the quadrilateral mesh to the    *
 * specified output stream in a format suitable for     *
 * plotting the grid with GNUPLOT.                      *
 *                                                      *
 ********************************************************/
void Output_Gnuplot(Grid2D_Quad_Block &Grid,
                    const int Block_Number,
                    const int Output_Title,
	            ostream &Out_File) {

    int i, j;

    Out_File << setprecision(14);
    if (Output_Title) {
       Out_File << "# " << CFFC_Name()
                << ": 2D Structured Curvilinear Grid Block (Node Locations)"
                << "\n"
	        << "# x(m), y(m)\n";
    } /* endif */

    for (j  = Grid.JNl ; j <= Grid.JNu ; ++j ) {
       for ( i = Grid.INl ; i <= Grid.INu ; ++i ) {
	   Out_File << " " << Grid.Node[i][j].X << "\n";
       } /* endfor */
       Out_File << "\n";
    } /* endfor */
    for (i  = Grid.INl ; i <= Grid.INu ; ++i ) {
       for ( j = Grid.JNl ; j <= Grid.JNu ; ++j ) {
	   Out_File << " " << Grid.Node[i][j].X << "\n";
       } /* endfor */
       Out_File << "\n";
    } /* endfor */
    Out_File << setprecision(6);

}

/********************************************************
 * Routine: Double_Mesh_Resolution                      *
 *                                                      *
 * Returns a new quadrilateral mesh block with twice    *
 * the mesh resolution of the input grid block.         *
 *                                                      *
 ********************************************************/
void Double_Mesh_Resolution(Grid2D_Quad_Block &Grid_Double,
                            Grid2D_Quad_Block &Grid_Original) {

    int i, j, double_resolution_permitted;
    double sp_l, sp_r, sp_m, ds_ratio;
 
    /* Allocate memory for the cells and nodes of the 
       quadrilateral mesh block with twice the resolution. */

    if ( (Grid_Original.NCi-2*Grid_Original.Nghost-2*((Grid_Original.NCi-2*Grid_Original.Nghost)/2) != 0) || 
         (Grid_Original.NCj-2*Grid_Original.Nghost-2*((Grid_Original.NCj-2*Grid_Original.Nghost)/2) != 0) ||
         (Grid_Original.NCi-2*Grid_Original.Nghost < 2*Grid_Original.Nghost) ||
         (Grid_Original.NCj-2*Grid_Original.Nghost < 2*Grid_Original.Nghost) ||
         (Grid_Original.Node == NULL) ) { 
      double_resolution_permitted = 0;
    } else {
      double_resolution_permitted = 1;
      if (Grid_Double.Node != NULL && Grid_Double.Cell != NULL) { 
          Grid_Double.deallocate();
       } else if (Grid_Double.Node != NULL) {
          Grid_Double.deallocateNodes();
       } else if (Grid_Double.Cell != NULL) {
          Grid_Double.deallocateCells();
       } /* endif */
       Grid_Double.allocate(2*(Grid_Original.NCi-2*Grid_Original.Nghost), 
                            2*(Grid_Original.NCj-2*Grid_Original.Nghost),
			    Grid_Original.Nghost);
    } /* endif */

    /* Copy boundary spline info to quadrilateral mesh block 
       with twice the resolution. */

    if (double_resolution_permitted) {

       if (Grid_Original.BndNorthSpline.np != 0) {
          Copy_Spline(Grid_Double.BndNorthSpline, 
                      Grid_Original.BndNorthSpline);
       } /* endif */

       if (Grid_Original.BndSouthSpline.np != 0) {
          Copy_Spline(Grid_Double.BndSouthSpline, 
                      Grid_Original.BndSouthSpline);
       } /* endif */

       if (Grid_Original.BndEastSpline.np != 0) {
          Copy_Spline(Grid_Double.BndEastSpline, 
                      Grid_Original.BndEastSpline);
       } /* endif */
  
       if (Grid_Original.BndWestSpline.np != 0) {
          Copy_Spline(Grid_Double.BndWestSpline, 
                      Grid_Original.BndWestSpline);
       } /* endif */

    /* Copy boundary spline pathlength info to quadrilateral mesh block 
       with twice the resolution. */

       Grid_Double.SminN = Grid_Original.SminN;
       Grid_Double.SmaxN = Grid_Original.SmaxN;
       Grid_Double.SminS = Grid_Original.SminS;
       Grid_Double.SmaxS = Grid_Original.SmaxS;
       Grid_Double.SminE = Grid_Original.SminE;
       Grid_Double.SmaxE = Grid_Original.SmaxE;
       Grid_Double.SminW = Grid_Original.SminW;
       Grid_Double.SmaxW = Grid_Original.SmaxW;

    /* Copy node stretching info to quadrilateral mesh block 
       with twice the resolution. */

       Grid_Double.StretchI = Grid_Original.StretchI;
       Grid_Double.BetaI = Grid_Original.BetaI;
       Grid_Double.TauI = Grid_Original.TauI;
       Grid_Double.StretchJ = Grid_Original.StretchJ;
       Grid_Double.BetaJ = Grid_Original.BetaJ;
       Grid_Double.TauJ = Grid_Original.TauJ;
       Grid_Double.OrthogonalN = Grid_Original.OrthogonalN;
       Grid_Double.OrthogonalS = Grid_Original.OrthogonalS;
       Grid_Double.OrthogonalE = Grid_Original.OrthogonalE;
       Grid_Double.OrthogonalW = Grid_Original.OrthogonalW;

    /* Determine the node locations of quadrilateral mesh block 
       with twice the resolution. */

       for (j  = Grid_Original.JCl ; j <= Grid_Original.JCu ; ++j ) {
           for ( i = Grid_Original.ICl ; i <= Grid_Original.ICu ; ++i ) {
    	       Grid_Double.Node[2*(i-Grid_Original.INl)+Grid_Original.INl  ]
                               [2*(j-Grid_Original.JNl)+Grid_Original.JNl  ].X 
                  = Grid_Original.nodeSW(i, j).X;
    	       Grid_Double.Node[2*(i-Grid_Original.INl)+Grid_Original.INl+1]
                               [2*(j-Grid_Original.JNl)+Grid_Original.JNl  ].X 
                  = Grid_Original.xfaceS(i, j);
    	       Grid_Double.Node[2*(i-Grid_Original.INl)+Grid_Original.INl  ]
                               [2*(j-Grid_Original.JNl)+Grid_Original.JNl+1].X 
                  = Grid_Original.xfaceW(i, j);
    	       Grid_Double.Node[2*(i-Grid_Original.INl)+Grid_Original.INl+1]
                               [2*(j-Grid_Original.JNl)+Grid_Original.JNl+1].X 
                  = Grid_Original.Cell[i][j].Xc;
               if (j == Grid_Original.JCu) {
                  Grid_Double.Node[2*(i-Grid_Original.INl)+Grid_Original.INl  ]
                                  [2*(j-Grid_Original.JNl)+Grid_Original.JNl+2].X 
                     = Grid_Original.nodeNW(i, j).X;
                  Grid_Double.Node[2*(i-Grid_Original.INl)+Grid_Original.INl+1]
                                  [2*(j-Grid_Original.JNl)+Grid_Original.JNl+2].X 
                     = Grid_Original.xfaceN(i, j);
               } /* endif */
               if (i == Grid_Original.ICu) {
                  Grid_Double.Node[2*(i-Grid_Original.INl)+Grid_Original.INl+2]
                                  [2*(j-Grid_Original.JNl)+Grid_Original.JNl  ].X 
                     = Grid_Original.nodeSE(i, j).X;
                  Grid_Double.Node[2*(i-Grid_Original.INl)+Grid_Original.INl+2]
                                  [2*(j-Grid_Original.JNl)+Grid_Original.JNl+1].X 
                     = Grid_Original.xfaceE(i, j);
               } /* endif */
               if (i == Grid_Original.ICu && j == Grid_Original.JCu) {
                  Grid_Double.Node[2*(i-Grid_Original.INl)+Grid_Original.INl+2]
                                  [2*(j-Grid_Original.JNl)+Grid_Original.JNl+2].X 
                  = Grid_Original.nodeNE(i, j).X;
               } /* endif */
           } /* endfor */
       } /* endfor */

       if (Grid_Double.BndWestSpline.np != 0) {
          for (j  = Grid_Double.JNl+1 ; j < Grid_Double.JNu ; j += 2 ) {
 	      sp_l = getS(Grid_Double.Node[Grid_Double.INl][j-1].X, 
                          Grid_Double.BndWestSpline);
 	      sp_r = getS(Grid_Double.Node[Grid_Double.INl][j+1].X, 
                          Grid_Double.BndWestSpline);
              ds_ratio = abs(Grid_Double.Node[Grid_Double.INl+1][j].X-
                             Grid_Double.Node[Grid_Double.INl+1][j-1].X)/
                         abs(Grid_Double.Node[Grid_Double.INl+1][j+1].X-
                             Grid_Double.Node[Grid_Double.INl+1][j-1].X);
              sp_m = sp_l + ds_ratio*(sp_r-sp_l);
	      Grid_Double.Node[Grid_Double.INl][j].X = 
                 Spline(sp_m, Grid_Double.BndWestSpline);
          } /* endfor */
       } /* endif */

       if (Grid_Double.BndEastSpline.np != 0) {
          for (j  = Grid_Double.JNl+1 ; j < Grid_Double.JNu ; j += 2 ) {
   	      sp_l = getS(Grid_Double.Node[Grid_Double.INu][j-1].X, 
                          Grid_Double.BndEastSpline);
 	      sp_r = getS(Grid_Double.Node[Grid_Double.INu][j+1].X, 
                          Grid_Double.BndEastSpline);
              ds_ratio = abs(Grid_Double.Node[Grid_Double.INu-1][j].X-
                             Grid_Double.Node[Grid_Double.INu-1][j-1].X)/
                         abs(Grid_Double.Node[Grid_Double.INu-1][j+1].X-
                             Grid_Double.Node[Grid_Double.INu-1][j-1].X);
              sp_m = sp_l + ds_ratio*(sp_r-sp_l);
    	      Grid_Double.Node[Grid_Double.INu][j].X = 
	         Spline(sp_m, Grid_Double.BndEastSpline);
          } /* endfor */
       } /* endif */

       if (Grid_Double.BndSouthSpline.np != 0) {
          for ( i = Grid_Double.INl+1 ; i < Grid_Double.INu ; i += 2 ) {
 	      sp_l = getS(Grid_Double.Node[i-1][Grid_Double.JNl].X, 
                          Grid_Double.BndSouthSpline);
              sp_r = getS(Grid_Double.Node[i+1][Grid_Double.JNl].X, 
                          Grid_Double.BndSouthSpline);
              ds_ratio = abs(Grid_Double.Node[i][Grid_Double.JNl+1].X-
                             Grid_Double.Node[i-1][Grid_Double.JNl+1].X)/
                         abs(Grid_Double.Node[i+1][Grid_Double.JNl+1].X-
                             Grid_Double.Node[i-1][Grid_Double.JNl+1].X);
              sp_m = sp_l + ds_ratio*(sp_r-sp_l);
    	      Grid_Double.Node[i][Grid_Double.JNl].X =  
	         Spline(sp_m, Grid_Double.BndSouthSpline);
          } /* endfor */
       } /* endif */

       if (Grid_Double.BndNorthSpline.np != 0) {
          for ( i = Grid_Double.INl+1 ; i < Grid_Double.INu ; i += 2 ) {
              sp_l = getS(Grid_Double.Node[i-1][Grid_Double.JNu].X, 
                          Grid_Double.BndNorthSpline);
 	      sp_r = getS(Grid_Double.Node[i+1][Grid_Double.JNu].X, 
                          Grid_Double.BndNorthSpline);
              ds_ratio = abs(Grid_Double.Node[i][Grid_Double.JNu-1].X-
                             Grid_Double.Node[i-1][Grid_Double.JNu-1].X)/
                         abs(Grid_Double.Node[i+1][Grid_Double.JNu-1].X-
                             Grid_Double.Node[i-1][Grid_Double.JNu-1].X);
              sp_m = sp_l + ds_ratio*(sp_r-sp_l);
    	      Grid_Double.Node[i][Grid_Double.JNu].X =  
	         Spline(sp_m, Grid_Double.BndNorthSpline);
          } /* endfor */
       } /* endif */

    /* Set the boundary condition types for quadrilateral mesh block 
       with twice the resolution. */

       Set_BCs(Grid_Double);

    /* Compute the exterior nodes for quadrilateral mesh block 
       with twice the resolution. */

       Update_Exterior_Nodes(Grid_Double);

    /* Compute the cells for quadrilateral mesh block 
       with twice the resolution. */

       Update_Cells(Grid_Double);

    } /* endif */

}

/********************************************************
 * Routine: Half_Mesh_Resolution                        *
 *                                                      *
 * Returns a new quadrilateral mesh block with half the *
 * mesh resolution of the input grid block.             *
 *                                                      *
 ********************************************************/
void Half_Mesh_Resolution(Grid2D_Quad_Block &Grid_Half,
                          Grid2D_Quad_Block &Grid_Original) {

    int i, j, half_resolution_permitted;
 
    /* Allocate memory for the cells and nodes of the 
       quadrilateral mesh block with half the resolution. */

    if ( (Grid_Original.NCi-2*Grid_Original.Nghost-2*((Grid_Original.NCi-2*Grid_Original.Nghost)/2) != 0) || 
         (Grid_Original.NCj-2*Grid_Original.Nghost-2*((Grid_Original.NCj-2*Grid_Original.Nghost)/2) != 0) ||
         (Grid_Original.NCi-2*Grid_Original.Nghost < 2*Grid_Original.Nghost) ||
         (Grid_Original.NCj-2*Grid_Original.Nghost < 2*Grid_Original.Nghost) ||
         (Grid_Original.Node == NULL) ) {
       half_resolution_permitted = 0;
    } else {
       half_resolution_permitted = 1;
       if (Grid_Half.Node != NULL && Grid_Half.Cell != NULL) { 
          Grid_Half.deallocate();
       } else if (Grid_Half.Node != NULL) {
          Grid_Half.deallocateNodes();
       } else if (Grid_Half.Cell != NULL) {
          Grid_Half.deallocateCells();
       } /* endif */
       Grid_Half.allocate((Grid_Original.NCi-2*Grid_Original.Nghost)/2, 
                          (Grid_Original.NCj-2*Grid_Original.Nghost)/2,
			  Grid_Original.Nghost);
    } /* endif */

    /* Copy boundary spline info to quadrilateral mesh block 
       with half the resolution. */

    if (half_resolution_permitted) {

       if (Grid_Original.BndNorthSpline.np != 0) {
          Copy_Spline(Grid_Half.BndNorthSpline, 
                      Grid_Original.BndNorthSpline);
       } /* endif */

       if (Grid_Original.BndSouthSpline.np != 0) {
          Copy_Spline(Grid_Half.BndSouthSpline, 
                      Grid_Original.BndSouthSpline);
       } /* endif */

       if (Grid_Original.BndEastSpline.np != 0) {
          Copy_Spline(Grid_Half.BndEastSpline, 
                      Grid_Original.BndEastSpline);
       } /* endif */
  
       if (Grid_Original.BndWestSpline.np != 0) {
          Copy_Spline(Grid_Half.BndWestSpline, 
                      Grid_Original.BndWestSpline);
       } /* endif */

    /* Copy boundary spline pathlength info to quadrilateral mesh block 
       with half the resolution. */

       Grid_Half.SminN = Grid_Original.SminN;
       Grid_Half.SmaxN = Grid_Original.SmaxN;
       Grid_Half.SminS = Grid_Original.SminS;
       Grid_Half.SmaxS = Grid_Original.SmaxS;
       Grid_Half.SminE = Grid_Original.SminE;
       Grid_Half.SmaxE = Grid_Original.SmaxE;
       Grid_Half.SminW = Grid_Original.SminW;
       Grid_Half.SmaxW = Grid_Original.SmaxW;

    /* Copy node stretching info to quadrilateral mesh block 
       with half the resolution. */

       Grid_Half.StretchI = Grid_Original.StretchI;
       Grid_Half.BetaI = Grid_Original.BetaI;
       Grid_Half.TauI = Grid_Original.TauI;
       Grid_Half.StretchJ = Grid_Original.StretchJ;
       Grid_Half.BetaJ = Grid_Original.BetaJ;
       Grid_Half.TauJ = Grid_Original.TauJ;
       Grid_Half.OrthogonalN = Grid_Original.OrthogonalN;
       Grid_Half.OrthogonalS = Grid_Original.OrthogonalS;
       Grid_Half.OrthogonalE = Grid_Original.OrthogonalE;
       Grid_Half.OrthogonalW = Grid_Original.OrthogonalW;

    /* Determine the node locations of quadrilateral mesh block 
       with half the resolution. */

       for (j  = Grid_Half.JNl ; j <= Grid_Half.JNu ; ++j ) {
           for ( i = Grid_Half.INl ; i <= Grid_Half.INu ; ++i ) {
    	       Grid_Half.Node[i][j].X = 
                  Grid_Original.Node[2*(i-Grid_Half.INl)+Grid_Half.INl]
                                    [2*(j-Grid_Half.JNl)+Grid_Half.JNl].X;
           } /* endfor */
       } /* endfor */

    /* Set the boundary condition types for quadrilateral mesh block 
       with half the resolution. */

       Set_BCs(Grid_Half);

    /* Compute the exterior nodes for quadrilateral mesh block 
       with half the resolution. */

       Update_Exterior_Nodes(Grid_Half);

    /* Compute the cells for quadrilateral mesh block 
       with half the resolution. */

       Update_Cells(Grid_Half);

    } /* endif */

}

/********************************************************
 * Routine: Refine_Mesh                                 *
 *                                                      *
 * Returns a new quadrilateral mesh block containing    *
 * one of the specified sectors of the original grid    *
 * block with twice the mesh resolution.                *
 *                                                      *
 ********************************************************/
void Refine_Mesh(Grid2D_Quad_Block &Grid_Fine,
                 Grid2D_Quad_Block &Grid_Original,
                 const int Sector) {

    int i, j, i_min, i_max, j_min, j_max, mesh_refinement_permitted;
    double sp_l, sp_r, sp_m, ds_ratio, dl, dr;

    /* Allocate memory for the cells and nodes for the 
       refined quadrilateral mesh block. */

    if ( (Grid_Original.NCi-2*Grid_Original.Nghost-2*((Grid_Original.NCi-2*Grid_Original.Nghost)/2) != 0) || 
         (Grid_Original.NCj-2*Grid_Original.Nghost-2*((Grid_Original.NCj-2*Grid_Original.Nghost)/2) != 0) ||
         (Grid_Original.NCi-2*Grid_Original.Nghost < 2*Grid_Original.Nghost) ||
         (Grid_Original.NCj-2*Grid_Original.Nghost < 2*Grid_Original.Nghost) ||
         (Grid_Original.Node == NULL) ) {
       mesh_refinement_permitted = 0;
    } else {
       mesh_refinement_permitted = 1;
       if (Grid_Fine.Node != NULL && Grid_Fine.Cell != NULL) { 
          Grid_Fine.deallocate();
       } else if (Grid_Fine.Node != NULL) {
          Grid_Fine.deallocateNodes();
       } else if (Grid_Fine.Cell != NULL) {
          Grid_Fine.deallocateCells();
       } /* endif */
       Grid_Fine.allocate(Grid_Original.NCi-2*Grid_Original.Nghost, 
                          Grid_Original.NCj-2*Grid_Original.Nghost,
			  Grid_Original.Nghost);
    } /* endif */

    /* Copy boundary spline info for the refined
       quadrilateral mesh block. */

    if (mesh_refinement_permitted) {

       if ((Sector == GRID2D_QUAD_BLOCK_SECTOR_NW ||
            Sector == GRID2D_QUAD_BLOCK_SECTOR_NE) &&
           Grid_Original.BndNorthSpline.np != 0) {
          Copy_Spline(Grid_Fine.BndNorthSpline, 
                      Grid_Original.BndNorthSpline);
       } /* endif */

       if ((Sector == GRID2D_QUAD_BLOCK_SECTOR_SE ||
            Sector == GRID2D_QUAD_BLOCK_SECTOR_SW) &&
           Grid_Original.BndSouthSpline.np != 0) {
          Copy_Spline(Grid_Fine.BndSouthSpline, 
                      Grid_Original.BndSouthSpline);
       } /* endif */

       if ((Sector == GRID2D_QUAD_BLOCK_SECTOR_NE ||
            Sector == GRID2D_QUAD_BLOCK_SECTOR_SE) &&
           Grid_Original.BndEastSpline.np != 0) {
          Copy_Spline(Grid_Fine.BndEastSpline, 
                      Grid_Original.BndEastSpline);
       } /* endif */
  
       if ((Sector == GRID2D_QUAD_BLOCK_SECTOR_NW ||
            Sector == GRID2D_QUAD_BLOCK_SECTOR_SW) &&
           Grid_Original.BndWestSpline.np != 0) {
          Copy_Spline(Grid_Fine.BndWestSpline, 
                      Grid_Original.BndWestSpline);
       } /* endif */

    /* Assign boundary spline pathlength info for the refined
       quadrilateral mesh block. */

       switch(Sector) {
         case GRID2D_QUAD_BLOCK_SECTOR_NW :
           if (Grid_Original.BndNorthSpline.np != 0) {
              Grid_Fine.SminN = Grid_Original.SminN;
              Grid_Fine.SmaxN = getS(Grid_Original.Node[Grid_Original.INl+
                                                        (Grid_Original.INu-
                                                         Grid_Original.INl)/2]
                                                       [Grid_Original.JNu].X, 
                                     Grid_Original.BndNorthSpline);
           } else {
              Grid_Fine.SminN = ZERO;
              Grid_Fine.SmaxN = ZERO;
           } /* endif */
           Grid_Fine.SminS = ZERO;
           Grid_Fine.SmaxS = ZERO;
           Grid_Fine.SminE = ZERO;
           Grid_Fine.SmaxE = ZERO;
           if (Grid_Original.BndWestSpline.np != 0) {
              Grid_Fine.SminW = getS(Grid_Original.Node[Grid_Original.INl]
                                                       [Grid_Original.JNl+
                                                        (Grid_Original.JNu-
                                                         Grid_Original.JNl)/2].X, 
                                     Grid_Original.BndWestSpline);
              Grid_Fine.SmaxW = Grid_Original.SmaxW;
           } else {
              Grid_Fine.SminW = ZERO;
              Grid_Fine.SmaxW = ZERO;
           } /* endif */
           break;
         case GRID2D_QUAD_BLOCK_SECTOR_NE :
           if (Grid_Original.BndNorthSpline.np != 0) {
              Grid_Fine.SminN = getS(Grid_Original.Node[Grid_Original.INl+
                                                        (Grid_Original.INu-
                                                         Grid_Original.INl)/2]
                                                       [Grid_Original.JNu].X, 
                                     Grid_Original.BndNorthSpline);
              Grid_Fine.SmaxN = Grid_Original.SmaxN;
           } else {
              Grid_Fine.SminN = ZERO;
              Grid_Fine.SmaxN = ZERO;
           } /* endif */
           Grid_Fine.SminS = ZERO;
           Grid_Fine.SmaxS = ZERO;
           if (Grid_Original.BndEastSpline.np != 0) {
              Grid_Fine.SminE = getS(Grid_Original.Node[Grid_Original.INu]
                                                       [Grid_Original.JNl+
                                                        (Grid_Original.JNu-
                                                         Grid_Original.JNl)/2].X, 
                                     Grid_Original.BndEastSpline);
              Grid_Fine.SmaxE = Grid_Original.SmaxE;
           } else {
              Grid_Fine.SminE = ZERO;
              Grid_Fine.SmaxE = ZERO;
           } /* endif */
           Grid_Fine.SminW = ZERO;
           Grid_Fine.SmaxW = ZERO;
           break;
         case GRID2D_QUAD_BLOCK_SECTOR_SE :
           Grid_Fine.SminN = ZERO;
           Grid_Fine.SmaxN = ZERO;
           if (Grid_Original.BndSouthSpline.np != 0) {
              Grid_Fine.SminS = getS(Grid_Original.Node[Grid_Original.INl+
                                                        (Grid_Original.INu-
                                                         Grid_Original.INl)/2]
                                                       [Grid_Original.JNl].X, 
                                     Grid_Original.BndSouthSpline);
              Grid_Fine.SmaxS = Grid_Original.SmaxS;
           } else {
              Grid_Fine.SminS = ZERO;
              Grid_Fine.SmaxS = ZERO;
           } /* endif */
           if (Grid_Original.BndEastSpline.np != 0) {
              Grid_Fine.SminE = Grid_Original.SminE;
              Grid_Fine.SmaxE = getS(Grid_Original.Node[Grid_Original.INu]
                                                       [Grid_Original.JNl+
                                                        (Grid_Original.JNu-
                                                         Grid_Original.JNl)/2].X, 
                                     Grid_Original.BndEastSpline);
              Grid_Fine.SmaxE = Grid_Original.SmaxE;
           } else {
              Grid_Fine.SminE = ZERO;
              Grid_Fine.SmaxE = ZERO;
           } /* endif */
           Grid_Fine.SminW = ZERO;
           Grid_Fine.SmaxW = ZERO;
           break;
         case GRID2D_QUAD_BLOCK_SECTOR_SW :
           Grid_Fine.SminN = ZERO;
           Grid_Fine.SmaxN = ZERO;
           if (Grid_Original.BndSouthSpline.np != 0) {
              Grid_Fine.SminS = Grid_Original.SminS;
              Grid_Fine.SmaxS = getS(Grid_Original.Node[Grid_Original.INl+
                                                        (Grid_Original.INu-
                                                         Grid_Original.INl)/2]
                                                       [Grid_Original.JNl].X, 
                                     Grid_Original.BndSouthSpline);
           } else {
              Grid_Fine.SminS = ZERO;
              Grid_Fine.SmaxS = ZERO;
           } /* endif */
           Grid_Fine.SminE = ZERO;
           Grid_Fine.SmaxE = ZERO;
           if (Grid_Original.BndWestSpline.np != 0) {
              Grid_Fine.SminW = Grid_Original.SminW;
              Grid_Fine.SmaxW = getS(Grid_Original.Node[Grid_Original.INl]
                                                       [Grid_Original.JNl+
                                                        (Grid_Original.JNu-
                                                         Grid_Original.JNl)/2].X, 
                                     Grid_Original.BndWestSpline);
           } else {
              Grid_Fine.SminW = ZERO;
              Grid_Fine.SmaxW = ZERO;
           } /* endif */
           break;
       } /* endswitch */

    /* Copy node stretching info to refined quadrilateral 
       mesh block. */

       Grid_Fine.StretchI = Grid_Original.StretchI;
       Grid_Fine.BetaI = Grid_Original.BetaI;
       Grid_Fine.TauI = Grid_Original.TauI;
       Grid_Fine.StretchJ = Grid_Original.StretchJ;
       Grid_Fine.BetaJ = Grid_Original.BetaJ;
       Grid_Fine.TauJ = Grid_Original.TauJ;
       if (Grid_Original.BndNorthSpline.np != 0) {
          Grid_Fine.OrthogonalN = Grid_Original.OrthogonalN;
       } else {
          Grid_Fine.OrthogonalN = ORTHOGONAL;
       } /* endif */
       if (Grid_Original.BndSouthSpline.np != 0) {
          Grid_Fine.OrthogonalS = Grid_Original.OrthogonalS;
       } else {
          Grid_Fine.OrthogonalS = ORTHOGONAL;
       } /* endif */
       if (Grid_Original.BndEastSpline.np != 0) {
          Grid_Fine.OrthogonalE = Grid_Original.OrthogonalE;
       } else {
          Grid_Fine.OrthogonalE = ORTHOGONAL;
       } /* endif */
       if (Grid_Original.BndWestSpline.np != 0) {
          Grid_Fine.OrthogonalW = Grid_Original.OrthogonalW;
       } else {
          Grid_Fine.OrthogonalW = ORTHOGONAL;
       } /* endif */

       // Force orthogonality at all refined mesh boundaries.
       Grid_Fine.OrthogonalN = ORTHOGONAL;
       Grid_Fine.OrthogonalS = ORTHOGONAL;
       Grid_Fine.OrthogonalE = ORTHOGONAL;
       Grid_Fine.OrthogonalW = ORTHOGONAL;

    /* Determine the node locations for refined 
       quadrilateral mesh block. */

       switch(Sector) {
         case GRID2D_QUAD_BLOCK_SECTOR_NW :
           i_min = Grid_Original.ICl;
	   i_max = Grid_Original.ICl+(Grid_Original.ICu-Grid_Original.ICl-1)/2;
           j_min = Grid_Original.JCl+(Grid_Original.JCu-Grid_Original.JCl+1)/2; 
           j_max = Grid_Original.JCu;
           break;
         case GRID2D_QUAD_BLOCK_SECTOR_NE :
           i_min = Grid_Original.ICl+(Grid_Original.ICu-Grid_Original.ICl+1)/2;
	   i_max = Grid_Original.ICu;
           j_min = Grid_Original.JCl+(Grid_Original.JCu-Grid_Original.JCl+1)/2; 
           j_max = Grid_Original.JCu;
           break;
         case GRID2D_QUAD_BLOCK_SECTOR_SE :
           i_min = Grid_Original.ICl+(Grid_Original.ICu-Grid_Original.ICl+1)/2;
	   i_max = Grid_Original.ICu;
           j_min = Grid_Original.JCl; 
           j_max = Grid_Original.JCl+(Grid_Original.JCu-Grid_Original.JCl-1)/2;
           break;
         case GRID2D_QUAD_BLOCK_SECTOR_SW :
           i_min = Grid_Original.ICl;
	   i_max = Grid_Original.ICl+(Grid_Original.ICu-Grid_Original.ICl-1)/2;
           j_min = Grid_Original.JCl; 
           j_max = Grid_Original.JCl+(Grid_Original.JCu-Grid_Original.JCl-1)/2;
           break;
         default:
           i_min = Grid_Original.ICl;
	   i_max = Grid_Original.ICl+(Grid_Original.ICu-Grid_Original.ICl-1)/2;
           j_min = Grid_Original.JCl+(Grid_Original.JCu-Grid_Original.JCl+1)/2; 
           j_max = Grid_Original.JCu;
           break;
       } /* endswitch */

       for ( j  = j_min; j <= j_max ; ++j ) {
	   for ( i = i_min ; i <= i_max ; ++i ) {
    	       Grid_Fine.Node[2*(i-i_min)+Grid_Fine.INl  ]
                             [2*(j-j_min)+Grid_Fine.JNl  ].X 
                  = Grid_Original.nodeSW(i, j).X;
    	       Grid_Fine.Node[2*(i-i_min)+Grid_Fine.INl+1]
                             [2*(j-j_min)+Grid_Fine.JNl  ].X 
                  = Grid_Original.xfaceS(i, j);
    	       Grid_Fine.Node[2*(i-i_min)+Grid_Fine.INl  ]
                             [2*(j-j_min)+Grid_Fine.JNl+1].X 
                  = Grid_Original.xfaceW(i, j);
//     	       Grid_Fine.Node[2*(i-i_min)+Grid_Fine.INl+1]
//                              [2*(j-j_min)+Grid_Fine.JNl+1].X 
//                   = Grid_Original.Cell[i][j].Xc;
    	       Grid_Fine.Node[2*(i-i_min)+Grid_Fine.INl+1]
                             [2*(j-j_min)+Grid_Fine.JNl+1].X 
                  = (Grid_Original.nodeSW(i,j).X +
		     Grid_Original.nodeSE(i,j).X +
		     Grid_Original.nodeNW(i,j).X +
		     Grid_Original.nodeNE(i,j).X)/FOUR;
               if (j == j_max) {
                  Grid_Fine.Node[2*(i-i_min)+Grid_Fine.INl  ]
                                [2*(j-j_min)+Grid_Fine.JNl+2].X 
                  = Grid_Original.nodeNW(i, j).X;
                  Grid_Fine.Node[2*(i-i_min)+Grid_Fine.INl+1]
                                [2*(j-j_min)+Grid_Fine.JNl+2].X 
                  = Grid_Original.xfaceN(i, j);
               } /* endif */
               if (i == i_max) {
                  Grid_Fine.Node[2*(i-i_min)+Grid_Fine.INl+2]
                                [2*(j-j_min)+Grid_Fine.JNl  ].X 
                  = Grid_Original.nodeSE(i, j).X;
                  Grid_Fine.Node[2*(i-i_min)+Grid_Fine.INl+2]
                                [2*(j-j_min)+Grid_Fine.JNl+1].X 
                  = Grid_Original.xfaceE(i, j);
               } /* endif */
               if (i == i_max && j == j_max) {
                  Grid_Fine.Node[2*(i-i_min)+Grid_Fine.INl+2]
                                [2*(j-j_min)+Grid_Fine.JNl+2].X 
                  = Grid_Original.nodeNE(i, j).X;
               } /* endif */
           } /* endfor */
       } /* endfor */

       if (Grid_Fine.BndWestSpline.np != 0) {
          for (j  = Grid_Fine.JNl+1 ; j < Grid_Fine.JNu ; j += 2 ) {
 	      sp_l = getS(Grid_Fine.Node[Grid_Fine.INl][j-1].X, 
                          Grid_Fine.BndWestSpline);
 	      sp_r = getS(Grid_Fine.Node[Grid_Fine.INl][j+1].X, 
                          Grid_Fine.BndWestSpline);
              dl = abs(Grid_Fine.Node[Grid_Fine.INl][j  ].X - 
                       Grid_Fine.Node[Grid_Fine.INl][j-1].X);
              dr = abs(Grid_Fine.Node[Grid_Fine.INl][j+1].X - 
                       Grid_Fine.Node[Grid_Fine.INl][j  ].X);
              ds_ratio = dl/(dl+dr);
              sp_m = sp_l + ds_ratio*(sp_r-sp_l);
	      Grid_Fine.Node[Grid_Fine.INl][j].X = 
                 Spline(sp_m, Grid_Fine.BndWestSpline);
          } /* endfor */
       } /* endif */

       if (Grid_Fine.BndEastSpline.np != 0) {
          for (j  = Grid_Fine.JNl+1 ; j < Grid_Fine.JNu ; j += 2 ) {
              sp_l = getS(Grid_Fine.Node[Grid_Fine.INu][j-1].X, 
                          Grid_Fine.BndEastSpline);
 	      sp_r = getS(Grid_Fine.Node[Grid_Fine.INu][j+1].X, 
                          Grid_Fine.BndEastSpline);
              dl = abs(Grid_Fine.Node[Grid_Fine.INu][j  ].X - 
                       Grid_Fine.Node[Grid_Fine.INu][j-1].X);
              dr = abs(Grid_Fine.Node[Grid_Fine.INu][j+1].X - 
                       Grid_Fine.Node[Grid_Fine.INu][j  ].X);
              ds_ratio = dl/(dl+dr);
              sp_m = sp_l + ds_ratio*(sp_r-sp_l);
    	      Grid_Fine.Node[Grid_Fine.INu][j].X = 
	         Spline(sp_m, Grid_Fine.BndEastSpline);
          } /* endfor */
       } /* endif */

       if (Grid_Fine.BndSouthSpline.np != 0) {
          for ( i = Grid_Fine.INl+1 ; i < Grid_Fine.INu ; i += 2 ) {
 	      sp_l = getS(Grid_Fine.Node[i-1][Grid_Fine.JNl].X, 
                          Grid_Fine.BndSouthSpline);
 	      sp_r = getS(Grid_Fine.Node[i+1][Grid_Fine.JNl].X, 
                          Grid_Fine.BndSouthSpline);
              dl = abs(Grid_Fine.Node[i  ][Grid_Fine.JNl].X - 
                       Grid_Fine.Node[i-1][Grid_Fine.JNl].X);
              dr = abs(Grid_Fine.Node[i+1][Grid_Fine.JNl].X - 
                       Grid_Fine.Node[i  ][Grid_Fine.JNl].X);
              ds_ratio = dl/(dl+dr);
              sp_m = sp_l + ds_ratio*(sp_r-sp_l);
    	      Grid_Fine.Node[i][Grid_Fine.JNl].X =  
	         Spline(sp_m, Grid_Fine.BndSouthSpline);
          } /* endfor */
       } /* endif */

       if (Grid_Fine.BndNorthSpline.np != 0) {
          for ( i = Grid_Fine.INl+1 ; i < Grid_Fine.INu ; i += 2 ) {
 	      sp_l = getS(Grid_Fine.Node[i-1][Grid_Fine.JNu].X, 
                          Grid_Fine.BndNorthSpline);
 	      sp_r = getS(Grid_Fine.Node[i+1][Grid_Fine.JNu].X, 
                          Grid_Fine.BndNorthSpline);
              dl = abs(Grid_Fine.Node[i  ][Grid_Fine.JNu].X - 
                       Grid_Fine.Node[i-1][Grid_Fine.JNu].X);
              dr = abs(Grid_Fine.Node[i+1][Grid_Fine.JNu].X - 
                       Grid_Fine.Node[i  ][Grid_Fine.JNu].X);
              ds_ratio = dl/(dl+dr);
              sp_m = sp_l + ds_ratio*(sp_r-sp_l);
    	      Grid_Fine.Node[i][Grid_Fine.JNu].X =  
	         Spline(sp_m, Grid_Fine.BndNorthSpline);
          } /* endfor */
// 	 for ( i = Grid_Fine.ICl; i <= Grid_Fine.ICu; i++) {
// 	   if (Grid_Fine.area(i,Grid_Fine.JCu) <= ZERO) {
// 	     Grid_Fine.Node[i][Grid_Fine.JNu-1].X = Grid_Fine.Node[i][Grid_Fine.JNu-3].X + (2.0/3.0)*(Grid_Fine.Node[i][Grid_Fine.JNu].X-
// 												      Grid_Fine.Node[i][Grid_Fine.JNu-3].X);
// 	     Grid_Fine.Node[i][Grid_Fine.JNu-2].X = Grid_Fine.Node[i][Grid_Fine.JNu-3].X + (1.0/3.0)*(Grid_Fine.Node[i][Grid_Fine.JNu].X-
// 												      Grid_Fine.Node[i][Grid_Fine.JNu-3].X);
// 	     Grid_Fine.Node[i+1][Grid_Fine.JNu-1].X = Grid_Fine.Node[i+1][Grid_Fine.JNu-3].X + (2.0/3.0)*(Grid_Fine.Node[i+1][Grid_Fine.JNu].X-
// 												      Grid_Fine.Node[i+1][Grid_Fine.JNu-3].X);
// 	     Grid_Fine.Node[i+1][Grid_Fine.JNu-2].X = Grid_Fine.Node[i+1][Grid_Fine.JNu-3].X + (1.0/3.0)*(Grid_Fine.Node[i+1][Grid_Fine.JNu].X-
// 												      Grid_Fine.Node[i+1][Grid_Fine.JNu-3].X);
// 	   }
// 	 }
       } /* endif */

    /* Set the boundary condition types for refined 
       quadrilateral mesh block. */

       Set_BCs(Grid_Fine);

    /* Compute the exterior nodes for refined 
       quadrilateral mesh block. */

       Update_Exterior_Nodes(Grid_Fine);

    /* Compute the cells for refined 
       quadrilateral mesh block. */

       Update_Cells(Grid_Fine);

    } /* endif */

}

/********************************************************
 * Routine: Coarsen_Mesh                                *
 *                                                      *
 * Returns a new quadrilateral mesh block resulting     *
 * from the coarsening of four original grid blocks     *
 * with half the original mesh resolution.              *
 *                                                      *
 ********************************************************/
void Coarsen_Mesh(Grid2D_Quad_Block &Grid_Coarse,
                  Grid2D_Quad_Block &Grid_Original_SW,
                  Grid2D_Quad_Block &Grid_Original_SE,
                  Grid2D_Quad_Block &Grid_Original_NW,
                  Grid2D_Quad_Block &Grid_Original_NE) {

    int i, j, i_coarse, j_coarse, mesh_coarsening_permitted;
 
    /* Allocate memory for the cells and nodes for the 
       coarsened quadrilateral mesh block. */

    if ( (Grid_Original_SW.NCi-2*Grid_Original_SW.Nghost-
	  2*((Grid_Original_SW.NCi-2*Grid_Original_SW.Nghost)/2) != 0) || 
         (Grid_Original_SW.NCj-2*Grid_Original_SW.Nghost-
	  2*((Grid_Original_SW.NCj-2*Grid_Original_SW.Nghost)/2) != 0) ||
         (Grid_Original_SW.NCi-2*Grid_Original_SW.Nghost < 2*Grid_Original_SW.Nghost) ||
         (Grid_Original_SW.NCj-2*Grid_Original_SW.Nghost < 2*Grid_Original_SW.Nghost) ||
         (Grid_Original_SE.NCi != Grid_Original_SW.NCi) ||
         (Grid_Original_SE.NCj != Grid_Original_SW.NCj) ||
         (Grid_Original_NW.NCi != Grid_Original_SW.NCi) ||
         (Grid_Original_NW.NCj != Grid_Original_SW.NCj) ||
         (Grid_Original_NE.NCi != Grid_Original_SW.NCi) ||
         (Grid_Original_NE.NCj != Grid_Original_SW.NCj) ||
         (Grid_Original_SW.Node == NULL) ||
         (Grid_Original_SE.Node == NULL) ||
         (Grid_Original_NW.Node == NULL) ||
         (Grid_Original_NE.Node == NULL) ) {
       mesh_coarsening_permitted = 0;
    } else {
       mesh_coarsening_permitted = 1;
       if (Grid_Coarse.Node != NULL && Grid_Coarse.Cell != NULL) { 
          Grid_Coarse.deallocate();
       } else if (Grid_Coarse.Node != NULL) {
          Grid_Coarse.deallocateNodes();
       } else if (Grid_Coarse.Cell != NULL) {
          Grid_Coarse.deallocateCells();
       } /* endif */
       Grid_Coarse.allocate(Grid_Original_SW.NCi-2*Grid_Original_SW.Nghost, 
                            Grid_Original_SW.NCj-2*Grid_Original_SW.Nghost,
			    Grid_Original_SW.Nghost);
    } /* endif */

    /* Copy boundary spline info for the coarsened
       quadrilateral mesh block. */

    if (mesh_coarsening_permitted) {

       if (Grid_Original_NW.BndNorthSpline.np != 0 &&
           Grid_Original_NE.BndNorthSpline.np != 0) {
          Copy_Spline(Grid_Coarse.BndNorthSpline, 
                      Grid_Original_NW.BndNorthSpline);
       } /* endif */

       if (Grid_Original_SW.BndSouthSpline.np != 0 &&
           Grid_Original_SE.BndSouthSpline.np != 0) {
          Copy_Spline(Grid_Coarse.BndSouthSpline, 
                      Grid_Original_SW.BndSouthSpline);
       } /* endif */

       if (Grid_Original_SE.BndEastSpline.np != 0 &&
           Grid_Original_NE.BndEastSpline.np != 0) {
          Copy_Spline(Grid_Coarse.BndEastSpline, 
                      Grid_Original_SE.BndEastSpline);
       } /* endif */
  
       if (Grid_Original_SW.BndWestSpline.np != 0 &&
           Grid_Original_NW.BndWestSpline.np != 0) {
          Copy_Spline(Grid_Coarse.BndWestSpline, 
                      Grid_Original_SW.BndWestSpline);
       } /* endif */

    /* Assign boundary spline pathlength info for the coarsened
       quadrilateral mesh block. */

       if (Grid_Original_NW.BndNorthSpline.np != 0 &&
           Grid_Original_NE.BndNorthSpline.np != 0) {
          Grid_Coarse.SminN = Grid_Original_NW.SminN;
          Grid_Coarse.SmaxN = Grid_Original_NE.SmaxN;
       } else {
          Grid_Coarse.SminN = ZERO;
          Grid_Coarse.SmaxN = ZERO;
       } /* endif */

       if (Grid_Original_SW.BndSouthSpline.np != 0 &&
           Grid_Original_SE.BndSouthSpline.np != 0) {
          Grid_Coarse.SminS = Grid_Original_SW.SminS;
          Grid_Coarse.SmaxS = Grid_Original_SE.SmaxS;
       } else {
          Grid_Coarse.SminS = ZERO;
          Grid_Coarse.SmaxS = ZERO;
       } /* endif */

       if (Grid_Original_SE.BndEastSpline.np != 0 &&
           Grid_Original_NE.BndEastSpline.np != 0) {
          Grid_Coarse.SminE = Grid_Original_SE.SmaxE;
          Grid_Coarse.SmaxE = Grid_Original_NE.SmaxE;
       } else {
          Grid_Coarse.SminE = ZERO;
          Grid_Coarse.SmaxE = ZERO;
       } /* endif */

       if (Grid_Original_SW.BndWestSpline.np != 0 &&
           Grid_Original_NW.BndWestSpline.np != 0) {
          Grid_Coarse.SminW = Grid_Original_SW.SminW;
          Grid_Coarse.SmaxW = Grid_Original_NW.SmaxW;
       } else {
          Grid_Coarse.SminW = ZERO;
          Grid_Coarse.SmaxW = ZERO;
       } /* endif */

    /* Copy node stretching info to coarsened quadrilateral 
       mesh block. */

       Grid_Coarse.StretchI = Grid_Original_SW.StretchI;
       Grid_Coarse.BetaI = Grid_Original_SW.BetaI;
       Grid_Coarse.TauI = Grid_Original_SW.TauI;
       Grid_Coarse.StretchJ = Grid_Original_SW.StretchJ;
       Grid_Coarse.BetaJ = Grid_Original_SW.BetaJ;
       Grid_Coarse.TauJ = Grid_Original_SW.TauJ;
       if (Grid_Original_NW.BndNorthSpline.np != 0 &&
           Grid_Original_NE.BndNorthSpline.np != 0) {
          Grid_Coarse.OrthogonalN = Grid_Original_NW.OrthogonalN;
       } else {
          Grid_Coarse.OrthogonalN = ORTHOGONAL;
       } /* endif */
       if (Grid_Original_SW.BndSouthSpline.np != 0 &&
           Grid_Original_SE.BndSouthSpline.np != 0) {
          Grid_Coarse.OrthogonalS = Grid_Original_SW.OrthogonalS;
       } else {
          Grid_Coarse.OrthogonalS = ORTHOGONAL;
       } /* endif */
       if (Grid_Original_SE.BndEastSpline.np != 0 &&
           Grid_Original_NE.BndEastSpline.np != 0) {
          Grid_Coarse.OrthogonalE = Grid_Original_SE.OrthogonalE;
       } else {
          Grid_Coarse.OrthogonalE = ORTHOGONAL;
       } /* endif */
       if (Grid_Original_SW.BndWestSpline.np != 0 &&
           Grid_Original_NW.BndWestSpline.np != 0) {
          Grid_Coarse.OrthogonalW = Grid_Original_SW.OrthogonalW;
       } else {
          Grid_Coarse.OrthogonalW = ORTHOGONAL;
       } /* endif */

       // Force orthogonality at all coarsened mesh boundaries.
       Grid_Coarse.OrthogonalN = ORTHOGONAL;
       Grid_Coarse.OrthogonalS = ORTHOGONAL;
       Grid_Coarse.OrthogonalE = ORTHOGONAL;
       Grid_Coarse.OrthogonalW = ORTHOGONAL;

    /* Determine the node locations for the coarsened
       quadrilateral mesh block. */

       for ( j = Grid_Original_SW.JNl; j <= Grid_Original_SW.JNu ; j += 2 ) {
	   for ( i = Grid_Original_SW.INl ; i <= Grid_Original_SW.INu ; i += 2 ) {
 	      i_coarse = (i-Grid_Original_SW.INl)/2+
                         Grid_Coarse.INl;
	      j_coarse = (j-Grid_Original_SW.JNl)/2+
                         Grid_Coarse.JNl;
              Grid_Coarse.Node[i_coarse][j_coarse].X = Grid_Original_SW.Node[i][j].X;
           } /* endfor */
       } /* endfor */

       for ( j = Grid_Original_SE.JNl; j <= Grid_Original_SE.JNu ; j += 2 ) {
	   for ( i = Grid_Original_SE.INl ; i <= Grid_Original_SE.INu ; i += 2 ) {
 	      i_coarse = (i-Grid_Original_SE.INl)/2+
                         (Grid_Coarse.INu-Grid_Coarse.INl)/2+Grid_Coarse.INl;
	      j_coarse = (j-Grid_Original_SE.JNl)/2+
                         Grid_Coarse.JNl;
              Grid_Coarse.Node[i_coarse][j_coarse].X = Grid_Original_SE.Node[i][j].X;
           } /* endfor */
       } /* endfor */

       for ( j = Grid_Original_NW.JNl; j <= Grid_Original_NW.JNu ; j += 2 ) {
	   for ( i = Grid_Original_NW.INl ; i <= Grid_Original_NW.INu ; i += 2 ) {
 	      i_coarse = (i-Grid_Original_NW.INl)/2+
                         Grid_Coarse.INl;
	      j_coarse = (j-Grid_Original_NW.JNl)/2+
                         (Grid_Coarse.JNu-Grid_Coarse.JNl)/2+Grid_Coarse.JNl;
              Grid_Coarse.Node[i_coarse][j_coarse].X = Grid_Original_NW.Node[i][j].X;
           } /* endfor */
       } /* endfor */

       for ( j = Grid_Original_NE.JNl; j <= Grid_Original_NE.JNu ; j += 2 ) {
	   for ( i = Grid_Original_NE.INl ; i <= Grid_Original_NE.INu ; i += 2 ) {
 	      i_coarse = (i-Grid_Original_NE.INl)/2+
                         (Grid_Coarse.INu-Grid_Coarse.INl)/2+Grid_Coarse.INl;
	      j_coarse = (j-Grid_Original_NE.JNl)/2+
                         (Grid_Coarse.JNu-Grid_Coarse.JNl)/2+Grid_Coarse.JNl;
              Grid_Coarse.Node[i_coarse][j_coarse].X = Grid_Original_NE.Node[i][j].X;
           } /* endfor */
       } /* endfor */

    /* Set the boundary condition types for newly coarsened
       quadrilateral mesh block. */

       Set_BCs(Grid_Coarse);

    /* Compute the exterior nodes for newly coarsened 
       quadrilateral mesh block. */

       Update_Exterior_Nodes(Grid_Coarse);

    /* Compute the cells for newly coarsened
       quadrilateral mesh block. */

       Update_Cells(Grid_Coarse);

    } /* endif */

}

/********************************************************
 * Routine: Fix_Refined_Mesh_Boundaries                 *
 *                                                      *
 * Adjusts the locations of the boundary nodes of a     *
 * quadrilateral grid block so that the new node        *
 * locations will match with cell volumes of adjacent   *
 * quadrilateral mesh blocks that have lower levels of  *
 * mesh refinement (i.e., are coarser mesh blocks).     *
 *                                                      *
 ********************************************************/
void Fix_Refined_Mesh_Boundaries(Grid2D_Quad_Block &Grid,
                                 const int Fix_North_Boundary,
                                 const int Fix_South_Boundary,
                                 const int Fix_East_Boundary,
                                 const int Fix_West_Boundary) {

    int i, j;
    double ds_ratio, dl, dr;
 
    /* Adjust the node locations of the north boundary
       of the quadrilateral mesh block. */

    if (Fix_North_Boundary) {
       for ( i = Grid.INl+1 ; i <= Grid.INu-1 ; i+=2 ) {
          dl = abs(Grid.Node[i  ][Grid.JNu].X - 
                   Grid.Node[i-1][Grid.JNu].X);
          dr = abs(Grid.Node[i+1][Grid.JNu].X - 
                   Grid.Node[i  ][Grid.JNu].X);
          ds_ratio = dl/(dl+dr);
          Grid.Node[i][Grid.JNu].X = 
             Grid.Node[i-1][Grid.JNu].X +
             ds_ratio*(Grid.Node[i+1][Grid.JNu].X-
                       Grid.Node[i-1][Grid.JNu].X);
       } /* endfor */
    } /* endif */

    /* Adjust the node locations of the south boundary
       of the quadrilateral mesh block. */

    if (Fix_South_Boundary) {
       for ( i = Grid.INl+1 ; i <= Grid.INu-1 ; i+=2 ) {
          dl = abs(Grid.Node[i  ][Grid.JNl].X - 
                   Grid.Node[i-1][Grid.JNl].X);
          dr = abs(Grid.Node[i+1][Grid.JNl].X - 
                   Grid.Node[i  ][Grid.JNl].X);
          ds_ratio = dl/(dl+dr);
          Grid.Node[i][Grid.JNl].X = 
             Grid.Node[i-1][Grid.JNl].X +
             ds_ratio*(Grid.Node[i+1][Grid.JNl].X-
                       Grid.Node[i-1][Grid.JNl].X);
       } /* endfor */
    } /* endif */

    /* Adjust the node locations of the east boundary
       of the quadrilateral mesh block. */

    if (Fix_East_Boundary) {
       for ( j  = Grid.JNl+1; j <= Grid.JNu-1; j+=2 ) {
          dl = abs(Grid.Node[Grid.INu][j  ].X - 
                   Grid.Node[Grid.INu][j-1].X);
          dr = abs(Grid.Node[Grid.INu][j+1].X - 
                   Grid.Node[Grid.INu][j  ].X);
          ds_ratio = dl/(dl+dr);
          Grid.Node[Grid.INu][j].X = 
             Grid.Node[Grid.INu][j-1].X +
             ds_ratio*(Grid.Node[Grid.INu][j+1].X-
                       Grid.Node[Grid.INu][j-1].X);
       } /* endfor */
    } /* endif */

    /* Adjust the node locations of the west boundary
       of the quadrilateral mesh block. */

    if (Fix_West_Boundary) {
       for ( j  = Grid.JNl+1; j <= Grid.JNu-1; j+=2 ) {
          dl = abs(Grid.Node[Grid.INl][j  ].X - 
                   Grid.Node[Grid.INl][j-1].X);
          dr = abs(Grid.Node[Grid.INl][j+1].X - 
                   Grid.Node[Grid.INl][j  ].X);
          ds_ratio = dl/(dl+dr);
          Grid.Node[Grid.INl][j].X = 
             Grid.Node[Grid.INl][j-1].X +
             ds_ratio*(Grid.Node[Grid.INl][j+1].X-
                       Grid.Node[Grid.INl][j-1].X);
       } /* endfor */
    } /* endif */

    /* Reset the boundary condition types at the grid block
       boundaries. */
 
    Set_BCs(Grid);

    /* Recompute the exterior nodes for the quadrilateral mesh block. */

    Update_Exterior_Nodes(Grid);

    /* Recompute the cells for the quadrilateral mesh block. */

    Update_Cells(Grid);

}

/********************************************************
 * Routine: Unfix_Refined_Mesh_Boundaries               *
 *                                                      *
 * Returns the adjusted the locations of the boundary   *
 * nodes of a quadrilateral grid block to their         *
 * original unmodified positions.                       *
 *                                                      *
 ********************************************************/
void Unfix_Refined_Mesh_Boundaries(Grid2D_Quad_Block &Grid) {

    int i, j;
    double sp_l, sp_r, sp_m, ds_ratio, dl, dr;
 
    /* Return the nodes at the north boundary
       to their original positions. */

    if (Grid.BndNorthSpline.np != 0) {
       for ( i = Grid.INl+1 ; i < Grid.INu ; i += 2 ) {
           sp_l = getS(Grid.Node[i-1][Grid.JNu].X, 
                       Grid.BndNorthSpline);
           sp_r = getS(Grid.Node[i+1][Grid.JNu].X, 
                       Grid.BndNorthSpline);
           dl = abs(Grid.Node[i  ][Grid.JNu].X - 
                    Grid.Node[i-1][Grid.JNu].X);
           dr = abs(Grid.Node[i+1][Grid.JNu].X - 
                    Grid.Node[i  ][Grid.JNu].X);
           ds_ratio = dl/(dl+dr);
           sp_m = sp_l + ds_ratio*(sp_r-sp_l);
           Grid.Node[i][Grid.JNu].X = Spline(sp_m, Grid.BndNorthSpline);
       } /* endfor */
    } /* endif */

    /* Return the nodes at the south boundary
       to their original positions. */

    if (Grid.BndSouthSpline.np != 0) {
       for ( i = Grid.INl+1 ; i < Grid.INu ; i += 2 ) {
           sp_l = getS(Grid.Node[i-1][Grid.JNl].X, 
                       Grid.BndSouthSpline);
           sp_r = getS(Grid.Node[i+1][Grid.JNl].X, 
                       Grid.BndSouthSpline);
           dl = abs(Grid.Node[i  ][Grid.JNl].X - 
                    Grid.Node[i-1][Grid.JNl].X);
           dr = abs(Grid.Node[i+1][Grid.JNl].X - 
                    Grid.Node[i  ][Grid.JNl].X);
           ds_ratio = dl/(dl+dr);
           sp_m = sp_l + ds_ratio*(sp_r-sp_l);
           Grid.Node[i][Grid.JNl].X = Spline(sp_m, Grid.BndSouthSpline);
       } /* endfor */
    } /* endif */

    /* Return the nodes at the east boundary
       to their original positions. */

    if (Grid.BndEastSpline.np != 0) {
       for (j  = Grid.JNl+1 ; j < Grid.JNu ; j += 2 ) {
           sp_l = getS(Grid.Node[Grid.INu][j-1].X, 
                       Grid.BndEastSpline);
           sp_r = getS(Grid.Node[Grid.INu][j+1].X, 
                       Grid.BndEastSpline);
           dl = abs(Grid.Node[Grid.INu][j  ].X - 
                    Grid.Node[Grid.INu][j-1].X);
           dr = abs(Grid.Node[Grid.INu][j+1].X - 
                    Grid.Node[Grid.INu][j  ].X);
           ds_ratio = dl/(dl+dr);
           sp_m = sp_l + ds_ratio*(sp_r-sp_l);
           Grid.Node[Grid.INu][j].X = Spline(sp_m, Grid.BndEastSpline);
       } /* endfor */
    } /* endif */

    /* Return the nodes at the west boundary
       to their original positions. */

    if (Grid.BndWestSpline.np != 0) {
       for (j  = Grid.JNl+1 ; j < Grid.JNu ; j += 2 ) {
           sp_l = getS(Grid.Node[Grid.INl][j-1].X, 
                       Grid.BndWestSpline);
           sp_r = getS(Grid.Node[Grid.INl][j+1].X, 
                       Grid.BndWestSpline);
           dl = abs(Grid.Node[Grid.INl][j  ].X - 
                    Grid.Node[Grid.INl][j-1].X);
           dr = abs(Grid.Node[Grid.INl][j+1].X - 
                    Grid.Node[Grid.INl][j  ].X);
           ds_ratio = dl/(dl+dr);
           sp_m = sp_l + ds_ratio*(sp_r-sp_l);
           Grid.Node[Grid.INl][j].X = Spline(sp_m, Grid.BndWestSpline);
       } /* endfor */
    } /* endif */

    /* Reset the boundary condition types at the grid block
       boundaries. */
 
    Set_BCs(Grid);

    /* Recompute the exterior nodes for the quadrilateral mesh block. */

    Update_Exterior_Nodes(Grid);

    /* Recompute the cells for the quadrilateral mesh block. */

    Update_Cells(Grid);

}

