/* TurbulenceModelling.h:  Header file for turbulence modelling parameters and subroutines. */

#ifndef _TURBULENCEMODELLING_INCLUDED
#define _TURBULENCEMODELLING_INCLUDED

/* Include 2D quadrilateral multiblock grid, adaptive block,
   and quadtree header files. */

#ifndef _GRID2D_QUAD_BLOCK_INCLUDED
#include "../Grid/Grid2DQuad.h"
#endif // _GRID2D_QUAD_BLOCK_INCLUDED

#ifndef _ADAPTIVBLOCK_INCLUDED
#include "../AMR/AdaptiveBlock.h"
#endif // _ADAPTIVEBLOCK_INCLUDED

#ifndef _QUADTREE_INCLUDED
#include "../AMR/QuadTree.h"
#endif // _QUADTREE_INCLUDED

/******************************************************************
 * TurbulenceModelling -- Templated subroutines.                  *
 ******************************************************************/

/**********************************************************
 * Routine: Distance_to_Wall                              *
 *                                                        *
 * Determines the normal distance to the wall for each    *
 * cell of a multiblock body-fitted quadrilateral mesh.   *
 * A direct search is performed.                          *
 *                                                        *
 **********************************************************/
template <class Quad_Soln_Block>
int Distance_to_Wall(Quad_Soln_Block             *Soln_ptr,
                     QuadTreeBlock_DataStructure &QuadTree,
                     AdaptiveBlock2D_List        &LocalSolnBlockList) {

    int iCPU, iBLK, buffer_size;
    int i, j, ilow, iup, jlow, jup, ni, nj;
    double y_wall;
    double *x_buffer;
    Vector2D x_cell;

    for ( iCPU = 0 ; iCPU <= QuadTree.Ncpu-1 ; ++iCPU ) { // Loop over available processors.
        for ( iBLK = 0 ; iBLK <= QuadTree.Nblk-1 ; ++iBLK ) { // Loop over available blocks.
            if (QuadTree.Blocks[iCPU][iBLK] != NULL) {
               if (QuadTree.Blocks[iCPU][iBLK]->block.used) {

		  // Determine dimensions of block of interest.
                  ilow = QuadTree.Blocks[iCPU][iBLK]->block.info.dimen.ghost;
                  jlow = QuadTree.Blocks[iCPU][iBLK]->block.info.dimen.ghost;
                  iup = ilow + QuadTree.Blocks[iCPU][iBLK]->block.info.dimen.i - 1;
                  jup = jlow + QuadTree.Blocks[iCPU][iBLK]->block.info.dimen.j - 1;

                  // Allocate buffer for storing cell centers.
                  ni = iup - ilow + 1;
                  nj = jup - jlow + 1;
                  x_buffer = new double[2*ni*nj];

                  // Get coordinates of cell centers.
                  if (LocalSolnBlockList.ThisCPU == iCPU) {
                     buffer_size = 0;
                     for ( j = jlow; j <= jup; ++j ) {
	  	       for (i = ilow; i <= iup; ++i ) {
     	                  x_buffer[buffer_size] = Soln_ptr[iBLK].Grid.Cell[i][j].Xc.x;
 	                  x_buffer[buffer_size+1] = Soln_ptr[iBLK].Grid.Cell[i][j].Xc.y;
                          buffer_size = buffer_size + 2;
                       } /* endfor */
                     } /* endfor */
                  } /* endif */

#ifdef _MPI_VERSION
                  buffer_size = 2*ni*nj;
                  MPI::COMM_WORLD.Bcast(x_buffer, buffer_size, MPI::DOUBLE, iCPU);
#endif

                  // Calculate normal distance to the wall.
                  buffer_size = 0;
                  for ( j = jlow; j <= jup; ++j ) {
		    for (i = ilow; i <= iup; ++i ) {
                       x_cell.x = x_buffer[buffer_size];
                       x_cell.y = x_buffer[buffer_size+1];
                       buffer_size = buffer_size + 2;

                       y_wall = Distance_to_Wall(Soln_ptr,
                                                 LocalSolnBlockList,
                                                 x_cell);

#ifdef _MPI_VERSION
                       MPI::COMM_WORLD.Allreduce(&y_wall, &y_wall, 1, MPI::DOUBLE, MPI::MIN);
#endif

                      double friction_vel = 1.587;
                      double viscosity = 1.83437E-05;
                      y_wall = friction_vel*y_wall/viscosity;

                      if (LocalSolnBlockList.ThisCPU == iCPU) {
                          Soln_ptr[iBLK].Grid.Cell[i][j].ywall = y_wall;
                          //cout << "\n " << i << " " << j << Soln_ptr[iBLK].Grid.Cell[i][j].Xc << " " << y_wall;
                       } /* endif */

                    } /* endfor */
                  } /* endfor */

                  // Deallocate buffer for storing cell centers.
                  delete []x_buffer;
                  x_buffer = NULL;

               } /* endif */
            } /* endif */
        } /* endfor */
    } /* endfor */    

    /* Distance to wall calculation complete.  Return zero value. */

    return(0);

}


/*!
 *
 * Class: SubfilterScaleModels
 *
 */
class SubfilterScaleModels{

  public:
    static double Smagorinsky_coef; //!< Smagorinsky coefficient.
    static double CI_Yoshizawa;     //!< Coefficient to calculate k based on Yoshizawa's model.
    static double CV_coef;          //!< Coefficient to estimate the eddy viscosity in the SFS k-equation.  
    static double CEPS_coef;        //!< Coefficient to calculate the dissipation rate of SFS k.


    // Two coefficients for Smagorinsky and Yoshizawa or two coefficients for k-equation.
    void set_SFSmodel_constants(const double &, const double &);

    // sqrt(2*S*S)
    double abs_strain_rate(const Tensor2D &strain_rate) const;

    //! Kinematic turbulent viscosity using Smagorinsky model
    double eddy_viscosity_Smagorinsky(const Tensor2D &strain_rate,
				      const double &filter_width) const;

    //! Eddy viscosity based on k
    double eddy_viscosity_k(const double &k, 
			    const double &filter_width) const;

    //! k based on Yoshizawa's model 
    double sfs_k_Yoshizawa(const Tensor2D &strain_rate,
			   const double &filter_width) const;
       
};


#endif /* _TURBULENCEMODELLING_INCLUDED  */
