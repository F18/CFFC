/**********************************************************************
 **********************************************************************
 **                                                                  **
 ** File: Medium2DQuad.h                                             **
 **                                                                  **
 ** Description: The header file defining the 2D Medium              **
 ** Quadrialateral mesh solution classes. This is a stripped down    **
 ** version of the regular Xxx2D_Quad_Block classes.  Its purpose    **
 ** is just to hold the state of the medium on a seperate grid.      **
 **                                                                  **
 ** Author: Marc "T-Bone" Charest                                    **
 **                                                                  **
 ** Revision:  Date        Initials   Change                         **
 **            04/03/2007  MRC        Original creation              **
 **                                                                  **
 **********************************************************************
 **********************************************************************/

#ifndef _MEDIUM2D_QUAD_INCLUDED
#define _MEDIUM2D_QUAD_INCLUDED

/********************************************************
 * Class Declaration                                    *
 ********************************************************/
class Medium2D_State;


/********************************************************
 * Required CFFC header files                           *
 ********************************************************/

#include "../Grid/Cell2D.h"
#include "../Grid/Grid2DQuad.h"
#include "../AMR/QuadTree.h"




/***********************************************************************/
/*!
 * Class: Medium2D_Quad_Block
 *
 * @brief Class definition of the 2D Medium solution blocks.
 *
 * \verbatim
 * Member functions
 *       U      -- Return conserved variable solution for the block
 *                 (cell average).
 *    Grid      -- Return the solution block quadrilateral grid or mesh.
 *   ownsGrid   -- True if the grid was allocated by this class, false if 
 *                 it just points to someone elses grid.
 *   allocate   -- Allocate memory for structured quadrilateral
 *                 solution block.
 *   deallocate -- Deallocate memory for structured quadrilateral
 *                 solution block.
 * \endverbatim
 */
/***********************************************************************/

class Medium2D_Quad_Block{

private:
  bool                ownsGrid; //!< Does the object own the grid. 


public:
  //@{ @name Solution state arrays:
  Medium2D_State           **U; //!< Participating medium state.
  //@}

  //@{ @name Grid block information:
  Grid2D_Quad_Block       Grid; //!< 2D quadrilateral grid geometry.
                                //!< (false if just pointing to another grid)
  //@}

  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  Medium2D_Quad_Block(void) : U(NULL), Grid() {}

  //! Copy constructor. (use default)
  // Rte2D_Quad_Block(const Rte2D_Quad_Block &Soln);

  //! Destructor. (use default)
  // ~Medium2D_Quad_Block(void);
  //@}

  /* Assignment operator. */
  // Medium2D_Quad_Block operator = (const Medium2D_Quad_Block &Soln);
  // Use automatically generated assignment operator.

  //@{ @name Allocate and deallocate functions.
  //! Allocate memory for structured quadrilateral solution block.
  void allocate(const int Ni, const int Nj, const int Ng);
  void allocate(const Grid2D_Quad_Block &Grid_ptr);

  //! Deallocate memory for structured quadrilateral solution block.
  void deallocate(void);
  //@}

};


/**************************************************************************
 * Rte2D_Quad_Block::allocate -- Allocate memory.                         *
 **************************************************************************/
inline void Medium2D_Quad_Block::allocate(const int Ni, const int Nj, const int Ng) {

  // allocate its own grid
  ownsGrid = true;
  Grid.allocate(Ni, Nj, Ng);
  U = new Medium2D_State*[Grid.NCi];
  for ( int i = 0; i < Grid.NCi ; ++i ) U[i] = new Medium2D_State[Grid.NCj];
  
  // Set values to zero.
  for (int j  = Grid.JCl-Grid.Nghost ; j <= Grid.JCu+Grid.Nghost ; ++j )
    for ( int i = Grid.ICl-Grid.Nghost ; i <= Grid.ICu+Grid.Nghost ; ++i )
      U[i][j].Zero();

}


inline void Medium2D_Quad_Block::allocate(const Grid2D_Quad_Block &Grid_ptr) {

  // copy someones grid (NOTE: this is a soft copy -> don't deallocate this)
  ownsGrid = false;
  Grid = Grid_ptr;
  U = new Medium2D_State*[Grid.NCi];
  for ( int i = 0; i < Grid.NCi ; ++i ) U[i] = new Medium2D_State[Grid.NCj];
  
  // Set values to zero.
  for (int j  = Grid.JCl-Grid.Nghost ; j <= Grid.JCu+Grid.Nghost ; ++j )
    for ( int i = Grid.ICl-Grid.Nghost ; i <= Grid.ICu+Grid.Nghost ; ++i )
      U[i][j].Zero();

}


/**************************************************************************
 * Rte2D_Quad_Block::deallocate -- Deallocate memory.                     *
 **************************************************************************/
inline void Medium2D_Quad_Block::deallocate(void) {
  if (ownsGrid) Grid.deallocate(); 
  for ( int i = 0; i <= Grid.NCi-1 ; ++i ) delete[] U[i];
  delete[] U; U = NULL;
}


/********************************************************
 * Required CFFC header files                           *
 ********************************************************/
#include "Medium2DState.h"


#endif /* _MEDIUM2D_QUAD_INCLUDED  */
