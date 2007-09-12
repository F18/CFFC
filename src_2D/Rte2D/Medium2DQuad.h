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
 *      Un      -- Return conserved variable solution 
 *                 at the specified node.              
 *    UnNW      -- Return conserved variable solution   
 *                 at the north-west node.              
 *    UnNE      -- Return conserved variable solution   
 *                 at the north-east node.              
 *    UnSW      -- Return conserved variable solution   
 *                 at the south-west node.              
 *    UnSE      -- Return conserved variable solution   
 *                 at the south-east node.              
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

  //@{ @name Bilinear interplation (Zingg & Yarrow).
  //! Return conserverd solution state at specified nodes.
  Medium2D_State Un(const int &ii, const int &jj);
  Medium2D_State UnNW(const int &ii, const int &jj);
  Medium2D_State UnNE(const int &ii, const int &jj);
  Medium2D_State UnSE(const int &ii, const int &jj);
  Medium2D_State UnSW(const int &ii, const int &jj);
  //@}

  //@{ @name Member functions required for message passing.
  //! Number of solution state variables.
  int NumVar(void) { return U[0][0].NUM_VAR_MEDIUM2D; }
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


/**************************************************************************
 * Medium2D_Quad_Block::Un -- Node conservative solution.                 *
 **************************************************************************/
inline Medium2D_State Medium2D_Quad_Block::Un(const int &ii, const int &jj) {
 
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y,
    eta1, zeta1, eta2, zeta2, eta, zeta;
 
  x=Grid.Node[ii][jj].X.x; y=Grid.Node[ii][jj].X.y;
  ax=Grid.Cell[ii-1][jj-1].Xc.x;
  bx=Grid.Cell[ii-1][jj].Xc.x-Grid.Cell[ii-1][jj-1].Xc.x;
  cx=Grid.Cell[ii][jj-1].Xc.x-Grid.Cell[ii-1][jj-1].Xc.x;
  dx=Grid.Cell[ii][jj].Xc.x+Grid.Cell[ii-1][jj-1].Xc.x-
    Grid.Cell[ii-1][jj].Xc.x-Grid.Cell[ii][jj-1].Xc.x;
  ay=Grid.Cell[ii-1][jj-1].Xc.y;
  by=Grid.Cell[ii-1][jj].Xc.y-Grid.Cell[ii-1][jj-1].Xc.y;
  cy=Grid.Cell[ii][jj-1].Xc.y-Grid.Cell[ii-1][jj-1].Xc.y;
  dy=Grid.Cell[ii][jj].Xc.y+Grid.Cell[ii-1][jj-1].Xc.y-
    Grid.Cell[ii-1][jj].Xc.y-Grid.Cell[ii][jj-1].Xc.y;
  aa=bx*dy-dx*by; bb=dy*(ax-x)+bx*cy-cx*by+dx*(y-ay); cc=cy*(ax-x)+cx*(y-ay);
  if (fabs(aa) < TOLER*TOLER) {
    if (fabs(bb) >= TOLER*TOLER) { zeta1=-cc/bb; }
    else { zeta1 = -cc/sgn(bb)*(TOLER*TOLER); }
    if (fabs(cy+dy*zeta1) >= TOLER*TOLER) { eta1=(y-ay-by*zeta1)/(cy+dy*zeta1); }
    else { eta1 = HALF; } zeta2=zeta1; eta2=eta1;
  } else {
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) { zeta1=HALF*(-bb+sqrt(bb*bb-FOUR*aa*cc))/aa; }
    else { zeta1 = -HALF*bb/aa; }
    if (fabs(cy+dy*zeta1) < TOLER*TOLER) { eta1=-ONE; }
    else { eta1=(y-ay-by*zeta1)/(cy+dy*zeta1); }
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) { zeta2=HALF*(-bb-sqrt(bb*bb-FOUR*aa*cc))/aa; }
    else { zeta2 = -HALF*bb/aa; }
    if (fabs(cy+dy*zeta2) < TOLER*TOLER) { eta2=-ONE; }
    else { eta2=(y-ay-by*zeta2)/(cy+dy*zeta2); }
  } /* end if */
  if (zeta1 > -TOLER && zeta1 < ONE + TOLER &&
      eta1  > -TOLER && eta1  < ONE + TOLER) {
    zeta=zeta1; eta=eta1;
  } else if (zeta2 > -TOLER && zeta2 < ONE + TOLER &&
	     eta2  > -TOLER && eta2  < ONE + TOLER) {
    zeta=zeta2; eta=eta2;
  } else {
    zeta=HALF; eta=HALF;
  } /* endif */
  
  return (U[ii-1][jj-1] +(U[ii-1][jj]-U[ii-1][jj-1])*zeta+  
	  (U[ii][jj-1]-U[ii-1][jj-1])*eta + 
	  (U[ii][jj]+U[ii-1][jj-1]-U[ii-1][jj]-U[ii][jj-1])*zeta*eta);

}

/**************************************************************************
 * Medium2D_Quad_Block::Un -- Get cell node conserved solution states.    *
 **************************************************************************/
inline Medium2D_State Medium2D_Quad_Block::UnNW(const int &ii, const int &jj) {
  return (Un(ii, jj+1));
}

inline Medium2D_State Medium2D_Quad_Block::UnNE(const int &ii, const int &jj) {
  return (Un(ii+1, jj+1));
}

inline Medium2D_State Medium2D_Quad_Block::UnSE(const int &ii, const int &jj) {
  return (Un(ii+1, jj));
}

inline Medium2D_State Medium2D_Quad_Block::UnSW(const int &ii, const int &jj) {
  return (Un(ii, jj));
}


/********************************************************
 * Required CFFC header files                           *
 ********************************************************/
#include "Medium2DState.h"


#endif /* _MEDIUM2D_QUAD_INCLUDED  */
