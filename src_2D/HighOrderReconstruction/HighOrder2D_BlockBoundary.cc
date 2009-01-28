/*!\file HighOrder2D_BlockBoundary.cc
   \brief File to implement member functions of HighOrder2D_BlockBoundary class. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "HighOrder2D_BlockBoundary.h"


/*! 
 * Default Constructor
 */ 
HighOrder2D_BlockBoundary::HighOrder2D_BlockBoundary(void):
  block_boundary_reconstruction(unconstrained),
  block_boundary_stencil(transparent)
{
  // Nothing
}

/*! 
 * Advanced Constructor
 */ 
HighOrder2D_BlockBoundary::HighOrder2D_BlockBoundary(const BoundaryReconstructionCalculationType& _boundary_reconstruction_,
						     const BoundaryInfluenceOnCellStencilType& _boundary_stencil_influence_):
  block_boundary_reconstruction(_boundary_reconstruction_),
  block_boundary_stencil(_boundary_stencil_influence_)
{
  // Nothing
}

/*! 
 * Reset the high-order boundary object.
 */ 
void HighOrder2D_BlockBoundary::Reset(void){
  block_boundary_reconstruction = unconstrained;
  block_boundary_stencil = transparent;
}
