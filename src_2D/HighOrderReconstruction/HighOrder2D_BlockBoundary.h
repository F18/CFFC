/*!\file HighOrder2D_BlockBoundary.h
   \brief Header file defining the block boundary class used in the HighOrder2D class. */

#ifndef _HIGHORDER_2D_BLOCKBOUNDARY_INCLUDED
#define _HIGHORDER_2D_BLOCKBOUNDARY_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
// None


/*!
 * \class HighOrder2D_BlockBoundary
 *
 * \brief Class for high-order block boundary in 2D (used at solution level).
 * \nosubgrouping
 *********************************************************/
class HighOrder2D_BlockBoundary{

public:
  //! @name Defined public types:
  //@{
  //!Type to classify the block boundary from the point of view of how the reconstruction is calculated.
  enum BoundaryReconstructionCalculationType { constrained, unconstrained};
  /*!
   * Type to classify the block boundary from the point of view of the influence that it has
   * on the stencils of nearby cells.
   * A transparent boundary can be crossed and thus, the cells in its proximity will have a central stencil,
   * whereas an opaque boundary cannot be crossed and thus, the cells will have special stencils.
   */
  enum BoundaryInfluenceOnCellStencilType { transparent, opaque };
  //@}

  //! @name Constructors:
  //@{
  HighOrder2D_BlockBoundary(void);		//!< Simple constructor 
  HighOrder2D_BlockBoundary(const BoundaryReconstructionCalculationType& _boundary_reconstruction_,
			    const BoundaryInfluenceOnCellStencilType& _boundary_stencil_influence_); //!< Advanced constructor
  //@}

  //! @name Get block boundary setup:
  //@{
  bool IsReconstructionConstrained(void) const { return (block_boundary_reconstruction == constrained)? true: false; }
  bool IsReconstructionUnconstrained(void) const { return (block_boundary_reconstruction == unconstrained)? true: false; }
  bool IsReconstructionStencilAffected(void) const { return (block_boundary_stencil == opaque)? true: false; }
  //@}

  //! @name Set block boundary variables:
  //@{
  void setTransparentBoundary(void){ block_boundary_stencil = transparent; }
  void setOpaqueBoundary(void){ block_boundary_stencil = opaque; }
  void setConstrainedBoundary(void){
    block_boundary_reconstruction = constrained;
    setOpaqueBoundary();	// a constrained boundary is also opaque!
  }
  void setUnconstrainedBoundary(void){ block_boundary_reconstruction = unconstrained; }
  //@}

  //! Reset object
  void Reset(void);

private:

  //! Flag to indicate the type of reconstruction for cells near this boundary
  BoundaryReconstructionCalculationType block_boundary_reconstruction;
  //! Flag to indicate the influence of this boundary on nearby cell stencils
  BoundaryInfluenceOnCellStencilType block_boundary_stencil;
  
};


#endif
