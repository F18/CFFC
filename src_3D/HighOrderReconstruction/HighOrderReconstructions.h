/*!\file HighOrderReconstructions.h
  \brief Header file implementing the member functions of the HighOrder class which reconstructs the solution.
  \note To use the functions defined in this file include 'HighOrder.h'!
*/

#ifndef _HIGHORDER_RECONSTRUCTIONS_INCLUDED
#define _HIGHORDER_RECONSTRUCTIONS_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
// None


/* -----------------------------------------------------------------
 * =============== BLOCK LEVEL 3D RECONSTRUCTIONS ==================
 * ----------------------------------------------------------------*/

/*! 
 * Compute the unlimited k-exact high-order reconstruction
 * proposed by Barth (1993) for all block computational cells.
 *
 * \param SolnBlk the quad block for which the solution reconstruction is done.
 * \param ReconstructedSoln member function of Soln_Block_Type which returns the solution.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type> inline
void HighOrder<SOLN_STATE>::ComputeUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk,
     								     const Soln_State & 
								     (Soln_Block_Type::*ReconstructedSoln)(const int &,
													   const int &,
													   const int &) const){

  int i,j,k;

  /***************************************************************************
   *    Perform unconstrained unlimited high-order solution reconstruction   *
   **************************************************************************/

  // Set the freeze limiter flag. It's value affects the reset monotonicity data!
  _freeze_limiter = SolnBlk.Freeze_Limiter;

  /* If reconstruction with pseudo-inverse is required,
     check if the associated grid encountered modifications
     since the pseudo-inverse was computed last time. */
  if ( IsPseudoInverseAllocated() || (IsPseudoInversePreComputed() == false) ) {
       /*&& ( ObserverInteriorCellGeometryState != Geom->getInteriorStateTracker() ||
	    ObserverGhostCellGeometryState != Geom->getGhostStateTracker() || 
	    ObserverCornerGhostCellGeometryState != Geom->getCornerGhostStateTracker() */


    // Require pseudo-inverse update
    MustUpdatePseudoInverse();
    ComputeReconstructionPseudoInverse();
  }

  // Check whether constrained reconstruction is required anywhere in the block.
  //  if ( !_constrained_block_reconstruction ){

  // Carry out the solution reconstruction using the central stencil for cells in the specified range.
  for ( k  = StartK ; k <= EndK ; ++k ) {
    for ( j  = StartJ ; j <= EndJ ; ++j ) {
      for ( i = StartI ; i <= EndI ; ++i ) {
	
	// Reset the monotonicity data
	ResetMonotonicityData(i,j,k);
	
	// Set the stencil of points used for reconstruction
	SetReconstructionStencil(i, j, k, i_index, j_index, k_index);
	
	// Compute the reconstruction for the current cell
	ComputeUnconstrainedUnlimitedSolutionReconstruction(SolnBlk, ReconstructedSoln,
							    i, j, k, i_index, j_index, k_index);
      } /* endfor */
    } /* endfor */
  }/* endfor */

// --> RR: comment out higher-level execution of constrained reconstruction
//  } else {
//    
//    /*****************************************************************************************
//     *    Depending on the ReconstructionTypeMap value for a given cell,                     *
//     *    perform constrained or unconstrained unlimited high-order solution reconstruction, *
//     *    or no reconstruction at all.                                                       *
//     ****************************************************************************************/
//
//    for (j = StartJ; j <= EndJ; ++j){
//      for (i = StartI; i <= EndI; ++i){
//	
//	// Reset the monotonicity data
//	ResetMonotonicityData(i,j);
//	
//	switch (ReconstructionTypeMap[i][j]){
//	case 'r':		// "Regular reconstruction" (i.e. uses the central stencil)
//	  // Set the stencil of points used for reconstruction
//	  SetReconstructionStencil(i, j, i_index, j_index);
//
//	  // Compute the reconstruction for the current cell
//	  ComputeUnconstrainedUnlimitedSolutionReconstruction(SolnBlk, ReconstructedSoln,
//							      i, j, i_index, j_index);
//	  break;
//	  
//	case 'm':		// "Modified reconstruction" (i.e. uses a deviated stencil but it has no constraints)
//	case 'c':		// "Constrained reconstruction" (i.e. uses a deviated stencil and it has constraints)
//	  // Set the biased stencil of points used for reconstruction
//	  SetDeviatedReconstructionStencil(i, j, i_index_ave, j_index_ave, rings);
//
//	  // Compute the constrained reconstruction for the current cell
//	  ComputeConstrainedUnlimitedSolutionReconstruction(SolnBlk, ReconstructedSoln,
//							    i, j, i_index_ave, j_index_ave);
//	  break;
//	  
//	case 'n':		// "No reconstruction" (i.e. cell for which no reconstruction should be performed)
//	  // Do nothing
//	  break;
//	}
//
//      }	// endfor
//    }// endfor
//    
//  } // endif (_constrained_block_reconstruction)
}

/*! 
 * Compute the pseudo-inverse of the left-hand-side term in the 
 * unlimited k-exact high-order reconstruction for all computational
 * cells based on the information provided by the associated grid.
 * For cells influenced by the presence of constrained boundary
 * conditions, either the pseudo-inverse of the noncentral reconstruction
 * is computed (if no constraints are imposed), or the LHS matrix 
 * generated by the conservation of mean quantities in the supporting stencil.
 * 
 */
template<class SOLN_STATE> inline
void HighOrder<SOLN_STATE>::ComputeReconstructionPseudoInverse(void){

  int i,j,k;

  // == Check if the pseudo-inverse has been allocated and it hasn't been precomputed
  if ( IsPseudoInverseAllocated() && !IsPseudoInversePreComputed() ){

    // == Check if the reconstruction polynomial is piecewise constant
    if (RecOrder() == 0){
      // There is no need to calculate pseudo-inverse
      // Confirm the pseudo-inverse calculation
      _calculated_psinv = true;

      // Memorise the corresponding state of the grid
//      ObserverInteriorCellGeometryState = Geom->getInteriorStateTracker();
//      ObserverGhostCellGeometryState = Geom->getGhostStateTracker();
//      ObserverCornerGhostCellGeometryState = Geom->getCornerGhostStateTracker();
      return;
    }


    // Check whether constrained reconstruction is required anywhere in the block.
    //    if ( !_constrained_block_reconstruction ){

      // Calculate the pseudo-inverse using the central stencil for cells in the specified range
    for ( k  = StartK ; k <= EndK ; ++k ) {
      for ( j  = StartJ ; j <= EndJ ; ++j ) {
	for ( i = StartI ; i <= EndI ; ++i ) {
	
	  // Set the stencil of points used for reconstruction
	  SetReconstructionStencil(i, j, k, i_index, j_index, k_index);

	  // Compute the pseudo-inverse for the current cell
	  ComputeCellReconstructionPseudoInverse(i, j, k, i_index, j_index, k_index);
	}/* endfor */
      }/* endfor */
    }/* endfor */

//    } else {
//
//
//      // --> RR: comment out constrained pseudo inverse stuff
//
//      // Calculate the pseudo-inverse for a block that has some constrained cells.
//      // The ReconstructionTypeMap variable indicates which cells are constrained, 
//      // which once have only the stencil modified and which ones use the central stencil.
//
//      /********************************************************************************************
//       * Calculate the matrices that need to be stored in order to speed up the high-order solution
//       * reconstruction of the cells affected by the presence of constrained boundaries
//       * (i.e. calculate the pseudo-inverse of the unconstrained unlimited reconstruction
//       * for those cells that don't have constraints but use a biased supporting stencil,
//       * calculate the LHS matrix associate with the least-squares problem and the mean conservation
//       * for those cells that have also constraint equations.
//       ********************************************************************************************/
//
//      for ( j  = StartJ ; j <= EndJ ; ++j ) {
//	for ( i = StartI ; i <= EndI ; ++i ) {
//	  
//	  switch (ReconstructionTypeMap[i][j]){
//	    
//	  case 'r':		// "Regular reconstruction" (i.e. uses the central stencil)
//	    // Set the stencil of points used for reconstruction
//	    SetReconstructionStencil(i, j, i_index, j_index);
//	    
//	    // Compute the pseudo-inverse for the current cell
//	    ComputeCellReconstructionPseudoInverse(i, j, i_index, j_index);
//	    break;
//
//	  case 'm':		// "Modified reconstruction" (i.e. uses a deviated stencil but it has no constraints)
//	  case 'c':		// "Constrained reconstruction" (i.e. uses a deviated stencil and it has constraints)
//	    // Get matrix for the current cell (i.e. pseudo-inverse or LHS)
//	    ComputeCellReconstructionPseudoInverseNearConstrainedBoundaries(i,j);
//	    break;
//
//	  case 'n':		// "No reconstruction" (i.e. cell for which no reconstruction should be performed)
//	    // Do nothing
//	    break;
//	  }
//
//	}/* endfor */
//      }/* endfor */
//      
//    } // endif (_constrained_block_reconstruction)

    
    // Confirm the pseudo-inverse calculation
    _calculated_psinv = true;

    // Memorise the corresponding state of the grid
//    ObserverInteriorCellGeometryState = Geom->getInteriorStateTracker();
//    ObserverGhostCellGeometryState = Geom->getGhostStateTracker();
//    ObserverCornerGhostCellGeometryState = Geom->getCornerGhostStateTracker();

  } // endif ( IsPseudoInverseAllocated() && !IsPseudoInversePreComputed() )
  
}

/*! 
 * Compute the limited linear least-squares reconstruction proposed
 * by Barth (1993) for those interpolants detected as non-smooth.
 * This reconstruction is carried out when the order or reconstruction
 * is required to be dropped since the high-order interpolant is detected
 * to be non-smooth. \n
 * The high-order interpolant is going to be overwritten by the low-order one.
 * \param [in] SolnBlk The solution block which provides solution data.
 * \param ReconstructedSoln member function of Soln_Block_Type which returns the solution.
 * \param [in] Limiter The limiter used during this limited reconstruction.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type> inline
void HighOrder<SOLN_STATE>::EnforceMonotonicityToNonSmoothInterpolants(Soln_Block_Type &SolnBlk,
									 const int &Limiter,
									 const Soln_State &
									 (Soln_Block_Type::*ReconstructedSoln)(const int &,
													       const int &,
													       const int &) const){
//  // --> RR: comment out contents (empty function) of ENFORCE MONOTONICITY (REDUCTION TO PEICEWISE CONSTANT) !!
// ------------- comment back in when you want to to drop order ---------------//
//  // Set local variables
//  int i,j,k;
//  
//  // == Check if the reconstruction polynomial is piecewise constant
//  if (RecOrder() == 0){
//    // There is no need to enforce any monotonicity
//    return;
//  }
//
//  if (CENO_Execution_Mode::CENO_DROP_ORDER){
//    // Carry on actions required to enforce monotonicity
//
//    // Switch to limited piecewise linear reconstruction those interior interpolants detected as non-smooth
//    for ( k  = KCl ; k <= KCu ; ++k ) {
//      for ( j  = JCl ; j <= JCu ; ++j ) {
//	for ( i = ICl ; i <= ICu ; ++i ) {
//	  
//	  if ( IsThereAnyNonSmoothHighOrderReconstruction(i,j,k) ){
//	    // One or more solution variables need to have the interpolant switched to a limited piecewise linear one.
//	    ComputeLimitedPiecewiseLinearSolutionReconstruction(SolnBlk,
//								i,j,k,
//								Limiter,
//								ReconstructedSoln);
//	  } // endif
// 
//	} /* endfor(i) */
//      } /* endfor(j) */
//    } /* endfor(k) */
//  } // endif(CENO_Execution_Mode::CENO_DROP_ORDER)
// --------------------------------------------------------------------------------//


    // --> RR: comment out constrained stuff in enforce monotonicity
//    // Check whether reconstruction based flux calculation is required anywhere in the block.
//    if ( !_constrained_block_reconstruction ){
//
//      // Switch to limited piecewise linear reconstruction those ghost cells 
//      // involved in flux calculation and detected as non-smooth
//
//      // == South and North boundaries ==
//      for (i = ICl; i<=ICu; ++i){
//	// == South bnd.
//	if ( IsThereAnyNonSmoothHighOrderReconstruction(i,JCl-1) ){
//	  // One or more solution variables need to have the interpolant switched to a limited piecewise linear one.
//	  ComputeLimitedPiecewiseLinearSolutionReconstruction(SolnBlk,
//							      i,JCl-1,
//							      Limiter,
//							      ReconstructedSoln);
//	}
//	
//	// == North bnd.
//	if ( IsThereAnyNonSmoothHighOrderReconstruction(i,JCu+1) ){
//	  // One or more solution variables need to have the interpolant switched to a limited piecewise linear one.
//	  ComputeLimitedPiecewiseLinearSolutionReconstruction(SolnBlk,
//							      i,JCu+1,
//							      Limiter,
//							      ReconstructedSoln);
//	}
//      }	// endfor
//
//      // == West and East boundaries ==
//      for (j = JCl; j<=JCu; ++j){
//	// == West bnd.
//	if ( IsThereAnyNonSmoothHighOrderReconstruction(ICl-1,j) ){
//	  // One or more solution variables need to have the interpolant switched to a limited piecewise linear one.
//	  ComputeLimitedPiecewiseLinearSolutionReconstruction(SolnBlk,
//							      ICl-1,j,
//							      Limiter,
//							      ReconstructedSoln);
//	}
//	
//	// == East bnd.
//	if ( IsThereAnyNonSmoothHighOrderReconstruction(ICu+1,j) ){
//	  // One or more solution variables need to have the interpolant switched to a limited piecewise linear one.
//	  ComputeLimitedPiecewiseLinearSolutionReconstruction(SolnBlk,
//							      ICu+1,j,
//							      Limiter,
//							      ReconstructedSoln);
//	}
//      }	// endfor
//      
//    } else {
//
//      // === Some boundaries require reconstruction based flux calculation
//
//      /* Motivation of the algorithm below:
//	 If reconstruction based flux calculation is required at some of the boundaries
//	 and non-smooth solution interpolants are detected near these boundaries the 
//	 flux is not going to be computed based on the high-order interpolant but on
//	 solving a Riemann problem at the interface.
//	 The purpose of the algorithm below is to ensure that a limited piecewise linear
//	 reconstruction is available in the first ghost cells that have interface with an
//	 interior cell detected with inadequate interpolant. Thus, when the flux calculation 
//	 is performed, the Riemann problem for those interfaces can be solved.
//	 Note that no interior cells are going to be affected by the code that follows!
//	 Note also that trying to obtain a high-order interpolant in these ghost cells is not
//	 justified based on accuracy and computational efficiency reasons.
//      */
//
//      // Check WEST boundary
//      if (WestBnd.IsReconstructionConstrained()){
//	for (j = JCl; j <= JCu; ++j){
//	  if ( IsThereAnyNonSmoothHighOrderReconstruction(ICl,j) ) { // check the interior cell
//	    // flag all reconstructions of the adjacent ghost cell as non-smooth
//	    FlagCellReconstructionsAsNonSmooth(ICl-1,j);
//	    // perform a limited piecewise linear reconstruction
//	    ComputeLimitedPiecewiseLinearSolutionReconstruction(SolnBlk,
//								ICl-1,j,
//								Limiter,
//								ReconstructedSoln);
//	  }// endif
//	}// enfor 
//
//      } else {
//	// Check the ghost cells near this boundary
//	for (j = JCl; j<=JCu; ++j){
//	  if ( IsThereAnyNonSmoothHighOrderReconstruction(ICl-1,j) ){
//	    // One or more solution variables need to have the interpolant switched to a limited piecewise linear one.
//	    ComputeLimitedPiecewiseLinearSolutionReconstruction(SolnBlk,
//								ICl-1,j,
//								Limiter,
//								ReconstructedSoln);
//	  }
//	} // endfor
//
//      }// endif (WestBnd)
//
//      // Check EAST boundary
//      if (EastBnd.IsReconstructionConstrained()){
//	for (j = JCl; j <= JCu; ++j){
//	  if ( IsThereAnyNonSmoothHighOrderReconstruction(ICu,j) ) { // check the interior cell
//	    // flag all reconstructions of the adjacent ghost cell as non-smooth
//	    FlagCellReconstructionsAsNonSmooth(ICu+1,j);
//	    // perform a limited piecewise linear reconstruction
//	    ComputeLimitedPiecewiseLinearSolutionReconstruction(SolnBlk,
//								ICu+1,j,
//								Limiter,
//								ReconstructedSoln);
//	  }// endif
//	}// enfor 
//
//      } else {
//	// Check the ghost cells near this boundary	
//	for (j = JCl; j <= JCu; ++j){
//	  if ( IsThereAnyNonSmoothHighOrderReconstruction(ICu+1,j) ){
//	    // One or more solution variables need to have the interpolant switched to a limited piecewise linear one.
//	    ComputeLimitedPiecewiseLinearSolutionReconstruction(SolnBlk,
//								ICu+1,j,
//								Limiter,
//								ReconstructedSoln);
//	  }
//	} // endfor
//
//      }// endif (EastBnd)
//
//      // Check NORTH boundary
//      if (NorthBnd.IsReconstructionConstrained()){
//	for (i = ICl; i <= ICu; ++i){
//	  if ( IsThereAnyNonSmoothHighOrderReconstruction(i,JCu) ) { // check the interior cell
//	    // flag all reconstructions of the adjacent ghost cell as non-smooth
//	    FlagCellReconstructionsAsNonSmooth(i,JCu+1);
//	    // perform a limited piecewise linear reconstruction
//	    ComputeLimitedPiecewiseLinearSolutionReconstruction(SolnBlk,
//								i,JCu+1,
//								Limiter,
//								ReconstructedSoln);
//	  }// endif
//	}// enfor 
//
//      } else {
//	// Check the ghost cells near this boundary	
//	for (i = ICl; i <= ICu; ++i){ 
//	  if ( IsThereAnyNonSmoothHighOrderReconstruction(i,JCu+1) ){
//	    // One or more solution variables need to have the interpolant switched to a limited piecewise linear one.
//	    ComputeLimitedPiecewiseLinearSolutionReconstruction(SolnBlk,
//								i,JCu+1,
//								Limiter,
//								ReconstructedSoln);
//	  }
//	} // endfor
//
//      }// endif (NorthBnd)
//
//      // Check SOUTH boundary
//      if (SouthBnd.IsReconstructionConstrained()){
//	for (i = ICl; i <= ICu; ++i){
//	  if ( IsThereAnyNonSmoothHighOrderReconstruction(i,JCl) ) { // check the interior cell
//	    // flag all reconstructions of the adjacent ghost cell as non-smooth
//	    FlagCellReconstructionsAsNonSmooth(i,JCl-1);
//	    // perform a limited piecewise linear reconstruction
//	    ComputeLimitedPiecewiseLinearSolutionReconstruction(SolnBlk,
//								i,JCl-1,
//								Limiter,
//								ReconstructedSoln);
//	  }// endif
//	}// enfor 
//
//      } else {
//	// Check the ghost cells near this boundary	
//	for (i = ICl; i <= ICu; ++i){ 
//	  if ( IsThereAnyNonSmoothHighOrderReconstruction(i,JCl-1) ){
//	    // One or more solution variables need to have the interpolant switched to a limited piecewise linear one.
//	    ComputeLimitedPiecewiseLinearSolutionReconstruction(SolnBlk,
//								i,JCl-1,
//								Limiter,
//								ReconstructedSoln);
//	  }
//	} // endfor
//
//      }// endif (SouthBnd)
//
//    } // endif (_constrained_block_reconstruction)
//
//
//  } else {
//    /* Reset monotonicity flags for boundary cells near splines which require reconstruction based flux calculation.
//       If Riemann based flux calculation is desired there is nothing to be done.
//    */
//
//    // Check whether reset of monotonicity flags is required anywhere in the block.
//    if ( !_constrained_block_reconstruction ){
//      // No need to reset
//      return;
//    }
//
//    // Check WEST boundary
//    if (WestBnd.IsReconstructionConstrained()){
//      for (j = JCl - 1; j <= JCu + 1; ++j){
//	// reset the monotonicity flag for the interior cell
//	ResetMonotonicityData(ICl,j);
//      }
//    } 
//
//    // Check EAST boundary
//    if (EastBnd.IsReconstructionConstrained()){
//      for (j = JCl - 1; j <= JCu + 1; ++j){
//	// reset the monotonicity flag for the interior cell
//	ResetMonotonicityData(ICu,j);
//      }
//    } 
//
//    // Check NORTH boundary
//    if (NorthBnd.IsReconstructionConstrained()){
//      for (i = ICl - 1; i <= ICu + 1; ++i){
//	// reset the monotonicity flag for the interior cell
//	ResetMonotonicityData(i,JCu);
//      }
//    } 
//
//    // Check SOUTH boundary
//    if (SouthBnd.IsReconstructionConstrained()){
//      for (i = ICl - 1; i <= ICu + 1; ++i){
//	// reset the monotonicity flag for the interior cell
//	ResetMonotonicityData(i,JCl);
//      }
//    }

//  } // endif(CENO_Execution_Mode::CENO_DROP_ORDER)
  
}

/*! 
 * Compute the high-order reconstruction for all SolnBlk cells 
 * using the CENO algorithm proposed by Ivan and Groth (AIAA-2007-4323-670).
 * This algorithm consists of three steps:
 *      --> Perform k-exact reconstruction in each computation cell for all solution variables.
 *      --> Compute the smoothness indicator for each computation cell and solution variable.
 *      --> Switch to a limited piecewise linear reconstruction those interpolants detected as non-smooth.
 * Slightly different variants of the same idea can be employed to 
 * compute the smoothness indicator. To choose between different algorithms
 * use the control execution flags provided by CENO_Execution_Mode class.
 *
 * \note This routine performs all three steps of the algorithm, 
 *       which might not be suited for equations with both elliptic and hyperbolic terms!
 *
 * \param SolnBlk the quad block for which the solution reconstruction is done.
 * \param ReconstructedSoln member function of Soln_Block_Type which returns the solution.
 * \param Limiter flag to indicate which limiter is used in the limited piecewise linear reconstruction.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type> inline
void HighOrder<SOLN_STATE>::ComputeHighOrderSolutionReconstruction(Soln_Block_Type &SolnBlk,
								     const int &Limiter,
								     const Soln_State & 
								     (Soln_Block_Type::*ReconstructedSoln)(const int &,
													   const int &,
													   const int &) const ){


  // Step 1. Compute the unlimited solution reconstruction in all required computational cells.
  ComputeUnlimitedSolutionReconstruction(SolnBlk,
					 ReconstructedSoln);


  // Step 2. Perform smoothness indicator calculation and analysis (i.e. flag those interpolants detected as non-smooth).
  ComputeSmoothnessIndicator(SolnBlk,
			     ReconstructedSoln);

  // Step 3. Enforce monotonicity to the non-smooth interpolants (i.e. switch to a limited linear interpolant)
  EnforceMonotonicityToNonSmoothInterpolants(SolnBlk,
					     Limiter,
					     ReconstructedSoln);

}

/* -----------------------------------------------------------------
 * ================ CELL LEVEL 3D RECONSTRUCTIONS ==================
 * ----------------------------------------------------------------*/

/*! 
 * Compute the unlimited k-exact high-order reconstruction
 * proposed by Barth (1993) for a specified computational cell.
 * This reconstruction doesn't account for any constrain
 * other than the mean quantity conservation.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type>
void HighOrder<SOLN_STATE>::
ComputeUnconstrainedUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk,
						    const Soln_State & 
						    (Soln_Block_Type::*ReconstructedSoln)(const int &,const int &,const int &) const,
						    const int &iCell, const int &jCell, const int &kCell,
						    const IndexType & i_index,
						    const IndexType & j_index,
						    const IndexType & k_index){

  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS

  int StencilSize(i_index.size()); 
  int cell, i, parameter;
  int P1, P2, P3;

  // *********  Assign the average solution to D00 ***********
  CellTaylorDerivState(iCell,jCell,kCell,0) = (SolnBlk.*ReconstructedSoln)(iCell,jCell,kCell);

  // == Check if the reconstruction polynomial is piecewise constant
  if (RecOrder() == 0){
    // There is no need to calculate the reconstruction
    return;
  }

  // Check if the pseudo-inverse has been allocated and pre-computed
  if ( IsPseudoInverseAllocated() && IsPseudoInversePreComputed() ){

    // Use the pseudo-inverse to calculate the least-squares reconstruction
    
    // Ensure that Delta_U vector has proper dimension
    if ( Delta_U.size() != (StencilSize - 1) ){
      Delta_U.newsize(StencilSize - 1);
    }

    // START: Compute for every parameter the high-order approximation
    // ***************************************************************
    for (parameter = 1; parameter <= NumberOfVariables() ; ++parameter){

      // Step 1. SET the vector Delta_U of the linear system (RHS) for the current parameter ***
      for (cell = 1 ; cell < StencilSize; ++cell) { // for each neighbour cell in the stencil

	// Compute Delta_U = U[neighbour] - U[cell] for each parameter
	Delta_U(cell-1) = ( (SolnBlk.*ReconstructedSoln)(i_index[cell],j_index[cell],k_index[cell])[parameter] -
			    (SolnBlk.*ReconstructedSoln)(iCell,jCell,kCell)[parameter] );

	// Apply the precomputed geometric weight to the Delta_U term
	Delta_U(cell-1) *= GeomWeightValue(iCell,jCell,kCell,cell);
      } 
     
      // Step 2. Find the solution of the linear-system for the current parameter
      X = Cell_LHS_Inv(iCell,jCell,kCell) * Delta_U;
      
      // Step 3. Update the high-order derivatives for the current parameter
      for (i = 1; i <= CellTaylorDeriv(iCell,jCell,kCell).LastElem(); ++i){

	// Identify 'x' and 'y' powers of the i-th derivative
	P1 = CellTaylorDeriv(iCell,jCell,kCell,i).P1();  // identify P1
	P2 = CellTaylorDeriv(iCell,jCell,kCell,i).P2();  // identify P2
	P3 = CellTaylorDeriv(iCell,jCell,kCell,i).P3();  // identify P3

	// Set the i-th derivative for the current parameter
	CellTaylorDeriv(iCell,jCell,kCell,i).D(parameter) = X(i-1);
	
	// This equation ensures the mean conservation of the current parameter inside the reconstructed cell.
	CellTaylorDeriv(iCell,jCell,kCell,0).D(parameter) -= Geom->CellGeomCoeffValue(iCell,jCell,kCell,P1,P2,P3) * X(i-1);
      }
    }
    // STOP: Reconstruction solution (i.e. derivatives) obtained.
    // ************************************************************

  } else {
    
    // Form both LHS and RHS in order to calculate the least-squares reconstruction

    // SET VARIABLES USED ONLY IN THIS RECONSTRUCTION PROCESS

    int krank;                        //< the final rank of matrix A is returned here
    int IndexSumZ, IndexSumY, IndexSumX, P1, P2, P3;
    double CombP1X, CombP2Y, CombP3Z;
    double PowDistanceZC,PowDistanceYC, PowDistanceXC;
    double MaxWeight(0.0);
    double IntSum1(0.0),IntSum2(0.0);

    // Ensure that enough memory is allocated for the current least-squares problem.
    if (A.size(0) != StencilSize-1) {
      // Resize the matrices accordingly
      A.newsize(StencilSize-1,NumberOfTaylorDerivatives()-1);
      All_Delta_U.newsize(StencilSize-1,NumberOfVariables());
    }

    // START:   Set the LHS and RHS of the linear system 
    // ***************************************************

    // ==== Set the geometric weight associated with the reconstructed cell
    GeometricWeights[0] = 1;

    // Step1. Compute the normalized geometric weights
    for (cell=1; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil

      /* Compute the X, Y, and Z component of the distance between
	 the cell centers of the neighbour and the reconstructed cell */
      DeltaCellCenters[cell] = CellCenter(i_index[cell],j_index[cell],k_index[cell]) - CellCenter(iCell,jCell,kCell);
    
      /* Compute the geometric weight based on the centroid distance */
      CENO_Geometric_Weighting(GeometricWeights[cell], DeltaCellCenters[cell].abs());

      /* Compute the maximum geometric weight (this is used for normalization) */
      MaxWeight = max(MaxWeight, GeometricWeights[cell]);
    }

    // Step2. Set the approximate equations
    for (cell=1 ; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil
    
      // compute the normalized geometric weight
      GeometricWeights[cell] /= MaxWeight;

      // *** SET the matrix A of the linear system (LHS) ***
      /* compute for each derivative the corresponding entry in the matrix of the linear system */
      for (i=1; i<=CellTaylorDeriv(iCell,jCell,kCell).LastElem(); ++i){
	// build the row of the matrix
	P1 = CellTaylorDeriv(iCell,jCell,kCell,i).P1();  // identify P1
	P2 = CellTaylorDeriv(iCell,jCell,kCell,i).P2();  // identify P2
	P3 = CellTaylorDeriv(iCell,jCell,kCell,i).P3();  // identify P3

	A(cell-1,i-1) = 0.0;  // set sumation variable to zero
	CombP3Z = 1.0;        // the binomial coefficient "nC k" for k=0 is 1
	PowDistanceZC = 1.0;  // initialize PowDistanceZC
	//	CombP2Y = 1.0;        // the binomial coefficient "nC k" for k=0 is 1
	//	PowDistanceYC = 1.0;  // initialize PowDistanceYC

	// Compute geometric integral over the neighbour's domain
	for (IndexSumZ = 0; IndexSumZ<=P3; ++IndexSumZ){
	  CombP2Y = 1.0;        // the binomial coefficient "nC k" for k=0 is 1
	  PowDistanceYC = 1.0;  // initialize PowDistanceYC
	  IntSum2 = 0.0;         // reset internal summation variable

	  for (IndexSumY = 0; IndexSumY<=P2; ++IndexSumY){
	    CombP1X = 1.0;       // the binomial coefficient "nC k" for k=0 is 1
	    PowDistanceXC = 1.0; // initialize PowDistanceXC
	    IntSum1 = 0.0;	     // reset internal sumation variable

	    for (IndexSumX = 0; IndexSumX<=P1; ++IndexSumX){
	      IntSum1 += ( CombP1X*PowDistanceXC*
			  Geom->CellGeomCoeffValue(i_index[cell],j_index[cell],k_index[cell],P1-IndexSumX,P2-IndexSumY,P3-IndexSumZ) );
	      
	      // update the binomial coefficients
	      CombP1X = (P1-IndexSumX)*CombP1X/(IndexSumX+1); // the index is still the old one => expression for "nC k+1"
	      PowDistanceXC *= DeltaCellCenters[cell].x;      // Update PowDistanceXC
	    }//endfor
	    
	    IntSum2 += CombP2Y*PowDistanceYC*IntSum1;
	    CombP2Y = (P2-IndexSumY)*CombP2Y/(IndexSumY+1); // the index is still the old one => expression for "nC k+1"
	    PowDistanceYC *= DeltaCellCenters[cell].y;      // Update PowDistanceYC
	  }//endfor

	  A(cell-1,i-1) += CombP3Z*PowDistanceZC*IntSum2;  // update the external sum
	  
	  CombP3Z = (P3-IndexSumZ)*CombP3Z/(IndexSumZ+1); // the index is still the old one => expression for "nC k+1"
	  PowDistanceZC *= DeltaCellCenters[cell].z;      // Update PowDistanceYC
	}//endfor

	// subtract the corresponding geometric moment of cell (iCell,jCell) 
	A(cell-1,i-1) -= Geom->CellGeomCoeffValue(iCell,jCell,kCell,P1,P2,P3);

	// apply geometric weighting
	A(cell-1,i-1) *= GeometricWeights[cell];
      } //endfor (i)

      // *** SET the matrix All_Delta_U of the linear system (RHS) ***
      for (parameter = 1; parameter <= NumberOfVariables(); ++parameter){

	// Compute Delta_U = U[neighbour] - U[cell] for each parameter
	All_Delta_U(cell-1,parameter-1) = ( (SolnBlk.*ReconstructedSoln)(i_index[cell],j_index[cell],k_index[cell])[parameter] -
					    (SolnBlk.*ReconstructedSoln)(iCell,jCell,kCell)[parameter] );
	
	// Apply geometric weighting
	All_Delta_U(cell-1,parameter-1) *= GeometricWeights[cell];
      }	// endfor (parameter)
      
    }//endfor (cell)

    // STOP:   Matrix A of the linear system (LHS) built.
    //         Matrix All_Delta_U of the linear system (RHS) built.
    // **********************************************************************

    if (CENO_Execution_Mode::USE_LAPACK_LEAST_SQUARES) {
      
      // Solve the least-squares system with Lapack subroutine
      /*********************************************************/
      Solve_LS_Householder_F77(A, All_Delta_U, krank, NumberOfVariables(), StencilSize-1, NumberOfTaylorDerivatives()-1);

      // Update the high-order derivatives
      //***********************************
      for (i = 1; i <= CellTaylorDeriv(iCell,jCell,kCell).LastElem(); ++i){

	// Identify 'x' and 'y' powers of the i-th derivative
	P1 = CellTaylorDeriv(iCell,jCell,kCell,i).P1();  // identify P1
	P2 = CellTaylorDeriv(iCell,jCell,kCell,i).P2();  // identify P2
	P3 = CellTaylorDeriv(iCell,jCell,kCell,i).P3();  // identify P3

	for (parameter = 1; parameter <= NumberOfVariables(); ++parameter){
	  // Set the i-th derivative for the current parameter
	  CellTaylorDeriv(iCell,jCell,kCell,i).D(parameter) = All_Delta_U(i-1,parameter-1);
	  
	  // This equation ensures the mean conservation of the current parameter inside the reconstructed cell.
	  CellTaylorDeriv(iCell,jCell,kCell,0).D(parameter) -= Geom->CellGeomCoeffValue(iCell,jCell,kCell,P1,P2,P3) * All_Delta_U(i-1,parameter-1);
	} // endfor (parameter)
      }	// endfor (i)

    }
    // --> RR: Comment out Solve_LS_Householder call in ComputeUnconstrainedUnlimitedSolutionReconstruction
//else { 
//
//  ColumnVector Rnorm(NumberOfVariables());       //< store the residual norm of the LS problem for each parameter.
//  DenseMatrix Xm(NumberOfTaylorDerivatives()-1, NumberOfVariables()); //< store the solution to the least-square problem
//
//  /* Solve the overdetermined linear system of equations using the internal least-squares procedure */
//  /**************************************************************************************************/
//  Solve_LS_Householder(A,All_Delta_U,Xm,krank,Rnorm);
//
//  // Update the high-order derivatives
//  //***********************************
//  for (i = 1; i <= CellTaylorDeriv(iCell,jCell,kCell).LastElem(); ++i){
//
//	// Identify 'x','y', and 'z' powers of the i-th derivative
//	P1 = CellTaylorDeriv(iCell,jCell,kCell,i).P1();  // identify P1
//	P2 = CellTaylorDeriv(iCell,jCell,kCell,i).P2();  // identify P2
//	P3 = CellTaylorDeriv(iCell,jCell,kCell,i).P3();  // identify P3
//
//	for (parameter = 1; parameter <= NumberOfVariables(); ++parameter){
//	  // Set the i-th derivative for the current parameter
//	  CellTaylorDeriv(iCell,jCell,kCell,i).D(parameter) = Xm(i-1,parameter-1);
//
//	  // This equation ensures the mean conservation of the current parameter inside the reconstructed cell.
//	  CellTaylorDeriv(iCell,jCell,kCell,0).D(parameter) -= Geom->CellGeomCoeffValue(iCell,jCell,kCell,P1,P2,P3) * Xm(i-1,parameter-1);
//
//	} // endfor (parameter)
//  } // endfor (i)
//
//} // endif (CENO_Execution_Mode::USE_LAPACK_LEAST_SQUARES)

  } // endif (IsPseudoInverseAllocated() && IsPseudoInversePreComputed())

}

/*! 
 * Compute the pseudo-inverse of the left-hand-side term in the 
 * unlimited k-exact high-order reconstruction for a specified
 * computational cell based on the information provided by the
 * associated grid.
 */
template<class SOLN_STATE>
void HighOrder<SOLN_STATE>::ComputeCellReconstructionPseudoInverse(const int &iCell, 
								   const int &jCell, 
								   const int &kCell,
								   const IndexType & i_index,
								   const IndexType & j_index,
								   const IndexType & k_index){
  
  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS

  int StencilSize(i_index.size());
  int IndexSumZ, IndexSumY, IndexSumX, P1, P2, P3;
  double CombP1X, CombP2Y, CombP3Z;
  double  PowDistanceZC, PowDistanceYC, PowDistanceXC;
  int cell, i;
  double MaxWeight(0.0);
  double IntSum1(0.0), IntSum2(0.0);

  // Ensure that the LHS matrix is formated correctly.
  // Memory shouldn't be allocated here, only the dimensions should be defined properly.
  Cell_LHS_Inv(iCell,jCell,kCell).newsize(StencilSize - 1, NumberOfTaylorDerivatives() - 1);
  GeomWeights(iCell,jCell,kCell).resize(StencilSize);


  // START:   Set the LHS of the linear system
  // ***************************************************

  // ==== Set the geometric weight associated with the reconstructed cell
  GeomWeightValue(iCell,jCell,kCell,0) = 1;

  // Step1. Compute the normalized geometric weights
  for (cell=1; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil

    /* Compute the X, Y, and Z components of the distance between
       the cell centers of the neighbour and the reconstructed cell */
    DeltaCellCenters[cell] = CellCenter(i_index[cell],j_index[cell],k_index[cell]) - CellCenter(iCell,jCell,kCell);
    
    /* Compute the geometric weight based on the centroid distance */
    CENO_Geometric_Weighting(GeomWeightValue(iCell,jCell,kCell,cell), DeltaCellCenters[cell].abs());

    /* Compute the maximum geometric weight (this is used for normalization) */
    MaxWeight = max(MaxWeight, GeomWeightValue(iCell,jCell,kCell,cell));
  }

  // Step2. Set the approximate equations
  for (cell=1 ; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil
    
    // compute the normalized geometric weight
    GeomWeightValue(iCell,jCell,kCell,cell) /= MaxWeight;

    // *** SET the matrix of the linear system (LHS) ***
    /* compute for each derivative the corresponding entry in the matrix of the linear system */
    for (i=1; i<=CellTaylorDeriv(iCell,jCell,kCell).LastElem(); ++i){
      // build the row of the matrix
      P1 = CellTaylorDeriv(iCell,jCell,kCell,i).P1();  // identify P1
      P2 = CellTaylorDeriv(iCell,jCell,kCell,i).P2();  // identify P2
      P3 = CellTaylorDeriv(iCell,jCell,kCell,i).P3();  // identify P3

      //----------------------------------------------
      Cell_LHS_Inv_Value(iCell,jCell,kCell,cell-1,i-1) = 0.0;  // set sumation variable to zero
      CombP3Z = 1.0;        // the binomial coefficient "nC k" for k=0 is 1
      PowDistanceZC = 1.0;  // initialize PowDistanceZC
      
      // Compute geometric integral over the neighbour's domain
      for (IndexSumZ = 0; IndexSumZ<=P3; ++IndexSumZ){
	CombP2Y = 1.0;        // the binomial coefficient "nC k" for k=0 is 1
	PowDistanceYC = 1.0;  // initialize PowDistanceYC
	IntSum2 = 0.0;         // reset internal summation variable
	
	for (IndexSumY = 0; IndexSumY<=P2; ++IndexSumY){
	  CombP1X = 1.0;       // the binomial coefficient "nC k" for k=0 is 1
	  PowDistanceXC = 1.0; // initialize PowDistanceXC
	  IntSum1 = 0.0;	     // reset internal sumation variable
	  
	  for (IndexSumX = 0; IndexSumX<=P1; ++IndexSumX){
	    IntSum1 += ( CombP1X*PowDistanceXC*
			 Geom->CellGeomCoeffValue(i_index[cell],j_index[cell],k_index[cell],P1-IndexSumX,P2-IndexSumY,P3-IndexSumZ) );
	    
	    // update the binomial coefficients
	    CombP1X = (P1-IndexSumX)*CombP1X/(IndexSumX+1); // the index is still the old one => expression for "nC k+1"
	    PowDistanceXC *= DeltaCellCenters[cell].x;      // Update PowDistanceXC
	  }//endfor
	  
	  IntSum2 += CombP2Y*PowDistanceYC*IntSum1;
	  CombP2Y = (P2-IndexSumY)*CombP2Y/(IndexSumY+1); // the index is still the old one => expression for "nC k+1"
	  PowDistanceYC *= DeltaCellCenters[cell].y;      // Update PowDistanceYC
	}//endfor

	Cell_LHS_Inv_Value(iCell,jCell,kCell,cell-1,i-1) += CombP3Z*PowDistanceZC*IntSum2;  // update the external sum
	  
	CombP3Z = (P3-IndexSumZ)*CombP3Z/(IndexSumZ+1); // the index is still the old one => expression for "nC k+1"
	PowDistanceZC *= DeltaCellCenters[cell].z;      // Update PowDistanceYC
      }//endfor

      // subtract the corresponding geometric moment of cell (iCell,jCell) 
      Cell_LHS_Inv_Value(iCell,jCell,kCell,cell-1,i-1) -= Geom->CellGeomCoeffValue(iCell,jCell,kCell,P1,P2,P3);
      
      // apply geometric weighting
      Cell_LHS_Inv_Value(iCell,jCell,kCell,cell-1,i-1) *= GeomWeightValue(iCell,jCell,kCell,cell);

    } // endfor (i)
  }//endfor (cell)

  // STOP:   Matrix of the linear system (LHS) built. 
  //         For kExact_Reconstruction away from some special curved boundaries 
  //         the same matrix is used for all variables (same geometry) and  
  //         at every time step as long as the mesh is the same.
  // **********************************************************************

  // Compute the pseudo-inverse and override the LHS term.
  // This operation will change the dimensions of the matrix.
  Cell_LHS_Inv(iCell,jCell,kCell).pseudo_inverse_override();
  
}
// --> RR: comment out defn of ComputeCellReconstructionPseudoInverseNearConstrainedBoundaries
//
///*! 
// * Compute the pseudo-inverse of the left-hand-side term in the 
// * unlimited k-exact high-order reconstruction for a specified
// * computational cell based on the information provided by the
// * associated grid. 
// * The reconstruction of this cell is influence by the presence
// * of curved boundaries.
// * If the pseudo-inverse cannot be computed because boundary condition
// * constraints must be added to the linear system, the LHS matrix of the 
// * k-exact least-squares reconstruction is stored instead.
// *
// * \note This routine must be called ONLY IF the block is constrained (i.e. ReconstructionTypeMap has been created)
// */
//template<class SOLN_STATE>
//void HighOrder<SOLN_STATE>::ComputeCellReconstructionPseudoInverseNearConstrainedBoundaries(const int &iCell,
//											      const int &jCell){
//
//  ostringstream ErrorMsg;
//
//  // Set the biased stencil of points used for the reconstruction of the current cell
//  SetDeviatedReconstructionStencil(iCell, jCell,
//				   i_index_ave, j_index_ave,
//				   rings);
//    
//
//  if (ReconstructionTypeMap[iCell][jCell] == 'm'){
//    /* There are NO constraints for this cell,
//       so a pseudo-inverse with the biased stencil may be computed.
//       The reconstruction of this cell has only the stencil MODIFIED (i.e. deviated from central). */
//
//    // Check overdeterminancy
//    if ( i_index_ave.size() > NumberOfTaylorDerivatives() ){
//      // The pseudo-inverse CAN be computed
//
//      // Compute the pseudo-inverse for the current cell
//      ComputeCellReconstructionPseudoInverse(iCell, jCell, i_index_ave, j_index_ave);
//	      
//    } else {
//      // The pseudo-inverse CANNOT be computed
//      // Throw an error
//      // Set the error message
//      ErrorMsg << "HighOrder<SOLN_STATE>::ComputeCellReconstructionPseudoInverseNearConstrainedBoundaries() ERROR!"
//	       << " The pseudo-inverse couldn't be computed for cell (" 
//	       << iCell << "," << jCell << ")";
//      throw runtime_error(ErrorMsg.str());
//    }
//
//  } else if (ReconstructionTypeMap[iCell][jCell] == 'c'){
//    /* There are constraints for this cell, so the pseudo-inverse CANNOT be computed.
//       Instead, the least-squares part of the reconstruction matrix will be stored.
//       That is, the mean conservation equation in the reconstructed cell (first equation)
//       and the approximate equations for the neighbouring cells. */
//	    
//    /********* Generate the exact and approximate mean conservation equations ***********/
//    /************************************************************************************/
//    Set_LHS_MeanValueConservation_Equations(iCell,jCell,
//					    i_index_ave, j_index_ave,
//					    Cell_LHS_Inv(iCell,jCell),
//					    GeomWeights(iCell,jCell));
//  } else {
//    // The pseudo-inverse SHOULDN'T be computed with this routine
//    // Throw an error
//    // Set the error message
//    ErrorMsg << "HighOrder<SOLN_STATE>::ComputeCellReconstructionPseudoInverseNearConstrainedBoundaries() ERROR!"
//	     << "This routine should be used to compute pseudo-inverse only for 'modified' and 'constrained' reconstructions"
//	     << "The reconstruction type of cell (" << iCell << "," << jCell << ") is different than those!"; 
//    throw runtime_error(ErrorMsg.str());    
//  } // endif (ReconstructionTypeMap[iCell][jCell])  
//}

/*! 
 * Performs the reconstruction of a limited piecewise 
 * linear solution state within a given cell (iCell,jCell) of   
 * the computational mesh for the specified             
 * quadrilateral solution block.  A least squares       
 * approach is used in the evaluation of the unlimited  
 * solution gradients.  Several slope limiters may be   
 * used.
 * The high-order derivatives of the variables flagged as unfit
 * are replaced with the first-order ones computed with
 * this subroutine. 
 */
template<class SOLN_STATE>
template<class Soln_Block_Type> 
void HighOrder<SOLN_STATE>::
ComputeLimitedPiecewiseLinearSolutionReconstruction(Soln_Block_Type &SolnBlk,
						    const int &iCell, const int &jCell, const int &kCell,
						    const int &Limiter,
						    const Soln_State & 
						    (Soln_Block_Type::*ReconstructedSoln)(const int &,const int &,const int &) const){


  // --> RR: VIP: comment out (empty function) for ComputeLimitedPiecewiseLinearSolutionReconstruction
//
//  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS
//
//  int n, parameter, n_pts;
//  double MaxWeight(0.0);
//  double DxDx_ave(0), DxDy_ave(0), DyDy_ave(0);
//  Soln_State DU, DUDx_ave(0), DUDy_ave(0);
//  double u0Min, u0Max;
//  double *uQuad(NULL);
//  Vector3D *GQP(NULL);
//  int NumGQP, GQP_North, GQP_South, GQP_East, GQP_West;
//  bool faceNorth, faceSouth, faceEast, faceWest;
//
//  // == Check if the reconstruction polynomial is piecewise constant
//  if (RecOrder() == 0){
//    // There is no need to calculate the reconstruction
//    return;
//  }
//
//  /* Carry out the limited solution reconstruction in
//     each cell of the computational mesh. */
//
//  /* Determine the number of neighbouring cells to
//     be used in the reconstruction procedure.
//     This stencil might be different near boundaries
//     and it is influenced by the boundary conditions. */
//  SolnBlk.SetPiecewiseLinearReconstructionStencil(iCell,jCell,
//						  I_Index,J_Index,
//						  n_pts);
//
//  // Perform reconstruction.
//  if (n_pts > 0) {
//    
//    // Perform piecewise linear reconstruction only if the reconstruction order is greater than 1
//    if (RecOrder() > 1){
//      
//      // Compute distance between centroids and the geometric weights
//      for ( n = 0 ; n < n_pts ; ++n ) {
//	/* Compute the X and Y component of the distance between
//	   the cell centers of the neighbour and the reconstructed cell */
//	dX[n] = Geom->Cell[ I_Index[n] ][ J_Index[n] ].Xc - Geom->Cell[iCell][jCell].Xc;
//
//	/* Compute the geometric weight based on the centroid distance */
//	CENO_Geometric_Weighting(geom_weights[n], dX[n].abs());
//
//	/* Compute the maximum geometric weight (this is used for normalization) */
//	MaxWeight = max(MaxWeight, geom_weights[n]);
//      }
//
//      for ( n = 0 ; n < n_pts ; ++n ) {
//	// compute the normalized geometric weight
//	geom_weights[n] /= MaxWeight;
//
//	// compute the square of the normalized geometric weight
//	geom_weights[n] *= geom_weights[n];
//
//	DU = (SolnBlk.*ReconstructedSoln)( I_Index[n], J_Index[n] ) - (SolnBlk.*ReconstructedSoln)(iCell,jCell);
//	DUDx_ave += DU*(geom_weights[n]*dX[n].x);
//	DUDy_ave += DU*(geom_weights[n]*dX[n].y);
//	DxDx_ave += geom_weights[n]*dX[n].x*dX[n].x;
//	DxDy_ave += geom_weights[n]*dX[n].x*dX[n].y;
//	DyDy_ave += geom_weights[n]*dX[n].y*dX[n].y;
//      } /* endfor */
//    					    
//      DUDx_ave /= double(n_pts);
//      DUDy_ave /= double(n_pts);
//      DxDx_ave /= double(n_pts);
//      DxDy_ave /= double(n_pts);
//      DyDy_ave /= double(n_pts);
//
//      // Calculate the first-order derivatives
//      dUdx = ( (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
//	       (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
//      dUdy = ( (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
//	       (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
//
//    }//endif(RecOrder())
//
//
//    // Calculate slope limiters or used the frozen ones.
//    if (!SolnBlk.Freeze_Limiter) {
//
//      // Calculate the new slope limiters.
//
//      // Get number of flux calculation points and edge type for each of the 4 cell faces.
//      GQP_North = Geom->NumOfFluxCalculationGaussQuadPoints_North(iCell,jCell,faceNorth);
//      GQP_South = Geom->NumOfFluxCalculationGaussQuadPoints_South(iCell,jCell,faceSouth);
//      GQP_East  = Geom->NumOfFluxCalculationGaussQuadPoints_East(iCell,jCell,faceEast);
//      GQP_West  = Geom->NumOfFluxCalculationGaussQuadPoints_West(iCell,jCell,faceWest);
//
//      // Get total number of points which are used to assess the slope limiter
//      NumGQP = GQP_North + GQP_South + GQP_East + GQP_West;
//
//      // Allocate memory for the location of the points and the solution
//      uQuad = new double [NumGQP];
//      GQP = new Vector3D [NumGQP];
//      
//      // Get North face GQPs
//      n = 0;
//      if (faceNorth){
//
//	if (jCell == JCu){ // this is an interior cell
//	  // used GQPs from BndNorthSplineInfo
//	  Geom->BndNorthSplineInfo[iCell].CopyGQPoints(&GQP[n]);
//	} else { // this is a ghost cell
//	  // used GQPs from BndSouthSplineInfo
//	  Geom->BndSouthSplineInfo[iCell].CopyGQPoints(&GQP[n]);
//	}
//	
//      }	else {
//	// used GQPs from straight edge
//	Geom->getGaussQuadPointsFaceN(iCell,jCell,&GQP[n],GQP_North);
//      }
//      // Update n
//      n += GQP_North;
//      
//      // Get West face GQPs
//      if (faceWest){
//
//	if (iCell == ICl){ // this is an interior cell
//	  // used GQPs from BndWestSplineInfo
//	  Geom->BndWestSplineInfo[jCell].CopyGQPoints(&GQP[n]);
//	} else { // this is a ghost cell
//	  // used GQPs from BndEastSplineInfo
//	  Geom->BndEastSplineInfo[jCell].CopyGQPoints(&GQP[n]);
//	}
//	
//	
//      } else {
//	// used GQPs from straight edge
//	Geom->getGaussQuadPointsFaceW(iCell,jCell,&GQP[n],GQP_West);
//      }
//      // Update n
//      n += GQP_West;
//      
//      // Get South face GQPs
//      if (faceSouth){
//
//	if (jCell == JCl){ // this is an interior cell
//	  // used GQPs from BndSouthSplineInfo
//	  Geom->BndSouthSplineInfo[iCell].CopyGQPoints(&GQP[n]);
//	} else { // this is a ghost cell
//	  // used GQPs from BndNorthSplineInfo
//	  Geom->BndNorthSplineInfo[iCell].CopyGQPoints(&GQP[n]);
//	}
//
//      } else {
//	// used GQPs from straight edge
//	Geom->getGaussQuadPointsFaceS(iCell,jCell,&GQP[n],GQP_South);
//      }
//      // Update n
//      n += GQP_South;
//      
//      // Get East face GQPs
//      if (faceEast){
//
//	if (iCell == ICu){ // this is an interior cell
//	  // used GQPs from BndEastSplineInfo
//	  Geom->BndEastSplineInfo[jCell].CopyGQPoints(&GQP[n]);
//	} else { // this is a ghost cell
//	  // used GQPs from BndWestSplineInfo
//	  Geom->BndWestSplineInfo[jCell].CopyGQPoints(&GQP[n]);
//	}
//
//      } else {
//	// used GQPs from straight edge
//	Geom->getGaussQuadPointsFaceE(iCell,jCell,&GQP[n],GQP_East);
//      }
//    
//      // Calculate the limiter for each solution variable (i.e. parameter)
//      for (parameter = 1; parameter <= NumberOfVariables(); ++parameter) {
//      
//	// Drop the order only for the variables that are flagged as unfit
//	if ( CellInadequateFitValue(iCell,jCell,parameter) ){
//	
//	  if (RecOrder() > 1){
//	    // Zero all derivatives but D00 associated with this parameter.
//	    for (n = 1; n < NumberOfTaylorDerivatives(); ++n){
//	      CellTaylorDerivState(iCell,jCell,n)[parameter] = 0.0;
//	    }
//	  
//	    // Set D00 and U_ave(see memory pool) of the current parameter to the correspondent average solution.
//	    // U_ave is used in UnlimitedLinearSolutionAtLocation() (see below).
//	    U_ave[parameter] = CellTaylorDerivState(iCell,jCell,0)[parameter] = (SolnBlk.*ReconstructedSoln)(iCell,jCell)[parameter];
//	    
//	    // Set D01 and D10 to the values of the first-order derivatives.
//	    CellTaylorDerivValue(iCell,jCell,0,1,parameter) = dUdy[parameter];
//	    CellTaylorDerivValue(iCell,jCell,1,0,parameter) = dUdx[parameter];
//	    
//	  } else {
//	    // Set U_ave of the current parameter
//	    U_ave[parameter] = CellTaylorDerivState(iCell,jCell,0)[parameter];
//	    dUdy[parameter] = CellTaylorDerivState(iCell,jCell,1)[parameter];
//	    dUdx[parameter] = CellTaylorDerivState(iCell,jCell,2)[parameter];
//	  }
//
//	  // Compute the minimum and maximum average solution in the stencil for the current parameter
//	  u0Min = U_ave[parameter];
//	  u0Max = u0Min;
//	  for (n = 0; n < n_pts; ++n) {
//	    u0Min = min(u0Min, (SolnBlk.*ReconstructedSoln)(I_Index[n],J_Index[n])[parameter]);
//	    u0Max = max(u0Max, (SolnBlk.*ReconstructedSoln)(I_Index[n],J_Index[n])[parameter]);
//	  }// endfor(n)
//	  
//	  // Evaluate the solution at all required points
//	  for (n = 0; n < NumGQP; ++n){
//	    uQuad[n] = UnlimitedLinearSolutionAtLocation(iCell,jCell,
//							 GQP[n],
//							 parameter);
//	  }// endfor(n)
//
//	  // Evaluate limiter for the current parameter
//	  CellTaylorDeriv(iCell,jCell).Limiter(parameter) = CalculateLimiter(SolnBlk,uQuad,NumGQP,U_ave[parameter],
//									     u0Min,u0Max,Limiter);
//
//	  // Save a copy of the limiter for the current parameter for later usage if limiter freezing is required
//	  CellTaylorDeriv(iCell,jCell).Make_Limiter_Copy(parameter);
//	
//	}// endif
//
//      }//endfor(parameter)
//
//      // Deallocate memory
//      delete [] uQuad; uQuad = NULL;
//      delete [] GQP; GQP = NULL;
//      NumGQP = 0;
//
//    } else {
//      // Set limiters to the frozen ones.
//
//      for (parameter = 1; parameter <= NumberOfVariables(); ++parameter) {
//      
//	// Drop the order only for the variables that are flagged as unfit
//	if ( CellInadequateFitValue(iCell,jCell,parameter) ){
//	
//	  if (RecOrder() > 1){
//	    // Zero all derivatives but D00 associated with this parameter.
//	    for (n = 1; n < NumberOfTaylorDerivatives(); ++n){
//	      CellTaylorDerivState(iCell,jCell,n)[parameter] = 0.0;
//	    }
//	  
//	    // Set D00 and U_ave(see memory pool) of the current parameter to the correspondent average solution.
//	    // U_ave is used in UnlimitedLinearSolutionAtLocation() (see below).
//	    U_ave[parameter] = CellTaylorDerivState(iCell,jCell,0)[parameter] = (SolnBlk.*ReconstructedSoln)(iCell,jCell)[parameter];
//	    
//	    // Set D01 and D10 to the values of the first-order derivatives.
//	    CellTaylorDerivValue(iCell,jCell,0,1,parameter) = dUdy[parameter];
//	    CellTaylorDerivValue(iCell,jCell,1,0,parameter) = dUdx[parameter];
//	    
//	  } else {
//	    // Set U_ave of the current parameter
//	    U_ave[parameter] = CellTaylorDerivState(iCell,jCell,0)[parameter];
//	    dUdy[parameter] = CellTaylorDerivState(iCell,jCell,1)[parameter];
//	    dUdx[parameter] = CellTaylorDerivState(iCell,jCell,2)[parameter];
//	  }
//
//	  // Set limiter for the current parameter to the frozen one
//	  CellTaylorDeriv(iCell,jCell).Limiter(parameter) = CellTaylorDeriv(iCell,jCell).Frozen_Limiter(parameter);
//
//	}// endif
//
//      }//endfor(parameter)      
//
//    } // endif (!SolnBlk.Freeze_Limiter)
//
//  } else {
//    // Use piecewise constant
//
//    // Set D00 to average solution
//    CellTaylorDerivState(iCell,jCell,0) = (SolnBlk.*ReconstructedSoln)(iCell,jCell);
//    // Zero all derivatives
//    for (n = 1; n < NumberOfTaylorDerivatives(); ++n){
//      CellTaylorDerivState(iCell,jCell,n).Vacuum();
//    }
//    // Reset limiter
//    CellTaylorDeriv(iCell,jCell).Limiter().Vacuum();
//    
//  } /* endif */

}

/*! 
 * Compute the unlimited k-exact high-order reconstruction
 * proposed by Barth (1993) for the given cell.
 * The reconstruction stencil in these routine is not necessarily
 * the central one. Therefore, this routines is designed more for
 * testing purposes than for regular reconstructions.
 *
 * \param SolnBlk the quad block for which the solution reconstruction is done.
 * \param ReconstructedSoln member function of Soln_Block_Type which returns the solution.
 * \param iCell i-index of the reconstructed cell
 * \param jCell j-index of the reconstructed cell
 */
template<class SOLN_STATE>
template<class Soln_Block_Type> inline
void HighOrder<SOLN_STATE>::
ComputeUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk, 
				       const int &iCell, const int &jCell, const int & kCell,
				       const bool & UseSpecialStencil,
				       const Soln_State & 
				       (Soln_Block_Type::*ReconstructedSoln)(const int &,const int &,const int &) const){

  IndexType Special_I_Index(i_index.size()), Special_J_Index(j_index.size()), Special_K_Index(k_index.size());

  if ( iCell >= StartI && iCell <= EndI && jCell >= StartJ && jCell <= EndJ && kCell >= StartK && kCell <= EndK){

    /***************************************************************************
     *    Perform unconstrained unlimited high-order solution reconstruction   *
     **************************************************************************/

 //   if (UseSpecialStencil){
 //
 //     // Set the stencil of points used for reconstruction with the special routine
 //     SetSpecialReconstructionStencil(iCell, jCell, Special_I_Index, Special_J_Index);
 //
 //   } else {

    // Set the stencil of points used for reconstruction with the regular routine
    SetReconstructionStencil(iCell, jCell, kCell, Special_I_Index, Special_J_Index, Special_K_Index);

//    }

    // Compute the reconstruction for the current cell
    ComputeUnconstrainedUnlimitedSolutionReconstruction(SolnBlk, ReconstructedSoln,
							iCell, jCell, kCell, Special_I_Index, Special_J_Index, Special_K_Index);
    
  } else {
    
    /***************************************************************************
     *    Perform constrained unlimited high-order solution reconstruction     *
     **************************************************************************/

    // Add constrained reconstruction here

    throw runtime_error("HighOrder<SOLN_STATE>::ComputeUnlimitedSolutionReconstruction() doesn't know how to handle constrained reconstruction!");

  }

}

// --> RR: comment out of definitions of Constrained reconstructions
//
///*! 
// * Compute the unlimited k-exact high-order reconstruction
// * proposed by Barth (1993) combined with equations satisfied 
// * exactly (i.e. constraints) which account for boundary conditions.
// * This reconstruction procedure should be used for computational cells
// * affected by the presence of a boundary condition enforced with constrained 
// * reconstruction.
// * The mean quantity conservation is also enforced as a constraint.
// */
//template<class SOLN_STATE>
//template<class Soln_Block_Type>
//void HighOrder<SOLN_STATE>::
//ComputeConstrainedUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk,
//						  const Soln_State & 
//						  (Soln_Block_Type::*ReconstructedSoln)(const int &,const int &) const,
//						  const int &iCell, const int &jCell,
//						  const IndexType & i_index,
//						  const IndexType & j_index) {
//
//  // == Check if the reconstruction polynomial is piecewise constant
//  if (RecOrder() == 0){
//    // There is no need to calculate the reconstruction
//    return;
//  }
//
//  int constrGQP_W, constrGQP_E, constrGQP_N, constrGQP_S; // number of constrained points on each edge
//  bool IC_Flag(false), RC_Flag(false);			  // flags to indicate the type of constraints encountered
//
//  /**********************************************************************************
//   *  STEP 1. DETERMINE THE NUMBER AND TYPE OF CONSTRAINTS THAT NEED TO BE IMPOSED  *
//   *********************************************************************************/
//  
//  // Determine the number of constraints for the West boundary.
//  constrGQP_W = SolnBlk.Grid.NumOfConstrainedGaussQuadPoints_West(iCell,jCell);
//  if (constrGQP_W > 0){
//    // Ensure physical constraints (i.e. proper boundary conditions)
//    SolnBlk.EnsurePhysicalBCsConstraints(WEST,jCell);
//    // Identify what type of constraints are imposed by this face
//    RC_Flag = RC_Flag || SolnBlk.BC_WestCell(jCell).IsThereAnyRelationalConstraintRequired();
//    IC_Flag = IC_Flag || SolnBlk.BC_WestCell(jCell).IsThereAnyIndividualConstraintRequired();
//  }
//
//  // Determine the number of constraints for the South boundary
//  constrGQP_S = SolnBlk.Grid.NumOfConstrainedGaussQuadPoints_South(iCell,jCell);
//  if (constrGQP_S > 0){
//    // Ensure physical constraints (i.e. proper boundary conditions)
//    SolnBlk.EnsurePhysicalBCsConstraints(SOUTH,iCell);
//    // Identify what type of constraints are imposed by this face
//    RC_Flag = RC_Flag || SolnBlk.BC_SouthCell(iCell).IsThereAnyRelationalConstraintRequired();
//    IC_Flag = IC_Flag || SolnBlk.BC_SouthCell(iCell).IsThereAnyIndividualConstraintRequired();
//  }
//
//  // Determine the number of constraints for the East boundary
//  constrGQP_E = SolnBlk.Grid.NumOfConstrainedGaussQuadPoints_East(iCell,jCell);
//  if (constrGQP_E > 0){
//    // Ensure physical constraints (i.e. proper boundary conditions)
//    SolnBlk.EnsurePhysicalBCsConstraints(EAST,jCell);
//    // Identify what type of constraints are imposed by this face
//    RC_Flag = RC_Flag || SolnBlk.BC_EastCell(jCell).IsThereAnyRelationalConstraintRequired();
//    IC_Flag = IC_Flag || SolnBlk.BC_EastCell(jCell).IsThereAnyIndividualConstraintRequired();
//  }
//    
//  // Determine the number of constraints for the North boundary
//  constrGQP_N = SolnBlk.Grid.NumOfConstrainedGaussQuadPoints_North(iCell,jCell);
//  if (constrGQP_N > 0){
//    // Ensure physical constraints (i.e. proper boundary conditions)
//    SolnBlk.EnsurePhysicalBCsConstraints(NORTH,iCell);
//    // Identify what type of constraints are imposed by this face
//    RC_Flag = RC_Flag || SolnBlk.BC_NorthCell(iCell).IsThereAnyRelationalConstraintRequired();
//    IC_Flag = IC_Flag || SolnBlk.BC_NorthCell(iCell).IsThereAnyIndividualConstraintRequired();
//  }
//
//
//  /***********************************************************************************************
//   * STEP 2. PROCEED WITH A DIFFERENT ALGORITHM DEPENDING ON THE TYPE OF ENCOUNTERED CONSTRAINTS *
//   **********************************************************************************************/
//  
//  if (RC_Flag && IC_Flag){
//    // There are both relational and individual constraints
//    ComputeRelationallyAndIndividuallyConstrainedUnlimitedSolutionReconstruction(SolnBlk,
//										 ReconstructedSoln,
//										 iCell, jCell,
//										 i_index, j_index,
//										 constrGQP_W, constrGQP_S,
//										 constrGQP_E, constrGQP_N);
//  } else if (RC_Flag){
//    // There are only relational constraints
//    ComputeRelationallyConstrainedUnlimitedSolutionReconstruction(SolnBlk,
//								  ReconstructedSoln,
//								  iCell, jCell,
//								  i_index, j_index,
//								  constrGQP_W, constrGQP_S,
//								  constrGQP_E, constrGQP_N);
//  } else if (IC_Flag){
//    // There are only individual constraints
//    ComputeIndividuallyConstrainedUnlimitedSolutionReconstruction(SolnBlk,
//								  ReconstructedSoln,
//								  iCell, jCell,
//								  i_index, j_index,
//								  constrGQP_W, constrGQP_S,
//								  constrGQP_E, constrGQP_N);
//  } else {
//    // There are no constraints
//    ComputeUnconstrainedUnlimitedSolutionReconstruction(SolnBlk,
//							ReconstructedSoln,
//							iCell, jCell,
//							i_index, j_index);
//  } // endif
//
//}


/*!
 * Generate the LHS and RHS of the least-squares problem associated 
 * with the reconstruction procedure and write the values at the specified
 * locations.
 * This matrix reflects the conservation of mean value in the control       
 * volumes of cells specified by (i_index,j_index,k_index).                         
 * The row associated with cell (iCell,jCell,kCell) represents a constraint and    
 * therefore it must be satisfied exactly.
 * The rest of the system is approximately solved in the least squares sense.
 *                                                                          
 * \note It is assumed that the (i_index[0],j_index[0],k_index[0]) is equal to          
 *       (iCell,jCell,kCell). That is, the first line of A is a constraint!!!            
 * The constraint is filled in the matrix A at the position (RowConstraint, StartCol).                                     
 * The approximate equations are filled in the matrix starting from (StartRow, StartCol) position. 
 *
 * \param SolnBlk the quad block for which the solution reconstruction is done.
 * \param ReconstructedSoln member function of Soln_Block_Type which returns the solution.
 * \param iCell i-index of the reconstructed cell
 * \param jCell j-index of the reconstructed cell
 * \param kCell k-index of the reconstructed cell
 * \param ParameterIndex related to the indexes of the solution
 * \param A the LHS assemble matrix 
 * \param B the RHS assemble matrix
 *
 * \note This routine doesn't normalize the geometric weights.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type> inline
void HighOrder<SOLN_STATE>::
Set_MeanValueConservation_Equations(Soln_Block_Type & SolnBlk,
				    const Soln_State & 
				    (Soln_Block_Type::*ReconstructedSoln)(const int &,const int &,const int &) const,
				    const int &iCell, const int &jCell, const int & kCell,
				    const IndexType & i_index, const IndexType & j_index, const IndexType & k_index,
				    DenseMatrix & A, DenseMatrix & All_U,
				    const IndexType & ParameterIndex,
				    const int &RowConstraint,
				    const int &StartRow, const int &StartCol){


  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS

  int StencilSize(i_index.size());
  int ParameterIndexSize(ParameterIndex.size());
  ColumnVector GeomWeights(StencilSize);   // The column vector of the geometric weights
  Vector3D *DeltaCellCenters;              /* stores the difference between the cell center of
					      neighbour cells and the one of (i,j,k) cell.
					      The first value is 0!!! */
  int IndexSumZ,IndexSumY, IndexSumX, P1, P2,P3;
  double CombP1X, CombP2Y, CombP3Z;
  double PowDistanceZC, PowDistanceYC, PowDistanceXC;
  int cell, i, parameter;
  double IntSum1(0.0), IntSum2(0.0);

  // Allocate memory
  DeltaCellCenters = new Vector3D [StencilSize];

  // START:   Set the bottom part of the LHS and RHS of the linear system 
  // *********************************************************************

  // Step1. Set the constraint equation
  for (i=0; i <= CellTaylorDeriv(iCell,jCell,kCell).LastElem(); ++i){
    A(RowConstraint,i+StartCol) = Geom->CellGeomCoeffValue(iCell,jCell,kCell,i);
  }

  for (parameter=0; parameter<ParameterIndexSize; ++parameter){
    All_U(RowConstraint,parameter) = (SolnBlk.*ReconstructedSoln)(iCell,jCell,kCell)[ParameterIndex[parameter]];
  }

  // Step3. Set the approximate equations
  for (cell=1; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil

    /* Compute the X, Y, and Z component of the distance between
       the cell center of the neighbours and the reconstructed cell */
    DeltaCellCenters[cell] = CellCenter(i_index[cell],j_index[cell],k_index[cell]) - CellCenter(iCell,jCell,kCell);

    /* Compute the geometric weight based on the centroid distance */
    CENO_Geometric_Weighting(GeomWeights[cell], DeltaCellCenters[cell].abs());

    // *** SET the matrix A of the linear system (LHS) ***
    /* compute for each derivative the corresponding entry in the matrix of the linear system */
    for (i=0; i<=CellTaylorDeriv(iCell,jCell,kCell).LastElem(); ++i){
      // build the row of the matrix
      P1 = CellTaylorDeriv(iCell,jCell,kCell,i).P1();  // identify P1
      P2 = CellTaylorDeriv(iCell,jCell,kCell,i).P2();  // identify P2
      P3 = CellTaylorDeriv(iCell,jCell,kCell,i).P3();  // identify P3

      A(cell+StartRow-1,i+StartCol) = 0.0;  // set sumation variable to zero
      CombP3Z = 1.0;        // the binomial coefficient "nC k" for k=0 is 1
      PowDistanceZC = 1.0;  // initialize PowDistanceZC

      // Compute geometric integral over the neighbour's domain
      for (IndexSumZ = 0; IndexSumZ<=P3; ++IndexSumZ){
	CombP2Y = 1.0;        // the binomial coefficient "nC k" for k=0 is 1
	PowDistanceYC = 1.0;  // initialize PowDistanceYC
	IntSum2 = 0.0;         // reset internal summation variable

	for (IndexSumY = 0; IndexSumY<=P2; ++IndexSumY){
	  CombP1X = 1.0;       // the binomial coefficient "nC k" for k=0 is 1
	  PowDistanceXC = 1.0; // initialize PowDistanceXC
	  IntSum1 = 0.0;	     // reset internal sumation variable

	  for (IndexSumX = 0; IndexSumX<=P1; ++IndexSumX){
	    IntSum1 += ( CombP1X*PowDistanceXC*
			 Geom->CellGeomCoeffValue(i_index[cell],j_index[cell],k_index[cell],P1-IndexSumX,P2-IndexSumY,P3-IndexSumZ) );
	      
	    // update the binomial coefficients
	    CombP1X = (P1-IndexSumX)*CombP1X/(IndexSumX+1); // the index is still the old one => expression for "nC k+1"
	    PowDistanceXC *= DeltaCellCenters[cell].x;      // Update PowDistanceXC
	  }//endfor
	    
	  IntSum2 += CombP2Y*PowDistanceYC*IntSum1;
	  CombP2Y = (P2-IndexSumY)*CombP2Y/(IndexSumY+1); // the index is still the old one => expression for "nC k+1"
	  PowDistanceYC *= DeltaCellCenters[cell].y;      // Update PowDistanceYC
	}//endfor
	
	A(cell-1,i-1) += CombP3Z*PowDistanceZC*IntSum2;  // update the external sum
	
	CombP3Z = (P3-IndexSumZ)*CombP3Z/(IndexSumZ+1); // the index is still the old one => expression for "nC k+1"
	PowDistanceZC *= DeltaCellCenters[cell].z;      // Update PowDistanceYC
      }//endfor

      // subtract the corresponding geometric moment of cell (iCell,jCell) 
      //A(cell-1,i-1) -= Geom->CellGeomCoeffValue(iCell,jCell,kCell,P1,P2,P3);

      // apply geometric weighting
      A(cell+StartRow-1,i+StartCol) *= GeomWeights(cell);

    } //endfor (i)

    // --> RR: from 2D implementation: delete later
//      // Compute geometric integral over the neighbour's domain
//      for (IndexSumY = 0; IndexSumY<=P2; ++IndexSumY){
//	CombP1X = 1.0;       // the binomial coefficient "nC k" for k=0 is 1
//	PowDistanceXC = 1.0; // initialize PowDistanceXC
//	IntSum = 0.0;	     // reset internal sumation variable
//
//	for (IndexSumX = 0; IndexSumX<=P1; ++IndexSumX){
//	  IntSum += ( CombP1X*PowDistanceXC*
//		      Geom->CellGeomCoeffValue(i_index[cell],j_index[cell],P1-IndexSumX,P2-IndexSumY) );
//	    
//	  // update the binomial coefficients
//	  CombP1X = (P1-IndexSumX)*CombP1X/(IndexSumX+1);  // The index is still the old one => expression for "nC k+1"
//	  PowDistanceXC *= DeltaCellCenters[cell].x;       // Update PowDistanceXC
//	}//endfor (IndexSumX)
//
//	A(cell+StartRow-1,i+StartCol) += CombP2Y*PowDistanceYC*IntSum; // update the external sum
//
//	CombP2Y = (P2-IndexSumY)*CombP2Y/(IndexSumY+1); // the index is still the old one => expression for "nC k+1"
//	PowDistanceYC *= DeltaCellCenters[cell].y;      // Update PowDistanceYC
//      }//endfor (IndexSumY)
//
//      // apply geometric weighting
//      A(cell+StartRow-1,i+StartCol) *= GeomWeights(cell);
//
//    }//endfor (i)
      
    // *** SET the matrix All_U of the linear system (RHS) ***
    for (parameter=0; parameter < ParameterIndexSize; ++parameter){
      All_U(cell+StartRow-1,parameter) = ( GeomWeights(cell)*
					   (SolnBlk.*ReconstructedSoln)(i_index[cell],j_index[cell],k_index[cell])[ParameterIndex[parameter]] );
    }
  } //endfor (cell)

  delete [] DeltaCellCenters;

}

/*!
 * Generate the LHS of the least-squares problem associated with the 
 * reconstruction procedure and write the values at the specified locations.
 * This matrix reflects the conservation of mean value in the control
 * volumes of cells specified by (i_index,j_index,k_index).
 * The row associated with cell (iCell,jCell,kCell) represents a constraint and    
 * therefore it must be satisfied exactly.
 * The rest of the system is approximately solved in the least squares sense.
 *                                                                          
 * \note It is assumed that the (i_index[0],j_index[0],k_index[0]) is equal to          
 *       (iCell,jCell,kCell). That is, the first line of A is a constraint!!!
 * The approximate equations are filled in the rest of the matrix.
 *
 * \param iCell i-index of the reconstructed cell
 * \param jCell j-index of the reconstructed cell
 * \param kCell k-index of the reconstructed cell
 * \param A the LHS matrix
 * \param CellGeometricWeights the array of geometric weights
 *
 */
template<class SOLN_STATE> inline
void HighOrder<SOLN_STATE>::
Set_LHS_MeanValueConservation_Equations(const int &iCell, const int &jCell, const int &kCell,
					const IndexType & i_index, const IndexType & j_index, const IndexType & k_index,
					DenseMatrix & A,
					DoubleArrayType & GeometricWeights){


  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS

  int StencilSize(i_index.size());
  Vector3D *DeltaCellCenters;              /* stores the difference between the cell center of
					      neighbour cells and the one of (i,j) cell.
					      The first value is 0!!! */
  int IndexSumZ, IndexSumY, IndexSumX, P1, P2, P3;
  double CombP1X, CombP2Y, CombP3Z;
  double PowDistanceZC, PowDistanceYC, PowDistanceXC;
  int cell, i;
  double IntSum1(0.0), IntSum2(0.0);
  double MaxWeight(0.0);

  // Allocate memory
  DeltaCellCenters = new Vector3D [StencilSize];

  // Ensure that enough memory is allocated for the current least-squares problem.
  if (A.size(0) != StencilSize || A.size(1) != NumberOfTaylorDerivatives()) {
    // Resize the matrix accordingly
    A.newsize(StencilSize,NumberOfTaylorDerivatives());
  }
  if (GeometricWeights.size() != StencilSize){
    // Resize the vector accordingly
    GeometricWeights.resize(StencilSize);
  }

  // START:   Set the LHS of the linear system 
  // *****************************************************

  // ==== Set the geometric weight associated with the reconstructed cell
  GeometricWeights[0] = 1;

  // Step1. Compute the normalized geometric weights
  for (cell=1; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil
    
    /* Compute the X, Y, and Z component of the distance between
       the cell centers of the neighbour and the reconstructed cell */
    DeltaCellCenters[cell] = CellCenter(i_index[cell],j_index[cell],k_index[cell]) - CellCenter(iCell,jCell,kCell);
    
    /* Compute the geometric weight based on the centroid distance */
    CENO_Geometric_Weighting(GeometricWeights[cell], DeltaCellCenters[cell].abs());
    
    /* Compute the maximum geometric weight (this is used for normalization) */
    MaxWeight = max(MaxWeight, GeometricWeights[cell]);
  }

  // Step 2. Set the constraint equation
  for (i=0; i <= CellTaylorDeriv(iCell,jCell,kCell).LastElem(); ++i){
    A(0,i) = Geom->CellGeomCoeffValue(iCell,jCell,kCell,i);
  }

  // Step 3. Set the approximate equations
  for (cell=1; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil

    // compute the normalized geometric weight
    GeometricWeights[cell] /= MaxWeight;

    // *** SET the matrix A of the linear system (LHS) ***
    /* compute for each derivative the corresponding entry in the matrix of the linear system */
    for (i=0; i<=CellTaylorDeriv(iCell,jCell,kCell).LastElem(); ++i){
      // build the row of the matrix
      P1 = CellTaylorDeriv(iCell,jCell,i).P1();  // identify P1
      P2 = CellTaylorDeriv(iCell,jCell,i).P2();  // identify P2
      P3 = CellTaylorDeriv(iCell,jCell,kCell,i).P3();  // identify P3
      
      A(cell,i) = 0.0;  // set sumation variable to zero
      CombP3Z = 1.0;        // the binomial coefficient "nC k" for k=0 is 1
      PowDistanceZC = 1.0;  // initialize PowDistanceZC

      // Compute geometric integral over the neighbour's domain
      for (IndexSumZ = 0; IndexSumZ<=P3; ++IndexSumZ){
	CombP2Y = 1.0;        // the binomial coefficient "nC k" for k=0 is 1
	PowDistanceYC = 1.0;  // initialize PowDistanceYC
	IntSum2 = 0.0;         // reset internal summation variable
	
	for (IndexSumY = 0; IndexSumY<=P2; ++IndexSumY){
	  CombP1X = 1.0;       // the binomial coefficient "nC k" for k=0 is 1
	  PowDistanceXC = 1.0; // initialize PowDistanceXC
	  IntSum1 = 0.0;	     // reset internal sumation variable
	  
	  for (IndexSumX = 0; IndexSumX<=P1; ++IndexSumX){
	    IntSum1 += ( CombP1X*PowDistanceXC*
			 Geom->CellGeomCoeffValue(i_index[cell],j_index[cell],k_index[cell],P1-IndexSumX,P2-IndexSumY,P3-IndexSumZ) );
	    
	    // update the binomial coefficients
	    CombP1X = (P1-IndexSumX)*CombP1X/(IndexSumX+1); // the index is still the old one => expression for "nC k+1"
	    PowDistanceXC *= DeltaCellCenters[cell].x;      // Update PowDistanceXC
	  }//endfor
	  
	  IntSum2 += CombP2Y*PowDistanceYC*IntSum1;
	  CombP2Y = (P2-IndexSumY)*CombP2Y/(IndexSumY+1); // the index is still the old one => expression for "nC k+1"
	  PowDistanceYC *= DeltaCellCenters[cell].y;      // Update PowDistanceYC
	}//endfor
	
	A(cell,i) += CombP3Z*PowDistanceZC*IntSum2;  // update the external sum
	
	CombP3Z = (P3-IndexSumZ)*CombP3Z/(IndexSumZ+1); // the index is still the old one => expression for "nC k+1"
	PowDistanceZC *= DeltaCellCenters[cell].z;      // Update PowDistanceYC
      }//endfor
      
      // apply geometric weighting
      A(cell,i) *= GeometricWeights[cell];
    
    } //endfor (i)

    // --> RR: from 2D implementation: delete later
//      // Compute geometric integral over the neighbour's domain
//      for (IndexSumY = 0; IndexSumY<=P2; ++IndexSumY){
//	CombP1X = 1.0;       // the binomial coefficient "nC k" for k=0 is 1
//	PowDistanceXC = 1.0; // initialize PowDistanceXC
//	IntSum = 0.0;	     // reset internal sumation variable
//
//	for (IndexSumX = 0; IndexSumX<=P1; ++IndexSumX){
//	  IntSum += ( CombP1X*PowDistanceXC*
//		      Geom->CellGeomCoeffValue(i_index[cell],j_index[cell],P1-IndexSumX,P2-IndexSumY) );
//	    
//	  // update the binomial coefficients
//	  CombP1X = (P1-IndexSumX)*CombP1X/(IndexSumX+1);  // The index is still the old one => expression for "nC k+1"
//	  PowDistanceXC *= DeltaCellCenters[cell].x;       // Update PowDistanceXC
//	}//endfor (IndexSumX)
//
//	A(cell,i) += CombP2Y*PowDistanceYC*IntSum; // update the external sum
//
//	CombP2Y = (P2-IndexSumY)*CombP2Y/(IndexSumY+1); // the index is still the old one => expression for "nC k+1"
//	PowDistanceYC *= DeltaCellCenters[cell].y;      // Update PowDistanceYC
//      }//endfor (IndexSumY)
//
//      // apply geometric weighting
//      A(cell,i) *= GeometricWeights[cell];
//
//    }//endfor (i)
      
  } //endfor (cell)


  delete [] DeltaCellCenters;

}

// --> RR: huge comment out of def'ns of contrained reconstructions 6 fcns
//
///*!
// * Set the constraints equations in the constrained least-square reconstruction.
// * These constraints are called individual because they are independent of
// * other solution parameters
// * (i.e they don't express any relationship between solution parameters).
// * The starting position for the entries in the LHS and RHS of the linear
// * system are specified by the (StartRow, StartCol) input parameters.   \n                                
// *                                                                          
// * The BCs handled by this subroutine have the mixed form:                  
// *   F(GQP, A_GQP, B_GQP) = A_GQP*F1(GQP) + B_GQP*F2(GQP)                   
// *                                                                          
// * where: GQP      -- Gauss Quadrature Point Locations (i.e. flux calculation points)
// *        A_GQP    -- The coefficient for the Dirichlet boundary condition at GQP                                                
// *        B_GQP    -- The coefficient for the Neumann boundary condition at GQP                                                
// *        F1(GQP)  -- The value of the Dirichlet boundary condition at GQP  
// *        F2(GQP)  -- The value of the Neumann boundary condition at GQP    
// *
// * \param SolnBlk the quad block for which the solution reconstruction is done.
// * \param iCell i-index of the reconstructed cell
// * \param jCell j-index of the reconstructed cell
// * \param Constraints_Loc GQP array
// * \param Constraints_Normals normal vectors at GQP locations
// * \param Constraints_BCs provide the boundary condition coefficients (i.e. A_GQP, B_GQP, F1, F2)
// * \param ParameterIndex related to the indexes of the solution
// * \param A the LHS assemble matrix 
// * \param B the RHS assemble matrix
// *
// * \note This routine uses the boundary condition coefficients (i.e. 'a' and 'b')
// *       from the first entry in the ParameterIndex!
// */
//template<class SOLN_STATE>
//template<class Soln_Block_Type> inline
//void HighOrder<SOLN_STATE>::
//Generalized_IndividualConstraints_Equations(Soln_Block_Type & SolnBlk,
//					    const int &iCell, const int &jCell,
//					    Vector3DArray & Constraints_Loc,
//					    Vector3DArray & Constraints_Normals,
//					    BC_Type_Array & Constraints_BCs,
//					    DenseMatrix & A, DenseMatrix & All_U,
//					    const IndexType & ParameterIndex,
//					    const int &StartRow, const int &StartCol) {
//  
//  int P1, P2, i;
//  int BCs, BCs_Entry, Eq;			// equation
//  double PowXC, PowYC;		/* PowXC = DistXi^(P1-1); PowYC = DistYi^(P2-1) */
//  double DistXi, DistYi;
//  int IndexP1, IndexP2;
//  double GeometricWeight;
//
//  
//  for (BCs_Entry = 0, Eq = 0; BCs_Entry < Constraints_BCs.size(); ++BCs_Entry){ // for each boundary condition entry
//
//    for (BCs = 1; BCs <= Constraints_BCs[BCs_Entry]->NumOfPoints(); ++BCs, ++Eq){ // for each flux calculation point
//
//      // Compute entrie in the LHS and RHS for the current equation
//      
//      // Determine distance between the current GQP and the centroid of cell (iCell,jCell)
//      DistXi = Constraints_Loc[Eq].x - XCellCenter(iCell,jCell);
//      DistYi = Constraints_Loc[Eq].y - YCellCenter(iCell,jCell);
//
//      // Step 1. Form the LHS  -- build the row of the matrix A associated with the current GQP
//      for (i=0; i<=CellTaylorDeriv(iCell,jCell).LastElem(); ++i){
//	// build the row of the matrix
//	P1 = CellTaylorDeriv(iCell,jCell,i).P1();  // identify P1
//	P2 = CellTaylorDeriv(iCell,jCell,i).P2();  // identify P2
//
//	/* Initialize PowXC & PowYC */
//	PowXC = 1.0/DistXi;
//	PowYC = 1.0/DistYi;
//
//	/* Update PowXC & PowYC */
//	for (IndexP1 = 1; IndexP1 <= P1; ++IndexP1){ PowXC *= DistXi; }
//	for (IndexP2 = 1; IndexP2 <= P2; ++IndexP2){ PowYC *= DistYi; }
//
//	// Use the first entry in the ParameterIndex to retrieve data for the a(BCs) and b(BCs).
//	A(Eq+StartRow,i+StartCol) = ( PowXC * PowYC * 
//				      ( Constraints_BCs[BCs_Entry]->a(BCs)[ParameterIndex[0]] * DistXi * DistYi + 
//					Constraints_BCs[BCs_Entry]->b(BCs)[ParameterIndex[0]] * (P1 * DistYi * 
//												 Constraints_Normals[Eq].x  + 
//												 P2 * DistXi * 
//												 Constraints_Normals[Eq].y )) );
//
//      } //endfor (i)
//       
//
//      // Step 2. Form the RHS  -- build the row of the matrix All_U associated with the current GQP
//      for (i=0; i<ParameterIndex.size(); ++i){
//	All_U(Eq+StartRow, i) = ( ( Constraints_BCs[BCs_Entry]->a(BCs)[ParameterIndex[0]] * 
//				    Constraints_BCs[BCs_Entry]->DirichletBC(BCs)[ParameterIndex[i]] ) + 
//				  ( Constraints_BCs[BCs_Entry]->b(BCs)[ParameterIndex[0]] * 
//				    Constraints_BCs[BCs_Entry]->NeumannBC(BCs)[ParameterIndex[i]] ) );
//
//      }//endfor (i)
//
//    } // endfor (BCs)
//
//  }// endfor (BCs_Entry) 
//
//}
//
//
///*!
// * Set the constraints equations in the constrained least-square reconstruction.
// * These constraints are called relational because they put in relationship
// * several solution parameters.
// * The starting position for the entries in the LHS and RHS of the linear
// * system are specified by the (StartRow, StartCol) input parameters.   \n                                
// *                                                                          
// * The BCs handled by this subroutine are specific to each solver,
// * so this subroutine MUST be specialized!
// *
// * The data passed:
// *  GQP      -- Gauss Quadrature Point Locations (i.e. flux calculation points)
// *  A_GQP    -- The coefficient for the Dirichlet boundary condition at GQP                                                
// *  B_GQP    -- The coefficient for the Neumann boundary condition at GQP                                                
// *  F1(GQP)  -- The value of the Dirichlet boundary condition at GQP  
// *  F2(GQP)  -- The value of the Neumann boundary condition at GQP    
// *
// * \param SolnBlk the quad block for which the solution reconstruction is done.
// * \param iCell i-index of the reconstructed cell
// * \param jCell j-index of the reconstructed cell
// * \param Constraints_Loc GQP array
// * \param Constraints_Normals normal vectors at GQP locations
// * \param Constraints_BCs provide the boundary condition coefficients (i.e. A_GQP, B_GQP, F1, F2)
// * \param BC_Type the type of the boundary condition for which the relational constraint is built
// * \param ParameterIndex related to the indexes of the solution
// * \param A the LHS assemble matrix 
// * \param B the RHS assemble matrix
// *
// */
//template<class SOLN_STATE>
//template<class Soln_Block_Type> inline
//void HighOrder<SOLN_STATE>::
//Generalized_RelationalConstraints_Equations(Soln_Block_Type & SolnBlk,
//					    const int &iCell, const int &jCell,
//					    Vector3DArray & Constraints_Loc,
//					    Vector3DArray & Constraints_Normals,
//					    BC_Type_Array & Constraints_BCs,
//					    const int & BC_Type,
//					    DenseMatrix & A, DenseMatrix & All_U,
//					    const IndexType & ParameterIndex,
//					    const int &StartRow, const int &StartCol) {  
//  // Specialize this routine!
//}
//
//
///*! 
// * Compute the unlimited k-exact high-order reconstruction
// * proposed by Barth (1993) combined with exactly satisfied equations 
// * (i.e. only relational constraints) which account for boundary conditions.
// * This reconstruction procedure should be used for computational cells
// * affected by the presence of a boundary condition enforced with relational constraints.
// * Relational constraints represent exact equations relating multiple solution parameters.
// * The mean quantity conservation is also enforced as a constraint for all solution parameters.
// *
// * \param SolnBlk the quad block for which the solution reconstruction is done.
// * \param ReconstructedSoln member function of SolnBlk which returns the average solution of (iCell,jCell)
// * \param iCell i-index of the reconstructed cell
// * \param jCell j-index of the reconstructed cell
// * \param i_index i-indexes of the cells that are part of the supporting stencil
// * \param j_index j-indexes of the cells that are part of the supporting stencil
// * \param ConstrainedGQPs_West number of Gauss quadrature points constrained on West cell edge
// * \param ConstrainedGQPs_South number of Gauss quadrature points constrained on South cell edge
// * \param ConstrainedGQPs_East number of Gauss quadrature points constrained on East cell edge
// * \param ConstrainedGQPs_North number of Gauss quadrature points constrained on North cell edge
// */
//template<class SOLN_STATE>
//template<class Soln_Block_Type> inline
//void HighOrder<SOLN_STATE>::
//ComputeRelationallyConstrainedUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk,
//							      const Soln_State & 
//							      (Soln_Block_Type::*ReconstructedSoln)(const int &,
//												    const int &) const,
//							      const int &iCell, const int &jCell,
//							      const IndexType & i_index,
//							      const IndexType & j_index,
//							      const int & ConstrainedGQPs_West,
//							      const int & ConstrainedGQPs_South,
//							      const int & ConstrainedGQPs_East,
//							      const int & ConstrainedGQPs_North) {
//
//  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS
//  IndexType ParametersWithPhysicalBCs,   //!< List of solution parameters that have physical BCs (i.e. have relational constraints)
//    ParametersWithNumericalBCs,		 //!< List of solution parameters that have numerical BCs (i.e. no constraints)
//    ParameterIndex(1);			 //!< List of constrained solution parameters (used for current calculations)
//  //! List for West edge of BCtype for each parameter (i.e. how many conditions per edge)
//  int *RelationalBCs_W = new int [NumberOfVariables()+1];
//  //! List for South edge of BCtype for each parameter (i.e. how many conditions per edge)
//  int *RelationalBCs_S = new int [NumberOfVariables()+1];
//  //! List for East edge of BCtype for each parameter (i.e. how many conditions per edge)
//  int *RelationalBCs_E = new int [NumberOfVariables()+1];
//  //! List for North edge of BCtype for each parameter (i.e. how many conditions per edge)
//  int *RelationalBCs_N = new int [NumberOfVariables()+1];
//  int parameter;
//  int TotalNumberOfEquations, TotalNumberOfExactlySatisfiedEquations, TotalNumberOfRelationalBCsEquations(0);
//  int ApproxMeanConservationRow, MeanConservationRow, MeanConservationCol(0);
//  int cell, n, iterator;
//  int BC_Type, TotalRelational;
//  int StencilSize(i_index.size());
//  DenseMatrix X;		//< the matrix of unknowns
//
//  
//  /********************************************************************************************************
//   *  STEP 1. SORT THE PARAMETERS INTO THE DESIGNATED CATEGORIES AND DETECT CONSTRAINTS ON EACH CELL EDGE *
//   *******************************************************************************************************/
//  for (parameter = 1; parameter <= NumberOfVariables(); ++parameter){
//
//    // (Re)set variables
//    TotalRelational = 0;
//
//    RelationalBCs_W[parameter] = 0;
//    RelationalBCs_S[parameter] = 0;
//    RelationalBCs_E[parameter] = 0;
//    RelationalBCs_N[parameter] = 0;
//
//    // West edge
//    if (ConstrainedGQPs_West > 0){
//      RelationalBCs_W[parameter] = ConstrainedGQPs_West * SolnBlk.BC_WestCell(jCell).NumberOfRelationalConstraints(parameter);
//      TotalRelational += RelationalBCs_W[parameter];
//    }
//
//    // South edge
//    if (ConstrainedGQPs_South > 0){
//      RelationalBCs_S[parameter] = ConstrainedGQPs_South * SolnBlk.BC_SouthCell(iCell).NumberOfRelationalConstraints(parameter);
//      TotalRelational += RelationalBCs_S[parameter];
//    }
//
//    // East edge
//    if (ConstrainedGQPs_East > 0){
//      RelationalBCs_E[parameter] = ConstrainedGQPs_East * SolnBlk.BC_EastCell(jCell).NumberOfRelationalConstraints(parameter);
//      TotalRelational += RelationalBCs_E[parameter];
//    }
//
//    // North edge
//    if (ConstrainedGQPs_North > 0){
//      RelationalBCs_N[parameter] = ConstrainedGQPs_North * SolnBlk.BC_NorthCell(iCell).NumberOfRelationalConstraints(parameter);
//      TotalRelational += RelationalBCs_N[parameter];
//    }
//    
//    if (TotalRelational){
//      // At least one edge has relational constraints for this parameter
//      ParametersWithPhysicalBCs.push_back(parameter);
//    } else {
//      // There are no constraints for this parameter.
//      ParametersWithNumericalBCs.push_back(parameter);
//    }
//  }// endfor
//
//
//  /****************************************************************************
//   *  STEP 2. SOLVE THE RECONSTRUCTION FOR PARAMETERS THAT HAVE NUMERICAL BCs *
//   ***************************************************************************/
//  ComputeUnconstrainedUnlimitedSolutionReconstructionInConstrainedCell(SolnBlk, ReconstructedSoln,
//								       iCell, jCell,
//								       i_index, j_index,
//								       ParametersWithNumericalBCs);
//  
//
//  /****************************************************************************************
//   *  STEP 3. SOLVE THE RECONSTRUCTION FOR PARAMETERS THAT HAVE PHYSICAL (RELATIONAL) BCs *
//   ***************************************************************************************/
//  /* Note: It is assumed that each GQP introduces one relational constraint 
//     that relates all parameters in the ParametersWithPhysicalBCs vector. */
//
//  if (ParametersWithPhysicalBCs.size() == 0){
//    throw runtime_error("HighOrder<SOLN_STATE>::ComputeRelationallyConstrainedUnlimitedSolutionReconstruction() ERROR! Reconstruction with relational constraints is required but there are no relational constraints present!!");
//  }
//
//  // === Reset and initialize the local variables ===
//  // Calculate TotalNumberOfExactlySatisfiedEquations
//  // Initialize with the number of mean conservation constraints
//  TotalNumberOfExactlySatisfiedEquations = ParametersWithPhysicalBCs.size();
//
//  // Add number of relational constraints on each cell face.
//  // Use the first solution variable with physical constraints to determine this number.
//  parameter = ParametersWithPhysicalBCs[0];
//  TotalNumberOfRelationalBCsEquations = (RelationalBCs_W[parameter] + RelationalBCs_S[parameter] +
//					 RelationalBCs_E[parameter] + RelationalBCs_N[parameter]);
//  TotalNumberOfExactlySatisfiedEquations += TotalNumberOfRelationalBCsEquations;
//
//  // Calculate TotalNumberOfEquations 
//  // Add TotalNumberOfExactlySatisfiedEquations and the number of approximate equations for all related variables
//  TotalNumberOfEquations = TotalNumberOfExactlySatisfiedEquations + (StencilSize-1) * ParametersWithPhysicalBCs.size();
//
//  /******** Set dimensions for the matrices of the least-squares problem **********/
//  /********************************************************************************/
//  A_Assembled.newsize(TotalNumberOfEquations, ParametersWithPhysicalBCs.size()*NumberOfTaylorDerivatives());
//  A_Assembled.zero();
//  All_U_Assembled.newsize(TotalNumberOfEquations, 1); // There is only ONE column
//  All_U_Assembled.zero();
//  X.newsize(ParametersWithPhysicalBCs.size()*NumberOfTaylorDerivatives(), 1); // There is only ONE column
//
//
//  /******** Fetch the data on each cell edge for the exactly satisfied RELATIONAL constraints ************/
//  /*******************************************************************************************************/
//  parameter = ParametersWithPhysicalBCs[0];
//  FetchDataConstraints(SolnBlk, iCell, jCell,
//		       BC_Type, parameter,
//		       ConstrainedGQPs_West,  RelationalBCs_W,
//		       ConstrainedGQPs_South, RelationalBCs_S,
//		       ConstrainedGQPs_East,  RelationalBCs_E,
//		       ConstrainedGQPs_North, RelationalBCs_N,
//		       Constraints_Loc, Constraints_Normals, Constraints_BCs);
//
//  /************ Generate exactly satisfied relational constraints for COUPLED variables *************/
//  /************ Only one type of relational constraint can be accommodated at a time *****************/
//  Generalized_RelationalConstraints_Equations(SolnBlk,
//  					      iCell, jCell,
//  					      Constraints_Loc,
//  					      Constraints_Normals,
//  					      Constraints_BCs,
//					      BC_Type,
//  					      A_Assembled, All_U_Assembled,
//  					      ParametersWithPhysicalBCs,
//  					      0, 0);
//
//  /*** Set the individual constraints matrix entries ****/
//  /******************************************************/
//  // === Initialize matrix entry indexes
//  MeanConservationRow = TotalNumberOfRelationalBCsEquations;  
//  ApproxMeanConservationRow = TotalNumberOfExactlySatisfiedEquations;
//
//  for (iterator = 0;  iterator < ParametersWithPhysicalBCs.size(); ++iterator){
//    
//    // === Reset and initialize the local variables ===
//    parameter = ParametersWithPhysicalBCs[iterator]; // index of current state variable
//    ParameterIndex[0] = parameter;
//
//    /** Generate the exact and approximate mean conservation equations for each parameter in the ParametersWithPhysicalBCs **/
//    /************************************************************************************************************************/
//    if ( IsPseudoInverseAllocated() && IsPseudoInversePreComputed() ){  
//
//      // Copy the matrix (i.e. the least-squares equations + the mean conservation)
//      // into the assembled matrix for each parameter.
//
//      // ===== Copy mean conservation equation and set RHS =====
//      A_Assembled.incorporate_matrix(MeanConservationRow,MeanConservationCol,
//				     Cell_LHS_Inv(iCell,jCell),
//				     0, StencilSize-1); //< copy only the first row in the Cell_LHS_Inv matrix
//      
//      All_U_Assembled(MeanConservationRow,0) = (SolnBlk.*ReconstructedSoln)(iCell,jCell)[parameter];
//      
//      // ===== Copy approximate mean conservation equations and set RHS =====
//      A_Assembled.incorporate_matrix(ApproxMeanConservationRow,MeanConservationCol,
//				     Cell_LHS_Inv(iCell,jCell),
//				     1); //< skip only the first row in the Cell_LHS_Inv matrix
//      
//      // Form the RHS for the least-squares problem
//      for (cell=1; cell<StencilSize; ++cell) { //for each cell in the stencil
//	All_U_Assembled(ApproxMeanConservationRow + cell-1, 0) = ( GeomWeightValue(iCell,jCell,cell)*
//								   (SolnBlk.*ReconstructedSoln)(i_index[cell],
//												j_index[cell])[parameter] );
//      } //endfor (cell)
//
//    } else {
//
//      // Set the mean conservation equation for the current parameter
//      Set_MeanValueConservation_Equations(SolnBlk,
//					  ReconstructedSoln,
//					  iCell,jCell,
//					  i_index, j_index,
//					  A_Assembled, All_U_Assembled,
//					  ParameterIndex,
//					  MeanConservationRow,
//					  ApproxMeanConservationRow, MeanConservationCol);
//    } // endif
//
//
//    // Update entry indexes for the next parameter
//    ++MeanConservationRow;
//    MeanConservationCol += NumberOfTaylorDerivatives();
//    ApproxMeanConservationRow += StencilSize - 1;
//  } // endfor (iterator)
//
//
//  /** Obtain solution to the linear equality-constrained least-square problem ********/
//  /***********************************************************************************/
//  Solve_Constrained_LS_Householder(A_Assembled,
//				   All_U_Assembled,
//				   X,
//				   TotalNumberOfExactlySatisfiedEquations);
//  
//  /** Update the coefficients D (derivatives) ********/
//  /***************************************************/
//  for (parameter = 0; parameter < ParametersWithPhysicalBCs.size(); ++parameter){
//    MeanConservationRow = parameter*NumberOfTaylorDerivatives(); // < calculate the shift in X vector
//    for (n=0; n<=CellTaylorDeriv(iCell,jCell).LastElem(); ++n){
//      CellTaylorDerivState(iCell,jCell,n)[ParametersWithPhysicalBCs[parameter]] = X(MeanConservationRow + n,0);
//    }//endfor
//  }//endfor
//
//  // Deallocate memory
//  delete [] RelationalBCs_W;
//  delete [] RelationalBCs_S;
//  delete [] RelationalBCs_E;
//  delete [] RelationalBCs_N;
//}
//
//
///*! 
// * Compute the unlimited k-exact high-order reconstruction
// * proposed by Barth (1993) combined with exactly satisfied equations 
// * (i.e. only individual constraints) which account for boundary conditions.
// * This reconstruction procedure should be used for computational cells
// * affected by the presence of a boundary condition enforced with individual constraints.
// * Individual constraints represent exact equations in only one solution parameter.
// * The mean quantity conservation is also enforced as a constraint.
// */
//template<class SOLN_STATE>
//template<class Soln_Block_Type> inline
//void HighOrder<SOLN_STATE>::
//ComputeIndividuallyConstrainedUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk,
//							      const Soln_State & 
//							      (Soln_Block_Type::*ReconstructedSoln)(const int &,
//												    const int &) const,
//							      const int &iCell, const int &jCell,
//							      const IndexType & i_index,
//							      const IndexType & j_index,
//							      const int & ConstrainedGQPs_West,
//							      const int & ConstrainedGQPs_South,
//							      const int & ConstrainedGQPs_East,
//							      const int & ConstrainedGQPs_North) {
//
//  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS
//  IndexType ParametersWithPhysicalBCs,   //!< List of solution parameters that have physical BCs (i.e. have individual constraints)
//    ParametersWithNumericalBCs,		 //!< List of solution parameters that have numerical BCs (i.e. no constraints)
//    ParameterIndex(1);			 //!< List of constrained solution parameters (used for current calculations)
//  //! List for West edge of BCtype for each parameter (i.e. how many conditions per edge)
//  int *IndividualBCs_W = new int [NumberOfVariables()+1];
//  //! List for South edge of BCtype for each parameter (i.e. how many conditions per edge)
//  int *IndividualBCs_S = new int [NumberOfVariables()+1];
//  //! List for East edge of BCtype for each parameter (i.e. how many conditions per edge)
//  int *IndividualBCs_E = new int [NumberOfVariables()+1];
//  //! List for North edge of BCtype for each parameter (i.e. how many conditions per edge)
//  int *IndividualBCs_N = new int [NumberOfVariables()+1];
//  int parameter;
//  int TotalNumberOfEquations, TotalNumberOfExactlySatisfiedEquations, StartRow;
//  int TotalIndividual;
//  int cell, n, iterator;
//  int BC_Type;			//< defined only for compatibility
//  
//
//  /********************************************************************************************************
//   *  STEP 1. SORT THE PARAMETERS INTO THE DESIGNATED CATEGORIES AND DETECT CONSTRAINTS ON EACH CELL EDGE *
//   *******************************************************************************************************/
//  for (parameter = 1; parameter <= NumberOfVariables(); ++parameter){
//
//    // (Re)set variables
//    TotalIndividual = 0;
//
//    IndividualBCs_W[parameter] = 0;
//    IndividualBCs_S[parameter] = 0;
//    IndividualBCs_E[parameter] = 0;
//    IndividualBCs_N[parameter] = 0;
//
//    // West edge
//    if (ConstrainedGQPs_West > 0){
//      IndividualBCs_W[parameter] = ConstrainedGQPs_West * SolnBlk.BC_WestCell(jCell).NumberOfIndividualConstraints(parameter);
//      TotalIndividual += IndividualBCs_W[parameter];
//    }
//
//    // South edge
//    if (ConstrainedGQPs_South > 0){
//      IndividualBCs_S[parameter] = ConstrainedGQPs_South * SolnBlk.BC_SouthCell(iCell).NumberOfIndividualConstraints(parameter);
//      TotalIndividual += IndividualBCs_S[parameter];
//    }
//
//    // East edge
//    if (ConstrainedGQPs_East > 0){
//      IndividualBCs_E[parameter] = ConstrainedGQPs_East * SolnBlk.BC_EastCell(jCell).NumberOfIndividualConstraints(parameter);
//      TotalIndividual += IndividualBCs_E[parameter];
//    }
//
//    // North edge
//    if (ConstrainedGQPs_North > 0){
//      IndividualBCs_N[parameter] = ConstrainedGQPs_North * SolnBlk.BC_NorthCell(iCell).NumberOfIndividualConstraints(parameter);
//      TotalIndividual += IndividualBCs_N[parameter];
//    }
//    
//    if (TotalIndividual){
//      // At least one edge has individual constraints for this parameter
//      ParametersWithPhysicalBCs.push_back(parameter);
//    } else {
//      // There are no constraints for this parameter.
//      ParametersWithNumericalBCs.push_back(parameter);
//    }
//  }// endfor
//
//
//  /****************************************************************************
//   *  STEP 2. SOLVE THE RECONSTRUCTION FOR PARAMETERS THAT HAVE NUMERICAL BCs *
//   ***************************************************************************/
//  ComputeUnconstrainedUnlimitedSolutionReconstructionInConstrainedCell(SolnBlk, ReconstructedSoln,
//								       iCell, jCell,
//								       i_index, j_index,
//								       ParametersWithNumericalBCs);
//  
//  
//  /***************************************************************************
//   *  STEP 3. SOLVE THE RECONSTRUCTION FOR PARAMETERS THAT HAVE PHYSICAL BCs *
//   **************************************************************************/
//  for (iterator = 0; iterator < ParametersWithPhysicalBCs.size(); ++iterator){
//    /*** Perform solution reconstruction for each parameter with individual constraints (i.e. physical BCs) *********/
//    /****************************************************************************************************************/
//
//    // === Reset and initialize the local variables ===
//    parameter = ParametersWithPhysicalBCs[iterator]; // index of current state variable    
//    // Calculate TotalNumberOfExactlySatisfiedEquations and account for average conservation in cell (iCell,jCell)
//    TotalNumberOfExactlySatisfiedEquations = (1 + IndividualBCs_W[parameter] + IndividualBCs_S[parameter] + 
//					      IndividualBCs_E[parameter] + IndividualBCs_N[parameter]);
//    // Calculate TotalNumberOfEquations (i.e. TotalNumberOfExactlySatisfiedEquations + approximate equations)
//    TotalNumberOfEquations = TotalNumberOfExactlySatisfiedEquations + (i_index.size() - 1);
//
//    /******** Determine dimensions of the least-squares problem and set matrices accordingly ************/
//    /****************************************************************************************************/    
//    A_Assembled.newsize(TotalNumberOfEquations, NumberOfTaylorDerivatives());
//    All_U_Assembled.newsize(TotalNumberOfEquations, 1); // There is only ONE column
//
//    /******** Fetch the data on each cell edge for the exactly satisfied individual constraints ************/
//    /*******************************************************************************************************/
//    FetchDataConstraints(SolnBlk, iCell, jCell,
//			 BC_Type, parameter,
//			 ConstrainedGQPs_West,  IndividualBCs_W,
//			 ConstrainedGQPs_South, IndividualBCs_S,
//			 ConstrainedGQPs_East,  IndividualBCs_E,
//			 ConstrainedGQPs_North, IndividualBCs_N,
//			 Constraints_Loc, Constraints_Normals, Constraints_BCs);
//
//    /************ Generate exactly satisfied individual constraints for UNCOUPLED variables *************/
//    /****************************************************************************************************/
//    ParameterIndex[0] = parameter;
//    Generalized_IndividualConstraints_Equations(SolnBlk,
//						iCell, jCell,
//						Constraints_Loc,
//						Constraints_Normals,
//						Constraints_BCs,
//						A_Assembled, All_U_Assembled,
//						ParameterIndex,
//						0, 0);
//
//    /******************** Generate the exact and approximate mean conservation equations ***************************/
//    /***************************************************************************************************************/
//    if ( IsPseudoInverseAllocated() && IsPseudoInversePreComputed() ){
//
//      // Copy the matrix (i.e. the least-squares equations + the mean conservation) into the assembled matrix
//      // starting at the right position
//      StartRow = TotalNumberOfExactlySatisfiedEquations-1;
//      A_Assembled.incorporate_matrix(StartRow,0,
//				     Cell_LHS_Inv(iCell,jCell));
//
//      
//      // Form the RHS for the least-squares problem
//      for (cell=0; cell<i_index.size(); ++cell) { //for each cell in the stencil
//	All_U_Assembled(cell+StartRow,0) = ( GeomWeightValue(iCell,jCell,cell)*
//					     (SolnBlk.*ReconstructedSoln)(i_index[cell],
//									  j_index[cell])[parameter] );
//      } //endfor (cell)
//
//    } else {
//      // Set the mean conservation equation for the current parameter
//      Set_MeanValueConservation_Equations(SolnBlk,
//					  ReconstructedSoln,
//					  iCell,jCell,
//					  i_index, j_index,
//					  A_Assembled, All_U_Assembled,
//					  ParameterIndex,
//					  TotalNumberOfExactlySatisfiedEquations-1,
//					  TotalNumberOfExactlySatisfiedEquations,
//					  0);
//    } // endif
//
//
//    /** Obtain solution to the linear equality-constrained least-square problem ********/
//    /***********************************************************************************/
//    // Note: X_Assembled might have more columns than necessary, but that's not a problem!
//    Solve_Constrained_LS_Householder(A_Assembled,
//				     All_U_Assembled,
//				     X_Assembled,
//				     TotalNumberOfExactlySatisfiedEquations);
//
//    /** Update the coefficients D (derivatives) ********/
//    /***************************************************/
//    for (n=0; n<=CellTaylorDeriv(iCell,jCell).LastElem(); ++n){
//      CellTaylorDerivState(iCell,jCell,n)[ParameterIndex[0]] = X_Assembled(n,0);
//    }//endfor
//
//  } // endfor (iterator)
//
//  // Deallocate memory
//  delete [] IndividualBCs_W;
//  delete [] IndividualBCs_S;
//  delete [] IndividualBCs_E;
//  delete [] IndividualBCs_N;
//}
//
///*! 
// * Compute the unlimited k-exact high-order reconstruction
// * proposed by Barth (1993) combined with exactly satisfied equations 
// * (i.e. relational and individual constraints) which account for boundary conditions.
// * The mean quantity conservation is also enforced as a constraint for all solution parameters.
// *
// * \param SolnBlk the quad block for which the solution reconstruction is done.
// * \param ReconstructedSoln member function of SolnBlk which returns the average solution of (iCell,jCell)
// * \param iCell i-index of the reconstructed cell
// * \param jCell j-index of the reconstructed cell
// * \param i_index i-indexes of the cells that are part of the supporting stencil
// * \param j_index j-indexes of the cells that are part of the supporting stencil
// * \param ConstrainedGQPs_West number of Gauss quadrature points constrained on West cell edge
// * \param ConstrainedGQPs_South number of Gauss quadrature points constrained on South cell edge
// * \param ConstrainedGQPs_East number of Gauss quadrature points constrained on East cell edge
// * \param ConstrainedGQPs_North number of Gauss quadrature points constrained on North cell edge
// */
//template<class SOLN_STATE>
//template<class Soln_Block_Type> inline
//void HighOrder<SOLN_STATE>::
//ComputeRelationallyAndIndividuallyConstrainedUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk,
//									     const Soln_State & 
//									     (Soln_Block_Type::
//									      *ReconstructedSoln)(const int &,
//												  const int &) const,
//									     const int &iCell, const int &jCell,
//									     const IndexType & i_index,
//									     const IndexType & j_index,
//									     const int & ConstrainedGQPs_West,
//									     const int & ConstrainedGQPs_South,
//									     const int & ConstrainedGQPs_East,
//									     const int & ConstrainedGQPs_North) {
//
//  /* Notes on the algorithm:
//   * 
//   * The solution variables are sorted into 3 categories:
//   *    A. Variables with physical BCs imposed with relational constraints and individual (if they exist) constraints
//   *    B. Variables with physical BCs imposed with only individual constraints
//   *    C. Variables with numerical BCs (i.e. what doens't belong to A or B)
//   */
//  
//  
//  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS
//  IndexType ParametersWithOnlyIndividualBCs,  /*!< List of solution parameters in category B */
//    ParametersWithRelationalBCs,         /*!< List of solution parameters in category A */
//    ParametersWithNumericalBCs,		 //!< List of solution parameters that have numerical BCs (i.e. no constraints)
//    ParameterIndex(1);			 //!< List of constrained solution parameters (used for current calculations)
//
//  //! List for West edge of relational BCtype for each parameter (i.e. how many conditions per edge)
//  int *RelationalBCs_W = new int [NumberOfVariables()+1];
//  //! List for South edge of relational BCtype for each parameter (i.e. how many conditions per edge)
//  int *RelationalBCs_S = new int [NumberOfVariables()+1];
//  //! List for East edge of relational BCtype for each parameter (i.e. how many conditions per edge)
//  int *RelationalBCs_E = new int [NumberOfVariables()+1];
//  //! List for North edge of relational BCtype for each parameter (i.e. how many conditions per edge)
//  int *RelationalBCs_N = new int [NumberOfVariables()+1];
//  //! List for West edge of individual BCtype for each parameter (i.e. how many conditions per edge)
//  int *IndividualBCs_W = new int [NumberOfVariables()+1];
//  //! List for South edge of individual BCtype for each parameter (i.e. how many conditions per edge)
//  int *IndividualBCs_S = new int [NumberOfVariables()+1];
//  //! List for East edge of individual BCtype for each parameter (i.e. how many conditions per edge)
//  int *IndividualBCs_E = new int [NumberOfVariables()+1];
//  //! List for North edge of individual BCtype for each parameter (i.e. how many conditions per edge)
//  int *IndividualBCs_N = new int [NumberOfVariables()+1];
//
//  int parameter, iterator, StartRow, StartCol;
//  int TotalNumberOfEquations, TotalNumberOfExactlySatisfiedEquations;
//  int ApproxMeanConservationRow, MeanConservationRow, MeanConservationCol(0);
//  int cell, n, BC_Type;
//  int StencilSize(i_index.size());
//  DenseMatrix X;		//< the matrix of unknowns
//  int TotalRelational, TotalIndividual;  
//  int TotalNumberOfRelationalBCsEquations(0), TotalNumberOfIndividualBCsEquations(0);
//
//
//  /********************************************************************************************************
//   *  STEP 1. SORT THE PARAMETERS INTO THE DESIGNATED CATEGORIES AND DETECT CONSTRAINTS ON EACH CELL EDGE *
//   *******************************************************************************************************/
//  for (parameter = 1; parameter <= NumberOfVariables(); ++parameter){
//
//    // (Re)set variables
//    TotalRelational = 0;
//    TotalIndividual = 0;
//
//    RelationalBCs_W[parameter] = 0;
//    RelationalBCs_S[parameter] = 0;
//    RelationalBCs_E[parameter] = 0;
//    RelationalBCs_N[parameter] = 0;
//
//    IndividualBCs_W[parameter] = 0;
//    IndividualBCs_S[parameter] = 0;
//    IndividualBCs_E[parameter] = 0;
//    IndividualBCs_N[parameter] = 0;
//
//    // West edge
//    if (ConstrainedGQPs_West > 0){
//      RelationalBCs_W[parameter] = ConstrainedGQPs_West * SolnBlk.BC_WestCell(jCell).NumberOfRelationalConstraints(parameter);
//      TotalRelational += RelationalBCs_W[parameter];
//      IndividualBCs_W[parameter] = ConstrainedGQPs_West * SolnBlk.BC_WestCell(jCell).NumberOfIndividualConstraints(parameter);
//      TotalIndividual += IndividualBCs_W[parameter];
//    }
//
//    // South edge
//    if (ConstrainedGQPs_South > 0){
//      RelationalBCs_S[parameter] = ConstrainedGQPs_South * SolnBlk.BC_SouthCell(iCell).NumberOfRelationalConstraints(parameter);
//      TotalRelational += RelationalBCs_S[parameter];
//      IndividualBCs_S[parameter] = ConstrainedGQPs_South * SolnBlk.BC_SouthCell(iCell).NumberOfIndividualConstraints(parameter);
//      TotalIndividual += IndividualBCs_S[parameter];
//    }
//
//    // East edge
//    if (ConstrainedGQPs_East > 0){
//      RelationalBCs_E[parameter] = ConstrainedGQPs_East * SolnBlk.BC_EastCell(jCell).NumberOfRelationalConstraints(parameter);
//      TotalRelational += RelationalBCs_E[parameter];
//      IndividualBCs_E[parameter] = ConstrainedGQPs_East * SolnBlk.BC_EastCell(jCell).NumberOfIndividualConstraints(parameter);
//      TotalIndividual += IndividualBCs_E[parameter];
//    }
//
//    // North edge
//    if (ConstrainedGQPs_North > 0){
//      RelationalBCs_N[parameter] = ConstrainedGQPs_North * SolnBlk.BC_NorthCell(iCell).NumberOfRelationalConstraints(parameter);
//      TotalRelational += RelationalBCs_N[parameter];
//      IndividualBCs_N[parameter] = ConstrainedGQPs_North * SolnBlk.BC_NorthCell(iCell).NumberOfIndividualConstraints(parameter);
//      TotalIndividual += IndividualBCs_N[parameter];
//    }
//    
//    if (TotalRelational){
//      // At least one edge has relational constraints for this parameter
//      ParametersWithRelationalBCs.push_back(parameter);
//    } else if (TotalIndividual){
//      // At least one edge has individual constraints for this parameter but no edge has relational constraints
//      ParametersWithOnlyIndividualBCs.push_back(parameter);
//    } else {
//      // There are no constraints for this parameter.
//      ParametersWithNumericalBCs.push_back(parameter);
//    }
//  }// endfor
//
//
//  /****************************************************************************
//   *  STEP 2. SOLVE THE RECONSTRUCTION FOR PARAMETERS THAT HAVE NUMERICAL BCs *
//   ***************************************************************************/
//  ComputeUnconstrainedUnlimitedSolutionReconstructionInConstrainedCell(SolnBlk, ReconstructedSoln,
//								       iCell, jCell,
//								       i_index, j_index,
//								       ParametersWithNumericalBCs);
//
//
//
//  /*********************************************************************************
//   *  STEP 3. SOLVE THE RECONSTRUCTION FOR PARAMETERS THAT HAVE ONY INDIVIDUAL BCs *
//   *********************************************************************************/
//  for (iterator = 0;  iterator < ParametersWithOnlyIndividualBCs.size(); ++iterator){
//    /*** Perform solution reconstruction for each parameter with only individual constraints  ****/
//    /*********************************************************************************************/
//
//    // === Reset and initialize the local variables ===
//    parameter = ParametersWithOnlyIndividualBCs[iterator]; // index of current state variable
//    // Calculate TotalNumberOfExactlySatisfiedEquations and account for average conservation in cell (iCell,jCell)
//    TotalNumberOfExactlySatisfiedEquations = (1 + IndividualBCs_W[parameter] + IndividualBCs_S[parameter] + 
//					      IndividualBCs_E[parameter] + IndividualBCs_N[parameter]);
//    // Calculate TotalNumberOfEquations (i.e. TotalNumberOfExactlySatisfiedEquations + approximate equations)
//    TotalNumberOfEquations = TotalNumberOfExactlySatisfiedEquations + (StencilSize - 1);
//
//    /******** Determine dimensions of the least-squares problem and set matrices accordingly ************/
//    /****************************************************************************************************/    
//    A_Assembled.newsize(TotalNumberOfEquations, NumberOfTaylorDerivatives());
//    All_U_Assembled.newsize(TotalNumberOfEquations, 1); // There is only ONE column
//
//    /******** Fetch the data on each cell edge for the exactly satisfied individual constraints ************/
//    /*******************************************************************************************************/
//    FetchDataConstraints(SolnBlk, iCell, jCell,
//			 BC_Type, parameter,
//			 ConstrainedGQPs_West,  IndividualBCs_W,
//			 ConstrainedGQPs_South, IndividualBCs_S,
//			 ConstrainedGQPs_East,  IndividualBCs_E,
//			 ConstrainedGQPs_North, IndividualBCs_N,
//			 Constraints_Loc, Constraints_Normals, Constraints_BCs);
//
//    /************ Generate exactly satisfied individual constraints for UNCOUPLED variables *************/
//    /****************************************************************************************************/
//    ParameterIndex[0] = parameter;
//    Generalized_IndividualConstraints_Equations(SolnBlk,
//						iCell, jCell,
//						Constraints_Loc,
//						Constraints_Normals,
//						Constraints_BCs,
//						A_Assembled, All_U_Assembled,
//						ParameterIndex,
//						0, 0);
//
//    /******************** Generate the exact and approximate mean conservation equations ***************************/
//    /***************************************************************************************************************/
//    if ( IsPseudoInverseAllocated() && IsPseudoInversePreComputed() ){
//
//      // Copy the matrix (i.e. the least-squares equations + the mean conservation) into the assembled matrix
//      // starting at the right position
//      StartRow = TotalNumberOfExactlySatisfiedEquations-1;
//      A_Assembled.incorporate_matrix(StartRow,0,
//				     Cell_LHS_Inv(iCell,jCell));
//
//      
//      // Form the RHS for the least-squares problem
//      for (cell=0; cell<i_index.size(); ++cell) { //for each cell in the stencil
//	All_U_Assembled(cell+StartRow,0) = ( GeomWeightValue(iCell,jCell,cell)*
//					     (SolnBlk.*ReconstructedSoln)(i_index[cell],
//									  j_index[cell])[ParameterIndex[0]] );
//      } //endfor (cell)
//
//    } else {
//      // Set the mean conservation equation for the current parameter
//      Set_MeanValueConservation_Equations(SolnBlk,
//					  ReconstructedSoln,
//					  iCell,jCell,
//					  i_index, j_index,
//					  A_Assembled, All_U_Assembled,
//					  ParameterIndex,
//					  TotalNumberOfExactlySatisfiedEquations-1,
//					  TotalNumberOfExactlySatisfiedEquations,
//					  0);
//    } // endif
//    
//    
//    /** Obtain solution to the linear equality-constrained least-square problem ********/
//    /***********************************************************************************/
//    // Note: X_Assembled might have more columns than necessary, but that's not a problem!
//    Solve_Constrained_LS_Householder(A_Assembled,
//				     All_U_Assembled,
//				     X_Assembled,
//				     TotalNumberOfExactlySatisfiedEquations);
//    
//    /** Update the coefficients D (derivatives) ********/
//    /***************************************************/
//    for (n=0; n<=CellTaylorDeriv(iCell,jCell).LastElem(); ++n){
//      CellTaylorDerivState(iCell,jCell,n)[ParameterIndex[0]] = X_Assembled(n,0);
//    }//endfor
//    
//  } // endfor (iterator)
//
//
//  /********************************************************************************************
//   *  STEP 4. SOLVE THE RECONSTRUCTION FOR PARAMETERS THAT HAVE RELATIONAL AND INDIVIDUAL BCs *
//   *******************************************************************************************/
//  /* Note: It is assumed that each GQP introduces only one relational constraint
//           that relates all parameters in the ParametersWithRelationalBCs vector.
//	   There is only one type of relational constraint in a given cell!
//  */
//
//  if (ParametersWithRelationalBCs.size() == 0){
//    throw runtime_error("HighOrder<SOLN_STATE>::ComputeRelationallyAndIndividuallyConstrainedUnlimitedSolutionReconstruction() ERROR! Reconstruction with relational constraints is required but there are no relational constraints present!!");
//  }
//
//  // === Reset and initialize the local variables ===
//  // Calculate TotalNumberOfExactlySatisfiedEquations
//  // Initialize with the number of mean conservation constraints
//  TotalNumberOfExactlySatisfiedEquations = ParametersWithRelationalBCs.size();
//
//  // Add number of relational constraints on each cell face.
//  // Use the first solution variable with relational constraints to determine this number.
//  parameter = ParametersWithRelationalBCs[0];
//  TotalNumberOfRelationalBCsEquations += (RelationalBCs_W[parameter] + RelationalBCs_S[parameter] +
//					  RelationalBCs_E[parameter] + RelationalBCs_N[parameter]);
//  TotalNumberOfExactlySatisfiedEquations += TotalNumberOfRelationalBCsEquations;
//  
//  // Add individual constraints required by each solution variable in ParametersWithRelationalBCs on each cell face.
//  for (iterator = 0; iterator < ParametersWithRelationalBCs.size(); ++iterator){
//    parameter = ParametersWithRelationalBCs[iterator];
//    
//    TotalNumberOfIndividualBCsEquations += (IndividualBCs_W[parameter] + IndividualBCs_S[parameter] +
//					    IndividualBCs_E[parameter] + IndividualBCs_N[parameter]);
//  }
//  TotalNumberOfExactlySatisfiedEquations += TotalNumberOfIndividualBCsEquations;
//  
//  // Calculate TotalNumberOfEquations 
//  // Add TotalNumberOfExactlySatisfiedEquations and the number of approximate equations for all related variables
//  TotalNumberOfEquations = TotalNumberOfExactlySatisfiedEquations + (StencilSize-1) * ParametersWithRelationalBCs.size();
//
//
//  /******** Set dimensions for the matrices of the least-squares problem **********/
//  /********************************************************************************/
//  A_Assembled.newsize(TotalNumberOfEquations, ParametersWithRelationalBCs.size()*NumberOfTaylorDerivatives());
//  A_Assembled.zero();
//  All_U_Assembled.newsize(TotalNumberOfEquations, 1); // There is only ONE column
//  All_U_Assembled.zero();
//  X.newsize(ParametersWithRelationalBCs.size()*NumberOfTaylorDerivatives(), 1); // There is only ONE column
//
//
//  /******** Fetch the data on each cell edge for the exactly satisfied RELATIONAL constraints ************/
//  /*******************************************************************************************************/
//  parameter = ParametersWithRelationalBCs[0];
//  FetchDataConstraints(SolnBlk, iCell, jCell,
//		       BC_Type, parameter,
//		       ConstrainedGQPs_West,  RelationalBCs_W,
//		       ConstrainedGQPs_South, RelationalBCs_S,
//		       ConstrainedGQPs_East,  RelationalBCs_E,
//		       ConstrainedGQPs_North, RelationalBCs_N,
//		       Constraints_Loc, Constraints_Normals, Constraints_BCs);
//
//  /************ Generate exactly satisfied relational constraints for COUPLED variables *************/
//  /************ Only one type of relational constraint can be accommodated at a time *****************/
//  Generalized_RelationalConstraints_Equations(SolnBlk,
//  					      iCell, jCell,
//  					      Constraints_Loc,
//  					      Constraints_Normals,
//  					      Constraints_BCs,
//					      BC_Type,
//  					      A_Assembled, All_U_Assembled,
//  					      ParametersWithRelationalBCs,
//  					      0, 0);
//
//  
//  /*** Set the individual constraints matrix entries ****/
//  /******************************************************/
//  StartRow = TotalNumberOfRelationalBCsEquations;
//  StartCol = 0; 
//  MeanConservationRow = TotalNumberOfRelationalBCsEquations + TotalNumberOfIndividualBCsEquations;
//  ApproxMeanConservationRow = TotalNumberOfExactlySatisfiedEquations;
//  for (iterator = 0;  iterator < ParametersWithRelationalBCs.size(); ++iterator){
//
//    // === Reset and initialize the local variables ===
//    parameter = ParametersWithRelationalBCs[iterator]; // index of current state variable
//    ParameterIndex[0] = parameter;
//
//    /******** Fetch the data on each cell edge for the exactly satisfied individual constraints ************/
//    /*******************************************************************************************************/
//    FetchDataConstraints(SolnBlk, iCell, jCell,
//			 BC_Type, parameter,
//			 ConstrainedGQPs_West,  IndividualBCs_W,
//			 ConstrainedGQPs_South, IndividualBCs_S,
//			 ConstrainedGQPs_East,  IndividualBCs_E,
//			 ConstrainedGQPs_North, IndividualBCs_N,
//			 Constraints_Loc, Constraints_Normals, Constraints_BCs);
//
//    /************ Generate exactly satisfied individual constraints for UNCOUPLED variables *************/
//    /****************************************************************************************************/
//    Generalized_IndividualConstraints_Equations(SolnBlk,
//						iCell, jCell,
//						Constraints_Loc,
//						Constraints_Normals,
//						Constraints_BCs,
//						A_Assembled, All_U_Assembled,
//						ParameterIndex,
//						StartRow, StartCol);
//    
//
//    /** Generate the exact and approximate mean conservation equations for each parameter **/
//    /***************************************************************************************/
//    if ( IsPseudoInverseAllocated() && IsPseudoInversePreComputed() ){  
//
//      // Copy the matrix (i.e. the least-squares equations + the mean conservation)
//      // into the assembled matrix for each parameter.
//
//      // ===== Copy mean conservation equation and set RHS =====
//      A_Assembled.incorporate_matrix(MeanConservationRow,StartCol,
//				     Cell_LHS_Inv(iCell,jCell),
//				     0, StencilSize-1); //< copy only the first row in the Cell_LHS_Inv matrix
//      
//      All_U_Assembled(MeanConservationRow,0) = (SolnBlk.*ReconstructedSoln)(iCell,jCell)[parameter];
//      
//      // ===== Copy approximate mean conservation equations and set RHS =====
//      A_Assembled.incorporate_matrix(ApproxMeanConservationRow,StartCol,
//				     Cell_LHS_Inv(iCell,jCell),
//				     1); //< skip only the first row in the Cell_LHS_Inv matrix
//      
//      // Form the RHS for the least-squares problem
//      for (cell=1; cell<StencilSize; ++cell) { //for each cell in the stencil
//	All_U_Assembled(ApproxMeanConservationRow + cell-1, 0) = ( GeomWeightValue(iCell,jCell,cell)*
//								   (SolnBlk.*ReconstructedSoln)(i_index[cell],
//												j_index[cell])[parameter] );
//      } //endfor (cell)
//      
//    } else {
//
//      // Set the mean conservation equation for the current parameter
//      Set_MeanValueConservation_Equations(SolnBlk,
//					  ReconstructedSoln,
//					  iCell,jCell,
//					  i_index, j_index,
//					  A_Assembled, All_U_Assembled,
//					  ParameterIndex,
//					  MeanConservationRow,
//					  ApproxMeanConservationRow, StartCol);
//    }// endif
//
//
//    // Update entry indexes for the next parameter
//    StartRow += ( IndividualBCs_W[parameter] + IndividualBCs_S[parameter] + 
//		  IndividualBCs_E[parameter] + IndividualBCs_N[parameter] );
//    StartCol += NumberOfTaylorDerivatives();
//    ++MeanConservationRow;    
//    ApproxMeanConservationRow += StencilSize - 1;
//
//  } // endfor (iterator)
//
//
//  /** Obtain solution to the linear equality-constrained least-square problem ********/
//  /***********************************************************************************/
//  Solve_Constrained_LS_Householder(A_Assembled,
//				   All_U_Assembled,
//				   X,
//				   TotalNumberOfExactlySatisfiedEquations);
//  
//  /** Update the coefficients D (derivatives) ********/
//  /***************************************************/
//  for (iterator = 0; iterator < ParametersWithRelationalBCs.size(); ++iterator){
//    MeanConservationRow = iterator*NumberOfTaylorDerivatives(); // < calculate the shift in X vector
//    for (n=0; n<=CellTaylorDeriv(iCell,jCell).LastElem(); ++n){
//      CellTaylorDerivState(iCell,jCell,n)[ParametersWithRelationalBCs[iterator]] = X(MeanConservationRow + n,0);
//    }//endfor
//  }//endfor
//
//
//  // Deallocate memory
//  delete [] RelationalBCs_W;
//  delete [] RelationalBCs_S;
//  delete [] RelationalBCs_E;
//  delete [] RelationalBCs_N;
//  delete [] IndividualBCs_W;
//  delete [] IndividualBCs_S;
//  delete [] IndividualBCs_E;
//  delete [] IndividualBCs_N;
//}
//
///*! 
// * Compute the unlimited k-exact high-order reconstruction
// * proposed by Barth (1993) for a set of solution parameters 
// * that have imposed only numerical boundary conditions.
// * This reconstruction procedure should be used for computational cells
// * next to the boundary and which have boundary conditions enforced with constraints.
// * Due to the physics of the BCs several parameters might not have any constraints enforced.
// * This routine is designed for these parameters.
// * The mean quantity conservation is enforced as a constraint.
// */
//template<class SOLN_STATE>
//template<class Soln_Block_Type> inline
//void HighOrder<SOLN_STATE>::
//ComputeUnconstrainedUnlimitedSolutionReconstructionInConstrainedCell(Soln_Block_Type &SolnBlk, 
//								     const Soln_State & 
//								     (Soln_Block_Type::*ReconstructedSoln)(const int &,
//													   const int &) const,
//								     const int &iCell, const int &jCell,
//								     const IndexType & i_index, const IndexType & j_index,
//								     const IndexType & ParameterIndex){
//
//  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS
//
//  int StencilSize(i_index.size()); 
//  int cell, i, parameter;
//  int P1, P2;
//  int ParameterIndexSize(ParameterIndex.size());
//  int n;
//
//  // === Check if there are any variables for which this reconstruction must be performed ===
//  if (ParameterIndexSize == 0){
//    return;
//  }
//
//  // Check if the full matrix has been allocated and pre-computed
//  if ( IsPseudoInverseAllocated() && IsPseudoInversePreComputed() ){
//    // Use the full matrix to speedup the least-squares reconstruction
//
//    // Step 1. Copy the full matrix (i.e. the least-squares equations + the mean conservation) into a local variable
//    A_Assembled = Cell_LHS_Inv(iCell,jCell);
//    
//    // Step 2. Set the RHS term
//    // Ensure that All_U_Assembled matrix has proper dimensions
//    if ( All_U_Assembled.size(0) != StencilSize || All_U_Assembled.size(1) != ParameterIndexSize ){
//      All_U_Assembled.newsize(StencilSize, ParameterIndexSize);
//    }
//    for (cell=0; cell<StencilSize; ++cell){ //for each cell in the stencil
//      for (parameter=0; parameter < ParameterIndexSize; ++parameter){
//	All_U_Assembled(cell,parameter) = ( GeomWeightValue(iCell,jCell,cell)*
//					    (SolnBlk.*ReconstructedSoln)(i_index[cell],
//									 j_index[cell])[ParameterIndex[parameter]] );
//      }
//    } //endfor (cell)
//
//    // Step 3. Obtain solution to the linear equality-constrained least-square problem (i.e. first equation is a constraint)
//    // X_Assembled might have more columns than necessary, but that's not a problem.
//    Solve_Constrained_LS_Householder(A_Assembled,
//				     All_U_Assembled,
//				     X_Assembled,
//				     1);
//
//    // Step 4. Update the coefficients D (derivatives)
//    //**************************************************
//    for (parameter=0; parameter < ParameterIndexSize; ++parameter){
//      for (n=0; n<=CellTaylorDeriv(iCell,jCell).LastElem(); ++n){
//	CellTaylorDerivState(iCell,jCell,n)[ParameterIndex[parameter]] = X_Assembled(n,parameter);
//      }//endfor
//    }//endfor
//    
//  } else {
//    // Form both LHS and RHS in order to calculate the least-squares reconstruction
//
//    // SET VARIABLES USED ONLY IN THIS RECONSTRUCTION PROCESS
//
//    int krank;                        //< the final rank of matrix A is returned here
//    int IndexSumY, IndexSumX, P1, P2;
//    double CombP1X, CombP2Y;
//    double PowDistanceYC, PowDistanceXC;
//    double MaxWeight(0.0);
//    double IntSum(0.0);
//    int CurrentParameter;
//
//    // Ensure that enough memory is allocated for the current least-squares problem.
//    if (    A_Assembled.size(0) != (StencilSize-1)  ||     A_Assembled.size(1) != (NumberOfTaylorDerivatives()-1) ||
//	All_U_Assembled.size(0) != (StencilSize-1)  || All_U_Assembled.size(1) != ParameterIndexSize) {
//      // Resize the matrices accordingly
//      A_Assembled.newsize(StencilSize-1,NumberOfTaylorDerivatives()-1);
//      All_U_Assembled.newsize(StencilSize-1,ParameterIndexSize);
//    }
//
//    // START:   Set the LHS and RHS of the linear system 
//    // ***************************************************
//
//    // ==== Set the geometric weight associated with the reconstructed cell
//    GeometricWeights[0] = 1;
//
//    // Step 1. Compute the normalized geometric weights
//    for (cell=1; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil
//
//      /* Compute the X and Y component of the distance between
//	 the cell centers of the neighbour and the reconstructed cell */
//      DeltaCellCenters[cell] = CellCenter(i_index[cell],j_index[cell]) - CellCenter(iCell,jCell);
//    
//      /* Compute the geometric weight based on the centroid distance */
//      CENO_Geometric_Weighting(GeometricWeights[cell], DeltaCellCenters[cell].abs());
//
//      /* Compute the maximum geometric weight (this is used for normalization) */
//      MaxWeight = max(MaxWeight, GeometricWeights[cell]);
//    }
//
//    // Step 2. Set the approximate equations
//    for (cell=1 ; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil
//    
//      // compute the normalized geometric weight
//      GeometricWeights[cell] /= MaxWeight;
//
//      // *** SET the matrix A_Assembled of the linear system (LHS) ***
//      /* compute for each derivative the corresponding entry in the matrix of the linear system */
//      for (i=1; i<=CellTaylorDeriv(iCell,jCell).LastElem(); ++i){
//	// build the row of the matrix
//	P1 = CellTaylorDeriv(iCell,jCell,i).P1();  // identify P1
//	P2 = CellTaylorDeriv(iCell,jCell,i).P2();  // identify P2
//	A_Assembled(cell-1,i-1) = 0.0;  // set sumation variable to zero
//	CombP2Y = 1.0;        // the binomial coefficient "nC k" for k=0 is 1
//	PowDistanceYC = 1.0;  // initialize PowDistanceYC
//
//	// Compute geometric integral over the neighbour's domain
//	for (IndexSumY = 0; IndexSumY<=P2; ++IndexSumY){
//	  CombP1X = 1.0;       // the binomial coefficient "nC k" for k=0 is 1
//	  PowDistanceXC = 1.0; // initialize PowDistanceXC
//	  IntSum = 0.0;	     // reset internal sumation variable
//
//	  for (IndexSumX = 0; IndexSumX<=P1; ++IndexSumX){
//	    IntSum += ( CombP1X*PowDistanceXC*
//			Geom->CellGeomCoeffValue(i_index[cell],j_index[cell],P1-IndexSumX,P2-IndexSumY) );
//	    
//	    // update the binomial coefficients
//	    CombP1X = (P1-IndexSumX)*CombP1X/(IndexSumX+1); // the index is still the old one => expression for "nC k+1"
//	    PowDistanceXC *= DeltaCellCenters[cell].x;      // Update PowDistanceXC
//	  }//endfor
//
//	  A_Assembled(cell-1,i-1) += CombP2Y*PowDistanceYC*IntSum; // update the external sum
//
//	  CombP2Y = (P2-IndexSumY)*CombP2Y/(IndexSumY+1); // the index is still the old one => expression for "nC k+1"
//	  PowDistanceYC *= DeltaCellCenters[cell].y;    // Update PowDistanceYC
//	}//endfor
//
//	// subtract the corresponding geometric moment of cell (iCell,jCell) 
//	A_Assembled(cell-1,i-1) -= Geom->CellGeomCoeffValue(iCell,jCell,P1,P2);
//
//	// apply geometric weighting
//	A_Assembled(cell-1,i-1) *= GeometricWeights[cell];
//      }
//
//      // *** SET the matrix All_U_Assembled (Delta_U) of the linear system (RHS) ***
//      for (parameter = 0; parameter < ParameterIndexSize; ++parameter){
//	
//	// Solution parameter for the current setup
//	CurrentParameter = ParameterIndex[parameter];
//
//	// Compute Delta_U = U[neighbour] - U[cell] for each parameter
//	All_U_Assembled(cell-1,parameter) = ( (SolnBlk.*ReconstructedSoln)(i_index[cell],j_index[cell])[CurrentParameter] -
//					      (SolnBlk.*ReconstructedSoln)(iCell,jCell)[CurrentParameter] );
//	
//	// Apply geometric weighting
//	All_U_Assembled(cell-1,parameter) *= GeometricWeights[cell];
//
//      }	// endfor (parameter)
//      
//    }//endfor (cell)
//
//    // STOP:   Matrix A_Assembled of the linear system (LHS) built.
//    //         Matrix All_U_Assembled of the linear system (RHS) built.
//    // **********************************************************************
//
//    // Assign the average solution to D00 for each parameter in ParameterIndex
//    for (parameter = 0; parameter < ParameterIndexSize; ++parameter){
//      CellTaylorDeriv(iCell,jCell,0).D(ParameterIndex[parameter]) =
//	(SolnBlk.*ReconstructedSoln)(iCell,jCell)[ParameterIndex[parameter]];
//    }
//    
//    if (CENO_Execution_Mode::USE_LAPACK_LEAST_SQUARES) {
//      
//      // Solve the least-squares system with Lapack subroutine
//      /*********************************************************/
//      Solve_LS_Householder_F77(A_Assembled, All_U_Assembled,
//			       krank, ParameterIndexSize, StencilSize-1, NumberOfTaylorDerivatives()-1);
//
//      // Update the high-order derivatives
//      //***********************************
//      for (i = 1; i <= CellTaylorDeriv(iCell,jCell).LastElem(); ++i){
//
//	// Identify 'x' and 'y' powers of the i-th derivative
//	P1 = CellTaylorDeriv(iCell,jCell,i).P1();  // identify P1
//	P2 = CellTaylorDeriv(iCell,jCell,i).P2();  // identify P2
//
//	for (parameter = 0; parameter < ParameterIndexSize; ++parameter){
//	  // Set the i-th derivative for the current parameter
//	  CellTaylorDeriv(iCell,jCell,i).D(ParameterIndex[parameter]) = All_U_Assembled(i-1,parameter);
//	  
//	  // This equation ensures the mean conservation of the current parameter inside the reconstructed cell.
//	  CellTaylorDeriv(iCell,jCell,0).D(ParameterIndex[parameter]) -= ( Geom->CellGeomCoeffValue(iCell,jCell,P1,P2) 
//									   * All_U_Assembled(i-1,parameter) );
//	} // endfor (parameter)
//      }	// endfor (i)
//
//    } else { 
//
//      ColumnVector Rnorm(ParameterIndexSize);       //< residual norm of the LS problem for each parameter.
//      DenseMatrix Xm(NumberOfTaylorDerivatives()-1, ParameterIndexSize); //< storage for the solution to the least-square problem
//
//      /* Solve the overdetermined linear system of equations using the internal least-squares procedure */
//      /**************************************************************************************************/
//      Solve_LS_Householder(A_Assembled,All_U_Assembled,Xm,krank,Rnorm);
//
//      // Update the high-order derivatives
//      //***********************************
//      for (i = 1; i <= CellTaylorDeriv(iCell,jCell).LastElem(); ++i){
//
//	// Identify 'x' and 'y' powers of the i-th derivative
//	P1 = CellTaylorDeriv(iCell,jCell,i).P1();  // identify P1
//	P2 = CellTaylorDeriv(iCell,jCell,i).P2();  // identify P2
//
//	for (parameter = 0; parameter < ParameterIndexSize; ++parameter){
//	  // Set the i-th derivative for the current parameter
//	  CellTaylorDeriv(iCell,jCell,i).D(ParameterIndex[parameter]) = Xm(i-1,parameter);
//    
//	  // This equation ensures the mean conservation of the current parameter inside the reconstructed cell.
//	  CellTaylorDeriv(iCell,jCell,0).D(ParameterIndex[parameter]) -= ( Geom->CellGeomCoeffValue(iCell,jCell,P1,P2)
//									   * Xm(i-1,parameter) );
//
//	} // endfor (parameter)
//      } // endfor (i)
//
//    } // endif (CENO_Execution_Mode::USE_LAPACK_LEAST_SQUARES)
//
//  } // endif (IsPseudoInverseAllocated() && IsPseudoInversePreComputed())
//  
//}
//
//
///*!
// * Gather the data on each cell edge for imposing the required constraints 
// */
//template<class SOLN_STATE>
//template<class Soln_Block_Type> inline
//void HighOrder<SOLN_STATE>::FetchDataConstraints(Soln_Block_Type & SolnBlk,
//						   const int &iCell, const int &jCell,
//						   int & BC_Type,
//						   const int & parameter,
//						   const int & ConstrainedGQPs_West,  const int * ConstraintBCs_W,
//						   const int & ConstrainedGQPs_South, const int * ConstraintBCs_S,
//						   const int & ConstrainedGQPs_East,  const int * ConstraintBCs_E,
//						   const int & ConstrainedGQPs_North, const int * ConstraintBCs_N,
//						   Vector3DArray & Constraints_Loc,
//						   Vector3DArray & Constraints_Normals,
//						   BC_Type_Array & Constraints_BCs){
//  int n;
//
//  // === Ensure clean containers
//  Constraints_Loc.clear();
//  Constraints_Normals.clear();
//  Constraints_BCs.clear();
//  
//  //=== West edge ===
//  if ( ConstraintBCs_W[parameter] ){
//    // Constraints detected on the West face
//    // Fetch the data for imposing the constraints
//    if (jCell<JCl && Geom->ExtendSouth_BndWestSplineInfo != NULL){
//      Geom->ExtendSouth_BndWestSplineInfo[jCell].CopyGQPoints(Constraints_Loc);
//      Geom->ExtendSouth_BndWestSplineInfo[jCell].CopyNormalGQPoints(Constraints_Normals);
//    } else if ( jCell>=JCl && jCell<=JCu && Geom->BndWestSplineInfo != NULL){
//      Geom->BndWestSplineInfo[jCell].CopyGQPoints(Constraints_Loc);
//      Geom->BndWestSplineInfo[jCell].CopyNormalGQPoints(Constraints_Normals);
//    } else if ( jCell>JCu && Geom->ExtendNorth_BndWestSplineInfo != NULL){
//      Geom->ExtendNorth_BndWestSplineInfo[jCell-(JCu+1)].CopyGQPoints(Constraints_Loc);
//      Geom->ExtendNorth_BndWestSplineInfo[jCell-(JCu+1)].CopyNormalGQPoints(Constraints_Normals);
//    } else {
//      Geom->addGaussQuadPointsFaceW(iCell,jCell,Constraints_Loc,ConstrainedGQPs_West);
//      for (n = 0; n < ConstrainedGQPs_West; ++n){
//	Constraints_Normals.push_back(Geom->nfaceW(iCell,jCell));
//      }
//    }
//      
//    Constraints_BCs.push_back(&SolnBlk.BC_WestCell(jCell));
//
//    // store the type of the boundary condition for this cell
//    BC_Type = Geom->BCtypeW[jCell];
//  }
//
//
//  //=== South edge ===
//  if ( ConstraintBCs_S[parameter] ){
//    // Constraints detected on the South face
//    // Fetch the data for imposing the constraints
//    if (iCell<ICl && Geom->ExtendWest_BndSouthSplineInfo != NULL){
//      Geom->ExtendWest_BndSouthSplineInfo[iCell].CopyGQPoints(Constraints_Loc);
//      Geom->ExtendWest_BndSouthSplineInfo[iCell].CopyNormalGQPoints(Constraints_Normals);      
//    } else if (iCell>=ICl && iCell<=ICu && Geom->BndSouthSplineInfo != NULL){
//      Geom->BndSouthSplineInfo[iCell].CopyGQPoints(Constraints_Loc);
//      Geom->BndSouthSplineInfo[iCell].CopyNormalGQPoints(Constraints_Normals);      
//    } else if (iCell>ICu && Geom->ExtendEast_BndSouthSplineInfo != NULL){
//      Geom->ExtendEast_BndSouthSplineInfo[iCell-(ICu+1)].CopyGQPoints(Constraints_Loc);
//      Geom->ExtendEast_BndSouthSplineInfo[iCell-(ICu+1)].CopyNormalGQPoints(Constraints_Normals);      
//    } else {
//      Geom->addGaussQuadPointsFaceS(iCell,jCell,Constraints_Loc,ConstrainedGQPs_South);
//      for (n = 0; n < ConstrainedGQPs_South; ++n){
//	Constraints_Normals.push_back(Geom->nfaceS(iCell,jCell));
//      }
//    }
//    
//    Constraints_BCs.push_back(&SolnBlk.BC_SouthCell(iCell));
//
//    // store the type of the boundary condition for this cell
//    BC_Type = Geom->BCtypeS[iCell];
//  }
//
//  //=== East edge ===
//  if ( ConstraintBCs_E[parameter] ){
//    // Constraints detected on the East face
//    // Fetch the data for imposing the constraints
//    if (jCell<JCl && Geom->ExtendSouth_BndEastSplineInfo != NULL){
//      Geom->ExtendSouth_BndEastSplineInfo[jCell].CopyGQPoints(Constraints_Loc);
//      Geom->ExtendSouth_BndEastSplineInfo[jCell].CopyNormalGQPoints(Constraints_Normals);
//    } else if (jCell>=JCl && jCell<=JCu && Geom->BndEastSplineInfo != NULL){
//      Geom->BndEastSplineInfo[jCell].CopyGQPoints(Constraints_Loc);
//      Geom->BndEastSplineInfo[jCell].CopyNormalGQPoints(Constraints_Normals);
//    } else if (jCell>JCu && Geom->ExtendNorth_BndEastSplineInfo != NULL){
//      Geom->ExtendNorth_BndEastSplineInfo[jCell-(JCu+1)].CopyGQPoints(Constraints_Loc);
//      Geom->ExtendNorth_BndEastSplineInfo[jCell-(JCu+1)].CopyNormalGQPoints(Constraints_Normals);
//    } else {
//      Geom->addGaussQuadPointsFaceE(iCell,jCell,Constraints_Loc,ConstrainedGQPs_East);
//      for (n = 0; n < ConstrainedGQPs_East; ++n){
//	Constraints_Normals.push_back(Geom->nfaceE(iCell,jCell));
//      }
//    }
//      
//    Constraints_BCs.push_back(&SolnBlk.BC_EastCell(jCell));
//
//    // store the type of the boundary condition for this cell
//    BC_Type = Geom->BCtypeE[jCell];
//  }
//
//  //=== North edge ===
//  if ( ConstraintBCs_N[parameter] ){
//    // Constraints detected on the North face
//    // Fetch the data for imposing the constraints
//    if (iCell<ICl && Geom->ExtendWest_BndNorthSplineInfo != NULL){
//      Geom->ExtendWest_BndNorthSplineInfo[iCell].CopyGQPoints(Constraints_Loc);
//      Geom->ExtendWest_BndNorthSplineInfo[iCell].CopyNormalGQPoints(Constraints_Normals);      
//    } else if (iCell>=ICl && iCell<=ICu && Geom->BndNorthSplineInfo != NULL){
//      Geom->BndNorthSplineInfo[iCell].CopyGQPoints(Constraints_Loc);
//      Geom->BndNorthSplineInfo[iCell].CopyNormalGQPoints(Constraints_Normals);
//    } else if (iCell>ICu && Geom->ExtendEast_BndNorthSplineInfo != NULL){
//      Geom->ExtendEast_BndNorthSplineInfo[iCell-(ICu+1)].CopyGQPoints(Constraints_Loc);
//      Geom->ExtendEast_BndNorthSplineInfo[iCell-(ICu+1)].CopyNormalGQPoints(Constraints_Normals);
//    } else {
//      Geom->addGaussQuadPointsFaceN(iCell,jCell,Constraints_Loc,ConstrainedGQPs_North);
//      for (n = 0; n < ConstrainedGQPs_North; ++n){
//	Constraints_Normals.push_back(Geom->nfaceN(iCell,jCell));
//      }
//    }
//      
//    Constraints_BCs.push_back(&SolnBlk.BC_NorthCell(iCell));
//
//    // store the type of the boundary condition for this cell
//    BC_Type = Geom->BCtypeN[iCell];
//  }
//
//}

#endif
