/*!\file HighOrder2D_Reconstructions.h
  \brief Header file implementing the member functions of HighOrder2D class which reconstruct the solution.
  \note To use the functions defined in this file include 'HighOrder2D.h'!
*/

#ifndef _HIGHORDER_2D_RECONSTRUCTIONS_INCLUDED
#define _HIGHORDER_2D_RECONSTRUCTIONS_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
// None


/* -----------------------------------------------------------------
 * =============== BLOCK LEVEL 2D RECONSTRUCTIONS ==================
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
void HighOrder2D<SOLN_STATE>::ComputeUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk,
     								     const Soln_State & 
								     (Soln_Block_Type::*ReconstructedSoln)(const int &,
													   const int &) const){

  int i,j;

  /***************************************************************************
   *    Perform unconstrained unlimited high-order solution reconstruction   *
   **************************************************************************/

  // Set the freeze limiter flag. It's value affects the reset monotonicity data!
  _freeze_limiter = SolnBlk.Freeze_Limiter;

  // Carry out the solution reconstruction for cells in the specified range.
  for ( j  = StartJ ; j <= EndJ ; ++j ) {
    for ( i = StartI ; i <= EndI ; ++i ) {

      // Reset the monotonicity data
      ResetMonotonicityData(i,j);
      
      // Set the stencil of points used for reconstruction
      SetReconstructionStencil(i, j, i_index, j_index);

      // Compute the reconstruction for the current cell
      ComputeUnconstrainedUnlimitedSolutionReconstruction(SolnBlk, ReconstructedSoln,
							  i, j, i_index, j_index);
      
    } /* endfor */
  }/* endfor */

  // Check whether constrained reconstruction is required anywhere in the block.
  if ( !_constrained_block_reconstruction ){
    // No need to perform constrained reconstructions
    return;
  }


  /***************************************************************************
   *    Perform constrained unlimited high-order solution reconstruction     *
   **************************************************************************/

  // Check WEST boundary
  if (_constrained_WEST_reconstruction){
    // Add constrained reconstruction here
    for (j = StartJ_ConstrWest; j <= EndJ_ConstrWest; ++j){
      for (i = StartI_ConstrWest; i <= EndI_ConstrWest; ++i){
	
	// Reset the monotonicity data
	ResetMonotonicityData(i,j);
	
	// Set the biased stencil of points used for reconstruction
	SetConstrainedReconstructionStencil(i, j, i_index_ave, j_index_ave);
	
	// Compute the constrained reconstruction for the current cell
	ComputeConstrainedUnlimitedSolutionReconstruction(SolnBlk, ReconstructedSoln,
							  i, j, i_index_ave, j_index_ave);

      }	// endfor
    }// endfor
  } 

  // Check EAST boundary
  if (_constrained_EAST_reconstruction){
    // Add constrained reconstruction here
    for (j = StartJ_ConstrEast; j <= EndJ_ConstrEast; ++j){
      for (i = StartI_ConstrEast; i <= EndI_ConstrEast; ++i){
	
	// Reset the monotonicity data
	ResetMonotonicityData(i,j);
	
	// Set the biased stencil of points used for reconstruction
	SetConstrainedReconstructionStencil(i, j, i_index_ave, j_index_ave);
	
	// Compute the constrained reconstruction for the current cell
	ComputeConstrainedUnlimitedSolutionReconstruction(SolnBlk, ReconstructedSoln,
							  i, j, i_index_ave, j_index_ave);

      }	// endfor
    }// endfor
  } 

  // Check NORTH boundary
  if (_constrained_NORTH_reconstruction){
    // Add constrained reconstruction here
    for (j = StartJ_ConstrNorth; j <= EndJ_ConstrNorth; ++j){
      for (i = StartI_ConstrNorth; i <= EndI_ConstrNorth; ++i){
	
	// Reset the monotonicity data
	ResetMonotonicityData(i,j);
	
	// Set the biased stencil of points used for reconstruction
	SetConstrainedReconstructionStencil(i, j, i_index_ave, j_index_ave);
	
	// Compute the constrained reconstruction for the current cell
	ComputeConstrainedUnlimitedSolutionReconstruction(SolnBlk, ReconstructedSoln,
							  i, j, i_index_ave, j_index_ave);

      }	// endfor
    }// endfor
  } 

  // Check SOUTH boundary
  if (_constrained_SOUTH_reconstruction){
    // Add constrained reconstruction here
    for (j = StartJ_ConstrSouth; j <= EndJ_ConstrSouth; ++j){
      for (i = StartI_ConstrSouth; i <= EndI_ConstrSouth; ++i){
	
	// Reset the monotonicity data
	ResetMonotonicityData(i,j);
	
	// Set the biased stencil of points used for reconstruction
	SetConstrainedReconstructionStencil(i, j, i_index_ave, j_index_ave);
	
	// Compute the constrained reconstruction for the current cell
	ComputeConstrainedUnlimitedSolutionReconstruction(SolnBlk, ReconstructedSoln,
							  i, j, i_index_ave, j_index_ave);

      }	// endfor
    }// endfor
  }

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
void HighOrder2D<SOLN_STATE>::ComputeReconstructionPseudoInverse(void){

  int i,j;

  // == Check if the pseudo-inverse has been allocated and it hasn't been precomputed
  if ( IsPseudoInverseAllocated() && !IsPseudoInversePreComputed() ){

    // == Check if the reconstruction polynomial is piecewise constant
    if (RecOrder() == 0){
      // There is no need to calculate pseudo-inverse
      // Confirm the pseudo-inverse calculation
      _calculated_psinv = true;
      return;
    }

    // Calculate the pseudo-inverse using the central stencil for cells in the specified range
    for ( j  = StartJ ; j <= EndJ ; ++j ) {
      for ( i = StartI ; i <= EndI ; ++i ) {
	
	// Set the stencil of points used for reconstruction
	SetReconstructionStencil(i, j, i_index, j_index);

	// Compute the pseudo-inverse for the current cell
	ComputeCellReconstructionPseudoInverse(i, j, i_index, j_index);
      }/* endfor */
    }/* endfor */

    
    // Check whether constrained reconstruction is required anywhere in the block.
    if ( !_constrained_block_reconstruction ){
      // No need to calculate anything for the constrained reconstruction
      
      // Confirm the pseudo-inverse calculation
      _calculated_psinv = true;
      
      return;
    }


    /********************************************************************************************
     * Calculate the matrices that need to be stored in order to speed up the high-order solution
     * reconstruction of the cells affected by the presence of constrained boundaries
     * (i.e. calculate the pseudo-inverse of the unconstrained unlimited reconstruction
     * for those cells that don't have constraints but use a biased supporting stencil,
     * calculate the LHS matrix associate with the least-squares problem and the mean conservation
     * for those cells that have also constraint equations.
     ********************************************************************************************/
    
    // Check WEST boundary
    if (_constrained_WEST_reconstruction){
      // Calculate pseudo-inverse or part of LHS assemble matrix here
      for (j = StartJ_ConstrWest; j <= EndJ_ConstrWest; ++j){
	for (i = StartI_ConstrWest; i <= EndI_ConstrWest; ++i){
	  
	  // Get matrix for the current cell (i.e. pseudo-inverse or LHS)
	  ComputeCellReconstructionPseudoInverseNearConstrainedBoundaries(WEST,i,j);

	} // endfor
      }// endfor
    } 

    // Check EAST boundary
    if (_constrained_EAST_reconstruction){
      // Calculate pseudo-inverse or part of LHS assemble matrix here
      for (j = StartJ_ConstrEast; j <= EndJ_ConstrEast; ++j){
	for (i = StartI_ConstrEast; i <= EndI_ConstrEast; ++i){

	  // Get matrix for the current cell (i.e. pseudo-inverse or LHS)
	  ComputeCellReconstructionPseudoInverseNearConstrainedBoundaries(EAST,i,j);

	} // endfor
      }// endfor
    } 

    // Check NORTH boundary
    if (_constrained_NORTH_reconstruction){
      // Calculate pseudo-inverse or part of LHS assemble matrix here
      for (j = StartJ_ConstrNorth; j <= EndJ_ConstrNorth; ++j){
	for (i = StartI_ConstrNorth; i <= EndI_ConstrNorth; ++i){

	  // Get matrix for the current cell (i.e. pseudo-inverse or LHS)
	  ComputeCellReconstructionPseudoInverseNearConstrainedBoundaries(NORTH,i,j);
	
	} // endfor
      }// endfor
    } 

    // Check SOUTH boundary
    if (_constrained_SOUTH_reconstruction){
      // Calculate pseudo-inverse or part of LHS assemble matrix here
      for (j = StartJ_ConstrSouth; j <= EndJ_ConstrSouth; ++j){
	for (i = StartI_ConstrSouth; i <= EndI_ConstrSouth; ++i){

	  // Get matrix for the current cell (i.e. pseudo-inverse or LHS)
	  ComputeCellReconstructionPseudoInverseNearConstrainedBoundaries(SOUTH,i,j);
	
	} // endfor
      }// endfor
    }
    
    // Confirm the pseudo-inverse calculation
    _calculated_psinv = true;

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
void HighOrder2D<SOLN_STATE>::EnforceMonotonicityToNonSmoothInterpolants(Soln_Block_Type &SolnBlk,
									 const int &Limiter,
									 const Soln_State &
									 (Soln_Block_Type::*ReconstructedSoln)(const int &,
													       const int &) const){

  // Set local variables
  int i,j;
  
  // == Check if the reconstruction polynomial is piecewise constant
  if (RecOrder() == 0){
    // There is no need to enforce any monotonicity
    return;
  }

  if (CENO_Execution_Mode::CENO_DROP_ORDER){
    // Carry on actions required to enforce monotonicity

    // Switch to limited piecewise linear reconstruction those interpolants detected as non-smooth
    for ( j  = StartJ_LPWL ; j <= EndJ_LPWL ; ++j ) {
      for ( i = StartI_LPWL ; i <= EndI_LPWL ; ++i ) {
      
	if ( IsThereAnyNonSmoothHighOrderReconstruction(i,j) ){
	  // One or more solution variables need to have the interpolant switched to a limited piecewise linear one.
	  ComputeLimitedPiecewiseLinearSolutionReconstruction(SolnBlk,
							      i,j,
							      Limiter,
							      ReconstructedSoln);
	} // endif
 
      } /* endfor(i) */
    } /* endfor(j) */

    // Check whether reconstruction based flux calculation is required anywhere in the block.
    if ( !_constrained_block_reconstruction ){
      // All cells involved in flux calculation have been checked for non-smooth interpolants.
      // No need to do anything more.
      return;
    }


    /* Motivation of the algorithm below:
       If reconstruction based flux calculation is required at some of the boundaries
       and non-smooth solution interpolants are detected near these boundaries the 
       flux is not going to be computed based on the high-order interpolant but on
       solving a Riemann problem at the interface.
       The purpose of the algorithm below is to ensure that a limited piecewise linear
       reconstruction is available in the first ghost cells that have interface with an
       interior cell detected with inadequate interpolant. Thus, when the flux calculation 
       is performed, the Riemann problem for those interfaces can be solved.
       Note that no interior cells are going to be affected by the code that follows!
       Note also that trying to obtain a high-order interpolant in these ghost cells is not
       justified based on accuracy and computational efficiency reasons.
    */

    // Check WEST boundary
    if (_constrained_WEST_reconstruction){
      for (j = JCl; j <= JCu; ++j){
	if ( IsThereAnyNonSmoothHighOrderReconstruction(ICl,j) ) { // check the interior cell
	  // flag all reconstructions of the adjacent ghost cell as non-smooth
	  FlagCellReconstructionsAsNonSmooth(ICl-1,j);
	  // perform a limited piecewise linear reconstruction
	  ComputeLimitedPiecewiseLinearSolutionReconstruction(SolnBlk,
							      ICl-1,j,
							      Limiter,
							      ReconstructedSoln);
	}// endif
      }// enfor 
    }// endif 

    // Check EAST boundary
    if (_constrained_EAST_reconstruction){
      for (j = JCl; j <= JCu; ++j){
	if ( IsThereAnyNonSmoothHighOrderReconstruction(ICu,j) ) { // check the interior cell
	  // flag all reconstructions of the adjacent ghost cell as non-smooth
	  FlagCellReconstructionsAsNonSmooth(ICu+1,j);
	  // perform a limited piecewise linear reconstruction
	  ComputeLimitedPiecewiseLinearSolutionReconstruction(SolnBlk,
							      ICu+1,j,
							      Limiter,
							      ReconstructedSoln);
	}// endif
      }// enfor 
    }// endif 

    // Check NORTH boundary
    if (_constrained_NORTH_reconstruction){
      for (i = ICl; i <= ICu; ++i){
	if ( IsThereAnyNonSmoothHighOrderReconstruction(i,JCu) ) { // check the interior cell
	  // flag all reconstructions of the adjacent ghost cell as non-smooth
	  FlagCellReconstructionsAsNonSmooth(i,JCu+1);
	  // perform a limited piecewise linear reconstruction
	  ComputeLimitedPiecewiseLinearSolutionReconstruction(SolnBlk,
							      i,JCu+1,
							      Limiter,
							      ReconstructedSoln);

	}// endif
      }// enfor 
    }// endif 

    // Check SOUTH boundary
    if (_constrained_SOUTH_reconstruction){
      for (i = ICl; i <= ICu; ++i){
	if ( IsThereAnyNonSmoothHighOrderReconstruction(i,JCl) ) { // check the interior cell
	  // flag all reconstructions of the adjacent ghost cell as non-smooth
	  FlagCellReconstructionsAsNonSmooth(i,JCl-1);
	  // perform a limited piecewise linear reconstruction
	  ComputeLimitedPiecewiseLinearSolutionReconstruction(SolnBlk,
							      i,JCl-1,
							      Limiter,
							      ReconstructedSoln);

	}// endif
      }// enfor 
    }// endif 

  } else {
    /* Reset monotonicity flags for boundary cells near splines which require reconstruction based flux calculation.
       If Riemann based flux calculation is desired there is nothing to be done.
     */

    // Check whether reset of monotonicity flags is required anywhere in the block.
    if ( !_constrained_block_reconstruction ){
      // No need to reset
      return;
    }

    // Check WEST boundary
    if (_constrained_WEST_reconstruction){
      for (j = JCl - 1; j <= JCu + 1; ++j){
	// reset the monotonicity flag for the interior cell
	ResetMonotonicityData(ICl,j);
      }
    } 

    // Check EAST boundary
    if (_constrained_EAST_reconstruction){
      for (j = JCl - 1; j <= JCu + 1; ++j){
	// reset the monotonicity flag for the interior cell
	ResetMonotonicityData(ICu,j);
      }
    } 

    // Check NORTH boundary
    if (_constrained_NORTH_reconstruction){
      for (i = ICl - 1; i <= ICu + 1; ++i){
	// reset the monotonicity flag for the interior cell
	ResetMonotonicityData(i,JCu);
      }
    } 

    // Check SOUTH boundary
    if (_constrained_SOUTH_reconstruction){
      for (i = ICl - 1; i <= ICu + 1; ++i){
	// reset the monotonicity flag for the interior cell
	ResetMonotonicityData(i,JCl);
      }
    }

  } // endif(CENO_Execution_Mode::CENO_DROP_ORDER)
  
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
void HighOrder2D<SOLN_STATE>::ComputeHighOrderSolutionReconstruction(Soln_Block_Type &SolnBlk,
								     const int &Limiter,
								     const Soln_State & 
								     (Soln_Block_Type::*ReconstructedSoln)(const int &,
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
 * ================ CELL LEVEL 2D RECONSTRUCTIONS ==================
 * ----------------------------------------------------------------*/

/*! 
 * Compute the unlimited k-exact high-order reconstruction
 * proposed by Barth (1993) for a specified computational cell.
 * This reconstruction doesn't account for any constrain
 * other than the mean quantity conservation.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type>
void HighOrder2D<SOLN_STATE>::
ComputeUnconstrainedUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk,
						    const Soln_State & 
						    (Soln_Block_Type::*ReconstructedSoln)(const int &,const int &) const,
						    const int &iCell, const int &jCell,
						    const IndexType & i_index,
						    const IndexType & j_index){

  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS

  int StencilSize(i_index.size()); 
  int cell, i, parameter;
  int P1, P2;

  // *********  Assign the average solution to D00 ***********
  CellTaylorDerivState(iCell,jCell,0) = (SolnBlk.*ReconstructedSoln)(iCell,jCell);

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
	Delta_U(cell-1) = ( (SolnBlk.*ReconstructedSoln)(i_index[cell],j_index[cell])[parameter] -
			    (SolnBlk.*ReconstructedSoln)(iCell,jCell)[parameter] );

	// Apply the precomputed geometric weight to the Delta_U term
	Delta_U(cell-1) *= GeomWeightValue(iCell,jCell,cell);
      } 
     
      // Step 2. Find the solution of the linear-system for the current parameter
      X = Cell_LHS_Inv(iCell,jCell) * Delta_U;
      
      // Step 3. Update the high-order derivatives for the current parameter
      for (i = 1; i <= CellTaylorDeriv(iCell,jCell).LastElem(); ++i){

	// Identify 'x' and 'y' powers of the i-th derivative
	P1 = CellTaylorDeriv(iCell,jCell,i).P1();  // identify P1
	P2 = CellTaylorDeriv(iCell,jCell,i).P2();  // identify P2

	// Set the i-th derivative for the current parameter
	CellTaylorDeriv(iCell,jCell,i).D(parameter) = X(i-1);
	
	// This equation ensures the mean conservation of the current parameter inside the reconstructed cell.
	CellTaylorDeriv(iCell,jCell,0).D(parameter) -= Geom->CellGeomCoeffValue(iCell,jCell,P1,P2) * X(i-1);
      }
    }
    // STOP: Reconstruction solution (i.e. derivatives) obtained.
    // ************************************************************

  } else {
    
    // Form both LHS and RHS in order to calculate the least-squares reconstruction

    // SET VARIABLES USED ONLY IN THIS RECONSTRUCTION PROCESS

    int krank;                        //< the final rank of matrix A is returned here
    int IndexSumY, IndexSumX, P1, P2;
    double CombP1X, CombP2Y;
    double PowDistanceYC, PowDistanceXC;
    double MaxWeight(0.0);
    double IntSum(0.0);

    // Ensure that enough memory is allocated for the current least-squares problem.
    if (A.size(0) != StencilSize) {
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

      /* Compute the X and Y component of the distance between
	 the cell centers of the neighbour and the reconstructed cell */
      DeltaCellCenters[cell] = CellCenter(i_index[cell],j_index[cell]) - CellCenter(iCell,jCell);
    
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
      for (i=1; i<=CellTaylorDeriv(iCell,jCell).LastElem(); ++i){
	// build the row of the matrix
	P1 = CellTaylorDeriv(iCell,jCell,i).P1();  // identify P1
	P2 = CellTaylorDeriv(iCell,jCell,i).P2();  // identify P2
	A(cell-1,i-1) = 0.0;  // set sumation variable to zero
	CombP2Y = 1.0;        // the binomial coefficient "nC k" for k=0 is 1
	PowDistanceYC = 1.0;  // initialize PowDistanceYC

	// Compute geometric integral over the neighbour's domain
	for (IndexSumY = 0; IndexSumY<=P2; ++IndexSumY){
	  CombP1X = 1.0;       // the binomial coefficient "nC k" for k=0 is 1
	  PowDistanceXC = 1.0; // initialize PowDistanceXC
	  IntSum = 0.0;	     // reset internal sumation variable

	  for (IndexSumX = 0; IndexSumX<=P1; ++IndexSumX){
	    IntSum += ( CombP1X*PowDistanceXC*
			Geom->CellGeomCoeffValue(i_index[cell],j_index[cell],P1-IndexSumX,P2-IndexSumY) );
	    
	    // update the binomial coefficients
	    CombP1X = (P1-IndexSumX)*CombP1X/(IndexSumX+1); // the index is still the old one => expression for "nC k+1"
	    PowDistanceXC *= DeltaCellCenters[cell].x;      // Update PowDistanceXC
	  }//endfor

	  A(cell-1,i-1) += CombP2Y*PowDistanceYC*IntSum; // update the external sum

	  CombP2Y = (P2-IndexSumY)*CombP2Y/(IndexSumY+1); // the index is still the old one => expression for "nC k+1"
	  PowDistanceYC *= DeltaCellCenters[cell].y;    // Update PowDistanceYC
	}//endfor

	// subtract the corresponding geometric moment of cell (iCell,jCell) 
	A(cell-1,i-1) -= Geom->CellGeomCoeffValue(iCell,jCell,P1,P2);

	// apply geometric weighting
	A(cell-1,i-1) *= GeometricWeights[cell];
      }

      // *** SET the matrix All_Delta_U of the linear system (RHS) ***
      for (parameter = 1; parameter <= NumberOfVariables(); ++parameter){

	// Compute Delta_U = U[neighbour] - U[cell] for each parameter
	All_Delta_U(cell-1,parameter-1) = ( (SolnBlk.*ReconstructedSoln)(i_index[cell],j_index[cell])[parameter] -
					    (SolnBlk.*ReconstructedSoln)(iCell,jCell)[parameter] );
	
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
      for (i = 1; i <= CellTaylorDeriv(iCell,jCell).LastElem(); ++i){

	// Identify 'x' and 'y' powers of the i-th derivative
	P1 = CellTaylorDeriv(iCell,jCell,i).P1();  // identify P1
	P2 = CellTaylorDeriv(iCell,jCell,i).P2();  // identify P2

	for (parameter = 1; parameter <= NumberOfVariables(); ++parameter){
	  // Set the i-th derivative for the current parameter
	  CellTaylorDeriv(iCell,jCell,i).D(parameter) = All_Delta_U(i-1,parameter-1);
	  
	  // This equation ensures the mean conservation of the current parameter inside the reconstructed cell.
	  CellTaylorDeriv(iCell,jCell,0).D(parameter) -= Geom->CellGeomCoeffValue(iCell,jCell,P1,P2) * All_Delta_U(i-1,
														   parameter-1);
	} // endfor (parameter)
      }	// endfor (i)

    } else { 

      ColumnVector Rnorm(NumberOfVariables());       //< store the residual norm of the LS problem for each parameter.
      DenseMatrix Xm(NumberOfTaylorDerivatives()-1, NumberOfVariables()); //< store the solution to the least-square problem

      /* Solve the overdetermined linear system of equations using the internal least-squares procedure */
      /**************************************************************************************************/
      Solve_LS_Householder(A,All_Delta_U,Xm,krank,Rnorm);

      // Update the high-order derivatives
      //***********************************
      for (i = 1; i <= CellTaylorDeriv(iCell,jCell).LastElem(); ++i){

	// Identify 'x' and 'y' powers of the i-th derivative
	P1 = CellTaylorDeriv(iCell,jCell,i).P1();  // identify P1
	P2 = CellTaylorDeriv(iCell,jCell,i).P2();  // identify P2

	for (parameter = 1; parameter <= NumberOfVariables(); ++parameter){
	  // Set the i-th derivative for the current parameter
	  CellTaylorDeriv(iCell,jCell,i).D(parameter) = Xm(i-1,parameter-1);
    
	  // This equation ensures the mean conservation of the current parameter inside the reconstructed cell.
	  CellTaylorDeriv(iCell,jCell,0).D(parameter) -= Geom->CellGeomCoeffValue(iCell,jCell,P1,P2) * Xm(i-1,parameter-1);

	} // endfor (parameter)
      } // endfor (i)

    } // endif (CENO_Execution_Mode::USE_LAPACK_LEAST_SQUARES)

  } // endif (IsPseudoInverseAllocated() && IsPseudoInversePreComputed())

}

/*! 
 * Compute the pseudo-inverse of the left-hand-side term in the 
 * unlimited k-exact high-order reconstruction for a specified
 * computational cell based on the information provided by the
 * associated grid.
 */
template<class SOLN_STATE>
void HighOrder2D<SOLN_STATE>::ComputeCellReconstructionPseudoInverse(const int &iCell, const int &jCell,
								     const IndexType & i_index,
								     const IndexType & j_index){
  
  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS

  int StencilSize(i_index.size());  
  int IndexSumY, IndexSumX, P1, P2;
  double CombP1X, CombP2Y;
  double PowDistanceYC, PowDistanceXC;
  int cell, i;
  double MaxWeight(0.0);
  double IntSum(0.0);

  // Ensure that the LHS matrix is formated correctly.
  // Memory shouldn't be allocated here, only the dimensions should be defined properly.
  Cell_LHS_Inv(iCell,jCell).newsize(StencilSize - 1, NumberOfTaylorDerivatives() - 1);
  GeomWeights(iCell,jCell).resize(StencilSize);


  // START:   Set the LHS of the linear system 
  // ***************************************************

  // ==== Set the geometric weight associated with the reconstructed cell
  GeomWeightValue(iCell,jCell,0) = 1;

  // Step1. Compute the normalized geometric weights
  for (cell=1; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil

    /* Compute the X and Y component of the distance between
       the cell centers of the neighbour and the reconstructed cell */
    DeltaCellCenters[cell] = CellCenter(i_index[cell],j_index[cell]) - CellCenter(iCell,jCell);
    
    /* Compute the geometric weight based on the centroid distance */
    CENO_Geometric_Weighting(GeomWeightValue(iCell,jCell,cell), DeltaCellCenters[cell].abs());

    /* Compute the maximum geometric weight (this is used for normalization) */
    MaxWeight = max(MaxWeight, GeomWeightValue(iCell,jCell,cell));
  }

  // Step2. Set the approximate equations
  for (cell=1 ; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil
    
    // compute the normalized geometric weight
    GeomWeightValue(iCell,jCell,cell) /= MaxWeight;

    // *** SET the matrix of the linear system (LHS) ***
    /* compute for each derivative the corresponding entry in the matrix of the linear system */
    for (i=1; i<=CellTaylorDeriv(iCell,jCell).LastElem(); ++i){
      // build the row of the matrix
      P1 = CellTaylorDeriv(iCell,jCell,i).P1();  // identify P1
      P2 = CellTaylorDeriv(iCell,jCell,i).P2();  // identify P2
      Cell_LHS_Inv_Value(iCell,jCell,cell-1,i-1) = 0.0;  // set sumation variable to zero
      CombP2Y = 1.0;        // the binomial coefficient "nC k" for k=0 is 1
      PowDistanceYC = 1.0;  // initialize PowDistanceYC

      // Compute geometric integral over the neighbour's domain
      for (IndexSumY = 0; IndexSumY<=P2; ++IndexSumY){
	CombP1X = 1.0;       // the binomial coefficient "nC k" for k=0 is 1
	PowDistanceXC = 1.0; // initialize PowDistanceXC
	IntSum = 0.0;	     // reset internal sumation variable

	for (IndexSumX = 0; IndexSumX<=P1; ++IndexSumX){
	  IntSum += ( CombP1X*PowDistanceXC*
		      Geom->CellGeomCoeffValue(i_index[cell],j_index[cell],P1-IndexSumX,P2-IndexSumY) );
	    
	  // update the binomial coefficients
	  CombP1X = (P1-IndexSumX)*CombP1X/(IndexSumX+1); // the index is still the old one => expression for "nC k+1"
	  PowDistanceXC *= DeltaCellCenters[cell].x;      // Update PowDistanceXC
	}//endfor

	Cell_LHS_Inv_Value(iCell,jCell,cell-1,i-1) += CombP2Y*PowDistanceYC*IntSum; // update the external sum

	CombP2Y = (P2-IndexSumY)*CombP2Y/(IndexSumY+1); // the index is still the old one => expression for "nC k+1"
	PowDistanceYC *= DeltaCellCenters[cell].y;    // Update PowDistanceYC
      }//endfor

      // subtract the corresponding geometric moment of cell (iCell,jCell) 
      Cell_LHS_Inv_Value(iCell,jCell,cell-1,i-1) -= Geom->CellGeomCoeffValue(iCell,jCell,P1,P2);

      // apply geometric weighting
      Cell_LHS_Inv_Value(iCell,jCell,cell-1,i-1) *= GeomWeightValue(iCell,jCell,cell);
    }
      
  }//endfor (cell)

  // STOP:   Matrix of the linear system (LHS) built. 
  //         For kExact_Reconstruction away from some special curved boundaries 
  //         the same matrix is used for all variables (same geometry) and  
  //         at every time step as long as the mesh is the same.
  // **********************************************************************

  // Compute the pseudo-inverse and override the LHS term.
  // This operation will change the dimensions of the matrix.
  Cell_LHS_Inv(iCell,jCell).pseudo_inverse_override();
  
}

/*! 
 * Compute the pseudo-inverse of the left-hand-side term in the 
 * unlimited k-exact high-order reconstruction for a specified
 * computational cell based on the information provided by the
 * associated grid. 
 * The reconstruction of this cell is influence by the presence
 * of curved boundaries.
 * If the pseudo-inverse cannot be computed because boundary condition
 * constraints must be added to the linear system, the LHS matrix of the 
 * k-exact least-squares reconstruction is stored instead.
 */
template<class SOLN_STATE>
void HighOrder2D<SOLN_STATE>::ComputeCellReconstructionPseudoInverseNearConstrainedBoundaries(const int& BOUNDARY,
											      const int &iCell,
											      const int &jCell){

  int constrGQP;
  string ErrorMsg;

  // Set the biased stencil of points used for the reconstruction of the current cell
  SetConstrainedReconstructionStencil(iCell, jCell, i_index_ave, j_index_ave);

  // Determine the number of constraints for the current cell and set the error message
  switch(BOUNDARY){
  case NORTH:
    constrGQP = Geom->NumOfConstrainedGaussQuadPoints_North(iCell,jCell);
    
    ErrorMsg = "HighOrder2D<SOLN_STATE>::ComputeCellReconstructionPseudoInverseNearConstrainedBoundaries() ERROR! The pseudo-inverse couldn't be computed for a cell affected by the North boundary.";
    break;

  case SOUTH:
    constrGQP = Geom->NumOfConstrainedGaussQuadPoints_South(iCell,jCell);

    ErrorMsg = "HighOrder2D<SOLN_STATE>::ComputeCellReconstructionPseudoInverseNearConstrainedBoundaries() ERROR! The pseudo-inverse couldn't be computed for a cell affected by the South boundary.";
    break;

  case EAST:
    constrGQP = Geom->NumOfConstrainedGaussQuadPoints_East(iCell,jCell);

    ErrorMsg = "HighOrder2D<SOLN_STATE>::ComputeCellReconstructionPseudoInverseNearConstrainedBoundaries() ERROR! The pseudo-inverse couldn't be computed for a cell affected by the East boundary.";
    break;

  case WEST:
    constrGQP = Geom->NumOfConstrainedGaussQuadPoints_West(iCell,jCell);

    ErrorMsg = "HighOrder2D<SOLN_STATE>::ComputeCellReconstructionPseudoInverseNearConstrainedBoundaries() ERROR! The pseudo-inverse couldn't be computed for a cell affected by the West boundary.";
    break;
  }

  if (constrGQP == 0){
    /* There are NO constraints for this cell,
       so a pseudo-inverse with the biased stencil may be computed. */

    // Check overdeterminancy
    if ( i_index_ave.size() > NumberOfTaylorDerivatives() ){
      // The pseudo-inverse CAN be computed

      // Compute the pseudo-inverse for the current cell
      ComputeCellReconstructionPseudoInverse(iCell, jCell, i_index_ave, j_index_ave);
	      
    } else {
      // The pseudo-inverse CANNOT be computed
      // Throw an error
      throw runtime_error(ErrorMsg);
    }

  } else {
    /* There are constraints for this cell, so the pseudo-inverse CANNOT be computed.
       Instead, the least-squares part of the reconstruction matrix will be stored.
       That is, the mean conservation equation in the reconstructed cell (first equation)
       and the approximate equations for the neighbouring cells. */
	    
    /********* Generate the exact and approximate mean conservation equations ***********/
    /************************************************************************************/
    Set_LHS_MeanValueConservation_Equations(iCell,jCell,
					    i_index_ave, j_index_ave,
					    Cell_LHS_Inv(iCell,jCell),
					    GeomWeights(iCell,jCell));
  } // endif (constrGQP == 0)
  
}

/*! 
 * Performs the reconstruction of a limited piecewise 
 * linear solution state within a given cell (iCell,jCell) of   
 * the computational mesh for the specified             
 * quadrilateral solution block.  A least squares       
 * approach is used in the evaluation of the unlimited  
 * solution gradients.  Several slope limiters may be   
 * used.
 * The high-order derivatives of the variable flag as unfit
 * are replaced with the first-order ones computed with
 * this subroutine. 
 */
template<class SOLN_STATE>
template<class Soln_Block_Type>
void HighOrder2D<SOLN_STATE>::
ComputeLimitedPiecewiseLinearSolutionReconstruction(Soln_Block_Type &SolnBlk,
						    const int &iCell, const int &jCell,
						    const int &Limiter,
						    const Soln_State & 
						    (Soln_Block_Type::*ReconstructedSoln)(const int &,const int &) const){

  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS

  int n, parameter, n_pts;
  double MaxWeight(0.0);
  double DxDx_ave(0), DxDy_ave(0), DyDy_ave(0);
  Soln_State DU, DUDx_ave(0), DUDy_ave(0);
  double u0Min, u0Max;
  double *uQuad(NULL);
  Vector2D *GQP(NULL);
  int NumGQP, GQP_North, GQP_South, GQP_East, GQP_West;
  bool faceNorth, faceSouth, faceEast, faceWest;

  // == Check if the reconstruction polynomial is piecewise constant
  if (RecOrder() == 0){
    // There is no need to calculate the reconstruction
    return;
  }

  /* Carry out the limited solution reconstruction in
     each cell of the computational mesh. */

  /* Determine the number of neighbouring cells to
     be used in the reconstruction procedure.
     This stencil might be different near boundaries
     and it is influenced by the boundary conditions. */
  SolnBlk.SetPiecewiseLinearReconstructionStencil(iCell,jCell,
						  I_Index,J_Index,
						  n_pts);

  // Perform reconstruction.
  if (n_pts > 0) {
    
    // Perform piecewise linear reconstruction only if the reconstruction order is greater than 1
    if (RecOrder() > 1){
      
      // Compute distance between centroids and the geometric weights
      for ( n = 0 ; n < n_pts ; ++n ) {
	/* Compute the X and Y component of the distance between
	   the cell centers of the neighbour and the reconstructed cell */
	dX[n] = Geom->Cell[ I_Index[n] ][ J_Index[n] ].Xc - Geom->Cell[iCell][jCell].Xc;

	/* Compute the geometric weight based on the centroid distance */
	CENO_Geometric_Weighting(geom_weights[n], dX[n].abs());

	/* Compute the maximum geometric weight (this is used for normalization) */
	MaxWeight = max(MaxWeight, geom_weights[n]);
      }

      for ( n = 0 ; n < n_pts ; ++n ) {
	// compute the normalized geometric weight
	geom_weights[n] /= MaxWeight;

	// compute the square of the normalized geometric weight
	geom_weights[n] *= geom_weights[n];

	DU = (SolnBlk.*ReconstructedSoln)( I_Index[n], J_Index[n] ) - (SolnBlk.*ReconstructedSoln)(iCell,jCell);
	DUDx_ave += DU*(geom_weights[n]*dX[n].x);
	DUDy_ave += DU*(geom_weights[n]*dX[n].y);
	DxDx_ave += geom_weights[n]*dX[n].x*dX[n].x;
	DxDy_ave += geom_weights[n]*dX[n].x*dX[n].y;
	DyDy_ave += geom_weights[n]*dX[n].y*dX[n].y;
      } /* endfor */
    					    
      DUDx_ave /= double(n_pts);
      DUDy_ave /= double(n_pts);
      DxDx_ave /= double(n_pts);
      DxDy_ave /= double(n_pts);
      DyDy_ave /= double(n_pts);

      // Calculate the first-order derivatives
      dUdx = ( (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
	       (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
      dUdy = ( (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
	       (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );

    }//endif(RecOrder())

    // Calculate slope limiters or used the frozen ones.
    if (!SolnBlk.Freeze_Limiter) {

      // Calculate the new slope limiters.

      // Get number of flux calculation points and edge type for each of the 4 cell faces.
      GQP_North = Geom->NumOfFluxCalculationGaussQuadPoints_North(iCell,jCell,faceNorth);
      GQP_South = Geom->NumOfFluxCalculationGaussQuadPoints_South(iCell,jCell,faceSouth);
      GQP_East  = Geom->NumOfFluxCalculationGaussQuadPoints_East(iCell,jCell,faceEast);
      GQP_West  = Geom->NumOfFluxCalculationGaussQuadPoints_West(iCell,jCell,faceWest);

      // Get total number of points which are used to assess the slope limiter
      NumGQP = GQP_North + GQP_South + GQP_East + GQP_West;

      // Allocate memory for the location of the points and the solution
      uQuad = new double [NumGQP];
      GQP = new Vector2D [NumGQP];
      
      // Get North face GQPs
      n = 0;
      if (faceNorth){
	// used GQPs from BndSplineInfo
	Geom->BndNorthSplineInfo[iCell].CopyGQPoints(&GQP[n]);
      } else {
	// used GQPs from straight edge
	Geom->getGaussQuadPointsFaceN(iCell,jCell,&GQP[n],GQP_North);
      }
      // Update n
      n += GQP_North;
      
      // Get West face GQPs
      if (faceWest){
	// used GQPs from BndSplineInfo
	Geom->BndWestSplineInfo[jCell].CopyGQPoints(&GQP[n]);
      } else {
	// used GQPs from straight edge
	Geom->getGaussQuadPointsFaceW(iCell,jCell,&GQP[n],GQP_West);
      }
      // Update n
      n += GQP_West;
      
      // Get South face GQPs
      if (faceSouth){
	// used GQPs from BndSplineInfo
	Geom->BndSouthSplineInfo[iCell].CopyGQPoints(&GQP[n]);
      } else {
	// used GQPs from straight edge
	Geom->getGaussQuadPointsFaceS(iCell,jCell,&GQP[n],GQP_South);
      }
      // Update n
      n += GQP_South;
      
      // Get East face GQPs
      if (faceEast){
	// used GQPs from BndSplineInfo
	Geom->BndEastSplineInfo[jCell].CopyGQPoints(&GQP[n]);
      } else {
	// used GQPs from straight edge
	Geom->getGaussQuadPointsFaceE(iCell,jCell,&GQP[n],GQP_East);
      }
    
      // Calculate the limiter for each solution variable (i.e. parameter)
      for (parameter = 1; parameter <= NumberOfVariables(); ++parameter) {
      
	// Drop the order only for the variables that are flagged as unfit
	if ( CellInadequateFitValue(iCell,jCell,parameter) ){
	
	  if (RecOrder() > 1){
	    // Zero all derivatives but D00 associated with this parameter.
	    for (n = 1; n < NumberOfTaylorDerivatives(); ++n){
	      CellTaylorDerivState(iCell,jCell,n)[parameter] = 0.0;
	    }
	  
	    // Set D00 and U_ave(see memory pool) of the current parameter to the correspondent average solution.
	    // U_ave is used in UnlimitedLinearSolutionAtLocation() (see below).
	    U_ave[parameter] = CellTaylorDerivState(iCell,jCell,0)[parameter] = (SolnBlk.*ReconstructedSoln)(iCell,jCell)[parameter];
	    
	    // Set D01 and D10 to the values of the first-order derivatives.
	    CellTaylorDerivValue(iCell,jCell,0,1,parameter) = dUdy[parameter];
	    CellTaylorDerivValue(iCell,jCell,1,0,parameter) = dUdx[parameter];
	    
	  } else {
	    // Set U_ave of the current parameter
	    U_ave[parameter] = CellTaylorDerivState(iCell,jCell,0)[parameter];
	    dUdy[parameter] = CellTaylorDerivState(iCell,jCell,1)[parameter];
	    dUdx[parameter] = CellTaylorDerivState(iCell,jCell,2)[parameter];
	  }

	  // Compute the minimum and maximum average solution in the stencil for the current parameter
	  u0Min = U_ave[parameter];
	  u0Max = u0Min;
	  for (n = 0; n < n_pts; ++n) {
	    u0Min = min(u0Min, (SolnBlk.*ReconstructedSoln)(I_Index[n],J_Index[n])[parameter]);
	    u0Max = max(u0Max, (SolnBlk.*ReconstructedSoln)(I_Index[n],J_Index[n])[parameter]);
	  }// endfor(n)
	  
	  // Evaluate the solution at all required points
	  for (n = 0; n < NumGQP; ++n){
	    uQuad[n] = UnlimitedLinearSolutionAtLocation(iCell,jCell,
							 GQP[n],
							 parameter);
	  }// endfor(n)

	  // Evaluate limiter for the current parameter
	  CellTaylorDeriv(iCell,jCell).Limiter(parameter) = CalculateLimiter(uQuad,NumGQP,U_ave[parameter],
									     u0Min,u0Max,Limiter);

	  // Save a copy of the limiter for the current parameter for later usage if limiter freezing is required
	  CellTaylorDeriv(iCell,jCell).Make_Limiter_Copy(parameter);
	
	}// endif

      }//endfor(parameter)

      // Deallocate memory
      delete [] uQuad; uQuad = NULL;
      delete [] GQP; GQP = NULL;
      NumGQP = 0;

    } else {
      // Set limiters to the frozen ones.

      for (parameter = 1; parameter <= NumberOfVariables(); ++parameter) {
      
	// Drop the order only for the variables that are flagged as unfit
	if ( CellInadequateFitValue(iCell,jCell,parameter) ){
	
	  if (RecOrder() > 1){
	    // Zero all derivatives but D00 associated with this parameter.
	    for (n = 1; n < NumberOfTaylorDerivatives(); ++n){
	      CellTaylorDerivState(iCell,jCell,n)[parameter] = 0.0;
	    }
	  
	    // Set D00 and U_ave(see memory pool) of the current parameter to the correspondent average solution.
	    // U_ave is used in UnlimitedLinearSolutionAtLocation() (see below).
	    U_ave[parameter] = CellTaylorDerivState(iCell,jCell,0)[parameter] = (SolnBlk.*ReconstructedSoln)(iCell,jCell)[parameter];
	    
	    // Set D01 and D10 to the values of the first-order derivatives.
	    CellTaylorDerivValue(iCell,jCell,0,1,parameter) = dUdy[parameter];
	    CellTaylorDerivValue(iCell,jCell,1,0,parameter) = dUdx[parameter];
	    
	  } else {
	    // Set U_ave of the current parameter
	    U_ave[parameter] = CellTaylorDerivState(iCell,jCell,0)[parameter];
	    dUdy[parameter] = CellTaylorDerivState(iCell,jCell,1)[parameter];
	    dUdx[parameter] = CellTaylorDerivState(iCell,jCell,2)[parameter];
	  }

	  // Set limiter for the current parameter to the frozen one
	  CellTaylorDeriv(iCell,jCell).Limiter(parameter) = CellTaylorDeriv(iCell,jCell).Frozen_Limiter(parameter);

	}// endif

      }//endfor(parameter)      

    } // endif (!SolnBlk.Freeze_Limiter)

  } else {
    // Use piecewise constant

    // Set D00 to average solution
    CellTaylorDerivState(iCell,jCell,0) = (SolnBlk.*ReconstructedSoln)(iCell,jCell);
    // Zero all derivatives
    for (n = 1; n < NumberOfTaylorDerivatives(); ++n){
      CellTaylorDerivState(iCell,jCell,n).Vacuum();
    }
    // Reset limiter
    CellTaylorDeriv(iCell,jCell).Limiter().Vacuum();
    
  } /* endif */
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
void HighOrder2D<SOLN_STATE>::
ComputeUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk, 
				       const int &iCell, const int &jCell,
				       const bool & UseSpecialStencil,
				       const Soln_State & 
				       (Soln_Block_Type::*ReconstructedSoln)(const int &,const int &) const){

  IndexType Special_I_Index(i_index.size()), Special_J_Index(j_index.size());

  if ( iCell >= StartI && iCell <= EndI && jCell >= StartJ && jCell <= EndJ ){

    /***************************************************************************
     *    Perform unconstrained unlimited high-order solution reconstruction   *
     **************************************************************************/

    if (UseSpecialStencil){

      // Set the stencil of points used for reconstruction with the special routine
      SetSpecialReconstructionStencil(iCell, jCell, Special_I_Index, Special_J_Index);

    } else {

      // Set the stencil of points used for reconstruction with the regular routine
      SetReconstructionStencil(iCell, jCell, Special_I_Index, Special_J_Index);
    }

    // Compute the reconstruction for the current cell
    ComputeUnconstrainedUnlimitedSolutionReconstruction(SolnBlk, ReconstructedSoln,
							iCell, jCell, Special_I_Index, Special_J_Index);
    
  } else {
    
    /***************************************************************************
     *    Perform constrained unlimited high-order solution reconstruction     *
     **************************************************************************/

    // Add constrained reconstruction here

    throw runtime_error("HighOrder2D<SOLN_STATE>::ComputeUnlimitedSolutionReconstruction() doesn't know how to handle constrained reconstruction!");

  }

}


/*! 
 * Compute the unlimited k-exact high-order reconstruction
 * proposed by Barth (1993) combined with equations satisfied 
 * exactly (i.e. constraints) which account for boundary conditions.
 * This reconstruction procedure should be used for computational cells
 * affected by the presence of a boundary condition enforced with constrained 
 * reconstruction.
 * The mean quantity conservation is also enforced as a constraint.
 *
 * \note This routine is customized for advection-diffusion state class!
 */
template<class SOLN_STATE>
template<class Soln_Block_Type>
void HighOrder2D<SOLN_STATE>::
ComputeConstrainedUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk,
						  const Soln_State & 
						  (Soln_Block_Type::*ReconstructedSoln)(const int &,const int &) const,
						  const int &iCell, const int &jCell,
						  const IndexType & i_index,
						  const IndexType & j_index) {

  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS
  int TotalNumberOfExactlySatisfiedConstraints(0),
    TotalNumberOfExactlySatisfiedEquations(1), // account for average conservation in cell (iCell,jCell)
    TotalNumberOfApproximatelySatisfiedConstraints(0),
    TotalNumberOfApproximatelySatisfiedEquations(i_index.size() - 1); // account for average conservation in neighbour cells

  int Temp, cell, n;
  int NumGQP(Geom->getNumGQP());

  IndexType ParameterIndex(1,1); //< for advection-diffusion!!!



  // == Check if the reconstruction polynomial is piecewise constant
  if (RecOrder() == 0){
    // *********  Assign the average solution to D00 ***********
    CellTaylorDerivState(iCell,jCell,0) = (SolnBlk.*ReconstructedSoln)(iCell,jCell);
    // There is no need to calculate the reconstruction
    return;
  }

  // Reset the memory pools
  Constraints_Loc.clear();
  Constraints_Normals.clear();
  Constraints_BCs.clear();
  Approx_Constraints_Loc.clear();
  Approx_Constraints_Normals.clear();
  Approx_Constraints_BCs.clear();
  

  // ========  Determine the number of exactly satisfied constraints and fetch the data  ==========
  // Check North cell face
  Temp = Geom->NumOfConstrainedGaussQuadPoints_North(iCell,jCell);
  if (Temp > 0){
    // Constraints detected on the North face
    TotalNumberOfExactlySatisfiedConstraints += Temp;

    // Fetch the data
    if (Geom->BndNorthSplineInfo != NULL){
      Geom->BndNorthSplineInfo[iCell].CopyGQPoints(Constraints_Loc);
      Geom->BndNorthSplineInfo[iCell].CopyNormalGQPoints(Constraints_Normals);      
    } else {
      Geom->addGaussQuadPointsFaceN(iCell,jCell,Constraints_Loc,NumGQP);
      for (n = 0; n < NumGQP; ++n){
	Constraints_Normals.push_back(Geom->nfaceN(iCell,jCell));
      }
    }
    
    Constraints_BCs.push_back(&SolnBlk.BC_NorthCell(iCell));
  }

  // Check South cell face
  Temp = Geom->NumOfConstrainedGaussQuadPoints_South(iCell,jCell);
  if (Temp > 0){
    // Constraints detected on the South face
    TotalNumberOfExactlySatisfiedConstraints += Temp;

    // Fetch the data
    if (Geom->BndSouthSplineInfo != NULL){
      Geom->BndSouthSplineInfo[iCell].CopyGQPoints(Constraints_Loc);
      Geom->BndSouthSplineInfo[iCell].CopyNormalGQPoints(Constraints_Normals);      
    } else {
      Geom->addGaussQuadPointsFaceS(iCell,jCell,Constraints_Loc,NumGQP);
      for (n = 0; n < NumGQP; ++n){
	Constraints_Normals.push_back(Geom->nfaceS(iCell,jCell));
      }
    }
    
    Constraints_BCs.push_back(&SolnBlk.BC_SouthCell(iCell));
  }
  
  // Check East cell face
  Temp = Geom->NumOfConstrainedGaussQuadPoints_East(iCell,jCell);
  if (Temp > 0){
    // Constraints detected on the East face
    TotalNumberOfExactlySatisfiedConstraints += Temp;

    // Fetch the data
    if (Geom->BndEastSplineInfo != NULL){
      Geom->BndEastSplineInfo[jCell].CopyGQPoints(Constraints_Loc);
      Geom->BndEastSplineInfo[jCell].CopyNormalGQPoints(Constraints_Normals);      
    } else {
      Geom->addGaussQuadPointsFaceE(iCell,jCell,Constraints_Loc,NumGQP);
      for (n = 0; n < NumGQP; ++n){
	Constraints_Normals.push_back(Geom->nfaceE(iCell,jCell));
      }
    }
    
    Constraints_BCs.push_back(&SolnBlk.BC_EastCell(jCell));
  }

  // Check West cell face
  Temp = Geom->NumOfConstrainedGaussQuadPoints_West(iCell,jCell);
  if (Temp > 0){
    // Constraints detected on the West face
    TotalNumberOfExactlySatisfiedConstraints += Temp;

    // Fetch the data
    if (Geom->BndWestSplineInfo != NULL){
      Geom->BndWestSplineInfo[jCell].CopyGQPoints(Constraints_Loc);
      Geom->BndWestSplineInfo[jCell].CopyNormalGQPoints(Constraints_Normals);      
    } else {
      Geom->addGaussQuadPointsFaceW(iCell,jCell,Constraints_Loc,NumGQP);
      for (n = 0; n < NumGQP; ++n){
	Constraints_Normals.push_back(Geom->nfaceW(iCell,jCell));
      }
    }
    
    Constraints_BCs.push_back(&SolnBlk.BC_WestCell(jCell));
  }


  if (CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_ADDITIONAL_APPROXIMATE_CONSTRAINTS == ON && 
      CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL == OFF) {

    // ======= Determine the number of approximately satisfied constraints and fetch the data ======
    for (cell = 1; cell < i_index.size(); ++cell){ // for each neighbour cell
    
      // Check North neighbour cell face
      Temp = Geom->NumOfConstrainedGaussQuadPoints_North(i_index[cell],j_index[cell]);
      if (Temp > 0){
	// Constraints detected on the North face
	TotalNumberOfApproximatelySatisfiedConstraints += Temp;

	// Fetch the data
	if (Geom->BndNorthSplineInfo != NULL){
	  Geom->BndNorthSplineInfo[i_index[cell]].CopyGQPoints(Approx_Constraints_Loc);
	  Geom->BndNorthSplineInfo[i_index[cell]].CopyNormalGQPoints(Approx_Constraints_Normals);      
	} else {
	  Geom->addGaussQuadPointsFaceN(i_index[cell],j_index[cell],Approx_Constraints_Loc,NumGQP);
	  for (n = 0; n < NumGQP; ++n){
	    Approx_Constraints_Normals.push_back(Geom->nfaceN(i_index[cell],j_index[cell]));
	  }
	}
    
	Approx_Constraints_BCs.push_back(&SolnBlk.BC_NorthCell(i_index[cell]));
      }

      // Check South neighbour cell face
      Temp = Geom->NumOfConstrainedGaussQuadPoints_South(i_index[cell],j_index[cell]);
      if (Temp > 0){
	// Constraints detected on the South face
	TotalNumberOfApproximatelySatisfiedConstraints += Temp;

	// Fetch the data
	if (Geom->BndSouthSplineInfo != NULL){
	  Geom->BndSouthSplineInfo[i_index[cell]].CopyGQPoints(Approx_Constraints_Loc);
	  Geom->BndSouthSplineInfo[i_index[cell]].CopyNormalGQPoints(Approx_Constraints_Normals);      
	} else {
	  Geom->addGaussQuadPointsFaceS(i_index[cell],j_index[cell],Approx_Constraints_Loc,NumGQP);
	  for (n = 0; n < NumGQP; ++n){
	    Approx_Constraints_Normals.push_back(Geom->nfaceS(i_index[cell],j_index[cell]));
	  }
	}
    
	Approx_Constraints_BCs.push_back(&SolnBlk.BC_SouthCell(i_index[cell]));
      }
  
      // Check East neighbour cell face
      Temp = Geom->NumOfConstrainedGaussQuadPoints_East(i_index[cell],j_index[cell]);
      if (Temp > 0){
	// Constraints detected on the East face
	TotalNumberOfApproximatelySatisfiedConstraints += Temp;

	// Fetch the data
	if (Geom->BndEastSplineInfo != NULL){
	  Geom->BndEastSplineInfo[j_index[cell]].CopyGQPoints(Approx_Constraints_Loc);
	  Geom->BndEastSplineInfo[j_index[cell]].CopyNormalGQPoints(Approx_Constraints_Normals);      
	} else {
	  Geom->addGaussQuadPointsFaceE(i_index[cell],j_index[cell],Approx_Constraints_Loc,NumGQP);
	  for (n = 0; n < NumGQP; ++n){
	    Approx_Constraints_Normals.push_back(Geom->nfaceE(i_index[cell],j_index[cell]));
	  }
	}
    
	Approx_Constraints_BCs.push_back(&SolnBlk.BC_EastCell(j_index[cell]));
      }

      // Check West neighbour cell face
      Temp = Geom->NumOfConstrainedGaussQuadPoints_West(i_index[cell],j_index[cell]);
      if (Temp > 0){
	// Constraints detected on the West face
	TotalNumberOfApproximatelySatisfiedConstraints += Temp;

	// Fetch the data
	if (Geom->BndWestSplineInfo != NULL){
	  Geom->BndWestSplineInfo[j_index[cell]].CopyGQPoints(Approx_Constraints_Loc);
	  Geom->BndWestSplineInfo[j_index[cell]].CopyNormalGQPoints(Approx_Constraints_Normals);      
	} else {
	  Geom->addGaussQuadPointsFaceW(i_index[cell],j_index[cell],Approx_Constraints_Loc,NumGQP);
	  for (n = 0; n < NumGQP; ++n){
	    Approx_Constraints_Normals.push_back(Geom->nfaceW(i_index[cell],j_index[cell]));
	  }
	}
    
	Approx_Constraints_BCs.push_back(&SolnBlk.BC_WestCell(j_index[cell]));
      }
    } // endfor (cell)

  } // endif

  /******** Determine dimensions of the least-squares problem and set matrices accordingly ************/
  /****************************************************************************************************/
  TotalNumberOfExactlySatisfiedEquations += TotalNumberOfExactlySatisfiedConstraints;
  TotalNumberOfApproximatelySatisfiedEquations += TotalNumberOfApproximatelySatisfiedConstraints;

  A_Assembled.newsize(TotalNumberOfExactlySatisfiedEquations + TotalNumberOfApproximatelySatisfiedEquations,
		      NumberOfTaylorDerivatives());
  All_U_Assembled.newsize(TotalNumberOfExactlySatisfiedEquations + TotalNumberOfApproximatelySatisfiedEquations,
			  NumberOfVariables());

  /************ Generate exactly satisfied individual constraints for UNCOUPLED variables *************/
  /****************************************************************************************************/
  Generalized_IndividualConstraints_Equations(SolnBlk,
					      iCell, jCell,
					      Constraints_Loc,
					      Constraints_Normals,
					      Constraints_BCs,
					      A_Assembled, All_U_Assembled,
					      ParameterIndex,
					      0, 0);


  /******** Generate approximately satisfied individual constraints for UNCOUPLED variables ***********/
  /****************************************************************************************************/
  Generalized_IndividualConstraints_Equations(SolnBlk,
					      iCell, jCell,
					      Approx_Constraints_Loc,
					      Approx_Constraints_Normals,
					      Approx_Constraints_BCs,
					      A_Assembled, All_U_Assembled,
					      ParameterIndex,
					      TotalNumberOfExactlySatisfiedEquations, 0);

  /******************** Generate the exact and approximate mean conservation equations ***************************/
  /***************************************************************************************************************/
  Set_MeanValueConservation_Equations(SolnBlk,
				      ReconstructedSoln,
				      iCell,jCell,
				      i_index, j_index,
				      A_Assembled, All_U_Assembled,
				      ParameterIndex,
				      TotalNumberOfExactlySatisfiedConstraints,
				      TotalNumberOfExactlySatisfiedEquations + TotalNumberOfApproximatelySatisfiedConstraints,
				      0);


  /************* Obtain solution to the constrained least-square problem *************/
  /***********************************************************************************/
  Solve_Constrained_LS_Householder(A_Assembled,
				   All_U_Assembled,
				   X_Assembled,
				   TotalNumberOfExactlySatisfiedEquations);

  // Update the coefficients D (derivatives)
  //**************************************************
  for (n=0; n<=CellTaylorDeriv(iCell,jCell).LastElem(); ++n){
    CellTaylorDerivState(iCell,jCell,n)[1] = X_Assembled(n,0);
  }//endfor

}

/*!
 * Generate the LHS and RHS of the least-squares problem associated 
 * with the reconstruction procedure and write the values at the specified
 * locations.
 * This matrix reflects the conservation of mean value in the control       
 * volumes of cells specified by (i_index,j_index).                         
 * The row associated with cell (iCell,jCell) represents a constraint and    
 * therefore it must be satisfied exactly.
 * The rest of the system is approximately solved in the least squares sense.
 *                                                                          
 * \note It is assumed that the (i_index[0],j_index[0]) is equal to          
 *       (iCell,jCell). That is, the first line of A is a constrain!!!            
 * The constraint is filled in the matrix A at the position (RowConstraint, StartCol).                                     
 * The approximate equations are filled in the matrix starting from (StartRow, StartCol) position. 
 *
 * \param SolnBlk the quad block for which the solution reconstruction is done.
 * \param ReconstructedSoln member function of Soln_Block_Type which returns the solution.
 * \param iCell i-index of the reconstructed cell
 * \param jCell j-index of the reconstructed cell
 * \param ParameterIndex related to the indexes of the solution
 * \param A the LHS assemble matrix 
 * \param B the RHS assemble matrix
 *
 * \note This routine doesn't normalize the geometric weights.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type> inline
void HighOrder2D<SOLN_STATE>::
Set_MeanValueConservation_Equations(Soln_Block_Type & SolnBlk,
				    const Soln_State & 
				    (Soln_Block_Type::*ReconstructedSoln)(const int &,const int &) const,
				    const int &iCell, const int &jCell,
				    const IndexType & i_index, const IndexType & j_index,
				    DenseMatrix & A, DenseMatrix & All_U,
				    const IndexType & ParameterIndex,
				    const int &RowConstraint,
				    const int &StartRow, const int &StartCol){


  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS

  int StencilSize(i_index.size());
  int ParameterIndexSize(ParameterIndex.size());
  ColumnVector GeomWeights(StencilSize);   // The column vector of the geometric weights
  Vector2D *DeltaCellCenters;              /* stores the difference between the cell center of
					      neighbour cells and the one of (i,j) cell.
					      The first value is 0!!! */
  int IndexSumY, IndexSumX, P1, P2;
  double CombP1X, CombP2Y;
  double PowDistanceYC, PowDistanceXC;
  int cell, i, parameter;
  double IntSum(0.0);

  // Allocate memory
  DeltaCellCenters = new Vector2D [StencilSize];

  // START:   Set the bottom part of the LHS and RHS of the linear system 
  // *********************************************************************

  // Step1. Set the constraint equation
  for (i=0; i <= CellTaylorDeriv(iCell,jCell).LastElem(); ++i){
    A(RowConstraint,i+StartCol) = Geom->CellGeomCoeffValue(iCell,jCell,i);
  }

  for (parameter=0; parameter<ParameterIndexSize; ++parameter){
    All_U(RowConstraint,parameter) = (SolnBlk.*ReconstructedSoln)(iCell,jCell)[ParameterIndex[parameter]];
  }

  // Step3. Set the approximate equations
  for (cell=1; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil

    /* Compute the X and Y component of the distance between
       the cell center of the neighbours and the reconstructed cell */
    DeltaCellCenters[cell] = CellCenter(i_index[cell],j_index[cell]) - CellCenter(iCell,jCell);

    /* Compute the geometric weight based on the centroid distance */
    CENO_Geometric_Weighting(GeomWeights[cell], DeltaCellCenters[cell].abs());

    // *** SET the matrix A of the linear system (LHS) ***
    /* compute for each derivative the corresponding entry in the matrix of the linear system */
    for (i=0; i<=CellTaylorDeriv(iCell,jCell).LastElem(); ++i){
      // build the row of the matrix
      P1 = CellTaylorDeriv(iCell,jCell,i).P1();  // identify P1
      P2 = CellTaylorDeriv(iCell,jCell,i).P2();  // identify P2
      A(cell+StartRow-1,i+StartCol) = 0.0;  // set sumation variable to zero
      CombP2Y = 1.0;                        // the binomial coefficient "nC k" for k=0 is 1
      PowDistanceYC = 1.0; 	            // initialize PowDistanceYC

      // Compute geometric integral over the neighbour's domain
      for (IndexSumY = 0; IndexSumY<=P2; ++IndexSumY){
	CombP1X = 1.0;       // the binomial coefficient "nC k" for k=0 is 1
	PowDistanceXC = 1.0; // initialize PowDistanceXC
	IntSum = 0.0;	     // reset internal sumation variable

	for (IndexSumX = 0; IndexSumX<=P1; ++IndexSumX){
	  IntSum += ( CombP1X*PowDistanceXC*
		      Geom->CellGeomCoeffValue(i_index[cell],j_index[cell],P1-IndexSumX,P2-IndexSumY) );
	    
	  // update the binomial coefficients
	  CombP1X = (P1-IndexSumX)*CombP1X/(IndexSumX+1);  // The index is still the old one => expression for "nC k+1"
	  PowDistanceXC *= DeltaCellCenters[cell].x;       // Update PowDistanceXC
	}//endfor (IndexSumX)

	A(cell+StartRow-1,i+StartCol) += CombP2Y*PowDistanceYC*IntSum; // update the external sum

	CombP2Y = (P2-IndexSumY)*CombP2Y/(IndexSumY+1); // the index is still the old one => expression for "nC k+1"
	PowDistanceYC *= DeltaCellCenters[cell].y;      // Update PowDistanceYC
      }//endfor (IndexSumY)

      // apply geometric weighting
      A(cell+StartRow-1,i+StartCol) *= GeomWeights(cell);

    }//endfor (i)
      
    // *** SET the matrix All_U of the linear system (RHS) ***
    for (parameter=0; parameter < ParameterIndexSize; ++parameter){
      All_U(cell+StartRow-1,parameter) = ( GeomWeights(cell)*
					   (SolnBlk.*ReconstructedSoln)(i_index[cell],j_index[cell])[ParameterIndex[parameter]] );
    }
  } //endfor (cell)


  delete [] DeltaCellCenters;

}

/*!
 * Generate the LHS of the least-squares problem associated with the 
 * reconstruction procedure and write the values at the specified locations.
 * This matrix reflects the conservation of mean value in the control
 * volumes of cells specified by (i_index,j_index).
 * The row associated with cell (iCell,jCell) represents a constraint and    
 * therefore it must be satisfied exactly.
 * The rest of the system is approximately solved in the least squares sense.
 *                                                                          
 * \note It is assumed that the (i_index[0],j_index[0]) is equal to          
 *       (iCell,jCell). That is, the first line of A is a constraint!!!
 * The approximate equations are filled in the rest of the matrix.
 *
 * \param iCell i-index of the reconstructed cell
 * \param jCell j-index of the reconstructed cell
 * \param A the LHS matrix
 * \param CellGeometricWeights the array of geometric weights
 *
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::
Set_LHS_MeanValueConservation_Equations(const int &iCell, const int &jCell,
					const IndexType & i_index, const IndexType & j_index,
					DenseMatrix & A,
					DoubleArrayType & GeometricWeights){


  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS

  int StencilSize(i_index.size());
  Vector2D *DeltaCellCenters;              /* stores the difference between the cell center of
					      neighbour cells and the one of (i,j) cell.
					      The first value is 0!!! */
  int IndexSumY, IndexSumX, P1, P2;
  double CombP1X, CombP2Y;
  double PowDistanceYC, PowDistanceXC;
  int cell, i;
  double IntSum(0.0);
  double MaxWeight(0.0);

  // Allocate memory
  DeltaCellCenters = new Vector2D [StencilSize];

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
    
    /* Compute the X and Y component of the distance between
       the cell centers of the neighbour and the reconstructed cell */
    DeltaCellCenters[cell] = CellCenter(i_index[cell],j_index[cell]) - CellCenter(iCell,jCell);
    
    /* Compute the geometric weight based on the centroid distance */
    CENO_Geometric_Weighting(GeometricWeights[cell], DeltaCellCenters[cell].abs());
    
    /* Compute the maximum geometric weight (this is used for normalization) */
    MaxWeight = max(MaxWeight, GeometricWeights[cell]);
  }

  // Step 2. Set the constraint equation
  for (i=0; i <= CellTaylorDeriv(iCell,jCell).LastElem(); ++i){
    A(0,i) = Geom->CellGeomCoeffValue(iCell,jCell,i);
  }

  // Step 3. Set the approximate equations
  for (cell=1; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil

    // compute the normalized geometric weight
    GeometricWeights[cell] /= MaxWeight;

    // *** SET the matrix A of the linear system (LHS) ***
    /* compute for each derivative the corresponding entry in the matrix of the linear system */
    for (i=0; i<=CellTaylorDeriv(iCell,jCell).LastElem(); ++i){
      // build the row of the matrix
      P1 = CellTaylorDeriv(iCell,jCell,i).P1();  // identify P1
      P2 = CellTaylorDeriv(iCell,jCell,i).P2();  // identify P2
      A(cell,i) = 0.0;                      // set sumation variable to zero
      CombP2Y = 1.0;                        // the binomial coefficient "nC k" for k=0 is 1
      PowDistanceYC = 1.0; 	            // initialize PowDistanceYC

      // Compute geometric integral over the neighbour's domain
      for (IndexSumY = 0; IndexSumY<=P2; ++IndexSumY){
	CombP1X = 1.0;       // the binomial coefficient "nC k" for k=0 is 1
	PowDistanceXC = 1.0; // initialize PowDistanceXC
	IntSum = 0.0;	     // reset internal sumation variable

	for (IndexSumX = 0; IndexSumX<=P1; ++IndexSumX){
	  IntSum += ( CombP1X*PowDistanceXC*
		      Geom->CellGeomCoeffValue(i_index[cell],j_index[cell],P1-IndexSumX,P2-IndexSumY) );
	    
	  // update the binomial coefficients
	  CombP1X = (P1-IndexSumX)*CombP1X/(IndexSumX+1);  // The index is still the old one => expression for "nC k+1"
	  PowDistanceXC *= DeltaCellCenters[cell].x;       // Update PowDistanceXC
	}//endfor (IndexSumX)

	A(cell,i) += CombP2Y*PowDistanceYC*IntSum; // update the external sum

	CombP2Y = (P2-IndexSumY)*CombP2Y/(IndexSumY+1); // the index is still the old one => expression for "nC k+1"
	PowDistanceYC *= DeltaCellCenters[cell].y;      // Update PowDistanceYC
      }//endfor (IndexSumY)

      // apply geometric weighting
      A(cell,i) *= GeometricWeights[cell];

    }//endfor (i)
      
  } //endfor (cell)


  delete [] DeltaCellCenters;

}

/*!
 * Set the constraints equations in the constrained least-square reconstruction.
 * This constraints are called individual because they are independent of
 * other solution parameters
 * (i.e they don't express any relationship between solution parameters).
 * The starting position for the entries in the LHS and RHS of the linear
 * system are specified by the (StartRow, StartCol) input parameters.   \n                                
 *                                                                          
 * The BCs handled by this subroutine have the mixed form:                  
 *   F(GQP, A_GQP, B_GQP) = A_GQP*F1(GQP) + B_GQP*F2(GQP)                   
 *                                                                          
 * where: GQP      -- Gauss Quadrature Point Locations (i.e. flux calculation points)
 *        A_GQP    -- The coefficient for the Dirichlet boundary condition at GQP                                                
 *        B_GQP    -- The coefficient for the Neumann boundary condition at GQP                                                
 *        F1(GQP)  -- The value of the Dirichlet boundary condition at GQP  
 *        F2(GQP)  -- The value of the Neumann boundary condition at GQP    
 *
 * \param SolnBlk the quad block for which the solution reconstruction is done.
 * \param iCell i-index of the reconstructed cell
 * \param jCell j-index of the reconstructed cell
 * \param Constraints_Loc GQP array
 * \param Constraints_Normals normal vectors at GQP locations
 * \param Constraints_BCs provide the boundary condition coefficients (i.e. A_GQP, B_GQP, F1, F2)
 * \param ParameterIndex related to the indexes of the solution
 * \param A the LHS assemble matrix 
 * \param B the RHS assemble matrix
 *
 * \note This routine is customized for advection-diffusion state class!
 */
template<class SOLN_STATE>
template<class Soln_Block_Type> inline
void HighOrder2D<SOLN_STATE>::
Generalized_IndividualConstraints_Equations(Soln_Block_Type & SolnBlk,
					    const int &iCell, const int &jCell,
					    Vector2DArray & Constraints_Loc,
					    Vector2DArray & Constraints_Normals,
					    BC_Type_Array & Constraints_BCs,
					    DenseMatrix & A, DenseMatrix & All_U,
					    const IndexType & ParameterIndex,
					    const int &StartRow, const int &StartCol) {
  
  int P1, P2, i;
  int BCs, BCs_Entry, Eq;			// equation
  double PowXC, PowYC;		/* PowXC = DistXi^(P1-1); PowYC = DistYi^(P2-1) */
  double DistXi, DistYi;
  int IndexP1, IndexP2;
  double GeometricWeight;

  
  for (BCs_Entry = 0, Eq = 0; BCs_Entry < Constraints_BCs.size(); ++BCs_Entry){ // for each boundary condition entry

    for (BCs = 1; BCs <= Constraints_BCs[BCs_Entry]->NumOfPoints(); ++BCs, ++Eq){ // for each flux calculation point

      // Compute entrie in the LHS and RHS for the current equation
      
      // Determine distance between the current GQP and the centroid of cell (iCell,jCell)
      DistXi = Constraints_Loc[Eq].x - XCellCenter(iCell,jCell);
      DistYi = Constraints_Loc[Eq].y - YCellCenter(iCell,jCell);

      /* Compute the geometric weight based on the centroid distance */
      CENO_Geometric_Weighting(GeometricWeight, Vector2D(DistXi,DistYi).abs());

      
      // Step 1. Form the LHS  -- build the row of the matrix A associated with the current GQP
      for (i=0; i<=CellTaylorDeriv(iCell,jCell).LastElem(); ++i){
	// build the row of the matrix
	P1 = CellTaylorDeriv(iCell,jCell,i).P1();  // identify P1
	P2 = CellTaylorDeriv(iCell,jCell,i).P2();  // identify P2

	/* Initialize PowXC & PowYC */
	PowXC = 1.0/DistXi;
	PowYC = 1.0/DistYi;

	/* Update PowXC & PowYC */
	for (IndexP1 = 1; IndexP1 <= P1; ++IndexP1){ PowXC *= DistXi; }
	for (IndexP2 = 1; IndexP2 <= P2; ++IndexP2){ PowYC *= DistYi; }

	A(Eq+StartRow,i+StartCol) = ( PowXC * PowYC * 
				      ( Constraints_BCs[BCs_Entry]->a(BCs)[ParameterIndex[0]] * DistXi * DistYi + 
					Constraints_BCs[BCs_Entry]->b(BCs)[ParameterIndex[0]] * (P1 * DistYi * 
												 Constraints_Normals[Eq].x  + 
												 P2 * DistXi * 
												 Constraints_Normals[Eq].y )) );

	// Apply geometric weighting
	A(Eq+StartRow,i+StartCol) *= GeometricWeight;
      } //endfor (i)
       

      // Step 2. Form the RHS  -- build the row of the matrix All_U associated with the current GQP
      for (i=0; i<ParameterIndex.size(); ++i){
	All_U(Eq+StartRow, i) = ( ( Constraints_BCs[BCs_Entry]->a(BCs)[ParameterIndex[0]] * 
				    Constraints_BCs[BCs_Entry]->DirichletBC(BCs)[ParameterIndex[0]] ) + 
				  ( Constraints_BCs[BCs_Entry]->b(BCs)[ParameterIndex[0]] * 
				    Constraints_BCs[BCs_Entry]->NeumannBC(BCs)[ParameterIndex[0]] ) );

	// Apply geometric weighting
	All_U(Eq+StartRow, i) *= GeometricWeight;

      }//endfor (i)

    } // endfor (BCs)

  }// endfor (BCs_Entry) 

}


/*! 
 * Compute the unlimited k-exact high-order reconstruction
 * proposed by Barth (1993) combined with equations satisfied 
 * exactly (i.e. constraints) which account for boundary conditions.
 * This reconstruction procedure should be used for computational cells
 * affected by the presence of a boundary condition enforced with constrained 
 * reconstruction.
 * The mean quantity conservation is also enforced as a constraint.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type>
void HighOrder2D<SOLN_STATE>::
HighLevel_ConstrainedUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk,
						     const Soln_State & 
						     (Soln_Block_Type::*ReconstructedSoln)(const int &,
											   const int &) const,
						     const int &iCell, const int &jCell,
						     const IndexType & i_index,
						     const IndexType & j_index) {

  
  int constrGQP_W, constrGQP_E, constrGQP_N, constrGQP_S; // number of constrained points on each edge
  bool IC_Flag(false), RC_Flag(false);			  // flags to indicate the type of constraints encountered

  /**********************************************************************************
   *  STEP 1. DETERMINE THE NUMBER AND TYPE OF CONSTRAINTS THAT NEED TO BE IMPOSED  *
   *********************************************************************************/
  
  // Determine the number of constraints for the West boundary.
  constrGQP_W = SolnBlk.Grid.NumOfConstrainedGaussQuadPoints_West(iCell,jCell);
  if (constrGQP_W > 0){
    // Identify what type of constraints are present
    if ( SolnBlk.BC_WestCell(jCell).IsThereAnyRelationalConstraintRequired() ){
      RC_Flag = true;
    }
    if ( SolnBlk.BC_WestCell(jCell).IsThereAnyIndividualConstraintRequired() ){
      IC_Flag = true;
    }
  }

  // Determine the number of constraints and for the South boundary
  constrGQP_S = SolnBlk.Grid.NumOfConstrainedGaussQuadPoints_South(iCell,jCell);
  if (constrGQP_S > 0){
    // Identify what type of constraints are present
    if ( SolnBlk.BC_SouthCell(iCell).IsThereAnyRelationalConstraintRequired() ){
      RC_Flag = true;
    }
    if ( SolnBlk.BC_SouthCell(iCell).IsThereAnyIndividualConstraintRequired() ){
      IC_Flag = true;
    }
  }    

  // Determine the number of constraints and for the East boundary
  constrGQP_E = SolnBlk.Grid.NumOfConstrainedGaussQuadPoints_East(iCell,jCell);
  if (constrGQP_E > 0){
    // Identify what type of constraints are present
    if ( SolnBlk.BC_EastCell(jCell).IsThereAnyRelationalConstraintRequired() ){
      RC_Flag = true;
    }
    if ( SolnBlk.BC_EastCell(jCell).IsThereAnyIndividualConstraintRequired() ){
      IC_Flag = true;
    }
  }
    
  // Determine the number of constraints and for the North boundary
  constrGQP_N = SolnBlk.Grid.NumOfConstrainedGaussQuadPoints_North(iCell,jCell);
  if (constrGQP_N > 0){
    // Identify what type of constraints are present
    if ( SolnBlk.BC_NorthCell(iCell).IsThereAnyRelationalConstraintRequired() ){
      RC_Flag = true;
    }
    if ( SolnBlk.BC_NorthCell(iCell).IsThereAnyIndividualConstraintRequired() ){
      IC_Flag = true;
    }
  }

  
  /***********************************************************************************************
   * STEP 2. PROCEED WITH A DIFFERENT ALGORITHM DEPENDING ON THE TYPE OF ENCOUNTERED CONSTRAINTS *
   **********************************************************************************************/
  
  if (RC_Flag && IC_Flag){
    // There are both relational and individual constraints
    throw runtime_error("HighOrder2D<SOLN_STATE>::HighLevel_ConstrainedUnlimitedSolutionReconstruction() ERROR! Both relational and individual constraints case cannot be handled yet!");
    
  } else if (RC_Flag){
    // There are only relational constraints
    throw runtime_error("HighOrder2D<SOLN_STATE>::HighLevel_ConstrainedUnlimitedSolutionReconstruction() ERROR! Only relational constraints case cannot be handled yet!");

  } else if (IC_Flag){
    // There are only individual constraints
    throw runtime_error("HighOrder2D<SOLN_STATE>::HighLevel_ConstrainedUnlimitedSolutionReconstruction() ERROR! Only individual constraints case cannot be handled yet!");

  } else {
    // There are no constraints
    ComputeUnconstrainedUnlimitedSolutionReconstruction(SolnBlk,
							ReconstructedSoln,
							iCell, jCell,
							i_index, j_index);
  } // endif

}


#endif
