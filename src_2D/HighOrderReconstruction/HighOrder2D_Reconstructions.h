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
  } 

  // Check EAST boundary
  if (_constrained_EAST_reconstruction){
    // Add constrained reconstruction here
  } 

  // Check NORTH boundary
  if (_constrained_NORTH_reconstruction){
    // Add constrained reconstruction here
  } 

  // Check SOUTH boundary
  if (_constrained_SOUTH_reconstruction){
    // Add constrained reconstruction here
  }

}

/*! 
 * Compute the pseudo-inverse of the left-hand-side term in the 
 * unlimited k-exact high-order reconstruction for all computational
 * cells based on the information provided by the associated grid.
 *
 * \note The pseudo-inverse in not determined for cells which use 
 * constrained reconstruction.
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

    // Calculate the pseudo-inverse for cells in the specified range
    for ( j  = StartJ ; j <= EndJ ; ++j ) {
      for ( i = StartI ; i <= EndI ; ++i ) {
	
	// Set the stencil of points used for reconstruction
	SetReconstructionStencil(i, j, i_index, j_index);

	// Compute the pseudo-inverse for the current cell
	ComputeCellReconstructionPseudoInverse(i, j, i_index, j_index);
      }/* endfor */
    }/* endfor */

    
    // Confirm the pseudo-inverse calculation
    _calculated_psinv = true;
  }

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


    // Calculate slope limiters.
    if (!SolnBlk.Freeze_Limiter) {

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
	
	}// endif

      }//endfor(parameter)

      // Deallocate memory
      delete [] uQuad;
      delete [] GQP;
      NumGQP = 0;
    }

  } else {
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


#endif
