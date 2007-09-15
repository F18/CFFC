#ifndef _RECONSTRUCTION_FUNCTION_1D_INCLUDED
#define _RECONSTRUCTION_FUNCTION_1D_INCLUDED

/* Include defined header file. */
#include <vector>
#include "../../../src_2D/Math/LinearSystems.h"
#include "../../../src_2D/Utilities/Utilities.h"
#include "include/TypeDefinition.h"
#include "Reconstruction/ReconstructionHelpers.h"


/**********************************************************************************************************************
 ----------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------------
 *                                                                                                                    *
 *                                         1D Reconstructions                                                         *
 *                                                                                                                    *
 ----------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------------
 **********************************************************************************************************************/

 /***************************************************************************
 * TEMPLATIZED Function: kExact_Reconstruction for 1D                       *
 *                                                                          *
 * This subroutine determines the coefficients of a Taylor series expansion *
 * which approximates the solution over the domain of the cell specified    *
 * by "iCell" index.                                                        *
 *                                                                          *
 ***************************************************************************/
template< class SolutionContainer> inline
void kExact_Reconstruction (SolutionContainer & SolnBlk, IndexType & i_index,
			    const int iCell, int ND, int DroppedOrderToPWL)
{

  /* Obs:
     ND --> number of Taylor expansion coefficients
     DroppedOrderToPWL == OFF: indicates that the kExact_Reconstruction is used for determining the high-order reconstruction
     DroppedOrderToPWL == ON: indicates that the kExact_Reconstruction is used for determining the PWL reconstruction
                              because the high-order reconstruction is not smooth.
                              The variables used for storing this new reconstruction are different than the ones used
                              for storing the high-order reconstruction.
  */

  typedef typename SolutionContainer::CompCellType  ComputationalCellType;
  static const int NumberOfParameters = ComputationalCellType::NumberOfVariables;

  // SET THE VARIABLES USED IN THE RECONSTRUCTION PROCESS
  int StencilSize(i_index.size());                            /* number of cells used in the reconstruction stencil */
  DenseMatrix  A(StencilSize-1,ND-1);                         /* the matrix which the linear system is solved for */
  DenseMatrix  All_Delta_U(StencilSize-1,NumberOfParameters); /* matrix for storing U[neighbour]-U[cell] */
  ColumnVector DistanceCellCenters(StencilSize);              /* stores the distance between the cell center of neighbour
  						                 cells and the I cell */
  ColumnVector GeomWeights(StencilSize);                      /* The column vector of the geometric weights */
  int krank;                                                  /* the final rank of A is returned here */

  int cell, i, IndexSumX, parameter, P1, IndexLastDerivative;
  double CombP1X, WeightsSum, PowDistanceXC;

  // *********  Assign the average solution to D0  ***********
  if (DroppedOrderToPWL == OFF){
    /* Use the variables for high-order */
    SolnBlk(iCell).CellDeriv(0,true,true,true).D( ) = SolnBlk(iCell).CellSolution( );
    IndexLastDerivative = SolnBlk(iCell).CellDeriv().LastElem();
  } else {
    /* Use the variables for PWL */
    SolnBlk(iCell).CellDerivFirstOrder(0,true,true,true).D( ) = SolnBlk(iCell).CellSolution( );
    IndexLastDerivative = SolnBlk(iCell).CellDerivFirstOrder().LastElem();
  }

  if (SolnBlk(iCell).NumberOfTaylorDerivatives() == 1){ // piecewise constant
    return;
  }

  // START:   Set the LHS and RHS of the linear system 
  // ***************************************************
  // Step1. Compute the normalized geometric weights
  for (cell=1; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil

    /* Compute the distance between the cell centers of the neighbour and the reconstructed cell */
    DistanceCellCenters(cell) = SolnBlk(i_index[cell]).CellCenter() - SolnBlk(iCell).CellCenter();

    /* Compute the geometric weights and their sum (this is used for normalization)
       based on the distance to each control volume */
    GeomWeights(cell) = DistanceCellCenters(cell)*DistanceCellCenters(cell);
    GeomWeights(cell) = 1.0/(1.0E-15 + GeomWeights(cell));
    WeightsSum += GeomWeights(cell);
  }

  // Step2. Set the approximate equations
  for (cell=1 ; cell<StencilSize; ++cell){ //for each cell in the stencil
    
    // compute the normalized geometric weight
    GeomWeights(cell) /= WeightsSum;

    // *** SET the matrix A of the linear system (LHS) ***
    /* compute for each derivative the corresponding entry in the matrix of the linear system */
    for (i=1; i<=IndexLastDerivative; ++i){
      // build the row of the matrix
      P1 = SolnBlk(iCell).CellDeriv(i,true,true,true).P1(); // identify P1
      A(cell-1,i-1) = 0.0;  // set sumation variable to zero
      CombP1X = 1.0;        // the binomial coefficient "nC k" for k=0 is 1
      PowDistanceXC = 1.0;  // initialize PowDistanceXC

      // ********* Compute geometric integral over the neighbour's domain **********
      for (IndexSumX = 0; IndexSumX<=P1; ++IndexSumX){
	A(cell-1,i-1) += ( CombP1X * PowDistanceXC *
			   SolnBlk(i_index[cell]).CellGeomCoeff()(P1-IndexSumX) );
	
	// update the binomial coefficients
	CombP1X = (P1-IndexSumX)*CombP1X/(IndexSumX+1);  // the index is still the old one => expression for "nC k+1"
	PowDistanceXC *= DistanceCellCenters(cell);      // Update PowDistanceXC
      } //endfor

      // subtract the corresponding geometric moment of (iCell) cell
      A(cell-1,i-1) -= SolnBlk(iCell).CellGeomCoeff(i,true,true,true);

      // apply geometric weighting
      A(cell-1,i-1) *= GeomWeights(cell);
    }
    
    // *** SET the matrix All_Delta_U of the linear system (RHS) ***
    for (parameter = 1; parameter <= NumberOfParameters; ++parameter){
      All_Delta_U(cell-1,parameter-1) = ( SolnBlk(i_index[cell]).CellSolution(parameter) -
					  SolnBlk(iCell).CellSolution(parameter) );
      All_Delta_U(cell-1,parameter-1) *= GeomWeights(cell);
    }
  }
  // Matrix A of the linear system (LHS) built. The same matrix is used for all the variables (same geometry).
  // Matrix All_Delta_U of the linear system (RHS) built.
  // **********************************************************************

  /* Solve the overdetermined linear system of equations using a least-squares procedure*/
  /**************************************************************************************/
   Solve_LS_Householder_F77(A, All_Delta_U, krank, NumberOfParameters, StencilSize-1, ND-1);

  // Update the coefficients D (derivatives)
  //**************************************************
  if (DroppedOrderToPWL == OFF){
    for (i=1; i<=IndexLastDerivative; ++i){
      for (parameter = 1; parameter <= NumberOfParameters; ++parameter){
	SolnBlk(iCell).CellDeriv(i,true,true,true).D(parameter) = All_Delta_U(i-1,parameter-1);
      
	/* this equation makes sure that the mean conservation of each parameter is satisfied inside the reconstructed cell */
	SolnBlk(iCell).CellDeriv(0,true,true,true).D(parameter) -= ( SolnBlk(iCell).CellGeomCoeff(i,true,true,true) * 
								     All_Delta_U(i-1,parameter-1) );
      }
    }
  } else {
    for (i=1; i<=IndexLastDerivative; ++i){
      for (parameter = 1; parameter <= NumberOfParameters; ++parameter){
	SolnBlk(iCell).CellDerivFirstOrder(i,true,true,true).D(parameter) = All_Delta_U(i-1,parameter-1);
      }
    }
  }

}

 /***************************************************************************
 * TEMPLATIZED Function: Essentially NonOscillatory Reconstruction for 1D   *
 *                                                                          *
 * This function determines the coefficients of a Taylor series expansion   *
 * which approximates the solution over the domain of the cell specified    *
 * by "i" index.                                                            *
 * The method used is presented in:                                         *
 * "Uniformly High Order Accurate Essentially Non-oscillatory Schemes, III" *
 * Ami Harten, Journal of Computational Physics, 71, (1987)                 *
 *                                                                          *
 ***************************************************************************/

template< class SolutionContainer> inline
void ENO_Reconstruction (SolutionContainer & SolnBlk, const int iCell)
{

  typedef typename SolutionContainer::CompCellType  ComputationalCellType;
  static const int NumberOfParameters = ComputationalCellType::NumberOfVariables;

  int StencilStartIndex = iCell;
  int StencilPoints = 1;
  int DivDiffOrder=1;
  int NumOfUnknowns = SolnBlk(iCell).NumberOfTaylorDerivatives();

  DenseMatrix A(NumOfUnknowns,NumOfUnknowns);
  ColumnVector U(NumOfUnknowns);
  ColumnVector X(NumOfUnknowns);
  double DeltaCellCenter;
  int DerivOrder, P1, IndexSumX;
  double CombP1X;
  double GeomIntegralOverCell;

  for(int parameter=1; parameter <= NumberOfParameters; ++parameter){

    /* Determine the "smoothest" stencil */
    while(StencilPoints < NumOfUnknowns){

      // determine the starting position of the stencil    
      if (fabs(DivDifference(SolnBlk,parameter,DivDiffOrder,StencilStartIndex-1)) <
	  fabs(DivDifference(SolnBlk,parameter,DivDiffOrder,StencilStartIndex))){

	// change starting position of the stencil
	--StencilStartIndex;
      }

      /* advance to the next level of finite differences */
      ++DivDiffOrder;

      /* add one more point to the stencil */
      ++StencilPoints;
    }

    /* Solve the reconstruction using the obtained stencil */

    /* Set the LHS and RHS */
    for (int cell=StencilStartIndex, Counter=0; Counter<NumOfUnknowns; ++cell, ++Counter){
      //for each cell in the stencil

      /*** SET the matrix A of the linear system (LHS) ***/

      /* Compute the distance between the cell center of the current cell and the reconstructed cell */
      DeltaCellCenter = SolnBlk(cell).CellCenter() - SolnBlk(iCell).CellCenter();
   
      /* compute for each derivative the corresponding entry in the matrix of the linear system */
      for(int i = SolnBlk(iCell).CellDeriv().FirstElem(); i<=SolnBlk(iCell).CellDeriv().LastElem(); ++i){

	// build the row of the matrix
	GeomIntegralOverCell = 0.0;  // set sumation variable to zero
	// Compute geometric integral over the neighbour's domain
	CombP1X = 1.0; // the binomial coefficient "nC k" for k=0 is 1
	for (IndexSumX = 0; IndexSumX<=SolnBlk(iCell).CellDeriv(i,true,true,true).P1(); ++IndexSumX){
	  GeomIntegralOverCell += (CombP1X * pow(DeltaCellCenter,IndexSumX)*
				       SolnBlk(cell).CellGeomCoeff()(SolnBlk(iCell).CellDeriv(i,true,true,true).P1()-IndexSumX));
	  // update the binomial coefficients
	  CombP1X = (SolnBlk(iCell).CellDeriv(i,true,true,true).P1()-IndexSumX)*CombP1X/(IndexSumX+1);
	  // the index is still the old one => expression for "nC k+1"
	}

	A(Counter,i) = GeomIntegralOverCell; 
      }
      
      // *** SET U of the linear system (RHS) ***
      U(Counter) = SolnBlk(cell).CellSolution(parameter);

    }

    /* Solve the linear system */
    Solve_LU_Decomposition(A,U,X);

    /* Update the coefficients D (derivatives) */
    for(int i = SolnBlk(iCell).CellDeriv().FirstElem(); i<=SolnBlk(iCell).CellDeriv().LastElem(); ++i){
      SolnBlk(iCell).CellDeriv(i,true,true,true) = X(i);
    }

    /* Reset stencil variables for the reconstruction of the next parameter */
    StencilStartIndex = iCell;
    StencilPoints = 1;
    DivDiffOrder = 1;
  }

}

template< class SolutionContainer> inline
double DivDifference (SolutionContainer & SolnBlk, int Var, int Order, int SI)
{

  // SI - StartIndex

  switch(Order){

  case 1:
    return (SolnBlk(SI).CellSolution(Var) - SolnBlk(SI+1).CellSolution(Var))/
              (SolnBlk(SI).CellCenter() - SolnBlk(SI+1).CellCenter());
    break;

  case 2: 
    return ((SolnBlk(SI).CellSolution(Var) - SolnBlk(SI+1).CellSolution(Var))*
            (SolnBlk(SI+1).CellCenter() - SolnBlk(SI+2).CellCenter()) - 
           (SolnBlk(SI+1).CellSolution(Var) - SolnBlk(SI+2).CellSolution(Var))*
            (SolnBlk(SI).CellCenter() - SolnBlk(SI+1).CellCenter()))/

           ((SolnBlk(SI).CellCenter() - SolnBlk(SI+1).CellCenter())*
           (SolnBlk(SI+1).CellCenter() - SolnBlk(SI+2).CellCenter())*
           (SolnBlk(SI).CellCenter() - SolnBlk(SI+2).CellCenter()));
    break;

  case 3:
    return
      ( (SolnBlk(SI).CellCenter() - SolnBlk(SI+1).CellCenter())*
	(SolnBlk(SI).CellCenter() - SolnBlk(SI+2).CellCenter())*
	( (SolnBlk(SI+1).CellCenter() - SolnBlk(SI+2).CellCenter())*
	  (SolnBlk(SI+2).CellSolution(Var) - SolnBlk(SI+3).CellSolution(Var)) -
	  (SolnBlk(SI+2).CellCenter() - SolnBlk(SI+3).CellCenter())*
	  (SolnBlk(SI+1).CellSolution(Var) - SolnBlk(SI+2).CellSolution(Var)) ) +

	(SolnBlk(SI+1).CellCenter() - SolnBlk(SI+3).CellCenter())*
	(SolnBlk(SI+2).CellCenter() - SolnBlk(SI+3).CellCenter())*
	( (SolnBlk(SI+1).CellCenter() - SolnBlk(SI+2).CellCenter())*
	  (SolnBlk(SI).CellSolution(Var) - SolnBlk(SI+1).CellSolution(Var)) -
	  (SolnBlk(SI).CellCenter() - SolnBlk(SI+1).CellCenter())*
	  (SolnBlk(SI+1).CellSolution(Var) - SolnBlk(SI+2).CellSolution(Var)) )
	) /

      ((SolnBlk(SI).CellCenter() - SolnBlk(SI+1).CellCenter())*
       (SolnBlk(SI).CellCenter() - SolnBlk(SI+2).CellCenter())*
       (SolnBlk(SI).CellCenter() - SolnBlk(SI+3).CellCenter())*
       (SolnBlk(SI+1).CellCenter() - SolnBlk(SI+2).CellCenter())*
       (SolnBlk(SI+1).CellCenter() - SolnBlk(SI+3).CellCenter())*
       (SolnBlk(SI+2).CellCenter() - SolnBlk(SI+3).CellCenter())
       );

  default:
    return (DivDifference(SolnBlk,Var,Order-1,SI) - DivDifference(SolnBlk,Var,Order-1,SI+1))/ \
              (SolnBlk(SI).CellCenter() - SolnBlk(SI+Order).CellCenter());

  }
}

 /***************************************************************************
 * TEMPLATIZED Function: DataDependent_ENO_Reconstruction for 1D            *
 *                                                                          *
 * This function determines the coefficients of a Taylor series expansion   *
 * which approximates the solution over the domain of the cell specified    *
 * by "iCell" index.                                                        *
 * The method used is following closely the technique proposed by           *
 * Carl F. Ollivier-Gooch in                                                *
 * "High-order ENO schemes for unstructured meshes based on least-squares   *
 * reconstruction", AIAA Meeting, 35th, Reno, NV, Jan 6-9, 1997             *
 *                                                                          *
 ***************************************************************************/
template< class SolutionContainer> inline
void DataDependent_ENO_Reconstruction (SolutionContainer & SolnBlk, IndexType & i_index,
				       const int iCell)
{
  std::cout << "Deprecated reconstruction.\n"
	    << "Still interested in this reconstruction?\n"
	    << "Checkout code version older than 8.\n";
}

#endif
