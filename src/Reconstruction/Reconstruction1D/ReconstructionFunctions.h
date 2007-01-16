#ifndef _RECONSTRUCTION_FUNCTION_1D_INCLUDED
#define _RECONSTRUCTION_FUNCTION_1D_INCLUDED

/* Include defined header file. */
#include <vector>
#include "Math/LinearSystems.h"
#include "include/require.h"
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

  // Obs: The subroutine assumes that one index in the i_index[] is equal to (iCell) !!!
  // check this condition
#ifndef __No_Checking__
  {
    bool result = false;
    for(int index = 0; index<i_index.size(); ++index){
      if ((i_index[index] == iCell))
	result = true;
    }
    require(result == true, "Weighted_ENO_Reconstruction 1D Error. The cell used for reconstruction is not included in"\
	    " the stencil. The Taylor series expansion is not true for this situation. Check the stencil please!\n");
  }
#endif

  typedef typename SolutionContainer::CompCellType  ComputationalCellType;
  //  typedef typename ComputationalCellType::DerivativesContainer DerivativesContainer;
  //  typedef typename ComputationalCellType::GeometricIntegrals GeometricIntegrals;
  static const int NumberOfParameters = ComputationalCellType::NumberOfVariables;

  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS
  int RecOrder = SolnBlk(0).CellRecOrder(); /* Reconstruction Order */
  int ND = SolnBlk.NumberOfTaylorDerivatives();	/* number of Taylor expansion coefficients */
  int StencilSize = i_index.size(); /* number of cells used in the reconstruction. It includes the cell I */
  DenseMatrix A(StencilSize-1,ND-1);      // The matrix which the linear system is
  DenseMatrix A_Original;                 // solved for and its copy
  DenseMatrix All_Delta_U (StencilSize-1,NumberOfParameters); //Matrix for storing U[neighbour]-U[cell].
  DenseMatrix All_Delta_U_Original;       //Matrix for storing a copy of All_Delta_U.
  ColumnVector Delta_U(StencilSize-1);    // Vector for storing U[neighbour] - U[cell] (it's only one column of All_Delta_U)
  ColumnVector Delta_U_Original(StencilSize-1);  // Vector for storing a copy of the
                                          // original Delta_U
  ColumnVector X(ND-1);                      // The column vector of the unknowns
  ColumnVector GeomWeights(StencilSize-1);   // The column vector of the geometric weights
  ColumnVector DistanceCellCenters(StencilSize-1); /* stores the distance between the cell center of neighbour
						      cells and the I cell */
  int    krank;                 // the rank of the matrix is returned here
  double Rnorm;                 // the residual norm of the LS problem
  //  int LowerOrderControlVolume = 0; // the number of cells that didn't get the nominal reconstruction order

  double GeomIntegralOverNeighbourCell; /* represents the geometric integral over the domain of the neigbour cell */
  double DistanceToPower;

  // START:   Set the LHS and RHS of the linear system 
  // ***************************************************
  // "counter" -> refers to neighbour
  // "DerivOrder" -> refers to the position of a derivative in the Taylor expanssion

  int DerivOrder,IndexSumX;
  int P1; //int P2;
  //  int sumation = 0;
  double CombP1X;
  //  double PowDistanceYC;
  for (int cell=0, CounterNeighbCell=0 ; cell<StencilSize; ++cell){ //for each cell in the stencil
    // check if the indexes are different from those of the reconstructed cell
    if( !(i_index[cell]==iCell) ){

      // *** SET the matrix A of the linear system (LHS) ***

      /* Compute the distance between the cell centers of the neighbours and the reconstructed cell */
      DistanceCellCenters(CounterNeighbCell) = SolnBlk(i_index[cell]).CellCenter() - SolnBlk(iCell).CellCenter();

      /* compute for each derivative but D0, the corresponding entry in the matrix of the linear system */
      for (int i = 1; i<=SolnBlk(iCell).CellDeriv().LastElem(); ++i){
	// build the row of the matrix
	GeomIntegralOverNeighbourCell = 0.0;  // set sumation variable to zero
	// Compute geometric integral over the neighbour's domain
	CombP1X = 1.0; // the binomial coefficient "nC k" for k=0 is 1
	// Initialize DistanceToPower
	DistanceToPower = 1.0;
	// Compute geometric integral over the neighbour's domain
	for (IndexSumX = 0; IndexSumX<=SolnBlk(iCell).CellDeriv(i,true,true,true).P1(); ++IndexSumX){
	  GeomIntegralOverNeighbourCell += CombP1X * DistanceToPower *
	                                   SolnBlk(i_index[cell]).CellGeomCoeff()(SolnBlk(iCell).CellDeriv(i,true,true,true).P1()-
										  IndexSumX);
	    
	  // update the binomial coefficients
	  CombP1X = (SolnBlk(iCell).CellDeriv(i,true,true,true).P1()-IndexSumX)*CombP1X/(IndexSumX+1);
	  // the index is still the old one => expression for "nC k+1"

	  // Update DistanceToPower
	  DistanceToPower *= DistanceCellCenters(CounterNeighbCell);
	}

	A(CounterNeighbCell,i-1) = ( GeomIntegralOverNeighbourCell - 
				     SolnBlk(iCell).CellGeomCoeff(SolnBlk(iCell).CellDeriv(i,true,true,true).P1()) );
      }

      // *** SET the matrix All_Delta_U of the linear system (RHS) ***
      for (int parameter = 1; parameter <= NumberOfParameters; ++parameter){
	All_Delta_U(CounterNeighbCell,parameter-1) = SolnBlk( i_index[cell] ).CellSolution(parameter) -
	                                             SolnBlk(iCell).CellSolution(parameter);
      }
      ++CounterNeighbCell;
    }
  }

  // STOP:   Matrix A of the linear system (LHS) built. The same matrix is used for all the variables (same geometry).
  //         Matrix All_Delta_U of the linear system (RHS) built.
  // **********************************************************************

  /* Apply geometric weighting */
  All_Delta_U_Original = All_Delta_U;  // make copy of the original RHS
  Weighting_LS_Problem(A,All_Delta_U,DistanceCellCenters,GeomWeights);

  /* Copy the matrix A */
  A_Original = A;

  /* Solve Reconstruction for each solution variable independently.*/
  for (int parameter = 1; parameter <= NumberOfParameters; ++parameter){

    // Set the RHS column vector of the linear system.
    // *******************************************************************
    for (int IndexRow = 0; IndexRow<All_Delta_U.size(0); ++IndexRow){
      Delta_U(IndexRow) = All_Delta_U(IndexRow,parameter-1);
      Delta_U_Original(IndexRow) = All_Delta_U_Original(IndexRow,parameter-1);
    }

    // Apply data-dependent weighting
    /**
     * Method: Determine directly the Data-Dependent Weights
     */
    /* Lucian Weighting */
    /*     Weighting_LS_Problem (A,GeomWeights,DistanceCellCenters,Delta_U, */
    /*     			  Delta_U_Original, MaxDeltaSolutionStencil(parameter-1), */
    /*     			  SolnBlk.MaxDeltaSolutionDomain(parameter-1), SolnBlk.GeomCharacteristicLength(), */
    /*     			  RecOrder, SolnBlk(iCell).CellFOrder(), SolnBlk.KnobCutoff()); */
    
    /* Carl Weighting */
    /*     Carl_Weighting_LS_Problem (A,GeomWeights,DistanceCellCenters,Delta_U, */
    /* 			       Delta_U_Original, SolnBlk.MaxDeltaSolutionDomain(parameter-1), */
    /*     			       SolnBlk.GeomCharacteristicLength(), */
    /* 			       RecOrder, SolnBlk(iCell).CellFOrder(), SolnBlk.KnobCutoff()); */
    

    // Solve the Least Squares Problem
    //*************************************************
    Solve_LS_Householder(A,Delta_U,X,krank,Rnorm);

    /************************************************
     * Solution to the LS_Problem obtained
     *________________________________________________*
     *************************************************/

    // Update the coefficients D (derivatives)
    //**************************************************
    for (int i=1; i<=SolnBlk(iCell).CellDeriv().LastElem(); ++i){
      SolnBlk(iCell).CellDeriv(i,true,true,true).D(parameter) = X(i-1);

      /* this equation makes sure that the conservation of parameters is satisfied inside the cell */
      SolnBlk(iCell).CellDeriv(0,true,true,true).D(parameter) -= SolnBlk(iCell).CellGeomCoeff(i,true,true,true) * X(i-1);
    }

    /* Restore the matrix A */
    A = A_Original;
  }

}

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
			    const int iCell)
{

  typedef typename SolutionContainer::CompCellType  ComputationalCellType;
  static const int NumberOfParameters = ComputationalCellType::NumberOfVariables;

  // SET THE VARIABLES USED IN THE RECONSTRUCTION PROCESS
  int ND(SolnBlk.NumberOfTaylorDerivatives());	              /* number of Taylor expansion coefficients */
  int StencilSize(i_index.size());                            /* number of cells used in the reconstruction stencil */
  DenseMatrix  A(StencilSize-1,ND-1);                         /* the matrix which the linear system is solved for */
  DenseMatrix  All_Delta_U(StencilSize-1,NumberOfParameters); /* matrix for storing U[neighbour]-U[cell] */
  DenseMatrix  X(ND-1,NumberOfParameters);                    /* the matrix of the unknowns for the unique reconstruction */
  ColumnVector DistanceCellCenters(StencilSize-1);            /* stores the distance between the cell center of neighbour
  						                 cells and the I cell */
  ColumnVector GeomWeights(StencilSize-1);                    /* The column vector of the geometric weights */
  int krank;                                                  /* the final rank of A is returned here */
  ColumnVector Rnorm(NumberOfParameters);                     /* the residual norm of the LS problem for each parameter */

  int cell, i, IndexSumX, parameter, P1;
  double CombP1X, WeightsSum, PowDistanceXC;

  // *********  Assign the average solution to D0  ***********
  SolnBlk(iCell).CellDeriv(0,true,true,true).D( ) = SolnBlk(iCell).CellSolution( );

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
    for (i=1; i<=SolnBlk(iCell).CellDeriv().LastElem(); ++i){
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
  Solve_LS_Householder(A,All_Delta_U,X,krank,Rnorm);
  
  // Update the coefficients D (derivatives)
  //**************************************************
  for (i=1; i<=SolnBlk(iCell).CellDeriv().LastElem(); ++i){
    for (parameter = 1; parameter <= NumberOfParameters; ++parameter){
      SolnBlk(iCell).CellDeriv(i,true,true,true).D(parameter) = X(i-1,parameter-1);
      
      /* this equation makes sure that the mean conservation of each parameter is satisfied inside the reconstructed cell */
      SolnBlk(iCell).CellDeriv(0,true,true,true).D(parameter) -= ( SolnBlk(iCell).CellGeomCoeff(i,true,true,true) * 
								   X(i-1,parameter-1) );
    }
  }
}

 /***************************************************************************
 * TEMPLATIZED Function: Reconstruction_FirstOrder for 1D                   *
 *                                                                          *
 * This function determines the coefficients of a Taylor series expansion   *
 * which approximates the solution over the domain of the cell specified    *
 * by "iCell" index using a first order polynomial function.                *
 ***************************************************************************/
template< class SolutionContainer> inline
void Reconstruction_FirstOrder (SolutionContainer & SolnBlk, IndexType & i_index,
				const int iCell)
{

  // Obs: The subroutine assumes that one index in the i_index[] is equal to (iCell) !!!
  // check this condition
#ifndef __No_Checking__
  {
    bool result = false;
    for(int index = 0; index<i_index.size(); ++index){
      if ((i_index[index] == iCell))
	result = true;
    }
    require(result == true, "Reconstruction_FirstOrder 1D Error. The cell used for reconstruction is not included in"\
	    " the stencil. The Taylor series expansion is not true for this situation. Check the stencil please!\n");
  }
#endif

  /* Initialize iterators for the containers of the Taylor derivatives
     and the Geometric Coefficients (geometric integrals).*/
  typedef typename SolutionContainer::CompCellType  ComputationalCellType;
  static const int NumberOfParameters = ComputationalCellType::NumberOfVariables;

  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS
  int RecOrder = 1; /* Reconstruction Order */
  int ND = 2;	    /* number of Taylor expansion coefficients */
  int StencilSize = i_index.size(); /* number of cells used in the reconstruction. It includes the cell I */
  DenseMatrix A(StencilSize-1,ND-1);      // The matrix which the linear system is
  DenseMatrix A_Original;                 // solved for and its copy
  DenseMatrix All_Delta_U (StencilSize-1,NumberOfParameters); //Matrix for storing U[neighbour]-U[cell].
  DenseMatrix All_Delta_U_Original;       //Matrix for storing a copy of All_Delta_U.
  ColumnVector Delta_U(StencilSize-1);    // Vector for storing U[neighbour] - U[cell] (it's only one column of All_Delta_U)
  ColumnVector Delta_U_Original(StencilSize-1);  // Vector for storing a copy of the
                                          // original Delta_U
  ColumnVector X(ND-1);                      // The column vector of the unknowns
  ColumnVector GeomWeights(StencilSize-1);   // The column vector of the geometric weights
  ColumnVector DistanceCellCenters(StencilSize-1); /* stores the distance between the cell center of neighbour
						      cells and the I cell */
  int    krank;                 // the rank of the matrix is returned here
  double Rnorm;                 // the residual norm of the LS problem
  double GeomIntegralOverNeighbourCell; /* represents the geometric integral over the domain of the neigbour cell */
  double DistanceToPower;

  // START:   Set the LHS and RHS of the linear system 
  // ***************************************************
  // "counter" -> refers to neighbour
  // "DerivOrder" -> refers to the position of a derivative in the Taylor expanssion

  int DerivOrder,IndexSumX;
  int P1;
  double CombP1X;

  for (int cell=0, CounterNeighbCell=0 ; cell<StencilSize; ++cell){ //for each cell in the stencil
    // check if the indexes are different from those of the reconstructed cell
    if( !(i_index[cell]==iCell) ){

      // *** SET the matrix A of the linear system (LHS) ***

      /* Compute the distance between the cell centers of the neighbours and the reconstructed cell */
      DistanceCellCenters(CounterNeighbCell) = SolnBlk(i_index[cell]).CellCenter() - SolnBlk(iCell).CellCenter();

      /* compute for each derivative but D0, the corresponding entry in the matrix of the linear system */
      for (int i = 1; i<=SolnBlk(iCell).CellDerivFirstOrder().LastElem(); ++i){
	// build the row of the matrix
	GeomIntegralOverNeighbourCell = 0.0;  // set sumation variable to zero
	// Compute geometric integral over the neighbour's domain
	CombP1X = 1.0; // the binomial coefficient "nC k" for k=0 is 1
	// Initialize DistanceToPower
	DistanceToPower = 1.0;
	// Compute geometric integral over the neighbour's domain
	for (IndexSumX = 0; IndexSumX<=SolnBlk(iCell).CellDerivFirstOrder(i,true,true,true).P1(); ++IndexSumX){
	  GeomIntegralOverNeighbourCell += CombP1X * DistanceToPower *
	                                   SolnBlk(i_index[cell]).CellGeomCoeff()(SolnBlk(iCell).CellDerivFirstOrder(i,true,true,true).P1() - IndexSumX);
	    
	  // update the binomial coefficients
	  CombP1X = (SolnBlk(iCell).CellDerivFirstOrder(i,true,true,true).P1()-IndexSumX)*CombP1X/(IndexSumX+1);
	  // the index is still the old one => expression for "nC k+1"

	  // Update DistanceToPower
	  DistanceToPower *= DistanceCellCenters(CounterNeighbCell);
	}

	A(CounterNeighbCell,i-1) = ( GeomIntegralOverNeighbourCell - 
				     SolnBlk(iCell).CellGeomCoeff(SolnBlk(iCell).CellDerivFirstOrder(i,true,true,true).P1()) );
      }

      // *** SET the matrix All_Delta_U of the linear system (RHS) ***
      for (int parameter = 1; parameter <= NumberOfParameters; ++parameter){
	All_Delta_U(CounterNeighbCell,parameter-1) = SolnBlk( i_index[cell] ).CellSolution(parameter) -
	                                             SolnBlk(iCell).CellSolution(parameter);
      }
      ++CounterNeighbCell;
    }
  }

  // STOP:   Matrix A of the linear system (LHS) built. The same matrix is used for all the variables (same geometry).
  //         Matrix All_Delta_U of the linear system (RHS) built.
  // **********************************************************************

  /* Apply geometric weighting */
  All_Delta_U_Original = All_Delta_U;  // make copy of the original RHS
  Weighting_LS_Problem(A,All_Delta_U,DistanceCellCenters,GeomWeights);

  /* Copy the matrix A */
  A_Original = A;

  /* Solve Reconstruction for each solution variable independently.*/
  for (int parameter = 1; parameter <= NumberOfParameters; ++parameter){

    // Set the RHS column vector of the linear system.
    // *******************************************************************
    for (int IndexRow = 0; IndexRow<All_Delta_U.size(0); ++IndexRow){
      Delta_U(IndexRow) = All_Delta_U(IndexRow,parameter-1);
      Delta_U_Original(IndexRow) = All_Delta_U_Original(IndexRow,parameter-1);
    }

    // Solve the Least Squares Problem
    //*************************************************
    Solve_LS_Householder(A,Delta_U,X,krank,Rnorm);

    /************************************************
     * Solution to the LS_Problem obtained
     *________________________________________________*
     *************************************************/

    // Update the coefficients D (derivatives)
    //**************************************************
    for (int i=1; i<=SolnBlk(iCell).CellDerivFirstOrder().LastElem(); ++i){
      SolnBlk(iCell).CellDerivFirstOrder(i,true,true,true).D(parameter) = X(i-1);

      /* this equation makes sure that the conservation of parameters is satisfied inside the cell */
      SolnBlk(iCell).CellDerivFirstOrder(0,true,true,true).D(parameter) -= SolnBlk(iCell).CellGeomCoeff(SolnBlk(iCell).CellDerivFirstOrder(i,true,true,true).P1()) * X(i-1);
    }

    /* Restore the matrix A */
    A = A_Original;
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

#endif
