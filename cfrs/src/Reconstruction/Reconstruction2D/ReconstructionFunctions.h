#ifndef _RECONSTRUCTION_FUNCTION_2D_INCLUDED
#define _RECONSTRUCTION_FUNCTION_2D_INCLUDED

/* Include defined header file. */
#include <vector>
#include "../../../src_2D/Math/LinearSystems.h"
#include "../../../src_2D/Utilities/Utilities.h"
#include "include/TypeDefinition.h"
#include "Reconstruction/ReconstructionHelpers.h"

/* Investigate matrix A -> parameters */
static const  int InvestCellI = 0;
static const  int InvestCellJ = 0;
static int FoundIndexCell = 0;

/**********************************************************************************************************************
 ----------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------------
 *                                                                                                                    *
 *                                         2D Reconstructions                                                         *
 *                                                                                                                    *
 ----------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------------
 **********************************************************************************************************************/

 /***************************************************************************
 * TEMPLATIZED Function: kExact_Reconstruction for 2D                       *
 *                                                                          *
 * This subroutine determines the coefficients of a Taylor series expansion *
 * which approximates the solution over the domain of the cell specified    *
 * by "i_index[0]" and "j_index[0]" indexes.                                *
 ***************************************************************************/
template< class SolutionContainer> inline
void kExact_Reconstruction (SolutionContainer & SolnBlk, const int *i_index, const int *j_index, const int & StencilSize)
{
  typedef typename SolutionContainer::CompCellType  ComputationalCellType;
  static const int NumberOfParameters = ComputationalCellType::NumberOfVariables;

  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS
  int ND(SolnBlk.NumberOfTaylorDerivatives());	/* number of Taylor expansion coefficients */
  DenseMatrix A(StencilSize-1,ND-1);            /* the matrix which the linear system is solved for */
  DenseMatrix All_Delta_U (StencilSize-1,NumberOfParameters); /* matrix for storing U[neighbour]-U[cell] */
  DenseMatrix X (StencilSize-1,NumberOfParameters); /* matrix for storing the solution to the linear system A*X = All_Delta_U */
  ColumnVector GeomWeights(StencilSize);        /* the column vector of the geometric weights */
  Vector2D* DeltaCellCenters;                   /* array for storing the X-distance and Y-distance between the cell center
						   of neighbour cells and the one of i,j cell */
  int krank;                                      /* the final rank of A matrix is returned here */
  int IndexSumY, IndexSumX, P1, P2;
  double CombP1X, CombP2Y;
  double PowDistanceYC, PowDistanceXC;
  int cell, i, parameter;
  double WeightsSum(0.0);
  double IntSum(0.0);

  // Allocate memory
  DeltaCellCenters = new Vector2D [StencilSize];

  // *********  Assign the average solution to D00 ***********
  SolnBlk(i_index[0],j_index[0]).CellDeriv(0,0) = SolnBlk(i_index[0],j_index[0]).CellSolution();

  if (ND == 1){ // piecewise constant
    return;
  }

  // START:   Set the LHS and RHS of the linear system 
  // ***************************************************

  // Step1. Compute the normalized geometric weights
  for (cell=1; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil

    /* Compute the X and Y component of the distance between
       the cell center of the neighbours and the reconstructed cell */
    DeltaCellCenters[cell] = SolnBlk(i_index[cell],j_index[cell]).CellCenter() - SolnBlk(i_index[0],j_index[0]).CellCenter();
    
    /* Compute the geometric weights and their sum (this is used for normalization)
       based on the distance to each control volume */
    GeomWeights(cell) = sqrt(DeltaCellCenters[cell].x*DeltaCellCenters[cell].x + 
			     DeltaCellCenters[cell].y*DeltaCellCenters[cell].y);
    GeomWeights(cell) *= GeomWeights(cell);
    GeomWeights(cell) = 1.0/(1.0E-15 + GeomWeights(cell));
    
    WeightsSum += GeomWeights(cell);
  }

  // Step2. Set the approximate equations
  for (cell=1 ; cell<StencilSize; ++cell){ //for each cell in the stencil
    
    // compute the normalized geometric weight
    GeomWeights(cell) /= WeightsSum;

    // *** SET the matrix A of the linear system (LHS) ***
    /* compute for each derivative the corresponding entry in the matrix of the linear system */
    for (i=1; i<=SolnBlk(i_index[0],j_index[0]).CellDeriv().LastElem(); ++i){
      // build the row of the matrix
      P1 = SolnBlk(i_index[0],j_index[0]).CellDeriv(i,true,true,true).P1();  // identify P1
      P2 = SolnBlk(i_index[0],j_index[0]).CellDeriv(i,true,true,true).P2();  // identify P2
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
		      SolnBlk(i_index[cell],j_index[cell]).CellGeomCoeff(P1-IndexSumX,P2-IndexSumY) );
	    
	  // update the binomial coefficients
	  CombP1X = (P1-IndexSumX)*CombP1X/(IndexSumX+1); // the index is still the old one => expression for "nC k+1"
	  PowDistanceXC *= DeltaCellCenters[cell].x;      // Update PowDistanceXC
	}//endfor

	A(cell-1,i-1) += CombP2Y*PowDistanceYC*IntSum; // update the external sum

	CombP2Y = (P2-IndexSumY)*CombP2Y/(IndexSumY+1); // the index is still the old one => expression for "nC k+1"
	PowDistanceYC *= DeltaCellCenters[cell].y;    // Update PowDistanceYC
      }//endfor

      // subtract the corresponding geometric moment of (i_index[0],j_index[0]) cell
      A(cell-1,i-1) -= SolnBlk(i_index[0],j_index[0]).CellGeomCoeff(i,true,true,true);

      // apply geometric weighting
      A(cell-1,i-1) *= GeomWeights(cell);
    }
      
    // *** SET the matrix All_Delta_U of the linear system (RHS) ***
    for (parameter = 1; parameter <= NumberOfParameters; ++parameter){
      All_Delta_U(cell-1,parameter-1) = (SolnBlk(i_index[cell],j_index[cell]).CellSolution(parameter) -
					 SolnBlk(i_index[0],j_index[0]).CellSolution(parameter));
      All_Delta_U(cell-1,parameter-1) *= GeomWeights(cell);
    }

  }//endfor (cell)

  // STOP:   Matrix A of the linear system (LHS) built. The same matrix is used for all variables (same geometry).
  //         Matrix All_Delta_U of the linear system (RHS) built.
  // **********************************************************************

  /* Solve the overdetermined linear system of equations: Two METHODS available */
  /******************************************************************************/

  if(CENO_CFRS_Execution_Mode::USE_PSEUDO_INVERSE == ON){

    // METHOD 1: Use the Pseudo Inverse of the LHS matrix
    // **************************************************
    
    // Ensure that the LHS matrix is formated correctly.
    // Memory shouldn't be allocated here, only the dimensions should be defined properly.
    SolnBlk.Cell_LHS_Inv(i_index[0],j_index[0]).newsize(StencilSize - 1, ND - 1);

    // Step 1. Compute and store the pseudo-inverse 
    // This operation will change the dimensions of the matrix and override the LHS term.
    SolnBlk.Cell_LHS_Inv(i_index[0],j_index[0]) = A;    
    SolnBlk.Cell_LHS_Inv(i_index[0],j_index[0]).pseudo_inverse_override();

     // Step 2. Find the solution of the linear-system for the current parameter
    // Note the matrix "A" used here is really "A_inverse" via Step 1 above.
    X = SolnBlk.Cell_LHS_Inv(i_index[0],j_index[0]) * All_Delta_U;

//   if (i_index[0] == 5 && j_index[0] == 5){
//      Print_(SolnBlk.Cell_LHS_Inv(i_index[0],j_index[0]));
//      Print_(X);
//      Print_(All_Delta_U);
//    }

    // Step 3. Update the coefficients D (derivatives)
    for (i=1; i<=SolnBlk(i_index[0],j_index[0]).CellDeriv().LastElem(); ++i){
      for (parameter = 1; parameter <= NumberOfParameters; ++parameter){
	SolnBlk(i_index[0],j_index[0]).CellDeriv(i,true,true,true).D(parameter) = X(i-1,parameter-1);
	/* this equation makes sure that the mean conservation of each parameter is satisfied inside the reconstructed cell */
	SolnBlk(i_index[0],j_index[0]).CellDeriv(0,true,true,true).D(parameter) -= 
	  (SolnBlk(i_index[0],j_index[0]).CellGeomCoeff(i,true,true,true)*X(i-1,parameter-1));
      }
    } 
  } else {

    // METHOD 2: Use a Least-Squares procedure with the original LHS matrix
    // ********************************************************************

    // Step 1. Find the solution of the linear-system for the current parameter
    Solve_LS_Householder_F77(A, All_Delta_U, krank, NumberOfParameters, StencilSize-1, ND-1);

    // Step 2. Update the coefficients D (derivatives)
    for (i=1; i<=SolnBlk(i_index[0],j_index[0]).CellDeriv().LastElem(); ++i){
      for (parameter = 1; parameter <= NumberOfParameters; ++parameter){
	SolnBlk(i_index[0],j_index[0]).CellDeriv(i,true,true,true).D(parameter) = All_Delta_U(i-1,parameter-1);
	/* this equation makes sure that the mean conservation of each parameter is satisfied inside the reconstructed cell */
	SolnBlk(i_index[0],j_index[0]).CellDeriv(0,true,true,true).D(parameter) -= 
	  (SolnBlk(i_index[0],j_index[0]).CellGeomCoeff(i,true,true,true)*All_Delta_U(i-1,parameter-1));
      }
    }

  } // endif UsePseudoInverse()
  
  // Deallocate memory
  delete [] DeltaCellCenters; DeltaCellCenters = NULL;

}

 /***************************************************************************
 * TEMPLATIZED Function: FirstOrder_kExact_Reconstruction for 2D            *
 *                                                                          *
 * This function determines the coefficients of a Taylor series expansion   *
 * which approximates the solution over the domain of the cell specified    *
 * by "i_index[0]" and "j_index[0]" indexes using a first order polynomial  *
 * function.                                                                *
 ***************************************************************************/
template< class SolutionContainer> inline
void FirstOrder_kExact_Reconstruction (SolutionContainer & SolnBlk, const int *i_index, 
				       const int *j_index, const int & StencilSize)
{

  // Obs. (i_index[0],j_index[0]) --> The indexes of the cell for which the reconstruction is carried out.

  typedef typename SolutionContainer::CompCellType  ComputationalCellType;
  static const int NumberOfParameters = ComputationalCellType::NumberOfVariables;

  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS
  int ND(3);	                          // Number of Taylor expansion coefficients (Derivatives)
  DenseMatrix A(StencilSize-1,ND-1);      // The matrix which the linear system is solved for
  DenseMatrix All_Delta_U (StencilSize-1,NumberOfParameters); //Matrix for storing U[neighbour]-U[cell].
  ColumnVector GeomWeights(StencilSize);     // The column vector of the geometric weights
  Vector2D *DeltaCellCenters;                /* array for storing the X-distance and Y-distance between the cell center
						of neighbour cells and the one of i,j cell */        
  int krank;
  int IndexSumY, IndexSumX, P1, P2;
  double CombP1X, CombP2Y;
  double PowDistanceYC, PowDistanceXC;
  int cell, i, parameter;
  double WeightsSum(0.0);

  // Allocate memory
  DeltaCellCenters = new Vector2D [StencilSize];

  // *********  Assign the average solution to D00 ***********
  SolnBlk(i_index[0],j_index[0]).CellDeriv(0,true,true,true) = SolnBlk(i_index[0],j_index[0]).CellSolution();

  // START:   Set the LHS and RHS of the linear system 
  // ***************************************************

  // Step1. Compute the normalized geometric weights
  for (cell=1; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil

    /* Compute the X and Y component of the distance between
       the cell center of the neighbours and the reconstructed cell */
    DeltaCellCenters[cell] = SolnBlk(i_index[cell],j_index[cell]).CellCenter() - SolnBlk(i_index[0],j_index[0]).CellCenter();
    
    /* Compute the geometric weights and their sum (this is used for normalization)
       based on the distance to each control volume */
    GeomWeights(cell) = sqrt(DeltaCellCenters[cell].x*DeltaCellCenters[cell].x + 
			     DeltaCellCenters[cell].y*DeltaCellCenters[cell].y);
    GeomWeights(cell) *= GeomWeights(cell);
    GeomWeights(cell) = 1.0/(1.0E-15 + GeomWeights(cell));
    
    WeightsSum += GeomWeights(cell);
  }

  // Step2. Set the approximate equations
  for (cell=1 ; cell<StencilSize; ++cell){ //for each cell in the stencil
    
    // compute the normalized geometric weight
    GeomWeights(cell) /= WeightsSum;
    
    // *** SET the matrix A of the linear system (LHS) ***
    A(cell-1,0) = GeomWeights(cell)*DeltaCellCenters[cell].y; // D01 coefficient
    A(cell-1,1) = GeomWeights(cell)*DeltaCellCenters[cell].x; // D10 coefficient
      
    // *** SET the matrix All_Delta_U of the linear system (RHS) ***
    for (parameter = 1; parameter <= NumberOfParameters; ++parameter){
      All_Delta_U(cell-1,parameter-1) = (SolnBlk(i_index[cell],j_index[cell]).CellSolution(parameter) -
					 SolnBlk(i_index[0],j_index[0]).CellSolution(parameter));
      All_Delta_U(cell-1,parameter-1) *= GeomWeights(cell);
    }
  
  }//endfor (cell)

  // STOP:   Matrix A of the linear system (LHS) built. The same matrix is used for all variables (same geometry).
  //         Matrix All_Delta_U of the linear system (RHS) built.
  // **********************************************************************
  
  /* Solve the overdetermined linear system of equations using a least-squares procedure*/
  /**************************************************************************************/
  Solve_LS_Householder_F77(A, All_Delta_U, krank, NumberOfParameters, StencilSize-1, ND-1);

  // Update the coefficients D (derivatives)
  //**************************************************
  for (parameter = 1; parameter <= NumberOfParameters; ++parameter){
    SolnBlk(i_index[0],j_index[0]).CellDeriv(0,1,parameter) = All_Delta_U(0,parameter-1);
    SolnBlk(i_index[0],j_index[0]).CellDeriv(1,0,parameter) = All_Delta_U(1,parameter-1);
  }

  // Deallocate memory
  delete [] DeltaCellCenters; DeltaCellCenters = NULL;

}

 /***************************************************************************
 * TEMPLATIZED Function: DataDependent_ENO_Reconstruction for 2D            *
 *                                                                          *
 * This function determines the coefficients of a Taylor series expansion   *
 * which approximates the solution over the domain of the cell specified    *
 * by "iCell" and "jCell" indexes.                                          *
 * The method used is following closely the technique proposed by           *
 * Carl F. Ollivier-Gooch in                                                *
 * "High-order ENO schemes for unstructured meshes based on least-squares   *
 * reconstruction", AIAA Meeting, 35th, Reno, NV, Jan 6-9, 1997             *
 *                                                                          *
 *                                                                          *
 ***************************************************************************/

template< class SolutionContainer> inline
void DataDependent_ENO_Reconstruction (SolutionContainer & SolnBlk, const int *i_index, const int *j_index,
				       const int & StencilSize,
				       const int iCell, const int jCell)
{
  std::cout << "Deprecated reconstruction.\n"
	    << "Still interested in this reconstruction?\n"
	    << "Checkout code version older than 8.\n";
}

#endif
