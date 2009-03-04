#ifndef _RECONSTRUCTION_FUNCTION_3D_INCLUDED
#define _RECONSTRUCTION_FUNCTION_3D_INCLUDED

/* Include defined header file. */
#include <vector>
#include "../../../src_3D/Math/LinearSystems.h"
#include "../../../src_3D/Math/Matrix.h"
#include "../../../src_3D/Utilities/Utilities.h"
#include "include/TypeDefinition.h"
#include "Reconstruction/ReconstructionHelpers.h"

/* Investigate matrix A -> parameters */
#ifndef _RECONSTRUCTION_FUNCTION_3D_INCLUDED
static const  int InvestCellI = 0;
static const  int InvestCellJ = 0;
static const  int InvestCellK = 0;
static int FoundIndexCell = 0;
#endif

/**********************************************************************************************************************
 ----------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------------
 *                                                                                                                    *
 *                                         3D Reconstructions                                                         *
 *                                                                                                                    *
 ----------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------------
 **********************************************************************************************************************/

 /***************************************************************************
 * TEMPLATIZED Function: kExact_Reconstruction for 3D                       *
 *                                                                          *
 * This subroutine determines the coefficients of a Taylor series expansion *
 * which approximates the solution over the domain of the cell specified    *
 * by "i_index[0]", "j_index[0]", and "k_index[0]" indexes.                 *
 ***************************************************************************/
template< class SolutionContainer> inline
void kExact_Reconstruction (SolutionContainer & SolnBlk, const int *i_index, const int *j_index, const int *k_index, 
			    const int & StencilSize)
{
  typedef typename SolutionContainer::CompCellType  ComputationalCellType;
  static const int NumberOfParameters = ComputationalCellType::NumberOfVariables;

  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS
  int ND(SolnBlk.NumberOfTaylorDerivatives());	/* number of Taylor expansion coefficients */
  DenseMatrix A(StencilSize-1,ND-1);            /* the matrix which the linear system is solved for */
  DenseMatrix All_Delta_U (StencilSize-1,NumberOfParameters); /* matrix for storing U[neighbour]-U[cell] */
  DenseMatrix X (StencilSize-1,NumberOfParameters); /* matrix for storing the solution to the linear system A*X = All_Delta_U */
  ColumnVector GeomWeights(StencilSize);        /* the column vector of the geometric weights */
  Vector3D* DeltaCellCenters;                   /* array for storing the X-distance and Y-distance between the cell center
						   of neighbour cells and the one of i,j,k cell */
  int krank;                                    /* the final rank of A matrix is returned here */
  int IndexSumZ, IndexSumY, IndexSumX, P1, P2, P3;
  double CombP1X, CombP2Y, CombP3Z;
  double PowDistanceXC, PowDistanceYC, PowDistanceZC;
  int cell, i, parameter;
  double WeightsSum(0.0);
  double IntSum1(0.0),IntSum2(0.0);

  // Allocate memory
  DeltaCellCenters = new Vector3D [StencilSize];

  // *********  Assign the average solution to D00 ***********
  SolnBlk(i_index[0],j_index[0],k_index[0]).CellDeriv(0,0,0) = SolnBlk(i_index[0],j_index[0],k_index[0]).CellSolution();

  if (ND == 1){ // piecewise constant
    return;
  }

  // START:   Set the LHS and RHS of the linear system 
  // **************************************************

  // Step1. Compute the normalized geometric weights
  for (cell=1; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil

    /* Compute the X, Y, and Z component of the distance between
       the cell center of the neighbours and the reconstructed cell */
    DeltaCellCenters[cell] = SolnBlk(i_index[cell],j_index[cell],k_index[cell]).CellCenter() - 
                             SolnBlk(i_index[0],j_index[0],k_index[0]).CellCenter();

    /* Compute the geometric weights and their sum (this is used for normalization)
       based on the distance to each control volume */
    GeomWeights(cell) = sqrt(DeltaCellCenters[cell].x*DeltaCellCenters[cell].x + 
			     DeltaCellCenters[cell].y*DeltaCellCenters[cell].y +
                             DeltaCellCenters[cell].z*DeltaCellCenters[cell].z);
    GeomWeights(cell) *= GeomWeights(cell);
    GeomWeights(cell) = 1.0/(1.0E-15 + GeomWeights(cell));
    
    WeightsSum += GeomWeights(cell);
  }


  // Step2. Set the approximate equations
  for (cell=1 ; cell<StencilSize; ++cell){ //for each cell in the stencil
    
    // compute the normalized geometric weight
    GeomWeights(cell) /= WeightsSum;
    // --> RR: The following line simply gets rid of geometric weighting for the time being
    GeomWeights(cell) = 1.0;

    // *** SET the matrix A of the linear system (LHS) ***
    /* compute for each derivative the corresponding entry in the matrix of the linear system */
    for (i=1; i<=SolnBlk(i_index[0],j_index[0],k_index[0]).CellDeriv().LastElem(); ++i){
      // build the row of the matrix
      P1 = SolnBlk(i_index[0],j_index[0],k_index[0]).CellDeriv(i).P1();  // identify P1
      P2 = SolnBlk(i_index[0],j_index[0],k_index[0]).CellDeriv(i).P2();  // identify P2
      P3 = SolnBlk(i_index[0],j_index[0],k_index[0]).CellDeriv(i).P3();  // identify P3

      A(cell-1,i-1) = 0.0;  // set sumation variable to zero
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
          IntSum1 = 0.0;        // reset internal summation variable
          
          for (IndexSumX = 0; IndexSumX<=P1; ++IndexSumX){
            IntSum1 += ( CombP1X*PowDistanceXC*
                        SolnBlk(i_index[cell],j_index[cell],k_index[cell]).CellGeomCoeff(P1-IndexSumX,P2-IndexSumY,P3-IndexSumZ) );
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

      // subtract the corresponding geometric moment of (i_index[0],j_index[0],k_index[0]) cell
      A(cell-1,i-1) -= SolnBlk(i_index[0],j_index[0],k_index[0]).CellGeomCoeff(i);

      // apply geometric weighting
      A(cell-1,i-1) *= GeomWeights(cell);

    } // endfor (i) - Number of Derivatives

    // STOP:   Matrix of the linear system (LHS) built. 
    // ************************************************
      
    // *** SET the matrix All_Delta_U of the linear system (RHS) ***
    for (parameter = 1; parameter <= NumberOfParameters; ++parameter){
      All_Delta_U(cell-1,parameter-1) = (SolnBlk(i_index[cell],j_index[cell],k_index[cell]).CellSolution(parameter) -
					 SolnBlk(i_index[0],j_index[0],k_index[0]).CellSolution(parameter));
      All_Delta_U(cell-1,parameter-1) *= GeomWeights(cell);
    }

  } //endfor (cell) - Stencil Size


  // STOP:   Matrix A of the linear system (LHS) built. The same matrix is used for all variables (same geometry).
  //         Matrix All_Delta_U of the linear system (RHS) built.
  // **********************************************************************

  //Print_(A);

  /* Solve the overdetermined linear system of equations: Two METHODS available */
  /******************************************************************************/

  if(SolnBlk.UsePseudoInverse() == ON){

    // METHOD 1: Use the Pseudo Inverse of the LHS matrix
    // **************************************************
    
    // Step 1. Compute the pseudo-inverse and override the LHS term.
    // This operation will change the dimensions of the matrix.
    A.pseudo_inverse_override();

    // Step 2. Find the solution of the linear-system for the current parameter
    // Note the matrix "A" used here is really "A_inverse" via Step 1 above.
    X = A * All_Delta_U;

    // Step 3. Update the coefficients D (derivatives)
    for (i=1; i<=SolnBlk(i_index[0],j_index[0],k_index[0]).CellDeriv().LastElem(); ++i){
      for (parameter = 1; parameter <= NumberOfParameters; ++parameter){
	SolnBlk(i_index[0],j_index[0],k_index[0]).CellDeriv(i).D(parameter) = X(i-1,parameter-1);
	/* this equation makes sure that the mean conservation of each parameter is satisfied inside the reconstructed cell */
	SolnBlk(i_index[0],j_index[0],k_index[0]).CellDeriv(0).D(parameter) -= 
	  (SolnBlk(i_index[0],j_index[0],k_index[0]).CellGeomCoeff(i)*X(i-1,parameter-1));
      }
    } 
  } else {

    // METHOD 2: Use a Least-Squares procedure with the original LHS matrix
    // ********************************************************************

    // Step 1. Find the solution of the linear-system for the current parameter
    Solve_LS_Householder_F77(A, All_Delta_U, krank, NumberOfParameters, StencilSize-1, ND-1);

    // Step 2. Update the coefficients D (derivatives)
    for (i=1; i<=SolnBlk(i_index[0],j_index[0],k_index[0]).CellDeriv().LastElem(); ++i){
      for (parameter = 1; parameter <= NumberOfParameters; ++parameter){
	SolnBlk(i_index[0],j_index[0],k_index[0]).CellDeriv(i).D(parameter) = All_Delta_U(i-1,parameter-1);
	/* this equation makes sure that the mean conservation of each parameter is satisfied inside the reconstructed cell */
	SolnBlk(i_index[0],j_index[0],k_index[0]).CellDeriv(0).D(parameter) -= 
	  (SolnBlk(i_index[0],j_index[0],k_index[0]).CellGeomCoeff(i)*All_Delta_U(i-1,parameter-1));
      }
    }

  } // endif UsePseudoInverse()

  // Deallocate memory
  delete [] DeltaCellCenters; DeltaCellCenters = NULL;

}



 /***************************************************************************
 * TEMPLATIZED Function: FirstOrder_kExact_Reconstruction for 3D            *
 *                                                                          *
 * This function determines the coefficients of a Taylor series expansion   *
 * which approximates the solution over the domain of the cell specified    *
 * by "i_index[0]","j_index[0]", and "k_index[0]" indexes using a first     *
 * order polynomial function                                                *
 ***************************************************************************/
template< class SolutionContainer> inline
void FirstOrder_kExact_Reconstruction (SolutionContainer & SolnBlk, const int *i_index, 
				       const int *j_index, const int *k_index, const int & StencilSize)
{

  // Obs. (i_index[0],j_index[0],k_index[0]) --> The indexes of the cell for which the reconstruction is carried out.

  typedef typename SolutionContainer::CompCellType  ComputationalCellType;
  static const int NumberOfParameters = ComputationalCellType::NumberOfVariables;
  
  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS
  int ND(4);	                          // Number of Taylor expansion coefficients (Derivatives)
  DenseMatrix A(StencilSize-1,ND-1);      // The matrix which the linear system is solved for
  DenseMatrix All_Delta_U (StencilSize-1,NumberOfParameters); //Matrix for storing U[neighbour]-U[cell].
  ColumnVector GeomWeights(StencilSize);     // The column vector of the geometric weights
  Vector3D *DeltaCellCenters;                /* array for storing the distances between the cell center
						of neighbour cells and the one of i,j,k cell */
  int krank;
  int IndexSumX, IndexSumY, IndexSumZ, P1, P2, P3;
  double CombP1X, CombP2Y, CombP3Z;
  double PowDistanceXC, PowDistanceYC, PowDistanceZC;
  int cell, i, parameter;
  double WeightsSum(0.0);

  // Allocate memory
  DeltaCellCenters = new Vector3D [StencilSize];

  // *********  Assign the average solution to D00 ***********
  SolnBlk(i_index[0],j_index[0],k_index[0]).CellDeriv(0,0,0) = SolnBlk(i_index[0],j_index[0],k_index[0]).CellSolution();

  // START:   Set the LHS and RHS of the linear system 
  // ***************************************************

  // Step1. Compute the normalized geometric weights
  for (cell=1; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil

    /* Compute the X, Y, and Z component of the distance between
       the cell center of the neighbours and the reconstructed cell */
    DeltaCellCenters[cell] = SolnBlk(i_index[cell],j_index[cell],k_index[cell]).CellCenter() - 
      SolnBlk(i_index[0],j_index[0],k_index[0]).CellCenter();
    
    /* Compute the geometric weights and their sum (this is used for normalization)
       based on the distance to each control volume */
    GeomWeights(cell) = sqrt(DeltaCellCenters[cell].x*DeltaCellCenters[cell].x + 
			     DeltaCellCenters[cell].y*DeltaCellCenters[cell].y +
                             DeltaCellCenters[cell].z*DeltaCellCenters[cell].z);
    GeomWeights(cell) *= GeomWeights(cell);
    GeomWeights(cell) = 1.0/(1.0E-15 + GeomWeights(cell));
    
    WeightsSum += GeomWeights(cell);
  }

  // Step2. Set the approximate equations
  for (cell=1 ; cell<StencilSize; ++cell){ //for each cell in the stencil
    
    // compute the normalized geometric weight
    GeomWeights(cell) /= WeightsSum;

    // RR: --> The following simply gets rid of geometric weighting for the time being
    GeomWeights(cell) = 1.0;
    
    // *** SET the matrix A of the linear system (LHS) ***
    A(cell-1,0) = GeomWeights(cell)*DeltaCellCenters[cell].z; // D001 coefficient    
    A(cell-1,1) = GeomWeights(cell)*DeltaCellCenters[cell].y; // D010 coefficient
    A(cell-1,2) = GeomWeights(cell)*DeltaCellCenters[cell].x; // D100 coefficient
    
      
    // *** SET the matrix All_Delta_U of the linear system (RHS) ***
    for (parameter = 1; parameter <= NumberOfParameters; ++parameter){
      All_Delta_U(cell-1,parameter-1) = (SolnBlk(i_index[cell],j_index[cell],k_index[cell]).CellSolution(parameter) -
					 SolnBlk(i_index[0],j_index[0],k_index[0]).CellSolution(parameter));
      All_Delta_U(cell-1,parameter-1) *= GeomWeights(cell);
    }
  
  }//endfor (cell) - Stencil Size

  // STOP:   Matrix A of the linear system (LHS) built. The same matrix is used for all variables (same geometry).
  //         Matrix All_Delta_U of the linear system (RHS) built.
  // **********************************************************************
  
  /* Solve the overdetermined linear system of equations using a least-squares procedure*/
  /**************************************************************************************/
  Solve_LS_Householder_F77(A, All_Delta_U, krank, NumberOfParameters, StencilSize-1, ND-1);
  
  // Update the coefficients D (derivatives)
  //**************************************************
  for (parameter = 1; parameter <= NumberOfParameters; ++parameter){
    SolnBlk(i_index[0],j_index[0],k_index[0]).CellDeriv(0,0,1,parameter) = All_Delta_U(0,parameter-1);
    SolnBlk(i_index[0],j_index[0],k_index[0]).CellDeriv(0,1,0,parameter) = All_Delta_U(1,parameter-1);
    SolnBlk(i_index[0],j_index[0],k_index[0]).CellDeriv(1,0,0,parameter) = All_Delta_U(2,parameter-1);
  }
  
  // Deallocate memory
  delete [] DeltaCellCenters; DeltaCellCenters = NULL;
  
}

#endif
