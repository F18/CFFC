#ifndef _RECONSTRUCTION_FUNCTION_2D_INCLUDED
#define _RECONSTRUCTION_FUNCTION_2D_INCLUDED

/* Include defined header file. */
#include <vector>
#include "Math/LinearSystems.h"
#include "include/require.h"
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
#define CentralStencil

  /* Initialize iterators for the containers of the Taylor derivatives
     and the Geometric Coefficients (geometric integrals).*/
  typedef typename SolutionContainer::CompCellType  ComputationalCellType;
  typedef typename ComputationalCellType::DerivativesContainer DerivativesContainer;
  typedef typename ComputationalCellType::GeometricIntegrals GeometricIntegrals;

  static const int NumberOfParameters = ComputationalCellType::NumberOfVariables;

  // SET VARIABLES USED IN THE RECONSTRUCTION PROCESS

  int RecOrder = SolnBlk(0,0).CellRecOrder(); /* Reconstruction Order */
  int ND = SolnBlk.NumberOfTaylorDerivatives();	/* number of Taylor expansion coefficients */
  DenseMatrix A(StencilSize-1,ND-1);      // The matrix which the linear system is solved for
  DenseMatrix All_Delta_U (StencilSize-1,NumberOfParameters); //Matrix for storing U[neighbour]-U[cell].

#ifndef CentralStencil
  DenseMatrix A_Original(StencilSize-1,ND-1);
  DenseMatrix All_Delta_U_Original(StencilSize-1,NumberOfParameters); //Matrix for storing a copy of All_Delta_U.
  ColumnVector Delta_U(StencilSize-1);  // Vector for storing U[neighbour] - U[cell] (it's only one column of All_Delta_U)
  ColumnVector Delta_U_Original(StencilSize-1); // Vector for storing a copy of the
                                                // original Delta_U
  ColumnVector MaxDeltaSolutionStencil(NumberOfParameters); //The maximum delta solution for each parameter
  ColumnVector MinSolutionStencil(NumberOfParameters);      //The minimum solution for each parameter
  double MaxSolution, MinSolution;  // The maximum and minimum solution values in the stencil
  vector<int> iCellMaxSolution, jCellMaxSolution, iCellMinSolution, jCellMinSolution; // The indeces of the cells having the
                                                                        //max and min solutions
  ColumnVector X(ND-1);                      // The column vector of the unknowns
  double Rnorm;                 // the residual norm of the LS problem
#else
  DenseMatrix X(ND-1,NumberOfParameters);    // The matrix of the unknowns
  ColumnVector Rnorm(NumberOfParameters);    // The residual norm of the LS problem for each parameter
#endif // CentralStencil

  ColumnVector GeomWeights(StencilSize-1);   // The column vector of the geometric weights
  vector<Vector2D> DeltaCellCenters; /* stores the difference on x and y between the cell center of neighbour cells
					and the one of i,j cell */
  ColumnVector DistanceCellCenters(StencilSize-1); /* stores the distance between the cell center of neighbour
						      cells and the i,j cell */
  int    krank;                 // the rank of the matrix is returned here
  double GeomIntegralOverNeighbourCell; /* represents the geometric integral over the domain of the neigbour cell */

  int IndexSumY, IndexSumX, P1, P2;
  double CombP1X, CombP2Y;
  double PowDistanceYC, PowDistanceXC;

  // START:   Set the LHS and RHS of the linear system 
  // ***************************************************
  // "counter" -> refers to neighbour
  // "DerivOrder" -> refers to the position of a derivative in the Taylor expanssion

  /* Reserve memory for DeltaCellCenters */
  DeltaCellCenters.reserve(StencilSize-1);

  for (int cell=0, CounterNeighbCell=0 ; cell<StencilSize; ++cell){ //for each cell in the stencil
    // check if the indexes are different from those of the reconstructed cell
    if( !((i_index[cell]==iCell) && (j_index[cell]==jCell)) ){

      // *** SET the matrix A of the linear system (LHS) ***
      {
	/* Compute the distance between the cell centers of the neighbours and the reconstructed cell */
	DeltaCellCenters.push_back(SolnBlk(i_index[cell],j_index[cell]).CellCenter() - SolnBlk(iCell,jCell).CellCenter());
	
	/* compute for each derivative the corresponding entry in the matrix of the linear system */
	for (int i=1; i<=SolnBlk(iCell,jCell).CellDeriv().LastElem(); ++i){
	  // build the row of the matrix
	  P1 = SolnBlk(iCell,jCell).CellDeriv(i,true,true,true).P1();  // identify P1
	  P2 = SolnBlk(iCell,jCell).CellDeriv(i,true,true,true).P2();  // identify P2
	  GeomIntegralOverNeighbourCell = 0.0;  // set sumation variable to zero
	  CombP2Y = 1.0; // the binomial coefficient "nC k" for k=0 is 1
	  PowDistanceYC = 1.0; 	  // Initialize PowDistanceYC
	  // Compute geometric integral over the neighbour's domain
	  for (IndexSumY = 0; IndexSumY<=P2; ++IndexSumY){
	    CombP1X = 1.0; // the binomial coefficient "nC k" for k=0 is 1
	    PowDistanceXC = 1.0; 	    // Initialize PowDistanceXC

	    for (IndexSumX = 0; IndexSumX<=P1; ++IndexSumX){
	      GeomIntegralOverNeighbourCell += CombP2Y*CombP1X*PowDistanceXC*
		                               (PowDistanceYC*
						SolnBlk(i_index[cell],j_index[cell]).CellGeomCoeff()(P1-IndexSumX,P2-IndexSumY));
	    
	      // update the binomial coefficients
	      CombP1X = (P1-IndexSumX)*CombP1X/(IndexSumX+1); // the index is still the old one => expression for "nC k+1"
	      PowDistanceXC *= DeltaCellCenters[CounterNeighbCell].x;	      // Update PowDistanceXC
	    }//endfor
	    CombP2Y = (P2-IndexSumY)*CombP2Y/(IndexSumY+1); // the index is still the old one => expression for "nC k+1"
	    PowDistanceYC *= DeltaCellCenters[CounterNeighbCell].y; 	    // Update PowDistanceYC
	  }//endfor
	  A(CounterNeighbCell,i-1) = GeomIntegralOverNeighbourCell - SolnBlk(iCell,jCell).CellGeomCoeff()(P1,P2);
	}
      }

      // *** SET the matrix All_Delta_U of the linear system (RHS) ***
      for (int parameter = 1; parameter <= NumberOfParameters; ++parameter){
	All_Delta_U(CounterNeighbCell,parameter-1) = SolnBlk( i_index[cell], j_index[cell] ).CellSolution(parameter) -
	                                             SolnBlk(iCell,jCell).CellSolution(parameter);

#if 0
	/* Compare the Min and Max Solution Stencil against the values of the current cell for each parameter */
	MaxSolution = MaxDeltaSolutionStencil(parameter-1);
	MaxDeltaSolutionStencil(parameter-1) = max(SolnBlk( i_index[cell], j_index[cell] ).CellSolution(parameter),
						   MaxDeltaSolutionStencil(parameter-1));
	if(MaxSolution != MaxDeltaSolutionStencil(parameter-1)){
	  iCellMaxSolution[parameter-1] = i_index[cell];
	  jCellMaxSolution[parameter-1] = j_index[cell];
	}
	MinSolution = MinSolutionStencil(parameter-1);
	MinSolutionStencil(parameter-1) = min(SolnBlk( i_index[cell], j_index[cell] ).CellSolution(parameter),
					      MinSolutionStencil(parameter-1));
	if(MinSolution != MinSolutionStencil(parameter-1)){
	  iCellMinSolution[parameter-1] = i_index[cell];
	  jCellMinSolution[parameter-1] = j_index[cell];
	}
#endif

      }//endfor (parameter)
      ++CounterNeighbCell;
    }//endif
  }//endfor (cell)

  // STOP:   Matrix A of the linear system (LHS) built. The same matrix is used for all the variables (same geometry).
  //         Matrix All_Delta_U of the linear system (RHS) built.
  // **********************************************************************

#if 0
  /* Determine the maximum delta solution for each parameter.*/
  DetermineMaxDeltaSolutionStencil(MaxDeltaSolutionStencil,MinSolutionStencil,NumberOfParameters);
  /* Try using max(All_Delta_U) for MaxDeltaSolutionStencil  */
  //DetermineMaxDeltaSolutionStencil(MaxDeltaSolutionStencil,All_Delta_U,NumberOfParameters);
#endif

  /* Determine the distance between the cell centers */
  DetermineDistanceCellCenters(DistanceCellCenters,DeltaCellCenters);




#ifdef CentralStencil
  /* This is the fast reconstruction -> uses the same geometric matrix (Left Hand Side) for all variables */

  /* Apply geometric weighting */
  Weighting_LS_Problem(A,All_Delta_U,DistanceCellCenters,GeomWeights);

  // Solve the Least Squares Problem
  //*************************************************
  Solve_LS_Householder(A,All_Delta_U,X,krank,Rnorm);

  // Update the coefficients D (derivatives)
  //**************************************************
  int parameter;
  for (int i=1; i<=SolnBlk(iCell,jCell).CellDeriv().LastElem(); ++i){
    for (parameter = 1; parameter <= NumberOfParameters; ++parameter){
      SolnBlk(iCell,jCell).CellDeriv(i,true,true,true).D(parameter) = X(i-1,parameter-1);
    
      /* this equation makes sure that the conservation of parameters is satisfied inside the cell */
      SolnBlk(iCell,jCell).CellDeriv(0,true,true,true).D(parameter) -= SolnBlk(iCell,jCell).CellGeomCoeff(i,true,true,true) * X(i-1,parameter-1);
    }
  }

#else 
  /* This is the slow reconstruction -> uses different geometric matrix (Left Hand Side) for each variable */

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

#if 0
    // Apply data-dependent weighting
    /**
     * Method: Determine directly the Data-Dependent Weights
     */
    Weighting_LS_Problem (A,GeomWeights,DistanceCellCenters,Delta_U,
    			  Delta_U_Original, MaxDeltaSolutionStencil(parameter-1),
    			  SolnBlk.MaxDeltaSolutionDomain(parameter-1), SolnBlk.GeomCharacteristicLength(),
    			  RecOrder, SolnBlk(iCell,jCell).CellFOrder(), SolnBlk.KnobCutoff());
    
    
    Carl_Weighting_LS_Problem (A,GeomWeights,DistanceCellCenters,Delta_U,
    			       Delta_U_Original, SolnBlk.MaxDeltaSolutionDomain(parameter-1),
			       SolnBlk.GeomCharacteristicLength(),
    			       RecOrder, SolnBlk(iCell,jCell).CellFOrder(), SolnBlk.KnobCutoff());
    
#endif

    // Solve the Least Squares Problem
    //*************************************************
    Solve_LS_Householder(A,Delta_U,X,krank,Rnorm);

#if 0
    if ((krank!=RecOrder) && (MaxDeltaSolutionStencil(parameter-1) > 0)){
      SolnBlk(iCell,jCell).CellFOrder() = krank;
    } else {
      SolnBlk(iCell,jCell).CellFOrder() = RecOrder;
    }
#endif

    /************************************************
     * Solution to the LS_Problem obtained
     *________________________________________________*
     *************************************************/
    // Update the coefficients D (derivatives)
    //**************************************************
    for (int i=1; i<=SolnBlk(iCell,jCell).CellDeriv().LastElem(); ++i){
      SolnBlk(iCell,jCell).CellDeriv(i,true,true,true).D(parameter) = X(i-1);

      /* this equation makes sure that the conservation of parameters is satisfied inside the cell */
      SolnBlk(iCell,jCell).CellDeriv(0,true,true,true).D(parameter) -= SolnBlk(iCell,jCell).CellGeomCoeff(i,true,true,true) * X(i-1);
    }

    /* Restore the matrix A */
    A = A_Original;
  }

#endif //CentralStencil

}

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
  DenseMatrix X(ND-1,NumberOfParameters);       /* the matrix of the unknowns */
  ColumnVector Rnorm(NumberOfParameters);       /* the residual norm of the LS problem for each parameter */
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
  
  /* Solve the overdetermined linear system of equations using a least-squares procedure*/
  /**************************************************************************************/
  Solve_LS_Householder(A,All_Delta_U,X,krank,Rnorm);

  // Update the coefficients D (derivatives)
  //**************************************************
  for (i=1; i<=SolnBlk(i_index[0],j_index[0]).CellDeriv().LastElem(); ++i){
    for (parameter = 1; parameter <= NumberOfParameters; ++parameter){
      SolnBlk(i_index[0],j_index[0]).CellDeriv(i,true,true,true).D(parameter) = X(i-1,parameter-1);
    
      /* this equation makes sure that the mean conservation of each parameter is satisfied inside the reconstructed cell */
      SolnBlk(i_index[0],j_index[0]).CellDeriv(0,true,true,true).D(parameter) -= (
			                                           SolnBlk(i_index[0],j_index[0]).CellGeomCoeff(i,true,true,true) * 
								   X(i-1,parameter-1));
    }
  }

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
  DenseMatrix X(ND-1,NumberOfParameters);    // The matrix of the unknowns
  ColumnVector Rnorm(NumberOfParameters);    // The residual norm of the LS problem for each parameter
  ColumnVector GeomWeights(StencilSize);     // The column vector of the geometric weights
  Vector2D *DeltaCellCenters;                /* array for storing the X-distance and Y-distance between the cell center
						of neighbour cells and the one of i,j cell */        
  int krank;
  int IndexSumY, IndexSumX, P1, P2;
  double CombP1X, CombP2Y;
  double PowDistanceYC, PowDistanceXC;
  int cell, i, parameter;
  double WeightsSum(0.0);
  double IntSum(0.0);

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
  Solve_LS_Householder(A,All_Delta_U,X,krank,Rnorm);

  // Update the coefficients D (derivatives)
  //**************************************************
  for (parameter = 1; parameter <= NumberOfParameters; ++parameter){
    SolnBlk(i_index[0],j_index[0]).CellDeriv(0,1,parameter) = X(0,parameter-1);
    SolnBlk(i_index[0],j_index[0]).CellDeriv(1,0,parameter) = X(1,parameter-1);
  }

  // Deallocate memory
  delete [] DeltaCellCenters; DeltaCellCenters = NULL;

}

#endif
