/*!\file ReconstructionSolvers1D.h
  \brief Header file defining template functions which are used 
  for performing piecewise solution reconstruction based 
  on cell average solutions in one-dimension. */

#ifndef _RECONSTRUCTIONSOLVERS1D_INCLUDED
#define _RECONSTRUCTIONSOLVERS1D_INCLUDED


/* Include required C++ libraries. */
// NONE

/* Using std namespace functions */
// None

/* Include CFFC header files. */
#include "../Math/LinearSystems.h"
#include "ReconstructionHelpers.h"


/**********************************************************************************************************************
 ----------------------------------------------------------------------------------------------------------------------
 *                                  1D Memory Pools Implementation                                                    *
 ----------------------------------------------------------------------------------------------------------------------
**********************************************************************************************************************/
class MemoryStorageENO_1D {

 public:

  /* Stencil detection variable */
  int StencilStartIndex;
  int StencilPoints;
  int DivDiffOrder;

  /* Linear system variable */
  DenseMatrix A;
  ColumnVector U;
  ColumnVector X;

  /* Helper variable for setting the linear system */
  int NumOfUnknowns;
  double DeltaCellCenter;
  double CombP1X;
  double GeomIntegralOverCell;

  /* Constructors */
  MemoryStorageENO_1D(void): StencilStartIndex(0), StencilPoints(1), DivDiffOrder(1),
			     NumOfUnknowns(0), DeltaCellCenter(0.0), CombP1X(0.0), GeomIntegralOverCell(0.0){
    A.newsize(NumOfUnknowns, NumOfUnknowns);
    U.newsize(NumOfUnknowns);
    X.newsize(NumOfUnknowns);
  }

  MemoryStorageENO_1D(const int NumberOfUnknowns): StencilStartIndex(0), StencilPoints(1), DivDiffOrder(1),
						   NumOfUnknowns(NumberOfUnknowns),DeltaCellCenter(0.0),
						   CombP1X(0.0), GeomIntegralOverCell(0.0){
    A.newsize(NumOfUnknowns, NumOfUnknowns);
    U.newsize(NumOfUnknowns);
    X.newsize(NumOfUnknowns);
  }
   
  void newsize(const int NumberOfUnknowns){
    NumOfUnknowns = NumberOfUnknowns;
    A.newsize(NumOfUnknowns, NumOfUnknowns);
    U.newsize(NumOfUnknowns);
    X.newsize(NumOfUnknowns);
  }

  // Use the default destructor
  // ~MemoryStorageENO_1D(void){};
  
  friend std::ostream & operator<< (std::ostream & os, const MemoryStorageENO_1D & Obj){
    
    os << "\t Memory Storage for ENO 1D Reconstruction\n";
    
    os << "\tStencil Detection Variables:\n"
	   << "StencilStartIndex=" << Obj.StencilStartIndex << std::endl
       << "StencilPoints=" << Obj.StencilPoints << std::endl
       << "DividedDifferenceOrder=" << Obj.DivDiffOrder << std::endl;
    
    os << "\t Linear System Variables:\n"
       << "A=" << Obj.A << std::endl
       << "U=" << Obj.U << std::endl
       << "X=" << Obj.X << std::endl;
    
    os << "\t Helper Variables:\n"
       << "NumOfUnknowns=" << Obj.NumOfUnknowns << std::endl
       << "DeltaCellCenter=" << Obj.DeltaCellCenter << std::endl
       << "Combination=" << Obj.CombP1X << std::endl
       << "GeomIntegralOverCell=" << Obj.GeomIntegralOverCell << std::endl;
    
    return os;
  }
  
};

/**********************************************************************************************************************
 ----------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------------
 *                                                                                                                    *
 *                                         1D Reconstructions                                                         *
 *                                                                                                                    *
 ----------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------------
**********************************************************************************************************************/

/****************************************************************************
 * TEMPLATE Function: kExact_Reconstruction for 1D                          *
 *                                                                          *
 * This subroutine determines the coefficients of a Taylor series expansion *
 * which approximates the solution over the domain of the cell specified    *
 * by "iCell" index.                                                        *
 *                                                                          *
 ***************************************************************************/
#ifdef _CENO_SPEED_EFFICIENT
template< class SolutionContainer> inline
void kExact_Reconstruction (SolutionContainer * SolnBlk, IndexType & i_index, const int iCell)
{

  static const int NumberOfParameters = SolutionContainer::NumberOfVariables;

  // SET THE VARIABLES USED IN THE RECONSTRUCTION PROCESS
  int ND(SolnBlk[iCell].NumberOfTaylorDerivatives());	      /* number of Taylor expansion coefficients */
  int StencilSize(i_index.size());                            /* number of cells used in the reconstruction stencil */
  ColumnVector Delta_U(StencilSize-1);                        /* vector for storing U[neighbour]-U[cell] for each parameter*/
  ColumnVector X(StencilSize-1);                              /* solution vector */
  int cell, i, parameter;

  // *********  Assign the average solution to D0  ***********
  SolnBlk[iCell].CellDeriv(0,true,true).D( ) = SolnBlk[iCell].CellSolutionPrimVar( );

  if (ND == 1){ // piecewise constant
    return;
  }

  // START: Compute for every parameter the high-order approximation
  // ***************************************************************
  for (parameter = 1; parameter <= NumberOfParameters; ++parameter){
    // Step 1. SET the vector Delta_U of the linear system (RHS) for the current parameter ***
    for (cell=1 ; cell<StencilSize; ++cell){ //for each cell in the stencil
      Delta_U(cell-1) = ( SolnBlk[i_index[cell]].CellSolutionPrimVar(parameter) -
			  SolnBlk[iCell].CellSolutionPrimVar(parameter) );
      Delta_U(cell-1) *= SolnBlk[iCell].GeomWeights(cell);
    }

    // Step 2. Find the solution of the linear-system for the current parameter
    X = SolnBlk[iCell].LHS() * Delta_U;

    // Step 3. Update the coefficients D (derivatives) for the current parameter
    for (i=1; i<=SolnBlk[iCell].CellDeriv().LastElem(); ++i){
      SolnBlk[iCell].CellDeriv(i,true,true).D(parameter) = X(i-1);

      /* this equation makes sure that the mean conservation of each parameter is satisfied inside the reconstructed cell */
      SolnBlk[iCell].CellDeriv(0,true,true).D(parameter) -= SolnBlk[iCell].CellGeomCoeff(i,true) * X(i-1);
    }
  }
}
//********************************************* ELSE ******************************************
#else
//********************************************* ELSE ******************************************
template< class SolutionContainer> inline
void kExact_Reconstruction (SolutionContainer * SolnBlk, IndexType & i_index, const int iCell)
{

  static const int NumberOfParameters = SolutionContainer::NumberOfVariables;

  // SET THE VARIABLES USED IN THE RECONSTRUCTION PROCESS
  int ND(SolnBlk[iCell].NumberOfTaylorDerivatives());	      /* number of Taylor expansion coefficients */
  int StencilSize(i_index.size());                            /* number of cells used in the reconstruction stencil */
  DenseMatrix  A(StencilSize-1,ND-1);                         /* the matrix which the linear system is solved for */
  DenseMatrix  All_Delta_U(StencilSize-1,NumberOfParameters); /* matrix for storing U[neighbour]-U[cell] */
  ColumnVector DistanceCellCenters(StencilSize);              /* stores the distance between the cell center of neighbour
  						                 cells and the I cell */
  ColumnVector GeomWeights(StencilSize);                      /* The column vector of the geometric weights */
  int krank;                                                  /* the final rank of A is returned here */

  int cell, i, IndexSumX, parameter, P1;
  double CombP1X, MaxWeight(0.0), PowDistanceXC;

  // *********  Assign the average solution to D0  ***********
  SolnBlk[iCell].CellDeriv(0,true,true).D( ) = SolnBlk[iCell].CellSolutionPrimVar( );

  if (ND == 1){ // piecewise constant
    return;
  }

  // START:   Set the LHS and RHS of the linear system 
  // ***************************************************

  // Step1. Compute the normalized geometric weights
  for (cell=1; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil

    /* Compute the distance between the cell centers of the neighbour and the reconstructed cell */
    DistanceCellCenters(cell) = SolnBlk[i_index[cell]].CellCenter() - SolnBlk[iCell].CellCenter();

    /* Compute the geometric weights and their maximum value (this is used for normalization)
       based on the distance to each control volume */
    GeomWeights(cell) = CENO_Geometric_Weighting(DistanceCellCenters(cell));
    MaxWeight = max(MaxWeight,GeomWeights(cell));
  }

  // Step2. Set the approximate equations
  for (cell=1 ; cell<StencilSize; ++cell){ //for each cell in the stencil
    
    // compute the normalized geometric weight
    GeomWeights(cell) /= MaxWeight;

    // *** SET the matrix A of the linear system (LHS) ***
    /* compute for each derivative the corresponding entry in the matrix of the linear system */
    for (i=1; i<=SolnBlk[iCell].CellDeriv().LastElem(); ++i){
      // build the row of the matrix
      P1 = SolnBlk[iCell].CellDeriv(i,true,true).P1(); // identify P1
      A(cell-1,i-1) = 0.0;  // set sumation variable to zero
      CombP1X = 1.0;        // the binomial coefficient "nC k" for k=0 is 1
      PowDistanceXC = 1.0;  // initialize PowDistanceXC

      // ********* Compute geometric integral over the neighbour's domain **********
      for (IndexSumX = 0; IndexSumX<=P1; ++IndexSumX){
	A(cell-1,i-1) += ( CombP1X * PowDistanceXC *
			   SolnBlk[i_index[cell]].CellGeomCoeff()(P1-IndexSumX) );
	
	// update the binomial coefficients
	CombP1X = (P1-IndexSumX)*CombP1X/(IndexSumX+1);  // the index is still the old one => expression for "nC k+1"
	PowDistanceXC *= DistanceCellCenters(cell);      // Update PowDistanceXC
      } //endfor


      // subtract the corresponding geometric moment of (iCell) cell
      A(cell-1,i-1) -= SolnBlk[iCell].CellGeomCoeff(i,true);

      // apply geometric weighting
      A(cell-1,i-1) *= GeomWeights(cell);
    }
    
    // *** SET the matrix All_Delta_U of the linear system (RHS) ***
    for (parameter = 1; parameter <= NumberOfParameters; ++parameter){
      All_Delta_U(cell-1,parameter-1) = ( SolnBlk[i_index[cell]].CellSolutionPrimVar(parameter) -
					  SolnBlk[iCell].CellSolutionPrimVar(parameter) );
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
  for (i=1; i<=SolnBlk[iCell].CellDeriv().LastElem(); ++i){
    for (parameter = 1; parameter <= NumberOfParameters; ++parameter){
      SolnBlk[iCell].CellDeriv(i,true,true).D(parameter) = All_Delta_U(i-1,parameter-1);
      
      /* this equation makes sure that the mean conservation of each parameter is satisfied inside the reconstructed cell */
      SolnBlk[iCell].CellDeriv(0,true,true).D(parameter) -= SolnBlk[iCell].CellGeomCoeff(i,true) * All_Delta_U(i-1,parameter-1);
    }
  }
}
#endif //_CENO_SPEED_EFFICIENT


/****************************************************************************
 * TEMPLATE Function: kExact_Reconstruction_LHS for 1D                      *
 *                                                                          *
 * This subroutine forms the left hand side (LHS) term that arise in        *
 * solving the CENO solution reconstruction of the cell specified           *
 * by "iCell" index.                                                        *
 * The exactly satisfied conservation equation for cell iCell is left out.  *
 *                                                                          *
 ***************************************************************************/
#ifdef _CENO_SPEED_EFFICIENT
template< class SolutionContainer> inline
void kExact_Reconstruction_LHS (SolutionContainer * SolnBlk, IndexType & i_index, const int iCell)
{

  // SET THE VARIABLES USED IN THE RECONSTRUCTION PROCESS
  int ND(SolnBlk[iCell].NumberOfTaylorDerivatives());	      /* number of Taylor expansion coefficients */
  int StencilSize(i_index.size());                            /* number of cells used in the reconstruction stencil */
  ColumnVector DistanceCellCenters(StencilSize);              /* stores the distance between the cell center of neighbour
  						                 cells and the I cell */
  int cell, i, IndexSumX, P1;
  double CombP1X, MaxWeight(0.0), PowDistanceXC;

  if (ND == 1){ // piecewise constant
    return;
  }

  // START:   Set the LHS of the linear system 
  // ***************************************************

  // Step1. Compute the normalized geometric weights
  for (cell=1; cell<StencilSize; ++cell){ //for each neighbour cell in the stencil

    /* Compute the distance between the cell centers of the neighbour and the reconstructed cell */
    DistanceCellCenters(cell) = SolnBlk[i_index[cell]].CellCenter() - SolnBlk[iCell].CellCenter();

    /* Compute the geometric weights and their maximum value (this is used for normalization)
       based on the distance to each control volume */
    SolnBlk[iCell].GeomWeights(cell) = CENO_Geometric_Weighting(DistanceCellCenters(cell));
    MaxWeight = max(MaxWeight,SolnBlk[iCell].GeomWeights(cell));
  }

  // Step2. Set the approximate equations
  for (cell=1 ; cell<StencilSize; ++cell){ //for each cell in the stencil
    
    // compute the normalized geometric weight
    SolnBlk[iCell].GeomWeights(cell) /= MaxWeight;

    // *** SET the matrix A of the linear system (LHS) ***
    /* compute for each derivative the corresponding entry in the matrix of the linear system */
    for (i=1; i<=SolnBlk[iCell].CellDeriv().LastElem(); ++i){
      // build the row of the matrix
      P1 = SolnBlk[iCell].CellDeriv(i,true,true).P1(); // identify P1
      SolnBlk[iCell].LHS(cell-1,i-1) = 0.0;  // set sumation variable to zero
      CombP1X = 1.0;        // the binomial coefficient "nC k" for k=0 is 1
      PowDistanceXC = 1.0;  // initialize PowDistanceXC

      // ********* Compute geometric integral over the neighbour's domain **********
      for (IndexSumX = 0; IndexSumX<=P1; ++IndexSumX){
	SolnBlk[iCell].LHS(cell-1,i-1) += ( CombP1X * PowDistanceXC *
					    SolnBlk[i_index[cell]].CellGeomCoeff()(P1-IndexSumX) );
	
	// update the binomial coefficients
	CombP1X = (P1-IndexSumX)*CombP1X/(IndexSumX+1);  // the index is still the old one => expression for "nC k+1"
	PowDistanceXC *= DistanceCellCenters(cell);      // Update PowDistanceXC
      } //endfor


      // subtract the corresponding geometric moment of (iCell) cell
      SolnBlk[iCell].LHS(cell-1,i-1) -= SolnBlk[iCell].CellGeomCoeff(i,true);

      // apply geometric weighting
      SolnBlk[iCell].LHS(cell-1,i-1) *= SolnBlk[iCell].GeomWeights(cell);
    }
  }

  // Matrix A of the linear system (LHS) built. The same matrix is used for all the variables (same geometry) and
  // at every time step as long as the mesh is the same.
  // **********************************************************************
}
#endif //_CENO_SPEED_EFFICIENT

/****************************************************************************
 * TEMPLATE Function: ENO_Reconstruction                                    *
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
void ENO_Reconstruction (SolutionContainer * SolnBlk, const int iCell, MemoryStorageENO_1D & Mem)
{

  int cell, Counter, i, P1;
  double PowDistanceXC;
  int IndexSumX;
  Mem.StencilStartIndex = iCell;
  Mem.StencilPoints = 1;
  Mem.DivDiffOrder=1;

  for(int parameter=1; parameter <= SolutionContainer::NumberOfVariables; ++parameter){

    /* Determine the "smoothest" stencil */
    while(Mem.StencilPoints < Mem.NumOfUnknowns){

      // determine the starting position of the stencil    
      if (fabs(DivDifference(SolnBlk,parameter,Mem.DivDiffOrder,Mem.StencilStartIndex-1)) <
	  fabs(DivDifference(SolnBlk,parameter,Mem.DivDiffOrder,Mem.StencilStartIndex))){

	// change starting position of the stencil
	--Mem.StencilStartIndex;
      }

      /* advance to the next level of finite differences */
      ++Mem.DivDiffOrder;

      /* add one more point to the stencil */
      ++Mem.StencilPoints;
    }

    /* Solve the reconstruction using the obtained stencil */

    /* Set the LHS and RHS */
    for (cell=Mem.StencilStartIndex, Counter=0; Counter<Mem.NumOfUnknowns; ++cell, ++Counter){
      //for each cell in the stencil

      /*** SET the matrix A of the linear system (LHS) ***/

      /* Compute the distance between the cell center of the current cell and the reconstructed cell */
      Mem.DeltaCellCenter = SolnBlk[cell].CellCenter() - SolnBlk[iCell].CellCenter();

      /* compute for each derivative the corresponding entry in the matrix of the linear system */
      for( i = SolnBlk[iCell].CellDeriv().FirstElem(); i<=SolnBlk[iCell].CellDeriv().LastElem(); ++i){

	// build the row of the matrix
	P1 = SolnBlk[iCell].CellDeriv(i,true,true).P1(); // identify P1
	Mem.GeomIntegralOverCell = 0.0;  // set sumation variable to zero
	Mem.CombP1X = 1.0;               // the binomial coefficient "nC k" for k=0 is 1
	PowDistanceXC = 1.0;             // initialize PowDistanceXC

	// ********* Compute geometric integral over the neighbour's domain **********
	for (IndexSumX = 0; IndexSumX<=P1; ++IndexSumX){
	  Mem.GeomIntegralOverCell += (Mem.CombP1X * PowDistanceXC *
				       SolnBlk[cell].CellGeomCoeff()(P1-IndexSumX) );
	  // update the binomial coefficients
	  Mem.CombP1X = (P1-IndexSumX)*Mem.CombP1X/(IndexSumX+1); // the index is still the old one => expression for "nC k+1"
	  PowDistanceXC *= Mem.DeltaCellCenter;      // Update PowDistanceXC
	}

	Mem.A(Counter,i) = Mem.GeomIntegralOverCell;
      }
      
      // *** SET U of the linear system (RHS) ***
      Mem.U(Counter) = SolnBlk[cell].CellSolutionPrimVar(parameter);

    }

    /* Solve the linear system */
    Solve_LU_Decomposition(Mem.A,Mem.U,Mem.X);

    /* Update the coefficients D (derivatives) */
    for( i = SolnBlk[iCell].CellDeriv().FirstElem(); i<=SolnBlk[iCell].CellDeriv().LastElem(); ++i){
      SolnBlk[iCell].CellDeriv(i,parameter) = Mem.X(i);
    }

    /* Reset stencil variables for the reconstruction of the next parameter */
    Mem.StencilStartIndex = iCell;
    Mem.StencilPoints = 1;
    Mem.DivDiffOrder = 1;
  }

}

/****************************************************************************
 * TEMPLATE Function: ENO_Characteristics_Reconstruction                    *
 *                                                                          *
 * This function determines the coefficients of a Taylor series expansion   *
 * which approximates the solution over the domain of the cell specified    *
 * by "i" index.                                                            *
 * The subroutine reconstruct the solution in characteristic variables      *
 * which are computed from either the primitive or conserved variables.     *
 * The method used is presented in:                                         *
 * "Uniformly High Order Accurate Essentially Non-oscillatory Schemes, III" *
 * Ami Harten, Journal of Computational Physics, 71, (1987)                 *
 *                                                                          *
 ***************************************************************************/
template< class SolutionContainer> inline
void ENO_Characteristics_Reconstruction (SolutionContainer * SolnBlk, const int iCell, MemoryStorageENO_1D & Mem)
{

  int cell, Counter, i;
  int IndexSumX;
  Mem.StencilStartIndex = iCell;
  Mem.StencilPoints = 1;
  Mem.DivDiffOrder=1;

  /* Compute the characteristic variables of the stencil candidates using the local left-eigenvector system */

  /* 1st choice: use primitive variables */
  #ifdef UsePrimitiveVar
  for( i = iCell-Mem.NumOfUnknowns+1; i<iCell+Mem.NumOfUnknowns; ++i){
    SolnBlk[i].CellSolutionCharactVar() = SolnBlk[i].CellSolutionPrimVar().PrimitiveVarToCharactVar(SolnBlk[iCell].CellSolutionPrimVar());
  }
  #else
  /* 2nd choice: use conserved variables */
  for( i = iCell-Mem.NumOfUnknowns+1; i<iCell+Mem.NumOfUnknowns; ++i){
    SolnBlk[i].CellSolutionCharactVar() = SolnBlk[i].CellSolutionConsVar().ConservedVarToCharactVar(SolnBlk[iCell].CellSolutionConsVar());
  }
  #endif

  for(int parameter=1; parameter <= SolutionContainer::NumberOfVariables; ++parameter){

    /* Determine the "smoothest" stencil */
    while(Mem.StencilPoints < Mem.NumOfUnknowns){

      // determine the starting position of the stencil    
      if (fabs(DivDifferenceCharacteristic(SolnBlk,parameter,Mem.DivDiffOrder,Mem.StencilStartIndex-1)) <
	  fabs(DivDifferenceCharacteristic(SolnBlk,parameter,Mem.DivDiffOrder,Mem.StencilStartIndex))){

	// change starting position of the stencil
	--Mem.StencilStartIndex;
      }

      /* advance to the next level of finite differences */
      ++Mem.DivDiffOrder;

      /* add one more point to the stencil */
      ++Mem.StencilPoints;
    }

    /* Solve the reconstruction using the obtained stencil */
    
    /* Set the LHS and RHS */
    for (cell=Mem.StencilStartIndex, Counter=0; Counter<Mem.NumOfUnknowns; ++cell, ++Counter){
      //for each cell in the stencil

      /*** SET the matrix A of the linear system (LHS) ***/

      /* Compute the distance between the cell center of the current cell and the reconstructed cell */
      Mem.DeltaCellCenter = SolnBlk[cell].CellCenter() - SolnBlk[iCell].CellCenter();

   
      /* compute for each derivative the corresponding entry in the matrix of the linear system */
      for( i = SolnBlk[iCell].CellDeriv().FirstElem(); i<=SolnBlk[iCell].CellDeriv().LastElem(); ++i){

	// build the row of the matrix
	Mem.GeomIntegralOverCell = 0.0;  // set sumation variable to zero
	// Compute geometric integral over the neighbour's domain
	Mem.CombP1X = 1.0; // the binomial coefficient "nC k" for k=0 is 1
	for (IndexSumX = 0; IndexSumX<=(int)SolnBlk[iCell].CellDeriv(i,true,true).P1(); ++IndexSumX){
	  Mem.GeomIntegralOverCell += (Mem.CombP1X * pow(Mem.DeltaCellCenter,IndexSumX)*
				       SolnBlk[cell].CellGeomCoeff()(SolnBlk[iCell].CellDeriv(i,true,true).P1()-IndexSumX));
	  // update the binomial coefficients
	  Mem.CombP1X = (SolnBlk[iCell].CellDeriv(i,true,true).P1()-IndexSumX)*Mem.CombP1X/(IndexSumX+1);
	  // the index is still the old one => expression for "nC k+1"
	}

	Mem.A(Counter,i) = Mem.GeomIntegralOverCell;
      }
      
      // *** SET U of the linear system (RHS) ***
      Mem.U(Counter) = SolnBlk[cell].CellSolutionCharactVar(parameter);

    }

    /* Solve the linear system */
    Solve_LU_Decomposition(Mem.A,Mem.U,Mem.X);

    /* Update the coefficients D (derivatives) */
    for( i = SolnBlk[iCell].CellDeriv().FirstElem(); i<=SolnBlk[iCell].CellDeriv().LastElem(); ++i){
      SolnBlk[iCell].CellDeriv(i,parameter) = Mem.X(i);
    }

    /* Reset stencil variables for the reconstruction of the next parameter */
    Mem.StencilStartIndex = iCell;
    Mem.StencilPoints = 1;
    Mem.DivDiffOrder = 1;
  }
}

/**********************************************************************************************************************
 ----------------------------------------------------------------------------------------------------------------------
 *                                  1D Devided Differences                                                            *
 ----------------------------------------------------------------------------------------------------------------------
**********************************************************************************************************************/
template< class SolutionContainer> inline
double DivDifference (SolutionContainer * SolnBlk, int Var, int Order, int SI)
{

  // SI - StartIndex

  switch(Order){

  case 1:
    return (SolnBlk[SI].CellSolutionPrimVar(Var) - SolnBlk[SI+1].CellSolutionPrimVar(Var))/
      (SolnBlk[SI].CellCenter() - SolnBlk[SI+1].CellCenter());
    break;

  case 2: 
    return ((SolnBlk[SI].CellSolutionPrimVar(Var) - SolnBlk[SI+1].CellSolutionPrimVar(Var))*
            (SolnBlk[SI+1].CellCenter() - SolnBlk[SI+2].CellCenter()) - 
	    (SolnBlk[SI+1].CellSolutionPrimVar(Var) - SolnBlk[SI+2].CellSolutionPrimVar(Var))*
            (SolnBlk[SI].CellCenter() - SolnBlk[SI+1].CellCenter()))/

      ((SolnBlk[SI].CellCenter() - SolnBlk[SI+1].CellCenter())*
       (SolnBlk[SI+1].CellCenter() - SolnBlk[SI+2].CellCenter())*
       (SolnBlk[SI].CellCenter() - SolnBlk[SI+2].CellCenter()));
    break;

  case 3:

    return
      ( (SolnBlk[SI].CellCenter() - SolnBlk[SI+1].CellCenter())*
	(SolnBlk[SI].CellCenter() - SolnBlk[SI+2].CellCenter())*
	( (SolnBlk[SI+1].CellCenter() - SolnBlk[SI+2].CellCenter())*
	  (SolnBlk[SI+2].CellSolutionPrimVar(Var) - SolnBlk[SI+3].CellSolutionPrimVar(Var)) -
	  (SolnBlk[SI+2].CellCenter() - SolnBlk[SI+3].CellCenter())*
	  (SolnBlk[SI+1].CellSolutionPrimVar(Var) - SolnBlk[SI+2].CellSolutionPrimVar(Var)) ) +

	(SolnBlk[SI+1].CellCenter() - SolnBlk[SI+3].CellCenter())*
	(SolnBlk[SI+2].CellCenter() - SolnBlk[SI+3].CellCenter())*
	( (SolnBlk[SI+1].CellCenter() - SolnBlk[SI+2].CellCenter())*
	  (SolnBlk[SI].CellSolutionPrimVar(Var) - SolnBlk[SI+1].CellSolutionPrimVar(Var)) -
	  (SolnBlk[SI].CellCenter() - SolnBlk[SI+1].CellCenter())*
	  (SolnBlk[SI+1].CellSolutionPrimVar(Var) - SolnBlk[SI+2].CellSolutionPrimVar(Var)) )
	) /

      ((SolnBlk[SI].CellCenter() - SolnBlk[SI+1].CellCenter())*
       (SolnBlk[SI].CellCenter() - SolnBlk[SI+2].CellCenter())*
       (SolnBlk[SI].CellCenter() - SolnBlk[SI+3].CellCenter())*
       (SolnBlk[SI+1].CellCenter() - SolnBlk[SI+2].CellCenter())*
       (SolnBlk[SI+1].CellCenter() - SolnBlk[SI+3].CellCenter())*
       (SolnBlk[SI+2].CellCenter() - SolnBlk[SI+3].CellCenter())
       );

  default:
    return (DivDifference(SolnBlk,Var,Order-1,SI) - DivDifference(SolnBlk,Var,Order-1,SI+1))/ \
      (SolnBlk[SI].CellCenter() - SolnBlk[SI+Order].CellCenter());

  }
}

template< class SolutionContainer> inline
double DivDifferenceCharacteristic (SolutionContainer * SolnBlk, int Var, int Order, int SI)
{

  // SI - StartIndex

  switch(Order){

  case 1:
    return (SolnBlk[SI].CellSolutionCharactVar(Var) - SolnBlk[SI+1].CellSolutionCharactVar(Var))/
      (SolnBlk[SI].CellCenter() - SolnBlk[SI+1].CellCenter());
    break;

  case 2: 
    return ((SolnBlk[SI].CellSolutionCharactVar(Var) - SolnBlk[SI+1].CellSolutionCharactVar(Var))*
            (SolnBlk[SI+1].CellCenter() - SolnBlk[SI+2].CellCenter()) - 
	    (SolnBlk[SI+1].CellSolutionCharactVar(Var) - SolnBlk[SI+2].CellSolutionCharactVar(Var))*
            (SolnBlk[SI].CellCenter() - SolnBlk[SI+1].CellCenter()))/

      ((SolnBlk[SI].CellCenter() - SolnBlk[SI+1].CellCenter())*
       (SolnBlk[SI+1].CellCenter() - SolnBlk[SI+2].CellCenter())*
       (SolnBlk[SI].CellCenter() - SolnBlk[SI+2].CellCenter()));
    break;

  case 3:

    return 
      ( (SolnBlk[SI].CellCenter() - SolnBlk[SI+1].CellCenter())*
	(SolnBlk[SI].CellCenter() - SolnBlk[SI+2].CellCenter())*
	( (SolnBlk[SI+1].CellCenter() - SolnBlk[SI+2].CellCenter())*
	  (SolnBlk[SI+2].CellSolutionCharactVar(Var) - SolnBlk[SI+3].CellSolutionCharactVar(Var)) -
	  (SolnBlk[SI+2].CellCenter() - SolnBlk[SI+3].CellCenter())*
	  (SolnBlk[SI+1].CellSolutionCharactVar(Var) - SolnBlk[SI+2].CellSolutionCharactVar(Var)) ) +

	(SolnBlk[SI+1].CellCenter() - SolnBlk[SI+3].CellCenter())*
	(SolnBlk[SI+2].CellCenter() - SolnBlk[SI+3].CellCenter())*
	( (SolnBlk[SI+1].CellCenter() - SolnBlk[SI+2].CellCenter())*
	  (SolnBlk[SI].CellSolutionCharactVar(Var) - SolnBlk[SI+1].CellSolutionCharactVar(Var)) -
	  (SolnBlk[SI].CellCenter() - SolnBlk[SI+1].CellCenter())*
	  (SolnBlk[SI+1].CellSolutionCharactVar(Var) - SolnBlk[SI+2].CellSolutionCharactVar(Var)) )
	) /

      ((SolnBlk[SI].CellCenter() - SolnBlk[SI+1].CellCenter())*
       (SolnBlk[SI].CellCenter() - SolnBlk[SI+2].CellCenter())*
       (SolnBlk[SI].CellCenter() - SolnBlk[SI+3].CellCenter())*
       (SolnBlk[SI+1].CellCenter() - SolnBlk[SI+2].CellCenter())*
       (SolnBlk[SI+1].CellCenter() - SolnBlk[SI+3].CellCenter())*
       (SolnBlk[SI+2].CellCenter() - SolnBlk[SI+3].CellCenter())
       );

  default:
    return (DivDifferenceCharacteristic(SolnBlk,Var,Order-1,SI) - DivDifferenceCharacteristic(SolnBlk,Var,Order-1,SI+1))/ \
      (SolnBlk[SI].CellCenter() - SolnBlk[SI+Order].CellCenter());

  }
}

#endif
