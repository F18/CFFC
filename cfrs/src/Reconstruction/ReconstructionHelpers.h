#ifndef _RECONSTRUCTIONHELPERS_INCLUDED
#define _RECONSTRUCTIONHELPERS_INCLUDED

#include <vector>
#include "Math/LinearSystems.h"
#include "include/require.h"
#include "include/TypeDefinition.h"

/* Function Definition */
void AnalyzeWeights(DenseMatrix & A, ColumnVector & Weights, const double & MinDistance, 
		    const int & ReconstructionOrder, int & FinalOrder, const double & CutoffKnob);

int InverseOrder2D(const int & ReconstructionOrder, const int & ComputedDerivatives);

inline int NumberOfDerivatives2D(const int & ReconstructionOrder){
  return (ReconstructionOrder+1)*(ReconstructionOrder+2)/2;
}

void DetermineMaxDeltaSolutionStencil(ColumnVector & MaxDeltaSolutionStencil,
				      DenseMatrix & All_Delta_U,
				      const int & NumberOfVariables);

void DetermineMaxDeltaSolutionStencil(ColumnVector & MaxDeltaSolutionStencil,
				      ColumnVector & MinSolutionStencil,
				      const int & NumberOfVariables);

void Weighting_LS_Problem (DenseMatrix & LHS, DenseMatrix & RHS,
			   ColumnVector & DistanceCellCenters,
			   ColumnVector & W);

void Weighting_LS_Problem (DenseMatrix & A, ColumnVector & Weights,
			   const ColumnVector & DistanceCellCenters,
			   ColumnVector & Delta_U, const ColumnVector & Delta_U_Original,
			   const double & MaxDeltaSolutionStencil,
			   const double &MaxDeltaSolutionOverDomain,
			   const double & CharacteristicLength,
			   int & ReconstructionOrder, int & FinalOrder, const double & CutoffKnob);

void Carl_Weighting_LS_Problem (DenseMatrix & A, ColumnVector & Weights,
				ColumnVector & DistanceCellCenters,
				ColumnVector & Delta_U, ColumnVector & Delta_U_Original,
				double & MaxDeltaSolutionOverDomain, double & CharacteristicLength,
				int & ReconstructionOrder, int & FinalOrder, const double & CutoffKnob);

void DetermineDistanceCellCenters(ColumnVector & DistanceCellCenters, vector<Vector2D> & DeltaCellCenters);

int InverseOrder2D(const int & ReconstructionOrder, const int & ComputedDerivatives);
  


#endif
