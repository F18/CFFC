#ifndef _RECONSTRUCTIONHELPERS_CFRS_INCLUDED
#define _RECONSTRUCTIONHELPERS_CFRS_INCLUDED

#include <vector>
#include "../../../src_2D/Math/LinearSystems.h"
#include "../../../src_2D/Math/Vector2D.h"
#include "../../../src_2D/Utilities/Utilities.h"
#include "../../../src_2D/HighOrderReconstruction/ReconstructionHelpers.h"
#include "include/TypeDefinition.h"

//#include "../../../src_3D/Math/LinearSystems.h"
//#include "../../../src_3D/Math/Vector3D.h"
//#include "../../../src_3D/Utilities/Utilities.h"

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
  
/*! Set the stencil for the 3D kExact reconstruction
******************************************************************************************/
void MakeReconstructionStencil(const int & rings, const int & iCell, const int & jCell, const int & kCell,
			       vector<int> & i_index, vector<int> & j_index, vector<int> & k_index);
void MakeReconstructionStencil(const int & rings, const int & iCell, const int & jCell, const int & kCell,
			       int *i_index, int *j_index, int *k_index);

#endif
