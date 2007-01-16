/* DataAnalysis.h: Header file defining subroutines used for data analysis */

#ifndef _DATAANALYSIS_INCLUDED
#define _DATAANALYSIS_INCLUDED

#include <vector>
#include "include/require.h"
#include "include/TypeDefinition.h"

template< class SolutionContainer>
void MultipleCorrelationCoefficient(SolutionContainer & SolnBlk, IndexType & i_index,
				    const int iCell){

  int StencilSize(i_index.size());

  double SS_Regression(0.0), SS_Residual(0.0), MeanSolution;
  double Temp;
  double MultipleCorrCoeff, MaxDeltaU(0.0);
  double Coordinate;
  int DOF = SolnBlk(iCell).NumberOfTaylorDerivatives(); 			/* degrees of freedom */
  int cell;

  /* Define the worst fit */
  MeanSolution = SolnBlk(iCell).CellSolution();

  // Get the maximum DeltaU over the stencil
  for( cell=1; cell<StencilSize; ++cell){ // for cell=0 the index is (iCell)
    MaxDeltaU = max(MaxDeltaU, fabs(SolnBlk(i_index[cell]).CellSolution() - MeanSolution));
  }

  // Decide if the smoothness indicator is computed or not
  if (MaxDeltaU/(1.0 + fabs(MeanSolution)) > 1.0e-6){ 
    /* compute S(quare)S(sum)_Regression and SS_Residual  */
    for(cell=0; cell<StencilSize; ++cell){
      /* compute SS_Regression */
      Temp = SolnBlk(i_index[cell]).CellDeriv(0) - MeanSolution;
      SS_Regression += Temp*Temp;

      /* compute SS_Residual */
      Coordinate = SolnBlk(i_index[cell]).CellCenter();
      Temp = SolnBlk(i_index[cell]).CellDeriv(0) - SolnBlk(iCell).SolutionAtCoordinates(Coordinate);
      SS_Residual += Temp*Temp;
    }

    /* Multiple-Correlation Coefficient */
    MultipleCorrCoeff = 1.0 - SS_Residual/SS_Regression;
  } else {

    // Assign the perfect fit value to the smoothness indicator
    MultipleCorrCoeff = 1.0;
  }

  // Compute the smoothness indicator
  SolnBlk(iCell).CellMCC() = (MultipleCorrCoeff*(StencilSize - DOF))/(max(1.0e-6,1.0 - MultipleCorrCoeff)*(DOF - 1));
}

template< class SolutionContainer>
void MultipleCorrelationCoefficient(SolutionContainer & SolnBlk, const int *i_index, const int *j_index, const int & StencilSize,
				    const int iCell, const int jCell){

  typedef typename SolutionContainer::CompCellType  ComputationalCellType;
  static const int NumberOfParameters = ComputationalCellType::NumberOfVariables;

  double SS_Regression, SS_Residual, MeanSolution;
  double Temp;
  double MultipleCorrCoeff, F_R;
  Vector2D Coordinate;
  int DOF = SolnBlk(iCell,jCell).NumberOfTaylorDerivatives(); 			/* degrees of freedom */
  int parameter, cell;
  double MaxDeltaU(0.0);

  for (parameter=1; parameter<=NumberOfParameters; ++parameter){
    // Mean Solution
    MeanSolution = SolnBlk(iCell,jCell).CellSolution(parameter);
    // Reset MaxDeltaU
    MaxDeltaU = 0.0;

    // Get the maximum DeltaU over the stencil
    for( cell=1; cell<StencilSize; ++cell){ // for cell=0 the indexes are (iCell,jCell)
      MaxDeltaU = max(MaxDeltaU, fabs(SolnBlk(i_index[cell],j_index[cell]).CellSolution(parameter) - MeanSolution));
    }

    // Decide if the smoothness indicator is computed or not
    if (MaxDeltaU/(1.0 + fabs(MeanSolution)) > 1.0e-15){ 
      SS_Regression = 0.0;
      SS_Residual = 0.0;

      for(cell=0; cell<StencilSize; ++cell){
	Temp = SolnBlk(i_index[cell],j_index[cell]).CellDeriv(0,0,parameter) - MeanSolution;
	SS_Regression += Temp*Temp;
	
	Coordinate = SolnBlk(i_index[cell],j_index[cell]).CellCenter();
	Temp = (SolnBlk(i_index[cell],j_index[cell]).CellDeriv(0,0,parameter) - 
		SolnBlk(iCell,jCell).SolutionAtCoordinates(Coordinate.x, Coordinate.y, parameter));
	SS_Residual += Temp*Temp;
      }

      MultipleCorrCoeff = 1.0 - SS_Residual/SS_Regression;
    } else {
      // Assign the perfect fit value to the smoothness indicator
      MultipleCorrCoeff = 1.0;
    }
    
    F_R = (MultipleCorrCoeff*(StencilSize - DOF))/(max(1.0e-6,1.0 - MultipleCorrCoeff)*(DOF - 1.0));

    SolnBlk(iCell,jCell).CellMCC() = F_R;
  }
}

#endif
