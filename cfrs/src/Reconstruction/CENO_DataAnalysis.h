/* DataAnalysis.h: Header file defining subroutines used for data analysis */

#ifndef _CENO_DATAANALYSIS_INCLUDED
#define _CENO_DATAANALYSIS_INCLUDED

#include <vector>
#include "../../../src_2D/Utilities/Utilities.h"
#include "include/TypeDefinition.h"
#include "../../../src_2D/Utilities/EpsilonTol.h"


/**********************************************************************************************************************
 ----------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------------
 *                                                                                                                    *
 *                                         1D ANALYSIS                                                                *
 *                                                                                                                    *
 ----------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------------
**********************************************************************************************************************/
template< class SolutionContainer>
void ComputeSmoothnessIndicator(SolutionContainer & SolnBlk, IndexType & i_index,
				const int iCell){

  static const int NumberOfParameters(SolutionContainer::NumberOfParameters);
  int StencilSize(i_index.size());
  double SS_Regression, SS_Residual, MeanSolution, alpha;
  double Temp, Coordinate, DeltaTol;
  int DOF(SolnBlk(iCell).NumberOfTaylorDerivatives()); 			/* degrees of freedom */
  int parameter, cell, ComputeSI;

  for (parameter=1; parameter<=NumberOfParameters; ++parameter){

    /* Mean Solution  */
    MeanSolution = SolnBlk(iCell).CellSolution(parameter);

    /* DeltaTol */
    DeltaTol = CENO_EpsilonTol::SquareToleranceAroundValue(MeanSolution);

    /* Compute the regression and residual sums */
    SS_Regression = 0.0;
    SS_Residual = 0.0;
    ComputeSI = OFF;		/* assume that the smoothness indicator is not computed but assigned */

    /* compute S(quare)S(sum)_Regression and SS_Residual  */
    for(cell=0; cell<StencilSize; ++cell){
      /* compute SS_Regression */
      Temp = SolnBlk(i_index[cell]).CellDeriv(0,parameter) - MeanSolution;
      Temp *= Temp;
      SS_Regression += Temp;

      /* Check if any of the Temp(s) is greater than DeltaTol */
      if ((ComputeSI==OFF) && (Temp > DeltaTol)){
	ComputeSI = ON; 	/* Decide to compute the smoothness indicator */
      }

      /* compute SS_Residual */
      Coordinate = SolnBlk(i_index[cell]).CellCenter();
      Temp = ( SolnBlk(i_index[cell]).CellDeriv(0,parameter) - 
	       SolnBlk(iCell).SolutionAtCoordinates(Coordinate) );
      SS_Residual += Temp*Temp;
    }

    /* Decide if the smoothness indicator is computed or not */
    if (ComputeSI){ 
      /* Compute alpha based on the ratio of the two sums */
      alpha = 1.0 - SS_Residual/SS_Regression;
    } else {
      /* Assign the perfect fit value to the smoothness indicator */
      alpha = 1.0;
    }

    // Compute the smoothness indicator
    SolnBlk(iCell).CellMCC(parameter) = (alpha*(StencilSize - DOF))/(max(CENO_EpsilonTol::epsilon,1.0 - alpha)*(DOF - 1));
  }
}


/**********************************************************************************************************************
 ----------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------------
 *                                                                                                                    *
 *                                         2D ANALYSIS                                                                *
 *                                                                                                                    *
 ----------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------------
**********************************************************************************************************************/
template< class SolutionContainer>
void ComputeSmoothnessIndicator(SolutionContainer & SolnBlk, const int *i_index, const int *j_index,
				const int & StencilSize,
				const int iCell, const int jCell){

  typedef typename SolutionContainer::CompCellType  ComputationalCellType;
  static const int NumberOfParameters = ComputationalCellType::NumberOfVariables;

  double SS_Regression, SS_Residual, MeanSolution, alpha;
  double Temp, DeltaTol;
  Vector2D Coordinate;
  int DOF = SolnBlk(iCell,jCell).NumberOfTaylorDerivatives(); 			/* degrees of freedom */
  int parameter, cell, ComputeSI;

  for (parameter=1; parameter<=NumberOfParameters; ++parameter){

    // Mean Solution
    MeanSolution = SolnBlk(iCell,jCell).CellSolution(parameter);

    /* DeltaTol */
    DeltaTol = CENO_EpsilonTol::SquareToleranceAroundValue(MeanSolution);

    /* Compute the regression and residual sums */
    SS_Regression = 0.0;
    SS_Residual = 0.0;
    ComputeSI = OFF;		/* assume that the smoothness indicator is not computed but assigned */

    /* compute S(quare)S(sum)_Regression and SS_Residual  */
    for(cell=0; cell<StencilSize; ++cell){
      /* compute SS_Regression */
      Temp = SolnBlk(i_index[cell],j_index[cell]).CellDeriv(0,0,parameter) - MeanSolution;
      Temp *= Temp;
      SS_Regression += Temp;
	
      /* Check if any of the Temp(s) is greater than DeltaTol */
      if ((ComputeSI==OFF) && (Temp > DeltaTol)){
	ComputeSI = ON; 	/* Decide to compute the smoothness indicator */
      }

      /* compute SS_Residual */
      Coordinate = SolnBlk(i_index[cell],j_index[cell]).CellCenter();
      Temp = (SolnBlk(i_index[cell],j_index[cell]).CellDeriv(0,0,parameter) - 
	      SolnBlk(iCell,jCell).SolutionAtCoordinates(Coordinate.x, Coordinate.y, parameter));
      SS_Residual += Temp*Temp;
    }

    /* Decide if the smoothness indicator is computed or not */
    if (ComputeSI){ 
      /* Compute alpha based on the ratio of the two sums */
      alpha = 1.0 - SS_Residual/SS_Regression;
    } else {
      /* Assign the perfect fit value to the smoothness indicator */
      alpha = 1.0;
    }
 
    // Compute the smoothness indicator
    SolnBlk(iCell,jCell).CellMCC(parameter) = (alpha*(StencilSize - DOF))/(max(CENO_EpsilonTol::epsilon,1.0 - alpha)*(DOF - 1.0));
  }
}

#endif
