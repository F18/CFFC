/* ReconstructionSolver.cc: Source file defining the Reconstruction functions */

#include "ReconstructionHelpers.h"

/**********************************************************
 * Routine: Weighting_LS_Problem                          *
 *                                                        *
 * Assigns and computes the geometric weights for the     *
 * linear system.                                         *
 *                                                        *
 * The weights are equal to 1/d^2                         *
 * where "d" is the distance between the cell centers     *
 *********************************************************/
void Weighting_LS_Problem (DenseMatrix & LHS, DenseMatrix & RHS,
			   ColumnVector & DistanceCellCenters,
			   ColumnVector & W){

  /******************************************************************
   *  -- geometric weighting                                          *
   * returns the geometric weights W()                                *
   *******************************************************************/

  int M = LHS.size(0);
  int N_LHS = LHS.size(1);
  int N_RHS = RHS.size(1); // the number of rows in All_Delta_U is equal to M
  double tol = 1.0E-15; // it should be much less than the mesh size
  double SumWeights(0.0);
  double Distance;

  // compute the geometric weights based on the distance to each control volume
  for (int i=0; i<=M-1; ++i){
    W(i) = fabs(DistanceCellCenters(i));
    W(i) = 1.0/(tol + W(i)*W(i));
    SumWeights += W(i);
  }

  for (int i=0; i<=M-1; ++i){

    // normalize weights by their total sum
    W(i) /= SumWeights;

    // apply geometric weights to the matrix LHS & the free term RHS
    for (int j=0; j<=N_LHS-1; ++j){
      LHS(i,j) *= W(i);
    }
    for (int j=0; j<=N_RHS-1; ++j){
      RHS(i,j) *= W(i);
    }
  }
}

/************************************************
 * Routine: DetermineMaxDeltaSolutionStencil    *
 *                                              *
 * This subroutine determines the value of      *
 * MaxDeltaSolutionStencil for each parameter   *
 *                                              *
 * Input: MaxDeltaSolutionStencil vector        *
 *        (when the function is called it has   *
 *         the values of the MaxSolutionStencil)*
 *        MinSolutionStencil vector             *
 * Output: the result is written in             *
 *         MaxSolutionStencil                   *
 ***********************************************/
void DetermineMaxDeltaSolutionStencil(ColumnVector & MaxDeltaSolutionStencil,
				      ColumnVector & MinSolutionStencil,
				      const int & NumberOfVariables){

  double tol = 1.0E-13;  // jumps smaller than that won't be considered
  /* Determine the maximum solution value for each parameter */
  for (int parameter=0; parameter<NumberOfVariables; ++parameter){
    MaxDeltaSolutionStencil(parameter) -= MinSolutionStencil(parameter);
    if (MaxDeltaSolutionStencil(parameter) <= tol)
      MaxDeltaSolutionStencil(parameter) = 0.0;
  }
}

/************************************************
 * Routine: DetermineMaxDeltaSolutionStencil    *
 *                                              *
 * This subroutine determines the value of      *
 * MaxDeltaSolutionStencil for each parameter   *
 *                                              *
 * Input: MaxDeltaSolutionStencil vector        *
 *        (when the function is called it has   *
 *         the values of the MaxSolutionStencil)*
 *        All_Delta_U matrix                    *
 * Output: the result is written in             *
 *         MaxSolutionStencil                   *
 ***********************************************/
void DetermineMaxDeltaSolutionStencil(ColumnVector & MaxDeltaSolutionStencil,
				      DenseMatrix & All_Delta_U,
				      const int & NumberOfVariables){

  double tol = 5.0E-13;
  /* Determine the maximum solution difference for each parameter */
  for (int parameter=0; parameter<NumberOfVariables; ++parameter){
    MaxDeltaSolutionStencil(parameter) = fabs(All_Delta_U(0,parameter));
    for(int CounterNeighbCell = 1; CounterNeighbCell<All_Delta_U.size(0); ++CounterNeighbCell)
      MaxDeltaSolutionStencil(parameter) = max(MaxDeltaSolutionStencil(parameter), fabs(All_Delta_U(CounterNeighbCell,parameter)));
  }
}


void Weighting_LS_Problem (DenseMatrix & A, ColumnVector & Weights,
			   const ColumnVector & DistanceCellCenters,
			   ColumnVector & Delta_U, const ColumnVector & Delta_U_Original,
			   const double & MaxDeltaSolutionStencil,
			   const double &MaxDeltaSolutionOverDomain,
			   const double & CharacteristicLength,
			   int & ReconstructionOrder, int & FinalOrder, const double & CutoffKnob){

  /*************************************************************
   *  -- data-dependent weighting WITHOUT the residual from DI proposed by Lucian*
   ************************************************************/

  int M = A.size(0)-1;
  int N = A.size(1)-1;
  double tol = 1.0E-13;
  double MaxWeight = 0.0;
  double MaxDistance = fabs(DistanceCellCenters(0));
  double MinDistance = fabs(DistanceCellCenters(0));
  double NormalizedMaxDeltaSolutionStencil = MaxDeltaSolutionStencil/MaxDeltaSolutionOverDomain;
  double NormalizedDeltaSolution;
  double NormalizedDistance;
  double IS; 			// smoothness indicator for a control volume
  double eps;

  switch(ReconstructionOrder){
  case 1:
    eps = 0.01;
    break;
  case 2:
    eps = 0.001;
    break;
  case 3:
    eps = 0.0001;
    break;
  default:
    eps = 1.0e-5;
  }

  // determine MinDistance and MaxDistance
  for (int i=1; i<=M; ++i){
    MaxDistance = max(MaxDistance, fabs(DistanceCellCenters(i)));
    MinDistance = min(MinDistance, fabs(DistanceCellCenters(i)));
  }

  // Compute the data-dependent weight associated with each control volume

  /*******************************************************************************
   *************************  DIFFERENT WEIGTINGS  *******************************
   *******************************************************************************/


  /* %%%%%%%%%%%%%%%%%%%%% FIRST WEIGHTING %%%%%%%%%%%%%%%%%%%%%%*/
  
  //  if (MaxDeltaSolutionStencil > 0){
  /*     for (int i=0; i<=M; ++i){ */
	
  /*       // Normalize the parameters */
  /*       NormalizedDeltaSolution = fabs(Delta_U_Original(i))/MaxDeltaSolutionOverDomain; */
      
  /*       NormalizedDistance = DistanceCellCenters(i)/CharacteristicLength; */
      
  /*       //Weights(i) = (NormalizedMaxDeltaSolutionStencil * NormalizedDistance)/(NormalizedDeltaSolution + NormalizedDistance); */

  /*       Weights(i) = 1.0/(NormalizedDeltaSolution/NormalizedDistance + 1); */
  /*       Weights(i) = pow(Weights(i),ReconstructionOrder+1); */
      
  /*       // Determine the maximum weight */
  /*       MaxWeight = max(MaxWeight,Weights(i)); */
  /*       //       cout << "Final W=" << Weights(i) << endl; */
  /*     } */
  /*   } else { */
  /*     for (int i=0; i<=M; ++i){ */
  /*       Weights(i) = 0.0; */
  /*     } */
  /*   } */
    
  /* %%%%%%%%%%%%%%%%%%%%% SECOND WEIGHTING %%%%%%%%%%%%%%%%%%%%%%*/
  /*     if (MaxDeltaSolutionStencil > tol){ */
      
  /*       // Normalize the parameters */
  /*       Weights(i) = MaxDeltaSolutionStencil*fabs(DistanceCellCenters(i))/CharacteristicLength; */
  /*       Weights(i) /= fabs(Delta_U_Original(i))+fabs(DistanceCellCenters(i))/CharacteristicLength; */
  /*       Weights(i) = pow(Weights(i),ReconstructionOrder+1); */
      
  /* %%%%%%%%%%%%%%%%%%%%% THIRD WEIGHTING %%%%%%%%%%%%%%%%%%%%%%*/

  /*   if (MaxDeltaSolutionStencil > 0){ */
  /*     for (int i=0; i<=M; ++i){ */
  /*       IS = (Delta_U_Original(i)/DistanceCellCenters(i))/(MaxDeltaSolutionStencil/MaxDistance); */
  /*       //      IS *= IS; */
  /*       //IS = fabs(IS); */
  /*       //             Print(MaxDeltaSolutionStencil); */
  /*       //             Print(Delta_U_copy(i)); */
  /*       //             Print(Delta_XC(i)); */
  /*       //      Print(IS); */
  /*       //      IS = max(IS,eps); */
  /*       Weights(i) = 1/pow(IS,ReconstructionOrder+1); */
  
  /*       //      Determine the maximum weight */
  /*       MaxWeight = max(MaxWeight,Weights(i)); */
  /*     } */
  /*   } else { */
  /*     for (int i=0; i<=M; ++i){ */
  /*       Weights(i) = 0.0; */
  /*     } */
  /*   } */
  

  /* %%%%%%%%%%%%%%%%%%%%% FOURTH WEIGHTING %%%%%%%%%%%%%%%%%%%%%%*/
  if(MaxDeltaSolutionOverDomain > 0){
    for (int i=0; i<=M; ++i){
      Weights(i) = pow(fabs(Delta_U_Original(i))/MaxDeltaSolutionOverDomain,2.0) + pow(DistanceCellCenters(i)/CharacteristicLength,2.0);
      Weights(i) = 1.0/sqrt(Weights(i));
      Weights(i) = pow(Weights(i),ReconstructionOrder+1);
	
      /*        Print(fabs(Delta_U_Original(i))/MaxDeltaSolutionOverDomain); */
      /*        Print(DistanceCellCenters(i)/CharacteristicLength); */
      
      // Determine the maximum weight
      MaxWeight = max(MaxWeight,Weights(i));
      //       cout << "Final W=" << Weights(i) << endl;
    }
  } else {
    for (int i=0; i<=M; ++i){
      Weights(i) = 1.0;
    }
  }


  /*******************************************************************************
   *************************  DIFFERENT WEIGTINGS  *******************************
   *******************************************************************************/


  /* Normalize the weights */
  if(MaxWeight > 0){
    for (int i=0; i<=M; ++i){
      Weights(i) /= MaxWeight;
    }
  }

  // Determine the order of the reconstruction (number of computed derivatives)
  // based on the analysis of the weights
  /*   if (MaxWeight > 0) */
  /*     AnalyzeWeights(A,Weights,MaxDistance/CharacteristicLength,ReconstructionOrder,FinalOrder, CutoffKnob); */


  // Assign the weights to control volumes
  for (int i=0; i<=M; ++i){
    for (int j=0; j<=N; ++j){
      A(i,j) *= Weights(i);
    }
    Delta_U(i) *= Weights(i);
  }
  
}

void Carl_Weighting_LS_Problem (DenseMatrix & A, ColumnVector & Weights,
				ColumnVector & DistanceCellCenters,
				ColumnVector & Delta_U, ColumnVector & Delta_U_Original,
				double & MaxDeltaSolutionOverDomain, double & CharacteristicLength,
				int & ReconstructionOrder, int & FinalOrder, const double & CutoffKnob){

  /***************************************************************************************************************
   *  -- data-dependent weighting proposed by Carl Ollivier-Gooch in                                             *
   *  "A high-order accurate unstructured mesh ENO scheme based on Data-Dependent Least-Squares Reconstruction"  *
   *  Obs. The weights are normalized by their sum.                                                              *
   **************************************************************************************************************/

  int i,j;
  double WeightsSum = 0.0;

  // Compute the data-dependent weight associated with each control volume
  for (i=0; i<A.size(0); ++i){

    /************************ FIRST WEIGHTING **************************/
    if(MaxDeltaSolutionOverDomain > 0){
      Weights(i)  = fabs(Delta_U_Original(i))/MaxDeltaSolutionOverDomain;
      Weights(i) += fabs(DistanceCellCenters(i))/CharacteristicLength;
      Weights(i)  = pow(Weights(i),ReconstructionOrder+1);
      Weights(i)  = 1.0/Weights(i);
    } else {
      Weights(i) = 1.0;
    }

    WeightsSum += Weights(i);
  }


  // Assign the weights to the control volume
  for (i=0; i<A.size(0); ++i){

    // Compute normalized weight
    Weights(i) /= WeightsSum;

    // LHS
    for (j=0; j<A.size(1); ++j){
      A(i,j) *= Weights(i);
    }

    // RHS
    Delta_U(i) *= Weights(i);
  }
}

/**********************************************************
 * Routine: AnalyzeWeights                                *
 *                                                        *
 * This subroutine determines the number of derivatives   *
 * that can be computed in the reconstruction process     *
 * based on the data-dependent weights.                   *
 * Criterium: Values less than a determined cutoff value  *
 *            are considered too small to support the     *
 *            reconstruction.                             *
 * The cutoff value is of the order of the truncation     *
 * error of the mesh. The truncation error is proportional*
 * to Delta_X^(Order+1)                                   *
 *                                                        *
 *********************************************************/
void AnalyzeWeights(DenseMatrix & A, ColumnVector & Weights, const double & MaxDistance, 
		    const int & ReconstructionOrder, int & FinalOrder, const double & CutoffKnob){

  double cutoff= 1.0e-5; 	/* minimum value */
  int WeightCounter = 0;
  int M = A.size(0)-1;
  int N = A.size(1)-1;

  /*************************************  The cutoff test     *****************************/

  // Determine the value of the "cutoff" (pow(MaxDistance, ReconstructionOrder+1))
  //  cutoff = CutoffKnob*max(cutoff,min(0.1, 2*pow((ReconstructionOrder)*MaxDistance,ReconstructionOrder+1)));
  cutoff = CutoffKnob*(ReconstructionOrder*5)*pow(MaxDistance,ReconstructionOrder+1);

  // Analyze the data-dependent weights and determine the max. number of derivatives that can be computed
  for (int i=0; i<=M; ++i){
    // Analyze weight
    if (Weights(i) >= cutoff){
      ++WeightCounter;
      Weights(i) = 1.0e-12;
    }
  }

  // Reduce the order of reconstruction if not all the derivatives can be computed
  if (WeightCounter < M+1){	// that is: there are CVs with weights smaller than the value of the cutoff

    switch(ReconstructionOrder){
    case 3:			/* eliminate the derivatives related to the 3rd Order */
      if(WeightCounter > 19){
	FinalOrder = min(FinalOrder,3);
	break;			/* there are still enough CVs to support the 3rd-order reconstruction */
      } else {
	/* there are not enough CVs and therefore drop the order of the reconstruction to 2 */
	for (int i=0; i<=M; ++i){ /* there are 4 derivatives related to the 3rd-order */
	    A(i,2) = 0.0; 	/* third column D_03 */
	    A(i,5) = 0.0;	/* sixth column D_12 */
	    A(i,7) = 0.0;	/* eighth column D_21 */
	    A(i,8) = 0.0;	/* ninth column D_30 */
	}
      }

      if(WeightCounter > 16){
	FinalOrder = min(FinalOrder,2);
	break;                  /* there are still enough CVs to support the 2nd-order reconstruction */
      } else {
	/* there are not enough CVs and therefore drop the order of the reconstruction to 1 */
	for (int i=0; i<=M; ++i){ /* there are 3 derivatives related to the 2nd-order */
	    A(i,1) = 0.0;	/* second column D_02 */
	    A(i,4) = 0.0; 	/* fifth column D_11 */
	    A(i,6) = 0.0;	/* seventh column D_20 */
	}
	/* drop the number of points to have a compact stencil */
/* 	Weights(0) = 0.0; */
/* 	Weights(1) = 0.0; */
/* 	Weights(2) = 0.0; */
/* 	Weights(3) = 0.0; */
/* 	Weights(4) = 0.0; */
/* 	Weights(5) = 0.0; */
/* 	Weights(9) = 0.0; */
/* 	Weights(10) = 0.0; */
/* 	Weights(13) = 0.0; */
/* 	Weights(14) = 0.0; */
/* 	Weights(18) = 0.0; */
/* 	Weights(19) = 0.0; */
/* 	Weights(20) = 0.0; */
/* 	Weights(21) = 0.0; */
/* 	Weights(22) = 0.0; */
/* 	Weights(23) = 0.0; */
      }

      if(WeightCounter > 16){
	FinalOrder = min(FinalOrder,1);
	break;                  /* there are still enough CVs to support the 1st-order reconstruction */
      } else {
	/* there are not enough CVs and therefore drop the order of the reconstruction to 0 */
	for (int i=0; i<=M; ++i){ /* there are 2 derivatives related to the 1st-order */
	    A(i,0) = 0.0;	/* first column D_01 */
	    A(i,3) = 0.0;	/* fourth column D_10 */
	  }
	FinalOrder = 0;
      }
      break;

    case 2:                     /* eliminate the derivatives related to the 2nd Order */

      if(WeightCounter > 16){
	FinalOrder = min(FinalOrder,2);
	break;                  /* there are still enough CVs to support the 2nd-order reconstruction */
      } else {
	/* there are not enough CVs and therefore drop the order of the reconstruction to 1 */
	for (int i=0; i<=M; ++i){ /* there are 3 derivatives related to the 2nd-order */
	    A(i,1) = 0.0;	/* second column D_02 */
	    A(i,3) = 0.0; 	/* fourth column D_11 */
	    A(i,4) = 0.0;	/* fifth column D_20 */
	  }
      }

      if(WeightCounter > 6){
	FinalOrder = min(FinalOrder,1);
	break;                  /* there are still enough CVs to support the 1st-order reconstruction */
      } else {
	/* there are not enough CVs and therefore drop the order of the reconstruction to 0 */
	for (int i=0; i<=M; ++i){ /* there are 2 derivatives related to the 1st-order */
	    A(i,0) = 0.0;	/* first column D_01 */
	    A(i,2) = 0.0;	/* fourth column D_10 */
	  }
	FinalOrder = 0;
      }

      break;

    case 1:                     /* eliminate the derivatives related to the 1st Order */

      if(WeightCounter > 6){
	FinalOrder = min(FinalOrder,1);
	break;                  /* there are still enough CVs to support the 1st-order reconstruction */
      } else {
	/* there are not enough CVs and therefore drop the order of the reconstruction to 0 */
	for (int i=0; i<=M; ++i){ /* there are 2 derivatives related to the 1st-order */
	    A(i,0) = 0.0;	/* first column D_01 */
	    A(i,1) = 0.0;	/* second column D_10 */
	  }
	FinalOrder = 0;
      }

      break;
    }
  }

}

/* Determines the distance between the cell centers */
void DetermineDistanceCellCenters(ColumnVector & DistanceCellCenters, vector<Vector2D> & DeltaCellCenters){

  for (int i=0; i<DistanceCellCenters.size(); ++i){
    DistanceCellCenters(i)= sqrt(DeltaCellCenters[i].x*DeltaCellCenters[i].x + DeltaCellCenters[i].y*DeltaCellCenters[i].y);
  }
}

/* Determines the attained order of reconstruction based on the number of computed derivatives (in 2D) */
int InverseOrder2D(const int & ReconstructionOrder, const int & ComputedDerivatives){

  int FinalOrder = ReconstructionOrder;
  int ND = NumberOfDerivatives2D(ReconstructionOrder);	/* number of derivatives */
  int Difference;

  while (ND > 1){
    Difference = ND - ComputedDerivatives;
    if (Difference > 0){
      --FinalOrder;
      ND = NumberOfDerivatives2D(FinalOrder);
    }
    else{
      break;
    }
  }

  return FinalOrder;
}
