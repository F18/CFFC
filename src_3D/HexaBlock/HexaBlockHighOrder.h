/*! \file HexaBlockHighOrder.h
  @brief High-order Subroutines for 3D Euler Equations on Hexahedral Grids */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#ifndef _HEXA_BLOCK_INCLUDED
#include "HexaBlock.h"
#endif //_HEXA_BLOCK_INCLUDED


//! Routine: dUdt_Multistage_Explicit_HighOrder
//  --------------------------------------------------------
/*!
 * This routine determines the solution residuals for a 
 * given stage of a variety of multi-stage explicit     
 * time integration schemes for the solution block.
 * The solution residuals are evaluated  
 * using a high-order CENO upwind finite-volume 
 * spatial discretization scheme.
 *
 *///-------------------------------------------------------
template<class SOLN_pSTATE, class SOLN_cSTATE>
int Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::
dUdt_Multistage_Explicit_HighOrder(const int i_stage,
				   Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs) {

  int i, j, k,  k_residual;
  double omega; 
  Vector3D dX;
   
  SOLN_pSTATE Wl, Wr;
  SOLN_cSTATE Flux;

  SOLN_cSTATE U_VACUUM;
  U_VACUUM.Vacuum();
  SOLN_pSTATE W_VACUUM;
  W_VACUUM.Vacuum();
 

  /* Evaluate the solution residual for stage 
   *       i_stage of an N stage scheme. */

  /* Evaluate the time step fraction and residual storage location 
   *       for the stage. */

  switch(IPs.i_Time_Integration) {
  case TIME_STEPPING_EXPLICIT_EULER :
    omega = Runge_Kutta(i_stage, IPs.N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
    omega = Runge_Kutta(i_stage, IPs.N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
    omega = Runge_Kutta(i_stage, IPs.N_Stage);
    k_residual = 0; 
    if (IPs.N_Stage == 4) { 
      if (i_stage == 4) { 
	k_residual = 0; 
      } else { 
	k_residual = i_stage - 1; 
      } // endif  
    } // endif 
    break;
  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
    omega = MultiStage_Optimally_Smoothing(i_stage, 
					   IPs.N_Stage,
					   IPs.i_Limiter);
    k_residual = 0;
    break;
  default: 
    omega = Runge_Kutta(i_stage, IPs.N_Stage);
    k_residual = 0;
    break;
  } /* endswitch */

  // ************* Step 1. (Re)-Set parameters in all affected cells based on the time integration scheme **************
  // *******************************************************************************************************************
  for ( k = KCl-1 ; k <= KCu+1 ; ++k ){
    for ( j = JCl-1 ; j <= JCu+1 ; ++j ){
      for ( i = ICl-1 ; i <= ICu+1 ; ++i ) {

	if ( i_stage == 1 ){
	  Uo[i][j][k] = U[i][j][k];
	  dUdt[i][j][k][k_residual].Vacuum();  // set to zero
	} else {
	  switch(IPs.i_Time_Integration) {
	  case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
	    // 
	    break;
	  case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
	    if (IPs.N_Stage == 2) {
	      // 
	    } else if (IPs.N_Stage == 4 && i_stage == 4) {
	      dUdt[i][j][k][k_residual] = ( dUdt[i][j][k][0] + 
					 TWO*dUdt[i][j][k][1] +
					 TWO*dUdt[i][j][k][2] );
	    } else {
	      dUdt[i][j][k][k_residual].Vacuum();  // set to zero
	    } /* endif */
	    break;
	  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
	    dUdt[i][j][k][k_residual].Vacuum(); // set to zero
	    break;
	  default:
	    dUdt[i][j][k][k_residual].Vacuum(); // set to zero
	    break;
	  } /* endswitch */
	}/* endif */

      } // endfor (i)
    } // endfor (j)
  } // endfor (k)

  // ** Step 2. Compute high-order spatial residual for the current time step fraction **
  // ************************************************************************************
  // return dUdt_Residual_HighOrder(IPs, k_residual, true);
  dUdt_Residual_HighOrder(IPs, k_residual, true);

  // Temporarily return the original Euler3D peicewise linear function:
  // ------------------------------------------------------------------
  return dUdt_Multistage_Explicit(i_stage, IPs);

}

/*!
 * Evaluate the residual for the solution block 
 * using the high-order CENO upwind finite-volume 
 * spatial discretization scheme.
 * The residual is stored in dUdt[][][][k_residual].
 *
 * \param IPs  input parameters object
 * \param k_residual index to identify the residual storage location
 *
 */
template<class SOLN_pSTATE, class SOLN_cSTATE>
int Hexa_Block<SOLN_pSTATE, SOLN_cSTATE>::
dUdt_Residual_HighOrder(Input_Parameters<SOLN_pSTATE, SOLN_cSTATE> &IPs,
			const int & k_residual,
			const bool & UseTimeStep) {

   // SET VARIABLES USED IN THE RESIDUAL CALCULATION PROCESS

  int i, j, k, GQPoint, Position, SplineSegment;
  bool IsNonSmoothHighOrderReconstruction;
  SOLN_pSTATE Wl, Wr, W_face;
  SOLN_cSTATE Flux, FaceFlux;
  int NumGQP(Grid.getNumGQP());	  // Number of Gauss quadrature points per face used to compute the flux integral

  Vector3D *GaussQuadPoints = new Vector3D [NumGQP]; // the GQPs at which a Riemann-like problem is solved
  double * GaussQuadWeights = new double [NumGQP];   // the Gauss integration weights for each Gauss quadrature

  /* Set the GaussQuadWeights. */
  GaussQuadratureData::getGaussQuadWeights(GaussQuadWeights, NumGQP);

  /* Evaluate the solution residual 
     and write it to dUdt[][][][k_residual]. */

  /***************************************************************************************
   *                 EVALUATE THE HIGH-ORDER SOLUTION RESIDUALS                          *
   *                                                                                     *
   * Algorithm Purpose: To evaluate solution residuals for solution blocks               *
   *                    characterized by a broad range of options.                       *
   *                                                                                     *
   * Important options to consider:                                                      *
   *         --> Geometry treatment: high-order or low-order                             *
   *         --> Spatial accuracy:   order of accuracy for flux calculation              *
   *         --> Boundary flux calculation: 'Riemann' problem or reconstruction based    *
   *                                                                                     *
   * In order to respond easier to all these parameter variations, the following         *
   * algorithm is adopted to sweep through the cell interfaces:                          *
   *         --> Compute all fluxes at interior inter-cellular faces.                    *
   *         --> Compute fluxes for North, South, East, West, Top, and Bottom            *
   *             block boundary faces.                                                   *
   *                                                                                     *
   ***************************************************************************************/

  /* Evaluate the time rate of change of the solution
     (i.e., the solution residuals) using a high-order
     CENO upwind finite-volume scheme. */

  /* Perform the high-order CENO reconstruction within
     each cell of the computational grid for this stage.
     NOTE: This solution reconstruction enforces monotonicity if required so!
  */
 








  // Deallocate memory
  delete [] GaussQuadPoints;
  delete [] GaussQuadWeights;
  
  /* residual for the stage successfully calculated. */
  return (0);
  
}
