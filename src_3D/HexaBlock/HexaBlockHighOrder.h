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


//  // SET VARIABLES USED IN THE RESIDUAL CALCULATION PROCESS
//
//  int i, j, k_residual;
//
//  /* Evaluate the solution residual for stage 
//     i_stage of an N stage scheme. */
//
//  /* Evaluate the time step fraction and residual storage location for the stage. */
//  
//  switch(IP.i_Time_Integration) {
//  case TIME_STEPPING_EXPLICIT_EULER :
//    k_residual = 0;
//    break;
//  case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
//    k_residual = 0;
//    break;
//  case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
//    k_residual = 0;
//    if (IP.N_Stage == 4) {
//      if (i_stage == 4) {
//	k_residual = 0;
//      } else {
//	k_residual = i_stage - 1;
//      } /* endif */
//    } /* endif */
//    break;
//  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
//    k_residual = 0;
//    break;
//  default:
//    k_residual = 0;
//  } /* endswitch */
//
//
//  // ************* Step 1. (Re)-Set parameters in all affected cells based on the time integration scheme **************
//  // *******************************************************************************************************************
//  for ( j = JCl-1 ; j <= JCu+1 ; ++j ){
//    for ( i = ICl-1 ; i <= ICu+1 ; ++i ) {
//
//      if ( i_stage == 1 ){
//	Uo[i][j] = U[i][j];
//	dUdt[i][j][k_residual].Vacuum();  // set to zero
//      } else {
//	switch(IP.i_Time_Integration) {
//	case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
//	  // 
//	  break;
//	case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
//	  if (IP.N_Stage == 2) {
//	    // 
//	  } else if (IP.N_Stage == 4 && i_stage == 4) {
//	    dUdt[i][j][k_residual] = ( dUdt[i][j][0] + 
//				       TWO*dUdt[i][j][1] +
//				       TWO*dUdt[i][j][2] );
//	  } else {
//	    dUdt[i][j][k_residual].Vacuum();  // set to zero
//	  } /* endif */
//	  break;
//	case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
//	  dUdt[i][j][k_residual].Vacuum(); // set to zero
//	  break;
//	default:
//	  dUdt[i][j][k_residual].Vacuum(); // set to zero
//	  break;
//	} /* endswitch */
//      }/* endif */
//
//    } // endfor (i)
//  } // endfor (j)
//
//
//  // ** Step 2. Compute high-order spatial residual for the current time step fraction **
//  // ************************************************************************************
//  return dUdt_Residual_HighOrder(IP, k_residual, true, Pos);

  // Temporarily return the original Euler3D peicewise linear function:
  // ------------------------------------------------------------------
  return dUdt_Multistage_Explicit(i_stage, IPs);

}
