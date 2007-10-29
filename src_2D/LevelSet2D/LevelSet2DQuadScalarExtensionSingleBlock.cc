/******************************************************************//**
 * \file LevelSet2DQuadScalarExtensionSingleBlock.cc                  
 *                                                                    
 * Scalar extension single-block versions of subroutines for 2D Level 
 * Set multi-block quadrilateral mesh solution classes.               
 *                                                                    
 **********************************************************************/

// Include 2D LevelSet quadrilateral mesh solution header file.

#ifndef _LEVELSET2D_QUAD_INCLUDED
#include "LevelSet2DQuad.h"
#endif // _LEVELSET2D_QUAD_INCLUDED

/**********************************************************************
 * LevelSet2D_Quad_Block -- Scalar Extension Single Block External    *
 *                          Subroutines.                              *
 **********************************************************************/

/******************************************************************//**
 * Routine: CFL_Scalar_Extension                                      
 *                                                                    
 * Determines the allowable global and local time steps (for explicit 
 * Euler time stepping scheme) for the specified quadrilateral        
 * solution block according to the Courant-Friedrichs-Lewy condition. 
 *                                                                    
 **********************************************************************/
double CFL_Scalar_Extension(LevelSet2D_Quad_Block &SolnBlk) {

  double dtMin = MILLION, d_i, d_j;

  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      if (i < SolnBlk.ICl || i > SolnBlk.ICu || j < SolnBlk.JCl || j > SolnBlk.JCu) {
	SolnBlk.dt[i][j] = ZERO;
      } else {
	d_i = TWO*(SolnBlk.Grid.Cell[i][j].A/(SolnBlk.Grid.lfaceE(i,j)+SolnBlk.Grid.lfaceW(i,j)));
	d_j = TWO*(SolnBlk.Grid.Cell[i][j].A/(SolnBlk.Grid.lfaceN(i,j)+SolnBlk.Grid.lfaceS(i,j)));
	SolnBlk.dt[i][j] = min(d_i,d_j);
	dtMin = min(dtMin,SolnBlk.dt[i][j]);
      }
    }
  }

  // Return the global time step.
  return dtMin;

}

/******************************************************************//**
 * Routine: dUdt_Multistage_Scalar_Extension                          
 *                                                                    
 * This routine evaulates the stage solution residual for the scalar  
 * (front speed) extension equation on the specified solution block.  
 *                                                                    
 **********************************************************************/
int dUdt_Multistage_Scalar_Extension(LevelSet2D_Quad_Block &SolnBlk,
				     const int i_stage,
				     LevelSet2D_Input_Parameters &IP) {

  int error_flag, k_residual;
  double omega;

  // Evaluate the solution residual for stage i_stage of n_stage 
  // scheme.

  // Evaluate the time step fraction and residual storage location
  // for the stage.   
  switch(IP.i_Time_Integration) {
  case TIME_STEPPING_EXPLICIT_EULER :
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    if (IP.N_Stage == 4) {
      if (i_stage == 4) {
	k_residual = 0;
      } else {
	k_residual = i_stage - 1;
      }
    }
    break;
  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
    omega = MultiStage_Optimally_Smoothing(i_stage,IP.N_Stage,IP.i_Limiter);
    k_residual = 0;
    break;
  default:
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    break;
  };

  // Perform the reconstruction within each cell of the
  // computational grid for this stage.
  Linear_Reconstruction_LeastSquares(SolnBlk,IP.i_Limiter);
  switch(IP.i_Reconstruction) {
  case RECONSTRUCTION_LINEAR_ESSENTIALLY_NON_OSCILLATORY :
  case RECONSTRUCTION_QUADRATIC_ESSENTIALLY_NON_OSCILLATORY :
  case RECONSTRUCTION_CUBIC_ESSENTIALLY_NON_OSCILLATORY :
    Reconstruction_EssentiallyNonOscillatory(SolnBlk,2,IP.i_Reconstruction);
    break;
  case RECONSTRUCTION_WEIGHTED_ESSENTIALLY_NON_OSCILLATORY :
    Reconstruction_WeightedEssentiallyNonOscillatory(SolnBlk,2);
    break;
  default:
    Reconstruction_EssentiallyNonOscillatory(SolnBlk,2,IP.i_Reconstruction);
    break;
  };

  // Determine the sign function.
  error_flag = Calculate_Sign_Function(SolnBlk,IP);
  if (error_flag) return error_flag;

  // Evaluate the time rate of change of the solution (i.e., the 
  // solution residuals) using a second-order upwind scheme.

  // Add i-direction (zeta-direction) fluxes.
  for (int j = SolnBlk.JCl-1; j <= SolnBlk.JCu+1; j++) {

    if (i_stage == 1) {
      SolnBlk.Uo[SolnBlk.ICl-1][j] = SolnBlk.U[SolnBlk.ICl-1][j];
      SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual] = LevelSet2D_ZERO;
    } else {
      SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual] = LevelSet2D_ZERO;
    }

    for (int i = SolnBlk.ICl-1; i <= SolnBlk.ICu; i++) {

      if (i_stage == 1) {
	SolnBlk.Uo[i+1][j] = SolnBlk.U[i+1][j];
	SolnBlk.dUdt[i+1][j][k_residual] = LevelSet2D_ZERO;
      } else if (j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
	switch(IP.i_Time_Integration) {
	case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
	  //SolnBlk.dUdt[i+1][j][k_residual] = SolnBlk.dUdt[i+1][j][k_residual];
	  break;
	case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
	  if (IP.N_Stage == 2) {
	    //SolnBlk.dUdt[i+1][j][k_residual] = SolnBlk.dUdt[i+1][j][k_residual];
	  } else if (IP.N_Stage == 4 && i_stage == 4) {
	    SolnBlk.dUdt[i+1][j][k_residual] = SolnBlk.dUdt[i+1][j][0] + 
                                               TWO*SolnBlk.dUdt[i+1][j][1] +
                                               TWO*SolnBlk.dUdt[i+1][j][2];
	  } else {
	    SolnBlk.dUdt[i+1][j][k_residual] = LevelSet2D_ZERO;
	  }
	  break;
	case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
	  SolnBlk.dUdt[i+1][j][k_residual] = LevelSet2D_ZERO;
	  break;
	default:
	  SolnBlk.dUdt[i+1][j][k_residual] = LevelSet2D_ZERO;
	  break;
	}
      }

      if (j >= SolnBlk.JCl && j <= SolnBlk.JCu &&
	  i >= SolnBlk.ICl && i <= SolnBlk.ICu) {

	// Update F solution residual.
	SolnBlk.dUdt[i][j][k_residual].F += (IP.Scalar_Extension_CFL_Number*SolnBlk.dt[i][j])*
                                            (max(SolnBlk.sign[i][j]*SolnBlk.dUdx[i][j].psi,ZERO)*SolnBlk.dUdxm[i][j].F +
					     min(SolnBlk.sign[i][j]*SolnBlk.dUdx[i][j].psi,ZERO)*SolnBlk.dUdxp[i][j].F +
					     max(SolnBlk.sign[i][j]*SolnBlk.dUdy[i][j].psi,ZERO)*SolnBlk.dUdym[i][j].F +
					     min(SolnBlk.sign[i][j]*SolnBlk.dUdy[i][j].psi,ZERO)*SolnBlk.dUdyp[i][j].F);

      }

      SolnBlk.dUdt[i][SolnBlk.JCl-1][k_residual] = LevelSet2D_ZERO;
      SolnBlk.dUdt[i][SolnBlk.JCu+1][k_residual] = LevelSet2D_ZERO;

    }

    if (j >= SolnBlk.JCl && j <= SolnBlk.JCu) {
      SolnBlk.dUdt[SolnBlk.ICl-1][j][k_residual] = LevelSet2D_ZERO;
      SolnBlk.dUdt[SolnBlk.ICu+1][j][k_residual] = LevelSet2D_ZERO;
    }

  }

  // Residual for the stage successfully calculated.
  return 0;

}

/******************************************************************//**
 * Routine: Update_Multistage_Scalar_Extension                        
 *                                                                    
 * This routine updates solution states of the given solution block   
 * for a variety of multi-stage explicit time integration schemes.    
 *                                                                    
 **********************************************************************/
int Update_Multistage_Scalar_Extension(LevelSet2D_Quad_Block &SolnBlk,
				       const int i_stage,
				       LevelSet2D_Input_Parameters &IP) {

  int k_residual;
  double omega;

  // Perform update of solution variables for stage i_stage of n_stage
  // scheme.

  // Evaluate the time step fraction and residual storage location for
  // the stage.
  switch(IP.i_Time_Integration) {
  case TIME_STEPPING_EXPLICIT_EULER :
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    if (IP.N_Stage == 4) {
      if (i_stage == 4) {
	k_residual = 0;
      } else {
	k_residual = i_stage - 1;
      }
    }
    break;
  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
    omega = MultiStage_Optimally_Smoothing(i_stage,IP.N_Stage,IP.i_Limiter);
    k_residual = 0;
    break;
  default:
    omega = Runge_Kutta(i_stage,IP.N_Stage);
    k_residual = 0;
    break;
  };

  // Update solution variables for this stage.
  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      SolnBlk.U[i][j].F = SolnBlk.Uo[i][j].F + omega*SolnBlk.dUdt[i][j][k_residual].F;
    }
  }

  // Solution successfully updated.
  return 0;

}
