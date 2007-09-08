/**********************************************************************
 * LevelSet2DQuadHamiltonJacobiSingleBlock.cc                         *
 *                                                                    *
 * Single-block versions of subroutines for the solution of the 2D    *
 * Hamilton-Jacobi-type equations for the 2D Level Set multi-block    *
 * quadrilateral mesh solution classes.                               *
 *                                                                    *
 **********************************************************************/

// Include 2D LevelSet quadrilateral mesh solution header file.

#ifndef _LEVELSET2D_QUAD_INCLUDED
#include "LevelSet2DQuad.h"
#endif // _LEVELSET2D_QUAD_INCLUDED

/**********************************************************************
 * LevelSet2D_Quad_Block -- Hamilton-Jacobi Single Block External     *
 *                          Subroutines.                              *
 **********************************************************************/

/**********************************************************************
 * Routine: CFL_Hamilton_Jacobi                                       *
 *                                                                    *
 * Determines the allowable global and local time steps (for explicit *
 * Euler time stepping scheme) for the specified quadrilateral        *
 * solution block according to the Courant-Friedrichs-Lewy condition. *
 *                                                                    *
 **********************************************************************/
double CFL_Hamilton_Jacobi(LevelSet2D_Quad_Block &SolnBlk) {

  double dtMin = MILLION, d_i, d_j, w_i, w_j, v_i, v_j, f;
  double dt_d, dt_v, dt_f;
  Vector2D W, dpsi;

  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      if (i < SolnBlk.ICl || i > SolnBlk.ICu || j < SolnBlk.JCl || j > SolnBlk.JCu) {
	SolnBlk.dt[i][j] = ZERO;
      } else {
	d_i = TWO*(SolnBlk.Grid.Cell[i][j].A/(SolnBlk.Grid.lfaceE(i,j)+SolnBlk.Grid.lfaceW(i,j)));
	d_j = TWO*(SolnBlk.Grid.Cell[i][j].A/(SolnBlk.Grid.lfaceN(i,j)+SolnBlk.Grid.lfaceS(i,j)));
	dt_d = min(d_i,d_j);
	f = fabs(SolnBlk.U[i][j].F);
	dt_f = min(d_i,d_j)/(fabs(f) + TOLER);
	v_i = HALF*(SolnBlk.U[i][j].V*(SolnBlk.Grid.nfaceE(i,j)-SolnBlk.Grid.nfaceW(i,j)));
	v_j = HALF*(SolnBlk.U[i][j].V*(SolnBlk.Grid.nfaceN(i,j)-SolnBlk.Grid.nfaceS(i,j)));
	dt_v = min(d_i/(fabs(v_i)+TOLER),d_j/(fabs(v_j)+TOLER));
	SolnBlk.dt[i][j] = min(dt_d,min(dt_f,dt_v));
	dtMin = min(dtMin,SolnBlk.dt[i][j]);
      }
    }
  }
  
  // Return the global time step.
  return dtMin;

}

/**********************************************************************
 * Routine: dUdt_Multistage_Hamilton_Jacobi                           *
 *                                                                    *
 * This routine determines the solution residuals for a given stage   *
 * of a variety of multi-stage explicit time integration schemes for  *
 * a given solution block.                                            *
 *                                                                    *
 **********************************************************************/
int dUdt_Multistage_Hamilton_Jacobi(LevelSet2D_Quad_Block &SolnBlk,
				    const int i_stage,
				    LevelSet2D_Input_Parameters &IP) {

  int k_residual;
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

  /* Perform the biased reconstruction within each cell of the
     computational grid for this stage. */
  switch(IP.i_Reconstruction) {
  case RECONSTRUCTION_LINEAR_ESSENTIALLY_NON_OSCILLATORY :
  case RECONSTRUCTION_QUADRATIC_ESSENTIALLY_NON_OSCILLATORY :
  case RECONSTRUCTION_CUBIC_ESSENTIALLY_NON_OSCILLATORY :
    Reconstruction_EssentiallyNonOscillatory(SolnBlk,1,IP.i_Reconstruction);
    break;
  case RECONSTRUCTION_WEIGHTED_ESSENTIALLY_NON_OSCILLATORY :
    Reconstruction_WeightedEssentiallyNonOscillatory(SolnBlk,1);
    break;
  default:
    Reconstruction_EssentiallyNonOscillatory(SolnBlk,1,IP.i_Reconstruction);
    break;
  };

  /* Perform the reconstruction of the curvature within each cell of the
     computational grid for this stage. */
  Reconstruction_Curvature(SolnBlk,1);

  /* Evaluate the time rate of change of the solution (i.e., the 
     solution residuals). */

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
	};
      }

      if (j >= SolnBlk.JCl && j <= SolnBlk.JCu &&
	  i >= SolnBlk.ICl && i <= SolnBlk.ICu) {

	// Set the gradient of the level set function:
	SolnBlk.dUdx[i][j].psi = sqrt(max(sqr(max(SolnBlk.dUdxm[i][j].psi,ZERO)),sqr(min(SolnBlk.dUdxp[i][j].psi,ZERO))) +
				      max(sqr(max(SolnBlk.dUdym[i][j].psi,ZERO)),sqr(min(SolnBlk.dUdyp[i][j].psi,ZERO))));
	SolnBlk.dUdy[i][j].psi = sqrt(max(sqr(min(SolnBlk.dUdxm[i][j].psi,ZERO)),sqr(max(SolnBlk.dUdxp[i][j].psi,ZERO))) +
				      max(sqr(min(SolnBlk.dUdym[i][j].psi,ZERO)),sqr(max(SolnBlk.dUdyp[i][j].psi,ZERO))));

	// Add the front propagation contribution to the level set function solution residual.
 	SolnBlk.dUdt[i][j][k_residual].psi -= (IP.Hamilton_Jacobi_CFL_Number*SolnBlk.dt[i][j])*
	                                      (max(SolnBlk.U[i][j].F,ZERO)*SolnBlk.dUdx[i][j].psi + 
					       min(SolnBlk.U[i][j].F,ZERO)*SolnBlk.dUdy[i][j].psi);

	// Add the curvature-based flow contribution to the level set function solution residual.
	SolnBlk.dUdt[i][j][k_residual].psi -= (IP.Hamilton_Jacobi_CFL_Number*SolnBlk.dt[i][j])*
	  (max(-IP.Curvature_Motion*SolnBlk.kappa[i][j].psi,ZERO)*SolnBlk.dUdx[i][j].psi + min(-IP.Curvature_Motion*SolnBlk.kappa[i][j].psi,ZERO)*SolnBlk.dUdy[i][j].psi);

	// Add the convective/bulk-velocity flow contribution to the
	// level set function solution residual.
 	SolnBlk.dUdt[i][j][k_residual].psi -= (IP.Hamilton_Jacobi_CFL_Number*SolnBlk.dt[i][j])*
	                                      (max(SolnBlk.U[i][j].V.x,ZERO)*SolnBlk.dUdxm[i][j].psi +
 					       min(SolnBlk.U[i][j].V.x,ZERO)*SolnBlk.dUdxp[i][j].psi +
 					       max(SolnBlk.U[i][j].V.y,ZERO)*SolnBlk.dUdym[i][j].psi +
 					       min(SolnBlk.U[i][j].V.y,ZERO)*SolnBlk.dUdyp[i][j].psi);

	// Add the curvature-based flow contribution to the
	// level set function solution residual.
// 	SolnBlk.dUdt[i][j][k_residual].psi += (IP.Hamilton_Jacobi_CFL_Number*SolnBlk.dt[i][j])*
// 	                                      IP.Curvature_Motion*SolnBlk.kappa[i][j].psi*
	
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

/**********************************************************************
 * Routine: Update_Solution_Multistage_Hamilton_Jacobi                *
 *                                                                    *
 * This routine updates solution states of the given solution block   *
 * for a variety of multi-stage explicit time integration schemes.    *
 *                                                                    *
 **********************************************************************/
int Update_Solution_Multistage_Hamilton_Jacobi(LevelSet2D_Quad_Block &SolnBlk,
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
      SolnBlk.U[i][j].psi = SolnBlk.Uo[i][j].psi + omega*SolnBlk.dUdt[i][j][k_residual].psi;
    }
  }

  // Solution successfully updated.
  return 0;
    
}
