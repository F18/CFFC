/**********************************************************************
 * LevelSet2DQuadEikonalSingleBlock.cc                                *
 *                                                                    *
 * Single-block versions of subroutines for the solution of the       *
 * Eikonal equations for the 2D Level Set multi-block quadrilateral   *
 * mesh solution classes.                                             *
 *                                                                    *
 **********************************************************************/

// Include 2D LevelSet quadrilateral mesh solution header file.

#ifndef _LEVELSET2D_QUAD_INCLUDED
#include "LevelSet2DQuad.h"
#endif // _LEVELSET2D_QUAD_INCLUDED

/**********************************************************************
 * LevelSet2D_Quad_Block -- Eikonal Single Block External Subroutines.*
 **********************************************************************/

/**********************************************************************
 * Routine: Eikonal_Error                                             *
 *                                                                    *
 * This routine calculates this block's area-weighted error in the    *
 * Eikonal equation solution. Local Error = | grad(phi) - 1 |         *
 *                                                                    *
 **********************************************************************/
int Eikonal_Error(LevelSet2D_Quad_Block &SolnBlk,
		  LevelSet2D_Input_Parameters &IP,
		  double &block_error,
		  double &block_area) {

  double grad;

  // Error in Eikonal equation at (i,j):
  double local_error;
  local_error = ZERO;

  // Perform the biased reconstruction within each cell of the
  // computational grid.
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

  // Compute the error.
  for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      grad = sqrt(max(sqr(max(SolnBlk.dUdxm[i][j].psi,ZERO)),
		      sqr(min(SolnBlk.dUdxp[i][j].psi,ZERO))) +
		  max(sqr(max(SolnBlk.dUdym[i][j].psi,ZERO)),
		      sqr(min(SolnBlk.dUdyp[i][j].psi,ZERO))));

      local_error = abs(grad-ONE);
      block_error = block_error + local_error*SolnBlk.Grid.Cell[i][j].A;
      block_area = block_area + SolnBlk.Grid.Cell[i][j].A;
    }
  }

  // Eikonal error calculation is successful.
  return 0;

}

/**********************************************************************
 * Routine: CFL_Eikonal                                               *
 *                                                                    *
 * Determines the allowable global and local time steps (for explicit *
 * Euler time stepping scheme) for the specified quadrilateral        *
 * solution block according to the Courant-Friedrichs-Lewy condition. *
 *                                                                    *
 **********************************************************************/
double CFL_Eikonal(LevelSet2D_Quad_Block &SolnBlk) {

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

/**********************************************************************
 * Routine: dUdt_Multistage_Eikonal                                   *
 *                                                                    *
 * This routine evaulates the stage solution residual for the Eikonal *
 * equation on the specified solution block.                          *
 *                                                                    *
 **********************************************************************/
int dUdt_Multistage_Eikonal(LevelSet2D_Quad_Block &SolnBlk,
			    const int i_stage,
			    LevelSet2D_Input_Parameters &IP) {

  int k_residual;
  double omega, dx;

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

  // Perform the reconstruction within each cell of the computational
  // grid for this stage.
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

  // Evaluate the time rate of change of the solution (i.e., the 
  // solution residuals).

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

	switch(IP.i_Eikonal_Scheme) {
	case EIKONAL_SCHEME_SUSSMAN :
	  /**********************************************************************
	   * Sussman, Smereka, Osher's reinitialization scheme.                 *
	   *                                                                    *
	   * See: Sussman, Smereka, Osher, "A level set approach for computing  *
	   *      solutions to incompressible two-phase flow", J. Comput. Phys. *
	   *      Vol. 114, pp. 146-159 (1994).                                 *
	   *                                                                    *
	   * The scheme is set up so that the gradient is calculated based on   *
	   * the user's choice of Reconstruction_Type (ie. WENO).               *
	   **********************************************************************/
	  // Update the Eikonal equation residual based on the Godunov
	  // scheme for a Hamilton-Jacobi type equation.
	  if (SolnBlk.sign[i][j] > ZERO) {
	    SolnBlk.dUdt[i][j][k_residual].psi += (SolnBlk.sign[i][j]*IP.Eikonal_CFL_Number*SolnBlk.dt[i][j])*
		                                  (ONE - sqrt(max(sqr(max(SolnBlk.dUdxm[i][j].psi,ZERO)),
								  sqr(min(SolnBlk.dUdxp[i][j].psi,ZERO))) +
							      max(sqr(max(SolnBlk.dUdym[i][j].psi,ZERO)),
								  sqr(min(SolnBlk.dUdyp[i][j].psi,ZERO)))));
	  } else {
	    SolnBlk.dUdt[i][j][k_residual].psi += (SolnBlk.sign[i][j]*IP.Eikonal_CFL_Number*SolnBlk.dt[i][j])*
		                                  (ONE - sqrt(max(sqr(min(SolnBlk.dUdxm[i][j].psi,ZERO)),
								  sqr(max(SolnBlk.dUdxp[i][j].psi,ZERO))) +
							      max(sqr(min(SolnBlk.dUdym[i][j].psi,ZERO)),
								  sqr(max(SolnBlk.dUdyp[i][j].psi,ZERO)))));
	  }
	  break;
	case EIKONAL_SCHEME_RUSSO_SMEREKA :
	  /**********************************************************************
	   * Russo and Smereka's new method for computing the signed distance   *
	   * function.                                                          *
	   *                                                                    *
	   * See: Russo, Smereka, "A Remark on Computing Distance Functions",   *
	   *      J. Comput. Phys. Vol. 163, pp. 51-67 (2000).                  *
	   **********************************************************************/
	  if ((SolnBlk.Uoo[i][j].psi*SolnBlk.Uoo[i-1][j].psi < 0) ||
	      (SolnBlk.Uoo[i][j].psi*SolnBlk.Uoo[i+1][j].psi < 0) ||
	      (SolnBlk.Uoo[i][j].psi*SolnBlk.Uoo[i][j-1].psi < 0) ||
	      (SolnBlk.Uoo[i][j].psi*SolnBlk.Uoo[i][j+1].psi < 0)) {
	    // The current cell is within one grid point from the level set.
	    dx = max(TOLER,min(fabs(SolnBlk.Grid.Cell[2][2].Xc.x-SolnBlk.Grid.Cell[1][2].Xc.x),
			       fabs(SolnBlk.Grid.Cell[2][2].Xc.y-SolnBlk.Grid.Cell[2][1].Xc.y)));


	    double D;
	    D = (TWO*SolnBlk.Uoo[i][j].psi)/sqrt(sqr(SolnBlk.Uoo[i+1][j].psi-SolnBlk.Uoo[i-1][j].psi) + sqr(SolnBlk.Uoo[i][j+1].psi-SolnBlk.Uoo[i][j-1].psi));

	    SolnBlk.dUdt[i][j][k_residual].psi -= IP.Eikonal_CFL_Number*SolnBlk.dt[i][j]*(SolnBlk.sign[i][j]*fabs(SolnBlk.U[i][j].psi)/dx-D);



// 	    SolnBlk.dUdt[i][j][k_residual].psi += (IP.Eikonal_CFL_Number*SolnBlk.dt[i][j])*(TWO*SolnBlk.Uoo[i][j].psi)/sqrt(sqr(SolnBlk.Uoo[i+1][j].psi-SolnBlk.Uoo[i-1][j].psi) + sqr(SolnBlk.Uoo[i][j+1].psi-SolnBlk.Uoo[i][j-1].psi));


// 	    if (SolnBlk.Uoo[i][j].psi < ZERO) {
// 	      SolnBlk.dUdt[i][j][k_residual].psi += (IP.Eikonal_CFL_Number*SolnBlk.dt[i][j])*fabs(SolnBlk.U[i][j].psi)/dx;
// 	    } else if (SolnBlk.Uoo[i][j].psi > ZERO) {
// 	      SolnBlk.dUdt[i][j][k_residual].psi -= (IP.Eikonal_CFL_Number*SolnBlk.dt[i][j])*fabs(SolnBlk.U[i][j].psi)/dx;
// 	    }

	  } else {
	    // The current cell is not within one grid point from the
	    // level set.  Update the Eikonal equation residual based on
	    // the Godunov scheme for a Hamilton-Jacobi type equation.
	    if (SolnBlk.sign[i][j] > ZERO) {
	      SolnBlk.dUdt[i][j][k_residual].psi += (SolnBlk.sign[i][j]*IP.Eikonal_CFL_Number*SolnBlk.dt[i][j])*
		                                    (ONE - sqrt(max(sqr(max(SolnBlk.dUdxm[i][j].psi,ZERO)),
								    sqr(min(SolnBlk.dUdxp[i][j].psi,ZERO))) +
								max(sqr(max(SolnBlk.dUdym[i][j].psi,ZERO)),
								    sqr(min(SolnBlk.dUdyp[i][j].psi,ZERO)))));
	    } else {
	      SolnBlk.dUdt[i][j][k_residual].psi += (SolnBlk.sign[i][j]*IP.Eikonal_CFL_Number*SolnBlk.dt[i][j])*
		                                    (ONE - sqrt(max(sqr(min(SolnBlk.dUdxm[i][j].psi,ZERO)),
								    sqr(max(SolnBlk.dUdxp[i][j].psi,ZERO))) +
								max(sqr(min(SolnBlk.dUdym[i][j].psi,ZERO)),
								    sqr(max(SolnBlk.dUdyp[i][j].psi,ZERO)))));
	    }
	  }
	  break;
	};

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
 * Routine: Update_Multistage_Eikonal                                 *
 *                                                                    *
 * This routine updates solution states of the given solution block   *
 * for a variety of multi-stage explicit time integration schemes.    *
 *                                                                    *
 **********************************************************************/
int Update_Multistage_Eikonal(LevelSet2D_Quad_Block &SolnBlk,
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

