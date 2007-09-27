/**********************************************************************
 * NavierStokes2DQuadTurbulenceSingleBlock.cc: Single-block versions  *
 *                                             of turbulence          *
 *                                             subroutines for 2D     *
 *                                             Navier-Stokes multi-   *
 *                                             block quadrilateral    *
 *                                             mesh solution classes. *
 **********************************************************************/

// Include 2D NavierStokes quadrilateral mesh solution header file.

#ifndef _NAVIERSTOKES2D_QUAD_INCLUDED
#include "NavierStokes2DQuad.h"
#endif // _NAVIERSTOKES2D_QUAD_INCLUDED

/**********************************************************************
 * NavierStokes2D_Quad_Block -- Turbulence Single Block External      *
 *                              Subroutines.                          *
 **********************************************************************/

/**********************************************************************
 * Routine: Turbulence_ICs                                            *
 *                                                                    *
 * Assigns initial turbulence conditions and data to the solution     *
 * variables of the specified quadrilateral solution block.           *
 *                                                                    *
 **********************************************************************/
int Turbulence_ICs(NavierStokes2D_Quad_Block &SolnBlk,
		   NavierStokes2D_Input_Parameters &IP,
		   NavierStokes2D_pState *Wo) {

  // Exit immediately if not a turbulent flow.
  if (SolnBlk.Flow_Type != FLOWTYPE_TURBULENT_RANS_K_OMEGA) return 0;

//   if (IP.i_ICs == IC_UNIFORM && IP.i_Grid == GRID_ROCKET_MOTOR) {
//     double length = IP.Chamber_Length+IP.Nozzle_Length;
//     Vector2D twall, that, v;
//     for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
//       for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
// 	SolnBlk.W[i][j] = TurbulentPipeFlow(IP.Wo,
// 					    SolnBlk.Grid.Cell[i][j].Xc,
// 					    IP.dp,
// 					    length,
// 					    SolnBlk.Wall[i][j].Xwall.y);
// 	v = SolnBlk.W[i][j].v;
// 	twall.x = -SolnBlk.Wall[i][j].nwall.y;
// 	twall.y =  SolnBlk.Wall[i][j].nwall.x;
// 	that = twall*(SolnBlk.Grid.Cell[i][j].Xc.y/SolnBlk.Wall[i][j].Xwall.y) +
// 	       ihat*(ONE - SolnBlk.Grid.Cell[i][j].Xc.y/SolnBlk.Wall[i][j].Xwall.y);
// 	SolnBlk.W[i][j].v.x = v.x*that.x + v.y*that.y;
// 	SolnBlk.W[i][j].v.y = v.x*that.y - v.y*that.x;
//   	//SolnBlk.W[i][j].p = IP.Wo.p - (2.0/3.0)*SolnBlk.W[i][j].rho*SolnBlk.W[i][j].k;
//   	SolnBlk.W[i][j].p += (2.0/3.0)*SolnBlk.W[i][j].rho*SolnBlk.W[i][j].k;
// 	//SolnBlk.W[i][j].omega *= 100.0;
// // 	SolnBlk.W[i][j].k = 0.01*SolnBlk.W[i][j].v.abs();
// // 	SolnBlk.W[i][j].p += (2.0/3.0)*SolnBlk.W[i][j].rho*SolnBlk.W[i][j].k;
// // 	SolnBlk.W[i][j].omega = MILLION;
// // 	//SolnBlk.W[i][j].omega = pow(10.0,1.5)/max(SolnBlk.W[i][j].lepsilon(SolnBlk.W[i][j].nu(),SolnBlk.Wall[i][j].ywall),TOLER));
// 	SolnBlk.U[i][j] = U(SolnBlk.W[i][j]);
//       }
//     }
//   }

  // Turbulence initial conditions applied successfully.
  return 0;

}

/**********************************************************************
 * Routine: Turbulence_BCs                                            *
 *                                                                    *
 * Apply any required turbulent boundary conditions such as standard, *
 * two-layer, or three-layer wall functions for the k-epsilon and     *
 * k-omega turbulence models or the sublayer value of omega for the   *
 * k-omega turbulence model when integrating directly to the wall.    *
 *                                                                    *
 **********************************************************************/
int Turbulence_BCs(NavierStokes2D_Quad_Block &SolnBlk,
		   NavierStokes2D_Input_Parameters &IP) {

  // Exit immediately if not a turbulent flow.
  if (SolnBlk.Flow_Type != FLOWTYPE_TURBULENT_RANS_K_OMEGA) return 0;

  double omega_di, omega_wf;

  // Apply the appropriate turblence boundary condition if required.
  if (IP.i_Turbulence_BCs == TURBULENT_BC_STANDARD_WALL_FUNCTION) {
    // Standard wall function:
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
 	if (SolnBlk.Wall[i][j].yplus <= IP.yplus_outer_layer &&
	    SolnBlk.Wall[i][j].BCwall != BC_BURNING_SURFACE &&
	    SolnBlk.Wall[i][j].BCwall != BC_MASS_INJECTION
	    ) {
// 	if (SolnBlk.Wall[i][j].yplus <= IP.yplus_buffer_layer ||
// 	    (SolnBlk.Wall[i][j].yplus <= IP.yplus_outer_layer &&
// 	     (SolnBlk.Wall[i+1][j].yplus <= IP.yplus_buffer_layer ||
// 	      SolnBlk.Wall[i-1][j].yplus <= IP.yplus_buffer_layer ||
// 	      SolnBlk.Wall[i][j+1].yplus <= IP.yplus_buffer_layer ||
// 	      SolnBlk.Wall[i][j-1].yplus <= IP.yplus_buffer_layer))) {
	  SolnBlk.W[i][j].k = sqr(SolnBlk.Wall[i][j].utau)/sqrt(SolnBlk.W[i][j].beta_k_o);
	  if (SolnBlk.Wall[i][j].utau == ZERO) {
	    SolnBlk.W[i][j].omega = SIX*SolnBlk.W[i][j].nu()/(SolnBlk.W[i][j].beta_omega_o*
							      sqr(SolnBlk.Wall[i][j].ywall));
	  } else {
	    SolnBlk.W[i][j].omega = SolnBlk.Wall[i][j].utau/(sqrt(SolnBlk.W[i][j].beta_k_o)*
							     SolnBlk.W[i][j].von_karman*SolnBlk.Wall[i][j].ywall);
	  }
	  SolnBlk.U[i][j].dk = SolnBlk.W[i][j].dk();
	  SolnBlk.U[i][j].domega = SolnBlk.W[i][j].domega();
	}
      }
    }

  } else if (IP.i_Turbulence_BCs == TURBULENT_BC_DIRECT_INTEGRATION) {
    // Sublayer for direct integration.
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
	if (SolnBlk.Wall[i][j].yplus <= IP.yplus_sublayer &&
	    SolnBlk.Wall[i][j].BCwall != BC_BURNING_SURFACE &&
	    SolnBlk.Wall[i][j].BCwall != BC_MASS_INJECTION
	    ) {
	  //if (SolnBlk.Wall[i][j].BCwall != BC_BURNING_SURFACE) {
	  SolnBlk.W[i][j].omega = SIX*SolnBlk.W[i][j].nu()/(SolnBlk.W[i][j].beta_omega_o*
							    sqr(SolnBlk.Wall[i][j].ywall));
	  //} else {
	  //nu = SolnBlk.W[i][j].nu();
	  //omega_w = (sqr(SolnBlk.Wall[i][j].utau)/nu)*25.0/(SolnBlk.Wall[i][j].vwplus*(ONE + FIVE*SolnBlk.Wall[i][j].vwplus));
	  //omega_di = omega_w/sqr(ONE + sqrt(SolnBlk.W[i][j].beta_omega_o*omega_w/(SIX*nu))*SolnBlk.Wall[i][j].ywall);
	  //}
	  SolnBlk.U[i][j].domega = SolnBlk.W[i][j].domega();
	}
	if (SolnBlk.Wall[i][j].yplus <= 6.0 &&
	    SolnBlk.Wall[i][j].BCwall != BC_BURNING_SURFACE &&
	    SolnBlk.Wall[i][j].BCwall != BC_MASS_INJECTION){
	  //double Vpr = 0.2 + 0.0125*sqr(SolnBlk.Wall[i][j].yplus);
	  // double Vpr = 0.01 + 0.0028*sqr(SolnBlk.Wall[i][j].yplus);
	  double Vpr = 0.2 + 0.00255*sqr(SolnBlk.Wall[i][j].yplus);
	  double CmuFmu = SolnBlk.W[i][j].muT()*SolnBlk.W[i][j].epsilon()/SolnBlk.W[i][j].rho/max(sqr(SolnBlk.W[i][j].k),sqr(TOLER));
	  double Coeff = sqr(Vpr*0.14*SolnBlk.W[i][j].f_lambda(SolnBlk.Wall[i][j].ywall)/max(CmuFmu,sqr(TOLER)));
	  SolnBlk.W[i][j].ee = Coeff*SolnBlk.W[i][j].epsilon()* SolnBlk.W[i][j].ke/ max(SolnBlk.W[i][j].k,TOLER);
	  //  SolnBlk.W[i][j].ee = SolnBlk.W[i][j].Alpha()*(ONE/FOUR/SolnBlk.W[i][j].ke*sqr(SolnBlk.dWdy[i][j].ke));
	  
	  SolnBlk.U[i][j].dee = SolnBlk.W[i][j].dee();
	}
      }
    }

  } else if (IP.i_Turbulence_BCs == TURBULENT_BC_AUTOMATIC_WALL_TREATMENT) {

    // Apply the automatic treatment for the wall boundary condition.
    // The automated algorithm is only used in the first cell off the
    // wall.  Direct integration is used on the other cells only if they
    // are within the laminar sublayer.  For details see Gao and Groth
    // AIAA-2006-1448.

    int Turbulent_BCtype = TURBULENT_BC_STANDARD_WALL_FUNCTION;
    int applied_flag;

    // NORTH Boundary:
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
//       if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
// 	  SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX) {
	if (SolnBlk.Wall[i][SolnBlk.JCu].yplus <= IP.yplus_sublayer &&
	    SolnBlk.Wall[i][SolnBlk.JCu-IP.n_cells_sublayer_support].yplus <= IP.yplus_sublayer) {
	  Turbulent_BCtype = TURBULENT_BC_DIRECT_INTEGRATION;
	} else if (SolnBlk.Wall[i][SolnBlk.JCu].yplus <= IP.yplus_buffer_layer) {
	  Turbulent_BCtype = TURBULENT_BC_AUTOMATIC_WALL_TREATMENT;
// 	} else if (SolnBlk.Wall[i][SolnBlk.JCu].yplus <= IP.yplus_outer_layer) {
	} else {
	  Turbulent_BCtype = TURBULENT_BC_STANDARD_WALL_FUNCTION;
// 	} else {
// 	  cout << endl << " NORTH:" << SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc << SolnBlk.W[i][SolnBlk.JCu] << " " << SolnBlk.Wall[i][SolnBlk.JCu].yplus; cout.flush();
// 	  return 1101;
	}
	for (int j = SolnBlk.JCu; j >= SolnBlk.JCl; j--) {
	  applied_flag = Apply_Turbulence_BCs(SolnBlk,IP,
					      Turbulent_BCtype,i,j);
	  if (!applied_flag) break;
	}
//       }
    }

    // SOUTH Boundary:
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	  SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX) {
	if (SolnBlk.Wall[i][SolnBlk.JCl].yplus <= IP.yplus_sublayer &&
	    SolnBlk.Wall[i][SolnBlk.JCl+IP.n_cells_sublayer_support].yplus <= IP.yplus_sublayer) {
	  Turbulent_BCtype = TURBULENT_BC_DIRECT_INTEGRATION;
	} else if (SolnBlk.Wall[i][SolnBlk.JCl].yplus <= IP.yplus_buffer_layer) {
	  Turbulent_BCtype = TURBULENT_BC_AUTOMATIC_WALL_TREATMENT;
	} else if (SolnBlk.Wall[i][SolnBlk.JCl].yplus <= IP.yplus_outer_layer) {
	  Turbulent_BCtype = TURBULENT_BC_STANDARD_WALL_FUNCTION;
	} else {
	  cout << endl << " SOUTH:" << SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc << " " << SolnBlk.Wall[i][SolnBlk.JCl].yplus; cout.flush();
	  return 1102;
	}
	for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
	  applied_flag = Apply_Turbulence_BCs(SolnBlk,IP,
					      Turbulent_BCtype,i,j);
	  if (!applied_flag) break;
	}
      }
    }

    // EAST Boundary:
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	  SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX) {
	if (SolnBlk.Wall[SolnBlk.ICu][j].yplus <= IP.yplus_sublayer &&
	    SolnBlk.Wall[SolnBlk.ICu-IP.n_cells_sublayer_support][j].yplus <= IP.yplus_sublayer) {
	  Turbulent_BCtype = TURBULENT_BC_DIRECT_INTEGRATION;
	} else if (SolnBlk.Wall[SolnBlk.ICu][j].yplus <= IP.yplus_buffer_layer) {
	  Turbulent_BCtype = TURBULENT_BC_AUTOMATIC_WALL_TREATMENT;
	} else if (SolnBlk.Wall[SolnBlk.ICu][j].yplus <= IP.yplus_outer_layer) {
	  Turbulent_BCtype = TURBULENT_BC_STANDARD_WALL_FUNCTION;
	} else {
	  cout << endl << " EAST:" << SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc << " " << SolnBlk.Wall[SolnBlk.ICu][j].yplus; cout.flush();
	  return 1103;
	}
	for (int i = SolnBlk.ICu; i >= SolnBlk.ICl; i--) {
	  applied_flag = Apply_Turbulence_BCs(SolnBlk,IP,
					      Turbulent_BCtype,i,j);
	  if (!applied_flag) break;
	}
      }
    }

    // WEST Boundary:
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	  SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX) {
	if (SolnBlk.Wall[SolnBlk.ICl][j].yplus <= IP.yplus_sublayer &&
	    SolnBlk.Wall[SolnBlk.ICl+IP.n_cells_sublayer_support][j].yplus <= IP.yplus_sublayer) {
	  Turbulent_BCtype = TURBULENT_BC_DIRECT_INTEGRATION;
	} else if (SolnBlk.Wall[SolnBlk.ICl][j].yplus <= IP.yplus_buffer_layer) {
	  Turbulent_BCtype = TURBULENT_BC_AUTOMATIC_WALL_TREATMENT;
	} else if (SolnBlk.Wall[SolnBlk.ICl][j].yplus <= IP.yplus_outer_layer) {
	  Turbulent_BCtype = TURBULENT_BC_STANDARD_WALL_FUNCTION;
	} else {
	  cout << endl << " WEST:" << SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc << " " << SolnBlk.Wall[SolnBlk.ICl][j].yplus; cout.flush();
	  return 1104;
	}
	for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
	  applied_flag = Apply_Turbulence_BCs(SolnBlk,IP,
					      Turbulent_BCtype,i,j);
	  if (!applied_flag) break;
	}
      }
    }

  } else {

    // Turbulence boundary condition not specified properly.
    cout << endl << " " << 1105;
    return 1105;

  }

  // The turbulence boundary conditions have been successfully applied.
  return 0;

}

/**********************************************************************
 **********************************************************************/
int Apply_Turbulence_BCs(NavierStokes2D_Quad_Block &SolnBlk,
			 NavierStokes2D_Input_Parameters &IP,
			 const int &Turbulent_BCtype,
			 const int &i, const int &j) {

  if (SolnBlk.Wall[i][j].BCwall == BC_BURNING_SURFACE) return 0;
  //if (SolnBlk.Wall[i][j].BCwall == BC_MASS_INJECTION) return 0;

  double omega_di, omega_wf;

  // The cell is within the laminar sublayer: use direct integration.
  if (Turbulent_BCtype == TURBULENT_BC_DIRECT_INTEGRATION &&
      SolnBlk.Wall[i][j].yplus <= IP.yplus_sublayer) {
    if (SolnBlk.Wall[i][j].BCwall != BC_MASS_INJECTION) {
      SolnBlk.W[i][j].omega = SIX*SolnBlk.W[i][j].nu()/(SolnBlk.W[i][j].beta_omega_o*sqr(SolnBlk.Wall[i][j].ywall));
    } else {
      SolnBlk.W[i][j].omega = sqrt(SolnBlk.W[i][j].k)/SolnBlk.W[i][j].lw;
    }
    SolnBlk.U[i][j].domega = SolnBlk.W[i][j].domega();

    SolnBlk.W[i][j].ee = SolnBlk.W[i][j].Alpha()*sqr(ONE/TWO/sqrt(SolnBlk.W[i][j].ke)*(SolnBlk.dWdx[i][j].ke*SolnBlk.Wall[i][j].nwall.x+SolnBlk.dWdy[i][j].ke*SolnBlk.Wall[i][j].nwall.y));
    SolnBlk.U[i][j].dee = SolnBlk.W[i][j].dee();

  // The cell is within the buffer layer: apply the blending function.
  } else if (Turbulent_BCtype == TURBULENT_BC_AUTOMATIC_WALL_TREATMENT &&
	     SolnBlk.Wall[i][j].yplus <= IP.yplus_buffer_layer) {
    // Determine the specific dissipation as specifed by direct integration.
    if (SolnBlk.Wall[i][j].BCwall != BC_MASS_INJECTION) {
      omega_di = SIX*SolnBlk.W[i][j].nu()/(SolnBlk.W[i][j].beta_omega_o*sqr(SolnBlk.Wall[i][j].ywall));
    } else {
      SolnBlk.W[i][j].omega = sqrt(SolnBlk.W[i][j].k)/SolnBlk.W[i][j].lw;
    }
    // Determine the specific dissipation as specifed by the
    // standard wall-functions.
    omega_wf = SolnBlk.Wall[i][j].utau/(sqrt(SolnBlk.W[i][j].beta_k_o)*
					SolnBlk.W[i][j].von_karman*SolnBlk.Wall[i][j].ywall);
    // Determine the specific dissipation by blending the two results.
    SolnBlk.W[i][j].omega = sqrt(sqr(omega_di) + sqr(omega_wf));
    // The turbulent kinetic energy behaves quadratically in the buffer layer.
    SolnBlk.W[i][j].k = sqr(SolnBlk.Wall[i][j].utau*SolnBlk.Wall[i][j].yplus/IP.yplus_buffer_layer)/sqrt(SolnBlk.W[i][j].beta_k_o);
    SolnBlk.U[i][j].dk = SolnBlk.W[i][j].dk();
    SolnBlk.U[i][j].domega = SolnBlk.W[i][j].domega();

  // The cell is within the outer layer: use standard wall-functions.
  } else if (Turbulent_BCtype == TURBULENT_BC_STANDARD_WALL_FUNCTION &&
	     SolnBlk.Wall[i][SolnBlk.JCu].yplus <= IP.yplus_outer_layer) {
    if (SolnBlk.Wall[i][j].BCwall != BC_MASS_INJECTION) {
      SolnBlk.W[i][j].k = sqr(SolnBlk.Wall[i][j].utau)/sqrt(SolnBlk.W[i][j].beta_k_o);
      if (SolnBlk.Wall[i][j].utau == ZERO) {
	SolnBlk.W[i][j].omega = SIX*SolnBlk.W[i][j].nu()/(SolnBlk.W[i][j].beta_omega_o*sqr(SolnBlk.Wall[i][j].ywall));
      } else {
	SolnBlk.W[i][j].omega = SolnBlk.Wall[i][j].utau/(sqrt(SolnBlk.W[i][j].beta_k_o)*
							 SolnBlk.W[i][j].von_karman*SolnBlk.Wall[i][j].ywall);
      }
    } else {
      SolnBlk.W[i][j].k = sqr(SolnBlk.W[i][j].sigmav)*sqr(3.1);//SolnBlk.W[i][j].v.abs());
      SolnBlk.W[i][j].omega = sqrt(SolnBlk.W[i][j].k)/SolnBlk.W[i][j].lw;
    }
    SolnBlk.U[i][j].dk = SolnBlk.W[i][j].dk();
    SolnBlk.U[i][j].domega = SolnBlk.W[i][j].domega();

  // No boundary condition to use.  Report that none were applied.
  } else {
    return 0;

  }

  // Appropriate turbulence boundary condition applied.
  return 1;

}

/**********************************************************************
 * Routine: Turbulence_Zero_Residual                                  *
 *                                                                    *
 **********************************************************************/
int Turbulence_Zero_Residual(NavierStokes2D_Quad_Block &SolnBlk,
                             const int i_stage,
			     NavierStokes2D_Input_Parameters &IP) {

  // Exit immediately if not a turbulent flow.
  if (!SolnBlk.Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) return 0;

  int error_flag, k_residual;
  double omega;

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
      if (i_stage == 4) k_residual = 0;
      else k_residual = i_stage - 1;
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

  if (IP.i_Turbulence_BCs == TURBULENT_BC_STANDARD_WALL_FUNCTION) {
    // Zero the residuals of the turbulence quantities if the current
    // cell is within the log layer (yplus is less than the specified
    // value of yplus_outer_layer).
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
 	if (SolnBlk.Wall[i][j].yplus <= IP.yplus_outer_layer &&
	    SolnBlk.Wall[i][j].BCwall != BC_BURNING_SURFACE// &&
	    //SolnBlk.Wall[i][j].BCwall != BC_MASS_INJECTION
	    ) {
// 	  if (SolnBlk.Wall[i][j].yplus <= IP.yplus_buffer_layer ||
// 	      (SolnBlk.Wall[i][j].yplus <= IP.yplus_outer_layer &&
// 	       (SolnBlk.Wall[i+1][j].yplus <= IP.yplus_buffer_layer ||
// 		SolnBlk.Wall[i-1][j].yplus <= IP.yplus_buffer_layer ||
// 		SolnBlk.Wall[i][j+1].yplus <= IP.yplus_buffer_layer ||
// 		SolnBlk.Wall[i][j-1].yplus <= IP.yplus_buffer_layer))) {
 	  SolnBlk.dUdt[i][j][k_residual].dk = ZERO;
 	  SolnBlk.dUdt[i][j][k_residual].domega = ZERO;
 	}
      }
    }
  } else if (IP.i_Turbulence_BCs == TURBULENT_BC_DIRECT_INTEGRATION) {
    // Zero the residuals of the specific dissipation if the current
    // cell is within the viscous sublayer (yplus is less than the
    // specified value of yplus_sublayer).
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
 	if (SolnBlk.Wall[i][j].yplus <= IP.yplus_sublayer &&
	    SolnBlk.Wall[i][j].BCwall != BC_BURNING_SURFACE// &&
	    //SolnBlk.Wall[i][j].BCwall != BC_MASS_INJECTION
	    ) {
	  SolnBlk.dUdt[i][j][k_residual].domega = ZERO;
 	}
	if (SolnBlk.Wall[i][j].yplus <= 6.0&&
	    SolnBlk.Wall[i][j].BCwall != BC_BURNING_SURFACE &&
	    SolnBlk.Wall[i][j].BCwall != BC_MASS_INJECTION){
	  SolnBlk.dUdt[i][j][k_residual].dee = ZERO;
	}
      }
    }
  } else if (IP.i_Turbulence_BCs == TURBULENT_BC_AUTOMATIC_WALL_TREATMENT) {
    // Zero the residuals of the turbulence quantities depending on the
    // automatic wall treatment.
    int Turbulent_BCtype = TURBULENT_BC_STANDARD_WALL_FUNCTION;
    int zeroed_flag;

    // NORTH Boundary:
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      if (SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	  SolnBlk.Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX) {
	if (SolnBlk.Wall[i][SolnBlk.JCu].yplus <= IP.yplus_sublayer &&
	    SolnBlk.Wall[i][SolnBlk.JCu-IP.n_cells_sublayer_support].yplus <= IP.yplus_sublayer) {
	  Turbulent_BCtype = TURBULENT_BC_DIRECT_INTEGRATION;
	} else if (SolnBlk.Wall[i][SolnBlk.JCu].yplus <= IP.yplus_buffer_layer) {
	  Turbulent_BCtype = TURBULENT_BC_AUTOMATIC_WALL_TREATMENT;
	} else if (SolnBlk.Wall[i][SolnBlk.JCu].yplus <= IP.yplus_outer_layer) {
	  Turbulent_BCtype = TURBULENT_BC_STANDARD_WALL_FUNCTION;
	} else {
	  cout << endl << " NORTH:" << SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc << " " << SolnBlk.Wall[i][SolnBlk.JCu].yplus; cout.flush();
	  return 2101;
	}
	for (int j = SolnBlk.JCu; j >= SolnBlk.JCl; j--) {
	  zeroed_flag = Zero_Turbulence_Residuals(SolnBlk,IP,
						  Turbulent_BCtype,i,j,k_residual);
	  if (!zeroed_flag) break;
	}
      }
    }

    // SOUTH Boundary:
    for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
      if (SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	  SolnBlk.Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX) {
	if (SolnBlk.Wall[i][SolnBlk.JCl].yplus <= IP.yplus_sublayer &&
	    SolnBlk.Wall[i][SolnBlk.JCl+IP.n_cells_sublayer_support].yplus <= IP.yplus_sublayer) {
	  Turbulent_BCtype = TURBULENT_BC_DIRECT_INTEGRATION;
	} else if (SolnBlk.Wall[i][SolnBlk.JCl].yplus <= IP.yplus_buffer_layer) {
	  Turbulent_BCtype = TURBULENT_BC_AUTOMATIC_WALL_TREATMENT;
	} else if (SolnBlk.Wall[i][SolnBlk.JCl].yplus <= IP.yplus_outer_layer) {
	  Turbulent_BCtype = TURBULENT_BC_STANDARD_WALL_FUNCTION;
	} else {
	  cout << endl << " SOUTH:" << SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc << " " << SolnBlk.Wall[i][SolnBlk.JCl].yplus; cout.flush();
	  return 2102;
	}
	for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
	  zeroed_flag = Zero_Turbulence_Residuals(SolnBlk,IP,
						  Turbulent_BCtype,i,j,k_residual);
	  if (!zeroed_flag) break;
	}
      }
    }

    // EAST Boundary:
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      if (SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	  SolnBlk.Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX) {
	if (SolnBlk.Wall[SolnBlk.ICu][j].yplus <= IP.yplus_sublayer &&
	    SolnBlk.Wall[SolnBlk.ICu-IP.n_cells_sublayer_support][j].yplus <= IP.yplus_sublayer) {
	  Turbulent_BCtype = TURBULENT_BC_DIRECT_INTEGRATION;
	} else if (SolnBlk.Wall[SolnBlk.ICu][j].yplus <= IP.yplus_buffer_layer) {
	  Turbulent_BCtype = TURBULENT_BC_AUTOMATIC_WALL_TREATMENT;
	} else if (SolnBlk.Wall[SolnBlk.ICu][j].yplus <= IP.yplus_outer_layer) {
	  Turbulent_BCtype = TURBULENT_BC_STANDARD_WALL_FUNCTION;
	} else {
	  cout << endl << " EAST:" << SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc << " " << SolnBlk.Wall[SolnBlk.ICu][j].yplus; cout.flush();
	  return 2103;
	}
	for (int i = SolnBlk.ICu; i >= SolnBlk.ICl; i--) {
	  zeroed_flag = Zero_Turbulence_Residuals(SolnBlk,IP,
						  Turbulent_BCtype,i,j,k_residual);
	  if (!zeroed_flag) break;
	}
      }
    }

    // WEST Boundary:
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      if (SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	  SolnBlk.Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX) {
	if (SolnBlk.Wall[SolnBlk.ICl][j].yplus <= IP.yplus_sublayer &&
	    SolnBlk.Wall[SolnBlk.ICl+IP.n_cells_sublayer_support][j].yplus <= IP.yplus_sublayer) {
	  Turbulent_BCtype = TURBULENT_BC_DIRECT_INTEGRATION;
	} else if (SolnBlk.Wall[SolnBlk.ICl][j].yplus <= IP.yplus_buffer_layer) {
	  Turbulent_BCtype = TURBULENT_BC_AUTOMATIC_WALL_TREATMENT;
	} else if (SolnBlk.Wall[SolnBlk.ICl][j].yplus <= IP.yplus_outer_layer) {
	  Turbulent_BCtype = TURBULENT_BC_STANDARD_WALL_FUNCTION;
	} else {
	  cout << endl << " WEST:" << SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc << " " << SolnBlk.Wall[SolnBlk.ICl][j].yplus; cout.flush();
	  return 2104;
	}
	for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
	  zeroed_flag = Zero_Turbulence_Residuals(SolnBlk,IP,
						  Turbulent_BCtype,i,j,k_residual);
	  if (!zeroed_flag) break;
	}
      }
    }

  } else {
    return 2105;

  }

  // The residuals for the turbulence parameters have been zeroed
  // according to the dimensionless wall distance and the local 
  // turbulence boundary condition type.
  return 0;

}

/**********************************************************************
 **********************************************************************/
int Zero_Turbulence_Residuals(NavierStokes2D_Quad_Block &SolnBlk,
			      NavierStokes2D_Input_Parameters &IP,
			      const int &Turbulent_BCtype,
			      const int &i, const int &j,
			      const int &k_residual) {

  if (SolnBlk.Wall[i][j].BCwall == BC_BURNING_SURFACE) return 0;
  //if (SolnBlk.Wall[i][j].BCwall == BC_MASS_INJECTION) return 0;

  // The cell is within the laminar sublayer: use direct integration.
  if (Turbulent_BCtype == TURBULENT_BC_DIRECT_INTEGRATION &&
      SolnBlk.Wall[i][j].yplus <= IP.yplus_sublayer) {
    SolnBlk.dUdt[i][j][k_residual].domega = ZERO;
    SolnBlk.dUdt[i][j][k_residual].dee = ZERO;

  // The cell is within the buffer layer: apply the blending function.
  } else if (Turbulent_BCtype == TURBULENT_BC_AUTOMATIC_WALL_TREATMENT &&
	     SolnBlk.Wall[i][j].yplus <= IP.yplus_buffer_layer) {
    SolnBlk.dUdt[i][j][k_residual].dk = ZERO;
    SolnBlk.dUdt[i][j][k_residual].domega = ZERO;

  // The cell is within the outer layer: use standard wall-functions.
  } else if (Turbulent_BCtype == TURBULENT_BC_STANDARD_WALL_FUNCTION &&
	     SolnBlk.Wall[i][SolnBlk.JCu].yplus <= IP.yplus_outer_layer) {
    SolnBlk.dUdt[i][j][k_residual].dk = ZERO;
    SolnBlk.dUdt[i][j][k_residual].domega = ZERO;

  // No boundary condition to use.  Report that no residuals were zeroed.
  } else {
    return 0;

  }

  // Appropriate turbulence residuals were zeroed.
  return 1;

}
