/**********************************************************************
 * HighTemp2DQuadTurbulenceSingleBlock.cc: Single-block versions      *
 *                                             of turbulence          *
 *                                             subroutines for 2D     *
 *                                             High-Temp multi-       *
 *                                             block quadrilateral    *
 *                                             mesh solution classes. *
 **********************************************************************/

#include "HighTemp2DQuad.h"

/**********************************************************************
 * HighTemp2D_Quad_Block -- Turbulence Single Block External          *
 *                          Subroutines.                              *
 **********************************************************************/

/**********************************************************************
 * Routine: Turbulent_ICs                                             *
 *                                                                    *
 * Assigns initial turbulence conditions and data to the solution     *
 * variables of the specified quadrilateral solution block.           *
 *                                                                    *
 **********************************************************************/
int Turbulent_ICs(HighTemp2D_Quad_Block &SolnBlk,
		  HighTemp2D_Input_Parameters &IP,
		  HighTemp2D_pState *Wo) {

  // Exit immediately if not a turbulent flow.
  if (SolnBlk.Flow_Type != FLOWTYPE_TURBULENT_RANS_K_OMEGA) return 0;

  // Turbulence initial conditions applied successfully.
  return 0;

}

/**********************************************************************
 * Routine: Turbulent_BCs                                             *
 *                                                                    *
 * Apply any required turbulent boundary conditions such as standard, *
 * two-layer, or three-layer wall functions for the k-epsilon and     *
 * k-omega turbulence models or the sublayer value of omega for the   *
 * k-omega turbulence model when integrating directly to the wall.    *
 *                                                                    *
 **********************************************************************/
int Turbulent_BCs(HighTemp2D_Quad_Block &SolnBlk,
		  const HighTemp2D_Input_Parameters &IP) {

  // Exit immediately if not a turbulent flow.
  if (SolnBlk.Flow_Type != FLOWTYPE_TURBULENT_RANS_K_OMEGA) return 0;

  // Apply the appropriate turblence boundary condition if required.
  if (IP.i_Turbulent_BCs == TURBULENT_BC_STANDARD_WALL_FUNCTION) {
    // Standard wall function:
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
	  if (SolnBlk.Wall[i][j].yplus <= IP.yplus_outer_layer &&
	    SolnBlk.Wall[i][j].BCwall != BC_BURNING_SURFACE) {
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

  } else if (IP.i_Turbulent_BCs == TURBULENT_BC_DIRECT_INTEGRATION) {
    // Sublayer for direct integration.
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
	if (SolnBlk.Wall[i][j].yplus <= IP.yplus_sublayer) {
// 	  if (SolnBlk.Wall[i][j].BCwall != BC_BURNING_SURFACE) {
	    SolnBlk.W[i][j].omega = SIX*SolnBlk.W[i][j].nu()/(SolnBlk.W[i][j].beta_omega_o*
							      sqr(SolnBlk.Wall[i][j].ywall));
// 	  } else {
// 	    SolnBlk.W[i][j].omega = (sqr(SolnBlk.Wall[i][j].utau)/SolnBlk.W[i][j].nu())*25.0/
//                                  (SolnBlk.Wall[i][j].vwplus*(ONE + FIVE*SolnBlk.Wall[i][j].vwplus));
// 	  }
	  SolnBlk.U[i][j].domega = SolnBlk.W[i][j].domega();
	}
      }
    }

  } else if (IP.i_Turbulent_BCs == TURBULENT_BC_AUTOMATIC_WALL_TREATMENT) {

    double omega_di, omega_wf;

    // Apply the automatic treatment for the wall boundary condition.
    // The automated algorithm is only used in the first cell off the
    // wall.  Direct integration is used on the other cells only if they
    // are within the laminar sublayer.  For details see Gao and Groth
    // AIAA-2006-1448.

     int Turbulent_BCtype = TURBULENT_BC_STANDARD_WALL_FUNCTION;
     int applied_flag;

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
	  cout << endl << " NORTH:" << SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc 
               << SolnBlk.W[i][SolnBlk.JCu] << " " << SolnBlk.Wall[i][SolnBlk.JCu].yplus; cout.flush();
	  return 1101;
	}
	for (int j = SolnBlk.JCu; j >= SolnBlk.JCl; j--) {
	  applied_flag = Apply_Turbulence_BCs(SolnBlk,IP,
					      Turbulent_BCtype,i,j);
	  if (!applied_flag) break;
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
	  cout << endl << " SOUTH:" << SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc << " " 
               << SolnBlk.Wall[i][SolnBlk.JCl].yplus; cout.flush();
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
	  cout << endl << " EAST:" << SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc << " " 
               << SolnBlk.Wall[SolnBlk.ICu][j].yplus; cout.flush();
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
	  cout << endl << " WEST:" << SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc << " " 
               << SolnBlk.Wall[SolnBlk.ICl][j].yplus; cout.flush();
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
    //return 1102;
    cout << endl << " " << 1105;
    return 1105;

  }

  // The turbulence boundary conditions have been successfully applied.
  return 0;

}

/**********************************************************************
 **********************************************************************/
int Apply_Turbulence_BCs(HighTemp2D_Quad_Block &SolnBlk,
			 const HighTemp2D_Input_Parameters &IP,
			 const int &Turbulent_BCtype,
			 const int &i, const int &j) {

  if (SolnBlk.Wall[i][j].BCwall == BC_BURNING_SURFACE) return 0;
  //if (SolnBlk.Wall[i][j].BCwall == BC_MASS_INJECTION) return 0;

  double omega_di, omega_wf;

  // The cell is within the laminar sublayer: use direct integration.
  if (Turbulent_BCtype == TURBULENT_BC_DIRECT_INTEGRATION &&
      SolnBlk.Wall[i][j].yplus <= IP.yplus_sublayer) {
    SolnBlk.W[i][j].omega = SIX*SolnBlk.W[i][j].nu()/(SolnBlk.W[i][j].beta_omega_o*sqr(SolnBlk.Wall[i][j].ywall));
    SolnBlk.U[i][j].domega = SolnBlk.W[i][j].domega();

  // The cell is within the buffer layer: apply the blending function.
  } else if (Turbulent_BCtype == TURBULENT_BC_AUTOMATIC_WALL_TREATMENT &&
	     SolnBlk.Wall[i][j].yplus <= IP.yplus_buffer_layer) {
    // Determine the specific dissipation as specifed by direct integration.
    omega_di = SIX*SolnBlk.W[i][j].nu()/(SolnBlk.W[i][j].beta_omega_o*sqr(SolnBlk.Wall[i][j].ywall));
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
    SolnBlk.W[i][j].k = sqr(SolnBlk.Wall[i][j].utau)/sqrt(SolnBlk.W[i][j].beta_k_o);
    if (SolnBlk.Wall[i][j].utau == ZERO) {
      SolnBlk.W[i][j].omega = SIX*SolnBlk.W[i][j].nu()/(SolnBlk.W[i][j].beta_omega_o*sqr(SolnBlk.Wall[i][j].ywall));
    } else {
      SolnBlk.W[i][j].omega = SolnBlk.Wall[i][j].utau/(sqrt(SolnBlk.W[i][j].beta_k_o)*
						       SolnBlk.W[i][j].von_karman*SolnBlk.Wall[i][j].ywall);
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
int Turbulence_Zero_Residual(HighTemp2D_Quad_Block &SolnBlk,
                             const int i_stage,
			     HighTemp2D_Input_Parameters &IP) {

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

  if (IP.i_Turbulent_BCs == TURBULENT_BC_STANDARD_WALL_FUNCTION) {
    // Zero the residuals of the turbulence quantities if the current
    // cell is within the log layer (yplus is less than the specified
    // value of yplus_outer_layer).
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
 	if (SolnBlk.Wall[i][j].yplus <= IP.yplus_outer_layer &&
	    SolnBlk.Wall[i][j].BCwall != BC_BURNING_SURFACE) {
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
  } else if (IP.i_Turbulent_BCs == TURBULENT_BC_DIRECT_INTEGRATION) {
    // Zero the residuals of the specific dissipation if the current
    // cell is within the viscous sublayer (yplus is less than the
    // specified value of yplus_sublayer).
    for (int j = SolnBlk.JCl; j <= SolnBlk.JCu; j++) {
      for (int i = SolnBlk.ICl; i <= SolnBlk.ICu; i++) {
 	if (SolnBlk.Wall[i][j].yplus <= IP.yplus_sublayer &&
	    SolnBlk.Wall[i][j].BCwall != BC_BURNING_SURFACE) {
	  SolnBlk.dUdt[i][j][k_residual].domega = ZERO;
 	}
      }
    }
  } else if (IP.i_Turbulent_BCs == TURBULENT_BC_AUTOMATIC_WALL_TREATMENT) {
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
	  cout << endl << " NORTH:" << SolnBlk.Grid.Cell[i][SolnBlk.JCu].Xc << " " 
               << SolnBlk.Wall[i][SolnBlk.JCu].yplus; cout.flush();
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
	  cout << endl << " SOUTH:" << SolnBlk.Grid.Cell[i][SolnBlk.JCl].Xc << " " 
               << SolnBlk.Wall[i][SolnBlk.JCl].yplus; cout.flush();
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
	  cout << endl << " EAST:" << SolnBlk.Grid.Cell[SolnBlk.ICu][j].Xc << " " 
               << SolnBlk.Wall[SolnBlk.ICu][j].yplus; cout.flush();
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
	  cout << endl << " WEST:" << SolnBlk.Grid.Cell[SolnBlk.ICl][j].Xc << " " 
               << SolnBlk.Wall[SolnBlk.ICl][j].yplus; cout.flush();
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
int Zero_Turbulence_Residuals(HighTemp2D_Quad_Block &SolnBlk,
			      HighTemp2D_Input_Parameters &IP,
			      const int &Turbulent_BCtype,
			      const int &i, const int &j,
			      const int &k_residual) {

  if (SolnBlk.Wall[i][j].BCwall == BC_BURNING_SURFACE) return 0;
  //  if (SolnBlk.Wall[i][j].BCwall == BC_MASS_INJECTION) return 0;

  // The cell is within the laminar sublayer: use direct integration.
  if (Turbulent_BCtype == TURBULENT_BC_DIRECT_INTEGRATION &&
      SolnBlk.Wall[i][j].yplus <= IP.yplus_sublayer) {
    SolnBlk.dUdt[i][j][k_residual].domega = ZERO;

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
