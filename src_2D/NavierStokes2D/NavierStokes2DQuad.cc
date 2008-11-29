/*! \file NavierStokes2DQuad.cc
  @brief Subroutines for 2D Navier-Stokes quadrilateral solution block class. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "NavierStokes2DQuad.h"  // Include 2D NavierStokes quadrilateral solution block class.


/**********************************************************************
 * NavierStokes2D_Quad_Block -- Create storage for the static         *
 *                              variables.                            *
 **********************************************************************/
int NavierStokes2D_Quad_Block::residual_variable = 1;
#ifdef _NS_PARALLEL_DEBUG_
ofstream NavierStokes2D_Quad_Block::dout;
#endif
int NavierStokes2D_Quad_Block::Number_of_Residual_Norms = 4;
// Initialize RefW
NavierStokes2D_pState NavierStokes2D_Quad_Block::RefW(1.0);
// Initialize ExactSoln
NavierStokes2D_ExactSolutions *NavierStokes2D_Quad_Block::ExactSoln = NULL;

/***********************************************************************
 * NavierStokes2D_Quad_Block -- Single Block Member Functions.         *
 **********************************************************************/

/*****************************************************//**
 * Copy the solution information of quadrilateral solution 
 * block SolnBlk to the current solution block.
 ********************************************************/
NavierStokes2D_Quad_Block & NavierStokes2D_Quad_Block::operator =(const NavierStokes2D_Quad_Block &Soln){

  int i, j, k;

  // Handle self-assignment:
  if (this == & Soln) return *this;

  // check if solution block Soln has memory allocated
  if (Soln.U != NULL){
    /* Allocate (re-allocate) memory for the solution
       of the quadrilateral solution block as necessary. */
    allocate(Soln.NCi-2*Soln.Nghost,
	     Soln.NCj-2*Soln.Nghost,
	     Soln.Nghost);
  } else {
    deallocate();
  }

  // Copy the viscous flow indicator.
  Flow_Type = Soln.Flow_Type;

  // Copy the axisymmetric/planar flow indicator.
  Axisymmetric = Soln.Axisymmetric;

  // Copy the compressibility effect correction indicator.
  Compressibility_Effect = Soln.Compressibility_Effect;

  // Copy the transition model correction indicator.
  Transition_Model = Soln.Transition_Model;

  // Copy the variable Prandtl number indicator.
  Variable_Prandtl = Soln.Variable_Prandtl;

  // Copy the wall velocity.
  Vwall = Soln.Vwall;

  // Copy the wall temperature.
  Twall = Soln.Twall;

  /* Copy the grid. */
  Grid = Soln.Grid;

  /* Copy the solution information from Soln. */
  if (Soln.U != NULL) {
    for ( j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
      for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
	U[i][j] = Soln.U[i][j];
	W[i][j] = Soln.W[i][j];
	for (k = 0; k < NUMBER_OF_RESIDUAL_VECTORS_NAVIERSTOKES2D; ++k){
	  dUdt[i][j][k] = Soln.dUdt[i][j][k];
	} /* endfor */
	dWdx[i][j] = Soln.dWdx[i][j];
	dWdy[i][j] = Soln.dWdy[i][j];
	for (k = 0; k < 5; ++k){
	  d_dWdx_dW[i][j][k] = Soln.d_dWdx_dW[i][j][k];
	  d_dWdy_dW[i][j][k] = Soln.d_dWdy_dW[i][j][k];
	}
	phi[i][j] = Soln.phi[i][j];
	Uo[i][j] = Soln.Uo[i][j];
	dt[i][j] = Soln.dt[i][j];
      } /* endfor */
    } /* endfor */

    for (j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
      WoW[j] = Soln.WoW[j];
      WoE[j] = Soln.WoE[j];
    }/* endfor */
    
    for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
      WoS[i] = Soln.WoS[i];
      WoN[i] = Soln.WoN[i];
    }/* endfor */
    
  }/* endif */

  // Copy high-order objects
  copy_HighOrder_Objects(Soln);

  // Copy boundary reference states
  Ref_State_BC_North = Soln.Ref_State_BC_North;
  Ref_State_BC_South = Soln.Ref_State_BC_South;
  Ref_State_BC_East = Soln.Ref_State_BC_East;
  Ref_State_BC_West = Soln.Ref_State_BC_West;  

  // Reset accuracy assessment flag
  AssessAccuracy.ResetForNewCalculation();

  return *this;
}

/***************************************************//**
 * Assigns boundary condition reference data based on
 * the boundary type for the specified quadrilateral 
 * solution block based on the input parameters.
 * This routine makes the link between user's specifications
 * and the values that are set as boundary reference states.
 ********************************************************/
void NavierStokes2D_Quad_Block::Set_Boundary_Reference_States_Based_On_Input(const NavierStokes2D_Input_Parameters &IP){

  // Set the reference values for boundary reference states
  Set_Reference_Values_For_Boundary_States(IP.Ref_State_BC_North,
					   IP.Ref_State_BC_South,
					   IP.Ref_State_BC_East,
					   IP.Ref_State_BC_West);
}

/***************************************************//**
 * Assigns reference values for the boundary condition 
 * reference states to what the user specified.
 * These values are used in Set_Boundary_Reference_States()
 * routine.
 ********************************************************/
void NavierStokes2D_Quad_Block::Set_Reference_Values_For_Boundary_States(const NavierStokes2D_pState & Ref_North,
									 const NavierStokes2D_pState & Ref_South,
									 const NavierStokes2D_pState & Ref_East,
									 const NavierStokes2D_pState & Ref_West){
  Ref_State_BC_North = Ref_North;
  Ref_State_BC_South = Ref_South;
  Ref_State_BC_East = Ref_East;
  Ref_State_BC_West = Ref_West;
}

/********************************************************************//**
 * This routine evaluates the residual for the specified solution     
 * block using a 2nd-order limited upwind finite-volume spatial       
 * discretization scheme with either the Godunov, Roe, Rusanov, HLLE, 
 * HLLL, or HLLC flux functions for the gas-phase and with either the 
 * Saurel or Equilibrium flux functions for the particle-phase.  The  
 * residual is stored in dUdt[][][0].                                 
 *                                                                    
 **********************************************************************/
int NavierStokes2D_Quad_Block::dUdt_Residual_Evaluation(NavierStokes2D_Input_Parameters &IP) {

  int error_flag;
  Vector2D dX;
  NavierStokes2D_pState Wl, Wr;
  NavierStokes2D_cState Flux;

  NavierStokes2D_pState Wu, Wd, dWdxl, dWdyl, dWdxr, dWdyr;
  NavierStokes2D_pState dWdx_face, dWdy_face;
  Vector2D Xl, Xr, Xu, Xd;
  int viscous_bc_flag;

  int flux_function = IP.i_Flux_Function;

  NavierStokes2D_cState NavierStokes2D_U_VACUUM; NavierStokes2D_U_VACUUM.Vacuum();
  NavierStokes2D_pState NavierStokes2D_W_STDATM; NavierStokes2D_W_STDATM.Standard_Atmosphere();

  // Perform the linear reconstruction within each cell of the
  // computational grid for this stage.
  switch(IP.i_Reconstruction) {
  case RECONSTRUCTION_GREEN_GAUSS :
    Linear_Reconstruction_GreenGauss(*this,
				     IP.i_Limiter);
    break;
  case RECONSTRUCTION_LINEAR_LEAST_SQUARES :
    Linear_Reconstruction_LeastSquares(*this,
				       IP.i_Limiter);
    break;
  default:
    Linear_Reconstruction_LeastSquares(*this,
				       IP.i_Limiter);
    break;
  };

  // Evaluate the time rate of change of the solution (i.e., the
  // solution residuals) using a second-order limited upwind scheme
  // with a variety of flux functions.

  // Add i-direction (zeta-direction) fluxes.
  for (int j = JCl-1; j <= JCu+1; j++) {

    dUdt[ICl-1][j][0] = NavierStokes2D_U_VACUUM;

    for (int i = ICl-1; i <= ICu; i++) {

      dUdt[i+1][j][0] = NavierStokes2D_U_VACUUM;

      if (j >= JCl && j <= JCu) {

	if (i == ICl-1 && 
	    (Grid.BCtypeW[j] == BC_REFLECTION ||
	     Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
	     Grid.BCtypeW[j] == BC_MASS_INJECTION ||
	     Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
	     Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
	     Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
	     Grid.BCtypeW[j] == BC_RINGLEB_FLOW ||
	     Grid.BCtypeW[j] == BC_CHARACTERISTIC)) {
	  
	  dX = Grid.xfaceW(i+1,j) - Grid.Cell[i+1][j].Xc;
	  Wr = (phi[i+1][j]^dWdx[i+1][j]);


	  Wr = W[i+1][j] + ( (phi[i+1][j]^dWdx[i+1][j])*dX.x +
			     (phi[i+1][j]^dWdy[i+1][j])*dX.y);

	  // WEST face of cell (i+1,j) is a normal boundary.
	  if (Grid.BCtypeW[j] == BC_REFLECTION) {
	    // WEST face of cell (i+1,j) is a REFLECTION boundary.
	    Wl = Reflect(Wr,Grid.nfaceW(i+1,j));
	  } else if (Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
	    // WEST face of cell (i+1,j) is a BURNING_SURFACE boundary.
	    Wl = BurningSurface(Wr,Grid.nfaceW(i+1,j));
	  } else if (Grid.BCtypeW[j] == BC_MASS_INJECTION) {
	    // WEST face of cell (i+1,j) is a MASS_INJECTION boundary.
	    //Wl = MassInjection(Wr,Grid.nfaceW(i+1,j));
	    Wl = MassInjection2(Wr,Grid.xfaceW(i+1,j),Grid.nfaceW(i+1,j),Twall);
	    flux_function = IP.i_Flux_Function;
	    IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV_WRS;
	  } else if (Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX) {
	    // WEST face of cell (i+1,j) is a WALL_VISCOUS_HEATFLUX boundary.
	    Wl = WallViscousHeatFlux(Wr,Grid.nfaceW(i+1,j));
	  } else if (Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    // WEST face of cell (i+1,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	    Wl = WallViscousIsothermal(Wr,Grid.nfaceW(i+1,j),Twall);
	  } else if (Grid.BCtypeW[j] == BC_MOVING_WALL) {
	    // WEST face of cell (i+1,j) is a MOVINGWALL boundary.
	    Wl = MovingWallHeatFlux(Wr,Grid.nfaceW(i+1,j),Vwall.x);
	  } else if (Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL) {
	    // WEST face of cell (i+1,j) is a MOVINGWALL_ISOTHERMAL boundary.
	    Wl = MovingWallIsothermal(Wr,Grid.nfaceW(i+1,j),Vwall.x,Twall);
	  } else if (Grid.BCtypeW[j] == BC_RINGLEB_FLOW) {
	    // WEST face of cell (i+1,j) is a RINGLEB_FLOW boundary.
	    Wl = RinglebFlow(Wl,Grid.xfaceW(i+1,j));
	  } else {
	    // WEST face of cell (i+1,j) is a CHARACTERISTIC boundary.
	    Wl = BC_Characteristic_Pressure(Wr,WoW[j],Grid.nfaceW(i+1,j));
	  }

	} else if (i == ICu &&
		   (Grid.BCtypeE[j] == BC_REFLECTION ||
		    Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
		    Grid.BCtypeE[j] == BC_MASS_INJECTION ||
		    Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		    Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		    Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
		    Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
		    Grid.BCtypeE[j] == BC_RINGLEB_FLOW ||
		    Grid.BCtypeE[j] == BC_CHARACTERISTIC)) {

	  dX = Grid.xfaceE(i,j) - Grid.Cell[i][j].Xc;
	  Wl = W[i][j] + ( (phi[i][j]^dWdx[i][j])*dX.x +
			   (phi[i][j]^dWdy[i][j])*dX.y);

	  // EAST face of cell (i,j) is a normal boundary.
	  if (Grid.BCtypeE[j] == BC_REFLECTION) {
	    // EAST face of cell (i,j) is a REFLECTION boundary.
	    Wr = Reflect(Wl,Grid.nfaceE(i,j));
	  } else if (Grid.BCtypeE[j] == BC_BURNING_SURFACE) {
	    // EAST face of cell (i,j) is a BURNING_SURFACE boundary.
	    Wr = BurningSurface(Wl,Grid.nfaceE(i,j));
	  } else if (Grid.BCtypeE[j] == BC_MASS_INJECTION) {
	    // EAST face of cell (i,j) is a MASS_INJECTION boundary.
	    //Wr = MassInjection(Wl,Grid.nfaceE(i,j));
	    Wr = MassInjection2(Wl,Grid.xfaceE(i,j),Grid.nfaceE(i,j),Twall);
	    flux_function = IP.i_Flux_Function;
	    IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV_WRS;
	  } else if (Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX) {
	    // EAST face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
	    Wr = WallViscousHeatFlux(Wl,Grid.nfaceE(i,j));
	  } else if (Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    // EAST face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	    Wr = WallViscousIsothermal(Wl,Grid.nfaceE(i,j),Twall);
	  } else if (Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX) {
	    // EAST face of cell (i,j) is a MOVINGWALL_HEATFLUX boundary.
	    Wr = MovingWallHeatFlux(Wl,Grid.nfaceE(i,j),Vwall.x);
	  } else if (Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL) {
	    // EAST face of cell (i,j) is a MOVINGWALL_ISOTHERMAL boundary.
	    Wr = MovingWallIsothermal(Wl,Grid.nfaceE(i,j),Vwall.x,Twall);
	  } else if (Grid.BCtypeE[j] == BC_RINGLEB_FLOW) {
	    // EAST face of cell (i,j) is a RINGLEB_FLOW boundary.
	    Wr = RinglebFlow(Wr,Grid.xfaceE(i,j));
	  } else {
	    // EAST face of cell (i,j) is a CHARACTERISTIC boundary.
	    Wr = BC_Characteristic_Pressure(Wl,WoE[j],Grid.nfaceE(i,j));
	  }

	} else {

	  // EAST face is either a normal cell or possibly a FIXED, 
	  // NONE or EXTRAPOLATION boundary.
	  dX = Grid.xfaceE(i  ,j) - Grid.Cell[i  ][j].Xc;
	  Wl = W[i  ][j] + ( (phi[i  ][j]^dWdx[i  ][j])*dX.x +
			     (phi[i  ][j]^dWdy[i  ][j])*dX.y);
	  dX = Grid.xfaceW(i+1,j) - Grid.Cell[i+1][j].Xc;
	  Wr = W[i+1][j] + ( (phi[i+1][j]^dWdx[i+1][j])*dX.x +
			     (phi[i+1][j]^dWdy[i+1][j])*dX.y);

	}

	// Determine EAST face INVISCID flux.
	switch(IP.i_Flux_Function) {
	case FLUX_FUNCTION_GODUNOV :
	  Flux = FluxGodunov_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_ROE :
	  Flux = FluxRoe_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_RUSANOV :
	  Flux = FluxRusanov_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_HLLE :
	  Flux = FluxHLLE_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_HLLL :
	  Flux = FluxHLLL_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_HLLC :
	  Flux = FluxHLLC_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_VANLEER :
	  Flux = FluxVanLeer_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_AUSM :
	  Flux = FluxAUSM_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_AUSMplus :
	  Flux = FluxAUSMplus_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_GODUNOV_WRS :
	  Flux = FluxGodunov_Wrs_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	default:
	  Flux = FluxRoe_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	};

	if (IP.i_Flux_Function == FLUX_FUNCTION_GODUNOV_WRS) {
	  IP.i_Flux_Function = flux_function;
	}

	// Compute the cell centred stress tensor and heat flux vector if required.
 	if (Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
	    (Flow_Type && Axisymmetric)) {
	  W[i][j].ComputeViscousTerms(dWdx[i][j],
				      dWdy[i][j],
				      Grid.Cell[i][j].Xc,
				      Axisymmetric,
				      OFF,
				      Wall[i][j].ywall,
				      Wall[i][j].yplus);
	  U[i][j].tau = W[i][j].tau;
	  U[i][j].q = W[i][j].q;
	}

	// Evaluate the cell interface i-direction VISCOUS flux if necessary.
	if (Flow_Type) {
	  // Determine the EAST face VISCOUS flux.
	  if (i == ICl-1 && 
	      (Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
	       Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	       Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
	       Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
	       Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
	       Grid.BCtypeW[j] == BC_MASS_INJECTION)) {
	    // WEST face of cell (i+1,j) is a normal boundary.
	    Xr = Grid.Cell[i+1][j].Xc; Wr = W[i+1][j];
	    if (Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX) {
	      // WEST face of cell (i+1,j) is a WALL_VISCOUS_HEATFLUX boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	      //Wu.rho = Wr.rho; Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	      Wu = HALF*(Wr+W[i+1][j+1]); Wu.v = Vector2D_ZERO;
	      Wd = HALF*(Wr+W[i+1][j-1]); Wd.v = Vector2D_ZERO;
	    } else if (Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	      // WEST face of cell (i+1,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	      //Wu.rho = Wr.p/(Wr.R*Twall); Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	      Wu = HALF*(Wr+W[i+1][j+1]); Wu.v = Vector2D_ZERO; Wu.rho = Wu.p/(Wr.R*Twall);
	      Wd = HALF*(Wr+W[i+1][j-1]); Wd.v = Vector2D_ZERO; Wd.rho = Wd.p/(Wr.R*Twall);
	    } else if (Grid.BCtypeW[j] == BC_MOVING_WALL) {
	      // WEST face of cell (i+1,j) is a MOVINGWALL_HEATFLUX boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	      //Wu.rho = Wr.rho; Wu.v.x = Vwall.x; Wu.v.y = Vwall.y; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	      Wu = HALF*(Wr+W[i+1][j+1]); Wu.v = Vwall;
	      Wd = HALF*(Wr+W[i+1][j-1]); Wd.v = Vwall;
	    } else if (Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL) {
	      // WEST face of cell (i+1,j) is a MOVINGWALL_ISOTHERMAL boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	      //Wu.rho = Wr.p/(Wr.R*Twall); Wu.v.x = Vwall.x; Wu.v.y = Vwall.y; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	      Wu = HALF*(Wr+W[i+1][j+1]); Wu.v = Vwall; Wu.rho = Wu.p/(Wr.R*Twall);
	      Wd = HALF*(Wr+W[i+1][j-1]); Wd.v = Vwall; Wd.rho = Wd.p/(Wr.R*Twall);
	    } else if (Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
	      // WEST face of cell (i+1,j) is a BURNING_SURFACE boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	      Wu = BurningSurface(Wr,Grid.nfaceW(i+1,j));
	      Wd = Wu;
	    } else if (Grid.BCtypeW[j] == BC_MASS_INJECTION) {
	      // WEST face of cell (i+1,j) is a MASS_INJECTION boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	      //Wu = MassInjection(Wr,Grid.nfaceW(i+1,j));
	      Wu = MassInjection2(Wr,Grid.xfaceW(i+1,j),Grid.nfaceW(i+1,j),Twall);
	      Wd = Wu;
	    }
	    switch(IP.i_Viscous_Reconstruction) {
	    case VISCOUS_RECONSTRUCTION_CARTESIAN :
	    case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	      Xu = Grid.Node[i+1][j+1].X;
	      Xd = Grid.Node[i+1][j  ].X;
	      break;
	    case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    case VISCOUS_RECONSTRUCTION_HYBRID :
	      Xu = Grid.xfaceW(i+1,j);
	      Xl = Xu; Wl = Wu;
	      dWdxr = dWdx[i+1][j]; dWdxl = dWdxr;
	      dWdyr = dWdy[i+1][j]; dWdyl = dWdyr;
	      break;
	    };
	  } else if (i == ICu &&
		     (Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		      Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		      Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
		      Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
		      Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
		      Grid.BCtypeE[j] == BC_MASS_INJECTION)) {
	    // EAST face of cell (i,j) is a normal boundary.
	    Xl = Grid.Cell[i][j].Xc; Wl = W[i][j];
	    if (Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX) {
	      // EAST face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	      //Wu.rho = Wl.rho; Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	      Wu = HALF*(Wl+W[i][j+1]); Wu.v = Vector2D_ZERO;
	      Wd = HALF*(Wl+W[i][j-1]); Wd.v = Vector2D_ZERO;
	    } else if (Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	      // EAST face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	      //Wu.rho = Wl.p/(Wl.R*Twall); Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	      Wu = HALF*(Wl+W[i][j+1]); Wu.v = Vector2D_ZERO; Wu.rho = Wu.p/(Wl.R*Twall);
	      Wd = HALF*(Wl+W[i][j-1]); Wd.v = Vector2D_ZERO; Wd.rho = Wd.p/(Wl.R*Twall);
	    } else if (Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX) {
	      // EAST face of cell (i,j) is a MOVINGWALL_HEATFLUX boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	      //Wu.rho = Wl.rho; Wu.v.x = Vwall.x; Wu.v.y = Vwall.y; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	      Wu = HALF*(Wl+W[i][j+1]); Wu.v = Vwall;
	      Wd = HALF*(Wl+W[i][j-1]); Wd.v = Vwall;
	    } else if (Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL) {
	      // EAST face of cell (i,j) is a MOVINGWALL_ISOTHERMAL boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	      //Wu.rho = Wl.p/(Wl.R*Twall); Wu.v.x = Vwall.x; Wu.v.y = Vwall.y; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	      Wu = HALF*(Wl+W[i][j+1]); Wu.v = Vwall; Wu.rho = Wu.p/(Wl.R*Twall);
	      Wd = HALF*(Wl+W[i][j-1]); Wd.v = Vwall; Wd.rho = Wd.p/(Wl.R*Twall);
	    } else if (Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL) {
	      // EAST face of cell (i,j) is a BURNING_SURFACE boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	      Wu = BurningSurface(Wr,Grid.nfaceE(i,j));
	      Wd = Wu;
	    } else if (Grid.BCtypeE[j] == BC_MASS_INJECTION) {
	      // EAST face of cell (i,j) is a MASS_INJECTION boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	      //Wu = MassInjection(Wr,Grid.nfaceE(i,j));
	      Wu = MassInjection2(Wr,Grid.xfaceE(i,j),Grid.nfaceE(i,j),Twall);
	      Wd = Wu;
	    }
	    switch(IP.i_Viscous_Reconstruction) {
	    case VISCOUS_RECONSTRUCTION_CARTESIAN :
	    case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	      Xu = Grid.Node[i+1][j+1].X;
	      Xd = Grid.Node[i+1][j  ].X;
	      break;
	    case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    case VISCOUS_RECONSTRUCTION_HYBRID :
	      Xu = Grid.xfaceE(i,j);
	      Xr = Xu; Wr = Wu;
	      dWdxl = dWdx[i][j]; dWdxr = dWdxl;
	      dWdyl = dWdy[i][j]; dWdyr = dWdyl;
	      break;
	    };
	  } else {
	    // EAST face is either a normal cell or possibly a non-
	    // viscous boundary condition.
	    Xl = Grid.Cell[i  ][j].Xc; Wl = W[i  ][j];
	    Xr = Grid.Cell[i+1][j].Xc; Wr = W[i+1][j];
	    switch(IP.i_Viscous_Reconstruction) {
	    case VISCOUS_RECONSTRUCTION_CARTESIAN :
	    case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	      viscous_bc_flag = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;
	      Xu = Grid.Node[i+1][j+1].X; Wu = WnNE(i,j);
	      Xd = Grid.Node[i+1][j  ].X; Wd = WnSE(i,j);
	      break;
	    case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    case VISCOUS_RECONSTRUCTION_HYBRID :
	      Xu = Grid.xfaceE(i,j); Wu = HALF*(Wl + Wr);
	      dWdxl = dWdx[i][j]; dWdxr = dWdy[i+1][j];
	      dWdyl = dWdy[i][j]; dWdyr = dWdy[i+1][j];
	      break;
	    };
	  }
	  // Compute the EAST face viscous flux.
	  switch(IP.i_Viscous_Reconstruction) {
	  case VISCOUS_RECONSTRUCTION_CARTESIAN :
	  case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	    Flux -= ViscousFluxDiamondPath_n(Grid.xfaceE(i,j),
					     Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
					     Grid.nfaceE(i,j),
					     Axisymmetric,
					     viscous_bc_flag,
					     Wall[i][j].ywall, 
					     Wall[i][j].yplus,
					     dWdx_face,dWdy_face);
	    if (IP.Solver_Type == IMPLICIT && face_grad_arrays_allocated) {
	      dWdx_faceE[i][j] = dWdx_face;
	      dWdy_faceE[i][j] = dWdy_face;
	    }
	    break;
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    Flux -= ViscousFlux_n(Xu,Wu,HALF*(dWdxl+dWdxr),HALF*(dWdyl+dWdyr),
				  Grid.nfaceE(i,j),Axisymmetric,
				  viscous_bc_flag,
				  Wall[i][j].ywall, 
				  Wall[i][j].yplus);
	    if (IP.Solver_Type == IMPLICIT && face_grad_arrays_allocated) {
	      dWdx_faceE[i][j] = HALF*(dWdxl+dWdxr);
	      dWdy_faceE[i][j] = HALF*(dWdyl+dWdyr);
	    }
	    break;
	  case VISCOUS_RECONSTRUCTION_HYBRID :
 	    Flux -= ViscousFluxHybrid_n(Xu,Wu,Xl,Wl,dWdxl,dWdyl,
 					Xr,Wr,dWdxr,dWdyr,
 					Grid.nfaceE(i,j),
 					Axisymmetric,
 					Wall[i][j].ywall,
 					Wall[i][j].yplus,
 					dWdx_face,dWdy_face);
	    if (IP.Solver_Type == IMPLICIT && face_grad_arrays_allocated) {
	      dWdx_faceE[i][j] = dWdx_face;
	      dWdy_faceE[i][j] = dWdy_face;
	    }
	    break;
	  };
	}

	// Evaluate cell-averaged solution changes.
	dUdt[i  ][j][0] -= Flux*Grid.lfaceE(i,j)/Grid.Cell[i][j].A;
	dUdt[i+1][j][0] += Flux*Grid.lfaceW(i+1,j)/Grid.Cell[i+1][j].A;

	// Include axisymmetric source terms if required.
	if (Axisymmetric) {
	  dUdt[i][j][0] += W[i][j].Si(Grid.Cell[i][j].Xc);
	  if (Flow_Type)
	    dUdt[i][j][0] += W[i][j].Sv(Grid.Cell[i][j].Xc,
					dWdy[i][j],
					Wall[i][j].ywall,
					Wall[i][j].yplus);
	}

	// Include turbulent production and destruction source term.
	if (Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
	  dUdt[i][j][0] += W[i][j].St(Grid.Cell[i][j].Xc,
				      W[i][j],
				      dWdx[i][j],
				      dWdy[i][j],
				      Axisymmetric,
				      Wall[i][j].ywall,
				      Wall[i][j].yplus);
	}

	// Save west and east face boundary flux.
	if (i == ICl-1) {
	  FluxW[j] = -Flux*Grid.lfaceW(i+1,j);
	} else if (i == ICu) {
	  FluxE[j] =  Flux*Grid.lfaceE(i,j);
	}

      }
    }

    if (j > JCl-1 && j < JCu+1) {
      dUdt[ICl-1][j][0] = NavierStokes2D_U_VACUUM;
      dUdt[ICu+1][j][0] = NavierStokes2D_U_VACUUM;
    }

  }

  // Add j-direction (eta-direction) fluxes.
  for (int i = ICl; i <= ICu; i++) {
    for (int j = JCl-1; j <= JCu; j++) {

      // Evaluate the cell interface j-direction INVISCID fluxes.
      if (j == JCl-1 && 
	  (Grid.BCtypeS[i] == BC_REFLECTION ||
	   Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
	   Grid.BCtypeS[i] == BC_MASS_INJECTION ||
	   Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
	   Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	   Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX ||
	   Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
	   Grid.BCtypeS[i] == BC_RINGLEB_FLOW ||
	   Grid.BCtypeS[i] == BC_CHARACTERISTIC)) {

	dX = Grid.xfaceS(i,j+1) - Grid.Cell[i][j+1].Xc;
	Wr = W[i][j+1] + ( (phi[i][j+1]^dWdx[i][j+1])*dX.x +
			   (phi[i][j+1]^dWdy[i][j+1])*dX.y);

	// SOUTH face of cell (i,j+1) is a normal boundary.
	if (Grid.BCtypeS[i] == BC_REFLECTION) {
	  // SOUTH face of cell (i,j+1) is a REFLECTION boundary.
	  Wl = Reflect(Wr,Grid.nfaceS(i,j+1));
	} else if (Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
	  // SOUTH face of cell (i,j+1) is a BURNING_SURFACE boundary.
	  Wl = BurningSurface(Wr,Grid.nfaceS(i,j+1));
	} else if (Grid.BCtypeS[i] == BC_MASS_INJECTION) {
	  // SOUTH face of cell (i,j+1) is a MASS_INJECTION boundary.
	  //Wl = MassInjection(Wr,Grid.nfaceS(i,j+1));
	  Wl = MassInjection2(Wr,Grid.xfaceS(i,j+1),Grid.nfaceS(i,j+1),Twall);
	  flux_function = IP.i_Flux_Function;
	  IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV_WRS;
	} else if (Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX) {
	  // SOUTH face of cell (i,j+1) is a WALL_VISCOUS_HEATFLUX boundary.
	  Wl = WallViscousHeatFlux(Wr,Grid.nfaceS(i,j+1));
	} else if (Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	  // SOUTH face of cell (i,j+1) is a WALL_VISCOUS_ISOTHERMAL boundary.
	  Wl = WallViscousIsothermal(Wr,Grid.nfaceS(i,j+1),Twall);
	} else if (Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX) {
	  // SOUTH face of cell (i,j+1) is a MOVINGWALL_HEATFLUX boundary.
	  Wl = MovingWallHeatFlux(Wr,Grid.nfaceS(i,j+1),Vwall.x);
	} else if (Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL) {
	  // SOUTH face of cell (i,j+1) is a MOVINGWALL_ISOTHERMAL boundary.
	  Wl = MovingWallIsothermal(Wr,Grid.nfaceS(i,j+1),Vwall.x,Twall);
	} else if (Grid.BCtypeS[i] == BC_RINGLEB_FLOW) {
	  // SOUTH face of cell (i,j+1) is a RINGLEB_FLOW boundary.
	  Wl = RinglebFlow(Wl,Grid.xfaceS(i,j+1));
	  Wl = BC_Characteristic_Pressure(Wr,Wl,Grid.nfaceS(i,j+1));
	} else if (Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
	  // SOUTH face of cell (i,j+1) is a CHARACTERISTIC boundary.
	  Wl = BC_Characteristic_Pressure(Wr,WoS[i],Grid.nfaceS(i,j+1));
	}

      } else if (j == JCu && 
		 (Grid.BCtypeN[i] == BC_REFLECTION ||
		  Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
		  Grid.BCtypeN[i] == BC_MASS_INJECTION ||
		  Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
		  Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		  Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
		  Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
		  Grid.BCtypeN[i] == BC_RINGLEB_FLOW ||
		  Grid.BCtypeN[i] == BC_CHARACTERISTIC)) {

	dX = Grid.xfaceN(i,j) - Grid.Cell[i][j].Xc;
	Wl = W[i][j] + ( (phi[i][j]^dWdx[i][j])*dX.x +
			 (phi[i][j]^dWdy[i][j])*dX.y);

	// NORTH face of cell (i,j) is a normal boundary.
	if (Grid.BCtypeN[i] == BC_REFLECTION) {
	  // NORTH face of cell (i,j) is a REFLECTION boundary.
	  Wr = Reflect(Wl,Grid.nfaceN(i,j));
	} else if (Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
	  // NORTH face of cell (i,j) is a BURNING_SURFACE boundary.
	  Wr = BurningSurface(Wl,Grid.nfaceN(i,j));
	} else if (Grid.BCtypeN[i] == BC_MASS_INJECTION) {
	  // NORTH face of cell (i,j) is a MASS_INJECTION boundary.
	  //Wr = MassInjection(Wl,Grid.nfaceN(i,j));
	  Wr = MassInjection2(Wl,Grid.xfaceN(i,j),Grid.nfaceN(i,j),Twall);
	  flux_function = IP.i_Flux_Function;
	  IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV_WRS;
	} else if (Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX) {
	  // NORTH face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
	  Wr = WallViscousHeatFlux(Wl,Grid.nfaceN(i,j));
	} else if (Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	  // NORTH face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	  Wr = WallViscousIsothermal(Wl,Grid.nfaceN(i,j),Twall);
	} else if (Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX) {
	  // NORTH face of cell (i,j) is a MOVINGWALL_HEATFLUX boundary.
	  Wr = MovingWallHeatFlux(Wl,Grid.nfaceN(i,j),Vwall.x);
	} else if (Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL) {
	  // NORTH face of cell (i,j) is a MOVINGWALL_ISOTHERMAL boundary.
	  Wr = MovingWallIsothermal(Wl,Grid.nfaceN(i,j),Vwall.x,Twall);
	} else if (Grid.BCtypeN[i] == BC_RINGLEB_FLOW) {
	  // NORTH face of cell (i,j) is a RINGLEB_FLOW boundary.
	  Wr = RinglebFlow(Wr,Grid.xfaceN(i,j));
	  Wr = BC_Characteristic_Pressure(Wl,Wr,Grid.nfaceN(i,j));
	} else {
	  // NORTH face of cell (i,j) is a CHARACTERISTIC boundary.
	  Wr = BC_Characteristic_Pressure(Wl,WoN[i],Grid.nfaceN(i,j));
	}

      } else {

	// NORTH face is either a normal cell or possibly a FIXED, 
	// NONE or EXTRAPOLATION boundary.
	dX = Grid.xfaceN(i,j  ) - Grid.Cell[i][j  ].Xc;
	Wl = W[i][j  ] + ( (phi[i][j  ]^dWdx[i][j  ])*dX.x +
			   (phi[i][j  ]^dWdy[i][j  ])*dX.y);
	dX = Grid.xfaceS(i,j+1) - Grid.Cell[i][j+1].Xc;
	Wr = W[i][j+1] + ( (phi[i][j+1]^dWdx[i][j+1])*dX.x +
			   (phi[i][j+1]^dWdy[i][j+1])*dX.y);

      }

      // Determine NORTH face inviscid flux.
      switch(IP.i_Flux_Function) {
      case FLUX_FUNCTION_GODUNOV :
	Flux = FluxGodunov_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_ROE :
	Flux = FluxRoe_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_RUSANOV :
	Flux = FluxRusanov_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_HLLE :
	Flux = FluxHLLE_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_HLLL :
	Flux = FluxHLLL_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_HLLC :
	Flux = FluxHLLC_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_VANLEER :
	Flux = FluxVanLeer_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_AUSM :
	Flux = FluxAUSM_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_AUSMplus :
	Flux = FluxAUSMplus_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_GODUNOV_WRS :
	Flux = FluxGodunov_Wrs_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      default:
	Flux = FluxRoe_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      };

      if (IP.i_Flux_Function == FLUX_FUNCTION_GODUNOV_WRS) {
	IP.i_Flux_Function = flux_function;
      }

      // Evaluate the cell interface j-direction VISCOUS flux if necessary.
      if (Flow_Type) {
	// Determine the NORTH face VISCOUS flux.
	if (j == JCl-1 && 
	    (Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
	     Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX ||
	     Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
	     Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
	     Grid.BCtypeS[i] == BC_MASS_INJECTION)) {
	  // SOUTH face of cell (i,j+1) is a normal boundary.
	  Xr = Grid.Cell[i][j+1].Xc; Wr = W[i][j+1];
	  if (Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX) {
	    // SOUTH face of cell (i,j+1) is a WALL_VISCOUS_HEATFLUX boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	    //Wu.rho = Wr.rho; Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	    Wu = HALF*(Wr+W[i-1][j+1]); Wu.v = Vector2D_ZERO;
	    Wd = HALF*(Wr+W[i+1][j+1]); Wd.v = Vector2D_ZERO;
	  } else if (Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    // SOUTH face of cell (i,j+1) is a WALL_VISCOUS_ISOTHERMAL boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	    //Wu.rho = Wr.p/(Wr.R*Twall); Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	    Wu = HALF*(Wr+W[i-1][j+1]); Wu.v = Vector2D_ZERO; Wu.rho = Wu.p/(Wr.R*Twall);
	    Wd = HALF*(Wr+W[i+1][j+1]); Wd.v = Vector2D_ZERO; Wd.rho = Wd.p/(Wr.R*Twall);
	  } else if (Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX) {
	    // SOUTH face of cell (i,j+1) is a MOVINGWALL_HEATFLUX boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	    //Wu.rho = Wr.rho; Wu.v.x = Vwall.x; Wu.v.y = Vwall.y; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	    Wu = HALF*(Wr+W[i-1][j+1]); Wu.v = Vwall;
	    Wd = HALF*(Wr+W[i+1][j+1]); Wd.v = Vwall;
	  } else if (Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL) {
	    // SOUTH face of cell (i,j+1) is a MOVINGWALL_ISOTHERMAL boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	    //Wu.rho = Wr.p/(Wr.R*Twall); Wu.v.x = Vwall.x; Wu.v.y = Vwall.y; Wu.p = Wr.p; Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	    Wu = HALF*(Wr+W[i-1][j+1]); Wu.v = Vwall; Wu.rho = Wu.p/(Wr.R*Twall);
	    Wd = HALF*(Wr+W[i+1][j+1]); Wd.v = Vwall; Wd.rho = Wd.p/(Wr.R*Twall);
	  } else if (Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
	    // SOUTH face of cell (i,j+1) is a BURNING_SURFACE boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	    Wu = BurningSurface(Wr,Grid.nfaceS(i,j+1));
	    Wd = Wu;
	  } else if (Grid.BCtypeS[i] == BC_MASS_INJECTION) {
	    // SOUTH face of cell (i,j+1) is a MASS_INJECTION boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	    //Wu = MassInjection(Wr,Grid.nfaceS(i,j+1));
	    Wu = MassInjection2(Wr,Grid.xfaceS(i,j+1),Grid.nfaceS(i,j+1),Twall);
	    Wd = Wu;
	  }
	  switch(IP.i_Viscous_Reconstruction) {
	  case VISCOUS_RECONSTRUCTION_CARTESIAN :
	  case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	    Xu = Grid.Node[i  ][j+1].X;
	    Xd = Grid.Node[i+1][j+1].X;
	    break;
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	  case VISCOUS_RECONSTRUCTION_HYBRID :
	    Xu = Grid.xfaceS(i,j+1);
	    Xl = Xu; Wl = Wu;
	    dWdxr = dWdx[i][j+1]; dWdxl = dWdxr;
	    dWdyr = dWdy[i][j+1]; dWdyl = dWdyr;
	    break;
	  };
	} else if (j == JCu && 
		   (Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
		    Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		    Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
		    Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
		    Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
		    Grid.BCtypeN[i] == BC_MASS_INJECTION)) {
	  // NORTH face of cell (i,j) is a normal boundary.
	  Xl = Grid.Cell[i][j].Xc; Wl = W[i][j];
	  if (Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX) {
	    // NORTH face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	    //Wu.rho = Wl.rho; Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	    Wu = HALF*(Wl+W[i-1][j]); Wu.v = Vector2D_ZERO;
	    Wd = HALF*(Wl+W[i+1][j]); Wd.v = Vector2D_ZERO;
	  } else if (Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    // NORTH face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	    //Wu.rho = Wl.p/(Wl.R*Twall); Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	    Wu = HALF*(Wl+W[i-1][j]); Wu.v = Vector2D_ZERO; Wu.rho = Wu.p/(Wl.R*Twall);
	    Wd = HALF*(Wl+W[i+1][j]); Wd.v = Vector2D_ZERO; Wd.rho = Wd.p/(Wl.R*Twall);
	  } else if (Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX) {
	    // NORTH face of cell (i,j) is a MOVINGWALL_HEATFLUX boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	    //Wu.rho = Wl.rho; Wu.v.x = Vwall.x; Wu.v.y = Vwall.y; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	    Wu = HALF*(Wl+W[i-1][j]); Wu.v = Vwall;
	    Wd = HALF*(Wl+W[i+1][j]); Wd.v = Vwall;
	  } else if (Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL) {
	    // NORTH face of cell (i,j) is a MOVINGWALL_ISOTHERMAL boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	    //Wu.rho = Wl.p/(Wl.R*Twall); Wu.v.x = Vwall.x; Wu.v.y = Vwall.y; Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	    Wu = HALF*(Wl+W[i-1][j]); Wu.v = Vwall; Wu.rho = Wu.p/(Wl.R*Twall);
	    Wd = HALF*(Wl+W[i+1][j]); Wd.v = Vwall; Wd.rho = Wd.p/(Wl.R*Twall);
	  } else if (Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
	    // NORTH face of cell (i,j) is a BURNING_SURFACE boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	    Wu = BurningSurface(Wl,Grid.nfaceN(i,j));
	    Wd = Wu;
	  } else if (Grid.BCtypeN[i] == BC_MASS_INJECTION) {
	    // NORTH face of cell (i,j) is a MASS_INJECTION boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	    //Wu = MassInjection(Wl,Grid.nfaceN(i,j));
	    Wu = MassInjection2(Wl,Grid.xfaceN(i,j),Grid.nfaceN(i,j),Twall);
	    Wd = Wu;
	  }
	  switch(IP.i_Viscous_Reconstruction) {
	  case VISCOUS_RECONSTRUCTION_CARTESIAN :
	  case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	    Xu = Grid.Node[i  ][j+1].X;
	    Xd = Grid.Node[i+1][j+1].X;
	    break;
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	  case VISCOUS_RECONSTRUCTION_HYBRID :
	    Xu = Grid.xfaceN(i,j);
	    Xr = Xu; Wr = Wu;
	    dWdxl = dWdx[i][j]; dWdxr = dWdxl;
	    dWdyl = dWdy[i][j]; dWdyr = dWdyl;
	    break;
	  };
	} else {
	  // NORTH face is either a normal cell or possibly a non-viscous
	  // boundary condition.
	  Xl = Grid.Cell[i][j  ].Xc; Wl = W[i][j  ];
	  Xr = Grid.Cell[i][j+1].Xc; Wr = W[i][j+1];
	  switch(IP.i_Viscous_Reconstruction) {
	  case VISCOUS_RECONSTRUCTION_CARTESIAN :
	  case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	    viscous_bc_flag = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;
	    Xu = Grid.Node[i  ][j+1].X; Wu = WnNW(i,j);
	    Xd = Grid.Node[i+1][j+1].X; Wd = WnNE(i,j);
	    break;
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	  case VISCOUS_RECONSTRUCTION_HYBRID :
	    Xu = Grid.xfaceN(i,j); Wu = HALF*(Wl + Wr);
	    dWdxl = dWdx[i][j]; dWdxr = dWdy[i][j+1];
	    dWdyl = dWdy[i][j]; dWdyr = dWdy[i][j+1];
	    break;
	  };
	}
	// Compute the NORTH face viscous flux.
	switch(IP.i_Viscous_Reconstruction) {
	case VISCOUS_RECONSTRUCTION_CARTESIAN :
	case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	  Flux -= ViscousFluxDiamondPath_n(Grid.xfaceN(i,j),
					   Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
					   Grid.nfaceN(i,j),
					   Axisymmetric,
					   viscous_bc_flag,
					   Wall[i][j].ywall,
					   Wall[i][j].yplus,
					   dWdx_face,dWdy_face);
	  if (IP.Solver_Type == IMPLICIT && face_grad_arrays_allocated) {
	    dWdx_faceN[i][j] = dWdx_face;
	    dWdy_faceN[i][j] = dWdy_face;
	  }
	  break;
	case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	  Flux -= ViscousFlux_n(Xu,Wu,HALF*(dWdxl+dWdxr),HALF*(dWdyl+dWdyr),
				Grid.nfaceN(i,j),Axisymmetric,OFF,
				Wall[i][j].ywall,
				Wall[i][j].yplus);
	  if (IP.Solver_Type == IMPLICIT && face_grad_arrays_allocated) {
	    dWdx_faceN[i][j] = HALF*(dWdxl+dWdxr);
	    dWdy_faceN[i][j] = HALF*(dWdyl+dWdyr);
	  }
	  break;
	case VISCOUS_RECONSTRUCTION_HYBRID :
 	  Flux -= ViscousFluxHybrid_n(Xu,Wu,Xl,Wl,dWdxl,dWdyl,
 				      Xr,Wr,dWdxr,dWdyr,
 				      Grid.nfaceN(i,j),
 				      Axisymmetric,
 				      Wall[i][j].ywall,
 				      Wall[i][j].yplus,
 				      dWdx_face,dWdy_face);
	  if (IP.Solver_Type == IMPLICIT && face_grad_arrays_allocated) {
	    dWdx_faceN[i][j] = dWdx_face;
	    dWdy_faceN[i][j] = dWdy_face;
	  }
	  break;
	};
      }

      // Evaluate cell-averaged solution changes.
      dUdt[i][j  ][0] -= Flux*Grid.lfaceN(i,j)/Grid.Cell[i][j].A;
      dUdt[i][j+1][0] += Flux*Grid.lfaceS(i,j+1)/Grid.Cell[i][j+1].A;

      // Save south and north face boundary flux.
      if (j == JCl-1) {
	FluxS[i] = -Flux*Grid.lfaceS(i,j+1);
      } else if (j == JCu) {
	FluxN[i] = Flux*Grid.lfaceN(i,j);
      }

    }

    dUdt[i][JCl-1][0] = NavierStokes2D_U_VACUUM;
    dUdt[i][JCu+1][0] = NavierStokes2D_U_VACUUM;

  }

  // Zero the residuals for the turbulence variables according to
  // the turbulence boundary condition.
  error_flag = Turbulence_Zero_Residual(*this,
					1,IP);
  if (error_flag) return error_flag;

  // Residual successfully evaluated.
  return 0;

}

/**********************************************************************
 * Routine: dUdt_Multistage_Explicit                                  *
 *                                                                    *
 * This routine determines the solution residuals for a given stage   *
 * of a variety of multi-stage explicit time integration schemes for  *
 * a given solution block.                                            *
 *                                                                    *
 **********************************************************************/
int NavierStokes2D_Quad_Block::dUdt_Multistage_Explicit(const int &i_stage,
							NavierStokes2D_Input_Parameters &IP) {

  int error_flag, k_residual;
  double omega;
  Vector2D dX;
  NavierStokes2D_pState Wl, Wr;
  NavierStokes2D_cState Flux;

  NavierStokes2D_pState Wu, Wd, dWdxl, dWdyl, dWdxr, dWdyr;
  Vector2D Xl, Xr, Xu, Xd;
  int viscous_bc_flag;

  int flux_function = IP.i_Flux_Function;

  NavierStokes2D_cState NavierStokes2D_U_VACUUM; NavierStokes2D_U_VACUUM.Vacuum();
  NavierStokes2D_pState NavierStokes2D_W_STDATM; NavierStokes2D_W_STDATM.Standard_Atmosphere();

  // Evaluate the solution residual for stage i_stage of n_stage scheme.

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

  // Perform the linear reconstruction within each cell of the
  // computational grid for this stage.
  switch(IP.i_Reconstruction) {
  case RECONSTRUCTION_GREEN_GAUSS :
    Linear_Reconstruction_GreenGauss(*this,
				     IP.i_Limiter);
    break;
  case RECONSTRUCTION_LINEAR_LEAST_SQUARES :
    Linear_Reconstruction_LeastSquares(*this,
				       IP.i_Limiter);
    break;
  default:
    Linear_Reconstruction_LeastSquares(*this,
				       IP.i_Limiter);
    break;
  };

  // Evaluate the time rate of change of the solution (i.e., the
  // solution residuals) using a second-order limited upwind scheme
  // with a variety of flux functions.

  // Add i-direction (zeta-direction) fluxes.
  for (int j = JCl-1; j <= JCu+1; j++) {
    if (i_stage == 1) {
      Uo[ICl-1][j] = U[ICl-1][j];
      dUdt[ICl-1][j][k_residual] = NavierStokes2D_U_VACUUM;
    } else {
      dUdt[ICl-1][j][k_residual] = NavierStokes2D_U_VACUUM;
    }

    for (int i = ICl-1; i <= ICu; i++) {
      if (i_stage == 1) {
	Uo[i+1][j] = U[i+1][j];
	dUdt[i+1][j][k_residual] = NavierStokes2D_U_VACUUM;
      } else if (j > JCl-1 && j < JCu+1) {
	switch(IP.i_Time_Integration) {
	case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
	  //dUdt[i+1][j][k_residual] = dUdt[i+1][j][k_residual];
	  break;
	case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
	  if (IP.N_Stage == 2) {
	    //dUdt[i+1][j][k_residual] = dUdt[i+1][j][k_residual];
	  } else if (IP.N_Stage == 4 && i_stage == 4) {
	    dUdt[i+1][j][k_residual] = dUdt[i+1][j][0] + 
	      TWO*dUdt[i+1][j][1] +
	      TWO*dUdt[i+1][j][2];
	  } else {
	    dUdt[i+1][j][k_residual] = NavierStokes2D_U_VACUUM;
	  }
	  break;
	case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
	  dUdt[i+1][j][k_residual] = NavierStokes2D_U_VACUUM;
	  break;
	default:
	  dUdt[i+1][j][k_residual] = NavierStokes2D_U_VACUUM;
	  break;
	};
      }

      if (j >= JCl && j <= JCu) {

	// Evaluate the cell interface i-direction INVISCID fluxes.
	if (i == ICl-1 && 
	    (Grid.BCtypeW[j] == BC_REFLECTION ||
	     Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
	     Grid.BCtypeW[j] == BC_MASS_INJECTION ||
	     Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
	     Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
	     Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
	     Grid.BCtypeW[j] == BC_RINGLEB_FLOW ||
	     Grid.BCtypeW[j] == BC_CHARACTERISTIC)) {

	  dX = Grid.xfaceW(i+1,j) - Grid.Cell[i+1][j].Xc;
	  Wr = W[i+1][j] + ( (phi[i+1][j]^dWdx[i+1][j])*dX.x +
			     (phi[i+1][j]^dWdy[i+1][j])*dX.y);

	  // WEST face of cell (i+1,j) is a normal boundary.
	  if (Grid.BCtypeW[j] == BC_REFLECTION) {
	    // WEST face of cell (i+1,j) is a REFLECTION boundary.
	    Wl = Reflect(Wr,Grid.nfaceW(i+1,j));
	  } else if (Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
	    // WEST face of cell (i+1,j) is a BURNING_SURFACE boundary.
	    Wl = BurningSurface(Wr,Grid.nfaceW(i+1,j));
	  } else if (Grid.BCtypeW[j] == BC_MASS_INJECTION) {
	    // WEST face of cell (i+1,j) is a MASS_INJECTION boundary.
	    //Wl = MassInjection(Wr,Grid.nfaceW(i+1,j));
	    Wl = MassInjection2(Wr,Grid.xfaceW(i+1,j),Grid.nfaceW(i+1,j),Twall);
	    flux_function = IP.i_Flux_Function;
	    IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV_WRS;
	  } else if (Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX) {
	    // WEST face of cell (i+1,j) is a WALL_VISCOUS_HEATFLUX boundary.
	    Wl = WallViscousHeatFlux(Wr,Grid.nfaceW(i+1,j));
	  } else if (Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    // WEST face of cell (i+1,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	    Wl = WallViscousIsothermal(Wr,Grid.nfaceW(i+1,j),Twall);
	  } else if (Grid.BCtypeW[j] == BC_MOVING_WALL) {
	    // WEST face of cell (i+1,j) is a MOVINGWALL_HEATFLUX boundary.
	    Wl = MovingWallHeatFlux(Wr,Grid.nfaceW(i+1,j),Vwall.x);
	  } else if (Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL) {
	    // WEST face of cell (i+1,j) is a MOVINGWALL_ISOTHERMAL boundary.
	    Wl = MovingWallIsothermal(Wr,Grid.nfaceW(i+1,j),Vwall.x,Twall);
	  } else if (Grid.BCtypeW[j] == BC_RINGLEB_FLOW) {
	    // WEST face of cell (i+1,j) is a RINGLEB_FLOW boundary.
	    Wl = RinglebFlow(Wl,Grid.xfaceW(i+1,j));
	  } else {
	    // WEST face of cell (i+1,j) is a CHARACTERISTIC boundary.
	    Wl = BC_Characteristic_Pressure(Wr,WoW[j],Grid.nfaceW(i+1,j));
	  }

	} else if (i == ICu &&
		   (Grid.BCtypeE[j] == BC_REFLECTION ||
		    Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
		    Grid.BCtypeE[j] == BC_MASS_INJECTION ||
		    Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		    Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		    Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
		    Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
		    Grid.BCtypeE[j] == BC_RINGLEB_FLOW ||
		    Grid.BCtypeE[j] == BC_CHARACTERISTIC)) {

	  dX = Grid.xfaceE(i,j) - Grid.Cell[i][j].Xc;
	  Wl = W[i][j] + ( (phi[i][j]^dWdx[i][j])*dX.x +
			   (phi[i][j]^dWdy[i][j])*dX.y);

	  // EAST face of cell (i,j) is a normal boundary.
	  if (Grid.BCtypeE[j] == BC_REFLECTION) {
	    // EAST face of cell (i,j) is a REFLECTION boundary.
	    Wr = Reflect(Wl,Grid.nfaceE(i,j));
	  } else if (Grid.BCtypeE[j] == BC_BURNING_SURFACE) {
	    // EAST face of cell (i,j) is a BURNING_SURFACE boundary.
	    Wr = BurningSurface(Wl,Grid.nfaceE(i,j));
	  } else if (Grid.BCtypeE[j] == BC_MASS_INJECTION) {
	    // EAST face of cell (i,j) is a MASS_INJECTION boundary.
	    //Wr = MassInjection(Wl,Grid.nfaceE(i,j));
	    Wr = MassInjection2(Wl,Grid.xfaceE(i,j),Grid.nfaceE(i,j),Twall);
	    flux_function = IP.i_Flux_Function;
	    IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV_WRS;
	  } else if (Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX) {
	    // EAST face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
	    Wr = WallViscousHeatFlux(Wl,Grid.nfaceE(i,j));
	  } else if (Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    // EAST face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	    Wr = WallViscousIsothermal(Wl,Grid.nfaceE(i,j),Twall);
	  } else if (Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX) {
	    // EAST face of cell (i,j) is a MOVINGWALL_HEATFLUX boundary.
	    Wr = MovingWallHeatFlux(Wl,Grid.nfaceE(i,j),Vwall.x);
	  } else if (Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL) {
	    // EAST face of cell (i,j) is a MOVINGWALL_ISOTHERMAL boundary.
	    Wr = MovingWallIsothermal(Wl,Grid.nfaceE(i,j),Vwall.x,Twall);
	  } else if (Grid.BCtypeE[j] == BC_RINGLEB_FLOW) {
	    // EAST face of cell (i,j) is a RINGLEB_FLOW boundary.
	    Wr = RinglebFlow(Wr,Grid.xfaceE(i,j));
	  } else {
	    // EAST face of cell (i,j) is a CHARACTERISTIC boundary.
	    Wr = BC_Characteristic_Pressure(Wl,WoE[j],Grid.nfaceE(i,j));
	  }

	} else {

	  // EAST face is either a normal cell or possibly a FIXED, 
	  // NONE or EXTRAPOLATION boundary.
	  dX = Grid.xfaceE(i  ,j) - Grid.Cell[i  ][j].Xc;
	  Wl = W[i  ][j] + ( (phi[i  ][j]^dWdx[i  ][j])*dX.x +
			     (phi[i  ][j]^dWdy[i  ][j])*dX.y);
	  dX = Grid.xfaceW(i+1,j) - Grid.Cell[i+1][j].Xc;
	  Wr = W[i+1][j] + ( (phi[i+1][j]^dWdx[i+1][j])*dX.x +
			     (phi[i+1][j]^dWdy[i+1][j])*dX.y);

	}

	// Determine EAST face INVISCID flux.
	switch(IP.i_Flux_Function) {
	case FLUX_FUNCTION_GODUNOV :
	  Flux = FluxGodunov_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_ROE :
	  Flux = FluxRoe_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_RUSANOV :
	  Flux = FluxRusanov_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_HLLE :
	  Flux = FluxHLLE_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_HLLL :
	  Flux = FluxHLLL_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_HLLC :
	  Flux = FluxHLLC_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_VANLEER :
	  Flux = FluxVanLeer_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_AUSM :
	  Flux = FluxAUSM_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_AUSMplus :
	  Flux = FluxAUSMplus_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_GODUNOV_WRS :
	  Flux = FluxGodunov_Wrs_n(Wl,Wr,Grid.nfaceN(i,j));
	  break;
	default:
	  Flux = FluxRoe_n(Wl,Wr,Grid.nfaceE(i,j));
	  break;
	};

	if (IP.i_Flux_Function == FLUX_FUNCTION_GODUNOV_WRS) {
	  IP.i_Flux_Function = flux_function;
	}

 	// Compute the cell centred stress tensor and heat flux vector if required.
 	if (Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
	    (Flow_Type && Axisymmetric)) {
 	  W[i][j].ComputeViscousTerms(dWdx[i][j],
				      dWdy[i][j],
				      Grid.Cell[i][j].Xc,
				      Axisymmetric,
				      OFF,
				      Wall[i][j].ywall,
				      Wall[i][j].yplus);
 	  U[i][j].tau = W[i][j].tau;
 	  U[i][j].q = W[i][j].q;
 	}

	// Evaluate the cell interface i-direction VISCOUS flux if necessary.
	if (Flow_Type) {
	  // Determine the EAST face VISCOUS flux.
	  if (i == ICl-1 && 
	      (Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX ||
	       Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
	       Grid.BCtypeW[j] == BC_MOVING_WALL_HEATFLUX ||
	       Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL ||
	       Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
	       Grid.BCtypeW[j] == BC_MASS_INJECTION)) {
	    // WEST face of cell (i+1,j) is a normal boundary.
	    Xr = Grid.Cell[i+1][j].Xc; Wr = W[i+1][j];
	    if (Grid.BCtypeW[j] == BC_WALL_VISCOUS_HEATFLUX) {
	      // WEST face of cell (i+1,j) is a WALL_VISCOUS_HEATFLUX boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	      //Wu.rho = Wr.rho; Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wr.p; 
	      //Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	      Wu = HALF*(Wr+W[i+1][j+1]); Wu.v = Vector2D_ZERO;
	      Wd = HALF*(Wr+W[i+1][j-1]); Wd.v = Vector2D_ZERO;
	    } else if (Grid.BCtypeW[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	      // WEST face of cell (i+1,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	      //Wu.rho = Wr.p/(Wr.R*Twall); Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wr.p;
	      //Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	      Wu = HALF*(Wr+W[i+1][j+1]); Wu.v = Vector2D_ZERO; Wu.rho = Wu.p/(Wr.R*Twall);
	      Wd = HALF*(Wr+W[i+1][j-1]); Wd.v = Vector2D_ZERO; Wd.rho = Wd.p/(Wr.R*Twall);
	    } else if (Grid.BCtypeW[j] == BC_MOVING_WALL) {
	      // WEST face of cell (i+1,j) is a MOVINGWALL_HEATFLUX boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	      //Wu.rho = Wr.rho; Wu.v.x = Vwall.x; Wu.v.y = Vwall.y; Wu.p = Wr.p;
	      //Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	      Wu = HALF*(Wr+W[i+1][j+1]); Wu.v = Vwall;
	      Wd = HALF*(Wr+W[i+1][j-1]); Wd.v = Vwall;
	    } else if (Grid.BCtypeW[j] == BC_MOVING_WALL_ISOTHERMAL) {
	      // WEST face of cell (i+1,j) is a MOVINGWALL_ISOTHERMAL boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	      //Wu.rho = Wr.p/(Wr.R*Twall); Wu.v.x = Vwall.x; Wu.v.y = Vwall.y; Wu.p = Wr.p;
	      //Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	      Wu = HALF*(Wr+W[i+1][j+1]); Wu.v = Vwall; Wu.rho = Wu.p/(Wr.R*Twall);
	      Wd = HALF*(Wr+W[i+1][j-1]); Wd.v = Vwall; Wd.rho = Wd.p/(Wr.R*Twall);
	    } else if (Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
	      // WEST face of cell (i+1,j) is a BURNING_SURFACE boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	      Wu = BurningSurface(Wr,Grid.nfaceW(i+1,j));
	      Wd = Wu;
	    } else if (Grid.BCtypeW[j] == BC_MASS_INJECTION) {
	      // WEST face of cell (i+1,j) is a MASS_INJECTION boundary.
	      viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	      //Wu = MassInjection(Wr,Grid.nfaceW(i+1,j));
	      Wu = MassInjection2(Wr,Grid.xfaceW(i+1,j),Grid.nfaceW(i+1,j),Twall);
	      Wd = Wu;
	    }
	    switch(IP.i_Viscous_Reconstruction) {
	    case VISCOUS_RECONSTRUCTION_CARTESIAN :
	    case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	      Xu = Grid.Node[i+1][j+1].X;
	      Xd = Grid.Node[i+1][j  ].X; Wd = Wu;
	      break;
	    case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    case VISCOUS_RECONSTRUCTION_HYBRID :
	      Xu = Grid.xfaceW(i+1,j);
	      Xl = Xu; Wl = Wu;
	      dWdxr = dWdx[i+1][j]; dWdxl = dWdxr;
	      dWdyr = dWdy[i+1][j]; dWdyl = dWdyr;
	      break;
	    };
	  } else if (i == ICu &&
		     (Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX ||
		      Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL ||
		      Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX ||
		      Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL ||
		      Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
		      Grid.BCtypeE[j] == BC_MASS_INJECTION)) {
	    // EAST face of cell (i,j) is a normal boundary.
	    Xl = Grid.Cell[i][j].Xc; Wl = W[i][j];
	    if (Grid.BCtypeE[j] == BC_WALL_VISCOUS_HEATFLUX) {
	      // EAST face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	      //Wu.rho = Wl.rho; Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wl.p;
	      //Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	      Wu = HALF*(Wl+W[i][j+1]); Wu.v = Vector2D_ZERO;
	      Wd = HALF*(Wl+W[i][j-1]); Wd.v = Vector2D_ZERO;
	    } else if (Grid.BCtypeE[j] == BC_WALL_VISCOUS_ISOTHERMAL) {
	      // EAST face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	      //Wu.rho = Wl.p/(Wl.R*Twall); Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wl.p;
	      //Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	      Wu = HALF*(Wl+W[i][j+1]); Wu.v = Vector2D_ZERO; Wu.rho = Wu.p/(Wl.R*Twall);
	      Wd = HALF*(Wl+W[i][j-1]); Wd.v = Vector2D_ZERO; Wd.rho = Wd.p/(Wl.R*Twall);
	    } else if (Grid.BCtypeE[j] == BC_MOVING_WALL_HEATFLUX) {
	      // EAST face of cell (i,j) is a MOVINGWALL_HEATFLUX boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	      //Wu.rho = Wl.rho; Wu.v.x = Vwall.x; Wu.v.y = Vwall.y; Wu.p = Wl.p;
	      //Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	      Wu = HALF*(Wl+W[i][j+1]); Wu.v = Vwall;
	      Wd = HALF*(Wl+W[i][j-1]); Wd.v = Vwall;
	    } else if (Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL) {
	      // EAST face of cell (i,j) is a MOVINGWALL_ISOTHERMAL boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	      //Wu.rho = Wl.p/(Wl.R*Twall); Wu.v.x = Vwall.x; Wu.v.y = Vwall.y; Wu.p = Wl.p;
	      //Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	      Wu = HALF*(Wl+W[i][j+1]); Wu.v = Vwall; Wu.rho = Wu.p/(Wl.R*Twall);
	      Wd = HALF*(Wl+W[i][j-1]); Wd.v = Vwall; Wd.rho = Wd.p/(Wl.R*Twall);
	    } else if (Grid.BCtypeE[j] == BC_MOVING_WALL_ISOTHERMAL) {
	      // EAST face of cell (i,j) is a BURNING_SURFACE boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	      Wu = BurningSurface(Wr,Grid.nfaceE(i,j));
	      Wd = Wu;
	    } else if (Grid.BCtypeE[j] == BC_MASS_INJECTION) {
	      // EAST face of cell (i,j) is a MASS_INJECTION boundary.
	      viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	      //Wu = MassInjection(Wr,Grid.nfaceE(i,j));
	      Wu = MassInjection2(Wr,Grid.xfaceE(i,j),Grid.nfaceE(i,j),Twall);
	      Wd = Wu;
	    }
	    switch(IP.i_Viscous_Reconstruction) {
	    case VISCOUS_RECONSTRUCTION_CARTESIAN :
	    case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	      Xu = Grid.Node[i+1][j+1].X;
	      Xd = Grid.Node[i+1][j  ].X;
	      break;
	    case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    case VISCOUS_RECONSTRUCTION_HYBRID :
	      Xu = Grid.xfaceE(i,j);
	      Xr = Xu; Wr = Wu;
	      dWdxl = dWdx[i][j]; dWdxr = dWdxl;
	      dWdyl = dWdy[i][j]; dWdyr = dWdyl;
	      break;
	    };
	  } else {
	    // EAST face is either a normal cell or possibly a non-
	    // viscous boundary condition.
	    Xl = Grid.Cell[i  ][j].Xc; Wl = W[i  ][j];
	    Xr = Grid.Cell[i+1][j].Xc; Wr = W[i+1][j];
	    switch(IP.i_Viscous_Reconstruction) {
	    case VISCOUS_RECONSTRUCTION_CARTESIAN :
	    case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	      viscous_bc_flag = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;
	      Xu = Grid.Node[i+1][j+1].X; Wu = WnNE(i,j);
	      Xd = Grid.Node[i+1][j  ].X; Wd = WnSE(i,j);
	      break;
	    case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    case VISCOUS_RECONSTRUCTION_HYBRID :
	      Xu = Grid.xfaceE(i,j); Wu = HALF*(Wl + Wr);
	      dWdxl = dWdx[i][j]; dWdxr = dWdy[i+1][j];
	      dWdyl = dWdy[i][j]; dWdyr = dWdy[i+1][j];
	      break;
	    };
	  }
	  // Compute the EAST face viscous flux.
	  switch(IP.i_Viscous_Reconstruction) {
	  case VISCOUS_RECONSTRUCTION_CARTESIAN :
	  case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	    Flux -= ViscousFluxDiamondPath_n(Grid.xfaceE(i,j),
					     Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
					     Grid.nfaceE(i,j),
					     Axisymmetric,
					     viscous_bc_flag,
					     Wall[i][j].ywall,
					     Wall[i][j].yplus);
	    break;
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	    Flux -= ViscousFlux_n(Xu,Wu,HALF*(dWdxl+dWdxr),
				  HALF*(dWdyl+dWdyr),
				  Grid.nfaceE(i,j),
				  Axisymmetric,OFF,
				  Wall[i][j].ywall,
				  Wall[i][j].yplus);
	    break;
	  case VISCOUS_RECONSTRUCTION_HYBRID :
 	    Flux -= ViscousFluxHybrid_n(Xu,Wu,Xl,Wl,dWdxl,dWdyl,
 					Xr,Wr,dWdxr,dWdyr,
 					Grid.nfaceE(i,j),
 					Axisymmetric,
 					Wall[i][j].ywall,
 					Wall[i][j].yplus);
	    break;
	  };
	}

	// Evaluate cell-averaged solution changes.
	dUdt[i  ][j][k_residual] -= (IP.CFL_Number*dt[i][j])*
	  Flux*Grid.lfaceE(i,j)/Grid.Cell[i][j].A;
	dUdt[i+1][j][k_residual] += (IP.CFL_Number*dt[i+1][j])*
	  Flux*Grid.lfaceW(i+1,j)/Grid.Cell[i+1][j].A;

   	// Include axisymmetric source terms if required.
       	if (Axisymmetric) {
   	  dUdt[i][j][k_residual] += (IP.CFL_Number*dt[i][j])*
	    W[i][j].Si(Grid.Cell[i][j].Xc);
	  if (Flow_Type)
  	    dUdt[i][j][k_residual] += (IP.CFL_Number*dt[i][j])*
	      W[i][j].Sv(Grid.Cell[i][j].Xc,
			 dWdy[i][j],
			 Wall[i][j].ywall,
			 Wall[i][j].yplus);
 	}

 	// Include turbulent production and destruction source term if required.
 	if (Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
 	  dUdt[i][j][k_residual] += (IP.CFL_Number*dt[i][j])*
	    W[i][j].St(Grid.Cell[i][j].Xc,
		       W[i][j],
		       dWdx[i][j],
		       dWdy[i][j],
		       Axisymmetric,
		       Wall[i][j].ywall,
		       Wall[i][j].yplus);
	}

	// Save west and east face boundary flux.
 	if (i == ICl-1) {
 	  FluxW[j] = -Flux*Grid.lfaceW(i+1,j);
 	} else if (i == ICu) {
 	  FluxE[j] =  Flux*Grid.lfaceE(i,j);
 	}

      }
    }

    if (j > JCl-1 && j < JCu+1) {
      dUdt[ICl-1][j][k_residual] = NavierStokes2D_U_VACUUM;
      dUdt[ICu+1][j][k_residual] = NavierStokes2D_U_VACUUM;
    }

  }

  // Add j-direction (eta-direction) fluxes.
  for (int i = ICl; i <= ICu; i++) {
    for (int j = JCl-1; j <= JCu; j++) {

      // Evaluate the cell interface j-direction INVISCID fluxes.
      if (j == JCl-1 && 
	  (Grid.BCtypeS[i] == BC_REFLECTION ||
	   Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
	   Grid.BCtypeS[i] == BC_MASS_INJECTION ||
	   Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
	   Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	   Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX ||
	   Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
	   Grid.BCtypeS[i] == BC_RINGLEB_FLOW ||
	   Grid.BCtypeS[i] == BC_CHARACTERISTIC)) {

	dX = Grid.xfaceS(i,j+1) - Grid.Cell[i][j+1].Xc;
	Wr = W[i][j+1] + ( (phi[i][j+1]^dWdx[i][j+1])*dX.x +
			   (phi[i][j+1]^dWdy[i][j+1])*dX.y);

	// SOUTH face of cell (i,j+1) is a normal boundary.
	if (Grid.BCtypeS[i] == BC_REFLECTION) {
	  // SOUTH face of cell (i,j+1) is a REFLECTION boundary.
	  Wl = Reflect(Wr,Grid.nfaceS(i,j+1));
	} else if (Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
	  // SOUTH face of cell (i,j+1) is a BURNING_SURFACE boundary.
	  Wl = BurningSurface(Wr,Grid.nfaceS(i,j+1));
	} else if (Grid.BCtypeS[i] == BC_MASS_INJECTION) {
	  // SOUTH face of cell (i,j+1) is a MASS_INJECTION boundary.
	  //Wl = MassInjection(Wr,Grid.nfaceS(i,j+1));
	  Wl = MassInjection2(Wr,Grid.xfaceS(i,j+1),Grid.nfaceS(i,j+1),Twall);
	  flux_function = IP.i_Flux_Function;
	  IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV_WRS;
	} else if (Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX) {
	  // SOUTH face of cell (i,j+1) is a WALL_VISCOUS_HEATFLUX boundary.
	  Wl = WallViscousHeatFlux(Wr,Grid.nfaceS(i,j+1));
	} else if (Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	  // SOUTH face of cell (i,j+1) is a WALL_VISCOUS_ISOTHERMAL boundary.
	  Wl = WallViscousIsothermal(Wr,Grid.nfaceS(i,j+1),Twall);
	} else if (Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX) {
	  // SOUTH face of cell (i,j+1) is a MOVINGWALL_HEATFLUX boundary.
	  Wl = MovingWallHeatFlux(Wr,Grid.nfaceS(i,j+1),Vwall.x);
	} else if (Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL) {
	  // SOUTH face of cell (i,j+1) is a MOVINGWALL_ISOTHERMAL boundary.
	  Wl = MovingWallIsothermal(Wr,Grid.nfaceS(i,j+1),Vwall.x,Twall);
	} else if (Grid.BCtypeS[i] == BC_RINGLEB_FLOW) {
	  // SOUTH face of cell (i,j+1) is a RINGLEB_FLOW boundary.
	  Wl = RinglebFlow(Wl,Grid.xfaceS(i,j+1));
	  Wl = BC_Characteristic_Pressure(Wr,Wl,Grid.nfaceS(i,j+1));
	} else if (Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
	  // SOUTH face of cell (i,j+1) is a CHARACTERISTIC boundary.
	  Wl = BC_Characteristic_Pressure(Wr,WoS[i],Grid.nfaceS(i,j+1));
	}

      } else if (j == JCu && 
		 (Grid.BCtypeN[i] == BC_REFLECTION ||
		  Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
		  Grid.BCtypeN[i] == BC_MASS_INJECTION ||
		  Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
		  Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		  Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
		  Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
		  Grid.BCtypeN[i] == BC_RINGLEB_FLOW ||
		  Grid.BCtypeN[i] == BC_CHARACTERISTIC)) {

	dX = Grid.xfaceN(i,j) - Grid.Cell[i][j].Xc;
	Wl = W[i][j] + ( (phi[i][j]^dWdx[i][j])*dX.x +
			 (phi[i][j]^dWdy[i][j])*dX.y );

	// NORTH face of cell (i,j) is a normal boundary.
	if (Grid.BCtypeN[i] == BC_REFLECTION) {
	  // NORTH face of cell (i,j) is a REFLECTION boundary.
	  Wr = Reflect(Wl,Grid.nfaceN(i,j));
	} else if (Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
	  // NORTH face of cell (i,j) is a BURNING_SURFACE boundary.
	  Wr = BurningSurface(Wl,Grid.nfaceN(i,j));
	} else if (Grid.BCtypeN[i] == BC_MASS_INJECTION) {
	  // NORTH face of cell (i,j) is a MASS_INJECTION boundary.
	  //Wr = MassInjection(Wl,Grid.nfaceN(i,j));
	  Wr = MassInjection2(Wl,Grid.xfaceN(i,j),Grid.nfaceN(i,j),Twall);
	  flux_function = IP.i_Flux_Function;
	  IP.i_Flux_Function = FLUX_FUNCTION_GODUNOV_WRS;
	} else if (Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX) {
	  // NORTH face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
	  Wr = WallViscousHeatFlux(Wl,Grid.nfaceN(i,j));
	} else if (Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	  // NORTH face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	  Wr = WallViscousIsothermal(Wl,Grid.nfaceN(i,j),Twall);
	} else if (Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX) {
	  // NORTH face of cell (i,j) is a MOVINGWALL_HEATFLUX boundary.
	  Wr = MovingWallHeatFlux(Wl,Grid.nfaceN(i,j),Vwall.x);
	} else if (Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL) {
	  // NORTH face of cell (i,j) is a MOVINGWALL_ISOTHERMAL boundary.
	  Wr = MovingWallIsothermal(Wl,Grid.nfaceN(i,j),Vwall.x,Twall);
	} else if (Grid.BCtypeN[i] == BC_RINGLEB_FLOW) {
	  // NORTH face of cell (i,j) is a RINGLEB_FLOW boundary.
	  Wr = RinglebFlow(Wr,Grid.xfaceN(i,j));
	  Wr = BC_Characteristic_Pressure(Wl,Wr,Grid.nfaceN(i,j));
	} else {
	  // NORTH face of cell (i,j) is a CHARACTERISTIC boundary.
	  Wr = BC_Characteristic_Pressure(Wl,WoN[i],Grid.nfaceN(i,j));
	}

      } else {
	
	// NORTH face is either a normal cell or possibly a FIXED, 
	// NONE or EXTRAPOLATION boundary.
	dX = Grid.xfaceN(i,j  ) - Grid.Cell[i][j  ].Xc;
	Wl = W[i][j  ] + ( (phi[i][j  ]^dWdx[i][j  ])*dX.x +
			   (phi[i][j  ]^dWdy[i][j  ])*dX.y);
	dX = Grid.xfaceS(i,j+1) - Grid.Cell[i][j+1].Xc;
	Wr = W[i][j+1] + ( (phi[i][j+1]^dWdx[i][j+1])*dX.x +
			   (phi[i][j+1]^dWdy[i][j+1])*dX.y);

      }

      // Determine NORTH face inviscid flux.
      switch(IP.i_Flux_Function) {
      case FLUX_FUNCTION_GODUNOV :
	Flux = FluxGodunov_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_ROE :
	Flux = FluxRoe_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_RUSANOV :
	Flux = FluxRusanov_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_HLLE :
	Flux = FluxHLLE_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_HLLL :
	Flux = FluxHLLL_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_HLLC :
	Flux = FluxHLLC_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_VANLEER :
	Flux = FluxVanLeer_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_AUSM :
	Flux = FluxAUSM_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_AUSMplus :
	Flux = FluxAUSMplus_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      case FLUX_FUNCTION_GODUNOV_WRS :
	Flux = FluxGodunov_Wrs_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      default:
	Flux = FluxRoe_n(Wl,Wr,Grid.nfaceN(i,j));
	break;
      };

      if (IP.i_Flux_Function == FLUX_FUNCTION_GODUNOV_WRS) {
	IP.i_Flux_Function = flux_function;
      }

      // Evaluate the cell interface j-direction VISCOUS flux if necessary.
      if (Flow_Type) {
	// Determine the NORTH face VISCOUS flux.
	if (j == JCl-1 && 
	    (Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX ||
	     Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
	     Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX ||
	     Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL ||
	     Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
	     Grid.BCtypeS[i] == BC_MASS_INJECTION)) {
	  // SOUTH face of cell (i,j+1) is a normal boundary.
	  Xr = Grid.Cell[i][j+1].Xc; Wr = W[i][j+1];
	  if (Grid.BCtypeS[i] == BC_WALL_VISCOUS_HEATFLUX) {
	    // SOUTH face of cell (i,j+1) is a WALL_VISCOUS_HEATFLUX boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	    //Wu.rho = Wr.rho; Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wr.p;
	    //Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	    Wu = HALF*(Wr+W[i-1][j+1]); Wu.v = Vector2D_ZERO;
	    Wd = HALF*(Wr+W[i+1][j+1]); Wd.v = Vector2D_ZERO;
	  } else if (Grid.BCtypeS[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    // SOUTH face of cell (i,j+1) is a WALL_VISCOUS_ISOTHERMAL boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	    //Wu.rho = Wr.p/(Wr.R*Twall); Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wr.p;
	    //Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	    Wu = HALF*(Wr+W[i-1][j+1]); Wu.v = Vector2D_ZERO; Wu.rho = Wu.p/(Wr.R*Twall);
	    Wd = HALF*(Wr+W[i+1][j+1]); Wd.v = Vector2D_ZERO; Wd.rho = Wd.p/(Wr.R*Twall);
	  } else if (Grid.BCtypeS[i] == BC_MOVING_WALL_HEATFLUX) {
	    // SOUTH face of cell (i,j+1) is a MOVINGWALL_HEATFLUX boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_HEATFLUX;
	    //Wu.rho = Wr.rho; Wu.v.x = Vwall.x; Wu.v.y = Vwall.y; Wu.p = Wr.p;
	    //Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	    Wu = HALF*(Wr+W[i-1][j+1]); Wu.v = Vwall;
	    Wd = HALF*(Wr+W[i+1][j+1]); Wd.v = Vwall;
	  } else if (Grid.BCtypeS[i] == BC_MOVING_WALL_ISOTHERMAL) {
	    // SOUTH face of cell (i,j+1) is a MOVINGWALL_ISOTHERMAL boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	    //Wu.rho = Wr.p/(Wr.R*Twall); Wu.v.x = Vwall.x; Wu.v.y = Vwall.y; Wu.p = Wr.p;
	    //Wu.k = Wr.k; Wu.omega = Wr.omega;Wu.ke = Wr.ke; Wu.ee = Wr.ee;
	    Wu = HALF*(Wr+W[i-1][j+1]); Wu.v = Vwall; Wu.rho = Wu.p/(Wr.R*Twall);
	    Wd = HALF*(Wr+W[i+1][j+1]); Wd.v = Vwall; Wd.rho = Wd.p/(Wr.R*Twall);
	  } else if (Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
	    // SOUTH face of cell (i,j+1) is a BURNING_SURFACE boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	    Wu = BurningSurface(Wr,Grid.nfaceS(i,j+1));
	    Wd = Wu;
	  } else if (Grid.BCtypeS[i] == BC_MASS_INJECTION) {
	    // SOUTH face of cell (i,j+1) is a MASS_INJECTION boundary.
	    viscous_bc_flag = DIAMONDPATH_RIGHT_TRIANGLE_ISOTHERMAL;
	    //Wu = MassInjection(Wr,Grid.nfaceS(i,j+1));
	    Wu = MassInjection2(Wr,Grid.xfaceS(i,j+1),Grid.nfaceS(i,j+1),Twall);
	    Wd = Wu;
	  }
	  switch(IP.i_Viscous_Reconstruction) {
	  case VISCOUS_RECONSTRUCTION_CARTESIAN :
	  case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	    Xu = Grid.Node[i  ][j+1].X;
	    Xd = Grid.Node[i+1][j+1].X;
	    break;
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	  case VISCOUS_RECONSTRUCTION_HYBRID :
	    Xu = Grid.xfaceS(i,j+1);
	    Xl = Xu; Wl = Wu;
	    dWdxr = dWdx[i][j+1]; dWdxl = dWdxr;
	    dWdyr = dWdy[i][j+1]; dWdyl = dWdyr;
	    break;
	  };
	} else if (j == JCu && 
		   (Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX ||
		    Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL ||
		    Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX ||
		    Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL ||
		    Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
		    Grid.BCtypeN[i] == BC_MASS_INJECTION)) {
	  // NORTH face of cell (i,j) is a normal boundary.
	  Xl = Grid.Cell[i][j].Xc; Wl = W[i][j];
	  if (Grid.BCtypeN[i] == BC_WALL_VISCOUS_HEATFLUX) {
	    // NORTH face of cell (i,j) is a WALL_VISCOUS_HEATFLUX boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	    //Wu.rho = Wl.rho; Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wl.p;
	    //Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	    Wu = HALF*(Wl+W[i-1][j]); Wu.v = Vector2D_ZERO;
	    Wd = HALF*(Wl+W[i+1][j]); Wd.v = Vector2D_ZERO;
	  } else if (Grid.BCtypeN[i] == BC_WALL_VISCOUS_ISOTHERMAL) {
	    // NORTH face of cell (i,j) is a WALL_VISCOUS_ISOTHERMAL boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	    //Wu.rho = Wl.p/(Wl.R*Twall); Wu.v.x = ZERO; Wu.v.y = ZERO; Wu.p = Wl.p; 
	    //Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	    Wu = HALF*(Wl+W[i-1][j]); Wu.v = Vector2D_ZERO; Wu.rho = Wu.p/(Wl.R*Twall);
	    Wd = HALF*(Wl+W[i+1][j]); Wd.v = Vector2D_ZERO; Wd.rho = Wd.p/(Wl.R*Twall);
	  } else if (Grid.BCtypeN[i] == BC_MOVING_WALL_HEATFLUX) {
	    // NORTH face of cell (i,j) is a MOVINGWALL_HEATFLUX boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_HEATFLUX;
	    //Wu.rho = Wl.rho; Wu.v.x = Vwall.x; Wu.v.y = Vwall.y; Wu.p = Wl.p; 
	    //Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	    Wu = HALF*(Wl+W[i-1][j]); Wu.v = Vwall;
	    Wd = HALF*(Wl+W[i+1][j]); Wd.v = Vwall;
	  } else if (Grid.BCtypeN[i] == BC_MOVING_WALL_ISOTHERMAL) {
	    // NORTH face of cell (i,j) is a MOVINGWALL_ISOTHERMAL boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	    //Wu.rho = Wl.p/(Wl.R*Twall); Wu.v.x = Vwall.x; Wu.v.y = Vwall.y; 
	    //Wu.p = Wl.p; Wu.k = Wl.k; Wu.omega = Wl.omega;Wu.ke = Wl.ke; Wu.ee = Wl.ee;
	    Wu = HALF*(Wl+W[i-1][j]); Wu.v = Vwall; Wu.rho = Wu.p/(Wl.R*Twall);
	    Wd = HALF*(Wl+W[i+1][j]); Wd.v = Vwall; Wd.rho = Wd.p/(Wl.R*Twall);
	  } else if (Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
	    // NORTH face of cell (i,j) is a BURNING_SURFACE boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	    Wu = BurningSurface(Wl,Grid.nfaceN(i,j));
	    Wd = Wu;
	  } else if (Grid.BCtypeN[i] == BC_MASS_INJECTION) {
	    // NORTH face of cell (i,j) is a MASS_INJECTION boundary.
	    viscous_bc_flag = DIAMONDPATH_LEFT_TRIANGLE_ISOTHERMAL;
	    //Wu = MassInjection(Wl,Grid.nfaceN(i,j));
	    Wu = MassInjection2(Wl,Grid.xfaceN(i,j),Grid.nfaceN(i,j),Twall);
	    Wd = Wu;
	  }
	  switch(IP.i_Viscous_Reconstruction) {
	  case VISCOUS_RECONSTRUCTION_CARTESIAN :
	  case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	    Xu = Grid.Node[i  ][j+1].X;
	    Xd = Grid.Node[i+1][j+1].X;
	    break;
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	  case VISCOUS_RECONSTRUCTION_HYBRID :
	    Xu = Grid.xfaceN(i,j);
	    Xr = Xu; Wr = Wu;
	    dWdxl = dWdx[i][j]; dWdxr = dWdxl;
	    dWdyl = dWdy[i][j]; dWdyr = dWdyl;
	    break;
	  };
	} else {
	  // NORTH face is either a normal cell or possibly a non-viscous
	  // boundary condition.
	  Xl = Grid.Cell[i][j  ].Xc; Wl = W[i][j  ];
	  Xr = Grid.Cell[i][j+1].Xc; Wr = W[i][j+1];
	  switch(IP.i_Viscous_Reconstruction) {
	  case VISCOUS_RECONSTRUCTION_CARTESIAN :
	  case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	    viscous_bc_flag = DIAMONDPATH_QUADRILATERAL_RECONSTRUCTION;
	    Xu = Grid.Node[i  ][j+1].X; Wu = WnNW(i,j);
	    Xd = Grid.Node[i+1][j+1].X; Wd = WnNE(i,j);
	    break;
	  case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	  case VISCOUS_RECONSTRUCTION_HYBRID :
	    Xu = Grid.xfaceN(i,j); Wu = HALF*(Wl + Wr);
	    dWdxl = dWdx[i][j]; dWdxr = dWdy[i][j+1];
	    dWdyl = dWdy[i][j]; dWdyr = dWdy[i][j+1];
	    break;
	  };
	}
	// Compute the NORTH face viscous flux.
	switch(IP.i_Viscous_Reconstruction) {
	case VISCOUS_RECONSTRUCTION_CARTESIAN :
	case VISCOUS_RECONSTRUCTION_DIAMOND_PATH :
	  Flux -= ViscousFluxDiamondPath_n(Grid.xfaceN(i,j),
					   Xl,Wl,Xd,Wd,Xr,Wr,Xu,Wu,
					   Grid.nfaceN(i,j),
					   Axisymmetric,
					   viscous_bc_flag,
					   Wall[i][j].ywall,
					   Wall[i][j].yplus);
	  break;
	case VISCOUS_RECONSTRUCTION_ARITHMETIC_AVERAGE :
	  Flux -= ViscousFlux_n(Xu,Wu,HALF*(dWdxl+dWdxr),
				HALF*(dWdyl+dWdyr),
				Grid.nfaceN(i,j),
				Axisymmetric,OFF,
				Wall[i][j].ywall,
				Wall[i][j].yplus);
	  break;
	case VISCOUS_RECONSTRUCTION_HYBRID :
	  Flux -= ViscousFluxHybrid_n(Xu,Wu,Xl,Wl,dWdxl,dWdyl,
				      Xr,Wr,dWdxr,dWdyr,
				      Grid.nfaceN(i,j),
				      Axisymmetric,
				      Wall[i][j].ywall,
 				      Wall[i][j].yplus);
	  break;
	};
      }

      // Evaluate cell-averaged solution changes.
      dUdt[i][j  ][k_residual] -= (IP.CFL_Number*dt[i][j  ])*
	Flux*Grid.lfaceN(i,j)/Grid.Cell[i][j].A;
      dUdt[i][j+1][k_residual] += (IP.CFL_Number*dt[i][j+1])*
	Flux*Grid.lfaceS(i,j+1)/Grid.Cell[i][j+1].A;

      // Save south and north face boundary flux.
      if (j == JCl-1) {
 	FluxS[i] = -Flux*Grid.lfaceS(i,j+1);
      } else if (j == JCu) {
	FluxN[i] = Flux*Grid.lfaceN(i,j);
      }

    }

    dUdt[i][JCl-1][k_residual] = NavierStokes2D_U_VACUUM;
    dUdt[i][JCu+1][k_residual] = NavierStokes2D_U_VACUUM;

  }

  // Zero the residuals for the turbulence variables according to
  // the turbulence boundary condition.
  error_flag = Turbulence_Zero_Residual(*this,
					i_stage,IP);
  if (error_flag) return error_flag;

  // Residual for the stage successfully calculated.
  return 0;

}
