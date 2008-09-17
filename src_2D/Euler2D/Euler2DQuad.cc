/*! \file Euler2DQuad.cc
  @brief Subroutines for 2D Euler Quadrilateral Mesh Solution Class. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "Euler2DQuad.h"                // Euler2D_Quad_Block class

/**********************************************************************
 * Euler2D_Quad_Block -- Create storage for the static variables.     *
 **********************************************************************/
// Initialize residual_variable
int Euler2D_Quad_Block::residual_variable = 1;
// Initialize Number_of_Residual_Norms
int Euler2D_Quad_Block::Number_of_Residual_Norms = 1;
// Initialize Flow_Type
int Euler2D_Quad_Block::Flow_Type = FLOWTYPE_INVISCID;
// Initialize RefW
Euler2D_pState Euler2D_Quad_Block::RefW(1.0);
// Initialize ExactSoln
Euler2D_ExactSolutions *Euler2D_Quad_Block::ExactSoln = NULL;

/***********************************************************************
 * Euler2D_Quad_Block -- Single Block Member Functions.                *
 **********************************************************************/

/*****************************************************//**
 * Copy the solution information of quadrilateral solution 
 * block SolnBlk to the current solution block.
 ********************************************************/
Euler2D_Quad_Block & Euler2D_Quad_Block::operator =(const Euler2D_Quad_Block &Soln){
  
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

    /* Set the same number of high-order objects
       as that of the rhs block. */
    allocate_HighOrder_Array(Soln.NumberOfHighOrderVariables);

  } else {
    deallocate();
  }

  /* Set the axisymmetric/planar flow indicator. */
  Axisymmetric = Soln.Axisymmetric;
  
  /* Copy the grid. */
  Grid = Soln.Grid;
  
  /* Copy the solution information from Soln. */
  if (Soln.U != NULL) {
    for ( j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
      for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
	U[i][j] = Soln.U[i][j];
	W[i][j] = Soln.W[i][j];
	for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_EULER2D-1 ; ++k ) {
	  dUdt[i][j][k] = Soln.dUdt[i][j][k];
	} /* endfor */
	dWdx[i][j] = Soln.dWdx[i][j];
	dWdy[i][j] = Soln.dWdy[i][j];
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

    // allocate memory for high-order boundary conditions.
    allocate_HighOrder_BoundaryConditions();

    for (j  = JCl ; j <= JCu ; ++j ) {
      // Copy West high-order BCs
      if (HO_WoW != NULL){
	HO_WoW[j] = Soln.HO_WoW[j];
      }

      // Copy East high-order BCs
      if (HO_WoE != NULL){
	HO_WoE[j] = Soln.HO_WoE[j];
      }
    }

    for ( i = ICl ; i <= ICu ; ++i ) {
      // Copy South high-order BCs
      if (HO_WoS != NULL){
	HO_WoS[i] = Soln.HO_WoS[i];
      }
      
      // Copy North high-order BCs
      if (HO_WoN != NULL){
	HO_WoN[i] = Soln.HO_WoN[i];
      }
    }
    
    // Copy the high-order objects
    for (k = 1; k <= NumberOfHighOrderVariables; ++k){
      HighOrderVariable(k-1) = Soln.HighOrderVariable(k-1);
    }/* endfor */
    
  }/* endif */

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
void Euler2D_Quad_Block::Set_Boundary_Reference_States_Based_On_Input(const Euler2D_Input_Parameters &IP){

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
void Euler2D_Quad_Block::Set_Reference_Values_For_Boundary_States(const Euler2D_pState & Ref_North,
								  const Euler2D_pState & Ref_South,
								  const Euler2D_pState & Ref_East,
								  const Euler2D_pState & Ref_West){
  Ref_State_BC_North = Ref_North;
  Ref_State_BC_South = Ref_South;
  Ref_State_BC_East = Ref_East;
  Ref_State_BC_West = Ref_West;
}

/*********************************************************//**
 * This routine evaluates the residual for the specified
 * solution block using a 2nd-order limited upwind      
 * finite-volume spatial discretization scheme with     
 * either the Godunov, Roe, Rusanov, HLLE, Linde, or    
 * HLLC flux functions.                                 
 * The residual is stored in dUdt[][][0].               
 *                                                      
 *************************************************************/
int Euler2D_Quad_Block::dUdt_Residual_Evaluation(const Euler2D_Input_Parameters &IP){

  int i, j;
  Vector2D dX;
  Euler2D_pState Wl, Wr;
  Euler2D_cState Flux;

  /* Perform the linear reconstruction within each cell
     of the computational grid for this stage. */
    
  switch(IP.i_Reconstruction) {
  case RECONSTRUCTION_GREEN_GAUSS :
    Linear_Reconstruction_GreenGauss(*this,
				     IP.i_Limiter);    
    break;
  case RECONSTRUCTION_LEAST_SQUARES :
    Linear_Reconstruction_LeastSquares(*this,
				       IP.i_Limiter);
    break;
  default:
    Linear_Reconstruction_LeastSquares(*this,
				       IP.i_Limiter);
    break;
  } /* endswitch */

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using a second-order
       limited upwind scheme with a variety of flux functions. */
    
    // Add i-direction (zeta-direction) fluxes.
  for ( j  = JCl-1 ; j <= JCu+1 ; ++j ) {
    dUdt[ICl-1][j][0] = Euler2D_U_VACUUM;
          
    for ( i = ICl-1 ; i <= ICu ; ++i ) {

      dUdt[i+1][j][0] = Euler2D_U_VACUUM;
    
      if ( j > JCl-1 && j < JCu+1 ) {
    
	/* Evaluate the cell interface i-direction fluxes. */
	  
	if (i == ICl-1 && 
	    (Grid.BCtypeW[j] == BC_REFLECTION ||
	     Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
	     Grid.BCtypeW[j] == BC_CHARACTERISTIC_VELOCITY ||
	     Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
	     Grid.BCtypeW[j] == BC_MASS_INJECTION ||
	     Grid.BCtypeW[j] == BC_RINGLEB_FLOW)) {
	  dX = Grid.xfaceW(i+1, j)-Grid.Cell[i+1][j].Xc;
	  Wr = W[i+1][j] + 
	    (phi[i+1][j]^dWdx[i+1][j])*dX.x +
	    (phi[i+1][j]^dWdy[i+1][j])*dX.y;
	  if (Grid.BCtypeW[j] == BC_REFLECTION) {
	    Wl = Reflect(Wr, Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
	    Wl = BurningSurface(Wr, Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_MASS_INJECTION) {
	    Wl = MassInjection(Wr,Grid.nfaceW(i+1,j),OFF);
	  } else if (Grid.BCtypeW[j] == BC_CHARACTERISTIC) {
	    Wl = BC_Characteristic_Pressure(Wr, 
					    WoW[j], 
					    Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_CHARACTERISTIC_VELOCITY) {
	    Wl = BC_Characteristic(Wr, 
				   WoW[j], 
				   Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_RINGLEB_FLOW) {
	    Wl = RinglebFlow(Wl,Grid.xfaceW(i+1,j));
	    //Wl = Reflect(Wr, Grid.nfaceW(i+1, j));
	  } /* endif */
	} else if (i == ICu && 
		   (Grid.BCtypeE[j] == BC_REFLECTION ||
		    Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
		    Grid.BCtypeE[j] == BC_CHARACTERISTIC_VELOCITY ||
		    Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
		    Grid.BCtypeE[j] == BC_MASS_INJECTION ||
		    Grid.BCtypeE[j] == BC_RINGLEB_FLOW)) {
	  dX = Grid.xfaceE(i, j)-Grid.Cell[i][j].Xc;
	  Wl = W[i][j] + 
	    (phi[i][j]^dWdx[i][j])*dX.x +
	    (phi[i][j]^dWdy[i][j])*dX.y;
	  if (Grid.BCtypeE[j] == BC_REFLECTION) {
	    Wr = Reflect(Wl, Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_BURNING_SURFACE) {
	    Wr = BurningSurface(Wl, Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_MASS_INJECTION) {
	    Wr = MassInjection(Wl,Grid.nfaceE(i,j),OFF);
	  } else if (Grid.BCtypeE[j] == BC_CHARACTERISTIC) {
	    Wr = BC_Characteristic_Pressure(Wl, 
					    WoE[j], 
					    Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_CHARACTERISTIC_VELOCITY) {
	    Wr = BC_Characteristic(Wl, 
				   WoE[j], 
				   Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_RINGLEB_FLOW) {
	    Wr = RinglebFlow(Wr,Grid.xfaceE(i,j));
	    //Wr = Reflect(Wl, Grid.nfaceE(i, j));
	  } /* endif */
	} else {            
	  dX = Grid.xfaceE(i, j)-Grid.Cell[i][j].Xc;
	  Wl = W[i][j] + 
	    (phi[i][j]^dWdx[i][j])*dX.x +
	    (phi[i][j]^dWdy[i][j])*dX.y;
	  dX = Grid.xfaceW(i+1, j)-Grid.Cell[i+1][j].Xc;
	  Wr = W[i+1][j] + 
	    (phi[i+1][j]^dWdx[i+1][j])*dX.x +
	    (phi[i+1][j]^dWdy[i+1][j])*dX.y;
	} /* endif */

	switch(IP.i_Flux_Function) {
	case FLUX_FUNCTION_GODUNOV :
	  Flux = FluxGodunov_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_ROE :
	  Flux = FluxRoe_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_RUSANOV :
	  Flux = FluxRusanov_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_HLLE :
	  Flux = FluxHLLE_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_LINDE :
	  Flux = FluxLinde_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_HLLC :
	  Flux = FluxHLLC_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_VANLEER :
	  Flux = FluxVanLeer_n(Wl, Wr, Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_AUSM :
	  Flux = FluxAUSM_n(Wl, Wr, Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_AUSMplus :
	  Flux = FluxAUSMplus_n(Wl, Wr, Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_ROE_PRECON_WS :
	  Flux = FluxRoe_n_Precon_WS(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_HLLE_PRECON_WS :
	  Flux = FluxHLLE_n_Precon_WS(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	default:
	  Flux = FluxRoe_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	} /* endswitch */

	  /* Evaluate cell-averaged solution changes. */
	  
	dUdt[i][j][0] -= 
	  Flux*Grid.lfaceE(i, j)/
	  Grid.Cell[i][j].A;
	dUdt[i+1][j][0] += 
	  Flux*Grid.lfaceW(i+1, j)/
	  Grid.Cell[i+1][j].A;

	/* Include axisymmetric source terms as required. */

	if (Axisymmetric) {
	  dUdt[i][j][0] += 
	    S(W[i][j], Grid.Cell[i][j].Xc);
	} /* endif */

	  /* Save west and east face boundary flux. */
	  
	if (i == ICl-1) {
	  FluxW[j] = -Flux*Grid.lfaceW(i+1, j);
	} else if (i == ICu) {
	  FluxE[j] = Flux*Grid.lfaceE(i, j);
	} /* endif */ 

      } /* endif */
    } /* endfor */
      
    if ( j > JCl-1 && j < JCu+1 ) {
      dUdt[ICl-1][j][0] = Euler2D_U_VACUUM;
      dUdt[ICu+1][j][0] = Euler2D_U_VACUUM;
    } /* endif */
  } /* endfor */
    
    // Add j-direction (eta-direction) fluxes.
  for ( i = ICl ; i <= ICu ; ++i ) {
    for ( j  = JCl-1 ; j <= JCu ; ++j ) {
	
      /* Evaluate the cell interface j-direction fluxes. */
	
      if (j == JCl-1 && 
	  (Grid.BCtypeS[i] == BC_REFLECTION ||
	   Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
	   Grid.BCtypeS[i] == BC_CHARACTERISTIC_VELOCITY ||
	   Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
	   Grid.BCtypeS[i] == BC_MASS_INJECTION ||
	   Grid.BCtypeS[i] == BC_RINGLEB_FLOW)) {
	dX = Grid.xfaceS(i, j+1)-Grid.Cell[i][j+1].Xc;
	Wr = W[i][j+1] +
	  (phi[i][j+1]^dWdx[i][j+1])*dX.x +
	  (phi[i][j+1]^dWdy[i][j+1])*dX.y;
	if (Grid.BCtypeS[i] == BC_REFLECTION) {
	  Wl = Reflect(Wr, Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
	  Wl = BurningSurface(Wr, Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_MASS_INJECTION) {
	  Wl = MassInjection(Wr,Grid.nfaceS(i,j+1),OFF);
	} else if (Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
	  Wl = BC_Characteristic_Pressure(Wr, 
					  WoS[i], 
					  Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_CHARACTERISTIC_VELOCITY) {
	  Wl = BC_Characteristic(Wr, 
				 WoS[i], 
				 Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_RINGLEB_FLOW) {
	  //Wl = RinglebFlow(Wl,Grid.xfaceS(i,j+1));
	  Wl = BC_Characteristic_Pressure(Wr,
					  RinglebFlow(Wr,Grid.xfaceS(i,j+1)), 
					  Grid.nfaceS(i, j+1));
	} /* endif */
      } else if (j == JCu && 
		 (Grid.BCtypeN[i] == BC_REFLECTION ||
		  Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
		  Grid.BCtypeN[i] == BC_CHARACTERISTIC_VELOCITY ||
		  Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
		  Grid.BCtypeN[i] == BC_RINGLEB_FLOW)) {
	dX = Grid.xfaceN(i, j)-Grid.Cell[i][j].Xc;
	Wl = W[i][j] + 
	  (phi[i][j]^dWdx[i][j])*dX.x +
	  (phi[i][j]^dWdy[i][j])*dX.y;
	if (Grid.BCtypeN[i] == BC_REFLECTION) {
	  Wr = Reflect(Wl, Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
	  Wr = BurningSurface(Wl, Grid.nfaceN(i, j));
	} else if (Grid.BCtypeS[i] == BC_MASS_INJECTION) {
	  Wr = MassInjection(Wr,Grid.nfaceN(i,j),OFF);
	} else if (Grid.BCtypeN[i] == BC_CHARACTERISTIC) {
	  Wr = BC_Characteristic_Pressure(Wl, 
					  WoN[i], 
					  Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_CHARACTERISTIC_VELOCITY) {
	  Wr = BC_Characteristic(Wl, 
				 WoN[i], 
				 Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_RINGLEB_FLOW) {
	  //Wr = RinglebFlow(Wr,Grid.xfaceN(i,j));
	  Wr = BC_Characteristic_Pressure(Wl, 
					  RinglebFlow(Wr,Grid.xfaceN(i,j)), 
					  Grid.nfaceN(i, j));
	} /* endif */
      } else {
	dX = Grid.xfaceN(i, j)-Grid.Cell[i][j].Xc;
	Wl = W[i][j] + 
	  (phi[i][j]^dWdx[i][j])*dX.x +
	  (phi[i][j]^dWdy[i][j])*dX.y;
	dX = Grid.xfaceS(i, j+1)-Grid.Cell[i][j+1].Xc;
	Wr = W[i][j+1] +
	  (phi[i][j+1]^dWdx[i][j+1])*dX.x +
	  (phi[i][j+1]^dWdy[i][j+1])*dX.y;
      } /* endif */
	
      switch(IP.i_Flux_Function) {
      case FLUX_FUNCTION_GODUNOV :
	Flux = FluxGodunov_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_ROE :
	Flux = FluxRoe_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_RUSANOV :
	Flux = FluxRusanov_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_HLLE :
	Flux = FluxHLLE_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_LINDE :
	Flux = FluxLinde_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_HLLC :
	Flux = FluxHLLC_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_VANLEER :
	Flux = FluxVanLeer_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_AUSM :
	Flux = FluxAUSM_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_AUSMplus :
	Flux = FluxAUSMplus_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_ROE_PRECON_WS :
	Flux = FluxRoe_n_Precon_WS(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_HLLE_PRECON_WS :
	Flux = FluxHLLE_n_Precon_WS(Wl, Wr, Grid.nfaceN(i, j));
	break;
      default:
	Flux = FluxRoe_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      } /* endswitch */
	
      /* Evaluate cell-averaged solution changes. */
	
      dUdt[i][j][0] -= 
	Flux*Grid.lfaceN(i, j)/
	Grid.Cell[i][j].A;
      dUdt[i][j+1][0] += 
	Flux*Grid.lfaceS(i, j+1)/
	Grid.Cell[i][j+1].A;

      /* Save south and north face boundary flux. */
	
      if (j == JCl-1) {
	FluxS[i] = -Flux*Grid.lfaceS(i, j+1);
      } else if (j == JCu) {
	FluxN[i] = Flux*Grid.lfaceN(i, j);
      } /* endif */
	
    } /* endfor */
      
    dUdt[i][JCl-1][0] = Euler2D_U_VACUUM;
    dUdt[i][JCu+1][0] = Euler2D_U_VACUUM;
  } /* endfor */
    
    /* Residual successfully evaluated. */
  return 0;

}

/*********************************************************//**
 * This routine determines the solution residuals for a 
 * given stage of a variety of multi-stage explicit     
 * time integration schemes for a given solution block. 
 *                                                      
 ************************************************************/
int Euler2D_Quad_Block::dUdt_Multistage_Explicit(const int &i_stage,
						 const Euler2D_Input_Parameters &IP){
  
  int i, j, k_residual;
  double omega;
  Vector2D dX;
  Euler2D_pState Wl, Wr;
  Euler2D_cState Flux;

  /* Evaluate the solution residual for stage 
     i_stage of an N stage scheme. */

  /* Evaluate the time step fraction and residual storage location for the stage. */
    
  switch(IP.i_Time_Integration) {
  case TIME_STEPPING_EXPLICIT_EULER :
    omega = Runge_Kutta(i_stage, IP.N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
    omega = Runge_Kutta(i_stage, IP.N_Stage);
    k_residual = 0;
    break;
  case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
    omega = Runge_Kutta(i_stage, IP.N_Stage);
    k_residual = 0;
    if (IP.N_Stage == 4) {
      if (i_stage == 4) {
	k_residual = 0;
      } else {
	k_residual = i_stage - 1;
      } /* endif */
    } /* endif */
    break;
  case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
    omega = MultiStage_Optimally_Smoothing(i_stage, 
					   IP.N_Stage,
					   IP.i_Limiter);
    k_residual = 0;
    break;
  default:
    omega = Runge_Kutta(i_stage, IP.N_Stage);
    k_residual = 0;
    break;
  } /* endswitch */
    
    /* Perform the linear reconstruction within each cell
       of the computational grid for this stage. */
    
  switch(IP.i_Reconstruction) {
  case RECONSTRUCTION_GREEN_GAUSS :
    Linear_Reconstruction_GreenGauss(*this,
				     IP.i_Limiter);    
    break;
  case RECONSTRUCTION_LEAST_SQUARES :
    Linear_Reconstruction_LeastSquares(*this,
				       IP.i_Limiter);
    break;
  default:
    Linear_Reconstruction_LeastSquares(*this,
				       IP.i_Limiter);
    break;
  } /* endswitch */

    /* Evaluate the time rate of change of the solution
       (i.e., the solution residuals) using a second-order
       limited upwind scheme with a variety of flux functions. */
    
    // Add i-direction (zeta-direction) fluxes.
  for ( j  = JCl-1 ; j <= JCu+1 ; ++j ) {
    if ( i_stage == 1 ) {
      Uo[ICl-1][j] = U[ICl-1][j];
      dUdt[ICl-1][j][k_residual] = Euler2D_U_VACUUM;
    } else {
      dUdt[ICl-1][j][k_residual] = Euler2D_U_VACUUM;
    } /* endif */
    
    for ( i = ICl-1 ; i <= ICu ; ++i ) {
      if ( i_stage == 1 ) {
	Uo[i+1][j] = U[i+1][j];
	dUdt[i+1][j][k_residual] = Euler2D_U_VACUUM;
      } else if ( j > JCl-1 && j < JCu+1 ) {
	switch(IP.i_Time_Integration) {
	case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
	  //dUdt[i+1][j][k_residual] = 
	  //   dUdt[i+1][j][k_residual];
	  break;
	case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
	  if (IP.N_Stage == 2) {
	    //dUdt[i+1][j][k_residual] = 
	    //   dUdt[i+1][j][k_residual];
	  } else if (IP.N_Stage == 4 && i_stage == 4) {
	    dUdt[i+1][j][k_residual] = 
	      dUdt[i+1][j][0] + 
	      TWO*dUdt[i+1][j][1] +
	      TWO*dUdt[i+1][j][2];
	  } else {
	    dUdt[i+1][j][k_residual] = Euler2D_U_VACUUM;
	  } /* endif */
	  break;
	case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
	  dUdt[i+1][j][k_residual] = Euler2D_U_VACUUM;
	  break;
	default:
	  dUdt[i+1][j][k_residual] = Euler2D_U_VACUUM;
	  break;
	} /* endswitch */
      } /* endif */
    
      if ( j > JCl-1 && j < JCu+1 ) {
    
	/* Evaluate the cell interface i-direction fluxes. */
    
	if (i == ICl-1 && 
	    (Grid.BCtypeW[j] == BC_REFLECTION ||
	     Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
	     Grid.BCtypeW[j] == BC_CHARACTERISTIC_VELOCITY ||
	     Grid.BCtypeW[j] == BC_BURNING_SURFACE ||
	     Grid.BCtypeW[j] == BC_MASS_INJECTION ||
	     Grid.BCtypeW[j] == BC_RINGLEB_FLOW)) {
	  dX = Grid.xfaceW(i+1, j)-Grid.Cell[i+1][j].Xc;
	  Wr = W[i+1][j] + 
	    (phi[i+1][j]^dWdx[i+1][j])*dX.x +
	    (phi[i+1][j]^dWdy[i+1][j])*dX.y;
	  if (Grid.BCtypeW[j] == BC_REFLECTION) {
	    Wl = Reflect(Wr, Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_BURNING_SURFACE) {
	    Wl = BurningSurface(Wr, Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_MASS_INJECTION) {
	    Wl = MassInjection(Wr,Grid.nfaceW(i+1,j),OFF);
	  } else if (Grid.BCtypeW[j] == BC_CHARACTERISTIC) {
	    Wl = BC_Characteristic_Pressure(Wr, 
					    WoW[j], 
					    Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_CHARACTERISTIC_VELOCITY) {
	    Wl = BC_Characteristic(Wr, 
				   WoW[j], 
				   Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_RINGLEB_FLOW) {
	    Wl = RinglebFlow(Wl,Grid.xfaceW(i+1,j));
	    //Wl = Reflect(Wr, Grid.nfaceW(i+1, j));
	  } /* endif */
	} else if (i == ICu && 
		   (Grid.BCtypeE[j] == BC_REFLECTION ||
		    Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
		    Grid.BCtypeE[j] == BC_CHARACTERISTIC_VELOCITY ||
		    Grid.BCtypeE[j] == BC_BURNING_SURFACE ||
		    Grid.BCtypeE[j] == BC_MASS_INJECTION ||
		    Grid.BCtypeE[j] == BC_RINGLEB_FLOW)) {
	  dX = Grid.xfaceE(i, j)-Grid.Cell[i][j].Xc;
	  Wl = W[i][j] + 
	    (phi[i][j]^dWdx[i][j])*dX.x +
	    (phi[i][j]^dWdy[i][j])*dX.y;
	  if (Grid.BCtypeE[j] == BC_REFLECTION) {
	    Wr = Reflect(Wl, Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_BURNING_SURFACE) {
	    Wr = BurningSurface(Wl, Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_MASS_INJECTION) {
	    Wr = MassInjection(Wl,Grid.nfaceE(i,j),OFF);
	  } else if (Grid.BCtypeE[j] == BC_CHARACTERISTIC) {
	    Wr = BC_Characteristic_Pressure(Wl, 
					    WoE[j], 
					    Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_CHARACTERISTIC_VELOCITY) {
	    Wr = BC_Characteristic(Wl, 
				   WoE[j], 
				   Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_RINGLEB_FLOW) {
	    Wr = RinglebFlow(Wr,Grid.xfaceE(i,j));
	    //Wr = Reflect(Wl, Grid.nfaceE(i, j));
	  } /* endif */
	} else {            
	  dX = Grid.xfaceE(i, j)-Grid.Cell[i][j].Xc;
	  Wl = W[i][j] + 
	    (phi[i][j]^dWdx[i][j])*dX.x +
	    (phi[i][j]^dWdy[i][j])*dX.y;
	  dX = Grid.xfaceW(i+1, j)-Grid.Cell[i+1][j].Xc;
	  Wr = W[i+1][j] + 
	    (phi[i+1][j]^dWdx[i+1][j])*dX.x +
	    (phi[i+1][j]^dWdy[i+1][j])*dX.y;
	} /* endif */

	switch(IP.i_Flux_Function) {
	case FLUX_FUNCTION_GODUNOV :
	  Flux = FluxGodunov_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_ROE :
	  Flux = FluxRoe_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_RUSANOV :
	  Flux = FluxRusanov_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_HLLE :
	  Flux = FluxHLLE_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_LINDE :
	  Flux = FluxLinde_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_HLLC :
	  Flux = FluxHLLC_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_VANLEER :
	  Flux = FluxVanLeer_n(Wl, Wr, Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_AUSM :
	  Flux = FluxAUSM_n(Wl, Wr, Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_AUSMplus :
	  Flux = FluxAUSMplus_n(Wl, Wr, Grid.nfaceE(i,j));
	  break;
	case FLUX_FUNCTION_ROE_PRECON_WS :
	  Flux = FluxRoe_n_Precon_WS(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	case FLUX_FUNCTION_HLLE_PRECON_WS :
	  Flux = FluxHLLE_n_Precon_WS(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	default:
	  Flux = FluxRoe_n(Wl, Wr, Grid.nfaceE(i, j));
	  break;
	} /* endswitch */
    
	/* Evaluate cell-averaged solution changes. */
    
	dUdt[i][j][k_residual] -= 
	  (IP.CFL_Number*dt[i][j])*
	  Flux*Grid.lfaceE(i, j)/
	  Grid.Cell[i][j].A;
	dUdt[i+1][j][k_residual] += 
	  (IP.CFL_Number*dt[i+1][j])*
	  Flux*Grid.lfaceW(i+1, j)/
	  Grid.Cell[i+1][j].A;

	/* Include axisymmetric source terms as required. */

	if (Axisymmetric) {
	  dUdt[i][j][k_residual] += 
	    (IP.CFL_Number*dt[i][j])*
	    S(W[i][j], Grid.Cell[i][j].Xc);
	} /* endif */

	/* Save west and east face boundary flux. */

	if (i == ICl-1) {
	  FluxW[j] = -Flux*Grid.lfaceW(i+1, j);
	} else if (i == ICu) {
	  FluxE[j] = Flux*Grid.lfaceE(i, j);
	} /* endif */ 

      } /* endif */
    } /* endfor */
    
    if ( j > JCl-1 && j < JCu+1 ) {
      dUdt[ICl-1][j][k_residual] = Euler2D_U_VACUUM;
      dUdt[ICu+1][j][k_residual] = Euler2D_U_VACUUM;
    } /* endif */
  } /* endfor */
    
    // Add j-direction (eta-direction) fluxes.
  for ( i = ICl ; i <= ICu ; ++i ) {
    for ( j  = JCl-1 ; j <= JCu ; ++j ) {
    
      /* Evaluate the cell interface j-direction fluxes. */
         
      if (j == JCl-1 && 
	  (Grid.BCtypeS[i] == BC_REFLECTION ||
	   Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
	   Grid.BCtypeS[i] == BC_CHARACTERISTIC_VELOCITY ||
	   Grid.BCtypeS[i] == BC_BURNING_SURFACE ||
	   Grid.BCtypeS[i] == BC_MASS_INJECTION ||
	   Grid.BCtypeS[i] == BC_RINGLEB_FLOW)) {
	dX = Grid.xfaceS(i, j+1)-Grid.Cell[i][j+1].Xc;
	Wr = W[i][j+1] +
	  (phi[i][j+1]^dWdx[i][j+1])*dX.x +
	  (phi[i][j+1]^dWdy[i][j+1])*dX.y;
	if (Grid.BCtypeS[i] == BC_REFLECTION) {
	  Wl = Reflect(Wr, Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_BURNING_SURFACE) {
	  Wl = BurningSurface(Wr, Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_MASS_INJECTION) {
	  Wl = MassInjection(Wr,Grid.nfaceS(i,j+1),OFF);
	} else if (Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
	  Wl = BC_Characteristic_Pressure(Wr, 
					  WoS[i], 
					  Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_CHARACTERISTIC_VELOCITY) {
	  Wl = BC_Characteristic(Wr, 
				 WoS[i], 
				 Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_RINGLEB_FLOW) {
	  //Wl = RinglebFlow(Grid.xfaceS(i,j+1));
	  Wl = BC_Characteristic_Pressure(Wr,
					  RinglebFlow(Wl,Grid.xfaceS(i,j+1)), 
					  Grid.nfaceS(i, j+1));
	} /* endif */
      } else if (j == JCu && 
		 (Grid.BCtypeN[i] == BC_REFLECTION ||
		  Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
		  Grid.BCtypeN[i] == BC_CHARACTERISTIC_VELOCITY ||
		  Grid.BCtypeN[i] == BC_BURNING_SURFACE ||
		  Grid.BCtypeN[i] == BC_MASS_INJECTION ||
		  Grid.BCtypeN[i] == BC_RINGLEB_FLOW)) {
	dX = Grid.xfaceN(i, j)-Grid.Cell[i][j].Xc;
	Wl = W[i][j] + 
	  (phi[i][j]^dWdx[i][j])*dX.x +
	  (phi[i][j]^dWdy[i][j])*dX.y;
	if (Grid.BCtypeN[i] == BC_REFLECTION) {
	  Wr = Reflect(Wl, Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_BURNING_SURFACE) {
	  Wr = BurningSurface(Wl, Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_MASS_INJECTION) {
	  Wr = MassInjection(Wl,Grid.nfaceN(i,j),OFF);
	} else if (Grid.BCtypeN[i] == BC_CHARACTERISTIC) {
	  Wr = BC_Characteristic_Pressure(Wl, 
					  WoN[i], 
					  Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_CHARACTERISTIC_VELOCITY) {
	  Wr = BC_Characteristic(Wl, 
				 WoN[i], 
				 Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_RINGLEB_FLOW) {
	  //Wr = RinglebFlow(Grid.xfaceN(i,j));
	  Wr = BC_Characteristic_Pressure(Wl, 
					  RinglebFlow(Wr,Grid.xfaceN(i,j)), 
					  Grid.nfaceN(i, j));
	} /* endif */
      } else {
	dX = Grid.xfaceN(i, j)-Grid.Cell[i][j].Xc;
	Wl = W[i][j] + 
	  (phi[i][j]^dWdx[i][j])*dX.x +
	  (phi[i][j]^dWdy[i][j])*dX.y;
	dX = Grid.xfaceS(i, j+1)-Grid.Cell[i][j+1].Xc;
	Wr = W[i][j+1] +
	  (phi[i][j+1]^dWdx[i][j+1])*dX.x +
	  (phi[i][j+1]^dWdy[i][j+1])*dX.y;
      } /* endif */

      switch(IP.i_Flux_Function) {
      case FLUX_FUNCTION_GODUNOV :
	Flux = FluxGodunov_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_ROE :
	Flux = FluxRoe_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_RUSANOV :
	Flux = FluxRusanov_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_HLLE :
	Flux = FluxHLLE_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_LINDE :
	Flux = FluxLinde_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_HLLC :
	Flux = FluxHLLC_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_VANLEER :
	Flux = FluxVanLeer_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_AUSM :
	Flux = FluxAUSM_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_AUSMplus :
	Flux = FluxAUSMplus_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_ROE_PRECON_WS :
	Flux = FluxRoe_n_Precon_WS(Wl, Wr, Grid.nfaceN(i, j));
	break;
      case FLUX_FUNCTION_HLLE_PRECON_WS :
	Flux = FluxHLLE_n_Precon_WS(Wl, Wr, Grid.nfaceN(i, j));
	break;
      default:
	Flux = FluxRoe_n(Wl, Wr, Grid.nfaceN(i, j));
	break;
      } /* endswitch */
    
      /* Evaluate cell-averaged solution changes. */
    
      dUdt[i][j][k_residual] -= 
	(IP.CFL_Number*dt[i][j])*
	Flux*Grid.lfaceN(i, j)/
	Grid.Cell[i][j].A;
      dUdt[i][j+1][k_residual] += 
	(IP.CFL_Number*dt[i][j+1])*
	Flux*Grid.lfaceS(i, j+1)/
	Grid.Cell[i][j+1].A;

      /* Save south and north face boundary flux. */

      if (j == JCl-1) {
	FluxS[i] = -Flux*Grid.lfaceS(i, j+1);
      } else if (j == JCu) {
	FluxN[i] = Flux*Grid.lfaceN(i, j);
      } /* endif */

    } /* endfor */
    
    dUdt[i][JCl-1][k_residual] = Euler2D_U_VACUUM;
    dUdt[i][JCu+1][k_residual] = Euler2D_U_VACUUM;
  } /* endfor */
    
    /* Residual for the stage successfully calculated. */

  return (0);

}
