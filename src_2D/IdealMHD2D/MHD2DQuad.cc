/*! \file MHD2DQuad.cc
  @brief Subroutines for 2D MHD Quadrilateral Mesh Solution Class. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "MHD2DQuad.h"                // MHD2D_Quad_Block class

/********************************************************************
 * MHD2D_Quad_Block -- Create storage for the static variables.     *
 *******************************************************************/
// Initialize residual_variable
int MHD2D_Quad_Block::residual_variable = 1;
// Initialize Number_of_Residual_Norms
int MHD2D_Quad_Block::Number_of_Residual_Norms = 1;
// Initialize Flow_Type
int MHD2D_Quad_Block::Flow_Type = FLOWTYPE_INVISCID;
// Initialize RefW
MHD3D_pState MHD2D_Quad_Block::RefW(1.0);
// Initialize ExactSoln
MHD2D_ExactSolutions *MHD2D_Quad_Block::ExactSoln = NULL;

/*********************************************************************
 * MHD2D_Quad_Block -- Single Block Member Functions.                *
 ********************************************************************/

/*!
 * Copy the solution information of quadrilateral solution 
 * block SolnBlk to the current solution block.
 */
MHD2D_Quad_Block & MHD2D_Quad_Block::operator =(const MHD2D_Quad_Block &Soln){
  
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
	for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_MHD2D-1 ; ++k ) {
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

/*!
 * Assigns boundary condition reference data based on
 * the boundary type for the specified quadrilateral 
 * solution block based on the input parameters.
 * This routine makes the link between user's specifications
 * and the values that are set as boundary reference states.
 */
void MHD2D_Quad_Block::Set_Boundary_Reference_States_Based_On_Input(const MHD2D_Input_Parameters &IP){

  // Set the reference values for boundary reference states
  Set_Reference_Values_For_Boundary_States(IP.Ref_State_BC_North,
					   IP.Ref_State_BC_South,
					   IP.Ref_State_BC_East,
					   IP.Ref_State_BC_West);
}

/*!
 * Assigns reference values for the boundary condition 
 * reference states to what the user specified.
 * These values are used in Set_Boundary_Reference_States()
 * routine.
 */
void MHD2D_Quad_Block::Set_Reference_Values_For_Boundary_States(const MHD3D_pState & Ref_North,
								const MHD3D_pState & Ref_South,
								const MHD3D_pState & Ref_East,
								const MHD3D_pState & Ref_West){
  Ref_State_BC_North = Ref_North;
  Ref_State_BC_South = Ref_South;
  Ref_State_BC_East = Ref_East;
  Ref_State_BC_West = Ref_West;
}

/*!
 * Performs the reconstruction of a limited piecewise  
 * linear solution state within a given cell (i,j) of  
 * the computational mesh for the specified            
 * quadrilateral solution block.  A Green-Gauss        
 * approach is used in the evaluation of the unlimited 
 * solution gradients.  Several slope limiters may be  
 * used.                                               
 */
void MHD2D_Quad_Block::Linear_Reconstruction_GreenGauss(const int &i, const int &j, 
							const int & LimiterType){

  int n, n2, n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4], phi_limiter;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  double l_north, l_south, l_east, l_west;
  Vector2D n_north, n_south, n_east, n_west, dX;
  MHD3D_pState W_nw, W_ne, W_sw, W_se, W_face, DU, DUDx_ave, DUDy_ave;

  /* Carry out the limited solution reconstruction in
     the specified cell of the computational mesh. */

  // === Set reconstruction stencil === 
  SetPiecewiseLinearReconstructionStencil(i,j,
					  i_index, j_index,
					  n_pts);


  // Perform reconstruction.
  if (n_pts > 0) {
    // If 8 neighbours are used, apply Green-Gauss reconstruction
    if (n_pts == 8) {
      W_nw = WnNW(i, j);
      W_ne = WnNE(i, j);
      W_sw = WnSW(i, j);
      W_se = WnSE(i, j);

      l_north = Grid.lfaceN(i, j);
      l_south = Grid.lfaceS(i, j);
      l_east = Grid.lfaceE(i, j);
      l_west = Grid.lfaceW(i, j);

      n_north = Grid.nfaceN(i, j);
      n_south = Grid.nfaceS(i, j);
      n_east = Grid.nfaceE(i, j);
      n_west = Grid.nfaceW(i, j);

      W_face = HALF*(W_nw+W_ne)*l_north; 
      dWdx[i][j] = W_face*n_north.x;
      dWdy[i][j] = W_face*n_north.y;

      W_face = HALF*(W_sw+W_se)*l_south; 
      dWdx[i][j] += W_face*n_south.x;
      dWdy[i][j] += W_face*n_south.y;

      W_face = HALF*(W_ne+W_se)*l_east; 
      dWdx[i][j] += W_face*n_east.x;
      dWdy[i][j] += W_face*n_east.y;

      W_face = HALF*(W_nw+W_sw)*l_west; 
      dWdx[i][j] += W_face*n_west.x;
      dWdy[i][j] += W_face*n_west.y;

      dWdx[i][j] = dWdx[i][j]/Grid.Cell[i][j].A;
      dWdy[i][j] = dWdy[i][j]/Grid.Cell[i][j].A;

      // If <8 neighbours are used, apply least-squares reconstruction
    } else {
      DUDx_ave.zero_all();
      DUDy_ave.zero_all();
      DxDx_ave = ZERO;
      DxDy_ave = ZERO;
      DyDy_ave = ZERO;
    
      for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
	dX = Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - Grid.Cell[i][j].Xc;
	DU = W[ i_index[n2] ][ j_index[n2] ] - W[i][j];
	DUDx_ave += DU*dX.x;
	DUDy_ave += DU*dX.y;
	DxDx_ave += dX.x*dX.x;
	DxDy_ave += dX.x*dX.y;
	DyDy_ave += dX.y*dX.y;
      } /* endfor */
    					    
      DUDx_ave = DUDx_ave/double(n_pts);
      DUDy_ave = DUDy_ave/double(n_pts);
      DxDx_ave = DxDx_ave/double(n_pts);
      DxDy_ave = DxDy_ave/double(n_pts);
      DyDy_ave = DyDy_ave/double(n_pts);
      dWdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
	(DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
      dWdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
	(DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    } /* endif */

    // Calculate slope limiters.    
    if (!Freeze_Limiter) {
      for ( n = 1 ; n <= NUM_VAR_MHD3D ; ++n ) {
	u0Min = W[i][j][n];
	u0Max = u0Min;
	for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
	  u0Min = min(u0Min, W[ i_index[n2] ][ j_index[n2] ][n]);
	  u0Max = max(u0Max, W[ i_index[n2] ][ j_index[n2] ][n]);
	} /* endfor */
    
	dX = Grid.xfaceE(i, j)-Grid.Cell[i][j].Xc;
	uQuad[0] = W[i][j][n] + dWdx[i][j][n]*dX.x + dWdy[i][j][n]*dX.y ;
	dX = Grid.xfaceW(i, j)-Grid.Cell[i][j].Xc;
	uQuad[1] = W[i][j][n] + dWdx[i][j][n]*dX.x + dWdy[i][j][n]*dX.y ;
	dX = Grid.xfaceN(i, j)-Grid.Cell[i][j].Xc;
	uQuad[2] = W[i][j][n] + dWdx[i][j][n]*dX.x + dWdy[i][j][n]*dX.y ;
	dX = Grid.xfaceS(i, j)-Grid.Cell[i][j].Xc;
	uQuad[3] = W[i][j][n] + dWdx[i][j][n]*dX.x + dWdy[i][j][n]*dX.y ;
    
	switch(LimiterType) {
	case LIMITER_ONE :
	  phi_limiter = ONE;
	  break;
	case LIMITER_ZERO :
	  phi_limiter = ZERO;
	  break;
	case LIMITER_BARTH_JESPERSEN :
	  phi_limiter = Limiter_BarthJespersen(uQuad, W[i][j][n], 
					       u0Min, u0Max, 4);
	  break;
	case LIMITER_VENKATAKRISHNAN :
	  phi_limiter = Limiter_Venkatakrishnan(uQuad, W[i][j][n], 
						u0Min, u0Max, 4);
	  break;
	case LIMITER_VANLEER :
	  phi_limiter = Limiter_VanLeer(uQuad, W[i][j][n], 
					u0Min, u0Max, 4);
	  break;
	case LIMITER_VANALBADA :
	  phi_limiter = Limiter_VanAlbada(uQuad, W[i][j][n], 
					  u0Min, u0Max, 4);
	  break;
	default:
	  phi_limiter = Limiter_BarthJespersen(uQuad, W[i][j][n], 
					       u0Min, u0Max, 4);
	  break;
	} /* endswitch */
	
	phi[i][j][n] = phi_limiter;
	
      } /* endfor */
    } /* endif */
  } else {
    dWdx[i][j].zero_all();
    dWdy[i][j].zero_all(); 
    phi[i][j].zero_all();
  } /* endif */
  
}

/*!
 * Performs the reconstruction of a limited piecewise  
 * linear solution state within a given cell (i,j) of  
 * the computational mesh for the specified            
 * quadrilateral solution block.  A least squares      
 * approach is used in the evaluation of the unlimited 
 * solution gradients.  Several slope limiters may be  
 * used.                                               
 */
void MHD2D_Quad_Block::Linear_Reconstruction_LeastSquares(const int &i, const int &j,
							  const int & LimiterType){

  int n, n2, n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4], phi_limiter;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  MHD3D_pState DU, DUDx_ave, DUDy_ave;

  /* Carry out the limited solution reconstruction in
     each cell of the computational mesh. */

  // === Set reconstruction stencil === 
  SetPiecewiseLinearReconstructionStencil(i,j,
					  i_index, j_index,
					  n_pts);

  // Perform reconstruction
  if (n_pts > 0) {
    DUDx_ave.zero_all();
    DUDy_ave.zero_all();
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;
    
    for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
      dX = Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - Grid.Cell[i][j].Xc;
      DU = W[ i_index[n2] ][ j_index[n2] ] - W[i][j];
      DUDx_ave += DU*dX.x;
      DUDy_ave += DU*dX.y;
      DxDx_ave += dX.x*dX.x;
      DxDy_ave += dX.x*dX.y;
      DyDy_ave += dX.y*dX.y;
    } /* endfor */
    					    
    DUDx_ave = DUDx_ave/double(n_pts);
    DUDy_ave = DUDy_ave/double(n_pts);
    DxDx_ave = DxDx_ave/double(n_pts);
    DxDy_ave = DxDy_ave/double(n_pts);
    DyDy_ave = DyDy_ave/double(n_pts);
    dWdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    dWdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    
    // Calculate slope limiters. 
    if (!Freeze_Limiter) {
      for ( n = 1 ; n <= NUM_VAR_MHD3D ; ++n ) {
	u0Min = W[i][j][n];
	u0Max = u0Min;
	for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
	  u0Min = min(u0Min, W[ i_index[n2] ][ j_index[n2] ][n]);
	  u0Max = max(u0Max, W[ i_index[n2] ][ j_index[n2] ][n]);
	} /* endfor */
    
	dX = Grid.xfaceE(i, j)-Grid.Cell[i][j].Xc;
	uQuad[0] = W[i][j][n] + dWdx[i][j][n]*dX.x + dWdy[i][j][n]*dX.y ;
	dX = Grid.xfaceW(i, j)-Grid.Cell[i][j].Xc;
	uQuad[1] = W[i][j][n] + dWdx[i][j][n]*dX.x + dWdy[i][j][n]*dX.y ;
	dX = Grid.xfaceN(i, j)-Grid.Cell[i][j].Xc;
	uQuad[2] = W[i][j][n] + dWdx[i][j][n]*dX.x + dWdy[i][j][n]*dX.y ;
	dX = Grid.xfaceS(i, j)-Grid.Cell[i][j].Xc;
	uQuad[3] = W[i][j][n] + dWdx[i][j][n]*dX.x + dWdy[i][j][n]*dX.y ;
    
	switch(LimiterType) {
	case LIMITER_ONE :
	  phi_limiter = ONE;
	  break;
	case LIMITER_ZERO :
	  phi_limiter = ZERO;
	  break;
	case LIMITER_BARTH_JESPERSEN :
	  phi_limiter = Limiter_BarthJespersen(uQuad, W[i][j][n], 
					       u0Min, u0Max, 4);
	  break;
	case LIMITER_VENKATAKRISHNAN :
	  phi_limiter = Limiter_Venkatakrishnan(uQuad, W[i][j][n], 
						u0Min, u0Max, 4);
	  break;
	case LIMITER_VANLEER :
	  phi_limiter = Limiter_VanLeer(uQuad, W[i][j][n], 
					u0Min, u0Max, 4);
	  break;
	case LIMITER_VANALBADA :
	  phi_limiter = Limiter_VanAlbada(uQuad, W[i][j][n], 
					  u0Min, u0Max, 4);
	  break;
	default:
	  phi_limiter = Limiter_BarthJespersen(uQuad, W[i][j][n], 
					       u0Min, u0Max, 4);
	  break;
	} /* endswitch */

	phi[i][j][n] = phi_limiter;

      } /* endfor */
    } /* endif */
  } else {
    dWdx[i][j].zero_all();
    dWdy[i][j].zero_all(); 
    phi[i][j].zero_all();
  } /* endif */

}

/*!
 * This routine evaluates the residual for the specified
 * solution block using a 2nd-order limited upwind      
 * finite-volume spatial discretization scheme with     
 * different flux functions.
 * The residual is stored in dUdt[][][0].               
 */
int MHD2D_Quad_Block::dUdt_Residual_Evaluation(const MHD2D_Input_Parameters &IP){

  int i, j;
  Vector2D dX;
  MHD3D_pState Wl, Wr, Wi;
  MHD3D_cState Flux;

  /* Perform the linear reconstruction within each cell
     of the computational grid for this stage. */
  ComputeLinearSolutionReconstruction(IP.i_Reconstruction,
				      IP.i_Limiter);

  /* Evaluate the time rate of change of the solution
     (i.e., the solution residuals) using a second-order
     limited upwind scheme with a variety of flux functions. */
    
    // Add i-direction (zeta-direction) fluxes.
  for ( j  = JCl-1 ; j <= JCu+1 ; ++j ) {
    dUdt[ICl-1][j][0].zero_all();
          
    for ( i = ICl-1 ; i <= ICu ; ++i ) {
      dUdt[i+1][j][0].zero_all();
    
      if ( j > JCl-1 && j < JCu+1 ) {
    
	/* Evaluate the cell interface i-direction fluxes. */
	  
	if (i == ICl-1 && 
	    (Grid.BCtypeW[j] == BC_REFLECTION ||
	     Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
	     Grid.BCtypeW[j] == BC_CHARACTERISTIC_VELOCITY ||
	     Grid.BCtypeW[j] == BC_FROZEN ||
	     Grid.BCtypeW[j] == BC_EXACT_SOLUTION ||
	     Grid.BCtypeW[j] == BC_WALL_INVISCID)) {
	  // Determine right state
	  Wr = PiecewiseLinearSolutionAtLocation(i+1, j,
						 Grid.xfaceW(i+1, j));
	  // Determine left state
	  if (Grid.BCtypeW[j] == BC_REFLECTION) {
	    Wl = Reflect(Wr, Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_CHARACTERISTIC) {
	    Wl = BC_Characteristic_Pressure(Wr, 
					    WoW[j], 
					    Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_CHARACTERISTIC_VELOCITY) {
	    Wl = BC_Characteristic(Wr, 
				   WoW[j], 
				   Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_FROZEN) {
	    // calculate Wl based on the ghost cell reconstruction
	    Wl = PiecewiseLinearSolutionAtLocation(i,j,
						   Grid.xfaceW(i+1,j));
	  } else if (Grid.BCtypeW[j] == BC_EXACT_SOLUTION) {
	    // calculate Wl using the exact solution at the interface and avoiding the problem to become ill-posed.
	    if (ExactSoln->IsExactSolutionSet()){
	      Wl = BC_Characteristic_Pressure(Wr,
					      ExactSoln->Solution(Grid.xfaceW(i+1,j).x,Grid.xfaceW(i+1,j).y),
					      Grid.nfaceW(i+1, j));
	    } else {
	      throw runtime_error("MHD2D_Quad_Block::dUdt_Residual_Evaluation() ERROR! There is no exact solution set for the Exact_Solution BC.");
	    }
	  } else if (Grid.BCtypeW[j] == BC_WALL_INVISCID) {
	    Wl = Reflect(Wr, Grid.nfaceW(i+1, j));
	    Wi = PiecewiseLinearSolutionAtLocation(i,j,
						  Grid.xfaceW(i+1,j));
	    Wl.d() = Wi.d();
	    Wl.p() = Wi.p();
	  } /* endif */
	} else if (i == ICu && 
		   (Grid.BCtypeE[j] == BC_REFLECTION ||
		    Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
		    Grid.BCtypeE[j] == BC_CHARACTERISTIC_VELOCITY ||
		    Grid.BCtypeE[j] == BC_FROZEN ||
		    Grid.BCtypeE[j] == BC_EXACT_SOLUTION ||
		    Grid.BCtypeE[j] == BC_WALL_INVISCID)) {
	  // Determine left state
	  Wl = PiecewiseLinearSolutionAtLocation(i,j,
						 Grid.xfaceE(i, j));
	  // Determine right state
	  if (Grid.BCtypeE[j] == BC_REFLECTION) {
	    Wr = Reflect(Wl, Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_CHARACTERISTIC) {
	    Wr = BC_Characteristic_Pressure(Wl, 
					    WoE[j], 
					    Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_CHARACTERISTIC_VELOCITY) {
	    Wr = BC_Characteristic(Wl, 
				   WoE[j], 
				   Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_FROZEN) {
	    // calculate Wr based on the ghost cell reconstruction
	    Wr = PiecewiseLinearSolutionAtLocation(i+1,j,
						   Grid.xfaceE(i,j));
	  } else if (Grid.BCtypeE[j] == BC_EXACT_SOLUTION) {
	    // calculate Wr using the exact solution at the interface and avoiding the problem to become ill-posed.
	    if (ExactSoln->IsExactSolutionSet()){
	      Wr = BC_Characteristic_Pressure(Wl,
					      ExactSoln->Solution(Grid.xfaceE(i,j).x,Grid.xfaceE(i,j).y),
					      Grid.nfaceE(i, j));
	    } else {
	      throw runtime_error("MHD2D_Quad_Block::dUdt_Residual_Evaluation() ERROR! There is no exact solution set for the Exact_Solution BC.");
	    }	    
	  } else if (Grid.BCtypeE[j] == BC_WALL_INVISCID) {
	    Wr = Reflect(Wl, Grid.nfaceE(i, j));
	    Wi = PiecewiseLinearSolutionAtLocation(i+1,j,
						   Grid.xfaceW(i+1, j));
	    Wl.d() = Wi.d();
	    Wl.p() = Wi.p();
	  } /* endif */
	} else { 
	  // Determine left state
	  Wl = PiecewiseLinearSolutionAtLocation(i,j,
						 Grid.xfaceE(i,j));
	  // Determine right state
	  Wr = PiecewiseLinearSolutionAtLocation(i+1,j,
						 Grid.xfaceW(i+1, j));	 
	} /* endif */

	// Determine the flux
	Flux = RiemannFlux_n(IP.i_Flux_Function,
			     Wl, Wr, Grid.nfaceE(i, j));

	/* Evaluate cell-averaged solution changes. */
	dUdt[i  ][j][0] -= Flux*Grid.lfaceE(i, j)/Grid.Cell[i][j].A;
	dUdt[i+1][j][0] += Flux*Grid.lfaceW(i+1, j)/Grid.Cell[i+1][j].A;

#if 0
	/* Include axisymmetric source terms as required. */
	if (Axisymmetric) {
	  dUdt[i][j][0] += S(W[i][j], Grid.Cell[i][j].Xc);
	} /* endif */
#endif

	/* Save west and east face boundary flux. */
	if (i == ICl-1) {
	  FluxW[j] = -Flux*Grid.lfaceW(i+1, j);
	} else if (i == ICu) {
	  FluxE[j] = Flux*Grid.lfaceE(i, j);
	} /* endif */ 

      } /* endif */
    } /* endfor */
      
    if ( j > JCl-1 && j < JCu+1 ) {
      dUdt[ICl-1][j][0].zero_all();
      dUdt[ICu+1][j][0].zero_all();
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
	   Grid.BCtypeS[i] == BC_FROZEN ||
	   Grid.BCtypeS[i] == BC_EXACT_SOLUTION ||
	   Grid.BCtypeS[i] == BC_WALL_INVISCID)) {
	// Determine right state
	Wr = PiecewiseLinearSolutionAtLocation(i, j+1,
					       Grid.xfaceS(i, j+1));	       
	// Determine left state
	if (Grid.BCtypeS[i] == BC_REFLECTION) {
	  Wl = Reflect(Wr, Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
	  Wl = BC_Characteristic_Pressure(Wr, 
					  WoS[i], 
					  Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_CHARACTERISTIC_VELOCITY) {
	  Wl = BC_Characteristic(Wr, 
				 WoS[i], 
				 Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_FROZEN) {
	  // calculate Wl based on the ghost cell reconstruction
	  Wl = PiecewiseLinearSolutionAtLocation(i,j,
						 Grid.xfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_EXACT_SOLUTION) {
	  // calculate Wl using the exact solution at the interface and avoiding the problem to become ill-posed.
	  if (ExactSoln->IsExactSolutionSet()){
	    Wl = BC_Characteristic_Pressure(Wr,
					    ExactSoln->Solution(Grid.xfaceS(i,j+1).x,Grid.xfaceS(i,j+1).y),
					    Grid.nfaceS(i, j+1));
	  } else {
	    throw runtime_error("MHD2D_Quad_Block::dUdt_Residual_Evaluation() ERROR! There is no exact solution set for the Exact_Solution BC.");
	  }
	} else if (Grid.BCtypeS[i] == BC_WALL_INVISCID) {
	  Wl = Reflect(Wr, Grid.nfaceS(i, j+1));
	  Wi = PiecewiseLinearSolutionAtLocation(i,j,
						 Grid.xfaceS(i, j+1));
	  Wl.d() = Wi.d();
	  Wl.p() = Wi.p();
	} /* endif */
      } else if (j == JCu && 
		 (Grid.BCtypeN[i] == BC_REFLECTION ||
		  Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
		  Grid.BCtypeN[i] == BC_CHARACTERISTIC_VELOCITY ||
		  Grid.BCtypeN[i] == BC_FROZEN ||
		  Grid.BCtypeN[i] == BC_EXACT_SOLUTION ||
		  Grid.BCtypeN[i] == BC_WALL_INVISCID )) {
	// Determine left state
	Wl = PiecewiseLinearSolutionAtLocation(i,j,
					       Grid.xfaceN(i, j));	       
	// Determine right state
	if (Grid.BCtypeN[i] == BC_REFLECTION) {
	  Wr = Reflect(Wl, Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_CHARACTERISTIC) {
	  Wr = BC_Characteristic_Pressure(Wl, 
					  WoN[i], 
					  Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_CHARACTERISTIC_VELOCITY) {
	  Wr = BC_Characteristic(Wl, 
				 WoN[i], 
				 Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_FROZEN) {
	  // calculate Wr based on the ghost cell reconstruction
	  Wr = PiecewiseLinearSolutionAtLocation(i,j+1,
						 Grid.xfaceN(i,j));
	} else if (Grid.BCtypeN[i] == BC_EXACT_SOLUTION) {
	  // calculate Wr using the exact solution at the interface and avoiding the problem to become ill-posed.
	  if (ExactSoln->IsExactSolutionSet()){
	    Wr = BC_Characteristic_Pressure(Wl,
					    ExactSoln->Solution(Grid.xfaceN(i,j).x,Grid.xfaceN(i,j).y),
					    Grid.nfaceN(i, j));
	  } else {
	    throw runtime_error("MHD2D_Quad_Block::dUdt_Residual_Evaluation() ERROR! There is no exact solution set for the Exact_Solution BC.");
	  }
	} else if (Grid.BCtypeN[i] == BC_WALL_INVISCID) {
	  Wr = Reflect(Wl, Grid.nfaceN(i, j));
	  Wi = PiecewiseLinearSolutionAtLocation(i,j+1,
						 Grid.xfaceN(i,j));
	  Wr.d() = Wi.d();
	  Wr.p() = Wi.p();
	} /* endif */
      } else {
	// Determine left state
	Wl = PiecewiseLinearSolutionAtLocation(i,j,
					       Grid.xfaceN(i, j));	       
	// Determine right state
	Wr = PiecewiseLinearSolutionAtLocation(i,j+1,
					       Grid.xfaceS(i, j+1));	       
      } /* endif */

      // Determine the flux
      Flux = RiemannFlux_n(IP.i_Flux_Function,
			   Wl, Wr, Grid.nfaceN(i, j));
	
      /* Evaluate cell-averaged solution changes. */
      dUdt[i][j  ][0] -= Flux*Grid.lfaceN(i, j)/Grid.Cell[i][j].A;
      dUdt[i][j+1][0] += Flux*Grid.lfaceS(i, j+1)/Grid.Cell[i][j+1].A;

      /* Save south and north face boundary flux. */
      if (j == JCl-1) {
	FluxS[i] = -Flux*Grid.lfaceS(i, j+1);
      } else if (j == JCu) {
	FluxN[i] = Flux*Grid.lfaceN(i, j);
      } /* endif */
	
    } /* endfor */
      
    dUdt[i][JCl-1][0].zero_all();
    dUdt[i][JCu+1][0].zero_all();
  } /* endfor */
    
    /* Residual successfully evaluated. */
  return 0;

}

/*!
 * This routine determines the solution residuals for a 
 * given stage of a variety of multi-stage explicit     
 * time integration schemes for a given solution block. 
 */
int MHD2D_Quad_Block::dUdt_Multistage_Explicit(const int &i_stage,
					       const MHD2D_Input_Parameters &IP){
  
  int i, j, k_residual;
  double omega;
  Vector2D dX;
  MHD3D_pState Wl, Wr, Wi;
  MHD3D_cState Flux;

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
  ComputeLinearSolutionReconstruction(IP.i_Reconstruction,
				      IP.i_Limiter);

  /* Evaluate the time rate of change of the solution
     (i.e., the solution residuals) using a second-order
     limited upwind scheme with a variety of flux functions. */
  
  // Add i-direction (zeta-direction) fluxes.
  for ( j  = JCl-1 ; j <= JCu+1 ; ++j ) {
    if ( i_stage == 1 ) {
      Uo[ICl-1][j] = U[ICl-1][j];
      dUdt[ICl-1][j][k_residual].zero_all();
    } else {
      dUdt[ICl-1][j][k_residual].zero_all();
    } /* endif */
    
    for ( i = ICl-1 ; i <= ICu ; ++i ) {
      if ( i_stage == 1 ) {
	Uo[i+1][j] = U[i+1][j];
	dUdt[i+1][j][k_residual].zero_all();
      } else if ( j > JCl-1 && j < JCu+1 ) {
	switch(IP.i_Time_Integration) {
	case TIME_STEPPING_EXPLICIT_PREDICTOR_CORRECTOR :
	  break;
	case TIME_STEPPING_EXPLICIT_RUNGE_KUTTA :
	  if (IP.N_Stage == 2) {
	  } else if (IP.N_Stage == 4 && i_stage == 4) {
	    dUdt[i+1][j][k_residual] = 
	      dUdt[i+1][j][0] + 
	      TWO*dUdt[i+1][j][1] +
	      TWO*dUdt[i+1][j][2];
	  } else {
	    dUdt[i+1][j][k_residual].zero_all();
	  } /* endif */
	  break;
	case TIME_STEPPING_MULTISTAGE_OPTIMAL_SMOOTHING :
	  dUdt[i+1][j][k_residual].zero_all();
	  break;
	default:
	  dUdt[i+1][j][k_residual].zero_all();
	  break;
	} /* endswitch */
      } /* endif */
    
      if ( j > JCl-1 && j < JCu+1 ) {
    
	/* Evaluate the cell interface i-direction fluxes. */
    
	if (i == ICl-1 && 
	    (Grid.BCtypeW[j] == BC_REFLECTION ||
	     Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
	     Grid.BCtypeW[j] == BC_CHARACTERISTIC_VELOCITY ||
	     Grid.BCtypeW[j] == BC_FROZEN ||
	     Grid.BCtypeW[j] == BC_EXACT_SOLUTION ||
	     Grid.BCtypeW[j] == BC_WALL_INVISCID)) {
	  // Determine right state
	  Wr = PiecewiseLinearSolutionAtLocation(i+1, j,
						 Grid.xfaceW(i+1, j));	 
	  // Determine left state
	  if (Grid.BCtypeW[j] == BC_REFLECTION) {
	    Wl = Reflect(Wr, Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_CHARACTERISTIC) {
	    Wl = BC_Characteristic_Pressure(Wr, 
					    WoW[j], 
					    Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_CHARACTERISTIC_VELOCITY) {
	    Wl = BC_Characteristic(Wr, 
				   WoW[j], 
				   Grid.nfaceW(i+1, j));
	  } else if (Grid.BCtypeW[j] == BC_FROZEN) {
	    // calculate Wl based on the ghost cell reconstruction
	    Wl = PiecewiseLinearSolutionAtLocation(i,j,
						   Grid.xfaceW(i+1,j));
	  } else if (Grid.BCtypeW[j] == BC_EXACT_SOLUTION) {
	    // calculate Wl using the exact solution at the interface and avoiding the problem to become ill-posed.
	    if (ExactSoln->IsExactSolutionSet()){
	      Wl = BC_Characteristic_Pressure(Wr,
					      ExactSoln->Solution(Grid.xfaceW(i+1,j).x,Grid.xfaceW(i+1,j).y),
					      Grid.nfaceW(i+1, j));
	    } else {
	      throw runtime_error("MHD2D_Quad_Block::dUdt_Residual_Evaluation() ERROR! There is no exact solution set for the Exact_Solution BC.");
	    }
	  } else if (Grid.BCtypeW[j] == BC_WALL_INVISCID) {
	    Wl = Reflect(Wr, Grid.nfaceW(i+1, j));
	    Wi = PiecewiseLinearSolutionAtLocation(i,j,
						  Grid.xfaceW(i+1,j));
	    Wl.d() = Wi.d();
	    Wl.p() = Wi.p();
	  } /* endif */
	} else if (i == ICu && 
		   (Grid.BCtypeE[j] == BC_REFLECTION ||
		    Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
		    Grid.BCtypeE[j] == BC_CHARACTERISTIC_VELOCITY ||
		    Grid.BCtypeE[j] == BC_FROZEN ||
		    Grid.BCtypeE[j] == BC_EXACT_SOLUTION ||
		    Grid.BCtypeE[j] == BC_WALL_INVISCID )) {
	  // Determine left state
	  Wl = PiecewiseLinearSolutionAtLocation(i,j,
						 Grid.xfaceE(i, j));
	  // Determine right state
	  if (Grid.BCtypeE[j] == BC_REFLECTION) {
	    Wr = Reflect(Wl, Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_CHARACTERISTIC) {
	    Wr = BC_Characteristic_Pressure(Wl, 
					    WoE[j], 
					    Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_CHARACTERISTIC_VELOCITY) {
	    Wr = BC_Characteristic(Wl, 
				   WoE[j], 
				   Grid.nfaceE(i, j));
	  } else if (Grid.BCtypeE[j] == BC_FROZEN) {
	    // calculate Wr based on the ghost cell reconstruction
	    Wr = PiecewiseLinearSolutionAtLocation(i+1,j,
						   Grid.xfaceE(i,j));
	  } else if (Grid.BCtypeE[j] == BC_EXACT_SOLUTION) {
	    // calculate Wr using the exact solution at the interface and avoiding the problem to become ill-posed.
	    if (ExactSoln->IsExactSolutionSet()){
	      Wr = BC_Characteristic_Pressure(Wl,
					      ExactSoln->Solution(Grid.xfaceE(i,j).x,Grid.xfaceE(i,j).y),
					      Grid.nfaceE(i, j));
	    } else {
	      throw runtime_error("MHD2D_Quad_Block::dUdt_Residual_Evaluation() ERROR! There is no exact solution set for the Exact_Solution BC.");
	    }	    
	  } else if (Grid.BCtypeE[j] == BC_WALL_INVISCID) {
	    Wr = Reflect(Wl, Grid.nfaceE(i, j));
	    Wi = PiecewiseLinearSolutionAtLocation(i+1,j,
						   Grid.xfaceW(i+1, j));
	    Wl.d() = Wi.d();
	    Wl.p() = Wi.p();
	  } /* endif */
	} else {            
	  // Determine left state
	  Wl = PiecewiseLinearSolutionAtLocation(i,j,
						 Grid.xfaceE(i, j));
	  // Determine right state
	  Wr = PiecewiseLinearSolutionAtLocation(i+1,j,
						 Grid.xfaceW(i+1, j));	 
	} /* endif */

	// Determine the flux 
	Flux = RiemannFlux_n(IP.i_Flux_Function,
			     Wl, Wr, Grid.nfaceE(i, j));
    
	/* Evaluate cell-averaged solution changes. */
	dUdt[i  ][j][k_residual] -= (IP.CFL_Number*dt[i  ][j])* Flux*Grid.lfaceE(i  , j)/ Grid.Cell[i  ][j].A;
	dUdt[i+1][j][k_residual] += (IP.CFL_Number*dt[i+1][j])* Flux*Grid.lfaceW(i+1, j)/ Grid.Cell[i+1][j].A;

#if 0
	/* Include axisymmetric source terms as required. */
	if (Axisymmetric) {
	  dUdt[i][j][k_residual] += (IP.CFL_Number*dt[i][j])*S(W[i][j], Grid.Cell[i][j].Xc);
	} /* endif */
#endif

	/* Save west and east face boundary flux. */
	if (i == ICl-1) {
	  FluxW[j] = -Flux*Grid.lfaceW(i+1, j);
	} else if (i == ICu) {
	  FluxE[j] = Flux*Grid.lfaceE(i, j);
	} /* endif */ 

      } /* endif */
    } /* endfor */
    
    if ( j > JCl-1 && j < JCu+1 ) {
      dUdt[ICl-1][j][k_residual].zero_all();
      dUdt[ICu+1][j][k_residual].zero_all();
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
	   Grid.BCtypeS[i] == BC_FROZEN ||
	   Grid.BCtypeS[i] == BC_EXACT_SOLUTION ||
	   Grid.BCtypeS[i] == BC_WALL_INVISCID )) {
	// Determine right state
	Wr = PiecewiseLinearSolutionAtLocation(i, j+1,
					       Grid.xfaceS(i, j+1));
	// Determine left state
	if (Grid.BCtypeS[i] == BC_REFLECTION) {
	  Wl = Reflect(Wr, Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
	  Wl = BC_Characteristic_Pressure(Wr, 
					  WoS[i], 
					  Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_CHARACTERISTIC_VELOCITY) {
	  Wl = BC_Characteristic(Wr, 
				 WoS[i], 
				 Grid.nfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_FROZEN) {
	  // calculate Wl based on the ghost cell reconstruction
	  Wl = PiecewiseLinearSolutionAtLocation(i,j,
						 Grid.xfaceS(i, j+1));
	} else if (Grid.BCtypeS[i] == BC_EXACT_SOLUTION) {
	  // calculate Wl using the exact solution at the interface and avoiding the problem to become ill-posed.
	  if (ExactSoln->IsExactSolutionSet()){
	    Wl = BC_Characteristic_Pressure(Wr,
					    ExactSoln->Solution(Grid.xfaceS(i,j+1).x,Grid.xfaceS(i,j+1).y),
					    Grid.nfaceS(i, j+1));
	  } else {
	    throw runtime_error("MHD2D_Quad_Block::dUdt_Residual_Evaluation() ERROR! There is no exact solution set for the Exact_Solution BC.");
	  }
	} else if (Grid.BCtypeS[i] == BC_WALL_INVISCID) {
	  Wl = Reflect(Wr, Grid.nfaceS(i, j+1));
	  Wi = PiecewiseLinearSolutionAtLocation(i,j,
						 Grid.xfaceS(i, j+1));
	  Wl.d() = Wi.d();
	  Wl.p() = Wi.p();
	} /* endif */
      } else if (j == JCu && 
		 (Grid.BCtypeN[i] == BC_REFLECTION ||
		  Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
		  Grid.BCtypeN[i] == BC_CHARACTERISTIC_VELOCITY ||
		  Grid.BCtypeN[i] == BC_FROZEN ||
		  Grid.BCtypeN[i] == BC_EXACT_SOLUTION ||
		  Grid.BCtypeN[i] == BC_WALL_INVISCID )) {
	// Determine left state
	Wl = PiecewiseLinearSolutionAtLocation(i,j,
					       Grid.xfaceN(i, j));	       
	// Determine right state
	if (Grid.BCtypeN[i] == BC_REFLECTION) {
	  Wr = Reflect(Wl, Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_CHARACTERISTIC) {
	  Wr = BC_Characteristic_Pressure(Wl, 
					  WoN[i], 
					  Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_CHARACTERISTIC_VELOCITY) {
	  Wr = BC_Characteristic(Wl, 
				 WoN[i], 
				 Grid.nfaceN(i, j));
	} else if (Grid.BCtypeN[i] == BC_FROZEN) {
	  // calculate Wr based on the ghost cell reconstruction
	  Wr = PiecewiseLinearSolutionAtLocation(i,j+1,
						 Grid.xfaceN(i,j));
	} else if (Grid.BCtypeN[i] == BC_EXACT_SOLUTION) {
	  // calculate Wr using the exact solution at the interface and avoiding the problem to become ill-posed.
	  if (ExactSoln->IsExactSolutionSet()){
	    Wr = BC_Characteristic_Pressure(Wl,
					    ExactSoln->Solution(Grid.xfaceN(i,j).x,Grid.xfaceN(i,j).y),
					    Grid.nfaceN(i, j));
	  } else {
	    throw runtime_error("MHD2D_Quad_Block::dUdt_Residual_Evaluation() ERROR! There is no exact solution set for the Exact_Solution BC.");
	  }
	} else if (Grid.BCtypeN[i] == BC_WALL_INVISCID) {
	  Wr = Reflect(Wl, Grid.nfaceN(i, j));
	  Wi = PiecewiseLinearSolutionAtLocation(i,j+1,
						 Grid.xfaceN(i,j));
	  Wr.d() = Wi.d();
	  Wr.p() = Wi.p();
	} /* endif */
      } else {
	// Determine left state
	Wl = PiecewiseLinearSolutionAtLocation(i,j,
					       Grid.xfaceN(i, j));	       
	// Determine right state
	Wr = PiecewiseLinearSolutionAtLocation(i,j+1,
					       Grid.xfaceS(i, j+1));	       
      } /* endif */

      // Determine the flux
      Flux = RiemannFlux_n(IP.i_Flux_Function,
			   Wl, Wr, Grid.nfaceN(i, j));
    
      /* Evaluate cell-averaged solution changes. */
      dUdt[i][j  ][k_residual] -= (IP.CFL_Number*dt[i][j  ])*Flux*Grid.lfaceN(i, j  )/Grid.Cell[i][j  ].A;
      dUdt[i][j+1][k_residual] += (IP.CFL_Number*dt[i][j+1])*Flux*Grid.lfaceS(i, j+1)/Grid.Cell[i][j+1].A;

      /* Save south and north face boundary flux. */
      if (j == JCl-1) {
	FluxS[i] = -Flux*Grid.lfaceS(i, j+1);
      } else if (j == JCu) {
	FluxN[i] = Flux*Grid.lfaceN(i, j);
      } /* endif */

    } /* endfor */
    
    dUdt[i][JCl-1][k_residual].zero_all();
    dUdt[i][JCu+1][k_residual].zero_all();
  } /* endfor */
    
  /* Residual for the stage successfully calculated. */

  return (0);

}
