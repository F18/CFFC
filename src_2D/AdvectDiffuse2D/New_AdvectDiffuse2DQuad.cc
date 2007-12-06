/*! \file AdvectDiffuse2DState.cc
  @brief Subroutines for 2D Advection Diffusion Equation Quadrilateral Mesh Solution Classes. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "New_AdvectDiffuse2DQuad.h"        // AdvectDiffuse2D_Quad_Block class

/*********************************************************************
 * Static Variable Initialization for AdvectDiffuse2D_Quad_Block_New *
 *********************************************************************/

// Initialize residual_variable
int AdvectDiffuse2D_Quad_Block_New::residual_variable = 1;
// Initialize ExactSoln
AdvectDiffuse2D_ExactSolutions *AdvectDiffuse2D_Quad_Block_New::ExactSoln = NULL;
// Initialize Flow_Type
int AdvectDiffuse2D_Quad_Block_New::Flow_Type = FLOWTYPE_INVISCID;
// Initialize Axisymmetric
int AdvectDiffuse2D_Quad_Block_New::Axisymmetric = OFF;
// Initialize Number_of_Residual_Norms
int AdvectDiffuse2D_Quad_Block_New::Number_of_Residual_Norms = 1;


/*******************************************************************************
 * AdvectDiffuse2D_Quad_Block -- Single Block Member Functions.                *
 ******************************************************************************/
/**********************
 * Default constructor.
 **********************/
AdvectDiffuse2D_Quad_Block_New::AdvectDiffuse2D_Quad_Block_New(void): AssessAccuracy(this) {
  Freeze_Limiter = OFF;
  // Grid size and variables:
  NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
  // Solution variables:
  U = NULL; dt = NULL; dUdt = NULL; 
  dUdx = NULL; dUdy = NULL; phi = NULL; Uo = NULL;
  FluxN = NULL; FluxS = NULL; FluxE = NULL; FluxW = NULL;
  UoN = NULL; UoS = NULL; UoE = NULL; UoW = NULL;

  // Set the pointers in the solution state to the velocity, diffusion and source term fields.
  // These pointes are static variables in the state class.
  (*U)->Set_Pointers_To_Fields();

  // Get access to the AdvectDiffuse2D_ExactSolutions object
  ExactSoln = &AdvectDiffuse2D_ExactSolutions::getInstance();
}

/******************************************
 * Private copy constructor. (shallow copy)
 *****************************************/
AdvectDiffuse2D_Quad_Block_New::AdvectDiffuse2D_Quad_Block_New(const AdvectDiffuse2D_Quad_Block_New &Soln):
  AssessAccuracy(this)
{
  NCi = Soln.NCi; ICl = Soln.ICl; ICu = Soln.ICu; 
  NCj = Soln.NCj; JCl = Soln.JCl; JCu = Soln.JCu; Nghost = Soln.Nghost;
  Grid = Soln.Grid; U = Soln.U; dt = Soln.dt; dUdt = Soln.dUdt; 
  dUdx = Soln.dUdx; dUdy = Soln.dUdy; phi = Soln.phi;
  Uo = Soln.Uo;
  FluxN = Soln.FluxN; FluxS = Soln.FluxS; FluxE = Soln.FluxE; FluxW = Soln.FluxW;
  UoN = Soln.UoN; UoS = Soln.UoS; UoE = Soln.UoE; UoW = Soln.UoW;
  Freeze_Limiter = Soln.Freeze_Limiter;
}

/**********************
 * Allocate memory.            
 **********************/
void AdvectDiffuse2D_Quad_Block_New::allocate(const int Ni, const int Nj, const int Ng) {
  int i, j, k; assert(Ni > 1 && Nj > 1 && Ng > 1 && Ng > 1); Grid.allocate(Ni, Nj, Ng);
  NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1;
  NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1; Nghost = Ng;
  U = new AdvectDiffuse2D_State_New*[NCi]; dt = new double*[NCi]; dUdt = new AdvectDiffuse2D_State_New**[NCi]; 
  dUdx = new AdvectDiffuse2D_State_New*[NCi]; dUdy = new AdvectDiffuse2D_State_New*[NCi]; 
  phi = new AdvectDiffuse2D_State_New*[NCi]; Uo = new AdvectDiffuse2D_State_New*[NCi];
  for ( i = 0; i <= NCi-1 ; ++i ) {
    U[i] = new AdvectDiffuse2D_State_New[NCj]; 
    dt[i] = new double[NCj]; dUdt[i] = new AdvectDiffuse2D_State_New*[NCj];
    for ( j = 0; j <= NCj-1 ; ++j ) 
      { dUdt[i][j] = new AdvectDiffuse2D_State_New[NUMBER_OF_RESIDUAL_VECTORS_ADVECTDIFFUSE2D]; }
    dUdx[i] = new AdvectDiffuse2D_State_New[NCj]; dUdy[i] = new AdvectDiffuse2D_State_New[NCj]; 
    phi[i] = new AdvectDiffuse2D_State_New[NCj];
    Uo[i] = new AdvectDiffuse2D_State_New[NCj];
  } /* endfor */
  FluxN = new AdvectDiffuse2D_State_New[NCi]; FluxS = new AdvectDiffuse2D_State_New[NCi];
  FluxE = new AdvectDiffuse2D_State_New[NCj]; FluxW = new AdvectDiffuse2D_State_New[NCj];
  UoN = new AdvectDiffuse2D_State_New[NCi]; UoS = new AdvectDiffuse2D_State_New[NCi];
  UoE = new AdvectDiffuse2D_State_New[NCj]; UoW = new AdvectDiffuse2D_State_New[NCj];
  // Set the solution residuals, gradients, limiters, and other values to zero.
  for (j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
    for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
      for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_ADVECTDIFFUSE2D-1 ; ++k ) {
	dUdt[i][j][k].Vacuum();
      }
      dUdx[i][j].Vacuum(); dUdy[i][j].Vacuum();
      phi[i][j].Vacuum();
      Uo[i][j].Vacuum();
      dt[i][j] = ZERO;
    } /* endfor */
  } /* endfor */
}

/***********************
 * Deallocate memory.   
 ***********************/
void AdvectDiffuse2D_Quad_Block_New::deallocate(void) {
  int i, j; Grid.deallocate();
  for ( i = 0; i <= NCi-1 ; ++i ) {
    delete []U[i]; U[i] = NULL;
    delete []dt[i]; dt[i] = NULL; 
    for ( j = 0; j <= NCj-1 ; ++j ) { delete []dUdt[i][j]; dUdt[i][j] = NULL; }
    delete []dUdt[i]; dUdt[i] = NULL;
    delete []dUdx[i]; dUdx[i] = NULL; delete []dUdy[i]; dUdy[i] = NULL;
    delete []phi[i]; phi[i] = NULL; delete []Uo[i]; Uo[i] = NULL;
  } /* endfor */
  delete []U; U = NULL; delete []dt; dt = NULL; delete []dUdt; dUdt = NULL;
  delete []dUdx; dUdx = NULL; delete []dUdy; dUdy = NULL; 
  delete []phi; phi = NULL; delete []Uo; Uo = NULL;
  delete []FluxN; FluxN = NULL; delete []FluxS; FluxS = NULL;
  delete []FluxE; FluxE = NULL; delete []FluxW; FluxW = NULL;
  delete []UoN; UoN = NULL; delete []UoS; UoS = NULL;
  delete []UoE; UoE = NULL; delete []UoW; UoW = NULL;
  NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
}

/***********************
 * Node solution state. 
 ***********************/
AdvectDiffuse2D_State_New AdvectDiffuse2D_Quad_Block_New::Un(const int ii, const int jj) {
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
    eta1, zeta1, eta2, zeta2, eta, zeta;
  double As, Bs, Cs, Ds;
  Vector2D Av, Bv, Cv, Dv;

  x=Grid.Node[ii][jj].X.x; y=Grid.Node[ii][jj].X.y;
  ax=Grid.Cell[ii-1][jj-1].Xc.x;
  bx=Grid.Cell[ii-1][jj].Xc.x-Grid.Cell[ii-1][jj-1].Xc.x; 
  cx=Grid.Cell[ii][jj-1].Xc.x-Grid.Cell[ii-1][jj-1].Xc.x; 
  dx=Grid.Cell[ii][jj].Xc.x+Grid.Cell[ii-1][jj-1].Xc.x-
    Grid.Cell[ii-1][jj].Xc.x-Grid.Cell[ii][jj-1].Xc.x;
  ay=Grid.Cell[ii-1][jj-1].Xc.y; 
  by=Grid.Cell[ii-1][jj].Xc.y-Grid.Cell[ii-1][jj-1].Xc.y; 
  cy=Grid.Cell[ii][jj-1].Xc.y-Grid.Cell[ii-1][jj-1].Xc.y; 
  dy=Grid.Cell[ii][jj].Xc.y+Grid.Cell[ii-1][jj-1].Xc.y-
    Grid.Cell[ii-1][jj].Xc.y-Grid.Cell[ii][jj-1].Xc.y;
  aa=bx*dy-dx*by; bb=dy*(ax-x)+bx*cy-cx*by+dx*(y-ay); cc=cy*(ax-x)+cx*(y-ay);
  if (fabs(aa) < TOLER*TOLER) {
    if (fabs(bb) >= TOLER*TOLER) { zeta1=-cc/bb; }
    else { zeta1 = -cc/sgn(bb)*(TOLER*TOLER); } 
    if (fabs(cy+dy*zeta1) >= TOLER*TOLER) { eta1=(y-ay-by*zeta1)/(cy+dy*zeta1); } 
    else { eta1 = HALF; } zeta2=zeta1; eta2=eta1;
  } else {
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) { zeta1=HALF*(-bb+sqrt(bb*bb-FOUR*aa*cc))/aa; }
    else { zeta1 = -HALF*bb/aa; } 
    if (fabs(cy+dy*zeta1) < TOLER*TOLER) { eta1=-ONE; } 
    else { eta1=(y-ay-by*zeta1)/(cy+dy*zeta1); } 
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) { zeta2=HALF*(-bb-sqrt(bb*bb-FOUR*aa*cc))/aa; }
    else { zeta2 = -HALF*bb/aa; }
    if (fabs(cy+dy*zeta2) < TOLER*TOLER) { eta2=-ONE; } 
    else { eta2=(y-ay-by*zeta2)/(cy+dy*zeta2); }
  } /* end if */
  if (zeta1 > -TOLER && zeta1 < ONE + TOLER &&
      eta1  > -TOLER && eta1  < ONE + TOLER) {
    zeta=zeta1; eta=eta1;
  } else if (zeta2 > -TOLER && zeta2 < ONE + TOLER &&
	     eta2  > -TOLER && eta2  < ONE + TOLER) {
    zeta=zeta2; eta=eta2;
  } else {
    zeta=HALF; eta=HALF;
  } /* endif */
  As=U[ii-1][jj-1].u; Bs=U[ii-1][jj].u-U[ii-1][jj-1].u; Cs=U[ii][jj-1].u-U[ii-1][jj-1].u;
  Ds=U[ii][jj].u+U[ii-1][jj-1].u-U[ii-1][jj].u-U[ii][jj-1].u;

  return (AdvectDiffuse2D_State_New(As+Bs*zeta+Cs*eta+Ds*zeta*eta));
}

/*********************************************//**
 * Compute the average source term for cell (ii,jj).
 * 
 * \todo Add integration with high-order reconstruction!
 *************************************************/
AdvectDiffuse2D_State_New AdvectDiffuse2D_Quad_Block_New::SourceTerm(const int & ii, const int & jj) const{
  
  double SourceTermIntegral(0.0), _dummy_(0.0);

  if (U[ii][jj].SourceTerm->FieldRequireIntegration()){
    // this source term field requires numerical integration

    // Integrate the non-linear source term field with the high-order reconstruction
    


    // Integrate the non-linear source term field with the low-order reconstruction

    // Define a source term functional that uses the piecewise linear reconstruction 
    // as a closure for the non-linear source variation.
    SourceTermFunctionalWithPiecewiseLinear NonLinearSourceVariation(ii,jj,this);

    // Integrate the source term functional over the cell (ii,jj) domain
    SourceTermIntegral = Grid.Integration.IntegrateFunctionOverCell(ii,jj,NonLinearSourceVariation,10,_dummy_);

    // compute average value and cast the result to AdvectDiffuse2D_State_New
    return AdvectDiffuse2D_State_New(SourceTermIntegral/Grid.Cell[ii][jj].A);

  } else {
    // this source term field is already integrated numerically and expressed as a function of cell average solution
    return AdvectDiffuse2D_State_New(U[ii][jj].s());
  }
}

/*********************************************//**
 * Compute the inviscid-flux ghost-cell state at a 
 * boundary interface. 
 * This state depends on the boundary conditions that 
 * need to be enforced.
 * 
 * \param BOUNDARY boundary position specifier (i.e. WEST, EAST, SOUTH or NORTH)
 * \param ii i-index of the cell in which the calculation is done
 * \param ii j-index of the cell in which the calculation is done
 * \param Ul the left interface state
 * \param Ur the right interface state
 *************************************************/
void AdvectDiffuse2D_Quad_Block_New::InviscidFluxStateAtBoundaryInterface(const int &BOUNDARY,
									  const int &ii, const int &jj,
									  AdvectDiffuse2D_State_New &Ul,
									  AdvectDiffuse2D_State_New &Ur){
  double Vn;

  switch(BOUNDARY){

    // *******************************
    // === WEST boundary interface ===
    // *******************************
  case WEST :			
    // Compute left interface state based on reconstruction or the particular boundary condition
    switch (Grid.BCtypeW[jj]){
      // Category I
    case BC_NONE :
    case BC_PERIODIC :
      // Calculate Ul based on the reconstruction in the ghost cell
      Ul = PiecewiseLinearSolutionAtLocation(ii  , jj,Grid.xfaceE(ii  , jj));
      break;
      
      // Category II
    case BC_INFLOW :
    case BC_DIRICHLET :
      Ul = UoW[jj];
      break;
      
      // Category III
    case BC_NEUMANN : 
    case BC_SYMMETRY_PLANE :
    case BC_EXTRAPOLATE :
    case BC_LINEAR_EXTRAPOLATION :
    case BC_OUTFLOW :
    case BC_CONSTANT_EXTRAPOLATION :
      Ul = Ur;
      break;
      
      // Category IV
    case BC_FARFIELD :
      /* Farfield BC is implemented differently for flows that 
	 enter the domain than for flows that leave the domain.
	 Whether the flow enters or leaves the domain is decided based on
	 the normal component of the velocity at the face midpoint.
	 --> If the flow enters the domain then the reference data is used.
	 --> If the flow leaves the domain then the interior value is used. */

      // Compute the normal velocity
      Vn = dot(VelocityAtLocation(Grid.xfaceW(ICl,jj)),
	       Grid.nfaceW(ICl,jj));
	  
      if (Vn <= ZERO){
	// The flow enters the domain
	// Use the boundary condition value
	Ul = UoW[jj];
      } else {
	// The flow leaves the domain
	// Use the interior domain value
	Ul = Ur;
      }
      break;

    default:
      Print_(Grid.BCtypeW[jj])
      throw runtime_error("AdvectDiffuse2D_Quad_Block_New::InviscidFluxStateAtBoundaryInterface() ERROR! No such West BCtype!");
    }// endswitch (Grid.BCtypeW[jj])
    break;


    // *******************************
    // === EAST boundary interface ===
    // *******************************
  case EAST :
    // Compute right interface state based on reconstruction or the particular boundary condition
    switch (Grid.BCtypeE[jj]){
      // Category I
    case BC_NONE :
    case BC_PERIODIC :
      // Calculate Ur based on the reconstruction in the ghost cell
      Ur = PiecewiseLinearSolutionAtLocation(ii+1, jj,Grid.xfaceW(ii+1, jj));
      break;
      
      // Category II
    case BC_INFLOW :
    case BC_DIRICHLET :
      Ur = UoE[jj];
      break;
      
      // Category III
    case BC_NEUMANN : 
    case BC_SYMMETRY_PLANE :
    case BC_EXTRAPOLATE :
    case BC_LINEAR_EXTRAPOLATION :
    case BC_OUTFLOW :
    case BC_CONSTANT_EXTRAPOLATION :
      Ur = Ul;
      break;
      
      // Category IV
    case BC_FARFIELD :
      /* Farfield BC is implemented differently for flows that 
	 enter the domain than for flows that leave the domain.
	 Whether the flow enters or leaves the domain is decided based on
	 the normal component of the velocity at the face midpoint.
	 --> If the flow enters the domain then the reference data is used.
	 --> If the flow leaves the domain then the interior value is used. */

      // Compute the normal velocity
      Vn = dot(VelocityAtLocation(Grid.xfaceE(ICu,jj)),
	       Grid.nfaceE(ICu,jj));
	  
      if (Vn <= ZERO){
	// The flow enters the domain
	// Use the boundary condition value
	Ur = UoE[jj];
      } else {
	// The flow leaves the domain
	// Use the interior domain value
	Ur = Ul;
      }
      break;

    default:
      throw runtime_error("AdvectDiffuse2D_Quad_Block_New::InviscidFluxStateAtBoundaryInterface() ERROR! No such East BCtype!");
    }// endswitch (Grid.BCtypeE[jj])
    break;

    // ********************************
    // === SOUTH boundary interface ===
    // ********************************
  case SOUTH :
    // Compute left interface state based on reconstruction or the particular boundary condition
    switch (Grid.BCtypeS[ii]){
      // Category I
    case BC_NONE :
    case BC_PERIODIC :
      // Calculate Ul based on the reconstruction in the ghost cell
      Ul = PiecewiseLinearSolutionAtLocation(ii ,jj  ,Grid.xfaceN(ii ,jj  ));
      break;
      
      // Category II
    case BC_INFLOW :
    case BC_DIRICHLET :
      Ul = UoS[ii];
      break;
      
      // Category III
    case BC_NEUMANN : 
    case BC_SYMMETRY_PLANE :
    case BC_EXTRAPOLATE :
    case BC_LINEAR_EXTRAPOLATION :
    case BC_OUTFLOW :
    case BC_CONSTANT_EXTRAPOLATION :
      Ul = Ur;
      break;
      
      // Category IV
    case BC_FARFIELD :
      /* Farfield BC is implemented differently for flows that 
	 enter the domain than for flows that leave the domain.
	 Whether the flow enters or leaves the domain is decided based on
	 the normal component of the velocity at the face midpoint.
	 --> If the flow enters the domain then the reference data is used.
	 --> If the flow leaves the domain then the interior value is used. */

      // Compute the normal velocity
      Vn = dot(VelocityAtLocation(Grid.xfaceS(ii,JCl)),
	       Grid.nfaceS(ii,JCl));
	  
      if (Vn <= ZERO){
	// The flow enters the domain
	// Use the boundary condition value
	Ul = UoS[ii];
      } else {
	// The flow leaves the domain
	// Use the interior domain value
	Ul = Ur;
      }
      break;

    default:
      throw runtime_error("AdvectDiffuse2D_Quad_Block_New::InviscidFluxStateAtBoundaryInterface() ERROR! No such South BCtype!");
    }// endswitch (Grid.BCtypeS[ii])
    break;

    // ********************************
    // === NORTH boundary interface ===
    // ********************************
  case NORTH :
    // Compute right interface state based on reconstruction or the particular boundary condition
    switch (Grid.BCtypeN[ii]){
      // Category I
    case BC_NONE :
    case BC_PERIODIC :
      // Calculate Ur based on the reconstruction in the ghost cell
      Ur = PiecewiseLinearSolutionAtLocation(ii ,jj+1,Grid.xfaceS(ii ,jj+1));
      break;
      
      // Category II
    case BC_INFLOW :
    case BC_DIRICHLET :
      Ur = UoN[ii];
      break;
      
      // Category III
    case BC_NEUMANN : 
    case BC_SYMMETRY_PLANE :
    case BC_EXTRAPOLATE :
    case BC_LINEAR_EXTRAPOLATION :
    case BC_OUTFLOW :
    case BC_CONSTANT_EXTRAPOLATION :
      Ur = Ul;
      break;
      
      // Category IV
    case BC_FARFIELD :
      /* Farfield BC is implemented differently for flows that 
	 enter the domain than for flows that leave the domain.
	 Whether the flow enters or leaves the domain is decided based on
	 the normal component of the velocity at the face midpoint.
	 --> If the flow enters the domain then the reference data is used.
	 --> If the flow leaves the domain then the interior value is used. */
      
      // Compute the normal velocity
      Vn = dot(VelocityAtLocation(Grid.xfaceN(ii,JCu)),
	       Grid.nfaceN(ii,JCu));
	  
      if (Vn <= ZERO){
	// The flow enters the domain
	// Use the boundary condition value
	Ur = UoN[ii];
      } else {
	// The flow leaves the domain
	// Use the interior domain value
	Ur = Ul;
      }
      break;

    default:
      throw runtime_error("AdvectDiffuse2D_Quad_Block_New::InviscidFluxStateAtBoundaryInterface() ERROR! No such North BCtype!");
    }// endswitch (Grid.BCtypeN[ii])
    break;

  default:
    throw runtime_error("AdvectDiffuse2D_Quad_Block_New::InviscidFluxStateAtBoundaryInterface() ERROR! No such boundary!");
  }

}

/********************
 * Output operator.
 ********************/
ostream &operator << (ostream &out_file,
		      const AdvectDiffuse2D_Quad_Block_New &SolnBlk) {
  int i, j; 
  out_file << SolnBlk.Grid;
  out_file << SolnBlk.NCi << " " << SolnBlk.ICl << " " << SolnBlk.ICu << "\n";
  out_file << SolnBlk.NCj << " " << SolnBlk.JCl << " " << SolnBlk.JCu << "\n";
  out_file << SolnBlk.Nghost << "\n";
  out_file << SolnBlk.Axisymmetric << "\n";
  if (SolnBlk.NCi == 0 || SolnBlk.NCj == 0) return(out_file);
  for ( j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
    for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
      out_file << SolnBlk.U[i][j] << "\n";
    } /* endfor */
  } /* endfor */
  for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
    out_file << SolnBlk.UoW[j] << "\n";
    out_file << SolnBlk.UoE[j] << "\n";
  } /* endfor */
  for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
    out_file << SolnBlk.UoS[i] << "\n";
    out_file << SolnBlk.UoN[i] << "\n";
  } /* endfor */
  return (out_file);
}

/********************
 * Input operator.
 ********************/
istream &operator >> (istream &in_file,
		      AdvectDiffuse2D_Quad_Block_New &SolnBlk) {
  int i, j, k, ni, il, iu, nj, jl, ju, ng;
  Grid2D_Quad_Block New_Grid;
  in_file >> New_Grid;
  in_file.setf(ios::skipws);
  in_file >> ni >> il >> iu;
  in_file >> nj >> jl >> ju;
  in_file >> ng;
  in_file >> SolnBlk.Axisymmetric;
  in_file.unsetf(ios::skipws);

  if (ni == 0 || nj == 0) {
    SolnBlk.deallocate(); return(in_file);
  } /* endif */
  if (SolnBlk.U == NULL || SolnBlk.NCi != ni || SolnBlk.NCj != nj || SolnBlk.Nghost != ng) {
    if (SolnBlk.U != NULL) SolnBlk.deallocate(); 
    SolnBlk.allocate(ni - 2*ng, nj - 2*ng, ng);
  } /* endif */

  // Copy the temporary mesh into the grid of the current solution block
  Copy_Quad_Block(SolnBlk.Grid,New_Grid); New_Grid.deallocate();

  // Read the solution & Initialize some data structures
  for ( j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
    for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
      in_file >> SolnBlk.U[i][j];
      for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_ADVECTDIFFUSE2D-1 ; ++k ) {
	SolnBlk.dUdt[i][j][k] = ZERO;
      } /* endfor */
      SolnBlk.dUdx[i][j].Vacuum();
      SolnBlk.dUdy[i][j].Vacuum();
      SolnBlk.phi[i][j].Vacuum();
      SolnBlk.Uo[i][j].Vacuum();
      SolnBlk.dt[i][j] = ZERO;
    } /* endfor */
  } /* endfor */

  for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
    in_file >> SolnBlk.UoW[j];
    in_file >> SolnBlk.UoE[j];
  } /* endfor */
  for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
    in_file >> SolnBlk.UoS[i];
    in_file >> SolnBlk.UoN[i];
  } /* endfor */

  in_file.setf(ios::skipws);
  return (in_file);
}

/*******************************************************************************
 *                                                                             *
 * MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING.                              *
 *                                                                             *
 *******************************************************************************/

/********************************
 * Loads send message buffer.
 ********************************/
int AdvectDiffuse2D_Quad_Block_New::LoadSendBuffer(double *buffer,
						   int &buffer_count,
						   const int buffer_size,
						   const int i_min, 
						   const int i_max,
						   const int i_inc,
						   const int j_min, 
						   const int j_max,
						   const int j_inc) {
  int i, j, k;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
    for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      for ( k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	buffer_count++;
	if (buffer_count >= buffer_size) return(1);
	buffer[buffer_count] = U[i][j][k];
      } /* endfor */
    } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * Loads send message buffer for fine to coarse block message passing.          
 *******************************************************************************/
int AdvectDiffuse2D_Quad_Block_New::LoadSendBuffer_F2C(double *buffer,
						       int &buffer_count,
						       const int buffer_size,
						       const int i_min, 
						       const int i_max,
						       const int i_inc,
						       const int j_min, 
						       const int j_max,
						       const int j_inc) {
  int i, j, k;
  for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
    for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
      for ( k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	buffer_count++;
	if (buffer_count >= buffer_size) return(1);
	buffer[buffer_count] = ( (Grid.Cell[i  ][j  ].A*U[i  ][j  ][k]+
				  Grid.Cell[i+1][j  ].A*U[i+1][j  ][k]+
				  Grid.Cell[i  ][j+1].A*U[i  ][j+1][k]+
				  Grid.Cell[i+1][j+1].A*U[i+1][j+1][k])/
				 (Grid.Cell[i  ][j  ].A+
				  Grid.Cell[i+1][j  ].A+
				  Grid.Cell[i  ][j+1].A+
				  Grid.Cell[i+1][j+1].A) );
      } /* endfor */
    } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * Loads send message buffer for coarse to fine block message passing.         
 *******************************************************************************/
int AdvectDiffuse2D_Quad_Block_New::LoadSendBuffer_C2F(double *buffer,
						       int &buffer_count,
						       const int buffer_size,
						       const int i_min, 
						       const int i_max,
						       const int i_inc,
						       const int j_min, 
						       const int j_max,
						       const int j_inc,
						       const int face,
						       const int sector) {
  int i, j, k;
  Vector2D dX;
  double ufine;
  AdvectDiffuse2D_State_New Ufine;
  int Limiter = LIMITER_VENKATAKRISHNAN;

  if (j_inc > 0) {
    if (i_inc > 0) {
      for ( j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i, j_min).
	  SubcellReconstruction(i, j, Limiter);
	  // Evaluate SW sub (fine) cell values if required.
	  if (!(face == NORTH && sector == WEST && Nghost%2 && j == j_min) &&
	      !(face == NORTH && sector == EAST && Nghost%2 && (i == i_min || j == j_min)) &&
	      !(face == SOUTH && sector == EAST && Nghost%2 && i == i_min) &&
	      !(face == EAST && sector == NORTH && Nghost%2 && (i == i_min || j == j_min)) &&
	      !(face == EAST && sector == SOUTH && Nghost%2 && i == i_min) &&
	      !(face == WEST && sector == NORTH && Nghost%2 && j == j_min) &&
	      !(face == NORTH_EAST && Nghost%2 && (i == i_min || j == j_min)) &&
	      !(face == NORTH_WEST && Nghost%2 && j == j_min) &&
	      !(face == SOUTH_EAST && Nghost%2 && i == i_min)) {
	    dX = Grid.centroidSW(i,j) - Grid.Cell[i][j].Xc;
	    Ufine = U[i][j] + ( (phi[i][j]^dUdx[i][j])*dX.x +
				(phi[i][j]^dUdy[i][j])*dX.y );
	    for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	      buffer_count++;
	      if (buffer_count >= buffer_size) return(1);
	      buffer[buffer_count] = Ufine[k];
	    } /* endfor (k) */
	  } /* endif */
	  // Evaluate SE sub (fine) cell values if required.
	  if (!(face == NORTH && sector == WEST && Nghost%2 && (i == i_max || j == j_min)) &&
	      !(face == NORTH && sector == EAST && Nghost%2 && j == j_min) &&
	      !(face == SOUTH && sector == WEST && Nghost%2 && i == i_max) &&
	      !(face == EAST && sector == NORTH && Nghost%2 && j == j_min) &&
	      !(face == WEST && sector == NORTH && Nghost%2 && (i == i_max || j == j_min)) &&
	      !(face == WEST && sector == SOUTH && Nghost%2 && i == i_max) &&
	      !(face == NORTH_EAST && Nghost%2 && j == j_min) &&
	      !(face == NORTH_WEST && Nghost%2 && (i == i_max || j == j_min)) &&
	      !(face == SOUTH_WEST && Nghost%2 && i == i_max)) {
	    dX = Grid.centroidSE(i,j) - Grid.Cell[i][j].Xc;
	    Ufine = U[i][j] + ( (phi[i][j]^dUdx[i][j])*dX.x +
				(phi[i][j]^dUdy[i][j])*dX.y );
	    for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	      buffer_count++;
	      if (buffer_count >= buffer_size) return(1);
	      buffer[buffer_count] = Ufine[k];
	    } /* endfor (k) */
	  } /* endif */
	} /* endfor */
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Evaluate NW sub (fine) cell values if required.
	  if (!(face == NORTH && sector == EAST && Nghost%2 && i == i_min) &&
	      !(face == SOUTH && sector == EAST && Nghost%2 && (i == i_min || j == j_max)) &&
	      !(face == SOUTH && sector == WEST && Nghost%2 && j == j_max) &&
	      !(face == EAST && sector == NORTH && Nghost%2 && i == i_min) &&
	      !(face == EAST && sector == SOUTH && Nghost%2 && (i == i_min || j == j_max)) &&
	      !(face == WEST && sector == SOUTH && Nghost%2 && j == j_max) &&
	      !(face == NORTH_EAST && Nghost%2 && i == i_min) &&
	      !(face == SOUTH_EAST && Nghost%2 && (i == i_min || j == j_max)) &&
	      !(face == SOUTH_WEST && Nghost%2 && j == j_max)) {
	    dX = Grid.centroidNW(i,j) - Grid.Cell[i][j].Xc;
	    Ufine = U[i][j] + ( (phi[i][j]^dUdx[i][j])*dX.x +
				(phi[i][j]^dUdy[i][j])*dX.y );
	    for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	      buffer_count++;
	      if (buffer_count >= buffer_size) return(1);
	      buffer[buffer_count] = Ufine[k];
	    } /* endfor (k) */
	  } /* endif */
	  // Evaluate NE sub (fine) cell values if required.
	  if (!(face == NORTH && sector == WEST && Nghost%2 && i == i_max) &&
	      !(face == SOUTH && sector == EAST && Nghost%2 && j == j_max) &&
	      !(face == SOUTH && sector == WEST && Nghost%2 && (i == i_max || j == j_max)) &&
	      !(face == EAST && sector == SOUTH && Nghost%2 && j == j_max) &&
	      !(face == WEST && sector == NORTH && Nghost%2 && i == i_max) &&
	      !(face == WEST && sector == SOUTH && Nghost%2 && (i == i_max || j == j_max)) &&
	      !(face == NORTH_WEST && Nghost%2 && i == i_max) &&
	      !(face == SOUTH_EAST && Nghost%2 && j == j_max) &&
	      !(face == SOUTH_WEST && Nghost%2 && (i == i_max || j == j_max))) {
	    dX = Grid.centroidNE(i,j) - Grid.Cell[i][j].Xc;
	    Ufine = U[i][j] + ( (phi[i][j]^dUdx[i][j])*dX.x +
				(phi[i][j]^dUdy[i][j])*dX.y );
	    for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	      buffer_count++;
	      if (buffer_count >= buffer_size) return(1);
	      buffer[buffer_count] = Ufine[k];
	    } /* endfor (k) */
	  } /* endif */
	} /* endfor */
      } /* endfor */

      return 0;

    } /* endif */
  } /* endif */

  // Load send message buffer for the coarse-to-fine grid for cases in
  // which one (or both) of the increments is negative.  Only for two
  // ghost cells.

  if (j_min == j_max) { // North or south boundary.
    // Four different orderings to consider depending on the value of i_inc & j_inc.
    if (j_inc > 0) {
      if (i_inc > 0) {
	return 1;
      } else {
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i, j_min).
	  SubcellReconstruction(i, j_min, Limiter);
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
		Grid.Node[i+1][j_min].X+
		Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	    Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + ( (phi[i][j_min]^dUdx[i][j_min])*dX.x +
				  (phi[i][j_min]^dUdy[i][j_min])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i][j_min].X+
		HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
		HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
		Grid.Cell[i][j_min].Xc)/FOUR - Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + ( (phi[i][j_min]^dUdx[i][j_min])*dX.x +
				  (phi[i][j_min]^dUdy[i][j_min])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	} /* endfor */
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
		Grid.Node[i+1][j_min+1].X)/FOUR - Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + ( (phi[i][j_min]^dUdx[i][j_min])*dX.x +
				  (phi[i][j_min]^dUdy[i][j_min])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
		Grid.Cell[i][j_min].Xc+
		Grid.Node[i][j_min+1].X+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	    Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + ( (phi[i][j_min]^dUdx[i][j_min])*dX.x +
				  (phi[i][j_min]^dUdy[i][j_min])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	} /* endfor */
      } /* endif */
    } else {
      if (i_inc > 0) {
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i, j_min).
	  SubcellReconstruction(i, j_min, Limiter);
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
		Grid.Cell[i][j_min].Xc+
		Grid.Node[i][j_min+1].X+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	    Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + ( (phi[i][j_min]^dUdx[i][j_min])*dX.x +
				  (phi[i][j_min]^dUdy[i][j_min])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
		Grid.Node[i+1][j_min+1].X)/FOUR -
	    Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + ( (phi[i][j_min]^dUdx[i][j_min])*dX.x +
				  (phi[i][j_min]^dUdy[i][j_min])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	} /* endfor */
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i][j_min].X+
		HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
		HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
		Grid.Cell[i][j_min].Xc)/FOUR - Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + ( (phi[i][j_min]^dUdx[i][j_min])*dX.x +
				  (phi[i][j_min]^dUdy[i][j_min])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
		Grid.Node[i+1][j_min].X+
		Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	    Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + ( (phi[i][j_min]^dUdx[i][j_min])*dX.x +
				  (phi[i][j_min]^dUdy[i][j_min])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	} /* endfor */
      } else {
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i, j_min).
	  SubcellReconstruction(i, j_min, Limiter);
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
		Grid.Node[i+1][j_min+1].X)/FOUR - Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + ( (phi[i][j_min]^dUdx[i][j_min])*dX.x +
				  (phi[i][j_min]^dUdy[i][j_min])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
		Grid.Cell[i][j_min].Xc+
		Grid.Node[i][j_min+1].X+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	    Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + ( (phi[i][j_min]^dUdx[i][j_min])*dX.x +
				  (phi[i][j_min]^dUdy[i][j_min])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	} /* endfor */
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
		Grid.Node[i+1][j_min].X+
		Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	    Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + ( (phi[i][j_min]^dUdx[i][j_min])*dX.x +
				  (phi[i][j_min]^dUdy[i][j_min])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i][j_min].X+
		HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
		HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
		Grid.Cell[i][j_min].Xc)/FOUR - Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + ( (phi[i][j_min]^dUdx[i][j_min])*dX.x +
				  (phi[i][j_min]^dUdy[i][j_min])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	} /* endfor */
      } /* endif */
    } /* endif */
  } else { // East or west boundary.
    // Four different orderings to consider depending on the value of i_inc & j_inc.
    if (j_inc > 0) {
      if (i_inc > 0) {
	return 1;
      } else {
	for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i_min, j).
	  SubcellReconstruction(i_min, j, Limiter);
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		Grid.Node[i_min+1][j].X+
		Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
	    Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + ( (phi[i_min][j]^dUdx[i_min][j])*dX.x +
				  (phi[i_min][j]*dUdy[i_min][j])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i_min][j].X+
		HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc)/FOUR - Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + ( (phi[i_min][j]^dUdx[i_min][j])*dX.x +
				  (phi[i_min][j]*dUdy[i_min][j])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
		Grid.Node[i_min+1][j+1].X)/FOUR - Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + ( (phi[i_min][j]^dUdx[i_min][j])*dX.x +
				  (phi[i_min][j]*dUdy[i_min][j])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc+
		Grid.Node[i_min][j+1].X+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
	    Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + ( (phi[i_min][j]^dUdx[i_min][j])*dX.x +
				  (phi[i_min][j]*dUdy[i_min][j])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	} /* endfor */
      } /* endif */
    } else {
      if (i_inc > 0) {
	for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i_min, j).
	  SubcellReconstruction(i_min, j, Limiter);
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc+
		Grid.Node[i_min][j+1].X+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
	    Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + ( (phi[i_min][j]^dUdx[i_min][j])*dX.x +
				  (phi[i_min][j]*dUdy[i_min][j])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
		Grid.Node[i_min+1][j+1].X)/FOUR - Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + ( (phi[i_min][j]^dUdx[i_min][j])*dX.x +
				  (phi[i_min][j]*dUdy[i_min][j])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i_min][j].X+
		HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc)/FOUR - Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + ( (phi[i_min][j]^dUdx[i_min][j])*dX.x +
				  (phi[i_min][j]*dUdy[i_min][j])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		Grid.Node[i_min+1][j].X+
		Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
	    Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + ( (phi[i_min][j]^dUdx[i_min][j])*dX.x +
				  (phi[i_min][j]*dUdy[i_min][j])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	} /* endfor */
      } else {
	for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i_min, j).
	  SubcellReconstruction(i_min, j, Limiter);
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
		Grid.Node[i_min+1][j+1].X)/FOUR - Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + ( (phi[i_min][j]^dUdx[i_min][j])*dX.x +
				  (phi[i_min][j]*dUdy[i_min][j])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc+
		Grid.Node[i_min][j+1].X+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
	    Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + ( (phi[i_min][j]^dUdx[i_min][j])*dX.x +
				  (phi[i_min][j]*dUdy[i_min][j])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		Grid.Node[i_min+1][j].X+
		Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
	    Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + ( (phi[i_min][j]^dUdx[i_min][j])*dX.x +
				  (phi[i_min][j]*dUdy[i_min][j])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i_min][j].X+
		HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc)/FOUR - Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + ( (phi[i_min][j]^dUdx[i_min][j])*dX.x +
				  (phi[i_min][j]*dUdy[i_min][j])*dX.y );
	  for (k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = Ufine[k];
	  } /* endfor (k) */
	} /* endfor */
      } /* endif */
    } /* endif */
  } /* endif */

  return(0);

}

/*******************************************************************************
 * Unloads receive buffer.
 *******************************************************************************/
int AdvectDiffuse2D_Quad_Block_New::UnloadReceiveBuffer(double *buffer,
							int &buffer_count,
							const int buffer_size,
							const int i_min, 
							const int i_max,
							const int i_inc,
							const int j_min, 
							const int j_max,
							const int j_inc) {
  int i, j, k;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
    for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      for ( k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	U[i][j][k] = buffer[buffer_count];
      } /* endfor (k) */
    } /* endfor */
  } /* endfor */
  return(0);
}

/***********************************************************************************
 * Unloads receive message buffer for fine to coarse block message passing.
 ***********************************************************************************/
int AdvectDiffuse2D_Quad_Block_New::UnloadReceiveBuffer_F2C(double *buffer,
							    int &buffer_count,
							    const int buffer_size,
							    const int i_min, 
							    const int i_max,
							    const int i_inc,
							    const int j_min, 
							    const int j_max,
							    const int j_inc) {
  int i, j, k;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
    for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      for ( k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	U[i][j][k] = buffer[buffer_count];
      } /* endfor (k) */
    } /* endfor */
  } /* endfor */
  return(0);
}

/***********************************************************************************
 * Unloads receive message buffer for coarse to fine block message passing.
 ***********************************************************************************/
int AdvectDiffuse2D_Quad_Block_New::UnloadReceiveBuffer_C2F(double *buffer,
							    int &buffer_count,
							    const int buffer_size,
							    const int i_min, 
							    const int i_max,
							    const int i_inc,
							    const int j_min, 
							    const int j_max,
							    const int j_inc) {
  int i, j, k;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
    for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      for ( k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	U[i][j][k] = buffer[buffer_count];
      } /* endfor (k) */
    } /* endfor */
  } /* endfor */
  return(0);
}

/**************************************************************
 * Performs the subcell reconstruction of solution state 
 * within a given cell (i,j) of the computational mesh for 
 * the specified quadrilateral solution block.            
 **************************************************************/
void AdvectDiffuse2D_Quad_Block_New::SubcellReconstruction(const int i, 
							   const int j,
							   const int Limiter) {

  int n, n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4], phi_k;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  double Du, DuDx_ave, DuDy_ave;

  /* Carry out the limited solution reconstruction in
     each cell of the computational mesh. */

  // Determine the number of neighbouring cells to
  // be used in the reconstruction procedure.  Away from
  // boundaries this 8 neighbours will be used.
  if (i <= ICl-2 || i >= ICu+2 || j <= JCl-2 || j >= JCu+2) {
    n_pts = 0;
  } else if ((i == ICl-1) && (Grid.BCtypeW[j] != BC_NONE)) {
    if (j == JCl-1 || j == JCu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeW[j] == BC_PERIODIC ||
               Grid.BCtypeW[j] == BC_NEUMANN ||
               Grid.BCtypeW[j] == BC_ROBIN) {
      if (j == JCl) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i+1; j_index[1] = j  ;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (j == JCu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      } /* endif */
    } else {
      if (j == JCl) {
	n_pts = 3;
	i_index[0] = i+1; j_index[0] = j  ;
	i_index[1] = i  ; j_index[1] = j+1;
	i_index[2] = i+1; j_index[2] = j+1;
      } else if (j == JCu) {
	n_pts = 3;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
      } else {
	n_pts = 5;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } /* endif */
    } /* endif */           
  } else if ((i == ICu+1) && (Grid.BCtypeE[j] != BC_NONE)) {
    if (j == JCl-1 || j == JCu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeE[j] == BC_PERIODIC ||
               Grid.BCtypeE[j] == BC_NEUMANN ||
               Grid.BCtypeE[j] == BC_ROBIN) {
      if (j == JCl) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i+1; j_index[1] = j  ;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (j == JCu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      } /* endif */
    } else {
      if (j == JCl) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i  ; j_index[2] = j+1;
      } else if (j == JCu) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
      } else {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
	i_index[3] = i-1; j_index[3] = j+1;
	i_index[4] = i  ; j_index[4] = j+1;
      } /* endif */
    } /* endif */
  } else if ((j == JCl-1) && (Grid.BCtypeS[i] != BC_NONE)) {
    if (i == ICl-1 || i == ICu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeS[i] == BC_PERIODIC ||
               Grid.BCtypeS[i] == BC_NEUMANN ||
               Grid.BCtypeS[i] == BC_ROBIN) {
      if (i == ICl) {
	n_pts = 5;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (i == ICu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
	i_index[3] = i-1; j_index[3] = j+1;
	i_index[4] = i  ; j_index[4] = j+1;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      } /* endif */
    } else {
      if (i == ICl) {
	n_pts = 3;
	i_index[0] = i+1; j_index[0] = j  ;
	i_index[1] = i  ; j_index[1] = j+1;
	i_index[2] = i+1; j_index[2] = j+1;
      } else if (i == ICu) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i  ; j_index[2] = j+1;
      } else {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j  ;
	i_index[1] = i+1; j_index[1] = j  ;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } /* endif */
    } /* endif */
  } else if ((j == JCu+1) && (Grid.BCtypeN[i] != BC_NONE)) {
    if (i == ICl-1 || i == ICu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeN[i] == BC_PERIODIC ||
               Grid.BCtypeN[i] == BC_NEUMANN ||
               Grid.BCtypeN[i] == BC_ROBIN) {
      if (i == ICl) {
	n_pts = 5;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
	i_index[3] = i  ; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
      } else if (i == ICu) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
	i_index[3] = i-1; j_index[3] = j+1;
	i_index[4] = i  ; j_index[4] = j+1;
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
	i_index[5] = i-1; j_index[5] = j+1;
	i_index[6] = i  ; j_index[6] = j+1;
	i_index[7] = i+1; j_index[7] = j+1;
      } /* endif */
    } else {
      if (i == ICl) {
	n_pts = 3;
	i_index[0] = i  ; j_index[0] = j-1;
	i_index[1] = i+1; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j  ;
      } else if (i == ICu) {
	n_pts = 3;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i-1; j_index[2] = j  ;
      } else {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i  ; j_index[1] = j-1;
	i_index[2] = i+1; j_index[2] = j-1;
	i_index[3] = i-1; j_index[3] = j  ;
	i_index[4] = i+1; j_index[4] = j  ;
      } /* endif */
    } /* endif */
  } else {
    n_pts = 8;
    i_index[0] = i-1; j_index[0] = j-1;
    i_index[1] = i  ; j_index[1] = j-1;
    i_index[2] = i+1; j_index[2] = j-1;
    i_index[3] = i-1; j_index[3] = j  ;
    i_index[4] = i+1; j_index[4] = j  ;
    i_index[5] = i-1; j_index[5] = j+1;
    i_index[6] = i  ; j_index[6] = j+1;
    i_index[7] = i+1; j_index[7] = j+1;
  } /* endif */
  
  // Perform reconstruction.
  if (n_pts > 0) {
    DuDx_ave = ZERO;
    DuDy_ave = ZERO;
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;
  
    for ( n = 0 ; n <= n_pts-1 ; ++n ) {
      dX = Grid.Cell[ i_index[n] ][ j_index[n] ].Xc - Grid.Cell[i][j].Xc;
      Du = U[ i_index[n] ][ j_index[n] ][1] - U[i][j][1];
      DuDx_ave += Du*dX.x;
      DuDy_ave += Du*dX.y;
      DxDx_ave += dX.x*dX.x;
      DxDy_ave += dX.x*dX.y;
      DyDy_ave += dX.y*dX.y;
    } /* endfor */
  					    
    DuDx_ave = DuDx_ave/double(n_pts);
    DuDy_ave = DuDy_ave/double(n_pts);
    DxDx_ave = DxDx_ave/double(n_pts);
    DxDy_ave = DxDy_ave/double(n_pts);
    DyDy_ave = DyDy_ave/double(n_pts);
    dUdx[i][j][1] = ( (DuDx_ave*DyDy_ave-DuDy_ave*DxDy_ave)/
		      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
    dUdy[i][j][1] = ( (DuDy_ave*DxDx_ave-DuDx_ave*DxDy_ave)/
		      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave) );
  
    // Calculate slope limiter.
    if (!Freeze_Limiter) {
      u0Min = U[i][j][1];
      u0Max = u0Min;
      for ( n = 0 ; n <= n_pts-1 ; ++n ) {
	u0Min = min(u0Min, U[ i_index[n] ][ j_index[n] ][1]);
	u0Max = max(u0Max, U[ i_index[n] ][ j_index[n] ][1]);
      } /* endfor */
  
      dX = Grid.xfaceE(i, j)-Grid.Cell[i][j].Xc;
      uQuad[0] = U[i][j][1] + dUdx[i][j][1]*dX.x + dUdy[i][j][1]*dX.y;
      dX = Grid.xfaceW(i, j)-Grid.Cell[i][j].Xc;
      uQuad[1] = U[i][j][1] + dUdx[i][j][1]*dX.x + dUdy[i][j][1]*dX.y;
      dX = Grid.xfaceN(i, j)-Grid.Cell[i][j].Xc;
      uQuad[2] = U[i][j][1] + dUdx[i][j][1]*dX.x + dUdy[i][j][1]*dX.y ;
      dX = Grid.xfaceS(i, j)-Grid.Cell[i][j].Xc;
      uQuad[3] = U[i][j][1] + dUdx[i][j][1]*dX.x + dUdy[i][j][1]*dX.y ;
  
      switch(Limiter) {
      case LIMITER_ONE :
	phi_k = ONE;
	break;
      case LIMITER_ZERO :
	phi_k = ZERO;
	break;
      case LIMITER_BARTH_JESPERSEN :
	phi_k = Limiter_BarthJespersen(uQuad, U[i][j][1], u0Min, u0Max, 4);
	break;
      case LIMITER_VENKATAKRISHNAN :
	phi_k = Limiter_Venkatakrishnan(uQuad, U[i][j][1], u0Min, u0Max, 4);
	break;
      case LIMITER_VANLEER :
	phi_k = Limiter_VanLeer(uQuad, U[i][j][1], u0Min, u0Max, 4);
	break;
      case LIMITER_VANALBADA :
	phi_k = Limiter_VanAlbada(uQuad, U[i][j][1], u0Min, u0Max, 4);
	break;
      default:
	phi_k = Limiter_BarthJespersen(uQuad, U[i][j][1], u0Min, u0Max, 4);
	break;
      } /* endswitch */
  
      phi[i][j][1] = phi_k;
    } /* endif */
  } else {
    dUdx[i][j].Vacuum();
    dUdy[i][j].Vacuum();
    phi[i][j].Vacuum();
  } /* endif */

}

/***************************************************************
 * Loads send message buffer for fine to coarse block message 
 * passing of conservative solution fluxes.
 ****************************************************************/
int AdvectDiffuse2D_Quad_Block_New::LoadSendBuffer_Flux_F2C(double *buffer,
							    int &buffer_count,
							    const int buffer_size,
							    const int i_min, 
							    const int i_max,
							    const int i_inc,
							    const int j_min, 
							    const int j_max,
							    const int j_inc) {
  int i, j, k;
  if (j_min == j_max && j_min == JCl) {
    for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
      for ( k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = (FluxS[i  ][k] + FluxS[i+1][k]);
      }
    } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
    for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
      for ( k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = (FluxN[i  ][k] + FluxN[i+1][k]);
      }
    } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
    for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
      for ( k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = (FluxW[j][k] + FluxW[j+1][k]);
      }
    } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
    for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
      for ( k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = (FluxE[j][k] + FluxE[j+1][k]);
      }
    } /* endfor */
  } /* endif */
  return(0);
}

/***********************************************************
 * Unloads receive message buffer for fine to coarse 
 * block message passing of conservative solution fluxes.  
 ***********************************************************/
int AdvectDiffuse2D_Quad_Block_New::UnloadReceiveBuffer_Flux_F2C(double *buffer,
								 int &buffer_count,
								 const int buffer_size,
								 const int i_min, 
								 const int i_max,
								 const int i_inc,
								 const int j_min, 
								 const int j_max,
								 const int j_inc) {
  int i, j, k;
  if (j_min == j_max && j_min == JCl) {
    for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      for ( k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	FluxS[i][k] = -buffer[buffer_count] - FluxS[i][k];
      }
    } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
    for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      for ( k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	buffer_count++;
	if (buffer_count >= buffer_size) return(1);
	FluxN[i][k] = -buffer[buffer_count] - FluxN[i][k];
      }
    } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
    for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
      for ( k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	buffer_count++;
	if (buffer_count >= buffer_size) return(1);
	FluxW[j][k] = -buffer[buffer_count] - FluxW[j][k];
      }
    } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
    for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
      for ( k = 1; k <= NUM_VAR_ADVECTDIFFUSE2D; ++k) {
	buffer_count++;
	if (buffer_count >= buffer_size) return(1);
	FluxE[j][k] = -buffer[buffer_count] - FluxE[j][k];
      }
    } /* endfor */
  } /* endif */
  return(0);
}
