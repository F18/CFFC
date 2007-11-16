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
// Initialize ExactGrad
AdvectDiffuse2D_Quad_Block_New::Exact_Gradient_Function AdvectDiffuse2D_Quad_Block_New::ExactGrad = NULL;
// Initialize ExactSoln
FunctionType2D AdvectDiffuse2D_Quad_Block_New::ExactSoln = NULL;


/*******************************************************************************
 * AdvectDiffuse2D_Quad_Block -- Single Block Member Functions.                *
 ******************************************************************************/
/**********************
 * Default constructor.
 **********************/
AdvectDiffuse2D_Quad_Block_New::AdvectDiffuse2D_Quad_Block_New(void) {
  NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
  U = NULL; dt = NULL; dudt = NULL; 
  dudx = NULL; dudy = NULL; phi = NULL; uo = NULL;
  FluxN = NULL; FluxS = NULL; FluxE = NULL; FluxW = NULL;
  UoN = NULL; UoS = NULL; UoE = NULL; UoW = NULL;
  Axisymmetric = 0; Freeze_Limiter = OFF;
}

/******************************************
 * Private copy constructor. (shallow copy)
 *****************************************/
AdvectDiffuse2D_Quad_Block_New::AdvectDiffuse2D_Quad_Block_New(const AdvectDiffuse2D_Quad_Block_New &Soln) {
  NCi = Soln.NCi; ICl = Soln.ICl; ICu = Soln.ICu; 
  NCj = Soln.NCj; JCl = Soln.JCl; JCu = Soln.JCu; Nghost = Soln.Nghost;
  Grid = Soln.Grid; U = Soln.U; dt = Soln.dt; dudt = Soln.dudt; 
  dudx = Soln.dudx; dudy = Soln.dudy; phi = Soln.phi;
  uo = Soln.uo;
  FluxN = Soln.FluxN; FluxS = Soln.FluxS; FluxE = Soln.FluxE; FluxW = Soln.FluxW;
  UoN = Soln.UoN; UoS = Soln.UoS; UoE = Soln.UoE; UoW = Soln.UoW;
  Axisymmetric = Soln.Axisymmetric; Freeze_Limiter = Soln.Freeze_Limiter;
}

/**********************
 * Allocate memory.            
 **********************/
void AdvectDiffuse2D_Quad_Block_New::allocate(const int Ni, const int Nj, const int Ng) {
  int i, j, k; assert(Ni > 1 && Nj > 1 && Ng > 1 && Ng > 1); Grid.allocate(Ni, Nj, Ng);
  NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1;
  NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1; Nghost = Ng;
  U = new AdvectDiffuse2D_State_New*[NCi]; dt = new double*[NCi]; dudt = new double**[NCi]; 
  dudx = new double*[NCi]; dudy = new double*[NCi]; 
  phi = new double*[NCi]; uo = new double*[NCi];
  for ( i = 0; i <= NCi-1 ; ++i ) {
    U[i] = new AdvectDiffuse2D_State_New[NCj]; 
    dt[i] = new double[NCj]; dudt[i] = new double*[NCj];
    for ( j = 0; j <= NCj-1 ; ++j ) 
      { dudt[i][j] = new double[NUMBER_OF_RESIDUAL_VECTORS_ADVECTDIFFUSE2D]; }
    dudx[i] = new double[NCj]; dudy[i] = new double[NCj]; 
    phi[i] = new double[NCj];
    uo[i] = new double[NCj];
  } /* endfor */
  FluxN = new double[NCi]; FluxS = new double[NCi];
  FluxE = new double[NCj]; FluxW = new double[NCj];
  UoN = new AdvectDiffuse2D_State_New[NCi]; UoS = new AdvectDiffuse2D_State_New[NCi];
  UoE = new AdvectDiffuse2D_State_New[NCj]; UoW = new AdvectDiffuse2D_State_New[NCj];
  // Set the solution residuals, gradients, limiters, and other values to zero.
  for (j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
    for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
      for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_ADVECTDIFFUSE2D-1 ; ++k ) { dudt[i][j][k] = ZERO; }
      dudx[i][j] = ZERO; dudy[i][j] = ZERO; phi[i][j] = ZERO; uo[i][j] = ZERO; dt[i][j] = ZERO;
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
    for ( j = 0; j <= NCj-1 ; ++j ) { delete []dudt[i][j]; dudt[i][j] = NULL; }
    delete []dudt[i]; dudt[i] = NULL;
    delete []dudx[i]; dudx[i] = NULL; delete []dudy[i]; dudy[i] = NULL;
    delete []phi[i]; phi[i] = NULL; delete []uo[i]; uo[i] = NULL;
  } /* endfor */
  delete []U; U = NULL; delete []dt; dt = NULL; delete []dudt; dudt = NULL;
  delete []dudx; dudx = NULL; delete []dudy; dudy = NULL; 
  delete []phi; phi = NULL; delete []uo; uo = NULL;
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

/***********************
 * Node solution value. 
 ***********************/
double AdvectDiffuse2D_Quad_Block_New::un(const int ii, const int jj) {
  return Un(ii,jj)[1];
}

/********************
 * Output operator.
 ********************/
ostream &operator << (ostream &out_file,
		      const AdvectDiffuse2D_Quad_Block_New &SolnBlk) {
  int i, j; 
  out_file << SolnBlk.Grid;
  out_file << SolnBlk.NCi << " " << SolnBlk.ICl << " " 
           << SolnBlk.ICu << " " << SolnBlk.Nghost << "\n";
  out_file << SolnBlk.NCj << " " << SolnBlk.JCl << " " << SolnBlk.JCu << "\n";
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
  in_file >> SolnBlk.Grid; 
  in_file.setf(ios::skipws);
  in_file >> ni >> il >> iu >> ng; in_file >> nj >> jl >> ju;
  in_file >> SolnBlk.Axisymmetric;
  in_file.unsetf(ios::skipws);
  if (ni == 0 || nj == 0) {
    SolnBlk.deallocate(); return(in_file);
  } /* endif */
  if (SolnBlk.U == NULL || SolnBlk.NCi != ni || SolnBlk.NCj != nj || SolnBlk.Nghost != ng) {
    if (SolnBlk.U != NULL) SolnBlk.deallocate(); 
    SolnBlk.allocate(ni - 2*ng, nj - 2*ng, ng);
  } /* endif */
  for ( j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
    for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
      in_file >> SolnBlk.U[i][j];
      for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_ADVECTDIFFUSE2D-1 ; ++k ) {
	SolnBlk.dudt[i][j][k] = ZERO;
      } /* endfor */
      SolnBlk.dudx[i][j] = ZERO;
      SolnBlk.dudy[i][j] = ZERO;
      SolnBlk.phi[i][j] = ZERO;
      SolnBlk.uo[i][j] = ZERO;
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
      buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
      if (buffer_count >= buffer_size) return(1);
      buffer[buffer_count] = U[i][j].u;
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
      buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
      if (buffer_count >= buffer_size) return(1);
      buffer[buffer_count] = (Grid.Cell[i  ][j  ].A*U[i  ][j  ].u+
			      Grid.Cell[i+1][j  ].A*U[i+1][j  ].u+
			      Grid.Cell[i  ][j+1].A*U[i  ][j+1].u+
			      Grid.Cell[i+1][j+1].A*U[i+1][j+1].u)/
	(Grid.Cell[i  ][j  ].A+
	 Grid.Cell[i+1][j  ].A+
	 Grid.Cell[i  ][j+1].A+
	 Grid.Cell[i+1][j+1].A);
      /*         buffer[buffer_count] = (Grid.Cell[i  ][j  ].A*U[i  ][j  ].u+ */
      /*                                 Grid.Cell[i+1][j  ].A*U[i+1][j  ].u+ */
      /*                                 Grid.Cell[i  ][j+1].A*U[i  ][j+1].u+ */
      /*                                 Grid.Cell[i+1][j+1].A*U[i+1][j+1].u); */
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
  int i, j;
  Vector2D dX;
  double ufine;

  if (j_inc > 0) {
    if (i_inc > 0) {
      for ( j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i, j_min).
	  SubcellReconstruction(i, j, LIMITER_VENKATAKRISHNAN);
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
	    ufine = U[i][j].u + (phi[i][j]*dudx[i][j])*dX.x +
	      (phi[i][j]*dudy[i][j])*dX.y;
	    buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = ufine;
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
	    ufine = U[i][j].u + (phi[i][j]*dudx[i][j])*dX.x +
	      (phi[i][j]*dudy[i][j])*dX.y;
	    buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = ufine;
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
	    ufine = U[i][j].u + (phi[i][j]*dudx[i][j])*dX.x +
	      (phi[i][j]*dudy[i][j])*dX.y;
	    buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = ufine;
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
	    ufine = U[i][j].u + (phi[i][j]*dudx[i][j])*dX.x +
	      (phi[i][j]*dudy[i][j])*dX.y;
	    buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	    if (buffer_count >= buffer_size) return(1);
	    buffer[buffer_count] = ufine;
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
	  SubcellReconstruction(i, j_min, LIMITER_VENKATAKRISHNAN);
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
		Grid.Node[i+1][j_min].X+
		Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	    Grid.Cell[i][j_min].Xc;
	  ufine = U[i][j_min].u +
	    (phi[i][j_min]*dudx[i][j_min])*dX.x +
	    (phi[i][j_min]*dudy[i][j_min])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i][j_min].X+
		HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
		HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
		Grid.Cell[i][j_min].Xc)/FOUR -
	    Grid.Cell[i][j_min].Xc;
	  ufine = U[i][j_min].u +
	    (phi[i][j_min]*dudx[i][j_min])*dX.x +
	    (phi[i][j_min]*dudy[i][j_min])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	} /* endfor */
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
		Grid.Node[i+1][j_min+1].X)/FOUR -
	    Grid.Cell[i][j_min].Xc;
	  ufine = U[i][j_min].u +
	    (phi[i][j_min]*dudx[i][j_min])*dX.x +
	    (phi[i][j_min]*dudy[i][j_min])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
		Grid.Cell[i][j_min].Xc+
		Grid.Node[i][j_min+1].X+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	    Grid.Cell[i][j_min].Xc;
	  ufine = U[i][j_min].u +
	    (phi[i][j_min]*dudx[i][j_min])*dX.x +
	    (phi[i][j_min]*dudy[i][j_min])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	} /* endfor */
      } /* endif */
    } else {
      if (i_inc > 0) {
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i, j_min).
	  SubcellReconstruction(i, j_min, LIMITER_VENKATAKRISHNAN);
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
		Grid.Cell[i][j_min].Xc+
		Grid.Node[i][j_min+1].X+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	    Grid.Cell[i][j_min].Xc;
	  ufine = U[i][j_min].u +
	    (phi[i][j_min]*dudx[i][j_min])*dX.x +
	    (phi[i][j_min]*dudy[i][j_min])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
		Grid.Node[i+1][j_min+1].X)/FOUR -
	    Grid.Cell[i][j_min].Xc;
	  ufine = U[i][j_min].u +
	    (phi[i][j_min]*dudx[i][j_min])*dX.x +
	    (phi[i][j_min]*dudy[i][j_min])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	} /* endfor */
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i][j_min].X+
		HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
		HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
		Grid.Cell[i][j_min].Xc)/FOUR -
	    Grid.Cell[i][j_min].Xc;
	  ufine = U[i][j_min].u +
	    (phi[i][j_min]*dudx[i][j_min])*dX.x +
	    (phi[i][j_min]*dudy[i][j_min])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
		Grid.Node[i+1][j_min].X+
		Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	    Grid.Cell[i][j_min].Xc;
	  ufine = U[i][j_min].u +
	    (phi[i][j_min]*dudx[i][j_min])*dX.x +
	    (phi[i][j_min]*dudy[i][j_min])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	} /* endfor */
      } else {
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i, j_min).
	  SubcellReconstruction(i, j_min, LIMITER_VENKATAKRISHNAN);
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
		Grid.Node[i+1][j_min+1].X)/FOUR -
	    Grid.Cell[i][j_min].Xc;
	  ufine = U[i][j_min].u +
	    (phi[i][j_min]*dudx[i][j_min])*dX.x +
	    (phi[i][j_min]*dudy[i][j_min])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
		Grid.Cell[i][j_min].Xc+
		Grid.Node[i][j_min+1].X+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	    Grid.Cell[i][j_min].Xc;
	  ufine = U[i][j_min].u +
	    (phi[i][j_min]*dudx[i][j_min])*dX.x +
	    (phi[i][j_min]*dudy[i][j_min])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	} /* endfor */
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
		Grid.Node[i+1][j_min].X+
		Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	    Grid.Cell[i][j_min].Xc;
	  ufine = U[i][j_min].u +
	    (phi[i][j_min]*dudx[i][j_min])*dX.x +
	    (phi[i][j_min]*dudy[i][j_min])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i][j_min].X+
		HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
		HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
		Grid.Cell[i][j_min].Xc)/FOUR -
	    Grid.Cell[i][j_min].Xc;
	  ufine = U[i][j_min].u +
	    (phi[i][j_min]*dudx[i][j_min])*dX.x +
	    (phi[i][j_min]*dudy[i][j_min])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
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
	  SubcellReconstruction(i_min, j, LIMITER_VENKATAKRISHNAN);
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		Grid.Node[i_min+1][j].X+
		Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
	    Grid.Cell[i_min][j].Xc;
	  ufine = U[i_min][j].u +
	    (phi[i_min][j]*dudx[i_min][j])*dX.x +
	    (phi[i_min][j]*dudy[i_min][j])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i_min][j].X+
		HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc)/FOUR -
	    Grid.Cell[i_min][j].Xc;
	  ufine = U[i_min][j].u +
	    (phi[i_min][j]*dudx[i_min][j])*dX.x +
	    (phi[i_min][j]*dudy[i_min][j])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
		Grid.Node[i_min+1][j+1].X)/FOUR -
	    Grid.Cell[i_min][j].Xc;
	  ufine = U[i_min][j].u +
	    (phi[i_min][j]*dudx[i_min][j])*dX.x +
	    (phi[i_min][j]*dudy[i_min][j])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc+
		Grid.Node[i_min][j+1].X+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
	    Grid.Cell[i_min][j].Xc;
	  ufine = U[i_min][j].u +
	    (phi[i_min][j]*dudx[i_min][j])*dX.x +
	    (phi[i_min][j]*dudy[i_min][j])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	} /* endfor */
      } /* endif */
    } else {
      if (i_inc > 0) {
	for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i_min, j).
	  SubcellReconstruction(i_min, j, LIMITER_VENKATAKRISHNAN);
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc+
		Grid.Node[i_min][j+1].X+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
	    Grid.Cell[i_min][j].Xc;
	  ufine = U[i_min][j].u +
	    (phi[i_min][j]*dudx[i_min][j])*dX.x +
	    (phi[i_min][j]*dudy[i_min][j])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
		Grid.Node[i_min+1][j+1].X)/FOUR -
	    Grid.Cell[i_min][j].Xc;
	  ufine = U[i_min][j].u +
	    (phi[i_min][j]*dudx[i_min][j])*dX.x +
	    (phi[i_min][j]*dudy[i_min][j])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i_min][j].X+
		HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc)/FOUR -
	    Grid.Cell[i_min][j].Xc;
	  ufine = U[i_min][j].u +
	    (phi[i_min][j]*dudx[i_min][j])*dX.x +
	    (phi[i_min][j]*dudy[i_min][j])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		Grid.Node[i_min+1][j].X+
		Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
	    Grid.Cell[i_min][j].Xc;
	  ufine = U[i_min][j].u +
	    (phi[i_min][j]*dudx[i_min][j])*dX.x +
	    (phi[i_min][j]*dudy[i_min][j])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	} /* endfor */
      } else {
	for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i_min, j).
	  SubcellReconstruction(i_min, j, LIMITER_VENKATAKRISHNAN);
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
		Grid.Node[i_min+1][j+1].X)/FOUR -
	    Grid.Cell[i_min][j].Xc;
	  ufine = U[i_min][j].u +
	    (phi[i_min][j]*dudx[i_min][j])*dX.x +
	    (phi[i_min][j]*dudy[i_min][j])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc+
		Grid.Node[i_min][j+1].X+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
	    Grid.Cell[i_min][j].Xc;
	  ufine = U[i_min][j].u +
	    (phi[i_min][j]*dudx[i_min][j])*dX.x +
	    (phi[i_min][j]*dudy[i_min][j])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		Grid.Node[i_min+1][j].X+
		Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
	    Grid.Cell[i_min][j].Xc;
	  ufine = U[i_min][j].u +
	    (phi[i_min][j]*dudx[i_min][j])*dX.x +
	    (phi[i_min][j]*dudy[i_min][j])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i_min][j].X+
		HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc)/FOUR -
	    Grid.Cell[i_min][j].Xc;
	  ufine = U[i_min][j].u +
	    (phi[i_min][j]*dudx[i_min][j])*dX.x +
	    (phi[i_min][j]*dudy[i_min][j])*dX.y;
	  buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
	  if (buffer_count >= buffer_size) return(1);
	  buffer[buffer_count] = ufine;
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
  int i, j;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
    for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
      if (buffer_count >= buffer_size) return(1);
      U[i][j].u = buffer[buffer_count];
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
  int i, j;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
    for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
      if (buffer_count >= buffer_size) return(1);
      U[i][j].u = buffer[buffer_count];
      /*         U[i][j].u = buffer[buffer_count]/Grid.Cell[i][j].A; */
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
  int i, j;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
    for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
      if (buffer_count >= buffer_size) return(1);
      U[i][j].u = buffer[buffer_count];
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
  if (i == ICl-Nghost || i == ICu+Nghost ||
      j == JCl-Nghost || j == JCu+Nghost) {
    n_pts = 0;
  } else if ((i == ICl-Nghost+1) && 
             (Grid.BCtypeW[j] != BC_NONE)) {
    if (j == JCl-Nghost+1 || j == JCu+Nghost-1) {
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
  } else if ((i == ICu+Nghost-1) && 
             (Grid.BCtypeE[j] != BC_NONE)) {
    if (j == JCl-Nghost+1 || j == JCu+Nghost-1) {
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
  } else if ((j == JCl-Nghost+1) && 
             (Grid.BCtypeS[i] != BC_NONE)) {
    if (i == ICl-Nghost+1 || i == ICu+Nghost-1) {
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
  } else if ((j == JCu+Nghost-1) && 
             (Grid.BCtypeN[i] != BC_NONE)) {
    if (i == ICl-Nghost+1 || i == ICu+Nghost-1) {
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
      dX = Grid.Cell[ i_index[n] ][ j_index[n] ].Xc - 
	Grid.Cell[i][j].Xc;
      Du = U[ i_index[n] ][ j_index[n] ].u - 
	U[i][j].u;
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
    dudx[i][j] = (DuDx_ave*DyDy_ave-DuDy_ave*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    dudy[i][j] = (DuDy_ave*DxDx_ave-DuDx_ave*DxDy_ave)/
      (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
  
    // Calculate slope limiter.
    if (!Freeze_Limiter) {
      u0Min = U[i][j].u;
      u0Max = u0Min;
      for ( n = 0 ; n <= n_pts-1 ; ++n ) {
	u0Min = min(u0Min, U[ i_index[n] ][ j_index[n] ].u);
	u0Max = max(u0Max, U[ i_index[n] ][ j_index[n] ].u);
      } /* endfor */
  
      dX = Grid.xfaceE(i, j)-Grid.Cell[i][j].Xc;
      uQuad[0] = U[i][j].u + 
	dudx[i][j]*dX.x +
	dudy[i][j]*dX.y ;
      dX = Grid.xfaceW(i, j)-Grid.Cell[i][j].Xc;
      uQuad[1] = U[i][j].u + 
	dudx[i][j]*dX.x +
	dudy[i][j]*dX.y ;
      dX = Grid.xfaceN(i, j)-Grid.Cell[i][j].Xc;
      uQuad[2] = U[i][j].u + 
	dudx[i][j]*dX.x +
	dudy[i][j]*dX.y ;
      dX = Grid.xfaceS(i, j)-Grid.Cell[i][j].Xc;
      uQuad[3] = U[i][j].u + 
	dudx[i][j]*dX.x +
	dudy[i][j]*dX.y ;
  
      switch(Limiter) {
      case LIMITER_ONE :
	phi_k = ONE;
	break;
      case LIMITER_ZERO :
	phi_k = ZERO;
	break;
      case LIMITER_BARTH_JESPERSEN :
	phi_k = Limiter_BarthJespersen(uQuad, U[i][j].u, 
				       u0Min, u0Max, 4);
	break;
      case LIMITER_VENKATAKRISHNAN :
	phi_k = Limiter_Venkatakrishnan(uQuad, U[i][j].u, 
					u0Min, u0Max, 4);
	break;
      case LIMITER_VANLEER :
	phi_k = Limiter_VanLeer(uQuad, U[i][j].u, 
				u0Min, u0Max, 4);
	break;
      case LIMITER_VANALBADA :
	phi_k = Limiter_VanAlbada(uQuad, U[i][j].u, 
				  u0Min, u0Max, 4);
	break;
      default:
	phi_k = Limiter_BarthJespersen(uQuad, U[i][j].u, 
				       u0Min, u0Max, 4);
	break;
      } /* endswitch */
  
      phi[i][j] = phi_k;
    } /* endif */
  } else {
    dudx[i][j] = ZERO;
    dudy[i][j] = ZERO; 
    phi[i][j]  = ZERO;
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
  int i, j;
  if (j_min == j_max && j_min == JCl) {
    for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
      buffer_count = buffer_count + 1;
      if (buffer_count >= buffer_size) return(1);
      buffer[buffer_count] = (FluxS[i  ]+
			      FluxS[i+1]);
    } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
    for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
      buffer_count = buffer_count + 1;
      if (buffer_count >= buffer_size) return(1);
      buffer[buffer_count] = (FluxN[i  ]+
			      FluxN[i+1]);
    } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
    for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
      buffer_count = buffer_count + 1;
      if (buffer_count >= buffer_size) return(1);
      buffer[buffer_count] = (FluxW[j]+
			      FluxW[j+1]);
    } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
    for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
      buffer_count = buffer_count + 1;
      if (buffer_count >= buffer_size) return(1);
      buffer[buffer_count] = (FluxE[j]+
			      FluxE[j+1]);
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
  int i, j;
  if (j_min == j_max && j_min == JCl) {
    for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      buffer_count = buffer_count + 1;
      if (buffer_count >= buffer_size) return(1);
      FluxS[i] = -buffer[buffer_count]
	-FluxS[i];
    } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
    for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      buffer_count = buffer_count + 1;
      if (buffer_count >= buffer_size) return(1);
      FluxN[i] = -buffer[buffer_count]
	-FluxN[i];
    } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
    for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
      buffer_count = buffer_count + 1;
      if (buffer_count >= buffer_size) return(1);
      FluxW[j] = -buffer[buffer_count]
	-FluxW[j];
    } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
    for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
      buffer_count = buffer_count + 1;
      if (buffer_count >= buffer_size) return(1);
      FluxE[j] = -buffer[buffer_count]
	-FluxE[j];
    } /* endfor */
  } /* endif */
  return(0);
}
