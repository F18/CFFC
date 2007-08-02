/**********************************************************************
 * Electrostatic2DQuad.h: Header file defining 2D Electrostatic       *
 *                        quadrilateral mesh solution classes.        *
 **********************************************************************/

#ifndef _ELECTROSTATIC2D_QUAD_INCLUDED
#define _ELECTROSTATIC2D_QUAD_INCLUDED

// Include 2D electrostatic state header file.

#ifndef _ELECTROSTATIC2D_STATE_INCLUDED
#include "Electrostatic2DState.h"
#endif // _ELECTROSTATIC2D_STATE_INCLUDED

// Include 2D electrostatic input header file.

#ifndef _ELECTROSTATIC2D_INPUT_INCLUDED
#include "Electrostatic2DInput.h"
#endif // _ELECTROSTATIC2D_INPUT_INCLUDED

// Include the 2D cell header file.

#ifndef _CELL2D_INCLUDED
#include "../Grid/Cell2D.h"
#endif // _CELL2D_INCLUDED

// Include the 2D adjusted quadrilateral multiblock grid header file.

#ifndef _GRID2D_QUAD_BLOCK_INCLUDED
#include "../Grid/Grid2DQuad.h"
#endif // _GRID2D_QUAD_BLOCK_INCLUDED

// Include the quadtree header file.

#ifndef _QUADTREE_INCLUDED
#include "../AMR/QuadTree.h"
#endif // _QUADTREE_INCLUDED

// Include the AMR header file.

#ifndef _AMR_INCLUDED
#include "../AMR/AMR.h"
#endif // _AMR_INCLUDED

// Include the linear systems header file.

#ifndef _LINEARSYSTEMS_INCLUDED
#include "../Math/LinearSystems.h"
#endif // _LINEARSYSTEMS_INCLUDED

// Include ICEMCFD input header file.

#ifndef _ICEMCFD_INCLUDED
#include "../ICEM/ICEMCFD.h"
#endif // _ICEMCFD_INCLUDED

#define	NUMBER_OF_RESIDUAL_VECTORS_ELECTROSTATIC2D  3

/*!
 * Class: Electrostatic2D_Quad_Block
 *
 * @brief Class definition of the 2D electrostatic solution blocks.
 *
 * \verbatim
 * Member functions
 *       W      -- Return primitive variable solution for the block
 *                 (cell average).
 *       U      -- Return conserved variable solution for the block
 *                 (cell average).
 *    Grid      -- Return the solution block quadrilateral grid or mesh.
 *      dt      -- Return local time step for the solution block.
 *    dUdt      -- Return the solution block residuals.
 *    dUdx      -- Return the unlimited primitive variable solution
 *                 gradients (x-direction) for the block.
 *    dUdy      -- Return the unlimited primitive variable solution
 *                 gradients (y-direction) for the block.
 *      Uo      -- Return initial solution state.
 *     NCi      -- Return number of cells in the i-direction 
 *                 (zeta-direction).
 *     ICl      -- Return lower index for cells in the i-direction 
 *                 (zeta-direction).
 *     ICu      -- Return upper index for cells in the i-direction 
 *                 (zeta-direction).
 *     NCj      -- Return number of cells in the j-direction 
 *                 (eta-direction).
 *     JCl      -- Return lower index for cells in the j-direction 
 *                 (eta-direction).
 *     JCu      -- Return upper index for cells in the j-direction 
 *                 (eta-direction).
 *  Nghost      -- Number of ghost cells.
 * Axisymmetric -- Return axisymmetric flow indicator (= 1 for 
 *                 axisymmetric flow, = 0 for planar flow).
 *     UoN      -- Return array of reference states for the
 *                 application of boundary conditions at the north
 *                 boundary of the solution block.
 *     UoS      -- Return array of reference states for the
 *                 application of boundary conditions at the south
 *                 boundary of the solution block.
 *     UoE      -- Return array of reference states for the
 *                 application of boundary conditions at the east
 *                 boundary of the solution block.
 *     UoW      -- Return array of reference states for the
 *                 application of boundary conditions at the west
 *                 boundary of the solution block.
 *   allocate   -- Allocate memory for structured quadrilateral
 *                 solution block.
 *   deallocate -- Deallocate memory for structured quadrilateral
 *                 solution block.
 *      Wn      -- Return primitive variable solution at the
 *                 specified node.
 *      Un      -- Return conserved variable solution at the
 *                 specified node.
 *    WnNW      -- Return primitive variable solution at the
 *                 north-west node.
 *    WnNE      -- Return primitive variable solution at the
 *                 north-east node.
 *    WnSW      -- Return primitive variable solution at the
 *                 south-west node.
 *    WnSE      -- Return primitive variable solution at the
 *                 south-east node.
 *    UnNW      -- Return conserved variable solution at the
 *                 north-west node.
 *    UnNE      -- Return conserved variable solution at the
 *                 north-east node.
 *    UnSW      -- Return conserved variable solution at the
 *                 south-west node.
 *    UnSE      -- Return conserved variable solution at the
 *                 south-east node.
 * evaluate_limiters -- Set flag to evaluate limiters.
 * freeze_limiters -- Set flag to freeze limiters.
 * Member functions required for message passing.
 *    NumVar    -- Returns number of solution variables in primitive
 *                 and conserved solution state vectors.
 * LoadSendBuffer -- Loads send buffer.
 * LoadSendBuffer_F2C -- Loads send buffer for fine to coarse block 
 *                       messages.
 * UnloadReceiveBuffer -- Unloads receive buffer.
 * UnloadReceiveBuffer_F2C -- Unloads receive buffer for fine to
 *                            coarse block messages.
 * LinearSubcellReconstruction -- Performs linear subcell solution
 *                                reconstruction used in adaptive mesh
 *                                refinement.
 * LoadSendBuffer_Flux_F2C -- Loads send buffer for sending
 *                            conservative flux corrections from fine
 *                            to coarse solution blocks.
 * UnLoadSendBuffer_Flux_F2C -- Loads send buffer for sending
 *                              conservative flux corrections from
 *                              fine to coarse solution blocks.
 *
 * Member operators
 *      S -- a 2D Electrostatic solution
 *
 * S = S;
 * cout << S; (output function)
 * cin  >> S; (input function)
 * \endverbatim
 */
class Electrostatic2D_Quad_Block{
private:
public:
  //@{ @name Solution state arrays:
  Electrostatic2DState      **U; //!< Solution state.
  //@}

  //@{ @name Grid block information:
  int                        NCi, //!< Total number of i-direction cells.
                             ICl, //!< First i-direction non-ghost cell counter.
                             ICu; //!< Final i-direction non-ghost cell counter.
  int                        NCj, //!< Total number of j-direction cells.
                             JCl, //!< First j-direction non-ghost cell counter.
                             JCu; //!< Final j-direction non-ghost cell counter.
  int                     Nghost; //!< Number of ghost cells.
  Grid2D_Quad_Block         Grid; //!< 2D quadrilateral grid geometry.
  //@}

  //@{ @name Residual and time-stepping arrays:
  double                    **dt; //!< Local time step.
  Electrostatic2DState   ***dUdt; //!< Solution residual.
  Electrostatic2DState      **Uo; //!< Initial solution state.
  static int   residual_variable; //!< Static integer that indicates which variable is used for residual calculations.
  //@}

  //@{ @name Solution gradient arrays:
  Electrostatic2DState    **dUdx, //!< Unlimited solution gradient (x-direction).
                          **dUdy; //!< Unlimited solution gradient (y-direction).
  Electrostatic2DState **ddUdxdx, //!< Unlimited solution gradient (x-direction).
                       **ddUdxdy, //!< Unlimited solution gradient (x-direction).
                       **ddUdydy; //!< Unlimited solution gradient (x-direction).
  //@}

  //@{ @name Boundary solution flux arrays:
  Electrostatic2DState    *FluxN, //!< North boundary solution flux.
                          *FluxS, //!< South boundary solution flux.
                          *FluxE, //!< East boundary solution flux.
                          *FluxW; //!< West boundary solution flux.
  //@}

  //@{ @name Problem indicator flags:
  int               Axisymmetric; //!< Axisymmetric flow indicator.
  static int           Flow_Type; //!< A dummy static integer indicating the 'flow-type.'
  //@}

  //@{ @name Boundary condtion reference states:
  Electrostatic2DState      *UoN, //!< Boundary condition reference states for north boundary.
                            *UoS, //!< Boundary condition reference states for south boundary.
                            *UoE, //!< Boundary condition reference states for east boundary.
                            *UoW; //!< Boundary condition reference states for west boundary.
  //@}

  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  Electrostatic2D_Quad_Block(void) {
    // Problem flags:
    Axisymmetric = OFF;
    // Grid size and variables:
    NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
    // Solution variables:
    U = NULL;
    dt = NULL; dUdt = NULL; Uo = NULL;
    dUdx = NULL; dUdy = NULL;
    ddUdxdx = NULL; ddUdxdy = NULL; ddUdydy = NULL;
    FluxN = NULL; FluxS = NULL; FluxE = NULL; FluxW = NULL;
    UoN = NULL; UoS = NULL; UoE = NULL; UoW = NULL;
  }

  //! Copy constructor.
  Electrostatic2D_Quad_Block(const Electrostatic2D_Quad_Block &Soln) {
    // Problem flags:
    Axisymmetric = Soln.Axisymmetric;
    // Grid size and variables:
    NCi = Soln.NCi; ICl = Soln.ICl; ICu = Soln.ICu;
    NCj = Soln.NCj; JCl = Soln.JCl; JCu = Soln.JCu;
    Nghost = Soln.Nghost;
    Grid  = Soln.Grid;
    // Solution variables:
    U = Soln.U;
    dt = Soln.dt; dUdt = Soln.dUdt; Uo = Soln.Uo;
    dUdx = Soln.dUdx; dUdy = Soln.dUdy;
    ddUdxdx = Soln.ddUdxdx; ddUdxdy = Soln.ddUdxdy; ddUdydy = Soln.ddUdydy;
    FluxN = Soln.FluxN; FluxS = Soln.FluxS; 
    FluxE = Soln.FluxE; FluxW = Soln.FluxW;
    UoN   = Soln.UoN;   UoS   = Soln.UoS;
    UoE   = Soln.UoE;   UoW   = Soln.UoW;
  }

  // Destructor.
  // ~Electrostatic2D_Quad_Block(void);
  // Use automatically generated destructor.
  //@}

  // Assignment operator.
  // Electrostatic2D_Quad_Block operator = (const Electrostatic2D_Quad_Block &Soln);
  // Use automatically generated assignment operator.

  //@{ @name Allocate and deallocate functions.
  //! Allocate memory for structured quadrilateral solution block.
  void allocate(const int Ni, const int Nj, const int Ng);
  //! Deallocate memory for structured quadrilateral solution block.
  void deallocate(void);
  //@}

  //@{ @name Bilinear interplation (Zingg & Yarrow).
  //! Retern solution state at specified node.
  Electrostatic2DState Uzy(const int ii, const int jj);

  Electrostatic2DState UzyNW(const int ii, const int jj); //!< Return solution state at cell NW node.
  Electrostatic2DState UzyNE(const int ii, const int jj); //!< Return solution state at cell NE node.
  Electrostatic2DState UzySE(const int ii, const int jj); //!< Return solution state at cell SE node.
  Electrostatic2DState UzySW(const int ii, const int jj); //!< Return solution state at cell SW node.
  //@}

  //@{ @name Bilinear interplation (Connell & Holmes).
  //! Retern conserverd solution state at specified node.
  Electrostatic2DState Un(const int ii, const int jj);

  //! Return solution state at cell nodes.
  Electrostatic2DState UnNW(const int ii, const int jj); //!< Return solution state at cell NW node.
  Electrostatic2DState UnNE(const int ii, const int jj); //!< Return solution state at cell NE node.
  Electrostatic2DState UnSE(const int ii, const int jj); //!< Return solution state at cell SE node.
  Electrostatic2DState UnSW(const int ii, const int jj); //!< Return solution state at cell SW node.
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const Electrostatic2D_Quad_Block &Soln);
  friend istream &operator >> (istream &in_file, Electrostatic2D_Quad_Block &Soln);
  //@}

  //@{ @name Member functions required for message passing.
  //! Number of solution state variables.
  int NumVar(void);
  //! Load send message passing buffer.
  int LoadSendBuffer(double *buffer,
		     int &buffer_count,
		     const int buffer_size,
		     const int i_min,
		     const int i_max,
		     const int i_inc,
		     const int j_min,
		     const int j_max,
		     const int j_inc);
  //! Load F2C send message passing buffer.
  int LoadSendBuffer_F2C(double *buffer,
			 int &buffer_count,
			 const int buffer_size,
			 const int i_min,
			 const int i_max,
			 const int i_inc,
			 const int j_min,
			 const int j_max,
			 const int j_inc);
  //! Load C2F send message passing buffer.
  int LoadSendBuffer_C2F(double *buffer,
			 int &buffer_count,
			 const int buffer_size,
			 const int i_min,
			 const int i_max,
			 const int i_inc,
			 const int j_min,
			 const int j_max,
			 const int j_inc,
			 const int face,
			 const int sector);
  //! Unload receive message passing buffer.
  int UnloadReceiveBuffer(double *buffer,
			  int &buffer_count,
			  const int buffer_size,
			  const int i_min,
			  const int i_max,
			  const int i_inc,
			  const int j_min,
			  const int j_max,
			  const int j_inc);
  //! Unload F2C receive message passing buffer.
  int UnloadReceiveBuffer_F2C(double *buffer,
			      int &buffer_count,
			      const int buffer_size,
			      const int i_min,
			      const int i_max,
			      const int i_inc,
			      const int j_min,
			      const int j_max,
			      const int j_inc);
  //! Unload C2F receive message passing buffer.
  int UnloadReceiveBuffer_C2F(double *buffer,
			      int &buffer_count,
			      const int buffer_size,
			      const int i_min,
			      const int i_max,
			      const int i_inc,
			      const int j_min,
			      const int j_max,
			      const int j_inc);
  //! Linear subcell solution reconstruction within given computational cell.
  void LinearSubcellReconstruction(const int i,
				   const int j);
  //! Load F2C conservative flux message passing buffer.
  int LoadSendBuffer_Flux_F2C(double *buffer,
			      int &buffer_count,
			      const int buffer_size,
			      const int i_min,
			      const int i_max,
			      const int i_inc,
			      const int j_min,
			      const int j_max,
			      const int j_inc);
  //! Unload F2C conservative flux message passing buffer.
  int UnloadReceiveBuffer_Flux_F2C(double *buffer,
				   int &buffer_count,
				   const int buffer_size,
				   const int i_min,
				   const int i_max,
				   const int i_inc,
				   const int j_min,
				   const int j_max,
				   const int j_inc);
  //@}

};

/**********************************************************************
 * Electrostatic2D_Quad_Block::allocate -- Allocate memory.           *
 **********************************************************************/
inline void Electrostatic2D_Quad_Block::allocate(const int Ni, const int Nj, const int Ng) {
  assert(Ni > 1 && Nj > 1 && Ng > 1);
  Grid.allocate(Ni,Nj,Ng);
  NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1; 
  NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1; Nghost = Ng;
  U = new Electrostatic2DState*[NCi];
  dt = new double*[NCi]; dUdt = new Electrostatic2DState**[NCi];
  dUdx = new Electrostatic2DState*[NCi];
  dUdy = new Electrostatic2DState*[NCi];
  ddUdxdx = new Electrostatic2DState*[NCi];
  ddUdxdy = new Electrostatic2DState*[NCi];
  ddUdydy = new Electrostatic2DState*[NCi];
  Uo = new Electrostatic2DState*[NCi];
  for (int i = 0; i < NCi; i++) {
    U[i] = new Electrostatic2DState[NCj];
    dt[i] = new double[NCj];
    dUdt[i] = new Electrostatic2DState*[NCj];
    for (int j = 0; j < NCj; j++)
      dUdt[i][j] = new Electrostatic2DState[NUMBER_OF_RESIDUAL_VECTORS_ELECTROSTATIC2D];
    dUdx[i] = new Electrostatic2DState[NCj];
    dUdy[i] = new Electrostatic2DState[NCj];
    ddUdxdx[i] = new Electrostatic2DState[NCj];
    ddUdxdy[i] = new Electrostatic2DState[NCj];
    ddUdydy[i] = new Electrostatic2DState[NCj];
    Uo[i] = new Electrostatic2DState[NCj];
  }
  UoN = new Electrostatic2DState[NCi]; UoS = new Electrostatic2DState[NCi];
  UoE = new Electrostatic2DState[NCj]; UoW = new Electrostatic2DState[NCj];
  FluxN = new Electrostatic2DState[NCi]; FluxS = new Electrostatic2DState[NCi];
  FluxE = new Electrostatic2DState[NCj]; FluxW = new Electrostatic2DState[NCj];
  // Set the solution residuals, gradients, limiters, and other values to zero.
  for (int j = JCl-Nghost; j <= JCu+Nghost; j++) {
    for (int i = ICl-Nghost; i <= ICu+Nghost; i++) {
      for (int k = 0; k < NUMBER_OF_RESIDUAL_VECTORS_ELECTROSTATIC2D; k++)
	dUdt[i][j][k].Vacuum();
      dUdx[i][j].Vacuum(); dUdy[i][j].Vacuum();
      ddUdxdx[i][j].Vacuum();
      ddUdxdy[i][j].Vacuum();
      ddUdydy[i][j].Vacuum();
      Uo[i][j].Vacuum();
      dt[i][j] = ZERO;
    }
  }
}

/**********************************************************************
 * Electrostatic2D_Quad_Block::deallocate -- Deallocate memory.       *
 **********************************************************************/
inline void Electrostatic2D_Quad_Block::deallocate(void) {
  Grid.deallocate();
  for (int i = 0; i < NCi; i++) {
    delete []U[i]; U[i] = NULL;
    delete []dt[i]; dt[i] = NULL; 
    for (int j = 0; j < NCj; j++) { delete []dUdt[i][j]; dUdt[i][j] = NULL; }
    delete []dUdt[i]; dUdt[i] = NULL;
    delete []dUdx[i]; dUdx[i] = NULL; delete []dUdy[i]; dUdy[i] = NULL;
    delete []ddUdxdx[i]; ddUdxdx[i] = NULL;
    delete []ddUdxdy[i]; ddUdxdy[i] = NULL;
    delete []ddUdydy[i]; ddUdydy[i] = NULL;
    delete []Uo[i]; Uo[i] = NULL;
  }
  delete []U; U = NULL;
  delete []dt; dt = NULL; delete []dUdt; dUdt = NULL;
  delete []dUdx; dUdx = NULL; delete []dUdy; dUdy = NULL; 
  delete []ddUdxdx; ddUdxdx = NULL;
  delete []ddUdxdy; ddUdxdy = NULL;
  delete []ddUdydy; ddUdydy = NULL;
  delete []Uo; Uo = NULL;
  delete []FluxN; FluxN = NULL; delete []FluxS; FluxS = NULL;
  delete []FluxE; FluxE = NULL; delete []FluxW; FluxW = NULL;
  delete []UoN; UoN = NULL; delete []UoS; UoS = NULL;
  delete []UoE; UoE = NULL; delete []UoW; UoW = NULL;
  NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
}

/**********************************************************************
 * Electrostatic2D_Quad_Block::Uzy -- Node variable solution state.   *
 *                                                                    *
 * Zingg and Yarrow (SIAM J. Sci. Stat. Comput. Vol. 13 No. 3 1992)   *
 *                                                                    *
 **********************************************************************/
inline Electrostatic2DState Electrostatic2D_Quad_Block::Uzy(const int ii, const int jj) {
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
         eta1, zeta1, eta2, zeta2, eta, zeta;
  Electrostatic2DState A, B, C, D;
  x  = Grid.Node[ii][jj].X.x;
  y  = Grid.Node[ii][jj].X.y;
  ax = Grid.Cell[ii-1][jj-1].Xc.x;
  bx = Grid.Cell[ii-1][jj  ].Xc.x - Grid.Cell[ii-1][jj-1].Xc.x; 
  cx = Grid.Cell[ii  ][jj-1].Xc.x - Grid.Cell[ii-1][jj-1].Xc.x; 
  dx = Grid.Cell[ii  ][jj  ].Xc.x + Grid.Cell[ii-1][jj-1].Xc.x -
       Grid.Cell[ii-1][jj  ].Xc.x - Grid.Cell[ii  ][jj-1].Xc.x;
  ay = Grid.Cell[ii-1][jj-1].Xc.y;
  by = Grid.Cell[ii-1][jj  ].Xc.y - Grid.Cell[ii-1][jj-1].Xc.y; 
  cy = Grid.Cell[ii  ][jj-1].Xc.y - Grid.Cell[ii-1][jj-1].Xc.y; 
  dy = Grid.Cell[ii  ][jj  ].Xc.y + Grid.Cell[ii-1][jj-1].Xc.y -
       Grid.Cell[ii-1][jj  ].Xc.y - Grid.Cell[ii  ][jj-1].Xc.y;
  aa = bx*dy - dx*by;
  bb = dy*(ax-x) + bx*cy - cx*by + dx*(y-ay);
  cc = cy*(ax-x) + cx*(y-ay);
  if (fabs(aa) < TOLER*TOLER) {
    if (fabs(bb) >= TOLER*TOLER) { zeta1=-cc/bb; }

    else { zeta1 = -cc/sgn(bb)*(TOLER*TOLER); } 
    if (fabs(cy+dy*zeta1) >= TOLER*TOLER) { eta1=(y-ay-by*zeta1)/(cy+dy*zeta1); } 
    else { eta1 = HALF; } zeta2=zeta1; eta2=eta1;
  } else {
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) { 
      zeta1 = HALF*(-bb+sqrt(bb*bb-FOUR*aa*cc))/aa; 
    } else { 
      zeta1 = -HALF*bb/aa; } 
    if (fabs(cy+dy*zeta1) < TOLER*TOLER) {
      eta1 = -ONE;
    } else { 
      eta1 = (y-ay-by*zeta1)/(cy+dy*zeta1); 
    } 
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) { 
      zeta2 = HALF*(-bb-sqrt(bb*bb-FOUR*aa*cc))/aa;
    } else {
      zeta2 = -HALF*bb/aa; 
    }
    if (fabs(cy+dy*zeta2) < TOLER*TOLER) {
      eta2 = -ONE; 
    } else {
      eta2 = (y-ay-by*zeta2)/(cy+dy*zeta2);
    }
  }
  if (zeta1 > -TOLER && zeta1 < ONE + TOLER &&
      eta1  > -TOLER && eta1  < ONE + TOLER) {
    zeta = zeta1;
    eta  = eta1;
  } else if (zeta2 > -TOLER && zeta2 < ONE + TOLER &&
	     eta2  > -TOLER && eta2  < ONE + TOLER) {
    zeta = zeta2;
    eta  = eta2;
  } else {
    zeta = HALF;
    eta  = HALF;
  }
  A = U[ii-1][jj-1];
  B = U[ii-1][jj  ] - U[ii-1][jj-1]; 
  C = U[ii  ][jj-1] - U[ii-1][jj-1];
  D = U[ii  ][jj  ] + U[ii-1][jj-1] - U[ii-1][jj  ] - U[ii  ][jj-1];
  return A+B*zeta+C*eta+D*zeta*eta;
}

/**********************************************************************
 * Electrostatic2D_Quad_Block::Uzy?? -- Cell node solution states.    *
 **********************************************************************/
inline Electrostatic2DState Electrostatic2D_Quad_Block::UzyNW(const int ii, const int jj) {
  return Uzy(ii,jj+1);
}

inline Electrostatic2DState Electrostatic2D_Quad_Block::UzyNE(const int ii, const int jj) {
  return Uzy(ii+1,jj+1);
}

inline Electrostatic2DState Electrostatic2D_Quad_Block::UzySE(const int ii, const int jj) {
  return Uzy(ii+1,jj);
}

inline Electrostatic2DState Electrostatic2D_Quad_Block::UzySW(const int ii, const int jj) {
  return Uzy(ii,jj);
}

/**********************************************************************
 * Electrostatic2D_Quad_Block::Un -- Node variable solution state.    *
 *                                                                    *
 * Holmes and Connell (AIAA Paper 1989-1932-CP)                       *
 *                                                                    *
 **********************************************************************/
inline Electrostatic2DState Electrostatic2D_Quad_Block::Un(const int ii, const int jj) {
  double w1, w2, w3, w4;
  Vector2D lambda, R, X0, X1, X2, X3, X4;
  Tensor2D I;
  Electrostatic2DState U1, U2, U3, U4;
  // Summarize cell-centres and states.
  X0 = Grid.Node[ii][jj].X;
  X1 = Grid.Cell[ii-1][jj-1].Xc; U1 = U[ii-1][jj-1];
  X2 = Grid.Cell[ii  ][jj-1].Xc; U2 = U[ii  ][jj-1];
  X3 = Grid.Cell[ii-1][jj  ].Xc; U3 = U[ii-1][jj  ];
  X4 = Grid.Cell[ii  ][jj  ].Xc; U4 = U[ii  ][jj  ];
  // Determine weighting coefficients:
  R = (X1 - X0) + (X2 - X0) + (X3 - X0) + (X4 - X0);
  I.xx = (X1.x - X0.x)*(X1.x - X0.x) + (X2.x - X0.x)*(X2.x - X0.x) +
         (X3.x - X0.x)*(X3.x - X0.x) + (X4.x - X0.x)*(X4.x - X0.x);
  I.xy = (X1.x - X0.x)*(X1.y - X0.y) + (X2.x - X0.x)*(X2.y - X0.y) +
         (X3.x - X0.x)*(X3.y - X0.y) + (X4.x - X0.x)*(X4.y - X0.y);
  I.yy = (X1.y - X0.y)*(X1.y - X0.y) + (X2.y - X0.y)*(X2.y - X0.y) +
         (X3.y - X0.y)*(X3.y - X0.y) + (X4.y - X0.y)*(X4.y - X0.y);
  lambda.x = (I.xy*R.y - I.yy*R.x)/(I.xx*I.yy - I.xy*I.xy);
  lambda.y = (I.xy*R.x - I.xx*R.y)/(I.xx*I.yy - I.xy*I.xy);
  // Determine the weights:
  w1 = 1 + lambda.x*(X1.x - X0.x) + lambda.y*(X1.y - X0.y);
  w2 = 1 + lambda.x*(X2.x - X0.x) + lambda.y*(X2.y - X0.y);
  w3 = 1 + lambda.x*(X3.x - X0.x) + lambda.y*(X3.y - X0.y);
  w4 = 1 + lambda.x*(X4.x - X0.x) + lambda.y*(X4.y - X0.y);
  // Determine the interpolated state:
  return (w1*U1 + w2*U2 + w3*U3+ w4*U4)/(w1 + w2 + w3 + w4);
}

/**********************************************************************
 * Electrostatic2D_Quad_Block::Un?? -- Cell node solution states.     *
 **********************************************************************/
inline Electrostatic2DState Electrostatic2D_Quad_Block::UnNW(const int ii, const int jj) {
  return Un(ii,jj+1);
}

inline Electrostatic2DState Electrostatic2D_Quad_Block::UnNE(const int ii, const int jj) {
  return Un(ii+1,jj+1);
}

inline Electrostatic2DState Electrostatic2D_Quad_Block::UnSE(const int ii, const int jj) {
  return Un(ii+1,jj);
}

inline Electrostatic2DState Electrostatic2D_Quad_Block::UnSW(const int ii, const int jj) {
  return Un(ii,jj);
}

/**********************************************************************
 * Electrostatic2D_Quad_Block -- Input-output operators.              *
 **********************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Electrostatic2D_Quad_Block &SolnBlk) {
  out_file << SolnBlk.Grid;
  out_file << SolnBlk.NCi << " " << SolnBlk.ICl << " " << SolnBlk.ICu << endl;
  out_file << SolnBlk.NCj << " " << SolnBlk.JCl << " " << SolnBlk.JCu << endl;
  out_file << SolnBlk.Nghost << endl;
  out_file << SolnBlk.Axisymmetric << endl;
  if (SolnBlk.NCi == 0 || SolnBlk.NCj == 0) return out_file;
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      out_file << SolnBlk.U[i][j];
      out_file << endl;
    }
  }
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    out_file << SolnBlk.UoW[j] << endl;
    out_file << SolnBlk.UoE[j] << endl;
  }
  for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
    out_file << SolnBlk.UoS[i] << endl;
    out_file << SolnBlk.UoN[i] << endl;
  }
  return out_file;
}

inline istream &operator >> (istream &in_file,
			     Electrostatic2D_Quad_Block &SolnBlk) {
  int ni, il, iu, nj, jl, ju, ng;
  Grid2D_Quad_Block New_Grid;
  in_file >> New_Grid;
  in_file.setf(ios::skipws);
  in_file >> ni >> il >> iu;
  in_file >> nj >> jl >> ju;
  in_file >> ng;
  in_file >> SolnBlk.Axisymmetric;
  in_file.unsetf(ios::skipws);
  if (ni == 0 || nj == 0) {
    SolnBlk.deallocate(); return in_file;
  }
  if (SolnBlk.U == NULL || SolnBlk.NCi != ni || SolnBlk.NCj != nj || SolnBlk.Nghost != ng) {
    if (SolnBlk.U != NULL) SolnBlk.deallocate();
    SolnBlk.allocate(ni-2*ng,nj-2*ng,ng);
  }
  Copy_Quad_Block(SolnBlk.Grid,New_Grid); New_Grid.deallocate();
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      in_file >> SolnBlk.U[i][j];
      for (int k = 0; k < NUMBER_OF_RESIDUAL_VECTORS_ELECTROSTATIC2D; k++) {
	SolnBlk.dUdt[i][j][k].Vacuum();
      }
      SolnBlk.dUdx[i][j].Vacuum();
      SolnBlk.dUdy[i][j].Vacuum();
      SolnBlk.ddUdxdx[i][j].Vacuum();
      SolnBlk.ddUdxdy[i][j].Vacuum();
      SolnBlk.ddUdydy[i][j].Vacuum();
      SolnBlk.Uo[i][j].Vacuum();
      SolnBlk.dt[i][j] = ZERO;
    }
  }
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    in_file >> SolnBlk.UoW[j];
    in_file >> SolnBlk.UoE[j];
  }
  for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
    in_file >> SolnBlk.UoS[i];
    in_file >> SolnBlk.UoN[i];
  }
  return in_file;
}

/**********************************************************************
 *                                                                    *
 * MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING.                     *
 *                                                                    *
 **********************************************************************/

/**********************************************************************
 * Electrostatic2D_Quad_Block::NumVar -- Returns number of state      *
 *                                       variables.                   *
 **********************************************************************/
inline int Electrostatic2D_Quad_Block::NumVar(void) {
  return NUM_VAR_ELECTROSTATIC2D;
}

/**********************************************************************
 * Electrostatic2D_Quad_Block::LoadSendBuffer -- Loads send message   *
 *                                               buffer.              *
 **********************************************************************/
inline int Electrostatic2D_Quad_Block::LoadSendBuffer(double *buffer,
						      int &buffer_count,
						      const int buffer_size,
						      const int i_min,
						      const int i_max,
						      const int i_inc,
						      const int j_min,
						      const int j_max,
						      const int j_inc) {

  for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
    for (int i = i_min;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = U[i][j][k];
      }
    }
  }
  return 0;
}

/**********************************************************************
 * Electrostatic2D_Quad_Block::LoadSendBuffer_F2C --                  *
 *                     Loads send message buffer for fine to coarse   *
 *                     block message passing.                         *
 **********************************************************************/
inline int Electrostatic2D_Quad_Block::LoadSendBuffer_F2C(double *buffer,
							  int &buffer_count,
							  const int buffer_size,
							  const int i_min,
							  const int i_max,
							  const int i_inc,
							  const int j_min,
							  const int j_max,
							  const int j_inc) {

  for (int j = j_min; ((j_inc+2)/4) ? (j < j_max):(j > j_max); j += j_inc) {
    for (int i = i_min; ((i_inc+2)/4) ? (i < i_max):(i > i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = (Grid.Cell[i  ][j  ].A*U[i  ][j  ][k]+
				Grid.Cell[i+1][j  ].A*U[i+1][j  ][k]+
				Grid.Cell[i  ][j+1].A*U[i  ][j+1][k]+
				Grid.Cell[i+1][j+1].A*U[i+1][j+1][k])/
	                       (Grid.Cell[i  ][j  ].A+
				Grid.Cell[i+1][j  ].A+
				Grid.Cell[i  ][j+1].A+
				Grid.Cell[i+1][j+1].A);
      }
    }
  }
  return 0;
}

/**********************************************************************
 * Electrostatic2D_Quad_Block::LoadSendBuffer_C2F --                  *
 *                     Loads send message buffer for coarse to fine   *
 *                     block message passing.                         *
 **********************************************************************/
inline int Electrostatic2D_Quad_Block::LoadSendBuffer_C2F(double *buffer,
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

  Vector2D dX;
  Electrostatic2DState Ufine;

  if (j_inc > 0) {
    if (i_inc > 0) {
      for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
	  // Perform limited linear least squares reconstruction in cell (i,j_min).
 	  LinearSubcellReconstruction(i,j);
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
	    Ufine = U[i][j] + dUdx[i][j]*dX.x + dUdy[i][j]*dX.y;
	    for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	      buffer_count++;
	      if (buffer_count >= buffer_size) return 1;
	      buffer[buffer_count] = Ufine[k];
	    }
	  }
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
	    Ufine = U[i][j] + dUdx[i][j]*dX.x + dUdy[i][j]*dX.y;
	    for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	      buffer_count++;
	      if (buffer_count >= buffer_size) return 1;
	      buffer[buffer_count] = Ufine[k];
	    }
	  }
	}
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
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
	    Ufine = U[i][j] + dUdx[i][j]*dX.x + dUdy[i][j]*dX.y;
	    for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) { 
	      buffer_count++;
	      if (buffer_count >= buffer_size) return 1;
	      buffer[buffer_count] = Ufine[k];
	    }
	  }
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
	    Ufine = U[i][j] + dUdx[i][j]*dX.x + dUdy[i][j]*dX.y;
	    for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	      buffer_count++;
	      if (buffer_count >= buffer_size) return 1;
	      buffer[buffer_count] = Ufine[k];
	    }
	  }
	}
      }

      return 0;

    }
  }

  // Load send message buffer for the coarse-to-fine grid for cases in
  // which one (or both) of the increments is negative.  Only for two
  // ghost cells.

  if (j_min == j_max) { // North or south boundary.
    // Four different orderings to consider depending on the value of i_inc & j_inc.
    if (j_inc > 0) {
      if (i_inc > 0) {
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
	  // Perform limited linear least squares reconstruction in cell (i,j_min).
	  LinearSubcellReconstruction(i,j_min);
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i][j_min].X +
		HALF*(Grid.Node[i][j_min].X + Grid.Node[i+1][j_min].X)+
		HALF*(Grid.Node[i][j_min].X + Grid.Node[i][j_min+1].X)+
		Grid.Cell[i][j_min].Xc)/FOUR - Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + (dUdx[i][j_min])*dX.x +
	                        (dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i  ][j_min].X + Grid.Node[i+1][j_min  ].X) +
		Grid.Node[i+1][j_min].X +
		Grid.Cell[i  ][j_min].Xc +
		HALF*(Grid.Node[i+1][j_min].X + Grid.Node[i+1][j_min+1].X))/FOUR -
	        Grid.Cell[i  ][j_min].Xc;
	  Ufine = U[i][j_min] + (dUdx[i][j_min])*dX.x +
	                        (dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
		Grid.Cell[i][j_min].Xc +
		Grid.Node[i][j_min+1].X +
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	        Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + (dUdx[i][j_min])*dX.x +
	                        (dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) { 
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
		Grid.Node[i+1][j_min+1].X)/FOUR -
	        Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + (dUdx[i][j_min])*dX.x +
	                        (dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
      } else {
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
	  // Perform limited linear least squares reconstruction in cell (i,j_min).
	  LinearSubcellReconstruction(i,j_min);
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
		Grid.Node[i+1][j_min].X+
		Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	        Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + (dUdx[i][j_min])*dX.x +
	                        (dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i][j_min].X+
		HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min  ].X)+
		HALF*(Grid.Node[i][j_min].X+Grid.Node[i  ][j_min+1].X)+
		Grid.Cell[i][j_min].Xc)/FOUR -
	        Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + (dUdx[i][j_min])*dX.x +
	                        (dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min  ].X+Grid.Node[i+1][j_min+1].X)+
		HALF*(Grid.Node[i  ][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
		Grid.Node[i+1][j_min+1].X)/FOUR -
	        Grid.Cell[i  ][j_min  ].Xc;
	  Ufine = U[i][j_min] + (dUdx[i][j_min])*dX.x +
	                        (dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
		Grid.Cell[i][j_min].Xc+
		Grid.Node[i][j_min+1].X+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	        Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + (dUdx[i][j_min])*dX.x +
	                        (dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) { 
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
      }
    } else {
      if (i_inc > 0) {
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
	  // Perform limited linear least squares reconstruction in cell (i,j_min).
	  LinearSubcellReconstruction(i,j_min);
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
		Grid.Cell[i][j_min].Xc+
		Grid.Node[i][j_min+1].X+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	        Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + (dUdx[i][j_min])*dX.x +
	                        (dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) { 
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
		Grid.Node[i+1][j_min+1].X)/FOUR -
	        Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + (dUdx[i][j_min])*dX.x +
	                        (dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i][j_min].X+
		HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
		HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
		Grid.Cell[i][j_min].Xc)/FOUR -
	        Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + (dUdx[i][j_min])*dX.x +
	                        (dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
		Grid.Node[i+1][j_min].X+
		Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	        Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + (dUdx[i][j_min])*dX.x +
	                        (dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
      } else {
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
	  // Perform limited linear least squares reconstruction in cell (i,j_min).
	  LinearSubcellReconstruction(i,j_min);
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
		Grid.Node[i+1][j_min+1].X)/FOUR -
	        Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + (dUdx[i][j_min])*dX.x +
	                        (dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
	        Grid.Cell[i][j_min  ].Xc+
	        Grid.Node[i][j_min+1].X+
	        HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	        Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + (dUdx[i][j_min])*dX.x +
	                        (dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) { 
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
		Grid.Node[i+1][j_min].X+
		Grid.Cell[i  ][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	        Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + (dUdx[i][j_min])*dX.x +
	                        (dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i][j_min].X+
		HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
		HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
		Grid.Cell[i][j_min].Xc)/FOUR -
	        Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + (dUdx[i][j_min])*dX.x +
  	                        (dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
      }
    }
  } else { // East or west boundary.
    // Four different orderings to consider depending on the value of i_inc & j_inc.
    if (j_inc > 0) {
      if (i_inc > 0) {
	for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
	  // Perform limited linear least squares reconstruction in cell (i_min,j).
	  LinearSubcellReconstruction(i_min,j);
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i_min][j].X+
		HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc)/FOUR -
	        Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + (dUdx[i_min][j])*dX.x +
                                (dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		Grid.Node[i_min+1][j].X+
		Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
	        Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + (dUdx[i_min][j])*dX.x +
                                (dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc+
		Grid.Node[i_min][j+1].X+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
	        Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + (dUdx[i_min][j])*dX.x +
	                        (dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) { 
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
		Grid.Node[i_min+1][j+1].X)/FOUR -
	        Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + (dUdx[i_min][j])*dX.x +
	                        (dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
      } else {
	for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
	  // Perform limited linear least squares reconstruction in cell (i_min,j).
	  LinearSubcellReconstruction(i_min,j);
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		Grid.Node[i_min+1][j].X+
		Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
	        Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + (dUdx[i_min][j])*dX.x +
	                        (dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i_min][j].X+
	        HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
	        HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
	        Grid.Cell[i_min][j].Xc)/FOUR -
                Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + (dUdx[i_min][j])*dX.x +
                                (dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
		Grid.Node[i_min+1][j+1].X)/FOUR -
                Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + (dUdx[i_min][j])*dX.x +
	                        (dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc+
		Grid.Node[i_min][j+1].X+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + (dUdx[i_min][j])*dX.x +
                                (dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) { 
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
      }
    } else {
      if (i_inc > 0) {
	for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
	  // Perform limited linear least squares reconstruction in cell (i_min,j).
	  LinearSubcellReconstruction(i_min,j);
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc+
		Grid.Node[i_min][j+1].X+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + (dUdx[i_min][j])*dX.x +
                                (dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) { 
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
		Grid.Node[i_min+1][j+1].X)/FOUR -
                Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + (dUdx[i_min][j])*dX.x +
	                        (dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i_min][j].X+
		HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc)/FOUR -
	        Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + (dUdx[i_min][j])*dX.x +
	                        (dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		Grid.Node[i_min+1][j].X+
		Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + (dUdx[i_min][j])*dX.x +
	                        (dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
      } else {
	for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
	  // Perform limited linear least squares reconstruction in cell (i_min,j).
	  LinearSubcellReconstruction(i_min,j);
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i_min][j].Xc+
	        HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
		Grid.Node[i_min+1][j+1].X)/FOUR -
	        Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + (dUdx[i_min][j])*dX.x +
	                        (dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc+
		Grid.Node[i_min][j+1].X+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + (dUdx[i_min][j])*dX.x +
                                (dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) { 
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		Grid.Node[i_min+1][j].X+
		Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
	        Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + (dUdx[i_min][j])*dX.x +
                                (dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i_min][j].X+
		HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc)/FOUR -
                Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + (dUdx[i_min][j])*dX.x +
                                (dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
      }
    }
  }
  return 0;
}

/**********************************************************************
 * Electrostatic2D_Quad_Block::UnloadReceiveBuffer --                 *
 *                     Unloads receive message buffer.                *
 **********************************************************************/
inline int Electrostatic2D_Quad_Block::UnloadReceiveBuffer(double *buffer,
							   int &buffer_count,
							   const int buffer_size,
							   const int i_min,
							   const int i_max,
							   const int i_inc,
							   const int j_min,
							   const int j_max,
							   const int j_inc) {

  for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
    for (int i = i_min;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	U[i][j][k] = buffer[buffer_count];
      }
    }
  }
  return 0;
}

/**********************************************************************
 * Electrostatic2D_Quad_Block::UnloadReceiveBuffer_F2C --             *
 *                     Unloads receive message buffer for fine to     *
 *                     coarse block message passing.                  *
 **********************************************************************/
inline int Electrostatic2D_Quad_Block::UnloadReceiveBuffer_F2C(double *buffer,
							       int &buffer_count,
							       const int buffer_size,
							       const int i_min,
							       const int i_max,
							       const int i_inc,
							       const int j_min,
							       const int j_max,
							       const int j_inc) {

  for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
    for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;    
	U[i][j][k] = buffer[buffer_count];
      }
    }
  }
  return 0;
}

/**********************************************************************
 * Electrostatic2D_Quad_Block::UnloadReceiveBuffer_C2F --             *
 *                     Unloads receive message buffer for coarse to   *
 *                     fine block message passing.                    *
 **********************************************************************/
inline int Electrostatic2D_Quad_Block::UnloadReceiveBuffer_C2F(double *buffer,
							       int &buffer_count,
							       const int buffer_size,
							       const int i_min,
							       const int i_max,
							       const int i_inc,
							       const int j_min,
							       const int j_max,
							       const int j_inc) {

  for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
    for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	U[i][j][k] = buffer[buffer_count];
      }
    }
  }
  return 0;
}

/**********************************************************************
 * Electrostatic2D_Quad_Block::LinearSubcellReconstruction --         *
 *                     Performs the subcell linear reconstruction of  *
 *                     solution state within a given cell (i,j) of    *
 *                     the computational mesh for the specified       *
 *                     quadrilateral solution block.                  *
 **********************************************************************/
inline void Electrostatic2D_Quad_Block::LinearSubcellReconstruction(const int i,
								    const int j) {

  int n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4];
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  Electrostatic2DState DU, DUDx_ave, DUDy_ave;

  // Carry out the limited solution reconstruction in each cell of the 
  // computational mesh.

  if (i == ICl-2 || i == ICu+2 || j == JCl-2 || j == JCu+2) {
    n_pts = 0;
  } else if (i == ICl-1 && Grid.BCtypeW[j] != BC_NONE) {
    if (j == JCl-1 || j == JCu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeW[j] == BC_PERIODIC ||
               Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
               Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
               Grid.BCtypeW[j] == BC_CHARACTERISTIC) {
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
      }
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
      }
    }
  } else if (i == ICu+1 && Grid.BCtypeE[j] != BC_NONE) {
    if (j == JCl-1 || j == JCu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeE[j] == BC_PERIODIC ||
               Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
               Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
               Grid.BCtypeE[j] == BC_CHARACTERISTIC) {
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
      }
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
      }
    }
  } else if (j == JCl-1 && Grid.BCtypeS[i] != BC_NONE) {
    if (i == ICl-1 || i == ICu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeS[i] == BC_PERIODIC ||
               Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
               Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
               Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
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
      }
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
      }
    }
  } else if (j == JCu+1 && Grid.BCtypeN[i] != BC_NONE) {
    if (i == ICl-1 || i == ICu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeN[i] == BC_PERIODIC ||
               Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
               Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
               Grid.BCtypeN[i] == BC_CHARACTERISTIC) {
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
      }
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
      }
    }
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
  }

  if (n_pts > 0) {
    DUDx_ave.Vacuum();
    DUDy_ave.Vacuum();
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;

    for (int n2 = 0; n2 < n_pts; n2++) {
      dX = Grid.Cell[i_index[n2]][j_index[n2]].Xc - Grid.Cell[i][j].Xc;
      DU = U[i_index[n2]][j_index[n2]] - U[i][j];
      DUDx_ave += DU*dX.x;
      DUDy_ave += DU*dX.y;
      DxDx_ave += dX.x*dX.x;
      DxDy_ave += dX.x*dX.y;
      DyDy_ave += dX.y*dX.y;
    }

    DUDx_ave /= double(n_pts);
    DUDy_ave /= double(n_pts);
    DxDx_ave /= double(n_pts);
    DxDy_ave /= double(n_pts);
    DyDy_ave /= double(n_pts);
    dUdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                 (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    dUdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                 (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);

  } else {
    dUdx[i][j].Vacuum();
    dUdy[i][j].Vacuum(); 

  }

}

/*******************************************************************************
 * Electrostatic2D_Quad_Block::LoadSendBuffer_Flux_F2C --                      *
 *                     Loads send message buffer for fine to coarse block      *
 *                     message passing of conservative solution fluxes.        *
 *******************************************************************************/
inline int Electrostatic2D_Quad_Block::LoadSendBuffer_Flux_F2C(double *buffer,
							       int &buffer_count,
							       const int buffer_size,
							       const int i_min,
							       const int i_max,
							       const int i_inc,
							       const int j_min,
							       const int j_max,
							       const int j_inc) {

  if (j_min == j_max && j_min == JCl) {
    for (int i = i_min;  ((i_inc+2)/4) ? (i < i_max):(i > i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = FluxS[i  ][k] + FluxS[i+1][k];
      }
    }
  } else if (j_min == j_max && j_min == JCu) {
    for (int i = i_min; ((i_inc+2)/4) ? (i < i_max):(i > i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = FluxN[i  ][k] + FluxN[i+1][k];
      }
    }
  } else if (i_min == i_max && i_min == ICl) {
    for (int j = j_min; ((j_inc+2)/4) ? (j < j_max):(j > j_max); j += j_inc) {
      for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = FluxW[j][k] + FluxW[j+1][k];
      }
    }
  } else if (i_min == i_max && i_min == ICu) {
    for (int j = j_min; ((j_inc+2)/4) ? (j < j_max):(j > j_max); j += j_inc) {
      for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = FluxE[j][k] + FluxE[j+1][k];
      }
    }
  }
  return 0;
}

/**********************************************************************
 * Electrostatic2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C --        *
 *                     Unloads receive message buffer for fine to     *
 *                     coarse block message passing of conservative   *
 *                     solution fluxes.                               *
 **********************************************************************/
inline int Electrostatic2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C(double *buffer,
								    int &buffer_count,
								    const int buffer_size,
								    const int i_min,
								    const int i_max,
								    const int i_inc,
								    const int j_min,
								    const int j_max,
								    const int j_inc) {

  if (j_min == j_max && j_min == JCl) {
    for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	FluxS[i][k] = - buffer[buffer_count] - FluxS[i][k];
      }
    }
  } else if (j_min == j_max && j_min == JCu) {
    for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	FluxN[i][k] = - buffer[buffer_count] - FluxN[i][k];
      }
    }
  } else if (i_min == i_max && i_min == ICl) {
    for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
      for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	FluxW[j][k] = - buffer[buffer_count] - FluxW[j][k];
      }
    }
  } else if (i_min == i_max && i_min == ICu) {
    for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
      for (int k = 1; k <= NUM_VAR_ELECTROSTATIC2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	FluxE[j][k] = - buffer[buffer_count] - FluxE[j][k];
      }
    }
  }
  return 0;
}

/**********************************************************************
 * Electrostatic2D_Quad_Block -- Single Block External Subroutines.   *
 **********************************************************************/

extern void Broadcast_Solution_Block(Electrostatic2D_Quad_Block &SolnBlk);

#ifdef _MPI_VERSION
extern void Broadcast_Solution_Block(Electrostatic2D_Quad_Block &SolnBlk,
                                     MPI::Intracomm &Communicator,
                                     const int Source_CPU);
#endif

extern void Copy_Solution_Block(Electrostatic2D_Quad_Block &SolnBlk1,
		                Electrostatic2D_Quad_Block &SolnBlk2);

extern int Prolong_Solution_Block(Electrostatic2D_Quad_Block &SolnBlk_Fine,
				  Electrostatic2D_Quad_Block &SolnBlk_Original,
				  const int Sector);

extern int Restrict_Solution_Block(Electrostatic2D_Quad_Block &SolnBlk_Coarse,
				   Electrostatic2D_Quad_Block &SolnBlk_Original_SW,
				   Electrostatic2D_Quad_Block &SolnBlk_Original_SE,
				   Electrostatic2D_Quad_Block &SolnBlk_Original_NW,
				   Electrostatic2D_Quad_Block &SolnBlk_Original_NE);

extern void ICs(Electrostatic2D_Quad_Block &SolnBlk,
		Electrostatic2D_Input_Parameters &IP,
                Electrostatic2DState *Uo);

extern void BCs(Electrostatic2D_Quad_Block &SolnBlk,
		Electrostatic2D_Input_Parameters &IP);

extern double CFL(Electrostatic2D_Quad_Block &SolnBlk,
                  Electrostatic2D_Input_Parameters &IP);

extern void Set_Global_TimeStep(Electrostatic2D_Quad_Block &SolnBlk,
                                const double &Dt_min);

extern double L1_Norm_Residual(Electrostatic2D_Quad_Block &SolnBlk);

extern double L2_Norm_Residual(Electrostatic2D_Quad_Block &SolnBlk);

extern double Max_Norm_Residual(Electrostatic2D_Quad_Block &SolnBlk);

extern void Linear_Reconstruction_GreenGauss(Electrostatic2D_Quad_Block &SolnBlk,
                                             const int i,
                                             const int j,
					     const int Limiter);

extern void Linear_Reconstruction_GreenGauss(Electrostatic2D_Quad_Block &SolnBlk,
					     const int Limiter);

extern void Linear_Reconstruction_LeastSquares(Electrostatic2D_Quad_Block &SolnBlk,
                                               const int i,
                                               const int j,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares_2(Electrostatic2D_Quad_Block &SolnBlk,
                                                 const int i,
                                                 const int j,
					         const int Limiter);

extern void Linear_Reconstruction_LeastSquares(Electrostatic2D_Quad_Block &SolnBlk,
					       const int Limiter);

extern void Residual_Smoothing(Electrostatic2D_Quad_Block &SolnBlk,
                               const int k_residual,
			       double &epsilon,
                               const int number_of_Gauss_Seidel_iterations);

extern void Calculate_Refinement_Criteria(double *refinement_criteria,
					  Electrostatic2D_Input_Parameters &IP,
                                          int &number_refinement_criteria,
                                          Electrostatic2D_Quad_Block &SolnBlk);

extern void Fix_Refined_Block_Boundaries(Electrostatic2D_Quad_Block &SolnBlk,
                                         const int Fix_North_Boundary,
                                         const int Fix_South_Boundary,
                                         const int Fix_East_Boundary,
                                         const int Fix_West_Boundary);

extern void Unfix_Refined_Block_Boundaries(Electrostatic2D_Quad_Block &SolnBlk);

extern void Apply_Boundary_Flux_Corrections(Electrostatic2D_Quad_Block &SolnBlk,
                                            const int Number_Neighbours_North_Boundary,
                                            const int Number_Neighbours_South_Boundary,
                                            const int Number_Neighbours_East_Boundary,
                                            const int Number_Neighbours_West_Boundary);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Electrostatic2D_Quad_Block &SolnBlk,
                                                                const int i_stage,
								Electrostatic2D_Input_Parameters &IP,
                                                                const int Number_Neighbours_North_Boundary,
                                                                const int Number_Neighbours_South_Boundary,
                                                                const int Number_Neighbours_East_Boundary,
                                                                const int Number_Neighbours_West_Boundary);

extern int dUdt_Residual_Evaluation(Electrostatic2D_Quad_Block &SolnBlk,
				    Electrostatic2D_Input_Parameters &IP);

extern int dUdt_Multistage_Explicit(Electrostatic2D_Quad_Block &SolnBlk,
   	                            const int i_stage,
                                    Electrostatic2D_Input_Parameters &IP);

extern int Update_Solution_Multistage_Explicit(Electrostatic2D_Quad_Block &SolnBlk,
   	                                       const int i_stage,
					       Electrostatic2D_Input_Parameters &IP);

extern int Determine_Electric_Field(Electrostatic2D_Quad_Block &SolnBlk,
				    Electrostatic2D_Input_Parameters &IP);

/**********************************************************************
 * Electrostatic2D_Quad_Block -- Multiple Block External Subroutines. *
 **********************************************************************/

extern Electrostatic2D_Quad_Block* Allocate(Electrostatic2D_Quad_Block *Soln_ptr,
					    Electrostatic2D_Input_Parameters &Input_Parameters);

extern Electrostatic2D_Quad_Block* Deallocate(Electrostatic2D_Quad_Block *Soln_ptr,
					      Electrostatic2D_Input_Parameters &Input_Parameters);

extern Electrostatic2D_Quad_Block* CreateInitialSolutionBlocks(Grid2D_Quad_Block **InitMeshBlk,
							       Electrostatic2D_Quad_Block *Soln_ptr,
							       Electrostatic2D_Input_Parameters &Input_Parameters,
							       QuadTreeBlock_DataStructure &QuadTree,
							       AdaptiveBlockResourceList &GlobalSolnBlockList,
							       AdaptiveBlock2D_List &LocalSolnBlockList);

extern void ICs(Electrostatic2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
                Electrostatic2D_Input_Parameters &Input_Parameters);

extern void BCs(Electrostatic2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
		Electrostatic2D_Input_Parameters &Input_Parameters);

extern double CFL(Electrostatic2D_Quad_Block *Soln_ptr,
                  AdaptiveBlock2D_List &Soln_Block_List,
		  Electrostatic2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(Electrostatic2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                const double &Dt_min);

extern double L1_Norm_Residual(Electrostatic2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double L2_Norm_Residual(Electrostatic2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double Max_Norm_Residual(Electrostatic2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List);

extern void Residual_Smoothing(Electrostatic2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Electrostatic2D_Input_Parameters &Input_Parameters,
   	                       const int I_Stage);

extern void Apply_Boundary_Flux_Corrections(Electrostatic2D_Quad_Block *Soln_ptr,
                                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Electrostatic2D_Quad_Block *Soln_ptr,
                                                                AdaptiveBlock2D_List &Soln_Block_List,
                                                                Electrostatic2D_Input_Parameters &Input_Parameters,
   	                                                        const int I_Stage);

extern int dUdt_Residual_Evaluation(Electrostatic2D_Quad_Block *Soln_ptr,
                                    AdaptiveBlock2D_List &Soln_Block_List,
                                    Electrostatic2D_Input_Parameters &Input_Parameters);

extern int dUdt_Multistage_Explicit(Electrostatic2D_Quad_Block *Soln_ptr,
                                    AdaptiveBlock2D_List &Soln_Block_List,
                                    Electrostatic2D_Input_Parameters &Input_Parameters,
   	                            const int I_Stage);

extern int Update_Solution_Multistage_Explicit(Electrostatic2D_Quad_Block *Soln_ptr,
                                               AdaptiveBlock2D_List &Soln_Block_List,
                                               Electrostatic2D_Input_Parameters &Input_Parameters,
   	                                       const int I_Stage);

extern int Determine_Electric_Field(Electrostatic2D_Quad_Block *Soln_ptr,
				    AdaptiveBlock2D_List &Soln_Block_List,
				    Electrostatic2D_Input_Parameters &Input_Parameters);

/**********************************************************************
 * Electrostatic2D_Quad_Block -- IO Single Block External Subroutines.*
 **********************************************************************/

extern void Write_Solution_Block(Electrostatic2D_Quad_Block &SolnBlk,
	                         ostream &Out_File);

extern void Read_Solution_Block(Electrostatic2D_Quad_Block &SolnBlk,
	                        istream &In_File);

extern void Output_Tecplot(Electrostatic2D_Quad_Block &SolnBlk,
			   Electrostatic2D_Input_Parameters &IP,
		           const int Number_of_Time_Steps,
                           const double &Time,
                           const int Block_Number,
                           const int Output_Title,
	                   ostream &Out_File);

extern void Output_Cells_Tecplot(Electrostatic2D_Quad_Block &SolnBlk,
				 Electrostatic2D_Input_Parameters &IP,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
				 const int CPU_Number,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_Quasi3D_Tecplot(Electrostatic2D_Quad_Block &SolnBlk,
				   Electrostatic2D_Input_Parameters &IP,
				   const int Number_of_Time_Steps,
				   const double &Time,
				   const int Block_Number,
				   const int Output_Title,
				   ostream &Out_File);

/**********************************************************************
 * Electrostatic2D_Quad_Block -- IO Multiple Block External           *
                                 Subroutines.                         *
 **********************************************************************/

extern int Read_Restart_Solution(Electrostatic2D_Quad_Block *Soln_ptr,
                                 AdaptiveBlock2D_List &Soln_Block_List,
                                 Electrostatic2D_Input_Parameters &IP,
		                 int &Number_of_Time_Steps,
                                 double &Time,
                                 CPUTime &CPU_Time);

extern int Write_Restart_Solution(Electrostatic2D_Quad_Block *Soln_ptr,
                                  AdaptiveBlock2D_List &Soln_Block_List,
                                  Electrostatic2D_Input_Parameters &IP,
		                  const int Number_of_Time_Steps,
                                  const double &Time,
                                  const CPUTime &CPU_Time);

extern int Output_Tecplot(Electrostatic2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          Electrostatic2D_Input_Parameters &IP,
		          const int Number_of_Time_Steps,
                          const double &Time);

extern int Output_Cells_Tecplot(Electrostatic2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                Electrostatic2D_Input_Parameters &IP,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Quasi3D_Tecplot(Electrostatic2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  Electrostatic2D_Input_Parameters &IP,
				  const int Number_of_Time_Steps,
				  const double &Time);

extern int Output_Mesh_Tecplot(Electrostatic2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Electrostatic2D_Input_Parameters &IP,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Mesh_Gnuplot(Electrostatic2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Electrostatic2D_Input_Parameters &IP,
		               const int Number_of_Time_Steps,
                               const double &Time);

/**********************************************************************
 * Electrostatic2D_Quad_Block -- Grid Multiple Block External         *
 *                               Subroutines.                         *
 **********************************************************************/

extern Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
					    Electrostatic2D_Input_Parameters &IP);

extern Grid2D_Quad_Block** Broadcast_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
						      Electrostatic2D_Input_Parameters &IP);

extern int Write_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                             Electrostatic2D_Input_Parameters &IP);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
							    Electrostatic2D_Input_Parameters &IP);

extern int Write_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                  Electrostatic2D_Input_Parameters &IP);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
						 Electrostatic2D_Input_Parameters &IP);

extern int Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                          Electrostatic2D_Input_Parameters &IP);

extern int Output_Nodes_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                Electrostatic2D_Input_Parameters &IP);

extern int Output_Cells_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                Electrostatic2D_Input_Parameters &IP);

/**********************************************************************
 * Electrostatic2D_Quad_Block -- Solvers.                             *
 **********************************************************************/

extern int Electrostatic2DQuadSolver(char *Input_File_Name_ptr,
				     int batch_flag);

extern Electrostatic2D_Quad_Block* Electrostatic2DQuadSolver(char *Input_File_Name,
							     int batch_flag,
							     Electrostatic2D_Quad_Block *Local_SolnBlk,
							     Electrostatic2D_Input_Parameters &Input_Parameters,
							     QuadTreeBlock_DataStructure &QuadTree,
							     AdaptiveBlockResourceList &List_of_Global_Solution_Blocks,
							     AdaptiveBlock2D_List &List_of_Local_Solution_Blocks);

#endif // _ELECTROSTATIC2D_QUAD_INCLUDED
