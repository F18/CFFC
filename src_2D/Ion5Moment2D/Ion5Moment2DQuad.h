/* Ion5Moment2DQuad.h:  Header file defining 
                        2D 5-Moment Ion Transport Model 
                        Quadrilateral Mesh Solution Classes. */

#ifndef _ION5MOMENT2D_QUAD_INCLUDED
#define _ION5MOMENT2D_QUAD_INCLUDED

/* Include 2D 5-moment ion transport state, 2D Euler state, 
   2D cell, 2D quadrilateral multiblock grid, quadtree, AMR, and
   2D 5-moment ion transport model input header files. */

#ifndef _ION5MOMENT2D_STATE_INCLUDED
#include "Ion5Moment2DState.h"
#endif // _ION5MOMENT2D_STATE_INCLUDED

#ifndef _EULER2D_STATE_INCLUDED
#include "../Euler2D/Euler2DState.h"
#endif // _EULER2D_STATE_INCLUDED

#ifndef _CELL2D_INCLUDED
#include "../Grid/Cell2D.h"
#endif // _CELL2D_INCLUDED

#ifndef _GRID2D_QUAD_BLOCK_INCLUDED
#include "../Grid/Grid2DQuad.h"
#endif // _GRID2D_QUAD_BLOCK_INCLUDED

#ifndef _QUADTREE_INCLUDED
#include "../AMR/QuadTree.h"
#endif // _QUADTREE_INCLUDED

#ifndef _AMR_INCLUDED
#include "../AMR/AMR.h"
#endif // _AMR_INCLUDED

#ifndef _ION5MOMENT2D_INPUT_INCLUDED
#include "Ion5Moment2DInput.h"
#endif // _ION5MOMENT2D_INPUT_INCLUDED

/* Also include multigrid and linear systems header files. */

#ifndef _FASMULTIGRID2D_INCLUDED
#include "../FASMultigrid2D/FASMultigrid2D.h"
#endif // _FASMULTIGRID2D_INCLUDED

#ifndef _LINEARSYSTEMS_INCLUDED
#include "../Math/LinearSystems.h"
#endif // _LINEARSYSTEMS_INCLUDED

/* Include ICEMCFD input header file. */

#ifndef _ICEMCFD_INCLUDED
#include "../ICEM/ICEMCFD.h"
#endif // _ICEMCFD_INCLUDED

/* Define the structures and classes. */

#define	NUMBER_OF_RESIDUAL_VECTORS_ION5MOMENT2D    3

/*!
 * Class: Ion5Moment2D_Quad_Block
 *
 * @brief Class definition of the 2D Ion5Moment solution blocks.
 *
 * \verbatim
 * Member functions
 *       W      -- Return ion primitive variable solution for block 
 *                 (cell average).
 *       U      -- Return ion conserved variable solution for block
 *                 (cell average).
 *   Wneut      -- Return neutral gas primitive variable solution for
 *                 block (cell average).
 *       E      -- Return electric field for the solution block (cell
 *                 centred).
 *       V      -- Return electric potential for the solution block 
 *                 (cell centred).
 *    Grid      -- Return the solution block quadrilateral grid or mesh.
 *      dt      -- Return local time step for the solution block.
 *    dUdt      -- Return the ion solution block residuals.
 *    dWdx      -- Return the unlimited ion primitive variable solution 
 *                 gradients (x-direction) for the block.
 *    dWdy      -- Return the unlimited ion primitive variable solution
 *                 gradients (y-direction) for the block.
 *     phi      -- Return the ion solution slope limiters.
 *      Uo      -- Return initial ion solution state.
 *    FluxN     -- Return array of north boundary solution fluxes.
 *    FluxS     -- Return array of south boundary solution fluxes.
 *    FluxE     -- Return array of east boundary solution fluxes.
 *    FluxW     -- Return array of west boundary solution fluxes.
 *     NCi      -- Return number of cells in the i-direction
 *                 (zeta-direction).
 *     ICl      -- Return lower index for cells in the i-direction
 *                 (zeta-direction).
 *     ICu      -- Return upper index for cells in the i-direction 
 *                 (zeta-direction).
 *     NCj      -- Return number of cells in the j-direction (eta-
 *                 direction).
 *     JCl      -- Return lower index for cells in the j-direction (eta-
 *                 direction).
 *     JCu      -- Return upper index for cells in the j-direction (eta-
 *                 direction).
 *  Nghost      -- Return number of ghost (halo or overlap) cells.
 * Axisymmetric -- Return axisymmetric flow indicator (=1 for
 *                 axisymmetric flow, =0 for planar flow).
 *     WoN      -- Return array of reference ion states for the 
 *                 application of boundary conditions at the north 
 *                 boundary of the solution block.
 *     WoS      -- Return array of reference ion states for the
 *                 application of boundary conditions at the south 
 *                 boundary of the solution block.
 *     WoE      -- Return array of reference ion states for the 
 *                 application of boundary conditions at the east
 *                 boundary of the solution block.
 *     WoW      -- Return array of reference ion states for the
 *                 application of boundary conditions at the west
 *                 boundary of the solution block.
 *   allocate   -- Allocate memory for structured quadrilateral
 *                 solution block.
 *   deallocate -- Deallocate memory for structured quadrilateral
 *                 solution block.
 *      Wn      -- Return ion primitive variable solution at the
 *                 specified node.
 *      Un      -- Return ion conserved variable solution at the
 *                 specified node.
 *  Wneutn      -- Return neutral gas primitive variable solution at the
 *                 specified node.
 *      En      -- Return electric field at specified node.
 *      Vn      -- Return electric potential at specified node.
 *    WnNW      -- Return ion primitive variable solution at the north-
 *                 west node.
 *    WnNE      -- Return ion primitive variable solution at the north-
 *                 east node.
 *    WnSW      -- Return ion primitive variable solution at the south-
 *                 west node.
 *    WnSE      -- Return ion primitive variable solution at the south-
 *                 east node.
 *    UnNW      -- Return ion conserved variable solution at the north-
 *                 west node.
 *    UnNE      -- Return ion conserved variable solution at the north-
 *                 east node.
 *    UnSW      -- Return ion conserved variable solution at the south-
 *                 west node.
 *    UnSE      -- Return ion conserved variable solution at the south-
 *                 east node.
 *    WneutnNW  -- Return neutral gas primitive variable solution at the 
 *                 north-west node.
 *    WneutnNE  -- Return neutral gas primitive variable solution at the
 *                 north-east node.
 *    WneutnSW  -- Return neutral gas primitive variable solution at the
 *                 south-west node.
 *    WneutnSE  -- Return neutral gas primitive variable solution at the
 *                 south-east node.
 *    EnNW      -- Return electric field at the north-west node.
 *    EnNE      -- Return electric field at the north-east node.
 *    EnSW      -- Return electric field at the south-west node.
 *    EnSE      -- Return electric field at the south-east node.
 *    VnNW      -- Return electric potential at the north-west node.
 *    VnNE      -- Return electric potential at the north-east node.
 *    VnSW      -- Return electric potential at the south-west node.
 *    VnSE      -- Return electric potential at the south-east node.
 * Member functions required for message passing.
 *    NumVar    -- Returns number of solution variables in primitive and
 *                 conserved solution state vectors.
 * LoadSendBuffer -- Loads send buffer.
 * LoadSendBuffer_F2C -- Loads send buffer for fine to coarse block 
 *                       messages.
 * UnloadReceiveBuffer -- Unloads receive buffer.
 * UnloadReceiveBuffer_F2C -- Unloads receive buffer for fine to coarse
 *                            block messages.
 * SubcellReconstruction -- Performs subcell solution reconstruction
 *                          used in adaptive mesh refinement.
 * LoadSendBuffer_Flux_F2C -- Loads send buffer for sending conservative 
 *                            flux corrections from fine to coarse
 *                            solution blocks.
 * UnLoadSendBuffer_Flux_F2C -- Loads send buffer for sending
 *                              conservative flux corrections from fine
 *                              to coarse solution blocks.
 *
 * Member operators
 *      S -- a 2D 5-Moment Ion Transport Model Solution 
 *
 * S = S;
 * cout << S; (output function)
 * cin  >> S; (input function)
 * \endverbatim
 */
class Ion5Moment2D_Quad_Block{
private:
public:
  //@{ @name Solution state arrays:
  Ion5Moment2D_pState           **W; //!< Ion primitive solution state.
  Ion5Moment2D_cState           **U; //!< Ion conserved solution state.
  Euler2D_pState            **Wneut; //!< Neutral gas primitive solution state.
  Vector2D                      **E; //!< Electric field.
  double                        **V; //!< Electric potential.
  //@}

  //@{ @name Grid block information:
  int                           NCi, //!< Total number of i-direction cells.
                                ICl, //!< First i-direction non-ghost cell counter.
                                ICu; //!< Final i-direction non-ghost cell counter.
  int                           NCj, //!< Total number of j-direction cells.
                                JCl, //!< First j-direction non-ghost cell counter.
                                JCu; //!< Final j-direction non-ghost cell counter.
  int                        Nghost; //!< Number of ghost cells.
  Grid2D_Quad_Block            Grid; //!< 2D quadrilateral grid geometry.
  //@}

  //@{ @name Residual and time-stepping arrays:
  double                       **dt; //!< Local time step.
  Ion5Moment2D_cState       ***dUdt; //!< Ion solution residual.
  Ion5Moment2D_cState          **Uo; //!< Initial ion solution state.
  static int      residual_variable; //!< Static integer that indicates which variable is used for residual calculations.
  //@}

  //@{ @name Solution gradient arrays:
  Ion5Moment2D_pState        **dWdx; //!< Unlimited ion solution gradient (x-direction).
  Ion5Moment2D_pState        **dWdy; //!< Unlimited ion solution gradient (y-direction).
  Ion5Moment2D_pState         **phi; //!< Ion solution slope limiter.
  //@}

  //@{ @name Boundary solution flux arrays:
  Ion5Moment2D_cState        *FluxN, //!< North boundary solution flux.
                             *FluxS, //!< South boundary solution flux.
                             *FluxE, //!< East boundary solution flux.
                             *FluxW; //!< West boundary solution flux.
  //@}

  //@{ @name Problem indicator flags:
  int                  Axisymmetric; //!< Axisymmetric flow indicator.
  //int                Freeze_Limiter; //!< Limiter freezing indicator.
  //@}

  //@{ @name Boundary condtion reference states:
  Ion5Moment2D_pState          *WoN, //!< Boundary condition reference states for north boundary.
                               *WoS, //!< Boundary condition reference states for south boundary.
                               *WoE, //!< Boundary condition reference states for east boundary.
                               *WoW; //!< Boundary condition reference states for west boundary.
  //@}
		      
  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  Ion5Moment2D_Quad_Block(void) {
    NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
    W = NULL; U = NULL; Wneut = NULL; E = NULL;  V = NULL;
    dt = NULL; dUdt = NULL; dWdx = NULL; dWdy = NULL; phi = NULL; Uo = NULL;
    FluxN = NULL; FluxS = NULL; FluxE = NULL; FluxW = NULL;
    WoN = NULL; WoS = NULL; WoE = NULL; WoW = NULL;
    Axisymmetric = 0;
  }

  //! Copy constructor.
  Ion5Moment2D_Quad_Block(const Ion5Moment2D_Quad_Block &Soln) {
    NCi = Soln.NCi; ICl = Soln.ICl; ICu = Soln.ICu; 
    NCj = Soln.NCj; JCl = Soln.JCl; JCu = Soln.JCu; Nghost = Soln.Nghost;
    Grid = Soln.Grid; W = Soln.W; U = Soln.U;
    Wneut = Soln.Wneut; E = Soln.E;  V = Soln.V; 
    dt = Soln.dt; dUdt = Soln.dUdt; 
    dWdx = Soln.dWdx; dWdy = Soln.dWdy; phi = Soln.phi; Uo = Soln.Uo;
    FluxN = Soln.FluxN; FluxS = Soln.FluxS; FluxE = Soln.FluxE; FluxW = Soln.FluxW;
    WoN = Soln.WoN; WoS = Soln.WoS; WoE = Soln.WoE; WoW = Soln.WoW;
    Axisymmetric = 0;
  }

  /* Destructor. */
  // ~Ion5Moment2D_Quad_Block(void);
  // Use automatically generated destructor.
  //@}

  /* Assignment operator. */
  // Ion5Moment2D_Quad_Block operator = (const Ion5Moment2D_Quad_Block &Soln);
  // Use automatically generated assignment operator.

  //@{ @name Allocate and deallocate functions.
  //! Allocate memory for structured quadrilateral solution block.
  void allocate(const int Ni, const int Nj, const int Ng);
  //! Deallocate memory for structured quadrilateral solution block.
  void deallocate(void);
  //@}

  //@{ @name Bilinear interplation (Zingg & Yarrow).
  //! Return primitive solution state at specified node.
  Ion5Moment2D_pState Wn(const int ii, const int jj);

  //! Return conserverd solution state at specified node.
  Ion5Moment2D_cState Un(const int ii, const int jj);

  //! Return neutral gas primitive solution state at specified node.
  Euler2D_pState Wneutn(const int ii, const int jj);

  //! Return electric field at specified node.
  Vector2D En(const int ii, const int jj);

  //! Return electric potential at specified node.
  double Vn(const int ii, const int jj);

  Ion5Moment2D_pState WnNW(const int ii, const int jj); //!< Return primitive solution state at cell NW node.
  Ion5Moment2D_pState WnNE(const int ii, const int jj); //!< Return primitive solution state at cell NE node.
  Ion5Moment2D_pState WnSE(const int ii, const int jj); //!< Return primitive solution state at cell SE node.
  Ion5Moment2D_pState WnSW(const int ii, const int jj); //!< Return primitive solution state at cell SW node.

  Ion5Moment2D_cState UnNW(const int ii, const int jj); //!< Return conserved solution state at cell NW node.
  Ion5Moment2D_cState UnNE(const int ii, const int jj); //!< Return conserved solution state at cell NE node.
  Ion5Moment2D_cState UnSE(const int ii, const int jj); //!< Return conserved solution state at cell SE node.
  Ion5Moment2D_cState UnSW(const int ii, const int jj); //!< Return conserved solution state at cell SW node.

  Euler2D_pState WneutnNW(const int ii, const int jj); //!< Return neutral gas primitive solution state at cell NW node.
  Euler2D_pState WneutnNE(const int ii, const int jj); //!< Return neutral gas primitive solution state at cell NE node.
  Euler2D_pState WneutnSE(const int ii, const int jj); //!< Return neutral gas primitive solution state at cell SE node.
  Euler2D_pState WneutnSW(const int ii, const int jj); //!< Return neutral gas primitive solution state at cell SW node.

  Vector2D EnNW(const int ii, const int jj); //!< Return electric field at cell NW node.
  Vector2D EnNE(const int ii, const int jj); //!< Return electric field at cell NE node.
  Vector2D EnSE(const int ii, const int jj); //!< Return electric field at cell SE node.
  Vector2D EnSW(const int ii, const int jj); //!< Return electric field at cell SW node.

  double VnNW(const int ii, const int jj); //!< Return electric potential at cell NW node.
  double VnNE(const int ii, const int jj); //!< Return electric potential at cell NE node.
  double VnSE(const int ii, const int jj); //!< Return electric potential at cell SE node.
  double VnSW(const int ii, const int jj); //!< Return electric potential at cell SW node.
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file,
			       const Ion5Moment2D_Quad_Block &Soln);
  friend istream &operator >> (istream &in_file,
			       Ion5Moment2D_Quad_Block &Soln);
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
  //! Subcell solution reconstruction within given computational cell.
  void SubcellReconstruction(const int i,
			     const int j,
			     const int Limiter);
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

/*******************************************************************************
 * Ion5Moment2D_Quad_Block::allocate -- Allocate memory.                       *
 *******************************************************************************/
inline void Ion5Moment2D_Quad_Block::allocate(const int Ni, const int Nj, const int Ng) {
   int i, j; assert(Ni > 1 && Nj > 1 && Ng > 1); Grid.allocate(Ni, Nj, Ng);
   NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1;
   NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1; Nghost = Ng;
   W = new Ion5Moment2D_pState*[NCi]; U = new Ion5Moment2D_cState*[NCi];
   Wneut = new Euler2D_pState*[NCi]; E = new Vector2D*[NCi]; V = new double*[NCi];
   dt = new double*[NCi]; dUdt = new Ion5Moment2D_cState**[NCi];
   dWdx = new Ion5Moment2D_pState*[NCi]; dWdy = new Ion5Moment2D_pState*[NCi];
   phi = new Ion5Moment2D_pState*[NCi]; Uo = new Ion5Moment2D_cState*[NCi];
   for ( i = 0; i <= NCi-1 ; ++i ) {
      W[i] = new Ion5Moment2D_pState[NCj]; U[i] = new Ion5Moment2D_cState[NCj];
      Wneut[i] = new Euler2D_pState[NCj]; E[i] = new Vector2D[NCj]; V[i] = new double[NCj];
      dt[i] = new double[NCj]; dUdt[i] = new Ion5Moment2D_cState*[NCj];
      for ( j = 0; j <= NCj-1 ; ++j ) 
        { dUdt[i][j] = new Ion5Moment2D_cState[NUMBER_OF_RESIDUAL_VECTORS_ION5MOMENT2D]; }
      dWdx[i] = new Ion5Moment2D_pState[NCj]; dWdy[i] = new Ion5Moment2D_pState[NCj];
      phi[i] = new Ion5Moment2D_pState[NCj]; Uo[i] = new Ion5Moment2D_cState[NCj];
   } /* endfor */
   FluxN = new Ion5Moment2D_cState[NCi]; FluxS = new Ion5Moment2D_cState[NCi];
   FluxE = new Ion5Moment2D_cState[NCj]; FluxW = new Ion5Moment2D_cState[NCj];
   WoN = new Ion5Moment2D_pState[NCi]; WoS = new Ion5Moment2D_pState[NCi];
   WoE = new Ion5Moment2D_pState[NCj]; WoW = new Ion5Moment2D_pState[NCj];
}

/*******************************************************************************
 * Ion5Moment2D_Quad_Block::deallocate -- Deallocate memory.                   *
 *******************************************************************************/
inline void Ion5Moment2D_Quad_Block::deallocate(void) {
   int i, j; Grid.deallocate();
   for ( i = 0; i <= NCi-1 ; ++i ) {
      delete []W[i]; W[i] = NULL; delete []U[i]; U[i] = NULL;
      delete []Wneut[i]; Wneut[i] = NULL; delete []E[i]; E[i] = NULL; delete []V[i]; V[i] = NULL;
      delete []dt[i]; dt[i] = NULL; 
      for ( j = 0; j <= NCj-1 ; ++j ) { delete []dUdt[i][j]; dUdt[i][j] = NULL; }
      delete []dUdt[i]; dUdt[i] = NULL;
      delete []dWdx[i]; dWdx[i] = NULL; delete []dWdy[i]; dWdy[i] = NULL;
      delete []phi[i]; phi[i] = NULL; delete []Uo[i]; Uo[i] = NULL;
   } /* endfor */
   delete []W; W = NULL; delete []U; U = NULL;
   delete []Wneut; Wneut = NULL; delete []E; E = NULL; delete []V; V = NULL;
   delete []dt; dt = NULL; delete []dUdt; dUdt = NULL;
   delete []dWdx; dWdx = NULL; delete []dWdy; dWdy = NULL; 
   delete []phi; phi = NULL; delete []Uo; Uo = NULL;
   delete []FluxN; FluxN = NULL; delete []FluxS; FluxS = NULL;
   delete []FluxE; FluxE = NULL; delete []FluxW; FluxW = NULL;
   delete []WoN; WoN = NULL; delete []WoS; WoS = NULL;
   delete []WoE; WoE = NULL; delete []WoW; WoW = NULL;
   NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
}

/*******************************************************************************
 * Ion5Moment2D_Quad_Block::Wn -- Node ion primitive solution state.           *
 *******************************************************************************/
inline Ion5Moment2D_pState Ion5Moment2D_Quad_Block::Wn(const int ii, const int jj) {
    double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
           eta1, zeta1, eta2, zeta2, eta, zeta;
    Ion5Moment2D_pState A, B, C, D;
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
    A=W[ii-1][jj-1]; B=W[ii-1][jj]-W[ii-1][jj-1]; C=W[ii][jj-1]-W[ii-1][jj-1];
    D=W[ii][jj]+W[ii-1][jj-1]-W[ii-1][jj]-W[ii][jj-1];
    return (A+B*zeta+C*eta+D*zeta*eta);
}

/*********************************************************************************
 * Ion5Moment2D_Quad_Block::Un -- Node ion conserved solution state.             *
 *********************************************************************************/
inline Ion5Moment2D_cState Ion5Moment2D_Quad_Block::Un(const int ii, const int jj) {
    double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
           eta1, zeta1, eta2, zeta2, eta, zeta;
    Ion5Moment2D_cState A, B, C, D;
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
    A=U[ii-1][jj-1]; B=U[ii-1][jj]-U[ii-1][jj-1]; C=U[ii][jj-1]-U[ii-1][jj-1];
    D=U[ii][jj]+U[ii-1][jj-1]-U[ii-1][jj]-U[ii][jj-1];
    return (A+B*zeta+C*eta+D*zeta*eta);
}

/*********************************************************************************
 * Ion5Moment2D_Quad_Block::Wneutn -- Node neutral gas primitive solution state. *
 *********************************************************************************/
inline Euler2D_pState Ion5Moment2D_Quad_Block::Wneutn(const int ii, const int jj) {
    double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
           eta1, zeta1, eta2, zeta2, eta, zeta;
    Euler2D_pState A, B, C, D;
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
    A=Wneut[ii-1][jj-1]; B=Wneut[ii-1][jj]-Wneut[ii-1][jj-1]; C=Wneut[ii][jj-1]-Wneut[ii-1][jj-1];
    D=Wneut[ii][jj]+Wneut[ii-1][jj-1]-Wneut[ii-1][jj]-Wneut[ii][jj-1];
    return (A+B*zeta+C*eta+D*zeta*eta);
}

/*********************************************************************************
 * Ion5Moment2D_Quad_Block::En -- Node electric field.                           *
 *********************************************************************************/
inline Vector2D Ion5Moment2D_Quad_Block::En(const int ii, const int jj) {
    double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
           eta1, zeta1, eta2, zeta2, eta, zeta;
    Vector2D A, B, C, D;
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
    A=E[ii-1][jj-1]; B=E[ii-1][jj]-E[ii-1][jj-1]; C=E[ii][jj-1]-E[ii-1][jj-1];
    D=E[ii][jj]+E[ii-1][jj-1]-E[ii-1][jj]-E[ii][jj-1];
    return (A+B*zeta+C*eta+D*zeta*eta);
}

/*********************************************************************************
 * Ion5Moment2D_Quad_Block::Vn -- Node electric potential.                       *
 *********************************************************************************/
inline double Ion5Moment2D_Quad_Block::Vn(const int ii, const int jj) {
    double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
           eta1, zeta1, eta2, zeta2, eta, zeta;
    double A, B, C, D;
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
    A=V[ii-1][jj-1]; B=V[ii-1][jj]-V[ii-1][jj-1]; C=V[ii][jj-1]-V[ii-1][jj-1];
    D=V[ii][jj]+V[ii-1][jj-1]-V[ii-1][jj]-V[ii][jj-1];
    return (A+B*zeta+C*eta+D*zeta*eta);
}

/*********************************************************************************
 * Ion5Moment2D_Quad_Block::Wn?? -- Get cell node ion primitive solution states. *
 *********************************************************************************/
inline Ion5Moment2D_pState Ion5Moment2D_Quad_Block::WnNW(const int ii, const int jj) {
  return (Wn(ii, jj+1));
}

inline Ion5Moment2D_pState Ion5Moment2D_Quad_Block::WnNE(const int ii, const int jj) {
  return (Wn(ii+1, jj+1));
}

inline Ion5Moment2D_pState Ion5Moment2D_Quad_Block::WnSE(const int ii, const int jj) {
  return (Wn(ii+1, jj));
}

inline Ion5Moment2D_pState Ion5Moment2D_Quad_Block::WnSW(const int ii, const int jj) {
  return (Wn(ii, jj));
}

/*********************************************************************************
 * Ion5Moment2D_Quad_Block::Un?? -- Get cell node ion conserved solution states. *
 *********************************************************************************/
inline Ion5Moment2D_cState Ion5Moment2D_Quad_Block::UnNW(const int ii, const int jj) {
  return (Un(ii, jj+1));
}

inline Ion5Moment2D_cState Ion5Moment2D_Quad_Block::UnNE(const int ii, const int jj) {
  return (Un(ii+1, jj+1));
}

inline Ion5Moment2D_cState Ion5Moment2D_Quad_Block::UnSE(const int ii, const int jj) {
  return (Un(ii+1, jj));
}

inline Ion5Moment2D_cState Ion5Moment2D_Quad_Block::UnSW(const int ii, const int jj) {
  return (Un(ii, jj));
}

/***********************************************************************************
 * Ion5Moment2D_Quad_Block::Wneutn?? -- Get cell node neutral gas solution states. *
 ***********************************************************************************/
inline Euler2D_pState Ion5Moment2D_Quad_Block::WneutnNW(const int ii, const int jj) {
  return (Wneutn(ii, jj+1));
}

inline Euler2D_pState Ion5Moment2D_Quad_Block::WneutnNE(const int ii, const int jj) {
  return (Wneutn(ii+1, jj+1));
}

inline Euler2D_pState Ion5Moment2D_Quad_Block::WneutnSE(const int ii, const int jj) {
  return (Wneutn(ii+1, jj));
}

inline Euler2D_pState Ion5Moment2D_Quad_Block::WneutnSW(const int ii, const int jj) {
  return (Wneutn(ii, jj));
}

/*********************************************************************************
 * Ion5Moment2D_Quad_Block::En?? -- Get cell node electric field.                *
 *********************************************************************************/
inline Vector2D Ion5Moment2D_Quad_Block::EnNW(const int ii, const int jj) {
  return (En(ii, jj+1));
}

inline Vector2D Ion5Moment2D_Quad_Block::EnNE(const int ii, const int jj) {
  return (En(ii+1, jj+1));
}

inline Vector2D Ion5Moment2D_Quad_Block::EnSE(const int ii, const int jj) {
  return (En(ii+1, jj));
}

inline Vector2D Ion5Moment2D_Quad_Block::EnSW(const int ii, const int jj) {
  return (En(ii, jj));
}

/*********************************************************************************
 * Ion5Moment2D_Quad_Block::Vn?? -- Get cell node electric potential.            *
 *********************************************************************************/
inline double Ion5Moment2D_Quad_Block::VnNW(const int ii, const int jj) {
  return (Vn(ii, jj+1));
}

inline double Ion5Moment2D_Quad_Block::VnNE(const int ii, const int jj) {
  return (Vn(ii+1, jj+1));
}

inline double Ion5Moment2D_Quad_Block::VnSE(const int ii, const int jj) {
  return (Vn(ii+1, jj));
}

inline double Ion5Moment2D_Quad_Block::VnSW(const int ii, const int jj) {
  return (Vn(ii, jj));
}

/*********************************************************************************
 * Ion5Moment2D_Quad_Block -- Input-output operators.                            *
 *********************************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Ion5Moment2D_Quad_Block &SolnBlk) {
  int i, j; 
  out_file << SolnBlk.Grid;
  out_file << SolnBlk.NCi << " " << SolnBlk.ICl << " " 
           << SolnBlk.ICu << " " << SolnBlk.Nghost << "\n";
  out_file << SolnBlk.NCj << " " << SolnBlk.JCl << " " << SolnBlk.JCu << "\n";
  out_file << SolnBlk.Axisymmetric << "\n";
  if (SolnBlk.NCi == 0 || SolnBlk.NCj == 0) return(out_file);
  for ( j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
     for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
         out_file << SolnBlk.U[i][j] << " " 
                  << SolnBlk.Wneut[i][j] << " " 
                  << SolnBlk.E[i][j] << " ";
         out_file.setf(ios::scientific);
         out_file << SolnBlk.V[i][j] << "\n";
         out_file.unsetf(ios::scientific);
     } /* endfor */
  } /* endfor */
  for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
     out_file << SolnBlk.WoW[j] << "\n";
     out_file << SolnBlk.WoE[j] << "\n";
  } /* endfor */
  for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
     out_file << SolnBlk.WoS[i] << "\n";
     out_file << SolnBlk.WoN[i] << "\n";
  } /* endfor */
  return (out_file);
}

inline istream &operator >> (istream &in_file,
			     Ion5Moment2D_Quad_Block &SolnBlk) {
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
         in_file >> SolnBlk.U[i][j] 
                 >> SolnBlk.Wneut[i][j] 
                 >> SolnBlk.E[i][j];
         in_file.setf(ios::skipws); 
         in_file >> SolnBlk.V[i][j];
         in_file.unsetf(ios::skipws); 
         SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);
         for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_ION5MOMENT2D-1 ; ++k ) {
	     SolnBlk.dUdt[i][j][k] = Ion5Moment2D_U_VACUUM;
         } /* endfor */
	 SolnBlk.dWdx[i][j] = Ion5Moment2D_W_VACUUM;
	 SolnBlk.dWdy[i][j] = Ion5Moment2D_W_VACUUM;
	 SolnBlk.phi[i][j] = Ion5Moment2D_W_VACUUM;
	 SolnBlk.Uo[i][j] = Ion5Moment2D_U_VACUUM;
	 SolnBlk.dt[i][j] = ZERO;
     } /* endfor */
  } /* endfor */
  for (j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
     in_file >> SolnBlk.WoW[j];
     in_file >> SolnBlk.WoE[j];
  } /* endfor */
  for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
     in_file >> SolnBlk.WoS[i];
     in_file >> SolnBlk.WoN[i];
  } /* endfor */
  return (in_file);
}

/***********************************************************************************
 *                                                                                 *
 * MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING.                                  *
 *                                                                                 *
 ***********************************************************************************/

/***********************************************************************************
 * Ion5Moment2D_Quad_Block::NumVar -- Returns number of state variables.           *
 ***********************************************************************************/
inline int Ion5Moment2D_Quad_Block::NumVar(void) {
  return (int(NUM_VAR_ION5MOMENT2D + NUM_VAR_EULER2D + NUM_COMP_VECTOR2D + 1));
}

/***********************************************************************************
 * Ion5Moment2D_Quad_Block::LoadSendBuffer -- Loads send message buffer.           *
 ***********************************************************************************/
inline int Ion5Moment2D_Quad_Block::LoadSendBuffer(double *buffer,
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
        for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = U[i][j][k];
        } /* endfor */
        for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = Wneut[i][j][k];
        } /* endfor */
        for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = E[i][j][k];
        } /* endfor */
        buffer_count = buffer_count + 1;
        if (buffer_count >= buffer_size) return(1);
        buffer[buffer_count] = V[i][j];
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * Ion5Moment2D_Quad_Block::LoadSendBuffer_F2C -- Loads send message buffer for*
 *                                                fine to coarse block message *
 *                                                passing.                     *
 *******************************************************************************/
inline int Ion5Moment2D_Quad_Block::LoadSendBuffer_F2C(double *buffer,
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
        for ( k = 1 ; k <=  NUM_VAR_ION5MOMENT2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (Grid.Cell[i  ][j  ].A*U[i  ][j  ][k]+
                                   Grid.Cell[i+1][j  ].A*U[i+1][j  ][k]+
                                   Grid.Cell[i  ][j+1].A*U[i  ][j+1][k]+
                                   Grid.Cell[i+1][j+1].A*U[i+1][j+1][k]);
                                  (Grid.Cell[i  ][j  ].A+
                                   Grid.Cell[i+1][j  ].A+
                                   Grid.Cell[i  ][j+1].A+
                                   Grid.Cell[i+1][j+1].A);
/*            buffer[buffer_count] = (Grid.Cell[i  ][j  ].A*U[i  ][j  ][k]+ */
/*                                    Grid.Cell[i+1][j  ].A*U[i+1][j  ][k]+ */
/*                                    Grid.Cell[i  ][j+1].A*U[i  ][j+1][k]+ */
/*                                    Grid.Cell[i+1][j+1].A*U[i+1][j+1][k]); */
        } /* endfor */
        for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (Grid.Cell[i  ][j  ].A*Wneut[i  ][j  ][k]+
                                   Grid.Cell[i+1][j  ].A*Wneut[i+1][j  ][k]+
                                   Grid.Cell[i  ][j+1].A*Wneut[i  ][j+1][k]+
                                   Grid.Cell[i+1][j+1].A*Wneut[i+1][j+1][k]);
                                  (Grid.Cell[i  ][j  ].A+
                                   Grid.Cell[i+1][j  ].A+
                                   Grid.Cell[i  ][j+1].A+
                                   Grid.Cell[i+1][j+1].A);
/*            buffer[buffer_count] = (Grid.Cell[i  ][j  ].A*Wneut[i  ][j  ][k]+ */
/*                                    Grid.Cell[i+1][j  ].A*Wneut[i+1][j  ][k]+ */
/*                                    Grid.Cell[i  ][j+1].A*Wneut[i  ][j+1][k]+ */
/*                                    Grid.Cell[i+1][j+1].A*Wneut[i+1][j+1][k]); */
        } /* endfor */
        for ( k = 1 ; k <= NUM_COMP_VECTOR2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (Grid.Cell[i  ][j  ].A*E[i  ][j  ][k]+
                                   Grid.Cell[i+1][j  ].A*E[i+1][j  ][k]+
                                   Grid.Cell[i  ][j+1].A*E[i  ][j+1][k]+
                                   Grid.Cell[i+1][j+1].A*E[i+1][j+1][k]);
                                  (Grid.Cell[i  ][j  ].A+
                                   Grid.Cell[i+1][j  ].A+
                                   Grid.Cell[i  ][j+1].A+
                                   Grid.Cell[i+1][j+1].A);
/*            buffer[buffer_count] = (Grid.Cell[i  ][j  ].A*E[i  ][j  ][k]+ */
/*                                    Grid.Cell[i+1][j  ].A*E[i+1][j  ][k]+ */
/*                                    Grid.Cell[i  ][j+1].A*E[i  ][j+1][k]+ */
/*                                    Grid.Cell[i+1][j+1].A*E[i+1][j+1][k]); */
        } /* endfor */
        buffer_count = buffer_count + 1;
        if (buffer_count >= buffer_size) return(1);
        buffer[buffer_count] = (Grid.Cell[i  ][j  ].A*V[i  ][j  ]+
                                Grid.Cell[i+1][j  ].A*V[i+1][j  ]+
                                Grid.Cell[i  ][j+1].A*V[i  ][j+1]+
                                Grid.Cell[i+1][j+1].A*V[i+1][j+1]);
                               (Grid.Cell[i  ][j  ].A+
                                Grid.Cell[i+1][j  ].A+
                                Grid.Cell[i  ][j+1].A+
                                Grid.Cell[i+1][j+1].A);
/*         buffer[buffer_count] = (Grid.Cell[i  ][j  ].A*V[i  ][j  ]+ */
/*                                 Grid.Cell[i+1][j  ].A*V[i+1][j  ]+ */
/*                                 Grid.Cell[i  ][j+1].A*V[i  ][j+1]+ */
/*                                 Grid.Cell[i+1][j+1].A*V[i+1][j+1]); */
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * Ion5Moment2D_Quad_Block::LoadSendBuffer_C2F -- Loads send message buffer for*
 *                                                coarse to fine block message *
 *                                                passing.                     *
 *******************************************************************************/
inline int Ion5Moment2D_Quad_Block::LoadSendBuffer_C2F(double *buffer,
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
  Ion5Moment2D_pState Wfine;
  Ion5Moment2D_cState Ufine;

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
	    Wfine = W[i][j] + (phi[i][j]^dWdx[i][j])*dX.x +
                              (phi[i][j]^dWdy[i][j])*dX.y;
	    Ufine = Wfine.U();
	    for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
	      buffer_count = buffer_count + 1;
	      if (buffer_count >= buffer_size) return(1);
	      buffer[buffer_count] = Ufine[k];
	    } /* endfor */
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
	    Wfine = W[i][j] + (phi[i][j]^dWdx[i][j])*dX.x +
                              (phi[i][j]^dWdy[i][j])*dX.y;
	    Ufine = Wfine.U();
	    for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
	      buffer_count = buffer_count + 1;
	      if (buffer_count >= buffer_size) return(1);
	      buffer[buffer_count] = Ufine[k];
	    } /* endfor */
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
	    Wfine = W[i][j] + (phi[i][j]^dWdx[i][j])*dX.x +
                              (phi[i][j]^dWdy[i][j])*dX.y;
	    Ufine = Wfine.U();
	    for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) { 
	      buffer_count = buffer_count + 1;
	      if (buffer_count >= buffer_size) return(1);
	      buffer[buffer_count] = Ufine[k];
	    } /* endfor */
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
	    Wfine = W[i][j] + (phi[i][j]^dWdx[i][j])*dX.x +
                              (phi[i][j]^dWdy[i][j])*dX.y;
	    Ufine = Wfine.U();
	    for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
	      buffer_count = buffer_count + 1;
	      if (buffer_count >= buffer_size) return(1);
	      buffer[buffer_count] = Ufine[k];
	    } /* endfor */
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
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (Grid.Node[i][j_min].X+
                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
                    Grid.Cell[i][j_min].Xc)/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate NE sub (fine) cell values.
              dX = (Grid.Cell[i][j_min].Xc+
                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
                    Grid.Node[i+1][j_min+1].X)/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
                    Grid.Cell[i][j_min].Xc+
                    Grid.Node[i][j_min+1].X+
                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) { 
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
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
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) { 
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (Grid.Cell[i][j_min].Xc+
                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
                    Grid.Node[i+1][j_min+1].X)/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate SW sub (fine) cell values.
              dX = (Grid.Node[i][j_min].X+
                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
                    Grid.Cell[i][j_min].Xc)/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
                    Grid.Node[i+1][j_min].X+
                    Grid.Cell[i][j_min].Xc+
                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
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
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
                    Grid.Cell[i][j_min].Xc+
                    Grid.Node[i][j_min+1].X+
                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) { 
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
                    Grid.Node[i+1][j_min].X+
                    Grid.Cell[i][j_min].Xc+
                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (Grid.Node[i][j_min].X+
                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
                    Grid.Cell[i][j_min].Xc)/FOUR -
                   Grid.Cell[i][j_min].Xc;
              Wfine = W[i][j_min] +
                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } /* endif */
  } else { // East or west boundary.
     // Four different orderings to consider depending on the value of i_inc & j_inc.
     if (j_inc > 0) {
        if (i_inc > 0) {
           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
              // Perform limited linear least squares reconstruction in cell (i_min, j).
              SubcellReconstruction(i_min, j, LIMITER_VENKATAKRISHNAN);
              // Evaluate SW sub (fine) cell values.
              dX = (Grid.Node[i_min][j].X+
                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
                    Grid.Cell[i_min][j].Xc)/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
                    Grid.Node[i_min+1][j].X+
                    Grid.Cell[i_min][j].Xc+
                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
                    Grid.Cell[i_min][j].Xc+
                    Grid.Node[i_min][j+1].X+
                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) { 
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (Grid.Cell[i_min][j].Xc+
                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
                    Grid.Node[i_min+1][j+1].X)/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
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
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (Grid.Node[i_min][j].X+
                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
                    Grid.Cell[i_min][j].Xc)/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
             dX = (Grid.Cell[i_min][j].Xc+
                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
                    Grid.Node[i_min+1][j+1].X)/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
                    Grid.Cell[i_min][j].Xc+
                    Grid.Node[i_min][j+1].X+
                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) { 
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
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
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) { 
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (Grid.Cell[i_min][j].Xc+
                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
                    Grid.Node[i_min+1][j+1].X)/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (Grid.Node[i_min][j].X+
                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
                    Grid.Cell[i_min][j].Xc)/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
                    Grid.Node[i_min+1][j].X+
                    Grid.Cell[i_min][j].Xc+
                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
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
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
                    Grid.Cell[i_min][j].Xc+
                    Grid.Node[i_min][j+1].X+
                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) { 
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
                    Grid.Node[i_min+1][j].X+
                    Grid.Cell[i_min][j].Xc+
                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (Grid.Node[i_min][j].X+
                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
                    Grid.Cell[i_min][j].Xc)/FOUR -
                   Grid.Cell[i_min][j].Xc;
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } /* endif */
  } /* endif */
  return(0);
}

/***********************************************************************************
 * Ion5Moment2D_Quad_Block::UnloadReceiveBuffer -- Unloads receive message buffer. *
 ***********************************************************************************/
inline int Ion5Moment2D_Quad_Block::UnloadReceiveBuffer(double *buffer,
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
  	buffer_count = buffer_count + NUM_VAR_ION5MOMENT2D;
        if (buffer_count >= buffer_size) return(1);
        U[i][j] = Ion5Moment2D_cState(buffer[buffer_count-3],
                                      buffer[buffer_count-2],
                                      buffer[buffer_count-1],
                                      buffer[buffer_count]);
        W[i][j] = U[i][j].W();
  	buffer_count = buffer_count + NUM_VAR_EULER2D;
        if (buffer_count >= buffer_size) return(1);
        Wneut[i][j] = Euler2D_pState(buffer[buffer_count-3],
                                     buffer[buffer_count-2],
                                     buffer[buffer_count-1],
                                     buffer[buffer_count]);
  	buffer_count = buffer_count + NUM_COMP_VECTOR2D;
        if (buffer_count >= buffer_size) return(1);
        E[i][j] = Vector2D(buffer[buffer_count-1],
                           buffer[buffer_count]);
  	buffer_count = buffer_count + 1;
        if (buffer_count >= buffer_size) return(1);
        V[i][j] = buffer[buffer_count];
     } /* endfor */
  } /* endfor */
  return(0);
}

/********************************************************************************
 * Ion5Moment2D_Quad_Block::UnloadReceiveBuffer_F2C -- Unloads receive message  *
 *                                                     buffer for fine to coarse*
 *                                                     block message passing.   *
 ********************************************************************************/
inline int Ion5Moment2D_Quad_Block::UnloadReceiveBuffer_F2C(double *buffer,
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
  	buffer_count = buffer_count + NUM_VAR_ION5MOMENT2D;
        if (buffer_count >= buffer_size) return(1);
        U[i][j] = Ion5Moment2D_cState(buffer[buffer_count-3],
                                      buffer[buffer_count-2],
                                      buffer[buffer_count-1],
                                      buffer[buffer_count]);
/*         U[i][j] = Ion5Moment2D_cState(buffer[buffer_count-3], */
/*                                       buffer[buffer_count-2], */
/*                                       buffer[buffer_count-1], */
/*                                       buffer[buffer_count])/Grid.Cell[i][j].A; */
        W[i][j] = U[i][j].W();
  	buffer_count = buffer_count + NUM_VAR_EULER2D;
        if (buffer_count >= buffer_size) return(1);
        Wneut[i][j] = Euler2D_pState(buffer[buffer_count-3],
                                     buffer[buffer_count-2],
                                     buffer[buffer_count-1],
                                     buffer[buffer_count]);
/*         Wneut[i][j] = Euler2D_pState(buffer[buffer_count-3], */
/*                                      buffer[buffer_count-2], */
/*                                      buffer[buffer_count-1], */
/*                                      buffer[buffer_count])/Grid.Cell[i][j].A; */
  	buffer_count = buffer_count + NUM_COMP_VECTOR2D;
        if (buffer_count >= buffer_size) return(1);
        E[i][j] = Vector2D(buffer[buffer_count-1],
                           buffer[buffer_count]);
/*         E[i][j] = Vector2D(buffer[buffer_count-1], */
/*                            buffer[buffer_count])/Grid.Cell[i][j].A; */
  	buffer_count = buffer_count + 1;
        if (buffer_count >= buffer_size) return(1);
        V[i][j] = buffer[buffer_count];
/*         V[i][j] = buffer[buffer_count]/Grid.Cell[i][j].A; */
     } /* endfor */
  } /* endfor */
  return(0);
}

/********************************************************************************
 * Ion5Moment2D_Quad_Block::UnloadReceiveBuffer_C2F -- Unloads receive message  *
 *                                                     buffer for coarse to fine*
 *                                                     block message passing.   *
 ********************************************************************************/
inline int Ion5Moment2D_Quad_Block::UnloadReceiveBuffer_C2F(double *buffer,
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
  	buffer_count = buffer_count + NUM_VAR_ION5MOMENT2D;
        if (buffer_count >= buffer_size) return(1);
        U[i][j] = Ion5Moment2D_cState(buffer[buffer_count-3],
                                      buffer[buffer_count-2],
                                      buffer[buffer_count-1],
                                      buffer[buffer_count]);
        W[i][j] = U[i][j].W();
  	buffer_count = buffer_count + NUM_VAR_EULER2D;
        if (buffer_count >= buffer_size) return(1);
        Wneut[i][j] = Euler2D_pState(buffer[buffer_count-3],
                                     buffer[buffer_count-2],
                                     buffer[buffer_count-1],
                                     buffer[buffer_count]);
  	buffer_count = buffer_count + NUM_COMP_VECTOR2D;
        if (buffer_count >= buffer_size) return(1);
        E[i][j] = Vector2D(buffer[buffer_count-1],
                           buffer[buffer_count]);
  	buffer_count = buffer_count + 1;
        if (buffer_count >= buffer_size) return(1);
        V[i][j] = buffer[buffer_count];
     } /* endfor */
  } /* endfor */
  return(0);
}

/**************************************************************************
 * Ion5Moment2D_Quad_Block::SubcellReconstruction --                      *
 *               Performs the subcell reconstruction of solution state    *
 *               within a given cell (i,j) of the computational mesh for  *
 *               the specified quadrilateral solution block.              *
 **************************************************************************/
inline void Ion5Moment2D_Quad_Block::SubcellReconstruction(const int i, 
                                                           const int j,
                                                           const int Limiter) {

  int n, n2, n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4], phi_n;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  double Tratio, TratioMin, TratioMax;
  Vector2D dX;
  Ion5Moment2D_pState DU, DUDx_ave, DUDy_ave;

  /* Carry out the limited solution reconstruction in
     each cell of the computational mesh. */

  if (i == ICl-Nghost || i == ICu+Nghost ||
      j == JCl-Nghost || j == JCu+Nghost) {
    n_pts = 0;
  } else if ((i == ICl-Nghost+1) && 
             (Grid.BCtypeW[j] != BC_NONE)) {
    if (j == JCl-Nghost+1 || j == JCu+Nghost-1) {
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
  
  if (n_pts > 0) {
      DUDx_ave = Ion5Moment2D_W_VACUUM;
      DUDy_ave = Ion5Moment2D_W_VACUUM;
      DxDx_ave = ZERO;
      DxDy_ave = ZERO;
      DyDy_ave = ZERO;
  
      for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
          dX = Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - 
               Grid.Cell[i][j].Xc;
          DU = W[ i_index[n2] ][ j_index[n2] ] - 
               W[i][j];
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
  
      for ( n = 1 ; n <= NUM_VAR_ION5MOMENT2D ; ++n ) {
         u0Min = W[i][j][n];
         u0Max = u0Min;
         for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
            u0Min = min(u0Min, W[ i_index[n2] ][ j_index[n2] ][n]);
            u0Max = max(u0Max, W[ i_index[n2] ][ j_index[n2] ][n]);
         } /* endfor */
  
         dX = Grid.xfaceE(i, j)-Grid.Cell[i][j].Xc;
         uQuad[0] = W[i][j][n] + 
                    dWdx[i][j][n]*dX.x +
                    dWdy[i][j][n]*dX.y ;
         dX = Grid.xfaceW(i, j)-Grid.Cell[i][j].Xc;
         uQuad[1] = W[i][j][n] + 
                    dWdx[i][j][n]*dX.x +
                    dWdy[i][j][n]*dX.y ;
         dX = Grid.xfaceN(i, j)-Grid.Cell[i][j].Xc;
         uQuad[2] = W[i][j][n] + 
                    dWdx[i][j][n]*dX.x +
                    dWdy[i][j][n]*dX.y ;
         dX = Grid.xfaceS(i, j)-Grid.Cell[i][j].Xc;
         uQuad[3] = W[i][j][n] + 
                    dWdx[i][j][n]*dX.x +
                    dWdy[i][j][n]*dX.y ;
  
         switch(Limiter) {
           case LIMITER_ONE :
             phi_n = ONE;
             break;
           case LIMITER_ZERO :
             phi_n = ZERO;
             break;
           case LIMITER_BARTH_JESPERSEN :
             phi_n = Limiter_BarthJespersen(uQuad, W[i][j][n], 
                                            u0Min, u0Max, 4);
             break;
           case LIMITER_VENKATAKRISHNAN :
             phi_n = Limiter_Venkatakrishnan(uQuad, W[i][j][n], 
                                             u0Min, u0Max, 4);
             break;
           case LIMITER_VANLEER :
             phi_n = Limiter_VanLeer(uQuad, W[i][j][n], 
                                     u0Min, u0Max, 4);
             break;
           case LIMITER_VANALBADA :
             phi_n = Limiter_VanAlbada(uQuad, W[i][j][n], 
                                       u0Min, u0Max, 4);
             break;
           default:
             phi_n = Limiter_BarthJespersen(uQuad, W[i][j][n], 
                                            u0Min, u0Max, 4);
             break;
         } /* endswitch */
  
         phi[i][j][n] = phi_n;
      } /* endfor */

      TratioMax = ONE;
      TratioMin = ONE;
      for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
          Tratio = W[ i_index[n2] ][ j_index[n2] ].T()/W[i][j].T();
          TratioMin = min(TratioMin, Tratio);
          TratioMax = max(TratioMax, Tratio);
      } /* endfor */
      if (TratioMin < 0.005 ||
          TratioMax > 200.00) {
	 //cout << i << " " << j << " " << TratioMin << " " << TratioMax << "\n";
	 if (TratioMax > 200.00) phi_n = max(TOLER, exp(-(TratioMax-200.00)/25.00));
         if (TratioMin < 0.005) phi_n = max(TOLER, exp(-(ONE/TratioMin-200.00)/25.00));
         phi[i][j]  = phi_n*phi[i][j];
      } /* endif */
  } else {
      dWdx[i][j] = Ion5Moment2D_W_VACUUM;
      dWdy[i][j] = Ion5Moment2D_W_VACUUM; 
      phi[i][j]  = Ion5Moment2D_W_VACUUM;
  } /* endif */

}

/************************************************************************************
 * Ion5Moment2D_Quad_Block::LoadSendBuffer_Flux_F2C -- Loads send message buffer for*
 *                                                     fine to coarse block message *
 *                                                     passing of conservative      *
 *                                                     solution fluxes.             *
 ************************************************************************************/
inline int Ion5Moment2D_Quad_Block::LoadSendBuffer_Flux_F2C(double *buffer,
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
        for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxS[i  ][k]+
                                   FluxS[i+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxN[i  ][k]+
                                   FluxN[i+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
     for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
        for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxW[j][k]+
                                   FluxW[j+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
     for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
        for ( k = 1 ; k <= NUM_VAR_ION5MOMENT2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxE[j][k]+
                                   FluxE[j+1][k]);
        } /* endfor */
     } /* endfor */
  } /* endif */
  return(0);
}

/*****************************************************************************************
 * Ion5Moment2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C -- Unloads receive message      *
 *                                                          buffer for fine to coarse    *
 *                                                          block message passing of     *
 *                                                          conservative solution fluxes.*
 *****************************************************************************************/
inline int Ion5Moment2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C(double *buffer,
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
  	buffer_count = buffer_count + NUM_VAR_ION5MOMENT2D;
        if (buffer_count >= buffer_size) return(1);
        FluxS[i] = -Ion5Moment2D_cState(buffer[buffer_count-3],
                                        buffer[buffer_count-2],
                                        buffer[buffer_count-1],
                                        buffer[buffer_count])
                   -FluxS[i];
     } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
  	buffer_count = buffer_count + NUM_VAR_ION5MOMENT2D;
        if (buffer_count >= buffer_size) return(1);
        FluxN[i] = -Ion5Moment2D_cState(buffer[buffer_count-3],
                                        buffer[buffer_count-2],
                                        buffer[buffer_count-1],
                                        buffer[buffer_count])
                   -FluxN[i];
     } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
     for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
  	buffer_count = buffer_count + NUM_VAR_ION5MOMENT2D;
        if (buffer_count >= buffer_size) return(1);
        FluxW[j] = -Ion5Moment2D_cState(buffer[buffer_count-3],
                                        buffer[buffer_count-2],
                                        buffer[buffer_count-1],
                                        buffer[buffer_count])
                   -FluxW[j];
     } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
     for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
  	buffer_count = buffer_count + NUM_VAR_ION5MOMENT2D;
        if (buffer_count >= buffer_size) return(1);
        FluxE[j] = -Ion5Moment2D_cState(buffer[buffer_count-3],
                                        buffer[buffer_count-2],
                                        buffer[buffer_count-1],
                                        buffer[buffer_count])
                   -FluxE[j];
     } /* endfor */
  } /* endif */
  return(0);
}

/********************************************************************************
 * Ion5Moment2D_Quad_Block -- Single Block External Subroutines.                *
 ********************************************************************************/

extern void Write_Solution_Block(Ion5Moment2D_Quad_Block &SolnBlk,
	                         ostream &Out_File);

extern void Read_Solution_Block(Ion5Moment2D_Quad_Block &SolnBlk,
	                        istream &In_File);

extern void Read_Neutral_Gas_Solution_Block(Ion5Moment2D_Quad_Block &SolnBlk,
		                            istream &In_File);

extern void Read_Electric_Field_Solution_Block(Ion5Moment2D_Quad_Block &SolnBlk,
		                               istream &In_File,
                                               const int Add_Initial_and_Solution_File_Electric_Fields);

extern void Broadcast_Solution_Block(Ion5Moment2D_Quad_Block &SolnBlk);

#ifdef _MPI_VERSION
extern void Broadcast_Solution_Block(Ion5Moment2D_Quad_Block &SolnBlk,
                                     MPI::Intracomm &Communicator, 
                                     const int Source_CPU);
#endif

extern void Copy_Solution_Block(Ion5Moment2D_Quad_Block &SolnBlk1,
		                Ion5Moment2D_Quad_Block &SolnBlk2);

extern int Prolong_Solution_Block(Ion5Moment2D_Quad_Block &SolnBlk_Fine,
				  Ion5Moment2D_Quad_Block &SolnBlk_Original,
				  const int Sector);

extern int Restrict_Solution_Block(Ion5Moment2D_Quad_Block &SolnBlk_Coarse,
				   Ion5Moment2D_Quad_Block &SolnBlk_Original_SW,
				   Ion5Moment2D_Quad_Block &SolnBlk_Original_SE,
				   Ion5Moment2D_Quad_Block &SolnBlk_Original_NW,
				   Ion5Moment2D_Quad_Block &SolnBlk_Original_NE);

extern void Output_Tecplot(Ion5Moment2D_Quad_Block &SolnBlk,
		           Ion5Moment2D_Input_Parameters &Input_Parameters,
                           const int Number_of_Time_Steps,
                           const double &Time,
                           const int Block_Number,
                           const int Output_Title,
	                   ostream &Out_File);

extern void Output_Cells_Tecplot(Ion5Moment2D_Quad_Block &SolnBlk,
		                 Ion5Moment2D_Input_Parameters &Input_Parameters,
                                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void ICs(Ion5Moment2D_Quad_Block &SolnBlk,
 	        const int i_ICtype,
                Ion5Moment2D_pState *Wo);

extern void ICs_Neutral_Gas(Ion5Moment2D_Quad_Block &SolnBlk,
 	                    const int i_ICtype,
                            Euler2D_pState &Wno);

extern void ICs_Electric_Field(Ion5Moment2D_Quad_Block &SolnBlk,
                               const int i_Electric_Field,
 	                       const double &Electric_Field_Strength,
 	                       const double &Electric_Field_Angle);

extern void BCs(Ion5Moment2D_Quad_Block &SolnBlk);

extern void BCs_Neutral_Gas(Ion5Moment2D_Quad_Block &SolnBlk);

extern void BCs_Electric_Field(Ion5Moment2D_Quad_Block &SolnBlk);

extern double CFL(Ion5Moment2D_Quad_Block &SolnBlk,
                  Ion5Moment2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(Ion5Moment2D_Quad_Block &SolnBlk, 
                                const double &Dt_min);

extern double L1_Norm_Residual(Ion5Moment2D_Quad_Block &SolnBlk);

extern double L2_Norm_Residual(Ion5Moment2D_Quad_Block &SolnBlk);

extern double Max_Norm_Residual(Ion5Moment2D_Quad_Block &SolnBlk);

extern void Linear_Reconstruction_GreenGauss(Ion5Moment2D_Quad_Block &SolnBlk,
                                             const int i,
                                             const int j,
					     const int Limiter);

extern void Linear_Reconstruction_GreenGauss(Ion5Moment2D_Quad_Block &SolnBlk,
					     const int Limiter);

extern void Linear_Reconstruction_LeastSquares(Ion5Moment2D_Quad_Block &SolnBlk,
                                               const int i,
                                               const int j,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares_2(Ion5Moment2D_Quad_Block &SolnBlk,
                                                 const int i,
                                                 const int j,
					         const int Limiter);

extern void Linear_Reconstruction_LeastSquares(Ion5Moment2D_Quad_Block &SolnBlk,
					       const int Limiter);

extern void Residual_Smoothing(Ion5Moment2D_Quad_Block &SolnBlk,
                               const int k_residual,
			       double &epsilon, 
                               const int number_of_Gauss_Seidel_iterations);

extern void Calculate_Refinement_Criteria(double *refinement_criteria,
					  Ion5Moment2D_Input_Parameters &IP,
                                          int &number_refinement_criteria,
                                          Ion5Moment2D_Quad_Block &SolnBlk);

extern void Fix_Refined_Block_Boundaries(Ion5Moment2D_Quad_Block SolnBlk,
                                         const int Fix_North_Boundary,
                                         const int Fix_South_Boundary,
                                         const int Fix_East_Boundary,
                                         const int Fix_West_Boundary);

extern void Unfix_Refined_Block_Boundaries(Ion5Moment2D_Quad_Block SolnBlk);

extern void Apply_Boundary_Flux_Corrections(Ion5Moment2D_Quad_Block SolnBlk,
                                            const int Number_Neighbours_North_Boundary,
                                            const int Number_Neighbours_South_Boundary,
                                            const int Number_Neighbours_East_Boundary,
                                            const int Number_Neighbours_West_Boundary);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Ion5Moment2D_Quad_Block SolnBlk,
                                                                const int i_stage,
                                                                const int n_stage,
                                                                const double &CFL_Number,
                                                                const int Time_Integration_Type,
                                                                const int Local_Time_Stepping,
                                                                const int Reconstruction_Type,
                                                                const int Limiter_Type,
                                                                const int Flux_Function_Type,
                                                                const int Number_Neighbours_North_Boundary,
                                                                const int Number_Neighbours_South_Boundary,
                                                                const int Number_Neighbours_East_Boundary,
                                                                const int Number_Neighbours_West_Boundary);

extern int dUdt_Residual_Evaluation(Ion5Moment2D_Quad_Block &SolnBlk,
				    Ion5Moment2D_Input_Parameters &Input_Parameters);

extern int dUdt_Multistage_Explicit(Ion5Moment2D_Quad_Block &SolnBlk,
   	                            const int i_stage,
                                    Ion5Moment2D_Input_Parameters &Input_Parameters);

extern int Update_Solution_Multistage_Explicit(Ion5Moment2D_Quad_Block &SolnBlk,
   	                                       const int i_stage,
                                               Ion5Moment2D_Input_Parameters &Input_Parameters);

extern void dUdt_Output_Cells_Tecplot(Ion5Moment2D_Quad_Block &SolnBlk,
                                      const int Number_of_Time_Steps,
                                      const double &Time,
                                      const int Block_Number,
                                      const int Output_Title,
                                      const int Reconstruction_Type,
                                      const int Limiter_Type,
	                              const int Flux_Function_Type,
	                              ostream &Out_File);

/********************************************************************************
 * Ion5Moment2D_Quad_Block -- Multiple Block External Subroutines.              *
 ********************************************************************************/

extern Ion5Moment2D_Quad_Block* Allocate(Ion5Moment2D_Quad_Block *Soln_ptr,
                                         Ion5Moment2D_Input_Parameters &Input_Parameters);

extern Ion5Moment2D_Quad_Block* Deallocate(Ion5Moment2D_Quad_Block *Soln_ptr,
                                           Ion5Moment2D_Input_Parameters &Input_Parameters);

extern void ICs(Ion5Moment2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
                Ion5Moment2D_Input_Parameters &Input_Parameters);

extern int Read_Restart_Solution(Ion5Moment2D_Quad_Block *Soln_ptr,
                                 AdaptiveBlock2D_List &Soln_Block_List,
                                 Ion5Moment2D_Input_Parameters &Input_Parameters,
		                 int &Number_of_Time_Steps,
                                 double &Time,
                                 CPUTime &CPU_Time);

extern int Write_Restart_Solution(Ion5Moment2D_Quad_Block *Soln_ptr,
                                  AdaptiveBlock2D_List &Soln_Block_List,
                                  Ion5Moment2D_Input_Parameters &Input_Parameters,
		                  const int Number_of_Time_Steps,
                                  const double &Time,
                                  const CPUTime &CPU_Time);

extern int Read_Neutral_Gas_Solution(Ion5Moment2D_Quad_Block *Soln_ptr,
                                     AdaptiveBlock2D_List &Soln_Block_List,
                                     Ion5Moment2D_Input_Parameters &Input_Parameters);

extern int Read_Electric_Field_Solution(Ion5Moment2D_Quad_Block *Soln_ptr,
                                        AdaptiveBlock2D_List &Soln_Block_List,
                                        Ion5Moment2D_Input_Parameters &Input_Parameters);

extern int Output_Tecplot(Ion5Moment2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
			  Ion5Moment2D_Input_Parameters &Input_Parameters,
		          const int Number_of_Time_Steps,
                          const double &Time);

extern int Output_Cells_Tecplot(Ion5Moment2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                Ion5Moment2D_Input_Parameters &Input_Parameters,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Mesh_Tecplot(Ion5Moment2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Ion5Moment2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Mesh_Gnuplot(Ion5Moment2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Ion5Moment2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern void BCs(Ion5Moment2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
		Ion5Moment2D_Input_Parameters &Input_Parameters);

extern void BCs_Neutral_Gas(Ion5Moment2D_Quad_Block *Soln_ptr,
                            AdaptiveBlock2D_List &Soln_Block_List);

extern void BCs_Electric_Field(Ion5Moment2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double CFL(Ion5Moment2D_Quad_Block *Soln_ptr,
                  AdaptiveBlock2D_List &Soln_Block_List,
                  Ion5Moment2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(Ion5Moment2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List, 
                                const double &Dt_min);

extern double L1_Norm_Residual(Ion5Moment2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
			       Ion5Moment2D_Input_Parameters &Input_Parameters);

extern double L2_Norm_Residual(Ion5Moment2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
			       Ion5Moment2D_Input_Parameters &Input_Parameters);

extern double Max_Norm_Residual(Ion5Moment2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
				Ion5Moment2D_Input_Parameters &Input_Parameters);

extern void Residual_Smoothing(Ion5Moment2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Ion5Moment2D_Input_Parameters &Input_Parameters,
   	                       const int I_Stage);

extern void Apply_Boundary_Flux_Corrections(Ion5Moment2D_Quad_Block *Soln_ptr,
                                            AdaptiveBlock2D_List &Soln_Block_List,
                                            Ion5Moment2D_Input_Parameters &Input_Parameters);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Ion5Moment2D_Quad_Block *Soln_ptr,
                                                                AdaptiveBlock2D_List &Soln_Block_List,
                                                                Ion5Moment2D_Input_Parameters &Input_Parameters,
   	                                                        const int I_Stage);

extern int dUdt_Multistage_Explicit(Ion5Moment2D_Quad_Block *Soln_ptr,
				    AdaptiveBlockResourceList &Global_Soln_Block_List,
                                    AdaptiveBlock2D_List &Local_Soln_Block_List,
                                    Ion5Moment2D_Input_Parameters &Input_Parameters,
   	                            const int I_Stage);

extern int Update_Solution_Multistage_Explicit(Ion5Moment2D_Quad_Block *Soln_ptr,
                                               AdaptiveBlock2D_List &Soln_Block_List,
                                               Ion5Moment2D_Input_Parameters &Input_Parameters,
   	                                       const int I_Stage);

extern int dUdt_Output_Cells_Tecplot(Ion5Moment2D_Quad_Block *Soln_ptr,
                                     AdaptiveBlock2D_List &Soln_Block_List,
                                     Ion5Moment2D_Input_Parameters &Input_Parameters,
		                     const int Number_of_Time_Steps,
                                     const double &Time);

/********************************************************************************
 * Ion5Moment2D_Quad_Block -- Multiple Block External Subroutines for Mesh.     *
 ********************************************************************************/

extern Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                            Ion5Moment2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Broadcast_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                      Ion5Moment2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                             Ion5Moment2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                                            Ion5Moment2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                  Ion5Moment2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                 Ion5Moment2D_Input_Parameters &Input_Parameters);

extern int Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                          Ion5Moment2D_Input_Parameters &Input_Parameters);

extern int Output_Nodes_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                Ion5Moment2D_Input_Parameters &Input_Parameters);

extern int Output_Cells_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                Ion5Moment2D_Input_Parameters &Input_Parameters);

/********************************************************************************
 * Ion5Moment2D_Quad_Block -- Solvers.                                          *
 ********************************************************************************/

extern int Ion5Moment2DQuadSolver(char *Input_File_Name_ptr,
                                  int batch_flag);

#endif /* _ION5MOMENT2D_QUAD_INCLUDED  */
