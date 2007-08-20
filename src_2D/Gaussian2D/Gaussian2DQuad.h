/* Gaussian2DQuad.h:  Header file defining 
                      2D Gaussian Quadrilateral Mesh Solution Classes. */

#ifndef _GAUSSIAN2D_QUAD_INCLUDED
#define _GAUSSIAN2D_QUAD_INCLUDED

/* Include 2D Gaussian state, 2D cell, 2D quadrilateral multiblock 
   grid, quadtree, AMR, and 2D Gaussian input header files. */

#ifndef _GAUSSIAN2D_STATE_INCLUDED
#include "Gaussian2DState.h"
#endif // _GAUSSIAN2D_STATE_INCLUDED

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

#ifndef _GAUSSIAN2D_INPUT_INCLUDED
#include "Gaussian2DInput.h"
#endif // _GAUSSIAN2D_INPUT_INCLUDED

#ifndef _GAUSSIAN2D_INPUT_INCLUDED
#include "Gaussian2DInput.h"
#endif // _GAUSSIAN2D_INPUT_INCLUDED

#ifndef _SYSTEM_LINUX_INCLUDED
#include "../System/System_Linux.h"
#endif // _GAUSSIAN2D_INPUT_INCLUDED

/* Define the structures and classes. */

#define	NUMBER_OF_RESIDUAL_VECTORS_GAUSSIAN2D    3

/********************************************************
 * Class: Gaussian2D_Quad_Block                         *
 *                                                      *
 * Member functions                                     *
 *       W      -- Return primitive variable solution   *
 *                 for the block (cell average).        *
 *       U      -- Return conserved variable solution   *
 *                 for the block (cell average).        *
 *    Grid      -- Return the solution block            *
 *                 quadrilateral grid or mesh.          *
 *      dt      -- Return local time step for the       *
 *                 solution block.                      *
 *    dUdt      -- Return the solution block residuals. *
 *    dWdx      -- Return the unlimited primitive       *
 *                 variable solution gradients          *
 *                 (x-direction) for the block.         *
 *    dWdy      -- Return the unlimited primitive       *
 *                 variable solution gradients          *
 *                 (y-direction) for the block.         *
 *     phi      -- Return the solution slope limiters.  *
 *      Uo      -- Return initial solution state.       *
 *    FluxN     -- Return array of north boundary       *
 *                 solution fluxes.                     *
 *    FluxS     -- Return array of south boundary       *
 *                 solution fluxes.                     *
 *    FluxE     -- Return array of east boundary        *
 *                 solution fluxes.                     *
 *    FluxW     -- Return array of west boundary        *
 *                 solution fluxes.                     *
 *     NCi      -- Return number of cells in            *
 *                 the i-direction (zeta-direction).    *
 *     ICl      -- Return lower index for cells in      *
 *                 the i-direction (zeta-direction).    *
 *     ICu      -- Return upper index for cells in      *
 *                 the i-direction (zeta-direction).    *
 *     NCj      -- Return number of cells in            *
 *                 the j-direction (eta-direction).     *
 *     JCl      -- Return lower index for cells in      *
 *                 the j-direction (eta-direction).     *
 *     JCu      -- Return upper index for cells in      *
 *                 the j-direction (eta-direction).     *
 *     Nghost   -- Return number of ghost (halo or      *
 *                 overlap) cells.                      *
 * Axisymmetric -- Return axisymmetric flow             *
 *                 indicator (=1 for axisymmetric flow, *
 *                            =0 for planar flow).      *
 * Freeze_Limiter -- Return limiter freezing indicator  *
 *                 (=1 for limiter freezing on,         *
 *                 (=0 for limiter freezing off).       * 
 *     WoN      -- Return array of reference states for *
 *                 the application of boundary          *
 *                 conditions at the north boundary of  *
 *                 the solution block.                  *
 *     WoS      -- Return array of reference states for *
 *                 the application of boundary          *
 *                 conditions at the south boundary of  *
 *                 the solution block.                  *
 *     WoE      -- Return array of reference states for *
 *                 the application of boundary          *
 *                 conditions at the east boundary of   *
 *                 the solution block.                  *
 *     WoW      -- Return array of reference states for *
 *                 the application of boundary          *
 *                 conditions at the west boundary of   *
 *                 the solution block.                  *
 *   allocate   -- Allocate memory for structured       *
 *                 quadrilateral solution block.        *
 *   deallocate -- Deallocate memory for structured     *
 *                 quadrilateral solution block.        *
 *      Wn      -- Return primitive variable solution   *
 *                 at the specified node.               *
 *      Un      -- Return conserved variable solution   *
 *                 at the specified node.               *
 *    WnNW      -- Return primitive variable solution   *
 *                 at the north-west node.              *
 *    WnNE      -- Return primitive variable solution   *
 *                 at the north-east node.              *
 *    WnSW      -- Return primitive variable solution   *
 *                 at the south-west node.              *
 *    WnSE      -- Return primitive variable solution   *
 *                 at the south-east node.              *
 *    UnNW      -- Return conserved variable solution   *
 *                 at the north-west node.              *
 *    UnNE      -- Return conserved variable solution   *
 *                 at the north-east node.              *
 *    UnSW      -- Return conserved variable solution   *
 *                 at the south-west node.              *
 *    UnSE      -- Return conserved variable solution   *
 *                 at the south-east node.              *
 * evaluate_limiters -- Set flag to evaluate limiters.  *
 * freeze_limiters -- Set flag to freeze limiters.      *
 * Member functions required for message passing.       *
 *    NumVar    -- Returns number of solution variables *
 *                 in primitive and conserved solution  *
 *                 state vectors.                       *
 * LoadSendBuffer -- Loads send buffer.                 *
 * LoadSendBuffer_F2C -- Loads send buffer for fine to  *
 *                        coarse block messages.        *
 * UnloadReceiveBuffer -- Unloads receive buffer.       *
 * UnloadReceiveBuffer_F2C -- Unloads receive buffer    *
 *                            for fine to coarse block  *
 *                            messages.                 *
 * SubcellReconstruction -- Performs subcell solution   *
 *                          reconstruction used in      *
 *                          adaptive mesh refinement.   *
 * LoadSendBuffer_Flux_F2C -- Loads send buffer for     *
 *                            sending conservative flux *
 *                            corrections from fine to  *
 *                            coarse solution blocks.   *
 * UnLoadSendBuffer_Flux_F2C -- Loads send buffer for   *
 *                            sending conservative flux *
 *                            corrections from fine to  *
 *                            coarse solution blocks.   *
 *                                                      *
 * Member operators                                     *
 *      S -- a 2D Gaussian solution                     *
 *                                                      *
 * S = S;                                               *
 * cout << S; (output function)                         *
 * cin  >> S; (input function)                          *
 *                                                      *
 ********************************************************/
class Gaussian2D_Quad_Block{
  private:
  public:
    Gaussian2D_pState           **W;  // Primitive solution state.
    Gaussian2D_cState           **U;  // Conserved solution state.
    Grid2D_Quad_Block          Grid;  // 2D quadrilateral grid geometry.
    double                     **dt;  // Local time step.
    Gaussian2D_cState       ***dUdt;  // Solution residual.
    Gaussian2D_pState        **dWdx;  // Unlimited solution gradient
                                      // (x-direction).
    Gaussian2D_pState        **dWdy;  // Unlimited solution gradient
                                      // (y-direction).
    Gaussian2D_pState         **phi;  // Solution slope limiter.
    Gaussian2D_cState          **Uo;  // Initial solution state.
    Gaussian2D_cState *FluxN,*FluxS,  // Boundary solution fluxes.
                      *FluxE,*FluxW;
    int                 NCi,ICl,ICu;  // i-direction cell counters.
    int                 NCj,JCl,JCu;  // j-direction cell counters.
    int                     Nghost;   // Number of ghost cells.
    int                Axisymmetric;  // Axisymmetric flow indicator.
    int              Freeze_Limiter;  // Limiter freezing indicator.
    Gaussian2D_pState     *WoN,*WoS,  // Boundary condition reference states for
                          *WoE,*WoW;  // north, south, east, & west boundaries.
    static int    residual_variable;  // Static integer that indicates which variable is used for residual calculations.
    static int            Flow_Type;  // Static integer identifying the flow type (inviscid?)
                                      // Made public so can access them.
		      
    /* Creation, copy, and assignment constructors. */
    Gaussian2D_Quad_Block(void) {
      NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
       W = NULL; U = NULL; dt = NULL; dUdt = NULL; 
       dWdx = NULL; dWdy = NULL; phi = NULL; Uo = NULL;
       FluxN = NULL; FluxS = NULL; FluxE = NULL; FluxW = NULL;
       WoN = NULL; WoS = NULL; WoE = NULL; WoW = NULL;
       Axisymmetric = 0; Freeze_Limiter = OFF;
    }

    Gaussian2D_Quad_Block(const Gaussian2D_Quad_Block &Soln) {
       NCi = Soln.NCi; ICl = Soln.ICl; ICu = Soln.ICu; 
       NCj = Soln.NCj; JCl = Soln.JCl; JCu = Soln.JCu; Nghost = Soln.Nghost;
       Grid = Soln.Grid; W = Soln.W; U = Soln.U; dt = Soln.dt; dUdt = Soln.dUdt; 
       dWdx = Soln.dWdx; dWdy = Soln.dWdy; phi = Soln.phi; Uo = Soln.Uo;
       FluxN = Soln.FluxN; FluxS = Soln.FluxS; FluxE = Soln.FluxE; FluxW = Soln.FluxW;
       WoN = Soln.WoN; WoS = Soln.WoS; WoE = Soln.WoE; WoW = Soln.WoW;
       Axisymmetric = 0; Freeze_Limiter = Soln.Freeze_Limiter;
    }

    /* Destructor. */
    // ~Gaussian2D_Quad_Block(void);
    // Use automatically generated destructor.

    /* Assignment operator. */
    // Gaussian2D_Quad_Block operator = (const Gaussian2D_Quad_Block &Soln);
    // Use automatically generated assignment operator.

    /* Allocate memory for structured quadrilateral solution block. */
    void allocate(const int Ni, const int Nj, const int Ng);

    /* Deallocate memory for structured quadrilateral solution block. */
    void deallocate(void);

    /* Return primitive solution state at specified node. */
    Gaussian2D_pState Wn(const int ii, const int jj);

    /* Retern conserverd solution state at specified node. */
    Gaussian2D_cState Un(const int ii, const int jj);

    /* Return primitive solution state at cell nodes. */
    Gaussian2D_pState WnNW(const int ii, const int jj);
    Gaussian2D_pState WnNE(const int ii, const int jj);
    Gaussian2D_pState WnSE(const int ii, const int jj);
    Gaussian2D_pState WnSW(const int ii, const int jj);

    /* Return conserved solution state at cell nodes. */
    Gaussian2D_cState UnNW(const int ii, const int jj);
    Gaussian2D_cState UnNE(const int ii, const int jj);
    Gaussian2D_cState UnSE(const int ii, const int jj);
    Gaussian2D_cState UnSW(const int ii, const int jj);

    /* Set flags for limiter evaluation. */
    void evaluate_limiters(void);
    void freeze_limiters(void);

    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file,
				 const Gaussian2D_Quad_Block &Soln);
    friend istream &operator >> (istream &in_file,
				 Gaussian2D_Quad_Block &Soln);

    /* MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING. */
    /* Number of solution state variables. */
    int NumVar(void);
    /* Load send message passing buffer. */
    int LoadSendBuffer(double *buffer,
                       int &buffer_count,
                       const int buffer_size,
                       const int i_min, 
                       const int i_max,
                       const int i_inc,
                       const int j_min, 
                       const int j_max,
                       const int j_inc);
    int LoadSendBuffer_F2C(double *buffer,
                           int &buffer_count,
                           const int buffer_size,
                           const int i_min, 
                           const int i_max,
                           const int i_inc,
                           const int j_min, 
                           const int j_max,
                           const int j_inc);
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
    /* Unload receive message passing buffer. */
    int UnloadReceiveBuffer(double *buffer,
                            int &buffer_count,
                            const int buffer_size,
                            const int i_min, 
                            const int i_max,
                            const int i_inc,
                            const int j_min, 
                            const int j_max,
                            const int j_inc);
    int UnloadReceiveBuffer_F2C(double *buffer,
                                int &buffer_count,
                                const int buffer_size,
                                const int i_min, 
                                const int i_max,
                                const int i_inc,
                                const int j_min, 
                                const int j_max,
                                const int j_inc);
    int UnloadReceiveBuffer_C2F(double *buffer,
                                int &buffer_count,
                                const int buffer_size,
                                const int i_min, 
                                const int i_max,
                                const int i_inc,
                                const int j_min, 
                                const int j_max,
                                const int j_inc);
    /* Subcell solution reconstruction within given computational cell. */
    void SubcellReconstruction(const int i,
                               const int j,
                               const int Limiter);
    /* Load and unload conservative flux message passing buffer. */
    int LoadSendBuffer_Flux_F2C(double *buffer,
                                int &buffer_count,
                                const int buffer_size,
                                const int i_min, 
                                const int i_max,
                                const int i_inc,
                                const int j_min, 
                                const int j_max,
                                const int j_inc);
    int UnloadReceiveBuffer_Flux_F2C(double *buffer,
                                     int &buffer_count,
                                     const int buffer_size,
                                     const int i_min, 
                                     const int i_max,
                                     const int i_inc,
                                     const int j_min, 
                                     const int j_max,
                                     const int j_inc);

};

/**************************************************************************
 * Gaussian2D_Quad_Block::allocate -- Allocate memory.                    *
 **************************************************************************/
inline void Gaussian2D_Quad_Block::allocate(const int Ni, const int Nj, const int Ng) {
   int i, j, k; 
   assert(Ni > 1 && Nj > 1 && Ng > 1); 
   Grid.allocate(Ni, Nj, Ng);
   NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1; 
   NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1; 
   Nghost = Ng;

   W = new Gaussian2D_pState*[NCi]; 
   U = new Gaussian2D_cState*[NCi];
   dt = new double*[NCi]; 
   dUdt = new Gaussian2D_cState**[NCi];
   dWdx = new Gaussian2D_pState*[NCi]; 
   dWdy = new Gaussian2D_pState*[NCi];
   phi = new Gaussian2D_pState*[NCi]; 
   Uo = new Gaussian2D_cState*[NCi];

   for ( i = 0; i <= NCi-1 ; ++i ) {
      W[i] = new Gaussian2D_pState[NCj]; 
      U[i] = new Gaussian2D_cState[NCj];
      dt[i] = new double[NCj]; 
      dUdt[i] = new Gaussian2D_cState*[NCj];
      for ( j = 0; j <= NCj-1 ; ++j ) { 
	dUdt[i][j] = new Gaussian2D_cState[NUMBER_OF_RESIDUAL_VECTORS_GAUSSIAN2D];
      }
      dWdx[i] = new Gaussian2D_pState[NCj]; 
      dWdy[i] = new Gaussian2D_pState[NCj];
      phi[i] = new Gaussian2D_pState[NCj]; 
      Uo[i] = new Gaussian2D_cState[NCj];
   } /* endfor */

   FluxN = new Gaussian2D_cState[NCi]; 
   FluxS = new Gaussian2D_cState[NCi];
   FluxE = new Gaussian2D_cState[NCj]; 
   FluxW = new Gaussian2D_cState[NCj];
   WoN = new Gaussian2D_pState[NCi]; 
   WoS = new Gaussian2D_pState[NCi];
   WoE = new Gaussian2D_pState[NCj]; 
   WoW = new Gaussian2D_pState[NCj];
   // Set the solution residuals, gradients, limiters, and other values to zero.
   for (j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
      for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
          for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_GAUSSIAN2D-1 ; ++k ) {
	     dUdt[i][j][k] = Gaussian2D_U_VACUUM;
          } /* endfor */
	  dWdx[i][j] = Gaussian2D_W_VACUUM; 
          dWdy[i][j] = Gaussian2D_W_VACUUM;
	  phi[i][j] = Gaussian2D_W_VACUUM; 
          Uo[i][j] = Gaussian2D_U_VACUUM;
	  dt[i][j] = ZERO;
      } /* endfor */
   } /* endfor */
}

/**************************************************************************
 * Gaussian2D_Quad_Block::deallocate -- Deallocate memory.                *
 **************************************************************************/
inline void Gaussian2D_Quad_Block::deallocate(void) {
   int i, j; 
   Grid.deallocate(); 
   for ( i = 0; i <= NCi-1 ; ++i ) {
      delete []W[i]; 
      W[i] = NULL; 
      delete []U[i]; 
      U[i] = NULL;
      delete []dt[i]; 
      dt[i] = NULL; 

      for ( j = 0; j <= NCj-1 ; ++j ) { 
	delete []dUdt[i][j]; 
        dUdt[i][j] = NULL;
      }

      delete []dUdt[i]; 
      dUdt[i] = NULL;
      delete []dWdx[i]; 
      dWdx[i] = NULL; 
      delete []dWdy[i]; 
      dWdy[i] = NULL;
      delete []phi[i]; 
      phi[i] = NULL; 
      delete []Uo[i]; 
      Uo[i] = NULL;
   } /* endfor */

   delete []W; 
   W = NULL; 
   delete []U; 
   U = NULL;
   delete []dt; 
   dt = NULL; 
   delete []dUdt; 
   dUdt = NULL;
   delete []dWdx; 
   dWdx = NULL; 
   delete []dWdy; 
   dWdy = NULL; 
   delete []phi; 
   phi = NULL; 
   delete []Uo; 
   Uo = NULL;
   delete []FluxN; 
   FluxN = NULL; 
   delete []FluxS; 
   FluxS = NULL;
   delete []FluxE; 
   FluxE = NULL; 
   delete []FluxW; 
   FluxW = NULL;
   delete []WoN; 
   WoN = NULL; 
   delete []WoS; 
   WoS = NULL;
   delete []WoE; 
   WoE = NULL; 
   delete []WoW; 
   WoW = NULL;

   NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
}

/**************************************************************************
 * Gaussian2D_Quad_Block::Wn -- Node primitive variable solution state.   *
 **************************************************************************/
inline Gaussian2D_pState Gaussian2D_Quad_Block::Wn(const int ii, const int jj) {
    double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
           eta1, zeta1, eta2, zeta2, eta, zeta;
    double micro_toler = TOLER/THOUSAND;
    Gaussian2D_pState A, B, C, D;
    x=Grid.Node[ii][jj].X.x; 
    y=Grid.Node[ii][jj].X.y;
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
    if (fabs(aa) < micro_toler*micro_toler) {
      if (fabs(bb) >= micro_toler*micro_toler) { zeta1=-cc/bb; }
      else { zeta1 = -cc/sgn(bb)*(micro_toler*micro_toler); } 
      if (fabs(cy+dy*zeta1) >= micro_toler*micro_toler) { eta1=(y-ay-by*zeta1)/(cy+dy*zeta1); } 
      else { eta1 = HALF; } zeta2=zeta1; eta2=eta1;
    } else {
      if (bb*bb-FOUR*aa*cc >= micro_toler*micro_toler) { zeta1=HALF*(-bb+sqrt(bb*bb-FOUR*aa*cc))/aa; }
      else { zeta1 = -HALF*bb/aa; } 
      if (fabs(cy+dy*zeta1) < micro_toler*micro_toler) { eta1=-ONE; } 
      else { eta1=(y-ay-by*zeta1)/(cy+dy*zeta1); } 
      if (bb*bb-FOUR*aa*cc >= micro_toler*micro_toler) { zeta2=HALF*(-bb-sqrt(bb*bb-FOUR*aa*cc))/aa; }
      else { zeta2 = -HALF*bb/aa; }
      if (fabs(cy+dy*zeta2) < micro_toler*micro_toler) { eta2=-ONE; } 
      else { eta2=(y-ay-by*zeta2)/(cy+dy*zeta2); }
    } /* end if */
    if (zeta1 > -micro_toler && zeta1 < ONE + micro_toler &&
        eta1  > -micro_toler && eta1  < ONE + micro_toler) {
      zeta=zeta1; eta=eta1;
    } else if (zeta2 > -micro_toler && zeta2 < ONE + micro_toler &&
               eta2  > -micro_toler && eta2  < ONE + micro_toler) {
      zeta=zeta2; eta=eta2;
    } else {
      zeta=HALF; eta=HALF;
    } /* endif */
    A=W[ii-1][jj-1]; B=W[ii-1][jj]-W[ii-1][jj-1]; C=W[ii][jj-1]-W[ii-1][jj-1];
    D=W[ii][jj]+W[ii-1][jj-1]-W[ii-1][jj]-W[ii][jj-1];
    return (A+B*zeta+C*eta+D*zeta*eta);
}

/**************************************************************************
 * Gaussian2D_Quad_Block::Un -- Node conserved variable solution state.   *
 **************************************************************************/
inline Gaussian2D_cState Gaussian2D_Quad_Block::Un(const int ii, const int jj) {
    double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
           eta1, zeta1, eta2, zeta2, eta, zeta;
    double micro_toler = TOLER/THOUSAND;
    Gaussian2D_cState A, B, C, D;
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
    if (fabs(aa) < micro_toler*micro_toler) {
      if (fabs(bb) >= micro_toler*micro_toler) { zeta1=-cc/bb; }
      else { zeta1 = -cc/sgn(bb)*(micro_toler*micro_toler); } 
      if (fabs(cy+dy*zeta1) >= micro_toler*micro_toler) { eta1=(y-ay-by*zeta1)/(cy+dy*zeta1); } 
      else { eta1 = HALF; } zeta2=zeta1; eta2=eta1;
    } else {
      if (bb*bb-FOUR*aa*cc >= micro_toler*micro_toler) { zeta1=HALF*(-bb+sqrt(bb*bb-FOUR*aa*cc))/aa; }
      else { zeta1 = -HALF*bb/aa; } 
      if (fabs(cy+dy*zeta1) < micro_toler*micro_toler) { eta1=-ONE; } 
      else { eta1=(y-ay-by*zeta1)/(cy+dy*zeta1); } 
      if (bb*bb-FOUR*aa*cc >= micro_toler*micro_toler) { zeta2=HALF*(-bb-sqrt(bb*bb-FOUR*aa*cc))/aa; }
      else { zeta2 = -HALF*bb/aa; }
      if (fabs(cy+dy*zeta2) < micro_toler*micro_toler) { eta2=-ONE; } 
      else { eta2=(y-ay-by*zeta2)/(cy+dy*zeta2); }
    } /* end if */
    if (zeta1 > -micro_toler && zeta1 < ONE + micro_toler &&
        eta1  > -micro_toler && eta1  < ONE + micro_toler) {
      zeta=zeta1; eta=eta1;
    } else if (zeta2 > -micro_toler && zeta2 < ONE + micro_toler &&
               eta2  > -micro_toler && eta2  < ONE + micro_toler) {
      zeta=zeta2; eta=eta2;
    } else {
      zeta=HALF; eta=HALF;
    } /* endif */
    A=U[ii-1][jj-1]; B=U[ii-1][jj]-U[ii-1][jj-1]; C=U[ii][jj-1]-U[ii-1][jj-1];
    D=U[ii][jj]+U[ii-1][jj-1]-U[ii-1][jj]-U[ii][jj-1];
    return (A+B*zeta+C*eta+D*zeta*eta);
}

/**************************************************************************
 * Gaussian2D_Quad_Block::Wn?? -- Get cell node primitive solution states.*
 **************************************************************************/
inline Gaussian2D_pState Gaussian2D_Quad_Block::WnNW(const int ii, const int jj) {
  return (Wn(ii, jj+1));
}

inline Gaussian2D_pState Gaussian2D_Quad_Block::WnNE(const int ii, const int jj) {
  return (Wn(ii+1, jj+1));
}

inline Gaussian2D_pState Gaussian2D_Quad_Block::WnSE(const int ii, const int jj) {
  return (Wn(ii+1, jj));
}

inline Gaussian2D_pState Gaussian2D_Quad_Block::WnSW(const int ii, const int jj) {
  return (Wn(ii, jj));
}

/**************************************************************************
 * Gaussian2D_Quad_Block::Un?? -- Get cell node conserved solution states.*
 **************************************************************************/
inline Gaussian2D_cState Gaussian2D_Quad_Block::UnNW(const int ii, const int jj) {
  return (Un(ii, jj+1));
}

inline Gaussian2D_cState Gaussian2D_Quad_Block::UnNE(const int ii, const int jj) {
  return (Un(ii+1, jj+1));
}

inline Gaussian2D_cState Gaussian2D_Quad_Block::UnSE(const int ii, const int jj) {
  return (Un(ii+1, jj));
}

inline Gaussian2D_cState Gaussian2D_Quad_Block::UnSW(const int ii, const int jj) {
  return (Un(ii, jj));
}

/*****************************************************************************
 * Gaussian2D_Quad_Block::evaluate_limiters -- Set flag to evaluate limiters.*
 *****************************************************************************/
inline void Gaussian2D_Quad_Block::evaluate_limiters(void) {
  Freeze_Limiter = OFF; 
}

/*****************************************************************************
 * Gaussian2D_Quad_Block::freeze_limiters -- Set flag to freeze limiters.    *
 *****************************************************************************/
inline void Gaussian2D_Quad_Block::freeze_limiters(void) {
  Freeze_Limiter = ON; 
}

/**************************************************************************
 * Gaussian2D_Quad_Block -- Input-output operators.                       *
 **************************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Gaussian2D_Quad_Block &SolnBlk) {
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
			     Gaussian2D_Quad_Block &SolnBlk) {
  int i, j, k, ni, il, iu, ng, nj, jl, ju;
  Grid2D_Quad_Block New_Grid; in_file >> New_Grid; 
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
  Copy_Quad_Block(SolnBlk.Grid, New_Grid); New_Grid.deallocate();
  for ( j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
     for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
         in_file >> SolnBlk.U[i][j];
         SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);
         for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_GAUSSIAN2D-1 ; ++k ) {
	     SolnBlk.dUdt[i][j][k] = Gaussian2D_U_VACUUM;
         } /* endfor */
	 SolnBlk.dWdx[i][j] = Gaussian2D_W_VACUUM;
	 SolnBlk.dWdy[i][j] = Gaussian2D_W_VACUUM;
	 SolnBlk.phi[i][j] = Gaussian2D_W_VACUUM;
	 SolnBlk.Uo[i][j] = Gaussian2D_U_VACUUM;
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

/*******************************************************************************
 *                                                                             *
 * MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING.                              *
 *                                                                             *
 *******************************************************************************/

/*******************************************************************************
 * Gaussian2D_Quad_Block::NumVar -- Returns number of state variables.         *
 *******************************************************************************/
inline int Gaussian2D_Quad_Block::NumVar(void) {
  return (int(NUM_VAR_GAUSSIAN2D));
}

/*******************************************************************************
 * Gaussian2D_Quad_Block::LoadSendBuffer -- Loads send message buffer.         *
 *******************************************************************************/
inline int Gaussian2D_Quad_Block::LoadSendBuffer(double *buffer,
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
        for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = U[i][j][k];
        } /* endfor */
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * Gaussian2D_Quad_Block::LoadSendBuffer_F2C -- Loads send message buffer for  *
 *                                           fine to coarse block message      *
 *                                           passing.                          *
 *******************************************************************************/
inline int Gaussian2D_Quad_Block::LoadSendBuffer_F2C(double *buffer,
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
        for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (Grid.Cell[i  ][j  ].A*U[i  ][j  ][k]+
                                   Grid.Cell[i+1][j  ].A*U[i+1][j  ][k]+
                                   Grid.Cell[i  ][j+1].A*U[i  ][j+1][k]+
                                   Grid.Cell[i+1][j+1].A*U[i+1][j+1][k])/
                                  (Grid.Cell[i  ][j  ].A+
                                   Grid.Cell[i+1][j  ].A+
                                   Grid.Cell[i  ][j+1].A+
                                   Grid.Cell[i+1][j+1].A);
/*            buffer[buffer_count] = (Grid.Cell[i  ][j  ].A*U[i  ][j  ][k]+ */
/*                                    Grid.Cell[i+1][j  ].A*U[i+1][j  ][k]+ */
/*                                    Grid.Cell[i  ][j+1].A*U[i  ][j+1][k]+ */
/*                                    Grid.Cell[i+1][j+1].A*U[i+1][j+1][k]); */

        } /* endfor */
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * Gaussian2D_Quad_Block::LoadSendBuffer_C2F -- Loads send message buffer for  *
 *                                           coarse to fine block message      *
 *                                           passing.                          *
 *******************************************************************************/
inline int Gaussian2D_Quad_Block::LoadSendBuffer_C2F(double *buffer,
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
  Gaussian2D_pState Wfine;
  Gaussian2D_cState Ufine;

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
	    dX = (Grid.Node[i][j].X+
		  HALF*(Grid.Node[i][j].X+Grid.Node[i+1][j].X)+
		  HALF*(Grid.Node[i][j].X+Grid.Node[i][j+1].X)+
		  Grid.Cell[i][j].Xc)/FOUR -
                 Grid.Cell[i][j].Xc;
	    Wfine = W[i][j] + (phi[i][j]^dWdx[i][j])*dX.x +
                              (phi[i][j]^dWdy[i][j])*dX.y;
	    Ufine = Wfine.U();
	    for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
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
	    dX = (HALF*(Grid.Node[i][j].X+Grid.Node[i+1][j].X)+
		  Grid.Node[i+1][j].X+
		  Grid.Cell[i][j].Xc+
		  HALF*(Grid.Node[i+1][j].X+Grid.Node[i+1][j+1].X))/FOUR -
                 Grid.Cell[i][j].Xc;
	    Wfine = W[i][j] + (phi[i][j]^dWdx[i][j])*dX.x +
                              (phi[i][j]^dWdy[i][j])*dX.y;
	    Ufine = Wfine.U();
	    for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
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
	    dX = (HALF*(Grid.Node[i][j].X+Grid.Node[i][j+1].X)+
		  Grid.Cell[i][j].Xc+
		  Grid.Node[i][j+1].X+
		  HALF*(Grid.Node[i][j+1].X+Grid.Node[i+1][j+1].X))/FOUR -
                 Grid.Cell[i][j].Xc;
	    Wfine = W[i][j] + (phi[i][j]^dWdx[i][j])*dX.x +
                              (phi[i][j]^dWdy[i][j])*dX.y;
	    Ufine = Wfine.U();
	    for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) { 
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
	    dX = (Grid.Cell[i][j].Xc+
		  HALF*(Grid.Node[i+1][j].X+Grid.Node[i+1][j+1].X)+
		  HALF*(Grid.Node[i][j+1].X+Grid.Node[i+1][j+1].X)+
		  Grid.Node[i+1][j+1].X)/FOUR -
	         Grid.Cell[i][j].Xc;
	    Wfine = W[i][j] + (phi[i][j]^dWdx[i][j])*dX.x +
                              (phi[i][j]^dWdy[i][j])*dX.y;
	    Ufine = Wfine.U();
	    for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
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
              Wfine = W[i_min][j] +
                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
              Ufine = Wfine.U();
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } /* endif */
  } /* endif */
  return(0);

//  if (j_min == j_max) { // North or south boundary.
//     // Four different orderings to consider depending on the value of i_inc & j_inc.
//     if (j_inc > 0) {
//        if (i_inc > 0) {
//           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
//              // Perform limited linear least squares reconstruction in cell (i, j_min).
//              SubcellReconstruction(i, j_min, LIMITER_VENKATAKRISHNAN);
//              // Evaluate SW sub (fine) cell values.
//	      if (!(face == NORTH && sector == WEST && Nghost%2 && j == j_min) &&
//		  !(face == NORTH && sector == EAST && Nghost%2 && (i == i_min || j == j_min)) &&
//		  !(face == SOUTH && sector == EAST && Nghost%2 && i == i_min) &&
//		  !(face == EAST && sector == NORTH && Nghost%2 && (i == i_min || j == j_min)) &&
//		  !(face == EAST && sector == SOUTH && Nghost%2 && i == i_min) &&
//		  !(face == WEST && sector == NORTH && Nghost%2 && j == j_min) &&
//		  !(face == NORTH_EAST && Nghost%2 && (i == i_min || j == j_min)) &&
//		  !(face == NORTH_WEST && Nghost%2 && j == j_min) &&
//		  !(face == SOUTH_EAST && Nghost%2 && i == i_min)) {
//		dX = (Grid.Node[i][j_min].X+
//		      HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
//		      HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
//		      Grid.Cell[i][j_min].Xc)/FOUR -
//		  Grid.Cell[i][j_min].Xc;
//		Wfine = W[i][j_min] +
//		  (phi[i][j_min]^dWdx[i][j_min])*dX.x +
//		  (phi[i][j_min]^dWdy[i][j_min])*dX.y;
//		Ufine = Wfine.U();
//		for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//		  buffer_count = buffer_count + 1;
//		  if (buffer_count >= buffer_size) return(1);
//		  buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//              // Evaluate SE sub (fine) cell values.
//	      if (!(face == NORTH && sector == WEST && Nghost%2 && (i == i_max || j == j_min)) &&
//		  !(face == NORTH && sector == EAST && Nghost%2 && j == j_min) &&
//		  !(face == SOUTH && sector == WEST && Nghost%2 && i == i_max) &&
//		  !(face == EAST && sector == NORTH && Nghost%2 && j == j_min) &&
//		  !(face == WEST && sector == NORTH && Nghost%2 && (i == i_max || j == j_min)) &&
//		  !(face == WEST && sector == SOUTH && Nghost%2 && i == i_max) &&
//		  !(face == NORTH_EAST && Nghost%2 && j == j_min) &&
//		  !(face == NORTH_WEST && Nghost%2 && (i == i_max || j == j_min)) &&
//		  !(face == SOUTH_WEST && Nghost%2 && i == i_max)) {
//              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
//                    Grid.Node[i+1][j_min].X+
//                    Grid.Cell[i][j_min].Xc+
//                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
//                   Grid.Cell[i][j_min].Xc;
//              Wfine = W[i][j_min] +
//                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
//                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//  	         buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//           } /* endfor */
//           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
//              // Evaluate NW sub (fine) cell values.
//	     if (!(face == NORTH && sector == EAST && Nghost%2 && i == i_min) &&
//		 !(face == SOUTH && sector == EAST && Nghost%2 && (i == i_min || j == j_max)) &&
//		 !(face == SOUTH && sector == WEST && Nghost%2 && j == j_max) &&
//		 !(face == EAST && sector == NORTH && Nghost%2 && i == i_min) &&
//		 !(face == EAST && sector == SOUTH && Nghost%2 && (i == i_min || j == j_max)) &&
//		 !(face == WEST && sector == SOUTH && Nghost%2 && j == j_max) &&
//		 !(face == NORTH_EAST && Nghost%2 && i == i_min) &&
//		 !(face == SOUTH_EAST && Nghost%2 && (i == i_min || j == j_max)) &&
//		 !(face == SOUTH_WEST && Nghost%2 && j == j_max)) {
//              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
//                    Grid.Cell[i][j_min].Xc+
//                    Grid.Node[i][j_min+1].X+
//                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
//                   Grid.Cell[i][j_min].Xc;
//              Wfine = W[i][j_min] +
//                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
//                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) { 
//                 buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//              // Evaluate NE sub (fine) cell values.
//	  if (!(face == NORTH && sector == WEST && Nghost%2 && i == i_max) &&
//	      !(face == SOUTH && sector == EAST && Nghost%2 && j == j_max) &&
//	      !(face == SOUTH && sector == WEST && Nghost%2 && (i == i_max || j == j_max)) &&
//	      !(face == EAST && sector == SOUTH && Nghost%2 && j == j_max) &&
//	      !(face == WEST && sector == NORTH && Nghost%2 && i == i_max) &&
//	      !(face == WEST && sector == SOUTH && Nghost%2 && (i == i_max || j == j_max)) &&
//	      !(face == NORTH_WEST && Nghost%2 && i == i_max) &&
//	      !(face == SOUTH_EAST && Nghost%2 && j == j_max) &&
//	      !(face == SOUTH_WEST && Nghost%2 && (i == i_max || j == j_max))) {
//              dX = (Grid.Cell[i][j_min].Xc+
//                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
//                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
//                    Grid.Node[i+1][j_min+1].X)/FOUR -
//                   Grid.Cell[i][j_min].Xc;
//              Wfine = W[i][j_min] +
//                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
//                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//  	         buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//           } /* endfor */
//        } else {
//           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
//              // Perform limited linear least squares reconstruction in cell (i, j_min).
//              SubcellReconstruction(i, j_min, LIMITER_VENKATAKRISHNAN);
//              // Evaluate SE sub (fine) cell values.
//              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
//                    Grid.Node[i+1][j_min].X+
//                    Grid.Cell[i][j_min].Xc+
//                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
//                   Grid.Cell[i][j_min].Xc;
//              Wfine = W[i][j_min] +
//                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
//                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//  	         buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//              // Evaluate SW sub (fine) cell values.
//              dX = (Grid.Node[i][j_min].X+
//                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
//                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
//                    Grid.Cell[i][j_min].Xc)/FOUR -
//                   Grid.Cell[i][j_min].Xc;
//              Wfine = W[i][j_min] +
//                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
//                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//  	         buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//           } /* endfor */
//           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
//              // Evaluate NE sub (fine) cell values.
//              dX = (Grid.Cell[i][j_min].Xc+
//                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
//                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
//                    Grid.Node[i+1][j_min+1].X)/FOUR -
//                   Grid.Cell[i][j_min].Xc;
//              Wfine = W[i][j_min] +
//                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
//                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//  	         buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//              // Evaluate NW sub (fine) cell values.
//              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
//                    Grid.Cell[i][j_min].Xc+
//                    Grid.Node[i][j_min+1].X+
//                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
//                   Grid.Cell[i][j_min].Xc;
//              Wfine = W[i][j_min] +
//                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
//                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) { 
//                 buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//           } /* endfor */
//        } /* endif */
//     } else {
//        if (i_inc > 0) {
//           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
//              // Perform limited linear least squares reconstruction in cell (i, j_min).
//              SubcellReconstruction(i, j_min, LIMITER_VENKATAKRISHNAN);
//              // Evaluate NW sub (fine) cell values.
//              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
//                    Grid.Cell[i][j_min].Xc+
//                    Grid.Node[i][j_min+1].X+
//                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
//                   Grid.Cell[i][j_min].Xc;
//              Wfine = W[i][j_min] +
//                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
//                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) { 
//                 buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//              // Evaluate NE sub (fine) cell values.
//              dX = (Grid.Cell[i][j_min].Xc+
//                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
//                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
//                    Grid.Node[i+1][j_min+1].X)/FOUR -
//                   Grid.Cell[i][j_min].Xc;
//              Wfine = W[i][j_min] +
//                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
//                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//  	         buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//           } /* endfor */
//           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
//              // Evaluate SW sub (fine) cell values.
//              dX = (Grid.Node[i][j_min].X+
//                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
//                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
//                    Grid.Cell[i][j_min].Xc)/FOUR -
//                   Grid.Cell[i][j_min].Xc;
//              Wfine = W[i][j_min] +
//                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
//                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//  	         buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//              // Evaluate SE sub (fine) cell values.
//              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
//                    Grid.Node[i+1][j_min].X+
//                    Grid.Cell[i][j_min].Xc+
//                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
//                   Grid.Cell[i][j_min].Xc;
//              Wfine = W[i][j_min] +
//                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
//                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//  	         buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//           } /* endfor */
//        } else {
//           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
//              // Perform limited linear least squares reconstruction in cell (i, j_min).
//              SubcellReconstruction(i, j_min, LIMITER_VENKATAKRISHNAN);
//              // Evaluate NE sub (fine) cell values.
//              dX = (Grid.Cell[i][j_min].Xc+
//                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
//                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
//                    Grid.Node[i+1][j_min+1].X)/FOUR -
//                   Grid.Cell[i][j_min].Xc;
//              Wfine = W[i][j_min] +
//                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
//                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//  	         buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//              // Evaluate NW sub (fine) cell values.
//              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
//                    Grid.Cell[i][j_min].Xc+
//                    Grid.Node[i][j_min+1].X+
//                    HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
//                   Grid.Cell[i][j_min].Xc;
//              Wfine = W[i][j_min] +
//                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
//                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) { 
//                 buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//           } /* endfor */
//           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
//              // Evaluate SE sub (fine) cell values.
//              dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
//                    Grid.Node[i+1][j_min].X+
//                    Grid.Cell[i][j_min].Xc+
//                    HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
//                   Grid.Cell[i][j_min].Xc;
//              Wfine = W[i][j_min] +
//                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
//                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//  	         buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//              // Evaluate SW sub (fine) cell values.
//              dX = (Grid.Node[i][j_min].X+
//                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
//                    HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
//                    Grid.Cell[i][j_min].Xc)/FOUR -
//                   Grid.Cell[i][j_min].Xc;
//              Wfine = W[i][j_min] +
//                      (phi[i][j_min]^dWdx[i][j_min])*dX.x +
//                      (phi[i][j_min]^dWdy[i][j_min])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//  	         buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//           } /* endfor */
//        } /* endif */
//     } /* endif */
//  } else { // East or west boundary.
//     // Four different orderings to consider depending on the value of i_inc & j_inc.
//     if (j_inc > 0) {
//        if (i_inc > 0) {
//           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
//              // Perform limited linear least squares reconstruction in cell (i_min, j).
//              SubcellReconstruction(i_min, j, LIMITER_VENKATAKRISHNAN);
//              // Evaluate SW sub (fine) cell values.
//              dX = (Grid.Node[i_min][j].X+
//                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
//                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
//                    Grid.Cell[i_min][j].Xc)/FOUR -
//                   Grid.Cell[i_min][j].Xc;
//              Wfine = W[i_min][j] +
//                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
//                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//  	         buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//              // Evaluate SE sub (fine) cell values.
//              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
//                    Grid.Node[i_min+1][j].X+
//                    Grid.Cell[i_min][j].Xc+
//                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
//                   Grid.Cell[i_min][j].Xc;
//              Wfine = W[i_min][j] +
//                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
//                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//  	         buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//              // Evaluate NW sub (fine) cell values.
//              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
//                    Grid.Cell[i_min][j].Xc+
//                    Grid.Node[i_min][j+1].X+
//                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
//                   Grid.Cell[i_min][j].Xc;
//              Wfine = W[i_min][j] +
//                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
//                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) { 
//                 buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//              // Evaluate NE sub (fine) cell values.
//              dX = (Grid.Cell[i_min][j].Xc+
//                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
//                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
//                    Grid.Node[i_min+1][j+1].X)/FOUR -
//                   Grid.Cell[i_min][j].Xc;
//              Wfine = W[i_min][j] +
//                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
//                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//                 buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//           } /* endfor */
//        } else {
//           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
//              // Perform limited linear least squares reconstruction in cell (i_min, j).
//              SubcellReconstruction(i_min, j, LIMITER_VENKATAKRISHNAN);
//              // Evaluate SE sub (fine) cell values.
//              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
//                    Grid.Node[i_min+1][j].X+
//                    Grid.Cell[i_min][j].Xc+
//                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
//                   Grid.Cell[i_min][j].Xc;
//              Wfine = W[i_min][j] +
//                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
//                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//  	         buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//              // Evaluate SW sub (fine) cell values.
//              dX = (Grid.Node[i_min][j].X+
//                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
//                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
//                    Grid.Cell[i_min][j].Xc)/FOUR -
//                   Grid.Cell[i_min][j].Xc;
//              Wfine = W[i_min][j] +
//                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
//                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//  	         buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//              // Evaluate NE sub (fine) cell values.
//             dX = (Grid.Cell[i_min][j].Xc+
//                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
//                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
//                    Grid.Node[i_min+1][j+1].X)/FOUR -
//                   Grid.Cell[i_min][j].Xc;
//              Wfine = W[i_min][j] +
//                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
//                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//                 buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//              // Evaluate NW sub (fine) cell values.
//              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
//                    Grid.Cell[i_min][j].Xc+
//                    Grid.Node[i_min][j+1].X+
//                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
//                   Grid.Cell[i_min][j].Xc;
//              Wfine = W[i_min][j] +
//                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
//                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) { 
//                 buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//           } /* endfor */
//        } /* endif */
//     } else {
//        if (i_inc > 0) {
//           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
//              // Perform limited linear least squares reconstruction in cell (i_min, j).
//              SubcellReconstruction(i_min, j, LIMITER_VENKATAKRISHNAN);
//              // Evaluate NW sub (fine) cell values.
//              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
//                    Grid.Cell[i_min][j].Xc+
//                    Grid.Node[i_min][j+1].X+
//                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
//                   Grid.Cell[i_min][j].Xc;
//              Wfine = W[i_min][j] +
//                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
//                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) { 
//                 buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//              // Evaluate NE sub (fine) cell values.
//              dX = (Grid.Cell[i_min][j].Xc+
//                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
//                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
//                    Grid.Node[i_min+1][j+1].X)/FOUR -
//                   Grid.Cell[i_min][j].Xc;
//              Wfine = W[i_min][j] +
//                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
//                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//                 buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//              // Evaluate SW sub (fine) cell values.
//              dX = (Grid.Node[i_min][j].X+
//                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
//                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
//                    Grid.Cell[i_min][j].Xc)/FOUR -
//                   Grid.Cell[i_min][j].Xc;
//              Wfine = W[i_min][j] +
//                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
//                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//  	         buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//              // Evaluate SE sub (fine) cell values.
//              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
//                    Grid.Node[i_min+1][j].X+
//                    Grid.Cell[i_min][j].Xc+
//                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
//                   Grid.Cell[i_min][j].Xc;
//              Wfine = W[i_min][j] +
//                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
//                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//  	         buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//           } /* endfor */
//        } else {
//           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
//              // Perform limited linear least squares reconstruction in cell (i_min, j).
//              SubcellReconstruction(i_min, j, LIMITER_VENKATAKRISHNAN);
//              // Evaluate NE sub (fine) cell values.
//              dX = (Grid.Cell[i_min][j].Xc+
//                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
//                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
//                    Grid.Node[i_min+1][j+1].X)/FOUR -
//                   Grid.Cell[i_min][j].Xc;
//              Wfine = W[i_min][j] +
//                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
//                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//                 buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//              // Evaluate NW sub (fine) cell values.
//              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
//                    Grid.Cell[i_min][j].Xc+
//                    Grid.Node[i_min][j+1].X+
//                    HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
//                   Grid.Cell[i_min][j].Xc;
//              Wfine = W[i_min][j] +
//                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
//                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) { 
//                 buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//              // Evaluate SE sub (fine) cell values.
//              dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
//                    Grid.Node[i_min+1][j].X+
//                    Grid.Cell[i_min][j].Xc+
//                    HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
//                   Grid.Cell[i_min][j].Xc;
//              Wfine = W[i_min][j] +
//                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
//                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//  	         buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//              // Evaluate SW sub (fine) cell values.
//              dX = (Grid.Node[i_min][j].X+
//                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
//                    HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
//                    Grid.Cell[i_min][j].Xc)/FOUR -
//                   Grid.Cell[i_min][j].Xc;
//              Wfine = W[i_min][j] +
//                      (phi[i_min][j]^dWdx[i_min][j])*dX.x +
//                      (phi[i_min][j]^dWdy[i_min][j])*dX.y;
//              Ufine = Wfine.U();
//              for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
//  	         buffer_count = buffer_count + 1;
//                 if (buffer_count >= buffer_size) return(1);
//                 buffer[buffer_count] = Ufine[k];
//              } /* endfor */
//           } /* endfor */
//        } /* endif */
//     } /* endif */
//  } /* endif */
//  return(0);
}

/**********************************************************************************
 * Gaussian2D_Quad_Block::UnloadReceiveBuffer -- Unloads receive message buffer.  *
 **********************************************************************************/
inline int Gaussian2D_Quad_Block::UnloadReceiveBuffer(double *buffer,
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
  	buffer_count = buffer_count + NUM_VAR_GAUSSIAN2D;
        if (buffer_count >= buffer_size) return(1);
        U[i][j] = Gaussian2D_cState(buffer[buffer_count-7],
                                    buffer[buffer_count-6],
                                    buffer[buffer_count-5],
                                    buffer[buffer_count-4],
                                    buffer[buffer_count-3],
                                    buffer[buffer_count-2],
                                    buffer[buffer_count-1],
                                    buffer[buffer_count]);
        W[i][j] = U[i][j].W();
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * Gaussian2D_Quad_Block::UnloadReceiveBuffer_F2C -- Unloads receive message   *
 *                                                buffer for fine to coarse    *
 *                                                block message passing.       *
 *******************************************************************************/
inline int Gaussian2D_Quad_Block::UnloadReceiveBuffer_F2C(double *buffer,
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
  	buffer_count = buffer_count + NUM_VAR_GAUSSIAN2D;
        if (buffer_count >= buffer_size) return(1);
        U[i][j] = Gaussian2D_cState(buffer[buffer_count-7],
                                    buffer[buffer_count-6],
                                    buffer[buffer_count-5],
                                    buffer[buffer_count-4],
                                    buffer[buffer_count-3],
                                    buffer[buffer_count-2],
                                    buffer[buffer_count-1],
                                    buffer[buffer_count]);
        W[i][j] = U[i][j].W();
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * Gaussian2D_Quad_Block::UnloadReceiveBuffer_C2F -- Unloads receive message   *
 *                                                buffer for coarse to fine    *
 *                                                block message passing.       *
 *******************************************************************************/
inline int Gaussian2D_Quad_Block::UnloadReceiveBuffer_C2F(double *buffer,
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
  	buffer_count = buffer_count + NUM_VAR_GAUSSIAN2D;
        if (buffer_count >= buffer_size) return(1);
        U[i][j] = Gaussian2D_cState(buffer[buffer_count-7],
                                    buffer[buffer_count-6],
                                    buffer[buffer_count-5],
                                    buffer[buffer_count-4],
                                    buffer[buffer_count-3],
                                    buffer[buffer_count-2],
                                    buffer[buffer_count-1],
                                    buffer[buffer_count]);
        W[i][j] = U[i][j].W();
     } /* endfor */
  } /* endfor */
  return(0);
}

/**************************************************************************
 * Gaussian2D_Quad_Block::SubcellReconstruction --                        *
 *               Performs the subcell reconstruction of solution state    *
 *               within a given cell (i,j) of the computational mesh for  *
 *               the specified quadrilateral solution block.              *
 **************************************************************************/
inline void Gaussian2D_Quad_Block::SubcellReconstruction(const int i, 
                                                         const int j,
                                                         const int Limiter) {

  int n, n2, n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4], phi_n;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  Gaussian2D_pState DU, DUDx_ave, DUDy_ave;

  /* Carry out the limited solution reconstruction in
     each cell of the computational mesh. */

  //If on outside ring of ghost cells, no reconstruction

  if (i == ICl-2 || i == ICu+2 ||
      j == JCl-2 || j == JCu+2 ) {
    n_pts = 0;

    //If on inside ring of ghost cells, reconstruction may be required

    /////////////left edge//////////////////

  } else if (i == ICl-1) {
    //bottom or top left corner (never used to evaluate fluxes, therefore
    //no reconstruction needed
    if(j == JCl-1 || j == JCu+1){
      n_pts = 0;
      //left edge, may need reconstruction
      //(only for BC_NONE and BC_PERIODIC)???
    } else if(Grid.BCtypeW[j] == BC_FIXED ||
	      Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
	      Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
	      Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
	      Grid.BCtypeW[j] == BC_REFLECTION ||
	      Grid.BCtypeW[j] == BC_ADIABATIC_WALL) {
      n_pts = 0;
    } else { //BC's are Periodic or "none"
      //Bottom left
      if(j == JCl &&
	 (Grid.BCtypeS[i] == BC_FIXED ||
	  Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
	  Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
	  Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
	  Grid.BCtypeS[i] == BC_REFLECTION ||
	  Grid.BCtypeS[i] == BC_ADIABATIC_WALL)) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i;   j_index[2] = j+1;
	i_index[3] = i+1; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j;
	//Top Left
      } else if (j == JCu &&
		 (Grid.BCtypeN[i] == BC_FIXED ||
		  Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
		  Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
		  Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
		  Grid.BCtypeN[i] == BC_REFLECTION ||
		  Grid.BCtypeN[i] == BC_ADIABATIC_WALL)) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j;
	i_index[1] = i-1; j_index[1] = j-1;
	i_index[2] = i;   j_index[2] = j-1;
	i_index[3] = i+1; j_index[3] = j-1;
	i_index[4] = i+1; j_index[4] = j;
	//middle left
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i-1; j_index[1] = j;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i;   j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
	i_index[5] = i+1; j_index[5] = j;
	i_index[6] = i+1; j_index[6] = j-1;
	i_index[7] = i;   j_index[7] = j-1;
      } //endif
    } //endif

      /////////////Right edge//////////////////

  } else if (i == ICu+-1) {
    //bottom or top left corner (never used to evaluate fluxes, therefore
    //no reconstruction needed
    if(j == JCl-1 || j == JCu+1){
      n_pts = 0;
      //right edge, may need reconstruction
      //(only for BC_NONE and BC_PERIODIC)???
    } else if(Grid.BCtypeE[j] == BC_FIXED ||
	      Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
	      Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
	      Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
	      Grid.BCtypeE[j] == BC_REFLECTION ||
	      Grid.BCtypeE[j] == BC_ADIABATIC_WALL) {
      n_pts = 0;
    } else { //BC's are Periodic or "none"
      //Bottom right
      if(j == JCl &&
	 (Grid.BCtypeS[i] == BC_FIXED ||
	  Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
	  Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
	  Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
	  Grid.BCtypeS[i] == BC_REFLECTION ||
	  Grid.BCtypeS[i] == BC_ADIABATIC_WALL)) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i;   j_index[2] = j+1;
	i_index[3] = i+1; j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j;
	//Top right
      } else if (j == JCu &&
		 (Grid.BCtypeN[i] == BC_FIXED ||
		  Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
		  Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
		  Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
		  Grid.BCtypeN[i] == BC_REFLECTION ||
		  Grid.BCtypeN[i] == BC_ADIABATIC_WALL)) {
	n_pts = 5;
	i_index[0] = i-1; j_index[0] = j;
	i_index[1] = i-1; j_index[1] = j-1;
	i_index[2] = i;   j_index[2] = j-1;
	i_index[3] = i+1; j_index[3] = j-1;
	i_index[4] = i+1; j_index[4] = j;
	//middle right or top or bottom with none or priodic bc's
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i-1; j_index[1] = j;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i;   j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
	i_index[5] = i+1; j_index[5] = j;
	i_index[6] = i+1; j_index[6] = j-1;
	i_index[7] = i;   j_index[7] = j-1;
      } //endif
    } //endif

      /////////////bottom edge/////////////
      // (corners already taken care of) //

  } else if (j == JCl-1) {
    if(Grid.BCtypeS[i] == BC_FIXED ||
       Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
       Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
       Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
       Grid.BCtypeS[i] == BC_REFLECTION ||
       Grid.BCtypeS[i] == BC_ADIABATIC_WALL) {
      n_pts = 0;
    } else { //BC's are Periodic or "none"
      //Bottom left
      if(i == ICl &&
	 (Grid.BCtypeW[j] == BC_FIXED ||
	  Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
	  Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
	  Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
	  Grid.BCtypeW[j] == BC_REFLECTION ||
	  Grid.BCtypeW[j] == BC_ADIABATIC_WALL)) {
	n_pts = 5;
	i_index[0] = i;   j_index[0] = j+1;
	i_index[1] = i+1; j_index[1] = j+1;
	i_index[2] = i+1; j_index[2] = j;
	i_index[3] = i+1; j_index[3] = j-1;
	i_index[4] = i;   j_index[4] = j-1;
	//Bottom right
      } else if (i == ICu &&
		 (Grid.BCtypeE[j] == BC_FIXED ||
		  Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
		  Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
		  Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
		  Grid.BCtypeE[j] == BC_REFLECTION ||
		  Grid.BCtypeE[j] == BC_ADIABATIC_WALL)) {
	n_pts = 5;
	i_index[0] = i;   j_index[0] = j+1;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i-1; j_index[2] = j;
	i_index[3] = i-1; j_index[3] = j-1;
	i_index[4] = i;   j_index[4] = j-1;
	//middle bottom or left or right with none or priodic bc's
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i-1; j_index[1] = j;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i;   j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
	i_index[5] = i+1; j_index[5] = j;
	i_index[6] = i+1; j_index[6] = j-1;
	i_index[7] = i;   j_index[7] = j-1;
      } //endif
    } //endif

      /////////////top edge/////////////
      // (corners already taken care of) //

  } else if (j == JCu+1) {
    if(Grid.BCtypeN[i] == BC_FIXED ||
       Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
       Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
       Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
       Grid.BCtypeN[i] == BC_REFLECTION ||
       Grid.BCtypeN[i] == BC_ADIABATIC_WALL) {
      n_pts = 0;
    } else { //BC's are Periodic or "none"
      //Top left
      if(i == ICl &&
	 (Grid.BCtypeW[j] == BC_FIXED ||
	  Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
	  Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
	  Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
	  Grid.BCtypeW[j] == BC_REFLECTION ||
	  Grid.BCtypeW[j] == BC_ADIABATIC_WALL)) {
	n_pts = 5;
	i_index[0] = i;   j_index[0] = j+1;
	i_index[1] = i+1; j_index[1] = j+1;
	i_index[2] = i+1; j_index[2] = j;
	i_index[3] = i+1; j_index[3] = j-1;
	i_index[4] = i;   j_index[4] = j-1;
	//Top right
      } else if (i == ICu &&
		 (Grid.BCtypeE[j] == BC_FIXED ||
		  Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
		  Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
		  Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
		  Grid.BCtypeE[j] == BC_REFLECTION ||
		  Grid.BCtypeE[j] == BC_ADIABATIC_WALL)) {
	n_pts = 5;
	i_index[0] = i;   j_index[0] = j+1;
	i_index[1] = i-1; j_index[1] = j+1;
	i_index[2] = i-1; j_index[2] = j;
	i_index[3] = i-1; j_index[3] = j-1;
	i_index[4] = i;   j_index[4] = j-1;
	//top bottom or left or right with none or priodic bc's
      } else {
	n_pts = 8;
	i_index[0] = i-1; j_index[0] = j-1;
	i_index[1] = i-1; j_index[1] = j;
	i_index[2] = i-1; j_index[2] = j+1;
	i_index[3] = i;   j_index[3] = j+1;
	i_index[4] = i+1; j_index[4] = j+1;
	i_index[5] = i+1; j_index[5] = j;
	i_index[6] = i+1; j_index[6] = j-1;
	i_index[7] = i;   j_index[7] = j-1;
      } //endif
    } //endif
    
      /////////First Cell Inside Computational Domain////////

      //bottom left
  } else if (i == ICl && j == JCl ) {
    if((Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
	Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
	Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
	Grid.BCtypeW[j] == BC_REFLECTION ||
	Grid.BCtypeW[j] == BC_ADIABATIC_WALL)
       &&
       (Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
	Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
	Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
	Grid.BCtypeS[i] == BC_REFLECTION ||
	Grid.BCtypeS[i] == BC_ADIABATIC_WALL)) {
      n_pts = 3;
      i_index[0] = i;   j_index[0] = j+1;
      i_index[1] = i+1; j_index[1] = j+1;
      i_index[2] = i+1; j_index[2] = j;
    } else if (Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
	       Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
	       Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
	       Grid.BCtypeS[i] == BC_REFLECTION ||
	       Grid.BCtypeS[i] == BC_ADIABATIC_WALL) {
      n_pts = 5;
      i_index[0] = i-1; j_index[0] = j;
      i_index[1] = i-1; j_index[1] = j+1;
      i_index[2] = i;   j_index[2] = j+1;
      i_index[3] = i+1; j_index[3] = j+1;
      i_index[4] = i+1; j_index[4] = j;
    } else if (Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
	       Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
	       Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
	       Grid.BCtypeW[j] == BC_REFLECTION ||
	       Grid.BCtypeW[j] == BC_ADIABATIC_WALL) {
      n_pts = 5;
      i_index[0] = i;   j_index[0] = j+1;
      i_index[1] = i+1; j_index[1] = j+1;
      i_index[2] = i+1; j_index[2] = j;
      i_index[3] = i+1; j_index[3] = j-1;
      i_index[4] = i;   j_index[4] = j-1;
    } else {
      n_pts = 8;
      i_index[0] = i-1; j_index[0] = j-1;
      i_index[1] = i-1; j_index[1] = j;
      i_index[2] = i-1; j_index[2] = j+1;
      i_index[3] = i;   j_index[3] = j+1;
      i_index[4] = i+1; j_index[4] = j+1;
      i_index[5] = i+1; j_index[5] = j;
      i_index[6] = i+1; j_index[6] = j-1;
      i_index[7] = i;   j_index[7] = j-1;
    } //endif

      //top left
  } else if (i == ICl && j == JCu ) {
    if((Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
	Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
	Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
	Grid.BCtypeW[j] == BC_REFLECTION ||
	Grid.BCtypeW[j] == BC_ADIABATIC_WALL)
       &&
       (Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
	Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
	Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
	Grid.BCtypeN[i] == BC_REFLECTION ||
	Grid.BCtypeN[i] == BC_ADIABATIC_WALL)) {
      n_pts = 3;
      i_index[0] = i;   j_index[0] = j-1;
      i_index[1] = i+1; j_index[1] = j-1;
      i_index[2] = i+1; j_index[2] = j;
    } else if (Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
	       Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
	       Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
	       Grid.BCtypeN[i] == BC_REFLECTION ||
	       Grid.BCtypeN[i] == BC_ADIABATIC_WALL) {
      n_pts = 5;
      i_index[0] = i-1; j_index[0] = j;
      i_index[1] = i-1; j_index[1] = j-1;
      i_index[2] = i;   j_index[2] = j-1;
      i_index[3] = i+1; j_index[3] = j-1;
      i_index[4] = i+1; j_index[4] = j;
    } else if (Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
	       Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
	       Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
	       Grid.BCtypeW[j] == BC_REFLECTION ||
	       Grid.BCtypeW[j] == BC_ADIABATIC_WALL) {
      n_pts = 5;
      i_index[0] = i;   j_index[0] = j+1;
      i_index[1] = i+1; j_index[1] = j+1;
      i_index[2] = i+1; j_index[2] = j;
      i_index[3] = i+1; j_index[3] = j-1;
      i_index[4] = i;   j_index[4] = j-1;
    } else {
      n_pts = 8;
      i_index[0] = i-1; j_index[0] = j-1;
      i_index[1] = i-1; j_index[1] = j;
      i_index[2] = i-1; j_index[2] = j+1;
      i_index[3] = i;   j_index[3] = j+1;
      i_index[4] = i+1; j_index[4] = j+1;
      i_index[5] = i+1; j_index[5] = j;
      i_index[6] = i+1; j_index[6] = j-1;
      i_index[7] = i;   j_index[7] = j-1;
    } //endif

      //top right
  } else if (i == ICu && j == JCu ) {
    if((Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
	Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
	Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
	Grid.BCtypeE[j] == BC_REFLECTION ||
	Grid.BCtypeE[j] == BC_ADIABATIC_WALL)
       &&
       (Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
	Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
	Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
	Grid.BCtypeN[i] == BC_REFLECTION ||
	Grid.BCtypeN[i] == BC_ADIABATIC_WALL)) {
      n_pts = 3;
      i_index[0] = i-1; j_index[0] = j;
      i_index[1] = i-1; j_index[1] = j-1;
      i_index[2] = i;   j_index[2] = j-1;
    } else if (Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
	       Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
	       Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
	       Grid.BCtypeN[i] == BC_REFLECTION ||
	       Grid.BCtypeN[i] == BC_ADIABATIC_WALL) {
      n_pts = 5;
      i_index[0] = i-1; j_index[0] = j;
      i_index[1] = i-1; j_index[1] = j-1;
      i_index[2] = i;   j_index[2] = j-1;
      i_index[3] = i+1; j_index[3] = j-1;
      i_index[4] = i+1; j_index[4] = j;
    } else if (Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
	       Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
	       Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
	       Grid.BCtypeE[j] == BC_REFLECTION ||
	       Grid.BCtypeE[j] == BC_ADIABATIC_WALL) {
      n_pts = 5;
      i_index[0] = i;   j_index[0] = j+1;
      i_index[1] = i-1; j_index[1] = j+1;
      i_index[2] = i-1; j_index[2] = j;
      i_index[3] = i-1; j_index[3] = j-1;
      i_index[4] = i;   j_index[4] = j-1;
    } else {
      n_pts = 8;
      i_index[0] = i-1; j_index[0] = j-1;
      i_index[1] = i-1; j_index[1] = j;
      i_index[2] = i-1; j_index[2] = j+1;
      i_index[3] = i;   j_index[3] = j+1;
      i_index[4] = i+1; j_index[4] = j+1;
      i_index[5] = i+1; j_index[5] = j;
      i_index[6] = i+1; j_index[6] = j-1;
      i_index[7] = i;   j_index[7] = j-1;
    } //endif

      //Bottom right
  } else if (i == ICu && j == JCl ) {
    if((Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
	Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
	Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
	Grid.BCtypeE[j] == BC_REFLECTION ||
	Grid.BCtypeE[j] == BC_ADIABATIC_WALL)
       &&
       (Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
	Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
	Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
	Grid.BCtypeS[i] == BC_REFLECTION ||
	Grid.BCtypeS[i] == BC_ADIABATIC_WALL)) {
      n_pts = 3;
      i_index[0] = i;   j_index[0] = j+1;
      i_index[1] = i-1; j_index[1] = j+1;
      i_index[2] = i-1; j_index[2] = j;
    } else if (Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
	       Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
	       Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
	       Grid.BCtypeS[i] == BC_REFLECTION ||
	       Grid.BCtypeS[i] == BC_ADIABATIC_WALL) {
      n_pts = 5;
      i_index[0] = i-1; j_index[0] = j;
      i_index[1] = i-1; j_index[1] = j+1;
      i_index[2] = i;   j_index[2] = j+1;
      i_index[3] = i+1; j_index[3] = j+1;
      i_index[4] = i+1; j_index[4] = j;
    } else if (Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
	       Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
	       Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
	       Grid.BCtypeE[j] == BC_REFLECTION ||
	       Grid.BCtypeE[j] == BC_ADIABATIC_WALL) {
      n_pts = 5;
      i_index[0] = i;   j_index[0] = j+1;
      i_index[1] = i-1; j_index[1] = j+1;
      i_index[2] = i-1; j_index[2] = j;
      i_index[3] = i-1; j_index[3] = j-1;
      i_index[4] = i;   j_index[4] = j-1;
    } else {
      n_pts = 8;
      i_index[0] = i-1; j_index[0] = j-1;
      i_index[1] = i-1; j_index[1] = j;
      i_index[2] = i-1; j_index[2] = j+1;
      i_index[3] = i;   j_index[3] = j+1;
      i_index[4] = i+1; j_index[4] = j+1;
      i_index[5] = i+1; j_index[5] = j;
      i_index[6] = i+1; j_index[6] = j-1;
      i_index[7] = i;   j_index[7] = j-1;
    } //endif

      //left
  } else if (i == ICl) {

    if(Grid.BCtypeW[j] == BC_CHARACTERISTIC ||
       Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
       Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
       Grid.BCtypeW[j] == BC_REFLECTION ||
       Grid.BCtypeW[j] == BC_ADIABATIC_WALL) {
      n_pts = 5;
      i_index[0] = i;   j_index[0] = j+1;
      i_index[1] = i+1; j_index[1] = j+1;
      i_index[2] = i+1; j_index[2] = j;
      i_index[3] = i+1; j_index[3] = j-1;
      i_index[4] = i;   j_index[4] = j-1;
    } else {
      n_pts = 8;
      i_index[0] = i-1; j_index[0] = j-1;
      i_index[1] = i-1; j_index[1] = j;
      i_index[2] = i-1; j_index[2] = j+1;
      i_index[3] = i;   j_index[3] = j+1;
      i_index[4] = i+1; j_index[4] = j+1;
      i_index[5] = i+1; j_index[5] = j;
      i_index[6] = i+1; j_index[6] = j-1;
      i_index[7] = i;   j_index[7] = j-1;
    } //endif

      //top
  } else if (j == JCu) {

    if(Grid.BCtypeN[i] == BC_CHARACTERISTIC ||
       Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
       Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
       Grid.BCtypeN[i] == BC_REFLECTION ||
       Grid.BCtypeN[i] == BC_ADIABATIC_WALL) {
      n_pts = 5;
      i_index[0] = i-1; j_index[0] = j;
      i_index[1] = i-1; j_index[1] = j-1;
      i_index[2] = i;   j_index[2] = j-1;
      i_index[3] = i+1; j_index[3] = j-1;
      i_index[4] = i+1; j_index[4] = j;
    } else {
      n_pts = 8;
      i_index[0] = i-1; j_index[0] = j-1;
      i_index[1] = i-1; j_index[1] = j;
      i_index[2] = i-1; j_index[2] = j+1;
      i_index[3] = i;   j_index[3] = j+1;
      i_index[4] = i+1; j_index[4] = j+1;
      i_index[5] = i+1; j_index[5] = j;
      i_index[6] = i+1; j_index[6] = j-1;
      i_index[7] = i;   j_index[7] = j-1;
    } //endif

      //Right
  } else if (i == ICu) {

    if(Grid.BCtypeE[j] == BC_CHARACTERISTIC ||
       Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
       Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
       Grid.BCtypeE[j] == BC_REFLECTION ||
       Grid.BCtypeE[j] == BC_ADIABATIC_WALL) {
      n_pts = 5;
      i_index[0] = i;   j_index[0] = j+1;
      i_index[1] = i-1; j_index[1] = j+1;
      i_index[2] = i-1; j_index[2] = j;
      i_index[3] = i-1; j_index[3] = j-1;
      i_index[4] = i;   j_index[4] = j-1;
    } else {
      n_pts = 8;
      i_index[0] = i-1; j_index[0] = j-1;
      i_index[1] = i-1; j_index[1] = j;
      i_index[2] = i-1; j_index[2] = j+1;
      i_index[3] = i;   j_index[3] = j+1;
      i_index[4] = i+1; j_index[4] = j+1;
      i_index[5] = i+1; j_index[5] = j;
      i_index[6] = i+1; j_index[6] = j-1;
      i_index[7] = i;   j_index[7] = j-1;
    } //endif

      //Bottom
  } else if (j == JCl) {

    if(Grid.BCtypeS[i] == BC_CHARACTERISTIC ||
       Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
       Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
       Grid.BCtypeS[i] == BC_REFLECTION ||
       Grid.BCtypeS[i] == BC_ADIABATIC_WALL) {
      n_pts = 5;
      i_index[0] = i-1; j_index[0] = j;
      i_index[1] = i-1; j_index[1] = j+1;
      i_index[2] = i;   j_index[2] = j+1;
      i_index[3] = i+1; j_index[3] = j+1;
      i_index[4] = i+1; j_index[4] = j;
    } else {
      n_pts = 8;
      i_index[0] = i-1; j_index[0] = j-1;
      i_index[1] = i-1; j_index[1] = j;
      i_index[2] = i-1; j_index[2] = j+1;
      i_index[3] = i;   j_index[3] = j+1;
      i_index[4] = i+1; j_index[4] = j+1;
      i_index[5] = i+1; j_index[5] = j;
      i_index[6] = i+1; j_index[6] = j-1;
      i_index[7] = i;   j_index[7] = j-1;
    } //endif
    
      ///////inside computational domain//////////

  } else {
    n_pts = 8;
    i_index[0] = i-1; j_index[0] = j-1;
    i_index[1] = i-1; j_index[1] = j;
    i_index[2] = i-1; j_index[2] = j+1;
    i_index[3] = i;   j_index[3] = j+1;
    i_index[4] = i+1; j_index[4] = j+1;
    i_index[5] = i+1; j_index[5] = j;
    i_index[6] = i+1; j_index[6] = j-1;
    i_index[7] = i;   j_index[7] = j-1;
  } //endif

  if (n_pts > 0) {
      DUDx_ave = Gaussian2D_W_VACUUM;
      DUDy_ave = Gaussian2D_W_VACUUM;
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
  
      // Calculate slope limiters. 
      if (!Freeze_Limiter) {
         for ( n = 1 ; n <= NUM_VAR_GAUSSIAN2D ; ++n ) {
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
      } /* endif */
  } else {
      dWdx[i][j] = Gaussian2D_W_VACUUM;
      dWdy[i][j] = Gaussian2D_W_VACUUM; 
      phi[i][j]  = Gaussian2D_W_VACUUM;
  } /* endif */
    
}

/**********************************************************************************
 * Gaussian2D_Quad_Block::LoadSendBuffer_Flux_F2C -- Loads send message buffer for*
 *                                                   fine to coarse block message *
 *                                                   passing of conservative      *
 *                                                   solution fluxes.             *
 **********************************************************************************/
inline int Gaussian2D_Quad_Block::LoadSendBuffer_Flux_F2C(double *buffer,
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
        for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxS[i  ][k]+
                                   FluxS[i+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxN[i  ][k]+
                                   FluxN[i+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
     for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
        for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxW[j][k]+
                                   FluxW[j+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
     for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
        for ( k = 1 ; k <= NUM_VAR_GAUSSIAN2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxE[j][k]+
                                   FluxE[j+1][k]);
        } /* endfor */
     } /* endfor */
  } /* endif */
  return(0);
}

/**********************************************************************************
 * Gaussian2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C -- Unloads receive message *
 *                                                   buffer for fine to coarse    *
 *                                                   block message passing of     *
 *                                                   conservative solution fluxes.*
 **********************************************************************************/
inline int Gaussian2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C(double *buffer,
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
  	buffer_count = buffer_count + NUM_VAR_GAUSSIAN2D;
        if (buffer_count >= buffer_size) return(1);
        FluxS[i] = -Gaussian2D_cState(buffer[buffer_count-7],
                                      buffer[buffer_count-6],
                                      buffer[buffer_count-5],
                                      buffer[buffer_count-4],
                                      buffer[buffer_count-3],
                                      buffer[buffer_count-2],
                                      buffer[buffer_count-1],
                                      buffer[buffer_count])
                   -FluxS[i];
     } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
  	buffer_count = buffer_count + NUM_VAR_GAUSSIAN2D;
        if (buffer_count >= buffer_size) return(1);
        FluxN[i] = -Gaussian2D_cState(buffer[buffer_count-7],
                                      buffer[buffer_count-6],
                                      buffer[buffer_count-5],
                                      buffer[buffer_count-4],
                                      buffer[buffer_count-3],
                                      buffer[buffer_count-2],
                                      buffer[buffer_count-1],
                                      buffer[buffer_count])
                   -FluxN[i];
     } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
     for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
  	buffer_count = buffer_count + NUM_VAR_GAUSSIAN2D;
        if (buffer_count >= buffer_size) return(1);
        FluxW[j] = -Gaussian2D_cState(buffer[buffer_count-7],
                                      buffer[buffer_count-6],
                                      buffer[buffer_count-5],
                                      buffer[buffer_count-4],
                                      buffer[buffer_count-3],
                                      buffer[buffer_count-2],
                                      buffer[buffer_count-1],
                                      buffer[buffer_count])
                   -FluxW[j];
     } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
     for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
  	buffer_count = buffer_count + NUM_VAR_GAUSSIAN2D;
        if (buffer_count >= buffer_size) return(1);
        FluxE[j] = -Gaussian2D_cState(buffer[buffer_count-7],
                                      buffer[buffer_count-6],
                                      buffer[buffer_count-5],
                                      buffer[buffer_count-4],
                                      buffer[buffer_count-3],
                                      buffer[buffer_count-2],
                                      buffer[buffer_count-1],
                                      buffer[buffer_count])
                   -FluxE[j];
     } /* endfor */
  } /* endif */
  return(0);
}

/**************************************************************************
 * Gaussian2D_Quad_Block -- Single Block External Subroutines.            *
 **************************************************************************/

extern void Write_Solution_Block(Gaussian2D_Quad_Block &SolnBlk,
	                         ostream &Out_File);

extern void Read_Solution_Block(Gaussian2D_Quad_Block &SolnBlk,
	                        istream &In_File);

extern void Broadcast_Solution_Block(Gaussian2D_Quad_Block &SolnBlk);

#ifdef _MPI_VERSION
extern void Broadcast_Solution_Block(Gaussian2D_Quad_Block &SolnBlk,
                                     MPI::Intracomm &Communicator, 
                                     const int Source_CPU);
#endif

extern void Copy_Solution_Block(Gaussian2D_Quad_Block &SolnBlk1,
		                Gaussian2D_Quad_Block &SolnBlk2);

extern int Prolong_Solution_Block(Gaussian2D_Quad_Block &SolnBlk_Fine,
		                   Gaussian2D_Quad_Block &SolnBlk_Original,
                                   const int Sector);

extern int Restrict_Solution_Block(Gaussian2D_Quad_Block &SolnBlk_Coarse,
		                    Gaussian2D_Quad_Block &SolnBlk_Original_SW,
                                    Gaussian2D_Quad_Block &SolnBlk_Original_SE,
                                    Gaussian2D_Quad_Block &SolnBlk_Original_NW,
                                    Gaussian2D_Quad_Block &SolnBlk_Original_NE);


extern void Output_Tecplot(Gaussian2D_Quad_Block &SolnBlk,
		           const int Number_of_Time_Steps,
                           const double &Time,
                           const int Block_Number,
                           const int Output_Title,
	                   ostream &Out_File);

extern void Output_Tecplot(Gaussian2D_Quad_Block &SolnBlk,
			   Gaussian2D_Input_Parameters &IP,
		           const int Number_of_Time_Steps,
                           const double &Time,
                           const int Block_Number,
                           const int Output_Title,
	                   ostream &Out_File);

extern void Output_Cells_Tecplot(Gaussian2D_Quad_Block &SolnBlk,
				 Gaussian2D_Input_Parameters &IP,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_Gradients_Tecplot(Gaussian2D_Quad_Block &SolnBlk,
				     const int Number_of_Time_Steps,
				     const double &Time,
				     const int Block_Number,
				     const int Output_Title,
				     ostream &Out_File);

extern void Output_Flat_Plate(Gaussian2D_Quad_Block &SolnBlk,
			      const int Number_of_Time_Steps,
			      const double &Time,
			      const int Block_Number,
			      const int Output_Title,
			      ostream &Out_File,
			      const Gaussian2D_pState &Winf);

extern void Output_Cylinder_Free_Molecular(Gaussian2D_Quad_Block &SolnBlk,
					   const int Number_of_Time_Steps,
					   const double &Time,
					   const int Block_Number,
					   const int Output_Title,
					   ostream &Out_File,
					   const Gaussian2D_pState &Winf,
					   Vector2D *nodes,
					   const int number_of_nodes);

extern void Append_nodes_to_send_buffer(Gaussian2D_Quad_Block &SolnBlk,
					double *buffer1, double *buffer2, 
					int &count);

extern void Output_Drag(Gaussian2D_Quad_Block &SolnBlk,
			double &drag, double &lift);

extern void ICs(Gaussian2D_Quad_Block &SolnBlk,
 	        Gaussian2D_Input_Parameters &Input_Parameters,
                Gaussian2D_pState *Wo);

extern void BCs(Gaussian2D_Quad_Block &SolnBlk,
		Gaussian2D_Input_Parameters &IP);

extern double CFL(Gaussian2D_Quad_Block &SolnBlk);

extern void Set_Global_TimeStep(Gaussian2D_Quad_Block &SolnBlk, 
                                const double &Dt_min);

extern double L1_Norm_Residual(Gaussian2D_Quad_Block &SolnBlk);

extern double L2_Norm_Residual(Gaussian2D_Quad_Block &SolnBlk);

extern double Max_Norm_Residual(Gaussian2D_Quad_Block &SolnBlk);

extern void Linear_Reconstruction_GreenGauss(Gaussian2D_Quad_Block &SolnBlk,
                                             const int i,
                                             const int j,
					     const int Limiter);

extern void Linear_Reconstruction_GreenGauss(Gaussian2D_Quad_Block &SolnBlk,
					     const int Limiter);

extern void Linear_Reconstruction_LeastSquares(Gaussian2D_Quad_Block &SolnBlk,
                                               const int i,
                                               const int j,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares(Gaussian2D_Quad_Block &SolnBlk,
					       const int Limiter);

extern void Calculate_Refinement_Criteria(double *refinement_criteria,
					  Gaussian2D_Input_Parameters &IP,
                                          int &number_refinement_criteria,
                                          Gaussian2D_Quad_Block &SolnBlk);

extern void Fix_Refined_Block_Boundaries(Gaussian2D_Quad_Block &SolnBlk,
                                         const int Fix_North_Boundary,
                                         const int Fix_South_Boundary,
                                         const int Fix_East_Boundary,
                                         const int Fix_West_Boundary);

extern void Unfix_Refined_Block_Boundaries(Gaussian2D_Quad_Block &SolnBlk);

extern void Apply_Boundary_Flux_Corrections(Gaussian2D_Quad_Block &SolnBlk,
                                            const int Number_Neighbours_North_Boundary,
                                            const int Number_Neighbours_South_Boundary,
                                            const int Number_Neighbours_East_Boundary,
                                            const int Number_Neighbours_West_Boundary);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Gaussian2D_Quad_Block &SolnBlk,
                                                                const int i_stage,
                                                                const int n_stage,
                                                                const double &CFL_Number,
                                                                const int Time_Integration_Type,
                                                                const int Reconstruction_Type,
                                                                const int Limiter_Type,
                                                                const int Flux_Function_Type,
                                                                const int Number_Neighbours_North_Boundary,
                                                                const int Number_Neighbours_South_Boundary,
                                                                const int Number_Neighbours_East_Boundary,
                                                                const int Number_Neighbours_West_Boundary);

extern void dUdt_Residual_Evaluation(Gaussian2D_Quad_Block &SolnBlk,
                                     Gaussian2D_Input_Parameters &Input_Parameters);

extern int dUdt_Multistage_Explicit(Gaussian2D_Quad_Block &SolnBlk,
   	                            const int i_stage,
                                    Gaussian2D_Input_Parameters &Input_Parameters);

extern int Update_Solution_Multistage_Explicit(Gaussian2D_Quad_Block &SolnBlk,
   	                                       const int i_stage,
                                               Gaussian2D_Input_Parameters &Input_Parameters);

extern void Ramp_up_Reference_Mach_Number(Gaussian2D_Quad_Block &SolnBlk,
					  Gaussian2D_Input_Parameters &Input_Parameters,
					  const int Number_of_Time_Steps);

extern int Count_south_adiabatic_wall_cells(Gaussian2D_Quad_Block &SolnBlk);

/**************************************************************************
 * Gaussian2D_Quad_Block -- Multiple Block External Subroutines.          *
 **************************************************************************/

extern Gaussian2D_Quad_Block* Allocate(Gaussian2D_Quad_Block *Soln_ptr,
                                       Gaussian2D_Input_Parameters &Input_Parameters);

extern Gaussian2D_Quad_Block* Deallocate(Gaussian2D_Quad_Block *Soln_ptr,
                                         Gaussian2D_Input_Parameters &Input_Parameters);

extern Gaussian2D_Quad_Block* CreateInitialSolutionBlocks(Grid2D_Quad_Block **InitMeshBlk,
							  Gaussian2D_Quad_Block *Soln_ptr,
							  Gaussian2D_Input_Parameters &Input_Parameters,
							  QuadTreeBlock_DataStructure &QuadTree,
							  AdaptiveBlockResourceList &GlobalSolnBlockList,
							  AdaptiveBlock2D_List &LocalSolnBlockList);

extern void ICs(Gaussian2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
                Gaussian2D_Input_Parameters &Input_Parameters);

extern int Read_Restart_Solution(Gaussian2D_Quad_Block *Soln_ptr,
                                 AdaptiveBlock2D_List &Soln_Block_List,
                                 Gaussian2D_Input_Parameters &Input_Parameters,
		                 int &Number_of_Time_Steps,
                                 double &Time,
                                 CPUTime &CPU_Time);

extern int Write_Restart_Solution(Gaussian2D_Quad_Block *Soln_ptr,
                                  AdaptiveBlock2D_List &Soln_Block_List,
                                  Gaussian2D_Input_Parameters &Input_Parameters,
		                  const int Number_of_Time_Steps,
                                  const double &Time,
                                  const CPUTime &CPU_Time);

extern int Output_Tecplot(Gaussian2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          Gaussian2D_Input_Parameters &Input_Parameters,
		          const int Number_of_Time_Steps,
                          const double &Time);

extern int Output_Cells_Tecplot(Gaussian2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                Gaussian2D_Input_Parameters &Input_Parameters,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Gradients_Tecplot(Gaussian2D_Quad_Block *Soln_ptr,
				    AdaptiveBlock2D_List &Soln_Block_List,
				    Gaussian2D_Input_Parameters &Input_Parameters,
				    const int Number_of_Time_Steps,
				    const double &Time);

extern int Output_Mesh_Tecplot(Gaussian2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Gaussian2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Mesh_Gnuplot(Gaussian2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Gaussian2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Flat_Plate(Gaussian2D_Quad_Block *Soln_ptr,
			     AdaptiveBlock2D_List &Soln_Block_List,
			     Gaussian2D_Input_Parameters &Input_Parameters,
			     const int Number_of_Time_Steps,
			     const double &Time);

extern int Output_Cylinder_Free_Molecular(Gaussian2D_Quad_Block *Soln_ptr,
					  AdaptiveBlock2D_List &Soln_Block_List,
					  Gaussian2D_Input_Parameters &Input_Parameters,
					  const int Number_of_Time_Steps,
					  const double &Time);

extern void Output_Drag(Gaussian2D_Quad_Block *Soln_ptr,
			AdaptiveBlock2D_List &Soln_Block_List,
			double &drag, double &lift);


extern void BCs(Gaussian2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
		Gaussian2D_Input_Parameters &IP);

extern double CFL(Gaussian2D_Quad_Block *Soln_ptr,
                  AdaptiveBlock2D_List &Soln_Block_List,
                  Gaussian2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(Gaussian2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List, 
                                const double &Dt_min);

extern double L1_Norm_Residual(Gaussian2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double L2_Norm_Residual(Gaussian2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double Max_Norm_Residual(Gaussian2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List);

extern void Evaluate_Limiters(Gaussian2D_Quad_Block *Soln_ptr,
                              AdaptiveBlock2D_List &Soln_Block_List);

extern void Freeze_Limiters(Gaussian2D_Quad_Block *Soln_ptr,
                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Apply_Boundary_Flux_Corrections(Gaussian2D_Quad_Block *Soln_ptr,
                                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Gaussian2D_Quad_Block *Soln_ptr,
                                                                AdaptiveBlock2D_List &Soln_Block_List,
                                                                Gaussian2D_Input_Parameters &Input_Parameters,
   	                                                        const int I_Stage);

extern int dUdt_Multistage_Explicit(Gaussian2D_Quad_Block *Soln_ptr,
                                    AdaptiveBlock2D_List &Soln_Block_List,
                                    Gaussian2D_Input_Parameters &Input_Parameters,
   	                            const int I_Stage);

extern int Update_Solution_Multistage_Explicit(Gaussian2D_Quad_Block *Soln_ptr,
                                               AdaptiveBlock2D_List &Soln_Block_List,
                                               Gaussian2D_Input_Parameters &Input_Parameters,
   	                                       const int I_Stage);

extern void Ramp_up_Reference_Mach_Number(Gaussian2D_Quad_Block *SolnBlk,
					  AdaptiveBlock2D_List &Soln_Block_List,
					  Gaussian2D_Input_Parameters &Input_Parameters,
					  const int Number_of_Time_Steps);


/**************************************************************************
 * Gaussian2D_Quad_Block -- Multiple Block External Subroutines for Mesh. *
 **************************************************************************/

extern Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                            Gaussian2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Broadcast_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                      Gaussian2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                             Gaussian2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                                            Gaussian2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                  Gaussian2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                 Gaussian2D_Input_Parameters &Input_Parameters);

extern int Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                          Gaussian2D_Input_Parameters &Input_Parameters);

extern int Output_Nodes_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                Gaussian2D_Input_Parameters &Input_Parameters);

extern int Output_Cells_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                Gaussian2D_Input_Parameters &Input_Parameters);

/**************************************************************************
 * Gaussian2D_Quad_Block -- Solvers.                                      *
 **************************************************************************/

extern int Gaussian2DQuadSolver(char *Input_File_Name_ptr,
                                int batch_flag);

#endif /* _GAUSSIAN2D_QUAD_INCLUDED  */
