/*!\file New_AdvectDiffuse2DQuad.h
  \brief Temporary header file defining 2D Advection Diffusion Equation Quadrilateral Mesh Solution Classes. */

#ifndef _NEW_ADVECTDIFFUSE2D_QUAD_INCLUDED
#define _NEW_ADVECTDIFFUSE2D_QUAD_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "AdvectDiffuse2DState.h"  /* Include 2D advection diffusion equation solution state header file */
#include "../Grid/Cell2D.h"        /* Include 2D cell header file */
#include "../Grid/Grid2DQuad.h"    /* Include 2D quadrilateral multiblock grid header file */
#include "../AMR/QuadTree.h"       /* Include quadtree header file */
#include "../AMR/AMR.h"            /* Include AMR header file */
#include "AdvectDiffuse2DInput.h"  /* Include 2D advection diffusion equation header file */
#include "../ICEM/ICEMCFD.h"       /* Include ICEMCFD header file. */

/* Define the structures and classes. */

#define	NUMBER_OF_RESIDUAL_VECTORS_ADVECTDIFFUSE2D    3

/*!
 * Class: AdvectDiffuse2D_Quad_Block_New
 *
 * @brief Class definition of the 2D advection-diffusion solution blocks.
 *
 * \verbatim
 * Member functions
 *       U      -- Return solution variable for the block (cell average).
 *    Grid      -- Return the solution block quadrilateral grid.
 *      dt      -- Return local time step for the solution block.
 *    dudt      -- Return the solution residuals.
 *    dudx      -- Return the unlimited solution gradients (x-direction)
 *                 for block.
 *    dudy      -- Return the unlimited solution gradients (y-direction)
 *                 for block.
 *     phi      -- Return the solution slope limiter.
 *    FluxN     -- Return array of north boundary solution fluxes.
 *    FluxS     -- Return array of south boundary solution fluxes.
 *    FluxE     -- Return array of east boundary solution fluxes.
 *    FluxW     -- Return array of west boundary solution fluxes.
 *     NCi      -- Return number of cells in
 *                 the i-direction (zeta-direction).
 *     ICl      -- Return lower index for cells in 
 *                 the i-direction (zeta-direction).
 *     ICu      -- Return upper index for cells in 
 *                 the i-direction (zeta-direction).
 *     NCj      -- Return number of cells in 
 *                 the j-direction (eta-direction).
 *     JCl      -- Return lower index for cells in 
 *                 the j-direction (eta-direction).
 *     JCu      -- Return upper index for cells in 
 *                 the j-direction (eta-direction).
 *  Nghost      -- Return number of ghost (halo or 
 *                 overlap) cells.
 * Axisymmetric -- Return axisymmetric geometry
 *                 indicator (=1 for axisymmetric,
 *                            =0 for planar geometry).
 * Freeze_Limiter -- Return limiter freezing indicator
 *                 (=1 for limiter freezing on,
 *                 (=0 for limiter freezing off).
 *     UoN      -- Return array of reference states for the application
 *                 of boundary conditions at the north boundary of the
 *                 solution block.
 *     UoS      -- Return array of reference states for the application
 *                 of boundary conditions at the south boundary of the
 *                 solution block.
 *     UoE      -- Return array of reference states for the application
 *                 of boundary conditions at the east boundary of the
 *                 solution block.
 *     UoW      -- Return array of reference states for the application
 *                 of boundary conditions at the west boundary of the
 *                 solution block.
 *   allocate   -- Allocate memory for structured quadrilateral
 *                 solution block.
 *   deallocate -- Deallocate memory for structured quadrilateral
 *                 solution block.
 *      Un      -- Return solution state at the specified node.
 *    UnNW      -- Return solution state at the north-west node.
 *    UnNE      -- Return solution state at the north-east node.
 *    UnSW      -- Return solution state at the south-west node.
 *    UnSE      -- Return solution state at the south-east node.
 *      un      -- Return solution value at the specified node.
 *    unNW      -- Return solution value at the north-west node.
 *    unNE      -- Return solution value at the north-east node.
 *    unSW      -- Return solution value at the south-west node.
 *    unSE      -- Return solution value at the south-east node.
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
 *      S -- a 2D advection diffusion equation solution
 *
 * S = S;
 * cout << S; (output function)
 * cin  >> S; (input function)
 * \endverbatim
 */
class AdvectDiffuse2D_Quad_Block_New{
private:
public:
  //! @name Solution state arrays:
  //@{
  AdvectDiffuse2D_State    **U; //!< Solution state.
  //@}

  //! @name Grid block information:
  //@{
  int                      NCi, //!< Total number of i-direction cells.
                           ICl, //!< First i-direction non-ghost cell counter.
                           ICu; //!< Final i-direction non-ghost cell counter.
  int                      NCj, //!< Total number of j-direction cells.
                           JCl, //!< First j-direction non-ghost cell counter.
                           JCu; //!< Final j-direction non-ghost cell counter.
  int                   Nghost; //!< Number of ghost cells.
  Grid2D_Quad_Block       Grid; //!< 2D quadrilateral grid geometry.
  //@}

  //! @name Residual and time-stepping arrays:
  //@{
  double                  **dt; //!< Local time step.
  double               ***dudt; //!< Solution residual.
  double                  **uo; //!< Initial solution.
  double                  **ue; //!< Exact solution.
  static int residual_variable; //!< Static integer that indicates which variable is used for residual calculations.
  //@}

  //! @name Solution gradient arrays:
  //@{
  double                **dudx; //!< Unlimited solution gradient (x-direction).
  double                **dudy; //!< Unlimited solution gradient (y-direction).
  double                 **phi; //!< Solution slope limiter.
  //@}

  //! @name Boundary solution flux arrays:
  //@{
  double                *FluxN, //!< North boundary solution flux.
                        *FluxS, //!< South boundary solution flux.
                        *FluxE, //!< East boundary solution flux.
                        *FluxW; //!< West boundary solution flux.
  //@}

  //! @name Problem indicator flags:
  //@{
  int             Axisymmetric; //!< Axisymmetric geometry indicator.
  int           Freeze_Limiter; //!< Limiter freezing indicator.
  static char*   solutionTitle; //!< Solution title info
  //@}

  //! @name Boundary condtion reference states:
  //@{
  AdvectDiffuse2D_State   *UoN, //!< Boundary condition reference states for north boundary.
                          *UoS, //!< Boundary condition reference states for south boundary.
                          *UoE, //!< Boundary condition reference states for east boundary.
                          *UoW; //!< Boundary condition reference states for west boundary.
  //@}

  //! @name Pointers to exact solutions
  //@{
  typedef Vector2D (* Exact_Gradient_Function) (const double, const double);
  static Exact_Gradient_Function ExactGrad; /* Exact gradient for equations with analytic solution */
  static FunctionType2D ExactSoln;          /* Exact solution for equations with analytic solution */
  //@}
	      
  //! @name Creation, copy, and assignment constructors.
  //@{
  //! Creation constructor.
  AdvectDiffuse2D_Quad_Block_New(void) {
    NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
    U = NULL; dt = NULL; dudt = NULL; 
    dudx = NULL; dudy = NULL; phi = NULL; uo = NULL; ue = NULL;
    FluxN = NULL; FluxS = NULL; FluxE = NULL; FluxW = NULL;
    UoN = NULL; UoS = NULL; UoE = NULL; UoW = NULL;
    Axisymmetric = 0; Freeze_Limiter = OFF;
  }

  //! Copy constructor.
  AdvectDiffuse2D_Quad_Block_New(const AdvectDiffuse2D_Quad_Block_New &Soln) {
    NCi = Soln.NCi; ICl = Soln.ICl; ICu = Soln.ICu; 
    NCj = Soln.NCj; JCl = Soln.JCl; JCu = Soln.JCu; Nghost = Soln.Nghost;
    Grid = Soln.Grid; U = Soln.U; dt = Soln.dt; dudt = Soln.dudt; 
    dudx = Soln.dudx; dudy = Soln.dudy; phi = Soln.phi;
    uo = Soln.uo; ue = Soln.ue;
    FluxN = Soln.FluxN; FluxS = Soln.FluxS; FluxE = Soln.FluxE; FluxW = Soln.FluxW;
    UoN = Soln.UoN; UoS = Soln.UoS; UoE = Soln.UoE; UoW = Soln.UoW;
    Axisymmetric = 0; Freeze_Limiter = Soln.Freeze_Limiter;
  }

  /* Destructor. */
  // ~AdvectDiffuse2D_Quad_Block_New(void);
  // Use automatically generated destructor.
  //@}

  /* Assignment operator. */
  // AdvectDiffuse2D_Quad_Block_New operator = (const AdvectDiffuse2D_Quad_Block_New &Soln);
  // Use automatically generated assignment operator.

  //! @name Allocate and deallocate functions.
  //@{
  //! Allocate memory for structured quadrilateral solution block.
    void allocate(const int Ni, const int Nj, const int Ng);

  //! Deallocate memory for structured quadrilateral solution block.
  void deallocate(void);
  //@}

  //! @name Bilinear interplation (Zingg & Yarrow).
  //@{
  //! Return solution state at specified node.
  AdvectDiffuse2D_State Un(const int ii, const int jj);

  AdvectDiffuse2D_State UnNW(const int ii, const int jj); //!< Return solution state at cell NW node.
  AdvectDiffuse2D_State UnNE(const int ii, const int jj); //!< Return solution state at cell NE node.
  AdvectDiffuse2D_State UnSE(const int ii, const int jj); //!< Return solution state at cell SE node.
  AdvectDiffuse2D_State UnSW(const int ii, const int jj); //!< Return solution state at cell SW node.

  //! Return solution state at specified node.
  double un(const int ii, const int jj);

  //! Return exact solution state at specified node.
  double uen(const int ii, const int jj);

  double unNW(const int ii, const int jj); //!< Return solution state at cell NW node.
  double unNE(const int ii, const int jj); //!< Return solution state at cell NE node.
  double unSE(const int ii, const int jj); //!< Return solution state at cell SE node.
  double unSW(const int ii, const int jj); //!< Return solution state at cell SW node.
  //@}

  //! @name Evaluate diffusive flux for the cell.
  //@{
  void evalDiffusiveFlux(const int ii, const int jj);
  //@}

  //! @name Member functions for limiter freezing.
  //@{
  void evaluate_limiters(void);
  void freeze_limiters(void);
  //@}

  //! @name Input-output operators.
  //@{
  friend ostream &operator << (ostream &out_file,
			       const AdvectDiffuse2D_Quad_Block_New &Soln);
  friend istream &operator >> (istream &in_file,
			       AdvectDiffuse2D_Quad_Block_New &Soln);
  //@}

  //! @name Member functions required for message passing.
  //@{
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

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::allocate -- Allocate memory.               *
 **************************************************************************/
inline void AdvectDiffuse2D_Quad_Block_New::allocate(const int Ni, const int Nj, const int Ng) {
   int i, j, k; assert(Ni > 1 && Nj > 1 && Ng > 1 && Ng > 1); Grid.allocate(Ni, Nj, Ng);
   NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1;
   NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1; Nghost = Ng;
   U = new AdvectDiffuse2D_State*[NCi]; dt = new double*[NCi]; dudt = new double**[NCi]; 
   dudx = new double*[NCi]; dudy = new double*[NCi]; 
   phi = new double*[NCi]; uo = new double*[NCi]; ue = new double*[NCi];
   for ( i = 0; i <= NCi-1 ; ++i ) {
      U[i] = new AdvectDiffuse2D_State[NCj]; 
      dt[i] = new double[NCj]; dudt[i] = new double*[NCj];
      for ( j = 0; j <= NCj-1 ; ++j ) 
        { dudt[i][j] = new double[NUMBER_OF_RESIDUAL_VECTORS_ADVECTDIFFUSE2D]; }
      dudx[i] = new double[NCj]; dudy[i] = new double[NCj]; 
      phi[i] = new double[NCj];
      uo[i] = new double[NCj]; ue[i] = new double[NCj];
   } /* endfor */
   FluxN = new double[NCi]; FluxS = new double[NCi];
   FluxE = new double[NCj]; FluxW = new double[NCj];
   UoN = new AdvectDiffuse2D_State[NCi]; UoS = new AdvectDiffuse2D_State[NCi];
   UoE = new AdvectDiffuse2D_State[NCj]; UoW = new AdvectDiffuse2D_State[NCj];
   // Set the solution residuals, gradients, limiters, and other values to zero.
   for (j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
      for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
         for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_ADVECTDIFFUSE2D-1 ; ++k ) { dudt[i][j][k] = ZERO; }
	 dudx[i][j] = ZERO; dudy[i][j] = ZERO; phi[i][j] = ZERO; uo[i][j] = ZERO; dt[i][j] = ZERO;
      } /* endfor */
   } /* endfor */
}

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::deallocate -- Deallocate memory.           *
 **************************************************************************/
inline void AdvectDiffuse2D_Quad_Block_New::deallocate(void) {
   int i, j; Grid.deallocate();
   for ( i = 0; i <= NCi-1 ; ++i ) {
      delete []U[i]; U[i] = NULL;
      delete []dt[i]; dt[i] = NULL; 
      for ( j = 0; j <= NCj-1 ; ++j ) { delete []dudt[i][j]; dudt[i][j] = NULL; }
      delete []dudt[i]; dudt[i] = NULL;
      delete []dudx[i]; dudx[i] = NULL; delete []dudy[i]; dudy[i] = NULL;
      delete []phi[i]; phi[i] = NULL; delete []uo[i]; uo[i] = NULL;
      delete []ue[i]; ue[i] = NULL;
   } /* endfor */
   delete []U; U = NULL; delete []dt; dt = NULL; delete []dudt; dudt = NULL;
   delete []dudx; dudx = NULL; delete []dudy; dudy = NULL; 
   delete []phi; phi = NULL; delete []uo; uo = NULL; delete []ue; ue = NULL;
   delete []FluxN; FluxN = NULL; delete []FluxS; FluxS = NULL;
   delete []FluxE; FluxE = NULL; delete []FluxW; FluxW = NULL;
   delete []UoN; UoN = NULL; delete []UoS; UoS = NULL;
   delete []UoE; UoE = NULL; delete []UoW; UoW = NULL;
   NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
}

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::Un -- Node solution state.                 *
 **************************************************************************/
inline AdvectDiffuse2D_State AdvectDiffuse2D_Quad_Block_New::Un(const int ii, const int jj) {
    double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
           eta1, zeta1, eta2, zeta2, eta, zeta;
    AdvectDiffuse2D_State Unew; double As, Bs, Cs, Ds; Vector2D Av, Bv, Cv, Dv;
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
    Unew.u = As+Bs*zeta+Cs*eta+Ds*zeta*eta;
    Av=U[ii-1][jj-1].V; Bv=U[ii-1][jj].V-U[ii-1][jj-1].V; Cv=U[ii][jj-1].V-U[ii-1][jj-1].V;
    Dv=U[ii][jj].V+U[ii-1][jj-1].V-U[ii-1][jj].V-U[ii][jj-1].V;
    Unew.V = Av+Bv*zeta+Cv*eta+Dv*zeta*eta;
    As=U[ii-1][jj-1].k; Bs=U[ii-1][jj].k-U[ii-1][jj-1].k; Cs=U[ii][jj-1].k-U[ii-1][jj-1].k;
    Ds=U[ii][jj].k+U[ii-1][jj-1].k-U[ii-1][jj].k-U[ii][jj-1].k;
    Unew.k = As+Bs*zeta+Cs*eta+Ds*zeta*eta;
    As=U[ii-1][jj-1].T; Bs=U[ii-1][jj].T-U[ii-1][jj-1].T; Cs=U[ii][jj-1].T-U[ii-1][jj-1].T;
    Ds=U[ii][jj].T+U[ii-1][jj-1].T-U[ii-1][jj].T-U[ii][jj-1].T;
    Unew.T = As+Bs*zeta+Cs*eta+Ds*zeta*eta;
    Av=U[ii-1][jj-1].Fd; Bv=U[ii-1][jj].Fd-U[ii-1][jj-1].Fd; Cv=U[ii][jj-1].Fd-U[ii-1][jj-1].Fd;
    Dv=U[ii][jj].Fd+U[ii-1][jj-1].Fd-U[ii-1][jj].Fd-U[ii][jj-1].Fd;
    Unew.Fd = Av+Bv*zeta+Cv*eta+Dv*zeta*eta;
    return (Unew);
}

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::Un?? -- Get cell node solution states.     *
 **************************************************************************/
inline AdvectDiffuse2D_State AdvectDiffuse2D_Quad_Block_New::UnNW(const int ii, const int jj) {
  return (Un(ii, jj+1));
}

inline AdvectDiffuse2D_State AdvectDiffuse2D_Quad_Block_New::UnNE(const int ii, const int jj) {
  return (Un(ii+1, jj+1));
}

inline AdvectDiffuse2D_State AdvectDiffuse2D_Quad_Block_New::UnSE(const int ii, const int jj) {
  return (Un(ii+1, jj));
}

inline AdvectDiffuse2D_State AdvectDiffuse2D_Quad_Block_New::UnSW(const int ii, const int jj) {
  return (Un(ii, jj));
}

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::un -- Node solution value.                 *
 **************************************************************************/
inline double AdvectDiffuse2D_Quad_Block_New::un(const int ii, const int jj) {
    double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
           eta1, zeta1, eta2, zeta2, eta, zeta;
    double As, Bs, Cs, Ds;
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
    return (As+Bs*zeta+Cs*eta+Ds*zeta*eta);
}

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::uen -- Exact solution node value.          *
 **************************************************************************/
inline double AdvectDiffuse2D_Quad_Block_New::uen(const int ii, const int jj) {
    double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
           eta1, zeta1, eta2, zeta2, eta, zeta;
    double As, Bs, Cs, Ds;
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
    As=ue[ii-1][jj-1]; Bs=ue[ii-1][jj]-ue[ii-1][jj-1]; Cs=ue[ii][jj-1]-ue[ii-1][jj-1];
    Ds=ue[ii][jj]+ue[ii-1][jj-1]-ue[ii-1][jj]-ue[ii][jj-1];
    return (As+Bs*zeta+Cs*eta+Ds*zeta*eta);
}

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::un?? -- Get cell node solution values.     *
 **************************************************************************/
inline double AdvectDiffuse2D_Quad_Block_New::unNW(const int ii, const int jj) {
  return (un(ii, jj+1));
}

inline double AdvectDiffuse2D_Quad_Block_New::unNE(const int ii, const int jj) {
  return (un(ii+1, jj+1));
}

inline double AdvectDiffuse2D_Quad_Block_New::unSE(const int ii, const int jj) {
  return (un(ii+1, jj));
}

inline double AdvectDiffuse2D_Quad_Block_New::unSW(const int ii, const int jj) {
  return (un(ii, jj));
}

/*****************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::evalDiffusiveFlux -- Evaluate diffusive flux. *
 *****************************************************************************/
inline void AdvectDiffuse2D_Quad_Block_New::evalDiffusiveFlux(const int ii, const int jj) {
  U[ii][jj].Fd = U[ii][jj].F_diff(dudx[ii][jj], dudy[ii][jj]);
}

/**********************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::evaluate_limiters -- Set flag to evaluate limiters.*
 **********************************************************************************/
inline void AdvectDiffuse2D_Quad_Block_New::evaluate_limiters(void) {
  Freeze_Limiter = OFF; 
}

/**********************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::freeze_limiters -- Set flag to freeze limiters.    *
 **********************************************************************************/
inline void AdvectDiffuse2D_Quad_Block_New::freeze_limiters(void) {
  Freeze_Limiter = ON; 
}

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block_New -- Input-output operators.                  *
 **************************************************************************/
inline ostream &operator << (ostream &out_file,
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

inline istream &operator >> (istream &in_file,
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
	 SolnBlk.ue[i][j] = ZERO;
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

/*******************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::NumVar -- Returns number of state variables.    *
 *******************************************************************************/
inline int AdvectDiffuse2D_Quad_Block_New::NumVar(void) {
  return (int(NUM_VAR_ADVECTDIFFUSE2D));
}

/*******************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::LoadSendBuffer -- Loads send message buffer.    *
 *******************************************************************************/
inline int AdvectDiffuse2D_Quad_Block_New::LoadSendBuffer(double *buffer,
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
 * AdvectDiffuse2D_Quad_Block_New::LoadSendBuffer_F2C -- Loads send message buffer *
 *                                                   for fine to coarse block  *
 *                                                   message passing.          *
 *******************************************************************************/
inline int AdvectDiffuse2D_Quad_Block_New::LoadSendBuffer_F2C(double *buffer,
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
 * AdvectDiffuse2D_Quad_Block_New::LoadSendBuffer_C2F -- Loads send message buffer *
 *                                                   for coarse to fine block  *
 *                                                   message passing.          *
 *******************************************************************************/
inline int AdvectDiffuse2D_Quad_Block_New::LoadSendBuffer_C2F(double *buffer,
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
 * AdvectDiffuse2D_Quad_Block_New::UnloadReceiveBuffer -- Unloads receive buffer.  *
 *******************************************************************************/
inline int AdvectDiffuse2D_Quad_Block_New::UnloadReceiveBuffer(double *buffer,
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
 * AdvectDiffuse2D_Quad_Block_New::UnloadReceiveBuffer_F2C -- Unloads receive message  *
 *                                                        buffer for fine to coarse*
 *                                                        block message passing.   *
 ***********************************************************************************/
inline int AdvectDiffuse2D_Quad_Block_New::UnloadReceiveBuffer_F2C(double *buffer,
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
 * AdvectDiffuse2D_Quad_Block_New::UnloadReceiveBuffer_C2F -- Unloads receive message  *
 *                                                        buffer for coarse to fine*
 *                                                        block message passing.   *
 ***********************************************************************************/
inline int AdvectDiffuse2D_Quad_Block_New::UnloadReceiveBuffer_C2F(double *buffer,
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

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::SubcellReconstruction --                   *
 *               Performs the subcell reconstruction of solution state    *
 *               within a given cell (i,j) of the computational mesh for  *
 *               the specified quadrilateral solution block.              *
 **************************************************************************/
inline void AdvectDiffuse2D_Quad_Block_New::SubcellReconstruction(const int i, 
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

/*******************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::LoadSendBuffer_Flux_F2C --                      *
 *                 Loads send message buffer for                               *
 *                 fine to coarse block message                                *
 *                 passing of conservative                                     *
 *                 solution fluxes.                                            *
 *******************************************************************************/
inline int AdvectDiffuse2D_Quad_Block_New::LoadSendBuffer_Flux_F2C(double *buffer,
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

/*******************************************************************************
 * AdvectDiffuse2D_Quad_Block_New::UnloadReceiveBuffer_Flux_F2C --                 *
 *                 Unloads receive message                                     *
 *                 buffer for fine to coarse                                   *
 *                 block message passing of                                    *
 *                 conservative solution fluxes.                               *
 *******************************************************************************/
inline int AdvectDiffuse2D_Quad_Block_New::UnloadReceiveBuffer_Flux_F2C(double *buffer,
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

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block_New -- Single Block External Subroutines.       *
 **************************************************************************/
#if 0

extern void Write_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
	                         ostream &Out_File);

extern void Read_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
	                        istream &In_File);

extern void Broadcast_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk);

#ifdef _MPI_VERSION
extern void Broadcast_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
                                     MPI::Intracomm &Communicator, 
                                     const int Source_CPU);
#endif

extern void Copy_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk1,
		                AdvectDiffuse2D_Quad_Block_New &SolnBlk2);

extern int Prolong_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk_Fine,
				  AdvectDiffuse2D_Quad_Block_New &SolnBlk_Original,
				  const int Sector);

extern int Restrict_Solution_Block(AdvectDiffuse2D_Quad_Block_New &SolnBlk_Coarse,
				   AdvectDiffuse2D_Quad_Block_New &SolnBlk_Original_SW,
				   AdvectDiffuse2D_Quad_Block_New &SolnBlk_Original_SE,
				   AdvectDiffuse2D_Quad_Block_New &SolnBlk_Original_NW,
				   AdvectDiffuse2D_Quad_Block_New &SolnBlk_Original_NE);

extern void Output_Tecplot(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
			   AdvectDiffuse2D_Input_Parameters &IP,
		           const int Number_of_Time_Steps,
                           const double &Time,
                           const int Block_Number,
                           const int Output_Title,
	                   ostream &Out_File);

extern void Output_Cells_Tecplot(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
				 AdvectDiffuse2D_Input_Parameters &IP,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_Nodes_Tecplot(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void ICs(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
 	        const int i_ICtype,
                AdvectDiffuse2D_State *Wo);

extern void Set_Boundary_Ref_State(AdvectDiffuse2D_Quad_Block_New &SolnBlk, 
				   const int GlobalSolnBlkNum,
				   const QuadTreeBlock_DataStructure &QuadTree,
				   const AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern void Set_Analytical_Solution(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
				    const AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern void Set_Advection_Velocity_Field(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
                                         const int i_Velocity_Field,
                                         const double &Vx,
                                         const double &Vy);

extern void BCs(AdvectDiffuse2D_Quad_Block_New &SolnBlk);

extern double CFL(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
                  AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(AdvectDiffuse2D_Quad_Block_New &SolnBlk, 
                                const double &Dt_min);

extern double L1_Norm_Residual(AdvectDiffuse2D_Quad_Block_New &SolnBlk);

extern double L2_Norm_Residual(AdvectDiffuse2D_Quad_Block_New &SolnBlk);

extern double Max_Norm_Residual(AdvectDiffuse2D_Quad_Block_New &SolnBlk);

extern void Linear_Reconstruction_GreenGauss(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
                                             const int i,
                                             const int j,
					     const int Limiter);

extern void Linear_Reconstruction_GreenGauss2(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
                                              const int i,
                                              const int j,
					      const int Limiter);

extern void Linear_Reconstruction_GreenGauss(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
					     const int Limiter);

extern void Linear_Reconstruction_LeastSquares(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
                                               const int i,
                                               const int j,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
					       const int Limiter);

extern void Diffusive_Flux(AdvectDiffuse2D_Quad_Block_New &SolnBlk);

extern void Residual_Smoothing(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
                               const int k_residual,
			       double &epsilon, 
                               const int number_of_Gauss_Seidel_iterations);

extern void Calculate_Refinement_Criteria(double *refinement_criteria,
					  AdvectDiffuse2D_Input_Parameters &IP,
                                          int &number_refinement_criteria,
                                          AdvectDiffuse2D_Quad_Block_New &SolnBlk);

extern void Fix_Refined_Block_Boundaries(AdvectDiffuse2D_Quad_Block_New SolnBlk,
                                         const int Fix_North_Boundary,
                                         const int Fix_South_Boundary,
                                         const int Fix_East_Boundary,
                                         const int Fix_West_Boundary);

extern void Unfix_Refined_Block_Boundaries(AdvectDiffuse2D_Quad_Block_New SolnBlk);

extern void Apply_Boundary_Flux_Corrections(AdvectDiffuse2D_Quad_Block_New,
                                            const int Number_Neighbours_North_Boundary,
                                            const int Number_Neighbours_South_Boundary,
                                            const int Number_Neighbours_East_Boundary,
                                            const int Number_Neighbours_West_Boundary);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(AdvectDiffuse2D_Quad_Block_New SolnBlk,
                                                                const int i_stage,
                                                                const int n_stage,
                                                                const double &CFL_Number,
                                                                const int Time_Integration_Type,
                                                                const int Reconstruction_Type,
                                                                const int Limiter_Type,
                                                                const int Number_Neighbours_North_Boundary,
                                                                const int Number_Neighbours_South_Boundary,
                                                                const int Number_Neighbours_East_Boundary,
                                                                const int Number_Neighbours_West_Boundary);

extern int dUdt_Residual_Evaluation(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
				    AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int dUdt_Multistage_Explicit(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
   	                            const int i_stage,
                                    AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Update_Solution_Multistage_Explicit(AdvectDiffuse2D_Quad_Block_New &SolnBlk,
   	                                       const int i_stage,
                                               AdvectDiffuse2D_Input_Parameters &Input_Parameters);

/**************************************************************************
 * AdvectDiffuse2D_Quad_Block_New -- Multiple Block External Subroutines.     *
 **************************************************************************/

extern AdvectDiffuse2D_Quad_Block_New* Allocate(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                            AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern AdvectDiffuse2D_Quad_Block_New* Deallocate(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                              AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern void ICs(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
		AdaptiveBlock2D_List &Soln_Block_List,
		AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		QuadTreeBlock_DataStructure &QuadTree);

extern void Set_Advection_Velocity_Field(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                         AdaptiveBlock2D_List &Soln_Block_List,
                                         AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern void Set_Analytical_Solution(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
				    const AdaptiveBlock2D_List &Soln_Block_List,
				    const AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Read_Restart_Solution(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                 AdaptiveBlock2D_List &Soln_Block_List,
                                 AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		                 int &Number_of_Time_Steps,
                                 double &Time,
                                 CPUTime &CPU_Time);

extern int Write_Restart_Solution(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                  AdaptiveBlock2D_List &Soln_Block_List,
                                  AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		                  const int Number_of_Time_Steps,
                                  const double &Time,
                                  const CPUTime &CPU_Time);

extern int Output_Tecplot(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		          const int Number_of_Time_Steps,
                          const double &Time);

extern int Output_Cells_Tecplot(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Nodes_Tecplot(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Mesh_Tecplot(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Mesh_Gnuplot(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern void BCs(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
		AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern double CFL(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                  AdaptiveBlock2D_List &Soln_Block_List,
                  AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List, 
                                const double &Dt_min);

extern double L1_Norm_Residual(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double L2_Norm_Residual(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double Max_Norm_Residual(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List);

extern void Evaluate_Limiters(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                              AdaptiveBlock2D_List &Soln_Block_List);

extern void Freeze_Limiters(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Residual_Smoothing(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               AdvectDiffuse2D_Input_Parameters &Input_Parameters,
   	                       const int I_Stage);

extern void Apply_Boundary_Flux_Corrections(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                                                AdaptiveBlock2D_List &Soln_Block_List,
                                                                AdvectDiffuse2D_Input_Parameters &Input_Parameters,
   	                                                        const int I_Stage);

extern int dUdt_Multistage_Explicit(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
				    AdaptiveBlockResourceList &Global_Soln_Block_List,
                                    AdaptiveBlock2D_List &Local_Soln_Block_List,
                                    AdvectDiffuse2D_Input_Parameters &Input_Parameters,
   	                            const int I_Stage);

extern int Update_Solution_Multistage_Explicit(AdvectDiffuse2D_Quad_Block_New *Soln_ptr,
                                               AdaptiveBlock2D_List &Soln_Block_List,
                                               AdvectDiffuse2D_Input_Parameters &Input_Parameters,
   	                                       const int I_Stage);

extern int Initial_AMR2(AdvectDiffuse2D_Quad_Block_New       *Soln_ptr,
                        AdvectDiffuse2D_Input_Parameters &InputParameters,
                        QuadTreeBlock_DataStructure      &QuadTree,
                        AdaptiveBlockResourceList        &GlobalSolnBlockList,
                        AdaptiveBlock2D_List             &LocalSolnBlockList);

extern int Uniform_AMR2(AdvectDiffuse2D_Quad_Block_New       *Soln_ptr,
                        AdvectDiffuse2D_Input_Parameters &InputParameters,
                        QuadTreeBlock_DataStructure      &QuadTree,
                        AdaptiveBlockResourceList        &GlobalSolnBlockList,
                        AdaptiveBlock2D_List             &LocalSolnBlockList);

extern int Boundary_AMR2(AdvectDiffuse2D_Quad_Block_New       *Soln_ptr,
                         AdvectDiffuse2D_Input_Parameters &InputParameters,
                         QuadTreeBlock_DataStructure      &QuadTree,
                         AdaptiveBlockResourceList        &GlobalSolnBlockList,
                         AdaptiveBlock2D_List             &LocalSolnBlockList);

/********************************************************************************
 * AdvectDiffuse2D_Quad_Block_New -- Multiple Block External Subroutines for Mesh.  *
 ********************************************************************************/

extern Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                            AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Broadcast_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                      AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                             AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                                            AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                  AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                 AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                          AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Output_Nodes_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                AdvectDiffuse2D_Input_Parameters &Input_Parameters);

extern int Output_Cells_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                AdvectDiffuse2D_Input_Parameters &Input_Parameters);

#endif


#endif /* _ADVECTDIFFUSE2D_QUAD_INCLUDED  */
