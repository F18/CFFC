/**********************************************************************
 * NavierStokes2DQuad.h: Header file defining 2D Navier-Stokes        *
 *                       quadrilateral mesh solution classes.         *
 **********************************************************************/

#ifndef _NAVIERSTOKES2D_QUAD_INCLUDED
#define _NAVIERSTOKES2D_QUAD_INCLUDED

// Include 2D Navier-Stokes state header file.

#ifndef _NAVIERSTOKES2D_STATE_INCLUDED
#include "NavierStokes2DState.h"
#endif // _NAVIERSTOKES2D_STATE_INCLUDED

// Include 2D Navier-Stokes input header file.

#ifndef _NAVIERSTOKES2D_INPUT_INCLUDED
#include "NavierStokes2DInput.h"
#endif // _NAVIERSTOKES2D_INPUT_INCLUDED

// Include the 2D cell header file.

#ifndef _CELL2D_INCLUDED
#include "../Grid/Cell2D.h"
#endif // _CELL2D_INCLUDED

// Include the 2D quadrilateral multiblock grid header file.

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

// Include the turbulent wall data header file.

#ifndef _TURBULENT2D_WALLDATA_INCLUDED
#include "../Turbulent2D/Turbulent2DWallData.h"
#endif // _TURBULENT2D_WALLDATA_INCLUDED

// Define the structures and classes.

#define	NUMBER_OF_RESIDUAL_VECTORS_NAVIERSTOKES2D  3

// (Parallel) Debugging option.

#define _NS_PARALLEL_DEBUG_

/*!
 * Class: NavierStokes2D_Quad_Block
 *
 * @brief Class definition of the 2D Navier-Stokes solution blocks.
 *
 * \verbatim
 * Member functions
 *       W      -- Return primitive variable solution for the block
 *                 (cell average).
 *       U      -- Return conserved variable solution for the block
 *                 (cell average).
 *    Wall      -- Return turbulent wall data.
 *    Grid      -- Return the solution block quadrilateral grid or mesh.
 *      dt      -- Return local time step for the solution block.
 *    dUdt      -- Return the solution block residuals.
 *    dWdx      -- Return the unlimited primitive variable solution
 *                 gradients (x-direction) for the block.
 *    dWdy      -- Return the unlimited primitive variable solution
 *                 gradients (y-direction) for the block.
 *     phi      -- Return the solution slope limiters.
 *      Uo      -- Return initial solution state.
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
 *     NCj      -- Return number of cells in the j-direction 
 *                 (eta-direction).
 *     JCl      -- Return lower index for cells in the j-direction 
 *                 (eta-direction).
 *     JCu      -- Return upper index for cells in the j-direction 
 *                 (eta-direction).
 *  Nghost      -- Number of ghost cells.
 * Axisymmetric -- Return axisymmetric flow indicator (= 1 for 
 *                 axisymmetric flow, = 0 for planar flow).
 * Freeze_Limiter -- Return limiter freezing indicator
 *                   (=1 for limiter freezing on,
 *                   (=0 for limiter freezing off).
 *     WoN      -- Return array of reference states for the
 *                 application of boundary conditions at the north
 *                 boundary of the solution block.
 *     WoS      -- Return array of reference states for the
 *                 application of boundary conditions at the south
 *                 boundary of the solution block.
 *     WoE      -- Return array of reference states for the
 *                 application of boundary conditions at the east
 *                 boundary of the solution block.
 *     WoW      -- Return array of reference states for the
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
 * SubcellReconstruction -- Performs linear subcell solution
 *                          reconstruction used in adaptive mesh
 *                          refinement.
 * LoadSendBuffer_Flux_F2C -- Loads send buffer for sending
 *                            conservative flux corrections from fine
 *                            to coarse solution blocks.
 * UnLoadSendBuffer_Flux_F2C -- Loads send buffer for sending
 *                              conservative flux corrections from
 *                              fine to coarse solution blocks.
 *
 * Member operators
 *      S -- a 2D NavierStokes solution
 *
 * S = S;
 * cout << S; (output function)
 * cin  >> S; (input function)
 * \endverbatim
 */
class NavierStokes2D_Quad_Block{
private:
public:
  //@{ @name Solution state arrays:
  NavierStokes2D_pState     **W; //!< Primitive solution state.
  NavierStokes2D_cState     **U; //!< Conserved solution state.
  //@}

  //@{ @name Grid block information:
  int                       NCi, //!< Total number of i-direction cells.
                            ICl, //!< First i-direction non-ghost cell counter.
                            ICu; //!< Final i-direction non-ghost cell counter.
  int                       NCj, //!< Total number of j-direction cells.
                            JCl, //!< First j-direction non-ghost cell counter.
                            JCu; //!< Final j-direction non-ghost cell counter.
  int                    Nghost; //!< Number of ghost cells.
  Grid2D_Quad_Block        Grid; //!< 2D quadrilateral grid geometry.
  //@}

  //@{ @name Residual and time-stepping arrays:
  double                   **dt; //!< Local time step.
  NavierStokes2D_cState ***dUdt; //!< Solution residual.
  NavierStokes2D_cState    **Uo; //!< Initial solution state.
  static int  residual_variable; //!< Static integer that indicates which variable is used for residual calculations.
  //@}

  //@{ @name Solution gradient arrays:
  NavierStokes2D_pState  **dWdx; //!< Unlimited solution gradient (x-direction).
  NavierStokes2D_pState  **dWdy; //!< Unlimited solution gradient (y-direction).
  NavierStokes2D_pState   **phi; //!< Solution slope limiter.
  //@}

  //@{ @name Boundary solution flux arrays:
  NavierStokes2D_cState  *FluxN, //!< North boundary solution flux.
                         *FluxS, //!< South boundary solution flux.
                         *FluxE, //!< East boundary solution flux.
                         *FluxW; //!< West boundary solution flux.
  //@}

  //@{ @name Problem indicator flags:
  int              Axisymmetric; //!< Axisymmetric flow indicator.
  int                 Flow_Type; //!< Flow-type flag (inviscid, laminar, or k-omega).
  int            Freeze_Limiter; //!< Limiter freezing indicator.
  //@}

  //@{ @name Boundary condtion reference states:
  NavierStokes2D_pState    *WoN, //!< Boundary condition reference states for north boundary.
                           *WoS, //!< Boundary condition reference states for south boundary.
                           *WoE, //!< Boundary condition reference states for east boundary.
                           *WoW; //!< Boundary condition reference states for west boundary.
  //@}

  //@{ @name Turbulence wall data arrays:
  Turbulent2DWallData    **Wall; //!< Turbulent wall data.
  //@}

  //@{ @name Flow constants:
  Vector2D                Vwall; //!< Specified wall velocity.
  double                  Twall; //!< Specified wall temperature.
  //@}

  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  NavierStokes2D_Quad_Block(void) {
    // Problem flags:
    Axisymmetric   = OFF;
    Flow_Type      = FLOWTYPE_INVISCID;
    Freeze_Limiter = OFF;
    Vwall = Vector2D_ZERO;
    Twall = ZERO;
    // Grid size and variables:
    NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
    // Solution variables:
    W = NULL; U = NULL;
    dt = NULL; dUdt = NULL; Uo = NULL;
    dWdx = NULL; dWdy = NULL; phi = NULL;
    FluxN = NULL; FluxS = NULL; FluxE = NULL; FluxW = NULL;
    WoN = NULL; WoS = NULL; WoE = NULL; WoW = NULL;
    // Turbulent wall data:
    Wall = NULL;
  }

  //! Copy constructor.
  NavierStokes2D_Quad_Block(const NavierStokes2D_Quad_Block &Soln) {
    // Problem flags:
    Axisymmetric   = Soln.Axisymmetric;
    Flow_Type      = Soln.Flow_Type;
    Freeze_Limiter = Soln.Freeze_Limiter;
    Vwall = Soln.Vwall;
    Twall = Soln.Twall;
    // Grid size and variables:
    NCi = Soln.NCi; ICl = Soln.ICl; ICu = Soln.ICu;
    NCj = Soln.NCj; JCl = Soln.JCl; JCu = Soln.JCu;
    Nghost = Soln.Nghost;
    Grid  = Soln.Grid;
    // Solution variables:
    W = Soln.W; U = Soln.U;
    dt = Soln.dt; dUdt = Soln.dUdt; Uo = Soln.Uo;
    dWdx = Soln.dWdx; dWdy = Soln.dWdy; phi = Soln.phi;
    FluxN = Soln.FluxN; FluxS = Soln.FluxS; 
    FluxE = Soln.FluxE; FluxW = Soln.FluxW;
    WoN   = Soln.WoN;   WoS   = Soln.WoS;
    WoE   = Soln.WoE;   WoW   = Soln.WoW;
    // Turbulent wall data:
    Wall = Soln.Wall;
  }

  // Destructor.
  // ~NavierStokes2D_Quad_Block(void);
  // Use automatically generated destructor.
  //@}

  // Assignment operator.
  // NavierStokes2D_Quad_Block operator = (const NavierStokes2D_Quad_Block &Soln);
  // Use automatically generated assignment operator.

  //@{ @name Allocate and deallocate functions.
  //! Allocate memory for structured quadrilateral solution block.
  void allocate(const int Ni, const int Nj, const int Ng);
  //! Deallocate memory for structured quadrilateral solution block.
  void deallocate(void);
  //@}

  //@{ @name Bilinear interplation (Zingg & Yarrow).
  //! Return primitive solution state at specified node.
  NavierStokes2D_pState Wn(const int ii, const int jj);

  //! Return conserverd solution state at specified node.
  NavierStokes2D_cState Un(const int ii, const int jj);

  NavierStokes2D_pState WnNW(const int ii, const int jj); //!< Return primitive solution state at cell NW node.
  NavierStokes2D_pState WnNE(const int ii, const int jj); //!< Return primitive solution state at cell NE node.
  NavierStokes2D_pState WnSE(const int ii, const int jj); //!< Return primitive solution state at cell SE node.
  NavierStokes2D_pState WnSW(const int ii, const int jj); //!< Return primitive solution state at cell SW node.

  NavierStokes2D_cState UnNW(const int ii, const int jj); //!< Return conserved solution state at cell NW node.
  NavierStokes2D_cState UnNE(const int ii, const int jj); //!< Return conserved solution state at cell NE node.
  NavierStokes2D_cState UnSE(const int ii, const int jj); //!< Return conserved solution state at cell SE node.
  NavierStokes2D_cState UnSW(const int ii, const int jj); //!< Return conserved solution state at cell SW node.
  //@}

  //@{ @name Bilinear interplation (Connell & Holmes).
  //! Return primitive solution state at specified node.
  NavierStokes2D_pState Ww(const int ii, const int jj);

  //! Retern conserverd solution state at specified node.
  NavierStokes2D_cState Uw(const int ii, const int jj);

  //! Return primitive solution state at cell nodes.
  NavierStokes2D_pState WwNW(const int ii, const int jj); //!< Return primitive solution state at cell NW node.
  NavierStokes2D_pState WwNE(const int ii, const int jj); //!< Return primitive solution state at cell NE node.
  NavierStokes2D_pState WwSE(const int ii, const int jj); //!< Return primitive solution state at cell SE node.
  NavierStokes2D_pState WwSW(const int ii, const int jj); //!< Return primitive solution state at cell SW node.

  //! Return conserved solution state at cell nodes.
  NavierStokes2D_cState UwNW(const int ii, const int jj); //!< Return conserved solution state at cell NW node.
  NavierStokes2D_cState UwNE(const int ii, const int jj); //!< Return conserved solution state at cell NE node.
  NavierStokes2D_cState UwSE(const int ii, const int jj); //!< Return conserved solution state at cell SE node.
  NavierStokes2D_cState UwSW(const int ii, const int jj); //!< Return conserved solution state at cell SW node.
  //@}

  //@{ @name Member functions for limiter freezing.
  void evaluate_limiters(void); //!< Set flags for limiter evaluation.
  void freeze_limiters(void);   //!< Set flags for limiter freezing.
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const NavierStokes2D_Quad_Block &Soln);
  friend istream &operator >> (istream &in_file, NavierStokes2D_Quad_Block &Soln);
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

#ifdef _NS_PARALLEL_DEBUG_
  char extension[256], output_file_name[256];
  char *output_file_name_ptr;
  static ofstream dout;

  int Open_Debug_Output_File(const int &ThisCPU) {
    strcpy(output_file_name,"debug");
    strcat(output_file_name,"_cpu");
    sprintf(extension,"%.6d",ThisCPU);
    strcat(extension,".txt");
    strcat(output_file_name,extension);
    output_file_name_ptr = output_file_name;
    dout.open(output_file_name_ptr,ios::out);
    if (dout.bad()) return 1;
    return 0;
  }

  void Close_Debug_Output_File(void) {
    dout.close();
  }
#endif

};

/**********************************************************************
 * NavierStokes2D_Quad_Block::allocate -- Allocate memory.            *
 **********************************************************************/
inline void NavierStokes2D_Quad_Block::allocate(const int Ni, const int Nj, const int Ng) {
  assert(Ni > 1 && Nj > 1 && Ng > 1);
  Grid.allocate(Ni,Nj,Ng);
  NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1; 
  NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1; Nghost = Ng;
  W = new NavierStokes2D_pState*[NCi]; U = new NavierStokes2D_cState*[NCi];
  dt = new double*[NCi]; dUdt = new NavierStokes2D_cState**[NCi];
  dWdx = new NavierStokes2D_pState*[NCi]; dWdy = new NavierStokes2D_pState*[NCi];
  phi = new NavierStokes2D_pState*[NCi]; Uo = new NavierStokes2D_cState*[NCi];
  Wall = new Turbulent2DWallData*[NCi];
  for (int i = 0; i < NCi; i++) {
    W[i] = new NavierStokes2D_pState[NCj]; U[i] = new NavierStokes2D_cState[NCj];
    dt[i] = new double[NCj]; dUdt[i] = new NavierStokes2D_cState*[NCj];
    for (int j = 0; j < NCj; j++)
      dUdt[i][j] = new NavierStokes2D_cState[NUMBER_OF_RESIDUAL_VECTORS_NAVIERSTOKES2D];
    dWdx[i] = new NavierStokes2D_pState[NCj]; dWdy[i] = new NavierStokes2D_pState[NCj];
    phi[i] = new NavierStokes2D_pState[NCj]; Uo[i] = new NavierStokes2D_cState[NCj];
    Wall[i] = new Turbulent2DWallData[NCj];
  }
  WoN = new NavierStokes2D_pState[NCi]; WoS = new NavierStokes2D_pState[NCi];
  WoE = new NavierStokes2D_pState[NCj]; WoW = new NavierStokes2D_pState[NCj];
  FluxN = new NavierStokes2D_cState[NCi]; FluxS = new NavierStokes2D_cState[NCi];
  FluxE = new NavierStokes2D_cState[NCj]; FluxW = new NavierStokes2D_cState[NCj];
  // Set the solution residuals, gradients, limiters, and other values to zero.
  for (int j = JCl-Nghost; j <= JCu+Nghost; j++) {
    for (int i = ICl-Nghost; i <= ICu+Nghost; i++) {
      for (int k = 0; k < NUMBER_OF_RESIDUAL_VECTORS_NAVIERSTOKES2D; k++)
	dUdt[i][j][k].Vacuum();
      dWdx[i][j].Vacuum(); dWdy[i][j].Vacuum();
      phi[i][j].Vacuum(); Uo[i][j].Vacuum();
      dt[i][j] = ZERO;
    }
  }
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::deallocate -- Deallocate memory.        *
 **********************************************************************/
inline void NavierStokes2D_Quad_Block::deallocate(void) {
  Grid.deallocate();
  for (int i = 0; i < NCi; i++) {
    delete []W[i]; W[i] = NULL; delete []U[i]; U[i] = NULL;
    delete []dt[i]; dt[i] = NULL; 
    for (int j = 0; j < NCj; j++) { delete []dUdt[i][j]; dUdt[i][j] = NULL; }
    delete []dUdt[i]; dUdt[i] = NULL;
    delete []dWdx[i]; dWdx[i] = NULL; delete []dWdy[i]; dWdy[i] = NULL;
    delete []phi[i];  phi[i]  = NULL; delete []Uo[i]; Uo[i] = NULL;
    delete []Wall[i]; Wall[i] = NULL; 
  }
  delete []W; W = NULL; delete []U; U = NULL;
  delete []dt; dt = NULL; delete []dUdt; dUdt = NULL;
  delete []dWdx; dWdx = NULL; delete []dWdy; dWdy = NULL; 
  delete []phi; phi = NULL; delete []Uo; Uo = NULL;
  delete []Wall; Wall = NULL;
  delete []FluxN; FluxN = NULL; delete []FluxS; FluxS = NULL;
  delete []FluxE; FluxE = NULL; delete []FluxW; FluxW = NULL;
  delete []WoN; WoN = NULL; delete []WoS; WoS = NULL;
  delete []WoE; WoE = NULL; delete []WoW; WoW = NULL;
  NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::Wn -- Node primitive variable solution  *
 *                                  state.                            *
 *                                                                    *
 * Zingg and Yarrow (SIAM J. Sci. Stat. Comput. Vol. 13 No. 3 1992)   *
 *                                                                    *
 **********************************************************************/
inline NavierStokes2D_pState NavierStokes2D_Quad_Block::Wn(const int ii, const int jj) {
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
         eta1, zeta1, eta2, zeta2, eta, zeta;
  NavierStokes2D_pState A, B, C, D;
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
  bb = dy*(ax-x) + bx*cy - cx*by+dx*(y-ay);
  cc = cy*(ax-x) + cx*(y-ay);
  if (fabs(aa) < TOLER*TOLER) {
    if (fabs(bb) >= TOLER*TOLER) {
      zeta1 = -cc/bb;
    } else { 
      zeta1 = -cc/sgn(bb)*(TOLER*TOLER); 
    } 
    if (fabs(cy+dy*zeta1) >= TOLER*TOLER) {
      eta1 = (y-ay-by*zeta1)/(cy+dy*zeta1); 
    } else { 
      eta1 = HALF;
    } 
    zeta2 = zeta1;
    eta2  = eta1;
  } else {
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) { 
      zeta1 = HALF*(-bb+sqrt(bb*bb-FOUR*aa*cc))/aa; 
    } else { zeta1 = -HALF*bb/aa;
    } 
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
  A = W[ii-1][jj-1]; 
  B = W[ii-1][jj  ] - W[ii-1][jj-1]; 
  C = W[ii  ][jj-1] - W[ii-1][jj-1];
  D = W[ii  ][jj  ] + W[ii-1][jj-1] - W[ii-1][jj  ] - W[ii  ][jj-1];
  return A + B*zeta + C*eta + D*zeta*eta;
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::Un -- Node conserved variable solution  *
 *                                  state.                            *
 *                                                                    *
 * Zingg and Yarrow (SIAM J. Sci. Stat. Comput. Vol. 13 No. 3 1992)   *
 *                                                                    *
 **********************************************************************/
inline NavierStokes2D_cState NavierStokes2D_Quad_Block::Un(const int ii, const int jj) {
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
         eta1, zeta1, eta2, zeta2, eta, zeta;
  NavierStokes2D_cState A, B, C, D;
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
 * NavierStokes2D_Quad_Block::Wn?? -- Cell node primitive solution    *
 *                                    states.                         *
 **********************************************************************/
inline NavierStokes2D_pState NavierStokes2D_Quad_Block::WnNW(const int ii, const int jj) {
  return Wn(ii,jj+1);
}

inline NavierStokes2D_pState NavierStokes2D_Quad_Block::WnNE(const int ii, const int jj) {
  return Wn(ii+1,jj+1);
}

inline NavierStokes2D_pState NavierStokes2D_Quad_Block::WnSE(const int ii, const int jj) {
  return Wn(ii+1,jj);
}

inline NavierStokes2D_pState NavierStokes2D_Quad_Block::WnSW(const int ii, const int jj) {
  return Wn(ii,jj);
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::Un?? -- Cell node conserved solution    *
 *                                    states.                         *
 **********************************************************************/
inline NavierStokes2D_cState NavierStokes2D_Quad_Block::UnNW(const int ii, const int jj) {
  return Un(ii,jj+1);
}

inline NavierStokes2D_cState NavierStokes2D_Quad_Block::UnNE(const int ii, const int jj) {
  return Un(ii+1,jj+1);
}

inline NavierStokes2D_cState NavierStokes2D_Quad_Block::UnSE(const int ii, const int jj) {
  return Un(ii+1,jj);
}

inline NavierStokes2D_cState NavierStokes2D_Quad_Block::UnSW(const int ii, const int jj) {
  return Un(ii,jj);
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::Ww -- Node primitive variable solution  *
 *                                  state.                            *
 *                                                                    *
 * Holmes and Connell (AIAA Paper 1989-1932-CP)                       *
 *                                                                    *
 **********************************************************************/
inline NavierStokes2D_pState NavierStokes2D_Quad_Block::Ww(const int ii, const int jj) {
  double w1, w2, w3, w4;
  Vector2D lambda, R, X0, X1, X2, X3, X4;
  Tensor2D I;
  NavierStokes2D_pState W1, W2, W3, W4;
  // Summarize cell-centres and states.
  X0 = Grid.Node[ii][jj].X;
  X1 = Grid.Cell[ii-1][jj-1].Xc; W1 = W[ii-1][jj-1];
  X2 = Grid.Cell[ii  ][jj-1].Xc; W2 = W[ii  ][jj-1];
  X3 = Grid.Cell[ii-1][jj  ].Xc; W3 = W[ii-1][jj  ];
  X4 = Grid.Cell[ii  ][jj  ].Xc; W4 = W[ii  ][jj  ];
  // Determine weighting coefficients:
  R = X1 + X2 + X3 + X4 - FOUR*X0;
  I.xx = sqr(X1.x - X0.x) + sqr(X2.x - X0.x) + sqr(X3.x - X0.x) + sqr(X4.x - X0.x);
  I.xy = (X1.x - X0.x)*(X1.y - X0.y) + (X2.x - X0.x)*(X2.y - X0.y) +
         (X3.x - X0.x)*(X3.y - X0.y) + (X4.x - X0.x)*(X4.y - X0.y);
  I.yy = sqr(X1.y - X0.y) + sqr(X2.y - X0.y) + sqr(X3.y - X0.y) + sqr(X4.y - X0.y);
  lambda.x = (I.xy*R.y - I.yy*R.x)/(I.xx*I.yy - I.xy*I.xy);
  lambda.y = (I.xy*R.x - I.xx*R.y)/(I.xx*I.yy - I.xy*I.xy);
  // Determine the weights:
  w1 = ONE + lambda*(X1 - X0);
  w2 = ONE + lambda*(X2 - X0);
  w3 = ONE + lambda*(X3 - X0);
  w4 = ONE + lambda*(X4 - X0);
  // Determine the interpolated state:
  return (w1*W1 + w2*W2 + w3*W3+ w4*W4)/(w1 + w2 + w3 + w4);
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::Uw -- Node conserved variable solution  *
 *                                  state.                            *
 *                                                                    *
 * Holmes and Connell (AIAA Paper 1989-1932-CP)                       *
 *                                                                    *
 **********************************************************************/
inline NavierStokes2D_cState NavierStokes2D_Quad_Block::Uw(const int ii, const int jj) {
  double w1, w2, w3, w4;
  Vector2D lambda, R, X0, X1, X2, X3, X4;
  Tensor2D I;
  NavierStokes2D_cState U1, U2, U3, U4;
  // Summarize cell-centres and state.
  X0 = Grid.Node[ii][jj].X;
  X1 = Grid.Cell[ii-1][jj-1].Xc; U1 = U[ii-1][jj-1];
  X2 = Grid.Cell[ii  ][jj-1].Xc; U2 = U[ii  ][jj-1];
  X3 = Grid.Cell[ii-1][jj  ].Xc; U3 = U[ii-1][jj  ];
  X4 = Grid.Cell[ii  ][jj  ].Xc; U4 = U[ii  ][jj  ];
  // Determine weighting coefficients:
  R = X1 + X2 + X3 + X4 - FOUR*X0;
  I.xx = sqr(X1.x - X0.x) + sqr(X2.x - X0.x) + sqr(X3.x - X0.x) + sqr(X4.x - X0.x);
  I.xy = (X1.x - X0.x)*(X1.y - X0.y) + (X2.x - X0.x)*(X2.y - X0.y) +
         (X3.x - X0.x)*(X3.y - X0.y) + (X4.x - X0.x)*(X4.y - X0.y);
  I.yy = sqr(X1.y - X0.y) + sqr(X2.y - X0.y) + sqr(X3.y - X0.y) + sqr(X4.y - X0.y);
  lambda.x = (I.xy*R.y - I.yy*R.x)/(I.xx*I.yy - I.xy*I.xy);
  lambda.y = (I.xy*R.x - I.xx*R.y)/(I.xx*I.yy - I.xy*I.xy);
  // Determine the weights:
  w1 = ONE + lambda*(X1 - X0);
  w2 = ONE + lambda*(X2 - X0);
  w3 = ONE + lambda*(X3 - X0);
  w4 = ONE + lambda*(X4 - X0);
  // Determine the interpolated state:
  return (w1*U1 + w2*U2 + w3*U3+ w4*U4)/(w1 + w2 + w3 + w4);
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::Ww?? -- Cell node primitive solution    *
 *                                    states.                         *
 **********************************************************************/
inline NavierStokes2D_pState NavierStokes2D_Quad_Block::WwNW(const int ii, const int jj) {
  return Ww(ii,jj+1);
}

inline NavierStokes2D_pState NavierStokes2D_Quad_Block::WwNE(const int ii, const int jj) {
  return Ww(ii+1,jj+1);
}

inline NavierStokes2D_pState NavierStokes2D_Quad_Block::WwSE(const int ii, const int jj) {
  return Ww(ii+1,jj);
}

inline NavierStokes2D_pState NavierStokes2D_Quad_Block::WwSW(const int ii, const int jj) {
  return Ww(ii,jj);
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::Uw?? -- Cell node conserved solution    *
 *                                    states.                         *
 **********************************************************************/
inline NavierStokes2D_cState NavierStokes2D_Quad_Block::UwNW(const int ii, const int jj) {
  return Uw(ii,jj+1);
}

inline NavierStokes2D_cState NavierStokes2D_Quad_Block::UwNE(const int ii, const int jj) {
  return Uw(ii+1,jj+1);
}

inline NavierStokes2D_cState NavierStokes2D_Quad_Block::UwSE(const int ii, const int jj) {
  return Uw(ii+1,jj);
}

inline NavierStokes2D_cState NavierStokes2D_Quad_Block::UwSW(const int ii, const int jj) {
  return Uw(ii,jj);
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::evaluate_limiters -- Set flag to        *
 *                                                 evaluate limiters. *
 **********************************************************************/
inline void NavierStokes2D_Quad_Block::evaluate_limiters(void) {
  Freeze_Limiter = OFF;
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::freeze_limiters -- Set flag to freeze   *
 *                                               limiters.            *
 **********************************************************************/
inline void NavierStokes2D_Quad_Block::freeze_limiters(void) {
  Freeze_Limiter = ON; 
}

/**********************************************************************
 * NavierStokes2D_Quad_Block -- Input-output operators.               *
 **********************************************************************/
inline ostream &operator << (ostream &out_file,
			     const NavierStokes2D_Quad_Block &SolnBlk) {
  out_file << SolnBlk.Grid;
  out_file << SolnBlk.NCi << " " << SolnBlk.ICl << " " << SolnBlk.ICu << endl;
  out_file << SolnBlk.NCj << " " << SolnBlk.JCl << " " << SolnBlk.JCu << endl;
  out_file << SolnBlk.Nghost << endl;
  out_file << SolnBlk.Axisymmetric << endl;
  out_file << SolnBlk.Flow_Type << endl;
  out_file << SolnBlk.Vwall << endl;
  out_file << SolnBlk.Twall << endl;
  if (SolnBlk.NCi == 0 || SolnBlk.NCj == 0) return out_file;
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      out_file << SolnBlk.U[i][j];
      out_file << endl;
    }
  }
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    out_file << SolnBlk.WoW[j] << endl;
    out_file << SolnBlk.WoE[j] << endl;
  }
  for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
    out_file << SolnBlk.WoS[i] << endl;
    out_file << SolnBlk.WoN[i] << endl;
  }
  return out_file;
}

inline istream &operator >> (istream &in_file,
			     NavierStokes2D_Quad_Block &SolnBlk) {
  int ni, il, iu, nj, jl, ju, ng;
  Grid2D_Quad_Block New_Grid;
  in_file >> New_Grid;
  in_file.setf(ios::skipws);
  in_file >> ni >> il >> iu;
  in_file >> nj >> jl >> ju;
  in_file >> ng;
  in_file >> SolnBlk.Axisymmetric;
  in_file >> SolnBlk.Flow_Type;
  in_file.unsetf(ios::skipws);
  in_file >> SolnBlk.Vwall;
  in_file.setf(ios::skipws);
  in_file >> SolnBlk.Twall;
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
      SolnBlk.W[i][j] = SolnBlk.U[i][j].W();
      for (int k = 0; k < NUMBER_OF_RESIDUAL_VECTORS_NAVIERSTOKES2D; k++) {
	SolnBlk.dUdt[i][j][k].Vacuum();
      }
      SolnBlk.dWdx[i][j].Vacuum();
      SolnBlk.dWdy[i][j].Vacuum();
      SolnBlk.phi[i][j].Vacuum();
      SolnBlk.Uo[i][j].Vacuum();
      SolnBlk.dt[i][j] = ZERO;
    }
  }
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    in_file >> SolnBlk.WoW[j];
    in_file >> SolnBlk.WoE[j];
  }
  for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
    in_file >> SolnBlk.WoS[i];
    in_file >> SolnBlk.WoN[i];
  }
  in_file.setf(ios::skipws);
  return in_file;
}

/**********************************************************************
 *                                                                    *
 * MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING.                     *
 *                                                                    *
 **********************************************************************/

/**********************************************************************
 * NavierStokes2D_Quad_Block::NumVar -- Returns number of state       *
 *                                      variables.                    *
 **********************************************************************/
inline int NavierStokes2D_Quad_Block::NumVar(void) {
  return NUM_VAR_NAVIERSTOKES2D;
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::LoadSendBuffer -- Loads send message    *
 *                                              buffer.               *
 **********************************************************************/
inline int NavierStokes2D_Quad_Block::LoadSendBuffer(double *buffer,
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
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = U[i][j][k];
      }
    }
  }
  return 0;
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::LoadSendBuffer_F2C --                   *
 *                     Loads send message buffer for fine to coarse   *
 *                     block message passing.                         *
 **********************************************************************/
inline int NavierStokes2D_Quad_Block::LoadSendBuffer_F2C(double *buffer,
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
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
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
 * NavierStokes2D_Quad_Block::LoadSendBuffer_C2F --                   *
 *                     Loads send message buffer for coarse to fine   *
 *                     block message passing.                         *
 **********************************************************************/
inline int NavierStokes2D_Quad_Block::LoadSendBuffer_C2F(double *buffer,
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
  NavierStokes2D_pState Wfine;
  NavierStokes2D_cState Ufine;
  int Limiter = LIMITER_VENKATAKRISHNAN;

  if (j_inc > 0) {
    if (i_inc > 0) {
      for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
	  // Perform limited linear least squares reconstruction in cell (i,j_min).
	  SubcellReconstruction(i,j,Limiter);
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
	    for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
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
	    Wfine = W[i][j] + (phi[i][j]^dWdx[i][j])*dX.x +
		              (phi[i][j]^dWdy[i][j])*dX.y;
	    Ufine = Wfine.U();
	    for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
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
	    Wfine = W[i][j] + (phi[i][j]^dWdx[i][j])*dX.x +
	                      (phi[i][j]^dWdy[i][j])*dX.y;
	    Ufine = Wfine.U();
	    for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) { 
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
	    Wfine = W[i][j] + (phi[i][j]^dWdx[i][j])*dX.x +
		              (phi[i][j]^dWdy[i][j])*dX.y;
	    Ufine = Wfine.U();
	    for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
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
	return 1;
      } else {
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
	  // Perform limited linear least squares reconstruction in cell (i,j_min).
	  SubcellReconstruction(i,j_min,Limiter);
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X + Grid.Node[i+1][j_min].X) +
		Grid.Node[i+1][j_min].X + Grid.Cell[i][j_min].Xc +
		HALF*(Grid.Node[i+1][j_min].X + Grid.Node[i+1][j_min+1].X))/FOUR -
	       Grid.Cell[i][j_min].Xc;
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i][j_min].X +
		HALF*(Grid.Node[i][j_min].X + Grid.Node[i+1][j_min].X) +
		HALF*(Grid.Node[i][j_min].X + Grid.Node[i][j_min+1].X) +
		Grid.Cell[i][j_min].Xc)/FOUR - Grid.Cell[i][j_min].Xc;
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i][j_min].Xc +
		HALF*(Grid.Node[i+1][j_min].X + Grid.Node[i+1][j_min+1].X) +
		HALF*(Grid.Node[i][j_min+1].X + Grid.Node[i+1][j_min+1].X) +
		Grid.Node[i+1][j_min+1].X)/FOUR - Grid.Cell[i][j_min].Xc;
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X + Grid.Node[i][j_min+1].X) +
		Grid.Cell[i][j_min].Xc + Grid.Node[i][j_min+1].X +
		HALF*(Grid.Node[i][j_min+1].X + Grid.Node[i+1][j_min+1].X))/FOUR -
	       Grid.Cell[i][j_min].Xc;
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) { 
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
	  SubcellReconstruction(i,j_min,Limiter);
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X + Grid.Node[i][j_min+1].X) +
		Grid.Cell[i][j_min].Xc + Grid.Node[i][j_min+1].X +
		HALF*(Grid.Node[i][j_min+1].X + Grid.Node[i+1][j_min+1].X))/FOUR -
	       Grid.Cell[i][j_min].Xc;
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) { 
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i][j_min].Xc +
		HALF*(Grid.Node[i+1][j_min].X + Grid.Node[i+1][j_min+1].X) +
		HALF*(Grid.Node[i][j_min+1].X + Grid.Node[i+1][j_min+1].X) +
		Grid.Node[i+1][j_min+1].X)/FOUR - Grid.Cell[i][j_min].Xc;
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i][j_min].X +
		HALF*(Grid.Node[i][j_min].X + Grid.Node[i+1][j_min].X) +
		HALF*(Grid.Node[i][j_min].X + Grid.Node[i][j_min+1].X) +
		Grid.Cell[i][j_min].Xc)/FOUR - Grid.Cell[i][j_min].Xc;
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X + Grid.Node[i+1][j_min].X) +
		Grid.Node[i+1][j_min].X + Grid.Cell[i][j_min].Xc +
		HALF*(Grid.Node[i+1][j_min].X + Grid.Node[i+1][j_min+1].X))/FOUR -
	       Grid.Cell[i][j_min].Xc;
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
      } else {
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
	  // Perform limited linear least squares reconstruction in cell (i,j_min).
	  SubcellReconstruction(i,j_min,Limiter);
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i][j_min].Xc +
		HALF*(Grid.Node[i+1][j_min].X + Grid.Node[i+1][j_min+1].X) +
		HALF*(Grid.Node[i][j_min+1].X + Grid.Node[i+1][j_min+1].X) +
		Grid.Node[i+1][j_min+1].X)/FOUR - Grid.Cell[i][j_min].Xc;
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X + Grid.Node[i][j_min+1].X) +
		Grid.Cell[i][j_min].Xc + Grid.Node[i][j_min+1].X +
		HALF*(Grid.Node[i][j_min+1].X + Grid.Node[i+1][j_min+1].X))/FOUR -
	       Grid.Cell[i][j_min].Xc;
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) { 
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X + Grid.Node[i+1][j_min].X) +
		Grid.Node[i+1][j_min].X + Grid.Cell[i][j_min].Xc +
		HALF*(Grid.Node[i+1][j_min].X + Grid.Node[i+1][j_min+1].X))/FOUR -
	       Grid.Cell[i][j_min].Xc;
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i][j_min].X +
		HALF*(Grid.Node[i][j_min].X + Grid.Node[i+1][j_min].X) +
		HALF*(Grid.Node[i][j_min].X + Grid.Node[i][j_min+1].X) +
		Grid.Cell[i][j_min].Xc)/FOUR - Grid.Cell[i][j_min].Xc;
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
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
	return 1;
      } else {
	for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
	  // Perform limited linear least squares reconstruction in cell (i_min,j).
	  SubcellReconstruction(i_min,j,Limiter);
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X + Grid.Node[i_min+1][j].X) +
		Grid.Node[i_min+1][j].X + Grid.Cell[i_min][j].Xc +
		HALF*(Grid.Node[i_min+1][j].X + Grid.Node[i_min+1][j+1].X))/FOUR -
	       Grid.Cell[i_min][j].Xc;
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i_min][j].X +
		HALF*(Grid.Node[i_min][j].X + Grid.Node[i_min+1][j].X) +
		HALF*(Grid.Node[i_min][j].X + Grid.Node[i_min][j+1].X) +
		Grid.Cell[i_min][j].Xc)/FOUR - Grid.Cell[i_min][j].Xc;
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i_min][j].Xc +
		HALF*(Grid.Node[i_min+1][j].X + Grid.Node[i_min+1][j+1].X) +
		HALF*(Grid.Node[i_min][j+1].X + Grid.Node[i_min+1][j+1].X) +
		Grid.Node[i_min+1][j+1].X)/FOUR - Grid.Cell[i_min][j].Xc;
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X + Grid.Node[i_min][j+1].X) +
		Grid.Cell[i_min][j].Xc + Grid.Node[i_min][j+1].X +
		HALF*(Grid.Node[i_min][j+1].X + Grid.Node[i_min+1][j+1].X))/FOUR -
	       Grid.Cell[i_min][j].Xc;
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) { 
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
	  SubcellReconstruction(i_min,j,Limiter);
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X + Grid.Node[i_min][j+1].X) +
		Grid.Cell[i_min][j].Xc + Grid.Node[i_min][j+1].X +
		HALF*(Grid.Node[i_min][j+1].X + Grid.Node[i_min+1][j+1].X))/FOUR -
	       Grid.Cell[i_min][j].Xc;
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) { 
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i_min][j].Xc +
		HALF*(Grid.Node[i_min+1][j].X + Grid.Node[i_min+1][j+1].X) +
		HALF*(Grid.Node[i_min][j+1].X + Grid.Node[i_min+1][j+1].X) +
		Grid.Node[i_min+1][j+1].X)/FOUR - Grid.Cell[i_min][j].Xc;
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i_min][j].X +
		HALF*(Grid.Node[i_min][j].X + Grid.Node[i_min+1][j].X) +
		HALF*(Grid.Node[i_min][j].X + Grid.Node[i_min][j+1].X) +
		Grid.Cell[i_min][j].Xc)/FOUR - Grid.Cell[i_min][j].Xc;
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X + Grid.Node[i_min+1][j].X) +
		Grid.Node[i_min+1][j].X + Grid.Cell[i_min][j].Xc +
		HALF*(Grid.Node[i_min+1][j].X + Grid.Node[i_min+1][j+1].X))/FOUR -
	       Grid.Cell[i_min][j].Xc;
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
      } else {
	for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
	  // Perform limited linear least squares reconstruction in cell (i_min,j).
	  SubcellReconstruction(i_min,j,Limiter);
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i_min][j].Xc +
		HALF*(Grid.Node[i_min+1][j].X + Grid.Node[i_min+1][j+1].X) +
		HALF*(Grid.Node[i_min][j+1].X + Grid.Node[i_min+1][j+1].X) +
		Grid.Node[i_min+1][j+1].X)/FOUR - Grid.Cell[i_min][j].Xc;
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X + Grid.Node[i_min][j+1].X) +
		Grid.Cell[i_min][j].Xc + Grid.Node[i_min][j+1].X +
		HALF*(Grid.Node[i_min][j+1].X + Grid.Node[i_min+1][j+1].X))/FOUR -
	       Grid.Cell[i_min][j].Xc;
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) { 
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X + Grid.Node[i_min+1][j].X) +
		Grid.Node[i_min+1][j].X + Grid.Cell[i_min][j].Xc +
		HALF*(Grid.Node[i_min+1][j].X + Grid.Node[i_min+1][j+1].X))/FOUR -
	       Grid.Cell[i_min][j].Xc;
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	  // Evaluate SW sub (fine) cell values.
	  dX = (Grid.Node[i_min][j].X +
		HALF*(Grid.Node[i_min][j].X + Grid.Node[i_min+1][j].X) +
		HALF*(Grid.Node[i_min][j].X + Grid.Node[i_min][j+1].X) +
		Grid.Cell[i_min][j].Xc)/FOUR - Grid.Cell[i_min][j].Xc;
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
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
 * NavierStokes2D_Quad_Block::UnloadReceiveBuffer -- Unloads receive  *
 *                                                   message buffer.  *
 **********************************************************************/
inline int NavierStokes2D_Quad_Block::UnloadReceiveBuffer(double *buffer,
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
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	U[i][j][k] = buffer[buffer_count];
      }
      W[i][j] = U[i][j].W();
    }
  }
  return 0;
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::UnloadReceiveBuffer_F2C --              *
 *                     Unloads receive message buffer for fine to     *
 *                     coarse block message passing.                  *
 **********************************************************************/
inline int NavierStokes2D_Quad_Block::UnloadReceiveBuffer_F2C(double *buffer,
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
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;    
	U[i][j][k] = buffer[buffer_count];
      }
      W[i][j] = U[i][j].W();
    }
  }
  return 0;
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::UnloadReceiveBuffer_C2F --              *
 *                     Unloads receive message buffer for coarse to   *
 *                     fine block message passing.                    *
 **********************************************************************/
inline int NavierStokes2D_Quad_Block::UnloadReceiveBuffer_C2F(double *buffer,
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
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	U[i][j][k] = buffer[buffer_count];
      }
      W[i][j] = U[i][j].W();
    }
  }
  return 0;
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::SubcellReconstruction --                *
 *                     Performs the subcell linear reconstruction of  *
 *                     solution state within a given cell (i,j) of    *
 *                     the computational mesh for the specified       *
 *                     quadrilateral solution block.                  *
 **********************************************************************/
inline void NavierStokes2D_Quad_Block::SubcellReconstruction(const int i,
							     const int j,
							     const int Limiter) {

  int n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4], phi_n;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  NavierStokes2D_pState DU, DUDx_ave, DUDy_ave;

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
      DU = W[i_index[n2]][j_index[n2]] - W[i][j];
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
    dWdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                 (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
    dWdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                 (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
  
    // Calculate slope limiters. 
    if (!Freeze_Limiter) {
      for (int n = 1; n <= NUM_VAR_NAVIERSTOKES2D; n++) {
	u0Min = W[i][j][n];
	u0Max = u0Min;
	for (int n2 = 0; n2 < n_pts; n2++) {
	  u0Min = min(u0Min,W[i_index[n2]][j_index[n2]][n]);
	  u0Max = max(u0Max,W[i_index[n2]][j_index[n2]][n]);
	}

	dX = Grid.xfaceE(i,j) - Grid.Cell[i][j].Xc;
	uQuad[0] = W[i][j][n] + dWdx[i][j][n]*dX.x + dWdy[i][j][n]*dX.y;
	dX = Grid.xfaceW(i,j) - Grid.Cell[i][j].Xc;
	uQuad[1] = W[i][j][n] + dWdx[i][j][n]*dX.x + dWdy[i][j][n]*dX.y;
	dX = Grid.xfaceN(i,j) - Grid.Cell[i][j].Xc;
	uQuad[2] = W[i][j][n] + dWdx[i][j][n]*dX.x + dWdy[i][j][n]*dX.y;
	dX = Grid.xfaceS(i,j) - Grid.Cell[i][j].Xc;
	uQuad[3] = W[i][j][n] + dWdx[i][j][n]*dX.x + dWdy[i][j][n]*dX.y;
    
	switch(Limiter) {
	case LIMITER_ONE :
	  phi_n = ONE;
	  break;
	case LIMITER_ZERO :
	  phi_n = ZERO;
	  break;
	case LIMITER_BARTH_JESPERSEN :
	  phi_n = Limiter_BarthJespersen(uQuad,W[i][j][n],u0Min,u0Max,4);
	  break;
	case LIMITER_VENKATAKRISHNAN :
	  phi_n = Limiter_Venkatakrishnan(uQuad,W[i][j][n],u0Min,u0Max,4);
	  break;
	case LIMITER_VANLEER :
	  phi_n = Limiter_VanLeer(uQuad,W[i][j][n],u0Min,u0Max,4);
	  break;
	case LIMITER_VANALBADA :
	  phi_n = Limiter_VanAlbada(uQuad,W[i][j][n],u0Min,u0Max,4);
	  break;
	default:
	  phi_n = Limiter_BarthJespersen(uQuad,W[i][j][n],u0Min,u0Max,4);
	  break;
	}

	phi[i][j][n] = phi_n;

      }
    }
  } else {
    dWdx[i][j].Vacuum();
    dWdy[i][j].Vacuum();
    phi[i][j].Vacuum();
  }

}

/*******************************************************************************
 * NavierStokes2D_Quad_Block::LoadSendBuffer_Flux_F2C --                       *
 *                     Loads send message buffer for fine to coarse block      *
 *                     message passing of conservative solution fluxes.        *
 *******************************************************************************/
inline int NavierStokes2D_Quad_Block::LoadSendBuffer_Flux_F2C(double *buffer,
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
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = FluxS[i  ][k] + FluxS[i+1][k];
      }
    }
  } else if (j_min == j_max && j_min == JCu) {
    for (int i = i_min; ((i_inc+2)/4) ? (i < i_max):(i > i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = FluxN[i  ][k] + FluxN[i+1][k];
      }
    }
  } else if (i_min == i_max && i_min == ICl) {
    for (int j = j_min; ((j_inc+2)/4) ? (j < j_max):(j > j_max); j += j_inc) {
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = FluxW[j][k] + FluxW[j+1][k];
      }
    }
  } else if (i_min == i_max && i_min == ICu) {
    for (int j = j_min; ((j_inc+2)/4) ? (j < j_max):(j > j_max); j += j_inc) {
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = FluxE[j][k] + FluxE[j+1][k];
      }
    }
  }
  return 0;
}

/**********************************************************************
 * NavierStokes2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C --         *
 *                     Unloads receive message buffer for fine to     *
 *                     coarse block message passing of conservative   *
 *                     solution fluxes.                               *
 **********************************************************************/
inline int NavierStokes2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C(double *buffer,
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
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	FluxS[i][k] = - buffer[buffer_count] - FluxS[i][k];
      }
    }
  } else if (j_min == j_max && j_min == JCu) {
    for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	FluxN[i][k] = - buffer[buffer_count] - FluxN[i][k];
      }
    }
  } else if (i_min == i_max && i_min == ICl) {
    for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	FluxW[j][k] = - buffer[buffer_count] - FluxW[j][k];
      }
    }
  } else if (i_min == i_max && i_min == ICu) {
    for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
      for (int k = 1; k <= NUM_VAR_NAVIERSTOKES2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	FluxE[j][k] = - buffer[buffer_count] - FluxE[j][k];
      }
    }
  }
  return 0;
}

/**********************************************************************
 * NavierStokes2D_Quad_Block -- Single Block External Subroutines.    *
 **********************************************************************/

extern void Broadcast_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk);

#ifdef _MPI_VERSION
extern void Broadcast_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk,
                                     MPI::Intracomm &Communicator,
                                     const int Source_CPU);
#endif

extern void Copy_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk1,
		                NavierStokes2D_Quad_Block &SolnBlk2);

extern int Prolong_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk_Fine,
				  NavierStokes2D_Quad_Block &SolnBlk_Original,
				  const int Sector);

extern int Restrict_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk_Coarse,
				   NavierStokes2D_Quad_Block &SolnBlk_Original_SW,
				   NavierStokes2D_Quad_Block &SolnBlk_Original_SE,
				   NavierStokes2D_Quad_Block &SolnBlk_Original_NW,
				   NavierStokes2D_Quad_Block &SolnBlk_Original_NE);

extern void ICs(NavierStokes2D_Quad_Block &SolnBlk,
		NavierStokes2D_Input_Parameters &IP,
                NavierStokes2D_pState *Wo);

extern void BCs(NavierStokes2D_Quad_Block &SolnBlk,
		NavierStokes2D_Input_Parameters &IP);

extern double CFL(NavierStokes2D_Quad_Block &SolnBlk,
                  NavierStokes2D_Input_Parameters &IP);

extern void Set_Global_TimeStep(NavierStokes2D_Quad_Block &SolnBlk,
                                const double &Dt_min);

extern double L1_Norm_Residual(NavierStokes2D_Quad_Block &SolnBlk);

extern double L2_Norm_Residual(NavierStokes2D_Quad_Block &SolnBlk);

extern double Max_Norm_Residual(NavierStokes2D_Quad_Block &SolnBlk);

extern void Linear_Reconstruction_GreenGauss(NavierStokes2D_Quad_Block &SolnBlk,
                                             const int i,
                                             const int j,
					     const int Limiter);

extern void Linear_Reconstruction_GreenGauss(NavierStokes2D_Quad_Block &SolnBlk,
					     const int Limiter);

extern void Linear_Reconstruction_LeastSquares(NavierStokes2D_Quad_Block &SolnBlk,
                                               const int i,
                                               const int j,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares_2(NavierStokes2D_Quad_Block &SolnBlk,
                                                 const int i,
                                                 const int j,
					         const int Limiter);

extern void Linear_Reconstruction_LeastSquares(NavierStokes2D_Quad_Block &SolnBlk,
					       const int Limiter);

extern void Residual_Smoothing(NavierStokes2D_Quad_Block &SolnBlk,
                               const int k_residual,
			       double &epsilon,
                               const int number_of_Gauss_Seidel_iterations);

extern void Calculate_Refinement_Criteria(double *refinement_criteria,
					  NavierStokes2D_Input_Parameters &IP,
                                          int &number_refinement_criteria,
                                          NavierStokes2D_Quad_Block &SolnBlk);

extern void Fix_Refined_Block_Boundaries(NavierStokes2D_Quad_Block &SolnBlk,
                                         const int Fix_North_Boundary,
                                         const int Fix_South_Boundary,
                                         const int Fix_East_Boundary,
                                         const int Fix_West_Boundary);

extern void Unfix_Refined_Block_Boundaries(NavierStokes2D_Quad_Block &SolnBlk);

extern void Apply_Boundary_Flux_Corrections(NavierStokes2D_Quad_Block &SolnBlk,
                                            const int Number_Neighbours_North_Boundary,
                                            const int Number_Neighbours_South_Boundary,
                                            const int Number_Neighbours_East_Boundary,
                                            const int Number_Neighbours_West_Boundary);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(NavierStokes2D_Quad_Block &SolnBlk,
                                                                const int i_stage,
								NavierStokes2D_Input_Parameters &IP,
                                                                const int Number_Neighbours_North_Boundary,
                                                                const int Number_Neighbours_South_Boundary,
                                                                const int Number_Neighbours_East_Boundary,
                                                                const int Number_Neighbours_West_Boundary);

extern int dUdt_Residual_Evaluation(NavierStokes2D_Quad_Block &SolnBlk,
				    NavierStokes2D_Input_Parameters &IP);

extern int dUdt_Multistage_Explicit(NavierStokes2D_Quad_Block &SolnBlk,
   	                            const int i_stage,
                                    NavierStokes2D_Input_Parameters &IP);

extern int Update_Solution_Multistage_Explicit(NavierStokes2D_Quad_Block &SolnBlk,
   	                                       const int i_stage,
					       NavierStokes2D_Input_Parameters &IP);

/**********************************************************************
 * NavierStokes2D_Quad_Block -- Multiple Block External Subroutines.  *
 **********************************************************************/

extern NavierStokes2D_Quad_Block* Allocate(NavierStokes2D_Quad_Block *Soln_ptr,
					   NavierStokes2D_Input_Parameters &Input_Parameters);

extern NavierStokes2D_Quad_Block* Deallocate(NavierStokes2D_Quad_Block *Soln_ptr,
					     NavierStokes2D_Input_Parameters &Input_Parameters);

extern NavierStokes2D_Quad_Block* CreateInitialSolutionBlocks(Grid2D_Quad_Block **InitMeshBlk,
							      NavierStokes2D_Quad_Block *Soln_ptr,
							      NavierStokes2D_Input_Parameters &Input_Parameters,
							      QuadTreeBlock_DataStructure &QuadTree,
							      AdaptiveBlockResourceList &GlobalSolnBlockList,
							      AdaptiveBlock2D_List &LocalSolnBlockList);

extern void ICs(NavierStokes2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
                NavierStokes2D_Input_Parameters &Input_Parameters);

extern void BCs(NavierStokes2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
		NavierStokes2D_Input_Parameters &Input_Parameters);

extern double CFL(NavierStokes2D_Quad_Block *Soln_ptr,
                  AdaptiveBlock2D_List &Soln_Block_List,
		  NavierStokes2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(NavierStokes2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                const double &Dt_min);

extern double L1_Norm_Residual(NavierStokes2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double L2_Norm_Residual(NavierStokes2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double Max_Norm_Residual(NavierStokes2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List);

extern void Evaluate_Limiters(NavierStokes2D_Quad_Block *Soln_ptr,
                              AdaptiveBlock2D_List &Soln_Block_List);

extern void Freeze_Limiters(NavierStokes2D_Quad_Block *Soln_ptr,
                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Residual_Smoothing(NavierStokes2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               NavierStokes2D_Input_Parameters &Input_Parameters,
   	                       const int I_Stage);

extern void Apply_Boundary_Flux_Corrections(NavierStokes2D_Quad_Block *Soln_ptr,
                                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(NavierStokes2D_Quad_Block *Soln_ptr,
                                                                AdaptiveBlock2D_List &Soln_Block_List,
                                                                NavierStokes2D_Input_Parameters &Input_Parameters,
   	                                                        const int I_Stage);

extern int dUdt_Residual_Evaluation(NavierStokes2D_Quad_Block *Soln_ptr,
				    AdaptiveBlockResourceList &Global_Soln_Block_List,
                                    AdaptiveBlock2D_List &Soln_Block_List,
                                    NavierStokes2D_Input_Parameters &Input_Parameters);

extern int dUdt_Multistage_Explicit(NavierStokes2D_Quad_Block *Soln_ptr,
				    AdaptiveBlockResourceList &Global_Soln_Block_List,
                                    AdaptiveBlock2D_List &Local_Soln_Block_List,
                                    NavierStokes2D_Input_Parameters &Input_Parameters,
   	                            const int I_Stage);

extern int Update_Solution_Multistage_Explicit(NavierStokes2D_Quad_Block *Soln_ptr,
                                               AdaptiveBlock2D_List &Soln_Block_List,
                                               NavierStokes2D_Input_Parameters &Input_Parameters,
   	                                       const int I_Stage);

extern int Adaptive_Mesh_Refinement(NavierStokes2D_Quad_Block *Soln_ptr,
				    NavierStokes2D_Input_Parameters &Input_Parameters,
				    QuadTreeBlock_DataStructure &QuadTree,
				    AdaptiveBlockResourceList &GlobalSolnBlockList,
				    AdaptiveBlock2D_List &LocalSolnBlockList);

extern int Initial_Adaptive_Mesh_Refinement(NavierStokes2D_Quad_Block *Soln_ptr,
					    NavierStokes2D_Input_Parameters &Input_Parameters,
					    QuadTreeBlock_DataStructure &QuadTree,
					    AdaptiveBlockResourceList &GlobalSolnBlockList,
					    AdaptiveBlock2D_List &LocalSolnBlockList);

extern int Uniform_Adaptive_Mesh_Refinement(NavierStokes2D_Quad_Block *Soln_ptr,
					    NavierStokes2D_Input_Parameters &Input_Parameters,
					    QuadTreeBlock_DataStructure &QuadTree,
					    AdaptiveBlockResourceList &GlobalSolnBlockList,
					    AdaptiveBlock2D_List &LocalSolnBlockList);

extern int Boundary_Adaptive_Mesh_Refinement(NavierStokes2D_Quad_Block *Soln_ptr,
					     NavierStokes2D_Input_Parameters &Input_Parameters,
					     QuadTreeBlock_DataStructure &QuadTree,
					     AdaptiveBlockResourceList &GlobalSolnBlockList,
					     AdaptiveBlock2D_List &LocalSolnBlockList);

extern int Flat_Plate_Adaptive_Mesh_Refinement(NavierStokes2D_Quad_Block *Soln_ptr,
					       NavierStokes2D_Input_Parameters &Input_Parameters,
					       QuadTreeBlock_DataStructure &QuadTree,
					       AdaptiveBlockResourceList &GlobalSolnBlockList,
					       AdaptiveBlock2D_List &LocalSolnBlockList);

/**********************************************************************
 * NavierStokes2D_Quad_Block -- IO Single Block External Subroutines. *
 **********************************************************************/

extern void Write_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk,
	                         ostream &Out_File);

extern void Read_Solution_Block(NavierStokes2D_Quad_Block &SolnBlk,
	                        istream &In_File);

extern void Output_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
			   NavierStokes2D_Input_Parameters &IP,
		           const int Number_of_Time_Steps,
                           const double &Time,
                           const int Block_Number,
                           const int Output_Title,
	                   ostream &Out_File);

extern void Output_Cells_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_Nozzleless_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
				      const int Number_of_Time_Steps,
				      const double &Time,
				      const int Block_Number,
				      const int Output_Title,
				      ostream &Out_File,
				      const double &po,
				      const double &rhoo);

extern void Output_Nodes_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_Gradients_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
				     const int Number_of_Time_Steps,
				     const double &Time,
				     const int Block_Number,
				     const int Output_Title,
				     ostream &Out_File);

extern void Output_Quasi3D_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
				   NavierStokes2D_Input_Parameters &IP,
				   const int Number_of_Time_Steps,
				   const double &Time,
				   const int Block_Number,
				   const int Output_Title,
				   ostream &Out_File);

extern void Output_Ringleb_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
				   const int Block_Number,
				   const int Output_Title,
				   ostream &Out_File,
				   double &l1_norm,
				   double &l2_norm,
				   double &max_norm,
				   double &area);

extern void Output_Viscous_Channel_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
					   NavierStokes2D_Input_Parameters &IP,
					   const int Block_Number,
					   const int Output_Title,
					   ostream &Out_File,
					   double &l1_norm,
					   double &l2_norm,
					   double &max_norm);

extern void Output_Viscous_Pipe_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
					NavierStokes2D_Input_Parameters &IP,
					const int Block_Number,
					const int Output_Title,
					ostream &Out_File,
					double &l1_norm,
					double &l2_norm,
					double &max_norm);

extern void Output_Turbulent_Pipe_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
					  const int Block_Number,
					  const int Output_Title,
					  const int Output_Data,
					  ostream &Out_File,
					  const double &Re,
					  const double &Pipe_Radius,
					  const int &variable_flag);

extern void Output_Flat_Plate_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
				      const int Block_Number,
				      const int Output_Title_Soln,
				      ostream &Out_File_Soln,
				      const int Output_Title_Skin,
				      ostream &Out_File_Skin,
				      const NavierStokes2D_pState &Winf,
				      const double &plate_length,
				      double &l1_norm,
				      double &l2_norm,
				      double &max_norm,
				      double &area,
				      int &numberofactivecells,
				      double &l1_norm_cf,
				      double &l2_norm_cf,
				      double &max_norm_cf,
				      double &area_cf,
				      int &numberofactivecells_cf);

extern void Output_Driven_Cavity_Flow_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
					      const int Block_Number,
					      const int Output_Title,
					      ostream &Out_File_u,
					      ostream &Out_File_v,
					      const double &Re,
					      const Vector2D &Vwall,
					      const double &length);

extern void Output_Backward_Facing_Step_Tecplot(NavierStokes2D_Quad_Block &SolnBlk,
						const int Block_Number,
						const int Output_Title,
						ostream &Out_File,
						const double &step_height,
						const double &top_wall_deflection);

/**********************************************************************
 * NavierStokes2D_Quad_Block -- IO Multiple Block External Subroutines.*
 **********************************************************************/

extern int Read_Restart_Solution(NavierStokes2D_Quad_Block *Soln_ptr,
                                 AdaptiveBlock2D_List &Soln_Block_List,
                                 NavierStokes2D_Input_Parameters &IP,
		                 int &Number_of_Time_Steps,
                                 double &Time,
                                 CPUTime &CPU_Time);

extern int Write_Restart_Solution(NavierStokes2D_Quad_Block *Soln_ptr,
                                  AdaptiveBlock2D_List &Soln_Block_List,
                                  NavierStokes2D_Input_Parameters &IP,
		                  const int Number_of_Time_Steps,
                                  const double &Time,
                                  const CPUTime &CPU_Time);

extern int Output_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          NavierStokes2D_Input_Parameters &IP,
		          const int Number_of_Time_Steps,
                          const double &Time);

extern int Output_Cells_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                NavierStokes2D_Input_Parameters &IP,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Nodes_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                NavierStokes2D_Input_Parameters &IP,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Gradients_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
				    AdaptiveBlock2D_List &Soln_Block_List,
				    NavierStokes2D_Input_Parameters &IP,
				    const int Number_of_Time_Steps,
				    const double &Time);

extern int Output_Quasi3D_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  NavierStokes2D_Input_Parameters &IP,
				  const int Number_of_Time_Steps,
				  const double &Time);

extern int Output_Mesh_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               NavierStokes2D_Input_Parameters &IP,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Mesh_Gnuplot(NavierStokes2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               NavierStokes2D_Input_Parameters &IP,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Ringleb_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  NavierStokes2D_Input_Parameters &IP);

extern int Output_Viscous_Channel_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
					  AdaptiveBlock2D_List &Soln_Block_List,
					  NavierStokes2D_Input_Parameters &IP);

extern int Output_Viscous_Pipe_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
				       AdaptiveBlock2D_List &Soln_Block_List,
				       NavierStokes2D_Input_Parameters &IP);

extern int Output_Turbulent_Pipe_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
					 AdaptiveBlock2D_List &Soln_Block_List,
					 NavierStokes2D_Input_Parameters &IP);

extern int Output_Flat_Plate_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
				     AdaptiveBlock2D_List &Soln_Block_List,
				     NavierStokes2D_Input_Parameters &IP);

extern int Output_Driven_Cavity_Flow_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
					     AdaptiveBlock2D_List &Soln_Block_List,
					     NavierStokes2D_Input_Parameters &IP);

extern int Output_Backward_Facing_Step_Tecplot(NavierStokes2D_Quad_Block *Soln_ptr,
					       AdaptiveBlock2D_List &Soln_Block_List,
					       NavierStokes2D_Input_Parameters &IP);

/**********************************************************************
 * NavierStokes2D_Quad_Block -- Grid Multiple Block External          *
 *                              Subroutines.                          *
 **********************************************************************/

extern Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
					    NavierStokes2D_Input_Parameters &IP);

extern Grid2D_Quad_Block** Broadcast_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
						      NavierStokes2D_Input_Parameters &IP);

extern int Write_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                             NavierStokes2D_Input_Parameters &IP);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
							    NavierStokes2D_Input_Parameters &IP);

extern int Write_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                  NavierStokes2D_Input_Parameters &IP);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
						 NavierStokes2D_Input_Parameters &IP);

extern int Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                          NavierStokes2D_Input_Parameters &IP);

extern int Output_Nodes_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                NavierStokes2D_Input_Parameters &IP);

extern int Output_Cells_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                NavierStokes2D_Input_Parameters &IP);

/**********************************************************************
 * NavierStokes2D_Quad_Block -- Turbulence Single Block External      *
 *                              Subroutines.                          *
 **********************************************************************/

extern int Turbulence_ICs(NavierStokes2D_Quad_Block &SolnBlk,
			  NavierStokes2D_Input_Parameters &IP,
			  NavierStokes2D_pState *Wo);

extern int Turbulence_BCs(NavierStokes2D_Quad_Block &SolnBlk,
			  NavierStokes2D_Input_Parameters &IP);

extern int Apply_Turbulence_BCs(NavierStokes2D_Quad_Block &SolnBlk,
				NavierStokes2D_Input_Parameters &IP,
				const int &Turbulent_BCtype,
				const int &i, const int &j);

extern int Turbulence_Zero_Residual(NavierStokes2D_Quad_Block &SolnBlk,
				    const int i_stage,
				    NavierStokes2D_Input_Parameters &IP);

extern int Zero_Turbulence_Residuals(NavierStokes2D_Quad_Block &SolnBlk,
				     NavierStokes2D_Input_Parameters &IP,
				     const int &Turbulent_BCtype,
				     const int &i, const int &j,
				     const int &k_residual);

/**********************************************************************
 * NavierStokes2D_Quad_Block -- Turbulence Multiple Block External    *
 *                              Subroutines.                          *
 **********************************************************************/

extern int Turbulence_BCs(NavierStokes2D_Quad_Block *Soln_ptr,
			  AdaptiveBlock2D_List &Soln_Block_List,
			  NavierStokes2D_Input_Parameters &Input_Parameters);

/**********************************************************************
 * NavierStokes2D_Quad_Block -- Solvers.                              *
 **********************************************************************/

extern int NavierStokes2DQuadSolver(char *Input_File_Name_ptr,
				    int batch_flag);

#endif // _NAVIERSTOKES2D_QUAD_INCLUDED
