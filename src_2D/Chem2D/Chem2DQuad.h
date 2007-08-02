/********************** Chem2DQuad.h **********************************
  Header file defining 2D Thermally Perfect Multiple Species Navier-Stokes
  Solution State Classes.

  - based on Euler2DQuad.h 
************************************************************************/

#ifndef _CHEM2D_QUAD_INCLUDED
#define _CHEM2D_QUAD_INCLUDED

/* Include 2D Chem state, 2D cell, and 2D quadrilateral block 
   grid, quadtree, and 2D Chem input header files. */
#ifndef _CELL2D_INCLUDED
#include "../Grid/Cell2D.h"
#endif // _CELL2D_INCLUDED

#ifndef _GRID2D_INCLUDED
#include "../Grid/Grid2DQuad.h"
#endif // _GRID2D_INCLUDED

#ifndef _QUADTREE_INCLUDED
#include "../AMR/QuadTree.h"
#endif // _QUADTREE_INCLUDED

#ifndef _AMR_INCLUDED
#include "../AMR/AMR.h"
#endif // _AMR_INCLUDED

/* Also include linear systems header files. */

#ifndef _LINEARSYSTEMS_INCLUDED
#include "../Math/LinearSystems.h"
#endif // _LINEARSYSTEMS_INCLUDED

//--------- Modified for Chem ------------------
// Added and Header file for Chem

#ifndef _CHEM2D_STATE_INCLUDED
#include "Chem2DState.h"
#endif // _CHEM2D_STATE_INCLUDED

#ifndef _CHEM2D_INPUT_INCLUDED
#include "Chem2DInput.h"
#endif // _CHEM2D_INPUT_INCLUDED

// Include the turbulent wall data header file.

#ifndef _TURBULENT2D_WALLDATA_INCLUDED
#include "../Turbulent2D/Turbulent2DWallData.h"
#endif // _TURBULENT2D_WALLDATA_INCLUDED

/* Define the structures and classes. */

//I think this is OK
#define	NUMBER_OF_RESIDUAL_VECTORS_CHEM2D    3

/*!
 * Class: Chem2D_Quad_Block
 *
 * @brief Class definition of the 2D Chemistry solution blocks.
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
 *  Viscous     -- Return Viscous(=1) or Inviscid(=0)
 *  Gravity     -- Return Gravity Flag (1 ON, 0 OFF)
 * debug_level  -- Level of Verboseness (0 none, 1,2..) 
 *  Moving_wall_veloicty -- Moving wall boundary
 *                          condition velocity           
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
 * SubcellReconstruction -- Performs subcell solution reconstruction
 *                          used in adaptive mesh refinement.
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
class Chem2D_Quad_Block{
  private:
  public:
  //@{ @name Solution state arrays:
  Chem2D_pState            **W; //!< Primitive solution state.
  Chem2D_cState            **U; //!< Conserved solution state.
  //@}

  //@{ @name Grid block information:
  int                      NCi, //!< Total number of i-direction cells.
                           ICl, //!< First i-direction non-ghost cell counter.
                           ICu; //!< Final i-direction non-ghost cell counter.
  int                      NCj, //!< Total number of j-direction cells.
                           JCl, //!< First j-direction non-ghost cell counter.
                           JCu; //!< Final j-direction non-ghost cell counter.
  int                   Nghost; //!< Number of ghost cells.
  Grid2D_Quad_Block       Grid; //!< 2D quadrilateral grid geometry.
  //@}

  //@{ @name Residual and time-stepping arrays:
  double                  **dt; //!< Local time step.
  Chem2D_cState        ***dUdt; //!< Solution residual.
  Chem2D_cState           **Uo; //!< Initial solution state.
  Chem2D_cState          **RHS; //!< Solution residual???
  Chem2D_cState        **Entry; //!< Solution residual???
  static int residual_variable; //!< Static integer that indicates which variable is used for residual calculations.
  //@}

  //@{ @name Solution gradient arrays:
  Chem2D_pState         **dWdx; //!< Unlimited solution gradient (x-direction).
  Chem2D_pState         **dWdy; //!< Unlimited solution gradient (y-direction).
  Chem2D_pState          **phi; //!< Solution slope limiter.
  Chem2D_pState   **dWdx_faceN; //!< North face unlimited solution gradient (x-direction).
  Chem2D_pState   **dWdx_faceS; //!< South face unlimited solution gradient (x-direction).
  Chem2D_pState   **dWdx_faceE; //!< East face unlimited solution gradient (x-direction).
  Chem2D_pState   **dWdx_faceW; //!< West face unlimited solution gradient (x-direction).
  Chem2D_pState   **dWdy_faceN; //!< North face unlimited solution gradient (y-direction).
  Chem2D_pState   **dWdy_faceS; //!< South face unlimited solution gradient (y-direction).
  Chem2D_pState   **dWdy_faceE; //!< East face unlimited solution gradient (y-direction).
  Chem2D_pState   **dWdy_faceW; //!< West face unlimited solution gradient (y-direction).
  double          ***d_dWdx_dW; //!< Derivative of primitive solution gradient (x-direction).
  double          ***d_dWdy_dW; //!< Derivative of primitive solution gradient (y-direction).
  //@}

  //@{ @name Boundary solution flux arrays:
  Chem2D_cState         *FluxN, //!< North boundary solution flux.
                        *FluxS, //!< South boundary solution flux.
                        *FluxE, //!< East boundary solution flux.
                        *FluxW; //!< West boundary solution flux.
  //@}

  //@{ @name Problem indicator flags:
  int             Axisymmetric; //!< Axisymmetric flow indicator.
  int                Flow_Type; //!< Flow-type flag (inviscid, laminar, or k-omega).
  int                  Gravity; //!< Gravity flag
  int           Freeze_Limiter; //!< Limiter freezing indicator.
  int              debug_level; //!< Debug level flag (0=none, 1,2,3,.. level of verboseness)
  //@}

  //@{ @name Boundary condtion reference states:
  Chem2D_pState           *WoN, //!< Boundary condition reference states for north boundary.
                          *WoS, //!< Boundary condition reference states for south boundary.
                          *WoE, //!< Boundary condition reference states for east boundary.
                          *WoW; //!< Boundary condition reference states for west boundary.
  //@}

  //@{ @name Turbulence wall data arrays:
  Turbulent2DWallData   **Wall; //!< Turbulent wall data.
  //@}

  //@{ @name Flow constants:
  double               Schmidt; //!< Schmidt Number.
  double  Moving_wall_velocity; //!< Moving wall BC velocity.
  double     Pressure_gradient; //!< Specified flow pressure gradient.
  //@}

  //@{ @name Creation, copy, and assignment constructors.
  Chem2D_Quad_Block(void) {
    NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
    W = NULL; U = NULL; dt = NULL; d_dWdx_dW = NULL; d_dWdy_dW = NULL; dUdt = NULL;
    dWdx = NULL; dWdy = NULL; dWdx_faceN = NULL; dWdy_faceN = NULL;
    dWdx_faceE = NULL; dWdy_faceE = NULL;
    dWdx_faceW = NULL; dWdy_faceW = NULL;
    dWdx_faceS = NULL; dWdy_faceS = NULL;
    phi = NULL; RHS = NULL; Entry = NULL; Uo = NULL;
    FluxN = NULL; FluxS = NULL; FluxE = NULL; FluxW = NULL;
    WoN = NULL; WoS = NULL; WoE = NULL; WoW = NULL;
    Axisymmetric = 0; Gravity = 0; Flow_Type = 0;
    Moving_wall_velocity = ZERO; Pressure_gradient = ZERO; debug_level=0; 
    // Turbulent wall data:
    Wall = NULL;
  }

  Chem2D_Quad_Block(const Chem2D_Quad_Block &Soln) {
    NCi = Soln.NCi; ICl = Soln.ICl; ICu = Soln.ICu;
    NCj = Soln.NCj; JCl = Soln.JCl; JCu = Soln.JCu; Nghost = Soln.Nghost;
    Grid = Soln.Grid; W = Soln.W; U = Soln.U;
    dt = Soln.dt; d_dWdx_dW = Soln.d_dWdx_dW; d_dWdy_dW = Soln.d_dWdy_dW;
    dUdt = Soln.dUdt; dWdx = Soln.dWdx; dWdy = Soln.dWdy;
    dWdx_faceN = Soln.dWdx_faceN; dWdy_faceN = Soln.dWdy_faceN;
    dWdx_faceE = Soln.dWdx_faceE; dWdy_faceE = Soln.dWdy_faceE;
    dWdx_faceW = Soln.dWdx_faceW; dWdy_faceW = Soln.dWdy_faceW;
    dWdx_faceS = Soln.dWdx_faceS; dWdy_faceS = Soln.dWdy_faceS;
    phi = Soln.phi; RHS = Soln.RHS; Entry = Soln.Entry; Uo = Soln.Uo;
    FluxN = Soln.FluxN; FluxS = Soln.FluxS; FluxE = Soln.FluxE; FluxW = Soln.FluxW;
    WoN = Soln.WoN; WoS = Soln.WoS; WoE = Soln.WoE; WoW = Soln.WoW;
    Axisymmetric = Soln.Axisymmetric; Gravity = Soln.Gravity;
    Schmidt = Soln.Schmidt;
    Freeze_Limiter = Soln.Freeze_Limiter; debug_level = 0; 
    Moving_wall_velocity = Soln.Moving_wall_velocity; 
    Pressure_gradient = Soln.Pressure_gradient;
    // Turbulent wall data:
    Wall = Soln.Wall;
  }

  /* Destructor. */
  // ~Chem2D_Quad_Block(void);
  // Use automatically generated destructor.
  //@}

  /* Assignment operator. */
  // Chem2D_Quad_Block operator = (const Chem2D_Quad_Block &Soln);
  // Use automatically generated assignment operator.

  //@{ @name Allocate and deallocate functions.
  //! Allocate memory for structured quadrilateral solution block.
  void allocate(const int Ni, const int Nj, const int Ng);
  //! Deallocate memory for structured quadrilateral solution block.
  void deallocate(void);
  //@}

  //@{ @name Bilinear interplation (Zingg & Yarrow).
  //! Return primitive solution state at specified node.
  Chem2D_pState Wn(const int &ii, const int &jj);

  //! Retern conserverd solution state at specified node.
  Chem2D_cState Un(const int &ii, const int &jj);

  //! Return primitive solution state at cell nodes.
  Chem2D_pState WnNW(const int &ii, const int &jj); //!< Return primitive solution state at cell NW node.
  Chem2D_pState WnNE(const int &ii, const int &jj); //!< Return primitive solution state at cell NE node.
  Chem2D_pState WnSE(const int &ii, const int &jj); //!< Return primitive solution state at cell SE node.
  Chem2D_pState WnSW(const int &ii, const int &jj); //!< Return primitive solution state at cell SW node.

  Chem2D_cState UnNW(const int &ii, const int &jj); //!< Return conserved solution state at cell NW node.
  Chem2D_cState UnNE(const int &ii, const int &jj); //!< Return conserved solution state at cell NE node.
  Chem2D_cState UnSE(const int &ii, const int &jj); //!< Return conserved solution state at cell SE node.
  Chem2D_cState UnSW(const int &ii, const int &jj); //!< Return conserved solution state at cell SW node.

  //! Bilinear interpolation coefficients for Jacobian calculations.
  int BiLinearInterpolationCoefficients(double &eta, double &zeta, const int &ii, const int &jj);
  //@}

  //@{ @name Viscous Jacobian calculation routines.
  /*****************************************************************************
     dWn_dWc is the derivative of node solution w.r.t. cell center solution
     which is used by calculating the viscous Jacobians
  ******************************************************************************/
  double dWn_dWc (const int &i, const int &j, const string &Orient);
  double Wn_mu(const int &i, const int &j);
  double Wn_dmudT(const int &i, const int &j);
  double Wn_hi(const int &i, const int &j, const int &SpeciesNum);
  double Wn_Cp(const int &i, const int &j);
  double Wn_kappa(const int &i, const int &j);
  double Wn_Rtot(const int &i, const int &j);
  double Wn_Dm(const int &i, const int &j, const int &SpeciesNum);
  double Wn_Cpi(const int &i, const int &j, const int &SpeciesNum); //dhdT

  double WnNW_mu(const int &i, const int &j);
  double WnNE_mu(const int &i, const int &j);
  double WnSE_mu(const int &i, const int &j);
  double WnSW_mu(const int &i, const int &j);

  double WnNW_dmudT(const int &i, const int &j);
  double WnNE_dmudT(const int &i, const int &j);
  double WnSE_dmudT(const int &i, const int &j);
  double WnSW_dmudT(const int &i, const int &j);
     
  double WnNW_hi(const int &i, const int &j, const int &SpeciesNum);
  double WnNE_hi(const int &i, const int &j, const int &SpeciesNum);
  double WnSE_hi(const int &i, const int &j, const int &SpeciesNum);
  double WnSW_hi(const int &i, const int &j, const int &SpeciesNum);
    
  double WnNW_Dm(const int &i, const int &j, const int &SpeciesNum);
  double WnNE_Dm(const int &i, const int &j, const int &SpeciesNum);
  double WnSE_Dm(const int &i, const int &j, const int &SpeciesNum);
  double WnSW_Dm(const int &i, const int &j, const int &SpeciesNum);
  
  double WnNW_Cpi(const int &i, const int &j, const int &SpeciesNum);
  double WnNE_Cpi(const int &i, const int &j, const int &SpeciesNum);
  double WnSE_Cpi(const int &i, const int &j, const int &SpeciesNum);
  double WnSW_Cpi(const int &i, const int &j, const int &SpeciesNum);
  
  double WnNW_Cp(const int &i, const int &j);
  double WnNE_Cp(const int &i, const int &j);
  double WnSE_Cp(const int &i, const int &j);
  double WnSW_Cp(const int &i, const int &j);
 
  double WnNW_kappa(const int &i, const int &j);
  double WnNE_kappa(const int &i, const int &j);
  double WnSE_kappa(const int &i, const int &j);
  double WnSW_kappa(const int &i, const int &j);

  double WnNW_Rtot(const int &i, const int &j);
  double WnNE_Rtot(const int &i, const int &j);
  double WnSE_Rtot(const int &i, const int &j);
  double WnSW_Rtot(const int &i, const int &j);
  //@}

  //@{ @name Member functions for limiter freezing.
  void evaluate_limiters(void); //!< Set flags for limiter evaluation.
  void freeze_limiters(void);   //!< Set flags for limiter freezing.
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const Chem2D_Quad_Block &Soln);
  friend istream &operator >> (istream &in_file, Chem2D_Quad_Block &Soln);
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

/**************************************************************************
 * Chem2D_Quad_Block::allocate -- Allocate memory.                        *
 **************************************************************************/
inline void Chem2D_Quad_Block::allocate(const int Ni, const int Nj, const int Ng) {
  assert(Ni > 1 && Nj > 1); Grid.allocate(Ni, Nj, Ng);
   NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1;
   NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1; Nghost = Ng;
   W = new Chem2D_pState*[NCi]; U = new Chem2D_cState*[NCi];
   dt = new double*[NCi];
   // 2D center 0, North 1, South 2, West 3, East 4  for the d_dWdx_dW and d_dWdx_dW;
   dUdt = new Chem2D_cState**[NCi];
   dWdx = new Chem2D_pState*[NCi]; dWdy = new Chem2D_pState*[NCi];
   dWdx_faceN = new Chem2D_pState*[NCi]; dWdy_faceN = new Chem2D_pState*[NCi];
   dWdx_faceE = new Chem2D_pState*[NCi]; dWdy_faceE = new Chem2D_pState*[NCi];
   dWdx_faceW = new Chem2D_pState*[NCi]; dWdy_faceW = new Chem2D_pState*[NCi];
   dWdx_faceS = new Chem2D_pState*[NCi]; dWdy_faceS = new Chem2D_pState*[NCi];
   phi = new Chem2D_pState*[NCi]; Uo = new Chem2D_cState*[NCi];
   RHS = new Chem2D_cState*[NCi];
   Entry = new Chem2D_cState*[NCi];
   d_dWdx_dW = new double  **[NCi];
   d_dWdy_dW = new double  **[NCi];
   Wall = new Turbulent2DWallData*[NCi];

   for (int i = 0; i <= NCi-1 ; ++i ) {
      W[i] = new Chem2D_pState[NCj]; U[i] = new Chem2D_cState[NCj];
      dt[i] = new double[NCj]; dUdt[i] = new Chem2D_cState*[NCj];
      d_dWdx_dW[i] = new double  *[NCj];
      d_dWdy_dW[i] = new double  *[NCj];
      for (int j = 0; j <= NCj-1 ; ++j ){
	dUdt[i][j] = new Chem2D_cState[NUMBER_OF_RESIDUAL_VECTORS_CHEM2D];
	d_dWdx_dW[i][j] = new double [5];
	d_dWdy_dW[i][j] = new double [5];
      }
      dWdx[i] = new Chem2D_pState[NCj]; dWdy[i] = new Chem2D_pState[NCj];
      dWdx_faceN[i] = new Chem2D_pState[NCj]; dWdy_faceN[i] = new Chem2D_pState[NCj];
      dWdx_faceE[i] = new Chem2D_pState[NCj]; dWdy_faceE[i] = new Chem2D_pState[NCj];
      dWdx_faceW[i] = new Chem2D_pState[NCj]; dWdy_faceW[i] = new Chem2D_pState[NCj];
      dWdx_faceS[i] = new Chem2D_pState[NCj]; dWdy_faceS[i] = new Chem2D_pState[NCj];
      phi[i] = new Chem2D_pState[NCj]; Uo[i] = new Chem2D_cState[NCj];
      RHS[i] = new Chem2D_cState[NCj];
      Entry[i] = new Chem2D_cState[NCj];
      Wall[i] = new Turbulent2DWallData[NCj];
   } /* endfor */
   FluxN = new Chem2D_cState[NCi]; FluxS = new Chem2D_cState[NCi];
   FluxE = new Chem2D_cState[NCj]; FluxW = new Chem2D_cState[NCj];
   WoN = new Chem2D_pState[NCi]; WoS = new Chem2D_pState[NCi];
   WoE = new Chem2D_pState[NCj]; WoW = new Chem2D_pState[NCj];
   // Set the solution residuals, gradients, limiters, and other values to zero.
   for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
      for ( int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
          for ( int k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_CHEM2D-1 ; ++k ) {
	     dUdt[i][j][k].Vacuum();
          } /* endfor */
	  dWdx[i][j].Vacuum() ; dWdy[i][j].Vacuum();
	  phi[i][j].Vacuum(); Uo[i][j].Vacuum();
	  dt[i][j] = ZERO;
      } /* endfor */
   } /* endfor */

}

/**************************************************************************
 * Chem2D_Quad_Block::deallocate -- Deallocate memory.                   *
 **************************************************************************/
inline void Chem2D_Quad_Block::deallocate(void) {
  int i, j; Grid.deallocate();
   for ( i = 0; i <= NCi-1 ; ++i ) {
      delete []W[i]; W[i] = NULL; delete []U[i]; U[i] = NULL;
      delete []dt[i]; dt[i] = NULL;
      for ( j = 0; j <= NCj-1 ; ++j ) { delete []dUdt[i][j]; dUdt[i][j] = NULL;
      delete []d_dWdx_dW[i][j]; d_dWdx_dW[i][j] =NULL;
      delete []d_dWdy_dW[i][j]; d_dWdy_dW[i][j] =NULL; }
      delete []dUdt[i]; dUdt[i] = NULL;
      delete []d_dWdx_dW[i];d_dWdx_dW[i] = NULL;
      delete []d_dWdy_dW[i];d_dWdy_dW[i] = NULL;
      delete []dWdx[i]; dWdx[i] = NULL; delete []dWdy[i]; dWdy[i] = NULL;
      delete []dWdx_faceN[i]; dWdx_faceN[i] = NULL; delete []dWdy_faceN[i]; dWdy_faceN[i] = NULL;
      delete []dWdx_faceE[i]; dWdx_faceE[i] = NULL; delete []dWdy_faceE[i]; dWdy_faceE[i] = NULL;
      delete []dWdx_faceW[i]; dWdx_faceW[i] = NULL; delete []dWdy_faceW[i]; dWdy_faceW[i] = NULL;
      delete []dWdx_faceS[i]; dWdx_faceS[i] = NULL; delete []dWdy_faceS[i]; dWdy_faceS[i] = NULL;
      delete []phi[i]; phi[i] = NULL; delete []Uo[i]; Uo[i] = NULL; 
      delete []RHS[i]; RHS[i] = NULL;
      delete []Entry[i]; Entry[i] = NULL;
      delete []Wall[i]; Wall[i] = NULL; 
   } /* endfor */
   delete []W; W = NULL; delete []U; U = NULL;
   delete []dt; dt = NULL; delete []dUdt; dUdt = NULL;
   delete [] d_dWdx_dW; d_dWdx_dW = NULL;
   delete [] d_dWdy_dW; d_dWdy_dW = NULL;
   delete []dWdx; dWdx = NULL; delete []dWdy; dWdy = NULL;
   delete []dWdx_faceN; dWdx_faceN = NULL; delete []dWdy_faceN; dWdy_faceN = NULL;
   delete []dWdx_faceE; dWdx_faceE = NULL; delete []dWdy_faceE; dWdy_faceE = NULL;
   delete []dWdx_faceW; dWdx_faceW = NULL; delete []dWdy_faceW; dWdy_faceW = NULL;
   delete []dWdx_faceS; dWdx_faceS = NULL; delete []dWdy_faceS; dWdy_faceS = NULL;

   delete []phi; phi = NULL; delete []Uo; Uo = NULL;
   delete []RHS; RHS = NULL;
   delete []Entry; Entry = NULL;
   delete []Wall; Wall = NULL;
   delete []FluxN; FluxN = NULL; delete []FluxS; FluxS = NULL;
   delete []FluxE; FluxE = NULL; delete []FluxW; FluxW = NULL;
   delete []WoN; WoN = NULL; delete []WoS; WoS = NULL;
   delete []WoE; WoE = NULL; delete []WoW; WoW = NULL;
   NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
}

// Coefficients used by bilinear interpolation
inline int Chem2D_Quad_Block::BiLinearInterpolationCoefficients(double &eta, double &zeta, const int &ii, const int &jj){

   double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
           eta1, zeta1, eta2, zeta2;

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

    return 0;


}


/**************************************************************************
 * Chem2D_Quad_Block::dWn_dWc -- Derivative of Node primitive solution 
                                 w.r.r Cell primitive solution *
 **************************************************************************/
inline double Chem2D_Quad_Block::dWn_dWc(const int &i, const int &j, const string &Orient) {
 
  int ii, jj;

  if (Orient == "NW"){
    ii = i; 
    jj=j+1;
  }
  if (Orient == "NE"){
    ii = i+1; 
    jj=j+1;
  }
    
  if (Orient == "SW"){
    ii = i; 
    jj=j;
  }
  if (Orient == "SE"){
    ii = i+1; 
    jj=j;
  }
 
  double eta, zeta;
  zeta = ZERO;
  zeta = ZERO;

  BiLinearInterpolationCoefficients(eta, zeta, ii, jj);
     
  if (Orient == "NW"){
    return (eta - zeta*eta);
  }
  if (Orient == "NE"){
    return (ONE - zeta - eta + zeta*eta);
  }
  if (Orient == "SW"){
    return (zeta*eta);
  }
  if (Orient == "SE"){
    return (zeta - zeta*eta);
  }
   
}//
//Gas property
inline double Chem2D_Quad_Block::Wn_mu(const int &ii, const int &jj) {

  double eta, zeta;
  zeta = ZERO;
  zeta = ZERO;

  BiLinearInterpolationCoefficients(eta, zeta, ii, jj);
     
    return (W[ii-1][jj-1].mu() +(W[ii-1][jj].mu()-W[ii-1][jj-1].mu())*zeta+  
    	    (W[ii][jj-1].mu()-W[ii-1][jj-1].mu())*eta + 
    	    (W[ii][jj].mu()+W[ii-1][jj-1].mu()-W[ii-1][jj].mu()-W[ii][jj-1].mu())*zeta*eta);

   
}
//dmudT ... used by viscous Jacobians
inline double Chem2D_Quad_Block::Wn_dmudT(const int &ii, const int &jj) {

  double eta, zeta;
  zeta = ZERO;
  zeta = ZERO;

  BiLinearInterpolationCoefficients(eta, zeta, ii, jj);
     
    return (W[ii-1][jj-1].dmudT() +(W[ii-1][jj].dmudT()-W[ii-1][jj-1].dmudT())*zeta+  
    	    (W[ii][jj-1].dmudT()-W[ii-1][jj-1].dmudT())*eta + 
    	    (W[ii][jj].dmudT()+W[ii-1][jj-1].dmudT()-W[ii-1][jj].dmudT()-W[ii][jj-1].dmudT())*zeta*eta);

   
}
//Cp  -- specific heat
inline double Chem2D_Quad_Block::Wn_Cp(const int &ii, const int &jj) {
  double eta, zeta;
  zeta = ZERO;
  zeta = ZERO;

  BiLinearInterpolationCoefficients(eta, zeta, ii, jj);
     
    return (W[ii-1][jj-1].Cp() +(W[ii-1][jj].Cp()-W[ii-1][jj-1].Cp())*zeta+  
    	    (W[ii][jj-1].Cp()-W[ii-1][jj-1].Cp())*eta + 
    	    (W[ii][jj].Cp()+W[ii-1][jj-1].Cp()-W[ii-1][jj].Cp()-W[ii][jj-1].Cp())*zeta*eta);

   
}

//enthapy 
inline double Chem2D_Quad_Block::Wn_hi(const int &ii, const int &jj, const int &SpeciesNum) {
  double eta, zeta;
    zeta = ZERO;
    zeta = ZERO;
    
    BiLinearInterpolationCoefficients(eta, zeta, ii, jj);
     
    double Temp1, Temp2, Temp3, Temp4;

    Temp1 = W[ii-1][jj-1].T();
    Temp2 = W[ii-1][jj].T();
    Temp3 = W[ii][jj-1].T();
    Temp4 = W[ii][jj].T();
 
    double h1, h2, h3, h4;

    h1 = W[ii-1][jj-1].specdata[SpeciesNum].Enthalpy(Temp1)+W[ii-1][jj-1].specdata[SpeciesNum].Heatofform();
    h2 = W[ii-1][jj].specdata[SpeciesNum].Enthalpy(Temp2)+ W[ii-1][jj].specdata[SpeciesNum].Heatofform();
    h3 = W[ii][jj-1].specdata[SpeciesNum].Enthalpy(Temp3)+ W[ii][jj-1].specdata[SpeciesNum].Heatofform();
    h4 = W[ii][jj].specdata[SpeciesNum].Enthalpy(Temp4)+ W[ii][jj].specdata[SpeciesNum].Heatofform();


    return (h1 +(h2- h1)*zeta+ (h3-h1)*eta + (h4+h1-h2-h3)*zeta*eta);

   
}
//dhdT for species
inline double Chem2D_Quad_Block::Wn_Cpi(const int &ii, const int &jj, const int &SpeciesNum) {
  double eta, zeta;
    zeta = ZERO;
    zeta = ZERO;
    
    BiLinearInterpolationCoefficients(eta, zeta, ii, jj);
     
    double Temp1, Temp2, Temp3, Temp4;

    Temp1 = W[ii-1][jj-1].T();
    Temp2 = W[ii-1][jj].T();
    Temp3 = W[ii][jj-1].T();
    Temp4 = W[ii][jj].T();
 
    double cpi1, cpi2, cpi3, cpi4;

    cpi1 = W[ii-1][jj-1].specdata[SpeciesNum].Enthalpy_prime(Temp1);
    cpi2 = W[ii-1][jj].specdata[SpeciesNum].Enthalpy_prime(Temp2);
    cpi3 = W[ii][jj-1].specdata[SpeciesNum].Enthalpy_prime(Temp3);
    cpi4 = W[ii][jj].specdata[SpeciesNum].Enthalpy_prime(Temp4);


    return (cpi1 +(cpi2- cpi1)*zeta+ (cpi3-cpi1)*eta + (cpi4+cpi1-cpi2-cpi3)*zeta*eta);

   
}
inline double Chem2D_Quad_Block::Wn_Dm(const int &ii, const int &jj, const int &SpeciesNum) {
  double eta, zeta;
    zeta = ZERO;
    zeta = ZERO;
    
    BiLinearInterpolationCoefficients(eta, zeta, ii, jj);
     
    double Dm1, Dm2, Dm3, Dm4;

    Dm1 = W[ii-1][jj-1].spec[SpeciesNum].diffusion_coef;
    Dm2= W[ii-1][jj].spec[SpeciesNum].diffusion_coef;
    Dm3 = W[ii][jj-1].spec[SpeciesNum].diffusion_coef;
    Dm4 = W[ii][jj].spec[SpeciesNum].diffusion_coef;

    return (Dm1 +(Dm2- Dm1)*zeta+ (Dm3-Dm1)*eta + (Dm4+Dm1-Dm2-Dm3)*zeta*eta);

   
}
inline double Chem2D_Quad_Block::Wn_kappa(const int &ii, const int &jj) {
  double eta, zeta;
  zeta = ZERO;
  zeta = ZERO;
  
  BiLinearInterpolationCoefficients(eta, zeta, ii, jj);
  
     
  return (W[ii-1][jj-1].kappa() +(W[ii-1][jj].kappa()-W[ii-1][jj-1].kappa())*zeta+  
    	    (W[ii][jj-1].kappa()-W[ii-1][jj-1].kappa())*eta + 
    	    (W[ii][jj].kappa()+W[ii-1][jj-1].kappa()-W[ii-1][jj].kappa()-W[ii][jj-1].kappa())*zeta*eta);

   
}

inline double Chem2D_Quad_Block::Wn_Rtot(const int &ii, const int &jj) {
  double eta, zeta;
  zeta = ZERO;
  zeta = ZERO;
  
  BiLinearInterpolationCoefficients(eta, zeta, ii, jj);
  
     
  return (W[ii-1][jj-1].Rtot() +(W[ii-1][jj].Rtot()-W[ii-1][jj-1].Rtot())*zeta+  
    	    (W[ii][jj-1].Rtot()-W[ii-1][jj-1].Rtot())*eta + 
    	    (W[ii][jj].Rtot()+W[ii-1][jj-1].Rtot()-W[ii-1][jj].Rtot()-W[ii][jj-1].Rtot())*zeta*eta);

   
}

/**************************************************************************
 * Chem2D_Quad_Block::Wn??-- Get cell node gas properties                *
 **************************************************************************/
inline double Chem2D_Quad_Block::WnNW_mu(const int &ii, const int &jj) {
  return (Wn_mu(ii, jj+1));
}

inline double Chem2D_Quad_Block::WnNE_mu(const int &ii, const int &jj) {
  return (Wn_mu(ii+1, jj+1));
}

inline double Chem2D_Quad_Block::WnSE_mu(const int &ii, const int &jj) {
  return (Wn_mu(ii+1, jj));
}

inline double Chem2D_Quad_Block::WnSW_mu(const int &ii, const int &jj) {
  return (Wn_mu(ii, jj));
}

/**************************************************************************
 * Chem2D_Quad_Block::Wn??-- Get cell node dmudT                          *
 **************************************************************************/
inline double Chem2D_Quad_Block::WnNW_dmudT(const int &ii, const int &jj) {
  return (Wn_dmudT(ii, jj+1));
}

inline double Chem2D_Quad_Block::WnNE_dmudT(const int &ii, const int &jj) {
  return (Wn_dmudT(ii+1, jj+1));
}

inline double Chem2D_Quad_Block::WnSE_dmudT(const int &ii, const int &jj) {
  return (Wn_dmudT(ii+1, jj));
}

inline double Chem2D_Quad_Block::WnSW_dmudT(const int &ii, const int &jj) {
  return (Wn_dmudT(ii, jj));
}

/**************************************************************************
 * Chem2D_Quad_Block::Wn??-- Get cell node gas properties                *
 **************************************************************************/
inline double Chem2D_Quad_Block::WnNW_Cp(const int &ii, const int &jj) {
  return (Wn_Cp(ii, jj+1));
}

inline double Chem2D_Quad_Block::WnNE_Cp(const int &ii, const int &jj) {
  return (Wn_Cp(ii+1, jj+1));
}

inline double Chem2D_Quad_Block::WnSE_Cp(const int &ii, const int &jj) {
  return (Wn_Cp(ii+1, jj));
}

inline double Chem2D_Quad_Block::WnSW_Cp(const int &ii, const int &jj) {
  return (Wn_Cp(ii, jj));
}

/**************************************************************************
 * Chem2D_Quad_Block::Wn??-- Get cell node gas properties                *
 **************************************************************************/
inline double Chem2D_Quad_Block::WnNW_hi(const int &ii, const int &jj, const int &SpeciesNum) {
  return (Wn_hi(ii, jj+1, SpeciesNum));
}

inline double Chem2D_Quad_Block::WnNE_hi(const int &ii, const int &jj, const int &SpeciesNum) {
  return (Wn_hi(ii+1, jj+1, SpeciesNum));
}

inline double Chem2D_Quad_Block::WnSE_hi(const int &ii, const int &jj, const int &SpeciesNum) {
  return (Wn_hi(ii+1, jj, SpeciesNum));
}

inline double Chem2D_Quad_Block::WnSW_hi(const int &ii, const int &jj, const int &SpeciesNum) {
  return (Wn_hi(ii, jj, SpeciesNum));
}

/**************************************************************************
 * Chem2D_Quad_Block::Wn??-- Get cell node gas properties                *
 **************************************************************************/
inline double Chem2D_Quad_Block::WnNW_Dm(const int &ii, const int &jj, const int &SpeciesNum) {
  return (Wn_Dm(ii, jj+1, SpeciesNum));
}

inline double Chem2D_Quad_Block::WnNE_Dm(const int &ii, const int &jj, const int &SpeciesNum) {
  return (Wn_Dm(ii+1, jj+1, SpeciesNum));
}

inline double Chem2D_Quad_Block::WnSE_Dm(const int &ii, const int &jj, const int &SpeciesNum) {
  return (Wn_Dm(ii+1, jj, SpeciesNum));
}

inline double Chem2D_Quad_Block::WnSW_Dm(const int &ii, const int &jj, const int &SpeciesNum) {
  return (Wn_Dm(ii, jj, SpeciesNum));
}
/**************************************************************************
 * Chem2D_Quad_Block::Wn??-- Get cell node gas properties                *
 **************************************************************************/
inline double Chem2D_Quad_Block::WnNW_Cpi(const int &ii, const int &jj, const int &SpeciesNum) {
  return (Wn_Cpi(ii, jj+1, SpeciesNum));
}

inline double Chem2D_Quad_Block::WnNE_Cpi(const int &ii, const int &jj, const int &SpeciesNum) {
  return (Wn_Cpi(ii+1, jj+1, SpeciesNum));
}

inline double Chem2D_Quad_Block::WnSE_Cpi(const int &ii, const int &jj, const int &SpeciesNum) {
  return (Wn_Cpi(ii+1, jj, SpeciesNum));
}

inline double Chem2D_Quad_Block::WnSW_Cpi(const int &ii, const int &jj, const int &SpeciesNum) {
  return (Wn_Cpi(ii, jj, SpeciesNum));
}
/**************************************************************************
 * Chem2D_Quad_Block::Wn??-- Get cell node gas properties                *
 **************************************************************************/
inline double Chem2D_Quad_Block::WnNW_kappa(const int &ii, const int &jj) {
  return (Wn_kappa(ii, jj+1));
}

inline double Chem2D_Quad_Block::WnNE_kappa(const int &ii, const int &jj) {
  return (Wn_kappa(ii+1, jj+1));
}

inline double Chem2D_Quad_Block::WnSE_kappa(const int &ii, const int &jj) {
  return (Wn_kappa(ii+1, jj));
}

inline double Chem2D_Quad_Block::WnSW_kappa(const int &ii, const int &jj) {
  return (Wn_kappa(ii, jj));
}

/**************************************************************************
 * Chem2D_Quad_Block::Wn??-- Get cell node gas properties                *
 **************************************************************************/
inline double Chem2D_Quad_Block::WnNW_Rtot(const int &ii, const int &jj) {
  return (Wn_Rtot(ii, jj+1));
}

inline double Chem2D_Quad_Block::WnNE_Rtot(const int &ii, const int &jj) {
  return (Wn_Rtot(ii+1, jj+1));
}

inline double Chem2D_Quad_Block::WnSE_Rtot(const int &ii, const int &jj) {
  return (Wn_Rtot(ii+1, jj));
}

inline double Chem2D_Quad_Block::WnSW_Rtot(const int &ii, const int &jj) {
  return (Wn_Rtot(ii, jj));
}

/**************************************************************************
 * Chem2D_Quad_Block::Wn -- Node primitive solution.                     *
 **************************************************************************/
inline Chem2D_pState Chem2D_Quad_Block::Wn(const int &ii, const int &jj) {
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y,
    eta1, zeta1, eta2, zeta2, eta, zeta;

  //Chem2D_pState A, B, C, D;
       
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
    //A=W[ii-1][jj-1]; B=W[ii-1][jj]-W[ii-1][jj-1]; C=W[ii][jj-1]-W[ii-1][jj-1];
    //D=W[ii][jj]+W[ii-1][jj-1]-W[ii-1][jj]-W[ii][jj-1];
    //cout<<"\n Wn "<<ii<<" "<<jj<<" "<<W[ii-1][jj-1]<<endl; cout.flush;
  
  return (W[ii-1][jj-1] +(W[ii-1][jj]-W[ii-1][jj-1])*zeta+  
	  (W[ii][jj-1]-W[ii-1][jj-1])*eta + 
	  (W[ii][jj]+W[ii-1][jj-1]-W[ii-1][jj]-W[ii][jj-1])*zeta*eta);
  
  // (A+B*zeta+C*eta+D*zeta*eta);
}
/**************************************************************************
/**************************************************************************
 * Chem2D_Quad_Block::Wn?? -- Get cell node primitive solution states.   *
 **************************************************************************/
inline Chem2D_pState Chem2D_Quad_Block::WnNW(const int &ii, const int &jj) {
  return (Wn(ii, jj+1));
}

inline Chem2D_pState Chem2D_Quad_Block::WnNE(const int &ii, const int &jj) {
  return (Wn(ii+1, jj+1));
}

inline Chem2D_pState Chem2D_Quad_Block::WnSE(const int &ii, const int &jj) {
  return (Wn(ii+1, jj));
}

inline Chem2D_pState Chem2D_Quad_Block::WnSW(const int &ii, const int &jj) {
  return (Wn(ii, jj));
}

/**************************************************************************
 * Chem2D_Quad_Block::Un -- Node primitive solution.                     *
 **************************************************************************/
inline Chem2D_cState Chem2D_Quad_Block::Un(const int &ii, const int &jj) {
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y,
    eta1, zeta1, eta2, zeta2, eta, zeta;
    //Chem2D_cState A, B, C, D;
 
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
    //A=U[ii-1][jj-1]; B=U[ii-1][jj]-U[ii-1][jj-1]; C=U[ii][jj-1]-U[ii-1][jj-1];
    // D=U[ii][jj]+U[ii-1][jj-1]-U[ii-1][jj]-U[ii][jj-1];
  
    return (U[ii-1][jj-1] +(U[ii-1][jj]-U[ii-1][jj-1])*zeta+  
	    (U[ii][jj-1]-U[ii-1][jj-1])*eta + 
	    (U[ii][jj]+U[ii-1][jj-1]-U[ii-1][jj]-U[ii][jj-1])*zeta*eta);
    //(A+B*zeta+C*eta+D*zeta*eta);
}
/**************************************************************************
 * Chem2D_Quad_Block::Un?? -- Get cell node conserved solution states.    *
 **************************************************************************/
inline Chem2D_cState Chem2D_Quad_Block::UnNW(const int &ii, const int &jj) {
  return (Un(ii, jj+1));
}

inline Chem2D_cState Chem2D_Quad_Block::UnNE(const int &ii, const int &jj) {
  return (Un(ii+1, jj+1));
}

inline Chem2D_cState Chem2D_Quad_Block::UnSE(const int &ii, const int &jj) {
  return (Un(ii+1, jj));
}

inline Chem2D_cState Chem2D_Quad_Block::UnSW(const int &ii, const int &jj) {
  return (Un(ii, jj));
}

/**************************************************************************
 * Chem2D_Quad_Block::evaluate_limiters -- Set flag to evaluate limiters.*
 **************************************************************************/
inline void Chem2D_Quad_Block::evaluate_limiters(void) {
  Freeze_Limiter = OFF; 
}

/**************************************************************************
 * Chem2D_Quad_Block::freeze_limiters -- Set flag to freeze limiters.    *
 **************************************************************************/
inline void Chem2D_Quad_Block::freeze_limiters(void) {
  Freeze_Limiter = ON; 
}

/**************************************************************************
 * Chem2D_Quad_Block -- Input-output operators.                           *
 **************************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Chem2D_Quad_Block &SolnBlk) {
  int i, j; 
  out_file << SolnBlk.Grid;
  out_file << SolnBlk.NCi << " " << SolnBlk.ICl << " " 
	   << SolnBlk.ICu << " " << SolnBlk.Nghost << "\n";
  out_file << SolnBlk.NCj << " " << SolnBlk.JCl << " " << SolnBlk.JCu << "\n";
  out_file << SolnBlk.Axisymmetric << "\n";
  out_file << SolnBlk.Flow_Type <<"\n";
  out_file << SolnBlk.Schmidt <<"\n";
  out_file << SolnBlk.Gravity <<"\n";
  out_file << SolnBlk.debug_level <<"\n";
  out_file << SolnBlk.Moving_wall_velocity <<"\n";
  out_file << SolnBlk.Pressure_gradient <<"\n";
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
			     Chem2D_Quad_Block &SolnBlk) {
  int i, j, k, ni, il, iu, nj, jl, ju, ng;
  Chem2D_pState Chem2D_W_VACUUM(ZERO, Vector2D_ZERO, ZERO, ZERO, ZERO);
  Chem2D_cState Chem2D_U_VACUUM(ZERO, Vector2D_ZERO, ZERO, ZERO, ZERO);
  Grid2D_Quad_Block New_Grid; in_file >> New_Grid;
  in_file.setf(ios::skipws);
  in_file >> ni >> il >> iu >> ng; in_file >> nj >> jl >> ju;
  in_file >> SolnBlk.Axisymmetric;
  in_file >> SolnBlk.Flow_Type;
  in_file >> SolnBlk.Schmidt;
  in_file >> SolnBlk.Gravity; 
  in_file >> SolnBlk.debug_level;
  in_file >> SolnBlk.Moving_wall_velocity;
  in_file >> SolnBlk.Pressure_gradient;
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
         for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_CHEM2D-1 ; ++k ) {
	     SolnBlk.dUdt[i][j][k] = Chem2D_U_VACUUM;
         } /* endfor */
	 SolnBlk.dWdx[i][j] = Chem2D_W_VACUUM;
	 SolnBlk.dWdy[i][j] = Chem2D_W_VACUUM;
	 SolnBlk.dWdx_faceN[i][j] = Chem2D_W_VACUUM;
	 SolnBlk.dWdy_faceN[i][j] = Chem2D_W_VACUUM;
	 SolnBlk.dWdx_faceE[i][j] = Chem2D_W_VACUUM;
	 SolnBlk.dWdy_faceE[i][j] = Chem2D_W_VACUUM;
	 SolnBlk.dWdx_faceW[i][j] = Chem2D_W_VACUUM;
	 SolnBlk.dWdy_faceW[i][j] = Chem2D_W_VACUUM;
	 SolnBlk.dWdx_faceS[i][j] = Chem2D_W_VACUUM;
	 SolnBlk.dWdy_faceS[i][j] = Chem2D_W_VACUUM;
	 SolnBlk.phi[i][j] = Chem2D_W_VACUUM;
	 SolnBlk.Uo[i][j] = Chem2D_U_VACUUM;
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
 * Chem2D_Quad_Block::NumVar -- Returns number of state variables.            *
 *******************************************************************************/
inline int Chem2D_Quad_Block::NumVar(void) {
  int num=0;
  int NUM_VAR_CHEM2D = W[num][num].NUM_VAR_CHEM2D;    
  return (int(NUM_VAR_CHEM2D));
}

/*******************************************************************************
 * Chem2D_Quad_Block::LoadSendBuffer -- Loads send message buffer.            *
 *******************************************************************************/
inline int Chem2D_Quad_Block::LoadSendBuffer(double *buffer,
                                              int &buffer_count,
                                              const int buffer_size,
                                              const int i_min, 
                                              const int i_max,
                                              const int i_inc,
                                              const int j_min, 
                                              const int j_max,
                                              const int j_inc) {
  int NUM_VAR_CHEM2D = NumVar();
  int i, j, k;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = U[i][j][k];
        } /* endfor */
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * Chem2D_Quad_Block::LoadSendBuffer_F2C -- Loads send message buffer for     *
 *                                           fine to coarse block message      *
 *                                           passing.                          *
 *******************************************************************************/
inline int Chem2D_Quad_Block::LoadSendBuffer_F2C(double *buffer,
                                                  int &buffer_count,
                                                  const int buffer_size,
                                                  const int i_min, 
                                                  const int i_max,
                                                  const int i_inc,
                                                  const int j_min, 
                                                  const int j_max,
                                                  const int j_inc) {
  int NUM_VAR_CHEM2D = NumVar();
  int i, j, k;
  for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
 * Chem2D_Quad_Block::LoadSendBuffer_C2F -- Loads send message buffer for      *
 *                                          coarse to fine block message       *
 *                                          passing.                           *
 *******************************************************************************/
inline int Chem2D_Quad_Block::LoadSendBuffer_C2F(double *buffer,
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
  int NUM_VAR_CHEM2D = NumVar();
  int i, j, k;
  Vector2D dX;
  Chem2D_pState Wfine;
  Chem2D_cState Ufine;

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
	    for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
	    for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
	    for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) { 
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
	    for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
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

/*******************************************************************************
 * Chem2D_Quad_Block::UnloadReceiveBuffer -- Unloads receive message buffer.  *
 *******************************************************************************/
inline int Chem2D_Quad_Block::UnloadReceiveBuffer(double *buffer,
                                                   int &buffer_count,
                                                   const int buffer_size,
                                                   const int i_min, 
                                                   const int i_max,
                                                   const int i_inc,
                                                   const int j_min, 
                                                   const int j_max,
                                                   const int j_inc) { 
  int NUM_VAR_CHEM2D = NumVar();
 
  for ( int j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
    for ( int i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      //Changed from Euler2D
      for ( int k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	U[i][j][k] = buffer[buffer_count];
      }
      W[i][j] = U[i][j].W();
    } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * Chem2D_Quad_Block::UnloadReceiveBuffer_F2C -- Unloads receive message      *
 *                                                buffer for fine to coarse    *
 *                                                block message passing.       *
 *******************************************************************************/
inline int Chem2D_Quad_Block::UnloadReceiveBuffer_F2C(double *buffer,
                                                       int &buffer_count,
                                                       const int buffer_size,
                                                       const int i_min, 
                                                       const int i_max,
                                                       const int i_inc,
                                                       const int j_min, 
                                                       const int j_max,
                                                       const int j_inc) { 
  int NUM_VAR_CHEM2D = NumVar();
  int i, j;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
       //Changed from Euler2D
       for ( int k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
	 buffer_count = buffer_count + 1;
	 if (buffer_count >= buffer_size) return(1);    
	 U[i][j][k] = buffer[buffer_count];
/* 	 U[i][j][k] = buffer[buffer_count]/Grid.Cell[i][j].A; */
       }
       W[i][j] = U[i][j].W();
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * Chem2D_Quad_Block::UnloadReceiveBuffer_C2F -- Unloads receive message      *
 *                                                buffer for coarse to fine    *
 *                                                block message passing.       *
 *******************************************************************************/
inline int Chem2D_Quad_Block::UnloadReceiveBuffer_C2F(double *buffer,
                                                       int &buffer_count,
                                                       const int buffer_size,
                                                       const int i_min, 
                                                       const int i_max,
                                                       const int i_inc,
                                                       const int j_min, 
                                                       const int j_max,
                                                       const int j_inc) { 
  int NUM_VAR_CHEM2D = NumVar();
  int i, j;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
       //Changed from Euler2D
       for ( int k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
	 buffer_count = buffer_count + 1;
	 if (buffer_count >= buffer_size) return(1);    
	 U[i][j][k] = buffer[buffer_count];
       }
       W[i][j] = U[i][j].W();
     } /* endfor */
  } /* endfor */
  return(0);
}

/**************************************************************************
 * Chem2D_Quad_Block::SubcellReconstruction --                           *
 *               Performs the subcell reconstruction of solution state    *
 *               within a given cell (i,j) of the computational mesh for  *
 *               the specified quadrilateral solution block.              *
 **************************************************************************/
inline void Chem2D_Quad_Block::SubcellReconstruction(const int i, 
                                                      const int j,
                                                      const int Limiter) {

  int n, n2, n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4], phi_n;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  Chem2D_pState DU, DUDx_ave, DUDy_ave;

  Chem2D_pState Chem2D_W_VACUUM;
  Chem2D_W_VACUUM.Vacuum();
 
  int NUM_VAR_CHEM2D = NumVar(); 
 
  /* Carry out the limited solution reconstruction in
     each cell of the computational mesh. */

  if (i == ICl-2 || i == ICu+2 ||
      j == JCl-2 || j == JCu+2) {
    n_pts = 0;
  } else if ((i == ICl-1) && 
             (Grid.BCtypeW[j] != BC_NONE)) {
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
  } else if ((i == ICu+1) && 
             (Grid.BCtypeE[j] != BC_NONE)) {
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
  } else if ((j == JCl-1) && 
             (Grid.BCtypeS[i] != BC_NONE)) {
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
  } else if ((j == JCu+1) && 
             (Grid.BCtypeN[i] != BC_NONE)) {
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
      DUDx_ave = Chem2D_W_VACUUM;
      DUDy_ave = Chem2D_W_VACUUM;
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
	for ( n = 1 ; n <= NUM_VAR_CHEM2D ; ++n ) {
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
      } // end limiter if
  } else {
      dWdx[i][j] = Chem2D_W_VACUUM;
      dWdy[i][j] = Chem2D_W_VACUUM; 
      phi[i][j]  = Chem2D_W_VACUUM;
  } /* endif */
    
}


/*******************************************************************************
 * Chem2D_Quad_Block::LoadSendBuffer_Flux_F2C -- Loads send message buffer for*
 *                                                fine to coarse block message *
 *                                                passing of conservative      *
 *                                                solution fluxes.             *
 *******************************************************************************/
inline int Chem2D_Quad_Block::LoadSendBuffer_Flux_F2C(double *buffer,
                                                       int &buffer_count,
                                                       const int buffer_size,
                                                       const int i_min, 
                                                       const int i_max,
                                                       const int i_inc,
                                                       const int j_min, 
                                                       const int j_max,
                                                       const int j_inc) { 
  int NUM_VAR_CHEM2D = NumVar(); 
  int i, j, k;
  if (j_min == j_max && j_min == JCl) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxS[i  ][k]+
                                   FluxS[i+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxN[i  ][k]+
                                   FluxN[i+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
     for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
        for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxW[j][k]+
                                   FluxW[j+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
     for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
        for ( k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxE[j][k]+
                                   FluxE[j+1][k]);
        } /* endfor */
     } /* endfor */
  } /* endif */
  return(0);
}

/*******************************************************************************
 * Chem2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C -- Unloads receive message *
 *                                                buffer for fine to coarse    *
 *                                                block message passing of     *
 *                                                conservative solution fluxes.*
 *******************************************************************************/
inline int Chem2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C(double *buffer,
                                                            int &buffer_count,
                                                            const int buffer_size,
                                                            const int i_min, 
                                                            const int i_max,
                                                            const int i_inc,
                                                            const int j_min, 
                                                            const int j_max,
                                                            const int j_inc) {
  int NUM_VAR_CHEM2D = NumVar();
  if (j_min == j_max && j_min == JCl) {
    for ( int i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      for ( int k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	FluxS[i][k] = - buffer[buffer_count] - FluxS[i][k];
      }
    }
  } else if (j_min == j_max && j_min == JCu) {
    for ( int i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      for ( int k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	FluxN[i][k] = - buffer[buffer_count] - FluxN[i][k];
      }
    } 
  } else if (i_min == i_max && i_min == ICl) {
    for ( int j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
      for ( int k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	FluxW[j][k] = - buffer[buffer_count] - FluxW[j][k];
      }
    } 
  } else if (i_min == i_max && i_min == ICu) {
    for ( int j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
      for ( int k = 1 ; k <= NUM_VAR_CHEM2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	FluxE[j][k] = - buffer[buffer_count] - FluxE[j][k];
      }
    } 
  } /* endif */
  return(0);
}

/**************************************************************************
 * Chem2D_Quad_Block -- Single Block External Subroutines.                * 
 * Chem2DQuadSingleBlock.cc                                               *
 **************************************************************************/

extern void Write_Solution_Block(Chem2D_Quad_Block &SolnBlk,
	                         ostream &Out_File);

extern void Read_Solution_Block(Chem2D_Quad_Block &SolnBlk,
	                        istream &In_File);

extern void Broadcast_Solution_Block(Chem2D_Quad_Block &SolnBlk);

#ifdef _MPI_VERSION
extern void Broadcast_Solution_Block(Chem2D_Quad_Block &SolnBlk,
                                     MPI::Intracomm &Communicator, 
                                     const int Source_CPU);
#endif

extern void Copy_Solution_Block(Chem2D_Quad_Block &SolnBlk1,
		                Chem2D_Quad_Block &SolnBlk2);

extern int Prolong_Solution_Block(Chem2D_Quad_Block &SolnBlk_Fine,
				  Chem2D_Quad_Block &SolnBlk_Original,
				  const int Sector);

extern int Restrict_Solution_Block(Chem2D_Quad_Block &SolnBlk_Coarse,
				   Chem2D_Quad_Block &SolnBlk_Original_SW,
				   Chem2D_Quad_Block &SolnBlk_Original_SE,
				   Chem2D_Quad_Block &SolnBlk_Original_NW,
				   Chem2D_Quad_Block &SolnBlk_Original_NE);

extern void Output_Tecplot(Chem2D_Quad_Block &SolnBlk,
			   Chem2D_Input_Parameters &Input_Parameters,
			   const int Number_of_Time_Steps,
                           const double &Time,
                           const int Block_Number,
                           const int Output_Title,
	                   ostream &Out_File);

extern void Output_Cells_Tecplot(Chem2D_Quad_Block &SolnBlk,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_Nodes_Tecplot(Chem2D_Quad_Block &SolnBlk,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_RHS(Chem2D_Quad_Block &SolnBlk,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);


extern void ICs(Chem2D_Quad_Block &SolnBlk,
 	        const int i_ICtype,
                Chem2D_pState *Wo, Chem2D_Input_Parameters &Input_Parameters);

extern void BCs(Chem2D_Quad_Block &SolnBlk, Chem2D_Input_Parameters &IP);

extern double CFL(Chem2D_Quad_Block &SolnBlk,
                  Chem2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(Chem2D_Quad_Block &SolnBlk, 
                                const double &Dt_min);

extern double L1_Norm_Residual(Chem2D_Quad_Block &SolnBlk);

extern double L2_Norm_Residual(Chem2D_Quad_Block &SolnBlk);

extern double Max_Norm_Residual(Chem2D_Quad_Block &SolnBlk);

extern void Linear_Reconstruction_GreenGauss(Chem2D_Quad_Block &SolnBlk,
                                             const int i,
                                             const int j,
					     const int Limiter);

extern void Linear_Reconstruction_GreenGauss(Chem2D_Quad_Block &SolnBlk,
					     const int Limiter);

extern void Linear_Reconstruction_LeastSquares(Chem2D_Quad_Block &SolnBlk,
                                               const int i,
                                               const int j,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares_Whatname(Chem2D_Quad_Block &SolnBlk,
				        const int i, 
                                        const int j,
				        const int Limiter);

extern void Linear_Reconstruction_LeastSquares_2(Chem2D_Quad_Block &SolnBlk,
                                                 const int i,
                                                 const int j,
					         const int Limiter);

extern void Linear_Reconstruction_LeastSquares_Diamond(Chem2D_Quad_Block &SolnBlk,
						       const int i, 
						       const int j,
						       const int Limiter);
extern void Linear_Reconstruction_LeastSquares(Chem2D_Quad_Block &SolnBlk,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares_Diamond(Chem2D_Quad_Block &SolnBlk,
					       const int Limiter);

extern void Residual_Smoothing(Chem2D_Quad_Block &SolnBlk,
                               const int k_residual,
			       double &epsilon, 
                               const int number_of_Gauss_Seidel_iterations);

extern void Calculate_Refinement_Criteria(double *refinement_criteria,
					  Chem2D_Input_Parameters &IP,
                                          int &number_refinement_criteria,
                                          Chem2D_Quad_Block &SolnBlk);

extern void Fix_Refined_Block_Boundaries(Chem2D_Quad_Block &SolnBlk,
                                         const int Fix_North_Boundary,
                                         const int Fix_South_Boundary,
                                         const int Fix_East_Boundary,
                                         const int Fix_West_Boundary);

extern void Unfix_Refined_Block_Boundaries(Chem2D_Quad_Block &SolnBlk);

extern void Apply_Boundary_Flux_Corrections(Chem2D_Quad_Block &SolnBlk,
                                            const int Number_Neighbours_North_Boundary,
                                            const int Number_Neighbours_South_Boundary,
                                            const int Number_Neighbours_East_Boundary,
                                            const int Number_Neighbours_West_Boundary);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Chem2D_Quad_Block &SolnBlk,
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

extern int dUdt_Residual_Evaluation(Chem2D_Quad_Block &SolnBlk,
				    Chem2D_Input_Parameters &Input_Parameters);


extern int dUdt_Multistage_Explicit(Chem2D_Quad_Block &SolnBlk,
   	                            const int i_stage,
                                    Chem2D_Input_Parameters &Input_Parameters);

extern int Update_Solution_Multistage_Explicit(Chem2D_Quad_Block &SolnBlk,
   	                                       const int i_stage,
                                               Chem2D_Input_Parameters &Input_Parameters);

extern void Viscous_Calculations(Chem2D_Quad_Block &SolnBlk);

/**************************************************************************
 * Chem2D_Quad_Block -- Multiple Block External Subroutines.              *
 * Chem2DQuadMultiBlock.cc                                                *
 **************************************************************************/

extern Chem2D_Quad_Block* Allocate(Chem2D_Quad_Block *Soln_ptr,
                                    Chem2D_Input_Parameters &Input_Parameters);

extern Chem2D_Quad_Block* Deallocate(Chem2D_Quad_Block *Soln_ptr,
                                      Chem2D_Input_Parameters &Input_Parameters);

extern void ICs(Chem2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
                Chem2D_Input_Parameters &Input_Parameters);

extern int Read_Restart_Solution(Chem2D_Quad_Block *Soln_ptr,
                                 AdaptiveBlock2D_List &Soln_Block_List,
                                 Chem2D_Input_Parameters &Input_Parameters,
		                 int &Number_of_Time_Steps,
                                 double &Time,
                                 CPUTime &CPU_Time);

extern int Write_Restart_Solution(Chem2D_Quad_Block *Soln_ptr,
                                  AdaptiveBlock2D_List &Soln_Block_List,
                                  Chem2D_Input_Parameters &Input_Parameters,
		                  const int Number_of_Time_Steps,
                                  const double &Time,
                                  const CPUTime &CPU_Time);

extern int Output_Tecplot(Chem2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          Chem2D_Input_Parameters &Input_Parameters,
		          const int Number_of_Time_Steps,
                          const double &Time);

extern int Output_RHS(Chem2D_Quad_Block *Soln_ptr,
		      AdaptiveBlock2D_List &Soln_Block_List,
		      Chem2D_Input_Parameters &Input_Parameters,
		      const int Number_of_Time_Steps,
		      const double &Time);

extern int Output_PERTURB(Chem2D_Quad_Block *Soln_ptr,
			  AdaptiveBlock2D_List &Soln_Block_List,
			  Chem2D_Input_Parameters &Input_Parameters,
			  const int Number_of_Time_Steps,
			  const double &Time,
			  const CPUTime &CPU_Time);

extern int Output_Cells_Tecplot(Chem2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                Chem2D_Input_Parameters &Input_Parameters,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Nodes_Tecplot(Chem2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                Chem2D_Input_Parameters &Input_Parameters,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Mesh_Tecplot(Chem2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Chem2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Mesh_Gnuplot(Chem2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Chem2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Ringleb(Chem2D_Quad_Block *Soln_ptr,
			  AdaptiveBlock2D_List &Soln_Block_List,
			  Chem2D_Input_Parameters &IP);

extern int Output_Viscous_Channel(Chem2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  Chem2D_Input_Parameters &IP);

extern int Output_Flat_Plate(Chem2D_Quad_Block *Soln_ptr,
			     AdaptiveBlock2D_List &Soln_Block_List,
			     Chem2D_Input_Parameters &IP);

extern void BCs(Chem2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
		Chem2D_Input_Parameters &Input_Parameters);

extern double CFL(Chem2D_Quad_Block *Soln_ptr,
                  AdaptiveBlock2D_List &Soln_Block_List,
                  Chem2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(Chem2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List, 
                                const double &Dt_min);

extern double L1_Norm_Residual(Chem2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double L2_Norm_Residual(Chem2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double Max_Norm_Residual(Chem2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List);

extern void Evaluate_Limiters(Chem2D_Quad_Block *Soln_ptr,
                              AdaptiveBlock2D_List &Soln_Block_List);

extern void Freeze_Limiters(Chem2D_Quad_Block *Soln_ptr,
                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Residual_Smoothing(Chem2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Chem2D_Input_Parameters &Input_Parameters,
   	                       const int I_Stage);

extern void Apply_Boundary_Flux_Corrections(Chem2D_Quad_Block *Soln_ptr,
                                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Chem2D_Quad_Block *Soln_ptr,
                                                                AdaptiveBlock2D_List &Soln_Block_List,
                                                                Chem2D_Input_Parameters &Input_Parameters,
   	                                                        const int I_Stage);

extern int dUdt_Multistage_Explicit(Chem2D_Quad_Block *Soln_ptr,
				    AdaptiveBlockResourceList &Global_Soln_Block_List,
                                    AdaptiveBlock2D_List &Soln_Block_List,
                                    Chem2D_Input_Parameters &Input_Parameters,
   	                            const int I_Stage);

extern int Update_Solution_Multistage_Explicit(Chem2D_Quad_Block *Soln_ptr,
                                               AdaptiveBlock2D_List &Soln_Block_List,
                                               Chem2D_Input_Parameters &Input_Parameters,
   	                                       const int I_Stage);

//************Viscous Stuff************************************/
//Viscous (Turbulence) Flux
extern void Viscous_Calculations(Chem2D_Quad_Block &SolnBlk);
extern Chem2D_cState Flux_Viscous_x(Chem2D_Quad_Block &SolnBlk, const int &i, const int &j);
extern Chem2D_cState Flux_Viscous_y(Chem2D_Quad_Block &SolnBlk, const int &i, const int &j);
extern Chem2D_cState Flux_Viscous_n(Chem2D_Quad_Block &SolnBlk, const string &Orient, const int &i, const int &j);
//Viscous (Turbulence) Flux (Diamond Path)
extern Chem2D_cState Flux_Viscous_EastFace_Diamond(Chem2D_Quad_Block &SolnBlk, const int &i, const int &j);
extern Chem2D_cState Flux_Viscous_NorthFace_Diamond(Chem2D_Quad_Block &SolnBlk, const int &i, const int &j);

/* Low-Reynolds number formulation;                                  */
/* Now temperarily set 5~7 nodes closest to the wall inside y+ <2.5  */
/*-----------------------------------------------------------------  */
extern void Low_ReynoldsNumber_Formulation(Chem2D_Quad_Block &SolnBlk, 
                                           Chem2D_Input_Parameters &Input_Parameters, 
                                           int i, int j);

extern double  Wall_Distance(Chem2D_Quad_Block &SolnBlk, Chem2D_Input_Parameters &Input_Parameters, int i, int j);
extern double Distance_to_Wall(Chem2D_Quad_Block &SolnBlk, const Vector2D X_cell);
/**************************************************************************
 * Chem2D_Quad_Block -- Multiple Block External Subroutines for Mesh.     * 
 * Chem2DQuadGrid.cc                                                      *
 **************************************************************************/

extern Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                            Chem2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Broadcast_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                      Chem2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                             Chem2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                                            Chem2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                  Chem2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                 Chem2D_Input_Parameters &Input_Parameters);

extern int Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                          Chem2D_Input_Parameters &Input_Parameters);

extern int Output_Nodes_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                Chem2D_Input_Parameters &Input_Parameters);

extern int Output_Cells_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                Chem2D_Input_Parameters &Input_Parameters);


/**************************************************************************
 * Chem2D_Quad_Block -- Solvers.                                          *
 * Chem2DQuadSolvers.cc                                                   *
 **************************************************************************/

extern int Chem2DQuadSolver(char *Input_File_Name_ptr, int batch_flag);

/**************************************************************************
 * Chem2D_Tools -- Tools and stuff                                         *
 * Chem2DTools.cc                                                          *
 **************************************************************************/
extern int Open_Time_Accurate_File(ofstream &Time_Accurate_File,
				   char *File_Name,
				   const int Append_to_Fileconst,
				   const Chem2D_pState &Soln);

extern int Close_Time_Accurate_File(ofstream &Time_Accurate_File);

extern void Output_to_Time_Accurate_File(ostream &Time_Accurate_File,
					 const double &Time,
					 const Chem2D_pState &Soln);

extern void Output_Ringleb_Solution(Chem2D_Quad_Block &SolnBlk,
				    const int Block_Number,
				    const int Output_Title,
				    ostream &Out_File);

extern void Output_Ringleb_Error(Chem2D_Quad_Block &SolnBlk,
				 double &l1_norm,
				 double &l2_norm,
				 double &max_norm,
				 int &numberofactivecells);

extern void Output_Viscous_Channel(Chem2D_Quad_Block &SolnBlk,
				   const int Block_Number,
				   const int Output_Title,
				   ostream &Out_File,
				   double &l1_norm,
				   double &l2_norm,
				   double &max_norm,
				   double &Vwall,
				   const double dp);

extern void Output_Flat_Plate(Chem2D_Quad_Block &SolnBlk,
			      const int Block_Number,
			      const int Output_Title_Soln,
			      ostream &Out_File_Soln,
			      const int Output_Title_Skin,
			      ostream &Out_File_Skin,
			      const Chem2D_pState &Winf,
			      double &l1_norm,
			      double &l2_norm,
			      double &max_norm);
			      


/*************** END CHEM2D ***************************************/



#endif /* _CHEM2D_QUAD_INCLUDED  */
