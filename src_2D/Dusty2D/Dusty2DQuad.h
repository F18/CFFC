/**********************************************************************
 * Dusty2DQuad.h: Header file defining 2D Dusty quadrilateral mesh    *
 *                solution classes.                                   *
 **********************************************************************/

#ifndef _DUSTY2D_QUAD_INCLUDED
#define _DUSTY2D_QUAD_INCLUDED

// Include 2D dusty state header file.

#ifndef _DUSTY2D_STATE_INCLUDED
#include "Dusty2DState.h"
#endif // _DUSTY2D_STATE_INCLUDED

// Include 2D dusty input header file.

#ifndef _DUSTY2D_INPUT_INCLUDED
#include "Dusty2DInput.h"
#endif // _DUSTY2D_INPUT_INCLUDED

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

// Include the NASA rotor 37 and 67 header files.

#ifndef _NASA_ROTOR37_INCLUDED
#include "../Grid/NASARotor37.h"
#endif // _NASA_ROTOR37_INCLUDED

#ifndef _NASA_ROTOR67_INCLUDED
#include "../Grid/NASARotor67.h"
#endif // _NASA_ROTOR67_INCLUDED

// Include ICEMCFD input header file.

#ifndef _ICEMCFD_INCLUDED
#include "../ICEM/ICEMCFD.h"
#endif // _ICEMCFD_INCLUDED

// Include the electrostatic quad block header file.

#ifndef _ELECTROSTATIC2D_QUAD_INCLUDED
#include "../Electrostatic2D/Electrostatic2DQuad.h"
#endif // _ELECTROSTATIC2D_QUAD_INCLUDED

// Include the turbulent wall data header file.

#ifndef _TURBULENT2D_WALLDATA_INCLUDED
#include "../Turbulent2D/Turbulent2DWallData.h"
#endif // _TURBULENT2D_WALLDATA_INCLUDED

// Define the structures and classes.

#define	NUMBER_OF_RESIDUAL_VECTORS_DUSTY2D  3

/*!
 * Class: Dusty2D_Quad_Block
 *
 * @brief Class definition of the 2D dusty solution blocks.
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
 *      S -- a 2D Dusty solution
 *
 * S = S;
 * cout << S; (output function)
 * cin  >> S; (input function)
 * \endverbatim
 */
class Dusty2D_Quad_Block{
private:
public:

  //! @name Defined public types:
  //@{
  typedef Grid2D_Quad_Block GridType;
  //@}

  //@{ @name Solution state arrays:
  Dusty2D_pState              **W; //!< Primitive solution state.
  Dusty2D_cState              **U; //!< Conserved solution state.
  static int      NUM_VAR_DUSTY2D; //!< Static integer that stores the number of variables.
  //@}

  //@{ @name Grid block information:
  int                         NCi, //!< Total number of i-direction cells.
                              ICl, //!< First i-direction non-ghost cell counter.
                              ICu; //!< Final i-direction non-ghost cell counter.
  int                         NCj, //!< Total number of j-direction cells.
                              JCl, //!< First j-direction non-ghost cell counter.
                              JCu; //!< Final j-direction non-ghost cell counter.
  int                      Nghost; //!< Number of ghost cells.
  Grid2D_Quad_Block          Grid; //!< 2D quadrilateral grid geometry.
  //@}

  //@{ @name Residual and time-stepping arrays:
  double                     **dt; //!< Local time step.
  Dusty2D_cState          ***dUdt; //!< Solution residual.
  Dusty2D_cState             **Uo; //!< Initial solution state.
  static int    residual_variable; //!< Static integer that indicates which variable is used for residual calculations.
  //@}

  //@{ @name Solution gradient arrays:
  Dusty2D_pState           **dWdx; //!< Unlimited solution gradient (x-direction).
  Dusty2D_pState           **dWdy; //!< Unlimited solution gradient (y-direction).
  Dusty2D_pState            **phi; //!< Solution slope limiter.
  //@}

  //@{ @name Boundary solution flux arrays:
  Dusty2D_cState           *FluxN, //!< North boundary solution flux.
                           *FluxS, //!< South boundary solution flux.
                           *FluxE, //!< East boundary solution flux.
                           *FluxW; //!< West boundary solution flux.
  //@}

  //@{ @name Problem indicator flags:
  int                Axisymmetric; //!< Axisymmetric flow indicator.
  int                   Flow_Type; //!< Flow-type indicator (inviscid, laminar, or k-omega).
  int                   Particles; //!< Particle-phase formulation flag (on or off).
  int            PhaseInteraction; //!< Phase interaction flag.
  int               Electrostatic; //!< Flag for electrostatic forces (on or off).
  int              Freeze_Limiter; //!< Limiter freezing indicator.
  //@}

  //@{ @name Boundary condtion reference states:
  Dusty2D_pState             *WoN, //!< Boundary condition reference states for north boundary.
                             *WoS, //!< Boundary condition reference states for south boundary.
                             *WoE, //!< Boundary condition reference states for east boundary.
                             *WoW; //!< Boundary condition reference states for west boundary.
  //@}

  //@{ @name Electrostatic force data (electric and potential fields):
  Electrostatic2DState       **E;
  //@}

  //@{ @name Turbulence wall data arrays:
  Turbulent2DWallData    **Wall; //!< Turbulent wall data.
  //@}

  //@{ @name Flow constants:
  Vector2D                 Vwall;
  double                   Twall;
  //@}

  //! Creation constructor.
  Dusty2D_Quad_Block(void) {
    // Problem flags:
    Axisymmetric     = OFF;
    Particles        = PARTICLE_PHASE_NONE;
    PhaseInteraction = OFF;
    Electrostatic    = OFF;
    Flow_Type        = FLOWTYPE_INVISCID;
    Freeze_Limiter   = OFF;
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
    // Electrostatic force data (electric and potential fields):
    E = NULL;
    // Turbulent wall data:
    Wall = NULL;
  }

  //! Copy constructor.
  Dusty2D_Quad_Block(const Dusty2D_Quad_Block &Soln) {
    // Problem flags:
    Axisymmetric     = Soln.Axisymmetric;
    Particles        = Soln.Particles;
    PhaseInteraction = Soln.PhaseInteraction;
    Electrostatic    = Soln.Electrostatic;
    Flow_Type        = Soln.Flow_Type;
    Freeze_Limiter   = Soln.Freeze_Limiter;
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
    // Electrostatic force data (electric and potential fields):
    E = Soln.E;
    // Turbulent wall data:
    Wall = Soln.Wall;
  }

  // Destructor.
  // ~Dusty2D_Quad_Block(void);
  // Use automatically generated destructor.

  // Assignment operator.
  // Dusty2D_Quad_Block operator = (const Dusty2D_Quad_Block &Soln);
  // Use automatically generated assignment operator.

  //@{ @name Allocate and deallocate functions.
  //! Allocate memory for structured quadrilateral solution block.
  void allocate(const int Ni, const int Nj, const int Ng);
  //! Deallocate memory for structured quadrilateral solution block.
  void deallocate(void);
  //! Allocate memory for the electrostatic variables.
  void allocate_electrostatic(void);
  //! Deallocate memory for the electrostatic variables.
  void deallocate_electrostatic(void);
  //@}

  //@{ @name Bilinear interplation (Zingg & Yarrow).
  //! Return primitive solution state at the specified node.
  Dusty2D_pState Wz(const int ii, const int jj);

  //! Return conserverd solution state at the specified node.
  Dusty2D_cState Uz(const int ii, const int jj);

  //! Return the electrostatic state at the specified node.
  Electrostatic2DState Ez(const int ii, const int jj);

  Dusty2D_pState WzNW(const int ii, const int jj); //!< Return primitive solution state at the NW node.
  Dusty2D_pState WzNE(const int ii, const int jj); //!< Return primitive solution state at the NE node.
  Dusty2D_pState WzSE(const int ii, const int jj); //!< Return primitive solution state at the SE node.
  Dusty2D_pState WzSW(const int ii, const int jj); //!< Return primitive solution state at the SW node.

  Dusty2D_cState UzNW(const int ii, const int jj); //!< Return conserved solution state at the NW node.
  Dusty2D_cState UzNE(const int ii, const int jj); //!< Return conserved solution state at the NE node.
  Dusty2D_cState UzSE(const int ii, const int jj); //!< Return conserved solution state at the SE node.
  Dusty2D_cState UzSW(const int ii, const int jj); //!< Return conserved solution state at the SW node.

  Electrostatic2DState EzNW(const int ii, const int jj); //!< Return electrostatic solution state at the NW node.
  Electrostatic2DState EzNE(const int ii, const int jj); //!< Return electrostatic solution state at the NE node.
  Electrostatic2DState EzSE(const int ii, const int jj); //!< Return electrostatic solution state at the SE node.
  Electrostatic2DState EzSW(const int ii, const int jj); //!< Return electrostatic solution state at the SW node.
  //@}

  //@{ @name Bilinear interplation.
  //! Return primitive solution state at the specified node.
  Dusty2D_pState Wn(const int ii, const int jj);

  //! Return primitive solution state at the specified node (Connell & Holmes).
  Dusty2D_pState Ww(const int ii, const int jj);

  //! Return conserverd solution state the at specified node.
  Dusty2D_cState Un(const int ii, const int jj);

  //! Return conserverd solution state the at specified node (Connell & Holmes).
  Dusty2D_cState Uw(const int ii, const int jj);

  //! Return the electrostatic solution state at the specified node.
  Electrostatic2DState En(const int ii, const int jj);

  Dusty2D_pState WnNW(const int ii, const int jj); //!< Return primitive solution state at the NW node.
  Dusty2D_pState WnNE(const int ii, const int jj); //!< Return primitive solution state at the NE node.
  Dusty2D_pState WnSE(const int ii, const int jj); //!< Return primitive solution state at the SE node.
  Dusty2D_pState WnSW(const int ii, const int jj); //!< Return primitive solution state at the SW node.

  Dusty2D_cState UnNW(const int ii, const int jj); //!< Return conserved solution state at the NW node.
  Dusty2D_cState UnNE(const int ii, const int jj); //!< Return conserved solution state at the NE node.
  Dusty2D_cState UnSE(const int ii, const int jj); //!< Return conserved solution state at the SE node.
  Dusty2D_cState UnSW(const int ii, const int jj); //!< Return conserved solution state at the SW node.

  Electrostatic2DState EnNW(const int ii, const int jj); //!< Return electrostatic solution state at the NW node.
  Electrostatic2DState EnNE(const int ii, const int jj); //!< Return electrostatic solution state at the NE node.
  Electrostatic2DState EnSE(const int ii, const int jj); //!< Return electrostatic solution state at the SE node.
  Electrostatic2DState EnSW(const int ii, const int jj); //!< Return electrostatic solution state at the SW node.
  //@}

  //@{ @name Member functions for limiter freezing.
  void evaluate_limiters(void); //!< Set flags for limiter evaluation.
  void freeze_limiters(void);   //!< Set flags for limiter freezing.
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const Dusty2D_Quad_Block &Soln);
  friend istream &operator >> (istream &in_file, Dusty2D_Quad_Block &Soln);
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

/**********************************************************************
 * Dusty2D_Quad_Block::allocate -- Allocate memory.                   *
 **********************************************************************/
inline void Dusty2D_Quad_Block::allocate(const int Ni, const int Nj, const int Ng) {
  assert(Ni > 1 && Nj > 1 && Ng > 1);
  Grid.allocate(Ni,Nj,Ng);
  NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1; 
  NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1; Nghost = Ng;
  W = new Dusty2D_pState*[NCi]; U = new Dusty2D_cState*[NCi];
  dt = new double*[NCi]; dUdt = new Dusty2D_cState**[NCi];
  dWdx = new Dusty2D_pState*[NCi]; dWdy = new Dusty2D_pState*[NCi];
  phi = new Dusty2D_pState*[NCi]; Uo = new Dusty2D_cState*[NCi];
  Wall = new Turbulent2DWallData*[NCi];
  for (int i = 0; i < NCi; i++) {
    W[i] = new Dusty2D_pState[NCj]; U[i] = new Dusty2D_cState[NCj];
    dt[i] = new double[NCj]; dUdt[i] = new Dusty2D_cState*[NCj];
    for (int j = 0; j < NCj; j++)
      dUdt[i][j] = new Dusty2D_cState[NUMBER_OF_RESIDUAL_VECTORS_DUSTY2D];
    dWdx[i] = new Dusty2D_pState[NCj]; dWdy[i] = new Dusty2D_pState[NCj];
    phi[i] = new Dusty2D_pState[NCj]; Uo[i] = new Dusty2D_cState[NCj];
    Wall[i] = new Turbulent2DWallData[NCj];
  }
  WoN = new Dusty2D_pState[NCi]; WoS = new Dusty2D_pState[NCi];
  WoE = new Dusty2D_pState[NCj]; WoW = new Dusty2D_pState[NCj];
  FluxN = new Dusty2D_cState[NCi]; FluxS = new Dusty2D_cState[NCi];
  FluxE = new Dusty2D_cState[NCj]; FluxW = new Dusty2D_cState[NCj];
  // Set the solution residuals, gradients, limiters, and other values to zero.
  for (int j = JCl-Nghost; j <= JCu+Nghost; j++) {
    for (int i = ICl-Nghost; i <= ICu+Nghost; i++) {
      for (int k = 0; k < NUMBER_OF_RESIDUAL_VECTORS_DUSTY2D; k++)
	dUdt[i][j][k].Vacuum();
      dWdx[i][j].Vacuum(); dWdy[i][j].Vacuum();
      phi[i][j].Vacuum(); Uo[i][j].Vacuum();
      dt[i][j] = ZERO;
    }
  }
}

/**********************************************************************
 * Dusty2D_Quad_Block::deallocate -- Deallocate memory.               *
 **********************************************************************/
inline void Dusty2D_Quad_Block::deallocate(void) {
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

  if (E != NULL) deallocate_electrostatic();
}

/**********************************************************************
 * Dusty2D_Quad_Block::allocate_electrostatic --                      *
 *                   Allocate memory for the electrostatic variables. *
 **********************************************************************/
inline void Dusty2D_Quad_Block::allocate_electrostatic(void) {
  E = new Electrostatic2DState*[NCi];
  for (int i = 0; i < NCi; i++) E[i] = new Electrostatic2DState[NCj];
}

/**********************************************************************
 * Dusty2D_Quad_Block::deallocate_electrostatic --                    *
 *                 Deallocate memory for the electrostatic variables. *
 **********************************************************************/
inline void Dusty2D_Quad_Block::deallocate_electrostatic(void) {
  for (int i = 0; i < NCi; i++) { delete []E[i]; E[i] = NULL; }
  delete []E; E = NULL;
}

/**********************************************************************
 * Dusty2D_Quad_Block::Wz -- Node primitive variable solution state.  *
 *                                                                    *
 * Zingg and Yarrow (SIAM J. Sci. Stat. Comput. Vol. 13 No. 3 1992)   *
 *                                                                    *
 **********************************************************************/
inline Dusty2D_pState Dusty2D_Quad_Block::Wz(const int ii, const int jj) {
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
         eta1, zeta1, eta2, zeta2, eta, zeta;
  Dusty2D_pState A, B, C, D;
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
 * Dusty2D_Quad_Block::Uz -- Node conserved variable solution state.  *
 *                                                                    *
 * Zingg and Yarrow (SIAM J. Sci. Stat. Comput. Vol. 13 No. 3 1992)   *
 *                                                                    *
 **********************************************************************/
inline Dusty2D_cState Dusty2D_Quad_Block::Uz(const int ii, const int jj) {
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
         eta1, zeta1, eta2, zeta2, eta, zeta;
  Dusty2D_cState A, B, C, D;
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
 * Dusty2D_Quad_Block::Ez -- Node electrostatic solution state.       *
 *                                                                    *
 * Zingg and Yarrow (SIAM J. Sci. Stat. Comput. Vol. 13 No. 3 1992)   *
 *                                                                    *
 **********************************************************************/
inline Electrostatic2DState Dusty2D_Quad_Block::Ez(const int ii, const int jj) {
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
  A = E[ii-1][jj-1];
  B = E[ii-1][jj  ] - E[ii-1][jj-1]; 
  C = E[ii  ][jj-1] - E[ii-1][jj-1];
  D = E[ii  ][jj  ] + E[ii-1][jj-1] - E[ii-1][jj  ] - E[ii  ][jj-1];
  return A+B*zeta+C*eta+D*zeta*eta;
}

/**********************************************************************
 * Dusty2D_Quad_Block::Wz?? -- Cell node primitive solution states.   *
 **********************************************************************/
inline Dusty2D_pState Dusty2D_Quad_Block::WzNW(const int ii, const int jj) {
  return Wz(ii,jj+1);
}

inline Dusty2D_pState Dusty2D_Quad_Block::WzNE(const int ii, const int jj) {
  return Wz(ii+1,jj+1);
}

inline Dusty2D_pState Dusty2D_Quad_Block::WzSE(const int ii, const int jj) {
  return Wz(ii+1,jj);
}

inline Dusty2D_pState Dusty2D_Quad_Block::WzSW(const int ii, const int jj) {
  return Wz(ii,jj);
}

/**********************************************************************
 * Dusty2D_Quad_Block::Uz?? -- Cell node conserved solution states.   *
 **********************************************************************/
inline Dusty2D_cState Dusty2D_Quad_Block::UzNW(const int ii, const int jj) {
  return Uz(ii,jj+1);
}

inline Dusty2D_cState Dusty2D_Quad_Block::UzNE(const int ii, const int jj) {
  return Uz(ii+1,jj+1);
}

inline Dusty2D_cState Dusty2D_Quad_Block::UzSE(const int ii, const int jj) {
  return Uz(ii+1,jj);
}

inline Dusty2D_cState Dusty2D_Quad_Block::UzSW(const int ii, const int jj) {
  return Uz(ii,jj);
}

/**********************************************************************
 * Dusty2D_Quad_Block::Ez?? -- Cell node conserved solution states.   *
 **********************************************************************/
inline Electrostatic2DState Dusty2D_Quad_Block::EzNW(const int ii, const int jj) {
  return Ez(ii,jj+1);
}

inline Electrostatic2DState Dusty2D_Quad_Block::EzNE(const int ii, const int jj) {
  return Ez(ii+1,jj+1);
}

inline Electrostatic2DState Dusty2D_Quad_Block::EzSE(const int ii, const int jj) {
  return Ez(ii+1,jj);
}

inline Electrostatic2DState Dusty2D_Quad_Block::EzSW(const int ii, const int jj) {
  return Ez(ii,jj);
}

/**********************************************************************
 * Dusty2D_Quad_Block::Wn -- Node primitive variable solution state.  *
 *                                                                    *
 * Zingg and Yarrow (SIAM J. Sci. Stat. Comput. Vol. 13 No. 3 1992)   *
 *                                                                    *
 **********************************************************************/
inline Dusty2D_pState Dusty2D_Quad_Block::Wn(const int ii, const int jj) {
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
         eta1, zeta1, eta2, zeta2, eta, zeta;
  Dusty2D_pState A, B, C, D;
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
 * Dusty2D_Quad_Block::Wn -- Node primitive variable solution state.  *
 *                                                                    *
 * Holmes and Connell (AIAA Paper 1989-1932-CP)                       *
 *                                                                    *
 **********************************************************************/
inline Dusty2D_pState Dusty2D_Quad_Block::Ww(const int ii, const int jj) {
  double wo, wi;
  Vector2D lambda, R;
  Tensor2D I;
  Dusty2D_pState Wo;
  // Determine the contribution to the weighting coefficients for the 
  // standard interpolation support set:
  R = Grid.Cell[ii-1][jj-1].Xc + Grid.Cell[ii  ][jj-1].Xc +
      Grid.Cell[ii-1][jj  ].Xc + Grid.Cell[ii  ][jj  ].Xc -
      FOUR*Grid.Node[ii][jj].X;
  I.xx = sqr(Grid.Cell[ii-1][jj-1].Xc.x - Grid.Node[ii][jj].X.x) +
         sqr(Grid.Cell[ii  ][jj-1].Xc.x - Grid.Node[ii][jj].X.x) +
         sqr(Grid.Cell[ii-1][jj  ].Xc.x - Grid.Node[ii][jj].X.x) +
         sqr(Grid.Cell[ii  ][jj  ].Xc.x - Grid.Node[ii][jj].X.x);
  I.xy = (Grid.Cell[ii-1][jj-1].Xc.x - Grid.Node[ii][jj].X.x)*
         (Grid.Cell[ii-1][jj-1].Xc.y - Grid.Node[ii][jj].X.y) +
         (Grid.Cell[ii  ][jj-1].Xc.x - Grid.Node[ii][jj].X.x)*
         (Grid.Cell[ii  ][jj-1].Xc.y - Grid.Node[ii][jj].X.y) +
         (Grid.Cell[ii-1][jj  ].Xc.x - Grid.Node[ii][jj].X.x)*
         (Grid.Cell[ii-1][jj  ].Xc.y - Grid.Node[ii][jj].X.y) +
         (Grid.Cell[ii  ][jj  ].Xc.x - Grid.Node[ii][jj].X.x)*
         (Grid.Cell[ii  ][jj  ].Xc.y - Grid.Node[ii][jj].X.y);
  I.yy = sqr(Grid.Cell[ii-1][jj-1].Xc.y - Grid.Node[ii][jj].X.y) +
         sqr(Grid.Cell[ii  ][jj-1].Xc.y - Grid.Node[ii][jj].X.y) +
         sqr(Grid.Cell[ii-1][jj  ].Xc.y - Grid.Node[ii][jj].X.y) +
         sqr(Grid.Cell[ii  ][jj  ].Xc.y - Grid.Node[ii][jj].X.y);
  // Determine the weighting coefficients:
  lambda.x = (I.xy*R.y - I.yy*R.x)/(I.xx*I.yy - I.xy*I.xy);
  lambda.y = (I.xy*R.x - I.xx*R.y)/(I.xx*I.yy - I.xy*I.xy);
  // Determine the weighted interpolated state from the standard
  // interpolation support set.
  wi = ONE + lambda*(Grid.Cell[ii-1][jj-1].Xc - Grid.Node[ii][jj].X);
  wo = wi;
  Wo = wi*W[ii-1][jj-1]; 
  wi = ONE + lambda*(Grid.Cell[ii  ][jj-1].Xc - Grid.Node[ii][jj].X);
  wo += wi;
  Wo += wi*W[ii  ][jj-1]; 
  wi = ONE + lambda*(Grid.Cell[ii-1][jj  ].Xc - Grid.Node[ii][jj].X);
  wo += wi;
  Wo += wi*W[ii-1][jj  ]; 
  wi = ONE + lambda*(Grid.Cell[ii  ][jj  ].Xc - Grid.Node[ii][jj].X);
  wo += wi;
  Wo += wi*W[ii  ][jj  ]; 
  // Return interpolated primitive node solution.
  return Wo/wo;
}

/**********************************************************************
 * Dusty2D_Quad_Block::Un -- Node conserved variable solution state.  *
 *                                                                    *
 * Zingg and Yarrow (SIAM J. Sci. Stat. Comput. Vol. 13 No. 3 1992)   *
 *                                                                    *
 **********************************************************************/
inline Dusty2D_cState Dusty2D_Quad_Block::Un(const int ii, const int jj) {
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
         eta1, zeta1, eta2, zeta2, eta, zeta;
  Dusty2D_cState A, B, C, D;
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
 * Dusty2D_Quad_Block::Uw -- Node conserved variable solution state.  *
 *                                                                    *
 * Holmes and Connell (AIAA Paper 1989-1932-CP)                       *
 *                                                                    *
 **********************************************************************/
inline Dusty2D_cState Dusty2D_Quad_Block::Uw(const int ii, const int jj) {
  double wo, wi;
  Vector2D lambda, R;
  Tensor2D I;
  Dusty2D_cState Uo;
  // Determine the contribution to the weighting coefficients for the 
  // standard interpolation support set:
  R = Grid.Cell[ii-1][jj-1].Xc + Grid.Cell[ii  ][jj-1].Xc +
      Grid.Cell[ii-1][jj  ].Xc + Grid.Cell[ii  ][jj  ].Xc - FOUR*Grid.Node[ii][jj].X;
  I.xx = sqr(Grid.Cell[ii-1][jj-1].Xc.x - Grid.Node[ii][jj].X.x) +
         sqr(Grid.Cell[ii  ][jj-1].Xc.x - Grid.Node[ii][jj].X.x) +
         sqr(Grid.Cell[ii-1][jj  ].Xc.x - Grid.Node[ii][jj].X.x) +
         sqr(Grid.Cell[ii  ][jj  ].Xc.x - Grid.Node[ii][jj].X.x);
  I.xy = (Grid.Cell[ii-1][jj-1].Xc.x - Grid.Node[ii][jj].X.x)*
         (Grid.Cell[ii-1][jj-1].Xc.y - Grid.Node[ii][jj].X.y) +
         (Grid.Cell[ii  ][jj-1].Xc.x - Grid.Node[ii][jj].X.x)*
         (Grid.Cell[ii  ][jj-1].Xc.y - Grid.Node[ii][jj].X.y) +
         (Grid.Cell[ii-1][jj  ].Xc.x - Grid.Node[ii][jj].X.x)*
         (Grid.Cell[ii-1][jj  ].Xc.y - Grid.Node[ii][jj].X.y) +
         (Grid.Cell[ii  ][jj  ].Xc.x - Grid.Node[ii][jj].X.x)*
         (Grid.Cell[ii  ][jj  ].Xc.y - Grid.Node[ii][jj].X.y);
  I.yy = sqr(Grid.Cell[ii-1][jj-1].Xc.y - Grid.Node[ii][jj].X.y) +
         sqr(Grid.Cell[ii  ][jj-1].Xc.y - Grid.Node[ii][jj].X.y) +
         sqr(Grid.Cell[ii-1][jj  ].Xc.y - Grid.Node[ii][jj].X.y) +
         sqr(Grid.Cell[ii  ][jj  ].Xc.y - Grid.Node[ii][jj].X.y);
  // Determine the contribution to the weighting coefficients for the 
  // distance-2 WEST pair of cells if requried:
  if (Grid.lfaceN(ii-1,jj-1) < NANO && ii-2 >= 0) {
    R += Grid.Cell[ii-2][jj-1].Xc + Grid.Cell[ii-2][jj  ].Xc - TWO*Grid.Node[ii][jj].X;
    I.xx += sqr(Grid.Cell[ii-2][jj-1].Xc.x - Grid.Node[ii][jj].X.x) +
            sqr(Grid.Cell[ii-2][jj  ].Xc.x - Grid.Node[ii][jj].X.x);
    I.xy += (Grid.Cell[ii-2][jj-1].Xc.x - Grid.Node[ii][jj].X.x)*
            (Grid.Cell[ii-2][jj-1].Xc.y - Grid.Node[ii][jj].X.y) +
            (Grid.Cell[ii-2][jj  ].Xc.x - Grid.Node[ii][jj].X.x)*
            (Grid.Cell[ii-2][jj  ].Xc.y - Grid.Node[ii][jj].X.y);
    I.yy += sqr(Grid.Cell[ii-2][jj-1].Xc.y - Grid.Node[ii][jj].X.y) +
            sqr(Grid.Cell[ii-2][jj  ].Xc.y - Grid.Node[ii][jj].X.y);
  }
  // Determine the contribution to the weighting coefficients for the 
  // distance-2 EAST pair of cells if requried:
  if (Grid.lfaceN(ii,jj-1) < NANO && ii+1 <= ICu+Nghost) {
    R += Grid.Cell[ii+1][jj-1].Xc + Grid.Cell[ii+1][jj  ].Xc - TWO*Grid.Node[ii][jj].X;
    I.xx += sqr(Grid.Cell[ii+1][jj-1].Xc.x - Grid.Node[ii][jj].X.x) +
            sqr(Grid.Cell[ii+1][jj  ].Xc.x - Grid.Node[ii][jj].X.x);
    I.xy += (Grid.Cell[ii+1][jj-1].Xc.x - Grid.Node[ii][jj].X.x)*
            (Grid.Cell[ii+1][jj-1].Xc.y - Grid.Node[ii][jj].X.y) +
            (Grid.Cell[ii+1][jj  ].Xc.x - Grid.Node[ii][jj].X.x)*
            (Grid.Cell[ii+1][jj  ].Xc.y - Grid.Node[ii][jj].X.y);
    I.yy += sqr(Grid.Cell[ii+1][jj-1].Xc.y - Grid.Node[ii][jj].X.y) +
            sqr(Grid.Cell[ii+1][jj  ].Xc.y - Grid.Node[ii][jj].X.y);
  }
  // Determine the contribution to the weighting coefficients for the 
  // distance-2 SOUTH pair of cells if requried:
  if (Grid.lfaceE(ii-1,jj-1) < NANO && jj-2 >= 0) {
    R += Grid.Cell[ii-1][jj-2].Xc + Grid.Cell[ii  ][jj-2].Xc - TWO*Grid.Node[ii][jj].X;
    I.xx += sqr(Grid.Cell[ii-1][jj-2].Xc.x - Grid.Node[ii][jj].X.x) +
            sqr(Grid.Cell[ii  ][jj-2].Xc.x - Grid.Node[ii][jj].X.x);
    I.xy += (Grid.Cell[ii-1][jj-2].Xc.x - Grid.Node[ii][jj].X.x)*
            (Grid.Cell[ii-1][jj-2].Xc.y - Grid.Node[ii][jj].X.y) +
            (Grid.Cell[ii  ][jj-2].Xc.x - Grid.Node[ii][jj].X.x)*
            (Grid.Cell[ii  ][jj-2].Xc.y - Grid.Node[ii][jj].X.y);
    I.yy += sqr(Grid.Cell[ii-1][jj-2].Xc.y - Grid.Node[ii][jj].X.y) +
            sqr(Grid.Cell[ii  ][jj-2].Xc.y - Grid.Node[ii][jj].X.y);
  }
  // Determine the contribution to the weighting coefficients for the 
  // distance-2 NORTH pair of cells if requried:
  if (Grid.lfaceE(ii-1,jj) < NANO && jj+1 <= JCu+Nghost) {
    R += Grid.Cell[ii-1][jj+1].Xc + Grid.Cell[ii  ][jj+1].Xc - TWO*Grid.Node[ii][jj].X;
    I.xx += sqr(Grid.Cell[ii-1][jj+1].Xc.x - Grid.Node[ii][jj].X.x) +
            sqr(Grid.Cell[ii  ][jj+1].Xc.x - Grid.Node[ii][jj].X.x);
    I.xy += (Grid.Cell[ii-1][jj+1].Xc.x - Grid.Node[ii][jj].X.x)*
            (Grid.Cell[ii-1][jj+1].Xc.y - Grid.Node[ii][jj].X.y) +
            (Grid.Cell[ii  ][jj+1].Xc.x - Grid.Node[ii][jj].X.x)*
            (Grid.Cell[ii  ][jj+1].Xc.y - Grid.Node[ii][jj].X.y);
    I.yy += sqr(Grid.Cell[ii-1][jj+1].Xc.y - Grid.Node[ii][jj].X.y) +
            sqr(Grid.Cell[ii  ][jj+1].Xc.y - Grid.Node[ii][jj].X.y);
  }
  // Determine the weighting coefficients:
  lambda.x = (I.xy*R.y - I.yy*R.x)/(I.xx*I.yy - I.xy*I.xy);
  lambda.y = (I.xy*R.x - I.xx*R.y)/(I.xx*I.yy - I.xy*I.xy);
  // Determine the weighted interpolated state from the standard
  // interpolation support set.
  wi = ONE + lambda*(Grid.Cell[ii-1][jj-1].Xc - Grid.Node[ii][jj].X);
  wo = wi;
  Uo = wi*U[ii-1][jj-1]; 
  wi = ONE + lambda*(Grid.Cell[ii  ][jj-1].Xc - Grid.Node[ii][jj].X);
  wo += wi;
  Uo += wi*U[ii  ][jj-1]; 
  wi = ONE + lambda*(Grid.Cell[ii-1][jj  ].Xc - Grid.Node[ii][jj].X);
  wo += wi;
  Uo += wi*U[ii-1][jj  ]; 
  wi = ONE + lambda*(Grid.Cell[ii  ][jj  ].Xc - Grid.Node[ii][jj].X);
  wo += wi;
  Uo += wi*U[ii  ][jj  ]; 
  // Determine the weighted interpolated state from the distance-2 WEST
  // pair of cells if required:
  if (Grid.lfaceN(ii-1,jj-1) < NANO && ii-2 >= 0) {
    wi = ONE + lambda*(Grid.Cell[ii-2][jj-1].Xc - Grid.Node[ii][jj].X);
    wo += wi;
    Uo += wi*U[ii-2][jj-1];
    wi = ONE + lambda*(Grid.Cell[ii-2][jj  ].Xc - Grid.Node[ii][jj].X);
    wo += wi;
    Uo += wi*U[ii-2][jj  ];
  }
  // Determine the weighted interpolated state from the distance-2 EAST
  // pair of cells if required:
  if (Grid.lfaceN(ii,jj-1) < NANO && ii+1 <= ICu+Nghost) {
    wi = ONE + lambda*(Grid.Cell[ii+1][jj-1].Xc - Grid.Node[ii][jj].X);
    wo += wi;
    Uo += wi*U[ii+1][jj-1];
    wi = ONE + lambda*(Grid.Cell[ii+1][jj  ].Xc - Grid.Node[ii][jj].X);
    wo += wi;
    Uo += wi*U[ii+1][jj  ];
  }
  // Determine the weighted interpolated state from the distance-2 SOUTH
  // pair of cells if required:
  if (Grid.lfaceE(ii-1,jj-1) < NANO && jj-2 >= 0) {
    wi = ONE + lambda*(Grid.Cell[ii-1][jj-2].Xc - Grid.Node[ii][jj].X);
    wo += wi;
    Uo += wi*U[ii-1][jj-2];
    wi = ONE + lambda*(Grid.Cell[ii  ][jj-2].Xc - Grid.Node[ii][jj].X);
    wo += wi;
    Uo += wi*U[ii  ][jj-2];
  }
  // Determine the weighted interpolated state from the distance-2 NORTH
  // pair of cells if required:
  if (Grid.lfaceE(ii-1,jj) < NANO && jj+1 <= JCu+Nghost) {
    wi = ONE + lambda*(Grid.Cell[ii-1][jj+1].Xc - Grid.Node[ii][jj].X);
    wo += wi;
    Uo += wi*U[ii-1][jj+1];
    wi = ONE + lambda*(Grid.Cell[ii  ][jj+1].Xc - Grid.Node[ii][jj].X);
    wo += wi;
    Uo += wi*U[ii  ][jj+1];
  }
  // Return interpolated primitive node solution.
  return Uo/wo;
}

/**********************************************************************
 * Dusty2D_Quad_Block::En -- Node variable solution state.            *
 *                                                                    *
 * Holmes and Connell (AIAA Paper 1989-1932-CP)                       *
 *                                                                    *
 **********************************************************************/
inline Electrostatic2DState Dusty2D_Quad_Block::En(const int ii, const int jj) {
  double w1, w2, w3, w4;
  Vector2D lambda, R, X0, X1, X2, X3, X4;
  Tensor2D I;
  // Summarize cell-centres and states.
  X0 = Grid.Node[ii][jj].X;
  X1 = Grid.Cell[ii-1][jj-1].Xc;
  X2 = Grid.Cell[ii  ][jj-1].Xc;
  X3 = Grid.Cell[ii-1][jj  ].Xc;
  X4 = Grid.Cell[ii  ][jj  ].Xc;
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
  return (w1*E[ii-1][jj-1] + w2*E[ii][jj-1] + w3*E[ii-1][jj]+ w4*E[ii][jj])/(w1 + w2 + w3 + w4);
}

/**********************************************************************
 * Dusty2D_Quad_Block::Wn?? -- Cell node primitive solution states.   *
 **********************************************************************/
inline Dusty2D_pState Dusty2D_Quad_Block::WnNW(const int ii, const int jj) {
  return Wn(ii,jj+1);
}

inline Dusty2D_pState Dusty2D_Quad_Block::WnNE(const int ii, const int jj) {
  return Wn(ii+1,jj+1);
}

inline Dusty2D_pState Dusty2D_Quad_Block::WnSE(const int ii, const int jj) {
  return Wn(ii+1,jj);
}

inline Dusty2D_pState Dusty2D_Quad_Block::WnSW(const int ii, const int jj) {
  return Wn(ii,jj);
}

/**********************************************************************
 * Dusty2D_Quad_Block::Un?? -- Cell node conserved solution states.   *
 **********************************************************************/
inline Dusty2D_cState Dusty2D_Quad_Block::UnNW(const int ii, const int jj) {
  return Un(ii,jj+1);
}

inline Dusty2D_cState Dusty2D_Quad_Block::UnNE(const int ii, const int jj) {
  return Un(ii+1,jj+1);
}

inline Dusty2D_cState Dusty2D_Quad_Block::UnSE(const int ii, const int jj) {
  return Un(ii+1,jj);
}

inline Dusty2D_cState Dusty2D_Quad_Block::UnSW(const int ii, const int jj) {
  return Un(ii,jj);
}

/**********************************************************************
 * Dusty2D_Quad_Block::Ew?? -- Cell node solution states.             *
 **********************************************************************/
inline Electrostatic2DState Dusty2D_Quad_Block::EnNW(const int ii, const int jj) {
  return En(ii,jj+1);
}

inline Electrostatic2DState Dusty2D_Quad_Block::EnNE(const int ii, const int jj) {
  return En(ii+1,jj+1);
}

inline Electrostatic2DState Dusty2D_Quad_Block::EnSE(const int ii, const int jj) {
  return En(ii+1,jj);
}

inline Electrostatic2DState Dusty2D_Quad_Block::EnSW(const int ii, const int jj) {
  return En(ii,jj);
}

/**********************************************************************
 * Dusty2D_Quad_Block::evaluate_limiters -- Set flag to evaluate      *
 *                                          limiters.                 *
 **********************************************************************/
inline void Dusty2D_Quad_Block::evaluate_limiters(void) {
  Freeze_Limiter = OFF; 
}

/**********************************************************************
 * Dusty2D_Quad_Block::freeze_limiters -- Set flag to freeze limiters.*
 **********************************************************************/
inline void Dusty2D_Quad_Block::freeze_limiters(void) {
  Freeze_Limiter = ON; 
}

/**********************************************************************
 * Dusty2D_Quad_Block -- Input-output operators.                      *
 **********************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Dusty2D_Quad_Block &SolnBlk) {
  out_file << SolnBlk.Grid;
  out_file << SolnBlk.NCi << " " << SolnBlk.ICl << " " << SolnBlk.ICu << endl;
  out_file << SolnBlk.NCj << " " << SolnBlk.JCl << " " << SolnBlk.JCu << endl;
  out_file << SolnBlk.Nghost            << endl;
  out_file << SolnBlk.Axisymmetric      << endl;
  out_file << SolnBlk.Particles         << endl;
  out_file << SolnBlk.PhaseInteraction  << endl;
  out_file << SolnBlk.Flow_Type           << endl;
  out_file << SolnBlk.Electrostatic     << endl;
  out_file << SolnBlk.Vwall << endl;
  out_file << SolnBlk.Twall << endl;
  if (SolnBlk.NCi == 0 || SolnBlk.NCj == 0) return out_file;
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      out_file << SolnBlk.U[i][j] << endl;
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
			     Dusty2D_Quad_Block &SolnBlk) {
  int ni, il, iu, nj, jl, ju, ng, dum;
  Grid2D_Quad_Block New_Grid;
  in_file >> New_Grid;
  in_file.setf(ios::skipws);
  in_file >> ni >> il >> iu;
  in_file >> nj >> jl >> ju;
  in_file >> ng;
  in_file >> SolnBlk.Axisymmetric;
  in_file >> SolnBlk.Particles;
  in_file >> SolnBlk.PhaseInteraction;
  in_file >> SolnBlk.Flow_Type;
  in_file >> SolnBlk.Electrostatic;
  in_file >> dum;
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
      for (int k = 0; k < NUMBER_OF_RESIDUAL_VECTORS_DUSTY2D; k++) {
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
  return in_file;
}

/**********************************************************************
 *                                                                    *
 * MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING.                     *
 *                                                                    *
 **********************************************************************/

/**********************************************************************
 * Dusty2D_Quad_Block::NumVar -- Returns number of state variables.   *
 **********************************************************************/
inline int Dusty2D_Quad_Block::NumVar(void) {
  return NUM_VAR_DUSTY2D;
}

/**********************************************************************
 * Dusty2D_Quad_Block::LoadSendBuffer -- Loads send message buffer.   *
 **********************************************************************/
inline int Dusty2D_Quad_Block::LoadSendBuffer(double *buffer,
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
      for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = U[i][j][k];
      }
    }
  }
  return 0;
}

/**********************************************************************
 * Dusty2D_Quad_Block::LoadSendBuffer_F2C --                          *
 *                     Loads send message buffer for fine to coarse   *
 *                     block message passing.                         *
 **********************************************************************/
inline int Dusty2D_Quad_Block::LoadSendBuffer_F2C(double *buffer,
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
      for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
 * Dusty2D_Quad_Block::LoadSendBuffer_C2F --                          *
 *                     Loads send message buffer for coarse to fine   *
 *                     block message passing.                         *
 **********************************************************************/
inline int Dusty2D_Quad_Block::LoadSendBuffer_C2F(double *buffer,
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
  Dusty2D_pState Wfine;
  Dusty2D_cState Ufine;
  int Limiter = LIMITER_VENKATAKRISHNAN;

  if (j_inc > 0) {
    if (i_inc > 0) {
      for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
	  // Perform limited linear least squares reconstruction in cell (i,j_min).
	  LinearSubcellReconstruction(i,j_min,Limiter);
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
	    for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
	    dX = Grid.centroidSE(i,j) - Grid.Cell[i  ][j].Xc;
	    Wfine = W[i][j] + (phi[i][j]^dWdx[i][j])*dX.x +
	                      (phi[i][j]^dWdy[i][j])*dX.y;
	    Ufine = Wfine.U();
	    for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
	    for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) { 
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
	    for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
	  LinearSubcellReconstruction(i,j_min,Limiter);
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
		Grid.Node[i+1][j_min].X+
		Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	        Grid.Cell[i][j_min].Xc;
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) { 
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
	  LinearSubcellReconstruction(i,j_min,Limiter);
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
		Grid.Cell[i][j_min].Xc+
		Grid.Node[i][j_min+1].X+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	        Grid.Cell[i][j_min].Xc;
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) { 
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
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
      } else {
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
	  // Perform limited linear least squares reconstruction in cell (i,j_min).
	  LinearSubcellReconstruction(i,j_min,Limiter);
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
		Grid.Node[i+1][j_min+1].X)/FOUR -
	        Grid.Cell[i][j_min].Xc;
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) { 
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
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
	  Wfine = W[i][j_min] + (phi[i][j_min]^dWdx[i][j_min])*dX.x +
  	                        (phi[i][j_min]^dWdy[i][j_min])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
	  LinearSubcellReconstruction(i_min,j,Limiter);
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		Grid.Node[i_min+1][j].X+
		Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
	        Grid.Cell[i_min][j].Xc;
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                                (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                                (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) { 
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
	  LinearSubcellReconstruction(i_min,j,Limiter);
	  // Evaluate NW sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc+
		Grid.Node[i_min][j+1].X+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                Grid.Cell[i_min][j].Xc;
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                                (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) { 
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
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
      } else {
	for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
	  // Perform limited linear least squares reconstruction in cell (i_min,j).
	  LinearSubcellReconstruction(i_min,j,Limiter);
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i_min][j].Xc+
	        HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
		Grid.Node[i_min+1][j+1].X)/FOUR -
	        Grid.Cell[i_min][j].Xc;
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                                (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) { 
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
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                                (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
	  Wfine = W[i_min][j] + (phi[i_min][j]^dWdx[i_min][j])*dX.x +
                                (phi[i_min][j]^dWdy[i_min][j])*dX.y;
	  Ufine = Wfine.U();
	  for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
 * Dusty2D_Quad_Block::UnloadReceiveBuffer --                         *
 *                     Unloads receive message buffer.                *
 **********************************************************************/
inline int Dusty2D_Quad_Block::UnloadReceiveBuffer(double *buffer,
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
      for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
 * Dusty2D_Quad_Block::UnloadReceiveBuffer_F2C --                     *
 *                     Unloads receive message buffer for fine to     *
 *                     coarse block message passing.                  *
 **********************************************************************/
inline int Dusty2D_Quad_Block::UnloadReceiveBuffer_F2C(double *buffer,
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
      for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
 * Dusty2D_Quad_Block::UnloadReceiveBuffer_C2F --                     *
 *                     Unloads receive message buffer for coarse to   *
 *                     fine block message passing.                    *
 **********************************************************************/
inline int Dusty2D_Quad_Block::UnloadReceiveBuffer_C2F(double *buffer,
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
      for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
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
 * Dusty2D_Quad_Block::LinearSubcellReconstruction --                 *
 *                     Performs the subcell linear reconstruction of  *
 *                     solution state within a given cell (i,j) of    *
 *                     the computational mesh for the specified       *
 *                     quadrilateral solution block.                  *
 **********************************************************************/
inline void Dusty2D_Quad_Block::LinearSubcellReconstruction(const int i,
							    const int j,
							    const int Limiter) {

  int n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4], phi_n;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  Dusty2D_pState DU, DUDx_ave, DUDy_ave;

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
      for (int n = 1; n <= NumVar(); n++) {
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

/**********************************************************************
 * Dusty2D_Quad_Block::LoadSendBuffer_Flux_F2C --                     *
 *                 Loads send message buffer for fine to coarse block *
 *                 message passing of conservative solution fluxes.   *
 **********************************************************************/
inline int Dusty2D_Quad_Block::LoadSendBuffer_Flux_F2C(double *buffer,
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
      for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = FluxS[i  ][k] + FluxS[i+1][k];
      }
    }
  } else if (j_min == j_max && j_min == JCu) {
    for (int i = i_min; ((i_inc+2)/4) ? (i < i_max):(i > i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = FluxN[i  ][k] + FluxN[i+1][k];
      }
    }
  } else if (i_min == i_max && i_min == ICl) {
    for (int j = j_min; ((j_inc+2)/4) ? (j < j_max):(j > j_max); j += j_inc) {
      for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = FluxW[j][k] + FluxW[j+1][k];
      }
    }
  } else if (i_min == i_max && i_min == ICu) {
    for (int j = j_min; ((j_inc+2)/4) ? (j < j_max):(j > j_max); j += j_inc) {
      for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = FluxE[j][k] + FluxE[j+1][k];
      }
    }
  }
  return 0;
}

/**********************************************************************
 * Dusty2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C --                *
 *                     Unloads receive message buffer for fine to     *
 *                     coarse block message passing of conservative   *
 *                     solution fluxes.                               *
 **********************************************************************/
inline int Dusty2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C(double *buffer,
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
      for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	FluxS[i][k] = - buffer[buffer_count] - FluxS[i][k];
      }
    }
  } else if (j_min == j_max && j_min == JCu) {
    for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	FluxN[i][k] = - buffer[buffer_count] - FluxN[i][k];
      }
    }
  } else if (i_min == i_max && i_min == ICl) {
    for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
      for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	FluxW[j][k] = - buffer[buffer_count] - FluxW[j][k];
      }
    }
  } else if (i_min == i_max && i_min == ICu) {
    for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
      for (int k = 1; k <= NUM_VAR_DUSTY2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	FluxE[j][k] = - buffer[buffer_count] - FluxE[j][k];
      }
    }
  }
  return 0;
}

/**********************************************************************
 * Dusty2D_Quad_Block -- Single Block External Subroutines.           *
 **********************************************************************/

extern void Broadcast_Solution_Block(Dusty2D_Quad_Block &SolnBlk);

#ifdef _MPI_VERSION
extern void Broadcast_Solution_Block(Dusty2D_Quad_Block &SolnBlk,
                                     MPI::Intracomm &Communicator,
                                     const int Source_CPU);
#endif

extern void Copy_Solution_Block(Dusty2D_Quad_Block &SolnBlk1,
		                Dusty2D_Quad_Block &SolnBlk2);

extern int Prolong_Solution_Block(Dusty2D_Quad_Block &SolnBlk_Fine,
				  Dusty2D_Quad_Block &SolnBlk_Original,
				  const int Sector);

extern int Restrict_Solution_Block(Dusty2D_Quad_Block &SolnBlk_Coarse,
				   Dusty2D_Quad_Block &SolnBlk_Original_SW,
				   Dusty2D_Quad_Block &SolnBlk_Original_SE,
				   Dusty2D_Quad_Block &SolnBlk_Original_NW,
				   Dusty2D_Quad_Block &SolnBlk_Original_NE);

extern void ICs(Dusty2D_Quad_Block &SolnBlk,
		Dusty2D_Input_Parameters &IP,
                Dusty2D_pState *Wo);

extern int Copy_Electrostatic_Field_Variables(Dusty2D_Quad_Block &SolnBlk,
					      Electrostatic2D_Quad_Block &ES_SolnBlk);

extern void BCs(Dusty2D_Quad_Block &SolnBlk,
		Dusty2D_Input_Parameters &IP);

extern double CFL(Dusty2D_Quad_Block &SolnBlk,
                  Dusty2D_Input_Parameters &IP);

extern void Set_Global_TimeStep(Dusty2D_Quad_Block &SolnBlk,
                                const double &Dt_min);

extern double L1_Norm_Residual(Dusty2D_Quad_Block &SolnBlk);

extern double L2_Norm_Residual(Dusty2D_Quad_Block &SolnBlk);

extern double Max_Norm_Residual(Dusty2D_Quad_Block &SolnBlk);

extern void Linear_Reconstruction_GreenGauss(Dusty2D_Quad_Block &SolnBlk,
                                             const int i,
                                             const int j,
					     const int Limiter);

extern void Linear_Reconstruction_GreenGauss(Dusty2D_Quad_Block &SolnBlk,
					     const int Limiter);

extern void Linear_Reconstruction_LeastSquares(Dusty2D_Quad_Block &SolnBlk,
                                               const int i,
                                               const int j,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares_2(Dusty2D_Quad_Block &SolnBlk,
                                                 const int i,
                                                 const int j,
					         const int Limiter);

extern void Linear_Reconstruction_LeastSquares(Dusty2D_Quad_Block &SolnBlk,
					       const int Limiter);

extern void Residual_Smoothing(Dusty2D_Quad_Block &SolnBlk,
                               const int k_residual,
			       double &epsilon,
                               const int number_of_Gauss_Seidel_iterations);

extern void Calculate_Refinement_Criteria(double *refinement_criteria,
					  Dusty2D_Input_Parameters &IP,
                                          int &number_refinement_criteria,
                                          Dusty2D_Quad_Block &SolnBlk);

extern void Fix_Refined_Block_Boundaries(Dusty2D_Quad_Block &SolnBlk,
                                         const int Fix_North_Boundary,
                                         const int Fix_South_Boundary,
                                         const int Fix_East_Boundary,
                                         const int Fix_West_Boundary);

extern void Unfix_Refined_Block_Boundaries(Dusty2D_Quad_Block &SolnBlk);

extern void Apply_Boundary_Flux_Corrections(Dusty2D_Quad_Block &SolnBlk,
                                            const int Number_Neighbours_North_Boundary,
                                            const int Number_Neighbours_South_Boundary,
                                            const int Number_Neighbours_East_Boundary,
                                            const int Number_Neighbours_West_Boundary);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Dusty2D_Quad_Block &SolnBlk,
                                                                const int i_stage,
								Dusty2D_Input_Parameters &IP,
                                                                const int Number_Neighbours_North_Boundary,
                                                                const int Number_Neighbours_South_Boundary,
                                                                const int Number_Neighbours_East_Boundary,
                                                                const int Number_Neighbours_West_Boundary);

extern int dUdt_Residual_Evaluation(Dusty2D_Quad_Block &SolnBlk,
				    Dusty2D_Input_Parameters &IP);

extern int dUdt_Multistage_Explicit(Dusty2D_Quad_Block &SolnBlk,
   	                            const int i_stage,
                                    Dusty2D_Input_Parameters &IP);

extern int Update_Solution_Multistage_Explicit(Dusty2D_Quad_Block &SolnBlk,
   	                                       const int i_stage,
					       Dusty2D_Input_Parameters &IP);

extern int Multi_Velocity_Component_Particle_Phase_Switch(Dusty2D_Quad_Block &SolnBlk,
							  Dusty2D_Input_Parameters &IP);

extern int Determine_Conservation_Properties(Dusty2D_Quad_Block &SolnBlk,
					     Dusty2D_cState &Umass);

/**********************************************************************
 * Dusty2D_Quad_Block -- Multiple Block External Subroutines.         *
 **********************************************************************/

extern Dusty2D_Quad_Block* Allocate(Dusty2D_Quad_Block *Soln_ptr,
                                    Dusty2D_Input_Parameters &Input_Parameters);

extern Dusty2D_Quad_Block* Deallocate(Dusty2D_Quad_Block *Soln_ptr,
                                      Dusty2D_Input_Parameters &Input_Parameters);

extern Dusty2D_Quad_Block* CreateInitialSolutionBlocks(Grid2D_Quad_Block **InitMeshBlk,
						       Dusty2D_Quad_Block *Soln_ptr,
						       Dusty2D_Input_Parameters &Input_Parameters,
						       QuadTreeBlock_DataStructure &QuadTree,
						       AdaptiveBlockResourceList &GlobalSolnBlockList,
						       AdaptiveBlock2D_List &LocalSolnBlockList);

extern void ICs(Dusty2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
                Dusty2D_Input_Parameters &Input_Parameters);

extern int Copy_Electrostatic_Field_Variables(Dusty2D_Quad_Block *Soln_ptr,
					      Electrostatic2D_Quad_Block *ES_Soln_ptr,
					      AdaptiveBlock2D_List &Soln_Block_List);

extern void BCs(Dusty2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
		Dusty2D_Input_Parameters &Input_Parameters);

extern double CFL(Dusty2D_Quad_Block *Soln_ptr,
                  AdaptiveBlock2D_List &Soln_Block_List,
		  Dusty2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(Dusty2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                const double &Dt_min);

extern double L1_Norm_Residual(Dusty2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double L2_Norm_Residual(Dusty2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double Max_Norm_Residual(Dusty2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List);

extern void Evaluate_Limiters(Dusty2D_Quad_Block *Soln_ptr,
                              AdaptiveBlock2D_List &Soln_Block_List);

extern void Freeze_Limiters(Dusty2D_Quad_Block *Soln_ptr,
                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Residual_Smoothing(Dusty2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Dusty2D_Input_Parameters &Input_Parameters,
   	                       const int I_Stage);

extern void Apply_Boundary_Flux_Corrections(Dusty2D_Quad_Block *Soln_ptr,
                                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Dusty2D_Quad_Block *Soln_ptr,
                                                                AdaptiveBlock2D_List &Soln_Block_List,
                                                                Dusty2D_Input_Parameters &Input_Parameters,
   	                                                        const int I_Stage);

extern int dUdt_Residual_Evaluation(Dusty2D_Quad_Block *Soln_ptr,
				    AdaptiveBlockResourceList &Global_Soln_Block_List,
                                    AdaptiveBlock2D_List &Local_Soln_Block_List,
                                    Dusty2D_Input_Parameters &Input_Parameters);

extern int dUdt_Multistage_Explicit(Dusty2D_Quad_Block *Soln_ptr,
				    AdaptiveBlockResourceList &Global_Soln_Block_List,
                                    AdaptiveBlock2D_List &Local_Soln_Block_List,
                                    Dusty2D_Input_Parameters &Input_Parameters,
   	                            const int I_Stage);

extern int Update_Solution_Multistage_Explicit(Dusty2D_Quad_Block *Soln_ptr,
                                               AdaptiveBlock2D_List &Soln_Block_List,
                                               Dusty2D_Input_Parameters &Input_Parameters,
   	                                       const int I_Stage);

extern int Adaptive_Mesh_Refinement(Dusty2D_Quad_Block *Soln_ptr,
				    Dusty2D_Input_Parameters &Input_Parameters,
				    QuadTreeBlock_DataStructure &QuadTree,
				    AdaptiveBlockResourceList &GlobalSolnBlockList,
				    AdaptiveBlock2D_List &LocalSolnBlockList,
				    const int Set_New_Refinement_Flags);

extern int Initial_Adaptive_Mesh_Refinement(Dusty2D_Quad_Block *Soln_ptr,
					    Dusty2D_Input_Parameters &Input_Parameters,
					    QuadTreeBlock_DataStructure &QuadTree,
					    AdaptiveBlockResourceList &GlobalSolnBlockList,
					    AdaptiveBlock2D_List &LocalSolnBlockList);

extern int Uniform_Adaptive_Mesh_Refinement(Dusty2D_Quad_Block *Soln_ptr,
					    Dusty2D_Input_Parameters &Input_Parameters,
					    QuadTreeBlock_DataStructure &QuadTree,
					    AdaptiveBlockResourceList &GlobalSolnBlockList,
					    AdaptiveBlock2D_List &LocalSolnBlockList);

extern int Boundary_Adaptive_Mesh_Refinement(Dusty2D_Quad_Block *Soln_ptr,
					     Dusty2D_Input_Parameters &Input_Parameters,
					     QuadTreeBlock_DataStructure &QuadTree,
					     AdaptiveBlockResourceList &GlobalSolnBlockList,
					     AdaptiveBlock2D_List &LocalSolnBlockList);

extern int Multi_Velocity_Component_Particle_Phase_Switch(Dusty2D_Quad_Block *Soln_ptr,
							  Dusty2D_Input_Parameters &Input_Parameters,
							  AdaptiveBlock2D_List &LocalSolnBlockList);

extern int Determine_Conservation_Properties(Dusty2D_Quad_Block *Soln_ptr,
					     AdaptiveBlock2D_List &LocalSolnBlockList,
					     Dusty2D_cState &Umass);

/**********************************************************************
 * Dusty2D_Quad_Block -- IO Single Block External Subroutines.        *
 **********************************************************************/

extern void Write_Solution_Block(Dusty2D_Quad_Block &SolnBlk,
	                         ostream &Out_File);

extern void Read_Solution_Block(Dusty2D_Quad_Block &SolnBlk,
	                        istream &In_File);

extern void Output_Tecplot(Dusty2D_Quad_Block &SolnBlk,
			   Dusty2D_Input_Parameters &IP,
		           const int Number_of_Time_Steps,
                           const double &Time,
                           const int Block_Number,
                           const int Output_Title,
	                   ostream &Out_File);

extern void Output_Cells_Tecplot(Dusty2D_Quad_Block &SolnBlk,
			         Dusty2D_Input_Parameters &IP,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_Nodes_Tecplot(Dusty2D_Quad_Block &SolnBlk,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
				 const int Block_Number,
				 const int Output_Title,
				 ostream &Out_File);

extern void Output_Gradients_Tecplot(Dusty2D_Quad_Block &SolnBlk,
				     const int Number_of_Time_Steps,
				     const double &Time,
				     const int Block_Number,
				     const int Output_Title,
				     ostream &Out_File);

extern void Output_Quasi3D_Tecplot(Dusty2D_Quad_Block &SolnBlk,
				   Dusty2D_Input_Parameters &IP,
				   const int Number_of_Time_Steps,
				   const double &Time,
				   const int Block_Number,
				   const int Output_Title,
				   ostream &Out_File);

extern void Output_Ringleb(Dusty2D_Quad_Block &SolnBlk,
			   const int Block_Number,
			   const int Output_Title,
			   ostream &Out_File,
			   double &l1_norm,
			   double &l2_norm,
			   double &max_norm,
			   double &area);

extern void Output_Viscous_Channel(Dusty2D_Quad_Block &SolnBlk,
				   const int Block_Number,
				   const int Output_Title,
				   ostream &Out_File,
				   double &l1_norm,
				   double &l2_norm,
				   double &max_norm,
				   const Vector2D Vwall,
				   const double dp,
				   const double length,
				   const double height);

extern void Output_Viscous_Pipe(Dusty2D_Quad_Block &SolnBlk,
				const int Block_Number,
				const int Output_Title,
				ostream &Out_File,
				double &l1_norm,
				double &l2_norm,
				double &max_norm,
				const double dp,
				const double length,
				const double radius);

extern void Output_Turbulent_Pipe_Tecplot(Dusty2D_Quad_Block &SolnBlk,
					  const int Block_Number,
					  const int Output_Title,
					  const int Output_Data,
					  ostream &Out_File,
					  const double &Re,
					  const double &Pipe_Radius,
					  const int &variable_flag);

extern void Output_Flat_Plate(Dusty2D_Quad_Block &SolnBlk,
			      const int Block_Number,
			      const int Output_Title_Soln,
			      ostream &Out_File_Soln,
			      const int Output_Title_Skin,
			      ostream &Out_File_Skin,
			      const Dusty2D_pState &Winf,
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

extern void Output_Driven_Cavity_Flow(Dusty2D_Quad_Block &SolnBlk,
				      const int Block_Number,
				      const int Output_Title,
				      ostream &Out_File_u,
				      ostream &Out_File_v,
				      const double &Re,
				      const Vector2D &Vwall,
				      const double &length);

/**********************************************************************
 * Dusty2D_Quad_Block -- IO Multiple Block External Subroutines.      *
 **********************************************************************/

extern int Read_Restart_Solution(Dusty2D_Quad_Block *Soln_ptr,
                                 AdaptiveBlock2D_List &Soln_Block_List,
                                 Dusty2D_Input_Parameters &IP,
		                 int &Number_of_Time_Steps,
                                 double &Time,
                                 CPUTime &CPU_Time);

extern int Write_Restart_Solution(Dusty2D_Quad_Block *Soln_ptr,
                                  AdaptiveBlock2D_List &Soln_Block_List,
                                  Dusty2D_Input_Parameters &IP,
		                  const int Number_of_Time_Steps,
                                  const double &Time,
                                  const CPUTime &CPU_Time);

extern int Output_Tecplot(Dusty2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          Dusty2D_Input_Parameters &IP,
		          const int Number_of_Time_Steps,
                          const double &Time);

extern int Output_Cells_Tecplot(Dusty2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                Dusty2D_Input_Parameters &IP,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Nodes_Tecplot(Dusty2D_Quad_Block *Soln_ptr,
				AdaptiveBlock2D_List &Soln_Block_List,
				Dusty2D_Input_Parameters &IP,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Gradients_Tecplot(Dusty2D_Quad_Block *Soln_ptr,
				    AdaptiveBlock2D_List &Soln_Block_List,
				    Dusty2D_Input_Parameters &IP,
				    const int Number_of_Time_Steps,
				    const double &Time);

extern int Output_Quasi3D_Tecplot(Dusty2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  Dusty2D_Input_Parameters &IP,
				  const int Number_of_Time_Steps,
				  const double &Time);

extern int Output_Mesh_Tecplot(Dusty2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Dusty2D_Input_Parameters &IP,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Mesh_Gnuplot(Dusty2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Dusty2D_Input_Parameters &IP,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Ringleb(Dusty2D_Quad_Block *Soln_ptr,
			  AdaptiveBlock2D_List &Soln_Block_List,
			  Dusty2D_Input_Parameters &IP);

extern int Output_Viscous_Channel(Dusty2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  Dusty2D_Input_Parameters &IP);

extern int Output_Viscous_Pipe(Dusty2D_Quad_Block *Soln_ptr,
			       AdaptiveBlock2D_List &Soln_Block_List,
			       Dusty2D_Input_Parameters &IP);

extern int Output_Turbulent_Pipe_Tecplot(Dusty2D_Quad_Block *Soln_ptr,
					 AdaptiveBlock2D_List &Soln_Block_List,
					 Dusty2D_Input_Parameters &IP);

extern int Output_Flat_Plate(Dusty2D_Quad_Block *Soln_ptr,
			     AdaptiveBlock2D_List &Soln_Block_List,
			     Dusty2D_Input_Parameters &IP);

extern int Output_Driven_Cavity_Flow(Dusty2D_Quad_Block *Soln_ptr,
				     AdaptiveBlock2D_List &Soln_Block_List,
				     Dusty2D_Input_Parameters &IP);

/**********************************************************************
 * Dusty2D_Quad_Block -- Grid Multiple Block External Subroutines.    *
 **********************************************************************/

extern Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
					    Dusty2D_Input_Parameters &IP);

extern Grid2D_Quad_Block** Broadcast_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
						      Dusty2D_Input_Parameters &IP);

extern int Write_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                             Dusty2D_Input_Parameters &IP);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
							    Dusty2D_Input_Parameters &IP);

extern int Write_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                  Dusty2D_Input_Parameters &IP);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
						 Dusty2D_Input_Parameters &IP);

extern int Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                          Dusty2D_Input_Parameters &IP);

extern int Output_Nodes_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                Dusty2D_Input_Parameters &IP);

extern int Output_Cells_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                Dusty2D_Input_Parameters &IP);

/**********************************************************************
 * Dusty2D_Quad_Block -- Turbulence Single Block External Subroutines.*
 **********************************************************************/

extern int Turbulence_ICs(Dusty2D_Quad_Block &SolnBlk,
			  Dusty2D_Input_Parameters &IP,
			  Dusty2D_pState *Wo);

extern int Turbulence_BCs(Dusty2D_Quad_Block &SolnBlk,
			  Dusty2D_Input_Parameters &IP);

extern int Apply_Turbulence_BCs(Dusty2D_Quad_Block &SolnBlk,
				Dusty2D_Input_Parameters &IP,
				const int &Turbulent_BCtype,
				const int &i, const int &j);

extern int Turbulence_Zero_Residual(Dusty2D_Quad_Block &SolnBlk,
				    const int i_stage,
				    Dusty2D_Input_Parameters &IP);

extern int Zero_Turbulence_Residuals(Dusty2D_Quad_Block &SolnBlk,
				     Dusty2D_Input_Parameters &IP,
				     const int &Turbulent_BCtype,
				     const int &i, const int &j,
				     const int &k_residual);

/**********************************************************************
 * Dusty2D_Quad_Block -- Turbulence Multiple Block External           *
 *                       Subroutines.                                 *
 **********************************************************************/

extern int Turbulence_BCs(Dusty2D_Quad_Block *Soln_ptr,
			  AdaptiveBlock2D_List &Soln_Block_List,
			  Dusty2D_Input_Parameters &Input_Parameters);

/**********************************************************************
 * Dusty2D_Quad_Block -- Solvers.                                     *
 **********************************************************************/

extern int Dusty2DQuadSolver(char *Input_File_Name_ptr,
			     int batch_flag);

extern int Open_Conservation_File(ofstream &Conservation_File,
				  char *File_Name,
				  const int Append_to_File,
				  const Dusty2D_cState &U);

extern int Close_Conservation_File(ofstream &Conservation_File);

extern void Output_Conservation_to_File(ostream &Conservation_File,
					const int Number_of_Time_Steps,
					const double &Time,
					const CPUTime &CPU_Time,
					const Dusty2D_cState &U);

#endif // _DUSTY2D_QUAD_INCLUDED
