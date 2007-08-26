/********************** LESPremixed2DQuad.h **********************************
  Header file defining 2D Thermally Perfect Multiple Species Favre-Filtered
  Navier-Stokes  Solution State Classes.

  - based on Chem2DQuad.h 
************************************************************************/

#ifndef _LESPREMIXED2D_QUAD_INCLUDED
#define _LESPREMIXED2D_QUAD_INCLUDED

/* Include 2D LES Premixed state, 2D cell, and 2D quadrilateral block 
   grid, quadtree, and 2D LES Premixed input header files. */
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

//--------- Modified for LES Premixed ------------------
// Added and Header file for LES Premixed

#ifndef _LESPREMIXED2D_STATE_INCLUDED
#include "LESPremixed2DState.h"
#endif // _LESPREMIXED2D_STATE_INCLUDED

#ifndef _LESPREMIXED2D_INPUT_INCLUDED
#include "LESPremixed2DInput.h"
#endif // _LESPREMIXED2D_INPUT_INCLUDED

#ifndef _FFTW_INCLUDED
#include "fftw3.h"
#endif //_FFTW_INCLUDED

// Include the turbulent wall data header file.

#ifndef _TURBULENT2D_WALLDATA_INCLUDED
#include "../Turbulent2D/Turbulent2DWallData.h"
#endif // _TURBULENT2D_WALLDATA_INCLUDED


/* Define the structures and classes. */

#define	NUMBER_OF_RESIDUAL_VECTORS_LESPREMIXED2D    3  //K for dUdt[i][j][K]

/*!
 * Class: LESPremixed2D_Quad_Block
 *
 * @brief Class definition of the 2D LES Premixed solution blocks.
 *
 * \verbatim
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
 *      Ut      -- Return solution state at the current *
 *                 physical time.                       *
 *                 (Used in the treatment of physical   *
 *                  time for dual time-stepping)        *
 *     Uold     -- Return solution state at the         *
 *                 previous physical time step.         *
 *                 (Used in the treatment of physical   *
 *                  time for dual time-stepping)        *
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
 *  Gravity     -- Return Gravity Flag (1 ON, 0 OFF)    *
 * debug_level  -- Level of Verboseness (0 none, 1,2..) * 
 *  Moving_wall_veloicty -- Moving wall boundary        *
 *                          condition velocity          *             
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
 *      S -- a 2D LES Premixed solution                 *
 *                                                      *
 * S = S;                                               *
 * cout << S; (output function)                         *
 * cin  >> S; (input function)                          *
 *                                                      *
 * \endverbatim
 */
class LESPremixed2D_Quad_Block{
  private:
  public:
  
  /***************** CLASS DATA *****************************************************/
  //@{ @name Solution state arrays:
  LESPremixed2D_pState             **W;   //!< Primitive solution state.
  LESPremixed2D_cState             **U;   //!< Conserved solution state. 
  //@} 
  
  LESPremixed2D_cState            **Uo;   //!< Initial solution state.
  LESPremixed2D_cState            **Ut;   //!< Solution state at the current time step.
  LESPremixed2D_cState          **Uold;   //!< Solution state at the previous time step.
 
  double                          **dt;   //!< Local time step.  
  LESPremixed2D_cState         ***dUdt;   //!< Solution residual.


  //@{ @name Boundary solution flux arrays:
  LESPremixed2D_cState   *FluxN,*FluxS,   //!< Boundary solution fluxes.
                         *FluxE,*FluxW; 
  //@}  

  //@{ @name Boundary condtion reference states:
  LESPremixed2D_pState       *WoN,*WoS,   //!< Boundary condition reference states for
                             *WoE,*WoW;   //!< north, south, east, & west boundaries.
  //@}  

  //@{ @name Solution gradient arrays:
  LESPremixed2D_pState        **dWdx;     //!< (x-direction).
  LESPremixed2D_pState        **dWdy;     //!< (y-direction).
  LESPremixed2D_pState         **phi;     //!< Solution slope limiter.
  //@}  

  //Store face gradients for Diamond Path & Jacobian formation
  LESPremixed2D_pState       **dWdx_faceN;   //!< north cell face(x-direction).
  LESPremixed2D_pState       **dWdx_faceE;   //!< east  cell face(x-direction).
  LESPremixed2D_pState       **dWdx_faceS;   //!< south cell face(x-direction).
  LESPremixed2D_pState       **dWdx_faceW;   //!< west  cell face(x-direction).
  LESPremixed2D_pState       **dWdy_faceN;   //!< north cell face(y-direction).
  LESPremixed2D_pState       **dWdy_faceE;   //!< east  cell face(y-direction).
  LESPremixed2D_pState       **dWdy_faceS;   //!< south cell face(y-direction).
  LESPremixed2D_pState       **dWdy_faceW;   //!< west  cell face(y-direction).
  //@}  

  //@{ @name Grid block information:
  Grid2D_Quad_Block       Grid;   //!< 2D quadrilateral grid geometry.
  int              NCi,ICl,ICu;   //!< i-direction cell counters.
  int              NCj,JCl,JCu;   //!< j-direction cell counters.
  int                   Nghost;   //!< Number of ghost cells.
  //@}


  //@{ @name Turbulence wall data arrays:
  Turbulent2DWallData   **Wall; //!< Turbulent wall data.
  //@}

  //@{ @name FLAGS ( These are all esentially "static" Input parameters put in the SolnBlk for ease of access)
  int             Axisymmetric;   //!< Axisymmetric flow indicator.
  int                  Gravity;   //!< Gravity flag
  int                Flow_Type;   //!< Laminar, Viscous, Turbulent etc..
  int           Wall_Functions;
  int              debug_level;   //Debug level flag (0=none, 1,2,3,.. level of verboseness)
  int           Freeze_Limiter;   //Limiter freezing indicator (Multigrid/NKS) 0/1
  static int residual_variable;          //!< 1 = rho, 2,3 = rhou, 4 = E
  static int Number_of_Residual_Norms;   //!< (default 4 )
  double  Moving_wall_velocity;   //Moving Wall BC velocity
  double     Pressure_gradient;   //Moving Wall BC velocity
  //@}

  /****************** MEMBER FUNCTIONS ****************************************************/  
  /* Creation, copy, and assignment constructors. */
  LESPremixed2D_Quad_Block(void) {
    NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
    W = NULL; U = NULL; dt = NULL; dUdt = NULL;
    dWdx = NULL; dWdy = NULL; dWdx_faceN = NULL; dWdy_faceN = NULL;
    dWdx_faceE = NULL; dWdy_faceE = NULL;
    dWdx_faceW = NULL; dWdy_faceW = NULL;
    dWdx_faceS = NULL; dWdy_faceS = NULL;
    phi = NULL;  Uo = NULL; Ut = NULL; Uold = NULL;
    FluxN = NULL; FluxS = NULL; FluxE = NULL; FluxW = NULL;
    WoN = NULL; WoS = NULL; WoE = NULL; WoW = NULL;
    Axisymmetric = 0; Gravity =0; Flow_Type = 0;
    Wall_Functions = 0; Freeze_Limiter = OFF;
    Moving_wall_velocity=ZERO; Pressure_gradient =ZERO; debug_level=0;
    Wall = NULL; 
  }
  
  LESPremixed2D_Quad_Block(const LESPremixed2D_Quad_Block &Soln) {
    NCi = Soln.NCi; ICl = Soln.ICl; ICu = Soln.ICu;
    NCj = Soln.NCj; JCl = Soln.JCl; JCu = Soln.JCu; Nghost = Soln.Nghost;
    Grid = Soln.Grid; W = Soln.W; U = Soln.U;
    dt = Soln.dt; 
    dUdt = Soln.dUdt; dWdx = Soln.dWdx; dWdy = Soln.dWdy;
    dWdx_faceN = Soln.dWdx_faceN; dWdy_faceN = Soln.dWdy_faceN;
    dWdx_faceE = Soln.dWdx_faceE; dWdy_faceE = Soln.dWdy_faceE;
    dWdx_faceW = Soln.dWdx_faceW; dWdy_faceW = Soln.dWdy_faceW;
    dWdx_faceS = Soln.dWdx_faceS; dWdy_faceS = Soln.dWdy_faceS;
    phi = Soln.phi; Uo = Soln.Uo; Ut = Soln.Ut; Uold = Soln.Uold; 
    FluxN = Soln.FluxN; FluxS = Soln.FluxS; FluxE = Soln.FluxE; FluxW = Soln.FluxW;
    WoN = Soln.WoN; WoS = Soln.WoS; WoE = Soln.WoE; WoW = Soln.WoW;
    Axisymmetric = Soln.Axisymmetric ; Gravity = Soln.Gravity; 
    Flow_Type = Soln.Flow_Type; 
    Wall_Functions = Soln.Wall_Functions; Freeze_Limiter = Soln.Freeze_Limiter;
    Pressure_gradient = Soln.Pressure_gradient;
    debug_level=0; Wall = Soln.Wall;
  }
   
  /* Assignment operator. */
  // LESPremixed2D_Quad_Block operator = (const LESPremixed2D_Quad_Block &Soln);
  // Use automatically generated assignment operator.  //THIS IS NOT GOOD DUE TO MEMORY ALLOCATION!!!!, POINTERS ONLY
  
  /* Allocate memory for structured quadrilateral solution block. */
  void allocate(const int Ni, const int Nj, const int Ng);
  
  /* Deallocate memory for structured quadrilateral solution block. */
  void deallocate(void);
  
  /* Return primitive solution state at specified node. */
  LESPremixed2D_pState Wn(const int &ii, const int &jj);
  
  /* Return conserved solution state at specified node. */
  LESPremixed2D_cState Un(const int &ii, const int &jj);

  /* Return conserved solution state at specified node. */
  LESPremixed2D_cState Uno(const int &ii, const int &jj);

  /* Return primitive solution state at cell nodes. */
  LESPremixed2D_pState WnNW(const int &ii, const int &jj);
  LESPremixed2D_pState WnNE(const int &ii, const int &jj);
  LESPremixed2D_pState WnSE(const int &ii, const int &jj);
  LESPremixed2D_pState WnSW(const int &ii, const int &jj);
  
  /* Return conserved solution state at cell nodes. */
  LESPremixed2D_cState UnNW(const int &ii, const int &jj);
  LESPremixed2D_cState UnNE(const int &ii, const int &jj);
  LESPremixed2D_cState UnSE(const int &ii, const int &jj);
  LESPremixed2D_cState UnSW(const int &ii, const int &jj);

  /* Return conserved solution state at cell nodes. */
  LESPremixed2D_cState UnoNW(const int &ii, const int &jj);
  LESPremixed2D_cState UnoNE(const int &ii, const int &jj);
  LESPremixed2D_cState UnoSE(const int &ii, const int &jj);
  LESPremixed2D_cState UnoSW(const int &ii, const int &jj);

  int BiLinearInterpolationCoefficients(double &eta, double &zeta, const int &ii, const int &jj);
  
  void set_v_zero(void);

  /*****************************************************************************
     dWn_dWc is the derivative of node solution w.r.t. cell center solution
     which is used by calculating the viscous Jacobians
  ******************************************************************************/
  double dWn_dWc (const int &i, const int &j, const int &Orient);


  /*****************************************************************************
     Vorticity and Laplacian of the vorticity related functions
  ******************************************************************************/

  /* Return vorticity at specified node. */
  double vorticity_n(const int &ii, const int &jj);
  double enstrophy_n(const int &ii, const int &jj);
  double vorticity(const int &i,const int &j);
  double enstrophy(const int &i,const int &j);        
    /****************************************************************************/   

  
  /* Set flags for limiter evaluation. */
  void evaluate_limiters(void) {Freeze_Limiter = OFF; } 
  void freeze_limiters(void) {Freeze_Limiter = ON; }
  
  /* Input-output operators. */
  friend ostream &operator << (ostream &out_file,
			       const LESPremixed2D_Quad_Block
			       &Soln);
  friend istream &operator >> (istream &in_file,
			       LESPremixed2D_Quad_Block
			       &Soln);

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
  /* Destructor. */
  // ~LESPremixed2D_Quad_Block(void);
  // Use automatically generated destructor
  // should maybe call the deallocate functions
  // from here ???
  //
  // Something like..... or will this mess with the AMR ???
  // Local_SolnBlk = Deallocate(Local_SolnBlk,Input_Parameters);
  // Deallocate_Message_Buffers(List_of_Local_Solution_Blocks); // Not necessary here!
  // List_of_Local_Solution_Blocks.deallocate();
  // List_of_Global_Solution_Blocks.deallocate();
  // QuadTree.deallocate();
  // MeshBlk = Deallocate_Multi_Block_Grid(MeshBlk,
  //					    Input_Parameters.Number_of_Blocks_Idir,
  //					    Input_Parameters.Number_of_Blocks_Jdir);
  //
};

/**************************************************************************
 * LESPremixed2D_Quad_Block::allocate -- Allocate memory.                 *
 **************************************************************************/
inline void LESPremixed2D_Quad_Block::allocate(const int Ni, const int Nj, const int Ng) {
  assert(Ni > 1 && Nj > 1); Grid.allocate(Ni, Nj, Ng);
   NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1;
   NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1; Nghost = Ng;
   W = new LESPremixed2D_pState*[NCi]; U = new LESPremixed2D_cState*[NCi];
   dt = new double*[NCi];
   dUdt = new LESPremixed2D_cState**[NCi];
   dWdx = new LESPremixed2D_pState*[NCi]; dWdy = new LESPremixed2D_pState*[NCi];
   dWdx_faceN = new LESPremixed2D_pState*[NCi]; dWdy_faceN = new LESPremixed2D_pState*[NCi];
   dWdx_faceE = new LESPremixed2D_pState*[NCi]; dWdy_faceE = new LESPremixed2D_pState*[NCi];
   dWdx_faceW = new LESPremixed2D_pState*[NCi]; dWdy_faceW = new LESPremixed2D_pState*[NCi];
   dWdx_faceS = new LESPremixed2D_pState*[NCi]; dWdy_faceS = new LESPremixed2D_pState*[NCi];
   phi = new LESPremixed2D_pState*[NCi]; Uo = new LESPremixed2D_cState*[NCi];
   Ut = new LESPremixed2D_cState*[NCi]; Uold = new LESPremixed2D_cState*[NCi];
   Wall = new Turbulent2DWallData*[NCi];
   
   for (int i = 0; i <= NCi-1 ; ++i ) {
      W[i] = new LESPremixed2D_pState[NCj]; U[i] = new LESPremixed2D_cState[NCj];
      dt[i] = new double[NCj]; dUdt[i] = new LESPremixed2D_cState*[NCj];
      for (int j = 0; j <= NCj-1 ; ++j ){
	dUdt[i][j] = new LESPremixed2D_cState[NUMBER_OF_RESIDUAL_VECTORS_LESPREMIXED2D];
      }
      dWdx[i] = new LESPremixed2D_pState[NCj]; dWdy[i] = new LESPremixed2D_pState[NCj];
      dWdx_faceN[i] = new LESPremixed2D_pState[NCj]; dWdy_faceN[i] = new LESPremixed2D_pState[NCj];
      dWdx_faceE[i] = new LESPremixed2D_pState[NCj]; dWdy_faceE[i] = new LESPremixed2D_pState[NCj];
      dWdx_faceW[i] = new LESPremixed2D_pState[NCj]; dWdy_faceW[i] = new LESPremixed2D_pState[NCj];
      dWdx_faceS[i] = new LESPremixed2D_pState[NCj]; dWdy_faceS[i] = new LESPremixed2D_pState[NCj];
      phi[i] = new LESPremixed2D_pState[NCj]; Uo[i] = new LESPremixed2D_cState[NCj];
      Ut[i] = new LESPremixed2D_cState[NCj]; Uold[i] = new LESPremixed2D_cState[NCj]; 
      Wall[i] = new Turbulent2DWallData[NCj];
   } /* endfor */
   FluxN = new LESPremixed2D_cState[NCi]; FluxS = new LESPremixed2D_cState[NCi];
   FluxE = new LESPremixed2D_cState[NCj]; FluxW = new LESPremixed2D_cState[NCj];
   WoN = new LESPremixed2D_pState[NCi]; WoS = new LESPremixed2D_pState[NCi];
   WoE = new LESPremixed2D_pState[NCj]; WoW = new LESPremixed2D_pState[NCj];

  // Set the solution residuals, gradients, limiters, and other values to zero.
   for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
      for ( int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
          for ( int k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_LESPREMIXED2D-1 ; ++k ) {
	     dUdt[i][j][k].Vacuum();
          } 
	  dWdx[i][j].Vacuum() ; dWdy[i][j].Vacuum();
	  phi[i][j].Vacuum(); Uo[i][j].Vacuum();
	  dt[i][j] = ZERO;
      } 
   }  
}

/**************************************************************************
 * LESPremixed2D_Quad_Block::deallocate -- Deallocate memory.             *
 **************************************************************************/
inline void LESPremixed2D_Quad_Block::deallocate(void) {
  int i, j; Grid.deallocate();
   for ( i = 0; i <= NCi-1 ; ++i ) {
      delete []W[i]; W[i] = NULL; delete []U[i]; U[i] = NULL;
      delete []dt[i]; dt[i] = NULL;
      for ( j = 0; j <= NCj-1 ; ++j ) { delete []dUdt[i][j]; dUdt[i][j] = NULL;
      }
      delete []dUdt[i]; dUdt[i] = NULL;
      delete []dWdx[i]; dWdx[i] = NULL; delete []dWdy[i]; dWdy[i] = NULL;
      delete []dWdx_faceN[i]; dWdx_faceN[i] = NULL; delete []dWdy_faceN[i]; dWdy_faceN[i] = NULL;
      delete []dWdx_faceE[i]; dWdx_faceE[i] = NULL; delete []dWdy_faceE[i]; dWdy_faceE[i] = NULL;
      delete []dWdx_faceW[i]; dWdx_faceW[i] = NULL; delete []dWdy_faceW[i]; dWdy_faceW[i] = NULL;
      delete []dWdx_faceS[i]; dWdx_faceS[i] = NULL; delete []dWdy_faceS[i]; dWdy_faceS[i] = NULL;
      delete []phi[i]; phi[i] = NULL; delete []Uo[i]; Uo[i] = NULL; 
      delete []Ut[i]; Ut[i] = NULL; delete []Uold[i]; Uold[i] = NULL;
      delete []Wall[i]; Wall[i] = NULL; 
   } /* endfor */
   delete []W; W = NULL; delete []U; U = NULL;
   delete []dt; dt = NULL; delete []dUdt; dUdt = NULL;
   delete []dWdx; dWdx = NULL; delete []dWdy; dWdy = NULL;
   delete []dWdx_faceN; dWdx_faceN = NULL; delete []dWdy_faceN; dWdy_faceN = NULL;
   delete []dWdx_faceE; dWdx_faceE = NULL; delete []dWdy_faceE; dWdy_faceE = NULL;
   delete []dWdx_faceW; dWdx_faceW = NULL; delete []dWdy_faceW; dWdy_faceW = NULL;
   delete []dWdx_faceS; dWdx_faceS = NULL; delete []dWdy_faceS; dWdy_faceS = NULL;

   delete []phi; phi = NULL; delete []Uo; Uo = NULL;  
   delete []Ut; Ut = NULL; delete []Uold; Uold = NULL;
   delete []Wall; Wall = NULL;
   delete []FluxN; FluxN = NULL; delete []FluxS; FluxS = NULL;
   delete []FluxE; FluxE = NULL; delete []FluxW; FluxW = NULL;
   delete []WoN; WoN = NULL; delete []WoS; WoS = NULL;
   delete []WoE; WoE = NULL; delete []WoW; WoW = NULL;
   NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
}


/**************************************************************************
 * LESPremixed2D_Quad_Block:BiLinearInterpolationCoefficients --          *
 *                     Coefficients used by bilinear interpolation        *
 **************************************************************************/
inline int LESPremixed2D_Quad_Block::BiLinearInterpolationCoefficients(double &eta, double &zeta, const int &ii, const int &jj){

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

    //DOES NOT CORRESPOND TO HOLMES BILINEAR INTERPOLATION SEE Wn, Un functions
}


/******************************************************************************
 * LESPremixed2D_Quad_Block::dWn_dWc -- Derivative of Node primitive solution *  
 *                                 w.r.t Cell primitive solution              *
 ******************************************************************************/
inline double LESPremixed2D_Quad_Block::dWn_dWc(const int &i, const int &j, const int &Orient) {
 
  double eta(ZERO), zeta(ZERO); 
  BiLinearInterpolationCoefficients(eta, zeta, i, j);

  switch(Orient) {
  case  NORTH_WEST:  
    return (eta - zeta*eta);
    break;
  case NORTH_EAST:
    return (ONE - zeta - eta + zeta*eta);
    break;
  case SOUTH_WEST:
    return (zeta*eta);
    break;
  case SOUTH_EAST:
    return (zeta - zeta*eta);
    break;
  default:
    cerr<<"\n Improper Orient in LESPremixed2D_Quad_Block::dWn_dWc\n";
    break;
  }
      
}

/**********************************************************
 * Vorticity:   curl(v)                                   *
 * Enstrophy:   0.5*sqr(vorticity)                        *
 **********************************************************/

inline double LESPremixed2D_Quad_Block::vorticity(const int &i,const int &j) {
  double dudy, dvdx, vort;
  dudy = dWdy[i][j].v.x;
  dvdx = dWdx[i][j].v.y;
  vort = dvdx - dudy;
 
  if(fabs(vort) < TOLER){
    return ZERO;
  } else {
    return vort;
  }  
}

inline double LESPremixed2D_Quad_Block::enstrophy(const int &i,const int &j) {
  double w;
  w = vorticity(i, j);

  return (0.5*w*w);
}
  
inline double LESPremixed2D_Quad_Block::vorticity_n(const int &ii, const int &jj) {
  double eta, zeta;
  eta = ZERO;
  zeta = ZERO;
  
  BiLinearInterpolationCoefficients(eta, zeta, ii, jj);
      
  return (vorticity(ii-1,jj-1) +(vorticity(ii-1,jj)-vorticity(ii-1,jj-1))*zeta+  
    	 (vorticity(ii,jj-1)-vorticity(ii-1,jj-1))*eta + 
    	 (vorticity(ii,jj)+vorticity(ii-1,jj-1)-vorticity(ii-1,jj)-vorticity(ii,jj-1))*zeta*eta);
}

inline double LESPremixed2D_Quad_Block::enstrophy_n(const int &ii,const int &jj) {
  double w;
  w = vorticity_n(ii, jj);

  return (0.5*w*w);
}




/**************************************************************************
 * LESPremixed2D_Quad_Block::Wn -- Node primitive solution.               *
 **************************************************************************/
inline LESPremixed2D_pState LESPremixed2D_Quad_Block::Wn(const int &ii, const int &jj) {
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y,
    eta1, zeta1, eta2, zeta2, eta, zeta;
  
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


  //THICKENING FACTOR CHECKS
/* #ifdef THICKENED_FLAME_ON */
/*   if (W[ii-1][jj-1].flame.TF<ONE) cout << "\nIn Soln.Wn, W[" << ii-1 << "][" <<jj-1 <<"] = " */
/* 				       << W[ii-1][jj-1].flame.TF; */

/*   if (W[ii-1][jj].flame.TF<ONE) cout << "\nIn Soln.Wn, W[" << ii-1 << "][" <<jj << "] = " */
/* 				     << W[ii-1][jj].flame.TF; */

/*   if (W[ii][jj-1].flame.TF<ONE) cout << "\nIn Soln.Wn, W[" << ii << "][" <<jj-1 << "] = " */
/* 				     << W[ii][jj-1].flame.TF; */

/*   if (W[ii][jj].flame.TF<ONE) cout << "\nIn Soln.Wn, W[" << ii << "][" << jj << "] = " */
/* 				   << W[ii][jj].flame.TF; */  
/* #endif */
  
  return (W[ii-1][jj-1] +(W[ii-1][jj]-W[ii-1][jj-1])*zeta+  
	  (W[ii][jj-1]-W[ii-1][jj-1])*eta + 
	  (W[ii][jj]+W[ii-1][jj-1]-W[ii-1][jj]-W[ii][jj-1])*zeta*eta);  
 

  //HOLMES et al.
//    double w1, w2, w3, w4;
//    LESPremixed2D_pState W1, W2, W3, W4;
//    Vector2D X0, X1, X2, X3, X4, lambda, R;
//    Tensor2D I;
//    // Summarize cell-centres and state:
//    X1 = Grid.Cell[ii-1][jj-1].Xc; W1 = W[ii-1][jj-1];
//    X2 = Grid.Cell[ii  ][jj-1].Xc; W2 = W[ii  ][jj-1];
//    X3 = Grid.Cell[ii-1][jj  ].Xc; W3 = W[ii-1][jj  ];
//    X4 = Grid.Cell[ii  ][jj  ].Xc; W4 = W[ii  ][jj  ];
//    X0 = Grid.Node[ii][jj].X;
//    // Determine weighting coefficients:
//    R = (X1 - X0) + (X2 - X0) + (X3 - X0) + (X4 - X0);
//    I.xx = (X1.x - X0.x)*(X1.x - X0.x) + (X2.x - X0.x)*(X2.x - X0.x) +
//       (X3.x - X0.x)*(X3.x - X0.x) + (X4.x - X0.x)*(X4.x - X0.x);
//    I.xy = (X1.x - X0.x)*(X1.y - X0.y) + (X2.x - X0.x)*(X2.y - X0.y) +
//       (X3.x - X0.x)*(X3.y - X0.y) + (X4.x - X0.x)*(X4.y - X0.y);
//    I.yy = (X1.y - X0.y)*(X1.y - X0.y) + (X2.y - X0.y)*(X2.y - X0.y) +
//       (X3.y - X0.y)*(X3.y - X0.y) + (X4.y - X0.y)*(X4.y - X0.y);
//    lambda.x = (I.xy*R.y - I.yy*R.x)/(I.xx*I.yy - I.xy*I.xy);
//    lambda.y = (I.xy*R.x - I.xx*R.y)/(I.xx*I.yy - I.xy*I.xy);
//    // Determine the weights:
//    w1 = 1 + lambda.x*(X1.x - X0.x) + lambda.y*(X1.y - X0.y);
//    w2 = 1 + lambda.x*(X2.x - X0.x) + lambda.y*(X2.y - X0.y);
//    w3 = 1 + lambda.x*(X3.x - X0.x) + lambda.y*(X3.y - X0.y);
//    w4 = 1 + lambda.x*(X4.x - X0.x) + lambda.y*(X4.y - X0.y);
//    // Return the interpolated state:
//    return (w1*W1 + w2*W2 + w3*W3+ w4*W4)/(w1 + w2 + w3 + w4);

}

/******************************************************************************
/******************************************************************************
 * LESPremixed2D_Quad_Block::Wn?? -- Get cell node primitive solution states. *
 ******************************************************************************/
inline LESPremixed2D_pState LESPremixed2D_Quad_Block::WnNW(const int &ii, const int &jj) {
  return (Wn(ii, jj+1));
}

inline LESPremixed2D_pState LESPremixed2D_Quad_Block::WnNE(const int &ii, const int &jj) {
  return (Wn(ii+1, jj+1));
}

inline LESPremixed2D_pState LESPremixed2D_Quad_Block::WnSE(const int &ii, const int &jj) {
  return (Wn(ii+1, jj));
}

inline LESPremixed2D_pState LESPremixed2D_Quad_Block::WnSW(const int &ii, const int &jj) {
  return (Wn(ii, jj));
}

/**************************************************************************
 * LESPremixed2D_Quad_Block::Un -- Node conservative solution.            *
 **************************************************************************/
inline LESPremixed2D_cState LESPremixed2D_Quad_Block::Un(const int &ii, const int &jj) {
 
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y,
    eta1, zeta1, eta2, zeta2, eta, zeta;
 
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
  
  return (U[ii-1][jj-1] +(U[ii-1][jj]-U[ii-1][jj-1])*zeta+  
	  (U[ii][jj-1]-U[ii-1][jj-1])*eta + 
	  (U[ii][jj]+U[ii-1][jj-1]-U[ii-1][jj]-U[ii][jj-1])*zeta*eta);
  
// /**********************************************************************
//  * LESPremixed2DQuad::Un -- Return the solution state                 *
//  *                              at the specified node using bilinear  *
//  *                              the interpolation presented by Holmes *
//  *                              and (AIAA Paper 1989-1932).           *
//  **********************************************************************/
//    double w1, w2, w3, w4;
//    LESPremixed2D_cState U1, U2, U3, U4;
//    Vector2D X0, X1, X2, X3, X4, lambda, R;
//    Tensor2D I;
//    // Summarize cell-centres and state:
//    X1 = Grid.Cell[ii-1][jj-1].Xc; U1 = U[ii-1][jj-1];
//    X2 = Grid.Cell[ii  ][jj-1].Xc; U2 = U[ii  ][jj-1];
//    X3 = Grid.Cell[ii-1][jj  ].Xc; U3 = U[ii-1][jj  ];
//    X4 = Grid.Cell[ii  ][jj  ].Xc; U4 = U[ii  ][jj  ];
//    X0 = Grid.Node[ii][jj].X;
//    // Determine weighting coefficients:
//    R = (X1 - X0) + (X2 - X0) + (X3 - X0) + (X4 - X0);
//    I.xx = (X1.x - X0.x)*(X1.x - X0.x) + (X2.x - X0.x)*(X2.x - X0.x) +
//       (X3.x - X0.x)*(X3.x - X0.x) + (X4.x - X0.x)*(X4.x - X0.x);
//    I.xy = (X1.x - X0.x)*(X1.y - X0.y) + (X2.x - X0.x)*(X2.y - X0.y) +
//       (X3.x - X0.x)*(X3.y - X0.y) + (X4.x - X0.x)*(X4.y - X0.y);
//    I.yy = (X1.y - X0.y)*(X1.y - X0.y) + (X2.y - X0.y)*(X2.y - X0.y) +
//       (X3.y - X0.y)*(X3.y - X0.y) + (X4.y - X0.y)*(X4.y - X0.y);
//    lambda.x = (I.xy*R.y - I.yy*R.x)/(I.xx*I.yy - I.xy*I.xy);
//    lambda.y = (I.xy*R.x - I.xx*R.y)/(I.xx*I.yy - I.xy*I.xy);
//    // Determine the weights:
//    w1 = 1 + lambda.x*(X1.x - X0.x) + lambda.y*(X1.y - X0.y);
//    w2 = 1 + lambda.x*(X2.x - X0.x) + lambda.y*(X2.y - X0.y);
//    w3 = 1 + lambda.x*(X3.x - X0.x) + lambda.y*(X3.y - X0.y);
//    w4 = 1 + lambda.x*(X4.x - X0.x) + lambda.y*(X4.y - X0.y);
//    // Return the interpolated state:
//    return (w1*U1 + w2*U2 + w3*U3+ w4*U4)/(w1 + w2 + w3 + w4);

}

/****************************************************************************
 * LESPremixed2D_Quad_Block::Un -- Get cell node conserved solution states. *
 ****************************************************************************/
inline LESPremixed2D_cState LESPremixed2D_Quad_Block::UnNW(const int &ii, const int &jj) {
  return (Un(ii, jj+1));
}

inline LESPremixed2D_cState LESPremixed2D_Quad_Block::UnNE(const int &ii, const int &jj) {
  return (Un(ii+1, jj+1));
}

inline LESPremixed2D_cState LESPremixed2D_Quad_Block::UnSE(const int &ii, const int &jj) {
  return (Un(ii+1, jj));
}

inline LESPremixed2D_cState LESPremixed2D_Quad_Block::UnSW(const int &ii, const int &jj) {
  return (Un(ii, jj));
}

/**************************************************************************
 * LESPremixed2D_Quad_Block::Un -- Node conservative solution.            *
 **************************************************************************/
inline LESPremixed2D_cState LESPremixed2D_Quad_Block::Uno(const int &ii, const int &jj) {
 
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y,
    eta1, zeta1, eta2, zeta2, eta, zeta;
 
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
  
  return (Uo[ii-1][jj-1] +(Uo[ii-1][jj]-Uo[ii-1][jj-1])*zeta+  
	  (Uo[ii][jj-1]-Uo[ii-1][jj-1])*eta + 
	  (Uo[ii][jj]+Uo[ii-1][jj-1]-Uo[ii-1][jj]-Uo[ii][jj-1])*zeta*eta);
}

/*****************************************************************************
 * LESPremixed2D_Quad_Block::Un -- Get cell node conserved solution states.  *
 *****************************************************************************/
inline LESPremixed2D_cState LESPremixed2D_Quad_Block::UnoNW(const int &ii, const int &jj) {
  return (Uno(ii, jj+1));
}

inline LESPremixed2D_cState LESPremixed2D_Quad_Block::UnoNE(const int &ii, const int &jj) {
  return (Uno(ii+1, jj+1));
}

inline LESPremixed2D_cState LESPremixed2D_Quad_Block::UnoSE(const int &ii, const int &jj) {
  return (Uno(ii+1, jj));
}

inline LESPremixed2D_cState LESPremixed2D_Quad_Block::UnoSW(const int &ii, const int &jj) {
  return (Uno(ii, jj));
}


/**************************************************************************
 * LESPremixed2D_Quad_Block -- Input-output operators.                    *
 **************************************************************************/
inline ostream &operator << (ostream &out_file,
			     const LESPremixed2D_Quad_Block &SolnBlk) {
  int i, j; 
  out_file << SolnBlk.Grid;
  out_file << SolnBlk.NCi << " " << SolnBlk.ICl << " " << SolnBlk.ICu << " " << SolnBlk.Nghost << "\n";
  out_file << SolnBlk.NCj << " " << SolnBlk.JCl << " " << SolnBlk.JCu << "\n";
  out_file << SolnBlk.Axisymmetric << "\n";
  out_file << SolnBlk.Flow_Type <<"\n";
  out_file << SolnBlk.Wall_Functions <<"\n";
  out_file << SolnBlk.Gravity <<"\n";
  out_file << SolnBlk.debug_level <<"\n";
  out_file << SolnBlk.Moving_wall_velocity <<"\n";
  out_file << SolnBlk.Pressure_gradient <<"\n";
  if (SolnBlk.NCi == 0 || SolnBlk.NCj == 0) return(out_file);
  for ( j  = SolnBlk.JCl-SolnBlk.Nghost ; j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
     for ( i = SolnBlk.ICl-SolnBlk.Nghost ; i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
         out_file << SolnBlk.U[i][j] << "\n"
	          << SolnBlk.Ut[i][j] << "\n"
		  << SolnBlk.Uold[i][j] << "\n";
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
			     LESPremixed2D_Quad_Block &SolnBlk) {

  int i, j, k, ni, il, iu, nj, jl, ju, ng;
  LESPremixed2D_pState LESPremixed2D_W_VACUUM(ZERO, Vector2D_ZERO, ZERO, ZERO);
  LESPremixed2D_cState LESPremixed2D_U_VACUUM(ZERO, Vector2D_ZERO, ZERO, ZERO);
  Grid2D_Quad_Block New_Grid; in_file >> New_Grid;
  in_file.setf(ios::skipws);
  in_file >> ni >> il >> iu >> ng; in_file >> nj >> jl >> ju;
  in_file >> SolnBlk.Axisymmetric;
  in_file >> SolnBlk.Flow_Type;
  in_file >> SolnBlk.Wall_Functions;
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
         in_file >> SolnBlk.U[i][j] >> SolnBlk.Ut[i][j] >> SolnBlk.Uold[i][j];
         SolnBlk.W[i][j] = W(SolnBlk.U[i][j]);
	 SolnBlk.Ut[i][j] = SolnBlk.U[i][j]; 
         for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_LESPREMIXED2D-1 ; ++k ) {
	     SolnBlk.dUdt[i][j][k] = LESPremixed2D_U_VACUUM;
         } /* endfor */
	 SolnBlk.dWdx[i][j] = LESPremixed2D_W_VACUUM;
	 SolnBlk.dWdy[i][j] = LESPremixed2D_W_VACUUM;
	 SolnBlk.dWdx_faceN[i][j] = LESPremixed2D_W_VACUUM;
	 SolnBlk.dWdy_faceN[i][j] = LESPremixed2D_W_VACUUM;
	 SolnBlk.dWdx_faceE[i][j] = LESPremixed2D_W_VACUUM;
	 SolnBlk.dWdy_faceE[i][j] = LESPremixed2D_W_VACUUM;
	 SolnBlk.dWdx_faceW[i][j] = LESPremixed2D_W_VACUUM;
	 SolnBlk.dWdy_faceW[i][j] = LESPremixed2D_W_VACUUM;
	 SolnBlk.dWdx_faceS[i][j] = LESPremixed2D_W_VACUUM;
	 SolnBlk.dWdy_faceS[i][j] = LESPremixed2D_W_VACUUM;
	 SolnBlk.phi[i][j] = LESPremixed2D_W_VACUUM;
	 SolnBlk.Uo[i][j] = LESPremixed2D_U_VACUUM;
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
 * LESPremixed2D_Quad_Block::NumVar -- Returns number of state variables.      *
 *******************************************************************************/
inline int LESPremixed2D_Quad_Block::NumVar(void) {
#ifdef THICKENED_FLAME_ON
  return (int(W[0][0].NUM_VAR_LESPREMIXED2D + 2));
#else   
  return (int(W[0][0].NUM_VAR_LESPREMIXED2D));    
#endif
}

/*******************************************************************************
 * LESPremixed2D_Quad_Block::LoadSendBuffer -- Loads send message buffer.      *
 *******************************************************************************/
inline int LESPremixed2D_Quad_Block::LoadSendBuffer(double *buffer,
                                              int &buffer_count,
                                              const int buffer_size,
                                              const int i_min, 
                                              const int i_max,
                                              const int i_inc,
                                              const int j_min, 
                                              const int j_max,
                                              const int j_inc) {
  int NUM_VAR_LESPREMIXED2D = NumVar();
  int i, j, k;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = U[i][j][k];
        } /* endfor */
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * LESPremixed2D_Quad_Block::LoadSendBuffer_F2C -- Loads send message buffer   *
 *                                           for fine to coarse block message  *
 *                                           passing.                          *
 *******************************************************************************/
inline int LESPremixed2D_Quad_Block::LoadSendBuffer_F2C(double *buffer,
                                                  int &buffer_count,
                                                  const int buffer_size,
                                                  const int i_min, 
                                                  const int i_max,
                                                  const int i_inc,
                                                  const int j_min, 
                                                  const int j_max,
                                                  const int j_inc) {
  int NUM_VAR_LESPREMIXED2D = NumVar();
  int i, j, k;
  for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
 * LESPremixed2D_Quad_Block::LoadSendBuffer_C2F -- Loads send message buffer   *
 *                                           for coarse to fine block message  *
 *                                           passing.                          *
 *******************************************************************************/
inline int LESPremixed2D_Quad_Block::LoadSendBuffer_C2F(double *buffer,
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
  int NUM_VAR_LESPREMIXED2D = NumVar();
  int i, j, k;
  Vector2D dX;
  LESPremixed2D_pState Wfine;
  LESPremixed2D_cState Ufine;

  if (j_min == j_max) { // North or south boundary.
     // Four different orderings to consider depending on the value of i_inc & j_inc.
     if (j_inc > 0) {
        if (i_inc > 0) {
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Perform limited linear least squares reconstruction in cell (i, j_min).
              SubcellReconstruction(i, j_min, LIMITER_VENKATAKRISHNAN);
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Ufine[k];
              } /* endfor */
           } /* endfor */
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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

/************************************************************************************
 * LESPremixed2D_Quad_Block::UnloadReceiveBuffer -- Unloads receive message buffer. *
 ************************************************************************************/
inline int LESPremixed2D_Quad_Block::UnloadReceiveBuffer(double *buffer,
                                                   int &buffer_count,
                                                   const int buffer_size,
                                                   const int i_min, 
                                                   const int i_max,
                                                   const int i_inc,
                                                   const int j_min, 
                                                   const int j_max,
                                                   const int j_inc) { 
  int NUM_VAR_LESPREMIXED2D = NumVar();
 
  for ( int j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
    for ( int i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      //Changed from Euler2D
      for ( int k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	U[i][j][k] = buffer[buffer_count];
      }
      W[i][j] = U[i][j].W();
    } /* endfor */
  } /* endfor */
  return(0);
}

/********************************************************************************
 * LESPremixed2D_Quad_Block::UnloadReceiveBuffer_F2C -- Unloads receive message *
 *                                                buffer for fine to coarse     *
 *                                                block message passing.        *
 ********************************************************************************/
inline int LESPremixed2D_Quad_Block::UnloadReceiveBuffer_F2C(double *buffer,
                                                       int &buffer_count,
                                                       const int buffer_size,
                                                       const int i_min, 
                                                       const int i_max,
                                                       const int i_inc,
                                                       const int j_min, 
                                                       const int j_max,
                                                       const int j_inc) { 
  int NUM_VAR_LESPREMIXED2D = NumVar();
  int i, j;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
       //Changed from Euler2D
       for ( int k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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

/********************************************************************************
 * LESPremixed2D_Quad_Block::UnloadReceiveBuffer_C2F -- Unloads receive message *
 *                                                buffer for coarse to fine     *
 *                                                block message passing.        *
 ********************************************************************************/
inline int LESPremixed2D_Quad_Block::UnloadReceiveBuffer_C2F(double *buffer,
                                                       int &buffer_count,
                                                       const int buffer_size,
                                                       const int i_min, 
                                                       const int i_max,
                                                       const int i_inc,
                                                       const int j_min, 
                                                       const int j_max,
                                                       const int j_inc) { 
  int NUM_VAR_LESPREMIXED2D = NumVar();
  int i, j;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
       //Changed from Euler2D
       for ( int k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
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
 * LESPremixed2D_Quad_Block::SubcellReconstruction --                     *
 *               Performs the subcell reconstruction of solution state    *
 *               within a given cell (i,j) of the computational mesh for  *
 *               the specified quadrilateral solution block.              *
 **************************************************************************/
inline void LESPremixed2D_Quad_Block::SubcellReconstruction(const int i, 
							    const int j,
							    const int Limiter) {

  int n, n2, n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4], phi_n;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  LESPremixed2D_pState DU, DUDx_ave, DUDy_ave;

  LESPremixed2D_pState LESPremixed2D_W_VACUUM;
  LESPremixed2D_W_VACUUM.Vacuum();
 
#ifdef THICKENED_FLAME_ON    
    int NUM_VAR_LESPREMIXED2D = NumVar() - 2;
#else
    int NUM_VAR_LESPREMIXED2D = NumVar();
#endif 
 
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
      DUDx_ave = LESPremixed2D_W_VACUUM;
      DUDy_ave = LESPremixed2D_W_VACUUM;
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
	for ( n = 1 ; n <= NUM_VAR_LESPREMIXED2D ; ++n ) {
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
      dWdx[i][j] = LESPremixed2D_W_VACUUM;
      dWdy[i][j] = LESPremixed2D_W_VACUUM; 
      phi[i][j]  = LESPremixed2D_W_VACUUM;
  } /* endif */
    
}


/***********************************************************************************
 * LESPremixed2D_Quad_Block::LoadSendBuffer_Flux_F2C -- Loads send message buffer  *
 *                                                for fine to coarse block message *
 *                                                passing of conservative          *
 *                                                solution fluxes.                 *
 ***********************************************************************************/
inline int LESPremixed2D_Quad_Block::LoadSendBuffer_Flux_F2C(double *buffer,
                                                       int &buffer_count,
                                                       const int buffer_size,
                                                       const int i_min, 
                                                       const int i_max,
                                                       const int i_inc,
                                                       const int j_min, 
                                                       const int j_max,
                                                       const int j_inc) { 
  int NUM_VAR_LESPREMIXED2D = NumVar(); 
  int i, j, k;
  if (j_min == j_max && j_min == JCl) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxS[i  ][k]+
                                   FluxS[i+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxN[i  ][k]+
                                   FluxN[i+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
     for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
        for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxW[j][k]+
                                   FluxW[j+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
     for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
        for ( k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxE[j][k]+
                                   FluxE[j+1][k]);
        } /* endfor */
     } /* endfor */
  } /* endif */
  return(0);
}

/*************************************************************************************
 * LESPremixed2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C -- Unloads receive message *
 *                                                buffer for fine to coarse          *
 *                                                block message passing of           *
 *                                                conservative solution fluxes.      *
 *************************************************************************************/
inline int LESPremixed2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C(double *buffer,
                                                            int &buffer_count,
                                                            const int buffer_size,
                                                            const int i_min, 
                                                            const int i_max,
                                                            const int i_inc,
                                                            const int j_min, 
                                                            const int j_max,
                                                            const int j_inc) {
  int NUM_VAR_LESPREMIXED2D = NumVar();
  if (j_min == j_max && j_min == JCl) {
    for ( int i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      for ( int k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	FluxS[i][k] = - buffer[buffer_count] - FluxS[i][k];
      }
    }
  } else if (j_min == j_max && j_min == JCu) {
    for ( int i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      for ( int k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	FluxN[i][k] = - buffer[buffer_count] - FluxN[i][k];
      }
    } 
  } else if (i_min == i_max && i_min == ICl) {
    for ( int j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
      for ( int k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	FluxW[j][k] = - buffer[buffer_count] - FluxW[j][k];
      }
    } 
  } else if (i_min == i_max && i_min == ICu) {
    for ( int j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
      for ( int k = 1 ; k <= NUM_VAR_LESPREMIXED2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	FluxE[j][k] = - buffer[buffer_count] - FluxE[j][k];
      }
    } 
  } /* endif */
  return(0);
}

/**************************************************************************
 * LESPremixed2D_Quad_Block -- Single Block External Subroutines.         * 
 * LESPremixed2DQuadSingleBlock.cc                                        *
 **************************************************************************/

extern void Write_Solution_Block(LESPremixed2D_Quad_Block &SolnBlk,
	                         ostream &Out_File);

extern void Read_Solution_Block(LESPremixed2D_Quad_Block &SolnBlk,
	                        istream &In_File);

extern void Broadcast_Solution_Block(LESPremixed2D_Quad_Block &SolnBlk);

#ifdef _MPI_VERSION
extern void Broadcast_Solution_Block(LESPremixed2D_Quad_Block &SolnBlk,
                                     MPI::Intracomm &Communicator, 
                                     const int Source_CPU);
#endif

extern void Copy_Solution_Block(LESPremixed2D_Quad_Block &SolnBlk1,
		                LESPremixed2D_Quad_Block &SolnBlk2);

extern int Prolong_Solution_Block(LESPremixed2D_Quad_Block &SolnBlk_Fine,
				  LESPremixed2D_Quad_Block &SolnBlk_Original,
				  const int Sector);

extern int Restrict_Solution_Block(LESPremixed2D_Quad_Block &SolnBlk_Coarse,
				   LESPremixed2D_Quad_Block &SolnBlk_Original_SW,
				   LESPremixed2D_Quad_Block &SolnBlk_Original_SE,
				   LESPremixed2D_Quad_Block &SolnBlk_Original_NW,
				   LESPremixed2D_Quad_Block &SolnBlk_Original_NE);

extern void Output_Tecplot(LESPremixed2D_Quad_Block &SolnBlk,
			   LESPremixed2D_Input_Parameters &IP,
			   const int Number_of_Time_Steps,
                           const double &Time,
                           const int Block_Number,
                           const int Output_Title,
	                   ostream &Out_File);

extern void Output_Cells_Tecplot(LESPremixed2D_Quad_Block &SolnBlk,
				 LESPremixed2D_Input_Parameters &IP,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_Nodes_Tecplot(LESPremixed2D_Quad_Block &SolnBlk,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_RHS(LESPremixed2D_Quad_Block &SolnBlk,
		       const int Number_of_Time_Steps,
		       const double &Time,
		       const int Block_Number,
		       const int Output_Title,
		       ostream &Out_File);


extern void ICs(LESPremixed2D_Quad_Block &SolnBlk,
 	        const int i_ICtype,
                LESPremixed2D_pState *Wo, 
		LESPremixed2D_Input_Parameters &Input_Parameters,
		const int Block_Number);

extern void Reset_Wo(LESPremixed2D_Quad_Block &SolnBlk );

extern void BCs(LESPremixed2D_Quad_Block &SolnBlk, 
		LESPremixed2D_Input_Parameters &IP);

extern double CFL(LESPremixed2D_Quad_Block &SolnBlk,
                  LESPremixed2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(LESPremixed2D_Quad_Block &SolnBlk, 
                                const double &Dt_min);

extern double L1_Norm_Residual(LESPremixed2D_Quad_Block &SolnBlk, const int &norm);

extern double L2_Norm_Residual(LESPremixed2D_Quad_Block &SolnBlk, const int &norm);

extern double Max_Norm_Residual(LESPremixed2D_Quad_Block &SolnBlk, const int &norm);


/**************** Reconstruction Methods ****************************/

extern void Linear_Reconstruction_LeastSquares(LESPremixed2D_Quad_Block &SolnBlk,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares(LESPremixed2D_Quad_Block &SolnBlk,
                                               const int i,
                                               const int j,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares_2(LESPremixed2D_Quad_Block &SolnBlk,
                                                 const int i,
                                                 const int j,
					         const int Limiter);

extern void Linear_Reconstruction_GreenGauss(LESPremixed2D_Quad_Block &SolnBlk,
					     const int Limiter);

extern void Linear_Reconstruction_GreenGauss(LESPremixed2D_Quad_Block &SolnBlk,
                                             const int i,
                                             const int j,
					     const int Limiter);

/******* Diamond Path Reconstruction Methods *************************/
extern void Linear_Reconstruction_LeastSquares_Diamond(LESPremixed2D_Quad_Block &SolnBlk,
						       const int Limiter);

extern void Linear_Reconstruction_LeastSquares_Diamond(LESPremixed2D_Quad_Block &SolnBlk,
						       const int i, 
						       const int j,
						       const int Limiter);


extern void Linear_Reconstruction_GreenGauss_Diamond(LESPremixed2D_Quad_Block &SolnBlk,
						     const int Limiter);

extern void Linear_Reconstruction_GreenGauss_Diamond(LESPremixed2D_Quad_Block &SolnBlk,
						     const int i,
						     const int j,
						     const int Limiter);

extern double Laplacian_of_Vorticity(LESPremixed2D_Quad_Block &SolnBlk,
                                     const int i,const int j);

extern void Residual_Smoothing(LESPremixed2D_Quad_Block &SolnBlk,
                               const int k_residual,
			       double &epsilon, 
                               const int number_of_Gauss_Seidel_iterations);

extern void Calculate_Refinement_Criteria(double *refinement_criteria,
					  LESPremixed2D_Input_Parameters &IP, 
                                          int &number_refinement_criteria,					  
                                          LESPremixed2D_Quad_Block &SolnBlk);

extern void Fix_Refined_Block_Boundaries(LESPremixed2D_Quad_Block &SolnBlk,
                                         const int Fix_North_Boundary,
                                         const int Fix_South_Boundary,
                                         const int Fix_East_Boundary,
                                         const int Fix_West_Boundary);

extern void Unfix_Refined_Block_Boundaries(LESPremixed2D_Quad_Block &SolnBlk);

extern void Apply_Boundary_Flux_Corrections(LESPremixed2D_Quad_Block &SolnBlk,
                                            const int Number_Neighbours_North_Boundary,
                                            const int Number_Neighbours_South_Boundary,
                                            const int Number_Neighbours_East_Boundary,
                                            const int Number_Neighbours_West_Boundary);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(LESPremixed2D_Quad_Block &SolnBlk,
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

extern int dUdt_Residual_Evaluation(LESPremixed2D_Quad_Block &SolnBlk,
				    LESPremixed2D_Input_Parameters &Input_Parameters);


extern int dUdt_Multistage_Explicit(LESPremixed2D_Quad_Block &SolnBlk,
   	                            const int i_stage,
                                    LESPremixed2D_Input_Parameters &Input_Parameters);

extern int Update_Solution_Multistage_Explicit(LESPremixed2D_Quad_Block &SolnBlk,
   	                                       const int i_stage,
                                               LESPremixed2D_Input_Parameters &Input_Parameters);

extern int Update_Dual_Solution_States(LESPremixed2D_Quad_Block &SolnBlk); 

extern void Viscous_Calculations(LESPremixed2D_Quad_Block &SolnBlk);

/**************************************************************************
 * LESPremixed2D_Quad_Block -- Multiple Block External Subroutines.       *
 * LESPremixed2DQuadMultiBlock.cc                                         *
 **************************************************************************/

extern LESPremixed2D_Quad_Block* Allocate(LESPremixed2D_Quad_Block *Soln_ptr,
					  LESPremixed2D_Input_Parameters &Input_Parameters);

extern LESPremixed2D_Quad_Block* Deallocate(LESPremixed2D_Quad_Block *Soln_ptr,
					    LESPremixed2D_Input_Parameters &Input_Parameters);

extern void ICs(LESPremixed2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
                LESPremixed2D_Input_Parameters &Input_Parameters);
  
extern int Read_Restart_Solution(LESPremixed2D_Quad_Block *Soln_ptr,
                                 AdaptiveBlock2D_List &Soln_Block_List,
                                 LESPremixed2D_Input_Parameters &Input_Parameters,
		                 int &Number_of_Time_Steps,
                                 double &Time,
                                 CPUTime &CPU_Time);

extern int Write_Restart_Solution(LESPremixed2D_Quad_Block *Soln_ptr,
                                  AdaptiveBlock2D_List &Soln_Block_List,
                                  LESPremixed2D_Input_Parameters &Input_Parameters,
		                  const int Number_of_Time_Steps,
                                  const double &Time,
                                  const CPUTime &CPU_Time);

extern int Output_Tecplot(LESPremixed2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          LESPremixed2D_Input_Parameters &Input_Parameters,
		          const int Number_of_Time_Steps,
                          const double &Time);

extern int Output_RHS(LESPremixed2D_Quad_Block *Soln_ptr,
		      AdaptiveBlock2D_List &Soln_Block_List,
		      LESPremixed2D_Input_Parameters &Input_Parameters,
		      const int Number_of_Time_Steps,
		      const double &Time);

int Output_PERTURB(LESPremixed2D_Quad_Block *Soln_ptr,
		   AdaptiveBlock2D_List &Soln_Block_List,
		   LESPremixed2D_Input_Parameters &Input_Parameters,
		   const int Number_of_Time_Steps,
		   const double &Time,
		   const CPUTime &CPU_Time);

extern int Output_Cells_Tecplot(LESPremixed2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                LESPremixed2D_Input_Parameters &Input_Parameters,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Mesh_Tecplot(LESPremixed2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               LESPremixed2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Mesh_Gnuplot(LESPremixed2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               LESPremixed2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Ringleb(LESPremixed2D_Quad_Block *Soln_ptr,
			  AdaptiveBlock2D_List &Soln_Block_List,
			  LESPremixed2D_Input_Parameters &IP);

extern int Output_Viscous_Channel(LESPremixed2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  LESPremixed2D_Input_Parameters &IP);

extern int Output_Flat_Plate(LESPremixed2D_Quad_Block *Soln_ptr,
			     AdaptiveBlock2D_List &Soln_Block_List,
			     LESPremixed2D_Input_Parameters &IP);

extern int Output_Driven_Cavity_Flow(LESPremixed2D_Quad_Block *Soln_ptr,
				     AdaptiveBlock2D_List &Soln_Block_List,
				     LESPremixed2D_Input_Parameters &IP);

extern void BCs(LESPremixed2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
		LESPremixed2D_Input_Parameters &Input_Parameters);

extern double CFL(LESPremixed2D_Quad_Block *Soln_ptr,
                  AdaptiveBlock2D_List &Soln_Block_List,
                  LESPremixed2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(LESPremixed2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List, 
                                const double &Dt_min);

extern double L1_Norm_Residual(LESPremixed2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);
		
extern double L2_Norm_Residual(LESPremixed2D_Quad_Block *Soln_ptr ,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double Max_Norm_Residual(LESPremixed2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List);

extern void L1_Norm_Residual(LESPremixed2D_Quad_Block *Soln_ptr, 
                               AdaptiveBlock2D_List &Soln_Block_List,
			       double *l1_norm);

extern void L2_Norm_Residual(LESPremixed2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
			       double *l2_norm);

extern void Max_Norm_Residual(LESPremixed2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
				double *max_norm);

extern void Evaluate_Limiters(LESPremixed2D_Quad_Block *Soln_ptr,
                              AdaptiveBlock2D_List &Soln_Block_List);

extern void Freeze_Limiters(LESPremixed2D_Quad_Block *Soln_ptr,
                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Change_Mref(LESPremixed2D_Quad_Block *Soln_ptr,
			AdaptiveBlock2D_List &Soln_Block_List,
			const double &Mr);

extern void Residual_Smoothing(LESPremixed2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               LESPremixed2D_Input_Parameters &Input_Parameters,
   	                       const int I_Stage);

extern void Apply_Boundary_Flux_Corrections(LESPremixed2D_Quad_Block *Soln_ptr,
                                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(LESPremixed2D_Quad_Block *Soln_ptr,
                                                                AdaptiveBlock2D_List &Soln_Block_List,
                                                                LESPremixed2D_Input_Parameters &Input_Parameters,
   	                                                        const int I_Stage);

extern int dUdt_Multistage_Explicit(LESPremixed2D_Quad_Block *Soln_ptr,
				    AdaptiveBlockResourceList &Global_Soln_Block_List,
                                    AdaptiveBlock2D_List &Soln_Block_List,
                                    LESPremixed2D_Input_Parameters &Input_Parameters,
   	                            const int I_Stage);

extern int Update_Solution_Multistage_Explicit(LESPremixed2D_Quad_Block *Soln_ptr,
                                               AdaptiveBlock2D_List &Soln_Block_List,
                                               LESPremixed2D_Input_Parameters &Input_Parameters,
   	                                       const int I_Stage);

extern int Update_Dual_Solution_States(LESPremixed2D_Quad_Block *Soln_ptr,
                                        AdaptiveBlock2D_List &Soln_Block_List); 



//************ Viscous Stuff************************************/
extern void Viscous_Calculations(LESPremixed2D_Quad_Block &SolnBlk);

// extern LESPremixed2D_cState Flux_Viscous_x(LESPremixed2D_Quad_Block &SolnBlk, const int &i, const int &j);
// extern LESPremixed2D_cState Flux_Viscous_y(LESPremixed2D_Quad_Block &SolnBlk, const int &i, const int &j);
// extern LESPremixed2D_cState Flux_Viscous_n(LESPremixed2D_Quad_Block &SolnBlk, const int &Orient, const int &i, const int &j);
// extern LESPremixed2D_cState Flux_Viscous_EastFace_Diamond(LESPremixed2D_Quad_Block &SolnBlk, const int &i, const int &j);
// extern LESPremixed2D_cState Flux_Viscous_NorthFace_Diamond(LESPremixed2D_Quad_Block &SolnBlk, const int &i, const int &j);


// extern double  Wall_Distance(LESPremixed2D_Quad_Block &SolnBlk, LESPremixed2D_Input_Parameters &Input_Parameters, int i, int j);
extern double Distance_to_Wall(LESPremixed2D_Quad_Block &SolnBlk, const Vector2D X_cell);

/******************************************************************************
 * LESPremixed2D_Quad_Block -- Multiple Block External Subroutines for Mesh.  * 
 * LESPremixed2DQuadGrid.cc                                                   *
 ******************************************************************************/

extern Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                            LESPremixed2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Broadcast_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                      LESPremixed2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                             LESPremixed2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                                            LESPremixed2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                  LESPremixed2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                 LESPremixed2D_Input_Parameters &Input_Parameters);

extern int Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                          LESPremixed2D_Input_Parameters &Input_Parameters);

extern int Output_Nodes_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                LESPremixed2D_Input_Parameters &Input_Parameters);

extern int Output_Cells_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                LESPremixed2D_Input_Parameters &Input_Parameters);


/**************************************************************************
 * LESPremixed2D_Quad_Block -- Solvers.                                   *
 * LESPremixed2DQuadSolvers.cc                                            *
 **************************************************************************/

extern int LESPremixed2DQuadSolver(char *Input_File_Name_ptr, int batch_flag);

/**************************************************************************
 * LESPremixed2D_Tools -- Tools and stuff                                 *
 * LESPremixed2DTools.cc                                                  *
 **************************************************************************/

extern int Open_Time_Accurate_File(ofstream &Time_Accurate_File,
				   char *File_Name,
				   const int Append_to_Fileconst,
				   const LESPremixed2D_pState &Soln);

extern int Close_Time_Accurate_File(ofstream &Time_Accurate_File);

extern void Output_to_Time_Accurate_File(ostream &Time_Accurate_File,
					 const double &Time,
					 const LESPremixed2D_pState &Soln);
  
extern void Output_Ringleb_Solution(LESPremixed2D_Quad_Block &SolnBlk,
				    const int Block_Number,
				    const int Output_Title,
				    ostream &Out_File);

extern void Output_Ringleb_Error(LESPremixed2D_Quad_Block &SolnBlk,
				 double &l1_norm,
				 double &l2_norm,
				 double &max_norm,
				 int &numberofactivecells);

extern void Output_Viscous_Channel(LESPremixed2D_Quad_Block &SolnBlk,
				   const int Block_Number,
				   const int Output_Title,
				   ostream &Out_File,
				   double &l1_norm,
				   double &l2_norm,
				   double &max_norm,
				   double &Vwall,
				   const double dp);

extern void Output_Flat_Plate(LESPremixed2D_Quad_Block &SolnBlk,
			      const int Block_Number,
			      const int Output_Title_Soln,
			      ostream &Out_File_Soln,
			      const int Output_Title_Skin,
			      ostream &Out_File_Skin,
			      const LESPremixed2D_pState &Winf,
			      double &l1_norm,
			      double &l2_norm,
			      double &max_norm);
			      
extern void Output_Driven_Cavity_Flow(LESPremixed2D_Quad_Block &SolnBlk,
				      const int Block_Number,
				      const int Output_Title,
				      ostream &Out_File_u,
				      ostream &Out_File_v,
				      const double &Re,
				      const double &Vwall,
				      const double &length);

extern int Output_Quasi3D_Tecplot(LESPremixed2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  LESPremixed2D_Input_Parameters &IP,
				  const int Number_of_Time_Steps,
				  const double &Time);

extern void Output_Quasi3D_Tecplot(LESPremixed2D_Quad_Block &SolnBlk,
				   LESPremixed2D_Input_Parameters &IP,
				   const int Number_of_Time_Steps,
				   const double &Time,
				   const int Block_Number,
				   const int Output_Title,
				   ostream &Out_File);


extern int Open_Energy_File(ofstream &Energy_File,
			    char *File_Name,
			    const int Append_to_File);

extern int Close_Energy_File(ofstream &Energy_File);

extern void Output_Energy_to_File(ostream &Energy_File,
				  const int Number_of_Time_Steps,
				  const double &Time,
				  const CPUTime &CPU_Time,
				  const double &Total_Energy,
				  const double &Total_Enstrophy,
                                  const int &n_inner);

extern int Open_Turbulence_Progress_File(ofstream &Turbulence_Progress_File,
					 char *File_Name,
					 const int Append_to_File);

extern int Close_Turbulence_Progress_File(ofstream &Turbulence_Progress_File);

extern void Output_Turbulence_Progress_to_File(ostream &Turbulence_Progress_File,
					       const int Number_of_Time_Steps,
					       const double &Time,
					       const CPUTime &CPU_Time,
					       const double &u_rms,
					       const double &Taylor_scale,
					       const double &viscosity,
					       const double &turbulent_burning_rate,
                                               const int &n_inner); 

extern void Write_Turbulent_Solution(LESPremixed2D_Quad_Block *Soln_ptr,
				     AdaptiveBlock2D_List &Soln_Block_List,
				     LESPremixed2D_Input_Parameters &Input_Parameters);

extern void Write_Slice_Solution(LESPremixed2D_Quad_Block *Soln_ptr,
				 AdaptiveBlock2D_List &Soln_Block_List,
				 LESPremixed2D_Input_Parameters &Input_Parameters,
				 int j_slice);

extern void Write_Wrinkled_Flame(LESPremixed2D_Quad_Block *Soln_ptr,
				 AdaptiveBlock2D_List &Soln_Block_List,
				 LESPremixed2D_Input_Parameters &Input_Parameters);

extern double Max_Grad_Temp(LESPremixed2D_Quad_Block *Soln_ptr,
			    AdaptiveBlock2D_List &Soln_Block_List);

extern double Max_Sensor(LESPremixed2D_Quad_Block *Soln_ptr,
			 AdaptiveBlock2D_List &Soln_Block_List);


/*****************************************************************
 *  Turbulence related functions  -> Single block                *
 *  LESPremixed2DTurbInit.cc                                     *
 *****************************************************************/  
extern double Total_Block_Area(LESPremixed2D_Quad_Block &SolnBlk);


/*****************************************************************
 *  Turbulence related functions  -> Multiple block              *
 *  LESPremixed2DTurbInit.cc                                     *
 *****************************************************************/ 
extern double Total_Area(LESPremixed2D_Quad_Block *Soln_ptr,
                         AdaptiveBlock2D_List &Soln_Block_List);

//extern double Energy_Spectrum(double abs_wave_num, const int &spectrum_flag);

extern double Max_k(LESPremixed2D_Input_Parameters &Input_Parameters);

extern double Average_viscosity(LESPremixed2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               LESPremixed2D_Input_Parameters &Input_Parameters);

extern void Average_V(LESPremixed2D_Quad_Block *Soln_ptr,
                      AdaptiveBlock2D_List &Soln_Block_List,
                      LESPremixed2D_Input_Parameters &Input_Parameters,
                      double &u_average, double &v_average);

extern double Total_TKE(LESPremixed2D_Quad_Block *Soln_ptr,
                        AdaptiveBlock2D_List &Soln_Block_List,
                        LESPremixed2D_Input_Parameters &Input_Parameters);

extern double Total_Enstrophy(LESPremixed2D_Quad_Block *Soln_ptr,
			      AdaptiveBlock2D_List &Soln_Block_List,
			      LESPremixed2D_Input_Parameters &Input_Parameters);

extern double u_rms(LESPremixed2D_Quad_Block *Soln_ptr,
		    AdaptiveBlock2D_List &Soln_Block_List,
		    LESPremixed2D_Input_Parameters &Input_Parameters);

extern double Taylor_Scale(LESPremixed2D_Quad_Block *Soln_ptr,
			   AdaptiveBlock2D_List &Soln_Block_List,
			   LESPremixed2D_Input_Parameters &Input_Parameters);

extern double Turbulent_Burning_Rate(LESPremixed2D_Quad_Block *Soln_ptr,
				     AdaptiveBlock2D_List &Soln_Block_List,
				     LESPremixed2D_Input_Parameters &Input_Parameters);

extern double Total_Species_Mass(LESPremixed2D_Quad_Block *Soln_ptr,
				 AdaptiveBlock2D_List &Soln_Block_List,
				 const int &k);

extern double Average_Flame_Position(LESPremixed2D_Quad_Block *Soln_ptr,
				     AdaptiveBlock2D_List &Soln_Block_List,
				     LESPremixed2D_Input_Parameters &Input_Parameters);

extern void Average(LESPremixed2D_Quad_Block *Soln_ptr,
                    AdaptiveBlock2D_List &Soln_Block_List,
                    LESPremixed2D_Input_Parameters &Input_Parameters,
                    const double &u_average, const double &v_average,
                    double &sqr_u);

extern void Velocity_Fluctuations(Grid2D_Quad_Block **Grid_ptr,
				  QuadTreeBlock_DataStructure &QuadTree,
				  LESPremixed2D_Input_Parameters &Input_Parameters);

extern void Compute_Spectrum(LESPremixed2D_Quad_Block *Soln_ptr,
                             AdaptiveBlock2D_List &Soln_Block_List, 
                             LESPremixed2D_Input_Parameters &Input_Parameters,
                             double &u_average, double &v_average, double &vis);

extern void Rescale_Velocity(LESPremixed2D_Quad_Block *Soln_ptr,
			     AdaptiveBlock2D_List &Soln_Block_List,
			     LESPremixed2D_Input_Parameters &Input_Parameters,
			     const double &u_average, const double &v_average,
			     double &ko);

extern void Rescale_Spectrum(/* LESPremixed2D_Quad_Block *Soln_ptr, */
			     /* AdaptiveBlock2D_List &Soln_Block_List, */
			     LESPremixed2D_Input_Parameters &Input_Parameters,
			     double &u_average, double &v_average, double &vis);

extern void Write_Rescaled_Velocity(LESPremixed2D_Input_Parameters &Input_Parameters,
				    double *u, double *v, const int &Ny);

extern void Write_Rescaled_Velocity(LESPremixed2D_Quad_Block *Soln_ptr,
				    AdaptiveBlock2D_List &Soln_Block_List,
				    LESPremixed2D_Input_Parameters &Input_Parameters);

extern void Write_Initial_Turbulence(LESPremixed2D_Input_Parameters &Input_Parameters,
                                     double *u, double *v, const int &Ny);

extern void Write_Initial_Turbulence(Grid2D_Quad_Block   **Grid_ptr,
				     QuadTreeBlock_DataStructure &QuadTree,
				     LESPremixed2D_Input_Parameters &Input_Parameters,
				     double *u, double *v, 
				     const int &Ny);

extern void Read_Turbulent_Fluctuations(LESPremixed2D_Input_Parameters &Input_Parameters,
					double *u, double *v, 
					double &u_ave, double &v_ave,
					int &Ny, double &scale);

extern int Longitudinal_Correlation(QuadTreeBlock_DataStructure &QuadTree,
                                    AdaptiveBlockResourceList  &Global_Soln_Block_List,
				    AdaptiveBlock2D_List &Soln_Block_List,
				    LESPremixed2D_Quad_Block  *Soln_ptr,
				    LESPremixed2D_Input_Parameters  &Input_Parameters, 
				    const double & u_ave, const double &sqr_u);


/*************** END LESPREMIXED2D ***************************************/



#endif /* _LESPREMIXED2D_QUAD_INCLUDED  */
