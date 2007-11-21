/********************** Flame2DQuad.h **********************************
  Header file defining 2D Thermally Perfect Multiple Species Navier-Stokes
  Solution State Classes.

  - based on Euler2DQuad.h 
************************************************************************/

#ifndef _FLAME2D_QUAD_INCLUDED
#define _FLAME2D_QUAD_INCLUDED

/* Include 2D Chem state, 2D cell, and 2D quadrilateral block 
   grid, quadtree, and 2D Chem input header files. */
#include "../Grid/Cell2D.h"
#include "../Grid/Grid2DQuad.h"
#include "../AMR/QuadTree.h"
#include "../AMR/AMR.h"

// Also include linear systems header files.
#include "../Math/LinearSystems.h"

// Added and Header file for Chem
#include "Flame2DState.h"
#include "Flame2DInput.h"

// include SNBCK for optically thin radiation 
#include "../Physics/SNBCK/PlanckMean.h"

/* Define the structures and classes. */

#define	NUMBER_OF_RESIDUAL_VECTORS_FLAME2D    3  //K for dUdt[i][j][K]

/*!
 * Class: Flame2D_Quad_Block
 *
 * @brief Class definition of the 2D Chemistry solution blocks.
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
 *      S -- a 2D Chem solution                         *
 *                                                      *
 * S = S;                                               *
 * cout << S; (output function)                         *
 * cin  >> S; (input function)                          *
 *                                                      *
 * \endverbatim
 */
class Flame2D_Quad_Block{
  private:
  public:
  
  /***************** CLASS DATA *****************************************************/
  //@{ @name Solution state arrays:
  Flame2D_pState             **W;   //!< Primitive solution state.
  Flame2D_cState             **U;   //!< Conserved solution state. 
//   Mixture                   **M;   //!< Thermodynamic and transport mixture properties
  //@} 
  
  Flame2D_cState            **Uo;   //!< Initial solution state.
 
  double                   **dt;   //!< Local time step.  
  Flame2D_cState         ***dUdt;   //!< Solution residual.


  //@{ @name Boundary solution flux arrays:
  Flame2D_cState   *FluxN,*FluxS,   //!< Boundary solution fluxes.
                  *FluxE,*FluxW; 
  //@}  

  //@{ @name Boundary condtion reference states:
  Flame2D_pState       *WoN,*WoS,   //!< Boundary condition reference states for
                      *WoE,*WoW;   //!< north, south, east, & west boundaries.
  //@}  

  //@{ @name Solution gradient arrays:
  Flame2D_pState        **dWdx;     //!< (x-direction).
  Flame2D_pState        **dWdy;     //!< (y-direction).
  Flame2D_pState         **phi;     //!< Solution slope limiter.
  //@}  

  //Store face gradients for Diamond Path & Jacobian formation
  Flame2D_pState       **dWdx_faceN;   //!< north cell face(x-direction).
  Flame2D_pState       **dWdx_faceE;   //!< east  cell face(x-direction).
  Flame2D_pState       **dWdx_faceS;   //!< south cell face(x-direction).
  Flame2D_pState       **dWdx_faceW;   //!< west  cell face(x-direction).
  Flame2D_pState       **dWdy_faceN;   //!< north cell face(y-direction).
  Flame2D_pState       **dWdy_faceE;   //!< east  cell face(y-direction).
  Flame2D_pState       **dWdy_faceS;   //!< south cell face(y-direction).
  Flame2D_pState       **dWdy_faceW;   //!< west  cell face(y-direction).
  //@}  

  //@{ @name Grid block information:
  Grid2D_Quad_Block       Grid;   //!< 2D quadrilateral grid geometry.
  int              NCi,ICl,ICu;   //!< i-direction cell counters.
  int              NCj,JCl,JCu;   //!< j-direction cell counters.
  int                   Nghost;   //!< Number of ghost cells.
  //@}


  //@{ @name Radiation source term and static data:
  double                       **Srad; //!< radiant source term (divergence of rad. heat flux)
  static PlanckMean  *PlanckMean_data; //!< planck mean data object
  //@}


  //@{ @name FLAGS ( These are all esentially "static" Input parameters put in the SolnBlk for ease of access)
  int             Axisymmetric;   //!< Axisymmetric flow indicator.
  int                  Gravity;   //!< Gravity flag
  int                Flow_Type;   //!< Laminar, Viscous, Turbulent etc..
  int              debug_level;   //Debug level flag (0=none, 1,2,3,.. level of verboseness)
  int           Freeze_Limiter;   //Limiter freezing indicator (Multigrid/NKS) 0/1
  static int residual_variable;          //!< 1 = rho, 2,3 = rhou, 4 = E
  static int Number_of_Residual_Norms;   //!< (default 4 )
  double  Moving_wall_velocity;   //Moving Wall BC velocity
  double     Pressure_gradient;   //Moving Wall BC velocity
  //@}

  /****************** MEMBER FUNCTIONS ****************************************************/  
  /* Creation, copy, and assignment constructors. */
  Flame2D_Quad_Block(void) {
    NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
    W = NULL; U = NULL; dt = NULL; dUdt = NULL;
    dWdx = NULL; dWdy = NULL; dWdx_faceN = NULL; dWdy_faceN = NULL;
    dWdx_faceE = NULL; dWdy_faceE = NULL;
    dWdx_faceW = NULL; dWdy_faceW = NULL;
    dWdx_faceS = NULL; dWdy_faceS = NULL;
    phi = NULL;  Uo = NULL;
    FluxN = NULL; FluxS = NULL; FluxE = NULL; FluxW = NULL;
    WoN = NULL; WoS = NULL; WoE = NULL; WoW = NULL;
    Axisymmetric = 0; Gravity =0; Flow_Type = 0;
    Freeze_Limiter = OFF;
    Moving_wall_velocity=ZERO; Pressure_gradient =ZERO; debug_level=0;
    Srad = NULL;
  }
  
   
  /* Allocate memory for structured quadrilateral solution block. */
  void allocate(const int Ni, const int Nj, const int Ng);
  
  /* Deallocate memory for structured quadrilateral solution block. */
  void deallocate(void);
  static void deallocate_static(void);

  // resize the species arrays
  void resize_species(const int n_old, const string* species_old);

  /* Return primitive solution state at specified node. */
  Flame2D_pState Wn(const int &ii, const int &jj);
  
  /* Return conserverd solution state at specified node. */
  Flame2D_cState Un(const int &ii, const int &jj);

  /* Return conserverd solution state at specified node. */
  Flame2D_cState Uno(const int &ii, const int &jj);

  /* Return primitive solution state at cell nodes. */
  Flame2D_pState WnNW(const int &ii, const int &jj);
  Flame2D_pState WnNE(const int &ii, const int &jj);
  Flame2D_pState WnSE(const int &ii, const int &jj);
  Flame2D_pState WnSW(const int &ii, const int &jj);
  
  /* Return conserved solution state at cell nodes. */
  Flame2D_cState UnNW(const int &ii, const int &jj);
  Flame2D_cState UnNE(const int &ii, const int &jj);
  Flame2D_cState UnSE(const int &ii, const int &jj);
  Flame2D_cState UnSW(const int &ii, const int &jj);

  /* Return conserved solution state at cell nodes. */
  Flame2D_cState UnoNW(const int &ii, const int &jj);
  Flame2D_cState UnoNE(const int &ii, const int &jj);
  Flame2D_cState UnoSE(const int &ii, const int &jj);
  Flame2D_cState UnoSW(const int &ii, const int &jj);

  int BiLinearInterpolationCoefficients(double &eta, double &zeta, const int &ii, const int &jj);
  
  void set_v_zero(void);

  /*****************************************************************************
     dWn_dWc is the derivative of node solution w.r.t. cell center solution
     which is used by calculating the viscous Jacobians
  ******************************************************************************/
  double dWn_dWc (const int &i, const int &j, const int &Orient);
  
  /* Set flags for limiter evaluation. */
  void evaluate_limiters(void) {Freeze_Limiter = OFF; } 
  void freeze_limiters(void) {Freeze_Limiter = ON; }
  
  /* Input-output operators. */
  friend ostream &operator << (ostream &out_file,
			       const Flame2D_Quad_Block
			       &Soln);
  friend istream &operator >> (istream &in_file,
			       Flame2D_Quad_Block
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
};

/**************************************************************************
 * Flame2D_Quad_Block::allocate -- Allocate memory.                        *
 **************************************************************************/
inline void Flame2D_Quad_Block::allocate(const int Ni, const int Nj, const int Ng) {
  assert(Ni > 1 && Nj > 1); Grid.allocate(Ni, Nj, Ng);
   NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1;
   NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1; Nghost = Ng;
   W = new Flame2D_pState*[NCi]; U = new Flame2D_cState*[NCi];
   dt = new double*[NCi];
   dUdt = new Flame2D_cState**[NCi];
   dWdx = new Flame2D_pState*[NCi]; dWdy = new Flame2D_pState*[NCi];
   dWdx_faceN = new Flame2D_pState*[NCi]; dWdy_faceN = new Flame2D_pState*[NCi];
   dWdx_faceE = new Flame2D_pState*[NCi]; dWdy_faceE = new Flame2D_pState*[NCi];
   dWdx_faceW = new Flame2D_pState*[NCi]; dWdy_faceW = new Flame2D_pState*[NCi];
   dWdx_faceS = new Flame2D_pState*[NCi]; dWdy_faceS = new Flame2D_pState*[NCi];
   phi = new Flame2D_pState*[NCi]; Uo = new Flame2D_cState*[NCi];
   Srad = new double*[NCi];
   
   for (int i = 0; i <= NCi-1 ; ++i ) {
      W[i] = new Flame2D_pState[NCj]; U[i] = new Flame2D_cState[NCj];
      dt[i] = new double[NCj]; dUdt[i] = new Flame2D_cState*[NCj];
      for (int j = 0; j <= NCj-1 ; ++j ){
	dUdt[i][j] = new Flame2D_cState[NUMBER_OF_RESIDUAL_VECTORS_FLAME2D];
      }
      dWdx[i] = new Flame2D_pState[NCj]; dWdy[i] = new Flame2D_pState[NCj];
      dWdx_faceN[i] = new Flame2D_pState[NCj]; dWdy_faceN[i] = new Flame2D_pState[NCj];
      dWdx_faceE[i] = new Flame2D_pState[NCj]; dWdy_faceE[i] = new Flame2D_pState[NCj];
      dWdx_faceW[i] = new Flame2D_pState[NCj]; dWdy_faceW[i] = new Flame2D_pState[NCj];
      dWdx_faceS[i] = new Flame2D_pState[NCj]; dWdy_faceS[i] = new Flame2D_pState[NCj];
      phi[i] = new Flame2D_pState[NCj]; Uo[i] = new Flame2D_cState[NCj];
      Srad[i] = new double[NCj];
   } /* endfor */
   FluxN = new Flame2D_cState[NCi]; FluxS = new Flame2D_cState[NCi];
   FluxE = new Flame2D_cState[NCj]; FluxW = new Flame2D_cState[NCj];
   WoN = new Flame2D_pState[NCi]; WoS = new Flame2D_pState[NCi];
   WoE = new Flame2D_pState[NCj]; WoW = new Flame2D_pState[NCj];

  // Set the solution residuals, gradients, limiters, and other values to zero.
   for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
      for ( int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
          for ( int k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_FLAME2D-1 ; ++k ) {
	     dUdt[i][j][k].Vacuum();
          } 
	  dWdx[i][j].Vacuum() ; dWdy[i][j].Vacuum();
	  phi[i][j].Vacuum(); Uo[i][j].Vacuum();
	  dt[i][j] = ZERO;
	  Srad[i][j] = ZERO;
      } 
   }  
}

/**************************************************************************
 * Flame2D_Quad_Block::deallocate -- Deallocate memory.                   *
 **************************************************************************/
inline void Flame2D_Quad_Block::deallocate(void) {
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
      delete []Srad[i]; Srad[i] = NULL; 
   } /* endfor */
   delete []W; W = NULL; delete []U; U = NULL;
   delete []dt; dt = NULL; delete []dUdt; dUdt = NULL;
   delete []dWdx; dWdx = NULL; delete []dWdy; dWdy = NULL;
   delete []dWdx_faceN; dWdx_faceN = NULL; delete []dWdy_faceN; dWdy_faceN = NULL;
   delete []dWdx_faceE; dWdx_faceE = NULL; delete []dWdy_faceE; dWdy_faceE = NULL;
   delete []dWdx_faceW; dWdx_faceW = NULL; delete []dWdy_faceW; dWdy_faceW = NULL;
   delete []dWdx_faceS; dWdx_faceS = NULL; delete []dWdy_faceS; dWdy_faceS = NULL;

   delete []phi; phi = NULL; delete []Uo; Uo = NULL;  
   delete []Srad; Srad = NULL;
   delete []FluxN; FluxN = NULL; delete []FluxS; FluxS = NULL;
   delete []FluxE; FluxE = NULL; delete []FluxW; FluxW = NULL;
   delete []WoN; WoN = NULL; delete []WoS; WoS = NULL;
   delete []WoE; WoE = NULL; delete []WoW; WoW = NULL;
   NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
}

/**************************************************************************
 * Flame2D_Quad_Block::deallocate_static -- Deallocate static memory.      *
 **************************************************************************/
inline void Flame2D_Quad_Block::deallocate_static(void) {
  if(PlanckMean_data != NULL) { delete PlanckMean_data; PlanckMean_data = NULL; }
}


/**************************************************************************
 * Flame2D_Quad_Block::resize_species                                      *
 *                                                                        *
 * Resize all the species arrays. This function is DISCUSTING!!!!.        *
 * Because the number of species is static, we have to carefully resize   *
 * all the species arrays in the solution block.  I hope that there are no*
 * Flame2D_pState and Flame2D_cState objects allocated anywhere else.       *
 * If you resize the reaction mechanism and this function is called,      *
 * and you get a SIGSEGV, you should look for leftover state objects.     *
 **************************************************************************/
inline void Flame2D_Quad_Block::resize_species(const int n_old, 
					      const string* species_old) {
  
  // temporary storage
  double* c_old = new double[n_old];
  int index;


  //
  // Loop over the grid
  //
  for ( int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
    for ( int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {

      //---------------------------------------------------
      // CONSERVED STATE
      U[i][j].resize_species(n_old, species_old); // keep old values
      //---------------------------------------------------
      // PRIMITIVE STATE
      W[i][j].resize_species(); // throw out old values
      W[i][j] = U[i][j].W();
      //---------------------------------------------------
      // REMAINING STATES
      for ( int k = 0 ; k < NUMBER_OF_RESIDUAL_VECTORS_FLAME2D; ++k ) {
	dUdt[i][j][k].resize_species(); // throw out old values
      } // endfor
      dWdx[i][j].resize_species(); // throw out old values
      dWdy[i][j].resize_species(); // throw out old values
      phi[i][j].resize_species(); // throw out old values
      Uo[i][j].resize_species(); // throw out old values
      //---------------------------------------------------

    } // endfor
  } // endfor
  

  //
  // resize boundary ref states
  //
  for ( int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {

    //---------------------------------------------------
    // WEST STATE
    WoW[j].resize_species(n_old, species_old); // keep old values
    //---------------------------------------------------
    // EAST STATE
    WoE[j].resize_species(n_old, species_old); // keep old values
    //---------------------------------------------------

  } // endfor


  //
  // resize boundary ref states
  //
  for ( int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {

    //---------------------------------------------------
    // SOUTH STATE
    WoS[i].resize_species(n_old, species_old); // keep old values
    //---------------------------------------------------
    // NORTH STATE
    WoN[i].resize_species(n_old, species_old); // keep old values
    //---------------------------------------------------

  } // endfor

  // delete temporary memory
  delete[] c_old;
}

/**************************************************************************
 * Flame2D_Quad_Block:BiLinearInterpolationCoefficients --                 *
 *                     Coefficients used by bilinear interpolation        *
 **************************************************************************/
inline int Flame2D_Quad_Block::BiLinearInterpolationCoefficients(double &eta, double &zeta, const int &ii, const int &jj){

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


/**************************************************************************
 * Flame2D_Quad_Block::dWn_dWc -- Derivative of Node primitive solution    *  
 *                                 w.r.t Cell primitive solution          *
 **************************************************************************/
inline double Flame2D_Quad_Block::dWn_dWc(const int &i, const int &j, const int &Orient) {
 
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
    cerr<<"\n Improper Orient in Flame2D_Quad_Block::dWn_dWc\n";
    break;
  }
      
}


/**************************************************************************
 * Flame2D_Quad_Block::Wn -- Node primitive solution.                     *
 **************************************************************************/
inline Flame2D_pState Flame2D_Quad_Block::Wn(const int &ii, const int &jj) {
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
   
  return (W[ii-1][jj-1] +(W[ii-1][jj]-W[ii-1][jj-1])*zeta+  
	  (W[ii][jj-1]-W[ii-1][jj-1])*eta + 
	  (W[ii][jj]+W[ii-1][jj-1]-W[ii-1][jj]-W[ii][jj-1])*zeta*eta);  

}

/**************************************************************************
/**************************************************************************
 * Flame2D_Quad_Block::Wn?? -- Get cell node primitive solution states.   *
 **************************************************************************/
inline Flame2D_pState Flame2D_Quad_Block::WnNW(const int &ii, const int &jj) {
  return (Wn(ii, jj+1));
}

inline Flame2D_pState Flame2D_Quad_Block::WnNE(const int &ii, const int &jj) {
  return (Wn(ii+1, jj+1));
}

inline Flame2D_pState Flame2D_Quad_Block::WnSE(const int &ii, const int &jj) {
  return (Wn(ii+1, jj));
}

inline Flame2D_pState Flame2D_Quad_Block::WnSW(const int &ii, const int &jj) {
  return (Wn(ii, jj));
}

/**************************************************************************
 * Flame2D_Quad_Block::Un -- Node conservative solution.                     *
 **************************************************************************/
inline Flame2D_cState Flame2D_Quad_Block::Un(const int &ii, const int &jj) {
 
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
  
}

/**************************************************************************
 * Flame2D_Quad_Block::Un -- Get cell node conserved solution states.    *
 **************************************************************************/
inline Flame2D_cState Flame2D_Quad_Block::UnNW(const int &ii, const int &jj) {
  return (Un(ii, jj+1));
}

inline Flame2D_cState Flame2D_Quad_Block::UnNE(const int &ii, const int &jj) {
  return (Un(ii+1, jj+1));
}

inline Flame2D_cState Flame2D_Quad_Block::UnSE(const int &ii, const int &jj) {
  return (Un(ii+1, jj));
}

inline Flame2D_cState Flame2D_Quad_Block::UnSW(const int &ii, const int &jj) {
  return (Un(ii, jj));
}

/**************************************************************************
 * Flame2D_Quad_Block::Un -- Node conservative solution.                     *
 **************************************************************************/
inline Flame2D_cState Flame2D_Quad_Block::Uno(const int &ii, const int &jj) {
 
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

/**************************************************************************
 * Flame2D_Quad_Block::Un -- Get cell node conserved solution states.    *
 **************************************************************************/
inline Flame2D_cState Flame2D_Quad_Block::UnoNW(const int &ii, const int &jj) {
  return (Uno(ii, jj+1));
}

inline Flame2D_cState Flame2D_Quad_Block::UnoNE(const int &ii, const int &jj) {
  return (Uno(ii+1, jj+1));
}

inline Flame2D_cState Flame2D_Quad_Block::UnoSE(const int &ii, const int &jj) {
  return (Uno(ii+1, jj));
}

inline Flame2D_cState Flame2D_Quad_Block::UnoSW(const int &ii, const int &jj) {
  return (Uno(ii, jj));
}


/**************************************************************************
 * Flame2D_Quad_Block -- Input-output operators.                           *
 **************************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Flame2D_Quad_Block &SolnBlk) {
  int i, j; 
  out_file << SolnBlk.Grid;
  out_file << SolnBlk.NCi << " " << SolnBlk.ICl << " " << SolnBlk.ICu << " " << SolnBlk.Nghost << "\n";
  out_file << SolnBlk.NCj << " " << SolnBlk.JCl << " " << SolnBlk.JCu << "\n";
  out_file << SolnBlk.Axisymmetric << "\n";
  out_file << SolnBlk.Flow_Type <<"\n";
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
			     Flame2D_Quad_Block &SolnBlk) {

  int i, j, k, ni, il, iu, nj, jl, ju, ng;
  Flame2D_pState Flame2D_W_VACUUM; Flame2D_W_VACUUM.Vacuum();
  Flame2D_cState Flame2D_U_VACUUM; Flame2D_U_VACUUM.Vacuum();
  Grid2D_Quad_Block New_Grid; in_file >> New_Grid;
  in_file.setf(ios::skipws);
  in_file >> ni >> il >> iu >> ng; in_file >> nj >> jl >> ju;
  in_file >> SolnBlk.Axisymmetric;
  in_file >> SolnBlk.Flow_Type;
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
         for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_FLAME2D-1 ; ++k ) {
	     SolnBlk.dUdt[i][j][k] = Flame2D_U_VACUUM;
         } /* endfor */
	 SolnBlk.dWdx[i][j] = Flame2D_W_VACUUM;
	 SolnBlk.dWdy[i][j] = Flame2D_W_VACUUM;
	 SolnBlk.dWdx_faceN[i][j] = Flame2D_W_VACUUM;
	 SolnBlk.dWdy_faceN[i][j] = Flame2D_W_VACUUM;
	 SolnBlk.dWdx_faceE[i][j] = Flame2D_W_VACUUM;
	 SolnBlk.dWdy_faceE[i][j] = Flame2D_W_VACUUM;
	 SolnBlk.dWdx_faceW[i][j] = Flame2D_W_VACUUM;
	 SolnBlk.dWdy_faceW[i][j] = Flame2D_W_VACUUM;
	 SolnBlk.dWdx_faceS[i][j] = Flame2D_W_VACUUM;
	 SolnBlk.dWdy_faceS[i][j] = Flame2D_W_VACUUM;
	 SolnBlk.phi[i][j] = Flame2D_W_VACUUM;
	 SolnBlk.Uo[i][j] = Flame2D_U_VACUUM;
	 SolnBlk.dt[i][j] = ZERO;
	 SolnBlk.Srad[i][j] = ZERO;
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
 * Flame2D_Quad_Block::NumVar -- Returns number of state variables.            *
 *******************************************************************************/
inline int Flame2D_Quad_Block::NumVar(void) {
  return (int(W[0][0].NUM_VAR_FLAME2D));    
}

/*******************************************************************************
 * Flame2D_Quad_Block::LoadSendBuffer -- Loads send message buffer.            *
 *******************************************************************************/
inline int Flame2D_Quad_Block::LoadSendBuffer(double *buffer,
                                              int &buffer_count,
                                              const int buffer_size,
                                              const int i_min, 
                                              const int i_max,
                                              const int i_inc,
                                              const int j_min, 
                                              const int j_max,
                                              const int j_inc) {
  int NUM_VAR_FLAME2D = NumVar();
  int i, j, k;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = U[i][j][k];
        } /* endfor */
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * Flame2D_Quad_Block::LoadSendBuffer_F2C -- Loads send message buffer for     *
 *                                           fine to coarse block message      *
 *                                           passing.                          *
 *******************************************************************************/
inline int Flame2D_Quad_Block::LoadSendBuffer_F2C(double *buffer,
                                                  int &buffer_count,
                                                  const int buffer_size,
                                                  const int i_min, 
                                                  const int i_max,
                                                  const int i_inc,
                                                  const int j_min, 
                                                  const int j_max,
                                                  const int j_inc) {
  int NUM_VAR_FLAME2D = NumVar();
  int i, j, k;
  for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
 * Flame2D_Quad_Block::LoadSendBuffer_C2F -- Loads send message buffer for     *
 *                                           coarse to fine block message      *
 *                                           passing.                          *
 *******************************************************************************/
inline int Flame2D_Quad_Block::LoadSendBuffer_C2F(double *buffer,
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
  int NUM_VAR_FLAME2D = NumVar();
  int i, j, k;
  Vector2D dX;
  Flame2D_pState Wfine;
  Flame2D_cState Ufine;

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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) { 
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
              for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
 * Flame2D_Quad_Block::UnloadReceiveBuffer -- Unloads receive message buffer.  *
 *******************************************************************************/
inline int Flame2D_Quad_Block::UnloadReceiveBuffer(double *buffer,
                                                   int &buffer_count,
                                                   const int buffer_size,
                                                   const int i_min, 
                                                   const int i_max,
                                                   const int i_inc,
                                                   const int j_min, 
                                                   const int j_max,
                                                   const int j_inc) { 
  int NUM_VAR_FLAME2D = NumVar();
 
  for ( int j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
    for ( int i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      //Changed from Euler2D
      for ( int k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
 * Flame2D_Quad_Block::UnloadReceiveBuffer_F2C -- Unloads receive message      *
 *                                                buffer for fine to coarse    *
 *                                                block message passing.       *
 *******************************************************************************/
inline int Flame2D_Quad_Block::UnloadReceiveBuffer_F2C(double *buffer,
                                                       int &buffer_count,
                                                       const int buffer_size,
                                                       const int i_min, 
                                                       const int i_max,
                                                       const int i_inc,
                                                       const int j_min, 
                                                       const int j_max,
                                                       const int j_inc) { 
  int NUM_VAR_FLAME2D = NumVar();
  int i, j;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
       //Changed from Euler2D
       for ( int k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
 * Flame2D_Quad_Block::UnloadReceiveBuffer_C2F -- Unloads receive message      *
 *                                                buffer for coarse to fine    *
 *                                                block message passing.       *
 *******************************************************************************/
inline int Flame2D_Quad_Block::UnloadReceiveBuffer_C2F(double *buffer,
                                                       int &buffer_count,
                                                       const int buffer_size,
                                                       const int i_min, 
                                                       const int i_max,
                                                       const int i_inc,
                                                       const int j_min, 
                                                       const int j_max,
                                                       const int j_inc) { 
  int NUM_VAR_FLAME2D = NumVar();
  int i, j;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
       //Changed from Euler2D
       for ( int k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
 * Flame2D_Quad_Block::SubcellReconstruction --                           *
 *               Performs the subcell reconstruction of solution state    *
 *               within a given cell (i,j) of the computational mesh for  *
 *               the specified quadrilateral solution block.              *
 **************************************************************************/
inline void Flame2D_Quad_Block::SubcellReconstruction(const int i, 
                                                      const int j,
                                                      const int Limiter) {

  int n, n2, n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4], phi_n;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  Flame2D_pState DU, DUDx_ave, DUDy_ave;

  Flame2D_pState Flame2D_W_VACUUM;
  Flame2D_W_VACUUM.Vacuum();
 
  int NUM_VAR_FLAME2D = NumVar(); 
 
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
      DUDx_ave = Flame2D_W_VACUUM;
      DUDy_ave = Flame2D_W_VACUUM;
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
	for ( n = 1 ; n <= NUM_VAR_FLAME2D ; ++n ) {
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
      dWdx[i][j] = Flame2D_W_VACUUM;
      dWdy[i][j] = Flame2D_W_VACUUM; 
      phi[i][j]  = Flame2D_W_VACUUM;
  } /* endif */
    
}


/*******************************************************************************
 * Flame2D_Quad_Block::LoadSendBuffer_Flux_F2C -- Loads send message buffer for*
 *                                                fine to coarse block message *
 *                                                passing of conservative      *
 *                                                solution fluxes.             *
 *******************************************************************************/
inline int Flame2D_Quad_Block::LoadSendBuffer_Flux_F2C(double *buffer,
                                                       int &buffer_count,
                                                       const int buffer_size,
                                                       const int i_min, 
                                                       const int i_max,
                                                       const int i_inc,
                                                       const int j_min, 
                                                       const int j_max,
                                                       const int j_inc) { 
  int NUM_VAR_FLAME2D = NumVar(); 
  int i, j, k;
  if (j_min == j_max && j_min == JCl) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxS[i  ][k]+
                                   FluxS[i+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxN[i  ][k]+
                                   FluxN[i+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
     for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
        for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxW[j][k]+
                                   FluxW[j+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
     for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
        for ( k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
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
 * Flame2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C -- Unloads receive message *
 *                                                buffer for fine to coarse    *
 *                                                block message passing of     *
 *                                                conservative solution fluxes.*
 *******************************************************************************/
inline int Flame2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C(double *buffer,
                                                            int &buffer_count,
                                                            const int buffer_size,
                                                            const int i_min, 
                                                            const int i_max,
                                                            const int i_inc,
                                                            const int j_min, 
                                                            const int j_max,
                                                            const int j_inc) {
  int NUM_VAR_FLAME2D = NumVar();
  if (j_min == j_max && j_min == JCl) {
    for ( int i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      for ( int k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	FluxS[i][k] = - buffer[buffer_count] - FluxS[i][k];
      }
    }
  } else if (j_min == j_max && j_min == JCu) {
    for ( int i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
      for ( int k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	FluxN[i][k] = - buffer[buffer_count] - FluxN[i][k];
      }
    } 
  } else if (i_min == i_max && i_min == ICl) {
    for ( int j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
      for ( int k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	FluxW[j][k] = - buffer[buffer_count] - FluxW[j][k];
      }
    } 
  } else if (i_min == i_max && i_min == ICu) {
    for ( int j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
      for ( int k = 1 ; k <= NUM_VAR_FLAME2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	FluxE[j][k] = - buffer[buffer_count] - FluxE[j][k];
      }
    } 
  } /* endif */
  return(0);
}

/**************************************************************************
 * Flame2D_Quad_Block -- Single Block External Subroutines.                * 
 * Flame2DQuadSingleBlock.cc                                               *
 **************************************************************************/

extern void Write_Solution_Block(Flame2D_Quad_Block &SolnBlk,
	                         ostream &Out_File);

extern void Read_Solution_Block(Flame2D_Quad_Block &SolnBlk,
	                        istream &In_File);

extern void Broadcast_Solution_Block(Flame2D_Quad_Block &SolnBlk);

#ifdef _MPI_VERSION
extern void Broadcast_Solution_Block(Flame2D_Quad_Block &SolnBlk,
                                     MPI::Intracomm &Communicator, 
                                     const int Source_CPU);
#endif

extern void Copy_Solution_Block(Flame2D_Quad_Block &SolnBlk1,
		                Flame2D_Quad_Block &SolnBlk2);

extern int Prolong_Solution_Block(Flame2D_Quad_Block &SolnBlk_Fine,
		                  Flame2D_Quad_Block &SolnBlk_Original,
                                  const int Sector);

extern int Restrict_Solution_Block(Flame2D_Quad_Block &SolnBlk_Coarse,
		                   Flame2D_Quad_Block &SolnBlk_Original_SW,
                                   Flame2D_Quad_Block &SolnBlk_Original_SE,
                                   Flame2D_Quad_Block &SolnBlk_Original_NW,
                                   Flame2D_Quad_Block &SolnBlk_Original_NE);

extern void Output_Tecplot(Flame2D_Quad_Block &SolnBlk,
			   Flame2D_Input_Parameters &IP,
			   const int Number_of_Time_Steps,
                           const double &Time,
                           const int Block_Number,
                           const int Output_Title,
	                   ostream &Out_File);

extern void Output_Cells_Tecplot(Flame2D_Quad_Block &SolnBlk,
				 Flame2D_Input_Parameters &IP,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_Nodes_Tecplot(Flame2D_Quad_Block &SolnBlk,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_RHS(Flame2D_Quad_Block &SolnBlk,
		       const int Number_of_Time_Steps,
                       const double &Time,
                       const int Block_Number,
                       const int Output_Title,
	               ostream &Out_File);

extern void ICs(Flame2D_Quad_Block &SolnBlk,
 	        const int i_ICtype,
                Flame2D_pState *Wo, Flame2D_Input_Parameters &Input_Parameters);

extern void Reset_Wo(Flame2D_Quad_Block &SolnBlk );

extern void BCs(Flame2D_Quad_Block &SolnBlk, 
		Flame2D_Input_Parameters &IP);

extern double CFL(Flame2D_Quad_Block &SolnBlk,
                  Flame2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(Flame2D_Quad_Block &SolnBlk, 
                                const double &Dt_min);

extern double L1_Norm_Residual(Flame2D_Quad_Block &SolnBlk, const int &norm);

extern double L2_Norm_Residual(Flame2D_Quad_Block &SolnBlk, const int &norm);

extern double Max_Norm_Residual(Flame2D_Quad_Block &SolnBlk, const int &norm);


/**************** Reconstruction Methods ****************************/

extern void Linear_Reconstruction_LeastSquares(Flame2D_Quad_Block &SolnBlk,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares(Flame2D_Quad_Block &SolnBlk,
                                               const int i,
                                               const int j,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares_2(Flame2D_Quad_Block &SolnBlk,
                                                 const int i,
                                                 const int j,
					         const int Limiter);

extern void Linear_Reconstruction_GreenGauss(Flame2D_Quad_Block &SolnBlk,
					     const int Limiter);

extern void Linear_Reconstruction_GreenGauss(Flame2D_Quad_Block &SolnBlk,
                                             const int i,
                                             const int j,
					     const int Limiter);

/******* Diamond Path Reconstruction Methods *************************/
extern void Linear_Reconstruction_LeastSquares_Diamond(Flame2D_Quad_Block &SolnBlk,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares_Diamond(Flame2D_Quad_Block &SolnBlk,
						       const int i, 
						       const int j,
						       const int Limiter);


extern void Linear_Reconstruction_GreenGauss_Diamond(Flame2D_Quad_Block &SolnBlk,
						     const int Limiter);

extern void Linear_Reconstruction_GreenGauss_Diamond(Flame2D_Quad_Block &SolnBlk,
						     const int i,
						     const int j,
						     const int Limiter);

extern void Residual_Smoothing(Flame2D_Quad_Block &SolnBlk,
                               const int k_residual,
			       double &epsilon, 
                               const int number_of_Gauss_Seidel_iterations);

extern void Calculate_Refinement_Criteria(double *refinement_criteria,
					  Flame2D_Input_Parameters &IP, 
                                          int &number_refinement_criteria,					  
                                          Flame2D_Quad_Block &SolnBlk);

extern void Fix_Refined_Block_Boundaries(Flame2D_Quad_Block &SolnBlk,
                                         const int Fix_North_Boundary,
                                         const int Fix_South_Boundary,
                                         const int Fix_East_Boundary,
                                         const int Fix_West_Boundary);

extern void Unfix_Refined_Block_Boundaries(Flame2D_Quad_Block &SolnBlk);

extern void Apply_Boundary_Flux_Corrections(Flame2D_Quad_Block &SolnBlk,
                                            const int Number_Neighbours_North_Boundary,
                                            const int Number_Neighbours_South_Boundary,
                                            const int Number_Neighbours_East_Boundary,
                                            const int Number_Neighbours_West_Boundary);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Flame2D_Quad_Block &SolnBlk,
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

extern int dUdt_Residual_Evaluation(Flame2D_Quad_Block &SolnBlk,
				    Flame2D_Input_Parameters &Input_Parameters);

extern int dUdt_Multistage_Explicit(Flame2D_Quad_Block &SolnBlk,
   	                            const int i_stage,
                                    Flame2D_Input_Parameters &Input_Parameters);

extern int Update_Solution_Multistage_Explicit(Flame2D_Quad_Block &SolnBlk,
   	                                       const int i_stage,
                                               Flame2D_Input_Parameters &Input_Parameters);

extern void Viscous_Calculations(Flame2D_Quad_Block &SolnBlk);

// optically thin radiation source term evaluation
extern void Radiation_Source_Eval( Flame2D_Quad_Block &SolnBlk,
				   Flame2D_Input_Parameters &Input_Parameters );

/**************************************************************************
 * Flame2D_Quad_Block -- Multiple Block External Subroutines.              *
 * Flame2DQuadMultiBlock.cc                                                *
 **************************************************************************/

extern Flame2D_Quad_Block* Allocate(Flame2D_Quad_Block *Soln_ptr,
                                    Flame2D_Input_Parameters &Input_Parameters);

extern Flame2D_Quad_Block* Deallocate(Flame2D_Quad_Block *Soln_ptr,
                                      Flame2D_Input_Parameters &Input_Parameters);

extern void ICs(Flame2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
                Flame2D_Input_Parameters &Input_Parameters);
  
extern int Read_Restart_Solution(Flame2D_Quad_Block *Soln_ptr,
                                 AdaptiveBlock2D_List &Soln_Block_List,
                                 Flame2D_Input_Parameters &Input_Parameters,
		                 int &Number_of_Time_Steps,
                                 double &Time,
                                 CPUTime &CPU_Time);

extern int Write_Restart_Solution(Flame2D_Quad_Block *Soln_ptr,
                                  AdaptiveBlock2D_List &Soln_Block_List,
                                  Flame2D_Input_Parameters &Input_Parameters,
		                  const int Number_of_Time_Steps,
                                  const double &Time,
                                  const CPUTime &CPU_Time);

extern int Output_Tecplot(Flame2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          Flame2D_Input_Parameters &Input_Parameters,
		          const int Number_of_Time_Steps,
                          const double &Time);

extern int Output_RHS(Flame2D_Quad_Block *Soln_ptr,
		      AdaptiveBlock2D_List &Soln_Block_List,
		      Flame2D_Input_Parameters &Input_Parameters,
		      const int Number_of_Time_Steps,
		      const double &Time);

int Output_PERTURB(Flame2D_Quad_Block *Soln_ptr,
		   AdaptiveBlock2D_List &Soln_Block_List,
		   Flame2D_Input_Parameters &Input_Parameters,
		   const int Number_of_Time_Steps,
		   const double &Time,
		   const CPUTime &CPU_Time);

extern int Output_Cells_Tecplot(Flame2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                Flame2D_Input_Parameters &Input_Parameters,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Mesh_Tecplot(Flame2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Flame2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Mesh_Gnuplot(Flame2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Flame2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Ringleb(Flame2D_Quad_Block *Soln_ptr,
			  AdaptiveBlock2D_List &Soln_Block_List,
			  Flame2D_Input_Parameters &IP);

extern int Output_Viscous_Channel(Flame2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  Flame2D_Input_Parameters &IP);

extern int Output_Flat_Plate(Flame2D_Quad_Block *Soln_ptr,
			     AdaptiveBlock2D_List &Soln_Block_List,
			     Flame2D_Input_Parameters &IP);

extern int Output_Driven_Cavity_Flow(Flame2D_Quad_Block *Soln_ptr,
				     AdaptiveBlock2D_List &Soln_Block_List,
				     Flame2D_Input_Parameters &IP);

extern void BCs(Flame2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
		Flame2D_Input_Parameters &Input_Parameters);

extern double CFL(Flame2D_Quad_Block *Soln_ptr,
                  AdaptiveBlock2D_List &Soln_Block_List,
                  Flame2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(Flame2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List, 
                                const double &Dt_min);

extern double L1_Norm_Residual(Flame2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);
		
extern double L2_Norm_Residual(Flame2D_Quad_Block *Soln_ptr ,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double Max_Norm_Residual(Flame2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List);

extern void L1_Norm_Residual(Flame2D_Quad_Block *Soln_ptr, 
                               AdaptiveBlock2D_List &Soln_Block_List,
			       double *l1_norm);

extern void L2_Norm_Residual(Flame2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
			       double *l2_norm);

extern void Max_Norm_Residual(Flame2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
				double *max_norm);

extern void Evaluate_Limiters(Flame2D_Quad_Block *Soln_ptr,
                              AdaptiveBlock2D_List &Soln_Block_List);

extern void Freeze_Limiters(Flame2D_Quad_Block *Soln_ptr,
                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Change_Mref(Flame2D_Quad_Block *Soln_ptr,
			AdaptiveBlock2D_List &Soln_Block_List,
			const double &Mr);

extern void Residual_Smoothing(Flame2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Flame2D_Input_Parameters &Input_Parameters,
   	                       const int I_Stage);

extern void Apply_Boundary_Flux_Corrections(Flame2D_Quad_Block *Soln_ptr,
                                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Flame2D_Quad_Block *Soln_ptr,
                                                                AdaptiveBlock2D_List &Soln_Block_List,
                                                                Flame2D_Input_Parameters &Input_Parameters,
   	                                                        const int I_Stage);

extern int dUdt_Multistage_Explicit(Flame2D_Quad_Block *Soln_ptr,
				    AdaptiveBlockResourceList &Global_Soln_Block_List,
                                    AdaptiveBlock2D_List &Soln_Block_List,
                                    Flame2D_Input_Parameters &Input_Parameters,
   	                            const int I_Stage);

extern int Update_Solution_Multistage_Explicit(Flame2D_Quad_Block *Soln_ptr,
                                               AdaptiveBlock2D_List &Soln_Block_List,
                                               Flame2D_Input_Parameters &Input_Parameters,
   	                                       const int I_Stage);



//************ Viscous Stuff************************************/
extern void Viscous_Calculations(Flame2D_Quad_Block &SolnBlk);

// extern Flame2D_cState Flux_Viscous_x(Flame2D_Quad_Block &SolnBlk, const int &i, const int &j);
// extern Flame2D_cState Flux_Viscous_y(Flame2D_Quad_Block &SolnBlk, const int &i, const int &j);
// extern Flame2D_cState Flux_Viscous_n(Flame2D_Quad_Block &SolnBlk, const int &Orient, const int &i, const int &j);
// extern Flame2D_cState Flux_Viscous_EastFace_Diamond(Flame2D_Quad_Block &SolnBlk, const int &i, const int &j);
// extern Flame2D_cState Flux_Viscous_NorthFace_Diamond(Flame2D_Quad_Block &SolnBlk, const int &i, const int &j);


/**************************************************************************
 * Flame2D_Quad_Block -- Multiple Block External Subroutines for Mesh.     * 
 * Flame2DQuadGrid.cc                                                      *
 **************************************************************************/

extern Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                            Flame2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Broadcast_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                      Flame2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                             Flame2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                                            Flame2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                  Flame2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                 Flame2D_Input_Parameters &Input_Parameters);

extern int Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                          Flame2D_Input_Parameters &Input_Parameters);

extern int Output_Nodes_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                Flame2D_Input_Parameters &Input_Parameters);

extern int Output_Cells_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                Flame2D_Input_Parameters &Input_Parameters);


/**************************************************************************
 * Flame2D_Quad_Block -- Solvers.                                          *
 * Flame2DQuadSolvers.cc                                                   *
 **************************************************************************/

extern int Flame2DQuadSolver(char *Input_File_Name_ptr, int batch_flag);

/**************************************************************************
 * Flame2D_Tools -- Tools and stuff                                         *
 * Flame2DTools.cc                                                          *
 **************************************************************************/

extern int Open_Time_Accurate_File(ofstream &Time_Accurate_File,
				   char *File_Name,
				   const int Append_to_Fileconst,
				   const Flame2D_pState &Soln);

extern int Close_Time_Accurate_File(ofstream &Time_Accurate_File);

extern void Output_to_Time_Accurate_File(ostream &Time_Accurate_File,
					 const double &Time,
					 const Flame2D_pState &Soln);
  
extern void Output_Ringleb_Solution(Flame2D_Quad_Block &SolnBlk,
				    const int Block_Number,
				    const int Output_Title,
				    ostream &Out_File);

extern void Output_Ringleb_Error(Flame2D_Quad_Block &SolnBlk,
				 double &l1_norm,
				 double &l2_norm,
				 double &max_norm,
				 int &numberofactivecells);

extern void Output_Viscous_Channel(Flame2D_Quad_Block &SolnBlk,
				   const int Block_Number,
				   const int Output_Title,
				   ostream &Out_File,
				   double &l1_norm,
				   double &l2_norm,
				   double &max_norm,
				   double &Vwall,
				   const double dp);

extern void Output_Flat_Plate(Flame2D_Quad_Block &SolnBlk,
			      const int Block_Number,
			      const int Output_Title_Soln,
			      ostream &Out_File_Soln,
			      const int Output_Title_Skin,
			      ostream &Out_File_Skin,
			      const Flame2D_pState &Winf,
			      double &l1_norm,
			      double &l2_norm,
			      double &max_norm);
			      
extern void Output_Driven_Cavity_Flow(Flame2D_Quad_Block &SolnBlk,
				      const int Block_Number,
				      const int Output_Title,
				      ostream &Out_File_u,
				      ostream &Out_File_v,
				      const double &Re,
				      const double &Vwall,
				      const double &length);

extern int Output_Quasi3D_Tecplot(Flame2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  Flame2D_Input_Parameters &IP,
				  const int Number_of_Time_Steps,
				  const double &Time);

extern void Output_Quasi3D_Tecplot(Flame2D_Quad_Block &SolnBlk,
				   Flame2D_Input_Parameters &IP,
				   const int Number_of_Time_Steps,
				   const double &Time,
				   const int Block_Number,
				   const int Output_Title,
				   ostream &Out_File);

int Set_Equilibrium_State(Flame2D_Quad_Block &SolnBlk);

/*************** END FLAME2D ***************************************/

#endif /* _FLAME2D_QUAD_INCLUDED  */
