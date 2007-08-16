/**********************************************************************
 **********************************************************************
 **                                                                  **
 ** File: RteQuad.h                                                  **
 **                                                                  **
 ** Description: The header file defining the 2D RTE Quadrialateral  **
 **              mesh solution classes.                              **
 **                                                                  **
 ** Author: Marc "T-Bone" Charest                                    **
 **                                                                  **
 ** Revision:  Date        Initials   Change                         **
 **            04/03/2007  MRC        Original creation              **
 **                                                                  **
 **********************************************************************
 **********************************************************************/

#ifndef _RTE2D_QUAD_INCLUDED
#define _RTE2D_QUAD_INCLUDED

/* Include 2D Rte state, 2D cell, 2D quadrilateral multiblock 
   grid, quadtree, AMR, and 2D Rte input 
   header files. */

#include "Rte2DState.h"
#include "Rte2DInput.h"
#include "Rte2DTools.h"
#include "../Grid/Cell2D.h"
#include "../Grid/Grid2DQuad.h"
#include "../AMR/QuadTree.h"
#include "../AMR/AMR.h"
#include "../Grid/NASARotor37.h"
#include "../Grid/NASARotor67.h"

/* Include the linear systems header file. */
#include "../Math/LinearSystems.h"

/* Include ICEMCFD input header file. */
#include "../ICEM/ICEMCFD.h"


/* Define the structures and classes. */

#define	NUMBER_OF_RESIDUAL_VECTORS_RTE2D    3

/*!
 * Class: Rte2D_Quad_Block
 *
 * @brief Class definition of the 2D Rte solution blocks.
 *
 * \verbatim
 * Member functions
 *       U      -- Return conserved variable solution for the block
 *                 (cell average).
 *    Grid      -- Return the solution block quadrilateral grid or mesh.
 *      dt      -- Return local time step for the solution block.
 *    dUdt      -- Return the solution block residuals.
 *    dUdx      -- Return the unlimited conservative variable solution
 *                 gradients (x-direction) for the block.
 *    dUdy      -- Return the unlimited conservative variable solution
 *                 gradients (y-direction) for the block.
 *     phi      -- Return the solution slope limiters.
 *   dUdpsi     -- Return the unlimited conservative variable solution
 *                 gradients (azimythal-direction) for the block.
 *   phi_psi    -- Return the solution slope limiter for azimuthal direction.
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
 *      Un      -- Return conserved variable solution at the
 *                 specified node.
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
 * QuadraticSubcellReconstruction -- Performs quadratic subcell solution
 *                                   reconstruction used in adaptive
 *                                   mesh refinement.
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

class Rte2D_Quad_Block{
private:
public:
  //@{ @name Solution state arrays:
  Rte2D_State             **U; //!< Conserved solution state.
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
  Rte2D_State         ***dUdt; //!< Solution residual.
  Rte2D_State            **Uo; //!< Initial solution state.
  static int residual_variable; //!< Static integer that indicates which variable is used for residual calculations.  
  static int Number_of_Residual_Norms; //!< How many Residual norms to plot?
  //@}

  //@{ @name Solution gradient arrays:
  Rte2D_State           **dUdx; //!< Unlimited solution gradient (x-direction).
  Rte2D_State           **dUdy; //!< Unlimited solution gradient (y-direction).
  Rte2D_State            **phi; //!< Solution slope limiter.
  Rte2D_State        **dUdpsi;  // Unlimited solution gradient
                                // (azimuthal-direction).
  Rte2D_State       **phi_psi;  // Solution slope limiter 
                                // (azimuthal-direction).
  //@}

  //@{ @name Boundary solution flux arrays:
  Rte2D_State           *FluxN, //!< North boundary solution flux.
                        *FluxS, //!< South boundary solution flux.
                        *FluxE, //!< East boundary solution flux.
                        *FluxW; //!< West boundary solution flux.
  //@}

  //@{ @name Problem indicator flags:
  int             Axisymmetric; //!< Axisymmetric flow indicator.
  int           Freeze_Limiter; //!< Limiter freezing indicator.
  static int         Flow_Type; //!< Flow type flag (always inviscid).
  static char*   solutionTitle; //!< Solution title info
  //@}

  //@{ @name Boundary condtion reference states:
  Rte2D_State             *UoN, //!< Boundary condition reference states for north boundary.
                          *UoS, //!< Boundary condition reference states for south boundary.
                          *UoE, //!< Boundary condition reference states for east boundary.
                          *UoW; //!< Boundary condition reference states for west boundary.
  //@}
  
  //@{ @name Grid scaling parameter (i.e. scale a 2D grid to a quasi-3D grid):
  double               **Sp;  // scale cell area to volume
  double              **SpN;  // scale north face length to area
  double              **SpS;  // scale south face length to area
  double              **SpE;  // scale east face length to area
  double              **SpW;  // scale west face length to area
  //@}
  

  //@{ @name Boundary conditions:
  double      NorthWallTemp;  // North Wall Temperature
  double      SouthWallTemp;  // South Wall Temperature
  double       EastWallTemp;  // East Wall Temperature  
  double       WestWallTemp;  // West Wall Temperature
  double     NorthWallEmiss;  // North Wall emissivity
  double     SouthWallEmiss;  // South Wall emissivity
  double      EastWallEmiss;  // East Wall emissivity  
  double      WestWallEmiss;  // West Wall emissivity
  //@}


  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  Rte2D_Quad_Block(void) {
      NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
      U = NULL; dt = NULL; dUdt = NULL; 
      dUdx = NULL; dUdy = NULL; phi = NULL; Uo = NULL;
      dUdpsi = NULL; phi_psi = NULL;
      FluxN = NULL; FluxS = NULL; FluxE = NULL; FluxW = NULL;
      UoN = NULL; UoS = NULL; UoE = NULL; UoW = NULL;
      Axisymmetric = 0; Freeze_Limiter = OFF;
      Sp = NULL; SpN = NULL; SpS = NULL; SpE = NULL; SpW = NULL; 
      NorthWallTemp  = ZERO;   SouthWallTemp  = ZERO;  
      EastWallTemp   = ZERO;   WestWallTemp   = ZERO;  
      NorthWallEmiss = ZERO;   SouthWallEmiss = ZERO;  
      EastWallEmiss  = ZERO;   WestWallEmiss  = ZERO;
  }

  //! Copy constructor.
  Rte2D_Quad_Block(const Rte2D_Quad_Block &Soln) {
    NCi = Soln.NCi; ICl = Soln.ICl; ICu = Soln.ICu; 
    NCj = Soln.NCj; JCl = Soln.JCl; JCu = Soln.JCu; Nghost = Soln.Nghost;
    Grid = Soln.Grid; U = Soln.U; dt = Soln.dt; dUdt = Soln.dUdt; 
    dUdx = Soln.dUdx; dUdy = Soln.dUdy; phi = Soln.phi; Uo = Soln.Uo;
    dUdpsi = Soln.dUdpsi; phi_psi = Soln.phi_psi;
    FluxN = Soln.FluxN; FluxS = Soln.FluxS; FluxE = Soln.FluxE; FluxW = Soln.FluxW;
    UoN = Soln.UoN; UoS = Soln.UoS; UoE = Soln.UoE; UoW = Soln.UoW;
    Axisymmetric = 0; Freeze_Limiter = Soln.Freeze_Limiter;
    Sp = Soln.Sp; SpN = Soln.SpN; SpS = Soln.SpS; SpE = Soln.SpE; SpW = Soln.SpW; 
    NorthWallTemp  = Soln.NorthWallTemp;    SouthWallTemp  = Soln.SouthWallTemp;  
    EastWallTemp   = Soln.EastWallTemp;     WestWallTemp   = Soln.WestWallTemp;  
    NorthWallEmiss = Soln.NorthWallEmiss;   SouthWallEmiss = Soln.SouthWallEmiss;  
    EastWallEmiss  = Soln.EastWallEmiss;    WestWallEmiss  = Soln.WestWallEmiss;
  }

  /* Destructor. */
  // ~Rte2D_Quad_Block(void);
  // Use automatically generated destructor.
  //@}

  /* Assignment operator. */
  // Rte2D_Quad_Block operator = (const Rte2D_Quad_Block &Soln);
  // Use automatically generated assignment operator.

  //@{ @name Allocate and deallocate functions.
  //! Allocate memory for structured quadrilateral solution block.
  void allocate(const int Ni, const int Nj, const int Ng);

  //! Deallocate memory for structured quadrilateral solution block.
  void deallocate(void);
  //@}

  //@{ @name Bilinear interplation (Zingg & Yarrow).
  //! Return conserverd solution state at specified node.
  Rte2D_State Un(const int &ii, const int &jj);

  Rte2D_State UnNW(const int &ii, const int &jj); //!< Return conserved solution state at cell NW node.
  Rte2D_State UnNE(const int &ii, const int &jj); //!< Return conserved solution state at cell NE node.
  Rte2D_State UnSE(const int &ii, const int &jj); //!< Return conserved solution state at cell SE node.
  Rte2D_State UnSW(const int &ii, const int &jj); //!< Return conserved solution state at cell SW node.
  //@}

  //@{ @name Member functions for limiter freezing.
  void evaluate_limiters(void); //!< Set flags for limiter evaluation.
  void freeze_limiters(void);   //!< Set flags for limiter freezing.
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const Rte2D_Quad_Block &Soln);
  friend istream &operator >> (istream &in_file, Rte2D_Quad_Block &Soln);
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


  //! Load send boundary ref state message passing buffer.
  int LoadSendBuffer_BC(double *buffer,
			int &buffer_count,
			const int buffer_size,
			const int i_min, 
			const int i_max,
			const int i_inc,
			const int j_min, 
			const int j_max,
			const int j_inc);
  //! Load F2C send boundary ref state message passing buffer.
  int LoadSendBuffer_BC_F2C(double *buffer,
			    int &buffer_count,
			    const int buffer_size,
			    const int i_min, 
			    const int i_max,
			    const int i_inc,
			    const int j_min, 
			    const int j_max,
			    const int j_inc);
  //! Load C2F send boundary ref state message passing buffer.
  int LoadSendBuffer_BC_C2F(double *buffer,
			    int &buffer_count,
			    const int buffer_size,
			    const int i_min, 
			    const int i_max,
			    const int i_inc,
			    const int j_min, 
			    const int j_max,
			    const int j_inc);
  //! Unload receive boundary ref state message passing buffer.
  int UnloadReceiveBuffer_BC(double *buffer,
			     int &buffer_count,
			     const int buffer_size,
			     const int i_min, 
			     const int i_max,
			     const int i_inc,
			     const int j_min, 
			     const int j_max,
			     const int j_inc);
  //! Unload F2C receive boundary ref state message passing buffer.
  int UnloadReceiveBuffer_BC_F2C(double *buffer,
				 int &buffer_count,
				 const int buffer_size,
				 const int i_min, 
				 const int i_max,
				 const int i_inc,
				 const int j_min, 
				 const int j_max,
				 const int j_inc);
  //! Unload C2F receive boundary ref state message passing buffer.
  int UnloadReceiveBuffer_BC_C2F(double *buffer,
				 int &buffer_count,
				 const int buffer_size,
				 const int i_min, 
				 const int i_max,
				 const int i_inc,
				 const int j_min, 
				 const int j_max,
				 const int j_inc);
  //@}

  //@{ @name Member functions required to scale the grid from 2D to quasi-3D for axisymmetric case
  //! Number of solution
  void ScaleGridTo3D( );
  void ScaleGridTo3D( const int AxisymmetricFlag );
  //! Compute scaling parameter for centriod (required for quasi-3D axisymmetric) */
  double Sp_c(const int i, const int j);
  //@}

};


/**************************************************************************
 * Rte2D_Quad_Block::allocate -- Allocate memory.                         *
 **************************************************************************/
inline void Rte2D_Quad_Block::allocate(const int Ni, const int Nj, const int Ng) {
   int i, j, k; /*assert(Ni > 1 && Nj > 1 && Ng > 1);*/ Grid.allocate(Ni, Nj, Ng);
   NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1; 
   NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1; Nghost = Ng;
   U = new Rte2D_State*[NCi];
   dt = new double*[NCi]; dUdt = new Rte2D_State**[NCi];
   dUdx = new Rte2D_State*[NCi]; dUdy = new Rte2D_State*[NCi];
   phi = new Rte2D_State*[NCi]; Uo = new Rte2D_State*[NCi];
   dUdpsi = new Rte2D_State*[NCi]; phi_psi = new Rte2D_State*[NCi];
   Sp = new double*[NCi]; 
   SpN = new double*[NCi]; SpS = new double*[NCi]; 
   SpE = new double*[NCi]; SpW = new double*[NCi]; 
   for ( i = 0; i <= NCi-1 ; ++i ) {
      U[i] = new Rte2D_State[NCj];
      dt[i] = new double[NCj]; dUdt[i] = new Rte2D_State*[NCj];
      for ( j = 0; j <= NCj-1 ; ++j ) 
        { dUdt[i][j] = new Rte2D_State[NUMBER_OF_RESIDUAL_VECTORS_RTE2D]; }
      dUdx[i] = new Rte2D_State[NCj]; dUdy[i] = new Rte2D_State[NCj];
      phi[i] = new Rte2D_State[NCj]; Uo[i] = new Rte2D_State[NCj];
      dUdpsi[i] = new Rte2D_State[NCj];  phi_psi[i] = new Rte2D_State[NCj];
      Sp[i] = new double[NCj]; 
      SpN[i] = new double[NCj]; SpS[i] = new double[NCj]; 
      SpE[i] = new double[NCj]; SpW[i] = new double[NCj]; 
   } /* endfor */
   FluxN = new Rte2D_State[NCi]; FluxS = new Rte2D_State[NCi];
   FluxE = new Rte2D_State[NCj]; FluxW = new Rte2D_State[NCj];
   UoN = new Rte2D_State[NCi]; UoS = new Rte2D_State[NCi];
   UoE = new Rte2D_State[NCj]; UoW = new Rte2D_State[NCj];
   // Set the solution residuals, gradients, limiters, and other values to zero.
   for (j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
      for ( i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
          for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_RTE2D-1 ; ++k ) {
	    dUdt[i][j][k].Zero();
          } /* endfor */
	  dUdx[i][j].Zero(); dUdy[i][j].Zero();
	  phi[i][j].Zero(); Uo[i][j].Zero();
	  dUdpsi[i][j].Zero();  phi_psi[i][j].Zero();
	  dt[i][j] = ZERO;
	  Sp[i][j] = ONE; 
	  SpN[i][j] = ONE; SpS[i][j] = ONE; 
	  SpE[i][j] = ONE; SpW[i][j] = ONE; 
      } /* endfor */
   } /* endfor */

}

/**************************************************************************
 * Rte2D_Quad_Block::deallocate -- Deallocate memory.                   *
 **************************************************************************/
inline void Rte2D_Quad_Block::deallocate(void) {
   int i, j; Grid.deallocate(); 
   for ( i = 0; i <= NCi-1 ; ++i ) {
      delete []U[i]; U[i] = NULL;
      delete []dt[i]; dt[i] = NULL; 
      for ( j = 0; j <= NCj-1 ; ++j ) { delete []dUdt[i][j]; dUdt[i][j] = NULL; }
      delete []dUdt[i]; dUdt[i] = NULL;
      delete []dUdx[i]; dUdx[i] = NULL; delete []dUdy[i]; dUdy[i] = NULL;
      delete []phi[i]; phi[i] = NULL; delete []Uo[i]; Uo[i] = NULL;
      delete []dUdpsi[i]; dUdpsi[i] = NULL; delete []phi_psi[i]; phi_psi[i] = NULL;
      delete[] Sp[i];  Sp[i] = NULL; 
      delete[] SpN[i]; SpN[i] = NULL; 
      delete[] SpS[i]; SpS[i] = NULL; 
      delete[] SpE[i]; SpE[i] = NULL; 
      delete[] SpW[i]; SpW[i] = NULL; 
   } /* endfor */
   delete []U; U = NULL;
   delete []dt; dt = NULL; delete []dUdt; dUdt = NULL;
   delete []dUdx; dUdx = NULL; delete []dUdy; dUdy = NULL; 
   delete []phi; phi = NULL; delete []Uo; Uo = NULL;
   delete []dUdpsi; dUdpsi = NULL; delete []phi_psi; phi_psi = NULL;
   delete []FluxN; FluxN = NULL; delete []FluxS; FluxS = NULL;
   delete []FluxE; FluxE = NULL; delete []FluxW; FluxW = NULL;
   delete []UoN; UoN = NULL; delete []UoS; UoS = NULL;
   delete []UoE; UoE = NULL; delete []UoW; UoW = NULL;
   delete[] Sp;  Sp  = NULL; 
   delete[] SpN; SpN = NULL; 
   delete[] SpS; SpS = NULL; 
   delete[] SpE; SpE = NULL; 
   delete[] SpW; SpW = NULL; 
   NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
}

/**************************************************************************
 * Rte2D_Quad_Block::Un -- Node conservative solution.                    *
 **************************************************************************/
inline Rte2D_State Rte2D_Quad_Block::Un(const int &ii, const int &jj) {
 
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
 * Rte2D_Quad_Block::Un?? -- Get cell node conserved solution states.     *
 **************************************************************************/
inline Rte2D_State Rte2D_Quad_Block::UnNW(const int &ii, const int &jj) {
  return (Un(ii, jj+1));
}

inline Rte2D_State Rte2D_Quad_Block::UnNE(const int &ii, const int &jj) {
  return (Un(ii+1, jj+1));
}

inline Rte2D_State Rte2D_Quad_Block::UnSE(const int &ii, const int &jj) {
  return (Un(ii+1, jj));
}

inline Rte2D_State Rte2D_Quad_Block::UnSW(const int &ii, const int &jj) {
  return (Un(ii, jj));
}

/**************************************************************************
 * Rte2D_Quad_Block::evaluate_limiters -- Set flag to evaluate limiters.*
 **************************************************************************/
inline void Rte2D_Quad_Block::evaluate_limiters(void) {
  Freeze_Limiter = OFF; 
}

/**************************************************************************
 * Rte2D_Quad_Block::freeze_limiters -- Set flag to freeze limiters.    *
 **************************************************************************/
inline void Rte2D_Quad_Block::freeze_limiters(void) {
  Freeze_Limiter = ON; 
}

/**************************************************************************
 * Rte2D_Quad_Block -- Input-output operators.                          *
 **************************************************************************/
inline ostream &operator << (ostream &out_file,
			     const Rte2D_Quad_Block &SolnBlk) {
  int i, j; 
  out_file << SolnBlk.Grid;
  out_file << SolnBlk.NCi << " " << SolnBlk.ICl << " " 
	   << SolnBlk.ICu << " " << SolnBlk.Nghost << "\n";
  out_file << SolnBlk.NCj << " " << SolnBlk.JCl << " " << SolnBlk.JCu << "\n";
  out_file << SolnBlk.Axisymmetric << "\n";
  out_file << SolnBlk.NorthWallTemp << "\n";
  out_file << SolnBlk.SouthWallTemp << "\n";
  out_file << SolnBlk.EastWallTemp << "\n";
  out_file << SolnBlk.WestWallTemp << "\n";
  out_file << SolnBlk.NorthWallEmiss << "\n";
  out_file << SolnBlk.SouthWallEmiss << "\n";
  out_file << SolnBlk.EastWallEmiss << "\n";
  out_file << SolnBlk.WestWallEmiss << "\n";
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
			     Rte2D_Quad_Block &SolnBlk) {
  int i, j, k, ni, il, iu, ng, nj, jl, ju;
  Grid2D_Quad_Block New_Grid; 
  
  in_file >> New_Grid; 
  in_file.setf(ios::skipws);
  in_file >> ni >> il >> iu >> ng; in_file >> nj >> jl >> ju;
  in_file >> SolnBlk.Axisymmetric;
  in_file >> SolnBlk.NorthWallTemp;
  in_file >> SolnBlk.SouthWallTemp;
  in_file >> SolnBlk.EastWallTemp;
  in_file >> SolnBlk.WestWallTemp;
  in_file >> SolnBlk.NorthWallEmiss;
  in_file >> SolnBlk.SouthWallEmiss;
  in_file >> SolnBlk.EastWallEmiss;
  in_file >> SolnBlk.WestWallEmiss;
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
         for ( k = 0 ; k <= NUMBER_OF_RESIDUAL_VECTORS_RTE2D-1 ; ++k ) {
	   SolnBlk.dUdt[i][j][k].Zero();
         } /* endfor */
	 SolnBlk.dUdx[i][j].Zero();
	 SolnBlk.dUdy[i][j].Zero();
	 SolnBlk.phi[i][j].Zero();
	 SolnBlk.Uo[i][j].Zero();
	 SolnBlk.dt[i][j] = ZERO;
	 SolnBlk.dUdpsi[i][j].Zero();
	 SolnBlk.phi_psi[i][j].Zero();
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
 * Rte2D_Quad_Block::NumVar -- Returns number of state variables.              *
 *******************************************************************************/
inline int Rte2D_Quad_Block::NumVar(void) {
  return (int(U[0][0].NUM_VAR_RTE2D));
}

/*******************************************************************************
 * Rte2D_Quad_Block::LoadSendBuffer -- Loads send message buffer.              *
 *******************************************************************************/
inline int Rte2D_Quad_Block::LoadSendBuffer(double *buffer,
                                              int &buffer_count,
                                              const int buffer_size,
                                              const int i_min, 
                                              const int i_max,
                                              const int i_inc,
                                              const int j_min, 
                                              const int j_max,
                                              const int j_inc) {
  int NUM_VAR_RTE2D = NumVar();
  int i, j, k;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = U[i][j][k];	   
        } /* endfor */
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * Rte2D_Quad_Block::LoadSendBuffer_F2C -- Loads send message buffer for       *
 *                                           fine to coarse block message      *
 *                                           passing.                          *
 *******************************************************************************/
inline int Rte2D_Quad_Block::LoadSendBuffer_F2C(double *buffer,
                                                  int &buffer_count,
                                                  const int buffer_size,
                                                  const int i_min, 
                                                  const int i_max,
                                                  const int i_inc,
                                                  const int j_min, 
                                                  const int j_max,
                                                  const int j_inc) {
  int NUM_VAR_RTE2D = NumVar();
  int i, j, k;
  for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (Grid.Cell[i  ][j  ].A*Sp[i  ][j  ]*U[i  ][j  ][k]+
                                   Grid.Cell[i+1][j  ].A*Sp[i+1][j  ]*U[i+1][j  ][k]+
                                   Grid.Cell[i  ][j+1].A*Sp[i  ][j+1]*U[i  ][j+1][k]+
                                   Grid.Cell[i+1][j+1].A*Sp[i+1][j+1]*U[i+1][j+1][k])/
                                  (Grid.Cell[i  ][j  ].A*Sp[i  ][j  ]+
                                   Grid.Cell[i+1][j  ].A*Sp[i+1][j  ]+
                                   Grid.Cell[i  ][j+1].A*Sp[i  ][j+1]+
                                   Grid.Cell[i+1][j+1].A*Sp[i+1][j+1]);
        } /* endfor */
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * Rte2D_Quad_Block::LoadSendBuffer_C2F -- Loads send message buffer for     *
 *                                           coarse to fine block message      *
 *                                           passing.                          *
 *******************************************************************************/
inline int Rte2D_Quad_Block::LoadSendBuffer_C2F(double *buffer,
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
  int NUM_VAR_RTE2D = NumVar();
  int i, j, k;
  Vector2D dX;
  Rte2D_State Ufine;

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
              Ufine = U[i][j_min] +
                      (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dUdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i][j_min] +
                      (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dUdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i][j_min] +
                      (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dUdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) { 
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
              Ufine = U[i][j_min] +
                      (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dUdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i][j_min] +
                      (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dUdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i][j_min] +
                      (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dUdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i][j_min] +
                      (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dUdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i][j_min] +
                      (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dUdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) { 
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
              Ufine = U[i][j_min] +
                      (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dUdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) { 
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
              Ufine = U[i][j_min] +
                      (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dUdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i][j_min] +
                      (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dUdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i][j_min] +
                      (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dUdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i][j_min] +
                      (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dUdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i][j_min] +
                      (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dUdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) { 
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
              Ufine = U[i][j_min] +
                      (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dUdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i][j_min] +
                      (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                      (phi[i][j_min]^dUdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i_min][j] +
                      (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dUdy[i_min][j])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i_min][j] +
                      (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dUdy[i_min][j])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i_min][j] +
                      (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dUdy[i_min][j])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) { 
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
              Ufine = U[i_min][j] +
                      (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dUdy[i_min][j])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i_min][j] +
                      (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dUdy[i_min][j])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i_min][j] +
                      (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dUdy[i_min][j])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i_min][j] +
                      (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dUdy[i_min][j])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i_min][j] +
                      (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dUdy[i_min][j])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) { 
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
              Ufine = U[i_min][j] +
                      (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dUdy[i_min][j])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) { 
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
              Ufine = U[i_min][j] +
                      (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dUdy[i_min][j])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i_min][j] +
                      (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dUdy[i_min][j])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i_min][j] +
                      (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dUdy[i_min][j])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i_min][j] +
                      (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dUdy[i_min][j])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i_min][j] +
                      (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dUdy[i_min][j])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) { 
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
              Ufine = U[i_min][j] +
                      (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dUdy[i_min][j])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
              Ufine = U[i_min][j] +
                      (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                      (phi[i_min][j]^dUdy[i_min][j])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
 * Rte2D_Quad_Block::UnloadReceiveBuffer -- Unloads receive message buffer.  *
 *******************************************************************************/
inline int Rte2D_Quad_Block::UnloadReceiveBuffer(double *buffer,
                                                   int &buffer_count,
                                                   const int buffer_size,
                                                   const int i_min, 
                                                   const int i_max,
                                                   const int i_inc,
                                                   const int j_min, 
                                                   const int j_max,
                                                   const int j_inc) {
  int NUM_VAR_RTE2D = NumVar();
  int i, j;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {

       /* Changed from Rte2D */
       for ( int k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
	 buffer_count = buffer_count + 1;
	 if (buffer_count >= buffer_size) return(1);    
	 U[i][j][k] = buffer[buffer_count];
       }
       /* End Changeed from Rte2D */

     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * Rte2D_Quad_Block::UnloadReceiveBuffer_F2C -- Unloads receive message      *
 *                                                buffer for fine to coarse    *
 *                                                block message passing.       *
 *******************************************************************************/
inline int Rte2D_Quad_Block::UnloadReceiveBuffer_F2C(double *buffer,
                                                       int &buffer_count,
                                                       const int buffer_size,
                                                       const int i_min, 
                                                       const int i_max,
                                                       const int i_inc,
                                                       const int j_min, 
                                                       const int j_max,
                                                       const int j_inc) {
  int NUM_VAR_RTE2D = NumVar();
  int i, j;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {

       /* Changed from Rte2D */
       for ( int k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
	 buffer_count = buffer_count + 1;
	 if (buffer_count >= buffer_size) return(1);    
	 U[i][j][k] = buffer[buffer_count];
       }
       /* End Changeed from Rte2D */

     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * Rte2D_Quad_Block::UnloadReceiveBuffer_C2F -- Unloads receive message      *
 *                                                buffer for coarse to fine    *
 *                                                block message passing.       *
 *******************************************************************************/
inline int Rte2D_Quad_Block::UnloadReceiveBuffer_C2F(double *buffer,
                                                       int &buffer_count,
                                                       const int buffer_size,
                                                       const int i_min, 
                                                       const int i_max,
                                                       const int i_inc,
                                                       const int j_min, 
                                                       const int j_max,
                                                       const int j_inc) {
  int NUM_VAR_RTE2D = NumVar();
  int i, j;
  for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {

       /* Changed from Rte2D */
       for ( int k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
	 buffer_count = buffer_count + 1;
	 if (buffer_count >= buffer_size) return(1);    
	 U[i][j][k] = buffer[buffer_count];
       }
       /* End Changeed from Rte2D */

     } /* endfor */
  } /* endfor */
  return(0);
}

/**************************************************************************
 * Rte2D_Quad_Block::SubcellReconstruction --                           *
 *               Performs the subcell reconstruction of solution state    *
 *               within a given cell (i,j) of the computational mesh for  *
 *               the specified quadrilateral solution block.              *
 **************************************************************************/
inline void Rte2D_Quad_Block::SubcellReconstruction(const int i, 
                                                      const int j,
                                                      const int Limiter) {

  int NUM_VAR_RTE2D = NumVar();
  int n, n2, n_pts, i_index[8], j_index[8];
  double u0Min, u0Max, uQuad[4], phi_n;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  Rte2D_State DU, DUDx_ave, DUDy_ave;

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
               Grid.BCtypeW[j] == BC_GRAY_WALL) {
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
               Grid.BCtypeE[j] == BC_GRAY_WALL) {
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
               Grid.BCtypeS[i] == BC_GRAY_WALL) {
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
               Grid.BCtypeN[i] == BC_GRAY_WALL) {
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
      DUDx_ave.Zero();
      DUDy_ave.Zero();
      DxDx_ave = ZERO;
      DxDy_ave = ZERO;
      DyDy_ave = ZERO;
  
      for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
          dX = Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - 
               Grid.Cell[i][j].Xc;
          DU = U[ i_index[n2] ][ j_index[n2] ] - 
               U[i][j];
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
      dUdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                   (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
      dUdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                   (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
  
      // Calculate slope limiters. 
      if (!Freeze_Limiter) {
         for ( n = 1 ; n <= NUM_VAR_RTE2D ; ++n ) {
            u0Min = U[i][j][n];
            u0Max = u0Min;
            for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
               u0Min = min(u0Min, U[ i_index[n2] ][ j_index[n2] ][n]);
               u0Max = max(u0Max, U[ i_index[n2] ][ j_index[n2] ][n]);
            } /* endfor */
    
            dX = Grid.xfaceE(i, j)-Grid.Cell[i][j].Xc;
            uQuad[0] = U[i][j][n] + 
                       dUdx[i][j][n]*dX.x +
                       dUdy[i][j][n]*dX.y ;
            dX = Grid.xfaceW(i, j)-Grid.Cell[i][j].Xc;
            uQuad[1] = U[i][j][n] + 
                       dUdx[i][j][n]*dX.x +
                       dUdy[i][j][n]*dX.y ;
            dX = Grid.xfaceN(i, j)-Grid.Cell[i][j].Xc;
            uQuad[2] = U[i][j][n] + 
                       dUdx[i][j][n]*dX.x +
                       dUdy[i][j][n]*dX.y ;
            dX = Grid.xfaceS(i, j)-Grid.Cell[i][j].Xc;
            uQuad[3] = U[i][j][n] + 
                       dUdx[i][j][n]*dX.x +
                       dUdy[i][j][n]*dX.y ;
    
            switch(Limiter) {
              case LIMITER_ONE :
                phi_n = ONE;
                break;
              case LIMITER_ZERO :
                phi_n = ZERO;
                break;
              case LIMITER_BARTH_JESPERSEN :
                phi_n = Limiter_BarthJespersen(uQuad, U[i][j][n], 
                                               u0Min, u0Max, 4);
                break;
              case LIMITER_VENKATAKRISHNAN :
                phi_n = Limiter_Venkatakrishnan(uQuad, U[i][j][n], 
                                                u0Min, u0Max, 4);
                break;
              case LIMITER_VANLEER :
                phi_n = Limiter_VanLeer(uQuad, U[i][j][n], 
                                        u0Min, u0Max, 4);
                break;
              case LIMITER_VANALBADA :
                phi_n = Limiter_VanAlbada(uQuad, U[i][j][n], 
                                          u0Min, u0Max, 4);
                break;
              default:
                phi_n = Limiter_BarthJespersen(uQuad, U[i][j][n], 
                                               u0Min, u0Max, 4);
                break;
            } /* endswitch */

	    phi[i][j][n] = phi_n;

         } /* endfor */
      } /* endif */
  } else {
      dUdx[i][j].Zero();
      dUdy[i][j].Zero(); 
      phi[i][j].Zero();
  } /* endif */
    
}

/*******************************************************************************
 * Rte2D_Quad_Block::LoadSendBuffer_Flux_F2C -- Loads send message buffer for*
 *                                                fine to coarse block message *
 *                                                passing of conservative      *
 *                                                solution fluxes.             *
 *******************************************************************************/
inline int Rte2D_Quad_Block::LoadSendBuffer_Flux_F2C(double *buffer,
						     int &buffer_count,
						     const int buffer_size,
						     const int i_min, 
						     const int i_max,
						     const int i_inc,
						     const int j_min, 
						     const int j_max,
						     const int j_inc) {
  int NUM_VAR_RTE2D = NumVar();
  int i, j, k;
  if (j_min == j_max && j_min == JCl) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxS[i  ][k]+
                                   FluxS[i+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxN[i  ][k]+
                                   FluxN[i+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
     for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
        for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (FluxW[j][k]+
                                   FluxW[j+1][k]);
        } /* endfor */
     } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
     for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
        for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
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
 * Rte2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C -- Unloads receive message   *
 *                                                buffer for fine to coarse    *
 *                                                block message passing of     *
 *                                                conservative solution fluxes.*
 *******************************************************************************/
inline int Rte2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C(double *buffer,
							  int &buffer_count,
							  const int buffer_size,
							  const int i_min, 
							  const int i_max,
							  const int i_inc,
							  const int j_min, 
							  const int j_max,
							  const int j_inc) {
  int NUM_VAR_RTE2D = NumVar();
  int i, j;
  if (j_min == j_max && j_min == JCl) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
       
       /* Changed for Rte2D */
       for ( int k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
	 buffer_count = buffer_count + 1;
	 if (buffer_count >= buffer_size) return(1);    
	 FluxS[i][k] = - buffer[buffer_count] - FluxS[i][k];
       }
       /* End Change */
 
     } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {

       /* Changed for Rte2D */
      for ( int k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	FluxN[i][k] = - buffer[buffer_count] - FluxN[i][k];
      }
       /* End Change */

     } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
     for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {

       /* Changed for Rte2D */
      for ( int k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	FluxW[j][k] = - buffer[buffer_count] - FluxW[j][k];
      }
       /* End Change */

     } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
     for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {

       /* Changed for Rte2D */
      for ( int k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	FluxE[j][k] = - buffer[buffer_count] - FluxE[j][k];
      }
       /* End Change */

     } /* endfor */
  } /* endif */
  return(0);
}

/*******************************************************************************
 * Rte2D_Quad_Block::LoadSendBuffer_BC -- Loads send message buffer for passing*
 *                                        of boundary ref state.               *
 *******************************************************************************/
inline int Rte2D_Quad_Block::LoadSendBuffer_BC(double *buffer,
					       int &buffer_count,
					       const int buffer_size,
					       const int i_min, 
					       const int i_max,
					       const int i_inc,
					       const int j_min, 
					       const int j_max,
					       const int j_inc) {
  int NUM_VAR_RTE2D = NumVar();
  int i, j, k;
  if (j_min == j_max && j_min == JCl) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = UoS[i][k];
        } /* endfor */
     } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = UoN[i][k];
        } /* endfor */
     } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
    for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
        for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = UoW[j][k];
        } /* endfor */
     } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
    for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
        for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = UoE[j][k];
        } /* endfor */
     } /* endfor */
  } /* endif */



//   cout<<"\n EULER i_min "<<i_min<<" i_max "<<i_max<<" i_inc "<<i_inc
//       <<" j_min "<<j_min<<" j_max "<<j_max<<" j_inc "<<j_inc
//       <<" buffersize "<<buffer_size;

  return(0);
}

/*******************************************************************************
 * Rte2D_Quad_Block::LoadSendBuffer_F2C_BC -- Loads send message buffer for    *
 *                                           fine to coarse block message      *
 *                                           passing of boundary ref states.   *
 *******************************************************************************/
inline int Rte2D_Quad_Block::LoadSendBuffer_BC_F2C(double *buffer,
						   int &buffer_count,
						   const int buffer_size,
						   const int i_min, 
						   const int i_max,
						   const int i_inc,
						   const int j_min, 
						   const int j_max,
						   const int j_inc) {
  int NUM_VAR_RTE2D = NumVar();
  int i, j, k;

  if (j_min == j_max && j_min == JCl) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (Grid.Cell[i  ][j_min].A*Sp[i  ][j_min]*UoS[i  ][k]+
                                   Grid.Cell[i+1][j_min].A*Sp[i+1][j_min]*UoS[i+1][k] ) /
                                  (Grid.Cell[i  ][j_min].A*Sp[i  ][j_min]+
                                   Grid.Cell[i+1][j_min].A*Sp[i+1][j_min]);
       } /* endfor */
     } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (Grid.Cell[i  ][j_min].A*Sp[i  ][j_min]*UoN[i  ][k]+
                                   Grid.Cell[i+1][j_min].A*Sp[i+1][j_min]*UoN[i+1][k] ) /
                                  (Grid.Cell[i  ][j_min].A*Sp[i  ][j_min]+
                                   Grid.Cell[i+1][j_min].A*Sp[i+1][j_min]);
        } /* endfor */
     } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
     for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
        for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (Grid.Cell[i_min][j  ].A*Sp[i_min][j  ]*UoW[j  ][k]+
                                   Grid.Cell[i_min][j+1].A*Sp[i_min][j+1]*UoW[j+1][k] ) /
                                  (Grid.Cell[i_min][j  ].A*Sp[i_min][j  ]+
                                   Grid.Cell[i_min][j+1].A*Sp[i_min][j+1]);
        } /* endfor */
     } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
     for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
        for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = (Grid.Cell[i_min][j  ].A*Sp[i_min][j  ]*UoE[j  ][k]+
                                   Grid.Cell[i_min][j+1].A*Sp[i_min][j+1]*UoE[j+1][k] ) /
                                  (Grid.Cell[i_min][j  ].A*Sp[i_min][j  ]+
                                   Grid.Cell[i_min][j+1].A*Sp[i_min][j+1]);
        } /* endfor */
     } /* endfor */
  } /* endif */



  return(0);
}

/*******************************************************************************
 * Rte2D_Quad_Block::LoadSendBuffer_BC_C2F -- Loads send message buffer for    *
 *                                           coarse to fine block message      *
 *                                           passing for boundary ref states.  *
 *******************************************************************************/
inline int Rte2D_Quad_Block::LoadSendBuffer_BC_C2F(double *buffer,
                                                  int &buffer_count,
                                                  const int buffer_size,
                                                  const int i_min, 
                                                  const int i_max,
                                                  const int i_inc,
                                                  const int j_min, 
                                                  const int j_max,
                                                  const int j_inc) {
  int NUM_VAR_RTE2D = NumVar();
  int i, j, k;
  if (j_min == j_max && j_min == JCl) {
    for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = UoS[i][k];
        } /* endfor */
     } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
    for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
        for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = UoN[i][k];
        } /* endfor */
     } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
    for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
        for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = UoW[j][k];
        } /* endfor */
     } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
    for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
        for ( k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
           buffer[buffer_count] = UoE[j][k];
        } /* endfor */
     } /* endfor */
  } /* endif */

  return(0);
}

/*******************************************************************************
 * Rte2D_Quad_Block::UnloadReceiveBuffer -- Unloads receive message buffer.  *
 *******************************************************************************/
inline int Rte2D_Quad_Block::UnloadReceiveBuffer_BC(double *buffer,
                                                   int &buffer_count,
                                                   const int buffer_size,
                                                   const int i_min, 
                                                   const int i_max,
                                                   const int i_inc,
                                                   const int j_min, 
                                                   const int j_max,
                                                   const int j_inc) {
  int NUM_VAR_RTE2D = NumVar();
  int i, j;

  if (j_min == j_max && j_min == JCl) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
       
       /* Changed for Rte2D */
       for ( int k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
	 buffer_count = buffer_count + 1;
	 if (buffer_count >= buffer_size) return(1);    
	 UoS[i][k] = buffer[buffer_count];
       }
       /* End Change */
 
     } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {

       /* Changed for Rte2D */
      for ( int k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	UoN[i][k] = buffer[buffer_count];
      }
       /* End Change */

     } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
    for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {

       /* Changed for Rte2D */
      for ( int k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	UoW[j][k] = buffer[buffer_count];
      }
       /* End Change */

     } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
    for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
    
       /* Changed for Rte2D */
      for ( int k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	UoE[j][k] = buffer[buffer_count];
      }
       /* End Change */

     } /* endfor */
  } /* endif */


  return(0);
}

/*******************************************************************************
 * Rte2D_Quad_Block::UnloadReceiveBuffer_BC_F2C -- Unloads receive message     *
 *                                                buffer for fine to coarse    *
 *                                                block message passing.       *
 *                                                for boundary ref states.     *
 *******************************************************************************/
inline int Rte2D_Quad_Block::UnloadReceiveBuffer_BC_F2C(double *buffer,
                                                       int &buffer_count,
                                                       const int buffer_size,
                                                       const int i_min, 
                                                       const int i_max,
                                                       const int i_inc,
                                                       const int j_min, 
                                                       const int j_max,
                                                       const int j_inc) {
  int NUM_VAR_RTE2D = NumVar();
  int i, j;
  if (j_min == j_max && j_min == JCl) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
       
       /* Changed for Rte2D */
       for ( int k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
	 buffer_count = buffer_count + 1;
	 if (buffer_count >= buffer_size) return(1);    
	 UoS[i][k] = buffer[buffer_count];
       }
       /* End Change */
 
     } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {

       /* Changed for Rte2D */
      for ( int k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	UoN[i][k] = buffer[buffer_count];
      }
       /* End Change */

     } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
     for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {

       /* Changed for Rte2D */
      for ( int k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	UoW[j][k] = buffer[buffer_count];
      }
       /* End Change */

     } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
     for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {

       /* Changed for Rte2D */
      for ( int k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	UoE[j][k] = buffer[buffer_count];
      }
       /* End Change */

     } /* endfor */
  } /* endif */

  return(0);
}

/*******************************************************************************
 * Rte2D_Quad_Block::UnloadReceiveBuffer_BC_C2F -- Unloads receive message      *
 *                                                buffer for coarse to fine    *
 *                                                block message passing.       *
 *                                                for boundary ref states.     *
 *******************************************************************************/
inline int Rte2D_Quad_Block::UnloadReceiveBuffer_BC_C2F(double *buffer,
                                                       int &buffer_count,
                                                       const int buffer_size,
                                                       const int i_min, 
                                                       const int i_max,
                                                       const int i_inc,
                                                       const int j_min, 
                                                       const int j_max,
                                                       const int j_inc) {
  int NUM_VAR_RTE2D = NumVar();
  int i, j;

  if (j_min == j_max && j_min == JCl) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
       
       /* Changed for Rte2D */
       for ( int k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
	 buffer_count = buffer_count + 1;
	 if (buffer_count >= buffer_size) return(1);    
	 UoS[i][k] = buffer[buffer_count];
       }
       /* End Change */
 
     } /* endfor */
  } else if (j_min == j_max && j_min == JCu) {
     for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {

       /* Changed for Rte2D */
      for ( int k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	UoN[i][k] = buffer[buffer_count];
      }
       /* End Change */

     } /* endfor */
  } else if (i_min == i_max && i_min == ICl) {
    for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {

       /* Changed for Rte2D */
      for ( int k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	UoW[j][k] = buffer[buffer_count];
      }
       /* End Change */

     } /* endfor */
  } else if (i_min == i_max && i_min == ICu) {
    for ( j  = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
    
       /* Changed for Rte2D */
      for ( int k = 1 ; k <= NUM_VAR_RTE2D; ++ k) {
	buffer_count = buffer_count + 1;
	if (buffer_count >= buffer_size) return(1);    
	UoE[j][k] = buffer[buffer_count];
      }
       /* End Change */

     } /* endfor */
  } /* endif */


  return(0);
}




/*******************************************************************************
 * Rte2D_Quad_Block::ScaleGrid -- Scale the grid to quasi-3D for axisymmetric  *
 *                                case.                                        *
 *******************************************************************************/
inline void Rte2D_Quad_Block::ScaleGridTo3D() {


  switch (Axisymmetric) {

  // PLANAR RADIATION
  case PLANAR:
  default:
    
    for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
      for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
	Sp[i][j] = ONE; 
	SpN[i][j] = ONE; SpS[i][j] = ONE; 
	SpE[i][j] = ONE; SpW[i][j] = ONE; 
      } /* endfor */
    } /* endfor */
    break;


  // AXISYMMETRIC X RADIATION
  case AXISYMMETRIC_X:
    for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
      for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
	Sp[i][j]  = fabs(Grid.Cell[i][j].Xc.x); 
	SpN[i][j] = fabs(Grid.xfaceN(i,j).x); 
	SpS[i][j] = fabs(Grid.xfaceS(i,j).x); 
	SpE[i][j] = fabs(Grid.xfaceE(i,j).x); 
	SpW[i][j] = fabs(Grid.xfaceW(i,j).x); 
      } /* endfor */
    } /* endfor */
    break;


  // AXISYMMETRIC Y RADIATION
  case AXISYMMETRIC_Y:
    for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
      for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
	Sp[i][j]  = fabs(Grid.Cell[i][j].Xc.y); 
	SpN[i][j] = fabs(Grid.xfaceN(i,j).y); 
	SpS[i][j] = fabs(Grid.xfaceS(i,j).y); 
	SpE[i][j] = fabs(Grid.xfaceE(i,j).y); 
	SpW[i][j] = fabs(Grid.xfaceW(i,j).y); 
      } /* endfor */
    } /* endfor */
    break;

  } /* endswitch Axisymmetric */


}

inline void Rte2D_Quad_Block::ScaleGridTo3D(const int AxisymmetricFlag) {


  switch (AxisymmetricFlag) {

  // PLANAR RADIATION
  case PLANAR:
  default:
    
    for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
      for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
	Sp[i][j] = ONE; 
	SpN[i][j] = ONE; SpS[i][j] = ONE; 
	SpE[i][j] = ONE; SpW[i][j] = ONE; 
      } /* endfor */
    } /* endfor */
    break;


  // AXISYMMETRIC X RADIATION
  case AXISYMMETRIC_X:
    for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
      for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
	Sp[i][j]  = fabs(Grid.Cell[i][j].Xc.x); 
	SpN[i][j] = fabs(Grid.xfaceN(i,j).x); 
	SpS[i][j] = fabs(Grid.xfaceS(i,j).x); 
	SpE[i][j] = fabs(Grid.xfaceE(i,j).x); 
	SpW[i][j] = fabs(Grid.xfaceW(i,j).x); 
      } /* endfor */
    } /* endfor */
    break;


  // AXISYMMETRIC Y RADIATION
  case AXISYMMETRIC_Y:
    for (int j  = JCl-Nghost ; j <= JCu+Nghost ; ++j ) {
      for (int i = ICl-Nghost ; i <= ICu+Nghost ; ++i ) {
	Sp[i][j]  = fabs(Grid.Cell[i][j].Xc.y); 
	SpN[i][j] = fabs(Grid.xfaceN(i,j).y); 
	SpS[i][j] = fabs(Grid.xfaceS(i,j).y); 
	SpE[i][j] = fabs(Grid.xfaceE(i,j).y); 
	SpW[i][j] = fabs(Grid.xfaceW(i,j).y); 
      } /* endfor */
    } /* endfor */
    break;

  } /* endswitch Axisymmetric */


}


/*******************************************************************************
 * Rte2D_Quad_Block::Sp_c -- Evaluate the current quasi-3D grid scaling        *
 *                           parameter for the cell centroid.                  *
 *******************************************************************************/
inline double Rte2D_Quad_Block::Sp_c( const int i, const int j ) {


  switch (Axisymmetric) {

  // PLANAR RADIATION
  case PLANAR:
  default:
    return (ONE);

  // AXISYMMETRIC X RADIATION
  case AXISYMMETRIC_X:
    return fabs(Grid.centroid(i,j).x);

  // AXISYMMETRIC Y RADIATION
  case AXISYMMETRIC_Y:
    return fabs(Grid.centroid(i,j).y);

  } /* endswitch Axisymmetric */


}



/**************************************************************************
 * Rte2D_Quad_Block -- Single Block External Subroutines.                 *
 **************************************************************************/

extern void Write_Solution_Block(Rte2D_Quad_Block &SolnBlk,
	                         ostream &Out_File);

extern void Read_Solution_Block(Rte2D_Quad_Block &SolnBlk,
	                        istream &In_File);

extern void Broadcast_Solution_Block(Rte2D_Quad_Block &SolnBlk);

#ifdef _MPI_VERSION
extern void Broadcast_Solution_Block(Rte2D_Quad_Block &SolnBlk,
                                     MPI::Intracomm &Communicator, 
                                     const int Source_CPU);
#endif

extern void Copy_Solution_Block(Rte2D_Quad_Block &SolnBlk1,
		                Rte2D_Quad_Block &SolnBlk2);

extern int Prolong_Solution_Block(Rte2D_Quad_Block &SolnBlk_Fine,
				  Rte2D_Quad_Block &SolnBlk_Original,
				  const int Sector);

extern int Restrict_Solution_Block(Rte2D_Quad_Block &SolnBlk_Coarse,
				   Rte2D_Quad_Block &SolnBlk_Original_SW,
				   Rte2D_Quad_Block &SolnBlk_Original_SE,
				   Rte2D_Quad_Block &SolnBlk_Original_NW,
				   Rte2D_Quad_Block &SolnBlk_Original_NE);

extern void Output_Tecplot(Rte2D_Quad_Block &SolnBlk,
			   Rte2D_Input_Parameters &IP,
		           const int Number_of_Time_Steps,
                           const double &Time,
                           const int Block_Number,
                           const int Output_Title,
	                   ostream &Out_File);

extern void Output_Cells_Tecplot(Rte2D_Quad_Block &SolnBlk,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_Nodes_Tecplot(Rte2D_Quad_Block &SolnBlk,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_Gradients_Tecplot(Rte2D_Quad_Block &SolnBlk,
				     const int Number_of_Time_Steps,
				     const double &Time,
				     const int Block_Number,
				     const int Output_Title,
				     ostream &Out_File);

extern void Output_Tecplot_Quasi3D(Rte2D_Quad_Block &SolnBlk,
				   Rte2D_Input_Parameters &IP,
				   const int Number_of_Time_Steps,
				   const double &Time,
				   const int Block_Number,
				   const int Output_Title,
				   ostream &Out_File);

extern void Output_Exact(Rte2D_Quad_Block &SolnBlk,
			 const int Block_Number,
			 const int Output_Title,
			 ostream &Out_File,
			 double &l1_norm,
			 double &l2_norm,
			 double &max_norm,
			 Rte2D_Input_Parameters &IP );


extern void ICs(Rte2D_Quad_Block &SolnBlk,
		Rte2D_Input_Parameters &IP,
                Rte2D_State *Wo);

extern void BCs(Rte2D_Quad_Block &SolnBlk,
		Rte2D_Input_Parameters &IP);

extern void BCs_Space_March(Rte2D_Quad_Block &SolnBlk, 
			    Rte2D_Input_Parameters &IP);

extern void Prescribe_NonSol(Rte2D_Quad_Block &SolnBlk,
			     Rte2D_Input_Parameters &Input_Parameters);

extern double CFL(Rte2D_Quad_Block &SolnBlk,
                  Rte2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(Rte2D_Quad_Block &SolnBlk, 
                                const double &Dt_min);

extern double L1_Norm_Residual(Rte2D_Quad_Block &SolnBlk, const int &norm);

extern double L2_Norm_Residual(Rte2D_Quad_Block &SolnBlk, const int &norm);

extern double Max_Norm_Residual(Rte2D_Quad_Block &SolnBlk, const int &norm);

extern void Linear_Reconstruction_GreenGauss(Rte2D_Quad_Block &SolnBlk,
                                             const int i,
                                             const int j,
					     const int Limiter);

extern void Linear_Reconstruction_GreenGauss(Rte2D_Quad_Block &SolnBlk,
					     const int Limiter);

extern void Linear_Reconstruction_LeastSquares(Rte2D_Quad_Block &SolnBlk,
                                               const int i,
                                               const int j,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares_2(Rte2D_Quad_Block &SolnBlk,
                                                 const int i,
                                                 const int j,
					         const int Limiter);

extern void Linear_Reconstruction_LeastSquares(Rte2D_Quad_Block &SolnBlk,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares_Angular(Rte2D_Quad_Block &SolnBlk,
						       const int Limiter);

extern void Linear_Reconstruction_LeastSquares_Angular(Rte2D_Quad_Block &SolnBlk,
						       const int i, 
						       const int j,
						       const int Limiter);


extern void Linear_Reconstruction_GreenGauss_Angular(Rte2D_Quad_Block &SolnBlk,
						     const int Limiter);

extern void Linear_Reconstruction_GreenGauss_Angular(Rte2D_Quad_Block &SolnBlk,
						     const int i,
						     const int j,
						     const int Limiter);


extern void Residual_Smoothing(Rte2D_Quad_Block &SolnBlk,
                               const int k_residual,
			       double &epsilon, 
                               const int number_of_Gauss_Seidel_iterations);

extern void Calculate_Refinement_Criteria(double *refinement_criteria,
					  Rte2D_Input_Parameters &IP,
                                          int &number_refinement_criteria,
                                          Rte2D_Quad_Block &SolnBlk);

extern void Fix_Refined_Block_Boundaries(Rte2D_Quad_Block SolnBlk,
                                         const int Fix_North_Boundary,
                                         const int Fix_South_Boundary,
                                         const int Fix_East_Boundary,
                                         const int Fix_West_Boundary);

extern void Unfix_Refined_Block_Boundaries(Rte2D_Quad_Block SolnBlk);

extern void Apply_Boundary_Flux_Corrections(Rte2D_Quad_Block SolnBlk,
                                            const int Number_Neighbours_North_Boundary,
                                            const int Number_Neighbours_South_Boundary,
                                            const int Number_Neighbours_East_Boundary,
                                            const int Number_Neighbours_West_Boundary);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Rte2D_Quad_Block SolnBlk,
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

extern int dUdt_Residual_Evaluation(Rte2D_Quad_Block &SolnBlk,
				    Rte2D_Input_Parameters &Input_Parameters);

extern int dUdt_Multistage_Explicit(Rte2D_Quad_Block &SolnBlk,
   	                            const int i_stage,
                                    Rte2D_Input_Parameters &Input_Parameters);

extern int dUdt_Space_March_Flux_Eval(Rte2D_Quad_Block &SolnBlk,
				      Rte2D_Input_Parameters &Input_Parameters,
				      const int i, const int j,
				      const int v, const int m, const int l,
				      double* &Ix_f, double &Iy_f,
				      const int xd, const int yd);

extern int dUdt_Space_March(Rte2D_Quad_Block &SolnBlk,
			    Rte2D_Input_Parameters &Input_Parameters);

extern int Update_Solution_Multistage_Explicit(Rte2D_Quad_Block &SolnBlk,
   	                                       const int i_stage,
                                               Rte2D_Input_Parameters &Input_Parameters);

/**************************************************************************
 * Rte2D_Quad_Block -- Multiple Block External Subroutines.               *
 **************************************************************************/

extern Rte2D_Quad_Block* Allocate(Rte2D_Quad_Block *Soln_ptr,
                                    Rte2D_Input_Parameters &Input_Parameters);

extern Rte2D_Quad_Block* Deallocate(Rte2D_Quad_Block *Soln_ptr,
                                      Rte2D_Input_Parameters &Input_Parameters);

extern Rte2D_Quad_Block* CreateInitialSolutionBlocks(Grid2D_Quad_Block **InitMeshBlk,
						       Rte2D_Quad_Block *Soln_ptr,
						       Rte2D_Input_Parameters &Input_Parameters,
						       QuadTreeBlock_DataStructure &QuadTree,
						       AdaptiveBlockResourceList &GlobalSolnBlockList,
						       AdaptiveBlock2D_List &LocalSolnBlockList);

extern void ICs(Rte2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
                Rte2D_Input_Parameters &Input_Parameters);

extern int Read_Restart_Solution(Rte2D_Quad_Block *Soln_ptr,
                                 AdaptiveBlock2D_List &Soln_Block_List,
                                 Rte2D_Input_Parameters &Input_Parameters,
		                 int &Number_of_Time_Steps,
                                 double &Time,
                                 CPUTime &CPU_Time);

extern int Write_Restart_Solution(Rte2D_Quad_Block *Soln_ptr,
                                  AdaptiveBlock2D_List &Soln_Block_List,
                                  Rte2D_Input_Parameters &Input_Parameters,
		                  const int Number_of_Time_Steps,
                                  const double &Time,
                                  const CPUTime &CPU_Time);

extern int Output_Tecplot(Rte2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          Rte2D_Input_Parameters &Input_Parameters,
		          const int Number_of_Time_Steps,
                          const double &Time);


extern int Output_Cells_Tecplot(Rte2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                Rte2D_Input_Parameters &Input_Parameters,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Nodes_Tecplot(Rte2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                Rte2D_Input_Parameters &Input_Parameters,
		                const int Number_of_Time_Steps,
                                const double &Time);

extern int Output_Gradients_Tecplot(Rte2D_Quad_Block *Soln_ptr,
				    AdaptiveBlock2D_List &Soln_Block_List,
				    Rte2D_Input_Parameters &Input_Parameters,
				    const int Number_of_Time_Steps,
				    const double &Time);

extern int Output_Tecplot_Quasi3D(Rte2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  Rte2D_Input_Parameters &Input_Parameters,
				  const int Number_of_Time_Steps,
				  const double &Time);

extern int Output_Mesh_Tecplot(Rte2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Rte2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Mesh_Gnuplot(Rte2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Rte2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Exact(Rte2D_Quad_Block *Soln_ptr,
			AdaptiveBlock2D_List &Soln_Block_List,
			Rte2D_Input_Parameters &IP);

extern void Prescribe_NonSol(Rte2D_Quad_Block *Soln_ptr,
			     AdaptiveBlock2D_List &Soln_Block_List,
			     Rte2D_Input_Parameters &Input_Parameters);

extern void BCs(Rte2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
		Rte2D_Input_Parameters &IP);

extern void BCs_Space_March(Rte2D_Quad_Block *Soln_ptr,
			    AdaptiveBlock2D_List &Soln_Block_List,
			    Rte2D_Input_Parameters &IP);


extern double CFL(Rte2D_Quad_Block *Soln_ptr,
                  AdaptiveBlock2D_List &Soln_Block_List,
                  Rte2D_Input_Parameters &Input_Parameters);

extern void Set_Global_TimeStep(Rte2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List, 
                                const double &Dt_min);

extern double L1_Norm_Residual(Rte2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double L2_Norm_Residual(Rte2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List);

extern double Max_Norm_Residual(Rte2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List);

extern void L1_Norm_Residual(Rte2D_Quad_Block *Soln_ptr, 
			     AdaptiveBlock2D_List &Soln_Block_List,
			     double *l1_norm);

extern void L2_Norm_Residual(Rte2D_Quad_Block *Soln_ptr,
			     AdaptiveBlock2D_List &Soln_Block_List,
			     double *l2_norm);

extern void Max_Norm_Residual(Rte2D_Quad_Block *Soln_ptr,
			      AdaptiveBlock2D_List &Soln_Block_List,
			      double *max_norm);

extern void Evaluate_Limiters(Rte2D_Quad_Block *Soln_ptr,
                              AdaptiveBlock2D_List &Soln_Block_List);

extern void Freeze_Limiters(Rte2D_Quad_Block *Soln_ptr,
                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Residual_Smoothing(Rte2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               Rte2D_Input_Parameters &Input_Parameters,
   	                       const int I_Stage);

extern void Apply_Boundary_Flux_Corrections(Rte2D_Quad_Block *Soln_ptr,
                                            AdaptiveBlock2D_List &Soln_Block_List);

extern void Apply_Boundary_Flux_Corrections_Multistage_Explicit(Rte2D_Quad_Block *Soln_ptr,
                                                                AdaptiveBlock2D_List &Soln_Block_List,
                                                                Rte2D_Input_Parameters &Input_Parameters,
   	                                                        const int I_Stage);

extern int dUdt_Multistage_Explicit(Rte2D_Quad_Block *Soln_ptr,
				    AdaptiveBlockResourceList &Global_Soln_Block_List,
                                    AdaptiveBlock2D_List &Soln_Block_List,
                                    Rte2D_Input_Parameters &Input_Parameters,
   	                            const int I_Stage);

extern int Update_Solution_Multistage_Explicit(Rte2D_Quad_Block *Soln_ptr,
                                               AdaptiveBlock2D_List &Soln_Block_List,
                                               Rte2D_Input_Parameters &Input_Parameters,
   	                                       const int I_Stage);

extern int dUdt_Space_March(Rte2D_Quad_Block *Soln_ptr,
			    AdaptiveBlockResourceList &Global_Soln_Block_List,
			    AdaptiveBlock2D_List &Soln_Block_List,
			    Rte2D_Input_Parameters &Input_Parameters);

extern void ScaleGridTo3D( Rte2D_Quad_Block *Soln_ptr,
			   AdaptiveBlock2D_List &Soln_Block_List );

/**************************************************************************
 * Rte2D_Quad_Block -- Multiple Block External Subroutines for Mesh.      *
 **************************************************************************/

extern Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                            Rte2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Broadcast_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                      Rte2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                             Rte2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                                            Rte2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                  Rte2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                 Rte2D_Input_Parameters &Input_Parameters);

extern int Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                          Rte2D_Input_Parameters &Input_Parameters);

extern int Output_Nodes_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                Rte2D_Input_Parameters &Input_Parameters);

extern int Output_Cells_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                Rte2D_Input_Parameters &Input_Parameters);

/**************************************************************************
 * Rte2D_Quad_Block -- Solvers.                                           *
 **************************************************************************/

extern int Rte2DQuadSolver(char *Input_File_Name_ptr,
			   int batch_flag);




#endif /* _RTE2D_QUAD_INCLUDED  */
