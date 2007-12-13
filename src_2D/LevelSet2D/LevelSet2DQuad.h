/******************************************************************//**
 * \file LevelSet2DQuad.h:
 *
 * Header file defining 2D Level Set quadrilateral mesh solution classes.
 **********************************************************************/

#ifndef _LEVELSET2D_QUAD_INCLUDED
#define _LEVELSET2D_QUAD_INCLUDED

// Include 2D LevelSet state, 2D cell, 2D quadrilateral multiblock 
// grid, quadtree, and 2D LevelSet input header files.

#ifndef _LEVELSET2D_STATE_INCLUDED
#include "LevelSet2DState.h"
#endif // _LEVELSET2D_STATE_INCLUDED

#ifndef _INTERFACE2D_INCLUDED
#include "../Interface2D/Interface2D.h"
#endif // _INTERFACE2D_INCLUDED

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

#ifndef _TENSOR2D_INCLUDED
#include "../Math//Tensor2D.h"
#endif // _TENSOR2D_INCLUDED

#ifndef _LEVELSET2D_INPUT_INCLUDED
#include "LevelSet2DInput.h"
#endif // _LEVELSET2D_INPUT_INCLUDED

#define	NUMBER_OF_RESIDUAL_VECTORS   3

#define CLOSEDCURVE     50
#define DOMAINBOUNDED   51

#define CLEAN           0  //!< Clean flag used in interface tracing.
#define INFECTED        1  //!< Infected flag used in interface tracing.
#define UNCHECKED       2  //!< Probably a deprecated flag.

//! Interface retrieval debugging option. Uncomment to print debugging files.
// #define _RETRIEVE_DEBUG_

/*!
 * Class: CutType
 *
 * @brief A class defined to contain face status values to faciltate 
 * the interface (zero contour) tracing algorithm.
 *
 * A struct defined to contain the face status values: INFECTED or
 * CLEAN.  The north face includes the NW corner, the south face
 * includes the SE corner, the east face includes the NE corner, and
 * the west face includes the SW corner.  Used for retrieving/capturing
 * the zero level set contour.
 */
class CutType {
 public:
  int     center, //!< Cell centre flag.
           north, //!< North face cut flag.
           south, //!< South face cut flag.
            east, //!< East face cut flag.
            west; //!< West face cut flag.
  int north_east, //!< North-east face cut flag.
      south_west; //!< South-west face cut flag.

  //! Default constructor.
  CutType(void) { clean(); }

  //! Reset/clean cut types.
  void clean(void) {
    center = CLEAN;
    north = CLEAN; south = CLEAN;
    east = CLEAN; west = CLEAN;
    north_east = CLEAN; south_west = CLEAN;
  }

};

/*!
 * Class: Trace_Data
 */
class Trace_Data {
 public:
  int    i, //!< Starting points.
         j, //!< Starting points.
      face; //!< Starting face.

  Trace_Data(void) { Zero(); }

  Trace_Data(int ii, int jj, int fface) {
    i = ii; j = jj; face = fface;
  }

  void Zero(void) { i = 0; j = 0; face = 0; }

  Trace_Data& operator =(const Trace_Data &T) {
    i = T.i; j = T.j; face = T.face;
  }

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const Trace_Data &T) {
    out_file << " " << T.i << " " << T.j << " " << T.face;
  }
  friend istream &operator >> (istream &in_file, Trace_Data &T) {
    in_file.setf(ios::skipws);
    in_file >> T.i >> T.j >> T.face;
    in_file.unsetf(ios::skipws);
  }
  //@}

};

/*!
 * Class: LevelSet2D_Quad_Block
 *
 * @brief Class definition of a 2D level set solution block.
 *
 * This class defines a single computational block for the 2D level set
 * method.  The level set method requires the solution of three scalar
 * equations: a Hamilton-Jacobi type equation, the Eikonal equations,
 * and the scalar extension equations.  The state variable, U, includes
 * the level set function, psi, the scalar front speed, F, and a bulk
 * flow=field, V = (u,v).  
 *
 * \verbatim
 * Member functions
 *          U   -- Return level set variable solution for the block.
 *       Grid   -- Return the solution block quadrilateral grid or mesh.
 *         dt   -- Return local time step for the solution block.
 *       dUdt   -- Return the solution block residuals.
 *       dUdx   -- Return the unlimited variable solution gradients
 *                 (x-direction) for the block.
 *      dUdxp   -- Return the unlimited variable solution gradients
 *                 (x-direction, right-biased stencil) for the block.
 *      dUdxm   -- Return the unlimited variable solution gradients
 *                 (x-direction, left-biased stencil) for the block.
 *     ddUdxx   -- Return the second-derivative of the variable solution
 *                 (x-direction) for the block.
 *       dUdy   -- Return the unlimited variable solution gradients 
 *                 (y-direction) for the block.
 *      dUdyp   -- Return the unlimited variable solution gradients
 *                 (y-direction, right-biased stencil) for the block.
 *      dUdym   -- Return the unlimited variable solution gradients
 *                 (y-direction, left-biased stencil) for the block.
 *     ddUdyy   -- Return the second-derivative of the variable solution
 *                 (y-direction) for the block.
 *      kappa   -- Return the curvature of the variable solution for the block.
 *    gradMag   -- Return the centered magnitude of the gradient of the
 *                 variable solution for the block.
 *        phi   -- Return the solution slope limiters.
 *         Uo   -- Return initial solution state.
 *        Uoo   -- Return initial solution state (before Eikonal iterations).
 *        NCi   -- Return number of cells in the i-direction.
 *        ICl   -- Return lower index for cells in the i-direction.
 *        ICu   -- Return upper index for cells in the i-direction.
 *        NCj   -- Return number of cells in the j-direction.
 *        JCl   -- Return lower index for cells in the j-direction.
 *        JCu   -- Return upper index for cells in the j-direction.
 *     Nghost   -- Number of ghost cells.
 * Number_of_Interfaces -- Number of arbitrary interfaces.
 *  Interface   -- Array of interfaces.
 *        cut   -- Cell cut types.
 *   allocate   -- Allocate memory for structured quadrilateral
 *                 solution block.
 * deallocate   -- Deallocate memory for structured quadrilateral
 *                 solution block.
 *         Un   -- Return level set solution state at the specified node.
 *       UnNW   -- Return level set solution state at the north-west node.
 *       UnNE   -- Return level set solution state at the north-east node.
 *       UnSW   -- Return level set solution state at the south-west node.
 *       UnSE   -- Return level set solution state at the south-east node.
 *
 * Member functions required for message passing.
 *   NumVar                    -- Returns number of solution variables
 *                                in the solution state vectors.
 *   LoadSendBuffer            -- Loads send buffer.
 *   LoadSendBuffer_F2C        -- Loads send buffer for fine to coarse  
 *                                block messages.
 *   UnloadReceiveBuffer       -- Unloads receive buffer.
 *   UnloadReceiveBuffer_F2C   -- Unloads receive buffer for fine to
 *                                coarse block messages.
 *   SubcellReconstruction     -- Performs subcell solution reconstruction 
 *                                used in adaptive mesh refinement.
 *   LoadSendBuffer_Flux_F2C   -- Loads send buffer for sending
 *                                conservative flux corrections from fine
 *                                to coarse solution blocks.
 *   UnLoadSendBuffer_Flux_F2C -- Loads send buffer for sending
 *                                conservative flux corrections from
 *                                fine to coarse solution blocks.
 *
 * Member operators
 *      S -- a 2D level set solution block.
 *
 * S = S;
 * cout << S; (output function)
 * cin  >> S; (input function)
 * \endverbatim
 */
class LevelSet2D_Quad_Block {
private:
public:
  //@{ @name Solution state arrays:
  LevelSet2DState            **U; //!< Level set solution state.
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
  LevelSet2DState         ***dUdt; //!< Solution residual.
  LevelSet2DState            **Uo, //!< Initial solution state.
                            **Uoo; //!< Initial solution state (before Eikonal iterations).
  double                   **sign; //!< Sign of distance function.
  static int    residual_variable; //!< Static integer that indicates which variable is used for residual calculations.
  //@}

  //@{ @name Solution gradient arrays:
  LevelSet2DState          **dUdx, //!< Unlimited solution gradient (x-direction).
                           **dUdy; //!< Unlimited solution gradient (y-direction).
  LevelSet2DState         **dUdxp, //!< Unlimited solution gradient (x-direction, right-biased stencil).
                          **dUdxm, //!< Unlimited solution gradient (x-direction, left-biased stencil).
                          **dUdyp, //!< Unlimited solution gradient (y-direction, right-biased stencil).
                          **dUdym; //!< Unlimited solution gradient (y-direction, left-biased stencil).
  LevelSet2DState        **ddUdxx, //!< Second-derivative of the solution (x-direction).
                         **ddUdyy, //!< Second-derivative of the solution (y-direction).
                          **kappa, //!< Curvature of the solution.
                        **gradMag; //!< Magnitude of the centered gradient of the solution.
  LevelSet2DState           **phi; //!< Solution slope limiter.
  //@}

  //@{ @name Boundary solution flux arrays:
  LevelSet2DState          *FluxN, //!< North boundary solution flux.
                           *FluxS, //!< South boundary solution flux.
                           *FluxE, //!< East boundary solution flux.
                           *FluxW; //!< West boundary solution flux.
  //@}

  //@{ @name Interface storage and retracing information:
  Interface2D_List Interface_List; //!< Interface list.
  CutType                   **cut; //!< Cell cut types.
  LinkedList<Trace_Data>    Trace; //!< Data to assist in interface tracing.
  //@}

  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  LevelSet2D_Quad_Block(void) {
    // Grid size and variables.
    NCi = 0; ICl = 0; ICu = 0;
    NCj = 0; JCl = 0; JCu = 0;
    Nghost = 0;
    // Solution variables.
    U      = NULL; Uo     = NULL; Uoo = NULL;
    dUdt   = NULL; dt     = NULL;
    dUdx   = NULL; dUdy   = NULL;
    dUdxp  = NULL; dUdyp  = NULL; dUdxm = NULL; dUdym = NULL;
    ddUdxx = NULL; ddUdyy = NULL; kappa = NULL; gradMag = NULL;
    phi    = NULL;
    sign   = NULL;
    FluxN  = NULL; FluxS = NULL; FluxE = NULL; FluxW = NULL;
    // Interfaces variables.
    //Interface_List = NULL;
    // Cut class.
    cut = NULL;
  }

  //! Copy constructor.
  LevelSet2D_Quad_Block(const LevelSet2D_Quad_Block &Soln) {
    // Grid size and variables.
    NCi = Soln.NCi; ICl = Soln.ICl; ICu = Soln.ICu; 
    NCj = Soln.NCj; JCl = Soln.JCl; JCu = Soln.JCu;
    Nghost = Soln.Nghost;
    Grid = Soln.Grid;
    // Solution variables.
    U     = Soln.U;       Uo    = Soln.Uo;    Uoo   = Soln.Uoo;
    dUdt  = Soln.dUdt;    dt    = Soln.dt;
    dUdx  = Soln.dUdx;    dUdy  = Soln.dUdy;
    dUdxp = Soln.dUdxp;   dUdyp = Soln.dUdyp;
    dUdxm = Soln.dUdxm;   dUdym = Soln.dUdym;
    ddUdxx = Soln.ddUdxx; ddUdyy = Soln.ddUdyy; kappa = Soln.kappa; gradMag = Soln.gradMag;
    phi   = Soln.phi;
    sign  = Soln.sign;
    FluxN = Soln.FluxN; FluxS = Soln.FluxS; FluxE = Soln.FluxE; FluxW = Soln.FluxW;
    // Interfaces variables.
    Interface_List = Soln.Interface_List;
    // Cut class.
    cut = Soln.cut;
    Trace = Soln.Trace;
  }

  // Destructor.
  // ~LevelSet2D_Quad_Block(void);
  // Use automatically generated destructor.
  //@}

  // Assignment operator.
  // LevelSet2D_Quad_Block operator = (const LevelSet2D_Quad_Block &Soln);
  // Use automatically generated assignment operator.

  //@{ @name Allocation and deallocation function.
  //! Allocate memory for structured quadrilateral solution block.
  void allocate(const int Ni, const int Nj, const int Ng);
  //! Deallocate memory for structured quadrilateral solution block.
  void deallocate(void);
  //@}

  //@{ @name Bilinear interplation (Holmes and Connell).
  //! Return primitive solution state at specified node.
  LevelSet2DState Un(const int ii, const int jj);

  //! Return solution state at cell nodes.
  LevelSet2DState UnNW(const int ii, const int jj); //!< Return level set solution state at cell NW node.
  LevelSet2DState UnNE(const int ii, const int jj); //!< Return level set solution state at cell NE node.
  LevelSet2DState UnSE(const int ii, const int jj); //!< Return level set solution state at cell SE node.
  LevelSet2DState UnSW(const int ii, const int jj); //!< Return level set solution state at cell SW node.
  //@}

  //@{ @name Bilinear interplation (Zingg and Yarrow).
  //! Return the level set solution state at the specified node.
  LevelSet2DState Uzy(const int ii, const int jj);

  //! Return the level set solution state at cell nodes.
  LevelSet2DState UzyNW(const int ii, const int jj); //!< Return level set solution state at cell NW node.
  LevelSet2DState UzyNE(const int ii, const int jj); //!< Return level set solution state at cell NE node.
  LevelSet2DState UzySE(const int ii, const int jj); //!< Return level set solution state at cell SE node.
  LevelSet2DState UzySW(const int ii, const int jj); //!< Return level set solution state at cell SW node.
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const LevelSet2D_Quad_Block &Soln);
  friend istream &operator >> (istream &in_file, LevelSet2D_Quad_Block &Soln);
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
  //! Load the fine grid solution message passing buffer.
  int LoadSendBuffer_FineGridSolution(double *buffer,
				      int &buffer_count,
				      const int buffer_size,
				      const int i_min,
				      const int i_max,
				      const int i_inc,
				      const int j_min,
				      const int j_max,
				      const int j_inc) { return 0; }
  //! Unload the fine grid solution message passing buffer.
  int UnloadReceiveBuffer_FineGridSolution(double *buffer,
					   int &buffer_count,
					   const int buffer_size,
					   const int i_min,
					   const int i_max,
					   const int i_inc,
					   const int j_min,
					   const int j_max,
					   const int j_inc) { return 0; }
  //! Subcell solution reconstruction within given computational cell.
  void SubcellReconstruction(const int i,
                             const int j,
                             const int Limiter);
  //! Load and unload conservative flux message passing buffer.
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
 * LevelSet2D_Quad_Block::allocate -- Allocate memory.                *
 **********************************************************************/
inline void LevelSet2D_Quad_Block::allocate(const int Ni, const int Nj, const int Ng) {
  assert(Ni >= 2*Ng && Nj > 2*Ng && Ng > 1); 
  Grid.allocate(Ni, Nj, Ng);
  NCi = Ni+2*Ng; ICl = Ng; ICu = Ni+Ng-1; 
  NCj = Nj+2*Ng; JCl = Ng; JCu = Nj+Ng-1; Nghost = Ng;
  U  = new LevelSet2DState*[NCi]; 
  dt = new double*[NCi]; 
  dUdt = new LevelSet2DState**[NCi];
  dUdx = new LevelSet2DState*[NCi];
  dUdy = new LevelSet2DState*[NCi];
  dUdxm = new LevelSet2DState*[NCi];
  dUdxp = new LevelSet2DState*[NCi];
  dUdym = new LevelSet2DState*[NCi];
  dUdyp = new LevelSet2DState*[NCi];
  ddUdxx = new LevelSet2DState*[NCi];
  ddUdyy = new LevelSet2DState*[NCi];
  kappa = new LevelSet2DState*[NCi];
  gradMag = new LevelSet2DState*[NCi];
  phi  = new LevelSet2DState*[NCi];
  sign = new double*[NCi];
  Uo   = new LevelSet2DState*[NCi];
  Uoo  = new LevelSet2DState*[NCi];
  cut  = new CutType*[NCi];
  for (int i = 0; i < NCi; i++) {
    U[i] = new LevelSet2DState[NCj];
    dt[i] = new double[NCj]; 
    dUdt[i] = new LevelSet2DState*[NCj];
    for (int j = 0; j < NCj; j++)
      dUdt[i][j] = new LevelSet2DState[NUMBER_OF_RESIDUAL_VECTORS];
    dUdx[i] = new LevelSet2DState[NCj];
    dUdy[i] = new LevelSet2DState[NCj];
    dUdxm[i] = new LevelSet2DState[NCj];
    dUdxp[i] = new LevelSet2DState[NCj];
    dUdym[i] = new LevelSet2DState[NCj];
    dUdyp[i] = new LevelSet2DState[NCj];
    ddUdxx[i] = new LevelSet2DState[NCj];
    ddUdyy[i] = new LevelSet2DState[NCj];
    kappa[i] = new LevelSet2DState[NCj];
    gradMag[i] = new LevelSet2DState[NCj];
    phi[i]  = new LevelSet2DState[NCj];
    sign[i] = new double[NCj];
    Uo[i]   = new LevelSet2DState[NCj];
    Uoo[i]  = new LevelSet2DState[NCj];
    cut[i]  = new CutType[NCj];
    for (int j = 0; j < NCj; j++) {
      U[i][j]    = LevelSet2D_ZERO;
      dUdx[i][j] = LevelSet2D_ZERO;
      dUdy[i][j] = LevelSet2D_ZERO;
      dUdxm[i][j] = LevelSet2D_ZERO;
      dUdxp[i][j] = LevelSet2D_ZERO;
      dUdym[i][j] = LevelSet2D_ZERO;
      dUdyp[i][j] = LevelSet2D_ZERO;
      ddUdxx[i][j] = LevelSet2D_ZERO;
      ddUdyy[i][j] = LevelSet2D_ZERO;
      kappa[i][j] = LevelSet2D_ZERO;
      gradMag[i][j] = LevelSet2D_ZERO;
      phi[i][j]  = LevelSet2D_ZERO;
      sign[i][j] = ZERO;
      Uo[i][j]   = LevelSet2D_ZERO;
      Uoo[i][j]  = LevelSet2D_ZERO;
      dt[i][j]   = ZERO;
      for (int k = 0; k < NUMBER_OF_RESIDUAL_VECTORS; k++)
	dUdt[i][j][k] = LevelSet2D_ZERO;
    }
  }
  FluxN = new LevelSet2DState[NCi]; FluxS = new LevelSet2DState[NCi];
  FluxE = new LevelSet2DState[NCj]; FluxW = new LevelSet2DState[NCj];
}

/**********************************************************************
 * LevelSet2D_Quad_Block::deallocate -- Deallocate memory.            *
 **********************************************************************/
inline void LevelSet2D_Quad_Block::deallocate(void) {
  Grid.deallocate();
  for (int i = 0; i < NCi; i++) {
    delete []U[i];  U[i]  = NULL;
    delete []dt[i]; dt[i] = NULL; 
    for (int j = 0; j < NCj; j++) {
      delete []dUdt[i][j]; dUdt[i][j] = NULL;
    }
    delete []dUdt[i];  dUdt[i] = NULL;
    delete []dUdx[i];  dUdx[i] = NULL;
    delete []dUdy[i];  dUdy[i] = NULL;
    delete []dUdxm[i]; dUdxm[i] = NULL;
    delete []dUdxp[i]; dUdxp[i] = NULL;
    delete []dUdym[i]; dUdym[i] = NULL;
    delete []dUdyp[i]; dUdyp[i] = NULL;
    delete []ddUdxx[i]; ddUdxx[i] = NULL;
    delete []ddUdyy[i]; ddUdyy[i] = NULL;
    delete []kappa[i]; kappa[i] = NULL;
    delete []gradMag[i]; gradMag[i] = NULL;
    delete []phi[i];   phi[i]  = NULL;
    delete []sign[i];  sign[i] = NULL;
    delete []Uo[i];    Uo[i]   = NULL;
    delete []Uoo[i];   Uoo[i]  = NULL;
    delete []cut[i];   cut[i]  = NULL;
  }
  delete []U;     U    = NULL;
  delete []dt;    dt   = NULL; 
  delete []dUdt;  dUdt = NULL;
  delete []dUdx;  dUdx = NULL;
  delete []dUdy;  dUdy = NULL;
  delete []dUdxm; dUdxm = NULL;
  delete []dUdxp; dUdxp = NULL;
  delete []dUdym; dUdym = NULL;
  delete []dUdyp; dUdyp = NULL;
  delete []ddUdxx; ddUdxx = NULL;
  delete []ddUdyy; ddUdyy = NULL;
  delete []kappa; kappa = NULL;
  delete []gradMag; gradMag = NULL;
  delete []phi;   phi  = NULL; 
  delete []sign;  sign = NULL; 
  delete []Uo;    Uo   = NULL;
  delete []Uoo;   Uoo  = NULL;
  delete []cut;   cut  = NULL;
  //delete []Trace; Trace = NULL;
  delete []FluxN; FluxN = NULL; delete []FluxS; FluxS = NULL;
  delete []FluxE; FluxE = NULL; delete []FluxW; FluxW = NULL;
  NCi = 0; ICl = 0; ICu = 0; NCj = 0; JCl = 0; JCu = 0; Nghost = 0;
  if (Interface_List.Interface != NULL) Interface_List.deallocate();
}

/**********************************************************************
 * LevelSet2D_Quad_Block::Un -- Return the level set solution state   *
 *                              at the specified node using the       *
 *                              bilinear interpolation of Holmes and  *
 *                              Connell (AIAA-1989-1932).             *
 **********************************************************************/
inline LevelSet2DState LevelSet2D_Quad_Block::Un(const int ii,
                                                 const int jj) {
  double w1, w2, w3, w4;
  LevelSet2DState U1, U2, U3, U4;
  Vector2D X0, X1, X2, X3, X4, lambda, R;
  Tensor2D I;
  // Summarize cell-centres and state:
  X1 = Grid.Cell[ii-1][jj-1].Xc;  U1 = U[ii-1][jj-1];
  X2 = Grid.Cell[ii  ][jj-1].Xc;  U2 = U[ii  ][jj-1];
  X3 = Grid.Cell[ii-1][jj  ].Xc;  U3 = U[ii-1][jj  ];
  X4 = Grid.Cell[ii  ][jj  ].Xc;  U4 = U[ii  ][jj  ];
  X0 = Grid.Node[ii][jj].X;
  // Determine weighting coefficients:
  R = X1 + X2 + X3 + X4 - 4.0*X0;
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
  // Return the interpolated state:
  return (w1*U1 + w2*U2 + w3*U3+ w4*U4)/(w1 + w2 + w3 + w4);
}

/**********************************************************************
 * LevelSet2D_Quad_Block::Un?? -- Get cell node solution states.      *
 **********************************************************************/
inline LevelSet2DState LevelSet2D_Quad_Block::UnNW(const int ii, const int jj) {
  return Un(ii  ,jj+1);
}

inline LevelSet2DState LevelSet2D_Quad_Block::UnNE(const int ii, const int jj) {
  return Un(ii+1,jj+1);
}

inline LevelSet2DState LevelSet2D_Quad_Block::UnSE(const int ii, const int jj) {
  return Un(ii+1,jj  );
}

inline LevelSet2DState LevelSet2D_Quad_Block::UnSW(const int ii, const int jj) {
  return Un(ii  ,jj  );
}

/**********************************************************************
 * LevelSet2D_Quad_Block::Uzy -- Return the level set solution state  *
 *                               at the specified node using bilinear *
 *                               interpolation, Zingg and Yarrow      *
 *                               (SIAM J. Sci. Stat. Comput. Vol. 13  *
 *                               No. 3 1992).                         *
 **********************************************************************/
inline LevelSet2DState LevelSet2D_Quad_Block::Uzy(const int ii, const int jj) {
  double ax, bx, cx, dx, ay, by, cy, dy, aa, bb, cc, x, y, 
         eta1, zeta1, eta2, zeta2, eta, zeta;
  LevelSet2DState A, B, C, D;
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
    if (fabs(bb) >= TOLER*TOLER) zeta1 = -cc/bb;
    else zeta1 = -cc/sgn(bb)*(TOLER*TOLER);
    if (fabs(cy+dy*zeta1) >= TOLER*TOLER) eta1 = (y-ay-by*zeta1)/(cy+dy*zeta1);
    else eta1 = HALF;
    zeta2 = zeta1;
    eta2  = eta1;
  } else {
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) zeta1 = HALF*(-bb+sqrt(bb*bb-FOUR*aa*cc))/aa;
    else zeta1 = -HALF*bb/aa;
    if (fabs(cy+dy*zeta1) < TOLER*TOLER) eta1 = -ONE;
    else eta1 = (y-ay-by*zeta1)/(cy+dy*zeta1);
    if (bb*bb-FOUR*aa*cc >= TOLER*TOLER) zeta2 = HALF*(-bb-sqrt(bb*bb-FOUR*aa*cc))/aa;
    else zeta2 = -HALF*bb/aa;
    if (fabs(cy+dy*zeta2) < TOLER*TOLER) eta2 = -ONE;
    else eta2 = (y-ay-by*zeta2)/(cy+dy*zeta2);
  }
  if (zeta1 > -TOLER && zeta1 < ONE + TOLER &&
      eta1  > -TOLER && eta1  < ONE + TOLER) {
    zeta = zeta1; eta = eta1;
  } else if (zeta2 > -TOLER && zeta2 < ONE + TOLER &&
	     eta2  > -TOLER && eta2  < ONE + TOLER) {
    zeta = zeta2; eta = eta2;
  } else {
    zeta = HALF;  eta = HALF;
  }
  A = U[ii-1][jj-1];
  B = U[ii-1][jj  ] - U[ii-1][jj-1];
  C = U[ii  ][jj-1] - U[ii-1][jj-1];
  D = U[ii  ][jj  ] + U[ii-1][jj-1] - U[ii-1][jj  ] - U[ii  ][jj-1];
  return A + B*zeta + C*eta + D*zeta*eta;
}

/**********************************************************************
 * LevelSet2D_Quad_Block::Uzy?? -- Get cell node solution states.     *
 **********************************************************************/
inline LevelSet2DState LevelSet2D_Quad_Block::UzyNW(const int ii, const int jj) {
  return Uzy(ii  ,jj+1);
}

inline LevelSet2DState LevelSet2D_Quad_Block::UzyNE(const int ii, const int jj) {
  return Uzy(ii+1,jj+1);
}

inline LevelSet2DState LevelSet2D_Quad_Block::UzySE(const int ii, const int jj) {
  return Uzy(ii+1,jj  );
}

inline LevelSet2DState LevelSet2D_Quad_Block::UzySW(const int ii, const int jj) {
  return Uzy(ii  ,jj  );
}

/**********************************************************************
 * LevelSet2D_Quad_Block -- Input-output operators.                   *
 **********************************************************************/
inline ostream &operator << (ostream &out_file, const LevelSet2D_Quad_Block &SolnBlk) {
  out_file << SolnBlk.Grid;
  out_file << SolnBlk.NCi << " " << SolnBlk.ICl << " " << SolnBlk.ICu << endl;
  out_file << SolnBlk.NCj << " " << SolnBlk.JCl << " " << SolnBlk.JCu << endl;
  out_file << SolnBlk.Nghost << endl;
  if (SolnBlk.NCi == 0 || SolnBlk.NCj == 0) return out_file;
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      out_file << SolnBlk.U[i][j] << endl;
    }
  }
  return out_file;
}

inline istream &operator >> (istream &in_file, LevelSet2D_Quad_Block &SolnBlk) {
  int ni, il, iu, nj, jl, ju, ng, nsp;
  Grid2D_Quad_Block New_Grid; in_file >> New_Grid; 
  in_file.setf(ios::skipws);
  in_file >> ni >> il >> iu;
  in_file >> nj >> jl >> ju;
  in_file >> ng;
  in_file.unsetf(ios::skipws);
  if (ni == 0 || nj == 0) {
    SolnBlk.deallocate(); return in_file;
  }
  if (SolnBlk.U == NULL || SolnBlk.NCi != ni || SolnBlk.NCj != nj) {
    if (SolnBlk.U != NULL) SolnBlk.deallocate(); 
    SolnBlk.allocate(ni-2*ng,nj-2*ng,ng);
  }
  Copy_Quad_Block(SolnBlk.Grid,New_Grid); New_Grid.deallocate();
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      in_file >> SolnBlk.U[i][j];
      for (int k = 0; k < NUMBER_OF_RESIDUAL_VECTORS; k++)
	SolnBlk.dUdt[i][j][k] = LevelSet2D_ZERO;
      SolnBlk.dUdx[i][j] = LevelSet2D_ZERO;
      SolnBlk.dUdy[i][j] = LevelSet2D_ZERO;
      SolnBlk.dUdxp[i][j] = LevelSet2D_ZERO;
      SolnBlk.dUdyp[i][j] = LevelSet2D_ZERO;
      SolnBlk.dUdxm[i][j] = LevelSet2D_ZERO;
      SolnBlk.dUdym[i][j] = LevelSet2D_ZERO;
      SolnBlk.phi[i][j]  = LevelSet2D_ZERO;
      SolnBlk.Uo[i][j]   = LevelSet2D_ZERO;
      SolnBlk.dt[i][j]   = ZERO;
    }
  }
  return in_file;
}

/**********************************************************************
 *                                                                    *
 * MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING.                     *
 *                                                                    *
 **********************************************************************/

/**********************************************************************
 * LevelSet2D_Quad_Block::NumVar -- Return number of state variables. *
 **********************************************************************/
inline int LevelSet2D_Quad_Block::NumVar(void) {
  return int(NUM_VAR_LEVELSET2D);
}

/**********************************************************************
 * LevelSet2D_Quad_Block::LoadSendBuffer -- Load send message buffer. *
 **********************************************************************/
inline int LevelSet2D_Quad_Block::LoadSendBuffer(double *buffer,
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
      for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = U[i][j][k];
      }
    }
  }

  return 0;

}

/**********************************************************************
 * LevelSet2D_Quad_Block::LoadSendBuffer_F2C -- Loads send message    *
 *                                              buffer for fine to    *
 *                                              coarse block message  *
 *                                              passing.              *
 **********************************************************************/
inline int LevelSet2D_Quad_Block::LoadSendBuffer_F2C(double *buffer,
                                                     int &buffer_count,
                                                     const int buffer_size,
                                                     const int i_min, 
                                                     const int i_max,
                                                     const int i_inc,
                                                     const int j_min, 
                                                     const int j_max,
                                                     const int j_inc) {

  for (int j = j_min; ((j_inc+2)/4) ? (j < j_max):(j > j_max); j += j_inc) {
    for (int i = i_min;  ((i_inc+2)/4) ? (i < i_max):(i > i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = UnNE(i,j)[k];
      }
    }
  }

  return 0;

}

/**********************************************************************
 * LevelSet2D_Quad_Block::LoadSendBuffer_C2F -- Loads send message    *
 *                                              buffer for coarse to  *
 *                                              fine block message    *
 *                                              passing.              *
 **********************************************************************/
inline int LevelSet2D_Quad_Block::LoadSendBuffer_C2F(double *buffer,
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
  LevelSet2DState Ufine;
  int Limiter = LIMITER_ONE;

  if (j_inc > 0) {
    if (i_inc > 0) {
      for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
	  // Perform limited linear least squares reconstruction in cell (i,j).
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
	    Ufine = U[i][j] + (phi[i][j]^dUdx[i][j])*dX.x +
	                      (phi[i][j]^dUdy[i][j])*dX.y;
	    for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
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
	    Ufine = U[i][j] + (phi[i][j]^dUdx[i][j])*dX.x +
	                      (phi[i][j]^dUdy[i][j])*dX.y;
	    for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
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
	    Ufine = U[i][j] + (phi[i][j]^dUdx[i][j])*dX.x +
                              (phi[i][j]^dUdy[i][j])*dX.y;
	    for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) { 
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
	    Ufine = U[i][j] + (phi[i][j]^dUdx[i][j])*dX.x +
	                      (phi[i][j]^dUdy[i][j])*dX.y;
	    for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
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
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
               Grid.Node[i+1][j_min].X+
               Grid.Cell[i][j_min].Xc+
               HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
               Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                                (phi[i][j_min]^dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
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
	  Ufine = U[i][j_min] + (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                                (phi[i][j_min]^dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
	  // Evaluate NE sub (fine) cell values.
	  dX = (Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
		Grid.Node[i+1][j_min+1].X)/FOUR -
	        Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                                (phi[i][j_min]^dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
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
	  Ufine = U[i][j_min] + (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                                (phi[i][j_min]^dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) { 
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
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i][j_min+1].X)+
		Grid.Cell[i][j_min].Xc+
		Grid.Node[i][j_min+1].X+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X))/FOUR -
                Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                                (phi[i][j_min]^dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) { 
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
	  Ufine = U[i][j_min] + (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                                (phi[i][j_min]^dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
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
	  Ufine = U[i][j_min] + (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                                (phi[i][j_min]^dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
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
	  Ufine = U[i][j_min] + (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                                (phi[i][j_min]^dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
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
	  dX = (Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X)+
		HALF*(Grid.Node[i][j_min+1].X+Grid.Node[i+1][j_min+1].X)+
		Grid.Node[i+1][j_min+1].X)/FOUR -
                Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                                (phi[i][j_min]^dUdy[i][j_min])*dX.y;
	  for (int k = 1 ; k <= NUM_VAR_LEVELSET2D; ++ k) {
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
	  Ufine = U[i][j_min] + (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                                (phi[i][j_min]^dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) { 
	    buffer_count++;
	    if (buffer_count >= buffer_size) return 1;
	    buffer[buffer_count] = Ufine[k];
	  }
	}
	for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
	  // Evaluate SE sub (fine) cell values.
	  dX = (HALF*(Grid.Node[i][j_min].X+Grid.Node[i+1][j_min].X)+
		Grid.Node[i+1][j_min].X+
		Grid.Cell[i][j_min].Xc+
		HALF*(Grid.Node[i+1][j_min].X+Grid.Node[i+1][j_min+1].X))/FOUR -
	        Grid.Cell[i][j_min].Xc;
	  Ufine = U[i][j_min] + (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                                (phi[i][j_min]^dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
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
	  Ufine = U[i][j_min] + (phi[i][j_min]^dUdx[i][j_min])*dX.x +
                                (phi[i][j_min]^dUdy[i][j_min])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
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
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min+1][j].X)+
		Grid.Node[i_min+1][j].X+
		Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                                (phi[i_min][j]^dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
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
	  Ufine = U[i_min][j] + (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                                (phi[i_min][j]^dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
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
	  Ufine = U[i_min][j] + (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                                (phi[i_min][j]^dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
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
	  Ufine = U[i_min][j] + (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                                (phi[i_min][j]^dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) { 
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
	  dX = (HALF*(Grid.Node[i_min][j].X+Grid.Node[i_min][j+1].X)+
		Grid.Cell[i_min][j].Xc+
		Grid.Node[i_min][j+1].X+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X))/FOUR -
                Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                                (phi[i_min][j]^dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) { 
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
	  Ufine = U[i_min][j] +
	          (phi[i_min][j]^dUdx[i_min][j])*dX.x +
	          (phi[i_min][j]^dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
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
	  Ufine = U[i_min][j] + (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                                (phi[i_min][j]^dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
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
	  Ufine = U[i_min][j] + (phi[i_min][j]^dUdx[i_min][j])*dX.x +
                                (phi[i_min][j]^dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
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
	  dX = (Grid.Cell[i_min][j].Xc+
		HALF*(Grid.Node[i_min+1][j].X+Grid.Node[i_min+1][j+1].X)+
		HALF*(Grid.Node[i_min][j+1].X+Grid.Node[i_min+1][j+1].X)+
		Grid.Node[i_min+1][j+1].X)/FOUR -
                Grid.Cell[i_min][j].Xc;
	  Ufine = U[i_min][j] + (phi[i_min][j]^dUdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
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
	  Ufine = U[i_min][j] + (phi[i_min][j]^dUdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) { 
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
	  Ufine = U[i_min][j] + (phi[i_min][j]^dUdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
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
	  Ufine = U[i_min][j] + (phi[i_min][j]^dUdx[i_min][j])*dX.x +
	                        (phi[i_min][j]^dUdy[i_min][j])*dX.y;
	  for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
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
 * LevelSet2D_Quad_Block::UnloadReceiveBuffer -- Unloads receive      *
                                                 message buffer.      *
 **********************************************************************/
inline int LevelSet2D_Quad_Block::UnloadReceiveBuffer(double *buffer,
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
      buffer_count += NUM_VAR_LEVELSET2D;
      if (buffer_count >= buffer_size) return 1;
      U[i][j] = LevelSet2DState(buffer[buffer_count-3],
				buffer[buffer_count-2],
				buffer[buffer_count-1],
				buffer[buffer_count  ]);
    }
  }
  return 0;

}

/**********************************************************************
 * LevelSet2D_Quad_Block::UnloadReceiveBuffer_F2C -- Unloads receive  *
 *                                       message buffer for fine to   *
 *                                       coarse block message passing.*
 **********************************************************************/
inline int LevelSet2D_Quad_Block::UnloadReceiveBuffer_F2C(double *buffer,
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
      buffer_count = buffer_count + NUM_VAR_LEVELSET2D;
      if (buffer_count >= buffer_size) return 1;
      U[i][j] = LevelSet2DState(buffer[buffer_count-3],
				buffer[buffer_count-2],
				buffer[buffer_count-1],
				buffer[buffer_count  ]);
    }
  }
  return 0;

}

/**********************************************************************
 * LevelSet2D_Quad_Block::UnloadReceiveBuffer_C2F -- Unloads receive  *
 *                                       message buffer for coarse to *
 *                                       fine block message passing.  *
 **********************************************************************/
inline int LevelSet2D_Quad_Block::UnloadReceiveBuffer_C2F(double *buffer,
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
      buffer_count = buffer_count + NUM_VAR_LEVELSET2D;
      if (buffer_count >= buffer_size) return 1;
      U[i][j] = LevelSet2DState(buffer[buffer_count-3],
				buffer[buffer_count-2],
				buffer[buffer_count-1],
				buffer[buffer_count  ]);
    }
  }
  return 0;

}

/**********************************************************************
 * LevelSet2D_Quad_Block::SubcellReconstruction --                    *
 *               Performs the subcell reconstruction of solution      *
 *               state within a given cell (i,j) of the computational *
 *               mesh for the specified quadrilateral solution block. *
 **********************************************************************/
inline void LevelSet2D_Quad_Block::SubcellReconstruction(const int i,
                                                         const int j,
                                                         const int Limiter) {

  int n, n2, n_pts, i_index[8], j_index[8];
  double U0Min, U0Max, UQuad[4], phi_n;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  LevelSet2DState DU, DUDx_ave, DUDy_ave;

  // Carry out the limited solution reconstruction in each cell of the
  // computational mesh.

  if (i == ICl-2 || i == ICu+2 || j == JCl-2 || j == JCu+2) {
    n_pts = 0;
  } else if (i == ICl-1 && Grid.BCtypeW[j] != BC_NONE) {
    if (j == JCl-1 || j == JCu+1) {
      n_pts = 0;
    } else if (Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
	       Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION) {
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
    } else if (Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
	       Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION) {
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
    } else if (Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
	       Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION) {
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
    } else if (Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
	       Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION) {
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
    DUDx_ave = LevelSet2D_ZERO;
    DUDy_ave = LevelSet2D_ZERO;
    DxDx_ave = ZERO;
    DxDy_ave = ZERO;
    DyDy_ave = ZERO;

    for (n2 = 0; n2 < n_pts; n2++) {
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

    for (n = 1; n <= NUM_VAR_LEVELSET2D; n++) {
      U0Min = U[i][j][n];
      U0Max = U0Min;
      for (n2 = 0; n2 <= n_pts-1; n2++) {
	U0Min = min(U0Min,U[i_index[n2]][j_index[n2]][n]);
	U0Max = max(U0Max,U[i_index[n2]][j_index[n2]][n]);
      }

      dX       = Grid.xfaceE(i,j) - Grid.Cell[i][j].Xc;
      UQuad[0] = U[i][j][n] + dUdx[i][j][n]*dX.x + dUdy[i][j][n]*dX.y;
      dX       = Grid.xfaceW(i,j) - Grid.Cell[i][j].Xc;
      UQuad[1] = U[i][j][n] + dUdx[i][j][n]*dX.x + dUdy[i][j][n]*dX.y;
      dX       = Grid.xfaceN(i,j) - Grid.Cell[i][j].Xc;
      UQuad[2] = U[i][j][n] + dUdx[i][j][n]*dX.x + dUdy[i][j][n]*dX.y;
      dX       = Grid.xfaceS(i,j) - Grid.Cell[i][j].Xc;
      UQuad[3] = U[i][j][n] + dUdx[i][j][n]*dX.x + dUdy[i][j][n]*dX.y;

      switch(Limiter) {
      case LIMITER_ONE :
	phi_n = ONE;
	break;
      case LIMITER_ZERO :
	phi_n = ZERO;
	break;
      case LIMITER_BARTH_JESPERSEN :
	phi_n = Limiter_BarthJespersen(UQuad,U[i][j][n],U0Min,U0Max,4);
	break;
      case LIMITER_VENKATAKRISHNAN :
	phi_n = Limiter_Venkatakrishnan(UQuad,U[i][j][n],U0Min,U0Max,4);
	break;
      case LIMITER_VANLEER :
	phi_n = Limiter_VanLeer(UQuad,U[i][j][n],U0Min,U0Max,4);
	break;
      case LIMITER_VANALBADA :
	phi_n = Limiter_VanAlbada(UQuad,U[i][j][n],U0Min,U0Max,4);
	break;
      default:
	phi_n = Limiter_BarthJespersen(UQuad,U[i][j][n],U0Min,U0Max,4);
	break;
      }

      phi[i][j][n] = phi_n;

    }

  } else {
    dUdx[i][j] = LevelSet2D_ZERO;
    dUdy[i][j] = LevelSet2D_ZERO; 
    phi[i][j]  = LevelSet2D_ZERO;
  }
    
}

/**********************************************************************
 * LevelSet2D_Quad_Block::LoadSendBuffer_Flux_F2C -- Loads send       *
 *                                       message buffer for fine to   *
 *                                       coarse block message passing *
 *                                       of conservative solution     *
 *                                       fluxes.                      *
 **********************************************************************/
inline int LevelSet2D_Quad_Block::LoadSendBuffer_Flux_F2C(double *buffer,
                                                          int &buffer_count,
                                                          const int buffer_size,
                                                          const int i_min,
                                                          const int i_max,
                                                          const int i_inc,
                                                          const int j_min,
                                                          const int j_max,
                                                          const int j_inc) {

  if (j_min == j_max && j_min == JCl) {
    for (int i = i_min; ((i_inc+2)/4) ? (i < i_max):(i > i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = FluxS[i  ][k] + FluxS[i+1][k];
      }
    }
  } else if (j_min == j_max && j_min == JCu) {
    for (int i = i_min; ((i_inc+2)/4) ? (i < i_max):(i > i_max); i += i_inc) {
      for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = FluxN[i  ][k] + FluxN[i+1][k];
      }
    }
  } else if (i_min == i_max && i_min == ICl) {
    for (int j = j_min; ((j_inc+2)/4) ? (j < j_max):(j > j_max); j += j_inc) {
      for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = FluxW[j  ][k] + FluxW[j+1][k];
      }
    }
  } else if (i_min == i_max && i_min == ICu) {
    for (int j = j_min; ((j_inc+2)/4) ? (j < j_max):(j > j_max); j += j_inc) {
      for (int k = 1; k <= NUM_VAR_LEVELSET2D; k++) {
	buffer_count++;
	if (buffer_count >= buffer_size) return 1;
	buffer[buffer_count] = FluxE[j  ][k] + FluxE[j+1][k];
      }
    }
  }

  return 0;

}

/**********************************************************************
 * LevelSet2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C -- Unloads     *
 *                                       receive message buffer for   *
 *                                       fine to coarse block message *
 *                                       passing of conservative      *
 *                                       solution fluxes.             *
 **********************************************************************/
inline int LevelSet2D_Quad_Block::UnloadReceiveBuffer_Flux_F2C(double *buffer,
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
      buffer_count = buffer_count + NUM_VAR_LEVELSET2D;
      if (buffer_count >= buffer_size) return 1;
      FluxS[i] = - LevelSet2DState(buffer[buffer_count-3],
				   buffer[buffer_count-2],
				   buffer[buffer_count-1],
				   buffer[buffer_count  ]) - FluxS[i];
    }
  } else if (j_min == j_max && j_min == JCu) {
    for (int i = i_min; ((i_inc+1)/2) ? (i <= i_max):(i >= i_max); i += i_inc) {
      buffer_count = buffer_count + NUM_VAR_LEVELSET2D;
      if (buffer_count >= buffer_size) return 1;
      FluxN[i] = - LevelSet2DState(buffer[buffer_count-3],
				   buffer[buffer_count-2],
				   buffer[buffer_count-1],
				   buffer[buffer_count  ]) - FluxN[i];
    }
  } else if (i_min == i_max && i_min == ICl) {
    for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
      buffer_count = buffer_count + NUM_VAR_LEVELSET2D;
      if (buffer_count >= buffer_size) return 1;
      FluxW[j] = - LevelSet2DState(buffer[buffer_count-3],
				   buffer[buffer_count-2],
				   buffer[buffer_count-1],
				   buffer[buffer_count  ]) - FluxW[j];
    }
  } else if (i_min == i_max && i_min == ICu) {
    for (int j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
      buffer_count = buffer_count + NUM_VAR_LEVELSET2D;
      if (buffer_count >= buffer_size) return 1;
      FluxE[j] = - LevelSet2DState(buffer[buffer_count-3],
				   buffer[buffer_count-2],
				   buffer[buffer_count-1],
				   buffer[buffer_count  ]) - FluxE[j];
    }
  }
  return 0;
}

/**********************************************************************
 * LevelSet2D_Quad_Block -- Single Block External Subroutines.        *
 **********************************************************************/

extern void Broadcast_Solution_Block(LevelSet2D_Quad_Block &SolnBlk);

#ifdef _MPI_VERSION
extern void Broadcast_Solution_Block(LevelSet2D_Quad_Block &SolnBlk,
                                     MPI::Intracomm &Communicator,
                                     const int Source_CPU);
#endif

extern void Copy_Solution_Block(LevelSet2D_Quad_Block &SolnBlk1,
		                LevelSet2D_Quad_Block &SolnBlk2);

extern int Prolong_Solution_Block(LevelSet2D_Quad_Block &SolnBlk_Fine,
				  LevelSet2D_Quad_Block &SolnBlk_Original,
				  const int Sector);

extern int Restrict_Solution_Block(LevelSet2D_Quad_Block &SolnBlk_Coarse,
				   LevelSet2D_Quad_Block &SolnBlk_Original_SW,
				   LevelSet2D_Quad_Block &SolnBlk_Original_SE,
				   LevelSet2D_Quad_Block &SolnBlk_Original_NW,
				   LevelSet2D_Quad_Block &SolnBlk_Original_NE);

extern int Construct_Bulk_Flow_Field(LevelSet2D_Quad_Block &SolnBlk,
				     LevelSet2D_Input_Parameters &IP);

extern void BCs(LevelSet2D_Quad_Block &SolnBlk);

extern void Set_Global_TimeStep(LevelSet2D_Quad_Block &SolnBlk,
                                const double &Dt_min);

extern double L1_Norm_Residual(LevelSet2D_Quad_Block &SolnBlk, const int &n);

extern double L2_Norm_Residual(LevelSet2D_Quad_Block &SolnBlk, const int &n);

extern double Max_Norm_Residual(LevelSet2D_Quad_Block &SolnBlk, const int &n);

extern int Store_Initial_Eikonal_Solution(LevelSet2D_Quad_Block &SolnBlk);

extern int Calculate_Sign_Function(LevelSet2D_Quad_Block &SolnBlk,
				   LevelSet2D_Input_Parameters &Input_Parameters);

extern void Linear_Reconstruction_GreenGauss(LevelSet2D_Quad_Block &SolnBlk,
                                             const int i,
                                             const int j,
					     const int Limiter);

extern void Linear_Reconstruction_GreenGauss(LevelSet2D_Quad_Block &SolnBlk,
					     const int Limiter);

extern void Linear_Reconstruction_LeastSquares(LevelSet2D_Quad_Block &SolnBlk,
                                               const int i,
                                               const int j,
					       const int Limiter);

extern void Linear_Reconstruction_LeastSquares(LevelSet2D_Quad_Block &SolnBlk,
					       const int Limiter);

extern void Laplacian_Reconstruction_GreenGauss(LevelSet2D_Quad_Block &SolnBlk,
						const int i,
						const int j,
						const int n,
						double &ddUdd);

extern double dUdx(LevelSet2D_Quad_Block &SolnBlk, const int &i, const int &j, const int &n);

extern double dUdy(LevelSet2D_Quad_Block &SolnBlk, const int &i, const int &j, const int &n);

extern double ddUdx2(LevelSet2D_Quad_Block &SolnBlk, const int &i, const int &j, const int &n);

extern double ddUdy2(LevelSet2D_Quad_Block &SolnBlk, const int &i, const int &j, const int &n);

extern double dddUdx3(LevelSet2D_Quad_Block &SolnBlk, const int &i, const int &j, const int &n);

extern double dddUdy3(LevelSet2D_Quad_Block &SolnBlk, const int &i, const int &j, const int &n);

extern void Reconstruction_EssentiallyNonOscillatory(LevelSet2D_Quad_Block &SolnBlk,
						     const int i,
						     const int j,
						     const int n,
						     const int i_Reconstruction);

extern void Reconstruction_EssentiallyNonOscillatory(LevelSet2D_Quad_Block &SolnBlk,
						     const int n,
						     const int i_Reconstruction);

extern void Reconstruction_WeightedEssentiallyNonOscillatory(LevelSet2D_Quad_Block &SolnBlk,
							     const int i,
							     const int j,
							     const int n);

extern double Reconstruction_WeightedEssentiallyNonOscillatory(const double v1,
							       const double v2,
							       const double v3,
							       const double v4,
							       const double v5);

extern void Reconstruction_WeightedEssentiallyNonOscillatory(LevelSet2D_Quad_Block &SolnBlk,
							     const int n);

extern void Reconstruction_Curvature(LevelSet2D_Quad_Block &SolnBlk,
				     LevelSet2D_Input_Parameters &IP,
				     const int n);

extern void Reconstruction_Curvature_Laplacian(LevelSet2D_Quad_Block &SolnBlk,
					       const int i,
					       const int j,
					       const int n);

extern void Reconstruction_Curvature_Regular(LevelSet2D_Quad_Block &SolnBlk,
					     const int i,
					     const int j,
					     const int n);

extern void Calculate_Refinement_Criteria(double *refinement_criteria,
					  LevelSet2D_Input_Parameters &IP,
					  int &number_refinement_criteria,
					  LevelSet2D_Quad_Block &SolnBlk);

extern void Fix_Refined_Block_Boundaries(LevelSet2D_Quad_Block &SolnBlk,
                                         const int Fix_North_Boundary,
                                         const int Fix_South_Boundary,
                                         const int Fix_East_Boundary,
                                         const int Fix_West_Boundary);

extern void Unfix_Refined_Block_Boundaries(LevelSet2D_Quad_Block &SolnBlk);

/**********************************************************************
 * LevelSet2D_Quad_Block -- Multiple Block External Subroutines.      *
 **********************************************************************/

extern LevelSet2D_Quad_Block* Allocate(LevelSet2D_Quad_Block *Soln_ptr,
				       LevelSet2D_Input_Parameters &Input_Parameters);

extern LevelSet2D_Quad_Block* Deallocate(LevelSet2D_Quad_Block *Soln_ptr,
					 LevelSet2D_Input_Parameters &Input_Parameters);

extern void ICs(LevelSet2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
                LevelSet2D_Input_Parameters &Input_Parameters);

extern int Construct_Bulk_Flow_Field(LevelSet2D_Quad_Block *Soln_ptr,
				     AdaptiveBlock2D_List &Soln_Block_List,
				     LevelSet2D_Input_Parameters &Input_Parameters);

extern void BCs(LevelSet2D_Quad_Block *Soln_ptr,
                AdaptiveBlock2D_List &Soln_Block_List,
		LevelSet2D_Input_Parameters &Input_Parameters);

// extern int Linear_Reconstruction(LevelSet2D_Quad_Block *Soln_ptr,
// 				 AdaptiveBlock2D_List &Soln_Block_List,
// 				 LevelSet2D_Input_Parameters &Input_Parameters);

extern int Reconstruction_EssentiallyNonOscillatory(LevelSet2D_Quad_Block *Soln_ptr,
						    AdaptiveBlock2D_List &Soln_Block_List,
						    LevelSet2D_Input_Parameters &Input_Parameters,
						    const int n);

extern int Reconstruction_WeightedEssentiallyNonOscillatory(LevelSet2D_Quad_Block *Soln_ptr,
							    AdaptiveBlock2D_List &Soln_Block_List,
							    const int n);

extern int Reconstruction_Curvature(LevelSet2D_Quad_Block *Soln_ptr,
				    AdaptiveBlock2D_List &Soln_Block_List,
				    LevelSet2D_Input_Parameters &Input_Parameters,
				    const int n);

extern void Set_Global_TimeStep(LevelSet2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                const double &Dt_min);

extern double L1_Norm_Residual(LevelSet2D_Quad_Block *Soln_ptr,
			       AdaptiveBlock2D_List &Soln_Block_List,
			       const int &n);

extern double L2_Norm_Residual(LevelSet2D_Quad_Block *Soln_ptr,
			       AdaptiveBlock2D_List &Soln_Block_List,
			       const int &n);

extern double Max_Norm_Residual(LevelSet2D_Quad_Block *Soln_ptr,
				AdaptiveBlock2D_List &Soln_Block_List,
				const int &n);

extern int Store_Initial_Eikonal_Solution(LevelSet2D_Quad_Block *Soln_ptr,
					  AdaptiveBlock2D_List &Soln_Block_List);

extern int Calculate_Sign_Function(LevelSet2D_Quad_Block *Soln_ptr,
				   AdaptiveBlock2D_List &Soln_Block_List,
				   LevelSet2D_Input_Parameters &Input_Parameters);

extern int Initial_Adaptive_Mesh_Refinement(LevelSet2D_Quad_Block *Soln_ptr,
					    LevelSet2D_Input_Parameters &Input_Parameters,
					    QuadTreeBlock_DataStructure &QuadTree,
					    AdaptiveBlockResourceList &GlobalSolnBlockList,
					    AdaptiveBlock2D_List &Soln_Block_List);

extern int Initial_Adaptive_Mesh_Refinement(LevelSet2D_Quad_Block *Soln_ptr,
					    LevelSet2D_Input_Parameters &Input_Parameters,
					    QuadTreeBlock_DataStructure &QuadTree,
					    AdaptiveBlockResourceList &GlobalSolnBlockList,
					    AdaptiveBlock2D_List &Soln_Block_List,
					    const Interface2D_List &Interface_List);

extern int Uniform_Adaptive_Mesh_Refinement(LevelSet2D_Quad_Block *Soln_ptr,
					    LevelSet2D_Input_Parameters &Input_Parameters,
					    QuadTreeBlock_DataStructure &QuadTree,
					    AdaptiveBlockResourceList &GlobalSolnBlockList,
					    AdaptiveBlock2D_List &Soln_Block_List);

extern int Uniform_Adaptive_Mesh_Refinement(LevelSet2D_Quad_Block *Soln_ptr,
					    LevelSet2D_Input_Parameters &Input_Parameters,
					    QuadTreeBlock_DataStructure &QuadTree,
					    AdaptiveBlockResourceList &GlobalSolnBlockList,
					    AdaptiveBlock2D_List &Soln_Block_List,
					    const Interface2D_List &Interface_List);

/**********************************************************************
 * LevelSet2D_Quad_Block -- Interface Single Block External           *
 *                          Subroutines                               *
 **********************************************************************/


extern int Set_Interface_List(LevelSet2D_Quad_Block &SolnBlk,
			      const Interface2D_List &Interface_List);

extern int Initialize_Interfaces(LevelSet2D_Quad_Block &SolnBlk,
				 LevelSet2D_Input_Parameters &IP);

extern int Initialize_Interfaces(LevelSet2D_Quad_Block &SolnBlk,
				 LevelSet2D_Input_Parameters &IP,
				 const Interface2D_List &Interface_List);

extern int Exact_Initial_Extension(LevelSet2D_Quad_Block &SolnBlk,
				   LevelSet2D_Input_Parameters &IP);

extern int Geometric_Extension_Problem(LevelSet2D_Quad_Block &SolnBlk,
				       LevelSet2D_Input_Parameters &IP);

extern int Scalar_Geometric_Extension_Problem(LevelSet2D_Quad_Block &SolnBlk,
					      LevelSet2D_Input_Parameters &IP);

#ifndef _RETRIEVE_DEBUG_
extern int Retrieve_Interface_Spline(LevelSet2D_Quad_Block &SolnBlk,
				     const int &gblknum);

extern int Trace_Interface_Spline(LevelSet2D_Quad_Block &SolnBlk,
				  LinkedList<Vector2D> &p,
				  LinkedList<Vector2D> &F,
				  const double &epsilon,
				  int &npts,
				  const int &start);
#endif

#ifdef _RETRIEVE_DEBUG_
extern int Retrieve_Interface_Spline(LevelSet2D_Quad_Block &SolnBlk,
				     const int &gblknum,
				     ofstream &dout);

extern int Trace_Interface_Spline(LevelSet2D_Quad_Block &SolnBlk,
				  LinkedList<Vector2D> &p,
				  LinkedList<Vector2D> &F,
				  const double &epsilon,
				  int &npts,
				  const int &start,
				  ofstream &dout);
#endif

extern void Update_Start_List(LevelSet2D_Quad_Block &SolnBlk,
			      LinkedList<int> &start);

extern void Flag_Infected_Cell(LevelSet2D_Quad_Block &SolnBlk,
			       const int &ic,
			       const int &jc,
			       const double &epsilon,
			       LinkedList<int> &start);

extern void Clean_Infected_Cell(LevelSet2D_Quad_Block &SolnBlk,
				const int &ic,
				const int &jc,
				const int &face);

/**********************************************************************
 * LevelSet2D_Quad_Block -- Interface Multi Block External            *
 *                          Subroutines                               *
 **********************************************************************/

extern int Set_Interface_List(LevelSet2D_Quad_Block *Soln_ptr,
			      AdaptiveBlock2D_List &Soln_Block_List,
			      const Interface2D_List &Interface_List);

extern int Initialize_Interfaces(LevelSet2D_Quad_Block *Soln_ptr,
				 AdaptiveBlock2D_List &Soln_Block_List,
				 LevelSet2D_Input_Parameters &Input_Parameters);

extern int Initialize_Interfaces(LevelSet2D_Quad_Block *Soln_ptr,
				 AdaptiveBlock2D_List &Soln_Block_List,
				 LevelSet2D_Input_Parameters &Input_Parameters,
				 const Interface2D_List &Interface_List);

extern int Exact_Initial_Extension(LevelSet2D_Quad_Block *Soln_ptr,
				   AdaptiveBlock2D_List &Soln_Block_List,
				   LevelSet2D_Input_Parameters &Input_Parameters);

extern int Geometric_Extension_Problem(LevelSet2D_Quad_Block *Soln_ptr,
				       AdaptiveBlock2D_List &Soln_Block_List,
				       LevelSet2D_Input_Parameters &Input_Parameters);

extern int Scalar_Geometric_Extension_Problem(LevelSet2D_Quad_Block *Soln_ptr,
					      AdaptiveBlock2D_List &Soln_Block_List,
					      LevelSet2D_Input_Parameters &Input_Parameters);

extern int Retrieve_Interface_Spline(LevelSet2D_Quad_Block *Soln_ptr,
				     AdaptiveBlock2D_List &Soln_Block_List);

extern int Share_Interface_Information(LevelSet2D_Quad_Block *Soln_ptr,
				       QuadTreeBlock_DataStructure &QuadTree,
				       AdaptiveBlockResourceList &GlobalSolnBlockList,
				       AdaptiveBlock2D_List &Soln_Block_List,
				       LevelSet2D_Input_Parameters &Input_Parameters);

/**********************************************************************
 * LevelSet2D_Quad_Block -- Hamilton-Jacobi Single Block External     *
 *                          Subroutines                               *
 **********************************************************************/

extern double CFL_Hamilton_Jacobi(LevelSet2D_Quad_Block &SolnBlk,
				  LevelSet2D_Input_Parameters &IP);

extern int dUdt_Multistage_Hamilton_Jacobi(LevelSet2D_Quad_Block &SolnBlk,
					   const int i_stage,
					   LevelSet2D_Input_Parameters &IP);

extern int Update_Solution_Multistage_Hamilton_Jacobi(LevelSet2D_Quad_Block &SolnBlk,
						      const int i_stage,
						      LevelSet2D_Input_Parameters &IP);

/**********************************************************************
 * LevelSet2D_Quad_Block -- Hamilton-Jacobi Multi Block External      *
 *                          Subroutines                               *
 **********************************************************************/

extern double CFL_Hamilton_Jacobi(LevelSet2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  LevelSet2D_Input_Parameters &Input_Parameters);

extern int dUdt_Multistage_Hamilton_Jacobi(LevelSet2D_Quad_Block *Soln_ptr,
					   AdaptiveBlock2D_List &Soln_Block_List,
					   LevelSet2D_Input_Parameters &Input_Parameters,
					   const int I_Stage);

extern int Update_Solution_Multistage_Hamilton_Jacobi(LevelSet2D_Quad_Block *Soln_ptr,
						      AdaptiveBlock2D_List &Soln_Block_List,
						      LevelSet2D_Input_Parameters &Input_Parameters,
						      const int I_Stage);

/**********************************************************************
 * LevelSet2D_Quad_Block -- Eikonal Single Block External Subroutines *
 **********************************************************************/

extern int Eikonal_Error(LevelSet2D_Quad_Block &SolnBlk,
			 LevelSet2D_Input_Parameters &IP,
			 double &block_error,
			 double &block_area);

extern double CFL_Eikonal(LevelSet2D_Quad_Block &SolnBlk);

extern int dUdt_Multistage_Eikonal(LevelSet2D_Quad_Block &SolnBlk,
				   const int i_stage,
				   LevelSet2D_Input_Parameters &IP);

extern int Update_Multistage_Eikonal(LevelSet2D_Quad_Block &SolnBlk,
				     const int i_stage,
				     LevelSet2D_Input_Parameters &IP);

/**********************************************************************
 * LevelSet2D_Quad_Block -- Eikonal Multi Block External Subroutines  *
 **********************************************************************/

extern int Eikonal_Error(LevelSet2D_Quad_Block *Soln_ptr,
			 LevelSet2D_Input_Parameters &IP,
			 AdaptiveBlock2D_List &Soln_Block_List,
			 double &global_error,
			 double &global_area);

extern double CFL_Eikonal(LevelSet2D_Quad_Block *Soln_ptr,
			  AdaptiveBlock2D_List &Soln_Block_List);

extern int Explicit_Eikonal_Equation(LevelSet2D_Quad_Block *Soln_ptr,
				     LevelSet2D_Input_Parameters &Input_Parameters,
				     QuadTreeBlock_DataStructure &QuadTree,
				     AdaptiveBlockResourceList &GlobalSolnBlockList,
				     AdaptiveBlock2D_List &Soln_Block_List,
				     const int Redistance);

extern int dUdt_Multistage_Eikonal(LevelSet2D_Quad_Block *Soln_ptr,
				   AdaptiveBlock2D_List &Soln_Block_List,
				   LevelSet2D_Input_Parameters &Input_Parameters,
				   const int I_Stage);

extern int Update_Multistage_Eikonal(LevelSet2D_Quad_Block *Soln_ptr,
				     AdaptiveBlock2D_List &Soln_Block_List,
				     LevelSet2D_Input_Parameters &Input_Parameters,
				     const int I_Stage);

/**********************************************************************
 * LevelSet2D_Quad_Block -- Scalar Extension Single Block External    *
 *                          Subroutines                               *
 **********************************************************************/

extern double CFL_Scalar_Extension(LevelSet2D_Quad_Block &SolnBlk);

extern int dUdt_Multistage_Scalar_Extension(LevelSet2D_Quad_Block &SolnBlk,
					    const int i_stage,
					    LevelSet2D_Input_Parameters &IP);

extern int Update_Multistage_Scalar_Extension(LevelSet2D_Quad_Block &SolnBlk,
					      const int i_stage,
					      LevelSet2D_Input_Parameters &IP);

/**********************************************************************
 * LevelSet2D_Quad_Block -- Scalar Extension Multi Block External     *
 *                          Subroutines                               *
 **********************************************************************/

extern double CFL_Scalar_Extension(LevelSet2D_Quad_Block *Soln_ptr,
				   AdaptiveBlock2D_List &Soln_Block_List);

extern int Explicit_Scalar_Extension_Equation(LevelSet2D_Quad_Block *Soln_ptr,
					      LevelSet2D_Input_Parameters &Input_Parameters,
					      QuadTreeBlock_DataStructure &QuadTree,
					      AdaptiveBlockResourceList &GlobalSolnBlockList,
					      AdaptiveBlock2D_List &Soln_Block_List);

extern int dUdt_Multistage_Scalar_Extension(LevelSet2D_Quad_Block *Soln_ptr,
					    AdaptiveBlock2D_List &Soln_Block_List,
					    LevelSet2D_Input_Parameters &Input_Parameters,
					    const int I_Stage);

extern int Update_Multistage_Scalar_Extension(LevelSet2D_Quad_Block *Soln_ptr,
					      AdaptiveBlock2D_List &Soln_Block_List,
					      LevelSet2D_Input_Parameters &Input_Parameters,
					      const int I_Stage);

/**********************************************************************
 * LevelSet2D_Quad_Block -- Input and Output Single Block External    *
 *                          Subroutines.                              *
 **********************************************************************/

extern void Write_Solution_Block(LevelSet2D_Quad_Block &SolnBlk,
	                         ostream &Out_File);

extern void Read_Solution_Block(LevelSet2D_Quad_Block &SolnBlk,
	                        istream &In_File);

extern void Output_Tecplot(LevelSet2D_Quad_Block &SolnBlk,
		           const int Number_of_Time_Steps,
                           const double &Time,
                           const int Block_Number,
                           const int Output_Title,
	                   ostream &Out_File);

extern void Output_Cells_Tecplot(LevelSet2D_Quad_Block &SolnBlk,
		                 const int Number_of_Time_Steps,
                                 const double &Time,
				 const int CPU_Number,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_Nodes_Tecplot(LevelSet2D_Quad_Block &SolnBlk,
                                 const int Block_Number,
                                 const int Output_Title,
	                         ostream &Out_File);

extern void Output_Interface_Tecplot(LevelSet2D_Quad_Block &SolnBlk,
				     const int Number_of_Time_Steps,
				     const double &Time,
				     const int Block_Number,
				     const int Output_Title,
				     ostream &Out_File);

extern void Output_Circle_Tecplot(LevelSet2D_Quad_Block &SolnBlk,
				  LevelSet2D_Input_Parameters &IP,
				  const double &Time,
				  const int Block_Number,
				  const int Output_Title,
				  ostream &Out_File,
				  double &l1_norm,
				  double &l2_norm,
				  double &max_norm);

extern void Output_Ellipse_Tecplot(LevelSet2D_Quad_Block &SolnBlk,
				   LevelSet2D_Input_Parameters &IP,
				   const double &Time,
				   const int Block_Number,
				   const int Output_Title,
				   ostream &Out_File,
				   double &l1_norm,
				   double &l2_norm,
				   double &max_norm);

extern void Output_Zalesaks_Disk_Tecplot(LevelSet2D_Quad_Block &SolnBlk,
					 LevelSet2D_Input_Parameters &IP,
					 const double &Time,
					 const int Block_Number,
					 const int Output_Title,
					 ostream &Out_File,
					 double &l1_norm,
					 double &l2_norm,
					 double &max_norm);

/**********************************************************************
 * LevelSet2D_Quad_Block -- Input and Output Multiple Block External  *
 *                          Subroutines.                              *
 **********************************************************************/

extern int Read_Restart_Solution(LevelSet2D_Quad_Block *Soln_ptr,
                                 AdaptiveBlock2D_List &Soln_Block_List,
                                 LevelSet2D_Input_Parameters &Input_Parameters,
		                 int &Number_of_Time_Steps,
                                 double &Time,
                                 CPUTime &CPU_Time);

extern int Write_Restart_Solution(LevelSet2D_Quad_Block *Soln_ptr,
                                  AdaptiveBlock2D_List &Soln_Block_List,
                                  LevelSet2D_Input_Parameters &Input_Parameters,
		                  const int Number_of_Time_Steps,
                                  const double &Time,
                                  const CPUTime &CPU_Time);

extern int Output_Tecplot(LevelSet2D_Quad_Block *Soln_ptr,
                          AdaptiveBlock2D_List &Soln_Block_List,
                          LevelSet2D_Input_Parameters &Input_Parameters,
		          const int Number_of_Time_Steps,
                          const double &Time);

extern int Output_Cells_Tecplot(LevelSet2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                LevelSet2D_Input_Parameters &Input_Parameters,
		                const int  Number_of_Time_Steps,
                                const double &Time);

extern int Output_Nodes_Tecplot(LevelSet2D_Quad_Block *Soln_ptr,
                                AdaptiveBlock2D_List &Soln_Block_List,
                                LevelSet2D_Input_Parameters &Input_Parameters);

extern int Output_Mesh_Tecplot(LevelSet2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               LevelSet2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Mesh_Gnuplot(LevelSet2D_Quad_Block *Soln_ptr,
                               AdaptiveBlock2D_List &Soln_Block_List,
                               LevelSet2D_Input_Parameters &Input_Parameters,
		               const int Number_of_Time_Steps,
                               const double &Time);

extern int Output_Interface_Tecplot(LevelSet2D_Quad_Block *Soln_ptr,
				    AdaptiveBlock2D_List &Soln_Block_List,
				    LevelSet2D_Input_Parameters &Input_Parameters,
				    const int Number_of_Time_Steps,
				    const double &Time);

extern int Output_Circle_Tecplot(LevelSet2D_Quad_Block *Soln_ptr,
				 AdaptiveBlock2D_List &Soln_Block_List,
				 LevelSet2D_Input_Parameters &Input_Parameters,
				 const int Number_of_Time_Steps,
				 const double &Time);

extern int Output_Ellipse_Tecplot(LevelSet2D_Quad_Block *Soln_ptr,
				  AdaptiveBlock2D_List &Soln_Block_List,
				  LevelSet2D_Input_Parameters &Input_Parameters,
				  const int Number_of_Time_Steps,
				  const double &Time);

extern int Output_Zalesaks_Disk_Tecplot(LevelSet2D_Quad_Block *Soln_ptr,
					AdaptiveBlock2D_List &Soln_Block_List,
					LevelSet2D_Input_Parameters &Input_Parameters,
					const int Number_of_Time_Steps,
					const double &Time);

/**********************************************************************
 * LevelSet2D_Quad_Block -- Multiple Block External Mesh Subroutines. *
 **********************************************************************/

extern Grid2D_Quad_Block** Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                            LevelSet2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Broadcast_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                                      LevelSet2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
                                             LevelSet2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid_Definition(Grid2D_Quad_Block **Grid_ptr,
							    LevelSet2D_Input_Parameters &Input_Parameters);

extern int Write_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
                                  LevelSet2D_Input_Parameters &Input_Parameters);

extern Grid2D_Quad_Block** Read_Multi_Block_Grid(Grid2D_Quad_Block **Grid_ptr,
						 LevelSet2D_Input_Parameters &Input_Parameters);

extern int Output_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                          LevelSet2D_Input_Parameters &Input_Parameters);

extern int Output_Cells_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                LevelSet2D_Input_Parameters &Input_Parameters);

extern int Output_Nodes_Tecplot(Grid2D_Quad_Block **Grid_ptr,
                                LevelSet2D_Input_Parameters &Input_Parameters);

/**********************************************************************
 * LevelSet2D_Quad_Block -- Solvers.                                  *
 **********************************************************************/

extern int LevelSet2DQuadSolver(char *Input_File_Name_ptr, int batch_flag);

extern LevelSet2D_Quad_Block* Initialize_Level_Set_Solution(char *Input_File_Name,
							    const int &batch_flag,
							    int &error_flag,
							    LevelSet2D_Quad_Block *Local_SolnBlk,
							    LevelSet2D_Input_Parameters &Input_Parameters,
							    QuadTreeBlock_DataStructure &QuadTree,
							    AdaptiveBlockResourceList &List_of_Global_Solution_Blocks,
							    AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
							    const Interface2D_List &Interface_List);

extern LevelSet2D_Quad_Block* Restart_Level_Set_Solution(char *Input_File_Name,
							 const int &batch_flag,
							 int &error_flag,
							 int &number_of_time_steps,
							 double &Time,
							 LevelSet2D_Quad_Block *Local_SolnBlk,
							 LevelSet2D_Input_Parameters &Input_Parameters,
							 QuadTreeBlock_DataStructure &QuadTree,
							 AdaptiveBlockResourceList &List_of_Global_Solution_Blocks,
							 AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
							 const Interface2D_List &Interface_List);

extern int Evolve_Level_Set_Solution(const int &batch_flag,
				     LevelSet2D_Quad_Block *Local_SolnBlk,
				     LevelSet2D_Input_Parameters &Input_Parameters,
				     QuadTreeBlock_DataStructure &QuadTree,
				     AdaptiveBlockResourceList &List_of_Global_Solution_Blocks,
				     AdaptiveBlock2D_List &List_of_Local_Solution_Blocks,
				     double &Time,
				     const double &Time_Max,
				     int &number_of_time_steps);

#endif // _LEVELSET2D_QUAD_INCLUDED
