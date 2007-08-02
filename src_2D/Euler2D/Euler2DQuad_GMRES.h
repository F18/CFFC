#ifndef _EULER2D_QUAD_GMRES_INCLUDED
#define _EULER2D_QUAD_GMRES_INCLUDED

#include <iostream>
#include <cstdio>
#include <cmath> 
using namespace std;

#include "BpMatrix.h"
#include "BpPrecon.h"
#include "BpResource.h"
#include "blas1.h"
#include "BlockMat.h"

#ifndef _BLOCKMAT_H_
#include "BlockMat.h"
#endif //_BLOCKMAT_H_

#ifndef _DENSEMAT_H_
#include "DenseMat.h"
#endif //_DENSEMAT_H_

#ifndef _HBTMAT_H_
#include "HBTMat.h"
#endif //_HBTMAT_H_

#ifndef _BILUK_H_
#include "BILUK.h"
#endif //_BILUK_H_

#ifndef _BRELAX_H_
#include "BRelax.h"
#endif //_BRELAX_H_

/* Include 2D Euler quadrilateral mesh solution header file. */

#ifndef _EULER2D_QUAD_INCLUDED
#include "Euler2DQuad.h"
#endif // _EULER2D_QUAD_INCLUDED

/***********************************************************
 * Class: GMRES_Block                                      *
 *                                                         *
 * Member functions                                        *
 *       s      -- Return ..........                       *
 *      cs      -- Return cos vector                       *
 *      sn      -- Return sin vector                       *
 *       W      -- Return Az vector                        *
 *       z      -- Return inversion of preconditioner      *
 *                  times v vector                         *
 *       b      -- Return RHS vector                       *
 *      Q1      -- Return matrix free vector               *
 *      Q2      -- Return matrix free vector               *
 *       x      -- Return solution vector, delta u         *
 *                 (Note: it also contains initial guess   * 
 *                        of zeros.)                       *
 *       V      -- Return Krylov search direction vector   *
 *     tmp      -- Return a temporary vector               *
 *       Q      -- Return a temporary vector for matrix-   *
 *                 free calculation                        *
 *       U      -- Return original solution vector         *
 *       H      -- Return Hessenberg matrix                *
 *  dt_vec      -- Return finite time step vector          *
 *   dTime      -- Return global time step                 *
 *       m      -- Return restart                          *
 *    xpts      -- Return number of nodes in x-direction   *
 *    ypts      -- Return number of nodes in y-direction   *   
 * overlap      -- Return level of overlap                 *
 * P_Switch     -- Return preconditioner switch            *
 * blocksize    -- Return number of variables              *
 * scalar_dim   -- Return xpts * ypts * blocksize          *
 * max_gmres_iter -- Return maximun number of gmres        *
 *                   iterations                            *
 *  copy_residual -- Copy residual from SolnBlk to an      *
 *                   1D array                              *
 *                                                         *
 * ** MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING **     *
 * vector_swtich -- Return vector switch                   *
 *                  Note: 0 for x vector (solution)        *
 *                        1 for z vector (z = Minv * v)    *
 * search_directions -- Return number of search            *
 *                      directions                         *
 *     NCi      -- Return number of cells in               *
 *                 the i-direction (zeta-direction)        *
 *     ICl      -- Return lower index for cells in         *
 *                 the i-direction (zeta-direction)        *
 *     ICu      -- Return upper index for cells in         *
 *                 the i-direction (zeta-direction)        *
 *     NCj      -- Return number of cells in               *
 *                 the j-direction (eta-direction)         *
 *     JCl      -- Return lower index for cells in         *
 *                 the j-direction (eta-direction)         *
 *     JCu      -- Return upper index for cells in         *
 *                 the j-direction (eta-direction)         *
 *    Grid      -- Return a dummy object.                  *
 *                 (Required by message passing routines)  *   
 *  Nghost      -- Return number of ghost (halo or         *
 *                 overlap) cells.                         *
 *  SolBlk_ptr  -- Return pointer to solution block.       *
 *                                                         *
 * Member operators                                        *
 *      G -- a GMRES_Block variable object                 *
 *                                                         *
 * G = G;                                                  *
 ***********************************************************/
class GMRES_Block{
private:
public:
  double *                      s; // 
  double *                     cs; // cos vector
  double *                     sn; // sin vector
  double *                      W; // Az
  double *                      b; // RHS vector
  double *                     Q1; // matrix free vector
  double *                     Q2; // matrix free vector
  double *                      x; // initial guess of delta u
  double *                      V; // Krylov search direction vector
  double *                    tmp; // temporary pointer
  double *                      Q; // temporary vector
  double *                      U; // original solution vector
  double *                      H; // Hessenberg matrix
  double *                 dt_vec; // finite time step vector
  double                    dTime; // global time step
  int                           m; // restart
  int                        xpts; // number of nodes in x-direction
  int                        ypts; // number of nodes in y-direction
  int                     overlap; // level of overlap
  int                    P_Switch; // preconditioner switch
  int                   blocksize; // number of variables
  int                  scalar_dim; // xpts * ypts * blocksize
  int              max_gmres_iter; // maximun number of gmres iterations
  
  /* MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING. */
  int               vector_switch; // to select the specified vector 
                                   // for message passing
  int           search_directions; // number of search directions
  int                 NCi,ICl,ICu; // i-direction cell counters.
  int                 NCj,JCl,JCu; // j-direction cell counters.
  int                      Nghost; // Number of ghost cells.
  Grid2D_Quad_Block          Grid; // dummy pointer.
  Euler2D_Quad_Block  *SolBlk_ptr; // Pointer to solution block. 
                                   // Made public so can access them.
  
  /* Creation and copy constructors. */
  GMRES_Block(void) {
    dTime          = 0;
    m              = 0; 
    xpts           = 0; 
    ypts           = 0;
    Nghost         = 0;
    overlap        = 0;
    P_Switch       = 1;
    blocksize      = 0;
    scalar_dim     = 0;
    max_gmres_iter = 0;
    tmp = NULL;

    vector_switch = 99;
    NCi = 0; 
    ICl = 0; ICu = 0; 
    NCj = 0; 
    JCl = 0; JCu = 0;
    SolBlk_ptr = NULL;

    search_directions = 0;
  }
  
  GMRES_Block(const GMRES_Block &G) {
    s  = G.s ; cs = G.cs; sn = G.sn;
    W  = G.W ;
    b  = G.b ;  
    V  = G.V ;
    Q1 = G.Q1; 
    Q2 = G.Q2;  
    x  = G.x ;

    tmp     = G.tmp;
    Q       = G.Q;
    H       = G.H;
    U       = G.U;
    dt_vec  = G.dt_vec; 

    dTime          = G.dTime;
    m              = G.m;  
    xpts           = G.xpts;
    ypts           = G.ypts;
    Nghost         = G.Nghost;
    overlap        = G.overlap;
    P_Switch       = G.P_Switch;
    blocksize      = G.blocksize;
    max_gmres_iter = G.max_gmres_iter;

    vector_switch = G.vector_switch;
    NCi = G.NCi; ICl = G.ICl; ICu = G.ICu; 
    NCj = G.NCj; JCl = G.JCl; JCu = G.JCu;
    Grid = G.Grid;
    SolBlk_ptr = G.SolBlk_ptr;

    search_directions = G.search_directions;
  }
  
  /* Destructor. */
  ~GMRES_Block() {} 
  
  /* Assignment operator. */
  // GMRES_Block operator = (const GMRES_Block &G);
  // Use automatically generated assignment operator.
  
  /* Allocate and deallocate memory for GMRES_Block variables. */
  void allocate(const int restart, const int blocksize, int xpts, int ypts);
  void deallocate();

  int  index(int i, int j, int num = 0) {return ((j*xpts+i)*blocksize+num);}
  void copy_residual(double * RHS, Euler2D_Quad_Block SolnBlk);
  
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
 * GMRES_Block::allocate -- Allocate memory.                              *
 **************************************************************************/
inline void GMRES_Block::allocate(const int restart, const int blocksize, 
			          const int xpts, const int ypts)
{
  assert(restart > 1);
  assert(xpts > 1);
  assert(ypts > 1);

  scalar_dim = xpts * ypts * blocksize;

  // Allocate Memory
   s = new double[restart+1];
  cs = new double[restart];
  sn = new double[restart];

   b         = new double[scalar_dim];
  Q1         = new double[scalar_dim];
  Q2         = new double[scalar_dim];
   x         = new double[scalar_dim]; 
   Q         = new double[scalar_dim];
   U         = new double[scalar_dim];
   dt_vec    = new double[scalar_dim];
   H         = new double[restart*(restart+1)];

   W = new double[restart * scalar_dim];
   V = new double[(restart + 1) * scalar_dim];

   // Initialization
   int i;

   for (i=0;i<restart+1;++i) s[i] = 0.0;
   for (i=0;i<restart;++i)  cs[i] = 0.0;
   for (i=0;i<restart;++i)  sn[i] = 0.0;

   for (i=0;i<scalar_dim;++i)         b[i] = 0.0;
   for (i=0;i<scalar_dim;++i)        Q1[i] = 0.0;
   for (i=0;i<scalar_dim;++i)        Q2[i] = 0.0;
   for (i=0;i<scalar_dim;++i)         x[i] = 0.0;
   for (i=0;i<scalar_dim;++i)         Q[i] = 0.0;
   for (i=0;i<scalar_dim;++i)         U[i] = 0.0;
   for (i=0;i<scalar_dim;++i)    dt_vec[i] = 0.0;

   for (i=0;i<(restart*(restart+1));++i)        H[i] = 0.0;
   for (i=0;i<restart* scalar_dim;++i)          W[i] = 0.0;
   for (i=0;i<((restart + 1) * scalar_dim);++i) V[i] = 0.0;

}

/**************************************************************************
 * GMRES_Block::deallocate -- Deallocate memory.                          *
 **************************************************************************/
inline void GMRES_Block::deallocate() 
{
  delete [] s;
  delete [] cs;
  delete [] sn;
 
  delete [] b;
  delete [] Q1;
  delete [] Q2;
  delete [] dt_vec; 
  delete [] x;
  delete [] Q;
  delete [] U;
  delete [] H;

  delete [] W;
  delete [] V;

}

/**************************************************************************
 * GMRES_Block::copy_residual -- copy residual from SolnBlk  to any       *
 *                               specified 1D vector.                     *
 **************************************************************************/
inline void GMRES_Block::copy_residual(double * RHS, Euler2D_Quad_Block SolnBlk)
{
  int index = 0;
  for (int ycount = SolnBlk.JCl-SolnBlk.Nghost ; 
       ycount <= SolnBlk.JCu+SolnBlk.Nghost ; ycount++) {
    for (int xcount = SolnBlk.ICl-SolnBlk.Nghost ; 
	 xcount <= SolnBlk.ICu+SolnBlk.Nghost ; xcount++) {
      index = (ycount * xpts + xcount) * blocksize;
      RHS[index] = SolnBlk.dUdt[xcount][ycount][0].d;    index++;
      RHS[index] = SolnBlk.dUdt[xcount][ycount][0].dv.x; index++;
      RHS[index] = SolnBlk.dUdt[xcount][ycount][0].dv.y; index++;
      RHS[index] = SolnBlk.dUdt[xcount][ycount][0].E;
    } /* endfor */
  } /* endfor */

}

/*******************************************************************************
 *                                                                             *
 * MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING.                              *
 *                                                                             *
 *******************************************************************************/

/*******************************************************************************
 * GMRES_Block::NumVar -- Returns number of state variables.                   *
 *******************************************************************************/
inline int GMRES_Block::NumVar(void) {
  return (int(NUM_VAR_EULER2D));
}

/*******************************************************************************
 * GMRES_Block::LoadSendBuffer -- Loads send message buffer.                   *
 *******************************************************************************/
inline int GMRES_Block::LoadSendBuffer(double *buffer,
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
        for ( k = 0 ; k < NUM_VAR_EULER2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
	   if (vector_switch) {
              buffer[buffer_count] = W[(search_directions)*scalar_dim + index(i,j,k)];
	   } else {
              buffer[buffer_count] = x[index(i,j,k)];
	   } /* endif */
        } /* endfor */
     } /* endfor */
  } /* endfor */
  return(0);
} 

/*******************************************************************************
 * GMRES_Block::LoadSendBuffer_F2C -- Loads send message buffer for            *
 *                                    fine to coarse block message             *
 *                                    passing.                                 *
 *******************************************************************************/
inline int GMRES_Block::LoadSendBuffer_F2C(double *buffer,
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
        for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
	   if (vector_switch) {
              buffer[buffer_count] = (SolBlk_ptr->Grid.Cell[i  ][j  ].A*
                                      W[(search_directions)*scalar_dim + index(i  ,j  ,k)]+
                                      SolBlk_ptr->Grid.Cell[i+1][j  ].A*
                                      W[(search_directions)*scalar_dim + index(i+1,j  ,k)]+
                                      SolBlk_ptr->Grid.Cell[i  ][j+1].A*
                                      W[(search_directions)*scalar_dim + index(i  ,j+1,k)]+
                                      SolBlk_ptr->Grid.Cell[i+1][j+1].A*
                                      W[(search_directions)*scalar_dim + index(i+1,j+1,k)])/
                                     (SolBlk_ptr->Grid.Cell[i  ][j  ].A+
                                      SolBlk_ptr->Grid.Cell[i+1][j  ].A+
                                      SolBlk_ptr->Grid.Cell[i  ][j+1].A+
                                      SolBlk_ptr->Grid.Cell[i+1][j+1].A);
/*               buffer[buffer_count] = (SolBlk_ptr->Grid.Cell[i  ][j  ].A* */
/*                                       W[(search_directions)*scalar_dim + index(i  ,j  ,k)]+ */
/*                                       SolBlk_ptr->Grid.Cell[i+1][j  ].A* */
/*                                       W[(search_directions)*scalar_dim + index(i+1,j  ,k)]+ */
/*                                       SolBlk_ptr->Grid.Cell[i  ][j+1].A* */
/*                                       W[(search_directions)*scalar_dim + index(i  ,j+1,k)]+ */
/*                                       SolBlk_ptr->Grid.Cell[i+1][j+1].A* */
/*                                       W[(search_directions)*scalar_dim + index(i+1,j+1,k)]); */
	   } else {
              buffer[buffer_count] = (SolBlk_ptr->Grid.Cell[i  ][j  ].A*x[index(i  ,j  ,k)]+
                                      SolBlk_ptr->Grid.Cell[i+1][j  ].A*x[index(i+1,j  ,k)]+
                                      SolBlk_ptr->Grid.Cell[i  ][j+1].A*x[index(i  ,j+1,k)]+
                                      SolBlk_ptr->Grid.Cell[i+1][j+1].A*x[index(i+1,j+1,k)])/
                                     (SolBlk_ptr->Grid.Cell[i  ][j  ].A+
                                      SolBlk_ptr->Grid.Cell[i+1][j  ].A+
                                      SolBlk_ptr->Grid.Cell[i  ][j+1].A+
                                      SolBlk_ptr->Grid.Cell[i+1][j+1].A);
/*               buffer[buffer_count] = (SolBlk_ptr->Grid.Cell[i  ][j  ].A*x[index(i  ,j  ,k)]+ */
/*                                       SolBlk_ptr->Grid.Cell[i+1][j  ].A*x[index(i+1,j  ,k)]+ */
/*                                       SolBlk_ptr->Grid.Cell[i  ][j+1].A*x[index(i  ,j+1,k)]+ */
/*                                       SolBlk_ptr->Grid.Cell[i+1][j+1].A*x[index(i+1,j+1,k)]); */
	   } /* endif */
        } /* endfor */
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * GMRES_Block::LoadSendBuffer_C2F -- Loads send message buffer for            *
 *                                    coarse to fine block message             *
 *                                    passing.                                 *
 *******************************************************************************/
inline int GMRES_Block::LoadSendBuffer_C2F(double *buffer,
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
  Euler2D_pState Wcoarse, Wfine;

  if (j_inc > 0) {
    if (i_inc > 0) {
      for ( j = j_min; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max); j += j_inc) {
	for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
	  // Perform limited linear least squares reconstruction in cell (i, j_min).
	  SubcellReconstruction(i, j, LIMITER_ZERO);
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
	    for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
	      if (vector_switch) {
		Wcoarse[k] = W[(search_directions)*scalar_dim + index(i,j_min,k-1)];
	      } else {
		Wcoarse[k] = x[index(i,j_min,k-1)];
	      } /* endif */
	    } /* endfor */
	    dX = (SolBlk_ptr->Grid.Node[i][j_min].X+
		  HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
		  HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
		  SolBlk_ptr->Grid.Cell[i][j_min].Xc)/FOUR -
                 SolBlk_ptr->Grid.Cell[i][j_min].Xc;
	    Wfine = Wcoarse +
                    (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdx[i][j_min])*dX.x +
                    (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdy[i][j_min])*dX.y;
	    for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
	      buffer_count = buffer_count + 1;
	      if (buffer_count >= buffer_size) return(1);
	      buffer[buffer_count] = Wfine[k];
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
	    dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
		  SolBlk_ptr->Grid.Node[i+1][j_min].X+
		  SolBlk_ptr->Grid.Cell[i][j_min].Xc+
		  HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                 SolBlk_ptr->Grid.Cell[i][j_min].Xc;
	    Wfine = Wcoarse +
                    (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdx[i][j_min])*dX.x +
                    (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdy[i][j_min])*dX.y;
	    for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
	      buffer_count = buffer_count + 1;
	      if (buffer_count >= buffer_size) return(1);
	      buffer[buffer_count] = Wfine[k];
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
	    for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
	      if (vector_switch) {
		Wcoarse[k] = W[(search_directions)*scalar_dim + index(i,j_min,k-1)];
	      } else {
		Wcoarse[k] = x[index(i,j_min,k-1)];
	      } /* endif */
	    } /* endfor */
	    dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
		  SolBlk_ptr->Grid.Cell[i][j_min].Xc+
		  SolBlk_ptr->Grid.Node[i][j_min+1].X+
		  HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                 SolBlk_ptr->Grid.Cell[i][j_min].Xc;
	    Wfine = Wcoarse +
                    (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdx[i][j_min])*dX.x +
                    (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdy[i][j_min])*dX.y;
	    for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
	      buffer_count = buffer_count + 1;
	      if (buffer_count >= buffer_size) return(1);
	      buffer[buffer_count] = Wfine[k];
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
	    dX = (SolBlk_ptr->Grid.Cell[i][j_min].Xc+
		  HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
		  HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
		  SolBlk_ptr->Grid.Node[i+1][j_min+1].X)/FOUR -
                 SolBlk_ptr->Grid.Cell[i][j_min].Xc;
	    Wfine = Wcoarse +
                    (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdx[i][j_min])*dX.x +
                    (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdy[i][j_min])*dX.y;
	    for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
	      buffer_count = buffer_count + 1;
	      if (buffer_count >= buffer_size) return(1);
	      buffer[buffer_count] = Wfine[k];
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
              SubcellReconstruction(i, j_min, LIMITER_ZERO);
              // Evaluate SE sub (fine) cell values.
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i,j_min,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i,j_min,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
                    SolBlk_ptr->Grid.Node[i+1][j_min].X+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdy[i][j_min])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Node[i][j_min].X+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate NE sub (fine) cell values.
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i,j_min,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i,j_min,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
                    SolBlk_ptr->Grid.Node[i+1][j_min+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdy[i][j_min])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    SolBlk_ptr->Grid.Node[i][j_min+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdy[i][j_min])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } else {
        if (i_inc > 0) {
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Perform limited linear least squares reconstruction in cell (i, j_min).
              SubcellReconstruction(i, j_min, LIMITER_ZERO);
              // Evaluate NW sub (fine) cell values.
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i,j_min,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i,j_min,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    SolBlk_ptr->Grid.Node[i][j_min+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdy[i][j_min])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
                    SolBlk_ptr->Grid.Node[i+1][j_min+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdy[i][j_min])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate SW sub (fine) cell values.
             for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i,j_min,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i,j_min,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (SolBlk_ptr->Grid.Node[i][j_min].X+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
                    SolBlk_ptr->Grid.Node[i+1][j_min].X+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdy[i][j_min])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
        } else {
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Perform limited linear least squares reconstruction in cell (i, j_min).
              SubcellReconstruction(i, j_min, LIMITER_ZERO);
              // Evaluate NE sub (fine) cell values.
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i,j_min,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i,j_min,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
                    SolBlk_ptr->Grid.Node[i+1][j_min+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdy[i][j_min])*dX.y;
              Wfine = Wfine.U();
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    SolBlk_ptr->Grid.Node[i][j_min+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdy[i][j_min])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate SE sub (fine) cell values.
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i,j_min,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i,j_min,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
                    SolBlk_ptr->Grid.Node[i+1][j_min].X+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdy[i][j_min])*dX.y;
              for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Node[i][j_min].X+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdx[i][j_min])*dX.x +
                      (SolBlk_ptr->phi[i][j_min]^SolBlk_ptr->dWdy[i][j_min])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
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
              SubcellReconstruction(i_min, j, LIMITER_ZERO);
              // Evaluate SW sub (fine) cell values.
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i_min,j,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i_min,j,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (SolBlk_ptr->Grid.Node[i_min][j].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdy[i_min][j])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j].X+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdy[i_min][j])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    SolBlk_ptr->Grid.Node[i_min][j+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdy[i_min][j])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdy[i_min][j])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
        } else {
           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
              // Perform limited linear least squares reconstruction in cell (i_min, j).
              SubcellReconstruction(i_min, j, LIMITER_ZERO);
              // Evaluate SE sub (fine) cell values.
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i_min,j,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i_min,j,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j].X+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdy[i_min][j])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Node[i_min][j].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdy[i_min][j])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdy[i_min][j])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    SolBlk_ptr->Grid.Node[i_min][j+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdy[i_min][j])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } else {
        if (i_inc > 0) {
           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
              // Perform limited linear least squares reconstruction in cell (i_min, j).
              SubcellReconstruction(i_min, j, LIMITER_ZERO);
              // Evaluate NW sub (fine) cell values.
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i_min,j,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i_min,j,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    SolBlk_ptr->Grid.Node[i_min][j+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdy[i_min][j])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate NE sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdy[i_min][j])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Node[i_min][j].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdy[i_min][j])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j].X+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdy[i_min][j])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
        } else {
           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
              // Perform limited linear least squares reconstruction in cell (i_min, j).
              SubcellReconstruction(i_min, j, LIMITER_ZERO);
              // Evaluate NE sub (fine) cell values.
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
   	         if (vector_switch) {
	            Wcoarse[k] = W[(search_directions)*scalar_dim + index(i_min,j,k-1)];
	         } else {
	            Wcoarse[k] = x[index(i_min,j,k-1)];
	         } /* endif */
	      } /* endfor */
              dX = (SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdy[i_min][j])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    SolBlk_ptr->Grid.Node[i_min][j+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdy[i_min][j])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
                 buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j].X+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdy[i_min][j])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
              // Evaluate SW sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Node[i_min][j].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              Wfine = Wcoarse +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdx[i_min][j])*dX.x +
                      (SolBlk_ptr->phi[i_min][j]^SolBlk_ptr->dWdy[i_min][j])*dX.y;
              for (k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
  	         buffer_count = buffer_count + 1;
                 if (buffer_count >= buffer_size) return(1);
                 buffer[buffer_count] = Wfine[k];
              } /* endfor */
           } /* endfor */
        } /* endif */
     } /* endif */
  } /* endif */
  return(0);
}

/*******************************************************************************
 * GMRES_Block::UnloadReceiveBuffer -- Unloads receive message buffer.         *
 *******************************************************************************/
inline int GMRES_Block::UnloadReceiveBuffer(double *buffer,
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
        for ( k = 0 ; k < NUM_VAR_EULER2D; ++ k) {
  	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
	   if (vector_switch) {
	      W[(search_directions)*scalar_dim + index(i,j,k)] = buffer[buffer_count];
	   } else {
	      x[index(i,j,k)] = buffer[buffer_count];
	   } /* endif */
        } /* endfor */
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * GMRES_Block::UnloadReceiveBuffer_F2C -- Unloads receive message             *
 *                                         buffer for fine to coarse           *
 *                                         block message passing.              *
 *******************************************************************************/
inline int GMRES_Block::UnloadReceiveBuffer_F2C(double *buffer,
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
        for ( k = 0 ; k < NUM_VAR_EULER2D; ++ k) {
           buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
	   if (vector_switch) {
	      W[(search_directions)*scalar_dim + index(i,j,k)] = 
                 buffer[buffer_count];
/* 	      W[(search_directions)*scalar_dim + index(i,j,k)] =  */
/*                  buffer[buffer_count]/SolBlk_ptr->Grid.Cell[i][j].A; */
	   } else {
	      x[index(i,j,k)] = buffer[buffer_count];
/* 	      x[index(i,j,k)] = buffer[buffer_count]/SolBlk_ptr->Grid.Cell[i][j].A; */
	   } /* endif */
        } /* endfor */
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * GMRES_Block::UnloadReceiveBuffer_C2F -- Unloads receive message             *
 *                                         buffer for coarse to fine           *
 *                                         block message passing.              *
 *******************************************************************************/
inline int GMRES_Block::UnloadReceiveBuffer_C2F(double *buffer,
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
        for ( k = 0 ; k < NUM_VAR_EULER2D; ++ k) {
    	   buffer_count = buffer_count + 1;
           if (buffer_count >= buffer_size) return(1);
	   if (vector_switch) {
	      W[(search_directions)*scalar_dim + index(i,j,k)] = buffer[buffer_count];
	   } else {
	      x[index(i,j,k)] = buffer[buffer_count];
	   } /* endif */
        } /* endfor */
     } /* endfor */
  } /* endfor */
  return(0);
}

/**************************************************************************
 * GMRES_Block::SubcellReconstruction --                                  *
 *              Performs the subcell reconstruction of solution state     *
 *              within a given cell (i,j) of the computational mesh for   *
 *              the specified quadrilateral solution block.               *
 **************************************************************************/
inline void GMRES_Block::SubcellReconstruction(const int i, 
					       const int j,
					       const int Limiter) {
    
  int n, n2, n_pts, i_index[8], j_index[8], k;
  double u0, u0Min, u0Max, uQuad[4], phi_n;
  double DxDx_ave, DxDy_ave, DyDy_ave;
  Vector2D dX;
  Euler2D_pState U0, DU, DUDx_ave, DUDy_ave;

  /* Carry out the limited solution reconstruction in
     each cell of the computational mesh. */

  if (i == ICl-2 || i == ICu+2 ||
      j == JCl-2 || j == JCu+2) {
    n_pts = 0;
  } else if ((i == ICl-1) && 
             (SolBlk_ptr->Grid.BCtypeW[j] != BC_NONE)) {
    if (j == JCl-1 || j == JCu+1) {
       n_pts = 0;
    } else if (SolBlk_ptr->Grid.BCtypeW[j] == BC_PERIODIC ||
               SolBlk_ptr->Grid.BCtypeW[j] == BC_CONSTANT_EXTRAPOLATION ||
               SolBlk_ptr->Grid.BCtypeW[j] == BC_LINEAR_EXTRAPOLATION ||
               SolBlk_ptr->Grid.BCtypeW[j] == BC_CHARACTERISTIC) {
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
             (SolBlk_ptr->Grid.BCtypeE[j] != BC_NONE)) {
    if (j == JCl-1 || j == JCu+1) {
       n_pts = 0;
    } else if (SolBlk_ptr->Grid.BCtypeE[j] == BC_PERIODIC ||
               SolBlk_ptr->Grid.BCtypeE[j] == BC_CONSTANT_EXTRAPOLATION ||
               SolBlk_ptr->Grid.BCtypeE[j] == BC_LINEAR_EXTRAPOLATION ||
               SolBlk_ptr->Grid.BCtypeE[j] == BC_CHARACTERISTIC) {
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
             (SolBlk_ptr->Grid.BCtypeS[i] != BC_NONE)) {
    if (i == ICl-1 || i == ICu+1) {
       n_pts = 0;
    } else if (SolBlk_ptr->Grid.BCtypeS[i] == BC_PERIODIC ||
               SolBlk_ptr->Grid.BCtypeS[i] == BC_CONSTANT_EXTRAPOLATION ||
               SolBlk_ptr->Grid.BCtypeS[i] == BC_LINEAR_EXTRAPOLATION ||
               SolBlk_ptr->Grid.BCtypeS[i] == BC_CHARACTERISTIC) {
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
             (SolBlk_ptr->Grid.BCtypeN[i] != BC_NONE)) {
    if (i == ICl-1 || i == ICu+1) {
       n_pts = 0;
    } else if (SolBlk_ptr->Grid.BCtypeN[i] == BC_PERIODIC ||
               SolBlk_ptr->Grid.BCtypeN[i] == BC_CONSTANT_EXTRAPOLATION ||
               SolBlk_ptr->Grid.BCtypeN[i] == BC_LINEAR_EXTRAPOLATION ||
               SolBlk_ptr->Grid.BCtypeN[i] == BC_CHARACTERISTIC) {
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
      DUDx_ave = Euler2D_W_VACUUM;
      DUDy_ave = Euler2D_W_VACUUM;
      DxDx_ave = ZERO;
      DxDy_ave = ZERO;
      DyDy_ave = ZERO;

      for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
         if (vector_switch) {
            U0[k] = W[(search_directions)*scalar_dim + index(i,j,k-1)];
         } else {
            U0[k] = x[index(i,j,k-1)];
         } /* endif */
      } /* endfor */

      for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
          dX = SolBlk_ptr->Grid.Cell[ i_index[n2] ][ j_index[n2] ].Xc - 
               SolBlk_ptr->Grid.Cell[i][j].Xc;
          for ( k = 1 ; k <= NUM_VAR_EULER2D; ++ k) {
             if (vector_switch) {
                DU[k] = W[(search_directions)*scalar_dim + index(i_index[n2] , j_index[n2] , k-1)] - 
                        U0[k];
             } else {
                DU[k] = x[index( i_index[n2] , j_index[n2] , k-1)] - 
                        U0[k];
             } /* endif */
          } /* endfor */
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
      SolBlk_ptr->dWdx[i][j] = (DUDx_ave*DyDy_ave-DUDy_ave*DxDy_ave)/
                               (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
      SolBlk_ptr->dWdy[i][j] = (DUDy_ave*DxDx_ave-DUDx_ave*DxDy_ave)/
                               (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
  
      // Calculate slope limiters. 
      if (!SolBlk_ptr->Freeze_Limiter) {
         for ( n = 1 ; n <= NUM_VAR_EULER2D ; ++n ) {
	    u0 = U0[n];
            u0Min = U0[n];
            u0Max = u0Min;
            for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
               if (vector_switch) {
                  u0Min = min(u0Min, W[(search_directions)*scalar_dim + index(i_index[n2] , j_index[n2] ,n-1)]);
                  u0Max = max(u0Max, W[(search_directions)*scalar_dim + index(i_index[n2] , j_index[n2] ,n-1)]);
               } else {
                  u0Min = min(u0Min, x[index(i_index[n2] , j_index[n2] , n-1)]);
                  u0Max = max(u0Max, x[index(i_index[n2] , j_index[n2] , n-1)]);
               } /* endif */
            } /* endfor */
    
            dX = SolBlk_ptr->Grid.xfaceE(i, j)-SolBlk_ptr->Grid.Cell[i][j].Xc;
            uQuad[0] = u0 + 
                       SolBlk_ptr->dWdx[i][j][n]*dX.x +
                       SolBlk_ptr->dWdy[i][j][n]*dX.y ;
            dX = SolBlk_ptr->Grid.xfaceW(i, j)-SolBlk_ptr->Grid.Cell[i][j].Xc;
            uQuad[1] = u0 + 
                       SolBlk_ptr->dWdx[i][j][n]*dX.x +
                       SolBlk_ptr->dWdy[i][j][n]*dX.y ;
            dX = SolBlk_ptr->Grid.xfaceN(i, j)-SolBlk_ptr->Grid.Cell[i][j].Xc;
            uQuad[2] = u0 + 
                       SolBlk_ptr->dWdx[i][j][n]*dX.x +
                       SolBlk_ptr->dWdy[i][j][n]*dX.y ;
            dX = SolBlk_ptr->Grid.xfaceS(i, j)-SolBlk_ptr->Grid.Cell[i][j].Xc;
            uQuad[3] = u0 + 
                       SolBlk_ptr->dWdx[i][j][n]*dX.x +
                       SolBlk_ptr->dWdy[i][j][n]*dX.y ;
    
            switch(Limiter) {
              case LIMITER_ONE :
                phi_n = ONE;
                break;
              case LIMITER_ZERO :
                phi_n = ZERO;
                break;
              case LIMITER_BARTH_JESPERSEN :
                phi_n = Limiter_BarthJespersen(uQuad, u0, 
                                               u0Min, u0Max, 4);
                break;
              case LIMITER_VENKATAKRISHNAN :
                phi_n = Limiter_Venkatakrishnan(uQuad, u0, 
                                                u0Min, u0Max, 4);
                break;
              case LIMITER_VANLEER :
                phi_n = Limiter_VanLeer(uQuad, u0, 
                                        u0Min, u0Max, 4);
                break;
              case LIMITER_VANALBADA :
                phi_n = Limiter_VanAlbada(uQuad, u0, 
                                          u0Min, u0Max, 4);
                break;
              default:
                phi_n = Limiter_BarthJespersen(uQuad, u0, 
                                               u0Min, u0Max, 4);
                break;
            } /* endswitch */

	    SolBlk_ptr->phi[i][j][n] = phi_n;

         } /* endfor */
      } /* endif */
  } else {
      SolBlk_ptr->dWdx[i][j] = Euler2D_W_VACUUM;
      SolBlk_ptr->dWdy[i][j] = Euler2D_W_VACUUM; 
      SolBlk_ptr->phi[i][j]  = Euler2D_W_VACUUM;
  } /* endif */

}

/*******************************************************************************
 * GMRES_Block::LoadSendBuffer_Flux_F2C -- Loads send message buffer for       *
 *                                         fine to coarse block message        *
 *                                         passing of conservative             *
 *                                         solution fluxes.                    *
 *******************************************************************************/
inline int GMRES_Block::LoadSendBuffer_Flux_F2C(double *buffer,
						int &buffer_count,
						const int buffer_size,
						const int i_min, 
						const int i_max,
						const int i_inc,
						const int j_min, 
						const int j_max,
						const int j_inc) {
  return(0);
}

/*******************************************************************************
 * GMRES_Block::UnloadReceiveBuffer_Flux_F2C -- Unloads receive message        *
 *                                              buffer for fine to coarse      *
 *                                              block message passing of       *
 *                                              conservative solution fluxes.  *
 *******************************************************************************/
inline int GMRES_Block::UnloadReceiveBuffer_Flux_F2C(double *buffer,
						     int &buffer_count,
						     const int buffer_size,
						     const int i_min, 
						     const int i_max,
						     const int i_inc,
						     const int j_min, 
						     const int j_max,
						     const int j_inc) {
  return(0);
}

/***********************************************************
 * Class: GMRES_Block_RightPrecon_MatrixFree               *
 *                                                         *
 * Member functions                                        *
 *              solve -- Solve linear system using GMRES   *
 *                       algorithm                         *
 *           get_iter -- Return total number of iterations *
 *                       performed in the GMRES algorithm  *
 * get_rel_resid_norm -- Return final relative residual    *
 *                                                         *
 ***********************************************************/
class GMRES_Block_RightPrecon_MatrixFree
{
private:
    int iter;
    double rel_resid;

public:
    int dim;
    int max_iter;
    double tol;

    GMRES_Block_RightPrecon_MatrixFree(int dim, int max_iter, double tol);
    ~GMRES_Block_RightPrecon_MatrixFree() {}

    void solve(Euler2D_Quad_Block *Soln_ptr, 
	       AdaptiveBlock2D_List &Soln_Block_List, 
               Euler2D_Input_Parameters &Input_Parameters,
	       GMRES_Block *G,
	       BpPrecon *M,
	       const int &normalization,
               const int Number_of_Newton_Steps);
    
    int    get_iter() {return iter;}
    double get_rel_resid_norm() {return rel_resid;}
    
};

/********************************************************
 * Euler2DQuad_GMRES -- External subroutines.           *
 ********************************************************/

extern void GMRES_Algorithm(Euler2D_Quad_Block *Soln_ptr, 
		            AdaptiveBlock2D_List &Soln_Block_List,
                            Euler2D_Input_Parameters &Input_Parameters,
		            GMRES_Block *G, double &gmrestol, 
		            BILUK *MBILUK, BJacobi *MBJacobi,
		            const int normalization,
                            const int Number_of_Newton_Steps);

extern void GMRES_BCs(Euler2D_Quad_Block &SolnBlk, 
		      double *v,
		      GMRES_Block &G,
		      int normalization);

extern ColumnVector GMRES_Reflect(const ColumnVector cv, 
				  const Vector2D &norm_dir,
				  int blocksize);

extern void normalize_RHS(double * RHS, 
			   Euler2D_Quad_Block &SolnBlk);

extern void denormalize_solution(double *RHS, 
			          Euler2D_Quad_Block &SolnBlk);

extern void normalize_dWdU(DenseMatrix &dWdU);

extern void normalize_dUdW(DenseMatrix &dUdW);

#endif // _EULER2D_QUAD_GMRES_INCLUDED
