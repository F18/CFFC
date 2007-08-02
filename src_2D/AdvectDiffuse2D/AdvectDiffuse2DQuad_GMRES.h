#ifndef _ADVECTDIFFUSE2D_QUAD_GMRES_INCLUDED
#define _ADVECTDIFFUSE2D_QUAD_GMRES_INCLUDED

/* Include header file. */
#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _MV_VECTOR_ALL_H_
#include "mvv.h" 
#endif // _MV_VECTOR_ALL_H_

#ifndef _MV_MATRIX_ALL_H_
#include "mvm.h"
#endif // _MV_MATRIX_ALL_H_

#ifndef _MVMTP_H
#include "mvmtp.h"
#endif // _MVMTP_H

#ifndef _MV_BLAS_H
#include "mvblas.h" 
#endif // _MV_BLAS_H

#ifndef _COMPROW_DOUBLE_H
#include "comprow_double.h"
#endif // _COMPROW_DOUBLE_H

#ifndef ILUPRE_H
#include "ilupre_double.h"
#endif // ILUPRE_H

#ifndef DIAGPRE_H
#include "diagpre_double.h"
#endif // DIAGPRE_H

/***********************************************************
 * Class: GMRES                                            *
 *                                                         *
 * Member functions                                        *
 *       s      -- Return ..........                       *
 *      cs      -- Return cos vector                       *
 *      sn      -- Return sin vector                       *
 *       w      -- Return Az vector                        *
 *       z      -- Return inversion of preconditioner      *
 *                  times v vector                         *
 *       r      -- Return residual vector                  *
 *       b      -- Return RHS vector                       *
 *      Ax      -- Return LHS vector                       *
 *      Q1      -- Return matrix free vector               *
 *      Q2      -- Return matrix free vector               *
 *       x      -- Return solution vector, delta u         *
 *                 (Note: it also contains initial guess   * 
 *                        of zeros.)                       *
 *       v      -- Return Krylov search direction vector   *
 *       H      -- Return upper Hessenberg H matrix        *
 *       A      -- Reutrn Jacobian matrix                  *
 *      BM      -- Return block ILU perconditioner matrix  *
 *      BD      -- Return block diagonal perconditioner    *
 *                 matrix                                  *
 *       m      -- Return restart                          *
 *    xpts      -- Return number of nodes in x-direction   *
 *    ypts      -- Return number of nodes in y-direction   *  
 * overlap      -- Return level of overlap                 *
 * P_Switch     -- Return preconditioner switch            *
 * max_gmres_iter -- Return maximun number of gmres        *
 *                   iterations                            * 
 *  copy_residual -- Copy residual from SolnBlk to an      *
 *                   1D array(VECTOR_double)               *
 *                                                         *
 * ** MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING **     *
 * vector_switch -- Return vector switch                   *
 *                  Note: 0 for x vector (solution)        *
 *                        1 for z vector (z = Minv * v)    *
 *     NCi      -- Return number of cells in               *
 *                 the i-direction (zeta-direction).       *
 *     ICl      -- Return lower index for cells in         *
 *                 the i-direction (zeta-direction).       *
 *     ICu      -- Return upper index for cells in         *
 *                 the i-direction (zeta-direction).       *
 *     NCj      -- Return number of cells in               *
 *                 the j-direction (eta-direction).        *
 *     JCl      -- Return lower index for cells in         *
 *                 the j-direction (eta-direction).        *
 *     JCu      -- Return upper index for cells in         *
 *                 the j-direction (eta-direction).        *
 *  Nghost      -- Return number of ghost (halo or         *
 *                 overlap) cells.                         *
 *    Grid      -- Return a dummy object.                  *
 *                 (Required by message passing routines)  *  
 *  SolBlk_ptr  -- Return pointer to solution block.       *
 *                                                         *
 * Member operators                                        *
 *      G -- a GMRES variable object                       *
 *                                                         *
 * G = G;                                                  *
 * cout << G; (output function)                            *
 * cin  >> G; (input function)                             *
 *                                                         *
 ***********************************************************/
class GMRES{
private:
public:
  VECTOR_double                      s; //
  VECTOR_double                     cs; // cos vector
  VECTOR_double                     sn; // sin vector
  VECTOR_double                      w; // Az
  VECTOR_double                      z; // inversion of preconditioner times v
  VECTOR_double                      r; // residual vector
  VECTOR_double                     Ax; // LHS vector
  VECTOR_double                      b; // RHS vector
  VECTOR_double                     Q1; // matrix free vector
  VECTOR_double                     Q2; // matrix free vector
  VECTOR_double                      x; // initial guess of delta u
  VECTOR_double                     *v; // Krylov search direction vector
  MATRIX_double                      H; // storage for upper Hessenberg H
  CompRow_Mat_double                 A; // Jacobian matrix
  CompRow_ILUPreconditioner_double  BM; // ILU perconditioner matrix
  DiagPreconditioner_double         BD; // diagonal perconditioner
  int                                m; // restart
  int                             xpts; // number of nodes in x-direction
  int                             ypts; // number of nodes in y-direction
  int                          overlap; // level of overlap
  int                         P_Switch; // preconditioner switch
  int                   max_gmres_iter; // maximun number of gmres 
                                        // iterations
  double                         dTime; // finite time step

  /* MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING. */
  int                      vector_switch; // to select the specified vector 
                                          // for message passing
  int                        NCi,ICl,ICu; // i-direction cell counters.
  int                        NCj,JCl,JCu; // j-direction cell counters.
  int                             Nghost; // Number of ghost cells.
  Grid2D_Quad_Block                 Grid; // dummy Grid variable.
  AdvectDiffuse2D_Quad_Block *SolBlk_ptr; // Pointer to solution block.   
                                          // Made public so can access them.
  
  /* Creation, copy, and assignment constructors. */
  GMRES(void) {
    m              = 0; 
    xpts           = 0; 
    ypts           = 0;
    overlap        = 0;
    P_Switch       = 1;
    max_gmres_iter = 0;
    dTime          = 0.0;

    vector_switch = 99;
    NCi = 0; 
    ICl = 0; ICu = 0; 
    NCj = 0; 
    JCl = 0; JCu = 0; Nghost = 0;
    SolBlk_ptr = NULL;
  }
  
  GMRES(const GMRES &G) {
    s  = G.s ; cs = G.cs; sn = G.sn;
    w  = G.w ;  z = G.z ;  r = G.r ;
    Ax = G.Ax;  b = G.b ;  v = G.v ;
    Q1 = G.Q1; Q2 = G.Q2;  x = G.x ;
    H  = G.H ;
    dTime = G.dTime;

    A  = G.A;
    BM = G.BM;
    BD = G.BD;

    m              = G.m;  
    xpts           = G.xpts;
    ypts           = G.ypts;
    overlap        = G.overlap;
    P_Switch       = G.P_Switch;
    max_gmres_iter = G.max_gmres_iter;

    vector_switch = G.vector_switch;
    NCi = G.NCi; ICl = G.ICl; ICu = G.ICu; 
    NCj = G.NCj; JCl = G.JCl; JCu = G.JCu; Nghost = G.Nghost;
    Grid = G.Grid; SolBlk_ptr = G.SolBlk_ptr;
  }
  
  /* Destructor. */
  // ~GMRES(void);
  // Use automatically generated destructor.
  
  /* Assignment operator. */
  // GMRES operator = (const GMRES &G);
  // Use automatically generated assignment operator.
  
  /* Allocate memory for GMRES variables. */
  void allocate(const int restart, int xpts, int ypts);
  
  /* Deallocate memory for GMRES variables. */
  void deallocate(void);
      
  /* Size of domain with ghost cells. */
  int Xcount_gc();
  int Ycount_gc();

  /* Size of domain with ghost cells. */
  int index(int i, int j);

  void copy_residual(VECTOR_double &RHS, AdvectDiffuse2D_Quad_Block SolnBlk);
  
  /* Input-output operators. */
  friend ostream &operator << (ostream &out_file, const GMRES &G);
  friend istream &operator >> (istream &in_file, GMRES &G);

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
 * GMRES::allocate -- Allocate memory.                                    *
 **************************************************************************/
inline void GMRES::allocate(const int restart, int xpts, int ypts) {
  assert(restart > 1);
  assert(xpts > 1);
  assert(ypts > 1);

  // Allocate memory
   s.newsize(restart+1);
  cs.newsize(restart+1);
  sn.newsize(restart+1);

   w.newsize(xpts*ypts);
   z.newsize(xpts*ypts);
   r.newsize(xpts*ypts);
   b.newsize(xpts*ypts);
  Q1.newsize(xpts*ypts);
  Q2.newsize(xpts*ypts);
  Ax.newsize(xpts*ypts);
  
   x.newsize(xpts*ypts); 
   H.newsize(restart+1, restart);

  v = new VECTOR_double[restart+1];
  for (int count = 0; count < restart+1; count++) {
    v[count].newsize(xpts*ypts);
  }

  // Initialization
  s = 0.0;
  cs = 0.0;
  sn = 0.0;

  w = 0.0;
  z = 0.0;
  r = 0.0;
  b = 0.0;
  Q1 = 0.0;
  Q2 = 0.0;
  Ax = 0.0;

  x = 0.0;
  H = 0.0;
  for (int count = 0; count < restart+1; count++) {
    v[count] = 0.0;
  }

}

/**************************************************************************
 * GMRES::deallocate -- Deallocate memory.                                *
 **************************************************************************/
inline void GMRES::deallocate(void) {

  delete []  v;
}

/**************************************************************************
 * GMRES::Xcount_gc -- Return number of nodes in x-dir with               *
 *                     ghost cells in both sides.                         *
 **************************************************************************/
inline int GMRES::Xcount_gc() {

  return (xpts+Nghost*2);
}

/**************************************************************************
 * GMRES::Ycount_gc -- Return number of nodes in y-dir with               *
 *                     ghost cells in both sides.                         *
 **************************************************************************/
inline int GMRES::Ycount_gc() {

  return (ypts+Nghost*2);
}

/**************************************************************************
 * GMRES::index -- Return index from zero to xpts*ypts-1                  *
 *                 given i and j positions in the domain                  *
 **************************************************************************/
inline int GMRES::index(int i, int j) {

  return (j*xpts+i);
}

/**************************************************************************
 * GMRES::copy_residual -- Copy residual from SolnBlk to an 1D array      *
 **************************************************************************/
inline void GMRES::copy_residual(VECTOR_double &RHS, AdvectDiffuse2D_Quad_Block SolnBlk)
{
  for (int j = SolnBlk.JCl-SolnBlk.Nghost; j <= SolnBlk.JCu+SolnBlk.Nghost; j++) {
    for (int i = SolnBlk.ICl-SolnBlk.Nghost; i <= SolnBlk.ICu+SolnBlk.Nghost; i++) {
      RHS[j*xpts+i] = SolnBlk.dudt[i][j][0];
    } /* endfor */
  } /* endfor */
}

/**************************************************************************
 * GMRES -- Input-output operators.                                       *
 **************************************************************************/
inline ostream &operator << (ostream &out_file, const GMRES &G) {

  int i, j; 
  out_file << G.m    << "\n";
  out_file << G.xpts << "\n";
  out_file << G.ypts << "\n";
  out_file << G.overlap << "\n";
  out_file << G.P_Switch << "\n";
  
  for ( j  = 0 ; j <= G.m ; ++j ) {
    out_file << G.s[j] << "\n";
  } /* endfor */
  for ( j  = 0 ; j <= G.m ; ++j ) {
    out_file << G.cs[j] << "\n";
  } /* endfor */
  for ( j  = 0 ; j <= G.m ; ++j ) {
    out_file << G.sn[j] << "\n";
  } /* endfor */
  for ( j  = 0 ; j <= G.xpts*G.ypts ; ++j ) {
    out_file << G.w[j] << "\n";
  } /* endfor */
  for ( j  = 0 ; j <= G.xpts*G.ypts ; ++j ) {
    out_file << G.z[j] << "\n";
  } /* endfor */
  for ( j  = 0 ; j <= G.xpts*G.ypts ; ++j ) {
    out_file << G.r[j] << "\n";
  } /* endfor */
  for ( j  = 0 ; j <= G.xpts*G.ypts ; ++j ) {
    out_file << G.Ax[j] << "\n";
  } /* endfor */
  for ( j  = 0 ; j <= G.xpts*G.ypts ; ++j ) {
    out_file << G.b[j] << "\n";
  } /* endfor */
  for ( j  = 0 ; j <= G.xpts*G.ypts ; ++j ) {
    out_file << G.Q1[j] << "\n";
  } /* endfor */
  for ( j  = 0 ; j <= G.xpts*G.ypts ; ++j ) {
    out_file << G.Q2[j] << "\n";
  } /* endfor */
  for ( i = 0 ; i <= G.m ; i++) {
    for ( j  = 0 ; j <= G.xpts*G.ypts ; ++j ) {
      out_file << "[" << i << "][" << j << "]" << G.v[i][j] << "\n";
    } /* endfor */ 
  } /* endfor */  
  return (out_file);
}

inline istream &operator>> (istream &in_file, GMRES &G) {
  
  int i, j;

  in_file >> G.m ; 
  in_file >> G.xpts;
  in_file >> G.ypts;
  in_file >> G.overlap;
  in_file >> G.P_Switch;
 
  for ( j  = 0 ; j <= G.m ; ++j ) {
    in_file >> G.s[j] ;
  } /* endfor */
  for ( j  = 0 ; j <= G.m ; ++j ) {
    in_file >> G.cs[j];
  } /* endfor */
  for ( j  = 0 ; j <= G.m ; ++j ) {
    in_file >> G.sn[j];
  } /* endfor */
  
  for ( j  = 0 ; j <= G.xpts*G.ypts ; ++j ) {
    in_file >> G.w[j] ;
  } /* endfor */
  for ( j  = 0 ; j <= G.xpts*G.ypts ; ++j ) {
    in_file >> G.z[j] ;
  } /* endfor */
  for ( j  = 0 ; j <= G.xpts*G.ypts ; ++j ) {
    in_file >> G.r[j] ;
  } /* endfor */
  for ( j  = 0 ; j <= G.xpts*G.ypts ; ++j ) {
    in_file >> G.Ax[j] ;
  } /* endfor */
  for ( j  = 0 ; j <= G.xpts*G.ypts ; ++j ) {
    in_file >> G.b[j] ;
  } /* endfor */
  for ( j  = 0 ; j <= G.xpts*G.ypts ; ++j ) {
    in_file >> G.Q1[j] ;
  } /* endfor */
  for ( j  = 0 ; j <= G.xpts*G.ypts ; ++j ) {
    in_file >> G.Q2[j] ;
  } /* endfor */
  for ( i = 0 ; i <= G.m ; i++) {
    for ( j  = 0 ; j <= G.xpts*G.ypts ; ++j ) {
      in_file >> G.v[i][j] ;
    } /* endfor */
  } /* endfor */
  
  return (in_file);
}

/*******************************************************************************
 *                                                                             *
 * MEMBER FUNCTIONS REQUIRED FOR MESSAGE PASSING.                              *
 *                                                                             *
 *******************************************************************************/

/*******************************************************************************
 * GMRES::NumVar -- Returns number of state variables.                         *
 *******************************************************************************/
inline int GMRES::NumVar(void) {
  return (int(NUM_VAR_ADVECTDIFFUSE2D));
}

/*******************************************************************************
 * GMRES::LoadSendBuffer -- Loads send message buffer.                         *
 *******************************************************************************/
inline int GMRES::LoadSendBuffer(double *buffer,
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
	   if (vector_switch) {
	     buffer[buffer_count] = z[index(i,j)];
	   } else {
	     buffer[buffer_count] = x[index(i,j)];
	   } /* endif */
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * GMRES::LoadSendBuffer_F2C -- Loads send message buffer for                  *
 *                              fine to coarse block message                   *
 *                              passing.                                       *
 *******************************************************************************/
inline int GMRES::LoadSendBuffer_F2C(double *buffer,
				     int &buffer_count,
				     const int buffer_size,
				     const int i_min, 
				     const int i_max,
				     const int i_inc,
				     const int j_min, 
				     const int j_max,
				     const int j_inc) {
  int i, j;
  for ( j  = j_min ; ((j_inc+2)/4) ? (j < j_max):(j > j_max) ; j += j_inc ) {
     for ( i = i_min ;  ((i_inc+2)/4) ? (i < i_max):(i > i_max) ; i += i_inc ) {
  	buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
        if (buffer_count >= buffer_size) return(1);
	if (vector_switch) {
           buffer[buffer_count] = (SolBlk_ptr->Grid.Cell[i  ][j  ].A*z[index(i  ,j  )]+
                                   SolBlk_ptr->Grid.Cell[i+1][j  ].A*z[index(i+1,j  )]+
                                   SolBlk_ptr->Grid.Cell[i  ][j+1].A*z[index(i  ,j+1)]+
                                   SolBlk_ptr->Grid.Cell[i+1][j+1].A*z[index(i+1,j+1)])/
                                  (SolBlk_ptr->Grid.Cell[i  ][j  ].A+
                                   SolBlk_ptr->Grid.Cell[i+1][j  ].A+
                                   SolBlk_ptr->Grid.Cell[i  ][j+1].A+
                                   SolBlk_ptr->Grid.Cell[i+1][j+1].A);
/*            buffer[buffer_count] = (SolBlk_ptr->Grid.Cell[i  ][j  ].A*z[index(i  ,j  )]+ */
/*                                    SolBlk_ptr->Grid.Cell[i+1][j  ].A*z[index(i+1,j  )]+ */
/*                                    SolBlk_ptr->Grid.Cell[i  ][j+1].A*z[index(i  ,j+1)]+ */
/*                                    SolBlk_ptr->Grid.Cell[i+1][j+1].A*z[index(i+1,j+1)]); */
	} else {
           buffer[buffer_count] = (SolBlk_ptr->Grid.Cell[i  ][j  ].A*x[index(i  ,j  )]+
                                   SolBlk_ptr->Grid.Cell[i+1][j  ].A*x[index(i+1,j  )]+
                                   SolBlk_ptr->Grid.Cell[i  ][j+1].A*x[index(i  ,j+1)]+
                                   SolBlk_ptr->Grid.Cell[i+1][j+1].A*x[index(i+1,j+1)])/
                                  (SolBlk_ptr->Grid.Cell[i  ][j  ].A+
                                   SolBlk_ptr->Grid.Cell[i+1][j  ].A+
                                   SolBlk_ptr->Grid.Cell[i  ][j+1].A+
                                   SolBlk_ptr->Grid.Cell[i+1][j+1].A);
/*            buffer[buffer_count] = (SolBlk_ptr->Grid.Cell[i  ][j  ].A*x[index(i  ,j  )]+ */
/*                                    SolBlk_ptr->Grid.Cell[i+1][j  ].A*x[index(i+1,j  )]+ */
/*                                    SolBlk_ptr->Grid.Cell[i  ][j+1].A*x[index(i  ,j+1)]+ */
/*                                    SolBlk_ptr->Grid.Cell[i+1][j+1].A*x[index(i+1,j+1)]); */
	} /* endif */
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * GMRES::LoadSendBuffer_C2F -- Loads send message buffer for                  *
 *                              coarse to fine block message                   *
 *                              passing.                                       *
 *******************************************************************************/
inline int GMRES::LoadSendBuffer_C2F(double *buffer,
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
	    dX = SolBlk_ptr->Grid.centroidSW(i,j) - SolBlk_ptr->Grid.Cell[i][j].Xc;
            if (vector_switch) {
              ufine = z[index(i,j)] +
                      (SolBlk_ptr->phi[i][j]*SolBlk_ptr->dudx[i][j])*dX.x +
                      (SolBlk_ptr->phi[i][j]*SolBlk_ptr->dudy[i][j])*dX.y;
            } else {
              ufine = x[index(i,j)] +
                      (SolBlk_ptr->phi[i][j]*SolBlk_ptr->dudx[i][j])*dX.x +
                      (SolBlk_ptr->phi[i][j]*SolBlk_ptr->dudy[i][j])*dX.y;
            } /* endif */
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
	    dX = SolBlk_ptr->Grid.centroidSE(i,j) - SolBlk_ptr->Grid.Cell[i][j].Xc;
            if (vector_switch) {
              ufine = z[index(i,j)] +
                      (SolBlk_ptr->phi[i][j]*SolBlk_ptr->dudx[i][j])*dX.x +
                      (SolBlk_ptr->phi[i][j]*SolBlk_ptr->dudy[i][j])*dX.y;
            } else {
              ufine = x[index(i,j)] +
                      (SolBlk_ptr->phi[i][j]*SolBlk_ptr->dudx[i][j])*dX.x +
                      (SolBlk_ptr->phi[i][j]*SolBlk_ptr->dudy[i][j])*dX.y;
            } /* endif */
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
	    dX = SolBlk_ptr->Grid.centroidNW(i,j) - SolBlk_ptr->Grid.Cell[i][j].Xc;
            if (vector_switch) {
              ufine = z[index(i,j)] +
                      (SolBlk_ptr->phi[i][j]*SolBlk_ptr->dudx[i][j])*dX.x +
                      (SolBlk_ptr->phi[i][j]*SolBlk_ptr->dudy[i][j])*dX.y;
            } else {
              ufine = x[index(i,j)] +
                      (SolBlk_ptr->phi[i][j]*SolBlk_ptr->dudx[i][j])*dX.x +
                      (SolBlk_ptr->phi[i][j]*SolBlk_ptr->dudy[i][j])*dX.y;
            } /* endif */
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
	    dX = SolBlk_ptr->Grid.centroidNE(i,j) - SolBlk_ptr->Grid.Cell[i][j].Xc;
            if (vector_switch) {
              ufine = z[index(i,j)] +
                      (SolBlk_ptr->phi[i][j]*SolBlk_ptr->dudx[i][j])*dX.x +
                      (SolBlk_ptr->phi[i][j]*SolBlk_ptr->dudy[i][j])*dX.y;
            } else {
              ufine = x[index(i,j)] +
                      (SolBlk_ptr->phi[i][j]*SolBlk_ptr->dudx[i][j])*dX.x +
                      (SolBlk_ptr->phi[i][j]*SolBlk_ptr->dudy[i][j])*dX.y;
            } /* endif */
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
              dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
                    SolBlk_ptr->Grid.Node[i+1][j_min].X+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              if (vector_switch) {
                 ufine = z[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } /* endif */
  	      buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
              if (buffer_count >= buffer_size) return(1);
              buffer[buffer_count] = ufine;
              // Evaluate SW sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Node[i][j_min].X+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              if (vector_switch) {
                 ufine = z[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } /* endif */
  	      buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
              if (buffer_count >= buffer_size) return(1);
              buffer[buffer_count] = ufine;
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate NE sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
                    SolBlk_ptr->Grid.Node[i+1][j_min+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              if (vector_switch) {
                 ufine = z[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } /* endif */
  	      buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
              if (buffer_count >= buffer_size) return(1);
              buffer[buffer_count] = ufine;
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    SolBlk_ptr->Grid.Node[i][j_min+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              if (vector_switch) {
                 ufine = z[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } /* endif */
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
              dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    SolBlk_ptr->Grid.Node[i][j_min+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              if (vector_switch) {
                 ufine = z[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } /* endif */
              buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
              if (buffer_count >= buffer_size) return(1);
              buffer[buffer_count] = ufine;
              // Evaluate NE sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
                    SolBlk_ptr->Grid.Node[i+1][j_min+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              if (vector_switch) {
                 ufine = z[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } /* endif */
  	      buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
              if (buffer_count >= buffer_size) return(1);
              buffer[buffer_count] = ufine;
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate SW sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Node[i][j_min].X+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              if (vector_switch) {
                 ufine = z[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } /* endif */
  	      buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
              if (buffer_count >= buffer_size) return(1);
              buffer[buffer_count] = ufine;
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
                    SolBlk_ptr->Grid.Node[i+1][j_min].X+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              if (vector_switch) {
                 ufine = z[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } /* endif */
  	      buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
              if (buffer_count >= buffer_size) return(1);
              buffer[buffer_count] = ufine;
           } /* endfor */
        } else {
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Perform limited linear least squares reconstruction in cell (i, j_min).
              SubcellReconstruction(i, j_min, LIMITER_VENKATAKRISHNAN);
              // Evaluate NE sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X)+
                    SolBlk_ptr->Grid.Node[i+1][j_min+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              if (vector_switch) {
                 ufine = z[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } /* endif */
  	      buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
              if (buffer_count >= buffer_size) return(1);
              buffer[buffer_count] = ufine;
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    SolBlk_ptr->Grid.Node[i][j_min+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min+1].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
             if (vector_switch) {
                 ufine = z[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } /* endif */
              buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
              if (buffer_count >= buffer_size) return(1);
              buffer[buffer_count] = ufine;
           } /* endfor */
           for ( i = i_min ;  ((i_inc+1)/2) ? (i <= i_max):(i >= i_max) ; i += i_inc ) {
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
                    SolBlk_ptr->Grid.Node[i+1][j_min].X+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i+1][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              if (vector_switch) {
                 ufine = z[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } /* endif */
  	      buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
              if (buffer_count >= buffer_size) return(1);
              buffer[buffer_count] = ufine;
              // Evaluate SW sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Node[i][j_min].X+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i+1][j_min].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i][j_min].X+SolBlk_ptr->Grid.Node[i][j_min+1].X)+
                    SolBlk_ptr->Grid.Cell[i][j_min].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i][j_min].Xc;
              if (vector_switch) {
                 ufine = z[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudx[i][j_min])*dX.x +
                         (SolBlk_ptr->phi[i][j_min]*SolBlk_ptr->dudy[i][j_min])*dX.y;
              } /* endif */
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
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j].X+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              if (vector_switch) {
                 ufine = z[index(i_min,j)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } /* endif */
  	      buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
              if (buffer_count >= buffer_size) return(1);
              buffer[buffer_count] = ufine;
              // Evaluate SW sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Node[i_min][j].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              if (vector_switch) {
                 ufine = z[index(i_min,j)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } /* endif */
  	      buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
              if (buffer_count >= buffer_size) return(1);
              buffer[buffer_count] = ufine;
              // Evaluate NE sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              if (vector_switch) {
                 ufine = z[index(i_min,j)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } /* endif */
              buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
              if (buffer_count >= buffer_size) return(1);
              buffer[buffer_count] = ufine;
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    SolBlk_ptr->Grid.Node[i_min][j+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              if (vector_switch) {
                 ufine = z[index(i_min,j)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } /* endif */
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
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    SolBlk_ptr->Grid.Node[i_min][j+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              if (vector_switch) {
                 ufine = z[index(i_min,j)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } /* endif */
              buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
              if (buffer_count >= buffer_size) return(1);
              buffer[buffer_count] = ufine;
              // Evaluate NE sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              if (vector_switch) {
                 ufine = z[index(i_min,j)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } /* endif */
              buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
              if (buffer_count >= buffer_size) return(1);
              buffer[buffer_count] = ufine;
              // Evaluate SW sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Node[i_min][j].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              if (vector_switch) {
                 ufine = z[index(i_min,j)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } /* endif */
  	      buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
              if (buffer_count >= buffer_size) return(1);
              buffer[buffer_count] = ufine;
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j].X+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              if (vector_switch) {
                 ufine = z[index(i_min,j)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } /* endif */
  	      buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
              if (buffer_count >= buffer_size) return(1);
              buffer[buffer_count] = ufine;
           } /* endfor */
        } else {
           for ( j = j_min ; ((j_inc+1)/2) ? (j <= j_max):(j >= j_max) ; j += j_inc ) {
              // Perform limited linear least squares reconstruction in cell (i_min, j).
              SubcellReconstruction(i_min, j, LIMITER_VENKATAKRISHNAN);
              // Evaluate NE sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j+1].X)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              if (vector_switch) {
                 ufine = z[index(i_min,j)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } /* endif */
              buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
              if (buffer_count >= buffer_size) return(1);
              buffer[buffer_count] = ufine;
              // Evaluate NW sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    SolBlk_ptr->Grid.Node[i_min][j+1].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j+1].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              if (vector_switch) {
                 ufine = z[index(i_min,j)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } /* endif */
              buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
              if (buffer_count >= buffer_size) return(1);
              buffer[buffer_count] = ufine;
              // Evaluate SE sub (fine) cell values.
              dX = (HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    SolBlk_ptr->Grid.Node[i_min+1][j].X+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc+
                    HALF*(SolBlk_ptr->Grid.Node[i_min+1][j].X+SolBlk_ptr->Grid.Node[i_min+1][j+1].X))/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              if (vector_switch) {
                 ufine = z[index(i_min,j)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } /* endif */
  	      buffer_count = buffer_count + NUM_VAR_ADVECTDIFFUSE2D;
              if (buffer_count >= buffer_size) return(1);
              buffer[buffer_count] = ufine;
              // Evaluate SW sub (fine) cell values.
              dX = (SolBlk_ptr->Grid.Node[i_min][j].X+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min+1][j].X)+
                    HALF*(SolBlk_ptr->Grid.Node[i_min][j].X+SolBlk_ptr->Grid.Node[i_min][j+1].X)+
                    SolBlk_ptr->Grid.Cell[i_min][j].Xc)/FOUR -
                   SolBlk_ptr->Grid.Cell[i_min][j].Xc;
              if (vector_switch) {
                 ufine = z[index(i_min,j)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } else {
                 ufine = x[index(i,j_min)] +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudx[i_min][j])*dX.x +
                         (SolBlk_ptr->phi[i_min][j]*SolBlk_ptr->dudy[i_min][j])*dX.y;
              } /* endif */
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
 * GMRES::UnloadReceiveBuffer -- Unloads receive message buffer.               *
 *******************************************************************************/
inline int GMRES::UnloadReceiveBuffer(double *buffer,
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
	if (vector_switch) {
	   z[index(i,j)] = buffer[buffer_count];
	} else {
	   x[index(i,j)] = buffer[buffer_count];
	} /* endif */
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * GMRES::UnloadReceiveBuffer_F2C -- Unloads receive message                   *
 *                                   buffer for fine to coarse                 *
 *                                   block message passing.                    *
 *******************************************************************************/
inline int GMRES::UnloadReceiveBuffer_F2C(double *buffer,
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
	if (vector_switch) {
	   z[index(i,j)] = buffer[buffer_count];
/* 	   z[index(i,j)] = buffer[buffer_count]/SolBlk_ptr->Grid.Cell[i][j].A; */
	} else {
	   x[index(i,j)] = buffer[buffer_count];
/* 	   x[index(i,j)] = buffer[buffer_count]/SolBlk_ptr->Grid.Cell[i][j].A; */
	} /* endif */
     } /* endfor */
  } /* endfor */
  return(0);
}

/*******************************************************************************
 * GMRES::UnloadReceiveBuffer_C2F -- Unloads receive message                   *
 *                                   buffer for coarse to fine                 *
 *                                   block message passing.                    *
 *******************************************************************************/
inline int GMRES::UnloadReceiveBuffer_C2F(double *buffer,
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
	if (vector_switch) {
	   z[index(i,j)] = buffer[buffer_count];
	} else {
	   x[index(i,j)] = buffer[buffer_count];
	} /* endif */
     } /* endfor */
  } /* endfor */
  return(0);
}

/**************************************************************************
 * GMRES::SubcellReconstruction --                                        *
 *        Performs the subcell reconstruction of solution state           *
 *        within a given cell (i,j) of the computational mesh for         *
 *        the specified quadrilateral solution block.                     *
 **************************************************************************/
inline void GMRES::SubcellReconstruction(const int i, 
					 const int j,
					 const int Limiter) {
    
  int n, n_pts, i_index[8], j_index[8];
  double u0, u0Min, u0Max, uQuad[4], phi_k;
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
             (SolBlk_ptr->Grid.BCtypeW[j] != BC_NONE)) {
    if (j == JCl-Nghost+1 || j == JCu+Nghost-1) {
       n_pts = 0;
    } else if (SolBlk_ptr->Grid.BCtypeW[j] == BC_PERIODIC ||
               SolBlk_ptr->Grid.BCtypeW[j] == BC_NEUMANN ||
               SolBlk_ptr->Grid.BCtypeW[j] == BC_ROBIN) {
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
             (SolBlk_ptr->Grid.BCtypeE[j] != BC_NONE)) {
    if (j == JCl-Nghost+1 || j == JCu+Nghost-1) {
       n_pts = 0;
    } else if (SolBlk_ptr->Grid.BCtypeE[j] == BC_PERIODIC ||
               SolBlk_ptr->Grid.BCtypeE[j] == BC_NEUMANN ||
               SolBlk_ptr->Grid.BCtypeE[j] == BC_ROBIN) {
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
             (SolBlk_ptr->Grid.BCtypeS[i] != BC_NONE)) {
    if (i == ICl-Nghost+1 || i == ICu+Nghost-1) {
       n_pts = 0;
    } else if (SolBlk_ptr->Grid.BCtypeS[i] == BC_PERIODIC ||
               SolBlk_ptr->Grid.BCtypeS[i] == BC_NEUMANN ||
               SolBlk_ptr->Grid.BCtypeS[i] == BC_ROBIN) {
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
             (SolBlk_ptr->Grid.BCtypeN[i] != BC_NONE)) {
    if (i == ICl-Nghost+1 || i == ICu+Nghost-1) {
       n_pts = 0;
    } else if (SolBlk_ptr->Grid.BCtypeN[i] == BC_PERIODIC ||
               SolBlk_ptr->Grid.BCtypeN[i] == BC_NEUMANN ||
               SolBlk_ptr->Grid.BCtypeN[i] == BC_ROBIN) {
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
          dX = SolBlk_ptr->Grid.Cell[ i_index[n] ][ j_index[n] ].Xc - 
               SolBlk_ptr->Grid.Cell[i][j].Xc;
          if (vector_switch) {
             Du = z(index(i_index[n],j_index[n])) - 
                  z(index(i,j));
          } else {
             Du = x(index(i_index[n],j_index[n])) - 
                  x(index(i,j));
          } /* endif */
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
      SolBlk_ptr->dudx[i][j] = (DuDx_ave*DyDy_ave-DuDy_ave*DxDy_ave)/
                               (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
      SolBlk_ptr->dudy[i][j] = (DuDy_ave*DxDx_ave-DuDx_ave*DxDy_ave)/
                               (DxDx_ave*DyDy_ave-DxDy_ave*DxDy_ave);
  
      // Calculate slope limiter.
      if (!SolBlk_ptr->Freeze_Limiter) {
         if (vector_switch) {
            u0 = z(index(i,j));
         } else {
            u0 = x(index(i,j));
         } /* endif */
         u0Min = u0;
         u0Max = u0Min;
         for ( n = 0 ; n <= n_pts-1 ; ++n ) {
            if (vector_switch) {
               u0Min = min(u0Min, z(index(i_index[n],j_index[n])));
               u0Max = max(u0Max, z(index(i_index[n],j_index[n])));
               u0 = z(index(i,j));
            } else {
               u0Min = min(u0Min, x(index(i_index[n],j_index[n])));
               u0Max = max(u0Max, x(index(i_index[n],j_index[n])));
               u0 = x(index(i,j));
            } /* endif */
         } /* endfor */
  
         dX = SolBlk_ptr->Grid.xfaceE(i, j)-SolBlk_ptr->Grid.Cell[i][j].Xc;
         uQuad[0] = u0 + 
                    SolBlk_ptr->dudx[i][j]*dX.x +
                    SolBlk_ptr->dudy[i][j]*dX.y ;
         dX = SolBlk_ptr->Grid.xfaceW(i, j)-SolBlk_ptr->Grid.Cell[i][j].Xc;
         uQuad[1] = u0 + 
                    SolBlk_ptr->dudx[i][j]*dX.x +
                    SolBlk_ptr->dudy[i][j]*dX.y ;
         dX = SolBlk_ptr->Grid.xfaceN(i, j)-SolBlk_ptr->Grid.Cell[i][j].Xc;
         uQuad[2] = u0 + 
                    SolBlk_ptr->dudx[i][j]*dX.x +
                    SolBlk_ptr->dudy[i][j]*dX.y ;
         dX = SolBlk_ptr->Grid.xfaceS(i, j)-SolBlk_ptr->Grid.Cell[i][j].Xc;
         uQuad[3] = u0 + 
                    SolBlk_ptr->dudx[i][j]*dX.x +
                    SolBlk_ptr->dudy[i][j]*dX.y ;
  
         switch(Limiter) {
           case LIMITER_ONE :
             phi_k = ONE;
             break;
           case LIMITER_ZERO :
             phi_k = ZERO;
             break;
           case LIMITER_BARTH_JESPERSEN :
             phi_k = Limiter_BarthJespersen(uQuad, u0, 
                                            u0Min, u0Max, 4);
             break;
           case LIMITER_VENKATAKRISHNAN :
             phi_k = Limiter_Venkatakrishnan(uQuad, u0, 
                                             u0Min, u0Max, 4);
             break;
           case LIMITER_VANLEER :
             phi_k = Limiter_VanLeer(uQuad, u0, 
                                     u0Min, u0Max, 4);
             break;
           case LIMITER_VANALBADA :
             phi_k = Limiter_VanAlbada(uQuad, u0, 
                                       u0Min, u0Max, 4);
             break;
           default:
             phi_k = Limiter_BarthJespersen(uQuad, u0, 
                                            u0Min, u0Max, 4);
             break;
         } /* endswitch */
  
         SolBlk_ptr->phi[i][j] = phi_k;
      } /* endif */
  } else {
      SolBlk_ptr->dudx[i][j] = ZERO;
      SolBlk_ptr->dudy[i][j] = ZERO; 
      SolBlk_ptr->phi[i][j]  = ZERO;
  } /* endif */

}

/*******************************************************************************
 * GMRES::LoadSendBuffer_Flux_F2C -- Loads send message buffer for             *
 *                                   fine to coarse block message              *
 *                                   passing of conservative                   *
 *                                   solution fluxes.                          *
 *******************************************************************************/
inline int GMRES::LoadSendBuffer_Flux_F2C(double *buffer,
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
 * GMRES::UnloadReceiveBuffer_Flux_F2C -- Unloads receive message              *
 *                                        buffer for fine to coarse            *
 *                                        block message passing of             *
 *                                        conservative solution fluxes.        *
 *******************************************************************************/
inline int GMRES::UnloadReceiveBuffer_Flux_F2C(double *buffer,
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

/************************************************************
 * FUNCTIONS REQUIRED BY GMRES ALGORITHM                    *
 ************************************************************/

/********************************************************
 * Routine: GMRES_RightPrecon_MatrixFree                *
 *                                                      *
 * Right-preconditioned Matrix-free GMRES rountine.     *
 *                                                      *
 ********************************************************/
template < class AdvectDiffuse2D_Quad_Block, 
           class AdaptiveBlock2D_List,
           class AdvectDiffuse2D_Input_Parameters,
           class GMRES, 
           class Preconditioner, 
           class Real>
int 
GMRES_RightPrecon_MatrixFree(AdvectDiffuse2D_Quad_Block *Soln_ptr,
		             AdaptiveBlock2D_List &Soln_Block_List,
                             AdvectDiffuse2D_Input_Parameters &Input_Parameters,
      		             GMRES *G, 
		             const Preconditioner *M, 
		             int &m, int &max_iter, Real &tol)

{
  int NBLK = Soln_Block_List.Nblk;
  int i, j = 1, k;
  int xcount, ycount, Bcount;
  int error_flag;

  double total_residual;
  double epsilon;
  double total_norm_uo= 0.0;
  double total_norm_r = 0.0;
  double total_norm_b = 0.0;
  double total_norm_z = 0.0;
  double total_norm_x = 0.0;
  double total_beta   = 0.0;
  double total_H_temp = 0.0;

  for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
    if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {

      /***********************************************************
       * Note: The initial solution block residual has been      *
       * calculated outside the GMRES routine and being stored   *
       * in each G[Bcount].b vector.                             * 
       *                                                         *
       * Also, since u = uo in the beginning the calculation,    *
       * the residual (RHS) can be obtained using uo. However,   *
       * because G.b is -R, therefore, G.Q2 [= R(u)] should be   *
       * negative of G.b.                                        *
       ***********************************************************/
      G[Bcount].Q2 = -1.0 *  G[Bcount].b; 
      
      /*************************/
      /** Matrix-Free Version **/
      /*************************/
      total_norm_b = total_norm_b + sqr(norm(G[Bcount].b));
      
      for (xcount = Soln_ptr[Bcount].ICl-Soln_ptr[Bcount].Nghost;  
	   xcount <= Soln_ptr[Bcount].ICu+Soln_ptr[Bcount].Nghost;  
	   ++xcount ) { 
	for (ycount = Soln_ptr[Bcount].JCu+Soln_ptr[Bcount].Nghost;  
	     ycount > Soln_ptr[Bcount].JCu+G[Bcount].overlap;  
	     --ycount ) { 
	  G[Bcount].r(G[Bcount].index(xcount, ycount)) =  
	    Soln_ptr[Bcount].uo[xcount][ycount]; 
	  
	} /* endfor */
      } /* endfor */
      total_norm_uo = total_norm_uo + sqr(norm(G[Bcount].r));
      
      G[Bcount].r = G[Bcount].b; // instead of calculating 
                                 // G.r = (G.b - G.A * G.x) 
                                 // since the initial guess (G.x) 
                                 // is always zero.
      
    } /* endif */
  } /* endfor */ 
  
  total_norm_uo = sqrt(CFDkit_Summation_MPI(total_norm_uo));
  total_norm_b = sqrt(CFDkit_Summation_MPI(total_norm_b));
  total_norm_r = total_norm_b;
  total_beta   = total_norm_r;

  /* Calculate residual for the entire domain. */ 
  if (total_norm_b == 0.0) {
    total_residual = 0;
  } else {
    total_residual = total_norm_r / total_norm_b;
  }
  
  if (total_residual <= tol) {
    tol = total_residual;
    max_iter = 0;
    if (CFDkit_Primary_MPI_Processor()) { 
      cout << "** norm(G.b)=ZERO **" << endl; cout.flush();
    }
    return 0;
  }
  
  while (j <= max_iter) {
    
    for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
      if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {  
	G[Bcount].v[0] = G[Bcount].r * (1.0 / total_beta);
	G[Bcount].s = 0.0;
	G[Bcount].s(0) =  total_beta;

      } /* endif */
    } /* endfor */ 
    
    for (i = 0; i < m && j <= max_iter; i++, j++) {
      
      for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	  
	  /* Set vector switch. */
	  G[Bcount].vector_switch = 1;
	  
	  G[Bcount].z = M[Bcount].solve(G[Bcount].v[i]);
	  
	} /* endif */
      } /* endfor */       
      
      ////////////////////////////////////////////////////////////////////
      ///////////////////* BEGIN MESSAGE PASSING *////////////////////////
      ////////////////////////////////////////////////////////////////////
      
      /* MPI barrier to ensure processor synchronization. */
      CFDkit_Barrier_MPI();  
      
      /* Send solution information between neighbouring blocks.*/
      error_flag = Send_All_Messages(G, 
                                     Soln_Block_List,
                                     NUM_VAR_ADVECTDIFFUSE2D, 
                                     OFF);
      if (error_flag) {
	cout << "\n AdvectDiffuse2D NKS ERROR: AdvectDiffuse2D message passing error on processor "
	     << Soln_Block_List.ThisCPU
	     << ".\n";
	cout.flush();
      } /* endif */
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
      
      ////////////////////////////////////////////////////////////////////
      ///////////////////*  END OF MESSAGE PASSING  */////////////////////    
      ////////////////////////////////////////////////////////////////////  
      
      for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) { 
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) { 
	  
	  /* Apply boundary condition on z vector. */
	  GMRES_BCs(Soln_ptr[Bcount], G[Bcount].z, G[Bcount]); 
	  
	} /* endif */
      } /* endfor */
      
      /* Calculate global 2-norm of z. */
      total_norm_z = 0.0;
      for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	  G[Bcount].w = 0;
	  total_norm_z = total_norm_z + sqr(norm(G[Bcount].z));
	  
	} /* endif */
      } /* endfor */
      total_norm_z = sqrt(CFDkit_Summation_MPI(total_norm_z));
      
      /* Determine epsilon based on 2-norm of z. */
      //epsilon = 1.e-08/total_norm_z; 
      epsilon = 1000;

      /* BEGIN MATRIX-FREE ************************************/

      // Determine u0+epsilon*z and evaluate residual.      
      for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	  
	  // Calculate U = uo + epsilon * z
	  for (ycount  = Soln_ptr[Bcount].JCl-Soln_ptr[Bcount].Nghost; 
	       ycount <= Soln_ptr[Bcount].JCu+Soln_ptr[Bcount].Nghost; 
	       ++ycount ) {
	    for (xcount = Soln_ptr[Bcount].ICl-Soln_ptr[Bcount].Nghost; 
		 xcount <= Soln_ptr[Bcount].ICu+Soln_ptr[Bcount].Nghost; 
		 ++xcount ) {
	      Soln_ptr[Bcount].U[xcount][ycount].u =  
		Soln_ptr[Bcount].uo[xcount][ycount] +  
		epsilon * G[Bcount].z[G[Bcount].index(xcount, ycount)]; 
	      
	    } /* endfor */ 
	  } /* endfor */ 

	  // Calculate residual R(U = uo + epsilon * z)
	  dUdt_Residual_Evaluation(Soln_ptr[Bcount], 
				   Input_Parameters);
	  
	} /* endif */ 
      } /* endfor */

      // Send boundary flux corrections at block interfaces with resolution changes.
      error_flag = Send_Conservative_Flux_Corrections(Soln_ptr, 
			        		      Soln_Block_List,
						      NUM_VAR_ADVECTDIFFUSE2D);
      if (error_flag) {
	cout << "\n AdvectDiffuse2D GMRES ERROR: AdvectDiffuse2D flux correction message passing error on processor "
	     << Soln_Block_List.ThisCPU
	     << ".\n";
	cout.flush();
      } /* endif */
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
	  
      // Apply boundary flux corrections to residual to ensure that method is conservative.
      Apply_Boundary_Flux_Corrections(Soln_ptr, 
			              Soln_Block_List);

      // Copy residual to Q1 vector and apply finite time step modification. 
      for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	  G[Bcount].copy_residual(G[Bcount].Q1, Soln_ptr[Bcount]);
	  
	  /* Apply finite time step */
	  G[Bcount].w = ZERO;
	  for (ycount  = Soln_ptr[Bcount].JCl; 
	       ycount <= Soln_ptr[Bcount].JCu; ++ycount ) {
	    for (xcount = Soln_ptr[Bcount].ICl; 
		 xcount <= Soln_ptr[Bcount].ICu; ++xcount ) {
	      G[Bcount].w[G[Bcount].index(xcount, ycount)] = 
		G[Bcount].z[G[Bcount].index(xcount, ycount)] / 
		G[Bcount].dTime;
		//Soln_ptr[Bcount].dt[xcount][ycount]; local time step instead of global?
	    } /* endfor */
	  } /* endfor */
	  
	  /* Calculate [R(uo-epsilon * z) - R(uo)]/epsilon - z/dtime. */
	  G[Bcount].w = (G[Bcount].Q1-G[Bcount].Q2)/epsilon - G[Bcount].w;
	  
	} /* endif */
      } /* endfor */
      
      /* END OF MATRIX-FREE ***********************************/
      
      for (k = 0; k <= i; k++) {
	/***************************************/
	/* Calculate global G[Bcount].H(k, i). */
	total_H_temp = 0.0;
	for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	  if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	    total_H_temp = total_H_temp + dot(G[Bcount].w, G[Bcount].v[k]);
	  } /* endif */
	} /* endfor */ 
	total_H_temp = CFDkit_Summation_MPI(total_H_temp);
	/***************************************/
	
	for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	  if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	    G[Bcount].H(k, i) = total_H_temp;
	    
	    /* Update w vector using global G[Bcount].H(k, i). */
	    G[Bcount].w -= G[Bcount].H(k, i) * G[Bcount].v[k];
	    
	  } /* endif */
	} /* endfor */ 
	
      } /* endfor */ 
      
      /*****************************************/
      /* Calculate global G[Bcount].H(i+1, i). */
      total_H_temp = 0.0;
      for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	  total_H_temp = total_H_temp + sqr(norm(G[Bcount].w));
	} /* endif */
      } /* endfor */ 
      total_H_temp = sqrt(CFDkit_Summation_MPI(total_H_temp));
      /*****************************************/
      
      /*******************************************/
      /* Calculate v[i+1] vector for all blocks. */
      for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {  
	  G[Bcount].H(i+1, i) = total_H_temp;   
	  G[Bcount].v[i+1] = G[Bcount].w * (1.0 / G[Bcount].H(i+1, i));
	} /* endif */
      } /* endfor */ 
      /*******************************************/

      for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {  
	    for (k = 0; k < i; k++) { 
	      RPMF_ApplyPlaneRotation(G[Bcount].H(k,i), G[Bcount].H(k+1,i), 
				      G[Bcount].cs(k), G[Bcount].sn(k));
	    }

	    RPMF_GeneratePlaneRotation(G[Bcount].H(i,i), G[Bcount].H(i+1,i), 
				       G[Bcount].cs(i), G[Bcount].sn(i));

	    RPMF_ApplyPlaneRotation(G[Bcount].H(i,i), G[Bcount].H(i+1,i), 
				    G[Bcount].cs(i), G[Bcount].sn(i));

	    RPMF_ApplyPlaneRotation(G[Bcount].s(i), G[Bcount].s(i+1), 
				    G[Bcount].cs(i), G[Bcount].sn(i));   
	} /* endif */
      } /* endfor */ 
      
      if ((total_residual = abs(G[0].s(i+1)) / total_norm_b) < tol) {

	error_flag = RPMF_Update(Soln_ptr, Soln_Block_List, 
				 G, M, i, tol, 
				 total_norm_b);

	tol = total_residual;
	max_iter = j;
	return 0;

      } /* endif */
    } /* endfor */
  
    if (CFDkit_Primary_MPI_Processor()) { 
      cout <<"###### RESTART " << int(j/m) << " #####" 
           << " total_residual = " << total_residual << " tol = " << tol << endl;
    } /* endif */

    error_flag = RPMF_Update(Soln_ptr, Soln_Block_List, 
			     G, M, m-1, tol, 
			     total_norm_b);

    /* Calculate global 2-norm of dU vector (x vector). */
    total_norm_x = 0.0;
    for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
      if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	total_norm_x = total_norm_x + sqr(norm(G[Bcount].x)); 
	
      } /* endif */
    } /* endfor */ 
    total_norm_x = sqrt(CFDkit_Summation_MPI(total_norm_x));

    /*************************/
    /** Matrix-Free Version **/
    /*************************/

    if (total_norm_x == 0.0) {
      for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	  G[Bcount].Ax=0;
	  
	} /* endif */
      } /* endfor */  
    } else {

      /* Global epsilon based on 2-norm of x. */
      //epsilon = 1.e-08*total_norm_uo/total_norm_x;
      epsilon = 1000;

      /* BEGIN MATRIX-FREE ************************************/

     // Determine u0+epsilon*x and evaluate residual. 
      for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	  for (ycount  = Soln_ptr[Bcount].JCl-Soln_ptr[Bcount].Nghost; 
	       ycount <= Soln_ptr[Bcount].JCu+Soln_ptr[Bcount].Nghost;
	       ++ycount ) {
	    for (xcount = Soln_ptr[Bcount].ICl-Soln_ptr[Bcount].Nghost; 
		 xcount <= Soln_ptr[Bcount].ICu+Soln_ptr[Bcount].Nghost;
		 ++xcount ) {
	      Soln_ptr[Bcount].U[xcount][ycount].u = 
		Soln_ptr[Bcount].uo[xcount][ycount] + 
		epsilon * G[Bcount].x[G[Bcount].index(xcount, ycount)];

	    } /* endfor */  
	  } /* endfor */  
	  
	  /* Calculate residual R(U = uo + epsilon * z) */
	  dUdt_Residual_Evaluation(Soln_ptr[Bcount],
				   Input_Parameters);
	  
	} /* endif */
      } /* endfor */  

      // Send boundary flux corrections at block interfaces with resolution changes.
      error_flag = Send_Conservative_Flux_Corrections(Soln_ptr, 
			        		      Soln_Block_List,
						      NUM_VAR_ADVECTDIFFUSE2D);
      if (error_flag) {
	 cout << "\n AdvectDiffuse2D GMRES ERROR: AdvectDiffuse2D flux correction message passing error on processor "
	      << Soln_Block_List.ThisCPU
	      << ".\n";
	 cout.flush();
      } /* endif */
      error_flag = CFDkit_OR_MPI(error_flag);
      if (error_flag) return (error_flag);
	  
      // Apply boundary flux corrections to residual to ensure that method is conservative.
      Apply_Boundary_Flux_Corrections(Soln_ptr, 
			              Soln_Block_List);

      // Copy residual to Q1 vector and apply finite time step modification.
      for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
	if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	  G[Bcount].copy_residual(G[Bcount].Q1, Soln_ptr[Bcount]);
	  
	  /* Apply finite time step */
	  G[Bcount].Ax = ZERO;
	  for (ycount = Soln_ptr[Bcount].JCl; 
	       ycount <= Soln_ptr[Bcount].JCu; ++ycount ) {
	    for (xcount = Soln_ptr[Bcount].ICl; 
		 xcount <= Soln_ptr[Bcount].ICu; ++xcount ) {
	      G[Bcount].Ax[G[Bcount].index(xcount, ycount)] = 
		G[Bcount].x[G[Bcount].index(xcount, ycount)] / 
		G[Bcount].dTime;
		//Soln_ptr[Bcount].dt[xcount][ycount];

	    } /* endfor */  
	  } /* endfor */  

	  G[Bcount].Ax = (G[Bcount].Q1-G[Bcount].Q2)/epsilon - G[Bcount].Ax;

	} /* endif */
      } /* endfor */   
      
      /* END OF MATRIX-FREE ***********************************/

    } /* endif */

    total_norm_r = 0.0;
    for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
      if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	G[Bcount].r = G[Bcount].b - G[Bcount].Ax;
	total_norm_r = total_norm_r + sqr(norm(G[Bcount].r));
	
      } /* endif */
    } /* endfor */ 
    total_norm_r = sqrt(CFDkit_Summation_MPI(total_norm_r)); 
    total_beta = total_norm_r;
    
    if ((total_residual = total_beta / total_norm_b) < tol) {
      tol = total_residual;
      max_iter = j;
      return 0;
    } /* endif */
  } /* endwhile */
  
  tol = total_residual;
  return 0;
  
}/* End of GMRES_RightPrecon_MatrixFree algorithm */

/********************************************************
 * Routine: RPMF_Update                                 *
 *                                                      *
 * Updates and message passes all solutions.            *
 *                                                      *
 ********************************************************/
template < class AdvectDiffuse2D_Quad_Block, 
           class AdaptiveBlock2D_List,
           class GMRES,
           class Preconditioner, 
           class Real>
int  
RPMF_Update(AdvectDiffuse2D_Quad_Block *Soln_ptr, 
	    AdaptiveBlock2D_List &Soln_Block_List, 
	    GMRES *G,
            const Preconditioner *M,
	    Real k, 
	    double tol,
	    double total_norm_b)

{
  int i,j=0;
  int NBLK = Soln_Block_List.Nblk;
  int Bcount, xcount, ycount, error_flag;

  VECTOR_double *y = new VECTOR_double[NBLK];
  for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
    if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      y[Bcount] = G[Bcount].s;
    }
  }

  for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
    if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      // Backsolve:  
      for (i = k; i >= 0; i--) {
	y[Bcount](i) /= G[Bcount].H(i,i);
	for (j = i - 1; j >= 0; j--) {
	  y[Bcount](j) -= G[Bcount].H(j,i) * y[Bcount](i);
	}
      }
      for (j = 0; j <= k; j++){
	G[Bcount].x += M[Bcount].solve(G[Bcount].v[j] * y[Bcount](j));
      }

      /* Set vector switch. */
      G[Bcount].vector_switch = 0;
      
    } /* endif */ 
  } /* endfor */

  ////////////////////////////////////////////////////////////////////
  ///////////////////* BEGIN MESSAGE PASSING *////////////////////////
  ////////////////////////////////////////////////////////////////////
  
  /* MPI barrier to ensure processor synchronization. */
  CFDkit_Barrier_MPI();
  
  /* Send solution information between neighbouring blocks. */
  error_flag = Send_All_Messages(G, 
                                 Soln_Block_List,
                                 NUM_VAR_ADVECTDIFFUSE2D, 
                                 OFF);
  if (error_flag) {
     cout << "\n AdvectDiffuse2D GMRES ERROR: AdvectDiffuse2D message passing error on processor "
          << Soln_Block_List.ThisCPU
          << ".\n";
     cout.flush();
  } /* endif */
  error_flag = CFDkit_OR_MPI(error_flag);
  if (error_flag) cout << " CFDkit_OR_MPI(error_flag) " << endl;

  ////////////////////////////////////////////////////////////////////
  ///////////////////*  END MESSAGE PASSING  *////////////////////////    
  //////////////////////////////////////////////////////////////////// 

  for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
    if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
      /* Apply boundary condition on dU. */ 
      GMRES_BCs(Soln_ptr[Bcount], G[Bcount].x, G[Bcount]);
    } /* endif */
  } /* endfor */

  return 0;

} /* End of RPMF_Update Function. */

template < class Real >
Real 
abs(Real x)
{
  return (x > 0 ? x : -x);
}

template<class Real> 
void RPMF_GeneratePlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
{
  if (dy == 0.0) {
    cs = 1.0;
    sn = 0.0;
  } else if (abs(dy) > abs(dx)) {
    Real temp = dx / dy;
    sn = 1.0 / sqrt( 1.0 + temp*temp );
    cs = temp * sn;
  } else {
    Real temp = dy / dx;
    cs = 1.0 / sqrt( 1.0 + temp*temp );
    sn = temp * cs;
  }
} /* End of RPMF_GeneratePlaneRotation Function. */

template<class Real> 
void RPMF_ApplyPlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
{
  Real temp  =  cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = temp;
} /* End of RPMF_ApplyPlaneRotation Function. */

/********************************************************
 * Routine: GMRES_BCs                                   *
 *                                                      *
 * Apply boundary conditions at boundaries of the       *
 * specified GMRES vector.                              *
 *                                                      *
 ********************************************************/
template < class AdvectDiffuse2D_Quad_Block, 
           class Vector, 
           class GMRES>
void
GMRES_BCs(AdvectDiffuse2D_Quad_Block &SolnBlk, 
	  Vector &v, 
	  GMRES &G) 
{
    int i, j;
    double dx_norm, du, dudx;
    Vector2D dX;

    for ( j = SolnBlk.JCl-SolnBlk.Nghost ; 
	  j <= SolnBlk.JCu+SolnBlk.Nghost ; ++j ) {
      if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
           (j < SolnBlk.JCl && 
            (SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_NONE ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_NEUMANN ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_DIRICHLET||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICl-1] == BC_ROBIN) ) ||
           (j > SolnBlk.JCu && 
            (SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_NONE ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_DIRICHLET||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_NEUMANN ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICl-1] == BC_ROBIN) ) ) {
        switch(SolnBlk.Grid.BCtypeW[j]) {
	case BC_NONE :
            break;
	case BC_DIRICHLET :
	  v[G.index(SolnBlk.ICl-1,j)] = 0.0;
	  v[G.index(SolnBlk.ICl-2,j)] = 0.0;
	  break;
	case BC_NEUMANN :
	  v[G.index(SolnBlk.ICl-1,j)] = v[G.index(SolnBlk.ICl,j)];
	  v[G.index(SolnBlk.ICl-2,j)] = v[G.index(SolnBlk.ICl,j)];
	  break;
	case BC_ROBIN :
	  break;
	default:
	  break;
        } /* endswitch */ 
      } /* endif */
	
      if ( (j >= SolnBlk.JCl && j <= SolnBlk.JCu) ||
           (j < SolnBlk.JCl && 
            (SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_NONE ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_NEUMANN ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_DIRICHLET ||
             SolnBlk.Grid.BCtypeS[SolnBlk.ICu+1] == BC_ROBIN) ) ||
           (j > SolnBlk.JCu && 
            (SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_NONE ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_PERIODIC ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_DIRICHLET ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_NEUMANN ||
             SolnBlk.Grid.BCtypeN[SolnBlk.ICu+1] == BC_ROBIN) ) ) {
        switch(SolnBlk.Grid.BCtypeE[j]) {
          case BC_NONE :
            break;
          case BC_DIRICHLET :
	    v[G.index(SolnBlk.ICu+1,j)] = 0.0;
	    v[G.index(SolnBlk.ICu+2,j)] = 0.0;
	    break;
          case BC_NEUMANN :
	    v[G.index(SolnBlk.ICu+1,j)] = v[G.index(SolnBlk.ICu,j)];
	    v[G.index(SolnBlk.ICu+2,j)] = v[G.index(SolnBlk.ICu,j)];
            break;
          case BC_ROBIN :
            break;
          default:
            break;
        } /* endswitch */
      } /* endif */
    } /* endfor */

    for ( i = SolnBlk.ICl-SolnBlk.Nghost ; 
	  i <= SolnBlk.ICu+SolnBlk.Nghost ; ++i ) {
      if ( (i >= SolnBlk.ICl && i <= SolnBlk.ICu) ||
	   (i < SolnBlk.ICl && 
	    (SolnBlk.Grid.BCtypeW[SolnBlk.JCl-1] == BC_NONE ||
	     SolnBlk.Grid.BCtypeW[SolnBlk.JCl-1] == BC_PERIODIC ||
	     SolnBlk.Grid.BCtypeW[SolnBlk.JCl-1] == BC_NEUMANN ||
	     SolnBlk.Grid.BCtypeW[SolnBlk.JCl-1] == BC_DIRICHLET||
	     SolnBlk.Grid.BCtypeW[SolnBlk.JCl-1] == BC_ROBIN) ) ||
	   (i > SolnBlk.ICu && 
	    (SolnBlk.Grid.BCtypeE[SolnBlk.JCl-1] == BC_NONE ||
	     SolnBlk.Grid.BCtypeE[SolnBlk.JCl-1] == BC_PERIODIC ||
	     SolnBlk.Grid.BCtypeE[SolnBlk.JCl-1] == BC_DIRICHLET ||
	     SolnBlk.Grid.BCtypeE[SolnBlk.JCl-1] == BC_NEUMANN ||
	     SolnBlk.Grid.BCtypeE[SolnBlk.JCl-1] == BC_ROBIN) ) ) {
	switch(SolnBlk.Grid.BCtypeS[i]) {
	case BC_NONE :
	  break;
	case BC_DIRICHLET :
	  v[G.index(i, SolnBlk.JCl-1)] = 0.0;
	  v[G.index(i, SolnBlk.JCl-2)] = 0.0;
	  break;
	case BC_NEUMANN :
	  v[G.index(i, SolnBlk.JCl-1)] = v[G.index(i, SolnBlk.JCl)];
	  v[G.index(i, SolnBlk.JCl-2)] = v[G.index(i, SolnBlk.JCl)];
	  break;
	case BC_ROBIN :	
	  break;
	default:

	  break;
        } /* endswitch */
      } /* endif */

      if ( (i >= SolnBlk.ICl && i <= SolnBlk.ICu) ||
	   (i < SolnBlk.ICl && 
	    (SolnBlk.Grid.BCtypeW[SolnBlk.JCu+1] == BC_NONE ||
	     SolnBlk.Grid.BCtypeW[SolnBlk.JCu+1] == BC_PERIODIC ||
	     SolnBlk.Grid.BCtypeW[SolnBlk.JCu+1] == BC_NEUMANN ||
	     SolnBlk.Grid.BCtypeW[SolnBlk.JCu+1] == BC_DIRICHLET||
	     SolnBlk.Grid.BCtypeW[SolnBlk.JCu+1] == BC_ROBIN) ) ||
	   (i > SolnBlk.ICu && 
           (SolnBlk.Grid.BCtypeE[SolnBlk.JCu+1] == BC_NONE ||
            SolnBlk.Grid.BCtypeE[SolnBlk.JCu+1] == BC_PERIODIC ||
	    SolnBlk.Grid.BCtypeE[SolnBlk.JCu+1] == BC_DIRICHLET ||
	    SolnBlk.Grid.BCtypeE[SolnBlk.JCu+1] == BC_NEUMANN ||
	    SolnBlk.Grid.BCtypeE[SolnBlk.JCu+1] == BC_ROBIN) ) ) {
	switch(SolnBlk.Grid.BCtypeN[i]) {
	case BC_NONE :
	  break;
	case BC_DIRICHLET :
	  v[G.index(i, SolnBlk.JCu+1)] = 0.0;
	  v[G.index(i, SolnBlk.JCu+2)] = 0.0;
	  break;
	case BC_NEUMANN :
	  v[G.index(i, SolnBlk.JCu+1)] = v[G.index(i, SolnBlk.JCu)];
	  v[G.index(i, SolnBlk.JCu+2)] = v[G.index(i, SolnBlk.JCu)];
	  break;
	case BC_ROBIN :
	  break;
	default:

	  break;
        } /* endswitch */ 
      } /* endif */
    } /* endfor */
    
} /* End of GMRES_BCs Function */ 

/***********************************************************************
 * Iterative template routine -- GMRES_Algorithm                       *
 *                               (Right-preconditioned matrix-free)    *
 *                                                                     * 
 *  GMRES solves the unsymmetric linear system Ax = b using the        *
 *  Generalized Minimum Residual method                                * 
 *                                                                     *
 * The return value indicates convergence within max_iter (input)      *
 * iterations (0), or no convergence within max_iter iterations (1).   *
 *                                                                     * 
 *  error_flag = GMRES_Algorithm                                       *
 *                    (Soln_ptr, Soln_Block_List, Input_Parameters,    *
 *		       G, gmrestol);                                   *
 *                                                                     *
 * List of neccessary inputs:                                          *
 *                                                                     *
 *         Soln_ptr --  solution block pointer                         *
 *  Soln_Block_List -- local solution block list                       *
 * Input_Parameters --  parameters describing solution                 *
 *                G --  GMRES object containing all function variables *
 *         P_Switch --  preconditioner switch                          *
 *          restart --  number of iterations for each restart          *    
 *         gmrestol --  the residual after the final iteration         *
 *                                                                     *
 * Note:  GMRES follows the algorithm described on p. 20 of the        *
 *        SIAM Templates book.                                         *
 *                                                                     *
 ***********************************************************************/
template < class AdvectDiffuse2D_Quad_Block,
           class AdaptiveBlock2D_List, 
           class GMRES, 
           class Real >
int 
GMRES_Algorithm(AdvectDiffuse2D_Quad_Block *Soln_ptr, 
		AdaptiveBlock2D_List &Soln_Block_List, 
                AdvectDiffuse2D_Input_Parameters &Input_Parameters,
		GMRES *G, 
	        Real &tol)
{
  double gmres_tol;

  int i, j;
  int Bcount;
  int xpts, ypts;
  int NBLK           = Soln_Block_List.Nblk;
  int error_flag     = 0;
  int max_gmres_iter = G[0].max_gmres_iter;

  CompRow_ILUPreconditioner_double *BM = 
    new CompRow_ILUPreconditioner_double[NBLK];
  DiagPreconditioner_double *BD = 
    new DiagPreconditioner_double[NBLK];

  for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
    if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {

      // Initialize vector x and matrix H.
      G[Bcount].x = 0.0;
      G[Bcount].H = 0.0;

      //RESET Maximum, restart iterations and tolerance
      gmres_tol  = tol;

      if (G[0].P_Switch == 1) {
	BM[Bcount] = G[Bcount].BM;
      } else {
	BD[Bcount] = G[Bcount].BD;
      }

    } /* endif */
  } /* endfor */

  if (G[0].P_Switch == 1) {
    error_flag = GMRES_RightPrecon_MatrixFree(Soln_ptr, 
					      Soln_Block_List, 
                                              Input_Parameters,
                                              G, BM,
				              G[0].m, 
					      max_gmres_iter, 
					      gmres_tol);
  } else {
    error_flag = GMRES_RightPrecon_MatrixFree(Soln_ptr, 
					      Soln_Block_List, 
                                              Input_Parameters,
                                              G, BD,
				              G[0].m, 
					      max_gmres_iter, 
					      gmres_tol); 
  } 

  if (CFDkit_Primary_MPI_Processor()) {  
    cout << "        GMRES iter = " << setw(5) << max_gmres_iter 
	 << "  gmrestol = " << setw(10) << gmres_tol << endl;
  } /* endif */

  if (error_flag) {
    if (CFDkit_Primary_MPI_Processor()) {  
      cout << "\nAdvectDiffuse2D GMRES ERROR: Unable to reach the specified convergence tolerance in the maximum number of iterations in block number " 
           << Bcount << ".\n";
      cout.flush();
    }
    return 1;
  } else {
    
    for ( Bcount = 0 ; Bcount < NBLK ; ++Bcount ) {
      if (Soln_Block_List.Block[Bcount].used == ADAPTIVEBLOCK2D_USED) {
	
	// Copy block solution (du) to dudt vector 
	for (j  = Soln_ptr[Bcount].JCl-Soln_ptr[Bcount].Nghost ; 
	     j <= Soln_ptr[Bcount].JCu+Soln_ptr[Bcount].Nghost ; ++j ) {
	  for ( i = Soln_ptr[Bcount].ICl-Soln_ptr[Bcount].Nghost ; 
		i <= Soln_ptr[Bcount].ICu+Soln_ptr[Bcount].Nghost ; ++i ) {
	    Soln_ptr[Bcount].dudt[i][j][1] = 
	      G[Bcount].x[G[Bcount].index(i, j)];

	  } /* endfor */
	} /* endfor */
      } /* endif */
    } /* endfor */
  } /* endif */
  
  delete []BM;
  delete []BD;
  
  return 0;
  
} /* End of GMRES_Algorithm. */

#endif // _ADVECTDIFFUSE2D_QUAD_GMRES_INCLUDED
