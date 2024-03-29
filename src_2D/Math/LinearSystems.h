/*!\file LinearSystems.h
  \brief Header file defining a variety of linear system classes. */

#ifndef _LINEARSYSTEMS_INCLUDED
#define _LINEARSYSTEMS_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include <cmath>

#include "lapackd.h"   /* Include lapack double precision library. */

using namespace std;

/* Include CFFC header files */
#include "Matrix.h"             /* Include the matrix header file. */
#include "../Utilities/EpsilonTol.h" /* Include numerical tolerances header file. */
#include "../Utilities/Utilities.h" /* Include utilities header file. */

/* Define some useful constants. */

#define THOMAS_ALGORITHM            0
#define	LU_DECOMPOSITION            1
#define POINT_JACOBI                2
#define GAUSS_SEIDEL                3
#define LS_Householder              4

/* External subroutines for dense systems. */

extern "C" 
{
  void dgells  (int*, void *, int *, void *, int*, void *, int*,
		double *, double*, int*, int*, int*, int*, double *, int*);
}

extern void Solve_LU_Decomposition(DenseMatrix &A,
                                   ColumnVector &b,
                                   ColumnVector &x);

extern void Solve_Point_Jacobi(DenseMatrix &A,
                               ColumnVector &b,
                               ColumnVector &x);

extern void Solve_Gauss_Seidel(DenseMatrix &A,
                               ColumnVector &b,
                               ColumnVector &x);

void Solve_LS_Householder(DenseMatrix &A,
			  ColumnVector &b,
			  ColumnVector &x,
			  int &krank, double &RNorm);

void Solve_LS_Householder(DenseMatrix &A,
			  DenseMatrix &B,
			  DenseMatrix &X,
			  int &krank, ColumnVector &RNorm);

void Solve_Constrained_LS_Householder(DenseMatrix &A,
				      DenseMatrix &B,
				      DenseMatrix &X,
				      const int NumberOfConstraints);

void HouseholderOrthogTransf_col(bool KnownTransformation, bool MaxAbsValueDetermined,
				 int PivotElement,int StartZeroElem, int MaxColumnIndex,
				 DenseMatrix &V, int PivotVector, double &UP, DenseMatrix &C,
				 int CLow, int CUp);

void HouseholderOrthogTransf_row(bool KnownTransformation, bool MaxAbsValueDetermined,
				 int PivotElement,int StartZeroElem,int MaxRowIndex,
				 DenseMatrix &V, int PivotVector, double &UP, DenseMatrix &C,
				 int RLow, int RUp);

void ComputeSquaredColumnLength(const int &j,const DenseMatrix &A,int &lmax,
				double &Hmax, double *H, const int &M, const int &N);

/* Define the n x n dense system of linear equations class. */

/********************************************************
 * Class: DenseSystemLinEqs                             *
 *                                                      *
 * Member functions                                     *
 *      A          -- Return LHS dense (full) matrix.   *
 *      b          -- Return RHS column vector.         *
 *      x          -- Return solution vector.           *
 *      allocate   -- Allocate memory for system.       *
 *      deallocate -- Deallocate memory for system.     *
 *      solve      -- Solve dense system of equations.  *
 *                                                      *
 * Member operators                                     *
 *      Ax_equal_b -- Dense system of equations.        *
 *                                                      *
 * Ax_equal_b = Ax_equal_b;                             *
 *                                                      *
 ********************************************************/
class DenseSystemLinEqs{
  private:
  public:
    DenseMatrix     A;  // LHS dense (full) matrix.
    ColumnVector  b,x;  // RHS and solution column vectors.
                        // Made public so can access them.
			
    /* Creation, copy, and assignment constructors. */
    DenseSystemLinEqs(void) {
    }

    DenseSystemLinEqs(const DenseSystemLinEqs &Ax_equal_b) {
       A = Ax_equal_b.A; b = Ax_equal_b.b; x = Ax_equal_b.x;
    }

    DenseSystemLinEqs(const DenseMatrix &AA,
   	              const ColumnVector &bb,
		      const ColumnVector &xx) {
       A = AA; b = bb; x = xx;
    }

    /* Destructor. */
    // ~DenseSystemLinEqs(void);
    // Use automatically generated destructor.

    /* Allocate memory for dense system. */
    void allocate(const int N);

    /* Deallocate memory for dense system. */
    void deallocate(void);

    /* Solve dense system of linear equations. */
    void solve(void);
    void solve(const int solver_type);

    /* Assignment operator. */
    // DenseSystemLinEqs operator = (const DenseSystemLinEqs &Ax_equal_b);
    // Use automatically generated assignment operator.

};

/**************************************************************
 * DenseSystemLinEqs::allocate -- Allocate memory.            *
 **************************************************************/
inline void DenseSystemLinEqs::allocate(const int N) {
   double temp = 0.0; // quick cludge to avoid copies of uninitialized b 
                      // in MV_Vector_double& MV_Vector_double::operator=(const MV_Vector_double & m)
   assert( N >= 1 ); A = DenseMatrix(N, N, temp); 
   b = ColumnVector(N,temp); x = ColumnVector(N, temp);
}

/**************************************************************
 * DenseSystemLinEqs::deallocate -- Deallocate memory.        *
 **************************************************************/
inline void DenseSystemLinEqs::deallocate(void) {
  //deallocation really done by MV_Vector_double, and MV_ColMat_double classes destructors
  //this just sets some NULL pointers.
  A = DenseMatrix(); b = ColumnVector(); x = ColumnVector();
}

/**************************************************************
 * DenseSystemLinEqs::solve -- Solve system of equations.     *
 **************************************************************/
inline void DenseSystemLinEqs::solve(void) {
   Solve_LU_Decomposition(A, b, x);
}

inline void DenseSystemLinEqs::solve(const int solver_type) {
   switch(solver_type) {
     case LU_DECOMPOSITION :
        Solve_LU_Decomposition(A, b, x); break;
     case POINT_JACOBI :
        Solve_Point_Jacobi(A, b, x); break;
     case GAUSS_SEIDEL :
        Solve_Gauss_Seidel(A, b, x); break;
     default:
        Solve_LU_Decomposition(A, b, x); break;
   };
}

/* External subroutines for tridiagonal systems. */

extern void Solve_Thomas_Algorithm(TriDiagonalMatrix &A,
                                   ColumnVector &b,
                                   ColumnVector &x);

extern void Solve_Point_Jacobi(TriDiagonalMatrix &A,
                               ColumnVector &b,
                               ColumnVector &x);

extern void Solve_Gauss_Seidel(TriDiagonalMatrix &A,
                               ColumnVector &b,
                               ColumnVector &x);

/* Define the n x n tridiagonal system of linear equations class. */

/********************************************************
 * Class: TriDiagonalSystemLinEqs                       *
 *                                                      *
 * Member functions                                     *
 *      A          -- Return LHS tridiagonal matrix.    *
 *      b          -- Return RHS column vector.         *
 *      x          -- Return solution vector.           *
 *      allocate   -- Allocate memory for system.       *
 *      deallocate -- Deallocate memory for system.     *
 *      solve      -- Solve tridiagonal system of       *
 *                    equations.                        *
 *                                                      *
 * Member operators                                     *
 *      Ax_equal_b -- Tridiagonal system of equations.  *
 *                                                      *
 * Ax_equal_b = Ax_equal_b;                             *
 *                                                      *
 ********************************************************/
class TriDiagonalSystemLinEqs{
  private:
  public:
    TriDiagonalMatrix   A;  // LHS tridiagonal matrix.
    ColumnVector      b,x;  // RHS and solution column vectors.
                            // Made public so can access them.
			
    /* Creation, copy, and assignment constructors. */
    TriDiagonalSystemLinEqs(void) {
       A = TriDiagonalMatrix();
    }

    TriDiagonalSystemLinEqs(const TriDiagonalSystemLinEqs &Ax_equal_b) {
       A = Ax_equal_b.A; b = Ax_equal_b.b; x = Ax_equal_b.x;
    }

    TriDiagonalSystemLinEqs(const TriDiagonalMatrix &AA,
   	                    const ColumnVector &bb,
		            const ColumnVector &xx) {
       A = AA; b = bb; x = xx;
    }

    /* Destructor. */
    // ~TriDiagonalSystemLinEqs(void);
    // Use automatically generated destructor.

    /* Allocate memory for tridiagonal system. */
    void allocate(const int N);

    /* Deallocate memory for tridiagonal system. */
    void deallocate(void);

    /* Solve tridiagonal system of equations. */
    void solve(void);
    void solve(const int solver_type);

    /* Assignment operator. */
    // TriDiagonalSystemLinEqs operator = (const TriDiagonalSystemLinEqs &Ax_equal_b);
    // Use automatically generated assignment operator.

};

/**************************************************************
 * TriDiagonalSystemLinEqs::allocate -- Allocate memory.      *
 **************************************************************/
inline void TriDiagonalSystemLinEqs::allocate(const int N) {
   assert( N >= 1 ); A.allocate(N); b = ColumnVector(N); x = ColumnVector(N);
}

/**************************************************************
 * TriDiagonalSystemLinEqs::deallocate -- Deallocate memory.  *
 **************************************************************/
inline void TriDiagonalSystemLinEqs::deallocate(void) {
   A.deallocate(); b = ColumnVector(); x = ColumnVector();
}

/**************************************************************
 * TriDiagonalSystemLinEqs::solve -- Solve system of eqs.     *
 **************************************************************/
inline void TriDiagonalSystemLinEqs::solve(void) {
   Solve_Thomas_Algorithm(A, b, x);
}

inline void TriDiagonalSystemLinEqs::solve(const int solver_type) {
   switch(solver_type) {
     case THOMAS_ALGORITHM :
     case LU_DECOMPOSITION :
        Solve_Thomas_Algorithm(A, b, x); break;
     case POINT_JACOBI :
        Solve_Point_Jacobi(A, b, x); break;
     case GAUSS_SEIDEL :
        Solve_Gauss_Seidel(A, b, x); break;
     default:
        Solve_Thomas_Algorithm(A, b, x); break;
   };
}

inline void Solve_LAPACK_dgesv(DenseMatrix &A,
			       ColumnVector &B) {

  static integer INFO;
  static integer N(A.size(0));
  static integer NRHS(1);
  int *IPIV = new int[N];
  
  F77NAME(dgesv)(&N, &NRHS, &A(0,0), &N, IPIV, &B(0), &N, &INFO);

  delete []IPIV; IPIV = NULL;
}

inline void Solve_LS_Householder_F77(DenseMatrix &A,
				     DenseMatrix &B,
				     int &krank, const int & NumberOfParameters,
				     int _NROW_, int _NCOL_){
  
#ifdef _USE_ESSL_LIB

  /* Initialize variables */
  static double RCOND(EpsilonTol::epsilon);
  static integer NRHS, NROW, NCOL, LWORK;
  static double *WORK;
  static double rn;			// it doesn't get used in the current setup
  static int k;			// number of columns of matrix A used in the solution.
  static int i,j;


  NRHS = NumberOfParameters;
  NROW = _NROW_;
  NCOL = _NCOL_;
  LWORK = 4*NCOL + max(NCOL, NRHS) + 3;
  WORK = new double[LWORK];
  DenseMatrix X(NCOL,NRHS);

  /* Call Fortran subroutine from IBM ESSL library */  
  dgells (0, &A(0,0), &NROW, &B(0,0), &NROW, &X(0,0), &NCOL,
	  &rn, &RCOND, &NROW, &NCOL, &NRHS, &k, WORK, &LWORK);

  // Copy solution X into the elements of B
  for (j = 0; j < NRHS; ++j){
    for (i = 0; i < NCOL; ++ i){
      B(i,j) = X(i,j);
    }
  }

  /* Free memory */
  delete [] WORK; WORK = NULL;

#else

  //  static char TRANS('N');
  static integer INFO;
  static double RCOND(EpsilonTol::epsilon);
  static integer NRHS, NROW, NCOL, LWORK;
  static double *WORK;
  static int *JPVT;

  /* Initialize variables */
  NRHS = NumberOfParameters;
  NROW = _NROW_;
  NCOL = _NCOL_;
  LWORK = max( min(NROW, NCOL) + 3*NCOL+1, 2*min(NROW,NCOL) +  NumberOfParameters);
  WORK = new double[LWORK];
  JPVT = new int[NCOL];
  for (int i=0; i<NCOL; ++i){
    JPVT[i] = 0;
  }

  /* Call Fortran subroutine from lapack */
  F77NAME(dgelsy)(&NROW, &NCOL, &NRHS, &A(0,0), &NROW,
		  &B(0,0), &NROW, JPVT, &RCOND, &krank, WORK, &LWORK, &INFO);

  /* Free memory */
  delete [] WORK; WORK = NULL;
  delete [] JPVT; JPVT = NULL;

#endif

}

inline void Solve_LS_Householder_F77(DenseMatrix &A,
				     ColumnVector &b,
				     int &krank,
				     int _NROW_, int _NCOL_){

  DenseMatrix B(&b(0),b.size(),1,MV_Matrix_::ref);
  Solve_LS_Householder_F77(A,B,krank,1,_NROW_,_NCOL_);
}

#endif /* _LINEARSYSTEMS_INCLUDED  */
