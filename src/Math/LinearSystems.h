/* LinearSystems.h:  Header file defining a variety of linear system classes. */

#ifndef _LINEARSYSTEMS_INCLUDED
#define _LINEARSYSTEMS_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include "include/require.h"

/* Include the matrix header file. */

#ifndef _MATRIX_INCLUDED
#include "Matrix.h"
#endif // _MATRIX_INCLUDED

using namespace std;

/* Define some useful constants. */

#define THOMAS_ALGORITHM            0
#define	LU_DECOMPOSITION            1
#define POINT_JACOBI                2
#define GAUSS_SEIDEL                3
#define LS_Householder              4

/* External subroutines for dense systems. */

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

void HouseholderOrthogTransf_col(bool KnownTransformation,int & PivotElement,int StartZeroElem,int &MaxColumnIndex,
				 DenseMatrix &V, int PivotVector, double &UP, DenseMatrix &C, int CLow,
				 int &CUp);

void HouseholderOrthogTransf_row(bool KnownTransformation,int & PivotElement,int StartZeroElem,int &MaxRowIndex,
				 DenseMatrix &V, int PivotVector, double &UP, DenseMatrix &C, int RLow,
				 int RUp);

void ComputeSquaredColumnLength(int j,const DenseMatrix &A,int &lmax,
				double &Hmax, double *H);



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
   assert( N >= 1 ); A = DenseMatrix(N, N); 
   b = ColumnVector(N); x = ColumnVector(N);
}

/**************************************************************
 * DenseSystemLinEqs::deallocate -- Deallocate memory.        *
 **************************************************************/
inline void DenseSystemLinEqs::deallocate(void) {
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

// ********************************** Fast Least Squares Solver ************************************

/********************************************************************
* Routine: Fast_SquaredColumnLength                                 *
*                                                                   *
* Determines the column of matrix A that has the sum of squares of  *
* components in rows j through Rows the greatest.                   *
********************************************************************/
template<int Rows, int Cols>
void ComputeSquaredColumnLength(int j, const double A[][Cols], int &lmax,
				double &Hmax, double H[]){


  int M,N;
  double factor = 0.001;
  int row,col;

  M = Rows - 1;
  N = Cols - 1;

  lmax = j;
  if (j==0){
    for (col=0; col<=N; ++col){	// for each column
      H[col] = 0.0;
      for (row=0; row<=M; ++row) // sum the elements of the rows
	H[col] += A(row,col)*A(row,col);
      if (H[col] > H[lmax])
	lmax = col;
    }
    Hmax = H[lmax];
  }
  else{
    for (col=j; col<=N; ++col){
      H[col] -= A(j-1,col)*A(j-1,col);
      if (H[col] > H[lmax])
	lmax = col;
    }
    if ((Hmax + factor*H[lmax])-Hmax<=0){
      lmax = j;
      for (col=j; col<=N; ++col){	// for each column
	H[col] = 0.0;
	for (row=j; row<=M; ++row){ // sum the elements of the rows
	  H[col] += A(row,col)*A(row,col);
	}
	if (H[col] > H[lmax])
	  lmax = col;
      }
    }
    Hmax = H[lmax];
  }
}



#endif /* _LINEARSYSTEMS_INCLUDED  */
