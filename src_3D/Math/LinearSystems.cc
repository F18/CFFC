/* LinearSystems.cc:  Subroutines for various linear system classes. */

/* Include linear systems header file. */

#ifndef _LINEARSYSTEMS_INCLUDED
#include "LinearSystems.h"
#endif // _LINEARSYSTEMS_INCLUDED

/********************************************************
 * Routine: Solve_LU_Decomposition                      *
 *                                                      *
 * Solves a dense system of equations of the form       *
 *                                                      *
 *               A x = L U x = b,                       *
 *                                                      *
 * where A is a dense (full) matrix, using Gaussian     *
 * elimination carried out via LU factorization.  The   *
 * matrix A is decomposed into a product of lower and   *
 * upper triangular matrices and then the solution x is *
 * found by forward and back substitutions.  The matrix *
 * A is overwritten and contains entries associated     *
 * with the lower and upper matrices, L and U, and the  *
 * soluton vector is returned in the column vector x.   *
 * Column vector b is also overwritten.                 *
 *                                                      *
 ********************************************************/
void Solve_LU_Decomposition(DenseMatrix &A,
                            ColumnVector &b,
                            ColumnVector &x) {

    int i, j, k, N;

    /* Perform the LU decomposition (factorization) of the 
       matrix A without pivoting such that A = LU. The entries
       of L and U are stored in place in the matrix A and L
       is assumed to have values of unity on the main diagonal. */

    N = A.size(0)-1;

    for (i = 0; i <= N-1; ++i) {
       for (j = i + 1; j <= N; ++j) { A(j,i) = A(j,i)/A(i,i); }
       for (j = i + 1; j <= N; ++j) { 
	  for (k = i + 1; k <= N; ++k) {
             A(j,k) -= A(j,i)*A(i,k);
          } /* endfor */
       } /* endfor */
    } /* endfor */

   /* Solve the lower triangular system Lz = b via forward
       substitution. */ 

    for (i = 0; i <= N-1; ++i) {
       x(i) = b(i);
       for (j = i+1; j <= N; ++j) { b(j) -= A(j,i)*x(i); }
    } /* endfor */
    x(N) = b(N);

    /* Solve the upper triangular system Ux = z via back
       substitution. */ 

    b = x;
    for (i = N; i >= 1; --i) {
      x(i) = b(i)/A(i,i);
      for (j = 0; j <= i-1; ++j) { b(j) -= x(i)*A(j,i); }
    } /* endfor */
    x(0) = b(0)/A(0,0);

}

/********************************************************
 * Routine: Solve_Point_Jacobi                          *
 *                                                      *
 * Solves a dense system of equations of the form       *
 *                                                      *
 *               A x = b,                               *
 *                                                      *
 * where A is a dense (full) matrix, using the point    *
 * Jacobi iterative method.  The soluton vector is      *
 * returned in the column vector x.                     *
 *                                                      *
 ********************************************************/
void Solve_Point_Jacobi(DenseMatrix &A,
                        ColumnVector &b,
                        ColumnVector &x) {

    int i, j, iterations, N;
    double xnorm, x0norm;
    ColumnVector x0;

    /* Determine system size. */

    N = A.size(0)-1;

    /* Create temporary solution variable. */

    x0 = ColumnVector(N+1);

    /* Determine initial estimate for x. */

    for (i = 0; i <= N; ++i) {
       x0(i) = b(i)/A(i,i);
    } /* endfor */
    x0norm = max(x0.norm(), TOLER);

    /* Apply point Jacobi iterative procedure. */

    for (iterations = 0; iterations <= 10*(N+1)*(N+1); ++iterations) {
       for (i = 0; i <= N; ++i) {
	  x(i) = b(i);
          for (j = 0; j <= N; ++j) { if (j!=i) x(i) -= A(i,j)*x0(j); }
          x(i) = x(i)/A(i,i);
       } /* endfor */

       xnorm = norm(x-x0);

       if (fabs(xnorm/x0norm) < TEN*TOLER*TOLER) break;

       x0 = x;
    } /* endfor */

    /* Delete storage for temporary solution variable. */

    x0 = ColumnVector();

}

/********************************************************
 * Routine: Solve_Gauss_Seidel                          *
 *                                                      *
 * Solves a dense system of equations of the form       *
 *                                                      *
 *               A x = b,                               *
 *                                                      *
 * where A is a dense (full) matrix, using the Gauss-   *
 * Seidel iterative method.  The soluton vector is      *
 * returned in the column vector x.                     *
 *                                                      *
 ********************************************************/
void Solve_Gauss_Seidel(DenseMatrix &A,
                        ColumnVector &b,
                        ColumnVector &x) {

    int i, j, iterations, N;
    double xnorm, x0norm;
    ColumnVector x0;

    /* Determine system size. */

    N = A.size(0)-1;

    /* Create temporary solution variable. */

    x0 = ColumnVector(N+1);

    /* Determine initial estimate for x. */

    for (i = 0; i <= N; ++i) {
       x0(i) = b(i)/A(i,i);
       x(i) = x0(i);
    } /* endfor */
    x0norm = max(x0.norm(), TOLER);

    /* Apply point Gauss-Seidel iterative procedure. */

    for (iterations = 0; iterations <= 10*(N+1)*(N+1); ++iterations) {
       for (i = 0; i <= N; ++i) {
	  x(i) = b(i);
          for (j = 0; j <= N; ++j) { if (j!=i) x(i) -= A(i,j)*x(j); }
          x(i) = x(i)/A(i,i);
       } /* endfor */

       xnorm = norm(x-x0);

       if (fabs(xnorm/x0norm) < TEN*TOLER*TOLER) break;

       x0 = x;
    } /* endfor */

    /* Delete storage for temporary solution variable. */

    x0 = ColumnVector();

}

/********************************************************
 * Routine: Solve_Thomas_Algorithm                      *
 *                                                      *
 * Solves a tridiagonal system of equations of the form *
 *                                                      *
 *               A x = b,                               *
 *                                                      *
 * where A is a tridiagonal matrix, using the           *
 * Thomas algorithm, which is a version of Gaussian     *
 * elimination carried out via LU decomposition.  The   *
 * matrix A is put into upper triangular form and then  *
 * the solution x is found by back substitution.        *
 * Column vector b is also overwritten.                 *
 *                                                      *
 ********************************************************/
void Solve_Thomas_Algorithm(TriDiagonalMatrix &A,
                            ColumnVector &b,
                            ColumnVector &x) {

    int i, N;
    double ratio;

    /* Obtain the number of unknowns or size of system
       to be solved. */

    N = A.n-1;

    /* Establish the upper triangular matrix. */

    // b(0) = b(0); not needed.
    for ( i = 1 ; i <= N ; ++i ) {
       ratio=A.B(i)/A.D(i-1);
       A.D(i)=A.D(i)-ratio*A.A(i-1);
       b(i)=b(i)-ratio*b(i-1);
    } /* endfor */

    /* Perform the back substitution. */

    x(N)=b(N)/A.D(N);
    for ( i = N-1 ; i >= 0 ; --i ) {
       x(i)=(b(i)-A.A(i)*x(i+1))/A.D(i);
    } /* endfor */

}

/********************************************************
 * Routine: Solve_Point_Jacobi                          *
 *                                                      *
 * Solves a tridiagonal system of equations of the form *
 *                                                      *
 *               A x = b,                               *
 *                                                      *
 * where A is a tridiagonal matrix, using the point     *
 * Jacobi iterative method.  The soluton vector is      *
 * returned in the column vector x.                     *
 *                                                      *
 ********************************************************/
void Solve_Point_Jacobi(TriDiagonalMatrix &A,
                        ColumnVector &b,
                        ColumnVector &x) {

    int i, iterations, N;
    double xnorm, x0norm;
    ColumnVector x0;

    /* Determine system size. */

    N = A.n-1;

    /* Create temporary solution variable. */

    x0 = ColumnVector(N+1);

    /* Determine initial estimate for x. */

    for (i = 0; i <= N; ++i) {
       x0(i) = b(i)/A.D(i);
    } /* endfor */
    x0norm = max(x0.norm(), TOLER);

    /* Apply point Jacobi iterative procedure. */

    for (iterations = 0; iterations <= 10*(N+1)*(N+1); ++iterations) {
       for (i = 0; i <= N; ++i) {
          x(i) = b(i);
          if (i < N) x(i) -= A.A(i)*x0(i+1);
          if (i > 0) x(i) -= A.B(i)*x0(i-1);
          x(i) = x(i)/A.D(i);
       } /* endfor */

       xnorm = norm(x-x0);

       if (fabs(xnorm/x0norm) < TEN*TOLER*TOLER) break;

       x0 = x;
    } /* endfor */

    /* Delete storage for temporary solution variable. */

    x0 = ColumnVector();

}

/********************************************************
 * Routine: Solve_Gauss_Seidel                          *
 *                                                      *
 * Solves a tridiagonal system of equations of the form *
 *                                                      *
 *               A x = b,                               *
 *                                                      *
 * where A is a tridiagonal matrix, using the Gauss-    *
 * Seidel iterative method.  The soluton vector is      *
 * returned in the column vector x.                     *
 *                                                      *
 ********************************************************/
void Solve_Gauss_Seidel(TriDiagonalMatrix &A,
                        ColumnVector &b,
                        ColumnVector &x) {

    int i, iterations, N;
    double xnorm, x0norm;
    ColumnVector x0;

    /* Determine system size. */

    N = A.n-1;

    /* Create temporary solution variable. */

    x0 = ColumnVector(N+1);

    /* Determine initial estimate for x. */

    for (i = 0; i <= N; ++i) {
       x0(i) = b(i)/A.D(i);
       x(i) = x0(i);
    } /* endfor */
    x0norm = max(x0.norm(), TOLER);

    /* Apply Gauss Seidel iterative procedure. */

    for (iterations = 0; iterations <= 10*(N+1)*(N+1); ++iterations) {
       for (i = 0; i <= N; ++i) {
          x(i) = b(i);
          if (i < N) x(i) -= A.A(i)*x(i+1);
          if (i > 0) x(i) -= A.B(i)*x(i-1);
          x(i) = x(i)/A.D(i);
       } /* endfor */

       xnorm = norm(x-x0);

       if (fabs(xnorm/x0norm) < TEN*TOLER*TOLER) break;

       x0 = x;
    } /* endfor */

    /* Delete storage for temporary solution variable. */

    x0 = ColumnVector();

}
