/*!\file LinearSystems.cc
  \brief Subroutines for various linear system classes. */

/* Include linear systems header file. */

#ifndef _LINEARSYSTEMS_INCLUDED
#include "LinearSystems.h"
#endif // _LINEARSYSTEMS_INCLUDED

#include "../Utilities/Utilities.h" // Include util functions & macros header file.

// Define macro to use Lapack library to solve least-squares problems
#define __LAPACK_LEAST_SQUARES__

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

void Solve_LS_Householder  (DenseMatrix &A,
                            ColumnVector &b,
                            ColumnVector &x,
			    int &krank, double &RNorm) {

  DenseMatrix B(&b(0),b.size(),1,MV_Matrix_::ref), X(&x(0),x.size(),1,MV_Matrix_::ref);
  ColumnVector Rnorm(1);

  // Solve the system
  Solve_LS_Householder (A, B, X, krank, Rnorm);

  // Transfer the residual  
  RNorm = Rnorm(0);

  return;
}

/********************************************************
 * Routine: Solve_LS_Householder                        *
 *                                                      *
 * Solves a linear least squares problem or a set of    *
 * linear least squares problems having the same LHS    *
 * matrix but different RHS vectors.                    *
 *                                                      *
 *               A X = B, (in the least squares sense)  *
 *                                                      *
 * where: A is a dense matrix MxN,                      *
 *        B is a dense matrix MxNB                      *
 *        X is a dense matrix NxNB                      *
 *          (its column vectors reprezents the solution *
 *           for each least squares problem)            *
 *                                                      *
 * The subrutine uses an absolute tolerance parameter:  *
 *  ->tau                                               *
 *                                                      *
 * Method:                                              *
 *  This subrutine first transforms the augmented matrix*
 * [A:B] to a matrix [R:C] using premultiplying         *
 * Householder transformations with column interchanges.*
 * All subdiagonal elements in the matrix R are zero and*
 * its diagonal elements satisfy                        *
 *   fabs(R[i,i]) >= fabs(R[i+1,i+1]), i=1,..,l-1,      *
 * where l=min{M,N}.                                    *
 *  The subrutine will set the pseudorank KRANK equal   *
 * to the number of diagonal elements of R exceeding    *
 * "tau" in magnitude.                                  *
 * Minimal length solution vectors x[j],                *
 * j=1,...,NB, will be computed for the problems defined*
 * by the first KRANK rows of [R:C].                    *
 * If the relative uncertainty in the data matrix B is  *
 * "rho", it is suggested that tau be set approximately *
 * equal to rho*||A||.                                  *
 * The residual for the least-squares problem is written*
 * in the RNorm vector.                                 *
 *                                                      *
 * The implemented method doesn't work with             *
 *          UNDERDETERMINED systems                     *
 *                                                      *
 * Some modifications may appear in the future          *
 ********************************************************/

void Solve_LS_Householder  (DenseMatrix &A,
                            DenseMatrix &B,
                            DenseMatrix &X,
			    int &krank,
			    ColumnVector &RNorm) {

  int M(A.size(0)-1), N(A.size(1)-1), JB(B.size(1)-1), ldiag;
  int K(0), lmax, pos;
  int index_K; // the index of the last diagonal element > tau
               // => A(index_K,index_K) > tau 
  double SM, Hmax;

  //  double CondNb;
  // A = Q*R*KQ_transpose
  //   DenseMatrix Q(A.size(0), A.size(0));
  //   DenseMatrix KQ(A.size(1), A.size(1));

  // Determine the value of tau
  static const double B_Uncertainty(1.0e-10);
  double tau(B_Uncertainty*A.NormFro());
  ldiag = min(M,N);

#ifndef __No_Checking__
  if (ldiag < N) {
    cout << "System not SOLVED!!!" << endl
	 << "The subroutine doesn't handle UNDERDETERMINED systems!" 
	 << endl;
    cout.flush();
    return;
  }
#endif

  double *H, *G;
  int *PermutationVector;

  H = new double [ldiag+1];	// helper array
  G = new double [ldiag+1];
  PermutationVector = new int[ldiag+1];	// stores the permutation indexes

  int i,j,jb,ix,jx,ii,l; 	// indices

  require((B.size(0)== (M+1)) && (X.size(0)==(N+1)) && (B.size(1)==X.size(1)) && (ldiag>=0), 
	  "Solve_LS_Householder() ERROR! The dimension of the RHS, LHS and Solution vector are inconsistant.\n");

  for (j=0; j<=ldiag; ++j){  // for each column of matrix A

    // Compute squared column lengths and find lmax
    ComputeSquaredColumnLength(j, A, lmax, Hmax, H, M, N);

    // Obs. lmax has been determined

    // Do column interchanges if needed
    PermutationVector[j] = lmax;
    if (PermutationVector[j]!=j){
      A.permute_col(j,lmax);
      H[lmax] = H[j];
    }

    // Determine the largest element in magnitude between row j to M
    // for the column j

    pos = j;
    for (i=j+1; i<=M; ++i){
      if (fabs(A(i,j)) > fabs(A(pos,j)))
	pos = i;
    } 

    // Do row interchange - based on the parameter pos
    if (pos != j){
      A.permute_row(j,pos);
      B.permute_row(j,pos);
    }

    // Compute the j-th transformation and apply it to A and B.
    HouseholderOrthogTransf_col(false,true,j,j+1,M,A,j,H[j],A,j+1,N);
    HouseholderOrthogTransf_col(true ,true,j,j+1,M,A,j,H[j],B,  0,JB);
  }

  // Determine the pseudorank, K, using the tolerance, tau.
  for (i=0; (i<=ldiag)&&(fabs(A(i,i))>tau); ++i)
    K = i+1;    		// K=0 means that all abs(a_ii) <= tau

  index_K = K-1;


// ************************ Compute the norms of the residual vectors ********************

  // only for speed
#if 0
  for (jb = 0; jb <= JB; ++jb){
    RNorm(jb) = 1.0;		       
  }
#endif

  // the true computation of RNorm
  for (jb = 0; jb <= JB; ++jb){
    SM = 0.0;
    if (index_K > M) 
      RNorm(jb) = 0.0;
    else {
      for (i=index_K+1; i<=M; ++i)
	SM += B(i,jb)*B(i,jb);
      RNorm(jb) = sqrt(SM);
    }
  } 
  //#endif // ************************ Compute the norms of the residual vectors ********************


  // Special for pseudorank = 0
  if (K == 0){
    for (jx = 0; jx <= (int)X.size(1)-1; ++jx)
      for (ix=0; ix <= (int)X.size(0)-1; ++ix){
	X(ix,jx) = 0.0;
      }
    krank = 0;
    delete [] H; H = NULL;
    delete [] G; G = NULL;
    delete [] PermutationVector; PermutationVector = NULL;
    return;
  }
 
  // If the pseudorank is less than N compute Householder decomposition
  // of the rows up to index_K
  if (index_K < N){
    for(ii=index_K; ii>=0; --ii){
      HouseholderOrthogTransf_row(false,true,ii,K,N,A,ii,G[ii],A,0,ii-1);
    }
  }


  for (jb=0; jb<=JB; ++jb){
    // Solve the K by K triangular system
    for (l=index_K; l>=0; l--){
      if (l==(index_K)){
	X(l,jb) = B(l,jb)/A(l,l);
      }
      else {
	SM = 0.0;
	for (j=l+1; j<=index_K; ++j)
	  SM += A(l,j)*X(j,jb);
	X(l,jb) = (B(l,jb) - SM)/A(l,l);
      }
    }

    // Complete computation of the solution vector
    if (index_K < N){
      DenseMatrix X_transposed(&X(0,0),X.dim(1),X.dim(0),MV_Matrix_::ref);
      for (l=K; l<=N; ++l){
	X(l,jb) = 0.0;
      }
      for (i=0; i<=index_K; ++i){
	HouseholderOrthogTransf_row(true,true,i,index_K+1,N,A,i,G[i],X_transposed,jb,jb);
      }
    }

    // Re-order the solution vector to compensate for the column interchanges
    for (j=ldiag; j>=0; --j){
      if (PermutationVector[j] != j){
	SM = X(PermutationVector[j],jb);
	X(PermutationVector[j],jb) = X(j,jb);
	X(j,jb) = SM;
      }
    }

  }//endfor (jb)

  // The solution vectors X are now computed
  krank = K;

  delete [] H; H = NULL;
  delete [] G; G = NULL;
  delete [] PermutationVector; PermutationVector = NULL;

  return;
}

/**************************************************************************
* Subroutine HouseholderOrthogonalTransformation_column                   *
*            HouseholderOrthogTransf_col                                  * 
*                                                                         *
* L.Ivan, UTIAS                                                           *
* The implementation follows the one introduced by C.L.Lawson and         *
* R.J.Hanson, Jet Propulsion Laboratory,                                  *
* in "Solving Least Squares Problems                                      *
*                                                                         *
* When "KnownTransformation == false", this subroutine constructs the     *
* orthogonal matrix "Q = I + (U*U_transpose)/B" and applies it to the     *
* specified column vectors of the dense matrix C.                         *
*                                                                         *
* When "KnownTransformation == true", the subroutine will only apply the  *
* transformation matrix to the specified column vectors of the dense      *
* matrix C.                                                               *
*                                                                         *
* PivotElement -> the index of the pivot element (counting of the column  *
*                                                  starts from 0)         *
* StartZeroElem -> if StartZeroElem <= ColumnSize, the transformation will*
*                  zero all the elements indexed from StartZeroElem       *
*                  through ColumnSize.                                    *
*               -> if StartZeroElem >= ColumnSize, the subroutine does an *
*                  identity transformation.                               *
* V    -> is the matrix that contains the pivot vector, which is specified*
*         by PivotVector.                                                 *
* PivotVector  -> the index of the column vector of matrix V.             *
* UP   -> contains quantities defining the vector U of the Householder    *
*         transformation. When the subroutine is called with              *
*         "KnownTransformation==true", this array should contain          *
*         quantities previously computed.                                 *
* C   ->  represents a matrix to which the determined Householder         *
*         transformation will be applied. On exit, C contains the set     *
*         of transformed vectors.                                         *
* CLow -> the lower index of the column vector of matrix C to which the   *
*         transformation is applied                                       *
* CUp ->  the upper index of the column vector of matrix C to which the   *
*         transformation is applied                                       *
* Q ->    for explicit computation of the orthogonal transformation matrix*
*         Obs. For fast execution, this computation must be commented out.* 
**************************************************************************/
void HouseholderOrthogTransf_col (bool KnownTransformation, bool MaxAbsValueDetermined,
				  int PivotElement, int StartZeroElem, 
				  int MaxColumnIndex , DenseMatrix &V, int PivotVector,
				  double &UP, DenseMatrix &C, int CLow, int CUp) {

  double SM, B;
  double MaxAbsValue, Inv_MaxAbsValue;
  int i,j;

  if( (PivotElement<0) || (PivotElement>=StartZeroElem) || (StartZeroElem>MaxColumnIndex) ){
    // These are situations when the transformation is equivalent to an indentity transformation
    return;
  }

  switch(KnownTransformation){

  case false:
    /*************************************************************************************
     ************************** Construct the Transformation *****************************
     *************************************************************************************/
    
    // Assigned initial value to MaxAbsValue
    MaxAbsValue = fabs( V(PivotElement,PivotVector) );

    if(MaxAbsValueDetermined == false){
      /************************************************************************************
       If this transformation is part of a subroutine which determines the largest absolute
       value for the pivot element, the "for" loop below is useless.
       (Pass "true" to "DeterminedMaxAbsValue" for faster execution)
       The "V_PivotElement" already has the maximum value in magnitude.
      *************************************************************************************/
      for (i=StartZeroElem; i<=MaxColumnIndex; ++i){
	// Determine the maximum absolute value of all elements between StartZeroElem and MaxColumnIndex
	MaxAbsValue = max(MaxAbsValue, fabs(V(i,PivotVector)) );
      }
    }

    if (MaxAbsValue == 0.0){	
      // The MaxAbsValue must be different than zero in order to apply the transformation
      // otherwise return the identity transformation
      return;
    }

    // For robustness, the computation of the square root of sum of squares can be made
    // resistant to underflow by using MaxAbsValue
    Inv_MaxAbsValue = 1.0/MaxAbsValue;

    // Compute the sum of squares
    SM = sqr(V(PivotElement,PivotVector)*Inv_MaxAbsValue);
    for(i=StartZeroElem; i<=MaxColumnIndex; ++i){
      SM += sqr(V(i,PivotVector)*Inv_MaxAbsValue);
    }
    // Take the square root of SM and adjust it with MaxAbsValue
    SM = MaxAbsValue * sqrt(SM);

    // Apply the correct sign to the sum based on the value of pivot element
    if ( V(PivotElement,PivotVector) >= 0.0){
      SM = -SM;
    }

    UP = V(PivotElement,PivotVector) - SM;

    // Store the value of SM in the location of V(PivotElement,PivotVector)
    V(PivotElement,PivotVector) = SM;

  case true:
    /*************************************************************************************
     **************************** Apply the Transformation *******************************
     *************************************************************************************/

    require( (CLow>=0) && (CUp<C.size(1)), "\nHouseholderOrthogTransf_col() ERROR! The transformation cannot be applied" \
	                                   " to the specified columns. The indexes are out of bounds!\n");

    // V(PivotElement,PivotVector) represents actually SM here !
    B = UP*V(PivotElement,PivotVector);

    if ( (CLow > CUp) || (B >= 0.0) ){
      // there is no column for which the transformation to be applied && B must be nonpositive
      return;
    }

    // Apply the transformation to the specified colomns of matrix C
    // Obs. Variable SM is reused
    for (j=CLow; j<=CUp; ++j){
      // for each column of matrix C between specified indexes
      // compute SM
      SM = C(PivotElement,j)*UP;
      for (i=StartZeroElem; i<=MaxColumnIndex; ++i){
	SM += C(i,j)*V(i,PivotVector);
      }

      // compute final values for the coefficients of matrix C
      if (SM==0.0)
	continue;
      else{
	SM /= B;
	C(PivotElement,j) += SM*UP;
	for(i=StartZeroElem; i<=MaxColumnIndex; ++i){
	  C(i,j) += SM*V(i,PivotVector);
	}
      }//endif
    }//endfor(j)
  }//endswitch
}

/**************************************************************************
* Subroutine HT12_row                                                     *
*                                                                         *
* L.Ivan, UTIAS                                                           *
* The implementation follows the one introduced by C.L.Lawson and         *
* R.J.Hanson, Jet Propulsion Laboratory,                                  *
* in "Solving Least Squares Problems                                      *
*                                                                         *
* Similar to HT12_col, the only difference being that the transformation  *
* is applied to rows instead of columns                                   *
**************************************************************************/
void HouseholderOrthogTransf_row (bool KnownTransformation, bool MaxAbsValueDetermined,
				  int PivotElement, int StartZeroElem,
				  int MaxRowIndex, DenseMatrix &V, int PivotVector,
				  double &UP, DenseMatrix &C, int RLow, int RUp) {
  
  double SM, B;
  double MaxAbsValue, Inv_MaxAbsValue;
  int i,j;

  if((PivotElement<0) || (PivotElement>=StartZeroElem) || (StartZeroElem>MaxRowIndex)){
    // These are situations when the transformation is equal to the indentity matrix
    return;
  }
  
  switch(KnownTransformation){
  case false:
    /*************************************************************************************
     ************************** Construct the Transformation *****************************
     *************************************************************************************/

    // Assigned initial value to MaxAbsValue
    MaxAbsValue = fabs( V(PivotVector,PivotElement) );

    if(MaxAbsValueDetermined == false){
      /************************************************************************************
       If this transformation is part of a subroutine which determines the biggest value in
       magnitude for the pivot element, the "for" loop below is useless.
       (Pass "true" to "DeterminedMaxAbsValue" for faster execution)
       The V_PivotElement already has the maximum value in magnitude.
      *************************************************************************************/
      for (j=StartZeroElem; j<=MaxRowIndex; ++j){
	// Determine the maximum absolute value of all the elements between StartZeroElem and MaxRowIndex
	MaxAbsValue = max( MaxAbsValue, fabs(V(PivotVector,j)) );
      }
    }

    if (MaxAbsValue == 0.0){	
      // The MaxAbsValue must be different than zero in order to apply the transformation
      // otherwise return the identity transformation
      return;
    }

    // For robustness, the computation of the square root of sum of squares can be made
    // resistant to underflow by using MaxAbsValue
    Inv_MaxAbsValue = 1.0/MaxAbsValue;

    // Compute the sum of squares
    SM = sqr(V(PivotVector,PivotElement)*Inv_MaxAbsValue);
    for(i=StartZeroElem; i<=MaxRowIndex; ++i){
      SM += sqr(V(PivotVector,i)*Inv_MaxAbsValue);
    }
    // Take the square root of SM and adjust it with MaxAbsValue
    SM = MaxAbsValue * sqrt(SM);

    // Apply the correct sign to the sum based on the value of pivot element
    if (V(PivotVector,PivotElement) >= 0.0) {
      SM = -SM;
    }
      
    UP = V(PivotVector,PivotElement) - SM;

    // Store the value of SM in the location of V(PivotVector,PivotElement)
    V(PivotVector,PivotElement) = SM;

  case true:
    /*************************************************************************************
     **************************** Apply the Transformation *******************************
     *************************************************************************************/

#if 0
    require( (RLow>=0) && (RUp<C.size(0)), "HouseholderOrthogTransf_row() ERROR! The transformation cannot be applied" \
	     " to the specified rows. The indexes are out of bounds!\n");
#endif

    // V(PivotVector,PivotElement) represents actually SM here !
    B = UP*V(PivotVector,PivotElement);

    if ( (RLow > RUp) || (B >= 0.0)){
      // there is no row for which the transformation to be applied && B must be nonpositive
      return;
    }

    // Apply the transformation to the specified colomns of matrix C
    // Obs. Variable SM is reused
    for (i=RLow; i<=RUp; ++i){
      // for each row of matrix C between the specified indexes
      SM = C(i,PivotElement)*UP;
      for (j=StartZeroElem; j<=MaxRowIndex; ++j){
	SM += C(i,j)*V(PivotVector,j);
      }

      // compute final values for the coefficients of matrix C
      if (SM==0.0)
	continue;
      else{
	SM /= B;
	C(i,PivotElement) += SM*UP;
	for(j=StartZeroElem; j<=MaxRowIndex; ++j){
	  C(i,j) += SM*V(PivotVector,j);
	}
      }//endif
    }//endfor(i)
  }//endswitch
}

/********************************************************************
* Routine: ComputeSquaredColumnLength                               *
*                                                                   *
* Determines the column of matrix A that has the sum of squares of  *
* components in rows j through dim0_ the greatest.                  *
********************************************************************/

void ComputeSquaredColumnLength(const int &j, const DenseMatrix &A, int &lmax,
				double &Hmax, double *H, const int &M, const int &N){


  static double factor = 0.001;
  int l,i;

  lmax = j;
  if (j==0){
    for (l=0; l<=N; ++l){	// for each column
      H[l] = 0.0;
      for (i=0; i<=M; ++i) // sum the elements of the rows
	H[l] += A(i,l)*A(i,l);
      if (H[l] > H[lmax])
	lmax = l;
    }
    Hmax = H[lmax];
  }
  else{
    for (l=j; l<=N; ++l){
      H[l] -= A(j-1,l)*A(j-1,l);
      if (H[l] > H[lmax])
	lmax = l;
    }
    if ((Hmax + factor*H[lmax])-Hmax<=0){
      lmax = j;
      for (l=j; l<=N; ++l){	// for each column
	H[l] = 0.0;
	for (i=j; i<=M; ++i){ // sum the elements of the rows
	  H[l] += A(i,l)*A(i,l);
	}
	if (H[l] > H[lmax])
	  lmax = l;
      }
    }
    Hmax = H[lmax];
  }
}


/********************************************************
 * Routine: Solve_Constrained_LS_Householder            *
 *                                                      *
 * Solves a linear constrained least squares problem    *
 * or a set of linear least squares problems having the *
 * same LHS matrix but different RHS vectors.           *
 * This is called a linear-equality-constrained least   *
 * squares problem (LSE).                               *
 * The first (NumberOfConstraints) equations represent  *
 * constraints and therefore must be satisfied exactly. *
 *                                                      *
 * The routine solves this problem on the assumption    *
 * that the matrix of the constrains has a full row     *
 * rank (no linear combination of constraints is        *
 * specified)!!!!                                       *
 *                                                      *
 *               A X = B, (in the least squares sense)  *
 *                                                      *
 * where: A is a dense matrix MxN,                      *
 *        B is a dense matrix MxNB                      *
 *        X is a dense matrix NxNB                      *
 *          (its column vectors reprezents the solution *
 *           for each least squares problem)            *
 *                                                      *
 * Method: Gauss elimination of the constaints followed *
 *         by solving a LS problem with Householder     *
 *         transformation.                              *
 *         Total pivoting procedure is carried out in   *
 *         the Gauss elimination step. The pivot is     *
 *         selected from the exactly satisfied equations*
 ********************************************************/

void Solve_Constrained_LS_Householder  (DenseMatrix &A,
					DenseMatrix &B,
					DenseMatrix &X,
					const int NumberOfConstraints) {

  int M(A.size(0)-1);		// maxim row indice of matrix A
  int N(A.size(1)-1);		// maxim column indice of matrix A
  int JB(B.size(1)-1);		// maxim column indice of matrix B

  int stage,row,col;
  int PivotR, PivotC;		// the row and the column indices of the pivot

  // Memory allocation for the Least-Squares problem
  DenseMatrix A_LS(M+1-NumberOfConstraints,N+1-NumberOfConstraints), B_LS(M+1-NumberOfConstraints,JB+1);
  int krank;

  int *PermutationVector;
  PermutationVector = new int[NumberOfConstraints];	// stores the permutation indexes

  // Apply Gauss elimination for the first NumberOfConstraints rows
  for (stage=0; stage<NumberOfConstraints; ++stage){

    // determine the coefficient with the maximum absolute value for the current stage
    // Obs. This coefficient is going to be the pivot
    PivotR = stage;		// initialize PivotR
    PivotC = stage;             // initialize PivotC
    for (col=stage; col<=N; ++col){
      for (row=stage; row<NumberOfConstraints; ++row){
	if ( fabs(A(row,col)) > fabs(A(PivotR,PivotC)) ){
	  PivotR = row;
	  PivotC = col;
	}
      }
    }

    // do row interchanges if needed
    if (PivotR != stage){
      A.permute_row(stage,PivotR);
      B.permute_row(stage,PivotR);
    }

    // do column interchanges if needed
    PermutationVector[stage] = PivotC;
    if (PivotC!=stage){
      A.permute_col(stage,PivotC);
    }


    // modify the pivot's equations such that the final pivot value is 1
    for (col=stage+1; col<=N; ++col){
      A(stage,col) /= A(stage,stage);
    }//endfor (col of A)
    for (col=0; col<=JB; ++col){ // modify the entries of the free term
      B(stage,col) /= A(stage,stage);
    }//endfor (col of B)
    A(stage,stage) = 1.0; 	// the diagonal element is 1

    // Gauss elimination
    for (row=stage+1; row<=M; ++row){ //for each row below the pivot's row
      for(col=stage+1; col<=N; ++col){ //for each column after the pivot's column
	A(row,col) -= A(row,stage)*A(stage,col); // update the value at the current elimination stage
      }//endfor (col of A)

      for(col=0; col<=JB; ++col){ //for each column of B
	B(row,col) -= A(row,stage)*B(stage,col); // update the value at the current elimination stage
      }//endfor (col of B)

    }//endfor (row)
    
  }//endfor (stage)

  if ( (int)A.size(1) > NumberOfConstraints){

    // copy the approximate part of the linear system
    for(row=NumberOfConstraints; row<=M; ++row){
      for(col=NumberOfConstraints; col<=N; ++col){
	A_LS(row-NumberOfConstraints,col-NumberOfConstraints) = A(row,col);
      }
    
      for(col=0; col<=JB; ++col){
	B_LS(row-NumberOfConstraints,col) = B(row,col);
      }
    }

    // Solve the Least Squares problem and copy solution to the X matrix
    /*******************************************************************/
#ifdef __LAPACK_LEAST_SQUARES__

    // Solve the least-squares system with Lapack subroutine
    Solve_LS_Householder_F77(A_LS, B_LS, krank, JB+1, M+1-NumberOfConstraints,N+1-NumberOfConstraints);

    // copy the solution of the least-squares problem (stored in B_LS) into X
    for(col=0; col<=JB; ++col){
      for(row=0; row<=N-NumberOfConstraints; ++row){
	X(row+NumberOfConstraints,col) = B_LS(row,col);
      }
    }//endfor

#else

    /* Solve the overdetermined linear system of equations using a least-squares procedure written by L. Ivan */
    /**********************************************************************************************************/
    DenseMatrix X_LS(N+1-NumberOfConstraints,JB+1);
    ColumnVector Rnorm(JB+1);

    Solve_LS_Householder(A_LS,B_LS,X_LS,krank,Rnorm);

    // copy X_LS into X
    for(col=0; col<=JB; ++col){
      for(row=0; row<=N-NumberOfConstraints; ++row){
	X(row+NumberOfConstraints,col) = X_LS(row,col);
      }
    }//endfor
#endif // __LAPACK_LEAST_SQUARES__

  } else {  // This part is for exactly solved systems of equations
    // set some derivatives to zero
    for(col=0; col<=JB; ++col){
      for(row=0; row<=N-NumberOfConstraints; ++row){
	X(row+NumberOfConstraints,col) = 0.0;
      }
    }//endfor
  }

  // back substitution for the constrained equations
  for(row=NumberOfConstraints-1; row>=0; --row){
    for(col=0; col<=JB; ++col){	// for each column of B
      A(row,0) = 0.0;		// reset the sum (Obs. This matrix coefficient is Zero anyway)
      for(stage=row+1; stage<=N; ++stage){
	A(row,0) += A(row,stage)*X(stage,col);
      }

      X(row,col) = B(row,col) - A(row,0); // compute final unknown
    }
  }

  // Re-order the solution vector to compensate for the column interchanges
  for (col=0; col<=JB; ++col){ // for each column of B
    for (stage=NumberOfConstraints-1; stage>=0; --stage){
      if (PermutationVector[stage] != stage){ // check the permutation vector at every stage
	A(0,0) = X(PermutationVector[stage],col); // A(0,0) is used as transfer variable
	X(PermutationVector[stage],col) = X(stage,col);
	X(stage,col) = A(0,0);
      }//endif
    }//endfor (stage)
  }//endfor (col)


  // Deallocate memory
  delete [] PermutationVector; PermutationVector = NULL;
}


