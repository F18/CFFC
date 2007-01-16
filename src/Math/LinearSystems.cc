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

  int M, N, JB, ldiag;
  double CondNb;
  // A = Q*R*KQ_transpose
  //   DenseMatrix Q(A.size(0), A.size(0));
  //   DenseMatrix KQ(A.size(1), A.size(1));

  // Determine the value of tau
  double B_Uncertainty = 1.0e-10;
  double tau = B_Uncertainty*A.NormFro();

  M = A.size(0)-1;
  N = A.size(1)-1;
  JB = B.size(1)-1;

  ldiag = min(M,N);


  if (ldiag < N) {
    cout << "System not SOLVED!!!" << endl
	 << "The subroutine doesn't handle UNDERDETERMINED systems!" 
	 << endl;
    cout.flush();
    return;
  }

  double *H, *G;
  int *PermutationVector;

  H = new double [ldiag+1];	// helper array
  G = new double [ldiag+1];
  PermutationVector = new int[ldiag+1];	// stores the permutation indexes

  double Hmax;
  int K, lmax, pos;

  int index_K; // the index of the last diagonal element > tau
               // => A(index_K,index_K) > tau 
  double SM;

  K=0;

  require((B.size(0)== (M+1)) && (X.size(0)==(N+1)) && (B.size(1)==X.size(1)) && (ldiag>=0), 
	  "Solve_LS_Householder() ERROR! The dimension of the RHS, LHS and Solution vector are inconsistant.\n");

  for (int j=0; j<=ldiag; ++j){  // for each column of matrix A

    // Compute squared column lengths and find lmax

    ComputeSquaredColumnLength(j, A, lmax, Hmax, H);
    // Lmax has been determined

    // Do column interchanges is needed
    PermutationVector[j] = lmax;
    if (PermutationVector[j]!=j){
      A.permute_col(j,lmax);
      H[lmax] = H[j];
    }

    // Determine the largest element in magnitude between row j to m
    // for the column j

    pos = j;
    for (int i=j+1; i<=M; ++i){
      if (fabs(A(i,j)) > fabs(A(pos,j)))
	pos = i;
    } 

    // Do row interchange - based on the parameter pos
    if (pos != j){
      A.permute_row(j,pos);
      B.permute_row(j,pos);
    }

    // Compute the j-th transformation and apply it to A and B.
    HouseholderOrthogTransf_col(false,j,j+1,M,A,j,H[j],A,j+1,N);
    HouseholderOrthogTransf_col(true ,j,j+1,M,A,j,H[j],B,  0,JB);
  }


  // Determine the pseudorank, K, using the tolerance, tau.
  for (int i=0; (i<=ldiag)&&(fabs(A(i,i))>tau); ++i)
    K = i+1;    		// K=0 means that all abs(a_ii) <= tau

  index_K = K-1;

  // Compute the norms of the residual vectors
  for (int jb = 0; jb <= JB; ++jb){
    SM = 0.0;
    if (index_K > M) 
      RNorm(jb) = 0.0;
    else {
      for (int i=index_K+1; i<=M; ++i)
	SM += B(i,jb)*B(i,jb);
      RNorm(jb) = sqrt(SM);
    }
  } 

  // Special for pseudorank = 0
  if (K == 0){
    for (int jx = 0; jx <= X.size(1)-1; ++jx)
      for (int ix=0; ix <= X.size(0)-1; ++ix){
	X(ix,jx) = 0.0;
      }
    krank = K;
    delete [] H; H = NULL;
    delete [] G; G = NULL;
    delete [] PermutationVector; PermutationVector = NULL;
    return;
  }
 
  // If the pseudorank is less than N compute Householder decomposition
  // of the rows up to index_K
  if (index_K < N){
    for(int ii=index_K; ii>=0; --ii){
      HouseholderOrthogTransf_row(false,ii,K,N,A,ii,G[ii],A,0,ii-1);
    }
  }

 
  for (int jb=0; jb<=JB; ++jb){
    // Solve the K by K triangular system
    for (int l=index_K; l>=0; l--){
      if (l==(index_K)){
	X(l,jb) = B(l,jb)/A(l,l);
      }
      else {
	SM = 0.0;
	for (int j=l+1; j<=index_K; ++j)
	  SM += A(l,j)*X(j,jb);
	X(l,jb) = (B(l,jb) - SM)/A(l,l);
      }
    }

    // Complete computation of the solution vector
    DenseMatrix X_transposed(&X(0,0),X.dim(1),X.dim(0),MV_Matrix_::ref);
    if (index_K < N){
      for (int l=K; l<=N; ++l)
	X(l,jb) = 0.0;
       for (int i=0; i<=index_K; ++i){
	HouseholderOrthogTransf_row(true,i,index_K+1,N,A,i,G[i],X_transposed,jb,jb);
       }
    }

    // Re-order the solution vector to compensate for the column interchanges

    for (int j=ldiag; j>=0; --j){
      if (PermutationVector[j] != j){
	SM = X(PermutationVector[j],jb);
	X(PermutationVector[j],jb) = X(j,jb);
	X(j,jb) = SM;
      }
    }
  }

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
* C   ->  represents a matrix on which the determined Householder         *
*         transformation will be applied to. On exit, C contains the set  *
*         of transformed vectors.                                         *
* CLow -> the lower index of the column vector of matrix C which the      *
*         transformation is applied to                                    *
* CUp ->  the upper index of the column vector of matrix C which the      *
*         transformation is applied to                                    *
* Q ->    for explicit computation of the orthogonal transformation matrix*
*         Obs. For fast execution, this computation must be commented out.* 
**************************************************************************/
inline void HouseholderOrthogTransf_col (bool KnownTransformation, int & PivotElement, int StartZeroElem, 
					 int &MaxColumnIndex , DenseMatrix &V, int PivotVector,
					 double &UP, DenseMatrix &C, int CLow, int &CUp) {

  double SM, B;
  double VarPower2;
  double V_PivotElement;
  //  ColumnVector u(MaxColumnIndex+1);

  require((PivotElement>=0) && (PivotElement<StartZeroElem) && (StartZeroElem<=MaxColumnIndex),
	  "HouseholderOrthogTransf_col() ERROR! The relation between the indexes of the PivotVector is not satisfied!\n");

  V_PivotElement = fabs( V(PivotElement,PivotVector) );

  // If V_PivotElement <= 0.0 one gets an identity transformation
  if (V_PivotElement <= 0.0){
    //     cout << "HouseholderOrthogTransf_col()::Abort1" << endl;
    //     cout.flush();
    return;
  }

  switch(KnownTransformation){
  case false:
    /********************* Construct The Transformation *****************************/

    /************************************************************************************
     If this transformation is part of a subroutine which determines the biggest value in
     magnitude for the pivot element, the "for" loop below is useless.
     The V_PivotElement already has the maximum value in magnitude.

    If V_PivotElement hasn't the maximum value in magnitude (no pivoting), the first 
    "if" statement from above must be moved in the switch expression. 
    *************************************************************************************/
    //     for (int i=StartZeroElem; i<=MaxColumnIndex; ++i){
    //       // Determine the maximum absolute value of all the elements
    //       V_PivotElement = max(V_PivotElement, fabs(V(i,PivotVector)) );
    //     }

    VarPower2 = V(PivotElement,PivotVector)/V_PivotElement;
    // initialize sum SM
    SM = VarPower2*VarPower2;
    for(int i=StartZeroElem; i<=MaxColumnIndex; ++i){
      VarPower2 = V(i,PivotVector)/V_PivotElement;
      SM += VarPower2 * VarPower2;
    }
    SM = V_PivotElement*sqrt(SM);
    
    if ( V(PivotElement,PivotVector) > 0.0){
      -SM;
    }
    UP = V(PivotElement,PivotVector) - SM;
    V(PivotElement,PivotVector) = SM;
    

    /************************** Compute The Transformation Matrix ****************************/

    //     B = UP*V(PivotElement,PivotVector);		   // B must be negative here
    //     // if B=0.0 , return
    //     if (B >= 0.0){
    //       cout << "HouseholderOrthogTransf_col()::Abort2" << endl;
    //       cout.flush();
    //       return;
    //     }
    //     B = 1.0/B;
    
    //     // construct the matrix Q
    //     Q.identity();
    //     u.zero();
    //     u(PivotElement) = UP;
    //     for (int i=StartZeroElem; i<=MaxColumnIndex; ++i)
    //       u(i)= V(i,PivotVector);
    //     Q += B*(u*u.transpose());
    
  case true:
    /********************* Apply The Transformation *****************************/


    require( (CLow>=0) && (CUp<C.size(1)), "HouseholderOrthogTransf_col() ERROR! The transformation cannot be applied" \
	     " to the specified columns. The indeces are out of bounds!\n");

    if ( CLow > CUp){
      // there is no column for which the transformation be applied to 
      //       cout << "HouseholderOrthogTransf_col()::Abort3" << endl;
      //        //    cout << CLow <<"  "<< CUp << "  "<< C.size(1)-1 << endl;
      //       cout.flush();
      return;
    }

    B = UP*V(PivotElement,PivotVector);		   // B must be nonpositive here
    // if B=0.0 , return
    if (B >= 0.0){
      //       cout << "HouseholderOrthogTransf_col()::Abort4" << endl;
      //       cout.flush();
      return;
    }

    for (int j=CLow; j<=CUp; ++j){	   // for each column of matrix C between specified indexes
      SM = C(PivotElement,j)*UP;
      for (int i=StartZeroElem; i<=MaxColumnIndex; ++i){
	SM += C(i,j)*V(i,PivotVector);
      }
      if (SM==0.0)
	continue;
      else{
	SM /= B;
	C(PivotElement,j) += SM*UP;
	for(int i=StartZeroElem; i<=MaxColumnIndex; ++i)
	  C(i,j) += SM*V(i,PivotVector);
      }
    }
  }
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
inline void HouseholderOrthogTransf_row (bool KnownTransformation,int & PivotElement, int StartZeroElem,
					 int & MaxRowIndex, DenseMatrix &V, int PivotVector,
					 double &UP, DenseMatrix &C, int RLow, int RUp) {

  double SM, B;
  double VarPower2;
  double V_PivotElement;
  ColumnVector u(MaxRowIndex+1);

  require((PivotElement>=0) && (PivotElement<StartZeroElem) && (StartZeroElem<=MaxRowIndex),
	  "HouseholderOrthogTransf_row() ERROR! The relation between the indexes of the PivotVector is not satisfied!\n");

  V_PivotElement = fabs( V(PivotVector,PivotElement) );

  // If V_PivotElement <= 0.0 one gets an identity transformation
  if (V_PivotElement <= 0.0){
    //     cout << "HouseholderOrthogTransf_row()::Abort1" << endl;
    //     cout.flush();
    return;
  }

  switch(KnownTransformation){
  case false:
    /********************* Construct The Transformation *****************************/

    /************************************************************************************
     If this transformation is part of a subroutine which determines the biggest value in
     magnitude for the pivot element, the "for" loop below is useless.
     The V_PivotElement already has the maximum value in magnitude.

    If V_PivotElement hasn't the maximum value in magnitude (no pivoting), the first 
    "if" statement from above must be moved in the switch expression. 
    *************************************************************************************/
    //     for (int j=StartZeroElem; j<=MaxRowIndex; ++j){
    //       // Determine the maximum absolute value of all the elements
    //       V_PivotElement = max( V_PivotElement, fabs(V(PivotVector,j)) );
    //     }

    VarPower2 = V(PivotVector,PivotElement)/V_PivotElement;
    // initialize sum SM
    SM = VarPower2*VarPower2;
    for(int i=StartZeroElem; i<=MaxRowIndex; ++i){
      VarPower2 = V(PivotVector,i)/V_PivotElement;
      SM += VarPower2 * VarPower2;
    }
    SM = V_PivotElement*sqrt(SM);

    if (V(PivotVector,PivotElement) > 0.0) {
      -SM;
    }
    UP = V(PivotVector,PivotElement) - SM;
    V(PivotVector,PivotElement) = SM;

    /************************** Compute The Transformation Matrix ****************************/

    //     B = UP*V_PivotElement;    		   // B must be nonpositive here
    // 					   // if B=0.0 , return
    //     if (B >= 0.0){
    //       cout << "HouseholderOrthogTransf_row()::Abort2" << endl;
    //       cout.flush();
    //       return;
    //     }
    //     B = 1.0/B;

    //     // construct the matrix K
    //     K.identity();
    //     u.zero();
    //     u(PivotElement) = UP;
    //     for (int i=StartZeroElem; i<=MaxRowIndex; ++i)
    //       u(i)= V(PivotVector,i);
    //     K += B*(u*u.transpose());

  case true:
    /********************* Apply The Transformation *****************************/

    int MaxSizeC = C.size(0);

    require( (RLow>=0) && (RUp<MaxSizeC), "HouseholderOrthogTransf_row() ERROR! The transformation cannot be applied" \
	     " to the specified rows. The indexes are out of bounds!\n");

    if ( RLow > RUp){
      // there is no row for which the transformation be applied to 
      //       cout << "HouseholderOrthogTransf_row()::Abort3" << endl;
      //       //        cout << RLow <<"  "<< RUp << "  "<< C.size(0)-1 << endl;
      //       cout.flush();
      return;
    }

    B = UP*V(PivotVector,PivotElement);     // B must be nonpositive here
					   // if B=0.0 , return
    if (B >= 0.0){
      //       cout << "HouseholderOrthogTransf_row()::Abort4" << endl;
      //       cout.flush();
      return;
    }

    for (int i=RLow; i<=RUp; ++i){	   // for each row of matrix C between the specified indexes
      SM = C(i,PivotElement)*UP;
      for (int j=StartZeroElem; j<=MaxRowIndex; ++j){
	SM += C(i,j)*V(PivotVector,j);
      }
      if (SM==0.0)
	continue;
      else{
	SM /= B;
	C(i,PivotElement) += SM*UP;
	for(int j=StartZeroElem; j<=MaxRowIndex; ++j)
	  C(i,j) += SM*V(PivotVector,j);
      }
    }
  }
}

/********************************************************************
* Routine: ComputeSquaredColumnLength                               *
*                                                                   *
* Determines the column of matrix A that has the sum of squares of  *
* components in rows j through dim0_ the greatest.                  *
********************************************************************/
void ComputeSquaredColumnLength(int j, const DenseMatrix &A, int &lmax,
				double &Hmax, double *H){


  int M,N;
  double factor = 0.001;

  M = A.size(0)-1;
  N = A.size(1)-1;

  lmax = j;
  if (j==0){
    for (int l=0; l<=N; ++l){	// for each column
      H[l] = 0.0;
      for (int i=0; i<=M; ++i) // sum the elements of the rows
	H[l] += A(i,l)*A(i,l);
      if (H[l] > H[lmax])
	lmax = l;
    }
    Hmax = H[lmax];
  }
  else{
    for (int l=j; l<=N; ++l){
      H[l] -= A(j-1,l)*A(j-1,l);
      if (H[l] > H[lmax])
	lmax = l;
    }
    if ((Hmax + factor*H[lmax])-Hmax<=0){
      lmax = j;
      for (int l=j; l<=N; ++l){	// for each column
	H[l] = 0.0;
	for (int i=j; i<=M; ++i){ // sum the elements of the rows
	  H[l] += A(i,l)*A(i,l);
	}
	if (H[l] > H[lmax])
	  lmax = l;
      }
    }
    Hmax = H[lmax];
  }
}
