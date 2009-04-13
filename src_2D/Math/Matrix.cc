/* Matrix.cc:  Subroutines for various matrix classes. */

/* Include matrix header file. */
#include "Matrix.h"
#include <lapackd.h>
#include <BPKIT.h>
#include "../Utilities/Utilities.h"    /* Include utilities header file. */
#include "../Utilities/EpsilonTol.h"   /* Include numerical tolerances header file. */

/*************************************************************
 * DenseMatrix -- Create storage for temp vector.            *
 *************************************************************/
RowVector    DenseMatrix::temp_RVec;

/*************************************************************
 * TriDiagonalMatrix -- Create storage for temp vector.      *
 *************************************************************/
RowVector    TriDiagonalMatrix::temp_RVec;


/*************************************************************
 *                   Fortran Wrappers                        *
 *************************************************************/
extern "C"
{
void F77NAME(dgeev) (char *, char *, integer *, doublereal *, integer *,
    doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *,
    doublereal *, integer *, integer *);

void F77NAME(dgetrf)(int *, int *, doublereal *, int *, int *, int *);

#ifdef _USE_ESSL_LIB
  int dgesvf  (int*, void *, int *, void *, int*, int*,
	       double *, int*, int*, double *, int*);
#endif
}

/*************************************************************
 * pseudo_inverse -- Compute the pseudo-inverse of a dense   *
 *                   matrix MxN and write it on top of the   *
 *                   initial matrix.                         *
 *************************************************************/
void DenseMatrix::pseudo_inverse_override(void){

#ifdef _USE_ESSL_LIB

  int iopt(12);		// require singular values, V and UtB to be computed. Sort the singular values in desceding order.
  integer NROW(size(0)), NCOL(size(1));
  integer LWORK( 10 * 2*NCOL + max(NROW, NCOL) );
  integer MN(min(NROW,NCOL));
  double *S, *WORK;
  WORK = new double[LWORK];	// work space
  S = new double[MN];		// array of singular values
  double LargestSV;		// the largest single value
  integer i,j,k;
  double tau(EpsilonTol::epsilon);

  // RHS matrix; ensure calculation of Ut (U transpose)
  DenseMatrix UtB(NROW,NROW); UtB.identity();

  /* Call Fortran subroutine from ESSL. 
     This routine calculates V matrix and Ut.
     The pseudo-inverse is the product of V, SIGMA-inverse (i.e. the inverse of the diagonal matrix), and Ut. */
  dgesvf(&iopt, &v_(0), &NROW, &UtB(0,0), &NROW, &NROW,
	 S, &NROW, &NCOL, WORK, &LWORK);

  // Store V matrix
  DenseMatrix V(*this);  

  // Get the reciprocal of each non-zero singular value of S
  LargestSV = S[0];
  if (LargestSV != 0){
    for (i=0; i<MN; ++i){
      if (S[i]/LargestSV > tau){
	S[i] = 1.0/S[i];
      } else {
	S[i] = 0.0;
      }
    }
  }
 
  // Resize the initial matrix to store the pseudo_inverse
  newsize(NCOL,NROW);

  // Compute the pseudo-inverse
  for (i=0; i<NCOL; ++i){
    for (j=0; j<NROW; ++j){
      operator()(i,j) = 0.0;
      for (k=0; k<NCOL; ++k){
	operator()(i,j) += V(i,k) * S[k] * UtB(k,j);
      }//endfor(k)
    }//endfor(j)
  }//endfor(i)

  // Deallocate locally allocated memory
  delete [] S; S = NULL;
  delete [] WORK; WORK = NULL;

  return;

#else

  // Get the SVD decomposition of matrix A using the 'dgesvd' Lapack subroutine

  // Set variables
  char JOBU('A'), JOBVT('A');
  integer NROW(size(0)), NCOL(size(1));
  integer LDA(NROW), LDU(NROW), LDVT(NCOL), INFO;
  double *S, *WORK;
  integer LWORK;
  integer MN(min(NROW,NCOL));
  integer i,j,k;
  double tau(1.0e-14);

  LWORK = 2*(max(1,max(3*MN+max(NROW,NCOL),5*MN)));
  WORK = new double[LWORK];
  S = new double[MN];
  double LargestSV;		// the largest single value

  DenseMatrix U(NROW,NROW), VT(NCOL,NCOL);

  /* Call Fortran subroutine */
  F77NAME(dgesvd)(&JOBU, &JOBVT, &NROW,
		  &NCOL,
		  &v_(0),
		  &LDA,
		  S, &U(0,0), &LDU, &VT(0,0), &LDVT, 
		  WORK, &LWORK, &INFO);

  if (INFO == 0){
    // Transpose the matrices
    U = U.transpose();
    VT = VT.transpose();

    // Get the reciprocal of each non-zero singular value of S
    LargestSV = S[0];
    if (LargestSV != 0){
      for (i=0; i<MN; ++i){
	if (S[i]/LargestSV > tau){
	  S[i] = 1.0/S[i];
	} else {
	  S[i] = 0.0;
	}
      }
    }

    // Resize the initial matrix to store the pseudo_inverse
    newsize(NCOL,NROW);

    // Compute the pseudo-inverse
    for (i=0; i<NCOL; ++i){
      for (j=0; j<NROW; ++j){
	operator()(i,j) = 0.0;
	for (k=0; k<NCOL; ++k){
	  operator()(i,j) += VT(i,k) * S[k] * U(k,j);
	}//endfor(k)
      }//endfor(j)
    }//endfor(i)

    // Deallocate locally allocated memory
    delete [] S; S = NULL;
    delete [] WORK; WORK = NULL;

    return;

  } else if (INFO < 0) {
    std::cout << "ERROR in Get_Matrix_PseudoInverse()! Error in the " << abs(INFO) << "-th argument.\n";
  } else {
    std::cout << "ERROR in Get_Matrix_PseudoInverse()! The DBDSQR did not converge for " << INFO
	      << " superdiagonals.\n";
  }

  // Deallocate locally allocated memory
  delete [] S; S = NULL;
  delete [] WORK; WORK = NULL;

#endif

}

/*************************************************************
 * eigenvalues -- Return a ColumnVector containng the        *
 *                eigenvalues of an NxN matrix.              *
 *             (Original Matrix is over-written!)            *
 *************************************************************/
ColumnVector DenseMatrix::eigenvalues_overwrite(void) {

  assert(size(0)==size(1));

  ColumnVector REAL(size(0)), IMAG(size(0));
  char JOBVL('N'), JOBVR('N');
  integer N(size(0));
  integer LDA(N), LDVL(N), LDVR(N);
  double dummy;
  double *WORK;
  integer LWORK;
  integer INFO;

  LWORK = 10*N;   //how big should this thing be?
  WORK = new double[LWORK];

  /* Call Fortran subroutine */
  F77NAME(dgeev)(&JOBVL, &JOBVR,
		 &N, &v_(0), &LDA,
		 &REAL(0), &IMAG(0),
		 &dummy, &LDVL, &dummy, &LDVR,
		 WORK, &LWORK, &INFO);


  delete [] WORK; WORK=NULL;

  return REAL;
}

/*************************************************************
 * inverse -- Return matrix inverse using LAPACK's           *
 *            dgetrf & dgetri functions.                     *
 *************************************************************/
void DenseMatrix::inverse_overwrite(void) {
  assert(size(0)==size(1));

  int N(size(0));
  int LDA(N);
  int *IPIV = new int[N];
  for(int i=0; i<N; ++i) {IPIV[i] = 0;}
  int IWORK(40*N); // how big should this be ?
  double *WORK = new double[IWORK];
  for(int i=0; i<N; ++i) {WORK[i] = 0.0;}
  int INFO(0);

  /* Call Fortran subroutines */
  F77NAME(dgetrf)(&N, &N, &v_(0), &LDA, IPIV, &INFO);
  if(INFO != 0) {
    delete [] IPIV; IPIV = NULL;
    delete [] WORK; WORK = NULL;
    bperror("Error in BPKIT:dgetrf",INFO);
  }
  F77NAME(dgetri)(&N, &v_(0), &LDA, IPIV,
		  WORK, &IWORK, &INFO);
  if(INFO != 0) {
    delete [] IPIV; IPIV = NULL;
    delete [] WORK; WORK = NULL;
    bperror("Error in BPKIT:dgetri",INFO);
  }

  delete [] IPIV; IPIV = NULL;
  delete [] WORK; WORK = NULL;

  return;
}
