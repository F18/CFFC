/* Matrix.cc:  Subroutines for various matrix classes. */

/* Include matrix header file. */
#include "Matrix.h"
#include <lapackd.h>
#include <BPKIT.h>

/*************************************************************
 * DenseMatrix -- Create storage for temp vector.            *
 *************************************************************/
RowVector    DenseMatrix::temp_RVec;

/*************************************************************
 * TriDiagonalMatrix -- Create storage for temp vector.      *
 *************************************************************/
RowVector    TriDiagonalMatrix::temp_RVec;

/*************************************************************
 * pseudo_inverse -- Compute the pseudo-inverse of a dense   *
 *                   matrix MxN and write it on top of the   *
 *                   initial matrix.                         *
 *************************************************************/
void DenseMatrix::pseudo_inverse_override(void){

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

}
