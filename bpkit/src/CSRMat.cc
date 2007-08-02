//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#include <cassert>

using namespace std;

#include "CSRMat.h"
#include "LocalPrecon.h"
#include "BPKIT.h"
#include "lapackd.h"
#include "spblas.h"

// get resources
int CSRMat_APINV0::transpose = BpGetInt(APINV0_TRANSPOSE, 1);
int CSRMat_APINV_BANDED::transpose = BpGetInt(APINV_BANDED_TRANSPOSE, 1);

int    CSRMat_RILUK::growth      = BpGetInt(RILUK_GROWTH, 2);
double CSRMat_RILUK::stab_thresh = BpGetDouble(RILUK_THRESH, 0.0);
double CSRMat_RILUK::stab_rel    = BpGetDouble(RILUK_REL, 1.0);
double CSRMat_RILUK::stab_abs    = BpGetDouble(RILUK_ABS, 0.0);
int    CSRMat_RILUK::stat        = BpGetInt(RILUK_STAT, 0);

void CSRMat::SetToZero(int r, int c)
{
    // always deallocate space
    delete [] a;
    delete [] ja;
    delete [] ia;

    nrow = r;
    ncol = c;
    nnz = 0;
    a = NULL;
    ja = NULL;
    ia = new int[nrow+1];

    for (int i=0; i<=nrow; i++)
        ia[i] = 0;
}

void CSRMat::MatCopy(const LocalMat& A_)
{
    CSRMat *A = (CSRMat *) &A_;

    nrow = A->numrow();
    ncol = A->numcol();
    nnz  = A->numnz();
    delete [] a;
    delete [] ja;
    delete [] ia;
    a = new double[nnz];
    ja = new int[nnz];
    ia = new int[nrow+1];

    int i;
    for (i=0; i<nnz; i++)
    {
        a[i] = A->a[i];
        ja[i] = A->ja[i];
    }

    for (i=0; i<=nrow; i++)
        ia[i] = A->ia[i];
}

// C = alpha A B + beta C

void CSRMat::Mat_Vec_Mult(const BlockVec& B, BlockVec& C,
    double alpha, double beta) const
{
    assert(B.v != C.v);

    int descra[9];
    descra[0] = 0;
    descra[1] = 0;
    descra[2] = 0;

    int M = nrow;
    int N = C.dim1;
    int K = ncol;
    int LDB = B.dim0;
    int LDC = C.dim0;

    // Sparselib++ version of dcsrmm does not actually use workspace.
    // When this is changed, do not dynamically allocate space repeatedly
    // for this function.

    // non-FORTRAN with 1-based indexing
    F77NAME(dcsrmm) (0, M, N, K, alpha,
        descra, a, ja, ia, B.v+1, LDB,
        beta, C.v, LDC, NULL, 0);
}

// C = alpha A' B + beta C

void CSRMat::Mat_Trans_Vec_Mult(const BlockVec& B, BlockVec& C,
    double alpha, double beta) const
{
    assert(B.v != C.v);

    int descra[9];
    descra[0] = 0;
    descra[1] = 0;
    descra[2] = 0;

    int M = ncol;
    int N = C.dim1;
    int K = nrow;
    int LDB = B.dim0;
    int LDC = C.dim0;

    // non-FORTRAN with 1-based indexing
    F77NAME(dcsrmm) (1, M, N, K, alpha,
        descra, a, ja, ia, B.v, LDB,
        beta, C.v+1, LDC, (double *)1, 0);
}

void CSRMat::Mat_Vec_Solve(const BlockVec&, BlockVec&) const
{
    bperror("CSRMat::Mat_Vec_Solve: to solve with a local matrix,\n"
            "a local preconditioner should be used.\n", 0);
}

void CSRMat::Mat_Trans_Vec_Solve(const BlockVec&, BlockVec&) const
{
    bperror("CSRMat::Mat_Trans_Vec_Solve: to solve with a local matrix,\n"
            "a local preconditioner should be used.\n", 0);
}

LocalMat *CSRMat::CreateInv(LocalPrecon& local_precon) const
{
    switch (local_precon.name)
    {
    case LP_INVERSE:
        return new CSRMat_INVERSE(*this);
    case LP_LU:
        return new CSRMat_LU(*this);
    case LP_RILUK:
        return new CSRMat_RILUK(*this, local_precon.iarg1, local_precon.darg1);
    case LP_ILUT:
        return new CSRMat_ILUT(*this, local_precon.iarg1, local_precon.darg1);
    case LP_APINV0:
        return new CSRMat_APINV0(*this);
    case LP_APINVS:
        return new CSRMat_APINVS(*this, local_precon.iarg1, local_precon.iarg2);
    case LP_DIAG:
        return new CSRMat_DIAG(*this);
    case LP_TRIDIAG:
        return new CSRMat_TRIDIAG(*this);
    case LP_SOR:
        return new CSRMat_SOR(*this, local_precon.darg1, local_precon.iarg1);
    case LP_SSOR:
        return new CSRMat_SSOR(*this, local_precon.darg1, local_precon.iarg1);
    case LP_GMRES:
        return new CSRMat_GMRES(*this, local_precon.iarg1, local_precon.darg1);
    case LP_APINV_BANDED:
        return new CSRMat_APINV_BANDED(*this, local_precon.iarg1);
    case LP_APINV_TRUNC:
        return new CSRMat_APINV_TRUNC(*this, local_precon.iarg1);
    default:
        bperror("The local preconditioner you have chosen is not available\n"
                "for the type of blocks (CSR) in the block matrix.", 0);
    }
    return NULL; // suppress compiler warning
}

// this simply calls RILUK for now

CSRMat_LU::CSRMat_LU(const CSRMat& A)
{
    nrow = A.nrow;
    ncol = A.ncol;

    assert (nrow == ncol);

    assert (nrow != 0);
    assert (A.a != NULL);

    int nzl, nzu;

    allocate_ilu(nrow, nrow, &nzl, &nzu, A.ia, A.ja, &ial, &jal, &iau, &jau, 0);

    symbolic_ilu(nrow, nrow, &nzl, &nzu, A.ia, A.ja, ial, jal, iau, jau);

    // allocate for values
    nnz = nzl + nzu;
    al = new double[nzl];
    au = new double[nzu];
 
    numeric_ilu(0.0, nrow, A.ia, A.ja, A.a, ial, jal, al, iau, jau, au,
        0.0, 1.0, 0.0);
}

CSRMat_LU::~CSRMat_LU()
{
    delete [] al;
    delete [] jal;
    delete [] ial;
    delete [] au;
    delete [] jau;
    delete [] iau;
}

void CSRMat_LU::Mat_Vec_Solve(const BlockVec& B, BlockVec& X) const
{
  // dcsrsm allows solves in place

  int M = X.size0;
  int N = X.dim1;
  double *work = NULL; // no actual workspace

  int descra[9];

  descra[0] = 0;

  // lower unit
  descra[1] = 1;
  descra[2] = 1;

  // non-FORTRAN with 1-based indexing
  F77NAME(dcsrsm) (0, M, N, 1, NULL, 1.0,
           descra, al, jal, ial,
           B.v, M, 0.0, X.v+1, M,
           work, 0);

  // upper diag
  descra[1] = 2;
  descra[2] = 0;

  // non-FORTRAN with 1-based indexing
  F77NAME(dcsrsm) (0, M, N, 1, NULL, 1.0,
           descra, au, jau, iau,
           X.v, M, 0.0, X.v+1, M,
           work, 0);
}

CSRMat_RILUK::CSRMat_RILUK(const CSRMat& A, int levfill, double omega)
{
    nrow = A.nrow;
    ncol = A.ncol;

    assert (nrow == ncol);

    assert (nrow != 0);
    assert (A.a != NULL);

    int nzl, nzu;

    allocate_ilu(levfill, nrow, &nzl, &nzu, A.ia, A.ja, &ial, &jal, &iau, &jau,
        growth);

    symbolic_ilu(levfill, nrow, &nzl, &nzu, A.ia, A.ja, ial, jal, iau, jau);

    // allocate for values
    nnz = nzl + nzu;
    al = new double[nzl];
    au = new double[nzu];
 
    numeric_ilu(omega, nrow, A.ia, A.ja, A.a, ial, jal, al, iau, jau, au,
        stab_thresh, stab_rel, stab_abs);

    if (stat)
        ilu_stat(nrow, A.ia, A.ja, A.a, ial, jal, al, iau, jau, au);
}

CSRMat_RILUK::~CSRMat_RILUK()
{
    delete [] al;
    delete [] jal;
    delete [] ial;
    delete [] au;
    delete [] jau;
    delete [] iau;
}

void CSRMat_RILUK::Mat_Vec_Solve(const BlockVec& B, BlockVec& X) const
{
  // dcsrsm allows solves in place

  int M = X.size0;
  int N = X.dim1;
  double *work = NULL; // no actual workspace

  int descra[9];

  descra[0] = 0;

  // lower unit
  descra[1] = 1;
  descra[2] = 1;

  // non-FORTRAN with 1-based indexing
  F77NAME(dcsrsm) (0, M, N, 1, NULL, 1.0,
           descra, al, jal, ial,
           B.v, M, 0.0, X.v+1, M,
           work, 0);

  // upper diag
  descra[1] = 2;
  descra[2] = 0;

  // non-FORTRAN with 1-based indexing
  F77NAME(dcsrsm) (0, M, N, 1, NULL, 1.0,
           descra, au, jau, iau,
           X.v, M, 0.0, X.v+1, M,
           work, 0);
}

// exact inverse with pivoting, implemented with LAPACK

CSRMat_INVERSE::CSRMat_INVERSE(const CSRMat& A)
{
    int i, j, k;

    nrow = A.nrow;
    ncol = A.ncol;
    assert (nrow == ncol);
    assert (A.a != NULL);

    // copy sparse matrix A into dense storage array
    a = new double[nrow*ncol];
    ja = new int[nrow*ncol];
    ia = new int[nrow+1];
    nnz = nrow*ncol;

    double *p = a;
    for (i=0; i<nrow*ncol; i++)
        *p++ = 0.0;

    for (i=0; i<nrow; i++)
        for (j=A.ia[i]; j<A.ia[i+1]; j++)
	    a[A.ja[j] + i*nrow] = A.a[j];

    integer M = nrow;
    integer N = ncol;
    integer LDA = nrow;
    integer INFO;
    integer *ipiv = new integer[nrow];
    integer LWORK = N*N;
    double *work = new double[LWORK];

    // LU factorize
    F77NAME(dgetrf)(&M, &N, a, &LDA, ipiv, &INFO);
    if (INFO != 0)
        bperror("CSRMat_INVERSE: dgetrf returned", INFO);

    // compute inverse
    F77NAME(dgetri)(&N, a, &LDA, ipiv, work, &LWORK, &INFO);
    if (INFO != 0)
        bperror("CSRMat_INVERSE: dgetri returned", INFO);

    delete [] ipiv;
    delete [] work;

    // assign ia and ja
    k = 0;
    for (i=0; i<nrow; i++)
    {
        ia[i] = i*ncol;
        for (j=0; j<nrow; j++)
        {
            ja[k] = j;
            k++;
        }
    }
    ia[nrow] = nrow*ncol;
}

CSRMat_APINV_TRUNC::CSRMat_APINV_TRUNC(const CSRMat& A, int semi)
{
    int i, j, k;

    nrow = A.nrow;
    ncol = A.ncol;
    assert (nrow == ncol);
    assert (A.a != NULL);

    // copy sparse matrix A into dense storage array
    double *b = new double[nrow*ncol];
    double *p = b;
    for (i=0; i<nrow*ncol; i++)
        *p++ = 0.0;

    for (i=0; i<nrow; i++)
        for (j=A.ia[i]; j<A.ia[i+1]; j++)
	    b[A.ja[j] + i*nrow] = A.a[j];

    integer M = nrow;
    integer N = ncol;
    integer LDA = nrow;
    integer INFO;
    integer *ipiv = new integer[nrow];
    integer LWORK = N*N;
    double *work = new double[LWORK];

    // LU factorize
    F77NAME(dgetrf)(&M, &N, b, &LDA, ipiv, &INFO);
    if (INFO != 0)
        bperror("CSRMat_APINV_TRUNC: dgetrf returned", INFO);

    // compute inverse
    F77NAME(dgetri)(&N, b, &LDA, ipiv, work, &LWORK, &INFO);
    if (INFO != 0)
        bperror("CSRMat_INVERSE: dgetri returned", INFO);

    delete [] ipiv;
    delete [] work;

    assert(nrow >= semi);
    nnz = nrow*(2*semi+1) - semi*(semi+1);
     a = new double[nnz];
    ja = new int[nnz];
    ia = new int[nrow+1];

    int col1, col2;
    k = 0;
    ia[0] = 0;
    for (i=0; i<nrow; i++)
    {
	col1 = MAX(i-semi,0);
	col2 = MIN(i+semi,nrow-1);

	for (j=col1; j<=col2; j++)
	{
	    a[k] = b[j + i*nrow];
	    ja[k] = j;
	    k++;
	}
	ia[i+1] = k;
    }
    assert(k == nnz);

    delete [] b;
}

CSRMat_APINV0::CSRMat_APINV0(const CSRMat& A)
{
    assert(A.nrow == A.ncol);
    int i, nind, colnnz, col, ii, jrow, j;
    double *ahat, *p;

    int n = nrow = ncol = A.nrow;
    nnz = A.nnz;

    // matrix for which right approx inverse is required
    double *ao;
    int *jao;
    int *iao;

    // transpose matrix if want right apinv
    if (transpose)
    {
        ao = new double[nnz];
        jao = new int[nnz];
        iao = new int[n+1];
        csrtrans(n, A.a, A.ja, A.ia, ao, jao, iao);
    }
    else
    {
	ao = A.a;
	jao = A.ja;
	iao = A.ia;
    }

    // allocate space for approximate inverse
    a = new double[nnz];
    ja = new int[nnz];
    ia = new int[n+1];

    // copy structure of ia and ja
    for (i=0; i<=n; i++)
        ia[i] = iao[i];
    for (i=0; i<ia[n]; i++)
        ja[i] = jao[i];

    // workspaces
    int *ind = new int[n];
    int *iwk = new int[n];
    for (i=0; i<n; i++)
        iwk[i] = -1;

    // b vector and workspace for dgels
    double *b = new double[n];
    integer LWORK = 10000;
    double *WORK = new double[LWORK];

    // main loop over columns
    for (i=0; i<n; i++)
    {
        // determine all nonzeros row indices for 
        // columns indicated in i-th column of A
        nind = 0;
        colnnz = ia[i+1] - ia[i];

        for (j=ia[i]; j<ia[i+1]; j++)
        {
            // scan column
            col = ja[j];
            for (ii=iao[col]; ii<iao[col+1]; ii++)
            {
                jrow = jao[ii];
                if (iwk[jrow] == -1)
                {
                    iwk[jrow] = nind;
                    ind[nind++] = jrow;
                }
            }
        }

        assert (nind > 0);
        assert (colnnz > 0);

        // form least-squares matrix
        ahat = new double[nind*colnnz];
        for (j=0; j<nind*colnnz; j++)
            ahat[j] = 0.0;

        p = ahat;
        for (j=ia[i]; j<ia[i+1]; j++)
        {
            col = jao[j];
            for (ii=iao[col]; ii<iao[col+1]; ii++)
            {
                jrow = ja[ii];
                p[iwk[jrow]] = ao[ii];
            }
            p = p + nind;
        }

        // form b
        for (j=0; j<nind; j++)
            b[j] = 0.0;
        if (iwk[i] >= 0)
            b[iwk[i]] = 1.0;

#if DEBUG
        if (iwk[i] < 0)
           cerr << "APINV0: no diag in col " << i << endl;
#endif

        // reset workspace
        for (j=0; j<nind; j++)
            iwk[ind[j]] = -1;

        // solve least squares system
        char TRANS = 'N';
        integer NROW = nind;
        integer NCOL = colnnz;
        assert (NROW >= NCOL);
        integer ONE = 1;
        integer INFO;

        F77NAME(dgels)(&TRANS, &NROW, &NCOL, &ONE, ahat, &NROW,
            b, &NROW, WORK, &LWORK, &INFO);
        if (INFO != 0)
            bperror("CSRMat_APINV0: dgels returned", INFO);

        delete [] ahat;

        // copy result into matrix
        j = 0;
        for (col=ia[i]; col<ia[i+1]; col++)
        {
            a[col] = b[j++];
#if 0
            if (a[col] == 0.0)
	    {
		cerr << "APINV0: 0 in col " << i << " of len " << 
		   ia[i+1]-ia[i] << endl;
	    }
#endif
        }
    }

    // transpose the solution in place
    if (transpose)
    {
        csrtrans(n, a, ja, ia);
        delete [] ao;
        delete [] jao;
        delete [] iao;
    }

    delete [] WORK;
    delete [] b;
    delete [] iwk;
    delete [] ind;
}

CSRMat_APINV_BANDED::CSRMat_APINV_BANDED(const CSRMat& A, int semi)
{
    assert(A.nrow == A.ncol);
    int i, nind, colnnz, col, ii, jrow, j;
    double *ahat, *p;

    int n = nrow = ncol = A.nrow;
    nnz = nrow*(2*semi+1) - semi*(semi+1);

    // matrix for which right approx inverse is required
    double *ao;
    int *jao;
    int *iao;

    // transpose matrix if want right apinv
    if (transpose)
    {
        ao = new double[A.nnz];
        jao = new int[A.nnz];
        iao = new int[n+1];
        csrtrans(n, A.a, A.ja, A.ia, ao, jao, iao);
    }
    else
    {
        ao = A.a;
        jao = A.ja;
        iao = A.ia;
    }

    // allocate space for approximate inverse
    a = new double[nnz];
    ja = new int[nnz];
    ia = new int[n+1];

    // structure of ia and ja
    int col1, col2;
    int k = 0;
    ia[0] = 0;
    for (i=0; i<nrow; i++)
    {
        col1 = MAX(i-semi,0);
        col2 = MIN(i+semi,nrow-1);

        for (j=col1; j<=col2; j++)
        {
            ja[k] = j;
            k++;
        }
        ia[i+1] = k;
    }

    // workspaces
    int *ind = new int[n];
    int *iwk = new int[n];
    for (i=0; i<n; i++)
        iwk[i] = -1;

    // b vector and workspace for dgels
    double *b = new double[n];
    integer LWORK = 10000;
    double *WORK = new double[LWORK];

    // main loop over columns
    for (i=0; i<n; i++)
    {
        // determine all nonzeros row indices for 
        // columns indicated in i-th column of A
        nind = 0;
        colnnz = ia[i+1] - ia[i];

        for (j=ia[i]; j<ia[i+1]; j++)
        {
            // scan column
            col = ja[j];
            for (ii=iao[col]; ii<iao[col+1]; ii++)
            {
                jrow = jao[ii];
                if (iwk[jrow] == -1)
                {
                    iwk[jrow] = nind;
                    ind[nind++] = jrow;
                }
            }
        }

        assert (nind > 0);
        assert (colnnz > 0);

        // form least-squares matrix
        ahat = new double[nind*colnnz];
        for (j=0; j<nind*colnnz; j++)
            ahat[j] = 0.0;

        p = ahat;
        for (j=ia[i]; j<ia[i+1]; j++)
        {
            col = ja[j];
            for (ii=iao[col]; ii<iao[col+1]; ii++)
            {
                jrow = jao[ii];
                p[iwk[jrow]] = ao[ii];
            }
            p = p + nind;
        }

        // form b
        for (j=0; j<nind; j++)
            b[j] = 0.0;
        if (iwk[i] >= 0)
            b[iwk[i]] = 1.0;

#if DEBUG
        if (iwk[i] < 0)
           cerr << "APINV_BANDED: no diag in col " << i << endl;
#endif

        // reset workspace
        for (j=0; j<nind; j++)
            iwk[ind[j]] = -1;

        // solve least squares system
        char TRANS = 'N';
        integer NROW = nind;
        integer NCOL = colnnz;
        assert (NROW >= NCOL);
        integer ONE = 1;
        integer INFO;

        F77NAME(dgels)(&TRANS, &NROW, &NCOL, &ONE, ahat, &NROW,
            b, &NROW, WORK, &LWORK, &INFO);
        if (INFO != 0)
            bperror("CSRMat_APINV_BANDED: dgels returned", INFO);

        delete [] ahat;

        // copy result into matrix
        j = 0;
        for (col=ia[i]; col<ia[i+1]; col++)
        {
            a[col] = b[j++];
#if 0
            if (a[col] == 0.0)
	    {
		cerr << "APINV0: 0 in col " << i << " of len " << 
		   ia[i+1]-ia[i] << endl;
	    }
#endif
        }
    }

    // transpose the solution in place
    if (transpose)
    {
        csrtrans(n, a, ja, ia);
        delete [] ao;
        delete [] jao;
        delete [] iao;
    }

    delete [] WORK;
    delete [] b;
    delete [] iwk;
    delete [] ind;
}

CSRMat_DIAG::CSRMat_DIAG(const CSRMat& A)
{
    nrow = A.nrow;
    ncol = A.nrow;
    nnz = A.nrow; // one for each row

    a = new double[nnz];
    ja = new int[nnz];
    ia = new int[nrow+1];

    int i, j;
    int got_diag;

    for (i=0; i<nrow; i++)
    {
        got_diag = FALSE;
        for (j=A.ia[i]; j<A.ia[i+1]; j++)
        {
            if (A.ja[j] == i)
            {
                a[i] = 1.0 / A.a[j];
                ja[i] = i;
                ia[i] = i;
                got_diag = TRUE;
            }
        }
        if (!got_diag)
            bperror("CSRMat_DIAG: missing diagonal entry in local matrix", 0);
    }
    ia[nrow] = nrow;
}

// discard everything outside tridiagonal structure
// matrix is stored in sparse diagonal format

CSRMat_TRIDIAG::CSRMat_TRIDIAG(const CSRMat& A)
{
    assert (A.nrow == A.ncol);
    int i, j;

    nrow = ncol = A.nrow;
    nnz = 3*nrow;

    val = new double[nnz];
    for (i=0; i<nnz; i++)
        val[i] = 0.0;

    for (i=0; i<nrow; i++)
    {
        for (j=A.ia[i]; j<A.ia[i+1]; j++)
        {
            if (A.ja[j] == i-1)
                val[i] = A.a[j];
            else if (A.ja[j] == i)
                val[nrow+i] = A.a[j];
            else if (A.ja[j] == i+1)
                val[2*nrow+i] = A.a[j];
        }
    }

    // Factor without pivoting.
    // First some aliases, for my own sanity.
    double *e = &val[1]; //val[0] is empty
    double *d = &val[nrow];
    double *f = &val[2*nrow];

    for (i=0; i<nrow-1; i++)
    {
        e[i] /= d[i];
        d[i+1] -= e[i]*f[i];
    }
}

void CSRMat_TRIDIAG::Mat_Vec_Solve(const BlockVec& B, BlockVec& X) const
{
    // Should call ddiasm in sparse blas.

    // First some aliases, for my own sanity.
    int i, n = nrow;
    double *e = &val[0]; // different indexing
    double *d = &val[nrow];
    double *f = &val[2*nrow];

    // Solve in place.
    X.BlockCopy(B);

    for (int j=0; j<X.dim1; j++)
    {
        double *x = X.v + j*X.dim0;

        // Lower triangular solve.
	for (i=1; i<n; i++)
	    x[i] -= e[i]*x[i-1];

	// Upper triangular solve.
	x[n-1] /= d[n-1];
	for (i=n-2; i>=0; i--)
	{
	    x[i] -= f[i]*x[i+1];
	    x[i] /= d[i];
	}
    }
}

// assumes matrix in sorted order

CSRMat_SOR::CSRMat_SOR(const CSRMat& A, 
    const double omega, const int iterations)
{
    int i, j, nrow;
    int got_diag;

    // Store parameters for use in iterations.

    Ap = &A;
    omega_ = omega;
    iterations_ = iterations;

    nrow = A.numrow();
    idiag = new int[nrow];

    // search for diagonal elements
    for (i=0; i<nrow; i++)
    {
        got_diag = FALSE;
        for (j=A.row_ptr(i); j<A.row_ptr(i+1); j++)
        {
            if (A.col_ind(j) == i)
            {
                idiag[i] = j;
                got_diag = TRUE;
            }
        }
        if (!got_diag)
            bperror("CSRMat_SOR: missing diagonal entry in local matrix", 0);
    }
}

// this does not handle blocks of vectors

void CSRMat_SOR::Mat_Vec_Solve(const BlockVec& U, BlockVec& V) const
{
    const int *ia = &Ap->row_ptr(0);
    const int *ja = &Ap->col_ind(0);
    const double *a = &Ap->val(0);
    int it, i, j;
    double temp;

    V.BlockCopy(U);
    double *v = V.v;
    const double *u = U.v;

    // Specialized code for first step.

    for (i=0; i<Ap->numrow(); i++)
    {
        for (j=ia[i]; j<idiag[i]; j++)
        {
            v[i] = v[i] - omega_ * a[j] * v[ja[j]];
        }
        v[i] /= a[idiag[i]];
    }

    for (i=0; i<Ap->numrow(); i++)
        v[i] *= omega_;

    // After first step....

    for (it=1; it<iterations_; it++)
    {
        for (i=0; i<Ap->numrow(); i++)
        {
            temp = u[i];

            for (j=ia[i]; j<idiag[i]; j++)
            {
                temp = temp - a[j] * v[ja[j]];
            }
            for (j=idiag[i]+1; j<ia[i+1]; j++)
            {
                temp = temp - a[j] * v[ja[j]];
            }

            temp /= a[idiag[i]];

            v[i] = v[i] + omega_ * (temp - v[i]);
        }
    }
}

CSRMat_SSOR::CSRMat_SSOR(const CSRMat& A, 
    const double omega, const int iterations)
{
    int i, j, nrow;
    int got_diag;

    // Store parameters for use in iterations.

    Ap = &A;
    omega_ = omega;
    iterations_ = iterations;

    nrow = A.numrow();
    idiag = new int[nrow];

    // search for diagonal elements
    for (i=0; i<nrow; i++)
    {
        got_diag = FALSE;
        for (j=A.row_ptr(i); j<A.row_ptr(i+1); j++)
        {
            if (A.col_ind(j) == i)
            {
                idiag[i] = j;
                got_diag = TRUE;
            }
        }
        if (!got_diag)
            bperror("CSRMat_SSOR: missing diagonal entry in local matrix", 0);
    }
}

void CSRMat_SSOR::Mat_Vec_Solve(const BlockVec& U, BlockVec& V) const
{
    const int *ia = &Ap->row_ptr(0);
    const int *ja = &Ap->col_ind(0);
    const double *a = &Ap->val(0);
    int it, i, j;
    double temp;

#if 0
    V.BlockCopy(U);
    double *v = V.v;
    const double *u = U.v;

    // Specialized code for first step.

    // lower sweep
    for (i=0; i<Ap->numrow(); i++)
    {
        for (j=ia[i]; j<idiag[i]; j++)
        {
            v[i] = v[i] - omega_ * a[j] * v[ja[j]];
        }
        v[i] /= a[idiag[i]];
    }

    // multiply by diagonal blocks
    for (i=0; i<Ap->numrow(); i++)
        v[i] *= a[idiag[i]];

    // upper sweep
    for (i=Ap->numrow()-1; i>=0; i--)
    {
        for (j=idiag[i]+1; j<ia[i+1]; j++)
        {
            v[i] = v[i] - omega_ * a[j] * v[ja[j]];
        }
        v[i] /= a[idiag[i]];
    }

    // After first step....
#endif

    V.BlockSetToZero();
    double *v = V.v;
    const double *u = U.v;

    for (it=1; it<=iterations_; it++)
    {
        for (i=0; i<Ap->numrow(); i++)
        {
            temp = u[i];

            for (j=ia[i]; j<idiag[i]; j++)
            {
                temp = temp - a[j] * v[ja[j]];
            }
            for (j=idiag[i]+1; j<ia[i+1]; j++)
            {
                temp = temp - a[j] * v[ja[j]];
            }

            temp /= a[idiag[i]];

            v[i] = v[i] + omega_ * (temp - v[i]);
        }
        for (i=Ap->numrow()-1; i>=0; i--)
        {
            temp = u[i];

            for (j=ia[i]; j<idiag[i]; j++)
            {
                temp = temp - a[j] * v[ja[j]];
            }
            for (j=idiag[i]+1; j<ia[i+1]; j++)
            {
                temp = temp - a[j] * v[ja[j]];
            }

            temp /= a[idiag[i]];

            v[i] = v[i] + omega_ * (temp - v[i]);
        }
    }
}

void CSRMat::Print(ostream& os) const
{
        assert (a != NULL);

        int M = numrow();
        int N = numcol();
        int rowp1, colp1;
        int flag = 0;

#ifdef _GNU_GCC_3
        ios::fmtflags olda = os.setf(ios::right,ios::adjustfield);
        ios::fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
#else
#ifdef _GNU_GCC_296
        long olda = os.setf(ios::right,ios::adjustfield);
        long oldf = os.setf(ios::scientific,ios::floatfield);
#else
        long olda = os.setf(ios::right,ios::adjustfield);
        long oldf = os.setf(ios::scientific,ios::floatfield);
#endif
#endif

        int oldp = os.precision(4);

        for (int i = 0; i < M ; i++)
           for (int j=row_ptr(i);j<row_ptr(i+1);j++)
           {
              rowp1 =  i + 1;
              colp1 =  col_ind(j) + 1;
              if ( rowp1 == M && colp1 == N ) flag = 1;
              os.width(8);
              os <<  rowp1;
              os.width(8);
              os <<  colp1;
              os.width(14);
              os <<  val(j) << endl;
           }

        if (flag == 0)
        {
           os.width(8);
           os <<  M ;
           os.width(8);
           os <<  N ;
           os.width(14);
           os << 0.0 << endl;
        }

        os.setf(olda,ios::adjustfield);
        os.setf(oldf,ios::floatfield);
        os.precision(oldp);
}

// non-member functions

ostream& operator << (ostream& os, const CSRMat& mat)
{
        assert (mat.a != NULL);

        int M = mat.numrow();
        int N = mat.numcol();
        int rowp1, colp1;
        int flag = 0;

#ifdef _GNU_GCC_3
        ios::fmtflags olda = os.setf(ios::right,ios::adjustfield);
        ios::fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
#else
#ifdef _GNU_GCC_296
        long olda = os.setf(ios::right,ios::adjustfield);
        long oldf = os.setf(ios::scientific,ios::floatfield);
#else
        long olda = os.setf(ios::right,ios::adjustfield);
        long oldf = os.setf(ios::scientific,ios::floatfield);
#endif
#endif

        int oldp = os.precision(4);

//      Loop through rows...
        for (int i = 0; i < M ; i++)
           for (int j=mat.row_ptr(i);j<mat.row_ptr(i+1);j++)
           {
              rowp1 =  i + 1;
              colp1 =  mat.col_ind(j) + 1;
              if ( rowp1 == M && colp1 == N ) flag = 1;
              os.width(8);
              os <<  rowp1;
              os.width(8);
              os <<  colp1;
              os.width(14);
              os <<  mat.val(j) << "\n";
           }

        if (flag == 0)
        {
           os.width(8);
           os <<  M ; os << "    " ;
           os.width(8);
           os <<  N ; os << "    " ;
           os.width(14);
           os <<  0.0 << "\n";
        }

        os.setf(olda,ios::adjustfield);
        os.setf(oldf,ios::floatfield);
        os.precision(oldp);

        return os;
}
