//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#ifndef _DENSEMAT_H_
#define _DENSEMAT_H_

#include <cmath>
#include <cassert>

using namespace std;

#include "LocalMat.h"
#include "BlockVec.h"
#include "LocalPrecon.h"
#include "BPKIT.h"
#include "blas3.h"
#include "lapackd.h"

// DenseMat objects typically allocate their own storage for matrix data.
// The default constructor is the only constructor that does not allocate
// storage.  Thus this constructor may be used to control memory
// allocation, for example by allocating an array of DenseMat using new.
// In this case, the data pointer must be set to NULL before the destructor
// is called.

class DenseMat : public LocalMat
{
protected:
    double *a;
    int    nrow;
    int    ncol;

public:
    DenseMat() {a = NULL; nrow = ncol = 0;}
   ~DenseMat() {delete [] a;}

    DenseMat(const int r, const int c)
        {nrow = r; ncol = c; a = new double[r*c];}

    DenseMat(const DenseMat& A);
// cannot inline on Cray
//      {nrow = A.nrow; 
//       ncol = A.ncol;
//       register double *p = a = new double[nrow*ncol];
//       register double *q = A.a;
//       for (int i=0; i<nrow*ncol; i++) *p++ = *q++;}

    void set(const int r, const int c, double *d)
        {nrow = r; ncol = c; a = d;}

    int numrow() const {return nrow;}
    int numcol() const {return ncol;}

    double operator=(double s)
        {register double *p = a;
         for (int i=0; i<nrow*ncol; i++) *p++ = s; return s;}

    // 0-based indexing
    const double& operator()(const unsigned int i, const unsigned int j) const
        {return a[j*nrow+i];}
    double& operator()(const unsigned int i, const unsigned int j)
        {return a[j*nrow+i];}

    // virtual functions

    double *& Data() {return a;}
    const double *Data() const {return a;}
    LocalMat *CreateEmpty() const {return new DenseMat();}
    inline LocalMat *CreateInv(LocalPrecon&) const;
    inline void SetToZero(int, int);
    void MatCopy(const LocalMat& A)  // UNDONE: realloc if necessary
        {register double *p = a;
         register double *q = ((DenseMat *)&A)->a;
         for (int i=0; i<nrow*ncol; i++) *p++ = *q++;}
    void Print(ostream&) const;

    // virtual mathematical functions

    inline void Mat_Trans(LocalMat *B) const;
    inline void Mat_Mat_Add(const LocalMat *B, LocalMat *C, double alpha) const;
    inline void Mat_Mat_Mult(const LocalMat *B, LocalMat *C, 
        double alpha, double beta) const;
    inline void Mat_Vec_Mult(const BlockVec& B, BlockVec& C, 
        double alpha, double beta) const;
    inline void Mat_Trans_Vec_Mult(const BlockVec& B, BlockVec& C,
        double alpha, double beta) const;
    inline void Mat_Vec_Solve(const BlockVec& B, BlockVec& C) const;
    inline void Mat_Trans_Vec_Solve(const BlockVec& B, BlockVec& C) const;
};

ostream& operator << (ostream& os, const DenseMat& mat);

class DenseMat_LU : public DenseMat
{
private:
    double  *lu;
    integer *ipiv;  // FORTRAN integer

public:
    inline DenseMat_LU(const DenseMat&);
    inline void Mat_Vec_Solve(const BlockVec&, BlockVec&) const;
    inline ~DenseMat_LU();
};

class DenseMat_INVERSE : public DenseMat
{
private:
public:
    inline DenseMat_INVERSE(const DenseMat&);
    void Mat_Vec_Solve(const BlockVec& B, BlockVec& X) const
        {solve_is_mult(B, X);}
};

class DenseMat_SVD : public DenseMat
{
private:
public:
    inline DenseMat_SVD(const DenseMat&, double, double);
    void Mat_Vec_Solve(const BlockVec& B, BlockVec& X) const
        {solve_is_mult(B, X);}
};

class DenseMat_DIAGDOM : public DenseMat
{
private:
public:
    inline DenseMat_DIAGDOM(const DenseMat&, double alpha);
    void Mat_Vec_Solve(const BlockVec& B, BlockVec& X) const
        {solve_is_mult(B, X);}
};

class DenseMat_GERSH : public DenseMat
{
private:
public:
    inline DenseMat_GERSH(const DenseMat&, double);
    void Mat_Vec_Solve(const BlockVec& B, BlockVec& X) const
        {solve_is_mult(B, X);}
};

// B = A'
// must copy data
inline void DenseMat::Mat_Trans(LocalMat *B) const
{
    DenseMat& b = *(DenseMat *) B;

    assert (this != &b);       // in-place not allowed

    if (b.a == NULL && b.nrow == 0 && b.ncol == 0)
    {
	b.nrow = ncol;
	b.ncol = nrow;
	b.a = new double[b.nrow*b.ncol];
    }

    assert (a != NULL);
    assert (b.a != NULL);
    assert (b.nrow == ncol);
    assert (b.ncol == nrow);

    int i, j;
    double *p, *q;
    
    // Traverse B and fill from A.

    p = b.a;
    for (i=0; i<nrow; i++)
    {
	q = a+i;
        for (j=0; j<ncol; j++)
	{
	    *p++ = *q;
	    q += nrow;
	}
    }
}

// C = A + alpha B
inline void DenseMat::Mat_Mat_Add(const LocalMat *B, LocalMat *C,
    double alpha) const
{
    DenseMat& b = *(DenseMat *) B;
    DenseMat& c = *(DenseMat *) C;

    if (c.a == NULL && c.nrow == 0 && c.ncol == 0)
    {
	c.nrow = nrow;
	c.ncol = ncol;
	c.a = new double[c.nrow*c.ncol];
    }

    assert (a != NULL);
    assert (b.a != NULL);
    assert (c.a != NULL);
    assert (nrow == b.nrow);
    assert (ncol == b.ncol);
    assert (nrow == c.nrow);
    assert (ncol == c.ncol);

    int i;
    double *ap = a;
    double *bp = b.a;
    double *cp = c.a;

    if (alpha == 1.0)
        for (i=0; i<nrow*ncol; i++)
            *cp++ = *ap++ + *bp++;

    else if (alpha == -1.0)
        for (i=0; i<nrow*ncol; i++)
	    *cp++ = *ap++ - *bp++;

    else
        for (i=0; i<nrow*ncol; i++)
	    *cp++ = *ap++ + alpha * *bp++;
}

// C = alpha A B + beta C
inline void DenseMat::Mat_Mat_Mult(const LocalMat *B, LocalMat *C, 
    double alpha, double beta) const
{
    assert (B != C);
    DenseMat& b = *(DenseMat *) B;
    DenseMat& c = *(DenseMat *) C;

    if (beta == 0.0 && c.a == NULL && c.nrow == 0 && c.ncol == 0)
    {
	c.nrow = nrow;
	c.ncol = b.ncol;
	c.a = new double[c.nrow*c.ncol];
    }

// added Feb 26, 1996
    if (beta == 0.0 && (c.nrow != nrow || c.ncol != b.ncol))
    {
        assert (c.a != NULL);
	delete c.a;

	c.nrow = nrow;
	c.ncol = b.ncol;
	c.a = new double[c.nrow*c.ncol];
    }

    assert (a != NULL);
    assert (b.a != NULL);
    assert (c.a != NULL);

    char transa = 'N';
    char transb = 'N';
    integer M = nrow;
    integer N = c.ncol;
    integer K = ncol;
    integer LDA = nrow;
    integer LDB = b.nrow;
    integer LDC = c.nrow;

    F77NAME(dgemm)(&transa, &transb, &M, &N, &K, &alpha, a, &LDA, 
        b.a, &LDB, &beta, c.a, &LDC);
}

// C = alpha A B + beta C
inline void DenseMat::Mat_Vec_Mult(const BlockVec& B, BlockVec& C,
    double alpha, double beta) const
{
    char trans = 'N';
    integer M = nrow;
    integer N = C.dim1;
    integer K = ncol;
    integer LDA = nrow;
    integer LDB = B.dim0;
    integer LDC = C.dim0;

    assert (a != NULL);
    assert (&B != &C);

    F77NAME(dgemm)(&trans, &trans, &M, &N, &K, &alpha, a, &LDA, 
        B.v, &LDB, &beta, C.v, &LDC);
}

inline void DenseMat::Mat_Trans_Vec_Mult(const BlockVec& B, 
    BlockVec& C, double alpha, double beta) const
{
    char transa = 'T';
    char transb = 'N';
    integer M = ncol;
    integer N = C.dim1;
    integer K = nrow;
    integer LDA = nrow;
    integer LDB = B.dim0;
    integer LDC = C.dim0;

    assert (a != NULL);
    assert (&B != &C);

    F77NAME(dgemm)(&transa, &transb, &M, &N, &K, &alpha, a, &LDA,
        B.v, &LDB, &beta, C.v, &LDC);
}

inline void DenseMat::Mat_Vec_Solve(const BlockVec&, BlockVec&) const
{
    bperror("DenseMat::Mat_Vec_Solve: to solve with a local matrix,\n"
            "a local preconditioner should be used.\n", 0);
}

inline void DenseMat::Mat_Trans_Vec_Solve(const BlockVec&, BlockVec&) const
{
    bperror("DenseMat::Mat_Trans_Vec_Solve: to solve with a local matrix,\n"
            "a local preconditioner should be used.\n", 0);
}

inline void DenseMat::SetToZero(int r, int c)
{
    if (a == NULL && nrow == 0 && ncol == 0)
    {
	nrow = r;
	ncol = c;
	a = new double[nrow*ncol];
    }

    assert (a != NULL);
    assert (nrow == r);
    assert (ncol == c);

    double *p = a;
    for (int i=0; i<nrow*ncol; i++)
	*p++ = 0.0;
}

inline LocalMat *DenseMat::CreateInv(LocalPrecon& local_precon) const
{
    assert (a != NULL);

    switch (local_precon.name)
    {
    case LP_LU:
        return new DenseMat_LU(*this);
    case LP_INVERSE:
        return new DenseMat_INVERSE(*this);
    case LP_SVD:
        return new DenseMat_SVD(*this, local_precon.darg1, local_precon.darg2);
    case LP_DIAGDOM:
        return new DenseMat_DIAGDOM(*this, local_precon.darg1);
    case LP_GERSH:
        return new DenseMat_GERSH(*this, local_precon.darg1);
    default:
        bperror("The local preconditioner you have chosen is not available\n"
                "for the type of blocks (DenseMat) in the block matrix.", 0);
    }
    return NULL;
}

inline DenseMat_LU::DenseMat_LU(const DenseMat& A)
{
    a = NULL; // indicate this is implied inverse
    nrow = A.numrow();
    ncol = A.numcol();

    assert (nrow == ncol);

    lu = new double[nrow*ncol];
    // copy to lu

    double *p = lu;
    const double *q = &A(0,0);
    for (int i=0; i<nrow*ncol; i++)
        *p++ = *q++;

    integer M = nrow;
    integer N = ncol;
    integer LDA = nrow;
    integer INFO;

    ipiv = new integer[nrow];

    F77NAME(dgetrf)(&M, &N, lu, &LDA, ipiv, &INFO);
    if (INFO != 0)
	bperror("DenseMat_LU: dgetrf error", INFO);
}

inline void DenseMat_LU::Mat_Vec_Solve(const BlockVec& b, BlockVec& x) const
{
    char trans = 'N';
    integer N = ncol;
    integer NRHS = b.dim1;
    integer LDA = nrow;
    integer LDB = b.dim0;
    integer INFO;

    if (b.v != x.v)
        x.BlockCopy(b);

    F77NAME(dgetrs)(&trans, &N, &NRHS, lu, &LDA, ipiv, 
        x.v, &LDB, &INFO);
    if (INFO != 0)
	bperror("DenseMat_LU: dgetrs error", INFO);
}

inline DenseMat_LU::~DenseMat_LU()
{
    delete [] ipiv;
    delete [] lu;
}

inline DenseMat_INVERSE::DenseMat_INVERSE(const DenseMat& A)
    : DenseMat(A)
{
    integer M = nrow;
    integer N = ncol;
    integer LDA = nrow;
    integer INFO;

    assert (nrow == ncol);

    integer *ipiv = new integer[nrow];
    integer LWORK = N*N;
    double *work = new double[LWORK];

    // LU factorize
    F77NAME(dgetrf)(&M, &N, a, &LDA, ipiv, &INFO);
    if (INFO != 0)
	bperror("DenseMat_INVERSE: dgetrf error", INFO);

    // compute inverse
    F77NAME(dgetri)(&N, a, &LDA, ipiv, work, &LWORK, &INFO);
    if (INFO != 0)
	bperror("DenseMat_INVERSE: dgetri error", INFO);

    delete [] ipiv;
    delete [] work;
}

inline DenseMat_SVD::DenseMat_SVD(const DenseMat& A, double rthresh, 
  double athresh)
    : DenseMat(A)
{

  cout<< "\n NO DenseMat_SVD due to no dgesvd from LAPACK "; exit(1);

  // ESSL does have dgesvf though 

//     assert (nrow == ncol);

//     int i, j;
//     double thresh;
//     double *u  = new double[nrow*ncol];
//     double *s  = new double[nrow];
//     double *vt = new double[nrow*ncol];

//     char job = 'A';
//     integer N = nrow;
//     integer LWORK = 5*N;
//     double *work = new double[LWORK];
//     integer INFO;

//     F77NAME(dgesvd) (&job, &job, &N, &N, a, &N, s, u, &N, vt, &N,
// 	work, &LWORK, &INFO);
//     if (INFO != 0)
// 	bperror("DenseMat_SVD: dgesvd error", INFO);

//     delete [] work;

//     // apply threshold
//     thresh = s[0]*rthresh + athresh;
//     for (i=0; i<nrow; i++)
// 	if (s[i] < thresh)
// 	    s[i] = thresh;

//     if (s[0] == 0.0)
// 	bperror("DenseMat_SVD: block is zero after thresholding", 0);

//     // scale the columns of u with reciprocal singular values
//     double *p = u;
//     for (i=0; i<ncol; i++)
//     {
// 	double scal = s[i];
// 	for (j=0; j<nrow; j++)
// 	    *p++ /= scal;
//     }

//     // multiply through and store the result
//     char trans = 'T';
//     double alpha = 1.0;
//     double beta = 0.0;
//     F77NAME(dgemm)(&trans, &trans, &N, &N, &N, &alpha, vt, &N,
//         u, &N, &beta, a, &N);

//     delete [] u;
//     delete [] s;
//     delete [] vt;
}

inline DenseMat_DIAGDOM::DenseMat_DIAGDOM(const DenseMat& A, double alpha)
    : DenseMat(A)
{
    int i, j;
    integer M = nrow;
    integer N = ncol;
    integer LDA = nrow;
    integer INFO;

    assert (nrow == ncol);

    integer *ipiv = new integer[nrow];
    integer LWORK = N*N;
    double *work = new double[LWORK];

    // compute sum of abs of off-diagonals and put in work array
    double *p = a;
    for (i=0; i<nrow; i++)
        work[i] = 0.0;

    for (j=0; j<ncol; j++)
    {
        for (i=0; i<nrow; i++)
        {
            if (i != j)
            {
                work[i] += ABS(*p);
            }
            p++;
        }
    }

    for (i=0; i<nrow; i++)
    {
	if (ABS(a[i*nrow+i]) < alpha*work[i])
	    a[i*nrow+i] = SGN(a[i*nrow+i])*alpha*work[i];
    }

    // LU factorize
    F77NAME(dgetrf)(&M, &N, a, &LDA, ipiv, &INFO);
    if (INFO != 0)
	bperror("DenseMat_DIAGDOM: dgetrf error", INFO);

    // compute inverse
    F77NAME(dgetri)(&N, a, &LDA, ipiv, work, &LWORK, &INFO);
    if (INFO != 0)
	bperror("DenseMat_DIAGDOM: dgetri error", INFO);

    delete [] ipiv;
    delete [] work;
}

inline DenseMat_GERSH::DenseMat_GERSH(const DenseMat& A, double alpha)
    : DenseMat(A)
{
    int i, j;
    integer M = nrow;
    integer N = ncol;
    integer LDA = nrow;
    integer INFO;

    assert (nrow == ncol);

    integer *ipiv = new integer[nrow];
    integer LWORK = N*N;
    double *work = new double[LWORK];

    // compute sum of abs of off-diagonals and put in work array
    double *p = a;
    for (i=0; i<nrow; i++)
        work[i] = 0.0;

    for (j=0; j<ncol; j++)
    {
        for (i=0; i<nrow; i++)
        {
            if (i != j)
            {
                work[i] += ABS(*p);
            }
            p++;
        }
    }

    double aii;
    for (i=0; i<nrow; i++)
    {
	aii = a[i*nrow+i];

	if (aii >= 0.0)
        {
	    if (aii - work[i] < alpha)
		a[i*nrow+i] += alpha - aii + work[i];
	}
	else
	{
	    if (aii + work[i] > -alpha)
		a[i*nrow+i] -= alpha + aii + work[i];
	}
    }

    // LU factorize
    F77NAME(dgetrf)(&M, &N, a, &LDA, ipiv, &INFO);
    if (INFO != 0)
	bperror("DenseMat_GERSH: dgetrf error", INFO);

    // compute inverse
    F77NAME(dgetri)(&N, a, &LDA, ipiv, work, &LWORK, &INFO);
    if (INFO != 0)
	bperror("DenseMat_GERSH: dgetri error", INFO);

    delete [] ipiv;
    delete [] work;
}

#endif // _DENSEMAT_H_
