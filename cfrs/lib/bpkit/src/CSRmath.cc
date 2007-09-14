//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#include "CSRMat.h"

// B = A', B is released

void CSRMat::Mat_Trans(LocalMat *B) const
{
    CSRMat& b = * (CSRMat *) B;
    int i, j, k, next;

    delete [] b.a;
    delete [] b.ja;
    delete [] b.ia;
    b.set(new double[nnz], new int[nnz], new int[ncol+1], ncol, nrow, nnz);

    // some aliases for speed and readability
    double *ao  = b.a;
    int    *jao = b.ja;
    int    *iao = b.ia;

    // initialize
    for (i=0; i<=ncol; i++)
        iao[i] = 0;

    // compute lengths of rows of transpose of A
    for (i=0; i<nrow; i++)
        for (k=ia[i]; k<ia[i+1]; k++)
            iao[ja[k]+1]++;

    // compute pointers from lengths
    iao[0] = 0;
    for (i=0; i<ncol; i++)
        iao[i+1] = iao[i] + iao[i+1];

    // now do the actual copying
    for (i=0; i<nrow; i++)
    {
        for (k=ia[i]; k<ia[i+1]; k++)
	{
            j = ja[k];
            next = iao[j];
            ao[next] = a[k];
            jao[next] = i;
            iao[j] = next+1;
        }
    }

    // reshift iao
    for (i=ncol-1; i>=0; i--)
        iao[i+1] = iao[i];
    iao[0] = 0;
}

// C = A + alpha B, C is released

void CSRMat::Mat_Mat_Add(const LocalMat *B_, LocalMat *C_,
    double alpha) const
{
    CSRMat& B = *(CSRMat *) B_;
    CSRMat& C = *(CSRMat *) C_;

    assert (nrow == B.nrow);
    assert (ncol == B.ncol);

    // semantics of nnz: it contains the ACTUAL number of nonzeros
    const int nzmax = nnz + B.nnz;

    // Release any existing data of C.
    delete [] C.a;
    delete [] C.ja;
    delete [] C.ia;

    C.nrow = nrow;
    C.ncol = ncol;
    C.nnz = -1; // unknown at the moment

    C.a = new double[nzmax];
    C.ja = new int[nzmax];
    C.ia = new int[nrow+1];

    int *iw = new int[ncol];
    for (int i=0; i<ncol; i++)
        iw[i] = -1;

    // some aliases for speed and readability
    const double *b  = B.a;
    const int    *jb = B.ja;
    const int    *ib = B.ia;
          double *c  = C.a;
          int    *jc = C.ja;
          int    *ic = C.ia;

    int ii, ka, kb, kc, jcol, jpos;

    kc = 0;
    ic[0] = 0;

    for (ii=0; ii<C.nrow; ii++)
    {
	for (ka=ia[ii]; ka<ia[ii+1]; ka++)
	{
            jcol = ja[ka];
            jc[kc] = jcol;
             c[kc] = a[ka];
            iw[jcol] = kc;
            kc++;
	}

        for (kb=ib[ii]; kb<ib[ii+1]; kb++)
        {
            jcol = jb[kb];
            jpos = iw[jcol];

            if (jpos < 0)
            {
                jc[kc] = jcol;
                 c[kc] = alpha*b[kb];
                iw[jcol] = kc;
                kc++;
            }
            else
            {
                c[jpos] += alpha*b[kb];
            }
        }
        ic[ii+1] = kc;
        for (int i=ic[ii]; i<kc; i++)
            iw[jc[i]] = -1;
    }
    delete [] iw;
    C.nnz = kc;
}

// D = alpha A B + beta D, D is released
// The variable C is used as a temporary inside the code.

void CSRMat::Mat_Mat_Mult(const LocalMat *B_, LocalMat *D_,
    double alpha, double beta) const
{
    CSRMat& B = *(CSRMat *) B_;
    CSRMat& D = *(CSRMat *) D_;

    assert (B.ia != NULL); // make sure they are not stored in factorized form

    assert (ncol == B.nrow);
    if (beta != 0.0)
    {
        assert (nrow == D.nrow);
        assert (B.ncol == D.ncol);
        assert (D.ia != NULL); // but a may be NULL
    }

    // some aliases for speed and readability
    const double *b  = B.a;
    const int    *jb = B.ja;
    const int    *ib = B.ia;
    const double *d  = D.a;
    const int    *jd = D.ja;
    const int    *id = D.ia;
          int    *jc = new int[B.ncol];

    int nzmax, i, ii, jj, ka, kb, kc, kd, jcol, jpos;
    double scal;

    int *iw = new int[B.ncol];
    for (i=0; i<B.ncol; i++)
	iw[i] = -1;

    //
    // first part: count number of nnz in product
    //

    nzmax = 0;

    for (ii=0; ii<nrow; ii++)                 // row ii
    {
        kc = 0;

	if (beta != 0.0)  // first add in row of D if necessary
	{
	    for (kd=id[ii]; kd<id[ii+1]; kd++)  // row ii of D
	    {
		jcol = jd[kd];
		jc[kc] = jcol;
		iw[jcol] = kc;
		kc++;
	    }
	}

	for (ka=ia[ii]; ka<ia[ii+1]; ka++)      // row ii of A
	{
	    jj = ja[ka];

	    for (kb=ib[jj]; kb<ib[jj+1]; kb++)  // row jj of B
	    {
		jcol = jb[kb];
		jpos = iw[jcol];

		if (jpos < 0)
		{
		    jc[kc] = jcol;
		    iw[jcol] = kc;
		    kc++;
		}
	    }
	}
	for (i=0; i<kc; i++)
	    iw[jc[i]] = -1;

	nzmax += kc;
    }

    delete [] jc;

    //
    // second part: form actual product with values
    //

    CSRMat C;
    C.set(new double[nzmax], new int[nzmax], new int[nrow+1],
	nrow, B.ncol, -1);

    // more aliases
    double *c  = C.a;
            jc = C.ja;
    int    *ic = C.ia;

    kc = 0;
    ic[0] = 0;

    for (ii=0; ii<C.nrow; ii++)                 // row ii
    {
	if (beta != 0.0)  // first add in row of D if necessary
	{
	    for (kd=id[ii]; kd<id[ii+1]; kd++)  // row ii of D
	    {
		jcol = jd[kd];
		jc[kc] = jcol;
		 c[kc] = beta*d[kd];
		iw[jcol] = kc;
		kc++;
	    }
	}

	for (ka=ia[ii]; ka<ia[ii+1]; ka++)      // row ii of A
	{
	    scal = alpha * a[ka];
	    jj = ja[ka];

	    for (kb=ib[jj]; kb<ib[jj+1]; kb++)  // row jj of B
	    {
		jcol = jb[kb];
		jpos = iw[jcol];

		if (jpos < 0)
		{
		    jc[kc] = jcol;
		     c[kc] = scal*b[kb];
		    iw[jcol] = kc;
		    kc++;
		}
		else
		{
		    c[jpos] += scal*b[kb];
		}
	    }
	}
	ic[ii+1] = kc;
	for (i=ic[ii]; i<kc; i++)
	    iw[jc[i]] = -1;
    }

    delete [] iw;

    // move result to matrix D (by copying pointers of temp C)

    D.nrow = nrow;
    D.ncol = B.ncol;
    D.nnz = kc;
    delete [] D.a;
    delete [] D.ja;
    delete [] D.ia;
    D.a  = C.a;
    D.ja = C.ja;
    D.ia = C.ia;
    C.a  = NULL; // need this, or destructor will free memory
    C.ja = NULL;
    C.ia = NULL;
}
