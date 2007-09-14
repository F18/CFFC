//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#include <cassert>

using namespace std;

#include "BPKIT.h"
#include "CSRMat.h"

double CSRMat_ILUT::permtol = BpGetDouble(ILUT_PERMTOL, 0.0);
int    CSRMat_ILUT::pblock  = BpGetInt(ILUT_PBLOCK, 0);
int    CSRMat_ILUT::stat    = BpGetInt(ILUT_STAT, 0);

CSRMat_ILUT::CSRMat_ILUT(const CSRMat& A,
    int ilut_lfil, double ilut_thresh)
{
    numrow = A.numrow();
    int ierr;

    // allocate space for factorization
    int iwk = 2*numrow*ilut_lfil+numrow+1;
    alu = new double[iwk];
    jlu = new int[iwk];
    ju  = new int[numrow];

    // allocate workspace
    double *w = new double[2*numrow];
    int *jw = new int[2*numrow];

    // The fortran has been modified to handle 0-based indexing.

    if (permtol == 0.0)
    {
        F77NAME(bpilut)(&numrow, A.a, A.ja, A.ia, &ilut_lfil, &ilut_thresh,
	    alu, jlu, ju, &iwk, w, jw, &ierr);
    }
    else
    {
        int *iperm = new int[2*numrow];
	int mbloc = pblock;

	if (mbloc == 0)
	    mbloc = numrow;

        F77NAME(bpilutp)(&numrow, A.a, A.ja, A.ia, &ilut_lfil, &ilut_thresh,
	    &permtol, &mbloc,
	    alu, jlu, ju, &iwk, w, jw, iperm, &ierr);

        delete [] iperm;
    }

    if (stat)
    {
        F77NAME(bpilutstat)(&numrow, A.a, A.ja, A.ia, alu, jlu, ju, jw);
    }

    // delete workspace
    delete [] jw;
    delete [] w;

    if (ierr != 0)
	bperror("CSRMat_ILUT: ILUT error", ierr);
}

void CSRMat_ILUT::Mat_Vec_Solve(const BlockVec& B, BlockVec& X)
    const
{
    assert (X.dim1 == B.dim1);

    double *b, *x;
    for (int j=0; j<X.dim1; j++)
    {
        b = B.v + j*B.dim0;
        x = X.v + j*X.dim0;

        F77NAME(bplusol)(&numrow, b, x, alu, jlu, ju);
    }
}

CSRMat_ILUT::~CSRMat_ILUT()
{
    delete [] alu;
    delete [] jlu;
    delete [] ju;
}
