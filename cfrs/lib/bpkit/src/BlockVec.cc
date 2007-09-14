//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#include "BlockVec.h"

BlockVec::BlockVec(const BlockVec& A)
{
    dim0 = A.dim0;
    dim1 = A.dim1;
    ld = A.dim0;
    size0 = A.dim0;
    base = new double[dim0*dim1];
    v = base;
    partit = A.partit;
    owndata = 1;
    VecCopy(A);
}

BlockVec::BlockVec(const BlockVec& A, int index)
{
if (index >= 0)
{
    dim0 = A.partit[index+1] - A.partit[index];;
    dim1 = A.dim1;
    ld = dim0;
    size0 = dim0;
    base = new double[dim0*dim1];
    v = base;
    partit = NULL;
    owndata = 1;

    // copy the values

    int i, j;
    double *p;
    const double *q;

    for (i=0; i<dim1; i++)
    {
	p = v + i*ld;
        q = A.base + A.partit[index] + i*A.ld;
        for (j=0; j<size0; j++)
            *p++ = *q++;
    }
}
else
{
    dim0 = A.size0;
    dim1 = A.dim1;
    ld = dim0;
    size0 = dim0;
    base = new double[dim0*dim1];
    v = base;
    partit = NULL;
    owndata = 1;

    // copy the values

    int i, j;
    double *p;
    const double *q;

    for (i=0; i<dim1; i++)
    {
	p = v + i*ld;
        q = A.v + i*A.ld;
        for (j=0; j<size0; j++)
            *p++ = *q++;
    }
}
}

void BlockVec::VecCopy(const BlockVec& A)
{
#ifdef DEBUG
    assert(dim0 == A.dim0 && dim1 == A.dim1);
#endif
    int i, j;
    double *p;
    const double *q;

    for (i=0; i<dim1; i++)
    {
	p = v + i*ld;
        q = A.v + i*A.ld;
        for (j=0; j<dim0; j++)
            *p++ = *q++;
    }
}

void BlockVec::VecSetToZero()
{
    int i, j;
    double *p;

    for (i=0; i<dim1; i++)
    {
	p = v + i*ld;
        for (j=0; j<dim0; j++)
            *p++ = 0.0;
    }
}

void BlockVec::BlockCopy(const BlockVec& A)
{
#ifdef DEBUG
    assert(size0 == A.size0 && dim1 == A.dim1);
#endif
    int i, j;
    double *p;
    const double *q;

    for (i=0; i<dim1; i++)
    {
	p = v + i*ld;
        q = A.v + i*A.ld;
        for (j=0; j<size0; j++)
            *p++ = *q++;
    }
}

void BlockVec::BlockSetToZero()
{
    int i, j;
    double *p;

    for (i=0; i<dim1; i++)
    {
	p = v + i*ld;
        for (j=0; j<size0; j++)
            *p++ = 0.0;
    }
}
