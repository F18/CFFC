//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#ifndef _BLOCKVEC_H_
#define _BLOCKVEC_H_

#include <cstdlib>
#ifndef NULL
#define NULL 0
#endif
#ifdef DEBUG
#include <cassert>
#endif

using namespace std;

class BlockVec
{
private:
    double *base;
    const int *partit;
    int owndata;

public:
    double *v;
    int     dim0;
    int     dim1;
    int     ld;
    int     size0;

    BlockVec& operator()(int i)
    {
#ifdef DEBUG
            assert(partit != NULL);
            assert(partit[i] < partit[i+1]);
            assert(partit[i+1] <= dim0);
#endif
	    v = base + partit[i];
	    size0 = partit[i+1] - partit[i];
	    return *this;
    }

    BlockVec(int nr, int nc, const double *a, int lda, const int *partitioning)
    {
        dim0 = nr;
        dim1 = nc;
        ld = lda;
	size0 = nr;
        base = (double *)a;
        v = (double *)a;
        partit = partitioning;
	owndata = 0;
    }

    BlockVec(const BlockVec& A);
    BlockVec(const BlockVec& A, int i);

   ~BlockVec()
    {
	if (owndata)
	    delete [] base;
	base = NULL;
    }

    void VecCopy(const BlockVec& A);
    void VecSetToZero();

    void BlockCopy(const BlockVec& A);
    void BlockSetToZero();
};

#endif // _BLOCKVEC_H_
