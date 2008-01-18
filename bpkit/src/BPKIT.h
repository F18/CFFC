//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#ifndef _BPKIT_H_
#define _BPKIT_H_

// added by Kalvin Tsang.
#include <cstdlib>
#include <iostream>

using namespace std;

#include "arch.h"

#ifdef CRAY
extern "C" { void exit(int); }
#endif

#define bperror(msg,code) \
{ cerr << "BPKIT: " << msg << "    CODE " << code << endl; exit(1); }

#define MIN(x,y) ((x)<(y) ? (x) : (y))
#define MAX(x,y) ((x)>(y) ? (x) : (y))
#define ABS(x)   ((x)<0 ? (-(x)) : (x))
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif
#define SGN(x) ((x)<0.0 ? -1.0 : 1.0)

// miscellaneous prototypes for FORTRAN functions

extern "C"
{
void F77NAME(dgetri)(integer *n, doublereal *A, integer *lda, integer *ipiv,
    doublereal *work, integer *lwork, integer *info);

void F77NAME(dgesvd) (char *, char *, integer *, integer *, const doublereal *,
    integer *, doublereal *, doublereal *, integer *, doublereal *, integer *,
    doublereal *, integer *, integer *);

void F77NAME(dgeev) (char *, char *, integer *, doublereal *, integer *,
    doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *,
    doublereal *, integer *, integer *);

void F77NAME(bpilut)(int *n, doublereal *a, int *ja, int *ia, int *lfil,
    doublereal *tol, doublereal *alu, int *jlu, int *ju,
    int *iwk, doublereal *w, int *jw, int *ierr);

void F77NAME(bpilutp)(int *n, doublereal *a, int *ja, int *ia, int *lfil,
    doublereal *tol, doublereal *permtol, int *pblock,
    doublereal *alu, int *jlu, int *ju,
    int *iwk, doublereal *w, int *jw, int *iperm, int *ierr);

void F77NAME(bplusol)(const int *n, const doublereal *y, doublereal *x,
    const doublereal *alu, const int *jlu, const int *ju);

void F77NAME(bpilutstat)(int *n, doublereal *a, int *ja, int *ia,
    doublereal *alu, int *jlu, int *ju, int *levnum);

void F77NAME(readhb)(const char *filename, int *flen,
    int *n, int *nnz, int *nrhs,
    doublereal *a, int *ja, int *ia,
    doublereal *rhs, doublereal *guess, doublereal *exact,
    doublereal *rhsdata);
}

#endif // _BPKIT_H_
