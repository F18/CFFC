/*
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota
*/

#ifndef _BPFORTRAN_H_
#define _BPFORTRAN_H_

#include "arch.h"

#ifdef __cplusplus
extern "C" {
#endif
void F77NAME(bpinitialize)();

void F77NAME(blockmatrix)(
  integer8 *bmat,
  integer *n,
  double  *a,
  int     *ja,
  int     *ia,
  integer *nb,
  int     *kvst,
  integer *type);

void F77NAME(blockmatrix2)(
  integer8 *bmat,
  integer *nrow,
  integer *nnz,
  int     *row,
  int     *col,
  double  *A,
  integer *nb);

void F77NAME(freeblockmatrix)(
        integer8 *bmat);

void F77NAME(preconditioner)(
        integer8 *precon,
  const integer8 *bmat,
  const integer *global,
  const double  *gparam1,
  const double  *gparam2,
  const integer *local,
  const double  *lparam1,
  const double  *lparam2);

void F77NAME(freepreconditioner)(
        integer8 *precon);

void F77NAME(flexgmres)(
  const integer8 *bmat,
        double  *x,
  const double  *rhs,
  const integer8 *precon,
  const integer *dim,
  const integer *max_iter,
  const double  *tol);

void F77NAME(matvec)(
  integer8 *bmat,
  integer *nr,
  integer *nc,
  const double *u,
  integer *ldu,
  double *v,
  integer *ldv);

void F77NAME(apply)(
  integer8 *prec,
  integer *nr,
  integer *nc,
  const double *u,
  integer *ldu,
  double *v,
  integer *ldv);

#ifdef __cplusplus
}
#endif

#endif /* _BPFORTRAN_H_ */
