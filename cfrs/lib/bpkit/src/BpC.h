/*
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota
*/

#ifndef _BPC_H_
#define _BPC_H_

#define BP_DENSE 	0
#define BP_SPARSE 	1

#define BP_NONE 	0
#define BP_BJACOBI 	1
#define BP_BSOR 	2
#define BP_BSSOR 	3
#define BP_BILUK 	4
#define BP_BTIF 	5

#define BP_LU 	 	1
#define BP_INVERSE 	2
#define BP_SVD 	 	3
#define BP_RILUK 	4
#define BP_ILUT 	5
#define BP_APINV_TRUNC 	6
#define BP_APINV_BANDED 7
#define BP_APINV0 	8
#define BP_APINVS 	9
#define BP_DIAG 	10
#define BP_TRIDIAG 	11
#define BP_SOR 	 	12
#define BP_SSOR 	13
#define BP_GMRES 	14

#ifdef __cplusplus
extern "C" {
#endif

void bpinitialize();

void bpblockmatrix(
  void    **bmat,
  int     n,
  double  *a,
  int     *ja,
  int     *ia,
  int     nb,
  int     *kvst,
  int     type);

void bpblockmatrix2(
  void   **bmat,
  int      nrow,
  int      nnz,
  int     *row,
  int     *col,
  double  *A,
  int      nb);

void bpfreeblockmatrix(
        void    *bmat);

void bppreconditioner(
        void    **precon,
  const void    *bmat,
  const int     global,
  const double  gparam1,
  const double  gparam2,
  const int     local,
  const double  lparam1,
  const double  lparam2);

void bpfreepreconditioner(
        void    *precon);

void bpflexgmres(
  const void    *bmat,
        double  *x,
  const double  *rhs,
  const void    *precon,
  const int     dim,
  const int     max_iter,
  const double  tol);

void bpmatvec(
  void    *bmat,
  int     nr,
  int     nc,
  const double *u,
  int     ldu,
  double *v,
  int     ldv);

void bpapply(
  void    *prec,
  int     nr,
  int     nc,
  const double *u,
  int     ldu,
  double *v,
  int     ldv);

#ifdef __cplusplus
}
#endif

#endif /* _BPC_H_ */
