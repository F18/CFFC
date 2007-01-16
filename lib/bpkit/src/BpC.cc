//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#include "BpC.h"
#include "BlockMat.h"
#include "BRelax.h"
#include "BILUK.h"
#include "BTIF.h"
#include "BPKIT.h"
#include "HBTMat.h"
#include "FGMRES.h"
#include "CSRMat.h"

#ifdef __cplusplus
extern "C" {
#endif

// prototypes for this file
static void set_localprecon(GlobalPrecon *M, int local, 
  double lparam1, double lparam2);

#if 0
// catch memory errors
static void freeStoreException()
{
    cout << endl; // flush standard output
    bperror("BpFortran: free store exhausted", 0);
}
#endif

void bpinitialize()
{
#if 0
    void (*_new_handler)();
    _new_handler = freeStoreException;
#endif

    // sun4 and Crays cannot initialize static members of class except this way
    HBTMat::rhs_resource = BpGetInt(HBTMAT_RHS, HBTMAT_RHS_RAND);
    HBTMat::transpose_resource = BpGetInt(HBTMAT_TRANSPOSE, 1);
    HBTMat::scale_resource = BpGetInt(HBTMAT_SCALE, 0);
    fgmres::dim           = BpGetInt(FGMRES_DIM, 20);
    fgmres::max_iter      = BpGetInt(FGMRES_MAX_ITER, 300);
    fgmres::tol           = BpGetDouble(FGMRES_TOL, 1.e-6);
    fgmres::print_history = BpGetInt(FGMRES_HISTORY, 0);
    CSRMat_APINVS::method     = BpGetInt(APINVS_METHOD,     1);
    CSRMat_APINVS::guess      = BpGetInt(APINVS_GUESS,      0);
    CSRMat_APINVS::nouter     = BpGetInt(APINVS_NOUTER,     1);
    CSRMat_APINVS::ninner     = BpGetInt(APINVS_NINNER,     1);
    CSRMat_APINVS::transpose  = BpGetInt(APINVS_TRANSPOSE,  1);
    CSRMat_APINVS::epsilon = BpGetDouble(APINVS_EPSILON, 0.0);
    CSRMat_APINV0::transpose = BpGetInt(APINV0_TRANSPOSE, 1);
    CSRMat_APINV_BANDED::transpose = BpGetInt(APINV_BANDED_TRANSPOSE, 1);
    CSRMat_RILUK::growth      = BpGetInt(RILUK_GROWTH, 2);
    CSRMat_RILUK::stab_thresh = BpGetDouble(RILUK_THRESH, 0.0);
    CSRMat_RILUK::stab_rel    = BpGetDouble(RILUK_REL, 1.0);
    CSRMat_RILUK::stab_abs    = BpGetDouble(RILUK_ABS, 0.0);
    CSRMat_RILUK::stat        = BpGetInt(RILUK_STAT, 0);
    BILUK::growth = BpGetInt(BILUK_GROWTH, 2);
    CSRMat_ILUT::permtol = BpGetDouble(ILUT_PERMTOL, 0.0);
    CSRMat_ILUT::pblock  = BpGetInt(ILUT_PBLOCK, 0);
    CSRMat_ILUT::stat    = BpGetInt(ILUT_STAT, 0);
}

void bpblockmatrix(
  void    **bmat,
  int     n,
  double  *a,
  int     *ja,
  int     *ia,
  int     nb,
  int     *kvst,
  int     type)
{
    // convert to 0-based indexing
    int i;
    for (i=0; i<=n; i++)
        ia[i]--;
    for (i=0; i<ia[n]; i++)
        ja[i]--;
    for (i=0; i<=nb; i++)
        kvst[i]--;

    HBTMat H(n, a, ja, ia, NULL, NULL, NULL);
    *bmat = (void *) new BlockMat(H, nb, kvst, (BlockType) type);

    // convert back to 1-based indexing
    for (i=0; i<=n; i++)
        ia[i]++;
    for (i=0; i<ia[n]; i++)
        ja[i]++;
    for (i=0; i<=nb; i++)
        kvst[i]++;
}

// C interface for conversion routine from block coordinate format
// assumes matrix is in 1-based indexing
void bpblockmatrix2(
  void   **bmat,
  int      nrow,
  int      nnz,
  int     *row,
  int     *col,
  double  *A,
  int      nb)
{
    *bmat = (void *) new BlockMat(nrow, nnz, row, col, A, nb);
}

void bpfreeblockmatrix(void *bmat)
{
    delete (BlockMat *) bmat;
}

void bppreconditioner(
        void    **precon,
  const void    *bmat,
  const int     global,
  const double  gparam1,
  const double  gparam2,
  const int     local,
  const double  lparam1,
  const double  lparam2)
{
    GlobalPrecon *M;
    BlockMat *B = (BlockMat *) bmat;

    switch (global)
    {
    case 0:
	M = new None;
	break;
    case 1:
	M = new BJacobi;
        set_localprecon(M, (int)local, lparam1, lparam2);
	((BJacobi *)M)->setup(*B);
        break;
    case 2:
	M = new BSOR;
        set_localprecon(M, (int)local, lparam1, lparam2);
        ((BSOR *)M)->setup(*B, (double)gparam1, (int)gparam2);
	break;
    case 3:
	M = new BSSOR;
	set_localprecon(M, (int)local, lparam1, lparam2);
	((BSSOR *)M)->setup(*B, (double)gparam1, (int)gparam2);
	break;
    case 4:
	M = new BILUK;
	set_localprecon(M, (int)local, lparam1, lparam2);
	((BILUK *)M)->setup(*B, (int)gparam1);
	break;
    case 5:
	M = new BTIF;
	set_localprecon(M, (int)local, lparam1, lparam2);
	((BTIF *)M)->setup(*B);
        break;
    default:
        bperror("BpC: no such global preconditioner", 0);
    }

    *precon = (void *) M;
}

void bpfreepreconditioner(void *precon)
{
    delete (GlobalPrecon *) precon;
}

void bpflexgmres(
  const void    *bmat,
        double  *x,
  const double  *rhs,
  const void    *precon,
  const int     dim,
  const int     max_iter,
  const double  tol)
{
    BlockMat *B = (BlockMat *) bmat;
    GlobalPrecon *M = (GlobalPrecon *) precon;

    fgmres f(dim, max_iter, tol);
    f.solve(*B, x, rhs, *M);
    f.status(cout);
}

// static functions

#define INT(d) ((int)(d+0.5)) // round to integer

static void set_localprecon(GlobalPrecon *M, int local_, 
  double lparam1, double lparam2)
{
    LocalPreconName local = (LocalPreconName) local_;

    switch (local)
    {
    case LP_LU:
	M->localprecon(local);
	break;
    case LP_INVERSE:
	M->localprecon(local);
	break;
    case LP_SVD:
	M->localprecon(local, lparam1, lparam2);
	break;
    case LP_RILUK:
	M->localprecon(local, INT(lparam1), lparam2);
	break;
    case LP_ILUT:
	M->localprecon(local, INT(lparam1), lparam2);
	break;
    case LP_APINV_TRUNC:
	M->localprecon(local, INT(lparam1));
	break;
    case LP_APINV_BANDED:
	M->localprecon(local, INT(lparam1));
	break;
    case LP_APINV0:
	M->localprecon(local);
	break;
    case LP_APINVS:
	M->localprecon(local, INT(lparam1), INT(lparam2));
	break;
    case LP_DIAG:
	M->localprecon(local);
	break;
    case LP_TRIDIAG:
	M->localprecon(local);
	break;
    case LP_SOR:
	M->localprecon(local, lparam1, INT(lparam2));
	break;
    case LP_SSOR:
	M->localprecon(local, lparam1, INT(lparam2));
	break;
    case LP_GMRES:
	M->localprecon(local, INT(lparam1), INT(lparam2));
	break;
    default:
        bperror("BpFortran: no such local preconditioner", 0);
    }
}

void bpmatvec(
  void    *bmat,
  int     nr,
  int     nc,
  const double *u,
  int     ldu,
  double *v,
  int     ldv)
{
    BlockMat *B = (BlockMat *) bmat;
    B->mult(nr, nc, u, ldu, v, ldv);
}

void bpapply(
  void    *prec,
  int     nr,
  int     nc,
  const double *u,
  int     ldu,
  double *v,
  int     ldv)
{
    GlobalPrecon *M = (GlobalPrecon *) prec;
    M->apply(nr, nc, u, ldu, v, ldv);
}

#ifdef __cplusplus
}
#endif
