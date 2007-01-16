//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

// A user designed global preconditioner for a 2x2 block matrix.
// The preconditioner is RILUK for the 1,1 block, and DIAG for the 2,2 block.

#include "LocalMat.h"
#include "BlockVec.h"
#include "BlockMat.h"
#include "B22.h"
#include "BPKIT.h"

B22::B22()
{
    diag0 = NULL;
    diag1 = NULL;
}

B22::~B22()
{
    delete diag0;
    delete diag1;
}

void B22::setup(const BlockMat& A)
{
    int j;
    LocalPrecon lp;

    Ap = &A;

    lp.name = LP_RILUK;
    lp.name = LP_DIAG;
    lp.iarg1 = 0;
    lp.darg1 = 0.0;
    j = (A.col_ind(0) == 0 ? 0:1);
    diag0 = A.val(j).CreateInv(lp);

    lp.name = LP_RILUK;
    lp.iarg1 = 0;
    lp.darg1 = 0.0;
    j = (A.col_ind(2) == 1 ? 2:3);
    diag1 = A.val(j).CreateInv(lp);
}

void B22::apply(int nr, int nc, const double *u, int ldu, 
    double *v, int ldv)
{
    BlockVec  V(nr, nc, v, ldv, &Ap->kvst_col(0));
    BlockVec  U(nr, nc, u, ldu, &Ap->kvst_col(0));

    diag0->Mat_Vec_Solve(U(0), V(0));
    diag1->Mat_Vec_Solve(U(1), V(1));
}

