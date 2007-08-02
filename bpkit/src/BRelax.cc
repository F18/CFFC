//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#include "LocalMat.h"
#include "BlockVec.h"
#include "BlockMat.h"
#include "BRelax.h"
#include "BPKIT.h"

void None::apply(int nr, int nc, const double *u, int ldu, 
    double *v, int ldv)
{
    int i, j;
    double *p;
    const double *q;

    for (i=0; i<nc; i++)
    {
        p = v + i*ldv;
        q = u + i*ldu;
        for (j=0; j<nr; j++)
            *p++ = *q++;
    }
}

BJacobi::BJacobi()
{
    Ap = (BlockMat *) NULL;
    diag = (LocalMat **) NULL;
}

BJacobi::~BJacobi()
{
    if (Ap != NULL)
        for (int i=0; i<Ap->numrow(); i++)
            delete diag[i];
    delete [] diag;
}

void BJacobi::setup(const BlockMat& A)
{
    int i, j;
    int got_diag;

    Ap = &A;
    diag = new LocalMatp[A.numrow()];

    // search for diagonal blocks
    for (i=0; i<A.numrow(); i++)
    {
        got_diag = FALSE;
        for (j=A.row_ptr(i); j<A.row_ptr(i+1); j++)
        {
            if (A.col_ind(j) == i)
            {
                diag[i] = A.val(j).CreateInv(local_precon);
                got_diag = TRUE;
            }
        }
        if (!got_diag)
            bperror("BRelax: matrix does not have diagonal block", i);
    }
}

void BJacobi::apply(int nr, int nc, const double *u, int ldu,
    double *v, int ldv)
{
    BlockVec U(nr, nc, u, ldu, &Ap->kvst_row(0));
    BlockVec V(nr, nc, v, ldv, &Ap->kvst_row(0));

    for (int i=0; i<Ap->numrow(); i++)
    {
        diag[i]->Mat_Vec_Solve(U(i), V(i));
    }
}

// note:
// assumes each matrix row is stored such that lower triangular elements,
// diagonal elements, upper triangular elements are stored in order

BSOR_Base::BSOR_Base()
{
    Ap = (BlockMat *) NULL;
    diag = (LocalMat **) NULL;
    idiag = (int *) NULL;
}

BSOR_Base::~BSOR_Base()
{
    if (Ap != NULL)
        for (int i=0; i<Ap->numrow(); i++)
            delete diag[i];
    delete [] diag;
    delete [] idiag;
}

void BSOR_Base::setup(const BlockMat& A, double omega, int iterations)
{
    int i, j, nrow;
    int got_diag;

    omega_ = omega;
    iterations_ = iterations;
    Ap = &A;
    nrow = A.numrow();
    diag = new LocalMatp[nrow];
    idiag = new int[nrow];

    // search for diagonal blocks
    for (i=0; i<nrow; i++)
    {
        got_diag = FALSE;
        for (j=A.row_ptr(i); j<A.row_ptr(i+1); j++)
        {
            if (A.col_ind(j) == i)
            {
                diag[i] = A.val(j).CreateInv(local_precon);
                idiag[i] = j;
                got_diag = TRUE;
            }
        }
        if (!got_diag)
            bperror("BRelax: matrix does not have diagonal block", i);
    }
}

// c = alpha * a + beta * b
// works on Blocks within a BlockVec

static void gaxpy(const double& alpha, const BlockVec& a,
    const double& beta, const BlockVec& b, BlockVec& c)
{
    double *ap, *bp, *cp;
    for (int j=0; j<a.dim1; j++)
    {
        ap = a.v + j * a.ld;
        bp = b.v + j * b.ld;
        cp = c.v + j * c.ld;

        for (int i=0; i<a.size0; i++) 
	    *cp++ = alpha * *ap++ + beta * *bp++;
            // c(i,j) = alpha*a(i,j) + beta*b(i,j);
    }
}

void BSOR::apply(int nr, int nc, const double *u, int ldu,
    double *v, int ldv)
{
    const int *ia = &Ap->row_ptr(0);
    const int *ja = &Ap->col_ind(0);
    int it, i, j;

#if 0
    BlockVec  V(nr, nc, v, ldv, &Ap->kvst_col(0));
    BlockVec V2(nr, nc, v, ldv, &Ap->kvst_col(0));
    BlockVec  U(nr, nc, u, ldu, &Ap->kvst_col(0));
    V.VecCopy(U);

    // Specialized code for first step.

    for (i=0; i<Ap->numrow(); i++)
    {
        for (j=ia[i]; j<idiag[i]; j++)
        {
            // V(i) = V(i) - omega_ a[j] * V(ja[j])
            Ap->val(j).Mat_Vec_Mult(V(ja[j]), V2(i), -omega_, 1.0);
        }
        diag[i]->Mat_Vec_Solve(V(i), V(i));
    }

    // After first step....
#endif
    BlockVec  V(nr, nc, v, ldv, &Ap->kvst_col(0));
    BlockVec  U(nr, nc, u, ldu, &Ap->kvst_col(0));
    V.VecSetToZero();

    for (it=1; it<=iterations_; it++)
    {
        for (i=0; i<Ap->numrow(); i++)
        {
            BlockVec temp(U,i);

            for (j=ia[i]; j<idiag[i]; j++)
            {
                // temp = temp - a[j] * v[ja[j]];
                Ap->val(j).Mat_Vec_Mult(V(ja[j]), temp, -1.0, 1.0);
            }
            for (j=idiag[i]+1; j<ia[i+1]; j++)
            {
                // temp = temp - a[j] * v[ja[j]];
                Ap->val(j).Mat_Vec_Mult(V(ja[j]), temp, -1.0, 1.0);
            }

            diag[i]->Mat_Vec_Solve(temp, temp);

            // v[i] = (1.0-omega_) * v[i] + omega_ * temp
            gaxpy(1.0-omega_, V(i), omega_, temp, V(i));
        }
    }
}

void BSSOR::apply(int nr, int nc, const double *u, int ldu,
    double *v, int ldv)
{
    const int *ia = &Ap->row_ptr(0);
    const int *ja = &Ap->col_ind(0);
    int it, i, j;

#if 0
    BlockVec  V(nr, nc, v, ldv, &Ap->kvst_col(0));
    BlockVec V2(nr, nc, v, ldv, &Ap->kvst_col(0));
    BlockVec  U(nr, nc, u, ldu, &Ap->kvst_col(0));
    V.VecCopy(U);

    // Specialized code for first step.

    // lower sweep
    for (i=0; i<Ap->numrow(); i++)
    {
        for (j=ia[i]; j<idiag[i]; j++)
        {
            // V(i) = V(i) - omega_ a[j] * V(ja[j])
            Ap->val(j).Mat_Vec_Mult(V(ja[j]), V2(i), -omega_, 1.0);
        }
        diag[i]->Mat_Vec_Solve(V(i), V(i));
    }

    // multiply by diagonal blocks
    for (i=0; i<Ap->numrow(); i++)
    {
        // V(i) = diag[i] * V(i)
        // this cannot be done in place
        BlockVec y(V,i); // make a copy of V(i)

        Ap->val(idiag[i]).Mat_Vec_Mult(y, V(i));
    }

    // upper sweep
    for (i=Ap->numrow()-1; i>=0; i--)
    {
        for (j=idiag[i]+1; j<ia[i+1]; j++)
        {
            // V(i) = V(i) - omega_ a[j] * V(ja[j])
            Ap->val(j).Mat_Vec_Mult(V(ja[j]), V2(i), -omega_, 1.0);
        }
        diag[i]->Mat_Vec_Solve(V(i), V(i));
    }

    // After first step....
#endif

    BlockVec  V(nr, nc, v, ldv, &Ap->kvst_col(0));
    BlockVec  U(nr, nc, u, ldu, &Ap->kvst_col(0));
    V.VecSetToZero();

    for (it=1; it<=iterations_; it++)
    {
        for (i=0; i<Ap->numrow(); i++)
        {
            BlockVec temp(U,i);

            for (j=ia[i]; j<idiag[i]; j++)
            {
                // temp = temp - a[j] * v[ja[j]];
                Ap->val(j).Mat_Vec_Mult(V(ja[j]), temp, -1.0, 1.0);
            }
            for (j=idiag[i]+1; j<ia[i+1]; j++)
            {
                // temp = temp - a[j] * v[ja[j]];
                Ap->val(j).Mat_Vec_Mult(V(ja[j]), temp, -1.0, 1.0);
            }

            diag[i]->Mat_Vec_Solve(temp, temp);

            // v[i] = (1.0-omega_) * v[i] + omega_ * temp
            gaxpy(1.0-omega_, V(i), omega_, temp, V(i));
        }

        for (i=Ap->numrow()-1; i>=0; i--)
        {
            BlockVec temp(U,i);

            for (j=ia[i]; j<idiag[i]; j++)
            {
                // temp = temp - a[j] * v[ja[j]];
                Ap->val(j).Mat_Vec_Mult(V(ja[j]), temp, -1.0, 1.0);
            }
            for (j=idiag[i]+1; j<ia[i+1]; j++)
            {
                // temp = temp - a[j] * v[ja[j]];
                Ap->val(j).Mat_Vec_Mult(V(ja[j]), temp, -1.0, 1.0);
            }

            diag[i]->Mat_Vec_Solve(temp, temp);

            // v[i] = (1.0-omega_) * v[i] + omega_ * temp
            gaxpy(1.0-omega_, V(i), omega_, temp, V(i));
        }
    }
}
