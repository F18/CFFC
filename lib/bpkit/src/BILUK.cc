//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#include "LocalMat.h"
#include "BlockVec.h"
#include "BlockMat.h"
#include "BILUK.h"
#include "BPKIT.h"
#include "SparseUtil.h"
#include "BpResource.h"

int BILUK::growth = BpGetInt(BILUK_GROWTH, 2);

BILUK::BILUK()
{
    diag = (LocalMat **) NULL;
     al = NULL;
    jal = NULL;
    ial = NULL;
     au = NULL;
    jau = NULL;
    iau = NULL;
}

BILUK::~BILUK()
{
    int i;

    if (diag != NULL)
        for (i=0; i<Ap->numrow(); i++)
            delete diag[i];
    delete [] diag;

    if (al != NULL)
        for (i=0; i<ial[Ap->numrow()]; i++)
            delete al[i];
    delete [] al;

    if (au != NULL)
        for (i=0; i<iau[Ap->numrow()]; i++)
            delete au[i];
    delete [] au;

    delete [] jal;
    delete [] ial;
    delete [] jau;
    delete [] iau;
}

// if levfill < 0, then use existing factorization for the pattern

void BILUK::setup(const BlockMat& A, int levfill)
{
    Ap = &A;

    if (levfill >= 0)
    {
	int nzl, nzu;

        delete jal;
        delete ial;
        delete jau;
        delete iau;

	allocate_ilu(levfill, Ap->numrow(), &nzl, &nzu, 
            &A.row_ptr(0), &A.col_ind(0), &ial, &jal, &iau, &jau, growth);

        symbolic_ilu(levfill, Ap->numrow(), &nzl, &nzu, 
            &A.row_ptr(0), &A.col_ind(0), ial, jal, iau, jau);

        // allocate for values
  	delete [] al;
  	delete [] au;
        al = new LocalMatp[nzl];
        au = new LocalMatp[nzu];
    }

    int i, j, k, kk, id, idd;

    int *marker = new int[Ap->numrow()];
    for (i=0; i<Ap->numrow(); i++)
        marker[i] = -1; // negative flag

    // full length work array of matrices
    LocalMat **row = new LocalMatp[A.numcol()];  // array of pointers to matrices

    diag = new LocalMatp[A.numrow()];

    for (i=0; i<A.numrow(); i++)
    {
	int neqr = A.kvst_row(i+1) - A.kvst_row(i);

        // scatter data structure of L and U
	j = 0;
        for (k=ial[i]; k<ial[i+1]; k++)
	{
            marker[jal[k]] = j;
	    row[j] = A.val(0).CreateEmpty();
	    row[j++]->SetToZero(neqr, A.kvst_col(jal[k]+1)-A.kvst_col(jal[k]));
	}
        for (k=iau[i]; k<iau[i+1]; k++)
	{
            marker[jau[k]] = j;
	    row[j] = A.val(0).CreateEmpty();
	    row[j++]->SetToZero(neqr, A.kvst_col(jau[k]+1)-A.kvst_col(jau[k]));
	}
        
        // scatter row of A
        for (k=A.row_ptr(i); k<A.row_ptr(i+1); k++)
            row[marker[A.col_ind(k)]]->MatCopy(A.val(k));

        LocalMat *mult = A.val(0).CreateEmpty();

        // eliminate the elements in L in order
        for (k=ial[i]; k<ial[i+1]; k++)
        {
            id = jal[k];

            // mult = row[id] / au[idiag[id]];
            row[marker[id]]->Mat_Mat_Mult(diag[id], mult);

            row[marker[id]]->MatCopy(*mult);

            for (kk=iau[id]+1; kk<iau[id+1]; kk++)
            {
                idd = jau[kk];
                if (marker[idd] >= 0) 
		{
                    // row[idd] = row[idd] - mult * au[kk];
                    mult->Mat_Mat_Mult(au[kk], row[marker[idd]], -1.0, 1.0);
		}
            }
        }

	delete mult;

        // gather resulting rows in L and U
        for (k=ial[i]; k<ial[i+1]; k++)
        {
            al[k] = row[marker[jal[k]]];
            marker[jal[k]] = -1;
        }
        for (k=iau[i]; k<iau[i+1]; k++)
        {
            au[k] = row[marker[jau[k]]];
            marker[jau[k]] = -1;
        }

        diag[i] = au[iau[i]]->CreateInv(local_precon);  // inverse of diagonal
    }

    delete [] row;
    delete [] marker;
}

void BILUK::apply(int nr, int nc, const double *u, int ldu, 
    double *v, int ldv)
{
    applyl(nr, nc, u, ldu, v, ldv);
    applyr(nr, nc, v, ldv, v, ldv);
}

void BILUK::applyl(int nr, int nc, const double *u, int ldu, 
    double *v, int ldv)
{
    int i, j;

    BlockVec  V(nr, nc, v, ldv, &Ap->kvst_col(0));
    BlockVec V2(nr, nc, v, ldv, &Ap->kvst_col(0));
    BlockVec  U(nr, nc, u, ldu, &Ap->kvst_col(0));
    V.VecCopy(U);

    // forward solve with lower triang factor (identity on diagonal assumed)

    for (i=0; i<Ap->numrow(); i++)
    {
        for (j=ial[i]; j<ial[i+1]; j++)
        {
            // V(i) = V(i) - al[j] * V(jal[j])
            al[j]->Mat_Vec_Mult(V(jal[j]), V2(i), -1.0, 1.0);
        }
    }

}

void BILUK::applyr(int nr, int nc, const double *u, int ldu, 
    double *v, int ldv)
{
    int i, j;

    BlockVec  V(nr, nc, v, ldv, &Ap->kvst_col(0));
    BlockVec V2(nr, nc, v, ldv, &Ap->kvst_col(0));
    BlockVec  U(nr, nc, u, ldu, &Ap->kvst_col(0));
    V.VecCopy(U);

    // backward solve with upper triang factor

    for (i=Ap->numrow()-1; i>=0; i--)
    {
        for (j=iau[i]+1; j<iau[i+1]; j++)
        {
            // V(i) = V(i) - au[j] * V(jau[j])
            au[j]->Mat_Vec_Mult(V(jau[j]), V2(i), -1.0, 1.0);
        }
        diag[i]->Mat_Vec_Solve(V(i), V(i));
    }
}

void BILUK::multiply(int nr, int nc, const double *u, int ldu, 
    double *v, int ldv)
{
    int i, j;

    BlockVec U(nr, nc, u, ldu, &Ap->kvst_col(0));
    BlockVec V(nr, nc, v, ldv, &Ap->kvst_col(0));
    BlockVec T(U);
    V.VecSetToZero();

    // multiply with upper triang factor

    for (i=0; i<Ap->numrow(); i++)
    {
        for (j=iau[i]; j<iau[i+1]; j++)
        {
            // V(i) = V(i) + au[j] * U(jau[j])
            au[j]->Mat_Vec_Mult(U(jau[j]), V(i), 1.0, 1.0);
        }
    }

    // multiply with lower triang factor (unit diagonal assumed)
    T.VecCopy(V);

    for (i=0; i<Ap->numrow(); i++)
    {
        for (j=ial[i]; j<ial[i+1]; j++)
        {
            // V(i) = V(i) + al[j] * temp(jal[j])
            al[j]->Mat_Vec_Mult(T(jal[j]), V(i), 1.0, 1.0);
        }
    }
}
