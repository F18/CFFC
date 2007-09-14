//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#include "LocalMat.h"
#include "BlockVec.h"
#include "BlockMat.h"
#include "BTIF.h"
#include "BPKIT.h"

BTIF::BTIF()
{
    diag = (LocalMat **) NULL;
}

BTIF::~BTIF()
{
    int i;

    if (diag != NULL)
        for (i=0; i<Ap->numrow(); i++)
            delete diag[i];
    delete [] diag;
}

void BTIF::setup(const BlockMat& A)
{
    int i;

    // Check that the matrix is block tridiagonal.
    int fail = 0;
    int k;
    const int *ja = &A.col_ind(0);
    if (ja[0] != 0 || ja[1] != 1) fail = 1;
    for (i=1, k=2; i<A.numrow()-1; i++, k+=3)
	if (ja[k] != i-1 || ja[k+1] != i || ja[k+2] != i+1) fail = 1;
    if (ja[k++] != A.numrow()-2 || ja[k] != A.numrow()-1) fail = 1;
    if (fail)
    {
	bperror("BTIF: matrix is not block tridiagonal", 0);
    }

    Ap = &A;
    diag = new LocalMatp[A.numrow()];

    diag[0] = A.val(0).CreateInv(local_precon);

    for (i=1; i<A.numrow(); i++)
    {
        LocalMat *temp1 = A.val(0).CreateEmpty();
        LocalMat *temp2 = A.val(0).CreateEmpty();
        LocalMat *temp3 = A.val(0).CreateEmpty();

        A.val(3*i-1).Mat_Mat_Mult(diag[i-1], temp1);   // mult with subdiag

        temp1->Mat_Mat_Mult(&A.val(3*i-2), temp2);     // mult with superdiag

        A.val(3*i).Mat_Mat_Add(temp2, temp3, -1.0);

        diag[i] = temp3->CreateInv(local_precon);

        delete temp3;
        delete temp2;
        delete temp1;
    }
}

void BTIF::apply(int nr, int nc, const double *b, int ldb,
    double *x, int ldx)
{
    // b and x may be the same memory space

    int i;
    BlockVec T(nr, nc, b, ldb, &Ap->kvst_col(0));
    BlockVec B(T);
    BlockVec X(nr, nc, x, ldx, &Ap->kvst_col(0));

    // X(0) = diag(0) * B(0)
    diag[0]->Mat_Vec_Mult(B(0), X(0));

    for (i=1; i<Ap->numrow(); i++)
    {
        // B(i) = B(i) - a(i,i-1) * X(i-1)
        Ap->val(3*i-1).Mat_Vec_Mult(X(i-1), B(i), -1.0, 1.0);

        // X(i) = diag(i) * B(i)
        diag[i]->Mat_Vec_Mult(B(i), X(i));
    }

    for (i=Ap->numrow()-2; i>=0; i--)
    {
	// B(i) = a(i,i+1) * X(i+1)
	Ap->val(3*i+1).Mat_Vec_Mult(X(i+1), B(i));

	// X(i) = X(i) - diag(i) * B(i)
	diag[i]->Mat_Vec_Mult(B(i), X(i), -1.0, 1.0);
    }
}
