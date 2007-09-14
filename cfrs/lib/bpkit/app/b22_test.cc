//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

// Example that uses a user-designed global preconditioner.
// The default output is:            600      4.933497030187e-05

#include "BlockMat.h"
#include "B22.h"
#include "HBTMat.h"
#include "FGMRES.h"

int main()
{
    int partit[3] = {0, 500, 1000};

    HBTMat H("SHERMAN1");
    double *rhs = H.get_rhs();
    double *x = H.get_guess();

    BlockMat B(H, 2, partit, CSR);
    B22 M;
    M.setup(B);

    fgmres f;
    f.solve(B, x, rhs, M);
    f.status(cout);

    return 0;
}
