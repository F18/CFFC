//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

// With the original .BpResource, the output should be: 24 9.994607604203e-09

#include "BlockMat.h"
#include "BRelax.h"
#include "HBTMat.h"
#include "FGMRES.h"

int main()
{
    HBTMat H("SHERMAN1");      // read Harwell-Boeing file from disk
    double *rhs = H.get_rhs();
    double *x = H.get_guess();

    BlockMat B(H, 100, CSR);   // block matrix with sparse 100x100 blocks

    BSSOR M;
    M.localprecon(LP_INVERSE);
    M.setup(B, 1.0, 1);        // BSSOR(omega=1, iterations=1)

    fgmres f(20, 100, 1.e-8);  // Kry_dim=20, max_iter=100, tol=1.e-8
    f.solve(B, x, rhs, M);
    f.status(cout);            // solution should be reached in 24 steps

    return 0;
}
