/*
 *    Example program that calls BPKIT library from C.
 *    Hardcoded for SHERMAN1 matrix from HB collection.
 *    With the original .BpResource, the output should be:
 *          88      9.520099361807e-09
 */

#include "BpC.h"
#include "BpFortran.h" // only used for FORTRAN readhb

#define nmax  1000
#define nzmax 3750

int main()
{
      int       ia[nmax+1], ja[nzmax], kvst[nmax+1];
      double    a[nzmax], rhs[nmax], sol[nmax], exact[nmax];
      double    temp[nmax*3];
      integer   i, n, nnz, nrhs;
      void      *bmat, *precon;

      integer   eight=8;
      int       dim = 20, maxits = 600;
      double    tol = 1.e-8;

      /* read matrix */
      n = nmax;
      nnz = nzmax;
      nrhs = nmax*3;

      F77NAME(readhb)("SHERMAN1", &eight, &n, &nnz, &nrhs, a, ja, ia,
          rhs, sol, exact, temp);

      for (i=0; i<n; i++)
      {
          rhs[i] = temp[i];  /* right-hand side */
          sol[i] = 0.0;
      }

      for (i=0; i<11; i++)
          kvst[i] = i*100 + 1;

      bpinitialize();
      bpblockmatrix(&bmat, n, a, ja, ia, 10, kvst, BP_SPARSE);
      bppreconditioner(&precon, bmat, BP_BJACOBI, 0.0, 0.0, 
				       BP_INVERSE, 0.0, 0.0);
      bpflexgmres(bmat, sol, rhs, precon, dim, maxits, tol);

      return 0;
}
