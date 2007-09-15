//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#include <cmath> 

using namespace std;

#include "CSRMat.h"
#include "blas1.h"
#include "spblas.h"

// simulate 2-D arrays at the cost of some arithmetic
#define V(i) (&V[(i)*n])
#define W(i) (&W[(i)*n])
#define H(i,j) (H[(j)*m1+(i)])
#define ABS(x)   ((x)<0 ? (-(x)) : (x))

static void 
GeneratePlaneRotation(double &dx, double &dy, double &cs, double &sn)
{
  if (dy == 0.0) {
    cs = 1.0;
    sn = 0.0;
  } else if (ABS(dy) > ABS(dx)) {
    double temp = dx / dy;
    sn = 1.0 / sqrt( 1.0 + temp*temp );
    cs = temp * sn;
  } else {
    double temp = dy / dx;
    cs = 1.0 / sqrt( 1.0 + temp*temp );
    sn = temp * cs;
  }
}

static void ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn)
{
  double temp  =  cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = temp;
}

CSRMat_GMRES::CSRMat_GMRES(const CSRMat& A,
    const int iterations, const double threshold)
{
    Ap = &A;
    iterations_ = iterations;
    threshold_ = threshold;
}

void CSRMat_GMRES::Mat_Vec_Solve(const BlockVec& B, BlockVec& X) const
{
    // this is single-vector GMRES
    assert (B.dim1 == 1);

    // make a copy of the right-hand side, in case B.v == X.v (solve in place)
    BlockVec B1(B, -1);
    const double *b = B1.v;

    // zero initial guess
    X.BlockSetToZero();
    double *x = X.v;

    // This is modified from the FGMRES object.

    int max_iter = iterations_;
    int dim = iterations_;
    double tol = threshold_;
    int iter;
    double rel_resid;

    int n = Ap->numrow();
    int m1 = dim+1; // used inside H macro
    int i, j, k;
    double beta, resid0;

    integer N = n;
    integer inc = 1;  // vector stride is always 1
    doublereal temp;

    double *s  = new double[dim+1];
    double *cs = new double[dim];
    double *sn = new double[dim];

    double *V  = new double[n*(dim+1)];
    double *H  = new double[dim*(dim+1)];

    int descra[9];
    descra[0] = 0;
    descra[1] = 0;
    descra[2] = 0;

    iter = 0;
    do
    {
        // compute initial residual and its norm

	F77NAME(dcsrmm) (0, n, 1, n, 1.0,
	    descra, Ap->a, Ap->ja, Ap->ia, x+1, n,
	    0.0, V(0), n, NULL, 0);                         // V(0) = A*x

        temp = -1.0;
        F77NAME(daxpy)(&N, &temp, b, &inc, V(0), &inc);     // V(0) = V(0) - b
        beta = F77NAME(dnrm2)(&N, V(0), &inc);              // beta = norm(V(0))
        temp = -1.0/beta;
        F77NAME(dscal)(&N, &temp, V(0), &inc);              // V(0) = -V(0)/beta

        // save very first residual norm
        if (iter == 0)
            resid0 = beta;

        for (i = 1; i < dim+1; i++)
            s[i] = 0.0;
        s[0] = beta;

        i = -1;
        do
        {
            i++;
            iter++;

	    F77NAME(dcsrmm) (0, n, 1, n, 1.0,
	        descra, Ap->a, Ap->ja, Ap->ia, V(i)+1, n,
	        0.0, V(i+1), n, NULL, 0);                    // V(i+1) = A*V(i)

            for (k = 0; k <= i; k++)
            {
                H(k, i) = F77NAME(ddot)(&N, V(i+1), &inc, V(k), &inc);
                temp = -H(k,i);
                // V(i+1) -= H(k, i) * V(k);
                F77NAME(daxpy)(&N, &temp, V(k), &inc, V(i+1), &inc);
            }

            H(i+1, i) = F77NAME(dnrm2)(&N, V(i+1), &inc);
            temp = 1.0 / H(i+1, i);
            // V(i+1) = V(i+1) / H(i+1, i)
            F77NAME(dscal)(&N, &temp, V(i+1), &inc);

            for (k = 0; k < i; k++)
                ApplyPlaneRotation(H(k,i), H(k+1,i), cs[k], sn[k]);
          
            GeneratePlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
            ApplyPlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
            ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);
          
            rel_resid = ABS(s[i+1]) / resid0;
            if (rel_resid <= tol)
                break;
        }
        while (i+1 < dim && iter+1 <= max_iter);

        // solve upper triangular system in place
        for (j = i; j >= 0; j--)
        {
            s[j] /= H(j,j);
            for (k = j-1; k >= 0; k--)
                s[k] -= H(k,j) * s[j];
        }

        // update the solution
        for (j = 0; j <= i; j++)
        {
            // x = x + s[j] * V(j)
            F77NAME(daxpy)(&N, &s[j], V(j), &inc, x, &inc); 
        }
    }
    while (rel_resid > tol && iter+1 <= max_iter);

    delete [] s;
    delete [] cs;
    delete [] sn;
    delete [] V;
    delete [] H;
}
