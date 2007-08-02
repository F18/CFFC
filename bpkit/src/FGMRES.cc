//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#include <cmath>
 
using namespace std;

#include "FGMRES.h"
#include "BpMatrix.h"
#include "BpPrecon.h"
#include "BpResource.h"
#include "blas1.h"

// simulate 2-D arrays at the cost of some arithmetic
#define V(i) (&V[(i)*n])
#define W(i) (&W[(i)*n])
#define H(i,j) (H[(j)*m1+(i)])
#define ABS(x)   ((x)<0 ? (-(x)) : (x))

int    fgmres::dim           = BpGetInt(FGMRES_DIM, 20);
int    fgmres::max_iter      = BpGetInt(FGMRES_MAX_ITER, 300);
double fgmres::tol           = BpGetDouble(FGMRES_TOL, 1.e-6);
int    fgmres::print_history = BpGetInt(FGMRES_HISTORY, 0);

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

fgmres::fgmres()
{
    H  = new double[dim*(dim+1)];
}

fgmres::fgmres(int dim_, int max_iter_, double tol_)
{
    // override defaults and resource values
    dim = dim_;
    max_iter = max_iter_;
    tol = tol_;

    H  = new double[dim*(dim+1)];
}

fgmres::~fgmres()
{
    delete [] H;
}

void fgmres::solve(const BpMatrix &A, double *x, const double *b, BpPrecon &M)
{
    int n = A.dimrow();
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
    double *W  = new double[n*dim];

    iter = 0;
    do
    {
        // compute initial residual and its norm
        A.mult(n, 1, x, n, V(0), n);                        // V(0) = A*x
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

            M.apply(n, 1, V(i), n, W(i), n);
            A.mult(n, 1, W(i), n, V(i+1), n);

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
            if (print_history)
                status(cout);
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
            // x = x + s[j] * W(j)
            F77NAME(daxpy)(&N, &s[j], W(j), &inc, x, &inc); 
        }
    }
    while (rel_resid > tol && iter+1 <= max_iter);

    // compute exact residual norm reduction
    A.mult(n, 1, x, n, V(0), n);                        // V(0) = A*x
    temp = -1.0;
    F77NAME(daxpy)(&N, &temp, b, &inc, V(0), &inc);     // V(0) = V(0) - b
    beta = F77NAME(dnrm2)(&N, V(0), &inc);              // beta = norm(V(0))
    rel_resid = beta / resid0;

    delete [] s;
    delete [] cs;
    delete [] sn;
    delete [] V;
    delete [] W;
}

// Modified by C. P. T. Groth
//#ifdef CRAY
#include <stdio.h>
//#endif

void fgmres::status(ostream& os)
{
#if 1
    //cout << iter << endl;
    fprintf(stderr, "%14d    %20.12g\n", iter, rel_resid);
#else
    long olda = os.setf(ios::right,ios::adjustfield);
    long oldf = os.setf(ios::scientific,ios::floatfield);
    int oldp = os.precision(12);

    os.width(14);
    os << iter << "    "; 
    os.width(20);
    os << rel_resid << endl;

    os.setf(olda,ios::adjustfield);
    os.setf(oldf,ios::floatfield);
    os.precision(oldp);
#endif
}
