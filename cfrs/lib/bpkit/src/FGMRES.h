//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#ifndef _FGMRES_H_
#define _FGMRES_H_

#include <iostream>

using namespace std;

class BpMatrix;
class BpPrecon;

class fgmres
{
private:
    int iter;
    double rel_resid;

    double *H; // Hessenberg matrix

public:
    fgmres();
    fgmres(int dim, int max_iter, double tol);
   ~fgmres();
    void   solve(const BpMatrix &A, double *x, const double *b, BpPrecon &M);
    int    get_iter() {return iter;}
    double get_rel_resid_norm() {return rel_resid;}
    void   status(ostream& os);

    static int dim;
    static int max_iter;
    static double tol;
    static int print_history;
};

// resource names
#define FGMRES_DIM         "fgmres.dim"
#define FGMRES_MAX_ITER    "fgmres.max_iter"
#define FGMRES_TOL         "fgmres.tol"
#define FGMRES_HISTORY     "fgmres.history"

#endif // _FGMRES_H_
