//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#ifndef _LOCALMAT_H_
#define _LOCALMAT_H_

#include <iostream>

using namespace std;

#include "BlockVec.h"
class LocalPrecon;

class LocalMat
{
protected:
    inline void solve_is_mult(const BlockVec& B, BlockVec& X) const;

public:
    virtual ~LocalMat() {} // virtual, to make sure we call derived destructor
    virtual double *& Data() = 0; // ref to ptr to double
    virtual const double *Data() const = 0;

    virtual LocalMat *CreateEmpty() const = 0;
    virtual LocalMat *CreateInv(LocalPrecon&) const = 0;
    virtual void SetToZero(int, int) = 0;
    virtual void MatCopy(const LocalMat& A) = 0;
    virtual void Print(ostream&) const = 0;

    virtual void Mat_Trans(LocalMat *B) const = 0;
    virtual void Mat_Mat_Add(const LocalMat *B, LocalMat *C, 
        double alpha = 1.0) const = 0;
    virtual void Mat_Mat_Mult(const LocalMat *B, LocalMat *C, 
        double alpha = 1.0, double beta = 0.0) const = 0;
    virtual void Mat_Vec_Mult(const BlockVec& B, BlockVec& C,
        double alpha = 1.0, double beta = 0.0) const = 0;
    virtual void Mat_Trans_Vec_Mult(const BlockVec& B, BlockVec&C,
        double alpha = 1.0, double beta = 0.0) const = 0;
    virtual void Mat_Vec_Solve(const BlockVec& b, 
        BlockVec& x) const = 0;
    virtual void Mat_Trans_Vec_Solve(const BlockVec& b, 
        BlockVec& x) const = 0;
};

inline void LocalMat::solve_is_mult(const BlockVec& B, BlockVec& X) const
{
    if (&B != &X)
    {
        Mat_Vec_Mult(B, X);
    }
    else
    {
        // solves must allow solves in place
        BlockVec T(B, -1);
        Mat_Vec_Mult(T, X);
    }
}

typedef LocalMat *LocalMatp;

#endif // _LOCALMAT_H_
