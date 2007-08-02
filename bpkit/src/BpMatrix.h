//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#ifndef _BPMATRIX_H_
#define _BPMATRIX_H_

class BpMatrix
{
public:
    virtual ~BpMatrix() {}

    virtual int dimrow() const = 0;
    virtual int dimcol() const = 0;

    virtual void mult(int, int, const double *, int, double *, int) const {}
    virtual void trans_mult (int, int, const double *, int, double *, int) const
      {}
};

#endif // _BPMATRIX_H_
