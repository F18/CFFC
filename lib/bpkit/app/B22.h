//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#ifndef _B22_H_
#define _B22_H_

#include "GlobalPrecon.h"

class LocalMat;
class BlockMat;

class B22 : public GlobalPrecon
{
private:
    const BlockMat *Ap;
    LocalMat *diag0;
    LocalMat *diag1;

public:
    B22();
   ~B22();
 
    void setup(const BlockMat& A);
    void apply (int, int, const double *, int, double *, int);
};

#endif // _B22_H_
