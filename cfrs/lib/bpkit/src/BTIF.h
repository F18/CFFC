//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#ifndef _BTIF_H_
#define _BTIF_H_

#include "GlobalPrecon.h"

class LocalMat;
class BlockMat;

class BTIF : public GlobalPrecon
{
private:
    const BlockMat *Ap;
    LocalMat **diag;

public:
    BTIF();
   ~BTIF();
 
    void setup(const BlockMat& A);
    void apply (int, int, const double *, int, double *, int);
};

#endif // _BTIF_H_
