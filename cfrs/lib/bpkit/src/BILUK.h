//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#ifndef _BILUK_H_
#define _BILUK_H_

#include "GlobalPrecon.h"

class LocalMat;
class BlockMat;

class BILUK : public GlobalPrecon
{
private:
    const BlockMat *Ap;

    LocalMat **diag;  // inverse or factors of diagonal blocks

    LocalMat **al;    // lower triangular factor (strict lower part stored)
    int      *jal;
    int      *ial;
    LocalMat **au;    // upper triangular factor
    int      *jau;
    int      *iau;

public:
    BILUK();
   ~BILUK();
 
    void setup(const BlockMat& A, int levfill);
    void apply (int, int, const double *, int, double *, int);
    void applyr(int, int, const double *, int, double *, int);
    void applyl(int, int, const double *, int, double *, int);

    void multiply(int, int, const double *, int, double *, int);

    static int growth;
};

// resource names

#define BILUK_GROWTH "BILUK.growth"

#endif // _BILUK_H_
