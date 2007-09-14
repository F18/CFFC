//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#ifndef _BRELAX_H_
#define _BRELAX_H_

#include "GlobalPrecon.h"

class LocalMat;
class BlockMat;

class None : public GlobalPrecon
{
public:
    void setup(const BlockMat&) {};
    void apply (int nr, int nc, const double *u, int ldu, double *v, int ldv);
};

class BJacobi : public GlobalPrecon
{
private:
    const BlockMat *Ap;
    LocalMat **diag;  // inverse or factors of diagonal blocks

public:
    BJacobi();
   ~BJacobi();
 
    void setup(const BlockMat& A);
    void apply (int, int, const double *, int, double *, int);
};

class BSOR_Base : public GlobalPrecon
{
protected:
    const BlockMat *Ap;
    LocalMat **diag;  // inverse or factors of diagonal blocks
    int    *idiag;

    double omega_;
    int iterations_;

public:
    BSOR_Base();
    virtual ~BSOR_Base();
 
    double& omega() {return omega_;}
    int& iterations() {return iterations_;}

    void setup(const BlockMat& A, double omega = 1.0, int iterations = 1);
};

class BSOR : public BSOR_Base
{
public:
    BSOR() {};
   ~BSOR() {};
    void apply (int, int, const double *, int, double *, int);
};

class BSSOR : public BSOR_Base
{
public:
    BSSOR() {};
   ~BSSOR() {};
    void apply (int, int, const double *, int, double *, int);
};

#endif // _BRELAX_H_
