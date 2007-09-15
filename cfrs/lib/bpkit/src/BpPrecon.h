//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#ifndef _BPPRECON_H_
#define _BPPRECON_H_

class BpPrecon
{
public:
    virtual ~BpPrecon() {}

    virtual void apply  (int, int, const double *, int, double *, int) {}
    virtual void applyt (int, int, const double *, int, double *, int) {}
    virtual void applyr (int, int, const double *, int, double *, int) {}
    virtual void applyrt(int, int, const double *, int, double *, int) {}
    virtual void applyl (int, int, const double *, int, double *, int) {}
    virtual void applylt(int, int, const double *, int, double *, int) {}
    virtual void applyc (int, int, const double *, int, double *, int) {}
    virtual void applyct(int, int, const double *, int, double *, int) {}
};

#endif // _BPPRECON_H_
