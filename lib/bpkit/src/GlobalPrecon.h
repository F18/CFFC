//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#ifndef _GLOBALPRECON_H_
#define _GLOBALPRECON_H_

#include "BpPrecon.h"
#include "LocalPrecon.h"

class GlobalPrecon : public BpPrecon
{
protected:
    LocalPrecon local_precon;

public:
    GlobalPrecon() {local_precon.name = (LocalPreconName) 0;}
    virtual ~GlobalPrecon() {}

    // set the method for solving or inverting the blocks

    void localprecon(LocalPreconName b) 
        {local_precon.name = b;}
    void localprecon(LocalPreconName b, int i1) 
        {local_precon.name = b; local_precon.iarg1 = i1;}
    void localprecon(LocalPreconName b, double d1) 
        {local_precon.name = b; local_precon.darg1 = d1;}
    void localprecon(LocalPreconName b, int i1, double d1) 
      {local_precon.name = b; local_precon.iarg1 = i1; local_precon.darg1 = d1;}
    void localprecon(LocalPreconName b, double d1, int i1) 
      {local_precon.name = b; local_precon.iarg1 = i1; local_precon.darg1 = d1;}
    void localprecon(LocalPreconName b, int i1, int i2) 
      {local_precon.name = b; local_precon.iarg1 = i1; local_precon.iarg2 = i2;}
    void localprecon(LocalPreconName b, double d1, double d2) 
      {local_precon.name = b; local_precon.darg1 = d1; local_precon.darg2 = d2;}
};

#endif // _GLOBALPRECON_H_
