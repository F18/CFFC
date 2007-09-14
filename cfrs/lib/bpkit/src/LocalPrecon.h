//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#ifndef _LOCALPRECON_H_
#define _LOCALPRECON_H_

// Make the values explicit, so we can match them with FORTRAN parameters.

enum LocalPreconName
{                          // arguments
    LP_LU           =  1,  //
    LP_INVERSE      =  2,  //
    LP_SVD          =  3,  // rthresh, athresh
    LP_RILUK        =  4,  // level, omega 
    LP_ILUT         =  5,  // lfil, thresh
    LP_APINV_TRUNC  =  6,  // semibw
    LP_APINV_BANDED =  7,  // semibw
    LP_APINV0       =  8,  //
    LP_APINVS       =  9,  // lfil, self_precon
    LP_DIAG         = 10,  //
    LP_TRIDIAG      = 11,  //
    LP_SOR          = 12,  // omega, iterations
    LP_SSOR         = 13,  // omega, iterations
    LP_GMRES        = 14,  // iterations, tol
    LP_DIAGDOM      = 15,  // 
    LP_GERSH        = 16   // alpha
};

class LocalPrecon
{
public:
    LocalPreconName name;
    int iarg1;
    int iarg2;
    double darg1;
    double darg2;
};

#endif // _LOCALPRECON_H_
