//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#include "DenseMat.h"

DenseMat::DenseMat(const DenseMat& A)
{
    nrow = A.nrow; 
    ncol = A.ncol;
    register double *p = a = new double[nrow*ncol];
    register double *q = A.a;
    for (int i=0; i<nrow*ncol; i++)
        *p++ = *q++;
}

void DenseMat::Print(ostream& os) const
{
        // check not an implicit inverse
        assert (a != NULL);

#ifdef _GNU_GCC_3
        ios::fmtflags olda = os.setf(ios::right,ios::adjustfield);
        ios::fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
#else
#ifdef _GNU_GCC_296
        long olda = os.setf(ios::right,ios::adjustfield);
        long oldf = os.setf(ios::scientific,ios::floatfield);
#else
        long olda = os.setf(ios::right,ios::adjustfield);
        long oldf = os.setf(ios::scientific,ios::floatfield);
#endif
#endif

        int oldp = os.precision(12);

        const double *p = a;
        for (int j=0; j<numcol(); j++)
        for (int i=0; i<numrow(); i++)
               os << i+1 << "  " << j+1 << "  " << *p++ << endl;

        os.setf(olda,ios::adjustfield);
        os.setf(oldf,ios::floatfield);
        os.precision(oldp);
}

// non-member functions

ostream& operator << (ostream& os, const DenseMat& mat)
{
        // should check not an implicit inverse
 
#ifdef _GNU_GCC_3
        ios::fmtflags olda = os.setf(ios::right,ios::adjustfield);
        ios::fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
#else
#ifdef _GNU_GCC_296
        long olda = os.setf(ios::right,ios::adjustfield);
        long oldf = os.setf(ios::scientific,ios::floatfield);
#else
        long olda = os.setf(ios::right,ios::adjustfield);
        long oldf = os.setf(ios::scientific,ios::floatfield);
#endif
#endif

        int oldp = os.precision(12);

        const double *a = &mat(0,0);
        for (int j=0; j<mat.numcol(); j++)
        for (int i=0; i<mat.numrow(); i++)
               os << i+1 << "  " << j+1 << "  " << *a++ << endl;

        os.setf(olda,ios::adjustfield);
        os.setf(oldf,ios::floatfield);
        os.precision(oldp);
        return os;
}
