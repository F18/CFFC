/*      LAPACK++ (V. 1.0 Beta)						*/
/*      (C) 1992-1994 All Rights Reserved.				*/
/*
              LAPACK++ 1.0: Linear Algebra Package 1.0
               University of Tennessee, Knoxvilee, TN.
            Oak Ridge National Laboratory, Oak Ridge, TN.
        Authors: J. J. Dongarra, E. Greaser, R. Pozo, D. Walker
                 (C) 1992-1993 All Rights Reserved

                             NOTICE

 Permission to use, copy, modify, and distribute this software and
 its documentation for any purpose and without fee is hereby granted
 provided that the above copyright notice appear in all copies and
 that both the copyright notice and this permission notice appear in
 supporting documentation.

 Neither the Institutions (University of Tennessee, and Oak Ridge National
 Laboratory) nor the Authors make any representations about the suitability
 of this software for any purpose.  This software is provided ``as is''
 without express or implied warranty.

 LAPACK++ was funded in part by the U.S. Department of Energy, the
 National Science Foundation and the State of Tennessee.
*/

#ifndef _BLAS1_H_
#define _BLAS1_H_

#include "arch.h"

extern "C"
{


    double F77NAME(dasum)(const integer *n, const double *dx, const integer *incx);


    void F77NAME(daxpy)(const integer *n, const double *da, const double *dx, 
			const integer *incx, double *dy, const integer *incy);

    void F77NAME(dcopy)(const integer *n, double *dx, const integer *incx, double *dy, 
                        const integer *incy);


    double F77NAME(ddot)(const integer *n, const double *dx, const integer *incx, 
                        const double *dy, const integer *incy);

    double F77NAME(dnrm2)(const integer *n, const double *dx, const integer *incx); 

    void F77NAME(drot)(const integer *n, double *dx, const integer *incx, double *dy, 
                        const integer *incy, const double *c, const double *s);

    void F77NAME(drotg)(double *da, double *db, double *c, double *s);

    void F77NAME(dscal)(const integer *n, double *da, double *dx, const integer *incx);

    void F77NAME(dswap)(const integer *n, double *dx, const integer *incx, double *dy, 
                        const integer *incy);

    integer F77NAME(idamax)(const integer *n, const double *dx, const integer *incx);


/*           AT&T CC        Gnu g++  */
#if defined(COMPLEXH)|| defined(_Complex_h )

    double F77NAME(zdotc)(complex *c, const integer *n, const complex *cx, 
			const integer *incx, const complex *cy, const integer *incy);

    double F77NAME(zdotu)(complex *c, const integer *n, const complex *cx, 
		const integer *incx, const complex *cy, const integer *incy);

    void F77NAME(zaxpy)(const integer *n, const complex *da, const complex *dx, 
			const integer *incx, complex *dy, const integer *incy);

    void F77NAME(zcopy)(const integer *n, complex *dx, const integer *incx, 
				complex *dy, const integer *incy);

    double	F77NAME(dzasum)(const integer *n, const complex *dx, const integer *incx);

    double	F77NAME(dznrm2)(const integer *n, const complex *dx, const integer *incx); 

    void F77NAME(zdscal)(const integer *n, const double *da, complex *dx, 
			const integer *incx);

    void F77NAME(zscal)(const integer *n, const complex *da, complex *dx, 
			const integer *incx);

    integer F77NAME(izamax)(const integer *n, const complex *dx, const integer *incx);

    void F77NAME(zswap)(const integer *n, complex *dx, const integer *incx, 
				complex *dy, integer *incy);

#endif
}

#endif

