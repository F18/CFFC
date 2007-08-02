//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#ifndef _ARCH_H_
#define _ARCH_H_

#define name2(a,b) a ## b

#ifdef CRAY
#define F77NAME(x) x
#else
#define F77NAME(x) name2(x,_)
#endif

typedef int      integer;
typedef long int integer8;    /* should be same size as a pointer */
typedef double   doublereal;

/* Note that C and Fortran integer sizes should match, since the storage
 * for Fortran arrays (for Harwell-Boeing matrices) is allocated by C. 
 */

#endif // _ARCH_H_
