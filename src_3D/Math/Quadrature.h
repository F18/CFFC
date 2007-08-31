/********************************************************************
 * Program:  Quadrature.h                                           *
 *                                                                  *
 * Description:  This program generates numerical quadrature values *
 *               and weights for several different types.           *
 *               See:                                               *   
 *                 Numerical Recipes in C, Second Edition (1992).   *
 *                 Available at http://www.nr.com/                  *
 *                                                                  *
 ********************************************************************/
#ifndef _QUADRATURE_INCLUDED
#define _QUADRATURE_INCLUDED

//--------------------------------------------------//
// SAMPLE USAGE
//--------------------------------------------------//
//--------------------------------------------------//
// END SAMPLE USAGE
//--------------------------------------------------//

// includes
#include <math.h>

// constants
#define EPS 3.0e-11


/********************************************************************
 * gauleg
 *                                                                  
 * Determines the locations in the interval x1<x<x2 and the weights 
 * of the Gauss-Lobatto integration scheme for n integration points.
 * The locations are reported in vector x and the weights in vector 
 * w.  The first vector term corresponds to x=x1.
 ********************************************************************/
inline void gauleg(double x1, double x2, double x[], double w[], int n)
{
        int m,j,i;
        double z1,z,xm,xl,pp,p3,p2,p1;

        m=(n+1)/2;
        xm=0.5*(x2+x1);
        xl=0.5*(x2-x1);
        for (i=1;i<=m;i++) {
                z=cos(3.141592654*(i-0.25)/(n+0.5));
                do {
                        p1=1.0;
                        p2=0.0;
                        for (j=1;j<=n;j++) {
                                p3=p2;
                                p2=p1;
                                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
                        }
                        pp=n*(z*p1-p2)/(z*z-1.0);
                        z1=z;
                        z=z1-p1/pp;
                } while (fabs(z-z1) > EPS);
                x[i-1]=xm-xl*z;
                x[n+1-i-1]=xm+xl*z; 
                w[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
                w[n+1-i-1]=w[i-1];
        }
}


#undef EPS



#endif // _QUADRATURE_INCLUDED
