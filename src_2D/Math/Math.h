/* Math.h:  Header file for some useful math macros. */

#ifndef _MATH_MACROS_INCLUDED
#define _MATH_MACROS_INCLUDED

/* Include C++ math library. */

#include <cmath>
#include <iostream>

using namespace std;

/* Define some useful constants. */

#define	ZERO	0.00
#define	ONE	1.00
#define	TWO	2.00
#define	THREE	3.00
#define	FOUR	4.00
#define	FIVE	5.00
#define	SIX	6.00
#define	SEVEN	7.00
#define	EIGHT	8.00
#define	NINE	9.00
#define	TEN    10.00

#define TWENTY   20.00
#define THIRTY   30.00
#define FOURTY   40.00
#define FIFTY    50.00
#define SIXTY    60.00
#define SEVENTY  70.00
#define EIGHTY   80.00
#define NINETY   90.00

#define QUARTER 0.25
#define	HALF	0.50

#define HUNDRED  100.00
#define THOUSAND 1000.00
#define MILLION  1000000.00
#define MILLI    0.001
#define MICRO    0.000001
#define TOLER    0.0000001
#define NANO     0.000000001

#ifndef PI
#define	PI	3.14159265358979323846
#endif
#define	SQRT_PI	1.7724538509055160273

/* Define additional functions. */

// max(x,y)
inline int    max(int x, int y)
                { return (x>y) ? x:y; }
inline short  max(short x, short y)
                { return (x>y) ? x:y; }
inline long   max(long x, long y)
                { return (x>y) ? x:y; }
inline float  max(float x, float y)
                { return (x>y) ? x:y; }
inline double max(const double &x, const double &y)
                { return (x>y) ? x:y; }

// min(x,y)
inline int    min(int x, int y)
                { return (x<y) ? x:y; }
inline short  min(short x, short y)
                { return (x<y) ? x:y; }
inline long   min(long x, long y)
                { return (x<y) ? x:y; }
inline float  min(float x, float y)
                { return (x<y) ? x:y; }
inline double min(const double &x, const double &y)
                { return (x<y) ? x:y; }

// x^2 (sqr(x) is defined in /usr/include/math.h)
# ifndef _HP_CC
inline int    sqr(int x)
                { return x*x; }
inline double sqr(const double &x)
                { return x*x; }
# endif
inline short  sqr(short x)
                { return x*x; }
inline long   sqr(long x)
                { return x*x; }
inline float  sqr(float x)
                { return x*x; }

// x^3
inline int    cube(int x)
                { return x*x*x; }
inline short  cube(short x)
                { return x*x*x; }
inline long   cube(long x)
                { return x*x*x; }
inline float  cube(float x)
                { return x*x*x; }
inline double cube(const double &x)
                { return x*x*x; }

// sgn(x) (note that sgn(0)=1)
inline int sgn(int x)
             { return (x<0) ? -1:1; }
inline int sgn(short x)
             { return (x<0) ? -1:1; }
inline int sgn(long x)
             { return (x<0) ? -1:1; }
inline int sgn(float x)
             {return (x<0.) ? -1:1; }
inline int sgn(const double &x)
             {return (x<0.) ? -1:1; }

// arctan(y, x) 
inline float  arctan(float y, float x) { 
   float z; z = atan2(y, x);
   return (z>=ZERO) ? z:TWO*PI+z;
}
inline double arctan(const double &y, const double &x) { 
   double z; z = atan2(y, x);
   return (z>=ZERO) ? z:TWO*PI+z;
}

// heaviside(x)
inline float heaviside(float &x) {
    return (x > ZERO) ? ONE : ZERO;
}

inline double heaviside(const double &x) {
    return (x > ZERO) ? ONE : ZERO;
}

// For use with qsort (used for sorting list of integers)
// e.g., qsort(list, n_list_size, sizeof(int), compare_integers);
inline int compare_integers(const void *p, const void *q) {
    return *(int *)p - *(int *)q;
}

// Returns the slope of the line of best fit (in the 
// least-squares sense) where the data sets are:
//   ( x, y ) = 
//   ( 0, values[position] ),
//   ( 1, values[position+1] ),
//   ...
//   (n-1, values[position+n-1] )
//
// "values" is a circular array so the actual indexing 
// is the remainder after division by n.
inline double linear_regression_slope(double *values, int n, int position) {
    double ssxx = ZERO, ssxy = ZERO;
    double ymean = ZERO, xmean = (n-ONE)/TWO;
    for (int x = 0; x < n; x++) { 
	ymean += values[x];
	ssxx += x * x;
	ssxy += x * values[ (position+x) % n ];
    }
    ymean *= ONE / n;
    ssxx -= n * xmean * xmean;
    ssxy -= n * xmean * ymean;
    return ssxy / ssxx;
}

// factorial(x)
inline int factorial(const int &x) {
  int factor = 1;
  if (x<0) {
    cerr << "\nError in factorial function: No factorial for negative integers! "<<endl;
    exit(1);
  } else if (x==0) {
    return 1;
  } else {
    for (int i=1; i<=x; i++) {
      factor *= i;
    }
    return factor;
  }
}

template <class T>
T LinearInterpolation(const double &x1, const double &x2, const double &x,
		      const T &T1, const T &T2) {
  return (T1 + (T2-T1)*(x-x1)/(x2-x1));  
}

template <class T>
T TwoPointFiniteDifference(const T &T1, const T &T2, const double &d2_d1) {
  return (T2-T1)/d2_d1;
}

#endif // _MATH_MACROS_INCLUDED
