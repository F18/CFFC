/* Math.h:  Header file for some useful math macros. */

#ifndef _MATH_MACROS_INCLUDED
#define _MATH_MACROS_INCLUDED

/* Include C++ math library. */

#include <cmath>
#include <iostream>
#include <cassert>

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

#define	FIFTEEN  15.00
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
#define ONETHIRD  0.3333333333333333
#define TWOTHIRDS 0.6666666666666667
#define FOURTHIRDS 1.333333333333333
#define TENTHIRDS  3.333333333333333 
#define ONESIXTH   0.1666666666666667

#define HUNDRED  100.00
#define THOUSAND 1000.00
#define MILLION  1000000.00
#define MILLI    0.001
#define MICRO    0.000001
#define TOLER    0.0000001
#define NANO     0.000000001
#define PICO     0.000000000001

#ifndef PI
#define	PI	3.14159265358979323846
#endif
#define	SQRT_PI	1.7724538509055160273
#define ONE_COMPLEX Complex(1.0,0.0)

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

inline int Factorial(int N){
  if (N==0) return 1;
  else return N*Factorial(N-1);
}

inline double ConvertDomainToMinusOneOne (double xmin, double xmax, double x){
  // convert the domain [xmin,xmax] to [-1:1]
  return (2*x-xmax-xmin)/(xmax-xmin);
}

inline double ConvertDomainToZeroOne (double xmin, double xmax, double x){
  // convert the domain [xmin,xmax] to [0:1]
  return (x-xmin)/(xmax-xmin);
}

inline double ConvertDomain (double DomainMin, double DomainMax, double NewDomainMin, double NewDomainMax, double x){
  // convert the domain [DomainMin,DomainMax] to [NewDomainMin,NewDomainMax]
  return (NewDomainMin*DomainMax - DomainMin*NewDomainMax + x*(NewDomainMax - NewDomainMin))/(DomainMax - DomainMin);
}

inline int Pascals_Triangle(int n, int k) {
  assert(n>=k);
  return Factorial(n)/(Factorial(k)*Factorial(n-k));
}

/*!
 * Calculate the centroid and the area of a polygon defined by its vertices.
 * The centroid determined by this function is for the case when 
 * the polygon is treated as a sheet of uniform density.
 * If the polygon is a quadrilateral, this centroid type has the property 
 * that the value of a linear function at this location is equal to the 
 * average value of the linear function over the quadrilateral domain.
 * Note that there are 2 more other possibilities of defining the 
 * centroid of a polygon and the resultant centroids are all 
 * different for a quadrilateral.
 * The first alternative case is to consider the point masses at 
 * the vertices of the polygon whereas the second one is to represent
 * the sides of the polygon as wire rods of uniform density.
 * To find centers of gravity of uniform density sheets, 
 * one can simply divide the polygon into non-overlapping triangles and 
 * treat the system as a set of point masses at the centroids of these 
 * triangles with a mass equal to the area of the triangle.
 * The current implementation is based on the one suggested in chapter I.1 
 * (Centroid of a Polygon) of "Graphics Gems IV" by Paul S. Heckbert.
 * It gives correct results for both convex and concave polygons.
 * If edge crossing occurs the calculated area might still be positive,
 * so this is not a very good way of checking edge crossing.
 *
 * \param [in]  Vertices the array of vertices in x-y plane.
 *              The PositionVectorType must contain variables x and y (e.g. Vector2D class)
 * \param [in]   n       the number of entries in the vertices array.
 * \param [out] Centroid the value of the centroid is written here.
 * \param [out] area     the value of the polygon area is written here.
 *                       The algebraic sign of the area is positive for counterclockwise
 *                       ordering of vertices in x-y plane; otherwise negative.
 *
 * \return 0 for normal execution;
 *         1 if the polygon is degenerated (i.e. number of vertices less than 3);
 *         2 if area = zero and the centroid is undefined.
 */
template<typename PositionVectorType>
inline int polyCentroid(const PositionVectorType * Vertices, const int &n,
			PositionVectorType &Centroid, double &area){

  int i,j;
  
  double ai, atmp(0), xtmp(0), ytmp(0);
  
  if (n < 3) return 1;
  for (i=n-1, j=0; j < n; i=j, ++j){
    ai = Vertices[i].x * Vertices[j].y - Vertices[j].x * Vertices[i].y;
    atmp += ai;
    xtmp += (Vertices[j].x + Vertices[i].x) * ai;
    ytmp += (Vertices[j].y + Vertices[i].y) * ai;
  }
  area = 0.5* atmp;
  
  if (atmp != 0){
    Centroid.x = xtmp / (3*atmp);
    Centroid.y = ytmp / (3*atmp);
    return 0;
  }
  
  return 2;
}

/*!
 * Calculate the centroid and the area of a quadrilateral polygon.
 * This subroutine gives correct results for both concave 
 * and convex quadrilaterals.
 * However, the subroutine doesn't check for successful
 * execution (i.e. occurrence of negative area)!
 *
 * X0, X1, X2, X3 are the vertices of the quad in counterclockwise order.
 * If the vertices are in clockwise order the area is going to be negative!!!
 *
 * \param [out] Centroid the value of the centroid is written here.
 * \param [out] area     the value of the polygon area is written here.
 *                       The algebraic sign of the area is positive for counterclockwise
 *                       ordering of vertices in x-y plane; otherwise negative.
 * \return 0 for normal execution;
 *         2 if area = zero (i.e. the centroid is undefined).
 */
template<typename PositionVectorType>
inline int quadCentroid(const PositionVectorType &X0, const PositionVectorType &X1,
			const PositionVectorType &X2, const PositionVectorType &X3,
			PositionVectorType &Centroid, double &area){

  // Local variables
  double ai, atmp(0), xtmp(0), ytmp(0);

  // First edge between X3 and X0
  ai = X3.x * X0.y - X0.x * X3.y;
  atmp += ai;
  xtmp += (X0.x + X3.x) * ai;
  ytmp += (X0.y + X3.y) * ai;

  // Second edge between X0 and X1
  ai = X0.x * X1.y - X1.x * X0.y;
  atmp += ai;
  xtmp += (X1.x + X0.x) * ai;
  ytmp += (X1.y + X0.y) * ai;

  // Third edge between X1  and X2
  ai = X1.x * X2.y - X2.x * X1.y;
  atmp += ai;
  xtmp += (X2.x + X1.x) * ai;
  ytmp += (X2.y + X1.y) * ai;

  // Fourth edge between X2 and X3
  ai = X2.x * X3.y - X3.x * X2.y;
  atmp += ai;
  xtmp += (X3.x + X2.x) * ai;
  ytmp += (X3.y + X2.y) * ai;

  // Calculate the centroid and area of cell (ii,jj)
  area = 0.5* atmp;

  if (atmp != 0){
    Centroid.x = xtmp / (3*atmp);
    Centroid.y = ytmp / (3*atmp);
    return 0;
  }

  return 2;
}

/*
 * \verbatim
 * Parameterization of two lines:                                    
 *                                                                   
 *   Ra = Xa1 + s*Xa2  and  Rb = Xb1 + t*Xb2                         
 *                                                                   
 * Set Ra = Rb:                                                      
 *                                                                   
 *   [ Xa2.x  - Xb2.x ] [ s ] = [ Xb1.x - Xa1.x ]                    
 *   [ Xa2.y  - Xb2.y ] [ t ] = [ Xb1.y - Xa1.y ]                    
 *                                                                   
 * By Cramer's Rule:                                                 
 *                                                                   
 *   s = ((Xb1.y - Xa1.y)*Xb2.x - (Xb1.x - Xa1.x)*Xb2.y)/det         
 *   t = ((Xb1.y - Xa1.y)*Xa2.x - (Xb1.x - Xa1.x)*Xa2.y)/det         
 *                                                                   
 * where det = Xa2.y*Xb2.x - Xa2.x*Xb2.y                             
 *                                                                   
 * If det = 0 then the lines are parallel (maybe coincident).        
 *                                                                   
 * Then the point of intersection is given by                        
 *                                                                   
 *   x = Xa1.x + s*Xa2.x = Xb1.x + t*Xb2.x                           
 *   y = Xa1.y + s*Xa2.y = Xb1.y + t*Xb2.y                           
 *                                                                   
 * Note that this point of intersection may not be contained within  
 * the two line segments.                                            
 *
 * \endverbatim 
 * The PositionVectorType must have defined x,y variables and scalar multiplication operator (e.g. Vector2D class)
 *
 * \param Xa1 the start point of the first line
 * \param Xa3 the end point of the first line
 * \param Xb1 the start point of the second line
 * \param Xb3 the end point of the second line
 *
 * \return 0 if the point couldn't be determined (i.e. parallel or coincident lines)
 * \return 1 for successful execution
 */
template<typename PositionVectorType>
inline int getIntersectionPointOfTwoLines(const PositionVectorType &Xa1,
					  const PositionVectorType &Xa3,
					  const PositionVectorType &Xb1,
					  const PositionVectorType &Xb3,
					  PositionVectorType &Xp) {

  double det, s, t;
  PositionVectorType Xa2(Xa3 - Xa1);
  PositionVectorType Xb2(Xb3 - Xb1);

  // Check for too small numbers
  if (abs(Xa2) < NANO || abs(Xb2) < NANO) return 0;

  // Calculate the determinant
  det = Xa2.y*Xb2.x - Xa2.x*Xb2.y;

  // Check for existence of solution
  if (fabs(det)/min(Xa2*Xa2,Xb2*Xb2) < NANO) return 0;

  // Determine the solution
  s = ((Xb1.y - Xa1.y)*Xb2.x - (Xb1.x - Xa1.x)*Xb2.y)/det;
  t = ((Xb1.y - Xa1.y)*Xa2.x - (Xb1.x - Xa1.x)*Xa2.y)/det;

  // Get the intersection point
  Xp = Xa1 + s*Xa2;

  // Return successful execution.
  return 1;
}

/*!
 * Determine the type of the quadrilateral
 * formed with the ordered vertices X1, X2, X3, X4
 * by analysing the relative position of the diagonal
 * intersection with respect to the diagonals.
 * \verbose  
 *                                            (4)
 *                   (4)     * (2)             *
 *     *(4)           *     /|               /  |
 *    / \             |\   / |              /   |
 *(1)*   \            | \ /  |             / (2)|
 *    \   *(3)        |  /   |            / *   |
 *     \  |           | / \  |           //  \  |
 *      \ |           *    \ |        (1)*    \ |
 *       *(2)        (1)    \|                  *
 *                           * (3)               (3)
 *      (a)             (b)                 (c)   
 *
 * \endverbose
 * Possible types are: \n
 *  - degenerated (i.e. impossible to be formed by the 4 vertexes because the intersection point doesn't exist) \n
 *  - convex (a) (i.e. the intersection point is contained by both diagonals) \n
 *  - concave (c) (i.e. the intersection point is contained by only one of the diagonals) \n
 *  - crossed (b) (i.e. the intersection point is not contained by any of the diagonals)  \n
 *
 * \return 0 for degenerated \n
 *         1 for convex \n
 *         2 for concave. The X2-X4 diagonal is contained by the quad.
 *         3 for concave. The X1-X3 diagonal is contained by the quad.
 *         4 for crossed.
 */
template<typename PositionVectorType>
inline int Find_Quadrilateral_Type(const PositionVectorType &X1, const PositionVectorType &X2,
				   const PositionVectorType &X3, const PositionVectorType &X4){

  // Determine the intersection point of the two diagonals
  // For details see getIntersectionPointOfTwoLines()
  
  double det, s, t;
  PositionVectorType Xa(X3 - X1);
  PositionVectorType Xb(X4 - X2);
  bool PointContainedByDiagonal1, PointContainedByDiagonal2;

  // Check for too small numbers
  if (abs(Xa) < NANO || abs(Xb) < NANO) return 0;

  // Calculate the determinant
  det = Xa.y*Xb.x - Xa.x*Xb.y;

  // Check for existence of solution
  if (fabs(det)/min(Xa*Xa,Xb*Xb) < NANO) return 0;

  // Determine the solution
  s = ((X2.y - X1.y)*Xb.x - (X2.x - X1.x)*Xb.y)/det;
  t = ((X2.y - X1.y)*Xa.x - (X2.x - X1.x)*Xa.y)/det;

  // Analyse the intersection parameters.
  PointContainedByDiagonal1 = (s > ZERO-NANO && s < ONE+NANO);
  PointContainedByDiagonal2 = (t > ZERO-NANO && t < ONE+NANO);

  if ( PointContainedByDiagonal1 && PointContainedByDiagonal2 ){
    // Interior intersection
    return 1;

  } else if (PointContainedByDiagonal1){
    // Second diagonal is contained by the quadrilateral
    return 2;

  } else if (PointContainedByDiagonal2){
    // First diagonal is contained by the quadrilateral
    return 3;

  } else {
    // No diagonal contains the intersection
    return 4;
  }  

  // Something went wrong!
  return 0;
}

#endif // _MATH_MACROS_INCLUDED
