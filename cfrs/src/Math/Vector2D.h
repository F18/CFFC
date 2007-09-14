/* Vector2D.h:  Header file defining 2D vector class. */

#ifndef _VECTOR2D_INCLUDED
#define _VECTOR2D_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

/* Include math macro header file. */

#ifndef _MATH_MACROS_INCLUDED
#include "Math.h"
#endif // _MATH_MACROS_INCLUDED

/* Define the class. */

#define	NUM_COMP_VECTOR2D    2

/********************************************************
 * Class: Vector2D                                      *
 *                                                      *
 * Member functions                                     *
 *      zero  -- Assign zero vector.                    *
 *      x     -- Return x-component.                    *
 *      y     -- Return y-component.                    *
 *      abs   -- Absolute value (magnitude) of vector.  *
 *      sqr   -- Square of vector.                      *
 *      dot   -- Inner product of two vectors.          *
 *      cross -- Cross product of two vectors (scalar). *
 *      vmin  -- Minimum of two vectors.                *
 *      vmax  -- Maximum of two vectors.                *
 *                                                      *
 * Member operators                                     *
 *      V -- a 2D vector                                *
 *      a -- a scalar (double)                          *
 *                                                      *
 * V = V;                                               *
 * a = V[i];                                            *
 * V = V + V;                                           *
 * V = V - V;                                           *
 * V = V * V; (inner product)                           *
 * V = a * V;                                           *
 * V = V * a;                                           *
 * V = V / a;                                           *
 * V = V ^ V; (cross product, scalar in this case)      *
 * V = +V;                                              *
 * V = -V;                                              *
 * V += V;                                              *
 * V -= V;                                              *
 * V == V;                                              *
 * V != V;                                              *
 * cout << V; (output function)                         *
 * cin  >> V; (input function)                          *
 *                                                      *
 ********************************************************/
class Vector2D{
  private:
  public:
    double x,y;  // x- and y-components of 2D vector.
	         // Made public so can access them.
			
    /* Creation, copy, and assignment constructors. */
    Vector2D(void) {
       x = ZERO; y = ZERO;
    }

    Vector2D(const Vector2D &V) {
       x = V.x; y = V.y;
    }

    Vector2D(const double &xx,
	     const double &yy) {
       x = xx; y = yy;
    }

    /* Destructor. */
    // ~Vector2D(void);
    // Use automatically generated destructor.

    /* Zero (re-initialize) the vector. */
    void zero(void);

    /* Absolute value (magnitude) of vector. */
    double abs(void);
    double abs(void) const;
    friend double abs(const Vector2D &V);

    /* Square of vector. */
    double sqr(void);
    double sqr(void) const;
    friend double sqr(const Vector2D &V);

    /* Inner (dot) product of two vectors. */
    double dot(const Vector2D &V);
    double dot(const Vector2D &V) const;
    friend double dot(const Vector2D &V1,
		      const Vector2D &V2);
    
    /* Cross (outer) product of two vectors. */
    double cross(const Vector2D &V);
    double cross(const Vector2D &V) const;
    friend double cross(const Vector2D &V1,
		    	const Vector2D &V2);
    
    /* Maximum of two vectors. */
    Vector2D vmax(const Vector2D &V);
    Vector2D vmax(const Vector2D &V) const;
    friend Vector2D vmax(const Vector2D &V1,
			 const Vector2D &V2);

    /* Minimum of two vectors. */
    Vector2D vmin(const Vector2D &V);
    Vector2D vmin(const Vector2D &V) const;
    friend Vector2D vmin(const Vector2D &V1,
 	    		 const Vector2D &V2);

    /* Assignment operator. */
    Vector2D& operator= (const Vector2D &V){
      if (this == &V) return *this;
      x = V.x; y = V.y;
      return *this;
    }

    /* Index operator. */
    double &operator[](int index) {
      assert( index >= 1 && index <= NUM_COMP_VECTOR2D );
      switch(index) {
      case 1 :
	return (x);
      case 2 :
	return (y);
      default:
	return (x);
      };
    }
    
    const double &operator[](int index) const {
      assert( index >= 1 && index <= NUM_COMP_VECTOR2D );
      switch(index) {
      case 1 :
	return (x);
      case 2 :
	return (y);
      default:
	return (x);
      };
    }

    /* Binary arithmetic operators. */
    friend Vector2D operator +(const Vector2D &V1, const Vector2D &V2);
    friend Vector2D operator -(const Vector2D &V1, const Vector2D &V2);
    friend double operator *(const Vector2D &V1, const Vector2D &V2);
    friend Vector2D operator *(const Vector2D &V, const double &a);
    friend Vector2D operator *(const double &a, const Vector2D &V);
    friend Vector2D operator /(const Vector2D &V, const double &a);
    friend double operator ^(const Vector2D &V1, const Vector2D &V2);

    /* Unary arithmetic operators. */
    friend Vector2D operator +(const Vector2D &V);
    friend Vector2D operator -(const Vector2D &V);

    /* Shortcut arithmetic operators. */
    friend Vector2D &operator +=(Vector2D &V1, const Vector2D &V2);
    friend Vector2D &operator -=(Vector2D &V1, const Vector2D &V2);
    
    /* Relational operators. */
    friend int operator ==(const Vector2D &V1, const Vector2D &V2);
    friend int operator !=(const Vector2D &V1, const Vector2D &V2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const Vector2D &V);
    friend istream &operator >> (istream &in_file, Vector2D &V);

};

/********************************************************
 * Vector2D::zero -- Assign zero vector.                *
 ********************************************************/
inline void Vector2D::zero(void) {
    x = ZERO; y = ZERO;
}

/********************************************************
 * Vector2D::abs -- Absolute value of a vector.         *
 ********************************************************/
inline double Vector2D::abs(void) {
    return (sqrt(x*x+y*y));
}

inline double Vector2D::abs(void) const {
    return (sqrt(x*x+y*y));
}

inline double abs(const Vector2D &V) {
    return (sqrt(V.x*V.x+V.y*V.y));  
}

/********************************************************
 * Vector2D::sqr -- Square of a vector.                 *
 ********************************************************/
inline double Vector2D::sqr(void) {
    return (x*x+y*y);
}

inline double Vector2D::sqr(void) const {
    return (x*x+y*y);
}

inline double sqr(const Vector2D &V) {
    return (V.x*V.x+V.y*V.y);  
}

/********************************************************
 * Vector2D::dot -- Inner product of two vectors.       *
 ********************************************************/
inline double Vector2D::dot(const Vector2D &V) {
    return (x*V.x+y*V.y);
}

inline double Vector2D::dot(const Vector2D &V) const {
    return (x*V.x+y*V.y);
}

inline double dot(const Vector2D &V1,
		  const Vector2D &V2) {
    return (V1.x*V2.x+V1.y*V2.y); 
}

/********************************************************
 * Vector2D::cross -- Cross product of two vectors.     *
 ********************************************************/
inline double Vector2D::cross(const Vector2D &V) {
    return (x*V.y-y*V.x);
}

inline double Vector2D::cross(const Vector2D &V) const {
    return (x*V.y-y*V.x);
}

inline double cross(const Vector2D &V1,
		    const Vector2D &V2) {
    return (V1.x*V2.y-V1.y*V2.x);
}

/********************************************************
 * Vector2D::max -- Maximum of two vectors.             *
 ********************************************************/
inline Vector2D Vector2D::vmax(const Vector2D &V) {
    return (Vector2D(max(x,V.x),
		     max(y,V.y)));
}

inline Vector2D Vector2D::vmax(const Vector2D &V) const {
    return (Vector2D(max(x,V.x),
		     max(y,V.y)));
}

inline Vector2D vmax(const Vector2D &V1,
		     const Vector2D &V2) {
    return (Vector2D(max(V1.x,V2.x),
		     max(V1.y,V2.y)));
}

/********************************************************
 * Vector2D::min -- Minimum of two vectors.             *
 ********************************************************/
inline Vector2D Vector2D::vmin(const Vector2D &V) {
    return (Vector2D(min(x,V.x),
		     min(y,V.y)));
}

inline Vector2D Vector2D::vmin(const Vector2D &V) const {
    return (Vector2D(min(x,V.x),
		     min(y,V.y)));
}

inline Vector2D vmin(const Vector2D &V1,
		     const Vector2D &V2) {
    return (Vector2D(min(V1.x,V2.x),
		     min(V1.y,V2.y)));
}

/********************************************************
 * Vector2D -- Binary arithmetic operators.             *
 ********************************************************/
inline Vector2D operator +(const Vector2D &V1, const Vector2D &V2) {
  return (Vector2D(V1.x+V2.x,V1.y+V2.y));
}

inline Vector2D operator -(const Vector2D &V1, const Vector2D &V2) {
  return (Vector2D(V1.x-V2.x,V1.y-V2.y));
}

// Inner product operator.
inline double operator *(const Vector2D &V1, const Vector2D &V2) {
   return (V1.x*V2.x+V1.y*V2.y);
}

inline Vector2D operator *(const Vector2D &V, const double &a) {
  return (Vector2D(a*V.x,a*V.y));
}

inline Vector2D operator *(const double &a, const Vector2D &V) {
  return (Vector2D(a*V.x,a*V.y));
}

inline Vector2D operator /(const Vector2D &V, const double &a) {
  return (Vector2D(V.x/a,V.y/a));
}

// Cross product operator.
inline double operator ^(const Vector2D &V1, const Vector2D &V2) {
  return (V1.x*V2.y-V1.y*V2.x);
}

/********************************************************
 * Vector2D -- Unary arithmetic operators.              *
 ********************************************************/
inline Vector2D operator +(const Vector2D &V) {
  return (Vector2D(V.x,V.y));
}

inline Vector2D operator -(const Vector2D &V) {
  return (Vector2D(-V.x,-V.y));
}

/********************************************************
 * Vector2D -- Shortcut arithmetic operators.           *
 ********************************************************/
inline Vector2D &operator +=(Vector2D &V1, const Vector2D &V2) {
  V1.x += V2.x;
  V1.y += V2.y;
  return (V1);
}

inline Vector2D &operator -=(Vector2D &V1, const Vector2D &V2) {
  V1.x -= V2.x;
  V1.y -= V2.y;
  return (V1);
}

/********************************************************
 * Vector2D -- Relational operators.                    *
 ********************************************************/
inline int operator ==(const Vector2D &V1, const Vector2D &V2) {
  return (V1.x == V2.x && V1.y == V2.y);
}

inline int operator !=(const Vector2D &V1, const Vector2D &V2) {
  return (V1.x != V2.x || V1.y != V2.y);
}

/*******************************************************
 * Vector2D -- Input-output operators.                  *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const Vector2D &V) {
  out_file.setf(ios::scientific);
  out_file << " " << V.x << " " << V.y;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, Vector2D &V) {
  in_file.setf(ios::skipws);
  in_file >> V.x >> V.y;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/********************************************************
 * Vector2D -- Useful 2D vector constants.              *
 ********************************************************/
const Vector2D Vector2D_ZERO(ZERO,ZERO);
const Vector2D Vector2D_NX(ONE,ZERO);
const Vector2D Vector2D_NY(ZERO,ONE);

/**********************************************************************
 * Routine --  Interpolate.                                           *
 *                                                                    *
 * Used linear interpolation to detrmine the location on a line where *
 * f = 0.                                                             *
 *                                                                    *
 **********************************************************************/
inline Vector2D Interpolate(Vector2D p1, 
			    Vector2D p2, 
			    double   f1, 
			    double   f2) {

  Vector2D p = p2-p1, pa;
  double COS = dot(p,Vector2D(ONE,ZERO))/(abs(p));
  double SIN = sin(acos(COS));

  if ((fabs(f1) < TOLER*TOLER) && (fabs(f2) < TOLER*TOLER)) {
    p = p1;
  } else if (fabs(f1) < TOLER*TOLER) {
    p = p1;
  } else if (fabs(f2) < TOLER*TOLER) {
    p = p2;
  } else {
    pa.x =   COS*p.x + SIN*p.y;
    pa.y = - SIN*p.x + COS*p.y;
    p.x = p1.x - COS*f1*pa.x/(f2-f1) + SIN*f1*pa.y/(f2-f1);
    p.y = p1.y - SIN*f1*pa.x/(f2-f1) - COS*f1*pa.y/(f2-f1);
  }

  return p;

}

/**********************************************************************
 * Routine: Intersection_Point_Between_Two_Lines                      *
 *                                                                    *
 * Ra = Xa1 + s*Xa2  and  Rb = Xb1 + t*Xb2                            *
 *                                                                    *
 * set Ra = Rb                                                        *
 *                                                                    *
 *   [ Xa2.x  - Xb2.x ] [ s ] = [ Xb1.x - Xa1.x ]                     *
 *   [ Xa2.y  - Xb2.y ] [ t ] = [ Xb1.y - Xa1.y ]                     *
 *                                                                    *
 * By Cramer's Rule:                                                  *
 *                                                                    *
 *   s = ((Xb1.y - Xa1.y)*Xb2.x - (Xb1.x - Xa1.x)*Xb2.y)/det          *
 *   t = ((Xb1.y - Xa1.y)*Xa2.x - (Xb1.x - Xa1.x)*Xa2.y)/det          *
 *                                                                    *
 * where det = Xa2.y*Xb2.x - Xa2.x*Xb2.y                              *
 *                                                                    *
 * If det = 0 then the lines do not intersect.                        *
 *                                                                    *
 * Then the point of intersection is given by                         *
 *                                                                    *
 *   x = Xa1.x + s*Xa2.x = Xb1.x + t*Xb2.x                            *
 *   y = Xa1.y + s*Xa2.y = Xb1.y + t*Xb2.y                            *
 *                                                                    *
 **********************************************************************/
inline Vector2D Intersection_Point_Between_Two_Lines(const Vector2D Xa1,
						     const Vector2D Xa3,
						     const Vector2D Xb1,
						     const Vector2D Xb3) {

  //int intersect_flag = 0;
  double det, s; // t;
  Vector2D Xa2 = Xa3 - Xa1;
  Vector2D Xb2 = Xb3 - Xb1;
  Vector2D Xp  = Vector2D_ZERO;
  det = Xa2.y*Xb2.x - Xa2.x*Xb2.y;
  if (fabs(det) > TOLER) {
    s = ((Xb1.y - Xa1.y)*Xb2.x - (Xb1.x - Xa1.x)*Xb2.y)/det;
    //    t = ((Xb1.y - Xa1.y)*Xa2.x - (Xb1.x - Xa1.x)*Xa2.y)/det;
    Xp = Xa1 + s*Xa2;
    //intersect_flag = 1;
  }
  if (fabs(Xp.x) < TOLER) Xp.x = ZERO;
  if (fabs(Xp.y) < TOLER) Xp.y = ZERO;
  //return intersect_flag;
  return Xp;

}

#endif /* _VECTOR2D_INCLUDED  */
