/* Vector3D.h:  Header file defining 3D vector class. */

#ifndef _VECTOR3D_INCLUDED
#define _VECTOR3D_INCLUDED

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

#define	NUM_COMP_VECTOR3D    3

/********************************************************
 * Class: Vector3D                                      *
 *                                                      *
 * Member functions                                     *
 *      zero  -- Assign zero vector.                    *
 *      x     -- Return x-component.                    *
 *      y     -- Return y-component.                    *
 *      z     -- Return z-component.                    *
 *      abs   -- Absolute value (magnitude) of vector.  *
 *      sqr   -- Square of vector.                      *
 *      dot   -- Inner product of two vectors.          *
 *      cross -- Cross product of two vectors.          *
 *      vmin  -- Minimum of two vectors.                *
 *      vmax  -- Maximum of two vectors.                *
 *                                                      *
 * Member operators                                     *
 *      V -- a 3D vector                                *
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
 * V = V ^ V; (cross product)                           *
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
class Vector3D{
  private:
  public:
    double x,y,z;  // x-, y-, and z-components of 3D vector.
		   // Made public so can access them.
			
    /* Creation, copy, and assignment constructors. */
    Vector3D(void) {
       x = ZERO; y = ZERO; z = ZERO;
    }

    Vector3D(const Vector3D &V) {
       x = V.x; y = V.y; z = V.z;
    }

    Vector3D(const double &xx,
	     const double &yy,
	     const double &zz) {
       x = xx; y = yy; z = zz;
    }

    /* Destructor. */
    // ~Vector3D(void);
    // Use automatically generated destructor.

    /* Zero (re-initialize) the vector. */
    void zero(void);

    /* Absolute value (magnitude) of vector. */
    double abs(void);
    double abs(void) const;
    friend double abs(const Vector3D &V);

    /* Square of vector. */
    double sqr(void);
    double sqr(void) const;
    friend double sqr(const Vector3D &V);

    /* Inner (dot) product of two vectors. */
    double dot(const Vector3D &V);
    double dot(const Vector3D &V) const;
    friend double dot(const Vector3D &V1,
		      const Vector3D &V2);
    
    /* Cross (outer) product of two vectors. */
    Vector3D cross(const Vector3D &V);
    Vector3D cross(const Vector3D &V) const;
    friend Vector3D cross(const Vector3D &V1,
		    	  const Vector3D &V2);
    
    /* Maximum of two vectors. */
    Vector3D vmax(const Vector3D &V);
    Vector3D vmax(const Vector3D &V) const;
    friend Vector3D vmax(const Vector3D &V1,
 			 const Vector3D &V2);

    /* Minimum of two vectors. */
    Vector3D vmin(const Vector3D &V);
    Vector3D vmin(const Vector3D &V) const;
    friend Vector3D vmin(const Vector3D &V1,
 	    		 const Vector3D &V2);

    /* Assignment operator. */
    // Vector3D operator = (const Vector3D &V);
    // Use automatically generated assignment operator.

    /* Index operator. */
    double &operator[](int index) {
      assert( index >= 1 && index <= NUM_COMP_VECTOR3D );
      switch(index) {
      case 1 :
	return (x);
      case 2 :
	return (y);
      case 3 :
	return (z);
      default:
	return (x);
      };
    }
    
    const double &operator[](int index) const {
      assert( index >= 1 && index <= NUM_COMP_VECTOR3D );
      switch(index) {
      case 1 :
	return (x);
      case 2 :
	return (y);
      case 3 :
	return (z);
      default:
	return (x);
      };
    }

    /* Binary arithmetic operators. */
    friend Vector3D operator +(const Vector3D &V1, const Vector3D &V2);
    friend Vector3D operator -(const Vector3D &V1, const Vector3D &V2);
    friend double operator *(const Vector3D &V1, const Vector3D &V2);
    friend Vector3D operator *(const Vector3D &V, const double &a);
    friend Vector3D operator *(const double &a, const Vector3D &V);
    friend Vector3D operator /(const Vector3D &V, const double &a);
    friend Vector3D operator ^(const Vector3D &V1, const Vector3D &V2);

    /* Unary arithmetic operators. */
    friend Vector3D operator +(const Vector3D &V);
    friend Vector3D operator -(const Vector3D &V);

    /* Shortcut arithmetic operators. */
    friend Vector3D &operator +=(Vector3D &V1, const Vector3D &V2);
    friend Vector3D &operator -=(Vector3D &V1, const Vector3D &V2);
    
    /* Relational operators. */
    friend int operator ==(const Vector3D &V1, const Vector3D &V2);
    friend int operator !=(const Vector3D &V1, const Vector3D &V2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const Vector3D &V);
    friend istream &operator >> (istream &in_file, Vector3D &V);

};

/********************************************************
 * Vector3D::zero -- Assign zero vector.                *
 ********************************************************/
inline void Vector3D::zero(void) {
    x = ZERO; y = ZERO; z = ZERO;
}

/********************************************************
 * Vector3D::abs -- Absolute value of a vector.         *
 ********************************************************/
inline double Vector3D::abs(void) {
    return (sqrt(x*x+y*y+z*z));
}

inline double Vector3D::abs(void) const {
    return (sqrt(x*x+y*y+z*z));
}

inline double abs(const Vector3D &V) {
    return (sqrt(V.x*V.x+V.y*V.y+V.z*V.z));  
}

/********************************************************
 * Vector3D::sqr -- Square of a vector.                 *
 ********************************************************/
inline double Vector3D::sqr(void) {
    return (x*x+y*y+z*z);
}

inline double Vector3D::sqr(void) const {
    return (x*x+y*y+z*z);
}

inline double sqr(const Vector3D &V) {
    return (V.x*V.x+V.y*V.y+V.z*V.z);  
}

/********************************************************
 * Vector3D::dot -- Inner product of two vectors.       *
 ********************************************************/
inline double Vector3D::dot(const Vector3D &V) {
    return (x*V.x+y*V.y+z*V.z);
}

inline double Vector3D::dot(const Vector3D &V) const {
    return (x*V.x+y*V.y+z*V.z);
}

inline double dot(const Vector3D &V1,
		  const Vector3D &V2) {
    return (V1.x*V2.x+V1.y*V2.y+V1.z*V2.z); 
}

/********************************************************
 * Vector3D::cross -- Cross product of two vectors.     *
 ********************************************************/
inline Vector3D Vector3D::cross(const Vector3D &V) {
    return (Vector3D(y*V.z-z*V.y,
		     z*V.x-x*V.z,
		     x*V.y-y*V.x));
}

inline Vector3D Vector3D::cross(const Vector3D &V) const {
    return (Vector3D(y*V.z-z*V.y,
		     z*V.x-x*V.z,
		     x*V.y-y*V.x));
}

inline Vector3D cross(const Vector3D &V1,
		      const Vector3D &V2) {
    return (Vector3D(V1.y*V2.z-V1.z*V2.y,
	             V1.z*V2.x-V1.x*V2.z,
	             V1.x*V2.y-V1.y*V2.x));
}

/********************************************************
 * Vector3D::max -- Maximum of two vectors.             *
 ********************************************************/
inline Vector3D Vector3D::vmax(const Vector3D &V) {
    return (Vector3D(max(x,V.x),
		     max(y,V.y),
		     max(z,V.z)));
}

inline Vector3D Vector3D::vmax(const Vector3D &V) const {
    return (Vector3D(max(x,V.x),
		     max(y,V.y),
		     max(z,V.z)));
}

inline Vector3D vmax(const Vector3D &V1,
		     const Vector3D &V2) {
    return (Vector3D(max(V1.x,V2.x),
		     max(V1.y,V2.y),
		     max(V1.z,V2.z)));
}

/********************************************************
 * Vector3D::min -- Minimum of two vectors.             *
 ********************************************************/
inline Vector3D Vector3D::vmin(const Vector3D &V) {
    return (Vector3D(min(x,V.x),
		     min(y,V.y),
		     min(z,V.z)));
}

inline Vector3D Vector3D::vmin(const Vector3D &V) const {
    return (Vector3D(min(x,V.x),
		     min(y,V.y),
		     min(z,V.z)));
}

inline Vector3D vmin(const Vector3D &V1,
		     const Vector3D &V2) {
    return (Vector3D(min(V1.x,V2.x),
		     min(V1.y,V2.y),
		     min(V1.z,V2.z)));
}

/********************************************************
 * Vector3D -- Binary arithmetic operators.             *
 ********************************************************/
inline Vector3D operator +(const Vector3D &V1, const Vector3D &V2) {
  return (Vector3D(V1.x+V2.x,V1.y+V2.y,V1.z+V2.z));
}

inline Vector3D operator -(const Vector3D &V1, const Vector3D &V2) {
  return (Vector3D(V1.x-V2.x,V1.y-V2.y,V1.z-V2.z));
}

// Inner product operator.
inline double operator *(const Vector3D &V1, const Vector3D &V2) {
   return (V1.x*V2.x+V1.y*V2.y+V1.z*V2.z);
}

inline Vector3D operator *(const Vector3D &V, const double &a) {
  return (Vector3D(a*V.x,a*V.y,a*V.z));
}

inline Vector3D operator *(const double &a, const Vector3D &V) {
  return (Vector3D(a*V.x,a*V.y,a*V.z));
}

inline Vector3D operator /(const Vector3D &V, const double &a) {
  return (Vector3D(V.x/a,V.y/a,V.z/a));
}

// Cross product operator.
inline Vector3D operator ^(const Vector3D &V1, const Vector3D &V2) {
  return (Vector3D(V1.y*V2.z-V1.z*V2.y,
	           V1.z*V2.x-V1.x*V2.z,
	           V1.x*V2.y-V1.y*V2.x));
}

/********************************************************
 * Vector3D -- Unary arithmetic operators.              *
 ********************************************************/
inline Vector3D operator +(const Vector3D &V) {
  return (Vector3D(V.x,V.y,V.z));
}

inline Vector3D operator -(const Vector3D &V) {
  return (Vector3D(-V.x,-V.y,-V.z));
}

/********************************************************
 * Vector3D -- Shortcut arithmetic operators.           *
 ********************************************************/
inline Vector3D &operator +=(Vector3D &V1, const Vector3D &V2) {
  V1.x += V2.x;
  V1.y += V2.y;
  V1.z += V2.z;
  return (V1);
}

inline Vector3D &operator -=(Vector3D &V1, const Vector3D &V2) {
  V1.x -= V2.x;
  V1.y -= V2.y;
  V1.z -= V2.z;
  return (V1);
}

/********************************************************
 * Vector3D -- Relational operators.                    *
 ********************************************************/
inline int operator ==(const Vector3D &V1, const Vector3D &V2) {
  return (V1.x == V2.x && V1.y == V2.y && V1.z == V2.z);
}

inline int operator !=(const Vector3D &V1, const Vector3D &V2) {
  return (V1.x != V2.x || V1.y != V2.y || V1.z != V2.z);
}

/********************************************************
 * Vector3D -- Input-output operators.                  *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const Vector3D &V) {
  out_file.setf(ios::scientific);
  out_file << " " << V.x << " " << V.y << " " << V.z;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, Vector3D &V) {
  in_file.setf(ios::skipws);
  in_file >> V.x >> V.y >> V.z;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/********************************************************
 * Vector3D -- Useful 3D vector constants.              *
 ********************************************************/
const Vector3D Vector3D_ZERO(ZERO,ZERO,ZERO);
const Vector3D Vector3D_NX(ONE,ZERO,ZERO);
const Vector3D Vector3D_NY(ZERO,ONE,ZERO);
const Vector3D Vector3D_NZ(ZERO,ZERO,ONE);

#endif /* _VECTOR3D_INCLUDED  */
