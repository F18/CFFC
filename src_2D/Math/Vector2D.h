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

/* Include math macro and linked list header files. */

#ifndef _MATH_MACROS_INCLUDED
#include "Math.h"
#endif // _MATH_MACROS_INCLUDED

/* Define the class. */

#define	NUM_COMP_VECTOR2D    2

/*!
 * Class: Vector2D
 *
 * @brief Defintion of a two-dimensional vector class.
 *
 * \verbatim
 * Member functions
 *      zero  -- Assign zero vector.
 *      x     -- Return x-component.
 *      y     -- Return y-component.
 *      abs   -- Absolute value (magnitude) of vector.
 *      sqr   -- Square of vector.
 *      dot   -- Inner product of two vectors.
 *      cross -- Cross product of two vectors (scalar).
 *      vmin  -- Minimum of two vectors.
 *      vmax  -- Maximum of two vectors.
 *
 * Member operators
 *      V -- a 2D vector
 *      a -- a scalar (double)
 *
 * V = V;
 * a = V[i];
 * V = V + V;
 * V = V - V;
 * V = V * V; (inner product)
 * V = a * V;
 * V = V * a;
 * V = V / a;
 * V = V ^ V; (cross product, scalar in this case)
 * V = +V;
 * V = -V;
 * V += V;
 * V -= V;
 * V *= a;
 * V /= a;
 * V == V;
 * V != V;
 * cout << V; (output function)
 * cin  >> V; (input function)
 * \endverbatim
 */
class Vector2D{
  private:
  public:
    double x,y;  //!< x- and y-components of 2D vector.
	         // Made public so can access them.
		
    //@{ @name Constructors.

    //! Creation constructor.
    Vector2D(void) {
       x = ZERO; y = ZERO;
    }
  
    //! Copy constructor.
    Vector2D(const double &V) {
       x = V; y = V;
    }
  
    //! Copy constructor.
    Vector2D(const Vector2D &V) {
       x = V.x; y = V.y;
    }

    //! Copy constructor.
    Vector2D(const double &xx,
	     const double &yy) {
       x = xx; y = yy;
    }

    /* Destructor. */
    // ~Vector2D(void);
    // Use automatically generated destructor.

    //@} 

    //! Zero (re-initialize) the vector.
    void zero(void);

    //@{ @name Absolute value (magnitude) of vector.
    double abs(void);
    double abs(void) const;
    friend double abs(const Vector2D &V);
    //@}

    //@{ @name Square of vector.
    double sqr(void);
    double sqr(void) const;
    friend double sqr(const Vector2D &V);
    //@}

    //@{ @name Inner (dot) product of two vectors.
    double dot(const Vector2D &V);
    double dot(const Vector2D &V) const;
    friend double dot(const Vector2D &V1,
		      const Vector2D &V2);
    //@}

    //@{ @name Cross (outer) product of two vectors.
    double cross(const Vector2D &V);
    double cross(const Vector2D &V) const;
    friend double cross(const Vector2D &V1,
		    	const Vector2D &V2);
    //@}

    //@{ @name Maximum of two vectors.
    Vector2D vmax(const Vector2D &V);
    Vector2D vmax(const Vector2D &V) const;
    friend Vector2D vmax(const Vector2D &V1,
			 const Vector2D &V2);
    //@}

    //@{ @name Minimum of two vectors.
    Vector2D vmin(const Vector2D &V);
    Vector2D vmin(const Vector2D &V) const;
    friend Vector2D vmin(const Vector2D &V1,
 	    		 const Vector2D &V2);

    //@}

    //@{ @name Set the vector to the given vector in the local rotated frame.
    void Rotate(const Vector2D &V, const Vector2D &norm_dir) {
      x =   V.x*norm_dir.x + V.y*norm_dir.y;
      y = - V.x*norm_dir.y + V.y*norm_dir.x;
    }

    /* Assignment operator. */
    // Vector2D operator = (const Vector2D &V);
    // Use automatically generated assignment operator.

    //@{ @name Index operator.
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
    //@}

    //@{ @name Binary arithmetic operators.
    friend Vector2D operator +(const Vector2D &V1, const Vector2D &V2);
    friend Vector2D operator -(const Vector2D &V1, const Vector2D &V2);
    friend double operator *(const Vector2D &V1, const Vector2D &V2);
    friend Vector2D operator *(const Vector2D &V, const double &a);
    friend Vector2D operator *(const double &a, const Vector2D &V);
    friend Vector2D operator /(const Vector2D &V, const double &a);
    friend double operator ^(const Vector2D &V1, const Vector2D &V2);
    //@}

    //@{ Unary arithmetic operators.
    friend Vector2D operator +(const Vector2D &V);
    friend Vector2D operator -(const Vector2D &V);
    //@}

    //@{ Shortcut arithmetic operators.
    friend Vector2D &operator +=(Vector2D &V1, const Vector2D &V2);
    friend Vector2D &operator -=(Vector2D &V1, const Vector2D &V2);
    Vector2D &operator *=(const double &a);
    Vector2D &operator /=(const double &a);
    //@}

    //@{ Relational operators.
    friend int operator ==(const Vector2D &V1, const Vector2D &V2);
    friend int operator !=(const Vector2D &V1, const Vector2D &V2);
    //@}
    
    //@{ @name Input-output operators.
    friend ostream &operator << (ostream &out_file, const Vector2D &V);
    friend istream &operator >> (istream &in_file, Vector2D &V);
    //@}

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

inline Vector2D& Vector2D::operator *=(const double &a) {
  x *= a; y *= a;
  return *this;
}

inline Vector2D& Vector2D::operator /=(const double &a) {
  x /= a; y /= a;
  return *this;
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

/********************************************************
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
const Vector2D ihat(ONE,ZERO);
const Vector2D jhat(ZERO,ONE);

/********************************************************
 * Vector2D -- External subroutines.                    *
 ********************************************************/

extern Vector2D nhat(const Vector2D X1, const Vector2D X2);

extern Vector2D Interpolate(const Vector2D X1, const double f1,
			    const Vector2D X2, const double f2);

extern Vector2D Interpolate(const Vector2D X1, const double f1,
			    const Vector2D X2, const double f2,
			    const double epsilon);

extern int Line_Intersection(const Vector2D Xa1,
			     const Vector2D Xa3,
			     const Vector2D Xb1,
			     const Vector2D Xb3,
			     Vector2D &Xp,
			     const double eps);

extern int Line_Intersection(const Vector2D Xa1,
			     const Vector2D Xa3,
			     const Vector2D Xb1,
			     const Vector2D Xb3,
			     Vector2D &Xp);

extern int Line_Intersection(const Vector2D Xa1,
			     const Vector2D Xa3,
			     const Vector2D Xb1,
			     const Vector2D Xb3);

extern int Point_On_Line(const Vector2D X1,
			 const Vector2D X2,
			 const Vector2D Xp);

extern double Triangle_Area(const Vector2D X1,
			    const Vector2D X2,
			    const Vector2D X3);

extern double Quadrilateral_Area(const Vector2D X1,
				 const Vector2D X2,
				 const Vector2D X3,
				 const Vector2D X4);

#endif /* _VECTOR2D_INCLUDED  */
