/* Tensor2D.h:  Header file defining 2D Tensor class. */

#ifndef _TENSOR2D_INCLUDED
#define _TENSOR2D_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

/* Include math macro, 2D vector, and 3D vector header files. */

#ifndef _MATH_MACROS_INCLUDED
#include "Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _VECTOR2D_INCLUDED
#include "Vector2D.h"
#endif //_VECTOR2D_INCLUDED

#ifndef _VECTOR3D_INCLUDED
#include "Vector3D.h"
#endif //_VECTOR3D_INCLUDED

/* Define the class. */

#define	NUM_COMP_TENSOR2D    3

/********************************************************
 * Class: Tensor2D (A symmetric tensor)                 *
 *                                                      *
 * Member functions                                     *
 *      zero      -- Assign zero tensor.                *
 *      identity  -- Assign identity tensor.            *
 *      xx        -- Return xx-entry.                   *
 *      xy        -- Return xy-entry.                   *
 *      xz        -- Return xz-entry.                   *
 *      yy        -- Return xx-entry.                   *
 *      yz        -- Return yz-entry.                   *
 *      zz        -- Return xy-entry.                   *
 *      trace     -- Return trace of tensor.            *
 *      det       -- Return determinant of tensor.      *
 *      inverse   -- Teturn inverse of tensor.          *
 *                                                      *
 * Member operators                                     *
 *      T -- a 2D tensor                                *
 *      V -- a 3D vector                                *
 *      a -- a scalar (double)                          *
 *                                                      *
 * T = T;                                               *
 * V = T[i];                                            *
 * a = T[i][j];                                         * 
 * T = T + T;                                           *
 * T = T - T;                                           *
 * T = a * T;                                           *
 * T = T * a;                                           *
 * V = T * V;                                           *
 * T = T / a;                                           *
 * T = +T;                                              *
 * T = -T;                                              *
 * T += T;                                              *
 * T -= T;                                              *
 * T == T;                                              *
 * T != T;                                              *
 * cout << T; (output function)                         *
 * cin  >> T; (input function)                          *
 *                                                      *
 ********************************************************/
class Tensor2D{
  private:
    static Vector3D temp_Vec;
  public:
    double xx,xy,yy,zz;  // xx, xy, yy, zz-entries of 2D tensor.
    static double xz,yz; // Other 3D tensor entries.
	                 // Made public so can access them.
			      
    /* Creation, copy, and assignment constructors. */
    Tensor2D(void) {
       xx = ZERO; xy = ZERO; yy = ZERO; zz = ZERO; 
    }

    Tensor2D(const Tensor2D &T) {
       xx = T.xx; xy = T.xy; yy = T.yy; zz = T.zz;
    }

    Tensor2D(const double &Txx,
	     const double &Txy,
	     const double &Tyy,
	     const double &Tzz) {
       xx = Txx; xy = Txy; yy = Tyy; zz = Tzz;
    }

    /* Destructor. */
    // ~Tensor2D(void);
    // Use automatically generated destructor.

    /* Zero (re-initialize) the tensor. */
    void zero(void);
    
    /* Assign identity tensor. */
    void identity(void);

    /* Compute the trace of a tensor. */
    double trace(void);
    double trace(void) const;
    double trace(const Tensor2D &T);
    friend double trace(const Tensor2D &T);
    
    /* Compute the determinant of a tensor. */
    double det(void);
    double det(void) const;
    double det(const Tensor2D &T);
    friend double det(const Tensor2D &T);
      
    /* Compute the inverse of a tensor. */
    Tensor2D inverse(void);
    Tensor2D inverse(void) const;
    Tensor2D inverse(const Tensor2D &T);
    friend Tensor2D inverse(const Tensor2D &T);
      
    /* Assignment operator. */
    // Tensor2D operator = (const Tensor2D &T);
    // Use automatically generated assignment operator.

    /* Index operator. */
    Vector3D &operator[](int index) {
      assert( index >= 1 && index <= NUM_COMP_TENSOR2D );
      switch(index) {
      case 1 :
	temp_Vec = Vector3D(xx,xy,ZERO);
	return (temp_Vec);
      case 2 :
	temp_Vec = Vector3D(xy,yy,ZERO);
	return (temp_Vec);
      case 3 :
	temp_Vec = Vector3D(ZERO,ZERO,zz);
	return (temp_Vec);
      default:
	temp_Vec = Vector3D(xx,xy,ZERO);
	return (temp_Vec);
      };
    }
    
    const Vector3D &operator[](int index) const {
      assert( index >= 1 && index <= NUM_COMP_TENSOR2D );
      switch(index) {
      case 1 :
	temp_Vec = Vector3D(xx,xy,ZERO);
	return (temp_Vec);
      case 2 :
	temp_Vec = Vector3D(xy,yy,ZERO);
	return (temp_Vec);
      case 3 :
	temp_Vec = Vector3D(ZERO,ZERO,zz);
	return (temp_Vec);
      default:
	temp_Vec = Vector3D(xx,xy,ZERO);
	return (temp_Vec);
      };
    }

    /* Binary arithmetic operators. */
    friend Tensor2D operator +(const Tensor2D &T1, const Tensor2D &T2);
    friend Tensor2D operator -(const Tensor2D &T1, const Tensor2D &T2);
    friend Tensor2D operator *(const Tensor2D &T, const double &a);
    friend Tensor2D operator *(const double &a, const Tensor2D &T);
    friend Vector3D operator *(const Tensor2D &T, const Vector3D &V);
    friend Tensor2D operator /(const Tensor2D &T, const double &a);

    /* Unary arithmetic operators. */
    friend Tensor2D operator +(const Tensor2D &T);
    friend Tensor2D operator -(const Tensor2D &T);

    /* Shortcut arithmetic operators. */
    friend Tensor2D &operator +=(Tensor2D &T1, const Tensor2D &T2);
    friend Tensor2D &operator -=(Tensor2D &T1, const Tensor2D &T2);
    
    /* Relational operators. */
    friend int operator ==(const Tensor2D &T1, const Tensor2D &T2);
    friend int operator !=(const Tensor2D &T1, const Tensor2D &T2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const Tensor2D &T);
    friend istream &operator >> (istream &in_file, Tensor2D &T);

};

/********************************************************
 * Tensor2D::zero -- Assign zero tensor.                *
 ********************************************************/
inline void Tensor2D::zero(void) {
    xx = ZERO; xy = ZERO; yy = ZERO; zz = ZERO;
}

/********************************************************
 * Tensor2D::identity -- Assign identity tensor.        *
 ********************************************************/
inline void Tensor2D::identity(void) {
    xx = ONE; xy = ZERO; yy = ONE; zz = ONE;
}

/********************************************************
 * Tensor2D::trace -- Trace of a tensor.                *
 ********************************************************/
inline double Tensor2D::trace(void) {
    return(xx+yy+zz);
}

inline double Tensor2D::trace(void) const {
    return(xx+yy+zz);
}

inline double Tensor2D::trace(const Tensor2D &T) {
    return(T.xx+T.yy+T.zz);
}

inline double trace(const Tensor2D &T) {
    return(T.xx+T.yy+T.zz);
}
    
/********************************************************
 * Tensor2D::det -- Determinant of a tensor.            *
 ********************************************************/
inline double Tensor2D::det(void) {
    return(zz*(xx*yy-xy*xy));
}

inline double Tensor2D::det(void) const {
    return(zz*(xx*yy-xy*xy));
}

inline double Tensor2D::det(const Tensor2D &T) {
    return(T.zz*(T.xx*T.yy-T.xy*T.xy));
}

inline double det(const Tensor2D &T) {
    return(T.zz*(T.xx*T.yy-T.xy*T.xy));
}
      
/********************************************************
 * Tensor2D::inverse -- Inverse of a tensor.            *
 ********************************************************/
inline Tensor2D Tensor2D::inverse(void) {
  double inv_det_T = ONE/det();
  return (Tensor2D(yy*zz*inv_det_T,
		   -xy*zz*inv_det_T,
		   xx*zz*inv_det_T,
		   ONE/zz));
}

inline Tensor2D Tensor2D::inverse(void) const {
  double inv_det_T = ONE/det();
  return (Tensor2D(yy*zz*inv_det_T,
		   -xy*zz*inv_det_T,
		   xx*zz*inv_det_T,
		   ONE/zz));
}

inline Tensor2D Tensor2D::inverse(const Tensor2D &T) {
  double inv_det_T = ONE/T.det();
  return (Tensor2D(T.yy*T.zz*inv_det_T,
		   -T.xy*T.zz*inv_det_T,
		   T.xx*T.zz*inv_det_T,
		   ONE/T.zz));
}

inline Tensor2D inverse(const Tensor2D &T) {
  double inv_det_T = ONE/T.det();
  return (Tensor2D(T.yy*T.zz*inv_det_T,
		   -T.xy*T.zz*inv_det_T,
		   T.xx*T.zz*inv_det_T,
		   ONE/T.zz));
}

/********************************************************
 * Tensor2D -- Binary arithmetic operators.             *
 ********************************************************/
inline Tensor2D operator +(const Tensor2D &T1, const Tensor2D &T2) {
  return (Tensor2D(T1.xx+T2.xx,T1.xy+T2.xy,T1.yy+T2.yy,T1.zz+T2.zz));
}

inline Tensor2D operator -(const Tensor2D &T1, const Tensor2D &T2) {
  return (Tensor2D(T1.xx-T2.xx,T1.xy-T2.xy,T1.yy-T2.yy,T1.zz-T2.zz));
}

inline Tensor2D operator *(const Tensor2D &T, const double &a) {
  return (Tensor2D(a*T.xx,a*T.xy,a*T.yy,a*T.zz));
}

inline Tensor2D operator *(const double &a, const Tensor2D &T) {
  return (Tensor2D(a*T.xx,a*T.xy,a*T.yy,a*T.zz));
}

inline Vector3D operator *(const Tensor2D &T, const Vector3D &V) {
  return (Vector3D(T.xx*V.x+T.xy*V.y,
                   T.xy*V.x+T.yy*V.y,
                   T.zz*V.z));
}

inline Tensor2D operator /(const Tensor2D &T, const double &a) {
  return (Tensor2D(T.xx/a,T.xy/a,T.yy/a,T.zz/a));
}

/********************************************************
 * Tensor2D -- Unary arithmetic operators.              *
 ********************************************************/
inline Tensor2D operator +(const Tensor2D &T) {
  return (Tensor2D(T.xx,T.xy,T.yy,T.zz));
}

inline Tensor2D operator -(const Tensor2D &T) {
  return (Tensor2D(-T.xx,-T.xy,-T.yy,-T.zz));
}

/********************************************************
 * Tensor2D -- Shortcut arithmetic operators.           *
 ********************************************************/
inline Tensor2D &operator +=(Tensor2D &T1, const Tensor2D &T2) {
  T1.xx += T2.xx;
  T1.xy += T2.xy;
  T1.yy += T2.yy;
  T1.zz += T2.zz;
  return (T1);
}

inline Tensor2D &operator -=(Tensor2D &T1, const Tensor2D &T2) {
  T1.xx -= T2.xx;
  T1.xy -= T2.xy;
  T1.yy -= T2.yy;
  T1.zz -= T2.zz;
  return (T1);
}

/********************************************************
 * Tensor2D -- Relational operators.                    *
 ********************************************************/
inline int operator ==(const Tensor2D &T1, const Tensor2D &T2) {
  return (T1.xx == T2.xx && T1.xy == T2.xy &&
	  T1.yy == T2.yy && T1.zz == T2.zz);
}

inline int operator !=(const Tensor2D &T1, const Tensor2D &T2) {
  return (T1.xx != T2.xx || T1.xy != T2.xy ||
	  T1.yy != T2.yy || T1.zz != T2.zz);
}

/********************************************************
 * Tensor2D -- Input-output operators.                  *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const Tensor2D &T) {
  out_file.setf(ios::scientific);
  out_file << " " << T.xx << " " << T.xy << " " << T.yy << " " << T.zz;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, Tensor2D &T) {
  in_file.setf(ios::skipws);
  in_file >> T.xx >> T.xy >> T.yy >> T.zz;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/********************************************************
 * Tensor2D -- Useful 2D tensor constants.              *
 ********************************************************/
const Tensor2D Tensor2D_ZERO(ZERO,ZERO,ZERO,ZERO);
const Tensor2D Tensor2D_IDENTITY(ONE,ZERO,ONE,ONE);

#endif /* _TENSOR2D_INCLUDED  */
