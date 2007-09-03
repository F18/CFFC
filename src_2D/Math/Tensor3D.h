/* Tensor3D.h:  Header file defining 3D Tensor class. */

#ifndef _TENSOR3D_INCLUDED
#define _TENSOR3D_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

/* Include math macro and 3D vector header file. */

#ifndef _MATH_MACROS_INCLUDED
#include "Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _VECTOR3D_INCLUDED
#include "Vector3D.h"
#endif //_VECTOR3D_INCLUDED

/* Define the class. */

#define	NUM_COMP_TENSOR3D    3

/********************************************************
 * Class: Tensor3D (A symmetric tensor)                 *
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
 *      T -- a 3D tensor                                *
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
 * V = T * T;                                           *
 * T = T / a;                                           *
 * T = +T;                                              *
 * T = -T;                                              *
 * T += T;                                              *
 * T -= T;                                              *
 * T *= a;                                              *
 * T /= a;                                              *
 * T == T;                                              *
 * T != T;                                              *
 * cout << T; (output function)                         *
 * cin  >> T; (input function)                          *
 *                                                      *
 ********************************************************/
class Tensor3D{
  private:
    static Vector3D temp_Vec;
  public:
    double xx,xy,xz,yy,yz,zz;  // Six entries of 3D tensor.
   	                       // Made public so can access them.
			      
    /* Creation, copy, and assignment constructors. */
    Tensor3D(void) {
       xx = ZERO; xy = ZERO; xz = ZERO; yy = ZERO; yz = ZERO; zz = ZERO; 
    }

    Tensor3D(const Tensor3D &T) {
       xx = T.xx; xy = T.xy; xz = T.xz; yy = T.yy; yz = T.yz; zz = T.zz;
    }

    Tensor3D(const double &Txx,
	     const double &Txy,
	     const double &Txz,
	     const double &Tyy,
	     const double &Tyz,
	     const double &Tzz) {
       xx = Txx; xy = Txy; xz = Txz, yy = Tyy; yz = Tyz; zz = Tzz;
    }

    /* Destructor. */
    // ~Tensor3D(void);
    // Use automatically generated destructor.

    /* Zero (re-initialize) the tensor. */
    void zero(void);
    
    /* Assign identity tensor. */
    void identity(void);

    /* Compute the trace of a tensor. */
    double trace(void);
    double trace(void) const;
    double trace(const Tensor3D &T);
    friend double trace(const Tensor3D &T);
    
    /* Compute the determinant of a tensor. */
    double det(void);
    double det(void) const;
    double det(const Tensor3D &T);
    friend double det(const Tensor3D &T);
      
    /* Compute the inverse of a tensor. */
    Tensor3D inverse(void);
    Tensor3D inverse(void) const;
    Tensor3D inverse(const Tensor3D &T);
    friend Tensor3D inverse(const Tensor3D &T);
      
    /* Assignment operator. */
    // Tensor3D operator = (const Tensor3D &T);
    // Use automatically generated assignment operator.

    /* Index operator. */
    Vector3D &operator[](int index) {
      assert( index >= 1 && index <= NUM_COMP_TENSOR3D );
      switch(index) {
      case 1 :
	temp_Vec = Vector3D(xx,xy,xz);
	return (temp_Vec);
      case 2 :
	temp_Vec = Vector3D(xy,yy,yz);
	return (temp_Vec);
      case 3 :
	temp_Vec = Vector3D(xz,yz,zz);
	return (temp_Vec);
      default:
	temp_Vec = Vector3D(xx,xy,xz);
	return (temp_Vec);
      };
    }
    
    const Vector3D &operator[](int index) const {
      assert( index >= 1 && index <= NUM_COMP_TENSOR3D );
      switch(index) {
      case 1 :
	temp_Vec = Vector3D(xx,xy,xz);
	return (temp_Vec);
      case 2 :
	temp_Vec = Vector3D(xy,yy,yz);
	return (temp_Vec);
      case 3 :
	temp_Vec = Vector3D(xz,yz,zz);
	return (temp_Vec);
      default:
	temp_Vec = Vector3D(xx,xy,xz);
	return (temp_Vec);
      };
    }

    /* Binary arithmetic operators. */
    friend Tensor3D operator +(const Tensor3D &T1, const Tensor3D &T2);
    friend Tensor3D operator -(const Tensor3D &T1, const Tensor3D &T2);
    friend Tensor3D operator *(const Tensor3D &T, const double &a);
    friend Tensor3D operator *(const double &a, const Tensor3D &T);
    friend Vector3D operator *(const Tensor3D &T, const Vector3D &V);
    friend Tensor3D operator *(const Tensor3D &T1, const Tensor3D &T2);
    friend Tensor3D operator /(const Tensor3D &T, const double &a);

    /* Unary arithmetic operators. */
    friend Tensor3D operator +(const Tensor3D &T);
    friend Tensor3D operator -(const Tensor3D &T);

    /* Shortcut arithmetic operators. */
    friend Tensor3D &operator +=(Tensor3D &T1, const Tensor3D &T2);
    friend Tensor3D &operator -=(Tensor3D &T1, const Tensor3D &T2);
    Tensor3D &operator *=(const double &a);
    Tensor3D &operator /=(const double &a);

    /* Relational operators. */
    friend int operator ==(const Tensor3D &T1, const Tensor3D &T2);
    friend int operator !=(const Tensor3D &T1, const Tensor3D &T2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const Tensor3D &T);
    friend istream &operator >> (istream &in_file, Tensor3D &T);

};

/********************************************************
 * Tensor3D::zero -- Assign zero tensor.                *
 ********************************************************/
inline void Tensor3D::zero(void) {
    xx = ZERO; xy = ZERO; xz = ZERO; yy = ZERO; yz = ZERO; zz = ZERO;
}

/********************************************************
 * Tensor3D::identity -- Assign identity tensor.        *
 ********************************************************/
inline void Tensor3D::identity(void) {
    xx = ONE; xy = ZERO; xz = ZERO; yy = ONE; yz = ZERO; zz = ONE;
}

/********************************************************
 * Tensor3D::trace -- Trace of a tensor.                *
 ********************************************************/
inline double Tensor3D::trace(void) {
    return(xx+yy+zz);
}

inline double Tensor3D::trace(void) const {
    return(xx+yy+zz);
}

inline double Tensor3D::trace(const Tensor3D &T) {
    return(T.xx+T.yy+T.zz);
}

inline double trace(const Tensor3D &T) {
    return(T.xx+T.yy+T.zz);
}
    
/********************************************************
 * Tensor3D::det -- Determinant of a tensor.            *
 ********************************************************/
inline double Tensor3D::det(void) {
    return(xx*yy*zz-xx*yz*yz-xy*xy*zz-xz*xz*yy+TWO*xy*xz*yz);
}

inline double Tensor3D::det(void) const {
    return(xx*yy*zz-xx*yz*yz-xy*xy*zz-xz*xz*yy+TWO*xy*xz*yz);
}

inline double Tensor3D::det(const Tensor3D &T) {
    return(T.xx*T.yy*T.zz-T.xx*T.yz*T.yz-T.xy*T.xy*T.zz-
	   T.xz*T.xz*T.yy+TWO*T.xy*T.xz*T.yz);
}

inline double det(const Tensor3D &T) {
    return(T.xx*T.yy*T.zz-T.xx*T.yz*T.yz-T.xy*T.xy*T.zz-
	   T.xz*T.xz*T.yy+TWO*T.xy*T.xz*T.yz);
}
      
/********************************************************
 * Tensor3D::inverse -- Inverse of a tensor.            *
 ********************************************************/
inline Tensor3D Tensor3D::inverse(void) {
  double inv_det_T = ONE/det();
  return (Tensor3D((yy*zz-yz*yz)*inv_det_T,
		   -(xy*zz-xz*yz)*inv_det_T,
		   -(xy*yz-xz*yy)*inv_det_T,
		   (xx*zz-xz*xz)*inv_det_T,
		   -(xx*yz-xy*xz)*inv_det_T,
		   (xx*yy-xy*xy)*inv_det_T));
}

inline Tensor3D Tensor3D::inverse(void) const {
  double inv_det_T = ONE/det();
  return (Tensor3D((yy*zz-yz*yz)*inv_det_T,
		   -(xy*zz-xz*yz)*inv_det_T,
		   -(xy*yz-xz*yy)*inv_det_T,
		   (xx*zz-xz*xz)*inv_det_T,
		   -(xx*yz-xy*xz)*inv_det_T,
		   (xx*yy-xy*xy)*inv_det_T));
}

inline Tensor3D Tensor3D::inverse(const Tensor3D &T) {
  double inv_det_T = ONE/T.det();
  return (Tensor3D((T.yy*T.zz-T.yz*T.yz)*inv_det_T,
		   -(T.xy*T.zz-T.xz*T.yz)*inv_det_T,
		   -(T.xy*T.yz-T.xz*T.yy)*inv_det_T,
		   (T.xx*T.zz-T.xz*T.xz)*inv_det_T,
		   -(T.xx*T.yz-T.xy*T.xz)*inv_det_T,
		   (T.xx*T.yy-T.xy*T.xy)*inv_det_T));
}

inline Tensor3D inverse(const Tensor3D &T) {
  double inv_det_T = ONE/T.det();
  return (Tensor3D((T.yy*T.zz-T.yz*T.yz)*inv_det_T,
		   -(T.xy*T.zz-T.xz*T.yz)*inv_det_T,
		   -(T.xy*T.yz-T.xz*T.yy)*inv_det_T,
		   (T.xx*T.zz-T.xz*T.xz)*inv_det_T,
		   -(T.xx*T.yz-T.xy*T.xz)*inv_det_T,
		   (T.xx*T.yy-T.xy*T.xy)*inv_det_T));
}

/********************************************************
 * Tensor3D -- Binary arithmetic operators.             *
 ********************************************************/
inline Tensor3D operator +(const Tensor3D &T1, const Tensor3D &T2) {
  return (Tensor3D(T1.xx+T2.xx,T1.xy+T2.xy,T1.xz+T2.xz,
		   T1.yy+T2.yy,T1.yz+T2.yz,T1.zz+T2.zz));
}

inline Tensor3D operator -(const Tensor3D &T1, const Tensor3D &T2) {
  return (Tensor3D(T1.xx-T2.xx,T1.xy-T2.xy,T1.xz-T2.xz,
		   T1.yy-T2.yy,T1.yz-T2.yz,T1.zz-T2.zz));
}

inline Tensor3D operator *(const Tensor3D &T, const double &a) {
  return (Tensor3D(a*T.xx,a*T.xy,a*T.xz,a*T.yy,a*T.yz,a*T.zz));
}

inline Tensor3D operator *(const double &a, const Tensor3D &T) {
  return (Tensor3D(a*T.xx,a*T.xy,a*T.xz,a*T.yy,a*T.yz,a*T.zz));
}

inline Vector3D operator *(const Tensor3D &T, const Vector3D &V) {
  return (Vector3D(T.xx*V.x+T.xy*V.y+T.xz*V.z,
                   T.xy*V.x+T.yy*V.y+T.yz*V.z,
                   T.xz*V.x+T.yz*V.y+T.zz*V.z));
}

inline Tensor3D operator *(const Tensor3D &T1, const Tensor3D &T2) {
   return (Tensor3D(T1.xx*T2.xx+T1.xy*T2.xy+T1.xz*T2.xz, 
                    T1.xx*T2.xy+T1.xy*T2.yy+T1.xz*T2.yz,
                    T1.xx*T2.xz+T1.xy*T2.yz+T1.xz*T2.zz,
                    T1.xy*T2.xy+T1.yy*T2.yy+T1.yz*T2.yz,
                    T1.xy*T2.xz+T1.yy*T2.yz+T1.yz*T2.zz,
                    T1.xz*T2.xz+T1.yz*T2.yz+T1.zz*T2.zz ));
}

inline Tensor3D operator /(const Tensor3D &T, const double &a) {
  return (Tensor3D(T.xx/a,T.xy/a,T.xz/a,T.yy/a,T.yz/a,T.zz/a));
}

/********************************************************
 * Tensor3D -- Unary arithmetic operators.              *
 ********************************************************/
inline Tensor3D operator +(const Tensor3D &T) {
  return (Tensor3D(T.xx,T.xy,T.xz,T.yy,T.yz,T.zz));
}

inline Tensor3D operator -(const Tensor3D &T) {
  return (Tensor3D(-T.xx,-T.xy,-T.xz,-T.yy,-T.yz,-T.zz));
}

/********************************************************
 * Tensor3D -- Shortcut arithmetic operators.           *
 ********************************************************/
inline Tensor3D &operator +=(Tensor3D &T1, const Tensor3D &T2) {
  T1.xx += T2.xx;
  T1.xy += T2.xy;
  T1.xz += T2.xz;
  T1.yy += T2.yy;
  T1.yz += T2.yz;
  T1.zz += T2.zz;
  return (T1);
}

inline Tensor3D &operator -=(Tensor3D &T1, const Tensor3D &T2) {
  T1.xx -= T2.xx;
  T1.xy -= T2.xy;
  T1.xz -= T2.xz;
  T1.yy -= T2.yy;
  T1.yz -= T2.yz;
  T1.zz -= T2.zz;
  return (T1);
}

inline Tensor3D& Tensor3D::operator *=(const double &a) {
  xx *= a; xy *= a; xz *= a; yy *= a; yz *= a; zz *= a;
  return *this;
}

inline Tensor3D& Tensor3D::operator /=(const double &a) {
  xx /= a; xy /= a; xz /= a; yy /= a; yz /= a; zz /= a;
  return *this;
}

/********************************************************
 * Tensor3D -- Relational operators.                    *
 ********************************************************/
inline int operator ==(const Tensor3D &T1, const Tensor3D &T2) {
  return (T1.xx == T2.xx && T1.xy == T2.xy && T1.xz == T2.xz &&
	  T1.yy == T2.yy && T1.yz == T2.yz && T1.zz == T2.zz);
}

inline int operator !=(const Tensor3D &T1, const Tensor3D &T2) {
  return (T1.xx != T2.xx || T1.xy != T2.xy || T1.xz != T2.xz ||
	  T1.yy != T2.yy || T1.yz != T2.yz || T1.zz != T2.zz);
}

/********************************************************
 * Tensor3D -- Input-output operators.                  *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const Tensor3D &T) {
  out_file.setf(ios::scientific);
  out_file << " " << T.xx << " " << T.xy << " " << T.xz
	   << " " << T.yy << " " << T.yz << " " << T.zz;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, Tensor3D &T) {
  in_file.setf(ios::skipws);
  in_file >> T.xx >> T.xy >> T.xz >> T.yy >> T.yz >> T.zz;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/********************************************************
 * Tensor3D -- Useful 3D tensor constants.              *
 ********************************************************/
const Tensor3D Tensor3D_ZERO(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO);
const Tensor3D Tensor3D_IDENTITY(ONE,ZERO,ZERO,ONE,ZERO,ONE);

#endif /* _TENSOR3D_INCLUDED  */
