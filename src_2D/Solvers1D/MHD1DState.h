
/* MHD1DState.h:  Header file defining 1D MHD Solution State Classes. */

#ifndef _MHD1D_STATE_INCLUDED
#define _MHD1D_STATE_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

/* Include Math macro, CFD, 3D vector, and gas constant header files. */

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _VECTOR3D_INCLUDED
#include "../Math/Vector3D.h"
#endif //_VECTOR3D_INCLUDED

#ifndef _GAS_CONSTANTS_INCLUDED
#include "../Physics/GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED

/* Define the classes. */

#define	NUM_VAR_MHD1D    8

class MHD1D_cState;

/********************************************************
 * Class: MHD1D_pState                                  *
 *                                                      *
 * Member functions                                     *
 *      d       -- Return density.                      *
 *      v       -- Return flow velocity.                *
 *      B1      -- Return perturbative magnetic field.  *
 *      B0      -- Return intrinsic magnetic field.     *
 *      p       -- Return pressure.                     *
 *      g       -- Return specific heat ratio.          *
 *      gm1     -- Return g-1                           *
 *      gm1i    -- Return 1/(g-1).                      *
 *      setgas  -- Set gas constants.                   *
 *      B       -- Return total magnetic field.         *
 *      T       -- Return temperature.                  *
 *      e       -- Return specific internal energy.     *
 *      E       -- Return total energy.                 *
 *      E1      -- Return total perturbative energy.    *
 *      h       -- Return specific enthalpy.            *
 *      h1      -- Return spec. perturbative enthalpy.  *
 *      H       -- Return total enthalpy.               *
 *      H1      -- Return total perturbative enthalpy.  *
 *      a       -- Return sound speed.                  *
 *      a2      -- Return sound speed square.           *
 *      Va      -- Return Alfven wave velocity.         *
 *      Va2     -- Return Alfven wave speed square.     *
 *      s       -- Return specific entropy.             *
 *      dv      -- Return momentum.                     *
 *      U       -- Return conserved solution state.     *
 *      F       -- Return solution flux.                *
 *      lambda  -- Return eigenvalue.                   *
 *      rp      -- Return primitive right eigenvector.  *
 *      rc      -- Return conserved right eigenvector.  *
 *      lp      -- Return primitive left eigenvector.   *
 *      lc      -- Return conserved left eigenvector.   * 
 *                                                      *
 * Member operators                                     *
 *      W -- a primitive solution state                 *
 *      c -- a scalar (double)                          *
 *                                                      *
 * W = W;                                               *
 * c = W[i];                                            *
 * W = W + W;                                           *
 * W = W - W;                                           *
 * c = W * W; (inner product)                           *
 * W = c * W;                                           *
 * W = W * c;                                           *
 * W = W / c;                                           *
 * W = +W;                                              *
 * W = -W;                                              *
 * W += W;                                              *
 * W -= W;                                              *
 * W == W;                                              *
 * W != W;                                              *
 * cout << W; (output function)                         *
 * cin  >> W; (input function)                          *
 *                                                      *
 ********************************************************/
class MHD1D_pState{
  private:
  public:
    double          d;   // Density.
    Vector3D        v;   // Flow velocity (3D vector).
    Vector3D       B1;   // Perturbative magnetic field (3D vector).
    Vector3D       B0;   // Intrinsic magnetic field (3D vector).
    double          p;   // Pressure.
    static double   g;   // Specific heat ratio.
    static double gm1;   // g-1
    static double gm1i;  // 1/(g-1)
		         // Made public so can access them.
		      
    /* Creation, copy, and assignment constructors. */
    MHD1D_pState(void) {
       d = ONE; v.zero(); B1.zero(); B0.zero(); p = ONE;
    }

    MHD1D_pState(const MHD1D_pState &W) {
       d = W.d; v = W.v; B1 = W.B1; B0 = W.B0; p = W.p;
    }

    MHD1D_pState(const MHD1D_cState &U);

    MHD1D_pState(const double &rho,
	         const Vector3D &V,
	         const Vector3D &Bm1,
	         const double &pre) {
       d = rho; v = V; B1 = Bm1; B0.zero(); p = pre;
    }

    MHD1D_pState(const double &rho,
	         const Vector3D &V,
	         const Vector3D &Bm1,
                 const Vector3D &Bm0,
	         const double &pre) {
       d = rho; v = V; B1 = Bm1; B0 = Bm0; p = pre;
    }

    MHD1D_pState(const double &rho,
	         const double &vx,
	         const double &vy,
	         const double &vz,
	         const double &B1x,
	         const double &B1y,
	         const double &B1z,
	         const double &pre) {
       d = rho; v.x = vx; v.y = vy; v.z = vz;
       B1.x = B1x; B1.y = B1y; B1.z = B1z; B0.zero(); p = pre;
    }

    MHD1D_pState(const double &rho,
	         const double &vx,
	         const double &vy,
	         const double &vz,
	         const double &B1x,
	         const double &B1y,
	         const double &B1z,
	         const double &B0x,
	         const double &B0y,
	         const double &B0z,
	         const double &pre) {
       d = rho; v.x = vx; v.y = vy; v.z = vz;
       B1.x = B1x; B1.y = B1y; B1.z = B1z;
       B0.x = B0x; B0.y = B0y; B0.z = B0z; p = pre;
    }

    /* Destructor. */
    // ~MHD1D_pState(void);
    // Use automatically generated destructor.

    /* Set gas constants. */
    void setgas(void);
    void setgas(char *string_ptr);

    /* Total Magnetic Field. */
    Vector3D B(void);
    Vector3D B(void) const;

    /* Temperature. */
    double T(void);
    double T(void) const;

    /* Specific internal energy. */
    double e(void);
    double e(void) const;

    /* Total energy. */
    double E(void);
    double E(void) const;

    /* Total perturbative energy. */
    double E1(void);
    double E1(void) const;

    /* Specific enthalpy. */
    double h(void);
    double h(void) const;

    /* Specific perturbative enthalpy. */
    double h1(void);
    double h1(void) const;

    /* Total enthalpy. */
    double H(void);
    double H(void) const;

    /* Total perturbative enthalpy. */
    double H1(void);
    double H1(void) const;

    /* Sound speed. */
    double a(void);
    double a(void) const;

    /* Sound speed squared. */
    double a2(void);
    double a2(void) const;

    /* Alfven wave velocity. */
    Vector3D Va(void);
    Vector3D Va(void) const;

    /* Alfven wave speed squared. */
    double Va2(void);
    double Va2(void) const;

    /* Specific entropy. */
    double s(void);
    double s(void) const;

    /* Momentum. */
    Vector3D dv(void);
    Vector3D dv(void) const;
    double dv(const Vector3D &n);
    double dv(const Vector3D &n) const;

    /* Conserved solution state. */
    MHD1D_cState U(void);
    MHD1D_cState U(void) const;
    MHD1D_cState U(const MHD1D_pState &W);
    friend MHD1D_cState U(const MHD1D_pState &W);
    
    /* Solution flux. */
    MHD1D_cState F(void);
    MHD1D_cState F(void) const;
    MHD1D_cState F(const MHD1D_pState &W);
    friend MHD1D_cState F(const MHD1D_pState &W);

    /* Eigenvalue(s). */
    MHD1D_pState lambda(void);
    MHD1D_pState lambda(void) const;
    MHD1D_pState lambda(const MHD1D_pState &W);
    friend MHD1D_pState lambda(const MHD1D_pState &W);
    double lambda(int index);
    double lambda(int index) const;
    friend double lambda(const MHD1D_pState &W, int index);

    /* Primitive right eigenvector. */
    MHD1D_pState rp(int index);
    MHD1D_pState rp(int index) const;
    friend MHD1D_pState rp(const MHD1D_pState &W, int index);

    /* Conserved right eigenvector. */
    MHD1D_cState rc(int index);
    MHD1D_cState rc(int index) const;
    friend MHD1D_cState rc(const MHD1D_pState &W, int index);

    /* Primitive left eigenvector. */
    MHD1D_pState lp(int index);
    MHD1D_pState lp(int index) const;
    friend MHD1D_pState lp(const MHD1D_pState &W, int index);

    /* Conserved left eigenvector. */
    MHD1D_cState lc(int index);
    MHD1D_cState lc(int index) const;
    friend MHD1D_cState lc(const MHD1D_pState &W, int index);

    /* Assignment operator. */
    // MHD1D_pState operator = (const MHD1D_pState &W);
    // Use automatically generated assignment operator.

    /* Index operator. */
    double &operator[](int index) {
      assert( index >= 1 && index <= NUM_VAR_MHD1D );
      switch(index) {
        case 1 :
	  return (d);
        case 2 :
	  return (v.x);
        case 3 :
	  return (v.y);
        case 4 :
	  return (v.z);
        case 5 :
	  return (B1.x);
        case 6 :
	  return (B1.y);
        case 7 :
	  return (B1.z);
        case 8 :
	  return (p);
        default:
	  return (d);
      };
    }
    
    const double &operator[](int index) const {
      assert( index >= 1 && index <= NUM_VAR_MHD1D );
      switch(index) {
        case 1 :
	  return (d);
        case 2 :
	  return (v.x);
        case 3 :
	  return (v.y);
        case 4 :
	  return (v.z);
        case 5 :
	  return (B1.x);
        case 6 :
	  return (B1.y);
        case 7 :
	  return (B1.z);
        case 8 :
	  return (p);
        default:
	  return (d);
      };
    }

    /* Binary arithmetic operators. */
    friend MHD1D_pState operator +(const MHD1D_pState &W1, const MHD1D_pState &W2);
    friend MHD1D_pState operator -(const MHD1D_pState &W1, const MHD1D_pState &W2);
    friend double operator *(const MHD1D_pState &W1, const MHD1D_pState &W2);
    friend MHD1D_pState operator *(const MHD1D_pState &W, const double &a);
    friend MHD1D_pState operator *(const double &a, const MHD1D_pState &W);
    friend MHD1D_pState operator /(const MHD1D_pState &W, const double &a);

    /* Unary arithmetic operators. */
    friend MHD1D_pState operator +(const MHD1D_pState &W);
    friend MHD1D_pState operator -(const MHD1D_pState &W);

    /* Shortcut arithmetic operators. */
    friend MHD1D_pState &operator +=(MHD1D_pState &W1, const MHD1D_pState &W2);
    friend MHD1D_pState &operator -=(MHD1D_pState &W1, const MHD1D_pState &W2);
    
    /* Relational operators. */
    friend int operator ==(const MHD1D_pState &W1, const MHD1D_pState &W2);
    friend int operator !=(const MHD1D_pState &W1, const MHD1D_pState &W2);

    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const MHD1D_pState &W);
    friend istream &operator >> (istream &in_file,  MHD1D_pState &W);
    
};

/********************************************************
 * Class: MHD1D_cState                                  *
 *                                                      *
 * Member functions                                     *
 *      d       -- Return density.                      *
 *      dv      -- Return momentum.                     *
 *      B1      -- Return perturbative magentic field.  *
 *      B0      -- Return intrinsic magnetic field.     *
 *      E1      -- Return total perturbative energy.    *
 *      g       -- Return specific heat ratio.          *
 *      gm1     -- Return g-1.                          *
 *      gm1i    -- Return 1/(g-1).                      *
 *      setgas  -- Set gas constants.                   *
 *      v       -- Return flow velocity.                *
 *      p       -- Return pressure.                     *
 *      B       -- Return total magnetic field.         *
 *      T       -- Return temperature.                  *
 *      e       -- Return specific internal energy.     *
 *      E       -- Return total energy.                 *
 *      h       -- Return specific enthalpy.            *
 *      h1      -- Return perturbative enthalpy.        *
 *      H       -- Return total enthalpy.               *
 *      H1      -- Return total perturbative enthalpy.  *
 *      a       -- Return sound speed.                  *
 *      a2      -- Return sound speed square.           *
 *      Va      -- Return Alfven wave velocity.         *
 *      Va2     -- Return Alfven wave speed square.     *
 *      s       -- Return specific entropy.             *
 *      W       -- Return primitive solution state.     *
 *                                                      *
 * Member operators                                     *
 *      U -- a primitive solution state                 *
 *      c -- a scalar (double)                          *
 *                                                      *
 * U = U;                                               *
 * c = U[i];                                            *
 * U = U + U;                                           *
 * U = U - U;                                           *
 * c = U * U; (inner product)                           *
 * U = c * U;                                           *
 * U = U * c;                                           *
 * U = U / c;                                           *
 * U = +U;                                              *
 * U = -U;                                              *
 * U += U;                                              *
 * U -= U;                                              *
 * U == U;                                              *
 * U != U;                                              *
 * cout << U; (output function)                         *
 * cin  >> U; (input function)                          *
 *                                                      *
 ********************************************************/
class MHD1D_cState{
  private:
  public:
    double          d;   // Density.
    Vector3D       dv;   // Momentum.
    Vector3D       B1;   // Perturbative magnetic field (3D vector).
    Vector3D       B0;   // Intrinsic magnetic field (3D vector).
    double         E1;   // Total perturbative Energy.
    static double   g;   // Specific heat ratio.
    static double gm1;   // g-1
    static double gm1i;  // 1/(g-1)
	                 // Made public so can access them.
		      
    /* Creation, copy, and assignment constructors. */
    MHD1D_cState(void) {
       d = ONE; dv.zero(); B1.zero(); B0.zero(); E1 = ONE/(GAMMA_MONATOMIC-ONE);
    }

    MHD1D_cState(const MHD1D_cState &U) {
       d = U.d; dv = U.dv; B1 = U.B1; B0 = U.B0; E1 = U.E1;
    }

    MHD1D_cState(const MHD1D_pState &W);

    MHD1D_cState(const double &rho,
	         const Vector3D &rhoV,
	         const Vector3D &Bm1,
	         const double &Etotal) {
       d = rho; dv = rhoV; B1 = Bm1; B0.zero(); E1 = Etotal;
    }

    MHD1D_cState(const double &rho,
	         const Vector3D &rhoV,
	         const Vector3D &Bm1,
                 const Vector3D &Bm0,
	         const double &Etotal) {
       d = rho; dv = rhoV; B1 = Bm1; B0 = Bm0; E1 = Etotal;
    }

    MHD1D_cState(const double &rho,
	         const double &rhovx,
	         const double &rhovy,
	         const double &rhovz,
	         const double &B1x,
	         const double &B1y,
	         const double &B1z,
	         const double &Etotal) {
       d = rho; dv.x = rhovx; dv.y = rhovy; dv.z = rhovz;
       B1.x = B1x; B1.y = B1y; B1.z = B1z; B0.zero(); E1 = Etotal;
    }

    MHD1D_cState(const double &rho,
	         const double &rhovx,
	         const double &rhovy,
	         const double &rhovz,
	         const double &B1x,
	         const double &B1y,
	         const double &B1z,
	         const double &B0x,
	         const double &B0y,
	         const double &B0z,
	         const double &Etotal) {
       d = rho; dv.x = rhovx; dv.y = rhovy; dv.z = rhovz;
       B1.x = B1x; B1.y = B1y; B1.z = B1z;
       B0.x = B0x; B0.y = B0y; B0.z = B0z; E1 = Etotal;
    }
    
    /* Destructor. */
    // ~MHD1D_cState(void);
    // Use automatically generated destructor.

    /* Set gas constants. */
    void setgas(void);
    void setgas(char *string_ptr);

    /* Flow velocity. */
    Vector3D v(void);
    Vector3D v(void) const;
    double v(const Vector3D &n);
    double v(const Vector3D &n) const;
    
    /* Pressure. */
    double p(void);
    double p(void) const;

    /* Total Magnetic Field. */
    Vector3D B(void);
    Vector3D B(void) const;

    /* Temperature. */
    double T(void);
    double T(void) const;

    /* Specific internal energy. */
    double e(void);
    double e(void) const;

    /* Total energy. */
    double E(void);
    double E(void) const;

    /* Specific enthalpy. */
    double h(void);
    double h(void) const;

    /* Specific perturbative enthalpy. */
    double h1(void);
    double h1(void) const;

    /* Total enthalpy. */
    double H(void);
    double H(void) const;

    /* Total perturbative enthalpy. */
    double H1(void);
    double H1(void) const;

    /* Sound speed. */
    double a(void);
    double a(void) const;

    /* Sound speed squared. */
    double a2(void);
    double a2(void) const;

    /* Alfven wave velocity. */
    Vector3D Va(void);
    Vector3D Va(void) const;

    /* Alfven wave speed squared. */
    double Va2(void);
    double Va2(void) const;

    /* Specific entropy. */
    double s(void);
    double s(void) const;

    /* Primitive solution state. */
    MHD1D_pState W(void);
    MHD1D_pState W(void) const;
    MHD1D_pState W(const MHD1D_cState &U);
    friend MHD1D_pState W(const MHD1D_cState &U);
    
    /* Assignment operator. */
    // MHD1D_cState operator = (const MHD1D_cState &U);
    // Use automatically generated assignment operator.

    /* Index operator. */
    double &operator[](int index) {
      assert( index >= 1 && index <= NUM_VAR_MHD1D );
      switch(index) {
        case 1 :
	  return (d);
        case 2 :
	  return (dv.x);
        case 3 :
	  return (dv.y);
        case 4 :
	  return (dv.z);
        case 5 :
	  return (B1.x);
        case 6 :
	  return (B1.y);
        case 7 :
	  return (B1.z);
        case 8 :
	  return (E1);
        default:
	  return (d);
      };
    }
    
    const double &operator[](int index) const {
      assert( index >= 1 && index <= NUM_VAR_MHD1D );
      switch(index) {
        case 1 :
	  return (d);
        case 2 :
	  return (dv.x);
        case 3 :
	  return (dv.y);
        case 4 :
	  return (dv.z);
        case 5 :
	  return (B1.x);
        case 6 :
	  return (B1.y);
        case 7 :
	  return (B1.z);
        case 8 :
	  return (E1);
        default:
	  return (d);
      };
    }

    /* Binary arithmetic operators. */
    friend MHD1D_cState operator +(const MHD1D_cState &U1, const MHD1D_cState &U2);
    friend MHD1D_cState operator -(const MHD1D_cState &U1, const MHD1D_cState &U2);
    friend double operator *(const MHD1D_cState &U1, const MHD1D_cState &U2);
    friend MHD1D_cState operator *(const MHD1D_cState &U, const double &a);
    friend MHD1D_cState operator *(const double &a, const MHD1D_cState &U);
    friend MHD1D_cState operator /(const MHD1D_cState &U, const double &a);

    /* Unary arithmetic operators. */
    friend MHD1D_cState operator +(const MHD1D_cState &U);
    friend MHD1D_cState operator -(const MHD1D_cState &U);

    /* Shortcut arithmetic operators. */
    friend MHD1D_cState &operator +=(MHD1D_cState &U1, const MHD1D_cState &U2);
    friend MHD1D_cState &operator -=(MHD1D_cState &U1, const MHD1D_cState &U2);
    
    /* Relational operators. */
    friend int operator ==(const MHD1D_cState &U1, const MHD1D_cState &U2);
    friend int operator !=(const MHD1D_cState &U1, const MHD1D_cState &U2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const MHD1D_cState &U);
    friend istream &operator >> (istream &in_file,  MHD1D_cState &U);
    
};

/********************************************************
 * MHD1D_pState::setgas -- Assign gas constants.        *
 ********************************************************/
inline void MHD1D_pState::setgas(void) {
    g = GAMMA_MONATOMIC;
    gm1 = g - ONE;
    gm1i = ONE/gm1;
}

inline void MHD1D_pState::setgas(char *string_ptr) {
   if (strcmp(string_ptr, "MONATOMIC") == 0) {
     g = GAMMA_MONATOMIC;
   } else if (strcmp(string_ptr, "DIATOMIC") == 0) {
     g = GAMMA_DIATOMIC;
   } else if (strcmp(string_ptr, "POLYATOMIC") == 0) {
     g = GAMMA_POLYATOMIC;
   } else if (strcmp(string_ptr, "BRIOWU") == 0) {
     g = GAMMA_TWO;
   } else if (strcmp(string_ptr, "ISOTHERMAL") == 0) {
     g = GAMMA_ISOTHERMAL;
   } /* endif */
   gm1 = g - ONE;
   gm1i = ONE/gm1;
}

/********************************************************
 * MHD1D_pState::B -- Total Magnetic Field.             *
 ********************************************************/
inline Vector3D MHD1D_pState::B(void) {
    return (B1+B0);
}

inline Vector3D MHD1D_pState::B(void) const {
    return (B1+B0);
}

/********************************************************
 * MHD1D_pState::T -- Temperature.                      *
 ********************************************************/
inline double MHD1D_pState::T(void) {
    return (p/d);
}

inline double MHD1D_pState::T(void) const {
    return (p/d);
}

/********************************************************
 * MHD1D_pState::e -- Specific internal energy.         *
 ********************************************************/
inline double MHD1D_pState::e(void) {
    return (p/(gm1*d));
}

inline double MHD1D_pState::e(void) const {
    return (p/(gm1*d));
}

/********************************************************
 * MHD1D_pState::E -- Total energy.                     *
 ********************************************************/
inline double MHD1D_pState::E(void) {
    return (p*gm1i + HALF*d*v.sqr() + HALF*B().sqr());
}

inline double MHD1D_pState::E(void) const {
    return (p*gm1i + HALF*d*v.sqr() + HALF*B().sqr());
}

/********************************************************
 * MHD1D_pState::E1 -- Total perturbative energy.       *
 ********************************************************/
inline double MHD1D_pState::E1(void) {
    return (p*gm1i + HALF*d*v.sqr() + HALF*B1.sqr());
}

inline double MHD1D_pState::E1(void) const {
    return (p*gm1i + HALF*d*v.sqr() + HALF*B1.sqr());
}

/********************************************************
 * MHD1D_pState::h -- Specific enthalpy.                *
 ********************************************************/
inline double MHD1D_pState::h(void) {
    return (g*p/(gm1*d) + HALF*v.sqr() + B().sqr()/d);
}

inline double MHD1D_pState::h(void) const {
    return (g*p/(gm1*d) + HALF*v.sqr() + B().sqr()/d);
}

/********************************************************
 * MHD1D_pState::h1 -- Specific perturbative enthalpy.  *
 ********************************************************/
inline double MHD1D_pState::h1(void) {
    return (g*p/(gm1*d) + HALF*v.sqr() + B1.sqr()/d);
}

inline double MHD1D_pState::h1(void) const {
    return (g*p/(gm1*d) + HALF*v.sqr() + B1.sqr()/d);
}

/********************************************************
 * MHD1D_pState::H -- Total enthalpy.                   *
 ********************************************************/
inline double MHD1D_pState::H(void) {
    return (g*gm1i*p + HALF*d*v.sqr() + B().sqr());
}

inline double MHD1D_pState::H(void) const {
    return (g*gm1i*p + HALF*d*v.sqr() + B().sqr());
}

/********************************************************
 * MHD1D_pState::H1 -- Total perturbative enthalpy.     *
 ********************************************************/
inline double MHD1D_pState::H1(void) {
    return (g*gm1i*p + HALF*d*v.sqr() + B1.sqr());
}

inline double MHD1D_pState::H1(void) const {
    return (g*gm1i*p + HALF*d*v.sqr() + B1.sqr());
}

/********************************************************
 * MHD1D_pState::a -- Sound speed.                      *
 ********************************************************/
inline double MHD1D_pState::a(void) {
    return (sqrt(g*p/d));
}

inline double MHD1D_pState::a(void) const {
    return (sqrt(g*p/d));
}

/********************************************************
 * MHD1D_pState::a2 -- Sound speed squared.             *
 ********************************************************/
inline double MHD1D_pState::a2(void) {
    return (g*p/d);
}

inline double MHD1D_pState::a2(void) const {
    return (g*p/d);
}

/********************************************************
 * MHD1D_pState::Va -- Alfven wave velocity.            *
 ********************************************************/
inline Vector3D MHD1D_pState::Va(void) {
    return (B()/sqrt(d));
}

inline Vector3D MHD1D_pState::Va(void) const {
    return (B()/sqrt(d));
}

/********************************************************
 * MHD1D_pState::Va2 -- Alfven wave speed squared.      *
 ********************************************************/
inline double MHD1D_pState::Va2(void) {
    return ((B().sqr())/d);
}

inline double MHD1D_pState::Va2(void) const {
    return ((B().sqr())/d);
}

/********************************************************
 * MHD1D_pState::s -- Specific entropy.                 *
 ********************************************************/
inline double MHD1D_pState::s(void) {
    return (-gm1i*log(p/pow(d, g)));
}

inline double MHD1D_pState::s(void) const {
    return (-gm1i*log(p/pow(d, g)));
}

/********************************************************
 * MHD1D_pState::dv -- Momentum.                        *
 ********************************************************/
inline Vector3D MHD1D_pState::dv(void) {
    return (d*v);
}

inline Vector3D MHD1D_pState::dv(void) const {
    return (d*v);
}

inline double MHD1D_pState::dv(const Vector3D &n) {
    return (d*(v*n));
}

inline double MHD1D_pState::dv(const Vector3D &n) const {
    return (d*(v*n));
}

/********************************************************
 * MHD1D_pState -- Binary arithmetic operators.         *
 ********************************************************/
inline MHD1D_pState operator +(const MHD1D_pState &W1, const MHD1D_pState &W2) {
  return (MHD1D_pState(W1.d+W2.d,W1.v+W2.v,W1.B1+W2.B1,W1.p+W2.p));
}

inline MHD1D_pState operator -(const MHD1D_pState &W1, const MHD1D_pState &W2) {
  return (MHD1D_pState(W1.d-W2.d,W1.v-W2.v,W1.B1-W2.B1,W1.p-W2.p));
}

// Inner product operator.
inline double operator *(const MHD1D_pState &W1, const MHD1D_pState &W2) {
   return (W1.d*W2.d+W1.v*W2.v+W1.B1*W2.B1+W1.p*W2.p);
}

inline MHD1D_pState operator *(const MHD1D_pState &W, const double &a) {
  return (MHD1D_pState(a*W.d,a*W.v,a*W.B1,a*W.p));
}

inline MHD1D_pState operator *(const double &a, const MHD1D_pState &W) {
  return (MHD1D_pState(a*W.d,a*W.v,a*W.B1,a*W.p));
}

inline MHD1D_pState operator /(const MHD1D_pState &W, const double &a) {
  return (MHD1D_pState(W.d/a,W.v/a,W.B1/a,W.p/a));
}

/********************************************************
 * MHD1D_pState -- Unary arithmetic operators.          *
 ********************************************************/
inline MHD1D_pState operator +(const MHD1D_pState &W) {
  return (MHD1D_pState(W.d,W.v,W.B1,W.p));
}

inline MHD1D_pState operator -(const MHD1D_pState &W) {
  return (MHD1D_pState(-W.d,-W.v,-W.B1,-W.p));
}

/********************************************************
 * MHD1D_pState -- Shortcut arithmetic operators.       *
 ********************************************************/
inline MHD1D_pState &operator +=(MHD1D_pState &W1, const MHD1D_pState &W2) {
  W1.d += W2.d;
  W1.v += W2.v;
  W1.B1 += W2.B1;
  W1.p += W2.p;
  return (W1);
}

inline MHD1D_pState &operator -=(MHD1D_pState &W1, const MHD1D_pState &W2) {
  W1.d -= W2.d;
  W1.v -= W2.v;
  W1.B1 -= W2.B1;
  W1.p -= W2.p;
  return (W1);
}

/********************************************************
 * MHD1D_pState -- Relational operators.                *
 ********************************************************/
inline int operator ==(const MHD1D_pState &W1, const MHD1D_pState &W2) {
  return (W1.d == W2.d && W1.v == W2.v && W1.B() == W2.B() && W1.p == W2.p);
}

inline int operator !=(const MHD1D_pState &W1, const MHD1D_pState &W2) {
  return (W1.d != W2.d || W1.v != W2.v || W1.B() != W2.B() || W1.p != W2.p);
}

/********************************************************
 * MHD1D_pState -- Input-output operators.              *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const MHD1D_pState &W) {
  Vector3D Bt; Bt = W.B();
  out_file.setf(ios::scientific);
  out_file << " " << W.d  << " " << W.v.x << " " << W.v.y << " " << W.v.z
	   << " " << Bt.x << " " << Bt.y << " " << Bt.z
           << " " << W.B0.x << " " << W.B0.y << " " << W.B0.z << " " << W.p;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, MHD1D_pState &W) {
  Vector3D Bt;
  in_file.setf(ios::skipws);
  in_file >> W.d >> W.v.x >> W.v.y >> W.v.z >> Bt.x >> Bt.y >> Bt.z
	  >> W.B0.x >> W.B0.y >> W.B0.z >> W.p;
  in_file.unsetf(ios::skipws); W.B1 = Bt - W.B0;
  return (in_file);
}

/********************************************************
 * MHD1D_cState::setgas -- Assign gas constants.        *
 ********************************************************/
inline void MHD1D_cState::setgas(void) {
    g = GAMMA_MONATOMIC;
    gm1 = g - ONE;
    gm1i = ONE/gm1;
}

inline void MHD1D_cState::setgas(char *string_ptr) {
   if (strcmp(string_ptr, "MONATOMIC") == 0) {
     g = GAMMA_MONATOMIC;
   } else if (strcmp(string_ptr, "DIATOMIC") == 0) {
     g = GAMMA_DIATOMIC;
   } else if (strcmp(string_ptr, "POLYATOMIC") == 0) {
     g = GAMMA_POLYATOMIC;
   } else if (strcmp(string_ptr, "BRIOWU") == 0) {
     g = GAMMA_TWO;
   } else if (strcmp(string_ptr, "ISOTHERMAL") == 0) {
     g = GAMMA_ISOTHERMAL;
   } /* endif */
   gm1 = g - ONE;
   gm1i = ONE/gm1;
}

/********************************************************
 * MHD1D_cState::v -- Flow velocity.                    *
 ********************************************************/
inline Vector3D MHD1D_cState::v(void) {
    return (dv/d);
}

inline Vector3D MHD1D_cState::v(void) const {
    return (dv/d);
}

inline double MHD1D_cState::v(const Vector3D &n) {
    return ((dv*n)/d);
}

inline double MHD1D_cState::v(const Vector3D &n) const {
    return ((dv*n)/d);
}

/********************************************************
 * MHD1D_cState::p -- Pressure.                         *
 ********************************************************/
inline double MHD1D_cState::p(void) {
    return (gm1*(E1 - HALF*dv.sqr()/d - HALF*B1.sqr()));
}

inline double MHD1D_cState::p(void) const {
    return (gm1*(E1 - HALF*dv.sqr()/d - HALF*B1.sqr()));
}

/********************************************************
 * MHD1D_cState::B -- Total Magnetic Field.             *
 ********************************************************/
inline Vector3D MHD1D_cState::B(void) {
    return (B1+B0);
}

inline Vector3D MHD1D_cState::B(void) const {
    return (B1+B0);
}

/********************************************************
 * MHD1D_cState::T -- Temperature.                      *
 ********************************************************/
inline double MHD1D_cState::T(void) {
    return (gm1*(E1 - HALF*dv.sqr()/d - HALF*B1.sqr())/d);
}

inline double MHD1D_cState::T(void) const {
    return (gm1*(E1 - HALF*dv.sqr()/d - HALF*B1.sqr())/d);
}

/********************************************************
 * MHD1D_cState::e -- Specific internal energy.         *
 ********************************************************/
inline double MHD1D_cState::e(void) {
    return ((E1 - HALF*dv.sqr()/d - HALF*B1.sqr())/d);
}

inline double MHD1D_cState::e(void) const {
    return ((E1 - HALF*dv.sqr()/d - HALF*B1.sqr())/d);
}

/********************************************************
 * MHD1D_cState::E -- Total energy.                     *
 ********************************************************/
inline double MHD1D_cState::E(void) {
    return (E1 + HALF*B0.sqr() + B0*B1);
}

inline double MHD1D_cState::E(void) const {
    return (E1 + HALF*B1.sqr() + B0*B1);
}

/********************************************************
 * MHD1D_cState::h -- Specific enthalpy.                *
 ********************************************************/
inline double MHD1D_cState::h(void) {
    return ((g*E1 - gm1*HALF*dv.sqr()/d - (g-TWO)*HALF*B1.sqr()
	     + B0.sqr() + TWO*B0*B1)/d);
}

inline double MHD1D_cState::h(void) const {
    return ((g*E1 - gm1*HALF*dv.sqr()/d - (g-TWO)*HALF*B1.sqr()
	     + B0.sqr() + TWO*B0*B1)/d);
}

/********************************************************
 * MHD1D_cState::h1 -- Specific perturbative enthalpy.  *
 ********************************************************/
inline double MHD1D_cState::h1(void) {
    return ((g*E1 - gm1*HALF*dv.sqr()/d - (g-TWO)*HALF*B1.sqr())/d);
}

inline double MHD1D_cState::h1(void) const {
    return ((g*E1 - gm1*HALF*dv.sqr()/d - (g-TWO)*HALF*B1.sqr())/d);
}

/********************************************************
 * MHD1D_cState::H -- Total enthalpy.                   *
 ********************************************************/
inline double MHD1D_cState::H(void) {
     return (g*E1 - gm1*HALF*dv.sqr()/d - (g-TWO)*HALF*B1.sqr()
	     + B0.sqr() + TWO*B0*B1);
}

inline double MHD1D_cState::H(void) const {
     return (g*E1 - gm1*HALF*dv.sqr()/d - (g-TWO)*HALF*B1.sqr()
	     + B0.sqr() + TWO*B0*B1);
}

/********************************************************
 * MHD1D_cState::H1 -- Total perturbative enthalpy.     *
 ********************************************************/
inline double MHD1D_cState::H1(void) {
    return (g*E1 - gm1*HALF*dv.sqr()/d - (g-TWO)*HALF*B1.sqr());
}

inline double MHD1D_cState::H1(void) const {
    return (g*E1 - gm1*HALF*dv.sqr()/d - (g-TWO)*HALF*B1.sqr());
}

/********************************************************
 * MHD1D_cState::a -- Sound speed.                      *
 ********************************************************/
inline double MHD1D_cState::a(void) {
    return (sqrt(g*gm1*(E1 - HALF*dv.sqr()/d - HALF*B1.sqr())/d));
}

inline double MHD1D_cState::a(void) const {
    return (sqrt(g*gm1*(E1 - HALF*dv.sqr()/d - HALF*B1.sqr())/d));
}

/********************************************************
 * MHD1D_cState::a2 -- Sound speed squared.             *
 ********************************************************/
inline double MHD1D_cState::a2(void) {
    return (g*gm1*(E1 - HALF*dv.sqr()/d - HALF*B1.sqr())/d);
}

inline double MHD1D_cState::a2(void) const {
    return (g*gm1*(E1 - HALF*dv.sqr()/d - HALF*B1.sqr())/d);
}

/********************************************************
 * MHD1D_cState::Va -- Alfven wave velocity.            *
 ********************************************************/
inline Vector3D MHD1D_cState::Va(void) {
    return (B()/sqrt(d));
}

inline Vector3D MHD1D_cState::Va(void) const {
    return (B()/sqrt(d));
}

/********************************************************
 * MHD1D_cState::Va2 -- Alfven wave speed squared.      *
 ********************************************************/
inline double MHD1D_cState::Va2(void) {
    return ((B().sqr())/d);
}

inline double MHD1D_cState::Va2(void) const {
    return ((B().sqr())/d);
}

/********************************************************
 * MHD1D_cState::s -- Specific entropy.                 *
 ********************************************************/
inline double MHD1D_cState::s(void) {
    return (-gm1i*log(gm1*(E1 - HALF*dv.sqr()/d - HALF*B1.sqr())/pow(d, g)));
}

inline double MHD1D_cState::s(void) const {
    return (-gm1i*log(gm1*(E1 - HALF*dv.sqr()/d - HALF*B1.sqr())/pow(d, g)));
}

/********************************************************
 * MHD1D_cState -- Binary arithmetic operators.         *
 ********************************************************/
inline MHD1D_cState operator +(const MHD1D_cState &U1, const MHD1D_cState &U2) {
  return (MHD1D_cState(U1.d+U2.d,U1.dv+U2.dv,U1.B1+U2.B1,U1.E1+U2.E1));
}

inline MHD1D_cState operator -(const MHD1D_cState &U1, const MHD1D_cState &U2) {
  return (MHD1D_cState(U1.d-U2.d,U1.dv-U2.dv,U1.B1-U2.B1,U1.E1-U2.E1));
}

// Inner product operator.
inline double operator *(const MHD1D_cState &U1, const MHD1D_cState &U2) {
   return (U1.d*U2.d+U1.dv*U2.dv+U1.B1*U2.B1+U1.E1*U2.E1);
}

inline MHD1D_cState operator *(const MHD1D_cState &U, const double &a) {
  return (MHD1D_cState(a*U.d,a*U.dv,a*U.B1,a*U.E1));
}

inline MHD1D_cState operator *(const double &a, const MHD1D_cState &U) {
  return (MHD1D_cState(a*U.d,a*U.dv,a*U.B1,a*U.E1));
}

inline MHD1D_cState operator /(const MHD1D_cState &U, const double &a) {
  return (MHD1D_cState(U.d/a,U.dv/a,U.B1/a,U.E1/a));
}

/********************************************************
 * MHD1D_cState -- Unary arithmetic operators.          *
 ********************************************************/
inline MHD1D_cState operator +(const MHD1D_cState &U) {
  return (MHD1D_cState(U.d,U.dv,U.B1,U.E1));
}

inline MHD1D_cState operator -(const MHD1D_cState &U) {
  return (MHD1D_cState(-U.d,-U.dv,-U.B1,-U.E1));
}

/********************************************************
 * MHD1D_cState -- Shortcut arithmetic operators.       *
 ********************************************************/
inline MHD1D_cState &operator +=(MHD1D_cState &U1, const MHD1D_cState &U2) {
  U1.d += U2.d;
  U1.dv += U2.dv;
  U1.B1 += U2.B1;
  U1.E1 += U2.E1;
  return (U1);
}

inline MHD1D_cState &operator -=(MHD1D_cState &U1, const MHD1D_cState &U2) {
  U1.d -= U2.d;
  U1.dv -= U2.dv;
  U1.B1 -= U2.B1;
  U1.E1 -= U2.E1;
  return (U1);
}

/********************************************************
 * MHD1D_cState -- Relational operators.                *
 ********************************************************/
inline int operator ==(const MHD1D_cState &U1, const MHD1D_cState &U2) {
  return (U1.d == U2.d && U1.dv == U2.dv && U1.B() == U2.B() && U1.E1 == U2.E1);
}

inline int operator !=(const MHD1D_cState &U1, const MHD1D_cState &U2) {
  return (U1.d != U2.d || U1.dv != U2.dv || U1.B() != U2.B() || U1.E1 != U2.E1);
}

/********************************************************
 * MHD1D_cState -- Input-output operators.              *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const MHD1D_cState &U) {
  Vector3D Bt; Bt = U.B();
  out_file.setf(ios::scientific);
  out_file << " " << U.d  << " " << U.dv.x << " " << U.dv.y << " " << U.dv.z
           << " " << Bt.x << " " << Bt.y << " " << Bt.z
	   << " " << U.B0.x << " " << U.B0.y << " " << U.B0.z << " " << U.E1;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, MHD1D_cState &U) {
  Vector3D Bt;
  in_file.setf(ios::skipws);
  in_file >> U.d >> U.dv.x >> U.dv.y >> U.dv.z >> Bt.x >> Bt.y >> Bt.z
	  >> U.B0.x >> U.B0.y >> U.B0.z >> U.E1;
  in_file.unsetf(ios::skipws); U.B1 = Bt - U.B0;
  return (in_file);
}

/********************************************************
 * MHD1D_pState::MHD1D_pState -- Constructor.           *
 ********************************************************/
inline MHD1D_pState::MHD1D_pState(const MHD1D_cState &U) {
  d = U.d; v = U.v(); B1 = U.B1; B0 = U.B0; p = U.p();
}

/********************************************************
 * MHD1D_pState::U -- Conserved solution state.         *
 ********************************************************/
inline MHD1D_cState MHD1D_pState::U(void) {
  return (MHD1D_cState(d, dv(), B1, B0, E()));
}

inline MHD1D_cState MHD1D_pState::U(void) const {
  return (MHD1D_cState(d, dv(), B1, B0, E1()));
}

inline MHD1D_cState MHD1D_pState::U(const MHD1D_pState &W) {
  return (MHD1D_cState(W.d, W.dv(), W.B1, W.B0, W.E1()));
}

inline MHD1D_cState U(const MHD1D_pState &W) {
  return (MHD1D_cState(W.d, W.dv(), W.B1, W.B0, W.E1()));
}

/********************************************************
 * MHD1D_pState::F -- Solution flux.                    *
 ********************************************************/
inline MHD1D_cState MHD1D_pState::F(void) {
  return (MHD1D_cState(d*v.x, d*v.x*v.x + p + HALF*B1.sqr() -
          B1.x*B1.x + B1*B0 - (B0.x*B1.x + B1.x*B0.x), d*v.x*v.y -
          B1.x*B1.y - (B0.x*B1.y + B1.x*B0.y), d*v.x*v.z -
          B1.x*B1.z - (B0.x*B1.z + B1.x*B0.z), v.x*B1.x - B1.x*v.x +
          v.x*B0.x - B0.x*v.x, v.x*B1.y - B1.x*v.y + v.x*B0.y - 
          B0.x*v.y, v.x*B1.z - B1.x*v.z + v.x*B0.z - B0.x*v.z,
          v.x*H1() - (B1*v)*B1.x + (B0*B1)*v.x - (B1*v)*B0.x));
}

inline MHD1D_cState MHD1D_pState::F(void) const {
  return (MHD1D_cState(d*v.x, d*v.x*v.x + p + HALF*B1.sqr() -
          B1.x*B1.x + B1*B0 - (B0.x*B1.x + B1.x*B0.x), d*v.x*v.y -
          B1.x*B1.y - (B0.x*B1.y + B1.x*B0.y), d*v.x*v.z -
          B1.x*B1.z - (B0.x*B1.z + B1.x*B0.z), v.x*B1.x - B1.x*v.x +
          v.x*B0.x - B0.x*v.x, v.x*B1.y - B1.x*v.y + v.x*B0.y - 
          B0.x*v.y, v.x*B1.z - B1.x*v.z + v.x*B0.z - B0.x*v.z,
          v.x*H1() - (B1*v)*B1.x + (B0*B1)*v.x - (B1*v)*B0.x));
}

inline MHD1D_cState MHD1D_pState::F(const MHD1D_pState &W) {
  return (MHD1D_cState(W.d*W.v.x, W.d*W.v.x*W.v.x + W.p + HALF*W.B1.sqr() -
          W.B1.x*W.B1.x + W.B1*W.B0 - (W.B0.x*W.B1.x + W.B1.x*W.B0.x), W.d*W.v.x*W.v.y -
          W.B1.x*W.B1.y - (W.B0.x*W.B1.y + W.B1.x*W.B0.y), W.d*W.v.x*W.v.z -
          W.B1.x*W.B1.z - (W.B0.x*W.B1.z + W.B1.x*W.B0.z), W.v.x*W.B1.x - W.B1.x*W.v.x +
          W.v.x*W.B0.x - W.B0.x*W.v.x, W.v.x*W.B1.y - W.B1.x*W.v.y + W.v.x*W.B0.y -
          W.B0.x*W.v.y, W.v.x*W.B1.z - W.B1.x*W.v.z + W.v.x*W.B0.z - W.B0.x*W.v.z,
          W.v.x*W.H1() - (W.B1*W.v)*W.B1.x + (W.B0*W.B1)*W.v.x - (W.B1*W.v)*W.B0.x));
}

inline MHD1D_cState F(const MHD1D_pState &W) {
  return (MHD1D_cState(W.d*W.v.x, W.d*W.v.x*W.v.x + W.p + HALF*W.B1.sqr() -
          W.B1.x*W.B1.x + W.B1*W.B0 - (W.B0.x*W.B1.x + W.B1.x*W.B0.x), W.d*W.v.x*W.v.y -
          W.B1.x*W.B1.y - (W.B0.x*W.B1.y + W.B1.x*W.B0.y), W.d*W.v.x*W.v.z -
          W.B1.x*W.B1.z - (W.B0.x*W.B1.z + W.B1.x*W.B0.z), W.v.x*W.B1.x - W.B1.x*W.v.x +
          W.v.x*W.B0.x - W.B0.x*W.v.x, W.v.x*W.B1.y - W.B1.x*W.v.y + W.v.x*W.B0.y -
          W.B0.x*W.v.y, W.v.x*W.B1.z - W.B1.x*W.v.z + W.v.x*W.B0.z - W.B0.x*W.v.z,
          W.v.x*W.H1() - (W.B1*W.v)*W.B1.x + (W.B0*W.B1)*W.v.x - (W.B1*W.v)*W.B0.x));
}

/********************************************************
 * MHD1D_pState::lambda -- Eigenvalue(s).               *
 ********************************************************/
inline MHD1D_pState MHD1D_pState::lambda(void) {
  double c = a(), c2 = a2(), v1, v2, cs, cf;
  Vector3D ca = Va();
  v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
  cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
  cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
  return (MHD1D_pState(v.x - cf, v.x - ca.x, v.x - cs, v.x, v.x,
		       v.x + cs, v.x + ca.x, v.x + cf));
}

inline MHD1D_pState MHD1D_pState::lambda(void) const {
  double c = a(), c2 = a2(), v1, v2, cs, cf;
  Vector3D ca = Va();
  v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
  cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
  cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
  return (MHD1D_pState(v.x - cf, v.x - ca.x, v.x - cs, v.x, v.x,
		       v.x + cs, v.x + ca.x, v.x + cf));
}

inline MHD1D_pState MHD1D_pState::lambda(const MHD1D_pState &W) {
  double c = W.a(), c2 = W.a2(), v1, v2, cs, cf;
  Vector3D ca = W.Va();
  v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
  cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
  cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
  return (MHD1D_pState(W.v.x - cf, W.v.x - ca.x, W.v.x - cs, W.v.x, W.v.x,
		       W.v.x + cs, W.v.x + ca.x, W.v.x + cf));
}

inline MHD1D_pState lambda(const MHD1D_pState &W) {
  double c = W.a(), c2 = W.a2(), v1, v2, cs, cf;
  Vector3D ca = W.Va();
  v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
  cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
  cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
  return (MHD1D_pState(W.v.x - cf, W.v.x - ca.x, W.v.x - cs, W.v.x, W.v.x,
		       W.v.x + cs, W.v.x + ca.x, W.v.x + cf));
}

inline double MHD1D_pState::lambda(int index) {
  double c, c2, v1, v2, cs, cf;
  assert( index >= 1 && index <= NUM_VAR_MHD1D );
  switch(index) {
    case 1 :
      c = a(); c2 = a2();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(Va().x)));
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      return (v.x-cf);
    case 2 :
      return (v.x-Va().x);
    case 3 :
      c = a(); c2 = a2();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(Va().x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      return (v.x-cs);
    case 4 :
      return (v.x);
    case 5 :
      return (v.x);
    case 6 :
      c = a(); c2 = a2();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(Va().x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      return (v.x+cs);
    case 7 :
      return (v.x+Va().x);
    case 8 :
      c = a(); c2 = a2();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(Va().x)));
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      return (v.x+cf);
    default:
      return (v.x);
  };
}

inline double MHD1D_pState::lambda(int index) const {
  double c, c2, v1, v2, cs, cf;
  assert( index >= 1 && index <= NUM_VAR_MHD1D );
  switch(index) {
    case 1 :
      c = a(); c2 = a2();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(Va().x)));
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      return (v.x-cf);
    case 2 :
      return (v.x-Va().x);
    case 3 :
      c = a(); c2 = a2();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(Va().x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      return (v.x-cs);
    case 4 :
      return (v.x);
    case 5 :
      return (v.x);
    case 6 :
      c = a(); c2 = a2();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(Va().x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      return (v.x+cs);
    case 7 :
      return (v.x+Va().x);
    case 8 :
      c = a(); c2 = a2();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(Va().x)));
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      return (v.x+cf);
    default:
      return (v.x);
  };
}

inline double lambda(const MHD1D_pState &W, int index) {
  double c, c2, v1, v2, cs, cf;
  assert( index >= 1 && index <= NUM_VAR_MHD1D );
  switch(index) {
    case 1 :
      c = W.a(); c2 = W.a2();
      v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(W.Va().x)));
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      return (W.v.x-cf);
    case 2 :
      return (W.v.x-W.Va().x);
    case 3 :
      c = W.a(); c2 = W.a2();
      v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(W.Va().x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      return (W.v.x-cs);
    case 4 :
      return (W.v.x);
    case 5 :
      return (W.v.x);
    case 6 :
      c = W.a(); c2 = W.a2();
      v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(W.Va().x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      return (W.v.x+cs);
    case 7 :
      return (W.v.x+W.Va().x);
    case 8 :
      c = W.a(); c2 = W.a2();
      v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(W.Va().x)));
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      return (W.v.x+cf);
    default:
      return (W.v.x);
  };
}

/********************************************************
 * MHD1D_pState::rp -- Primitive right eigenvector.     *
 ********************************************************/
inline MHD1D_pState MHD1D_pState::rp(int index) {
  double c, c2, v1, v2, cs, cf, alpha_f, alpha_s, beta_y, beta_z;
  Vector3D ca;
  assert( index >= 1 && index <= NUM_VAR_MHD1D );
  switch(index) {
    case 1 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(alpha_f*d, -alpha_f*cf, alpha_s*cs*beta_y*sgn(ca.x),
                           alpha_s*cs*beta_z*sgn(ca.x), ZERO,
                           alpha_s*c*beta_y*sqrt(d), alpha_s*c*beta_z*sqrt(d),
                           alpha_f*g*p));
    case 2 :
      ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, ZERO, -beta_z/sqrt(TWO), beta_y/sqrt(TWO), ZERO,
			   -beta_z*sqrt(d/TWO), beta_y*sqrt(d/TWO), ZERO));
    case 3 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(alpha_s*d, -alpha_s*cs, -alpha_f*cf*beta_y*sgn(ca.x),
                           -alpha_f*cf*beta_z*sgn(ca.x), ZERO,
                           -alpha_f*c*beta_y*sqrt(d), -alpha_f*c*beta_z*sqrt(d),
                           alpha_s*g*p));
    case 4 :
      return (MHD1D_pState(ONE, Vector3D_ZERO, Vector3D_ZERO, ZERO));
    case 5 :
      return (MHD1D_pState(ZERO, Vector3D_ZERO, Vector3D_NX, ZERO));
    case 6 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(alpha_s*d, alpha_s*cs, alpha_f*cf*beta_y*sgn(ca.x),
                           alpha_f*cf*beta_z*sgn(ca.x), ZERO,
                           -alpha_f*c*beta_y*sqrt(d), -alpha_f*c*beta_z*sqrt(d),
                           alpha_s*g*p));
    case 7 :
      ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, ZERO, -beta_z/sqrt(TWO), beta_y/sqrt(TWO), ZERO,
			   beta_z*sqrt(d/TWO), -beta_y*sqrt(d/TWO), ZERO));
    case 8 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(alpha_f*d, alpha_f*cf, -alpha_s*cs*beta_y*sgn(ca.x),
                           -alpha_s*cs*beta_z*sgn(ca.x), ZERO,
                           alpha_s*c*beta_y*sqrt(d), alpha_s*c*beta_z*sqrt(d),
                           alpha_f*g*p));
    default: 
      return (MHD1D_pState(ONE, Vector3D_ZERO, Vector3D_ZERO, ZERO));
  };
}

inline MHD1D_pState MHD1D_pState::rp(int index) const {
  double c, c2, v1, v2, cs, cf, alpha_f, alpha_s, beta_y, beta_z;
  Vector3D ca;
  assert( index >= 1 && index <= NUM_VAR_MHD1D );
  switch(index) {
    case 1 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(alpha_f*d, -alpha_f*cf, alpha_s*cs*beta_y*sgn(ca.x),
                           alpha_s*cs*beta_z*sgn(ca.x), ZERO,
                           alpha_s*c*beta_y*sqrt(d), alpha_s*c*beta_z*sqrt(d),
                           alpha_f*g*p));
    case 2 :
      ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, ZERO, -beta_z/sqrt(TWO), beta_y/sqrt(TWO), ZERO,
			   -beta_z*sqrt(d/TWO), beta_y*sqrt(d/TWO), ZERO));
    case 3 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(alpha_s*d, -alpha_s*cs, -alpha_f*cf*beta_y*sgn(ca.x),
                           -alpha_f*cf*beta_z*sgn(ca.x), ZERO,
                           -alpha_f*c*beta_y*sqrt(d), -alpha_f*c*beta_z*sqrt(d),
                           alpha_s*g*p));
    case 4 :
      return (MHD1D_pState(ONE, Vector3D_ZERO, Vector3D_ZERO, ZERO));
    case 5 :
      return (MHD1D_pState(ZERO, Vector3D_ZERO, Vector3D_NX, ZERO));
    case 6 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(alpha_s*d, alpha_s*cs, alpha_f*cf*beta_y*sgn(ca.x),
                           alpha_f*cf*beta_z*sgn(ca.x), ZERO,
                           -alpha_f*c*beta_y*sqrt(d), -alpha_f*c*beta_z*sqrt(d),
                           alpha_s*g*p));
    case 7 :
      ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, ZERO, -beta_z/sqrt(TWO), beta_y/sqrt(TWO), ZERO,
			   beta_z*sqrt(d/TWO), -beta_y*sqrt(d/TWO), ZERO));
    case 8 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(alpha_f*d, alpha_f*cf, -alpha_s*cs*beta_y*sgn(ca.x),
                           -alpha_s*cs*beta_z*sgn(ca.x), ZERO,
                           alpha_s*c*beta_y*sqrt(d), alpha_s*c*beta_z*sqrt(d),
                           alpha_f*g*p));
    default: 
      return (MHD1D_pState(ONE, Vector3D_ZERO, Vector3D_ZERO, ZERO));
  };
}

inline MHD1D_pState rp(const MHD1D_pState &W, int index) {
  double c, c2, v1, v2, cs, cf, alpha_f, alpha_s, beta_y, beta_z;
  Vector3D ca;
  assert( index >= 1 && index <= NUM_VAR_MHD1D );
  switch(index) {
    case 1 :
      c = W.a(); c2 = W.a2(); ca = W.Va();
      v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(alpha_f*W.d, -alpha_f*cf, alpha_s*cs*beta_y*sgn(ca.x),
                           alpha_s*cs*beta_z*sgn(ca.x), ZERO,
                           alpha_s*c*beta_y*sqrt(W.d), alpha_s*c*beta_z*sqrt(W.d),
                           alpha_f*W.g*W.p));
    case 2 :
      ca = W.Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, ZERO, -beta_z/sqrt(TWO), beta_y/sqrt(TWO), ZERO,
			   -beta_z*sqrt(W.d/TWO), beta_y*sqrt(W.d/TWO), ZERO));
    case 3 :
      c = W.a(); c2 = W.a2(); ca = W.Va();
      v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(alpha_s*W.d, -alpha_s*cs, -alpha_f*cf*beta_y*sgn(ca.x),
                           -alpha_f*cf*beta_z*sgn(ca.x), ZERO,
                           -alpha_f*c*beta_y*sqrt(W.d), -alpha_f*c*beta_z*sqrt(W.d),
                           alpha_s*W.g*W.p));
    case 4 :
      return (MHD1D_pState(ONE, Vector3D_ZERO, Vector3D_ZERO, ZERO));
    case 5 :
      return (MHD1D_pState(ZERO, Vector3D_ZERO, Vector3D_NX, ZERO));
    case 6 :
      c = W.a(); c2 = W.a2(); ca = W.Va();
      v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(alpha_s*W.d, alpha_s*cs, alpha_f*cf*beta_y*sgn(ca.x),
                           alpha_f*cf*beta_z*sgn(ca.x), ZERO,
                           -alpha_f*c*beta_y*sqrt(W.d), -alpha_f*c*beta_z*sqrt(W.d),
                           alpha_s*W.g*W.p));
    case 7 :
      ca = W.Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, ZERO, -beta_z/sqrt(TWO), beta_y/sqrt(TWO), ZERO,
			   beta_z*sqrt(W.d/TWO), -beta_y*sqrt(W.d/TWO), ZERO));
    case 8 :
      c = W.a(); c2 = W.a2(); ca = W.Va();
      v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(alpha_f*W.d, alpha_f*cf, -alpha_s*cs*beta_y*sgn(ca.x),
                           -alpha_s*cs*beta_z*sgn(ca.x), ZERO,
                           alpha_s*c*beta_y*sqrt(W.d), alpha_s*c*beta_z*sqrt(W.d),
                           alpha_f*W.g*W.p));
    default: 
      return (MHD1D_pState(ONE, Vector3D_ZERO, Vector3D_ZERO, ZERO));
  };
}

/********************************************************
 * MHD1D_pState::rc -- Conserved right eigenvector.     *
 ********************************************************/
inline MHD1D_cState MHD1D_pState::rc(int index) {
  double c, c2, v1, v2, cs, cf,
         alpha_f, alpha_s, beta_y, beta_z, gamma;
  Vector3D ca;
  assert( index >= 1 && index <= NUM_VAR_MHD1D );
  switch(index) {
    case 1 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = alpha_s*(c*sqrt(d)*(beta_y*B1.y+beta_z*B1.z)+
		       d*cs*sgn(ca.x)*(v.y*beta_y+v.z*beta_z));
      return (MHD1D_cState(alpha_f*d, alpha_f*d*(v.x-cf),
			   d*(alpha_f*v.y+alpha_s*cs*beta_y*sgn(ca.x)),
                           d*(alpha_f*v.z+alpha_s*cs*beta_z*sgn(ca.x)), ZERO,
                           alpha_s*c*beta_y*sqrt(d), alpha_s*c*beta_z*sqrt(d),
                           alpha_f*(HALF*d*sqr(v)+g*p*gm1i-d*v.x*cf)+gamma));
    case 2 :
      ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_cState(ZERO, ZERO, -beta_z*d/sqrt(TWO), beta_y*d/sqrt(TWO), ZERO,
			   -beta_z*sqrt(d/TWO), beta_y*sqrt(d/TWO),
			   ((v.z*beta_y-v.y*beta_z)-(B1.y*beta_z-B1.z*beta_y))*sqrt(d/TWO)));
    case 3 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = alpha_f*(c*sqrt(d)*(beta_y*B1.y+beta_z*B1.z)+
		       d*cf*sgn(ca.x)*(v.y*beta_y+v.z*beta_z));
      return (MHD1D_cState(alpha_s*d, alpha_s*d*(v.x-cs),
			   d*(alpha_s*v.y-alpha_f*cf*beta_y*sgn(ca.x)),
                           d*(alpha_s*v.z-alpha_f*cf*beta_z*sgn(ca.x)), ZERO,
                           -alpha_f*c*beta_y*sqrt(d), -alpha_f*c*beta_z*sqrt(d),
                           alpha_s*(HALF*d*sqr(v)+g*p*gm1i-d*v.x*cs)-gamma));
    case 4 :
      return (MHD1D_cState(ONE, v, Vector3D_ZERO, HALF*sqr(v)));
    case 5 :
      return (MHD1D_cState(ZERO, Vector3D_ZERO, Vector3D_NX, B1.x));
    case 6 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = alpha_f*(c*sqrt(d)*(beta_y*B1.y+beta_z*B1.z)-
		       d*cf*sgn(ca.x)*(v.y*beta_y+v.z*beta_z));
      return (MHD1D_cState(alpha_s*d, alpha_s*d*(v.x+cs),
			   d*(alpha_s*v.y+alpha_f*cf*beta_y*sgn(ca.x)),
                           d*(alpha_s*v.z+alpha_f*cf*beta_z*sgn(ca.x)), ZERO,
                           -alpha_f*c*beta_y*sqrt(d), -alpha_f*c*beta_z*sqrt(d),
                           alpha_s*(HALF*d*sqr(v)+g*p*gm1i+d*v.x*cs)-gamma));
    case 7 :
      ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_cState(ZERO, ZERO, -beta_z*d/sqrt(TWO), beta_y*d/sqrt(TWO), ZERO,
			   beta_z*sqrt(d/TWO), -beta_y*sqrt(d/TWO),
			   ((v.z*beta_y-v.y*beta_z)+(B1.y*beta_z-B1.z*beta_y))*sqrt(d/TWO)));
    case 8 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = alpha_s*(c*sqrt(d)*(beta_y*B1.y+beta_z*B1.z)-
		       d*cs*sgn(ca.x)*(v.y*beta_y+v.z*beta_z));
      return (MHD1D_cState(alpha_f*d, alpha_f*d*(v.x+cf),
			   d*(alpha_f*v.y-alpha_s*cs*beta_y*sgn(ca.x)),
                           d*(alpha_f*v.z-alpha_s*cs*beta_z*sgn(ca.x)), ZERO,
                           alpha_s*c*beta_y*sqrt(d), alpha_s*c*beta_z*sqrt(d),
                           alpha_f*(HALF*d*sqr(v)+g*p*gm1i+d*v.x*cf)+gamma));
    default: 
      return (MHD1D_cState(ONE, Vector3D_ZERO, Vector3D_ZERO, ZERO));
  };
}

inline MHD1D_cState MHD1D_pState::rc(int index) const {
  double c, c2, v1, v2, cs, cf,
         alpha_f, alpha_s, beta_y, beta_z, gamma;
  Vector3D ca;
  assert( index >= 1 && index <= NUM_VAR_MHD1D );
  switch(index) {
    case 1 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = alpha_s*(c*sqrt(d)*(beta_y*B1.y+beta_z*B1.z)+
		       d*cs*sgn(ca.x)*(v.y*beta_y+v.z*beta_z));
      return (MHD1D_cState(alpha_f*d, alpha_f*d*(v.x-cf),
			   d*(alpha_f*v.y+alpha_s*cs*beta_y*sgn(ca.x)),
                           d*(alpha_f*v.z+alpha_s*cs*beta_z*sgn(ca.x)), ZERO,
                           alpha_s*c*beta_y*sqrt(d), alpha_s*c*beta_z*sqrt(d),
                           alpha_f*(HALF*d*sqr(v)+g*p*gm1i-d*v.x*cf)+gamma));
    case 2 :
      ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_cState(ZERO, ZERO, -beta_z*d/sqrt(TWO), beta_y*d/sqrt(TWO), ZERO,
			   -beta_z*sqrt(d/TWO), beta_y*sqrt(d/TWO),
			   ((v.z*beta_y-v.y*beta_z)-(B1.y*beta_z-B1.z*beta_y))*sqrt(d/TWO)));
    case 3 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = alpha_f*(c*sqrt(d)*(beta_y*B1.y+beta_z*B1.z)+
		       d*cf*sgn(ca.x)*(v.y*beta_y+v.z*beta_z));
      return (MHD1D_cState(alpha_s*d, alpha_s*d*(v.x-cs),
			   d*(alpha_s*v.y-alpha_f*cf*beta_y*sgn(ca.x)),
                           d*(alpha_s*v.z-alpha_f*cf*beta_z*sgn(ca.x)), ZERO,
                           -alpha_f*c*beta_y*sqrt(d), -alpha_f*c*beta_z*sqrt(d),
                           alpha_s*(HALF*d*sqr(v)+g*p*gm1i-d*v.x*cs)-gamma));
    case 4 :
      return (MHD1D_cState(ONE, v, Vector3D_ZERO, HALF*sqr(v)));
    case 5 :
      return (MHD1D_cState(ZERO, Vector3D_ZERO, Vector3D_NX, B1.x));
    case 6 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = alpha_f*(c*sqrt(d)*(beta_y*B1.y+beta_z*B1.z)-
		       d*cf*sgn(ca.x)*(v.y*beta_y+v.z*beta_z));
      return (MHD1D_cState(alpha_s*d, alpha_s*d*(v.x+cs),
			   d*(alpha_s*v.y+alpha_f*cf*beta_y*sgn(ca.x)),
                           d*(alpha_s*v.z+alpha_f*cf*beta_z*sgn(ca.x)), ZERO,
                           -alpha_f*c*beta_y*sqrt(d), -alpha_f*c*beta_z*sqrt(d),
                           alpha_s*(HALF*d*sqr(v)+g*p*gm1i+d*v.x*cs)-gamma));
    case 7 :
      ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_cState(ZERO, ZERO, -beta_z*d/sqrt(TWO), beta_y*d/sqrt(TWO), ZERO,
			   beta_z*sqrt(d/TWO), -beta_y*sqrt(d/TWO),
			   ((v.z*beta_y-v.y*beta_z)+(B1.y*beta_z-B1.z*beta_y))*sqrt(d/TWO)));
    case 8 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = alpha_s*(c*sqrt(d)*(beta_y*B1.y+beta_z*B1.z)-
		       d*cs*sgn(ca.x)*(v.y*beta_y+v.z*beta_z));
      return (MHD1D_cState(alpha_f*d, alpha_f*d*(v.x+cf),
			   d*(alpha_f*v.y-alpha_s*cs*beta_y*sgn(ca.x)),
                           d*(alpha_f*v.z-alpha_s*cs*beta_z*sgn(ca.x)), ZERO,
                           alpha_s*c*beta_y*sqrt(d), alpha_s*c*beta_z*sqrt(d),
                           alpha_f*(HALF*d*sqr(v)+g*p*gm1i+d*v.x*cf)+gamma));
    default: 
      return (MHD1D_cState(ONE, Vector3D_ZERO, Vector3D_ZERO, ZERO));
  };
}

inline MHD1D_cState rc(const MHD1D_pState &W, int index) {
  double c, c2, v1, v2, cs, cf,
         alpha_f, alpha_s, beta_y, beta_z, gamma;
  Vector3D ca;
  assert( index >= 1 && index <= NUM_VAR_MHD1D );
  switch(index) {
    case 1 :
      c = W.a(); c2 = W.a2(); ca = W.Va();
      v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = alpha_s*(c*sqrt(W.d)*(beta_y*W.B1.y+beta_z*W.B1.z)+
		       W.d*cs*sgn(ca.x)*(W.v.y*beta_y+W.v.z*beta_z));
      return (MHD1D_cState(alpha_f*W.d, alpha_f*W.d*(W.v.x-cf),
			   W.d*(alpha_f*W.v.y+alpha_s*cs*beta_y*sgn(ca.x)),
                           W.d*(alpha_f*W.v.z+alpha_s*cs*beta_z*sgn(ca.x)), ZERO,
                           alpha_s*c*beta_y*sqrt(W.d), alpha_s*c*beta_z*sqrt(W.d),
                           alpha_f*(HALF*W.d*sqr(W.v)+W.g*W.p*W.gm1i-W.d*W.v.x*cf)+gamma));
    case 2 :
      ca = W.Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_cState(ZERO, ZERO, -beta_z*W.d/sqrt(TWO), beta_y*W.d/sqrt(TWO), ZERO,
			   -beta_z*sqrt(W.d/TWO), beta_y*sqrt(W.d/TWO),
			   ((W.v.z*beta_y-W.v.y*beta_z)-(W.B1.y*beta_z-W.B1.z*beta_y))*sqrt(W.d/TWO)));
    case 3 :
      c = W.a(); c2 = W.a2(); ca = W.Va();
      v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = alpha_f*(c*sqrt(W.d)*(beta_y*W.B1.y+beta_z*W.B1.z)+
		       W.d*cf*sgn(ca.x)*(W.v.y*beta_y+W.v.z*beta_z));
      return (MHD1D_cState(alpha_s*W.d, alpha_s*W.d*(W.v.x-cs),
			   W.d*(alpha_s*W.v.y-alpha_f*cf*beta_y*sgn(ca.x)),
                           W.d*(alpha_s*W.v.z-alpha_f*cf*beta_z*sgn(ca.x)), ZERO,
                           -alpha_f*c*beta_y*sqrt(W.d), -alpha_f*c*beta_z*sqrt(W.d),
                           alpha_s*(HALF*W.d*sqr(W.v)+W.g*W.p*W.gm1i-W.d*W.v.x*cs)-gamma));
    case 4 :
      return (MHD1D_cState(ONE, W.v, Vector3D_ZERO, HALF*sqr(W.v)));
    case 5 :
      return (MHD1D_cState(ZERO, Vector3D_ZERO, Vector3D_NX, W.B1.x));
    case 6 :
      c = W.a(); c2 = W.a2(); ca = W.Va();
      v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = alpha_f*(c*sqrt(W.d)*(beta_y*W.B1.y+beta_z*W.B1.z)-
		       W.d*cf*sgn(ca.x)*(W.v.y*beta_y+W.v.z*beta_z));
      return (MHD1D_cState(alpha_s*W.d, alpha_s*W.d*(W.v.x+cs),
			   W.d*(alpha_s*W.v.y+alpha_f*cf*beta_y*sgn(ca.x)),
                           W.d*(alpha_s*W.v.z+alpha_f*cf*beta_z*sgn(ca.x)), ZERO,
                           -alpha_f*c*beta_y*sqrt(W.d), -alpha_f*c*beta_z*sqrt(W.d),
                           alpha_s*(HALF*W.d*sqr(W.v)+W.g*W.p*W.gm1i+W.d*W.v.x*cs)-gamma));
    case 7 :
      ca = W.Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_cState(ZERO, ZERO, -beta_z*W.d/sqrt(TWO), beta_y*W.d/sqrt(TWO), ZERO,
			   beta_z*sqrt(W.d/TWO), -beta_y*sqrt(W.d/TWO),
			   ((W.v.z*beta_y-W.v.y*beta_z)+(W.B1.y*beta_z-W.B1.z*beta_y))*sqrt(W.d/TWO)));
    case 8 :
      c = W.a(); c2 = W.a2(); ca = W.Va();
      v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = alpha_s*(c*sqrt(W.d)*(beta_y*W.B1.y+beta_z*W.B1.z)-
		       W.d*cs*sgn(ca.x)*(W.v.y*beta_y+W.v.z*beta_z));
      return (MHD1D_cState(alpha_f*W.d, alpha_f*W.d*(W.v.x+cf),
			   W.d*(alpha_f*W.v.y-alpha_s*cs*beta_y*sgn(ca.x)),
                           W.d*(alpha_f*W.v.z-alpha_s*cs*beta_z*sgn(ca.x)), ZERO,
                           alpha_s*c*beta_y*sqrt(W.d), alpha_s*c*beta_z*sqrt(W.d),
                           alpha_f*(HALF*W.d*sqr(W.v)+W.g*W.p*W.gm1i+W.d*W.v.x*cf)+gamma));
    default: 
      return (MHD1D_cState(ONE, Vector3D_ZERO, Vector3D_ZERO, ZERO));
  };
}

/********************************************************
 * MHD1D_pState::lp -- Primitive left eigenvector.      *
 ********************************************************/
inline MHD1D_pState MHD1D_pState::lp(int index) {
  double c, c2, v1, v2, cs, cf, alpha_f, alpha_s, beta_y, beta_z;
  Vector3D ca;
  assert( index >= 1 && index <= NUM_VAR_MHD1D );
  switch(index) {
    case 1 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, -alpha_f*cf/(TWO*c2),
			   alpha_s*cs*beta_y*sgn(ca.x)/(TWO*c2),
                           alpha_s*cs*beta_z*sgn(ca.x)/(TWO*c2), ZERO,
                           alpha_s*beta_y/(TWO*c*sqrt(d)),
			   alpha_s*beta_z/(TWO*c*sqrt(d)), alpha_f/(TWO*d*c2)));
    case 2 :
      ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, ZERO, -beta_z/sqrt(TWO), beta_y/sqrt(TWO), ZERO,
			   -beta_z/sqrt(TWO*d), beta_y/sqrt(TWO*d), ZERO));
    case 3 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, -alpha_s*cs/(TWO*c2),
                           -alpha_f*cf*beta_y*sgn(ca.x)/(TWO*c2),
                           -alpha_f*cf*beta_z*sgn(ca.x)/(TWO*c2), ZERO,
                           -alpha_f*beta_y/(TWO*c*sqrt(d)),
                           -alpha_f*beta_z/(TWO*c*sqrt(d)), alpha_s/(TWO*d*c2)));
    case 4 :
      return (MHD1D_pState(ONE, Vector3D_ZERO, Vector3D_ZERO, -ONE/a2()));
    case 5 :
      return (MHD1D_pState(ZERO, Vector3D_ZERO, Vector3D_NX, ZERO));
    case 6 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, alpha_s*cs/(TWO*c2),
                           alpha_f*cf*beta_y*sgn(ca.x)/(TWO*c2),
                           alpha_f*cf*beta_z*sgn(ca.x)/(TWO*c2), ZERO,
                           -alpha_f*beta_y/(TWO*c*sqrt(d)),
                           -alpha_f*beta_z/(TWO*c*sqrt(d)), alpha_s/(TWO*d*c2)));
    case 7 :
      ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, ZERO, -beta_z/sqrt(TWO), beta_y/sqrt(TWO), ZERO,
			   beta_z/sqrt(TWO*d), -beta_y/sqrt(TWO*d), ZERO));
    case 8 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, alpha_f*cf/(TWO*c2),
			   -alpha_s*cs*beta_y*sgn(ca.x)/(TWO*c2),
                           -alpha_s*cs*beta_z*sgn(ca.x)/(TWO*c2), ZERO,
                           alpha_s*beta_y/(TWO*c*sqrt(d)),
			   alpha_s*beta_z/(TWO*c*sqrt(d)), alpha_f/(TWO*d*c2)));
    default: 
      return (MHD1D_pState(ONE, Vector3D_ZERO, Vector3D_ZERO, ZERO));
  };
}

inline MHD1D_pState MHD1D_pState::lp(int index) const {
  double c, c2, v1, v2, cs, cf, alpha_f, alpha_s, beta_y, beta_z;
  Vector3D ca;
  assert( index >= 1 && index <= NUM_VAR_MHD1D );
  switch(index) {
    case 1 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, -alpha_f*cf/(TWO*c2),
			   alpha_s*cs*beta_y*sgn(ca.x)/(TWO*c2),
                           alpha_s*cs*beta_z*sgn(ca.x)/(TWO*c2), ZERO,
                           alpha_s*beta_y/(TWO*c*sqrt(d)),
			   alpha_s*beta_z/(TWO*c*sqrt(d)), alpha_f/(TWO*d*c2)));
    case 2 :
      ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, ZERO, -beta_z/sqrt(TWO), beta_y/sqrt(TWO), ZERO,
			   -beta_z/sqrt(TWO*d), beta_y/sqrt(TWO*d), ZERO));
    case 3 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, -alpha_s*cs/(TWO*c2),
                           -alpha_f*cf*beta_y*sgn(ca.x)/(TWO*c2),
                           -alpha_f*cf*beta_z*sgn(ca.x)/(TWO*c2), ZERO,
                           -alpha_f*beta_y/(TWO*c*sqrt(d)),
                           -alpha_f*beta_z/(TWO*c*sqrt(d)), alpha_s/(TWO*d*c2)));
    case 4 :
      return (MHD1D_pState(ONE, Vector3D_ZERO, Vector3D_ZERO, -ONE/a2()));
    case 5 :
      return (MHD1D_pState(ZERO, Vector3D_ZERO, Vector3D_NX, ZERO));
    case 6 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, alpha_s*cs/(TWO*c2),
                           alpha_f*cf*beta_y*sgn(ca.x)/(TWO*c2),
                           alpha_f*cf*beta_z*sgn(ca.x)/(TWO*c2), ZERO,
                           -alpha_f*beta_y/(TWO*c*sqrt(d)),
                           -alpha_f*beta_z/(TWO*c*sqrt(d)), alpha_s/(TWO*d*c2)));
    case 7 :
      ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, ZERO, -beta_z/sqrt(TWO), beta_y/sqrt(TWO), ZERO,
			   beta_z/sqrt(TWO*d), -beta_y/sqrt(TWO*d), ZERO));
    case 8 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, alpha_f*cf/(TWO*c2),
			   -alpha_s*cs*beta_y*sgn(ca.x)/(TWO*c2),
                           -alpha_s*cs*beta_z*sgn(ca.x)/(TWO*c2), ZERO,
                           alpha_s*beta_y/(TWO*c*sqrt(d)),
			   alpha_s*beta_z/(TWO*c*sqrt(d)), alpha_f/(TWO*d*c2)));
    default: 
      return (MHD1D_pState(ONE, Vector3D_ZERO, Vector3D_ZERO, ZERO));
  };
}

inline MHD1D_pState lp(const MHD1D_pState &W, int index) {
  double c, c2, v1, v2, cs, cf, alpha_f, alpha_s, beta_y, beta_z;
  Vector3D ca;
  assert( index >= 1 && index <= NUM_VAR_MHD1D );
  switch(index) {
    case 1 :
      c = W.a(); c2 = W.a2(); ca = W.Va();
      v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, -alpha_f*cf/(TWO*c2),
			   alpha_s*cs*beta_y*sgn(ca.x)/(TWO*c2),
                           alpha_s*cs*beta_z*sgn(ca.x)/(TWO*c2), ZERO,
                           alpha_s*beta_y/(TWO*c*sqrt(W.d)),
			   alpha_s*beta_z/(TWO*c*sqrt(W.d)), alpha_f/(TWO*W.d*c2)));
    case 2 :
      ca = W.Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, ZERO, -beta_z/sqrt(TWO), beta_y/sqrt(TWO), ZERO,
			   -beta_z/sqrt(TWO*W.d), beta_y/sqrt(TWO*W.d), ZERO));
    case 3 :
      c = W.a(); c2 = W.a2(); ca = W.Va();
      v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, -alpha_s*cs/(TWO*c2),
                           -alpha_f*cf*beta_y*sgn(ca.x)/(TWO*c2),
                           -alpha_f*cf*beta_z*sgn(ca.x)/(TWO*c2), ZERO,
                           -alpha_f*beta_y/(TWO*c*sqrt(W.d)),
                           -alpha_f*beta_z/(TWO*c*sqrt(W.d)), alpha_s/(TWO*W.d*c2)));
    case 4 :
      return (MHD1D_pState(ONE, Vector3D_ZERO, Vector3D_ZERO, -ONE/W.a2()));
    case 5 :
      return (MHD1D_pState(ZERO, Vector3D_ZERO, Vector3D_NX, ZERO));
    case 6 :
      c = W.a(); c2 = W.a2(); ca = W.Va();
      v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, alpha_s*cs/(TWO*c2),
                           alpha_f*cf*beta_y*sgn(ca.x)/(TWO*c2),
                           alpha_f*cf*beta_z*sgn(ca.x)/(TWO*c2), ZERO,
                           -alpha_f*beta_y/(TWO*c*sqrt(W.d)),
                           -alpha_f*beta_z/(TWO*c*sqrt(W.d)), alpha_s/(TWO*W.d*c2)));
    case 7 :
      ca = W.Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, ZERO, -beta_z/sqrt(TWO), beta_y/sqrt(TWO), ZERO,
			   beta_z/sqrt(TWO*W.d), -beta_y/sqrt(TWO*W.d), ZERO));
    case 8 :
      c = W.a(); c2 = W.a2(); ca = W.Va();
      v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_pState(ZERO, alpha_f*cf/(TWO*c2),
			   -alpha_s*cs*beta_y*sgn(ca.x)/(TWO*c2),
                           -alpha_s*cs*beta_z*sgn(ca.x)/(TWO*c2), ZERO,
                           alpha_s*beta_y/(TWO*c*sqrt(W.d)),
			   alpha_s*beta_z/(TWO*c*sqrt(W.d)), alpha_f/(TWO*W.d*c2)));
    default: 
      return (MHD1D_pState(ONE, Vector3D_ZERO, Vector3D_ZERO, ZERO));
  };
}

/********************************************************
 * MHD1D_pState::lc -- Conserved left eigenvector.      *
 ********************************************************/
inline MHD1D_cState MHD1D_pState::lc(int index) {
  double c, c2, v1, v2, cs, cf,
         alpha_f, alpha_s, beta_y, beta_z, gamma;
  Vector3D ca;
  assert( index >= 1 && index <= NUM_VAR_MHD1D );
  switch(index) {
    case 1 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = -HALF*alpha_s*cs*sgn(ca.x)*(v.y*beta_y+v.z*beta_z)/(d*c2);
      return (MHD1D_cState(HALF*alpha_f*(HALF*gm1*sqr(v)+v.x*cf)/(d*c2)+gamma,
			   -HALF*alpha_f*(gm1*v.x+cf)/(d*c2),
			   -HALF*(gm1*alpha_f*v.y-alpha_s*cs*beta_y*sgn(ca.x))/(d*c2),
                           -HALF*(gm1*alpha_f*v.z-alpha_s*cs*beta_z*sgn(ca.x))/(d*c2),
			   -HALF*gm1*alpha_f*B1.x/(d*c2),
                           HALF*(alpha_s*c*beta_y*sqrt(d)-gm1*alpha_f*B1.y)/(d*c2),
			   HALF*(alpha_s*c*beta_z*sqrt(d)-gm1*alpha_f*B1.z)/(d*c2),
                           HALF*gm1*alpha_f/(d*c2)));
    case 2 :
      ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_cState((v.y*beta_z - v.z*beta_y)/(sqrt(TWO)*d),
			   ZERO, -beta_z/(sqrt(TWO)*d), beta_y/(sqrt(TWO)*d), ZERO,
			   -beta_z/sqrt(TWO*d), beta_y/sqrt(TWO*d), ZERO));
    case 3 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = HALF*alpha_f*cf*sgn(ca.x)*(v.y*beta_y+v.z*beta_z)/(d*c2);
      return (MHD1D_cState(HALF*alpha_s*(HALF*gm1*sqr(v)+v.x*cs)/(d*c2)+gamma,
			   -HALF*alpha_s*(gm1*v.x+cs)/(d*c2),
			   -HALF*(gm1*alpha_s*v.y+alpha_f*cf*beta_y*sgn(ca.x))/(d*c2),
                           -HALF*(gm1*alpha_s*v.z+alpha_f*cf*beta_z*sgn(ca.x))/(d*c2),
			   -HALF*gm1*alpha_s*B1.x/(d*c2),
                           -HALF*(alpha_f*c*beta_y*sqrt(d)+gm1*alpha_s*B1.y)/(d*c2),
			   -HALF*(alpha_f*c*beta_z*sqrt(d)+gm1*alpha_s*B1.z)/(d*c2),
                           HALF*gm1*alpha_s/(d*c2)));
    case 4 :
      c2 = a2();
      return (MHD1D_cState(ONE-HALF*gm1*sqr(v)/c2, gm1*v.x/c2,
			   gm1*v.y/c2, gm1*v.z/c2, gm1*B1.x/c2,
			   gm1*B1.y/c2, gm1*B1.z/c2, -gm1/c2));
    case 5 :
      return (MHD1D_cState(ZERO, Vector3D_ZERO, Vector3D_NX, ZERO));
    case 6 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = -HALF*alpha_f*cf*sgn(ca.x)*(v.y*beta_y+v.z*beta_z)/(d*c2);
      return (MHD1D_cState(HALF*alpha_s*(HALF*gm1*sqr(v)-v.x*cs)/(d*c2)+gamma,
			   -HALF*alpha_s*(gm1*v.x-cs)/(d*c2),
			   -HALF*(gm1*alpha_s*v.y-alpha_f*cf*beta_y*sgn(ca.x))/(d*c2),
                           -HALF*(gm1*alpha_s*v.z-alpha_f*cf*beta_z*sgn(ca.x))/(d*c2),
			   -HALF*gm1*alpha_s*B1.x/(d*c2),
                           -HALF*(alpha_f*c*beta_y*sqrt(d)+gm1*alpha_s*B1.y)/(d*c2),
			   -HALF*(alpha_f*c*beta_z*sqrt(d)+gm1*alpha_s*B1.z)/(d*c2),
                           HALF*gm1*alpha_s/(d*c2)));
    case 7 :
      ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_cState((v.y*beta_z - v.z*beta_y)/(sqrt(TWO)*d),
			   ZERO, -beta_z/(sqrt(TWO)*d), beta_y/(sqrt(TWO)*d), ZERO,
			   beta_z/sqrt(TWO*d), -beta_y/sqrt(TWO*d), ZERO));
    case 8 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = HALF*alpha_s*cs*sgn(ca.x)*(v.y*beta_y+v.z*beta_z)/(d*c2);
      return (MHD1D_cState(HALF*alpha_f*(HALF*gm1*sqr(v)-v.x*cf)/(d*c2)+gamma,
			   -HALF*alpha_f*(gm1*v.x-cf)/(d*c2),
			   -HALF*(gm1*alpha_f*v.y+alpha_s*cs*beta_y*sgn(ca.x))/(d*c2),
                           -HALF*(gm1*alpha_f*v.z+alpha_s*cs*beta_z*sgn(ca.x))/(d*c2),
			   -HALF*gm1*alpha_f*B1.x/(d*c2),
                           HALF*(alpha_s*c*beta_y*sqrt(d)-gm1*alpha_f*B1.y)/(d*c2),
			   HALF*(alpha_s*c*beta_z*sqrt(d)-gm1*alpha_f*B1.z)/(d*c2),
                           HALF*gm1*alpha_f/(d*c2)));
    default: 
      return (MHD1D_cState(ONE, Vector3D_ZERO, Vector3D_ZERO, ZERO));
  };
}

inline MHD1D_cState MHD1D_pState::lc(int index) const {
  double c, c2, v1, v2, cs, cf,
         alpha_f, alpha_s, beta_y, beta_z, gamma;
  Vector3D ca;
  assert( index >= 1 && index <= NUM_VAR_MHD1D );
  switch(index) {
    case 1 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = -HALF*alpha_s*cs*sgn(ca.x)*(v.y*beta_y+v.z*beta_z)/(d*c2);
      return (MHD1D_cState(HALF*alpha_f*(HALF*gm1*sqr(v)+v.x*cf)/(d*c2)+gamma,
			   -HALF*alpha_f*(gm1*v.x+cf)/(d*c2),
			   -HALF*(gm1*alpha_f*v.y-alpha_s*cs*beta_y*sgn(ca.x))/(d*c2),
                           -HALF*(gm1*alpha_f*v.z-alpha_s*cs*beta_z*sgn(ca.x))/(d*c2),
			   -HALF*gm1*alpha_f*B1.x/(d*c2),
                           HALF*(alpha_s*c*beta_y*sqrt(d)-gm1*alpha_f*B1.y)/(d*c2),
			   HALF*(alpha_s*c*beta_z*sqrt(d)-gm1*alpha_f*B1.z)/(d*c2),
                           HALF*gm1*alpha_f/(d*c2)));
    case 2 :
      ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_cState((v.y*beta_z - v.z*beta_y)/(sqrt(TWO)*d),
			   ZERO, -beta_z/(sqrt(TWO)*d), beta_y/(sqrt(TWO)*d), ZERO,
			   -beta_z/sqrt(TWO*d), beta_y/sqrt(TWO*d), ZERO));
    case 3 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = HALF*alpha_f*cf*sgn(ca.x)*(v.y*beta_y+v.z*beta_z)/(d*c2);
      return (MHD1D_cState(HALF*alpha_s*(HALF*gm1*sqr(v)+v.x*cs)/(d*c2)+gamma,
			   -HALF*alpha_s*(gm1*v.x+cs)/(d*c2),
			   -HALF*(gm1*alpha_s*v.y+alpha_f*cf*beta_y*sgn(ca.x))/(d*c2),
                           -HALF*(gm1*alpha_s*v.z+alpha_f*cf*beta_z*sgn(ca.x))/(d*c2),
			   -HALF*gm1*alpha_s*B1.x/(d*c2),
                           -HALF*(alpha_f*c*beta_y*sqrt(d)+gm1*alpha_s*B1.y)/(d*c2),
			   -HALF*(alpha_f*c*beta_z*sqrt(d)+gm1*alpha_s*B1.z)/(d*c2),
                           HALF*gm1*alpha_s/(d*c2)));
    case 4 :
      c2 = a2();
      return (MHD1D_cState(ONE-HALF*gm1*sqr(v)/c2, gm1*v.x/c2,
			   gm1*v.y/c2, gm1*v.z/c2, gm1*B1.x/c2,
			   gm1*B1.y/c2, gm1*B1.z/c2, -gm1/c2));
    case 5 :
      return (MHD1D_cState(ZERO, Vector3D_ZERO, Vector3D_NX, ZERO));
    case 6 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = -HALF*alpha_f*cf*sgn(ca.x)*(v.y*beta_y+v.z*beta_z)/(d*c2);
      return (MHD1D_cState(HALF*alpha_s*(HALF*gm1*sqr(v)-v.x*cs)/(d*c2)+gamma,
			   -HALF*alpha_s*(gm1*v.x-cs)/(d*c2),
			   -HALF*(gm1*alpha_s*v.y-alpha_f*cf*beta_y*sgn(ca.x))/(d*c2),
                           -HALF*(gm1*alpha_s*v.z-alpha_f*cf*beta_z*sgn(ca.x))/(d*c2),
			   -HALF*gm1*alpha_s*B1.x/(d*c2),
                           -HALF*(alpha_f*c*beta_y*sqrt(d)+gm1*alpha_s*B1.y)/(d*c2),
			   -HALF*(alpha_f*c*beta_z*sqrt(d)+gm1*alpha_s*B1.z)/(d*c2),
                           HALF*gm1*alpha_s/(d*c2)));
    case 7 :
      ca = Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_cState((v.y*beta_z - v.z*beta_y)/(sqrt(TWO)*d),
			   ZERO, -beta_z/(sqrt(TWO)*d), beta_y/(sqrt(TWO)*d), ZERO,
			   beta_z/sqrt(TWO*d), -beta_y/sqrt(TWO*d), ZERO));
    case 8 :
      c = a(); c2 = a2(); ca = Va();
      v1 = c2 + Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = HALF*alpha_s*cs*sgn(ca.x)*(v.y*beta_y+v.z*beta_z)/(d*c2);
      return (MHD1D_cState(HALF*alpha_f*(HALF*gm1*sqr(v)-v.x*cf)/(d*c2)+gamma,
			   -HALF*alpha_f*(gm1*v.x-cf)/(d*c2),
			   -HALF*(gm1*alpha_f*v.y+alpha_s*cs*beta_y*sgn(ca.x))/(d*c2),
                           -HALF*(gm1*alpha_f*v.z+alpha_s*cs*beta_z*sgn(ca.x))/(d*c2),
			   -HALF*gm1*alpha_f*B1.x/(d*c2),
                           HALF*(alpha_s*c*beta_y*sqrt(d)-gm1*alpha_f*B1.y)/(d*c2),
			   HALF*(alpha_s*c*beta_z*sqrt(d)-gm1*alpha_f*B1.z)/(d*c2),
                           HALF*gm1*alpha_f/(d*c2)));
    default: 
      return (MHD1D_cState(ONE, Vector3D_ZERO, Vector3D_ZERO, ZERO));
  };
}

inline MHD1D_cState lc(const MHD1D_pState &W, int index) {
  double c, c2, v1, v2, cs, cf,
         alpha_f, alpha_s, beta_y, beta_z, gamma;
  Vector3D ca;
  assert( index >= 1 && index <= NUM_VAR_MHD1D );
  switch(index) {
    case 1 :
      c = W.a(); c2 = W.a2(); ca = W.Va();
      v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = -HALF*alpha_s*cs*sgn(ca.x)*(W.v.y*beta_y+W.v.z*beta_z)/(W.d*c2);
      return (MHD1D_cState(HALF*alpha_f*(HALF*W.gm1*sqr(W.v)+W.v.x*cf)/(W.d*c2)+gamma,
			   -HALF*alpha_f*(W.gm1*W.v.x+cf)/(W.d*c2),
			   -HALF*(W.gm1*alpha_f*W.v.y-alpha_s*cs*beta_y*sgn(ca.x))/(W.d*c2),
                           -HALF*(W.gm1*alpha_f*W.v.z-alpha_s*cs*beta_z*sgn(ca.x))/(W.d*c2),
			   -HALF*W.gm1*alpha_f*W.B1.x/(W.d*c2),
                           HALF*(alpha_s*c*beta_y*sqrt(W.d)-W.gm1*alpha_f*W.B1.y)/(W.d*c2),
			   HALF*(alpha_s*c*beta_z*sqrt(W.d)-W.gm1*alpha_f*W.B1.z)/(W.d*c2),
                           HALF*W.gm1*alpha_f/(W.d*c2)));
    case 2 :
      ca = W.Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_cState((W.v.y*beta_z - W.v.z*beta_y)/(sqrt(TWO)*W.d),
			   ZERO, -beta_z/(sqrt(TWO)*W.d), beta_y/(sqrt(TWO)*W.d), ZERO,
			   -beta_z/sqrt(TWO*W.d), beta_y/sqrt(TWO*W.d), ZERO));
    case 3 :
      c = W.a(); c2 = W.a2(); ca = W.Va();
      v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = HALF*alpha_f*cf*sgn(ca.x)*(W.v.y*beta_y+W.v.z*beta_z)/(W.d*c2);
      return (MHD1D_cState(HALF*alpha_s*(HALF*W.gm1*sqr(W.v)+W.v.x*cs)/(W.d*c2)+gamma,
			   -HALF*alpha_s*(W.gm1*W.v.x+cs)/(W.d*c2),
			   -HALF*(W.gm1*alpha_s*W.v.y+alpha_f*cf*beta_y*sgn(ca.x))/(W.d*c2),
                           -HALF*(W.gm1*alpha_s*W.v.z+alpha_f*cf*beta_z*sgn(ca.x))/(W.d*c2),
			   -HALF*W.gm1*alpha_s*W.B1.x/(W.d*c2),
                           -HALF*(alpha_f*c*beta_y*sqrt(W.d)+W.gm1*alpha_s*W.B1.y)/(W.d*c2),
			   -HALF*(alpha_f*c*beta_z*sqrt(W.d)+W.gm1*alpha_s*W.B1.z)/(W.d*c2),
                           HALF*W.gm1*alpha_s/(W.d*c2)));
    case 4 :
      c2 = W.a2();
      return (MHD1D_cState(ONE-HALF*W.gm1*sqr(W.v)/c2, W.gm1*W.v.x/c2,
			   W.gm1*W.v.y/c2, W.gm1*W.v.z/c2, W.gm1*W.B1.x/c2,
			   W.gm1*W.B1.y/c2, W.gm1*W.B1.z/c2, -W.gm1/c2));
    case 5 :
      return (MHD1D_cState(ZERO, Vector3D_ZERO, Vector3D_NX, ZERO));
    case 6 :
      c = W.a(); c2 = W.a2(); ca = W.Va();
      v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = -HALF*alpha_f*cf*sgn(ca.x)*(W.v.y*beta_y+W.v.z*beta_z)/(W.d*c2);
      return (MHD1D_cState(HALF*alpha_s*(HALF*W.gm1*sqr(W.v)-W.v.x*cs)/(W.d*c2)+gamma,
			   -HALF*alpha_s*(W.gm1*W.v.x-cs)/(W.d*c2),
			   -HALF*(W.gm1*alpha_s*W.v.y-alpha_f*cf*beta_y*sgn(ca.x))/(W.d*c2),
                           -HALF*(W.gm1*alpha_s*W.v.z-alpha_f*cf*beta_z*sgn(ca.x))/(W.d*c2),
			   -HALF*W.gm1*alpha_s*W.B1.x/(W.d*c2),
                           -HALF*(alpha_f*c*beta_y*sqrt(W.d)+W.gm1*alpha_s*W.B1.y)/(W.d*c2),
			   -HALF*(alpha_f*c*beta_z*sqrt(W.d)+W.gm1*alpha_s*W.B1.z)/(W.d*c2),
                           HALF*W.gm1*alpha_s/(W.d*c2)));
    case 7 :
      ca = W.Va(); v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      return (MHD1D_cState((W.v.y*beta_z - W.v.z*beta_y)/(sqrt(TWO)*W.d),
			   ZERO, -beta_z/(sqrt(TWO)*W.d), beta_y/(sqrt(TWO)*W.d), ZERO,
			   beta_z/sqrt(TWO*W.d), -beta_y/sqrt(TWO*W.d), ZERO));
    case 8 :
      c = W.a(); c2 = W.a2(); ca = W.Va();
      v1 = c2 + W.Va2(); v2 = sqrt(max(ZERO, v1*v1-FOUR*c2*sqr(ca.x)));
      cs = sqrt(max(ZERO, HALF*(v1-v2))); cs = min(cs, c);
      cf = sqrt(HALF*(v1+v2)); cf = max(cf, c);
      v1 = cf*cf-cs*cs; v2 = sqrt(sqr(ca.y) + sqr(ca.z));
      if (v1 >= TOLER) {
         alpha_f = sqrt(max(ZERO, c2-cs*cs)/v1); alpha_s = sqrt(max(ZERO, cf*cf-c2)/v1);
      } else if (ca.x <= c) {
         alpha_f = ONE; alpha_s = ZERO;
      } else {
         alpha_f = ZERO; alpha_s = ONE;
      } /* endif */
      if (v2 >= TOLER) {
         beta_y = ca.y/v2; beta_z = ca.z/v2;
      } else {
         beta_y = ONE/sqrt(TWO); beta_z = ONE/sqrt(TWO);
      } /* endif */
      gamma = HALF*alpha_s*cs*sgn(ca.x)*(W.v.y*beta_y+W.v.z*beta_z)/(W.d*c2);
      return (MHD1D_cState(HALF*alpha_f*(HALF*W.gm1*sqr(W.v)-W.v.x*cf)/(W.d*c2)+gamma,
			   -HALF*alpha_f*(W.gm1*W.v.x-cf)/(W.d*c2),
			   -HALF*(W.gm1*alpha_f*W.v.y+alpha_s*cs*beta_y*sgn(ca.x))/(W.d*c2),
                           -HALF*(W.gm1*alpha_f*W.v.z+alpha_s*cs*beta_z*sgn(ca.x))/(W.d*c2),
			   -HALF*W.gm1*alpha_f*W.B1.x/(W.d*c2),
                           HALF*(alpha_s*c*beta_y*sqrt(W.d)-W.gm1*alpha_f*W.B1.y)/(W.d*c2),
			   HALF*(alpha_s*c*beta_z*sqrt(W.d)-W.gm1*alpha_f*W.B1.z)/(W.d*c2),
                           HALF*W.gm1*alpha_f/(W.d*c2)));
    default: 
      return (MHD1D_cState(ONE, Vector3D_ZERO, Vector3D_ZERO, ZERO));
  };
}

/********************************************************
 * MHD1D_cState::MHD1D_cState -- Constructor.           *
 ********************************************************/
inline MHD1D_cState::MHD1D_cState(const MHD1D_pState &W) {
  d = W.d; dv = W.dv(); B1 = W.B1; B0 = W.B0; E1 = W.E1();
}

/********************************************************
 * MHD1D_cState::W -- Primitive solution state.         *
 ********************************************************/
inline MHD1D_pState MHD1D_cState::W(void) {
  return (MHD1D_pState(d, v(), B1, B0, p()));
}

inline MHD1D_pState MHD1D_cState::W(void) const {
  return (MHD1D_pState(d, v(), B1, B0, p()));
}

inline MHD1D_pState MHD1D_cState::W(const MHD1D_cState &U) {
  return (MHD1D_pState(U.d, U.v(), U.B1, U.B0, U.p()));
}

inline MHD1D_pState W(const MHD1D_cState &U) {
  return (MHD1D_pState(U.d, U.v(), U.B1, U.B0, U.p()));
}

/********************************************************
 * Useful 1D MHD state constants.                       *
 ********************************************************/
const MHD1D_pState MHD1D_W_REF(ONE, Vector3D_ZERO, Vector3D_ZERO,
			            Vector3D_ZERO, ONE);
const MHD1D_pState MHD1D_W_ZERO(ZERO, Vector3D_ZERO, Vector3D_ZERO,
				      Vector3D_ZERO, ZERO);
const MHD1D_cState MHD1D_U_REF(MHD1D_W_REF);
const MHD1D_cState MHD1D_U_ZERO(MHD1D_W_ZERO);

/********************************************************
 * MHD1DState -- External subroutines.                  *
 ********************************************************/

extern MHD1D_pState RoeAverage(const MHD1D_pState &Wl,
	      	               const MHD1D_pState &Wr);

extern MHD1D_pState WaveSpeedPos(const MHD1D_pState &lambda_a,
                                 const MHD1D_pState &lambda_l,
                                 const MHD1D_pState &lambda_r);

extern MHD1D_pState WaveSpeedNeg(const MHD1D_pState &lambda_a,
                                 const MHD1D_pState &lambda_l,
                                 const MHD1D_pState &lambda_r);

extern MHD1D_pState WaveSpeedAbs(const MHD1D_pState &lambda_a,
                                 const MHD1D_pState &lambda_l,
                                 const MHD1D_pState &lambda_r);

extern MHD1D_pState HartenFixPos(const MHD1D_pState &lambda_a,
                                 const MHD1D_pState &lambda_l,
                                 const MHD1D_pState &lambda_r);

extern MHD1D_pState HartenFixNeg(const MHD1D_pState &lambda_a,
                                 const MHD1D_pState &lambda_l,
                                 const MHD1D_pState &lambda_r);

extern MHD1D_pState HartenFixAbs(const MHD1D_pState &lambda_a,
                                 const MHD1D_pState &lambda_l,
                                 const MHD1D_pState &lambda_r);

extern MHD1D_cState FluxRoe(const MHD1D_pState &Wl,
	      	            const MHD1D_pState &Wr);

extern MHD1D_cState FluxRoe(const MHD1D_cState &Ul,
	      	            const MHD1D_cState &Ur);

extern MHD1D_cState FluxRusanov(const MHD1D_pState &Wl,
	      	            const MHD1D_pState &Wr);

extern MHD1D_cState FluxRusanov(const MHD1D_cState &Ul,
	      	            const MHD1D_cState &Ur);

extern MHD1D_cState FluxHLLE(const MHD1D_pState &Wl,
	      	             const MHD1D_pState &Wr);

extern MHD1D_cState FluxHLLE(const MHD1D_cState &Ul,
	      	             const MHD1D_cState &Ur);

extern MHD1D_cState FluxLinde(const MHD1D_pState &Wl,
	      	              const MHD1D_pState &Wr);

extern MHD1D_cState FluxLinde(const MHD1D_cState &Ul,
	      	              const MHD1D_cState &Ur);

#endif /* _MHD1D_STATE_INCLUDED  */
