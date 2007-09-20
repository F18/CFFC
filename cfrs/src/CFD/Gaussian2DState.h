/* Gaussian2DState.h:  Header file defining 2D Euler Solution State Classes. */

#ifndef _GAUSSIAN2D_STATE_INCLUDED
#define _GAUSSIAN2D_STATE_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

/* Include math macro, CFD, 2D vector, and gas constant header files. */

#ifndef _CFD_INCLUDED
#include "CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _MATH_MACROS_INCLUDED
#include "../../../src_2D/Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _VECTOR2D_INCLUDED
#include "../../../src_2D/Math/Vector2D.h"
#endif //_VECTOR2D_INCLUDED

#ifndef _TENSOR2D_INCLUDED
#include "../../../src_2D/Math/Tensor2D.h"
#endif //_TENSOR2D_INCLUDED

#ifndef _GAS_CONSTANTS_INCLUDED
#include "CFD/GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED

#ifndef _HEADERDATA_INCLUDED
#include "include/HeaderData.h"
#endif

/* Define the classes. */

#define	NUM_VAR_GAUSSIAN2D    8
#define GAUSSIAN_MONATOMIC    1
#define GAUSSIAN_DIATOMIC     2


class Gaussian2D_cState;

/********************************************************
 * Class: Gaussian2D_pState                             *
 *                                                      *
 * Member functions                                     *
 *     d        -- Return density.                      *
 *     v        -- Return flow velocity.                *
 *     p        -- Return pressure.                     *
 *     atoms    -- Return number of atoms in a molecule *
 *     M        -- Return Molar mass                    *
 *     setgas   -- Set gas constants.                   *
 *     T        -- Return temperature.                  *
 *     axx      -- Return pxx/rho.                      *
 *     ayy      -- Return pyy/rho.                      *
 *     sound    -- Return sound speed.                  *
 *     dv       -- Return momentum.                     *
 *     U        -- Return conserved solution state.     *
 *     F        -- Return x-direction solution flux.    *
 *     Fx       -- Return x-direction solution flux.    *
 *     Fy       -- Return y-direction solution flux.    *
 *     Fn       -- Return n-direction solution flux.    *
 *     lambda   -- Return x-direction eigenvalue(s).    *
 *     lambda_x -- Return x-direction eigenvalue(s).    *
 *     lambda_y -- Return y-direction eigenvalue(s).    *
 *     rp       -- Return primitive right eigenvector   *
 *                 (x-direction).                       *
 *     rp_x     -- Return primitive right eigenvector   *
 *                 (x-direction).                       *
 *     rp_y     -- Return primitive right eigenvector   *
 *                 (y-direction).                       *
 *     rc       -- Return conserved right eigenvector   *
 *                 (x-direction).                       *
 *     rc_x     -- Return conserved right eigenvector   *
 *                 (x-direction).                       *
 *     rc_y     -- Return conserved right eigenvector   *
 *                 (y-direction).                       *
 *     lp       -- Return primitive left eigenvector    *
 *                 (x-direction).                       *
 *     lp_x     -- Return primitive left eigenvector    *
 *                 (x-direction).                       *
 *     lp_y     -- Return primitive left eigenvector    *
 *                 (y-direction).                       *
 *     relax    -- relaxes state towards thermodynamic  *
 *                 equilibrium (BGK source terms)       *
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
 * W = W ^ W; (my useful product)                       *
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
class Gaussian2D_pState{
  private:
  public:
    double            d;   // Density.
    Vector2D          v;   // Flow velocity (2D vector).
    Tensor2D          p;   // Pressure tensor.
    double         erot;   // rotational energy
    static double     M;   // Molar Weight
    static int    atoms;   // Monatomic or Diatomic
    static int      gas;   // Gas type
    static double alpha;   // Accommodation coefficient for gas/solid boundary
	                   // Made public so can access them.

      
    /* Creation, copy, and assignment constructors. */
    Gaussian2D_pState(void) {
       d = DENSITY_STDATM; v.zero(); p.xx = PRESSURE_STDATM;
       p.xy = ZERO; p.yy = PRESSURE_STDATM; p.zz = PRESSURE_STDATM;
       erot = PRESSURE_STDATM;
    }

    Gaussian2D_pState(const double & Val){
      d = Val; v.x = Val, v.y = Val; p.xx = Val;
      p.xy = Val; p.yy = Val; p.zz = Val;
      erot = Val;
    }

    Gaussian2D_pState(const Gaussian2D_pState &W) {
       d = W.d; v = W.v; p = W.p; erot = W.erot;
    }

    Gaussian2D_pState(const Gaussian2D_cState &U);

    Gaussian2D_pState(const double &rho,
	           const Vector2D &V,
	           const double &pre) {
       d = rho; v = V; p.xx = pre;
       p.xy = ZERO; p.yy = pre; p.zz = pre;
       erot = pre;
    }

    Gaussian2D_pState(const double &rho,
	           const double &vx,
	           const double &vy,
	           const double &pre) {
       d = rho; v.x = vx; v.y = vy; p.xx = pre;
       p.xy = ZERO; p.yy = pre; p.zz = pre;
       erot = pre;
    }

    Gaussian2D_pState(const double &rho,
		      const Vector2D &V,
		      const Tensor2D &pre,
		      const double &energyrot) {
      d = rho; v = V; p = pre; erot = energyrot;
    }

    Gaussian2D_pState(const double &rho,
		      const double &vx,
		      const double &vy,
		      const Tensor2D &pre,
		      const double &energyrot) {
      d = rho; v.x = vx; v.y = vy; p = pre; 
      erot = energyrot;
    }


    Gaussian2D_pState(const double &rho,
	           const double &vx,
	           const double &vy,
	           const double &pxx,
                   const double &pxy,
                   const double &pyy,
                   const double &pzz,
                   const double &energyrot) {
       d = rho; v.x = vx; v.y = vy; p.xx = pxx;
       p.xy = pxy; p.yy = pyy; p.zz = pzz;
       erot = energyrot;
    }

    
    /* Destructor. */
    // ~Gaussian2D_pState(void);
    // Use automatically generated destructor.

    /* Source terms */
    void relax(double deltat, int stage, const Gaussian2D_pState &W);

    /* Set gas constants. */
    void setgas(void);
    void setgas(char *string_ptr);

    /* Thermodynamic Pressure */

    double pressure(void);
    double pressure(void) const;

    /* Determinant of Pressure tensor */

    double DetP(void);
    double DetP(void) const;

    /* Validity check and adjustment */
    int invalid();
    void make_valid();

    /* Temperature. */
    double T(void);
    double T(void) const;

    /* Total energy. */
    Tensor2D E(void);
    Tensor2D E(void) const;

    /* Specific enthalpy. */
    Tensor2D h(void);
    Tensor2D h(void) const;

    /* Sound speeds. */
    double axx(void);
    double axx(void) const;

    double ayy(void);
    double ayy(void) const;

    double sound(void);
    double sound(void) const;

    /* Momentum. */
    Vector2D dv(void);
    Vector2D dv(void) const;
    double dv(const Vector2D &n);
    double dv(const Vector2D &n) const;

    /* Visocsities. */
    double viscosity(void);
    double viscosity(void) const;
    double nu(void);
    double nu(void) const;
    double bulk_viscosity(void);
    double bulk_viscosity(void) const;

    /* Mean Free Path */
    double mfp(void);

    /* Conserved solution state. */
    Gaussian2D_cState U(void);
    Gaussian2D_cState U(void) const;
    Gaussian2D_cState U(const Gaussian2D_pState &W);
    friend Gaussian2D_cState U(const Gaussian2D_pState &W);

    /* Absolute value */
    friend Gaussian2D_pState fabs(const Gaussian2D_pState &W);

    
    /* Solution flux and Jacobian (x-direction). */
    Gaussian2D_cState F(void);
    Gaussian2D_cState F(void) const;
    Gaussian2D_cState F(const Gaussian2D_pState &W);
    friend Gaussian2D_cState F(const Gaussian2D_pState &W);

    Gaussian2D_cState Fx(void);
    Gaussian2D_cState Fx(void) const;
    Gaussian2D_cState Fx(const Gaussian2D_pState &W);
    friend Gaussian2D_cState Fx(const Gaussian2D_pState &W);

    /* Solution flux and Jacobian (y-direction). */
    Gaussian2D_cState Fy(void);
    Gaussian2D_cState Fy(void) const;
    Gaussian2D_cState Fy(const Gaussian2D_pState &W);
    friend Gaussian2D_cState Fy(const Gaussian2D_pState &W);

    /* Solution flux and Jacobian (n-direction). */
    Gaussian2D_cState Fn(void);
    Gaussian2D_cState Fn(void) const;
    Gaussian2D_cState Fn(const Gaussian2D_pState &W);
    friend Gaussian2D_cState Fn(const Gaussian2D_pState &W);

    /* Eigenvalue(s) (x-direction). */
    Gaussian2D_pState lambda(void);
    Gaussian2D_pState lambda(void) const;
    Gaussian2D_pState lambda(const Gaussian2D_pState &W);
    friend Gaussian2D_pState lambda(const Gaussian2D_pState &W);
    double lambda(int index);
    double lambda(int index) const;
    friend double lambda(const Gaussian2D_pState &W, int index);

    Gaussian2D_pState lambda_x(void);
    Gaussian2D_pState lambda_x(void) const;
    Gaussian2D_pState lambda_x(const Gaussian2D_pState &W);
    friend Gaussian2D_pState lambda_x(const Gaussian2D_pState &W);
    double lambda_x(int index);
    double lambda_x(int index) const;
    friend double lambda_x(const Gaussian2D_pState &W, int index);

    /* Eigenvalue(s) (y-direction). */
    Gaussian2D_pState lambda_y(void);
    Gaussian2D_pState lambda_y(void) const;
    Gaussian2D_pState lambda_y(const Gaussian2D_pState &W);
    friend Gaussian2D_pState lambda_y(const Gaussian2D_pState &W);
    double lambda_y(int index);
    double lambda_y(int index) const;
    friend double lambda_y(const Gaussian2D_pState &W, int index);

    /* Conserved right eigenvector (x-direction). */
    Gaussian2D_cState rc(int index);
    Gaussian2D_cState rc(int index) const;
    friend Gaussian2D_cState rc(const Gaussian2D_pState &W, int index);

    Gaussian2D_cState rc_x(int index);
    Gaussian2D_cState rc_x(int index) const;
    friend Gaussian2D_cState rc_x(const Gaussian2D_pState &W, int index);

   /* Conserved right eigenvector (y-direction). */
    Gaussian2D_cState rc_y(int index);
    Gaussian2D_cState rc_y(int index) const;
    friend Gaussian2D_cState rc_y(const Gaussian2D_pState &W, int index);

    /* Primitive left eigenvector (x-direction). */
    Gaussian2D_pState lp(int index);
    Gaussian2D_pState lp(int index) const;
    friend Gaussian2D_pState lp(const Gaussian2D_pState &W, int index);

    Gaussian2D_pState lp_x(int index);
    Gaussian2D_pState lp_x(int index) const;
    friend Gaussian2D_pState lp_x(const Gaussian2D_pState &W, int index);

    /* Primitive left eigenvector (y-direction). */
    Gaussian2D_pState lp_y(int index);
    Gaussian2D_pState lp_y(int index) const;
    friend Gaussian2D_pState lp_y(const Gaussian2D_pState &W, int index);

    /* Index operator. */
    double &operator[](int index) {
      assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
      switch(index) {
        case 1 :
	  return (d);
        case 2 :
	  return (v.x);
        case 3 :
	  return (v.y);
        case 4 :
	  return (p.xx);
        case 5 :
	  return (p.xy);
        case 6 :
	  return (p.yy);
        case 7 :
	  return (p.zz);
        case 8 :
	  return (erot);
        default:
	  return (d);
      };
    }
    
    const double &operator[](int index) const {
      assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
      switch(index) {
        case 1 :
	  return (d);
        case 2 :
	  return (v.x);
        case 3 :
	  return (v.y);
        case 4 :
	  return (p.xx);
        case 5 :
	  return (p.xy);
        case 6 :
	  return (p.yy);
        case 7 :
	  return (p.zz);
        case 8 :
	  return (erot);
        default:
	  return (d);
      };
    }

    /* Binary arithmetic operators. */
    friend Gaussian2D_pState operator +(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2);
    friend Gaussian2D_pState operator -(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2);
    friend double operator *(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2);
    friend Gaussian2D_pState operator *(const Gaussian2D_pState &W, const double &a);
    friend Gaussian2D_pState operator *(const double &a, const Gaussian2D_pState &W);
    friend Gaussian2D_pState operator /(const Gaussian2D_pState &W, const double &a);
    friend Gaussian2D_pState operator /(const Gaussian2D_pState &up, const Gaussian2D_pState &down );
    friend Gaussian2D_pState operator ^(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2);
    
    /* Unary arithmetic operators. */
    friend Gaussian2D_pState operator +(const Gaussian2D_pState &W);
    friend Gaussian2D_pState operator -(const Gaussian2D_pState &W);

    /* Shortcut arithmetic operators. */
    friend Gaussian2D_pState &operator +=(Gaussian2D_pState &W1, const Gaussian2D_pState &W2);
    friend Gaussian2D_pState &operator -=(Gaussian2D_pState &W1, const Gaussian2D_pState &W2);
    
    /* Relational operators. */
    friend int operator ==(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2);
    friend int operator !=(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const Gaussian2D_pState &W);
    friend istream &operator >> (istream &in_file,  Gaussian2D_pState &W);

    /* Comparison operators */
    friend bool operator < (const Gaussian2D_pState &W1, const Gaussian2D_pState &W2);
};

/********************************************************
 * Class: Gaussian2D_cState                             *
 *                                                      *
 * Member functions                                     *
 *      d       -- Return density.                      *
 *      dv      -- Return momentum.                     *
 *      E       -- Return total energy.                 *
 *      atoms   -- Return number of atoms in a molecule *
 *      M       -- Return molar mass                    *
 *      setgas  -- Set gas constants.                   *
 *      v       -- Return flow velocity.                *
 *      p       -- Return pressure.                     *
 *      W       -- Return primitive solution state.     *
 *      F       -- Return x-direction solution flux.    *
 *      Fx      -- Return x-direction solution flux.    *
 *      Fy      -- Return y-direction solution flux.    *
 *      Fn      -- Return n-direction solution flux.    *
 *                                                      *
 * Member operators                                     *
 *      U -- a conserved solution state                 *
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
 * U = U ^ U; (my useful product)                       *
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
class Gaussian2D_cState{
  private:
  public:
    double            d;   // Density.
    Vector2D         dv;   // Momentum.
    Tensor2D          E;   // Total Energy.
    double         erot;   // Rotational Energy
    static double     M;   // Molar Weight
    static int    atoms;   // Monatomic or Diatomic
    static int      gas;   // gas type
    static double alpha;   // Accommodation coefficient for gas/solid boundary
	                   // Made public so can access them.
		      
    /* Creation, copy, and assignment constructors. */
    Gaussian2D_cState(void) {
       d = DENSITY_STDATM; dv.zero(); E.xx = PRESSURE_STDATM;
       E.xy = ZERO; E.yy = PRESSURE_STDATM; E.zz = PRESSURE_STDATM;
       erot = PRESSURE_STDATM;
    }

    Gaussian2D_cState(double & Val){
      d = Val; dv.x = Val, dv.y = Val; E.xx = Val;
      E.xy = Val; E.yy = Val; E.zz = Val;
      erot = Val;
    }

    Gaussian2D_cState(const Gaussian2D_cState &U) {
       d = U.d; dv = U.dv; E = U.E; erot = U.erot;
    }

    Gaussian2D_cState(const Gaussian2D_pState &W);

    Gaussian2D_cState(const double &rho,
	           const Vector2D &rhoV,
	           const Tensor2D &Etotal,
                   const double &energyrot) {
       d = rho; dv = rhoV; E = Etotal;
       erot = energyrot;
    }

    Gaussian2D_cState(const double &rho,
	           const double &rhovx,
	           const double &rhovy,
	           const double &exx,
                   const double &exy,
                   const double &eyy,
                   const double &ezz,
                   const double &energyrot) {
       d = rho; dv.x = rhovx; dv.y = rhovy; E.xx = exx;
       E.xy = exy, E.yy = eyy; E.zz = ezz; erot = energyrot;
    }
    
    /* Destructor. */
    // ~Gaussian2D_cState(void);
    // Use automatically generated destructor.

    /* Set gas constants. */
    void setgas(void);
    void setgas(char *string_ptr);

    /* Flow velocity. */
    Vector2D v(void);
    Vector2D v(void) const;
    double v(const Vector2D &n);
    double v(const Vector2D &n) const;

    /* Pressure. */
    Tensor2D p(void);
    Tensor2D p(void) const;

    /* Validity check */
    int invalid();

    /* Primitive solution state. */
    Gaussian2D_pState W(void);
    Gaussian2D_pState W(void) const;
    Gaussian2D_pState W(const Gaussian2D_cState &U);
    friend Gaussian2D_pState W(const Gaussian2D_cState &U);
    
    /* Solution flux and Jacobian (x-direction). */
    Gaussian2D_cState F(void);
    Gaussian2D_cState F(void) const;
    Gaussian2D_cState F(const Gaussian2D_cState &U);
    friend Gaussian2D_cState F(const Gaussian2D_cState &U);

    Gaussian2D_cState Fx(void);
    Gaussian2D_cState Fx(void) const;
    Gaussian2D_cState Fx(const Gaussian2D_cState &U);
    friend Gaussian2D_cState Fx(const Gaussian2D_cState &U);

    /* Solution flux and Jacobian (y-direction). */
    Gaussian2D_cState Fy(void);
    Gaussian2D_cState Fy(void) const;
    Gaussian2D_cState Fy(const Gaussian2D_cState &U);
    friend Gaussian2D_cState Fy(const Gaussian2D_cState &U);

    /* Solution flux and Jacobian (n-direction). */
    Gaussian2D_cState Fn(void);
    Gaussian2D_cState Fn(void) const;
    Gaussian2D_cState Fn(const Gaussian2D_cState &U);
    friend Gaussian2D_cState Fn(const Gaussian2D_cState &U);

    /* Assignment operator. */
    // Gaussian2D_cState operator = (const Gaussian2D_cState &U);
    // Use automatically generated assignment operator.

    /* Index operator. */
    double &operator[](int index) {
      assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
      switch(index) {
        case 1 :
	  return (d);
        case 2 :
	  return (dv.x);
        case 3 :
	  return (dv.y);
        case 4 :
	  return (E.xx);
        case 5 :
	  return (E.xy);
        case 6 :
	  return (E.yy);
        case 7 :
	  return (E.zz);
        case 8 :
	  return (erot);
        default:
	  return (d);
      };
    }
    
    const double &operator[](int index) const {
      assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
      switch(index) {
        case 1 :
	  return (d);
        case 2 :
	  return (dv.x);
        case 3 :
	  return (dv.y);
        case 4 :
	  return (E.xx);
        case 5 :
	  return (E.xy);
        case 6 :
	  return (E.yy);
        case 7 :
	  return (E.zz);
        case 8 :
	  return (erot);
        default:
	  return (d);
      };
    }

    /* Binary arithmetic operators. */
    friend Gaussian2D_cState operator +(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2);
    friend Gaussian2D_cState operator -(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2);
    friend double operator *(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2);
    friend Gaussian2D_cState operator *(const Gaussian2D_cState &U, const double &a);
    friend Gaussian2D_cState operator *(const double &a, const Gaussian2D_cState &U);
    friend Gaussian2D_cState operator /(const Gaussian2D_cState &U, const double &a);
    friend Gaussian2D_cState operator ^(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2);

    /* Unary arithmetic operators. */
    friend Gaussian2D_cState operator +(const Gaussian2D_cState &U);
    friend Gaussian2D_cState operator -(const Gaussian2D_cState &U);

    /* Shortcut arithmetic operators. */
    friend Gaussian2D_cState &operator +=(Gaussian2D_cState &U1, const Gaussian2D_cState &U2);
    friend Gaussian2D_cState &operator -=(Gaussian2D_cState &U1, const Gaussian2D_cState &U2);
    
    /* Relational operators. */
    friend int operator ==(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2);
    friend int operator !=(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const Gaussian2D_cState &U);
    friend istream &operator >> (istream &in_file,  Gaussian2D_cState &U);
    
};

/***********************************************************************
 * Gaussian2D_pState::fabs -- Return the absolute value of the object  *
 ***********************************************************************/
inline Gaussian2D_pState fabs (const Gaussian2D_pState &W)
{
  Gaussian2D_pState temp;
  for (int i=1; i<=NUM_VAR_GAUSSIAN2D; i++)
    temp[i] = W[i];

  return temp;
}

/********************************************************
 * Gaussian2D_pState::setgas -- Assign gas constants.   *
 ********************************************************/
inline void Gaussian2D_pState::setgas(void) {
  M = MOLE_WT_AIR;
  atoms = GAUSSIAN_DIATOMIC;
  gas = GAS_AIR;
}

inline void Gaussian2D_pState::setgas(char *string_ptr) {
   if (strcmp(string_ptr, "AIR") == 0) {
     M = MOLE_WT_AIR;
     atoms = GAUSSIAN_DIATOMIC;
     gas = GAS_AIR;
   } else if (strcmp(string_ptr, "A") == 0) {
     M = MOLE_WT_A;
     atoms = GAUSSIAN_MONATOMIC;
     gas = GAS_A;
   } else if (strcmp(string_ptr, "CO") == 0) {
     M = MOLE_WT_CO;
     atoms = GAUSSIAN_DIATOMIC;
     gas = GAS_CO;
   } else if (strcmp(string_ptr, "H2") == 0) {
     M = MOLE_WT_H2;
     atoms = GAUSSIAN_DIATOMIC;
     gas = GAS_H2;
   } else if (strcmp(string_ptr, "HE") == 0) {
     M = MOLE_WT_HE;
     atoms = GAUSSIAN_MONATOMIC;
     gas = GAS_HE;
   } else if (strcmp(string_ptr, "N2") == 0) {
     M = MOLE_WT_N2;
     atoms = GAUSSIAN_DIATOMIC;
     gas = GAS_N2;
   } else if (strcmp(string_ptr, "O2") == 0) {
     M = MOLE_WT_O2;
     atoms = GAUSSIAN_DIATOMIC;
     gas = GAS_O2;
   } else {
     M = MOLE_WT_AIR;
     atoms = GAUSSIAN_DIATOMIC;
     gas = GAS_AIR;
   } /* endif */
   if(atoms == GAUSSIAN_MONATOMIC){
     erot = 0.0;
   }
}

/*********************************************************
 * Gaussian2D_pState::pressure -- Thermodynamic pressure.*
 *********************************************************/

inline double Gaussian2D_pState::pressure(void) {
  assert( atoms == GAUSSIAN_MONATOMIC || atoms == GAUSSIAN_DIATOMIC );
  if(atoms==GAUSSIAN_MONATOMIC){
    return (p.xx+p.yy+p.zz)/3.0;
  }else{
    return (p.xx+p.yy+p.zz+2.0*erot)/5.0;
  }
}

inline double Gaussian2D_pState::pressure(void) const {
  assert( atoms == GAUSSIAN_MONATOMIC || atoms == GAUSSIAN_DIATOMIC );
  if(atoms==GAUSSIAN_MONATOMIC){
    return (p.xx+p.yy+p.zz)/3.0;
  }else{
    return (p.xx+p.yy+p.zz+2.0*erot)/5.0;
  }
}

/*********************************************************
 * Gaussian2D_pState::DetP -- Determinant of Press tensor*
 *********************************************************/

inline double Gaussian2D_pState::DetP(void) {
  return (p.xx*p.yy-sqr(p.xy));
}

inline double Gaussian2D_pState::DetP(void) const {
  return (p.xx*p.yy-sqr(p.xy));
}

/************************************************************
 * Gaussian2D_pState::invalid, checks for physical validity *
 ************************************************************/

inline int Gaussian2D_pState::invalid() {

  int check(0);
  double det = DetP();

  //ensure density > 0
  //and all relevant pressures > 0

  if(   d <= 0.0   ||
        det <= 0.0 ||
	p.xx <= 0.0 ||
	p.yy <= 0.0 ||
        p.zz <= 0.0 ){
    check = 1;
  }

  if((atoms == GAUSSIAN_DIATOMIC) &&
     (erot <= 0.0))
    { check = 1; }

  return check;
}

inline void Gaussian2D_pState::make_valid() {

  if(p.xx<0.0){p.xx = TOLER;}
  if(p.yy<0.0){p.yy = TOLER;}

  if(DetP() < TOLER/MILLION*(p.xx*p.yy)) {
    p.xy = sgn(p.xy)*(1.0-TOLER)*sqrt(p.xx*p.yy);
  }
  return;
}


/********************************************************
 * Gaussian2D_pState::T -- Temperature.                 *
 ********************************************************/

inline double Gaussian2D_pState::T(void) {
  if(pressure()<=0.0) {
    cout << "P=" << pressure();
    return(THOUSAND);
  }
  if(d<=0.0) {
    cout << "d=" << d;
    return (THOUSAND);
  }
    return (pressure()/(d*AVOGADRO*THOUSAND/M*BOLTZMANN));
}

inline double Gaussian2D_pState::T(void) const {
  if(pressure()<=0.0) {
    cout << "P=" << pressure();
    return(THOUSAND);
  }
  if(d<=0.0) {
    cout << "d=" << d;
    return (THOUSAND);
  }
    return (pressure()/(d*AVOGADRO*THOUSAND/M*BOLTZMANN));
}

/********************************************************
 * Euler2D_pState::E -- Total energy.                   *
 ********************************************************/

inline Tensor2D Gaussian2D_pState::E(void) {

    Tensor2D temp;

    temp.xx = p.xx + d*v.x*v.x;
    temp.xy = p.xy + d*v.x*v.y;
    temp.yy = p.yy + d*v.y*v.y;
    temp.zz = p.zz;

    return (temp); 
}

inline Tensor2D Gaussian2D_pState::E(void) const {

    Tensor2D temp;

    temp.xx = p.xx + d*v.x*v.x;
    temp.xy = p.xy + d*v.x*v.y;
    temp.yy = p.yy + d*v.y*v.y;
    temp.zz = p.zz;

    return (temp); 

}

/********************************************************
 * Euler2D_pState::h -- Specific enthalpy.              *
 ********************************************************/

inline Tensor2D Gaussian2D_pState::h(void) {

    double hxx, hxy, hyy, hzz;

    hxx = 1.5*p.xx/d+0.5*v.x*v.x;
    hxy = 1.5*p.xy/d+0.5*v.x*v.y;
    hyy = 1.5*p.yy/d+0.5*v.y*v.y;
    hzz = 1.5*p.zz/d;

    return (Tensor2D(hxx,hxy,hyy,hzz));
}

inline Tensor2D Gaussian2D_pState::h(void) const {

    double hxx, hxy, hyy, hzz;

    hxx = 1.5*p.xx/d+0.5*v.x*v.x;
    hxy = 1.5*p.xy/d+0.5*v.x*v.y;
    hyy = 1.5*p.yy/d+0.5*v.y*v.y;
    hzz = 1.5*p.zz/d;

    return (Tensor2D(hxx,hxy,hyy,hzz));
}

/********************************************************
 * Euler2D_pState::axx & ayy -- Sound speeds. (sort of) *
 ********************************************************/

inline double Gaussian2D_pState::axx(void) {
  assert( p.xx>0.0 && d>0 );
  double a(p.xx/d);
  //double c(5.0);
  //if(a>5){
  return (sqrt(a));
    //} else {
    //return (sqrt(c)/exp(0.5)*exp(a/(2.0*c)));
    //}
}

inline double Gaussian2D_pState::axx(void) const {
  assert( p.xx>0.0 && d>0 );
  double a(p.xx/d);
  //  double c(5.0);
  //if(a>5){
  return (sqrt(a));
    //} else {
    //return (sqrt(c)/exp(0.5)*exp(a/(2.0*c)));
    //}
}

inline double Gaussian2D_pState::ayy(void) {
    assert( p.yy>0.0 && d>0 );
    return (sqrt(p.yy/d));
}

inline double Gaussian2D_pState::ayy(void) const {
    assert( p.yy>0.0 && d>0 );
    return (sqrt(p.yy/d));
}

inline double Gaussian2D_pState::sound(void) {
    if(atoms == 1){
      return (sqrt((5.0/3.0)*(p.xx+p.yy+p.zz)/3.0/d));
    }else if (atoms == 2){
      return (sqrt(1.4*(p.xx+p.yy+p.zz+2.0*erot)/5.0/d));
    }else{
      cout << "Error....atoms value of " << atoms << "is not allowed.\n";
      cout << "assuming atoms = 1\n";
      return (sqrt((5.0/3.0)*(p.xx+p.yy+p.zz)/3/d));
    }
}

inline double Gaussian2D_pState::sound(void) const {
    if(atoms == 1){
      return (sqrt((5.0/3.0)*(p.xx+p.yy+p.zz)/3.0/d));
    }else if (atoms == 2){
      return (sqrt(1.4*(p.xx+p.yy+p.zz+2.0*erot)/5.0/d));
    }else{
      cout << "Error....atoms value of " << atoms << "is not allowed.\n";
      cout << "assuming atoms = 1\n";
      return (sqrt((5.0/3.0)*(p.xx+p.yy+p.zz)/3/d));
    }
}

/********************************************************
 * Gaussian2D_pState::dv -- Momentum.                   *
 ********************************************************/
inline Vector2D Gaussian2D_pState::dv(void) {
    return (d*v);
}

inline Vector2D Gaussian2D_pState::dv(void) const {
    return (d*v);
}

inline double Gaussian2D_pState::dv(const Vector2D &n) {
    return (d*(v*n));
}

inline double Gaussian2D_pState::dv(const Vector2D &n) const {
    return (d*(v*n));
}

/********************************************************
 * Gaussian2D_pState::dv -- Viscosity.                  *
 ********************************************************/
inline double Gaussian2D_pState::viscosity(void){
  
  double omega, mu_not;

  switch(gas) {
  case GAS_AIR:
    omega = OMEGA_AIR;
    mu_not = MU_NOT_AIR;
    break;
  case GAS_A:
    omega = OMEGA_A;
    mu_not = MU_NOT_A;
    break;
  case GAS_CO:
    omega = OMEGA_CO;
    mu_not = MU_NOT_CO;
    break;
  case GAS_H2:
    omega = OMEGA_H2;
    mu_not = MU_NOT_H2;
    break;
  case GAS_HE:
    omega = OMEGA_HE;
    mu_not = MU_NOT_HE;
    break;
  case GAS_N2:
    omega = OMEGA_N2;
    mu_not = MU_NOT_N2;
    break;
  case GAS_O2:
    omega = OMEGA_O2;
    mu_not = MU_NOT_O2;
    break;
  default:
    omega = OMEGA_AIR;
    mu_not = MU_NOT_AIR;
    break;
  }  

  return (mu_not*pow((T()/273.0),omega));

}

inline double Gaussian2D_pState::viscosity(void) const{
  
  double omega, mu_not;

  switch(gas) {
  case GAS_AIR:
    omega = OMEGA_AIR;
    mu_not = MU_NOT_AIR;
    break;
  case GAS_A:
    omega = OMEGA_A;
    mu_not = MU_NOT_A;
    break;
  case GAS_CO:
    omega = OMEGA_CO;
    mu_not = MU_NOT_CO;
    break;
  case GAS_H2:
    omega = OMEGA_H2;
    mu_not = MU_NOT_H2;
    break;
  case GAS_HE:
    omega = OMEGA_HE;
    mu_not = MU_NOT_HE;
    break;
  case GAS_N2:
    omega = OMEGA_N2;
    mu_not = MU_NOT_N2;
    break;
  case GAS_O2:
    omega = OMEGA_O2;
    mu_not = MU_NOT_O2;
    break;
  default:
    omega = OMEGA_AIR;
    mu_not = MU_NOT_AIR;
    break;
  }  

  return (mu_not*pow((T()/273.0),omega));

}

inline double Gaussian2D_pState::nu(void){
  return (viscosity()/d);
}

inline double Gaussian2D_pState::nu(void) const{
  return (viscosity()/d);
}

inline double Gaussian2D_pState::bulk_viscosity(void){
  return (3.0*viscosity());
}

inline double Gaussian2D_pState::bulk_viscosity(void) const{
  return (3.0*viscosity());
}

/********************************************************
 * Gaussian2D_pState -- Mean Free Path (hard spheres)   *
 ********************************************************/
inline double Gaussian2D_pState::mfp(void){
  return (16.0*viscosity()/(5.0*sqrt(2.0*PI*d*pressure())));
}


/********************************************************
 * Gaussian2D_pState -- Binary arithmetic operators.    *
 ********************************************************/
inline Gaussian2D_pState operator +(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2) {
  return (Gaussian2D_pState(W1.d+W2.d,W1.v+W2.v,W1.p+W2.p,W1.erot+W2.erot));
}

inline Gaussian2D_pState operator -(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2) {
  return (Gaussian2D_pState(W1.d-W2.d,W1.v-W2.v,W1.p-W2.p,W1.erot-W2.erot));
}

// Inner product operator.
inline double operator *(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2) {
   return (W1.d*W2.d+W1.v*W2.v+W1.p.xx*W2.p.xx+W1.p.xy*W2.p.xy+W1.p.yy*W2.p.yy+
           W1.p.zz*W2.p.zz+W1.erot*W2.erot);
}

inline Gaussian2D_pState operator *(const Gaussian2D_pState &W, const double &a) {
  return (Gaussian2D_pState(a*W.d,a*W.v,a*W.p,a*W.erot));
}

inline Gaussian2D_pState operator *(const double &a, const Gaussian2D_pState &W) {
  return (Gaussian2D_pState(a*W.d,a*W.v,a*W.p,a*W.erot));
}

inline Gaussian2D_pState operator /(const Gaussian2D_pState &W, const double &a) {
  return (Gaussian2D_pState(W.d/a,W.v/a,W.p/a,W.erot/a));
}

inline Gaussian2D_pState operator /(const Gaussian2D_pState &up, const Gaussian2D_pState &down) {

  if((down[1]==0)||(down[2]==0)||(down[3]==0)||(down[4]==0)||(down[5]==0)
     ||(down[6]==0)||(down[7]==0)||(down[8]==0))
    return Gaussian2D_pState(0.0);
  
  return (Gaussian2D_pState(up.d/down.d, up.v.x/down.v.x, up.v.y/down.v.y, up.p.xx/down.p.xx,
			    up.p.xy/down.p.xy, up.p.yy/down.p.yy, up.p.zz/down.p.zz, up.erot/down.erot));
}

inline bool operator < (const Gaussian2D_pState &W1, const Gaussian2D_pState &W2){

  bool result = true;
  for (int i=1; i<=NUM_VAR_GAUSSIAN2D; ++i)
    result = result && (W1[i]<W2[i]);
  return result;
}

// My useful solution state product operator.
inline Gaussian2D_pState operator ^(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2) {
   return (Gaussian2D_pState(W1.d*W2.d,W1.v.x*W2.v.x,W1.v.y*W2.v.y,W1.p.xx*W2.p.xx,
                          W1.p.xy*W2.p.xy,W1.p.yy*W2.p.yy,W1.p.zz*W2.p.zz,W1.erot*W2.erot));
}

/********************************************************
 * Gaussian2D_pState -- Unary arithmetic operators.     *
 ********************************************************/
inline Gaussian2D_pState operator +(const Gaussian2D_pState &W) {
  return (Gaussian2D_pState(W.d,W.v,W.p,W.erot));
}

inline Gaussian2D_pState operator -(const Gaussian2D_pState &W) {
  return (Gaussian2D_pState(-W.d,-W.v,-W.p,-W.erot));
}

/********************************************************
 * Gaussian2D_pState -- Shortcut arithmetic operators.  *
 ********************************************************/
inline Gaussian2D_pState &operator +=(Gaussian2D_pState &W1, const Gaussian2D_pState &W2) {
  W1.d    += W2.d;
  W1.v    += W2.v;
  W1.p    += W2.p;
  W1.erot += W2.erot;
  return (W1);
}

inline Gaussian2D_pState &operator -=(Gaussian2D_pState &W1, const Gaussian2D_pState &W2) {
  W1.d    -= W2.d;
  W1.v    -= W2.v;
  W1.p    -= W2.p;
  W1.erot -= W2.erot;
  return (W1);
}

/********************************************************
 * Gaussian2D_pState -- Relational operators.           *
 ********************************************************/
inline int operator ==(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2) {
  return (W1.d == W2.d && W1.v == W2.v && W1.p == W2.p && W1.erot == W2.erot);
}

inline int operator !=(const Gaussian2D_pState &W1, const Gaussian2D_pState &W2) {
  return (W1.d != W2.d || W1.v != W2.v || W1.p != W2.p || W1.erot != W2.erot);
}

/********************************************************
 * Gaussian2D_pState -- Input-output operators.         *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const Gaussian2D_pState &W) {
  out_file.setf(ios::scientific);
  out_file << " " << W.d  << " " << W.v.x << " " << W.v.y << " " << W.p.xx
           << " " << W.p.xy << " " << W.p.yy << " " << W.p.zz << " " << W.erot;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, Gaussian2D_pState &W) {
  in_file.setf(ios::skipws);
  in_file >> W.d >> W.v.x >> W.v.y >> W.p.xx >> W.p.xy
          >> W.p.yy >> W.p.zz >> W.erot;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/********************************************************
 * Gaussian2D_cState::setgas -- Assign gas constants.   *
 ********************************************************/
inline void Gaussian2D_cState::setgas(void) {
  M = MOLE_WT_AIR;
  atoms = GAUSSIAN_DIATOMIC;
  gas = GAS_AIR;
}

inline void Gaussian2D_cState::setgas(char *string_ptr) {
   if (strcmp(string_ptr, "AIR") == 0) {
     M = MOLE_WT_AIR;
     atoms = GAUSSIAN_DIATOMIC;
     gas = GAS_AIR;
    } else if (strcmp(string_ptr, "A") == 0) {
     M = MOLE_WT_A;
     atoms = GAUSSIAN_MONATOMIC;
     gas = GAS_A;
   } else if (strcmp(string_ptr, "CO") == 0) {
     M = MOLE_WT_CO;
     atoms = GAUSSIAN_DIATOMIC;
     gas = GAS_CO;
   } else if (strcmp(string_ptr, "H2") == 0) {
     M = MOLE_WT_H2;
     atoms = GAUSSIAN_DIATOMIC;
     gas = GAS_H2;
   } else if (strcmp(string_ptr, "HE") == 0) {
     M = MOLE_WT_HE;
     atoms = GAUSSIAN_MONATOMIC;
     gas = GAS_HE;
   } else if (strcmp(string_ptr, "N2") == 0) {
     M = MOLE_WT_N2;
     atoms = GAUSSIAN_DIATOMIC;
     gas = GAS_N2;
   } else if (strcmp(string_ptr, "O2") == 0) {
     M = MOLE_WT_O2;
     atoms = GAUSSIAN_DIATOMIC;
     gas = GAS_O2;
   } else {
     M = MOLE_WT_AIR;
     atoms = GAUSSIAN_DIATOMIC;
     gas = GAS_AIR;
   } /* endif */
   if(atoms == GAUSSIAN_MONATOMIC){
     erot = 0.0;
   }
}

/********************************************************
 * Gaussian2D_cState::v -- Flow velocity.               *
 ********************************************************/
inline Vector2D Gaussian2D_cState::v(void) {
    return (dv/d);
}

inline Vector2D Gaussian2D_cState::v(void) const {
    return (dv/d);
}

inline double Gaussian2D_cState::v(const Vector2D &n) {
    return ((dv*n)/d);
}

inline double Gaussian2D_cState::v(const Vector2D &n) const {
    return ((dv*n)/d);
}

/********************************************************
 * Gaussian2D_cState::p -- Pressure.                    *
 ********************************************************/

inline Tensor2D Gaussian2D_cState::p(void) {

    Tensor2D temp;

    temp.xx = E.xx - dv.x*dv.x/d;
    temp.xy = E.xy - dv.x*dv.y/d;
    temp.yy = E.yy - dv.y*dv.y/d;
    temp.zz = E.zz;

    return (temp); 

}

inline Tensor2D Gaussian2D_cState::p(void) const {

    Tensor2D temp;

    temp.xx = E.xx - dv.x*dv.x/d;
    temp.xy = E.xy - dv.x*dv.y/d;
    temp.yy = E.yy - dv.y*dv.y/d;
    temp.zz = E.zz;

    return (temp); 

}

/************************************************************
 * Gaussian2D_cState::invalid, checks for physical validity *
 ************************************************************/

inline int Gaussian2D_cState::invalid()
{
  int check(0);
  double det;
  Tensor2D pressure;

  pressure = p();

  //find determinant of P tensor

  det = (pressure.xx*pressure.yy-sqr(pressure.xy));

  //ensure density > 0
  //and all relevant pressures > 0

  if(   d <= 0.0   ||
        det <= 0.0 ||
	pressure.xx <= 0.0 ||
	pressure.yy <= 0.0 ||
        pressure.zz <= 0.0 ){
    check = 1;
  }

  if((atoms == GAUSSIAN_DIATOMIC) &&
     (erot <= 0.0))
    { check = 1; }

  return check;
}

/********************************************************
 * Gaussian2D_cState -- Binary arithmetic operators.       *
 ********************************************************/

inline Gaussian2D_cState operator +(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2) {
  return (Gaussian2D_cState(U1.d+U2.d,U1.dv+U2.dv,U1.E+U2.E,U1.erot+U2.erot));
}

inline Gaussian2D_cState operator -(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2) {
  return (Gaussian2D_cState(U1.d-U2.d,U1.dv-U2.dv,U1.E-U2.E,U1.erot-U2.erot));
}

// Inner product operator.
inline double operator *(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2) {
   return (U1.d*U2.d+U1.dv*U2.dv+U1.E.xx*U2.E.xx+U1.E.xy*U2.E.xy+U1.E.yy*U2.E.yy
           +U1.E.zz*U2.E.zz+U1.erot*U2.erot);
}

inline Gaussian2D_cState operator *(const Gaussian2D_cState &U, const double &a) {
  return (Gaussian2D_cState(a*U.d,a*U.dv,a*U.E,a*U.erot));
}

inline Gaussian2D_cState operator *(const double &a, const Gaussian2D_cState &U) {
  return (Gaussian2D_cState(a*U.d,a*U.dv,a*U.E,a*U.erot));
}

inline Gaussian2D_cState operator /(const Gaussian2D_cState &U, const double &a) {
  return (Gaussian2D_cState(U.d/a,U.dv/a,U.E/a,U.erot/a));
}

// My useful solution state product operator.
inline Gaussian2D_cState operator ^(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2) {
  return (Gaussian2D_cState(U1.d*U2.d,U1.dv.x*U2.dv.x,U1.dv.y*U2.dv.y,U1.E.xx*U2.E.xx,
                             U1.E.xy*U2.E.xy,U1.E.yy*U2.E.yy,U1.E.zz*U2.E.zz,U1.erot*U2.erot));
}

/********************************************************
 * Gaussian2D_cState -- Unary arithmetic operators.        *
 ********************************************************/
inline Gaussian2D_cState operator +(const Gaussian2D_cState &U) {
  return (Gaussian2D_cState(U.d,U.dv,U.E,U.erot));
}

inline Gaussian2D_cState operator -(const Gaussian2D_cState &U) {
  return (Gaussian2D_cState(-U.d,-U.dv,-U.E,-U.erot));
}

/********************************************************
 * Gaussian2D_cState -- Shortcut arithmetic operators.     *
 ********************************************************/
inline Gaussian2D_cState &operator +=(Gaussian2D_cState &U1, const Gaussian2D_cState &U2) {
  U1.d    += U2.d;
  U1.dv   += U2.dv;
  U1.E    += U2.E;
  U1.erot += U2.erot;
  return (U1);
}

inline Gaussian2D_cState &operator -=(Gaussian2D_cState &U1, const Gaussian2D_cState &U2) {
  U1.d    -= U2.d;
  U1.dv   -= U2.dv;
  U1.E    -= U2.E;
  U1.erot -= U2.erot;
  return (U1);
}

/********************************************************
 * Gaussian2D_cState -- Relational operators.              *
 ********************************************************/
inline int operator ==(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2) {
  return (U1.d == U2.d && U1.dv == U2.dv && U1.E == U2.E && U1.erot == U2.erot);
}

inline int operator !=(const Gaussian2D_cState &U1, const Gaussian2D_cState &U2) {
  return (U1.d != U2.d || U1.dv != U2.dv || U1.E != U2.E || U1.erot != U2.erot);
}

/********************************************************
 * Gaussian2D_cState -- Input-output operators.            *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const Gaussian2D_cState &U) {
  out_file.setf(ios::scientific);
  out_file << " " << U.d  << " " << U.dv.x << " " << U.dv.y << " " << U.E.xx
           << " " << U.E.xy << " " << U.E.yy << " " << U.E.zz << " " << U.erot;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, Gaussian2D_cState &U) {
  in_file.setf(ios::skipws);
  in_file >> U.d >> U.dv.x >> U.dv.y >> U.E.xx >> U.E.xy 
          >> U.E.yy >> U.E.zz >> U.erot;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/********************************************************
 * Gaussian2D_pState::Euler2D_pState -- Constructor.    *
 ********************************************************/

inline Gaussian2D_pState::Gaussian2D_pState(const Gaussian2D_cState &U) {
  d = U.d; v = U.v(); p = U.p(); erot = U.erot;
}

/********************************************************
 * Gaussian2D_pState::U -- Conserved solution state.    *
 ********************************************************/

inline Gaussian2D_cState Gaussian2D_pState::U(void) {
  return (Gaussian2D_cState(d, dv(), E(), erot));
}

inline Gaussian2D_cState Gaussian2D_pState::U(void) const {
  return (Gaussian2D_cState(d, dv(), E(), erot));
}

inline Gaussian2D_cState Gaussian2D_pState::U(const Gaussian2D_pState &W) {
  return (Gaussian2D_cState(W.d, W.dv(), W.E(), W.erot));
}

inline Gaussian2D_cState U(const Gaussian2D_pState &W) {
  return (Gaussian2D_cState(W.d, W.dv(), W.E(), W.erot));
}

/********************************************************
 * Gaussian2D_pState::F -- Solution flux (x-direction). *
 ********************************************************/

inline Gaussian2D_cState Gaussian2D_pState::F(void) {
  return (Gaussian2D_cState(d*v.x, 
                            d*v.x*v.x + p.xx, 
                            d*v.x*v.y + p.xy, 
                            d*v.x*v.x*v.x + 3.0*v.x*p.xx,
                            d*v.x*v.x*v.y + 2.0*v.x*p.xy + v.y*p.xx,
                            d*v.x*v.y*v.y + v.x*p.yy + 2.0*v.y*p.xy,
                            v.x*p.zz,
                            v.x*erot));
}

inline Gaussian2D_cState Gaussian2D_pState::F(void) const {
  return (Gaussian2D_cState(d*v.x, 
                            d*v.x*v.x + p.xx, 
                            d*v.x*v.y + p.xy, 
                            d*v.x*v.x*v.x + 3.0*v.x*p.xx,
                            d*v.x*v.x*v.y + 2.0*v.x*p.xy + v.y*p.xx,
                            d*v.x*v.y*v.y + v.x*p.yy + 2.0*v.y*p.xy,
                            v.x*p.zz,
                            v.x*erot));
}

inline Gaussian2D_cState Gaussian2D_pState::F(const Gaussian2D_pState &W) {
  return (Gaussian2D_cState(W.d*W.v.x, 
                            W.d*W.v.x*W.v.x + W.p.xx, 
                            W.d*W.v.x*W.v.y + W.p.xy, 
                            W.d*W.v.x*W.v.x*W.v.x + 3.0*W.v.x*W.p.xx,
                            W.d*W.v.x*W.v.x*W.v.y + 2.0*W.v.x*W.p.xy + W.v.y*W.p.xx,
                            W.d*W.v.x*W.v.y*W.v.y + W.v.x*W.p.yy + 2.0*W.v.y*W.p.xy,
                            W.v.x*W.p.zz,
                            W.v.x*W.erot));
}

inline Gaussian2D_cState F(const Gaussian2D_pState &W) {
  return (Gaussian2D_cState(W.d*W.v.x, 
                            W.d*W.v.x*W.v.x + W.p.xx, 
                            W.d*W.v.x*W.v.y + W.p.xy, 
                            W.d*W.v.x*W.v.x*W.v.x + 3.0*W.v.x*W.p.xx,
                            W.d*W.v.x*W.v.x*W.v.y + 2.0*W.v.x*W.p.xy + W.v.y*W.p.xx,
                            W.d*W.v.x*W.v.y*W.v.y + W.v.x*W.p.yy + 2.0*W.v.y*W.p.xy,
                            W.v.x*W.p.zz,
                            W.v.x*W.erot));
}

/********************************************************
 * Gaussian2D_pState::Fx -- Solution flux (x-direction).*
 ********************************************************/

inline Gaussian2D_cState Gaussian2D_pState::Fx(void) {
  return (Gaussian2D_cState(d*v.x, 
                            d*v.x*v.x + p.xx, 
                            d*v.x*v.y + p.xy, 
                            d*v.x*v.x*v.x + 3.0*v.x*p.xx,
                            d*v.x*v.x*v.y + 2.0*v.x*p.xy + v.y*p.xx,
                            d*v.x*v.y*v.y + v.x*p.yy + 2.0*v.y*p.xy,
                            v.x*p.zz,
                            v.x*erot));
}

inline Gaussian2D_cState Gaussian2D_pState::Fx(void) const {
  return (Gaussian2D_cState(d*v.x, 
                            d*v.x*v.x + p.xx, 
                            d*v.x*v.y + p.xy, 
                            d*v.x*v.x*v.x + 3.0*v.x*p.xx,
                            d*v.x*v.x*v.y + 2.0*v.x*p.xy + v.y*p.xx,
                            d*v.x*v.y*v.y + v.x*p.yy + 2.0*v.y*p.xy,
                            v.x*p.zz,
                            v.x*erot));
}

inline Gaussian2D_cState Gaussian2D_pState::Fx(const Gaussian2D_pState &W) {
  return (Gaussian2D_cState(W.d*W.v.x, 
                            W.d*W.v.x*W.v.x + W.p.xx, 
                            W.d*W.v.x*W.v.y + W.p.xy, 
                            W.d*W.v.x*W.v.x*W.v.x + 3.0*W.v.x*W.p.xx,
                            W.d*W.v.x*W.v.x*W.v.y + 2.0*W.v.x*W.p.xy + W.v.y*W.p.xx,
                            W.d*W.v.x*W.v.y*W.v.y + W.v.x*W.p.yy + 2.0*W.v.y*W.p.xy,
                            W.v.x*W.p.zz,
                            W.v.x*W.erot));
}

inline Gaussian2D_cState Fx(const Gaussian2D_pState &W) {
  return (Gaussian2D_cState(W.d*W.v.x, 
                            W.d*W.v.x*W.v.x + W.p.xx, 
                            W.d*W.v.x*W.v.y + W.p.xy, 
                            W.d*W.v.x*W.v.x*W.v.x + 3.0*W.v.x*W.p.xx,
                            W.d*W.v.x*W.v.x*W.v.y + 2.0*W.v.x*W.p.xy + W.v.y*W.p.xx,
                            W.d*W.v.x*W.v.y*W.v.y + W.v.x*W.p.yy + 2.0*W.v.y*W.p.xy,
                            W.v.x*W.p.zz,
                            W.v.x*W.erot));
}

/********************************************************
 * Euler2D_pState::Fy -- Solution flux (y-direction).   *
 ********************************************************/

inline Gaussian2D_cState Gaussian2D_pState::Fy(void) {
  return (Gaussian2D_cState(d*v.y, 
                            d*v.x*v.y + p.xy, 
                            d*v.y*v.y + p.yy, 
                            d*v.x*v.x*v.y + 2.0*v.x*p.xy + v.y*p.xx,
                            d*v.x*v.y*v.y + 2.0*v.y*p.xy + v.x*p.yy,
                            d*v.y*v.y*v.y + 3.0*v.y*p.yy,
                            v.y*p.zz,
                            v.y*erot));
}

inline Gaussian2D_cState Gaussian2D_pState::Fy(void) const {
  return (Gaussian2D_cState(d*v.y, 
                            d*v.x*v.y + p.xy, 
                            d*v.y*v.y + p.yy, 
                            d*v.x*v.x*v.y + 2.0*v.x*p.xy + v.y*p.xx,
                            d*v.x*v.y*v.y + 2.0*v.y*p.xy + v.x*p.yy,
                            d*v.y*v.y*v.y + 3.0*v.y*p.yy,
                            v.y*p.zz,
                            v.y*erot));
}

inline Gaussian2D_cState Gaussian2D_pState::Fy(const Gaussian2D_pState &W) {
  return (Gaussian2D_cState(W.d*W.v.y, 
                            W.d*W.v.x*W.v.y + W.p.xy, 
                            W.d*W.v.y*W.v.y + W.p.yy, 
                            W.d*W.v.x*W.v.x*W.v.y + 2.0*W.v.x*W.p.xy + W.v.y*W.p.xx,
                            W.d*W.v.x*W.v.y*W.v.y + 2.0*W.v.y*W.p.xy + W.v.x*W.p.yy,
                            W.d*W.v.y*W.v.y*W.v.y + 3.0*W.v.y*W.p.yy,
                            W.v.y*W.p.zz,
                            W.v.y*W.erot));
}

inline Gaussian2D_cState Fy(const Gaussian2D_pState &W) {
  return (Gaussian2D_cState(W.d*W.v.y, 
                            W.d*W.v.x*W.v.y + W.p.xy, 
                            W.d*W.v.y*W.v.y + W.p.yy, 
                            W.d*W.v.x*W.v.x*W.v.y + 2.0*W.v.x*W.p.xy + W.v.y*W.p.xx,
                            W.d*W.v.x*W.v.y*W.v.y + 2.0*W.v.y*W.p.xy + W.v.x*W.p.yy,
                            W.d*W.v.y*W.v.y*W.v.y + 3.0*W.v.y*W.p.yy,
                            W.v.y*W.p.zz,
                            W.v.y*W.erot));
}

/********************************************************
 * Euler2D_pState::Fn -- Solution flux (n-direction).   *
 ********************************************************/
/*
inline Euler2D_cState Euler2D_pState::Fn(void) {
  return (Euler2D_cState(d*v.x, d*sqr(v.x) + p, d*v.x*v.y, v.x*H()));
}

inline Euler2D_cState Euler2D_pState::Fn(void) const {
  return (Euler2D_cState(d*v.x, d*sqr(v.x) + p, d*v.x*v.y, v.x*H()));
}

inline Euler2D_cState Euler2D_pState::Fn(const Euler2D_pState &W) {
  return (Euler2D_cState(W.d*W.v.x, W.d*sqr(W.v.x) + W.p,
                         W.d*W.v.x*W.v.y, W.v.x*W.H()));
}

inline Euler2D_cState Fn(const Euler2D_pState &W) {
  return (Euler2D_cState(W.d*W.v.x, W.d*sqr(W.v.x) + W.p,
                         W.d*W.v.x*W.v.y, W.v.x*W.H()));
}
*/
/************************************************************
 * Gaussian2D_pState::lambda -- Eigenvalue(s) (x-direction).*
 ************************************************************/

inline Gaussian2D_pState Gaussian2D_pState::lambda(void) {
  double c = axx();
  return (Gaussian2D_pState(v.x-sqrt(3.0)*c, v.x-c, v.x, v.x, v.x, v.x, v.x+c, v.x+sqrt(3.0)*c));
}

inline Gaussian2D_pState Gaussian2D_pState::lambda(void) const {
  double c = axx();
  return (Gaussian2D_pState(v.x-sqrt(3.0)*c, v.x-c, v.x, v.x, v.x, v.x, v.x+c, v.x+sqrt(3.0)*c));
}

inline Gaussian2D_pState Gaussian2D_pState::lambda(const Gaussian2D_pState &W) {
  double c = W.axx();
  return (Gaussian2D_pState(W.v.x-sqrt(3.0)*c, W.v.x-c, W.v.x, W.v.x, W.v.x, W.v.x, W.v.x+c, W.v.x+sqrt(3.0)*c));
}

inline Gaussian2D_pState lambda(const Gaussian2D_pState &W) {
  double c = W.axx();
  return (Gaussian2D_pState(W.v.x-sqrt(3.0)*c, W.v.x-c, W.v.x, W.v.x, W.v.x, W.v.x, W.v.x+c, W.v.x+sqrt(3.0)*c));
}

inline double Gaussian2D_pState::lambda(int index) {
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  switch(index) {
    case 1 :
      return (v.x-sqrt(3.0)*axx());
    case 2 :
      return (v.x-axx());
    case 3 :
      return (v.x);
    case 4 :
      return (v.x);
    case 5 :
      return (v.x);
    case 6 :
      return (v.x);
    case 7 :
      return (v.x+axx());
    case 8 :
      return (v.x+sqrt(3.0)*axx());
    default:
      return (v.x);
  };
}

inline double Gaussian2D_pState::lambda(int index) const {
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  switch(index) {
    case 1 :
      return (v.x-sqrt(3.0)*axx());
    case 2 :
      return (v.x-axx());
    case 3 :
      return (v.x);
    case 4 :
      return (v.x);
    case 5 :
      return (v.x);
    case 6 :
      return (v.x);
    case 7 :
      return (v.x+axx());
    case 8 :
      return (v.x+sqrt(3.0)*axx());
    default:
      return (v.x);
  };
}

inline double lambda(const Gaussian2D_pState &W, int index) {
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  switch(index) {
    case 1 :
      return (W.v.x-sqrt(3.0)*W.axx());
    case 2 :
      return (W.v.x-W.axx());
    case 3 :
      return (W.v.x);
    case 4 :
      return (W.v.x);
    case 5 :
      return (W.v.x);
    case 6 :
      return (W.v.x);
    case 7 :
      return (W.v.x+W.axx());
    case 8 :
      return (W.v.x+sqrt(3.0)*W.axx());
    default:
      return (W.v.x);
  };
}

/************************************************************
 * Euler2D_pState::lambda_x -- Eigenvalue(s) (x-direction). *
 ************************************************************/
inline Gaussian2D_pState Gaussian2D_pState::lambda_x(void) {
  double c = axx();
  return (Gaussian2D_pState(v.x-sqrt(3.0)*c, v.x-c, v.x, v.x, v.x, v.x, v.x+c, v.x+sqrt(3.0)*c));
}

inline Gaussian2D_pState Gaussian2D_pState::lambda_x(void) const {
  double c = axx();
  return (Gaussian2D_pState(v.x-sqrt(3.0)*c, v.x-c, v.x, v.x, v.x, v.x, v.x+c, v.x+sqrt(3.0)*c));
}

inline Gaussian2D_pState Gaussian2D_pState::lambda_x(const Gaussian2D_pState &W) {
  double c = W.axx();
  return (Gaussian2D_pState(W.v.x-sqrt(3.0)*c, W.v.x-c, W.v.x, W.v.x, W.v.x, W.v.x, W.v.x+c, W.v.x+sqrt(3.0)*c));
}

inline Gaussian2D_pState lambda_x(const Gaussian2D_pState &W) {
  double c = W.axx();
  return (Gaussian2D_pState(W.v.x-sqrt(3.0)*c, W.v.x-c, W.v.x, W.v.x, W.v.x, W.v.x, W.v.x+c, W.v.x+sqrt(3.0)*c));
}

inline double Gaussian2D_pState::lambda_x(int index) {
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  switch(index) {
    case 1 :
      return (v.x-sqrt(3.0)*axx());
    case 2 :
      return (v.x-axx());
    case 3 :
      return (v.x);
    case 4 :
      return (v.x);
    case 5 :
      return (v.x);
    case 6 :
      return (v.x);
    case 7 :
      return (v.x+axx());
    case 8 :
      return (v.x+sqrt(3.0)*axx());
    default:
      return (v.x);
  };
}

inline double Gaussian2D_pState::lambda_x(int index) const {
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  switch(index) {
    case 1 :
      return (v.x-sqrt(3.0)*axx());
    case 2 :
      return (v.x-axx());
    case 3 :
      return (v.x);
    case 4 :
      return (v.x);
    case 5 :
      return (v.x);
    case 6 :
      return (v.x);
    case 7 :
      return (v.x+axx());
    case 8 :
      return (v.x+sqrt(3.0)*axx());
    default:
      return (v.x);
  };
}

inline double lambda_x(const Gaussian2D_pState &W, int index) {
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  switch(index) {
    case 1 :
      return (W.v.x-sqrt(3.0)*W.axx());
    case 2 :
      return (W.v.x-W.axx());
    case 3 :
      return (W.v.x);
    case 4 :
      return (W.v.x);
    case 5 :
      return (W.v.x);
    case 6 :
      return (W.v.x);
    case 7 :
      return (W.v.x+W.axx());
    case 8 :
      return (W.v.x+sqrt(3.0)*W.axx());
    default:
      return (W.v.x);
  };
}

/************************************************************
 * Euler2D_pState::lambda_y -- Eigenvalue(s) (y-direction). *
 ************************************************************/
inline Gaussian2D_pState Gaussian2D_pState::lambda_y(void) {
  double c = ayy();
  return (Gaussian2D_pState(v.y-sqrt(3.0)*c, v.y-c, v.y, v.y, v.y, v.y, v.y+c, v.y+sqrt(3.0)*c));
}

inline Gaussian2D_pState Gaussian2D_pState::lambda_y(void) const {
  double c = ayy();
  return (Gaussian2D_pState(v.y-sqrt(3.0)*c, v.y-c, v.y, v.y, v.y, v.y, v.y+c, v.y+sqrt(3.0)*c));
}

inline Gaussian2D_pState Gaussian2D_pState::lambda_y(const Gaussian2D_pState &W) {
  double c = W.ayy();
  return (Gaussian2D_pState(W.v.y-sqrt(3.0)*c, W.v.y-c, W.v.y, W.v.y, W.v.y, W.v.y, W.v.y+c, W.v.y+sqrt(3.0)*c));
}

inline Gaussian2D_pState lambda_y(const Gaussian2D_pState &W) {
  double c = W.ayy();
  return (Gaussian2D_pState(W.v.y-sqrt(3.0)*c, W.v.y-c, W.v.y, W.v.y, W.v.y, W.v.y, W.v.y+c, W.v.y+sqrt(3.0)*c));
}

inline double Gaussian2D_pState::lambda_y(int index) {
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  switch(index) {
    case 1 :
      return (v.y-sqrt(3.0)*ayy());
    case 2 :
      return (v.y-ayy());
    case 3 :
      return (v.y);
    case 4 :
      return (v.y);
    case 5 :
      return (v.y);
    case 6 :
      return (v.y);
    case 7 :
      return (v.y+ayy());
    case 8 :
      return (v.y+sqrt(3.0)*ayy());
    default:
      return (v.y);
  };
}

inline double Gaussian2D_pState::lambda_y(int index) const {
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  switch(index) {
    case 1 :
      return (v.y-sqrt(3.0)*ayy());
    case 2 :
      return (v.y-ayy());
    case 3 :
      return (v.y);
    case 4 :
      return (v.y);
    case 5 :
      return (v.y);
    case 6 :
      return (v.y);
    case 7 :
      return (v.y+ayy());
    case 8 :
      return (v.y+sqrt(3.0)*ayy());
    default:
      return (v.y);
  };
}

inline double lambda_y(const Gaussian2D_pState &W, int index) {
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  switch(index) {
    case 1 :
      return (W.v.y-sqrt(3.0)*W.ayy());
    case 2 :
      return (W.v.y-W.ayy());
    case 3 :
      return (W.v.y);
    case 4 :
      return (W.v.y);
    case 5 :
      return (W.v.y);
    case 6 :
      return (W.v.y);
    case 7 :
      return (W.v.y+W.ayy());
    case 8 :
      return (W.v.y+sqrt(3.0)*W.ayy());
    default:
      return (W.v.y);
  };
}

/********************************************************
 * Gaussian2D_pState::rc -- Conserved right eigenvector *
 *                       (x-direction).                 *
 ********************************************************/

inline Gaussian2D_cState Gaussian2D_pState::rc(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*c, v.y-sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c-2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c-v.x*sqrt(3.0)*p.xy-sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c-2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
    case 2 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,v.x-c,2.0*(v.y-p.xy/(c*d)),ZERO,ZERO));
    case 3 :
      return (Gaussian2D_cState(ONE,v.x,v.y,v.x*v.x,v.x*v.y,v.y*v.y,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,v.x+c,2.0*(v.y+p.xy/(c*d)),ZERO,ZERO));
    case 8 :
      return (Gaussian2D_cState(ONE, v.x+sqrt(3.0)*c, v.y+sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c+2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c+v.x*sqrt(3.0)*p.xy+sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c+2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
    default:
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*c, v.y-sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c-2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c-v.x*sqrt(3.0)*p.xy-sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c-2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
  };
}

inline Gaussian2D_cState Gaussian2D_pState::rc(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*c, v.y-sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c-2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c-v.x*sqrt(3.0)*p.xy-sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c-2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
    case 2 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,v.x-c,2.0*(v.y-p.xy/(c*d)),ZERO,ZERO));
    case 3 :
      return (Gaussian2D_cState(ONE,v.x,v.y,v.x*v.x,v.x*v.y,v.y*v.y,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,v.x+c,2.0*(v.y+p.xy/(c*d)),ZERO,ZERO));
    case 8 :
      return (Gaussian2D_cState(ONE, v.x+sqrt(3.0)*c, v.y+sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c+2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c+v.x*sqrt(3.0)*p.xy+sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c+2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
    default:
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*c, v.y-sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c-2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c-v.x*sqrt(3.0)*p.xy-sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c-2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
  };
}

inline Gaussian2D_cState rc(const Gaussian2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = W.axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_cState(ONE, W.v.x-sqrt(3.0)*c, W.v.y-sqrt(3.0)*W.p.xy/(c*W.d),
               3.0*c*c-2.0*sqrt(3.0)*W.v.x*c+W.v.x*W.v.x,
               (W.v.x*W.d*W.v.y*c-W.v.x*sqrt(3.0)*W.p.xy-sqrt(3.0)*c*c*W.d*W.v.y+3.0*c*W.p.xy)/(c*W.d),
               (W.d*W.d*W.v.y*W.v.y*c*c-2.0*sqrt(3.0)*W.p.xy*c*W.d*W.v.y+W.d*c*c*W.p.yy+2.0*W.p.xy*W.p.xy)/(W.d*W.d*c*c),
               W.p.zz/W.d,W.erot/W.d));
    case 2 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,W.v.x-c,2.0*(W.v.y-W.p.xy/(c*W.d)),ZERO,ZERO));
    case 3 :
      return (Gaussian2D_cState(ONE,W.v.x,W.v.y,W.v.x*W.v.x,W.v.x*W.v.y,W.v.y*W.v.y,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,W.v.x+c,2.0*(W.v.y+W.p.xy/(c*W.d)),ZERO,ZERO));
    case 8 :
      return (Gaussian2D_cState(ONE, W.v.x+sqrt(3.0)*c, W.v.y+sqrt(3.0)*W.p.xy/(c*W.d),
               3.0*c*c+2.0*sqrt(3.0)*W.v.x*c+W.v.x*W.v.x,
               (W.v.x*W.d*W.v.y*c+W.v.x*sqrt(3.0)*W.p.xy+sqrt(3.0)*c*c*W.d*W.v.y+3.0*c*W.p.xy)/(c*W.d),
               (W.d*W.d*W.v.y*W.v.y*c*c+2.0*sqrt(3.0)*W.p.xy*c*W.d*W.v.y+W.d*c*c*W.p.yy+2.0*W.p.xy*W.p.xy)/(W.d*W.d*c*c),
               W.p.zz/W.d,W.erot/W.d));
    default:
      return (Gaussian2D_cState(ONE, W.v.x-sqrt(3.0)*c, W.v.y-sqrt(3.0)*W.p.xy/(c*W.d),
               3.0*c*c-2.0*sqrt(3.0)*W.v.x*c+W.v.x*W.v.x,
               (W.v.x*W.d*W.v.y*c-W.v.x*sqrt(3.0)*W.p.xy-sqrt(3.0)*c*c*W.d*W.v.y+3.0*c*W.p.xy)/(c*W.d),
               (W.d*W.d*W.v.y*W.v.y*c*c-2.0*sqrt(3.0)*W.p.xy*c*W.d*W.v.y+W.d*c*c*W.p.yy+2.0*W.p.xy*W.p.xy)/(W.d*W.d*c*c),
               W.p.zz/W.d,W.erot/W.d));
  };
}

/*********************************************************
 * Gaussian2D_pState::rc_x -- Conserved right eigenvector*
 *                         (x-direction).                *
 *********************************************************/
inline Gaussian2D_cState Gaussian2D_pState::rc_x(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*c, v.y-sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c-2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c-v.x*sqrt(3.0)*p.xy-sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c-2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
    case 2 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,v.x-c,2.0*(v.y-p.xy/(c*d)),ZERO,ZERO));
    case 3 :
      return (Gaussian2D_cState(ONE,v.x,v.y,v.x*v.x,v.x*v.y,v.y*v.y,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,v.x+c,2.0*(v.y+p.xy/(c*d)),ZERO,ZERO));
    case 8 :
      return (Gaussian2D_cState(ONE, v.x+sqrt(3.0)*c, v.y+sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c+2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c+v.x*sqrt(3.0)*p.xy+sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c+2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
    default:
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*c, v.y-sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c-2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c-v.x*sqrt(3.0)*p.xy-sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c-2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
  };
}

inline Gaussian2D_cState Gaussian2D_pState::rc_x(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*c, v.y-sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c-2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c-v.x*sqrt(3.0)*p.xy-sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c-2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
    case 2 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,v.x-c,2.0*(v.y-p.xy/(c*d)),ZERO,ZERO));
    case 3 :
      return (Gaussian2D_cState(ONE,v.x,v.y,v.x*v.x,v.x*v.y,v.y*v.y,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,v.x+c,2.0*(v.y+p.xy/(c*d)),ZERO,ZERO));
    case 8 :
      return (Gaussian2D_cState(ONE, v.x+sqrt(3.0)*c, v.y+sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c+2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c+v.x*sqrt(3.0)*p.xy+sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c+2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
    default:
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*c, v.y-sqrt(3.0)*p.xy/(c*d),
                                3.0*c*c-2.0*sqrt(3.0)*v.x*c+v.x*v.x,
                                (v.x*d*v.y*c-v.x*sqrt(3.0)*p.xy-sqrt(3.0)*c*c*d*v.y+3.0*c*p.xy)/(c*d),
                                (d*d*v.y*v.y*c*c-2.0*sqrt(3.0)*p.xy*c*d*v.y+d*c*c*p.yy+2.0*p.xy*p.xy)/(d*d*c*c),
                                p.zz/d,erot/d));
  };
}

inline Gaussian2D_cState rc_x(const Gaussian2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = W.axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_cState(ONE, W.v.x-sqrt(3.0)*c, W.v.y-sqrt(3.0)*W.p.xy/(c*W.d),
               3.0*c*c-2.0*sqrt(3.0)*W.v.x*c+W.v.x*W.v.x,
               (W.v.x*W.d*W.v.y*c-W.v.x*sqrt(3.0)*W.p.xy-sqrt(3.0)*c*c*W.d*W.v.y+3.0*c*W.p.xy)/(c*W.d),
               (W.d*W.d*W.v.y*W.v.y*c*c-2.0*sqrt(3.0)*W.p.xy*c*W.d*W.v.y+W.d*c*c*W.p.yy+2.0*W.p.xy*W.p.xy)/(W.d*W.d*c*c),
               W.p.zz/W.d,W.erot/W.d));
    case 2 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,W.v.x-c,2.0*(W.v.y-W.p.xy/(c*W.d)),ZERO,ZERO));
    case 3 :
      return (Gaussian2D_cState(ONE,W.v.x,W.v.y,W.v.x*W.v.x,W.v.x*W.v.y,W.v.y*W.v.y,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_cState(ZERO,ZERO,ONE,ZERO,W.v.x+c,2.0*(W.v.y+W.p.xy/(c*W.d)),ZERO,ZERO));
    case 8 :
      return (Gaussian2D_cState(ONE, W.v.x+sqrt(3.0)*c, W.v.y+sqrt(3.0)*W.p.xy/(c*W.d),
               3.0*c*c+2.0*sqrt(3.0)*W.v.x*c+W.v.x*W.v.x,
               (W.v.x*W.d*W.v.y*c+W.v.x*sqrt(3.0)*W.p.xy+sqrt(3.0)*c*c*W.d*W.v.y+3.0*c*W.p.xy)/(c*W.d),
               (W.d*W.d*W.v.y*W.v.y*c*c+2.0*sqrt(3.0)*W.p.xy*c*W.d*W.v.y+W.d*c*c*W.p.yy+2.0*W.p.xy*W.p.xy)/(W.d*W.d*c*c),
               W.p.zz/W.d,W.erot/W.d));
    default:
      return (Gaussian2D_cState(ONE, W.v.x-sqrt(3.0)*c, W.v.y-sqrt(3.0)*W.p.xy/(c*W.d),
               3.0*c*c-2.0*sqrt(3.0)*W.v.x*c+W.v.x*W.v.x,
               (W.v.x*W.d*W.v.y*c-W.v.x*sqrt(3.0)*W.p.xy-sqrt(3.0)*c*c*W.d*W.v.y+3.0*c*W.p.xy)/(c*W.d),
               (W.d*W.d*W.v.y*W.v.y*c*c-2.0*sqrt(3.0)*W.p.xy*c*W.d*W.v.y+W.d*c*c*W.p.yy+2.0*W.p.xy*W.p.xy)/(W.d*W.d*c*c),
               W.p.zz/W.d,W.erot/W.d));
  };
}

/*********************************************************
 * Gaussian2D_pState::rc_y -- Conserved right eigenvector*
 *                         (y-direction).                *
 *********************************************************/
inline Gaussian2D_cState Gaussian2D_pState::rc_y(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = ayy();
  switch(index) {
    case 1 :
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*p.xy/(c*d),v.y-sqrt(3.0)*c,
                                (d*d*v.x*v.x*c*c-2.0*sqrt(3.0)*p.xy*d*v.x*c+p.xx*d*c*c+2.0*p.xy*p.xy)/(c*c*d*d),
                                (v.y*d*v.x*c-sqrt(3.0)*v.y*p.xy-sqrt(3.0)*c*c*d*v.x+3.0*c*p.xy)/(d*c),
                                3.0*c*c-2.0*sqrt(3.0)*v.y*c+v.y*v.y,
                                p.zz/d,erot/d));
    case 2 :
      return (Gaussian2D_cState(ZERO,ONE,ZERO,2.0*(v.x-p.xy/(c*d)),v.y-c,ZERO,ZERO,ZERO));
    case 3 :
      return (Gaussian2D_cState(ONE,v.x,v.y,v.x*v.x,v.x*v.y,v.y*v.y,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ONE,ZERO,ZERO,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_cState(ZERO,ONE,ZERO,2.0*(v.x+p.xy/(c*d)),v.y+c,ZERO,ZERO,ZERO));
    case 8 :
      return (Gaussian2D_cState(ONE, v.x+sqrt(3.0)*p.xy/(c*d),v.y+sqrt(3.0)*c,
                                (d*d*v.x*v.x*c*c+2.0*sqrt(3.0)*p.xy*d*v.x*c+p.xx*d*c*c+2.0*p.xy*p.xy)/(c*c*d*d),
                                (v.y*d*v.x*c+sqrt(3.0)*v.y*p.xy+sqrt(3.0)*c*c*d*v.x+3.0*c*p.xy)/(d*c),
                                3.0*c*c+2.0*sqrt(3.0)*v.y*c+v.y*v.y,
                                p.zz/d,erot/d));
    default:
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*p.xy/(c*d),v.y-sqrt(3.0)*c,
                                (d*d*v.x*v.x*c*c-2.0*sqrt(3.0)*p.xy*d*v.x*c+p.xx*d*c*c+2.0*p.xy*p.xy)/(c*c*d*d),
                                (v.y*d*v.x*c-sqrt(3.0)*v.y*p.xy-sqrt(3.0)*c*c*d*v.x+3.0*c*p.xy)/(d*c),
                                3.0*c*c-2.0*sqrt(3.0)*v.y*c+v.y*v.y,
                                p.zz/d,erot/d));
  };
}

inline Gaussian2D_cState Gaussian2D_pState::rc_y(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = ayy();
  switch(index) {
    case 1 :
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*p.xy/(c*d),v.y-sqrt(3.0)*c,
                                (d*d*v.x*v.x*c*c-2.0*sqrt(3.0)*p.xy*d*v.x*c+p.xx*d*c*c+2.0*p.xy*p.xy)/(c*c*d*d),
                                (v.y*d*v.x*c-sqrt(3.0)*v.y*p.xy-sqrt(3.0)*c*c*d*v.x+3.0*c*p.xy)/(d*c),
                                3.0*c*c-2.0*sqrt(3.0)*v.y*c+v.y*v.y,
                                p.zz/d,erot/d));
    case 2 :
      return (Gaussian2D_cState(ZERO,ONE,ZERO,2.0*(v.x-p.xy/(c*d)),v.y-c,ZERO,ZERO,ZERO));
    case 3 :
      return (Gaussian2D_cState(ONE,v.x,v.y,v.x*v.x,v.x*v.y,v.y*v.y,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ONE,ZERO,ZERO,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_cState(ZERO,ONE,ZERO,2.0*(v.x+p.xy/(c*d)),v.y+c,ZERO,ZERO,ZERO));
    case 8 :
      return (Gaussian2D_cState(ONE, v.x+sqrt(3.0)*p.xy/(c*d),v.y+sqrt(3.0)*c,
                                (d*d*v.x*v.x*c*c+2.0*sqrt(3.0)*p.xy*d*v.x*c+p.xx*d*c*c+2.0*p.xy*p.xy)/(c*c*d*d),
                                (v.y*d*v.x*c+sqrt(3.0)*v.y*p.xy+sqrt(3.0)*c*c*d*v.x+3.0*c*p.xy)/(d*c),
                                3.0*c*c+2.0*sqrt(3.0)*v.y*c+v.y*v.y,
                                p.zz/d,erot/d));
    default:
      return (Gaussian2D_cState(ONE, v.x-sqrt(3.0)*p.xy/(c*d),v.y-sqrt(3.0)*c,
                                (d*d*v.x*v.x*c*c-2.0*sqrt(3.0)*p.xy*d*v.x*c+p.xx*d*c*c+2.0*p.xy*p.xy)/(c*c*d*d),
                                (v.y*d*v.x*c-sqrt(3.0)*v.y*p.xy-sqrt(3.0)*c*c*d*v.x+3.0*c*p.xy)/(d*c),
                                3.0*c*c-2.0*sqrt(3.0)*v.y*c+v.y*v.y,
                                p.zz/d,erot/d));
  };
}

inline Gaussian2D_cState rc_y(const Gaussian2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = W.ayy();
  switch(index) {
    case 1 :
      return (Gaussian2D_cState(ONE, W.v.x-sqrt(3.0)*W.p.xy/(c*W.d),W.v.y-sqrt(3.0)*c,
               (W.d*W.d*W.v.x*W.v.x*c*c-2.0*sqrt(3.0)*W.p.xy*W.d*W.v.x*c+W.p.xx*W.d*c*c+2.0*W.p.xy*W.p.xy)/(c*c*W.d*W.d),
               (W.v.y*W.d*W.v.x*c-sqrt(3.0)*W.v.y*W.p.xy-sqrt(3.0)*c*c*W.d*W.v.x+3.0*c*W.p.xy)/(W.d*c),
               3.0*c*c-2.0*sqrt(3.0)*W.v.y*c+W.v.y*W.v.y,
               W.p.zz/W.d,W.erot/W.d));
    case 2 :
      return (Gaussian2D_cState(ZERO,ONE,ZERO,2.0*(W.v.x-W.p.xy/(c*W.d)),W.v.y-c,ZERO,ZERO,ZERO));
    case 3 :
      return (Gaussian2D_cState(ONE,W.v.x,W.v.y,W.v.x*W.v.x,W.v.x*W.v.y,W.v.y*W.v.y,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ONE,ZERO,ZERO,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_cState(ZERO,ONE,ZERO,2.0*(W.v.x+W.p.xy/(c*W.d)),W.v.y+c,ZERO,ZERO,ZERO));
    case 8 :
      return (Gaussian2D_cState(ONE, W.v.x+sqrt(3.0)*W.p.xy/(c*W.d),W.v.y+sqrt(3.0)*c,
               (W.d*W.d*W.v.x*W.v.x*c*c+2.0*sqrt(3.0)*W.p.xy*W.d*W.v.x*c+W.p.xx*W.d*c*c+2.0*W.p.xy*W.p.xy)/(c*c*W.d*W.d),
               (W.v.y*W.d*W.v.x*c+sqrt(3.0)*W.v.y*W.p.xy+sqrt(3.0)*c*c*W.d*W.v.x+3.0*c*W.p.xy)/(W.d*c),
               3.0*c*c+2.0*sqrt(3.0)*W.v.y*c+W.v.y*W.v.y,
               W.p.zz/W.d,W.erot/W.d));
    default:
      return (Gaussian2D_cState(ONE, W.v.x-sqrt(3.0)*W.p.xy/(c*W.d),W.v.y-sqrt(3.0)*c,
               (W.d*W.d*W.v.x*W.v.x*c*c-2.0*sqrt(3.0)*W.p.xy*W.d*W.v.x*c+W.p.xx*W.d*c*c+2.0*W.p.xy*W.p.xy)/(c*c*W.d*W.d),
               (W.v.y*W.d*W.v.x*c-sqrt(3.0)*W.v.y*W.p.xy-sqrt(3.0)*c*c*W.d*W.v.x+3.0*c*W.p.xy)/(W.d*c),
               3.0*c*c-2.0*sqrt(3.0)*W.v.y*c+W.v.y*W.v.y,
               W.p.zz/W.d,W.erot/W.d));
  };
}


/********************************************************
 * Gaussian2D_pState::lp -- Primitive left eigenvector     *
 *                       (x-direction).                 *
 ********************************************************/

inline Gaussian2D_pState Gaussian2D_pState::lp(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 2 :
      return (Gaussian2D_pState(ZERO,-p.xy/(2.0*c*c),d/2.0,p.xy/(2.0*d*c*c*c),-ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 3 :
      return (Gaussian2D_pState(ONE,ZERO,ZERO,-1.0/(3.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,(4.0*p.xy*p.xy-d*c*c*p.yy)/(3.0*d*d*c*c*c*c),
                                -2.0*p.xy/(c*c*d),ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-p.zz/(3.0*c*c*d),ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-erot/(3.0*c*c*d),ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_pState(ZERO,-p.xy/(2.0*c*c),d/2.0,-p.xy/(2.0*d*c*c*c),ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 8 :
      return (Gaussian2D_pState(ZERO, sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    default:
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
  };
}

inline Gaussian2D_pState Gaussian2D_pState::lp(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 2 :
      return (Gaussian2D_pState(ZERO,-p.xy/(2.0*c*c),d/2.0,p.xy/(2.0*d*c*c*c),-ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 3 :
      return (Gaussian2D_pState(ONE,ZERO,ZERO,-1.0/(3.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,(4.0*p.xy*p.xy-d*c*c*p.yy)/(3.0*d*d*c*c*c*c),
                                -2.0*p.xy/(c*c*d),ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-p.zz/(3.0*c*c*d),ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-erot/(3.0*c*c*d),ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_pState(ZERO,-p.xy/(2.0*c*c),d/2.0,-p.xy/(2.0*d*c*c*c),ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 8 :
      return (Gaussian2D_pState(ZERO, sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    default:
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
  };
}

inline Gaussian2D_pState lp(const Gaussian2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = W.axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*W.d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 2 :
      return (Gaussian2D_pState(ZERO,-W.p.xy/(2.0*c*c),W.d/2.0,W.p.xy/(2.0*W.d*c*c*c),-ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 3 :
      return (Gaussian2D_pState(ONE,ZERO,ZERO,-1.0/(3.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,(4.0*W.p.xy*W.p.xy-W.d*c*c*W.p.yy)/(3.0*W.d*W.d*c*c*c*c),
                                -2.0*W.p.xy/(c*c*W.d),ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-W.p.zz/(3.0*c*c*W.d),ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-W.erot/(3.0*c*c*W.d),ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_pState(ZERO,-W.p.xy/(2.0*c*c),W.d/2.0,-W.p.xy/(2.0*W.d*c*c*c),ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 8 :
      return (Gaussian2D_pState(ZERO, sqrt(3.0)/6.0*W.d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    default:
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*W.d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
  };
}

/********************************************************
 * Gaussian2D_pState::lp_x -- Primitive left eigenvector   *
 *                         (x-direction).               *
 ********************************************************/
inline Gaussian2D_pState Gaussian2D_pState::lp_x(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 2 :
      return (Gaussian2D_pState(ZERO,-p.xy/(2.0*c*c),d/2.0,p.xy/(2.0*d*c*c*c),-ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 3 :
      return (Gaussian2D_pState(ONE,ZERO,ZERO,-1.0/(3.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,(4.0*p.xy*p.xy-d*c*c*p.yy)/(3.0*d*d*c*c*c*c),
                                -2.0*p.xy/(c*c*d),ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-p.zz/(3.0*c*c*d),ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-erot/(3.0*c*c*d),ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_pState(ZERO,-p.xy/(2.0*c*c),d/2.0,-p.xy/(2.0*d*c*c*c),ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 8 :
      return (Gaussian2D_pState(ZERO, sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    default:
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
  };
}

inline Gaussian2D_pState Gaussian2D_pState::lp_x(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 2 :
      return (Gaussian2D_pState(ZERO,-p.xy/(2.0*c*c),d/2.0,p.xy/(2.0*d*c*c*c),-ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 3 :
      return (Gaussian2D_pState(ONE,ZERO,ZERO,-1.0/(3.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,(4.0*p.xy*p.xy-d*c*c*p.yy)/(3.0*d*d*c*c*c*c),
                                -2.0*p.xy/(c*c*d),ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-p.zz/(3.0*c*c*d),ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-erot/(3.0*c*c*d),ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_pState(ZERO,-p.xy/(2.0*c*c),d/2.0,-p.xy/(2.0*d*c*c*c),ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 8 :
      return (Gaussian2D_pState(ZERO, sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    default:
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
  };
}

inline Gaussian2D_pState lp_x(const Gaussian2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = W.axx();
  switch(index) {
    case 1 :
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*W.d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 2 :
      return (Gaussian2D_pState(ZERO,-W.p.xy/(2.0*c*c),W.d/2.0,W.p.xy/(2.0*W.d*c*c*c),-ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 3 :
      return (Gaussian2D_pState(ONE,ZERO,ZERO,-1.0/(3.0*c*c),ZERO,ZERO,ZERO,ZERO));
    case 4 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,(4.0*W.p.xy*W.p.xy-W.d*c*c*W.p.yy)/(3.0*W.d*W.d*c*c*c*c),
                                -2.0*W.p.xy/(c*c*W.d),ONE,ZERO,ZERO));
    case 5 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-W.p.zz/(3.0*c*c*W.d),ZERO,ZERO,ONE,ZERO));
    case 6 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,-W.erot/(3.0*c*c*W.d),ZERO,ZERO,ZERO,ONE));
    case 7 :
      return (Gaussian2D_pState(ZERO,-W.p.xy/(2.0*c*c),W.d/2.0,-W.p.xy/(2.0*W.d*c*c*c),ONE/(2.0*c),ZERO,ZERO,ZERO));
    case 8 :
      return (Gaussian2D_pState(ZERO, sqrt(3.0)/6.0*W.d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
    default:
      return (Gaussian2D_pState(ZERO, -sqrt(3.0)/6.0*W.d/c, ZERO, ONE/(6.0*c*c),ZERO,ZERO,ZERO,ZERO));
  };
}

/********************************************************
 * Gaussian2D_pState::lp_y -- Primitive left eigenvector*
 *                         (y-direction).               *
 ********************************************************/
inline Gaussian2D_pState Gaussian2D_pState::lp_y(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = ayy();
  switch(index) {
    case 1 :
      return (Gaussian2D_pState(ZERO, ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ZERO, ONE/(6.0*c*c),ZERO,ZERO));
    case 2 :
      return (Gaussian2D_pState(ZERO,d/2.0,-p.xy/(2.0*c*c),ZERO,-ONE/(2*c),p.xy/(2*d*c*c*c),ZERO,ZERO));
    case 3 :
      return (Gaussian2D_pState(ONE,ZERO,ZERO,ZERO,ZERO,-1.0/(3.0*c*c),ZERO,ZERO));
    case 4 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,ONE,-2.0*p.xy/(d*c*c),(4.0*p.xy*p.xy-p.xx*d*c*c)/(3.0*d*d*c*c*c*c),
                                ZERO,ZERO));
    case 5 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,-p.zz/(3.0*c*c*d),ONE,ZERO));
    case 6 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,-erot/(3.0*c*c*d),ZERO,ONE));
    case 7 :
      return (Gaussian2D_pState(ZERO,d/2.0,-p.xy/(2.0*c*c),ZERO,ONE/(2*c),-p.xy/(2*d*c*c*c),ZERO,ZERO));
    case 8 :
      return (Gaussian2D_pState(ZERO, ZERO, sqrt(3.0)/6.0*d/c, ZERO, ZERO, ONE/(6.0*c*c),ZERO,ZERO));
    default:
      return (Gaussian2D_pState(ZERO, ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ZERO, ONE/(6.0*c*c),ZERO,ZERO));
  };
}

inline Gaussian2D_pState Gaussian2D_pState::lp_y(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = ayy();
  switch(index) {
    case 1 :
      return (Gaussian2D_pState(ZERO, ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ZERO, ONE/(6.0*c*c),ZERO,ZERO));
    case 2 :
      return (Gaussian2D_pState(ZERO,d/2.0,-p.xy/(2.0*c*c),ZERO,-ONE/(2*c),p.xy/(2*d*c*c*c),ZERO,ZERO));
    case 3 :
      return (Gaussian2D_pState(ONE,ZERO,ZERO,ZERO,ZERO,-1.0/(3.0*c*c),ZERO,ZERO));
    case 4 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,ONE,-2.0*p.xy/(d*c*c),(4.0*p.xy*p.xy-p.xx*d*c*c)/(3.0*d*d*c*c*c*c),
                                ZERO,ZERO));
    case 5 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,-p.zz/(3.0*c*c*d),ONE,ZERO));
    case 6 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,-erot/(3.0*c*c*d),ZERO,ONE));
    case 7 :
      return (Gaussian2D_pState(ZERO,d/2.0,-p.xy/(2.0*c*c),ZERO,ONE/(2*c),-p.xy/(2*d*c*c*c),ZERO,ZERO));
    case 8 :
      return (Gaussian2D_pState(ZERO, ZERO, sqrt(3.0)/6.0*d/c, ZERO, ZERO, ONE/(6.0*c*c),ZERO,ZERO));
    default:
      return (Gaussian2D_pState(ZERO, ZERO, -sqrt(3.0)/6.0*d/c, ZERO, ZERO, ONE/(6.0*c*c),ZERO,ZERO));
  };
}

inline Gaussian2D_pState lp_y(const Gaussian2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_GAUSSIAN2D );
  c = W.ayy();
  switch(index) {
    case 1 :
      return (Gaussian2D_pState(ZERO, ZERO, -sqrt(3.0)/6.0*W.d/c, ZERO, ZERO, ONE/(6.0*c*c),ZERO,ZERO));
    case 2 :
      return (Gaussian2D_pState(ZERO,W.d/2.0,-W.p.xy/(2.0*c*c),ZERO,-ONE/(2*c),W.p.xy/(2*W.d*c*c*c),ZERO,ZERO));
    case 3 :
      return (Gaussian2D_pState(ONE,ZERO,ZERO,ZERO,ZERO,-1.0/(3.0*c*c),ZERO,ZERO));
    case 4 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,ONE,-2.0*W.p.xy/(W.d*c*c),
                                (4.0*W.p.xy*W.p.xy-W.p.xx*W.d*c*c)/(3.0*W.d*W.d*c*c*c*c),ZERO,ZERO));
    case 5 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,-W.p.zz/(3.0*c*c*W.d),ONE,ZERO));
    case 6 :
      return (Gaussian2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,-W.erot/(3.0*c*c*W.d),ZERO,ONE));
    case 7 :
      return (Gaussian2D_pState(ZERO,W.d/2.0,-W.p.xy/(2.0*c*c),ZERO,ONE/(2*c),-W.p.xy/(2*W.d*c*c*c),ZERO,ZERO));
    case 8 :
      return (Gaussian2D_pState(ZERO, ZERO, sqrt(3.0)/6.0*W.d/c, ZERO, ZERO, ONE/(6.0*c*c),ZERO,ZERO));
    default:
      return (Gaussian2D_pState(ZERO, ZERO, -sqrt(3.0)/6.0*W.d/c, ZERO, ZERO, ONE/(6.0*c*c),ZERO,ZERO));
  };
}

/********************************************************
 * Gaussian2D_cState::Gaussian2D_cState -- Constructor. *
 ********************************************************/

inline Gaussian2D_cState::Gaussian2D_cState(const Gaussian2D_pState &W) {
  d = W.d; dv = W.dv(); E = W.E(); erot = W.erot;
}

/********************************************************
 * Gaussian2D_cState::W -- Primitive solution state.    *
 ********************************************************/

inline Gaussian2D_pState Gaussian2D_cState::W(void) {
  return (Gaussian2D_pState(d, v(), p(), erot));
}

inline Gaussian2D_pState Gaussian2D_cState::W(void) const {
  return (Gaussian2D_pState(d, v(), p(), erot));
}

inline Gaussian2D_pState Gaussian2D_cState::W(const Gaussian2D_cState &U) {
  return (Gaussian2D_pState(U.d, U.v(), U.p(), U.erot));
}

inline Gaussian2D_pState W(const Gaussian2D_cState &U) {
  return (Gaussian2D_pState(U.d, U.v(), U.p(), U.erot));
}

/********************************************************
 * Gaussian2D_cState::F -- Solution flux (x-direction).    *
 ********************************************************/

inline Gaussian2D_cState Gaussian2D_cState::F(void) {
  return (Gaussian2D_cState(dv.x, 
                            dv.x*dv.x/d + p().xx,
                            dv.x*dv.y/d + p().xy, 
                            E.xx*dv.x/d + 2.0*dv.x/d*p().xx,
                            E.xy*dv.x/d + dv.x/d*p().xy + dv.y/d*p().xx,
                            E.yy*dv.x/d + 2.0*dv.y/d*p().xy,
                            E.zz*dv.x/d,
                            erot*dv.x/d));
}

inline Gaussian2D_cState Gaussian2D_cState::F(void) const {
  return (Gaussian2D_cState(dv.x, 
                            dv.x*dv.x/d + p().xx,
                            dv.x*dv.y/d + p().xy, 
                            E.xx*dv.x/d + 2.0*dv.x/d*p().xx,
                            E.xy*dv.x/d + dv.x/d*p().xy + dv.y/d*p().xx,
                            E.yy*dv.x/d + 2.0*dv.y/d*p().xy,
                            E.zz*dv.x/d,
                            erot*dv.x/d));
}

inline Gaussian2D_cState Gaussian2D_cState::F(const Gaussian2D_cState &U) {
  return (Gaussian2D_cState(U.dv.x, 
                            U.dv.x*U.dv.x/U.d + U.p().xx,
                            U.dv.x*U.dv.y/U.d + U.p().xy, 
                            U.E.xx*U.dv.x/U.d + 2.0*U.dv.x/U.d*U.p().xx,
                            U.E.xy*U.dv.x/U.d + U.dv.x/U.d*U.p().xy + U.dv.y/U.d*U.p().xx,
                            U.E.yy*U.dv.x/U.d + 2.0*U.dv.y/U.d*U.p().xy,
                            U.E.zz*U.dv.x/U.d,
                            U.erot*U.dv.x/U.d));
}

inline Gaussian2D_cState F(const Gaussian2D_cState &U) {
  return (Gaussian2D_cState(U.dv.x, 
                            U.dv.x*U.dv.x/U.d + U.p().xx,
                            U.dv.x*U.dv.y/U.d + U.p().xy, 
                            U.E.xx*U.dv.x/U.d + 2.0*U.dv.x/U.d*U.p().xx,
                            U.E.xy*U.dv.x/U.d + U.dv.x/U.d*U.p().xy + U.dv.y/U.d*U.p().xx,
                            U.E.yy*U.dv.x/U.d + 2.0*U.dv.y/U.d*U.p().xy,
                            U.E.zz*U.dv.x/U.d,
                            U.erot*U.dv.x/U.d));
}

/********************************************************
 * Euler2D_cState::Fx -- Solution flux (x-direction).   *
 ********************************************************/

inline Gaussian2D_cState Gaussian2D_cState::Fx(void) {
  return (Gaussian2D_cState(dv.x, 
                            dv.x*dv.x/d + p().xx,
                            dv.x*dv.y/d + p().xy, 
                            E.xx*dv.x/d + 2.0*dv.x/d*p().xx,
                            E.xy*dv.x/d + dv.x/d*p().xy + dv.y/d*p().xx,
                            E.yy*dv.x/d + 2.0*dv.y/d*p().xy,
                            E.zz*dv.x/d,
                            erot*dv.x/d));
}

inline Gaussian2D_cState Gaussian2D_cState::Fx(void) const {
  return (Gaussian2D_cState(dv.x, 
                            dv.x*dv.x/d + p().xx,
                            dv.x*dv.y/d + p().xy, 
                            E.xx*dv.x/d + 2.0*dv.x/d*p().xx,
                            E.xy*dv.x/d + dv.x/d*p().xy + dv.y/d*p().xx,
                            E.yy*dv.x/d + 2.0*dv.y/d*p().xy,
                            E.zz*dv.x/d,
                            erot*dv.x/d));
}

inline Gaussian2D_cState Gaussian2D_cState::Fx(const Gaussian2D_cState &U) {
  return (Gaussian2D_cState(U.dv.x, 
                            U.dv.x*U.dv.x/U.d + U.p().xx,
                            U.dv.x*U.dv.y/U.d + U.p().xy, 
                            U.E.xx*U.dv.x/U.d + 2.0*U.dv.x/U.d*U.p().xx,
                            U.E.xy*U.dv.x/U.d + U.dv.x/U.d*U.p().xy + U.dv.y/U.d*U.p().xx,
                            U.E.yy*U.dv.x/U.d + 2.0*U.dv.y/U.d*U.p().xy,
                            U.E.zz*U.dv.x/U.d,
                            U.erot*U.dv.x/U.d));
}

inline Gaussian2D_cState Fx(const Gaussian2D_cState &U) {
  return (Gaussian2D_cState(U.dv.x, 
                            U.dv.x*U.dv.x/U.d + U.p().xx,
                            U.dv.x*U.dv.y/U.d + U.p().xy, 
                            U.E.xx*U.dv.x/U.d + 2.0*U.dv.x/U.d*U.p().xx,
                            U.E.xy*U.dv.x/U.d + U.dv.x/U.d*U.p().xy + U.dv.y/U.d*U.p().xx,
                            U.E.yy*U.dv.x/U.d + 2.0*U.dv.y/U.d*U.p().xy,
                            U.E.zz*U.dv.x/U.d,
                            U.erot*U.dv.x/U.d));
}

/********************************************************
 * Gaussian2D_cState::Fy -- Solution flux (y-direction).*
 ********************************************************/

inline Gaussian2D_cState Gaussian2D_cState::Fy(void) {
  return (Gaussian2D_cState(dv.y,
                            dv.y*dv.x/d + p().xy,
                            dv.y*dv.y/d + p().yy, 
                            E.xx*dv.y/d + 2.0*dv.x/d*p().xy,
                            E.xy*dv.y/d + dv.y/d*p().xy + dv.x/d*p().yy,
                            E.yy*dv.y/d + 2.0*dv.y/d*p().yy,
                            E.zz*dv.y/d,
                            erot*dv.y/d));
}

inline Gaussian2D_cState Gaussian2D_cState::Fy(void) const {
  return (Gaussian2D_cState(dv.y,
                            dv.y*dv.x/d + p().xy,
                            dv.y*dv.y/d + p().yy, 
                            E.xx*dv.y/d + 2.0*dv.x/d*p().xy,
                            E.xy*dv.y/d + dv.y/d*p().xy + dv.x/d*p().yy,
                            E.yy*dv.y/d + 2.0*dv.y/d*p().yy,
                            E.zz*dv.y/d,
                            erot*dv.y/d));
}

inline Gaussian2D_cState Gaussian2D_cState::Fy(const Gaussian2D_cState &U) {
  return (Gaussian2D_cState(U.dv.y,
                            U.dv.y*U.dv.x/U.d + U.p().xy,
                            U.dv.y*U.dv.y/U.d + U.p().yy, 
                            U.E.xx*U.dv.y/U.d + 2.0*U.dv.x/U.d*U.p().xy,
                            U.E.xy*U.dv.y/U.d + U.dv.y/U.d*U.p().xy + U.dv.x/U.d*U.p().yy,
                            U.E.yy*U.dv.y/U.d + 2.0*U.dv.y/U.d*U.p().yy,
                            U.E.zz*U.dv.y/U.d,
                            U.erot*U.dv.y/U.d));
}

inline Gaussian2D_cState Fy(const Gaussian2D_cState &U) {
  return (Gaussian2D_cState(U.dv.y,
                            U.dv.y*U.dv.x/U.d + U.p().xy,
                            U.dv.y*U.dv.y/U.d + U.p().yy, 
                            U.E.xx*U.dv.y/U.d + 2.0*U.dv.x/U.d*U.p().xy,
                            U.E.xy*U.dv.y/U.d + U.dv.y/U.d*U.p().xy + U.dv.x/U.d*U.p().yy,
                            U.E.yy*U.dv.y/U.d + 2.0*U.dv.y/U.d*U.p().yy,
                            U.E.zz*U.dv.y/U.d,
                            U.erot*U.dv.y/U.d));
}

/********************************************************
 * Euler2D_cState::Fn -- Solution flux (n-direction).   *
 ********************************************************/
/*
inline Euler2D_cState Euler2D_cState::Fn(void) {
  return (Euler2D_cState(dv.x, sqr(dv.x)/d + p(),
                         dv.x*dv.y/d, dv.x*H()/d));
}

inline Euler2D_cState Euler2D_cState::Fn(void) const {
  return (Euler2D_cState(dv.x, sqr(dv.x)/d + p(),
                         dv.x*dv.y/d, dv.x*H()/d));
}

inline Euler2D_cState Euler2D_cState::Fn(const Euler2D_cState &U) {
  return (Euler2D_cState(U.dv.x, sqr(U.dv.x)/U.d + U.p(),
                         U.dv.x*U.dv.y/U.d, U.dv.x*U.H()/U.d));
}

inline Euler2D_cState Fn(const Euler2D_cState &U) {
  return (Euler2D_cState(U.dv.x, sqr(U.dv.x)/U.d + U.p(),
                         U.dv.x*U.dv.y/U.d, U.dv.x*U.H()/U.d));
}

inline void Euler2D_cState::dFndU(DenseMatrix &dFndU) {
  W().dFndU(dFndU);
}

inline void Euler2D_cState::dFndU(DenseMatrix &dFndU) const {
  W().dFndU(dFndU);
}

inline void dFndU(DenseMatrix &dFndU, const Euler2D_cState &U) {
  U.W().dFndU(dFndU);
}
*/
/********************************************************
 * Useful 2D Gaussian2D_cState::relax                   *
 ********************************************************/

inline void Gaussian2D_pState::relax(double deltat, int stage, const Gaussian2D_pState &W) {
  return;
  double tau_trans, tau_rot;
  double a, b, c;
  double pxx2(0.0), pxy2(0.0), pyy2(0.0), pzz2(0.0), erot2(0.0);
  double implicit_coeff1(0.0), implicit_coeff2(0.0), implicit_coeff3(0.0);
  double implicit_coeff4(0.0), implicit_coeff5(0.0), implicit_coeff6(0.0);

  tau_trans = viscosity()/pressure();
  tau_rot = 15.0/4.0*bulk_viscosity()/pressure();

  if(atoms==GAUSSIAN_MONATOMIC){
    a = 3.0*(double)stage*tau_trans+deltat;
    b = deltat;
    c = 3.0*(deltat+(double)stage*tau_trans);

    implicit_coeff1 = a/c;
    implicit_coeff2 = b/c;
    implicit_coeff3 = 0.0;
    implicit_coeff4 = (double)stage*tau_trans/(deltat+(double)stage*tau_trans);
    implicit_coeff5 = 0.0;
    implicit_coeff6 = 0.0;
  }else{
    a = 15.0*(double)stage*(double)stage*tau_trans*tau_rot+13.0*(double)stage*tau_trans*deltat;
    a = a+5.0*deltat*tau_rot*(double)stage+3.0*deltat*deltat;
    b = deltat*(5.0*(double)stage*tau_rot+3.0*deltat-2.0*(double)stage*tau_trans);
    c = 15.0*((double)stage*tau_trans+deltat)*((double)stage*tau_rot+deltat);

    implicit_coeff1 = a/c;
    implicit_coeff2 = b/c;
    implicit_coeff3 = 2.0*deltat/(5.0*(deltat+(double)stage*tau_rot));
    implicit_coeff4 = (double)stage*tau_trans/(deltat+(double)stage*tau_trans);
    implicit_coeff5 = deltat/(5.0*(deltat+(double)stage*tau_rot));
    implicit_coeff6 = (2.0*deltat+5.0*(double)stage*tau_rot)/(5.0*(deltat+(double)stage*tau_rot));
  }

  if(stage == ONE){

    pxx2  = implicit_coeff1*p.xx+implicit_coeff2*p.yy+implicit_coeff2*p.zz+implicit_coeff3*erot;
    pxy2  = implicit_coeff4*p.xy;
    pyy2  = implicit_coeff1*p.yy+implicit_coeff2*p.xx+implicit_coeff2*p.zz+implicit_coeff3*erot;
    pzz2  = implicit_coeff1*p.zz+implicit_coeff2*p.xx+implicit_coeff2*p.yy+implicit_coeff3*erot;
    erot2 = implicit_coeff5*(p.xx+p.yy+p.zz)+implicit_coeff6*erot;

  }else{

    p.xx = p.xx-deltat*(2.0*W.p.xx-W.p.yy-W.p.zz)/(3.0*tau_trans)/2.0;
    p.xy = p.xy-deltat*W.p.xy/tau_trans/2.0;
    p.yy = p.yy-deltat*(2.0*W.p.yy-W.p.xx-W.p.zz)/(3.0*tau_trans)/2.0;
    p.zz = p.zz-deltat*(2.0*W.p.zz-W.p.xx-W.p.yy)/(3.0*tau_trans)/2.0;

    if(atoms == GAUSSIAN_DIATOMIC){

      p.xx = p.xx-2.0*deltat*(W.p.xx+W.p.yy+W.p.zz-3.0*W.erot)/(15.0*tau_rot)/2.0;
      p.yy = p.yy-2.0*deltat*(W.p.xx+W.p.yy+W.p.zz-3.0*W.erot)/(15.0*tau_rot)/2.0;
      p.zz = p.zz-2.0*deltat*(W.p.xx+W.p.yy+W.p.zz-3.0*W.erot)/(15.0*tau_rot)/2.0;
      erot = erot-deltat*(3.0*W.erot-W.p.xx-W.p.yy-W.p.zz)/(5.0*tau_rot)/2.0;
    }

    pxx2  = implicit_coeff1*p.xx+implicit_coeff2*p.yy+implicit_coeff2*p.zz+implicit_coeff3*erot;
    pxy2  = implicit_coeff4*p.xy;
    pyy2  = implicit_coeff1*p.yy+implicit_coeff2*p.xx+implicit_coeff2*p.zz+implicit_coeff3*erot;
    pzz2  = implicit_coeff1*p.zz+implicit_coeff2*p.xx+implicit_coeff2*p.yy+implicit_coeff3*erot;
    erot2 = implicit_coeff5*(p.xx+p.yy+p.zz)+implicit_coeff6*erot;

  }
  
  p.xx = pxx2;
  p.xy = pxy2;
  p.yy = pyy2;
  p.zz = pzz2;
  erot = erot2;

}

/********************************************************
 * Useful 2D Gaussian state constants.                  *
 ********************************************************/

const Gaussian2D_pState Gaussian2D_W_STDATM(DENSITY_STDATM,
				      Vector2D_ZERO, PRESSURE_STDATM);
const Gaussian2D_pState Gaussian2D_W_VACUUM(ZERO, Vector2D_ZERO, ZERO);
const Gaussian2D_pState Gaussian2D_W_ONE(ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE);
const Gaussian2D_cState Gaussian2D_U_STDATM(Gaussian2D_W_STDATM);
const Gaussian2D_cState Gaussian2D_U_VACUUM(Gaussian2D_W_VACUUM);
const Gaussian2D_cState Gaussian2D_U_ONE(ONE, ONE, ONE, ONE, ONE, ONE, ONE , ONE);

/********************************************************
 * Euler2DState -- External subroutines.                *
 ********************************************************/

extern Gaussian2D_pState RoeAverage(const Gaussian2D_pState &Wl,
	      	                 const Gaussian2D_pState &Wr);

extern Gaussian2D_pState Rotate(const Gaussian2D_pState &W,
	      	                const Vector2D &norm_dir);

extern Gaussian2D_pState Reflect(const Gaussian2D_pState &W,
	      	                 const Vector2D &norm_dir);

extern Gaussian2D_pState NoSlip(const Gaussian2D_pState &W,
	      	                const Vector2D &norm_dir);

extern Gaussian2D_pState Adiabatic_Wall(const Gaussian2D_pState &W,
					const Gaussian2D_pState &Wo,
					const Vector2D &norm_dir);

extern Gaussian2D_pState BC_Characteristic_Pressure(const Gaussian2D_pState &Wi,
                                                    const Gaussian2D_pState &Wo,
	      	                                    const Vector2D &norm_dir);

extern Gaussian2D_pState BC_Characteristic_Velocity(const Gaussian2D_pState &Wi,
                                                    const Gaussian2D_pState &Wo,
						    const Vector2D &norm_dir);

extern Gaussian2D_pState BCs(const Gaussian2D_pState &Wb,
			  const Gaussian2D_pState &Wi,
			  const Gaussian2D_pState &dWdx,
			  const double &dx,
			  const int BC_type,
			  const int End_type);

extern Gaussian2D_pState BCs_x(const Gaussian2D_pState &Wb,
			    const Gaussian2D_pState &Wi,
			    const Gaussian2D_pState &dWdx,
			    const double &dx,
			    const int BC_type,
			    const int End_type);

extern Gaussian2D_pState BCs_y(const Gaussian2D_pState &Wb,
			    const Gaussian2D_pState &Wi,
			    const Gaussian2D_pState &dWdy,
			    const double &dy,
			    const int BC_type,
			    const int End_type);

extern Gaussian2D_pState WaveSpeedPos(const Gaussian2D_pState &lambda_a,
                                   const Gaussian2D_pState &lambda_l,
                                   const Gaussian2D_pState &lambda_r);

extern Gaussian2D_pState WaveSpeedNeg(const Gaussian2D_pState &lambda_a,
                                   const Gaussian2D_pState &lambda_l,
                                   const Gaussian2D_pState &lambda_r);

extern Gaussian2D_pState WaveSpeedAbs(const Gaussian2D_pState &lambda_a,
                                   const Gaussian2D_pState &lambda_l,
                                   const Gaussian2D_pState &lambda_r);

extern Gaussian2D_pState HartenFixPos(const Gaussian2D_pState &lambda_a,
                                   const Gaussian2D_pState &lambda_l,
                                   const Gaussian2D_pState &lambda_r);

extern Gaussian2D_pState HartenFixNeg(const Gaussian2D_pState &lambda_a,
                                   const Gaussian2D_pState &lambda_l,
                                   const Gaussian2D_pState &lambda_r);

extern Gaussian2D_pState HartenFixAbs(const Gaussian2D_pState &lambda_a,
                                   const Gaussian2D_pState &lambda_l,
                                   const Gaussian2D_pState &lambda_r);

extern Gaussian2D_cState FluxRoe(const Gaussian2D_pState &Wl,
	      	              const Gaussian2D_pState &Wr);

extern Gaussian2D_cState FluxRoe(const Gaussian2D_cState &Ul,
	      	              const Gaussian2D_cState &Ur);

extern Gaussian2D_cState FluxRoe_x(const Gaussian2D_pState &Wl,
	      	                const Gaussian2D_pState &Wr);

extern Gaussian2D_cState FluxRoe_x(const Gaussian2D_cState &Ul,
	      	                const Gaussian2D_cState &Ur);

extern Gaussian2D_cState FluxRoe_y(const Gaussian2D_pState &Wl,
	      	                const Gaussian2D_pState &Wr);

extern Gaussian2D_cState FluxRoe_y(const Gaussian2D_cState &Ul,
	      	                const Gaussian2D_cState &Ur);

extern Gaussian2D_cState FluxRoe_n(const Gaussian2D_pState &Wl,
	      	                const Gaussian2D_pState &Wr,
                                const Vector2D &norm_dir);

extern Gaussian2D_cState FluxRoe_n(const Gaussian2D_cState &Ul,
	      	                const Gaussian2D_cState &Ur,
                                const Vector2D &norm_dir);

extern Gaussian2D_cState FluxHLLE(const Gaussian2D_pState &Wl,
	      	               const Gaussian2D_pState &Wr);

extern Gaussian2D_cState FluxHLLE(const Gaussian2D_cState &Ul,
	      	               const Gaussian2D_cState &Ur);

extern Gaussian2D_cState FluxHLLE_x(const Gaussian2D_pState &Wl,
	      	                 const Gaussian2D_pState &Wr);

extern Gaussian2D_cState FluxHLLE_x(const Gaussian2D_cState &Ul,
	      	                 const Gaussian2D_cState &Ur);

extern Gaussian2D_cState FluxHLLE_y(const Gaussian2D_pState &Wl,
	      	                 const Gaussian2D_pState &Wr);

extern Gaussian2D_cState FluxHLLE_y(const Gaussian2D_cState &Ul,
	      	                 const Gaussian2D_cState &Ur);
  
extern Gaussian2D_cState FluxHLLE_n(const Gaussian2D_pState &Wl,
	      	                 const Gaussian2D_pState &Wr,
                                 const Vector2D &norm_dir);

extern Gaussian2D_cState FluxHLLE_n(const Gaussian2D_cState &Ul,
	      	                 const Gaussian2D_cState &Ur,
                                 const Vector2D &norm_dir);

extern Vector2D HLLE_wavespeeds(const Gaussian2D_pState &Wl,
                                const Gaussian2D_pState &Wr,
                                const Vector2D &norm_dir);
/*
extern Gaussian2D_cState FluxHLLC(const Gaussian2D_pState &Wl,
	      	               const Gaussian2D_pState &Wr);

extern Gaussian2D_cState FluxHLLC(const Gaussian2D_cState &Ul,
	      	               const Gaussian2D_cState &Ur);

extern Gaussian2D_cState FluxHLLC_x(const Gaussian2D_pState &Wl,
	      	                 const Gaussian2D_pState &Wr);

extern Gaussian2D_cState FluxHLLC_x(const Gaussian2D_cState &Ul,
	      	                 const Gaussian2D_cState &Ur);

extern Gaussian2D_cState FluxHLLC_y(const Gaussian2D_pState &Wl,
	      	                 const Gaussian2D_pState &Wr);

extern Gaussian2D_cState FluxHLLC_y(const Gaussian2D_cState &Ul,
	      	                 const Gaussian2D_cState &Ur);

extern Gaussian2D_cState FluxHLLC_n(const Gaussian2D_pState &Wl,
	      	                 const Gaussian2D_pState &Wr,
                                 const Vector2D &norm_dir);

extern Gaussian2D_cState FluxHLLC_n(const Gaussian2D_cState &Ul,
	      	                 const Gaussian2D_cState &Ur,
                                 const Vector2D &norm_dir);
*/

extern Gaussian2D_cState FluxKinetic_x(const Gaussian2D_pState &W1,
				       const Gaussian2D_pState &W2);

extern Gaussian2D_pState FlatPlate(const Gaussian2D_pState &Winf,
				   const Vector2D X,
				   double &eta,
				   double &f,
				   double &fp,
				   double &fpp);

//Test function
Gaussian2D_pState GaussianTestFunction (const double x);

#endif /* _GAUSSIAN2D_STATE_INCLUDED  */