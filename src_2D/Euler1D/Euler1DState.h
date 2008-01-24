/* Euler1DState.h:  Header file defining 1D Euler Solution State Classes. */

#ifndef _EULER1D_STATE_INCLUDED
#define _EULER1D_STATE_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

/* Include math macro, CFD, and gas constant header files. */

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _GAS_CONSTANTS_INCLUDED
#include "../Physics/GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED

/* Define the classes. */

#define	NUM_VAR_EULER1D    3

class Euler1D_cState;

/********************************************************
 * Class: Euler1D_pState                                *
 *                                                      *
 * Member functions                                     *
 *      d       -- Return density.                      *
 *      v       -- Return flow velocity.                *
 *      p       -- Return pressure.                     *
 *      g       -- Return specific heat ratio.          *
 *      gm1     -- Return g-1.                          *
 *      gm1i    -- Return 1/(g-1).                      *
 *      R       -- Return gas constant.                 *
 *      setgas  -- Set gas constants.                   *
 *      T       -- Return temperature.                  *
 *      e       -- Return specific internal energy.     *
 *      E       -- Return total energy.                 *
 *      h       -- Return specific enthalpy.            *
 *      H       -- Return total enthalpy.               *
 *      a       -- Return sound speed.                  *
 *      a2      -- Return sound speed square.           *
 *      M       -- Return Mach number.                  *
 *      s       -- Return specific entropy.             *
 *      dv      -- Return momentum.                     *
 *      To      -- Return stagnation temperature.       *
 *      po      -- Return stagnation pressure.          *
 *      ao      -- Return stagnation sound speed.       *
 *      ho      -- Return stagnation enthalpy.          * 
 *      U       -- Return conserved solution state.     *
 *      F       -- Return solution flux.                *
 *      C       -- Return characteristic variable       *
 *                 solution state.                      *
 *      CtoW    -- Return primitive solution state      *
 *                 given the characteristic variables.  *
 *      lambda  -- Return eigenvalue(s).                *
 *      rp      -- Return primitive right eigenvector.  *
 *      rc      -- Return conserved right eigenvector.  *
 *      lp      -- Return primitive left eigenvectors.  *
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
class Euler1D_pState{
  private:
  public:
    double          d;   // Density.
    double          v;   // Flow velocity.
    double          p;   // Pressure.
    static double   g;   // Specific heat ratio.
    static double gm1;   // g-1
    static double gm1i;  // 1/(g-1)
    static double   R;   // Gas constant.
	                 // Made public so can access them.

    static const int NumberOfVariables = NUM_VAR_EULER1D;
		      
    /* Creation, copy, and assignment constructors. */
    Euler1D_pState(void) {
       d = DENSITY_STDATM; v = ZERO; p = PRESSURE_STDATM;
    }

    Euler1D_pState(const Euler1D_pState &W) {
       d = W.d; v = W.v; p = W.p;
    }

    Euler1D_pState(const Euler1D_cState &U);

    Euler1D_pState(const double &rho,
	           const double &u,
	           const double &pre) {
       d = rho; v = u; p = pre;
    }
    
    Euler1D_pState(const double & Val){
      d = Val; v = Val; p = Val;
    }
  
    /* Vacuum operator. */
    void Vacuum() {
      d = ZERO; v = ZERO; p = ZERO;
    }

    /* One operator. */
    void One() {
      d = ONE; v = ONE; p = ONE;
    }
   
    /* Destructor. */
    // ~Euler1D_pState(void);
    // Use automatically generated destructor.

    /* Set gas constants. */
    void setgas(void);
    void setgas(char *string_ptr);

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

    /* Total enthalpy. */
    double H(void);
    double H(void) const;

    /* Sound speed. */
    double a(void);
    double a(void) const;

    /* Sound speed squared. */
    double a2(void);
    double a2(void) const;

    /* Mach number. */
    double M(void);
    double M(void) const;

    /* Specific entropy. */
    double s(void);
    double s(void) const;

    /* Momentum. */
    double dv(void);
    double dv(void) const;

    /* Stagnation temperature. */
    double To(void);
    double To(void) const;

    /* Stagnation pressure. */
    double po(void);
    double po(void) const;

    /* Stagnation sound speed. */
    double ao(void);
    double ao(void) const;

    /* Stagnation enthalpy. */
    double ho(void);
    double ho(void) const;

    /* Conserved solution state. */
    Euler1D_cState U(void);
    Euler1D_cState U(void) const;
    Euler1D_cState U(const Euler1D_pState &W);
    friend Euler1D_cState U(const Euler1D_pState &W);

    /* Solution flux. */
    Euler1D_cState F(void);
    Euler1D_cState F(void) const;
    Euler1D_cState F(const Euler1D_pState &W);
    friend Euler1D_cState F(const Euler1D_pState &W);

    /* Characteristic variables. */
    Euler1D_pState C(void);
    Euler1D_pState C(void) const;
    Euler1D_pState C(const Euler1D_pState &W);
    friend Euler1D_pState C(const Euler1D_pState &W);
    friend Euler1D_pState CtoW(const Euler1D_pState &W);

    /* Functions for transforming to the locally defined characteristic variables from 
       the conserved variables.
       The left and the right eigenvectors used for the transformation are based on the values of 
       a reference state. */
    Euler1D_cState CharactVarToConservedVar(const double SoundSpeed, const double Vel);
    Euler1D_cState CharactVarToConservedVar(const Euler1D_pState &W_Ref);
    Euler1D_cState CharactVarToConservedVar(const Euler1D_cState &U_Ref);

    /* Eigenvalue(s). */
    Euler1D_pState lambda(void);
    Euler1D_pState lambda(void) const;
    Euler1D_pState lambda(const Euler1D_pState &W);
    friend Euler1D_pState lambda(const Euler1D_pState &W);
    double lambda(int index);
    double lambda(int index) const;
    friend double lambda(const Euler1D_pState &W, int index);

    /* Primitive right eigenvector. */
    Euler1D_pState rp(int index);
    Euler1D_pState rp(int index) const;
    friend Euler1D_pState rp(const Euler1D_pState &W, int index);

    /* Conserved right eigenvector. */
    Euler1D_cState rc(int index);
    Euler1D_cState rc(int index) const;
    friend Euler1D_cState rc(const Euler1D_pState &W, int index);

    /* Primitive left eigenvector. */
    Euler1D_pState lp(int index);
    Euler1D_pState lp(int index) const;
    friend Euler1D_pState lp(const Euler1D_pState &W, int index);
    
    /* Assignment operator. */
    // Euler1D_pState operator = (const Euler1D_pState &W);
    // Use automatically generated assignment operator.

    /* Index operator. */
    double &operator[](int index) {
      assert( index >= 1 && index <= NUM_VAR_EULER1D );
      switch(index) {
        case 1 :
	  return (d);
        case 2 :
	  return (v);
        case 3 :
	  return (p);
        default:
	  return (d);
      };
    }
    
    const double &operator[](int index) const {
      assert( index >= 1 && index <= NUM_VAR_EULER1D );
      switch(index) {
        case 1 :
	  return (d);
        case 2 :
	  return (v);
        case 3 :
	  return (p);
        default:
	  return (d);
      };
    }

    /* Binary arithmetic operators. */
    friend Euler1D_pState operator +(const Euler1D_pState &W1, const Euler1D_pState &W2);
    friend Euler1D_pState operator -(const Euler1D_pState &W1, const Euler1D_pState &W2);
    friend double operator *(const Euler1D_pState &W1, const Euler1D_pState &W2);
    friend Euler1D_pState operator *(const Euler1D_pState &W, const double &a);
    friend Euler1D_pState operator *(const double &a, const Euler1D_pState &W);
    friend Euler1D_pState operator /(const Euler1D_pState &W, const double &a);
    friend Euler1D_pState operator ^(const Euler1D_pState &W1, const Euler1D_pState &W2);
    
    /* Unary arithmetic operators. */
    friend Euler1D_pState operator +(const Euler1D_pState &W);
    friend Euler1D_pState operator -(const Euler1D_pState &W);

    /* Shortcut arithmetic operators. */
    friend Euler1D_pState &operator +=(Euler1D_pState &W1, const Euler1D_pState &W2);
    friend Euler1D_pState &operator -=(Euler1D_pState &W1, const Euler1D_pState &W2);
    
    /* Relational operators. */
    friend int operator ==(const Euler1D_pState &W1, const Euler1D_pState &W2);
    friend int operator !=(const Euler1D_pState &W1, const Euler1D_pState &W2);
    friend int operator <=(const Euler1D_pState &W1, const Euler1D_pState &W2);
    friend int operator >=(const Euler1D_pState &W1, const Euler1D_pState &W2);
    friend int operator <(const Euler1D_pState &W1, const Euler1D_pState &W2);
    friend int operator >(const Euler1D_pState &W1, const Euler1D_pState &W2);

    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const Euler1D_pState &W);
    friend istream &operator >> (istream &in_file,  Euler1D_pState &W);
    
};

/********************************************************
 * Class: Euler1D_cState                                *
 *                                                      *
 * Member functions                                     *
 *      d       -- Return density.                      *
 *      dv      -- Return momentum.                     *
 *      E       -- Return total energy.                 *
 *      g       -- Return specific heat ratio.          *
 *      gm1     -- Return g-1.                          *
 *      gm1i    -- Return 1/(g-1).                      *
 *      R       -- Return gas constant.                 *
 *      setgas  -- Set gas constants.                   *
 *      v       -- Return flow velocity.                *
 *      p       -- Return pressure.                     *
 *      T       -- Return temperature.                  *
 *      e       -- Return specific internal energy.     *
 *      h       -- Return specific enthalpy.            *
 *      H       -- Return total enthalpy.               *
 *      a       -- Return sound speed.                  *
 *      a2      -- Return sound speed square.           *
 *      M       -- Return Mach number.                  *
 *      s       -- Return specific entropy.             *
 *      To      -- Return stagnation temperature.       *
 *      po      -- Return stagnation pressure.          *
 *      ao      -- Return stagnation sound speed.       *
 *      ho      -- Return stagnation enthalpy.          *
 *      W       -- Return primitive solution state.     *
 *      F       -- Return solution flux.                *
 *      C       -- Return characteristic variable       *
 *                 solution state.                      *
 *      CtoU    -- Return conserved solution state      *
 *                 given the characteristic variables.  *
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
class Euler1D_cState{
  private:
  public:
    double          d;   // Density.
    double         dv;   // Momentum.
    double          E;   // Total Energy.
    static double   g;   // Specific heat ratio.
    static double gm1;   // g-1
    static double gm1i;  // 1/(g-1)
    static double   R;   // Gas constant.
	                 // Made public so can access them.

    static const int NumberOfVariables = NUM_VAR_EULER1D;
		      
    /* Creation, copy, and assignment constructors. */
    Euler1D_cState(void) {
       d = DENSITY_STDATM; dv = ZERO; E = PRESSURE_STDATM/(GAMMA_AIR-ONE);
    }

    Euler1D_cState(const Euler1D_cState &U) {
       d = U.d; dv = U.dv; E = U.E;
    }

    Euler1D_cState(const Euler1D_pState &W);

    Euler1D_cState(const double &rho,
	           const double &rhou,
	           const double &Etotal) {
       d = rho; dv = rhou; E = Etotal;
    }
    
    /* Destructor. */
    // ~Euler1D_cState(void);
    // Use automatically generated destructor.

    /* Set gas constants. */
    void setgas(void);
    void setgas(char *string_ptr);

    /* Flow velocity. */
    double v(void);
    double v(void) const;

    /* Pressure. */
    double p(void);
    double p(void) const;

    /* Temperature. */
    double T(void);
    double T(void) const;

    /* Specific internal energy. */
    double e(void);
    double e(void) const;

    /* Specific enthalpy. */
    double h(void);
    double h(void) const;

    /* Total enthalpy. */
    double H(void);
    double H(void) const;

    /* Sound speed. */
    double a(void);
    double a(void) const;

    /* Sound speed squared. */
    double a2(void);
    double a2(void) const;

    /* Mach number. */
    double M(void);
    double M(void) const;

    /* Specific entropy. */
    double s(void);
    double s(void) const;

    /* Stagnation temperature. */
    double To(void);
    double To(void) const;

    /* Stagnation pressure. */
    double po(void);
    double po(void) const;

    /* Stagnation sound speed. */
    double ao(void);
    double ao(void) const;

    /* Stagnation enthalpy. */
    double ho(void);
    double ho(void) const;

    /* Primitive solution state. */
    Euler1D_pState W(void);
    Euler1D_pState W(void) const;
    Euler1D_pState W(const Euler1D_cState &U);
    friend Euler1D_pState W(const Euler1D_cState &U);
    
    /* Solution flux. */
    Euler1D_cState F(void);
    Euler1D_cState F(void) const;
    Euler1D_cState F(const Euler1D_cState &U);
    friend Euler1D_cState F(const Euler1D_cState &U);

   /* Characteristic variables. */
    Euler1D_pState C(void);
    Euler1D_pState C(void) const;
    Euler1D_pState C(const Euler1D_cState &U);
    friend Euler1D_pState C(const Euler1D_cState &U);
    friend Euler1D_cState CtoU(const Euler1D_pState &W);
    
   /* Characteristic variable using the left and right conserved eigenvectors */
    Euler1D_pState ConservedVarToCharactVar(const double SoundSpeed, const double Vel) const;
    Euler1D_pState ConservedVarToCharactVar(const double SoundSpeed, const double Vel);
    Euler1D_pState ConservedVarToCharactVar(void);
    Euler1D_pState ConservedVarToCharactVar(void) const;
   /* the transformation below is relative to the variables of other state */
    Euler1D_pState ConservedVarToCharactVar(const Euler1D_cState &U);
    Euler1D_pState ConservedVarToCharactVar(const Euler1D_cState &U) const;
    Euler1D_pState ConservedVarToCharactVar(const Euler1D_pState &W);
    Euler1D_pState ConservedVarToCharactVar(const Euler1D_pState &W) const;
    friend Euler1D_cState CharactVarToConservedVar(const Euler1D_pState &Wc, const Euler1D_cState &U_Ref);
     
    /* Assignment operator. */
    // Euler1D_cState operator = (const Euler1D_cState &U);
    // Use automatically generated assignment operator.

    /* Index operator. */
    double &operator[](int index) {
      assert( index >= 1 && index <= NUM_VAR_EULER1D );
      switch(index) {
        case 1 :
	  return (d);
        case 2 :
	  return (dv);
        case 3 :
	  return (E);
        default:
	  return (d);
      };
    }
    
    const double &operator[](int index) const {
      assert( index >= 1 && index <= NUM_VAR_EULER1D );
      switch(index) {
        case 1 :
	  return (d);
        case 2 :
	  return (dv);
        case 3 :
	  return (E);
        default:
	  return (d);
      };
    }

    /* Binary arithmetic operators. */
    friend Euler1D_cState operator +(const Euler1D_cState &U1, const Euler1D_cState &U2);
    friend Euler1D_cState operator -(const Euler1D_cState &U1, const Euler1D_cState &U2);
    friend double operator *(const Euler1D_cState &U1, const Euler1D_cState &U2);
    friend Euler1D_cState operator *(const Euler1D_cState &U, const double &a);
    friend Euler1D_cState operator *(const double &a, const Euler1D_cState &U);
    friend Euler1D_cState operator /(const Euler1D_cState &U, const double &a);
    friend Euler1D_cState operator ^(const Euler1D_cState &U1, const Euler1D_cState &U2);
    
    /* Unary arithmetic operators. */
    friend Euler1D_cState operator +(const Euler1D_cState &U);
    friend Euler1D_cState operator -(const Euler1D_cState &U);

    /* Shortcut arithmetic operators. */
    friend Euler1D_cState &operator +=(Euler1D_cState &U1, const Euler1D_cState &U2);
    friend Euler1D_cState &operator -=(Euler1D_cState &U1, const Euler1D_cState &U2);
    
    /* Relational operators. */
    friend int operator ==(const Euler1D_cState &U1, const Euler1D_cState &U2);
    friend int operator !=(const Euler1D_cState &U1, const Euler1D_cState &U2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const Euler1D_cState &U);
    friend istream &operator >> (istream &in_file,  Euler1D_cState &U);
    
};

/********************************************************
 * Euler1D_pState::setgas -- Assign gas constants.      *
 ********************************************************/
inline void Euler1D_pState::setgas(void) {
    g = GAMMA_AIR;
    R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    gm1 = g - ONE;
    gm1i = ONE/gm1;
}

inline void Euler1D_pState::setgas(char *string_ptr) {
   if (strcmp(string_ptr, "AIR") == 0) {
     g = GAMMA_AIR;
     R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
   } else if (strcmp(string_ptr, "A") == 0) {
     g = GAMMA_A;
     R = R_UNIVERSAL/(MOLE_WT_A*MILLI);
   } else if (strcmp(string_ptr, "CO") == 0) {
     g = GAMMA_CO;
     R = R_UNIVERSAL/(MOLE_WT_CO*MILLI);
   } else if (strcmp(string_ptr, "CO2") == 0) {
     g = GAMMA_CO2;
     R = R_UNIVERSAL/(MOLE_WT_CO2*MILLI);
   } else if (strcmp(string_ptr, "CH4") == 0) {
     g = GAMMA_CH4;
     R = R_UNIVERSAL/(MOLE_WT_CH4*MILLI);
   } else if (strcmp(string_ptr, "H") == 0) {
     g = GAMMA_H;
     R = R_UNIVERSAL/(MOLE_WT_H*MILLI);
   } else if (strcmp(string_ptr, "H2") == 0) {
     g = GAMMA_H2;
     R = R_UNIVERSAL/(MOLE_WT_H2*MILLI);
   } else if (strcmp(string_ptr, "HE") == 0) {
     g = GAMMA_HE;
     R = R_UNIVERSAL/(MOLE_WT_HE*MILLI);
   } else if (strcmp(string_ptr, "H2O") == 0) {
     g = GAMMA_H2O;
     R = R_UNIVERSAL/(MOLE_WT_H2O*MILLI);
   } else if (strcmp(string_ptr, "N2") == 0) {
     g = GAMMA_N2;
     R = R_UNIVERSAL/(MOLE_WT_N2*MILLI);
   } else if (strcmp(string_ptr, "O") == 0) {
     g = GAMMA_O;
     R = R_UNIVERSAL/(MOLE_WT_O*MILLI);
   } else if (strcmp(string_ptr, "O2") == 0) {
     g = GAMMA_O2;
     R = R_UNIVERSAL/(MOLE_WT_O2*MILLI);
   } else if (strcmp(string_ptr, "e") == 0) {
     g = GAMMA_e;
     R = R_UNIVERSAL/(MOLE_WT_e*MILLI);
   } else if (strcmp(string_ptr, "Zb") == 0) {
     g = GAMMA_ZB;
     R = R_UNIVERSAL/(MOLE_WT_ZB*MILLI);
   } else {
     g = GAMMA_AIR;
     R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
   } /* endif */
   gm1 = g - ONE;
   gm1i = ONE/gm1;
}

/********************************************************
 * Euler1D_pState::T -- Temperature.                    *
 ********************************************************/
inline double Euler1D_pState::T(void) {
    return (p/(d*R));
}

inline double Euler1D_pState::T(void) const {
    return (p/(d*R));
}

/********************************************************
 * Euler1D_pState::e -- Specific internal energy.       *
 ********************************************************/
inline double Euler1D_pState::e(void) {
    return (p/(gm1*d));
}

inline double Euler1D_pState::e(void) const {
    return (p/(gm1*d));
}

/********************************************************
 * Euler1D_pState::E -- Total energy.                   *
 ********************************************************/
inline double Euler1D_pState::E(void) {
    return (p*gm1i + HALF*d*sqr(v));
}

inline double Euler1D_pState::E(void) const {
    return (p*gm1i + HALF*d*sqr(v));
}

/********************************************************
 * Euler1D_pState::h -- Specific enthalpy.              *
 ********************************************************/
inline double Euler1D_pState::h(void) {
    return (g*p/(gm1*d) + HALF*sqr(v));
}

inline double Euler1D_pState::h(void) const {
    return (g*p/(gm1*d) + HALF*sqr(v));
}

/********************************************************
 * Euler1D_pState::H -- Total enthalpy.                 *
 ********************************************************/
inline double Euler1D_pState::H(void) {
    return (g*gm1i*p + HALF*d*sqr(v));
}

inline double Euler1D_pState::H(void) const {
    return (g*gm1i*p + HALF*d*sqr(v));
}

/********************************************************
 * Euler1D_pState::a -- Sound speed.                    *
 ********************************************************/
inline double Euler1D_pState::a(void) {
    return (sqrt(g*p/d));
}

inline double Euler1D_pState::a(void) const {
    return (sqrt(g*p/d));
}

/********************************************************
 * Euler1D_pState::a2 -- Sound speed squared.           *
 ********************************************************/
inline double Euler1D_pState::a2(void) {
    return (g*p/d);
}

inline double Euler1D_pState::a2(void) const {
    return (g*p/d);
}

/********************************************************
 * Euler1D_pState::M -- Mach number.                    *
 ********************************************************/
inline double Euler1D_pState::M(void) {
    return (v/sqrt(g*p/d));
}

inline double Euler1D_pState::M(void) const {
    return (v/sqrt(g*p/d));
}

/********************************************************
 * Euler1D_pState::s -- Specific entropy.               *
 ********************************************************/
inline double Euler1D_pState::s(void) {
    return (R*gm1i*log(p/pow(d, g)));
}

inline double Euler1D_pState::s(void) const {
    return (R*gm1i*log(p/pow(d, g)));
}

/********************************************************
 * Euler1D_pState::dv -- Momentum.                      *
 ********************************************************/
inline double Euler1D_pState::dv(void) {
    return (d*v);
}

inline double Euler1D_pState::dv(void) const {
    return (d*v);
}

/********************************************************
 * Euler1D_pState::To -- Stagnation temperature.        *
 ********************************************************/
inline double Euler1D_pState::To(void) {
    return ((p/(d*R))*(ONE+HALF*gm1*v*v/(g*p/d)));
}

inline double Euler1D_pState::To(void) const {
    return ((p/(d*R))*(ONE+HALF*gm1*v*v/(g*p/d)));
}

/********************************************************
 * Euler1D_pState::po -- Stagnation pressure.           *
 ********************************************************/
inline double Euler1D_pState::po(void) {
    return (p*pow(ONE+HALF*gm1*v*v/(g*p/d), g*gm1i));
}

inline double Euler1D_pState::po(void) const {
    return (p*pow(ONE+HALF*gm1*v*v/(g*p/d), g*gm1i));
}

/********************************************************
 * Euler1D_pState::ao -- Stagnation sound speed.        *
 ********************************************************/
inline double Euler1D_pState::ao(void) {
    return (sqrt((g*p/d)*(ONE+HALF*gm1*v*v/(g*p/d))));
}

inline double Euler1D_pState::ao(void) const {
    return (sqrt((g*p/d)*(ONE+HALF*gm1*v*v/(g*p/d))));
}

/********************************************************
 * Euler1D_pState::ho -- Stagnation enthalpy.           *
 ********************************************************/
inline double Euler1D_pState::ho(void) {
    return ((g*p/(gm1*d) + HALF*sqr(v))*(ONE+HALF*gm1*v*v/(g*p/d)));
}

inline double Euler1D_pState::ho(void) const {
    return ((g*p/(gm1*d) + HALF*sqr(v))*(ONE+HALF*gm1*v*v/(g*p/d)));
}

/********************************************************
 * Euler1D_pState -- Binary arithmetic operators.       *
 ********************************************************/
inline Euler1D_pState operator +(const Euler1D_pState &W1, const Euler1D_pState &W2) {
  return (Euler1D_pState(W1.d+W2.d,W1.v+W2.v,W1.p+W2.p));
}

inline Euler1D_pState operator -(const Euler1D_pState &W1, const Euler1D_pState &W2) {
  return (Euler1D_pState(W1.d-W2.d,W1.v-W2.v,W1.p-W2.p));
}

// Inner product operator.
inline double operator *(const Euler1D_pState &W1, const Euler1D_pState &W2) {
   return (W1.d*W2.d+W1.v*W2.v+W1.p*W2.p);
}

inline Euler1D_pState operator *(const Euler1D_pState &W, const double &a) {
  return (Euler1D_pState(a*W.d,a*W.v,a*W.p));
}

inline Euler1D_pState operator *(const double &a, const Euler1D_pState &W) {
  return (Euler1D_pState(a*W.d,a*W.v,a*W.p));
}

inline Euler1D_pState operator /(const Euler1D_pState &W, const double &a) {
  return (Euler1D_pState(W.d/a,W.v/a,W.p/a));
}

// My useful solution state product operator.
inline Euler1D_pState operator ^(const Euler1D_pState &W1, const Euler1D_pState &W2) {
   return (Euler1D_pState(W1.d*W2.d,W1.v*W2.v,W1.p*W2.p));
}

/********************************************************
 * Euler1D_pState -- Unary arithmetic operators.        *
 ********************************************************/
inline Euler1D_pState operator +(const Euler1D_pState &W) {
  return (Euler1D_pState(W.d,W.v,W.p));
}

inline Euler1D_pState operator -(const Euler1D_pState &W) {
  return (Euler1D_pState(-W.d,-W.v,-W.p));
}

/********************************************************
 * Euler1D_pState -- Shortcut arithmetic operators.     *
 ********************************************************/
inline Euler1D_pState &operator +=(Euler1D_pState &W1, const Euler1D_pState &W2) {
  W1.d += W2.d;
  W1.v += W2.v;
  W1.p += W2.p;
  return (W1);
}

inline Euler1D_pState &operator -=(Euler1D_pState &W1, const Euler1D_pState &W2) {
  W1.d -= W2.d;
  W1.v -= W2.v;
  W1.p -= W2.p;
  return (W1);
}

/********************************************************
 * Euler1D_pState -- Relational operators.              *
 ********************************************************/
inline int operator ==(const Euler1D_pState &W1, const Euler1D_pState &W2) {
  return (W1.d == W2.d && W1.v == W2.v && W1.p == W2.p);
}

inline int operator !=(const Euler1D_pState &W1, const Euler1D_pState &W2) {
  return (W1.d != W2.d || W1.v != W2.v || W1.p != W2.p);
}

inline int operator <=(const Euler1D_pState &W1, const Euler1D_pState &W2) {
  return (W1[1]<=W2[1] && W1[2]<=W2[2] && W1[3]<=W2[3]);
}

inline int operator >=(const Euler1D_pState &W1, const Euler1D_pState &W2) {
  return (W1[1]>=W2[1] && W1[2]>=W2[2] && W1[3]>=W2[3]);
}

inline int operator <(const Euler1D_pState &W1, const Euler1D_pState &W2) {
  return (W1[1]<W2[1] && W1[2]<W2[2] && W1[3]<W2[3]);
}

inline int operator >(const Euler1D_pState &W1, const Euler1D_pState &W2) {
  return (W1[1]>W2[1] && W1[2]>W2[2] && W1[3]>W2[3]);
}

/********************************************************
 * Euler1D_pState -- Input-output operators.            *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const Euler1D_pState &W) {
  out_file.setf(ios::scientific);
  out_file << " " << W.d << " " << W.v << " " << W.p;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, Euler1D_pState &W) {
  in_file.setf(ios::skipws);
  in_file >> W.d >> W.v >> W.p;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/********************************************************
 * Euler1D_cState::setgas -- Assign gas constants.      *
 ********************************************************/
inline void Euler1D_cState::setgas(void) {
    g = GAMMA_AIR;
    R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    gm1 = g - ONE;
    gm1i = ONE/gm1;
}

inline void Euler1D_cState::setgas(char *string_ptr) {
   if (strcmp(string_ptr, "AIR") == 0) {
     g = GAMMA_AIR;
     R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
   } else if (strcmp(string_ptr, "A") == 0) {
     g = GAMMA_A;
     R = R_UNIVERSAL/(MOLE_WT_A*MILLI);
   } else if (strcmp(string_ptr, "CO") == 0) {
     g = GAMMA_CO;
     R = R_UNIVERSAL/(MOLE_WT_CO*MILLI);
   } else if (strcmp(string_ptr, "CO2") == 0) {
     g = GAMMA_CO2;
     R = R_UNIVERSAL/(MOLE_WT_CO2*MILLI);
   } else if (strcmp(string_ptr, "CH4") == 0) {
     g = GAMMA_CH4;
     R = R_UNIVERSAL/(MOLE_WT_CH4*MILLI);
   } else if (strcmp(string_ptr, "H") == 0) {
     g = GAMMA_H;
     R = R_UNIVERSAL/(MOLE_WT_H*MILLI);
   } else if (strcmp(string_ptr, "H2") == 0) {
     g = GAMMA_H2;
     R = R_UNIVERSAL/(MOLE_WT_H2*MILLI);
   } else if (strcmp(string_ptr, "HE") == 0) {
     g = GAMMA_HE;
     R = R_UNIVERSAL/(MOLE_WT_HE*MILLI);
   } else if (strcmp(string_ptr, "H2O") == 0) {
     g = GAMMA_H2O;
     R = R_UNIVERSAL/(MOLE_WT_H2O*MILLI);
   } else if (strcmp(string_ptr, "N2") == 0) {
     g = GAMMA_N2;
     R = R_UNIVERSAL/(MOLE_WT_N2*MILLI);
   } else if (strcmp(string_ptr, "O") == 0) {
     g = GAMMA_O;
     R = R_UNIVERSAL/(MOLE_WT_O*MILLI);
   } else if (strcmp(string_ptr, "O2") == 0) {
     g = GAMMA_O2;
     R = R_UNIVERSAL/(MOLE_WT_O2*MILLI);
   } else if (strcmp(string_ptr, "e") == 0) {
     g = GAMMA_e;
     R = R_UNIVERSAL/(MOLE_WT_e*MILLI);
   } else if (strcmp(string_ptr, "Zb") == 0) {
     g = GAMMA_ZB;
     R = R_UNIVERSAL/(MOLE_WT_ZB*MILLI);
   } else {
     g = GAMMA_AIR;
     R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
   } /* endif */
   gm1 = g - ONE;
   gm1i = ONE/gm1;
}

/********************************************************
 * Euler1D_cState::v -- Flow velocity.                  *
 ********************************************************/
inline double Euler1D_cState::v(void) {
    return (dv/d);
}

inline double Euler1D_cState::v(void) const {
    return (dv/d);
}

/********************************************************
 * Euler1D_cState::p -- Pressure.                       *
 ********************************************************/
inline double Euler1D_cState::p(void) {
    return (gm1*(E - HALF*sqr(dv)/d));
}

inline double Euler1D_cState::p(void) const {
    return (gm1*(E - HALF*sqr(dv)/d));
}

/********************************************************
 * Euler1D_cState::T -- Temperature.                    *
 ********************************************************/
inline double Euler1D_cState::T(void) {
    return (gm1*(E - HALF*sqr(dv)/d)/(d*R));
}

inline double Euler1D_cState::T(void) const {
    return (gm1*(E - HALF*sqr(dv)/d)/(d*R));
}

/********************************************************
 * Euler1D_cState::e -- Specific internal energy.       *
 ********************************************************/
inline double Euler1D_cState::e(void) {
    return (E/d - HALF*sqr(dv/d));
}

inline double Euler1D_cState::e(void) const {
    return (E/d - HALF*sqr(dv/d));
}

/********************************************************
 * Euler1D_cState::h -- Specific enthalpy.              *
 ********************************************************/
inline double Euler1D_cState::h(void) {
    return (g*E/d - gm1*HALF*sqr(dv/d));
}

inline double Euler1D_cState::h(void) const {
    return (g*E/d - gm1*HALF*sqr(dv/d));
}

/********************************************************
 * Euler1D_cState::H -- Total enthalpy.                 *
 ********************************************************/
inline double Euler1D_cState::H(void) {
    return (g*E - gm1*HALF*sqr(dv)/d);
}

inline double Euler1D_cState::H(void) const {
    return (g*E - gm1*HALF*sqr(dv)/d);
}

/********************************************************
 * Euler1D_cState::a -- Sound speed.                    *
 ********************************************************/
inline double Euler1D_cState::a(void) {
    return (sqrt(g*gm1*(E/d - HALF*sqr(dv/d))));
}

inline double Euler1D_cState::a(void) const {
    return (sqrt(g*gm1*(E/d - HALF*sqr(dv/d))));
}

/********************************************************
 * Euler1D_cState::a2 -- Sound speed squared.           *
 ********************************************************/
inline double Euler1D_cState::a2(void) {
    return (g*gm1*(E/d - HALF*sqr(dv/d)));
}

inline double Euler1D_cState::a2(void) const {
    return (g*gm1*(E/d - HALF*sqr(dv/d)));
}

/********************************************************
 * Euler1D_cState::M -- Mach number.                    *
 ********************************************************/
inline double Euler1D_cState::M(void) {
    return (dv/(d*sqrt(g*gm1*(E/d - HALF*sqr(dv/d)))));
}

inline double Euler1D_cState::M(void) const {
    return (dv/(d*sqrt(g*gm1*(E/d - HALF*sqr(dv/d)))));
}

/********************************************************
 * Euler1D_cState::s -- Specific entropy.               *
 ********************************************************/
inline double Euler1D_cState::s(void) {
    return (R*gm1i*log(gm1*(E - HALF*sqr(dv)/d)/pow(d, g)));
}

inline double Euler1D_cState::s(void) const {
    return (R*gm1i*log(gm1*(E - HALF*sqr(dv)/d)/pow(d, g)));
}

/********************************************************
 * Euler1D_cState::To -- Stagnation temperature.        *
 ********************************************************/
inline double Euler1D_cState::To(void) {
    return ((gm1*(E - HALF*sqr(dv)/d)/(d*R))*
	    (ONE+HALF*gm1*dv*dv/(d*d*g*gm1*(E/d - HALF*sqr(dv/d)))));
}

inline double Euler1D_cState::To(void) const {
    return ((gm1*(E - HALF*sqr(dv)/d)/(d*R))*
	    (ONE+HALF*gm1*dv*dv/(d*d*g*gm1*(E/d - HALF*sqr(dv/d)))));
}

/********************************************************
 * Euler1D_cState::po -- Stagnation pressure.           *
 ********************************************************/
inline double Euler1D_cState::po(void) {
    return ((gm1*(E - HALF*sqr(dv)/d))*
	    pow(ONE+HALF*gm1*dv*dv/(d*d*g*gm1*(E/d - HALF*sqr(dv/d))), g*gm1i));
}

inline double Euler1D_cState::po(void) const {
    return ((gm1*(E - HALF*sqr(dv)/d))*
	    pow(ONE+HALF*gm1*dv*dv/(d*d*g*gm1*(E/d - HALF*sqr(dv/d))), g*gm1i));
}

/********************************************************
 * Euler1D_cState::ao -- Stagnation sound speed.        *
 ********************************************************/
inline double Euler1D_cState::ao(void) {
    return (sqrt((g*gm1*(E/d - HALF*sqr(dv/d)))*
	         (ONE+HALF*gm1*dv*dv/(d*d*g*gm1*(E/d - HALF*sqr(dv/d))))));
}

inline double Euler1D_cState::ao(void) const {
    return (sqrt((g*gm1*(E/d - HALF*sqr(dv/d)))*
	         (ONE+HALF*gm1*dv*dv/(d*d*g*gm1*(E/d - HALF*sqr(dv/d))))));
}

/********************************************************
 * Euler1D_cState::ho -- Stagnation enthalpy.           *
 ********************************************************/
inline double Euler1D_cState::ho(void) {
    return ((g*E/d - gm1*HALF*sqr(dv/d))*
	    (ONE+HALF*gm1*dv*dv/(d*d*g*gm1*(E/d - HALF*sqr(dv/d)))));
}

inline double Euler1D_cState::ho(void) const {
    return ((g*E/d - gm1*HALF*sqr(dv/d))*
	    (ONE+HALF*gm1*dv*dv/(d*d*g*gm1*(E/d - HALF*sqr(dv/d)))));
}

/********************************************************
 * Euler1D_cState -- Binary arithmetic operators.       *
 ********************************************************/
inline Euler1D_cState operator +(const Euler1D_cState &U1, const Euler1D_cState &U2) {
  return (Euler1D_cState(U1.d+U2.d,U1.dv+U2.dv,U1.E+U2.E));
}

inline Euler1D_cState operator -(const Euler1D_cState &U1, const Euler1D_cState &U2) {
  return (Euler1D_cState(U1.d-U2.d,U1.dv-U2.dv,U1.E-U2.E));
}

// Inner product operator.
inline double operator *(const Euler1D_cState &U1, const Euler1D_cState &U2) {
   return (U1.d*U2.d+U1.dv*U2.dv+U1.E*U2.E);
}

inline Euler1D_cState operator *(const Euler1D_cState &U, const double &a) {
  return (Euler1D_cState(a*U.d,a*U.dv,a*U.E));
}

inline Euler1D_cState operator *(const double &a, const Euler1D_cState &U) {
  return (Euler1D_cState(a*U.d,a*U.dv,a*U.E));
}

inline Euler1D_cState operator /(const Euler1D_cState &U, const double &a) {
  return (Euler1D_cState(U.d/a,U.dv/a,U.E/a));
}

// My useful solution state product operator.
inline Euler1D_cState operator ^(const Euler1D_cState &U1, const Euler1D_cState &U2) {
   return (Euler1D_cState(U1.d*U2.d,U1.dv*U2.dv,U1.E*U2.E));
}

/********************************************************
 * Euler1D_cState -- Unary arithmetic operators.        *
 ********************************************************/
inline Euler1D_cState operator +(const Euler1D_cState &U) {
  return (Euler1D_cState(U.d,U.dv,U.E));
}

inline Euler1D_cState operator -(const Euler1D_cState &U) {
  return (Euler1D_cState(-U.d,-U.dv,-U.E));
}

/********************************************************
 * Euler1D_cState -- Shortcut arithmetic operators.     *
 ********************************************************/
inline Euler1D_cState &operator +=(Euler1D_cState &U1, const Euler1D_cState &U2) {
  U1.d += U2.d;
  U1.dv += U2.dv;
  U1.E += U2.E;
  return (U1);
}

inline Euler1D_cState &operator -=(Euler1D_cState &U1, const Euler1D_cState &U2) {
  U1.d -= U2.d;
  U1.dv -= U2.dv;
  U1.E -= U2.E;
  return (U1);
}

/********************************************************
 * Euler1D_cState -- Relational operators.              *
 ********************************************************/
inline int operator ==(const Euler1D_cState &U1, const Euler1D_cState &U2) {
  return (U1.d == U2.d && U1.dv == U2.dv && U1.E == U2.E);
}

inline int operator !=(const Euler1D_cState &U1, const Euler1D_cState &U2) {
  return (U1.d != U2.d || U1.dv != U2.dv || U1.E != U2.E);
}

/********************************************************
 * Euler1D_cState -- Input-output operators.            *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const Euler1D_cState &U) {
  out_file.setf(ios::scientific);
  out_file << " " << U.d << " " << U.dv << " " << U.E;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, Euler1D_cState &U) {
  in_file.setf(ios::skipws);
  in_file >> U.d >> U.dv >> U.E;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/********************************************************
 * Euler1D_pState::Euler1D_pState -- Constructor.       *
 ********************************************************/
inline Euler1D_pState::Euler1D_pState(const Euler1D_cState &U) {
  d = U.d; v = U.v(); p = U.p();
}

/********************************************************
 * Euler1D_pState::U -- Conserved solution state.       *
 ********************************************************/
inline Euler1D_cState Euler1D_pState::U(void) {
  return (Euler1D_cState(d, dv(), E()));
}

inline Euler1D_cState Euler1D_pState::U(void) const {
  return (Euler1D_cState(d, dv(), E()));
}

inline Euler1D_cState Euler1D_pState::U(const Euler1D_pState &W) {
  return (Euler1D_cState(W.d, W.dv(), W.E()));
}

inline Euler1D_cState U(const Euler1D_pState &W) {
  return (Euler1D_cState(W.d, W.dv(), W.E()));
}

/********************************************************
 * Euler1D_pState::F -- Solution flux.                  *
 ********************************************************/
inline Euler1D_cState Euler1D_pState::F(void) {
  return (Euler1D_cState(dv(), d*sqr(v) + p, v*H()));
}

inline Euler1D_cState Euler1D_pState::F(void) const {
  return (Euler1D_cState(dv(), d*sqr(v) + p, v*H()));
}

inline Euler1D_cState Euler1D_pState::F(const Euler1D_pState &W) {
  return (Euler1D_cState(W.dv(), W.d*sqr(W.v) + W.p, W.v*W.H()));
}

inline Euler1D_cState F(const Euler1D_pState &W) {
  return (Euler1D_cState(W.dv(), W.d*sqr(W.v) + W.p, W.v*W.H()));
}

/********************************************************
 * Euler1D_pState::C -- Characteristic variables.       *
 ********************************************************/
inline Euler1D_pState Euler1D_pState::C(void) {
  double c = a();
  return (Euler1D_pState(TWO*gm1i*c-v, p/pow(d, g),
                         TWO*gm1i*c+v));
}

inline Euler1D_pState Euler1D_pState::C(void) const {
  double c = a();
  return (Euler1D_pState(TWO*gm1i*c-v, p/pow(d, g),
                         TWO*gm1i*c+v));
}

inline Euler1D_pState Euler1D_pState::C(const Euler1D_pState &W) {
  double c = W.a();
  return (Euler1D_pState(TWO*W.gm1i*c-W.v, W.p/pow(W.d, W.g),
                         TWO*W.gm1i*c+W.v));
}

inline Euler1D_pState C(const Euler1D_pState &W) {
  double c = W.a();
  return (Euler1D_pState(TWO*W.gm1i*c-W.v, W.p/pow(W.d, W.g),
                         TWO*W.gm1i*c+W.v));
}

inline Euler1D_pState CtoW(const Euler1D_pState &W) {
  double c = (W[3]+W[1])*W.gm1/FOUR; double u = (W[3]-W[1])/TWO;
  double dd = pow(sqr(c)/(W.g*W[2]), W.gm1i);
  return (Euler1D_pState(dd, u, sqr(c)*dd/W.g));
}

/* Transformation of the Characteristic Variables To Conserved Variables */
inline Euler1D_cState Euler1D_pState::CharactVarToConservedVar(const double SoundSpeed, const double Vel){
  /**** The characteristic variables ****
	Wc[1] = d;
	Wc[2] = v;
	Wc[3] = p; 
	Euler1D_cState( Wc[1]+Wc[2]+Wc[3], (v-c)*Wc[1]+v*Wc[2]+(v+c)*Wc[3] , (H-c*v)*Wc[1]+0.5*v^2*Wc[2]+(H+c*v)*Wc[3])
  ***************************************/
  double H = sqr(SoundSpeed)*gm1i + 0.5*sqr(Vel);
  /*   return Euler1D_cState( d + v + p, */
  /* 			    (Vel-SoundSpeed)*d + Vel*v + (Vel+SoundSpeed)*p, */
  /* 			    (H-SoundSpeed*Vel)*d + 0.5*sqr(Vel)*v + (H+SoundSpeed*Vel)*p ); */
  // The form below is only more compact than the form above
  return Euler1D_cState( d + v + p,
			 Vel*(d+v+p) - SoundSpeed*(d-p),
			 H*(d+p) - SoundSpeed*Vel*(d-p) + 0.5*sqr(Vel)*v );
}

inline Euler1D_cState Euler1D_pState::CharactVarToConservedVar(const Euler1D_pState &W_Ref){
  return CharactVarToConservedVar(W_Ref.a(), W_Ref.v);
}

inline Euler1D_cState Euler1D_pState::CharactVarToConservedVar(const Euler1D_cState &U_Ref){
  /* U -> the reference state that is used for computing the right eigenvectors
     "Wc" is a state storing primitive variables!!! */

  return CharactVarToConservedVar(U_Ref.a(), U_Ref.v());
  /*   double c = U_Ref.a(); */
  /*   double vel = U_Ref.v(); */
  /*   double H = U_Ref.h(); */
  
  /*   return Euler1D_cState( d+v+p, (vel-c)*d+vel*v+(vel+c)*p, (H-c*vel)*d+0.5*vel*vel*v+(H+c*vel)*p); */
}

inline Euler1D_cState CharactVarToConservedVar(const Euler1D_pState &Wc, const Euler1D_cState &U_Ref){
  /* U -> the reference state that is used for computing the right eigenvectors
     Wc -> the characteristic variables that are used for the transformation
     "Wc" is a state storing primitive variables!!! */

  double SoundSpeed = U_Ref.a();
  double Vel = U_Ref.v();
  double H = U_Ref.h();
  
  return Euler1D_cState( Wc.d+Wc.v+Wc.p,
			 Vel*(Wc.d+Wc.v+Wc.p) - SoundSpeed*(Wc.d-Wc.p),
			 H*(Wc.d+Wc.p) - SoundSpeed*Vel*(Wc.d-Wc.p) + 0.5*sqr(Vel)*Wc.v );
}

inline Euler1D_pState Euler1D_cState::ConservedVarToCharactVar(const double SoundSpeed, const double Vel){
  double omega = gm1*(0.5*sqr(Vel)*d - Vel*dv + E)/(SoundSpeed*SoundSpeed);
  return Euler1D_pState( 0.5*(omega + (Vel*d - dv)/SoundSpeed), d-omega, 0.5*(omega - (Vel*d - dv)/SoundSpeed) );
}

inline Euler1D_pState Euler1D_cState::ConservedVarToCharactVar(const double SoundSpeed, const double Vel) const{
  double omega = gm1*(0.5*sqr(Vel)*d - Vel*dv + E)/(SoundSpeed*SoundSpeed);
  return Euler1D_pState( 0.5*(omega + (Vel*d - dv)/SoundSpeed), d-omega, 0.5*(omega - (Vel*d - dv)/SoundSpeed) );
}

inline Euler1D_pState Euler1D_cState::ConservedVarToCharactVar(void){
  return ConservedVarToCharactVar(a(),v());
}

inline Euler1D_pState Euler1D_cState::ConservedVarToCharactVar(void) const{
  return ConservedVarToCharactVar(a(),v());
}

inline Euler1D_pState Euler1D_cState::ConservedVarToCharactVar(const Euler1D_cState &U){
  return ConservedVarToCharactVar(U.a(),U.v());
}

inline Euler1D_pState Euler1D_cState::ConservedVarToCharactVar(const Euler1D_cState &U) const{
  return ConservedVarToCharactVar(U.a(),U.v());
}

inline Euler1D_pState Euler1D_cState::ConservedVarToCharactVar(const Euler1D_pState &W){
  return ConservedVarToCharactVar(W.a(),W.v);
}

inline Euler1D_pState Euler1D_cState::ConservedVarToCharactVar(const Euler1D_pState &W) const{
  return ConservedVarToCharactVar(W.a(),W.v);
}

/********************************************************
 * Euler1D_pState::lambda -- Eigenvalue(s).             *
 ********************************************************/
inline Euler1D_pState Euler1D_pState::lambda(void) {
  double c = a();
  return (Euler1D_pState(v - c, v, v + c));
}

inline Euler1D_pState Euler1D_pState::lambda(void) const {
  double c = a();
  return (Euler1D_pState(v - c, v, v + c));
}

inline Euler1D_pState Euler1D_pState::lambda(const Euler1D_pState &W) {
  double c = W.a();
  return (Euler1D_pState(W.v - c, W.v, W.v + c));
}

inline Euler1D_pState lambda(const Euler1D_pState &W) {
  double c = W.a();
  return (Euler1D_pState(W.v - c, W.v, W.v + c));
}

inline double Euler1D_pState::lambda(int index) {
  assert( index >= 1 && index <= NUM_VAR_EULER1D );
  switch(index) {
    case 1 :
      return (v-a());
    case 2 :
      return (v);
    case 3 :
      return (v+a());
    default:
      return (v);
  };
}

inline double Euler1D_pState::lambda(int index) const {
  assert( index >= 1 && index <= NUM_VAR_EULER1D );
  switch(index) {
    case 1 :
      return (v-a());
    case 2 :
      return (v);
    case 3 :
      return (v+a());
    default:
      return (v);
  };
}

inline double lambda(const Euler1D_pState &W, int index) {
  assert( index >= 1 && index <= NUM_VAR_EULER1D );
  switch(index) {
    case 1 :
      return (W.v-W.a());
    case 2 :
      return (W.v);
    case 3 :
      return (W.v+W.a());
    default:
      return (W.v);
  };
}

/********************************************************
 * Euler1D_pState::rp -- Primitive right eigenvector.   *
 ********************************************************/
inline Euler1D_pState Euler1D_pState::rp(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER1D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler1D_pState(ONE, -c/d, sqr(c)));
    case 2 :
      return (Euler1D_pState(ONE, ZERO, ZERO));
    case 3 :
      c = a();
      return (Euler1D_pState(ONE, c/d, sqr(c)));
    default:
      c = a();
      return (Euler1D_pState(ONE, -c/d, sqr(c)));
  };
}

inline Euler1D_pState Euler1D_pState::rp(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER1D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler1D_pState(ONE, -c/d, sqr(c)));
    case 2 :
      return (Euler1D_pState(ONE, ZERO, ZERO));
    case 3 :
      c = a();
      return (Euler1D_pState(ONE, c/d, sqr(c)));
    default:
      c = a();
      return (Euler1D_pState(ONE, -c/d, sqr(c)));
  };
}

inline Euler1D_pState rp(const Euler1D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER1D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Euler1D_pState(ONE, -c/W.d, sqr(c)));
    case 2 :
      return (Euler1D_pState(ONE, ZERO, ZERO));
    case 3 :
      c = W.a();
      return (Euler1D_pState(ONE, c/W.d, sqr(c)));
    default:
      c = W.a();
      return (Euler1D_pState(ONE, -c/W.d, sqr(c)));
  };
}

/********************************************************
 * Euler1D_pState::rc -- Conserved right eigenvector.   *
 ********************************************************/
inline Euler1D_cState Euler1D_pState::rc(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER1D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler1D_cState(ONE, v-c, h()-v*c));
    case 2 :
      return (Euler1D_cState(ONE, v, HALF*sqr(v)));
    case 3 :
      c = a();
      return (Euler1D_cState(ONE, v+c, h()+v*c));
    default:
      c = a();
      return (Euler1D_cState(ONE, v-c, h()-v*c));
  };
}

inline Euler1D_cState Euler1D_pState::rc(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER1D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler1D_cState(ONE, v-c, h()-v*c));
    case 2 :
      return (Euler1D_cState(ONE, v, HALF*sqr(v)));
    case 3 :
      c = a();
      return (Euler1D_cState(ONE, v+c, h()+v*c));
    default:
      c = a();
      return (Euler1D_cState(ONE, v-c, h()-v*c));
  };
}

inline Euler1D_cState rc(const Euler1D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER1D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Euler1D_cState(ONE, W.v-c, W.h()-W.v*c));
    case 2 :
      return (Euler1D_cState(ONE, W.v, HALF*sqr(W.v)));
    case 3 :
      c = W.a();
      return (Euler1D_cState(ONE, W.v+c, W.h()+W.v*c));
    default:
      c = W.a();
      return (Euler1D_cState(ONE, W.v-c, W.h()-W.v*c));
  };
}

/********************************************************
 * Euler1D_pState::lp -- Primitive left eigenvector.    *
 ********************************************************/
inline Euler1D_pState Euler1D_pState::lp(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER1D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler1D_pState(ZERO, -HALF*d/c, HALF/sqr(c)));
    case 2 :
      return (Euler1D_pState(ONE, ZERO, -ONE/a2()));
    case 3 :
      c = a();
      return (Euler1D_pState(ZERO, HALF*d/c, HALF/sqr(c)));
    default:
      c = a();
      return (Euler1D_pState(ZERO, -HALF*d/c, HALF/sqr(c)));
  };
}

inline Euler1D_pState Euler1D_pState::lp(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER1D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler1D_pState(ZERO, -HALF*d/c, HALF/sqr(c)));
    case 2 :
      return (Euler1D_pState(ONE, ZERO, -ONE/a2()));
    case 3 :
      c = a();
      return (Euler1D_pState(ZERO, HALF*d/c, HALF/sqr(c)));
    default:
      c = a();
      return (Euler1D_pState(ZERO, -HALF*d/c, HALF/sqr(c)));
  };
}

inline Euler1D_pState lp(const Euler1D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER1D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Euler1D_pState(ZERO, -HALF*W.d/c, HALF/sqr(c)));
    case 2 :
      return (Euler1D_pState(ONE, ZERO, -ONE/W.a2()));
    case 3 :
      c = W.a();
      return (Euler1D_pState(ZERO, HALF*W.d/c, HALF/sqr(c)));
    default:
      c = W.a();
      return (Euler1D_pState(ZERO, -HALF*W.d/c, HALF/sqr(c)));
  };
}

/********************************************************
 * Euler1D_cState::Euler1D_cState -- Constructor.       *
 ********************************************************/
inline Euler1D_cState::Euler1D_cState(const Euler1D_pState &W) {
  d = W.d; dv = W.dv(); E = W.E();
}

/********************************************************
 * Euler1D_cState::W -- Primitive solution state.       *
 ********************************************************/
inline Euler1D_pState Euler1D_cState::W(void) {
  return (Euler1D_pState(d, v(), p()));
}

inline Euler1D_pState Euler1D_cState::W(void) const {
  return (Euler1D_pState(d, v(), p()));
}

inline Euler1D_pState Euler1D_cState::W(const Euler1D_cState &U) {
  return (Euler1D_pState(U.d, U.v(), U.p()));
}

inline Euler1D_pState W(const Euler1D_cState &U) {
  return (Euler1D_pState(U.d, U.v(), U.p()));
}

/********************************************************
 * Euler1D_cState::F -- Solution flux.                  *
 ********************************************************/
inline Euler1D_cState Euler1D_cState::F(void) {
  return (Euler1D_cState(dv, sqr(dv)/d + p(), dv*H()/d));
}

inline Euler1D_cState Euler1D_cState::F(void) const {
  return (Euler1D_cState(dv, sqr(dv)/d + p(), dv*H()/d));
}

inline Euler1D_cState Euler1D_cState::F(const Euler1D_cState &U) {
  return (Euler1D_cState(U.dv, sqr(U.dv)/U.d + U.p(), U.dv*U.H()/U.d));
}

inline Euler1D_cState F(const Euler1D_cState &U) {
  return (Euler1D_cState(U.dv, sqr(U.dv)/U.d + U.p(), U.dv*U.H()/U.d));
}

/********************************************************
 * Euler1D_cState::C -- Characteristics variables.      *
 ********************************************************/
inline Euler1D_pState Euler1D_cState::C(void) {
  double c = a();
  return (Euler1D_pState(TWO*gm1i*c-dv/d, p()/pow(d, g),
			 TWO*gm1i*c+dv/d));
}

inline Euler1D_pState Euler1D_cState::C(void) const {
  double c = a();
  return (Euler1D_pState(TWO*gm1i*c-dv/d, p()/pow(d, g),
			 TWO*gm1i*c+dv/d));
}

inline Euler1D_pState Euler1D_cState::C(const Euler1D_cState &U) {
  double c = U.a();
  return (Euler1D_pState(TWO*U.gm1i*c-U.dv/U.d, U.p()/pow(U.d, U.g),
			 TWO*U.gm1i*c+U.dv/U.d));
}

inline Euler1D_pState C(const Euler1D_cState &U) {
  double c = U.a();
  return (Euler1D_pState(TWO*U.gm1i*c-U.dv/U.d, U.p()/pow(U.d, U.g),
			 TWO*U.gm1i*c+U.dv/U.d));
}

inline Euler1D_cState CtoU(const Euler1D_pState &W) {
  double c = (W[3]+W[1])*W.gm1/FOUR; double u = (W[3]-W[1])/TWO;
  double dd = pow(sqr(c)/(W.g*W[2]), W.gm1i);
  return (Euler1D_cState(dd, dd*u, dd*sqr(u)+W.gm1i*sqr(c)*dd/W.g));
}

/********************************************************
 * Useful 1D Euler state constants.                     *
 ********************************************************/
const Euler1D_pState Euler1D_W_STDATM(DENSITY_STDATM, ZERO, PRESSURE_STDATM);
const Euler1D_pState Euler1D_W_VACUUM(ZERO, ZERO, ZERO);
const Euler1D_pState Euler1D_W_ONE(ONE, ONE, ONE);
const Euler1D_cState Euler1D_U_STDATM(Euler1D_W_STDATM);
const Euler1D_cState Euler1D_U_VACUUM(Euler1D_W_VACUUM);
const Euler1D_cState Euler1D_U_ONE(ONE, ONE, ONE);

/********************************************************
 * Euler1DState -- External subroutines.                *
 ********************************************************/

extern Euler1D_pState Riemann(const Euler1D_pState &Wl,
	      	              const Euler1D_pState &Wr);

extern Euler1D_pState RoeAverage(const Euler1D_pState &Wl,
	      	                 const Euler1D_pState &Wr);

extern Euler1D_pState WaveSpeedPos(const Euler1D_pState &lambda_a,
                                   const Euler1D_pState &lambda_l,
                                   const Euler1D_pState &lambda_r);

extern Euler1D_pState WaveSpeedNeg(const Euler1D_pState &lambda_a,
                                   const Euler1D_pState &lambda_l,
                                   const Euler1D_pState &lambda_r);

extern Euler1D_pState WaveSpeedAbs(const Euler1D_pState &lambda_a,
                                   const Euler1D_pState &lambda_l,
                                   const Euler1D_pState &lambda_r);

extern Euler1D_pState HartenFixPos(const Euler1D_pState &lambda_a,
                                   const Euler1D_pState &lambda_l,
                                   const Euler1D_pState &lambda_r);

extern Euler1D_pState HartenFixNeg(const Euler1D_pState &lambda_a,
                                   const Euler1D_pState &lambda_l,
                                   const Euler1D_pState &lambda_r);

extern Euler1D_pState HartenFixAbs(const Euler1D_pState &lambda_a,
                                   const Euler1D_pState &lambda_l,
                                   const Euler1D_pState &lambda_r);

extern Euler1D_cState FluxGodunov(const Euler1D_pState &Wl,
	      	                  const Euler1D_pState &Wr);

extern Euler1D_cState FluxGodunov(const Euler1D_cState &Ul,
	      	                  const Euler1D_cState &Ur);

extern Euler1D_cState FluxRoe(const Euler1D_pState &Wl,
	      	              const Euler1D_pState &Wr);

extern Euler1D_cState FluxRoe(const Euler1D_cState &Ul,
	      	              const Euler1D_cState &Ur);

extern Euler1D_cState FluxRusanov(const Euler1D_pState &Wl,
	      	                  const Euler1D_pState &Wr);

extern Euler1D_cState FluxRusanov(const Euler1D_cState &Ul,
	      	                  const Euler1D_cState &Ur);

extern Euler1D_cState FluxHLLE(const Euler1D_pState &Wl,
	      	               const Euler1D_pState &Wr);

extern Euler1D_cState FluxHLLE(const Euler1D_cState &Ul,
	      	               const Euler1D_cState &Ur);
  
extern Euler1D_cState FluxLinde(const Euler1D_pState &Wl,
	      	                const Euler1D_pState &Wr);

extern Euler1D_cState FluxLinde(const Euler1D_cState &Ul,
	      	                const Euler1D_cState &Ur);

extern Euler1D_cState FluxHLLC(const Euler1D_pState &Wl,
	      	               const Euler1D_pState &Wr);

extern Euler1D_cState FluxHLLC(const Euler1D_cState &Ul,
	      	               const Euler1D_cState &Ur);

extern Euler1D_cState FluxOsher(const Euler1D_pState &Wl,
	      	                const Euler1D_pState &Wr);

extern Euler1D_cState FluxOsher(const Euler1D_cState &Ul,
	      	                const Euler1D_cState &Ur);

extern Euler1D_cState RiemannFlux(const int & Flux_Function,
				  const Euler1D_pState &Wl,
				  const Euler1D_pState &Wr);

extern Euler1D_pState BC_Reflection(const Euler1D_pState &Wi,
                                    const double Vbnd,
                                    const int End_Type);

extern Euler1D_pState BC_Characteristic(const Euler1D_pState &Wi,
	      	                        const Euler1D_pState &Wo,
                                        const int End_Type);

extern Euler1D_pState BC_Characteristic_Pressure(const Euler1D_pState &Wi,
	      	                                 const Euler1D_pState &Wo,
                                                 const int End_Type);

extern Euler1D_pState BC_Characteristic_Mach_Number(const Euler1D_pState &Wi,
	      	                                    const Euler1D_pState &Wo,
                                                    const int End_Type);

extern Euler1D_pState BC_Open_End(const Euler1D_pState &Wi,
	      	                  const Euler1D_pState &Wo,
                                  const int End_Type);

#endif /* _EULER1D_STATE_INCLUDED  */
