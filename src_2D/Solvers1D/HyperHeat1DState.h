/* HyperHeat1DState.h:  Header file defining 1D Hyperbolic Heat Equations
                        (Maxwell-Cattaneo Equation) Solution State Class. */

#ifndef _HYPERHEAT1D_STATE_INCLUDED
#define _HYPERHEAT1D_STATE_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

/* Include math macro and CFD header files. */

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

/* Define the classes. */

#define	NUM_VAR_HYPERHEAT1D    2

/********************************************************
 * Class: HyperHeat1D_State                             *
 *                                                      *
 * Member functions                                     *
 *      T       -- Return temperature.                  *
 *      qx      -- Return heat flux.                    *
 *      kappa   -- Return thermal conductivity.         *
 *      tau     -- Return relaxation time.              *
 *      setcon  -- Set conductivity & relaxation time.  *
 *      a       -- Return wave speed.                   *
 *      a2      -- Return wave speed square.            *
 *      F       -- Return solution flux.                *
 *      S       -- Return solution sources.             *
 *      lambda  -- Return eigenvalue(s).                *
 *      r       -- Return right eigenvectors.           *
 *      l       -- Return left eigenvectors.            *
 *                                                      *
 * Member operators                                     *
 *      U -- a solution state                           *
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
class HyperHeat1D_State{
  private:
  public:
    double            T;   // Temperature.
    double           qx;   // Heat flux.
    static double kappa;   // Thermal conductivity.
    static double   tau;   // Relaxation time.
	                   // Made public so can access them.
		      
    /* Creation, copy, and assignment constructors. */
    HyperHeat1D_State(void) {
       T = ZERO; qx = ZERO;
    }

    HyperHeat1D_State(const HyperHeat1D_State &U) {
       T = U.T; qx = U.qx;
    }

    HyperHeat1D_State(const double &Temp,
	              const double &hx) {
       T = Temp; qx = hx;
    }
    
    /* Destructor. */
    // ~HyperHeat1D_State(void);
    // Use automatically generated destructor.

    /* Set thermal conductivity and relaxatime time. */
    void setcon(void);
    void setcon(const double &kk,
	        const double &tt);

    /* Wave speed. */
    double a(void);
    double a(void) const;

    /* Wave  speed squared. */
    double a2(void);
    double a2(void) const;

    /* Solution flux. */
    HyperHeat1D_State F(void);
    HyperHeat1D_State F(void) const;
    HyperHeat1D_State F(const HyperHeat1D_State &U);
    friend HyperHeat1D_State F(const HyperHeat1D_State &U);

    /* Solution source. */
    HyperHeat1D_State S(void);
    HyperHeat1D_State S(void) const;
    HyperHeat1D_State S(const HyperHeat1D_State &U);
    friend HyperHeat1D_State S(const HyperHeat1D_State &U);

    /* Eigenvalue(s). */
    HyperHeat1D_State lambda(void);
    HyperHeat1D_State lambda(void) const;
    HyperHeat1D_State lambda(const HyperHeat1D_State &U);
    friend HyperHeat1D_State lambda(const HyperHeat1D_State &U);
    double lambda(int index);
    double lambda(int index) const;
    friend double lambda(const HyperHeat1D_State &U, int index);

    /* Right eigenvector. */
    HyperHeat1D_State r(int index);
    HyperHeat1D_State r(int index) const;
    friend HyperHeat1D_State r(const HyperHeat1D_State &U, int index);

    /* Left eigenvector. */
    HyperHeat1D_State l(int index);
    HyperHeat1D_State l(int index) const;
    friend HyperHeat1D_State l(const HyperHeat1D_State &U, int index);
    
    /* Assignment operator. */
    // HyperHeat1D_State operator = (const HyperHeat1D_State &W);
    // Use automatically generated assignment operator.

    /* Index operator. */
    double &operator[](int index) {
      assert( index >= 1 && index <= NUM_VAR_HYPERHEAT1D );
      switch(index) {
        case 1 :
	  return (T);
        case 2 :
	  return (qx);
        default:
	  return (T);
      };
    }
    
    const double &operator[](int index) const {
      assert( index >= 1 && index <= NUM_VAR_HYPERHEAT1D );
      switch(index) {
        case 1 :
	  return (T);
        case 2 :
	  return (qx);
        default:
	  return (T);
      };
    }

    /* Binary arithmetic operators. */
    friend HyperHeat1D_State operator +(const HyperHeat1D_State &U1, const HyperHeat1D_State &U2);
    friend HyperHeat1D_State operator -(const HyperHeat1D_State &U1, const HyperHeat1D_State &U2);
    friend double operator *(const HyperHeat1D_State &U1, const HyperHeat1D_State &U2);
    friend HyperHeat1D_State operator *(const HyperHeat1D_State &U, const double &a);
    friend HyperHeat1D_State operator *(const double &a, const HyperHeat1D_State &U);
    friend HyperHeat1D_State operator /(const HyperHeat1D_State &U, const double &a);
    friend HyperHeat1D_State operator ^(const HyperHeat1D_State &U1,
					const HyperHeat1D_State &U2);

    /* Unary arithmetic operators. */
    friend HyperHeat1D_State operator +(const HyperHeat1D_State &U);
    friend HyperHeat1D_State operator -(const HyperHeat1D_State &U);

    /* Shortcut arithmetic operators. */
    friend HyperHeat1D_State &operator +=(HyperHeat1D_State &U1, const HyperHeat1D_State &U2);
    friend HyperHeat1D_State &operator -=(HyperHeat1D_State &U1, const HyperHeat1D_State &U2);
    
    /* Relational operators. */
    friend int operator ==(const HyperHeat1D_State &U1, const HyperHeat1D_State &U2);
    friend int operator !=(const HyperHeat1D_State &U1, const HyperHeat1D_State &U2);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const HyperHeat1D_State &U);
    friend istream &operator >> (istream &in_file,  HyperHeat1D_State &U);
    
};

/********************************************************
 * HyperHeat1D_State::setcon -- Assign thermal          *
 *                              conductivity and        *
 *                              relaxation time.        *
 ********************************************************/
inline void HyperHeat1D_State::setcon(void) {
    kappa = TWO/HUNDRED; tau = ONE/HUNDRED;
}

inline void HyperHeat1D_State::setcon(const double &kk,
	                              const double &tt) {
    kappa = kk; tau = tt;
}

/********************************************************
 * HyperHeat1D_State::a -- Wave speed.                  *
 ********************************************************/
inline double HyperHeat1D_State::a(void) {
    return (sqrt(kappa/tau));
}

inline double HyperHeat1D_State::a(void) const {
    return (sqrt(kappa/tau));
}

/********************************************************
 * HyperHeat1D_State::a2 -- Wave speed squared.         *
 ********************************************************/
inline double HyperHeat1D_State::a2(void) {
    return (kappa/tau);
}

inline double HyperHeat1D_State::a2(void) const {
    return (kappa/tau);
}

/********************************************************
 * HyperHeat1D_State -- Binary arithmetic operators.    *
 ********************************************************/
inline HyperHeat1D_State operator +(const HyperHeat1D_State &U1, const HyperHeat1D_State &U2) {
  return (HyperHeat1D_State(U1.T+U2.T,U1.qx+U2.qx));
}

inline HyperHeat1D_State operator -(const HyperHeat1D_State &U1, const HyperHeat1D_State &U2) {
  return (HyperHeat1D_State(U1.T-U2.T,U1.qx-U2.qx));
}

// Inner product operator.
inline double operator *(const HyperHeat1D_State &U1, const HyperHeat1D_State &U2) {
   return (U1.T*U2.T+U1.qx*U2.qx);
}

inline HyperHeat1D_State operator *(const HyperHeat1D_State &U, const double &a) {
  return (HyperHeat1D_State(a*U.T,a*U.qx));
}

inline HyperHeat1D_State operator *(const double &a, const HyperHeat1D_State &U) {
  return (HyperHeat1D_State(a*U.T,a*U.qx));
}

inline HyperHeat1D_State operator /(const HyperHeat1D_State &U, const double &a) {
  return (HyperHeat1D_State(U.T/a,U.qx/a));
}

// My useful solution state product operator.
inline HyperHeat1D_State operator ^(const HyperHeat1D_State &U1,
				    const HyperHeat1D_State &U2) {
   return (HyperHeat1D_State(U1.T*U2.T,U1.qx*U2.qx));
}

/********************************************************
 * HyperHeat1D_State -- Unary arithmetic operators.     *
 ********************************************************/
inline HyperHeat1D_State operator +(const HyperHeat1D_State &U) {
  return (HyperHeat1D_State(U.T,U.qx));
}

inline HyperHeat1D_State operator -(const HyperHeat1D_State &U) {
  return (HyperHeat1D_State(-U.T,-U.qx));
}

/********************************************************
 * HyperHeat1D_State -- Shortcut arithmetic operators.  *
 ********************************************************/
inline HyperHeat1D_State &operator +=(HyperHeat1D_State &U1, const HyperHeat1D_State &U2) {
  U1.T  += U2.T; U1.qx += U2.qx;
  return (U1);
}

inline HyperHeat1D_State &operator -=(HyperHeat1D_State &U1, const HyperHeat1D_State &U2) {
  U1.T  -= U2.T; U1.qx -= U2.qx;
  return (U1);
}

/********************************************************
 * HyperHeat1D_State -- Relational operators.           *
 ********************************************************/
inline int operator ==(const HyperHeat1D_State &U1, const HyperHeat1D_State &U2) {
  return (U1.T == U2.T && U1.qx == U2.qx);
}

inline int operator !=(const HyperHeat1D_State &U1, const HyperHeat1D_State &U2) {
  return (U1.T != U2.T || U1.qx != U2.qx);
}

/********************************************************
 * HyperHeat1D_State -- Input-output operators.         *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const HyperHeat1D_State &U) {
  out_file.setf(ios::scientific);
  out_file << " " << U.T << " " << U.qx;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, HyperHeat1D_State &U) {
  in_file.setf(ios::skipws);
  in_file >> U.T >> U.qx;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/********************************************************
 * HyperHeat1D_State::F -- Solution flux.               *
 ********************************************************/
inline HyperHeat1D_State HyperHeat1D_State::F(void) {
  return (HyperHeat1D_State(qx, kappa*T/tau));
}

inline HyperHeat1D_State HyperHeat1D_State::F(void) const {
  return (HyperHeat1D_State(qx, kappa*T/tau));
}

inline HyperHeat1D_State HyperHeat1D_State::F(const HyperHeat1D_State &U) {
  return (HyperHeat1D_State(U.qx, U.kappa*U.T/U.tau));
}

inline HyperHeat1D_State F(const HyperHeat1D_State &U) {
  return (HyperHeat1D_State(U.qx, U.kappa*U.T/U.tau));
}

/********************************************************
 * HyperHeat1D_State::S -- Solution source.               *
 ********************************************************/
inline HyperHeat1D_State HyperHeat1D_State::S(void) {
  return (HyperHeat1D_State(ZERO, -qx/tau));
}

inline HyperHeat1D_State HyperHeat1D_State::S(void) const {
  return (HyperHeat1D_State(ZERO, -qx/tau));
}

inline HyperHeat1D_State HyperHeat1D_State::S(const HyperHeat1D_State &U) {
  return (HyperHeat1D_State(ZERO, -U.qx/U.tau));
}

inline HyperHeat1D_State S(const HyperHeat1D_State &U) {
  return (HyperHeat1D_State(ZERO, -U.qx/U.tau));
}

/********************************************************
 * HyperHeat1D_State::lambda -- Eigenvalue(s).          *
 ********************************************************/
inline HyperHeat1D_State HyperHeat1D_State::lambda(void) {
  double c = a();
  return (HyperHeat1D_State(-c, c));
}

inline HyperHeat1D_State HyperHeat1D_State::lambda(void) const {
  double c = a();
  return (HyperHeat1D_State(-c, c));
}

inline HyperHeat1D_State HyperHeat1D_State::lambda(const HyperHeat1D_State &U) {
  double c = U.a();
  return (HyperHeat1D_State(-c, c));
}

inline HyperHeat1D_State lambda(const HyperHeat1D_State &U) {
  double c = U.a();
  return (HyperHeat1D_State(-c, c));
}

inline double HyperHeat1D_State::lambda(int index) {
  assert( index >= 1 && index <= NUM_VAR_HYPERHEAT1D );
  switch(index) {
    case 1 :
      return (-a());
    case 2 :
      return (a());
    default:
      return (a());
  };
}

inline double HyperHeat1D_State::lambda(int index) const {
  assert( index >= 1 && index <= NUM_VAR_HYPERHEAT1D );
  switch(index) {
    case 1 :
      return (-a());
    case 2 :
      return (a());
    default:
      return (a());
  };
}

inline double lambda(const HyperHeat1D_State &U, int index) {
  assert( index >= 1 && index <= NUM_VAR_HYPERHEAT1D );
  switch(index) {
    case 1 :
      return (-U.a());
    case 2 :
      return (U.a());
    default:
      return (U.a());
  };
}

/********************************************************
 * HyperHeat1D_State::r -- Right eigenvector.           *
 ********************************************************/
inline HyperHeat1D_State HyperHeat1D_State::r(int index) {
  assert( index >= 1 && index <= NUM_VAR_HYPERHEAT1D );
  switch(index) {
    case 1 :
      return (HyperHeat1D_State(ONE, -a()));
    case 2 :
      return (HyperHeat1D_State(ONE, a()));
    default:
      return (HyperHeat1D_State(ONE, a()));
  };
}

inline HyperHeat1D_State HyperHeat1D_State::r(int index) const {
  assert( index >= 1 && index <= NUM_VAR_HYPERHEAT1D );
  switch(index) {
    case 1 :
      return (HyperHeat1D_State(ONE, -a()));
    case 2 :
      return (HyperHeat1D_State(ONE, a()));
    default:
      return (HyperHeat1D_State(ONE, a()));
  };
}

inline HyperHeat1D_State r(const HyperHeat1D_State &U, int index) {
  assert( index >= 1 && index <= NUM_VAR_HYPERHEAT1D );
  switch(index) {
    case 1 :
      return (HyperHeat1D_State(ONE, -U.a()));
    case 2 :
      return (HyperHeat1D_State(ONE, U.a()));
    default:
      return (HyperHeat1D_State(ONE, U.a()));
  };
}

/********************************************************
 * HyperHeat1D_State::l -- Left eigenvector.            *
 ********************************************************/
inline HyperHeat1D_State HyperHeat1D_State::l(int index) {
  assert( index >= 1 && index <= NUM_VAR_HYPERHEAT1D );
  switch(index) {
    case 1 :
      return (HyperHeat1D_State(HALF, -HALF/a()));
    case 2 :
      return (HyperHeat1D_State(HALF, HALF/a()));
    default:
      return (HyperHeat1D_State(HALF, HALF/a()));
  };
}

inline HyperHeat1D_State HyperHeat1D_State::l(int index) const {
  assert( index >= 1 && index <= NUM_VAR_HYPERHEAT1D );
  switch(index) {
    case 1 :
      return (HyperHeat1D_State(HALF, -HALF/a()));
    case 2 :
      return (HyperHeat1D_State(HALF, HALF/a()));
    default:
      return (HyperHeat1D_State(HALF, HALF/a()));
  };
}

inline HyperHeat1D_State l(const HyperHeat1D_State &U, int index) {
  assert( index >= 1 && index <= NUM_VAR_HYPERHEAT1D );
  switch(index) {
    case 1 :
      return (HyperHeat1D_State(HALF, -HALF/U.a()));
    case 2 :
      return (HyperHeat1D_State(HALF, HALF/U.a()));
    default:
      return (HyperHeat1D_State(HALF, HALF/U.a()));
  };
}

/******************************************************** 
 * Useful 1D hyperbolic heat equations solution states. *
 ********************************************************/
const HyperHeat1D_State HyperHeat1D_U_ZERO(ZERO, ZERO);
const HyperHeat1D_State HyperHeat1D_U_ONE(ONE, ONE);

/********************************************************
 * HyperHeat1DState -- External subroutines.            *
 ********************************************************/

extern HyperHeat1D_State Riemann(const HyperHeat1D_State &Ul,
	      	                 const HyperHeat1D_State &Ur);

extern HyperHeat1D_State RiemannHomo(const HyperHeat1D_State &Ul,
	      	                     const HyperHeat1D_State &Ur);

extern HyperHeat1D_State BCs(const HyperHeat1D_State &Ub,
			     const HyperHeat1D_State &Ui,
			     const HyperHeat1D_State &dUdx,
			     const double &dx,
			     int BC_type,
			     int End_type);

extern HyperHeat1D_State WaveSpeedPos(const HyperHeat1D_State &lambdas);

extern HyperHeat1D_State WaveSpeedNeg(const HyperHeat1D_State &lambdas);

extern HyperHeat1D_State WaveSpeedAbs(const HyperHeat1D_State &lambdas);

extern HyperHeat1D_State FluxGodunov(const HyperHeat1D_State &Ul,
	      	                     const HyperHeat1D_State &Ur);

extern HyperHeat1D_State FluxRoe(const HyperHeat1D_State &Ul,
	      	                 const HyperHeat1D_State &Ur);

extern HyperHeat1D_State FluxRusanov(const HyperHeat1D_State &Ul,
	      	                     const HyperHeat1D_State &Ur);

extern HyperHeat1D_State FluxHLLE(const HyperHeat1D_State &Ul,
	      	                  const HyperHeat1D_State &Ur);

#endif /* _HYPERHEAT1D_STATE_INCLUDED  */
