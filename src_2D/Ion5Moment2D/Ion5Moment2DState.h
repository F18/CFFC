/* Ion5Moment2DState.h:  Header file defining 
                         2D 5-Moment Ion Transport Model
                         Solution State Classes. */

#ifndef _ION5MOMENT2D_STATE_INCLUDED
#define _ION5MOMENT2D_STATE_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

/* Include math macro, matrix, CFD, 2D vector, gas constant,
   and 2D Euler state header files. */

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _MATRIX_INCLUDED
#include "../Math/Matrix.h"
#endif // _MATRIX_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _VECTOR2D_INCLUDED
#include "../Math/Vector2D.h"
#endif //_VECTOR2D_INCLUDED

#ifndef _GAS_CONSTANTS_INCLUDED
#include "../Physics/GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED

#ifndef _EULER2D_STATE_INCLUDED
#include "../Euler2D/Euler2DState.h"
#endif // _EULER2D_STATE_INCLUDED

/* Define the classes. */

#define	NUM_VAR_ION5MOMENT2D    4

class Ion5Moment2D_cState;

/*!
 * Class: Ion5Moment2D_pState
 *
 * @brief Primitive variable solution state class definition for an
 *        inviscid Ion5moment flow.
 *
 * Primitive variable solution state class definition for an inviscid
 * Ion5moment flow.
 *
 * \verbatim
 * Member functions
 *     d        -- Return ion mass density.
 *     v        -- Return ion flow velocity.
 *     p        -- Return ion pressure.
 *     g        -- Return ion specific heat ratio.
 *     gm1      -- Return g-1.
 *     gm1i     -- Return 1/(g-1).
 *     R        -- Return ion gas constant.
 *     q        -- Return ion charge.
 *     m        -- Return ion particle mass.
 *     ion      -- Return ion type indicator.
 *     gas      -- Return neutral gas type indicator.
 *     setion   -- Set ion constants.
 *     setgas   -- Set neutral gas type.
 *     n        -- Return ion number density.
 *     T        -- Return ion temperature.
 *     e        -- Return specific internal energy.
 *     E        -- Return total energy.
 *     h        -- Return specific enthalpy.
 *     H        -- Return total enthalpy.
 *     a        -- Return sound speed.
 *     a2       -- Return sound speed square.
 *     M        -- Return Mach number.
 *     s        -- Return specific entropy.
 *     dv       -- Return momentum.
 *     To       -- Return stagnation temperature.
 *     po       -- Return stagnation pressure.
 *     ao       -- Return stagnation sound speed.
 *     ho       -- Return stagnation enthalpy.
 *     nu       -- Return ion-neutral collision frequency.
 *     mn       -- Return neutral gas particle mass.
 *     U        -- Return conserved solution state.
 *     F        -- Return x-direction solution flux.
 *     Fx       -- Return x-direction solution flux.
 *     Fy       -- Return y-direction solution flux.
 *     Fn       -- Return n-direction solution flux.
 *     lambda   -- Return x-direction eigenvalue(s).
 *     lambda_x -- Return x-direction eigenvalue(s).
 *     lambda_y -- Return y-direction eigenvalue(s).
 *     rp       -- Return primitive right eigenvector (x-direction).
 *     rp_x     -- Return primitive right eigenvector (x-direction).
 *     rp_y     -- Return primitive right eigenvector (y-direction).
 *     rc       -- Return conserved right eigenvector (x-direction).
 *     rc_x     -- Return conserved right eigenvector (x-direction).
 *     rc_y     -- Return conserved right eigenvector (y-direction).
 *     lp       -- Return primitive left eigenvector (x-direction).
 *     lp_x     -- Return primitive left eigenvector (x-direction).
 *     lp_y     -- Return primitive left eigenvector (y-direction).
 *     Sa       -- Return axisymmetric source term vector.
 *     Se       -- Return Coulomb force source vector produced by
 *                 external electric field.
 *     Sc       -- Return source terms associated with ion-neutral
 *                 collision processes
 *     dSadU    -- Return axisymmetric source term Jacobian.
 *     dSedU    -- Return Coulomb force source Jacobian.
 *     dScdU    -- Return source Jacobian for ion-neutral collision
 *                 collision processes.
 *
 * Member operators
 *      W -- a primitive solution state
 *      c -- a scalar (double)
 *
 * W = W;
 * c = W[i];
 * W = W + W;
 * W = W - W;
 * c = W * W; (inner product)
 * W = c * W;
 * W = W * c;
 * W = W / c;
 * W = W ^ W; (a useful product)
 * W = +W;
 * W = -W;
 * W += W;
 * W -= W;
 * W == W;
 * W != W;
 * cout << W; (output function)
 * cin  >> W; (input function)
 * \endverbatim
 */
class Ion5Moment2D_pState{
private:
public:
  //@{ @name Primitive variables and associated constants:
  double           d; //!< Ion mass density.
  Vector2D         v; //!< Ion flow velocity (2D vector).
  double           p; //!< Ion pressure.
  static double    g; //!< Ion specific heat ratio.
  static double  gm1; //!< g-1
  static double gm1i; //!< 1/(g-1)
  static double    R; //!< Ion gas constant.
  static double    q; //!< Ion charge.
  static double    m; //!< Ion particle mass.
  static int     ion; //!< Ion type indicator.
  static int     gas; //!< Neutral gas type indicator.
  //@}
		      
  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  Ion5Moment2D_pState(void) {
    d = DENSITY_STDATM; v.zero(); p = PRESSURE_STDATM;
  }

  //! Copy constructor.
  Ion5Moment2D_pState(const Ion5Moment2D_pState &W) {
    d = W.d; v = W.v; p = W.p;
  }

  //! Copy constructor.
  Ion5Moment2D_pState(const Ion5Moment2D_cState &U);

  //! Assignment constructor.
  Ion5Moment2D_pState(const double &rho,
		      const Vector2D &V,
		      const double &pre) {
    d = rho; v = V; p = pre;
  }

  //! Assignment constructor.
  Ion5Moment2D_pState(const double &rho,
		      const double &vx,
		      const double &vy,
		      const double &pre) {
    d = rho; v.x = vx; v.y = vy; p = pre;
  }
  //@}
    
  /* Destructor. */
  // ~Ion5Moment2D_pState(void);
  // Use automatically generated destructor.
  //@}

  //@{ @name Set static variables.
  //! Set ion constants.
  void setion(void);
  void setion(char *string_ptr);

  //! Set neutral gas type.
  void setgas(void);
  void setgas(char *string_ptr);
  //@}

  //@{ @name State functions.
  //! Number density.
  double n(void);
  double n(void) const;

  //! Temperature.
  double T(void);
  double T(void) const;

  //! Specific internal energy.
  double e(void);
  double e(void) const;

  //! Total energy.
  double E(void);
  double E(void) const;

  //! Specific enthalpy.
  double h(void);
  double h(void) const;

  //! Total enthalpy.
  double H(void);
  double H(void) const;

  //! Sound speed.
  double a(void);
  double a(void) const;

  //! Sound speed squared.
  double a2(void);
  double a2(void) const;

  //! Mach number.
  double M(void);
  double M(void) const;

  //! Specific entropy.
  double s(void);
  double s(void) const;

  //! Momentum.
  Vector2D dv(void);
  Vector2D dv(void) const;
  double dv(const Vector2D &n);
  double dv(const Vector2D &n) const;

  //! Stagnation temperature.
  double To(void);
  double To(void) const;

  //! Stagnation pressure.
  double po(void);
  double po(void) const;

  //! Stagnation sound speed.
  double ao(void);
  double ao(void) const;

  //! Stagnation enthalpy.
  double ho(void);
  double ho(void) const;

  //! Ion-neutral collision frequency.
  double nu(const Euler2D_pState &Wneut);
  double nu(const Euler2D_pState &Wneut) const;
  friend double nu(const Ion5Moment2D_pState &W, const Euler2D_pState &Wneut);

  //! Neutral gas particle mass.
  double mn(void);
  double mn(void) const;
  friend double mn(const Ion5Moment2D_pState &W);
  //@}

  //@{ @name Conserved solution state.
  Ion5Moment2D_cState U(void);
  Ion5Moment2D_cState U(void) const;
  Ion5Moment2D_cState U(const Ion5Moment2D_pState &W);
  friend Ion5Moment2D_cState U(const Ion5Moment2D_pState &W);
  //@}

  //@{ @name Solution flux and Jacobian (x-direction).
  Ion5Moment2D_cState F(void);
  Ion5Moment2D_cState F(void) const;
  Ion5Moment2D_cState F(const Ion5Moment2D_pState &W);
  friend Ion5Moment2D_cState F(const Ion5Moment2D_pState &W);

  Ion5Moment2D_cState Fx(void);
  Ion5Moment2D_cState Fx(void) const;
  Ion5Moment2D_cState Fx(const Ion5Moment2D_pState &W);
  friend Ion5Moment2D_cState Fx(const Ion5Moment2D_pState &W);
  //@}

  //@{ @name Solution flux and Jacobian (y-direction).
  Ion5Moment2D_cState Fy(void);
  Ion5Moment2D_cState Fy(void) const;
  Ion5Moment2D_cState Fy(const Ion5Moment2D_pState &W);
  friend Ion5Moment2D_cState Fy(const Ion5Moment2D_pState &W);
  //@}

  //@{ @name Solution flux and Jacobian (n-direction).
  Ion5Moment2D_cState Fn(void);
  Ion5Moment2D_cState Fn(void) const;
  Ion5Moment2D_cState Fn(const Ion5Moment2D_pState &W);
  friend Ion5Moment2D_cState Fn(const Ion5Moment2D_pState &W);
  //@}

  //@{ @name Eigenvalue(s) (x-direction).
  Ion5Moment2D_pState lambda(void);
  Ion5Moment2D_pState lambda(void) const;
  Ion5Moment2D_pState lambda(const Ion5Moment2D_pState &W);
  friend Ion5Moment2D_pState lambda(const Ion5Moment2D_pState &W);
  double lambda(int index);
  double lambda(int index) const;
  friend double lambda(const Ion5Moment2D_pState &W, int index);

  Ion5Moment2D_pState lambda_x(void);
  Ion5Moment2D_pState lambda_x(void) const;
  Ion5Moment2D_pState lambda_x(const Ion5Moment2D_pState &W);
  friend Ion5Moment2D_pState lambda_x(const Ion5Moment2D_pState &W);
  double lambda_x(int index);
  double lambda_x(int index) const;
  friend double lambda_x(const Ion5Moment2D_pState &W, int index);
  //@}

  //@{ @name Eigenvalue(s) (y-direction).
  Ion5Moment2D_pState lambda_y(void);
  Ion5Moment2D_pState lambda_y(void) const;
  Ion5Moment2D_pState lambda_y(const Ion5Moment2D_pState &W);
  friend Ion5Moment2D_pState lambda_y(const Ion5Moment2D_pState &W);
  double lambda_y(int index);
  double lambda_y(int index) const;
  friend double lambda_y(const Ion5Moment2D_pState &W, int index);
  //@}

  //@{ @name Primitive right eigenvector (x-direction).
  Ion5Moment2D_pState rp(int index);
  Ion5Moment2D_pState rp(int index) const;
  friend Ion5Moment2D_pState rp(const Ion5Moment2D_pState &W, int index);

  Ion5Moment2D_pState rp_x(int index);
  Ion5Moment2D_pState rp_x(int index) const;
  friend Ion5Moment2D_pState rp_x(const Ion5Moment2D_pState &W, int index);
  //@}

  //@{ @name Primitive right eigenvector (y-direction).
  Ion5Moment2D_pState rp_y(int index);
  Ion5Moment2D_pState rp_y(int index) const;
  friend Ion5Moment2D_pState rp_y(const Ion5Moment2D_pState &W, int index);
  //@}

  //@{ @name Conserved right eigenvector (x-direction).
  Ion5Moment2D_cState rc(int index);
  Ion5Moment2D_cState rc(int index) const;
  friend Ion5Moment2D_cState rc(const Ion5Moment2D_pState &W, int index);

  Ion5Moment2D_cState rc_x(int index);
  Ion5Moment2D_cState rc_x(int index) const;
  friend Ion5Moment2D_cState rc_x(const Ion5Moment2D_pState &W, int index);
  //@}

  //@{ @name Conserved right eigenvector (y-direction).
  Ion5Moment2D_cState rc_y(int index);
  Ion5Moment2D_cState rc_y(int index) const;
  friend Ion5Moment2D_cState rc_y(const Ion5Moment2D_pState &W, int index);
  //@}

  //@{ @name Primitive left eigenvector (x-direction).
  Ion5Moment2D_pState lp(int index);
  Ion5Moment2D_pState lp(int index) const;
  friend Ion5Moment2D_pState lp(const Ion5Moment2D_pState &W, int index);

  Ion5Moment2D_pState lp_x(int index);
  Ion5Moment2D_pState lp_x(int index) const;
  friend Ion5Moment2D_pState lp_x(const Ion5Moment2D_pState &W, int index);
  //@}

  //@{ @name Primitive left eigenvector (y-direction).
  Ion5Moment2D_pState lp_y(int index);
  Ion5Moment2D_pState lp_y(int index) const;
  friend Ion5Moment2D_pState lp_y(const Ion5Moment2D_pState &W, int index);
  //@}

  //@{ @name Axisymmetric flow source vector and Jacobian.
  Ion5Moment2D_cState Sa(const Vector2D &X);
  Ion5Moment2D_cState Sa(const Vector2D &X) const;
  friend Ion5Moment2D_cState Sa(const Ion5Moment2D_pState &W, const Vector2D &X);
  //@}

  //@{ @name Source vector (axisymmetric terms).
  void dSadU(DenseMatrix &dSadU, const Vector2D &X);
  void dSadU(DenseMatrix &dSadU, const Vector2D &X) const;
  friend void dSadU(DenseMatrix &dSadU, const Ion5Moment2D_pState &W, const Vector2D &X);
  //@}

  //@{ @name Electric field source vector and Jacobian.
  Ion5Moment2D_cState Se(const Vector2D &E);
  Ion5Moment2D_cState Se(const Vector2D &E) const;
  friend Ion5Moment2D_cState Se(const Ion5Moment2D_pState &W, const Vector2D &E);

  void dSedU(DenseMatrix &dSedU, const Vector2D &E);
  void dSedU(DenseMatrix &dSedU, const Vector2D &E) const;
  friend void dSedU(DenseMatrix &dSedU, const Ion5Moment2D_pState &W, const Vector2D &E);
  //@}

  //@{ @name Ion-neutral collisions source vector and Jacobian.
  Ion5Moment2D_cState Sc(const Euler2D_pState &Wneut);
  Ion5Moment2D_cState Sc(const Euler2D_pState &Wneut) const;
  friend Ion5Moment2D_cState Sc(const Ion5Moment2D_pState &W, const Euler2D_pState &Wneut);

  void dScdU(DenseMatrix &dScdU, const Euler2D_pState &Wneut);
  void dScdU(DenseMatrix &dScdU, const Euler2D_pState &Wneut) const;
  friend void dScdU(DenseMatrix &dScdU, const Ion5Moment2D_pState &W, const Euler2D_pState &Wneut);
  //@}

  /* Assignment operator. */
  // Ion5Moment2D_pState operator = (const Ion5Moment2D_pState &W);
  // Use automatically generated assignment operator.

  //@{ @name Index operator.
  double &operator[](int index) {
    assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
    switch(index) {
    case 1 :
      return (d);
    case 2 :
      return (v.x);
    case 3 :
      return (v.y);
    case 4 :
      return (p);
    default:
      return (d);
    };
  }
    
  const double &operator[](int index) const {
    assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
    switch(index) {
    case 1 :
      return (d);
    case 2 :
      return (v.x);
    case 3 :
      return (v.y);
    case 4 :
      return (p);
    default:
      return (d);
    };
  }
  //@}

  //@{ @name Binary arithmetic operators.
  friend Ion5Moment2D_pState operator +(const Ion5Moment2D_pState &W1, const Ion5Moment2D_pState &W2);
  friend Ion5Moment2D_pState operator -(const Ion5Moment2D_pState &W1, const Ion5Moment2D_pState &W2);
  friend double operator *(const Ion5Moment2D_pState &W1, const Ion5Moment2D_pState &W2);
  friend Ion5Moment2D_pState operator *(const Ion5Moment2D_pState &W, const double &a);
  friend Ion5Moment2D_pState operator *(const double &a, const Ion5Moment2D_pState &W);
  friend Ion5Moment2D_pState operator /(const Ion5Moment2D_pState &W, const double &a);
  friend Ion5Moment2D_pState operator ^(const Ion5Moment2D_pState &W1, const Ion5Moment2D_pState &W2);
  //@}

  //@{ @name Unary arithmetic operators.
  friend Ion5Moment2D_pState operator +(const Ion5Moment2D_pState &W);
  friend Ion5Moment2D_pState operator -(const Ion5Moment2D_pState &W);
  //@}

  //@{ @name Shortcut arithmetic operators.
  friend Ion5Moment2D_pState &operator +=(Ion5Moment2D_pState &W1, const Ion5Moment2D_pState &W2);
  friend Ion5Moment2D_pState &operator -=(Ion5Moment2D_pState &W1, const Ion5Moment2D_pState &W2);
  //@}

  //@{ @name Relational operators.
  friend int operator ==(const Ion5Moment2D_pState &W1, const Ion5Moment2D_pState &W2);
  friend int operator !=(const Ion5Moment2D_pState &W1, const Ion5Moment2D_pState &W2);
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const Ion5Moment2D_pState &W);
  friend istream &operator >> (istream &in_file,  Ion5Moment2D_pState &W);
  //@}

};

/*!
 * Class: Ion5Moment2D_cState
 *
 * @brief Conserved variable solution state class definition for an
 *        inviscid Ion5moment flow.
 *
 * Conserved variable solution state class definition for an inviscid
 * Ion5moment flow.
 *
 * \verbatim
 * Member functions
 *     d        -- Return ion mass density.
 *     dv       -- Return momentum.
 *     E        -- Return total energy.
 *     g        -- Return ion specific heat ratio.
 *     gm1      -- Return g-1.
 *     gm1i     -- Return 1/(g-1).
 *     R        -- Return ion gas constant.
 *     q        -- Return ion charge.
 *     m        -- Return ion particle mass.
 *     ion      -- Return ion type indicator.
 *     gas      -- Return neutral gas type indicator.
 *     setion   -- Set ion constants.
 *     setgas   -- Set neutral gas type.
 *     n        -- Return ion number density.
 *     v        -- Return ion flow velocity.
 *     p        -- Return ion pressure.
 *     T        -- Return ion temperature.
 *     e        -- Return specific internal energy.
 *     h        -- Return specific enthalpy.
 *     H        -- Return total enthalpy.
 *     a        -- Return sound speed.
 *     a2       -- Return sound speed square.
 *     M        -- Return Mach number.
 *     s        -- Return specific entropy.
 *     To       -- Return stagnation temperature.
 *     po       -- Return stagnation pressure.
 *     ao       -- Return stagnation sound speed.
 *     ho       -- Return stagnation enthalpy.
 *     nu       -- Return ion-neutral collision frequency.
 *     mn       -- Return neutral gas particle mass.
 *     W        -- Return primitive solution state.
 *     F        -- Return x-direction solution flux.
 *     Fx       -- Return x-direction solution flux.
 *     Fy       -- Return y-direction solution flux.
 *     Fn       -- Return n-direction solution flux.
 *     Sa       -- Return axisymmetric source term vector.
 *     Se       -- Return Coulomb force source vector produced by
 *                 external electric field.
 *     Sc       -- Return source terms associated with ion-neutral
 *                 collision processes.
 *     dSadU    -- Return axisymmetric source term Jacobian.
 *     dSedU    -- Return Coulomb force source Jacobian.
 *     dScdU    -- Return source Jacobian for ion-neutral collision
 *                 processes.
 *
 * Member operators
 *      U -- a primitive solution state
 *      c -- a scalar (double)
 *
 * U = U;
 * c = U[i];
 * U = U + U;
 * U = U - U;
 * c = U * U; (inner product)
 * U = c * U;
 * U = U * c;
 * U = U / c;
 * U = U ^ U; (a useful product)
 * U = +U;
 * U = -U;
 * U += U;
 * U -= U;
 * U == U;
 * U != U;
 * cout << U; (output function)
 * cin  >> U; (input function)
 * \endverbatim
 */
class Ion5Moment2D_cState{
private:
public:
  //@{ @name Conserved variables and associated constants:
  double           d; //!< Ion mass density.
  Vector2D        dv; //!< Ion momentum.
  double           E; //!< Ion total Energy.
  static double    g; //!< Ion specific heat ratio.
  static double  gm1; //!< g-1
  static double gm1i; //!< 1/(g-1)
  static double    R; //!< Ion gas constant.
  static double    q; //!< Ion charge.
  static double    m; //!< Ion particle mass.
  static int     ion; //!< Ion type indicator.
  static int     gas; //!< Neutral gas type indicator.
  //@}

  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  Ion5Moment2D_cState(void) {
    d = DENSITY_STDATM; dv.zero(); E = PRESSURE_STDATM/(GAMMA_AIR-ONE);
  }

  //! Copy constructor.
  Ion5Moment2D_cState(const Ion5Moment2D_cState &U) {
    d = U.d; dv = U.dv; E = U.E;
  }

  //! Copy constructor.
  Ion5Moment2D_cState(const Ion5Moment2D_pState &W);

  //! Assignment constructor.
  Ion5Moment2D_cState(const double &rho,
		      const Vector2D &rhoV,
		      const double &Etotal) {
    d = rho; dv = rhoV; E = Etotal;
  }

  //! Assignment constructor.
  Ion5Moment2D_cState(const double &rho,
		      const double &rhovx,
		      const double &rhovy,
		      const double &Etotal) {
    d = rho; dv.x = rhovx; dv.y = rhovy; E = Etotal;
  }
    
  /* Destructor. */
  // ~Ion5Moment2D_cState(void);
  // Use automatically generated destructor.
  //@}

  //@{ @name Set static variables.
  //! Set ion constants.
  void setion(void);
  void setion(char *string_ptr);

  //! Set neutral gas type.
  void setgas(void);
  void setgas(char *string_ptr);
  //@}

  //@{ @name State functions.
  //! Number density.
  double n(void);
  double n(void) const;

  //! Flow velocity.
  Vector2D v(void);
  Vector2D v(void) const;
  double v(const Vector2D &n);
  double v(const Vector2D &n) const;

  //! Pressure.
  double p(void);
  double p(void) const;

  //! Temperature.
  double T(void);
  double T(void) const;

  //! Specific internal energy.
  double e(void);
  double e(void) const;

  //! Specific enthalpy.
  double h(void);
  double h(void) const;

  //! Total enthalpy.
  double H(void);
  double H(void) const;

  //! Sound speed.
  double a(void);
  double a(void) const;

  //! Sound speed squared.
  double a2(void);
  double a2(void) const;

  //! Mach number.
  double M(void);
  double M(void) const;

  //! Specific entropy.
  double s(void);
  double s(void) const;

  //! Stagnation temperature.
  double To(void);
  double To(void) const;

  //! Stagnation pressure.
  double po(void);
  double po(void) const;

  //! Stagnation sound speed.
  double ao(void);
  double ao(void) const;

  //! Stagnation enthalpy.
  double ho(void);
  double ho(void) const;

  //! Ion-neutral collision frequency.
  double nu(const Euler2D_pState &Wneut);
  double nu(const Euler2D_pState &Wneut) const;
  friend double nu(const Ion5Moment2D_cState &W, const Euler2D_pState &Wneut);

  //! Neutral gas particle mass.
  double mn(void);
  double mn(void) const;
  friend double mn(const Ion5Moment2D_cState &U);
  //@}

  //@{ @name Conserved solution state.
  Ion5Moment2D_pState W(void);
  Ion5Moment2D_pState W(void) const;
  Ion5Moment2D_pState W(const Ion5Moment2D_cState &U);
  friend Ion5Moment2D_pState W(const Ion5Moment2D_cState &U);
  //@}

  //@{ @name Solution flux and Jacobian (x-direction).
  Ion5Moment2D_cState F(void);
  Ion5Moment2D_cState F(void) const;
  Ion5Moment2D_cState F(const Ion5Moment2D_cState &U);
  friend Ion5Moment2D_cState F(const Ion5Moment2D_cState &U);

  Ion5Moment2D_cState Fx(void);
  Ion5Moment2D_cState Fx(void) const;
  Ion5Moment2D_cState Fx(const Ion5Moment2D_cState &U);
  friend Ion5Moment2D_cState Fx(const Ion5Moment2D_cState &U);
  //@}

  //@{ @name Solution flux and Jacobian (y-direction).
  Ion5Moment2D_cState Fy(void);
  Ion5Moment2D_cState Fy(void) const;
  Ion5Moment2D_cState Fy(const Ion5Moment2D_cState &U);
  friend Ion5Moment2D_cState Fy(const Ion5Moment2D_cState &U);
  //@}

  //@{ @name Solution flux and Jacobian (n-direction).
  Ion5Moment2D_cState Fn(void);
  Ion5Moment2D_cState Fn(void) const;
  Ion5Moment2D_cState Fn(const Ion5Moment2D_cState &U);
  friend Ion5Moment2D_cState Fn(const Ion5Moment2D_cState &U);
  //@}

  //@{ @name Axisymmetric flow source vector and Jacobian.
  Ion5Moment2D_cState Sa(const Vector2D &X);
  Ion5Moment2D_cState Sa(const Vector2D &X) const;
  friend Ion5Moment2D_cState Sa(const Ion5Moment2D_cState &U, const Vector2D &X);

  void dSadU(DenseMatrix &dSadU, const Vector2D &X);
  void dSadU(DenseMatrix &dSadU, const Vector2D &X) const;
  friend void dSadU(DenseMatrix &dSadU, const Ion5Moment2D_cState &U, const Vector2D &X);
  //@}

  //@{ @name Electric field source vector and Jacobian.
  Ion5Moment2D_cState Se(const Vector2D &E);
  Ion5Moment2D_cState Se(const Vector2D &E) const;
  friend Ion5Moment2D_cState Se(const Ion5Moment2D_cState &U, const Vector2D &E);

  void dSedU(DenseMatrix &dSedU, const Vector2D &E);
  void dSedU(DenseMatrix &dSedU, const Vector2D &E) const;
  friend void dSedU(DenseMatrix &dSedU, const Ion5Moment2D_cState &U, const Vector2D &E);
  //@}

  //@{ @name Ion-neutral collisions source vector and Jacobian.
  Ion5Moment2D_cState Sc(const Euler2D_pState &Wneut);
  Ion5Moment2D_cState Sc(const Euler2D_pState &Wneut) const;
  friend Ion5Moment2D_cState Sc(const Ion5Moment2D_cState &U, const Euler2D_pState &Wneut);

  void dScdU(DenseMatrix &dScdU, const Euler2D_pState &Wneut);
  void dScdU(DenseMatrix &dScdU, const Euler2D_pState &Wneut) const;
  friend void dScdU(DenseMatrix &dScdU, const Ion5Moment2D_cState &U, const Euler2D_pState &Wneut);
  //@}

  /* Assignment operator. */
  // Ion5Moment2D_cState operator = (const Ion5Moment2D_cState &U);
  // Use automatically generated assignment operator.

  //@{ @name Index operator.
  double &operator[](int index) {
    assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
    switch(index) {
    case 1 :
      return (d);
    case 2 :
      return (dv.x);
    case 3 :
      return (dv.y);
    case 4 :
      return (E);
    default:
      return (d);
    };
  }
    
  const double &operator[](int index) const {
    assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
    switch(index) {
    case 1 :
      return (d);
    case 2 :
      return (dv.x);
    case 3 :
      return (dv.y);
    case 4 :
      return (E);
    default:
      return (d);
    };
  }
  //@}

  //@{ @name Binary arithmetic operators.
  friend Ion5Moment2D_cState operator +(const Ion5Moment2D_cState &U1, const Ion5Moment2D_cState &U2);
  friend Ion5Moment2D_cState operator -(const Ion5Moment2D_cState &U1, const Ion5Moment2D_cState &U2);
  friend double operator *(const Ion5Moment2D_cState &U1, const Ion5Moment2D_cState &U2);
  friend Ion5Moment2D_cState operator *(const Ion5Moment2D_cState &U, const double &a);
  friend Ion5Moment2D_cState operator *(const double &a, const Ion5Moment2D_cState &U);
  friend Ion5Moment2D_cState operator /(const Ion5Moment2D_cState &U, const double &a);
  friend Ion5Moment2D_cState operator ^(const Ion5Moment2D_cState &U1, const Ion5Moment2D_cState &U2);
  //@}

  //@{ @name Unary arithmetic operators.
  friend Ion5Moment2D_cState operator +(const Ion5Moment2D_cState &U);
  friend Ion5Moment2D_cState operator -(const Ion5Moment2D_cState &U);
  //@}

  //@{ @name Shortcut arithmetic operators.
  friend Ion5Moment2D_cState &operator +=(Ion5Moment2D_cState &U1, const Ion5Moment2D_cState &U2);
  friend Ion5Moment2D_cState &operator -=(Ion5Moment2D_cState &U1, const Ion5Moment2D_cState &U2);
  //@}

  //@{ @name Relational operators.
  friend int operator ==(const Ion5Moment2D_cState &U1, const Ion5Moment2D_cState &U2);
  friend int operator !=(const Ion5Moment2D_cState &U1, const Ion5Moment2D_cState &U2);
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const Ion5Moment2D_cState &U);
  friend istream &operator >> (istream &in_file,  Ion5Moment2D_cState &U);
  //@}

};

/*************************************************************
 * Ion5Moment2D_pState::setion -- Assign ion constants.      *
 *************************************************************/
inline void Ion5Moment2D_pState::setion(void) {
    g = GAMMA_H;
    R = R_UNIVERSAL/(MOLE_WT_H*MILLI);
    gm1 = g - ONE;
    gm1i = ONE/gm1;
    q = E_CHARGE;
    m = MOLE_WT_H/(AVOGADRO*THOUSAND);
    ion = ION_H;
}

inline void Ion5Moment2D_pState::setion(char *string_ptr) {
   if (strcmp(string_ptr, "H+") == 0) {
     g = GAMMA_H;
     R = R_UNIVERSAL/(MOLE_WT_H*MILLI);
     q = E_CHARGE;
     m = MOLE_WT_H/(AVOGADRO*THOUSAND);
     ion = ION_H;  
   } else if (strcmp(string_ptr, "HE+") == 0) {
     g = GAMMA_HE;
     R = R_UNIVERSAL/(MOLE_WT_HE*MILLI);
     q = E_CHARGE;
     m = MOLE_WT_HE/(AVOGADRO*THOUSAND);
     ion = ION_HE;
   } else if (strcmp(string_ptr, "O+") == 0) {
     g = GAMMA_O;
     R = R_UNIVERSAL/(MOLE_WT_O*MILLI);
     q = E_CHARGE;
     m = MOLE_WT_O/(AVOGADRO*THOUSAND);
     ion = ION_O;
   } else if (strcmp(string_ptr, "ION_LOW_MASS") == 0) {
     g = GAMMA_ION_LOW_MASS;
     R = R_UNIVERSAL/(MOLE_WT_ION_LOW_MASS*MILLI);
     q = E_CHARGE;
     m = MOLE_WT_ION_LOW_MASS/(AVOGADRO*THOUSAND);
     ion = ION_LOW_MASS;
   } else if (strcmp(string_ptr, "ION_MED_MASS") == 0) {
     g = GAMMA_ION_MED_MASS;
     R = R_UNIVERSAL/(MOLE_WT_ION_MED_MASS*MILLI);
     q = E_CHARGE;
     m = MOLE_WT_ION_MED_MASS/(AVOGADRO*THOUSAND);
     ion = ION_MED_MASS;
   } else if (strcmp(string_ptr, "ION_HI_MASS") == 0) {
     g = GAMMA_ION_HI_MASS;
     R = R_UNIVERSAL/(MOLE_WT_ION_HI_MASS*MILLI);
     q = E_CHARGE;
     m = MOLE_WT_ION_HI_MASS/(AVOGADRO*THOUSAND);
     ion = ION_HI_MASS;
   } else if (strcmp(string_ptr, "ION_HUGE_MASS") == 0) {
     g = GAMMA_ION_HUGE_MASS;
     R = R_UNIVERSAL/(MOLE_WT_ION_HUGE_MASS*MILLI);
     q = NINE*E_CHARGE;
     m = MOLE_WT_ION_HUGE_MASS/(AVOGADRO*THOUSAND);
     ion = ION_HUGE_MASS;
   } else if (strcmp(string_ptr, "ION_HUGE_MASS2") == 0) {
     g = GAMMA_ION_HUGE_MASS2;
     R = R_UNIVERSAL/(MOLE_WT_ION_HUGE_MASS2*MILLI);
     q = ONE*E_CHARGE;
     m = MOLE_WT_ION_HUGE_MASS2/(AVOGADRO*THOUSAND);
     ion = ION_HUGE_MASS2;
   } else {
     g = GAMMA_H;
     R = R_UNIVERSAL/(MOLE_WT_H*MILLI);
     q = E_CHARGE;
     m = MOLE_WT_H/(AVOGADRO*THOUSAND);
     ion = ION_H;
   } /* endif */
   gm1 = g - ONE;
   gm1i = ONE/gm1;
}

/*************************************************************
 * Ion5Moment2D_pState::setgas -- Set neutral gas type.      *
 *************************************************************/
inline void Ion5Moment2D_pState::setgas(void) {
   gas = GAS_AIR;
}

inline void Ion5Moment2D_pState::setgas(char *string_ptr) {
   if (strcmp(string_ptr, "AIR") == 0) {
     gas = GAS_AIR;
   } else if (strcmp(string_ptr, "A") == 0) {
     gas = GAS_A;
   } else if (strcmp(string_ptr, "CO") == 0) {
     gas = GAS_CO;
   } else if (strcmp(string_ptr, "CO2") == 0) {
     gas = GAS_CO2;
   } else if (strcmp(string_ptr, "CH4") == 0) {
     gas = GAS_CH4;
   } else if (strcmp(string_ptr, "H") == 0) {
     gas = GAS_H;
   } else if (strcmp(string_ptr, "H2") == 0) {
     gas = GAS_H2;
   } else if (strcmp(string_ptr, "HE") == 0) {
     gas = GAS_HE;
   } else if (strcmp(string_ptr, "H2O") == 0) {
     gas = GAS_H2O;
   } else if (strcmp(string_ptr, "N2") == 0) {
     gas = GAS_N2;
   } else if (strcmp(string_ptr, "O") == 0) {
     gas = GAS_O;
   } else if (strcmp(string_ptr, "O2") == 0) {
     gas = GAS_O2;
   } else if (strcmp(string_ptr, "e") == 0) {
     gas = GAS_e;
   } else if (strcmp(string_ptr, "AP-HTPB") == 0) {
     gas = GAS_APHTPB;
   } else {
     gas = GAS_AIR;
   } /* endif */
}

/*************************************************************
 * Ion5Moment2D_pState::n -- Number density.                 *
 *************************************************************/
inline double Ion5Moment2D_pState::n(void) {
    return (d/m);
}

inline double Ion5Moment2D_pState::n(void) const {
    return (d/m);
}

/*************************************************************
 * Ion5Moment2D_pState::T -- Temperature.                    *
 *************************************************************/
inline double Ion5Moment2D_pState::T(void) {
    return (p/(d*R));
}

inline double Ion5Moment2D_pState::T(void) const {
    return (p/(d*R));
}

/*************************************************************
 * Ion5Moment2D_pState::e -- Specific internal energy.       *
 *************************************************************/
inline double Ion5Moment2D_pState::e(void) {
    return (p/(gm1*d));
}

inline double Ion5Moment2D_pState::e(void) const {
    return (p/(gm1*d));
}

/*************************************************************
 * Ion5Moment2D_pState::E -- Total energy.                   *
 *************************************************************/
inline double Ion5Moment2D_pState::E(void) {
    return (p*gm1i + HALF*d*v.sqr());
}

inline double Ion5Moment2D_pState::E(void) const {
    return (p*gm1i + HALF*d*v.sqr());
}

/*************************************************************
 * Ion5Moment2D_pState::h -- Specific enthalpy.              *
 *************************************************************/
inline double Ion5Moment2D_pState::h(void) {
    return (g*p/(gm1*d) + HALF*v.sqr());
}

inline double Ion5Moment2D_pState::h(void) const {
    return (g*p/(gm1*d) + HALF*v.sqr());
}

/*************************************************************
 * Ion5Moment2D_pState::H -- Total enthalpy.                 *
 *************************************************************/
inline double Ion5Moment2D_pState::H(void) {
    return (g*gm1i*p + HALF*d*v.sqr());
}

inline double Ion5Moment2D_pState::H(void) const {
    return (g*gm1i*p + HALF*d*v.sqr());
}

/*************************************************************
 * Ion5Moment2D_pState::a -- Sound speed.                    *
 *************************************************************/
inline double Ion5Moment2D_pState::a(void) {
    return (sqrt(g*p/d));
}

inline double Ion5Moment2D_pState::a(void) const {
    return (sqrt(g*p/d));
}

/*************************************************************
 * Ion5Moment2D_pState::a2 -- Sound speed squared.           *
 *************************************************************/
inline double Ion5Moment2D_pState::a2(void) {
    return (g*p/d);
}

inline double Ion5Moment2D_pState::a2(void) const {
    return (g*p/d);
}

/*************************************************************
 * Ion5Moment2D_pState::M -- Mach number.                    *
 *************************************************************/
inline double Ion5Moment2D_pState::M(void) {
    return (abs(v)/sqrt(g*p/d));
}

inline double Ion5Moment2D_pState::M(void) const {
    return (abs(v)/sqrt(g*p/d));
}

/*************************************************************
 * Ion5Moment2D_pState::s -- Specific entropy.               *
 *************************************************************/
inline double Ion5Moment2D_pState::s(void) {
    return (R*gm1i*log(p/pow(d, g)));
}

inline double Ion5Moment2D_pState::s(void) const {
    return (R*gm1i*log(p/pow(d, g)));
}

/*************************************************************
 * Ion5Moment2D_pState::dv -- Momentum.                      *
 *************************************************************/
inline Vector2D Ion5Moment2D_pState::dv(void) {
    return (d*v);
}

inline Vector2D Ion5Moment2D_pState::dv(void) const {
    return (d*v);
}

inline double Ion5Moment2D_pState::dv(const Vector2D &n) {
    return (d*(v*n));
}

inline double Ion5Moment2D_pState::dv(const Vector2D &n) const {
    return (d*(v*n));
}

/*************************************************************
 * Ion5Moment2D_pState::To -- Stagnation temperature.        *
 *************************************************************/
inline double Ion5Moment2D_pState::To(void) {
    return ((p/(d*R))*(ONE+HALF*gm1*v.sqr()/(g*p/d)));
}

inline double Ion5Moment2D_pState::To(void) const {
    return ((p/(d*R))*(ONE+HALF*gm1*v.sqr()/(g*p/d)));
}

/*************************************************************
 * Ion5Moment2D_pState::po -- Stagnation pressure.           *
 *************************************************************/
inline double Ion5Moment2D_pState::po(void) {
    return (p*pow(ONE+HALF*gm1*v.sqr()/(g*p/d), g*gm1i));
}

inline double Ion5Moment2D_pState::po(void) const {
    return (p*pow(ONE+HALF*gm1*v.sqr()/(g*p/d), g*gm1i));
}

/*************************************************************
 * Ion5Moment2D_pState::ao -- Stagnation sound speed.        *
 *************************************************************/
inline double Ion5Moment2D_pState::ao(void) {
    return (sqrt((g*p/d)*(ONE+HALF*gm1*v.sqr()/(g*p/d))));
}

inline double Ion5Moment2D_pState::ao(void) const {
    return (sqrt((g*p/d)*(ONE+HALF*gm1*v.sqr()/(g*p/d))));
}

/*************************************************************
 * Ion5Moment2D_pState::ho -- Stagnation enthalpy.           *
 *************************************************************/
inline double Ion5Moment2D_pState::ho(void) {
    return ((g*p/(gm1*d) + HALF*v.sqr())*(ONE+HALF*gm1*v.sqr()/(g*p/d)));
}

inline double Ion5Moment2D_pState::ho(void) const {
    return ((g*p/(gm1*d) + HALF*v.sqr())*(ONE+HALF*gm1*v.sqr()/(g*p/d)));
}

/*************************************************************
 * Ion5Moment2D_pState -- Binary arithmetic operators.       *
 *************************************************************/
inline Ion5Moment2D_pState operator +(const Ion5Moment2D_pState &W1, const Ion5Moment2D_pState &W2) {
  return (Ion5Moment2D_pState(W1.d+W2.d,W1.v+W2.v,W1.p+W2.p));
}

inline Ion5Moment2D_pState operator -(const Ion5Moment2D_pState &W1, const Ion5Moment2D_pState &W2) {
  return (Ion5Moment2D_pState(W1.d-W2.d,W1.v-W2.v,W1.p-W2.p));
}

// Inner product operator.
inline double operator *(const Ion5Moment2D_pState &W1, const Ion5Moment2D_pState &W2) {
   return (W1.d*W2.d+W1.v*W2.v+W1.p*W2.p);
}

inline Ion5Moment2D_pState operator *(const Ion5Moment2D_pState &W, const double &a) {
  return (Ion5Moment2D_pState(a*W.d,a*W.v,a*W.p));
}

inline Ion5Moment2D_pState operator *(const double &a, const Ion5Moment2D_pState &W) {
  return (Ion5Moment2D_pState(a*W.d,a*W.v,a*W.p));
}

inline Ion5Moment2D_pState operator /(const Ion5Moment2D_pState &W, const double &a) {
  return (Ion5Moment2D_pState(W.d/a,W.v/a,W.p/a));
}

// My useful solution state product operator.
inline Ion5Moment2D_pState operator ^(const Ion5Moment2D_pState &W1, const Ion5Moment2D_pState &W2) {
   return (Ion5Moment2D_pState(W1.d*W2.d,W1.v.x*W2.v.x,W1.v.y*W2.v.y,W1.p*W2.p));
}

/*************************************************************
 * Ion5Moment2D_pState -- Unary arithmetic operators.        *
 *************************************************************/
inline Ion5Moment2D_pState operator +(const Ion5Moment2D_pState &W) {
  return (Ion5Moment2D_pState(W.d,W.v,W.p));
}

inline Ion5Moment2D_pState operator -(const Ion5Moment2D_pState &W) {
  return (Ion5Moment2D_pState(-W.d,-W.v,-W.p));
}

/*************************************************************
 * Ion5Moment2D_pState -- Shortcut arithmetic operators.     *
 *************************************************************/
inline Ion5Moment2D_pState &operator +=(Ion5Moment2D_pState &W1, const Ion5Moment2D_pState &W2) {
  W1.d += W2.d;
  W1.v += W2.v;
  W1.p += W2.p;
  return (W1);
}

inline Ion5Moment2D_pState &operator -=(Ion5Moment2D_pState &W1, const Ion5Moment2D_pState &W2) {
  W1.d -= W2.d;
  W1.v -= W2.v;
  W1.p -= W2.p;
  return (W1);
}

/*************************************************************
 * Ion5Moment2D_pState -- Relational operators.              *
 *************************************************************/
inline int operator ==(const Ion5Moment2D_pState &W1, const Ion5Moment2D_pState &W2) {
  return (W1.d == W2.d && W1.v == W2.v && W1.p == W2.p);
}

inline int operator !=(const Ion5Moment2D_pState &W1, const Ion5Moment2D_pState &W2) {
  return (W1.d != W2.d || W1.v != W2.v || W1.p != W2.p);
}

/*************************************************************
 * Ion5Moment2D_pState -- Input-output operators.            *
 *************************************************************/
inline ostream &operator << (ostream &out_file, const Ion5Moment2D_pState &W) {
  out_file.setf(ios::scientific);
  out_file << " " << W.d  << " " << W.v.x << " " << W.v.y << " " << W.p;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, Ion5Moment2D_pState &W) {
  in_file.setf(ios::skipws);
  in_file >> W.d >> W.v.x >> W.v.y >> W.p;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/*************************************************************
 * Ion5Moment2D_cState::setion -- Assign ion constants.      *
 *************************************************************/
inline void Ion5Moment2D_cState::setion(void) {
    g = GAMMA_AIR;
    R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    gm1 = g - ONE;
    gm1i = ONE/gm1;
    q = E_CHARGE;
    m = MOLE_WT_H/(AVOGADRO*THOUSAND);
    ion = ION_H;
}

inline void Ion5Moment2D_cState::setion(char *string_ptr) {
   if (strcmp(string_ptr, "H+") == 0) {
     g = GAMMA_H;
     R = R_UNIVERSAL/(MOLE_WT_H*MILLI);
     q = E_CHARGE;
     m = MOLE_WT_H/(AVOGADRO*THOUSAND);
     ion = ION_H;
   } else if (strcmp(string_ptr, "HE+") == 0) {
     g = GAMMA_HE;
     R = R_UNIVERSAL/(MOLE_WT_HE*MILLI);
     q = E_CHARGE;
     m = MOLE_WT_HE/(AVOGADRO*THOUSAND);
     ion = ION_HE;
   } else if (strcmp(string_ptr, "O+") == 0) {
     g = GAMMA_O;
     R = R_UNIVERSAL/(MOLE_WT_O*MILLI);
     q = E_CHARGE;
     m = MOLE_WT_O/(AVOGADRO*THOUSAND);
     ion = ION_O;
   } else if (strcmp(string_ptr, "ION_LOW_MASS") == 0) {
     g = GAMMA_ION_LOW_MASS;
     R = R_UNIVERSAL/(MOLE_WT_ION_LOW_MASS*MILLI);
     q = E_CHARGE;
     m = MOLE_WT_ION_LOW_MASS/(AVOGADRO*THOUSAND);
     ion = ION_LOW_MASS;
   } else if (strcmp(string_ptr, "ION_MED_MASS") == 0) {
     g = GAMMA_ION_MED_MASS;
     R = R_UNIVERSAL/(MOLE_WT_ION_MED_MASS*MILLI);
     q = E_CHARGE;
     m = MOLE_WT_ION_MED_MASS/(AVOGADRO*THOUSAND);
     ion = ION_MED_MASS;
   } else if (strcmp(string_ptr, "ION_HI_MASS") == 0) {
     g = GAMMA_ION_HI_MASS;
     R = R_UNIVERSAL/(MOLE_WT_ION_HI_MASS*MILLI);
     q = E_CHARGE;
     m = MOLE_WT_ION_HI_MASS/(AVOGADRO*THOUSAND);
     ion = ION_HI_MASS;
   } else if (strcmp(string_ptr, "ION_HUGE_MASS") == 0) {
     g = GAMMA_ION_HUGE_MASS;
     R = R_UNIVERSAL/(MOLE_WT_ION_HUGE_MASS*MILLI);
     q = NINE*E_CHARGE;
     m = MOLE_WT_ION_HUGE_MASS/(AVOGADRO*THOUSAND);
     ion = ION_HUGE_MASS;
   } else if (strcmp(string_ptr, "ION_HUGE_MASS2") == 0) {
     g = GAMMA_ION_HUGE_MASS2;
     R = R_UNIVERSAL/(MOLE_WT_ION_HUGE_MASS2*MILLI);
     q = ONE*E_CHARGE;
     m = MOLE_WT_ION_HUGE_MASS2/(AVOGADRO*THOUSAND);
     ion = ION_HUGE_MASS2;
   } else {
     g = GAMMA_H;
     R = R_UNIVERSAL/(MOLE_WT_H*MILLI);
     q = E_CHARGE;
     m = MOLE_WT_H/(AVOGADRO*THOUSAND);
     ion = ION_H;
   } /* endif */
   gm1 = g - ONE;
   gm1i = ONE/gm1;
}

/*************************************************************
 * Ion5Moment2D_cState::setgas -- Set neutral gas type.      *
 *************************************************************/
inline void Ion5Moment2D_cState::setgas(void) {
   gas = GAS_AIR;
}

inline void Ion5Moment2D_cState::setgas(char *string_ptr) {
   if (strcmp(string_ptr, "AIR") == 0) {
     gas = GAS_AIR;
   } else if (strcmp(string_ptr, "A") == 0) {
     gas = GAS_A;
   } else if (strcmp(string_ptr, "CO") == 0) {
     gas = GAS_CO;
   } else if (strcmp(string_ptr, "CO2") == 0) {
     gas = GAS_CO2;
   } else if (strcmp(string_ptr, "CH4") == 0) {
     gas = GAS_CH4;
   } else if (strcmp(string_ptr, "H") == 0) {
     gas = GAS_H;
   } else if (strcmp(string_ptr, "H2") == 0) {
     gas = GAS_H2;
   } else if (strcmp(string_ptr, "HE") == 0) {
     gas = GAS_HE;
   } else if (strcmp(string_ptr, "H2O") == 0) {
     gas = GAS_H2O;
   } else if (strcmp(string_ptr, "N2") == 0) {
     gas = GAS_N2;
   } else if (strcmp(string_ptr, "O") == 0) {
     gas = GAS_O;
   } else if (strcmp(string_ptr, "O2") == 0) {
     gas = GAS_O2;
   } else if (strcmp(string_ptr, "e") == 0) {
     gas = GAS_e;
   } else if (strcmp(string_ptr, "AP-HTPB") == 0) {
     gas = GAS_APHTPB;
   } else {
     gas = GAS_AIR;
   } /* endif */
}

/*************************************************************
 * Ion5Moment2D_cState::n -- Number density.                 *
 *************************************************************/
inline double Ion5Moment2D_cState::n(void) {
    return (d/m);
}

inline double Ion5Moment2D_cState::n(void) const {
    return (d/m);
}

/*************************************************************
 * Ion5Moment2D_cState::v -- Flow velocity.                  *
 *************************************************************/
inline Vector2D Ion5Moment2D_cState::v(void) {
    return (dv/d);
}

inline Vector2D Ion5Moment2D_cState::v(void) const {
    return (dv/d);
}

inline double Ion5Moment2D_cState::v(const Vector2D &n) {
    return ((dv*n)/d);
}

inline double Ion5Moment2D_cState::v(const Vector2D &n) const {
    return ((dv*n)/d);
}

/*************************************************************
 * Ion5Moment2D_cState::p -- Pressure.                       *
 *************************************************************/
inline double Ion5Moment2D_cState::p(void) {
    return (gm1*(E - HALF*dv.sqr()/d));
}

inline double Ion5Moment2D_cState::p(void) const {
    return (gm1*(E - HALF*dv.sqr()/d));
}

/*************************************************************
 * Ion5Moment2D_cState::T -- Temperature.                    *
 *************************************************************/
inline double Ion5Moment2D_cState::T(void) {
    return (gm1*(E - HALF*dv.sqr()/d)/(d*R));
}

inline double Ion5Moment2D_cState::T(void) const {
    return (gm1*(E - HALF*dv.sqr()/d)/(d*R));
}

/*************************************************************
 * Ion5Moment2D_cState::e -- Specific internal energy.       *
 *************************************************************/
inline double Ion5Moment2D_cState::e(void) {
    return (E/d - HALF*dv.sqr()/sqr(d));
}

inline double Ion5Moment2D_cState::e(void) const {
    return (E/d - HALF*dv.sqr()/sqr(d));
}

/*************************************************************
 * Ion5Moment2D_cState::h -- Specific enthalpy.              *
 *************************************************************/
inline double Ion5Moment2D_cState::h(void) {
    return (g*E/d - gm1*HALF*dv.sqr()/sqr(d));
}

inline double Ion5Moment2D_cState::h(void) const {
    return (g*E/d - gm1*HALF*dv.sqr()/sqr(d));
}

/*************************************************************
 * Ion5Moment2D_cState::H -- Total enthalpy.                 *
 *************************************************************/
inline double Ion5Moment2D_cState::H(void) {
     return (g*E - gm1*HALF*dv.sqr()/d);
}

inline double Ion5Moment2D_cState::H(void) const {
     return (g*E - gm1*HALF*dv.sqr()/d);
}

/*************************************************************
 * Ion5Moment2D_cState::a -- Sound speed.                    *
 *************************************************************/
inline double Ion5Moment2D_cState::a(void) {
    return (sqrt(g*gm1*(E/d - HALF*dv.sqr()/sqr(d))));
}

inline double Ion5Moment2D_cState::a(void) const {
    return (sqrt(g*gm1*(E/d - HALF*dv.sqr()/sqr(d))));
}

/*************************************************************
 * Ion5Moment2D_cState::a2 -- Sound speed squared.           *
 *************************************************************/
inline double Ion5Moment2D_cState::a2(void) {
    return (g*gm1*(E/d - HALF*dv.sqr()/sqr(d)));
}

inline double Ion5Moment2D_cState::a2(void) const {
    return (g*gm1*(E/d - HALF*dv.sqr()/sqr(d)));
}

/*************************************************************
 * Ion5Moment2D_cState::M -- Mach number.                    *
 *************************************************************/
inline double Ion5Moment2D_cState::M(void) {
    return (abs(dv)/(d*sqrt(g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))));
}

inline double Ion5Moment2D_cState::M(void) const {
    return (abs(dv)/(d*sqrt(g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))));
}

/*************************************************************
 * Ion5Moment2D_cState::s -- Specific entropy.               *
 *************************************************************/
inline double Ion5Moment2D_cState::s(void) {
    return (R*gm1i*log(gm1*(E - HALF*dv.sqr()/d)/pow(d, g)));
}

inline double Ion5Moment2D_cState::s(void) const {
    return (R*gm1i*log(gm1*(E - HALF*dv.sqr()/d)/pow(d, g)));
}

/*************************************************************
 * Ion5Moment2D_cState::To -- Stagnation temperature.        *
 *************************************************************/
inline double Ion5Moment2D_cState::To(void) {
    return ((gm1*(E - HALF*dv.sqr()/d)/(d*R))*
	    (ONE+HALF*gm1*dv.sqr()/(d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))));
}

inline double Ion5Moment2D_cState::To(void) const {
    return ((gm1*(E - HALF*dv.sqr()/d)/(d*R))*
	    (ONE+HALF*gm1*dv.sqr()/(d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))));
}

/*************************************************************
 * Ion5Moment2D_cState::po -- Stagnation pressure.           *
 *************************************************************/
inline double Ion5Moment2D_cState::po(void) {
    return ((gm1*(E - HALF*dv.sqr()/d))*
	    pow(ONE+HALF*gm1*dv.sqr()/(d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d))), g*gm1i));
}

inline double Ion5Moment2D_cState::po(void) const {
    return ((gm1*(E - HALF*dv.sqr()/d))*
	    pow(ONE+HALF*gm1*dv.sqr()/(d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d))), g*gm1i));
}

/*************************************************************
 * Ion5Moment2D_cState::ao -- Stagnation sound speed.        *
 *************************************************************/
inline double Ion5Moment2D_cState::ao(void) {
    return (sqrt((g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))*
	         (ONE+HALF*gm1*dv.sqr()/(d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d))))));
}

inline double Ion5Moment2D_cState::ao(void) const {
    return (sqrt((g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))*
	         (ONE+HALF*gm1*dv.sqr()/(d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d))))));
}

/*************************************************************
 * Ion5Moment2D_cState::ho -- Stagnation enthalpy.           *
 *************************************************************/
inline double Ion5Moment2D_cState::ho(void) {
    return ((g*E/d - gm1*HALF*dv.sqr()/sqr(d))*
	    (ONE+HALF*gm1*dv.sqr()/(d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))));
}

inline double Ion5Moment2D_cState::ho(void) const {
    return ((g*E/d - gm1*HALF*dv.sqr()/sqr(d))*
	    (ONE+HALF*gm1*dv.sqr()/(d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))));
}

/*************************************************************
 * Ion5Moment2D_cState -- Binary arithmetic operators.       *
 *************************************************************/
inline Ion5Moment2D_cState operator +(const Ion5Moment2D_cState &U1, const Ion5Moment2D_cState &U2) {
  return (Ion5Moment2D_cState(U1.d+U2.d,U1.dv+U2.dv,U1.E+U2.E));
}

inline Ion5Moment2D_cState operator -(const Ion5Moment2D_cState &U1, const Ion5Moment2D_cState &U2) {
  return (Ion5Moment2D_cState(U1.d-U2.d,U1.dv-U2.dv,U1.E-U2.E));
}

// Inner product operator.
inline double operator *(const Ion5Moment2D_cState &U1, const Ion5Moment2D_cState &U2) {
   return (U1.d*U2.d+U1.dv*U2.dv+U1.E*U2.E);
}

inline Ion5Moment2D_cState operator *(const Ion5Moment2D_cState &U, const double &a) {
  return (Ion5Moment2D_cState(a*U.d,a*U.dv,a*U.E));
}

inline Ion5Moment2D_cState operator *(const double &a, const Ion5Moment2D_cState &U) {
  return (Ion5Moment2D_cState(a*U.d,a*U.dv,a*U.E));
}

inline Ion5Moment2D_cState operator /(const Ion5Moment2D_cState &U, const double &a) {
  return (Ion5Moment2D_cState(U.d/a,U.dv/a,U.E/a));
}

// My useful solution state product operator.
inline Ion5Moment2D_cState operator ^(const Ion5Moment2D_cState &U1, const Ion5Moment2D_cState &U2) {
   return (Ion5Moment2D_cState(U1.d*U2.d,U1.dv.x*U2.dv.x,U1.dv.y*U2.dv.y,U1.E*U2.E));
}

/*************************************************************
 * Ion5Moment2D_cState -- Unary arithmetic operators.        *
 *************************************************************/
inline Ion5Moment2D_cState operator +(const Ion5Moment2D_cState &U) {
  return (Ion5Moment2D_cState(U.d,U.dv,U.E));
}

inline Ion5Moment2D_cState operator -(const Ion5Moment2D_cState &U) {
  return (Ion5Moment2D_cState(-U.d,-U.dv,-U.E));
}

/*************************************************************
 * Ion5Moment2D_cState -- Shortcut arithmetic operators.     *
 *************************************************************/
inline Ion5Moment2D_cState &operator +=(Ion5Moment2D_cState &U1, const Ion5Moment2D_cState &U2) {
  U1.d += U2.d;
  U1.dv += U2.dv;
  U1.E += U2.E;
  return (U1);
}

inline Ion5Moment2D_cState &operator -=(Ion5Moment2D_cState &U1, const Ion5Moment2D_cState &U2) {
  U1.d -= U2.d;
  U1.dv -= U2.dv;
  U1.E -= U2.E;
  return (U1);
}

/*************************************************************
 * Ion5Moment2D_cState -- Relational operators.              *
 *************************************************************/
inline int operator ==(const Ion5Moment2D_cState &U1, const Ion5Moment2D_cState &U2) {
  return (U1.d == U2.d && U1.dv == U2.dv && U1.E == U2.E);
}

inline int operator !=(const Ion5Moment2D_cState &U1, const Ion5Moment2D_cState &U2) {
  return (U1.d != U2.d || U1.dv != U2.dv || U1.E != U2.E);
}

/*************************************************************
 * Ion5Moment2D_cState -- Input-output operators.            *
 *************************************************************/
inline ostream &operator << (ostream &out_file, const Ion5Moment2D_cState &U) {
  out_file.setf(ios::scientific);
  out_file << " " << U.d  << " " << U.dv.x << " " << U.dv.y << " " << U.E;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, Ion5Moment2D_cState &U) {
  in_file.setf(ios::skipws);
  in_file >> U.d >> U.dv.x >> U.dv.y >> U.E;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/******************************************************************
 * Ion5Moment2D_pState::Ion5Moment2D_pState -- Constructor.       *
 ******************************************************************/
inline Ion5Moment2D_pState::Ion5Moment2D_pState(const Ion5Moment2D_cState &U) {
  d = U.d; v = U.v(); p = U.p();
}

/*************************************************************
 * Ion5Moment2D_pState::U -- Conserved solution state.       *
 *************************************************************/
inline Ion5Moment2D_cState Ion5Moment2D_pState::U(void) {
  return (Ion5Moment2D_cState(d, dv(), E()));
}

inline Ion5Moment2D_cState Ion5Moment2D_pState::U(void) const {
  return (Ion5Moment2D_cState(d, dv(), E()));
}

inline Ion5Moment2D_cState Ion5Moment2D_pState::U(const Ion5Moment2D_pState &W) {
  return (Ion5Moment2D_cState(W.d, W.dv(), W.E()));
}

inline Ion5Moment2D_cState U(const Ion5Moment2D_pState &W) {
  return (Ion5Moment2D_cState(W.d, W.dv(), W.E()));
}

/*************************************************************
 * Ion5Moment2D_pState::F -- Solution flux (x-direction).    *
 *************************************************************/
inline Ion5Moment2D_cState Ion5Moment2D_pState::F(void) {
  return (Ion5Moment2D_cState(d*v.x, d*sqr(v.x) + p, d*v.x*v.y, v.x*H()));
}

inline Ion5Moment2D_cState Ion5Moment2D_pState::F(void) const {
  return (Ion5Moment2D_cState(d*v.x, d*sqr(v.x) + p, d*v.x*v.y, v.x*H()));
}

inline Ion5Moment2D_cState Ion5Moment2D_pState::F(const Ion5Moment2D_pState &W) {
  return (Ion5Moment2D_cState(W.d*W.v.x, W.d*sqr(W.v.x) + W.p,
                              W.d*W.v.x*W.v.y, W.v.x*W.H()));
}

inline Ion5Moment2D_cState F(const Ion5Moment2D_pState &W) {
  return (Ion5Moment2D_cState(W.d*W.v.x, W.d*sqr(W.v.x) + W.p,
                              W.d*W.v.x*W.v.y, W.v.x*W.H()));
}

/*************************************************************
 * Ion5Moment2D_pState::Fx -- Solution flux (x-direction).   *
 *************************************************************/
inline Ion5Moment2D_cState Ion5Moment2D_pState::Fx(void) {
  return (Ion5Moment2D_cState(d*v.x, d*sqr(v.x) + p, d*v.x*v.y, v.x*H()));
}

inline Ion5Moment2D_cState Ion5Moment2D_pState::Fx(void) const {
  return (Ion5Moment2D_cState(d*v.x, d*sqr(v.x) + p, d*v.x*v.y, v.x*H()));
}

inline Ion5Moment2D_cState Ion5Moment2D_pState::Fx(const Ion5Moment2D_pState &W) {
  return (Ion5Moment2D_cState(W.d*W.v.x, W.d*sqr(W.v.x) + W.p,
                              W.d*W.v.x*W.v.y, W.v.x*W.H()));
}

inline Ion5Moment2D_cState Fx(const Ion5Moment2D_pState &W) {
  return (Ion5Moment2D_cState(W.d*W.v.x, W.d*sqr(W.v.x) + W.p,
                              W.d*W.v.x*W.v.y, W.v.x*W.H()));
}

/*************************************************************
 * Ion5Moment2D_pState::Fy -- Solution flux (y-direction).   *
 *************************************************************/
inline Ion5Moment2D_cState Ion5Moment2D_pState::Fy(void) {
  return (Ion5Moment2D_cState(d*v.y, d*v.x*v.y, d*sqr(v.y) + p, v.y*H()));
}

inline Ion5Moment2D_cState Ion5Moment2D_pState::Fy(void) const {
  return (Ion5Moment2D_cState(d*v.y, d*v.x*v.y, d*sqr(v.y) + p, v.y*H()));
}

inline Ion5Moment2D_cState Ion5Moment2D_pState::Fy(const Ion5Moment2D_pState &W) {
  return (Ion5Moment2D_cState(W.d*W.v.y, W.d*W.v.x*W.v.y,
                              W.d*sqr(W.v.y) + W.p, W.v.y*W.H()));
}

inline Ion5Moment2D_cState Fy(const Ion5Moment2D_pState &W) {
  return (Ion5Moment2D_cState(W.d*W.v.y, W.d*W.v.x*W.v.y,
                         W.d*sqr(W.v.y) + W.p, W.v.y*W.H()));
}

/*************************************************************
 * Ion5Moment2D_pState::Fn -- Solution flux (n-direction).   *
 *************************************************************/
inline Ion5Moment2D_cState Ion5Moment2D_pState::Fn(void) {
  return (Ion5Moment2D_cState(d*v.x, d*sqr(v.x) + p, d*v.x*v.y, v.x*H()));
}

inline Ion5Moment2D_cState Ion5Moment2D_pState::Fn(void) const {
  return (Ion5Moment2D_cState(d*v.x, d*sqr(v.x) + p, d*v.x*v.y, v.x*H()));
}

inline Ion5Moment2D_cState Ion5Moment2D_pState::Fn(const Ion5Moment2D_pState &W) {
  return (Ion5Moment2D_cState(W.d*W.v.x, W.d*sqr(W.v.x) + W.p,
                         W.d*W.v.x*W.v.y, W.v.x*W.H()));
}

inline Ion5Moment2D_cState Fn(const Ion5Moment2D_pState &W) {
  return (Ion5Moment2D_cState(W.d*W.v.x, W.d*sqr(W.v.x) + W.p,
                              W.d*W.v.x*W.v.y, W.v.x*W.H()));
}

/*****************************************************************
 * Ion5Moment2D_pState::lambda -- Eigenvalue(s) (x-direction).   *
 *****************************************************************/
inline Ion5Moment2D_pState Ion5Moment2D_pState::lambda(void) {
  double c = a();
  return (Ion5Moment2D_pState(v.x - c, v.x, v.x, v.x + c));
}

inline Ion5Moment2D_pState Ion5Moment2D_pState::lambda(void) const {
  double c = a();
  return (Ion5Moment2D_pState(v.x - c, v.x, v.x, v.x + c));
}

inline Ion5Moment2D_pState Ion5Moment2D_pState::lambda(const Ion5Moment2D_pState &W) {
  double c = W.a();
  return (Ion5Moment2D_pState(W.v.x - c, W.v.x, W.v.x, W.v.x + c));
}

inline Ion5Moment2D_pState lambda(const Ion5Moment2D_pState &W) {
  double c = W.a();
  return (Ion5Moment2D_pState(W.v.x - c, W.v.x, W.v.x, W.v.x + c));
}

inline double Ion5Moment2D_pState::lambda(int index) {
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      return (v.x-a());
    case 2 :
      return (v.x);
    case 3 :
      return (v.x);
    case 4 :
      return (v.x+a());
    default:
      return (v.x);
  };
}

inline double Ion5Moment2D_pState::lambda(int index) const {
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      return (v.x-a());
    case 2 :
      return (v.x);
    case 3 :
      return (v.x);
    case 4 :
      return (v.x+a());
    default:
      return (v.x);
  };
}

inline double lambda(const Ion5Moment2D_pState &W, int index) {
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      return (W.v.x-W.a());
    case 2 :
      return (W.v.x);
    case 3 :
      return (W.v.x);
    case 4 :
      return (W.v.x+W.a());
    default:
      return (W.v.x);
  };
}

/*****************************************************************
 * Ion5Moment2D_pState::lambda_x -- Eigenvalue(s) (x-direction). *
 *****************************************************************/
inline Ion5Moment2D_pState Ion5Moment2D_pState::lambda_x(void) {
  double c = a();
  return (Ion5Moment2D_pState(v.x - c, v.x, v.x, v.x + c));
}

inline Ion5Moment2D_pState Ion5Moment2D_pState::lambda_x(void) const {
  double c = a();
  return (Ion5Moment2D_pState(v.x - c, v.x, v.x, v.x + c));
}

inline Ion5Moment2D_pState Ion5Moment2D_pState::lambda_x(const Ion5Moment2D_pState &W) {
  double c = W.a();
  return (Ion5Moment2D_pState(W.v.x - c, W.v.x, W.v.x, W.v.x + c));
}

inline Ion5Moment2D_pState lambda_x(const Ion5Moment2D_pState &W) {
  double c = W.a();
  return (Ion5Moment2D_pState(W.v.x - c, W.v.x, W.v.x, W.v.x + c));
}

inline double Ion5Moment2D_pState::lambda_x(int index) {
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      return (v.x-a());
    case 2 :
      return (v.x);
    case 3 :
      return (v.x);
    case 4 :
      return (v.x+a());
    default:
      return (v.x);
  };
}

inline double Ion5Moment2D_pState::lambda_x(int index) const {
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      return (v.x-a());
    case 2 :
      return (v.x);
    case 3 :
      return (v.x);
    case 4 :
      return (v.x+a());
    default:
      return (v.x);
  };
}

inline double lambda_x(const Ion5Moment2D_pState &W, int index) {
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      return (W.v.x-W.a());
    case 2 :
      return (W.v.x);
    case 3 :
      return (W.v.x);
    case 4 :
      return (W.v.x+W.a());
    default:
      return (W.v.x);
  };
}

/*****************************************************************
 * Ion5Moment2D_pState::lambda_y -- Eigenvalue(s) (y-direction). *
 *****************************************************************/
inline Ion5Moment2D_pState Ion5Moment2D_pState::lambda_y(void) {
  double c = a();
  return (Ion5Moment2D_pState(v.y - c, v.y, v.y, v.y + c));
}

inline Ion5Moment2D_pState Ion5Moment2D_pState::lambda_y(void) const {
  double c = a();
  return (Ion5Moment2D_pState(v.y - c, v.y, v.y, v.y + c));
}

inline Ion5Moment2D_pState Ion5Moment2D_pState::lambda_y(const Ion5Moment2D_pState &W) {
  double c = W.a();
  return (Ion5Moment2D_pState(W.v.y - c, W.v.y, W.v.y, W.v.y + c));
}

inline Ion5Moment2D_pState lambda_y(const Ion5Moment2D_pState &W) {
  double c = W.a();
  return (Ion5Moment2D_pState(W.v.y - c, W.v.y, W.v.y, W.v.y + c));
}

inline double Ion5Moment2D_pState::lambda_y(int index) {
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      return (v.y-a());
    case 2 :
      return (v.y);
    case 3 :
      return (v.y);
    case 4 :
      return (v.y+a());
    default:
      return (v.y);
  };
}

inline double Ion5Moment2D_pState::lambda_y(int index) const {
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      return (v.y-a());
    case 2 :
      return (v.y);
    case 3 :
      return (v.y);
    case 4 :
      return (v.y+a());
    default:
      return (v.y);
  };
}

inline double lambda_y(const Ion5Moment2D_pState &W, int index) {
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      return (W.v.y-W.a());
    case 2 :
      return (W.v.y);
    case 3 :
      return (W.v.y);
    case 4 :
      return (W.v.y+W.a());
    default:
      return (W.v.y);
  };
}

/*************************************************************
 * Ion5Moment2D_pState::rp -- Primitive right eigenvector    *
 *                            (x-direction).                 *
 *************************************************************/
inline Ion5Moment2D_pState Ion5Moment2D_pState::rp(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = a();
      return (Ion5Moment2D_pState(ONE, -c/d, ZERO, sqr(c)));
    case 2 :
      return (Ion5Moment2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Ion5Moment2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = a();
      return (Ion5Moment2D_pState(ONE, c/d, ZERO, sqr(c)));
    default:
      c = a();
      return (Ion5Moment2D_pState(ONE, -c/d, ZERO, sqr(c)));
  };
}

inline Ion5Moment2D_pState Ion5Moment2D_pState::rp(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = a();
      return (Ion5Moment2D_pState(ONE, -c/d, ZERO, sqr(c)));
    case 2 :
      return (Ion5Moment2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Ion5Moment2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = a();
      return (Ion5Moment2D_pState(ONE, c/d, ZERO, sqr(c)));
    default:
      c = a();
      return (Ion5Moment2D_pState(ONE, -c/d, ZERO, sqr(c)));
  };
}

inline Ion5Moment2D_pState rp(const Ion5Moment2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Ion5Moment2D_pState(ONE, -c/W.d, ZERO, sqr(c)));
    case 2 :
      return (Ion5Moment2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Ion5Moment2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = W.a();
      return (Ion5Moment2D_pState(ONE, c/W.d, ZERO, sqr(c)));
    default:
      c = W.a();
      return (Ion5Moment2D_pState(ONE, -c/W.d, ZERO, sqr(c)));
  };
}

/*************************************************************
 * Ion5Moment2D_pState::rp_x -- Primitive right eigenvector  *
 *                              (x-direction).               *
 *************************************************************/
inline Ion5Moment2D_pState Ion5Moment2D_pState::rp_x(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = a();
      return (Ion5Moment2D_pState(ONE, -c/d, ZERO, sqr(c)));
    case 2 :
      return (Ion5Moment2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Ion5Moment2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = a();
      return (Ion5Moment2D_pState(ONE, c/d, ZERO, sqr(c)));
    default:
      c = a();
      return (Ion5Moment2D_pState(ONE, -c/d, ZERO, sqr(c)));
  };
}

inline Ion5Moment2D_pState Ion5Moment2D_pState::rp_x(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = a();
      return (Ion5Moment2D_pState(ONE, -c/d, ZERO, sqr(c)));
    case 2 :
      return (Ion5Moment2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Ion5Moment2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = a();
      return (Ion5Moment2D_pState(ONE, c/d, ZERO, sqr(c)));
    default:
      c = a();
      return (Ion5Moment2D_pState(ONE, -c/d, ZERO, sqr(c)));
  };
}

inline Ion5Moment2D_pState rp_x(const Ion5Moment2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Ion5Moment2D_pState(ONE, -c/W.d, ZERO, sqr(c)));
    case 2 :
      return (Ion5Moment2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Ion5Moment2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = W.a();
      return (Ion5Moment2D_pState(ONE, c/W.d, ZERO, sqr(c)));
    default:
      c = W.a();
      return (Ion5Moment2D_pState(ONE, -c/W.d, ZERO, sqr(c)));
  };
}

/*************************************************************
 * Ion5Moment2D_pState::rp_y -- Primitive right eigenvector  *
 *                              (y-direction).               *
 *************************************************************/
inline Ion5Moment2D_pState Ion5Moment2D_pState::rp_y(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = a();
      return (Ion5Moment2D_pState(ONE, ZERO, -c/d, sqr(c)));
    case 2 :
      return (Ion5Moment2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Ion5Moment2D_pState(ZERO, ONE, ZERO, ZERO));
    case 4 :
      c = a();
      return (Ion5Moment2D_pState(ONE, ZERO, c/d, sqr(c)));
    default:
      c = a();
      return (Ion5Moment2D_pState(ONE, ZERO, -c/d, sqr(c)));
  };
}

inline Ion5Moment2D_pState Ion5Moment2D_pState::rp_y(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = a();
      return (Ion5Moment2D_pState(ONE, ZERO, -c/d, sqr(c)));
    case 2 :
      return (Ion5Moment2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Ion5Moment2D_pState(ZERO, ONE, ZERO, ZERO));
    case 4 :
      c = a();
      return (Ion5Moment2D_pState(ONE, ZERO, c/d, sqr(c)));
    default:
      c = a();
      return (Ion5Moment2D_pState(ONE, ZERO, -c/d, sqr(c)));
  };
}

inline Ion5Moment2D_pState rp_y(const Ion5Moment2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Ion5Moment2D_pState(ONE, ZERO, -c/W.d, sqr(c)));
    case 2 :
      return (Ion5Moment2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Ion5Moment2D_pState(ZERO, ONE, ZERO, ZERO));
    case 4 :
      c = W.a();
      return (Ion5Moment2D_pState(ONE, ZERO, c/W.d, sqr(c)));
    default:
      c = W.a();
      return (Ion5Moment2D_pState(ONE, ZERO, -c/W.d, sqr(c)));
  };
}

/*************************************************************
 * Ion5Moment2D_pState::rc -- Conserved right eigenvector    *
 *                            (x-direction).                 *
 *************************************************************/
inline Ion5Moment2D_cState Ion5Moment2D_pState::rc(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = a();
      return (Ion5Moment2D_cState(ONE, v.x-c, v.y, h()-v.x*c));
    case 2 :
      return (Ion5Moment2D_cState(ONE, v.x, v.y, HALF*v.sqr()));
    case 3 :
      return (Ion5Moment2D_cState(ZERO, ZERO, d, d*v.y));
    case 4 :
      c = a();
      return (Ion5Moment2D_cState(ONE, v.x+c, v.y, h()+v.x*c));
    default:
      c = a();
      return (Ion5Moment2D_cState(ONE, v.x-c, v.y, h()-v.x*c));
  };
}

inline Ion5Moment2D_cState Ion5Moment2D_pState::rc(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = a();
      return (Ion5Moment2D_cState(ONE, v.x-c, v.y, h()-v.x*c));
    case 2 :
      return (Ion5Moment2D_cState(ONE, v.x, v.y, HALF*v.sqr()));
    case 3 :
      return (Ion5Moment2D_cState(ZERO, ZERO, d, d*v.y));
    case 4 :
      c = a();
      return (Ion5Moment2D_cState(ONE, v.x+c, v.y, h()+v.x*c));
    default:
      c = a();
      return (Ion5Moment2D_cState(ONE, v.x-c, v.y, h()-v.x*c));
  };
}

inline Ion5Moment2D_cState rc(const Ion5Moment2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Ion5Moment2D_cState(ONE, W.v.x-c, W.v.y, W.h()-W.v.x*c));
    case 2 :
      return (Ion5Moment2D_cState(ONE, W.v.x, W.v.y, HALF*W.v.sqr()));
    case 3 :
      return (Ion5Moment2D_cState(ZERO, ZERO, W.d, W.d*W.v.y));
    case 4 :
      c = W.a();
      return (Ion5Moment2D_cState(ONE, W.v.x+c, W.v.y, W.h()+W.v.x*c));
    default:
      c = W.a();
      return (Ion5Moment2D_cState(ONE, W.v.x-c, W.v.y, W.h()-W.v.x*c));
    };
}

/*************************************************************
 * Ion5Moment2D_pState::rc_x -- Conserved right eigenvector  *
 *                              (x-direction).               *
 *************************************************************/
inline Ion5Moment2D_cState Ion5Moment2D_pState::rc_x(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = a();
      return (Ion5Moment2D_cState(ONE, v.x-c, v.y, h()-v.x*c));
    case 2 :
      return (Ion5Moment2D_cState(ONE, v.x, v.y, HALF*v.sqr()));
    case 3 :
      return (Ion5Moment2D_cState(ZERO, ZERO, d, d*v.y));
    case 4 :
      c = a();
      return (Ion5Moment2D_cState(ONE, v.x+c, v.y, h()+v.x*c));
    default:
      c = a();
      return (Ion5Moment2D_cState(ONE, v.x-c, v.y, h()-v.x*c));
  };
}

inline Ion5Moment2D_cState Ion5Moment2D_pState::rc_x(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = a();
      return (Ion5Moment2D_cState(ONE, v.x-c, v.y, h()-v.x*c));
    case 2 :
      return (Ion5Moment2D_cState(ONE, v.x, v.y, HALF*v.sqr()));
    case 3 :
      return (Ion5Moment2D_cState(ZERO, ZERO, d, d*v.y));
    case 4 :
      c = a();
      return (Ion5Moment2D_cState(ONE, v.x+c, v.y, h()+v.x*c));
    default:
      c = a();
      return (Ion5Moment2D_cState(ONE, v.x-c, v.y, h()-v.x*c));
  };
}

inline Ion5Moment2D_cState rc_x(const Ion5Moment2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Ion5Moment2D_cState(ONE, W.v.x-c, W.v.y, W.h()-W.v.x*c));
    case 2 :
      return (Ion5Moment2D_cState(ONE, W.v.x, W.v.y, HALF*W.v.sqr()));
    case 3 :
      return (Ion5Moment2D_cState(ZERO, ZERO, W.d, W.d*W.v.y));
    case 4 :
      c = W.a();
      return (Ion5Moment2D_cState(ONE, W.v.x+c, W.v.y, W.h()+W.v.x*c));
    default:
      c = W.a();
      return (Ion5Moment2D_cState(ONE, W.v.x-c, W.v.y, W.h()-W.v.x*c));
    };
}

/*************************************************************
 * Ion5Moment2D_pState::rc_y -- Conserved right eigenvector  *
 *                              (y-direction).               *
 *************************************************************/
inline Ion5Moment2D_cState Ion5Moment2D_pState::rc_y(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = a();
      return (Ion5Moment2D_cState(ONE, v.x, v.y-c, h()-v.y*c));
    case 2 :
      return (Ion5Moment2D_cState(ONE, v.x, v.y, HALF*v.sqr()));
    case 3 :
      return (Ion5Moment2D_cState(ZERO, d, ZERO, d*v.x));
    case 4 :
      c = a();
      return (Ion5Moment2D_cState(ONE, v.x, v.y+c, h()+v.y*c));
    default:
      c = a();
      return (Ion5Moment2D_cState(ONE, v.x, v.y-c, h()-v.y*c));
  };
}

inline Ion5Moment2D_cState Ion5Moment2D_pState::rc_y(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = a();
      return (Ion5Moment2D_cState(ONE, v.x, v.y-c, h()-v.y*c));
    case 2 :
      return (Ion5Moment2D_cState(ONE, v.x, v.y, HALF*v.sqr()));
    case 3 :
      return (Ion5Moment2D_cState(ZERO, d, ZERO, d*v.x));
    case 4 :
      c = a();
      return (Ion5Moment2D_cState(ONE, v.x, v.y+c, h()+v.y*c));
    default:
      c = a();
      return (Ion5Moment2D_cState(ONE, v.x, v.y-c, h()-v.y*c));
  };
}

inline Ion5Moment2D_cState rc_y(const Ion5Moment2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Ion5Moment2D_cState(ONE, W.v.x, W.v.y-c, W.h()-W.v.y*c));
    case 2 :
      return (Ion5Moment2D_cState(ONE, W.v.x, W.v.y, HALF*W.v.sqr()));
    case 3 :
      return (Ion5Moment2D_cState(ZERO, W.d, ZERO, W.d*W.v.x));
    case 4 :
      c = W.a();
      return (Ion5Moment2D_cState(ONE, W.v.x, W.v.y+c, W.h()+W.v.y*c));
    default:
      c = W.a();
      return (Ion5Moment2D_cState(ONE, W.v.x, W.v.y-c, W.h()-W.v.y*c));
    };
}

/*************************************************************
 * Ion5Moment2D_pState::lp -- Primitive left eigenvector     *
 *                            (x-direction).                 *
 *************************************************************/
inline Ion5Moment2D_pState Ion5Moment2D_pState::lp(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = a();
      return (Ion5Moment2D_pState(ZERO, -HALF*d/c, ZERO, HALF/sqr(c)));
    case 2 :
      return (Ion5Moment2D_pState(ONE, ZERO, ZERO, -ONE/a2()));
    case 3 :
      return (Ion5Moment2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = a();
      return (Ion5Moment2D_pState(ZERO, HALF*d/c, ZERO, HALF/sqr(c)));
    default:
      c = a();
      return (Ion5Moment2D_pState(ZERO, -HALF*d/c, ZERO, HALF/sqr(c)));
  };
}

inline Ion5Moment2D_pState Ion5Moment2D_pState::lp(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = a();
      return (Ion5Moment2D_pState(ZERO, -HALF*d/c, ZERO, HALF/sqr(c)));
    case 2 :
      return (Ion5Moment2D_pState(ONE, ZERO, ZERO, -ONE/a2()));
    case 3 :
      return (Ion5Moment2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = a();
      return (Ion5Moment2D_pState(ZERO, HALF*d/c, ZERO, HALF/sqr(c)));
    default:
      c = a();
      return (Ion5Moment2D_pState(ZERO, -HALF*d/c, ZERO, HALF/sqr(c)));
  };
}

inline Ion5Moment2D_pState lp(const Ion5Moment2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Ion5Moment2D_pState(ZERO, -HALF*W.d/c, ZERO, HALF/sqr(c)));
    case 2 :
      return (Ion5Moment2D_pState(ONE, ZERO, ZERO, -ONE/W.a2()));
    case 3 :
      return (Ion5Moment2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = W.a();
      return (Ion5Moment2D_pState(ZERO, HALF*W.d/c, ZERO, HALF/sqr(c)));
    default:
      c = W.a();
      return (Ion5Moment2D_pState(ZERO, -HALF*W.d/c, ZERO, HALF/sqr(c)));
  };
}

/*************************************************************
 * Ion5Moment2D_pState::lp_x -- Primitive left eigenvector   *
 *                              (x-direction).               *
 *************************************************************/
inline Ion5Moment2D_pState Ion5Moment2D_pState::lp_x(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = a();
      return (Ion5Moment2D_pState(ZERO, -HALF*d/c, ZERO, HALF/sqr(c)));
    case 2 :
      return (Ion5Moment2D_pState(ONE, ZERO, ZERO, -ONE/a2()));
    case 3 :
      return (Ion5Moment2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = a();
      return (Ion5Moment2D_pState(ZERO, HALF*d/c, ZERO, HALF/sqr(c)));
    default:
      c = a();
      return (Ion5Moment2D_pState(ZERO, -HALF*d/c, ZERO, HALF/sqr(c)));
  };
}

inline Ion5Moment2D_pState Ion5Moment2D_pState::lp_x(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = a();
      return (Ion5Moment2D_pState(ZERO, -HALF*d/c, ZERO, HALF/sqr(c)));
    case 2 :
      return (Ion5Moment2D_pState(ONE, ZERO, ZERO, -ONE/a2()));
    case 3 :
      return (Ion5Moment2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = a();
      return (Ion5Moment2D_pState(ZERO, HALF*d/c, ZERO, HALF/sqr(c)));
    default:
      c = a();
      return (Ion5Moment2D_pState(ZERO, -HALF*d/c, ZERO, HALF/sqr(c)));
  };
}

inline Ion5Moment2D_pState lp_x(const Ion5Moment2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Ion5Moment2D_pState(ZERO, -HALF*W.d/c, ZERO, HALF/sqr(c)));
    case 2 :
      return (Ion5Moment2D_pState(ONE, ZERO, ZERO, -ONE/W.a2()));
    case 3 :
      return (Ion5Moment2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = W.a();
      return (Ion5Moment2D_pState(ZERO, HALF*W.d/c, ZERO, HALF/sqr(c)));
    default:
      c = W.a();
      return (Ion5Moment2D_pState(ZERO, -HALF*W.d/c, ZERO, HALF/sqr(c)));
  };
}

/*************************************************************
 * Ion5Moment2D_pState::lp_y -- Primitive left eigenvector   *
 *                              (y-direction).               *
 *************************************************************/
inline Ion5Moment2D_pState Ion5Moment2D_pState::lp_y(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = a();
      return (Ion5Moment2D_pState(ZERO, ZERO, -HALF*d/c, HALF/sqr(c)));
    case 2 :
      return (Ion5Moment2D_pState(ONE, ZERO, ZERO, -ONE/a2()));
    case 3 :
      return (Ion5Moment2D_pState(ZERO, ONE, ZERO, ZERO));
    case 4 :
      c = a();
      return (Ion5Moment2D_pState(ZERO, ZERO, HALF*d/c, HALF/sqr(c)));
    default:
      c = a();
      return (Ion5Moment2D_pState(ZERO, ZERO, -HALF*d/c, HALF/sqr(c)));
  };
}

inline Ion5Moment2D_pState Ion5Moment2D_pState::lp_y(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = a();
      return (Ion5Moment2D_pState(ZERO, ZERO, -HALF*d/c, HALF/sqr(c)));
    case 2 :
      return (Ion5Moment2D_pState(ONE, ZERO, ZERO, -ONE/a2()));
    case 3 :
      return (Ion5Moment2D_pState(ZERO, ONE, ZERO, ZERO));
    case 4 :
      c = a();
      return (Ion5Moment2D_pState(ZERO, ZERO, HALF*d/c, HALF/sqr(c)));
    default:
      c = a();
      return (Ion5Moment2D_pState(ZERO, ZERO, -HALF*d/c, HALF/sqr(c)));
  };
}

inline Ion5Moment2D_pState lp_y(const Ion5Moment2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_ION5MOMENT2D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Ion5Moment2D_pState(ZERO, ZERO, -HALF*W.d/c, HALF/sqr(c)));
    case 2 :
      return (Ion5Moment2D_pState(ONE, ZERO, ZERO, -ONE/W.a2()));
    case 3 :
      return (Ion5Moment2D_pState(ZERO, ONE, ZERO, ZERO));
    case 4 :
      c = W.a();
      return (Ion5Moment2D_pState(ZERO, ZERO, HALF*W.d/c, HALF/sqr(c)));
    default:
      c = W.a();
      return (Ion5Moment2D_pState(ZERO, ZERO, -HALF*W.d/c, HALF/sqr(c)));
  };
}

/***************************************************************
 * Ion5Moment2D_pState::Sa -- Axisymmetric flow source terms.  *
 ***************************************************************/
inline Ion5Moment2D_cState Ion5Moment2D_pState::Sa(const Vector2D &X) {
  return (Ion5Moment2D_cState(-d*v.y/X.y, -d*v.x*v.y/X.y, 
                              -d*sqr(v.y)/X.y, -v.y*H()/X.y));
}

inline Ion5Moment2D_cState Ion5Moment2D_pState::Sa(const Vector2D &X) const {
  return (Ion5Moment2D_cState(-d*v.y/X.y, -d*v.x*v.y/X.y, 
                              -d*sqr(v.y)/X.y, -v.y*H()/X.y));
}

inline Ion5Moment2D_cState Sa(const Ion5Moment2D_pState &W, const Vector2D &X) {
  return (Ion5Moment2D_cState(-W.d*W.v.y/X.y, -W.d*W.v.x*W.v.y/X.y, 
                              -W.d*sqr(W.v.y)/X.y, -W.v.y*W.H()/X.y));
}

inline void Ion5Moment2D_pState::dSadU(DenseMatrix &dSadU, const Vector2D &X) {
  dSadU(0,2) -= ONE/X.y;
  dSadU(1,0) += v.x*v.y/X.y;
  dSadU(1,1) -= v.y/X.y;
  dSadU(1,2) -= v.x/X.y;
  dSadU(2,0) += v.y*v.y/X.y;
  dSadU(2,2) -= TWO*v.y/X.y;
  dSadU(3,0) += HALF*gm1i*v.y*(TWO*g*p/d+(TWO-g)*gm1*(v.x*v.x+v.y*v.y))/X.y;
  dSadU(3,1) += gm1*v.x*v.y/X.y;
  dSadU(3,2) -= HALF*gm1i*(TWO*g*p/d+gm1*v.x*v.x-(TWO*g-THREE)*gm1*v.y*v.y)/X.y;
  dSadU(3,3) -= g*v.y/X.y;
}

inline void Ion5Moment2D_pState::dSadU(DenseMatrix &dSadU, const Vector2D &X) const {
  dSadU(0,2) -= ONE/X.y;
  dSadU(1,0) += v.x*v.y/X.y;
  dSadU(1,1) -= v.y/X.y;
  dSadU(1,2) -= v.x/X.y;
  dSadU(2,0) += v.y*v.y/X.y;
  dSadU(2,2) -= TWO*v.y/X.y;
  dSadU(3,0) += HALF*gm1i*v.y*(TWO*g*p/d+(TWO-g)*gm1*(v.x*v.x+v.y*v.y))/X.y;
  dSadU(3,1) += gm1*v.x*v.y/X.y;
  dSadU(3,2) -= HALF*gm1i*(TWO*g*p/d+gm1*v.x*v.x-(TWO*g-THREE)*gm1*v.y*v.y)/X.y;
  dSadU(3,3) -= g*v.y/X.y;
}

inline void dSadU(DenseMatrix &dSadU, const Ion5Moment2D_pState &W, const Vector2D &X) {
  dSadU(0,2) -= ONE/X.y;
  dSadU(1,0) += W.v.x*W.v.y/X.y;
  dSadU(1,1) -= W.v.y/X.y;
  dSadU(1,2) -= W.v.x/X.y;
  dSadU(2,0) += W.v.y*W.v.y/X.y;
  dSadU(2,2) -= TWO*W.v.y/X.y;
  dSadU(3,0) += HALF*W.gm1i*W.v.y*(TWO*W.g*W.p/W.d+(TWO-W.g)*W.gm1*(W.v.x*W.v.x+W.v.y*W.v.y))/X.y;
  dSadU(3,1) += W.gm1*W.v.x*W.v.y/X.y;
  dSadU(3,2) -= HALF*W.gm1i*(TWO*W.g*W.p/W.d+W.gm1*W.v.x*W.v.x-(TWO*W.g-THREE)*W.gm1*W.v.y*W.v.y)/X.y;
  dSadU(3,3) -= W.g*W.v.y/X.y;
}

/***************************************************************
 * Ion5Moment2D_pState::Se -- Electric field source terms.     *
 ***************************************************************/
inline Ion5Moment2D_cState Ion5Moment2D_pState::Se(const Vector2D &E) {
  double cm; cm = d*(q/m);
  return (Ion5Moment2D_cState(ZERO, cm*E.x, cm*E.y, cm*(v*E)));
}

inline Ion5Moment2D_cState Ion5Moment2D_pState::Se(const Vector2D &E) const {
  double cm; cm = d*(q/m);
  return (Ion5Moment2D_cState(ZERO, cm*E.x, cm*E.y, cm*(v*E)));
}

inline Ion5Moment2D_cState Se(const Ion5Moment2D_pState &W, const Vector2D &E) {
  double cm; cm = W.d*(W.q/W.m);
  return (Ion5Moment2D_cState(ZERO, cm*E.x, cm*E.y, cm*(W.v*E)));
}

inline void Ion5Moment2D_pState::dSedU(DenseMatrix &dSedU, const Vector2D &E) {
  dSedU(1,0) += (q/m)*E.x;
  dSedU(2,0) += (q/m)*E.y;
  dSedU(3,1) += (q/m)*E.x;
  dSedU(3,2) += (q/m)*E.y;
}

inline void Ion5Moment2D_pState::dSedU(DenseMatrix &dSedU, const Vector2D &E) const {
  dSedU(1,0) += (q/m)*E.x;
  dSedU(2,0) += (q/m)*E.y;
  dSedU(3,1) += (q/m)*E.x;
  dSedU(3,2) += (q/m)*E.y;
}

inline void dSedU(DenseMatrix &dSedU, const Ion5Moment2D_pState &W, const Vector2D &E) {
  dSedU(1,0) += (W.q/W.m)*E.x;
  dSedU(2,0) += (W.q/W.m)*E.y;
  dSedU(3,1) += (W.q/W.m)*E.x;
  dSedU(3,2) += (W.q/W.m)*E.y;
}

/***************************************************************
 * Ion5Moment2D_pState::nu -- Ion-neutral collision frequency. *
 ***************************************************************/
inline double Ion5Moment2D_pState::nu(const Euler2D_pState &Wneut) {
  double nu_in, mneut, sigma;
  mneut = mn();
  switch(ion) {
    case ION_H :
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
    case ION_HE :
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
    case ION_O :
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
    case ION_LOW_MASS :
      sigma = 105*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*T())/(m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    case ION_MED_MASS :
      sigma = 150*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*T())/(m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    case ION_HI_MASS :
      sigma = 280*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*T())/(m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    case ION_HUGE_MASS :
      sigma = 2560*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*T())/(m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    case ION_HUGE_MASS2 :
      sigma = 433*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*T())/(m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    default:
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
  };
  return (nu_in);
}

inline double Ion5Moment2D_pState::nu(const Euler2D_pState &Wneut) const {
  double nu_in, mneut, sigma;
  mneut = mn();
  switch(ion) {
    case ION_H :
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
    case ION_HE :
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
    case ION_O :
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
    case ION_LOW_MASS :
      sigma = 105*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*T())/(m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    case ION_MED_MASS :
      sigma = 150*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*T())/(m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    case ION_HI_MASS :
      sigma = 280*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*T())/(m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    case ION_HUGE_MASS :
      sigma = 2560*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*T())/(m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    case ION_HUGE_MASS2 :
      sigma = 433*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*T())/(m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    default:
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
  };
  return (nu_in);
}

inline double nu(const Ion5Moment2D_pState &W, const Euler2D_pState &Wneut) {
  double nu_in, mneut, sigma;
  mneut = W.mn();
  switch(W.ion) {
    case ION_H :
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
    case ION_HE :
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
    case ION_O :
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
    case ION_LOW_MASS :
      sigma = 105*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*W.T())/(W.m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    case ION_MED_MASS :
      sigma = 150*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*W.T())/(W.m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    case ION_HI_MASS :
      sigma = 280*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*W.T())/(W.m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    case ION_HUGE_MASS :
      sigma = 2560*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*W.T())/(W.m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    case ION_HUGE_MASS2 :
      sigma = 433*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*W.T())/(W.m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    default:
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
  };
  return (nu_in);
}

/***************************************************************
 * Ion5Moment2D_pState::mn -- Neutral gas particle mass.       *
 ***************************************************************/
inline double Ion5Moment2D_pState::mn(void) {
  switch(gas) {
    case GAS_AIR :
      return (MOLE_WT_AIR/(AVOGADRO*THOUSAND));
    case GAS_A :
      return (MOLE_WT_A/(AVOGADRO*THOUSAND));
    case GAS_CO :
      return (MOLE_WT_CO/(AVOGADRO*THOUSAND));
    case GAS_CO2 :
      return (MOLE_WT_CO2/(AVOGADRO*THOUSAND));
    case GAS_CH4 :
      return (MOLE_WT_CH4/(AVOGADRO*THOUSAND));
    case GAS_H :
      return (MOLE_WT_H/(AVOGADRO*THOUSAND));
    case GAS_H2 :
      return (MOLE_WT_H2/(AVOGADRO*THOUSAND));
    case GAS_HE :
      return (MOLE_WT_HE/(AVOGADRO*THOUSAND));
    case GAS_H2O :
      return (MOLE_WT_H2O/(AVOGADRO*THOUSAND));
    case GAS_N2 :
      return (MOLE_WT_N2/(AVOGADRO*THOUSAND));
    case GAS_O :
      return (MOLE_WT_O/(AVOGADRO*THOUSAND));
    case GAS_O2 :
      return (MOLE_WT_O2/(AVOGADRO*THOUSAND));
    case GAS_e :
      return (MOLE_WT_e/(AVOGADRO*THOUSAND));
    case GAS_APHTPB :
      return (MOLE_WT_APHTPB/(AVOGADRO*THOUSAND));
    default:
      return (MOLE_WT_AIR/(AVOGADRO*THOUSAND));
  };
}

inline double Ion5Moment2D_pState::mn(void) const {
  switch(gas) {
    case GAS_AIR :
      return (MOLE_WT_AIR/(AVOGADRO*THOUSAND));
    case GAS_A :
      return (MOLE_WT_A/(AVOGADRO*THOUSAND));
    case GAS_CO :
      return (MOLE_WT_CO/(AVOGADRO*THOUSAND));
    case GAS_CO2 :
      return (MOLE_WT_CO2/(AVOGADRO*THOUSAND));
    case GAS_CH4 :
      return (MOLE_WT_CH4/(AVOGADRO*THOUSAND));
    case GAS_H :
      return (MOLE_WT_H/(AVOGADRO*THOUSAND));
    case GAS_H2 :
      return (MOLE_WT_H2/(AVOGADRO*THOUSAND));
    case GAS_HE :
      return (MOLE_WT_HE/(AVOGADRO*THOUSAND));
    case GAS_H2O :
      return (MOLE_WT_H2O/(AVOGADRO*THOUSAND));
    case GAS_N2 :
      return (MOLE_WT_N2/(AVOGADRO*THOUSAND));
    case GAS_O :
      return (MOLE_WT_O/(AVOGADRO*THOUSAND));
    case GAS_O2 :
      return (MOLE_WT_O2/(AVOGADRO*THOUSAND));
    case GAS_e :
      return (MOLE_WT_e/(AVOGADRO*THOUSAND));
    case GAS_APHTPB :
      return (MOLE_WT_APHTPB/(AVOGADRO*THOUSAND));
    default:
      return (MOLE_WT_AIR/(AVOGADRO*THOUSAND));
  };
}

inline double mn(const Ion5Moment2D_pState &W) {
  switch(W.gas) {
    case GAS_AIR :
      return (MOLE_WT_AIR/(AVOGADRO*THOUSAND));
    case GAS_A :
      return (MOLE_WT_A/(AVOGADRO*THOUSAND));
    case GAS_CO :
      return (MOLE_WT_CO/(AVOGADRO*THOUSAND));
    case GAS_CO2 :
      return (MOLE_WT_CO2/(AVOGADRO*THOUSAND));
    case GAS_CH4 :
      return (MOLE_WT_CH4/(AVOGADRO*THOUSAND));
    case GAS_H :
      return (MOLE_WT_H/(AVOGADRO*THOUSAND));
    case GAS_H2 :
      return (MOLE_WT_H2/(AVOGADRO*THOUSAND));
    case GAS_HE :
      return (MOLE_WT_HE/(AVOGADRO*THOUSAND));
    case GAS_H2O :
      return (MOLE_WT_H2O/(AVOGADRO*THOUSAND));
    case GAS_N2 :
      return (MOLE_WT_N2/(AVOGADRO*THOUSAND));
    case GAS_O :
      return (MOLE_WT_O/(AVOGADRO*THOUSAND));
    case GAS_O2 :
      return (MOLE_WT_O2/(AVOGADRO*THOUSAND));
    case GAS_e :
      return (MOLE_WT_e/(AVOGADRO*THOUSAND));
    case GAS_APHTPB :
      return (MOLE_WT_APHTPB/(AVOGADRO*THOUSAND));
    default:
      return (MOLE_WT_AIR/(AVOGADRO*THOUSAND));
  };
}

/*******************************************************************
 * Ion5Moment2D_pState::Sc -- Ion-neutral collisions source terms. *
 *******************************************************************/
inline Ion5Moment2D_cState Ion5Moment2D_pState::Sc(const Euler2D_pState &Wneut) {
  double nu_in, mneut, dpdt; Vector2D dvdt;
  mneut = mn();
  nu_in = nu(Wneut);
  dvdt = (nu_in*mneut/(m+mneut))*(Wneut.v-v);
  dpdt = (nu_in*d*mneut/sqr(m+mneut))*(TWO*BOLTZMANN*(Wneut.T()-T())+TWO*mneut*sqr(Wneut.v-v)/THREE);
  return (Ion5Moment2D_cState(ZERO, d*dvdt.x, d*dvdt.y, gm1i*dpdt+d*(v*dvdt)));
}

inline Ion5Moment2D_cState Ion5Moment2D_pState::Sc(const Euler2D_pState &Wneut) const {
  double nu_in, mneut, dpdt; Vector2D dvdt;
  mneut = mn();
  nu_in = nu(Wneut);
  dvdt = (nu_in*mneut/(m+mneut))*(Wneut.v-v);
  dpdt = (nu_in*d*mneut/sqr(m+mneut))*(TWO*BOLTZMANN*(Wneut.T()-T())+TWO*mneut*sqr(Wneut.v-v)/THREE);
  return (Ion5Moment2D_cState(ZERO, d*dvdt.x, d*dvdt.y, gm1i*dpdt+d*(v*dvdt)));
}

inline Ion5Moment2D_cState Sc(const Ion5Moment2D_pState &W, const Euler2D_pState &Wneut) {
  double nu_in, mneut, dpdt; Vector2D dvdt;
  mneut = W.mn();
  nu_in = W.nu(Wneut);
  dvdt = (nu_in*mneut/(W.m+mneut))*(Wneut.v-W.v);
  dpdt = (nu_in*W.d*mneut/sqr(W.m+mneut))*(TWO*BOLTZMANN*(Wneut.T()-W.T())+TWO*mneut*sqr(Wneut.v-W.v)/THREE);
  return (Ion5Moment2D_cState(ZERO, W.d*dvdt.x, W.d*dvdt.y, W.gm1i*dpdt+W.d*(W.v*dvdt)));
}

inline void Ion5Moment2D_pState::dScdU(DenseMatrix &dScdU, const Euler2D_pState &Wneut) {
  double nu_in, mneut; mneut = mn(); nu_in = nu(Wneut);
  dScdU(1,0) += nu_in*mneut*Wneut.v.x/(m+mneut);
  dScdU(1,1) -= nu_in*mneut/(m+mneut);
  dScdU(2,0) += nu_in*mneut*Wneut.v.y/(m+mneut);
  dScdU(2,2) -= nu_in*mneut/(m+mneut);
  dScdU(3,0) += nu_in*mneut*(SIX*BOLTZMANN*Wneut.T()+TWO*mneut*(sqr(Wneut.v)-sqr(v))+
                THREE*gm1*mneut*sqr(v))/(THREE*gm1*sqr(m+mneut));
  dScdU(3,1) += nu_in*mneut*(THREE*gm1*(m+mneut)*Wneut.v.x-
                SIX*gm1*mneut*v.x-FOUR*mneut*(Wneut.v.x-v.x))/(THREE*gm1*sqr(m+mneut));
  dScdU(3,2) += nu_in*mneut*(THREE*gm1*(m+mneut)*Wneut.v.y-
                SIX*gm1*mneut*v.y-FOUR*mneut*(Wneut.v.y-v.y))/(THREE*gm1*sqr(m+mneut));
  dScdU(3,3) -= TWO*nu_in*m*mneut/sqr(m+mneut);
}

inline void Ion5Moment2D_pState::dScdU(DenseMatrix &dScdU, const Euler2D_pState &Wneut) const {
  double nu_in, mneut; mneut = mn(); nu_in = nu(Wneut);
  dScdU(1,0) += nu_in*mneut*Wneut.v.x/(m+mneut);
  dScdU(1,1) -= nu_in*mneut/(m+mneut);
  dScdU(2,0) += nu_in*mneut*Wneut.v.y/(m+mneut);
  dScdU(2,2) -= nu_in*mneut/(m+mneut);
  dScdU(3,0) += nu_in*mneut*(SIX*BOLTZMANN*Wneut.T()+TWO*mneut*(sqr(Wneut.v)-sqr(v))+
                THREE*gm1*mneut*sqr(v))/(THREE*gm1*sqr(m+mneut));
  dScdU(3,1) += nu_in*mneut*(THREE*gm1*(m+mneut)*Wneut.v.x-
                SIX*gm1*mneut*v.x-FOUR*mneut*(Wneut.v.x-v.x))/(THREE*gm1*sqr(m+mneut));
  dScdU(3,2) += nu_in*mneut*(THREE*gm1*(m+mneut)*Wneut.v.y-
                SIX*gm1*mneut*v.y-FOUR*mneut*(Wneut.v.y-v.y))/(THREE*gm1*sqr(m+mneut));
  dScdU(3,3) -= TWO*nu_in*m*mneut/sqr(m+mneut);
}

inline void dScdU(DenseMatrix &dScdU, const Ion5Moment2D_pState &W, const Euler2D_pState &Wneut) {
  double nu_in, mneut; mneut = W.mn(); nu_in = W.nu(Wneut);
  dScdU(1,0) += nu_in*mneut*Wneut.v.x/(W.m+mneut);
  dScdU(1,1) -= nu_in*mneut/(W.m+mneut);
  dScdU(2,0) += nu_in*mneut*Wneut.v.y/(W.m+mneut);
  dScdU(2,2) -= nu_in*mneut/(W.m+mneut);
  dScdU(3,0) += nu_in*mneut*(SIX*BOLTZMANN*Wneut.T()+TWO*mneut*(sqr(Wneut.v)-sqr(W.v))+
                THREE*W.gm1*mneut*sqr(W.v))/(THREE*W.gm1*sqr(W.m+mneut));
  dScdU(3,1) += nu_in*mneut*(THREE*W.gm1*(W.m+mneut)*Wneut.v.x-
                SIX*W.gm1*mneut*W.v.x-FOUR*mneut*(Wneut.v.x-W.v.x))/(THREE*W.gm1*sqr(W.m+mneut));
  dScdU(3,2) += nu_in*mneut*(THREE*W.gm1*(W.m+mneut)*Wneut.v.y-
                SIX*W.gm1*mneut*W.v.y-FOUR*mneut*(Wneut.v.y-W.v.y))/(THREE*W.gm1*sqr(W.m+mneut));
  dScdU(3,3) -= TWO*nu_in*W.m*mneut/sqr(W.m+mneut);
}

/******************************************************************
 * Ion5Moment2D_cState::Ion5Moment2D_cState -- Constructor.       *
 ******************************************************************/
inline Ion5Moment2D_cState::Ion5Moment2D_cState(const Ion5Moment2D_pState &W) {
  d = W.d; dv = W.dv(); E = W.E();
}

/*************************************************************
 * Ion5Moment2D_cState::W -- Primitive solution state.       *
 *************************************************************/
inline Ion5Moment2D_pState Ion5Moment2D_cState::W(void) {
  return (Ion5Moment2D_pState(d, v(), p()));
}

inline Ion5Moment2D_pState Ion5Moment2D_cState::W(void) const {
  return (Ion5Moment2D_pState(d, v(), p()));
}

inline Ion5Moment2D_pState Ion5Moment2D_cState::W(const Ion5Moment2D_cState &U) {
  return (Ion5Moment2D_pState(U.d, U.v(), U.p()));
}

inline Ion5Moment2D_pState W(const Ion5Moment2D_cState &U) {
  return (Ion5Moment2D_pState(U.d, U.v(), U.p()));
}

/*************************************************************
 * Ion5Moment2D_cState::F -- Solution flux (x-direction).    *
 *************************************************************/
inline Ion5Moment2D_cState Ion5Moment2D_cState::F(void) {
  return (Ion5Moment2D_cState(dv.x, sqr(dv.x)/d + p(),
                              dv.x*dv.y/d, dv.x*H()/d));
}

inline Ion5Moment2D_cState Ion5Moment2D_cState::F(void) const {
  return (Ion5Moment2D_cState(dv.x, sqr(dv.x)/d + p(),
                              dv.x*dv.y/d, dv.x*H()/d));
}

inline Ion5Moment2D_cState Ion5Moment2D_cState::F(const Ion5Moment2D_cState &U) {
  return (Ion5Moment2D_cState(U.dv.x, sqr(U.dv.x)/U.d + U.p(),
                              U.dv.x*U.dv.y/U.d, U.dv.x*U.H()/U.d));
}

inline Ion5Moment2D_cState F(const Ion5Moment2D_cState &U) {
  return (Ion5Moment2D_cState(U.dv.x, sqr(U.dv.x)/U.d + U.p(),
                              U.dv.x*U.dv.y/U.d, U.dv.x*U.H()/U.d));
}

/*************************************************************
 * Ion5Moment2D_cState::Fx -- Solution flux (x-direction).   *
 *************************************************************/
inline Ion5Moment2D_cState Ion5Moment2D_cState::Fx(void) {
  return (Ion5Moment2D_cState(dv.x, sqr(dv.x)/d + p(),
                              dv.x*dv.y/d, dv.x*H()/d));
}

inline Ion5Moment2D_cState Ion5Moment2D_cState::Fx(void) const {
  return (Ion5Moment2D_cState(dv.x, sqr(dv.x)/d + p(),
                              dv.x*dv.y/d, dv.x*H()/d));
}

inline Ion5Moment2D_cState Ion5Moment2D_cState::Fx(const Ion5Moment2D_cState &U) {
  return (Ion5Moment2D_cState(U.dv.x, sqr(U.dv.x)/U.d + U.p(),
                              U.dv.x*U.dv.y/U.d, U.dv.x*U.H()/U.d));
}

inline Ion5Moment2D_cState Fx(const Ion5Moment2D_cState &U) {
  return (Ion5Moment2D_cState(U.dv.x, sqr(U.dv.x)/U.d + U.p(),
                              U.dv.x*U.dv.y/U.d, U.dv.x*U.H()/U.d));
}

/*************************************************************
 * Ion5Moment2D_cState::Fy -- Solution flux (y-direction).   *
 *************************************************************/
inline Ion5Moment2D_cState Ion5Moment2D_cState::Fy(void) {
  return (Ion5Moment2D_cState(dv.y, dv.x*dv.y/d,
		  	      sqr(dv.y)/d + p(), dv.y*H()/d));
}

inline Ion5Moment2D_cState Ion5Moment2D_cState::Fy(void) const {
  return (Ion5Moment2D_cState(dv.y, dv.x*dv.y/d,
		  	      sqr(dv.y)/d + p(), dv.y*H()/d));
}

inline Ion5Moment2D_cState Ion5Moment2D_cState::Fy(const Ion5Moment2D_cState &U) {
  return (Ion5Moment2D_cState(U.dv.y, U.dv.x*U.dv.y/U.d,
		  	      sqr(U.dv.y)/U.d + U.p(), U.dv.y*U.H()/U.d));
}

inline Ion5Moment2D_cState Fy(const Ion5Moment2D_cState &U) {
  return (Ion5Moment2D_cState(U.dv.y, U.dv.x*U.dv.y/U.d,
		 	      sqr(U.dv.y)/U.d + U.p(), U.dv.y*U.H()/U.d));
}

/*************************************************************
 * Ion5Moment2D_cState::Fn -- Solution flux (n-direction).   *
 *************************************************************/
inline Ion5Moment2D_cState Ion5Moment2D_cState::Fn(void) {
  return (Ion5Moment2D_cState(dv.x, sqr(dv.x)/d + p(),
                              dv.x*dv.y/d, dv.x*H()/d));
}

inline Ion5Moment2D_cState Ion5Moment2D_cState::Fn(void) const {
  return (Ion5Moment2D_cState(dv.x, sqr(dv.x)/d + p(),
                              dv.x*dv.y/d, dv.x*H()/d));
}

inline Ion5Moment2D_cState Ion5Moment2D_cState::Fn(const Ion5Moment2D_cState &U) {
  return (Ion5Moment2D_cState(U.dv.x, sqr(U.dv.x)/U.d + U.p(),
                              U.dv.x*U.dv.y/U.d, U.dv.x*U.H()/U.d));
}

inline Ion5Moment2D_cState Fn(const Ion5Moment2D_cState &U) {
  return (Ion5Moment2D_cState(U.dv.x, sqr(U.dv.x)/U.d + U.p(),
                              U.dv.x*U.dv.y/U.d, U.dv.x*U.H()/U.d));
}

/***************************************************************
 * Ion5Moment2D_cState::Sa -- Axisymmetric flow source terms.  *
 ***************************************************************/
inline Ion5Moment2D_cState Ion5Moment2D_cState::Sa(const Vector2D &X) {
  return (Ion5Moment2D_cState(-dv.y/X.y, -dv.x*dv.y/(d*X.y), 
                              -sqr(dv.y)/(d*X.y), -dv.y*H()/(d*X.y)));
}

inline Ion5Moment2D_cState Ion5Moment2D_cState::Sa(const Vector2D &X) const {
  return (Ion5Moment2D_cState(-dv.y/X.y, -dv.x*dv.y/(d*X.y), 
                              -sqr(dv.y)/(d*X.y), -dv.y*H()/(d*X.y)));
}

inline Ion5Moment2D_cState Sa(const Ion5Moment2D_cState &U, const Vector2D &X) {
  return (Ion5Moment2D_cState(-U.dv.y/X.y, -U.dv.x*U.dv.y/(U.d*X.y), 
                              -sqr(U.dv.y)/(U.d*X.y), -U.dv.y*U.H()/(U.d*X.y)));
}

/***************************************************************
 * Ion5Moment2D_cState::Se -- Electric field source terms.     *
 ***************************************************************/
inline Ion5Moment2D_cState Ion5Moment2D_cState::Se(const Vector2D &E) {
  double cm; cm = d*(q/m);
  return (Ion5Moment2D_cState(ZERO, cm*E.x, cm*E.y, cm*(dv*E)/d));
}

inline Ion5Moment2D_cState Ion5Moment2D_cState::Se(const Vector2D &E) const {
  double cm; cm = d*(q/m);
  return (Ion5Moment2D_cState(ZERO, cm*E.x, cm*E.y, cm*(dv*E)/d));
}

inline Ion5Moment2D_cState Se(const Ion5Moment2D_cState &U, const Vector2D &E) {
  double cm; cm = U.d*(U.q/U.m);
  return (Ion5Moment2D_cState(ZERO, cm*E.x, cm*E.y, cm*(U.dv*E)/U.d));
}

/***************************************************************
 * Ion5Moment2D_cState::nu -- Ion-neutral collision frequency. *
 ***************************************************************/
inline double Ion5Moment2D_cState::nu(const Euler2D_pState &Wneut) {
  double nu_in, mneut, sigma;
  mneut = mn();
  switch(ion) {
    case ION_H :
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
    case ION_HE :
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
    case ION_O :
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
    case ION_LOW_MASS :
      sigma = 105*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*T())/(m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    case ION_MED_MASS :
      sigma = 150*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*T())/(m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    case ION_HI_MASS :
      sigma = 280*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*T())/(m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
   case ION_HUGE_MASS :
      sigma = 2560*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*T())/(m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
   case ION_HUGE_MASS2 :
      sigma = 433*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*T())/(m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    default:
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
  };
  return (nu_in);
}

inline double Ion5Moment2D_cState::nu(const Euler2D_pState &Wneut) const {
  double nu_in, mneut, sigma;
  mneut = mn();
  switch(ion) {
    case ION_H :
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
    case ION_HE :
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
    case ION_O :
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
    case ION_LOW_MASS :
      sigma = 105*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*T())/(m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    case ION_MED_MASS :
      sigma = 150*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*T())/(m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    case ION_HI_MASS :
      sigma = 280*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*T())/(m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
   case ION_HUGE_MASS :
      sigma = 2560*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*T())/(m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
   case ION_HUGE_MASS2 :
      sigma = 433*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*T())/(m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    default:
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
  };
  return (nu_in);
}

inline double nu(const Ion5Moment2D_cState &U, const Euler2D_pState &Wneut) {
  double nu_in, mneut, sigma;
  mneut = U.mn();
  switch(U.ion) {
    case ION_H :
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
    case ION_HE :
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
    case ION_O :
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
    case ION_LOW_MASS :
      sigma = 105*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              (ONE+sqrt((mneut*U.T())/(U.m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    case ION_MED_MASS :
      sigma = 150*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*U.T())/(U.m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    case ION_HI_MASS :
      sigma = 280*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*U.T())/(U.m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    case ION_HUGE_MASS :
      sigma = 2560*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*U.T())/(U.m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    case ION_HUGE_MASS2 :
      sigma = 433*1.0e-20;
      nu_in = sigma*(Wneut.d/mneut)*
              sqrt(ONE+((mneut*U.T())/(U.m*Wneut.T())))*
              sqrt(EIGHT*BOLTZMANN*Wneut.T()/(PI*mneut));
      break;
    default:
      nu_in = 3.36e-09*(Wneut.d/(mneut*MILLION));
      break;
  };
  return (nu_in);
}

/***************************************************************
 * Ion5Moment2D_cState::mn -- Neutral gas particle mass.       *
 ***************************************************************/
inline double Ion5Moment2D_cState::mn(void) {
  switch(gas) {
    case GAS_AIR :
      return (MOLE_WT_AIR/(AVOGADRO*THOUSAND));
    case GAS_A :
      return (MOLE_WT_A/(AVOGADRO*THOUSAND));
    case GAS_CO :
      return (MOLE_WT_CO/(AVOGADRO*THOUSAND));
    case GAS_CO2 :
      return (MOLE_WT_CO2/(AVOGADRO*THOUSAND));
    case GAS_CH4 :
      return (MOLE_WT_CH4/(AVOGADRO*THOUSAND));
    case GAS_H :
      return (MOLE_WT_H/(AVOGADRO*THOUSAND));
    case GAS_H2 :
      return (MOLE_WT_H2/(AVOGADRO*THOUSAND));
    case GAS_HE :
      return (MOLE_WT_HE/(AVOGADRO*THOUSAND));
    case GAS_H2O :
      return (MOLE_WT_H2O/(AVOGADRO*THOUSAND));
    case GAS_N2 :
      return (MOLE_WT_N2/(AVOGADRO*THOUSAND));
    case GAS_O :
      return (MOLE_WT_O/(AVOGADRO*THOUSAND));
    case GAS_O2 :
      return (MOLE_WT_O2/(AVOGADRO*THOUSAND));
    case GAS_e :
      return (MOLE_WT_e/(AVOGADRO*THOUSAND));
    case GAS_APHTPB :
      return (MOLE_WT_APHTPB/(AVOGADRO*THOUSAND));
    default:
      return (MOLE_WT_AIR/(AVOGADRO*THOUSAND));
  };
}

inline double Ion5Moment2D_cState::mn(void) const {
  switch(gas) {
    case GAS_AIR :
      return (MOLE_WT_AIR/(AVOGADRO*THOUSAND));
    case GAS_A :
      return (MOLE_WT_A/(AVOGADRO*THOUSAND));
    case GAS_CO :
      return (MOLE_WT_CO/(AVOGADRO*THOUSAND));
    case GAS_CO2 :
      return (MOLE_WT_CO2/(AVOGADRO*THOUSAND));
    case GAS_CH4 :
      return (MOLE_WT_CH4/(AVOGADRO*THOUSAND));
    case GAS_H :
      return (MOLE_WT_H/(AVOGADRO*THOUSAND));
    case GAS_H2 :
      return (MOLE_WT_H2/(AVOGADRO*THOUSAND));
    case GAS_HE :
      return (MOLE_WT_HE/(AVOGADRO*THOUSAND));
    case GAS_H2O :
      return (MOLE_WT_H2O/(AVOGADRO*THOUSAND));
    case GAS_N2 :
      return (MOLE_WT_N2/(AVOGADRO*THOUSAND));
    case GAS_O :
      return (MOLE_WT_O/(AVOGADRO*THOUSAND));
    case GAS_O2 :
      return (MOLE_WT_O2/(AVOGADRO*THOUSAND));
    case GAS_e :
      return (MOLE_WT_e/(AVOGADRO*THOUSAND));
    case GAS_APHTPB :
      return (MOLE_WT_APHTPB/(AVOGADRO*THOUSAND));
    default:
      return (MOLE_WT_AIR/(AVOGADRO*THOUSAND));
  };
}

inline double mn(const Ion5Moment2D_cState &U) {
  switch(U.gas) {
    case GAS_AIR :
      return (MOLE_WT_AIR/(AVOGADRO*THOUSAND));
    case GAS_A :
      return (MOLE_WT_A/(AVOGADRO*THOUSAND));
    case GAS_CO :
      return (MOLE_WT_CO/(AVOGADRO*THOUSAND));
    case GAS_CO2 :
      return (MOLE_WT_CO2/(AVOGADRO*THOUSAND));
    case GAS_CH4 :
      return (MOLE_WT_CH4/(AVOGADRO*THOUSAND));
    case GAS_H :
      return (MOLE_WT_H/(AVOGADRO*THOUSAND));
    case GAS_H2 :
      return (MOLE_WT_H2/(AVOGADRO*THOUSAND));
    case GAS_HE :
      return (MOLE_WT_HE/(AVOGADRO*THOUSAND));
    case GAS_H2O :
      return (MOLE_WT_H2O/(AVOGADRO*THOUSAND));
    case GAS_N2 :
      return (MOLE_WT_N2/(AVOGADRO*THOUSAND));
    case GAS_O :
      return (MOLE_WT_O/(AVOGADRO*THOUSAND));
    case GAS_O2 :
      return (MOLE_WT_O2/(AVOGADRO*THOUSAND));
    case GAS_e :
      return (MOLE_WT_e/(AVOGADRO*THOUSAND));
    case GAS_APHTPB :
      return (MOLE_WT_APHTPB/(AVOGADRO*THOUSAND));
    default:
      return (MOLE_WT_AIR/(AVOGADRO*THOUSAND));
  };
}

/*******************************************************************
 * Ion5Moment2D_cState::Sc -- Ion-neutral collisions source terms. *
 *******************************************************************/
inline Ion5Moment2D_cState Ion5Moment2D_cState::Sc(const Euler2D_pState &Wneut) {
  return (W().Sc(Wneut));
}

inline Ion5Moment2D_cState Ion5Moment2D_cState::Sc(const Euler2D_pState &Wneut) const {
  return (W().Sc(Wneut));
}

inline Ion5Moment2D_cState Sc(const Ion5Moment2D_cState &U, const Euler2D_pState &Wneut) {
  return (U.W().Sc(Wneut));
}

inline void Ion5Moment2D_cState::dScdU(DenseMatrix &dScdU, const Euler2D_pState &Wneut) {
  W().dScdU(dScdU, Wneut);
}

inline void Ion5Moment2D_cState::dScdU(DenseMatrix &dScdU, const Euler2D_pState &Wneut) const {
  W().dScdU(dScdU, Wneut);
}

inline void dScdU(DenseMatrix &dScdU, const Ion5Moment2D_cState &U, const Euler2D_pState &Wneut) {
  U.W().dScdU(dScdU, Wneut);
}

/*************************************************************
 * Useful 2D Ion 5-moment state constants.                   *
 *************************************************************/
const Ion5Moment2D_pState Ion5Moment2D_W_STDATM(DENSITY_STDATM,
 				                Vector2D_ZERO, PRESSURE_STDATM);
const Ion5Moment2D_pState Ion5Moment2D_W_VACUUM(ZERO, Vector2D_ZERO, ZERO);
const Ion5Moment2D_pState Ion5Moment2D_W_ONE(ONE, ONE, ONE, ONE);
const Ion5Moment2D_cState Ion5Moment2D_U_STDATM(Ion5Moment2D_W_STDATM);
const Ion5Moment2D_cState Ion5Moment2D_U_VACUUM(Ion5Moment2D_W_VACUUM);
const Ion5Moment2D_cState Ion5Moment2D_U_ONE(ONE, ONE, ONE, ONE);

/*************************************************************
 * Ion5Moment2DState -- External subroutines.                *
 *************************************************************/

extern Ion5Moment2D_pState Riemann(const Ion5Moment2D_pState &Wl,
                                   const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_pState Riemann_x(const Ion5Moment2D_pState &Wl,
                                     const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_pState Riemann_y(const Ion5Moment2D_pState &Wl,
                                     const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_pState RoeAverage(const Ion5Moment2D_pState &Wl,
                                      const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_pState Reflect(const Ion5Moment2D_pState &W,
	      	                   const Vector2D &norm_dir);

extern Ion5Moment2D_pState BC_Characteristic(const Ion5Moment2D_pState &Wi,
                                             const Ion5Moment2D_pState &Wo,
	      	                             const Vector2D &norm_dir);

extern Ion5Moment2D_pState BC_Characteristic_Pressure(const Ion5Moment2D_pState &Wi,
                                                      const Ion5Moment2D_pState &Wo,
	      	                                      const Vector2D &norm_dir);

extern Ion5Moment2D_pState BC_Characteristic_Mach_Number(const Ion5Moment2D_pState &Wi,
                                                         const Ion5Moment2D_pState &Wo,
	      	                                         const Vector2D &norm_dir);

extern Ion5Moment2D_pState WaveSpeedPos(const Ion5Moment2D_pState &lambda_a,
                                        const Ion5Moment2D_pState &lambda_l,
                                        const Ion5Moment2D_pState &lambda_r);

extern Ion5Moment2D_pState WaveSpeedNeg(const Ion5Moment2D_pState &lambda_a,
                                        const Ion5Moment2D_pState &lambda_l,
                                        const Ion5Moment2D_pState &lambda_r);

extern Ion5Moment2D_pState WaveSpeedAbs(const Ion5Moment2D_pState &lambda_a,
                                        const Ion5Moment2D_pState &lambda_l,
                                        const Ion5Moment2D_pState &lambda_r);

extern Ion5Moment2D_pState HartenFixPos(const Ion5Moment2D_pState &lambda_a,
                                        const Ion5Moment2D_pState &lambda_l,
                                        const Ion5Moment2D_pState &lambda_r);

extern Ion5Moment2D_pState HartenFixNeg(const Ion5Moment2D_pState &lambda_a,
                                        const Ion5Moment2D_pState &lambda_l,
                                        const Ion5Moment2D_pState &lambda_r);

extern Ion5Moment2D_pState HartenFixAbs(const Ion5Moment2D_pState &lambda_a,
                                        const Ion5Moment2D_pState &lambda_l,
                                        const Ion5Moment2D_pState &lambda_r);

extern Ion5Moment2D_cState FluxGodunov(const Ion5Moment2D_pState &Wl,
	      	                       const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_cState FluxGodunov(const Ion5Moment2D_cState &Ul,
	      	                       const Ion5Moment2D_cState &Ur);

extern Ion5Moment2D_cState FluxGodunov_x(const Ion5Moment2D_pState &Wl,
	      	                         const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_cState FluxGodunov_x(const Ion5Moment2D_cState &Ul,
	      	                         const Ion5Moment2D_cState &Ur);

extern Ion5Moment2D_cState FluxGodunov_y(const Ion5Moment2D_pState &Wl,
	      	                         const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_cState FluxGodunov_y(const Ion5Moment2D_cState &Ul,
	      	                         const Ion5Moment2D_cState &Ur);

extern Ion5Moment2D_cState FluxGodunov_n(const Ion5Moment2D_pState &Wl,
	      	                         const Ion5Moment2D_pState &Wr,
                                         const Vector2D &norm_dir);

extern Ion5Moment2D_cState FluxGodunov_n(const Ion5Moment2D_cState &Ul,
	      	                         const Ion5Moment2D_cState &Ur,
                                         const Vector2D &norm_dir);

extern Ion5Moment2D_cState FluxRoe(const Ion5Moment2D_pState &Wl,
	      	                   const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_cState FluxRoe(const Ion5Moment2D_cState &Ul,
	      	                   const Ion5Moment2D_cState &Ur);

extern Ion5Moment2D_cState FluxRoe_x(const Ion5Moment2D_pState &Wl,
	      	                     const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_cState FluxRoe_x(const Ion5Moment2D_cState &Ul,
	      	                     const Ion5Moment2D_cState &Ur);

extern Ion5Moment2D_cState FluxRoe_y(const Ion5Moment2D_pState &Wl,
	      	                     const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_cState FluxRoe_y(const Ion5Moment2D_cState &Ul,
	      	                     const Ion5Moment2D_cState &Ur);

extern Ion5Moment2D_cState FluxRoe_n(const Ion5Moment2D_pState &Wl,
	      	                     const Ion5Moment2D_pState &Wr,
                                     const Vector2D &norm_dir);

extern Ion5Moment2D_cState FluxRoe_n(const Ion5Moment2D_cState &Ul,
	      	                     const Ion5Moment2D_cState &Ur,
                                     const Vector2D &norm_dir);

extern Ion5Moment2D_cState FluxRusanov(const Ion5Moment2D_pState &Wl,
	      	                       const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_cState FluxRusanov(const Ion5Moment2D_cState &Ul,
	      	                       const Ion5Moment2D_cState &Ur);

extern Ion5Moment2D_cState FluxRusanov_x(const Ion5Moment2D_pState &Wl,
	      	                         const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_cState FluxRusanov_x(const Ion5Moment2D_cState &Ul,
	      	                         const Ion5Moment2D_cState &Ur);

extern Ion5Moment2D_cState FluxRusanov_y(const Ion5Moment2D_pState &Wl,
	      	                         const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_cState FluxRusanov_y(const Ion5Moment2D_cState &Ul,
	      	                         const Ion5Moment2D_cState &Ur);

extern Ion5Moment2D_cState FluxRusanov_n(const Ion5Moment2D_pState &Wl,
	      	                         const Ion5Moment2D_pState &Wr,
                                         const Vector2D &norm_dir);

extern Ion5Moment2D_cState FluxRusanov_n(const Ion5Moment2D_cState &Ul,
	      	                         const Ion5Moment2D_cState &Ur,
                                         const Vector2D &norm_dir);

extern Ion5Moment2D_cState FluxHLLE(const Ion5Moment2D_pState &Wl,
	      	                    const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_cState FluxHLLE(const Ion5Moment2D_cState &Ul,
	      	                    const Ion5Moment2D_cState &Ur);

extern Ion5Moment2D_cState FluxHLLE_x(const Ion5Moment2D_pState &Wl,
	      	                      const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_cState FluxHLLE_x(const Ion5Moment2D_cState &Ul,
	      	                      const Ion5Moment2D_cState &Ur);

extern Ion5Moment2D_cState FluxHLLE_y(const Ion5Moment2D_pState &Wl,
	      	                      const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_cState FluxHLLE_y(const Ion5Moment2D_cState &Ul,
	      	                      const Ion5Moment2D_cState &Ur);
  
extern Ion5Moment2D_cState FluxHLLE_n(const Ion5Moment2D_pState &Wl,
	      	                      const Ion5Moment2D_pState &Wr,
                                      const Vector2D &norm_dir);

extern Ion5Moment2D_cState FluxHLLE_n(const Ion5Moment2D_cState &Ul,
	      	                      const Ion5Moment2D_cState &Ur,
                                      const Vector2D &norm_dir);

extern Ion5Moment2D_cState FluxLinde(const Ion5Moment2D_pState &Wl,
	      	                     const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_cState FluxLinde(const Ion5Moment2D_cState &Ul,
	      	                     const Ion5Moment2D_cState &Ur);

extern Ion5Moment2D_cState FluxLinde_x(const Ion5Moment2D_pState &Wl,
	      	                       const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_cState FluxLinde_x(const Ion5Moment2D_cState &Ul,
	      	                       const Ion5Moment2D_cState &Ur);

extern Ion5Moment2D_cState FluxLinde_y(const Ion5Moment2D_pState &Wl,
	      	                       const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_cState FluxLinde_y(const Ion5Moment2D_cState &Ul,
	      	                       const Ion5Moment2D_cState &Ur);

extern Ion5Moment2D_cState FluxLinde_n(const Ion5Moment2D_pState &Wl,
	      	                       const Ion5Moment2D_pState &Wr,
                                       const Vector2D &norm_dir);

extern Ion5Moment2D_cState FluxLinde_n(const Ion5Moment2D_cState &Ul,
	         	               const Ion5Moment2D_cState &Ur,
                                       const Vector2D &norm_dir);

extern Ion5Moment2D_cState FluxHLLC(const Ion5Moment2D_pState &Wl,
	      	                    const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_cState FluxHLLC(const Ion5Moment2D_cState &Ul,
	      	                    const Ion5Moment2D_cState &Ur);

extern Ion5Moment2D_cState FluxHLLC_x(const Ion5Moment2D_pState &Wl,
	      	                      const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_cState FluxHLLC_x(const Ion5Moment2D_cState &Ul,
	      	                      const Ion5Moment2D_cState &Ur);

extern Ion5Moment2D_cState FluxHLLC_y(const Ion5Moment2D_pState &Wl,
	      	                      const Ion5Moment2D_pState &Wr);

extern Ion5Moment2D_cState FluxHLLC_y(const Ion5Moment2D_cState &Ul,
	      	                      const Ion5Moment2D_cState &Ur);

extern Ion5Moment2D_cState FluxHLLC_n(const Ion5Moment2D_pState &Wl,
	      	                      const Ion5Moment2D_pState &Wr,
                                      const Vector2D &norm_dir);

extern Ion5Moment2D_cState FluxHLLC_n(const Ion5Moment2D_cState &Ul,
	      	                      const Ion5Moment2D_cState &Ur,
                                      const Vector2D &norm_dir);

#endif /* _ION5MOMENT2D_STATE_INCLUDED  */
