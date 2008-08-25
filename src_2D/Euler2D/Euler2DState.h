/*!\file Euler2DState.h
  \brief Header file defining 2D Euler Solution State Classes. */

#ifndef _EULER2D_STATE_INCLUDED
#define _EULER2D_STATE_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

/* Using std namespace functions */
using namespace std;

/* Include CFFC header files */
#include "../Math/Math.h"  /* Include math macro header files. */
#include "../CFD/CFD.h"    /* Include CFD header files. */
#include "../Math/Matrix.h"  /* Include matrix header files. */
#include "../Math/Vector2D.h" /* Include vector 2D header files. */
#include "../Physics/GasConstants.h" /* Include gas constant header files. */
#include "../Physics/SolidConstants.h" /* Include solid constant header files. */

/* Define the classes. */

#define	NUM_VAR_EULER2D    4

class Euler2D_cState;

/*!
 * Class: Euler2D_pState
 *
 * @brief Primitive variable solution state class definition for an
 *        inviscid compressible gas-flow.
 *
 * Primitive variable solution state class definition for an inviscid
 * compressible gas-flow.
 *
 * \verbatim
 * Member functions
 *     d        -- Return density.
 *     v        -- Return flow velocity.
 *     p        -- Return pressure.
 *     g        -- Return specific heat ratio.
 *     gm1      -- Return g-1.
 *     gm1i     -- Return 1/(g-1).
 *     R        -- Return gas constant.
 *     setgas   -- Set gas constants.
 *     T        -- Return temperature.
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
 *     Mr_min   -- Reference Mach number used in Weiss-Smith low-Mach-
 *                 number local preconditioner.
 *     U        -- Return conserved solution state.
 *     F        -- Return x-direction solution flux.
 *     Fx       -- Return x-direction solution flux.
 *     Fy       -- Return y-direction solution flux.
 *     Fn       -- Return n-direction solution flux.
 *     dFdU     -- Return x-direction flux Jacobian.
 *     dFxdU    -- Return x-direction flux Jacobian.
 *     dFydU    -- Return y-direction flux Jacobian.
 *     dFndU    -- Return n-direction flux Jacobian.
 *     lambda   -- Return x-direction eigenvalue(s).
 *     lambda_x -- Return x-direction eigenvalue(s).
 *     lambda_y -- Return y-direction eigenvalue(s).
 *     lambda_precon_WS -- Return x-direction eigenvalue(s) of 
 *                 Weiss-Smith preconditioned system.
 *     rp       -- Return primitive right eigenvector (x-direction).
 *     rp_x     -- Return primitive right eigenvector (x-direction).
 *     rp_y     -- Return primitive right eigenvector (y-direction).
 *     rp_precon_WS -- Return primitive right eigenvector of 
 *                    Weiss-Smith preconditioned system (x-direction).
 *     rc       -- Return conserved right eigenvector (x-direction).
 *     rc_x     -- Return conserved right eigenvector (x-direction).
 *     rc_y     -- Return conserved right eigenvector (y-direction).
 *     rc_precon_WS -- Return conservative right eigenvector of
 *                     Weiss-Smith preconditioned system (x-direction).
 *     lp       -- Return primitive left eigenvector (x-direction).
 *     lp_x     -- Return primitive left eigenvector (x-direction).
 *     lp_y     -- Return primitive left eigenvector (y-direction).
 *     lp_precon_WS -- Return primitive left eigenvector of
 *                     Weiss-Smith preconditioned system (x-direction).
 *     S        -- Return axisymmetric source term vector.
 *     dUdW     -- Return Jacobian of conserved solution variables wrt
 *                 primitive solution variables.
 *     dWdU     -- Return Jacobian of primitive solution variables wrt
 *                 conserved solution variables.
 *     P_U_WS   -- Weiss-Smith preconditioner for the conservative
 *                 solution variables.
 *     P_U_WS_inv -- Inverse of Weiss-Smith preconditioner for the
 *                   conservative solution variables.
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
class Euler2D_pState{
private:
public:
  //@{ @name Primitive variables and associated constants:
  double             d; //!< Density.
  Vector2D           v; //!< Flow velocity (2D vector).
  double             p; //!< Pressure.
  static double      g; //!< Specific heat ratio.
  static double    gm1; //!< g-1
  static double   gm1i; //!< 1/(g-1)
  static double      R; //!< Gas constant.
  static double Mr_min; //!< Reference Mach number used in Weiss-Smith
                        //!< low-Mach-number local preconditioner.
  //@}

  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  Euler2D_pState(void) {
    d = DENSITY_STDATM; v.zero(); p = PRESSURE_STDATM;
  }

  //! Copy constructor.
  Euler2D_pState(const Euler2D_pState &W) {
    d = W.d; v = W.v; p = W.p;
  }

  //! Assignment constructor.
  Euler2D_pState(const Euler2D_cState &U);

  //! Assignment constructor.
  Euler2D_pState(const double &rho,
		 const Vector2D &V,
		 const double &pre) {
    d = rho; v = V; p = pre;
  }

  //! Assignment constructor.
  Euler2D_pState(const double &rho,
		 const double &vx,
		 const double &vy,
		 const double &pre) {
    d = rho; v.x = vx; v.y = vy; p = pre;
  }

  //! Value Constructor
  Euler2D_pState(const double &Val);

  /* Destructor. */
  // ~Euler2D_pState(void);
  // Use automatically generated destructor.
  //@}
 
  //@{ @name Useful operators.
  //! Return the number of variables.
  static int NumVar(void) { return NUM_VAR_EULER2D; }

  //! Copy operator.
  void Copy(const Euler2D_pState &W) {
    d = W.d; v = W.v; p = W.p;
  }

  //! Vacuum operator.
  void Vacuum(void) {
    d = ZERO; v = Vector2D_ZERO; p = ZERO;
  }

  //! One operator. Set the solution to ONE.
  void One(void) { 
    d = ONE; v.x = ONE; v.y = ONE; p = ONE;
  }

  //! Standard atmosphere operator.
  void Standard_Atmosphere(void) {
    d = DENSITY_STDATM; v.zero(); p = PRESSURE_STDATM;
  }

  //! Check for unphysical state properties.
  int Unphysical_Properties(void) const {
    if (d <= ZERO || p <= ZERO || E() <= ZERO) return 1;
    return 0;
  }
  //@}

  //@{ @name Set static variables.
  void setgas(void);
  void setgas(char *string_ptr);
  //@}

  //@{ @name State functions.
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

  //! Reference Mach number used by Weiss-Smith low-Mach-number local
  //! preconditioner.
  double Mr(void);
  double Mr(void) const;

  //! Burning rate.
  double burningrate(void) const;
  //@}

  //@{ @name Conserved solution state.
  Euler2D_cState U(void);
  Euler2D_cState U(void) const;
  Euler2D_cState U(const Euler2D_pState &W);
  friend Euler2D_cState U(const Euler2D_pState &W);
  //@}

  //@{ @name Solution flux and Jacobian (x-direction).
  Euler2D_cState F(void);
  Euler2D_cState F(void) const;
  Euler2D_cState F(const Euler2D_pState &W);
  friend Euler2D_cState F(const Euler2D_pState &W);
  Euler2D_cState F(const Vector2D &V) const;
  void dFdU(DenseMatrix &dFdU);
  void dFdU(DenseMatrix &dFdU) const;
  friend void dFdU(DenseMatrix &dFdU, const Euler2D_pState &W);

  Euler2D_cState Fx(void);
  Euler2D_cState Fx(void) const;
  Euler2D_cState Fx(const Euler2D_pState &W);
  friend Euler2D_cState Fx(const Euler2D_pState &W);
  void dFxdU(DenseMatrix &dFxdU);
  void dFxdU(DenseMatrix &dFxdU) const;
  friend void dFxdU(DenseMatrix &dFxdU, const Euler2D_pState &W);
  //@}

  //@{ @name Solution flux and Jacobian (y-direction).
  Euler2D_cState Fy(void);
  Euler2D_cState Fy(void) const;
  Euler2D_cState Fy(const Euler2D_pState &W);
  friend Euler2D_cState Fy(const Euler2D_pState &W);
  void dFydU(DenseMatrix &dFydU);
  void dFydU(DenseMatrix &dFydU) const;
  friend void dFydU(DenseMatrix &dFydU, const Euler2D_pState &W);
  //@}

  //@{ @name Solution flux and Jacobian (n-direction).
  Euler2D_cState Fn(void);
  Euler2D_cState Fn(void) const;
  Euler2D_cState Fn(const Euler2D_pState &W);
  friend Euler2D_cState Fn(const Euler2D_pState &W);
  void dFndU(DenseMatrix &dFndU);
  void dFndU(DenseMatrix &dFndU) const;
  friend void dFndU(DenseMatrix &dFndU, const Euler2D_pState &W);
  //@}

  //@{ @name Eigenvalue(s) (x-direction).
  Euler2D_pState lambda(void);
  Euler2D_pState lambda(void) const;
  Euler2D_pState lambda(const Euler2D_pState &W);
  friend Euler2D_pState lambda(const Euler2D_pState &W);
  double lambda(int index);
  double lambda(int index) const;
  friend double lambda(const Euler2D_pState &W, int index);

  Euler2D_pState lambda_x(void);
  Euler2D_pState lambda_x(void) const;
  Euler2D_pState lambda_x(const Euler2D_pState &W);
  friend Euler2D_pState lambda_x(const Euler2D_pState &W);
  double lambda_x(int index);
  double lambda_x(int index) const;
  friend double lambda_x(const Euler2D_pState &W, int index);
  Euler2D_pState lambda_x(const Vector2D &V) const;

  Euler2D_pState lambda_precon_WS(void);
  Euler2D_pState lambda_precon_WS(void) const;
  Euler2D_pState lambda_precon_WS(const Euler2D_pState &W);
  friend Euler2D_pState lambda_precon_WS(const Euler2D_pState &W);
  double lambda_precon_WS(int index);
  double lambda_precon_WS(int index) const;
  friend double lambda_precon_WS(const Euler2D_pState &W, int index);
  //@}

  //@{ @name Eigenvalue(s) (y-direction).
  Euler2D_pState lambda_y(void);
  Euler2D_pState lambda_y(void) const;
  Euler2D_pState lambda_y(const Euler2D_pState &W);
  friend Euler2D_pState lambda_y(const Euler2D_pState &W);
  double lambda_y(int index);
  double lambda_y(int index) const;
  friend double lambda_y(const Euler2D_pState &W, int index);
  //@}

  //@{ @name Primitive right eigenvector (x-direction).
  Euler2D_pState rp(int index);
  Euler2D_pState rp(int index) const;
  friend Euler2D_pState rp(const Euler2D_pState &W, int index);

  Euler2D_pState rp_x(int index);
  Euler2D_pState rp_x(int index) const;
  friend Euler2D_pState rp_x(const Euler2D_pState &W, int index);

  Euler2D_pState rp_precon_WS(int index);
  Euler2D_pState rp_precon_WS(int index) const;
  friend Euler2D_pState rp_precon_WS(const Euler2D_pState &W, int index);
  //@}

  //@{ @name Primitive right eigenvector (y-direction).
  Euler2D_pState rp_y(int index);
  Euler2D_pState rp_y(int index) const;
  friend Euler2D_pState rp_y(const Euler2D_pState &W, int index);
  //@}

  //@{ @name Conserved right eigenvector (x-direction).
  Euler2D_cState rc(int index);
  Euler2D_cState rc(int index) const;
  friend Euler2D_cState rc(const Euler2D_pState &W, int index);

  Euler2D_cState rc_x(int index);
  Euler2D_cState rc_x(int index) const;
  friend Euler2D_cState rc_x(const Euler2D_pState &W, int index);

  Euler2D_cState rc_precon_WS(int index);
  Euler2D_cState rc_precon_WS(int index) const;
  friend Euler2D_cState rc_precon_WS(const Euler2D_pState &W, int index);
  //@}

  //@{ @name Conserved right eigenvector (y-direction).
  Euler2D_cState rc_y(int index);
  Euler2D_cState rc_y(int index) const;
  friend Euler2D_cState rc_y(const Euler2D_pState &W, int index);
  //@}

  //@{ @name Primitive left eigenvector (x-direction).
  Euler2D_pState lp(int index);
  Euler2D_pState lp(int index) const;
  friend Euler2D_pState lp(const Euler2D_pState &W, int index);

  Euler2D_pState lp_x(int index);
  Euler2D_pState lp_x(int index) const;
  friend Euler2D_pState lp_x(const Euler2D_pState &W, int index);

  Euler2D_pState lp_precon_WS(int index);
  Euler2D_pState lp_precon_WS(int index) const;
  friend Euler2D_pState lp_precon_WS(const Euler2D_pState &W, int index);
  //@}

  //@{ @name Primitive left eigenvector (y-direction).
  Euler2D_pState lp_y(int index);
  Euler2D_pState lp_y(int index) const;
  friend Euler2D_pState lp_y(const Euler2D_pState &W, int index);
  //@}

  //@{ @name Source vector (axisymmetric terms).
  Euler2D_cState S(const Vector2D &X);
  Euler2D_cState S(const Vector2D &X) const;
  friend Euler2D_cState S(const Euler2D_pState &W, const Vector2D &X);

  void dSdU(DenseMatrix &dSdU,
	    const Vector2D &X,
	    const Euler2D_pState &dWdx,
	    const Euler2D_pState &dWdy,
	    const int &Axisymmetric) const {}
  //@}

  //@{ @name Solution variable Jacobian.
  void dUdW(DenseMatrix &dUdW);
  void dUdW(DenseMatrix &dUdW) const;
  friend void dUdW(DenseMatrix &dUdW, const Euler2D_pState &W);

  void dWdU(DenseMatrix &dWdU);
  void dWdU(DenseMatrix &dWdU) const;
  friend void dWdU(DenseMatrix &dWdU, const Euler2D_pState &W);
  //@}

  //@{ @name Weiss-Smith low-Mach-number local preconditioner.
  void P_U_WS(DenseMatrix &dUdW);
  void P_U_WS(DenseMatrix &dUdW) const;
  friend void P_U_WS(DenseMatrix &dUdW, const Euler2D_pState &W);

  void P_U_WS_inv(DenseMatrix &dWdU);
  void P_U_WS_inv(DenseMatrix &dWdU) const;
  friend void P_U_WS_inv(DenseMatrix &dWdU, const Euler2D_pState &W);
  //@}

  /* Assignment operator. */
  // Euler2D_pState operator = (const Euler2D_pState &W);
  // Use automatically generated assignment operator.

  //@{ @name Index operator.
  double &operator[](int index) {
    assert( index >= 1 && index <= NUM_VAR_EULER2D );
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
    assert( index >= 1 && index <= NUM_VAR_EULER2D );
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
  friend Euler2D_pState operator +(const Euler2D_pState &W1, const Euler2D_pState &W2);
  friend Euler2D_pState operator -(const Euler2D_pState &W1, const Euler2D_pState &W2);
  friend double operator *(const Euler2D_pState &W1, const Euler2D_pState &W2);
  friend Euler2D_pState operator *(const Euler2D_pState &W, const double &a);
  friend Euler2D_pState operator *(const double &a, const Euler2D_pState &W);
  friend Euler2D_pState operator /(const Euler2D_pState &W, const double &a);
  friend Euler2D_pState operator /(const Euler2D_pState &W1, const Euler2D_pState &W2);
  friend Euler2D_pState operator ^(const Euler2D_pState &W1, const Euler2D_pState &W2);
  friend Euler2D_pState max(const Euler2D_pState &W1, const Euler2D_pState &W2 );
  //@}

  //@{ @name Unary arithmetic operators.
  friend Euler2D_pState operator +(const Euler2D_pState &W);
  friend Euler2D_pState operator -(const Euler2D_pState &W);
  friend Euler2D_pState fabs(const Euler2D_pState &W);
  friend Euler2D_pState sqr(const Euler2D_pState &W);

  //@}

  //@{ @name Shortcut arithmetic operators.
  Euler2D_pState &operator +=(const Euler2D_pState &W);
  Euler2D_pState &operator -=(const Euler2D_pState &W);
  Euler2D_pState &operator /=(const Euler2D_pState &W);
  Euler2D_pState &operator *=(const Euler2D_pState &W);
  Euler2D_pState &operator *=(const double &a);
  Euler2D_pState &operator /=(const double &a);
  //@}

  //@{ @name Relational operators.
  friend int operator ==(const Euler2D_pState &W1, const Euler2D_pState &W2);
  friend int operator !=(const Euler2D_pState &W1, const Euler2D_pState &W2);
  friend bool operator >=(const Euler2D_pState& W1, const Euler2D_pState& W2);
  friend bool operator <=(const Euler2D_pState& W1, const Euler2D_pState& W2);
  friend bool operator <(const Euler2D_pState &W1, const Euler2D_pState &W2);
  friend bool operator >(const Euler2D_pState &W1, const Euler2D_pState &W2);
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const Euler2D_pState &W);
  friend istream &operator >> (istream &in_file,  Euler2D_pState &W);
  //@}

  //@{ @name Output functions.
  void output_labels(ostream &out_file) {
    out_file << "\"rho\" \\ \n"
	     << "\"u\" \\ \n"
	     << "\"v\" \\ \n"
	     << "\"p\" \\ \n"
	     << "\"T\" \\ \n"
	     << "\"M\" \\ \n"
	     << "\"H\" \\ \n"
	     << "\"s\" \\ \n";
  }

  void output_data(ostream &out_file, double dummy1, double dummy2) {  //dummies needed for compatibility with NS turbulent
    out_file << " " << d << " " << v.x << " " << v.y << " " << p
	     << " " << T() << " " << M() << " " << H() << " " << s();
  }
  //@}

  int analytically_inverted_relaxation() { //this is needed for embeddedboundaries with gaussian2D
    return 0;
  }
  void relax(double deltat, int stage, const Euler2D_pState &W) {return;} //this is needed for embeddedboundaries with gaussian2D

  double pressure() const {return p;} //added for compatibility with embeddedboundaries2D
};

/*!
 * Class: Euler2D_cState
 *
 * @brief Conserved variable solution state class definition for an
 *        inviscid compressible gas-flow.
 *
 * Conserved variable solution state class definition for an inviscid
 * compressible gas-flow.
 *
 * \verbatim
 * Member functions
 *     d        -- Return density.
 *     dv       -- Return momentum.
 *     E        -- Return total energy.
 *     g        -- Return specific heat ratio.
 *     gm1      -- Return g-1.
 *     gm1i     -- Return 1/(g-1).
 *     R        -- Return gas constant.
 *     setgas   -- Set gas constants.
 *     v        -- Return flow velocity.
 *     p        -- Return pressure.
 *     T        -- Return temperature.
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
 *     Mr_min   -- Reference Mach number used in Weiss-Smith low-Mach-
 *                 number local preconditioner.
 *     W        -- Return primitive solution state.
 *     F        -- Return x-direction solution flux.
 *     Fx       -- Return x-direction solution flux.
 *     Fy       -- Return y-direction solution flux.
 *     Fn       -- Return n-direction solution flux.
 *     dFdU     -- Return x-direction flux Jacobian.
 *     dFxdU    -- Return x-direction flux Jacobian.
 *     dFydU    -- Return y-direction flux Jacobian.
 *     dFndU    -- Return n-direction flux Jacobian.
 *     S        -- Return axisymmetric source term vector.
 *     dUdW     -- Return Jacobian of conserved solution variables wrt
 *                 primitive solution variables.
 *     dWdU     -- Return Jacobian of primitive solution variables wrt
 *                 conserved solution variables.
 *     P_U_WS   -- Weiss-Smith preconditioner for the conservative
 *                 solution variables.
 *     P_U_WS_inv -- Inverse of Weiss-Smith preconditioner for the
 *                   conservative solution variables.
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
class Euler2D_cState{
private:
public:
  //@{ @name Conserved variables and associated constants:
  double             d; //!< Density.
  Vector2D          dv; //!< Momentum.
  double             E; //!< Total Energy.
  static double      g; //!< Specific heat ratio.
  static double    gm1; //!< g-1
  static double   gm1i; //!< 1/(g-1)
  static double      R; //!< Gas constant.
  static double Mr_min; //!< Reference Mach number used in Weiss-Smith
                        //!< low-Mach-number local preconditioner.
  //@}
		      
  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  Euler2D_cState(void) {
    d = DENSITY_STDATM; dv.zero(); E = PRESSURE_STDATM/(GAMMA_AIR-ONE);
  }

  //! Copy constructor.
  Euler2D_cState(const Euler2D_cState &U) {
    d = U.d; dv = U.dv; E = U.E;
  }

  //! Copy constructor.
  Euler2D_cState(const Euler2D_pState &W);

  //! Assignment constructor.
  Euler2D_cState(const double &rho,
		 const Vector2D &rhoV,
		 const double &Etotal) {
    d = rho; dv = rhoV; E = Etotal;
  }

  //! Assignment constructor.
  Euler2D_cState(const double &rho,
		 const double &rhovx,
		 const double &rhovy,
		 const double &Etotal) {
    d = rho; dv.x = rhovx; dv.y = rhovy; E = Etotal;
  }

  /* Destructor. */
  // ~Euler2D_cState(void);
  // Use automatically generated destructor.
  //@}
  
  //@{ @name Useful operators.
  //! Return the number of variables.
  int NumVar(void) { return NUM_VAR_EULER2D; }

  //! Copy operator.
  void Copy(const Euler2D_cState &U) {
    d = U.d; dv = U.dv; E = U.E;
  }

  //! Vacuum operator.
  void Vacuum(void) {
    d = ZERO; dv = Vector2D_ZERO; E = ZERO;
  }

  //! Standard atmosphere operator.
  void Standard_Atmosphere(void) {
    d = DENSITY_STDATM; dv.zero(); E = PRESSURE_STDATM/(GAMMA_AIR-ONE);
  }

  //! Check for unphysical state properties.
  int Unphysical_Properties(void) const {
    if (d <= ZERO || E <= ZERO || e() <= ZERO) return 1;
    return 0;
  }

  //! Copy variables solved by multigrid only.
  void Copy_Multigrid_State_Variables(const Euler2D_cState &Ufine) {
    Copy(Ufine);
  }

  //! Zero variables not-solved by multigrid.
  void Zero_Non_Multigrid_State_Variables(void) { }

  //@{ @name State functions.
  //! Set gas constants.
  void setgas(void);
  void setgas(char *string_ptr);

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

  //! Reference Mach number used by Weiss-Smith low-Mach-number local
  //! preconditioner.
  double Mr(void);
  double Mr(void) const;

  //! Burning rate.
  double burningrate(void) const;
  //@}

  //@{ @name Primitive solution state.
  Euler2D_pState W(void);
  Euler2D_pState W(void) const;
  Euler2D_pState W(const Euler2D_cState &U);
  friend Euler2D_pState W(const Euler2D_cState &U);
  //@}

  //@{ @name Solution flux and Jacobian (x-direction).
  Euler2D_cState F(void);
  Euler2D_cState F(void) const;
  Euler2D_cState F(const Euler2D_cState &U);
  friend Euler2D_cState F(const Euler2D_cState &U);
  Euler2D_cState F(const Vector2D &V) const;
  void dFdU(DenseMatrix &dFdU);
  void dFdU(DenseMatrix &dFdU) const;
  friend void dFdU(DenseMatrix &dFdU, const Euler2D_cState &U);

  Euler2D_cState Fx(void);
  Euler2D_cState Fx(void) const;
  Euler2D_cState Fx(const Euler2D_cState &U);
  friend Euler2D_cState Fx(const Euler2D_cState &U);
  void dFxdU(DenseMatrix &dFxdU);
  void dFxdU(DenseMatrix &dFxdU) const;
  friend void dFxdU(DenseMatrix &dFxdU, const Euler2D_cState &U);
  //@}

  //@{ @name Solution flux and Jacobian (y-direction).
  Euler2D_cState Fy(void);
  Euler2D_cState Fy(void) const;
  Euler2D_cState Fy(const Euler2D_cState &U);
  friend Euler2D_cState Fy(const Euler2D_cState &U);
  void dFydU(DenseMatrix &dFydU);
  void dFydU(DenseMatrix &dFydU) const;
  friend void dFydU(DenseMatrix &dFydU, const Euler2D_cState &U);
  //@}

  //@{ @name Solution flux and Jacobian (n-direction).
  Euler2D_cState Fn(void);
  Euler2D_cState Fn(void) const;
  Euler2D_cState Fn(const Euler2D_cState &U);
  friend Euler2D_cState Fn(const Euler2D_cState &U);
  void dFndU(DenseMatrix &dFndU);
  void dFndU(DenseMatrix &dFndU) const;
  friend void dFndU(DenseMatrix &dFndU, const Euler2D_cState &U);
  //@}

  //@{ @name Source vector (axisymmetric terms).
  Euler2D_cState S(const Vector2D &X);
  Euler2D_cState S(const Vector2D &X) const;
  friend Euler2D_cState S(const Euler2D_cState &U, const Vector2D &X);
  //@}

  //@{ @name Solution variable Jacobian.
  void dUdW(DenseMatrix &dUdW);
  void dUdW(DenseMatrix &dUdW) const;
  friend void dUdW(DenseMatrix &dUdW, const Euler2D_cState &U);

  void dWdU(DenseMatrix &dWdU);
  void dWdU(DenseMatrix &dWdU) const;
  friend void dWdU(DenseMatrix &dWdU, const Euler2D_cState &U);
  //@}

  //@{ @name Weiss-Smith low-Mach-number local preconditioner.
  void P_U_WS(DenseMatrix &dUdW);
  void P_U_WS(DenseMatrix &dUdW) const;
  friend void P_U_WS(DenseMatrix &dUdW, const Euler2D_cState &W);

  void P_U_WS_inv(DenseMatrix &dWdU);
  void P_U_WS_inv(DenseMatrix &dWdU) const;
  friend void P_U_WS_inv(DenseMatrix &dWdU, const Euler2D_cState &W);
  //@}

  /* Assignment operator. */
  // Euler2D_cState operator = (const Euler2D_cState &U);
  // Use automatically generated assignment operator.

  //@{ @name Index operator.
  double &operator[](int index) {
    assert( index >= 1 && index <= NUM_VAR_EULER2D );
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
    assert( index >= 1 && index <= NUM_VAR_EULER2D );
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
  friend Euler2D_cState operator +(const Euler2D_cState &U1, const Euler2D_cState &U2);
  friend Euler2D_cState operator -(const Euler2D_cState &U1, const Euler2D_cState &U2);
  friend double operator *(const Euler2D_cState &U1, const Euler2D_cState &U2);
  friend Euler2D_cState operator *(const Euler2D_cState &U, const double &a);
  friend Euler2D_cState operator *(const double &a, const Euler2D_cState &U);
  friend Euler2D_cState operator /(const Euler2D_cState &U, const double &a);
  friend Euler2D_cState operator ^(const Euler2D_cState &U1, const Euler2D_cState &U2);
  //@}

  //@{ @name Unary arithmetic operators.
  friend Euler2D_cState operator +(const Euler2D_cState &U);
  friend Euler2D_cState operator -(const Euler2D_cState &U);
  //@}

  //@{ @name Shortcut arithmetic operators.
  Euler2D_cState &operator +=(const Euler2D_cState &U);
  Euler2D_cState &operator -=(const Euler2D_cState &U);
  Euler2D_cState &operator *=(const double &a);
  Euler2D_cState &operator /=(const double &a);
  //@}

  //@{ @name Relational operators.
  friend int operator ==(const Euler2D_cState &U1, const Euler2D_cState &U2);
  friend int operator !=(const Euler2D_cState &U1, const Euler2D_cState &U2);
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const Euler2D_cState &U);
  friend istream &operator >> (istream &in_file,  Euler2D_cState &U);
  //@}

};

/********************************************
 * Euler2D_pState Value Constructor.        *
 *******************************************/
inline Euler2D_pState::Euler2D_pState(const double &Val):
  d(Val), v(Val), p(Val){
}

/********************************************************
 * Euler2D_pState::setgas -- Assign gas constants.      *
 ********************************************************/
inline void Euler2D_pState::setgas(void) {
    g = GAMMA_AIR;
    R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    gm1 = g - ONE;
    gm1i = ONE/gm1;
}

inline void Euler2D_pState::setgas(char *string_ptr) {
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
   } else if (strcmp(string_ptr, "AP-HTPB") == 0) {
     g = GAMMA_APHTPB;
     R = R_UNIVERSAL/(MOLE_WT_APHTPB*MILLI);
   } else {
     g = GAMMA_AIR;
     R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
   } /* endif */
   gm1 = g - ONE;
   gm1i = ONE/gm1;
}

/********************************************************
 * Euler2D_pState::T -- Temperature.                    *
 ********************************************************/
inline double Euler2D_pState::T(void) {
    return (p/(d*R));
}

inline double Euler2D_pState::T(void) const {
    return (p/(d*R));
}

/********************************************************
 * Euler2D_pState::e -- Specific internal energy.       *
 ********************************************************/
inline double Euler2D_pState::e(void) {
    return (p/(gm1*d));
}

inline double Euler2D_pState::e(void) const {
    return (p/(gm1*d));
}

/********************************************************
 * Euler2D_pState::E -- Total energy.                   *
 ********************************************************/
inline double Euler2D_pState::E(void) {
    return (p*gm1i + HALF*d*v.sqr());
}

inline double Euler2D_pState::E(void) const {
    return (p*gm1i + HALF*d*v.sqr());
}

/********************************************************
 * Euler2D_pState::h -- Specific enthalpy.              *
 ********************************************************/
inline double Euler2D_pState::h(void) {
    return (g*p/(gm1*d) + HALF*v.sqr());
}

inline double Euler2D_pState::h(void) const {
    return (g*p/(gm1*d) + HALF*v.sqr());
}

/********************************************************
 * Euler2D_pState::H -- Total enthalpy.                 *
 ********************************************************/
inline double Euler2D_pState::H(void) {
    return (g*gm1i*p + HALF*d*v.sqr());
}

inline double Euler2D_pState::H(void) const {
    return (g*gm1i*p + HALF*d*v.sqr());
}

/********************************************************
 * Euler2D_pState::a -- Sound speed.                    *
 ********************************************************/
inline double Euler2D_pState::a(void) {
    return (sqrt(g*p/d));
}

inline double Euler2D_pState::a(void) const {
    return (sqrt(g*p/d));
}

/********************************************************
 * Euler2D_pState::a2 -- Sound speed squared.           *
 ********************************************************/
inline double Euler2D_pState::a2(void) {
    return (g*p/d);
}

inline double Euler2D_pState::a2(void) const {
    return (g*p/d);
}

/********************************************************
 * Euler2D_pState::M -- Mach number.                    *
 ********************************************************/
inline double Euler2D_pState::M(void) {
    return (abs(v)/sqrt(g*p/d));
}

inline double Euler2D_pState::M(void) const {
    return (abs(v)/sqrt(g*p/d));
}

/********************************************************
 * Euler2D_pState::s -- Specific entropy.               *
 ********************************************************/
inline double Euler2D_pState::s(void) {
    return (R*gm1i*log(p/pow(d, g)));
}

inline double Euler2D_pState::s(void) const {
    return (R*gm1i*log(p/pow(d, g)));
}

/********************************************************
 * Euler2D_pState::dv -- Momentum.                      *
 ********************************************************/
inline Vector2D Euler2D_pState::dv(void) {
    return (d*v);
}

inline Vector2D Euler2D_pState::dv(void) const {
    return (d*v);
}

inline double Euler2D_pState::dv(const Vector2D &n) {
    return (d*(v*n));
}

inline double Euler2D_pState::dv(const Vector2D &n) const {
    return (d*(v*n));
}

/********************************************************
 * Euler2D_pState::To -- Stagnation temperature.        *
 ********************************************************/
inline double Euler2D_pState::To(void) {
    return ((p/(d*R))*(ONE+HALF*gm1*v.sqr()/(g*p/d)));
}

inline double Euler2D_pState::To(void) const {
    return ((p/(d*R))*(ONE+HALF*gm1*v.sqr()/(g*p/d)));
}

/********************************************************
 * Euler2D_pState::po -- Stagnation pressure.           *
 ********************************************************/
inline double Euler2D_pState::po(void) {
    return (p*pow(ONE+HALF*gm1*v.sqr()/(g*p/d), g*gm1i));
}

inline double Euler2D_pState::po(void) const {
    return (p*pow(ONE+HALF*gm1*v.sqr()/(g*p/d), g*gm1i));
}

/********************************************************
 * Euler2D_pState::ao -- Stagnation sound speed.        *
 ********************************************************/
inline double Euler2D_pState::ao(void) {
    return (sqrt((g*p/d)*(ONE+HALF*gm1*v.sqr()/(g*p/d))));
}

inline double Euler2D_pState::ao(void) const {
    return (sqrt((g*p/d)*(ONE+HALF*gm1*v.sqr()/(g*p/d))));
}

/********************************************************
 * Euler2D_pState::ho -- Stagnation enthalpy.           *
 ********************************************************/
inline double Euler2D_pState::ho(void) {
    return ((g*p/(gm1*d) + HALF*v.sqr())*(ONE+HALF*gm1*v.sqr()/(g*p/d)));
}

inline double Euler2D_pState::ho(void) const {
    return ((g*p/(gm1*d) + HALF*v.sqr())*(ONE+HALF*gm1*v.sqr()/(g*p/d)));
}

/********************************************************
 * Euler2D_pState::Mr -- Reference Mach number.         *
 ********************************************************/
inline double Euler2D_pState::Mr(void) {
    return (min(ONE, max(Mr_min, abs(v)/sqrt(g*p/d))));
}

inline double Euler2D_pState::Mr(void) const {
    return (min(ONE, max(Mr_min, abs(v)/sqrt(g*p/d))));
}

/**********************************************************************
 * Euler2D_pState::burningrate -- Solid propellent burning rate.      *
 **********************************************************************/
inline double Euler2D_pState::burningrate(void) const {
  return -BETA_APHTPB*pow(p,N_APHTPB);
}

/********************************************************
 * Euler2D_pState -- Binary arithmetic operators.       *
 ********************************************************/
inline Euler2D_pState operator +(const Euler2D_pState &W1, const Euler2D_pState &W2) {
  return (Euler2D_pState(W1.d+W2.d,W1.v+W2.v,W1.p+W2.p));
}

inline Euler2D_pState operator -(const Euler2D_pState &W1, const Euler2D_pState &W2) {
  return (Euler2D_pState(W1.d-W2.d,W1.v-W2.v,W1.p-W2.p));
}

// Inner product operator.
inline double operator *(const Euler2D_pState &W1, const Euler2D_pState &W2) {
   return (W1.d*W2.d+W1.v*W2.v+W1.p*W2.p);
}

inline Euler2D_pState operator *(const Euler2D_pState &W, const double &a) {
  return (Euler2D_pState(a*W.d,a*W.v,a*W.p));
}

inline Euler2D_pState operator *(const double &a, const Euler2D_pState &W) {
  return (Euler2D_pState(a*W.d,a*W.v,a*W.p));
}

inline Euler2D_pState operator /(const Euler2D_pState &W, const double &a) {
  return (Euler2D_pState(W.d/a,W.v/a,W.p/a));
}

inline Euler2D_pState operator /(const Euler2D_pState &W1, const Euler2D_pState &W2) {
  return (Euler2D_pState(W1[1]/W2[1], W1[2]/W2[2], W1[3]/W2[3], W1[4]/W2[4]));
}

// My useful solution state product operator.
inline Euler2D_pState operator ^(const Euler2D_pState &W1, const Euler2D_pState &W2) {
   return (Euler2D_pState(W1.d*W2.d,W1.v.x*W2.v.x,W1.v.y*W2.v.y,W1.p*W2.p));
}

/*!
 * Compute maximum between 2 states. 
 * Return the state of maximum values.
 */
inline Euler2D_pState max(const Euler2D_pState &W1, const Euler2D_pState &W2 ){
  return Euler2D_pState(max(W1.d,W2.d),max(W1.v.x,W2.v.x),max(W1.v.y,W2.v.y),max(W1.p,W2.p));
}

/********************************************************
 * Euler2D_pState -- Unary arithmetic operators.        *
 ********************************************************/
inline Euler2D_pState operator +(const Euler2D_pState &W) {
  return (Euler2D_pState(W.d,W.v,W.p));
}

inline Euler2D_pState operator -(const Euler2D_pState &W) {
  return (Euler2D_pState(-W.d,-W.v,-W.p));
}

inline Euler2D_pState fabs(const Euler2D_pState &W){
  return Euler2D_pState(fabs(W[1]),fabs(W[2]), fabs(W[3]),fabs(W[4]));
}

inline Euler2D_pState sqr(const Euler2D_pState &W){
  return Euler2D_pState(sqr(W.d),sqr(W.v.x),sqr(W.v.y),sqr(W.p));
}


/********************************************************
 * Euler2D_pState -- Shortcut arithmetic operators.     *
 ********************************************************/
inline Euler2D_pState& Euler2D_pState::operator +=(const Euler2D_pState &W) {
  d += W.d; v.x += W.v.x; v.y += W.v.y; p += W.p;
  return *this;
}

inline Euler2D_pState& Euler2D_pState::operator -=(const Euler2D_pState &W) {
  d -= W.d; v.x -= W.v.x; v.y -= W.v.y; p -= W.p;
  return *this;
}

inline Euler2D_pState& Euler2D_pState::operator /=(const Euler2D_pState &W){
  d /= W.d;
  v.x /= W.v.x;
  v.y /= W.v.y;
  p /= W.p;
  return *this;
}

inline Euler2D_pState& Euler2D_pState::operator *=(const Euler2D_pState &W) {
  d *= W.d; v.x *= W.v.x; v.y *= W.v.y; p *= W.p;
  return *this;
}

inline Euler2D_pState& Euler2D_pState::operator *=(const double &a) {
  d *= a; v.x *= a; v.y *= a; p *= a;
  return *this;
}

inline Euler2D_pState& Euler2D_pState::operator /=(const double &a) {
  d /= a; v.x /= a; v.y /= a; p /= a;
  return *this;
}

/********************************************************
 * Euler2D_pState -- Relational operators.              *
 ********************************************************/
inline int operator ==(const Euler2D_pState &W1, const Euler2D_pState &W2) {
  return (W1.d == W2.d && W1.v == W2.v && W1.p == W2.p);
}

inline int operator !=(const Euler2D_pState &W1, const Euler2D_pState &W2) {
  return (W1.d != W2.d || W1.v != W2.v || W1.p != W2.p);
}

inline bool operator <=(const Euler2D_pState &W1, const Euler2D_pState &W2) {
  return (W1[1]<=W2[1] && W1[2]<=W2[2] && W1[3]<=W2[3] && W1[4]<=W2[4] );
}

inline bool operator >=(const Euler2D_pState &W1, const Euler2D_pState &W2) {
  return (W1[1]>=W2[1] && W1[2]>=W2[2] && W1[3]>=W2[3] && W1[4]>=W2[4] );
}

inline bool operator <(const Euler2D_pState &W1, const Euler2D_pState &W2) {
  return (W1[1]<W2[1] && W1[2]<W2[2] && W1[3]<W2[3] && W1[4]<W2[4] );
}

inline bool operator >(const Euler2D_pState &W1, const Euler2D_pState &W2) {
  return (W1[1]>W2[1] && W1[2]>W2[2] && W1[3]>W2[3] && W1[4]>W2[4] );
}

/********************************************************
 * Euler2D_pState -- Input-output operators.            *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const Euler2D_pState &W) {
  out_file.setf(ios::scientific);
  out_file << " " << W.d  << " " << W.v.x << " " << W.v.y << " " << W.p;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, Euler2D_pState &W) {
  in_file.setf(ios::skipws);
  in_file >> W.d >> W.v.x >> W.v.y >> W.p;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/********************************************************
 * Euler2D_cState::setgas -- Assign gas constants.      *
 ********************************************************/
inline void Euler2D_cState::setgas(void) {
    g = GAMMA_AIR;
    R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    gm1 = g - ONE;
    gm1i = ONE/gm1;
}

inline void Euler2D_cState::setgas(char *string_ptr) {
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
   } else if (strcmp(string_ptr, "AP-HTPB") == 0) {
     g = GAMMA_APHTPB;
     R = R_UNIVERSAL/(MOLE_WT_APHTPB*MILLI);
   } else {
     g = GAMMA_AIR;
     R = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
   } /* endif */
   gm1 = g - ONE;
   gm1i = ONE/gm1;
}

/********************************************************
 * Euler2D_cState::v -- Flow velocity.                  *
 ********************************************************/
inline Vector2D Euler2D_cState::v(void) {
    return (dv/d);
}

inline Vector2D Euler2D_cState::v(void) const {
    return (dv/d);
}

inline double Euler2D_cState::v(const Vector2D &n) {
    return ((dv*n)/d);
}

inline double Euler2D_cState::v(const Vector2D &n) const {
    return ((dv*n)/d);
}

/********************************************************
 * Euler2D_cState::p -- Pressure.                       *
 ********************************************************/
inline double Euler2D_cState::p(void) {
    return (gm1*(E - HALF*dv.sqr()/d));
}

inline double Euler2D_cState::p(void) const {
    return (gm1*(E - HALF*dv.sqr()/d));
}

/********************************************************
 * Euler2D_cState::T -- Temperature.                    *
 ********************************************************/
inline double Euler2D_cState::T(void) {
    return (gm1*(E - HALF*dv.sqr()/d)/(d*R));
}

inline double Euler2D_cState::T(void) const {
    return (gm1*(E - HALF*dv.sqr()/d)/(d*R));
}

/********************************************************
 * Euler2D_cState::e -- Specific internal energy.       *
 ********************************************************/
inline double Euler2D_cState::e(void) {
    return (E/d - HALF*dv.sqr()/sqr(d));
}

inline double Euler2D_cState::e(void) const {
    return (E/d - HALF*dv.sqr()/sqr(d));
}

/********************************************************
 * Euler2D_cState::h -- Specific enthalpy.              *
 ********************************************************/
inline double Euler2D_cState::h(void) {
    return (g*E/d - gm1*HALF*dv.sqr()/sqr(d));
}

inline double Euler2D_cState::h(void) const {
    return (g*E/d - gm1*HALF*dv.sqr()/sqr(d));
}

/********************************************************
 * Euler2D_cState::H -- Total enthalpy.                 *
 ********************************************************/
inline double Euler2D_cState::H(void) {
     return (g*E - gm1*HALF*dv.sqr()/d);
}

inline double Euler2D_cState::H(void) const {
     return (g*E - gm1*HALF*dv.sqr()/d);
}

/********************************************************
 * Euler2D_cState::a -- Sound speed.                    *
 ********************************************************/
inline double Euler2D_cState::a(void) {
    return (sqrt(g*gm1*(E/d - HALF*dv.sqr()/sqr(d))));
}

inline double Euler2D_cState::a(void) const {
    return (sqrt(g*gm1*(E/d - HALF*dv.sqr()/sqr(d))));
}

/********************************************************
 * Euler2D_cState::a2 -- Sound speed squared.           *
 ********************************************************/
inline double Euler2D_cState::a2(void) {
    return (g*gm1*(E/d - HALF*dv.sqr()/sqr(d)));
}

inline double Euler2D_cState::a2(void) const {
    return (g*gm1*(E/d - HALF*dv.sqr()/sqr(d)));
}

/********************************************************
 * Euler2D_cState::M -- Mach number.                    *
 ********************************************************/
inline double Euler2D_cState::M(void) {
    return (abs(dv)/(d*sqrt(g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))));
}

inline double Euler2D_cState::M(void) const {
    return (abs(dv)/(d*sqrt(g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))));
}

/********************************************************
 * Euler2D_cState::s -- Specific entropy.               *
 ********************************************************/
inline double Euler2D_cState::s(void) {
    return (R*gm1i*log(gm1*(E - HALF*dv.sqr()/d)/pow(d, g)));
}

inline double Euler2D_cState::s(void) const {
    return (R*gm1i*log(gm1*(E - HALF*dv.sqr()/d)/pow(d, g)));
}

/********************************************************
 * Euler2D_cState::To -- Stagnation temperature.        *
 ********************************************************/
inline double Euler2D_cState::To(void) {
    return ((gm1*(E - HALF*dv.sqr()/d)/(d*R))*
	    (ONE+HALF*gm1*dv.sqr()/(d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))));
}

inline double Euler2D_cState::To(void) const {
    return ((gm1*(E - HALF*dv.sqr()/d)/(d*R))*
	    (ONE+HALF*gm1*dv.sqr()/(d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))));
}

/********************************************************
 * Euler2D_cState::po -- Stagnation pressure.           *
 ********************************************************/
inline double Euler2D_cState::po(void) {
    return ((gm1*(E - HALF*dv.sqr()/d))*
	    pow(ONE+HALF*gm1*dv.sqr()/(d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d))), g*gm1i));
}

inline double Euler2D_cState::po(void) const {
    return ((gm1*(E - HALF*dv.sqr()/d))*
	    pow(ONE+HALF*gm1*dv.sqr()/(d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d))), g*gm1i));
}

/********************************************************
 * Euler2D_cState::ao -- Stagnation sound speed.        *
 ********************************************************/
inline double Euler2D_cState::ao(void) {
    return (sqrt((g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))*
	         (ONE+HALF*gm1*dv.sqr()/(d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d))))));
}

inline double Euler2D_cState::ao(void) const {
    return (sqrt((g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))*
	         (ONE+HALF*gm1*dv.sqr()/(d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d))))));
}

/********************************************************
 * Euler2D_cState::ho -- Stagnation enthalpy.           *
 ********************************************************/
inline double Euler2D_cState::ho(void) {
    return ((g*E/d - gm1*HALF*dv.sqr()/sqr(d))*
	    (ONE+HALF*gm1*dv.sqr()/(d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))));
}

inline double Euler2D_cState::ho(void) const {
    return ((g*E/d - gm1*HALF*dv.sqr()/sqr(d))*
	    (ONE+HALF*gm1*dv.sqr()/(d*d*g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))));
}

/********************************************************
 * Euler2D_cState::Mr -- Reference Mach number.         *
 ********************************************************/
inline double Euler2D_cState::Mr(void) {
    return (min(ONE, max(Mr_min, abs(dv)/(d*sqrt(g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))) )));
}

inline double Euler2D_cState::Mr(void) const {
    return (min(ONE, max(Mr_min, abs(dv)/(d*sqrt(g*gm1*(E/d - HALF*dv.sqr()/sqr(d)))) )));
}

/**********************************************************************
 * Euler2D_cState::burningrate -- Solid propellent burning rate.      *
 **********************************************************************/
inline double Euler2D_cState::burningrate(void) const {
  return -BETA_APHTPB*pow(p(),N_APHTPB);
}

/********************************************************
 * Euler2D_cState -- Binary arithmetic operators.       *
 ********************************************************/
inline Euler2D_cState operator +(const Euler2D_cState &U1, const Euler2D_cState &U2) {
  return (Euler2D_cState(U1.d+U2.d,U1.dv+U2.dv,U1.E+U2.E));
}

inline Euler2D_cState operator -(const Euler2D_cState &U1, const Euler2D_cState &U2) {
  return (Euler2D_cState(U1.d-U2.d,U1.dv-U2.dv,U1.E-U2.E));
}

// Inner product operator.
inline double operator *(const Euler2D_cState &U1, const Euler2D_cState &U2) {
   return (U1.d*U2.d+U1.dv*U2.dv+U1.E*U2.E);
}

inline Euler2D_cState operator *(const Euler2D_cState &U, const double &a) {
  return (Euler2D_cState(a*U.d,a*U.dv,a*U.E));
}

inline Euler2D_cState operator *(const double &a, const Euler2D_cState &U) {
  return (Euler2D_cState(a*U.d,a*U.dv,a*U.E));
}

inline Euler2D_cState operator /(const Euler2D_cState &U, const double &a) {
  return (Euler2D_cState(U.d/a,U.dv/a,U.E/a));
}

// My useful solution state product operator.
inline Euler2D_cState operator ^(const Euler2D_cState &U1, const Euler2D_cState &U2) {
   return (Euler2D_cState(U1.d*U2.d,U1.dv.x*U2.dv.x,U1.dv.y*U2.dv.y,U1.E*U2.E));
}

/********************************************************
 * Euler2D_cState -- Unary arithmetic operators.        *
 ********************************************************/
inline Euler2D_cState operator +(const Euler2D_cState &U) {
  return (Euler2D_cState(U.d,U.dv,U.E));
}

inline Euler2D_cState operator -(const Euler2D_cState &U) {
  return (Euler2D_cState(-U.d,-U.dv,-U.E));
}

/********************************************************
 * Euler2D_cState -- Shortcut arithmetic operators.     *
 ********************************************************/
inline Euler2D_cState& Euler2D_cState::operator +=(const Euler2D_cState &U) {
  d += U.d; dv.x += U.dv.x; dv.y += U.dv.y; E += U.E;
  return *this;
}

inline Euler2D_cState& Euler2D_cState::operator -=(const Euler2D_cState &U) {
  d -= U.d; dv.x -= U.dv.x; dv.y -= U.dv.y; E -= U.E;
  return *this;
}

inline Euler2D_cState& Euler2D_cState::operator *=(const double &a) {
  d *= a; dv.x *= a; dv.y *= a; E *= a;
  return *this;
}

inline Euler2D_cState& Euler2D_cState::operator /=(const double &a) {
  d /= a; dv.x /= a; dv.y /= a; E /= a;
  return *this;
}

/********************************************************
 * Euler2D_cState -- Relational operators.              *
 ********************************************************/
inline int operator ==(const Euler2D_cState &U1, const Euler2D_cState &U2) {
  return (U1.d == U2.d && U1.dv == U2.dv && U1.E == U2.E);
}

inline int operator !=(const Euler2D_cState &U1, const Euler2D_cState &U2) {
  return (U1.d != U2.d || U1.dv != U2.dv || U1.E != U2.E);
}

/********************************************************
 * Euler2D_cState -- Input-output operators.            *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const Euler2D_cState &U) {
  out_file.setf(ios::scientific);
  out_file << " " << U.d  << " " << U.dv.x << " " << U.dv.y << " " << U.E;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, Euler2D_cState &U) {
  in_file.setf(ios::skipws);
  in_file >> U.d >> U.dv.x >> U.dv.y >> U.E;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/********************************************************
 * Euler2D_pState::Euler2D_pState -- Constructor.       *
 ********************************************************/
inline Euler2D_pState::Euler2D_pState(const Euler2D_cState &U) {
  d = U.d; v = U.v(); p = U.p();
}

/********************************************************
 * Euler2D_pState::U -- Conserved solution state.       *
 ********************************************************/
inline Euler2D_cState Euler2D_pState::U(void) {
  return (Euler2D_cState(d, dv(), E()));
}

inline Euler2D_cState Euler2D_pState::U(void) const {
  return (Euler2D_cState(d, dv(), E()));
}

inline Euler2D_cState Euler2D_pState::U(const Euler2D_pState &W) {
  return (Euler2D_cState(W.d, W.dv(), W.E()));
}

inline Euler2D_cState U(const Euler2D_pState &W) {
  return (Euler2D_cState(W.d, W.dv(), W.E()));
}

/********************************************************
 * Euler2D_pState::F -- Solution flux (x-direction).    *
 ********************************************************/
inline Euler2D_cState Euler2D_pState::F(void) {
  return (Euler2D_cState(d*v.x, d*sqr(v.x) + p, d*v.x*v.y, v.x*H()));
}

inline Euler2D_cState Euler2D_pState::F(void) const {
  return (Euler2D_cState(d*v.x, d*sqr(v.x) + p, d*v.x*v.y, v.x*H()));
}

inline Euler2D_cState Euler2D_pState::F(const Euler2D_pState &W) {
  return (Euler2D_cState(W.d*W.v.x, W.d*sqr(W.v.x) + W.p,
                         W.d*W.v.x*W.v.y, W.v.x*W.H()));
}

inline Euler2D_cState F(const Euler2D_pState &W) {
  return (Euler2D_cState(W.d*W.v.x, W.d*sqr(W.v.x) + W.p,
                         W.d*W.v.x*W.v.y, W.v.x*W.H()));
}

inline Euler2D_cState Euler2D_pState::F(const Vector2D &V) const {
  return Euler2D_cState(d*(v.x-V.x),
			d*(v.x-V.x)*v.x + p,
			d*(v.x-V.x)*v.y,
			(v.x-V.x)*E() + v.x*p);
}

inline void Euler2D_pState::dFdU(DenseMatrix &dFdU) {
  dFdU(0,1) += ONE;
  dFdU(1,0) += HALF*(sqr(v.x)*(g-THREE)+sqr(v.y)*gm1);
  dFdU(1,1) -= v.x*(g-THREE);
  dFdU(1,2) -= v.y*gm1;
  dFdU(1,3) += gm1;
  dFdU(2,0) -= v.x*v.y;
  dFdU(2,1) += v.y;
  dFdU(2,2) += v.x;
  dFdU(3,0) += HALF*gm1i*v.x*(d*(TWO+g*g-THREE*g)*(sqr(v.x)+sqr(v.y))-TWO*p*g)/d;
  dFdU(3,1) += HALF*gm1i*
    (sqr(v.x)*d*(FIVE*g-THREE-TWO*g*g)+sqr(v.y)*d*gm1+TWO*p*g)/d;
  dFdU(3,2) -= v.x*v.y*gm1;
  dFdU(3,3) += g*v.x;
}

inline void Euler2D_pState::dFdU(DenseMatrix &dFdU) const {
  dFdU(0,1) += ONE;
  dFdU(1,0) += HALF*(sqr(v.x)*(g-THREE)+sqr(v.y)*gm1);
  dFdU(1,1) -= v.x*(g-THREE);
  dFdU(1,2) -= v.y*gm1;
  dFdU(1,3) += gm1;
  dFdU(2,0) -= v.x*v.y;
  dFdU(2,1) += v.y;
  dFdU(2,2) += v.x;
  dFdU(3,0) += HALF*gm1i*v.x*(d*(TWO+g*g-THREE*g)*(sqr(v.x)+sqr(v.y))-TWO*p*g)/d;
  dFdU(3,1) += HALF*gm1i*
    (sqr(v.x)*d*(FIVE*g-THREE-TWO*g*g)+sqr(v.y)*d*gm1+TWO*p*g)/d;
  dFdU(3,2) -= v.x*v.y*gm1;
  dFdU(3,3) += g*v.x;
}

inline void dFdU(DenseMatrix &dFdU, const Euler2D_pState &W) {
  dFdU(0,1) += ONE;
  dFdU(1,0) += HALF*(sqr(W.v.x)*(W.g-THREE)+sqr(W.v.y)*(W.g-1));
  dFdU(1,1) -= W.v.x*(W.g-THREE);
  dFdU(1,2) -= W.v.y*W.gm1;
  dFdU(1,3) += W.gm1;
  dFdU(2,0) -= W.v.x*W.v.y;
  dFdU(2,1) += W.v.y;
  dFdU(2,2) += W.v.x;
  dFdU(3,0) += HALF*W.gm1i*W.v.x*
    (W.d*(TWO+W.g*W.g-THREE*W.g)*(sqr(W.v.x)+sqr(W.v.y))-TWO*W.p*W.g)/W.d;
  dFdU(3,1) += HALF*W.gm1i*
    (sqr(W.v.x)*W.d*(FIVE*W.g-THREE-TWO*W.g*W.g)
     +sqr(W.v.y)*W.d*W.gm1+TWO*W.p*W.g)/W.d;
  dFdU(3,2) -= W.v.x*W.v.y*W.gm1;
  dFdU(3,3) += W.g*W.v.x;
}

/********************************************************
 * Euler2D_pState::Fx -- Solution flux (x-direction).   *
 ********************************************************/
inline Euler2D_cState Euler2D_pState::Fx(void) {
  return (Euler2D_cState(d*v.x, d*sqr(v.x) + p, d*v.x*v.y, v.x*H()));
}

inline Euler2D_cState Euler2D_pState::Fx(void) const {
  return (Euler2D_cState(d*v.x, d*sqr(v.x) + p, d*v.x*v.y, v.x*H()));
}

inline Euler2D_cState Euler2D_pState::Fx(const Euler2D_pState &W) {
  return (Euler2D_cState(W.d*W.v.x, W.d*sqr(W.v.x) + W.p,
                         W.d*W.v.x*W.v.y, W.v.x*W.H()));
}

inline Euler2D_cState Fx(const Euler2D_pState &W) {
  return (Euler2D_cState(W.d*W.v.x, W.d*sqr(W.v.x) + W.p,
                         W.d*W.v.x*W.v.y, W.v.x*W.H()));
}

inline void Euler2D_pState::dFxdU(DenseMatrix &dFxdU) {
  dFxdU(0,1) += ONE;
  dFxdU(1,0) += HALF*(sqr(v.x)*(g-THREE)+sqr(v.y)*gm1);
  dFxdU(1,1) -= v.x*(g-THREE);
  dFxdU(1,2) -= v.y*gm1;
  dFxdU(1,3) += gm1;
  dFxdU(2,0) -= v.x*v.y;
  dFxdU(2,1) += v.y;
  dFxdU(2,2) += v.x;
  dFxdU(3,0) -= HALF*gm1i*v.x*(d*(THREE*g-TWO-g*g)*(sqr(v.x)+sqr(v.y))+TWO*p*g)/d;
  dFxdU(3,1) += HALF*gm1i*
    (sqr(v.x)*d*(FIVE*g-THREE-TWO*g*g)+sqr(v.y)*d*gm1+TWO*p*g)/d;
  dFxdU(3,2) -= v.x*v.y*gm1;
  dFxdU(3,3) += g*v.x;
}

inline void Euler2D_pState::dFxdU(DenseMatrix &dFxdU) const {
  dFxdU(0,1) += ONE;
  dFxdU(1,0) += HALF*(sqr(v.x)*(g-THREE)+sqr(v.y)*gm1);
  dFxdU(1,1) -= v.x*(g-THREE);
  dFxdU(1,2) -= v.y*gm1;
  dFxdU(1,3) += gm1;
  dFxdU(2,0) -= v.x*v.y;
  dFxdU(2,1) += v.y;
  dFxdU(2,2) += v.x;
  dFxdU(3,0) -= HALF*gm1i*v.x*(d*(THREE*g-TWO-g*g)*(sqr(v.x)+sqr(v.y))+TWO*p*g)/d;
  dFxdU(3,1) += HALF*gm1i*
    (sqr(v.x)*d*(FIVE*g-THREE-TWO*g*g)+sqr(v.y)*d*gm1+TWO*p*g)/d;
  dFxdU(3,2) -= v.x*v.y*gm1;
  dFxdU(3,3) += g*v.x;
}

inline void dFxdU(DenseMatrix &dFxdU, const Euler2D_pState &W) {
  dFxdU(0,1) += ONE;
  dFxdU(1,0) += HALF*(sqr(W.v.x)*(W.g-THREE)+sqr(W.v.y)*(W.g-1));
  dFxdU(1,1) -= W.v.x*(W.g-THREE);
  dFxdU(1,2) -= W.v.y*W.gm1;
  dFxdU(1,3) += W.gm1;
  dFxdU(2,0) -= W.v.x*W.v.y;
  dFxdU(2,1) += W.v.y;
  dFxdU(2,2) += W.v.x;
  dFxdU(3,0) -= HALF*W.gm1i*W.v.x*(W.d*(THREE*W.g-TWO-W.g*W.g)*(sqr(W.v.x)+sqr(W.v.y))+TWO*W.p*W.g)/W.d;
  dFxdU(3,1) += HALF*W.gm1i*
    (sqr(W.v.x)*W.d*(FIVE*W.g-THREE-TWO*W.g*W.g)
     +sqr(W.v.y)*W.d*W.gm1+TWO*W.p*W.g)/W.d;
  dFxdU(3,2) -= W.v.x*W.v.y*W.gm1;
  dFxdU(3,3) += W.g*W.v.x;
}

/********************************************************
 * Euler2D_pState::Fy -- Solution flux (y-direction).   *
 ********************************************************/
inline Euler2D_cState Euler2D_pState::Fy(void) {
  return (Euler2D_cState(d*v.y, d*v.x*v.y, d*sqr(v.y) + p, v.y*H()));
}

inline Euler2D_cState Euler2D_pState::Fy(void) const {
  return (Euler2D_cState(d*v.y, d*v.x*v.y, d*sqr(v.y) + p, v.y*H()));
}

inline Euler2D_cState Euler2D_pState::Fy(const Euler2D_pState &W) {
  return (Euler2D_cState(W.d*W.v.y, W.d*W.v.x*W.v.y,
                         W.d*sqr(W.v.y) + W.p, W.v.y*W.H()));
}

inline Euler2D_cState Fy(const Euler2D_pState &W) {
  return (Euler2D_cState(W.d*W.v.y, W.d*W.v.x*W.v.y,
                         W.d*sqr(W.v.y) + W.p, W.v.y*W.H()));
}

inline void Euler2D_pState::dFydU(DenseMatrix &dFydU) {
  dFydU(0,2) += ONE;
  dFydU(1,0) -= v.x*v.y;
  dFydU(1,1) += v.y;
  dFydU(1,2) += v.x;
  dFydU(2,0) += HALF*(sqr(v.y)*(g-THREE)+sqr(v.x)*gm1);
  dFydU(2,1) -= v.x*gm1;
  dFydU(2,2) -= v.y*(g-THREE);
  dFydU(2,3) += gm1;
  dFydU(3,0) += HALF*gm1i*v.y*(d*(TWO+g*g-THREE*g)*(sqr(v.x)+sqr(v.y))-TWO*p*g)/d;
  dFydU(3,1) -= v.x*v.y*gm1;
  dFydU(3,2) += HALF*gm1i*
    (sqr(v.y)*d*(FIVE*g-THREE-TWO*g*g)+sqr(v.x)*d*gm1+TWO*p*g)/d;
  dFydU(3,3) += g*v.y;
}

inline void Euler2D_pState::dFydU(DenseMatrix &dFydU) const {
  dFydU(0,2) += ONE;
  dFydU(1,0) -= v.x*v.y;
  dFydU(1,1) += v.y;
  dFydU(1,2) += v.x;
  dFydU(2,0) += HALF*(sqr(v.y)*(g-THREE)+sqr(v.x)*gm1);
  dFydU(2,1) -= v.x*gm1;
  dFydU(2,2) -= v.y*(g-THREE);
  dFydU(2,3) += gm1;
  dFydU(3,0) += HALF*gm1i*v.y*(d*(TWO+g*g-THREE*g)*(sqr(v.x)+sqr(v.y))-TWO*p*g)/d;
  dFydU(3,1) -= v.x*v.y*gm1;
  dFydU(3,2) += HALF*gm1i*
    (sqr(v.y)*d*(FIVE*g-THREE-TWO*g*g)+sqr(v.x)*d*gm1+TWO*p*g)/d;
  dFydU(3,3) += g*v.y;
}

inline void dFydU(DenseMatrix &dFydU, const Euler2D_pState &W) {
  dFydU(0,2) += ONE;
  dFydU(1,0) -= W.v.x*W.v.y;
  dFydU(1,1) += W.v.y;
  dFydU(1,2) += W.v.x;
  dFydU(2,0) += HALF*(sqr(W.v.y)*(W.g-THREE)+sqr(W.v.x)*W.gm1);
  dFydU(2,1) -= W.v.x*W.gm1;
  dFydU(2,2) -= W.v.y*(W.g-THREE);
  dFydU(2,3) += W.gm1;
  dFydU(3,0) += HALF*W.gm1i*W.v.y*
    (W.d*(TWO+W.g*W.g-THREE*W.g)*(sqr(W.v.x)+sqr(W.v.y))-TWO*W.p*W.g)/W.d;
  dFydU(3,1) -= W.v.x*W.v.y*W.gm1;
  dFydU(3,2) += HALF*W.gm1i*
    (sqr(W.v.y)*W.d*(FIVE*W.g-THREE-TWO*W.g*W.g)
     +sqr(W.v.x)*W.d*W.gm1+TWO*W.p*W.g)/W.d;
  dFydU(3,3) += W.g*W.v.y;
}

/********************************************************
 * Euler2D_pState::Fn -- Solution flux (n-direction).   *
 ********************************************************/
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

inline void Euler2D_pState::dFndU(DenseMatrix &dFndU) {
  dFndU(0,1) += ONE;
  dFndU(1,0) += HALF*(sqr(v.x)*(g-THREE)+sqr(v.y)*gm1);
  dFndU(1,1) -= v.x*(g-THREE);
  dFndU(1,2) -= v.y*gm1;
  dFndU(1,3) += gm1;
  dFndU(2,0) -= v.x*v.y;
  dFndU(2,1) += v.y;
  dFndU(2,2) += v.x;
  dFndU(3,0) += HALF*gm1i*v.x*(d*(TWO+g*g-THREE*g)*(sqr(v.x)+sqr(v.y))-TWO*p*g)/d;
  dFndU(3,1) += HALF*gm1i*
    (sqr(v.x)*d*(FIVE*g-THREE-TWO*g*g)+sqr(v.y)*d*gm1+TWO*p*g)/d;
  dFndU(3,2) -= v.x*v.y*gm1;
  dFndU(3,3) += g*v.x;
}

inline void Euler2D_pState::dFndU(DenseMatrix &dFndU) const {
  dFndU(0,1) += ONE;
  dFndU(1,0) += HALF*(sqr(v.x)*(g-THREE)+sqr(v.y)*gm1);
  dFndU(1,1) -= v.x*(g-THREE);
  dFndU(1,2) -= v.y*gm1;
  dFndU(1,3) += gm1;
  dFndU(2,0) -= v.x*v.y;
  dFndU(2,1) += v.y;
  dFndU(2,2) += v.x;
  dFndU(3,0) += HALF*gm1i*v.x*(d*(TWO+g*g-THREE*g)*(sqr(v.x)+sqr(v.y))-TWO*p*g)/d;
  dFndU(3,1) += HALF*gm1i*
    (sqr(v.x)*d*(FIVE*g-THREE-TWO*g*g)+sqr(v.y)*d*gm1+TWO*p*g)/d;
  dFndU(3,2) -= v.x*v.y*gm1;
  dFndU(3,3) += g*v.x;
}

inline void dFndU(DenseMatrix &dFndU, const Euler2D_pState &W) {
  dFndU(0,1) += ONE;
  dFndU(1,0) += HALF*(sqr(W.v.x)*(W.g-THREE)+sqr(W.v.y)*(W.g-1));
  dFndU(1,1) -= W.v.x*(W.g-THREE);
  dFndU(1,2) -= W.v.y*W.gm1;
  dFndU(1,3) += W.gm1;
  dFndU(2,0) -= W.v.x*W.v.y;
  dFndU(2,1) += W.v.y;
  dFndU(2,2) += W.v.x;
  dFndU(3,0) += HALF*W.gm1i*W.v.x*
    (W.d*(TWO+W.g*W.g-THREE*W.g)*(sqr(W.v.x)+sqr(W.v.y))-TWO*W.p*W.g)/W.d;
  dFndU(3,1) += HALF*W.gm1i*
    (sqr(W.v.x)*W.d*(FIVE*W.g-THREE-TWO*W.g*W.g)
     +sqr(W.v.y)*W.d*W.gm1+TWO*W.p*W.g)/W.d;
  dFndU(3,2) -= W.v.x*W.v.y*W.gm1;
  dFndU(3,3) += W.g*W.v.x;
}

/************************************************************
 * Euler2D_pState::lambda -- Eigenvalue(s) (x-direction).   *
 ************************************************************/
inline Euler2D_pState Euler2D_pState::lambda(void) {
  double c = a();
  return (Euler2D_pState(v.x - c, v.x, v.x, v.x + c));
}

inline Euler2D_pState Euler2D_pState::lambda(void) const {
  double c = a();
  return (Euler2D_pState(v.x - c, v.x, v.x, v.x + c));
}

inline Euler2D_pState Euler2D_pState::lambda(const Euler2D_pState &W) {
  double c = W.a();
  return (Euler2D_pState(W.v.x - c, W.v.x, W.v.x, W.v.x + c));
}

inline Euler2D_pState lambda(const Euler2D_pState &W) {
  double c = W.a();
  return (Euler2D_pState(W.v.x - c, W.v.x, W.v.x, W.v.x + c));
}

inline double Euler2D_pState::lambda(int index) {
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
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

inline double Euler2D_pState::lambda(int index) const {
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
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

inline double lambda(const Euler2D_pState &W, int index) {
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
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

inline Euler2D_pState Euler2D_pState::lambda_x(const Vector2D &V) const {
  double c = a();
  return Euler2D_pState(v.x-V.x-c,v.x-V.x,v.x-V.x,v.x-V.x+c);
}

/************************************************************
 * Euler2D_pState::lambda_x -- Eigenvalue(s) (x-direction). *
 ************************************************************/
inline Euler2D_pState Euler2D_pState::lambda_x(void) {
  double c = a();
  return (Euler2D_pState(v.x - c, v.x, v.x, v.x + c));
}

inline Euler2D_pState Euler2D_pState::lambda_x(void) const {
  double c = a();
  return (Euler2D_pState(v.x - c, v.x, v.x, v.x + c));
}

inline Euler2D_pState Euler2D_pState::lambda_x(const Euler2D_pState &W) {
  double c = W.a();
  return (Euler2D_pState(W.v.x - c, W.v.x, W.v.x, W.v.x + c));
}

inline Euler2D_pState lambda_x(const Euler2D_pState &W) {
  double c = W.a();
  return (Euler2D_pState(W.v.x - c, W.v.x, W.v.x, W.v.x + c));
}

inline double Euler2D_pState::lambda_x(int index) {
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
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

inline double Euler2D_pState::lambda_x(int index) const {
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
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

inline double lambda_x(const Euler2D_pState &W, int index) {
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
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

/************************************************************
 * Euler2D_pState::lambda_precon_WS -- Eigenvalue(s) of     *
 * Weiss-Smith preconditioned system (x-direction).         *
 ************************************************************/
inline Euler2D_pState Euler2D_pState::lambda_precon_WS(void) {
  double c = a(), mr = Mr(), cc;
  cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
  return (Euler2D_pState(v.x*HALF*(ONE+mr*mr) - cc, v.x, v.x, v.x*HALF*(ONE+mr*mr) + cc));
}

inline Euler2D_pState Euler2D_pState::lambda_precon_WS(void) const {
  double c = a(), mr = Mr(), cc;
  cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
  return (Euler2D_pState(v.x*HALF*(ONE+mr*mr) - cc, v.x, v.x, v.x*HALF*(ONE+mr*mr) + cc));
}

inline Euler2D_pState Euler2D_pState::lambda_precon_WS(const Euler2D_pState &W) {
  double c = W.a(), mr = W.Mr(), cc;
  cc = sqrt(QUARTER*sqr(W.v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
  return (Euler2D_pState(W.v.x*HALF*(ONE+mr*mr) - cc, W.v.x, W.v.x, W.v.x*HALF*(ONE+mr*mr) + cc));
}

inline Euler2D_pState lambda_precon_WS(const Euler2D_pState &W) {
  double c = W.a(), mr = W.Mr(), cc;
  cc = sqrt(QUARTER*sqr(W.v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
  return (Euler2D_pState(W.v.x*HALF*(ONE+mr*mr) - cc, W.v.x, W.v.x, W.v.x*HALF*(ONE+mr*mr) + cc));
}

inline double Euler2D_pState::lambda_precon_WS(int index) {
  double c, mr, cc;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (v.x*HALF*(ONE+mr*mr)-cc);
    case 2 :
      return (v.x);
    case 3 :
      return (v.x);
    case 4 :
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (v.x*HALF*(ONE+mr*mr)+cc);
    default:
      return (v.x);
  };
}

inline double Euler2D_pState::lambda_precon_WS(int index) const {
  double c, mr, cc;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (v.x*HALF*(ONE+mr*mr)-cc);
    case 2 :
      return (v.x);
    case 3 :
      return (v.x);
    case 4 :
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (v.x*HALF*(ONE+mr*mr)+cc);
    default:
      return (v.x);
  };
}

inline double lambda_precon_WS(const Euler2D_pState &W, int index) {
  double c, mr, cc;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = W.a(); mr = W.Mr(); cc = sqrt(QUARTER*sqr(W.v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (W.v.x*HALF*(ONE+mr*mr)-cc);
    case 2 :
      return (W.v.x);
    case 3 :
      return (W.v.x);
    case 4 :
      c = W.a(); mr = W.Mr(); cc = sqrt(QUARTER*sqr(W.v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (W.v.x*HALF*(ONE+mr*mr)+cc);
    default:
      return (W.v.x);
  };
}

/************************************************************
 * Euler2D_pState::lambda_y -- Eigenvalue(s) (y-direction). *
 ************************************************************/
inline Euler2D_pState Euler2D_pState::lambda_y(void) {
  double c = a();
  return (Euler2D_pState(v.y - c, v.y, v.y, v.y + c));
}

inline Euler2D_pState Euler2D_pState::lambda_y(void) const {
  double c = a();
  return (Euler2D_pState(v.y - c, v.y, v.y, v.y + c));
}

inline Euler2D_pState Euler2D_pState::lambda_y(const Euler2D_pState &W) {
  double c = W.a();
  return (Euler2D_pState(W.v.y - c, W.v.y, W.v.y, W.v.y + c));
}

inline Euler2D_pState lambda_y(const Euler2D_pState &W) {
  double c = W.a();
  return (Euler2D_pState(W.v.y - c, W.v.y, W.v.y, W.v.y + c));
}

inline double Euler2D_pState::lambda_y(int index) {
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
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

inline double Euler2D_pState::lambda_y(int index) const {
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
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

inline double lambda_y(const Euler2D_pState &W, int index) {
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
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

/********************************************************
 * Euler2D_pState::rp -- Primitive right eigenvector    *
 *                       (x-direction).                 *
 ********************************************************/
inline Euler2D_pState Euler2D_pState::rp(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler2D_pState(ONE, -c/d, ZERO, sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Euler2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = a();
      return (Euler2D_pState(ONE, c/d, ZERO, sqr(c)));
    default:
      c = a();
      return (Euler2D_pState(ONE, -c/d, ZERO, sqr(c)));
  };
}

inline Euler2D_pState Euler2D_pState::rp(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler2D_pState(ONE, -c/d, ZERO, sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Euler2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = a();
      return (Euler2D_pState(ONE, c/d, ZERO, sqr(c)));
    default:
      c = a();
      return (Euler2D_pState(ONE, -c/d, ZERO, sqr(c)));
  };
}

inline Euler2D_pState rp(const Euler2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Euler2D_pState(ONE, -c/W.d, ZERO, sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Euler2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = W.a();
      return (Euler2D_pState(ONE, c/W.d, ZERO, sqr(c)));
    default:
      c = W.a();
      return (Euler2D_pState(ONE, -c/W.d, ZERO, sqr(c)));
  };
}

/********************************************************
 * Euler2D_pState::rp_x -- Primitive right eigenvector  *
 *                         (x-direction).               *
 ********************************************************/
inline Euler2D_pState Euler2D_pState::rp_x(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler2D_pState(ONE, -c/d, ZERO, sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Euler2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = a();
      return (Euler2D_pState(ONE, c/d, ZERO, sqr(c)));
    default:
      c = a();
      return (Euler2D_pState(ONE, -c/d, ZERO, sqr(c)));
  };
}

inline Euler2D_pState Euler2D_pState::rp_x(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler2D_pState(ONE, -c/d, ZERO, sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Euler2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = a();
      return (Euler2D_pState(ONE, c/d, ZERO, sqr(c)));
    default:
      c = a();
      return (Euler2D_pState(ONE, -c/d, ZERO, sqr(c)));
  };
}

inline Euler2D_pState rp_x(const Euler2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Euler2D_pState(ONE, -c/W.d, ZERO, sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Euler2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = W.a();
      return (Euler2D_pState(ONE, c/W.d, ZERO, sqr(c)));
    default:
      c = W.a();
      return (Euler2D_pState(ONE, -c/W.d, ZERO, sqr(c)));
  };
}

/********************************************************
 * Euler2D_pState::rp_precon_WS -- Primitive right      *
 * eigenvector of Weiss-Smith preconditioned system     *
 * (x-direction).                                       *
 ********************************************************/
inline Euler2D_pState Euler2D_pState::rp_precon_WS(int index) {
  double c, mr, cc;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_pState(ONE, -(cc+v.x*HALF*(mr*mr-ONE))/(d*mr*mr), ZERO, sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Euler2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_pState(ONE, (cc-v.x*HALF*(mr*mr-ONE))/(d*mr*mr), ZERO, sqr(c)));
    default:
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_pState(ONE, -(cc+v.x*HALF*(mr*mr-ONE))/(d*mr*mr), ZERO, sqr(c)));
  };
}

inline Euler2D_pState Euler2D_pState::rp_precon_WS(int index) const {
  double c, mr, cc;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_pState(ONE, -(cc+v.x*HALF*(mr*mr-ONE))/(d*mr*mr), ZERO, sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Euler2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_pState(ONE, (cc-v.x*HALF*(mr*mr-ONE))/(d*mr*mr), ZERO, sqr(c)));
    default:
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_pState(ONE, -(cc+v.x*HALF*(mr*mr-ONE))/(d*mr*mr), ZERO, sqr(c)));
  };
}

inline Euler2D_pState rp_precon_WS(const Euler2D_pState &W, int index) {
  double c, mr, cc;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = W.a(); mr = W.Mr(); cc = sqrt(QUARTER*sqr(W.v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_pState(ONE, -(cc+W.v.x*HALF*(mr*mr-ONE))/(W.d*mr*mr), ZERO, sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Euler2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = W.a(); mr = W.Mr(); cc = sqrt(QUARTER*sqr(W.v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_pState(ONE, (cc-W.v.x*HALF*(mr*mr-ONE))/(W.d*mr*mr), ZERO, sqr(c)));
    default:
      c = W.a(); mr = W.Mr(); cc = sqrt(QUARTER*sqr(W.v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_pState(ONE, -(cc+W.v.x*HALF*(mr*mr-ONE))/(W.d*mr*mr), ZERO, sqr(c)));
  };
}

/********************************************************
 * Euler2D_pState::rp_y -- Primitive right eigenvector  *
 *                         (y-direction).               *
 ********************************************************/
inline Euler2D_pState Euler2D_pState::rp_y(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler2D_pState(ONE, ZERO, -c/d, sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Euler2D_pState(ZERO, ONE, ZERO, ZERO));
    case 4 :
      c = a();
      return (Euler2D_pState(ONE, ZERO, c/d, sqr(c)));
    default:
      c = a();
      return (Euler2D_pState(ONE, ZERO, -c/d, sqr(c)));
  };
}

inline Euler2D_pState Euler2D_pState::rp_y(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler2D_pState(ONE, ZERO, -c/d, sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Euler2D_pState(ZERO, ONE, ZERO, ZERO));
    case 4 :
      c = a();
      return (Euler2D_pState(ONE, ZERO, c/d, sqr(c)));
    default:
      c = a();
      return (Euler2D_pState(ONE, ZERO, -c/d, sqr(c)));
  };
}

inline Euler2D_pState rp_y(const Euler2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Euler2D_pState(ONE, ZERO, -c/W.d, sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, ZERO));
    case 3 :
      return (Euler2D_pState(ZERO, ONE, ZERO, ZERO));
    case 4 :
      c = W.a();
      return (Euler2D_pState(ONE, ZERO, c/W.d, sqr(c)));
    default:
      c = W.a();
      return (Euler2D_pState(ONE, ZERO, -c/W.d, sqr(c)));
  };
}

/********************************************************
 * Euler2D_pState::rc -- Conserved right eigenvector    *
 *                       (x-direction).                 *
 ********************************************************/
inline Euler2D_cState Euler2D_pState::rc(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler2D_cState(ONE, v.x-c, v.y, h()-v.x*c));
    case 2 :
      return (Euler2D_cState(ONE, v.x, v.y, HALF*v.sqr()));
    case 3 :
      return (Euler2D_cState(ZERO, ZERO, d, d*v.y));
    case 4 :
      c = a();
      return (Euler2D_cState(ONE, v.x+c, v.y, h()+v.x*c));
    default:
      c = a();
      return (Euler2D_cState(ONE, v.x-c, v.y, h()-v.x*c));
  };
}

inline Euler2D_cState Euler2D_pState::rc(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler2D_cState(ONE, v.x-c, v.y, h()-v.x*c));
    case 2 :
      return (Euler2D_cState(ONE, v.x, v.y, HALF*v.sqr()));
    case 3 :
      return (Euler2D_cState(ZERO, ZERO, d, d*v.y));
    case 4 :
      c = a();
      return (Euler2D_cState(ONE, v.x+c, v.y, h()+v.x*c));
    default:
      c = a();
      return (Euler2D_cState(ONE, v.x-c, v.y, h()-v.x*c));
  };
}

inline Euler2D_cState rc(const Euler2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Euler2D_cState(ONE, W.v.x-c, W.v.y, W.h()-W.v.x*c));
    case 2 :
      return (Euler2D_cState(ONE, W.v.x, W.v.y, HALF*W.v.sqr()));
    case 3 :
      return (Euler2D_cState(ZERO, ZERO, W.d, W.d*W.v.y));
    case 4 :
      c = W.a();
      return (Euler2D_cState(ONE, W.v.x+c, W.v.y, W.h()+W.v.x*c));
    default:
      c = W.a();
      return (Euler2D_cState(ONE, W.v.x-c, W.v.y, W.h()-W.v.x*c));
    };
}

/********************************************************
 * Euler2D_pState::rc_x -- Conserved right eigenvector  *
 *                         (x-direction).               *
 ********************************************************/
inline Euler2D_cState Euler2D_pState::rc_x(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler2D_cState(ONE, v.x-c, v.y, h()-v.x*c));
    case 2 :
      return (Euler2D_cState(ONE, v.x, v.y, HALF*v.sqr()));
    case 3 :
      return (Euler2D_cState(ZERO, ZERO, d, d*v.y));
    case 4 :
      c = a();
      return (Euler2D_cState(ONE, v.x+c, v.y, h()+v.x*c));
    default:
      c = a();
      return (Euler2D_cState(ONE, v.x-c, v.y, h()-v.x*c));
  };
}

inline Euler2D_cState Euler2D_pState::rc_x(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler2D_cState(ONE, v.x-c, v.y, h()-v.x*c));
    case 2 :
      return (Euler2D_cState(ONE, v.x, v.y, HALF*v.sqr()));
    case 3 :
      return (Euler2D_cState(ZERO, ZERO, d, d*v.y));
    case 4 :
      c = a();
      return (Euler2D_cState(ONE, v.x+c, v.y, h()+v.x*c));
    default:
      c = a();
      return (Euler2D_cState(ONE, v.x-c, v.y, h()-v.x*c));
  };
}

inline Euler2D_cState rc_x(const Euler2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Euler2D_cState(ONE, W.v.x-c, W.v.y, W.h()-W.v.x*c));
    case 2 :
      return (Euler2D_cState(ONE, W.v.x, W.v.y, HALF*W.v.sqr()));
    case 3 :
      return (Euler2D_cState(ZERO, ZERO, W.d, W.d*W.v.y));
    case 4 :
      c = W.a();
      return (Euler2D_cState(ONE, W.v.x+c, W.v.y, W.h()+W.v.x*c));
    default:
      c = W.a();
      return (Euler2D_cState(ONE, W.v.x-c, W.v.y, W.h()-W.v.x*c));
    };
}

/********************************************************
 * Euler2D_pState::rc_precon_WS -- Conservative right   *
 * eigenvector of Weiss-Smith preconditioned system     *
 * (x-direction).                                       *
 ********************************************************/
inline Euler2D_cState Euler2D_pState::rc_precon_WS(int index) {
  double c, mr, cc;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_cState(ONE, (v.x*HALF*(ONE+mr*mr) - cc)/(mr*mr), 
                             v.y, h()-v.x*(cc+v.x*HALF*(mr*mr-ONE))/(mr*mr)));
    case 2 :
      return (Euler2D_cState(ONE, v.x, v.y, HALF*v.sqr()));
    case 3 :
      return (Euler2D_cState(ZERO, ZERO, d, d*v.y));
    case 4 :
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_cState(ONE, (v.x*HALF*(ONE+mr*mr) + cc)/(mr*mr), 
                             v.y, h()+v.x*(cc-v.x*HALF*(mr*mr-ONE))/(mr*mr)));
    default:
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_cState(ONE, (v.x*HALF*(ONE+mr*mr) - cc)/(mr*mr), 
                             v.y, h()-v.x*(cc+v.x*HALF*(mr*mr-ONE))/(mr*mr)));
  };
}

inline Euler2D_cState Euler2D_pState::rc_precon_WS(int index) const {
  double c, mr, cc;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_cState(ONE, (v.x*HALF*(ONE+mr*mr) - cc)/(mr*mr), 
                             v.y, h()-v.x*(cc+v.x*HALF*(mr*mr-ONE))/(mr*mr)));
    case 2 :
      return (Euler2D_cState(ONE, v.x, v.y, HALF*v.sqr()));
    case 3 :
      return (Euler2D_cState(ZERO, ZERO, d, d*v.y));
    case 4 :
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_cState(ONE, (v.x*HALF*(ONE+mr*mr) + cc)/(mr*mr), 
                             v.y, h()+v.x*(cc-v.x*HALF*(mr*mr-ONE))/(mr*mr)));
    default:
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_cState(ONE, (v.x*HALF*(ONE+mr*mr) - cc)/(mr*mr), 
                             v.y, h()-v.x*(cc+v.x*HALF*(mr*mr-ONE))/(mr*mr)));
  };
}

inline Euler2D_cState rc_precon_WS(const Euler2D_pState &W, int index) {
  double c, mr, cc;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = W.a(); mr = W.Mr(); cc = sqrt(QUARTER*sqr(W.v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_cState(ONE, (W.v.x*HALF*(ONE+mr*mr) - cc)/(mr*mr), 
                             W.v.y, W.h()-W.v.x*(cc+W.v.x*HALF*(mr*mr-ONE))/(mr*mr)));
    case 2 :
      return (Euler2D_cState(ONE, W.v.x, W.v.y, HALF*W.v.sqr()));
    case 3 :
      return (Euler2D_cState(ZERO, ZERO, W.d, W.d*W.v.y));
    case 4 :
      c = W.a(); mr = W.Mr(); cc = sqrt(QUARTER*sqr(W.v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_cState(ONE, (W.v.x*HALF*(ONE+mr*mr) + cc)/(mr*mr), 
                             W.v.y, W.h()+W.v.x*(cc-W.v.x*HALF*(mr*mr-ONE))/(mr*mr)));
    default:
      c = W.a(); mr = W.Mr(); cc = sqrt(QUARTER*sqr(W.v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_cState(ONE, (W.v.x*HALF*(ONE+mr*mr) - cc)/(mr*mr), 
                             W.v.y, W.h()-W.v.x*(cc+W.v.x*HALF*(mr*mr-ONE))/(mr*mr)));
    };
}

/********************************************************
 * Euler2D_pState::rc_y -- Conserved right eigenvector  *
 *                         (y-direction).               *
 ********************************************************/
inline Euler2D_cState Euler2D_pState::rc_y(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler2D_cState(ONE, v.x, v.y-c, h()-v.y*c));
    case 2 :
      return (Euler2D_cState(ONE, v.x, v.y, HALF*v.sqr()));
    case 3 :
      return (Euler2D_cState(ZERO, d, ZERO, d*v.x));
    case 4 :
      c = a();
      return (Euler2D_cState(ONE, v.x, v.y+c, h()+v.y*c));
    default:
      c = a();
      return (Euler2D_cState(ONE, v.x, v.y-c, h()-v.y*c));
  };
}

inline Euler2D_cState Euler2D_pState::rc_y(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler2D_cState(ONE, v.x, v.y-c, h()-v.y*c));
    case 2 :
      return (Euler2D_cState(ONE, v.x, v.y, HALF*v.sqr()));
    case 3 :
      return (Euler2D_cState(ZERO, d, ZERO, d*v.x));
    case 4 :
      c = a();
      return (Euler2D_cState(ONE, v.x, v.y+c, h()+v.y*c));
    default:
      c = a();
      return (Euler2D_cState(ONE, v.x, v.y-c, h()-v.y*c));
  };
}

inline Euler2D_cState rc_y(const Euler2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Euler2D_cState(ONE, W.v.x, W.v.y-c, W.h()-W.v.y*c));
    case 2 :
      return (Euler2D_cState(ONE, W.v.x, W.v.y, HALF*W.v.sqr()));
    case 3 :
      return (Euler2D_cState(ZERO, W.d, ZERO, W.d*W.v.x));
    case 4 :
      c = W.a();
      return (Euler2D_cState(ONE, W.v.x, W.v.y+c, W.h()+W.v.y*c));
    default:
      c = W.a();
      return (Euler2D_cState(ONE, W.v.x, W.v.y-c, W.h()-W.v.y*c));
    };
}

/********************************************************
 * Euler2D_pState::lp -- Primitive left eigenvector     *
 *                       (x-direction).                 *
 ********************************************************/
inline Euler2D_pState Euler2D_pState::lp(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler2D_pState(ZERO, -HALF*d/c, ZERO, HALF/sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, -ONE/a2()));
    case 3 :
      return (Euler2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = a();
      return (Euler2D_pState(ZERO, HALF*d/c, ZERO, HALF/sqr(c)));
    default:
      c = a();
      return (Euler2D_pState(ZERO, -HALF*d/c, ZERO, HALF/sqr(c)));
  };
}

inline Euler2D_pState Euler2D_pState::lp(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler2D_pState(ZERO, -HALF*d/c, ZERO, HALF/sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, -ONE/a2()));
    case 3 :
      return (Euler2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = a();
      return (Euler2D_pState(ZERO, HALF*d/c, ZERO, HALF/sqr(c)));
    default:
      c = a();
      return (Euler2D_pState(ZERO, -HALF*d/c, ZERO, HALF/sqr(c)));
  };
}

inline Euler2D_pState lp(const Euler2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Euler2D_pState(ZERO, -HALF*W.d/c, ZERO, HALF/sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, -ONE/W.a2()));
    case 3 :
      return (Euler2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = W.a();
      return (Euler2D_pState(ZERO, HALF*W.d/c, ZERO, HALF/sqr(c)));
    default:
      c = W.a();
      return (Euler2D_pState(ZERO, -HALF*W.d/c, ZERO, HALF/sqr(c)));
  };
}

/********************************************************
 * Euler2D_pState::lp_x -- Primitive left eigenvector   *
 *                         (x-direction).               *
 ********************************************************/
inline Euler2D_pState Euler2D_pState::lp_x(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler2D_pState(ZERO, -HALF*d/c, ZERO, HALF/sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, -ONE/a2()));
    case 3 :
      return (Euler2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = a();
      return (Euler2D_pState(ZERO, HALF*d/c, ZERO, HALF/sqr(c)));
    default:
      c = a();
      return (Euler2D_pState(ZERO, -HALF*d/c, ZERO, HALF/sqr(c)));
  };
}

inline Euler2D_pState Euler2D_pState::lp_x(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler2D_pState(ZERO, -HALF*d/c, ZERO, HALF/sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, -ONE/a2()));
    case 3 :
      return (Euler2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = a();
      return (Euler2D_pState(ZERO, HALF*d/c, ZERO, HALF/sqr(c)));
    default:
      c = a();
      return (Euler2D_pState(ZERO, -HALF*d/c, ZERO, HALF/sqr(c)));
  };
}

inline Euler2D_pState lp_x(const Euler2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Euler2D_pState(ZERO, -HALF*W.d/c, ZERO, HALF/sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, -ONE/W.a2()));
    case 3 :
      return (Euler2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = W.a();
      return (Euler2D_pState(ZERO, HALF*W.d/c, ZERO, HALF/sqr(c)));
    default:
      c = W.a();
      return (Euler2D_pState(ZERO, -HALF*W.d/c, ZERO, HALF/sqr(c)));
  };
}

/********************************************************
 * Euler2D_pState::lp_precon_WS -- Primitive left       *
 * eigenvector of Weiss-Smith preconditioned system     * 
 * (x-direction).                                       *
 ********************************************************/
inline Euler2D_pState Euler2D_pState::lp_precon_WS(int index) {
  double c, mr, cc;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_pState(ZERO, -HALF*d*mr*mr/cc, 
                             ZERO, HALF*(cc-v.x*HALF*(mr*mr-ONE))/(sqr(c)*cc)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, -ONE/a2()));
    case 3 :
      return (Euler2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_pState(ZERO, HALF*d*mr*mr/cc, 
                             ZERO, HALF*(cc+v.x*HALF*(mr*mr-ONE))/(sqr(c)*cc)));
    default:
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_pState(ZERO, -HALF*d*mr*mr/cc, 
                             ZERO, HALF*(cc-v.x*HALF*(mr*mr-ONE))/(sqr(c)*cc)));
  };
}

inline Euler2D_pState Euler2D_pState::lp_precon_WS(int index) const {
  double c, mr, cc;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_pState(ZERO, -HALF*d*mr*mr/cc, 
                             ZERO, HALF*(cc-v.x*HALF*(mr*mr-ONE))/(sqr(c)*cc)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, -ONE/a2()));
    case 3 :
      return (Euler2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_pState(ZERO, HALF*d*mr*mr/cc, 
                             ZERO, HALF*(cc+v.x*HALF*(mr*mr-ONE))/(sqr(c)*cc)));
    default:
      c = a(); mr = Mr(); cc = sqrt(QUARTER*sqr(v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_pState(ZERO, -HALF*d*mr*mr/cc, 
                             ZERO, HALF*(cc-v.x*HALF*(mr*mr-ONE))/(sqr(c)*cc)));
  };
}

inline Euler2D_pState lp_precon_WS(const Euler2D_pState &W, int index) {
  double c, mr, cc;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = W.a(); mr = W.Mr(); cc = sqrt(QUARTER*sqr(W.v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_pState(ZERO, -HALF*W.d*mr*mr/cc, 
                             ZERO, HALF*(cc-W.v.x*HALF*(mr*mr-ONE))/(sqr(c)*cc)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, -ONE/W.a2()));
    case 3 :
      return (Euler2D_pState(ZERO, ZERO, ONE, ZERO));
    case 4 :
      c = W.a(); mr = W.Mr(); cc = sqrt(QUARTER*sqr(W.v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_pState(ZERO, HALF*W.d*mr*mr/cc, 
                             ZERO, HALF*(cc+W.v.x*HALF*(mr*mr-ONE))/(sqr(c)*cc)));
    default:
      c = W.a(); mr = W.Mr(); cc = sqrt(QUARTER*sqr(W.v.x)*(mr*mr*mr*mr-TWO*mr*mr+ONE)+c*c*mr*mr);
      return (Euler2D_pState(ZERO, -HALF*W.d*mr*mr/cc, 
                             ZERO, HALF*(cc-W.v.x*HALF*(mr*mr-ONE))/(sqr(c)*cc)));
  };
}

/********************************************************
 * Euler2D_pState::lp_y -- Primitive left eigenvector   *
 *                         (y-direction).               *
 ********************************************************/
inline Euler2D_pState Euler2D_pState::lp_y(int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler2D_pState(ZERO, ZERO, -HALF*d/c, HALF/sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, -ONE/a2()));
    case 3 :
      return (Euler2D_pState(ZERO, ONE, ZERO, ZERO));
    case 4 :
      c = a();
      return (Euler2D_pState(ZERO, ZERO, HALF*d/c, HALF/sqr(c)));
    default:
      c = a();
      return (Euler2D_pState(ZERO, ZERO, -HALF*d/c, HALF/sqr(c)));
  };
}

inline Euler2D_pState Euler2D_pState::lp_y(int index) const {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = a();
      return (Euler2D_pState(ZERO, ZERO, -HALF*d/c, HALF/sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, -ONE/a2()));
    case 3 :
      return (Euler2D_pState(ZERO, ONE, ZERO, ZERO));
    case 4 :
      c = a();
      return (Euler2D_pState(ZERO, ZERO, HALF*d/c, HALF/sqr(c)));
    default:
      c = a();
      return (Euler2D_pState(ZERO, ZERO, -HALF*d/c, HALF/sqr(c)));
  };
}

inline Euler2D_pState lp_y(const Euler2D_pState &W, int index) {
  double c;
  assert( index >= 1 && index <= NUM_VAR_EULER2D );
  switch(index) {
    case 1 :
      c = W.a();
      return (Euler2D_pState(ZERO, ZERO, -HALF*W.d/c, HALF/sqr(c)));
    case 2 :
      return (Euler2D_pState(ONE, ZERO, ZERO, -ONE/W.a2()));
    case 3 :
      return (Euler2D_pState(ZERO, ONE, ZERO, ZERO));
    case 4 :
      c = W.a();
      return (Euler2D_pState(ZERO, ZERO, HALF*W.d/c, HALF/sqr(c)));
    default:
      c = W.a();
      return (Euler2D_pState(ZERO, ZERO, -HALF*W.d/c, HALF/sqr(c)));
  };
}

/**********************************************************
 * Euler2D_pState::S -- Source terms (axisymmetric flow). *
 **********************************************************/
inline Euler2D_cState Euler2D_pState::S(const Vector2D &X) {
  return (Euler2D_cState(-d*v.y/X.y, -d*v.x*v.y/X.y, 
                         -d*sqr(v.y)/X.y, -v.y*H()/X.y));
}

inline Euler2D_cState Euler2D_pState::S(const Vector2D &X) const {
  return (Euler2D_cState(-d*v.y/X.y, -d*v.x*v.y/X.y, 
                         -d*sqr(v.y)/X.y, -v.y*H()/X.y));
}

inline Euler2D_cState S(const Euler2D_pState &W, const Vector2D &X) {
  return (Euler2D_cState(-W.d*W.v.y/X.y, -W.d*W.v.x*W.v.y/X.y, 
                         -W.d*sqr(W.v.y)/X.y, -W.v.y*W.H()/X.y));
}

/**********************************************************
 * Euler2D_pState::dUdW -- Solution Jacobian.             *
 **********************************************************/
inline void Euler2D_pState::dUdW(DenseMatrix &dUdW) {
  dUdW(0,0) += ONE;
  dUdW(1,0) += v.x;
  dUdW(1,1) += d;
  dUdW(2,0) += v.y;
  dUdW(2,2) += d;
  dUdW(3,0) += HALF*(sqr(v.x)+sqr(v.y));
  dUdW(3,1) += d*v.x;
  dUdW(3,2) += d*v.y;
  dUdW(3,3) += gm1i;
}

inline void Euler2D_pState::dUdW(DenseMatrix &dUdW) const {
  dUdW(0,0) += ONE;
  dUdW(1,0) += v.x;
  dUdW(1,1) += d;
  dUdW(2,0) += v.y;
  dUdW(2,2) += d;
  dUdW(3,0) += HALF*(sqr(v.x)+sqr(v.y));
  dUdW(3,1) += d*v.x;
  dUdW(3,2) += d*v.y;
  dUdW(3,3) += gm1i;
}

inline void dUdW(DenseMatrix &dUdW, const Euler2D_pState &W) {
  dUdW(0,0) += ONE;
  dUdW(1,0) += W.v.x;
  dUdW(1,1) += W.d;
  dUdW(2,0) += W.v.y;
  dUdW(2,2) += W.d;
  dUdW(3,0) += HALF*(sqr(W.v.x)+sqr(W.v.y));
  dUdW(3,1) += W.d*W.v.x;
  dUdW(3,2) += W.d*W.v.y;
  dUdW(3,3) += W.gm1i;
}

/**********************************************************
 * Euler2D_pState::dWdU -- Solution Jacobian.             *
 **********************************************************/
inline void Euler2D_pState::dWdU(DenseMatrix &dWdU) {
  dWdU(0,0) += ONE;
  dWdU(1,0) -= v.x/d;
  dWdU(1,1) += ONE/d;
  dWdU(2,0) -= v.y/d;
  dWdU(2,2) += ONE/d;
  dWdU(3,0) += HALF*gm1*(sqr(v.x)+sqr(v.y));
  dWdU(3,1) -= v.x*gm1;
  dWdU(3,2) -= v.y*gm1;
  dWdU(3,3) += gm1;
}

inline void Euler2D_pState::dWdU(DenseMatrix &dWdU) const {
  dWdU(0,0) += ONE;
  dWdU(1,0) -= v.x/d;
  dWdU(1,1) += ONE/d;
  dWdU(2,0) -= v.y/d;
  dWdU(2,2) += ONE/d;
  dWdU(3,0) += HALF*gm1*(sqr(v.x)+sqr(v.y));
  dWdU(3,1) -= v.x*gm1;
  dWdU(3,2) -= v.y*gm1;
  dWdU(3,3) += gm1;
}

inline void dWdU(DenseMatrix &dWdU, const Euler2D_pState &W) {
  dWdU(0,0) += ONE;
  dWdU(1,0) -= W.v.x/W.d;
  dWdU(1,1) += ONE/W.d;
  dWdU(2,0) -= W.v.y/W.d;
  dWdU(2,2) += ONE/W.d;
  dWdU(3,0) += HALF*W.gm1*(sqr(W.v.x)+sqr(W.v.y));
  dWdU(3,1) -= W.v.x*W.gm1;
  dWdU(3,2) -= W.v.y*W.gm1;
  dWdU(3,3) += W.gm1;
}

/**********************************************************
 * Euler2D_pState::P_U_WS -- Weiss-Smith low-Mach-number  *
 * local preconditioner.                                  *
 **********************************************************/
inline void Euler2D_pState::P_U_WS(DenseMatrix &P_U_WS) {
  double uu = v.x, vv = v.y, c = a(), mr = Mr();
  P_U_WS(0,0) = -(uu*uu*g*mr*mr-uu*uu*g-uu*uu*mr*mr+uu*uu+vv*vv*g*mr*mr-
                  vv*vv*g-vv*vv*mr*mr+vv*vv-TWO*c*c*mr*mr)/(c*c)/(mr*mr)/TWO;
  P_U_WS(0,1) = uu*(mr-ONE)*(mr+ONE)*gm1/(c*c)/(mr*mr);
  P_U_WS(0,2) = vv*(mr-ONE)*(mr+ONE)*gm1/(c*c)/(mr*mr);
  P_U_WS(0,3) = -(mr-ONE)*(mr+ONE)*gm1/(c*c)/(mr*mr);
  P_U_WS(1,0) = -uu*(uu*uu+vv*vv)*(mr-ONE)*(mr+ONE)*gm1/(c*c)/(mr*mr)/TWO;
  P_U_WS(1,1) = (uu*uu*g*mr*mr-uu*uu*g-uu*uu*mr*mr+uu*uu+c*c*mr*mr)/(c*c)/(mr*mr);
  P_U_WS(1,2) = uu*vv*(mr-ONE)*(mr+ONE)*gm1/(c*c)/(mr*mr);
  P_U_WS(1,3) = -uu*(mr-ONE)*(mr+ONE)*gm1/(c*c)/(mr*mr);
  P_U_WS(2,0) = -vv*(uu*uu+vv*vv)*(mr-ONE)*(mr+ONE)*gm1/(c*c)/(mr*mr)/TWO;
  P_U_WS(2,1) = uu*vv*(mr-ONE)*(mr+ONE)*gm1/(c*c)/(mr*mr);
  P_U_WS(2,2) = (vv*vv*g*mr*mr-vv*vv*g-vv*vv*mr*mr+vv*vv+c*c*mr*mr)/(c*c)/(mr*mr);
  P_U_WS(2,3) = -vv*(mr-ONE)*(mr+ONE)*gm1/(c*c)/(mr*mr);
  P_U_WS(3,0) = -(uu*uu+vv*vv)*(mr-ONE)*(mr+ONE)*(uu*uu*g-uu*uu+vv*vv*g-vv*vv+TWO*c*c)/
                 (c*c)/(mr*mr)/FOUR;
  P_U_WS(3,1) = uu*(mr-ONE)*(mr+ONE)*(uu*uu*g-uu*uu+vv*vv*g-vv*vv+TWO*c*c)/(c*c)/
                (mr*mr)/TWO;
  P_U_WS(3,2) = vv*(mr-ONE)*(mr+ONE)*(uu*uu*g-uu*uu+vv*vv*g-vv*vv+TWO*c*c)/(c*c)/(mr*mr)/TWO;
  P_U_WS(3,3) = -(uu*uu*g*mr*mr-uu*uu*g-uu*uu*mr*mr+uu*uu+vv*vv*g*mr*mr-
                  vv*vv*g-vv*vv*mr*mr+vv*vv-TWO*c*c)/(c*c)/(mr*mr)/TWO;
}

inline void Euler2D_pState::P_U_WS(DenseMatrix &P_U_WS) const {
  double uu = v.x, vv = v.y, c = a(), mr = Mr();
  P_U_WS(0,0) = -(uu*uu*g*mr*mr-uu*uu*g-uu*uu*mr*mr+uu*uu+vv*vv*g*mr*mr-
                  vv*vv*g-vv*vv*mr*mr+vv*vv-TWO*c*c*mr*mr)/(c*c)/(mr*mr)/TWO;
  P_U_WS(0,1) = uu*(mr-ONE)*(mr+ONE)*gm1/(c*c)/(mr*mr);
  P_U_WS(0,2) = vv*(mr-ONE)*(mr+ONE)*gm1/(c*c)/(mr*mr);
  P_U_WS(0,3) = -(mr-ONE)*(mr+ONE)*gm1/(c*c)/(mr*mr);
  P_U_WS(1,0) = -uu*(uu*uu+vv*vv)*(mr-ONE)*(mr+ONE)*gm1/(c*c)/(mr*mr)/TWO;
  P_U_WS(1,1) = (uu*uu*g*mr*mr-uu*uu*g-uu*uu*mr*mr+uu*uu+c*c*mr*mr)/(c*c)/(mr*mr);
  P_U_WS(1,2) = uu*vv*(mr-ONE)*(mr+ONE)*gm1/(c*c)/(mr*mr);
  P_U_WS(1,3) = -uu*(mr-ONE)*(mr+ONE)*gm1/(c*c)/(mr*mr);
  P_U_WS(2,0) = -vv*(uu*uu+vv*vv)*(mr-ONE)*(mr+ONE)*gm1/(c*c)/(mr*mr)/TWO;
  P_U_WS(2,1) = uu*vv*(mr-ONE)*(mr+ONE)*gm1/(c*c)/(mr*mr);
  P_U_WS(2,2) = (vv*vv*g*mr*mr-vv*vv*g-vv*vv*mr*mr+vv*vv+c*c*mr*mr)/(c*c)/(mr*mr);
  P_U_WS(2,3) = -vv*(mr-ONE)*(mr+ONE)*gm1/(c*c)/(mr*mr);
  P_U_WS(3,0) = -(uu*uu+vv*vv)*(mr-ONE)*(mr+ONE)*(uu*uu*g-uu*uu+vv*vv*g-vv*vv+TWO*c*c)/
                 (c*c)/(mr*mr)/FOUR;
  P_U_WS(3,1) = uu*(mr-ONE)*(mr+ONE)*(uu*uu*g-uu*uu+vv*vv*g-vv*vv+TWO*c*c)/(c*c)/
                (mr*mr)/TWO;
  P_U_WS(3,2) = vv*(mr-ONE)*(mr+ONE)*(uu*uu*g-uu*uu+vv*vv*g-vv*vv+TWO*c*c)/(c*c)/(mr*mr)/TWO;
  P_U_WS(3,3) = -(uu*uu*g*mr*mr-uu*uu*g-uu*uu*mr*mr+uu*uu+vv*vv*g*mr*mr-
                  vv*vv*g-vv*vv*mr*mr+vv*vv-TWO*c*c)/(c*c)/(mr*mr)/TWO;
}

inline void P_U_WS(DenseMatrix &P_U_WS, const Euler2D_pState &W) {
  double uu = W.v.x, vv = W.v.y, c = W.a(), mr = W.Mr();
  P_U_WS(0,0) = -(uu*uu*W.g*mr*mr-uu*uu*W.g-uu*uu*mr*mr+uu*uu+vv*vv*W.g*mr*mr-
                  vv*vv*W.g-vv*vv*mr*mr+vv*vv-TWO*c*c*mr*mr)/(c*c)/(mr*mr)/TWO;
  P_U_WS(0,1) = uu*(mr-ONE)*(mr+ONE)*W.gm1/(c*c)/(mr*mr);
  P_U_WS(0,2) = vv*(mr-ONE)*(mr+ONE)*W.gm1/(c*c)/(mr*mr);
  P_U_WS(0,3) = -(mr-ONE)*(mr+ONE)*W.gm1/(c*c)/(mr*mr);
  P_U_WS(1,0) = -uu*(uu*uu+vv*vv)*(mr-ONE)*(mr+ONE)*W.gm1/(c*c)/(mr*mr)/TWO;
  P_U_WS(1,1) = (uu*uu*W.g*mr*mr-uu*uu*W.g-uu*uu*mr*mr+uu*uu+c*c*mr*mr)/(c*c)/(mr*mr);
  P_U_WS(1,2) = uu*vv*(mr-ONE)*(mr+ONE)*W.gm1/(c*c)/(mr*mr);
  P_U_WS(1,3) = -uu*(mr-ONE)*(mr+ONE)*W.gm1/(c*c)/(mr*mr);
  P_U_WS(2,0) = -vv*(uu*uu+vv*vv)*(mr-ONE)*(mr+ONE)*W.gm1/(c*c)/(mr*mr)/TWO;
  P_U_WS(2,1) = uu*vv*(mr-ONE)*(mr+ONE)*W.gm1/(c*c)/(mr*mr);
  P_U_WS(2,2) = (vv*vv*W.g*mr*mr-vv*vv*W.g-vv*vv*mr*mr+vv*vv+c*c*mr*mr)/(c*c)/(mr*mr);
  P_U_WS(2,3) = -vv*(mr-ONE)*(mr+ONE)*W.gm1/(c*c)/(mr*mr);
  P_U_WS(3,0) = -(uu*uu+vv*vv)*(mr-ONE)*(mr+ONE)*(uu*uu*W.g-uu*uu+vv*vv*W.g-vv*vv+TWO*c*c)/
                 (c*c)/(mr*mr)/FOUR;
  P_U_WS(3,1) = uu*(mr-ONE)*(mr+ONE)*(uu*uu*W.g-uu*uu+vv*vv*W.g-vv*vv+TWO*c*c)/(c*c)/
                (mr*mr)/TWO;
  P_U_WS(3,2) = vv*(mr-ONE)*(mr+ONE)*(uu*uu*W.g-uu*uu+vv*vv*W.g-vv*vv+TWO*c*c)/(c*c)/(mr*mr)/TWO;
  P_U_WS(3,3) = -(uu*uu*W.g*mr*mr-uu*uu*W.g-uu*uu*mr*mr+uu*uu+vv*vv*W.g*mr*mr-
                  vv*vv*W.g-vv*vv*mr*mr+vv*vv-TWO*c*c)/(c*c)/(mr*mr)/TWO;
}

/**********************************************************
 * Euler2D_pState::P_U_WS_inv -- Inverse of Weiss-Smith   *
 * low-Mach-number local preconditioner.                  *
 **********************************************************/
inline void Euler2D_pState::P_U_WS_inv(DenseMatrix &P_U_WS_inv) {
  double uu = v.x, vv = v.y, c = a(), mr = Mr();
  P_U_WS_inv(0,0) = (TWO*c*c+vv*vv*g*mr*mr-vv*vv*g+uu*uu*g*mr*mr-
                      uu*uu*g-uu*uu*mr*mr+uu*uu-vv*vv*mr*mr+vv*vv)/(c*c)/TWO;
  P_U_WS_inv(0,1) = -uu*(mr-ONE)*(mr+ONE)*gm1/(c*c);
  P_U_WS_inv(0,2) = -vv*(mr-ONE)*(mr+ONE)*gm1/(c*c);
  P_U_WS_inv(0,3) = (mr-ONE)*(mr+ONE)*gm1/(c*c);
  P_U_WS_inv(1,0) = uu*(uu*uu+vv*vv)*(mr-ONE)*(mr+ONE)*gm1/(c*c)/TWO;
  P_U_WS_inv(1,1) = -(-c*c+uu*uu*g*mr*mr-uu*uu*g-uu*uu*mr*mr+uu*uu)/(c*c);
  P_U_WS_inv(1,2) = -uu*vv*(mr-ONE)*(mr+ONE)*gm1/(c*c);
  P_U_WS_inv(1,3) = uu*(mr-ONE)*(mr+ONE)*gm1/(c*c);
  P_U_WS_inv(2,0) = vv*(uu*uu+vv*vv)*(mr-ONE)*(mr+ONE)*gm1/(c*c)/TWO;
  P_U_WS_inv(2,1) = -uu*vv*(mr-ONE)*(mr+ONE)*gm1/(c*c);
  P_U_WS_inv(2,2) = -(-c*c+vv*vv*g*mr*mr-vv*vv*g-vv*vv*mr*mr+vv*vv)/(c*c);
  P_U_WS_inv(2,3) = vv*(mr-ONE)*(mr+ONE)*gm1/(c*c);
  P_U_WS_inv(3,0) = (uu*uu+vv*vv)*(mr-ONE)*(mr+ONE)*(uu*uu*g-uu*uu+vv*vv*g-vv*vv+TWO*c*c)/
                     (c*c)/FOUR;
  P_U_WS_inv(3,1) = -uu*(mr-ONE)*(mr+ONE)*(uu*uu*g-uu*uu+vv*vv*g-vv*vv+TWO*c*c)/(c*c)/TWO;
  P_U_WS_inv(3,2) = -vv*(mr-ONE)*(mr+ONE)*(uu*uu*g-uu*uu+vv*vv*g-vv*vv+TWO*c*c)/(c*c)/TWO;
  P_U_WS_inv(3,3) = (TWO*c*c*mr*mr+vv*vv*g*mr*mr-vv*vv*g+uu*uu*g*mr*mr-
                      uu*uu*g-uu*uu*mr*mr+uu*uu-vv*vv*mr*mr+vv*vv)/(c*c)/TWO;
}

inline void Euler2D_pState::P_U_WS_inv(DenseMatrix &P_U_WS_inv) const {
  double uu = v.x, vv = v.y, c = a(), mr = Mr();
  P_U_WS_inv(0,0) = (TWO*c*c+vv*vv*g*mr*mr-vv*vv*g+uu*uu*g*mr*mr-
                      uu*uu*g-uu*uu*mr*mr+uu*uu-vv*vv*mr*mr+vv*vv)/(c*c)/TWO;
  P_U_WS_inv(0,1) = -uu*(mr-ONE)*(mr+ONE)*gm1/(c*c);
  P_U_WS_inv(0,2) = -vv*(mr-ONE)*(mr+ONE)*gm1/(c*c);
  P_U_WS_inv(0,3) = (mr-ONE)*(mr+ONE)*gm1/(c*c);
  P_U_WS_inv(1,0) = uu*(uu*uu+vv*vv)*(mr-ONE)*(mr+ONE)*gm1/(c*c)/TWO;
  P_U_WS_inv(1,1) = -(-c*c+uu*uu*g*mr*mr-uu*uu*g-uu*uu*mr*mr+uu*uu)/(c*c);
  P_U_WS_inv(1,2) = -uu*vv*(mr-ONE)*(mr+ONE)*gm1/(c*c);
  P_U_WS_inv(1,3) = uu*(mr-ONE)*(mr+ONE)*gm1/(c*c);
  P_U_WS_inv(2,0) = vv*(uu*uu+vv*vv)*(mr-ONE)*(mr+ONE)*gm1/(c*c)/TWO;
  P_U_WS_inv(2,1) = -uu*vv*(mr-ONE)*(mr+ONE)*gm1/(c*c);
  P_U_WS_inv(2,2) = -(-c*c+vv*vv*g*mr*mr-vv*vv*g-vv*vv*mr*mr+vv*vv)/(c*c);
  P_U_WS_inv(2,3) = vv*(mr-ONE)*(mr+ONE)*gm1/(c*c);
  P_U_WS_inv(3,0) = (uu*uu+vv*vv)*(mr-ONE)*(mr+ONE)*(uu*uu*g-uu*uu+vv*vv*g-vv*vv+TWO*c*c)/
                     (c*c)/FOUR;
  P_U_WS_inv(3,1) = -uu*(mr-ONE)*(mr+ONE)*(uu*uu*g-uu*uu+vv*vv*g-vv*vv+TWO*c*c)/(c*c)/TWO;
  P_U_WS_inv(3,2) = -vv*(mr-ONE)*(mr+ONE)*(uu*uu*g-uu*uu+vv*vv*g-vv*vv+TWO*c*c)/(c*c)/TWO;
  P_U_WS_inv(3,3) = (TWO*c*c*mr*mr+vv*vv*g*mr*mr-vv*vv*g+uu*uu*g*mr*mr-
                      uu*uu*g-uu*uu*mr*mr+uu*uu-vv*vv*mr*mr+vv*vv)/(c*c)/TWO;
}

inline void P_U_WS_inv(DenseMatrix &P_U_WS_inv, const Euler2D_pState &W) {
  double uu = W.v.x, vv = W.v.y, c = W.a(), mr = W.Mr();
  P_U_WS_inv(0,0) = (TWO*c*c+vv*vv*W.g*mr*mr-vv*vv*W.g+uu*uu*W.g*mr*mr-
                      uu*uu*W.g-uu*uu*mr*mr+uu*uu-vv*vv*mr*mr+vv*vv)/(c*c)/TWO;
  P_U_WS_inv(0,1) = -uu*(mr-ONE)*(mr+ONE)*W.gm1/(c*c);
  P_U_WS_inv(0,2) = -vv*(mr-ONE)*(mr+ONE)*W.gm1/(c*c);
  P_U_WS_inv(0,3) = (mr-ONE)*(mr+ONE)*W.gm1/(c*c);
  P_U_WS_inv(1,0) = uu*(uu*uu+vv*vv)*(mr-ONE)*(mr+ONE)*W.gm1/(c*c)/TWO;
  P_U_WS_inv(1,1) = -(-c*c+uu*uu*W.g*mr*mr-uu*uu*W.g-uu*uu*mr*mr+uu*uu)/(c*c);
  P_U_WS_inv(1,2) = -uu*vv*(mr-ONE)*(mr+ONE)*W.gm1/(c*c);
  P_U_WS_inv(1,3) = uu*(mr-ONE)*(mr+ONE)*W.gm1/(c*c);
  P_U_WS_inv(2,0) = vv*(uu*uu+vv*vv)*(mr-ONE)*(mr+ONE)*W.gm1/(c*c)/TWO;
  P_U_WS_inv(2,1) = -uu*vv*(mr-ONE)*(mr+ONE)*W.gm1/(c*c);
  P_U_WS_inv(2,2) = -(-c*c+vv*vv*W.g*mr*mr-vv*vv*W.g-vv*vv*mr*mr+vv*vv)/(c*c);
  P_U_WS_inv(2,3) = vv*(mr-ONE)*(mr+ONE)*W.gm1/(c*c);
  P_U_WS_inv(3,0) = (uu*uu+vv*vv)*(mr-ONE)*(mr+ONE)*(uu*uu*W.g-uu*uu+vv*vv*W.g-vv*vv+TWO*c*c)/
                     (c*c)/FOUR;
  P_U_WS_inv(3,1) = -uu*(mr-ONE)*(mr+ONE)*(uu*uu*W.g-uu*uu+vv*vv*W.g-vv*vv+TWO*c*c)/(c*c)/TWO;
  P_U_WS_inv(3,2) = -vv*(mr-ONE)*(mr+ONE)*(uu*uu*W.g-uu*uu+vv*vv*W.g-vv*vv+TWO*c*c)/(c*c)/TWO;
  P_U_WS_inv(3,3) = (TWO*c*c*mr*mr+vv*vv*W.g*mr*mr-vv*vv*W.g+uu*uu*W.g*mr*mr-
                      uu*uu*W.g-uu*uu*mr*mr+uu*uu-vv*vv*mr*mr+vv*vv)/(c*c)/TWO;
}

/********************************************************
 * Euler2D_cState::Euler2D_cState -- Constructor.       *
 ********************************************************/
inline Euler2D_cState::Euler2D_cState(const Euler2D_pState &W) {
  d = W.d; dv = W.dv(); E = W.E();
}

/********************************************************
 * Euler2D_cState::W -- Primitive solution state.       *
 ********************************************************/
inline Euler2D_pState Euler2D_cState::W(void) {
  return (Euler2D_pState(d, v(), p()));
}

inline Euler2D_pState Euler2D_cState::W(void) const {
  return (Euler2D_pState(d, v(), p()));
}

inline Euler2D_pState Euler2D_cState::W(const Euler2D_cState &U) {
  return (Euler2D_pState(U.d, U.v(), U.p()));
}

inline Euler2D_pState W(const Euler2D_cState &U) {
  return (Euler2D_pState(U.d, U.v(), U.p()));
}

/********************************************************
 * Euler2D_cState::F -- Solution flux (x-direction).    *
 ********************************************************/
inline Euler2D_cState Euler2D_cState::F(void) {
  return (Euler2D_cState(dv.x, sqr(dv.x)/d + p(),
                         dv.x*dv.y/d, dv.x*H()/d));
}

inline Euler2D_cState Euler2D_cState::F(void) const {
  return (Euler2D_cState(dv.x, sqr(dv.x)/d + p(),
                         dv.x*dv.y/d, dv.x*H()/d));
}

inline Euler2D_cState Euler2D_cState::F(const Euler2D_cState &U) {
  return (Euler2D_cState(U.dv.x, sqr(U.dv.x)/U.d + U.p(),
                         U.dv.x*U.dv.y/U.d, U.dv.x*U.H()/U.d));
}

inline Euler2D_cState F(const Euler2D_cState &U) {
  return (Euler2D_cState(U.dv.x, sqr(U.dv.x)/U.d + U.p(),
                         U.dv.x*U.dv.y/U.d, U.dv.x*U.H()/U.d));
}

inline Euler2D_cState Euler2D_cState::F(const Vector2D &V) const {
  double vx = v().x, pp = p();
  return Euler2D_cState(d*(vx-V.x),
			(vx-V.x)*dv.x + pp,
			(vx-V.x)*dv.y,
			(vx-V.x)*E + vx*pp);
}

inline void Euler2D_cState::dFdU(DenseMatrix &dFdU) {
  W().dFdU(dFdU);
}

inline void Euler2D_cState::dFdU(DenseMatrix &dFdU) const {
  W().dFdU(dFdU);
}

inline void dFdU(DenseMatrix &dFdU, const Euler2D_cState &U) {
  U.W().dFdU(dFdU);
}

/********************************************************
 * Euler2D_cState::Fx -- Solution flux (x-direction).   *
 ********************************************************/
inline Euler2D_cState Euler2D_cState::Fx(void) {
  return (Euler2D_cState(dv.x, sqr(dv.x)/d + p(),
                         dv.x*dv.y/d, dv.x*H()/d));
}

inline Euler2D_cState Euler2D_cState::Fx(void) const {
  return (Euler2D_cState(dv.x, sqr(dv.x)/d + p(),
                         dv.x*dv.y/d, dv.x*H()/d));
}

inline Euler2D_cState Euler2D_cState::Fx(const Euler2D_cState &U) {
  return (Euler2D_cState(U.dv.x, sqr(U.dv.x)/U.d + U.p(),
                         U.dv.x*U.dv.y/U.d, U.dv.x*U.H()/U.d));
}

inline Euler2D_cState Fx(const Euler2D_cState &U) {
  return (Euler2D_cState(U.dv.x, sqr(U.dv.x)/U.d + U.p(),
                         U.dv.x*U.dv.y/U.d, U.dv.x*U.H()/U.d));
}

inline void Euler2D_cState::dFxdU(DenseMatrix &dFxdU) {
  W().dFxdU(dFxdU);
}

inline void Euler2D_cState::dFxdU(DenseMatrix &dFxdU) const {
  W().dFxdU(dFxdU);
}

inline void dFxdU(DenseMatrix &dFxdU, const Euler2D_cState &U) {
  U.W().dFxdU(dFxdU);
}

/********************************************************
 * Euler2D_cState::Fy -- Solution flux (y-direction).   *
 ********************************************************/
inline Euler2D_cState Euler2D_cState::Fy(void) {
  return (Euler2D_cState(dv.y, dv.x*dv.y/d,
			 sqr(dv.y)/d + p(), dv.y*H()/d));
}

inline Euler2D_cState Euler2D_cState::Fy(void) const {
  return (Euler2D_cState(dv.y, dv.x*dv.y/d,
			 sqr(dv.y)/d + p(), dv.y*H()/d));
}

inline Euler2D_cState Euler2D_cState::Fy(const Euler2D_cState &U) {
  return (Euler2D_cState(U.dv.y, U.dv.x*U.dv.y/U.d,
			 sqr(U.dv.y)/U.d + U.p(), U.dv.y*U.H()/U.d));
}

inline Euler2D_cState Fy(const Euler2D_cState &U) {
  return (Euler2D_cState(U.dv.y, U.dv.x*U.dv.y/U.d,
			 sqr(U.dv.y)/U.d + U.p(), U.dv.y*U.H()/U.d));
}

inline void Euler2D_cState::dFydU(DenseMatrix &dFydU) {
  W().dFydU(dFydU);
}

inline void Euler2D_cState::dFydU(DenseMatrix &dFydU) const {
  W().dFydU(dFydU);
}

inline void dFydU(DenseMatrix &dFydU, const Euler2D_cState &U) {
  U.W().dFydU(dFydU);
}

/********************************************************
 * Euler2D_cState::Fn -- Solution flux (n-direction).   *
 ********************************************************/
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

/**********************************************************
 * Euler2D_cState::S -- Source terms (axisymmetric flow). *
 **********************************************************/
inline Euler2D_cState Euler2D_cState::S(const Vector2D &X) {
  return (Euler2D_cState(-dv.y/X.y, -dv.x*dv.y/(d*X.y), 
                         -sqr(dv.y)/(d*X.y), -dv.y*H()/(d*X.y)));
}

inline Euler2D_cState Euler2D_cState::S(const Vector2D &X) const {
  return (Euler2D_cState(-dv.y/X.y, -dv.x*dv.y/(d*X.y), 
                         -sqr(dv.y)/(d*X.y), -dv.y*H()/(d*X.y)));
}

inline Euler2D_cState S(const Euler2D_cState &U, const Vector2D &X) {
  return (Euler2D_cState(-U.dv.y/X.y, -U.dv.x*U.dv.y/(U.d*X.y), 
                         -sqr(U.dv.y)/(U.d*X.y), -U.dv.y*U.H()/(U.d*X.y)));
}

/**********************************************************
 * Euler2D_cState::dUdW -- Solution Jacobian.             *
 **********************************************************/
inline void Euler2D_cState::dUdW(DenseMatrix &dUdW) {
  W().dUdW(dUdW);
}

inline void Euler2D_cState::dUdW(DenseMatrix &dUdW) const {
  W().dUdW(dUdW);
}

inline void dUdW(DenseMatrix &dUdW, const Euler2D_cState &U) {
  U.W().dUdW(dUdW);
}

/**********************************************************
 * Euler2D_cState::dWdU -- Solution Jacobian.             *
 **********************************************************/
inline void Euler2D_cState::dWdU(DenseMatrix &dWdU) {
  W().dWdU(dWdU);
}

inline void Euler2D_cState::dWdU(DenseMatrix &dWdU) const {
  W().dWdU(dWdU);
}

inline void dWdU(DenseMatrix &dWdU, const Euler2D_cState &U) {
  U.W().dWdU(dWdU);
}

/**********************************************************
 * Euler2D_cState::P_U_WS -- Weiss-Smith low-Mach-number  *
 * local preconditioner.                                  *
 **********************************************************/
inline void Euler2D_cState::P_U_WS(DenseMatrix &P_U_WS) {
  W().P_U_WS(P_U_WS);
}

inline void Euler2D_cState::P_U_WS(DenseMatrix &P_U_WS) const {
  W().P_U_WS(P_U_WS);
}

inline void P_U_WS(DenseMatrix &P_U_WS, const Euler2D_cState &U) {
  U.W().P_U_WS(P_U_WS);
}

/**********************************************************
 * Euler2D_cState::P_U_WS_inv -- Inverse of Weiss-Smith   *
 * low-Mach-number local preconditioner.                  *
 **********************************************************/
inline void Euler2D_cState::P_U_WS_inv(DenseMatrix &P_U_WS_inv) {
  W().P_U_WS_inv(P_U_WS_inv);
}

inline void Euler2D_cState::P_U_WS_inv(DenseMatrix &P_U_WS_inv) const {
  W().P_U_WS_inv(P_U_WS_inv);
}

inline void P_U_WS_inv(DenseMatrix &P_U_WS_inv, const Euler2D_cState &U) {
  U.W().P_U_WS_inv(P_U_WS_inv);
}

/********************************************************
 * Useful 2D Euler state constants.                     *
 ********************************************************/
const Euler2D_pState Euler2D_W_STDATM(DENSITY_STDATM,
				      Vector2D_ZERO, PRESSURE_STDATM);
const Euler2D_pState Euler2D_W_VACUUM(ZERO, Vector2D_ZERO, ZERO);
const Euler2D_pState Euler2D_W_ONE(ONE, ONE, ONE, ONE);
const Euler2D_cState Euler2D_U_STDATM(Euler2D_W_STDATM);
const Euler2D_cState Euler2D_U_VACUUM(Euler2D_W_VACUUM);
const Euler2D_cState Euler2D_U_ONE(ONE, ONE, ONE, ONE);

/********************************************************
 * Euler2DState -- External subroutines.                *
 ********************************************************/

extern Euler2D_pState Riemann(const Euler2D_pState &Wl,
	      	              const Euler2D_pState &Wr);

extern Euler2D_pState Riemann_x(const Euler2D_pState &Wl,
	      	                const Euler2D_pState &Wr);

extern Euler2D_pState Riemann_y(const Euler2D_pState &Wl,
	      	                const Euler2D_pState &Wr);

extern Euler2D_pState RoeAverage(const Euler2D_pState &Wl,
	      	                 const Euler2D_pState &Wr);

extern Euler2D_pState Rotate(const Euler2D_pState &W,
	      	             const Vector2D &norm_dir);

extern Euler2D_cState Rotate(const Euler2D_cState &U,
	      	             const Vector2D &norm_dir);

extern Euler2D_pState Translate(const Euler2D_pState &W,
				const Vector2D &V);

extern Euler2D_pState Reflect(const Euler2D_pState &W,
	      	              const Vector2D &norm_dir);

extern Euler2D_pState Reflect(const Euler2D_pState &W,
	      	              const Vector2D &norm_dir,
			      const Vector2D &V);

extern Euler2D_pState NoSlip(const Euler2D_pState &W,
	      	             const Vector2D &norm_dir);

extern Euler2D_pState BurningSurface(const Euler2D_pState &W,
				     const Vector2D &norm_dir);

extern Euler2D_pState MassInjection(const Euler2D_pState &W,
				    const Vector2D &norm_dir,
				    const int &bc_flag);

extern Euler2D_pState RinglebFlow(const Euler2D_pState &Wdum,
				  const Vector2D &X);

extern Euler2D_pState RinglebFlow(const Euler2D_pState &Wdum,
				  const Vector2D &X,
				  double &q, double &k);

extern Euler2D_pState RinglebFlowAverageState(const Euler2D_pState &Wdum,
					      const Vector2D &Y1,
					      const Vector2D &Y2,
					      const Vector2D &Y3,
					      const Vector2D &Y4);

extern Euler2D_pState BC_Characteristic(const Euler2D_pState &Wi,
                                        const Euler2D_pState &Wo,
	      	                        const Vector2D &norm_dir);

extern Euler2D_pState BC_Characteristic_Pressure(const Euler2D_pState &Wi,
                                                 const Euler2D_pState &Wo,
	      	                                 const Vector2D &norm_dir);

extern Euler2D_pState BC_Characteristic_Mach_Number(const Euler2D_pState &Wi,
                                                    const Euler2D_pState &Wo,
	      	                                    const Vector2D &norm_dir);

extern Euler2D_pState BCs(const Euler2D_pState &Wb,
			  const Euler2D_pState &Wi,
			  const Euler2D_pState &dWdx,
			  const double &dx,
			  const int BC_type,
			  const int End_type);

extern Euler2D_pState BCs_x(const Euler2D_pState &Wb,
			    const Euler2D_pState &Wi,
			    const Euler2D_pState &dWdx,
			    const double &dx,
			    const int BC_type,
			    const int End_type);

extern Euler2D_pState BCs_y(const Euler2D_pState &Wb,
			    const Euler2D_pState &Wi,
			    const Euler2D_pState &dWdy,
			    const double &dy,
			    const int BC_type,
			    const int End_type);

extern Euler2D_pState WaveSpeedPos(const Euler2D_pState &lambda_a,
                                   const Euler2D_pState &lambda_l,
                                   const Euler2D_pState &lambda_r);

extern Euler2D_pState WaveSpeedNeg(const Euler2D_pState &lambda_a,
                                   const Euler2D_pState &lambda_l,
                                   const Euler2D_pState &lambda_r);

extern Euler2D_pState WaveSpeedAbs(const Euler2D_pState &lambda_a,
                                   const Euler2D_pState &lambda_l,
                                   const Euler2D_pState &lambda_r);

extern Euler2D_pState HartenFixPos(const Euler2D_pState &lambda_a,
                                   const Euler2D_pState &lambda_l,
                                   const Euler2D_pState &lambda_r);

extern Euler2D_pState HartenFixNeg(const Euler2D_pState &lambda_a,
                                   const Euler2D_pState &lambda_l,
                                   const Euler2D_pState &lambda_r);

extern Euler2D_pState HartenFixAbs(const Euler2D_pState &lambda_a,
                                   const Euler2D_pState &lambda_l,
                                   const Euler2D_pState &lambda_r);

extern Euler2D_cState FluxGodunov(const Euler2D_pState &Wl,
	      	                  const Euler2D_pState &Wr);

extern Euler2D_cState FluxGodunov(const Euler2D_cState &Ul,
	      	                  const Euler2D_cState &Ur);

extern Euler2D_cState FluxGodunov_x(const Euler2D_pState &Wl,
	      	                    const Euler2D_pState &Wr);

extern Euler2D_cState FluxGodunov_x(const Euler2D_cState &Ul,
	      	                    const Euler2D_cState &Ur);

extern Euler2D_cState FluxGodunov_y(const Euler2D_pState &Wl,
	      	                    const Euler2D_pState &Wr);

extern Euler2D_cState FluxGodunov_y(const Euler2D_cState &Ul,
	      	                    const Euler2D_cState &Ur);

extern Euler2D_cState FluxGodunov_n(const Euler2D_pState &Wl,
	      	                    const Euler2D_pState &Wr,
                                    const Vector2D &norm_dir);

extern Euler2D_cState FluxGodunov_n(const Euler2D_cState &Ul,
	      	                    const Euler2D_cState &Ur,
                                    const Vector2D &norm_dir);

extern Euler2D_cState FluxGodunov_MB_n(const Euler2D_pState &Wl,
				       const Euler2D_pState &Wr,
				       const Vector2D &V,
				       const Vector2D &norm_dir);

extern Euler2D_pState StateGodunov_n(const Euler2D_pState &Wl,
				     const Euler2D_pState &Wr,
				     const Vector2D &norm_dir);

extern Euler2D_cState FluxRoe(const Euler2D_pState &Wl,
	      	              const Euler2D_pState &Wr);

extern Euler2D_cState FluxRoe(const Euler2D_cState &Ul,
	      	              const Euler2D_cState &Ur);

extern Euler2D_cState FluxRoe_x(const Euler2D_pState &Wl,
	      	                const Euler2D_pState &Wr);

extern Euler2D_cState FluxRoe_x(const Euler2D_cState &Ul,
	      	                const Euler2D_cState &Ur);

extern Euler2D_cState FluxRoe_y(const Euler2D_pState &Wl,
	      	                const Euler2D_pState &Wr);

extern Euler2D_cState FluxRoe_y(const Euler2D_cState &Ul,
	      	                const Euler2D_cState &Ur);

extern Euler2D_cState FluxRoe_n(const Euler2D_pState &Wl,
	      	                const Euler2D_pState &Wr,
                                const Vector2D &norm_dir);

extern Euler2D_cState FluxRoe_n(const Euler2D_cState &Ul,
	      	                const Euler2D_cState &Ur,
                                const Vector2D &norm_dir);

extern Euler2D_cState FluxRoe_MB(const Euler2D_pState &Wl,
				 const Euler2D_pState &Wr,
				 const Vector2D &V);

extern Euler2D_cState FluxRoe_MB(const Euler2D_cState &Ul,
				 const Euler2D_cState &Ur,
				 const Vector2D &V);

extern Euler2D_cState FluxRoe_MB_n(const Euler2D_pState &Wl,
				   const Euler2D_pState &Wr,
				   const Vector2D &V,
				   const Vector2D &norm_dir);

extern Euler2D_cState FluxRoe_MB_n(const Euler2D_cState &Ul,
				   const Euler2D_cState &Ur,
				   const Vector2D &V,
				   const Vector2D &norm_dir);

extern Euler2D_cState FluxRusanov(const Euler2D_pState &Wl,
	      	                  const Euler2D_pState &Wr);

extern Euler2D_cState FluxRusanov(const Euler2D_cState &Ul,
	      	                  const Euler2D_cState &Ur);

extern Euler2D_cState FluxRusanov_x(const Euler2D_pState &Wl,
	      	                    const Euler2D_pState &Wr);

extern Euler2D_cState FluxRusanov_x(const Euler2D_cState &Ul,
	      	                    const Euler2D_cState &Ur);

extern Euler2D_cState FluxRusanov_y(const Euler2D_pState &Wl,
	      	                    const Euler2D_pState &Wr);

extern Euler2D_cState FluxRusanov_y(const Euler2D_cState &Ul,
	      	                    const Euler2D_cState &Ur);

extern Euler2D_cState FluxRusanov_n(const Euler2D_pState &Wl,
	      	                    const Euler2D_pState &Wr,
                                    const Vector2D &norm_dir);

extern Euler2D_cState FluxRusanov_n(const Euler2D_cState &Ul,
	      	                    const Euler2D_cState &Ur,
                                    const Vector2D &norm_dir);

extern Euler2D_cState FluxHLLE(const Euler2D_pState &Wl,
	      	               const Euler2D_pState &Wr);

extern Euler2D_cState FluxHLLE(const Euler2D_cState &Ul,
	      	               const Euler2D_cState &Ur);

extern Euler2D_cState FluxHLLE_x(const Euler2D_pState &Wl,
	      	                 const Euler2D_pState &Wr);

extern Euler2D_cState FluxHLLE_x(const Euler2D_cState &Ul,
	      	                 const Euler2D_cState &Ur);

extern Euler2D_cState FluxHLLE_y(const Euler2D_pState &Wl,
	      	                 const Euler2D_pState &Wr);

extern Euler2D_cState FluxHLLE_y(const Euler2D_cState &Ul,
	      	                 const Euler2D_cState &Ur);
  
extern Euler2D_cState FluxHLLE_n(const Euler2D_pState &Wl,
	      	                 const Euler2D_pState &Wr,
                                 const Vector2D &norm_dir);

extern Euler2D_cState FluxHLLE_n(const Euler2D_cState &Ul,
	      	                 const Euler2D_cState &Ur,
                                 const Vector2D &norm_dir);

extern Euler2D_cState FluxHLLE_MB(const Euler2D_pState &Wl,
				  const Euler2D_pState &Wr,
				  const Vector2D &V);

extern Euler2D_cState FluxHLLE_MB(const Euler2D_cState &Ul,
				  const Euler2D_cState &Ur,
				  const Vector2D &V);

extern Euler2D_cState FluxHLLE_MB_n(const Euler2D_pState &Wl,
				    const Euler2D_pState &Wr,
				    const Vector2D &V,
				    const Vector2D &norm_dir);

extern Euler2D_cState FluxHLLE_MB_n(const Euler2D_cState &Ul,
				    const Euler2D_cState &Ur,
				    const Vector2D &V,
				    const Vector2D &norm_dir);

extern Vector2D HLLE_wavespeeds(const Euler2D_pState &Wl,
                                const Euler2D_pState &Wr,
                                const Vector2D &norm_dir);

extern Euler2D_cState FluxLinde(const Euler2D_pState &Wl,
	      	                const Euler2D_pState &Wr);

extern Euler2D_cState FluxLinde(const Euler2D_cState &Ul,
	      	                const Euler2D_cState &Ur);

extern Euler2D_cState FluxLinde_x(const Euler2D_pState &Wl,
	      	                  const Euler2D_pState &Wr);

extern Euler2D_cState FluxLinde_x(const Euler2D_cState &Ul,
	      	                  const Euler2D_cState &Ur);

extern Euler2D_cState FluxLinde_y(const Euler2D_pState &Wl,
	      	                  const Euler2D_pState &Wr);

extern Euler2D_cState FluxLinde_y(const Euler2D_cState &Ul,
	      	                  const Euler2D_cState &Ur);

extern Euler2D_cState FluxLinde_n(const Euler2D_pState &Wl,
	      	                  const Euler2D_pState &Wr,
                                  const Vector2D &norm_dir);

extern Euler2D_cState FluxLinde_n(const Euler2D_cState &Ul,
	      	                  const Euler2D_cState &Ur,
                                  const Vector2D &norm_dir);

extern Euler2D_cState FluxHLLC(const Euler2D_pState &Wl,
	      	               const Euler2D_pState &Wr);

extern Euler2D_cState FluxHLLC(const Euler2D_cState &Ul,
	      	               const Euler2D_cState &Ur);

extern Euler2D_cState FluxHLLC_x(const Euler2D_pState &Wl,
	      	                 const Euler2D_pState &Wr);

extern Euler2D_cState FluxHLLC_x(const Euler2D_cState &Ul,
	      	                 const Euler2D_cState &Ur);

extern Euler2D_cState FluxHLLC_y(const Euler2D_pState &Wl,
	      	                 const Euler2D_pState &Wr);

extern Euler2D_cState FluxHLLC_y(const Euler2D_cState &Ul,
	      	                 const Euler2D_cState &Ur);

extern Euler2D_cState FluxHLLC_n(const Euler2D_pState &Wl,
	      	                 const Euler2D_pState &Wr,
                                 const Vector2D &norm_dir);

extern Euler2D_cState FluxHLLC_n(const Euler2D_cState &Ul,
	      	                 const Euler2D_cState &Ur,
                                 const Vector2D &norm_dir);

extern Euler2D_cState FluxVanLeer(const Euler2D_pState &Wl,
				  const Euler2D_pState &Wr);

extern Euler2D_cState FluxVanLeer(const Euler2D_cState &Wl,
				  const Euler2D_cState &Wr);

extern Euler2D_cState FluxVanLeer_n(const Euler2D_pState &Wl,
				    const Euler2D_pState &Wr,
				    const Vector2D &norm_dir);

extern Euler2D_cState FluxVanLeer_n(const Euler2D_cState &Wl,
				    const Euler2D_cState &Wr,
				    const Vector2D &norm_dir);

extern Euler2D_cState FluxVanLeer_MB(const Euler2D_pState &Wl,
				     const Euler2D_pState &Wr,
				     const Vector2D &V);

extern Euler2D_cState FluxVanLeer_MB_n(const Euler2D_pState &Wl,
				       const Euler2D_pState &Wr,
				       const Vector2D &V,
				       const Vector2D &norm_dir);

extern Euler2D_cState FluxAUSM(const Euler2D_pState &Wl,
			       const Euler2D_pState &Wr);

extern Euler2D_cState FluxAUSM(const Euler2D_cState &Wl,
			       const Euler2D_cState &Wr);

extern Euler2D_cState FluxAUSM_n(const Euler2D_pState &Wl,
				 const Euler2D_pState &Wr,
				 const Vector2D &norm_dir);

extern Euler2D_cState FluxAUSM_n(const Euler2D_cState &Wl,
				 const Euler2D_cState &Wr,
				 const Vector2D &norm_dir);

extern Euler2D_cState FluxAUSMplus(const Euler2D_pState &Wl,
				   const Euler2D_pState &Wr);

extern Euler2D_cState FluxAUSMplus(const Euler2D_cState &Wl,
				   const Euler2D_cState &Wr);

extern Euler2D_cState FluxAUSMplus_n(const Euler2D_pState &Wl,
				     const Euler2D_pState &Wr,
				     const Vector2D &norm_dir);

extern Euler2D_cState FluxAUSMplus_n(const Euler2D_cState &Wl,
				     const Euler2D_cState &Wr,
				     const Vector2D &norm_dir);

extern Euler2D_cState FluxRoe_Precon_WS(const Euler2D_cState &Ul,
	      	                        const Euler2D_cState &Ur);

extern Euler2D_cState FluxRoe_Precon_WS(const Euler2D_cState &Ul,
	      	                        const Euler2D_cState &Ur);

extern Euler2D_cState FluxRoe_n_Precon_WS(const Euler2D_pState &Wl,
	      	                          const Euler2D_pState &Wr,
                                          const Vector2D &norm_dir);

extern Euler2D_cState FluxRoe_n_Precon_WS(const Euler2D_cState &Ul,
	      	                          const Euler2D_cState &Ur,
                                          const Vector2D &norm_dir);

extern Euler2D_cState FluxHLLE_Precon_WS(const Euler2D_pState &Wl,
	      	                         const Euler2D_pState &Wr);

extern Euler2D_cState FluxHLLE_Precon_WS(const Euler2D_cState &Ul,
	      	                         const Euler2D_cState &Ur);

extern Euler2D_cState FluxHLLE_n_Precon_WS(const Euler2D_pState &Wl,
	      	                           const Euler2D_pState &Wr,
                                           const Vector2D &norm_dir);

extern Euler2D_cState FluxHLLE_n_Precon_WS(const Euler2D_cState &Ul,
	      	                           const Euler2D_cState &Ur,
                                           const Vector2D &norm_dir);

#endif /* _EULER2D_STATE_INCLUDED  */
