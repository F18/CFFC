/**********************************************************************
 * Particle2DState.h: Header file defining 2D particle solution state *
 *                    classes.                                        *
 **********************************************************************/

#ifndef _PARTICLE2D_STATE_INCLUDED
#define _PARTICLE2D_STATE_INCLUDED

// Include required C++ libraries.

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

// Include math macro, CFD, 2D vector, and gas constant header files.

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

#ifndef _SOLID_CONSTANTS_INCLUDED
#include "../Physics/SolidConstants.h"
#endif // _SOLID_CONSTANTS_INCLUDED

#define	NUM_VAR_PARTICLE2D  4

// Define the classes.

class Particle2D_cState;

/*!
 * Class: Particle2D_pState
 *
 * @brief Primitive variable solution state class definition for an
 * Eulerian description of an inert, disperse and dilute particle-phase.
 * This class is used by Dusty2D_pState.
 *
 * \verbatim
 * Member functions
 *     sigma  -- Particle concentration.
 *     u      -- Particle velocity.
 *     Tp     -- Particle temperature.
 *     nd     -- Particle-phase number density.
 *     ep     -- Particle-phase total thermal/internal energy.
 *     Ep     -- Particle-phase total energy.
 *     du     -- Particle-phase momentum.
 *
 *     W      -- Return primitive solution state.
 *     F      -- Return x-direction solution flux.
 *     Sa     -- Return axisymmetric source term vector.
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
class Particle2D_pState {
 private:
 public:
  double    sigma; //!< Particle-phase concentration.   
  Vector2D      u; //!< Particle-phase velocity.
  double       Tp; //!< Particle-phase temperature.

  //@{ @name Constructors and desctructors.
  //! Creation constructor.
  Particle2D_pState(void) {
    sigma = ZERO; u.zero(); Tp = ZERO;
  }

  //! Copy constructor.
  Particle2D_pState(const Particle2D_pState &W) {
    sigma = W.sigma; u = W.u; Tp = W.Tp;
  }

  //! Copy constructor.
  Particle2D_pState(const Particle2D_cState &U, const double &cm);

  //! Assignment constructors.
  Particle2D_pState(const double &sig, const Vector2D &U, const double &Tpart) {
    sigma = sig; u = U; Tp = Tpart;
  }

  //! Assignment constructors.
  Particle2D_pState(const double &sig, const double &ux, const double &uy, const double &Tpart) {
    sigma = sig; u.x = ux; u.y = uy; Tp = Tpart;
  }

  //! Destructor.
  // ~Particle2D_pState(void);
  // Use automatically generated destructor.
  //@}

  //@{ @name State operators.
  //! Copy function.
  void Copy(const Particle2D_pState &W) {
    sigma = W.sigma; u = W.u; Tp = W.Tp;
  }

  //! Copy function.
  void Copy(const double &sig, const double &ux, const double &uy, const double &Tpart) {
    sigma = sig; u.x = ux; u.y = uy; Tp = Tpart;
  }

  //! Vacuum operator.
  void Vacuum(void) {
    sigma = ZERO; u = Vector2D_ZERO; Tp = ZERO;
  }

  //! Constant operator.
  void Constant(const double &val) {
    sigma = val; u = Vector2D(val,val); Tp = val;
  }

  //! Check for unphysical state properties.
  int Unphysical_Properties(void) const {
    if (sigma < ZERO || Tp < ZERO) return 1;
    return 0;
  }

  //! Reset unphysical state properties.
  void Reset_Unphysical_Properties(void) {
    if (fabs(sigma) < NANO || fabs(Tp) < NANO) Vacuum();
  }
  //@}

  //! Return the number of variables.
  int NumVar(void) { return NUM_VAR_PARTICLE2D; }

  //! Particle-phase number density.
  double nd(const double &mp) const;

  //! Particle-phase momentum.
  Vector2D du(void) const;

  //! Particle-phase total thermal/internal energy.
  double ep(const double &cm) const;
  
  //! Particle-phase total energy.
  double Ep(const double &cm) const;

  //! Particle-phase current.
  Vector2D jc(const double &qe, const double &mp) const;

  //@{ @name Conserved solution state.
  Particle2D_cState U(const double &cm) const;
  Particle2D_cState U(const Particle2D_pState &W, const double &cm) const;
  //@}

  //@{ @name Transformation matrices.
  //! Jacobian of the conserved solution variables with respect to the
  //! primitive solution variables.
  void dUdW(DenseMatrix &dUdW, const int &n, const double &cm) const;

  //! Jacobian of the primitive solution variables with respect to the
  //! conserved solution variables.
  void dWdU(DenseMatrix &dWdU, const int &n, const double &cm) const;
  //@}

  //@{ @name Solution flux and x-direction Jacobian.
  Particle2D_cState F(const double &cm) const;
  Particle2D_cState F(const Vector2D &V, const double &cm) const;
  void dFdU(DenseMatrix &dFdU, const int &n, const double &cm) const;
  void dFdU(DenseMatrix &dFdU, const int &n, const Vector2D &V, const double &cm) const;
  //@}

  //@{ @name Eigenstructure (x-direction).
  Particle2D_pState lambda_x(void) const;
  Particle2D_pState lambda_x(const Vector2D &V) const;

  //! Primitive right eigenvector (x-direction).
  Particle2D_pState rp_x(int index) const;

  //! Conserved right eigenvector (x-direction).
  Particle2D_cState rc_x(int index) const;

  //! Primitive left eigenvector (x-direction).
  //Particle2D_pState lp_x(int index) const;
  //@}

  //@{ @name Axisymmetric flow source term vector and Jacobian.
  Particle2D_cState Sa(const Vector2D &X, const double &cm) const;
  void dSadU(DenseMatrix &dSadU, const Vector2D &X, const int &n, const double &cm) const;
  //@}

  //@{ @name Electrostatic-force source term vector and Jacobian.
  Particle2D_cState Se(const Vector2D &E, const double &Cm) const;
  void dSedU(DenseMatrix &dSedU, const int &n, const Vector2D &E, const double &Cm) const;
  //@}

  //@{ @name State tranformations, boundary conditions, and flux functions.
  void Rotate(const Particle2D_pState &W, const Vector2D &norm_dir);

  void Translate(const Particle2D_pState &W, const Vector2D &V);

  void Reflect(const Particle2D_pState &W, const Vector2D &norm_dir);

  void Mirror(const Particle2D_pState &W, const Vector2D &norm_dir);

  void Absorb(const Particle2D_pState &W, const Vector2D &norm_dir);

  void FluxSaurel(const Particle2D_pState &Wl, const Particle2D_pState &Wr);
  //@}

  //@{ @name Index operator.
  double &operator[](int index) {
    //assert(index >= 1 && index <= NUM_VAR_PARTICLE2D);
    switch(index) {
    case 1 :
      return sigma;
    case 2 :
      return u.x;
    case 3 :
      return u.y;
    case 4 :
      return Tp;
    default:
      return sigma;
    };
  }  
  const double &operator[](int index) const {
    //assert(index >= 1 && index <= NUM_VAR_PARTICLE2D);
    switch(index) {
    case 1 :
      return sigma;
    case 2 :
      return u.x;
    case 3 :
      return u.y;
    case 4 :
      return Tp;
    default:
      return sigma;
    };
  }
  //@}

  //@{ Binary arithmetic operators.
  Particle2D_pState operator +(const Particle2D_pState &W) const;
  Particle2D_pState operator -(const Particle2D_pState &W) const;
  double operator *(const Particle2D_pState &W) const;
  Particle2D_pState operator *(const double &a) const;
  friend Particle2D_pState operator *(const double &a, const Particle2D_pState &W);
  Particle2D_pState operator /(const double &a) const;
  Particle2D_pState operator ^(const Particle2D_pState &W) const;
  //@}

  //! Assignment operator.
  Particle2D_pState &operator =(const Particle2D_pState &W);

  //@{ Unary arithmetic operators.
  //Particle2D_pState operator +(const Particle2D_pState &W);
  Particle2D_pState operator -(const Particle2D_pState &W);
  //@}

  //@{ Shortcut arithmetic operators.
  Particle2D_pState &operator +=(const Particle2D_pState &W);
  Particle2D_pState &operator -=(const Particle2D_pState &W);
  Particle2D_pState &operator *=(const double &a);
  Particle2D_pState &operator /=(const double &a);
  //@}

  //@{ Relational operators.
  friend int operator ==(const Particle2D_pState &W1, const Particle2D_pState &W2);
  friend int operator !=(const Particle2D_pState &W1, const Particle2D_pState &W2);
  //@}

  //@{ Input-output operators.
  friend ostream &operator << (ostream &out_file, const Particle2D_pState &W);
  friend istream &operator >> (istream &in_file,  Particle2D_pState &W);
  //@}

};

/*!
 * Class: Particle2D_cState
 *
 * @brief Conserved variable solution state class definition for an
 * Eulerian description of an inert, disperse and dilute particle-phase.
 * This class is used by Dusty2D_pState.
 *
 * \verbatim
 * Member functions
 *     sigma  -- Particle concentration.
 *     du     -- Particle momentum.
 *     ep     -- Particle-phase total thermal/internal energy.
 *     nd     -- Particle-phase number density.
 *     u      -- Particle-phase velocity.
 *     Tp     -- Particle-phase temperature.
 *     Ep     -- Particle total energy.
 *
 *     W      -- Return primitive solution state.
 *     F      -- Return x-direction solution flux.
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
class Particle2D_cState {
  private:
  public:
  double    sigma; //!< Particle-phase concentration.
  Vector2D     du; //!< Particle-phase momentum.
  double       ep; //!< Particle-phase total thermal/internal energy.

  //@{ @name Constructors and desctructors.
  //! Creation constructor.
  Particle2D_cState(void) {
    sigma = ZERO; du.zero(); ep = ZERO;
  }

  //! Copy constructor.
  Particle2D_cState(const Particle2D_cState &U) {
    sigma = U.sigma; du = U.du; ep = U.ep;
  }

  //! Copy constructor.
  Particle2D_cState(const Particle2D_pState &W, const double &cm);
    
  //! Assignment constructors.
  Particle2D_cState(const double &sig,
		    const Vector2D &sigU,
		    const double &epart) {
    sigma = sig; du = sigU; ep = epart;
  }
  
  //! Assignment constructors.
  Particle2D_cState(const double &sig, const double &sigux, const double &siguy, const double &epart) {
    sigma = sig; du.x = sigux; du.y = siguy; ep = epart;
  }

  //! Destructor.
  // ~Particle2D_cState(void);
  // Use automatically generated destructor.
  //@}

  //@{ @name State operators.
  //! Copy function.
  void Copy(const Particle2D_cState &U) {
    sigma = U.sigma; du = U.du; ep = U.ep;
  }
  
  //! Copy function.
  void Copy(const double &sig, const double &sigux, const double &siguy, const double &epart) {
    sigma = sig; du.x = sigux; du.y = siguy; ep = epart;
  }

  //! Vacuum operator.
  void Vacuum(void) {
    sigma = ZERO; du = Vector2D_ZERO; ep = ZERO;
  }

  //! Constant operator.
  void Constant(const double &val) {
    sigma = val; du = Vector2D(val,val); ep = val;
  }

  //! Check for unphysical state properties.
  int Unphysical_Properties(void) const {
    if (sigma < ZERO || ep < ZERO) return 1;
    return 0;
  }

  //! Reset unphysical state properties.
  void Reset_Unphysical_Properties(void) {
    if (fabs(sigma) < NANO || ep < fabs(NANO)) Vacuum();
  }
  //@}

  //! Return the number of variables.
  int NumVar(void) { return NUM_VAR_PARTICLE2D; }

  //! Particle-phase number density.
  double nd(const double &mp) const;

  //! Particle-phase velocity.
  Vector2D u(void) const;

  //! Particle-phase temperature.
  double Tp(const double &cm) const;

  //! Particle-phase total energy.
  double Ep(void) const;

  //! Particle-phase current.
  Vector2D jc(const double &qe, const double &mp) const;

  //@{ @name Primitive solution state.
  Particle2D_pState W(const double &cm) const;
  Particle2D_pState W(const Particle2D_cState &U, const double &cm) const;
  //@}

  //! Jacobian of the conserved solution variables with respect to the
  //! primitive solution variables.
  void dUdW(DenseMatrix &dUdW, const int &n, const double &cm) const;

  //! Jacobian of the primitive solution variables with respect to the
  //! conserved solution variables.
  void dWdU(DenseMatrix &dWdU, const int &n, const double &cm) const;

  //@{ @name Solution flux and x-direction Jacobian.
  Particle2D_cState F(void) const;
  Particle2D_cState F(const Vector2D &V) const;
  void F(const Particle2D_pState &W, const double &cm);
  void F(const Particle2D_pState &W, const Vector2D &V, const double &cm);
  void dFdU(DenseMatrix &dFdU, const int &n) const;
  void dFdU(DenseMatrix &dFdU, const int &n, const Vector2D &V) const;
  //@}

  //@{ @name Eigenstructure (x-direction).
  Particle2D_pState lambda_x(void) const;
  Particle2D_pState lambda_x(const Vector2D &V) const;

  //! Primitive right eigenvector (x-direction).
  Particle2D_pState rp_x(int index) const;

  //! Conserved right eigenvector (x-direction).
  Particle2D_cState rc_x(int index) const;

  //! Primitive left eigenvector (x-direction).
  //Particle2D_pState lp_x(int index) const;
  //@}

  //@{ @name State tranformations, boundary conditions, and flux functions.
  void Rotate(const Particle2D_cState &U, const Vector2D &norm_dir);

  void FluxSaurel_n(const Particle2D_pState &Wl, const Particle2D_pState &Wr,
		    const Vector2D &norm_dir, const double &cm);

  void FluxSaurel_MB_n(const Particle2D_pState &Wl, const Particle2D_pState &Wr,
		       const Vector2D &V, const Vector2D &norm_dir, const double &cm);
  //@}

  //@{ @name Index operator.
  double &operator[](int index) {
    //assert(index >= 1 && index <= NUM_VAR_PARTICLE2D);
    switch(index) {
    case 1 :
      return sigma;
    case 2 :
      return du.x;
    case 3 :
      return du.y;
    case 4 :
      return ep;
    default:
      return sigma;
    };
  }    
  const double &operator[](int index) const {
    //assert(index >= 1 && index <= NUM_VAR_PARTICLE2D);
    switch(index) {
    case 1 :
      return sigma;
    case 2 :
      return du.x;
    case 3 :
      return du.y;
    case 4 :
      return ep;
    default:
      return sigma;
    };
  }
  //@}

  //@{ Binary arithmetic operators.
  Particle2D_cState operator +(const Particle2D_cState &U) const;
  Particle2D_cState operator -(const Particle2D_cState &U) const;
  double operator *(const Particle2D_cState &U) const;
  Particle2D_cState operator *(const double &a) const;
  friend Particle2D_cState operator *(const double &a, const Particle2D_cState &U);
  Particle2D_cState operator /(const double &a) const;
  Particle2D_cState operator ^(const Particle2D_cState &U) const;
  //@}

  // Assignment operator.
  Particle2D_cState &operator =(const Particle2D_cState &U);
  
  //@{ Unary arithmetic operators.
  //Particle2D_cState operator +(const Particle2D_cState &U);
  Particle2D_cState operator -(const Particle2D_cState &U);
  //@}

  //@{ Shortcut arithmetic operators.
  Particle2D_cState &operator +=(const Particle2D_cState &U);
  Particle2D_cState &operator -=(const Particle2D_cState &U);
  Particle2D_cState &operator *=(const double &a);
  Particle2D_cState &operator /=(const double &a);
  //@}

  //@{ Relational operators.
  friend int operator ==(const Particle2D_cState &U1, const Particle2D_cState &U2);
  friend int operator !=(const Particle2D_cState &U1, const Particle2D_cState &U2);
  //@}

  //@{ Input-output operators.
  friend ostream &operator << (ostream &out_file, const Particle2D_cState &U);
  friend istream &operator >> (istream &in_file, Particle2D_cState &U);
  //@}

};

/**********************************************************************
 * Particle2D_pState::nd -- Particle-phase number density.            *
 **********************************************************************/
inline double Particle2D_pState::nd(const double &mp) const {
  return sigma/mp;
}

/**********************************************************************
 * Particle2D_pState::du -- Particle-phase momentum.                  *
 **********************************************************************/
inline Vector2D Particle2D_pState::du(void) const {
  return sigma*u;
}

/**********************************************************************
 * Particle2D_pState::ep -- Particle-phase total thermal/internal     *
 *                       -- energy.                                   *
 **********************************************************************/
inline double Particle2D_pState::ep(const double &cm) const {
  return sigma*cm*Tp;
}

/**********************************************************************
 * Particle2D_pState::Ep -- Particle-phase total energy.              *
 **********************************************************************/
inline double Particle2D_pState::Ep(const double &cm) const {
  return sigma*cm*Tp + HALF*sigma*u.sqr();
}

/**********************************************************************
 * Particle2D_pState::jc -- Particle-phase current.                   *
 **********************************************************************/
inline Vector2D Particle2D_pState::jc(const double &qe, const double &mp) const {
  return qe*nd(mp)*u;
}

/**********************************************************************
 * Particle2D_pState -- Binary arithmetic operators.                  *
 **********************************************************************/
inline Particle2D_pState Particle2D_pState::operator +(const Particle2D_pState &W) const {
  return Particle2D_pState(sigma+W.sigma,u.x+W.u.x,u.y+W.u.y,Tp+W.Tp);
}

inline Particle2D_pState Particle2D_pState::operator -(const Particle2D_pState &W) const {
  return Particle2D_pState(sigma-W.sigma,u.x-W.u.x,u.y-W.u.y,Tp-W.Tp);
}

// Inner product operator.
inline double Particle2D_pState::operator *(const Particle2D_pState &W) const {
  return (sigma*W.sigma + u.x*W.u.x + u.y*W.u.y + Tp*W.Tp);
}

inline Particle2D_pState Particle2D_pState::operator *(const double &a) const {
  return Particle2D_pState(sigma*a,u.x*a,u.y*a,Tp*a);
}

inline Particle2D_pState operator *(const double &a, const Particle2D_pState &W) {
  return Particle2D_pState(W.sigma*a,W.u.x*a,W.u.y*a,W.Tp*a);
}

inline Particle2D_pState Particle2D_pState::operator /(const double &a) const {
  return Particle2D_pState(sigma/a,u.x/a,u.y/a,Tp/a);
}

// A useful solution state product operator.
inline Particle2D_pState Particle2D_pState::operator ^(const Particle2D_pState &W) const {
  return Particle2D_pState(sigma*W.sigma,u.x*W.u.x,u.y*W.u.y,Tp*W.Tp);
}

/**********************************************************************
 * Particle2D_pState -- Assignment operator.                          *
 **********************************************************************/
inline Particle2D_pState& Particle2D_pState::operator =(const Particle2D_pState &W) {
  //if (this != &W) {
  sigma = W.sigma; u.x = W.u.x; u.y = W.u.y; Tp = W.Tp;
  //}
  return *this;
}

/**********************************************************************
 * Particle2D_pState -- Unary arithmetic operators.                   *
 **********************************************************************/
//inline Particle2D_pState operator +(const Particle2D_pState &W) {
//return Particle2D_pState(W.sigma,W.u.x,W.u.y,W.Tp);
//}

inline Particle2D_pState operator -(const Particle2D_pState &W) {
  return Particle2D_pState(-W.sigma,-W.u.x,-W.u.y,-W.Tp);
}

/**********************************************************************
 * Particle2D_pState -- Shortcut arithmetic operators.                *
 **********************************************************************/
inline Particle2D_pState& Particle2D_pState::operator +=(const Particle2D_pState &W) {
  sigma += W.sigma; u.x += W.u.x; u.y += W.u.y; Tp += W.Tp;
  return *this;
}

inline Particle2D_pState& Particle2D_pState::operator -=(const Particle2D_pState &W) {
  sigma -= W.sigma; u.x -= W.u.x; u.y -= W.u.y; Tp -= W.Tp;
  return *this;
}

inline Particle2D_pState& Particle2D_pState::operator *=(const double &a) {
  sigma *= a; u.x *= a; u.y *= a; Tp *= a;
  return *this;
}

inline Particle2D_pState& Particle2D_pState::operator /=(const double &a) {
  sigma /= a; u.x /= a; u.y /= a; Tp /= a;
  return *this;
}

/**********************************************************************
 * Particle2D_pState -- Relational operators.                         *
 **********************************************************************/
inline int operator ==(const Particle2D_pState &W1, const Particle2D_pState &W2) {
  return (W1.sigma == W2.sigma && W1.u == W2.u && W1.Tp == W2.Tp);
}

inline int operator !=(const Particle2D_pState &W1, const Particle2D_pState &W2) {
  return (W1.sigma != W2.sigma || W1.u != W2.u || W1.Tp != W2.Tp);
}

/**********************************************************************
 * Particle2D_pState -- Input-output operators.                       *
 **********************************************************************/
inline ostream &operator << (ostream &out_file, const Particle2D_pState &W) {
  out_file.setf(ios::scientific);
  out_file << " " << W.sigma << " " << W.u.x << " " << W.u.y << " " << W.Tp;
  out_file.unsetf(ios::scientific);
  return out_file;
}

inline istream &operator >> (istream &in_file, Particle2D_pState &W) {
  in_file.setf(ios::skipws);
  in_file >> W.sigma >> W.u.x >> W.u.y >> W.Tp;
  in_file.unsetf(ios::skipws);
  return in_file;
}

/**********************************************************************
 * Particle2D_cState::nd -- Particle-phase number density.            *
 **********************************************************************/
inline double Particle2D_cState::nd(const double &mp) const {
  return sigma/mp;
}

/**********************************************************************
 * Particle2D_cState::u -- Particle-phase flow velocity.              *
 **********************************************************************/
inline Vector2D Particle2D_cState::u(void) const {
  if (sigma < NANO) return Vector2D_ZERO;
  return du/sigma;
}

/**********************************************************************
 * Particle2D_cState::Tp -- Particle-phase temperature.               *
 **********************************************************************/
inline double Particle2D_cState::Tp(const double &cm) const {
  if (sigma < NANO || ep < NANO) return ZERO;
  return ep/(cm*sigma);
}
  
/**********************************************************************
 * Particle2D_cState::Ep -- Particle-phase total energy.              *
 **********************************************************************/
inline double Particle2D_cState::Ep(void) const {
  return ep + HALF*sigma*u().sqr();
}

/**********************************************************************
 * Particle2D_cState::jc -- Particle-phase current.                   *
 **********************************************************************/
inline Vector2D Particle2D_cState::jc(const double &qe, const double &mp) const {
  return qe*nd(mp)*u();
}

/**********************************************************************
 * Particle2D_cState -- Binary arithmetic operators.                  *
 **********************************************************************/
inline Particle2D_cState Particle2D_cState::operator +(const Particle2D_cState &U) const {
  return Particle2D_cState(sigma+U.sigma,du.x+U.du.x,du.y+U.du.y,ep+U.ep);
}

inline Particle2D_cState Particle2D_cState::operator -(const Particle2D_cState &U) const {
  return Particle2D_cState(sigma-U.sigma,du.x-U.du.x,du.y-U.du.y,ep-U.ep);
}

// Inner product operator.
inline double Particle2D_cState::operator *(const Particle2D_cState &U) const {
  return (sigma*U.sigma + du.x*U.du.x + du.y*U.du.y + ep*U.ep);
}

inline Particle2D_cState Particle2D_cState::operator *(const double &a) const {
  return Particle2D_cState(sigma*a,du.x*a,du.y*a,ep*a);
}

inline Particle2D_cState operator *(const double &a, const Particle2D_cState &U) {
  return Particle2D_cState(U.sigma*a,U.du.x*a,U.du.y*a,U.ep*a);
}

inline Particle2D_cState Particle2D_cState::operator /(const double &a) const {
  return Particle2D_cState(sigma/a,du.x/a,du.y/a,ep/a);
}

// A useful solution state product operator.
inline Particle2D_cState Particle2D_cState::operator ^(const Particle2D_cState &U) const {
  return Particle2D_cState(sigma*U.sigma,du.x*U.du.x,du.y*U.du.y,ep*U.ep);
}

/**********************************************************************
 * Particle2D_cState -- Assignment operator.                          *
 **********************************************************************/
inline Particle2D_cState& Particle2D_cState::operator =(const Particle2D_cState &U) {
  //if (this != &U) {
  sigma = U.sigma; du = U.du; ep = U.ep;
  //}
  return *this;
}

/**********************************************************************
 * Particle2D_cState -- Unary arithmetic operators.                   *
 **********************************************************************/
//inline Particle2D_cState operator +(const Particle2D_cState &U) {
//return Particle2D_cState(U.sigma,U.du.x,U.du.y,U.ep);
//}

inline Particle2D_cState operator -(const Particle2D_cState &U) {
  return Particle2D_cState(-U.sigma,-U.du.x,-U.du.y,-U.ep);
}

/**********************************************************************
 * Particle2D_cState -- Shortcut arithmetic operators.                *
 **********************************************************************/
inline Particle2D_cState& Particle2D_cState::operator +=(const Particle2D_cState &U) {
  sigma += U.sigma; du.x += U.du.x; du.y += U.du.y; ep += U.ep;
  return *this;
}

inline Particle2D_cState& Particle2D_cState::operator -=(const Particle2D_cState &U) {
  sigma -= U.sigma; du.x -= U.du.x; du.y -= U.du.y; ep -= U.ep;
  return *this;
}

inline Particle2D_cState& Particle2D_cState::operator *=(const double &a) {
  sigma *= a; du.x *= a; du.y *= a; ep *= a;
  return *this;
}

inline Particle2D_cState& Particle2D_cState::operator /=(const double &a) {
  sigma /= a; du.x /= a; du.y /= a; ep /= a;
  return *this;
}

/**********************************************************************
 * Particle2D_cState -- Relational operators.                         *
 **********************************************************************/
inline int operator ==(const Particle2D_cState &U1, const Particle2D_cState &U2) {
  return (U1.sigma == U2.sigma && U1.du == U2.du && U1.ep == U2.ep);
}

inline int operator !=(const Particle2D_cState &U1, const Particle2D_cState &U2) {
  return (U1.sigma != U2.sigma || U1.du != U2.du || U1.ep != U2.ep);
}

/**********************************************************************
 * Particle2D_cState -- Input-output operators.                       *
 **********************************************************************/
inline ostream &operator << (ostream &out_file, const Particle2D_cState &U) {
  out_file.setf(ios::scientific);
  out_file << " " << U.sigma << " " << U.du.x << " " << U.du.y << " " << U.ep;
  out_file.unsetf(ios::scientific);
  return out_file;
}

inline istream &operator >> (istream &in_file, Particle2D_cState &U) {
  in_file.setf(ios::skipws);
  in_file >> U.sigma >> U.du.x >> U.du.y >> U.ep;
  in_file.unsetf(ios::skipws);
  return in_file;
}

/**********************************************************************
 * Particle2D_pState::Particle2D_pState -- Constructor.               *
 **********************************************************************/
inline Particle2D_pState::Particle2D_pState(const Particle2D_cState &U, const double &cm) {
  sigma = U.sigma; u = U.u(); Tp = U.Tp(cm);
}

/**********************************************************************
 * Particle2D_pState::U -- Conserved solution state.                  *
 **********************************************************************/
inline Particle2D_cState Particle2D_pState::U(const double &cm) const {
  return Particle2D_cState(sigma,du(),ep(cm));
}

inline Particle2D_cState Particle2D_pState::U(const Particle2D_pState &W, const double &cm) const {
  return Particle2D_cState(W.sigma,W.du(),W.ep(cm));
}

/**********************************************************************
 * Particle_pState::dUdW -- Jacobian of the conserved solution        *
 *                          variables with respect to the primitive   *
 *                          solution variables.                       *
 **********************************************************************/
inline void Particle2D_pState::dUdW(DenseMatrix &dUdW, const int &n, const double &cm) const {
  if (sigma > ZERO) {
    dUdW(n  ,n  ) += ONE;
    dUdW(n+1,n  ) += u.x;
    dUdW(n+1,n+1) += sigma;
    dUdW(n+2,n  ) += u.y;
    dUdW(n+2,n+2) += sigma;
    dUdW(n+3,n  ) += ep(cm)/sigma;
    dUdW(n+3,n+3) += sigma*cm;
  }
}

/**********************************************************************
 * Particle2D_pState::dWdU -- Jacobian of the primitive solution      *
 *                            variables with respect to the conserved *
 *                            solution variables.                     *
 **********************************************************************/
inline void Particle2D_pState::dWdU(DenseMatrix &dWdU, const int &n, const double &cm) const {
  if (sigma > ZERO) {
    dWdU(n  ,n  ) += ONE;
    dWdU(n+1,n  ) -= u.x/sigma;
    dWdU(n+1,n+1) += ONE/sigma;
    dWdU(n+2,n  ) -= u.y/sigma;
    dWdU(n+2,n+2) += ONE/sigma;
    dWdU(n+3,n  ) -= Tp/sigma;
    dWdU(n+3,n+3) += ONE/(sigma*cm);
  }
}

/**********************************************************************
 * Particle2D_pState::F -- Solution flux and x-direction Jacobian.    *
 **********************************************************************/
inline Particle2D_cState Particle2D_pState::F(const double &cm) const {
  return Particle2D_cState(sigma*u.x,sigma*sqr(u.x),sigma*u.x*u.y,u.x*ep(cm));
}

inline Particle2D_cState Particle2D_pState::F(const Vector2D &V, const double &cm) const {
  return Particle2D_cState(sigma*(u.x-V.x),sigma*u.x*(u.x-V.x),
 			   sigma*u.y*(u.x-V.x),ep(cm)*(u.x-V.x));
}

inline void Particle2D_pState::dFdU(DenseMatrix &dFdU, const int &n, const double &cm) const {
  if (sigma > ZERO) {
    dFdU(n  ,n+1) += ONE;
    dFdU(n+1,n  ) -= u.x*u.x;
    dFdU(n+1,n+1) += TWO*u.x;
    dFdU(n+2,n  ) -= u.x*u.y;
    dFdU(n+2,n+1) += u.y;
    dFdU(n+2,n+2) += u.x;
    dFdU(n+3,n  ) -= u.x*ep(cm)/sigma;
    dFdU(n+3,n+1) += ep(cm)/sigma;
    dFdU(n+3,n+3) += u.x;
  }
}

inline void Particle2D_pState::dFdU(DenseMatrix &dFdU, const int &n,
				    const Vector2D &V, const double &cm) const {
  if (sigma > ZERO) {
    dFdU(n  ,n  ) -= V.x;
    dFdU(n  ,n+1) += ONE;
    dFdU(n+1,n  ) -= sqr(u.x);
    dFdU(n+1,n+1) += TWO*u.x-V.x;
    dFdU(n+2,n  ) -= u.y*u.x;
    dFdU(n+2,n+1) += u.y;
    dFdU(n+2,n+2) += u.x-V.x;
    dFdU(n+3,n  ) -= u.x*ep(cm)/sigma;
    dFdU(n+3,n+1) += ep(cm)/sigma;
    dFdU(n+3,n+3) += u.x-V.x;
  }
}

/**********************************************************************
 * Particle2D_pState::lambda_x -- Eigenvalue(s) (x-direction).        *
 **********************************************************************/
inline Particle2D_pState Particle2D_pState::lambda_x(void) const {
  return Particle2D_pState(u.x,u.x,u.x,u.x);
}

inline Particle2D_pState Particle2D_pState::lambda_x(const Vector2D &V) const {
  return Particle2D_pState(u.x-V.x,u.x-V.x,u.x-V.x,u.x-V.x);
}

/**********************************************************************
 * Particle2D_pState::rp_x -- Primitive right eigenvector (x-direction).*
 **********************************************************************/
inline Particle2D_pState Particle2D_pState::rp_x(int index) const {
  //assert(index >= 1 && index < NUM_VAR_PARTICLE2D);
  switch(index) {
  case 1 :
    return Particle2D_pState(ONE,ZERO,ZERO,ZERO);
  case 2 :
    return Particle2D_pState(ZERO,ZERO,ONE,ZERO);
  case 3 :
    return Particle2D_pState(ZERO,ZERO,ZERO,ONE);
  default:
    return Particle2D_pState(ONE,ZERO,ZERO,ZERO);
  };
}

/**********************************************************************
 * Particle2D_pState::rc_x -- Conserved right eigenvector (x-direction).*
 **********************************************************************/
inline Particle2D_cState Particle2D_pState::rc_x(int index) const {
  //assert(index >= 1 && index < NUM_VAR_PARTICLE2D);
  switch(index) {
  case 1 :
    return Particle2D_cState(ONE,u.x,ZERO,ZERO);
  case 2 :
    return Particle2D_cState(ZERO,ZERO,ONE,ZERO);
  case 3 :
    return Particle2D_cState(ZERO,ZERO,ZERO,ONE);
  default:
    return Particle2D_cState(ONE,u.x,ZERO,ZERO);
  };
}

/**********************************************************************
 * Particle2D_pState::lp_x -- Primitive left eigenvector (x-direction).*
 **********************************************************************/
// inline Particle2D_pState Particle2D_pState::lp_x(int index) const {
//   return Particle2D_pState(ZERO,ZERO,ZERO,ZERO);
// }

/**********************************************************************
 * Particle2D_pState::Sa -- Axisymmetric source term vector and       *
 *                          Jacobian.                                 *
 **********************************************************************/
inline Particle2D_cState Particle2D_pState::Sa(const Vector2D &X, const double &cm) const {
 if (sigma < NANO) return Particle2D_cState(ZERO,ZERO,ZERO,ZERO);
 return Particle2D_cState(-sigma*u.y/X.y,-sigma*u.x*u.y/X.y,-sigma*sqr(u.y)/X.y,-u.y*ep(cm)/X.y);
}

inline void Particle2D_pState::dSadU(DenseMatrix &dSadU, const Vector2D &X, const int &n, const double &cm) const {
  if (sigma > ZERO) {
    dSadU(n  ,n+1) -= ONE/X.y;
    dSadU(n+1,n  ) += u.x*u.y/X.y;
    dSadU(n+1,n+1) -= u.y/X.y;
    dSadU(n+1,n+2) -= u.x/X.y;
    dSadU(n+2,n  ) += u.y*u.y/X.y;
    dSadU(n+2,n+2) -= TWO*u.y/X.y;
    dSadU(n+3,n  ) += u.y*ep(cm)/(sigma*X.y);
    dSadU(n+3,n+2) -= ep(cm)/(sigma*X.y);
    dSadU(n+3,n+3) -= u.y/X.y;
  }
}

/**********************************************************************
 * Particle2D_pState::Se-- Electrostatic-force source term vector and *
 *                         Jacobian.                                  *
 **********************************************************************/
inline Particle2D_cState Particle2D_pState::Se(const Vector2D &E, const double &Cm) const {
 if (sigma < NANO) return Particle2D_cState(ZERO,ZERO,ZERO,ZERO);
 return Particle2D_cState(ZERO,sigma*Cm*E.x,sigma*Cm*E.y,sigma*Cm*(u*E));
}

inline void Particle2D_pState::dSedU(DenseMatrix &dSedU, const int &n,
				     const Vector2D &E, const double &Cm) const {
  if (sigma > ZERO) {
    dSedU(n+1,n  ) += Cm*E.x;
    dSedU(n+2,n  ) += Cm*E.y;
    dSedU(n+3,n+1) += Cm*E.x;
    dSedU(n+3,n+2) += Cm*E.y;
  }
}

/**********************************************************************
 * Particle2D_pState::Rotate -- Set the solution in the local rotated *
 *                              frame.                                *
 **********************************************************************/
inline void Particle2D_pState::Rotate(const Particle2D_pState &W,
				      const Vector2D &norm_dir) {
  sigma = W.sigma;
  u.x   =   W.u.x*norm_dir.x + W.u.y*norm_dir.y;
  u.y   = - W.u.x*norm_dir.y + W.u.y*norm_dir.x;
  Tp    = W.Tp;
}

/**********************************************************************
 * Particle2D_cState::Rotate -- Set the solution in the local rotated *
 *                              frame.                                *
 **********************************************************************/
inline void Particle2D_cState::Rotate(const Particle2D_cState &U,
				      const Vector2D &norm_dir) {
  sigma = U.sigma;
  du.x  =   U.du.x*norm_dir.x + U.du.y*norm_dir.y;
  du.y  = - U.du.x*norm_dir.y + U.du.y*norm_dir.x;
  ep    = U.ep;
}

/**********************************************************************
 * Particle2D_cState::Translate -- Set the solution in a stationary   *
 *                                 frame.                             *
 **********************************************************************/
inline void Particle2D_pState::Translate(const Particle2D_pState &W,
					 const Vector2D &V) {
  sigma = W.sigma;
  u.x   = W.u.x - V.x;
  u.y   = W.u.y - V.y;
  Tp    = W.Tp;
}

/**********************************************************************
 * Particle2D_cState::Reflect -- Set the reflected solution state in  *
 *                               a given direction given the          *
 *                               primitive solution variables and the *
 *                               unit normal vector in the direction  *
 *                               of interest.                         *
 **********************************************************************/
inline void Particle2D_pState::Reflect(const Particle2D_pState &W,
				       const Vector2D &norm_dir) {
  
  Particle2D_pState Wr;

  // Apply the frame rotation and calculate the primitive solution 
  // state variables in the local rotated frame defined by the unit 
  // normal vector.
  Wr.Rotate(W,norm_dir);

  // Reflect the gas-phase normal velocity in the rotated frame.
  Wr.u.x = -Wr.u.x;

  // Rotate back to the original Cartesian reference frame.
  Rotate(Wr,Vector2D(norm_dir.x,-norm_dir.y));

}

/**********************************************************************
 * Particle2D_pState::Mirror -- Set the mirrored solution state in a  *
 *                              given direction given the primitive   *
 *                              solution variables and the unit       *
 *                              normal vector in the direction of     *
 *                              interest.                             *
 **********************************************************************/
inline void Particle2D_pState::Mirror(const Particle2D_pState &W,
				      const Vector2D &norm_dir) {
  
  Particle2D_pState Wr;

  // Apply the frame rotation and calculate the primitive solution 
  // state variables in the local rotated frame defined by the unit 
  // normal vector.
  Wr.Rotate(W,norm_dir);

  // Mirror the gas-phase normal velocity in the rotated frame.
  Wr.u.x = -Wr.u.x;
  Wr.u.y = -Wr.u.y;

  // Rotate back to the original Cartesian reference frame.
  Rotate(Wr,Vector2D(norm_dir.x,-norm_dir.y));

}

/**********************************************************************
 * Particle2D_pState::Absorb -- Set the absorbed (zero) solution      *
 *                              state in a given direction given the  *
 *                              primitive solution variables and the  *
 *                              unit normal vector in the direction   *
 *                              of interest.                          *
 **********************************************************************/
inline void Particle2D_pState::Absorb(const Particle2D_pState &W,
				      const Vector2D &norm_dir) {
  
  Particle2D_pState Wr;
  
  // Apply the frame rotation and calculate the primitive solution 
  // state variables in the local rotated frame defined by the unit 
  // normal vector.
  Wr.Rotate(W,norm_dir);

  // Absorb the particle-phase.
  if (Wr.sigma > ZERO && Wr.u.x < ZERO) Wr.Vacuum();

  // Rotate back to the original Cartesian reference frame.
  Rotate(Wr,Vector2D(norm_dir.x,-norm_dir.y));

}

/**********************************************************************
 * Particle2D_cState::Particle2D_cState -- Constructor.               *
 **********************************************************************/
inline Particle2D_cState::Particle2D_cState(const Particle2D_pState &W, const double &cm) {
  sigma = W.sigma; du = W.du(); ep = W.ep(cm);
}

/**********************************************************************
 * Particle2D_cState::W -- Primitive solution state.                  *
 **********************************************************************/
inline Particle2D_pState Particle2D_cState::W(const double &cm) const {
  return Particle2D_pState(sigma,u(),Tp(cm));
}

inline Particle2D_pState Particle2D_cState::W(const Particle2D_cState &U, const double &cm) const {
  return Particle2D_pState(U.sigma,U.u(),U.Tp(cm));
}

/**********************************************************************
 * Particle_cState::dUdW -- Jacobian of the conserved solution        *
 *                          variables with respect to the primitive   *
 *                          solution variables.                       *
 **********************************************************************/
inline void Particle2D_cState::dUdW(DenseMatrix &dUdW, const int &n, const double &cm) const {
  if (sigma > ZERO) {
    dUdW(n  ,n  ) += ONE;
    dUdW(n+1,n  ) += u().x;
    dUdW(n+1,n+1) += sigma;
    dUdW(n+2,n  ) += u().y;
    dUdW(n+2,n+2) += sigma;
    dUdW(n+3,n  ) += ep/sigma;
    dUdW(n+3,n+3) += sigma*cm;
  }
}

/**********************************************************************
 * Particle2D_pState::dWdU -- Jacobian of the primitive solution      *
 *                            variables with respect to the conserved *
 *                            solution variables.                     *
 **********************************************************************/
inline void Particle2D_cState::dWdU(DenseMatrix &dWdU, const int &n, const double &cm) const {
  if (sigma > ZERO) {
    dWdU(n  ,n  ) += ONE;
    dWdU(n+1,n  ) -= u().x/sigma;
    dWdU(n+1,n+1) += ONE/sigma;
    dWdU(n+2,n  ) -= u().y/sigma;
    dWdU(n+2,n+2) += ONE/sigma;
    dWdU(n+3,n  ) -= Tp(cm)/sigma;
    dWdU(n+3,n+3) += ONE/(sigma*cm);
  }
}

/**********************************************************************
 * Particle2D_cState::F -- Solution flux and x-direction Jacobian.    *
 **********************************************************************/
inline Particle2D_cState Particle2D_cState::F(void) const {
  if (sigma < NANO) return Particle2D_cState(ZERO,ZERO,ZERO,ZERO);
  return Particle2D_cState(du.x,sqr(du.x)/sigma,du.x*du.y/sigma,du.x*ep/sigma);
}

inline Particle2D_cState Particle2D_cState::F(const Vector2D &V) const {
  if (sigma < NANO) return Particle2D_cState(ZERO,ZERO,ZERO,ZERO);
  return Particle2D_cState(sigma*(u().x-V.x),du.x*(u().x-V.x),
 			   du.y*(u().x-V.x),ep*(u().x-V.x));
}

inline void Particle2D_cState::F(const Particle2D_pState &W, const double &cm) {
  if (W.sigma < NANO) Vacuum();
  else Copy(W.sigma*W.u.x,W.sigma*sqr(W.u.x),W.sigma*W.u.x*W.u.y,W.u.x*W.ep(cm));
}

inline void Particle2D_cState::F(const Particle2D_pState &W, const Vector2D &V, const double &cm) {
  if (W.sigma < NANO) Vacuum();
  else Copy(W.sigma*(W.u.x-V.x),W.sigma*W.u.x*(W.u.x-V.x),W.sigma*W.u.y*(W.u.x-V.x),W.ep(cm)*(W.u.x-V.x));
}

inline void Particle2D_cState::dFdU(DenseMatrix &dFdU, const int &n) const {
  if (sigma > ZERO) {
    double ux = u().x, uy = u().y;
    dFdU(n  ,n+1) += ONE;
    dFdU(n+1,n  ) -= sqr(ux);
    dFdU(n+1,n+1) += TWO*ux;
    dFdU(n+2,n  ) -= ux*uy;
    dFdU(n+2,n+1) += uy;
    dFdU(n+2,n+2) += ux;
    dFdU(n+3,n  ) -= ux*ep/sigma;
    dFdU(n+3,n+1) += ep/sigma;
    dFdU(n+3,n+3) += ux;
  }
}

inline void Particle2D_cState::dFdU(DenseMatrix &dFdU, const int &n,
				    const Vector2D &V) const {
  if (sigma > ZERO) {
    double ux = u().x, uy = u().y;
    dFdU(n  ,n  ) -= V.x;
    dFdU(n  ,n+1) += ONE;
    dFdU(n+1,n  ) -= sqr(ux);
    dFdU(n+1,n+1) += TWO*ux-V.x;
    dFdU(n+2,n  ) -= uy*ux;
    dFdU(n+2,n+1) += uy;
    dFdU(n+2,n+2) += ux-V.x;
    dFdU(n+3,n  ) -= ux*ep/sigma;
    dFdU(n+3,n+1) += ep/sigma;
    dFdU(n+3,n+3) += ux-V.x;
  }
}

/**********************************************************************
 * Particle2D_cState::lambda_x -- Eigenvalue(s) (x-direction).        *
 **********************************************************************/
inline Particle2D_pState Particle2D_cState::lambda_x(void) const {
  return Particle2D_pState(u().x,u().x,u().x,u().x);
}

inline Particle2D_pState Particle2D_cState::lambda_x(const Vector2D &V) const {
  return Particle2D_pState(u().x-V.x,u().x-V.x,u().x-V.x,u().x-V.x);
}

/**********************************************************************
 * Particle2D_cState::rp_x -- Primitive right eigenvector (x-direction).*
 **********************************************************************/
inline Particle2D_pState Particle2D_cState::rp_x(int index) const {
  //assert(index >= 1 && index < NUM_VAR_PARTICLE2D);
  switch(index) {
  case 1 :
    return Particle2D_pState(ONE,ZERO,ZERO,ZERO);
  case 2 :
    return Particle2D_pState(ZERO,ZERO,ONE,ZERO);
  case 3 :
    return Particle2D_pState(ZERO,ZERO,ZERO,ONE);
  default:
    return Particle2D_pState(ONE,ZERO,ZERO,ZERO);
  };
}

/**********************************************************************
 * Particle2D_cState::rc_x -- Conserved right eigenvector (x-direction).*
 **********************************************************************/
inline Particle2D_cState Particle2D_cState::rc_x(int index) const {
  //assert(index >= 1 && index < NUM_VAR_PARTICLE2D);
  switch(index) {
  case 1 :
    return Particle2D_cState(ONE,u().x,ZERO,ZERO);
  case 2 :
    return Particle2D_cState(ZERO,ZERO,ONE,ZERO);
  case 3 :
    return Particle2D_cState(ZERO,ZERO,ZERO,ONE);
  default:
    return Particle2D_cState(ONE,u().x,ZERO,ZERO);
  };
}

/**********************************************************************
 * Particle2D_cState::lp_x -- Primitive left eigenvector (x-direction).*
 **********************************************************************/
// inline Particle2D_pState Particle2D_cState::lp_x(int index) const {
//   return Particle2D_pState(ZERO,ZERO,ZERO,ZERO);
// }

/**********************************************************************
 * Particle2D_pState::FluxSaurel -- Particle-phase Riemann solver     *
 *                                  proposed by Saurel, Daniel, and   *
 *                                  Loraud (AIAA J. 32:6 1992).  This *
 *                                  function returns the x-direction  *
 *                                  intermediate state solution flux  *
 *                                  for the given left and right      *
 *                                  solution states.                  *
 **********************************************************************/
inline void Particle2D_pState::FluxSaurel(const Particle2D_pState &Wl,
					  const Particle2D_pState &Wr) {

  // Determine the particle-phase intermediate solution state.
  if (Wl.sigma < NANO && Wr.sigma < NANO) {
    Vacuum();
  } else {
    if (Wl.u.x < Wr.u.x) {
      if (Wl.u.x >= ZERO && Wr.u.x > ZERO) {
	Copy(Wl);
      } else if (Wl.u.x < ZERO && Wr.u.x > ZERO) {
	Vacuum();
      } else if (Wl.u.x < ZERO && Wr.u.x <= ZERO) {
	Copy(Wr);
      }
    } else if (Wl.u.x >= Wr.u.x) {
      if (Wl.u.x > ZERO && Wr.u.x >= ZERO) {
	Copy(Wl);
      } else if (Wl.u.x >= ZERO && Wr.u.x <= ZERO) {
	sigma = Wl.sigma + Wr.sigma;
	if (sigma < NANO) {
	  Vacuum();
	} else {
	  u.x = (Wl.sigma*Wl.u.x + Wr.sigma*Wr.u.x)/sigma;
	  u.y = (Wl.sigma*Wl.u.y + Wr.sigma*Wr.u.y)/sigma;
	  Tp  = (Wl.sigma*Wl.Tp  + Wr.sigma*Wr.Tp )/sigma;
	}
      } else if (Wl.u.x <= ZERO && Wr.u.x < ZERO) {
	Copy(Wr);
      }
    }
    // Error checking.
    if (sigma < NANO || Tp < NANO) Vacuum();
  }

}

/**********************************************************************
 * Particle2D_pState::FluxSaurel_n -- Particle-phase Riemann solver   *
 *                                    proposed by Saurel, Daniel, and *
 *                                    Loraud (AIAA J. 32:6 1992).     *
 *                                    This function returns the n-    *
 *                                    direction intermediate state    *
 *                                    solution flux for the given     *
 *                                    left and right solution states. *
 **********************************************************************/
inline void Particle2D_cState::FluxSaurel_n(const Particle2D_pState &Wl,
					    const Particle2D_pState &Wr,
					    const Vector2D &norm_dir,
					    const double &cm) {

  Particle2D_pState Wl_rotated, Wr_rotated, W_rotated;
  Particle2D_cState Flux_rotated;

  // Apply the frame rotation and evaluate left and right solution 
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated.Rotate(Wl,norm_dir);
  Wr_rotated.Rotate(Wr,norm_dir);

  // Evaluate the intermediate solution state in the rotated frame.
  W_rotated.FluxSaurel(Wl_rotated,Wr_rotated);

  // Evaluate the intermediate state solution flux in the rotated 
  // frame.
  Flux_rotated.F(W_rotated,cm);

  // Rotate back to the original Cartesian reference frame and return 
  // the solution flux.
  Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

}

/**********************************************************************
 * Particle2D_pState::FluxSaurel_MB_n -- Particle-phase Riemann       *
 *                                       solution proposed by Saurel, *
 *                                       Daniel, and Loraud (AIAA J.  *
 *                                       32:6 1992).  This function   *
 *                                       returns the n-direction      *
 *                                       intermediate state solution  *
 *                                       flux for the given left and  *
 *                                       right solution states at a   *
 *                                       nonstationary interface.     *
 **********************************************************************/
inline void Particle2D_cState::FluxSaurel_MB_n(const Particle2D_pState &Wl,
					       const Particle2D_pState &Wr,
					       const Vector2D &V,
					       const Vector2D &norm_dir,
					       const double &cm) {

  Particle2D_pState Wl_rotated, Wr_rotated, W_rotated;
  Particle2D_cState Flux_rotated;
  Vector2D V_rotated;

  // Apply the frame rotation and evaluate left and right solution 
  // states in the local rotated frame defined by the unit normal 
  // vector.
  Wl_rotated.Rotate(Wl,norm_dir);
  Wr_rotated.Rotate(Wr,norm_dir);
  V_rotated.Rotate(V,norm_dir);

  // Determine the intermediate solution state in the rotated frame.
  W_rotated.FluxSaurel(Wl_rotated,Wr_rotated);

  // Transform back into the moving body frame of reference.
  W_rotated.u.x += V_rotated.x;

  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated.F(W_rotated,V_rotated,cm);

  // Rotate back to the original Cartesian reference frame and return 
  // the solution flux.
  Rotate(Flux_rotated,Vector2D(norm_dir.x,-norm_dir.y));

}

#endif // _PARTICLE2D_STATE_INCLUDED
