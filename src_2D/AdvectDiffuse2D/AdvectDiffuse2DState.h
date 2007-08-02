/* AdvectDiffuse2DState.h:  Header file defining 
                            2D Advection Diffusion Equation Solution State Class. */

#ifndef _ADVECTDIFFUSE2D_STATE_INCLUDED
#define _ADVECTDIFFUSE2D_STATE_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

/* Include math macro, CFD, and 2D vector header files. */

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _VECTOR2D_INCLUDED
#include "../Math/Vector2D.h"
#endif //_VECTOR2D_INCLUDED

/* Define the classes. */

#define	NUM_VAR_ADVECTDIFFUSE2D    1

/*!
 * Class: AdvectDiffuse2D_State
 *
 * @brief Solution state class definition for the 2D advection-diffusion
 *        equation.
 *
 * Solution state class definition for the 2D advection-diffusion
 * equation.
 *
 * \verbatim
 * Member functions
 *     u        -- Return solution.
 *     V        -- Return advection velocity.
 *     k        -- Return diffusion coefficient.
 *     T        -- Return relaxation time.
 *    Fd        -- Return diffusive flux.
 *     F        -- Return advective flux.
 *     s        -- Return regular source term.
 *     s_axi    -- Return axisymmetric source term.
 *    F_diff    -- Evaluates diffusive flux.
 *
 * Member operators
 *      U -- a solution state
 *      c -- a scalar (double)
 *
 * U = U;
 * U = U + U;
 * U = U - U;
 * U = U * U;
 * U = c * U;
 * U = U * c;
 * U = U / c;
 * U = +U;
 * U = -U;
 * U += U;
 * U -= U;
 * U *= a;
 * U /= a;
 * U == U;
 * U != U;
 * cout << U; (output function)
 * cin  >> U; (input function)
 * \endverbatim
 */
class AdvectDiffuse2D_State{
private:
public:
  //@{ @name Solution state variables and associated constants:
  double          u;   // Solution.
  Vector2D        V;   // Advection velocity (2D vector).
  double          k;   // Diffusion coefficient.
  double          T;   // Relaxation time for source terms.
  Vector2D       Fd;   // Diffusive flux (2D vector).
  //@}

  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  AdvectDiffuse2D_State(void) {
    u = ZERO; V.zero(); k = ONE; T = ONE; Fd.zero();
  }

  //! Copy constructor.
  AdvectDiffuse2D_State(const AdvectDiffuse2D_State &U) {
    u = U.u; V = U.V; k = U.k; T = U.T; Fd = U.Fd;
  }

  //! Assignment constructor.
  AdvectDiffuse2D_State(const double &uu,
			const Vector2D &VV,
			const double &kappa,
			const double &tau,
			const Vector2D &FF) {
    u = uu; V = VV; k = kappa; T = tau; Fd = FF;
  }

  //! Assignment constructor.
  AdvectDiffuse2D_State(const double &uu,
			const double &Vx,
			const double &Vy,
			const double &kappa,
			const double &tau) {
    u = uu; V.x = Vx; V.y = Vy; k = kappa; T = tau; Fd.zero();
  }
    
  //! Assignment constructor.
  AdvectDiffuse2D_State(const double &uu,
			const double &Vx,
			const double &Vy,
			const double &kappa,
			const double &tau,
			const double &Fdx,
			const double &Fdy) {
    u = uu; V.x = Vx; V.y = Vy; k = kappa; T = tau; Fd.x = Fdx; Fd.y = Fdy;
  }

  /* Destructor. */
  // ~AdvectDiffuse2D_State(void);
  // Use automatically generated destructor.
  //@}

  //@{ @name Useful operators.
  //! Vacuum/zero operator.
  void Vacuum(void) { u = ZERO; }
  //@}

  //@{ @name Advective Flux.
  Vector2D F(void);
  Vector2D F(void) const;
  friend Vector2D F(const AdvectDiffuse2D_State &U);
  //@}

  //@{ @name Regular source term.
  double s(void);
  double s(void) const;
  friend double s(const AdvectDiffuse2D_State &U);
  //@}

  //@{ @name Axisymmetric source term.
  double s_axi(const Vector2D &X);
  double s_axi(const Vector2D &X) const;
  friend double s_axi(const AdvectDiffuse2D_State &U, const Vector2D &X);
  //@}

  //@{ @name Evaluates diffusive flux.
  Vector2D F_diff(const double &dudx, const double &dudy);
  Vector2D F_diff(const double &dudx, const double &dudy) const;
  friend Vector2D F_diff(const AdvectDiffuse2D_State &U, const double &dudx, const double &dudy);
  //@}

  /* @name Assignment operator. */
  // AdvectDiffuse2D_State operator = (const AdvectDiffuse2D_State &W);
  // Use automatically generated assignment operator.

  //@{ @name Binary arithmetic operators.
  friend AdvectDiffuse2D_State operator +(const AdvectDiffuse2D_State &U1, const AdvectDiffuse2D_State &U2);
  friend AdvectDiffuse2D_State operator -(const AdvectDiffuse2D_State &U1, const AdvectDiffuse2D_State &U2);
  friend AdvectDiffuse2D_State operator *(const AdvectDiffuse2D_State &U1, const AdvectDiffuse2D_State &U2);
  friend AdvectDiffuse2D_State operator *(const AdvectDiffuse2D_State &U, const double &a);
  friend AdvectDiffuse2D_State operator *(const double &a, const AdvectDiffuse2D_State &U);
  friend AdvectDiffuse2D_State operator /(const AdvectDiffuse2D_State &U, const double &a);
  //@}

  //@{ @name Unary arithmetic operators.
  friend AdvectDiffuse2D_State operator +(const AdvectDiffuse2D_State &U);
  friend AdvectDiffuse2D_State operator -(const AdvectDiffuse2D_State &U);
  //@}

  //@{ @name Shortcut arithmetic operators.
  friend AdvectDiffuse2D_State &operator +=(AdvectDiffuse2D_State &U1, const AdvectDiffuse2D_State &U2);
  friend AdvectDiffuse2D_State &operator -=(AdvectDiffuse2D_State &U1, const AdvectDiffuse2D_State &U2);
  AdvectDiffuse2D_State &operator *=(const double &a);
  AdvectDiffuse2D_State &operator /=(const double &a);
  //@}

  //@{ @name Relational operators.
  friend int operator ==(const AdvectDiffuse2D_State &U1, const AdvectDiffuse2D_State &U2);
  friend int operator !=(const AdvectDiffuse2D_State &U1, const AdvectDiffuse2D_State &U2);
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const AdvectDiffuse2D_State &U);
  friend istream &operator >> (istream &in_file,  AdvectDiffuse2D_State &U);
  //@}

};

/*************************************************************
 * AdvectDiffuse2D_State::F -- Advective Flux.               *
 *************************************************************/
inline Vector2D AdvectDiffuse2D_State::F(void) {
    return (u*V);
}

inline Vector2D AdvectDiffuse2D_State::F(void) const {
    return (u*V);
}

/*************************************************************
 * AdvectDiffuse2D_State::s -- Regular source term.          *
 *************************************************************/
inline double AdvectDiffuse2D_State::s(void) {
  return (-u/T);
}

inline double AdvectDiffuse2D_State::s(void) const {
  return (-u/T);
}

inline double s(const AdvectDiffuse2D_State &U) {
  return (-U.u/U.T);
}

/*************************************************************
 * AdvectDiffuse2D_State::s_axi -- Axisymmetric source term. *
 *************************************************************/
inline double AdvectDiffuse2D_State::s_axi(const Vector2D &X) {
    return (ZERO);
}

inline double AdvectDiffuse2D_State::s_axi(const Vector2D &X) const {
    return (ZERO);
}

inline double s_axi(const AdvectDiffuse2D_State &U, const Vector2D &X) {
    return (ZERO);
}

/*************************************************************
 * AdvectDiffuse2D_State::F_diff - Evaluate diffusive flux.  *
 *************************************************************/
inline Vector2D AdvectDiffuse2D_State::F_diff(const double &dudx, const double &dudy) {
    return (Vector2D(-k*dudx, -k*dudy));   
}

inline Vector2D AdvectDiffuse2D_State::F_diff(const double &dudx, const double &dudy) const {
    return (Vector2D(-k*dudx, -k*dudy));   
}

inline Vector2D F_diff(const AdvectDiffuse2D_State &U, const double &dudx, const double &dudy) {
    return (Vector2D(-U.k*dudx, -U.k*dudy));
}

/*************************************************************
 * AdvectDiffuse2D_State -- Binary arithmetic operators.     *
 *************************************************************/
inline AdvectDiffuse2D_State operator +(const AdvectDiffuse2D_State &U1, const AdvectDiffuse2D_State &U2) {
  return (AdvectDiffuse2D_State(U1.u+U2.u,U1.V,U1.k,U1.T,U1.Fd));
}

inline AdvectDiffuse2D_State operator -(const AdvectDiffuse2D_State &U1, const AdvectDiffuse2D_State &U2) {
  return (AdvectDiffuse2D_State(U1.u-U2.u,U1.V,U1.k,U1.T,U1.Fd));
}

inline AdvectDiffuse2D_State operator *(const AdvectDiffuse2D_State &U1, const AdvectDiffuse2D_State &U2) {
  return (AdvectDiffuse2D_State(U1.u*U2.u,U1.V,U1.k,U1.T,U1.Fd));
}

inline AdvectDiffuse2D_State operator *(const AdvectDiffuse2D_State &U, const double &a) {
  return (AdvectDiffuse2D_State(a*U.u,U.V,U.k,U.T,U.Fd));
}

inline AdvectDiffuse2D_State operator *(const double &a, const AdvectDiffuse2D_State &U) {
  return (AdvectDiffuse2D_State(a*U.u,U.V,U.k,U.T,U.Fd));
}

inline AdvectDiffuse2D_State operator /(const AdvectDiffuse2D_State &U, const double &a) {
  return (AdvectDiffuse2D_State(U.u/a,U.V,U.k,U.T,U.Fd));
}

/*************************************************************
 * AdvectDiffuse2D_State -- Unary arithmetic operators.      *
 *************************************************************/
inline AdvectDiffuse2D_State operator +(const AdvectDiffuse2D_State &U) {
  return (AdvectDiffuse2D_State(U.u,U.V,U.k,U.T,U.Fd));
}

inline AdvectDiffuse2D_State operator -(const AdvectDiffuse2D_State &U) {
  return (AdvectDiffuse2D_State(-U.u,U.V,U.k,U.T,U.Fd));
}

/*************************************************************
 * AdvectDiffuse2D_State -- Shortcut arithmetic operators.   *
 *************************************************************/
inline AdvectDiffuse2D_State &operator +=(AdvectDiffuse2D_State &U1, const AdvectDiffuse2D_State &U2) {
  U1.u += U2.u;  
  return (U1);
}

inline AdvectDiffuse2D_State &operator -=(AdvectDiffuse2D_State &U1, const AdvectDiffuse2D_State &U2) {
  U1.u -= U2.u;
  return (U1);
}

inline AdvectDiffuse2D_State& AdvectDiffuse2D_State::operator *=(const double &a) {
  u *= a;
  return *this;
}

inline AdvectDiffuse2D_State& AdvectDiffuse2D_State::operator /=(const double &a) {
  u /= a;
  return *this;
}

/*************************************************************
 * AdvectDiffuse2D_State -- Relational operators.            *
 *************************************************************/
inline int operator ==(const AdvectDiffuse2D_State &U1, const AdvectDiffuse2D_State &U2) {
  return (U1.u == U2.u);
}

inline int operator !=(const AdvectDiffuse2D_State &U1, const AdvectDiffuse2D_State &U2) {
  return (U1.u != U2.u);
}

/*************************************************************
 * AdvectDiffuse2D_State -- Input-output operators.          *
 *************************************************************/
inline ostream &operator << (ostream &out_file, const AdvectDiffuse2D_State &U) {
  out_file.setf(ios::scientific);
  out_file << " " << U.u  << " " << U.Fd.x << " " << U.Fd.y 
           << " " << U.V.x << " " << U.V.y << " " << U.k << " " << U.T;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, AdvectDiffuse2D_State &U) {
  in_file.setf(ios::skipws);
  in_file >> U.u >> U.Fd.x >> U.Fd.y >> U.V.x >> U.V.y >> U.k >> U.T;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/*************************************************************
 * AdvectDiffuse2D_State -- External subroutines.            *
 *************************************************************/

extern AdvectDiffuse2D_State Riemann(const AdvectDiffuse2D_State &Ul,
	      	                     const AdvectDiffuse2D_State &Ur);

extern AdvectDiffuse2D_State Riemann_x(const AdvectDiffuse2D_State &Ul,
	      	                       const AdvectDiffuse2D_State &Ur);

extern AdvectDiffuse2D_State Riemann_y(const AdvectDiffuse2D_State &Ul,
	      	                       const AdvectDiffuse2D_State &Ur);

extern double Flux(const AdvectDiffuse2D_State &Ul,
                   const AdvectDiffuse2D_State &Ur);
       
extern double Flux_x(const AdvectDiffuse2D_State &Ul,
                     const AdvectDiffuse2D_State &Ur);
       
extern double Flux_y(const AdvectDiffuse2D_State &Ul,
                     const AdvectDiffuse2D_State &Ur);
       
extern double Flux_n(const AdvectDiffuse2D_State &Ul,
                     const AdvectDiffuse2D_State &Ur,
                     const Vector2D &norm_dir);

extern AdvectDiffuse2D_State BC_Dirichlet(const AdvectDiffuse2D_State &Ui,
                                          const AdvectDiffuse2D_State &Uo,
	      	                          const Vector2D &norm_dir);

extern AdvectDiffuse2D_State BC_Neumann(const AdvectDiffuse2D_State &Ui,
                                        const AdvectDiffuse2D_State &Uo,
	      	                        const Vector2D &norm_dir);

extern AdvectDiffuse2D_State BC_Robin(const AdvectDiffuse2D_State &Ui,
                                      const AdvectDiffuse2D_State &Uo,
	      	                      const Vector2D &norm_dir);

#endif /* _DVECTDIFFUSE2D_STATE_INCLUDED  */
