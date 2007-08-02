/**********************************************************************
 * Electrostatic2DState.h: Header file defining 2D electrostatic      *
 *                         solution state class.                      *
 **********************************************************************/

#ifndef _ELECTROSTATIC2D_STATE_INCLUDED
#define _ELECTROSTATIC2D_STATE_INCLUDED

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

// Define the number of electrostatic variables.
#define NUM_VAR_ELECTROSTATIC2D  3

/*!
 * Class: Electrostatic2DState
 *
 * @brief Electrostatic solution state class definition.
 *
 * Definition of a electrostatic solution class containing the electric
 * and potential field values.  Included by the Dusty2D_pState and 
 * Dusty2D_cState classes when considering problems involving the
 * tranport of charged solid due to an electric field.
 *
 * \verbatim
 * Member functions
 *      E -- Electric field vector.
 *      V -- Electric potential.
 *
 * Member operators
 *      U -- a solution state
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
 * U = U ^ U; (my useful product)
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
class Electrostatic2DState {
 private:
 public:
  Vector2D                E; //!< Electric field.
  double                  V; //!< Electric potential.

  //@{ @name Constructors and destructor.
  //! Creation constructor.
  Electrostatic2DState(void) { E = Vector2D_ZERO; V = ZERO; }

  //! Assignment constructor.
  Electrostatic2DState(const double &Ex, const double &Ey, const double &v) {
    E.x = Ex; E.y = Ey; V = v;
  }

  //! Assignment constructor.
  Electrostatic2DState(const Vector2D &e, const double &v) {
    E = e; V = v;
  }

  //! Copy constructor.
  Electrostatic2DState(const Electrostatic2DState &U) { E = U.E; V = U.V; }

  //! Destructor.
  ~Electrostatic2DState(void) { E = Vector2D_ZERO; V = ZERO; }
  //@}

  //! Return the number of variables.
  int NumVar(void) { return NUM_VAR_ELECTROSTATIC2D; }

  //! Copy operator.
  void Copy(const Electrostatic2DState &U) { E = U.E; V = U.V; }

  //! Vacuum operator.
  void Vacuum(void) { E = Vector2D_ZERO; V = ZERO; }

  //! Constant operator.
  void Constant(const double &val) { E = Vector2D(val,val); V = val; }

  //@{ @name Index operator.
  double &operator[](int index) {
    assert(index >= 1 && index <= NUM_VAR_ELECTROSTATIC2D);
    switch(index) {
    case 1 :
      return E.x;
    case 2 :
      return E.y;
    case 3 :
      return V;
    };
    return E.x;
  }
  const double &operator[](int index) const {
    assert(index >= 1 && index <= NUM_VAR_ELECTROSTATIC2D);
    switch(index) {
    case 1 :
      return E.x;
    case 2 :
      return E.y;
    case 3 :
      return V;
    };
    return E.x;
  }
  //@}

  //@{ @name Axisymmetric source vector and Jacobian.
  Electrostatic2DState S(const Vector2D &X,
			 const Electrostatic2DState &dUdy) const;
  void dSdU(DenseMatrix &dSdU, const Vector2D &X, const Electrostatic2DState &dUdy) const;
  //@}

  //@{ @name Binary arithmetic operators.
  Electrostatic2DState operator +(const Electrostatic2DState &U) const;
  Electrostatic2DState operator -(const Electrostatic2DState &U) const;
  double operator *(const Electrostatic2DState &U) const;
  Electrostatic2DState operator *(const double &a) const;
  friend Electrostatic2DState operator *(const double &a, const Electrostatic2DState &U);
  Electrostatic2DState operator /(const double &a) const;
  Electrostatic2DState operator ^(const Electrostatic2DState &U1) const;
  //@}

  //! Assignment Operator.
  Electrostatic2DState& operator =(const Electrostatic2DState &U);

  //@{ @name Unary arithmetic operators.
  //Electrostatic2DState operator +(const Electrostatic2DState &U);
  friend Electrostatic2DState operator -(const Electrostatic2DState &U);
  //@}

  //@{ @name Shortcut arithmetic operators.
  Electrostatic2DState &operator +=(const Electrostatic2DState &U);
  Electrostatic2DState &operator -=(const Electrostatic2DState &U);
  Electrostatic2DState &operator *=(const double &a);
  Electrostatic2DState &operator /=(const double &a);
  //@}

  //@{ @name Relational operators.
  friend int operator ==(const Electrostatic2DState &U1, const Electrostatic2DState &U2);
  friend int operator !=(const Electrostatic2DState &U1, const Electrostatic2DState &U2);
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const Electrostatic2DState &U);
  friend istream &operator >> (istream &in_file, Electrostatic2DState &U);
  //@}

  void output_labels(ostream &out_file) {
    out_file << "\"V\" \\ \n"
	     << "\"Ex\" \\ \n"
	     << "\"Ey\" \\ \n";
  }

  void output_data(ostream &out_file) {
    out_file << " " << V << " " << E.x << " " << E.y;
  }

};

/**********************************************************************
 * Electrostatic2DState::S -- Axisymmetric source term and Jacobian.  *
 **********************************************************************/
inline Electrostatic2DState Electrostatic2DState::S(const Vector2D &X,
						    const Electrostatic2DState &dUdy) const {
  return Electrostatic2DState(ZERO,ZERO,dUdy.V/X.y);
}

inline void Electrostatic2DState::dSdU(DenseMatrix &dSdU,
				       const Vector2D &X,
				       const Electrostatic2DState &dUdy) const {
  dSdU(0,0) -= ZERO/X.y;
  dSdU(1,0) -= ZERO/X.y;
  dSdU(0,1) -= ZERO/X.y;
  dSdU(1,1) -= ZERO/X.y;
}

/**********************************************************************
 * Electrostatic2DState -- Binary arithmetic operators.               *
 **********************************************************************/
inline Electrostatic2DState Electrostatic2DState::operator +(const Electrostatic2DState &U) const {
  return Electrostatic2DState(E.x+U.E.x,E.y+U.E.y,V+U.V);
}

inline Electrostatic2DState Electrostatic2DState::operator -(const Electrostatic2DState &U) const {
  return Electrostatic2DState(E.x-U.E.x,E.y-U.E.y,V-U.V);
}

// Inner product operator.
inline double Electrostatic2DState::operator *(const Electrostatic2DState &U) const {
  return E.x*U.E.x + E.y*U.E.y + V*U.V;
}

inline Electrostatic2DState Electrostatic2DState::operator *(const double &a) const {
  return Electrostatic2DState(E.x*a,E.y*a,V*a);
}

inline Electrostatic2DState operator *(const double &a, const Electrostatic2DState &U) {
  return Electrostatic2DState(U.E.x*a,U.E.y*a,U.V*a);
}

inline Electrostatic2DState Electrostatic2DState::operator /(const double &a) const {
  return Electrostatic2DState(E.x/a,E.y/a,V/a);
}

// A useful solution state product operator.
inline Electrostatic2DState Electrostatic2DState::operator ^(const Electrostatic2DState &U) const {
  return Electrostatic2DState(E.x*U.E.x,E.y*U.E.y,V*U.V);
}

/**********************************************************************
 * Electrostatic2DState -- Assignment operator.                       *
 **********************************************************************/
inline Electrostatic2DState& Electrostatic2DState::operator =(const Electrostatic2DState &U) {
  E.x = U.E.x; E.y = U.E.y; V = U.V;
  return *this;
}

/**********************************************************************
 * Electrostatic2DState -- Unary arithmetic operators.                *
 **********************************************************************/
//inline Electrostatic2DState operator +(const Electrostatic2DState &U) {
//return U;
//}

inline Electrostatic2DState operator -(const Electrostatic2DState &U) {
  return Electrostatic2DState(U.E.x,-U.E.y,-U.V);
}

/**********************************************************************
 * Electrostatic2DState -- Shortcut arithmetic operators.             *
 **********************************************************************/
inline Electrostatic2DState& Electrostatic2DState::operator +=(const Electrostatic2DState &U) {
  E.x += U.E.x; E.y += U.E.y; V += U.V;
  return *this;
}

inline Electrostatic2DState& Electrostatic2DState::operator -=(const Electrostatic2DState &U) {
  E.x -= U.E.x; E.y -= U.E.y; V -= U.V;
  return *this;
}

inline Electrostatic2DState& Electrostatic2DState::operator *=(const double &a) {
  E.x *= a; E.y *= a; V *= a;
  return *this;
}

inline Electrostatic2DState& Electrostatic2DState::operator /=(const double &a) {
  E.x /= a; E.y /= a; V /= a;
  return *this;
}

/**********************************************************************
 * Electrostatic2DState -- Relational operators.                      *
 **********************************************************************/
inline int operator ==(const Electrostatic2DState &U1, const Electrostatic2DState &U2) {
  return (U1.E == U2.E && U1.V == U2.V);
}

inline int operator !=(const Electrostatic2DState &U1, const Electrostatic2DState &U2) {
  return (U1.E != U2.E || U1.V != U2.V);
}

/**********************************************************************
 * Electrostatic2DState -- Input-output operators.                    *
 **********************************************************************/
inline ostream &operator << (ostream &out_file, const Electrostatic2DState &U) {
  out_file.setf(ios::scientific);
  out_file << " " << U.E.x << " " << U.E.y << " " << U.V;
  out_file.unsetf(ios::scientific);
  return out_file;
}

inline istream &operator >> (istream &in_file, Electrostatic2DState &U) {
  in_file.setf(ios::skipws);
  in_file >> U.E.x >> U.E.y >> U.V;
  in_file.unsetf(ios::skipws);
  return in_file;
}

/**********************************************************************
 * Useful 2D electrostatic state constants.                           *
 **********************************************************************/
const Electrostatic2DState ELECTROSTATIC2D_VACUUM(ZERO,ZERO,ZERO);
const Electrostatic2DState ELECTROSTATIC2D_ONE(ONE,ONE,ONE);

/**********************************************************************
 * Electrostatic2DState -- External subroutines.                      *
 **********************************************************************/

extern Electrostatic2DState Rotate(const Electrostatic2DState &U,
				   const Vector2D &norm_dir);

extern Electrostatic2DState Translate(const Electrostatic2DState &U,
				      const Vector2D &V);

extern Electrostatic2DState Reflect(const Electrostatic2DState &U,
				    const Vector2D &norm_dir);

extern Electrostatic2DState Mirror(const Electrostatic2DState &U,
				   const Vector2D &norm_dir);

extern Electrostatic2DState FluxArithmetic_n(const Vector2D &X,
					     const Vector2D &X1,
					     const Electrostatic2DState &U1,
					     const Electrostatic2DState &dU1dx,
					     const Electrostatic2DState &dU1dy,
					     const Vector2D &X2,
					     const Electrostatic2DState &U2,
					     const Electrostatic2DState &dU2dx,
					     const Electrostatic2DState &dU2dy,
					     const Vector2D &norm_dir);

extern Electrostatic2DState FluxDiamondPath_n(const Vector2D &X,
					      const Vector2D &X1, const Electrostatic2DState &U1,
					      const Vector2D &X2, const Electrostatic2DState &U2,
					      const Vector2D &X3, const Electrostatic2DState &U3,
					      const Vector2D &X4, const Electrostatic2DState &U4,
					      const Vector2D &norm_dir,
					      const int &diamond_path);

#endif // _ELECTROSTATIC2D_STATE_INCLUDED
