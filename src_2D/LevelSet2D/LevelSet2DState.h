/**********************************************************************
 * LevelSet2DState.h: Header file defining 2D LevelSet Solution State *
 *                    Classes.                                        *
 **********************************************************************/

#ifndef _LEVELSET2D_STATE_INCLUDED
#define _LEVELSET2D_STATE_INCLUDED

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

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _VECTOR2D_INCLUDED
#include "../Math/Vector2D.h"
#endif //_VECTOR2D_INCLUDED

// Define the classes.

#define	NUM_VAR_LEVELSET2D    4

/*!
 * Class: LevelSet2DState
 *
 * @brief Level set variable solution state class definition.
 *
 * \verbatim
 * Member functions
 *     psi   -- Level set function.
 *     F     -- Extended front speed.
 *     V     -- Bulk velocity flow-field.
 *
 * Member operators
 *     U -- the level set solution state
 *     c -- a scalar (double)
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
 * U *= a;
 * U /= a;
 * U == U;
 * U != U;
 * cout << U; (output function)
 * cin  >> U; (input function)
 * \endverbatim
 */
class LevelSet2DState {
 private:
 public:
  double   psi; //!< Level set function.
  double     F; //!< Extended front speed.
  Vector2D   V; //!< Convection/bulk velocity flow.

  //@{ @name Constructors
  //! Creation constructor.
  LevelSet2DState(void) {
    psi = ZERO; F = ONE; V = Vector2D_ZERO;
  }
  
  //! Copy constructor.
  LevelSet2DState(const LevelSet2DState &U) {
    psi = U.psi; F = U.F; V = U.V;
  }
  
  //! Copy constructor.
  LevelSet2DState(LevelSet2DState &U) {
    psi = U.psi; F = U.F; V = U.V;
  }
  
  //! Copy constructor.
  LevelSet2DState(const double &psipsi,
		  const double &FF,
		  const Vector2D &VV) {
    psi = psipsi; F = FF; V = VV;
  }

  //! Copy constructor.
  LevelSet2DState(const double &psipsi,
		  const double &FF,
		  const double &Vx,
		  const double &Vy) {
    psi = psipsi; F = FF; V = Vector2D(Vx,Vy);
  }

  // Destructor.
  // ~LevelSet2DState(void);
  // Use automatically generated destructor.

  //@}

  // Assignment operator.
  // LevelSet2DState operator = (const LevelSet2DState &U);
  // Use automatically generated assignment operator.
  
  //! Copy functions.
  void Copy(const LevelSet2DState &U);

  //! Vacuum operator.
  void Vacuum(void);
 
  //! Constant operator.
  void Constant(const double &val);

  //@{ @name Index operator.
  double &operator[](int index) {
    //assert(index >= 1 && index <= NUM_VAR_LEVELSET2D);
    switch(index) {
    case 1 :
      return psi;
    case 2 :
      return F;
    case 3 :
      return V.x;
    case 4 :
      return V.y;
    default:
      return psi;
    };
  }
  
  const double &operator[](int index) const {
    //assert(index >= 1 && index <= NUM_VAR_LEVELSET2D);
    switch(index) {
    case 1 :
      return psi;
    case 2 :
      return F;
    case 3 :
      return V.x;
    case 4 :
      return V.y;
    default:
      return psi;
    };
  }
  //@}

  //@{ Binary arithmetic operators.
  LevelSet2DState operator +(const LevelSet2DState &U) const;
  LevelSet2DState operator -(const LevelSet2DState &U) const;
  double operator *(const LevelSet2DState &U) const;
  LevelSet2DState operator *(const double &a) const;
  friend LevelSet2DState operator *(const double &a, const LevelSet2DState &U);
  LevelSet2DState operator /(const double &a) const;
  LevelSet2DState operator ^(const LevelSet2DState &U) const;
  //@}

  //@{ Assignment Operator.
  LevelSet2DState& operator =(const LevelSet2DState &U);
  //@}

  //@{ Unary arithmetic operators.
  //LevelSet2DState operator +(const LevelSet2DState &U);
  friend LevelSet2DState operator -(const LevelSet2DState &U);
  //@}

  //@{ Shortcut arithmetic operators.
  LevelSet2DState &operator +=(const LevelSet2DState &U);
  LevelSet2DState &operator -=(const LevelSet2DState &U);
  LevelSet2DState &operator *=(const double &a);
  LevelSet2DState &operator /=(const double &a);
  //@}

  //@{ Relational operators.
  friend int operator ==(const LevelSet2DState &U1, const LevelSet2DState &U2);
  friend int operator !=(const LevelSet2DState &U1, const LevelSet2DState &U2);
  //@}

  //@{ Input-output operators.
  friend ostream &operator << (ostream &out_file, const LevelSet2DState &U);
  friend istream &operator >> (istream &in_file,  LevelSet2DState &U);
  //@}

};

/**********************************************************************
 * LevelSet2DState::Copy -- Copy the level state.                     *
 **********************************************************************/
inline void LevelSet2DState::Copy(const LevelSet2DState &U) {
  psi = U.psi; F = U.F; V = U.V;
}

/**********************************************************************
 * LevelSet2DState::Vacuum -- Set vacuum state.                       *
 **********************************************************************/
inline void LevelSet2DState::Vacuum(void) {
  psi = ZERO; F = ZERO; V = Vector2D_ZERO;
}

/**********************************************************************
 * LevelSet2DState::Constant -- Set constant state.                   *
 **********************************************************************/
inline void LevelSet2DState::Constant(const double &val) {
  psi = val; F = val; V.x = val; V.y = val;
}

/**********************************************************************
 * LevelSet2DState -- Useful operators.                               *
 **********************************************************************/
//inline LevelSet2DState fabs(LevelSet2DState &U) {
//  return LevelSet2DState(fabs(U.psi),fabs(U.F),fabs(V));
//}

//inline LevelSet2DState max(LevelSet2DState &U1, LevelSet2DState &U2) {
//  return LevelSet2DState(max(U1.psi,U2.psi),max(U1.F,U2.F),max(U1.V,U2.V));
//}

//inline LevelSet2DState sqr(LevelSet2DState &U) {
//  return LevelSet2DState(sqr(U.psi),sqr(U.F),);
//}

//inline LevelSet2DState sqrt(LevelSet2DState &U) {
//  return LevelSet2DState(sqrt(U.psi),sqrt(U.F));
//}

/**********************************************************************
 * LevelSet2DState -- Binary arithmetic operators.                    *
 **********************************************************************/
inline LevelSet2DState LevelSet2DState::operator +(const LevelSet2DState &U) const {
  LevelSet2DState Utemp;
  Utemp.Copy(*this);
  Utemp += U;
  return Utemp;
}

inline LevelSet2DState LevelSet2DState::operator -(const LevelSet2DState &U) const {
  LevelSet2DState Utemp;
  Utemp.Copy(*this);
  Utemp -= U;
  return Utemp;
}

// Inner product operator.
inline double LevelSet2DState::operator *(const LevelSet2DState &U) const {
  return psi*U.psi + F*U.F + V*U.V;
}

inline LevelSet2DState LevelSet2DState::operator *(const double &a) const {
  LevelSet2DState Utemp;
  Utemp.Copy(*this);
  Utemp.psi = psi*a; Utemp.F = F*a; Utemp.V = V*a;
  return Utemp;
}

inline LevelSet2DState operator *(const double &a, const LevelSet2DState &U) {
  return LevelSet2DState(a*U.psi,a*U.F,a*U.V);
}

inline LevelSet2DState LevelSet2DState::operator /(const double &a) const {
  LevelSet2DState Utemp;
  Utemp.Copy(*this);
  Utemp.psi = psi/a; Utemp.F = F/a; Utemp.V = V/a;
  return Utemp;
}

inline LevelSet2DState LevelSet2DState::operator ^(const LevelSet2DState &U) const {
  LevelSet2DState Utemp;
  Utemp.Copy(*this);
  Utemp.psi = psi*U.psi; Utemp.F   = F*U.F; Utemp.V.x = V.x*U.V.x; Utemp.V.y = V.y*U.V.y;
  return Utemp;
}

/**********************************************************************
 * LevelSet2DState -- Assignment operator.                            *
 **********************************************************************/
inline LevelSet2DState& LevelSet2DState::operator =(const LevelSet2DState &U) {
  if (this != &U) { psi = U.psi; F = U.F; V = U.V; }
  return *this;
}

/**********************************************************************
 * LevelSet2DState -- Unary arithmetic operators.                     *
 **********************************************************************/
//inline LevelSet2DState operator +(const LevelSet2DState &U) {
//return U;
//}

inline LevelSet2DState operator -(const LevelSet2DState &U) {
  LevelSet2DState Utemp;
  Utemp.psi = -U.psi; Utemp.F = -U.F; Utemp.V = -U.V;
  return Utemp;
}

/**********************************************************************
 * LevelSet2DState -- Shortcut arithmetic operators.                  *
 **********************************************************************/
inline LevelSet2DState& LevelSet2DState::operator +=(const LevelSet2DState &U) {
  psi += U.psi; F += U.F; V += U.V;
  return *this;
}

inline LevelSet2DState& LevelSet2DState::operator -=(const LevelSet2DState &U) {
  psi -= U.psi; F -= U.F; V -= U.V;
  return *this;
}

inline LevelSet2DState& LevelSet2DState::operator *=(const double &a) {
  psi *= a; F *= a; V.x *= a; V.y *= a;
  return *this;
}

inline LevelSet2DState& LevelSet2DState::operator /=(const double &a) {
  psi /= a; F /= a; V.x /= a; V.y /= a;
  return *this;
}

/**********************************************************************
 * LevelSet2DState -- Relational operators.                           *
 **********************************************************************/
inline int operator ==(const LevelSet2DState &U1, const LevelSet2DState &U2) {
  return (U1.psi == U2.psi && U1.F == U2.F && U1.V == U2.V);
}

inline int operator !=(const LevelSet2DState &U1, const LevelSet2DState &U2) {
  return (U1.psi != U2.psi || U1.F != U2.F || U1.V != U2.V);
}

/**********************************************************************
 * LevelSet2DState -- Input-output operators.                         *
 **********************************************************************/
inline ostream &operator << (ostream &out_file, const LevelSet2DState &U) {
  out_file.setf(ios::scientific);
  out_file << " " << U.psi   << " " << U.F << " " << U.V.x << " " << U.V.y;
  out_file.unsetf(ios::scientific);
  return out_file;
}

inline istream &operator >> (istream &in_file, LevelSet2DState &U) {
  in_file.setf(ios::skipws);
  in_file >> U.psi >> U.F >> U.V.x >> U.V.y;
  in_file.unsetf(ios::skipws);
  return in_file;
}

/**********************************************************************
 * Useful 2D Level Set state constants.                               *
 **********************************************************************/
const LevelSet2DState LevelSet2D_ZERO(ZERO,ZERO,Vector2D_ZERO);

#endif // _LEVELSET2D_STATE_INCLUDED
