/*!\file New_AdvectDiffuse2DState.h
  \brief Temporary header file defining 2D Advection Diffusion Equation Solution State Class. */

#ifndef _NEW_ADVECTDIFFUSE2D_STATE_INCLUDED
#define _NEW_ADVECTDIFFUSE2D_STATE_INCLUDED

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
#include "../Math/Math.h" /* Include math macro header files. */
#include "../CFD/CFD.h"   /* Include CFD header files. */
#include "../Math/Vector2D.h" /* Include vector header files. */
#include "AdvectDiffuse2DExactSolutions.h" /* Include 2D advection diffusion exact solutions header file */
#include "AdvectDiffuse2DParameterFields.h" /* Include 2D advection diffusion parameter fields */

/* Define the classes. */

#define	NUM_VAR_ADVECTDIFFUSE2D    1

/*!
 * \class AdvectDiffuse2D_State_New
 *
 * @brief Solution state class definition for the 2D advection-diffusion
 *        equation.
 *
 * Solution state class definition for the 2D advection-diffusion
 * equation.
 *
 * \verbatim
 * Member variables
 *     u        -- Return solution.
 *     V        -- Return advection velocity field
 *     k        -- Return diffusion coefficient field
 *
 * Member functions
 *    Fa        -- Return advective flux.
 *    Fd        -- Return diffusive flux.
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
 * U <= U;
 * U >= U;
 * cout << U; (output function)
 * cin  >> U; (input function)
 * \endverbatim
 */
class AdvectDiffuse2D_State_New{
public:
  //! @name Defined datatypes
  //@{
  typedef VelocityFields::VelocityFieldType AdvectionVelocityType;
  typedef DiffusionFields::NonlinearDiffusionFieldType DiffusionFieldType;
  //@}

  //! @name Solution state variables and associated constants:
  //@{
  double  u;   // Solution.
  static AdvectionVelocityType  V;      // Advection velocity as function of 'x' and 'y'.
  static DiffusionFieldType k;		// Non-linear diffusion coefficient as function of 'x', 'y' and 'u'
  //@}

  //! @name Creation and copy constructors.
  //@{
  //! Default constructor.
  AdvectDiffuse2D_State_New(void): u(0.0){ }

  //! Value Constructor
  AdvectDiffuse2D_State_New(const double &uu): u(uu){ };

  //! Copy constructor.
  AdvectDiffuse2D_State_New(const AdvectDiffuse2D_State_New &U){ u = U.u; }
  //@}

  //! Destructor.
  ~AdvectDiffuse2D_State_New(void){}

  //! @name Useful operators.
  //@{
  //! Vacuum/zero operator.
  void Vacuum(void) { u = ZERO; }
  //@}

  //! @name Advective Flux.
  //@{
  Vector2D Fa(const double &x, const double &y) const;
  friend Vector2D Fa(const double &x, const double &y, const AdvectDiffuse2D_State_New &U);
  //@}

  //! @name Diffusive flux.
  //@{
  Vector2D Fd(const double &x, const double &y,
	      const Vector2D &GradU) const;
  friend Vector2D Fd(const double &x, const double &y,
		     const AdvectDiffuse2D_State_New &U, const Vector2D &GradU);
  Vector2D Fd(const double &x, const double &y,
	      const double &dudx, const double &dudy) const;
  friend Vector2D Fd(const double &x, const double &y,
		     const AdvectDiffuse2D_State_New &U, const double &dudx, const double &dudy);
  //@}

  //! @name Regular source term.
  //@{
  //@}

  //! @name Axisymmetric source term.
  //@{
  //@}

  /*! @name Assignment operator. */
  //@{
  AdvectDiffuse2D_State_New operator= (const AdvectDiffuse2D_State_New &rhs);
  //@}

  /*! @name Index operators */
  //@{
  double &operator[](const int & index);
  const double &operator[](const int & index) const;
  //@}

  //! @name Binary arithmetic operators.
  //@{
  AdvectDiffuse2D_State_New operator +(const AdvectDiffuse2D_State_New &U) const;
  AdvectDiffuse2D_State_New operator -(const AdvectDiffuse2D_State_New &U) const;
  AdvectDiffuse2D_State_New operator *(const AdvectDiffuse2D_State_New &U) const;
  AdvectDiffuse2D_State_New operator *(const double &a) const;
  friend AdvectDiffuse2D_State_New operator *(const double &a, const AdvectDiffuse2D_State_New &U){ return U*a;}
  AdvectDiffuse2D_State_New operator /(const double &a){ return *this * (1.0/a); }
  //@}

  //! @name Unary arithmetic operators.
  //@{
  friend AdvectDiffuse2D_State_New & operator +(const AdvectDiffuse2D_State_New &U);
  friend AdvectDiffuse2D_State_New & operator -(const AdvectDiffuse2D_State_New &U);
  //@}

  //! @name Shortcut arithmetic operators.
  //@{
  AdvectDiffuse2D_State_New & operator +=(AdvectDiffuse2D_State_New &U);
  AdvectDiffuse2D_State_New & operator -=(AdvectDiffuse2D_State_New &U);
  AdvectDiffuse2D_State_New & operator *=(const double &a);
  AdvectDiffuse2D_State_New & operator /=(const double &a);
  //@}

  //! @name Relational operators.
  //@{
  friend bool operator ==(const AdvectDiffuse2D_State_New& lhs, const AdvectDiffuse2D_State_New& rhs);
  friend bool operator !=(const AdvectDiffuse2D_State_New& lhs, const AdvectDiffuse2D_State_New& rhs);
  friend bool operator >=(const AdvectDiffuse2D_State_New& lhs, const AdvectDiffuse2D_State_New& rhs);
  friend bool operator <=(const AdvectDiffuse2D_State_New& lhs, const AdvectDiffuse2D_State_New& rhs);
  //@}

  //! @name Input-output operators.
  //@{
  friend ostream &operator << (ostream &out_file, const AdvectDiffuse2D_State_New &U);
  friend istream &operator >> (istream &in_file,  AdvectDiffuse2D_State_New &U);
  //@}
};

/*!
 * Compute the advective flux \f$Fa=V(x,y)*u\f$
 * 
 * \param [in] x x-coordinate location
 * \param [in] y y-coordinate location
 * \return advective flux at a given location
 */
inline Vector2D AdvectDiffuse2D_State_New::Fa(const double &x, const double &y) const {
  return V(x,y)*u;
}

/*!
 * Compute the advective flux \f$Fa=V(x,y)*u\f$
 *
 * \param [in] x x-coordinate location
 * \param [in] y y-coordinate location
 * \param [in] U solution state used to compute the advective flux
 * \return advective flux at a given location based on the solution state U
 */
inline Vector2D Fa(const double &x, const double &y,
		   const AdvectDiffuse2D_State_New &U) {
  return U.V(x,y)*U.u;
}

/*!
 * Compute the diffusive flux \f$Fd=-k(x,y,u)*\nabla u\f$, where 'u' is the solution scalar.
 * 
 * \param [in] x     x-coordinate location
 * \param [in] y     y-coordinate location
 * \param [in] GradU solution gradient vector
 * \return diffusive flux at a given location
 */
inline Vector2D AdvectDiffuse2D_State_New::Fd(const double &x, const double &y,
					      const Vector2D &GradU) const {
  return -k(x,y,u)*GradU;
}

/*!
 * Compute the diffusive flux \f$Fd=-k(x,y,u)*\nabla u\f$, where 'u' is the solution scalar.
 * 
 * \param [in] x     x-coordinate location
 * \param [in] y     y-coordinate location
 * \param [in] GradU solution gradient vector
 * \param [in] U     solution state used to compute the advective flux
 * \return diffusive flux at a given location
 */
inline Vector2D Fd(const double &x, const double &y,
		   const AdvectDiffuse2D_State_New &U,
		   const Vector2D &GradU) {
  return -U.k(x,y,U.u)*GradU;
}

/*!
 * Compute the diffusive flux \f$Fd=-k(x,y,u)*\nabla u\f$, where 'u' is the solution scalar.
 * 
 * \param [in] x    x-coordinate location
 * \param [in] y    y-coordinate location
 * \param [in] dudx solution gradient with respect to x-direction
 * \param [in] dudy solution gradient with respect to y-direction
 * \return diffusive flux at a given location
 */
inline Vector2D AdvectDiffuse2D_State_New::Fd(const double &x, const double &y,
					      const double &dudx, const double &dudy) const {
  return Fd(x,y,Vector2D(dudx,dudy));
}

/*!
 * Compute the diffusive flux \f$Fd=-k(x,y)*\nabla u\f$.
 * 
 * \param [in] x    x-coordinate location
 * \param [in] y    y-coordinate location
 * \param [in] dudx solution gradient with respect to x-direction
 * \param [in] dudy solution gradient with respect to y-direction
 * \param [in] U    solution state used to compute the advective flux
 * \return diffusive flux at a given location
 */
inline Vector2D Fd(const double &x, const double &y,
		   const AdvectDiffuse2D_State_New &U,
		   const double &dudx, const double &dudy) {
  return Fd(x,y,U,Vector2D(dudx,dudy));
}

/*!
 * Assignment operator
 */
inline AdvectDiffuse2D_State_New AdvectDiffuse2D_State_New::operator= (const AdvectDiffuse2D_State_New &rhs){
  if(this == &rhs) return *this;
  u = rhs.u;
  return *this;
}

//! index operator 
inline double & AdvectDiffuse2D_State_New::operator[](const int & index){
  assert( index >= 1 && index <= NUM_VAR_ADVECTDIFFUSE2D );
  switch(index) {
  case 1 :
    return (u);
  default:
    return (u);
  };
}

//! constant index operator 
inline const double & AdvectDiffuse2D_State_New::operator[](const int & index) const {
  assert( index >= 1 && index <= NUM_VAR_ADVECTDIFFUSE2D );
  switch(index) {
  case 1 :
    return (u);
  default:
    return (u);
  };
}

//! summation between states
inline AdvectDiffuse2D_State_New AdvectDiffuse2D_State_New::operator +(const AdvectDiffuse2D_State_New &U) const {
  
}

//! subtraction between states
inline AdvectDiffuse2D_State_New AdvectDiffuse2D_State_New::operator -(const AdvectDiffuse2D_State_New &U) const {

}

//! multiplication between states
inline AdvectDiffuse2D_State_New AdvectDiffuse2D_State_New::operator *(const AdvectDiffuse2D_State_New &U) const {

}

//! multiplication with scalar
inline AdvectDiffuse2D_State_New AdvectDiffuse2D_State_New::operator *(const double &a) const{

}

//! unary arithmetic + operator
inline AdvectDiffuse2D_State_New & operator +(const AdvectDiffuse2D_State_New &U){
  
}

//! unary arithmetic - operator
inline AdvectDiffuse2D_State_New & operator -(const AdvectDiffuse2D_State_New &U){
  
}

//! shortcut summation
inline AdvectDiffuse2D_State_New & AdvectDiffuse2D_State_New::operator +=(AdvectDiffuse2D_State_New &U){

}

//! shortcut subtraction
inline AdvectDiffuse2D_State_New & AdvectDiffuse2D_State_New::operator -=(AdvectDiffuse2D_State_New &U){
  
}

//! shortcut multiplication with scalar
inline AdvectDiffuse2D_State_New & AdvectDiffuse2D_State_New::operator *=(const double &a){
  
}

//! shortcut division with scalar
inline AdvectDiffuse2D_State_New & AdvectDiffuse2D_State_New::operator /=(const double &a){

}

//! equal operator 
inline bool operator ==(const AdvectDiffuse2D_State_New &lhs, const AdvectDiffuse2D_State_New &rhs) {
  return (lhs.u == rhs.u);
}

//! not-equal operator
inline bool operator !=(const AdvectDiffuse2D_State_New& lhs, const AdvectDiffuse2D_State_New& rhs){
  return !(lhs == rhs);
}

//! greater or equal operator
inline bool operator >=(const AdvectDiffuse2D_State_New& lhs, const AdvectDiffuse2D_State_New& rhs){
  return (lhs.u >= rhs.u);
}

//! less or equal operator 
inline bool operator <=(const AdvectDiffuse2D_State_New& lhs, const AdvectDiffuse2D_State_New& rhs){
  return (lhs.u <= rhs.u);
}

//! output operator
inline ostream &operator << (ostream &out_file, const AdvectDiffuse2D_State_New &U) {
  out_file.setf(ios::scientific);
  out_file << " " << U.u ;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, AdvectDiffuse2D_State_New &U) {
  in_file.setf(ios::skipws);
  in_file >> U.u ;
  in_file.unsetf(ios::skipws);
  return (in_file);
}


#endif /* _DVECTDIFFUSE2D_STATE_INCLUDED  */
