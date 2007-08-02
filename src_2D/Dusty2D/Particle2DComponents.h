/**********************************************************************
 * Particle2DComponents.h: Header file defining the 2D particle       *
 *                         components solution classes.               *
 **********************************************************************/

#ifndef _PARTICLE2D_COMPONENTS_INCLUDED
#define _PARTICLE2D_COMPONENTS_INCLUDED

// Include required C++ libraries.
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

// Include CFD header file.

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

// Include the particle state class header file.

#ifndef _PARTICLE2D_STATE_INCLUDED
#include "Particle2DState.h"
#endif // _PARTICLE2D_STATE_INCLUDED

#define PARTICLE2D_SINGLE_VELOCITY_FORMULATION  1
#define PARTICLE2D_MULTI_VELOCITY_FORMULATION   4

// Define the classes.

class Particle2D_cComponents;

/*!
 * Class: Particle2D_pComponents
 *
 * @brief Primitive variable component definition for the particle-phase
 * of a gas-particle flow.
 *
 * Primitive variable component definition for the particle-phase of a
 * coupled gas-particle flow.  This class aids in the organization and
 * manipulation of the desired number of components (1 or 4) to be used 
 * in the solution of the coupled gas-particle flow.  Standard methods
 * correspond to the standard Eulerian approach.  The use of four 
 * components corresponds to the multi-velcoity solution method for an
 * Eulerian formulation of a particle-phase.
 *
 * \verbatim
 * Member functions
 *     W -- Array of primitive solution state components.
 *     NUM_CMP_PARTICLE2D -- Number of particle component groups.
 *     NUM_VAR_PARTICLE2D -- Total number of variables.
 *
 *     U        -- Return the conserved solution state.
 *     dUdW     -- Return the Jacobian of the conserved solution
 *                 variables with respect to the primitive solution
 *                 variables.
 *     dWdU     -- Return the Jacobian of the primitive solution
 *                 variables with respect to the conserved solution
 *                 variables.
 *
 *     F        -- Return x-direction solution flux.
 *     dFdU     -- Return the Jacobian of the solution flux vector with
 *                 respect to the conserved solution variables.
 *
 *     lambda_x -- Return x-direction eigenvalue(s).
 *     rp_x     -- Return primitive right eigenvector (x-direction).
 *     rc_x     -- Return conserved right eigenvector (x-direction).
 *     lp_x     -- Return primitive left eigenvector (x-direction).
 *
 *     S        -- Return axisymmetric source term vector.
 *     dSadU    -- Return the Jacobian of the axisymmetric source term
 *                 vector with respect to the conserved solution variables.
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
class Particle2D_pComponents {
 private:
 public:
  //@{ @name Particle-phase primitive variables:
  Particle2D_pState          *W; //!< Array of particle primitive states.
  static int NUM_CMP_PARTICLE2D; //!< Number of particle component groups.
  //@}

  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  Particle2D_pComponents(void) { W = NULL; allocate(); }

  //! Copy constructor.
  Particle2D_pComponents(const Particle2D_pComponents &Wc);

  //! Copy constructor.
  Particle2D_pComponents(const Particle2D_cComponents &Uc, const double &cm);

  //! Assignment constructor.
  Particle2D_pComponents(const Particle2D_pState &Wp);
  
  //! Assignment constructor.
  Particle2D_pComponents(const Particle2D_cState &Up, const double &cm);

  //! Destructor.
  ~Particle2D_pComponents(void) { deallocate(); }
  //@}

  //@{ @name Allocation and deallocation functions.
  //! Allocation function.
  void allocate(void) { deallocate(); W = new Particle2D_pState[NUM_CMP_PARTICLE2D]; }

  //! Deallocation function.
  void deallocate(void) { if (W != NULL) { delete []W; W = NULL; } }
  //@}

  //@{ @name Static variable functions.
  void set_particle_components(void);
  void set_particle_components(const int Nc);
  //@}

  //@{ @name Useful operators.
  //! Copy operator.
  void Copy(const Particle2D_pComponents &Wc);

  //! Vacuum operator.
  void Vacuum(void);

  //! Standard atmosphere operator.
  void Standard_Atmosphere(void);

  //! Check for unphysical state properties.
  int Unphysical_Properties(void) const;

  //! Reset unphysical state properties.
  void Reset_Unphysical_Properties(void);
  //@}

  //@{ @name Mass-averaged/total particle-phase quantities.
  //! Particle-phase concentration (total over the components).
  double sigma(void) const;

  //! Particle-phase number density (total over the components).
  double nd(const double &mp) const;

  //! Particle-phase velocity (mass-averaged over the components).
  Vector2D u(void) const;

  //! Particle-phase momentum (mass-averaged over the components).
  Vector2D du(void) const;

  //! Particle-phase temperature (mass-averaged over the components).
  double Tp(void) const;

  //! Particle-phase total thermal/internal energy (mass-averaged over the components).
  double ep(const double &cm) const;
  
  //! Particle-phase total energy (mass-averaged over the components).
  double Ep(const double &cm) const;

  //! Particle-phase current (mass-averaged over the components).
  Vector2D jc(const double &qe, const double &mp) const;
  //@}

  //@{ @name Conserved solution state.
  Particle2D_cComponents U(const double &cm) const;
  Particle2D_cComponents U(const Particle2D_pComponents &Wc, const double &cm) const;
  //@}

  //! Jacobian of the conserved solution variables with respect to the
  //! primitive solution variables.
  void dUdW(DenseMatrix &dUdW, const int &n, const double &cm) const;

  //! Jacobian of the primitive solution variables with respect to the
  //! conserved solution variables.
  void dWdU(DenseMatrix &dWdU, const int &n, const double &cm) const;

  //@{ @name Solution flux (x-direction) and Jacobian.
  Particle2D_cComponents F(const double &cm) const;
  Particle2D_cComponents F(const Vector2D &V, const double &cm) const;
  void dFdU(DenseMatrix &dFdU, const int &n, const double &cm) const;
  void dFdU(DenseMatrix &dFdU, const int &n, const Vector2D &V, const double &cm) const;
  //@}

  //@{ @name Eigenstructure.
  //! Eigenvalue(s) (x-direction).
  Particle2D_pComponents lambda_x(void) const;

  //! Eigenvalue(s) for a moving reference frame (x-direction).
  Particle2D_pComponents lambda_x(const Vector2D &V) const;

  //! Primitive right eigenvector (x-direction).
  Particle2D_pComponents rp_x(int index) const;

  //! Conserved right eigenvector (x-direction).
  Particle2D_cComponents rc_x(int index) const;

  //! Primitive left eigenvector (x-direction).
  //Particle2D_pComponents lp_x(int index) const;
  //@}

  //@{ Axisymmetric flow source vector and Jacobian.
  Particle2D_cComponents Sa(const Vector2D &X, const double &cm) const;
  void dSadU(DenseMatrix &dSadU, const Vector2D &X, const int &n, const double &cm) const;
  //@}

  //@{ @name Electrostatic-force source term vector and Jacobian.
  Particle2D_cComponents Se(const Vector2D &E, const double &Cm) const;
  void dSedU(DenseMatrix &dSedU, const int &n, const Vector2D &E, const double &Cm) const;
  //@}

  //@{ @name Index operator.
  Particle2D_pState &operator[](int index) {
    //assert(index >= 0 && index < NUM_CMP_PARTICLE2D);
    return W[index];
  }
  
  const Particle2D_pState &operator[](int index) const {
    //assert(index >= 0 && index < NUM_CMP_PARTICLE2D);
    return W[index];
  }

  void get_indices(int index, int &cindex, int &pindex) {
    if (NUM_CMP_PARTICLE2D == PARTICLE2D_SINGLE_VELOCITY_FORMULATION) cindex = 0;
    else cindex = (index%(NUM_CMP_PARTICLE2D*NUM_VAR_PARTICLE2D))/NUM_CMP_PARTICLE2D;
    pindex = index%NUM_VAR_PARTICLE2D + 1;
  }

  void get_indices(int index, int &cindex, int &pindex) const {
    if (NUM_CMP_PARTICLE2D == PARTICLE2D_SINGLE_VELOCITY_FORMULATION) cindex = 0;
    else cindex = (index%(NUM_CMP_PARTICLE2D*NUM_VAR_PARTICLE2D))/NUM_CMP_PARTICLE2D;
    pindex = index%NUM_VAR_PARTICLE2D + 1;
  }
  //@}

  //@{ @name Binary arithmetic operators.
  Particle2D_pComponents operator +(const Particle2D_pComponents &Wc) const;
  Particle2D_pComponents operator -(const Particle2D_pComponents &Wc) const;
  double operator *(const Particle2D_pComponents &W1) const;
  Particle2D_pComponents operator *(const double &a) const;
  friend Particle2D_pComponents operator *(const double &a, const Particle2D_pComponents &Wc);
  Particle2D_pComponents operator /(const double &a) const;
  Particle2D_pComponents operator ^(const Particle2D_pComponents &Wc) const;
  //@}

  //@{ @name Assignment operator.
  Particle2D_pComponents& operator =(const Particle2D_pComponents &Wc);
  //@}

  //@{ @name Unary arithmetic operators.
  //Particle2D_pComponents operator +(const Particle2D_pComponents &W);
  friend Particle2D_pComponents operator -(const Particle2D_pComponents &Wc);
  //@}

  //@{ @name Shortcut arithmetic operators.
  Particle2D_pComponents &operator +=(const Particle2D_pComponents &Wc);
  Particle2D_pComponents &operator -=(const Particle2D_pComponents &Wc);
  Particle2D_pComponents &operator *=(const double &a);
  Particle2D_pComponents &operator /=(const double &a);
  //@}

  //@{ @name Relational operators.
  friend int operator ==(const Particle2D_pComponents &W1, const Particle2D_pComponents &W2);
  friend int operator !=(const Particle2D_pComponents &W1, const Particle2D_pComponents &W2);
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const Particle2D_pComponents &Wc);
  friend istream &operator >> (istream &in_file, Particle2D_pComponents &Wc);
  //@}

};

/*!
 * Class: Particle2D_cComponents
 *
 * @brief Conserved variable component definition for the particle-phase
 * of a gas-particle flow.
 *
 * Conserved variable component definition for the particle-phase of a
 * coupled gas-particle flow.  This class aids in the organization and
 * manipulation of the desired number of components (1 or 4) to be used 
 * in the solution of the coupled gas-particle flow.  Standard methods
 * correspond to the standard Eulerian approach.  The use of four 
 * components corresponds to the multi-velcoity solution method for an
 * Eulerian formulation of a particle-phase.
 *
 * \verbatim
 * Member functions
 *     U -- Array of conservative solution state components.
 *     NUM_CMP_PARTICLE2D -- Number of particle component groups.
 *
 *     W        -- Return primitive solution state.
 *     dUdW     -- Return the Jacobian of the conserved solution
 *                 variables with respect to the primitive solution
 *                 variables.
 *     dWdU     -- Return the Jacobian of the primitive solution
 *                 variables with respect to the conserved solution
 *                 variables.
 *
 *     F        -- Return x-direction solution flux.
 *     dFdU     -- Return the Jacobian of the solution flux vector with
 *                 respect to the conserved solution variables.
 *
 *     lambda_x -- Return x-direction eigenvalue(s).
 *     rp_x     -- Return primitive right eigenvector (x-direction).
 *     rc_x     -- Return conserved right eigenvector (x-direction).
 *     lp_x     -- Return primitive left eigenvector (x-direction).
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
class Particle2D_cComponents {
  private:
  public:
  //@{ @name Particle-phase conservative variables:
  Particle2D_cState          *U; //!< Array of particle conserved states.
  static int NUM_CMP_PARTICLE2D; //!< Number of particle component groups.
  //@}

  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  Particle2D_cComponents(void) { U = NULL; allocate(); }

  //! Copy constructor.
  Particle2D_cComponents(const Particle2D_cComponents &Uc);

  //! Copy constructor.
  Particle2D_cComponents(const Particle2D_pComponents &Wc, const double &cm);

  //! Assignment constructor.
  Particle2D_cComponents(const Particle2D_cState &Up);

  //! Assignment constructor.
  Particle2D_cComponents(const Particle2D_pState &Wp, const double &cm);

  //! Destructor.
  ~Particle2D_cComponents(void) { deallocate(); };
  //@}

  //@{ @name Allocation and deallocation functions.
  //! Allocation function.
  void allocate(void) { deallocate(); U = new Particle2D_cState[NUM_CMP_PARTICLE2D]; }

  //! Deallocation function.
  void deallocate(void) { if (U != NULL) { delete []U; U = NULL; } }
  //@}

  //@{ @name Static variable functions.
  void set_particle_components(void);
  void set_particle_components(const int Nc);
  //@}

  //@{ @name Useful operators.
  //! Copy operator.
  void Copy(const Particle2D_cComponents &Uc);

  //! Vacuum operator.
  void Vacuum(void);

  //! Standard atmosphere operator.
  void Standard_Atmosphere(void);

  //! Check for unphysical state properties.
  int Unphysical_Properties(void) const;

  //! Reset unphysical state properties.
  void Reset_Unphysical_Properties(void);
  //@}

  //@{ @name Mass-averaged/total particle-phase quantities.
  //! Particle-phase concentration (total over the components).
  double sigma(void) const;

  //! Particle-phase number density (mass-averaged over the components).
  double nd(const double &mp) const;

  //! Particle-phase velocity (mass-averaged over the components).
  Vector2D u(void) const;

  //! Particle-phase momentum (mass-averaged over the components).
  Vector2D du(void) const;

  //! Particle-phase temperature (mass-averaged over the components).
  double Tp(const double &cm) const;

  //! Particle-phase total thermal/internal energy (mass-averaged over the components).
  double ep(void) const;

  //! Particle-phase total energy (mass-averaged over the components).
  double Ep(void) const;

  //! Particle-phase current (mass-averaged over the components).
  Vector2D jc(const double &qe, const double &mp) const;
  //@}

  //@{ @name Primitive solution state.
  Particle2D_pComponents W(const double &cm) const;
  Particle2D_pComponents W(const Particle2D_cComponents &U, const double &cm) const;
  //@}

  //! Jacobian of the conserved solution variables with respect to the
  //! primitive solution variables.
  void dUdW(DenseMatrix &dUdW, const int &n, const double &cm) const;

  //! Jacobian of the primitive solution variables with respect to the
  //! conserved solution variables.
  void dWdU(DenseMatrix &dWdU, const int &n, const double &cm) const;

  //@{ @name Solution flux (x-direction) and Jacobian.
  Particle2D_cComponents F(void) const;
  Particle2D_cComponents F(const Vector2D &V) const;
  void dFdU(DenseMatrix &dFdU, const int &n) const;
  void dFdU(DenseMatrix &dFdU, const int &n, const Vector2D &V) const;
  //@}

  //@{ @name Eigenstructure.
  //! Eigenvalue(s) (x-direction).
  Particle2D_pComponents lambda_x(void) const;

  //! Eigenvalue(s) for a moving reference frame (x-direction).
  Particle2D_pComponents lambda_x(const Vector2D &V) const;

  //! Primitive right eigenvector (x-direction).
  Particle2D_pComponents rp_x(int index) const;

  //! Conserved right eigenvector (x-direction).
  Particle2D_cComponents rc_x(int index) const;

  //! Primitive left eigenvector (x-direction).
  //Particle2D_pState lp_x(int index) const;
  //@}

  //@{ @name Index operator.
  Particle2D_cState &operator[](int index) {
    //assert(index >= 0 && index < NUM_CMP_PARTICLE2D);
    return U[index];
  }
    
  const Particle2D_cState &operator[](int index) const {
    //assert(index >= 0 && index < NUM_CMP_PARTICLE2D);
    return U[index];
  }

  void get_indices(int index, int &cindex, int &pindex) {
    if (NUM_CMP_PARTICLE2D == PARTICLE2D_SINGLE_VELOCITY_FORMULATION) cindex = 0;
    else cindex = (index%(NUM_CMP_PARTICLE2D*NUM_VAR_PARTICLE2D))/NUM_CMP_PARTICLE2D;
    pindex = index%NUM_VAR_PARTICLE2D + 1;
  }

  void get_indices(int index, int &cindex, int &pindex) const {
    if (NUM_CMP_PARTICLE2D == PARTICLE2D_SINGLE_VELOCITY_FORMULATION) cindex = 0;
    else cindex = (index%(NUM_CMP_PARTICLE2D*NUM_VAR_PARTICLE2D))/NUM_CMP_PARTICLE2D;
    pindex = index%NUM_VAR_PARTICLE2D + 1;
  }
  //@}

  //@{ @name Binary arithmetic operators.
  Particle2D_cComponents operator +(const Particle2D_cComponents &Uc) const;
  Particle2D_cComponents operator -(const Particle2D_cComponents &Uc) const;
  double operator *(const Particle2D_cComponents &Uc) const;
  Particle2D_cComponents operator *(const double &a) const;
  friend Particle2D_cComponents operator *(const double &a, const Particle2D_cComponents &Uc);
  Particle2D_cComponents operator /(const double &a) const;
  Particle2D_cComponents operator ^(const Particle2D_cComponents &Uc) const;
  //@}

  //@{ @name Assignment operator.
  Particle2D_cComponents& operator =(const Particle2D_cComponents &Uc);
  //@}

  //@{ @name Unary arithmetic operators.
  //Particle2D_cComponents operator +(const Particle2D_cComponents &Uc);
  friend Particle2D_cComponents operator -(const Particle2D_cComponents &Uc);
  //@}

  //@{ @name Shortcut arithmetic operators.
  Particle2D_cComponents &operator +=(const Particle2D_cComponents &Uc);
  Particle2D_cComponents &operator -=(const Particle2D_cComponents &Uc);
  Particle2D_cComponents &operator *=(const double &a);
  Particle2D_cComponents &operator /=(const double &a);
  //@}

  //@{ @name Relational operators.
  friend int operator ==(const Particle2D_cComponents &U1, const Particle2D_cComponents &U2);
  friend int operator !=(const Particle2D_cComponents &U1, const Particle2D_cComponents &U2);
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const Particle2D_cComponents &Uc);
  friend istream &operator >> (istream &in_file, Particle2D_cComponents &Uc);
  //@}

};

/**********************************************************************
 * Particle2D_pComponents::Copy -- Copy operator.                     *
 **********************************************************************/
inline void Particle2D_pComponents::Copy(const Particle2D_pComponents &Wc) {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) W[nc].Copy(Wc[nc]);
}

/**********************************************************************
 * Particle2D_pComponents::Vacuum -- Vacuum operator.                 *
 **********************************************************************/
inline void Particle2D_pComponents::Vacuum(void) {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) W[nc].Vacuum();
}

/**********************************************************************
 * Particle2D_pComponents::Standard_Atmosphere -- Standard atmosphere *
 *                                                operator.           *
 **********************************************************************/
inline void Particle2D_pComponents::Standard_Atmosphere(void) {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) W[nc].Vacuum();
}

/**********************************************************************
 * Particle2D_pComponents::Unphysical_Properties -- Check for         *
 *                                                  unphysical state  *
 *                                                  properties.       *
 **********************************************************************/
inline int Particle2D_pComponents::Unphysical_Properties(void) const {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) if (W[nc].Unphysical_Properties()) return 1;
  return 0;
}

/**********************************************************************
 * Particle2D_pComponents::Reset_Unphysical_Properties --             *
 *                                 Reset unphysical state properties. *
 **********************************************************************/
inline void Particle2D_pComponents::Reset_Unphysical_Properties(void) {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) W[nc].Reset_Unphysical_Properties();
}

/**********************************************************************
 * Particle2D_pComponents::set_particle_components                    *
 **********************************************************************/
inline void Particle2D_pComponents::set_particle_components(void) {
  // Deallocate component list.
  deallocate();
  // Set the number of particle components.
  NUM_CMP_PARTICLE2D = PARTICLE2D_SINGLE_VELOCITY_FORMULATION;
  // Allocate component list.
  allocate();
}

inline void Particle2D_pComponents::set_particle_components(const int Nc) {
  //assert(Nc == PARTICLE2D_SINGLE_VELOCITY_FORMULATION ||
  // Nc == PARTICLE2D_MULTI_VELOCITY_FORMULATION);
  // Deallocate component list.
  deallocate();
  // Set the number of particle components.
  NUM_CMP_PARTICLE2D = Nc;
  // Allocate component list.
  allocate();
}

/**********************************************************************
 * Particle2D_pComponents::sigma -- Particle-phase concentration      *
 *                                  (total over the components).      *
 **********************************************************************/
inline double Particle2D_pComponents::sigma(void) const {
  double c = ZERO;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) c += W[nc].sigma;
  return c;
}

/**********************************************************************
 * Particle2D_pComponents::nd -- Particle-phase number density (total *
 *                               over the components).                *
 **********************************************************************/
inline double Particle2D_pComponents::nd(const double &mp) const {
  double n = ZERO, c = sigma();
  if (c < TOLER*TOLER) return ZERO;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) n += W[nc].nd(mp);
  return n/c;
}

/**********************************************************************
 * Particle2D_pComponents::u -- Particle-phase velocity (mass-        *
 *                              averaged over the components).        *
 **********************************************************************/
inline Vector2D Particle2D_pComponents::u(void) const {
  Vector2D v = Vector2D_ZERO;
  double c = sigma();
  if (c < TOLER*TOLER) return Vector2D_ZERO;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) v += W[nc].sigma*W[nc].u;
  return v/c;
}

/**********************************************************************
 * Particle2D_pComponents::du -- Particle-phase momentum (mass-       *
 *                               averaged over the components).       *
 **********************************************************************/
inline Vector2D Particle2D_pComponents::du(void) const {
  Vector2D dv = Vector2D_ZERO;
  double c = sigma();
  if (c < TOLER*TOLER) return Vector2D_ZERO;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) dv += W[nc].sigma*W[nc].du();
  return dv/c;
}

/**********************************************************************
 * Particle2D_pComponents::Tp -- Particle-phase temperature (mass-    *
 *                               averaged over the components).       *
 **********************************************************************/
inline double Particle2D_pComponents::Tp(void) const {
  double T = ZERO, c = sigma();
  if (c < TOLER*TOLER) return ZERO;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) T += W[nc].sigma*W[nc].Tp;
  return T/c;
}

/**********************************************************************
 * Particle2D_pComponents::ep -- Particle-phase total thermal/        *
 *                               internal energy (mass-averaged over  *
 *                               the components.                      *
 **********************************************************************/
inline double Particle2D_pComponents::ep(const double &cm) const {
  double e = ZERO, c = sigma();
  if (c < TOLER*TOLER) return ZERO;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) e += W[nc].sigma*W[nc].ep(cm);
  return e/c;
}

/**********************************************************************
 * Particle2D_pComponents::Ep -- Particle-phase total energy (mass-   *
 *                               averaged over the components).       *
 **********************************************************************/
inline double Particle2D_pComponents::Ep(const double &cm) const {
  double E = ZERO, c = sigma();
  if (c < TOLER*TOLER) return ZERO;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) E += W[nc].sigma*W[nc].Ep(cm);
  return E/c;
}

/**********************************************************************
 * Particle2D_pComponents::jc -- Particle-phase current (mass-        *
 *                               averaged over the components).       *
 **********************************************************************/
inline Vector2D Particle2D_pComponents::jc(const double &qe, const double &mp) const {
  Vector2D j = Vector2D_ZERO;
  double c = sigma();
  if (c < TOLER*TOLER) return Vector2D_ZERO;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) j += W[nc].sigma*W[nc].jc(qe,mp);
  return j/c;
}

/**********************************************************************
 * Particle2D_pComponents -- Binary arithmetic operators.             *
 **********************************************************************/
 inline Particle2D_pComponents Particle2D_pComponents::operator +(const Particle2D_pComponents &Wc) const {
  Particle2D_pComponents Wtemp;
  Wtemp = *this;
  Wtemp += Wc;
  return Wtemp;
}

inline Particle2D_pComponents Particle2D_pComponents::operator -(const Particle2D_pComponents &Wc) const {
  Particle2D_pComponents Wtemp;
  Wtemp = *this;
  Wtemp -= Wc;
  return Wtemp;
}

// Inner product operator.
inline double Particle2D_pComponents::operator *(const Particle2D_pComponents &Wc) const {
  double sum;
  sum = ZERO;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) sum += W[nc]*Wc[nc];
  return sum;
}

inline Particle2D_pComponents Particle2D_pComponents::operator *(const double &a) const {
  Particle2D_pComponents Wtemp;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) Wtemp[nc] = W[nc]*a;
  return Wtemp;
}

inline Particle2D_pComponents operator *(const double &a, const Particle2D_pComponents &Wc) {
  Particle2D_pComponents Wtemp;
  for (int nc = 0; nc < Wc.NUM_CMP_PARTICLE2D; nc++) Wtemp[nc] = Wc[nc]*a;
  return Wtemp;
}

inline Particle2D_pComponents Particle2D_pComponents::operator /(const double &a) const {
  Particle2D_pComponents Wtemp;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) Wtemp.W[nc] = W[nc]/a;
  return Wtemp;
}

// A useful solution state product operator.
inline Particle2D_pComponents Particle2D_pComponents::operator ^(const Particle2D_pComponents &Wc) const {
  Particle2D_pComponents Wtemp;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) Wtemp[nc] = W[nc]^Wc[nc];  
  return Wtemp;
}

/**********************************************************************
 * Particle2D_pComponents -- Assignment operator.                     *
 **********************************************************************/
inline Particle2D_pComponents& Particle2D_pComponents::operator =(const Particle2D_pComponents &Wc) {
  //if (this != &Wc) {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) W[nc] = Wc[nc];
  //}
  return *this;
}

/**********************************************************************
 * Particle2D_pComponents -- Unary arithmetic operators.              *
 **********************************************************************/
//inline Particle2D_pComponents operator +(const Particle2D_pComponents &Wc) {
//return Wc;
//}

inline Particle2D_pComponents operator -(const Particle2D_pComponents &Wc) {
  Particle2D_pComponents Wtemp;
  for (int nc = 0; nc < Wc.NUM_CMP_PARTICLE2D; nc++) Wtemp[nc] = -Wc[nc];  
  return Wtemp;
}

/**********************************************************************
 * Particle2D_pComponents -- Shortcut arithmetic operators.           *
 **********************************************************************/
inline Particle2D_pComponents& Particle2D_pComponents::operator +=(const Particle2D_pComponents &Wc) {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) W[nc] += Wc[nc];
  return *this;
}

inline Particle2D_pComponents& Particle2D_pComponents::operator -=(const Particle2D_pComponents &Wc) {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) W[nc] -= Wc[nc];
  return *this;
}

inline Particle2D_pComponents& Particle2D_pComponents::operator *=(const double &a) {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) W[nc] *= a;
  return *this;
}

inline Particle2D_pComponents& Particle2D_pComponents::operator /=(const double &a) {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) W[nc] /= a;
  return *this;
}

/**********************************************************************
 * Particle2D_pComponents -- Relational operators.                    *
 **********************************************************************/
inline int operator ==(const Particle2D_pComponents &W1, const Particle2D_pComponents &W2) {
  int check = 1;
  for (int nc = 0; nc < W1.NUM_CMP_PARTICLE2D; nc++) if (W1[nc] != W2[nc]) check = 0;
  return check;
}

inline int operator !=(const Particle2D_pComponents &W1, const Particle2D_pComponents &W2) {
  int check = 0;
  for (int nc = 0; nc < W1.NUM_CMP_PARTICLE2D; nc++) if (W1[nc] == W2[nc]) check = 1;
  return check;
}

/**********************************************************************
 * Particle2D_pComponents -- Input-output operators.                  *
 **********************************************************************/
inline ostream &operator << (ostream &out_file, const Particle2D_pComponents &Wc) {
  for (int nc = 0; nc < Wc.NUM_CMP_PARTICLE2D; nc++) out_file << Wc[nc];
  return out_file;
}

inline istream &operator >> (istream &in_file, Particle2D_pComponents &Wc) {
  for (int nc = 0; nc < Wc.NUM_CMP_PARTICLE2D; nc++) in_file >> Wc[nc];
  return in_file;
}

/**********************************************************************
 * Particle2D_cComponents::Copy -- Copy operator.                     *
 **********************************************************************/
inline void Particle2D_cComponents::Copy(const Particle2D_cComponents &Uc) {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) U[nc].Copy(Uc[nc]);
}

/**********************************************************************
 * Particle2D_cComponents::Vacuum -- Vacuum operator.                 *
 **********************************************************************/
inline void Particle2D_cComponents::Vacuum(void) {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) U[nc].Vacuum();
}

/**********************************************************************
 * Particle2D_cComponents::Standard_Atmosphere -- Standard atmosphere *
 *                                                operator.           *
 **********************************************************************/
inline void Particle2D_cComponents::Standard_Atmosphere(void) {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) U[nc].Vacuum();
}

/**********************************************************************
 * Particle2D_cComponents::Unphysical_Properties -- Check for         *
 *                                                  unphysical state  *
 *                                                  properties.       *
 **********************************************************************/
inline int Particle2D_cComponents::Unphysical_Properties(void) const {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) if (U[nc].Unphysical_Properties()) return 1;
  return 0;
}

/**********************************************************************
 * Particle2D_cComponents::Reset_Unphysical_Properties --             *
 *                                 Reset unphysical state properties. *
 **********************************************************************/
inline void Particle2D_cComponents::Reset_Unphysical_Properties(void) {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) U[nc].Reset_Unphysical_Properties();
}

/**********************************************************************
 * Particle2D_cComponents::set_particle_components                    *
 **********************************************************************/
inline void Particle2D_cComponents::set_particle_components(void) {
  // Deallocate component list.
  deallocate();
  // Set the number of particle components.
  NUM_CMP_PARTICLE2D = 1;
  // Allocate component list.
  allocate();
}

inline void Particle2D_cComponents::set_particle_components(const int Nc) {
  //assert(Nc == PARTICLE2D_SINGLE_VELOCITY_FORMULATION ||
  // Nc == PARTICLE2D_MULTI_VELOCITY_FORMULATION);
  // Deallocate component list.
  deallocate();
  // Set the number of particle components.
  NUM_CMP_PARTICLE2D = Nc;
  // Allocate component list.
  allocate();
}

/**********************************************************************
 * Particle2D_cComponents::sigma -- Particle-phase concentration      *
 *                                  (total over the components).      *
 **********************************************************************/
inline double Particle2D_cComponents::sigma(void) const {
  double c = ZERO;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) c += U[nc].sigma;
  return c;
}

/**********************************************************************
 * Particle2D_cComponents::nd -- Particle-phase number density (total *
 *                               over the components).                *
 **********************************************************************/
inline double Particle2D_cComponents::nd(const double &mp) const {
  double n = ZERO, c = sigma();
  if (c < TOLER*TOLER) return ZERO;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) n += U[nc].nd(mp);
  return n/c;
}

/**********************************************************************
 * Particle2D_cComponents::u -- Particle-phase velocity (mass-        *
 *                              averaged over the components).        *
 **********************************************************************/
inline Vector2D Particle2D_cComponents::u(void) const {
  Vector2D v = Vector2D_ZERO;
  double c = sigma();
  if (c < TOLER*TOLER) return Vector2D_ZERO;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) v += U[nc].sigma*U[nc].u();
  return v/c;
}

/**********************************************************************
 * Particle2D_cComponents::du -- Particle-phase momentum (mass-       *
 *                               averaged over the components).       *
 **********************************************************************/
inline Vector2D Particle2D_cComponents::du(void) const {
  Vector2D dv = Vector2D_ZERO;
  double c = sigma();
  if (c < TOLER*TOLER) return Vector2D_ZERO;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) dv += U[nc].sigma*U[nc].du;
  return dv/c;
}

/**********************************************************************
 * Particle2D_cComponents::Tp -- Particle-phase temperature (mass-    *
 *                               averaged over the components).       *
 **********************************************************************/
inline double Particle2D_cComponents::Tp(const double &cm) const {
  double T = ZERO, c = sigma();
  if (c < TOLER*TOLER) return ZERO;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) T += U[nc].sigma*U[nc].Tp(cm);
  return T/c;
}

/**********************************************************************
 * Particle2D_cComponents::ep -- Particle-phase specific internal     *
 *                               energy (mass-averaged over the       *
 *                               components.                          *
 **********************************************************************/
inline double Particle2D_cComponents::ep(void) const {
  double e = ZERO, c = sigma();
  if (c < TOLER*TOLER) return ZERO;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) e += U[nc].sigma*U[nc].ep;
  return e/c;
}

/**********************************************************************
 * Particle2D_cComponents::Ep -- Particle-phase total energy (mass-   *
 *                               averaged over the components).       *
 **********************************************************************/
inline double Particle2D_cComponents::Ep(void) const {
  double E = ZERO, c = sigma();
  if (c < TOLER*TOLER) return ZERO;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) E += U[nc].sigma*U[nc].Ep();
  return E/c;
}

/**********************************************************************
 * Particle2D_cComponents::jc -- Particle-phase current (mass-        *
 *                               averaged over the components).       *
 **********************************************************************/
inline Vector2D Particle2D_cComponents::jc(const double &qe, const double &mp) const {
  Vector2D j = Vector2D_ZERO;
  double c = sigma();
  if (c < TOLER*TOLER) return Vector2D_ZERO;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) j += U[nc].sigma*U[nc].jc(qe,mp);
  return j/c;
}

/**********************************************************************
 * Particle2D_cComponents -- Binary arithmetic operators.             *
 **********************************************************************/
inline Particle2D_cComponents Particle2D_cComponents::operator +(const Particle2D_cComponents &Uc) const {
  Particle2D_cComponents Utemp;
  Utemp = *this;
  Utemp += Uc;
  return Utemp;
}

inline Particle2D_cComponents Particle2D_cComponents::operator -(const Particle2D_cComponents &Uc) const {
  Particle2D_cComponents Utemp;
  Utemp = *this;
  Utemp -= Uc;
  return Utemp;
}

// Inner product operator.
inline double Particle2D_cComponents::operator *(const Particle2D_cComponents &Uc) const {
  double sum;
  sum = ZERO;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) sum += U[nc]*Uc[nc];
  return sum;
}

inline Particle2D_cComponents Particle2D_cComponents::operator *(const double &a) const {
  Particle2D_cComponents Utemp;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) Utemp[nc] = U[nc]*a;
  return Utemp;
}

inline Particle2D_cComponents operator *(const double &a, const Particle2D_cComponents &Uc) {
  Particle2D_cComponents Utemp;
  for (int nc = 0; nc < Uc.NUM_CMP_PARTICLE2D; nc++) Utemp[nc] = Uc[nc]*a;
  return Utemp;
}

inline Particle2D_cComponents Particle2D_cComponents::operator /(const double &a) const {
  Particle2D_cComponents Utemp;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) Utemp[nc] = U[nc]/a;
  return Utemp;
}

// A useful solution state product operator.
inline Particle2D_cComponents Particle2D_cComponents::operator ^(const Particle2D_cComponents &Uc) const {
  Particle2D_cComponents Utemp;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) Utemp[nc] = U[nc]^Uc[nc];  
  return Utemp;
}

/**********************************************************************
 * Particle2D_cComponents -- Assignment operator.                     *
 **********************************************************************/
inline Particle2D_cComponents& Particle2D_cComponents::operator =(const Particle2D_cComponents &Uc) {
  //if (this != &Uc) {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) U[nc] = Uc[nc];
  //}
  return *this;
}

/**********************************************************************
 * Particle2D_cComponents -- Unary arithmetic operators.              *
 **********************************************************************/
//inline Particle2D_cComponents operator +(const Particle2D_cComponents &Uc) {
//return Uc;
//}

inline Particle2D_cComponents operator -(const Particle2D_cComponents &Uc) {
  Particle2D_cComponents Utemp;
  for (int nc = 0; nc < Uc.NUM_CMP_PARTICLE2D; nc++) Utemp[nc] = -Uc[nc];
  return Utemp;
}

/**********************************************************************
 * Particle2D_cComponents -- Shortcut arithmetic operators.           *
 **********************************************************************/
inline Particle2D_cComponents& Particle2D_cComponents::operator +=(const Particle2D_cComponents &Uc) {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) U[nc] += Uc[nc];
  return *this;
}

inline Particle2D_cComponents& Particle2D_cComponents::operator -=(const Particle2D_cComponents &Uc) {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) U[nc] -= Uc[nc];
  return *this;
}

inline Particle2D_cComponents& Particle2D_cComponents::operator *=(const double &a) {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) U[nc] *= a;
  return *this;
}

inline Particle2D_cComponents& Particle2D_cComponents::operator /=(const double &a) {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) U[nc] /= a;
  return *this;
}

/**********************************************************************
 * Particle2D_cComponents -- Relational operators.                    *
 **********************************************************************/
inline int operator ==(const Particle2D_cComponents &U1, const Particle2D_cComponents &U2) {
  int check = 1;
  for (int nc = 0; nc < U1.NUM_CMP_PARTICLE2D; nc++) if (U1[nc] != U2[nc]) check = 0;
  return check;
}

inline int operator !=(const Particle2D_cComponents &U1, const Particle2D_cComponents &U2) {
  int check = 0;
  for (int nc = 0; nc < U1.NUM_CMP_PARTICLE2D; nc++) if (U1[nc] == U2[nc]) check = 1;
  return check;
}

/**********************************************************************
 * Particle2D_cComponents -- Input-output operators.                  *
 **********************************************************************/
inline ostream &operator << (ostream &out_file, const Particle2D_cComponents &Uc) {
  for (int nc = 0; nc < Uc.NUM_CMP_PARTICLE2D; nc++) out_file << Uc[nc];
  return out_file;
}

inline istream &operator >> (istream &in_file, Particle2D_cComponents &Uc) {
  for (int nc = 0; nc < Uc.NUM_CMP_PARTICLE2D; nc++) in_file >> Uc[nc];
  return in_file;
}

/**********************************************************************
 * Particle2D_pComponents::Particle2D_pComponents -- Constructors.    *
 **********************************************************************/
inline Particle2D_pComponents::Particle2D_pComponents(const Particle2D_pComponents &Wc) {
  W = NULL; allocate();
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) W[nc] = Wc[nc];
}

inline Particle2D_pComponents::Particle2D_pComponents(const Particle2D_cComponents &Uc,
						      const double &cm) {
  W = NULL; allocate();
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) W[nc] = Uc[nc].W(cm);
}

inline Particle2D_pComponents::Particle2D_pComponents(const Particle2D_pState &Wp) {
  W = NULL; allocate();
  if (NUM_CMP_PARTICLE2D == PARTICLE2D_SINGLE_VELOCITY_FORMULATION) {
    // Single-velocity formulation.
    W[0] = Wp;
  } else {
    // Multi-velocity formulation.
    //if (dot(Wp.u,ihat) > ZERO && dot(Wp.u,jhat) >= ZERO) {
    if (dot(Wp.u,ihat) >= ZERO && dot(Wp.u,jhat) >= ZERO) {
      W[0] = Wp;
    } else if (dot(Wp.u,ihat) <= ZERO && dot(Wp.u,jhat) > ZERO) {
      W[1] = Wp;
    } else if (dot(Wp.u,ihat) < ZERO && dot(Wp.u,jhat) <= ZERO) {
      W[2] = Wp;
    } else if (dot(Wp.u,ihat) >= ZERO && dot(Wp.u,jhat) < ZERO) {
      W[3] = Wp;
    }
  }
}

inline Particle2D_pComponents::Particle2D_pComponents(const Particle2D_cState &Up,
						      const double &cm) {
  W = NULL; allocate();
  Particle2D_pState Wp;
  Wp = Up.W(cm);
  if (NUM_CMP_PARTICLE2D == PARTICLE2D_SINGLE_VELOCITY_FORMULATION) {
    // Single-velocity formulation.
    W[0] = Wp;
  } else {
    // Multi-velocity formulation.
    //if (dot(Wp.u,ihat) > ZERO && dot(Wp.u,jhat) >= ZERO) {
    if (dot(Wp.u,ihat) >= ZERO && dot(Wp.u,jhat) >= ZERO) {
      W[0] = Wp;
    } else if (dot(Wp.u,ihat) <= ZERO && dot(Wp.u,jhat) > ZERO) {
      W[1] = Wp;
    } else if (dot(Wp.u,ihat) < ZERO && dot(Wp.u,jhat) <= ZERO) {
      W[2] = Wp;
    } else if (dot(Wp.u,ihat) >= ZERO && dot(Wp.u,jhat) < ZERO) {
      W[3] = Wp;
    }
  }
}

/**********************************************************************
 * Particle2D_pComponents::U -- Conserved solution state.             *
 **********************************************************************/
inline Particle2D_cComponents Particle2D_pComponents::U(const double &cm) const {
  return U(*this,cm);
}

inline Particle2D_cComponents Particle2D_pComponents::U(const Particle2D_pComponents &Wc,
							const double &cm) const {
  Particle2D_cComponents Utemp;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) Utemp[nc] = Wc[nc].U(cm);
  return Utemp;
}

/**********************************************************************
 * Particle2D_pComponents::dUdW -- Jacobian of the conserved solution *
 *                                 variables with respect to the      *
 *                                 primitive solution variables.      *
 **********************************************************************/
inline void Particle2D_pComponents::dUdW(DenseMatrix &dUdW, const int &n, const double &cm) const {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++)
    W[nc].dUdW(dUdW,NUM_VAR_PARTICLE2D*nc+n,cm);
}

/**********************************************************************
 * Particle2D_pComponents::dWdU -- Jacobian of the primitive solution *
 *                                 variables with respect to the      *
 *                                 conserved solution variables.      *
 **********************************************************************/
inline void Particle2D_pComponents::dWdU(DenseMatrix &dWdU, const int &n, const double &cm) const {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++)
    W[nc].dWdU(dWdU,NUM_VAR_PARTICLE2D*nc+n,cm);
}

/**********************************************************************
 * Particle2D_pComponents::F -- Solution flux (x-direction).          *
 **********************************************************************/
inline Particle2D_cComponents Particle2D_pComponents::F(const double &cm) const {
  Particle2D_cComponents Flux;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) Flux[nc] = W[nc].F(cm);
  return Flux;
}

inline Particle2D_cComponents Particle2D_pComponents::F(const Vector2D &V,
							const double &cm) const {
  Particle2D_cComponents Flux;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) Flux[nc] = W[nc].F(V,cm);
  return Flux;
}

/**********************************************************************
 * Particle2D_pComponents::dFdU -- Jacobian of the solution flux with *
 *                                 respect to the conserved solution  *
 *                                 variables.                         *
 **********************************************************************/
inline void Particle2D_pComponents::dFdU(DenseMatrix &dFdU, const int &n, const double &cm) const {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++)
    W[nc].dFdU(dFdU,NUM_VAR_PARTICLE2D*nc+n,cm);
}

inline void Particle2D_pComponents::dFdU(DenseMatrix &dFdU, const int &n,
					 const Vector2D &V, const double &cm) const {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++)
    W[nc].dFdU(dFdU,NUM_VAR_PARTICLE2D*nc+n,V,cm);
}

/**********************************************************************
 * Particle2D_pComponents::lambda_x -- Eigenvalue(s) (x-direction).   *
 **********************************************************************/
inline Particle2D_pComponents Particle2D_pComponents::lambda_x(void) const {
  Particle2D_pComponents lambda;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) lambda = W[nc].lambda_x();
  return lambda;
}

inline Particle2D_pComponents Particle2D_pComponents::lambda_x(const Vector2D &V) const {
  Particle2D_pComponents lambda;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) lambda = W[nc].lambda_x(V);
  return lambda;
}

/**********************************************************************
 * Particle2D_pComponents::rp_x -- Primitive right eigenvectors       *
 *                                 (x-direction).                     *
 **********************************************************************/
inline Particle2D_pComponents Particle2D_pComponents::rp_x(int index) const {
  Particle2D_pComponents rpx;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) rpx[nc] = W[nc].rp_x(index);
  return rpx;
}

/**********************************************************************
 * Particle2D_pComponents::rc_x -- Conserved right eigenvectors       *
 *                                 (x-direction).                     *
 **********************************************************************/
inline Particle2D_cComponents Particle2D_pComponents::rc_x(int index) const {
  Particle2D_cComponents rcx;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) rcx[nc] = W[nc].rc_x(index);
  return rcx;
}

/**********************************************************************
 * Particle2D_pComponents::lp_x -- Primitive left eigenvector         *
 *                                 (x-direction).                     *
 **********************************************************************/
// inline Particle2D_pComponents Particle2D_pComponents::lp_x(int index) const {
//   Particle2D_pComponents lpx;
//   for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) lpx[nc] = W[nc].lp_x(index);
//   return lpx;
// }

/**********************************************************************
 * Particle2D_pComponents::Sa -- Axisymmetric source term vector and  *
 *                               Jacobian.                            *
 **********************************************************************/
inline Particle2D_cComponents Particle2D_pComponents::Sa(const Vector2D &X,
							 const double &cm) const {
  Particle2D_cComponents S;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) S[nc] = W[nc].Sa(X,cm);
  return S;
}

inline void Particle2D_pComponents::dSadU(DenseMatrix &dSadU, const Vector2D &X,
					  const int &n, const double &cm) const {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++)
    W[nc].dSadU(dSadU,X,NUM_VAR_PARTICLE2D*nc+n,cm);
}

/**********************************************************************
 * Particle2D_pComponents::Se -- Electrostatic-force source term      *
 *                               vector and Jacobian.                 *
 **********************************************************************/
inline Particle2D_cComponents Particle2D_pComponents::Se(const Vector2D &E,
							 const double &Cm) const {
  Particle2D_cComponents S;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) S[nc] = W[nc].Se(E,Cm);
  return S;
}

inline void Particle2D_pComponents::dSedU(DenseMatrix &dSedU, const int &n,
					  const Vector2D &E, const double &Cm) const {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++)
    W[nc].dSedU(dSedU,NUM_VAR_PARTICLE2D*nc+n,E,Cm);
}

/**********************************************************************
 * Particle2D_cComponents::Particle2D_cComponents -- Constructor.     *
 **********************************************************************/
inline Particle2D_cComponents::Particle2D_cComponents(const Particle2D_cComponents &Uc) {
  U = NULL; allocate();
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) U[nc] = Uc[nc];
}

inline Particle2D_cComponents::Particle2D_cComponents(const Particle2D_pComponents &Wc,
						      const double &cm) {
  U = NULL; allocate();
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) U[nc] = Wc[nc].U(cm);
}

inline Particle2D_cComponents::Particle2D_cComponents(const Particle2D_cState &Up) {
  U = NULL; allocate();
  if (NUM_CMP_PARTICLE2D == PARTICLE2D_SINGLE_VELOCITY_FORMULATION) {
    // Single-velocity formulation.
    U[0] = Up;
  } else {
    // Multi-velocity formulation.
    //if (dot(Up.du,ihat) > ZERO && dot(Up.du,jhat) >= ZERO) {
    if (dot(Up.du,ihat) >= ZERO && dot(Up.du,jhat) >= ZERO) {
      U[0] = Up;
    } else if (dot(Up.du,ihat) <= ZERO && dot(Up.du,jhat) > ZERO) {
      U[1] = Up;
    } else if (dot(Up.du,ihat) < ZERO && dot(Up.du,jhat) <= ZERO) {
      U[2] = Up;
    } else if (dot(Up.du,ihat) >= ZERO && dot(Up.du,jhat) < ZERO) {
      U[3] = Up;
    }
  }
}

inline Particle2D_cComponents::Particle2D_cComponents(const Particle2D_pState &Wp,
						      const double &cm) {
  U = NULL; allocate();
  Particle2D_cState Up;
  Up = Wp.U(cm);
  if (NUM_CMP_PARTICLE2D == PARTICLE2D_SINGLE_VELOCITY_FORMULATION) {
    // Single-velocity formulation.
    U[0] = Up;
  } else {
    // Multi-velocity formulation.
    //if (dot(Up.du,ihat) > ZERO && dot(Up.du,jhat) >= ZERO) {
    if (dot(Up.du,ihat) >= ZERO && dot(Up.du,jhat) >= ZERO) {
      U[0] = Up;
    } else if (dot(Up.du,ihat) <= ZERO && dot(Up.du,jhat) > ZERO) {
      U[1] = Up;
    } else if (dot(Up.du,ihat) < ZERO && dot(Up.du,jhat) <= ZERO) {
      U[2] = Up;
    } else if (dot(Up.du,ihat) >= ZERO && dot(Up.du,jhat) < ZERO) {
      U[3] = Up;
    }
  }
}

/**********************************************************************
 * Particle2D_cComponents::W -- Primitive solution state.             *
 **********************************************************************/
inline Particle2D_pComponents Particle2D_cComponents::W(const double &cm) const {
  return W(*this,cm);
}

inline Particle2D_pComponents Particle2D_cComponents::W(const Particle2D_cComponents &Uc,
							const double &cm) const {
  Particle2D_pComponents Wtemp;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) Wtemp[nc] = Uc[nc].W(cm);
  return Wtemp;
}

/**********************************************************************
 * Particle2D_cComponents::dUdW -- Jacobian of the conserved solution *
 *                                 variables with respect to the      *
 *                                 primitive solution variables.      *
 **********************************************************************/
inline void Particle2D_cComponents::dUdW(DenseMatrix &dUdW, const int &n, const double &cm) const {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++)
    U[nc].dUdW(dUdW,NUM_VAR_PARTICLE2D*nc+n,cm);
}

/**********************************************************************
 * Particle2D_cComponents::dWdU -- Jacobian of the primitive solution *
 *                                 variables with respect to the      *
 *                                 conserved solution variables.      *
 **********************************************************************/
inline void Particle2D_cComponents::dWdU(DenseMatrix &dWdU, const int &n, const double &cm) const {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++)
    U[nc].dWdU(dWdU,NUM_VAR_PARTICLE2D*nc+n,cm);
}

/**********************************************************************
 * Particle2D_cComponents::F -- Solution flux (x-direction).          *
 **********************************************************************/
inline Particle2D_cComponents Particle2D_cComponents::F(void) const {
  Particle2D_cComponents Flux;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) Flux[nc] = U[nc].F();
  return Flux;
}

inline Particle2D_cComponents Particle2D_cComponents::F(const Vector2D &V) const {
  Particle2D_cComponents Flux;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) Flux[nc] = U[nc].F(V);
  return Flux;
}

/**********************************************************************
 * Particle2D_pComponents::dFdU -- Jacobian of the solution flux with *
 *                                 respect to the conserved solution  *
 *                                 variables.                         *
 **********************************************************************/
inline void Particle2D_cComponents::dFdU(DenseMatrix &dFdU, const int &n) const {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++)
    U[nc].dFdU(dFdU,NUM_VAR_PARTICLE2D*nc+n);
}

inline void Particle2D_cComponents::dFdU(DenseMatrix &dFdU, const int &n,
					 const Vector2D &V) const {
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++)
    U[nc].dFdU(dFdU,NUM_VAR_PARTICLE2D*nc+n,V);
}

/**********************************************************************
 * Particle2D_cComponents::lambda_x -- Eigenvalue(s) (x-direction).   *
 **********************************************************************/
inline Particle2D_pComponents Particle2D_cComponents::lambda_x(void) const {
  Particle2D_pComponents lambda;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) lambda[nc] = U[nc].lambda_x();
  return lambda;
}

inline Particle2D_pComponents Particle2D_cComponents::lambda_x(const Vector2D &V) const {
  Particle2D_pComponents lambda;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) lambda[nc] = U[nc].lambda_x(V);
  return lambda;
}

/**********************************************************************
 * Particle2D_pComponents::rp_x -- Primitive right eigenvectors       *
 *                                 (x-direction).                     *
 **********************************************************************/
inline Particle2D_pComponents Particle2D_cComponents::rp_x(int index) const {
  Particle2D_pComponents rpx;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) rpx[nc] = U[nc].rp_x(index);
  return rpx;
}

/**********************************************************************
 * Particle2D_pComponents::rc_x -- Conserved right eigenvectors       *
 *                                 (x-direction).                     *
 **********************************************************************/
inline Particle2D_cComponents Particle2D_cComponents::rc_x(int index) const {
  Particle2D_cComponents rcx;
  for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) rcx[nc] = U[nc].rc_x(index);
  return rcx;
}

/**********************************************************************
 * Particle2D_pComponents::lp_x -- Primitive left eigenvector         *
 *                                 (x-direction).                     *
 **********************************************************************/
// inline Particle2D_pComponents Particle2D_cComponents::lp_x(int index) const {
//   Particle2D_pComponents lpx;
//   for (int nc = 0; nc < NUM_CMP_PARTICLE2D; nc++) lpx[nc] = U[nc].lp_x(index);
//   return lpx;
// }

/**********************************************************************
 * Particle2DComponents -- External subroutines.                      *
 **********************************************************************/

extern Particle2D_pComponents Rotate(const Particle2D_pComponents &W,
				     const Vector2D &norm_dir);

extern Particle2D_cComponents Rotate(const Particle2D_cComponents &U,
				     const Vector2D &norm_dir);

extern Particle2D_pComponents Translate(const Particle2D_pComponents &W,
					const Vector2D &V);

extern Particle2D_pComponents Reflect(const Particle2D_pComponents &W,
				      const Vector2D &norm_dir);

extern Particle2D_pComponents Mirror(const Particle2D_pComponents &W,
				     const Vector2D &norm_dir);

extern Particle2D_pComponents Absorb(const Particle2D_pComponents &W,
				     const Vector2D &norm_dir);

extern Particle2D_cComponents FluxSaurel_n(const Particle2D_pComponents &Wl,
					   const Particle2D_pComponents &Wr,
					   const Vector2D &norm_dir,
					   const double &cm);

extern Particle2D_cComponents FluxSaurel_MB_n(const Particle2D_pComponents &Wl,
					      const Particle2D_pComponents &Wr,
					      const Vector2D &V,
					      const Vector2D &norm_dir,
					      const double &cm);

extern Particle2D_pComponents FluxMultiVelocity(const Particle2D_pComponents &Wl,
						const Particle2D_pComponents &Wr);

extern Particle2D_pComponents FluxMultiVelocity(const Particle2D_cComponents &Ul,
						const Particle2D_cComponents &Ur,
						const double &cm);

extern Particle2D_cComponents FluxMultiVelocity_n(const Particle2D_pComponents &Wl,
						  const Particle2D_pComponents &Wr,
						  const Vector2D &norm_dir,
						  const double &cm);

extern Particle2D_cComponents FluxMultiVelocity_n(const Particle2D_cComponents &Ul,
						  const Particle2D_cComponents &Ur,
						  const Vector2D &norm_dir,
						  const double &cm);

extern Particle2D_cComponents FluxMultiVelocity_MB_n(const Particle2D_pComponents &Wl,
						     const Particle2D_pComponents &Wr,
						     const Vector2D &V,
						     const Vector2D &norm_dir,
						     const double &cm);

#endif // _PARTICLE2D_COMPONENTS_INCLUDED
