/* Levermore1DState.h:  Header file defining 1D Levermore Solution State Classes. */

#ifndef _LEVERMORE1D_STATE_INCLUDED
#define _LEVERMORE1D_STATE_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cmath>

using namespace std;

/* Include Headers */
#ifndef _LEVERMORE1D_VECTOR_INCLUDED
#include "Levermore1DVector.h"
#endif //_LEVERMORE1D_VECTOR_INCLUDED

#ifndef _MATRIX_INCLUDED
#include "../Math/Matrix.h"
#endif //_MATRIX_INCLUDED

#ifndef _GAS_CONSTANTS_INCLUDED
#include "../Physics/GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED

/* Define the classes. */
class Levermore1D_cState;
class Levermore1D_weights;

/********************************************************
 * Class: Levermore1D_pState                            *
 ********************************************************/
class Levermore1D_pState : public Levermore1D_Vector{

 public:

  Levermore1D_pState(void){}
  Levermore1D_pState(const Levermore1D_pState &W) : Levermore1D_Vector(W) {}
  explicit Levermore1D_pState(const Levermore1D_cState &U) {set_from_U(U);}
  explicit Levermore1D_pState(const Levermore1D_weights &A) {set_from_A(A);}

  /* Functions. */
  Levermore1D_Vector& operator=(const Levermore1D_Vector &V) {return Levermore1D_Vector::operator=(V);}
  void Vacuum() {Levermore1D_Vector::zero();}
  void set_from_U(const Levermore1D_cState &U);
  void set_from_A(const Levermore1D_weights &A);
  double conserved_extras(int i) const;

 protected:
  double conserved_extras_recursive(int i, int &pf, int pf_num, int pf_den) const;

};

/********************************************************
 * Class: Levermore1D_cState                            *
 ********************************************************/
class Levermore1D_cState : public Levermore1D_Vector{

  protected:
  public:

  Levermore1D_cState(void){}
  Levermore1D_cState(const Levermore1D_cState &U) : Levermore1D_Vector(U) {}
  explicit Levermore1D_cState(const Levermore1D_pState &W) {set_from_W(W);}
  explicit Levermore1D_cState(const Levermore1D_weights &A) {set_from_A(A);}

  /* Functions. */
  Levermore1D_Vector& operator=(const Levermore1D_Vector &V) {return Levermore1D_Vector::operator=(V);}
  void Vacuum() {Levermore1D_Vector::zero();}
  void set_from_W(const Levermore1D_pState &W);
  void set_from_A(const Levermore1D_weights &A);
  double moment(int n, const Levermore1D_weights &A) const;
  DenseMatrix d2hda2(const Levermore1D_weights &A) const;

};

/********************************************************
 * Class: Levermore1D_weights                           *
 ********************************************************/
class Levermore1D_weights : public Levermore1D_Vector{

  public:

  Levermore1D_weights(void){}
  Levermore1D_weights(const Levermore1D_weights &A) : Levermore1D_Vector(A) {}
  Levermore1D_weights(const double &d, const double &u, const double &p) {MaxBoltz(d,u,p);}
  explicit Levermore1D_weights(const Levermore1D_pState &W) {MaxBoltz(W); set_from_W(W);}
  explicit Levermore1D_weights(const Levermore1D_cState &U) {MaxBoltz(U); set_from_U(U);}

  /* Functions. */
  Levermore1D_Vector& operator=(const Levermore1D_Vector &V) {return Levermore1D_Vector::operator=(V);}
  void set_from_W(const Levermore1D_pState &W);
  void set_from_U(const Levermore1D_cState &U);
  double integrate_conserved_moment(int i) const;
  double integrate_random_moment(int i, double u) const;

  /* Inline Functions. */
  double value_at(double v) const {return exp(exponent_value_recursive(v,0));}
  double velocity_weighted_value_at(double v, int i) const {
    return pow(v,(double)i)*particle_mass*exp(exponent_value_recursive(v,0));
  }
  double random_velocity_weighted_value_at(double v, double u, int i) const {
    return pow(v-u,(double)i)*particle_mass*exp(exponent_value_recursive(v,0));
  }

  /* Static Functions. */
  static double m() {return particle_mass;}
  static void set_particle_mass(double m) {particle_mass=m;}

  protected:

  static double particle_mass;

  /* Inline Functions */
  double exponent_value_recursive(double v, int i) const {
    if(i<length-1) {
      return m_values[i] + v*exponent_value_recursive(v,i+1);
    } else {
      return m_values[i];
    }
  }

  void MaxBoltz(const Levermore1D_pState &W) {MaxBoltz(W[1], W[2], W[3]);}
  void MaxBoltz(const Levermore1D_cState &U) {MaxBoltz(U[1], U[2]/U[1], U[3]-U[2]*U[2]/U[1]);}
  void MaxBoltz(double rho, double u, double p) {
    double n(rho/m()); //number density
    double B(rho/(2.0*p));
    zero();
    m_values[0] = -B*u*u+log(n*sqrt(B/PI));
    m_values[1] = 2.0*B*u;
    m_values[2] = -B;
  }

};

/********************************************************
 *              External  Functions                     *
 ********************************************************/
inline ostream& operator<<(ostream &out, const Levermore1D_pState &W) {
  W.output(out);
  return out;
}
inline ostream& operator<<(ostream &out, const Levermore1D_cState &U) {
  U.output(out);
  return out;
}
inline ostream& operator<<(ostream &out, const Levermore1D_weights &A) {
  A.output(out);
  return out;
}
inline istream& operator<<(istream &in, const Levermore1D_pState &W) {
  W.input(in);
  return in;
}
inline istream& operator<<(istream &in, const Levermore1D_cState &U) {
  U.input(in);
  return in;
}
inline istream& operator<<(istream &in, const Levermore1D_weights &A) {
  A.input(in);
  return in;
}

#endif
