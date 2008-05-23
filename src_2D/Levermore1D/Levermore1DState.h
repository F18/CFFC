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

#ifndef _LINEARSYSTEMS_INCLUDED
#include "../Math/LinearSystems.h"
#endif //_LINEARSYSTEMS_INCLUDED

#ifndef _GAS_CONSTANTS_INCLUDED
#include "../Physics/GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED

#define STABILIZATION 1.0e-18  //positive numbers only please.
                               //otherwise it's anti-stabilization.

#define EXP_LIMIT                250.0

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
  explicit Levermore1D_pState(const Levermore1D_Vector &V) : Levermore1D_Vector(V) {}
  explicit Levermore1D_pState(const Levermore1D_cState &U) {set_from_U(U);}
  explicit Levermore1D_pState(const Levermore1D_weights &A, double us) {set_from_A(A,us);}
  Levermore1D_pState(double rho, double u, double p) { //Maxwell-Boltzmann
    int coef;
    zero();
    m_values[0] = rho; m_values[1] = u; m_values[2] = p;
    for(int i=4; i<Levermore1D_Vector::get_length(); i=i+2) {
      m_values[i] = (double)Double_Factorial(i-1) * pow(p,i/2) / pow(rho,i/2-1);
    }
  }

  /* Functions. */
  Levermore1D_Vector& operator=(const Levermore1D_Vector &V) {return Levermore1D_Vector::operator=(V);}
  void Vacuum() {Levermore1D_Vector::zero();}
  void set_from_U(const Levermore1D_cState &U);
  void set_from_A(const Levermore1D_weights &A, double us);
  static void set_relaxation_time(double tau) {m_relaxation_time = tau;}
  double conserved_extras(int i) const;
  double relaxation_time() const;
  DenseMatrix dUdW(void) const;
  DenseMatrix dU_MBdW(void) const;
  DenseMatrix dSdW(void) const;
  DenseMatrix dSdU(void) const;
  int valid(void) const {
    //this is hard-coded for the length == 5 case....must be changed later!!!
    return (m_values[0] > 0 && m_values[2] > 0
	    && m_values[4] > m_values[2]*m_values[2]/m_values[0]+m_values[3]*m_values[3]/m_values[2]);
  }

 protected:
  double conserved_extras_recursive(int i, int &pf, int pf_num, int pf_den) const;
  static double m_relaxation_time;
};

/********************************************************
 * Class: Levermore1D_cState                            *
 ********************************************************/
class Levermore1D_cState : public Levermore1D_Vector{

  protected:
  public:

  Levermore1D_cState(void){}
  Levermore1D_cState(const Levermore1D_cState &U) : Levermore1D_Vector(U) {}
  explicit Levermore1D_cState(const Levermore1D_Vector &V) : Levermore1D_Vector(V) {}
  explicit Levermore1D_cState(const Levermore1D_pState &W) {set_from_W(W);}
  explicit Levermore1D_cState(const Levermore1D_weights &A, double us) {set_from_A(A,us);}
  Levermore1D_cState(double rho, double u, double p) { //Maxwell-Boltzmann
    MaxBoltz(rho,u,p);
  }

  /* Functions. */
  Levermore1D_Vector& operator=(const Levermore1D_Vector &V) {return Levermore1D_Vector::operator=(V);}
  void Vacuum() {Levermore1D_Vector::zero();}
  void set_from_W(const Levermore1D_pState &W);
  void set_from_A(const Levermore1D_weights &A, double us);
  static void set_relaxation_time(double tau) {m_relaxation_time = tau;}
  void MaxBoltz(double rho, double u, double p) {
    set_from_W(Levermore1D_pState(rho,u,p)); // change to something better later.
  }
  void MaxBoltz() {
    MaxBoltz(rho(),u(),p());
  }

  double rho() const {return m_values[0];}
  double u() const {return m_values[1]/m_values[0];}
  double p() const {return (m_values[2] - m_values[1]*m_values[1]/m_values[0]);}
  double moment(int n, const Levermore1D_weights &A, const double &us) const;
  double moment_series(int n, const Levermore1D_weights &A, const double &us) const;
  double moment_series_L(const Levermore1D_weights &A, const double &us) const;
  DenseMatrix dUdW(void) const;
  DenseMatrix dU_MBdW(void) const;
  DenseMatrix dSdW(void) const;
  DenseMatrix dSdU(void) const;
  DenseMatrix d2hda2(const Levermore1D_weights &A, const double &us) const;
  DenseMatrix d2jda2(const Levermore1D_weights &A, const double &us) const;
  Levermore1D_Vector F(const Levermore1D_weights &A) const;
  int in_sync_with(const Levermore1D_weights &A) const;
  double detector_value(const Levermore1D_weights &A, double predicted) const;
  double relative_error(const Levermore1D_cState &U2) const;
  double relaxation_time() const;
  int valid(void) const {
    return Levermore1D_pState(*this).valid();
  }

  static double m_resync_tol;

 protected:
  static double m_relaxation_time;
};

/********************************************************
 * Class: Levermore1D_weights                           *
 ********************************************************/
class Levermore1D_weights : public Levermore1D_Vector{

  public:

  Levermore1D_weights(void){}
  Levermore1D_weights(const Levermore1D_weights &A) : Levermore1D_Vector(A) {}
  explicit Levermore1D_weights(const Levermore1D_Vector &V) : Levermore1D_Vector(V) {}
  Levermore1D_weights(const double &d, const double &u, const double &p) {MaxBoltz(d,u,p);}
  explicit Levermore1D_weights(const Levermore1D_pState &W) {MaxBoltz(W); set_from_W(W);}
  explicit Levermore1D_weights(const Levermore1D_cState &U) {MaxBoltz(U); set_from_U(U);}

  /* Functions. */
  Levermore1D_Vector& operator=(const Levermore1D_Vector &V) {return Levermore1D_Vector::operator=(V);}
  void set_from_W(const Levermore1D_pState &W);
  int set_from_U(const Levermore1D_cState &U);
  double integrate_conserved_moment(int i, double us) const;
  double integrate_conserved_moment_pos(int i, double us) const;
  double integrate_conserved_moment_neg(int i, double us) const;
  double integrate_random_moment(int i, double u, double us) const;  //could us ever not equal u?
  double integrate_random_moment_pos(int i, double u, double us) const;
  double integrate_random_moment_neg(int i, double u, double us) const;

  /* Inline Functions. */
  double value_at(double v, double us) const {
    return exp(min(exponent_value_recursive(v,0)-STABILIZATION*pow((v-us),(double)(length+1)),EXP_LIMIT));
  }
  double velocity_weighted_value_at(double v, int i, double us) const {
    return pow(v,(double)i)*exp(min(exponent_value_recursive(v,0)-STABILIZATION*pow((v-us),(double)(length+1)),EXP_LIMIT));
  }
  double random_velocity_weighted_value_at(double v, double u, int i, double us) const {
    return pow(v-u,(double)i)*exp(min(exponent_value_recursive(v,0)-STABILIZATION*pow((v-us),(double)(length+1)),EXP_LIMIT));
  }

  void MaxBoltz(const Levermore1D_pState &W) {MaxBoltz(W[1], W[2], W[3]);}
  void MaxBoltz(const Levermore1D_cState &U) {MaxBoltz(U[1], U[2]/U[1], U[3]-U[2]*U[2]/U[1]);}
  void MaxBoltz(double rho, double u, double p) {
    double B(rho/(2.0*p));
    zero();
    m_values[0] = -B*u*u+log(rho*sqrt(B/PI));
    m_values[1] = 2.0*B*u;
    m_values[2] = -B;
  }

  /* Static Functions. */
  static double m() {return particle_mass;}
  static void set_particle_mass(double m) {particle_mass=m;}
  static void setgas(char* gas);

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
inline istream& operator>>(istream &in, const Levermore1D_pState &W) {
  W.input(in);
  return in;
}
inline istream& operator>>(istream &in, const Levermore1D_cState &U) {
  U.input(in);
  return in;
}
inline istream& operator>>(istream &in, const Levermore1D_weights &A) {
  A.input(in);
  return in;
}

/********************************************************
 *                Inline  Functions                     *
 ********************************************************/
inline Levermore1D_Vector Levermore1D_cState::F(const Levermore1D_weights &A) const {
  Levermore1D_Vector Flux;
  double us = m_values[1]/m_values[0];
  for(int i=1; i<=length; ++i) {
    Flux[i] = moment(i,A,us);
  }
  return Flux;
}

inline void Levermore1D_weights::setgas(char* gas) {
   if (strcmp(gas, "ZB") != 0) {
     cout << endl << "levermore1D cannot use gas: " << gas
	  << ". Using Zoidbergium instead." << endl;
   }
   set_particle_mass(MOLE_WT_ZB/(AVOGADRO*THOUSAND));
}


/********************************************************
 *          External Function Declarations              *
 ********************************************************/
extern Levermore1D_Vector FluxHLLE(const Levermore1D_cState &Ul,
				   const Levermore1D_weights &Al,
				   const double &wavespeed_l,
				   const Levermore1D_cState &Ur,
				   const Levermore1D_weights &Ar,
				   const double &wavespeed_r);

extern Levermore1D_Vector FluxKinetic(const Levermore1D_weights &Al,
				      const double &us_l,
				      const Levermore1D_weights &Ar,
				      const double &us_r);

extern Levermore1D_Vector Collision_RHS(const Levermore1D_cState &U);

extern double relaxation_time(const Levermore1D_cState &U);

inline int detector_below_tolerance(const double &detector) {
  return (detector < Levermore1D_cState::m_resync_tol);
}


#endif
