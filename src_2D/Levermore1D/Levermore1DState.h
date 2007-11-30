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
  Levermore1D_pState(const Levermore1D_cState &U) {set_from_U(U);}
  Levermore1D_pState(const Levermore1D_weights &A) {set_from_A(A);}

  /* Functions. */
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
  Levermore1D_cState(const Levermore1D_pState &W) {set_from_W(W);}
  Levermore1D_cState(const Levermore1D_weights &A) {set_from_A(A);}

  /* Functions. */
  void Vacuum() {Levermore1D_Vector::zero();}
  void set_from_W(const Levermore1D_pState &W);
  void set_from_A(const Levermore1D_weights &A);

};

/********************************************************
 * Class: Levermore1D_weights                           *
 ********************************************************/
class Levermore1D_weights : public Levermore1D_Vector{

  public:

  Levermore1D_weights(void){}
  Levermore1D_weights(const Levermore1D_weights &A) : Levermore1D_Vector(A) {}
  Levermore1D_weights(const Levermore1D_pState &W) {set_from_W(W);}
  Levermore1D_weights(const Levermore1D_cState &U) {set_from_U(U);}

  /* Functions. */
  void Vacuum() {Levermore1D_Vector::zero();}
  void set_from_W(const Levermore1D_pState &W);
  void set_from_U(const Levermore1D_cState &U);
  double value_at(double v) {return exp(exponent_value_recursive(v,0));}

  protected:

  double exponent_value_recursive(double v, int i) {
    if(i<length) {
      return m_values[i] + v*exponent_value_recursive(v,i+1);
    } else { 
      return 1.0;
    }
  }

};

#endif
