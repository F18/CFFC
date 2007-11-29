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

  friend class Levermore1D_cState;
  friend class Levermore1D_weights;

 protected:
 public:

  Levermore1D_pState(void){}
  Levermore1D_pState(const Levermore1D_pState &W) : Levermore1D_Vector(W) {}
  Levermore1D_pState(const Levermore1D_cState &U) {set_from_U(U);}
  Levermore1D_pState(const Levermore1D_weights &A) {set_from_A(A);}

  /* Functions. */
  void Vacuum() {Levermore1D_Vector::zero();}
  void set_from_U(const Levermore1D_cState &U);
  void set_from_A(const Levermore1D_weights &A);

};

/********************************************************
 * Class: Levermore1D_cState                            *
 ********************************************************/
class Levermore1D_cState : public Levermore1D_Vector{

  friend class Levermore1D_pState;
  friend class Levermore1D_weights;

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

  friend class Levermore1D_pState;
  friend class Levermore1D_cState;

  protected:
  public:

  Levermore1D_weights(void){}
  Levermore1D_weights(const Levermore1D_weights &A) : Levermore1D_Vector(A) {}
  Levermore1D_weights(const Levermore1D_pState &W) {set_from_W(W);}
  Levermore1D_weights(const Levermore1D_cState &U) {set_from_U(U);}

  /* Functions. */
  void Vacuum() {Levermore1D_Vector::zero();}
  void set_from_W(const Levermore1D_pState &W);
  void set_from_U(const Levermore1D_cState &U);

};

/********************************************************
 * Class: Levermore1D_pState: Inline functions          *
 ********************************************************/
inline void Levermore1D_pState::set_from_U(const Levermore1D_cState &U) {
  return;
}
 
inline void Levermore1D_pState::set_from_A(const Levermore1D_weights &A) {
  return;
}

/********************************************************
 * Class: Levermore1D_cState: Inline functions          *
 ********************************************************/
inline void Levermore1D_cState::set_from_W(const Levermore1D_pState &W) {
  return;
}

inline void Levermore1D_cState::set_from_A(const Levermore1D_weights &A) {
  return;
}


/********************************************************
 * Class: Levermore1D_weights: Inline functions         *
 ********************************************************/
inline void Levermore1D_weights::set_from_W(const Levermore1D_pState &W) {
  return;
}

inline void Levermore1D_weights::set_from_U(const Levermore1D_cState &U) {
  return;
}

#endif
