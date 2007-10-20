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

#ifndef _LEVERMORE1D_VECTOR_INCLUDED
#include "Levermore1DVector.h"
#endif //_LEVERMORE1D_VECTOR_INCLUDED

/* Define the classes. */
template<int N_moments>
class Levermore1D_cState;

template<int N_moments>
class Levermore1D_weights;

/********************************************************
 * Class: Levermore1D_pState                            *
 ********************************************************/
template<int N_moments>
class Levermore1D_pState : public Levermore1D_Vector<N_moments>{
 private:
 public:

  Levermore1D_pState(void){}
  Levermore1D_pState(const Levermore1D_pState &W) {copy_from(W);}
  Levermore1D_pState(const Levermore1D_cState<N_moments> &U) {from_U(U);}
  Levermore1D_pState(const Levermore1D_weights<N_moments> &A) {from_A(A);}

  /* Functions. */
  void Vacuum() {Levermore1D_Vector<N_moments>::zero();}
  void from_U(const Levermore1D_cState<N_moments> &U);
  void from_A(const Levermore1D_weights<N_moments> &A);

};

/********************************************************
 * Class: Levermore1D_cState                            *
 ********************************************************/
template<int N_moments>
class Levermore1D_cState : public Levermore1D_Vector<N_moments>{
  private:
  public:

  Levermore1D_cState(void){}
  Levermore1D_cState(const Levermore1D_cState &U) {copy_from(U);}
  Levermore1D_cState(const Levermore1D_pState<N_moments> &W) {from_W(W);}
  Levermore1D_cState(const Levermore1D_weights<N_moments> &A) {from_A(A);}

  /* Functions. */
  void Vacuum() {Levermore1D_Vector<N_moments>::zero();}
  void from_W(const Levermore1D_pState<N_moments> &W);
  void from_A(const Levermore1D_weights<N_moments> &A);

};

/********************************************************
 * Class: Levermore1D_weights                           *
 ********************************************************/
template<int N_moments>
class Levermore1D_weights : public Levermore1D_Vector<N_moments>{
  private:
  public:

  Levermore1D_weights(void){}
  Levermore1D_weights(const Levermore1D_weights &A) {copy_from(A);}
  Levermore1D_weights(const Levermore1D_pState<N_moments> &W) {from_W(W);}
  Levermore1D_weights(const Levermore1D_cState<N_moments> &U) {from_U(U);}

  /* Functions. */
  void Vacuum() {Levermore1D_Vector<N_moments>::zero();}
  void from_W(const Levermore1D_pState<N_moments> &W);
  void from_U(const Levermore1D_cState<N_moments> &U);

};

#endif
