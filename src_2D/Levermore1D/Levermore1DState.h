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

#ifndef _LEVERMORE1D_TEMPLATED_METAPROGRAMMING
#include "Levermore1D_Templated_Metaprogramming.h"
#endif //_LEVERMORE1D_TEMPLATED_METAPROGRAMMING

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

  friend class Levermore1D_cState<N_moments>;
  friend class Levermore1D_weights<N_moments>;

 protected:
 public:

  Levermore1D_pState(void){}
  Levermore1D_pState(const Levermore1D_pState &W) {copy_from(W);}
  Levermore1D_pState(const Levermore1D_cState<N_moments> &U) {set_from_U(U);}
  Levermore1D_pState(const Levermore1D_weights<N_moments> &A) {set_from_A(A);}

  /* Functions. */
  void Vacuum() {Levermore1D_Vector<N_moments>::zero();}
  void set_from_U(const Levermore1D_cState<N_moments> &U);
  void set_from_A(const Levermore1D_weights<N_moments> &A);

};

/********************************************************
 * Class: Levermore1D_cState                            *
 ********************************************************/
template<int N_moments>
class Levermore1D_cState : public Levermore1D_Vector<N_moments>{

  friend class Levermore1D_pState<N_moments>;
  friend class Levermore1D_weights<N_moments>;

  protected:
  public:

  Levermore1D_cState(void){}
  Levermore1D_cState(const Levermore1D_cState &U) {copy_from(U);}
  Levermore1D_cState(const Levermore1D_pState<N_moments> &W) {set_from_W(W);}
  Levermore1D_cState(const Levermore1D_weights<N_moments> &A) {set_from_A(A);}

  /* Functions. */
  void Vacuum() {Levermore1D_Vector<N_moments>::zero();}
  void set_from_W(const Levermore1D_pState<N_moments> &W);
  void set_from_A(const Levermore1D_weights<N_moments> &A);

};

/********************************************************
 * Class: Levermore1D_weights                           *
 ********************************************************/
template<int N_moments>
class Levermore1D_weights : public Levermore1D_Vector<N_moments>{

  friend class Levermore1D_pState<N_moments>;
  friend class Levermore1D_cState<N_moments>;

  protected:
  public:

  Levermore1D_weights(void){}
  Levermore1D_weights(const Levermore1D_weights &A) {copy_from(A);}
  Levermore1D_weights(const Levermore1D_pState<N_moments> &W) {set_from_W(W);}
  Levermore1D_weights(const Levermore1D_cState<N_moments> &U) {set_from_U(U);}

  /* Functions. */
  void Vacuum() {Levermore1D_Vector<N_moments>::zero();}
  void set_from_W(const Levermore1D_pState<N_moments> &W);
  void set_from_U(const Levermore1D_cState<N_moments> &U);

};

/********************************************************
 * Class: Levermore1D_pState: Inline functions          *
 ********************************************************/
template<int N_moments>
inline void Levermore1D_pState<N_moments>::set_from_U(const Levermore1D_cState<N_moments> &U) {
  templated_metaprogramming::Levermore1D_convert_U_to_W<N_moments-1>(Levermore1D_Vector<N_moments>::m_values, U.m_values);
  return;
}

template<int N_moments>
inline void Levermore1D_pState<N_moments>::set_from_A(const Levermore1D_weights<N_moments> &A) {
  return;
}

/********************************************************
 * Class: Levermore1D_cState: Inline functions          *
 ********************************************************/
template<int N_moments>
inline void Levermore1D_cState<N_moments>::set_from_W(const Levermore1D_pState<N_moments> &W) {
  templated_metaprogramming::Levermore1D_convert_W_to_U<N_moments-1>(Levermore1D_Vector<N_moments>::m_values, W.m_values);
  return;
}

template<int N_moments>
inline void Levermore1D_cState<N_moments>::set_from_A(const Levermore1D_weights<N_moments> &A) {
  return;
}


/********************************************************
 * Class: Levermore1D_weights: Inline functions         *
 ********************************************************/
template<int N_moments>
inline void Levermore1D_weights<N_moments>::set_from_W(const Levermore1D_pState<N_moments> &W) {
  return;
}

template<int N_moments>
inline void Levermore1D_weights<N_moments>::set_from_U(const Levermore1D_cState<N_moments> &U) {
  return;
}

#endif
