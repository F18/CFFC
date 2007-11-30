/* Levermore1DState.cc:  Source file for Levermore1D_XState classes. */

#include "Levermore1DState.h"

/********************************************************
 * Function: Levermore1D_pState::set_from_U             *
 *                                                      *
 * Convert a converved vector into a primitive one.     *
 *                                                      *
 ********************************************************/
void Levermore1D_pState::set_from_U(const Levermore1D_cState &U) {

  for(int i=0; i<length; ++i) {
    if(i>1) {
      m_values[i] = U[i+1]-conserved_extras(i);
    } else if(i==1) {
      m_values[i] = U[2]/U[1];
    } else {  //i=0
      m_values[i] = U[1];
    }
  }

  return;
}

/********************************************************
 * Function: Levermore1D_cState::set_from_W             *
 *                                                      *
 * Convert a primitive vector into a conserved one.     *
 *                                                      *
 ********************************************************/
void Levermore1D_cState::set_from_W(const Levermore1D_pState &W) {

  for(int i=0; i<length; ++i) {
    if(i>1) {
      m_values[i] = W[i+1]+W.conserved_extras(i);
    } else if(i==1) {
      m_values[i] = W[1]*W[2];
    } else {  //i=0
      m_values[i] = W[1];
    }
  }

  return;
}

/********************************************************
 * Function: Levermore1D_pState::conserved_extras       *
 *                                                      *
 * Find the difference between the "i" entry in the     *
 * primitive and conserved vectors.  This requires all  *
 * previous entries in the primitive vector be up to    *
 * date.                                                *
 *                                                      *
 ********************************************************/
double Levermore1D_pState::conserved_extras(int i) const {
  int pf(1);
  return m_values[1]*conserved_extras_recursive(i-1,pf,i,1);
}

/********************************************************
 * Function: Levermore1D_pState::                       *
 *              conserved_extras_recursive              *
 *                                                      *
 * Used for the conserved_extras function.              *
 *                                                      *
 ********************************************************/
double Levermore1D_pState::conserved_extras_recursive(int i, int &pf, int pf_num, int pf_den) const {
  if(i==1) {
    return m_values[1]*m_values[0];
  } else {
    pf = (pf*pf_num)/pf_den;
    return (double)pf * m_values[i]  + m_values[1] * conserved_extras_recursive(i-1,pf,pf_num-1,pf_den+1);
  }
}
