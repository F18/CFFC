/**********************************************
***********************************************
**        NumericalLibrary_Wrapper.h         **
**                                           **
** I need this wrapper to use Lucian's       **
** numerical library to integrate my moments **
**                                           **
***********************************************
***********************************************/
#ifndef _LEVERMORE1D_NUMERICALLIBRARY_WRAPPER_H
#define _LEVERMORE1D_NUMERICALLIBRARY_WRAPPER_H

#ifndef _LEVERMORE1D_STATE_INCLUDED
#include "Levermore1DState.h"
#endif //_LEVERMORE1D_STATE_INCLUDED

typedef double (Levermore1D_weights::*L1D_weights_memfunc)(double, int, double) const;

class Levermore1D_Wrapper{

 private:
  Levermore1D_Wrapper();
  const Levermore1D_weights *p_A;
  L1D_weights_memfunc p_func;
  double us;
  int n;

 public:

  Levermore1D_Wrapper(const Levermore1D_weights *const p_A_in,L1D_weights_memfunc p_func_in, int n_in, double us_in):
    p_A(p_A_in), p_func(p_func_in), n(n_in), us(us_in){}
    
  double operator() (double v) const{
    return (p_A->*p_func)(v,n,us);
  }

};


#endif
