#ifndef _LEVERMORE1D_TEMPLATED_METAPROGRAMMING
#define _LEVERMORE1D_TEMPLATED_METAPROGRAMMING


namespace templated_metaprogramming {

  /*************************************************
   * Calculates extra stuff in conserved vector    *
   *************************************************/
  template<int N_t> inline double Levermore1D_conserved_extras(const double *primitive,
							       unsigned int pf,
							       unsigned int pf_num,
							       unsigned int pf_den) {
    pf = (pf*pf_num)/pf_den;
    return (double)(pf)*primitive[N_t] + primitive[1]*Levermore1D_conserved_extras<N_t-1>(primitive, pf, pf_num-1,pf_den+1);
  }

  template<> inline double Levermore1D_conserved_extras<0>(const double *primitive,
							   unsigned int pf,
							   unsigned int pf_num,
							   unsigned int pf_den) {
    return primitive[0];
  }

  template<> inline double Levermore1D_conserved_extras<1>(const double *primitive,
							   unsigned int pf,
							   unsigned int pf_num,
							   unsigned int pf_den) {
    return primitive[1]*Levermore1D_conserved_extras<0>(primitive, pf, pf_num-1,pf_den+1);
  }

  /*************************************************
   * Convert conserved (U) vector to primitive (W) *
   *************************************************/
  template<int N> inline void Levermore1D_convert_U_to_W(double *primitive, const double *conserved) {
    //convert moments from lowest- to highest-order
    Levermore1D_convert_U_to_W<N-1>(primitive,conserved);
    primitive[N] = conserved[N] - primitive[1]*Levermore1D_conserved_extras<N-1>(primitive,1,N,1);
  }

  template<> inline void Levermore1D_convert_U_to_W<0>(double *primitive, const double *conserved) {
    primitive[0] = conserved[0];
  }

  template<> inline void Levermore1D_convert_U_to_W<1>(double *primitive, const double *conserved) {
    Levermore1D_convert_U_to_W<0>(primitive,conserved);
    primitive[1] = conserved[1]/conserved[0];
  }

  /*************************************************
   * Convert primitive (W) vector to conserved (U) *
   *************************************************/
  template<int N> inline void Levermore1D_convert_W_to_U(double *conserved, const double *primitive) {
    //convert moments from lowest- to highest-order
    Levermore1D_convert_W_to_U<N-1>(conserved,primitive);
    conserved[N] = primitive[N] + primitive[1]*Levermore1D_conserved_extras<N-1>(primitive,1,N,1);
  }

  template<> inline void Levermore1D_convert_W_to_U<0>(double *conserved, const double *primitive) {
    conserved[0] = primitive[0];
  }

  template<> inline void Levermore1D_convert_W_to_U<1>(double *conserved, const double *primitive) {
    Levermore1D_convert_W_to_U<0>(conserved,primitive);
    conserved[1] = primitive[1]*primitive[0];
  }




}

#endif //_LEVERMORE1D_TEMPLATED_METAPROGRAMMING
