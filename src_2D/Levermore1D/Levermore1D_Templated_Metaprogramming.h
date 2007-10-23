#ifndef _LEVERMORE1D_TEMPLATED_METAPROGRAMMING
#define _LEVERMORE1D_TEMPLATED_METAPROGRAMMING


namespace templated_metaprogramming {

  /*************************************************
   * Convert one entry in primitive (W) vector     *
   * to its coresponding conserved value (U)       *
   *************************************************/
  template<int N_t> inline double Levermore1D_add_primitive_contribution(const double *in,
									 int pf,
									 int pf_num,
									 int pf_den) {
    pf = (pf*pf_num)/pf_den;
    return (double)(pf)*in[N_t] + in[1]*Levermore1D_add_primitive_contribution<N_t-1>(in, pf, pf_num-1,pf_den+1);
  }

  template<> inline double Levermore1D_add_primitive_contribution<0>(const double *in,
								     int pf,
								     int pf_num,
								     int pf_den) {
    return in[0];
  }

  template<> inline double Levermore1D_add_primitive_contribution<1>(const double *in,
								     int pf,
								     int pf_num,
								     int pf_den) {
    return in[1]*Levermore1D_add_primitive_contribution<0>(in, pf, pf_num-1,pf_den+1);
  }

  /*************************************************
   * Convert primitive (W) vector to conserved (U) *
   *************************************************/
  template<int N> inline void Levermore1D_convert_W_to_U(double *out, const double *in) {
    out[N] = in[N] + in[1]*Levermore1D_add_primitive_contribution<N-1>(in,1,N,1);
    Levermore1D_convert_W_to_U<N-1>(out,in);
  }

  template<> inline void Levermore1D_convert_W_to_U<0>(double *out, const double *in) {
    out[0] = in[0];
  }

  template<> inline void Levermore1D_convert_W_to_U<1>(double *out, const double *in) {
    out[1] = in[1]*in[0];
    Levermore1D_convert_W_to_U<0>(out,in);
  }




}

#endif //_LEVERMORE1D_TEMPLATED_METAPROGRAMMING
