
namespace templated_metaprogramming {


  template<int N_term> inline void copy_double_array(double *out, const double *in_1) {
    out[N_term] = in_1[N_term];
    copy_double_array<N_term-1>(out,in_1);
  }

  template<> inline void copy_double_array<0>(double *out, const double *in_1) {
    out[0] = in_1[0];
  }

  template<int N_term> inline void set_double_array(double *out, const double &val) {
    out[N_term] = val;
    set_double_array<N_term-1>(out,val);
  }

  template<> inline void set_double_array<0>(double *out, const double &val) {
    out[0] = val;
  }

  template<int N_term> inline void add_double_array(double *out,
						    const double *in_1,
						    const double *in_2) {
    out[N_term] = in_1[N_term] + in_2[N_term];
    add_double_array<N_term-1>(out,in_1,in_2);
  }

  template<> inline void add_double_array<0>(double *out,
					     const double *in_1,
					     const double *in_2) {
    out[0] = in_1[0] + in_2[0];
  }

  template<int N_term> inline void addto_double_array(double *out, const double *in_1) {
    out[N_term] += in_1[N_term];
    addto_double_array<N_term-1>(out,in_1);
  }

  template<> inline void addto_double_array<0>(double *out, const double *in_1) {
    out[0] += in_1[0];
  }

  template<int N_term> inline void subtract_double_array(double *out,
							 const double *in_1,
							 const double *in_2) {
    out[N_term] = in_1[N_term] - in_2[N_term];
    subtract_double_array<N_term-1>(out,in_1,in_2);
  }

  template<> inline void subtract_double_array<0>(double *out,
						  const double *in_1,
						  const double *in_2) {
    out[0] = in_1[0] - in_2[0];
  }

  template<int N_term> inline void subtractfrom_double_array(double *out, const double *in_1) {
    out[N_term] -= in_1[N_term];
    subtractfrom_double_array<N_term-1>(out,in_1);
  }

  template<> inline void subtractfrom_double_array<0>(double *out, const double *in_1) {
    out[0] -= in_1[0];
  }


  template<int N_term> inline double dot_product_double_array(double *in_1, double *in_2) {
    return ( in_1[N_term]*in_2[N_term] + dot_product_double_array<N_term-1>(in_1,in_2) );
  }

  template<> inline double dot_product_double_array<0>(double *in_1, double *in_2) {
    return ( in_1[0]*in_2[0] );
  }

  template<int N_term> inline void termwise_multiply_double_array(double *out,
								  const double *in_1,
								  const double *in_2) {
    out[N_term] = in_1[N_term] * in_2[N_term];
    termwise_multiply_double_array<N_term-1>(out,in_1,in_2);
  }

  template<> inline void termwise_multiply_double_array<0>(double *out,
							   const double *in_1,
							   const double *in_2) {
    out[0] = in_1[0] * in_2[0];
  }


}
