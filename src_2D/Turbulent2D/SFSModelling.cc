#ifndef _SFS_MODELLING_INCLUDED
#include "SFSModelling.h"
#endif



/****************************************************
 *    Members of the SubfilterScaleModels class     *
 ****************************************************/


// Initialize static members
double SubfilterScaleModels::Smagorinsky_coef = 0.18;
double SubfilterScaleModels::CI_Yoshizawa = 0.005;
double SubfilterScaleModels::CV_coef = 0.086; 
double SubfilterScaleModels::CEPS_coef = 0.845;


void SubfilterScaleModels::set_SFSmodel_constants(const double &Smagorinsky_c, 
						  const double &Yoshizawa_c) {
  Smagorinsky_coef = Smagorinsky_c;
  CI_Yoshizawa = Yoshizawa_c;
}

void SubfilterScaleModels::set_SFSmodel_constants_k(const double &CV_keq, 
						    const double &CEPS_keq) {
  CV_coef = CV_keq;
  CEPS_coef = CEPS_keq;
} 

double SubfilterScaleModels::abs_strain_rate(const Tensor2D &strain_rate) const {
  
  // S[i,j]*S[i,j]
  double SS = sqr(strain_rate.xx) + sqr(strain_rate.yy) + 2.0*sqr(strain_rate.xy);

  // sqrt(2*S*S)
  return sqrt(2.0*SS);

}

double SubfilterScaleModels::eddy_viscosity_Smagorinsky(const Tensor2D &strain_rate,
							const double &filter_width) const {

  double nu_smagorinsky = sqr(Smagorinsky_coef*filter_width)*abs_strain_rate(strain_rate);
  return nu_smagorinsky;  
}

double SubfilterScaleModels::eddy_viscosity_k(const double &k, 
					      const double &filter_width) const {

  return CV_coef*sqrt(k)*filter_width; 

}

double SubfilterScaleModels::sfs_k_Yoshizawa(const Tensor2D &strain_rate,
					     const double &filter_width) const {

  double trace_no_tworho = CI_Yoshizawa*sqr(filter_width*abs_strain_rate(strain_rate));
  return (trace_no_tworho);
}
