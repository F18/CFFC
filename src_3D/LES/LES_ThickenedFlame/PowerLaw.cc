#ifndef _POWER_LAW_INCLUDED
#include "PowerLaw.h"
#endif // _POWER_LAW_INCLUDED



//  Efficiency function of Charlette, Meneveau and Veynante
double PowerLaw::efficiency_function(const double &lam_speed,
				     const double &lam_thickness,
				     const double &cell_size,
				     const double &lapl_vor,
				     const double &rho,
				     const double &mu) const {

  double fu, fdelta, fRe, a, b, Delta, u_p, sfs_Re, nu;
  nu = (mu/rho)/*/(TF*WF)*/;  // Original viscosity, unmodified by the thickening 
  Delta = filter_size(lam_thickness);
  u_p = u_prime( cell_size /*lam_thickness*/, lapl_vor);
  sfs_Re = u_p*Delta/nu; 

  b = 1.4;
  if (u_p == 0.0){
    a = 0.8 - 0.2*exp(-0.01*Delta/lam_thickness);
    fu = 0.0;
    fRe = 0.0;
  } else {
    a = 0.6 + 0.2*exp(-0.1*u_p/lam_speed) - 0.2*exp(-0.01*Delta/lam_thickness);
    fu = FOUR * sqrt(27.0*KOLMOGOROV_CONSTANT/110.0) * (18.0*KOLMOGOROV_CONSTANT/55.0)
         * sqr(u_p/lam_speed);
    fRe = sqrt( (9.0/55.0) * exp(-3.0*KOLMOGOROV_CONSTANT*pow(PI, 4.0/3.0)/(2.0*sfs_Re)) ) 
          * sqrt(sfs_Re);
  }

  fdelta = sqrt( (27.0*KOLMOGOROV_CONSTANT*pow(PI,4.0/3.0)/110.0) 
		 * (pow(Delta/lam_thickness, 4.0/3.0) - 1.0) );

  double pow_fu, pow_fdelta, pow_fRe, pow_u_delta;
  if(Delta/lam_thickness < ONE){
    fRe = 0.0;
    cout << "\n Delta over laminar thickness less than unity.";
    cout << "\nDelta: " << Delta <<"\tlaminar_thickness: " << lam_thickness;
  }

  
  pow_fu = (fu == 0.0) ?  0.0 : pow(fu,-a);
  pow_fdelta = (fdelta == 0.0) ?  0.0 : pow(fdelta,-a);
  pow_fRe = (fRe == 0.0) ?  0.0 : pow(fRe,-b);
  pow_u_delta = pow_fu + pow_fdelta;
  pow_u_delta = (pow_u_delta > 0.0) ?  pow(pow_u_delta, b/a) : 0.0;

  return pow(pow_u_delta + pow_fRe, -1/b);
}


// Mean curvature (inverse inner cutoff scale)
double PowerLaw::mean_curvature(const double &lam_speed,
				const double &lam_thickness,
				const double &cell_size,
				const double &lapl_vor,
				const double &rho,
				const double &mu) const {

  double u_p, effunc;
  u_p = u_prime( cell_size /*lam_thickness*/ , lapl_vor);
  effunc = efficiency_function(lam_speed, lam_thickness, cell_size,
			       lapl_vor, rho, mu);

  if (u_p <= 0.0  ||  effunc <= 0.0) {
    return 0.0;
  } else {
    double mean_c;
    mean_c = u_p*effunc/(lam_speed*filter_size(lam_thickness)); 
    return mean_c;
  }
}


// Wrinkling factor
void PowerLaw::wrinkling_factor(const double &lam_speed,
				const double &lam_thickness,
				const double &cell_size,
				const double &lapl_vor,
				const double &rho,
				const double &mu) {

  const double beta = 0.5;  // exponent of the power-law expression
  double mean_c;
  mean_c = mean_curvature(lam_speed, lam_thickness, cell_size, lapl_vor, rho, mu);

  WF = pow(1.0 + min(filter_size(lam_thickness)/lam_thickness, 
		     filter_size(lam_thickness)*mean_c), beta );  
}  
