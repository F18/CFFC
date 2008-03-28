/*! \file  NavierStokes3DThermallyPerfectState.cc
\brief  Definition of member functions for 3D Navier Stokes solution classes
        associated with solution of compressible viscous flows of a thermally 
        perfect non-reactive or combusting mixture.
*/

/* Include NavierStokes3DThermallyPerfectState header file. */

#ifndef _NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED
#include "NavierStokes3DThermallyPerfectState.h"
#endif // NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED   

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState member functions                                  *
 ********************************************************************************************/

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::mu -- Return mixture viscosity.                  *
 ********************************************************************************************/
// Based on Wilke mixing rule
double NavierStokes3D_ThermallyPerfect_pState::mu(void) {
  double sum =0.0;
  double Temp = T();
  for (int i = 0; i < ns; i++) {
    double phi = 0.0;
    for (int j = 0; j < ns; j++) {
      if (i == 0) {
        _temp_values[j] = specdata[j].Viscosity(Temp);
      } /* endif */
      phi += (spec[j].c / specdata[j].Mol_mass())*pow(ONE + 
             sqrt(_temp_values[i]/_temp_values[j])*
             pow(specdata[j].Mol_mass()/specdata[i].Mol_mass(),0.25),2.0)/
             sqrt(EIGHT*(ONE +specdata[i].Mol_mass()/specdata[j].Mol_mass()));
    } /* endfor */
    sum += (spec[i].c * _temp_values[i]) / (specdata[i].Mol_mass() * phi);
  } /* endfor */
  return sum;
}

double NavierStokes3D_ThermallyPerfect_pState::nu(void) {
  return mu()/rho;
}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::kappa -- Return mixture thermal conductivity.    *
 ********************************************************************************************/
// Based on Mason & Saxena mixing rule
double NavierStokes3D_ThermallyPerfect_pState::kappa(void) {
  double sum = 0.0;  
  double Temp = T();
  for (int i=0; i<ns; i++) {
     double phi = 0.0;
     for (int j=0; j<ns; j++){
       if (i == 0) {
 	  _temp_values[j] = specdata[j].Viscosity(Temp);
       } /* endif */
       if (i != j) {
 	  phi += (spec[j].c / specdata[j].Mol_mass())*pow(ONE + 
                 sqrt(_temp_values[i]/_temp_values[j])*
 	         pow(specdata[j].Mol_mass()/specdata[i].Mol_mass(),0.25),2.0)/
 	         sqrt(EIGHT*(ONE +specdata[i].Mol_mass()/specdata[j].Mol_mass()));
       } /* endif */
     } /* endfor */
     sum += (specdata[i].ThermalConduct(Temp)*spec[i].c) / 
            (spec[i].c + (specdata[i].Mol_mass()) * 1.065 * phi);
  } /* endfor */
  return sum;
}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::Ds -- Species mass diffusion coefficient.        *
 ********************************************************************************************/
double NavierStokes3D_ThermallyPerfect_pState::Ds(const int &i) {
  return mu()/(rho*Schmidt[i]);
}

double NavierStokes3D_ThermallyPerfect_pState::Ds(const int &i,
                                                  const double &mu_temp) {
  return mu_temp/(rho*Schmidt[i]);
}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::Pr -- Return mixture Prandtl number.             *
 ********************************************************************************************/
double NavierStokes3D_ThermallyPerfect_pState::Pr(void) {
  //Pr = Cp*mu/k
  return Cp()*mu()/kappa();
}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::Sc -- Return species Schmidt number.             *
 ********************************************************************************************/
double NavierStokes3D_ThermallyPerfect_pState::Sc(const int &i) {
  return Schmidt[i];
}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::Le -- Return species Lewis number.               *
 ********************************************************************************************/
double NavierStokes3D_ThermallyPerfect_pState::Le(const int &i) {
  return kappa()/(rho*Cp()*Ds(i));
}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::thermal_diffusion -- Returns thermal diffusion   *
 *                                                              flux vector (due to species *
 *                                                              diffusion processes)        *
 ********************************************************************************************/
// Calculates the thermal diffusion component of the heat flux vector
// sum( hs * Ds * grad cs)
Vector3D NavierStokes3D_ThermallyPerfect_pState::
thermal_diffusion(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                  const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                  const NavierStokes3D_ThermallyPerfect_pState &dWdz) const {

   Vector3D sum, gradc;
   double Temp = T();

   sum.zero();

   for (int index = 0; index < ns; index++) {
      gradc.x = dWdx.spec[index].c;
      gradc.y = dWdy.spec[index].c;
      gradc.z = dWdz.spec[index].c;
      sum += (specdata[index].Enthalpy(Temp) + specdata[index].Heatofform())*_diff_coeff[index]*gradc;
   } /* endfor */

   return sum;

}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::thermal_diffusion_x -- Returns component of      *
 *                                                              thermal diffusion flux      *
 *                                                              vector in x-direction (due  *
 *                                                              to species diffusion)       *
 ********************************************************************************************/
Vector3D NavierStokes3D_ThermallyPerfect_pState::
thermal_diffusion_x(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                    const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                    const NavierStokes3D_ThermallyPerfect_pState &dWdz) const {

   Vector3D sum, gradc;
   double Temp = T();

   sum.zero();

   for (int index = 0; index < ns; index++) {
      gradc.x = dWdx.spec[index].c;
      sum.x += (specdata[index].Enthalpy(Temp) + specdata[index].Heatofform())*_diff_coeff[index]*gradc.x;
   } /* endfor */

   return sum;

}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::thermal_diffusion_y -- Returns component of      *
 *                                                              thermal diffusion flux      *
 *                                                              vector in y-direction (due  *
 *                                                              to species diffusion)       *
 ********************************************************************************************/
Vector3D NavierStokes3D_ThermallyPerfect_pState::
thermal_diffusion_y(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                    const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                    const NavierStokes3D_ThermallyPerfect_pState &dWdz) const {

   Vector3D sum, gradc;
   double Temp = T();

   sum.zero();

   for (int index = 0; index < ns; index++) {
      gradc.y = dWdy.spec[index].c;
      sum.y += (specdata[index].Enthalpy(Temp) + specdata[index].Heatofform())*_diff_coeff[index]*gradc.y;
   } /* endfor */

   return sum;

}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::thermal_diffusion_z -- Returns component of      *
 *                                                              thermal diffusion flux      *
 *                                                              vector in z-direction (due  *
 *                                                              to species diffusion)       *
 ********************************************************************************************/
Vector3D NavierStokes3D_ThermallyPerfect_pState::
thermal_diffusion_z(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                    const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                    const NavierStokes3D_ThermallyPerfect_pState &dWdz) const {

   Vector3D sum, gradc;
   double Temp = T();

   sum.zero();

   for (int index = 0; index < ns; index++) {
      gradc.z = dWdz.spec[index].c;
      sum.z += (specdata[index].Enthalpy(Temp) + specdata[index].Heatofform())*_diff_coeff[index]*gradc.z;
   } /* endfor */

   return sum;

}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::strain_rate -- Returns strain rate tensor.       *
 ********************************************************************************************/
Tensor3D NavierStokes3D_ThermallyPerfect_pState::
strain_rate(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
            const NavierStokes3D_ThermallyPerfect_pState &dWdy, 
            const NavierStokes3D_ThermallyPerfect_pState &dWdz) const {

   Tensor3D strain_rate_tensor;

   strain_rate_tensor.xx = dWdx.v.x;
   strain_rate_tensor.xy = (dWdy.v.x + dWdx.v.y)/TWO;
   strain_rate_tensor.xz = (dWdz.v.x + dWdx.v.z)/TWO;
   strain_rate_tensor.yy = dWdy.v.y;
   strain_rate_tensor.yz = (dWdz.v.y + dWdy.v.z)/TWO;
   strain_rate_tensor.zz = dWdz.v.z;

   return (strain_rate_tensor);

}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::rotation -- Returns rotation tensor.             *
 ********************************************************************************************/
// Note that this is an antisymmetric tensor that is stored as a symmetric tensor.
// Care must be exercised in it's use.
Tensor3D NavierStokes3D_ThermallyPerfect_pState::
rotation(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
         const NavierStokes3D_ThermallyPerfect_pState &dWdy, 
         const NavierStokes3D_ThermallyPerfect_pState &dWdz) const {

   Tensor3D rotation_tensor;

   rotation_tensor.xx = ZERO;
   rotation_tensor.xy = (dWdy.v.x - dWdx.v.y)/TWO;
   rotation_tensor.xz = (dWdz.v.x - dWdx.v.z)/TWO;
   rotation_tensor.yy = ZERO;
   rotation_tensor.yz = (dWdz.v.y - dWdy.v.z)/TWO;
   rotation_tensor.zz = ZERO;

   return (rotation_tensor);

}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::vorticity -- Returns vorticity vector.           *
 ********************************************************************************************/
Vector3D NavierStokes3D_ThermallyPerfect_pState::
vorticity(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
          const NavierStokes3D_ThermallyPerfect_pState &dWdy, 
          const NavierStokes3D_ThermallyPerfect_pState &dWdz) const {

   Vector3D vorticity_vector;

   vorticity_vector.x = dWdy.v.z - dWdz.v.y;
   vorticity_vector.y = -(dWdx.v.z - dWdz.v.x);
   vorticity_vector.z = dWdx.v.y - dWdy.v.x;

   return (vorticity_vector);

}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::tau -- Return (molecular) fluid stress tensor.   *
 ********************************************************************************************/
Tensor3D NavierStokes3D_ThermallyPerfect_pState::
tau(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
    const NavierStokes3D_ThermallyPerfect_pState &dWdy,
    const NavierStokes3D_ThermallyPerfect_pState &dWdz) {

   Tensor3D molecular_stress;
   double mu_temp = mu();
   
   molecular_stress.xx = mu_temp*(FOUR*dWdx.v.x - TWO*dWdy.v.y - TWO*dWdz.v.z)/THREE;
   molecular_stress.xy = mu_temp*(dWdx.v.y + dWdy.v.x);
   molecular_stress.xz = mu_temp*(dWdx.v.z + dWdz.v.x);
   molecular_stress.yy = mu_temp*(FOUR*dWdy.v.y - TWO*dWdx.v.x - TWO*dWdz.v.z)/THREE;
   molecular_stress.yz = mu_temp*(dWdy.v.z + dWdz.v.y);
   molecular_stress.zz = mu_temp*(FOUR*dWdz.v.z - TWO*dWdx.v.x - TWO*dWdy.v.y)/THREE;
    
   return (molecular_stress);
   
}

Tensor3D NavierStokes3D_ThermallyPerfect_pState::
tau(const double &mu_temp,
    const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
    const NavierStokes3D_ThermallyPerfect_pState &dWdy,
    const NavierStokes3D_ThermallyPerfect_pState &dWdz) {

   Tensor3D molecular_stress;
   
   molecular_stress.xx = mu_temp*(FOUR*dWdx.v.x - TWO*dWdy.v.y - TWO*dWdz.v.z)/THREE;
   molecular_stress.xy = mu_temp*(dWdx.v.y + dWdy.v.x);
   molecular_stress.xz = mu_temp*(dWdx.v.z + dWdz.v.x);
   molecular_stress.yy = mu_temp*(FOUR*dWdy.v.y - TWO*dWdx.v.x - TWO*dWdz.v.z)/THREE;
   molecular_stress.yz = mu_temp*(dWdy.v.z + dWdz.v.y);
   molecular_stress.zz = mu_temp*(FOUR*dWdz.v.z - TWO*dWdx.v.x - TWO*dWdy.v.y)/THREE;

   return (molecular_stress);
   
}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::tau_x -- Return components of (molecular) fluid  *
 *                                                  stress tensor in the x-direction.       *
 ********************************************************************************************/
Tensor3D NavierStokes3D_ThermallyPerfect_pState::
tau_x(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
      const NavierStokes3D_ThermallyPerfect_pState &dWdy,
      const NavierStokes3D_ThermallyPerfect_pState &dWdz) {

   Tensor3D molecular_stress;
   double mu_temp = mu();
   
   molecular_stress.xx = mu_temp*(FOUR*dWdx.v.x - TWO*dWdy.v.y - TWO*dWdz.v.z)/THREE;
   molecular_stress.xy = mu_temp*(dWdx.v.y + dWdy.v.x);
   molecular_stress.xz = mu_temp*(dWdx.v.z + dWdz.v.x);
    
   return (molecular_stress);
   
}

Tensor3D NavierStokes3D_ThermallyPerfect_pState::
tau_x(const double &mu_temp,
      const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
      const NavierStokes3D_ThermallyPerfect_pState &dWdy,
      const NavierStokes3D_ThermallyPerfect_pState &dWdz) {

   Tensor3D molecular_stress;
   
   molecular_stress.xx = mu_temp*(FOUR*dWdx.v.x - TWO*dWdy.v.y - TWO*dWdz.v.z)/THREE;
   molecular_stress.xy = mu_temp*(dWdx.v.y + dWdy.v.x);
   molecular_stress.xz = mu_temp*(dWdx.v.z + dWdz.v.x);

   return (molecular_stress);
   
}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::tau_y -- Return components of (molecular) fluid  *
 *                                                  stress tensor in the y-direction.       *
 ********************************************************************************************/
Tensor3D NavierStokes3D_ThermallyPerfect_pState::
tau_y(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
      const NavierStokes3D_ThermallyPerfect_pState &dWdy,
      const NavierStokes3D_ThermallyPerfect_pState &dWdz) {

   Tensor3D molecular_stress;
   double mu_temp = mu();
   
   molecular_stress.xy = mu_temp*(dWdx.v.y + dWdy.v.x);
   molecular_stress.yy = mu_temp*(FOUR*dWdy.v.y - TWO*dWdx.v.x - TWO*dWdz.v.z)/THREE;
   molecular_stress.yz = mu_temp*(dWdy.v.z + dWdz.v.y);
    
   return (molecular_stress);
   
}

Tensor3D NavierStokes3D_ThermallyPerfect_pState::
tau_y(const double &mu_temp,
      const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
      const NavierStokes3D_ThermallyPerfect_pState &dWdy,
      const NavierStokes3D_ThermallyPerfect_pState &dWdz) {

   Tensor3D molecular_stress;
   
   molecular_stress.xy = mu_temp*(dWdx.v.y + dWdy.v.x);
   molecular_stress.yy = mu_temp*(FOUR*dWdy.v.y - TWO*dWdx.v.x - TWO*dWdz.v.z)/THREE;
   molecular_stress.yz = mu_temp*(dWdy.v.z + dWdz.v.y);

   return (molecular_stress);
   
}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::tau_z -- Return components of (molecular) fluid  *
 *                                                  stress tensor in the z-direction.       *
 ********************************************************************************************/
Tensor3D NavierStokes3D_ThermallyPerfect_pState::
tau_z(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
      const NavierStokes3D_ThermallyPerfect_pState &dWdy,
      const NavierStokes3D_ThermallyPerfect_pState &dWdz) {

   Tensor3D molecular_stress;
   double mu_temp = mu();
   
   molecular_stress.xz = mu_temp*(dWdx.v.z + dWdz.v.x);
   molecular_stress.yz = mu_temp*(dWdy.v.z + dWdz.v.y);
   molecular_stress.zz = mu_temp*(FOUR*dWdz.v.z - TWO*dWdx.v.x - TWO*dWdy.v.y)/THREE;
    
   return (molecular_stress);
   
}

Tensor3D NavierStokes3D_ThermallyPerfect_pState::
tau_z(const double &mu_temp,
      const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
      const NavierStokes3D_ThermallyPerfect_pState &dWdy,
      const NavierStokes3D_ThermallyPerfect_pState &dWdz) {

   Tensor3D molecular_stress;
   
   molecular_stress.xz = mu_temp*(dWdx.v.z + dWdz.v.x);
   molecular_stress.yz = mu_temp*(dWdy.v.z + dWdz.v.y);
   molecular_stress.zz = mu_temp*(FOUR*dWdz.v.z - TWO*dWdx.v.x - TWO*dWdy.v.y)/THREE;

   return (molecular_stress);
   
}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::q -- Return (molecular) heat flux vector.        *
 ********************************************************************************************/
// Heat flux calculation based on Fourier's law of heat conduction
Vector3D NavierStokes3D_ThermallyPerfect_pState::
q(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
  const NavierStokes3D_ThermallyPerfect_pState &dWdy,
  const NavierStokes3D_ThermallyPerfect_pState &dWdz) {
  
    double Rmix = Rtot();
    double kappa_temp = kappa();    
    Vector3D heat_flux;
    
    heat_flux.x = -kappa_temp*(ONE/(rho*Rmix)) * (dWdx.p -(p/rho)*dWdx.rho);
    heat_flux.y = -kappa_temp*(ONE/(rho*Rmix)) * (dWdy.p -(p/rho)*dWdy.rho);
    heat_flux.z = -kappa_temp*(ONE/(rho*Rmix)) * (dWdz.p -(p/rho)*dWdz.rho);
        
    return (heat_flux);
   
}

Vector3D NavierStokes3D_ThermallyPerfect_pState::
q(const double &kappa_temp,
  const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
  const NavierStokes3D_ThermallyPerfect_pState &dWdy,
  const NavierStokes3D_ThermallyPerfect_pState &dWdz) {
  
    double Rmix = Rtot();
    Vector3D heat_flux;
    
    heat_flux.x = -kappa_temp*(ONE/(rho*Rmix)) * (dWdx.p -(p/rho)*dWdx.rho);
    heat_flux.y = -kappa_temp*(ONE/(rho*Rmix)) * (dWdy.p -(p/rho)*dWdy.rho);
    heat_flux.z = -kappa_temp*(ONE/(rho*Rmix)) * (dWdz.p -(p/rho)*dWdz.rho);
        
    return (heat_flux);

}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::q_x -- Return component of (molecular) heat flux *
 *                                                vector in the x-direction.                *
 ********************************************************************************************/
Vector3D NavierStokes3D_ThermallyPerfect_pState::
q_x(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
    const NavierStokes3D_ThermallyPerfect_pState &dWdy,
    const NavierStokes3D_ThermallyPerfect_pState &dWdz) {
  
    double Rmix = Rtot();
    double kappa_temp = kappa();    
    Vector3D heat_flux;
    
    heat_flux.x = -kappa_temp*(ONE/(rho*Rmix)) * (dWdx.p -(p/rho)*dWdx.rho);
        
    return (heat_flux);
   
}

Vector3D NavierStokes3D_ThermallyPerfect_pState::
q_x(const double &kappa_temp,
    const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
    const NavierStokes3D_ThermallyPerfect_pState &dWdy,
    const NavierStokes3D_ThermallyPerfect_pState &dWdz) {
  
    double Rmix = Rtot();
    Vector3D heat_flux;
    
    heat_flux.x = -kappa_temp*(ONE/(rho*Rmix)) * (dWdx.p -(p/rho)*dWdx.rho);
        
    return (heat_flux);

}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::q_y -- Return component of (molecular) heat flux *
 *                                                vector in the y-direction.                *
 ********************************************************************************************/
Vector3D NavierStokes3D_ThermallyPerfect_pState::
q_y(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
    const NavierStokes3D_ThermallyPerfect_pState &dWdy,
    const NavierStokes3D_ThermallyPerfect_pState &dWdz) {
  
    double Rmix = Rtot();
    double kappa_temp = kappa();    
    Vector3D heat_flux;
    
    heat_flux.y = -kappa_temp*(ONE/(rho*Rmix)) * (dWdy.p -(p/rho)*dWdy.rho);
        
    return (heat_flux);
   
}

Vector3D NavierStokes3D_ThermallyPerfect_pState::
q_y(const double &kappa_temp,
    const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
    const NavierStokes3D_ThermallyPerfect_pState &dWdy,
    const NavierStokes3D_ThermallyPerfect_pState &dWdz) {
  
    double Rmix = Rtot();
    Vector3D heat_flux;
    
    heat_flux.y = -kappa_temp*(ONE/(rho*Rmix)) * (dWdy.p -(p/rho)*dWdy.rho);
        
    return (heat_flux);

}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::q_z -- Return component of (molecular) heat flux *
 *                                                vector in the z-direction.                *
 ********************************************************************************************/
Vector3D NavierStokes3D_ThermallyPerfect_pState::
q_z(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
    const NavierStokes3D_ThermallyPerfect_pState &dWdy,
    const NavierStokes3D_ThermallyPerfect_pState &dWdz) {
  
    double Rmix = Rtot();
    double kappa_temp = kappa();    
    Vector3D heat_flux;
    
    heat_flux.z = -kappa_temp*(ONE/(rho*Rmix)) * (dWdz.p -(p/rho)*dWdz.rho);
        
    return (heat_flux);
   
}

Vector3D NavierStokes3D_ThermallyPerfect_pState::
q_z(const double &kappa_temp,
    const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
    const NavierStokes3D_ThermallyPerfect_pState &dWdy,
    const NavierStokes3D_ThermallyPerfect_pState &dWdz) {
  
    double Rmix = Rtot();
    Vector3D heat_flux;
    
    heat_flux.z = -kappa_temp*(ONE/(rho*Rmix)) * (dWdz.p -(p/rho)*dWdz.rho);
        
    return (heat_flux);

}

/**********************************************************************
 ***************** VISCOUS FLUX VECTORS *******************************
 **********************************************************************/

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::Fv -- Viscous flux (x-direction).                *
/********************************************************************************************/
NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_pState::
Fv(const NavierStokes3D_ThermallyPerfect_pState &dWdx,
   const NavierStokes3D_ThermallyPerfect_pState &dWdy,
   const NavierStokes3D_ThermallyPerfect_pState &dWdz) {
   
   NavierStokes3D_ThermallyPerfect_cState Temp;
   
   double mu_temp = mu();
   double kappa_temp = kappa();
   
   Tensor3D molecular_stress;
   Vector3D heat_flux;
   
   molecular_stress = tau_x(mu_temp, dWdx, dWdy, dWdz);
   heat_flux = q_x(kappa_temp, dWdx, dWdy, dWdz);
   
   for (int index = 0; index < ns; ++index) {
      _diff_coeff[index] = Ds(index, mu_temp);
   } /* endfor */

   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux -= rho*thermal_diffusion_x(dWdx, dWdy, dWdz);
   
   Temp.rho= ZERO;
   Temp.rhov.x = molecular_stress.xx;
   Temp.rhov.y = molecular_stress.xy;
   Temp.rhov.z = molecular_stress.xz;
   Temp.E = v.x*molecular_stress.xx + v.y*molecular_stress.xy +
            v.z*molecular_stress.xz - heat_flux.x;
   
   for (int index = 0; index < ns; ++index) {
      Temp.rhospec[index] = rho*_diff_coeff[index]*dWdx.spec[index].c;
   } /* endfor */
   
   return (Temp);
   
}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::Fvx -- Viscous flux (x-direction).               *
/********************************************************************************************/
NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_pState::
Fvx(const NavierStokes3D_ThermallyPerfect_pState &dWdx,
    const NavierStokes3D_ThermallyPerfect_pState &dWdy,
    const NavierStokes3D_ThermallyPerfect_pState &dWdz) {
   
   NavierStokes3D_ThermallyPerfect_cState Temp;
   
   double mu_temp = mu();
   double kappa_temp = kappa();
   
   Tensor3D molecular_stress;
   Vector3D heat_flux;
   
   molecular_stress = tau_x(mu_temp, dWdx, dWdy, dWdz);
   heat_flux = q_x(kappa_temp, dWdx, dWdy, dWdz);
   
   for (int index = 0; index < ns; ++index) {
      _diff_coeff[index] = Ds(index, mu_temp);
   } /* endfor */

   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux -= rho*thermal_diffusion_x(dWdx, dWdy, dWdz);
   
   Temp.rho= ZERO;
   Temp.rhov.x = molecular_stress.xx;
   Temp.rhov.y = molecular_stress.xy;
   Temp.rhov.z = molecular_stress.xz;
   Temp.E = v.x*molecular_stress.xx + v.y*molecular_stress.xy +
            v.z*molecular_stress.xz - heat_flux.x;
   
   for (int index = 0; index < ns; ++index) {
      Temp.rhospec[index] = rho*_diff_coeff[index]*dWdx.spec[index].c;
   } /* endfor */
   
   return (Temp);
   
}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::Fvy -- Viscous flux (y-direction).               *
/********************************************************************************************/
NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_pState::
Fvy(const NavierStokes3D_ThermallyPerfect_pState &dWdx,
    const NavierStokes3D_ThermallyPerfect_pState &dWdy,
    const NavierStokes3D_ThermallyPerfect_pState &dWdz) {

   NavierStokes3D_ThermallyPerfect_cState Temp;
   
   double mu_temp = mu();
   double kappa_temp = kappa();
   
   Tensor3D molecular_stress;
   Vector3D heat_flux;
   
   molecular_stress = tau_y(mu_temp, dWdx, dWdy, dWdz);
   heat_flux = q_y(kappa_temp, dWdx, dWdy, dWdz);
   
   for (int index = 0; index < ns; ++index) {
      _diff_coeff[index] = Ds(index, mu_temp);
   } /* endfor */

   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux -= rho*thermal_diffusion_y(dWdx, dWdy, dWdz);
   
   Temp.rho= ZERO;
   Temp.rhov.x = molecular_stress.xy;
   Temp.rhov.y = molecular_stress.yy;
   Temp.rhov.z = molecular_stress.yz;
   Temp.E = v.x*molecular_stress.xy+ v.y*molecular_stress.yy +
            v.z*molecular_stress.yz - heat_flux.y;
   
   for (int index = 0; index < ns; ++index) {
      Temp.rhospec[index] = rho*_diff_coeff[index]*dWdy.spec[index].c;
   } /* endfor */
   
   return (Temp);
}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::Fvz -- Viscous flux (z-direction).               *
/********************************************************************************************/
NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_pState::
Fvz(const NavierStokes3D_ThermallyPerfect_pState &dWdx,
    const NavierStokes3D_ThermallyPerfect_pState &dWdy,
    const NavierStokes3D_ThermallyPerfect_pState &dWdz) {
   
   NavierStokes3D_ThermallyPerfect_cState Temp;
   
   double mu_temp = mu();
   double kappa_temp = kappa();
   
   Tensor3D molecular_stress;
   Vector3D heat_flux;
   
   molecular_stress = tau_z(mu_temp, dWdx, dWdy, dWdz);
   heat_flux = q_z(kappa_temp, dWdx, dWdy, dWdz);
   
   for (int index = 0; index < ns; ++index) {
      _diff_coeff[index] = Ds(index, mu_temp);
   } /* endfor */

   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux -=  rho*thermal_diffusion_z(dWdx, dWdy, dWdz);
   
   Temp.rho= ZERO;
   Temp.rhov.x = molecular_stress.xz;
   Temp.rhov.y = molecular_stress.yz;
   Temp.rhov.z = molecular_stress.zz;
   Temp.E = v.x*molecular_stress.xz+ v.y*molecular_stress.yz +
            v.z*molecular_stress.zz - heat_flux.z;
   
   for (int index = 0; index < ns; ++index) {
      Temp.rhospec[index] = rho*_diff_coeff[index]*dWdz.spec[index].c;
   } /* endfor */

   return (Temp);
   
}

/************************************************************************
 *************** NUMERICAL EVALUATION OF VISCOUS FLUXES *****************
 ************************************************************************/

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::FluxViscous_n -- Viscous flux (n-direction).     *
/********************************************************************************************/
NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_pState::
FluxViscous_n(const NavierStokes3D_ThermallyPerfect_pState &Wl,
              const NavierStokes3D_ThermallyPerfect_pState &Wr,
              const NavierStokes3D_ThermallyPerfect_pState &Wc,
              const NavierStokes3D_ThermallyPerfect_pState &Wc_Neigbor,
              const NavierStokes3D_ThermallyPerfect_pState &dWdx,
              const NavierStokes3D_ThermallyPerfect_pState &dWdy,
              const NavierStokes3D_ThermallyPerfect_pState &dWdz,
              const NavierStokes3D_ThermallyPerfect_pState &dWdx_Neigbor,
              const NavierStokes3D_ThermallyPerfect_pState &dWdy_Neigbor,
              const NavierStokes3D_ThermallyPerfect_pState &dWdz_Neigbor,
              const Vector3D &norm, 
              const Vector3D &ts, 
              const double &deltad,
              const double &Volume, 
              const double &Volume_Neigbor) {
   
   // construct the gradients on the cell interface (surface) 
   // based on Hybrid Average Gradient-Diamond-Path Approach
   // Grad Phi = (Phi_c_neigbor - Phi_c)/ds *norm/(norm . ts)
   //            + [bar (Grad Phi) - bar (Grad Phi) . ts norm/(norm.ts)]

   // weighted factor based on volume
   double alpha = Volume/(Volume + Volume_Neigbor);
   
   NavierStokes3D_ThermallyPerfect_pState dWdx_Weighted, 
                                          dWdy_Weighted, 
                                          dWdz_Weighted, 
                                          dWdx_face, 
                                          dWdy_face, 
                                          dWdz_face, 
                                          Grad_middle_term,
                                          W_face;

   dWdx_Weighted = alpha*dWdx + (1.0 - alpha)*dWdx_Neigbor;
   dWdy_Weighted = alpha*dWdy + (1.0 - alpha)*dWdy_Neigbor;
   dWdz_Weighted = alpha*dWdz + (1.0 - alpha)*dWdz_Neigbor;
   
   // a weighted term  
   Grad_middle_term = dWdx_Weighted*ts.x + dWdy_Weighted*ts.y + dWdz_Weighted*ts.z;
      
   // gradients of primitive variables on the face
   dWdx_face = (Wc_Neigbor - Wc)/deltad *norm.x/dot(norm, ts) + 
               (dWdx_Weighted -  Grad_middle_term*norm.x/dot(norm, ts));
   dWdy_face = (Wc_Neigbor - Wc)/deltad *norm.y/dot(norm, ts) + 
               (dWdy_Weighted -  Grad_middle_term*norm.y/dot(norm, ts));
   dWdz_face = (Wc_Neigbor - Wc)/deltad *norm.z/dot(norm, ts) + 
               (dWdz_Weighted -  Grad_middle_term*norm.z/dot(norm, ts));
   
   W_face = HALF*(Wl + Wr);

   if (fabs(norm.y) < TOLER && fabs(norm.z) < TOLER) {
      return (W_face.Fvx(dWdx_face, dWdy_face, dWdz_face)*norm.x);

   } else if (fabs(norm.x) < TOLER && fabs(norm.z) < TOLER) {
      return (W_face.Fvy(dWdx_face, dWdy_face, dWdz_face)*norm.y);

   } else if (fabs(norm.x) < TOLER && fabs(norm.y) < TOLER) {
      return (W_face.Fvz(dWdx_face, dWdy_face, dWdz_face)*norm.z);

   } else {
      return (W_face.Fvx(dWdx_face, dWdy_face, dWdz_face)*norm.x +
              W_face.Fvy(dWdx_face, dWdy_face, dWdz_face)*norm.y +
              W_face.Fvz(dWdx_face, dWdy_face, dWdz_face)*norm.z);

   } /* endif */
   
}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_cState member functions                                  *
 ********************************************************************************************/

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_cState::mu -- Return mixture viscosity.                  *
 ********************************************************************************************/
// Based on Wilke mixing rule
double NavierStokes3D_ThermallyPerfect_cState::mu(void) {
   double sum = ZERO;
   double Temp = T();
   for (int i = 0; i < ns; i++) {
     double phi = 0.0;
     for (int j = 0; j < ns; j++) {
       phi += ((rhospec[j].c/rho) / specdata[j].Mol_mass())*
              pow(1.0 + sqrt(specdata[i].Viscosity(Temp)/specdata[j].Viscosity(Temp))*
              pow(specdata[j].Mol_mass()/specdata[i].Mol_mass(),0.25),2.0)/
              sqrt(8.0*(1.0 +specdata[i].Mol_mass()/specdata[j].Mol_mass()));
     } /* endfor */
     sum += ((rhospec[i].c/rho)* specdata[i].Viscosity(Temp))/(specdata[i].Mol_mass()*phi);
   } /* endfor */
   return sum;
}

double NavierStokes3D_ThermallyPerfect_cState::nu(void) {
  return mu()/rho;
}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_cState::kappa -- Return mixture thermal conductivity.    *
 ********************************************************************************************/
// Based on Mason & Saxena mixing rule
double NavierStokes3D_ThermallyPerfect_cState::kappa(void) {
  return ZERO;
}

/********************************************************************************************
 * NavierStokes3D_ThermallyPerfect_cState::Ds -- Species mass diffusion coefficient.        *
 ********************************************************************************************/
double NavierStokes3D_ThermallyPerfect_cState::Ds(const int &i) {
  return mu()/(rho*Schmidt[i]);
}

double NavierStokes3D_ThermallyPerfect_cState::Ds(const int &i,
                                                  const double &mu_temp) {
  return mu_temp/(rho*Schmidt[i]);
}
