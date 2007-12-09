/*! \file  FANS3DThermallyPerfectState.cc
\brief  Definition of member functions for 3D Favre-averaged Navier-Stokes solution classes
        associated with solution of compressible turbulent flows of a thermally 
        perfect non-reactive or combusting mixture.
*/

/* Include FANS3DThermallyPerfectState header file. */

#ifndef _FANS3D_THERMALLYPERFECT_STATE_INCLUDED
#include "FANS3DThermallyPerfectState.h"
#endif // _FANS3D_THERMALLYPERFECT_STATE_INCLUDED   

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState member functions                                *
 *****************************************************************************************/

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::Realizable_Solution_Check -- Check physical    *
 *                                                                      validity of      *
 *                                                                      solution state.  *
 *****************************************************************************************/
bool FANS3D_ThermallyPerfect_KOmega_pState::Realizable_Solution_Check(void) {
   if (rho <= ZERO || !negative_speccheck() || p <= ZERO ||
       k <= ZERO || omega <= ZERO) {    
      cout << "\n " << CFFC_Name() 
           << " ERROR: Primitive solution state has a negative density, pressure, mass fractions,"
           << " turbulent kinetic energy, and/or dissipation rate.\n";
      return false;
   } else {
      return true;
   } /* endif */
} 

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::E -- Return total energy of the mixture.       *
 *****************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_pState::E(void) const {   
   return (Euler3D_ThermallyPerfect_pState::E() + rho*k);
}

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::H -- Return total enthalpy of the mixture.     *
 *****************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_pState::H(void) const {
   return (Euler3D_ThermallyPerfect_pState::H() + rho*k);
}

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::Hs -- Return total mixture sensible enthalpy.  *
 *****************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_pState::Hs(void) const {
   return (Euler3D_ThermallyPerfect_pState::Hs() + rho*k);
}

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::a_t -- Return mixture sound speed (including   *
 *                                               turbulent kinetic energy).              *
 *****************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_pState::a_t(void) {
   double aa = sqr(a());
   aa += (TWO/THREE)*k*g();
   return sqrt(aa);
}

double FANS3D_ThermallyPerfect_KOmega_pState::a_t(void) const {
   double aa = sqr(a());
   aa += (TWO/THREE)*k*g();
   return sqrt(aa);
}

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::mu_t -- Return eddy (turbulent) viscosity.     *
 *****************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_pState::mu_t(void) {
   return (rho*k/max(TOLER, omega));
}

/********************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::kappa_t -- Return turbulent thermal conductivity. *
 ********************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_pState::kappa_t(void) {
  return (mu_t()*Cp()/Pr_t()); 
}

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::Ds_t -- Return species turbulent diffusion     *
 *                                                coefficient.                           *
 *****************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_pState::Ds_t(const int &i) {
   return (mu_t()/(rho*Sc_t()));
}

double FANS3D_ThermallyPerfect_KOmega_pState::Ds_t(const int &i,
                                                   const double &mu_t_temp) {
   return (mu_t_temp/(rho*Sc_t()));
}

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::Pr_t -- Return turbulent Prandtl number.       *
 *****************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_pState::Pr_t(void) {
   return (0.9); 
}

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::Sc_t -- Return turbulent Schmidt number.       *
 *****************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_pState::Sc_t(void) {
   return (1.0);
}

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::Le_t -- Return turbulent Lewis number.         *
 *****************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_pState::Le_t(void) {
  return (Sc_t()/Pr_t());
}

/**********************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::tau_t -- Return Reynolds (turbulent) stress tensor. *
 **********************************************************************************************/
Tensor3D FANS3D_ThermallyPerfect_KOmega_pState::
tau_t(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {
   
   Tensor3D Reynolds_stress;
   double mu_t_temp = mu_t();

   Reynolds_stress.xx = mu_t_temp*(FOUR*dWdx.v.x - TWO*dWdy.v.y - TWO*dWdz.v.z)/THREE - (TWO/THREE)*rho*k;
   Reynolds_stress.xy = mu_t_temp*(dWdx.v.y + dWdy.v.x);
   Reynolds_stress.xz = mu_t_temp*(dWdx.v.z + dWdz.v.x);
   Reynolds_stress.yy = mu_t_temp*(FOUR*dWdy.v.y - TWO*dWdx.v.x - TWO*dWdz.v.z)/THREE - (TWO/THREE)*rho*k;
   Reynolds_stress.yz = mu_t_temp*(dWdy.v.z + dWdz.v.y);
   Reynolds_stress.zz = mu_t_temp*(FOUR*dWdz.v.z - TWO*dWdx.v.x - TWO*dWdy.v.y)/THREE - (TWO/THREE)*rho*k;

   return (Reynolds_stress);

}

Tensor3D FANS3D_ThermallyPerfect_KOmega_pState::
tau_t(const double &mu_t_temp,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {
   
   Tensor3D Reynolds_stress;

   Reynolds_stress.xx = mu_t_temp*(FOUR*dWdx.v.x - TWO*dWdy.v.y - TWO*dWdz.v.z)/THREE - (TWO/THREE)*rho*k;
   Reynolds_stress.xy = mu_t_temp*(dWdx.v.y + dWdy.v.x);
   Reynolds_stress.xz = mu_t_temp*(dWdx.v.z + dWdz.v.x);
   Reynolds_stress.yy = mu_t_temp*(FOUR*dWdy.v.y - TWO*dWdx.v.x - TWO*dWdz.v.z)/THREE - (TWO/THREE)*rho*k;
   Reynolds_stress.yz = mu_t_temp*(dWdy.v.z + dWdz.v.y);
   Reynolds_stress.zz = mu_t_temp*(FOUR*dWdz.v.z - TWO*dWdx.v.x - TWO*dWdy.v.y)/THREE - (TWO/THREE)*rho*k;

   return (Reynolds_stress);

}

/**********************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::tau_t_x -- Return components of Reynolds (turbulent)*
 *                                                   stress tensor in the x-direction.        *
 **********************************************************************************************/
Tensor3D FANS3D_ThermallyPerfect_KOmega_pState::
tau_t_x(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
        const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
        const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {
   
   Tensor3D Reynolds_stress;
   double mu_t_temp = mu_t();

   Reynolds_stress.xx = mu_t_temp*(FOUR*dWdx.v.x - TWO*dWdy.v.y - TWO*dWdz.v.z)/THREE - (TWO/THREE)*rho*k;
   Reynolds_stress.xy = mu_t_temp*(dWdx.v.y + dWdy.v.x);
   Reynolds_stress.xz = mu_t_temp*(dWdx.v.z + dWdz.v.x);

   return (Reynolds_stress);

}

Tensor3D FANS3D_ThermallyPerfect_KOmega_pState::
tau_t_x(const double &mu_t_temp,
        const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
        const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
        const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {
   
   Tensor3D Reynolds_stress;

   Reynolds_stress.xx = mu_t_temp*(FOUR*dWdx.v.x - TWO*dWdy.v.y - TWO*dWdz.v.z)/THREE - (TWO/THREE)*rho*k;
   Reynolds_stress.xy = mu_t_temp*(dWdx.v.y + dWdy.v.x);
   Reynolds_stress.xz = mu_t_temp*(dWdx.v.z + dWdz.v.x);

   return (Reynolds_stress);

}

/**********************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::tau_t_y -- Return components of Reynolds (turbulent)*
 *                                                   stress tensor in the y-direction.        *
 **********************************************************************************************/
Tensor3D FANS3D_ThermallyPerfect_KOmega_pState::
tau_t_y(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
        const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
        const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {
   
   Tensor3D Reynolds_stress;
   double mu_t_temp = mu_t();

   Reynolds_stress.xy = mu_t_temp*(dWdx.v.y + dWdy.v.x);
   Reynolds_stress.yy = mu_t_temp*(FOUR*dWdy.v.y - TWO*dWdx.v.x - TWO*dWdz.v.z)/THREE - (TWO/THREE)*rho*k;
   Reynolds_stress.yz = mu_t_temp*(dWdy.v.z + dWdz.v.y);

   return (Reynolds_stress);

}

Tensor3D FANS3D_ThermallyPerfect_KOmega_pState::
tau_t_y(const double &mu_t_temp,
        const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
        const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
        const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {
   
   Tensor3D Reynolds_stress;

   Reynolds_stress.xy = mu_t_temp*(dWdx.v.y + dWdy.v.x);
   Reynolds_stress.yy = mu_t_temp*(FOUR*dWdy.v.y - TWO*dWdx.v.x - TWO*dWdz.v.z)/THREE - (TWO/THREE)*rho*k;
   Reynolds_stress.yz = mu_t_temp*(dWdy.v.z + dWdz.v.y);

   return (Reynolds_stress);

}

/**********************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::tau_t_z -- Return components of Reynolds (turbulent)*
 *                                                   stress tensor in the z-direction.        *
 **********************************************************************************************/
Tensor3D FANS3D_ThermallyPerfect_KOmega_pState::
tau_t_z(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
        const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
        const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {
   
   Tensor3D Reynolds_stress;
   double mu_t_temp = mu_t();

   Reynolds_stress.xz = mu_t_temp*(dWdx.v.z + dWdz.v.x);
   Reynolds_stress.yz = mu_t_temp*(dWdy.v.z + dWdz.v.y);
   Reynolds_stress.zz = mu_t_temp*(FOUR*dWdz.v.z - TWO*dWdx.v.x - TWO*dWdy.v.y)/THREE - (TWO/THREE)*rho*k;

   return (Reynolds_stress);

}

Tensor3D FANS3D_ThermallyPerfect_KOmega_pState::
tau_t_z(const double &mu_t_temp,
        const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
        const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
        const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {
   
   Tensor3D Reynolds_stress;

   Reynolds_stress.xz = mu_t_temp*(dWdx.v.z + dWdz.v.x);
   Reynolds_stress.yz = mu_t_temp*(dWdy.v.z + dWdz.v.y);
   Reynolds_stress.zz = mu_t_temp*(FOUR*dWdz.v.z - TWO*dWdx.v.x - TWO*dWdy.v.y)/THREE - (TWO/THREE)*rho*k;

   return (Reynolds_stress);

}

/**********************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::q_t -- Return turbulent heat flux vector.           *
 **********************************************************************************************/
Vector3D FANS3D_ThermallyPerfect_KOmega_pState::
q_t(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
    const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
    const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {
  
    double Rmix = Rtot();
    double kappa_t_temp = kappa_t();
    Vector3D heat_flux;
     
    heat_flux.x = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdx.p -(p/rho)*dWdx.rho);
    heat_flux.y = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdy.p -(p/rho)*dWdy.rho);
    heat_flux.z = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdz.p -(p/rho)*dWdz.rho);
        
    return (heat_flux);

}

Vector3D FANS3D_ThermallyPerfect_KOmega_pState::
q_t(const double &kappa_t_temp,
    const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
    const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
    const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {
  
    double Rmix = Rtot();
    Vector3D heat_flux;
     
    heat_flux.x = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdx.p -(p/rho)*dWdx.rho);
    heat_flux.y = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdy.p -(p/rho)*dWdy.rho);
    heat_flux.z = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdz.p -(p/rho)*dWdz.rho);
        
    return (heat_flux);

}

/**********************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::q_t_x -- Return component of turbulent heat flux    *
 *                                                 vector in the x-direction.                 *
 **********************************************************************************************/
Vector3D FANS3D_ThermallyPerfect_KOmega_pState::
q_t_x(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {
  
    double Rmix = Rtot();
    double kappa_t_temp = kappa_t();
    Vector3D heat_flux;
     
    heat_flux.x = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdx.p -(p/rho)*dWdx.rho);

    return (heat_flux);

}

Vector3D FANS3D_ThermallyPerfect_KOmega_pState::
q_t_x(const double &kappa_t_temp,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {
  
    double Rmix = Rtot();
    Vector3D heat_flux;
     
    heat_flux.x = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdx.p -(p/rho)*dWdx.rho);

    return (heat_flux);

}

/**********************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::q_t_y -- Return component of turbulent heat flux    *
 *                                                 vector in the y-direction.                 *
 **********************************************************************************************/
Vector3D FANS3D_ThermallyPerfect_KOmega_pState::
q_t_y(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {
  
    double Rmix = Rtot();
    double kappa_t_temp = kappa_t();
    Vector3D heat_flux;
     
    heat_flux.y = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdy.p -(p/rho)*dWdy.rho);
        
    return (heat_flux);

}

Vector3D FANS3D_ThermallyPerfect_KOmega_pState::
q_t_y(const double &kappa_t_temp,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {
  
    double Rmix = Rtot();
    Vector3D heat_flux;
     
    heat_flux.y = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdy.p -(p/rho)*dWdy.rho);
        
    return (heat_flux);

}

/**********************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::q_t_z -- Return component of turbulent heat flux    *
 *                                                 vector in the z-direction.                 *
 **********************************************************************************************/
Vector3D FANS3D_ThermallyPerfect_KOmega_pState::
q_t_z(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {
  
    double Rmix = Rtot();
    double kappa_t_temp = kappa_t();
    Vector3D heat_flux;
     
    heat_flux.z = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdz.p -(p/rho)*dWdz.rho);
        
    return (heat_flux);

}

Vector3D FANS3D_ThermallyPerfect_KOmega_pState::
q_t_z(const double &kappa_t_temp,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
      const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {
  
    double Rmix = Rtot();
    Vector3D heat_flux;
     
    heat_flux.z = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdz.p -(p/rho)*dWdz.rho);
        
    return (heat_flux);

}

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::U -- Return conserved solution state vector.   *
 *****************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::U(void) {
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.rho = rho;
   Temp.rhov = rhov();
   Temp.E = E();
   Temp.rhok = rho*k;
   Temp.rhoomega = rho*omega;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i] = rho*spec[i];
   } /* endfor */
   return Temp;
}

FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::U(void) const {
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.rho = rho;
   Temp.rhov = rhov();
   Temp.E = E();
   Temp.rhok = rho*k;
   Temp.rhoomega = rho*omega;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i] = rho*spec[i];
   } /* endfor */
   return Temp;
}

FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::
U(const FANS3D_ThermallyPerfect_KOmega_pState &W) {
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.rho = W.rho;
   Temp.rhov = W.rhov();
   Temp.E = W.E();
   Temp.rhok = W.rho*W.k;
   Temp.rhoomega = W.rho*W.omega;
   for (int i = 0; i < W.ns; i++) {
      Temp.rhospec[i] = W.rho*W.spec[i];
   } /* endfor */
   return Temp;
}

/***********************************************************************
 ***************** INVISCID FLUX VECTORS *******************************
 ***********************************************************************/

/***************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::F -- Inviscid flux (x-direction).            *
 ***************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::F(void) {
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.rho = rho*v.x;
   Temp.rhov.x = rho*sqr(v.x) + p + (TWO/THREE)*rho*k;
   Temp.rhov.y = rho*v.x*v.y;
   Temp.rhov.z = rho*v.x*v.z;
   Temp.E = v.x*(H() + (TWO/THREE)*rho*k);
   Temp.rhok = rho*v.x*k;
   Temp.rhoomega = rho*v.x*omega;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i].c = rho*v.x*spec[i].c;
   } /* endfor */
   return (Temp);
}

FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::F(void) const {
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.rho = rho*v.x;
   Temp.rhov.x = rho*sqr(v.x) + p + (TWO/THREE)*rho*k;
   Temp.rhov.y = rho*v.x*v.y;
   Temp.rhov.z = rho*v.x*v.z;
   Temp.E = v.x*(H() + (TWO/THREE)*rho*k);
   Temp.rhok = rho*v.x*k;
   Temp.rhoomega = rho*v.x*omega;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i].c = rho*v.x*spec[i].c;
   } /* endfor */
   return (Temp);
}

/***************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::Fx -- Inviscid flux (x-direction).           *
 ***************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::Fx(void) {
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.rho = rho*v.x;
   Temp.rhov.x = rho*sqr(v.x) + p + (TWO/THREE)*rho*k;
   Temp.rhov.y = rho*v.x*v.y;
   Temp.rhov.z = rho*v.x*v.z;
   Temp.E = v.x*(H() + (TWO/THREE)*rho*k);
   Temp.rhok = rho*v.x*k;
   Temp.rhoomega = rho*v.x*omega;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i].c = rho*v.x*spec[i].c;
   } /* endfor */
   return (Temp);
}

FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::Fx(void) const {
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.rho = rho*v.x;
   Temp.rhov.x = rho*sqr(v.x) + p + (TWO/THREE)*rho*k;
   Temp.rhov.y = rho*v.x*v.y;
   Temp.rhov.z = rho*v.x*v.z;
   Temp.E = v.x*(H() + (TWO/THREE)*rho*k);
   Temp.rhok = rho*v.x*k;
   Temp.rhoomega = rho*v.x*omega;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i].c = rho*v.x*spec[i].c;
   } /* endfor */
   return (Temp);
}

/***************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::Fy -- Inviscid flux (y-direction).           *
 ***************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::Fy(void) {
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.rho = rho*v.y;
   Temp.rhov.x = rho*v.x*v.y;
   Temp.rhov.y = rho*sqr(v.y) + p + (TWO/THREE)*rho*k;
   Temp.rhov.z = rho*v.y*v.z;
   Temp.E = v.y*(H() + (TWO/THREE)*rho*k);
   Temp.rhok = rho*v.y*k;
   Temp.rhoomega = rho*v.y*omega;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i].c = rho*v.y*spec[i].c;
   } /* endfor */
   return (Temp);
}

FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::Fy(void) const {
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.rho = rho*v.y;
   Temp.rhov.x = rho*v.x*v.y;
   Temp.rhov.y = rho*sqr(v.y) + p + (TWO/THREE)*rho*k;
   Temp.rhov.z = rho*v.y*v.z;
   Temp.E = v.y*(H() + (TWO/THREE)*rho*k);
   Temp.rhok = rho*v.y*k;
   Temp.rhoomega = rho*v.y*omega;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i].c = rho*v.y*spec[i].c;
   } /* endfor */
   return (Temp);
}

/***************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::Fz -- Inviscid flux (z-direction).           *
 ***************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::Fz(void) {
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.rho = rho*v.y;
   Temp.rhov.x = rho*v.x*v.y;
   Temp.rhov.y = rho*sqr(v.y) + p + (TWO/THREE)*rho*k;
   Temp.rhov.z = rho*v.y*v.z;
   Temp.E = v.y*(H() + (TWO/THREE)*rho*k);
   Temp.rhok = rho*v.y*k;
   Temp.rhoomega = rho*v.y*omega;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i].c = rho*v.y*spec[i].c;
   } /* endfor */
   return (Temp);
}

FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::Fz(void) const {
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.rho = rho*v.y;
   Temp.rhov.x = rho*v.x*v.y;
   Temp.rhov.y = rho*sqr(v.y) + p + (TWO/THREE)*rho*k;
   Temp.rhov.z = rho*v.y*v.z;
   Temp.E = v.y*(H() + (TWO/THREE)*rho*k);
   Temp.rhok = rho*v.y*k;
   Temp.rhoomega = rho*v.y*omega;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i].c = rho*v.y*spec[i].c;
   } /* endfor */
   return (Temp);
}

/**********************************************************************
 ***************** VISCOUS FLUX VECTORS *******************************
 **********************************************************************/

/***************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::Fv -- Viscous flux (x-direction).            *
 ***************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::
Fv(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx,
   const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
   const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {
   
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   
   double mu_temp = mu();
   double mu_t_temp = mu_t();
   double kappa_temp = kappa();
   double kappa_t_temp = kappa_t();
   
   Tensor3D fluid_stress;
   Vector3D heat_flux;
   
   fluid_stress = tau_t_x(mu_temp + mu_t_temp, dWdx, dWdy, dWdz);
   heat_flux = q_t_x(kappa_temp + kappa_t_temp, dWdx, dWdy, dWdz);
   
   for (int index = 0; index < ns; ++index) {
      _diff_coeff[index] = Ds(index, mu_temp)+Ds_t(index, mu_t_temp);
   } /* endfor */
   
   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux -= rho*thermal_diffusion_x(dWdx, dWdy, dWdz);
   
   Temp.rho= ZERO;
   Temp.rhov.x = fluid_stress.xx + (TWO/THREE)*rho*k;
   Temp.rhov.y = fluid_stress.xy;
   Temp.rhov.z = fluid_stress.xz;
   Temp.E = v.x*(fluid_stress.xx + (TWO/THREE)*rho*k) +
            v.y*fluid_stress.xy  + 
            v.z*fluid_stress.xz - 
            heat_flux.x +
            (mu_temp + mu_t_temp*k_omega_model.sigma_star)*dWdx.k;
   Temp.rhok = (mu_temp + mu_t_temp*k_omega_model.sigma_star)*dWdx.k;
   Temp.rhoomega = (mu_temp + mu_t_temp*k_omega_model.sigma)*dWdx.omega;
   
   for (int index = 0; index < ns; ++index) {
      Temp.rhospec[index] = rho*_diff_coeff[index]*dWdx.spec[index].c;
   } /* endfor */

   return (Temp);

}

/***************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::Fvx -- Viscous flux (x-direction).           *
 ***************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::
Fvx(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx,
    const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
    const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {
   
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   
   double mu_temp = mu();
   double mu_t_temp = mu_t();
   double kappa_temp = kappa();
   double kappa_t_temp = kappa_t();
   
   Tensor3D fluid_stress;
   Vector3D heat_flux;
   
   fluid_stress = tau_t_x(mu_temp + mu_t_temp, dWdx, dWdy, dWdz);
   heat_flux = q_t_x(kappa_temp + kappa_t_temp, dWdx, dWdy, dWdz);
   
   for (int index = 0; index < ns; ++index) {
      _diff_coeff[index] = Ds(index, mu_temp)+Ds_t(index, mu_t_temp);
   } /* endfor */
   
   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux -= rho*thermal_diffusion_x(dWdx, dWdy, dWdz);
   
   Temp.rho= ZERO;
   Temp.rhov.x = fluid_stress.xx + (TWO/THREE)*rho*k;
   Temp.rhov.y = fluid_stress.xy;
   Temp.rhov.z = fluid_stress.xz;
   Temp.E = v.x*(fluid_stress.xx + (TWO/THREE)*rho*k) +
            v.y*fluid_stress.xy  + 
            v.z*fluid_stress.xz - 
            heat_flux.x +
            (mu_temp + mu_t_temp*k_omega_model.sigma_star)*dWdx.k;
   Temp.rhok = (mu_temp + mu_t_temp*k_omega_model.sigma_star)*dWdx.k;
   Temp.rhoomega = (mu_temp + mu_t_temp*k_omega_model.sigma)*dWdx.omega;
   
   for (int index = 0; index < ns; ++index) {
      Temp.rhospec[index] = rho*_diff_coeff[index]*dWdx.spec[index].c;
   } /* endfor */

   return (Temp);

}

/***************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::Fvy -- Viscous flux (y-direction).           *
 ***************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::
Fvy(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx,
    const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
    const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {

   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   
   double mu_temp = mu();
   double mu_t_temp = mu_t();
   double kappa_temp = kappa();
   double kappa_t_temp = kappa_t();
   
   Tensor3D fluid_stress;
   Vector3D heat_flux;
   
   fluid_stress = tau_t_y(mu_temp + mu_t_temp, dWdx, dWdy, dWdz);
   heat_flux = q_t_y(kappa_temp + kappa_t_temp, dWdx, dWdy, dWdz);
   
   for (int index = 0; index < ns; ++index) {
      _diff_coeff[index] = Ds(index, mu_temp)+Ds_t(index, mu_t_temp);
   } /* endfor */
   
   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux -= rho*thermal_diffusion_y(dWdx, dWdy, dWdz);
   
   Temp.rho= ZERO;
   Temp.rhov.x = fluid_stress.xy;
   Temp.rhov.y = fluid_stress.yy + (TWO/THREE)*rho*k ;
   Temp.rhov.z = fluid_stress.yz;
   Temp.E = v.x*fluid_stress.xy +
            v.y*(fluid_stress.yy + (TWO/THREE)*rho*k) +
            v.z*fluid_stress.yz - 
            heat_flux.y +
            (mu_temp + mu_t_temp*k_omega_model.sigma_star)*dWdy.k;
   Temp.rhok = (mu_temp + mu_t_temp*k_omega_model.sigma_star)*dWdy.k;
   Temp.rhoomega = (mu_temp + mu_t_temp*k_omega_model.sigma)*dWdy.omega;
   
   for (int index = 0; index<ns; ++index) {
      Temp.rhospec[index] = rho*_diff_coeff[index]*dWdy.spec[index].c;
   } /* endif */

   return (Temp);

}

/***************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::Fvz -- Viscous flux (z-direction).           *
 ***************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::
Fvz(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx,
    const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
    const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {

   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   
   double mu_temp = mu();
   double mu_t_temp = mu_t();
   double kappa_temp = kappa();
   double kappa_t_temp = kappa_t();
   
   Tensor3D fluid_stress;
   Vector3D heat_flux;
   
   fluid_stress = tau_t_z(mu_temp + mu_t_temp, dWdx, dWdy, dWdz);
   heat_flux = q_t_z(kappa_temp + kappa_t_temp, dWdx, dWdy, dWdz);

   for (int index = 0; index < ns; ++index) {
      _diff_coeff[index] = Ds(index, mu_temp)+Ds_t(index, mu_t_temp);
   } /* endfor */
   
   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux -= rho*thermal_diffusion_z(dWdx, dWdy, dWdz);
   
   Temp.rho= ZERO;
   Temp.rhov.x = fluid_stress.xz;
   Temp.rhov.y = fluid_stress.yz;
   Temp.rhov.z = fluid_stress.zz + (TWO/THREE)*rho*k;
   Temp.E = v.x*fluid_stress.xz +
            v.y*fluid_stress.yz +
            v.z*(fluid_stress.zz + (TWO/THREE)*rho*k) - 
            heat_flux.z +
            (mu_temp + mu_t_temp*k_omega_model.sigma_star)*dWdz.k;
   Temp.rhok = (mu_temp + mu_t_temp*k_omega_model.sigma_star)*dWdz.k;
   Temp.rhoomega = (mu_temp + mu_t_temp*k_omega_model.sigma)*dWdz.omega;
   
   for (int index = 0; index<ns; ++index) {
      Temp.rhospec[index] = rho*_diff_coeff[index]*dWdz.spec[index].c;
   } /* endfor */

   return (Temp);

}

/************************************************************************
 ********************** EIGENVALUES *************************************
 ************************************************************************/

/************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::lambda -- Eigenvalue(s) (x-direction).    *
 ************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::lambda(void){
  double cc = a_t();
  FANS3D_ThermallyPerfect_KOmega_pState Temp;
  Temp.rho = v.x - cc;
  Temp.v.x = v.x;
  Temp.v.y = v.x;
  Temp.v.z = v.x;
  Temp.p = v.x + cc;
  Temp.k = v.x;
  Temp.omega = v.x;
  for (int i = 0; i < ns; i++) {
     Temp.spec[i].c = v.x;
  } /* endfor */
  return (Temp);
}

FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::lambda(void) const {
  double cc = a_t();
  FANS3D_ThermallyPerfect_KOmega_pState Temp;
  Temp.rho = v.x - cc;
  Temp.v.x = v.x;
  Temp.v.y = v.x;
  Temp.v.z = v.x;
  Temp.p = v.x + cc;
  Temp.k = v.x;
  Temp.omega = v.x;
  for (int i = 0; i < ns; i++) {
     Temp.spec[i].c = v.x;
  } /* endfor */
  return (Temp);
}

/************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::lambda_x -- Eigenvalue(s) (x-direction).  *
 ************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::lambda_x(void){
  double cc = a_t();
  FANS3D_ThermallyPerfect_KOmega_pState Temp;
  Temp.rho = v.x - cc;
  Temp.v.x = v.x;
  Temp.v.y = v.x;
  Temp.v.z = v.x;
  Temp.p = v.x + cc;
  Temp.k = v.x;
  Temp.omega = v.x;
  for (int i = 0; i < ns; i++) {
     Temp.spec[i].c = v.x;
  } /* endfor */
  return (Temp);
}

FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::lambda_x(void) const {
  double cc = a_t();
  FANS3D_ThermallyPerfect_KOmega_pState Temp;
  Temp.rho = v.x - cc;
  Temp.v.x = v.x;
  Temp.v.y = v.x;
  Temp.v.z = v.x;
  Temp.p = v.x + cc;
  Temp.k = v.x;
  Temp.omega = v.x;
  for (int i = 0; i < ns; i++) {
     Temp.spec[i].c = v.x;
  } /* endfor */
  return (Temp);
}

/************************************************************************
 ********************* EIGENVECTORS *************************************
 ************************************************************************/

/*********************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::rc -- Conserved right eigenvector (x-direction).   *
 *********************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::rc(const int &index) {
   double cc;
   switch(index){  
   case 1:
      cc = a_t();
      return (FANS3D_ThermallyPerfect_KOmega_cState(ONE, 
                                                    v.x-cc, v.y, v.z, 
                                                    (H() + (TWO/THREE)*rho*k)/rho-v.x*cc, 
                                                    k, omega, spec));
   case 2:
      return (FANS3D_ThermallyPerfect_KOmega_cState(ONE, 
                                                    v.x, v.y, v.z, 
                                                    (H() + (TWO/THREE)*rho*k)/rho-Cp()*T(), 
                                                    k, omega, spec)); 
   case 3:
      return (FANS3D_ThermallyPerfect_KOmega_cState(ZERO, 
                                                    ZERO, rho, ZERO, 
                                                    rho*v.y, 
                                                    ZERO, ZERO, ZERO));
   case 4:
      return (FANS3D_ThermallyPerfect_KOmega_cState(ZERO, 
                                                    ZERO, ZERO, rho, 
                                                    rho*v.z, 
                                                    ZERO, ZERO, ZERO));
   case 5:
      cc = a_t();
      return (FANS3D_ThermallyPerfect_KOmega_cState(ONE, 
                                                    v.x+cc, v.y, v.z, 
                                                    (H() + (TWO/THREE)*rho*k)/rho+v.x*cc, 
                                                    k, omega, spec));
   case 6:
      return (FANS3D_ThermallyPerfect_KOmega_cState(ZERO, 
                                                    ZERO, ZERO, ZERO, 
                                                    rho*(ONE-TWO/THREE*(ONE/(g()-ONE))), 
                                                    rho, ZERO, ZERO));
   case 7:
      return (FANS3D_ThermallyPerfect_KOmega_cState(ZERO, 
                                                    ZERO, ZERO, ZERO, 
                                                    ZERO, 
                                                    ZERO, rho, ZERO));
   default :
      FANS3D_ThermallyPerfect_KOmega_cState rr(ZERO);
      double RTOT = Rtot();
      double TEMP = p/(rho*RTOT);      
      rr.E = rho*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Enthalpy(TEMP) + 
                   specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Heatofform() 
                   - Cp(TEMP)*TEMP*specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_FANS3D_VAR_EXTRA+1)].Rs()/RTOT);
      double PHI = specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Enthalpy(TEMP)-
                   specdata[num_vars-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Enthalpy(TEMP)-
                   (Cp(TEMP) -RTOT)*TEMP*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_FANS3D_VAR_EXTRA+1)].Rs() - specdata[num_vars - (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_FANS3D_VAR_EXTRA+1)].Rs())/RTOT;
      rr.rhospec[index - (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].c = rho*PHI; 
      return rr;
   };
}

FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::rc(const int &index) const {
   double cc;
   switch(index){  
   case 1:
      cc = a_t();
      return (FANS3D_ThermallyPerfect_KOmega_cState(ONE, 
                                                    v.x-cc, v.y, v.z, 
                                                    (H() + (TWO/THREE)*rho*k)/rho-v.x*cc, 
                                                    k, omega, spec));
   case 2:
      return (FANS3D_ThermallyPerfect_KOmega_cState(ONE, 
                                                    v.x, v.y, v.z, 
                                                    (H() + (TWO/THREE)*rho*k)/rho-Cp()*T(), 
                                                    k, omega, spec)); 
   case 3:
      return (FANS3D_ThermallyPerfect_KOmega_cState(ZERO, 
                                                    ZERO, rho, ZERO, 
                                                    rho*v.y, 
                                                    ZERO, ZERO, ZERO));
   case 4:
      return (FANS3D_ThermallyPerfect_KOmega_cState(ZERO, 
                                                    ZERO, ZERO, rho, 
                                                    rho*v.z, 
                                                    ZERO, ZERO, ZERO));
   case 5:
      cc = a_t();
      return (FANS3D_ThermallyPerfect_KOmega_cState(ONE, 
                                                    v.x+cc, v.y, v.z, 
                                                    (H() + (TWO/THREE)*rho*k)/rho+v.x*cc, 
                                                    k, omega, spec));
   case 6:
      return (FANS3D_ThermallyPerfect_KOmega_cState(ZERO, 
                                                    ZERO, ZERO, ZERO, 
                                                    rho*(ONE-TWO/THREE*(ONE/(g()-ONE))), 
                                                    rho, ZERO, ZERO));
   case 7:
      return (FANS3D_ThermallyPerfect_KOmega_cState(ZERO, 
                                                    ZERO, ZERO, ZERO, 
                                                    ZERO, 
                                                    ZERO, rho, ZERO));
   default :
      FANS3D_ThermallyPerfect_KOmega_cState rr(ZERO);
      double RTOT = Rtot();
      double TEMP = p/(rho*RTOT);      
      rr.E = rho*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Enthalpy(TEMP) + 
                   specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Heatofform() 
                   - Cp(TEMP)*TEMP*specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_FANS3D_VAR_EXTRA+1)].Rs()/RTOT);
      double PHI = specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Enthalpy(TEMP)-
                   specdata[num_vars-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Enthalpy(TEMP)-
                   (Cp(TEMP) -RTOT)*TEMP*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_FANS3D_VAR_EXTRA+1)].Rs() - specdata[num_vars - (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_FANS3D_VAR_EXTRA+1)].Rs())/RTOT;
      rr.rhospec[index - (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].c = rho*PHI; 
      return rr;
   };
}

/*********************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::rc_x -- Conserved right eigenvector (x-direction). *
 *********************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::rc_x(const int &index) {
   double cc;
   switch(index){  
   case 1:
      cc = a_t();
      return (FANS3D_ThermallyPerfect_KOmega_cState(ONE, 
                                                    v.x-cc, v.y, v.z, 
                                                    (H() + (TWO/THREE)*rho*k)/rho-v.x*cc, 
                                                    k, omega, spec));
   case 2:
      return (FANS3D_ThermallyPerfect_KOmega_cState(ONE, 
                                                    v.x, v.y, v.z, 
                                                    (H() + (TWO/THREE)*rho*k)/rho-Cp()*T(), 
                                                    k, omega, spec)); 
   case 3:
      return (FANS3D_ThermallyPerfect_KOmega_cState(ZERO, 
                                                    ZERO, rho, ZERO, 
                                                    rho*v.y, 
                                                    ZERO, ZERO, ZERO));
   case 4:
      return (FANS3D_ThermallyPerfect_KOmega_cState(ZERO, 
                                                    ZERO, ZERO, rho, 
                                                    rho*v.z, 
                                                    ZERO, ZERO, ZERO));
   case 5:
      cc = a_t();
      return (FANS3D_ThermallyPerfect_KOmega_cState(ONE, 
                                                    v.x+cc, v.y, v.z, 
                                                    (H() + (TWO/THREE)*rho*k)/rho+v.x*cc, 
                                                    k, omega, spec));
   case 6:
      return (FANS3D_ThermallyPerfect_KOmega_cState(ZERO, 
                                                    ZERO, ZERO, ZERO, 
                                                    rho*(ONE-TWO/THREE*(ONE/(g()-ONE))), 
                                                    rho, ZERO, ZERO));
   case 7:
      return (FANS3D_ThermallyPerfect_KOmega_cState(ZERO, 
                                                    ZERO, ZERO, ZERO, 
                                                    ZERO, 
                                                    ZERO, rho, ZERO));
   default :
      FANS3D_ThermallyPerfect_KOmega_cState rr(ZERO);
      double RTOT = Rtot();
      double TEMP = p/(rho*RTOT);      
      rr.E = rho*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Enthalpy(TEMP) + 
                   specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Heatofform() 
                   - Cp(TEMP)*TEMP*specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_FANS3D_VAR_EXTRA+1)].Rs()/RTOT);
      double PHI = specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Enthalpy(TEMP)-
                   specdata[num_vars-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Enthalpy(TEMP)-
                   (Cp(TEMP) -RTOT)*TEMP*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_FANS3D_VAR_EXTRA+1)].Rs() - specdata[num_vars - (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_FANS3D_VAR_EXTRA+1)].Rs())/RTOT;
      rr.rhospec[index - (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].c = rho*PHI; 
      return rr;
   };
}

FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::rc_x(const int &index) const {
   double cc;
   switch(index){  
   case 1:
      cc = a_t();
      return (FANS3D_ThermallyPerfect_KOmega_cState(ONE, 
                                                    v.x-cc, v.y, v.z, 
                                                    (H() + (TWO/THREE)*rho*k)/rho-v.x*cc, 
                                                    k, omega, spec));
   case 2:
      return (FANS3D_ThermallyPerfect_KOmega_cState(ONE, 
                                                    v.x, v.y, v.z, 
                                                    (H() + (TWO/THREE)*rho*k)/rho-Cp()*T(), 
                                                    k, omega, spec)); 
   case 3:
      return (FANS3D_ThermallyPerfect_KOmega_cState(ZERO, 
                                                    ZERO, rho, ZERO, 
                                                    rho*v.y, 
                                                    ZERO, ZERO, ZERO));
   case 4:
      return (FANS3D_ThermallyPerfect_KOmega_cState(ZERO, 
                                                    ZERO, ZERO, rho, 
                                                    rho*v.z, 
                                                    ZERO, ZERO, ZERO));
   case 5:
      cc = a_t();
      return (FANS3D_ThermallyPerfect_KOmega_cState(ONE, 
                                                    v.x+cc, v.y, v.z, 
                                                    (H() + (TWO/THREE)*rho*k)/rho+v.x*cc, 
                                                    k, omega, spec));
   case 6:
      return (FANS3D_ThermallyPerfect_KOmega_cState(ZERO, 
                                                    ZERO, ZERO, ZERO, 
                                                    rho*(ONE-TWO/THREE*(ONE/(g()-ONE))), 
                                                    rho, ZERO, ZERO));
   case 7:
      return (FANS3D_ThermallyPerfect_KOmega_cState(ZERO, 
                                                    ZERO, ZERO, ZERO, 
                                                    ZERO, 
                                                    ZERO, rho, ZERO));
   default :
      FANS3D_ThermallyPerfect_KOmega_cState rr(ZERO);
      double RTOT = Rtot();
      double TEMP = p/(rho*RTOT);      
      rr.E = rho*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Enthalpy(TEMP) + 
                   specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Heatofform() 
                   - Cp(TEMP)*TEMP*specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_FANS3D_VAR_EXTRA+1)].Rs()/RTOT);
      double PHI = specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Enthalpy(TEMP)-
                   specdata[num_vars-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Enthalpy(TEMP)-
                   (Cp(TEMP) -RTOT)*TEMP*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_FANS3D_VAR_EXTRA+1)].Rs() - specdata[num_vars - (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_FANS3D_VAR_EXTRA+1)].Rs())/RTOT;
      rr.rhospec[index - (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].c = rho*PHI; 
      return rr;
   };
}

/*********************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::lp -- Primitive left eigenvector (x-direction).    *
 *********************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::lp(const int &index) {
   double cc;
   switch(index){  
   case 1:
      cc = a_t();
      return (FANS3D_ThermallyPerfect_KOmega_pState(k/(THREE*cc*cc), 
                                                    -HALF*rho/cc, ZERO, ZERO, 
                                                    HALF/(cc*cc), 
                                                    rho/(THREE*cc*cc), ZERO, ZERO));
   case 2:
      cc = a_t();
      return (FANS3D_ThermallyPerfect_KOmega_pState(ONE - TWO*k/(THREE*cc*cc), 
                                                    ZERO, ZERO, ZERO, 
                                                    -ONE/(cc*cc), 
                                                    -TWO*rho/(THREE*cc*cc), ZERO, ZERO));
   case 3 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(ZERO, 
                                                    ZERO, ONE, ZERO, 
                                                    ZERO, 
                                                    ZERO, ZERO, ZERO));
   case 4 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(ZERO, 
                                                    ZERO, ZERO, ONE, 
                                                    ZERO, 
                                                    ZERO, ZERO, ZERO));
   case 5 :
      cc = a_t();
      return (FANS3D_ThermallyPerfect_KOmega_pState(k/(THREE*cc*cc), 
                                                    HALF*rho/cc, ZERO, ZERO, 
                                                    HALF/(cc*cc), 
                                                    rho/(THREE*cc*cc), ZERO, ZERO));
   case 6 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(ZERO, 
                                                    ZERO, ZERO, ZERO, 
                                                    ZERO, 
                                                    ONE, ZERO, ZERO));
   case 7 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(ZERO, 
                                                    ZERO, ZERO, ZERO, 
                                                    ZERO, 
                                                    ZERO, ONE, ZERO));
   default :
      FANS3D_ThermallyPerfect_KOmega_pState NEW(ZERO);
      NEW.spec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].c = ONE;
      return NEW;
   };
}

FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::lp(const int &index) const {
   double cc;
   switch(index){  
   case 1:
      cc = a_t();
      return (FANS3D_ThermallyPerfect_KOmega_pState(k/(THREE*cc*cc), 
                                                    -HALF*rho/cc, ZERO, ZERO, 
                                                    HALF/(cc*cc), 
                                                    rho/(THREE*cc*cc), ZERO, ZERO));
   case 2:
      cc = a_t();
      return (FANS3D_ThermallyPerfect_KOmega_pState(ONE - TWO*k/(THREE*cc*cc), 
                                                    ZERO, ZERO, ZERO, 
                                                    -ONE/(cc*cc), 
                                                    -TWO*rho/(THREE*cc*cc), ZERO, ZERO));
   case 3 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(ZERO, 
                                                    ZERO, ONE, ZERO, 
                                                    ZERO, 
                                                    ZERO, ZERO, ZERO));
   case 4 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(ZERO, 
                                                    ZERO, ZERO, ONE, 
                                                    ZERO, 
                                                    ZERO, ZERO, ZERO));
   case 5 :
      cc = a_t();
      return (FANS3D_ThermallyPerfect_KOmega_pState(k/(THREE*cc*cc), 
                                                    HALF*rho/cc, ZERO, ZERO, 
                                                    HALF/(cc*cc), 
                                                    rho/(THREE*cc*cc), ZERO, ZERO));
   case 6 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(ZERO, 
                                                    ZERO, ZERO, ZERO, 
                                                    ZERO, 
                                                    ONE, ZERO, ZERO));
   case 7 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(ZERO, 
                                                    ZERO, ZERO, ZERO, 
                                                    ZERO, 
                                                    ZERO, ONE, ZERO));
   default :
      FANS3D_ThermallyPerfect_KOmega_pState NEW(ZERO);
      NEW.spec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].c = ONE;
      return NEW;
   };
}

/*********************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::lp_x -- Primitive left eigenvector (x-direction).  *
 *********************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::lp_x(const int &index) {
   double cc;
   switch(index){  
   case 1:
      cc = a_t();
      return (FANS3D_ThermallyPerfect_KOmega_pState(k/(THREE*cc*cc), 
                                                    -HALF*rho/cc, ZERO, ZERO, 
                                                    HALF/(cc*cc), 
                                                    rho/(THREE*cc*cc), ZERO, ZERO));
   case 2:
      cc = a_t();
      return (FANS3D_ThermallyPerfect_KOmega_pState(ONE - TWO*k/(THREE*cc*cc), 
                                                    ZERO, ZERO, ZERO, 
                                                    -ONE/(cc*cc), 
                                                    -TWO*rho/(THREE*cc*cc), ZERO, ZERO));
   case 3 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(ZERO, 
                                                    ZERO, ONE, ZERO, 
                                                    ZERO, 
                                                    ZERO, ZERO, ZERO));
   case 4 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(ZERO, 
                                                    ZERO, ZERO, ONE, 
                                                    ZERO, 
                                                    ZERO, ZERO, ZERO));
   case 5 :
      cc = a_t();
      return (FANS3D_ThermallyPerfect_KOmega_pState(k/(THREE*cc*cc), 
                                                    HALF*rho/cc, ZERO, ZERO, 
                                                    HALF/(cc*cc), 
                                                    rho/(THREE*cc*cc), ZERO, ZERO));
   case 6 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(ZERO, 
                                                    ZERO, ZERO, ZERO, 
                                                    ZERO, 
                                                    ONE, ZERO, ZERO));
   case 7 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(ZERO, 
                                                    ZERO, ZERO, ZERO, 
                                                    ZERO, 
                                                    ZERO, ONE, ZERO));
   default :
      FANS3D_ThermallyPerfect_KOmega_pState NEW(ZERO);
      NEW.spec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].c = ONE;
      return NEW;
   };
}

FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::lp_x(const int &index) const {
   double cc;
   switch(index){  
   case 1:
      cc = a_t();
      return (FANS3D_ThermallyPerfect_KOmega_pState(k/(THREE*cc*cc), 
                                                    -HALF*rho/cc, ZERO, ZERO, 
                                                    HALF/(cc*cc), 
                                                    rho/(THREE*cc*cc), ZERO, ZERO));
   case 2:
     cc = a_t();
      return (FANS3D_ThermallyPerfect_KOmega_pState(ONE - TWO*k/(THREE*cc*cc), 
                                                    ZERO, ZERO, ZERO, 
                                                    -ONE/(cc*cc), 
                                                    -TWO*rho/(THREE*cc*cc), ZERO, ZERO));
   case 3 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(ZERO, 
                                                    ZERO, ONE, ZERO, 
                                                    ZERO, 
                                                    ZERO, ZERO, ZERO));
   case 4 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(ZERO, 
                                                    ZERO, ZERO, ONE, 
                                                    ZERO, 
                                                    ZERO, ZERO, ZERO));
   case 5 :
      cc = a_t();
      return (FANS3D_ThermallyPerfect_KOmega_pState(k/(THREE*cc*cc), 
                                                    HALF*rho/cc, ZERO, ZERO, 
                                                    HALF/(cc*cc), 
                                                    rho/(THREE*cc*cc), ZERO, ZERO));
   case 6 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(ZERO, 
                                                    ZERO, ZERO, ZERO, 
                                                    ZERO, 
                                                    ONE, ZERO, ZERO));
   case 7 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(ZERO, 
                                                    ZERO, ZERO, ZERO, 
                                                    ZERO, 
                                                    ZERO, ONE, ZERO));
   default :
      FANS3D_ThermallyPerfect_KOmega_pState NEW(ZERO);
      NEW.spec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].c = ONE;
      return NEW;
   };
}

/************************************************************************
 *************** NUMERICAL FLUX FUNCTIONS *******************************
 ************************************************************************/

/************************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::RoeAverage -- Roe-averaged primitive solution state.  *
 ************************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::
RoeAverage(const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
           const FANS3D_ThermallyPerfect_KOmega_pState &Wr) {

   double Hl, Hr, srhol, srhor;
   double Ha, ha;
   FANS3D_ThermallyPerfect_KOmega_pState Temp;

   /* Determine the left and right state specific enthalpies
      and square roots of the density. */

   Hl = Wl.Hs()/Wl.rho;
   Hr = Wr.Hs()/Wr.rho;
   srhol = sqrt(Wl.rho);
   srhor = sqrt(Wr.rho);
   
   /* Determine the appropriate Roe averages. */

   Temp.rho = srhol*srhor;
   Temp.v.x = (srhol*Wl.v.x+srhor*Wr.v.x)/(srhol+srhor);
   Temp.v.y = (srhol*Wl.v.y+srhor*Wr.v.y)/(srhol+srhor);
   Temp.v.z = (srhol*Wl.v.z+srhor*Wr.v.z)/(srhol+srhor);
   Temp.k = (srhol*Wl.k+srhor*Wr.k)/(srhol+srhor);
   Temp.omega = (srhol*Wl.omega+srhor*Wr.omega)/(srhol+srhor);

   for (int i=0; i<Wl.ns; i++) {
      Temp.spec[i].c = (srhol*Wl.spec[i].c + srhor*Wr.spec[i].c)/(srhol+srhor);
   } /* endif */

   Ha = (srhol*Hl+srhor*Hr)/(srhol+srhor);
   ha = Ha - HALF*(sqr(Temp.v.x)+sqr(Temp.v.y)+sqr(Temp.v.z)) - Temp.k;
   Temp.p = Temp.rho*Temp.T(ha)*Temp.Rtot();
   
   /* Return the Roe-averged state. */

   return (Temp);     

}

/************************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::FluxHLLE_x -- HLLE flux function, x-direction flux.   *
 ************************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::
FluxHLLE_x(const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
           const FANS3D_ThermallyPerfect_KOmega_pState &Wr) {
   
   double wavespeed_l, wavespeed_r;
 
   FANS3D_ThermallyPerfect_KOmega_pState Wa, lambdas_l, lambdas_r, lambdas_a;
   FANS3D_ThermallyPerfect_KOmega_cState Flux, dUrl;
   
   int num_vars = Wl.num_vars;

   /* Evaluate the Roe-average primitive solution state. */
   
   Wa = RoeAverage(Wl, Wr);
   
   /* Evaluate the jumps in the conserved solution states. */
 
   dUrl = Wr.U() - Wl.U();
 
   /* Evaluate the left, right, and average state eigenvalues. */
    
   lambdas_l = Wl.lambda_x();
   lambdas_r = Wr.lambda_x();
   lambdas_a = Wa.lambda_x();

   /* Determine the intermediate state flux. */
   
   wavespeed_l = min(lambdas_l[1],
                     lambdas_a[1]);
   wavespeed_r = max(lambdas_r[num_vars-NUM_FANS3D_VAR_EXTRA-lambdas_r.ns],
                     lambdas_a[num_vars-NUM_FANS3D_VAR_EXTRA-lambdas_a.ns]);
   
   wavespeed_l = min(wavespeed_l, ZERO);
   wavespeed_r = max(wavespeed_r, ZERO); 

   if (wavespeed_l >= ZERO) {
      Flux = Wl.F();
   } else if (wavespeed_r <= ZERO) {
      Flux = Wr.F();
   } else {
      Flux = ( (wavespeed_r*Wl.F()-wavespeed_l*Wr.F())
              +(wavespeed_l*wavespeed_r)*dUrl)/
              (wavespeed_r-wavespeed_l);
   } /* endif */
   
   /* Return solution flux. */
   
   return (Flux);
   
}

FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::
FluxHLLE_x(const FANS3D_ThermallyPerfect_KOmega_cState &Ul,
           const FANS3D_ThermallyPerfect_KOmega_cState &Ur) {

   return (FluxHLLE_x(Ul.W(), Ur.W()));

}

/************************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::FluxHLLE_n -- HLLE flux function, n-direction flux.   *
 ************************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::
FluxHLLE_n(const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
           const FANS3D_ThermallyPerfect_KOmega_pState &Wr,
           const Vector3D &norm_dir) {

   // Determine the left and right solution states in the rotate frame.
   FANS3D_ThermallyPerfect_KOmega_pState Wl_rot(Wl.Rotate(norm_dir));
   FANS3D_ThermallyPerfect_KOmega_pState Wr_rot(Wr.Rotate(norm_dir));
    
   // Evaluate the intermediate state solution flux in the rotated frame.
   FANS3D_ThermallyPerfect_KOmega_cState Flux_rot = FluxHLLE_x(Wl_rot, Wr_rot);
 
   // Return numerical flux in un-rotated frame.
   return (Flux_rot.RotateBack(norm_dir));

}

FANS3D_ThermallyPerfect_KOmega_cState  FANS3D_ThermallyPerfect_KOmega_pState::
FluxHLLE_n(const FANS3D_ThermallyPerfect_KOmega_cState &Ul,
           const FANS3D_ThermallyPerfect_KOmega_cState &Ur,
           const Vector3D &norm_dir) {

    return (FluxHLLE_n(Ul.W(), Ur.W(), norm_dir));

}

/************************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::FluxRoe_x -- Roe flux function, x-direction flux.     *
 ************************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::
FluxRoe_x(const  FANS3D_ThermallyPerfect_KOmega_pState &Wl,  
          const  FANS3D_ThermallyPerfect_KOmega_pState &Wr) {
   
   FANS3D_ThermallyPerfect_KOmega_pState Wa, dWrl, wavespeeds, 
                                         lambdas_l, lambdas_r, lambdas_a;
   FANS3D_ThermallyPerfect_KOmega_cState Flux;

   /* Evaluate the Roe-average primitive solution state. */ 
   
   Wa = RoeAverage(Wl, Wr);
 
   /* Evaluate the jumps in the primitive solution states. */

   dWrl = Wr-Wl;
   
   /* Evaluate the left, right, and average state eigenvalues. */

   lambdas_l = Wl.lambda_x();
   lambdas_r = Wr.lambda_x();
   lambdas_a = Wa.lambda_x();

   /* Determine the intermediate state flux. */

   if (Wa.v.x >= ZERO) {
      Flux = Wl.F();   
      wavespeeds = Wa.lambda_minus(lambdas_a,
                                   lambdas_l,
                                   lambdas_r);
      for (int i=1 ; i <= Wl.num_vars; i++) {
         if (wavespeeds[i] < ZERO) {
            Flux += wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
         } /* endif */
      } /* endfor */

   } else {
      Flux = Wr.F();
      wavespeeds = Wa.lambda_plus(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
      for (int i=1; i <= Wl.num_vars; i++) {
         if (wavespeeds[i] > ZERO) {
            Flux -= wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i); 
         } /* endif */
      } /* endfor */
   } /* endif */
    
   /* Return solution flux. */    

   return (Flux);    
   
}
   
/************************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::FluxRoe_n -- Roe flux function, n-direction flux.     *
 ************************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::
FluxRoe_n(const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
          const FANS3D_ThermallyPerfect_KOmega_pState &Wr,
          const Vector3D &norm_dir) {
   
   // Determine the left and right solution states in the rotate frame.
   FANS3D_ThermallyPerfect_KOmega_pState Wl_rot(Wl.Rotate(norm_dir));
   FANS3D_ThermallyPerfect_KOmega_pState Wr_rot(Wr.Rotate(norm_dir));
    
   // Evaluate the intermediate state solution flux in the rotated frame.
   FANS3D_ThermallyPerfect_KOmega_cState Flux_rot = FluxRoe_x(Wl_rot, Wr_rot);
 
   // Return numerical flux in un-rotated frame.
   return (Flux_rot.RotateBack(norm_dir));

}

/************************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::lambda_minus -- Negative wave speeds determined using *
 *                                                        Harten entropy fix.                   *
 ************************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::
lambda_minus(const FANS3D_ThermallyPerfect_KOmega_pState &lambdas_a,
             const FANS3D_ThermallyPerfect_KOmega_pState &lambdas_l,
             const FANS3D_ThermallyPerfect_KOmega_pState &lambdas_r) {
  
   FANS3D_ThermallyPerfect_KOmega_pState W_temp;

   W_temp.rho = HartenFixNeg(lambdas_a[1], lambdas_l[1], lambdas_r[1]);
   W_temp.v.x = HALF*(lambdas_a[2] - fabs(lambdas_a[2]));
   W_temp.v.y = HALF*(lambdas_a[3] - fabs(lambdas_a[3]));
   W_temp.v.z = HALF*(lambdas_a[4] - fabs(lambdas_a[4]));
   W_temp.p = HartenFixNeg(lambdas_a[5], lambdas_l[5], lambdas_r[5]);
   W_temp.k = HALF*(lambdas_a[6] - fabs(lambdas_a[6]));
   W_temp.omega = HALF*(lambdas_a[7] - fabs(lambdas_a[7]));
   
   for (int i = (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA + 1); i <= W_temp.num_vars; i++) {
      W_temp.spec[i-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].c = 
         HALF*(lambdas_a[i] - fabs(lambdas_a[i]));
   } /* endfor */
   
   return (W_temp);

}

/************************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::lambda_plus -- Positive wave speeds determined using  *
 *                                                       Harten entropy fix.                    *
 ************************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::
lambda_plus(const FANS3D_ThermallyPerfect_KOmega_pState &lambdas_a,
            const FANS3D_ThermallyPerfect_KOmega_pState &lambdas_l,
            const FANS3D_ThermallyPerfect_KOmega_pState &lambdas_r) {
   
   FANS3D_ThermallyPerfect_KOmega_pState W_temp;

   W_temp.rho = HartenFixPos(lambdas_a[1], lambdas_l[1], lambdas_r[1]);
   W_temp.v.x = HALF*(lambdas_a[2] + fabs(lambdas_a[2]));
   W_temp.v.y = HALF*(lambdas_a[3] + fabs(lambdas_a[3]));
   W_temp.v.z = HALF*(lambdas_a[4] + fabs(lambdas_a[4]));
   W_temp.p = HartenFixPos(lambdas_a[5], lambdas_l[5], lambdas_r[5]);
   W_temp.k = HALF*(lambdas_a[6] + fabs(lambdas_a[6]));
   W_temp.omega = HALF*(lambdas_a[7] + fabs(lambdas_a[7]));
  
   for (int i = (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA +1); i <= W_temp.num_vars; i++) {
      W_temp.spec[i-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].c = 
         HALF*(lambdas_a[i] + fabs(lambdas_a[i]));
   } /* endfor */
   
   return (W_temp);

}

/************************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::HLLE_wavespeeds -- Returns the lambda plus and lambda *
 *                    minus wave speeds for rotated Riemann problem aligned with norm_dir       *    
 *                    given unroated solution states Wl and Wr.                                 *
 ************************************************************************************************/
Vector2D FANS3D_ThermallyPerfect_KOmega_pState::HLLE_wavespeeds(const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
                                                                const FANS3D_ThermallyPerfect_KOmega_pState &Wr,
                                                                const Vector3D &norm_dir) {

    Vector2D wavespeed;
    FANS3D_ThermallyPerfect_KOmega_pState Wa_n, lambdas_l, lambdas_r, lambdas_a;  //Lots of TEMPS
    FANS3D_ThermallyPerfect_KOmega_pState Wl_rotated(Wl.Rotate(norm_dir));
    FANS3D_ThermallyPerfect_KOmega_pState Wr_rotated(Wr.Rotate(norm_dir));

    /* Evaluate the Roe-average primitive solution state. */                           
    Wa_n = Wa_n.RoeAverage(Wl_rotated, Wr_rotated);
    
    /* Evaluate the left, right, and average state eigenvalues. */
    lambdas_l = Wl_rotated.lambda_x();
    lambdas_r = Wr_rotated.lambda_x();
    lambdas_a = Wa_n.lambda_x();

    /* Determine the intermediate state flux. */
    wavespeed.x = min(lambdas_l[1],
                      lambdas_a[1]);
    wavespeed.y = max(lambdas_r[lambdas_r.NumVarSansSpecies()],
                      lambdas_a[lambdas_a.NumVarSansSpecies()]);
 
    wavespeed.x = min(wavespeed.x, ZERO); //lambda minus
    wavespeed.y = max(wavespeed.y, ZERO); //lambda plus 

    return (wavespeed);

}

/************************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::Rotate -- Returns a rotated primitive state aligned   *
 *                                            with a local x-axis in the norm_dir.              *
 ************************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::Rotate(const Vector3D &norm_dir) const {

  // for a 3D unit normal rotated to align with the x-axis
  double Ct = norm_dir.x;  //cos_angle
  double St = sqrt( norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z); //sin_angle
  Vector3D rt(0,norm_dir.z,-norm_dir.y);  //rotation axis  

  return FANS3D_ThermallyPerfect_KOmega_pState(rho,
                                               v.x*Ct - v.y*rt.z*St + v.z*rt.y*St,
                                               v.x*rt.z*St + v.y*(rt.y*rt.y*(ONE-Ct)+Ct) + v.z*(rt.y*rt.z*(ONE-Ct)),
                                               -v.x*rt.y*St +  v.y*(rt.y*rt.z*(ONE-Ct)) + v.z*(rt.z*rt.z*(ONE-Ct)+Ct),
                                               p,
                                               k,
                                               omega,
                                               spec);

}

/************************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::Rotate -- Returns an un-rotated primitive state       *
 *                                            re-alinged from the x-axis of the global          *
 *                                            problem.                                          *
 ************************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::RotateBack(const Vector3D &norm_dir) const {

  // for a 3D unit normal rotated to align with the x-axis
  double Ct = norm_dir.x;  //cos_angle
  double St = sqrt( norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z); //sin_angle
  Vector3D rt(0,norm_dir.z,-norm_dir.y);  //rotation axis  

  return FANS3D_ThermallyPerfect_KOmega_pState(rho,
                                               v.x*Ct + v.y*rt.z*St - v.z*rt.y*St,
                                               -v.x*rt.z*St + v.y*(rt.y*rt.y*(ONE-Ct)+Ct) + v.z*(rt.y*rt.z*(ONE-Ct)),
                                               + v.x*rt.y*St +  v.y*(rt.y*rt.z*(ONE-Ct)) + v.z*(rt.z*rt.z*(ONE-Ct)+Ct),
                                               p,
                                               k,
                                               omega,
                                               spec);

}

/************************************************************************
 *************** NUMERICAL EVALUATION OF VISCOUS FLUXES *****************
 ************************************************************************/

/********************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::FluxViscous_n -- Viscous flux (n-direction).      *
/********************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::
FluxViscous_n(const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
              const FANS3D_ThermallyPerfect_KOmega_pState &Wr,
              const FANS3D_ThermallyPerfect_KOmega_pState &Wc,
              const FANS3D_ThermallyPerfect_KOmega_pState &Wc_Neigbor,
              const FANS3D_ThermallyPerfect_KOmega_pState &dWdx,
              const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
              const FANS3D_ThermallyPerfect_KOmega_pState &dWdz,
              const FANS3D_ThermallyPerfect_KOmega_pState &dWdx_Neigbor,
              const FANS3D_ThermallyPerfect_KOmega_pState &dWdy_Neigbor,
              const FANS3D_ThermallyPerfect_KOmega_pState &dWdz_Neigbor,
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
   
   FANS3D_ThermallyPerfect_KOmega_pState dWdx_Weighted, 
      dWdy_Weighted, dWdz_Weighted, dWdx_face, 
      dWdy_face, dWdz_face, Grad_middle_term;

   FANS3D_ThermallyPerfect_KOmega_pState W_face;
   
   dWdx_Weighted = alpha*dWdx + (1.0 - alpha)*dWdx_Neigbor;
   dWdy_Weighted = alpha*dWdy + (1.0 - alpha)*dWdy_Neigbor;
   dWdz_Weighted = alpha*dWdz + (1.0 - alpha)*dWdz_Neigbor;

   //dWdx_face  = dWdx_Weighted;
   //dWdy_face  = dWdy_Weighted;
   //dWdz_face  = dWdz_Weighted; 
   
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
   //W_face = HALF*(Wc + Wc_Neigbor);

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

/************************************************************************
 ******************** BOUNDARY CONDITIONS *******************************
 ************************************************************************/

/************************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::Reflect -- Return reflected solution state.           *
 ************************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::
Reflect(const FANS3D_ThermallyPerfect_KOmega_pState &W,
        const Vector3D &norm_dir) {
   
   Vector3D ur_norm, ur_tang, vr_tot;
   FANS3D_ThermallyPerfect_KOmega_pState Temp; Temp.Copy(W);

   ur_norm = dot(W.v, norm_dir)*norm_dir;
   ur_tang = W.v - ur_norm;

   ur_norm = -ur_norm;
   vr_tot = ur_norm + ur_tang;

   Temp.v = vr_tot;
   
   return (Temp);
       
}

/************************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::MovingWall -- Return moving wall boundary state.      *
 ************************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::
MovingWall(const FANS3D_ThermallyPerfect_KOmega_pState &Win,
           const FANS3D_ThermallyPerfect_KOmega_pState &Wout,
           const Vector3D &norm_dir, 
           const Vector3D &wall_velocity,
           const Vector3D &pressure_gradient,
           const int &TEMPERATURE_BC_FLAG) {
   
   FANS3D_ThermallyPerfect_KOmega_pState Temp; Temp.Copy(Win);
   
   if (wall_velocity == Vector3D_ZERO) {
     Temp.v = -Win.v;

   } else {
     double  Wall_velocity_tang ;
     Vector3D ur_norm, ur_tang, vr_tot, uw_tang;
     ur_norm = dot(Win.v, norm_dir)*norm_dir;
     ur_tang = Win.v - ur_norm;
     
     uw_tang = wall_velocity - dot(norm_dir,wall_velocity)*norm_dir;
     
     ur_norm = -ur_norm;
     ur_tang = 2.0*uw_tang - ur_tang;
     vr_tot = ur_norm + ur_tang;
     
     Temp.v = vr_tot;
   } /* endif */
  
  Temp.k = ZERO; 

  /* Fixed Wall Temperature or constant extrapolation for Adiabatic */

  if (TEMPERATURE_BC_FLAG == FIXED_TEMPERATURE_WALL) {
     if (pressure_gradient != Vector3D_ZERO) {
        Temp.rho = Wout.p/(Temp.Rtot()*Wout.T());
     } else {
        Temp.rho = Temp.p/(Temp.Rtot()*Wout.T());
     } /* endif */
  } /* endif */

  return (Temp);

}

/************************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::NoSlip -- Return no-slip wall boundary state.         *
 ************************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::
NoSlip(const FANS3D_ThermallyPerfect_KOmega_pState &Win,
       const FANS3D_ThermallyPerfect_KOmega_pState &Wout,
       const Vector3D &norm_dir,
       const Vector3D &pressure_gradient,
       const int &TEMPERATURE_BC_FLAG) {
   
   return (MovingWall(Win, 
                      Wout, 
                      norm_dir, 
                      Vector3D_ZERO, 
                      pressure_gradient, 
                      TEMPERATURE_BC_FLAG));
   
}

/************************************************************************
 *************************** SOURCE TERMS *******************************
 ************************************************************************/

/***************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::Sturbulene -- Turbulence model source terms. *
 ***************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::
Sturbulence(FANS3D_ThermallyPerfect_KOmega_pState &Wc,
            const FANS3D_ThermallyPerfect_KOmega_pState &dWdx,
            const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
            const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {
   
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   double mu_t_temp, production;
   Tensor3D Reynolds_stress, strain_rate_tensor, vorticity_tensor;
   
   mu_t_temp = Wc.mu_t();
   Reynolds_stress = Wc.tau_t(mu_t_temp, dWdx, dWdy, dWdz);
  
   production = Reynolds_stress.xx*dWdx.v.x + 
                Reynolds_stress.xy*(dWdy.v.x + dWdx.v.y) +
                Reynolds_stress.yy*dWdy.v.y +
                Reynolds_stress.xz*(dWdz.v.x + dWdx.v.z) +
                Reynolds_stress.yz*(dWdz.v.y + dWdy.v.z) +
                Reynolds_stress.zz*dWdz.v.z;
   
   vorticity_tensor = Wc.vorticity(dWdx, dWdy, dWdz);
   strain_rate_tensor = Wc.strain_rate(dWdx, dWdy, dWdz);
   
   Temp.Vacuum();
   Temp.rhok = production - 
               Wc.k_omega_model.beta_star*
               Wc.k_omega_model.f_betastar(dWdx.k, dWdy.k, dWdz.k, 
                                           dWdx.omega, dWdy.omega, dWdz.omega, Wc.omega)*
               Wc.rho*Wc.k*Wc.omega;

   Temp.rhoomega = Wc.k_omega_model.alpha*(Wc.omega/max(Wc.k, TOLER))*production - 
                   Wc.k_omega_model.beta*
                   Wc.k_omega_model.f_beta(vorticity_tensor, strain_rate_tensor, Wc.omega)*
                   Wc.rho*Wc.omega*Wc.omega;
   
   return (Temp);

} 

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState member functions                                *
 *****************************************************************************************/

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState::Realizable_Solution_Check -- Check physical    *
 *                                                                      validity of      *
 *                                                                      solution state.  *
 *****************************************************************************************/
bool FANS3D_ThermallyPerfect_KOmega_cState::Realizable_Solution_Check(void) {
   if (rho <= ZERO || !negative_speccheck() || es() <= ZERO ||
       rhok <= ZERO || rhoomega <= ZERO) {    
      cout << "\n " << CFFC_Name() 
           << " ERROR: Conserved solution state has a negative density, energy, mass fractions,"
           << " turbulent kinetic energy, and/or dissipation rate.\n";
      return false;
   } else {
      return true;
   } /* endif */
} 

/**************************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState::e -- Return mixture absolute internal energy.           *
 **************************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_cState::e(void) const {
  return ((E - HALF*rhov.sqr()/rho-rhok)/rho);
}

/**************************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState::es -- Return sensible internal energy.                  *
 **************************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_cState::es(void) const {
  return ((E - HALF*rhov.sqr()/rho-rhok)/rho-HeatofFormation());
}

/**************************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState::h -- Return mixture absolute internal enthalpy.         *
 **************************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_cState::h(void) const {
  return (e()+p()/rho);
}

double FANS3D_ThermallyPerfect_KOmega_cState::h(const double &Temp) const {
  return (Euler3D_ThermallyPerfect_cState::h(Temp));
}

/**************************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState::hs -- Return sensible internal enthalpy.                *
 **************************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_cState::hs(void) const {
  return (es()+p()/rho);
}

double FANS3D_ThermallyPerfect_KOmega_cState::hs(const double &Temp) const {
  return (Euler3D_ThermallyPerfect_cState::hs());
}

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState::T -- Return mixture temperature.               *
 *****************************************************************************************/
/********************* Temperature ***********************************
 This gets complicated when trying to obtain the primitive variables 
  from the conserved as E is a function of Temperature taken from 
  a nonlinear equation (ie. a polynomial).   Thus is can't be rearranged
  to find T and p, so a iterative newtons method has to used. Simple 
  enough, but can be extremely expensive in terms of computation time.
  Oh well ce la vie...

  Also note the "E" is actually rho*E is. E = (rho *(e + HALF*v^2+k))
  in the conserved state.

  If more that 20 iterations are taken than the look jumps out 
  and gives a warning but will continue.  On average it almost always 
  converges in less than 5 iterations with "tolerance" which is defined
  in the header.
**********************************************************************/
double FANS3D_ThermallyPerfect_KOmega_cState::T(void) const {

   double T = ZERO;
   double RTOT = Rtot();  
   //--------- Initial Guess ------------------------------//
   //using a polytropic gas assumption with gamma@200;
   double Tguess = (gamma_guess() - ONE)*(E - HALF*rhov.sqr()/rho-rhok)/(rho*RTOT);
   //--------- global newtons method to get T ---------------//
   double A = (E - HALF*rhov*rhov/rho -rhok)/rho;
 
   int numit =0;
   double Tmin = low_temp_range;
   double Tmax = high_temp_range;
   
   //check for start value
   if (Tguess > Tmin && Tguess < Tmax) {
      T=Tguess;
   } else {
      T=Tmin;
   } /* endif */
   
   double fa = h(Tmin) - Tmin*RTOT - A;
   double fn = h(T) - T*RTOT - A;
   double dfn = hprime(T) - RTOT;
   while (fabs(Tmax-Tmin) > CONV_TOLERANCE && 
          fabs(fn) > CONV_TOLERANCE && 
          numit < 20 && 
          T >= low_temp_range) {    
      // Newton 
      if (T >= Tmin && T <= Tmax) {
         T = T - fn/dfn;
         if(T >= Tmax) T = HALF*(Tmax - Tmin);	
         //Bisection
      } else {
         T = HALF*(Tmax - Tmin);
      } /* endif */
      //evaluate function and derivative
      fn = h(T) - T*RTOT - A;
      dfn = hprime(T) - RTOT;  
      //change bisection range
      if ( fa*fn <=ZERO) {
         Tmax = T;
      } else {
         Tmin = T;
         fa = fn;
      } /* endif */
      numit++;
   } /* endwhile */

   if (numit>=19 || T <= low_temp_range) {
      T = max(Tguess,low_temp_range); 	
      if (debug_level) { 
         cout << "\nTemperature didn't converge in FANS3D_ThermallyPerfect_KOmega_cState::T(void)";
         cout << " with polytopic Tguess " << Tguess << ", or lower than Tmin " 
              << low_temp_range << " using " << T;
      } /* endif */
   } /* endif */

   return T;

} 

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState::a_t -- Return mixture sound speed (including   *
 *                                               turbulent kinetic energy).              *
 *****************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_cState::a_t(void) const {
   double aa = sqr(a());
   aa += (TWO/THREE)*(rhok/rho)*g();
   return sqrt(aa);
}

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState::k -- Return turbulent kinetic energy.          *
 *****************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_cState::k() const {
   return (rhok/rho);
}

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState::omega -- Return specific dissipation rate for  *
 *                                                 turbulent energy.                     *
 *****************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_cState::omega() const {
   return (rhoomega/rho);
}

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState::mu_t -- Return eddy (turbulent) viscosity.     *
 *****************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_cState::mu_t(void) {
   return (rho*rhok/max(TOLER,rhoomega));
}

/********************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState::kappa_t -- Return turbulent thermal conductivity. *
 ********************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_cState::kappa_t(void) {
  return (mu_t()*Cp()/Pr_t()); 
}

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState::Ds_t -- Return species turbulent diffusion     *
 *                                                coefficient.                           *
 *****************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_cState::Ds_t(const int &i) {
   return (mu_t()/(rho*Sc_t()));
}

double FANS3D_ThermallyPerfect_KOmega_cState::Ds_t(const int &i,
                                                   const double &mu_t_temp) {
   return (mu_t_temp/(rho*Sc_t()));
}

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState::Pr_t -- Return turbulent Prandtl number.       *
 *****************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_cState::Pr_t(void) {
   return (0.9);
}

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState::Sc_t -- Return turbulent Schmidt number.       *
 *****************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_cState::Sc_t(void) {
   return (1.0);
}

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState::Le_t -- Return turbulent Lewis number.         *
 *****************************************************************************************/
double FANS3D_ThermallyPerfect_KOmega_cState::Le_t(void) {
  return (Sc_t()/Pr_t());
}

/*****************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState::W -- Return primitive solution state vector.   *
 *****************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_cState::W(void) {
   FANS3D_ThermallyPerfect_KOmega_pState Temp;
   Temp.rho = rho;
   Temp.v = v();  
   Temp.p = p();
   Temp.k = k();
   Temp.omega = omega();
   for (int i = 0; i < ns; i++) {
      Temp.spec[i] = rhospec[i]/rho;
   } /* endfor */
   return Temp;
}

FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_cState::W(void) const {
   FANS3D_ThermallyPerfect_KOmega_pState Temp;
   Temp.rho = rho;
   Temp.v = v();  
   Temp.p = p();
   Temp.k = k();  
   Temp.omega = omega();
   for (int i=0; i<ns; i++) {
      Temp.spec[i] = rhospec[i]/rho;
   } /* endfor */
   return Temp;
}

FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_cState::
W(const FANS3D_ThermallyPerfect_KOmega_cState &U) const {
   FANS3D_ThermallyPerfect_KOmega_pState Temp;
   Temp.rho = U.rho;
   Temp.v = U.v();  
   Temp.p = U.p();
   Temp.k = U.k();  
   Temp.omega = U.omega();
   for (int i = 0; i < U.ns; i++) {
      Temp.spec[i] = U.rhospec[i]/U.rho;
   } /* endfor */
   return Temp;
}

/************************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState::Rotate -- Returns a rotated primitive state aligned   *
 *                                            with a local x-axis in the norm_dir.              *
 ************************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_cState::Rotate(const Vector3D &norm_dir) const {

  // for a 3D unit normal rotated to align with the x-axis
  double Ct = norm_dir.x;  //cos_angle
  double St = sqrt( norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z); //sin_angle
  Vector3D rt(0,norm_dir.z,-norm_dir.y);  //rotation axis  

  return FANS3D_ThermallyPerfect_KOmega_cState(rho,
                                               rhov.x*Ct - rhov.y*rt.z*St + rhov.z*rt.y*St,
                                               rhov.x*rt.z*St + rhov.y*(rt.y*rt.y*(ONE-Ct)+Ct) + rhov.z*(rt.y*rt.z*(ONE-Ct)),
                                               -rhov.x*rt.y*St +  rhov.y*(rt.y*rt.z*(ONE-Ct)) + rhov.z*(rt.z*rt.z*(ONE-Ct)+Ct),
                                               E,
                                               rhok,
                                               rhoomega,
                                               rhospec);

}

/************************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState::Rotate -- Returns an un-rotated primitive state       *
 *                                            re-alinged from the x-axis of the global          *
 *                                            problem.                                          *
 ************************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_cState::RotateBack(const Vector3D &norm_dir) const {

  // for a 3D unit normal rotated to align with the x-axis
  double Ct = norm_dir.x;  //cos_angle
  double St = sqrt( norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z); //sin_angle
  Vector3D rt(0,norm_dir.z,-norm_dir.y);  //rotation axis  

  return FANS3D_ThermallyPerfect_KOmega_cState(rho,
                                               rhov.x*Ct + rhov.y*rt.z*St - rhov.z*rt.y*St,
                                               -rhov.x*rt.z*St + rhov.y*(rt.y*rt.y*(ONE-Ct)+Ct) + rhov.z*(rt.y*rt.z*(ONE-Ct)),
                                               + rhov.x*rt.y*St +  rhov.y*(rt.y*rt.z*(ONE-Ct)) + rhov.z*(rt.z*rt.z*(ONE-Ct)+Ct),
                                               E,
                                               rhok,
                                               rhoomega,
                                               rhospec);

}
