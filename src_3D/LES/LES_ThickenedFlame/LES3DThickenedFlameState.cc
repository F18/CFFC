/*! \file  LES3DThickenedFlameState.cc
\brief  Definition of member functions for 3D Favre-filtered Navier-Stokes solution classes
        associated with solution of premixed compressible turbulent combusting flows of a 
        thermally perfect using a large eddy simulation (LES) technique in conjunction with 
        thickened flame subfilter scale model.
*/

/* Include LES3DThickenedFlameState header file. */

#ifndef _LES3DTF_STATE_INCLUDED
#include "LES3DThickenedFlameState.h"
#endif // LES3DTF_STATE_INCLUDED   

/**************************************************************************************
 * LES3DTF_pState -- Create storage and assign various static values.                 *
 **************************************************************************************/
double LES3DTF_pState::Mref = 0.1;
double LES3DTF_pState::_laminar_flame_speed=0.3837;
double LES3DTF_pState::_laminar_flame_thickness=4.4E-04;
double LES3DTF_pState::_TFactor = 1.0;

/**************************************************************************************
 * LES3DTF_cState -- Create storage and assign various static values.                 *
 **************************************************************************************/
double LES3DTF_cState::Mref = 0.1;
double LES3DTF_cState::_laminar_flame_speed=0.3837;
double LES3DTF_cState::_laminar_flame_thickness=4.4E-04;
double LES3DTF_cState::_TFactor = 1.0;

/**************************************************************************************
 * LESS3DFsd_pState member functions                                                  *
 **************************************************************************************/

/**************************************************************************************
 * LES3DTF_pState::Copy -- Makes a copy of solution state vector.                     *
 **************************************************************************************/
void LES3DTF_pState::Copy(const LES3DTF_pState &W) {
  rho = W.rho;
  v = W.v; 
  p = W.p;  
  k = W.k; 
  for (int i=0; i<ns; i++) {
    spec[i] = W.spec[i];
  } /* endfor */
  flame = W.flame;
}

/********************************************************************************************
 * LES3DTF_pState::Realizable_Solution_Check -- Check physical validity of solution state.  *
 ********************************************************************************************/
bool LES3DTF_pState::Realizable_Solution_Check(void) {
   if (rho <= ZERO || !negative_speccheck() || p <= ZERO || k < ZERO ) {    
      cout << "\n " << CFFC_Name() 
           << " ERROR: Primitive solution state has a negative density, pressure, mass fractions,"
           << " and/or turbulent kinetic energy.\n";
      return false;
   } else {
      return true;
   } /* endif */
} 

/****************************************************************************************
 * LES3DTF_pState::E -- Return total energy of the mixture.                             *
 ****************************************************************************************/
double LES3DTF_pState::E(void) const {   
   return (rho*(e()+HALF*v.sqr()+k));   
}

/****************************************************************************************
 * LES3DTF_pState::H -- Return total enthalpy of the mixture.                           *
 ****************************************************************************************/
double LES3DTF_pState::H(void) const{
   return (rho*(h()+HALF*v.sqr()+FIVE*k/THREE));   
}

/****************************************************************************************
 * LES3DTF_pState::Hs -- Return total mixture sensible enthalpy.                        *
 ****************************************************************************************/
double LES3DTF_pState::Hs(void) const{
   return (rho*(hs()+HALF*v.sqr()+FIVE*k/THREE));   
}

/*****************************************************************************************
 * LES3DTF_pState::p_t -- Return turbulence modified pressure.                           *
 *****************************************************************************************/
double LES3DTF_pState::p_t(void) const { 
   return (p + 2.0*rho*k/3.0);
}

/*******************************************************************************************
 * LES3DTF_pState::a_t -- Return mixture sound speed (including turbulent kinetic energy). *
 *******************************************************************************************/
double LES3DTF_pState::a_t(void) {
   double aa = sqr(a());
   aa += (TWO/THREE)*k*g();
   return sqrt(aa);
}

double LES3DTF_pState::a_t(void) const {
   double aa = sqr(a());
   aa += (TWO/THREE)*k*g();
   return sqrt(aa);
}

/******************************************************************
 * LES3DTF_pState::mu_t -- Return eddy (turbulent) viscosity.     *
 ******************************************************************/
double LES3DTF_pState::mu_t(const LES3DTF_pState &dWdx,
 			     const LES3DTF_pState &dWdy,
			     const LES3DTF_pState &dWdz,
                             const int Flow_Type, 
                             const double &Volume) {
  double filter = filter_width(Volume);
  if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY) {
    double Cs = 0.18;
    return(rho*sqr(Cs*filter)*abs_strain_rate(dWdx,dWdy,dWdz));
  } else if(Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {
    double Cv = 0.086;
    return(rho*Cv*sqrt(k)*filter);
  }
}

/*************************************************************************
 * LES3DTF_pState::kappa_t -- Return turbulent thermal conductivity.     *
 *************************************************************************/
double LES3DTF_pState::kappa_t(const LES3DTF_pState &dWdx,
 			        const LES3DTF_pState &dWdy,
			        const LES3DTF_pState &dWdz,
                                const int Flow_Type, 
                                const double &Volume) {
  return (mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume)*Cp()/Pr_t()); 
}

/*******************************************************************************
 * LES3DTF_pState::Ds_t -- Return species turbulent diffusion coefficient.     *
 *******************************************************************************/
double LES3DTF_pState::Ds_t(const int i,
                             const LES3DTF_pState &dWdx,
 			     const LES3DTF_pState &dWdy,
			     const LES3DTF_pState &dWdz,
                             const int Flow_Type, 
                             const double &Volume) {
  return (mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume)/(rho*Sc_t())); 
}

double LES3DTF_pState::Ds_t(const int i,
                             const double &mu_t_temp) {
  return (mu_t_temp/(rho*Sc_t())); 
}

/******************************************************************
 * LES3DTF_pState::Pr_t -- Return turbulent Prandtl number.       *
 ******************************************************************/
double LES3DTF_pState::Pr_t(void) {
   return (0.9); 
}

/*******************************************************************
 * LES3DTF_pState::Sc_t -- Return turbulent Schmidt number.       *
 *******************************************************************/
double LES3DTF_pState::Sc_t(void) {
   return (1.0);
}

/******************************************************************
 * LES3DTF_pState::Le_t -- Return turbulent Lewis number.         *
 ******************************************************************/
double LES3DTF_pState::Le_t(void) {
  return (Sc_t()/Pr_t());
}

/************************************************************************
 * LES3DTF_pState::filter_width -- LES characteristic filter width      *
 ************************************************************************/
double LES3DTF_pState::filter_width(const double &Volume) const {
  return (2.0*pow(Volume,1.0/3.0)); 
}

/******************************************************************************
 * LES3DTF_pState::tau_t -- Return subfilter scale (turbulent) stress tensor. *
 ******************************************************************************/
Tensor3D LES3DTF_pState::tau_t(const LES3DTF_pState &dWdx, 
                                const LES3DTF_pState &dWdy,
                                const LES3DTF_pState &dWdz,
                                const int Flow_Type, 
                                const double &Volume) {
   
   Tensor3D SFS_stress;
   double mu_t_temp = mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume);

   SFS_stress.xx = mu_t_temp*(FOUR*dWdx.v.x - TWO*dWdy.v.y - TWO*dWdz.v.z)/THREE - 
                   (TWO/THREE)*rho*k;
   SFS_stress.xy = mu_t_temp*(dWdx.v.y + dWdy.v.x);
   SFS_stress.xz = mu_t_temp*(dWdx.v.z + dWdz.v.x);
   SFS_stress.yy = mu_t_temp*(FOUR*dWdy.v.y - TWO*dWdx.v.x - TWO*dWdz.v.z)/THREE - 
                   (TWO/THREE)*rho*k;
   SFS_stress.yz = mu_t_temp*(dWdy.v.z + dWdz.v.y);
   SFS_stress.zz = mu_t_temp*(FOUR*dWdz.v.z - TWO*dWdx.v.x - TWO*dWdy.v.y)/THREE - 
                   (TWO/THREE)*rho*k;

   return (SFS_stress);

}

Tensor3D LES3DTF_pState::tau_t(const double &mu_t_temp,
                                const LES3DTF_pState &dWdx, 
                                const LES3DTF_pState &dWdy,
                                const LES3DTF_pState &dWdz,
                                const int Flow_Type, 
                                const double &Volume) {
   
   Tensor3D SFS_stress;

   SFS_stress.xx = mu_t_temp*(FOUR*dWdx.v.x - TWO*dWdy.v.y - TWO*dWdz.v.z)/THREE - 
                   (TWO/THREE)*rho*k;
   SFS_stress.xy = mu_t_temp*(dWdx.v.y + dWdy.v.x);
   SFS_stress.xz = mu_t_temp*(dWdx.v.z + dWdz.v.x);
   SFS_stress.yy = mu_t_temp*(FOUR*dWdy.v.y - TWO*dWdx.v.x - TWO*dWdz.v.z)/THREE - 
                   (TWO/THREE)*rho*k;
   SFS_stress.yz = mu_t_temp*(dWdy.v.z + dWdz.v.y);
   SFS_stress.zz = mu_t_temp*(FOUR*dWdz.v.z - TWO*dWdx.v.x - TWO*dWdy.v.y)/THREE - 
                   (TWO/THREE)*rho*k;

   return (SFS_stress);

}

/**********************************************************************************************
 * LES3DTF_pState::tau_t_x -- Return components of subfilter scale (turbulent)               *
 *                             stress tensor in the x-direction.                              *
 **********************************************************************************************/
Tensor3D LES3DTF_pState::tau_t_x(const LES3DTF_pState &dWdx, 
                                  const LES3DTF_pState &dWdy,
                                  const LES3DTF_pState &dWdz,
                                  const int Flow_Type, 
                                  const double &Volume) {
   Tensor3D SFS_stress;
   double mu_t_temp = mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume);

   SFS_stress.xx = mu_t_temp*(FOUR*dWdx.v.x - TWO*dWdy.v.y - TWO*dWdz.v.z)/THREE - 
                   (TWO/THREE)*rho*k;
   SFS_stress.xy = mu_t_temp*(dWdx.v.y + dWdy.v.x);
   SFS_stress.xz = mu_t_temp*(dWdx.v.z + dWdz.v.x);

   return (SFS_stress);

}

Tensor3D LES3DTF_pState::tau_t_x(const double &mu_t_temp,
                                  const LES3DTF_pState &dWdx, 
                                  const LES3DTF_pState &dWdy,
                                  const LES3DTF_pState &dWdz,
                                  const int Flow_Type, 
                                  const double &Volume) {
   Tensor3D SFS_stress;

   SFS_stress.xx = mu_t_temp*(FOUR*dWdx.v.x - TWO*dWdy.v.y - TWO*dWdz.v.z)/THREE - 
                   (TWO/THREE)*rho*k;
   SFS_stress.xy = mu_t_temp*(dWdx.v.y + dWdy.v.x);
   SFS_stress.xz = mu_t_temp*(dWdx.v.z + dWdz.v.x);

   return (SFS_stress);

}

/**********************************************************************************************
 * LES3DTF_pState::tau_t_y -- Return components of subfilter scale (turbulent)               *
 *                             stress tensor in the y-direction.                              *
 **********************************************************************************************/
Tensor3D LES3DTF_pState::tau_t_y(const LES3DTF_pState &dWdx, 
                                  const LES3DTF_pState &dWdy,
                                  const LES3DTF_pState &dWdz,
                                  const int Flow_Type, 
                                  const double &Volume) {
   Tensor3D SFS_stress;
   double mu_t_temp = mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume);

   SFS_stress.xy = mu_t_temp*(dWdx.v.y + dWdy.v.x);
   SFS_stress.yy = mu_t_temp*(FOUR*dWdy.v.y - TWO*dWdx.v.x - TWO*dWdz.v.z)/THREE - 
                   (TWO/THREE)*rho*k;
   SFS_stress.yz = mu_t_temp*(dWdy.v.z + dWdz.v.y);

   return (SFS_stress);

}

Tensor3D LES3DTF_pState::tau_t_y(const double &mu_t_temp,
                                  const LES3DTF_pState &dWdx, 
                                  const LES3DTF_pState &dWdy,
                                  const LES3DTF_pState &dWdz,
                                  const int Flow_Type, 
                                  const double &Volume) {
   Tensor3D SFS_stress;

   SFS_stress.xy = mu_t_temp*(dWdx.v.y + dWdy.v.x);
   SFS_stress.yy = mu_t_temp*(FOUR*dWdy.v.y - TWO*dWdx.v.x - TWO*dWdz.v.z)/THREE - 
                   (TWO/THREE)*rho*k;
   SFS_stress.yz = mu_t_temp*(dWdy.v.z + dWdz.v.y);

   return (SFS_stress);

}

/**********************************************************************************************
 * LES3DTF_pState::tau_t_z -- Return components of subfilter scale (turbulent)               *
 *                             stress tensor in the z-direction.                              *
 **********************************************************************************************/
Tensor3D LES3DTF_pState::tau_t_z(const LES3DTF_pState &dWdx, 
                                  const LES3DTF_pState &dWdy,
                                  const LES3DTF_pState &dWdz,
                                  const int Flow_Type, 
                                  const double &Volume) {
   Tensor3D SFS_stress;
   double mu_t_temp = mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume);

   SFS_stress.xz = mu_t_temp*(dWdx.v.z + dWdz.v.x);
   SFS_stress.yz = mu_t_temp*(dWdy.v.z + dWdz.v.y);
   SFS_stress.zz = mu_t_temp*(FOUR*dWdz.v.z - TWO*dWdx.v.x - TWO*dWdy.v.y)/THREE - 
                   (TWO/THREE)*rho*k;

   return (SFS_stress);

}

Tensor3D LES3DTF_pState::tau_t_z(const double &mu_t_temp,
                                  const LES3DTF_pState &dWdx, 
                                  const LES3DTF_pState &dWdy,
                                  const LES3DTF_pState &dWdz,
                                  const int Flow_Type, 
                                  const double &Volume) {
   Tensor3D SFS_stress;

   SFS_stress.xz = mu_t_temp*(dWdx.v.z + dWdz.v.x);
   SFS_stress.yz = mu_t_temp*(dWdy.v.z + dWdz.v.y);
   SFS_stress.zz = mu_t_temp*(FOUR*dWdz.v.z - TWO*dWdx.v.x - TWO*dWdy.v.y)/THREE - 
                   (TWO/THREE)*rho*k;

   return (SFS_stress);

}

/************************************************************************
 * LES3DTF_pState::q_t -- Return turbulent heat flux vector.           *
 ************************************************************************/
Vector3D LES3DTF_pState::q_t(const LES3DTF_pState &dWdx, 
                              const LES3DTF_pState &dWdy,
                              const LES3DTF_pState &dWdz,
                              const int Flow_Type, 
                              const double &Volume) {
  
    double Rmix = Rtot();
    double kappa_t_temp = kappa_t(dWdx,dWdy,dWdz,Flow_Type,Volume);
    Vector3D heat_flux;
     
    heat_flux.x = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdx.p -(p/rho)*dWdx.rho);
    heat_flux.y = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdy.p -(p/rho)*dWdy.rho);
    heat_flux.z = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdz.p -(p/rho)*dWdz.rho);
        
    return (heat_flux);

}

Vector3D LES3DTF_pState::q_t(const double &kappa_t_temp,
                              const LES3DTF_pState &dWdx, 
                              const LES3DTF_pState &dWdy,
                              const LES3DTF_pState &dWdz,
                              const int Flow_Type, 
                              const double &Volume) {
  
    double Rmix = Rtot();
    Vector3D heat_flux;
     
    heat_flux.x = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdx.p -(p/rho)*dWdx.rho);
    heat_flux.y = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdy.p -(p/rho)*dWdy.rho);
    heat_flux.z = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdz.p -(p/rho)*dWdz.rho);
        
    return (heat_flux);

}

/************************************************************************
 * LES3DTF_pState::q_t_x -- Return component of turbulent heat flux    *
 *                           vector in the x-direction.                 *
 ************************************************************************/
Vector3D LES3DTF_pState::q_t_x(const LES3DTF_pState &dWdx, 
                                const LES3DTF_pState &dWdy,
                                const LES3DTF_pState &dWdz,
                                const int Flow_Type, 
                                const double &Volume) {
  
    double Rmix = Rtot();
    double kappa_t_temp = kappa_t(dWdx,dWdy,dWdz,Flow_Type,Volume);
    Vector3D heat_flux;
     
    heat_flux.x = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdx.p -(p/rho)*dWdx.rho);
    return (heat_flux);

}

Vector3D LES3DTF_pState::q_t_x(const double &kappa_t_temp,
                                const LES3DTF_pState &dWdx, 
                                const LES3DTF_pState &dWdy,
                                const LES3DTF_pState &dWdz,
                                const int Flow_Type, 
                                const double &Volume) {
  
    double Rmix = Rtot();
    Vector3D heat_flux;
     
    heat_flux.x = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdx.p -(p/rho)*dWdx.rho);
    return (heat_flux);

}

/************************************************************************
 * LES3DTF_pState::q_t_y -- Return component of turbulent heat flux    *
 *                           vector in the y-direction.                 *
 ************************************************************************/
Vector3D LES3DTF_pState::q_t_y(const LES3DTF_pState &dWdx, 
                                const LES3DTF_pState &dWdy,
                                const LES3DTF_pState &dWdz,
                                const int Flow_Type, 
                                const double &Volume) {
  
    double Rmix = Rtot();
    double kappa_t_temp = kappa_t(dWdx,dWdy,dWdz,Flow_Type,Volume);
    Vector3D heat_flux;
     
    heat_flux.y = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdy.p -(p/rho)*dWdy.rho);
    return (heat_flux);

}

Vector3D LES3DTF_pState::q_t_y(const double &kappa_t_temp,
                                const LES3DTF_pState &dWdx, 
                                const LES3DTF_pState &dWdy,
                                const LES3DTF_pState &dWdz,
                                const int Flow_Type, 
                                const double &Volume) {
  
    double Rmix = Rtot();
    Vector3D heat_flux;
     
    heat_flux.y = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdy.p -(p/rho)*dWdy.rho);
    return (heat_flux);

}

/************************************************************************
 * LES3DTF_pState::q_t_z -- Return component of turbulent heat flux    *
 *                           vector in the z-direction.                 *
 ************************************************************************/
Vector3D LES3DTF_pState::q_t_z(const LES3DTF_pState &dWdx, 
                                const LES3DTF_pState &dWdy,
                                const LES3DTF_pState &dWdz,
                                const int Flow_Type, 
                                const double &Volume) {
  
    double Rmix = Rtot();
    double kappa_t_temp = kappa_t(dWdx,dWdy,dWdz,Flow_Type,Volume);
    Vector3D heat_flux;
     
    heat_flux.z = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdz.p -(p/rho)*dWdz.rho);
    return (heat_flux);

}

Vector3D LES3DTF_pState::q_t_z(const double &kappa_t_temp,
                                const LES3DTF_pState &dWdx, 
                                const LES3DTF_pState &dWdy,
                                const LES3DTF_pState &dWdz,
                                const int Flow_Type, 
                                const double &Volume) {
  
    double Rmix = Rtot();
    Vector3D heat_flux;
     
    heat_flux.z = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdz.p -(p/rho)*dWdz.rho);
    return (heat_flux);

}

/*******************************************************************
 * LES3DTF_pState::U -- Return conserved solution state vector.   *
 *******************************************************************/
LES3DTF_cState LES3DTF_pState::U(void) {
   LES3DTF_cState Temp;
   Temp.rho = rho;
   Temp.rhov = rhov();
   Temp.E = E();
   Temp.rhok = rho*k;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i] = rho*spec[i];
   } /* endfor */
   Temp.flame = flame;
   return Temp;
}

LES3DTF_cState LES3DTF_pState::U(void) const {
   LES3DTF_cState Temp;
   Temp.rho = rho;
   Temp.rhov = rhov();
   Temp.E = E();
   Temp.rhok = rho*k;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i] = rho*spec[i];
   } /* endfor */
   Temp.flame = flame;
   return Temp;
}

LES3DTF_cState LES3DTF_pState::U(const LES3DTF_pState &W) {
    LES3DTF_cState Temp;
    Temp.rho = W.rho;
    Temp.rhov = W.rhov();
    Temp.E = W.E();
    Temp.rhok = W.rho*W.k;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i] = rho*spec[i];
   } /* endfor */
    Temp.flame = W.flame;
    return Temp;
}

/***********************************************************************
 ***************** INVISCID FLUX VECTORS *******************************
 ***********************************************************************/

/***********************************************************
 *    LES3DTF_pState::F -- Inviscid flux (x-direction).   *
 ***********************************************************/
LES3DTF_cState LES3DTF_pState::F(void) {
   LES3DTF_cState Temp;
   Temp.rho = rho*v.x;
   Temp.rhov.x = rho*sqr(v.x) + p + (TWO/THREE)*rho*k;
   Temp.rhov.y = rho*v.x*v.y;
   Temp.rhov.z = rho*v.x*v.z;
   Temp.E = v.x*H();
   Temp.rhok = rho*v.x*k;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i].c = rho*v.x*spec[i].c;
   } /* endfor */
   return (Temp);
}

LES3DTF_cState LES3DTF_pState::F(void) const {
   LES3DTF_cState Temp;
   Temp.rho = rho*v.x;
   Temp.rhov.x = rho*sqr(v.x) + p + (TWO/THREE)*rho*k;
   Temp.rhov.y = rho*v.x*v.y;
   Temp.rhov.z = rho*v.x*v.z;
   Temp.E = v.x*H();
   Temp.rhok = rho*v.x*k;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i].c = rho*v.x*spec[i].c;
   } /* endfor */
   return (Temp);
}

/************************************************************
 *    LES3DTF_pState::Fx -- Inviscid flux (x-direction).   *
 ************************************************************/
LES3DTF_cState LES3DTF_pState::Fx(void) {
   LES3DTF_cState Temp;
   Temp.rho = rho*v.x;
   Temp.rhov.x = rho*sqr(v.x) + p + (TWO/THREE)*rho*k;
   Temp.rhov.y = rho*v.x*v.y;
   Temp.rhov.z = rho*v.x*v.z;
   Temp.E = v.x*H();
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i].c = rho*v.x*spec[i].c;
   } /* endfor */
   Temp.rhok = rho*v.x*k;
   return (Temp);
}

LES3DTF_cState LES3DTF_pState::Fx(void) const {
   LES3DTF_cState Temp;
   Temp.rho = rho*v.x;
   Temp.rhov.x = rho*sqr(v.x) + p + (TWO/THREE)*rho*k;
   Temp.rhov.y = rho*v.x*v.y;
   Temp.rhov.z = rho*v.x*v.z;
   Temp.E = v.x*H();
   Temp.rhok = rho*v.x*k;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i].c = rho*v.x*spec[i].c;
   } /* endfor */
   return (Temp);
}

/************************************************************
 *    LES3DTF_pState::Fy -- Inviscid flux (y-direction).   *
 ************************************************************/
LES3DTF_cState LES3DTF_pState::Fy(void) {
   LES3DTF_cState Temp;
   Temp.rho = rho*v.y;
   Temp.rhov.x = rho*v.x*v.y;
   Temp.rhov.y = rho*sqr(v.y) + p + (TWO/THREE)*rho*k;
   Temp.rhov.z = rho*v.y*v.z;
   Temp.E = v.y*H();
   Temp.rhok = rho*v.y*k;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i].c = rho*v.y*spec[i].c;
   } /* endfor */
   return (Temp);
}

LES3DTF_cState LES3DTF_pState::Fy(void) const {
   LES3DTF_cState Temp;
   Temp.rho = rho*v.y;
   Temp.rhov.x = rho*v.x*v.y;
   Temp.rhov.y = rho*sqr(v.y) + p + (TWO/THREE)*rho*k;
   Temp.rhov.z = rho*v.y*v.z;
   Temp.E = v.y*H();
   Temp.rhok = rho*v.y*k;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i].c = rho*v.y*spec[i].c;
   } /* endfor */
   return (Temp);
}

/************************************************************
 *    LES3DTF_pState::Fz -- Inviscid flux (z-direction).   *
 ************************************************************/
LES3DTF_cState LES3DTF_pState::Fz(void) {
   LES3DTF_cState Temp;
   Temp.rho = rho*v.z;
   Temp.rhov.x = rho*v.x*v.z;
   Temp.rhov.y = rho*v.y*v.z;
   Temp.rhov.z = rho*sqr(v.z) + p + (TWO/THREE)*rho*k;
   Temp.E = v.z*H();
   Temp.rhok = rho*v.z*k;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i].c = rho*v.z*spec[i].c;
   } /* endfor */
   return (Temp);
}

LES3DTF_cState LES3DTF_pState::Fz(void) const {
   LES3DTF_cState Temp;
   Temp.rho = rho*v.z;
   Temp.rhov.x = rho*v.x*v.z;
   Temp.rhov.y = rho*v.y*v.z;
   Temp.rhov.z = rho*sqr(v.z) + p + (TWO/THREE)*rho*k;
   Temp.E = v.z*H();
   Temp.rhok = rho*v.z*k;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i].c = rho*v.z*spec[i].c;
   } /* endfor */
   return (Temp);
}

/**********************************************************************
 ***************** VISCOUS FLUX VECTORS *******************************
 **********************************************************************/

/********************************************************
 *  LES3DTF_pState::Fv -- Viscous flux (x-direction).  * 
 ********************************************************/
LES3DTF_cState LES3DTF_pState::Fv(const LES3DTF_pState &dWdx,
                                    const LES3DTF_pState &dWdy,
                                    const LES3DTF_pState &dWdz,
                                    const int Flow_Type,
                                    const double &Volume) {

   LES3DTF_cState Temp;

   double mu_temp = flame.WF*flame.TF*mu();
   double mu_t_temp = flame.WF*flame.TF*mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume);
   double kappa_temp = flame.WF*flame.TF*kappa();
   double kappa_t_temp = flame.WF*flame.TF*kappa_t(dWdx,dWdy,dWdz,Flow_Type,Volume);

   Tensor3D fluid_stress;
   Vector3D heat_flux;
   
   fluid_stress = tau_t_x(mu_temp + mu_t_temp, dWdx, dWdy, dWdz, Flow_Type, Volume);
   heat_flux = q_t_x(kappa_temp + kappa_t_temp, dWdx, dWdy, dWdz, Flow_Type, Volume);
   
   for (int index = 0; index < ns; ++index) {
     _diff_coeff[index] = flame.WF*flame.TF*(Ds(index, mu_temp)+Ds_t(index, mu_t_temp));
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
            heat_flux.x;
   Temp.rhok = (mu_temp+mu_t_temp/0.25)*dWdx.k;

   for (int index = 0; index<ns; ++index) {
   Temp.rhospec[index] = rho*_diff_coeff[index]*dWdx.spec[index].c;
   } /* endfor */

   return (Temp);

}

/*********************************************************
 *  LES3DTF_pState::Fvx -- Viscous flux (x-direction).  * 
 *********************************************************/
LES3DTF_cState LES3DTF_pState::Fvx(const LES3DTF_pState &dWdx,
                                     const LES3DTF_pState &dWdy,
                                     const LES3DTF_pState &dWdz,
                                     const int Flow_Type,
                                     const double &Volume) {

   LES3DTF_cState Temp;

   double mu_temp = flame.WF*flame.TF*mu();
   double mu_t_temp = flame.WF*flame.TF*mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume);
   double kappa_temp = flame.WF*flame.TF*kappa();
   double kappa_t_temp = flame.WF*flame.TF*kappa_t(dWdx,dWdy,dWdz,Flow_Type,Volume);

   Tensor3D fluid_stress;
   Vector3D heat_flux;
   
   fluid_stress = tau_t_x(mu_temp + mu_t_temp, dWdx, dWdy, dWdz, Flow_Type, Volume);
   heat_flux = q_t_x(kappa_temp + kappa_t_temp, dWdx, dWdy, dWdz, Flow_Type, Volume);
   
   for (int index = 0; index < ns; ++index) {
     _diff_coeff[index] = flame.WF*flame.TF*(Ds(index, mu_temp)+Ds_t(index, mu_t_temp));
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
            heat_flux.x;
   Temp.rhok = (mu_temp+mu_t_temp/0.25)*dWdx.k;
   for (int index = 0; index<ns; ++index) {
   Temp.rhospec[index] = rho*_diff_coeff[index]*dWdx.spec[index].c;
   } /* endfor */

   return (Temp);

}

/*********************************************************
 *  LES3DTF_pState::Fvy -- Viscous flux (y-direction).  * 
 *********************************************************/
LES3DTF_cState LES3DTF_pState::Fvy(const LES3DTF_pState &dWdx,
                                     const LES3DTF_pState &dWdy,
                                     const LES3DTF_pState &dWdz,
                                     const int Flow_Type,
                                     const double &Volume) {

   LES3DTF_cState Temp;

   double mu_temp = flame.WF*flame.TF*mu();
   double mu_t_temp = flame.WF*flame.TF*mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume);
   double kappa_temp = flame.WF*flame.TF*kappa();
   double kappa_t_temp = flame.WF*flame.TF*kappa_t(dWdx,dWdy,dWdz,Flow_Type,Volume);

   Tensor3D fluid_stress;
   Vector3D heat_flux;
   
   fluid_stress = tau_t_y(mu_temp + mu_t_temp, dWdx, dWdy, dWdz, Flow_Type, Volume);
   heat_flux = q_t_y(kappa_temp + kappa_t_temp, dWdx, dWdy, dWdz, Flow_Type, Volume);
   
   for (int index = 0; index < ns; ++index) {
     _diff_coeff[index] = flame.WF*flame.TF*(Ds(index, mu_temp)+Ds_t(index, mu_t_temp));
   } /* endfor */

   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux -= rho*thermal_diffusion_x(dWdx, dWdy, dWdz);
   
   Temp.rho= ZERO;
   Temp.rhov.x = fluid_stress.xy;
   Temp.rhov.y = fluid_stress.yy + (TWO/THREE)*rho*k ;
   Temp.rhov.z = fluid_stress.yz;
   Temp.E = v.x*fluid_stress.xy +
            v.y*(fluid_stress.yy + (TWO/THREE)*rho*k) +
            v.z*fluid_stress.yz - 
            heat_flux.y;
   Temp.rhok = (mu_temp+mu_t_temp/0.25)*dWdy.k;
   for (int index = 0; index<ns; ++index) {
   Temp.rhospec[index] = rho*_diff_coeff[index]*dWdy.spec[index].c;
   } /* endfor */

   return (Temp);

}

/*********************************************************
 *  LES3DTF_pState::Fvz -- Viscous flux (z-direction).  * 
 *********************************************************/
LES3DTF_cState LES3DTF_pState::Fvz(const LES3DTF_pState &dWdx,
                                     const LES3DTF_pState &dWdy,
                                     const LES3DTF_pState &dWdz,
                                     const int Flow_Type,
                                     const double &Volume) {

   LES3DTF_cState Temp;

   double mu_temp = flame.WF*flame.TF*mu();
   double mu_t_temp = flame.WF*flame.TF*mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume);
   double kappa_temp = flame.WF*flame.TF*kappa();
   double kappa_t_temp = flame.WF*flame.TF*kappa_t(dWdx,dWdy,dWdz,Flow_Type,Volume);

   Tensor3D fluid_stress;
   Vector3D heat_flux;
   
   fluid_stress = tau_t_z(mu_temp + mu_t_temp, dWdx, dWdy, dWdz, Flow_Type, Volume);
   heat_flux = q_t_z(kappa_temp + kappa_t_temp, dWdx, dWdy, dWdz, Flow_Type, Volume);
   
   for (int index = 0; index < ns; ++index) {
     _diff_coeff[index] = flame.WF*flame.TF*(Ds(index, mu_temp)+Ds_t(index, mu_t_temp));
   } /* endfor */

   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux -= rho*thermal_diffusion_x(dWdx, dWdy, dWdz);
   
   Temp.rho= ZERO;
   Temp.rhov.x = fluid_stress.xz;
   Temp.rhov.y = fluid_stress.yz;
   Temp.rhov.z = fluid_stress.zz + (TWO/THREE)*rho*k;
   Temp.E = v.x*fluid_stress.xz +
            v.y*fluid_stress.yz +
            v.z*(fluid_stress.zz + (TWO/THREE)*rho*k) - 
            heat_flux.z;
   Temp.rhok = (mu_temp+mu_t_temp/0.25)*dWdz.k;
   for (int index = 0; index<ns; ++index) {
   Temp.rhospec[index] = rho*_diff_coeff[index]*dWdz.spec[index].c;
   } /* endfor */

   return (Temp);

}

/************************************************************************
 ********************** EIGENVALUES *************************************
 ************************************************************************/

/**************************************************************
 * LES3DTF_pState::lambda -- Eigenvalue(s) (x-direction).    *
 **************************************************************/
LES3DTF_pState LES3DTF_pState::lambda(void) {
  double cc = a_t();
  LES3DTF_pState Temp;
  Temp.rho = v.x - cc;
  Temp.v.x = v.x;
  Temp.v.y = v.x;
  Temp.v.z = v.x;
  Temp.p = v.x + cc;
  Temp.k = v.x;
  for (int i = 0; i < ns; i++) {
     Temp.spec[i].c = v.x;
  } /* endfor */
  return (Temp);
}

LES3DTF_pState LES3DTF_pState::lambda(void) const {
  double cc = a_t();
  LES3DTF_pState Temp;
  Temp.rho = v.x - cc;
  Temp.v.x = v.x;
  Temp.v.y = v.x;
  Temp.v.z = v.x;
  Temp.p = v.x + cc;
  Temp.k = v.x;
  for (int i = 0; i < ns; i++) {
     Temp.spec[i].c = v.x;
  } /* endfor */
  return (Temp);
}

/**************************************************************
 * LES3DTF_pState::lambda_x -- Eigenvalue(s) (x-direction).  *
 **************************************************************/
LES3DTF_pState LES3DTF_pState::lambda_x(void) {
  double cc = a_t();
  LES3DTF_pState Temp;
  Temp.rho = v.x - cc;
  Temp.v.x = v.x;
  Temp.v.y = v.x;
  Temp.v.z = v.x;
  Temp.p = v.x + cc;
  Temp.k = v.x;
  for (int i = 0; i < ns; i++) {
     Temp.spec[i].c = v.x;
  } /* endfor */
  return (Temp);
}

LES3DTF_pState LES3DTF_pState::lambda_x(void) const {
  double cc = a_t();
  LES3DTF_pState Temp;
  Temp.rho = v.x - cc;
  Temp.v.x = v.x;
  Temp.v.y = v.x;
  Temp.v.z = v.x;
  Temp.p = v.x + cc;
  Temp.k = v.x;
  for (int i = 0; i < ns; i++) {
     Temp.spec[i].c = v.x;
  } /* endfor */
  return (Temp);
}

/************************************************************************
 ********************* EIGENVECTORS *************************************
 ************************************************************************/

/***********************************************************************
 * LES3DTF_pState::rc -- Conserved right eigenvector (x-direction).   *
 ***********************************************************************/
LES3DTF_cState LES3DTF_pState::rc(const int &index) {
   double cc;
   switch(index){  
   case 1:
     cc = a_t();
     return (LES3DTF_cState(ONE, v.x-cc, v.y, v.z, H()/rho+2.0/3.0*k-v.x*cc, k, spec));
   case 2:
     cc = a_t();
     return (LES3DTF_cState(ONE, v.x, v.y, v.z, H()/rho+2.0/3.0*k-cc*cc/(g()-1.0), k, spec)); 
   case 3:
     return (LES3DTF_cState(ZERO, ZERO, rho, ZERO, rho*v.y, ZERO, ZERO));
   case 4:
     return (LES3DTF_cState(ZERO, ZERO, ZERO, rho, rho*v.z, ZERO, ZERO));
   case 5:
     cc = a_t();
     return (LES3DTF_cState(ONE, v.x+cc, v.y, v.z, H()/rho+2.0/3.0*k+v.x*cc, k, spec));
   case 6:
     return (LES3DTF_cState(ZERO, ZERO, ZERO, ZERO, 5.0*rho/3.0, rho, ZERO));
   default :
      LES3DTF_cState rr(ZERO);
      double RTOT = Rtot();
      double TEMP = p/(rho*RTOT);      
      rr.E = rho*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].Enthalpy(TEMP) + 
                   specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].Heatofform() 
                   - Cp(TEMP)*TEMP*specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_LES3DTF_VAR_EXTRA+1)].Rs()/RTOT);
      double PHI = specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].Enthalpy(TEMP)-
                   specdata[num_vars-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].Enthalpy(TEMP)-
                   (Cp(TEMP) -RTOT)*TEMP*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_LES3DTF_VAR_EXTRA+1)].Rs() - specdata[num_vars - (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_LES3DTF_VAR_EXTRA+1)].Rs())/RTOT;
      rr.rhospec[index - (NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].c = rho*PHI; 
      return rr;
    };
}

LES3DTF_cState LES3DTF_pState::rc(const int &index) const{
   double cc;
   switch(index){  
   case 1:
     cc = a_t();
     return (LES3DTF_cState(ONE, v.x-cc, v.y, v.z, H()/rho+2.0/3.0*k-v.x*cc, k, spec));
   case 2:
     cc = a_t();
     return (LES3DTF_cState(ONE, v.x, v.y, v.z, H()/rho+2.0/3.0*k-cc*cc/(g()-1.0), k, spec)); 
   case 3:
     return (LES3DTF_cState(ZERO, ZERO, rho, ZERO, rho*v.y, ZERO, ZERO));
   case 4:
     return (LES3DTF_cState(ZERO, ZERO, ZERO, rho, rho*v.z, ZERO, ZERO));
   case 5:
     cc = a_t();
     return (LES3DTF_cState(ONE, v.x+cc, v.y, v.z, H()/rho+2.0/3.0*k+v.x*cc, k, spec));
   case 6:
     return (LES3DTF_cState(ZERO, ZERO, ZERO, ZERO, 5.0*rho/3.0, rho, ZERO));
   default :
      LES3DTF_cState rr(ZERO);
      double RTOT = Rtot();
      double TEMP = p/(rho*RTOT);      
      rr.E = rho*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].Enthalpy(TEMP) + 
                   specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].Heatofform() 
                   - Cp(TEMP)*TEMP*specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_LES3DTF_VAR_EXTRA+1)].Rs()/RTOT);
      double PHI = specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].Enthalpy(TEMP)-
                   specdata[num_vars-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].Enthalpy(TEMP)-
                   (Cp(TEMP) -RTOT)*TEMP*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_LES3DTF_VAR_EXTRA+1)].Rs() - specdata[num_vars - (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_LES3DTF_VAR_EXTRA+1)].Rs())/RTOT;
      rr.rhospec[index - (NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].c = rho*PHI; 
      return rr;
    };
}

/*************************************************************************
 * LES3DTF_pState::rc_x -- Conserved right eigenvector (x-direction).   *
 *************************************************************************/
LES3DTF_cState LES3DTF_pState::rc_x(const int &index) const{
   double cc;
   switch(index){  
   case 1:
     cc = a_t();
     return (LES3DTF_cState(ONE, v.x-cc, v.y, v.z, H()/rho+2.0/3.0*k-v.x*cc, k, spec));
   case 2:
     cc = a_t();
     return (LES3DTF_cState(ONE, v.x, v.y, v.z, H()/rho+2.0/3.0*k-cc*cc/(g()-1.0), k, spec)); 
   case 3:
     return (LES3DTF_cState(ZERO, ZERO, rho, ZERO, rho*v.y, ZERO, ZERO));
   case 4:
     return (LES3DTF_cState(ZERO, ZERO, ZERO, rho, rho*v.z, ZERO, ZERO));
   case 5:
     cc = a_t();
     return (LES3DTF_cState(ONE, v.x+cc, v.y, v.z, H()/rho+2.0/3.0*k+v.x*cc, k, spec));
   case 6:
     return (LES3DTF_cState(ZERO, ZERO, ZERO, ZERO, 5.0*rho/3.0, rho, ZERO));
   default :
      LES3DTF_cState rr(ZERO);
      double RTOT = Rtot();
      double TEMP = p/(rho*RTOT);      
      rr.E = rho*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].Enthalpy(TEMP) + 
                   specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].Heatofform() 
                   - Cp(TEMP)*TEMP*specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_LES3DTF_VAR_EXTRA+1)].Rs()/RTOT);
      double PHI = specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].Enthalpy(TEMP)-
                   specdata[num_vars-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].Enthalpy(TEMP)-
                   (Cp(TEMP) -RTOT)*TEMP*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_LES3DTF_VAR_EXTRA+1)].Rs() - specdata[num_vars - (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_LES3DTF_VAR_EXTRA+1)].Rs())/RTOT;
      rr.rhospec[index - (NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].c = rho*PHI; 
      return rr;
    };
}

LES3DTF_cState LES3DTF_pState::rc_x(const int &index) {
   double cc;
   switch(index){  
   case 1:
     cc = a_t();
     return (LES3DTF_cState(ONE, v.x-cc, v.y, v.z, H()/rho+2.0/3.0*k-v.x*cc, k, spec));
   case 2:
     cc = a_t();
     return (LES3DTF_cState(ONE, v.x, v.y, v.z, H()/rho+2.0/3.0*k-cc*cc/(g()-1.0), k, spec)); 
   case 3:
     return (LES3DTF_cState(ZERO, ZERO, rho, ZERO, rho*v.y, ZERO, ZERO));
   case 4:
     return (LES3DTF_cState(ZERO, ZERO, ZERO, rho, rho*v.z, ZERO, ZERO));
   case 5:
     cc = a_t();
     return (LES3DTF_cState(ONE, v.x+cc, v.y, v.z, H()/rho+2.0/3.0*k+v.x*cc, k, spec));
   case 6:
     return (LES3DTF_cState(ZERO, ZERO, ZERO, ZERO, 5.0*rho/3.0, rho, ZERO));
   default :
      LES3DTF_cState rr(ZERO);
      double RTOT = Rtot();
      double TEMP = p/(rho*RTOT);      
      rr.E = rho*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].Enthalpy(TEMP) + 
                   specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].Heatofform() 
                   - Cp(TEMP)*TEMP*specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_LES3DTF_VAR_EXTRA+1)].Rs()/RTOT);
      double PHI = specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].Enthalpy(TEMP)-
                   specdata[num_vars-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].Enthalpy(TEMP)-
                   (Cp(TEMP) -RTOT)*TEMP*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_LES3DTF_VAR_EXTRA+1)].Rs() - specdata[num_vars - (NUM_EULER3D_VAR_SANS_SPECIES + 
                   NUM_LES3DTF_VAR_EXTRA+1)].Rs())/RTOT;
      rr.rhospec[index - (NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].c = rho*PHI; 
      return rr;
    };
}

/***********************************************************************
 * LES3DTF_pState::lp -- Primitive left eigenvector (x-direction).    *
 ***********************************************************************/
LES3DTF_pState LES3DTF_pState::lp(const int &index) {
   double cc;  
   switch(index){  
   case 1:
      cc = a_t();
      return (LES3DTF_pState(ZERO, -HALF*rho/cc, ZERO, ZERO, HALF/(cc*cc), ZERO, ZERO));
   case 2:
      cc = a_t();      
      return (LES3DTF_pState(ONE, ZERO, ZERO, ZERO, -ONE/(cc*cc), ZERO, ZERO));
   case 3 :
      return (LES3DTF_pState(ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO));
   case 4 :
      return (LES3DTF_pState(ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO));
   case 5 :
      cc = a_t();
      return (LES3DTF_pState(ZERO, HALF*rho/cc, ZERO, ZERO, HALF/(cc*cc), ZERO, ZERO));
   case 6 :
     return (LES3DTF_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO));
   default :
      LES3DTF_pState NEW(ZERO);
      NEW.spec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].c = ONE;
      return NEW;
   };
}

LES3DTF_pState LES3DTF_pState::lp(const int &index) const {
   double cc;  
   switch(index){  
   case 1:
      cc = a_t();
      return (LES3DTF_pState(ZERO, -HALF*rho/cc, ZERO, ZERO, HALF/(cc*cc), ZERO, ZERO));
   case 2:
      cc = a_t();      
      return (LES3DTF_pState(ONE, ZERO, ZERO, ZERO, -ONE/(cc*cc), ZERO, ZERO));
   case 3 :
      return (LES3DTF_pState(ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO));
   case 4 :
      return (LES3DTF_pState(ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO));
   case 5 :
      cc = a_t();
      return (LES3DTF_pState(ZERO, HALF*rho/cc, ZERO, ZERO, HALF/(cc*cc), ZERO, ZERO));
   case 6 :
     return (LES3DTF_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO));
   default :
      LES3DTF_pState NEW(ZERO);
      NEW.spec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].c = ONE;
      return NEW;
   };
}

/*************************************************************************
 * LES3DTF_pState::lp_x -- Primitive left eigenvector (x-direction).    *
 *************************************************************************/
LES3DTF_pState LES3DTF_pState::lp_x(const int &index) {
   double cc;  
   switch(index){  
   case 1:
      cc = a_t();
      return (LES3DTF_pState(ZERO, -HALF*rho/cc, ZERO, ZERO, HALF/(cc*cc), ZERO, ZERO));
   case 2:
      cc = a_t();      
      return (LES3DTF_pState(ONE, ZERO, ZERO, ZERO, -ONE/(cc*cc), ZERO, ZERO));
   case 3 :
      return (LES3DTF_pState(ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO));
   case 4 :
      return (LES3DTF_pState(ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO));
   case 5 :
      cc = a_t();
      return (LES3DTF_pState(ZERO, HALF*rho/cc, ZERO, ZERO, HALF/(cc*cc), ZERO, ZERO));
   case 6 :
     return (LES3DTF_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO));
   default :
      LES3DTF_pState NEW(ZERO);
      NEW.spec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].c = ONE;
      return NEW;
   };
}

LES3DTF_pState LES3DTF_pState::lp_x(const int &index) const {
   double cc;  
   switch(index){  
   case 1:
      cc = a_t();
      return (LES3DTF_pState(ZERO, -HALF*rho/cc, ZERO, ZERO, HALF/(cc*cc), ZERO, ZERO));
   case 2:
      cc = a_t();      
      return (LES3DTF_pState(ONE, ZERO, ZERO, ZERO, -ONE/(cc*cc), ZERO, ZERO));
   case 3 :
      return (LES3DTF_pState(ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO));
   case 4 :
      return (LES3DTF_pState(ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO));
   case 5 :
      cc = a_t();
      return (LES3DTF_pState(ZERO, HALF*rho/cc, ZERO, ZERO, HALF/(cc*cc), ZERO, ZERO));
   case 6 :
     return (LES3DTF_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO));
   default :
      LES3DTF_pState NEW(ZERO);
      NEW.spec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].c = ONE;
      return NEW;
   };
}

/***********************************************************
 * LES3DTF_pState::Mr2 -- Square of Mach number           *
 ***********************************************************/
double LES3DTF_pState::Mr2(const double &deltax, 
                            const double &lengthx, 
                            const double &dTime) {

    double cc = a_t();
    double MR2 = min(max((v.sqr()/(cc*cc)),Mref*Mref),ONE);
    MR2 = pow(max(sqrt(MR2*cc*cc), mu()/(rho*deltax)),2.0)/(cc*cc);
    double MR_uns = (lengthx/(PI*dTime))/cc;  // ZERO;
// double MR_vis;
//     if (flow_type_flag == FLOWTYPE_INVISCID) {
//       MR_vis = ZERO;
//     } else {
//       MR_vis = (mu()/(rho*deltax))/c;
//     }  
   //  MR2 = max(MR_inv*MR_inv, MR_vis*MR_vis);
    MR2 = max(MR2, MR_uns*MR_uns);
    MR2 = min(max(MR2, Mref*Mref), ONE);

    return (MR2);

}

/***********************************************************************************
 * LES3DTF_pState::u_plus_aprecon -- Preconditioned velocity                      *
 ***********************************************************************************/
double LES3DTF_pState::u_plus_aprecon(const double &u, 
                                       const double &deltax, 
                                       const double &lengthx, 
                                       const double &dTime) {

  double Temp = T();
  double cc = a_t();
  double UR2 = Mr2(deltax, lengthx, dTime)*cc*cc;
  double alpha = HALF*( ONE - (ONE/(Rtot()*Temp) - ONE/(Cp()*Temp))*UR2);
  // uprime + cprime

  return ( u*(ONE - alpha) + sqrt(alpha*alpha*u*u + UR2) );

}

/***********************************************************************************
 * LES3DTF_pState::u_a_precon -- Preconditioned velocity and soundspeed           *
 ***********************************************************************************/
void LES3DTF_pState::u_a_precon(const double &UR2, 
                                 double &uprimed, 
                                 double &cprimed) const {

  double alpha = HALF*(ONE - ONE/(a_t()*a_t())*UR2);
  uprimed = v.x*(ONE - alpha);
  cprimed = sqrt(alpha*alpha*v.x*v.x + UR2); 

}

/******************************************************************************
 * LES3DTF_pState::rc_x_precon -- Conserved right eigenvector (x-direction)  * 
 *                                 for low Mach number preconditioner.        *
 ******************************************************************************/
LES3DTF_cState LES3DTF_pState::rc_x_precon(const int &index, const double &MR2) const {

   double cc, uprimed, cprimed;

   switch(index){  
   case 1:
     cc = a_t();
     u_a_precon(MR2*cc*cc,uprimed,cprimed);
     return (LES3DTF_cState(ONE, 
                           (uprimed-cprimed)/MR2, v.y, v.z, 
                            h()+(v.sqr()/MR2)/TWO+(FIVE*k)/THREE-(v.x*cprimed)/MR2, 
		            k,
                            spec));
   case 2:
     cc = a_t();
     return (LES3DTF_cState(ONE, v.x, v.y, v.z, (H()/rho-cc*cc/(g()-ONE)), k, spec)); 
   case 3:
     return (LES3DTF_cState(ZERO, ZERO, rho, ZERO, rho*v.y, ZERO, ZERO));
   case 4:
     return (LES3DTF_cState(ZERO, ZERO, ZERO, rho, rho*v.z, ZERO, ZERO));
   case 5:
     cc = a_t();
     u_a_precon(MR2*cc*cc,uprimed,cprimed);
     return (LES3DTF_cState(ONE, 
			   (uprimed+cprimed)/MR2, v.y, v.z, 
                            h()+(v.sqr()/MR2)/TWO+(FIVE*k)/THREE+(v.x*cprimed)/MR2, 
                            k,
                            spec));
   case 6:
     return (LES3DTF_cState(ZERO, ZERO, ZERO, ZERO, FIVE*rho/THREE, rho, ZERO));
   };

}

/******************************************************************************
 * LES3DTF_pState::lp_x_precon -- Primitive left eigenvector (x-direction)   * 
 *                                 for low Mach number preconditioner.        *
 ******************************************************************************/
LES3DTF_pState LES3DTF_pState::lp_x_precon(const int &index,const double &MR2) const {

   double cc, uprimed,cprimed;

   switch(index){  
   case 1:
      cc = a_t();
      u_a_precon(MR2*cc*cc,uprimed,cprimed);
      return (LES3DTF_pState(ZERO, -HALF*rho*MR2/cprimed, ZERO, ZERO, 
                             (-uprimed+cprimed+v.x)/(TWO*cprimed*cc*cc), 
                              ZERO, ZERO));
   case 2:
      cc = a_t();
      return (LES3DTF_pState(ONE, ZERO, ZERO, ZERO, -ONE/(cc*cc), ZERO, ZERO));
   case 3 :
      return (LES3DTF_pState(ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO));
   case 4 :
      return (LES3DTF_pState(ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO));
   case 5 :
      cc = a_t();
      u_a_precon(MR2*cc*cc,uprimed,cprimed);
      return (LES3DTF_pState(ZERO, HALF*rho*MR2/cprimed, ZERO, ZERO, (
                              uprimed+cprimed-v.x)/(TWO*cprimed*cc*cc), 
                              ZERO, ZERO));
   case 6 :
     return (LES3DTF_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO));
   };

}

/******************************************************************************
 * LES3DTF_pState::Low_Mach_Number_Preconditioner -- Preconditoner matrix    * 
 ******************************************************************************/
void LES3DTF_pState::Low_Mach_Number_Preconditioner(DenseMatrix &P, 
                                                     const double &deltax, 
                                                     const double &lengthx,
                                                     const double &dTime) {  
  double Temp = T();
  double pt = p_t();
  double Rmix = Rtot();
  double enthalpy = h();
  double CP = Cp();
  double c = a();
  double theta = ONE/(Mr2(deltax, lengthx, dTime)*c*c) 
                 + (g()-ONE)/(c*c);  
 
 
  double phi = ZERO;
  for(int j=0; j<ns-1; ++j){   
    phi += spec[j].c*(specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() 
		      - CP*Temp*specdata[j].Rs()/Rmix);   
  }


  double alpha = theta*pt/rho;
  double alpham1 = alpha - ONE;
  double omega = (Rmix - CP)*pt/(rho*Rmix);
  double beta = enthalpy - CP*pt/(rho*Rmix) - phi + 5.0/3.0*k;
  double V = HALF*v.sqr();


  P.zero();


  P(0,0) = (alpha*(beta-V)+V+pt/rho-enthalpy+phi-5.0/3.0*k)/omega;    
  P(0,1) = v.x*alpham1/omega;                               
  P(0,2) = v.y*alpham1/omega;                               
  P(0,3) = v.z*alpham1/omega;
  P(0,4) = -alpham1/omega;

  P(1,0) = v.x*(beta-V)*alpham1/omega;                      
  P(1,1) = v.x*v.x*alpham1/omega + 1.0;                     
  P(1,2) = v.x*v.y*alpham1/omega;                           
  P(1,3) = v.x*v.z*alpham1/omega;
  P(1,4) = -v.x*alpham1/omega;

  P(2,0) = v.y*(beta-V)*alpham1/omega;                      
  P(2,1) = v.x*v.y*alpham1/omega;                           
  P(2,2) = v.y*v.y*alpham1/omega + 1.0;                     
  P(2,3) = v.z*v.y*alpham1/omega;
  P(2,4) = -v.y*alpham1/omega;

  P(3,0) = v.z*(beta-V)*alpham1/omega;  
  P(3,1) = v.x*v.z*alpham1/omega;       
  P(3,2) = v.y*v.z*alpham1/omega;       
  P(3,3) = v.z*v.z*alpham1/omega + 1.0; 
  P(3,4) = -v.z*alpham1/omega;

  P(4,0) = (enthalpy+V+5.0/3.0*k)*(beta-V)*alpham1/omega;  
  P(4,1) = v.x*(enthalpy+V+5.0/3.0*k)*alpham1/omega;       
  P(4,2) = v.y*(enthalpy+V+5.0/3.0*k)*alpham1/omega;       
  P(4,3) = v.z*(enthalpy+V+5.0/3.0*k)*alpham1/omega; 
  P(4,4) = -(alpha*(enthalpy+V+5.0/3.0*k)-V-5.0/3.0*k-pt/rho-beta-phi)/omega;
  

  int NUM_VAR = NUM_EULER3D_VAR_SANS_SPECIES;
  

  // Multispecies
  for (int j=0; j<ns-1; ++j) {       

    double enth_j = specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/Rmix;
    P(0,j+NUM_VAR) = enth_j*alpham1/omega;       
    P(1,j+NUM_VAR) = v.x*enth_j*alpham1/omega;   
    P(2,j+NUM_VAR) = v.y*enth_j*alpham1/omega;   
    P(3,j+NUM_VAR) = v.z*enth_j*alpham1/omega;  
    P(4,j+NUM_VAR) = enth_j*(enthalpy+V+5.0/3.0*k)*alpham1/omega;

    for (int i=0; i<ns-1; ++i) { 
      if (i==j) { 
	P(i+NUM_VAR,0) = (spec[i].c)*(beta-V)*alpham1/omega;   
	P(i+NUM_VAR,1) = (spec[i].c)*v.x*alpham1/omega;        
	P(i+NUM_VAR,2) = (spec[i].c)*v.y*alpham1/omega;        
	P(i+NUM_VAR,3) = (spec[i].c)*v.z*alpham1/omega;
	P(i+NUM_VAR,4) = -(spec[i].c)*alpham1/omega;
	//diagonal
	P(i+NUM_VAR,j+NUM_VAR) = (spec[i].c)*enth_j*alpham1/omega + 1.0;    
      } else {
	P(i+NUM_VAR,j+NUM_VAR) = (spec[i].c)*enth_j*alpham1/omega;     
      }
    }       

  }

} // end  Low_Mach_Number_Preconditioner

/*********************************************************************************************
 * LES3DTF_pState::Low_Mach_Number_Preconditioner_Inverse -- Preconditoner matrix inverse   * 
 *********************************************************************************************/
void LES3DTF_pState::Low_Mach_Number_Preconditioner_Inverse(DenseMatrix &Pinv, 
                                                             const double &deltax, 
                                                             const double &lengthx,
                                                             const double &dTime) {

  double Temp = T();
  double pt = p_t();
  double Rmix = Rtot();
  double enthalpy = h();
  double CP = Cp();
  double c = a();
  double theta = ONE/(Mr2(deltax, lengthx, dTime)*c*c) 
                 + (g()-ONE)/(c*c);  


  double phi = ZERO;
  for(int j=0; j<ns-1; ++j){   
    phi += spec[j].c*(specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() 
		      - CP*Temp*specdata[j].Rs()/Rmix);   
  }

  double AA = pt*(rho*Rmix-theta*pt*CP);
  double BB = Rmix*rho*(theta*pt-rho);
  double EE = HALF*v.sqr() - enthalpy + phi;
  double CC = EE + CP*pt/(rho*Rmix); 
  double DD = HALF*v.sqr() + enthalpy + 5.0*k/3.0;


  Pinv.zero();


  Pinv(0,0) = rho*Rmix/AA*(theta*pt*EE-rho*CC+pt); 
  Pinv(0,1) = -v.x*BB/AA;                          
  Pinv(0,2) = -v.y*BB/AA;                          
  Pinv(0,3) = -v.z*BB/AA;
  Pinv(0,4) = BB/AA;

                               
  Pinv(1,0) = v.x*CC*BB/AA;                        
  Pinv(1,1) = rho*Rmix/AA*(pt+rho*v.x*v.x-theta*pt*(v.x*v.x+CP*pt/(rho*Rmix)));  
  Pinv(1,2) = -v.x*v.y*BB/AA;              
  Pinv(1,3) = -v.x*v.z*BB/AA;                   
  Pinv(1,4) = v.x*BB/AA;


  Pinv(2,0) = v.y*CC*BB/AA;                
  Pinv(2,1) = -v.x*v.y*BB/AA;              
  Pinv(2,2) = rho*Rmix/AA*(pt+v.y*v.y*rho-theta*pt*(v.y*v.y+CP*pt/(rho*Rmix)));  
  Pinv(2,3) = -v.z*v.y*BB/AA;                            
  Pinv(2,4) = v.y*BB/AA;


  Pinv(3,0) = v.z*CC*BB/AA;                
  Pinv(3,1) = -v.x*v.z*BB/AA;              
  Pinv(3,2) = -v.y*v.z*BB/AA;  
  Pinv(3,3) = rho*Rmix/AA*(pt+v.z*v.z*rho-theta*pt*(v.z*v.z+CP*pt/(rho*Rmix)));                            
  Pinv(3,4) = v.z*BB/AA;


  Pinv(4,0) = DD*CC*BB/AA;                          
  Pinv(4,1) = -v.x*DD*BB/AA;                       
  Pinv(4,2) = -v.y*DD*BB/AA;                       
  Pinv(4,3) = -v.z*DD*BB/AA;  
  Pinv(4,4) = rho*Rmix/AA*(theta*pt*(DD-CP*pt/(rho*Rmix))-rho*DD+pt);


  int NUM_VAR = NUM_EULER3D_VAR_SANS_SPECIES;


  // Multispecies
  for (int j=0; j<ns-1; ++j) {       
    double enth_j = specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() 
                    - CP*Temp*specdata[j].Rs()/Rmix;
    Pinv(0,j+NUM_VAR) = -enth_j*BB/AA;      
    Pinv(1,j+NUM_VAR) = -v.x*enth_j*BB/AA;  
    Pinv(2,j+NUM_VAR) = -v.y*enth_j*BB/AA;  
    Pinv(3,j+NUM_VAR) = -v.z*enth_j*BB/AA;   
    Pinv(4,j+NUM_VAR) = -enth_j*BB*DD/AA;

    for (int i=0; i<ns-1; ++i) { 
      if (i==j) { 
	Pinv(i+NUM_VAR,0) = (spec[i].c)*CC*BB/AA;  
	Pinv(i+NUM_VAR,1) = -(spec[i].c)*v.x*BB/AA;
	Pinv(i+NUM_VAR,2) = -(spec[i].c)*v.y*BB/AA;
	Pinv(i+NUM_VAR,3) = -(spec[i].c)*v.z*BB/AA;     
	Pinv(i+NUM_VAR,4) = (spec[i].c)*BB/AA;
	//diagonal
	Pinv(i+NUM_VAR,j+NUM_VAR) = 1.0 - spec[i].c*enth_j*BB/AA; 
      } else {
	Pinv(i+NUM_VAR,j+NUM_VAR) = -spec[i].c*enth_j*BB/AA;     
      }
    }       
  }


} // end Low_Mach_Number_Preconditioner_Inverse

/************************************************************************
 *************** NUMERICAL FLUX FUNCTIONS *******************************
 ************************************************************************/

/**************************************************************************
 * LES3DTF_pState::RoeAverage -- Roe-averaged primitive solution state.  *
 **************************************************************************/
LES3DTF_pState LES3DTF_pState::RoeAverage(const LES3DTF_pState &Wl,
                                            const LES3DTF_pState &Wr) {

   double Hl, Hr, srhol, srhor;
   double Ha, ha;
   LES3DTF_pState Temp;

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

   /* Set the species mass fractions. */
   for (int i=0; i<Wl.ns; i++) {
      Temp.spec[i].c = (srhol*Wl.spec[i].c + srhor*Wr.spec[i].c)/(srhol+srhor);
   } /* endif */

   /* Determine the pressure. */

   Ha = (srhol*Hl+srhor*Hr)/(srhol+srhor);
   ha = Ha - HALF*(sqr(Temp.v.x)+sqr(Temp.v.y)+sqr(Temp.v.z)) - FIVE*Temp.k/THREE;
   Temp.p = Temp.rho*Temp.T(ha)*Temp.Rtot();
   
   /* Return the Roe-averged state. */

   return (Temp);     

}

/**************************************************************************
 * LES3DTF_pState::FluxHLLE_x -- HLLE flux function, x-direction flux.   *
 **************************************************************************/
LES3DTF_cState  LES3DTF_pState::FluxHLLE_x(const LES3DTF_pState &Wl,
                                             const LES3DTF_pState &Wr) {
   
   double wavespeed_l, wavespeed_r;
   LES3DTF_pState Wa, lambdas_l, lambdas_r, lambdas_a;
   LES3DTF_cState Flux, dUrl;

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
   wavespeed_r = max(lambdas_r[5],
                     lambdas_a[5]);
   
   wavespeed_l = min(wavespeed_l, ZERO);
   wavespeed_r = max(wavespeed_r, ZERO); 

   if (wavespeed_l >= ZERO) {
      Flux = Wl.F();
   } else if (wavespeed_r <= ZERO) {
      Flux = Wr.F();
   } else {
      Flux = ((wavespeed_r*Wl.F()-wavespeed_l*Wr.F())
             +(wavespeed_l*wavespeed_r)*dUrl)/
              (wavespeed_r-wavespeed_l);
   } /* endif */
   
   /* Return solution flux. */

   return (Flux);

}

/**************************************************************************
 * LES3DTF_pState::FluxHLLE_n -- HLLE flux function, n-direction flux.   *
 **************************************************************************/
LES3DTF_cState LES3DTF_pState::FluxHLLE_n(const LES3DTF_pState &Wl,
                                            const LES3DTF_pState &Wr,
                                            const Vector3D &norm_dir) {

   // Determine the left and right solution states in the rotate frame.
   LES3DTF_pState Wl_rot(Wl.Rotate(norm_dir));
   LES3DTF_pState Wr_rot(Wr.Rotate(norm_dir));
    
   // Evaluate the intermediate state solution flux in the rotated frame.
   LES3DTF_cState Flux_rot = FluxHLLE_x(Wl_rot, Wr_rot);
 
   // Return numerical flux in un-rotated frame.
   return (Flux_rot.RotateBack(norm_dir));

}

/**************************************************************************
 * LES3DTF_pState::FluxRoe_x -- Roe flux function, x-direction flux.     *
 **************************************************************************/
LES3DTF_cState LES3DTF_pState::FluxRoe_x(const  LES3DTF_pState &Wl,  
                                           const  LES3DTF_pState &Wr) {
   
   LES3DTF_pState Wa, dWrl, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
   LES3DTF_cState Flux;

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
      for (int i=1; i < Wl.num_vars; i++) {
         if (wavespeeds[i] < ZERO) {
            Flux += wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
         } /* endif */
      } /* endfor */
    } else {
      Flux = Wr.F();
      wavespeeds = Wa.lambda_plus(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
      for (int i=1; i < Wl.num_vars; i++) {
         if (wavespeeds[i] > ZERO) {
            Flux -= wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i); 
         } /* endif */
      } /* endfor */
   } /* endif */
    
   /* Return solution flux. */    

   return (Flux);

}

/**************************************************************************
 * LES3DTF_pState::FluxRoe_n -- Roe flux function, n-direction flux.     *
 **************************************************************************/
LES3DTF_cState LES3DTF_pState::FluxRoe_n(const LES3DTF_pState &Wl,
                                           const LES3DTF_pState &Wr,
                                           const Vector3D &norm_dir) {
   
   // Determine the left and right solution states in the rotate frame.
   LES3DTF_pState Wl_rot(Wl.Rotate(norm_dir));
   LES3DTF_pState Wr_rot(Wr.Rotate(norm_dir));
    
   // Evaluate the intermediate state solution flux in the rotated frame.
   LES3DTF_cState Flux_rot = FluxRoe_x(Wl_rot, Wr_rot);
 
   // Return numerical flux in un-rotated frame.
   return (Flux_rot.RotateBack(norm_dir));

}

/**************************************************************************************
 * LES3DTF_pState::AUSMplus_up_x -- AUSMplus_up flux function, x-direction flux.     *
 **************************************************************************************/
LES3DTF_cState LES3DTF_pState::FluxAUSMplus_up_x(const LES3DTF_pState &Wl,
                                                   const LES3DTF_pState &Wr) {
 
  LES3DTF_cState Flux;
  double beta = 0.125, sigma = 0.75, Kp =0.25, Ku = 0.75;
  double alpha, rhohalf, mass_flux_half;
  double ahalf, Ml, Mr, Mplus, Mminus, Mhalf, pplus, pminus, phalf;
  //double al, ar, atilde_l, atilde_r;
 
  // Determine the intermediate state sound speed and density:
  //   al = sqrt(Wl.H()/Wl.rho)*sqrt(TWO*(Wl.g() - ONE)/(Wl.g() + ONE));
  //   ar = sqrt(Wr.H()/Wr.rho)*sqrt(TWO*(Wr.g() - ONE)/(Wr.g() + ONE));
  //   atilde_l = sqr(al)/max(al, Wl.v.x);    
  //   atilde_r = sqr(ar)/max(ar, -Wr.v.x);
  //   ahalf = min(atilde_l, atilde_r);
 
  ahalf = HALF*(Wl.a_t() + Wr.a_t());
  rhohalf = HALF*(Wl.rho + Wr.rho); 
   
 
  // Determine the left and right state Mach numbers based on the
  // intermediate state sound speed:
  Ml = Wl.v.x/ahalf;
  Mr = Wr.v.x/ahalf;
 
  // Determine the reference Mach number, scaling function and coefficient
  double M2_bar, M2_ref, fa;
  M2_bar = (Wl.v.x*Wl.v.x + Wr.v.x*Wr.v.x)/(TWO*ahalf*ahalf);
  M2_ref = min(ONE, max(M2_bar, Wl.Mref*Wl.Mref));
  if (M2_ref > ONE || M2_ref < 0.0) cout << "\nM2_ref out of range in AUSM+-up.";
  //fa = sqrt(M2_ref)*(TWO - sqrt(M2_ref));
  fa = sqrt(sqr(ONE - M2_ref)*M2_bar + FOUR*M2_ref)/(ONE + M2_ref);
  if (fa > ONE || fa <= ZERO) cout << "\nfa out of range in AUSM+-up.";
  alpha = (3.0/16.0)*(-4.0 + 5.0*fa*fa);
  if (alpha < (-3.0/4.0)  ||  alpha > (3.0/16.0)) cout << "\nalpha out of range in AUSM+-up.";
 
 
  // Determine the left state split Mach number:
  if (fabs(Ml) >= ONE) {
    Mplus = 0.5*(Ml+fabs(Ml));
    pplus = 0.5*(Ml+fabs(Ml))/Ml;
  } else {
    Mplus = 0.25*sqr(Ml+1.0) * (1.0 - 16.0*beta*(-0.25*sqr(Ml-1.0)));
    pplus = 0.25*sqr(Ml+1.0) * ((2.0 - Ml) - 16.0*alpha*Ml*(-0.25*sqr(Ml-1.0)));
  } /* endif */
  
  // Determine the right state split Mach number:
  if (fabs(Mr) >= ONE) {
    Mminus = 0.5*(Mr-fabs(Mr));
    pminus = 0.5*(Mr-fabs(Mr))/Mr;        
  } else {
    Mminus = -0.25*sqr(Mr-1.0) * (1.0 + 16.0*beta*0.25*sqr(Mr+1.0));
    pminus = -0.25*sqr(Mr-1.0) * ((-2.0 - Mr) + 16.0*alpha*Mr*0.25*sqr(Mr+1.0));
  } /* endif */

  // Determine the left state split Mach number:
//   if (fabs(Ml) >= ONE) {
//     Mplus = Mplus_1(Ml);
//     pplus = Mplus_1(Ml)/Ml;
//   } else {
//     Mplus = Mplus_2(Ml) * (1.0 - 16.0*beta*Mminus_2(Ml));
//     pplus = Mplus_2(Ml) * ((2.0 - Ml) - 16.0*alpha*Ml*Mminus_2(Ml));
//   }
  
  // Determine the right state split Mach number:
//   if (fabs(Mr) >= ONE) {
//     Mminus = Mminus_1(Mr);
//     pminus = Mminus_1(Mr)/Mr;        
//   } else {
//     Mminus = Mminus_2(Mr) * (1.0 + 16.0*beta*Mplus_2(Mr));
//     pminus = Mminus_2(Mr) * ((-2.0 - Mr) + 16.0*alpha*Mr*Mplus_2(Mr));
//   } 
 
   // Determine the intermediate state Mach number, pressure and mass flux:
  Mhalf = Mplus + Mminus -
          (Kp/fa)*max((ONE - sigma*M2_bar), ZERO)*(Wr.p - Wl.p)/(rhohalf*ahalf*ahalf);
 
  phalf = pplus*Wl.p + pminus*Wr.p -
          Ku*pplus*pminus*TWO*rhohalf*(fa*ahalf)*(Wr.v.x - Wl.v.x);
 
  mass_flux_half = (Mhalf > ZERO) ? ahalf*Mhalf*Wl.rho : ahalf*Mhalf*Wr.rho; 
  
  // Determine the intermediate state convective solution flux:
  if (mass_flux_half > ZERO) {
    Flux.rho = ONE;
    Flux.rhov.x = Wl.v.x; 
    Flux.rhov.y = Wl.v.y; 
    Flux.E = Wl.H()/Wl.rho;
    Flux.rhok = Wl.k;
    for(int i=0; i<Wl.ns; ++i){
      Flux.rhospec[i].c = Wl.spec[i].c;
    }
  } else {
    Flux.rho = ONE;
    Flux.rhov.x = Wr.v.x; 
    Flux.rhov.y = Wr.v.y; 
    Flux.E = Wr.H()/Wr.rho;
    Flux.rhok = Wr.k;
    for(int i=0; i<Wr.ns; ++i){
      Flux.rhospec[i].c = Wr.spec[i].c;
    }
  } /* endif */

  Flux = mass_flux_half*Flux;
 
  // Add the pressure contribution to the intermediate state solution flux:
  Flux[2] += phalf;
  
  // Return solution flux.
  return Flux;

}

/******************************************************************************************
 * LES3DTF_pState::FluxAUSMplus_up_n -- AUSMplus_up flux function, n-direction flux.     *
 ******************************************************************************************/
LES3DTF_cState LES3DTF_pState::FluxAUSMplus_up_n(const LES3DTF_pState &Wl,
                                                   const LES3DTF_pState &Wr,
                                                   const Vector3D &norm_dir) {
   
   // Determine the left and right solution states in the rotate frame.
   LES3DTF_pState Wl_rot(Wl.Rotate(norm_dir));
   LES3DTF_pState Wr_rot(Wr.Rotate(norm_dir));
    
   // Evaluate the intermediate state solution flux in the rotated frame.
   LES3DTF_cState Flux_rot = FluxAUSMplus_up_x(Wl_rot, Wr_rot);
 
   // Return numerical flux in un-rotated frame.
   return (Flux_rot.RotateBack(norm_dir));

}

/**************************************************************************
 * LES3DTF_pState::lambda_minus -- Negative wave speeds determined using *
 *                                  Harten entropy fix.                   *
 **************************************************************************/
LES3DTF_pState LES3DTF_pState::lambda_minus(const LES3DTF_pState &lambdas_a,
                                              const LES3DTF_pState &lambdas_l,
                                              const LES3DTF_pState &lambdas_r) {

   LES3DTF_pState W_temp;

   W_temp.rho = HartenFixNeg(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
   W_temp.v.x = HALF*(lambdas_a[2]-fabs(lambdas_a[2]));
   W_temp.v.y = HALF*(lambdas_a[3]-fabs(lambdas_a[3]));
   W_temp.v.z = HALF*(lambdas_a[4]-fabs(lambdas_a[4]));
   W_temp.p = HartenFixNeg(lambdas_a[5],lambdas_l[5],lambdas_r[5]);
   W_temp.k = HALF*(lambdas_a[8] -fabs(lambdas_a[8]));
   
   for (int i = (NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA + 1); i <= W_temp.num_vars; i++) {
      W_temp.spec[i-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].c = 
         HALF*(lambdas_a[i] - fabs(lambdas_a[i]));
   } /* endfor */

   return (W_temp);

}

/**************************************************************************
 * LES3DTF_pState::lambda_plus -- Positive wave speeds determined using  *
 *                                 Harten entropy fix.                    *
 **************************************************************************/
LES3DTF_pState LES3DTF_pState::lambda_plus(const LES3DTF_pState &lambdas_a,
                                             const LES3DTF_pState &lambdas_l,
                                             const LES3DTF_pState &lambdas_r) {
   
   LES3DTF_pState W_temp;

   W_temp.rho = HartenFixPos(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
   W_temp.v.x = HALF*(lambdas_a[2]+fabs(lambdas_a[2]));
   W_temp.v.y = HALF*(lambdas_a[3]+fabs(lambdas_a[3]));
   W_temp.v.z = HALF*(lambdas_a[4]+fabs(lambdas_a[4]));
   W_temp.p = HartenFixPos(lambdas_a[5],lambdas_l[5],lambdas_r[5]);
   W_temp.k = HALF*(lambdas_a[8] + fabs(lambdas_a[8]));

   for (int i = (NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA +1); i <= W_temp.num_vars; i++) {
      W_temp.spec[i-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].c = 
         HALF*(lambdas_a[i] + fabs(lambdas_a[i]));
   } /* endfor */

   return (W_temp);

}

/*************************************************************************************************
 * LES3DTF_pState::HLLE_wavespeeds -- Returns the lambda plus and lambda minus wave speeds for  *
 *                                     rotated Riemann problem aligned with norm_dir given       * 
 *                                     unroated solution states Wl and Wr.                       *
 ************************************************************************************************/
Vector2D LES3DTF_pState::HLLE_wavespeeds(const LES3DTF_pState &Wl,
                                          const LES3DTF_pState &Wr,
                                          const Vector3D &norm_dir) {

    Vector2D wavespeed;
    LES3DTF_pState Wa_n, lambdas_l, lambdas_r, lambdas_a;  //Lots of TEMPS

    /* Determine the left and right states. */

    LES3DTF_pState Wl_rotated(Wl.Rotate(norm_dir));
    LES3DTF_pState Wr_rotated(Wr.Rotate(norm_dir));

    /* Evaluate the Roe-average primitive solution state. */                           

    Wa_n = Wa_n.RoeAverage(Wl_rotated, Wr_rotated);
    
    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl_rotated.lambda_x();
    lambdas_r = Wr_rotated.lambda_x();
    lambdas_a = Wa_n.lambda_x();

    /* Determine the intermediate state flux. */

    wavespeed.x = min(lambdas_l[1],
                      lambdas_a[1]);
    wavespeed.y = max(lambdas_r[5],
                      lambdas_a[5]);
 
    wavespeed.x = min(wavespeed.x, ZERO); //lambda minus
    wavespeed.y = max(wavespeed.y, ZERO); //lambda plus 

    /* Return the required wave speeds. */

    return (wavespeed);

}

/**************************************************************************
 * LES3DTF_pState::Rotate -- Returns a rotated primitive state aligned   *
 *                            with a local x-axis in the norm_dir.        *
 **************************************************************************/
LES3DTF_pState LES3DTF_pState::Rotate(const Vector3D &norm_dir) const {

   double cos_psi, sin_psi, cos_phi, sin_phi, cos_theta, sin_theta;
   cos_phi = ONE;
   sin_phi = ZERO;
   if (fabs(fabs(norm_dir.x)-ONE) < TOLER) {
      cos_psi = norm_dir.x/fabs(norm_dir.x);
      sin_psi = ZERO;
      cos_theta = ONE;
      sin_theta = ZERO;
   } else {
      cos_psi = norm_dir.x;
      sin_psi = sqrt(norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z);
      cos_theta = norm_dir.y/sqrt(norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z);
      sin_theta = norm_dir.z/sqrt(norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z);
   } /* endif */

   return LES3DTF_pState(rho,
                          (cos_psi*cos_phi-cos_theta*sin_phi*sin_psi)*v.x + 
                          (cos_psi*sin_phi+cos_theta*cos_phi*sin_psi)*v.y + 
                          (sin_psi*sin_theta)*v.z,
                          (-sin_psi*cos_phi-cos_theta*sin_phi*cos_psi)*v.x + 
                          (-sin_psi*sin_phi+cos_theta*cos_phi*cos_psi)*v.y + 
                          (cos_psi*sin_theta)*v.z,
                          (sin_theta*sin_phi)*v.x + 
                          (-sin_theta*cos_phi)*v.y + 
                          (cos_theta)*v.z,
                          p,
                          k);

}

/**************************************************************************
 * LES3DTF_pState::Rotate -- Returns an un-rotated primitive state       *
 *                            re-alinged from the x-axis of the global    *
 *                            problem.                                    *
 **************************************************************************/
LES3DTF_pState LES3DTF_pState::RotateBack(const Vector3D &norm_dir) const {

   double cos_psi, sin_psi, cos_phi, sin_phi, cos_theta, sin_theta;
   cos_phi = ONE;
   sin_phi = ZERO;
   if (fabs(fabs(norm_dir.x)-ONE) < TOLER) {
      cos_psi = norm_dir.x/fabs(norm_dir.x);
      sin_psi = ZERO;
      cos_theta = ONE;
      sin_theta = ZERO;
   } else {
      cos_psi = norm_dir.x;
      sin_psi = sqrt(norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z);
      cos_theta = norm_dir.y/sqrt(norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z);
      sin_theta = norm_dir.z/sqrt(norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z);
   } /* endif */

   return LES3DTF_pState(rho,
                          (cos_psi*cos_phi-cos_theta*sin_phi*sin_psi)*v.x + 
                          (-sin_psi*cos_phi-cos_theta*sin_phi*cos_psi)*v.y + 
                          (sin_theta*sin_phi)*v.z,
                          (cos_psi*sin_phi+cos_theta*cos_phi*sin_psi)*v.x + 
                          (-sin_psi*sin_phi+cos_theta*cos_phi*cos_psi)*v.y + 
                          (-sin_theta*cos_phi)*v.z,
                          (sin_theta*sin_psi)*v.x + 
                          (sin_theta*cos_psi)*v.y + 
                          (cos_theta)*v.z,
                          p,
                          k);

}

/************************************************************************
 *************** NUMERICAL EVALUATION OF VISCOUS FLUXES *****************
 ************************************************************************/

/**********************************************************************
 * LES3DTF_pState::FluxViscous_n -- Viscous flux (n-direction).      *
/**********************************************************************/
LES3DTF_cState LES3DTF_pState::FluxViscous_n(const LES3DTF_pState &Wl,
                                               const LES3DTF_pState &Wr,
                                               const LES3DTF_pState &Wc,
                                               const LES3DTF_pState &Wc_Neigbor,
                                               const LES3DTF_pState &dWdx,
                                               const LES3DTF_pState &dWdy,
                                               const LES3DTF_pState &dWdz,
                                               const LES3DTF_pState &dWdx_Neigbor,
                                               const LES3DTF_pState &dWdy_Neigbor,
                                               const LES3DTF_pState &dWdz_Neigbor,
                                               const Vector3D &norm,
                                               const Vector3D &ts, 
                                               const double &deltad, 
                                               const double &Volume, 
                                               const double &Volume_Neigbor, 
                                               const int Flow_Type) {
   
   // construct the gradients on the cell interface (surface) 
   // based on Hybrid Average Gradient-Diamond-Path Approach
   // Grad Phi = (Phi_c_neigbor - Phi_c)/ds *norm/(norm . ts)
   //            + [bar (Grad Phi) - bar (Grad Phi) . ts norm/(norm.ts)]
  
   // weighted factor based on volume
   double alpha = Volume/(Volume + Volume_Neigbor);
   
   LES3DTF_pState dWdx_Weighted, dWdy_Weighted, dWdz_Weighted, 
                   dWdx_face, dWdy_face, dWdz_face, Grad_middle_term;

   LES3DTF_pState W_face;
   
   dWdx_Weighted = alpha*dWdx + (1.0 - alpha)*dWdx_Neigbor;
   dWdy_Weighted = alpha*dWdy + (1.0 - alpha)*dWdy_Neigbor;
   dWdz_Weighted = alpha*dWdz + (1.0 - alpha)*dWdz_Neigbor;

//    dWdx_face = dWdx_Weighted ;
//    dWdy_face = dWdy_Weighted ;
//    dWdz_face = dWdz_Weighted ; 
   
   // Evaluate a weighted term for solution gradients  
   Grad_middle_term = dWdx_Weighted*ts.x + dWdy_Weighted*ts.y + dWdz_Weighted*ts.z;
      
   // Evaluate gradients of primitive variables on the face
   dWdx_face = (Wc_Neigbor - Wc)/deltad *norm.x/dot(norm, ts) + 
               (dWdx_Weighted -  Grad_middle_term*norm.x/dot(norm, ts));
   dWdy_face = (Wc_Neigbor - Wc)/deltad *norm.y/dot(norm, ts) + 
               (dWdy_Weighted -  Grad_middle_term*norm.y/dot(norm, ts));
   dWdz_face = (Wc_Neigbor - Wc)/deltad *norm.z/dot(norm, ts) + 
               (dWdz_Weighted -  Grad_middle_term*norm.z/dot(norm, ts));

   // Determine face solution state.   
   W_face = HALF*(Wl + Wr);
   //W_face = HALF*(Wc + Wc_Neigbor);

   // Evaluate viscous flux
   if (fabs(norm.y) < TOLER && fabs(norm.z) < TOLER) {
      return (W_face.Fvx(dWdx_face, dWdy_face, dWdz_face, Flow_Type, Volume_Neigbor)*norm.x);

   } else if (fabs(norm.x) < TOLER && fabs(norm.z) < TOLER) {
      return (W_face.Fvy(dWdx_face, dWdy_face, dWdz_face, Flow_Type, Volume_Neigbor)*norm.y);

   } else if (fabs(norm.x) < TOLER && fabs(norm.y) < TOLER) {
      return (W_face.Fvz(dWdx_face, dWdy_face, dWdz_face, Flow_Type, Volume_Neigbor)*norm.z);

   } else {
      return (W_face.Fvx(dWdx_face, dWdy_face, dWdz_face, Flow_Type, Volume_Neigbor)*norm.x +
              W_face.Fvy(dWdx_face, dWdy_face, dWdz_face, Flow_Type, Volume_Neigbor)*norm.y +
              W_face.Fvz(dWdx_face, dWdy_face, dWdz_face, Flow_Type, Volume_Neigbor)*norm.z);
   } /* endif */

}

/************************************************************************
 ******************** BOUNDARY CONDITIONS *******************************
 ************************************************************************/

/**************************************************************************
 * LES3DTF_pState::Reflect -- Return reflected solution state.           *
 **************************************************************************/
LES3DTF_pState LES3DTF_pState::Reflect(const LES3DTF_pState &W,
                                         const Vector3D &norm_dir) {
   
   Vector3D ur_norm, ur_tang, vr_tot;
   LES3DTF_pState Temp; Temp.Copy(W);

   ur_norm = dot(W.v, norm_dir)*norm_dir;
   ur_tang = W.v - ur_norm;

   ur_norm = -ur_norm;
   vr_tot = ur_norm + ur_tang;

   Temp.v = vr_tot;
   
   return (Temp);
       
}

/**************************************************************************
 * LES3DTF_pState::MovingWall -- Return moving wall boundary state.      *
 **************************************************************************/
LES3DTF_pState LES3DTF_pState::MovingWall(const LES3DTF_pState &Win,
                                            const LES3DTF_pState &Wout,
                                            const Vector3D &norm_dir, 
                                            const Vector3D &wall_velocity,
                                            const Vector3D &pressure_gradient,
                                            const int &TEMPERATURE_BC_FLAG) {
   
   LES3DTF_pState Temp; Temp.Copy(Win);
   
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
        Temp.rho = Wout.p/(Win.Rtot()*Wout.T());
     } else {
        Temp.rho = Temp.p/(Win.Rtot()*Wout.T());
     } /* endif */
  } /* endif */

  return (Temp);

}

/**************************************************************************
 * LES3DTF_pState::NoSlip -- Return no-slip wall boundary state.         *
 **************************************************************************/
LES3DTF_pState LES3DTF_pState::NoSlip(const LES3DTF_pState &Win,
                                        const LES3DTF_pState &Wout,
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

/*********************************************************************************
 * LES3DTF_pState::Enstrophy                                                    *
 *********************************************************************************/
double LES3DTF_pState::Enstrophy(const LES3DTF_pState &dWdx, 
                                  const LES3DTF_pState &dWdy, 
                                  const LES3DTF_pState &dWdz) const{
   return sqr(vorticity(dWdx,dWdy,dWdz));
}

/*********************************************************************************
 * LES3DTF_pState::abs_strain_rate -- Absolute value of strain rate tensor      *
 *********************************************************************************/
double LES3DTF_pState::abs_strain_rate(const LES3DTF_pState &dWdx, 
                                        const LES3DTF_pState &dWdy, 
                                        const LES3DTF_pState &dWdz) const{

   Tensor3D strain_rate_tensor = strain_rate(dWdx, dWdy, dWdz);
   // S[i,j]*S[i,j]
   double SS = sqr(strain_rate_tensor.xx) + sqr(strain_rate_tensor.yy) + 
               sqr(strain_rate_tensor.zz) +
               TWO*(sqr(strain_rate_tensor.xy) + sqr(strain_rate_tensor.yz) + 
                    sqr(strain_rate_tensor.xz));
   // sqrt(2*S*S)
   return sqrt(TWO*SS);

}

/*********************************************************************************
 * LES3DTF_pState::SFS_Kinetic_Energy -- Subfilter scale kinetic energy     *
 *********************************************************************************/
double LES3DTF_pState::SFS_Kinetic_Energy(const LES3DTF_pState &dWdx,
                                               const LES3DTF_pState &dWdy,
                                               const LES3DTF_pState &dWdz,
                                               const int Flow_Type,
                                               const double &Volume) {

  if ( Flow_Type == FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY ) {
    double CI = 0.005;
    return (CI*sqr(filter_width(Volume)*abs_strain_rate(dWdx,dWdy,dWdz)));
  } else if ( Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K ) {
    return (k);
  } /* endif */

}

/*****************************************************************
 * LES3DTF_pState::K_equ_sources -- Source term for k-equation  *
 *****************************************************************/
double LES3DTF_pState::K_equ_sources(const LES3DTF_pState &dWdx,
                                      const LES3DTF_pState &dWdy,
                                      const LES3DTF_pState &dWdz,
                                      const int Flow_Type,
                                      const double &Volume) {

     double production, dissipation, source, mu_t;
     Tensor3D subfilter_stress;
     double Ceps = 0.845;
 
     subfilter_stress = tau_t(dWdx, dWdy, dWdz, Flow_Type, Volume);
  
       production = subfilter_stress.xx*dWdx.v.x + 
       subfilter_stress.xy*(dWdy.v.x + dWdx.v.y) +
       subfilter_stress.yy*dWdy.v.y +
       subfilter_stress.xz*(dWdz.v.x + dWdx.v.z) +
       subfilter_stress.yz*(dWdz.v.y + dWdy.v.z) +
       subfilter_stress.zz*dWdz.v.z;

       dissipation = Ceps*rho*pow(k, 1.5)/filter_width(Volume);
       source = production - dissipation;
       return(source);

}

/*****************************************************************
 * LES3DTF_pState::Sturbulene -- Turbulence model source terms. *
 *****************************************************************/
LES3DTF_cState LES3DTF_pState::Sturbchem(LES3DTF_pState &Wc,
                                           const LES3DTF_pState &dWdx,
                                           const LES3DTF_pState &dWdy,
                                           const LES3DTF_pState &dWdz,
                                           const int Flow_Type,
                                           const double &Volume) {
   
  LES3DTF_cState Temp;

      if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {
         // k-equation Source Term
         Temp.rhok = Wc.K_equ_sources(dWdx,dWdy,dWdz,Flow_Type,Volume);
  } /* endif */

  return (Temp);

} 

/*****************************************************************************************
 * LESS3DFsd_cState member functions                                                     *
 *****************************************************************************************/

/***************************************************************************************
 * LES3DTF_cState::Copy -- Makes a copy of solution state vector.                     *
 ***************************************************************************************/
void LES3DTF_cState::Copy(const LES3DTF_cState &U) {
  rho = U.rho;
  rhov = U.rhov; 
  E = U.E;  
  rhok = U.rhok; 
  for (int i = 0; i < ns; i++) { 
     rhospec[i] = U.rhospec[i];
  } /* endfor */
  flame = U.flame;
}

/*********************************************************************************************
 * LES3DTF_cState::Realizable_Solution_Check -- Check physical validity of solution state.  *
 *********************************************************************************************/
bool LES3DTF_cState::Realizable_Solution_Check(void) {
   if (rho <= ZERO || !negative_speccheck() || es() <= ZERO || rhok < ZERO ) {    
      cout << "\n " << CFFC_Name() 
           << " ERROR: Conservative solution state has a negative density, energy, mass fractions,"
           << " and/or turbulent kinetic energy.\n";
      return false;
   } else {
      return true;
   } /* endif */
} 

/*************************************************************************
 * LES3DTF_cState::e -- Return mixture absolute internal energy.        *
 *************************************************************************/
double LES3DTF_cState::e(void) const {
  return ((E - HALF*rhov.sqr()/rho-rhok)/rho);
}

/**************************************************************************
 * LES3DTF_cState::es -- Return sensible internal energy.                *
 **************************************************************************/
double LES3DTF_cState::es(void) const {
  return ((E - HALF*rhov.sqr()/rho-rhok)/rho-HeatofFormation());
}

/**************************************************************************
 * LES3DTF_cState::h -- Return mixture absolute internal enthalpy.       *
 **************************************************************************/
double LES3DTF_cState::h(void) const {
  return (e()+p()/rho);
}

double LES3DTF_cState::h(const double &Temp) const {
  return (Euler3D_ThermallyPerfect_cState::h(Temp));
}

/**************************************************************************
 * LES3DTF_cState::hs -- Return sensible internal enthalpy.              *
 **************************************************************************/
double LES3DTF_cState::hs(void) const {
  return (es()+p()/rho);
}

double LES3DTF_cState::hs(const double &Temp) const {
  return (Euler3D_ThermallyPerfect_cState::hs());
}

/*******************************************************************
 * LES3DTF_cState::T -- Return mixture temperature.               *
 *******************************************************************/
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
double LES3DTF_cState::T(void) const {

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
         cout << "\nTemperature didn't converge in LES3DTF_cState::T(void)";
         cout << " with polytopic Tguess " << Tguess << ", or lower than Tmin " 
              << low_temp_range << " using " << T;
      } /* endif */
   } /* endif */

   return T;

} 

/*******************************************************************
 * LES3DTF_cState::p_t -- Return turbulence modified pressure.     *
 *******************************************************************/
double LES3DTF_cState::p_t(void) const { 
  return (p() + 2.0*rhok/3.0);
}

/*******************************************************************
 * LES3DTF_cState::a_t -- Return mixture sound speed (including   *
 *                         turbulent kinetic energy).              *
 *******************************************************************/
double LES3DTF_cState::a_t(void) const {
   double aa = sqr(a());
   aa += (TWO/THREE)*(rhok/rho)*g();
   return sqrt(aa);
}

/*******************************************************************
 * LES3DTF_cState::k -- Return turbulent kinetic energy.          *
 *******************************************************************/
double LES3DTF_cState::k() const {
   return (rhok/rho);
}

/*******************************************************************
 * LES3DTF_cState::mu_t -- Return eddy (turbulent) viscosity.     *
 *******************************************************************/
double LES3DTF_cState::mu_t(const LES3DTF_pState &dWdx,
 			     const LES3DTF_pState &dWdy,
			     const LES3DTF_pState &dWdz,
                             const int Flow_Type, 
                             const double &Volume) {

  double filter = filter_width(Volume);
  if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY) {
    double Cs = 0.18;
    return(rho*sqr(Cs*filter)*abs_strain_rate(dWdx,dWdy,dWdz));
  } else if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {
    double Cv = 0.086;
    return(rho*Cv*sqrt(k())*filter);
  }

}

/**************************************************************************
 * LES3DTF_cState::kappa_t -- Return turbulent thermal conductivity.     *
 **************************************************************************/
double LES3DTF_cState::kappa_t(const LES3DTF_pState &dWdx,
 			        const LES3DTF_pState &dWdy,
			        const LES3DTF_pState &dWdz,
                                const int Flow_Type, 
                                const double &Volume) {

  return (mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume)*Cp()/Pr_t()); 

}

/********************************************************************************
 * LES3DTF_cState::Ds_t -- Return species turbulent diffusion coefficient.     *
 ********************************************************************************/
double LES3DTF_cState::Ds_t(const int i,
                             const LES3DTF_pState &dWdx,
 			     const LES3DTF_pState &dWdy,
			     const LES3DTF_pState &dWdz,
                             const int Flow_Type, 
                             const double &Volume) {

  return (mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume)/(rho*Sc_t())); 

}

double LES3DTF_cState::Ds_t(const int i,
                             const double &mu_t_temp) {

  return (mu_t_temp/(rho*Sc_t())); 

}

/*******************************************************************
 * LES3DTF_cState::Pr_t -- Return turbulent Prandtl number.       *
 *******************************************************************/
double LES3DTF_cState::Pr_t(void) {
   return (0.9); 
}

/*******************************************************************
 * LES3DTF_cState::Sc_t -- Return turbulent Schmidt number.       *
 *******************************************************************/
double LES3DTF_cState::Sc_t(void) {
   return (1.0);
}

/*******************************************************************
 * LES3DTF_cState::Le_t -- Return turbulent Lewis number.         *
 *******************************************************************/
double LES3DTF_cState::Le_t(void) {
  return (Sc_t()/Pr_t());
}

/*************************************************************************
 * LES3DTF_cState::filter_width -- LES characteristic filter width      *
 *************************************************************************/
double LES3DTF_cState::filter_width(const double &Volume) const {
  return (2.0*pow(Volume,1.0/3.0)); 
}

/*********************************************************************************
 * LES3DTF_cState::abs_strain_rate -- Absolute value of strain rate tensor      *
 *********************************************************************************/
double LES3DTF_cState::abs_strain_rate(const LES3DTF_pState &dWdx, 
                                        const LES3DTF_pState &dWdy, 
                                        const LES3DTF_pState &dWdz) const{

   Tensor3D strain_rate; 
   strain_rate.zero();
   strain_rate.xx = dWdx.v.x;
   strain_rate.yy = dWdy.v.y;
   strain_rate.zz = dWdz.v.z;
   strain_rate.xy = 0.5*(dWdy.v.x + dWdx.v.y);
   strain_rate.yz = 0.5*(dWdz.v.y + dWdy.v.z);
   strain_rate.xz = 0.5*(dWdx.v.z + dWdz.v.x);
   // S[i,j]*S[i,j]
   double SS = sqr(strain_rate.xx) + sqr(strain_rate.yy) + sqr(strain_rate.zz) +
               2.0*(sqr(strain_rate.xy) + sqr(strain_rate.yz) + sqr(strain_rate.xz));
   // sqrt(2*S*S)
   return sqrt(2.0*SS);

}

/*****************************************************************
 * LES3DTF_pState::Sw -- Chemical Reaction Rate Source Terms.    *
 * Using the Reaction class to get the source terms for the      * 
 * specific "Reaction_set".                                      * 
 *****************************************************************/
LES3DTF_cState  LES3DTF_pState::Sw(const int &REACT_SET_FLAG) {
 
   LES3DTF_cState NEW; 
   LES3DTF_pState Temp; 
   NEW.Vacuum();
//    bool test = negative_speccheck();

   //Adds concentration rate of change for species 1->N
   if( REACT_SET_FLAG != NO_REACTIONS) {
      React.omega(NEW,Temp);  
   }

   return (NEW.flame.WF/NEW.flame.TF)*NEW;
}


/**************************************************************** 
 * LES3DTF_pState::dSwdU_max_diagonal -- Max Diagonal of        *
 *                                       Jacobian for CFL       * 
 ****************************************************************/
// double LES3DTF_pState::dSwdU_max_diagonal(void) const {
   
//    double max_diagonal = ONE;
//    DenseMatrix dSwdU(num_vars-1, num_vars -1);
//    dSwdU.zero();

// //    bool test = negative_speccheck();

//    React.dSwdU(dSwdU,*this,true);
   
//    for(int i=0; i < num_vars -1; ++i){
//      max_diagonal = max(max_diagonal, fabs( (flame.WF/flame.TF)*dSwdU(i,i)) );
//    }

//    return max_diagonal;
// }

/*******************************************************************
 * LES3DTF_cState::W -- Return primitive solution state vector.   *
 *******************************************************************/
LES3DTF_pState LES3DTF_cState::W(void){
   LES3DTF_pState Temp;
   Temp.rho = rho;
   Temp.v = v();  
   Temp.p = p();
   Temp.k = k();  
   for (int i = 0; i < ns; i++) {
      Temp.spec[i] = rhospec[i]/rho;
   } /* endfor */
   Temp.flame = flame;
   return Temp;
}

LES3DTF_pState LES3DTF_cState::W(void) const{
   LES3DTF_pState Temp;
   Temp.rho = rho;
   Temp.v = v();  
   Temp.p = p();
   Temp.k = k();  
   for (int i = 0; i < ns; i++) {
      Temp.spec[i] = rhospec[i]/rho;
   } /* endfor */
   Temp.flame = flame;
   return Temp;
}

LES3DTF_pState LES3DTF_cState::W(const LES3DTF_cState &U) const {
   LES3DTF_pState Temp;
   Temp.rho = U.rho;
   Temp.v = U.v();
   Temp.p = U.p();
   Temp.k = U.k();
   for (int i = 0; i < ns; i++) {
      Temp.spec[i] = rhospec[i]/rho;
   } /* endfor */
   Temp.flame = U.flame;
   return Temp;
}


/**************************************************************************
 * LES3DTF_cState::Rotate -- Returns a rotated primitive state aligned   *
 *                            with a local x-axis in the norm_dir.        *
 **************************************************************************/
LES3DTF_cState LES3DTF_cState::Rotate(const Vector3D &norm_dir) const {

  // for a 3D unit normal rotated to align with the x-axis
  double Ct = norm_dir.x;  //cos_angle
  double St = sqrt( norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z); //sin_angle
  Vector3D rt(0,norm_dir.z,-norm_dir.y);  //rotation axis  

  return LES3DTF_cState(rho,
                         rhov.x*Ct - rhov.y*rt.z*St + rhov.z*rt.y*St,
                         rhov.x*rt.z*St + rhov.y*(rt.y*rt.y*(ONE-Ct)+Ct) + 
                         rhov.z*(rt.y*rt.z*(ONE-Ct)),
                         -rhov.x*rt.y*St +  rhov.y*(rt.y*rt.z*(ONE-Ct)) + 
                         rhov.z*(rt.z*rt.z*(ONE-Ct)+Ct),
                         E,
			rhok,
                        rhospec);

}

/**************************************************************************
 * LES3DTF_cState::Rotate -- Returns an un-rotated primitive state       *
 *                            re-alinged from the x-axis of the global    *
 *                            problem.                                    *
 **************************************************************************/
LES3DTF_cState LES3DTF_cState::RotateBack(const Vector3D &norm_dir) const {

  // for a 3D unit normal rotated to align with the x-axis
  double Ct = norm_dir.x;  //cos_angle
  double St = sqrt( norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z); //sin_angle
  Vector3D rt(0,norm_dir.z,-norm_dir.y);  //rotation axis  

  return LES3DTF_cState(rho,
                         rhov.x*Ct + rhov.y*rt.z*St - rhov.z*rt.y*St,
                         -rhov.x*rt.z*St + rhov.y*(rt.y*rt.y*(ONE-Ct)+Ct) + 
                         rhov.z*(rt.y*rt.z*(ONE-Ct)),
                         + rhov.x*rt.y*St +  rhov.y*(rt.y*rt.z*(ONE-Ct)) + 
                         rhov.z*(rt.z*rt.z*(ONE-Ct)+Ct),
                         E,
			rhok,
                        rhospec);

}

