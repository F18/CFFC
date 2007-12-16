/*! \file  LES3DFsdState.cc
\brief  Definition of member functions for 3D Favre-filtered Navier-Stokes solution classes
        associated with solution of premixed compressible turbulent combusting flows of a 
        thermally perfect using a large eddy simulation (LES) technique in conjunction with 
        a flame surface density (FSD) subfilter scale model.
*/

/* Include LES3DFsdState header file. */

#ifndef _LES3DFSD_STATE_INCLUDED
#include "LES3DFsdState.h"
#endif // LES3DFSD_STATE_INCLUDED   

/***************************************************************************************
 * LES3DFsd_pState -- Create storage and assign various static values.                 *
 ***************************************************************************************/
double LES3DFsd_pState::Mref = 0.1;
double LES3DFsd_pState::_fuel_equivalence_ratio=1.0;
double LES3DFsd_pState::_unburnt_fuel_mass_fraction=0.05518;
double LES3DFsd_pState::_reactants_density=1.13;
double LES3DFsd_pState::_laminar_flame_speed=0.3738;
double LES3DFsd_pState::_laminar_flame_thickness=4.4E-04;
double LES3DFsd_pState::_adiabatic_flame_temperature=2218.0;

/***************************************************************************************
 * LES3DFsd_cState -- Create storage and assign various static values.                 *
 ***************************************************************************************/
double LES3DFsd_cState::Mref = 0.1;
double LES3DFsd_cState::_fuel_equivalence_ratio=1.0;
double LES3DFsd_cState::_unburnt_fuel_mass_fraction=0.05518;
double LES3DFsd_cState::_reactants_density=1.13;
double LES3DFsd_cState::_laminar_flame_speed=0.3738;
double LES3DFsd_cState::_laminar_flame_thickness=4.4E-04;
double LES3DFsd_cState::_adiabatic_flame_temperature=2218.0;

/*****************************************************************************************
 * LESS3DFsd_pState member functions                                                     *
 *****************************************************************************************/

/***************************************************************************************
 * LES3DFsd_pState::Copy -- Makes a copy of solution state vector.                     *
 ***************************************************************************************/
void LES3DFsd_pState::Copy(const LES3DFsd_pState &W) {
  rho = W.rho;
  v = W.v; 
  p = W.p;  
  C = W.C; 
  Fsd = W.Fsd; 
  k = W.k; 
}

/*********************************************************************************************
 * LES3DFsd_pState::Realizable_Solution_Check -- Check physical validity of solution state.  *
 *********************************************************************************************/
bool LES3DFsd_pState::Realizable_Solution_Check(void) {
   if (rho <= ZERO || p <= ZERO || C <= ZERO || Fsd <= ZERO || k <= ZERO ) {    
      cout << "\n " << CFFC_Name() 
           << " ERROR: Primitive solution state has a negative density, pressure, progress variable,"
           << " FSD, and/or turbulent kinetic energy.\n";
      return false;
   } else {
      return true;
   } /* endif */
} 

/*****************************************************************************************
 * LES3DFsd_pState::E -- Return total energy of the mixture.                             *
 *****************************************************************************************/
double LES3DFsd_pState::E(void) const {   
   return (rho*(e()+HALF*v.sqr()+k));   
}

/*****************************************************************************************
 * LES3DFsd_pState::H -- Return total enthalpy of the mixture.                           *
 *****************************************************************************************/
double LES3DFsd_pState::H(void) const{
   return (rho*(h()+HALF*v.sqr()+FIVE*k/THREE));   
}

/*****************************************************************************************
 * LES3DFsd_pState::Hs -- Return total mixture sensible enthalpy.                        *
 *****************************************************************************************/
double LES3DFsd_pState::Hs(void) const{
   return (rho*(hs()+HALF*v.sqr()+FIVE*k/THREE));   
}

/*****************************************************************************************
 * LES3DFsd_pState::p_t -- Return turbulence modified pressure.                           *
 *****************************************************************************************/
double LES3DFsd_pState::p_t(void) const { 
   return (p + 2.0*rho*k/3.0);
}

/********************************************************************************************
 * LES3DFsd_pState::a_t -- Return mixture sound speed (including turbulent kinetic energy). *
 ********************************************************************************************/
double LES3DFsd_pState::a_t(void) {
   double aa = sqr(a());
   aa += (TWO/THREE)*k*g();
   return sqrt(aa);
}

double LES3DFsd_pState::a_t(void) const {
   double aa = sqr(a());
   aa += (TWO/THREE)*k*g();
   return sqrt(aa);
}

/**********************************************************************************************
 * LES3DFsd_pState::premixed_mfrac -- Species mass fractions for one-step reaction mechanism  *
 **********************************************************************************************/
void LES3DFsd_pState::premixed_mfrac(void) {
  double unburnt_oxygen_c, burnt_fuel_c, burnt_oxygen_c, stoich_ratio;
  double c_products, products_ratio;
  
  // Check realizability of progress variable.
  Realizable_C_Check();  

  // Determine stoichiometric ratio
  if (React.reactset_flag == CH4_1STEP) {
    stoich_ratio = 2.0*specdata[1].Mol_mass()/
                   specdata[0].Mol_mass();     // stoichiometric O2/CH4 mass ratio
  } else if (React.reactset_flag == C3H8_1STEP) {
    stoich_ratio = 5.0*specdata[1].Mol_mass()/
                   specdata[0].Mol_mass();     // stoichiometric O2/C3H8 mass ratio
  } else if (React.reactset_flag == H2O2_1STEP){
    stoich_ratio = specdata[1].Mol_mass()/
                   2.0/specdata[0].Mol_mass(); // stoichiometric O2/H2 mass ratio
  } /* endif */

  // Calculate fuel mass fraction
  spec[0].c = (ONE - C)*_unburnt_fuel_mass_fraction;
  if (spec[0].c < MICRO) {
      spec[0].c = ZERO;
  } /* endif */

  // Calculate oxygen mass fraction
  // Lean mixture(_fuel_equivalence_ratio < 1) => excessive oxygen
  if (_fuel_equivalence_ratio < ONE) {  
    burnt_fuel_c = ZERO;
    burnt_oxygen_c = (ONE/_fuel_equivalence_ratio - ONE)*_unburnt_fuel_mass_fraction*stoich_ratio;
    if (spec[0].c == ZERO) {
      spec[1].c = burnt_oxygen_c;      
    } else {
      spec[1].c = spec[0].c*stoich_ratio + burnt_oxygen_c;
    } /* endif */
    unburnt_oxygen_c = _unburnt_fuel_mass_fraction*stoich_ratio + burnt_oxygen_c;

  // Rich mixture(_fuel_equivalence_ratio > 1) => excessive fuel
  } else if (_fuel_equivalence_ratio > ONE) {  
    burnt_oxygen_c = ZERO;
    burnt_fuel_c = (ONE - ONE/_fuel_equivalence_ratio)*_unburnt_fuel_mass_fraction;
    spec[0].c = C*(burnt_fuel_c - _unburnt_fuel_mass_fraction) + _unburnt_fuel_mass_fraction;
    if (spec[0].c <= burnt_fuel_c) {
      spec[1].c = burnt_oxygen_c;
    } else {
      spec[1].c = (spec[0].c - burnt_fuel_c)*stoich_ratio;
    }
    unburnt_oxygen_c = (_unburnt_fuel_mass_fraction - burnt_fuel_c)*stoich_ratio;

  // Stoichiometric mixture(_fuel_equivalence_ratio = 1)
  } else { 
    burnt_fuel_c = ZERO;
    burnt_oxygen_c = ZERO;
    spec[1].c = spec[0].c*stoich_ratio; 
    unburnt_oxygen_c = _unburnt_fuel_mass_fraction*stoich_ratio; 
  } /* endif */

  // Determine mass fraction of nitrogen, N2
  if (React.reactset_flag == CH4_1STEP ||
      React.reactset_flag == C3H8_1STEP) {
     spec[4].c = ONE - _unburnt_fuel_mass_fraction - unburnt_oxygen_c;
  } else if (React.reactset_flag == H2O2_1STEP) {
     spec[3].c = ONE - _unburnt_fuel_mass_fraction - unburnt_oxygen_c;
  } /* endif */

  // Determine mass fractions of products
  if (React.reactset_flag == CH4_1STEP ||
      React.reactset_flag == C3H8_1STEP) {
    // mass fractions of products
    if (C <= MICRO) {
       c_products = ZERO;
    } else {
       c_products = ONE - (spec[0].c + spec[1].c + spec[4].c);
    } /* endif */
    if (React.reactset_flag == CH4_1STEP) {
       products_ratio = specdata[2].Mol_mass()/
                        (specdata[2].Mol_mass()+2.0*specdata[3].Mol_mass());
    } else if (React.reactset_flag == C3H8_1STEP) {
       products_ratio = 3.0*specdata[2].Mol_mass()/
                        (3.0*specdata[2].Mol_mass()+4.0*specdata[3].Mol_mass());
    } /* endif */ 
    spec[2].c = products_ratio*c_products; // CO2 mass fraction
    spec[3].c = c_products-spec[2].c;  // H2O mass fraction   
    spec[1].c = (ONE - c_products - spec[0].c - spec[4].c) ;  // O2 mass fraction
  } else if (React.reactset_flag == H2O2_1STEP) {
       spec[2].c = ONE - spec[0].c - spec[1].c - spec[3].c; // H2O mass fraction 
  } /* endif */

  // Check for realizability of each species mass fraction
  for (int i=0; i<ns; i++){
    if (spec[i].c < ZERO) spec[i].c = ZERO;
  } /* endfor */

  double suma = ZERO;
  for (int i=0; i < ns; i++) {
    suma = suma + spec[i].c;
  } /* endfor */

  if (suma > ZERO) {
     for (int i=0; i < ns; i++) {
       spec[i].c = spec[i].c*(ONE/suma);
     } /* endfor */
  } /* endif */

}

/*******************************************************************
 * LES3DFsd_pState::mu_t -- Return eddy (turbulent) viscosity.     *
 *******************************************************************/
double LES3DFsd_pState::mu_t(const LES3DFsd_pState &dWdx,
 			     const LES3DFsd_pState &dWdy,
			     const LES3DFsd_pState &dWdz,
                             const int Flow_Type, 
                             const double &Volume) {
  double filter = filter_width(Volume);
  if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY) {
    double Cs = 0.18;
    return(rho*Cs*sqr(filter)*abs_strain_rate(dWdx,dWdy,dWdz));
  } else if(Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K) {
    double Cv = 0.086;
    return(rho*Cv*sqrt(k)*filter);
  }
}

/**************************************************************************
 * LES3DFsd_pState::kappa_t -- Return turbulent thermal conductivity.     *
 **************************************************************************/
double LES3DFsd_pState::kappa_t(const LES3DFsd_pState &dWdx,
 			        const LES3DFsd_pState &dWdy,
			        const LES3DFsd_pState &dWdz,
                                const int Flow_Type, 
                                const double &Volume) {
  return (mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume)*Cp()/Pr_t()); 
}

/********************************************************************************
 * LES3DFsd_pState::Ds_t -- Return species turbulent diffusion coefficient.     *
 ********************************************************************************/
double LES3DFsd_pState::Ds_t(const int i,
                             const LES3DFsd_pState &dWdx,
 			     const LES3DFsd_pState &dWdy,
			     const LES3DFsd_pState &dWdz,
                             const int Flow_Type, 
                             const double &Volume) {
  return (mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume)/(rho*Sc_t())); 
}

double LES3DFsd_pState::Ds_t(const int i,
                             const double &mu_t_temp) {
  return (mu_t_temp/(rho*Sc_t())); 
}

/*******************************************************************
 * LES3DFsd_pState::Pr_t -- Return turbulent Prandtl number.       *
 *******************************************************************/
double LES3DFsd_pState::Pr_t(void) {
   return (0.9); 
}

/*******************************************************************
 * LES3DFsd_pState::Sc_t -- Return turbulent Schmidt number.       *
 *******************************************************************/
double LES3DFsd_pState::Sc_t(void) {
   return (1.0);
}

/*******************************************************************
 * LES3DFsd_pState::Le_t -- Return turbulent Lewis number.         *
 *******************************************************************/
double LES3DFsd_pState::Le_t(void) {
  return (Sc_t()/Pr_t());
}

/*************************************************************************
 * LES3DFsd_pState::filter_width -- LES characteristic filter width      *
 *************************************************************************/
double LES3DFsd_pState::filter_width(const double &Volume) const {
  return (0.0366/60.0);//pow(Volume,1.0/3.0); 
}

/*******************************************************************************
 * LES3DFsd_pState::tau_t -- Return subfilter scale (turbulent) stress tensor. *
 *******************************************************************************/
Tensor3D LES3DFsd_pState::tau_t(const LES3DFsd_pState &dWdx, 
                                const LES3DFsd_pState &dWdy,
                                const LES3DFsd_pState &dWdz,
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

Tensor3D LES3DFsd_pState::tau_t(const double &mu_t_temp,
                                const LES3DFsd_pState &dWdx, 
                                const LES3DFsd_pState &dWdy,
                                const LES3DFsd_pState &dWdz,
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
 * LES3DFsd_pState::tau_t_x -- Return components of subfilter scale (turbulent)               *
 *                             stress tensor in the x-direction.                              *
 **********************************************************************************************/
Tensor3D LES3DFsd_pState::tau_t_x(const LES3DFsd_pState &dWdx, 
                                  const LES3DFsd_pState &dWdy,
                                  const LES3DFsd_pState &dWdz,
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

Tensor3D LES3DFsd_pState::tau_t_x(const double &mu_t_temp,
                                  const LES3DFsd_pState &dWdx, 
                                  const LES3DFsd_pState &dWdy,
                                  const LES3DFsd_pState &dWdz,
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
 * LES3DFsd_pState::tau_t_y -- Return components of subfilter scale (turbulent)               *
 *                             stress tensor in the y-direction.                              *
 **********************************************************************************************/
Tensor3D LES3DFsd_pState::tau_t_y(const LES3DFsd_pState &dWdx, 
                                  const LES3DFsd_pState &dWdy,
                                  const LES3DFsd_pState &dWdz,
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

Tensor3D LES3DFsd_pState::tau_t_y(const double &mu_t_temp,
                                  const LES3DFsd_pState &dWdx, 
                                  const LES3DFsd_pState &dWdy,
                                  const LES3DFsd_pState &dWdz,
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
 * LES3DFsd_pState::tau_t_z -- Return components of subfilter scale (turbulent)               *
 *                             stress tensor in the z-direction.                              *
 **********************************************************************************************/
Tensor3D LES3DFsd_pState::tau_t_z(const LES3DFsd_pState &dWdx, 
                                  const LES3DFsd_pState &dWdy,
                                  const LES3DFsd_pState &dWdz,
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

Tensor3D LES3DFsd_pState::tau_t_z(const double &mu_t_temp,
                                  const LES3DFsd_pState &dWdx, 
                                  const LES3DFsd_pState &dWdy,
                                  const LES3DFsd_pState &dWdz,
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
 * LES3DFsd_pState::q_t -- Return turbulent heat flux vector.           *
 ************************************************************************/
Vector3D LES3DFsd_pState::q_t(const LES3DFsd_pState &dWdx, 
                              const LES3DFsd_pState &dWdy,
                              const LES3DFsd_pState &dWdz,
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

Vector3D LES3DFsd_pState::q_t(const double &kappa_t_temp,
                              const LES3DFsd_pState &dWdx, 
                              const LES3DFsd_pState &dWdy,
                              const LES3DFsd_pState &dWdz,
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
 * LES3DFsd_pState::q_t_x -- Return component of turbulent heat flux    *
 *                           vector in the x-direction.                 *
 ************************************************************************/
Vector3D LES3DFsd_pState::q_t_x(const LES3DFsd_pState &dWdx, 
                                const LES3DFsd_pState &dWdy,
                                const LES3DFsd_pState &dWdz,
                                const int Flow_Type, 
                                const double &Volume) {
  
    double Rmix = Rtot();
    double kappa_t_temp = kappa_t(dWdx,dWdy,dWdz,Flow_Type,Volume);
    Vector3D heat_flux;
     
    heat_flux.x = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdx.p -(p/rho)*dWdx.rho);
    return (heat_flux);

}

Vector3D LES3DFsd_pState::q_t_x(const double &kappa_t_temp,
                                const LES3DFsd_pState &dWdx, 
                                const LES3DFsd_pState &dWdy,
                                const LES3DFsd_pState &dWdz,
                                const int Flow_Type, 
                                const double &Volume) {
  
    double Rmix = Rtot();
    Vector3D heat_flux;
     
    heat_flux.x = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdx.p -(p/rho)*dWdx.rho);
    return (heat_flux);

}

/************************************************************************
 * LES3DFsd_pState::q_t_y -- Return component of turbulent heat flux    *
 *                           vector in the y-direction.                 *
 ************************************************************************/
Vector3D LES3DFsd_pState::q_t_y(const LES3DFsd_pState &dWdx, 
                                const LES3DFsd_pState &dWdy,
                                const LES3DFsd_pState &dWdz,
                                const int Flow_Type, 
                                const double &Volume) {
  
    double Rmix = Rtot();
    double kappa_t_temp = kappa_t(dWdx,dWdy,dWdz,Flow_Type,Volume);
    Vector3D heat_flux;
     
    heat_flux.y = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdy.p -(p/rho)*dWdy.rho);
    return (heat_flux);

}

Vector3D LES3DFsd_pState::q_t_y(const double &kappa_t_temp,
                                const LES3DFsd_pState &dWdx, 
                                const LES3DFsd_pState &dWdy,
                                const LES3DFsd_pState &dWdz,
                                const int Flow_Type, 
                                const double &Volume) {
  
    double Rmix = Rtot();
    Vector3D heat_flux;
     
    heat_flux.y = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdy.p -(p/rho)*dWdy.rho);
    return (heat_flux);

}

/************************************************************************
 * LES3DFsd_pState::q_t_z -- Return component of turbulent heat flux    *
 *                           vector in the z-direction.                 *
 ************************************************************************/
Vector3D LES3DFsd_pState::q_t_z(const LES3DFsd_pState &dWdx, 
                                const LES3DFsd_pState &dWdy,
                                const LES3DFsd_pState &dWdz,
                                const int Flow_Type, 
                                const double &Volume) {
  
    double Rmix = Rtot();
    double kappa_t_temp = kappa_t(dWdx,dWdy,dWdz,Flow_Type,Volume);
    Vector3D heat_flux;
     
    heat_flux.z = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdz.p -(p/rho)*dWdz.rho);
    return (heat_flux);

}

Vector3D LES3DFsd_pState::q_t_z(const double &kappa_t_temp,
                                const LES3DFsd_pState &dWdx, 
                                const LES3DFsd_pState &dWdy,
                                const LES3DFsd_pState &dWdz,
                                const int Flow_Type, 
                                const double &Volume) {
  
    double Rmix = Rtot();
    Vector3D heat_flux;
     
    heat_flux.z = -kappa_t_temp*(1.0/(rho*Rmix)) * (dWdz.p -(p/rho)*dWdz.rho);
    return (heat_flux);

}

/*******************************************************************
 * LES3DFsd_pState::U -- Return conserved solution state vector.   *
 *******************************************************************/
LES3DFsd_cState LES3DFsd_pState::U(void) {
   LES3DFsd_cState Temp;
   Temp.rho = rho;
   Temp.rhov = rhov();
   Temp.E = E();
   Temp.rhoC = rho*C;
   Temp.rhoFsd = rho*Fsd;
   Temp.rhok = rho*k;
   return Temp;
}

LES3DFsd_cState LES3DFsd_pState::U(void) const {
   LES3DFsd_cState Temp;
   Temp.rho = rho;
   Temp.rhov = rhov();
   Temp.E = E();
   Temp.rhoC = rho*C;
   Temp.rhoFsd = rho*Fsd;
   Temp.rhok = rho*k;
   return Temp;
}

LES3DFsd_cState LES3DFsd_pState::U(const LES3DFsd_pState &W) {
    LES3DFsd_cState Temp;
    Temp.rho = W.rho;
    Temp.rhov = W.rhov();
    Temp.E = W.E();
    Temp.rhoC = W.rho*W.C;
    Temp.rhoFsd = W.rho*W.Fsd;
    Temp.rhok = W.rho*W.k;
    return Temp;
}

/***********************************************************************
 ***************** INVISCID FLUX VECTORS *******************************
 ***********************************************************************/

/***********************************************************
 *    LES3DFsd_pState::F -- Inviscid flux (x-direction).   *
 ***********************************************************/
LES3DFsd_cState LES3DFsd_pState::F(void) {
   LES3DFsd_cState Temp;
   Temp.rho = rho*v.x;
   Temp.rhov.x = rho*sqr(v.x) + p + (TWO/THREE)*rho*k;
   Temp.rhov.y = rho*v.x*v.y;
   Temp.rhov.z = rho*v.x*v.z;
   Temp.E = v.x*(H() + (TWO/THREE)*rho*k);
   Temp.rhoC = rho*v.x*C;
   Temp.rhoFsd = rho*v.x*Fsd;
   Temp.rhok = rho*v.x*k;
   return (Temp);
}

LES3DFsd_cState LES3DFsd_pState::F(void) const {
   LES3DFsd_cState Temp;
   Temp.rho = rho*v.x;
   Temp.rhov.x = rho*sqr(v.x) + p + (TWO/THREE)*rho*k;
   Temp.rhov.y = rho*v.x*v.y;
   Temp.rhov.z = rho*v.x*v.z;
   Temp.E = v.x*(H() + (TWO/THREE)*rho*k);
   Temp.rhoC = rho*v.x*C;
   Temp.rhoFsd = rho*v.x*Fsd;
   Temp.rhok = rho*v.x*k;
   return (Temp);
}

/************************************************************
 *    LES3DFsd_pState::Fx -- Inviscid flux (x-direction).   *
 ************************************************************/
LES3DFsd_cState LES3DFsd_pState::Fx(void) {
   LES3DFsd_cState Temp;
   Temp.rho = rho*v.x;
   Temp.rhov.x = rho*sqr(v.x) + p + (TWO/THREE)*rho*k;
   Temp.rhov.y = rho*v.x*v.y;
   Temp.rhov.z = rho*v.x*v.z;
   Temp.E = v.x*(H() + (TWO/THREE)*rho*k);
   Temp.rhoC = rho*v.x*C;
   Temp.rhoFsd = rho*v.x*Fsd;
   Temp.rhok = rho*v.x*k;
   return (Temp);
}

LES3DFsd_cState LES3DFsd_pState::Fx(void) const {
   LES3DFsd_cState Temp;
   Temp.rho = rho*v.x;
   Temp.rhov.x = rho*sqr(v.x) + p + (TWO/THREE)*rho*k;
   Temp.rhov.y = rho*v.x*v.y;
   Temp.rhov.z = rho*v.x*v.z;
   Temp.E = v.x*(H() + (TWO/THREE)*rho*k);
   Temp.rhoC = rho*v.x*C;
   Temp.rhoFsd = rho*v.x*Fsd;
   Temp.rhok = rho*v.x*k;
   return (Temp);
}

/************************************************************
 *    LES3DFsd_pState::Fy -- Inviscid flux (y-direction).   *
 ************************************************************/
LES3DFsd_cState LES3DFsd_pState::Fy(void) {
   LES3DFsd_cState Temp;
   Temp.rho = rho*v.y;
   Temp.rhov.x = rho*v.x*v.y;
   Temp.rhov.y = rho*sqr(v.y) + p + (TWO/THREE)*rho*k;
   Temp.rhov.z = rho*v.y*v.z;
   Temp.E = v.y*(H() + (TWO/THREE)*rho*k);
   Temp.rhoC = rho*v.y*C;
   Temp.rhoFsd = rho*v.y*Fsd;
   Temp.rhok = rho*v.y*k;
   return (Temp);
}

LES3DFsd_cState LES3DFsd_pState::Fy(void) const {
   LES3DFsd_cState Temp;
   Temp.rho = rho*v.y;
   Temp.rhov.x = rho*v.x*v.y;
   Temp.rhov.y = rho*sqr(v.y) + p + (TWO/THREE)*rho*k;
   Temp.rhov.z = rho*v.y*v.z;
   Temp.E = v.y*(H() + (TWO/THREE)*rho*k);
   Temp.rhoC = rho*v.y*C;
   Temp.rhoFsd = rho*v.y*Fsd;
   Temp.rhok = rho*v.y*k;
   return (Temp);
}

/************************************************************
 *    LES3DFsd_pState::Fz -- Inviscid flux (z-direction).   *
 ************************************************************/
LES3DFsd_cState LES3DFsd_pState::Fz(void) {
   LES3DFsd_cState Temp;
   Temp.rho = rho*v.z;
   Temp.rhov.x = rho*v.x*v.z;
   Temp.rhov.y = rho*v.y*v.z;
   Temp.rhov.z = rho*sqr(v.z) + p + (TWO/THREE)*rho*k;
   Temp.E = v.z*(H() + (TWO/THREE)*rho*k);
   Temp.rhoC = rho*v.z*C;
   Temp.rhoFsd = rho*v.z*Fsd;
   Temp.rhok = rho*v.z*k;
   return (Temp);
}

LES3DFsd_cState LES3DFsd_pState::Fz(void) const {
   LES3DFsd_cState Temp;
   Temp.rho = rho*v.z;
   Temp.rhov.x = rho*v.x*v.z;
   Temp.rhov.y = rho*v.y*v.z;
   Temp.rhov.z = rho*sqr(v.z) + p + (TWO/THREE)*rho*k;
   Temp.E = v.z*(H() + (TWO/THREE)*rho*k);
   Temp.rhoC = rho*v.z*C;
   Temp.rhoFsd = rho*v.z*Fsd;
   Temp.rhok = rho*v.z*k;
   return (Temp);
}

/**********************************************************************
 ***************** VISCOUS FLUX VECTORS *******************************
 **********************************************************************/

/********************************************************
 *  LES3DFsd_pState::Fv -- Viscous flux (x-direction).  * 
 ********************************************************/
LES3DFsd_cState LES3DFsd_pState::Fv(const LES3DFsd_pState &dWdx,
                                    const LES3DFsd_pState &dWdy,
                                    const LES3DFsd_pState &dWdz,
                                    const int Flow_Type,
                                    const double &Volume) {

   LES3DFsd_cState Temp;

   double mu_temp = mu();
   double mu_t_temp = mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume);
   double kappa_temp = kappa();
   double kappa_t_temp = kappa_t(dWdx,dWdy,dWdz,Flow_Type,Volume);

   Tensor3D fluid_stress;
   Vector3D heat_flux;
   
   fluid_stress = tau_t_x(mu_temp + mu_t_temp, dWdx, dWdy, dWdz, Flow_Type, Volume);
   heat_flux = q_t_x(kappa_temp + kappa_t_temp, dWdx, dWdy, dWdz, Flow_Type, Volume);
   
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
            heat_flux.x;
   Temp.rhoC = mu_t_temp*dWdx.C/Sc_t();
   Temp.rhoFsd = mu_t_temp*dWdx.Fsd/Sc_t();
   Temp.rhok = (mu_temp+mu_t_temp/0.25)*dWdx.k;

   return (Temp);

}

/*********************************************************
 *  LES3DFsd_pState::Fvx -- Viscous flux (x-direction).  * 
 *********************************************************/
LES3DFsd_cState LES3DFsd_pState::Fvx(const LES3DFsd_pState &dWdx,
                                     const LES3DFsd_pState &dWdy,
                                     const LES3DFsd_pState &dWdz,
                                     const int Flow_Type,
                                     const double &Volume) {

   LES3DFsd_cState Temp;

   double mu_temp = mu();
   double mu_t_temp = mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume);
   double kappa_temp = kappa();
   double kappa_t_temp = kappa_t(dWdx,dWdy,dWdz,Flow_Type,Volume);

   Tensor3D fluid_stress;
   Vector3D heat_flux;
   
   fluid_stress = tau_t_x(mu_temp + mu_t_temp, dWdx, dWdy, dWdz, Flow_Type, Volume);
   heat_flux = q_t_x(kappa_temp + kappa_t_temp, dWdx, dWdy, dWdz, Flow_Type, Volume);
   
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
            heat_flux.x;
   Temp.rhoC = mu_t_temp*dWdx.C/Sc_t();
   Temp.rhoFsd = mu_t_temp*dWdx.Fsd/Sc_t();
   Temp.rhok = (mu_temp+mu_t_temp/0.25)*dWdx.k;

   return (Temp);

}

/*********************************************************
 *  LES3DFsd_pState::Fvy -- Viscous flux (y-direction).  * 
 *********************************************************/
LES3DFsd_cState LES3DFsd_pState::Fvy(const LES3DFsd_pState &dWdx,
                                     const LES3DFsd_pState &dWdy,
                                     const LES3DFsd_pState &dWdz,
                                     const int Flow_Type,
                                     const double &Volume) {

   LES3DFsd_cState Temp;

   double mu_temp = mu();
   double mu_t_temp = mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume);
   double kappa_temp = kappa();
   double kappa_t_temp = kappa_t(dWdx,dWdy,dWdz,Flow_Type,Volume);

   Tensor3D fluid_stress;
   Vector3D heat_flux;
   
   fluid_stress = tau_t_y(mu_temp + mu_t_temp, dWdx, dWdy, dWdz, Flow_Type, Volume);
   heat_flux = q_t_y(kappa_temp + kappa_t_temp, dWdx, dWdy, dWdz, Flow_Type, Volume);
   
   for (int index = 0; index < ns; ++index) {
      _diff_coeff[index] = Ds(index, mu_temp)+Ds_t(index, mu_t_temp);
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
   Temp.rhoC = mu_t_temp*dWdy.C/Sc_t();
   Temp.rhoFsd = mu_t_temp*dWdy.Fsd/Sc_t();
   Temp.rhok = (mu_temp+mu_t_temp/0.25)*dWdy.k;

   return (Temp);

}

/*********************************************************
 *  LES3DFsd_pState::Fvz -- Viscous flux (z-direction).  * 
 *********************************************************/
LES3DFsd_cState LES3DFsd_pState::Fvz(const LES3DFsd_pState &dWdx,
                                     const LES3DFsd_pState &dWdy,
                                     const LES3DFsd_pState &dWdz,
                                     const int Flow_Type,
                                     const double &Volume) {

   LES3DFsd_cState Temp;

   double mu_temp = mu();
   double mu_t_temp = mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume);
   double kappa_temp = kappa();
   double kappa_t_temp = kappa_t(dWdx,dWdy,dWdz,Flow_Type,Volume);

   Tensor3D fluid_stress;
   Vector3D heat_flux;
   
   fluid_stress = tau_t_z(mu_temp + mu_t_temp, dWdx, dWdy, dWdz, Flow_Type, Volume);
   heat_flux = q_t_z(kappa_temp + kappa_t_temp, dWdx, dWdy, dWdz, Flow_Type, Volume);
   
   for (int index = 0; index < ns; ++index) {
      _diff_coeff[index] = Ds(index, mu_temp)+Ds_t(index, mu_t_temp);
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
   Temp.rhoC = mu_t_temp*dWdz.C/Sc_t();
   Temp.rhoFsd = mu_t_temp*dWdz.Fsd/Sc_t();
   Temp.rhok = (mu_temp+mu_t_temp/0.25)*dWdz.k;

   return (Temp);

}

/************************************************************************
 ********************** EIGENVALUES *************************************
 ************************************************************************/

/**************************************************************
 * LES3DFsd_pState::lambda -- Eigenvalue(s) (x-direction).    *
 **************************************************************/
LES3DFsd_pState LES3DFsd_pState::lambda(void) {
  double cc = a_t();
  LES3DFsd_pState Temp;
  Temp.rho = v.x - cc;
  Temp.v.x = v.x;
  Temp.v.y = v.x;
  Temp.v.z = v.x;
  Temp.p = v.x + cc;
  Temp.C = v.x;
  Temp.Fsd = v.x;
  Temp.k = v.x;
  return (Temp);
}

LES3DFsd_pState LES3DFsd_pState::lambda(void) const {
  double cc = a_t();
  LES3DFsd_pState Temp;
  Temp.rho = v.x - cc;
  Temp.v.x = v.x;
  Temp.v.y = v.x;
  Temp.v.z = v.x;
  Temp.p = v.x + cc;
  Temp.C = v.x;
  Temp.Fsd = v.x;
  Temp.k = v.x;
  return (Temp);
}

/**************************************************************
 * LES3DFsd_pState::lambda_x -- Eigenvalue(s) (x-direction).  *
 **************************************************************/
LES3DFsd_pState LES3DFsd_pState::lambda_x(void) {
  double cc = a_t();
  LES3DFsd_pState Temp;
  Temp.rho = v.x - cc;
  Temp.v.x = v.x;
  Temp.v.y = v.x;
  Temp.v.z = v.x;
  Temp.p = v.x + cc;
  Temp.C = v.x;
  Temp.Fsd = v.x;
  Temp.k = v.x;
  return (Temp);
}

LES3DFsd_pState LES3DFsd_pState::lambda_x(void) const {
  double cc = a_t();
  LES3DFsd_pState Temp;
  Temp.rho = v.x - cc;
  Temp.v.x = v.x;
  Temp.v.y = v.x;
  Temp.v.z = v.x;
  Temp.p = v.x + cc;
  Temp.C = v.x;
  Temp.Fsd = v.x;
  Temp.k = v.x;
  return (Temp);
}

/************************************************************************
 ********************* EIGENVECTORS *************************************
 ************************************************************************/

/***********************************************************************
 * LES3DFsd_pState::rc -- Conserved right eigenvector (x-direction).   *
 ***********************************************************************/
LES3DFsd_cState LES3DFsd_pState::rc(const int &index) {
   double cc, eta_fsd;  
   switch(index){  
   case 1:
     cc = a_t();
     return (LES3DFsd_cState(ONE, v.x-cc, v.y, v.z, H()/rho-v.x*cc, C, Fsd, k));
   case 2:
     cc = a_t();
     return (LES3DFsd_cState(ONE, v.x, v.y, v.z, H()/rho-cc*cc/(g()-1.0), C, Fsd, k)); 
   case 3:
     return (LES3DFsd_cState(ZERO, ZERO, rho, ZERO, rho*v.y, ZERO, ZERO, ZERO));
   case 4:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, rho, rho*v.z, ZERO, ZERO, ZERO));
   case 5:
     cc = a_t();
     return (LES3DFsd_cState(ONE, v.x+cc, v.y, v.z, H()/rho+v.x*cc, C, Fsd, k));
   case 6:
     eta_fsd = Progvar_Species_Grad();
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, rho*eta_fsd, rho, ZERO, ZERO));
   case 7:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, rho, ZERO));
   case 8:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, 5.0*rho/3.0, ZERO, ZERO, rho));
    };
}

LES3DFsd_cState LES3DFsd_pState::rc(const int &index) const {
   double cc, eta_fsd;  
   switch(index){  
   case 1:
     cc = a_t();
     return (LES3DFsd_cState(ONE, v.x-cc, v.y, v.z, H()/rho-v.x*cc, C, Fsd, k));
   case 2:
     cc = a_t();
     return (LES3DFsd_cState(ONE, v.x, v.y, v.z, H()/rho-cc*cc/(g()-1.0), C, Fsd, k)); 
   case 3:
     return (LES3DFsd_cState(ZERO, ZERO, rho, ZERO, rho*v.y, ZERO, ZERO, ZERO));
   case 4:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, rho, rho*v.z, ZERO, ZERO, ZERO));
   case 5:
     cc = a_t();
     return (LES3DFsd_cState(ONE, v.x+cc, v.y, v.z, H()/rho+v.x*cc, C, Fsd, k));
   case 6:
     eta_fsd = Progvar_Species_Grad();
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, rho*eta_fsd, rho, ZERO, ZERO));
   case 7:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, rho, ZERO));
   case 8:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, 5.0*rho/3.0, ZERO, ZERO, rho));
    };
}

/*************************************************************************
 * LES3DFsd_pState::rc_x -- Conserved right eigenvector (x-direction).   *
 *************************************************************************/
LES3DFsd_cState LES3DFsd_pState::rc_x(const int &index) const {
   double cc, eta_fsd;  
   switch(index){  
   case 1:
     cc = a_t();
     return (LES3DFsd_cState(ONE, v.x-cc, v.y, v.z, H()/rho-v.x*cc, C, Fsd, k));
   case 2:
     cc = a_t();
     return (LES3DFsd_cState(ONE, v.x, v.y, v.z, H()/rho-cc*cc/(g()-1.0), C, Fsd, k)); 
   case 3:
     return (LES3DFsd_cState(ZERO, ZERO, rho, ZERO, rho*v.y, ZERO, ZERO, ZERO));
   case 4:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, rho, rho*v.z, ZERO, ZERO, ZERO));
   case 5:
     cc = a_t();
     return (LES3DFsd_cState(ONE, v.x+cc, v.y, v.z, H()/rho+v.x*cc, C, Fsd, k));
   case 6:
     eta_fsd = Progvar_Species_Grad();
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, rho*eta_fsd, rho, ZERO, ZERO));
   case 7:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, rho, ZERO));
   case 8:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, 5.0*rho/3.0, ZERO, ZERO, rho));
    };
}

LES3DFsd_cState LES3DFsd_pState::rc_x(const int &index) {
   double cc, eta_fsd;  
   switch(index){  
   case 1:
     cc = a_t();
     return (LES3DFsd_cState(ONE, v.x-cc, v.y, v.z, H()/rho-v.x*cc, C, Fsd, k));
   case 2:
     cc = a_t();
     return (LES3DFsd_cState(ONE, v.x, v.y, v.z, H()/rho-cc*cc/(g()-1.0), C, Fsd, k)); 
   case 3:
     return (LES3DFsd_cState(ZERO, ZERO, rho, ZERO, rho*v.y, ZERO, ZERO, ZERO));
   case 4:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, rho, rho*v.z, ZERO, ZERO, ZERO));
   case 5:
     cc = a_t();
     return (LES3DFsd_cState(ONE, v.x+cc, v.y, v.z, H()/rho+v.x*cc, C, Fsd, k));
   case 6:
     eta_fsd = Progvar_Species_Grad();
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, rho*eta_fsd, rho, ZERO, ZERO));
   case 7:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, rho, ZERO));
   case 8:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, 5.0*rho/3.0, ZERO, ZERO, rho));
    };
}

/***********************************************************************
 * LES3DFsd_pState::lp -- Primitive left eigenvector (x-direction).    *
 ***********************************************************************/
LES3DFsd_pState LES3DFsd_pState::lp(const int &index) {
   double cc;  
   switch(index){  
   case 1:
      cc = a_t();
      return (LES3DFsd_pState(ZERO, -HALF*rho/cc, ZERO, ZERO, HALF/(cc*cc), ZERO, ZERO, ZERO));
   case 2:
      cc = a_t();      
      return (LES3DFsd_pState(ONE, ZERO, ZERO, ZERO, -ONE/(cc*cc), ZERO, ZERO, ZERO));
   case 3 :
      return (LES3DFsd_pState(ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO, ZERO));
   case 4 :
      return (LES3DFsd_pState(ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO));
   case 5 :
      cc = a_t();
      return (LES3DFsd_pState(ZERO, HALF*rho/cc, ZERO, ZERO, HALF/(cc*cc), ZERO, ZERO, ZERO));
   case 6 :
     return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO, ZERO));
   case 7 :
     return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO));
   case 8 :
     return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE));
   };
}

LES3DFsd_pState LES3DFsd_pState::lp(const int &index) const {
   double cc;  
   switch(index){  
   case 1:
      cc = a_t();
      return (LES3DFsd_pState(ZERO, -HALF*rho/cc, ZERO, ZERO, HALF/(cc*cc), ZERO, ZERO, ZERO));
   case 2:
      cc = a_t();      
      return (LES3DFsd_pState(ONE, ZERO, ZERO, ZERO, -ONE/(cc*cc), ZERO, ZERO, ZERO));
   case 3 :
      return (LES3DFsd_pState(ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO, ZERO));
   case 4 :
      return (LES3DFsd_pState(ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO));
   case 5 :
      cc = a_t();
      return (LES3DFsd_pState(ZERO, HALF*rho/cc, ZERO, ZERO, HALF/(cc*cc), ZERO, ZERO, ZERO));
   case 6 :
     return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO, ZERO));
   case 7 :
     return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO));
   case 8 :
     return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE));
   };
}

/*************************************************************************
 * LES3DFsd_pState::lp_x -- Primitive left eigenvector (x-direction).    *
 *************************************************************************/
LES3DFsd_pState LES3DFsd_pState::lp_x(const int &index) {
   double cc;  
   switch(index){  
   case 1:
      cc = a_t();
      return (LES3DFsd_pState(ZERO, -HALF*rho/cc, ZERO, ZERO, HALF/(cc*cc), ZERO, ZERO, ZERO));
   case 2:
      cc = a_t();      
      return (LES3DFsd_pState(ONE, ZERO, ZERO, ZERO, -ONE/(cc*cc), ZERO, ZERO, ZERO));
   case 3 :
      return (LES3DFsd_pState(ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO, ZERO));
   case 4 :
      return (LES3DFsd_pState(ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO));
   case 5 :
      cc = a_t();
      return (LES3DFsd_pState(ZERO, HALF*rho/cc, ZERO, ZERO, HALF/(cc*cc), ZERO, ZERO, ZERO));
   case 6 :
     return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO, ZERO));
   case 7 :
     return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO));
   case 8 :
     return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE));
   };
}

LES3DFsd_pState LES3DFsd_pState::lp_x(const int &index) const {
   double cc;  
   switch(index){  
   case 1:
      cc = a_t();
      return (LES3DFsd_pState(ZERO, -HALF*rho/cc, ZERO, ZERO, HALF/(cc*cc), ZERO, ZERO, ZERO));
   case 2:
      cc = a_t();      
      return (LES3DFsd_pState(ONE, ZERO, ZERO, ZERO, -ONE/(cc*cc), ZERO, ZERO, ZERO));
   case 3 :
      return (LES3DFsd_pState(ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO, ZERO));
   case 4 :
      return (LES3DFsd_pState(ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO));
   case 5 :
      cc = a_t();
      return (LES3DFsd_pState(ZERO, HALF*rho/cc, ZERO, ZERO, HALF/(cc*cc), ZERO, ZERO, ZERO));
   case 6 :
     return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO, ZERO));
   case 7 :
     return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO));
   case 8 :
     return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE));
   };
}

/***********************************************************
 * LES3DFsd_pState::Mr2 -- Square of Mach number           *
 ***********************************************************/
double LES3DFsd_pState::Mr2(const double &deltax, 
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
 * LES3DFsd_pState::u_plus_aprecon -- Preconditioned velocity                      *
 ***********************************************************************************/
double LES3DFsd_pState::u_plus_aprecon(const double &u, 
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
 * LES3DFsd_pState::u_a_precon -- Preconditioned velocity and soundspeed           *
 ***********************************************************************************/
void LES3DFsd_pState::u_a_precon(const double &UR2, 
                                 double &uprimed, 
                                 double &cprimed) const {

  double alpha = HALF*(ONE - ONE/(a_t()*a_t())*UR2);
  uprimed = v.x*(ONE - alpha);
  cprimed = sqrt(alpha*alpha*v.x*v.x + UR2); 

}

/******************************************************************************
 * LES3DFsd_pState::rc_x_precon -- Conserved right eigenvector (x-direction)  * 
 *                                 for low Mach number preconditioner.        *
 ******************************************************************************/
LES3DFsd_cState LES3DFsd_pState::rc_x_precon(const int &index, const double &MR2) const {

   double cc, uprimed, cprimed, eta_fsd;

   switch(index){  
   case 1:
     cc = a_t();
     u_a_precon(MR2*cc*cc,uprimed,cprimed);
     return (LES3DFsd_cState(ONE, 
                             (uprimed-cprimed)/MR2, v.y, v.z, 
                             h()+(v.sqr()/MR2)/TWO+(FIVE*k)/THREE-(v.x*cprimed)/MR2, 
                             C, Fsd, k));
   case 2:
     cc = a_t();
     return (LES3DFsd_cState(ONE, v.x, v.y, v.z, (H()/rho-cc*cc/(g()-ONE)), C, Fsd, k)); 
   case 3:
     return (LES3DFsd_cState(ZERO, ZERO, rho, ZERO, rho*v.y, ZERO, ZERO, ZERO));
   case 4:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, rho, rho*v.z, ZERO, ZERO, ZERO));
   case 5:
     cc = a_t();
     u_a_precon(MR2*cc*cc,uprimed,cprimed);
     return (LES3DFsd_cState(ONE, 
                             (uprimed+cprimed)/MR2, v.y, v.z, 
                             h()+(v.sqr()/MR2)/TWO+(FIVE*k)/THREE+(v.x*cprimed)/MR2, 
                             C, Fsd, k));
   case 6:
     eta_fsd = Progvar_Species_Grad();
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, rho*eta_fsd, rho, ZERO, ZERO));
   case 7:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, rho, ZERO));
   case 8:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, FIVE*rho/THREE, ZERO, ZERO, rho));
   };

}

/******************************************************************************
 * LES3DFsd_pState::lp_x_precon -- Primitive left eigenvector (x-direction)   * 
 *                                 for low Mach number preconditioner.        *
 ******************************************************************************/
LES3DFsd_pState LES3DFsd_pState::lp_x_precon(const int &index,const double &MR2) const {

   double cc, uprimed,cprimed;

   switch(index){  
   case 1:
      cc = a_t();
      u_a_precon(MR2*cc*cc,uprimed,cprimed);
      return (LES3DFsd_pState(ZERO, -HALF*rho*MR2/cprimed, ZERO, ZERO, 
                             (-uprimed+cprimed+v.x)/(TWO*cprimed*cc*cc), 
                              ZERO, ZERO, ZERO));
   case 2:
      cc = a_t();
      return (LES3DFsd_pState(ONE, ZERO, ZERO, ZERO, -ONE/(cc*cc), ZERO, ZERO, ZERO));
   case 3 :
      return (LES3DFsd_pState(ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO, ZERO));
   case 4 :
      return (LES3DFsd_pState(ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO));
   case 5 :
      cc = a_t();
      u_a_precon(MR2*cc*cc,uprimed,cprimed);
      return (LES3DFsd_pState(ZERO, HALF*rho*MR2/cprimed, ZERO, ZERO, (
                              uprimed+cprimed-v.x)/(TWO*cprimed*cc*cc), 
                              ZERO, ZERO, ZERO));
   case 6 :
      return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO, ZERO));
   case 7 :
      return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO));
   case 8 :
      return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE));
   };

}

/******************************************************************************
 * LES3DFsd_pState::Low_Mach_Number_Preconditioner -- Preconditoner matrix    * 
 ******************************************************************************/
void LES3DFsd_pState::Low_Mach_Number_Preconditioner(DenseMatrix &P, 
                                                     const double &deltax, 
                                                     const double &lengthx,
                                                     const double &dTime) {  
  double Temp = T();
  double Rmix = Rtot();
  double enthalpy = h();
  double CP = Cp();
  double cc = a_t();
  double pt = p_t();
  double theta = (ONE/(Mr2(deltax,lengthx,dTime)*cc*cc) + (g()-ONE)/(cc*cc));// + ONE/(CP*Temp));
  double eta_fsd = Progvar_Species_Grad();
  double phi = C*eta_fsd;
  double alpha = theta*pt/rho;
  double alpham1 = alpha - ONE;
  double Omega = (Rmix - CP)*pt/(rho*Rmix);
  double beta = enthalpy - CP*pt/(rho*Rmix) - phi;
  double V = HALF*v.sqr();

  P.zero();

  P(0,0) = (alpha*(beta-V)+pt/rho+V+Rmix*Temp-enthalpy+phi)/Omega;
  P(0,1) = v.x*alpham1/Omega;
  P(0,2) = v.y*alpham1/Omega;
  P(0,3) = v.z*alpham1/Omega;
  P(0,4) = -alpham1/Omega;
  P(0,5) = eta_fsd*alpham1/Omega;
  P(0,7) = 5.0*alpham1/(3.0*Omega);
  P(1,0) = v.x*(beta-V)*alpham1/Omega;
  P(1,1) = v.x*v.x*alpham1/Omega+1.0;
  P(1,2) = v.x*v.y*alpham1/Omega;
  P(1,3) = v.x*v.z*alpham1/Omega;
  P(1,4) = -v.x*alpham1/Omega;
  P(1,5) = v.x*eta_fsd*alpham1/Omega;
  P(1,7) = 5.0*v.x*alpham1/(3.0*Omega);
  P(2,0) = v.y*(beta-V)*alpham1/Omega;
  P(2,1) = v.x*v.y*alpham1/Omega;
  P(2,2) = v.y*v.y*alpham1/Omega+1.0;
  P(2,3) = v.y*v.z*alpham1/Omega;
  P(2,4) = -v.y*alpham1/Omega;
  P(2,5) = v.y*eta_fsd*alpham1/Omega;
  P(2,7) = 5.0*v.y*alpham1/(3.0*Omega);
  P(3,0) = v.z*(beta-V)*alpham1/Omega;
  P(3,1) = v.x*v.z*alpham1/Omega;
  P(3,2) = v.y*v.z*alpham1/Omega+1.0;
  P(3,3) = v.z*v.z*alpham1/Omega;
  P(3,4) = -v.z*alpham1/Omega;
  P(3,5) = v.z*eta_fsd*alpham1/Omega;
  P(3,7) = 5.0*v.z*alpham1/(3.0*Omega);
  P(4,0) = (enthalpy+V+5.0*k/3.0)*(beta-V)*alpham1/Omega;
  P(4,1) = v.x*(enthalpy+V+5.0*k/3.0)*alpham1/Omega;
  P(4,2) = v.y*(enthalpy+V+5.0*k/3.0)*alpham1/Omega;
  P(4,3) = v.z*(enthalpy+V+5.0*k/3.0)*alpham1/Omega;
  P(4,4) = -(alpha*(enthalpy+V+5.0*k/3.0)-V-5.0*k/3.0-Rmix*Temp-beta-phi)/Omega;
  P(4,5) = eta_fsd*(enthalpy+V+5.0*k/3.0)*alpham1/Omega;
  P(4,7) = 5.0*(enthalpy+V+5.0*k/3.0)*alpham1/(3.0*Omega);
  P(5,0) = C*(beta-V)*alpham1/Omega;
  P(5,1) = C*v.x*alpham1/Omega;
  P(5,2) = C*v.y*alpham1/Omega;
  P(5,3) = C*v.z*alpham1/Omega;
  P(5,4) = -C*alpham1/Omega;
  P(5,5) = C*eta_fsd*alpham1/Omega+1.0;
  P(5,7) = 5.0*C*alpham1/(3.0*Omega);
  P(6,0) = Fsd*(beta-V)*alpham1/Omega;
  P(6,1) = Fsd*v.x*alpham1/Omega;
  P(6,2) = Fsd*v.y*alpham1/Omega;
  P(6,3) = Fsd*v.z*alpham1/Omega;
  P(6,4) = -Fsd*alpham1/Omega;
  P(6,5) = Fsd*eta_fsd*alpham1/Omega;
  P(6,6) = ONE;
  P(6,7) = 5.0*Fsd*alpham1/(3.0*Omega);
  P(7,0) = k*(beta-V)*alpham1/Omega;
  P(7,1) = k*v.x*alpham1/Omega;
  P(7,2) = k*v.y*alpham1/Omega;
  P(7,3) = k*v.z*alpham1/Omega;
  P(7,4) = -k*alpham1/Omega;
  P(7,5) = k*eta_fsd*alpham1/Omega;
  P(7,7) = 5.0*k*alpham1/(3.0*Omega);

}

/*********************************************************************************************
 * LES3DFsd_pState::Low_Mach_Number_Preconditioner_Inverse -- Preconditoner matrix inverse   * 
 *********************************************************************************************/
void LES3DFsd_pState::Low_Mach_Number_Preconditioner_Inverse(DenseMatrix &Pinv, 
                                                             const double &deltax, 
                                                             const double &lengthx,
                                                             const double &dTime) {

  double Temp = T();
  double Rmix = Rtot();
  double enthalpy = h();
  double CP = Cp();
  double cc = a_t();
  double pt = p_t();
  double theta = (ONE/(Mr2(deltax,lengthx,dTime)*cc*cc) + (g()-ONE)/(cc*cc));// + ONE/(CP*Temp));  

  double eta_fsd = Progvar_Species_Grad();
  double phi = C*eta_fsd;
  double AA = pt*(rho*Rmix-theta*pt*CP);
  double BB = Rmix*rho*(theta*pt-rho);
  double EE = HALF*v.sqr() - enthalpy + phi;
  double CC = EE + CP*pt/(rho*Rmix); 
  double DD = HALF*v.sqr() + enthalpy;
  Pinv.zero();    

  Pinv(0,0) = rho*Rmix*(theta*pt*EE-rho*CC+pt)/AA;
  Pinv(0,1) = -v.x*BB/AA;
  Pinv(0,2) = -v.y*BB/AA;
  Pinv(0,3) = -v.z*BB/AA;
  Pinv(0,4) = BB/AA;
  Pinv(0,5) = -eta_fsd*BB/AA;
  Pinv(0,7) = -5.0*BB/(3.0*AA);
  Pinv(1,0) = v.x*CC*BB/AA;
  Pinv(1,1) = rho*Rmix/AA*(pt+rho*v.x*v.x-theta*pt*(v.x*v.x+CP*Temp));
  Pinv(1,2) = -v.x*v.y*BB/AA;
  Pinv(1,3) = -v.x*v.z*BB/AA;
  Pinv(1,4) = v.x*BB/AA;    
  Pinv(1,5) = -v.x*eta_fsd*BB/AA;
  Pinv(1,7) = -5.0*v.x*BB/(3.0*AA);
  Pinv(2,0) = v.y*CC*BB/AA;
  Pinv(2,1) = -v.x*v.y*BB/AA;
  Pinv(2,2) = rho*Rmix/AA*(pt+v.y*v.y*rho-theta*pt*(v.y*v.y+CP*Temp));
  Pinv(2,3) = -v.z*v.y*BB/AA;
  Pinv(2,4) = v.y*BB/AA;  
  Pinv(2,5) = -v.y*eta_fsd*BB/AA;
  Pinv(2,7) = -5.0*v.y*BB/(3.0*AA);
  Pinv(3,0) = v.z*CC*BB/AA;
  Pinv(3,1) = -v.x*v.z*BB/AA;
  Pinv(3,2) = -v.y*v.z*BB/AA;
  Pinv(3,3) = rho*Rmix/AA*(pt+v.z*v.z*rho-theta*pt*(v.z*v.z+CP*Temp));
  Pinv(3,4) = v.z*BB/AA;  
  Pinv(3,5) = -v.z*eta_fsd*BB/AA;
  Pinv(3,7) = -5.0*v.z*BB/(3.0*AA);
  Pinv(4,0) = DD*CC*BB/AA;
  Pinv(4,1) = -v.x*DD*BB/AA;
  Pinv(4,2) = -v.y*DD*BB/AA;
  Pinv(4,3) = -v.z*DD*BB/AA;
  Pinv(4,4) = rho*Rmix/AA*(theta*pt*(DD-CP*Temp)-rho*DD+pt);
  Pinv(4,5) = -DD*eta_fsd*BB/AA;
  Pinv(4,7) = -5.0*DD*BB/(3.0*AA);
  Pinv(5,0) = C*CC*BB/AA;
  Pinv(5,1) = -C*v.x*BB/AA;
  Pinv(5,2) = -C*v.y*BB/AA;
  Pinv(5,3) = -C*v.z*BB/AA;
  Pinv(5,4) = C*BB/AA;
  Pinv(5,5) = 1.0 - C*eta_fsd*BB/AA;
  Pinv(5,7) = -5.0*C*BB/(3.0*AA);
  Pinv(6,0) = Fsd*CC*BB/AA;
  Pinv(6,1) = -Fsd*v.x*BB/AA;
  Pinv(6,2) = -Fsd*v.y*BB/AA;
  Pinv(6,3) = -Fsd*v.z*BB/AA;
  Pinv(6,4) = Fsd*BB/AA;
  Pinv(6,5) = -Fsd*eta_fsd*BB/AA;
  Pinv(6,6) = ONE;
  Pinv(6,7) = -5.0*Fsd*BB/(3.0*AA);
  Pinv(7,0) = k*CC*BB/AA;
  Pinv(7,1) = -k*v.x*BB/AA;
  Pinv(7,2) = -k*v.y*BB/AA;
  Pinv(7,3) = -k*v.z*BB/AA;
  Pinv(7,4) = k*BB/AA;
  Pinv(7,5) = -k*eta_fsd*BB/AA;
  Pinv(7,7) = ONE-5.0*k/(3.0*AA);

}

/*********************************************************************************************
 * LES3DFsd_pState::dWdU -- Transform Jacobin for primitive state to conservetive state      * 
 *********************************************************************************************/
void LES3DFsd_pState::dWdU(DenseMatrix &dWdU){

  dWdU.zero();
  double alpha = Rtot()/(Cp()-Rtot());
  double eta_fsd = Progvar_Species_Grad();      

  dWdU(0,0) = ONE;
  dWdU(1,0) = -v.x/rho;
  dWdU(1,1) = ONE/rho;
  dWdU(2,0) = -v.y/rho;
  dWdU(2,2) = -ONE/rho;
  dWdU(3,0) = -v.z/rho;
  dWdU(3,3) = -ONE/rho;
  dWdU(4,0) = alpha*(v.sqr()/TWO-h()+Cp()*T()+C*eta_fsd);
  dWdU(4,1) = -v.x*alpha;
  dWdU(4,2) = -v.y*alpha;
  dWdU(4,3) = -v.z*alpha;
  dWdU(4,4) = alpha;
  dWdU(4,5) = -eta_fsd*alpha;
  dWdU(4,7) = -FIVE*alpha/THREE;
  dWdU(5,0) = -C/rho;
  dWdU(5,5) = ONE/rho;
  dWdU(6,0) = -Fsd/rho;
  dWdU(6,6) = ONE/rho;
  dWdU(7,0) = -k/rho;
  dWdU(7,7) = ONE/rho;

}

/*********************************************************************************
 * LES3DFsd_pState::SemiImplicitSourceJacobi -- Source term Jacobin for          *
 *                                              Semi-Implicit time marching      * 
 *********************************************************************************/
void LES3DFsd_pState::SemiImplicitSourceJacobi(const LES3DFsd_pState &dWdx, 
                                               const LES3DFsd_pState &dWdy,
                                               const LES3DFsd_pState &dWdz,
                                               const double &d_dWdx_dW, 
                                               const double &d_dWdy_dW,
                                               const double &d_dWdz_dW,
                                               DenseMatrix &dStdW,
                                               const int Flow_Type,
                                               const double &Volume) {

     dStdW.zero();
     double filter = filter_width(Volume);
     double k_fsd = SFS_Kinetic_Energy_Fsd(dWdx,dWdy,dWdz,Flow_Type,Volume); 
     double kappa_fsd = Efficiency_Function_Fsd(dWdx,dWdy,dWdz,Flow_Type,Volume);   
     double tau_fsd = HeatRelease_Parameter();
     double Cv = 0.086, Cs = 0.018;
     double beta_fsd = 1.0;

     double t1,t4,t5,t6,t7,t8,t9;
     double t10,t11,t12,t14,t15,t16,t17,t18,t19;
     double t22,t23,t24,t27,t28;
     double t30,t31,t32,t33,t34,t36,t37,t38,t39;
     double t41,t42,t43,t44,t45,t47;
     double t50,t51,t52,t53,t54,t55,t59;
     double t60,t61,t65,t66,t67;
     double t73,t77;
     double t81,t82,t83,t86,t87,t89;
     double t90,t95,t97;
     double t100,t106,t108;
     double t111,t117,t119;
     double t122,t128,t129;
     double t131,t132,t135,t137,t139;
     double t140,t142,t143,t145,t146,t148,t149;
     double t151,t152,t157,t158;
     double t165,t167;
     double t184,t188;
     double t190,t193,t195,t198;
     double t200,t207;
     double t246,t248;
     double t254,t258;
     double t262;
     double t279;
     double t291,t295,t296;
     double t301,t302,t304,t308;
     double t317;
     double t321,t326;
     double t334,t336;
     double t347;
     double t352;
     double t364;
     double t376,t379;
     double t390,t394,t398; 
     double t403;

      t1 = _reactants_density*_laminar_flame_speed;
      t4 = dWdx.C;//cx(c);
      t5 = t4*t4;
      t6 = dWdy.C;//cy(c);
      t7 = t6*t6;
      t8 = dWdz.C;//cz(c);
      t9 = t8*t8;
      t10 = t5+t7+t9;
      t11 = 1/t10;
      t12 = t5*t11;
      t14 = t7*t11;
      t15 = t14/3.0;
      t16 = t9*t11;
      t17 = t16/3.0;
      t18 = 2.0/3.0-2.0/3.0*t12+t15+t17;
      t19 = dWdx.v.x;//Ux(U);
      t22 = t12/3.0;
      t23 = 2.0/3.0-2.0/3.0*t14+t22+t17;
      t24 = dWdy.v.y;//Vy(V);
      t27 = 2.0/3.0-2.0/3.0*t16+t22+t15;
      t28 = dWdz.v.z;//Gz(G);
      t30 = t4*t11;
      if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY){
      t31 = dWdy.v.x;//Uy(U);
      t32 = dWdx.v.y;//Vx(V);
      }else{
      t31 = dWdx.v.y;//Vx(V);
      t32 = dWdy.v.x;//Uy(U);
      }
      t33 = t31+t32;
      t34 = t6*t33;
      if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY){
      t36 = dWdz.v.x;//Uz(U);
      t37 = dWdx.v.z;//Gx(G);
      }else{
      t36 = dWdx.v.z;//Gx(G);
      t37 = dWdz.v.x;//Uz(U);
      }
      t38 = t36+t37;
      t39 = t8*t38;
      t41 = t6*t11;
      if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY){
      t42 = dWdz.v.y;//Vz(V);
      t43 = dWdy.v.z;//Gy(G);
      }else{
      t42 = dWdy.v.z;//Gy(G);
      t43 = dWdz.v.y;//Vz(V);
      }
      t44 = t42+t43;
      t45 = t8*t44;
      t47 = t18*t19+t23*t24+t27*t28-t30*t34-t30*t39-t41*t45;
      t50 = 1.0+tau_fsd*C;//c;
      t51 = sqrt(t10);
      t52 = 1/t51;
      t53 = t4*t52;
      t54 = dWdx.Fsd;//Fsdx(Fsd);
      t55 = d_dWdx_dW;//diff(rhox(rho),rho);
      t59 = t6*t52;
      t60 = dWdy.Fsd;//Fsdy(Fsd);
      t61 = d_dWdy_dW;//diff(rhoy(rho),rho);
      t65 = t8*t52;
      t66 = dWdz.Fsd;//Fsdz(Fsd);
      t67 = d_dWdz_dW;//diff(rhoz(rho),rho);
      t73 = tau_fsd*Fsd;
      t77 = -t5*t52-t7*t52-t9*t52;
      t81 = sqrt(k_fsd);
      t82 = kappa_fsd*t81;
      t83 = 1/filter;
      t86 = beta_fsd*_laminar_flame_speed;
      t87 = Fsd*Fsd;
      t89 = 1.0-C;
      t90 = 1/t89;
      t95 = d_dWdx_dW;//diff(Ux(U),U);
      t97 = d_dWdy_dW;//diff(Uy(U),U);
      t100 = d_dWdz_dW;//diff(Uz(U),U);
      t106 = d_dWdy_dW;//diff(Vy(V),V);
      t108 = d_dWdx_dW;//diff(Vx(V),V);
      t111 = d_dWdz_dW;//diff(Vz(V),V);
      t117 = d_dWdz_dW;//diff(Gz(G),G);
      t119 = d_dWdx_dW;//diff(Gx(G),G);
      t122 = d_dWdy_dW;//diff(Gy(G),G);
      t128 = d_dWdx_dW;//diff(cx(c),c);
      t129 = t30*t128;
      t131 = t10*t10;
      t132 = 1/t131;
      t135 = d_dWdy_dW;//diff(cy(c),c);
      t137 = d_dWdz_dW;//diff(cz(c),c);
      t139 = t4*t128+t6*t135+t8*t137;
      t140 = 2.0*t5*t132*t139;
      t142 = t41*t135;
      t143 = 2.0/3.0*t142;
      t145 = 2.0*t7*t132*t139;
      t146 = t145/3.0;
      t148 = t8*t11*t137;
      t149 = 2.0/3.0*t148;
      t151 = 2.0*t9*t132*t139;
      t152 = t151/3.0;
      t157 = 2.0/3.0*t129;
      t158 = t140/3.0;
      t165 = t128*t11;
      t167 = t4*t132;
      t184 = (-4.0/3.0*t129+2.0/3.0*t140+t143-t146+t149-t152)*t19+(-4.0/3.0*
t142+2.0/3.0*t145+t157-t158+t149-t152)*t24+(-4.0/3.0*t148+2.0/3.0*t151+t157-
t158+t143-t146)*t28-t165*t34+2.0*t167*t34*t139-t30*t135*t33-t165*t39+2.0*t167*
t39*t139-t30*t137*t38-t135*t11*t45+2.0*t6*t132*t45*t139-t41*t137*t44;
      t188 = dWdx.rho;//rhox(rho);
      t190 = rho*t54+Fsd*t188;
      t193 = dWdy.rho;//rhoy(rho);
      t195 = rho*t60+Fsd*t193;
      t198 = dWdz.rho;//rhoz(rho);
      t200 = rho*t66+Fsd*t198;
      t207 = 1/t51/t10;
      t246 = rho*rho;
      t248 = t89*t89;
      t254 = d_dWdx_dW;//diff(Fsdx(Fsd),Fsd);
      t258 = d_dWdy_dW;//diff(Fsdy(Fsd),Fsd);
      t262 = d_dWdz_dW;//diff(Fsdz(Fsd),Fsd);
      if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K){
      t279 = 1/t81;
      t291 = pow(k_fsd,0.15);
      t295 = Cv*t81;
      t296 = filter*t95;
      t301 = 0.2*t24;
      t302 = 0.2*t28;
      t304 = filter*(0.4*t19-t301-t302);
      t308 = 0.6666666667*rho*k_fsd;
      t317 = filter*t33;
      t321 = filter*t38;
      t326 = filter*t106;
      t334 = 0.2E1*t19;
      t336 = filter*(0.4*t24-t334-t302);
      t347 = filter*t44;
      t352 = filter*t117;
      t364 = filter*(0.4*t28-t334-t301);
      t376 = Cv*t279;
      t379 = 0.6666666667*rho;
      t390 = t33*t33;
      t394 = t38*t38;
      t398 = t44*t44;
      t403 = pow(k_fsd,0.5);
      }
      dStdW(5,0) = t1*Fsd;
      dStdW(5,6) = t1*rho;
      dStdW(6,0) = t47*Fsd-_laminar_flame_speed*(t50*(-t53*(t54+Fsd*t55)-t59*(t60+Fsd
*t61)-t65*(t66+Fsd*t67))+t73*t77)+t82*Fsd*t83-2.0*t86*t87*rho*t90;
      dStdW(6,1) = (t18*t95-t30*t6*t97-t30*t8*t100)*Fsd*rho;
      dStdW(6,2) = (t23*t106-t30*t6*t108-t41*t8*t111)*Fsd*rho;
      dStdW(6,3) = (t27*t117-t30*t8*t119-t41*t8*t122)*Fsd*rho;
      dStdW(6,5) = t184*Fsd*rho-_laminar_flame_speed*(tau_fsd*(-t53*t190-t59*t195-t65
*t200)+t50*(-t128*t52*t190+t4*t207*t190*t139-t135*t52*t195+t6*t207*t195*t139-
t137*t52*t200+t8*t207*t200*t139)+t73*rho*(-2.0*t53*t128+t5*t207*t139-2.0*t59*
t135+t7*t207*t139-2.0*t65*t137+t9*t207*t139))-t86*t87*t246/t248;
      dStdW(6,6) = t47*rho-_laminar_flame_speed*(t50*(-t53*(rho*t254+t188)-t59*(rho*
t258+t193)-t65*(rho*t262+t198))+tau_fsd*rho*t77)+t82*rho*t83-2.0*t86*Fsd*t246*
t90;
      if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K){
      dStdW(6,7) = kappa_fsd*t279*Fsd*rho*t83/2.0;
      dStdW(7,0) = -0.6666666667*k_fsd*t19-0.6666666667*k_fsd*t24-0.6666666667*k_fsd*t28-
Cs*t291*t83;
      dStdW(7,1) = 0.2666666667E1*t295*t296*t19+(0.6666666667*t295*t304-t308)*
t95-0.1333333333E1*t295*t296*t24-0.1333333333E1*t295*t296*t28+0.4E1*t295*t317*
t97+0.4E1*t295*t321*t100;
      dStdW(7,2) = -0.1333333333E1*t295*t326*t19+0.2666666667E1*t295*t326*t24+
(0.6666666667*t295*t336-t308)*t106-0.1333333333E1*t295*t326*t28+0.4E1*t295*t317
*t108+0.4E1*t295*t347*t111;
      dStdW(7,3) = -0.1333333333E1*t295*t352*t19-0.1333333333E1*t295*t352*t24+
0.2666666667E1*t295*t352*t28+(0.6666666667*t295*t364-t308)*t117+0.4E1*t295*t321
*t119+0.4E1*t295*t347*t122;
      dStdW(7,7) = (0.3333333334*t376*t304-t379)*t19+(0.3333333334*t376*t336-
t379)*t24+(0.3333333334*t376*t364-t379)*t28+0.1E1*t376*filter*t390+0.1E1*t376*
filter*t394+0.1E1*t376*filter*t398-0.15E1*Cs*rho*t403*t83;
      }

}

/************************************************************************
 *************** NUMERICAL FLUX FUNCTIONS *******************************
 ************************************************************************/

/**************************************************************************
 * LES3DFsd_pState::RoeAverage -- Roe-averaged primitive solution state.  *
 **************************************************************************/
LES3DFsd_pState LES3DFsd_pState::RoeAverage(const LES3DFsd_pState &Wl,
                                            const LES3DFsd_pState &Wr) {

   double Hl, Hr, srhol, srhor;
   double Ha, ha;
   LES3DFsd_pState Temp;

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
   Temp.C = (srhol*Wl.C+srhor*Wr.C)/(srhol+srhor);
   Temp.Fsd = (srhol*Wl.Fsd+srhor*Wr.Fsd)/(srhol+srhor);
   Temp.k = (srhol*Wl.k+srhor*Wr.k)/(srhol+srhor);

   /* Set the species mass fractions. */

   Temp.premixed_mfrac();

   /* Determine the pressure. */

   Ha = (srhol*Hl+srhor*Hr)/(srhol+srhor);
   ha = Ha - HALF*(sqr(Temp.v.x)+sqr(Temp.v.y)+sqr(Temp.v.z)) - FIVE*Temp.k/THREE;
   Temp.p = Temp.rho*Temp.T(ha)*Temp.Rtot();
   
   /* Return the Roe-averged state. */

   return (Temp);     

}

/**************************************************************************
 * LES3DFsd_pState::FluxHLLE_x -- HLLE flux function, x-direction flux.   *
 **************************************************************************/
LES3DFsd_cState  LES3DFsd_pState::FluxHLLE_x(const LES3DFsd_pState &Wl,
                                             const LES3DFsd_pState &Wr) {
   
   double wavespeed_l, wavespeed_r;
   LES3DFsd_pState Wa, lambdas_l, lambdas_r, lambdas_a;
   LES3DFsd_cState Flux, dUrl;

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
 * LES3DFsd_pState::FluxHLLE_n -- HLLE flux function, n-direction flux.   *
 **************************************************************************/
LES3DFsd_cState LES3DFsd_pState::FluxHLLE_n(const LES3DFsd_pState &Wl,
                                            const LES3DFsd_pState &Wr,
                                            const Vector3D &norm_dir) {

   // Determine the left and right solution states in the rotate frame.
   LES3DFsd_pState Wl_rot(Wl.Rotate(norm_dir));
   LES3DFsd_pState Wr_rot(Wr.Rotate(norm_dir));
    
   // Evaluate the intermediate state solution flux in the rotated frame.
   Wl_rot.premixed_mfrac();
   Wr_rot.premixed_mfrac();
   LES3DFsd_cState Flux_rot = FluxHLLE_x(Wl_rot, Wr_rot);
 
   // Return numerical flux in un-rotated frame.
   return (Flux_rot.RotateBack(norm_dir));

}

/**************************************************************************
 * LES3DFsd_pState::FluxRoe_x -- Roe flux function, x-direction flux.     *
 **************************************************************************/
LES3DFsd_cState LES3DFsd_pState::FluxRoe_x(const  LES3DFsd_pState &Wl,  
                                           const  LES3DFsd_pState &Wr) {
   
   LES3DFsd_pState Wa, dWrl, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
   LES3DFsd_cState Flux;

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
 * LES3DFsd_pState::FluxRoe_n -- Roe flux function, n-direction flux.     *
 **************************************************************************/
LES3DFsd_cState LES3DFsd_pState::FluxRoe_n(const LES3DFsd_pState &Wl,
                                           const LES3DFsd_pState &Wr,
                                           const Vector3D &norm_dir) {
   
   // Determine the left and right solution states in the rotate frame.
   LES3DFsd_pState Wl_rot(Wl.Rotate(norm_dir));
   LES3DFsd_pState Wr_rot(Wr.Rotate(norm_dir));
    
   // Evaluate the intermediate state solution flux in the rotated frame.
   Wl_rot.premixed_mfrac();
   Wr_rot.premixed_mfrac();
   LES3DFsd_cState Flux_rot = FluxRoe_x(Wl_rot, Wr_rot);
 
   // Return numerical flux in un-rotated frame.
   return (Flux_rot.RotateBack(norm_dir));

}

/**************************************************************************************
 * LES3DFsd_pState::AUSMplus_up_x -- AUSMplus_up flux function, x-direction flux.     *
 **************************************************************************************/
LES3DFsd_cState LES3DFsd_pState::FluxAUSMplus_up_x(const LES3DFsd_pState &Wl,
                                                   const LES3DFsd_pState &Wr) {
 
  LES3DFsd_cState Flux;
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
    Flux.rhoC = Wl.C;
    Flux.rhoFsd = Wl.Fsd;
    Flux.rhok = Wl.k;
  } else {
    Flux.rho = ONE;
    Flux.rhov.x = Wr.v.x; 
    Flux.rhov.y = Wr.v.y; 
    Flux.E = Wr.H()/Wr.rho;
    Flux.rhoC = Wr.C;
    Flux.rhoFsd = Wr.Fsd;
    Flux.rhok = Wr.k;
  } /* endif */

  Flux = mass_flux_half*Flux;
 
  // Add the pressure contribution to the intermediate state solution flux:
  Flux[2] += phalf;
  
  // Return solution flux.
  return Flux;

}

/******************************************************************************************
 * LES3DFsd_pState::FluxAUSMplus_up_n -- AUSMplus_up flux function, n-direction flux.     *
 ******************************************************************************************/
LES3DFsd_cState LES3DFsd_pState::FluxAUSMplus_up_n(const LES3DFsd_pState &Wl,
                                                   const LES3DFsd_pState &Wr,
                                                   const Vector3D &norm_dir) {
   
   // Determine the left and right solution states in the rotate frame.
   LES3DFsd_pState Wl_rot(Wl.Rotate(norm_dir));
   LES3DFsd_pState Wr_rot(Wr.Rotate(norm_dir));
    
   // Evaluate the intermediate state solution flux in the rotated frame.
   Wl_rot.premixed_mfrac();
   Wr_rot.premixed_mfrac();
   LES3DFsd_cState Flux_rot = FluxAUSMplus_up_x(Wl_rot, Wr_rot);
 
   // Return numerical flux in un-rotated frame.
   return (Flux_rot.RotateBack(norm_dir));

}

/**************************************************************************
 * LES3DFsd_pState::lambda_minus -- Negative wave speeds determined using *
 *                                  Harten entropy fix.                   *
 **************************************************************************/
LES3DFsd_pState LES3DFsd_pState::lambda_minus(const LES3DFsd_pState &lambdas_a,
                                              const LES3DFsd_pState &lambdas_l,
                                              const LES3DFsd_pState &lambdas_r) {

   LES3DFsd_pState W_temp;

   W_temp.rho = HartenFixNeg(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
   W_temp.v.x = HALF*(lambdas_a[2]-fabs(lambdas_a[2]));
   W_temp.v.y = HALF*(lambdas_a[3]-fabs(lambdas_a[3]));
   W_temp.v.z = HALF*(lambdas_a[4]-fabs(lambdas_a[4]));
   W_temp.p = HartenFixNeg(lambdas_a[5],lambdas_l[5],lambdas_r[5]);
   W_temp.C = HALF*(lambdas_a[6] -fabs(lambdas_a[6]));
   W_temp.Fsd = HALF*(lambdas_a[7] - fabs(lambdas_a[7]));
   W_temp.k = HALF*(lambdas_a[8] -fabs(lambdas_a[8]));
   
   return (W_temp);

}

/**************************************************************************
 * LES3DFsd_pState::lambda_plus -- Positive wave speeds determined using  *
 *                                 Harten entropy fix.                    *
 **************************************************************************/
LES3DFsd_pState LES3DFsd_pState::lambda_plus(const LES3DFsd_pState &lambdas_a,
                                             const LES3DFsd_pState &lambdas_l,
                                             const LES3DFsd_pState &lambdas_r) {
   
   LES3DFsd_pState W_temp;

   W_temp.rho = HartenFixPos(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
   W_temp.v.x = HALF*(lambdas_a[2]+fabs(lambdas_a[2]));
   W_temp.v.y = HALF*(lambdas_a[3]+fabs(lambdas_a[3]));
   W_temp.v.z = HALF*(lambdas_a[4]+fabs(lambdas_a[4]));
   W_temp.p = HartenFixPos(lambdas_a[5],lambdas_l[5],lambdas_r[5]);
   W_temp.C = HALF*(lambdas_a[6] + fabs(lambdas_a[6]));
   W_temp.Fsd = HALF*(lambdas_a[7] +fabs(lambdas_a[7]));
   W_temp.k = HALF*(lambdas_a[8] + fabs(lambdas_a[8]));

   return (W_temp);

}

/*************************************************************************************************
 * LES3DFsd_pState::HLLE_wavespeeds -- Returns the lambda plus and lambda minus wave speeds for  *
 *                                     rotated Riemann problem aligned with norm_dir given       * 
 *                                     unroated solution states Wl and Wr.                       *
 ************************************************************************************************/
Vector2D LES3DFsd_pState::HLLE_wavespeeds(const LES3DFsd_pState &Wl,
                                          const LES3DFsd_pState &Wr,
                                          const Vector3D &norm_dir) {

    Vector2D wavespeed;
    LES3DFsd_pState Wa_n, lambdas_l, lambdas_r, lambdas_a;  //Lots of TEMPS

    /* Determine the left and right states. */

    LES3DFsd_pState Wl_rotated(Wl.Rotate(norm_dir));
    LES3DFsd_pState Wr_rotated(Wr.Rotate(norm_dir));
    Wl_rotated.premixed_mfrac();
    Wr_rotated.premixed_mfrac();

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
 * LES3DFsd_pState::Rotate -- Returns a rotated primitive state aligned   *
 *                            with a local x-axis in the norm_dir.        *
 **************************************************************************/
LES3DFsd_pState LES3DFsd_pState::Rotate(const Vector3D &norm_dir) const {

  // for a 3D unit normal rotated to align with the x-axis
  double Ct = norm_dir.x;  //cos_angle
  double St = sqrt( norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z); //sin_angle
  Vector3D rt(0,norm_dir.z,-norm_dir.y);  //rotation axis  

  return LES3DFsd_pState(rho,
                         v.x*Ct - v.y*rt.z*St + v.z*rt.y*St,
                         v.x*rt.z*St + v.y*(rt.y*rt.y*(ONE-Ct)+Ct) + v.z*(rt.y*rt.z*(ONE-Ct)),
                         - v.x*rt.y*St + v.y*(rt.y*rt.z*(ONE-Ct)) + v.z*(rt.z*rt.z*(ONE-Ct)+Ct),
                         p,
                         C,
                         Fsd,  
                         k);

}

/**************************************************************************
 * LES3DFsd_pState::Rotate -- Returns an un-rotated primitive state       *
 *                            re-alinged from the x-axis of the global    *
 *                            problem.                                    *
 **************************************************************************/
LES3DFsd_pState LES3DFsd_pState::RotateBack(const Vector3D &norm_dir) const {

  // for a 3D unit normal rotated to align with the x-axis
  double Ct = norm_dir.x;  //cos_angle
  double St = sqrt( norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z); //sin_angle
  Vector3D rt(0,norm_dir.z,-norm_dir.y);  //rotation axis  

  return LES3DFsd_pState(rho,
                         v.x*Ct + v.y*rt.z*St - v.z*rt.y*St,
                         - v.x*rt.z*St + v.y*(rt.y*rt.y*(ONE-Ct)+Ct) + v.z*(rt.y*rt.z*(ONE-Ct)),
                         + v.x*rt.y*St + v.y*(rt.y*rt.z*(ONE-Ct)) + v.z*(rt.z*rt.z*(ONE-Ct)+Ct),
                         p,
                         C,
                         Fsd,  
                         k);

}

/************************************************************************
 *************** NUMERICAL EVALUATION OF VISCOUS FLUXES *****************
 ************************************************************************/

/**********************************************************************
 * LES3DFsd_pState::FluxViscous_n -- Viscous flux (n-direction).      *
/**********************************************************************/
LES3DFsd_cState LES3DFsd_pState::FluxViscous_n(const LES3DFsd_pState &Wl,
                                               const LES3DFsd_pState &Wr,
                                               const LES3DFsd_pState &Wc,
                                               const LES3DFsd_pState &Wc_Neigbor,
                                               const LES3DFsd_pState &dWdx,
                                               const LES3DFsd_pState &dWdy,
                                               const LES3DFsd_pState &dWdz,
                                               const LES3DFsd_pState &dWdx_Neigbor,
                                               const LES3DFsd_pState &dWdy_Neigbor,
                                               const LES3DFsd_pState &dWdz_Neigbor,
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
   
   LES3DFsd_pState dWdx_Weighted, dWdy_Weighted, dWdz_Weighted, 
                   dWdx_face, dWdy_face, dWdz_face, Grad_middle_term;

   LES3DFsd_pState W_face;
   
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

   // Ensure that we have correct species mass fractions
   W_face.premixed_mfrac();

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
 * LES3DFsd_pState::Reflect -- Return reflected solution state.           *
 **************************************************************************/
LES3DFsd_pState LES3DFsd_pState::Reflect(const LES3DFsd_pState &W,
                                         const Vector3D &norm_dir) {
   
   Vector3D ur_norm, ur_tang, vr_tot;
   LES3DFsd_pState Temp; Temp.Copy(W);

   ur_norm = dot(W.v, norm_dir)*norm_dir;
   ur_tang = W.v - ur_norm;

   ur_norm = -ur_norm;
   vr_tot = ur_norm + ur_tang;

   Temp.v = vr_tot;
   
   return (Temp);
       
}

/**************************************************************************
 * LES3DFsd_pState::MovingWall -- Return moving wall boundary state.      *
 **************************************************************************/
LES3DFsd_pState LES3DFsd_pState::MovingWall(const LES3DFsd_pState &Win,
                                            const LES3DFsd_pState &Wout,
                                            const Vector3D &norm_dir, 
                                            const Vector3D &wall_velocity,
                                            const Vector3D &pressure_gradient,
                                            const int &TEMPERATURE_BC_FLAG) {
   
   LES3DFsd_pState Temp; Temp.Copy(Win);
   
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
 * LES3DFsd_pState::NoSlip -- Return no-slip wall boundary state.         *
 **************************************************************************/
LES3DFsd_pState LES3DFsd_pState::NoSlip(const LES3DFsd_pState &Win,
                                        const LES3DFsd_pState &Wout,
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
 * LES3DFsd_pState::abs_strain_rate -- Absolute value of strain rate tensor      *
 *********************************************************************************/
double LES3DFsd_pState::abs_strain_rate(const LES3DFsd_pState &dWdx, 
                                        const LES3DFsd_pState &dWdy, 
                                        const LES3DFsd_pState &dWdz) const{

   Tensor3D strain_rate; 
   strain_rate.zero();
   strain_rate.xx = dWdx.v.x;
   strain_rate.yy = dWdy.v.y;
   strain_rate.zz = dWdz.v.z;
   strain_rate.xy = 0.5*(dWdy.v.x + dWdx.v.y);
   strain_rate.yz = 0.5*(dWdz.v.y + dWdy.v.z);
   strain_rate.xz = 0.5*(dWdx.v.z + dWdz.v.x);
  // S[i,j]*S[i,j]
  double SS = sqr(strain_rate.xx) + sqr(strain_rate.yy) + sqr(strain_rate.zz)
            + 2.0*(sqr(strain_rate.xy) + sqr(strain_rate.yz) + sqr(strain_rate.xz));
  // sqrt(2*S*S)
  return sqrt(2.0*SS);

}

/*********************************************************************************
 * LES3DFsd_pState::grad_abs_strain_rate -- gradient of strain rate tensor       *
 *********************************************************************************/
Vector3D LES3DFsd_pState::grad_abs_strain_rate(const LES3DFsd_pState &dWdx, 
                                               const LES3DFsd_pState &dWdy, 
                                               const LES3DFsd_pState &dWdz,
                                               const LES3DFsd_pState &d_dWdx_dx,
                                               const LES3DFsd_pState &d_dWdy_dy,
                                               const LES3DFsd_pState &d_dWdz_dz,
                                               const LES3DFsd_pState &d_dWdx_dy,
                                               const LES3DFsd_pState &d_dWdx_dz,
			  	               const LES3DFsd_pState &d_dWdy_dz) const {

  Vector3D temp;

  temp.zero();
  temp.x = dWdx.v.x*d_dWdx_dx.v.x+dWdy.v.y*d_dWdx_dy.v.y+dWdz.v.z*d_dWdx_dz.v.z
            +(dWdy.v.x+dWdx.v.y)*(d_dWdx_dy.v.x+d_dWdx_dx.v.y)
            +(dWdz.v.y+dWdy.v.z)*(d_dWdx_dz.v.y+d_dWdx_dy.v.z) 
            +(dWdx.v.z+dWdz.v.x)*(d_dWdx_dx.v.z+d_dWdx_dz.v.x);
  temp.y = dWdx.v.x*d_dWdx_dy.v.x+dWdy.v.y*d_dWdy_dy.v.y+dWdz.v.z*d_dWdy_dz.v.z
            +(dWdy.v.x+dWdx.v.y)*(d_dWdy_dy.v.x+d_dWdx_dy.v.y)
            +(dWdz.v.y+dWdy.v.z)*(d_dWdy_dz.v.y+d_dWdy_dy.v.z) 
            +(dWdx.v.z+dWdz.v.x)*(d_dWdx_dy.v.z+d_dWdy_dz.v.x);
  temp.z = dWdx.v.x*d_dWdx_dz.v.x+dWdy.v.y*d_dWdy_dz.v.y+dWdz.v.z*d_dWdz_dz.v.z
            +(dWdy.v.x+dWdx.v.y)*(d_dWdy_dz.v.x+d_dWdx_dz.v.y)
            +(dWdz.v.y+dWdy.v.z)*(d_dWdz_dz.v.y+d_dWdy_dz.v.z) 
            +(dWdx.v.z+dWdz.v.x)*(d_dWdx_dz.v.z+d_dWdz_dz.v.x);

  return (temp);

}

/*********************************************************************************
 * LES3DFsd_pState::HeatRelease_Parameter -- Heat release parameter              *
 *********************************************************************************/
double LES3DFsd_pState::HeatRelease_Parameter(void) const {
  return (_adiabatic_flame_temperature/298.0-1.0);
}

/*********************************************************************************
 * LES3DFsd_pState::SFS_Kinetic_Energy_Fsd -- Subfilter scale kinetic energy     *
 *********************************************************************************/
double LES3DFsd_pState::SFS_Kinetic_Energy_Fsd(const LES3DFsd_pState &dWdx,
                                               const LES3DFsd_pState &dWdy,
                                               const LES3DFsd_pState &dWdz,
                                               const int Flow_Type,
                                               const double &Volume) {

  if ( Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY ) {
    double CI = 0.005;
    return (CI*sqr(filter_width(Volume)*abs_strain_rate(dWdx,dWdy,dWdz)));
  } else if ( Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K ) {
    return (k);
  } /* endif */

}

/**************************************************************************************************
 * LES3DFsd_pState::Efficiency_Function_Fsd -- Efficiency function for subfilter scale strain     *
 **************************************************************************************************/
double LES3DFsd_pState::Efficiency_Function_Fsd(const LES3DFsd_pState &dWdx,
                                                const LES3DFsd_pState &dWdy,
                                                const LES3DFsd_pState &dWdz,
                                                const int Flow_Type,
                                                const double &Volume) {

  double filter, kappa_fsd, k_fsd;
  k_fsd = SFS_Kinetic_Energy_Fsd(dWdx,dWdy,dWdz,Flow_Type,Volume);
  filter = filter_width(Volume);
  kappa_fsd = 0.75*exp(-1.2/pow(sqrt(k_fsd)/_laminar_flame_speed,0.3))*
              pow(filter/_laminar_flame_thickness,2.0/3.0);
  return(kappa_fsd);

}

/********************************************************************************************************
 * LES3DFsd_pState::Progvar_Species_Grad -- Gradient of progress variable to species mass fractions     *
 ********************************************************************************************************/
double LES3DFsd_pState::Progvar_Species_Grad(void) const {

  double Temp, stoich, ratio, f_ub, eta_fsd;

  Temp = p/(rho*Rtot());

  if ( React.reactset_flag == CH4_1STEP ||
       React.reactset_flag == C3H8_1STEP ) {
     if ( React.reactset_flag == CH4_1STEP ){
       stoich = 2.0*specdata[1].Mol_mass()/specdata[0].Mol_mass();
       ratio = specdata[2].Mol_mass()/(specdata[2].Mol_mass()+2.0*specdata[3].Mol_mass());
       f_ub = specdata[0].Mol_mass()/(specdata[0].Mol_mass()+2.0*specdata[1].Mol_mass()+7.52*specdata[4].Mol_mass());
     }else if ( React.reactset_flag == C3H8_1STEP ){
       stoich = 5.0*specdata[1].Mol_mass()/specdata[0].Mol_mass();
       ratio = 3.0*specdata[2].Mol_mass()/(3.0*specdata[2].Mol_mass()+4.0*specdata[3].Mol_mass());
       f_ub = specdata[0].Mol_mass()/(specdata[0].Mol_mass()+5.0*specdata[1].Mol_mass()+18.8*specdata[4].Mol_mass());
     } /* endif */
    eta_fsd = (specdata[0].Enthalpy(Temp)+specdata[0].Heatofform()-Cp(Temp)*Temp*specdata[0].Rs()/Rtot())*(-f_ub)
             +(specdata[1].Enthalpy(Temp)+specdata[1].Heatofform()-Cp(Temp)*Temp*specdata[1].Rs()/Rtot())*(-stoich*f_ub/_fuel_equivalence_ratio)
	     +(specdata[2].Enthalpy(Temp)+specdata[2].Heatofform()-Cp(Temp)*Temp*specdata[2].Rs()/Rtot())*((1.0+stoich/_fuel_equivalence_ratio)*f_ub*ratio)
	     +(specdata[3].Enthalpy(Temp)+specdata[3].Heatofform()-Cp(Temp)*Temp*specdata[3].Rs()/Rtot())*((1.0+stoich/_fuel_equivalence_ratio)*f_ub*(1.0-ratio));
  }else if ( React.reactset_flag == H2O2_1STEP ){
    stoich = specdata[1].Mol_mass()/2.0/specdata[0].Mol_mass();
    f_ub = 2.0*specdata[0].Mol_mass()/(2.0*specdata[0].Mol_mass()+specdata[1].Mol_mass()+3.76*specdata[3].Mol_mass());
    eta_fsd = (specdata[0].Enthalpy(Temp)+specdata[0].Heatofform()-Cp(Temp)*Temp*specdata[0].Rs()/Rtot())*(-f_ub)
             +(specdata[1].Enthalpy(Temp)+specdata[1].Heatofform()-Cp(Temp)*Temp*specdata[1].Rs()/Rtot())*(-stoich*f_ub/_fuel_equivalence_ratio)
             +(specdata[2].Enthalpy(Temp)+specdata[2].Heatofform()-Cp(Temp)*Temp*specdata[2].Rs()/Rtot())*((1.0+stoich/_fuel_equivalence_ratio)*f_ub);
  } /* endif */

  return (eta_fsd);

}

/***********************************************************
 * LES3DFsd_pState::Reaction_Rate_Fsd -- Reaction rate     *
 ***********************************************************/
double LES3DFsd_pState::Reaction_Rate_Fsd(const LES3DFsd_pState &dWdx,
                                          const LES3DFsd_pState &dWdy,
                                          const LES3DFsd_pState &dWdz) {

   double tau_fsd = HeatRelease_Parameter();

   if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO ) {
     return ( _reactants_density*_laminar_flame_speed*Fsd*rho ); //-tau_fsd*_laminar_flame_speed*(rho*(1-2*C)*(dWdx.C+dWdy.C+dWdz.C)+C*(1-C)*(dWdx.rho+dWdy.rho+dWdz.rho)) );
   } else {
     return (ZERO);
   } /* endif */

}

/********************************************************************
 * LES3DFsd_pState::M_x -- x-direction surface averaged normal      *
 ********************************************************************/
double LES3DFsd_pState::M_x(const LES3DFsd_pState &dWdx,
                            const LES3DFsd_pState &dWdy,
                            const LES3DFsd_pState &dWdz) const {
   double Mx, tau_fsd;
   tau_fsd = HeatRelease_Parameter();

   if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO) { //&& dWdy.C != ZERO && dWdz.C != ZERO ) {//&& Fsd != ZERO ) {
     Mx = dWdx.C/sqrt(sqr(dWdx.C)+sqr(dWdy.C)+sqr(dWdz.C));
//  Mx = -((1.0+tau_fsd)*(1-exp(-0.2*1.28))/sqr(1.0+tau_fsd*C)+exp(-0.2*1.28))*dWdx.C/Fsd/rho;
//  Mx = -(1.0+tau_fsd)*dWdx.C/sqr(1.0+tau_fsd*C)/Fsd/rho;
//  Mx = -dWdx.C/Fsd/rho;
// if ( Mx < -1.0 ) { Mx = -1.0; }
     return Mx;
   }else {
    return (0.0);
   } /* endif */

}

/********************************************************************
 * LES3DFsd_pState::M_y -- y-direction surface averaged normal      *
 ********************************************************************/
double LES3DFsd_pState::M_y(const LES3DFsd_pState &dWdx,
                            const LES3DFsd_pState &dWdy,
                            const LES3DFsd_pState &dWdz) const {

   double My, tau_fsd;
   tau_fsd = HeatRelease_Parameter();
   if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ){//&& dWdy.C != ZERO && dWdz.C != ZERO && Fsd != ZERO ) {
     My = -dWdy.C/sqrt(sqr(dWdx.C)+sqr(dWdy.C)+sqr(dWdz.C));
//  My = -((1.0+tau_fsd)*(1-exp(-0.2*1.28))/sqr(1.0+tau_fsd*C)+exp(-0.2*1.28))*dWdy.C/Fsd/rho;
//  My = -(1.0+tau_fsd)*dWdy.C/sqr(1.0+tau_fsd*C)/Fsd/rho;
//  My = -dWdy.C/Fsd/rho;
//  if ( My < -1.0 ) { My = -1.0; }
     return My;
   } else {
     return (0.0);
   } /* endif */

}

/********************************************************************
 * LES3DFsd_pState::M_z -- z-direction surface averaged normal      *
 ********************************************************************/
double LES3DFsd_pState::M_z(const LES3DFsd_pState &dWdx,
                            const LES3DFsd_pState &dWdy,
                            const LES3DFsd_pState &dWdz) const {

    double Mz, tau_fsd;
    tau_fsd = HeatRelease_Parameter();

    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ){//&& dWdy.C != ZERO && dWdz.C != ZERO && Fsd != ZERO ) {
      Mz = -dWdz.C/sqrt(sqr(dWdx.C)+sqr(dWdy.C)+sqr(dWdz.C));
//  Mx = -((1.0+tau_fsd)*(1-exp(-0.2*1.28))/sqr(1.0+tau_fsd*C)+exp(-0.2*1.28))*dWdx.C/Fsd/rho;
//  Mx = -(1.0+tau_fsd)*dWdx.C/sqr(1.0+tau_fsd*C)/Fsd/rho;
//  Mx = -dWdx.C/Fsd/rho;
//  if ( Mx < -1.0 ) { Mx = -1.0; }
      return Mz;
    }else {
      return (0.0);
   } /* endif */

}

/***************************************************************
 * LES3DFsd_pState::Resolved_Strain -- Resolved strain term    *
 ***************************************************************/
double LES3DFsd_pState::Resolved_Strain(const LES3DFsd_pState &dWdx,
                                        const LES3DFsd_pState &dWdy,
                                        const LES3DFsd_pState &dWdz) const {

   double Mx, My, Mz, n_xx, n_yy, n_zz, n_xy, n_xz, n_yz, alpha_fsd;
   double resolved_strain_xx, resolved_strain_yy, resolved_strain_zz, resolved_strain_xy, resolved_strain_xz, resolved_strain_yz;

   Mx = M_x(dWdx,dWdy,dWdz);
   My = M_y(dWdx,dWdy,dWdz);
   Mz = M_z(dWdx,dWdy,dWdz);
   alpha_fsd = ONE - sqr(Mx) - sqr(My) - sqr(Mz);
   n_xx = sqr(Mx)+ONE/THREE*alpha_fsd;
   n_yy = sqr(My)+ONE/THREE*alpha_fsd;
   n_zz = sqr(Mz)+ONE/THREE*alpha_fsd;
   n_xy = Mx*My;
   n_xz = Mx*Mz;
   n_yz = My*Mz;

   if ( C < 0.999 && C > 0.001 &&  dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO ) {
     resolved_strain_xx = (ONE - n_xx)*dWdx.v.x*Fsd*rho;
     resolved_strain_yy = (ONE - n_yy)*dWdy.v.y*Fsd*rho;
     resolved_strain_zz = (ONE - n_zz)*dWdz.v.z*Fsd*rho;
     resolved_strain_xy = -n_xy*(dWdx.v.y + dWdy.v.x)*Fsd*rho;
     resolved_strain_xz = -n_xz*(dWdx.v.z + dWdz.v.x)*Fsd*rho;
     resolved_strain_yz = -n_yz*(dWdz.v.y + dWdy.v.y)*Fsd*rho;
     return (resolved_strain_xx + resolved_strain_yy + resolved_strain_zz + 
             resolved_strain_xy + resolved_strain_xz + resolved_strain_yz);
   } else {
     return (ZERO);
   } /* endif */

}

/**************************************************************************************************
 * LES3DFsd_pState::Resolved_Propagation_Curvature -- Resolved propagation and curvature term     *
 **************************************************************************************************/
double LES3DFsd_pState::Resolved_Propagation_Curvature(const LES3DFsd_pState &dWdx,
                                                       const LES3DFsd_pState &dWdy,
                                                       const LES3DFsd_pState &dWdz) const {
   double tau_fsd, Mx, My, Mz, 
          resolved_propagation_curvature_x, resolved_propagation_curvature_y, resolved_propagation_curvature_z;

   tau_fsd = HeatRelease_Parameter();
   Mx = M_x(dWdx,dWdy,dWdz);
   My = M_y(dWdx,dWdy,dWdz);
   Mz = M_z(dWdx,dWdy,dWdz);
   if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO ) {
      resolved_propagation_curvature_x = -_laminar_flame_speed*(ONE+tau_fsd*C)*Mx*(rho*dWdx.Fsd+Fsd*dWdx.rho)-
                                          _laminar_flame_speed*tau_fsd*Fsd*rho*Mx*dWdx.C;
      resolved_propagation_curvature_y = -_laminar_flame_speed*(ONE+tau_fsd*C)*My*(rho*dWdy.Fsd+Fsd*dWdy.rho)-
                                          _laminar_flame_speed*tau_fsd*Fsd*rho*My*dWdy.C;
      resolved_propagation_curvature_z = -_laminar_flame_speed*(ONE+tau_fsd*C)*Mz*(rho*dWdz.Fsd+Fsd*dWdz.rho)-
                                          _laminar_flame_speed*tau_fsd*Fsd*rho*Mz*dWdz.C;
      return ( resolved_propagation_curvature_x + resolved_propagation_curvature_y + resolved_propagation_curvature_z );
    } else {
      return (ZERO);
    } /* endif */

}

/******************************************************************
 * LES3DFsd_pState::SFS_Strain -- Subfilter scale strain term     *
 ******************************************************************/
double LES3DFsd_pState::SFS_Strain(const LES3DFsd_pState &dWdx,
                                   const LES3DFsd_pState &dWdy,
                                   const LES3DFsd_pState &dWdz,
                                   const int Flow_Type,
                                   const double &Volume) {

   double filter, kappa_fsd, k_fsd;
   filter = filter_width(Volume);
   k_fsd = SFS_Kinetic_Energy_Fsd(dWdx,dWdy,dWdz,Flow_Type,Volume); 
   kappa_fsd = Efficiency_Function_Fsd(dWdx,dWdy,dWdz,Flow_Type,Volume);   

   if (C < 0.999 && C > 0.001) {
     return ( kappa_fsd*sqrt(k_fsd)*Fsd*rho/filter );
   } else {
    return (0.0);
   } /* endif */

}

/***********************************************************************
 * LES3DFsd_pState::SFS_Curvature -- Subfilter scale curvature term    *
 ***********************************************************************/
double LES3DFsd_pState::SFS_Curvature(const LES3DFsd_pState &dWdx,
                                      const LES3DFsd_pState &dWdy,
                                      const LES3DFsd_pState &dWdz) const {

    double Mx, My, Mz, alpha_fsd, beta_fsd;
    beta_fsd = 1.0;
    Mx = M_x(dWdx,dWdy,dWdz);
    My = M_y(dWdx,dWdy,dWdz);
    Mz = M_z(dWdx,dWdy,dWdz);
    alpha_fsd = ONE - sqr(Mx) - sqr(My) - sqr(Mz);

    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO && Fsd != ZERO) {
      return ( -beta_fsd*_laminar_flame_speed*sqr(Fsd*rho)/(ONE-C) );

//     if ( Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_CHARLETTE ){
//     double tau_fsd = HeatRelease_Parameter();
//     double c_bar = (1.0+tau_fsd)*C/(1.0+tau_fsd*C);
//     return(-beta_fsd*_laminar_flame_speed*(Fsd-(1+tau_fsd)*sqrt(sqr(dWdx.C)+sqr(dWdy.C))/sqr(1+tau_fsd*C))*Fsd/c_bar/(1-c_bar));
   } else {
     return (0.0);
   } /* endif */

}

/*****************************************************************
 * LES3DFsd_pState::M_xx -- Gradient of surface averaged normal  *
 *****************************************************************/
double LES3DFsd_pState::M_xx(const LES3DFsd_pState &dWdx,
                             const LES3DFsd_pState &dWdy,
                             const LES3DFsd_pState &dWdz,
                             const LES3DFsd_pState &d_dWdx_dx,
                             const LES3DFsd_pState &d_dWdx_dy,
                             const LES3DFsd_pState &d_dWdx_dz) const {

    double Mxx, magnitude_C;
    magnitude_C = sqrt(sqr(dWdx.C)+sqr(dWdy.C)+sqr(dWdz.C));

    if (C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO ) { 
      Mxx = d_dWdx_dx.C/magnitude_C-dWdx.C*(dWdx.C*d_dWdx_dx.C+dWdy.C*d_dWdx_dy.C+dWdz.C*d_dWdx_dz.C)/pow(magnitude_C,3.0);
      return ( Mxx );
   }else{
      return (0.0);
   } /* endif */

}

/******************************************************************
 * LES3DFsd_pState::M_yy -- Gradient of surface averaged normal   *
 ******************************************************************/
double LES3DFsd_pState::M_yy(const LES3DFsd_pState &dWdx,
                             const LES3DFsd_pState &dWdy,
                             const LES3DFsd_pState &dWdz,
                             const LES3DFsd_pState &d_dWdy_dy,
                             const LES3DFsd_pState &d_dWdx_dy,
                             const LES3DFsd_pState &d_dWdy_dz) const {

    double Myy, magnitude_C;
    magnitude_C = sqrt(sqr(dWdx.C)+sqr(dWdy.C)+sqr(dWdz.C));
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO ) { 
      Myy = -d_dWdy_dy.C/magnitude_C+dWdy.C*(dWdx.C*d_dWdx_dy.C+dWdy.C*d_dWdy_dy.C+dWdz.C*d_dWdy_dz.C)/pow(magnitude_C,3.0);
      return ( Myy );
    } else {
      return (ZERO);
    } /* endif */

}

/*******************************************************************
 * LES3DFsd_pState::M_zz -- Gradient of surface averaged normal    *
 *******************************************************************/
double LES3DFsd_pState::M_zz(const LES3DFsd_pState &dWdx,
                             const LES3DFsd_pState &dWdy,
                             const LES3DFsd_pState &dWdz,
                             const LES3DFsd_pState &d_dWdz_dz,
                             const LES3DFsd_pState &d_dWdx_dz,
                             const LES3DFsd_pState &d_dWdy_dz) const {

    double Mzz, magnitude_C;
    magnitude_C = sqrt(sqr(dWdx.C)+sqr(dWdy.C)+sqr(dWdz.C));
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO ) { 
      Mzz = -d_dWdz_dz.C/magnitude_C+dWdz.C*(dWdx.C*d_dWdx_dz.C+dWdy.C*d_dWdy_dz.C+dWdz.C*d_dWdz_dz.C)/pow(magnitude_C,3.0);
      return ( Mzz );
    } else {
      return (ZERO);
    } /* endif */

}

/*******************************************************************
 * LES3DFsd_pState::M_xy -- Gradient of surface averaged normal    *
 *******************************************************************/
double LES3DFsd_pState::M_xy(const LES3DFsd_pState &dWdx,
                             const LES3DFsd_pState &dWdy,
                             const LES3DFsd_pState &dWdz,
                             const LES3DFsd_pState &d_dWdy_dy,
                             const LES3DFsd_pState &d_dWdx_dy,
                             const LES3DFsd_pState &d_dWdy_dz) const {

    double Mxy, magnitude_C;
    magnitude_C = sqrt(sqr(dWdx.C)+sqr(dWdy.C)+sqr(dWdz.C));
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO ) { 
      Mxy = -d_dWdx_dy.C/magnitude_C+dWdx.C*(dWdx.C*d_dWdx_dy.C+dWdy.C*d_dWdy_dy.C+dWdz.C*d_dWdy_dz.C)/pow(magnitude_C,3.0);
      return (Mxy);
    } else {
      return (ZERO);
    } /* endif */

}

/*******************************************************************
 * LES3DFsd_pState::M_xz -- Gradient of surface averaged normal    *
 *******************************************************************/
double LES3DFsd_pState::M_xz(const LES3DFsd_pState &dWdx,
                             const LES3DFsd_pState &dWdy,
                             const LES3DFsd_pState &dWdz,
                             const LES3DFsd_pState &d_dWdz_dz,
                             const LES3DFsd_pState &d_dWdx_dz,
                             const LES3DFsd_pState &d_dWdy_dz) const {

    double Mxz, magnitude_C;
    magnitude_C = sqrt(sqr(dWdx.C)+sqr(dWdy.C)+sqr(dWdz.C));
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO ) { 
      Mxz = -d_dWdx_dz.C/magnitude_C+dWdx.C*(dWdx.C*d_dWdx_dz.C+dWdy.C*d_dWdy_dz.C+dWdz.C*d_dWdz_dz.C)/pow(magnitude_C,3.0);
      return ( Mxz );
   }else{
      return (0.0);
   } /* endif */

}

/*******************************************************************
 * LES3DFsd_pState::M_yz -- Gradient of surface averaged normal    *
 *******************************************************************/
double LES3DFsd_pState::M_yz(const LES3DFsd_pState &dWdx,
                             const LES3DFsd_pState &dWdy,
                             const LES3DFsd_pState &dWdz,
                             const LES3DFsd_pState &d_dWdz_dz,
                             const LES3DFsd_pState &d_dWdx_dz,
                             const LES3DFsd_pState &d_dWdy_dz) const {

    double Myz, magnitude_C;
    magnitude_C = sqrt(sqr(dWdx.C)+sqr(dWdy.C)+sqr(dWdz.C));
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO ) { 
      Myz = -d_dWdy_dz.C/magnitude_C+dWdy.C*(dWdx.C*d_dWdx_dz.C+dWdy.C*d_dWdy_dz.C+dWdz.C*d_dWdz_dz.C)/pow(magnitude_C,3.0);
      return ( Myz );
   }else{
      return (0.0);
   } /* endif */

}

/*******************************************************************
 * LES3DFsd_pState::Resolved_Curvature -- Resolved curvature term  *
 *******************************************************************/
double LES3DFsd_pState::Resolved_Curvature(const LES3DFsd_pState &dWdx,
                                           const LES3DFsd_pState &dWdy,
                                           const LES3DFsd_pState &dWdz,
                                           const LES3DFsd_pState &d_dWdx_dx,
                                           const LES3DFsd_pState &d_dWdy_dy,
                                           const LES3DFsd_pState &d_dWdz_dz,
                                           const LES3DFsd_pState &d_dWdx_dy,
                                           const LES3DFsd_pState &d_dWdx_dz,
                                           const LES3DFsd_pState &d_dWdy_dz) const {

   double tau_fsd, Mxx, Myy, Mzz, 
          resolved_curvature_xx, resolved_curvature_yy, resolved_curvature_zz;

   tau_fsd = HeatRelease_Parameter();
   Mxx = M_xx(dWdx,dWdy,dWdz,d_dWdx_dx,d_dWdx_dy,d_dWdx_dz);
   Myy = M_yy(dWdx,dWdy,dWdz,d_dWdy_dy,d_dWdx_dy,d_dWdy_dz);
   Mzz = M_zz(dWdx,dWdy,dWdz,d_dWdz_dz,d_dWdx_dz,d_dWdy_dz);

   if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C!= ZERO && dWdz.C != ZERO ) {
     resolved_curvature_xx = _laminar_flame_speed*(1.0+tau_fsd*C)*Fsd*rho*Mxx;
     resolved_curvature_yy = _laminar_flame_speed*(1.0+tau_fsd*C)*Fsd*rho*Myy;
     resolved_curvature_zz = _laminar_flame_speed*(1.0+tau_fsd*C)*Fsd*rho*Mzz;
     return ( resolved_curvature_xx + resolved_curvature_yy + resolved_curvature_zz);
   } else {
     return (0.0);
   } /* endif */

}

/***********************************************************************
 * LES3DFsd_pState::Resolved_Propagation -- Resolved propagation term  *
 ***********************************************************************/
double LES3DFsd_pState::Resolved_Propagation(const LES3DFsd_pState &dWdx,
                                             const LES3DFsd_pState &dWdy,
                                             const LES3DFsd_pState &dWdz,
                                             const LES3DFsd_pState &d_dWdx_dx,
                                             const LES3DFsd_pState &d_dWdy_dy,
                                             const LES3DFsd_pState &d_dWdz_dz,
                                             const LES3DFsd_pState &d_dWdx_dy,
                                             const LES3DFsd_pState &d_dWdx_dz,
                                             const LES3DFsd_pState &d_dWdy_dz) const {

   if (C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C!= ZERO && dWdz.C != ZERO ) {
      return (Resolved_Propagation_Curvature(dWdx,dWdy,dWdz)
              -Resolved_Curvature(dWdx,dWdy,dWdz,d_dWdx_dx,d_dWdy_dy,d_dWdz_dz,
                                 d_dWdx_dy,d_dWdx_dz,d_dWdy_dz) );
   } else {
      return (0.0);
   } /* endif */

}

/***************************************************************************************************
 * LES3DFsd_pState::Resolved_Convection_Progvar -- Resolved convection term for progress variable  *
 ***************************************************************************************************/
double LES3DFsd_pState::Resolved_Convection_Progvar(const LES3DFsd_pState &dWdx,
                                                    const LES3DFsd_pState &dWdy,
                                                    const LES3DFsd_pState &dWdz) const {

    double resolved_convection_progvar_x, resolved_convection_progvar_y, 
           resolved_convection_progvar_z;

   if (C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C!= ZERO && dWdz.C != ZERO ) {
     resolved_convection_progvar_x = -(dWdx.rho*v.x*C+rho*dWdx.v.x*C+rho*v.x*dWdx.C);
     resolved_convection_progvar_y = -(dWdy.rho*v.y*C+rho*dWdy.v.y*C+rho*v.y*dWdy.C);
     resolved_convection_progvar_z = -(dWdz.rho*v.z*C+rho*dWdz.v.z*C+rho*v.z*dWdz.C);
     return(resolved_convection_progvar_x+resolved_convection_progvar_y+
             resolved_convection_progvar_z );
   } else {
      return (0.0);
   } /* endif */

}

/**************************************************************************************************
 * LES3DFsd_pState::Resolved_Convection_Fsd -- Resolved convection term for flame surface density.*
 **************************************************************************************************/
double LES3DFsd_pState::Resolved_Convection_Fsd(const LES3DFsd_pState &dWdx,
                                                const LES3DFsd_pState &dWdy,
                                                const LES3DFsd_pState &dWdz) const {

    double resolved_convection_fsd_x, resolved_convection_fsd_y, 
           resolved_convection_fsd_z;

    if (C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C!= ZERO && dWdz.C != ZERO ) {
      resolved_convection_fsd_x = -(dWdx.rho*v.x*Fsd+rho*dWdx.v.x*Fsd+rho*v.x*dWdx.Fsd);
      resolved_convection_fsd_y = -(dWdy.rho*v.y*Fsd+rho*dWdy.v.y*Fsd+rho*v.y*dWdy.Fsd);
      resolved_convection_fsd_z = -(dWdz.rho*v.z*Fsd+rho*dWdz.v.z*Fsd+rho*v.z*dWdz.Fsd);
      return( resolved_convection_fsd_x+resolved_convection_fsd_y+resolved_convection_fsd_z );
   }else{
      return (0.0);
   } /* endif */

}

/****************************************************************************
 * LES3DFsd_pState::NGT_Progvar -- Non-gradient term for progress variable  *
 ****************************************************************************/
double LES3DFsd_pState::NGT_Progvar(const LES3DFsd_pState &dWdx,
                                    const LES3DFsd_pState &dWdy,
                                    const LES3DFsd_pState &dWdz) const {

   double tau_fsd, NGT_progvar_x, NGT_progvar_y, NGT_progvar_z;
   tau_fsd = HeatRelease_Parameter();

   if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C!= ZERO && dWdz.C != ZERO ) {
     NGT_progvar_x = -tau_fsd*_laminar_flame_speed*(rho*(1-2*C)*dWdx.C+C*(1-C)*dWdx.rho);
     NGT_progvar_y = -tau_fsd*_laminar_flame_speed*(rho*(1-2*C)*dWdy.C+C*(1-C)*dWdy.rho);
     NGT_progvar_z = -tau_fsd*_laminar_flame_speed*(rho*(1-2*C)*dWdz.C+C*(1-C)*dWdz.rho);
     return ( NGT_progvar_x+NGT_progvar_y+NGT_progvar_z );
   } else {
     return (0.0);
   } /* endif */

}

/****************************************************************************
 * LES3DFsd_pState::NGT_Fsd -- Non-gradient term for flame surface density  *
 ****************************************************************************/
double LES3DFsd_pState::NGT_Fsd(const LES3DFsd_pState &dWdx,
                                const LES3DFsd_pState &dWdy,
                                const LES3DFsd_pState &dWdz,
                                const LES3DFsd_pState &d_dWdx_dx,
                                const LES3DFsd_pState &d_dWdy_dy,
                                const LES3DFsd_pState &d_dWdz_dz,
                                const LES3DFsd_pState &d_dWdx_dy,
                                const LES3DFsd_pState &d_dWdx_dz,
                                const LES3DFsd_pState &d_dWdy_dz) const {

   double tau_fsd, Mx, My, Mz, Mxx, Myy, Mzz, NGT_fsd_x, NGT_fsd_y, NGT_fsd_z;

   tau_fsd = HeatRelease_Parameter();
   Mx = M_x(dWdx,dWdy,dWdz);
   My = M_y(dWdx,dWdy,dWdz);
   Mz = M_z(dWdx,dWdy,dWdz);
   Mxx = M_xx(dWdx,dWdy,dWdz,d_dWdx_dx,d_dWdx_dy,d_dWdx_dz);
   Myy = M_yy(dWdx,dWdy,dWdz,d_dWdy_dy,d_dWdx_dy,d_dWdy_dz);
   Mzz = M_zz(dWdx,dWdy,dWdz,d_dWdz_dz,d_dWdx_dz,d_dWdy_dz);

   if (C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C!= ZERO && dWdz.C != ZERO ) {
      NGT_fsd_x = -tau_fsd*_laminar_flame_speed*((0.5-C)*
                  (Fsd*Mx*dWdx.rho+rho*Mx*dWdx.Fsd+rho*Fsd*Mxx)-rho*Fsd*Mx*dWdx.C);
      NGT_fsd_y = -tau_fsd*_laminar_flame_speed*((0.5-C)*
                  (Fsd*My*dWdy.rho+rho*My*dWdy.Fsd+rho*Fsd*Myy)-rho*Fsd*My*dWdy.C);
      NGT_fsd_z = -tau_fsd*_laminar_flame_speed*((0.5-C)*
                  (Fsd*Mz*dWdz.rho+rho*Mz*dWdz.Fsd+rho*Fsd*Mzz)-rho*Fsd*Mz*dWdz.C);
      return (NGT_fsd_x+NGT_fsd_y+NGT_fsd_z);
   } else {
      return (0.0);
   } /* endif */ 

}

/***************************************************************************************************
 * LES3DFsd_pState::SFS_Diffusion_Progvar -- Subfilter scale diffusion term for progress variable  *
 ***************************************************************************************************/
double LES3DFsd_pState::SFS_Diffusion_Progvar(const LES3DFsd_pState &dWdx,
                                              const LES3DFsd_pState &dWdy,
                                              const LES3DFsd_pState &dWdz,
                                              const LES3DFsd_pState &d_dWdx_dx,
                                              const LES3DFsd_pState &d_dWdy_dy,
                                              const LES3DFsd_pState &d_dWdz_dz,
                                              const LES3DFsd_pState &d_dWdx_dy,
                                              const LES3DFsd_pState &d_dWdx_dz,
                                              const LES3DFsd_pState &d_dWdy_dz,
                                              const int Flow_Type,
                                              const double &Volume) {

   double grad_eddyviscosity_x, grad_eddyviscosity_y, grad_eddyviscosity_z, 
          sfs_diffusion_progvar_x, sfs_diffusion_progvar_y, sfs_diffusion_progvar_z;
   double Cv = 0.0184, Cs = 0.086, Schmidt_sfs = 1.0;

   if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY) {
      double SS = abs_strain_rate(dWdx,dWdy,dWdz);
      Vector3D grad_rate = grad_abs_strain_rate(dWdx,dWdy,dWdz,d_dWdx_dx,
                                                d_dWdy_dy,d_dWdz_dz,d_dWdx_dy,
                                                d_dWdx_dz,d_dWdy_dz);
      grad_eddyviscosity_x = sqr(Cs*filter_width(Volume))*TWO*grad_rate.x/SS;
      grad_eddyviscosity_y = sqr(Cs*filter_width(Volume))*TWO*grad_rate.y/SS;
      grad_eddyviscosity_z = sqr(Cs*filter_width(Volume))*TWO*grad_rate.z/SS;
   } else if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K) {
      grad_eddyviscosity_x = HALF*Cv*filter_width(Volume)*dWdx.k/sqrt(k);
      grad_eddyviscosity_y = HALF*Cv*filter_width(Volume)*dWdy.k/sqrt(k);
      grad_eddyviscosity_z = HALF*Cv*filter_width(Volume)*dWdz.k/sqrt(k);
   } /* endif */

   if (C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C!= ZERO && dWdz.C != ZERO ) {
      sfs_diffusion_progvar_x = (mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume)*
                                 dWdx.C*dWdx.rho+rho*dWdx.C*grad_eddyviscosity_x+
                                 mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume)*d_dWdx_dx.C)/
                                Schmidt_sfs;
      sfs_diffusion_progvar_y = (mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume)*dWdy.C*
                                 dWdy.rho+rho*dWdy.C*grad_eddyviscosity_y+
                                 mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume)*d_dWdy_dy.C)/
                                Schmidt_sfs;
      sfs_diffusion_progvar_z = (mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume)*
                                 dWdz.C*dWdz.rho+rho*dWdz.C*grad_eddyviscosity_z+
                                 mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume)*d_dWdz_dz.C)/
                                Schmidt_sfs;
      return (sfs_diffusion_progvar_x+sfs_diffusion_progvar_y+sfs_diffusion_progvar_z);
   } else {
      return (0.0);
   } /* endif */

}

/**************************************************************************************************
 * LES3DFsd_pState::SFS_Diffusion_Fsd -- Subfilter scale diffusion term for flame surface density.*
 **************************************************************************************************/
double LES3DFsd_pState::SFS_Diffusion_Fsd(const LES3DFsd_pState &dWdx,
                                          const LES3DFsd_pState &dWdy,
                                          const LES3DFsd_pState &dWdz,
                                          const LES3DFsd_pState &d_dWdx_dx,
                                          const LES3DFsd_pState &d_dWdy_dy,
                                          const LES3DFsd_pState &d_dWdz_dz,
                                          const LES3DFsd_pState &d_dWdx_dy,
                                          const LES3DFsd_pState &d_dWdx_dz,
                                          const LES3DFsd_pState &d_dWdy_dz,
                                          const int Flow_Type,
                                          const double &Volume) {

    double grad_eddyviscosity_x, grad_eddyviscosity_y, grad_eddyviscosity_z, 
           sfs_diffusion_fsd_x, sfs_diffusion_fsd_y, sfs_diffusion_fsd_z;
    double Cv = 0.0184, Cs = 0.086, Schmidt_sfs = 1.0;

    if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY) {
      double SS = abs_strain_rate(dWdx,dWdy,dWdz);
      Vector3D grad_rate = grad_abs_strain_rate(dWdx,dWdy,dWdz,d_dWdx_dx,d_dWdy_dy,d_dWdz_dz,d_dWdx_dy,d_dWdx_dz,d_dWdy_dz);
      grad_eddyviscosity_x = sqr(Cs*filter_width(Volume))*TWO*grad_rate.x/SS;
      grad_eddyviscosity_y = sqr(Cs*filter_width(Volume))*TWO*grad_rate.y/SS;
      grad_eddyviscosity_z = sqr(Cs*filter_width(Volume))*TWO*grad_rate.z/SS;
    }else if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K) {
      grad_eddyviscosity_x = HALF*Cv*filter_width(Volume)*dWdx.k/sqrt(k);
      grad_eddyviscosity_y = HALF*Cv*filter_width(Volume)*dWdy.k/sqrt(k);
      grad_eddyviscosity_z = HALF*Cv*filter_width(Volume)*dWdz.k/sqrt(k);
    } /* endif */
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C!= ZERO && dWdz.C != ZERO ) {
      sfs_diffusion_fsd_x = (mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume)*dWdx.Fsd*dWdx.rho+rho*dWdx.Fsd*grad_eddyviscosity_x+mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume)*d_dWdx_dx.Fsd)/Schmidt_sfs;
      sfs_diffusion_fsd_y = (mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume)*dWdy.Fsd*dWdy.rho+rho*dWdy.Fsd*grad_eddyviscosity_y+mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume)*d_dWdy_dy.Fsd)/Schmidt_sfs;
      sfs_diffusion_fsd_z = (mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume)*dWdx.Fsd*dWdx.rho+rho*dWdx.Fsd*grad_eddyviscosity_x+mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume)*d_dWdx_dx.Fsd)/Schmidt_sfs;
      return ( sfs_diffusion_fsd_x+sfs_diffusion_fsd_y+sfs_diffusion_fsd_z );
    }else{
      return (0.0);
    } /* endif */

}

/*********************************************************************
 * LES3DFsd_pState::Heat_Release_Strain -- Heat release strain term  *
 *********************************************************************/
double LES3DFsd_pState::Heat_Release_Strain(const LES3DFsd_pState &dWdx,
                                            const LES3DFsd_pState &dWdy,
                                            const LES3DFsd_pState &dWdz,
                                            const LES3DFsd_pState &d_dWdx_dx,
                                            const LES3DFsd_pState &d_dWdy_dy,
                                            const LES3DFsd_pState &d_dWdz_dz,
                                            const LES3DFsd_pState &d_dWdx_dy,
                                            const LES3DFsd_pState &d_dWdx_dz,
 	                                    const LES3DFsd_pState &d_dWdy_dz) {

     double tau_fsd, Mxx, Myy, Mzz, heat_release_strain_xx, 
            heat_release_strain_yy, heat_release_strain_zz;
     tau_fsd = HeatRelease_Parameter();

     Mxx = M_xx(dWdx,dWdy,dWdz,d_dWdx_dx,d_dWdx_dy,d_dWdx_dz);
     Myy = M_yy(dWdx,dWdy,dWdz,d_dWdy_dy,d_dWdx_dy,d_dWdy_dz);
     Mzz = M_zz(dWdx,dWdy,dWdz,d_dWdz_dz,d_dWdx_dz,d_dWdy_dz);

     if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C!= ZERO && dWdz.C != ZERO ) {
       heat_release_strain_xx = (0.5-C)*tau_fsd*_laminar_flame_speed*Fsd*rho*Mxx;
       heat_release_strain_yy = (0.5-C)*tau_fsd*_laminar_flame_speed*Fsd*rho*Myy;
       heat_release_strain_zz = (0.5-C)*tau_fsd*_laminar_flame_speed*Fsd*rho*Mzz;
       return ( heat_release_strain_xx + heat_release_strain_yy + heat_release_strain_zz );
    }else{
       return (0.0);
    } /* endif */

}

/**************************************************************************************
 * LES3DFsd_pState::Net_Rate_Change_Progvar -- Net rate change for progress variable  *
 **************************************************************************************/
double LES3DFsd_pState::Net_Rate_Change_Progvar(const LES3DFsd_pState &dWdx,
                                                const LES3DFsd_pState &dWdy,
                                                const LES3DFsd_pState &dWdz,
                                                const LES3DFsd_pState &d_dWdx_dx,
                                                const LES3DFsd_pState &d_dWdy_dy,
                                                const LES3DFsd_pState &d_dWdz_dz,
                                                const LES3DFsd_pState &d_dWdx_dy,
                                                const LES3DFsd_pState &d_dWdx_dz,
						const LES3DFsd_pState &d_dWdy_dz,
                                                const int Flow_Type,
                                                const double &Volume) {

  return(Resolved_Convection_Progvar(dWdx,dWdy,dWdz)+
         SFS_Diffusion_Progvar(dWdx,dWdy,dWdz,d_dWdx_dx,d_dWdy_dy,d_dWdz_dz,
                               d_dWdx_dy,d_dWdx_dz,d_dWdy_dz,Flow_Type,Volume)+
         Reaction_Rate_Fsd(dWdx,dWdy,dWdz) );

}

/**************************************************************************************
 * LES3DFsd_pState::Net_Rate_Change_Fsd -- Net rate change for flame surface density  *
 **************************************************************************************/
double LES3DFsd_pState::Net_Rate_Change_Fsd(const LES3DFsd_pState &dWdx,
                                            const LES3DFsd_pState &dWdy,
                                            const LES3DFsd_pState &dWdz,
                                            const LES3DFsd_pState &d_dWdx_dx,
                                            const LES3DFsd_pState &d_dWdy_dy,
                                            const LES3DFsd_pState &d_dWdz_dz,
                                            const LES3DFsd_pState &d_dWdx_dy,
                                            const LES3DFsd_pState &d_dWdx_dz,
			  		    const LES3DFsd_pState &d_dWdy_dz,
                                            const int Flow_Type,
                                            const double &Volume) {

  return(Resolved_Convection_Fsd(dWdx,dWdy,dWdz)+
         SFS_Diffusion_Fsd(dWdx,dWdy,dWdz,d_dWdx_dx,d_dWdy_dy,d_dWdz_dz,
                           d_dWdx_dy,d_dWdx_dz,d_dWdy_dz,Flow_Type,Volume)
	 +Resolved_Strain(dWdx,dWdy,dWdz)
	 +Resolved_Propagation_Curvature(dWdx,dWdy,dWdz)
	 +SFS_Strain(dWdx,dWdy,dWdz,Flow_Type,Volume)
	 +SFS_Curvature(dWdx,dWdy,dWdz) );

}

/*****************************************************************
 * LES3DFsd_pState::K_equ_sources -- Source term for k-equation  *
 *****************************************************************/
double LES3DFsd_pState::K_equ_sources(const LES3DFsd_pState &dWdx,
                                      const LES3DFsd_pState &dWdy,
                                      const LES3DFsd_pState &dWdz,
                                      const int Flow_Type,
                                      const double &Volume) {

     double production, dissipation, source, mu_t;
     Tensor3D subfilter_stress;
     double Ceps = 0.845;
 
     subfilter_stress = tau_t(dWdx, dWdy, dWdz, Flow_Type, Volume);
  
     if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C!= ZERO && dWdz.C != ZERO ) {
       production = subfilter_stress.xx*dWdx.v.x + 
       subfilter_stress.xy*(dWdy.v.x + dWdx.v.y) +
       subfilter_stress.yy*dWdy.v.y +
       subfilter_stress.xz*(dWdz.v.x + dWdx.v.z) +
       subfilter_stress.yz*(dWdz.v.y + dWdy.v.z) +
       subfilter_stress.zz*dWdz.v.z;

       dissipation = Ceps*rho*pow(k, 1.5)/filter_width(Volume);
       source = production - dissipation;
       return(source);
     }else{
       return (0.0);
    } /* endif */

}

/*****************************************************************
 * LES3DFsd_pState::Sturbulene -- Turbulence model source terms. *
 *****************************************************************/
LES3DFsd_cState LES3DFsd_pState::Sturbchem(LES3DFsd_pState &Wc,
                                           const LES3DFsd_pState &dWdx,
                                           const LES3DFsd_pState &dWdy,
                                           const LES3DFsd_pState &dWdz,
                                           const int Flow_Type,
                                           const double &Volume) {
   
  double resolved_strain, resolved_propagation_curvature, 
         sfs_strain, sfs_curvature;
  LES3DFsd_cState Temp;

  if (Wc.C < 0.999 && Wc.C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO ) {
    // Reaction Rate for Progress Variable Equation --- Source term
      Temp.rhoC = Wc.Reaction_Rate_Fsd(dWdx,dWdy,dWdz);
    // FSD Equation Source Term
      resolved_strain = Wc.Resolved_Strain(dWdx,dWdy,dWdz);

      resolved_propagation_curvature = Wc.Resolved_Propagation_Curvature(dWdx,dWdy,dWdz);

      sfs_strain = Wc.SFS_Strain(dWdx,dWdy,dWdz,Flow_Type,Volume);

      sfs_curvature = Wc.SFS_Curvature(dWdx,dWdy,dWdz);

      Temp.rhoFsd = resolved_strain + resolved_propagation_curvature + 
                    sfs_strain + sfs_curvature;

      if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K) {
         // k-equation Source Term
         Temp.rhok = Wc.K_equ_sources(dWdx,dWdy,dWdz,Flow_Type,Volume);
      } /* endif */
  } /* endif */

  return (Temp);

} 

/*****************************************************************************************
 * LESS3DFsd_cState member functions                                                     *
 *****************************************************************************************/

/***************************************************************************************
 * LES3DFsd_cState::Copy -- Makes a copy of solution state vector.                     *
 ***************************************************************************************/
void LES3DFsd_cState::Copy(const LES3DFsd_cState &U) {
  rho = U.rho;
  rhov = U.rhov; 
  E = U.E;  
  rhoC = U.rhoC; 
  rhoFsd = U.rhoFsd; 
  rhok = U.rhok; 
}

/*********************************************************************************************
 * LES3DFsd_cState::Realizable_Solution_Check -- Check physical validity of solution state.  *
 *********************************************************************************************/
bool LES3DFsd_cState::Realizable_Solution_Check(void) {
  if (rho <= ZERO || es() <= ZERO || rhoC <= ZERO || rhoFsd <= ZERO || rhok <= ZERO ) {    
      cout << "\n " << CFFC_Name() 
           << " ERROR: Conservative solution state has a negative density, energy, progress variable,"
           << " FSD, and/or turbulent kinetic energy.\n";
      return false;
   } else {
      return true;
   } /* endif */
} 

/*************************************************************************
 * LES3DFsd_cState::e -- Return mixture absolute internal energy.        *
 *************************************************************************/
double LES3DFsd_cState::e(void) const {
  return ((E - HALF*rhov.sqr()/rho-rhok)/rho);
}

/**************************************************************************
 * LES3DFsd_cState::es -- Return sensible internal energy.                *
 **************************************************************************/
double LES3DFsd_cState::es(void) const {
  return ((E - HALF*rhov.sqr()/rho-rhok)/rho-HeatofFormation());
}

/**************************************************************************
 * LES3DFsd_cState::h -- Return mixture absolute internal enthalpy.       *
 **************************************************************************/
double LES3DFsd_cState::h(void) const {
  return (e()+p()/rho);
}

double LES3DFsd_cState::h(const double &Temp) const {
  return (Euler3D_ThermallyPerfect_cState::h(Temp));
}

/**************************************************************************
 * LES3DFsd_cState::hs -- Return sensible internal enthalpy.              *
 **************************************************************************/
double LES3DFsd_cState::hs(void) const {
  return (es()+p()/rho);
}

double LES3DFsd_cState::hs(const double &Temp) const {
  return (Euler3D_ThermallyPerfect_cState::hs());
}

/*******************************************************************
 * LES3DFsd_cState::T -- Return mixture temperature.               *
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
double LES3DFsd_cState::T(void) const {

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
         cout << "\nTemperature didn't converge in LES3DFsd_cState::T(void)";
         cout << " with polytopic Tguess " << Tguess << ", or lower than Tmin " 
              << low_temp_range << " using " << T;
      } /* endif */
   } /* endif */

   return T;

} 

/*******************************************************************
 * LES3DFsd_cState::p_t -- Return turbulence modified pressure.     *
 *******************************************************************/
double LES3DFsd_cState::p_t(void) const { 
  return (p() + 2.0*rhok/3.0);
}

/*******************************************************************
 * LES3DFsd_cState::a_t -- Return mixture sound speed (including   *
 *                         turbulent kinetic energy).              *
 *******************************************************************/
double LES3DFsd_cState::a_t(void) const {
   double aa = sqr(a());
   aa += (TWO/THREE)*(rhok/rho)*g();
   return sqrt(aa);
}

/*******************************************************************
 * LES3DFsd_cState::C -- Return progress variable.                 *
 *******************************************************************/
double LES3DFsd_cState::C() const {
   return (rhoC/rho);
}

/*******************************************************************
 * LES3DFsd_cState::Fsd -- Return flame surface density.           *
 *******************************************************************/
double LES3DFsd_cState::Fsd() const {
   return (rhoFsd/rho);
}

/*******************************************************************
 * LES3DFsd_cState::k -- Return turbulent kinetic energy.          *
 *******************************************************************/
double LES3DFsd_cState::k() const {
   return (rhok/rho);
}

/**********************************************************************************************
 * LES3DFsd_cState::premixed_mfrac -- Species mass fractions for one-step reaction mechanism  *
 **********************************************************************************************/
void LES3DFsd_cState::premixed_mfrac(void) {
  double unburnt_oxygen_c, burnt_fuel_c, burnt_oxygen_c, stoich_ratio;
  double c_products, products_ratio;

  LES3DFsd_pState W_temp;
  
  // Check realizability of progress variable.
  Realizable_C_Check();  

  // Determine stoichiometric ratio
  if (W_temp.React.reactset_flag == CH4_1STEP) {
    stoich_ratio = 2.0*specdata[1].Mol_mass()/
                   specdata[0].Mol_mass();     // stoichiometric O2/CH4 mass ratio
  } else if (W_temp.React.reactset_flag == C3H8_1STEP) {
    stoich_ratio = 5.0*specdata[1].Mol_mass()/
                   specdata[0].Mol_mass();     // stoichiometric O2/C3H8 mass ratio
  } else if (W_temp.React.reactset_flag == H2O2_1STEP){
    stoich_ratio = specdata[1].Mol_mass()/
                   2.0/specdata[0].Mol_mass(); // stoichiometric O2/H2 mass ratio
  } /* endif */

  // Calculate fuel mass fraction
  rhospec[0].c = (ONE - rhoC/rho)*rho*_unburnt_fuel_mass_fraction;
  if (rhospec[0].c < MICRO) {
      rhospec[0].c = ZERO;
  } /* endif */

  // Calculate oxygen mass fraction
  // Lean mixture(_fuel_equivalence_ratio < 1) => excessive oxygen
  if (_fuel_equivalence_ratio < ONE) {  
    burnt_fuel_c = ZERO;
    burnt_oxygen_c = (ONE/_fuel_equivalence_ratio - ONE)*_unburnt_fuel_mass_fraction*stoich_ratio;
    if (rhospec[0].c == ZERO) {
      rhospec[1].c = rho*burnt_oxygen_c;      
    } else {
      rhospec[1].c = rhospec[0].c*stoich_ratio + rho*burnt_oxygen_c;
    } /* endif */
    unburnt_oxygen_c = _unburnt_fuel_mass_fraction*stoich_ratio + burnt_oxygen_c;

  // Rich mixture(_fuel_equivalence_ratio > 1) => excessive fuel
  } else if (_fuel_equivalence_ratio > ONE) {  
    burnt_oxygen_c = ZERO;
    burnt_fuel_c = (ONE - ONE/_fuel_equivalence_ratio)*_unburnt_fuel_mass_fraction;
    rhospec[0].c = rhoC*(burnt_fuel_c - _unburnt_fuel_mass_fraction) + rho*_unburnt_fuel_mass_fraction;
    if (rhospec[0].c <= rho*burnt_fuel_c) {
      rhospec[1].c = rho*burnt_oxygen_c;
    } else {
      rhospec[1].c = (rhospec[0].c - rho*burnt_fuel_c)*stoich_ratio;
    } /* endif */
    unburnt_oxygen_c = (_unburnt_fuel_mass_fraction - burnt_fuel_c)*stoich_ratio;

  // Stoichiometric mixture(_fuel_equivalence_ratio = 1)
  } else { 
    burnt_fuel_c = ZERO;
    burnt_oxygen_c = ZERO;
    rhospec[1].c = rhospec[0].c*stoich_ratio; 
    unburnt_oxygen_c = _unburnt_fuel_mass_fraction*stoich_ratio; 
  } /* endif */

  // Determine mass fraction of nitrogen, N2
  if (W_temp.React.reactset_flag == CH4_1STEP ||
      W_temp.React.reactset_flag == C3H8_1STEP) {
     rhospec[4].c = rho*(ONE - _unburnt_fuel_mass_fraction - unburnt_oxygen_c);
  } else if (W_temp.React.reactset_flag == H2O2_1STEP) {
     rhospec[3].c = rho*(ONE - _unburnt_fuel_mass_fraction - unburnt_oxygen_c);
  } /* endif */

  // Determine mass fractions of products
  if (W_temp.React.reactset_flag == CH4_1STEP ||
      W_temp.React.reactset_flag == C3H8_1STEP) {
    // mass fractions of products
    if (rhoC/rho <= MICRO) {
       c_products = ZERO;
    } else {
       c_products = ONE - (rhospec[0].c + rhospec[1].c + rhospec[4].c)/rho;
    } /* endif */
    if (W_temp.React.reactset_flag == CH4_1STEP) {
       products_ratio = specdata[2].Mol_mass()/
                        (specdata[2].Mol_mass()+2.0*specdata[3].Mol_mass());
    } else if (W_temp.React.reactset_flag == C3H8_1STEP) {
       products_ratio = 3.0*specdata[2].Mol_mass()/
                        (3.0*specdata[2].Mol_mass()+4.0*specdata[3].Mol_mass());
    } /* endif */ 
    rhospec[2].c = products_ratio*c_products*rho; // CO2 mass fraction
    rhospec[3].c = rho*c_products-rhospec[2].c;  // H2O mass fraction   
    rhospec[1].c = rho*(ONE - c_products - (rhospec[0].c + rhospec[4].c)/rho) ;  // O2 mass fraction
  } else if (W_temp.React.reactset_flag == H2O2_1STEP) {
    rhospec[2].c = c_products*rho; // H2O mass fraction 
  } /* endif */

  // Check for realizability of each species mass fraction
  for (int i=0; i<ns; i++){
    if (rhospec[i].c < ZERO) rhospec[i].c = ZERO;
  } /* endfor */

  double suma = ZERO;
  for (int i=0; i < ns; i++) {
    suma = suma + rhospec[i].c/rho;
  } /* endfor */

  if (suma > ZERO) {
     for (int i=0; i < ns; i++) {
        rhospec[i].c = rhospec[i].c*(ONE/suma);
     } /* endfor */
  } /* endif */

}

/*******************************************************************
 * LES3DFsd_cState::mu_t -- Return eddy (turbulent) viscosity.     *
 *******************************************************************/
double LES3DFsd_cState::mu_t(const LES3DFsd_pState &dWdx,
 			     const LES3DFsd_pState &dWdy,
			     const LES3DFsd_pState &dWdz,
                             const int Flow_Type, 
                             const double &Volume) {

  double filter = filter_width(Volume);
  if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY) {
    double Cs = 0.18;
    return(rho*Cs*sqr(filter)*abs_strain_rate(dWdx,dWdy,dWdz));
  } else if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K) {
    double Cv = 0.086;
    return(rho*Cv*sqrt(k())*filter);
  }

}

/**************************************************************************
 * LES3DFsd_cState::kappa_t -- Return turbulent thermal conductivity.     *
 **************************************************************************/
double LES3DFsd_cState::kappa_t(const LES3DFsd_pState &dWdx,
 			        const LES3DFsd_pState &dWdy,
			        const LES3DFsd_pState &dWdz,
                                const int Flow_Type, 
                                const double &Volume) {

  return (mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume)*Cp()/Pr_t()); 

}

/********************************************************************************
 * LES3DFsd_cState::Ds_t -- Return species turbulent diffusion coefficient.     *
 ********************************************************************************/
double LES3DFsd_cState::Ds_t(const int i,
                             const LES3DFsd_pState &dWdx,
 			     const LES3DFsd_pState &dWdy,
			     const LES3DFsd_pState &dWdz,
                             const int Flow_Type, 
                             const double &Volume) {

  return (mu_t(dWdx,dWdy,dWdz,Flow_Type,Volume)/(rho*Sc_t())); 

}

double LES3DFsd_cState::Ds_t(const int i,
                             const double &mu_t_temp) {

  return (mu_t_temp/(rho*Sc_t())); 

}

/*******************************************************************
 * LES3DFsd_cState::Pr_t -- Return turbulent Prandtl number.       *
 *******************************************************************/
double LES3DFsd_cState::Pr_t(void) {
   return (0.9); 
}

/*******************************************************************
 * LES3DFsd_cState::Sc_t -- Return turbulent Schmidt number.       *
 *******************************************************************/
double LES3DFsd_cState::Sc_t(void) {
   return (1.0);
}

/*******************************************************************
 * LES3DFsd_cState::Le_t -- Return turbulent Lewis number.         *
 *******************************************************************/
double LES3DFsd_cState::Le_t(void) {
  return (Sc_t()/Pr_t());
}

/*************************************************************************
 * LES3DFsd_cState::filter_width -- LES characteristic filter width      *
 *************************************************************************/
double LES3DFsd_cState::filter_width(const double &Volume) const {
  return (0.0366/60.0);//pow(Volume,1.0/3.0); 
}

/*********************************************************************************
 * LES3DFsd_cState::abs_strain_rate -- Absolute value of strain rate tensor      *
 *********************************************************************************/
double LES3DFsd_cState::abs_strain_rate(const LES3DFsd_pState &dWdx, 
                                        const LES3DFsd_pState &dWdy, 
                                        const LES3DFsd_pState &dWdz) const{

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

/*******************************************************************
 * LES3DFsd_cState::W -- Return primitive solution state vector.   *
 *******************************************************************/
LES3DFsd_pState LES3DFsd_cState::W(void){
   LES3DFsd_pState Temp;
   Temp.rho = rho;
   Temp.v = v();  
   Temp.p = p();
   Temp.C = C();  
   Temp.Fsd = Fsd();
   Temp.k = k();  
   return Temp;
}

LES3DFsd_pState LES3DFsd_cState::W(void) const{
   LES3DFsd_pState Temp;
   Temp.rho = rho;
   Temp.v = v();  
   Temp.p = p();
   Temp.C = C();  
   Temp.Fsd = Fsd();
   Temp.k = k();  
   return Temp;
}

LES3DFsd_pState LES3DFsd_cState::W(const LES3DFsd_cState &U) const {
   LES3DFsd_pState Temp;
   Temp.rho = U.rho;
   Temp.v = U.v();
   Temp.p = U.p();
   Temp.C = U.C();  
   Temp.Fsd = U.Fsd();
   Temp.k = U.k();
   return Temp;
}


/**************************************************************************
 * LES3DFsd_cState::Rotate -- Returns a rotated primitive state aligned   *
 *                            with a local x-axis in the norm_dir.        *
 **************************************************************************/
LES3DFsd_cState LES3DFsd_cState::Rotate(const Vector3D &norm_dir) const {

  // for a 3D unit normal rotated to align with the x-axis
  double Ct = norm_dir.x;  //cos_angle
  double St = sqrt( norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z); //sin_angle
  Vector3D rt(0,norm_dir.z,-norm_dir.y);  //rotation axis  

  return LES3DFsd_cState(rho,
                         rhov.x*Ct - rhov.y*rt.z*St + rhov.z*rt.y*St,
                         rhov.x*rt.z*St + rhov.y*(rt.y*rt.y*(ONE-Ct)+Ct) + 
                         rhov.z*(rt.y*rt.z*(ONE-Ct)),
                         -rhov.x*rt.y*St +  rhov.y*(rt.y*rt.z*(ONE-Ct)) + 
                         rhov.z*(rt.z*rt.z*(ONE-Ct)+Ct),
                         E,
                         rhoC,
                         rhoFsd,
			 rhok);

}

/**************************************************************************
 * LES3DFsd_cState::Rotate -- Returns an un-rotated primitive state       *
 *                            re-alinged from the x-axis of the global    *
 *                            problem.                                    *
 **************************************************************************/
LES3DFsd_cState LES3DFsd_cState::RotateBack(const Vector3D &norm_dir) const {

  // for a 3D unit normal rotated to align with the x-axis
  double Ct = norm_dir.x;  //cos_angle
  double St = sqrt( norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z); //sin_angle
  Vector3D rt(0,norm_dir.z,-norm_dir.y);  //rotation axis  

  return LES3DFsd_cState(rho,
                         rhov.x*Ct + rhov.y*rt.z*St - rhov.z*rt.y*St,
                         -rhov.x*rt.z*St + rhov.y*(rt.y*rt.y*(ONE-Ct)+Ct) + 
                         rhov.z*(rt.y*rt.z*(ONE-Ct)),
                         + rhov.x*rt.y*St +  rhov.y*(rt.y*rt.z*(ONE-Ct)) + 
                         rhov.z*(rt.z*rt.z*(ONE-Ct)+Ct),
                         E,
                         rhoC,
                         rhoFsd,
			 rhok);

}

