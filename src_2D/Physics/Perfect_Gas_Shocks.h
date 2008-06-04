/*! \file Perfect_Gas_Shocks.h
  \brief Normal-shock jump conditions for perfect gases. */

#ifndef _PERFECT_GAS_SHOCKS_INCLUDED
#define _PERFECT_GAS_SHOCKS_INCLUDED

/******************************************************//**
 * Function: normal_shock_pressure_ratio
 *
 * Returns pressure ratio accross a normal shock
 * (downstream / upstream),
 * given a mach number (M) and gamma (G)
 *
 ********************************************************/
inline double normal_shock_pressure_ratio(double M, double G) {
  return 1.0/(G+1.0)*(2.0*G*M*M-(G-1.0));
}

/******************************************************//**
 * Function: normal_shock_density_ratio
 *
 * Returns density ratio accross a normal shock
 * (downstream / upstream),
 * given a mach number (M) and gamma (G)
 *
 ********************************************************/
inline double normal_shock_density_ratio(double M, double G) {
  return ((G+1.0)*M*M) / ((G-1.0)*M*M+2.0);
}

/******************************************************//**
 * Function: normal_shock_velocity_ratio
 *
 * Returns velocity ratio accross a normal shock
 * (downstream / upstream),
 * given a mach number (M) and gamma (G)
 *
 ********************************************************/
inline double normal_shock_velocity_ratio(double M, double G) {
  return ((G-1.0)*M*M+2.0) / ((G+1.0)*M*M);
}

/******************************************************//**
 * Function: normal_shock_temperature_ratio
 *
 * Returns temperature ratio accross a normal shock
 * (downstream / upstream),
 * given a mach number (M) and gamma (G)
 *
 ********************************************************/
inline double normal_shock_temperature_ratio(double M, double G) {
  return (2.0 +(G-1.0)*M*M) * (2.0*G*M*M-(G-1.0)) / ((G+1.0)*(G+1.0)*M*M);
}

#endif //_PERFECT_GAS_SHOCKS_INCLUDED
