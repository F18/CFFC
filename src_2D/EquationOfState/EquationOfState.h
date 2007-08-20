/**********************************************************************
 * EquationOfState.h:      Header file defining curve-fits            *
 * of Srinivasan et al. (1987) for a high temperature gas             *
 * equation of state                                                  *
 **********************************************************************/

#ifndef _EQUATIONOFSTATE_INCLUDED
#define _EQUATIONOFSTATE_INCLUDED

#include <cstdlib>

#define F(x) x ## _

//Datatype correlations
// C++      ref     Fortran
typedef int            Fortran_Integer;
typedef double         Fortran_Real;

extern "C"
{

  /* TGAS subroutines for thermodynamic properties. */

  // TGAS1 returns equilibrium pressure, sound speed, and temperature 
  // given the specific internal energy and density
  void F(tgas1)(Fortran_Real *E, 
                Fortran_Real *R, 
                Fortran_Real *P, 
                Fortran_Real *A, 
                Fortran_Real *T, 
                Fortran_Integer *MFLAG, 
                Fortran_Integer *IERROR);

  // values for MFLAG
#define TGAS1_RETURN_PRESSURE  0
#define TGAS1_RETURN_PRESSURE_AND_SOUND_SPEED  1
#define TGAS1_RETURN_PRESSURE_AND_TEMPERATURE  2
#define TGAS1_RETURN_PRESSURE_TEMPERATURE_AND_SOUND_SPEED  3

  // TGAS2 returns equilibrium entropy
  void F(tgas2)(Fortran_Real *E,
                Fortran_Real *R, 
                Fortran_Real *S, 
                Fortran_Integer *IERROR);
  
  // TGAS3 returns equilibrium temperature
  void F(tgas3)(Fortran_Real *P, 
                Fortran_Real *RHO, 
                Fortran_Real *T, 
                Fortran_Integer *IERROR);
  
  // TGAS4 returns equilibrium specific enthalpy
  void F(tgas4)(Fortran_Real *P, 
                Fortran_Real *RHO, 
                Fortran_Real *H, 
                Fortran_Integer *IERROR);
  
  // TGAS5 returns density, specific internal energy, and sound speed
  void F(tgas5)(Fortran_Real *PRE, 
                Fortran_Real *TEMP, 
                Fortran_Real *DEN, 
                Fortran_Real *EINT, 
                Fortran_Real *SSP, 
                Fortran_Integer *IERROR);
  
  // TGAS6 returns specific heats at constant pressure and volume and specific heat ratio
  void F(tgas6)(Fortran_Real *PRE, 
                Fortran_Real *TEMP, 
                Fortran_Real *CP, 
                Fortran_Real *CV, 
                Fortran_Real *GAM, 
                Fortran_Integer *IERROR);

  /* UGAS subroutines for transport quantities. */

  void F(ugas1)(Fortran_Real *T, 
                Fortran_Real *RHO, 
                Fortran_Real *MU, 
                Fortran_Integer *IERROR);

  void F(ugas2)(Fortran_Real *T, 
                Fortran_Real *RHO, 
                Fortran_Real *PR, 
                Fortran_Integer *IERROR);

  void F(ugas4)(Fortran_Real *E, 
                Fortran_Real *RHO, 
                Fortran_Real *K, 
                Fortran_Integer *IERROR);

}

/**********************************************************************
 * Tgas_p -- Pressure.                                                *
 **********************************************************************/
inline double Tgas_p(double inte, double rho) {
  int MFLAG = TGAS1_RETURN_PRESSURE, IERROR;
  double ss, p, temp;
  F(tgas1)(&inte, &rho, &p, &ss, &temp, &MFLAG, &IERROR);
  assert (IERROR==0);
  return p; 
}

inline double Tgas_a_from_e_rho(double inte, double rho) {
  int MFLAG = TGAS1_RETURN_PRESSURE_AND_SOUND_SPEED, IERROR;
  double ss, p, temp;
  F(tgas1)(&inte, &rho, &p, &ss, &temp, &MFLAG, &IERROR);
  assert (IERROR==0);
  return ss;
}

inline double Tgas_T_from_e_rho(double inte, double rho) {
  int MFLAG = TGAS1_RETURN_PRESSURE_AND_TEMPERATURE, IERROR;
  double ss, p, temp;
  F(tgas1)(&inte, &rho, &p, &ss, &temp, &MFLAG, &IERROR);
  assert (IERROR==0);
  return temp;
}

/**********************************************************************
 * Tgas_temp -- Temperature.                                          *
 **********************************************************************/
inline double Tgas_temp(double p, double rho) {
  int IERROR;
  double temp;
  F(tgas3)(&p, &rho, &temp, &IERROR);
  assert (IERROR==0);
  return temp;
}

/**********************************************************************
 * Tgas_e -- Specific internal energy.                                *
 **********************************************************************/
inline double Tgas_e(double p, double temp) {
  int IERROR;
  double inte, rho, ss;
  F(tgas5)(&p, &temp, &rho, &inte, &ss, &IERROR);
  assert (IERROR==0);
  return inte;
}

/**********************************************************************
 * Tgas_rho -- Density.                                               *
 **********************************************************************/
inline double Tgas_rho(double p, double temp) {
  int IERROR;
  double inte, rho, ss;
  F(tgas5)(&p, &temp, &rho, &inte, &ss, &IERROR);
  assert (IERROR==0);
  return rho;
}

/**********************************************************************
 * Tgas_a -- Sound speed.                                             *
 **********************************************************************/
inline double Tgas_a(double p, double temp) {
  int IERROR;
  double inte, rho, ss;
  F(tgas5)(&p, &temp, &rho, &inte, &ss, &IERROR);
  assert (IERROR==0);
  return ss;
}

/**********************************************************************
 * Tgas_s -- Entropy.                                                 *
 **********************************************************************/
inline double Tgas_s(double inte, double rho) {
  int IERROR;
  double s;
  F(tgas2)(&inte, &rho, &s, &IERROR);
  assert (IERROR==0);
  return s;
}

/**********************************************************************
 * Tgas_h -- Specific Enthalpy.                                       *
 **********************************************************************/
inline double Tgas_h(double p, double rho) {
  int IERROR;
  double h;
  F(tgas4)(&p, &rho, &h, &IERROR);
  return h;
}

/**********************************************************************
 * Tgas_cp -- Specific heat at constant pressure.                     *
 **********************************************************************/
inline double Tgas_cp(double p, double temp) {
  int IERROR;
  double cp, cv, g;
  F(tgas6)(&p, &temp, &cp, &cv, &g, &IERROR);
  assert (IERROR==0);
  return cp;
}

/**********************************************************************
 * Tgas_cv -- Specific heat at constant volume                        *
 **********************************************************************/
inline double Tgas_cv(double p, double temp) {
  int IERROR;
  double cp, cv, g;
  F(tgas6)(&p, &temp, &cp, &cv, &g, &IERROR); 
  assert (IERROR==0);
  return cv;
}

/**********************************************************************
 * Tgas_gamma -- Specific heat ratio.                                 *
 **********************************************************************/
inline double Tgas_gamma(double p, double temp) {
  int IERROR;
  double cp, cv, g;
  F(tgas6)(&p, &temp, &cp, &cv, &g, &IERROR); 
  assert (IERROR==0);
  return g;
}

/**********************************************************************
 * Tgas_mu -- Dynamic viscosity.                                      *
 **********************************************************************/
inline double Tgas_mu(double temp, double rho) {
  int IERROR;
  double mu;
  F(ugas1)(&temp, &rho, &mu, &IERROR);
  assert (IERROR==0);
  return mu;
}

/**********************************************************************
 * Tgas_Pr -- Prandtl number.                                         *
 **********************************************************************/
inline double Tgas_Pr(double temp, double rho) {
  int IERROR;
  double pr;
  F(ugas2)(&temp, &rho, &pr, &IERROR);
  assert (IERROR==0);
  return pr;
}

/**********************************************************************
 * Tgas_kappa -- Coefficient of thermal conductivity.                 *
 **********************************************************************/
inline double Tgas_kappa( double inte, double rho) {
  int IERROR;
  double kappa;
  F(ugas4)(&inte, &rho, &kappa, &IERROR);
  assert (IERROR==0);
  return kappa;
}

#endif // _EQUATIONOFSTATE_INCLUDED
