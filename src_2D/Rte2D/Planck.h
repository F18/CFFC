/**********************************************************************
 **********************************************************************
 **                                                                  **
 ** File: Planck.h                                                   **
 **                                                                  **
 ** Description: This file defines some necessary functions for      **
 **              the description and evaluation of Plancks law.      **
 **                                                                  **
 ** Author: Marc "T-Bone" Charest                                    **
 **                                                                  **
 ** Revision:  Date        Initials   Change                         **
 **            12/08/2007  MRC        Original creation              **
 **                                                                  **
 **********************************************************************
 **********************************************************************/
#ifndef _PLANCK_INCLUDED
#define _PLANCK_INCLUDED 

/* Required CFDkit header files */
#ifndef _GAS_CONSTANTS_INCLUDED
#include "../Physics/GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED


/********************************************************
 * Blackbody emmissive power at specified wavenumber.   *
 *                                                      *
 * wavenumber wn in [cm^-1]                             *
 * temperature T in [K]                                 *
 ********************************************************/
inline double Planck(const double T, 
	      const double wn ) { 

  // Plancks constants
  static const double C1 = 1.1909E-08; // [W/(m2 ster cm-4)]
  static const double C2 = 1.4388E+00; // [cm*K]

  // check for div by zero
  if (T<NANO) return ZERO;

  // X = eta/T (cm^-1/K)
  // V = C_2*eta/T
  // EN = E_beta/nT^3
  double V=C2*wn/T;
  double EX=exp(V);
  double EN=C1*wn*wn*wn/(EX-1.0); // W/m2cm-1K3
  return EN;                      //W/m2cm-1
}




/********************************************************
 * This subroutine calculates the fractional blackbody	*
 * emissive power f(n*lambda*T), where n*lambda*T in    *
 * (micro-m*K) and refractive index, n, approx 1.       *
 * See:                                                 *
 *   M.F. Modest, "Radiative Heat Transfer," 2nd ed,    *
 *   New York: Academic Press. 2003.                    *
 ********************************************************/
inline double Planck(const double lambdaT) { // micro-m K

  // constants
  double C2 = 1.4388E+04; // [micro-m K]
  double CC = FIFTEEN/pow(PI,4);
  double EPS = 1E-16;


  // V = C_2/lambdaT = C_2*eta/T
  double V  = C2/(lambdaT);
  double EX = exp(V);

  // Evaluation of f(n*lambda*T) in terms of an infinite series
  // func = f(n*lambda*T)
  double func = ZERO;
  double EM = ONE;
  int M = 0;
  bool converged = false;
  double VM, BM;
  while( !converged ) {
    M++;
    VM=M*V;
    BM=(SIX+VM*(SIX+VM*(THREE+VM)))/pow(double(M),4);
    EM=EM/EX;
    func = func + BM*EM;
    if((pow(VM,3)*EM)<EPS) converged = true;
  }
  func *= CC;

  return (func);

}

/********************************************************
 * Blackbody intensity                                  *
 ********************************************************/
inline double Ib(const double T) { return STEFFAN_BOLTZMANN*pow(T,4.0)/PI; }

/********************************************************
 * Blackbody spectral intensity at specified wavenumber.*
 * Radiative energy flow / time / area normal to rays / *
 * solid angle / wavelength [cm].                       *
 *                                                      *
 * wavenumber wn in [cm^-1]                             *
 * temperature T in [K]                                 *
 ********************************************************/
inline double Ib_v(const double T, 
	    const double wn ) { 
  return Planck( T, wn ) * PI;

}


#endif //end _PLANCK_INCLUDED 
