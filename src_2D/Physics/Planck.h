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

/* Required CFFC header files */
#include "GasConstants.h"


/********************************************************
 * Blackbody emmissive power at specified wavenumber.   *
 *                                                      *
 * wavenumber wn in [cm^-1]                             *
 * temperature T in [K]                                 *
 ********************************************************/
inline double BlackBody(const double T, 
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
 * Blackbody intensity                                  *
 ********************************************************/
inline double BlackBody(const double T) 
{ return STEFFAN_BOLTZMANN*pow(T,4.0)/PI; }


#endif //end _PLANCK_INCLUDED 
