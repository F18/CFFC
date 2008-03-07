/////////////////////////////////////////////////////////////////////
///
/// \file SootState.h
/// 
/// \author Marc R.J. Charest
/// 
/// \brief This header file contains the class member function 
///        definitions for describing the state of soot.
/// 
/// This header file contains the class member function definitions
/// for describing the state of soot.  A simplified reaction mechanism 
/// for the formation, growth, and combustion of soot particles in 
/// laminar nonpremixed flames is used.
///
/// \see K.M. Leung, R.P. Lindstedt, and W.P. Jones, "A simplified
///      reaction mechanism for soot formation in nonpremixed
///      flames," Combustion and Flame 87:289-305 (1991).
///
/////////////////////////////////////////////////////////////////////
#include "Soot2DState.h"
#include "Flame2DState.h"

/**
 * Initialization of static variables.
 */
const int Soot2D_State :: iN;
const int Soot2D_State :: iY;

/////////////////////////////////////////////////////////////////////
/// Member Functions
/////////////////////////////////////////////////////////////////////

/*! *****************************************************************
 * Computation of the soot mean particle diameter..
 *
 * \see K.M. Leung, R.P. Lindstedt, and W.P. Jones, "A simplified
 *      reaction mechanism for soot formation in nonpremixed
 *      flames," Combustion and Flame 87:289-305 (1991).
 *
 * \return Mean diameter in [m].
 ********************************************************************/
double Soot2D_State :: meanDiameter( const Flame2D_pState &W ) const
{

  // soot properties
  const double &N = W.sc(iN);  // number density (particles / kg)
  const double &Y = W.sc(iY);  // mass fraction

  // diameter
  if (Y<TOLER && N<TOLER) return 0.0;
  else return pow( (6.0/PI) * Y/(N*rho), 1.0/3.0);  
}


/*! ****************************************************************
 * Computation of the soot growth source terms. Note that the gas phase
 * mass fractions (excluding soot) sum to unity.  The soot mass fraction
 * is the mass fraction of soot in the whole gas/soot mixture.
 *
 * \see K.M. Leung, R.P. Lindstedt, and W.P. Jones, "A simplified
 *      reaction mechanism for soot formation in nonpremixed
 *      flames," Combustion and Flame 87:289-305 (1991).
 *
 * \param T The gas temperature in [K]
 * \param rho_gas The gas mixture density [kg/m^3]
 * \param iC2H2 The gas phase species index of C2H2
 * \param iO2 The gas phase species index of O2
 * \param iH2 The gas phase species index of H2
 * \param iCO The gas phase species index of CO
 * \param MW_gas Array of molecular weights for gas phase species [kg/kmol]
 * \param y_gas Array of mass fractions for gas phase species [kg/kg gas mix]
 * \param ydot_gas Array of mass fraction production rates for gas phase species [1/s]
 * \param ydot_soot Array of mass fraction production rates for soot [1/s]
 *
 * \return Return error flag, 0 for success, !=0 for error.
 ********************************************************************/
int Soot2D_State :: getRates( const Flame2D_pState &W,
			      Flame2D_State &dUdt,
			      const double &mult ) const
{
  //------------------------------------------------
  // declarations
  //------------------------------------------------
  // Check for necessary species. If one does not exist, error.
  // if ( iO2<0 || iC2H2<0 || iH2<0 || iCO<0 ) return -1;

  // temperature [K]
  const double &T = W.T();

  // soot properties
  const double &N = W.sc(iN);  // number density (particles / kg)
  const double &Y = W.sc(iY);  // mass fraction

  // approximate mixture density with gas phase density
  const double& rho_mix = W.rho(); // [kg/m^3]

  // all gasesous rates are corrected for gas mixture without soot
  const double fact = 1.0 / (1.0-Y);

  // Compute relevant concentrations (note we must correct gas phase 
  // concentrations to include soot mass fraction)
  double X_C2H2 = W.c(iC2H2)*(1.0-Y) * rho_mix / MW_gas[iC2H2];// [kmol/m^3]
  double X_O2   = W.c(iO2)  *(1.0-Y) * rho_mix / MW_gas[iO2];  // [kmol/m^3]
  double X_SOOT = Y                  * rho_mix / MW;           // [kmol/m^3]
 
  //------------------------------------------------
  // Nucleation rate
  //------------------------------------------------
  // >> C2H2 -> 2C(s) + H2
  // rate constant [1/s]
  double k1 = 0.1E5 * exp( -21000.0 / T );
  // rate [kmol/m^3/s]
  double R1 = k1 * X_C2H2;
  // C(s)
  dUdt.rhosc(iY)   +=               mult * 2.0*MW*R1; // [kg/(m^3 s)]
  // H2
  dUdt.rhoc(iH2)   +=   mult * MW_gas[iH2]*R1 * fact; // [kg/(m^3 s)]
  // C2H2
  dUdt.rhoc(iC2H2) -= mult * MW_gas[iC2H2]*R1 * fact; // [kg/(m^3 s)]


  //------------------------------------------------
  // Surface Growth
  //------------------------------------------------
  // Due to the absorbsion of C2H2 on the surface of the particles
  // particle diameter [m]
  // >> C2H2 + nC(s) -> (n+2)C(s) + H2
  double dp;
  if (Y==0 && N==0) dp = 0.0;
  else dp = pow( (6.0/PI) * Y/(N*rho), 1.0/3.0);  
  // Surface area [m^2]
  double S = PI * pow(dp, 2) * rho_mix * N;
  // soot size constant
  double ratio= 6.0*MW/(PI*rho);
  // rate constant [ m^(3/2) / m-soot / s]
  double k2 = 0.6E4 * exp( -12100.0 / T );
  // rate [kmol/m^3/s]
  double R2 = k2 * X_C2H2 * sqrt( PI*pow(ratio, 2.0/3.0) ) *
              pow(X_SOOT, 1.0/3.0) * pow(rho_mix*N, 1.0/6.0);
  // C(s)
  dUdt.rhosc(iY)   +=               mult * 2.0*MW*R2; // [kg/(m^3 s)]
  // H2
  dUdt.rhoc(iH2)   +=   mult * MW_gas[iH2]*R2 * fact; // [kg/(m^3 s)]
  // C2H2
  dUdt.rhoc(iC2H2) -= mult * MW_gas[iC2H2]*R2 * fact; // [kg/(m^3 s)]

  //------------------------------------------------
  // Oxidation
  //------------------------------------------------
  // Oxidation of soot due to O2
  // >> C(s) + 1/2 O2 -> CO
  // rate constant [ m^3 / m^2-soot / s]
  double k3 = 0.1E5 * sqrt(T) * exp( -19680.0 / T );
  // rate [kmol/m^3/s]
  double R3 = k3 * S * X_O2;
  // CO
  dUdt.rhoc(iCO) += mult * 2.0*MW_gas[iCO]*R3 * fact; // [kg/(m^3 s)]
  // C(s)
  dUdt.rhosc(iY) -=                 mult * 2.0*MW*R3; // [kg/(m^3 s)]
  // O2
  dUdt.rhoc(iO2) -=     mult * MW_gas[iO2]*R3 * fact; // [kg/(m^3 s)]

  //------------------------------------------------
  // Agglomeration
  //------------------------------------------------
  // Reduction in number density due to agglomeration of soot 
  // particles.
  // >> nC(s) -> C_n(s)
  // rate constant
  double k4 = 2.0 * Ca * pow(ratio, 1.0/6.0) *
              sqrt(6.0*BOLTZMANN*T/rho);
  // rate [kmol/m^3/s]
  double R4 = k4 * X_SOOT * pow(rho_mix*N, 11.0/6.0);
  // dN/dt
  dUdt.rhosc(iN) = mult * ( 2.0*AVOGADRO*R1/Cmin - R4 ); // [particles/m^3/s]

  // success
  return 0;
} 

/////////////////////////////////////////////////////////////////////
/// Auxiliary Functions
/////////////////////////////////////////////////////////////////////

/**
 * Computation of the soot oxidation rates per unit
 * surface area [kg m-2 s-1].  The model assumes oxidation by
 * O, OH, and O2.  O2 oxidation rates based on the 
 * Nagle–Strickland-Constable model with rate constants for both 
 * O2 and OH from Moss et al.  A constant collision 
 * efficiency of 0.2 was assumed for OH and O.  Rate constant
 * for O was from Bradley et al.
 *
 * \see J. Nagle, R.F. Strickland-Constable, "Oxidation of carbon 
 *      between 1000 and 2000 °C", in Proceedings of the Fifth 
 *      Conference on Carbon, 1962, pp. 154–164.
 *
 * \see J.B. Moss, C.D. Stewart and K.J. Young, 
 *      Combust. Flame 101 (1995), pp. 491–500.
 *
 * \see D. Bradley, G. Dixon-Lewis, S.E. Habik and E.M.J. Mushi, 
 *      Proc. Combust. Inst. 20 (1984), pp. 931–940.
 *
 * \param T The gas temperature in [K]
 * \param rho_gas The gas mixture density [kg/m^3]
 * \param iC2H2 The gas phase species index of C2H2
 * \param iO2 The gas phase species index of O2
 * \param iH2 The gas phase species index of H2
 * \param iCO The gas phase species index of CO
 * \param MW_gas Array of molecular weights for gas phase species [kg/kmol]
 * \param y_gas Array of mass fractions for gas phase species [kg/kg gas mix]
 * \param ydot_gas Array of mass fraction production rates for gas phase species [1/s]
 * \param ydot_soot Array of mass fraction production rates for soot [1/s]
 *
 * \return Total oxidation rate in kg m-2 s-1.
 */
/*
double getRates( const double T, 
		 const double Pa, 
		 const double rho_gas,
		 const int iC2H2,
		 const int iO2,
		 const int iCO,
		 const int iOH,
		 const int iH,
		 const int iO,
		 const double* MW_gas,
		 double* y_gas,
		 double* ydot_gas,
		 double* ydot_soot)
{
  //------------------------------------------------
  // O2 Oxidation
  //------------------------------------------------
  // >> O2 + 0.5C(S) -> CO,
  // rate constants
  double kA = 20.0 * exp(-15098.0/T);
  double kB = 4.46E-03 * exp(-7650.0/T);
  double kT = 1.51E+05 * exp(-48817.0/T);
  double kZ = 21.3 * exp(2063.0/T);
  double chi = 1.0 / (1.0 + kT/(kB*X_O2*Pa));
  // rate [kg/m^2/s]
  double R3 = 120.0 * ( (kA*X_O2*Pa*chi)/(1.0+kZ*X_O2*Pa) +
			kB*X_O2*Pa*(1.0-chi) );
  // CO
  ydot_gas[iCO] += 2.0*R3*As/rho_mix; // [1/s]
  // C(s)
  ydot_soot[0]  -= 2.0*R3*As/rho_mix; // [1/s]
  // O2
  ydot_gas[iO2] -=     R3*As/rho_mix; // [1/s]

  //------------------------------------------------
  // OH Oxidation
  //------------------------------------------------
  // >> OH + C(S) -> CO + H,
  // collision efficiency
  static const double theta_OH = 0.2;
  // rate constant [kg K^(l/2) m^(2) s^(-l) atm^(-1)]
  double k4 = 1.27E+03;
  // rate [kg/m^2/s]
  double R4 = theta_OH * k4 * X_OH * Pa / sqrt(T);
  // CO
  ydot_gas[iCO] += R4*As/rho_mix; // [1/s]
  // H
  ydot_gas[iH ] += R4*As/rho_mix; // [1/s]
  // C(s)
  ydot_soot[0]  -= R4*As/rho_mix; // [1/s]
  // OH
  ydot_gas[iOH] -= R4*As/rho_mix; // [1/s]

  //------------------------------------------------
  // O Oxidation
  //------------------------------------------------
  // >> O + C(S) -> CO.
  // collision efficiency
  static const double theta_O = 0.2;
  // rate constant [kg K^(l/2) m^(2) s^(-l) atm^(-1)]
  double k5 = 665.5;
  // rate [kg/m^2/s]
  double R4 = theta_O * k5 * X_O * Pa / sqrt(T);
  // CO
  ydot_gas[iCO] += R4*As/rho_mix; // [1/s]
  // O
  ydot_gas[iO ] -= R4*As/rho_mix; // [1/s]
  // C(s)
  ydot_soot[0]  -= R4*As/rho_mix; // [1/s]



} 
*/
