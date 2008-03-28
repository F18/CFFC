/*******************************************************************
  File: PlanckMean.cpp

  Description:  This file defines the PlanckMean class member 
                functions for the computation of the planck mean 
                absorbsion coefficient using the SNBCK.

  Author:  Marc R.J. Charest

  Date:    Dec 25th, 2006
*******************************************************************/
#include "PlanckMean.h"


/*********************************************************************
 ****************** PLANCK MEAN CLASS MEMBER FUNCTIONS ***************
 *********************************************************************/

/*********************************************************************
 * NOTE: Model parameters are stored based on the following units:   *
 *         length - cm                                               *
 *         press  - atm                                              *
 *         temp   - K                                                *
 *********************************************************************/


/*********************************************************************
 * Array allocators/deallocators                                     *
 *********************************************************************/
void PlanckMean :: Allocate() {

  // deallocate just in case
  Deallocate();

  // intialize
  if (Ninterp>0) {
    Tn      = new double[Ninterp];
    kp_CO   = new double[Ninterp];
    kp_CO2  = new double[Ninterp];
    kp_H2O  = new double[Ninterp];
    kp2_CO  = new double[Ninterp];
    kp2_CO2 = new double[Ninterp];
    kp2_H2O = new double[Ninterp];
  } /* endif */

}

void PlanckMean :: Deallocate() {
  if (kp_CO != NULL) { 
    delete[]   kp_CO;  kp_CO  =NULL;
    delete[]  kp_CO2;  kp_CO2 =NULL;
    delete[]  kp_H2O;  kp_H2O =NULL;
    delete[]  kp2_CO;  kp2_CO =NULL;
    delete[] kp2_CO2;  kp2_CO2=NULL;
    delete[] kp2_H2O;  kp2_H2O=NULL;
  }/* endif */

  if (Tn != NULL)  delete[] Tn;  Tn = NULL;
}


/*********************************************************************
 * SNBCK :: Setup                                                    *
 *                                                                   *
 * Load the dataset EM2C of Soufiani and Taine (1997) for CO, CO2,   *
 * and H2O which is used to compute the SNB model parameters.        *
 * Additionally, setup the statisctical narrow band model and        *
 * generate spline fits for the planck mean absorbsion coefficient.  *
 *********************************************************************/
void PlanckMean :: Setup( const char *CFFC_PATH ) // Current path
{

  //-----------------------------------------------
  // Setup SNB data object
  //-----------------------------------------------

  // get SNB model parameters for gas mixture (units: cm, atm, K)
  SNB.LoadParams(CFFC_PATH);


  //------------------------------------------------
  // Initialize
  //------------------------------------------------

  // for uniform spacing
  Tmin = EM2C::Tmin;
  Tmax = EM2C::Tmax;
  dT = (Tmax-Tmin)/double(Ninterp-1);

  //
  // loop over each temperature point
  //
  for (int n=0; n<Ninterp; n++) {

    // compute temperautre
    Tn[n] = Tmin + double(n)*dT;

    // compute planck mean for each species, kp [atm^-1 m^-1]
    SNB.PlanckMean( Tn[n], kp_CO[n], kp_H2O[n], kp_CO2[n] );

  } // endfor - T

  // Build a natural cubic spline
#ifdef PLANCKMEAN_USE_SPLINES
  spline( Tn,  kp_CO, Ninterp, 1.1E30, 1.1E30,  kp2_CO );
  spline( Tn, kp_CO2, Ninterp, 1.1E30, 1.1E30, kp2_CO2 );
  spline( Tn, kp_H2O, Ninterp, 1.1E30, 1.1E30, kp2_H2O );

  //Polynomial
#else //PLANCKMEAN_USE_POLY
  polcof( Tn,  kp_CO, Ninterp,  kp2_CO );
  polcof( Tn, kp_CO2, Ninterp, kp2_CO2 );
  polcof( Tn, kp_H2O, Ninterp, kp2_H2O );
#endif

}


/*********************************************************************
 * SNBCK :: EvalPlanckMean                                           *
 *                                                                   *
 * Calculate the planck mean absorbsion coefficient using the        *
 * optically thin limit. Thus, self-reabsorption of radiation of hot *
 * burned gas is neglected. See:                                     *
 *   Y. Ju, H. Guo, F. Liu, and K. Maruta, J Fluid Mech (1999), vol. *
 *   379, pp. 165-190.                                               *
 *                                                                   *
 * Assuming the ambient is cold, we have:                            *
 *  qr = \int_{-1}^{-1} \int_{0}^{\infty} k_v Ib_v dv d\mu           *
 *     = 4 K_p \sigma T^4                                            *
 *                                                                   *
 * Here, the spline fits are used to interpolate for the value of    *
 * the planck mean at the specified temperature.                     *
 *                                                                   *
 * Placnk mean in [m^-1]                                             *
 *********************************************************************/
double PlanckMean :: EvalPlanckMean( const RadiatingGas &gas )
{
  // declares
  double kp, kkp_CO, kkp_CO2, kkp_H2O;
  static int err_cnt(0);


  // check to make sure everything allocated
  if (kp_CO == NULL) {
    cerr << "\nError in PlanckMean::EvalPlanckMean() : Need to "
	 << "allocate first.\n";
    exit(-1);
  }

  // check to make sure T within bounds
  if ( (gas.T<Tmin || gas.T>Tmax) && err_cnt<5 ) {
    err_cnt++;
    cerr << "\nError in PlanckMean::EvalPlanckMean() : Temperature "
	 << "out of bounds. Tmin=" << Tmin << " < T=" << gas.T 
	 << " < Tmax=" << Tmax << "\n";
    if (err_cnt==5) 
      cerr << "PlanckMean::EvalPlanckMean(): Suppressing error output from now on.\n";
    //exit(-1);
  }

  // interpolate using the cubic spline
#ifdef PLANCKMEAN_USE_SPLINES
  splint( Tn,  kp_CO,  kp2_CO, dT, Ninterp, gas.T,  kkp_CO );
  splint( Tn, kp_CO2, kp2_CO2, dT, Ninterp, gas.T, kkp_CO2 );
  splint( Tn, kp_H2O, kp2_H2O, dT, Ninterp, gas.T, kkp_H2O );

  //Polynomial
#else //PLANCKMEAN_USE_POLY
  polval(  kp2_CO, Ninterp, gas.T,  kkp_CO );
  polval( kp2_CO2, Ninterp, gas.T, kkp_CO2 );
  polval( kp2_H2O, Ninterp, gas.T, kkp_H2O );
#endif

  // compute planck mean using uncorrelated approximation
  // See Ju et al. in J Fluid Mech (1999)
  kp = kkp_CO*gas.xco + kkp_CO2*gas.xco2 + kkp_H2O*gas.xh2o;
  kp *= gas.p;

  // return value
  return kp;

}


/*********************************************************************
 * SNBCK :: RadSourceOptThin                                         *
 *                                                                   *
 * Calculate the divergence of the radiative heat flux in the        *
 * optically thin limit. Thus, self-reabsorption of radiation of hot *
 * burned gas is neglected. See:                                     *
 *   Y. Ju, H. Guo, F. Liu, and K. Maruta, J Fluid Mech (1999), vol. *
 *   379, pp. 165-190.                                               *
 *                                                                   *
 * Assuming the ambient is cold, we have:                            *
 *  qr = \int_{-1}^{-1} \int_{0}^{\infty} k_v Ib_v dv d\mu           *
 *     = 4 K_p \sigma T^4                                            *
 *                                                                   *
 * Radiation source in [W/m^3]                                       *
 *********************************************************************/
double PlanckMean :: RadSourceOptThin( const RadiatingGas &gas )
{

  //
  // compute planck mean for gas band radiation
  // See Ju et al. in J Fluid Mech (1999)
  //
  double kp = EvalPlanckMean( gas );

  //
  // Add gas band radiation.
  // See Liu, Guo, Smallwood, Gulder, J QSRT 73 (2002) pp. 409-421. 
  //
  double Srad = kp*FOUR*STEFFAN_BOLTZMANN*pow(gas.T,4);

  //
  // Add soot radiation.
  // See Liu, Guo, Smallwood, Gulder, J QSRT 73 (2002) pp. 409-421. 
  //
  static const double C = 3.337E-4; // [W/(m^3 K^5)]
  Srad += C * gas.fsoot * pow(gas.T,5);

  // return value
  return -Srad;

}
