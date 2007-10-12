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

  //------------------------------------------------
  // declares
  //------------------------------------------------
  // referece composition (mole fractions)
  double xco(1.0), xh2o(1.0), xco2(1.0), xo2(0.0), xsoot(0.0);

  // reference pressure
  double p(1.0); // 1 atm

  //-----------------------------------------------
  // Setup SNBCK data object
  //-----------------------------------------------

  // setup temporary input parameters object
  SNBCK_Input_Parameters IP;
  
  // use 'full' model with no band lumping
  IP.EvaluationType   = SNBCK_EVAL_ONLINE;
  IP.QuadType         = GAUSS_LEGENDRE;
  IP.QuadPoints       = 12;
  IP.LumpedBands      = 1;
  IP.OptimizedLumping = false;
  IP.OverlapModel     = SNBCK_OVERLAP_UNCORRELATED;

  // perform main setup
  SNBCKdata.Setup(IP, CFFC_PATH);


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
  cout << endl << "Planck Mean Abs. Coeff" << endl;
  cout << setw(12) << "T [K]"
       << setw(12) << "Kp (H2O)" 
       << setw(12) << "Kp (CO)" 
       << setw(12) << "Kp (CO2)"
       << endl;
  for (int n=0; n<Ninterp; n++) {

    // compute temperautre
    Tn[n] = Tmin + double(n)*dT;

    // compute planck mean for each species, kp [atm^-1 m^-1]
    kp_CO[n] = SNBCKdata.PlanckMean( p, Tn[n],  xco,  0.0,  0.0, xo2, xsoot )
      / ( xco * p );
    kp_H2O[n] = SNBCKdata.PlanckMean( p, Tn[n], 0.0, xh2o,  0.0, xo2, xsoot )
      / ( xh2o * p );
    kp_CO2[n] = SNBCKdata.PlanckMean( p, Tn[n], 0.0,  0.0, xco2, xo2, xsoot )
      / ( xco2 * p );

    cout << setw(12) << Tn[n] 
	 << setw(12) << kp_H2O[n] 
	 << setw(12) << kp_CO[n] 
	 << setw(12) << kp_CO2[n] 
	 << endl;
  } // endfor - T

  // Build a natural cubic spline
  spline( Tn,  kp_CO, Ninterp, 1.1E30, 1.1E30,  kp2_CO );
  spline( Tn, kp_CO2, Ninterp, 1.1E30, 1.1E30, kp2_CO2 );
  spline( Tn, kp_H2O, Ninterp, 1.1E30, 1.1E30, kp2_H2O );

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
double PlanckMean :: EvalPlanckMean( const double p,        // pressure [atm]
				     const double T,        // temperature [K]
				     const double xco,      // mole fraction of CO
				     const double xh2o,     // mole fraction of H2O
				     const double xco2 )    // mole fraction of CO2
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
  if ( (T<Tmin || T>Tmax) && err_cnt<5 ) {
    err_cnt++;
    cerr << "\nError in PlanckMean::EvalPlanckMean() : Temperature "
	 << "out of bounds. Tmin=" << Tmin << " < T=" << T 
	 << " < Tmax=" << Tmax << "\n";
    if (err_cnt==5) 
      cerr << "PlanckMean::EvalPlanckMean(): Suppressing error output from now on.\n";
    //exit(-1);
  }

  // interpolate using the cubic spline
  splint( Tn,  kp_CO,  kp2_CO, dT, Ninterp, T,  kkp_CO );
  splint( Tn, kp_CO2, kp2_CO2, dT, Ninterp, T, kkp_CO2 );
  splint( Tn, kp_H2O, kp2_H2O, dT, Ninterp, T, kkp_H2O );

  // compute planck mean using uncorrelated approximation
  // See Ju et al. in J Fluid Mech (1999)
  kp = kkp_CO*xco + kkp_CO2*xco2 + kkp_H2O*xh2o;
  kp *= p;

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
double PlanckMean :: RadSourceOptThin( const double p,        // pressure [atm]
				       const double T,        // temperature [K]
				       const double xco,      // mole fraction of CO
				       const double xh2o,     // mole fraction of H2O
				       const double xco2,     // mole fraction of CO2
				       const double xo2,      // mole fraction of O2
				       const double xsoot )   // volume fraction of soot 
{

  //
  // compute planck mean for gas band radiation
  // See Ju et al. in J Fluid Mech (1999)
  //
  double kp = EvalPlanckMean( p, T, xco, xh2o, xco2 );

  //
  // Add gas band radiation.
  // See Liu, Guo, Smallwood, Gulder, J QSRT 73 (2002) pp. 409-421. 
  //
  double Srad = kp*FOUR*STEFFAN_BOLTZMANN*pow(T,4);

  //
  // Add soot radiation.
  // See Liu, Guo, Smallwood, Gulder, J QSRT 73 (2002) pp. 409-421. 
  //
  static const double C = 3.337E-4; // [W/(m^3 K^5)]
  Srad += C * xsoot * pow(T,5);

  // return value
  return -Srad;

}
