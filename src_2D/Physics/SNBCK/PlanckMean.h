/***********************************************************************
 ***********************************************************************
 **  File: PlanckMean.h                                               **
 **                                                                   **
 **  Description:  This header defines the function prototypes,       **
 **                constants, and any class definitions required for  **
 **                the computation of the planck mean absorbsion      **
 **                coefficient.  The SNBCK absorbsion coefficient     **
 **                model is used here. This allows the computation    **
 **                of radiant heat transfer using the optically       **
 **                thin assumption.                                   **
 **                                                                   **
 ** References: Y. Ju, H. Guo, F. Liu, and K. Maruta, "Effects of the **
 **             Lewis number and radiative heat loss on the           **
 **             on the bifurcation and extinction of CH4/O2-N2-He     **
 **             flames."  in J Fluid Mech (1999), vol 379,            **
 **             pp. 165-190.                                          **
 **                                                                   **
 **  Author:  Marc R.J. Charest                                       **
 **                                                                   **
 **  Date:    Dec 17th, 2006                                          **
 ***********************************************************************
 ***********************************************************************/


/*********************************************************************
 ***************************** INCLUDES ******************************
 *********************************************************************/

#ifndef _PLANCK_MEAN_INCLUDED
#define _PLANCK_MEAN_INCLUDED

// Required C++ libraries
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

// Required header files
#include "../../Math/Math.h"
#include "../../Math/SplineFit.h"
#include "../Planck.h"
#include "SNBCK.h"


/*********************************************************************
 *********************** CLASS DEFINITIONS ***************************
 *********************************************************************/

/******************* PLANCK MEAN DEFINTION *********************************
 *                                                                         *
 * Class: PLANCK MEAN                                                      *
 *                                                                         *
 * Description:  The planck mena class is used to store the table of       *
 *               Planck mean absorbsion coefficients computed using        *
 *               the Statistical Narrow Band Correlated-K model for the    *
 *               absorbsion coefficient.  The planck mean absorbsion       *
 *               coefficients are computed for h2o, co2, and co at         *
 *               specified temperature values.  The results are fit with   *
 *               a spline.  This allows the planck mean absorbsion         *
 *               coefficient to be quickly computed during expensive flow  *
 *               calculations.  Significant savings in time can be acheived*
 *               by sing this class and  fitting the average values rather *
 *               than using the SNBCK class and computing the mean         *
 *               absorbsion coefficients using the precalculated method .  *
 *               This is permissable since we do not need to know the      *
 *               values of the abs. coefficient at each band and quadrature*
 *               point.                                                    *
 *                                                                         *
 * Usage:   1 -> Declare object                                            *
 *               >> EM2C SNB;                                              *
 *          2 -> Load the database into memory (given the PATH to the data)*
 *               >> SNB.LoadParams(PATH);                                  *
 *          3a-> Given T, P, and species conc., compute the SNB parameters *
 *               'B' and 'S' for each band, for each gas                   *
 *               >> SNB.ComputeSNB( p, T, xco, xh2o, xco2, xo2 );          *
 *          3b-> Alternatively, to compute gamma at a reference state and  *
 *               generate the SNB model parameter S per unit mole X and    *
 *               per unit atm,                                             *
 *               >> SNB.ComputeRefSNB( p0, T0, xco0, xh2o0, xco20, xo20 ); *
 *          4 -> You can now use S and B as you wish.                      *
 *                                                                         *
 * NOTE: Model parameters are stored based on the following units:         *
 *         length - m                                                      *
 *         press  - atm                                                    *
 *         temp   - K                                                      *
 ***************************************************************************/
class PlanckMean {

  //------------------------------------------------
  // Objects
  //------------------------------------------------

public:
  
private:

  // SNB data object
  EM2C SNB;

  // Interpolation parameters required for cubic spline fitting of
  // the absorbsion coefficient wrt temperature.  Uniform spacing in 
  // temperature is used.
  int Ninterp;       // number of tempeature interpolation points
  double Tmin, Tmax; // temperature limits
  double dT;         // stepsize between temperature interpolation points (uniform spacing)
  double *Tn;        // temperature points
  double *kp_CO,      // kp [atm^-1 m^-1] wrt temp for interpolation
         *kp_CO2,
         *kp_H2O;
  double *kp2_CO,     // second derivative of k wrt temp for interpolation
         *kp2_CO2,
         *kp2_H2O;

  //------------------------------------------------
  // Member Functions
  //------------------------------------------------

public:
  
  //
  // Consutructor 1 - components uninitialized
  //
  PlanckMean( const char *CFFC_PATH ) :  // Current path
    Tn(NULL), kp_CO(NULL), kp_CO2(NULL), kp_H2O(NULL), 
    kp2_CO(NULL), kp2_CO2(NULL), kp2_H2O(NULL), dT(ZERO),
    Ninterp(EM2C::Npoints), Tmin(ZERO), Tmax(ZERO)
  { 
    // allocate
    Allocate(); 

    // Get SNB model parameters for gas mixture (units: cm, atm, K),
    // compute curve fits
    Setup(CFFC_PATH);
  }


  //
  // Consutructor 2 - components uninitialized
  //
  PlanckMean( const char *CFFC_PATH,  // Current path
	      const int N ) :         // number of interpolation points
    Tn(NULL), kp_CO(NULL), kp_CO2(NULL), kp_H2O(NULL), 
    kp2_CO(NULL), kp2_CO2(NULL), kp2_H2O(NULL), dT(ZERO),
    Tmin(ZERO), Tmax(ZERO) 
  { 
    // Check to make sure there are enough points.
    // If there aren't, default to number EM2C data points.
    if ( N<3 ) {
      cerr << "\nPlanckMean::PlanckMean() - Warning, number of interpolation"
	   << " points N=" << N << " < 3, using " << EM2C::Npoints
	   << endl;
      Ninterp = EM2C::Npoints;
    } else {
      Ninterp = N;
    } // endif

    // allocate
    Allocate(); 

    // Get SNB model parameters for gas mixture (units: cm, atm, K),
    // compute curve fits
    Setup(CFFC_PATH);
  }


  //
  // destructor
  //
  ~PlanckMean() {  Deallocate(); }

  // calculate the planck mean absorbsion coefficient [m^-1]
  double EvalPlanckMean( const double p,        // pressure [atm]
			 const double T,        // temperature [K]
			 const double xco,      // mole fraction of CO
			 const double xh2o,     // mole fraction of H2O
			 const double xco2 );   // mole fraction of CO2
  
  // calculate the volumetric radiation source based on the 
  // optically thin approximation in W/m^3
  double RadSourceOptThin( const double p,        // pressure [atm]
			   const double T,        // temperature [K]
			   const double xco,      // mole fraction of CO
			   const double xh2o,     // mole fraction of H2O
			   const double xco2,     // mole fraction of CO2
			   const double xo2,      // mole fraction of O2
			   const double xsoot );  // volume fraction of soot 
    

private:

  // allocate and deallocate the arrays
  void Allocate();
  void Deallocate();

  // Get SNB model parameters for gas mixture (units: cm, atm, K)
  // and compute curve fits
  void Setup(const char *CFFC_PATH);


};


#endif //_PLANCK_MEAN_INCLUDED
