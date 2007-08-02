/**********************************************************************
 * SolidConstants.h:  Header file defining various particle and       *
 *                    propellant constants.                           *
 **********************************************************************/

#ifndef _SOLID_CONSTANTS_INCLUDED
#define _SOLID_CONSTANTS_INCLUDED

// Include the gas constants header file.

#ifndef _GAS_CONSTANTS_INCLUDED
#include "GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED

/**********************************************************************
 * Definition of particle constants.                                  *
 **********************************************************************/

// Constants for AP-HPTB.
#define CM_APHTPB           900.0               // particle specific heat
#define DP_APHTPB           MICRO*TEN           // particle diameter
#define RHOP_APHTPB         2700.0              // particle solid density
#define MP_APHTPB           (PI/SIX)*RHOP_APHTPB*cube(DP_APHTPB)
#define QE_APHTPB           ZERO                // particle charge

// Constants for glass beads.
#define CM_GLASS_BEADS      840.0
#define DP_GLASS_BEADS      0.000010
#define RHOP_GLASS_BEADS    2400.0
#define MP_GLASS_BEADS      (PI/SIX)*RHOP_GLASS_BEADS*cube(DP_GLASS_BEADS)
#define QE_GLASS_BEADS      ZERO                // particle charge

// Constants for the tiny SCIEX particle.
#define CM_SCIEX_TINY       720.0               // particle specific heat
#define DP_SCIEX_TINY       30.0e-9             // particle diameter
#define RHOP_SCIEX_TINY     1000.0              // particle solid density
#define MP_SCIEX_TINY       (PI/SIX)*RHOP_SCIEX_TINY*cube(DP_SCIEX_TINY)
#define QE_SCIEX_TINY       208.0*E_CHARGE      // particle charge

// Constants for the small SCIEX particle.
#define CM_SCIEX_SMALL      720.0               // particle specific heat
#define DP_SCIEX_SMALL      300.0e-9            // particle diameter
#define RHOP_SCIEX_SMALL    1000.0              // particle solid density
#define MP_SCIEX_SMALL      (PI/SIX)*RHOP_SCIEX_SMALL*cube(DP_SCIEX_SMALL)
#define QE_SCIEX_SMALL      6587.0*E_CHARGE     // particle charge

// Constants for the medium SCIEX particle.
#define CM_SCIEX_MEDIUM     720.0               // particle specific heat
#define DP_SCIEX_MEDIUM     1.0e-6              // particle diameter
#define RHOP_SCIEX_MEDIUM   1000.0              // particle solid density
#define MP_SCIEX_MEDIUM     (PI/SIX)*RHOP_SCIEX_MEDIUM*cube(DP_SCIEX_MEDIUM)
#define QE_SCIEX_MEDIUM     40090.0*E_CHARGE    // particle charge

// Constants for the large SCIEX particle.
#define CM_SCIEX_LARGE      720.0               // particle specific heat
#define DP_SCIEX_LARGE      45.0e-6             // particle diameter
#define RHOP_SCIEX_LARGE    1000.0              // particle solid density
#define MP_SCIEX_LARGE      (PI/SIX)*RHOP_SCIEX_LARGE*cube(DP_SCIEX_LARGE)
#define QE_SCIEX_LARGE      (1.210e7)*E_CHARGE  // particle charge

// Constants for HEAVY_APHTPB.
#define CM_HEAVY_APHTPB     900.0               // particle specific heat
#define DP_HEAVY_APHTPB     MICRO*TEN           // particle diameter
#define RHOP_HEAVY_APHTPB   5400.0              // particle solid density
#define MP_HEAVY_APHTPB     (PI/SIX)*RHOP_HEAVY_APHTPB*cube(DP_HEAVY_APHTPB)
#define QE_HEAVY_APHTPB     ZERO                // particle charge

/**********************************************************************
 * Definition of propellant constants.                                *
 **********************************************************************/

// Constants for AP-HPTB.
#define RHOS_APHTPB         1740.0              // propellant density
#define N_APHTPB            0.33                // burning constant
#define BETA_APHTPB         0.0000511646        // burning coefficient
#define TF_APHTPB           3060.0              // propellant flame temperature
#define TS_APHTPB           1130.0              // propellant surface temperature
#define ALPHAS_APHTPB       0.03                // percent concentraion of particles in propellant
#define CS_APHTPB           1510.0              // propellant specific heat

// Constants for QUICK-AP-HPTB.
#define RHOS_QUICK          17.4                // propellant density
#define N_QUICK             0.50                // burning constant
#define BETA_QUICK          0.000511646         // burning coefficient
#define TF_QUICK            3060.0              // propellant flame temperature
#define TS_QUICK            1130.0              // propellant surface temperature
#define ALPHAS_QUICK        0.03                // percent concentraion of particles in propellant
#define CS_QUICK            1510.0              // propellant specific heat

// Constants for PLAID-AP-HPTB.
#define RHOS_PLAID          8.70                // propellant density
#define N_PLAID             0.50                // burning constant
#define BETA_PLAID          0.001023292         // burning coefficient
#define TF_PLAID            3060.0              // propellant flame temperature
#define TS_PLAID            1130.0              // propellant surface temperature
#define ALPHAS_PLAID        0.03                // percent concentraion of particles in propellant
#define CS_PLAID            1510.0              // propellant specific heat

#endif // _SOLID_CONSTANTS_INCLUDED
