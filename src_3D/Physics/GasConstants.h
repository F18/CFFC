/* GasConstants.h:  Header file defining useful gas-dynmic constants. */

#ifndef _GAS_CONSTANTS_INCLUDED
#define _GAS_CONSTANTS_INCLUDED

/* Define some useful physical constants. */

#define AVOGADRO            6.022137e23
#define BOLTZMANN           1.380658e-23
#define R_UNIVERSAL         8.314511
#define E_CHARGE            1.602189246e-19
#define ESU_CONSTANT        4.8e-13
#define GRAVITY_CONSTANT    6.672041e-11
#define VON_KARMAN_CONSTANT 0.41

/* Define standard atmospheric conditions
   (Sea Level US Standard Atmosphere). */

#define	PRESSURE_STDATM     101325.00
#define	DENSITY_STDATM      1.225
#define TEMPERATURE_STDATM  288.1600

/* Define various gas and ion types. */

#define GAS_AIR           1001
#define GAS_A             1002
#define GAS_CO            1003
#define GAS_CO2           1004
#define GAS_CH4           1005
#define GAS_H             1006
#define GAS_H2            1007
#define GAS_HE            1008
#define GAS_H2O           1009
#define GAS_N2            1010
#define GAS_O             1011
#define GAS_O2            1012
#define GAS_APHTPB        1013        // For solid propellent rocket motor gas.
#define GAS_e             2000

#define ION_H             2001
#define ION_HE            2002
#define ION_O             2003
#define ION_LOW_MASS      2004
#define ION_MED_MASS      2005
#define ION_HI_MASS       2006
#define ION_HUGE_MASS     2007
#define ION_HUGE_MASS2    2008
#define ION_e             3000

/* Define some specific heat ratios for various gases species
   (Zucrow and Hoffman, 1976). */

#define GAMMA_MONATOMIC     1.6666666667
#define GAMMA_DIATOMIC      1.40
#define GAMMA_POLYATOMIC    1.3333333333
#define GAMMA_TWO           2.00
#define GAMMA_ISOTHERMAL    1.01

#define GAMMA_AIR           1.40
#define GAMMA_A             1.6666666667
#define GAMMA_CO            1.398
#define GAMMA_CO2           1.288
#define GAMMA_CH4           1.304
#define GAMMA_H             1.6666666667
#define GAMMA_H2            1.405
#define GAMMA_HE            1.6666666667
#define GAMMA_H2O           1.329
#define GAMMA_N2            1.40
#define GAMMA_O             1.6666666667
#define GAMMA_O2            1.395
#define GAMMA_e             1.6666666667
#define GAMMA_APHTPB        1.208        // For solid propellent rocket motor gas.

#define GAMMA_ION_LOW_MASS   1.20
#define GAMMA_ION_MED_MASS   1.10
#define GAMMA_ION_HI_MASS    1.10
#define GAMMA_ION_HUGE_MASS  1.05
#define GAMMA_ION_HUGE_MASS2 1.05

/* Define some molecular weights of various gaseous species
   (Zucrow and Hoffman, 1976). */

#define MOLE_WT_AIR         28.9640
#define MOLE_WT_A           39.9440
#define MOLE_WT_CO          28.0100
#define MOLE_WT_CO2         44.0100
#define MOLE_WT_CH4         16.0430
#define MOLE_WT_H           1.00790
#define MOLE_WT_H2          2.01600
#define MOLE_WT_HE          4.00260
#define MOLE_WT_H2O         18.0160
#define MOLE_WT_N2          28.0130
#define MOLE_WT_O           15.9994
#define MOLE_WT_O2          32.0000
#define MOLE_WT_e           5.485802621e-04
#define MOLE_WT_APHTPB      26.14626101

#define MOLE_WT_ION_LOW_MASS   105.00
#define MOLE_WT_ION_MED_MASS   228.00
#define MOLE_WT_ION_HI_MASS    609.00
#define MOLE_WT_ION_HUGE_MASS  16951.00
#define MOLE_WT_ION_HUGE_MASS2 16241.00

/* Define some transport properties of various gases species. */

#define AIR_c1              0.000000521920
#define AIR_c2             -3.311320
#define AIR_c3              0.865351
#define AIR_c4           2365.270000
#define AIR_c5              0.000000

#define HE_c1               0.000000395509
#define HE_c2               0.000000
#define HE_c3               0.813299
#define HE_c4               0.000000
#define HE_c5               0.000000

#define H2_c1               0.000000172256
#define H2_c2               0.000000
#define H2_c3               0.808749
#define H2_c4               0.000000
#define H2_c5               0.000000

#define N2_c1               0.000000484270
#define N2_c2              -0.829421
#define N2_c3               0.851275
#define N2_c4            1219.260000
#define N2_c5               0.000000

#define O2_c1               0.000000540279
#define O2_c2              -0.164235
#define O2_c3               0.851096
#define O2_c4            2049.970000
#define O2_c5               0.000000

#define MU_APHTPB           0.0000819
#define KAPPA_APHTPB        0.1840
#define APHTPB_c1           0.000000
#define APHTPB_c2           0.000000
#define APHTPB_c3           0.000000
#define APHTPB_c4           0.000000
#define APHTPB_c5           MU_APHTPB

// Euken's formula for the Prandtl number (Pr = mu cp/k).
inline double Pr(double g) {
  return 20.0*g/(39.0*g - 15.0);  
}

// JJG's mu correlation.
inline double mu_gottlieb(double c1, double c2, double c3, double c4, double c5, double T) {
  return c1*pow(T,1.5)/(c2 + pow(T,c3) + c4/T) + c5;
}

// JJG's kappa correlation.
inline double kappa_gottlieb(double c1, double c2, double c3, double c4, double c5,
			     double T, double g, double cp) {
  return mu_gottlieb(c1,c2,c3,c4,c5,T)*cp/Pr(g);  
}

inline double kappa_gottlieb(double mu, double g, double cp) {
  return mu*cp/Pr(g);  
}

#endif // _GAS_CONSTANTS_INCLUDED
