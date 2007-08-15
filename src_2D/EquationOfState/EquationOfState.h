/**********************************************************************
 * EquationOfState.h: Header file defining Curve-fits                 *
 * of Srinivasan et al. (1987) for a High Temperature Gas             *
 * Equation of State                                                  *
 **********************************************************************/

#ifndef EquationOfState_h
#define EquationOfState_h

#include <cstdlib>

#define F(x) x ## _

//Datatype correlations
// C++      ref     Fortran
typedef int            integer;
typedef double         Freal;

//  The Fotran code in tgas.f and ugas.f is single precision given that
//  they date from the 1980's. Both the Alpha Fortran compiler (fort) and
//  Itanium Fortran compiler (ifort) have a -r8 option which can compile
//  the Fotran code in double precision.  Unfortunately I could not find
//  such an option for the GNU Fotran 77 compiler (g77). Professor Groth
//  says that we should convert the Fortran code to double precision by
//  hand. This would involve changing variable declarations to real*8,
//  changing 1.0E00 to 1.0D00, changing functions such as AMIN() to
//  DMIN(), changing the upper limit from 1.00E37 to the upper limit for
//  double precision, and so on. Professor Groth volunteered to do this
//  given his knowledge of Fortran. In the meantime I attempted a quick
//  fix along the lines of:
//  
//  	inline double Tgas_temp(double p, double rho) 
//  	{
//  		int IERROR;
//  	#ifdef _i386
//  		float temp;
//  		float p_single = static_cast<float>(p);
//  		float rho_single = static_cast<float>(rho);
//  		F(tgas3)(&p_single, &rho_single, &temp, &IERROR);
//  	#else
//  		double temp;
//  		F(tgas3)(&p, &rho, &temp, &IERROR);
//  	#endif
//  		return temp; // auto conversion to double for i386
//  	}
//  
//  I discovered that while the "forward" functions (such as tgas3)
//  worked, Professor Groth's "reverse" function tgas5 returned 901 in
//  IERROR where 901 is caused by:
//  
//  	DETMNT=DG1DD*DG2DE-DG1DE*DG2DD
//  	IF (ABS(DETMNT).LT.TOLER2) THEN
//  		IERROR=-901
//  		RETURN
//  	END IF
//   
//  I don't know why. Since I only wanted to run on i386 to see the output
//  of the profiler I gave up here. 
//  
//  But there is a fundamental issue that all the polynomial constants in
//  the tgas and ugas functions are single precision (most have exactly
//  six significant digits). To really make the code double precision, one
//  would need to redo the of Srinivasan et al (which was to fit
//  polynomials to the actual thermodynamic behaviour of air) but this
//  time in double precision from the ground up.
//
//    -- Alistair Wood. Wed Aug 08 2007.

extern "C"
{

  //TGAS subroutines

	// values for MFLAG
#define TGAS1_RETURN_PRESSURE  0
#define TGAS1_RETURN_PRESSURE_AND_SOUND_SPEED  1
#define TGAS1_RETURN_PRESSURE_AND_TEMPERATURE  2
#define TGAS1_RETURN_PRESSURE_TEMPERATURE_AND_SOUND_SPEED  3
  void F(tgas1)(Freal *E,Freal *R,Freal *P,Freal *A,Freal *T,integer *MFLAG,integer *IERROR);

  // TGAS2 returns equilibrium entropy
  void F(tgas2)(Freal *E,Freal *R, Freal *S,integer *IERROR);
  // TGAS3 returns equilibrium temperature
  void F(tgas3)(Freal *P, Freal *RHO, Freal *T, integer *IERROR);
  // TGAS4 returns equilibrium specific enthalpy
  void F(tgas4)(Freal *P, Freal *RHO, Freal *H, integer *IERROR);
  // TGAS5 returns density, specific internal energy, and sound speed
  void F(tgas5)(Freal *PRE, Freal *TEMP, Freal *DEN, Freal *EINT, Freal *SSP, integer *IERROR);
  // TGAS6 returns specific heats at constant pressure and volume and specific heat ratio
  void F(tgas6)(Freal *PRE, Freal *TEMP, Freal *CP, Freal *CV, Freal *GAM, integer *IERROR);

  //UGAS subroutines for transport quantities
  void F(ugas1)(Freal *T,Freal *RHO,Freal *MU,integer *IERROR);

  void F(ugas2)(Freal *T,Freal *RHO,Freal *PR,integer *IERROR);

  void F(ugas4)(Freal *E,Freal *RHO,Freal *K,integer *IERROR);

}

/**********************************************************************
 * Tgas_p -- Pressure.                         *
 **********************************************************************/
inline double Tgas_p(double inte, double rho)
{

#ifdef _i386
	cerr << "\nCurve fits do not work on i386. Returning.\n";
	exit(1);
#endif

  int MFLAG = TGAS1_RETURN_PRESSURE, IERROR;
  double ss, p, temp;
  F(tgas1)(&inte, &rho, &p, &ss, &temp, &MFLAG, &IERROR);
  return p; 
}

inline double Tgas_a_from_e_rho(double inte, double rho)
{
  int MFLAG = TGAS1_RETURN_PRESSURE_AND_SOUND_SPEED, IERROR;
  double ss, p, temp;
  F(tgas1)(&inte, &rho, &p, &ss, &temp, &MFLAG, &IERROR);
  return ss;
}

inline double Tgas_T_from_e_rho(double inte, double rho)
{
  int MFLAG = TGAS1_RETURN_PRESSURE_AND_TEMPERATURE, IERROR;
  double ss, p, temp;
  F(tgas1)(&inte, &rho, &p, &ss, &temp, &MFLAG, &IERROR);
  return temp;
}

/**********************************************************************
 * Tgas_temp -- Gas temperature.                                      *
 **********************************************************************/
inline double Tgas_temp(double p, double rho) 
{
	int IERROR;
  double temp;
  F(tgas3)(&p, &rho, &temp, &IERROR);
  return temp;
}

//  I eliminated this function!! Hooray!! (Once I realised that e = h - p/rho ...)
//  -- Alistair Wood. Thu Aug 09 2007. 01:30.
//
//  /**********************************************************************
//   * Tgas_e -- Gas specific internal energy.          *
//   **********************************************************************/
//  inline double Tgas_e(double p, double temp) 
//  {
//  	int IERROR;
//    double inte, rho, ss;
//    F(tgas5)(&p, &temp, &rho, &inte, &ss, &IERROR);
//    return inte;
//  }

/**********************************************************************
 * Tgas_rho -- Gas Density.                           *
 **********************************************************************/
inline double Tgas_rho(double p, double temp) 
{
	int IERROR;
  double inte, rho, ss;
  F(tgas5)(&p, &temp, &rho, &inte, &ss, &IERROR);
  return rho;
}

//  I eliminated this function!! (Once I realised that e = h - p/rho ...)
//  -- Alistair Wood. Thu Aug 09 2007. 00:25.
//
//  /**********************************************************************
//   * Tgas_a -- Equilibrium Sound Speed.                 *
//   **********************************************************************/
//  inline double Tgas_a(double p, double temp) 
//  {
//  	int IERROR;
//    double inte, rho, ss;
//    F(tgas5)(&p, &temp, &rho, &inte, &ss, &IERROR);
//    return ss;
//  }


/**********************************************************************
 * Tgas_s -- Equilibrium Entropy.                     *
 **********************************************************************/
inline double Tgas_s(double inte, double rho) 
{
	int IERROR;
  double s;
  F(tgas2)(&inte,&rho,&s,&IERROR);
  return s;
}

/**********************************************************************
 * Tgas_h -- Equilibrium Specific Enthalpy.           *
 **********************************************************************/
inline double Tgas_h(double p, double rho) 
{
	int IERROR;
  double h;
  F(tgas4)(&p, &rho, &h, &IERROR);
  return h;
}

/**********************************************************************
 * Tgas_cp -- Specific Heat at const. P.              *
 **********************************************************************/
inline double Tgas_cp(double p, double temp) 
{
	int IERROR;
  double cp, cv, g;
  F(tgas6)(&p, &temp, &cp, &cv, &g, &IERROR);
  return cp;
}

//  /**********************************************************************
//   * Tgas_cv -- Specific Heat at cost. V.              *
//   **********************************************************************/
//  inline double Tgas_cv(double p, double temp) 
//  {
//    double cp, cv, g, IERROR;
//    F(tgas6)(&p, &temp, &cp, &cv, &g, &IERROR); 
//    return cv;
//  }
//
//  Tgas_cv is not used so I removed it. Alistair Wood. Wed Aug 01 2007.


//  /**********************************************************************
//   * Tgas_gamma -- Specific Heat RATIO.                                 *
//   **********************************************************************/
//  inline double Tgas_gamma(double p, double temp) 
//  {
//    double cp, cv, g, IERROR;
//    F(tgas6)(&p, &temp, &cp, &cv, &g, &IERROR); 
//    return g;
//  }
//
//  Tgas_gamma is not used so I removed it. Alistair Wood. Wed Aug 01 2007.

/***************UGAS SUBROUTINES***************************************/

/**********************************************************************
 * Tgas_mu -- Equilibrium Dynamic Viscosity.          *
 **********************************************************************/
inline double Tgas_mu( double temp, double rho) 
{
	int IERROR;
  double mu;
  F(ugas1)(&temp, &rho, &mu, &IERROR);
  return mu;
}

/**********************************************************************
 * Tgas_Pr -- Prandtl Number.          *
 **********************************************************************/
inline double Tgas_Pr( double temp, double rho) 
{
	int IERROR;
  double pr;
  F(ugas2)(&temp, &rho, &pr, &IERROR);
  return pr;
}

/**********************************************************************
 * Tgas_kappa -- Coefficient of Thermal Conductivity. *
 **********************************************************************/
inline double Tgas_kappa( double inte, double rho) 
{
	int IERROR;
  double kappa;
  F(ugas4)(&inte, &rho, &kappa, &IERROR);
  return kappa;
}

//**********************************************************************


#endif
