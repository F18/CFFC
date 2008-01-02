/////////////////////////////////////////////////////////////////////
///
/// \file Mixture.h
/// 
/// \author Marc R.J. Charest
/// 
/// \brief This header file contains the class definition for 
///        describing the thermodynamic and transport properties
///        of a ideal gas mixture.
///
/////////////////////////////////////////////////////////////////////
#ifndef _MIXTURE_INCLUDED 
#define _MIXTURE_INCLUDED

//! include densematrix
#include "../Math/Matrix.h"

//! Cantera libraries
#include <cantera/Cantera.h>      // main include
#include <cantera/IdealGasMix.h>  // reacting, ideal gas mixture class
#include <cantera/transport.h>
#include <cantera/equilibrium.h>

//! Standard Includes
#include <math.h>
#include <iostream>
using namespace std;

/////////////////////////////////////////////////////////////////////
/// Defines
/////////////////////////////////////////////////////////////////////
//Reference temperature for polytropic heat ratio mixture gamma [K]
#undef TREF
#define TREF 298.15
#undef PREF
#define PREF 101325.0

// Properties of air at reference state (TREF and )
#undef  T_SDATM
#define T_SDATM       288.1600      // K
#undef  H_AIR_SDATM
#define H_AIR_SDATM  -10044.3039815 // J/kg ( 0 @ TREF, PREF )
#undef  HF_AIR_SDATM
#define HF_AIR_SDATM  0.0           // J/kg ( 0 @ TREF, PREF )
#undef  MW_AIR_SDATM
#define MW_AIR_SDATM  28.8507321008 // kg/kmole
#undef  CP_AIR_SDATM
#define CP_AIR_SDATM  1008.8401778  // J/(kg K)


// Temperature convergence tolerance
// For explicit -> TOL ~ 1.0E-08 is good
// For implicit -> TOL ~ 1.0E-10 is necessary 
#undef CONV_TOLERANCE
#define CONV_TOLERANCE  1e-8
#undef NUM_ITERATIONS
#define NUM_ITERATIONS  25

// If you define this variable, the number of species will be
// predetermined for faster calculations.., however it is not as general 
#define STATIC_NUMBER_OF_SPECIES 5 //2 AIR, 6 2STEP_CH4

/////////////////////////////////////////////////////////////////////
/// CLASS DEFINITIONS
/////////////////////////////////////////////////////////////////////
/**
 * \class Mixture
 *
 * This class defines the thermodynamic and transport properties 
 * for an ideal gas mixture.  By storing the thermodynamic properties
 * in a seperate class, we can control how much we recalculate the 
 * mixture properties.
 *
 * The relevant parameters include:
 * \arg Number density [particles / kg gaseous mixture]
 * \arg Mass fraction [kg soot / kg gaseous mixture]
 *
 * Currently, only a two parameter model of Leung et al. (1991) 
 * is implemented.  Thus: 
 * \arg  scalar[0] -> soot mass frac
 * \arg  scalar[1] -> number density
 */

class Mixture {

  /**
   * Public members
   */
public:


  /************** Constructors/Destructors ***********/

  //@{ @name constructors/destructors
  Mixture() : T(T_SDATM), MW(MW_AIR_SDATM), 
	      Cp(CP_AIR_SDATM), hs(H_AIR_SDATM), 
	      hf(HF_AIR_SDATM), mu_OK(false), kappa_OK(false),
	      diff_OK(false), dhdy_OK(false)
  { Nullify(); Allocate(); }
  
  Mixture(const Mixture &M) : T(M.temperature()),
			      MW(M.molarMass()), Cp(M.heatCapacity_p()), 
			      hs(M.enthalpySens()), hf(M.heatFormation()), 
			      mu_OK(false), kappa_OK(false),
			      diff_OK(false), dhdy_OK(false)
  { Nullify(); Allocate(); }

  ~Mixture() { Deallocate(); }

  void Copy( const Mixture &M ) {
    T = M.temperature();
    MW = M.molarMass();
    Cp = M.heatCapacity_p();
    hs = M.enthalpySens();
    hf = M.heatFormation();
    mu_OK = false;
    kappa_OK = false;
    diff_OK = false;
    dhdy_OK = false;
  }
  //@}

  /************ Accessors ****************************/

  //@{ @name Public accessors
  double temperature(void) const { return T; };
  double molarMass(void) const { return MW; };
  double heatCapacity_p(void) const { return Cp; };
  double enthalpySens(void) const { return hs; };
  double heatFormation(void) const { return hf; };
  double viscosity(void) const { assert(mu_OK); return mu; };
  double thermalCond(void) const { assert(kappa_OK); return kappa; };
  double speciesDiffCoef( const int &i ) const 
  { assert(diff_OK); return diff[i]; };
  void getDiffCoefs( double *d ) const 
  { assert(diff_OK); for (int i=0; i<ns; i++) d[i] = diff[i]; };
  double DihdDiy( const int &i ) const 
  { assert(dhdy_OK); return dhdy[i]; };
  void getDihdDiy( double* dh ) const 
  { assert(dhdy_OK); for (int i=0; i<ns; i++) dh[i] = dhdy[i]; };
  double Phi(void) const { assert(dhdy_OK); return phi; };
  //@}


  //@{ @name Static accessors
  static int nSpecies(void) {return ns;};
  static int nReactions(void) {return nr;};
  static double molarMass(const int&i) {return M[i];};
  static double heatFormation(const int&i) {return Hform[i];};
  static double LowTempRange(void) {return Tmin;};
  static double HighTempRange(void) {return Tmax;};
  //@}


  /************ Setup Functions **********************/
  //! allocate / deallocate for static memory
  static void AllocateStatic();
  static void DeallocateStatic();

  //! returns the number of species for a mechanism
  static int getNumSpecies(const string &mech_name,
			   const string &mech_file);

  //! Static setup function
  static void setMixture(const string &mech_name,
			 const string &mech_file);
  static void setConstantSchmidt(const double* Sc);

  //! static cantera setup functions
  static void parse_mass_string( const string& massFracStr, 
				 double* massFracs);
  static void parse_schmidt_string( const string& schmidtStr, 
				    double* schmidt);
  static void parse_mole_string( const string& moleFracStr, 
				 double* moleFracs);
  static int speciesIndex(const string &sp);
  static string speciesName(const int &i) { return names[i]; };
  static void composition( const string& fuel_species, 
			   const double &phi,
			   double* massFracs);
  static string mechName(void) { return ct_mech_name; };

  //! set the mixture state
  void setState_TPY(const double &T, const double &Press, const double* Y);
  void setState_DPY(const double &rho, const double &Press, const double* Y);
  void setState_DEY(const double &rho, const double &e, const double* Y);
  void setState_DHY(const double &rho, const double &h, const double* y);
private:
  void setState_DH(const double &rho, const double& h, const double tol = 1.e-4);
  void setState_DE(const double &rho, const double& e, const double tol = 1.e-4);

public:
  //! update viscosity
  void updateViscosity(const double &rho, const double* y);
  //! update thermal conductivity and diffusion coefficients
  void updateTransport(const double &rho, const double* y);
  //! update dihdiy and phi
  void updateDihdDic( const double &rho, 
		      const double* y, 
		      const int NSm1 );

  /***************** Mixing Rules ********************
    The following constructors return "total" physical
    parameters based on mixture rules for each.
  ****************************************************/
  static double molarMass( const double* Y );
  static double heatFormation( const double* Y );
  double gasConstant(void) const;   
  static double heatCapacity_p(const double &Temp, const double &Press, 
			       const double* y);
  void getHeatCapacity_p( const double &Press, 
			  const double* y, 
			  double*cp ) const;
  double heatCapacity_v(void) const;
  double heatRatio(void) const;
  double internalEnergy() const;
  double internalEnergySens() const;
  double enthalpy(void) const;
  static double enthalpy(const double &Temp, const double &Press, 
			 const double* y);
  void getEnthalpy( const double &Press, const double* y, double*h ) const;
  void get_cp_and_h( const double &Press, 
		     const double* y, 
		     double*cp,
		     double*h) const;
  double schmidt(const double &rho, const int &i) const;
  double prandtl(void) const;
  double lewis(const double &rho, const int &i) const;

  /***************** Reaction Rates ******************
    The following functions use CANTERA to compute    
    chemical reaction rates and equilibrium states.
  ****************************************************/
  void equilibrate_HP( double &rho, double &Press, double* y );
  void equilibrate_TP( double &rho, double &Press, double* y );
  void getRates( const double &Press, 
		 const double* y, 
		 double* rr ) const;
  void dSwdU( DenseMatrix &dSdU,
	      const double &rho,
	      const double &Press,
	      const double* y,
	      const int offset,
	      const int NSm1) const;
  double dSwdU_max_diagonal( const double &rho,
			     const double &Press,
			     const double* y) const;

  /**************** Output Functions *************************/
  void output (ostream &out) const;

  /**************** Operators Overloading ********************/
  // Input-output operators.
  friend ostream& operator << (ostream &out_file, const Mixture &Mix) {
    Mix.output(out_file);
    return out_file;
  }

  /**
   * Private members
   */
private:
  
  //! allocate / deallocate for objects
  void Allocate();
  void Deallocate();
  void Nullify();


  /**
   * Private Objects
   */
private:
  
  double   T;      //!< Temperature [K]
  double  MW;      //!< molecular weight [kg/kmole]
  double  Cp;      //!< Heat Capacity (const Pressure) [J/(kg*K)]
  double  hs;      //!< Sensible enthalpy [J/kg]
  double  hf;      //!< heat of formation [J/kg]
  double  mu;      //!< Viscosity [kg/(m*s), N*s/m^2]
  double  kappa;   //!< Thermal Conductivity [N/(s*K), W.(m*K)]
  double  phi;     //!< The quantity phi
  bool mu_OK;      //!< Flag indicating whether viscosity is up-to-date
  bool kappa_OK;   //!< Flag indicating whether transport properties are up-to-date
  bool diff_OK;    //!< Flag indicating whether transport properties are up-to-date
  bool dhdy_OK;    //!< Flag indicating whether dihdic is up-to-date

#ifdef STATIC_NUMBER_OF_SPECIES
  double  diff[STATIC_NUMBER_OF_SPECIES]; //!< species diffusion coefficient [m^2/s]
  double  dhdy[STATIC_NUMBER_OF_SPECIES]; //!< derivative of enthalpy wrt to mass fracs
#else 
  double* diff;
  double* dhdy;
#endif

  //@{ @name Static Variaables
#ifdef STATIC_NUMBER_OF_SPECIES
  static const int ns = STATIC_NUMBER_OF_SPECIES;
#else
  static int       ns; //!< Number of species
#endif
  static int      nr;              //!< Number of reactions
  static double Tmin;              //!< min temperature [k]
  static double Tmax;              //!< max temperature [k]
  static double* Sc_ref;           //!< reference species schmidt numbers
  static double* M;                //!< species molar masses [kg/kmole]
  static double* Hform;            //!< species heats of formation [J/kg]
  static string* names;            //!< species names
  static bool isConstSchmidt;      //!< flag indicating whether constant schmidt

  //! Cantera objects
  static string ct_mech_name;     //!< Reaction mechanism file path
  static string ct_mech_file;     //!< Reaction mechanism file path
  static Cantera::IdealGasMix* ct_gas;     //!< the Cantera IdealGasMix object
  static Cantera::Transport* ct_trans;     //!< the Cantera transport object

  //! Static storage
  static double *r, *r0, *c; 
  //@}


};


/////////////////////////////////////////////////////////////////////
/// Member Functions
/////////////////////////////////////////////////////////////////////

/****************************************************
 * Allocator/deallocator for object memory.
 ****************************************************/
inline void Mixture :: Allocate() {
#ifndef STATIC_NUMBER_OF_SPECIES
  Deallocate();
  if (ns>0) { 
    diff = new double[ns]; 
    dhdy = new double[ns];
  }
#endif
};

inline void Mixture :: Deallocate() { 
#ifndef STATIC_NUMBER_OF_SPECIES
  if (diff!=NULL) { delete[] diff; diff = NULL; } 
  if (dhdy!=NULL) { delete[] dhdy; dhdy = NULL; } 
#endif
};

inline void Mixture :: Nullify() { 
#ifndef STATIC_NUMBER_OF_SPECIES
  diff = NULL;
  dhdy = NULL;
#endif
};

/****************************************************
 * Allocator/deallocator for static memory.
 ****************************************************/
inline void Mixture :: AllocateStatic() {
  if (ns>0) { 
    Sc_ref = new double[ns]; 
    names = new string[ns]; 
    M = new double[ns]; 
    Hform = new double[ns]; 
    r = new double[ns]; 
    r0 = new double[ns]; 
    c = new double[ns]; 
  }
};

inline void Mixture :: DeallocateStatic() { 
  if (ct_gas!=NULL) { delete ct_gas; ct_gas = NULL; } 
  if (ct_trans!=NULL) { delete ct_trans; ct_trans = NULL; } 
  if (Sc_ref!=NULL) { delete[] Sc_ref; Sc_ref = NULL; } 
  if (names!=NULL) { delete[] names; names = NULL; } 
  if (M!=NULL) { delete[] M; M = NULL; } 
  if (Hform!=NULL) { delete[] Hform; Hform = NULL; } 
  if (r!=NULL) { delete[] r; r = NULL; } 
  if (r0!=NULL) { delete[] r0; r0 = NULL; } 
  if (c!=NULL) { delete[] c; c = NULL; } 
};

/****************************************************
 * Output functions.
 ****************************************************/
inline void Mixture :: output (ostream &out) const {
  out << "Temperature = " << T << endl;
  out << "molecular weight = " << MW << endl;
  out << "Heat Capacity = " << Cp << endl;
  out << "Sensible enthalpy = " << hs << endl;
  out << "heat of formation = " << hf << endl;
  out << "Viscosity = " << mu << endl;
  out << "Thermal Conductivity = " << kappa << endl;
  out << "Diffusivities:" << endl;
  for (int i=0; i<ns; i++)
    out << speciesName(i) << " = " << diff[i] << endl;
}


/////////////////////////////////////////////////////////////////////
/// LOCAL MEMBER FUNCTIONS
/////////////////////////////////////////////////////////////////////

/****************************************************
 * Set the state of the mixture and update all the properties.
 * State is set based on temperature, pressure, and
 * species mass fractions
 ****************************************************/
inline void Mixture :: setState_TPY(const double &Temp, 
				    const double &Press, 
				    const double* y)
{
  // initialize
  T = Temp;
  hf = heatFormation(y);
  
  // set cantera state
  ct_gas->setMassFractions_NoNorm(y);
  ct_gas->setTemperature(Temp);
  ct_gas->setPressure(Press);

  // compute thermo / transport properties
  MW = ct_gas->meanMolecularWeight();
  hs = ct_gas->enthalpy_mass() - hf;
  Cp = ct_gas->cp_mass();

  // set out of date flags
  mu_OK = false;
  kappa_OK = false;
  diff_OK = false;
  dhdy_OK = false;

  //--------------------------------------------------------
  // NO NEED TO STORE THESE
  // kappa = ct_trans->thermalConductivity();
  // mu = ct_trans->viscosity();

  // compute diffusion coefficient
  // if(isConstSchmidt){
  //   double rho( Press / (gasConstant()*T) );
  //   for (int i=0; i<ns; i++) diff[i] = mu/(rho*Sc_ref[i]);
  // } else {
  //   ct_trans->getMixDiffCoeffs(diff);
  // }
  //--------------------------------------------------------


}

/****************************************************
 * Set the state of the mixture and update all the properties.
 * State is set based on density, pressure, and
 * species mass fractions
 ****************************************************/
inline void Mixture :: setState_DPY(const double &rho, 
				    const double &Press, 
				    const double* y)
{
  // compute mixture molecular weight and heat of formation
  MW = 0.0; hf = 0.0;
  for(int i=0; i<ns; i++){
    MW += y[i]/M[i];
    hf += y[i]*Hform[i];
  }
  MW = 1.0/MW;

  // compute the temperature
  T = Press / (gasConstant()*rho);

  // set the cantera gas state
  ct_gas->setMassFractions_NoNorm(y);
  ct_gas->setTemperature(T);
  ct_gas->setPressure(Press);

  // compute thermo / trans properties
  hs = ct_gas->enthalpy_mass() - hf;
  Cp = ct_gas->cp_mass();

  // set out of date flags
  mu_OK = false;
  kappa_OK = false;
  diff_OK = false;
  dhdy_OK = false;

  //--------------------------------------------------------
  // NO NEED TO STORE THESE
  // kappa = ct_trans->thermalConductivity();
  // mu = ct_trans->viscosity();

  // compute diffusion coefficient
  // if(isConstSchmidt){
  //   for (int i=0; i<ns; i++) diff[i] = mu/(rho*Sc_ref[i]);
  // } else {
  //   ct_trans->getMixDiffCoeffs(diff);
  // }
  //--------------------------------------------------------


}

/****************************************************
 * Set the state of the mixture and update all the properties.
 * State is set based on density, energy, and
 * species mass fractions
 ****************************************************/
inline void Mixture :: setState_DEY(const double &rho, const double &e, 
				    const double* y)
{
  // compute heat of formation
  hf = heatFormation(y);

  // set the cantera gas state
  ct_gas->setMassFractions_NoNorm(y);
  MW = ct_gas->meanMolecularWeight();
  setState_DE(rho, e, CONV_TOLERANCE);

  // compute the thermo / trans properties
  T = ct_gas->temperature();
  hs = ct_gas->enthalpy_mass() - hf;
  Cp = ct_gas->cp_mass();

  // set out of date flags
  mu_OK = false;
  kappa_OK = false;
  diff_OK = false;
  dhdy_OK = false;

  //--------------------------------------------------------
  // NO NEED TO STORE THESE
  // kappa = ct_trans->thermalConductivity();
  // mu = ct_trans->viscosity();

  // compute diffusion coefficient
  // if(isConstSchmidt){
  //   for (int i=0; i<ns; i++) diff[i] = mu/(rho*Sc_ref[i]);
  // } else {
  //   ct_trans->getMixDiffCoeffs(diff);
  // }
  //--------------------------------------------------------


}


/****************************************************
 * Set the state of the mixture and update all the properties.
 * State is set based on density, enthalpy, and
 * species mass fractions
 ****************************************************/
inline void Mixture :: setState_DHY(const double &rho, const double &h, 
				    const double* y)
{
  // compute heat of formation
  hf = heatFormation(y);

  // set the cantera gas state
  ct_gas->setMassFractions_NoNorm(y);
  MW = ct_gas->meanMolecularWeight();
  setState_DH(rho, h, CONV_TOLERANCE);

  // compute the thermo / trans properties
  T = ct_gas->temperature();
  hs = ct_gas->enthalpy_mass() - hf;
  Cp = ct_gas->cp_mass();

  // set out of date flags
  mu_OK = false;
  kappa_OK = false;
  diff_OK = false;
  dhdy_OK = false;

  //--------------------------------------------------------
  // NO NEED TO STORE THESE
  // kappa = ct_trans->thermalConductivity();
  // mu = ct_trans->viscosity();
  // 
  // // compute diffusion coefficient
  // if(isConstSchmidt){
  //   for (int i=0; i<ns; i++) diff[i] = mu/(rho*Sc_ref[i]);
  // } else {
  //   ct_trans->getMixDiffCoeffs(diff);
  // }
  //--------------------------------------------------------


}


/****************************************************
 * Set the state of the mixture and update all the properties.
 * State is set based on density, and enthalpy.
 * Note, mass fractions should already be set.
 ****************************************************/
inline void Mixture :: setState_DH(const double &rho, const double& h, 
				   const double tol) {
  double dt;
  ct_gas->setDensity(rho);
  ct_gas->setTemperature(T);

  // Newton iteration
  int n;
  for (n = 0; n < NUM_ITERATIONS; n++) {
    // time step
    dt = (h - ct_gas->enthalpy_mass())/ct_gas->cp_mass();
    // limit step size to 100 K
    // if (dt > 100.0) dt = 100.0;
    // else if (dt < -100.0) dt = -100.0;
    // increment T
    ct_gas->setTemperature(ct_gas->temperature() + dt);
    //check convergence
    if (fabs(dt) < tol) return;
  }
  cerr << endl 
       << "Mixture::setState_DH() - No convergence. dt = " << dt 
       << ". H = " << h 
       << ". n = " << n 
       << ". T = " << ct_gas->temperature()
       << endl 
       << flush;
}

/****************************************************
 * Set the state of the mixture and update all the properties.
 * State is set based on density, and energy.
 * Note, mass fractions should already be set.
 ****************************************************/
inline void Mixture :: setState_DE(const double &rho, const double& e, 
				   const double tol) {
  double dt;
  ct_gas->setDensity(rho);
  ct_gas->setTemperature(T);
  double Rmix( gasConstant() );

  // Newton iteration
  int n;
  for (n = 0; n < NUM_ITERATIONS; n++) {
    // time step
    dt = (e - ct_gas->intEnergy_mass())/(ct_gas->cp_mass()-Rmix);
    // limit step size to 100 K
    // if (dt > 100.0) dt = 100.0;
    // else if (dt < -100.0) dt = -100.0;
    // increment T
    ct_gas->setTemperature(ct_gas->temperature() + dt);
    //check convergence
    if (fabs(dt) < tol) return;
  }
  cerr << endl 
       << "Mixture::setState_DE() - No convergence. dt = " << dt 
       << ". E = " << e
       << ". n = " << n 
       << ". T = " << ct_gas->temperature()
       << endl
       << flush;
}

/****************************************************
 * Compute the remaining transport properties, i.e.
 * get the thermal conductivity and the species 
 * diffusivities.
 ****************************************************/
inline void Mixture :: updateViscosity(const double &rho, 
				       const double* y) {
  if (!mu_OK) {

    // set the state
    ct_gas->setMassFractions_NoNorm(y);
    ct_gas->setTemperature(T);
    ct_gas->setDensity(rho);

    mu = ct_trans->viscosity();
    mu_OK = true;
      
  } // endif

}

/****************************************************
 * Compute the remaining transport properties, i.e.
 * get the thermal conductivity and the species 
 * diffusivities.
 ****************************************************/
inline void Mixture :: updateTransport(const double &rho, 
				       const double* y) {
  if (!kappa_OK || !diff_OK || !mu_OK) {

    // set the state
    ct_gas->setMassFractions_NoNorm(y);
    ct_gas->setTemperature(T);
    ct_gas->setDensity(rho);

    // viscosity
    if (!mu_OK) {
      mu = ct_trans->viscosity();
      mu_OK = true;
    }

    // thermal conductivity
    if (!kappa_OK) {
      kappa = ct_trans->thermalConductivity();
      kappa_OK = true;
    }

    // compute diffusion coefficient
    if (!diff_OK) {
      if(isConstSchmidt){
	for (int i=0; i<ns; i++) diff[i] = mu/(rho*Sc_ref[i]);
      } else {
	ct_trans->getMixDiffCoeffs(diff);
      }
      diff_OK = true;
    }

  } // endif

}

/****************************************************
 * Compute the derivative of the species enthalpies
 * wrt the mass fractions.
 ****************************************************/
inline void Mixture ::  updateDihdDic( const double &rho, 
				       const double* y, 
				       const int NSm1 ) {
    
  if (!dhdy_OK) {
	
    // set the state
    ct_gas->setMassFractions_NoNorm(y);
    ct_gas->setTemperature(T);
    ct_gas->setDensity(rho);
    ct_gas->getEnthalpy_RT(dhdy); // -> h = hs + hf

    // get last species value
    // for h_k - h_N
    double hN( 0.0 );
    if (NSm1) hN = ( dhdy[ns-1]*(Cantera::GasConstant/M[ns-1])*T -
		     Cp*T*MW/M[ns-1] );
	
    // compute the rest and phi at the same time
    phi = 0.0;
    for(int i=0; i<ns-NSm1; i++) {
      dhdy[i] = ( ( dhdy[i]*(Cantera::GasConstant/M[i])*T -
		    Cp*T*MW/M[i] ) - hN );
      phi += y[i] * dhdy[i];
    }

    dhdy_OK = true;
	
  } // endif

}



/****************************************************
 * Mixture molecular mass [kg/mol]
 ****************************************************/
inline double Mixture :: molarMass( const double* y ) {
  double sum( 0.0 );
  for(int i=0; i<ns; i++){
    sum += y[i]/M[i];
  }
  return 1.0/sum;
}

/****************************************************
 * Mixture heat of formation [J/kg]
 ****************************************************/
inline double Mixture :: heatFormation( const double* y ) {
  double sum( 0.0 );
  for(int i=0; i<ns; i++){
    sum += y[i]*Hform[i];
  }
  return sum;
}

/****************************************************
 * Mixture gas constant [kg/mol]
 ****************************************************/
inline double Mixture :: gasConstant(void) const {
  return (Cantera::GasConstant/MW);
}

/****************************************************
 * Mixture Heat Capacity (const pressure) J/(kg*K)
 ****************************************************/
inline double Mixture :: heatCapacity_p(const double &Temp, 
					const double &Press, 
					const double* y) {
  ct_gas->setMassFractions_NoNorm(y);
  ct_gas->setTemperature(Temp);
  ct_gas->setPressure(Press);
  return ct_gas->cp_mass();
}

/****************************************************
 * Mixture Heat Capacity (const volume) J/(kg*K)
 ****************************************************/
inline double Mixture :: heatCapacity_v(void) const {
  return (Cp - gasConstant());
}

/****************************************************
 * Mixture Heat Ratio gamma
 ****************************************************/
inline double Mixture :: heatRatio(void) const {
  return (heatCapacity_p()/heatCapacity_v());
}

/****************************************************
 * Mixture Specific Internal Energy J/(kg)
 ****************************************************/
//! etotal = sensible & chemical
inline double Mixture :: internalEnergy(void) const {
  //(Enthalpy(Temp) - (R/mol_mass)*Temp)
  return (hs + hf - gasConstant()*T);
}

//! internal energy with no heat of formation included (sensible)
inline double Mixture :: internalEnergySens(void) const {
  //(Enthalpy(Temp) - (R/mol_mass)*Temp)
  return (hs - gasConstant()*T);
}

/****************************************************
 * Specific absolute enthalpy J/(kg)
 ****************************************************/
//! htotal = mass fraction * (hsensible + heatofform)
inline double Mixture :: enthalpy(void) const {
  return ( hs + hf );
}

//! htotal = mass fraction * (hsensible + heatofform)
inline double Mixture :: enthalpy(const double &Temp, 
				  const double &Press, 
				  const double* y) {
  ct_gas->setMassFractions_NoNorm(y);
  ct_gas->setTemperature(Temp);
  ct_gas->setPressure(Press);
  return ct_gas->enthalpy_mass();
}


/****************************************************
 * Individual Species enthalpy
 ****************************************************/
inline void Mixture :: getEnthalpy( const double &Press, 
				    const double* y, 
				    double*h ) const {
  ct_gas->setMassFractions_NoNorm(y);
  ct_gas->setTemperature(T);
  ct_gas->setPressure(Press);
  ct_gas->getEnthalpy_RT(h); // -> h = hs + hf
  for(int i=0; i<ns; i++) h[i] *= (Cantera::GasConstant/M[i])*T;
}

/****************************************************
 * Individual Species Specific Heat
 ****************************************************/
inline void Mixture :: getHeatCapacity_p( const double &Press, 
					  const double* y, 
					  double*cp ) const {
  ct_gas->setMassFractions_NoNorm(y);
  ct_gas->setTemperature(T);
  ct_gas->setPressure(Press);
  ct_gas->getCp_R(cp);
  for(int i=0; i<ns; i++) cp[i] *= (Cantera::GasConstant/M[i]);
}

/****************************************************
 * How about both at once!!!
 ****************************************************/
inline void Mixture :: get_cp_and_h( const double &Press, 
				     const double* y, 
				     double*cp,
				     double*h) const {
  ct_gas->setMassFractions_NoNorm(y);
  ct_gas->setTemperature(T);
  ct_gas->setPressure(Press);
  ct_gas->getEnthalpy_RT(h); // -> h = hs + hf
  ct_gas->getCp_R(cp);

  double R;
  for(int i=0; i<ns; i++) {
    R = Cantera::GasConstant/M[i];
    h[i] *= R*T;
    cp[i] *= R;
  }
}


/****************************************************
 * Schmidt
 ****************************************************/
inline double Mixture :: schmidt(const double &rho,
				 const int &i) const {
  if(isConstSchmidt){
    return Sc_ref[i];
  } else {
    return mu/(rho*diff[i]);
  }
}

/****************************************************
 * Prandtl
 ****************************************************/
inline double Mixture :: prandtl(void) const {
  //Pr = Cp*mu/k
  return Cp*mu/kappa;
}

/****************************************************
 * Lewis
 ****************************************************/
inline double Mixture :: lewis(const double &rho,
			       const int &i) const {
  if(isConstSchmidt) {
    return (kappa*Sc_ref[i])/( Cp*mu );
  } else {
    return kappa / ( Cp*rho*diff[i] );
  }
}

/************************************************************************
  Calculates the concentration time rate of change of species from
  primitive state W using the general law of mass action.
  U is the conserved state container for passing back the 
  source terms. ie. U.rhospec[i].c 

  W.SpecCon:  is the  species mass fractions concentrations
              of Chem2D_pState. (c_i*rho/M_i)   mol/m^3

  Return units are  kg/m^3*s ie. rho*omega (kg/m^3)*(1/s)

************************************************************************/
inline void Mixture :: getRates( const double &Press, 
				 const double* y, double* rr ) const {
  ct_gas->setMassFractions_NoNorm(y);
  ct_gas->setTemperature(T);
  ct_gas->setPressure(Press);
  //ct_gas->setState_TPY(T, Press, y);
  ct_gas->getNetProductionRates(rr);
  for (int i=0; i<ns; i++){
    rr[i] = rr[i]*M[i];
  }
}


/************************************************************************
  Calculates the equilibrium composition given an unburnt mixture
  using CANTERA.  This is only for CANTERA reaction types.  Here
  we hold enthalpy and pressure fixed.
************************************************************************/
inline void Mixture::equilibrate_HP( double &rho, double &Press, double* y ) {

  // set state and equilibrate
  ct_gas->setMassFractions_NoNorm(y);
  ct_gas->setTemperature(T);
  ct_gas->setDensity(rho);
  equilibrate( *ct_gas, "HP" );

  //get burnt mass fractions
  ct_gas->getMassFractions(y);

  // the temperature
  T = ct_gas->temperature();

  // density and pressure
  Press = ct_gas->pressure();
  rho = ct_gas->density();
  
} // end of ct_equilibrate


/************************************************************************
  Calculates the equilibrium composition given an unburnt mixture
  using CANTERA.  This is only for CANTERA reaction types.  Here 
  we hold temperature and pressure fixed.
************************************************************************/
inline void Mixture::equilibrate_TP( double &rho, double &Press, double* y ) {

  // set state and equilibrate
  ct_gas->setMassFractions_NoNorm(y);
  ct_gas->setTemperature(T);
  ct_gas->setDensity(rho);
  equilibrate( *ct_gas, "TP" );

  //get burnt mass fractions
  ct_gas->getMassFractions(y);

  // density and pressure
  Press = ct_gas->pressure();
  rho = ct_gas->density();

} // end of ct_equilibrate


#endif // _MIXTURE_INCLUDED
