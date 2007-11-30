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

//! Cantera libraries
#include <cantera/Cantera.h>      // main include
#include <cantera/IdealGasMix.h>  // reacting, ideal gas mixture class
#include <cantera/transport.h>
#include <cantera/equilibrium.h>

//! Standard Includes
#include <math.h>
#include <iostream>
using namespace std;

//Reference temperature for polytropic heat ratio mixture gamma [K]
#define TREF 298.15
#define PREF 101325.0

// If you define this variable, the number of species will be
// predetermined for faster calculations.., however it is not as general 
#define STATIC_NUMBER_OF_SPECIES 2 //2 AIR, 6 2STEP_CH4

/////////////////////////////////////////////////////////////////////
/// FUNCTION ROTOTYPES
/////////////////////////////////////////////////////////////////////

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
  Mixture() : T(0.0), MW(0.0), Cp(0.0), hs(0.0), 
	      hf(0.0), mu(0.0), kappa(0.0)
  { Nullify(); Allocate(); for (int i=0; i<ns; i++) diff[i] = 0.0; }
  
  Mixture(const Mixture &M) : T(M.temperature()),
    MW(M.molarMass()), Cp(M.heatCapacity_p()), 
    hs(M.enthalpySens()), hf(M.heatFormation()), mu(M.viscosity()), 
    kappa(M.thermalCond())
  { Nullify(); Allocate(); M.getDiffCoefs(diff); }

  ~Mixture() { Deallocate(); }

  void Copy( const Mixture &M ) {
    MW = M.temperature();
    MW = M.molarMass();
    Cp = M.heatCapacity_p();
    hs = M.enthalpySens();
    hf = M.heatFormation();
    mu = M.viscosity();
    kappa = M.thermalCond();
    M.getDiffCoefs(diff);
  }
  //@}

  /************ Accessors ****************************/

  //@{ @name Public accessors
  double temperature(void) const { return T; };
  double molarMass(void) const { return MW; };
  double heatCapacity_p(void) const { return Cp; };
  double enthalpySens(void) const { return hs; };
  double heatFormation(void) const { return hf; };
  double viscosity(void) const { return mu; };
  double thermalCond(void) const { return kappa; };
  double speciesDiffCoef( const int &i ) const 
  { return diff[i]; };
  void getDiffCoefs( double *d ) const
  { for (int i=0; i<ns; i++) d[i] = diff[i]; };
  //@}


  //@{ @name Static accessors
  static int nSpecies(void) {return ns;};
  static double molarMass(const int&i) {return M[i];};
  static double heatFormation(const int&i) {return Hform[i];};
  static double LowTempRange(void) {return Tmin;};
  static double HighTempRange(void) {return Tmax;};
  //@}


  /************ Setup Functions **********************/
  //! allocate / deallocate for static memory
  static void AllocateStatic();
  static void DeallocateStatic();

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
  
public:
  /***************** Mixing Rules ********************
    The following constructors return "total" physical
    parameters based on mixture rules for each.
  ****************************************************/
  static double molarMass( const double* Y );
  static double heatFormation( const double* Y );
  double gasConstant(void) const;   
  static double gasConstant( const double* Y );   
  static double heatCapacity_p( const double &Temp, const double* y );
  double heatCapacity_v(void) const;
  static double heatCapacity_v( const double &Temp, const double* y );
  double heatRatio(void) const;
  static double heatRatio( const double &Temp, const double* y );
  double internalEnergy() const;
  static double internalEnergy( const double &Temp, const double* y );
  double internalEnergySens() const;
  static double internalEnergySens( const double &Temp, const double* y );
  double enthalpy(void) const;
  static double enthalpy( const double &Temp, const double* y );
  static double enthalpySens( const double &Temp, const double* y );
  double enthalpyPrime(void) const;
  static double enthalpyPrime( const double &Temp, const double* y );
  static double viscosity(const double &Temp, const double* y);
  static double thermalCond(const double &Temp, const double* y);
  static double speciesGibbsFree(const double &Temp, const int &species);
  static double speciesDiffCoeff(const double &Temp, const double* y, const int &i);
  double schmidt(const double &rho, const int &i) const;
  static double schmidt(const double &Temp, const double &rho, const double* y, const int &i);
  double prandtl(void) const;
  static double prandtl(const double &Temp, const double* y);
  double lewis(const double &rho, const int &i) const;
  static double lewis(const double &Temp, const double &rho, const double* y, const int &i);

  static double temperature(const double& Press, const double* y);

  void getDihdDc( const double &Press, 
		  const double* y, 
		  double* dh ) const;
  static void getDihdDc(const double &Temp, 
			const double &Press, 
			const double* y, 
			double* dh);

  /***************** Reaction Rates ******************
    The following functions use CANTERA to compute    
    chemical reaction rates and equilibrium states.
  ****************************************************/
  static void equilibrate_HP( double &Temp, double &Press, double* y );
  static void equilibrate_TP( double &Temp, double &Press, double* y );
  void getRates( const double &Press, 
		 const double* y, 
		 double* rr ) const;
  static void getRates( const double &Temp, 
			const double &Press, 
			const double* y, 
			double* rr );

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
  
  double   T;    //!< Temperature [K]
  double  MW;    //!< molecular weight [kg/kmole]
  double  Cp;    //!< Heat Capacity (const Pressure) [J/(kg*K)]
  double  hs;    //!< Sensible enthalpy [J/kg]
  double  hf;    //!< heat of formation [J/kg]
  double  mu;    //!< Viscosity [kg/(m*s), N*s/m^2]
  double  kappa; //!< Thermal Conductivity [N/(s*K), W.(m*K)]

  //!< species diffusion coefficient [m^2/s]
#ifdef STATIC_NUMBER_OF_SPECIES
  double  diff[STATIC_NUMBER_OF_SPECIES];
#else 
  double* diff;
#endif

  //@{ @name Static Variaables
#ifdef STATIC_NUMBER_OF_SPECIES
  static const int ns = STATIC_NUMBER_OF_SPECIES;
#else
  static int       ns; //!< Number of species
#endif
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
  static IdealGasMix* ct_gas;     //!< the Cantera IdealGasMix object
  static Transport* ct_trans;     //!< the Cantera transport object

  //! Static storage
  static double *r, *r0; 
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
  if (ns>0) { diff = new double[ns]; }
#endif
  };

inline void Mixture :: Deallocate() { 
#ifndef STATIC_NUMBER_OF_SPECIES
  if (diff!=NULL) { delete[] diff; diff = NULL; } 
#endif
};

inline void Mixture :: Nullify() { 
#ifndef STATIC_NUMBER_OF_SPECIES
  diff = NULL;
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
 *
 * \param Temp Mixture temperature [k]
 * \param y Mixture species mass fractions.
 ****************************************************/
inline void Mixture :: setState_TPY(const double &Temp, 
				    const double &Press, 
				    const double* y)
{
  T = Temp;
  hf = heatFormation(y);

  ct_gas->setMassFractions_NoNorm(y);
  ct_gas->setPressure(Press);
  ct_gas->setTemperature(Temp);
  MW = ct_gas->meanMolecularWeight();
  hs = ct_gas->enthalpy_mass() - hf;
  Cp = ct_gas->cp_mass();
  mu = ct_trans->viscosity();
  kappa = ct_trans->thermalConductivity();
  if(isConstSchmidt){
    double rho( Press / (gasConstant()*T) );
    for (int i=0; i<ns; i++) diff[i] = mu/(rho*Sc_ref[i]);
  } else {
    ct_trans->getMixDiffCoeffs(diff);
  }


}

inline void Mixture :: setState_DPY(const double &rho, 
				    const double &Press, 
				    const double* y)
{
  MW = 0.0; hf = 0.0;
  for(int i=0; i<ns; i++){
    MW += y[i]/M[i];
    hf += y[i]*Hform[i];
  }
  MW = 1.0/MW;
  // MW = molarMass(y);
  // hf = heatFormation(y);

  T = Press / (gasConstant()*rho);

  ct_gas->setMassFractions_NoNorm(y);
  ct_gas->setPressure(Press);
  ct_gas->setTemperature(T);
  hs = ct_gas->enthalpy_mass() - hf;
  Cp = ct_gas->cp_mass();
  mu = ct_trans->viscosity();
  kappa = ct_trans->thermalConductivity();
  if(isConstSchmidt){
    for (int i=0; i<ns; i++) diff[i] = mu/(rho*Sc_ref[i]);
  } else {
    ct_trans->getMixDiffCoeffs(diff);
  }


}

inline void Mixture :: setState_DEY(const double &rho, const double &e, 
				    const double* y)
{
  hf = heatFormation(y);

  ct_gas->setMassFractions_NoNorm(y);
  ct_gas->setState_UV(e, 1.0/rho);
  T = ct_gas->temperature();
  MW = ct_gas->meanMolecularWeight();
  hs = ct_gas->enthalpy_mass() - hf;
  Cp = ct_gas->cp_mass();
  mu = ct_trans->viscosity();
  kappa = ct_trans->thermalConductivity();
  if(isConstSchmidt){
    for (int i=0; i<ns; i++) diff[i] = mu/(rho*Sc_ref[i]);
  } else {
    ct_trans->getMixDiffCoeffs(diff);
  }


}


inline void Mixture :: setState_DHY(const double &rho, const double &h, 
				    const double* y)
{
  hf = heatFormation(y);

  ct_gas->setMassFractions_NoNorm(y);
  setState_DH(rho, h);
  T = ct_gas->temperature();
  MW = ct_gas->meanMolecularWeight();
  hs = ct_gas->enthalpy_mass() - hf;
  Cp = ct_gas->cp_mass();
  mu = ct_trans->viscosity();
  kappa = ct_trans->thermalConductivity();
  if(isConstSchmidt){
    for (int i=0; i<ns; i++) diff[i] = mu/(rho*Sc_ref[i]);
  } else {
    ct_trans->getMixDiffCoeffs(diff);
  }


}


inline void Mixture :: setState_DH(const double &rho, const double& h, 
				   const double tol) {
  double dt, h0;
  ct_gas->setDensity(rho);
  
  // Newton iteration
  int n;
  for (n = 0; n < 500; n++) {
    h0 = ct_gas->enthalpy_mass();
    dt = (h - h0)/ct_gas->cp_mass();
    // limit step size to 100 K
    if (dt > 100.0) dt = 100.0;
    else if (dt < -100.0) dt = -100.0; 
    ct_gas->setTemperature(ct_gas->temperature() + dt);
    if (fabs(dt) < tol) {
      return;
    }
  }
  cerr << endl 
       << "Mixture::setState_DH() - No convergence. dt = " << dt 
       << ". H = " << h 
       << ". n = " << n 
       << endl;
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
  ct_gas->setState_TPY(T, Press, y);
  ct_gas->getNetProductionRates(rr);
  for (int i=0; i<ns; i++){
    rr[i] = rr[i]*M[i];
  }
}

/****************************************************
 * Derivative of species h wrt to mass fraction
 ****************************************************/
inline void Mixture :: getDihdDc(const double &Press, 
				 const double* y, 
				 double* dh) const {
  ct_gas->setState_TPY(T, Press, y);
  ct_gas->getEnthalpy_RT(dh); // -> h = hs + hf
  for(int i=0; i<ns; i++)
    dh[i] = ( dh[i]*(Cantera::GasConstant/M[i])*T -
	      Cp*T*MW/M[i] );
}

/****************************************************
 * Mixture gas constant [kg/mol]
 ****************************************************/
inline double Mixture :: gasConstant(void) const {
  return (Cantera::GasConstant/MW);
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

/****************************************************
 * Derivative of specific enthalpy dh/dT
 * actually is just Cp as Cp = (dh/dT)_p
 ****************************************************/
inline double Mixture :: enthalpyPrime(void) const {
  return Cp;
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


#endif // _MIXTURE_INCLUDED
