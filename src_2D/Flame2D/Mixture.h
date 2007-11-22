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

//! Standard Includes
#include <math.h>
#include <iostream>
using namespace std;

// FLAME2D Specific headers
#include "../Physics/NASAData/NASARP1311data.h"

//! \def STATIC_NUMBER_OF_SPECIES
//! If you define this variable, the number of species will be
//! predetermined for faster calculations.., however it is not as general 
//#define STATIC_NUMBER_OF_SPECIES 36

//Temperature convergence tolerance
#define CONV_TOLERANCE  1e-8

//Reference temperature for polytropic heat ratio mixture gamma [K]
#define REFERENCE_TEMPERATURE 200.0

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

  //@{ @name constructors/destructors
  Mixture() : T(0.0), MW(0.0), Cp(0.0), hs(0.0), 
	      hf(0.0), mu(0.0), kappa(0.0), g_ref(0.0)
  { Allocate(); for (int i=0; i<ns; i++) diff[i] = 0.0; }
  
  Mixture(const Mixture &M) : 
    T(M.temperature()), MW(M.molarMass()), Cp(M.heatCapacity_p()), 
    hs(M.enthalpySens()), hf(M.heatFormation()), mu(M.viscosity()), 
    kappa(M.thermalCond()), g_ref(M.heatRatioRef())
  { Allocate(); M.getDiffCoefs(diff); }

  ~Mixture() { Deallocate(); }

  void Copy( const Mixture &M ) {
    T = M.temperature();
    MW = M.molarMass();
    Cp = M.heatCapacity_p();
    hs = M.enthalpySens();
    hf = M.heatFormation();
    g_ref = M.heatRatioRef();
    mu = M.viscosity();
    kappa = M.thermalCond();
    M.getDiffCoefs(diff);
  }
  //@}

  //! Static setup function
  static void setMixture(const int &n, 
			 const string *names,
			 const char *PATH,
			 const double* Sc, 
			 const int &trans_data_flag,
			 const bool &constant_schmidt);

  //! set the mixture state
  void setState_TY(const double &T, const double* Y);
  

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
  double heatRatioRef(void) const {return g_ref;};
  //@}


  //@{ @name Static accessors
  static int nSpecies(void) {return ns;};
  static double LowTempRange(void) {return Tmin;};
  static double HighTempRange(void) {return Tmax;};
  //@}


  /***************** Mixing Rules ********************
    The following constructors return "total" physical
    parameters based on mixture rules for each.
  ****************************************************/
  //@{ @name Static versions
  static double molarMass( const double* Y );
  static double gasConstant( const double* Y );   
  static double heatCapacity_p( const double &Temp, const double* y );
  static double heatCapacity_v( const double &Temp, const double* y );
  static double heatRatio( const double &Temp, const double* y );
  static double heatRatioRef( const double* y );
  static double internalEnergy( const double &Temp, const double* y );
  static double internalEnergyRef( const double &Temp, const double* y );
  static double internalEnergySens( const double &Temp, const double* y );
  static double enthalpy( const double &Temp, const double* y );
  static double enthalpyRef( const double &Temp, const double* y );
  static double enthalpySens( const double &Temp, const double* y );
  static double enthalpyPrime( const double &Temp, const double* y );
  static double viscosity(const double &Temp, const double* y);
  static double viscosityPrime(const double &Temp, const double* y);
  static double thermalCond(const double &Temp, const double* y);
  static double gammaGuess(const double* y);
  static double temperature(double &h_s, const double* y);
  static double speciesGibbsFree(const double &Temp, const int &species);
  static double speciesDiffCoeff(const double &Temp, const double* y, const int &i);
  static double schmidt(const double &Temp, const double &rho, const double* y, const int &i);
  static double prandtl(const double &Temp, const double* y);
  static double lewis(const double &Temp, const double &rho, const double* y, const int &i);
//    double a(void);                 //speed of sound
//    double diedip() const;
//    double diedirho() const;
//    Vector2D thermal_diffusion(void) const;
  //@}

  //@{ @name local versions
  double gasConstant(void) const;   
  double heatCapacity_v(void) const;
  double heatRatio(void) const;
  double internalEnergy(void) const;
  double internalEnergyRef(void) const;
  double internalEnergySens(void) const;
  double enthalpy(void) const;
  double enthalpyRef(void) const;
  double enthalpyPrime(void) const;
  double gammaGuess(void) const;
  double schmidt(const double &rho, const int &i) const;
  double prandtl(void) const;
  double lewis(const double &rho, const int &i) const;
  //@}


  /**
   * Private members
   */
private:
  
  //! allocate / deallocate for objects
  void Allocate();
  void Deallocate();

  //! allocate / deallocate for static memory
  static void AllocateStatic();
  static void DeallocateStatic();



  /**
   * Private Objects
   */
private:
  
  double  T;     //!< Temperature [K]
  double  MW;    //!< molecular weight [kg/kmole]
  double  Cp;    //!< Heat Capacity (const Pressure) [J/(kg*K)]
  double  hs;    //!< Sensible enthalpy [J/kg]
  double  hf;    //!< heat of formation [J/kg]
  double  mu;    //!< Viscosity [kg/(m*s), N*s/m^2]
  double  kappa; //!< Thermal Conductivity [N/(s*K), W.(m*K)]
  double  g_ref; //!< Polytropc mixture heat ratio evaluated @ T_ref

  //!< species diffusion coefficient [m^2/s]
#ifdef STATIC_NUMBER_OF_SPECIES
  double  diff[STATIC_NUMBER_OF_SPECIES];
#else 
  double* diff;
#endif

  //@{ @name Static Variaables
  static int ns;                   //!< number of species
  static double Tmin;              //!< min temperature [k]
  static double Tmax;              //!< max temperature [k]
  static NASARP1311data *specdata; //!< Global Species Data
  static double* Sc_ref;           //!< reference species schmidt numbers
  static bool isConstSchmidt;   //!< flag indicating whether constant schmidt
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

/****************************************************
 * Allocator/deallocator for static memory.
 ****************************************************/
inline void Mixture :: AllocateStatic() {
  DeallocateStatic();
  if (ns>0) { 
    specdata = new NASARP1311data[ns]; 
    Sc_ref = new double[ns]; 
  }
};

inline void Mixture :: DeallocateStatic() { 
  if (specdata!=NULL) { delete[] specdata; specdata = NULL; } 
  if (Sc_ref!=NULL) { delete[] Sc_ref; Sc_ref = NULL; } 
};



/****************************************************
 * Mixture gas constant [kg/mol]
 ****************************************************/
inline double Mixture :: gasConstant(void) const {
  return (R_UNIVERSAL/MW);
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
    cerr << "Mixture::schmidt() - Not implemented yet.";
    exit(-1);
    //return mu/(rho*diff[i]);
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
    cerr << "Mixture::lewis() - Not implemented yet.";
    exit(-1);
    //return kappa / ( Cp*rho*diff[i] );
  }
}


#endif // _MIXTURE_INCLUDED
