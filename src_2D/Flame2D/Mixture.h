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
//! Standard Includes
#include <math.h>
#include <iostream>
using namespace std;

//! \def STATIC_NUMBER_OF_SPECIES
//! If you define this variable, the number of species will be
//! predetermined for faster calculations.., however it is not as general 
//#define STATIC_NUMBER_OF_SPECIES 36

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

  //! Public members
public:

  /**
   * Default constructors / destructors
   */
  //@{ 
  //! constructor
  Mixture() : T(0.0), MW(0.0), Cp(0.0), h(0.0), 
	      mu(0.0), kappa(0.0)
  { Allocate(); ZeroDiff(); }
  
  //! destructor
  ~Mixture() { Deallocate(); }

  //! copy function
  void Copy( const Mixture &M ) {
    T = M.temperature();
    MW = M.molMass();
    Cp = M.specificHeat_p();
    h = M.enthalpy();
    mu = M.viscosity();
    kappa = M.thermalCond();
    M.getDiffCoefs(diff);
  }
  //@}

  /**
   * Public accessors
   */
  //@{ 
  double temperature() { return T; };
  double molarMass() { return MW; };
  double specificHeat_p() { return Cp; };
  double enthalpy() { return h; };
  double viscosity() { return mu; };
  double thermalCond() { return kappa; };
  double speciesDiffCoef( const int &i ) { return diff[i]; };
  void getDiffCoefs( const double *d ) 
  { for (int i=0; i<ns; i++) d[i] = diff[i]; };
  //@}



  //! Private members
private:
  
  //! allocate / deallocate
  void Allocate() {
#ifdef STATIC_NUMBER_OF_SPECIES
    if( STATIC_NUMBER_OF_SPECIES < ns) {
      cerr << "\n ERROR, Mixture::Allocate() - Built using static species with "
	   << STATIC_NUMBER_OF_SPECIES 
	   << " species predefined, asking for " << ns 
	   << endl; 
      exit(1); 
    }
#else 
    deallocate();
    if (ns>0) { diff = new double[ns]; }
#endif
  };

  void Deallocate() 
  { 
#ifdef STATIC_NUMBER_OF_SPECIES
#else 
    if (diff!=NULL) { delete[] diff; diff = NULL; } 
#endif
  };



  //! Private Objects
private:
  
  double  T;     //!< Temperature [K]
  double  MW;    //!< molecular weight [kg/kmole]
  double  Cp;    //!< Heat Capacity (const Pressure) [J/(kg*K)]
  double  h;     //!< Specific enthalpy [J/kg]
  double  mu;    //!< Viscosity [kg/(m*s), N*s/m^2]
  double  kappa; //!< Thermal Conductivity [N/(s*K), W.(m*K)]

  //!< species diffusion coefficient [m^2/s]
#ifdef STATIC_NUMBER_OF_SPECIES
  double  diff[STATIC_NUMBER_OF_SPECIES];
#else 
  double* diff;
#endif

  //@{ @name Static Variaables
  static int ns;      //!< number of species
  static double Tmin; //!< min temperature [k]
  static double Tmax; //!< max temperature [k]
  static NASARP1311data *specdata; //!< Global Species Data
  //@}


}
