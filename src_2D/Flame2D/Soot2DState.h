/////////////////////////////////////////////////////////////////////
///
/// \file SootState.h
/// 
/// \author Marc R.J. Charest
/// 
/// \brief This header file contains the class definition for 
///        describing the state of soot.
/// 
/// This header file contains the class definition for describing
/// the state of soot.  A simplified reaction mechanism for the 
/// formation, growth, and combustion of soot particles in laminar
/// nonpremixed flames is used.
///
/// \see K.M. Leung, R.P. Lindstedt, and W.P. Jones, "A simplified
///      reaction mechanism for soot formation in nonpremixed
///      flames," Combustion and Flame 87:289-305 (1991).
///
/// \todo Add check for soot important species.
///
/////////////////////////////////////////////////////////////////////
#ifndef _SOOT2D_STATE_INCLUDED 
#define _SOOT2D_STATE_INCLUDED

//! Standard Includes
#include <math.h>
#include <iostream>
using namespace std;

// CFFC includes
#include "../Math/Math.h"
#include "../Physics/GasConstants.h"


/////////////////////////////////////////////////////////////////////
/// Defines
/////////////////////////////////////////////////////////////////////

// Soot model Flags
enum SootFlags { SOOT2D_NONE = 0,                // No soot model
		 SOOT2D_LEUNG_ET_AL_1991 = 1 };  // Leung, Lindstedt, and Jones (1991)

const string SOOT_MODEL_NAMES[] = { "NONE",
				    "LEUNG_ET_AL_1991" };

const int NUM_SOOT_MODELS = 2;


// If you define this variable, the number of soot scalars will be
// predetermined for faster calculations.., however it is not as general 
#define STATIC_NUM_SOOT2D_VAR 2


/////////////////////////////////////////////////////////////////////
/// CLASS DEFINITIONS
/////////////////////////////////////////////////////////////////////

/**
 * \class Soot_State
 *
 * This class defines the state of soot and member functions used
 * to describe its evolution in time.  It is to be used with the 
 * well stirred and plug flow reactors and currently tracks its 
 * mass fraction and number density in a gaseous mixture.
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
class Soot2D_State {

private:
  
  /*********************** Objects **********************************/
  int      model; //! soot model flag

  //! Soot phase objects
  double      MW; //! Soot molar mass [kg / kmol]
  double     rho; //! Soot density [kg soot / m^3 soot]
  int       Cmin; //! Number of carbon atoms in the incipient carbon particles
  double      Ca; //! Aglomeration rate constant
  int         nv; //! number of state variables

  //! gaseous phase objects
  double* MW_gas; //! array of gas phase species mass fractions [kg/kmole]
  int     ns_gas; //! number of gasesous species

  //! species array indices
  int iC2H2, iO2, iH2, iCO;


public:

  /**************** Constructors / Destructors **********************/
  //! Default constructors
  Soot2D_State() : MW_gas(NULL) {
    
    // setup model
    setModelParams(SOOT2D_LEUNG_ET_AL_1991);

    // setup a N2 gas mixture
    double WT[] = {MOLE_WT_N2};
    setGasPhase( WT, 1);

  };

  //! Destructor
  ~Soot2D_State() { Deallocate(); };

  /****************** Allocators / Deallocators *********************/
  void Allocate();
  void Deallocate();

  /************************* Accessors ******************************/
  //! return number of variables
  int NumVar() { return nv; }

  /******************* Miscillaneous Functions **********************/
  //! mean soot diameter
  double meanDiameter( const double*sc ) const;

  //! set model params
  void setModelParams(const int model_flag);

  //! Set gas mixture
  void setGasPhase( const double* MWgas, const int nsp, const string* names=NULL);
  

  /*********************** Source Terms *****************************/
  //! Computation of source terms
  int getRates( const double T, 
		const double rho_gas,
		const double* y_gas,
		const double*sc,
		double* ydot_gas,
		double* ydot_soot,
		const double& mult=1.0 ) const;

  /******************** Operator Overloads **************************/
  //! Input-output operators.
  friend ostream& operator << (ostream &out_file, const Soot2D_State &U);


};

/////////////////////////////////////////////////////////////////////
/// INLINE FUNCTION DEFINITIONS
/////////////////////////////////////////////////////////////////////

/*********************************************************************
 * Allocate / Deallocate dynamic memory
 *********************************************************************/
inline void Soot2D_State::Allocate() {
  Deallocate();
  if (ns_gas>0) { MW_gas = new double[ns_gas]; }
}
inline void Soot2D_State::Deallocate() {
  if (MW_gas != NULL) { delete[] MW_gas; MW_gas = NULL; }
}

/*********************************************************************
 * Set relevant model parameters
 *********************************************************************/
inline void Soot2D_State::setModelParams(const int model_flag) {
    
  model = model_flag;

  //-----------------------------------------------------------------
  // Soot Model parameters from Leung et al. (1991)
  //-----------------------------------------------------------------
  if (model == SOOT2D_LEUNG_ET_AL_1991) {
    nv = 2;        // two parameters (mass frac, number density)
    MW = 12.011;   // [kg/kmol]
    rho = 2000.0;  // [kg/m^3]]
    Cmin = 100;
    Ca = 9.0;
    
  //-----------------------------------------------------------------
  // No model - do nothing
  //-----------------------------------------------------------------
  } else if (model == SOOT2D_NONE) {
    
  //-----------------------------------------------------------------
  // Invalid model
  //-----------------------------------------------------------------
  } else {
      cerr << "\n Error in Soot2D_State::setModelParams() - "
	   << "Invalid model flag.\n";
      exit(-1);
  }
}

/*********************************************************************
 * Set relevant model parameters
 *********************************************************************/
inline void Soot2D_State::setGasPhase( const double* MWgas, 
				       const int nsp, 
				       const string* names) {

  // allocate memory
  ns_gas = nsp;
  Allocate();

  // copy molar weights
  for ( int i=0; i<nsp; i++) MW_gas[i] = MWgas[i];

  // determine array indices
  iC2H2 = iO2 = iH2 = iCO = -1;
  if (names != NULL) {

    for ( int i=0; i<nsp; i++) {
      if (names[i] == "C2H2") iC2H2 = i;
      else if (names[i] == "O2") iO2 = i;
      else if (names[i] == "H2") iH2 = i;
      else if (names[i] == "CO") iCO = i;      
    }

    if ( iC2H2==-1 || iO2==-1 || iH2==-1 || iCO==-1) {
      cerr << "\n Error in Soot2D_State::setGasPhase() - "
	   << "Could not find C2H2, O2, H2, CO in species list.\n";
      exit(-1);
    }

  }
}


#endif // _SOOT2D_STATE_INCLUDED
