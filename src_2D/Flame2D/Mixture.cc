/////////////////////////////////////////////////////////////////////
///
/// \file Mixture.cc
/// 
/// \author Marc R.J. Charest
/// 
/// \brief This header file contains the class definition for 
///        describing the thermodynamic and transport properties
///        of a ideal gas mixture.
///
/////////////////////////////////////////////////////////////////////
#include "Mixture.h"

/**
 * Initialization of static variables.
 */
#ifndef STATIC_NUMBER_OF_SPECIES
int     Mixture :: ns = 0;
#endif
double  Mixture :: Tmin = 0.0;
double  Mixture :: Tmax = 0.0;
double* Mixture :: Sc_ref = NULL;
double* Mixture :: M = NULL;
double* Mixture :: Hform = NULL;
string* Mixture :: names = NULL;
double* Mixture :: r = NULL;
double* Mixture :: r0 = NULL;
bool    Mixture :: isConstSchmidt = false;
string  Mixture :: ct_mech_name = "gri30";
string  Mixture :: ct_mech_file = "gri30.xml";
IdealGasMix* Mixture :: ct_gas = NULL;
Transport*   Mixture :: ct_trans = NULL;

/////////////////////////////////////////////////////////////////////
/// STATIC MEMBER FUNCTIONS
/////////////////////////////////////////////////////////////////////

/****************************************************
 * Main setup of mixture properties.  Load the specific gas data for
 * each species, and set some static flags.
 *
 * \param n     Number of species.
 * \param names Array of species names.
 * \param PATH  Path to data files.
 * \param Schmidt Array of species schmidt numbers.
 * \param trans_data_flag Flag indicating whether NASA of Lennard-Jones data
 * \param constant_schmidt Boolean indicating whether constant schmidt values.
 ****************************************************/
void Mixture :: setMixture(const string &mech_name,
			   const string &mech_file) {
  
  // deallocate first
  DeallocateStatic();

  //create a new ideal gas mixture class
  ct_mech_name = mech_name;
  ct_mech_file = mech_file;
  try {
    ct_gas = new IdealGasMix(ct_mech_file, ct_mech_name);
    ct_trans = newTransportMgr("Mix", ct_gas);
  }
  catch (CanteraError) {
    Cantera::showErrors();
  }

  // ascertain the number of species
#ifdef STATIC_NUMBER_OF_SPECIES
  if( STATIC_NUMBER_OF_SPECIES != ct_gas->nSpecies()) {
      cerr << "\n ERROR, Mixture::setMixture() - Built using static species with "
	   << STATIC_NUMBER_OF_SPECIES 
	   << " species predefined, asking for " << ns
	   << endl; 
      exit(-1); 
    }
#else
  ns = ct_gas->nSpecies();
#endif

  //set data temperature ranges for mixture
  Tmin = ct_gas->minTemp();
  Tmax = ct_gas->maxTemp();
  
  //allocate static memory  
  AllocateStatic();

  // get non-dimensional heats of formation
  ct_gas->setTemperature(TREF);
  ct_gas->setPressure(PREF);
  ct_gas->getEnthalpy_RT( Hform );

  // load the species data
  for(int i=0; i<ns; i++){
    Sc_ref[i] = 1.0;
    M[i] = ct_gas->molarMass(i);
    names[i] = ct_gas->speciesName(i);
    Hform[i] *= (Cantera::GasConstant/M[i]) * TREF;
  }
  
}   

/****************************************************
 * Set a constant schmidt number
 ****************************************************/
void Mixture :: setConstantSchmidt(const double* Sc) {
  isConstSchmidt = true;
  for(int i=0; i<ns; i++){
    Sc_ref[i] = Sc[i];
  }  
}
/****************************************************
 * Mixture molecular mass [kg/mol]
 ****************************************************/
double Mixture :: molarMass( const double* y ) {
  double sum( 0.0 );
  for(int i=0; i<ns; i++){
    sum += y[i]/M[i];
  }
  return 1.0/sum;
}

/****************************************************
 * Mixture heat of formation [J/kg]
 ****************************************************/
double Mixture :: heatFormation( const double* y ) {
  double sum( 0.0 );
  for(int i=0; i<ns; i++){
    sum += y[i]*Hform[i];
  }
  return sum;
}

// /****************************************************
//  * Mixture gas constant [kg/mol]
//  ****************************************************/
// double Mixture :: gasConstant( const double* y ){
//   return Cantera::GasConstant/molarMass(y);
// }

// /****************************************************
//  * Mixture Heat Capacity (const pressure) J/(kg*K)
//  ****************************************************/
// double Mixture :: heatCapacity_p( const double &Temp, const double* y ){
//   ct_gas->setMassFractions_NoNorm(y);
//   ct_gas->setTemperature(Temp);
//   return ct_gas->cp_mass();
// }

// /****************************************************
//  * Mixture Heat Capacity (const volume) J/(kg*K)
//  ****************************************************/
// double Mixture :: heatCapacity_v( const double &Temp, const double* y ){
//   ct_gas->setMassFractions_NoNorm(y);
//   ct_gas->setTemperature(Temp);
//   return ct_gas->cv_mass();
// }

// /****************************************************
//  * Mixture Heat Ratio gamma
//  ****************************************************/
// double Mixture :: heatRatio( const double &Temp, const double* y ){
//   ct_gas->setMassFractions_NoNorm(y);
//   ct_gas->setTemperature(Temp);
//   return ct_gas->cp_mole()/ct_gas->cv_mole();
// }

// /****************************************************
//  * Mixture Specific Internal Energy J/(kg)
//  ****************************************************/
// //! etotal = sensible & chemical
// double Mixture :: internalEnergy( const double &Temp, const double* y ){
//   ct_gas->setMassFractions_NoNorm(y);
//   ct_gas->setTemperature(Temp);
//   return ct_gas->intEnergy_mass();
// }

// //! internal energy with no heat of formation included (sensible)
// double Mixture :: internalEnergySens( const double &Temp, const double* y ){
//   ct_gas->setMassFractions_NoNorm(y);
//   ct_gas->setTemperature(Temp);
//   return (ct_gas->intEnergy_mass() - heatFormation(y));
// }

// /****************************************************
//  * Specific absolute enthalpy J/(kg)
//  ****************************************************/
// //! htotal = mass fraction * (hsensible + heatofform)
// double Mixture :: enthalpy( const double &Temp, const double* y ){
//   ct_gas->setMassFractions_NoNorm(y);
//   ct_gas->setTemperature(Temp);
//   return ct_gas->enthalpy_mass();
// }

// //! absolute enthalpy with no heat of formation included (sensible)
// double Mixture :: enthalpySens( const double &Temp, const double* y ) {
//   ct_gas->setMassFractions_NoNorm(y);
//   ct_gas->setTemperature(Temp);
//   return (ct_gas->enthalpy_mass() - heatFormation(y));
// }

// /****************************************************
//  * Derivative of specific enthalpy dh/dT
//  * actually is just Cp as Cp = (dh/dT)_p
//  ****************************************************/
// double Mixture :: enthalpyPrime( const double &Temp, const double* y ){
//   return heatCapacity_p( Temp, y );
// }

// /****************************************************
//  * Viscosity using Wilke [1950] formulation
//  ****************************************************/
// double Mixture :: viscosity(const double &Temp, const double* y) {
//   ct_gas->setMassFractions_NoNorm(y);
//   ct_gas->setTemperature(Temp);
//   return ct_trans->viscosity();
// }

// /****************************************************
//  * Thermal Conductivity - Mason & Saxena (1958)  W/(m*K)
//  ****************************************************/
// double Mixture :: thermalCond(const double &Temp, const double* y) {
//   ct_gas->setMassFractions_NoNorm(y);
//   ct_gas->setTemperature(Temp);
//   return ct_trans->thermalConductivity();
// }

// /****************************************************
//  * Sepcies mixture diffusion coefficient (not implemented)
//  ****************************************************/
// double Mixture :: speciesDiffCoeff(const double &Temp, 
// 				   const double* y, 
// 				   const int &i){
//   cerr << "Mixture::speciesDiffCoeff() - Not implemented yet.";
//   exit(-1);
// }

// /****************************************************
//  * Schmidt
//  ****************************************************/
// double Mixture :: schmidt(const double &Temp, 
// 			  const double &rho,
// 			  const double* y, 
// 			  const int &i) {
//   if(isConstSchmidt){
//     return Sc_ref[i];
//   } else {
//     return viscosity(Temp, y)/(rho*speciesDiffCoeff(Temp,y,i));
//   }
// }

// /****************************************************
//  * Prandtl
//  ****************************************************/
// double Mixture :: prandtl(const double &Temp, const double* y) {
//   //Pr = Cp*mu/k
//   return heatCapacity_p(Temp,y)*viscosity(Temp,y)/thermalCond(Temp,y);
// }

// /****************************************************
//  * Lewis
//  ****************************************************/
// double Mixture :: lewis(const double &Temp, 
// 			const double &rho,
// 			const double* y,
// 			const int &i) {
//    if(isConstSchmidt)
//      return (thermalCond(Temp,y)*Sc_ref[i])/
//        ( heatCapacity_p(Temp,y)*viscosity(Temp,y) );
//    else
//      return thermalCond(Temp,y) / 
//        ( heatCapacity_p(Temp,y)*rho*speciesDiffCoeff(Temp,y,i) );
// }


// /****************************************************
//  * Derivative of species h wrt to mass fraction
//  ****************************************************/
// void Mixture :: getDihdDc(const double &Temp, 
// 			  const double &Press, 
// 			  const double* y, 
// 			  double* dh) {
//   ct_gas->setState_TPY(Temp, Press, y);
//   ct_gas->getEnthalpy_RT(dh); // -> h = hs + hf
//   double cp( ct_gas->cp_mass() );
//   double M_mix( molarMass(y) );
//   for(int i=0; i<ns; i++)
//     dh[i] = ( dh[i]*(Cantera::GasConstant/M[i])*Temp -
// 	      cp*Temp*M_mix/M[i] );
// }


/***********************************************************************
  Use cantera to parse the input mass fraction string of the form
      CH4:0.5, O2:0.5
  All other species will be assumed to have 0 mass fractions.  Cantera
  also normalizes the mass fractions to sum to unity.  Returns them
  in an array.
***********************************************************************/
void Mixture::parse_mass_string(const string& massFracStr, 
				double* massFracs) {

  //
  // for a string of 'space' delimited mass fractions
  //
  if ( massFracStr.find(":",0) == string::npos ) {
    
    // get first position
    string delimiters = " ";
    string::size_type lastPos = massFracStr.find_first_not_of(delimiters, 0);
    string::size_type pos     = massFracStr.find_first_of(delimiters, lastPos);

    // parse the string and set mass fractions
    int index = 0;
    while (string::npos != pos || string::npos != lastPos) {
      massFracs[index] = atof(massFracStr.substr(lastPos, pos - lastPos).c_str());
      massFracs[index] = max(massFracs[index], 0.0);
      lastPos = massFracStr.find_first_not_of(delimiters, pos);
      pos = massFracStr.find_first_of(delimiters, lastPos);
      index++;
    } // endwhile

    // normalize to unity
    double sum(0.0);
    for(int index=0; index<ns; index++) sum += massFracs[index];
    if (sum>1.e-7) for(int index=0; index<ns; index++) massFracs[index] /= sum;

  //
  // for a composition map (eg. CH4:0.5, O2:0.5)
  //
  } else {

    // set mass fractions and make sure they sum to unity
    ct_gas->setMassFractionsByName(massFracStr);
    for(int index =0; index<ns; index++)
      massFracs[index] = ct_gas->massFraction(index);

  } // endif

} // end of ct_parse_mass_string

/***********************************************************************
  Use cantera to parse the input mole fraction string of the form
      CH4:0.5, O2:0.5
  All other species will be assumed to have 0 mole fractions.  Cantera
  also normalizes the mole fractions to sum to unity.  Returns them in
  an array.
***********************************************************************/
void Mixture::parse_mole_string(const string& moleFracStr, 
				double* moleFracs) {

  //
  // for a string of 'space' delimited mole fractions
  //
  if ( moleFracStr.find(":",0) == string::npos ) {
    
    // get first position
    string delimiters = " ";
    string::size_type lastPos = moleFracStr.find_first_not_of(delimiters, 0);
    string::size_type pos     = moleFracStr.find_first_of(delimiters, lastPos);

    // parse the string and set molar fractions
    int index = 0;
    while (string::npos != pos || string::npos != lastPos) {
      moleFracs[index] = atof(moleFracStr.substr(lastPos, pos - lastPos).c_str());
      moleFracs[index] = max(moleFracs[index], 0.0);
      lastPos = moleFracStr.find_first_not_of(delimiters, pos);
      pos = moleFracStr.find_first_of(delimiters, lastPos);
      index++;
    } // endwhile

    // normalize to unity
    double sum(0.0);
    for(int index=0; index<ns; index++) sum += moleFracs[index];
    if (sum>1.e-7) for(int index=0; index<ns; index++) moleFracs[index] /= sum;

  //
  // for a composition map (eg. CH4:0.5, O2:0.5)
  //
  } else {

    // set mole fractions and make sure they sum to unity
    ct_gas->setMoleFractionsByName(moleFracStr);
    for(int index =0; index<ns; index++)
      moleFracs[index] = ct_gas->moleFraction(index);

  } // endif

} // end of ct_parse_mass_string


/***********************************************************************
  Use cantera to parse the input schmidt number string of the form
      CH4:0.5, O2:0.5
  All other species will be assumed to have unity Schmidt number.  
  Returns them in an array.
***********************************************************************/
void Mixture::parse_schmidt_string( const string& schmidtStr, 
				    double* schmidt) {

  //
  // for a string of 'space' delimited schmidt numbers
  //
  if ( schmidtStr.find(":",0) == string::npos ) {
    
    // get first position
    string delimiters = " ";
    string::size_type lastPos = schmidtStr.find_first_not_of(delimiters, 0);
    string::size_type pos     = schmidtStr.find_first_of(delimiters, lastPos);

    // parse the string and set schmidt numbers
    int index = 0;
    while (string::npos != pos || string::npos != lastPos) {
      schmidt[index] = atof(schmidtStr.substr(lastPos, pos - lastPos).c_str());
      if(schmidt[index]<0.0) schmidt[index] = 1.0;
      lastPos = schmidtStr.find_first_not_of(delimiters, pos);
      pos = schmidtStr.find_first_of(delimiters, lastPos);
      index++;
    } // endwhile

  //
  // for a composition map (eg. CH4:0.5, O2:0.5)
  //
  } else {

    // declares
    compositionMap xx;
    int kk = ct_gas->nSpecies();
    double s;
    
    // build composition map and initialize
    for (int k = 0; k < kk; k++) xx[ct_gas->speciesName(k)] = -1.0;
    
    // parse map
    parseCompString(schmidtStr, xx);
    
    // set schmidt numbers
    for (int k = 0; k < kk; k++) { 
      s = xx[ct_gas->speciesName(k)];
      if (s > 0.0) schmidt[k] = s;
      else schmidt[k] = 1.0;
    }
    
  } // endif

} // end of ct_parse_mass_string

/***********************************************************************
  Use cantera to return the species index.  Exits in error if an
  unidentified species is requested.
***********************************************************************/
int Mixture::speciesIndex(const string &sp) {

  // get index
  int index = ct_gas->speciesIndex(sp);

  // if index==-1, species not found
  if (index<0) {
    cerr<<"\n Reaction_set::ct_get_species_index() - Index of unkown species, "
	<<sp<<", requested.\n";
    exit(-1);
  } // endif

  // return index
  return index;


} // end of ct_get_species_index



/************************************************************************
  Calculates the equilibrium composition given an unburnt mixture
  using CANTERA.  This is only for CANTERA reaction types.  Here
  we hold enthalpy and pressure fixed.

  Wu - unburnt state
  Wb - burnt state

************************************************************************/

void Mixture::equilibrate_HP( double &Temp, double &Press, double* y ) {

  // set state and equilibrate
  ct_gas->setState_TPY(Temp, Press, y);
  equilibrate( *ct_gas, "HP" );

  //get burnt mass fractions
  ct_gas->getMassFractions(y);

  // the temperature
  Temp = ct_gas->temperature();
  
} // end of ct_equilibrate


/************************************************************************
  Calculates the equilibrium composition given an unburnt mixture
  using CANTERA.  This is only for CANTERA reaction types.  Here 
  we hold temperature and pressure fixed.

  Wu - unburnt state
  Wb - burnt state

************************************************************************/
void Mixture::equilibrate_TP( double &Temp, double &Press, double* y ) {

  // set state and equilibrate
  ct_gas->setState_TPY(Temp, Press, y);
  equilibrate( *ct_gas, "TP" );

  //get burnt mass fractions
  ct_gas->getMassFractions(y);

} // end of ct_equilibrate



/***********************************************************************

  This function computes the composition of a hydrocarbon (CxHy), O2, 
  and N2 mixture.  Cantera is used to convert molar fraction to mass
  fraction.

***********************************************************************/
void Mixture::composition( const string& fuel_species, 
			   const double &phi,
			   double* massFracs ) {

  // first, compute the stoichiometric fuel air ratio
  double C_atoms=ct_gas->nAtoms(ct_gas->speciesIndex(fuel_species), ct_gas->elementIndex("C"));
  double H_atoms=ct_gas->nAtoms(ct_gas->speciesIndex(fuel_species), ct_gas->elementIndex("H"));
  double ax=C_atoms+H_atoms/4.0;
  double fa_stoic=1.0/(4.76*ax);

  // determine the composition
  int nsp = ct_gas->nSpecies();  
  double* x = new double[nsp];
  for(int k=0; k<nsp; k++){
    if(k==ct_gas->speciesIndex(fuel_species)){ x[k]=phi; }
    else if(k==ct_gas->speciesIndex("O2")){    x[k]=1.00*ax; }
    else if(k==ct_gas->speciesIndex("N2")){    x[k]=3.76*ax; }
    else{ x[k]=0.0; }
  }

  // compute composition 
  // -> why do it yourself when you can get someone else to do it
  ct_gas->setMoleFractions(x);
  for(int k=0;k<nsp;k++) massFracs[k] = ct_gas->massFraction(k);

  // clean up memory
  delete[] x;

} // end of ct_composition


// /************************************************************************
//   Calculates the concentration time rate of change of species from
//   primitive state W using the general law of mass action.
//   U is the conserved state container for passing back the 
//   source terms. ie. U.rhospec[i].c 

//   W.SpecCon:  is the  species mass fractions concentrations
//               of Chem2D_pState. (c_i*rho/M_i)   mol/m^3

//   Return units are  kg/m^3*s ie. rho*omega (kg/m^3)*(1/s)

// ************************************************************************/
// void Mixture :: getRates( const double &Temp, const double &Press, 
// 			  const double* y, double* rr ) {
//   ct_gas->setState_TPY(Temp, Press, y);
//   ct_gas->getNetProductionRates(rr);
// }


/////////////////////////////////////////////////////////////////////
/// LOCAL MEMBER FUNCTIONS
/////////////////////////////////////////////////////////////////////
