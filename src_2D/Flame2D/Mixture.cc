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
int     Mixture :: ns = 0;
double  Mixture :: Tmin = 0.0;
double  Mixture :: Tmax = 0.0;
double* Mixture :: Sc_ref = NULL;
bool    Mixture :: isConstSchmidt = true;
NASARP1311data* Mixture :: specdata = NULL;

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
void Mixture :: setMixture(const int &n, 
			   const string *names,
			   const char *PATH,
			   const double* Schmidt, 
			   const int &trans_data_flag,
			   const bool &constant_schmidt) {
  
  // set new number of species
#ifdef STATIC_NUMBER_OF_SPECIES
  if( STATIC_NUMBER_OF_SPECIES < n) {
      cerr << "\n ERROR, Mixture::setMixture() - Built using static species with "
	   << STATIC_NUMBER_OF_SPECIES 
	   << " species predefined, asking for " << n 
	   << endl; 
      exit(1); 
    }
#endif
  ns = n;

  //allocate static memory and load the species data  
  AllocateStatic();
  for(int i=0; i<ns; i++){
    specdata[i].Getdata(names[i],PATH,trans_data_flag);  
    Sc_ref[i] = Schmidt[i];
  }

  // are we using constant schmidt
  isConstSchmidt = constant_schmidt;
  
  //set data temperature ranges for mixture
  Tmin = specdata[0].Low_range();
  Tmax = specdata[0].High_range();
  for(int i=1; i<ns; i++){
    Tmin = max(specdata[i].Low_range(), Tmin);
    Tmax = min(specdata[i].High_range(), Tmax);
  }
  
}   

/****************************************************
 * Mixture molecular mass [kg/mol]
 ****************************************************/
double Mixture :: molarMass( const double* y ) {
  double sum( 0.0 );
  for(int i=0; i<ns; i++) sum += y[i]/specdata[i].Mol_mass();
  return 1.0/sum;
}

/****************************************************
 * Mixture gas constant [kg/mol]
 ****************************************************/
double Mixture :: gasConstant( const double* y ){
  // = sum ( mass fraction * species gas constant)
  double sum ( 0.0 );
  for(int i=0; i<ns; i++) sum += y[i] * specdata[i].Rs();
  return sum;
}

/****************************************************
 * Mixture Heat Capacity (const pressure) J/(kg*K)
 ****************************************************/
double Mixture :: heatCapacity_p( const double &Temp, const double* y ){
  double sum( 0.0 );
  for(int i=0; i<ns; i++) sum += y[i]*specdata[i].HeatCapacity_p(Temp);
  return sum;
}

/****************************************************
 * Mixture Heat Capacity (const volume) J/(kg*K)
 ****************************************************/
double Mixture :: heatCapacity_v( const double &Temp, const double* y ){
  double sum( 0.0 );
  for(int i=0; i<ns; i++) sum += y[i]*specdata[i].HeatCapacity_v(Temp);
  return sum;
}

/****************************************************
 * Mixture Heat Ratio gamma
 ****************************************************/
double Mixture :: heatRatio( const double &Temp, const double* y ){
  // = Cp / Cv  
  return heatCapacity_p(Temp, y)/heatCapacity_v(Temp, y);
}

/****************************************************
 * polytropic heat ratio mixture gamma assuming 
 * T=200K as the temperature.
 ****************************************************/
double Mixture :: heatRatioRef( const double* y ){
  double sum1(ZERO); 
  double sum2(ZERO);
  double gamma_s(ZERO);

  for(int i=0; i<ns; i++){
    sum1 += y[i]*specdata[i].Rs();
    gamma_s = ONE/(ONE - specdata[i].Rs()/ specdata[i].HeatCapacity_p(REFERENCE_TEMPERATURE));
    sum2 += ((y[i]*specdata[i].Rs()) / (gamma_s - ONE)); 
  }
  return (ONE + sum1/sum2);
}

/****************************************************
 * Mixture Specific Internal Energy J/(kg)
 ****************************************************/
//! etotal = sensible & chemical
double Mixture :: internalEnergy( const double &Temp, const double* y ){
  // = sum (mass fraction * species e) 
  double sum( 0.0 );
  for(int i=0; i<ns; i++) {
    //(Enthalpy(Temp) - (R/mol_mass)*Temp)
    sum += y[i]*(specdata[i].Enthalpy(Temp) + 
		 specdata[i].Heatofform() -
		 specdata[i].Rs()*Temp);
  }
  return sum;
}

//! reference internal energy e + heat of formation - offset
double Mixture :: internalEnergyRef( const double &Temp, const double* y ){
  // = sum (mass fraction * species e) 
  double sum( 0.0 );
  for(int i=0; i<ns; i++){ 
    //(Enthalpy(Temp) - (R/mol_mass)*Temp)
    sum += y[i]*(specdata[i].Enthalpy(Temp) + 
		 specdata[i].Heatofform() - 
		 specdata[i].DeltaHref() - 
		 specdata[i].Rs()*Temp);
  }
  return sum;
}

//! internal energy with no heat of formation included (sensible)
double Mixture :: internalEnergySens( const double &Temp, const double* y ){
  // = sum (mass fraction * species e) 
  double sum( 0.0 );
  for(int i=0; i<ns; i++){ 
    //(Enthalpy(Temp) - (R/mol_mass)*Temp)
    sum += y[i]*(specdata[i].Enthalpy(Temp) - specdata[i].Rs()*Temp);
  }
  return sum;
}

/****************************************************
 * Specific absolute enthalpy J/(kg)
 ****************************************************/
//! htotal = mass fraction * (hsensible + heatofform)
double Mixture :: enthalpy( const double &Temp, const double* y ){
  // = sum (mass fraction * species h) 
  double sum( 0.0 );
  for(int i=0; i<ns; i++){
    sum += y[i]*(specdata[i].Enthalpy(Temp) + 
		 specdata[i].Heatofform());
  }
  return sum;
}

//! reference absolute enthalpy h - offset
double Mixture :: enthalpyRef( const double &Temp, const double* y ){
  // = sum (mass fraction * species h) 
  double sum( 0.0 );
  for(int i=0; i<ns; i++){ 
    sum += y[i]*(specdata[i].Enthalpy(Temp) + 
		 specdata[i].Heatofform() - 
		 specdata[i].DeltaHref() );
  }
  return sum;
}

//! absolute enthalpy with no heat of formation included (sensible)
double Mixture :: enthalpySens( const double &Temp, const double* y ) {
  // = sum (mass fraction * species h) 
  double sum( 0.0 );
  for(int i=0; i<ns; i++) sum += y[i]*(specdata[i].Enthalpy(Temp));
  return sum;
}

/****************************************************
 * Derivative of specific enthalpy dh/dT
 * actually is just Cp as Cp = (dh/dT)_p
 ****************************************************/
double Mixture :: enthalpyPrime( const double &Temp, const double* y ){
  return heatCapacity_p( Temp, y );
}

/****************************************************
 * Viscosity using Wilke [1950] formulation
 ****************************************************/
double Mixture :: viscosity(const double &Temp, const double* y) {
  double sum(0.0);
  double phi;
#ifdef STATIC_NUMBER_OF_SPECIES
  double  vis[ns];
#else 
  double* vis = new double[ns];
#endif

  for(int i=0; i<ns; i++){
    phi = 0.0;
    for (int j=0; j<ns; j++){
      if(i == 0) vis[j] = specdata[j].Viscosity(Temp);
      phi += (y[j] / specdata[j].Mol_mass())*
	sqr( ONE + 
	     sqrt(vis[i]/vis[j])*pow(specdata[j].Mol_mass()/
				     specdata[i].Mol_mass(),0.25) ) /
	sqrt(EIGHT*(ONE +specdata[i].Mol_mass()/specdata[j].Mol_mass()));
    }
    sum += (y[i] * vis[i]) / (specdata[i].Mol_mass() * phi);
  }  

#ifndef STATIC_NUMBER_OF_SPECIES
  delete[] vis;
#endif

  return sum;

}

/****************************************************
 * Molecular Viscosity using Wilke [1950] formulation
 * as a function of temperature.
 ****************************************************/
double Mixture :: viscosityPrime(const double &Temp, const double* y) {
  double sum(0.0);
  double phi;
#ifdef STATIC_NUMBER_OF_SPECIES
  double  vis[ns];
#else 
  double* vis = new double[ns];
#endif

  for(int i=0; i<ns; i++){
    phi = 0.0;
    for (int j=0; j<ns; j++){
      if(i == 0) vis[j] = specdata[j].dViscositydT(Temp);
      phi += (y[j] / specdata[j].Mol_mass())*
	sqr(1.0 + 
	    sqrt(vis[i]/vis[j])*pow(specdata[j].Mol_mass()/
				    specdata[i].Mol_mass(),0.25) ) /
	sqrt(8.0*(1.0 +specdata[i].Mol_mass()/specdata[j].Mol_mass()));
    }
    sum += (y[i] * vis[i] ) / (specdata[i].Mol_mass() * phi);
  }  


#ifndef STATIC_NUMBER_OF_SPECIES
  delete[] vis;
#endif

  return sum;

}

/****************************************************
 * Thermal Conductivity - Mason & Saxena (1958)  W/(m*K)
 ****************************************************/
double Mixture :: thermalCond(const double &Temp, const double* y) {
  double sum(0.0);
  double phi;
#ifdef STATIC_NUMBER_OF_SPECIES
  double  vis[ns];
#else 
  double* vis = new double[ns];
#endif

  for(int i=0; i<ns; i++){
    phi = 0.0;
    for (int j=0; j<ns; j++){
      if(i == 0) vis[j] = specdata[j].Viscosity(Temp);
      if(i != j){
	phi += (y[j] / specdata[j].Mol_mass())*
	  sqr(ONE + 
	      sqrt(vis[i]/vis[j])*pow(specdata[j].Mol_mass()/
				      specdata[i].Mol_mass(),0.25)) /
	  sqrt(EIGHT*(ONE +specdata[i].Mol_mass()/specdata[j].Mol_mass()));
      }
    }

    sum += (specdata[i].ThermalConduct(Temp)*y[i]) / 
      (y[i] + (specdata[i].Mol_mass()) * 1.065 * phi);
  }  


#ifndef STATIC_NUMBER_OF_SPECIES
  delete[] vis;
#endif

  return sum;

}

/****************************************************
 * Sepcies mixture diffusion coefficient (not implemented)
 ****************************************************/
double Mixture :: speciesDiffCoeff(const double &Temp, 
				   const double* y, 
				   const int &i){
  cerr << "Mixture::speciesDiffCoeff() - Not implemented yet.";
  exit(-1);
}

/****************************************************
 * Schmidt
 ****************************************************/
double Mixture :: schmidt(const double &Temp, 
			  const double &rho,
			  const double* y, 
			  const int &i) {
  if(isConstSchmidt){
    return Sc_ref[i];
  } else {
    return viscosity(Temp, y)/(rho*speciesDiffCoeff(Temp,y,i));
  }
}

/****************************************************
 * Prandtl
 ****************************************************/
double Mixture :: prandtl(const double &Temp, const double* y) {
  //Pr = Cp*mu/k
  return heatCapacity_p(Temp,y)*viscosity(Temp,y)/thermalCond(Temp,y);
}

/****************************************************
 * Lewis
 ****************************************************/
double Mixture :: lewis(const double &Temp, 
			const double &rho,
			const double* y,
			const int &i) {
   if(isConstSchmidt)
     return (thermalCond(Temp,y)*Sc_ref[i])/
       ( heatCapacity_p(Temp,y)*viscosity(Temp,y) );
   else
     return thermalCond(Temp,y) / 
       ( heatCapacity_p(Temp,y)*rho*speciesDiffCoeff(Temp,y,i) );
}


/****************************************************
 * polytropic heat ratio mixture gamma assuming 
 * T=200K as the temperature.
 ****************************************************/
double Mixture :: gammaGuess(const double* y) {  
  double sum1(ZERO), sum2(ZERO);
  double gamma_s(ZERO), Temp(200.0);

  for(int i=0; i<ns; i++){
    sum1 += y[i]*specdata[i].Rs();
    gamma_s = ONE/(ONE - specdata[i].Rs()/ specdata[i].HeatCapacity_p(Temp));
    sum2 += ((y[i]*specdata[i].Rs()) / (gamma_s - ONE)); 
  }
  return (ONE + sum1/sum2);
}

/**************************************************
 * Temperature derived from specific sensible
 * enthalpy
 **************************************************/
double Mixture :: temperature(double &h_s, const double* y) {

  // declares
  double RTOT(gasConstant(y));

  // set iteration parameters
  static const int Nmax = 20;

  // determine limits
#ifdef TLOWERBOUNDS
  double T_min(TLOWERBOUNDS);
#else
  double T_min(Tmin);
#endif
  double T_max(Tmax);
  //double fmin(hs(T_min) - h_s);
  //double fmax(hs(T_max) - h_s);


  //------------------------------------------------
  // Initial Guess
  //------------------------------------------------
  //using a polytropic gas assumption with gamma@200;
  double Tguess( (h_s/RTOT)*(ONE/(ONE/(gammaGuess(y) - ONE) +ONE)) );

  //check for start value
  double Tn;
  if (Tguess < T_min )
    Tn=T_min;
  else if (Tguess > T_max)
    Tn=T_max;
  else
    Tn=Tguess;

  // compute function values
  double fn( enthalpySens(Tn,y) - h_s );
  double dfn( enthalpyPrime(Tn,y) );
  double dTn( fabs(T_max - T_min) );
  double dT0( dTn );

  // No need to orient search such that f(T_min)<0,
  // already oriented.


  //------------------------------------------------
  // Newton-Raphson iteration
  //------------------------------------------------
  int i;
  for ( i=1; i<=Nmax; i++){

    // use bisection if Newton out of range or not decreasing 
    // fast enough
    if ( (((Tn-T_max)*dfn-fn)*((Tn-T_min)*dfn-fn) >= ZERO) || 
	 (fabs(TWO*fn) > fabs(dT0*dfn)) ) {
      dT0 = dTn;
      dTn = HALF*(T_max-T_min);
      Tn = T_min + dTn;

    // Newton acceptable 
    } else {
      dT0 = dTn;
      dTn = fn/dfn;
      Tn -= dTn;
    }

    // evaluate new guess
    fn = enthalpySens(Tn,y) - h_s;  
    dfn = enthalpyPrime(Tn,y); 

    // Convergence test
    if ( fabs(dTn)<CONV_TOLERANCE || fabs(fn)<CONV_TOLERANCE ) break;

    // change bisection bracket
    if ( fn < ZERO)  T_min=Tn;
    else T_max = Tn;

  } // end Newton-Raphson


  if (i>Nmax){
    cout<<"\nTemperature didn't converge in Flame2D_pState::T(double &h_s)";
    cout<<" with polytopic Tguess "<<Tguess<<" using "<<Tn;
  }

  //return value
  return Tn;
}

/****************************************************
 * GIBBS Free Energy.
 *
 * Eqn. (10.84) (10.87) Anderson
 * Gs = Hs - TS
 * Gs(ps=1) = Gs - R_UNIVERSAL*T*ln(ps)  //for data not at 1atm
 * ps = cs(M/Ms)*p
 ****************************************************/
double Mixture::speciesGibbsFree(const double &Temp, 
				 const int &species) {
  return specdata[species].Enthalpy_mol(Temp) - 
    Temp*specdata[species].Entropy_mol(Temp);
}


/////////////////////////////////////////////////////////////////////
/// MEMBER FUNCTIONS
/////////////////////////////////////////////////////////////////////


/****************************************************
 * Set the state of the mixture and update all the properties.
 *
 * \param Temp Mixture temperature [k]
 * \param y Mixture species mass fractions.
 ****************************************************/
void Mixture :: setState_TY(const double &Temp, const double* y)
{

  // declares
  double Cp_ref(0.0);
  double phi;
#ifdef STATIC_NUMBER_OF_SPECIES
  double  vis[ns];
#else 
  double* vis = new double[ns];
#endif


  // initialize
  T = Temp;
  MW = 0.0;  hs = 0.0;
  hf = 0.0;  Cp = 0.0;
  mu = 0.0;  kappa = 0.0;
  g_ref = 0.0;

  //
  // Main loop over species
  //
  for(int i=0; i<ns; i++){

    // thermodynamic coefficients
    MW += y[i]/specdata[i].Mol_mass();
    Cp +=  y[i]*specdata[i].HeatCapacity_p(Temp);
    hs += y[i]*specdata[i].Enthalpy(Temp);
    hf += y[i]*specdata[i].Heatofform();
    Cp_ref +=  y[i]*specdata[i].HeatCapacity_p(REFERENCE_TEMPERATURE);

    // compute mixture rule
    phi = 0.0;
    for (int j=0; j<ns; j++){
      if(i == 0) vis[j] = specdata[j].Viscosity(Temp);
      phi += (y[j] / specdata[j].Mol_mass())*
	sqr( ONE + 
	     sqrt(vis[i]/vis[j])*pow(specdata[j].Mol_mass()/
				     specdata[i].Mol_mass(),0.25) ) /
	sqrt(EIGHT*(ONE +specdata[i].Mol_mass()/specdata[j].Mol_mass()));
    }

    // transport coeffcieints
    mu += (y[i] * vis[i]) / (specdata[i].Mol_mass() * phi);
    kappa += (specdata[i].ThermalConduct(Temp)*y[i]) / 
      (y[i] + (specdata[i].Mol_mass()) * 1.065 * (phi-1.0));

  } // endfor - species

  // finalize
  MW = 1.0/MW;
  g_ref = Cp_ref / (Cp_ref - (R_UNIVERSAL/MW));

  // clean up memory
#ifndef STATIC_NUMBER_OF_SPECIES
  delete[] vis;
#endif

}

