/********************** LESPremixed2DState.cc ******************************
  Constructors for LESPremixed2Dstate which handles all the physical
  variables and mixture rules associated with multi-species
  chemically reacting flows.

   associated files:
           LESPremixed2DState.h    

*********************************************************************/

//DEBUGGING FLAG FOR FIGUREING OUT proper NS-1 setup.
//#define _NS_MINUS_ONE

#ifndef _LESPREMIXED2D_STATE_INCLUDED 
#include "LESPremixed2DState.h"
#endif // _LESPREMIXED2D_STATE_INCLUDED    

/***********************************************************/
//Static members initialization

int LESPremixed2D_pState::ns = 1;
int LESPremixed2D_pState::nscal = 0;
int LESPremixed2D_pState::NUM_VAR_LESPREMIXED2D = NUM_LESPREMIXED2D_VAR_SANS_SPECIES; 
NASARP1311data* LESPremixed2D_pState::specdata=NULL;
Reactionset LESPremixed2D_pState::React;
Set_scalar LESPremixed2D_pState::Scal_sys;
double LESPremixed2D_pState::low_temp_range = 200.0;
double LESPremixed2D_pState::high_temp_range = 300.0;
double LESPremixed2D_pState::Mref=0.5;
double* LESPremixed2D_pState::Schmidt=NULL;
SubfilterScaleModels LESPremixed2D_pState::SFSmodel;
double LESPremixed2D_pState::filter_width = 0.0;


int LESPremixed2D_cState::ns = 1;
int LESPremixed2D_cState::nscal = 0; 
int LESPremixed2D_cState::NUM_VAR_LESPREMIXED2D = NUM_LESPREMIXED2D_VAR_SANS_SPECIES;   
NASARP1311data* LESPremixed2D_cState::specdata=NULL;   
double LESPremixed2D_cState::low_temp_range = 200.0;
double LESPremixed2D_cState::high_temp_range = 300.0;
double LESPremixed2D_cState::Mref=0.5;
double* LESPremixed2D_cState::Schmidt=NULL;
SubfilterScaleModels LESPremixed2D_cState::SFSmodel;
double LESPremixed2D_cState::filter_width = 0.0;



/***********************************************************************
 * LESPremixed2D_pState -- Create storage and assign premixed flame    *
 *                         static variables.                           *
 ***********************************************************************/
double LESPremixed2D_pState::laminar_speed = 0.38;
double LESPremixed2D_pState::laminar_thickness = 0.446E-3;
double LESPremixed2D_pState::TFactor = 1.0;

/***********************************************************************
 * LESPremixed2D_cState -- Create storage and assign premixed flame    *
 *                         static variables.                           *
 ***********************************************************************/
double LESPremixed2D_cState::laminar_speed = 0.38;
double LESPremixed2D_cState::laminar_thickness = 0.446E-3;
double LESPremixed2D_cState::TFactor = 1.0;



/***********************************************************/
/************* set_species_data ***************************
**********************************************************/
//set Global data for Species (STATIC, ie. so only call once! in LESPremixed2Dinput)
void LESPremixed2D_pState::set_species_data(const int &n,
					    const int &n2,
					    const string *S,
					    const char *PATH,
					    const double &Mr, 
					    const double* Sc,
					    const int &trans_data){ 
 
#ifdef STATIC_LESPREMIXED2D_SPECIES
  if( STATIC_LESPREMIXED2D_SPECIES < n) {                    //fix for mpi as called by each processor!!!!!!
    cerr<<"\n WARNING USING STATIC LESPREMIXED2D BUILT WITH "
	<< STATIC_LESPREMIXED2D_SPECIES 
        <<" SPECIES PREDEFINED, HOWEVER ASKING FOR "<< n <<endl; 
    exit(1); 
  }
#else 
  Deallocate(); //Clean up memory before changing ns
#endif

  ns = n;
  nscal = n2;
  NUM_VAR_LESPREMIXED2D = ns + nscal + NUM_LESPREMIXED2D_VAR_SANS_SPECIES;

  //read in NASA data for each species to be used  
  Deallocate_static();
  specdata = new NASARP1311data[ns]; 
  Schmidt = new double[ns];
  for(int i=0; i<ns; ++i){
    //overwrite default data  
    specdata[i].Getdata(S[i],PATH,trans_data);  
    Schmidt[i] = Sc[i];
  } 
  
  //set data temperature ranges for mixture
  Temp_low_range(); 
  Temp_high_range(); 
  
  //Preconditioning Reference Mach Number
  Mref = Mr;

  //setup initial arrays for mass fractions and scalars
  set_initial_values(); 
  set_initial_values_scal();
}   

//exact same thing with conserved variables, 
//should just point to pState data but for now duplication will work 
void LESPremixed2D_cState::set_species_data(const int &n,
					    const int &n2,
					    const string *S, 
					    const char *PATH,
					    const double &Mr, 
					    const double* Sc,
					    const int &trans_data){ 
#ifdef STATIC_LESPREMIXED2D_SPECIES
  if( STATIC_LESPREMIXED2D_SPECIES < n ) { 
    cerr<<"\n WARNING USING STATIC LESPREMIXED2D BUILT WITH "
	<< STATIC_LESPREMIXED2D_SPECIES 
        <<" SPECIES PREDEFINED, HOWEVER ASKING FOR "<< n <<endl; 
    exit(1); 
  }
#else 
  Deallocate(); //Clean up memory before changing ns
#endif 
 
  ns = n;
  nscal = n2;
  NUM_VAR_LESPREMIXED2D = ns + nscal + NUM_LESPREMIXED2D_VAR_SANS_SPECIES;

  //read in NASA data for each species to be used
  Deallocate_static();
  specdata = new NASARP1311data[ns];
  Schmidt = new double[ns];
  for(int i=0; i<ns; ++i){
    //overwrite default data  
    specdata[i].Getdata(S[i],PATH,trans_data);
    Schmidt[i] = Sc[i];  
  }  

  //set data temperature ranges for mixture
  Temp_low_range();
  Temp_high_range();

  //Preconditioning Reference Mach Number
  Mref = Mr;
  
  //setup initial arrays for mass fractions and scalars
  set_initial_values();
  set_initial_values_scal();
}      

/**************************************************************************
***************** LESPREMIXED2D_PSTATE CONSTRUCTORS ***********************
***************************************************************************/

/***********************************************************
 ****************** Mixture Rules **************************
 ***********************************************************/

/**************************************************
  mixture molecular mass (kg/mol)
***************************************************/
double LESPremixed2D_pState::Mass() const{
  // = 1 / sum(mass fration/mol_mass)
  double sum = 0.0;
  for(int i=0; i<ns; ++i){
    sum += spec[i].c/specdata[i].Mol_mass();
  }
  return 1.0/sum;
}

/**************************************************
  mixture gas constant  J/(kg*K)
***************************************************/
double LESPremixed2D_pState::Rtot(){
  // = sum ( mass fraction * species gas constant)
  double sum = 0.0;
  for(int i=0; i<ns; ++i){
    sum += spec[i].c * specdata[i].Rs();
  }
  return sum;
}

double LESPremixed2D_pState::Rtot() const{
  // = sum ( mass fraction * species gas constant)
  double sum = 0.0;
  for(int i=0; i<ns; ++i){
    sum += spec[i].c * specdata[i].Rs();
  }
  return sum;
}

/**************************************************
  mixture Heat Capacity (const pressure) J/(kg*K)
***************************************************/
double LESPremixed2D_pState::Cp(void) const{
  // = sum ( mass fraction * species Cp) 
  double Temp = T();
  double sum = 0.0;
  for(int i=0; i<ns; ++i){
    sum += spec[i].c*specdata[i].HeatCapacity_p(Temp);
  }
  return sum;
}

double LESPremixed2D_pState::Cp(const double& TEMP) const{
  // = sum ( mass fraction * species Cp) 
  double sum = 0.0;
  for(int i=0; i<ns; ++i){
    sum += spec[i].c*specdata[i].HeatCapacity_p(TEMP);
  }
  return sum;
}

/**************************************************
  mixture Heat Capacity (const volume) J/(kg*K)
***************************************************/
double LESPremixed2D_pState::Cv(void) const{
  // = sum ( mass fraction * species Cv)  
  double Temp = T();
  double sum = 0.0;
  for(int i=0; i<ns; ++i){
    sum += spec[i].c*specdata[i].HeatCapacity_v(Temp);
  }
  return sum;
}

/**************************************************
  mixture Heat Ratio gamma J/(kg*K)
***************************************************/
double LESPremixed2D_pState::g(void) const{
  // = Cp / Cv  
  return Cp()/Cv();
}

/**************************************************
  Specific Internal Energy
***************************************************/
// etotal = sensible & chemical
double LESPremixed2D_pState::e(void) const{
  // = sum (mass fraction * species e) 
  double sum = 0.0;
  double Temp = T();
  for(int i=0; i<ns; ++i){ //(Enthalpy(Temp) - (R/mol_mass)*Temp)
    sum += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform()
    		      - specdata[i].Rs()*Temp);
  }
  return sum;
}

// double LESPremixed2D_pState::e(void) const{
//   // = sum (mass fraction * species e) 
//   double sum = 0.0;
//   double cs = 0.0;
//   double Temp = T();
//   for(int i=0; i<ns-1; ++i){ //(Enthalpy(Temp) - (R/mol_mass)*Temp)
//     sum += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform()
//     		      - specdata[i].Rs()*Temp);
//     cs += spec[i].c;
//   }
//   sum += (ONE - cs )*(specdata[ns-1].Enthalpy(Temp) + specdata[ns-1].Heatofform() - specdata[ns-1].Rs()*Temp);
//   return sum;
// }

// reference internal energy e + heat of formation - offset
double LESPremixed2D_pState::eref(void)const{
  // = sum (mass fraction * species e) 
  double sum = 0.0;
  double Temp = T();
  for(int i=0; i<ns; ++i){ //(Enthalpy(Temp) - (R/mol_mass)*Temp)
    sum += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform()
		      	 - specdata[i].DeltaHref() - specdata[i].Rs()*Temp);
  }
  return sum;
}

// internal energy with no heat of formation included (sensible)
double LESPremixed2D_pState::es(void) const{
  // = sum (mass fraction * species e) 
  double sum = 0.0;
  double Temp = T();
  for(int i=0; i<ns; ++i){ //(Enthalpy(Temp) - (R/mol_mass)*Temp)
    sum += spec[i].c*(specdata[i].Enthalpy(Temp) - specdata[i].Rs()*Temp);
  }
  return sum;
}

/************************************************** 
   Specific absolute enthalpy = 
          mass fraction * (hsensible + heatofform)
***************************************************/
double LESPremixed2D_pState::h(void) const{
  // = sum (mass fraction * species h) 
  double sum = 0.0;  
  double Temp = T();
  for(int i=0; i<ns; ++i){
    sum += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform());
  }
  return sum;
}

double LESPremixed2D_pState::h(const double &Temp) const{
  // = sum (mass fraction * species h) 
  double sum = 0.0;  
  for(int i=0; i<ns; ++i){
    sum += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform());
  }
  return sum;
}


double LESPremixed2D_pState::href(void)const{
  // = sum (mass fraction * species h) 
  double sum = 0.0;  
  double Temp = T();
  for(int i=0; i<ns; ++i){ 
    sum += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform()
		      - specdata[i].DeltaHref() );
  }
  return sum;
}

double LESPremixed2D_pState::hs(void) const{
  // = sum (mass fraction * species h) 
  double sum = 0.0;  
  double Temp = T();
  for(int i=0; i<ns; ++i){
    sum += spec[i].c*(specdata[i].Enthalpy(Temp));
  }
  return sum;
}

double LESPremixed2D_pState::hs(const double &Temp) const{
  // = sum (mass fraction * species h) 
  double sum = 0.0;  
  for(int i=0; i<ns; ++i){
    sum += spec[i].c*(specdata[i].Enthalpy(Temp));
  }
  return sum;
}

/**************************************************
   Derivative of specific enthalpy dh/dT
   actually is just Cp as Cp = (dh/dT)_p
***************************************************/
double LESPremixed2D_pState::hprime() const{
 double sum = 0.0; 
 double Temp = T();
 for(int i=0; i<ns; ++i){
   sum += spec[i].c*specdata[i].Enthalpy_prime(Temp);
 }
 return (sum);
}

double LESPremixed2D_pState::hprime(double &Temp) const{
 double sum = 0.0;  
 for(int i=0; i<ns; ++i){
   sum += spec[i].c*specdata[i].Enthalpy_prime(Temp);
 }
 return (sum);
}

/**************************************************
  Total Energy
***************************************************/
double LESPremixed2D_pState::E(void) const{
  // E = rho*(e + velocity^2 +k)  
  return (rho * (e() + HALF*v.sqr() + k()) ); 
}

/**************************************************
  Total Enthalpy
***************************************************/
double LESPremixed2D_pState::H(void) const{
  // H = h + velocity^2 
  return(rho * (h() + HALF*v.sqr() + 5.0*k()/3.0) ); 
}

double LESPremixed2D_pState::Hs(void) const{
  // H = h + velocity^2+k
  return (rho * (hs() + HALF*v.sqr() + 5.0*k()/3.0) );
}

/**************************************************
  Viscosity 
  using Wilke [1950] formulation
***************************************************/
double LESPremixed2D_pState::mu() const{
  double sum =0.0; 
  double Temp = T();
#ifdef STATIC_LESPREMIXED2D_SPECIES
  double  vis[STATIC_LESPREMIXED2D_SPECIES];
#else
  double  *vis = new double[ns];
#endif

  for(int i=0; i<ns; ++i){
    double phi = 0.0;
    for (int j=0; j<ns; ++j){
      if(i == 0){
	vis[j] = specdata[j].Viscosity(Temp);
      }
      phi += (spec[j].c / specdata[j].Mol_mass())*
	pow(ONE + sqrt(vis[i]/vis[j])*
	    pow(specdata[j].Mol_mass()/specdata[i].Mol_mass(),0.25),2.0)/
	sqrt(EIGHT*(ONE +specdata[i].Mol_mass()/specdata[j].Mol_mass()));
    }
    sum += (spec[i].c * vis[i]) / 
      (specdata[i].Mol_mass() * phi);
  }  

#ifndef STATIC_LESPREMIXED2D_SPECIES  
  delete[] vis;
#endif

#ifdef THICKENED_FLAME_ON
  return (flame.WF*flame.TF)*sum;
#else
  return sum;
#endif
}

/***********************************************************
  Molecular viscosity is function of Temperature
  This derivative is needed by Jacobian
 ***********************************************************/
double LESPremixed2D_pState::dmudT(void) const{
  double sum =0.0; 
  double Temp = T();
  for(int i=0; i<ns; ++i){
    double phi = 0.0;
    for (int j=0; j<ns; ++j){
      phi += (spec[j].c / specdata[j].Mol_mass())*
	pow(1.0 + sqrt(specdata[i].dViscositydT(Temp)/specdata[j].dViscositydT(Temp))*
	    pow(specdata[j].Mol_mass()/specdata[i].Mol_mass(),0.25),2.0)/
       sqrt(8.0*(1.0 +specdata[i].Mol_mass()/specdata[j].Mol_mass()));
    }
    sum += (spec[i].c * specdata[i].dViscositydT(Temp) ) / 
      (specdata[i].Mol_mass() * phi);
  }

#ifdef THICKENED_FLAME_ON
  return (flame.WF*flame.TF)*sum;
#else  
  return sum;
#endif
}

/****************************************************
  Thermal Conductivity - Mason & Saxena (1958)  W/(m*K)
****************************************************/
double LESPremixed2D_pState::kappa(void) const{
  double sum = 0.0;  
  double Temp = T();
#ifdef STATIC_LESPREMIXED2D_SPECIES
  double  vis[STATIC_LESPREMIXED2D_SPECIES];
#else
  double  *vis = new double[ns];
#endif

  for(int i=0; i<ns; ++i){
    double phi = 0.0;
    for (int j=0; j<ns; ++j){
      if(i == 0){
	vis[j] = specdata[j].Viscosity(Temp);
      }
      if(i != j){
	phi += (spec[j].c / specdata[j].Mol_mass())*
	  pow(ONE + sqrt(vis[i]/vis[j])*
	      pow(specdata[j].Mol_mass()/specdata[i].Mol_mass(),0.25),2.0)/
	  sqrt(EIGHT*(ONE +specdata[i].Mol_mass()/specdata[j].Mol_mass()));
      }
    }
 
    sum += (specdata[i].ThermalConduct(Temp)*spec[i].c) / 
      (spec[i].c + (specdata[i].Mol_mass()) * 1.065 * phi);
  }  
  
  //or Coffee and Heimerl (1981)
//   double one =0.0; double two=0.0;
//   for (int j=0; j<ns; ++j){
//     one += spec[i].c*specdata[i].ThermalConduct(Temp);
//     two += spec[i].c/specdata[i].ThermalConduct(Temp);
//   }
//   sum = HALF*(one + ONE/two);
#ifndef STATIC_LESPREMIXED2D_SPECIES
  delete[] vis;
#endif

#ifdef THICKENED_FLAME_ON
  return (flame.WF*flame.TF)*sum;
#else
  return sum;
#endif
}

/**************************************************
  Turbulence model related parameters
***************************************************/
double LESPremixed2D_pState::mu_t(const Tensor2D &strain_rate) const{
  double mut = rho*SFSmodel.eddy_viscosity_Smagorinsky(strain_rate, filter_width);
#ifdef THICKENED_FLAME_ON
  return (flame.WF*flame.TF)*mut;
#else
  return mut;
#endif
}

double LESPremixed2D_pState::Pr_turb(void) const{
  return (0.75); 
}

double LESPremixed2D_pState::Sc_turb(void) const{
  return (0.9);
}

double LESPremixed2D_pState::Kappa_turb(const double &mut) const{
  return (mut*Cp()/Pr_turb()); 
}

double LESPremixed2D_pState::Dm_turb(const double &mut) const{
  return (mut/(rho*Sc_turb()));
}     



/**************************************************
  polytropic heat ratio mixture gamma J/(kg*K)
  assuming T=200K as the temperature.
***************************************************/
double LESPremixed2D_pState::gamma_guess(void) const{  
  double sum1 = ZERO;     double sum2 = ZERO;
  double gamma_s = ZERO;  double Temp = 200.0;

  for(int i=0; i<ns; ++i){
    sum1 += spec[i].c*specdata[i].Rs();
    gamma_s = ONE/(ONE - specdata[i].Rs()/ specdata[i].HeatCapacity_p(Temp));
    sum2 += ((spec[i].c*specdata[i].Rs()) / (gamma_s - ONE)); 
  }
  return (ONE + sum1/sum2);
}

/**************************************************
   Temperature derived from specific sensible
   enthalpy
***************************************************/
double LESPremixed2D_pState::T(double &h_s) const{
  double T = ZERO;
  //--------- Initial Guess ------------------------------//
  //using a polytropic gas assumption with gamma@200;
  double Tguess = (h_s/Rtot())*(ONE/(ONE/(gamma_guess() - ONE) +ONE));
  //--------- newtons method to get T --------------------//
  int numit =0;
  double Tmin = low_temp_range;
  double Tmax = high_temp_range;
  //check for start value
  if(Tguess > Tmin && Tguess < Tmax){
    T=Tguess;
  } else { 
   T=Tmin;
  }

  double fa = hs(Tmin) - h_s;
  double fn = hs(T) - h_s; 
  double dfn = hprime(T);
  
  while( fabs(Tmax-Tmin) > CONV_TOLERANCE && fabs(fn) 
	 > CONV_TOLERANCE && numit<20 && T >= low_temp_range){  
    if(T >= Tmin && T <= Tmax){
      T = T - fn/dfn;
      if(T >= Tmax) T = HALF*(Tmax + Tmin);	
      //Bisection
    } else {
      T = HALF*(Tmax + Tmin);
    } 
    fn = hs(T) - h_s;  
    dfn = hprime(T); 
    //change bisection range
    if ( fa*fn <=ZERO){
      Tmax = T;
    } else {
      Tmin = T;
      fa = fn;
    }
    numit++;
  }
  if (numit>=19 || T <= low_temp_range){
    T = max(Tguess,low_temp_range); 
    cout<<"\nTemperature didn't converge in LESPremixed2D_cState::T(void)";
    cout<<" with polytopic Tguess "<<Tguess<<", or lower than Tmin "
	<<low_temp_range<<" using "<<T;
  }
  return T;
}
    
/****************************************************
  Speed of sound using 
  a^2 = dip/dirho + p/rho^2( die/dip)^-1
  from eigenvalue analysis using e =f(p,rho)
****************************************************/
double LESPremixed2D_pState::a(void){
  double sum;
  double RTOT= Rtot();
 
  sum = (p/rho)*(RTOT/( hprime() - RTOT) + ONE);
  //could also just call sqrt(g()*Rtot()*T());
  return sqrt(sum);
}

double LESPremixed2D_pState::a(void) const{
  double sum;
  double RTOT= Rtot();
  sum = (p/rho)*(RTOT/( hprime() - RTOT) + ONE);
  //could also just call sqrt(g()*Rtot()*T());
  return sqrt(sum);
}

//turbulence modified  speed of sound
double LESPremixed2D_pState::amodified(void) const{
  return sqrt( g()*pmodified()/rho );
}



/******************************************************
 Calculating the thermal diffusion component of 
 the heat flux vector (qflux)
  sum( hs * Ds * grad cs)
*******************************************************/
Vector2D LESPremixed2D_pState::thermal_diffusion(void) const{
  Vector2D sum;
  sum.zero();
  double Temp = T();
  //problems with Species overloaded operators
  for(int i=0; i<ns; ++i){
#ifdef THICKENED_FLAME_ON
    sum  +=  (specdata[i].Enthalpy(Temp) + specdata[i].Heatofform())
      * (flame.WF*flame.TF)*spec[i].diffusion_coef * spec[i].gradc;
#else 
    sum  +=  (specdata[i].Enthalpy(Temp) + specdata[i].Heatofform())
      * spec[i].diffusion_coef * spec[i].gradc;
#endif
  }
  return sum;
}

/*******************************************************************
 *************** INVISCID FLUXES ***********************************
 *******************************************************************/

/*************************************************************
 * LESPremixed2D_pState::Fx -- Inviscid flux (x-direction).  *
 *************************************************************/
LESPremixed2D_cState LESPremixed2D_pState::Fx(void) const{
  LESPremixed2D_cState Temp;
  
  Temp.rho = rho*v.x;
  Temp.rhov.x = rho*sqr(v.x) + pmodified();
  Temp.rhov.y = rho*v.x*v.y;
  Temp.E = v.x*H();

  if(nscal) for(int i=0; i<nscal; ++i) Temp.rhoscalar[i] = rho*v.x*scalar[i];
   
  //multispecies transport
  for(int i=0; i<ns; ++i){
    Temp.rhospec[i].c = rho*v.x*spec[i].c;
  }

  return (Temp);
}

/***************************************************************
 *************** INVISCID FLUX JACOBIANS ***********************
 ***************************************************************/

/**********************************************************
 * LESPremixed2D_pState::dFIdU -- Invisicd Flux Jacobian  * 
 **********************************************************/
void dFIdU(DenseMatrix &dFdU, const LESPremixed2D_pState &W, const int Flow_Type) {
  
  //cout<<"\n USING DFIDU \n";
  
  //int num_species = dFdU.get_n() - NUM_LESPREMIXED2D_VAR_SANS_SPECIES; 
  int num_species=W.ns-1;

  double Temp = W.T();
  double pt = W.pmodified();
  double Rt = W.Rtot();
  double C_p = W.Cp();
  double ht = W.h();
  double kk = W.k();
  double denominator = (C_p/Rt - ONE);

  double phi = ZERO;   
  for(int i=0; i<num_species; ++i){ 
#ifdef _NS_MINUS_ONE
    phi += W.spec[i].c*(W.specdata[i].Enthalpy(Temp) + W.specdata[i].Heatofform() - C_p*Temp*W.specdata[i].Rs()/Rt
			-(W.specdata[W.ns-1].Enthalpy(Temp)+W.specdata[W.ns-1].Heatofform()  -C_p*Temp*W.specdata[W.ns-1].Rs()/Rt));   
#else
    phi += W.spec[i].c*(W.specdata[i].Enthalpy(Temp) + W.specdata[i].Heatofform() - C_p*Temp*W.specdata[i].Rs()/Rt);
#endif 
  }

  dFdU(0,1) += ONE;
  
  if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY) {
    dFdU(1,0) += ( (C_p/Rt)*( - W.v.x*W.v.x) + HALF*(THREE*W.v.x*W.v.x + W.v.y*W.v.y) - ht - 5.0*kk/3.0 + C_p*pt/(W.rho*Rt) + phi )/denominator;
  } else if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {
    dFdU(1,0) += ( (C_p/Rt)*( - W.v.x*W.v.x) + HALF*(THREE*W.v.x*W.v.x + W.v.y*W.v.y) - ht + C_p*pt/(W.rho*Rt) + phi )/denominator;
  }

  dFdU(1,1) += W.v.x*(TWO*C_p/Rt-THREE)/denominator; 
  dFdU(1,2) -= W.v.y/denominator;
  dFdU(1,3) += ONE/denominator;
  dFdU(2,0) -= W.v.x*W.v.y;
  dFdU(2,1) += W.v.y;
  dFdU(2,2) += W.v.x;

  if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY) {
    dFdU(3,0) += W.v.x*( W.v.x*W.v.x + W.v.y*W.v.y + C_p*pt/(W.rho*Rt) - (C_p/Rt)*( HALF*(W.v.x*W.v.x + W.v.y*W.v.y) + 5.0*kk/3.0 + ht) + phi)/denominator;
  } else if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {
    dFdU(3,0) += W.v.x*( W.v.x*W.v.x + W.v.y*W.v.y + 5.0*kk/3.0 + C_p*pt/(W.rho*Rt) - (C_p/Rt)*( HALF*(W.v.x*W.v.x + W.v.y*W.v.y) + 5.0*kk/3.0 + ht) + phi)/denominator;
  }

  dFdU(3,1) += ht + HALF*(W.v.x*W.v.x + W.v.y*W.v.y) + 5.0*kk/3.0 - W.v.x*W.v.x/denominator;
  dFdU(3,2) -= W.v.x*W.v.y/denominator;
  dFdU(3,3) += W.v.x*C_p/denominator/Rt;


  int NUM_VAR = NUM_LESPREMIXED2D_VAR_SANS_SPECIES + W.nscal;

  //Species
  for(int i = 0; i<num_species; ++i){ 
#ifdef _NS_MINUS_ONE
    dFdU(1,NUM_VAR+i) -= (W.specdata[i].Enthalpy(Temp) + W.specdata[i].Heatofform() - C_p*Temp*W.specdata[i].Rs()/Rt
			  -(W.specdata[W.ns-1].Enthalpy(Temp)+W.specdata[W.ns-1].Heatofform() - C_p*Temp*W.specdata[W.ns-1].Rs()/Rt))/denominator; 
#else
    dFdU(1,NUM_VAR+i) -= (W.specdata[i].Enthalpy(Temp) + W.specdata[i].Heatofform() 
			  - C_p*Temp*W.specdata[i].Rs()/Rt)/denominator;	
#endif 
	
    dFdU(3,NUM_VAR+i) = W.v.x*dFdU(1,NUM_VAR+i);    
    dFdU(NUM_VAR+i, 0) -= W.spec[i].c*W.v.x ;
    dFdU(NUM_VAR+i, 1) += W.spec[i].c ;
    dFdU(NUM_VAR+i,NUM_VAR+i) += W.v.x ;        
  }

  
  //Scalars
  NUM_VAR = NUM_LESPREMIXED2D_VAR_SANS_SPECIES;  
  for(int i = 0; i<W.nscal; ++i){ 
    // column corresponding to k, I don't know about other scalars...
    if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {
      dFdU(1,NUM_VAR+i) -= (5.0/3.0)/denominator;
      dFdU(3,NUM_VAR+i) += W.v.x*dFdU(1,NUM_VAR+i);
    }
    
    // row corresponding to the scalar
    dFdU(NUM_VAR+i, 0) -= W.scalar[i]*W.v.x;
    dFdU(NUM_VAR+i, 1) += W.scalar[i];
    // diagonal
    dFdU(NUM_VAR+i,NUM_VAR+i) += W.v.x;
  }
 
}

// Finite difference check of dFxdU
void dFIdU_FD(DenseMatrix &dFdU, const LESPremixed2D_pState &WW, const int Flow_Type) {

  LESPremixed2D_cState UU = U(WW);
  LESPremixed2D_cState A,C;
  LESPremixed2D_pState B,D;
  double perturb = 5e-6;
  double a;

  for(int jcol=0; jcol<(dFdU.get_n()); ++jcol){    
    A = UU;  C = UU;
    if( jcol <NUM_LESPREMIXED2D_VAR_SANS_SPECIES) {
      A[jcol+1] += perturb*max(ONE,UU[jcol+1]); 
      C[jcol+1] -= perturb*max(ONE,UU[jcol+1]);
    } else {                                       //enforce sum(ci) = 1;
      a = perturb*max(ONE,UU[jcol+1]); 
      A[jcol+1] += a;
      A[WW.NUM_VAR_LESPREMIXED2D] -= a;      
      C[jcol+1] -= a;
      C[WW.NUM_VAR_LESPREMIXED2D] += a;
    }   
    B = W(A);  D = W(C);    
    A = B.Fx();  C = D.Fx();
    for(int irow=0; irow<(dFdU.get_n()); ++irow){
      dFdU(irow,jcol) = ( A[irow+1] - C[irow+1])/(TWO*perturb*max(ONE, UU[jcol+1]));      
    }
  } 

}

/*********************************************************
 * LESPremixed2D_pState::dFIdW -- Invisicd Flux Jacobian *
 *********************************************************/
void dFIdW(DenseMatrix &dFdW, const LESPremixed2D_pState &W, const int Flow_Type) {

  //int num_species = dFdW.get_n() - NUM_LESPREMIXED2D_VAR_SANS_SPECIES; 
  int num_species=W.ns-1;

  double Temp = W.T();
  double pt = W.pmodified();
  double Rt = W.Rtot();
  double C_p = W.Cp();
  double ht = W.h();
  double kk = W.k();

  dFdW(0,0) = W.v.x;
  dFdW(0,1) = W.rho;
  dFdW(1,0) = W.v.x*W.v.x;
  dFdW(1,1) = TWO*W.rho*W.v.x; 
  dFdW(1,3) = ONE;
  dFdW(2,0) = W.v.x*W.v.y;
  dFdW(2,1) = W.rho*W.v.y;
  dFdW(2,2) = W.rho*W.v.x;
  dFdW(3,0) = (HALF*(W.v.x*W.v.x+W.v.y*W.v.y) + ht + 5.0*kk/3.0)*W.v.x - C_p*pt/(W.rho*Rt)*W.v.x;
  dFdW(3,1) = W.rho*W.v.x*W.v.x+ W.rho*(ht + 5.0*kk/3.0 + HALF*(W.v.x*W.v.x+W.v.y*W.v.y));
  dFdW(3,2) = W.rho*W.v.x*W.v.y;
  dFdW(3,3) = W.v.x*C_p/Rt;

  int NUM_VAR = NUM_LESPREMIXED2D_VAR_SANS_SPECIES + W.nscal;

  //Species
  for(int i = 0; i<num_species; ++i){ 
#ifdef _NS_MINUS_ONE
    dFdW(3,NUM_VAR+i) = W.rho*W.v.x*((W.specdata[i].Enthalpy(Temp) + W.specdata[i].Heatofform() - C_p*Temp*W.specdata[i].Rs()/Rt) - 
				     (W.specdata[W.ns-1].Enthalpy(Temp) + W.specdata[W.ns-1].Heatofform() - C_p*Temp*W.specdata[W.ns-1].Rs()/Rt));    
#else
    dFdW(3,NUM_VAR+i) =  W.rho*W.v.x*(W.specdata[i].Enthalpy(Temp) + W.specdata[i].Heatofform() - C_p*Temp*W.specdata[i].Rs()/Rt); 
#endif 

     dFdW(NUM_VAR+i, 0) = W.spec[i].c*W.v.x ;
     dFdW(NUM_VAR+i, 1) = W.rho*W.spec[i].c ;
     dFdW(NUM_VAR+i,NUM_VAR+i) =W.rho*W.v.x ;        
  }

  //Scalars
  NUM_VAR = NUM_LESPREMIXED2D_VAR_SANS_SPECIES;  
  for(int i = 0; i<W.nscal; ++i){ 
    // column corresponding to k, I don't know about other scalars...
    if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {
      dFdW(3,NUM_VAR+i) = (5.0/3.0)*W.rho*W.v.x;
    }
    
    // row corresponding to the scalar
    dFdW(NUM_VAR+i, 0) = W.scalar[i]*W.v.x;
    dFdW(NUM_VAR+i, 1) = W.rho*W.scalar[i];
    // diagonal
    dFdW(NUM_VAR+i,NUM_VAR+i) = W.rho*W.v.x;
  }  
  
}

// Finite difference check of dFxdW
void dFIdW_FD(DenseMatrix &dFdW, const LESPremixed2D_pState &W, const int Flow_Type) {

  LESPremixed2D_pState A,C;
  LESPremixed2D_cState B,D;
  double perturb = 5e-6;
  double a;

  for(int jcol=0; jcol<(dFdW.get_n()); ++jcol){    
    A.Copy(W);  C.Copy(W);
    if( jcol <NUM_LESPREMIXED2D_VAR_SANS_SPECIES) {
      A[jcol+1] += perturb*max(ONE,W[jcol+1]); 
      C[jcol+1] -= perturb*max(ONE,W[jcol+1]);
    } else {                                       //enforce sum(ci) = 1;
      a = perturb*max(ONE,W[jcol+1]); 
      A[jcol+1] += a;
      A[W.NUM_VAR_LESPREMIXED2D] -= a;      
      C[jcol+1] -= a;
      C[W.NUM_VAR_LESPREMIXED2D] += a;
    }   
    B = A.Fx();   D = C.Fx();
    for(int irow=0; irow<(dFdW.get_n()); ++irow){
      dFdW(irow,jcol) = ( B[irow+1] - D[irow+1])/(TWO*perturb*max(ONE, W[jcol+1]));      
    }
  } 

}

/*****************************************************************************
 * LESPremixed2D_pState::dWdU -- Primitive/Conserved transformation Jacobian *
 *****************************************************************************/
void LESPremixed2D_pState::dWdU(DenseMatrix &dWdQ, const int Flow_Type) const{

  //int num_species = dWdQ.get_n() - NUM_LESPREMIXED2D_VAR_SANS_SPECIES; //ns -1 or ns
  int num_species=ns-1;

  double Temp = T();
  double pt = pmodified();
  double Rt = Rtot();
  double C_p = Cp();
  double denominator = (C_p/Rt - ONE);
  double kk = k();

  dWdQ(0,0) = ONE;
  dWdQ(1,0) = -v.x/rho;
  dWdQ(1,1) = ONE/rho;
  dWdQ(2,0) = -v.y/rho;
  dWdQ(2,2) = ONE/rho; 

  double phi = ZERO;   
  for(int i=0; i<num_species; ++i){  
#ifdef _NS_MINUS_ONE
    phi += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform() - C_p*Temp*specdata[i].Rs()/Rt
    - (specdata[ns-1].Enthalpy(Temp)+specdata[ns-1].Heatofform() - C_p*Temp*specdata[ns-1].Rs()/Rt));
#else
    phi += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform() - C_p*Temp*specdata[i].Rs()/Rt);
#endif
  }

  if(Flow_Type == FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY) { 
    dWdQ(3,0) = (HALF*(v.x*v.x+v.y*v.y) - h() - 5.0*kk/3.0 + C_p*pt/(rho*Rt) + phi)/denominator;
  } else if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {
    dWdQ(3,0) = (HALF*(v.x*v.x+v.y*v.y) - h() + C_p*pt/(rho*Rt) + phi)/denominator;
  }
  dWdQ(3,1) = -v.x/denominator;
  dWdQ(3,2) = -v.y/denominator;
  dWdQ(3,3) = ONE/denominator;


  //Species
  int NUM_VAR = NUM_LESPREMIXED2D_VAR_SANS_SPECIES + nscal;
  for(int i=0; i<num_species; ++i){  
#ifdef _NS_MINUS_ONE  
    dWdQ(3, NUM_VAR+i) = - (specdata[i].Enthalpy(Temp) + specdata[i].Heatofform() - C_p*Temp*specdata[i].Rs()/Rt
			    -(specdata[ns-1].Enthalpy(Temp)+specdata[ns-1].Heatofform() - C_p*Temp*specdata[ns-1].Rs()/Rt))/denominator;
#else
    dWdQ(3, NUM_VAR+i) = - (specdata[i].Enthalpy(Temp) + specdata[i].Heatofform() - C_p*Temp*specdata[i].Rs()/Rt)/denominator;
#endif
    dWdQ(NUM_VAR+i, 0) = - spec[i].c/rho;
    dWdQ(NUM_VAR+i, NUM_VAR+i) = ONE/rho;
  }

  //Scalars
  NUM_VAR = NUM_LESPREMIXED2D_VAR_SANS_SPECIES;
  for(int i = 0; i<nscal; ++i){ 
    // column corresponding to k, I don't know about other scalars...
    if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {
      dWdQ(3,NUM_VAR+i) = -(5.0/3.0)/denominator;
    }

    // row corresponding to the scalar
    dWdQ(NUM_VAR+i, 0) = - scalar[i]/rho;
    // diagonal
    dWdQ(NUM_VAR+i,NUM_VAR+i) = ONE/rho;
  }

}

// Finite difference check of dWdU
// shows error in (3,0) ie dp/rho due to pertubing rho and cState T() calc.
void LESPremixed2D_pState::dWdU_FD(DenseMatrix &dWdQ, const int Flow_Type){

  LESPremixed2D_cState UU = U(*this);
  LESPremixed2D_cState A,C;
  LESPremixed2D_pState B,D;
  double perturb = 5e-6;
  double a;

  for(int jcol=0; jcol<dWdQ.get_n(); ++jcol){    
    A = UU;    C = UU;
    if( jcol <NUM_LESPREMIXED2D_VAR_SANS_SPECIES) {
      A[jcol+1] += perturb*max(ONE,A[jcol+1]); 
      C[jcol+1] -= perturb*max(ONE,C[jcol+1]);
    } else {                                       //enforce sum(ci) = 1;
      a =  perturb*max(ONE,A[jcol+1]); 
      A[jcol+1] += a;
      A[NUM_VAR_LESPREMIXED2D] -= a;      
      C[jcol+1] -= a;
      C[NUM_VAR_LESPREMIXED2D] += a;
    }
    B = W(A);    D = W(C);
    for(int irow=0; irow<(dWdQ.get_n()); ++irow){
      dWdQ(irow,jcol) = ( B[irow+1] - D[irow+1])/(TWO*perturb*max(ONE,UU[jcol+1]));     
    }
  } 

}

/*****************************************************************************
 * LESPremixed2D_pState::dUdW -- Conserved/Primitive transformation Jacobian *
 *****************************************************************************/
void LESPremixed2D_pState::dUdW(DenseMatrix &dQdW, const int Flow_Type){

  //int num_species = dQdW.get_n() - NUM_LESPREMIXED2D_VAR_SANS_SPECIES; //ns -1 or ns
  int num_species=ns-1;
  
  double Temp = T();
  double pt = pmodified();
  double Rt = Rtot();
  double C_p = Cp();
  double kk = k();

  dQdW(0,0) =  ONE;
  dQdW(1,0) =  v.x;
  dQdW(1,1) =  rho;  
  dQdW(2,0) =  v.y;
  dQdW(2,2) =  rho;  
  dQdW(3,0) =  HALF*(v.x*v.x + v.y*v.y) + h() + 5.0*kk/3.0 - C_p*pt/(rho*Rt);
  dQdW(3,1) =  rho*v.x;
  dQdW(3,2) =  rho*v.y;
  dQdW(3,3) =  C_p/Rt - ONE;

  //Species
  int NUM_VAR = NUM_LESPREMIXED2D_VAR_SANS_SPECIES + nscal;
  for(int i=0; i<(num_species); ++i){  
#ifdef _NS_MINUS_ONE   
     dQdW(3,NUM_VAR+i) = rho*(specdata[i].Enthalpy(Temp)+specdata[i].Heatofform() - C_p*Temp*specdata[i].Rs()/Rt
			      - (specdata[ns-1].Enthalpy(Temp)+specdata[ns-1].Heatofform() - C_p*Temp*specdata[ns-1].Rs()/Rt));
#else
    dQdW(3,NUM_VAR+i) = rho*(specdata[i].Enthalpy(Temp)+specdata[i].Heatofform() - C_p*Temp*specdata[i].Rs()/Rt);
#endif
    dQdW(NUM_VAR+i,0) = spec[i].c;
    dQdW(NUM_VAR+i,NUM_VAR+i) =  rho;
  }

  //Scalars
  NUM_VAR = NUM_LESPREMIXED2D_VAR_SANS_SPECIES;
  for(int i = 0; i<nscal; ++i){
    // column corresponding to k, I don't know about other scalars...
    if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {
      dQdW(3,NUM_VAR+i) = (5.0/3.0)*rho;
    }

    // row corresponding to the scalar
    dQdW(NUM_VAR+i, 0) = scalar[i];
    // diagonal
    dQdW(NUM_VAR+i,NUM_VAR+i) = rho;
  }

}

// Finite difference check of dWdU
void LESPremixed2D_pState::dUdW_FD(DenseMatrix &dUdW, const int Flow_Type){

  LESPremixed2D_pState W(*this);
  LESPremixed2D_pState A,C;
  LESPremixed2D_cState B,D;
  double perturb = 5e-6;
  double a;

  for(int jcol=0; jcol<(dUdW.get_n()); ++jcol){    
    A.Copy(*this);  C.Copy(*this);
    if( jcol <NUM_LESPREMIXED2D_VAR_SANS_SPECIES) {
      A[jcol+1] += perturb*max(ONE,A[jcol+1]); 
      C[jcol+1] -= perturb*max(ONE,C[jcol+1]);
    } else {                                       //enforce sum(ci) = 1;
      a =  perturb*max(ONE,A[jcol+1]); 
      A[jcol+1] += a;
      A[NUM_VAR_LESPREMIXED2D] -= a;      
      C[jcol+1] -= a;
      C[NUM_VAR_LESPREMIXED2D] += a;
    }   
    B = U(A);   D = U(C);
    for(int irow=0; irow<(dUdW.get_n()); ++irow){
      dUdW(irow,jcol) = ( B[irow+1] - D[irow+1])/(TWO*perturb*max(ONE, W[jcol+1]));      
    }
  } 

}


/*******************************************************************
 ***************** EIGENVALUES *************************************
 *******************************************************************/

/****************************************************************
 * LESPremixed2D_pState::lambda -- Eigenvalue(s) (x-direction). *
 ****************************************************************/
LESPremixed2D_pState LESPremixed2D_pState::lambda_x(void) const {
  double c = amodified();
  LESPremixed2D_pState Temp;
  Temp.rho = v.x - c;
  Temp.v.x = v.x;
  Temp.v.y = v.x;
  Temp.p = v.x + c;
  if(nscal) for(int i=0; i<nscal; ++i) Temp.scalar[i] = v.x;

  for(int i=0; i<ns; ++i){
    Temp.spec[i].c = v.x;
  }
     
  return (Temp);
}


/************************************************************
 * Low Mach Number Preconditioned Eigenvalue(s)             *
 * LESPremixed2D_pState::lambda_preconditioned_x            *
 ************************************************************/
LESPremixed2D_pState LESPremixed2D_pState::lambda_preconditioned_x(const double &MR2) const {

  LESPremixed2D_pState NEW;
  double uprimed, cprimed;
  double c = amodified();
  u_a_precon(MR2*c*c,uprimed,cprimed);

  NEW.rho = uprimed - cprimed;
  NEW.v.x = v.x;
  NEW.v.y = v.x;
  NEW.p = uprimed + cprimed;
  if(nscal) for(int i=0; i<nscal; ++i) NEW.scalar[i] = v.x;

  for(int i=0; i<ns; ++i){
    NEW.spec[i].c = v.x;
  }
    
  return (NEW);
}

/*******************************************************************
 ***************** EIGENVECTORS ************************************
 *******************************************************************/
// Conserved Right Eigenvector -- (x-direction)
LESPremixed2D_cState LESPremixed2D_pState::rc_x(const int &index) const {

    if(index == 1){
      double c = amodified(); 
      return (LESPremixed2D_cState(ONE, v.x-c, v.y, H()/rho-v.x*c, scalar, spec));
    } else if(index == 2) {
      double c = amodified();
      return (LESPremixed2D_cState(ONE, v.x, v.y, H()/rho-c*c/(g()-ONE), scalar, spec)); 
    } else if(index == 3) {
      return (LESPremixed2D_cState(ZERO, ZERO, rho, rho*v.y, ZERO));
    } else if(index == 4) {
      double c = amodified();
      return (LESPremixed2D_cState(ONE, v.x+c, v.y, H()/rho+v.x*c, scalar, spec));
    } else if( nscal  &&  index >=5 && index<=(NUM_VAR_LESPREMIXED2D-ns)){
      for(int i=5; i<=(NUM_VAR_LESPREMIXED2D-ns); ++i){
        if(index == 5){
          LESPremixed2D_cState NEW(ZERO);
          NEW.E = FIVE*rho/THREE;   // For k equation
          NEW.rhoscalar[i-5] = rho; //rho*scalar[i-5];
          return NEW;
	} else {  
          LESPremixed2D_cState NEW(ZERO);
          NEW.rhoscalar[i-5] = rho*scalar[i-5]; // ????
          return NEW;
        }
      }
    } else { 
      LESPremixed2D_cState NEW(ZERO);
      double RTOT = Rtot();
      double TEMP = p/(rho*RTOT);     
      int count = index-(NUM_VAR_LESPREMIXED2D-ns+1);
#ifdef _NS_MINUS_ONE
      NEW.E = rho*((specdata[count].Enthalpy(TEMP) + specdata[count].Heatofform() - Cp(TEMP)*TEMP*specdata[count].Rs()/RTOT) -
		   (specdata[ns-1].Enthalpy(TEMP)  + specdata[ns-1].Heatofform() - Cp(TEMP)*TEMP*specdata[ns-1].Rs()/RTOT)); 
#else
      NEW.E = rho*(specdata[count].Enthalpy(TEMP) + specdata[count].Heatofform() - Cp(TEMP)*TEMP*specdata[count].Rs()/RTOT);      
#endif
      NEW.rhospec[count].c = rho;
      return NEW;
    }
}

// Primitive Left Eigenvector -- (x-direction)
LESPremixed2D_pState LESPremixed2D_pState::lp_x(const int &index) const {
 
   if(index == 1){
      double c = amodified(); 
      return (LESPremixed2D_pState(ZERO, -HALF*rho/c, ZERO, HALF/(c*c), ZERO));
   } else if(index == 2) {
      double c = amodified(); 
      return (LESPremixed2D_pState(ONE, ZERO, ZERO, -ONE/(c*c), ZERO));
   } else if(index == 3) {
      return  (LESPremixed2D_pState(ZERO, ZERO, ONE, ZERO, ZERO));
   } else if(index == 4) {  
      double c = amodified(); 
      return (LESPremixed2D_pState(ZERO, HALF*rho/c, ZERO, HALF/(c*c), ZERO));

   } else if(nscal  &&  index >=5  &&  index<=(NUM_VAR_LESPREMIXED2D-ns) ){
     for(int i=5; i<=(NUM_VAR_LESPREMIXED2D-ns); ++i){
       if(index == i){
	 LESPremixed2D_pState NEW(ZERO);
	 NEW.scalar[i-5] = ONE; 
	 return NEW;
       }
     }

   } else{ 
     LESPremixed2D_pState NEW(ZERO);
     NEW.spec[index-(NUM_VAR_LESPREMIXED2D-ns+1)].c = ONE;
     return NEW;
   } 

}

/************************************************************
 ************** PRECONDITIONED EIGENVECTORS *****************
 ************************************************************/
// Conserved Right Eigenvector -- (x-direction)
LESPremixed2D_cState LESPremixed2D_pState::rc_x_precon(const int &index, const double &MR2) const {

  if(index == 1){
    double c = amodified(); 
    double uprimed,cprimed;
    u_a_precon(MR2*c*c,uprimed,cprimed);
    return (LESPremixed2D_cState(ONE, 
				 (uprimed-cprimed)/MR2,
				 v.y,
				 h() + HALF*(v.y*v.y + v.x*v.x/MR2) + FIVE*k()/THREE - (v.x*cprimed)/MR2,
				 scalar, spec));
   
  } else if(index == 2) {
    double c = amodified();
    return (LESPremixed2D_cState(ONE, v.x, v.y, (H()/rho - c*c/(g()-ONE)), scalar,spec));
  
  } else if(index == 3) {
    return (LESPremixed2D_cState(ZERO, ZERO, rho, rho*v.y,ZERO));
   
  } else if(index == 4) { 
    double c = amodified(); 
    double uprimed,cprimed;
    u_a_precon(MR2*c*c,uprimed,cprimed);
    return (LESPremixed2D_cState(ONE,
				 (uprimed+cprimed)/MR2,
				 v.y, 
				 h() + HALF*(v.y*v.y + v.x*v.x/MR2) + FIVE*k()/THREE + (v.x*cprimed)/MR2,
				 scalar, spec));

  } else if(nscal  && index >=5  && index<=(NUM_VAR_LESPREMIXED2D-ns) ){
    for(int i=5; i<=(NUM_VAR_LESPREMIXED2D-ns); ++i){
      if(index == 5){
	LESPremixed2D_cState NEW(ZERO);
	NEW.E = FIVE*rho/THREE;   // For k equation
	NEW.rhoscalar[i-5] = rho; 
	return NEW;
      } else {  
	LESPremixed2D_cState NEW(ZERO);
	NEW.rhoscalar[i-5] = rho*scalar[i-5]; // ????
	return NEW;
      }        
    }  

  } else { 
    LESPremixed2D_cState NEW(ZERO);
    double RTOT = Rtot();
    double TEMP = p/(rho*RTOT);    
    int count = index-(NUM_VAR_LESPREMIXED2D-ns+1);
#ifdef _NS_MINUS_ONE
    NEW.E = rho*((specdata[count].Enthalpy(TEMP) + specdata[count].Heatofform() - Cp(TEMP)*TEMP*specdata[count].Rs()/RTOT) -
 		 (specdata[ns-1].Enthalpy(TEMP)  + specdata[ns-1].Heatofform()  - Cp(TEMP)*TEMP*specdata[ns-1].Rs()/RTOT)); 
#else
    NEW.E = rho*(specdata[count].Enthalpy(TEMP) + specdata[count].Heatofform() - Cp(TEMP)*TEMP*specdata[count].Rs()/RTOT);
#endif

    NEW.rhospec[count].c = rho;
    return NEW;    
  }

}

// Primitive Left Eigenvector -- (x-direction)
LESPremixed2D_pState LESPremixed2D_pState::lp_x_precon(const int &index, const double &MR2) const {
  
  if(index == 1){
    double c = amodified();   
    double uprimed,cprimed;
    u_a_precon(MR2*c*c,uprimed,cprimed);
    return (LESPremixed2D_pState(ZERO, 
				 -HALF*rho*MR2/cprimed, 
				 ZERO,
				 (-uprimed+cprimed + v.x)/(TWO*cprimed*c*c),
				 ZERO));
  } else if(index == 2) {
    double c = amodified(); 
    return (LESPremixed2D_pState(ONE, ZERO, ZERO, -ONE/(c*c), ZERO));
  } else if(index == 3) {
    return  (LESPremixed2D_pState(ZERO, ZERO, ONE, ZERO,ZERO));
  } else if(index == 4) {  
    double c = amodified(); 
    double uprimed,cprimed;
    u_a_precon(MR2*c*c,uprimed,cprimed);
    return (LESPremixed2D_pState(ZERO, 
				 HALF*rho*MR2/cprimed, 
				 ZERO,
				 (uprimed+cprimed - v.x)/(TWO*cprimed*c*c),
				 ZERO));

  } else if(nscal  &&  index >=5 && index<=(NUM_VAR_LESPREMIXED2D-ns) ){
    for(int i=5; i<=(NUM_VAR_LESPREMIXED2D-ns); ++i){
      if(index == i){
	LESPremixed2D_pState NEW(ZERO);
	NEW.scalar[i-5] = ONE; // scalar[i-5]; ?????
	return NEW;
      }
    }
  
  } else { 
    LESPremixed2D_pState NEW(ZERO);
    NEW.spec[index-(NUM_VAR_LESPREMIXED2D-ns+1)].c = ONE;
    return NEW;
  }
 
}

/*******************************************************************
 *******************************************************************
 *  Low Mach Number Preconditioner Based on Weiss & Smith (1995)   *
 *                                                                 *
 *******************************************************************
 *******************************************************************/
// For CFL Calculation
double LESPremixed2D_pState::u_plus_aprecon(const double &u, 
					    const int    &flow_type_flag, 
					    const double &deltax) const {
  double Temp = T();
  double c = amodified();
  double UR2 = Mr2(flow_type_flag,deltax)*c*c;
  double alpha = HALF*(ONE - ONE/(c*c) * UR2);
  
  // uprime + cprime
  return ( u*(ONE - alpha) + sqrt(alpha*alpha*u*u + UR2) );
}

/************************************************************/
// For eignvalues and eigenvectors return preconditioned
// velocity and soundspeed ie. u' and c'
void LESPremixed2D_pState::u_a_precon(const double &UR2, double &uprimed, double &cprimed) const{
  
  double Temp = T();
  double c = amodified();
  double alpha = HALF*(ONE - ONE/(c*c) * UR2);

  uprimed = v.x*(ONE - alpha);
  cprimed = sqrt(alpha*alpha*v.x*v.x + UR2); 

} 

/************************************************************/
// as defined by E.Turkel (1999)
double LESPremixed2D_pState::Mr2(const int    &flow_type_flag, 
				 const double &deltax) const {
  
  double c = amodified();
  double MR2 = min(max((v.sqr()/(c*c)),Mref*Mref),ONE);
  
  // Need deltax which is based on cell spacing 
  if(flow_type_flag){ 
    MR2 = pow(max(sqrt(MR2*c*c), mu()/(rho*deltax)),2.0)/(c*c);
  }

  return (MR2);
}


/*******************************************************************
  *******************************************************************
  *  Low Mach Number Preconditioner Based on Weiss & Smith (1995)   *
  *  - Includes the SFS turbulence kinetic energy                   *
  *******************************************************************
  *******************************************************************/
void LESPremixed2D_pState::Low_Mach_Number_Preconditioner(DenseMatrix &P,
							  const int &Viscous_flag, 
							  const double &deltax) const{
  
  double Temp = T();
  double pt = pmodified();
  double Rmix = Rtot();
  double enthalpy = h();
  double CP = Cp();
  double c = amodified();
  double theta = ONE/(Mr2(Viscous_flag,deltax)*c*c) + (g()-ONE)/(c*c);  
  double kk = k();
 
  double phi = ZERO;   
  for(int j=0; j<ns-1; ++j){   
#ifdef _NS_MINUS_ONE
    phi += spec[j].c*(specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/Rmix -    
		      (specdata[ns-1].Enthalpy(Temp) + specdata[ns-1].Heatofform() - CP*Temp*specdata[ns-1].Rs()/Rmix));
#else
    phi += spec[j].c*(specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() 
		      - CP*Temp*specdata[j].Rs()/Rmix);   
#endif
  }

  double V = HALF*v.sqr();
  double alpha = theta*pt/rho;
  double alpham1 = alpha - ONE;
  double omega = (Rmix - CP)*pt/(rho*Rmix);
  double beta = (Scal_sys.scalar_flag == LES_TF_K)  ?  (enthalpy - CP*pt/(rho*Rmix) - phi)
    : (enthalpy - CP*pt/(rho*Rmix) - phi + 5.0*kk/3.0);
  
  P.zero(); //RESET MATRIX TO ZERO!!!

  if(Scal_sys.scalar_flag == LES_TF_K) {
    P(0,0) = (alpha*(beta-V) + V + pt/rho - enthalpy + phi)/omega; 
  } else {
    P(0,0) = (alpha*(beta-V) + V + pt/rho - enthalpy + phi - 5.0*kk/3.0)/omega; 
  }
  P(0,1) = v.x*alpham1/omega;                              
  P(0,2) = v.y*alpham1/omega;                           
  P(0,3) = -alpham1/omega;                              
  P(1,0) = v.x*(beta-V)*alpham1/omega;                  
  P(1,1) = v.x*v.x*alpham1/omega + 1.0;                 
  P(1,2) = v.x*v.y*alpham1/omega;                       
  P(1,3) = -v.x*alpham1/omega;                            
  P(2,0) = v.y*(beta-V)*alpham1/omega;                     
  P(2,1) = v.x*v.y*alpham1/omega;                          
  P(2,2) = v.y*v.y*alpham1/omega + 1.0;                    
  P(2,3) = -v.y*alpham1/omega;                             
  P(3,0) = (enthalpy+V+5.0/3.0*kk)*(beta-V)*alpham1/omega; 
  P(3,1) = v.x*(enthalpy+V+5.0/3.0*kk)*alpham1/omega;      
  P(3,2) = v.y*(enthalpy+V+5.0/3.0*kk)*alpham1/omega;      
  if(Scal_sys.scalar_flag == LES_TF_K) {
    P(3,3) = -(alpha*(enthalpy+V+5.0/3.0*kk) - V - pt/rho - beta - phi - 5.0/3.0*kk)/omega;
  } else {
    P(3,3) = -(alpha*(enthalpy+V+5.0/3.0*kk) - V - pt/rho - beta - phi)/omega;
  }
      
  int NUM_VAR = NUM_LESPREMIXED2D_VAR_SANS_SPECIES + nscal;

  // Scalar (SFS turbulence kinetic energy)
  if (nscal  &&  Scal_sys.scalar_flag == LES_TF_K) {
    P(0,4) = 5.0/3.0*alpham1/omega;                   
    P(1,4) = 5.0/3.0*v.x*alpham1/omega;               
    P(2,4) = 5.0/3.0*v.y*alpham1/omega;               
    P(3,4) = 5.0/3.0*(enthalpy+V+5.0/3.0*kk)*alpham1/omega;   

    P(4,0) = kk*(beta-V)*alpham1/omega;   
    P(4,1) = kk*v.x*alpham1/omega;        
    P(4,2) = kk*v.y*alpham1/omega;        
    P(4,3) = -kk*alpham1/omega;           
    P(4,4) = 5.0/3.0*kk*alpham1/omega + 1.0;

    for (int j=0; j<ns-1; ++j) {
#ifdef _NS_MINUS_ONE     
      double enth_j = specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/Rmix-
	( specdata[ns-1].Enthalpy(Temp) + specdata[ns-1].Heatofform() - CP*Temp*specdata[ns-1].Rs()/Rmix);
#else
      double enth_j = specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/Rmix;
#endif

      P(NUM_VAR-1, j+NUM_VAR) = kk*enth_j*alpham1/omega;  
      P(j+NUM_VAR, NUM_VAR-1) = 5.0/3.0*(spec[j].c)*alpham1/omega;   
    }
  }
  

  // Multispecies
  for (int j=0; j<ns-1; ++j) {
#ifdef _NS_MINUS_ONE     
    double enth_j = specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/Rmix-
      ( specdata[ns-1].Enthalpy(Temp) + specdata[ns-1].Heatofform() - CP*Temp*specdata[ns-1].Rs()/Rmix);
#else
    double enth_j = specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/Rmix;
#endif       
    P(0,j+NUM_VAR) = enth_j*alpham1/omega;       
    P(1,j+NUM_VAR) = v.x*enth_j*alpham1/omega;   
    P(2,j+NUM_VAR) = v.y*enth_j*alpham1/omega;   
    P(3,j+NUM_VAR) = enth_j*(enthalpy+V+5.0/3.0*kk)*alpham1/omega; 	
    for (int i=0; i<ns-1; ++i) { 
      if (i==j) { 
	P(i+NUM_VAR,0) = (spec[i].c)*(beta-V)*alpham1/omega;   
	P(i+NUM_VAR,1) = (spec[i].c)*v.x*alpham1/omega;        
	P(i+NUM_VAR,2) = (spec[i].c)*v.y*alpham1/omega;        
	P(i+NUM_VAR,3) = -(spec[i].c)*alpham1/omega;           
	//diagonal
	P(i+NUM_VAR,j+NUM_VAR) = (spec[i].c)*enth_j*alpham1/omega + ONE;   
      } else {
	P(i+NUM_VAR,j+NUM_VAR) = (spec[i].c)*enth_j*alpham1/omega;  
      }
    }       
  }

} // end  Low_Mach_Number_Preconditioner


/*******************************************************************
 *******************************************************************
 *  Inverse of the Low Mach Number Preconditioner                  *
 *  - Includes the SFS turbulence kinetic energy                   *
 *******************************************************************
 *******************************************************************/
void LESPremixed2D_pState::Low_Mach_Number_Preconditioner_Inverse(DenseMatrix &Pinv,	
								  const int &Viscous_flag, 
								  const double &deltax) const{
  
  double Temp = T();
  double pt = pmodified();
  double Rmix = Rtot();
  double enthalpy = h();
  double CP = Cp();
  double c = amodified();  
  double theta = ONE/(Mr2(Viscous_flag,deltax)*c*c) + (g()-ONE)/(c*c);  
  double kk = k();

  double phi = ZERO;
  for(int j=0; j<ns-1; ++j){ 
#ifdef _NS_MINUS_ONE       
    phi += spec[j].c*(specdata[j].Enthalpy(Temp)+ specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/Rmix
		      - (specdata[ns-1].Enthalpy(Temp)+ specdata[ns-1].Heatofform() - CP*Temp*specdata[ns-1].Rs()/Rmix));
#else
    phi += spec[j].c*(specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/Rmix);   
#endif
  }
  
  double AA = pt*(rho*Rmix-theta*pt*CP);
  double BB = Rmix*rho*(theta*pt-rho);
  double EE = (Scal_sys.scalar_flag == LES_TF_K)  ?  (HALF*v.sqr() - enthalpy + phi)
    : (HALF*v.sqr() - enthalpy + phi - 5.0*kk/3.0);
      
  double CC = EE + CP*pt/(rho*Rmix); 
  double DD = HALF*v.sqr() + enthalpy + 5.0*kk/3.0;

  Pinv.zero();  //RESET MATRIX TO ZERO!!!

  Pinv(0,0) = rho*Rmix/AA*(theta*pt*EE-rho*CC+pt);  
  Pinv(0,1) = -v.x*BB/AA;                           
  Pinv(0,2) = -v.y*BB/AA;                           
  Pinv(0,3) = BB/AA;                                
  Pinv(1,0) = v.x*CC*BB/AA;                         
  Pinv(1,1) = rho*Rmix/AA*(pt+rho*v.x*v.x-theta*pt*(v.x*v.x+CP*pt/(rho*Rmix)));  
  Pinv(1,2) = -v.x*v.y*BB/AA;                      
  Pinv(1,3) = v.x*BB/AA;                           
  Pinv(2,0) = v.y*CC*BB/AA;                        
  Pinv(2,1) = -v.x*v.y*BB/AA;                      
  Pinv(2,2) = rho*Rmix/AA*(pt+v.y*v.y*rho-theta*pt*(v.y*v.y+CP*pt/(rho*Rmix))); 
  Pinv(2,3) = v.y*BB/AA;                           
  Pinv(3,0) = DD*CC*BB/AA;                         
  Pinv(3,1) = -v.x*DD*BB/AA;                       
  Pinv(3,2) = -v.y*DD*BB/AA;                       
  Pinv(3,3) = rho*Rmix/AA*(theta*pt*(DD-CP*pt/(rho*Rmix))-rho*DD+pt); 

  int NUM_VAR = NUM_LESPREMIXED2D_VAR_SANS_SPECIES + nscal;

  // Scalar (SFS turbulence kinetic energy)
  if (nscal  &&  Scal_sys.scalar_flag == LES_TF_K) {
    Pinv(0,4) = -5.0/3.0*BB/AA;         
    Pinv(1,4) = -5.0/3.0*v.x*BB/AA;     
    Pinv(2,4) = -5.0/3.0*v.y*BB/AA;     
    Pinv(3,4) = -5.0/3.0*DD*BB/AA;      

    Pinv(4,0) = kk*CC*BB/AA;            
    Pinv(4,1) = -kk*v.x*BB/AA;          
    Pinv(4,2) = -kk*v.y*BB/AA;          
    Pinv(4,3) = kk*BB/AA;               
    Pinv(4,4) = 1.0 - 5.0/3.0*kk*BB/AA; 

    for (int j=0; j<ns-1; ++j) {
#ifdef _NS_MINUS_ONE    
      double enth_j =  specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/Rmix -
	(specdata[ns-1].Enthalpy(Temp) + specdata[ns-1].Heatofform() - CP*Temp*specdata[ns-1].Rs()/Rmix);
#else
      double enth_j =  specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/Rmix;
#endif

      Pinv(NUM_VAR-1, j+NUM_VAR) = -kk*enth_j*BB/AA;   
      Pinv(j+NUM_VAR, NUM_VAR-1) = -5.0/3.0*(spec[j].c)*BB/AA; 

    }    
  }


  // Multispecies
  for (int j=0; j<ns-1; ++j) {
#ifdef _NS_MINUS_ONE    
    double enth_j =  specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/Rmix -
      (specdata[ns-1].Enthalpy(Temp) + specdata[ns-1].Heatofform() - CP*Temp*specdata[ns-1].Rs()/Rmix);
#else
    double enth_j =  specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/Rmix;
#endif
    Pinv(0,j+NUM_VAR) = -enth_j*BB/AA;         
    Pinv(1,j+NUM_VAR) = -v.x*enth_j*BB/AA;     
    Pinv(2,j+NUM_VAR) = -v.y*enth_j*BB/AA;     
    Pinv(3,j+NUM_VAR) = -enth_j*BB*DD/AA;      	
    for (int i=0; i<ns-1; ++i) { 
      if (i==j) { 
	Pinv(i+NUM_VAR,0) = (spec[i].c)*CC*BB/AA;        
	Pinv(i+NUM_VAR,1) = -(spec[i].c)*v.x*BB/AA;    
	Pinv(i+NUM_VAR,2) = -(spec[i].c)*v.y*BB/AA;    
	Pinv(i+NUM_VAR,3) = (spec[i].c)*BB/AA;            
	//diagonal
	Pinv(i+NUM_VAR,j+NUM_VAR) = ONE - spec[i].c*enth_j*BB/AA;    
      } else {
	Pinv(i+NUM_VAR,j+NUM_VAR) = -spec[i].c*enth_j*BB/AA;      
      }
    }       
  }

} // end Low_Mach_Number_Preconditioner_Inverse



/***********************************************************
 *************** Operator Overloading **********************
 ***********************************************************/
 //Assuming # of species are equal to avoid overhead of checking
 //ie.  if(ns == W.ns)

/********************************************************
 * LESPremixed2D_pState -- Binary arithmetic operators. *
 ********************************************************/
//----------------- Addition -----------------------------//
LESPremixed2D_pState LESPremixed2D_pState::operator +(const LESPremixed2D_pState &W) const{    
  LESPremixed2D_pState Temp(*this);
  Temp += W;
  return Temp;
}

//------------------ Subtraction ------------------------//
LESPremixed2D_pState LESPremixed2D_pState::operator -(const LESPremixed2D_pState &W) const{
  LESPremixed2D_pState Temp(*this);
  Temp -= W;
  return Temp;
}

//---------------- Scalar Multiplication ------------------//
LESPremixed2D_pState LESPremixed2D_pState::operator *(const double &a) const{
  LESPremixed2D_pState Temp(*this);
  Temp.rho = rho*a;  Temp.v = v*a; Temp.p = p*a;
  if(nscal) for(int i=0; i<nscal; ++i) Temp.scalar[i] = scalar[i]*a;
  for( int i=0; i<ns; ++i) Temp.spec[i] = spec[i]*a; 
  Temp.tau = tau*a;
  Temp.qflux = qflux*a;
  Temp.lambda= lambda*a;
  Temp.theta = theta*a;
#ifdef THICKENED_FLAME_ON
  Temp.flame = flame*a; 
#endif
  return(Temp);
}

LESPremixed2D_pState operator *(const double &a, const LESPremixed2D_pState &W){
  LESPremixed2D_pState Temp;
  Temp.rho = W.rho*a;  Temp.v = W.v*a; Temp.p = W.p*a;
  if(W.nscal) for(int i=0; i<W.nscal; ++i) Temp.scalar[i] = W.scalar[i]*a; 
  for( int i=0; i<W.ns; ++i) Temp.spec[i] = W.spec[i]*a;
  Temp.tau = W.tau*a;
  Temp.qflux = W.qflux*a;
  Temp.lambda= W.lambda*a;
  Temp.theta = W.theta*a;
#ifdef THICKENED_FLAME_ON
  Temp.flame = W.flame*a;
#endif
  return(Temp);
}

//--------------- Scalar Division ------------------------//
LESPremixed2D_pState LESPremixed2D_pState::operator /(const double &a) const {
  LESPremixed2D_pState Temp(*this);
  Temp.rho = rho/a; Temp.v = v/a; Temp.p = p/a;
  if(nscal) for(int i=0; i<nscal; ++i) Temp.scalar[i] = scalar[i]/a;
  for(int i=0; i<ns; ++i) Temp.spec[i] = spec[i]/a;
  Temp.tau = tau/a;
  Temp.qflux = qflux/a;
  Temp.lambda= lambda/a;
  Temp.theta = theta/a;
#ifdef THICKENED_FLAME_ON
  Temp.flame = flame/a;
#endif
  return(Temp);
}

//----------------- Inner Product ------------------------//
double LESPremixed2D_pState::operator *(const LESPremixed2D_pState &W) const{
  double sum=0.0;
  double sumscal=0.0;
  if(nscal) for(int i=0; i<nscal; ++i) sumscal += scalar[i]*W.scalar[i];
  for(int i=0; i<ns; ++i) sum += spec[i]*W.spec[i];
  return (rho*W.rho + v*W.v + p*W.p + sumscal + sum);
}

//----------- solution state product operator ------------//
LESPremixed2D_pState LESPremixed2D_pState::operator ^( const LESPremixed2D_pState &W) const {
    LESPremixed2D_pState Temp(*this);
    Temp.rho = rho*W.rho;
    Temp.v.x = v.x*W.v.x;
    Temp.v.y = v.y*W.v.y;
    Temp.p = p*W.p;
    if(nscal) for(int i=0; i<nscal; ++i) Temp.scalar[i] = scalar[i]*W.scalar[i];
    for(int i=0; i<ns; ++i) Temp.spec[i] = spec[i]*W.spec[i];
#ifdef THICKENED_FLAME_ON
  Temp.flame.WF = flame.WF*W.flame.WF;
  Temp.flame.TF = flame.TF*W.flame.TF;
#endif
    return(Temp);
}

//----------------- Assignment ----------------------------//
LESPremixed2D_pState& LESPremixed2D_pState::operator =(const LESPremixed2D_pState &W){
  //self assignment protection
  if( this != &W){ 
    Copy(W);
  }
  return (*this);
}

/**********************************************************
 * LESPremixed2D_pState -- Shortcut arithmetic operators. *
 **********************************************************/
LESPremixed2D_pState& LESPremixed2D_pState::operator +=(const LESPremixed2D_pState &W){
  rho += W.rho;
  v += W.v; 
  p += W.p; 
  if(nscal) for(int i=0; i<nscal; ++i) scalar[i] += W.scalar[i];
  for( int i=0; i<ns; ++i)  spec[i].c += W.spec[i].c;
  tau += W.tau;
  qflux += W.qflux;
  lambda += W.lambda;
  theta += W.theta;
#ifdef THICKENED_FLAME_ON
  flame += W.flame;
#endif
  return (*this);
}

LESPremixed2D_pState& LESPremixed2D_pState::operator -=(const LESPremixed2D_pState &W) {
  rho -= W.rho;
  v -= W.v;
  p -= W.p;
  if(nscal) for(int i=0; i<nscal; ++i) scalar[i] -= W.scalar[i];
  for(int i=0; i<ns; ++i) spec[i].c -= W.spec[i].c;
  tau -= W.tau;
  qflux -= W.qflux;
  lambda -= W.lambda;
  theta -= W.theta;
#ifdef THICKENED_FLAME_ON
  flame -= W.flame;
#endif
  return (*this); 
}

/********************************************************
 * LESPremixed2D_pState -- Unary arithmetic operators.  *
 ********************************************************/
//  LESPremixed2D_pState operator +(const LESPremixed2D_pState &W) {  
//   return (LESPremixed2D_pState(W.rho,W.v,W.p,W.spec));
// }

LESPremixed2D_pState operator -(const LESPremixed2D_pState &W) {
#ifdef STATIC_LESPREMIXED2D_SPECIES
  Species spt[STATIC_LESPREMIXED2D_SPECIES];
#else
  Species *spt= new Species[W.ns];
#endif

  double *scalpt = new double[W.nscal];
  for(int i=0; i<W.nscal; ++i) scalpt[i] = -W.scalar[i];

  for(int i=0; i<W.ns; ++i)  spt[i] = -W.spec[i]; 
  LESPremixed2D_pState Temp(-W.rho,-W.v,-W.p, scalpt, spt);
  Temp.tau = -W.tau;
  Temp.qflux = -W.qflux;
  Temp.lambda = -W.lambda;
  Temp.theta= -W.theta;
#ifdef THICKENED_FLAME_ON
  Temp.flame = -W.flame;
#endif

#ifndef STATIC_LESPREMIXED2D_SPECIES
  delete[] spt;
#endif

  delete[] scalpt;

  return(Temp);
}

/********************************************************
 * LESPremixed2D_pState -- Relational operators.        *
 ********************************************************/
int operator ==(const LESPremixed2D_pState &W1, const LESPremixed2D_pState &W2) {
  if(W1.ns == W2.ns){ //check that species are equal
    bool Temp, Temp2;
    //species
    for(int i=0; i<W1.ns; ++i){
      if( W1.spec[i] == W2.spec[i]){
	Temp = true;
      } else {
	Temp = false;
	break;
      }
    }
    //scalars
    if(W1.nscal){        
      for(int i=0; i<W1.nscal; ++i){
	if(W1.scalar[i] == W2.scalar[i]){
	  Temp2 = true;
	} else {
	  Temp2 = false;
	  break;
	}
      }
    }  
 
    return (W1.rho == W2.rho && W1.v == W2.v && W1.p == W2.p 
	    && Temp2 == true && Temp == true  && W1.tau == W2.tau &&
	    W1.qflux == W2.qflux&& W1.lambda == W2.lambda &&
	    W1.theta == W2.theta);
    
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
}

int operator !=(const LESPremixed2D_pState &W1, const LESPremixed2D_pState &W2) {
   if(W1.ns == W2.ns){ //check that species are equal
     bool Temp = true, Temp2 = true;
     //species
     for(int i=0; i<W1.ns; ++i){
       if( W1.spec[i] != W2.spec[i] ){
	 Temp = false;
	 break;
       }
     }    
     //scalars
     for(int i=0; i<W1.nscal; ++i){
       if( W1.scalar[i] != W2.scalar[i] ){
	 Temp2 = false;
	 break;
       }
     }
     
     return (W1.rho != W2.rho || W1.v != W2.v || W1.p != W2.p 
	     || Temp2 != true || Temp != true || W1.tau != W2.tau ||
	     W1.qflux != W2.qflux|| W1.lambda != W2.lambda ||
	     W1.theta != W2.theta);
     
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
}

/********************************************************
 * LESPremixed2D_pState -- Input-output operators.      *
 ********************************************************/
ostream &operator << (ostream &out_file, const LESPremixed2D_pState &W) {
  out_file.precision(10);
  out_file.setf(ios::scientific);
  out_file << " " << W.rho  << " " << W.v.x << " " << W.v.y << " " << W.p;

  if(W.nscal) for(int i=0; i<W.nscal; ++i) out_file <<" "<<W.scalar[i];

  for( int i=0; i<W.ns; ++i) out_file<<" "<<W.spec[i];
 
#ifdef THICKENED_FLAME_ON
  out_file << " " << W.flame;
#endif

  out_file.unsetf(ios::scientific);
  return (out_file);
}

istream &operator >> (istream &in_file, LESPremixed2D_pState &W) {
  in_file.setf(ios::skipws);
  in_file >> W.rho >> W.v.x >> W.v.y >> W.p;

  if(W.nscal) for(int i=0; i<W.nscal; ++i) in_file >> W.scalar[i];

  for( int i=0; i<W.ns; ++i) in_file>>W.spec[i];

#ifdef THICKENED_FLAME_ON
  in_file >> W.flame;
#endif

  in_file.unsetf(ios::skipws);
   return (in_file);
}


/***************************************************************
 ***************************************************************
 ********** SOURCE TERMS & JACOBIANS ***************************
 ***************************************************************
 ***************************************************************/

/***************************************************************
 * Axisymmetric flow source terms (Inviscid)                   *
 ***************************************************************/
LESPremixed2D_cState LESPremixed2D_pState::Sa_inviscid(const Vector2D &X,
						       const int Flow_Type,
						       const int Axisymmetric) const{

  LESPremixed2D_cState Temp; Temp.Vacuum();

  //x is radial
  if (Axisymmetric == AXISYMMETRIC_X) {
     Temp.rho = -rho*v.x/X.x;
     Temp.rhov.x = -rho*v.x*v.x/X.x; 
     Temp.rhov.y = -rho*v.x*v.y/X.x;
     Temp.E =  -v.x*H()/X.x;
     //species contributions
     for(int i=0; i<ns; ++i){
       Temp.rhospec[i].c = -rho*v.x*spec[i].c/X.x;         //correct for ns-1 ????
     }
     //scalars contributions
     if(nscal) for(int i=0; i<nscal; ++i) Temp.rhoscalar[i] =  -v.x*rho*scalar[i]/X.x;
 

     //y is radial
  } else if (Axisymmetric == AXISYMMETRIC_Y) {
     Temp.rho = -rho*v.y/X.y;
     Temp.rhov.x = -rho*v.x*v.y/X.y; 
     Temp.rhov.y = -rho*v.y*v.y/X.y;
     Temp.E =  -v.y*H()/X.y;
     //species contributions
     for(int i=0; i<ns; ++i){
       Temp.rhospec[i].c = -rho*v.y*spec[i].c/X.y;
     }
     //scalars contributions
     if(nscal) for(int i=0; i<nscal; ++i) Temp.rhoscalar[i] =  -v.y*rho*scalar[i]/X.y;   
  }

  return (Temp);

}


// Finite difference check of dSa_idU
void LESPremixed2D_pState::dSa_idU_FD(DenseMatrix &dSa_IdU, 
				      const Vector2D &X, 
				      const int Flow_Type,
				      const int Axisymmetric ) const {

  LESPremixed2D_cState UU = U(*this);
  LESPremixed2D_cState A,C;
  LESPremixed2D_pState B,D;
  double perturb = 5e-6;
  double a;

  for(int jcol=0; jcol<(dSa_IdU.get_n()); ++jcol){    
    A = UU; C = UU;
    if( jcol <NUM_LESPREMIXED2D_VAR_SANS_SPECIES) {
      A[jcol+1] += perturb*max(ONE,UU[jcol+1]); 
      C[jcol+1] -= perturb*max(ONE,UU[jcol+1]);
    } else {                                       //enforce sum(ci) = 1;
      a =  perturb*max(ONE,UU[jcol+1]); 
      A[jcol+1] += a;
      A[NUM_VAR_LESPREMIXED2D] -= a;      
      C[jcol+1] -= a;
      C[NUM_VAR_LESPREMIXED2D] += a;
    }   
    B = W(A);  D = W(C);
    A = B.Sa_inviscid(X,Flow_Type,Axisymmetric);   C = D.Sa_inviscid(X,Flow_Type,Axisymmetric);
    for(int irow=0; irow<(dSa_IdU.get_n()); ++irow){
      dSa_IdU(irow,jcol) = ( A[irow+1] - C[irow+1])/(TWO*perturb*max(ONE, UU[jcol+1]));      
    }
  } 

}

/****************************************************************
 * Axisymmetric Source Term Jacobian (Inviscid)                 * 
 ****************************************************************/
void LESPremixed2D_pState::dSa_idU(DenseMatrix &dSa_IdU, 
				   const Vector2D &X, 
				   const int Flow_Type,
				   const int Axisymmetric ) const {

  double enthalpy = h();
  double CP = Cp();
  double RTOT = Rtot();
  double phi = ZERO;
  double Temp = p/(rho*RTOT);

  for(int j=0; j<ns-1; ++j){   
#ifdef _NS_MINUS_ONE
    phi += spec[j].c*(specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/RTOT
		      -(specdata[ns-1].Enthalpy(Temp)+specdata[ns-1].Heatofform() - CP*Temp*specdata[ns-1].Rs()/RTOT)); 		      
#else
    phi += spec[j].c*(specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/RTOT);
#endif
  } 
 
  if(Axisymmetric == AXISYMMETRIC_Y){
    cout<<"\n ISSUES IN dSa_idU AXISYMMETRIC_Y \n"; exit(1);
    //end of axisymmetric case 1  -- y radial direction 
  } else if(Axisymmetric == AXISYMMETRIC_X){ 
    
    dSa_IdU(0,1) -= ONE/X.x;
    dSa_IdU(1,0) += v.x*v.x/X.x;
    dSa_IdU(1,1) -= TWO*v.x/X.x;
    dSa_IdU(2,0) += v.x*v.y/X.x;
    dSa_IdU(2,1) -= v.y/X.x;
    dSa_IdU(2,2) -= v.x/X.x;
   
//     if((Flow_Type ==FLOWTYPE_INVISCID)||(Flow_Type ==FLOWTYPE_LAMINAR)){    
//     dSa_IdU(3,0) -= v.x*( - CP*(HALF*(v.x*v.x+v.y*v.y)-enthalpy+phi-CP*Temp)/(RTOT-CP) - (phi +v.x*v.x +v.y*v.y + CP*Temp))/X.x;
//     dSa_IdU(3,1) -= (v.x*v.x + H()/rho + CP*v.x*v.x/(RTOT-CP))/X.x;
//     dSa_IdU(3,2) -= (v.x*v.y + CP*v.x*v.y/(RTOT-CP))/X.x ;    
//     dSa_IdU(3,3) -= CP/(CP-RTOT)*v.x/X.x;     

    // XINFENGS ORIGINAL
    dSa_IdU(3,0) -= RTOT/(CP-RTOT)*(v.x*(phi + v.sqr()) + v.x*CP*( p/rho - enthalpy - HALF*v.sqr())/RTOT)/X.x;
    dSa_IdU(3,1) -= (enthalpy + CP/(CP-RTOT)*HALF*v.sqr() - RTOT/(CP-RTOT)*HALF*(THREE*v.x*v.x +v.y*v.y))/X.x;
    dSa_IdU(3,2) += RTOT/(CP-RTOT)*v.x*v.y/X.x;
    dSa_IdU(3,3) -= CP/(CP-RTOT)*v.x/X.x;      
    //}//end of laminar case
    
    //Multispecies terms
    int NUM_VAR = NUM_LESPREMIXED2D_VAR_SANS_SPECIES + nscal;
    for(int i=0; i<(ns-1); ++i){
#ifdef _NS_MINUS_ONE     
      dSa_IdU(3,i+NUM_VAR) += v.x*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform() - CP*Temp*specdata[i].Rs()/RTOT
				   -(specdata[ns-1].Enthalpy(Temp)+specdata[ns-1].Heatofform() - CP*Temp*specdata[ns-1].Rs()/RTOT))/(CP/RTOT - ONE)/X.x; 
#else
      dSa_IdU(3,i+NUM_VAR) += v.x*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform() - CP*Temp*specdata[i].Rs()/RTOT)/(CP/RTOT - ONE)/X.x; 
#endif
      dSa_IdU(NUM_VAR+i,0) += v.x*spec[i].c/X.x;
      dSa_IdU(NUM_VAR+i,1) -= spec[i].c/X.x;
      dSa_IdU(NUM_VAR+i,NUM_VAR+i) -= v.x/X.x;
    }
    

    if(Flow_Type == FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY ||
       Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {
      cout <<"\nYOU ARE NOT SUPPOSED TO REACH THIS LINE, dSa_idU";
    }

  }//end of axisymmetric case   -- x radial 
}

/***************************************************************
 * Axisymmetric flow source terms (Viscous)                    *  
 ***************************************************************/
LESPremixed2D_cState LESPremixed2D_pState::Sa_viscous(const LESPremixed2D_pState &dWdx,
						      const LESPremixed2D_pState &dWdy,
						      const Vector2D &X,
						      const int Flow_Type,
						      const int Axisymmetric){
  double Mu, Temperature, Rmix;
  double mut, Dm_t;
  double rhohsDs;
  Vector2D grad_T;
  LESPremixed2D_cState Temp; Temp.Vacuum();

  //Transport and thermodynamic properties
  Mu = mu();
  Temperature = T();
  Rmix = Rtot();

  //Temperature gradient
  //dT/dx = 1/rho*R *( dP/dx - P/rho * drho/dx)
  grad_T.x = (ONE/(rho*Rmix))*(dWdx.p - (p/rho)*dWdx.rho);
  grad_T.y = (ONE/(rho*Rmix))*(dWdy.p - (p/rho)*dWdy.rho);

  
  //Molecular (laminar) stress tensor
  Tensor2D strain_rate;
  strain_rate = Strain_Rate(dWdx,dWdy, Flow_Type, Axisymmetric, X);  
  Laminar_Stress(strain_rate);


  //Molecular (laminar) heat flux
  //Thermal conduction, q = - kappa * grad(T)
  qflux = - kappa()*grad_T;
  //Thermal diffusion, q -= rho * sum ( hs * Ds *gradcs)
  for(int i=0; i<ns; ++i){ 
    rhohsDs = rho*(Mu/(rho*Schmidt[i]))*(specdata[i].Enthalpy(Temperature) + specdata[i].Heatofform());
    qflux.x -= rhohsDs*dWdx.spec[i].c;
    qflux.y -= rhohsDs*dWdy.spec[i].c;
  }
  
  //Turbulent heat flux
  //Thermal conduction, q = - kappa * grad(T)
  if(Flow_Type == FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY ||
     Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {

    //Turbulence model eddy viscosity
    mut = mu_t(strain_rate); 
    Dm_t = Dm_turb(mut);

    theta = - mut*Cp()/Pr_turb()*grad_T;
    //Thermal Diffusion, q -= rho * sum ( hs * Ds *gradcs)  
    for (int i=0; i<ns; ++i) {
      rhohsDs = rho*Dm_t*(/*specdata[i].Enthalpy(Temperature)+*/ specdata[i].Heatofform());
      theta.x -= rhohsDs*dWdx.spec[i].c;
      theta.y -= rhohsDs*dWdy.spec[i].c;
    }

    SFS_Stress(strain_rate);
  
  } 

  /////////////////////////////////////////////////////////////
  //Determine axisymmetric source terms
  /////////////////////////////////////////////////////////////
  if (Axisymmetric == AXISYMMETRIC_X) {
    Temp.rhov.x = (tau.xx - tau.zz)/X.x;
    Temp.rhov.y = tau.xy/X.x;
    Temp.E = (- qflux.x + v.x*tau.xx + v.y*tau.xy)/X.x;
    for(int i=0; i<ns; ++i){
      Temp.rhospec[i].c = rho*(Mu/(rho*Schmidt[i]))*dWdx.spec[i].c/X.x;
    }
    
    if(Flow_Type == FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY ||
       Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {

      Temp.rhov.x += (lambda.xx - lambda.zz)/X.x;
      Temp.rhov.y += lambda.xy/X.x;
      Temp.E += (- theta.x + v.x*lambda.xx + v.y*lambda.xy)/X.x; // +(Mu + mu_t)*dWdx.k/X.x; ??????

      if (nscal) for(int i=0; i<nscal; ++i) Temp.rhoscalar[i] = (Mu + mut)*dWdx.scalar[i]/X.x;
      for(int i=0; i<ns; ++i){
	Temp.rhospec[i].c = rho*Dm_t*dWdx.spec[i].c/X.x;
      }
    }
    
  } else if (Axisymmetric == AXISYMMETRIC_Y) {
    Temp.rhov.x = tau.xy/X.y;
    Temp.rhov.y = (tau.xx - tau.zz)/X.y;
    Temp.E = (- qflux.y + v.x*tau.xy + v.y*tau.yy)/X.y;
    for(int i=0; i<ns; ++i){
      Temp.rhospec[i].c = rho*(Mu/(rho*Schmidt[i]))*dWdy.spec[i].c/X.y;
    }
    
    if(Flow_Type == FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY ||
       Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {
      Temp.rhov.x += lambda.xy/X.y;
      Temp.rhov.y += (lambda.xx - lambda.zz)/X.y;
      Temp.E += (- theta.y + v.x*lambda.xy + v.y*lambda.yy)/X.y; //+(Mu + mu_t*sigma_star)*dWdy.k/X.y; ?????

      if (nscal) for(int i=0; i<nscal; ++i) Temp.rhoscalar[i] = (Mu + mut)*dWdy.scalar[i]/X.y;
      for(int i=0; i<ns; ++i){
	Temp.rhospec[i].c = rho*Dm_t*dWdy.spec[i].c/X.y;
      }
    }
    
  } /* endif */
  
  return (Temp);
  
}

/**************************************************************** 
 * Axisymmetric Source Term Jacobian (Viscous)                  *  
 ****************************************************************/
 void LESPremixed2D_pState::dSa_vdW(DenseMatrix &dSa_VdW,
				    const LESPremixed2D_pState &dWdx,
				    const LESPremixed2D_pState &dWdy,const Vector2D &X, 
				    const int Flow_Type,const int Axisymmetric,
				    const double d_dWdx_dW, const double d_dWdy_dW) const {

  int NUM_VAR = NUM_LESPREMIXED2D_VAR_SANS_SPECIES;
  double Rmix = Rtot();
  double Temp = T();
  double Mu = mu();
  double Kappa = kappa();
  
  if(Axisymmetric == AXISYMMETRIC_X){   
    double radius = X.x;

    dSa_VdW(2,2) += TWO*Mu*(d_dWdx_dW - ONE/radius)/radius;    
    dSa_VdW(2,1) += Mu*d_dWdy_dW/radius;
    dSa_VdW(2,2) += Mu*d_dWdx_dW/radius ;
    
    //energy
    double Sum_q(ZERO), Sum_dq(ZERO), Sum_dhi(ZERO);
    for(int Num = 0; Num<ns; ++Num){
      Sum_q += specdata[Num].Enthalpy_prime(Temp)*Mu*dWdx.spec[Num].c/(Schmidt[Num]*rho*Rmix);      // - ns-1 ???
      Sum_dq -= specdata[Num].Enthalpy_prime(Temp)*Mu*dWdx.spec[Num].c*Temp/(rho*Schmidt[Num]);
      //Sum_dhi += specdata[Num].Enthalpy_prime(Temp)*Temp*Mu*dWdx.spec[Num].c/(Rmix*Schmidt[Num]); 
    } 

    dSa_VdW(3,0) += (Kappa*( -(dWdx.p - p/rho*dWdx.rho)/rho + (p*dWdx.rho/(rho*rho) - p*d_dWdx_dW/rho))/(rho*Rmix) + Sum_dq)/radius;
    dSa_VdW(3,1) += (Mu*v.y*d_dWdy_dW + TWO*Mu*(TWO/THREE*dWdx.v.x - dWdx.v.y/THREE - v.x/(THREE*radius)) 
		     + TWO*v.x*Mu*(TWO/THREE*d_dWdx_dW-ONE/(THREE*radius)))/radius;    
    dSa_VdW(3,2) += (Mu*(dWdy.v.x +dWdx.v.y) + v.y*Mu*d_dWdx_dW - TWO/THREE*Mu*v.x*d_dWdy_dW)/radius;
    dSa_VdW(3,3) += (Kappa*(d_dWdx_dW - dWdx.rho/rho)/(rho*Rmix)+Sum_q)/radius;
   
    //Multispecies
    for(int Num = 0; Num<(ns-1); ++Num){ 
      dSa_VdW(3,NUM_VAR+Num) += (Mu/Schmidt[Num]*(specdata[Num].Enthalpy(Temp)+specdata[Num].Heatofform())*d_dWdx_dW)/radius;  
      //- specdata[Num].Rs()*Sum_dhi)/radius; 
      dSa_VdW(NUM_VAR+Num,NUM_VAR+Num) += (Mu/Schmidt[Num]*d_dWdx_dW)/radius;
    }
    
    if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY ||
	Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {
      cout<<"\n Missing Turbulent dSa_vdU components for AXISYMMETRIC_X ";
    }

  } else if(Axisymmetric == AXISYMMETRIC_Y){ 
    cout<<"\n AXISYMMETRIC_Y NOT DONE YET ";
  }
       
}


/***************************************************************
 * Turbulence model source terms                               *  
 ***************************************************************/
LESPremixed2D_cState LESPremixed2D_pState::S_turbulence_model(const LESPremixed2D_pState &dWdx,
							      const LESPremixed2D_pState &dWdy,
							      const Vector2D &X,
							      const int Flow_Type,
							      const int Axisymmetric){
  double radius;
  double production, dissipation;
  LESPremixed2D_cState Temp; Temp.Vacuum();

  //Tensor2D strain_rate;
  //strain_rate = Strain_Rate(dWdx,dWdy, Flow_Type, Axisymmetric, X);

  //Turbulence model eddy viscosity
  //mut = mu_t(strain_rate); 
 
  if (Axisymmetric == AXISYMMETRIC_X) {
    if(X.x !=0){
       radius = X.x;
     }
  } else if (Axisymmetric == AXISYMMETRIC_Y) {
     if(X.y !=0){
        radius = X.y;
     }
  } /* endif */
  

  //Production of SFS k
  production = lambda.xx*dWdx.v.x + lambda.xy*(dWdy.v.x + dWdx.v.y) + 
               lambda.yy*dWdy.v.y;


  if (Axisymmetric == AXISYMMETRIC_X) {
    production += lambda.zz*v.x/radius;
  } else if (Axisymmetric == AXISYMMETRIC_Y) {
    production += lambda.zz*v.y/radius;    
  } /* endif */
 

  //Dissipation of SFS K
  dissipation = rho*(SFSmodel.CEPS_coef)*pow(k(), 1.5)/filter_width;
 
  Temp.rhoscalar[0] = production - dissipation;
  
  return (Temp);

}


/***********************************************************************
 ***********************************************************************
 ** LESPremixed2D_pState::Sw -- Chemical Reaction Rate Source Terms.  **
 **                                                                   **
 ** Using the Reaction class to get the source terms for the          ** 
 ** specific "Reactionset".                                           ** 
 ***********************************************************************
 ***********************************************************************/
LESPremixed2D_cState LESPremixed2D_pState::Sw(int &REACT_SET_FLAG, 
					      const int Flow_Type) const {
  LESPremixed2D_cState NEW;     
  NEW.Vacuum();

  //Adds concentration rate of change for species 1->N
  if( REACT_SET_FLAG != NO_REACTIONS){
    //bool test = negative_speccheck();            //FOR TESTING 
    React.omega(NEW,*this,Flow_Type );  
  }

#ifdef THICKENED_FLAME_ON
  return (flame.WF/flame.TF)*NEW;
#else     
  return NEW;
#endif
}

/************* Chemical Source Term Jacobian ****************************/
void LESPremixed2D_pState::dSwdU(DenseMatrix &dSwdU, const int &Flow_Type,const int &solver_type) const {
  React.dSwdU(dSwdU,*this,false, Flow_Type,solver_type);

#ifdef THICKENED_FLAME_ON
  dSwdU = (flame.WF/flame.TF)*dSwdU;
#endif
}

void LESPremixed2D_pState::dSwdU_FD(DenseMatrix &dSwdU, const int Flow_Type) const{

  LESPremixed2D_cState A,C;
  LESPremixed2D_cState B,D;
  double perturb = 1e-6; //5e-6;
  double a;

  for(int jcol=0; jcol<NUM_VAR_LESPREMIXED2D-1; ++jcol){    
    A =U(*this);  C=U(*this);
    a =  perturb*max(ONE,A[jcol+1]);

    if( jcol <NUM_LESPREMIXED2D_VAR_SANS_SPECIES) {
      A[jcol+1] += a;
 //      C[jcol+1] -= perturb*max(ONE,C[jcol+1]);
    } else {                                       //enforce sum(ci) = 1;
      A[jcol+1] += a;
      A[NUM_VAR_LESPREMIXED2D] -= a;            
//       if(C[jcol+1] - a > ZERO) {  C[jcol+1] -= a; C[NUM_VAR_LESPREMIXED2D] += a;}      
    }   
    B = A.W().Sw(React.reactset_flag,Flow_Type);   D = C.W().Sw(React.reactset_flag,Flow_Type);
    for(int irow=0; irow<NUM_VAR_LESPREMIXED2D-1; ++irow){
      //dSwdU(irow,jcol) += ( B[irow+1] - D[irow+1])/(TWO*perturb*max(ONE, U(*this)[jcol+1]));      
      dSwdU(irow,jcol) += ( B[irow+1] - D[irow+1])/a;      
    }
  }

#ifdef THICKENED_FLAME_ON
  dSwdU = (flame.WF/flame.TF)*dSwdU;
#endif
}

/************* Max Diagonal of Jacobian for CFL **************************/
double LESPremixed2D_pState::dSwdU_max_diagonal(const int &Preconditioned,
						const int &flow_type_flag,
						const double &delta_n,
						const int &solver_type) const {

  //this is expensive as I am recalculating the whole Jacobian
  //as above, but its really easy to setup.
  //should change later to only calculate the diagonal terms!!!!

  double max_diagonal =ONE;
  DenseMatrix dSwdU(NUM_VAR_LESPREMIXED2D-1,NUM_VAR_LESPREMIXED2D-1,ZERO);
  React.dSwdU(dSwdU,*this,true,flow_type_flag,solver_type);

#ifdef THICKENED_FLAME_ON
  dSwdU = (flame.WF/flame.TF)*dSwdU;  //???????
#endif

//   if(Preconditioned){
//     DenseMatrix Pinv(NUM_VAR_LESPREMIXED2D-1,NUM_VAR_LESPREMIXED2D-1);
//     Low_Mach_Number_Preconditioner_Inverse(Pinv,flow_type_flag,delta_n);
//     dSwdU = Pinv*dSwdU;
//   }

  for(int i=0; i < NUM_VAR_LESPREMIXED2D-1; ++i){
    max_diagonal = max(max_diagonal,fabs(dSwdU(i,i)));
    // cout<<"\n "<<i<<" "<<dSwdU(i,i);
  }
  return max_diagonal;
}


/*********************************************************************
 *********************************************************************
 ** LESPremixed2D_pState::Sg -- Source Terms for gravitational body **
 **                             force.                              **
 **                                                                 **
 *********************************************************************
 *********************************************************************/
LESPremixed2D_cState LESPremixed2D_pState::Sg(void) const {
  LESPremixed2D_cState NEW;     
  NEW.Vacuum();
  //Gravity only in the z or axial direction
  //  z|   
  //   |___ r
  //   
  NEW[1] = ZERO;
  NEW[2] = ZERO;  //rho*g_r
  NEW[3] = rho*gravity_z;
  NEW[4] = rho*gravity_z*v.y;

  return NEW;

}

/************* Gravitational Term Jacobian ****************************/
void LESPremixed2D_pState::dSgdU(DenseMatrix &dSgdU) const {
  dSgdU(2,0) += gravity_z;
  dSgdU(3,2) += gravity_z;
}

/*****************************************************************
 *****************************************************************
 ** LESPremixed2D_pState::S_dual_time_stepping                  **
 **                                                             **
 **                -- Source Terms for dual time stepping.      **
 **                                                             **
 *****************************************************************
 *****************************************************************/
LESPremixed2D_cState LESPremixed2D_pState::S_dual_time_stepping(const LESPremixed2D_cState &U,
								const LESPremixed2D_cState &Ut,
								const LESPremixed2D_cState &Uold,
								const double &dTime,
								const int &first_step) const {
LESPremixed2D_cState NEW;
 if (first_step) {
   // Implicit Euler
   //cout << "\n First step, Implicit Euler" << endl;
   NEW = (U - Ut)/dTime;
 } else {
   // Second-order backward implicit
   //cout << "\n Second-order backward implicit" << endl;
   NEW = (THREE*U - FOUR*Ut + Uold)/(TWO*dTime);
 }
 
  return NEW; 
}


/**************************************************************************
 **************** LESPREMIXED2D_CSTATE CONSTRUCTORS ***********************
 **************************************************************************/

/**************************************************
  mixture gas constant  J/(kg*K)
***************************************************/
double LESPremixed2D_cState::Rtot() const{
  // = sum ( mass fraction * species gas constant)
  double sum = 0.0;
  for(int i=0; i<ns; ++i){
    sum += rhospec[i].c * specdata[i].Rs();
  }
  return (sum/rho);
}

// /**************************************************
//   mixture Heat Capacity (const pressure) J/(kg*K)
// ***************************************************/
// double LESPremixed2D_cState::Cp(void) const{
//   // = sum ( mass fraction * species Cp) 
//   double Temp = T();
//   double sum = 0.0;
//   for(int i=0; i<ns; ++i){
//     sum += rhospec[i].c*specdata[i].HeatCapacity_p(Temp);
//   }
//   return (sum/rho);
// }

// /**************************************************
//   mixture Heat Capacity (const volume) J/(kg*K)
// ***************************************************/
// double LESPremixed2D_cState::Cv(void) const{
//   // = sum ( mass fraction * species Cv)  
//   double Temp = T();
//   double sum = 0.0;
//   for(int i=0; i<ns; ++i){
//     sum += rhospec[i].c*specdata[i].HeatCapacity_v(Temp);
//   }
//   return (sum/rho);
// }

// /**************************************************
//   mixture Heat Ratio gamma J/(kg*K)
// ***************************************************/
// double LESPremixed2D_cState::g(void) const{
//   // = Cp / Cv  
//   return Cp()/Cv();
// }

/**************************************************
  Specific Internal Energy
 ***************************************************/
double LESPremixed2D_cState::e(void) const{
  // = sum (mass fraction * species e) 
  double sum = 0.0;
  double Temp = T();
  for(int i=0; i<ns; ++i){ //(Enthalpy(Temp) - (R/mol_mass)*Temp)
    sum += rhospec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform() 
			 - specdata[i].Rs()*Temp);
  }
  return (sum/rho);
}

double LESPremixed2D_cState::es(void) const{
  // = sum (mass fraction * species e) 
  double sum = 0.0;
  double Temp = T();
  for(int i=0; i<ns; ++i){ //(Enthalpy(Temp) - (R/mol_mass)*Temp)
    sum += rhospec[i].c*(specdata[i].Enthalpy(Temp) - specdata[i].Rs()*Temp);
  }
  return (sum/rho);
}

/**************************************************
 Specific Absolute enthalpy
***************************************************/
double LESPremixed2D_cState::h(const double &Temp) const{
  // = sum (mass fraction * species h) 
 double sum = 0.0;  
 for(int i=0; i<ns; ++i){
   sum += rhospec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform());
 }
 return (sum/rho);
}

double LESPremixed2D_cState::hs(const double &Temp) const{
  // = sum (mass fraction * species h) 
 double sum = 0.0;  
 for(int i=0; i<ns; ++i){
   sum += rhospec[i].c*(specdata[i].Enthalpy(Temp));
 }
 return (sum/rho);
}

/**************************************************
   Derivative of specific enthalpy dh/dT
   actually is just Cp as Cp = (dh/dT)_p
***************************************************/
double LESPremixed2D_cState::hprime(const double &Temp) const{
 double sum = 0.0;  
 for(int i=0; i<ns; ++i){
   sum += rhospec[i].c*specdata[i].Enthalpy_prime(Temp);
 }
 return (sum/rho);
}

/************* Mixture Heats of Formation *********/
double LESPremixed2D_cState::heatofform(void) const{ 
  double sum = 0.0;
  for(int i=0; i<ns; ++i){ 
    sum += rhospec[i].c*specdata[i].Heatofform();
  }
  return (sum);
}


/**************************************************
  polytropic heat ratio mixture gamma J/(kg*K)
  assuming T=273K as the temperature.
 **************************************************/
double LESPremixed2D_cState::gamma_guess(void) const{
  double sum1 = ZERO;     double sum2 = ZERO;
  double gamma_s = ZERO;  double Temp = 200.0;
  for(int i=0; i<ns; ++i){
    sum1 += (rhospec[i].c/rho)*specdata[i].Rs();
    gamma_s = ONE/(ONE - specdata[i].Rs()/ specdata[i].HeatCapacity_p(Temp));
    sum2 += (((rhospec[i].c/rho)*specdata[i].Rs()) / (gamma_s - ONE)); 
  }

  return (ONE + sum1/sum2);
}

/********************* Temperature ***********************************
  This gets complicated when trying to obtain the primitive variables 
  from the conserved as E is a function of Temperature taken from 
  a nonlinear equation (ie. a polynomial).   Thus is can't be rearranged
  to find T and p, so a iterative newtons method has to used. Simple 
  enough, but can be extremely expensive in terms of computation time.
  Oh well ce la vie...

  Also note the "E" is actually rho*E is. E = (rho *(e + HALF*v^2))
  in the conserved state.

  If more that 20 iterations are taken than the look jumps out 
  and gives a warning but will continue.  On average it almost always 
  converges in less than 5 iterations with "tolerance" which is defined
  in the header.
**********************************************************************/
double LESPremixed2D_cState::T(void) const{
  double T = ZERO;
  double RTOT = Rtot();  
  //--------- Initial Guess ------------------------------//
  //using a polytropic gas assumption with gamma@200;
  double Tguess = (gamma_guess() - ONE)*(E - HALF*rhov.sqr()/rho-rhok())/(rho*RTOT);
  //--------- global newtons method to get T ---------------//
  double A = (E - HALF*rhov*rhov/rho -rhok())/rho;
  //Note that there is k (turbulent kinetic energy) in above both variables 
  //Need to set a flag for the choice of laminar and turbulent flow 
  int numit =0;
  double Tmin = low_temp_range;
  double Tmax = high_temp_range;

  //check for start value
  if(Tguess > Tmin && Tguess < Tmax){
    T=Tguess;
  } else {
    T=Tmin;
  }
  
  double fa = h(Tmin) - Tmin*RTOT - A;
  double fn = h(T) - T*RTOT - A;
  double dfn = hprime(T) - RTOT;
  while( fabs(Tmax-Tmin) > CONV_TOLERANCE && fabs(fn) > CONV_TOLERANCE && numit<20 && T >= low_temp_range){    
    // Newton 
    if(T >= Tmin && T <= Tmax){
      T = T - fn/dfn;
      if(T >= Tmax) T = HALF*(Tmax + Tmin);	
      //Bisection
    } else {
      T = HALF*(Tmax + Tmin);
    } 
    //evaluate function and derivative
    fn = h(T) - T*RTOT - A;
    dfn = hprime(T) - RTOT;  
    //change bisection range
    if ( fa*fn <=ZERO){
      Tmax = T;
    } else {
      Tmin = T;
      fa = fn;
    }
    numit++;
  }  
  if (numit>=19 || T <= low_temp_range){
    T = max(Tguess,low_temp_range); 	
    cout<<"\nTemperature didn't converge in LESPremixed2D_cState::T(void)";
    cout<<" with polytopic Tguess "<<Tguess<<", or lower than Tmin "<<low_temp_range<<" using "<<T;
  }

  return T;
} 


/****************************************************
  Speed of sound using 
  a^2 = dip/dirho + p/rho^2( die/dip)^-1
  from eigenvalue analysis using e =f(p,rho)
****************************************************/
double LESPremixed2D_cState::a(void) const{
  double sum;
  double RTOT= Rtot();
  double Temp= T();
  sum = RTOT*Temp*(RTOT/( hprime(Temp) - RTOT) + ONE);
  //could also just call sqrt(g()*Rtot()*T());
  return sqrt(sum);
}

/**************************************************
  Turbulence model related parameters
***************************************************/
double LESPremixed2D_cState::mu_t(const Tensor2D &strain_rate) const{
  double mut = rho*SFSmodel.eddy_viscosity_Smagorinsky(strain_rate, filter_width);
#ifdef THICKENED_FLAME_ON
  return (flame.WF*flame.TF)*mut;
#else
  return mut;
#endif
}

double LESPremixed2D_cState::Pr_turb(void) const{
  return (0.75);
}

double LESPremixed2D_cState::Sc_turb(void) const{
  return (0.9);
}

double LESPremixed2D_cState::Dm_turb(const double &mut) const{
  return (mut/(rho*Sc_turb()));
}


/**************************************************
  Viscosity 
  using Wilke [1950] formulation
***************************************************/
double LESPremixed2D_cState::mu(void) const{
  double sum =0.0; 
  double Temp = T();

  for(int i=0; i<ns; ++i){
    double phi = 0.0;
    for (int j=0; j<ns; ++j){
      phi += ((rhospec[j].c/rho) / specdata[j].Mol_mass())*
	pow(1.0 + sqrt(specdata[i].Viscosity(Temp)/specdata[j].Viscosity(Temp))*
	    pow(specdata[j].Mol_mass()/specdata[i].Mol_mass(),0.25),2.0)/
       sqrt(8.0*(1.0 +specdata[i].Mol_mass()/specdata[j].Mol_mass()));
    }
    sum += ((rhospec[i].c/rho)* specdata[i].Viscosity(Temp) ) / 
      (specdata[i].Mol_mass() * phi);
  }  

#ifdef THICKENED_FLAME_ON
  return (flame.WF*flame.TF)*sum;
#else
  return sum;
#endif
}

double LESPremixed2D_cState::dmudT(void) const{
  double sum =0.0; 
  double Temp = T();

  for(int i=0; i<ns; ++i){
    double phi = 0.0;
    for (int j=0; j<ns; ++j){
      phi += (rhospec[j].c/rho / specdata[j].Mol_mass())*
	pow(1.0 + sqrt(specdata[i].dViscositydT(Temp)/specdata[j].dViscositydT(Temp))*
	    pow(specdata[j].Mol_mass()/specdata[i].Mol_mass(),0.25),2.0)/
       sqrt(8.0*(1.0 +specdata[i].Mol_mass()/specdata[j].Mol_mass()));
    }
    sum += (rhospec[i].c/rho * specdata[i].dViscositydT(Temp) ) / 
      (specdata[i].Mol_mass() * phi);
  }
#ifdef THICKENED_FLAME_ON
  return (flame.WF*flame.TF)*sum;
#else  
  return sum;
#endif
}

/******************************************************
 Calculating the thermal diffusion component of 
 the heat flux vector (qflux)

  sum( hs * Ds * grad cs)
*******************************************************/
Vector2D LESPremixed2D_cState::thermal_diffusion(const double &Temp) const{
  Vector2D sum;
  sum.zero();
  //double Temp = T();
  //problems with Species overloaded operators
  for(int i=0; i<ns; ++i){
#ifdef THICKENED_FLAME_ON
    sum  += (specdata[i].Enthalpy(Temp) + specdata[i].Heatofform())
      * (flame.WF*flame.TF)*rhospec[i].diffusion_coef * rhospec[i].gradc;
#else 
    sum  += (specdata[i].Enthalpy(Temp) + specdata[i].Heatofform())
      * rhospec[i].diffusion_coef*rhospec[i].gradc;
#endif
  }
  return sum/(rho*rho);
}

/*****************************************************************
 * Viscous fluxes  (laminar flow)                                * 
 * Viscous fluxes  (turbulent flows) are defined in single block * 
 ****************************************************************/
LESPremixed2D_cState LESPremixed2D_cState::Viscous_Flux_x(const LESPremixed2D_pState &dWdx,
							  const LESPremixed2D_pState &dWdy,
							  const int Flow_Type,
							  const int Axisymmetric,
							  const Vector2D X) const{
 
  LESPremixed2D_cState temp;

  temp[1] = ZERO;
  temp[2] = tau.xx;
  temp[3] = tau.xy;
  temp[4] = - qflux.x + v().x*tau.xx + v().y*tau.xy;		
  //rho * Diffusion_Coef * grad cn 
  for( int i=0; i<ns; ++i){
#ifdef THICKENED_FLAME_ON
    temp.rhospec[i].c = ( (flame.WF*flame.TF)*rhospec[i].diffusion_coef * rhospec[i].gradc.x)/rho;
#else
    temp.rhospec[i].c = (rhospec[i].diffusion_coef * rhospec[i].gradc.x)/rho;
#endif 
  }

  if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY ||
      Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {
    Tensor2D strain_rate = Strain_Rate(dWdx,dWdy, Flow_Type, Axisymmetric, X);
    double mut = mu_t(strain_rate); 
    double Dm_t = Dm_turb(mut);
    
    temp[2] += lambda.xx + 2.0*rhok()/3.0; 
    temp[3] += lambda.xy;
    temp[4] += - theta.x + v().x*(lambda.xx + 2.0*rhok()/3.0) + v().y*lambda.xy;
    // For the SFS k-equation
    if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {
      temp.rhoscalar[0] = (mu() + mut)*dWdx.scalar[0];
    }
    //species    
    for(int i=0; i<ns; ++i){
      temp.rhospec[i].c += Dm_t*rhospec[i].gradc.x; 
    }
  }
 
  return(temp);  
}

LESPremixed2D_cState LESPremixed2D_cState::Viscous_Flux_y(const LESPremixed2D_pState &dWdx,
							  const LESPremixed2D_pState &dWdy,
							  const int Flow_Type,
							  const int Axisymmetric,
							  const Vector2D X) const {
  LESPremixed2D_cState temp;

  temp[1] = ZERO;
  temp[2] = tau.xy; 
  temp[3] = tau.yy;
  temp[4] = - qflux.y + v().x*tau.xy + v().y*tau.yy;		
  //rho * Diffusion_Coef * grad cn 
  for( int i=0; i<ns; ++i){
#ifdef THICKENED_FLAME_ON
    temp.rhospec[i].c = ( (flame.WF*flame.TF)*rhospec[i].diffusion_coef * rhospec[i].gradc.y)/rho;
#else
    temp.rhospec[i].c = (rhospec[i].diffusion_coef * rhospec[i].gradc.y)/rho;
#endif     
  }

  if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY ||
      Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {
    Tensor2D strain_rate = Strain_Rate(dWdx,dWdy, Flow_Type, Axisymmetric, X);
    double mut = mu_t(strain_rate); 
    double Dm_t = Dm_turb(mut);

    temp[2] += lambda.xy; 
    temp[3] += lambda.yy + 2.0*rhok()/3.0;
    temp[4] += - theta.y + v().x*lambda.xy + v().y*(lambda.yy + 2.0*rhok()/3.0);

    // For the SFS k-equation
    if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {
      temp.rhoscalar[0] = (mu() + mut)*dWdy.scalar[0];
    }

    // species
    for(int i=0; i<ns; ++i){
      temp.rhospec[i].c += Dm_t*rhospec[i].gradc.y; 
    }
  }

  return(temp);  
}

/*********************************************************************
 *************** Overloaded Operators ********************************
 *********************************************************************/

//----------------- Addition -----------------------------//
LESPremixed2D_cState LESPremixed2D_cState::operator +(const LESPremixed2D_cState &U) const{ 
  LESPremixed2D_cState Temp(*this);
  Temp += U;
  return Temp;
}

//------------------ Subtraction ------------------------//
LESPremixed2D_cState LESPremixed2D_cState::operator -(const LESPremixed2D_cState &U) const{
  LESPremixed2D_cState Temp(*this);
  Temp -= U;
  return Temp;
}

//---------------- Scalar Multiplication ------------------//
LESPremixed2D_cState LESPremixed2D_cState::operator *(const double &a) const{
  LESPremixed2D_cState Temp(*this);
  Temp.rho = rho*a;  Temp.rhov = rhov*a; Temp.E = E*a;
  if(nscal) for(int i=0; i<nscal; ++i) Temp.rhoscalar[i] = rhoscalar[i]*a;
  for( int i=0; i<ns; ++i) Temp.rhospec[i] = rhospec[i]*a;
  Temp.tau = tau*a;
  Temp.qflux = qflux*a;
  Temp.lambda = lambda*a;
  Temp.theta = theta*a;
#ifdef THICKENED_FLAME_ON
  Temp.flame = flame*a;
#endif
  return(Temp);
}

LESPremixed2D_cState operator *(const double &a, const LESPremixed2D_cState &U){
  LESPremixed2D_cState Temp;
  Temp.rho = U.rho*a;  Temp.rhov = U.rhov*a; Temp.E = U.E*a;
  if(U.nscal) for(int i=0; i<U.nscal; ++i) Temp.rhoscalar[i] = U.rhoscalar[i]*a;
  for(int i=0; i<U.ns; ++i) Temp.rhospec[i] = U.rhospec[i]*a;
  Temp.tau = U.tau*a;
  Temp.qflux = U.qflux*a;
  Temp.lambda = U.lambda*a;
  Temp.theta = U.theta*a;
#ifdef THICKENED_FLAME_ON
  Temp.flame = U.flame*a;
#endif
  return(Temp);
}

//--------------- Scalar Division ------------------------//
LESPremixed2D_cState LESPremixed2D_cState::operator /(const double &a) const {
  LESPremixed2D_cState Temp(*this);
  Temp.rho = rho/a; Temp.rhov = rhov/a; Temp.E = E/a;
  if(nscal) for(int i=0; i<nscal; ++i) Temp.rhoscalar[i] = rhoscalar[i]/a;
  for(int i=0; i<ns; ++i) Temp.rhospec[i] = rhospec[i]/a; 
  Temp.tau = tau/a;
  Temp.qflux = qflux/a;
  Temp.lambda = lambda/a;
  Temp.theta = theta/a;
#ifdef THICKENED_FLAME_ON
  Temp.flame = flame/a;
#endif
  return(Temp);
}

//----------------- Inner Product ------------------------//
double LESPremixed2D_cState::operator *(const LESPremixed2D_cState &U) const{
  double sum=0.0; 
  double sumscal=0.0;
  if(nscal) for(int i=0; i<nscal; ++i) sumscal += rhoscalar[i]*U.rhoscalar[i];
  for(int i=0; i<ns; ++i)  sum += rhospec[i]*U.rhospec[i];
  return (rho*U.rho + rhov*U.rhov + E*U.E + sumscal + sum);
}

//----------- solution state product operator ------------//
LESPremixed2D_cState LESPremixed2D_cState::operator ^( const LESPremixed2D_cState &U) const {
  LESPremixed2D_cState Temp(*this);
  Temp.rho = rho*U.rho;
  Temp.rhov.x = rhov.x*U.rhov.x;
  Temp.rhov.y = rhov.y*U.rhov.y;
  Temp.E = E*U.E;
  if(nscal) for(int i=0;  i<nscal; ++i) Temp.rhoscalar[i] = rhoscalar[i]*U.rhoscalar[i];
  for(int i=0; i<ns; ++i) Temp.rhospec[i] = rhospec[i]*U.rhospec[i];
#ifdef THICKENED_FLAME_ON
  Temp.flame.WF = flame.WF*U.flame.WF;
  Temp.flame.TF = flame.TF*U.flame.TF;
#endif
  return(Temp);
}

//----------------- Assignment ----------------------------//
LESPremixed2D_cState& LESPremixed2D_cState::operator =(const LESPremixed2D_cState &U){
  //self assignment protection
  if( this != &U){   
    //copy assignment
    Copy(U);
  }
  return (*this);
}


/**********************************************************
 * LESPremixed2D_cState -- Shortcut arithmetic operators. *
 **********************************************************/
LESPremixed2D_cState& LESPremixed2D_cState::operator +=(const LESPremixed2D_cState &U){
  rho += U.rho;
  rhov += U.rhov; 
  E += U.E;
  if(nscal) for(int i=0; i<nscal; ++i) rhoscalar[i] += U.rhoscalar[i];
  for( int i=0; i<ns; ++i)  rhospec[i].c += U.rhospec[i].c;
  tau += U.tau;
  qflux += U.qflux; 
  lambda += U.lambda;
  theta += U.theta;
#ifdef THICKENED_FLAME_ON
  flame += U.flame;
#endif
  return (*this);
}

LESPremixed2D_cState& LESPremixed2D_cState::operator -=(const LESPremixed2D_cState &U) {
  rho -= U.rho;
  rhov -= U.rhov;
  E -= U.E;
  if(nscal) for(int i=0; i<nscal; ++i) rhoscalar[i] -= U.rhoscalar[i];
  for(int i=0; i<ns; ++i) rhospec[i].c -= U.rhospec[i].c;
  tau -= U.tau;
  qflux -= U.qflux; 
  lambda -= U.lambda;
  theta -= U.theta; 
#ifdef THICKENED_FLAME_ON
  flame -= U.flame;
#endif
  return (*this); 
}

LESPremixed2D_cState& LESPremixed2D_cState::operator *=(const double &a) {
  rho *= a;
  rhov.x *= a;
  rhov.y *= a;
  E *= a;
  if(nscal) for(int i=0; i<nscal; ++i) rhoscalar[i] *= a;
  for (int i = 0; i < ns; ++i) rhospec[i] *= a;
  tau *= a;
  qflux *= a;
  lambda *= a;
  theta *= a;
#ifdef THICKENED_FLAME_ON
  flame.WF *= a;
  flame.TF *= a;
#endif
  return (*this);
}   

LESPremixed2D_cState& LESPremixed2D_cState::operator /=(const double &a) {
  rho /= a;
  rhov.x /= a;
  rhov.y /= a;
  E /= a;
  if(nscal) for(int i=0; i<nscal; ++i) rhoscalar[i] /= a;
  for (int i = 0; i < ns; ++i) rhospec[i] /= a;
  tau /= a;
  qflux /= a;
  lambda /= a;
  theta /= a;
#ifdef THICKENED_FLAME_ON
  flame.WF /= a;
  flame.TF /= a;
#endif
  return (*this);
}

/********************************************************
 * LESPremixed2D_cState -- Unary arithmetic operators.  *
 ********************************************************/
LESPremixed2D_cState operator -(const LESPremixed2D_cState &U) {
#ifdef STATIC_LESPREMIXED2D_SPECIES
  Species spt[STATIC_LESPREMIXED2D_SPECIES];
#else
  Species *spt= new Species[U.ns];
#endif
  
  double *scalpt = new double[U.nscal];
  for(int i=0; i<U.nscal; ++i) scalpt[i] = -U.rhoscalar[i];

  for(int i=0; i<U.ns; ++i) spt[i] = -U.rhospec[i];
  LESPremixed2D_cState Temp(-U.rho,-U.rhov, -U.E, scalpt, spt);
  Temp.tau = -U.tau;
  Temp.qflux = -U.qflux;
  Temp.lambda = -U.lambda;
  Temp.theta = -U.theta;
#ifdef THICKENED_FLAME_ON
  Temp.flame = -U.flame;
#endif

#ifndef STATIC_LESPREMIXED2D_SPECIES
  delete[] spt;
#endif
  delete[] scalpt;

  return(Temp);
}

/********************************************************
 * LESPremixed2D_cState -- Relational operators.        *
 ********************************************************/
int operator ==(const LESPremixed2D_cState &U1, const LESPremixed2D_cState &U2) {
  if(U1.ns == U2.ns){ //check that species are equal
    bool Temp, Temp2;
    //species
    for(int i=0; i<U1.ns; ++i){
      if( U1.rhospec[i] == U2.rhospec[i] ){
	Temp = true;
      } else {
	Temp = false;
	break;
      }
    }
    //scalars
    if(U1.nscal){ 
      for(int i=0; i<U1.nscal; ++i){      
	if(U1.rhoscalar[i] == U2.rhoscalar[i] ){
	  Temp2 = true;
	} else {
	  Temp2 = false;
	  break;
	}
      }
    }
  
    return (U1.rho == U2.rho && U1.rhov == U2.rhov && U1.E == U2.E 
	    && Temp2 == true && Temp == true &&
	    U1.tau == U2.tau && U1.qflux == U2.qflux&&
	    U1.lambda == U2.lambda && U1.theta == U2.theta);
  
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
}

int operator !=(const LESPremixed2D_cState &U1, const LESPremixed2D_cState &U2) {
   if(U1.ns == U2.ns){ //check that species are equal
     bool Temp = true, Temp2 = true;
    // species
    for(int i=0; i<U1.ns; ++i){
      if( U1.rhospec[i] != U2.rhospec[i] ){
	Temp = false;
	break;
      }
    }
    // scalars
    if(U1.nscal){ 
      for(int i=0; i<U1.nscal; ++i){      
	if(U1.rhoscalar[i] != U2.rhoscalar[i] ){
	  Temp2 = false;
	  break;
	}
      }
    }

    return (U1.rho != U2.rho || U1.rhov != U2.rhov || U1.E != U2.E 
	    || Temp2 != true || Temp != true ||
	    U1.tau != U2.tau || U1.qflux != U2.qflux
	    || U1.lambda != U2.lambda || U1.theta != U2.theta);
    
   } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
}

/********************************************************
 * LESPremixed2D_cState -- Input-output operators.      *
 ********************************************************/
ostream &operator << (ostream &out_file, const LESPremixed2D_cState &U) {
  //out_file.precision(20);
  out_file.setf(ios::scientific);
  out_file << " " << U.rho  << " " << U.rhov.x << " " << U.rhov.y << " " << U.E;

  if(U.nscal) for(int i=0; i<U.nscal; ++i) out_file <<" "<< U.rhoscalar[i];

  for( int i=0; i<U.ns; ++i) out_file<<" "<<U.rhospec[i];

#ifdef THICKENED_FLAME_ON
  out_file << " " << U.flame;
#endif

  out_file.unsetf(ios::scientific);
  return (out_file);
}

istream &operator >> (istream &in_file, LESPremixed2D_cState &U) {
  in_file.setf(ios::skipws);
  in_file >> U.rho >> U.rhov.x >> U.rhov.y >> U.E;

  if(U.nscal) for(int i=0; i<U.nscal; ++i) in_file >> U.rhoscalar[i];

  for( int i=0; i<U.ns; ++i) in_file >> U.rhospec[i]; 

#ifdef THICKENED_FLAME_ON
  in_file >> U.flame;
#endif

  in_file.unsetf(ios::skipws);
  return (in_file);
}

/*******************************************************************
 *******************************************************************
 *  Low Mach Number Preconditioner Based on Weiss & Smith (1995)   *
 *                                                                 *
 *******************************************************************
 *******************************************************************/
void LESPremixed2D_cState::Low_Mach_Number_Preconditioner(DenseMatrix  &P,
							  const int    &flow_type_flag, 
							  const double &deltax) const {
  LESPremixed2D_pState NEW = W();
  NEW.Low_Mach_Number_Preconditioner(P,flow_type_flag,deltax);
}

void LESPremixed2D_cState::Low_Mach_Number_Preconditioner_Inverse(DenseMatrix  &Pinv,
								  const int    &flow_type_flag, 
								  const double &deltax) const {
 
  LESPremixed2D_pState NEW = W();
  NEW.Low_Mach_Number_Preconditioner_Inverse(Pinv,flow_type_flag,deltax);
}



/*********************************************************************************
 *********************************************************************************
 *********************************************************************************
 *********************** EXTERNAL FUNCTIONS **************************************
 *********************************************************************************
 *********************************************************************************
 *********************************************************************************/

/********************************************************
 * Routine: Reflect                                     *
 *                                                      *
 * This function returns the reflected solution state   *
 * in a given direction given the primitive solution    *
 * variables and the unit normal vector in the          *
 * direction of interest.                               *
 *                                                      *
 ********************************************************/
LESPremixed2D_pState Reflect(const LESPremixed2D_pState &W,
			     const Vector2D &norm_dir) {

    double ur, vr, u, v;
    double cos_angle, sin_angle;
    LESPremixed2D_pState Temp(W);
 
    /* Determine the direction cosine's for the frame
       rotation. */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    /* Apply the frame rotation and calculate the primitive
       solution state variables in the local rotated frame
       defined by the unit normal vector. */
    ur = W.v.x*cos_angle +
         W.v.y*sin_angle;
    vr = - W.v.x*sin_angle +
           W.v.y*cos_angle;

    /* Reflect the normal velocity in the rotated frame. */
    ur = -ur;
       //cout<<"\n Reflect";
    /* Rotate back to the original Cartesian reference frame. */

    Temp.v.x = ur*cos_angle - vr*sin_angle;
    Temp.v.y = ur*sin_angle + vr*cos_angle;
  
    /* Return the reflected state. */
    return (Temp);
}
/*********** FREE SLIP  ***********************************/
LESPremixed2D_pState Free_Slip(const LESPremixed2D_pState &Win,
			       const LESPremixed2D_pState &Wout,
			       const Vector2D &norm_dir,
			       const int &TEMPERATURE_BC_FLAG) { 
  
  /***** MOHAMMED et al *******/
  //LESPremixed2D_pState Temp(Wout);
  //Temp.v = Win.v;
  
  /**** DAY & BELL ********/
  LESPremixed2D_pState Temp;
  Temp = Reflect(Win,norm_dir);
 
  //Fixed Temperature
  if(TEMPERATURE_BC_FLAG == FIXED_TEMPERATURE_WALL){
   Temp.rho = Temp.p/(Temp.Rtot()*Wout.T());
  }
  
  return Temp;
}

/*********** NO SLIP - ************************************/
LESPremixed2D_pState No_Slip(const LESPremixed2D_pState &Win,
			     const LESPremixed2D_pState &Wout,
			     const Vector2D &norm_dir,
			     const int &TEMPERATURE_BC_FLAG) {  
  
  return(Moving_Wall(Win,Wout,norm_dir,ZERO,TEMPERATURE_BC_FLAG));

}

/************ Moving_Wall **********************************/
LESPremixed2D_pState Moving_Wall(const LESPremixed2D_pState &Win,
				 const LESPremixed2D_pState &Wout,
				 const Vector2D &norm_dir, 
				 const double &wall_velocity,
				 const int &TEMPERATURE_BC_FLAG) {

  double ur, vr;
  double cos_angle, sin_angle;
  LESPremixed2D_pState Temp(Win);
  
  /* Determine the direction cosine's for the frame
     rotation. */
  
  cos_angle = norm_dir.x; 
  sin_angle = norm_dir.y;

  /* Apply the frame rotation and calculate the primitive
     solution state variables in the local rotated frame
     defined by the unit normal vector. */
  ur = Win.v.x*cos_angle +  Win.v.y*sin_angle;
  vr = - Win.v.x*sin_angle +  Win.v.y*cos_angle;
  
  /* Use the Roeaveraged value to find the correct tangential velocity */
  ur = -ur;
  vr = -TWO*wall_velocity - vr;

  /* Rotate back to the original Cartesin reference frame. */  
  Temp.v.x = ur*cos_angle - vr*sin_angle;
  Temp.v.y = ur*sin_angle + vr*cos_angle;
   
 //   /* Fixed Wall Temperature */
//    if(TEMPERATURE_BC_FLAG == FIXED_TEMPERATURE_WALL){
//      Temp.rho = Temp.p/(Temp.Rtot()*Wout.T());
//    }
   
   //ADIABATIC_WALL -> constant extrapolation for Adiabatic */


   return (Temp);

}


/********************************************************
 * Routine: BC_Flame_Inflow                             *
 *        - 1DFlame Inflow conditions                   *
 *                                                      *
 ********************************************************/
LESPremixed2D_pState BC_1DFlame_Inflow(const LESPremixed2D_pState &Wi,
				       const LESPremixed2D_pState &Wo,
				       const LESPremixed2D_pState &Woutlet,
				       const Vector2D &norm_dir){ 
  //fixed Wo
  LESPremixed2D_pState Wnew(Wo);
  Wnew.v.y = ZERO;

  //Constant extrapolate the wrinkling and thickening factor
#ifdef THICKENED_FLAME_ON
  Wnew.flame = Wi.flame;
#endif

  //Constant Extrapolate Wi 
  //Calculate upstream velocity to balance flame ie. mass flow rate 
  // rho_1*u_1 = rho_2*u_2  
  //relaxation
  //Wnew.v.x = Wi.v.x + 0.5*(Woutlet.rho*Woutlet.v.x/Wo.rho - Wi.v.x);

  //no relaxation
  Wnew.v.x = Woutlet.rho*Woutlet.v.x/Wo.rho;

  if(Wnew.v.x  < 0.1 || Wnew.v.x > (Wo.v.x + 0.5)){
    Wnew.v.x = Wo.v.x;
  }
   
  return Wnew;

}
/********************************************************
 * Routine: BC_Flame_Outflow                            *
 *        - 1DFlame outflow conditions                  *
 *                                                      *
 ********************************************************/
LESPremixed2D_pState BC_1DFlame_Outflow(const LESPremixed2D_pState &Wi,       //ICu
					const LESPremixed2D_pState &Wo,       //Ghost
					const LESPremixed2D_pState &Winlet,   //ICl
					const Vector2D &norm_dir){

  //Constant extrapolate ( zero gradient)
  LESPremixed2D_pState Wnew(Wi);
  Wnew.v.y = ZERO;

#ifdef THICKENED_FLAME_ON
  Wnew.flame = Wi.flame;
#endif

  //Calculate Pressure assuming constant mass flow rate
  //and Wo.p == Winput.p (constant pressure initial condition)
  double sum = Wi.rho*Wi.v.x*(Wi.v.x - Winlet.v.x);
  if( sum < ZERO){
    Wnew.p = Wo.p;
  } else {
    //no relaxation
    Wnew.p = Wo.p - sum;

    // relaxation
    //Wnew.p = Wi.p + 0.5*( (Wo.p-sum) - Wi.p);
  }

  return Wnew;
 
}

/********************************************************
 * Routine: BC_2DFlame_Inflow                           *
 *        - 2DFlame Inflow conditions                   *
 *                                                      *
 ********************************************************/
LESPremixed2D_pState BC_2DFlame_Inflow(const LESPremixed2D_pState &Wi,
				       const LESPremixed2D_pState &Wo,
				       const Vector2D &norm_dir){ 

  //fixed rho, v, p, and species
  LESPremixed2D_pState Wnew(Wo);
  Wnew.v.x = Wi.v.x;

//   LESPremixed2D_pState Wnew(Wi);  
//   Wnew.p = Wo.p;         //fix pressure & V velocity
//   Wnew.v.y = Wo.v.y;

  return Wnew;
 
}

/********************************************************
 * Routine: BC_2DFlame_Outflow                          *
 *        - 2DFlame outflow conditions                  *
 *                                                      *
 ********************************************************/
LESPremixed2D_pState BC_2DFlame_Outflow(const LESPremixed2D_pState &Wi, 
					const LESPremixed2D_pState &Wo,
					const Vector2D &norm_dir){
  //2D Coreflame OUTFLOW hold pressure
  LESPremixed2D_pState Wnew(Wi);  
  Wnew.p = Wo.p;  
  if(Wnew.v.y < ZERO){ 
    Wnew.v.y = ZERO;
  }
  return Wnew;

//   LESPremixed2D_pState Wi_rotated(Wi), Wo_rotated(Wo), Wnew;
//   double ab, ub_rotated, vb_rotated;
//   double cos_angle = norm_dir.x; 
//   double sin_angle = norm_dir.y;

//   Wi_rotated.v.x = Wi.v.x*cos_angle + Wi.v.y*sin_angle;
//   Wi_rotated.v.y = - Wi.v.x*sin_angle +  Wi.v.y*cos_angle; 
//   Wo_rotated.v.x = Wo.v.x*cos_angle + Wo.v.y*sin_angle;
//   Wo_rotated.v.y = - Wo.v.x*sin_angle + Wo.v.y*cos_angle;
    
//   Wnew.Copy(Wi_rotated);
//   Wnew.p = Wo_rotated.p;
//   Wnew.rho = Wi_rotated.rho*pow(Wnew.p/Wi_rotated.p, ONE/Wi_rotated.g());

//   ab = Wnew.a();
//   ub_rotated = Wi_rotated.v.x + TWO*(Wi_rotated.a()-ab)/(Wi_rotated.g() -ONE);
//   vb_rotated = Wi_rotated.v.y;

//   Wnew.v.x = ub_rotated*cos_angle - vb_rotated*sin_angle;
//   Wnew.v.y = ub_rotated*sin_angle + vb_rotated*cos_angle;

//   return (Wnew);

}

/********************************************************
 * Routine: BC_Characteristic_Pressure                  *
 *   (Characteristic-Based Boundary Condition with      *
 *    Static Pressure Specified Whenever Possible)      *
 *                                                      *
 * This function returns the boundary solution state    *
 * for a given direction given the primitive solution   * 
 * state on the interior of the boundary, Wi, the       *
 * desired flow state to be imposed at the boundary,    *
 * Wo, and the unit normal vector in the direction of   *
 * interest. A simplified characteristic analysis is    *
 * used to specify the boundary flow state in which the *
 * static pressure is specified whenever possible.  The *
 * imposition of the boundary-data respects the         *
 * directions of propogation for the solution           *
 * characteristics at the boundary and thereby ensures  *
 * that the boundary conditions are not ill-posed (i.e.,*
 * the boundary data is not under- or over-determined). *
 * The following procedure is adopted:                  *
 *                                                      *
 * 1) for supersonic outflow: constant extrapolation    *
 *    is employed to specify boundary conditions        *
 *    using known interior solution state,              *
 * 2) for subsonic outflow: boundary conditions         *
 *    specified by employing a 1D unsteady isentropic   *
 *    wave approximation to match boundary pressure,    *
 * 3) for subsonic inflow: boundary conditions          *
 *    specified by employing a 1D unsteady isentropic   *
 *    wave approximation to match both boundary state   *
 *    pressure and sound speed,                         *
 * 4) for supersonic inflow: the known boundary state   *
 *    is used to specify the boundary state.            *    
 *                                                      *
 ********************************************************/
LESPremixed2D_pState BC_Characteristic_Pressure(const LESPremixed2D_pState &Wi,
						const LESPremixed2D_pState &Wo,
						const Vector2D &norm_dir) {

  //USED FOR BUMP FLOW EXIT
  LESPremixed2D_pState Wi_rotated(Wi), Wo_rotated(Wo), Wnew;
  double mi, ab, ub_rotated, vb_rotated;
  double cos_angle, sin_angle;
  
  /* Determine the direction cosine's for the frame
     rotation. */
  
  cos_angle = norm_dir.x; 
  sin_angle = norm_dir.y;
  
  /* Apply the frame rotation and evaluate interior and 
     imposed boundary solution states in the local rotated 
       frame defined by the unit normal vector. */

  Wi_rotated.v.x = Wi.v.x*cos_angle  +  Wi.v.y*sin_angle;
  Wi_rotated.v.y = - Wi.v.x*sin_angle + Wi.v.y*cos_angle;
 
  Wo_rotated.v.x = Wo.v.x*cos_angle +  Wo.v.y*sin_angle;
  Wo_rotated.v.y = - Wo.v.x*sin_angle +  Wo.v.y*cos_angle;
 
  /* Determine the Mach number at the interior node. */
  mi = Wi_rotated.v.x/Wi_rotated.a();
  
  /* Boundary condition for supersonic outflow. */
  if (mi >= ONE) {
    Wnew.Copy(Wi);
    
    /* Boundary condition for subsonic outflow. 
       Pressure specified. */
  } else if (mi >= ZERO) {
    Wnew.Copy(Wi_rotated);
    Wnew.p = Wo_rotated.p;
    Wnew.rho = Wi_rotated.rho*pow(Wnew.p/Wi_rotated.p, ONE/Wi_rotated.g());
    
    ab = Wnew.a();
    ub_rotated = Wi_rotated.v.x + TWO*(Wi_rotated.a()-ab)/(Wi_rotated.g() -ONE);
    vb_rotated = Wi_rotated.v.y;
    
    /* Boundary condition for subsonic inflow. 
       Pressure specified. */
  } else if (mi >= -ONE) {
    Wnew.Copy(Wo_rotated);
    
    ab = Wnew.a();
    ub_rotated = Wi_rotated.v.x + TWO*(Wi_rotated.a()-ab)*(Wo_rotated.g() -ONE);
    vb_rotated = Wo_rotated.v.y;
     
    /* Boundary condition for supersonic inflow.  */
  } else {
    Wnew.Copy(Wo);
  } 

  if( fabs(mi) <=ONE){
    //rotate frame back to cartesian
    Wnew.v.x = ub_rotated*cos_angle - vb_rotated*sin_angle;
    Wnew.v.y = ub_rotated*sin_angle + vb_rotated*cos_angle;
  }

  /* Return boundary solution state. */
  return (Wnew);
  
}

/************************************************************************
 ************************************************************************
 ************** Exact Test Case Solution Functions **********************     
 ************************************************************************
 ************************************************************************/
/**********************************************************************
 * Routine: RinglebFlow                                               *
 *                                                                    *
 * This function returns the exact solution to Ringleb's flow for the *
 * location X.                                                        *
 *                                                                    *
 **********************************************************************/
LESPremixed2D_pState RinglebFlow(const LESPremixed2D_pState &Wdum,
				 const Vector2D X) {

//   LESPremixed2D_pState W;
//   double sin_theta, cos_theta, theta;
//   double f_a, f_ab;
//   double J, J_a, J_ab;
//   double rho_a, rho_ab;
//   double q, q_a, q_ab, k;
//   double c, c_a = 0.70, c_b = 0.99, c_ab;
//   double g = GAMMA_AIR, po = PRESSURE_STDATM, rhoo = DENSITY_STDATM;

//   // Use bisection method to solve for the sound speed, c.
//   while (fabs(c_a - c_b) > NANO) {
//     // Determine f_a.
//     rho_a = pow(c_a,TWO/(g-ONE));
//     J_a = ONE/c_a + ONE/(THREE*c_a*c_a*c_a) + ONE/(FIVE*c_a*c_a*c_a*c_a*c_a) - HALF*log((ONE+c_a)/(ONE-c_a));
//     q_a = sqrt((TWO/(g-ONE))*(ONE-c_a*c_a));
//     f_a = (X.x + HALF*J_a)*(X.x + HALF*J_a) + X.y*X.y - ONE/(FOUR*rho_a*rho_a*q_a*q_a*q_a*q_a);
//     // Determine f_ab.
//     c_ab = HALF*(c_a + c_b);
//     rho_ab = pow(c_ab,TWO/(g-ONE));
//     J_ab = ONE/c_ab + ONE/(THREE*c_ab*c_ab*c_ab) + ONE/(FIVE*c_ab*c_ab*c_ab*c_ab*c_ab) - HALF*log((ONE+c_ab)/(ONE-c_ab));
//     q_ab = sqrt((TWO/(g-ONE))*(ONE-c_ab*c_ab));
//     f_ab = (X.x + HALF*J_ab)*(X.x + HALF*J_ab) + X.y*X.y - ONE/(FOUR*rho_ab*rho_ab*q_ab*q_ab*q_ab*q_ab);
//     if (f_a*f_ab <= ZERO) {
//       c_b = HALF*(c_a + c_b);
//     } else {
//       c_a = HALF*(c_a + c_b);
//     }
//   }

//   // Final sound speed, density, and total velocity (speed).
//   c = HALF*(c_a + c_b);
//   q = sqrt((TWO/(g-ONE))*(ONE-c*c));
//   W.rho = pow(c,TWO/(g-ONE));
//   J = ONE/c + ONE/(THREE*c*c*c) + ONE/(FIVE*c*c*c*c*c) - HALF*log((ONE+c)/(ONE-c));
//   k = sqrt(TWO/(TWO*W.rho*(X.x+HALF*J)+ONE/(q*q)));
//   //if (k > 5.0/3.0) cout << "k = " << k << " > 5/3 @ " << X << endl;
//   sin_theta = max(ZERO,min(ONE,q/k));
//   theta = TWO*PI-asin(sin_theta);
//   sin_theta = sin(theta);
//   cos_theta = cos(theta);

//   W.rho = rhoo*W.rho;
//   W.v.x = sqrt(g*po/rhoo)*q*cos_theta;
//   if (X.y < ZERO) W.v.x = -ONE*W.v.x;
//   W.v.y = sqrt(g*po/rhoo)*q*sin_theta;
//   W.p   = po*(W.rho/rhoo)*c*c;
  
//   W.g   = g;

//   // Return W state.
//   return W;

}

/**********************************************************************
 * Routine: ViscousChannelFlow                                        *
 *                                                                    *
 * This function will return the exact viscous channel flow solution  *
 * (Couette or Poiseuille flows) given an (x,y)-coordinate where x =  *
 * [0,0.2] and y = [0,0.02], an upper wall speed, and the imposed     *
 * pressure gradient.                                                 *
 *                                                                    *
 **********************************************************************/
LESPremixed2D_pState ViscousChannelFlow(const LESPremixed2D_pState &Wdum,
					const Vector2D X,
					const double Vwall,
					const double dp) {

  LESPremixed2D_pState W(Wdum);
  //W.Copy(Wdum);  //use same species and viscosity

  // Compute the exact viscous channel solution.
 //  W.rho = Wdum.rho;
//   W.p = Wdum.p;

  W.rho = DENSITY_STDATM; 
  W.p   = PRESSURE_STDATM + dp*(ONE - X.x/0.20);
 
  W.v.x =  (HALF/Wdum.mu())*(-dp/0.2)*
    (pow(X.y,TWO) -(0.001*0.001/4.0)) + Vwall*(X.y/0.001 + 0.5); 
  W.v.y = ZERO; 

  // Return W state.
  return W;

}

/**********************************************************************
 * Routine: FlatPlate                                                 *
 *                                                                    *
 * This function returns the exact solution for the flow over a flat  *
 * (adiabatic) plate (Blasius solution) at a given the position and   *
 * the freestream flow velocity.                                      *
 *                                                                    *
 **********************************************************************/
LESPremixed2D_pState FlatPlate(const LESPremixed2D_pState &Winf,
			       const Vector2D X,
			       double &eta,
			       double &f,
			       double &fp,
			       double &fpp) {

  LESPremixed2D_pState W(Winf);
  double fo, dn, sign, k1, k2, k3, k4;
  
  W.v.zero();

  // Initialize variables.
//   W.Vacuum(); W.rho = Winf.rho; W.p = Winf.p;
//   for(int i=0; i<W.ns; ++i){
//     W.spec[i].c = Winf.spec[i].c;
//   } 
  

  eta = ZERO;
  f = ZERO; fo = ZERO; fp = ZERO; fpp = 0.33206;
  sign = ONE;
  dn = 0.000005;

  // Return upstream conditions before flat plate.
  if (X.x < ZERO) return Winf;

  // Return upstream conditions with zero velocity at the leading edge
  // of the plate.
  if (X.x < TOLER) return W;

  // Determine the dimensionless similarity coordinate, eta:
  eta = X.y*sqrt(Winf.v.x/(X.x*Winf.mu()/Winf.rho));

  // If eta is greater than 8.4, for the sake of expediency, use linear
  // extrapolation to determine the value of f (fp = ONE and fpp = ZERO)
  // given the tabulated value at 8.4 (note, the analytic solution is 
  // linear in this region).
  if (eta > 8.4) {
    fp = ONE; fpp = ZERO; f = 6.67923 + fp*(eta - 8.4);
    W.v.x = fp*Winf.v.x;
    W.v.y = HALF*(eta*fp-f);
    return W;
  }

  // Compute the Blasius solution using a fourth-order Runge-Kutta method.
  for (double n = ZERO; n < eta; n += dn) {

    fo = f;

    // Increment f:
    k1 = dn*fp;
    k2 = dn*(fp+k1/2.0);
    k3 = dn*(fp+k2/2.0);
    k4 = dn*(fp+k3);
    f = f + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;
 
    // Increment fp:
    k1 = dn*fpp;
    k2 = dn*(fpp+k1/2.0);
    k3 = dn*(fpp+k2/2.0);
    k4 = dn*(fpp+k3);
    fp = fp + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;

    // Increment fpp:
    k1 = -dn*fo*fpp/2.0;
    k2 = -dn*(fo+dn/2.0)*(fpp+k1/2.0)/2.0;
    k3 = -dn*(fo+dn/2.0)*(fpp+k2/2.0)/2.0;
    k4 = -dn*(fo+dn)*(fpp+k3)/2.0;
    fpp = fpp + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;

  }

  // Compute the velocity vector at point X.
  W.v.x = fp*Winf.v.x;
  //W.v.y = HALF*(eta*fp-f);

  //if (eta <= 8.8) cout << eta << " " << f << " " << fp << " " << fpp << endl;

  // Return W state.
  return W;
}

/************************************************************************
 ************************************************************************
 ************** External Flux Function Functions     ********************     
 ************************************************************************
 ************************************************************************/
/********************************************************
 * Routine: RoeAverage (Roe Averages)                   *
 *                                                      *
 * This function returns the Roe-averaged primitive     *
 * solution state given left and right primitive        *
 * solution variables.  See Roe (1981).                 *
 *                                                      *
 *  NOTE: -> Differing species result in incorrect      *               
 *           Roeaveraged energy/pressure values, thus   *
 *           causing a energy flux inconsistency when   *
 *           used with the Roe Flux function.           *
 ********************************************************/
LESPremixed2D_pState RoeAverage(const LESPremixed2D_pState &Wl,
				const LESPremixed2D_pState &Wr) {

    double Hl, Hr, srhol, srhor;
    double Ha, ha;
    LESPremixed2D_pState Temp;

    /* Determine the left and right state specific enthalpies
       and square roots of the density. */
    Hl = Wl.Hs()/Wl.rho;
    Hr = Wr.Hs()/Wr.rho;
    srhol = sqrt(Wl.rho);
    srhor = sqrt(Wr.rho);

    /* Determine the appropriate Roe averages. */
    Temp.rho = srhol*srhor;
    Temp.v.x = (srhol*Wl.v.x+srhor*Wr.v.x)/(srhol+srhor);
    Temp.v.y = (srhol*Wl.v.y+srhor*Wr.v.y)/(srhol+srhor);

    // Scalars
    if(Wl.nscal) { 
      for(int i=0; i<Wl.nscal; i++) Temp.scalar[i] = (srhol*Wl.scalar[i] + srhor*Wr.scalar[i])/(srhol+srhor);
    }
    
    for(int i=0; i<Wl.ns; ++i){
      Temp.spec[i].c = (srhol*Wl.spec[i].c + srhor*Wr.spec[i].c)/(srhol+srhor);
    }
 
    Ha = (srhol*Hl+srhor*Hr)/(srhol+srhor);
    ha = Ha - HALF*(sqr(Temp.v.x)+sqr(Temp.v.y));
    ha -= 5.0*Temp.k()/3.0;

    //double TEMP = Temp.T(ha);
    Temp.p = Temp.rho*Temp.T(ha)*Temp.Rtot();
   
    /* Return the Roe-averged state. */
    return (Temp);     
  
}

/*********************************************************
 * Routine: FluxHLLE_x (Harten-Lax-van Leer flux         *
 *                      function, x-direction)           *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states by using the so-called Harten-Lax-    *
 * van Leer approximation for the fluxes.  See Harten,   *
 * Lax, van Leer (1983).                                 *
 *                                                       *
 *********************************************************/
LESPremixed2D_cState FluxHLLE_x(const LESPremixed2D_pState &Wl,
				const LESPremixed2D_pState &Wr,
				const int &Preconditioning,
				const int Flow_Type) {

    double wavespeed_l, wavespeed_r;
    LESPremixed2D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
    LESPremixed2D_cState Flux, dUrl;
    int NUM_VAR_LESPREMIXED2D = Wl.NUM_VAR_LESPREMIXED2D;

    /* Evaluate the Roe-average primitive solution state. */   
    Wa = RoeAverage(Wl, Wr);
 
    /* Evaluate the jumps in the conserved solution states. */
    dUrl = Wr.U()-Wl.U();

    /* Evaluate the left, right, and average state eigenvalues. */     
    lambdas_l = Wl.lambda_x();
    lambdas_r = Wr.lambda_x();
    lambdas_a = Wa.lambda_x();     
    
    /* Determine the intermediate state flux. */
    wavespeed_l = min(lambdas_l[1], lambdas_a[1]);

    //wavespeed_r = max(lambdas_r[NUM_VAR_LESPREMIXED2D],
    //                  lambdas_a[NUM_VAR_LESPREMIXED2D]);
    wavespeed_r = max(lambdas_r[4],lambdas_a[4]); //MARKTHIS XINFENG
    //   wavespeed_r = max(lambdas_r[NUM_VAR_LESPREMIXED2D-lambdas_r.ns],
    //                       lambdas_a[NUM_VAR_LESPREMIXED2D-lambdas_a.ns]); 
    wavespeed_l = min(wavespeed_l, ZERO);
    wavespeed_r = max(wavespeed_r, ZERO);
  

    if (wavespeed_l >= ZERO) {
      Flux = Wl.Fx();
    } else if (wavespeed_r <= ZERO) {
      Flux = Wr.Fx();
    } else { 

      Flux =   ( (wavespeed_r*Wl.Fx()-wavespeed_l*Wr.Fx())
		 +(wavespeed_l*wavespeed_r)*dUrl)/
	(wavespeed_r-wavespeed_l);
    } /* endif */

    /* Return solution flux. */
    return (Flux);

}

LESPremixed2D_cState FluxHLLE_x(const LESPremixed2D_cState &Ul,
				const LESPremixed2D_cState &Ur,
				const int &Preconditioning,
				const int Flow_Type) {
   return (FluxHLLE_x(Ul.W(), Ur.W(),Preconditioning, Flow_Type));
}


/*********************************************************
 * Routine: FluxHLLE_n (Harten-Lax-van Leer flux         *
 *                      function, n-direction)           *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for an arbitrary direction defined by a unit     *
 * normal vector in the direction of interest, given     *
 * left and right solution states.  The problem is       *
 * solved by first applying a frame rotation to rotate   *
 * the problem to a local frame aligned with the unit    *
 * normal vector and then by using the so-called         *
 * Harten-Lax-van Leer approximation to specify the      *
 * intermediate state fluxes in terms of the rotated     *
 * solution states.  See Harten, Lax, van Leer (1983).   *
 *                                                       *
 *********************************************************/
LESPremixed2D_cState FluxHLLE_n(const LESPremixed2D_pState &Wl,
				const LESPremixed2D_pState &Wr,
				const Vector2D &norm_dir,
				const int &Preconditioning,
				const int Flow_Type) {

    double cos_angle, sin_angle;
    LESPremixed2D_pState Wl_rotated(Wl), Wr_rotated(Wr);
    LESPremixed2D_cState Flux, Flux_rotated;
    Flux.Vacuum();
    /* Determine the direction cosine's for the frame
       rotation. */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    /* Apply the frame rotation and evaluate left and right
       solution states in the local rotated frame defined
       by the unit normal vector. */
 //    Wl_rotated.Copy(Wl);
//     Wr_rotated.Copy(Wr);


    Wl_rotated.v.x = Wl.v.x*cos_angle +
                     Wl.v.y*sin_angle;
    Wl_rotated.v.y = - Wl.v.x*sin_angle +
                       Wl.v.y*cos_angle;
   
    Wr_rotated.v.x = Wr.v.x*cos_angle +
                     Wr.v.y*sin_angle;
    Wr_rotated.v.y = - Wr.v.x*sin_angle +
                       Wr.v.y*cos_angle;

    /* Evaluate the intermediate state solution 
       flux in the rotated frame. */

    Flux_rotated = FluxHLLE_x(Wl_rotated, Wr_rotated, Preconditioning, Flow_Type);

    /* Rotate back to the original Cartesian reference
       frame and return the solution flux. */
    Flux.Copy(Flux_rotated);    

    Flux.rhov.x = Flux_rotated.rhov.x*cos_angle -
                Flux_rotated.rhov.y*sin_angle;
    Flux.rhov.y = Flux_rotated.rhov.x*sin_angle +
                Flux_rotated.rhov.y*cos_angle;
     
    Flux.zero_non_sol();
    
    return (Flux);
}

LESPremixed2D_cState FluxHLLE_n(const LESPremixed2D_cState &Ul,
				const LESPremixed2D_cState &Ur,
				const Vector2D &norm_dir,
				const int &Preconditioning,
				const int Flow_Type) {
  return (FluxHLLE_n(Ul.W(), Ur.W(),norm_dir,Preconditioning, Flow_Type));
}

/*********************************************************
 * Routine: FluxLinde (Timur Linde's flux function,      *
 *                     x-direction)                      *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states by using the Linde approximation for  *
 * the fluxes.  See Linde (1998).                        *
 *                                                       *
 *********************************************************/
LESPremixed2D_cState FluxLinde(const LESPremixed2D_pState &Wl,
			       const LESPremixed2D_pState &Wr, 
			       const int Flow_Type) {

    double wavespeed_l, wavespeed_r, wavespeed_m, rhoa, ca, dU, alpha;
    LESPremixed2D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
    LESPremixed2D_cState Flux, dFrl, dUrl, dFwave;
    int NUM_VAR_LESPREMIXED2D = Wl.NUM_VAR_LESPREMIXED2D;

    /* Evaluate the Roe-average primitive solution state. */   
    Wa = RoeAverage(Wl, Wr);

    /* Evaluate the left, right, and average state eigenvalues. */
    lambdas_l = Wl.lambda_x();
    lambdas_r = Wr.lambda_x();
    lambdas_a = Wa.lambda_x();
   
    /* Determine the intermediate state flux. */
    wavespeed_l = min(lambdas_l[1],
                      lambdas_a[1]);
  //   wavespeed_r = max(lambdas_r[NUM_VAR_LESPREMIXED2D-lambdas_r.ns],
//                       lambdas_a[NUM_VAR_LESPREMIXED2D-lambdas_a.ns]);

  wavespeed_r = max(lambdas_r[4],lambdas_a[4]);
    if (wavespeed_l >= ZERO) {
        Flux = Wl.Fx();
    } else if (wavespeed_r <= ZERO) {
        Flux = Wr.Fx();
    } else {
        dUrl = Wr.U()-Wl.U();
        dFrl = Wr.Fx()-Wl.Fx();
        wavespeed_m = Wa.v.x;
        rhoa = Wa.rho;
        ca = Wa.amodified();  //Wa.a();
        dU = fabs(dUrl.rho)/rhoa+
             fabs(dUrl.rhov.x)/(rhoa*ca)+
             fabs(dUrl.rhov.y)/(rhoa*ca)+
             fabs(dUrl.E)/(rhoa*ca*ca);
        if (dU <= TOLER) {
            alpha = ZERO;
        } else {
            dU = ONE/dU;
            dFwave = dFrl - wavespeed_m*dUrl;
            alpha = ONE - (fabs(dFwave.rho)/(rhoa*ca)+
                           fabs(dFwave.rhov.x)/(rhoa*ca*ca)+
                           fabs(dFwave.rhov.y)/(rhoa*ca*ca)+
                           fabs(dFwave.E)/(rhoa*ca*ca*ca))*dU;
 	    alpha = max(ZERO, alpha);

        } /* endif */
 
        Flux =   ( (wavespeed_r*Wl.Fx()-wavespeed_l*Wr.Fx())
                  +(wavespeed_l*wavespeed_r)*
		   (ONE-(ONE-max(wavespeed_m/wavespeed_r,
                                 wavespeed_m/wavespeed_l))*alpha)*dUrl)/
                   (wavespeed_r-wavespeed_l);
    } /* endif */
    
    /* Return solution flux. */
    return (Flux);

}

LESPremixed2D_cState FluxLinde(const LESPremixed2D_cState &Ul,
			       const LESPremixed2D_cState &Ur,
			       const int Flow_Type) {
   return (FluxLinde(Ul.W(), Ur.W(), Flow_Type));
}

/*********************************************************
 * Routine: FluxLinde_n (Timur Linde's flux function,    *
 *                       n-direction)                    *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for an arbitrary direction defined by a unit     *
 * normal vector in the direction of interest, given     *
 * left and right solution states.  The problem is       *
 * solved by first applying a frame rotation to rotate   *
 * the problem to a local frame aligned with the unit    *
 * normal vector and then by using the Linde             *
 * approximation to specify the intermediate state flux  *
 * in terms of the rotated solution states.  See Linde   *
 * (1998).                                               *
 *                                                       *
 *********************************************************/
LESPremixed2D_cState FluxLinde_n(const LESPremixed2D_pState &Wl,
				 const LESPremixed2D_pState &Wr,
				 const Vector2D &norm_dir,
				 const int Flow_Type) {

    double cos_angle, sin_angle;
    LESPremixed2D_pState Wl_rotated(Wl), Wr_rotated(Wr);
    LESPremixed2D_cState Flux, Flux_rotated;

    /* Determine the direction cosine's for the frame
       rotation. */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    /* Apply the frame rotation and evaluate left and right
       solution states in the local rotated frame defined
       by the unit normal vector. */
//     Wl_rotated.Copy(Wl);
//     Wr_rotated.Copy(Wr);

    Wl_rotated.v.x = Wl.v.x*cos_angle +
                     Wl.v.y*sin_angle;
    Wl_rotated.v.y = - Wl.v.x*sin_angle +
                       Wl.v.y*cos_angle;

    Wr_rotated.v.x = Wr.v.x*cos_angle +
                     Wr.v.y*sin_angle;
    Wr_rotated.v.y = - Wr.v.x*sin_angle +
                       Wr.v.y*cos_angle;

    /* Evaluate the intermediate state solution 
       flux in the rotated frame. */ 

    Flux_rotated = FluxLinde(Wl_rotated, Wr_rotated, Flow_Type);

    /* Rotate back to the original Cartesian reference
       frame and return the solution flux. */
    Flux.Copy(Flux_rotated);
 
    Flux.rhov.x = Flux_rotated.rhov.x*cos_angle -
                Flux_rotated.rhov.y*sin_angle;
    Flux.rhov.y = Flux_rotated.rhov.x*sin_angle +
                Flux_rotated.rhov.y*cos_angle;

    Flux.zero_non_sol();

    return (Flux);

}

LESPremixed2D_cState FluxLinde_n(const LESPremixed2D_cState &Ul,
				 const LESPremixed2D_cState &Ur,
				 const Vector2D &norm_dir,
				 const int Flow_Type) {
  return (FluxLinde_n(Ul.W(), Ur.W(), norm_dir, Flow_Type));
}


/*********************************************************
 * Routine: Rotate                                       *
 *                                                       *
 * This function returns the solution in the lcoal       *
 * rotated frame.                                        *
 *                                                       *
 *********************************************************/
LESPremixed2D_pState Rotate(const LESPremixed2D_pState &W,
			    const Vector2D &norm_dir) {

  LESPremixed2D_pState W_rotated(W);
  W_rotated.v.x =   W.v.x*norm_dir.x + W.v.y*norm_dir.y;
  W_rotated.v.y = - W.v.x*norm_dir.y + W.v.y*norm_dir.x;
  
  return (W_rotated);

}

/*********************************************************
 * Routine: HLLE_wavespeeds                              *
 *                                                       *
 * This function returns lambda plus and lambda minus    *
 * for rotated Riemann problem aligned with norm_dir     *
 * given unroated solution states Wl and Wr.             *
 * Note: wavespeed.x = wavespeed_l = lambda minus.       *
 *       wavespeed.y = wavespeed_r = lambda plus.        *
 *                                                       *
 *********************************************************/
Vector2D HLLE_wavespeeds(const LESPremixed2D_pState &Wl,
                         const LESPremixed2D_pState &Wr,
			 const Vector2D &norm_dir) {

    Vector2D wavespeed;
    LESPremixed2D_pState Wa_n, lambdas_l, lambdas_r, lambdas_a, Wl_n, Wr_n;  
    int NUM_VAR_LESPREMIXED2D = (Wl.NUM_VAR_LESPREMIXED2D );
    /* Use rotated values to calculate eignvalues */
    Wl_n = Rotate(Wl, norm_dir);
    Wr_n = Rotate(Wr, norm_dir);

    /* Evaluate the Roe-average primitive solution state. */
    Wa_n = RoeAverage(Wl_n, Wr_n);
    
    /* Evaluate the left, right, and average state eigenvalues. */

    lambdas_l = Wl_n.lambda_x();
    lambdas_r = Wr_n.lambda_x();
    lambdas_a = Wa_n.lambda_x();

    /* Determine the intermediate state flux. */
    wavespeed.x = min(lambdas_l[1],
                      lambdas_a[1]);   //u-a
    wavespeed.y = max(lambdas_r[4],
                      lambdas_a[4]);   //u+a
 
  //   wavespeed.y = max(lambdas_r[NUM_VAR_LESPREMIXED2D],
//                       lambdas_a[NUM_VAR_LESPREMIXED2D]);  //THIS IS u! not u+a WTF!!!
 
    wavespeed.x = min(wavespeed.x, ZERO); //lambda minus
    wavespeed.y = max(wavespeed.y, ZERO); //lambda plus 

 //    cout<< "lambda - = "<<wavespeed.x<< "lambda + = "<<wavespeed.y<<endl;
    return (wavespeed);

}


/********************************************************
 * Routine: WaveSpeedPos                                *
 *                                                      *
 * This function returns the positive parts of the      *
 * elemental wave speeds or eigenvalues.                *
 *                                                      *
 ********************************************************/
LESPremixed2D_pState WaveSpeedPos(const LESPremixed2D_pState &lambdas_a,
				  const LESPremixed2D_pState &lambdas_l,
				  const LESPremixed2D_pState &lambdas_r) {
   LESPremixed2D_pState NEW;   
   for(int i=1; i<=lambdas_a.NUM_VAR_LESPREMIXED2D; ++i){
     NEW[i] = HALF*(lambdas_a[i]+fabs(lambdas_a[i]));
   }
   return(NEW);
}

/********************************************************
 * Routine: WaveSpeedNeg                                *
 *                                                      *
 * This function returns the negative parts of the      *
 * elemental wave speeds or eigenvalues.                *
 *                                                      *
 ********************************************************/
LESPremixed2D_pState WaveSpeedNeg(const LESPremixed2D_pState &lambdas_a,
				  const LESPremixed2D_pState &lambdas_l,
				  const LESPremixed2D_pState &lambdas_r) {
  LESPremixed2D_pState NEW;   
  for(int i=1; i<=lambdas_a.NUM_VAR_LESPREMIXED2D; ++i){
     NEW[i] = HALF*(lambdas_a[i]-fabs(lambdas_a[i]));
   }
   return(NEW);
}

/********************************************************
 * Routine: WaveSpeedAbs                                *
 *                                                      *
 * This function returns the absolute values of the     *
 * elemental wave speeds or eigenvalues.                *
 *                                                      *
 ********************************************************/
LESPremixed2D_pState WaveSpeedAbs(const LESPremixed2D_pState &lambdas_a,
				  const LESPremixed2D_pState &lambdas_l,
				  const LESPremixed2D_pState &lambdas_r) {
   LESPremixed2D_pState NEW;   
   for(int i=1; i<=lambdas_a.NUM_VAR_LESPREMIXED2D; ++i){
     NEW[i] = fabs(lambdas_a[i]);
   }
   return(NEW);
}

/********************************************************
 * Routine: HartenFixPos (Harten Entropy Fix)           *
 *                                                      *
 * This function returns the positive parts of the      *
 * corrected elemental wave speeds or eigenvalues       *
 * according to the entropy fix of Harten (1983).       *
 *                                                      *
 ********************************************************/
LESPremixed2D_pState HartenFixPos(const LESPremixed2D_pState &lambdas_a,
				  const LESPremixed2D_pState &lambdas_l,
				  const LESPremixed2D_pState &lambdas_r) {
  LESPremixed2D_pState NEW;
  NEW.rho = HartenFixPos(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
  NEW.v.x = HALF*(lambdas_a[2]+fabs(lambdas_a[2]));
  NEW.v.y = HALF*(lambdas_a[3]+fabs(lambdas_a[3]));
  NEW.p = HartenFixPos(lambdas_a[4],lambdas_l[4],lambdas_r[4]);

  if(NEW.nscal){
    for(int i=5; i<=(NEW.NUM_VAR_LESPREMIXED2D-NEW.ns); ++i){
      NEW.scalar[i-5] = HALF*(lambdas_a[i]+fabs(lambdas_a[i])); //fabs(lambdas_a[i]); ??????
    }
  }
    
  for( int i=(NEW.NUM_VAR_LESPREMIXED2D-NEW.ns+1); i<=NEW.NUM_VAR_LESPREMIXED2D; ++i){
    NEW.spec[i-(NEW.NUM_VAR_LESPREMIXED2D-NEW.ns+1)].c = HALF*(lambdas_a[i]+fabs(lambdas_a[i]));
  }
  
  return (NEW);
}

/********************************************************
 * Routine: HartenFixNeg (Harten Entropy Fix)           *
 *                                                      *
 * This function returns the negative parts of the      *
 * corrected elemental wave speeds or eigenvalues       *
 * according to the entropy fix of Harten (1983).       *
 *                                                      *
 ********************************************************/
LESPremixed2D_pState HartenFixNeg(const LESPremixed2D_pState &lambdas_a,
				  const LESPremixed2D_pState &lambdas_l,
				  const LESPremixed2D_pState &lambdas_r) {
  LESPremixed2D_pState NEW;
  NEW.rho = HartenFixNeg(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
  NEW.v.x = HALF*(lambdas_a[2]-fabs(lambdas_a[2]));
  NEW.v.y = HALF*(lambdas_a[3]-fabs(lambdas_a[3]));
  NEW.p = HartenFixNeg(lambdas_a[4],lambdas_l[4],lambdas_r[4]);

  if(NEW.nscal){
    for(int i=5; i<=(NEW.NUM_VAR_LESPREMIXED2D-NEW.ns); ++i){
      NEW.scalar[i-5] = HALF*(lambdas_a[i]-fabs(lambdas_a[i]));  // fabs(lambdas_a[i]) ???????
    }
  }
  
  for( int i=(NEW.NUM_VAR_LESPREMIXED2D-NEW.ns+1); i<=NEW.NUM_VAR_LESPREMIXED2D; ++i){
    NEW.spec[i-(NEW.NUM_VAR_LESPREMIXED2D-NEW.ns+1)].c = HALF*(lambdas_a[i]-fabs(lambdas_a[i]));
  }
  return (NEW);
}
/********************************************************
 * Routine: HartenFixAbs (Harten Entropy Fix)           *
 *                                                      *
 * This function returns the absolute values of the     *
 * corrected elemental wave speeds or eigenvalues       *
 * according to the entropy fix of Harten (1983).       *
 *                                                      *
 ********************************************************/
LESPremixed2D_pState HartenFixAbs(const LESPremixed2D_pState &lambdas_a,
				  const LESPremixed2D_pState &lambdas_l,
				  const LESPremixed2D_pState &lambdas_r) {
  LESPremixed2D_pState NEW;
  NEW.rho = HartenFixAbs(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
  NEW.v.x = fabs(lambdas_a[2]);
  NEW.v.y = fabs(lambdas_a[3]);
  NEW.p = HartenFixAbs(lambdas_a[4],lambdas_l[4],lambdas_r[4]);

  if(NEW.nscal){
    for(int i=5; i<=(NEW.NUM_VAR_LESPREMIXED2D-NEW.ns); ++i){
      NEW.scalar[i-5] = fabs(lambdas_a[i]);
    }
  }
  
  for( int i=(NEW.NUM_VAR_LESPREMIXED2D-NEW.ns+1); i<=NEW.NUM_VAR_LESPREMIXED2D; ++i){
    NEW.spec[i-(NEW.NUM_VAR_LESPREMIXED2D-NEW.ns+1)].c = fabs(lambdas_a[i]);
  }
  return (NEW);
}
/*********************************************************
 * Routine: FluxRoe_x (Roe's flux function, x-direction) *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for the x-direction given left and right         *
 * solution states by using the "linearized" approximate *
 * Riemann solver of Roe for the two states.  See Roe    *
 * (1981).                                               *
 *                                                       *
 *********************************************************/
LESPremixed2D_cState FluxRoe_x(const LESPremixed2D_pState &Wl,
			       const LESPremixed2D_pState &Wr,
			       const int &Preconditioning,
			       const int &flow_type_flag,
			       const double &deltax) {



    LESPremixed2D_pState Wa, dWrl, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
    LESPremixed2D_cState Flux;

     /* Evaluate the Roe-average primitive solution state. */    
    Wa = RoeAverage(Wl, Wr);

    /* Evaluate the jumps in the primitive solution states. */
    dWrl = Wr-Wl;

    /* Evaluate the left, right, and average state eigenvalues. */
    if(!Preconditioning){
      lambdas_l = Wl.lambda_x();
      lambdas_r = Wr.lambda_x();
      lambdas_a = Wa.lambda_x();

//       //////TEST FOR DRDU///////////////
//       wavespeeds = HartenFixAbs(lambdas_a,
// 				lambdas_l,
// 				lambdas_r);
//       Flux = HALF*(Wl.Fx()+Wr.Fx()); 
//       LESPremixed2D_cState Flux_dissipation(ZERO);     
//       for ( int i = 1 ; i <= Wa.NUM_VAR_LESPREMIXED2D ; ++i ) {
// 	   Flux_dissipation -= HALF*wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
//       }
//       Flux += Flux_dissipation;
//       /////////////////////////////////

      /* Determine the intermediate state flux. */
      if (Wa.v.x >= ZERO) {
        Flux = Wl.Fx();   
        wavespeeds = HartenFixNeg(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
	
        for (int i=1 ; i < Wl.NUM_VAR_LESPREMIXED2D; ++i) {
	  if (wavespeeds[i] < ZERO) {
 	    Flux += wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
	  }
        } 
      } else {
        Flux = Wr.Fx();
        wavespeeds = HartenFixPos(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
        for (int i=1; i < Wl.NUM_VAR_LESPREMIXED2D; ++i) {
	  if (wavespeeds[i] > ZERO) {
	    Flux -= wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
          }
        } 
      } 
   
      /******* LOW MACH NUMBER PRECONDITIONING ********************/
      /* Evaluate the left, right, and average state eigenvalues. */
    } else if(Preconditioning){
	
      //calculating Mr^2 and passing to save computation,
      //not conceptually nice but saves from recalculating

      double MR2a = Wa.Mr2(flow_type_flag,deltax); 
      lambdas_l = Wl.lambda_preconditioned_x(Wl.Mr2(flow_type_flag,deltax)); 
      lambdas_r = Wr.lambda_preconditioned_x(Wr.Mr2(flow_type_flag,deltax));
      lambdas_a = Wa.lambda_preconditioned_x(MR2a);
                                                                    
      /* Evaluate the jumps in the primitive solution states. */
//       wavespeeds = WaveSpeedAbs(lambdas_a,
// 				lambdas_l,
// 				lambdas_r);
      wavespeeds = HartenFixAbs(lambdas_a,
				lambdas_l,
				lambdas_r);
          
      DenseMatrix P(Wa.NUM_VAR_LESPREMIXED2D-1,Wa.NUM_VAR_LESPREMIXED2D-1);     //COULD BE STORED IN CLASS AS STATIC AND REUSED REDUCING OVERHEAD???
      /* Evaluate the low-Mach-number local preconditioner for the Roe-averaged state. */  
    
      Wa.Low_Mach_Number_Preconditioner(P,flow_type_flag,deltax);
                                                                                     
      /* Determine the intermediate state flux. */                                                
      Flux = HALF*(Wl.Fx()+Wr.Fx()); 
      LESPremixed2D_cState Flux_dissipation(ZERO);   
    
      for ( int i = 1 ; i < Wa.NUM_VAR_LESPREMIXED2D ; ++i ) {
	Flux_dissipation -= HALF*wavespeeds[i]*(Wa.lp_x_precon(i,MR2a)*dWrl)*Wa.rc_x_precon(i,MR2a);
      }
  
      for ( int i = 1 ; i < Wa.NUM_VAR_LESPREMIXED2D ; ++i ) {
	for ( int j = 1 ; j < Wa.NUM_VAR_LESPREMIXED2D ; ++j ) {
	  Flux[i] += P(i-1,j-1)*Flux_dissipation[j]; // Add preconditioned upwind dissipation flux.
	} 
      } 

    } else {
      cerr<<"\n Not a valid Preconditioner "<<Preconditioning;
      exit(1);
    }  
 
    /* Return solution flux. */    
    return (Flux);    
}

LESPremixed2D_cState FluxRoe_x(const LESPremixed2D_cState &Ul,
			       const LESPremixed2D_cState &Ur,
			       const int &Preconditioning, 
			       const int &flow_type_flag,
			       const double &deltax) {
   return (FluxRoe_x(Ul.W(), Ur.W(),Preconditioning,flow_type_flag,deltax));
}

/*********************************************************
 * Routine: FluxRoe_n (Roe's flux function, n-direction) *
 *                                                       *
 * This function returns the intermediate state solution *
 * flux for an arbitrary direction defined by a unit     *
 * normal vector in the direction of interest, given     *
 * left and right solution states.  The problem is       *
 * solved by first applying a frame rotation to rotate   *
 * the problem to a local frame aligned with the unit    *
 * normal vector and then by using the "linearized"      *
 * approximate Riemann solver of Roe to specify the flux *
 * in terms of the rotated solution states.  See Roe     *
 * (1981).                                               *
 *                                                       *
 *********************************************************/
LESPremixed2D_cState FluxRoe_n(const LESPremixed2D_pState &Wl,
			       const LESPremixed2D_pState &Wr,
			       const Vector2D &norm_dir,
			       const int &Preconditioning,
			       const int &flow_type_flag,
			       const double &delta_n ) {

 
  double cos_angle, sin_angle, delta_rotated;
  LESPremixed2D_pState Wl_rotated(Wl), Wr_rotated(Wr);
  LESPremixed2D_cState Flux, Flux_rotated;

  /* Determine the direction cosine's for the frame rotation. */

  cos_angle = norm_dir.x; 
  sin_angle = norm_dir.y;

  /* Apply the frame rotation and evaluate left and right
     solution states in the local rotated frame defined
     by the unit normal vector. */
  
  Wl_rotated.v.x = Wl.v.x*cos_angle +  Wl.v.y*sin_angle;
  Wl_rotated.v.y = - Wl.v.x*sin_angle + Wl.v.y*cos_angle;
  
  Wr_rotated.v.x = Wr.v.x*cos_angle + Wr.v.y*sin_angle;
  Wr_rotated.v.y = - Wr.v.x*sin_angle + Wr.v.y*cos_angle;
  
  /* Evaluate the intermediate state solution 
     flux in the rotated frame. */
  
  Flux_rotated = FluxRoe_x(Wl_rotated, Wr_rotated, 
			   Preconditioning,flow_type_flag,delta_n);
  
  /* Rotate back to the original Cartesian reference
     frame and return the solution flux. */
  Flux.Copy(Flux_rotated);
  
  Flux.rhov.x = Flux_rotated.rhov.x*cos_angle - Flux_rotated.rhov.y*sin_angle;
  Flux.rhov.y = Flux_rotated.rhov.x*sin_angle + Flux_rotated.rhov.y*cos_angle;
 
  Flux.zero_non_sol();
 
  return (Flux);
  
}

LESPremixed2D_cState FluxRoe_n(const LESPremixed2D_cState &Ul,
			       const LESPremixed2D_cState &Ur,
			       const Vector2D &norm_dir, 
			       const int &Preconditioning,
			       const int &flow_type_flag,
			       const double &delta_n ) {
    return (FluxRoe_n(Ul.W(), Ur.W(), norm_dir,
		      Preconditioning,flow_type_flag,delta_n));
}


/**********************************************************************
 * Routine: FluxAUSMplus_up (Liou's updated Advection Upstream        * 
 *                           Splitting Method flux function for all   *
 *                           speeds,  x-direction)                    *
 *                                                                    *
 * This function returns the intermediate state solution flux for the *
 * x-direction given left and right solution states by using the      *
 * AUSM+-up (updated AUSM scheme) approximation for the fluxes.  See  *
 * M.-S. Liou (J. Comp. Physics 2006).                                *
 *                                                                    *
 **********************************************************************/
LESPremixed2D_cState FluxAUSMplus_up(const LESPremixed2D_pState &Wl,
				     const LESPremixed2D_pState &Wr) {

  LESPremixed2D_cState Flux, Convected_Quantities;
  double beta = 0.125, sigma = 0.75, Kp =0.25, Ku = 0.75;
  double alpha, rhohalf, mass_flux_half;
  double ahalf, Ml, Mr, Mplus, Mminus, Mhalf, pplus, pminus, phalf;
  //double al, ar, atilde_l, atilde_r;

  // Determine the intermediate state sound speed and density:
  // al = sqrt(Wl.H()/Wl.rho)*sqrt(TWO*(Wl.g() - ONE)/(Wl.g() + ONE));
//   ar = sqrt(Wr.H()/Wr.rho)*sqrt(TWO*(Wr.g() - ONE)/(Wr.g() + ONE));
//   atilde_l = sqr(al)/max(al, Wl.v.x);    
//   atilde_r = sqr(ar)/max(ar, -Wr.v.x);
//   ahalf = min(atilde_l, atilde_r);

  ahalf = HALF*(Wl.amodified() + Wr.amodified());
  rhohalf = HALF*(Wl.rho + Wr.rho); 
  

  // Determine the left and right state Mach numbers based on the
  // intermediate state sound speed:
  Ml = Wl.v.x/ahalf;
  Mr = Wr.v.x/ahalf;

  // Determine the reference Mach number, scaling function and coefficient
  double M2_bar, M2_ref, fa;
  M2_bar = (Wl.v.x*Wl.v.x + Wr.v.x*Wr.v.x)/(TWO*ahalf*ahalf);
  M2_ref = min(ONE, max(M2_bar, Wl.Mref*Wl.Mref));
  if (M2_ref > ONE || M2_ref < 0.0) cout << "\nM2_ref out of range, FluxAUSMplus_up.";
  //fa = sqrt(M2_ref)*(TWO - sqrt(M2_ref));
  fa = sqrt(sqr(ONE - M2_ref)*M2_bar + FOUR*M2_ref)/(ONE + M2_ref);
  if (fa > ONE || fa <= ZERO) cout << "\nfa out of range, FluxAUSMplus_up.";
  alpha = (3.0/16.0)*(-4.0 + 5.0*fa*fa);
  if (alpha < (-3.0/4.0)  ||  alpha > (3.0/16.0)) cout << "\nalpha out of range, FluxAUSMplus_up.";


  // Determine the left state split Mach number:
  if (fabs(Ml) >= ONE) {
    Mplus = Mplus_1(Ml);
    pplus = Mplus_1(Ml)/Ml;
  } else {
    Mplus = Mplus_2(Ml) * (1.0 - 16.0*beta*Mminus_2(Ml));
    pplus = Mplus_2(Ml) * ((2.0 - Ml) - 16.0*alpha*Ml*Mminus_2(Ml));
  }


  // Determine the right state split Mach number:
  if (fabs(Mr) >= ONE) {
    Mminus = Mminus_1(Mr);
    pminus = Mminus_1(Mr)/Mr;        
  } else {
    Mminus = Mminus_2(Mr) * (1.0 + 16.0*beta*Mplus_2(Mr));
    pminus = Mminus_2(Mr) * ((-2.0 - Mr) + 16.0*alpha*Mr*Mplus_2(Mr));
  } 


  // Determine the intermediate state Mach number, pressure and mass flux:
  Mhalf = Mplus + Mminus
    - (Kp/fa)*max((ONE - sigma*M2_bar), ZERO)*(Wr.pmodified() - Wl.pmodified())/(rhohalf*ahalf*ahalf);

  phalf = pplus*Wl.pmodified() + pminus*Wr.pmodified()
    - Ku*pplus*pminus*TWO*rhohalf*(fa*ahalf)*(Wr.v.x - Wl.v.x);

  mass_flux_half = (Mhalf > ZERO) ? ahalf*Mhalf*Wl.rho : ahalf*Mhalf*Wr.rho; 


  // Determine the intermediate state convective solution flux:
  if (mass_flux_half  > ZERO) {
    Convected_Quantities.rho = ONE;
    Convected_Quantities.rhov.x = Wl.v.x; 
    Convected_Quantities.rhov.y = Wl.v.y; 
    Convected_Quantities.E = Wl.H()/Wl.rho;

    if(Wl.nscal){
      for(int i=0; i<Wl.nscal; ++i) Convected_Quantities.rhoscalar[i] = Wl.scalar[i];
    }

    for(int i=0; i<Wl.ns; ++i){
      Convected_Quantities.rhospec[i].c = Wl.spec[i].c;
    }
    
  } else {
    Convected_Quantities.rho = ONE;
    Convected_Quantities.rhov.x = Wr.v.x; 
    Convected_Quantities.rhov.y = Wr.v.y; 
    Convected_Quantities.E = Wr.H()/Wr.rho;

    if(Wr.nscal){
      for(int i=0; i<Wr.nscal; ++i) Convected_Quantities.rhoscalar[i] = Wr.scalar[i];
    }

    for(int i=0; i<Wr.ns; ++i){
      Convected_Quantities.rhospec[i].c = Wr.spec[i].c;
    }

  } //end if

  Flux = mass_flux_half*Convected_Quantities;

  // Add the pressure contribution to the intermediate state solution flux:
  Flux[2] += phalf;


  // Return solution flux.
  return Flux;

}

LESPremixed2D_cState FluxAUSMplus_up(const LESPremixed2D_cState &Ul,
				     const LESPremixed2D_cState &Ur) {
  return FluxAUSMplus_up(Ul.W(),Ur.W());
}

/**********************************************************************
 * Routine: FluxAUSMplus_up_n (M.-S. Liou's Advection Upstream        *
 *                             Splitting Method flux function for     *
 *                             all speeds, n-direction)               *
 *                                                                    *
 * This function returns the intermediate state solution flux for an  *
 * arbitrary direction defined by a unit normal vector in the         *
 * direction of interest, given left and right solution states.  The  *
 * problem is solved by first applying a frame rotation to rotate the *
 * problem to a local frame aligned with the unit normal vector and   *
 * then by using the AUSM+-up approximation to specify the            * 
 * intermediate  state in terms of the rotated solution states.       *
 * See M.-S. Liou (J. Comp. Physics 2006).                            *
 **********************************************************************/
LESPremixed2D_cState FluxAUSMplus_up_n(const LESPremixed2D_pState &Wl,
				       const LESPremixed2D_pState &Wr,
				       const Vector2D &norm_dir) {

  LESPremixed2D_pState Wl_rotated(Wl), Wr_rotated(Wr);
  LESPremixed2D_cState Flux_rotated;

  /* Determine the direction cosine's for the frame rotation. */
   double cos_angle = norm_dir.x; 
   double sin_angle = norm_dir.y;

  /* Apply the frame rotation and evaluate left and right
     solution states in the local rotated frame defined
     by the unit normal vector. */
//   Wl_rotated.Copy(Wl);
//   Wr_rotated.Copy(Wr);
   
  Wl_rotated.v.x = Wl.v.x*cos_angle + Wl.v.y*sin_angle;
  Wl_rotated.v.y = - Wl.v.x*sin_angle +  Wl.v.y*cos_angle;

  Wr_rotated.v.x = Wr.v.x*cos_angle +  Wr.v.y*sin_angle;
  Wr_rotated.v.y = - Wr.v.x*sin_angle + Wr.v.y*cos_angle;
  
  // Evaluate the intermediate state solution flux in the rotated frame.
  Flux_rotated = FluxAUSMplus_up(Wl_rotated,Wr_rotated); 

  // Rotate back to the original Cartesian reference frame and return
  // the solution flux.
  LESPremixed2D_cState Flux(Flux_rotated);

  Flux.rhov.x = Flux_rotated.rhov.x*cos_angle - Flux_rotated.rhov.y*sin_angle;
  Flux.rhov.y = Flux_rotated.rhov.x*sin_angle + Flux_rotated.rhov.y*cos_angle;

  Flux.zero_non_sol();

  // Return the solution flux.  
  return Flux;

}

LESPremixed2D_cState FluxAUSMplus_up_n(const LESPremixed2D_cState &Ul,
				       const LESPremixed2D_cState &Ur,
				       const Vector2D &norm_dir) {
  return FluxAUSMplus_up_n(Ul.W(),Ur.W(),norm_dir);
}



/************************************************************************
 ************************************************************************
 ************** VISCOUS RECONSTRUCTION  ********************************* 
 ************************************************************************
 ************************************************************************/

/**********************************************************************
 * Routine: ViscousFluxArithmetic_n                                   *
 *                                                                    *
 * This function returns the intermediate state solution viscous flux *
 * calculated by the arithmetic mean of the cell-centred flux terms   *
 * of the neighbouring cells.                                         *
 *                                                                    *
 **********************************************************************/
LESPremixed2D_cState Viscous_FluxArithmetic_n(const LESPremixed2D_cState &Ul,
					      const LESPremixed2D_pState &dWdx_l,
					      const LESPremixed2D_pState &dWdy_l,
					      const LESPremixed2D_cState &Ur,
					      const LESPremixed2D_pState &dWdx_r,
					      const LESPremixed2D_pState &dWdy_r,
					      const int Flow_Type,
					      const int Axisymmetric,
					      const Vector2D X,
					      const Vector2D &norm_dir) {

  LESPremixed2D_cState Fx, Fy;
  Vector2D i = Vector2D(1,0), j = Vector2D(0,1);

  Fx = HALF*(Ul.Viscous_Flux_x(dWdx_l, dWdy_l, Flow_Type, Axisymmetric, X) + Ur.Viscous_Flux_x(dWdx_r, dWdy_r, Flow_Type, Axisymmetric, X));
  Fy = HALF*(Ul.Viscous_Flux_y(dWdx_l, dWdy_l, Flow_Type, Axisymmetric, X) + Ur.Viscous_Flux_y(dWdx_r, dWdy_r, Flow_Type, Axisymmetric, X));

  return (Fx*dot(i,norm_dir) + Fy*dot(j,norm_dir));

}

/**********************************************************************
 * Routine: Viscous_Flux_n                                            *
 *                                                                    *
 * This function returns the intermediate state solution viscous flux *
 * given the primitive variable solution state and the gradients of   *
 * the primitive variables.                                           *
 *                                                                    *
 **********************************************************************/
LESPremixed2D_cState Viscous_Flux_n(LESPremixed2D_pState &W,
				    const LESPremixed2D_pState &dWdx,
				    const LESPremixed2D_pState &dWdy,
				    const int Flow_Type,
				    const int Axisymmetric,
				    const Vector2D X,
				    const Vector2D &norm_dir) {

  LESPremixed2D_cState U;
  double Temperature, Rmix;
  Vector2D grad_T;
  Vector2D i = Vector2D(1,0), j = Vector2D(0,1);

  //Molecular transport properties
  Temperature = W.T();
  Rmix = W.Rtot();
 
  //Conserved variables from primitive
  U.rho = W.rho;
  U.rhov = W.rhov();
  U.E = W.E();
  for(int i=0; i<W.nscal; ++i) U.rhoscalar[i] = W.rho*W.scalar[i];
  
  for(int i=0; i<W.ns; ++i){
    U.rhospec[i].c = W.rho*W.spec[i].c;
  } 
 
  //Temperature gradient
  //dT/dx = 1/rho*R *( dP/dx - P/rho * drho/dx)
  grad_T.x = (ONE/(W.rho*Rmix))*(dWdx.p - (W.p/W.rho)*dWdx.rho);
  grad_T.y = (ONE/(W.rho*Rmix))*(dWdy.p - (W.p/W.rho)*dWdy.rho);
 
  // Molecular (Laminar) diffusion of species
  // for each of the "n" species
  for( int k=0; k<U.ns; ++k){
    /***************** Diffusion coefficients **********************/
    // using global Schmidt number relation Scs = mu/rho*Ds
#ifdef THICKENED_FLAME_ON
    U.rhospec[k].diffusion_coef = (W.mu()/(W.flame.WF*W.flame.TF))/U.Schmidt[k];
#else
    U.rhospec[k].diffusion_coef = W.mu()/U.Schmidt[k];
#endif
    /***************** mass fraction gradients *********************/
    U.rhospec[k].gradc.x = U.rho * dWdx.spec[k].c;
    U.rhospec[k].gradc.y = U.rho * dWdy.spec[k].c;
  }
  
  //Molecular (laminar) stress tensor
  Tensor2D strain_rate = W.Strain_Rate(dWdx,dWdy, Flow_Type,Axisymmetric, X); 
  W.Laminar_Stress(strain_rate);
  U.tau = W.tau;

  //Molecular (laminar) heat flux
  //Thermal conduction, q = - kappa * grad(T)
  U.qflux = - W.kappa()*grad_T;
  //Thermal diffusion, q -= rho * sum ( hs * Ds *gradcs)
  U.qflux -= U.rho*U.thermal_diffusion(Temperature);  
  
  //Turbulent heat flux
  //Thermal conduction, q = - kappa * grad(T)
  if (Flow_Type == FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY ||
      Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) {
 
    double mut = W.mu_t(strain_rate);
    double Dm_t = W.Dm_turb(mut);

    U.theta = - mut*W.Cp()/W.Pr_turb()*grad_T;
    //Thermal Diffusion, q -= rho * sum ( hs * Ds *gradcs)   
    for (int k=0; k<U.ns; ++k) {
      U.theta -= Dm_t*U.rhospec[k].gradc*
	(/*U.specdata[k].Enthalpy(Temperature)+*/U.specdata[k].Heatofform());
    }
    
    //SFS stresses 
    W.SFS_Stress(strain_rate); 
    U.lambda = W.lambda;

  } 

  return (U.Viscous_Flux_x(dWdx, dWdy, Flow_Type, Axisymmetric, X)*dot(i,norm_dir) 
	  + U.Viscous_Flux_y(dWdx, dWdy, Flow_Type, Axisymmetric, X)*dot(j,norm_dir));

}

/**********************************************************************
 * Routine: WallShearStress                                           *
 *                                                                    *
 * This routine computes and returns the shear stress at a wall.      *
 *                                                                    *
 **********************************************************************/
double WallShearStress(const LESPremixed2D_pState &W1,
		       const Vector2D &X1,
		       const Vector2D &X2,
		       const Vector2D &X3,
		       const Vector2D &norm_dir) {

  double l21, l32, l13, A, dWdn;
  Vector2D n21, n32, n13;
  LESPremixed2D_pState W2, W3, W_face, dWdx, dWdy;

  // Initialze W2 and W3.
  W2.Vacuum(); W2.rho = W1.rho; W2.p = W1.p;
  W3.Vacuum(); W3.rho = W1.rho; W3.p = W1.p;
  for(int i=0; i<W1.ns; ++i){
    W2.spec[i].c = W1.spec[i].c;
    W3.spec[i].c = W1.spec[i].c;
  }

  // Determine the lengths and normals of te faces and the 
  // areas of the regions of Green-Gauss integration.
  l21 = abs(X2-X1); n21 = Vector2D((X2.y-X1.y),-(X2.x-X1.x))/l21;
  l32 = abs(X3-X2); n32 = Vector2D((X3.y-X2.y),-(X3.x-X2.x))/l32;
  l13 = abs(X1-X3); n13 = Vector2D((X1.y-X3.y),-(X1.x-X3.x))/l13;
  A = HALF*((X2-X1)^(X3-X1));
  // Compute Green-Gauss integration on left triangle.
  W_face = HALF*(W2+W1)*l21;
  dWdx = W_face*n21.x;
  dWdy = W_face*n21.y;
  W_face = HALF*(W3+W2)*l32;
  dWdx += W_face*n32.x;
  dWdy += W_face*n32.y;
  W_face = HALF*(W1+W3)*l13;
  dWdx += W_face*n13.x;
  dWdy += W_face*n13.y;
  dWdx = dWdx/A;
  dWdy = dWdy/A;
  // Determine the normal gradient.
  dWdn = dWdy.v.x;//dot(Vector2D(dWdx.v.x,dWdy.v.y),norm_dir);

  // Return the wall shear stress.
  return W1.mu()*dWdn;

}
