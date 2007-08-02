/********************** Chem2DState.cc ******************************
  Constructors for Chem2Dstate which handles all the physical
  variables and mixture rules associated with multi-species
  chemically reacting flows.

   assosicated files:
           Chem2DState.h    

*********************************************************************/

#ifndef _CHEM2D_STATE_INCLUDED
#include "Chem2DState.h"
#endif // _CHEM2D_STATE_INCLUDED   

/***********************************************************/
//Static members initialization
int Chem2D_pState::ns = 1; 
int Chem2D_pState::NUM_VAR_CHEM2D = NUM_CHEM2D_VAR_SANS_SPECIES; 
NASARP1311data* Chem2D_pState::specdata=NULL;
Reaction_set Chem2D_pState::React;
double Chem2D_pState::low_temp_range = 200.0;
double Chem2D_pState::high_temp_range = 300.0;
int Chem2D_pState::debug_level = 0;            // debug level flag (0=none, 1,2,..)
double Chem2D_pState::Mref=0.5;
double* Chem2D_pState::Schmidt=NULL;
//total duplication -- find a better way   
int Chem2D_cState::ns = 1; 
int Chem2D_cState::NUM_VAR_CHEM2D = NUM_CHEM2D_VAR_SANS_SPECIES;   
NASARP1311data* Chem2D_cState::specdata=NULL;   
double Chem2D_cState::low_temp_range = 200.0;
double Chem2D_cState::high_temp_range = 300.0;

int Chem2D_pState::flow_type = FLOWTYPE_LAMINAR;
int Chem2D_cState::flow_type = FLOWTYPE_LAMINAR;

//k-omega Turbulence model coefficients
double Chem2D_pState::alpha = FIVE/NINE; //13.0/15.0; 
double Chem2D_pState::sigma = HALF;
double Chem2D_pState::sigma_star = HALF;
double Chem2D_pState::beta = THREE/40.00; //9.0/125.0;
double Chem2D_pState::f_beta = ONE;
double Chem2D_pState::beta_star = NINE/100.0;
double Chem2D_pState::f_beta_star = ONE;
double Chem2D_pState::Coeff_edm = ZERO;

double Chem2D_cState::alpha = FIVE/NINE; //13.0/15.0; 
double Chem2D_cState::sigma = HALF;
double Chem2D_cState::sigma_star = HALF;
double Chem2D_cState::beta = THREE/40.00; //9.0/125.0;
double Chem2D_cState::f_beta = ONE;
double Chem2D_cState::beta_star = NINE/100.0;
double Chem2D_cState::f_beta_star = ONE;
double Chem2D_cState::Coeff_edm = ZERO;

int Chem2D_cState::debug_level = 0;           // debug level flag (0=none, 1,2,..)
double Chem2D_cState::Mref=0.5;
double* Chem2D_cState::Schmidt=NULL;

/**********************************************************************
 * Chem2D_pState -- Create storage and assign turbulence static       *
 *                  variables.                                        *
 **********************************************************************/
// Turbulent boundary-layer constants:
double Chem2D_pState::yplus_o = 10.0;
double Chem2D_pState::C = 5.0;
double Chem2D_pState::von_karman = 0.41;
double Chem2D_pState::yplus_sublayer = 2.5;
double Chem2D_pState::yplus_buffer_layer = 30.0;
double Chem2D_pState::yplus_outer_layer = 100.0;

/**********************************************************************
 * Chem2D_cState -- Create storage and assign turbulence static       *
 *                  variables.                                        *
 **********************************************************************/
// Turbulent boundary-layer constants:
double Chem2D_cState::yplus_o = 10.0;
double Chem2D_cState::C = 5.0;
double Chem2D_cState::von_karman = 0.41;
double Chem2D_cState::yplus_sublayer = 2.5;
double Chem2D_cState::yplus_buffer_layer = 30.0;
double Chem2D_cState::yplus_outer_layer = 100.0;

/***********************************************************/
/************* set_species_data ***************************
**********************************************************/
//set Global data for Species (STATIC, ie. so only call once! in Chem2Dinput)
void Chem2D_pState::set_species_data(const int &n,const string *S,const char *PATH,
				     const int &debug, const double &Mr, const double* Sc){ 
  //Deallocate_static();
  Deallocate(); //Clean up memory before changing ns

  ns =n;
  NUM_VAR_CHEM2D = ns + NUM_CHEM2D_VAR_SANS_SPECIES;
  //read in NASA data for each species to be used 
  specdata = new NASARP1311data[ns]; 
  Schmidt = new double[ns];
  for(int i=0; i<ns; i++){
    //overwrite default data  
    specdata[i].Getdata(S[i],PATH);  
    Schmidt[i] = Sc[i];
  } 
  
  //set data temperature ranges for mixture
  Temp_low_range(); 
  Temp_high_range(); 
  
  //Set Debug Information level
  debug_level = debug;
  Mref = Mr;

  //setup initial array for mass fractions
  set_initial_values(); 
}   

//exact same thing with conserved variables, 
//should just point to pState data but for now duplication will work 
void Chem2D_cState::set_species_data(const int &n, const string *S, const char *PATH,
				     const int &debug, const double &Mr, const double* Sc){ 
  //Deallocate_static();
  Deallocate(); //Clean up memory before changing ns
 
  ns =n; 
  NUM_VAR_CHEM2D = ns + NUM_CHEM2D_VAR_SANS_SPECIES;
  //read in NASA data for each species to be used
  specdata = new NASARP1311data[ns];
  Schmidt = new double[ns];
  for(int i=0; i<ns; i++){
    //overwrite default data  
    specdata[i].Getdata(S[i],PATH);
    Schmidt[i] = Sc[i];  
  }  

  //set data temperature ranges for mixture
  Temp_low_range();
  Temp_high_range();

  //Set Debug Information level
  debug_level = debug;
  Mref = Mr;
  
  //setup initial array for mass fractions
  set_initial_values();
}      

/**************************************************************************
********************* CHEM2D_PSTATE CONSTRUCTORS **************************
***************************************************************************/

/***********************************************************
 ****************** Mixture Rules **************************
 ***********************************************************/

/**************************************************
  mixture molecular mass (kg/mol)
***************************************************/
double Chem2D_pState::Mass() const{
  // = 1 / sum(mass fration/mol_mass)
  double sum = 0.0;
  for(int i=0; i<ns; i++){
    sum += spec[i].c/specdata[i].Mol_mass();
  }
  return 1.0/sum;
}

/**************************************************
  mixture gas constant  J/(kg*K)
***************************************************/
double Chem2D_pState::Rtot(){
  // = sum ( mass fraction * species gas constant)
  double sum = 0.0;
  for(int i=0; i<ns; i++){
    sum += spec[i].c * specdata[i].Rs();
  }
  return sum;
}

double Chem2D_pState::Rtot() const{
  // = sum ( mass fraction * species gas constant)
  double sum = 0.0;
  for(int i=0; i<ns; i++){
    sum += spec[i].c * specdata[i].Rs();
  }
  return sum;
}

/**************************************************
  mixture Heat Capacity (const pressure) J/(kg*K)
***************************************************/
double Chem2D_pState::Cp(void) const{
  // = sum ( mass fraction * species Cp) 
  double Temp = T();
  double sum = 0.0;
  for(int i=0; i<ns; i++){
    sum += spec[i].c*specdata[i].HeatCapacity_p(Temp);
  }
  return sum;
}

double Chem2D_pState::Cp(const double& TEMP) const{
  // = sum ( mass fraction * species Cp) 
  double sum = 0.0;
  for(int i=0; i<ns; i++){
    sum += spec[i].c*specdata[i].HeatCapacity_p(TEMP);
  }
  return sum;
}

/**************************************************
  mixture Heat Capacity (const volume) J/(kg*K)
***************************************************/
double Chem2D_pState::Cv(void) const{
  // = sum ( mass fraction * species Cv)  
  double Temp = T();
  double sum = 0.0;
  for(int i=0; i<ns; i++){
    sum += spec[i].c*specdata[i].HeatCapacity_v(Temp);
  }
  return sum;
}

/**************************************************
  mixture Heat Ratio gamma J/(kg*K)
***************************************************/
double Chem2D_pState::g(void) const{
  // = Cp / Cv  
  return Cp()/Cv();
}

/**************************************************
  Specific Internal Energy
***************************************************/
// etotal = sensible & chemical
double Chem2D_pState::e(void) const{
  // = sum (mass fraction * species e) 
  double sum = 0.0;
  double Temp = T();
  for(int i=0; i<ns; i++){ //(Enthalpy(Temp) - (R/mol_mass)*Temp)
    sum += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform()
    		      - specdata[i].Rs()*Temp);
  }
  return sum;
}

// reference internal energy e + heat of formation - offset
double Chem2D_pState::eref(void)const{
  // = sum (mass fraction * species e) 
  double sum = 0.0;
  double Temp = T();
  for(int i=0; i<ns; i++){ //(Enthalpy(Temp) - (R/mol_mass)*Temp)
    sum += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform()
		      	 - specdata[i].DeltaHref() - specdata[i].Rs()*Temp);
  }
  return sum;
}

// internal energy with no heat of formation included (sensible)
double Chem2D_pState::es(void) const{
  // = sum (mass fraction * species e) 
  double sum = 0.0;
  double Temp = T();
  for(int i=0; i<ns; i++){ //(Enthalpy(Temp) - (R/mol_mass)*Temp)
    sum += spec[i].c*(specdata[i].Enthalpy(Temp) - specdata[i].Rs()*Temp);
  }
  return sum;
}

/************************************************** 
   Specific absolute enthalpy = 
          mass fraction * (hsensible + heatofform)
***************************************************/
double Chem2D_pState::h(void) const{
  // = sum (mass fraction * species h) 
  double sum = 0.0;  
  double Temp = T();
  for(int i=0; i<ns; i++){
    sum += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform());
  }
  return sum;
}

double Chem2D_pState::h(const double &Temp) const{
  // = sum (mass fraction * species h) 
  double sum = 0.0;  
  for(int i=0; i<ns; i++){
    sum += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform());
  }
  return sum;
}


double Chem2D_pState::href(void)const{
  // = sum (mass fraction * species h) 
  double sum = 0.0;  
  double Temp = T();
  for(int i=0; i<ns; i++){ 
    sum += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform()
		      - specdata[i].DeltaHref() );
  }
  return sum;
}

double Chem2D_pState::hs(void) const{
  // = sum (mass fraction * species h) 
  double sum = 0.0;  
  double Temp = T();
  for(int i=0; i<ns; i++){
    sum += spec[i].c*(specdata[i].Enthalpy(Temp));
  }
  return sum;
}

double Chem2D_pState::hs(const double &Temp) const{
  // = sum (mass fraction * species h) 
  double sum = 0.0;  
  for(int i=0; i<ns; i++){
    sum += spec[i].c*(specdata[i].Enthalpy(Temp));
  }
  return sum;
}

/**************************************************
   Derivative of specific enthalpy dh/dT
   actually is just Cp as Cp = (dh/dT)_p
***************************************************/
double Chem2D_pState::hprime() const{
 double sum = 0.0; 
 double Temp = T();
 for(int i=0; i<ns; i++){
   sum += spec[i].c*specdata[i].Enthalpy_prime(Temp);
 }
 return (sum);
}

double Chem2D_pState::hprime(double &Temp) const{
 double sum = 0.0;  
 for(int i=0; i<ns; i++){
   sum += spec[i].c*specdata[i].Enthalpy_prime(Temp);
 }
 return (sum);
}

/**************************************************
  Total Energy
***************************************************/
double Chem2D_pState::E(void) const{
  // E = rho*(e + velocity^2 +k)  
  return (rho*(e() + HALF*v.sqr()+ k)); 
}

/**************************************************
  Total Enthalpy
***************************************************/
double Chem2D_pState::H(void) const{
  // H = h + velocity^2 
  double TotalEnthalpy = ZERO;
  
  
  TotalEnthalpy = (rho*(h() + HALF*v.sqr() +k));
  

  return (TotalEnthalpy);
}

double Chem2D_pState::Hs(void) const{
  // H = h + velocity^2+k
  return (rho*(hs() + HALF*v.sqr() + k));
}

/**************************************************
  Viscosity 
  using Wilke [1950] formulation
***************************************************/
double Chem2D_pState::mu() const{
  double sum =0.0; 
  double Temp = T();
  double  *vis = new double[ns];

  for(int i=0; i<ns; i++){
    double phi = 0.0;
    for (int j=0; j<ns; j++){
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
  
  delete[] vis;
  return sum;
}

/****************************************************
  Thermal Conductivity - Mason & Saxena (1958)  W/(m*K)
****************************************************/
double Chem2D_pState::kappa(void) const{
  double sum = 0.0;  
  double Temp = T();
  double  *vis = new double[ns];

  for(int i=0; i<ns; i++){
    double phi = 0.0;
    for (int j=0; j<ns; j++){
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
//   for (int j=0; j<ns; j++){
//     one += spec[i].c*specdata[i].ThermalConduct(Temp);
//     two += spec[i].c/specdata[i].ThermalConduct(Temp);
//   }
//   sum = HALF*(one + ONE/two);

  delete[] vis;
  return sum;
}

/**************************************************
  Turbulence model related parameters
***************************************************/
double Chem2D_pState::eddy_viscosity(void) const{
  return (rho*k/omega);
}

double Chem2D_pState::Pr_turb(void) const{
  return (0.9); 
}

double Chem2D_pState::Sc_turb(void) const{
  return (1.1);
}

double Chem2D_pState::Kappa_turb(void) const{
  return (eddy_viscosity()*Cp()/Pr_turb()); 
}

double Chem2D_pState::Dm_turb(void) const{
  return (eddy_viscosity()/(rho*Sc_turb()));
}     

double Chem2D_pState::omega_sublayer_BC(const double &y) const {
  return (SIX*mu()/(rho*beta*y*y));
}

/**************************************************
  polytropic heat ratio mixture gamma J/(kg*K)
  assuming T=200K as the temperature.
***************************************************/
double Chem2D_pState::gamma_guess(void) const{  
  double sum1 = ZERO;     double sum2 = ZERO;
  double gamma_s = ZERO;  double Temp = 200.0;

  for(int i=0; i<ns; i++){
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
double Chem2D_pState::T(double &h_s) const{
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
      if(T >= Tmax) T = HALF*(Tmax - Tmin);	
      //Bisection
    } else {
      T = HALF*(Tmax - Tmin);
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
    if(debug_level){ 
      cout<<"\nTemperature didn't converge in Chem2D_cState::T(void)";
      cout<<" with polytopic Tguess "<<Tguess<<", or lower than Tmin "
	  <<low_temp_range<<" using "<<T;
    }
  }
  return T;
}
    
/****************************************************
  Speed of sound using 
  a^2 = dip/dirho + p/rho^2( die/dip)^-1
  from eigenvalue analysis using e =f(p,rho)
****************************************************/
double Chem2D_pState::a(void){
  double sum;
  double RTOT= Rtot();
 
  sum = (p/rho)*(RTOT/( hprime() - RTOT) + ONE);
  //could also just call sqrt(g()*Rtot()*T());
  return sqrt(sum);
}

double Chem2D_pState::a(void) const{
  double sum;
  double RTOT= Rtot();
  sum = (p/rho)*(RTOT/( hprime() - RTOT) + ONE);
  //could also just call sqrt(g()*Rtot()*T());
  return sqrt(sum);
}

// /******************************************************
//  Calculating the thermal diffusion component of 
//  the heat flux vector (qflux)
//   sum( hs * Ds * grad cs)
// *******************************************************/
Vector2D Chem2D_pState::thermal_diffusion(void) const{
  Vector2D sum;
  sum.zero();
  double Temp = T();
  //problems with Species overloaded operators
  for(int i=0; i<ns; i++){ 
    sum  +=  (specdata[i].Enthalpy(Temp) + specdata[i].Heatofform())
      * spec[i].diffusion_coef * spec[i].gradc;
  }
  return sum;
}
/*******************************************************************
 ***************** FLUXES ******************************************
 *******************************************************************/
/********************************************************
 * Chem2D_pState::Fx -- Inviscid flux (x-direction).   *
 ********************************************************/
Chem2D_cState Chem2D_pState::Fx(void) const{
  Chem2D_cState Temp;
  
  Temp.rho = rho*v.x;
  Temp.rhov.x = rho*sqr(v.x) + p;
  Temp.rhov.y = rho*v.x*v.y;
  Temp.E = v.x*H();
  
  Temp.rhok = rho*k*v.x;
  Temp.rhoomega = rho*omega*v.x;
  
  //multispecies transport
  for(int i=0; i<ns; i++){
    Temp.rhospec[i].c = rho*v.x*spec[i].c;
  }
  
  return (Temp);
}

/*******************************************************************
 ***************** EIGENVALUES *************************************
 *******************************************************************/

/************************************************************
 * Chem2D_pState::lambda -- Eigenvalue(s) (x-direction).    *
 ************************************************************/
Chem2D_pState Chem2D_pState::lambda_x(void) const {
  double c = a();
  Chem2D_pState Temp;
  Temp.rho = v.x - c;
  Temp.v.x = v.x;
  Temp.v.y = v.x;
  Temp.p = v.x + c;
  for(int i=0; i<ns; i++){
    Temp.spec[i].c = v.x;
  }
 
  Temp.k = v.x;
  Temp.omega = v.x;
  
  return (Temp);
}

/***********************************************************
Molecular viscosity is function of Temperature
This derivative is needed by Jacobian
 ***********************************************************/
double Chem2D_pState::dmudT(void) const{
  double sum =0.0; 
  double Temp = T();
  for(int i=0; i<ns; i++){
    double phi = 0.0;
    for (int j=0; j<ns; j++){
      phi += (spec[j].c / specdata[j].Mol_mass())*
	pow(1.0 + sqrt(specdata[i].dViscositydT(Temp)/specdata[j].dViscositydT(Temp))*
	    pow(specdata[j].Mol_mass()/specdata[i].Mol_mass(),0.25),2.0)/
       sqrt(8.0*(1.0 +specdata[i].Mol_mass()/specdata[j].Mol_mass()));
    }
    sum += (spec[i].c * specdata[i].dViscositydT(Temp) ) / 
      (specdata[i].Mol_mass() * phi);
  }  
  return sum;
}

double Chem2D_cState::dmudT(void) const{
  double sum =0.0; 
  double Temp = T();

  for(int i=0; i<ns; i++){
    double phi = 0.0;
    for (int j=0; j<ns; j++){
      phi += (rhospec[j].c/rho / specdata[j].Mol_mass())*
	pow(1.0 + sqrt(specdata[i].dViscositydT(Temp)/specdata[j].dViscositydT(Temp))*
	    pow(specdata[j].Mol_mass()/specdata[i].Mol_mass(),0.25),2.0)/
       sqrt(8.0*(1.0 +specdata[i].Mol_mass()/specdata[j].Mol_mass()));
    }
    sum += (rhospec[i].c/rho * specdata[i].dViscositydT(Temp) ) / 
      (specdata[i].Mol_mass() * phi);
  }  
  return sum;
}
/************************************************************
 * Low Mach Number Preconditioned Eigenvalue(s)             *
 * Chem2D_pState::lambda_preconditioned_x                   *
 ************************************************************/
Chem2D_pState Chem2D_pState::lambda_preconditioned_x(const double &MR2) const {

  Chem2D_pState NEW;
  double uprimed, cprimed;
  double c = a();
  u_a_precon(MR2*c*c,uprimed,cprimed);

  NEW.rho = uprimed - cprimed;
  NEW.v.x = v.x;
  NEW.v.y = v.x;
  NEW.p = uprimed + cprimed;
  for(int i=0; i<ns; i++){
    NEW.spec[i].c = v.x;
  }
  
  NEW.k = v.x;
  NEW.omega = v.x;
  
  
  return (NEW);
}

/*******************************************************************
 ***************** EIGENVECTORS ************************************
 *******************************************************************/
// Conserved Right Eigenvector -- (x-direction)
Chem2D_cState Chem2D_pState::rc_x(const int &index) const {
     assert( index >= 1 && index <= NUM_VAR_CHEM2D );
    if(index == 1){
      double c = a(); 
      return (Chem2D_cState(ONE, v.x-c, v.y, H()/rho-v.x*c, k, omega, spec));
    } else if(index == 2) {
      return (Chem2D_cState(ONE, v.x, v.y, H()/rho-Cp()*T(), k, omega, spec)); 
    } else if(index == 3) {
      return (Chem2D_cState(ZERO, ZERO, rho, rho*v.y, ZERO));
    } else if(index == 4) {
      double c = a(); 
      return (Chem2D_cState(ONE, v.x+c, v.y, H()/rho+v.x*c, k, omega, spec));
    } else if(index == 5) { 
      return (Chem2D_cState(ZERO, ZERO, ZERO, rho, rho, ZERO, ZERO) );
    } else if(index == 6) {
      return (Chem2D_cState(ZERO, ZERO, ZERO, ZERO, ZERO, rho, ZERO) ) ;
    } else{ 
      for(int i=7; i<=NUM_VAR_CHEM2D; i++){
	if(index == i){
	  Chem2D_cState NEW(ZERO);
	  double RTOT = Rtot();
	  double TEMP = p/(rho*RTOT);      
	  NEW.E = rho*(specdata[i-7].Enthalpy(TEMP) + specdata[i-7].Heatofform() 
		       - Cp(TEMP)*TEMP*specdata[i-7].Rs()/RTOT);
	  NEW.rhospec[i-7].c = rho;
	  return NEW;
	  break;
	}
      }
    }
}

// Primitive Left Eigenvector -- (x-direction)
Chem2D_pState Chem2D_pState::lp_x(const int &index) const {
     assert( index >= 1 && index <= NUM_VAR_CHEM2D );
    if(index == 1){
      double c = a(); 
      return (Chem2D_pState(ZERO, -HALF*rho/c, ZERO, HALF/(c*c), ZERO));
    } else if(index == 2) {
      double c = a(); 
      return (Chem2D_pState(ONE, ZERO, ZERO, -ONE/(c*c), ZERO));
    } else if(index == 3) {
      return  (Chem2D_pState(ZERO, ZERO, ONE, ZERO, ZERO));
    } else if(index == 4) {  
      double c = a(); 
      return (Chem2D_pState(ZERO, HALF*rho/c, ZERO, HALF/(c*c), ZERO));
    }else if(index == 5) {  
      return (Chem2D_pState(ZERO, ZERO, ZERO, ZERO, ONE,ZERO, ZERO));
    }else if(index == 6) {  
      return (Chem2D_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO));
    } else{ 
      for(int i=7; i<=NUM_VAR_CHEM2D; i++){
	if(index == i){
	  Chem2D_pState NEW(ZERO);
	  NEW.spec[i-7].c = ONE;
	  return NEW;
	  break;
	}
      }
    } 

}

/************************************************************
 ************** PRECONDITIONED EIGENVECTORS *****************
 ************************************************************/
// Conserved Right Eigenvector -- (x-direction)
Chem2D_cState Chem2D_pState::rc_x_precon(const int &index, const double &MR2) const {
  assert( index >= 1 && index <= NUM_VAR_CHEM2D );
  if(index == 1){
    double c = a(); 
    double uprimed,cprimed;
    u_a_precon(MR2*c*c,uprimed,cprimed);
    return (Chem2D_cState(ONE, 
			  (uprimed-cprimed)/MR2,
			  v.y,			
			  h()+HALF*(v.sqr()/MR2) - (v.x*cprimed)/MR2, k,omega,
			  spec));
   
  } else if(index == 2) {
    return (Chem2D_cState(ONE, v.x, v.y, (h()-Cp()*T()) + HALF*v.sqr(), k,omega,spec));
  
  } else if(index == 3) {
    return (Chem2D_cState(ZERO, ZERO, rho, rho*v.y,ZERO));
   
  } else if(index == 4) { 
    double c = a(); 
    double uprimed,cprimed;
    u_a_precon(MR2*c*c,uprimed,cprimed);
    return (Chem2D_cState(ONE,
			  (uprimed+cprimed)/MR2,
			  v.y, 
			  h()+HALF*(v.sqr()/MR2) + (v.x*cprimed)/MR2, k, omega,
			  spec));
   
  } else if(index == 5) { 
    return (Chem2D_cState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO) );
  } else if(index == 6) {
    return (Chem2D_cState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO) ) ;
  } else{ 
    for(int i=7; i<=NUM_VAR_CHEM2D; i++){
      if(index == i){
	Chem2D_cState NEW(ZERO);
	NEW.E = rho*(specdata[i-7].Enthalpy(T()) + specdata[i-7].Heatofform() 
		     - Cp()*T()*specdata[i-7].Rs()/Rtot());

	NEW.rhospec[i-7].c = rho;

	return NEW;
	break;
      }
    }
  }

}

// Primitive Left Eigenvector -- (x-direction)
Chem2D_pState Chem2D_pState::lp_x_precon(const int &index, const double &MR2) const {
  
  assert( index >= 1 && index <= NUM_VAR_CHEM2D );
  if(index == 1){
    double c = a();   
    double uprimed,cprimed;
    u_a_precon(MR2*c*c,uprimed,cprimed);
    return (Chem2D_pState(ZERO, 
			  -HALF*rho*MR2/cprimed, 
			  ZERO,
			  (-uprimed+cprimed + v.x)/(TWO*cprimed*c*c),k,omega,
			  ZERO));
  } else if(index == 2) {
    double c = a(); 
    return (Chem2D_pState(ONE, ZERO, ZERO, -ONE/(c*c),k,omega,ZERO));
  } else if(index == 3) {
    return  (Chem2D_pState(ZERO, ZERO, ONE, ZERO,ZERO));
  } else if(index == 4) {  
    double c = a(); 
    double uprimed,cprimed;
    u_a_precon(MR2*c*c,uprimed,cprimed);
    return (Chem2D_pState(ZERO, 
			  HALF*rho*MR2/cprimed, 
			  ZERO,
			  (uprimed+cprimed - v.x)/(TWO*cprimed*c*c), k, omega,
			  ZERO));
  }else if(index == 5) {  
    return (Chem2D_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO));
  }else if(index == 6) {  
    return (Chem2D_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO));
  } else{ 
    for(int i=7; i<=NUM_VAR_CHEM2D; i++){
      if(index == i){
	Chem2D_pState NEW(ZERO);
	NEW.spec[i-7].c = ONE;
	return NEW;
	break;
      }
    }
  } 
}

/*******************************************************************
 *******************************************************************
 *  Low Mach Number Preconditioner Based on Weiss & Smith (1995)   *
 *                                                                 *
 *******************************************************************
 *******************************************************************/
// For CFL Calculation
double Chem2D_pState::u_plus_aprecon(const double &u, 
				     const double &deltax) const {
  double Temp = T();
  double c = a();
  double UR2 = Mr2(deltax)*c*c;
  double alpha = HALF*( ONE - (ONE/(Rtot()*Temp) - ONE/(Cp()*Temp))*UR2);
  
  // uprime + cprime
  return ( u*(ONE - alpha) + sqrt(alpha*alpha*u*u + UR2) );
}

/************************************************************/
// For eignvalues and eigenvectors return preconditioned
// velocity and soundspeed ie. u' and c'
void Chem2D_pState::u_a_precon(const double &UR2, double &uprimed, double &cprimed) const{
  
  double Temp = T();
  double alpha = HALF*( ONE - (ONE/(Rtot()*Temp) - ONE/(Cp()*Temp))*UR2);

  uprimed = v.x*(ONE - alpha);
  cprimed = sqrt(alpha*alpha*v.x*v.x + UR2); 

} 

/************************************************************/
// as defined by E.Turkel (1999)
double Chem2D_pState::Mr2(const double &deltax) const {
  
  double c = a();
  double MR2 = min(max((v.sqr()/(c*c)),Mref*Mref),ONE);
  
  // Need deltax which is based on cell spacing 
  if (flow_type) {
    MR2 = pow(max(sqrt(MR2*c*c), mu()/(rho*deltax)),2.0)/(c*c);
  }

  return (MR2);
}

/************************************************************/
/************ ORIGINAL LAMINAR ONLY *************************/
/* (will be replaced by version commented out below when corrected) */
/************************************************************/
void Chem2D_pState::Low_Mach_Number_Preconditioner(DenseMatrix &P,
						   const double &deltax ) const{  
  double Temp = T();
  double Rmix = Rtot();
  double enthalpy = h();
  double CP = Cp();
  double c = a();
  double theta = (ONE/(Mr2(deltax)*c*c) + ONE/(CP*Temp));
 
  double phi = ZERO;   
  for(int j=0; j<ns-1; j++){   
    phi += spec[j].c*(specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() 
		      - CP*Temp*specdata[j].Rs()/Rmix);   
  }

  double alpha = theta*p/rho;
  double alpham1 = alpha - ONE;
  double Omega = (Rmix - CP)*p/(rho*Rmix);
  double beta = enthalpy - CP*p/(rho*Rmix) - phi;
  double V = HALF*v.sqr();
  P.zero();

  P(0,0) = (alpha*(beta-V)+V+Rmix*Temp-enthalpy+phi)/Omega;
  P(0,1) = v.x*alpham1/Omega;
  P(0,2) = v.y*alpham1/Omega;
  P(0,3) = -alpham1/Omega;
  P(1,0) = v.x*(beta-V)*alpham1/Omega;
  P(1,1) = v.x*v.x*alpham1/Omega+1.0;
  P(1,2) = v.x*v.y*alpham1/Omega;
  P(1,3) = -v.x*alpham1/Omega;
  P(2,0) = v.y*(beta-V)*alpham1/Omega;
  P(2,1) = v.x*v.y*alpham1/Omega;
  P(2,2) = v.y*v.y*alpham1/Omega+1.0;
  P(2,3) = -v.y*alpham1/Omega;
  P(3,0) = (enthalpy+V)*(beta-V)*alpham1/Omega;
  P(3,1) = v.x*(enthalpy+V)*alpham1/Omega;
  P(3,2) = v.y*(enthalpy+V)*alpham1/Omega;
  P(3,3) = -(alpha*(enthalpy+V)-V-Rmix*Temp-beta-phi)/Omega;

  //fixes so it can work without Turbulence !!!
  P(4,4) = ONE;
  P(5,5) = ONE;

  int NUM_VAR =   NUM_CHEM2D_VAR_SANS_SPECIES;

  //Multispecies
  for(int j=0; j<ns-1; j++){       
    double enth_j = specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() 
      - CP*Temp*specdata[j].Rs()/Rmix;
    P(0,j+NUM_VAR) = enth_j*alpham1/Omega;
    P(1,j+NUM_VAR) = v.x*enth_j*alpham1/Omega;
    P(2,j+NUM_VAR) = v.y*enth_j*alpham1/Omega;
    P(3,j+NUM_VAR) = enth_j*(V+enthalpy)*alpham1/Omega;	
    for(int i=0; i<ns-1; i++){ 
      if(i==j){ 
	P(i+NUM_VAR,0) = (spec[i].c)*(beta-V)*alpham1/Omega;
	P(i+NUM_VAR,1) = (spec[i].c)*v.x*alpham1/Omega;
	P(i+NUM_VAR,2) = (spec[i].c)*v.y*alpham1/Omega;
	P(i+NUM_VAR,3) = -(spec[i].c)*alpham1/Omega;
	//diagonal
	P(i+NUM_VAR,j+NUM_VAR) = spec[i].c*enth_j*alpham1/Omega+1.0;
      }
      else {
	P(i+NUM_VAR,j+NUM_VAR) = spec[i].c*enth_j*alpham1/Omega;
      }
    }       
  }
}

/************************************************************/
void Chem2D_pState::Low_Mach_Number_Preconditioner_Inverse(DenseMatrix &Pinv,	
							   const double &deltax ) const{  
  double Temp = T();
  double Rmix = Rtot();
  double enthalpy = h();
  double CP = Cp();
  double c = a();
  double theta = (ONE/(Mr2(deltax)*c*c) + ONE/(CP*Temp));  

  double phi = ZERO;
  for(int j=0; j<ns-1; j++){   
    phi += spec[j].c*(specdata[j].Enthalpy(Temp) 
		      + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/Rmix);   
  }

  double AA = p*(rho*Rmix-theta*p*CP);
  double BB = Rmix*rho*(theta*p-rho);
  double EE = HALF*v.sqr() - enthalpy + phi;
  double CC = EE + CP*Temp; 
  double DD = HALF*v.sqr() + enthalpy;
  Pinv.zero();    
  
  Pinv(0,0) = rho*Rmix/AA*(theta*p*EE-rho*CC+p);
  Pinv(0,1) = -v.x*BB/AA;
  Pinv(0,2) = -v.y*BB/AA;
  Pinv(0,3) = BB/AA;
  Pinv(1,0) = v.x*CC*BB/AA;
  Pinv(1,1) = rho*Rmix/AA*(p+rho*v.x*v.x-theta*p*(v.x*v.x+CP*Temp));
  Pinv(1,2) = -v.x*v.y*BB/AA;
  Pinv(1,3) = v.x*BB/AA;    
  Pinv(2,0) = v.y*CC*BB/AA;
  Pinv(2,1) = -v.x*v.y*BB/AA;
  Pinv(2,2) = rho*Rmix/AA*(p+v.y*v.y*rho-theta*p*(v.y*v.y+CP*Temp));
  Pinv(2,3) = v.y*BB/AA;  
  Pinv(3,0) = DD*CC*BB/AA;
  Pinv(3,1) = -v.x*DD*BB/AA;
  Pinv(3,2) = -v.y*DD*BB/AA;
  Pinv(3,3) = rho*Rmix/AA*(theta*p*(DD-CP*Temp)-rho*DD+p);

  Pinv(4,4) = ONE;
  Pinv(5,5) = ONE;
  
  int NUM_VAR =   NUM_CHEM2D_VAR_SANS_SPECIES;
  //Multispecies
  for(int j=0; j<ns-1; j++){   
    double enth_j =  specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() 
      - CP*Temp*specdata[j].Rs()/Rmix;
    Pinv(0,j+NUM_VAR) = -enth_j*BB/AA;
    Pinv(1,j+NUM_VAR) = -v.x*enth_j*BB/AA;
    Pinv(2,j+NUM_VAR) = -v.y*enth_j*BB/AA;
    Pinv(3,j+NUM_VAR) = -enth_j*BB*DD/AA;
    for(int i=0; i<ns-1; i++){  
      if(i==j){
	Pinv(i+NUM_VAR,0) = (spec[i].c)*CC*BB/AA;
	Pinv(i+NUM_VAR,1) = -(spec[i].c)*v.x*BB/AA;
	Pinv(i+NUM_VAR,2) = -(spec[i].c)*v.y*BB/AA;
	Pinv(i+NUM_VAR,3) = (spec[i].c)*BB/AA;
	//diagonal	
	Pinv(i+NUM_VAR,j+NUM_VAR) = 1.0 - spec[i].c*enth_j*BB/AA ;
      }
      else {
	Pinv(i+NUM_VAR,j+NUM_VAR) = -spec[i].c*enth_j*BB/AA;
      } 
    }   
  }  
}

// Still some problems with this form
// /*******************************************************************
//  *******************************************************************
//  *  Low Mach Number Preconditioner Based on Weiss & Smith (1995)   *
//  *                                                                 *
//  *******************************************************************
//  *******************************************************************/
// void Chem2D_pState::Low_Mach_Number_Preconditioner(DenseMatrix  &P,
// 						   const double &deltax ) const{ 
//   double Temp = T();
//   double Rmix = Rtot();
//   double enthalpy = h();
//   double CP = Cp();
//   double c = a();
//   double Theta = (ONE/(Mr2deltax)*c*c) + ONE/(CP*Temp));
 
//   double phi = ZERO;   
//   for(int j=0; j<ns-1; j++){   
//     phi += spec[j].c*(specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/Rmix);   
//   }
//   double V = HALF*v.sqr();
//   double Omega = (Rmix - CP)*p/(rho*Rmix);
//   double psi = Theta*p/rho; 
//   double chi =  enthalpy - CP*p/(rho*Rmix) - phi - V;
//   double coefficient = (psi -ONE)/Omega;
//   double H = (enthalpy+V);
 
//   P.zero();
 
//   P(0,0) = coefficient*chi +ONE;
//   P(0,1) = coefficient*v.x;
//   P(0,2) = coefficient*v.y;
//   P(0,3) =  -coefficient;

//   P(1,0) = coefficient*v.x*chi;
//   P(1,1) = coefficient*v.x*v.x+ONE;
//   P(1,2) = coefficient*v.x*v.y;
//   P(1,3) = -coefficient*v.x;

//   P(2,0) = coefficient*v.y*chi;
//   P(2,1) = coefficient*v.x*v.y;
//   P(2,2) = coefficient*v.y*v.y+ONE;
//   P(2,3) = -coefficient*v.y;

//   P(3,0) = coefficient*H*chi;
//   P(3,1) = coefficient*H*v.x;
//   P(3,2) = coefficient*H*v.y;
//   P(3,3) = -coefficient*H +ONE;

//   P(4,4) = coefficient*k+ONE;
//   P(5,5) = ONE;

//   if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
//       flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
//     P(0,4) = coefficient;
//     P(1,4) = coefficient*v.x;
//     P(2,4) = coefficient*v.y;
//     P(3,4) = coefficient*H;
//     //k-omega ...
//     P(4,0) = coefficient*k*chi;
//     P(4,1) = coefficient*k*v.x;
//     P(4,2) = coefficient*k*v.y;
//     P(4,3) = -coefficient*k;
  
//     P(5,0) = coefficient*omega*chi;
//     P(5,1) = coefficient*omega*v.x;
//     P(5,2) = coefficient*omega*v.y;
//     P(5,3) = -coefficient*omega;
//     P(5,4) = coefficient*omega;
 

//   }
//   int NUM_VAR =   NUM_CHEM2D_VAR_SANS_SPECIES;
//   //Multispecies (column entries)...
//   for(int j=0; j<ns-1; j++){ 
   
//     double enth_j = specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/Rmix;
//     P(0,j+NUM_VAR) = coefficient*enth_j;
//     P(1,j+NUM_VAR) = coefficient*v.x*enth_j;
//     P(2,j+NUM_VAR) = coefficient*v.y*enth_j;
//     P(3,NUM_VAR) =coefficient* H*enth_j;	
    
//     if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
// 	flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
   
//       P(4,j+NUM_VAR) = coefficient*k*enth_j;
//       P(5,j+NUM_VAR) = coefficient*omega*enth_j;
//     }
//   }
 
//   //Multispecies (row entries)
//   for(int i=0; i<ns-1; i++){ 
//     double enth_i = specdata[i].Enthalpy(Temp) + specdata[i].Heatofform() - CP*Temp*specdata[i].Rs()/Rmix;
//     P(i+NUM_VAR,0) = coefficient*(spec[i].c)*chi;
//     P(i+NUM_VAR,1) = coefficient*(spec[i].c)*v.x;
//     P(i+NUM_VAR,2) = coefficient*(spec[i].c)*v.y;
//     P(i+NUM_VAR,3) = -coefficient*(spec[i].c);
    
//     if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
// 	flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
     
//       P(i+NUM_VAR,4) = coefficient*spec[i].c;
//     }
    
//     for(int j=0; j<ns-1; j++){ 
//       if(j==i){
// 	//diagonal
// 	P(i+NUM_VAR,j+NUM_VAR) = coefficient*spec[i].c*enth_i+ONE;
//       }else{
// 	double enth_j = specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/Rmix;
// 	P(i+NUM_VAR,j+NUM_VAR) = coefficient*spec[i].c*enth_j;
//       }
//     }//multispecies  
    
//   }//end of row entries 


// }

// //Inverse of low mach number preconditioner
// void Chem2D_pState::Low_Mach_Number_Preconditioner_Inverse(DenseMatrix  &Pinv,	
// 							   const double &deltax ) const{ 
//   double Temp = T();
//   double Rmix = Rtot();
//   double enthalpy = h();
//   double CP = Cp();
//   double c = a();
//   double Theta = (ONE/(Mr2(deltax)*c*c) + ONE/(CP*Temp));
//   double phi = ZERO;   
//   for(int j=0; j<ns-1; j++){   
//     phi += spec[j].c*(specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/Rmix);   
//   }
//   double V = HALF*v.sqr();
//   double Omega = (Rmix - CP)*p/(rho*Rmix);
//   double psi = Theta*p/rho; 
//   double chi =  enthalpy - CP*p/(rho*Rmix) - phi - V;
//   double H = (enthalpy+V);
  
//   double AA = -H + TWO*V + k +chi +phi;
//   double BB = (psi-ONE);
//   double coefficient = ONE/(BB*AA+Omega);
 
//   Pinv.zero();
//   Pinv(0,0) =  BB*(AA-chi)+Omega;
//   Pinv(0,1) =  - BB*v.x;
//   Pinv(0,2) =  - BB*v.y;
//   Pinv(0,3) =   BB;

//   Pinv(1,0) = - BB*v.x*chi;
//   Pinv(1,1) =  BB*(AA-v.x*v.x)+Omega;
//   Pinv(1,2) = - BB*v.x*v.y;
//   Pinv(1,3) =  BB*v.x;

//   Pinv(2,0) = - BB*v.y*chi;
//   Pinv(2,1) = - BB*v.x*v.y;
//   Pinv(2,2) =  BB*(AA-v.y*v.y) +Omega; 
//   Pinv(2,3) =  BB*v.y;

//   Pinv(3,0) = - BB*H*chi;
//   Pinv(3,1) = - BB*H*v.x;
//   Pinv(3,2) = - BB*H*v.y;
//   Pinv(3,3) = BB*(AA+ H)+Omega;

//   Pinv(4,4) = BB*(AA-k)+Omega;
//   Pinv(5,5) = BB*AA+Omega;
  
//   if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
//       flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
//     Pinv(0,4) = -BB;
//     Pinv(1,4) = -BB*v.x;
//     Pinv(2,4) = -BB*v.y;
//     Pinv(3,4) = -BB*H;
//     //k-omega ...
//     Pinv(4,0) = -BB*k*chi;
//     Pinv(4,1) = -BB*k*v.x;
//     Pinv(4,2) = -BB*k*v.y;
//     Pinv(4,3) = BB*k;
    
//     Pinv(5,0) = -BB*omega*chi;
//     Pinv(5,1) = -BB*omega*v.x;
//     Pinv(5,2) = -BB*omega*v.y;
//     Pinv(5,3) = BB*omega;
//     Pinv(5,4) = -BB*omega;
    
    
//   }
//   //Multispecies (column entries)...
//   int NUM_VAR = NUM_CHEM2D_VAR_SANS_SPECIES;
//   for(int j=0; j<ns-1; j++){  
//     double enth_j = specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/Rmix;
//     Pinv(0,j+NUM_VAR) = -BB*enth_j;
//     Pinv(1,j+NUM_VAR) = -BB*v.x*enth_j;
//     Pinv(2,j+NUM_VAR) = -BB*v.y*enth_j;
//     Pinv(3,j+NUM_VAR) = -BB*H*enth_j;	
    
//     if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
// 	flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
      
//       Pinv(4,j+NUM_VAR) = -BB*k*enth_j;
//       Pinv(5,j+NUM_VAR) = -BB*omega*enth_j;
//     }
//   }
//   //Multispecies (row entries)
//   for(int i=0; i<ns-1; i++){ 
//     double enth_i = specdata[i].Enthalpy(Temp) + specdata[i].Heatofform() - CP*Temp*specdata[i].Rs()/Rmix;
//     Pinv(i+NUM_VAR,0) = -BB*(spec[i].c)*chi;
//     Pinv(i+NUM_VAR,1) = -BB*(spec[i].c)*v.x;
//     Pinv(i+NUM_VAR,2) = -BB*(spec[i].c)*v.y;
//     Pinv(i+NUM_VAR,3) = BB*(spec[i].c);
    
//     if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
// 	flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
//       Pinv(i+NUM_VAR,4) = -BB*(spec[i].c);
//     }
    
//     for(int j=0; j<ns-1; j++){ 
//       if(j==i){
// 	//diagonal
// 	Pinv(i+NUM_VAR,j+NUM_VAR) =BB*(AA - spec[i].c*enth_i)+ Omega;
//       }else{
// 	double enth_j = specdata[j].Enthalpy(Temp) + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/Rmix;
// 	Pinv(i+NUM_VAR,j+NUM_VAR) = -BB*spec[i].c*enth_j;
//       }
//     }//multispecies  
    
//   }//end of row entries 

 
//   Pinv = coefficient *Pinv; 
  
// }


/***********************************************************
 *************** Operator Overloading **********************
 ***********************************************************/

/********************************************************
 * Chem2D_pState -- Binary arithmetic operators.        *
 ********************************************************/
//----------------- Addition -----------------------------//
Chem2D_pState Chem2D_pState::operator +(const Chem2D_pState &W) const{ 
  
  if(ns == W.ns){ //check that species are equal   
    Chem2D_pState Temp(W.rho,W.v,W.p, W.k, W.omega);
    Temp.Copy(*this);
    Temp += W;
    return Temp;
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
}

//------------------ Subtraction ------------------------//
 Chem2D_pState Chem2D_pState::operator -(const Chem2D_pState &W) const{
  if(ns == W.ns){ //check that species are equal
    Chem2D_pState Temp(W.rho,W.v,W.p, W.k, W.omega);
    Temp.Copy(*this);
    Temp -= W;
    return Temp;
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
}

//---------------- Scalar Multiplication ------------------//
Chem2D_pState Chem2D_pState::operator *(const double &a) const{
  Chem2D_pState Temp(rho,v,p,k,omega);
  Temp.Copy(*this);
  Temp.rho = rho*a;  Temp.v = v*a; Temp.p = p*a;
  Temp.k = k*a; Temp.omega = omega*a;
  for( int i=0; i<ns; i++){
    Temp.spec[i] = spec[i]*a;
  } 
  Temp.tau = tau*a;
  Temp.qflux = qflux*a;
  Temp.lambda= lambda*a;
  Temp.theta = theta*a;

  return(Temp);
}

Chem2D_pState operator *(const double &a, const Chem2D_pState &W){
  Chem2D_pState Temp;
  //Temp.Copy(W);
  Temp.rho = W.rho*a;  Temp.v = W.v*a; Temp.p = W.p*a;
 Temp.k = W.k*a; Temp.omega = W.omega*a;
  for( int i=0; i<W.ns; i++){
    Temp.spec[i] = W.spec[i]*a;
  } 

  Temp.tau = W.tau*a;
  Temp.qflux = W.qflux*a;
  Temp.lambda= W.lambda*a;
  Temp.theta = W.theta*a;
  return(Temp);
}

//--------------- Scalar Division ------------------------//
Chem2D_pState Chem2D_pState::operator /(const double &a) const {
  Chem2D_pState Temp(rho,v,p,k,omega);
  Temp.Copy(*this);
  Temp.rho = rho/a; Temp.v = v/a; Temp.p = p/a; 
  Temp.k = k/a; Temp.omega = omega/a;
  for(int i=0; i<ns; i++){
     Temp.spec[i] = spec[i]/a; 
  } 

  Temp.tau = tau/a;
  Temp.qflux = qflux/a;
  Temp.lambda= lambda/a;
  Temp.theta = theta/a;

  return(Temp);
}

//----------------- Inner Product ------------------------//
double Chem2D_pState::operator *(const Chem2D_pState &W) const{
  double sum=0.0;
  if(ns == W.ns){ //check that species are equal
    for(int i=0; i<ns; i++){
      sum += spec[i]*W.spec[i];
    }  
    return (rho*W.rho + v*W.v + p*W.p + k*W.k+omega*W.omega+ sum);
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
}

//----------- solution state product operator ------------//
Chem2D_pState Chem2D_pState::operator ^( const Chem2D_pState &W) const {
  if(ns == W.ns){ //check that species are equal
    Chem2D_pState Temp(rho,v,p,k, omega);
    Temp.Copy(*this);
    Temp.rho = rho*W.rho;
    Temp.v.x = v.x*W.v.x;
    Temp.v.y = v.y*W.v.y;
    Temp.p = p*W.p;
    Temp.k = k*W.k;
    Temp.omega = omega*W.omega;
    for(int i=0; i<ns; i++){
      Temp.spec[i] = spec[i]*W.spec[i];
    }  
    return(Temp);
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
}

//----------------- Assignment ----------------------------//
Chem2D_pState& Chem2D_pState::operator =(const Chem2D_pState &W){
  //self assignment protection
  if( this != &W){   
    //copy assignment
    rho = W.rho;
    v = W.v; 
    p = W.p; 
    k = W.k;
    omega = W.omega;
    if ( ns == W.ns){
      for(int i=0; i<ns; i++){
	spec[i] = W.spec[i];
      }   
      tau = W.tau;
      qflux = W.qflux;
      lambda= W.lambda;
      theta = W.theta;
    }   
    else {
      cerr<<"\n Mismatch in number of species \n ";
    }  
  }
  return (*this);

}
/********************************************************
 * Chem2D_pState -- Shortcut arithmetic operators.     *
 ********************************************************/
Chem2D_pState& Chem2D_pState::operator +=(const Chem2D_pState &W){
  rho += W.rho;
  v += W.v; 
  p += W.p; 
  k += W.k;
  omega += W.omega;
  for( int i=0; i<ns; i++){
    spec[i] += W.spec[i];
  } 

  tau += W.tau;
  qflux += W.qflux;
  lambda += W.lambda;
  theta += W.theta;

  return (*this);
}

Chem2D_pState& Chem2D_pState::operator -=(const Chem2D_pState &W) {
  rho -= W.rho;
  v -= W.v;
  p -= W.p;
  k -= W.k;
  omega -= W.omega;
  for(int i=0; i<ns; i++){
    spec[i] -= W.spec[i];
  }

  tau -= W.tau;
  qflux -= W.qflux;
  lambda -= W.lambda;
  theta -= W.theta;

  return (*this); 
}

inline Chem2D_pState& Chem2D_pState::operator *=(const double &a) {
  rho *= a;
  v.x *= a;
  v.y *= a;
  p *= a;
  k *= a;
  omega *= a;
  for (int i = 0; i < ns; i++) spec[i] *= a;
  tau *= a;
  qflux *= a;
  lambda *= a;
  theta *= a;
  return *this;
}

inline Chem2D_pState& Chem2D_pState::operator /=(const double &a) {
  rho /= a;
  v.x /= a;
  v.y /= a;
  p /= a;
  k /= a;
  omega /= a;
  for (int i = 0; i < ns; i++) spec[i] /= a;
  tau /= a;
  qflux /= a;
  lambda /= a;
  theta /= a;
  return *this;
}

/********************************************************
 * Chem2D_pState -- Unary arithmetic operators.         *
 ********************************************************/
//  Chem2D_pState operator +(const Chem2D_pState &W) {  
//   return (Chem2D_pState(W.rho,W.v,W.p,W.spec));
// }

Chem2D_pState operator -(const Chem2D_pState &W) {
  Species *spt= new Species[W.ns];
  for(int i=0; i<W.ns; i++){
    spt[i] = -W.spec[i]; 
  }  

  Chem2D_pState Temp(-W.rho,-W.v,-W.p, -W.k, -W.omega ,spt);
  Temp.tau = -W.tau;
  Temp.qflux = -W.qflux;
  Temp.lambda = -W.lambda;
  Temp.theta= -W.theta;

  delete[] spt;
  return(Temp);
}

/********************************************************
 * Chem2D_pState -- Relational operators.               *
 ********************************************************/
int operator ==(const Chem2D_pState &W1, const Chem2D_pState &W2) {
  if(W1.ns == W2.ns){ //check that species are equal
    bool Temp;
    for(int i=0; i<W1.ns; i++){
      if( W1.spec[i] == W2.spec[i] ){
	Temp = true;
      } else {
	Temp = false;
	break;
      }  
      return (W1.rho == W2.rho && W1.v == W2.v && W1.p == W2.p&& W1.k == W2.k && W1.omega == W2.omega
	      && Temp == true && W1.tau == W2.tau &&
	      W1.qflux == W2.qflux&& W1.lambda == W2.lambda &&
	      W1.theta == W2.theta);
    }
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
}

int operator !=(const Chem2D_pState &W1, const Chem2D_pState &W2) {
   if(W1.ns == W2.ns){ //check that species are equal
    bool Temp = true;
    for(int i=0; i<W1.ns; i++){
      if( W1.spec[i] != W2.spec[i] ){
	Temp = false;
	break;
      } 
      return (W1.rho != W2.rho || W1.v != W2.v || W1.p != W2.p || W1.k != W2.k || W1.omega != W2.omega
	      || Temp != true || W1.tau != W2.tau ||
	      W1.qflux != W2.qflux|| W1.lambda != W2.lambda ||
	      W1.theta != W2.theta);
    }
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
}
/********************************************************
 * Chem2D_pState -- Input-output operators.            *
 ********************************************************/
ostream &operator << (ostream &out_file, const Chem2D_pState &W) {
  out_file.precision(10);
  out_file.setf(ios::scientific);
  out_file << " " << W.rho  << " " << W.v.x << " " << W.v.y << " " << W.p<< " " << W.k << " " << W.omega;
  for( int i=0; i<W.ns; i++){
    out_file<<" "<<W.spec[i];
  }
  out_file << " " << W.qflux << " " <<W.tau << " " << W.theta << " " << W.lambda;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

istream &operator >> (istream &in_file, Chem2D_pState &W) {
  in_file.setf(ios::skipws);
  in_file >> W.rho >> W.v.x >> W.v.y >> W.p >> W.k >> W.omega;
  //W.set_initial_values();
  for( int i=0; i<W.ns; i++){
    in_file>>W.spec[i];
  }
  in_file >>W.qflux >>W.tau >>W.theta >>W.lambda;
  in_file.unsetf(ios::skipws);
   return (in_file);
}
/***************************************************************
 ***************************************************************
 *************** AXISYMMETRIC SOURCE Jacobians *****************
 ***************************************************************
 ***************************************************************/

/****************************************************************
 * Axisymmetric Source Term Jacboian (Inviscid)                 * 
 ****************************************************************/
void Chem2D_pState::dSa_idU(DenseMatrix &dSa_IdU,
			    const Vector2D &X,
			    const int Axisymmetric) const {

  double enthalpy = h();
  double CP = Cp();
  double CV = Cv();
  double Gamma = g();
  double RTOT = Rtot();
  double phi = ZERO;
  double Temp = p/(rho*RTOT);
  for(int j=0; j<ns-1; j++){   
    phi += spec[j].c*(specdata[j].Enthalpy(Temp) 
 		      + specdata[j].Heatofform() - CP*Temp*specdata[j].Rs()/RTOT);   
  }
 
  if(Axisymmetric == 1){
    dSa_IdU(0,2) -= ONE;
    dSa_IdU(1,0) += v.x*v.y;
    dSa_IdU(1,1) -= v.y;
    dSa_IdU(1,2) -= v.x;
    dSa_IdU(2,0) += v.y*v.y;
    dSa_IdU(2,2) -= TWO*v.y;
    
    if((flow_type == FLOWTYPE_INVISCID) || (flow_type == FLOWTYPE_LAMINAR)) {
      double ee = e();
      double dedp = rho*rho*diedip();
      dSa_IdU(3,0) += v.y*(rho*ee +HALF*(v.x*v.x +v.y*v.y) + 
			   (-HALF*(v.x*v.x +v.y*v.y) + TWO*rho*ee)/(dedp));
      dSa_IdU(3,1) += ((v.x*v.y)/dedp);
      dSa_IdU(3,2) -= (rho*ee +HALF*(v.x*v.x +v.y*v.y) + p/rho - (v.y*v.y)/(dedp));
      dSa_IdU(3,3) -= v.y*(ONE + ONE/dedp);
      
    } //end of laminar or inviscid cases
    //multispecies terms 
    int NUM_VAR = NUM_CHEM2D_VAR_SANS_SPECIES;
    for(int i=0; i<(ns-1);i++){
      
      dSa_IdU(3,i+NUM_VAR) += v.y*(specdata[i].Enthalpy(Temp) 
				   + specdata[i].Heatofform() - CP*Temp*specdata[i].Rs()/RTOT)/(CP/RTOT - ONE); 
      //combine turbulence and laminar cases together
      dSa_IdU(NUM_VAR+i,0) += v.y*spec[i].c;
      dSa_IdU(NUM_VAR+i,2) -= spec[i].c;
      dSa_IdU(NUM_VAR+i,NUM_VAR+i) -= v.y;
    }
    
    if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
	flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
      
      dSa_IdU(3,0) -= -(v.x*v.x +v.y*v.y)*v.y+ Gamma*v.y*(v.x*v.x+v.y*v.y+2.0*RTOT*Temp)/2.0-k*v.y -v.y*h();
      dSa_IdU(3,1) -= v.x*v.y- Gamma *v.y*v.x;
      dSa_IdU(3,2) -= THREE/TWO*v.y*v.y +HALF*v.x*v.x +h() +k -Gamma*v.y*v.y;
      dSa_IdU(3,3) -= v.y*Gamma;
      dSa_IdU(3,4) -= -dSa_IdU(3,3)+v.y;
      
      dSa_IdU(4,0) = k*v.y;
      dSa_IdU(5,0) = omega*v.y;
      dSa_IdU(4,2) -= k;
      dSa_IdU(5,2) -= omega;
      dSa_IdU(4,4) -= v.y;
      dSa_IdU(5,5) -= v.y;
      
    }//end of turbulent case
    
    dSa_IdU = dSa_IdU/X.y;

    //end of axisymmetric case 1  -- y radial direction 
  } else if(Axisymmetric ==2){ 
    
    dSa_IdU(0,1) -= ONE;
    dSa_IdU(1,0) += v.x*v.x;
    dSa_IdU(1,1) -= TWO*v.x;
    dSa_IdU(2,0) += v.x*v.y;
    dSa_IdU(2,1) -= v.y;
    dSa_IdU(2,2) -= v.x;

    if((flow_type == FLOWTYPE_INVISCID) || (flow_type ==FLOWTYPE_LAMINAR)) {
    
      dSa_IdU(3,0) -= RTOT/(CP-RTOT)*(v.x*(phi + v.sqr()) + v.x*CP*( p/rho - enthalpy - HALF*v.sqr())/RTOT);
      dSa_IdU(3,1) -= (enthalpy + CP/(CP-RTOT)*HALF*v.sqr() - RTOT/(CP-RTOT)*HALF*(THREE*v.x*v.x +v.y*v.y));
      dSa_IdU(3,2) += RTOT/(CP-RTOT)*v.x*v.y;
      dSa_IdU(3,3) -= CP/(CP-RTOT)*v.x;
      
    }//end of laminar case
    
    //Multispecies terms
    int NUM_VAR = NUM_CHEM2D_VAR_SANS_SPECIES;
    for(int i=0; i<(ns-1);i++){
      dSa_IdU(3,i+NUM_VAR) += v.x*(specdata[i].Enthalpy(Temp) 
				   + specdata[i].Heatofform() - CP*Temp*specdata[i].Rs()/RTOT)/(CP/RTOT - ONE); 
      dSa_IdU(NUM_VAR+i,0) += v.x*spec[i].c;
      dSa_IdU(NUM_VAR+i,1) -= spec[i].c;
      dSa_IdU(NUM_VAR+i,NUM_VAR+i) -= v.x;
    }
    
    if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
	flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
       
      dSa_IdU(3,0) -= -(v.x*v.x +v.y*v.y)*v.x+ Gamma*v.x*(v.x*v.x+v.y*v.y+2.0*RTOT*Temp)/2.0-k*v.x -v.x*h();
      dSa_IdU(3,1) -= THREE/TWO*v.x*v.x +HALF*v.y*v.y +h() +k -Gamma*v.x*v.x;
      dSa_IdU(3,2) -= v.x*v.y- Gamma*v.y*v.x;
      
      dSa_IdU(3,3) -= Gamma*v.x;
      dSa_IdU(3,4) -= -dSa_IdU(3,3)+v.x;
    
      dSa_IdU(4,0) = k*v.x;
      dSa_IdU(4,1) -= k;
      dSa_IdU(4,4) -= v.x;
      dSa_IdU(5,0) = omega*v.x;
      dSa_IdU(5,1) = omega;
      dSa_IdU(5,5) -= v.x;  
    }//end of turbulent case

    dSa_IdU = dSa_IdU/X.x;

    //    cout<<"\n Axisymmetric-x dSa_IdU";
  }//end of axisymmetric case   -- x radial 
}
 
/****************************************************************
 * Axisymmetric Source Term Jacboian (Viscous)                  * 
 ****************************************************************/
void Chem2D_pState::dSa_vdU(DenseMatrix &dSa_VdU,
			    DenseMatrix &dWdQ,
			    const Chem2D_pState &dWdx,
			    const Chem2D_pState &dWdy,
			    const Vector2D &X, 
			    const int Axisymmetric, 
			    const double d_dWdx_dW,
			    const double d_dWdy_dW) const {

  DenseMatrix dSa_VdW(NUM_VAR_CHEM2D-1,NUM_VAR_CHEM2D-1); 
  dSa_VdW.zero();
  
  double enthalpy = h();
  double CP = Cp();
  double RTOT = Rtot();
  double phi = ZERO;
  double Temp = p/(rho*RTOT);
  
  if((flow_type ==FLOWTYPE_INVISCID)||(flow_type ==FLOWTYPE_LAMINAR)){
    if(Axisymmetric == 1){  
      double t1,t2,t4,t6, t7, t8, t9, t11;
      double t14,t15,t23,t25,t26,t27,t28,t32;
      double t54,t56,t62,t71;

      double r =X.y; 
      t1 = mu();
      t2 = d_dWdy_dW ;
      t4 = d_dWdx_dW;
      t6 = dmudT();
      t7 = dWdy.v.x;
      t8 = dWdx.v.y; 
      t9 = t7+t8;
      t11 = d_dWdx_dW ;
      t14 = d_dWdy_dW;
      t15 = 1/r;
      t23 = dWdx.v.x; 
      t25 = dWdy.v.y; 
      t26 = t25/3.0;
      t27 = v.y*t15;
      t28 = t27/3.0;
      t32 = t23/3.0;
      t54 = v.x*t1;
      t56 = v.y*t1;
      t62 = 2.0/3.0*t25-t32-t28;
      t71 =  d_dWdy_dW ;

      dSa_VdW(1,1) = t1*t2;
      dSa_VdW(1,2) = t1*t4;
      dSa_VdW(1,3) = t6*t9;

      dSa_VdW(2,1) = 2.0*t1*t11;
      dSa_VdW(2,2) = 2.0/3.0*t1*(-t14-t15)-2.0*t1*(2.0/3.0*t15-t14/3.0);
      dSa_VdW(2,3) = 2.0*t6*(2.0/3.0*t23-t26-t28)-2.0*t6*(2.0/3.0*t27-t32-t26);

      dSa_VdW(3,0) = thermal_diffusion().y;
      dSa_VdW(3,1) = t1*t9+t54*t2-2.0/3.0*t56*t11;
      dSa_VdW(3,2) = t54*t4+2.0*t1*t62+2.0*t56*(2.0/3.0*t14-t15/3.0);

      double Sum_dh = 0.0;
      for(int Num = 0; Num<ns; Num++)
	{
	  // sum of (dhidT *(Dmi)*gradc
	  Sum_dh +=  specdata[Num].Enthalpy_prime(Temp)
	    *spec[Num].diffusion_coef*dWdy.spec[Num].c;
	  
	}
      dSa_VdW(3,3)=  kappa()*t71+rho*Sum_dh +v.x*t6*t9+2.0*v.y*t6*t62;
      int NUM_VAR =  NUM_CHEM2D_VAR_SANS_SPECIES;

      for(int Num = 0; Num<(ns-1); Num++)
	{ 
	  
	  dSa_VdW(3,NUM_VAR+Num)= rho*spec[Num].diffusion_coef*(specdata[Num].Enthalpy(Temp)+specdata[Num].Heatofform())*d_dWdy_dW;
	  dSa_VdW(NUM_VAR+Num,0)= spec[Num].diffusion_coef*dWdy.spec[Num].c;
	  dSa_VdW(NUM_VAR+Num,NUM_VAR+Num) = rho*spec[Num].diffusion_coef*d_dWdy_dW;
	}
      //then devided by the radial distance
      
      dSa_VdU = dSa_VdW*dWdQ;
      dSa_VdU = dSa_VdU/r;
      
    }//end of axisymmetric case -- y

    if(Axisymmetric == 2){  
      double t1,t2,t4,t6, t9, t12;
      double t13,t15,t19,t22,t24,t26,t27,t28;
      double t49,t52,t60;
      
      double r =X.x; 
    
      t1 = mu();
      t2 = d_dWdx_dW;
      t4 = 1/r;
      t6 = 2.0/3.0*t2-t4/3.0;
      t9 = d_dWdy_dW;
      t12 = dmudT();
      t13 = dWdx.v.x;
      t15 = dWdy.v.y;
      t19 = 2.0/3.0*t13-t15/3.0-v.x*t4/3.0;
      t22 = d_dWdy_dW;
      t24 = d_dWdx_dW;
      t26 = dWdy.v.x;
      t27 = dWdx.v.y;
      t28 = t26+t27;

      t49 = v.x*t1;
      t52 = v.y*t1;

      t60 = d_dWdx_dW;

      dSa_VdW(1,1) = 2.0*t1*t6;
      dSa_VdW(1,2) = -2.0/3.0*t1*t9;
      dSa_VdW(1,3) = 2.0*t12*t19;
   
      dSa_VdW(2,1) = t1*t22;
      dSa_VdW(2,2) = t1*t24;
      dSa_VdW(2,3) = t12*t28;
    
      dSa_VdW(3,0) = thermal_diffusion().x;
      dSa_VdW(3,1) = 2.0*t1*t19+2.0*t49*t6+t52*t22;
      dSa_VdW(3,2) = -2.0/3.0*t49*t9+t1*t28+t52*t24;
      double Sum_dh = 0.0;
      for(int Num = 0; Num<ns; Num++)
	{
	  // sum of (dhidT *(Dmi)*gradc
	  Sum_dh +=  specdata[Num].Enthalpy_prime(Temp)
	    *spec[Num].diffusion_coef*dWdx.spec[Num].c;
	  
	}
      dSa_VdW(3,3) = kappa()*t60+rho*Sum_dh+2.0*v.x*t12*t19+v.y*t12*t28;
	  
      int NUM_VAR = NUM_CHEM2D_VAR_SANS_SPECIES;
      for(int Num = 0; Num<(ns-1); Num++)
	{ 
	  //rho*Di*hi*d_dWdx_dW
	  dSa_VdW(3,NUM_VAR+Num)= rho*spec[Num].diffusion_coef*(specdata[Num].Enthalpy(Temp)+specdata[Num].Heatofform())*d_dWdx_dW;
	  dSa_VdW(NUM_VAR+Num,0)= spec[Num].diffusion_coef*dWdx.spec[Num].c;
	  dSa_VdW(NUM_VAR+Num,NUM_VAR+Num) = rho*spec[Num].diffusion_coef*d_dWdx_dW;
	}
    
      dSa_VdU = dSa_VdW*dWdQ;
      dSa_VdU = dSa_VdU/r;
      
    }//x-radius 
    
  }//end of laminar case
  
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
      flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
    
    double Temp =T();
    double r;

    if(Axisymmetric == 1){
      double t1,t2,t3,t4,t5,t7,t8,t10;
      double t13,t14,t18,t19,t21,t23,t24,t27,t29;
      double t31,t34,t38,t39,t43,t44,t47,t51,t57; 
      double  t62, t63 ,t91, t92, t93, t95, t116;
      double t146, t155, t156, t159, t160, t161; 
      double t166, t170, t200, t201, t204, t210; 
    
      r = X.y;
      t1 = 1/max(omega,TOLER);
      t2 = k*t1;
      t3 =  dWdy.v.x; 
      t4 =  dWdx.v.y;
      t5 = t3+t4;
      t7 = mu();
      t8 =  d_dWdy_dW;
      t10 = rho*k;
      t13 = t7*t8+t10*t1*t8;
      t14 = d_dWdx_dW;
      t18 = t7*t14+t10*t1*t14;
      t19 = dmudT();
      t21 = rho*t1;
      t23 = max(omega,TOLER)*max(omega,TOLER);
      t24 = 1/t23;
      t27 = dWdy.v.y;
      t29 = dWdx.v.x; 
    
      t31 = 1/r;
      t34 = 2.0/3.0*t27-t29/3.0-v.x*t31/3.0;
      t38 = 2.0*t2*t34-2.0/3.0*k;
      t39 =  d_dWdx_dW;
      t43 = -t7*t39-t10*t1*t39;
      t44 =  d_dWdy_dW;
      t47 = 2.0/3.0*t44-t31/3.0;
      t51 = t7*t47+t10*t1*t47;
      t57 = 2.0*t21*t34-2.0/3.0*rho;
      t62 = Cp()/Pr_turb();
      // t63 = dTdy;
      t63 = (ONE/(rho*Rtot()))*(dWdy.p - (p/rho)*dWdy.rho);
      
      t91 = dWdy.k;
      t92 = t1*t91;
      t93 = sigma_star*k*t92;
      t95 = t1*t5;
      t116 = d_dWdy_dW;
      t146 = t19*t91;
      t155 = sigma_star*rho;
      t156 = t155*t92;
      t159 = d_dWdy_dW;
      t160 = (t7+t155*t2)*t159;
      t161 = v.x*rho;
      t166 = k*t24;
      t170 = t155*t166*t91;
    
      t200 = dWdy.omega;
      t201 = t1*t200;
      t204 = sigma*rho;
      t210 = d_dWdy_dW;
      
      dSa_VdW(1,0) = t2*t5;
      dSa_VdW(1,1) = t13;
      dSa_VdW(1,2) = t18;
      dSa_VdW(1,3) = t19*t5;
      dSa_VdW(1,4) = t21*t5;
      dSa_VdW(1,5) = -t10*t24*t5;
 
      dSa_VdW(2,0) = t38;
      dSa_VdW(2,1) = 2.0/3.0*t43;
      dSa_VdW(2,2) = 2.0*t51;
      dSa_VdW(2,3) = 2.0*t19*t34;
      dSa_VdW(2,4) = t57;
      dSa_VdW(2,5) = -2.0*t10*t24*t34;
      
      double Sum_q = 0.0;
      double Sum_dq = 0.0;
      // 6 represents rho, vr,vz, p, k, omega in 2D axisymmetric turbulent flows
      
      for(int Num = 0; Num<(ns); Num++)
	{
	  //for each species	// h*(Dm+Dt)*gradc  
	  Sum_q +=  (specdata[Num].Enthalpy(Temp)+specdata[Num].Heatofform())*(spec[Num].diffusion_coef+Dm_turb())*dWdy.spec[Num].c;
	  //dhdT *(Dm+Dmt)*gradc
	  Sum_dq +=  specdata[Num].Enthalpy_prime(Temp)*(spec[Num].diffusion_coef+Dm_turb())*dWdy.spec[Num].c;
	  
	}
      dSa_VdW(3,0) = t62*t2*t63+Sum_q+t93+v.x*k*t95+v.y*t38;
      dSa_VdW(3,1) = t7*t5+t10*t95+v.x*t13+2.0/3.0*v.y*t43;
      dSa_VdW(3,2) = v.x*t18+2.0*t7*t34+2.0*t10*t1*t34-2.0/3.0*t10+2.0*v.y*t51;
      dSa_VdW(3,3) = (kappa()+t62*t10*t1)*t116+ rho*Sum_dq+t146+v.x*t19*t5+2.0*v.y*t19*t34;
      dSa_VdW(3,4) = t62*t21*t63+t156+t160+t161*t95+v.y*t57;
      dSa_VdW(3,5) = -t62*rho*t166*t63-t170-t161*t166*t5-2.0*v.y*rho*t166*t34;

      int NUM_VAR = NUM_CHEM2D_VAR_SANS_SPECIES;
     for(int Num = 0; Num<(ns-1); Num++)
       {//h*rho*(D+Dt)*d_dWdy_dW
	 
	  dSa_VdW(3, NUM_VAR+Num) = (specdata[Num].Enthalpy(Temp)+specdata[Num].Heatofform())*rho*(spec[Num].diffusion_coef+Dm_turb())*d_dWdy_dW;
	  //(D+Dt)*gradc
	  dSa_VdW(NUM_VAR+Num, 0) =(spec[Num].diffusion_coef+Dm_turb())*dWdy.spec[Num].c;
	  //rho*(D+Dt)*dcydc
	  dSa_VdW(NUM_VAR+Num,NUM_VAR+Num) = rho*(spec[Num].diffusion_coef+Dm_turb())*d_dWdy_dW;
	}// end the center ones
      
     dSa_VdW(4,0) = t93;
     dSa_VdW(4,3) = t146;
     dSa_VdW(4,4) = t156+t160;
     dSa_VdW(4,5) = -t170;
     dSa_VdW(5,0) = sigma_star*k*t201;
     dSa_VdW(5,3) = t19*t200;
     dSa_VdW(5,4) = t204*t201;
     dSa_VdW(5,5) = -t204*t166*t200+(t7+t204*t2)*t210;
     
     dSa_VdU = dSa_VdW*dWdQ;
     dSa_VdU = dSa_VdU/r;
    }  
    if(Axisymmetric == 2){
   
      double t1,t2,t3,t5,t7,t10;
      double t14,t15,t16,t19,t21,t24,t25,t29,t30;
      double t33,t37,t38,t39,t43,t44,t45,t47,t51; 
      double  t52, t56 ,t62, t63 ,t91, t92, t93, t96 ;
      double t146, t155, t156, t159, t160, t162; 
      double t166, t170, t200, t201, t204, t210; 
      r = X.x;
      t1 = 1/max(omega,TOLER);
      t2 = k*t1;
      t3 = dWdx.v.x;
      t5 = dWdy.v.y;
      t7 = 1/r;
      t10 = 2.0/3.0*t3-t5/3.0-v.x*t7/3.0;
      t14 = 2.0*t2*t10-2.0/3.0*k;
      t15 = mu();
      t16 = d_dWdx_dW;
      t19 = 2.0/3.0*t16-t7/3.0;
      t21 = rho*k;
      t24 = t15*t19+t21*t1*t19;
      t25 = d_dWdy_dW;
      t29 = -t15*t25-t21*t1*t25;
      t30 = dmudT();
      t33 = rho*t1;
      t37 = 2.0*t33*t10-2.0/3.0*rho;
      t38 = max(omega,TOLER)*max(omega,TOLER);
      t39 = 1/t38;
      t43 = dWdy.v.x;
      t44 = dWdx.v.y;
      t45 = t43+t44;
      t47 = d_dWdy_dW;
      t51 = t15*t47+t21*t1*t47;
      t52 = d_dWdx_dW;
      t56 = t15*t52+t21*t1*t52;
      t62 = CP/Pr_turb();
      t63 = (ONE/(rho*Rtot()))*(dWdx.p - (p/rho)*dWdx.rho);//dTdx;

      t91 = dWdx.k;
      t92 = t1*t91;
      t93 = sigma*k*t92;
      t96 = t1*t45;
   
      t146 = t30*t91;
      t155 = sigma*rho;
      t156 = t155*t92;
      t159 = d_dWdx_dW;
      t160 = (t15+t155*t2)*t159;
      t162 = v.y*rho;
      t166 = k*t39;
      t170 = t155*t166*t91;

      t200 = dWdx.omega;
      t201 = t1*t200;
      t204 = sigma_star*rho;
      t210 = d_dWdx_dW;
   
      dSa_VdW(1,0) = t14;
      dSa_VdW(1,1) = 2.0*t24;
      dSa_VdW(1,2) = 2.0/3.0*t29;
      dSa_VdW(1,3) = 2.0*t30*t10;
      dSa_VdW(1,4) = t37;
      dSa_VdW(1,5) = -2.0*t21*t39*t10;
  
      dSa_VdW(2,0) = t2*t45;
      dSa_VdW(2,1) = t51;
      dSa_VdW(2,2) = t56;
      dSa_VdW(2,3) = t30*t45;
      dSa_VdW(2,4) = t33*t45;
      dSa_VdW(2,5) = -t21*t39*t45;
      double Sum_q = 0.0;
      double Sum_dq = 0.0;
      for(int Num = 0; Num<(ns); Num++)
	{
	  //for each species	// h*(Dm+Dt)*gradc  
	  Sum_q +=  (specdata[Num].Enthalpy(Temp)+specdata[Num].Heatofform())*(spec[Num].diffusion_coef+Dm_turb())*dWdx.spec[Num].c;
	  //dhdT *(Dm+Dmt)*gradc
	  Sum_dq +=  specdata[Num].Enthalpy_prime(Temp)*(spec[Num].diffusion_coef+Dm_turb())*dWdx.spec[Num].c;
	  
	}

    
      dSa_VdW(3,0) = t62*t2*t63 + Sum_dq +t93+v.x*t14+v.y*k*t96;
      dSa_VdW(3,1) = 2.0*t15*t10+2.0*t21*t1*t10-2.0/3.0*t21+2.0*v.x*t24+v.y*t51;
      dSa_VdW(3,2) = 2.0/3.0*v.x*t29+t15*t45+t21*t96+v.y*t56;
      dSa_VdW(3,3) =(kappa()+t62*t21*t1)*t52+rho*Sum_dq+t146+2.0*v.x*t30*t10+v.y*t30*t45;
      dSa_VdW(3,4) = t62*t33*t63+t156+t160+v.x*t37+t162*t96;
      dSa_VdW(3,5) = -t62*rho*t166*t63-t170-2.0*v.x*rho*t166*t10-t162*t166*t45;
      
      int NUM_VAR =  NUM_CHEM2D_VAR_SANS_SPECIES;
 
      for(int Num = 0; Num<(ns-1); Num++)
	{//h*rho*(D+Dt)*d_dWdy_dW
	  dSa_VdW(3, NUM_VAR+Num) = (specdata[Num].Enthalpy(Temp)+specdata[Num].Heatofform())*rho*(spec[Num].diffusion_coef+Dm_turb())*d_dWdx_dW;
	  //(D+Dt)*gradc
	  dSa_VdW(NUM_VAR+Num, 0) =(spec[Num].diffusion_coef+Dm_turb())*dWdx.spec[Num].c;
	  //rho*(D+Dt)*dcydc
	  dSa_VdW(NUM_VAR+Num,NUM_VAR+Num) = rho*(spec[Num].diffusion_coef+Dm_turb())*d_dWdx_dW;
	}// end the center ones
      dSa_VdW(4,0) = t93;
      dSa_VdW(4,3) = t146;
      dSa_VdW(4,4) = t156+t160;
      dSa_VdW(4,5) = -t170;
      
      dSa_VdW(5,0) = sigma_star*k*t201;
      dSa_VdW(5,3) = t30*t200;
      dSa_VdW(5,4) = t204*t201;
      dSa_VdW(5,5) = -t204*t166*t200+(t15+t204*t2)*t210;
      
      dSa_VdU = dSa_VdW*dWdQ;
      dSa_VdU = dSa_VdU/r;
  
    }//end of axisymmetric x

  }//end of turbulent case 
  
}

/***************************************************************
 ***************************************************************
 *************** SOURCE TERMS **********************************
 ***************************************************************
 ***************************************************************/

/***************************************************************
 * Axisymmetric flow source terms (Inviscid)                   *
 ***************************************************************/
Chem2D_cState Chem2D_pState::Sa_inviscid(const Vector2D &X,
                                         const int Axisymmetric) const{

  Chem2D_cState Temp; Temp.Vacuum();

  if (Axisymmetric == 2) {
     Temp.rho = -rho*v.x/X.x;
     Temp.rhov.x = -rho*v.x*v.x/X.x; 
     Temp.rhov.y = -rho*v.x*v.y/X.x;
     Temp.E =  -v.x*H()/X.x;
     //species contributions
     for(int i=0; i<ns;i++){
       Temp.rhospec[i].c = -rho*v.x*spec[i].c/X.x;
     }
     if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
         flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
        Temp.rhok =  -v.x*rho*k/X.x;
        Temp.rhoomega = -v.x*rho*omega/X.x;
     } 
//      cout<<"\n Axisymmetric-x Sa_inviscid";
  } else if (Axisymmetric == 1) {
     Temp.rho = -rho*v.y/X.y;
     Temp.rhov.x = -rho*v.x*v.y/X.y; 
     Temp.rhov.y = -rho*v.y*v.y/X.y;
     Temp.E =  -v.y*H()/X.y;
     //species contributions
     for(int i=0; i<ns;i++){
       Temp.rhospec[i].c = -rho*v.y*spec[i].c/X.y;
     }
     if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
         flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
        Temp.rhok =  -v.y*rho*k/X.y;
        Temp.rhoomega = -v.y*rho*omega/X.y;
     }
  }

  return (Temp);

}

/***************************************************************
 * Axisymmetric flow source terms (Viscous)                    *  
 ***************************************************************/
Chem2D_cState Chem2D_pState::Sa_viscous(const Chem2D_pState &dWdx,
					const Chem2D_pState &dWdy,
                                        const Vector2D &X,
                                        const int Axisymmetric){

  double div_v, radius;
  double Mu, Kappa, Temperature, Rmix, Cpmix;
  double mu_t, kappa_t, Dm_t;
  double rhohsDs;
  Vector2D grad_T;
  Tensor2D strain_rate;
  Chem2D_cState Temp; Temp.Vacuum();

  //Transport and thermodynamic properties
  Mu = mu();
  Kappa = kappa();
  Temperature = T();
  Rmix = Rtot();
  Cpmix = Cp();

  //Turbulence model transport properties
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
      flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
      mu_t = eddy_viscosity();
      kappa_t =  mu_t*Cpmix/Pr_turb();
      Dm_t = Dm_turb();
  } /* endif */
  //Strain rate (+dilatation)
  div_v = dWdx.v.x + dWdy.v.y;
  if (Axisymmetric == 2) {
     radius = X.x;
     div_v += v.x/radius;
  } else if (Axisymmetric == 1) {
     radius = X.y;
     div_v += v.y/radius;
  } /* endif */
  strain_rate.xx = dWdx.v.x-div_v/THREE;
  strain_rate.xy = HALF*(dWdx.v.y + dWdy.v.x);
  strain_rate.yy = dWdy.v.y-div_v/THREE;
  if (Axisymmetric == 0) {
     strain_rate.zz = -(strain_rate.xx + strain_rate.yy); 
  } else if (Axisymmetric == 2) {
     strain_rate.zz = v.x/radius-div_v/THREE;
  } else if (Axisymmetric == 1) {
     strain_rate.zz = v.y/radius-div_v/THREE;
  } /* endif */

  //Temperature gradient
  //dT/dx = 1/rho*R *( dP/dx - P/rho * drho/dx)
  grad_T.x = (ONE/(rho*Rmix))*(dWdx.p - (p/rho)*dWdx.rho);
  grad_T.y = (ONE/(rho*Rmix))*(dWdy.p - (p/rho)*dWdy.rho);
 //Molecular (laminar) stress tensor
  tau = TWO*Mu*strain_rate;

  //Turbulent (Reynolds) Stresses
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
      flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
     lambda = (TWO*mu_t)*strain_rate;
     lambda.xx -= (TWO/THREE)*rho*k; 
     lambda.yy -= (TWO/THREE)*rho*k; 
     lambda.zz -= (TWO/THREE)*rho*k; 
  } /* endif */

  //Molecular (laminar) heat flux
  //Thermal conduction, q = - kappa * grad(T)
  qflux = - Kappa*grad_T;
  //Thermal diffusion, q -= rho * sum ( hs * Ds *gradcs)
  for(int i=0; i<ns; i++){ 
    rhohsDs = rho*(Mu/(rho*Schmidt[i]))*(specdata[i].Enthalpy(Temperature) + specdata[i].Heatofform());
    qflux.x -= rhohsDs*dWdx.spec[i].c;
    qflux.y -= rhohsDs*dWdy.spec[i].c;
  }
 
  //Turbulent heat flux
  //Thermal conduction, q = - kappa * grad(T)
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
      flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {

    theta = - kappa_t*grad_T;
    //Thermal Diffusion, q -= rho * sum ( hs * Ds *gradcs)  
    for (int i=0; i<ns; i++) {
      rhohsDs = rho*Dm_t*(specdata[i].Enthalpy(Temperature)+specdata[i].Heatofform());
      theta.x -= rhohsDs*dWdx.spec[i].c;
      theta.y -= rhohsDs*dWdy.spec[i].c;


    }
  } /* endif */

  //Determine axisymmetric source terms
  if (Axisymmetric == 2) {
    Temp.rhov.x = (tau.xx - tau.zz)/X.x;
    Temp.rhov.y = tau.xy/X.x;
    Temp.E = (- qflux.x + v.x*tau.xx + v.y*tau.xy)/X.x;
    for(int i=0; i<ns;i++){
      Temp.rhospec[i].c = rho*(Mu/(rho*Schmidt[i]))*dWdx.spec[i].c/X.x;
    }
    if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
	flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
      Temp.rhov.x += (lambda.xx - lambda.zz)/X.x;
      Temp.rhov.y += lambda.xy/X.x;
      Temp.E += (- theta.x + v.x*lambda.xx + v.y*lambda.xy)/X.x 
	+(Mu + mu_t*sigma_star)*dWdx.k/X.x;
      Temp.rhok = (Mu + mu_t*sigma_star)*dWdx.k/X.x;
      Temp.rhoomega = (Mu + mu_t*sigma)*dWdx.omega/X.x;
      for(int i=0; i<ns;i++){
	Temp.rhospec[i].c = rho*Dm_t*dWdx.spec[i].c/X.x;
      }
    }
    
  } else if (Axisymmetric == 1) {
    Temp.rhov.x = tau.xy/X.y;
    Temp.rhov.y = (tau.xx - tau.zz)/X.y;
    Temp.E = (- qflux.y + v.x*tau.xy + v.y*tau.yy)/X.y;
    for(int i=0; i<ns;i++){
      Temp.rhospec[i].c = rho*(Mu/(rho*Schmidt[i]))*dWdy.spec[i].c/X.y;
    }
    if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
	flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
      Temp.rhov.x += lambda.xy/X.y;
      Temp.rhov.y += (lambda.xx - lambda.zz)/X.y;
      Temp.E += (- theta.y + v.x*lambda.xy + v.y*lambda.yy)/X.y 
	+(Mu + mu_t*sigma_star)*dWdy.k/X.y ;
      Temp.rhok = (Mu + mu_t*sigma_star)*dWdy.k/X.y;
      Temp.rhoomega = (Mu + mu_t*sigma)*dWdy.omega/X.y;
      for(int i=0; i<ns;i++){
	Temp.rhospec[i].c = rho*Dm_t*dWdy.spec[i].c/X.y;
      }
    }
    
  } /* endif */
  
 
  return (Temp);
  
}

/***************************************************************
 * Turbulence model source terms                               *  
 ***************************************************************/
Chem2D_cState Chem2D_pState::S_turbulence_model(const Chem2D_pState &dWdx,
					        const Chem2D_pState &dWdy,
                                                const Vector2D &X,
                                                const int Axisymmetric){

  double div_v, radius;
  double mu_t, production;
  Tensor2D strain_rate;
  Chem2D_cState Temp; Temp.Vacuum();

  //Turbulence model eddy viscosity
  mu_t = eddy_viscosity();

  //Strain rate (+dilatation)
  div_v = dWdx.v.x + dWdy.v.y;
  if (Axisymmetric == 2) {
     radius = X.x;
     div_v += v.x/radius;
  } else if (Axisymmetric == 1) {
     radius = X.y;
     div_v += v.y/radius;
  } /* endif */
  strain_rate.xx = dWdx.v.x-div_v/THREE;
  strain_rate.xy = HALF*(dWdx.v.y + dWdy.v.x);
  strain_rate.yy = dWdy.v.y-div_v/THREE;
  if (Axisymmetric == 0) {
     strain_rate.zz = -(strain_rate.xx + strain_rate.yy); 
  } else if (Axisymmetric == 2) {
     strain_rate.zz = v.x/radius-div_v/THREE;
  } else if (Axisymmetric == 1) {
     strain_rate.zz = v.y/radius-div_v/THREE;
  } /* endif */

  //Turbulent (Reynolds) Stresses
  lambda = (TWO*mu_t)*strain_rate;
  lambda.xx -= (TWO/THREE)*rho*k; 
  lambda.yy -= (TWO/THREE)*rho*k; 
  lambda.zz -= (TWO/THREE)*rho*k; 

  //Production term
  production = lambda.xx*dWdx.v.x + 
               lambda.xy*(dWdy.v.x + dWdx.v.y) + 
               lambda.yy*dWdy.v.y;
  if (Axisymmetric == 2) {
     production += lambda.zz*v.x/radius;
  } else if (Axisymmetric == 1) {
     production += lambda.zz*v.y/radius;
  } /* endif */

  //Determine axisymmetric source terms
  Temp.rhok = production - 
              f_beta_star*beta_star*rho*k*omega;
  Temp.rhoomega = alpha*(omega/max(k, TOLER))*production -
                  f_beta*beta*rho*omega*omega;

  return (Temp);

//   double radius;
//   double  production;
//   Chem2D_cState Temp; Temp.Vacuum();

// //Xinfeng Note: the lambda = Reynolds_Stress(... ... ) may be computed/assigned outside
// //To see which way is more general ...
// // maybe modifying the viscous calculations some sort...
//   lambda = Reynolds_Stress(dWdx,dWdy,Axisymmetric, X);
					       
					      
//   //Production term
//   production = lambda.xx*dWdx.v.x + 
//                lambda.xy*(dWdy.v.x + dWdx.v.y) + 
//                lambda.yy*dWdy.v.y;
//   if (Axisymmetric == 2) {
//     production += lambda.zz*v.x/radius;
//   } else if (Axisymmetric == 1) {
//     production += lambda.zz*v.y/radius;
//   } /* endif */

//   //Determine axisymmetric source terms
//   Temp.rhok = production - 
//               f_beta_star*beta_star*rho*k*omega;
//   Temp.rhoomega = alpha*(omega/max(k, TOLER))*production -
//                   f_beta*beta*rho*omega*omega;

//   return (Temp);


}

/*****************************************************************
 *****************************************************************
 ** Chem2D_pState::Sw -- Chemical Reaction Rate Source Terms.   **
 **                                                             **
 ** Using the Reaction class to get the source terms for the    ** 
 ** specific "Reaction_set".                                    ** 
 *****************************************************************
 *****************************************************************/
Chem2D_cState Chem2D_pState::Sw(int &REACT_SET_FLAG, const int &Flow_Type) const {
  Chem2D_cState NEW;     
  NEW.Vacuum();

  //Adds concentration rate of change for species 1->N
  if (REACT_SET_FLAG != NO_REACTIONS) {
    bool test = negative_speccheck();
    React.omega(NEW,*this,Flow_Type);  
  }
     
  return NEW;

}

/************* Chemical Source Term Jacobian ****************************/
void Chem2D_pState::dSwdU(DenseMatrix &dSwdU) const {
  React.dSwdU(dSwdU,*this,false,flow_type);
}

/*****************************************************************
 *****************************************************************
 ** Chem2D_pState::Sg -- Source Terms for gravitational body    **
 **                      force.                                 **
 **                                                             **
 *****************************************************************
 *****************************************************************/
Chem2D_cState Chem2D_pState::Sg(void) const {
  Chem2D_cState NEW;     
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
void Chem2D_pState::dSgdU(DenseMatrix &dSgdU) const {
  dSgdU(2,0) += gravity_z;
  dSgdU(3,2) += gravity_z;
}

/************* Max Diagonal of Jacobian for CFL **************************/
double Chem2D_pState::dSwdU_max_diagonal(const int &Preconditioned,
					 const double &delta_n) const {

  //this is expensive as I am recalculating the whole Jacobain
  //as above, but its really easy to setup.
  //should change later to only calculate the diagonal terms!!!!

  double max_diagonal =ONE;
  DenseMatrix dSwdU(NUM_VAR_CHEM2D-1,NUM_VAR_CHEM2D-1);
  dSwdU.zero();
  React.dSwdU(dSwdU,*this,true,flow_type);

  if(Preconditioned == 1){
    DenseMatrix Pinv(NUM_VAR_CHEM2D-1,NUM_VAR_CHEM2D-1);
    Low_Mach_Number_Preconditioner_Inverse(Pinv,delta_n);
    dSwdU = Pinv*dSwdU;
  }

  for(int i=0; i < NUM_VAR_CHEM2D-1; i++){
    max_diagonal = max(max_diagonal,fabs(dSwdU(i,i)));
  }
  
  return max_diagonal;
}

/**************************************************************************
********************* CHEM2D_CSTATE CONSTRUCTORS **************************
***************************************************************************/

/**************************************************
  mixture gas constant  J/(kg*K)
***************************************************/
double Chem2D_cState::Rtot() const{
  // = sum ( mass fraction * species gas constant)
  double sum = 0.0;
  for(int i=0; i<ns; i++){
    sum += rhospec[i].c * specdata[i].Rs();
  }
  return (sum/rho);
}

// /**************************************************
//   mixture Heat Capacity (const pressure) J/(kg*K)
// ***************************************************/
// double Chem2D_cState::Cp(void) const{
//   // = sum ( mass fraction * species Cp) 
//   double Temp = T();
//   double sum = 0.0;
//   for(int i=0; i<ns; i++){
//     sum += rhospec[i].c*specdata[i].HeatCapacity_p(Temp);
//   }
//   return (sum/rho);
// }

// /**************************************************
//   mixture Heat Capacity (const volume) J/(kg*K)
// ***************************************************/
// double Chem2D_cState::Cv(void) const{
//   // = sum ( mass fraction * species Cv)  
//   double Temp = T();
//   double sum = 0.0;
//   for(int i=0; i<ns; i++){
//     sum += rhospec[i].c*specdata[i].HeatCapacity_v(Temp);
//   }
//   return (sum/rho);
// }

// /**************************************************
//   mixture Heat Ratio gamma J/(kg*K)
// ***************************************************/
// double Chem2D_cState::g(void) const{
//   // = Cp / Cv  
//   return Cp()/Cv();
// }

/**************************************************
  Specific Internal Energy
 ***************************************************/
double Chem2D_cState::e(void) const{
  // = sum (mass fraction * species e) 
  double sum = 0.0;
  double Temp = T();
  for(int i=0; i<ns; i++){ //(Enthalpy(Temp) - (R/mol_mass)*Temp)
    sum += rhospec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform() 
			 - specdata[i].Rs()*Temp);
  }
  return (sum/rho);
}

double Chem2D_cState::es(void) const{
  // = sum (mass fraction * species e) 
  double sum = 0.0;
  double Temp = T();
  for(int i=0; i<ns; i++){ //(Enthalpy(Temp) - (R/mol_mass)*Temp)
    sum += rhospec[i].c*(specdata[i].Enthalpy(Temp) - specdata[i].Rs()*Temp);
  }
  return (sum/rho);
}

/**************************************************
 Specific Absolute enthalpy
***************************************************/
double Chem2D_cState::h(const double &Temp) const{
  // = sum (mass fraction * species h) 
 double sum = 0.0;  
 for(int i=0; i<ns; i++){
   sum += rhospec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform());
 }
 return (sum/rho);
}

double Chem2D_cState::hs(const double &Temp) const{
  // = sum (mass fraction * species h) 
 double sum = 0.0;  
 for(int i=0; i<ns; i++){
   sum += rhospec[i].c*(specdata[i].Enthalpy(Temp));
 }
 return (sum/rho);
}

/**************************************************
   Derivative of specific enthalpy dh/dT
   actually is just Cp as Cp = (dh/dT)_p
***************************************************/
double Chem2D_cState::hprime(const double &Temp) const{
 double sum = 0.0;  
 for(int i=0; i<ns; i++){
   sum += rhospec[i].c*specdata[i].Enthalpy_prime(Temp);
 }
 return (sum/rho);
}

/************* Mixture Heats of Formation *********/
double Chem2D_cState::heatofform(void) const{ 
  double sum = 0.0;
  for(int i=0; i<ns; i++){ 
    sum += rhospec[i].c*specdata[i].Heatofform();
  }
  return (sum);
}


/**************************************************
  polytropic heat ratio mixture gamma J/(kg*K)
  assuming T=273K as the temperature.
 **************************************************/
double Chem2D_cState::gamma_guess(void) const{
  double sum1 = ZERO;     double sum2 = ZERO;
  double gamma_s = ZERO;  double Temp = 200.0;
  for(int i=0; i<ns; i++){
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
double Chem2D_cState::T(void) const{
  double T = ZERO;
  double RTOT = Rtot();  
  //--------- Initial Guess ------------------------------//
  //using a polytropic gas assumption with gamma@200;
  double Tguess = (gamma_guess() - ONE)*(E - HALF*rhov.sqr()/rho-rhok)/(rho*RTOT);
  //--------- global newtons method to get T ---------------//
  double A = (E - HALF*rhov*rhov/rho -rhok)/rho;
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
      if(T >= Tmax) T = HALF*(Tmax - Tmin);	
      //Bisection
    } else {
      T = HALF*(Tmax - Tmin);
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
    if(debug_level){ 
      cout<<"\nTemperature didn't converge in Chem2D_cState::T(void)";
      cout<<" with polytopic Tguess "<<Tguess<<", or lower than Tmin "<<low_temp_range<<" using "<<T;
    }   
  }

  return T;
} 


/****************************************************
  Speed of sound using 
  a^2 = dip/dirho + p/rho^2( die/dip)^-1
  from eigenvalue analysis using e =f(p,rho)
****************************************************/
double Chem2D_cState::a(void) const{
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
double Chem2D_cState::eddy_viscosity(void) const{
  return (rho*rhok/rhoomega);
}

double Chem2D_cState::Pr_turb(void) const{
  return (0.9);
}

double Chem2D_cState::Sc_turb(void) const{
  return (1.01);
}

double Chem2D_cState::Dm_turb(void) const{
  return (eddy_viscosity()/(rho*Sc_turb()));
}

double Chem2D_cState::omega_sublayer_BC(const double &y) const {
  return (SIX*mu()/(rho*beta*y*y));
}

/**************************************************
  Viscosity 
  using Wilke [1950] formulation
***************************************************/
double Chem2D_cState::mu(void) const{
  double sum =0.0; 
  double Temp = T();

  for(int i=0; i<ns; i++){
    double phi = 0.0;
    for (int j=0; j<ns; j++){
      phi += ((rhospec[j].c/rho) / specdata[j].Mol_mass())*
	pow(1.0 + sqrt(specdata[i].Viscosity(Temp)/specdata[j].Viscosity(Temp))*
	    pow(specdata[j].Mol_mass()/specdata[i].Mol_mass(),0.25),2.0)/
       sqrt(8.0*(1.0 +specdata[i].Mol_mass()/specdata[j].Mol_mass()));
    }
    sum += ((rhospec[i].c/rho)* specdata[i].Viscosity(Temp) ) / 
      (specdata[i].Mol_mass() * phi);
  }  
  
  return sum;
}

/******************************************************
 Calculating the thermal diffusion component of 
 the heat flux vector (qflux)

  sum( hs * Ds * grad cs)
*******************************************************/
Vector2D Chem2D_cState::thermal_diffusion(const double &Temp) const{
  Vector2D sum;
  sum.zero();
  //double Temp = T();
  //problems with Species overloaded operators
  for(int i=0; i<ns; i++){ 
    sum  +=   (specdata[i].Enthalpy(Temp) + specdata[i].Heatofform())
            * rhospec[i].diffusion_coef*rhospec[i].gradc;
  }
  return sum/(rho*rho);
}

/*****************************************************************
 * Viscous fluxes  (laminar flow)                                * 
 * Viscous fluxes  (turbulent flows) are defined in single block * 
 ****************************************************************/
Chem2D_cState Chem2D_cState::Viscous_Flux_x(const Chem2D_pState &dWdx) const{
 
  Chem2D_cState temp;

  temp[1] = ZERO;
  temp[2] = tau.xx;
  temp[3] = tau.xy;
  temp[4] = - qflux.x + v().x*tau.xx + v().y*tau.xy;		

  //rho * Diffusion_Coef * grad cn 
  for( int i=0; i<ns; i++){
    temp.rhospec[i].c = (rhospec[i].diffusion_coef * rhospec[i].gradc.x)/rho; 
  }

  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
      flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
    
    temp[2] += lambda.xx; 
    temp[3] += lambda.xy;
    temp[4] += - theta.x + v().x*lambda.xx + v().y*lambda.xy
      + (mu()+eddy_viscosity()*sigma_star)*dWdx.k;
    
    temp[5] = (mu()+eddy_viscosity()*sigma_star)*dWdx.k;
    temp[6] = (mu()+eddy_viscosity()*sigma_star)*dWdx.omega;
    
    double Dm_t = Dm_turb();
    for( int i=0; i<ns; i++){
      temp.rhospec[i].c += Dm_t*rhospec[i].gradc.x; 
    }
  }
 
  return(temp);  
}

Chem2D_cState Chem2D_cState::Viscous_Flux_y(const Chem2D_pState &dWdy) const {
  Chem2D_cState temp;

  temp[1] = ZERO;
  temp[2] = tau.xy; 
  temp[3] = tau.yy;
  temp[4] = - qflux.y + v().x*tau.xy + v().y*tau.yy;		
  //rho * Diffusion_Coef * grad cn 
  for( int i=0; i<ns; i++){
    temp.rhospec[i].c = (rhospec[i].diffusion_coef * rhospec[i].gradc.y)/rho;     
  }

  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
      flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
     temp[2] += lambda.xy; 
     temp[3] += lambda.xy;
     temp[4] += - theta.y + v().x*lambda.xy + v().y*lambda.yy
                + (mu()+eddy_viscosity()*sigma_star)*dWdy.k;
     temp[5] = (mu()+eddy_viscosity()*sigma_star)*dWdy.k;
     temp[6] = (mu()+eddy_viscosity()*sigma_star)*dWdy.omega;
     double Dm_t = Dm_turb();
     for( int i=0; i<ns; i++){
        temp.rhospec[i].c += Dm_t*rhospec[i].gradc.y; 
     }
  }

  return(temp);  
}

/*********************************************************************
 *************** Overloaded Operators ********************************
 *********************************************************************/

//----------------- Addition -----------------------------//
Chem2D_cState Chem2D_cState::operator +(const Chem2D_cState &U) const{ 
  
  if(ns == U.ns){ //check that species are equal   
    Chem2D_cState Temp(U.rho,U.rhov,U.E, U.rhok, U.rhoomega);
    Temp.Copy(*this);
    Temp += U;
    return Temp;
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
}

//------------------ Subtraction ------------------------//
Chem2D_cState Chem2D_cState::operator -(const Chem2D_cState &U) const{
  if(ns == U.ns){ //check that species are equal
    Chem2D_cState Temp(U.rho,U.rhov,U.E, U.rhok, U.rhoomega);
    Temp.Copy(*this);
    Temp -= U;
    return Temp;
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
}

//---------------- Scalar Multiplication ------------------//
Chem2D_cState Chem2D_cState::operator *(const double &a) const{
  Chem2D_cState Temp(rho,rhov,E, rhok, rhoomega);
  Temp.Copy(*this);
  Temp.rho = rho*a;  Temp.rhov = rhov*a; Temp.E = E*a;
  Temp.rhok = rhok*a; Temp.rhoomega = rhoomega*a;
  for( int i=0; i<ns; i++){
    Temp.rhospec[i] = rhospec[i]*a;
  } 

  Temp.tau = tau*a;
  Temp.qflux = qflux*a;
  Temp.lambda = lambda*a;
  Temp.theta = theta*a;

  return(Temp);
}

Chem2D_cState operator *(const double &a, const Chem2D_cState &U){
  Chem2D_cState Temp;
  Temp.rho = U.rho*a;  Temp.rhov = U.rhov*a; Temp.E = U.E*a;
  Temp.rhok = U.rhok*a;Temp.rhoomega = U.rhoomega*a;
  for( int i=0; i<U.ns; i++){
    Temp.rhospec[i] = U.rhospec[i]*a;
  } 
  return(Temp);
}
//--------------- Scalar Division ------------------------//
Chem2D_cState Chem2D_cState::operator /(const double &a) const {
  Chem2D_cState Temp(rho,rhov,E,rhok,rhoomega);
  Temp.Copy(*this);
  Temp.rho = rho/a; Temp.rhov = rhov/a; Temp.E = E/a;
  Temp.rhok = rhok/a; Temp.rhoomega = rhoomega/a;
  for(int i=0; i<ns; i++){
     Temp.rhospec[i] = rhospec[i]/a; 
  } 

  Temp.tau = tau/a;
  Temp.qflux = qflux/a;
  Temp.lambda = lambda/a;
  Temp.theta = theta/a;

  return(Temp);
}

//----------------- Inner Product ------------------------//
double Chem2D_cState::operator *(const Chem2D_cState &U) const{
  double sum=0.0;
  if(ns == U.ns){ //check that species are equal
    for(int i=0; i<ns; i++){
      sum += rhospec[i]*U.rhospec[i];
    }  
    return (rho*U.rho + rhov*U.rhov + E*U.E +rhok*U.rhok+rhoomega*U.rhoomega+ sum);
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
}

//----------- solution state product operator ------------//
Chem2D_cState Chem2D_cState::operator ^( const Chem2D_cState &U) const {
  if(ns == U.ns){ //check that species are equal
    Chem2D_cState Temp(rho,rhov,E,rhok, rhoomega);
    Temp.Copy(*this);
    Temp.rho = rho*U.rho;
    Temp.rhov.x = rhov.x*U.rhov.x;
    Temp.rhov.y = rhov.y*U.rhov.y;
    Temp.E = E*U.E;
    Temp.rhok= rhok*U.rhok;
    Temp.rhoomega = rhoomega*U.rhoomega;
    for(int i=0; i<ns; i++){
      Temp.rhospec[i] = rhospec[i]*U.rhospec[i];
    }  
    return(Temp);
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
}

//----------------- Assignment ----------------------------//
Chem2D_cState& Chem2D_cState::operator =(const Chem2D_cState &U){
  //self assignment protection
  if( this != &U){   
    //copy assignment
    rho = U.rho;
    rhov = U.rhov; 
    E = U.E; 
    rhok = U.rhok;
    rhoomega = U.rhoomega;
    if ( ns == U.ns){
      for(int i=0; i<ns; i++){
	rhospec[i] = U.rhospec[i];
      }

      tau = U.tau;
      qflux = U.qflux; 
      lambda = U.lambda;
      theta = U.theta; 

    } else {
      cerr<<"\n Mismatch in number of species \n ";
    }  
  }
  return (*this);
}


/********************************************************
 * Chem2D_cState -- Shortcut arithmetic operators.     *
 ********************************************************/
Chem2D_cState& Chem2D_cState::operator +=(const Chem2D_cState &U){
  rho += U.rho;
  rhov += U.rhov; 
  E += U.E;
  rhok += U.rhok;
  rhoomega += U.rhoomega;
  for( int i=0; i<ns; i++){
    rhospec[i] += U.rhospec[i];
  } 

  tau += U.tau;
  qflux += U.qflux; 
  lambda += U.lambda;
  theta += U.theta; 

  return (*this);
}

Chem2D_cState& Chem2D_cState::operator -=(const Chem2D_cState &U) {
  rho -= U.rho;
  rhov -= U.rhov;
  E -= U.E;
  rhok -= U.rhok;
  rhoomega -= U.rhoomega;
  for(int i=0; i<ns; i++){
    rhospec[i] -= U.rhospec[i];
  }  

  tau -= U.tau;
  qflux -= U.qflux; 
  lambda -= U.lambda;
  theta -= U.theta; 

  return (*this); 
}

Chem2D_cState& Chem2D_cState::operator *=(const double &a) {
  rho *= a;
  rhov.x *= a;
  rhov.y *= a;
  E *= a;
  rhok *= a;
  rhoomega *= a;
  for (int i = 0; i < ns; i++) rhospec[i] *= a;
  tau *= a;
  qflux *= a;
  lambda *= a;
  theta *= a;
  return *this;
}

Chem2D_cState& Chem2D_cState::operator /=(const double &a) {
  rho /= a;
  rhov.x /= a;
  rhov.y /= a;
  E /= a;
  rhok /= a;
  rhoomega /= a;
  for (int i = 0; i < ns; i++) rhospec[i] /= a;
  tau /= a;
  qflux /= a;
  lambda /= a;
  theta /= a;
  return *this;
}

/********************************************************
 * Chem2D_cState -- Unary arithmetic operators.        *
 ********************************************************/
Chem2D_cState operator -(const Chem2D_cState &U) {
  Species *spt= new Species[U.ns];
  for(int i=0; i<U.ns; i++){
    spt[i] = -U.rhospec[i]; 
  }  

  Chem2D_cState Temp(-U.rho,-U.rhov, -U.E, -U.rhok, -U.rhoomega, spt);
  Temp.tau = -U.tau;
  Temp.qflux = -U.qflux;
  Temp.lambda = -U.lambda;
  Temp.theta = -U.theta;

  delete[] spt;
  return(Temp);
}

/********************************************************
 * Chem2D_cState -- Relational operators.              *
 ********************************************************/
int operator ==(const Chem2D_cState &U1, const Chem2D_cState &U2) {
  if(U1.ns == U2.ns){ //check that species are equal
    bool Temp;
    for(int i=0; i<U1.ns; i++){
      if( U1.rhospec[i] == U2.rhospec[i] ){
	Temp = true;
      } else {
	Temp = false;
	break;
      }  
      return (U1.rho == U2.rho && U1.rhov == U2.rhov && U1.E == U2.E 
	      && U1.rhok == U2.rhok && U1.rhoomega == U2.rhoomega &&
	      U1.tau == U2.tau && U1.qflux == U2.qflux&&
	      U1.lambda == U2.lambda && U1.theta == U2.theta
	      &&Temp == true);
    }
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
}

int operator !=(const Chem2D_cState &U1, const Chem2D_cState &U2) {
   if(U1.ns == U2.ns){ //check that species are equal
    bool Temp = true;
    for(int i=0; i<U1.ns; i++){
      if( U1.rhospec[i] != U2.rhospec[i] ){
	Temp = false;
	break;
      } 
     return (U1.rho != U2.rho || U1.rhov != U2.rhov || U1.E != U2.E 
	     ||  U1.rhok != U2.rhok ||  U1.rhoomega != U2.rhoomega ||
	     U1.tau != U2.tau || U1.qflux != U2.qflux
	     ||U1.lambda != U2.lambda || U1.theta != U2.theta
	     || Temp != true);
    }
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
}

/********************************************************
 * Chem2D_cState -- Input-output operators.            *
 ********************************************************/
ostream &operator << (ostream &out_file, const Chem2D_cState &U) {
  //out_file.precision(20);
  out_file.setf(ios::scientific);
  out_file << " " << U.rho  << " " << U.rhov.x << " " << U.rhov.y << " " << U.E<< " " << U.rhok<< " " << U.rhoomega;
  for( int i=0; i<U.ns; i++){
    out_file<<" "<<U.rhospec[i];
  } 
  out_file << " " <<U.qflux<< " " <<U.tau << " " <<U.theta<< " " <<U.lambda;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

istream &operator >> (istream &in_file, Chem2D_cState &U) {
  in_file.setf(ios::skipws);
  in_file >> U.rho >> U.rhov.x >> U.rhov.y >> U.E >> U.rhok >> U.rhoomega;
  //U.set_initial_values();
  for( int i=0; i<U.ns; i++){
    in_file>>U.rhospec[i]; 
  } 
  in_file >>U.qflux>>U.tau >>U.theta>>U.lambda;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/*******************************************************************
 *******************************************************************
 *  Low Mach Number Preconditioner Based on Weiss & Smith (1995)   *
 *                                                                 *
 *******************************************************************
 *******************************************************************/
void Chem2D_cState::Low_Mach_Number_Preconditioner(DenseMatrix  &P,
						   const double &deltax) const {
  Chem2D_pState NEW = W();
  NEW.Low_Mach_Number_Preconditioner(P,deltax);
}

void Chem2D_cState::Low_Mach_Number_Preconditioner_Inverse(DenseMatrix  &Pinv,
							   const double &deltax) const {
 
  Chem2D_pState NEW = W();
  NEW.Low_Mach_Number_Preconditioner_Inverse(Pinv,deltax);
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
Chem2D_pState Reflect(const Chem2D_pState &W,
		      const Vector2D &norm_dir) {

    double ur, vr, u, v;
    double cos_angle, sin_angle;
    Chem2D_pState Temp(W);
 
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
Chem2D_pState Free_Slip(const Chem2D_pState &Win,
			const Chem2D_pState &Wout,
			const Vector2D &norm_dir,
			const int &TEMPERATURE_BC_FLAG) { 
  
  /***** MOHAMMED et al *******/
  //Chem2D_pState Temp(Wout);
  //Temp.v = Win.v;
  
  /**** DAY & BELL ********/
  Chem2D_pState Temp;
  Temp = Reflect(Win,norm_dir);
 
  //Fixed Temperature
  if(TEMPERATURE_BC_FLAG == FIXED_TEMPERATURE_WALL){
   Temp.rho = Temp.p/(Temp.Rtot()*Wout.T());
  }
  
  return Temp;
}

/*********** NO SLIP - ************************************/
Chem2D_pState No_Slip(const Chem2D_pState &Win,
		      const Chem2D_pState &Wout,
		      const Vector2D &norm_dir,
		      const int &TEMPERATURE_BC_FLAG) {  
  
  return(Moving_Wall(Win,Wout,norm_dir,ZERO,TEMPERATURE_BC_FLAG));

}

/************ Moving_Wall **********************************/
Chem2D_pState Moving_Wall(const Chem2D_pState &Win,
			  const Chem2D_pState &Wout,
			  const Vector2D &norm_dir, 
			  const double &wall_velocity,
			  const int &TEMPERATURE_BC_FLAG) {

  double ur, vr;
  double q_xr,q_yr, theta_xr,theta_yr;  
  double cos_angle, sin_angle;
  Chem2D_pState Temp;
  Temp.Copy(Win);
  
  /* Determine the direction cosine's for the frame
     rotation. */
  
  cos_angle = norm_dir.x; 
  sin_angle = norm_dir.y;

  /* Apply the frame rotation and calculate the primitive
     solution state variables in the local rotated frame
     defined by the unit normal vector. */
  ur = Win.v.x*cos_angle +
    Win.v.y*sin_angle;
  vr = - Win.v.x*sin_angle +
    Win.v.y*cos_angle;
  
   /* Use the Roeaveraged value to find the correct tangential velocity */
   ur = -ur;
   vr = -TWO*wall_velocity - vr;

   /* Rotate back to the original Cartesin reference frame. */  
   Temp.v.x = ur*cos_angle - vr*sin_angle;
   Temp.v.y = ur*sin_angle + vr*cos_angle;
   //   Temp.k = -Win.k; //turbulent kinetic energy  vanishes on the wall (no-slip)
   
//   /* Fixed Wall Temperature ISOTHERMAL */
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
Chem2D_pState BC_Flame_Inflow(const Chem2D_pState &Wi,
			      const Chem2D_pState &Wo,
			      const Chem2D_pState &Woutlet,
			      const Vector2D &norm_dir){
  Chem2D_pState Wnew;

  //fixed
  Wnew.Copy(Wo);

  //Constant Extrapolate 
  //Wnew.Copy(Wi);
 
  //Calculate upstream velocity to balance
  //flame ie. mass flow rate 
  // rho_1*u_1 = rho_2*u_2
  //Wnew.v.x =  Woutlet.rho*Woutlet.v.x/Wo.rho;
  Wnew.v.y = ZERO;
  
  Wnew.v.x = Wi.v.x + 0.5*(Woutlet.rho*Woutlet.v.x/Wo.rho - Wi.v.x);

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
Chem2D_pState BC_Flame_Outflow(const Chem2D_pState &Wi, 
			       const Chem2D_pState &Wo,
			       const Chem2D_pState &Winlet,
			       const Vector2D &norm_dir){

  //Constant extrapolate ( zero gradient)
  Chem2D_pState Wnew(Wi);
 
  //Calculate Pressure assuming constant mass flow rate
  //and Wo.p == Winput.p (constant pressure initial condition)
  double sum = Wi.rho*Wi.v.x*(Wi.v.x - Winlet.v.x);
  if( sum < ZERO){
    Wnew.p = Wo.p;
  } else {
    //Wnew.p = Wo.p - sum;
    Wnew.p = Wi.p + 0.5*( (Wo.p-sum) - Wi.p);
  }

   return Wnew;
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
Chem2D_pState BC_Characteristic_Pressure(const Chem2D_pState &Wi,
					 const Chem2D_pState &Wo,
					 const Vector2D &norm_dir) {

  Chem2D_pState Wnew(Wi);  
  Wnew.p = Wo.p;
  
  if(Wnew.v.y < ZERO){ 
    Wnew.v.y = ZERO;
  }
  return Wnew;

//   Chem2D_pState Wi_rotated, Wo_rotated, Wnew;
//   double mi, ab, ub_rotated, vb_rotated;
//   double cos_angle, sin_angle;
  
//   /* Determine the direction cosine's for the frame
//      rotation. */
  
//   cos_angle = norm_dir.x; 
//   sin_angle = norm_dir.y;
  
//   /* Apply the frame rotation and evaluate interior and 
//      imposed boundary solution states in the local rotated 
//        frame defined by the unit normal vector. */
//   Wi_rotated.Copy(Wi);
//   Wo_rotated.Copy(Wo);


//   Wi_rotated.v.x = Wi.v.x*cos_angle +
//     Wi.v.y*sin_angle;
//   Wi_rotated.v.y = - Wi.v.x*sin_angle +
//     Wi.v.y*cos_angle;

 
//   Wo_rotated.v.x = Wo.v.x*cos_angle +
//     Wo.v.y*sin_angle;
//   Wo_rotated.v.y = - Wo.v.x*sin_angle +
//     Wo.v.y*cos_angle;
 

//   /* Determine the Mach number at the interior node. */
//   mi = Wi_rotated.v.x/Wi_rotated.a();
  
//   /* Boundary condition for supersonic outflow. */
//   if (mi >= ONE) {
//     Wnew.Copy(Wi);
    
//     /* Boundary condition for subsonic outflow. 
//        Pressure specified. */
//   } else if (mi >= ZERO) {
//     Wnew.Copy(Wi_rotated);
//     Wnew.p = Wo_rotated.p;
//     Wnew.rho = Wi_rotated.rho*pow(Wnew.p/Wi_rotated.p, ONE/Wi_rotated.g());
    
//     ab = Wnew.a();
//     ub_rotated = Wi_rotated.v.x + TWO*(Wi_rotated.a()-ab)/(Wi_rotated.g() -ONE);
//     vb_rotated = Wi_rotated.v.y;
    
//     /* Boundary condition for subsonic inflow. 
//        Pressure specified. */
//   } else if (mi >= -ONE) {
//     Wnew.Copy(Wo_rotated);
    
//     ab = Wnew.a();
//     ub_rotated = Wi_rotated.v.x + TWO*(Wi_rotated.a()-ab)*(Wo_rotated.g() -ONE);
//     vb_rotated = Wo_rotated.v.y;
     
//     /* Boundary condition for supersonic inflow.  */
//   } else {
//     Wnew.Copy(Wo);
//   } 

//   if( fabs(mi) <=ONE){
//     //rotate frame back to cartesian
//     Wnew.v.x = ub_rotated*cos_angle - vb_rotated*sin_angle;
//     Wnew.v.y = ub_rotated*sin_angle + vb_rotated*cos_angle;
//   }

//   /* Return boundary solution state. */
//   return (Wnew);
  
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
Chem2D_pState RinglebFlow(const Chem2D_pState &Wdum,
			   const Vector2D X) {

//   Chem2D_pState W;
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
Chem2D_pState ViscousChannelFlow(const Chem2D_pState &Wdum,
				  const Vector2D X,
				  const double Vwall,
				  const double dp) {

  Chem2D_pState W;
  W.Copy(Wdum);  //use same species and viscosity

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
Chem2D_pState FlatPlate(const Chem2D_pState &Winf,
			 const Vector2D X,
			 double &eta,
			 double &f,
			 double &fp,
			 double &fpp) {

  Chem2D_pState W;
  double fo, dn, sign, k1, k2, k3, k4;

  // Initialize variables.
  W.Vacuum(); W.rho = Winf.rho; W.p = Winf.p;
  for(int i=0; i<W.ns; i++){
    W.spec[i].c = Winf.spec[i].c;
  } 

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
  W.v.y = HALF*(eta*fp-f);

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
Chem2D_pState RoeAverage(const Chem2D_pState &Wl,
			 const Chem2D_pState &Wr) {

    double Hl, Hr, srhol, srhor;
    double Ha, ha;
    Chem2D_pState Temp;

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
    //Turbulence quantities
    Temp.k = (srhol*Wl.k + srhor*Wr.k)/(srhol+srhor);
    Temp.omega= (srhol*Wl.omega + srhor*Wr.omega)/(srhol+srhor);
    for(int i=0; i<Wl.ns; i++){
      Temp.spec[i].c = (srhol*Wl.spec[i].c + srhor*Wr.spec[i].c)/(srhol+srhor);
    }
 
    Ha = (srhol*Hl+srhor*Hr)/(srhol+srhor);
    ha = Ha - HALF*(sqr(Temp.v.x)+sqr(Temp.v.y));

    double TEMP = Temp.T(ha);
    Temp.p = Temp.rho*TEMP*Temp.Rtot();
   
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
Chem2D_cState FluxHLLE_x(const Chem2D_pState &Wl,
			 const Chem2D_pState &Wr,
			 const int &Preconditioning) {

    double wavespeed_l, wavespeed_r;
    Chem2D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
    Chem2D_cState Flux, dUrl;
    int NUM_VAR_CHEM2D = Wl.NUM_VAR_CHEM2D;

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

    //wavespeed_r = max(lambdas_r[NUM_VAR_CHEM2D],
    //                  lambdas_a[NUM_VAR_CHEM2D]);
    wavespeed_r = max(lambdas_r[4],lambdas_a[4]); //MARKTHIS XINFENG
    //   wavespeed_r = max(lambdas_r[NUM_VAR_CHEM2D-lambdas_r.ns],
    //                       lambdas_a[NUM_VAR_CHEM2D-lambdas_a.ns]); 
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

Chem2D_cState FluxHLLE_x(const Chem2D_cState &Ul,
			 const Chem2D_cState &Ur,
			 const int &Preconditioning) {
   return (FluxHLLE_x(Ul.W(), Ur.W(),Preconditioning));
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
Chem2D_cState FluxHLLE_n(const Chem2D_pState &Wl,
			 const Chem2D_pState &Wr,
			 const Vector2D &norm_dir,
			 const int &Preconditioning) {

    double cos_angle, sin_angle;
    Chem2D_pState Wl_rotated, Wr_rotated;
    Chem2D_cState Flux, Flux_rotated;
    Flux.Vacuum();
    /* Determine the direction cosine's for the frame
       rotation. */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    /* Apply the frame rotation and evaluate left and right
       solution states in the local rotated frame defined
       by the unit normal vector. */
    Wl_rotated.Copy(Wl);
    Wr_rotated.Copy(Wr);


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

    Flux_rotated = FluxHLLE_x(Wl_rotated, Wr_rotated, Preconditioning);

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

Chem2D_cState FluxHLLE_n(const Chem2D_cState &Ul,
			 const Chem2D_cState &Ur,
			 const Vector2D &norm_dir,
			 const int &Preconditioning) {
  return (FluxHLLE_n(Ul.W(), Ur.W(),norm_dir,Preconditioning));
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
Chem2D_cState FluxLinde(const Chem2D_pState &Wl,
			const Chem2D_pState &Wr) {

    double wavespeed_l, wavespeed_r, wavespeed_m, rhoa, ca, dU, alpha;
    Chem2D_pState Wa, lambdas_l, lambdas_r, lambdas_a;
    Chem2D_cState Flux, dFrl, dUrl, dFwave;
    int NUM_VAR_CHEM2D = Wl.NUM_VAR_CHEM2D;

    /* Evaluate the Roe-average primitive solution state. */   
    Wa = RoeAverage(Wl, Wr);

    /* Evaluate the left, right, and average state eigenvalues. */
    lambdas_l = Wl.lambda_x();
    lambdas_r = Wr.lambda_x();
    lambdas_a = Wa.lambda_x();
   
    /* Determine the intermediate state flux. */
    wavespeed_l = min(lambdas_l[1],
                      lambdas_a[1]);
  //   wavespeed_r = max(lambdas_r[NUM_VAR_CHEM2D-lambdas_r.ns],
//                       lambdas_a[NUM_VAR_CHEM2D-lambdas_a.ns]);

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
        ca = Wa.a();
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

Chem2D_cState FluxLinde(const Chem2D_cState &Ul,
	      	         const Chem2D_cState &Ur) {
   return (FluxLinde(Ul.W(), Ur.W()));
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
Chem2D_cState FluxLinde_n(const Chem2D_pState &Wl,
	      	           const Chem2D_pState &Wr,
                           const Vector2D &norm_dir) {

    double cos_angle, sin_angle;
    Chem2D_pState Wl_rotated, Wr_rotated;
    Chem2D_cState Flux, Flux_rotated;

    /* Determine the direction cosine's for the frame
       rotation. */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    /* Apply the frame rotation and evaluate left and right
       solution states in the local rotated frame defined
       by the unit normal vector. */
    Wl_rotated.Copy(Wl);
    Wr_rotated.Copy(Wr);

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

    Flux_rotated = FluxLinde(Wl_rotated, Wr_rotated);

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

Chem2D_cState FluxLinde_n(const Chem2D_cState &Ul,
			  const Chem2D_cState &Ur,
			  const Vector2D &norm_dir) {
  return (FluxLinde_n(Ul.W(), Ur.W(), norm_dir));
}


/*********************************************************
 * Routine: Rotate                                       *
 *                                                       *
 * This function returns the solution in the lcoal       *
 * rotated frame.                                        *
 *                                                       *
 *********************************************************/
Chem2D_pState Rotate(const Chem2D_pState &W,
                      const Vector2D &norm_dir) {
  Chem2D_pState W_rotated;
  double cos_angle = norm_dir.x;  
  double sin_angle = norm_dir.y;

  W_rotated.rho   = W.rho;
  W_rotated.v.x =   W.v.x*cos_angle + W.v.y*sin_angle;
  W_rotated.v.y = - W.v.x*sin_angle + W.v.y*cos_angle;
  W_rotated.p   = W.p;
  W_rotated.k   = W.k;
  W_rotated.omega   = W.omega;
  for( int i=0; i<W.ns; i++){
     W_rotated.spec[i] = W.spec[i];
  } 
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
Vector2D HLLE_wavespeeds(const Chem2D_pState &Wl,
                         const Chem2D_pState &Wr,
			 const Vector2D &norm_dir) {

    Vector2D wavespeed;
    Chem2D_pState Wa_n, lambdas_l, lambdas_r, lambdas_a, Wl_n, Wr_n;  
    int NUM_VAR_CHEM2D = (Wl.NUM_VAR_CHEM2D );
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
                      lambdas_a[1]);
    wavespeed.y = max(lambdas_r[NUM_VAR_CHEM2D],
                      lambdas_a[NUM_VAR_CHEM2D]);
 
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
Chem2D_pState WaveSpeedPos(const Chem2D_pState &lambdas_a,
			   const Chem2D_pState &lambdas_l,
			   const Chem2D_pState &lambdas_r) {
   Chem2D_pState NEW;   
   for(int i=1; i<=lambdas_a.NUM_VAR_CHEM2D; i++){
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
Chem2D_pState WaveSpeedNeg(const Chem2D_pState &lambdas_a,
			   const Chem2D_pState &lambdas_l,
			   const Chem2D_pState &lambdas_r) {
  Chem2D_pState NEW;   
  for(int i=1; i<=lambdas_a.NUM_VAR_CHEM2D; i++){
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
Chem2D_pState WaveSpeedAbs(const Chem2D_pState &lambdas_a,
			   const Chem2D_pState &lambdas_l,
			   const Chem2D_pState &lambdas_r) {
   Chem2D_pState NEW;   
   for(int i=1; i<=lambdas_a.NUM_VAR_CHEM2D; i++){
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
 *                                                       *
 ********************************************************/
Chem2D_pState HartenFixPos(const Chem2D_pState &lambdas_a,
			   const Chem2D_pState &lambdas_l,
			   const Chem2D_pState &lambdas_r) {
  Chem2D_pState NEW;
  NEW.rho = HartenFixPos(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
  NEW.v.x = HALF*(lambdas_a[2]+fabs(lambdas_a[2]));
  NEW.v.y = HALF*(lambdas_a[3]+fabs(lambdas_a[3]));
  NEW.p = HartenFixPos(lambdas_a[4],lambdas_l[4],lambdas_r[4]);
  NEW.k = fabs(lambdas_a[5]);
  NEW.omega = fabs(lambdas_a[6]);
  
  for( int i=7; i<=NEW.NUM_VAR_CHEM2D; i++){
    NEW.spec[i-7].c = HALF*(lambdas_a[i]+fabs(lambdas_a[i]));
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
Chem2D_pState HartenFixNeg(const Chem2D_pState &lambdas_a,
			   const Chem2D_pState &lambdas_l,
			   const Chem2D_pState &lambdas_r) {
  Chem2D_pState NEW;
  NEW.rho = HartenFixNeg(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
  NEW.v.x = HALF*(lambdas_a[2]-fabs(lambdas_a[2]));
  NEW.v.y = HALF*(lambdas_a[3]-fabs(lambdas_a[3]));
  NEW.p = HartenFixNeg(lambdas_a[4],lambdas_l[4],lambdas_r[4]);
 
  NEW.k = fabs(lambdas_a[5]);
  NEW.omega = fabs(lambdas_a[6]);
  
  for( int i=7; i<=NEW.NUM_VAR_CHEM2D; i++){
    NEW.spec[i-7].c = HALF*(lambdas_a[i]-fabs(lambdas_a[i]));
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
Chem2D_pState HartenFixAbs(const Chem2D_pState &lambdas_a,
			   const Chem2D_pState &lambdas_l,
			   const Chem2D_pState &lambdas_r) {
  Chem2D_pState NEW;
  NEW.rho = HartenFixAbs(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
  NEW.v.x = fabs(lambdas_a[2]);
  NEW.v.y = fabs(lambdas_a[3]);
  NEW.p = HartenFixAbs(lambdas_a[4],lambdas_l[4],lambdas_r[4]);
 
  NEW.k = fabs(lambdas_a[5]);
  NEW.omega = fabs(lambdas_a[6]);
  
  for( int i=7; i<=NEW.NUM_VAR_CHEM2D; i++){
    NEW.spec[i-7].c = fabs(lambdas_a[i]);
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
Chem2D_cState FluxRoe_x(const Chem2D_pState &Wl,
			const Chem2D_pState &Wr,
			const int &Preconditioning,
			const double &deltax) {



    Chem2D_pState Wa, dWrl, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
    Chem2D_cState Flux;

 
    /* Evaluate the Roe-average primitive solution state. */    
    Wa = RoeAverage(Wl, Wr);

    /* Evaluate the jumps in the primitive solution states. */
    dWrl = Wr-Wl;

    /* Evaluate the left, right, and average state eigenvalues. */
    if(Preconditioning == 0){
      lambdas_l = Wl.lambda_x();
      lambdas_r = Wr.lambda_x();
      lambdas_a = Wa.lambda_x();

    
      /* Determine the intermediate state flux. */
      if (Wa.v.x >= ZERO) {
        Flux = Wl.Fx();   

// 	wavespeeds = WaveSpeedNeg(lambdas_a,
// 				  lambdas_l,
// 				  lambdas_r);

        wavespeeds = HartenFixNeg(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
	

        for (int i=1 ; i <= Wl.NUM_VAR_CHEM2D; i++) {
	  if (wavespeeds[i] < ZERO) {
 	    Flux += wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);

	  }
        } 
      } else {
        Flux = Wr.Fx();
//  	wavespeeds = WaveSpeedPos(lambdas_a,
// 				  lambdas_l,
// 				  lambdas_r);
        wavespeeds = HartenFixPos(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
        for (int i=1; i <= Wl.NUM_VAR_CHEM2D; i++) {
	  if (wavespeeds[i] > ZERO) {
	    Flux -= wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
          }
        } 
      } 
   
      /******* LOW MACH NUMBER PRECONDITIONING ********************/
      /* Evaluate the left, right, and average state eigenvalues. */
      } else if(Preconditioning == 1){
	
      //calculating Mr^2 and passing to save computation,
      //not conceptually nice but saves from recalculating

      double MR2a = Wa.Mr2(deltax); 
      lambdas_l = Wl.lambda_preconditioned_x(Wl.Mr2(deltax)); 
      lambdas_r = Wr.lambda_preconditioned_x(Wr.Mr2(deltax));
      lambdas_a = Wa.lambda_preconditioned_x(MR2a);
                                                                    
      /* Evaluate the jumps in the primitive solution states. */
//       wavespeeds = WaveSpeedAbs(lambdas_a,
// 				lambdas_l,
// 				lambdas_r);
      wavespeeds = HartenFixAbs(lambdas_a,
				lambdas_l,
				lambdas_r);
          
      DenseMatrix P(Wa.NUM_VAR_CHEM2D-1,Wa.NUM_VAR_CHEM2D-1);        
      /* Evaluate the low-Mach-number local preconditioner for the Roe-averaged state. */  
    
      Wa.Low_Mach_Number_Preconditioner(P,deltax);
                                                                                     
      /* Determine the intermediate state flux. */                                                
      Flux = HALF*(Wl.Fx()+Wr.Fx()); 
      Chem2D_cState Flux_dissipation(ZERO);   
    
      for ( int i = 1 ; i <= Wa.NUM_VAR_CHEM2D ; i++ ) {

	
	Flux_dissipation -= HALF*wavespeeds[i]*(Wa.lp_x_precon(i,MR2a)*dWrl)*Wa.rc_x_precon(i,MR2a);
    
      }
  
      for ( int i = 1 ; i < Wa.NUM_VAR_CHEM2D ; i++ ) {
	for ( int j = 1 ; j < Wa.NUM_VAR_CHEM2D ; j++ ) {
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

Chem2D_cState FluxRoe_x(const Chem2D_cState &Ul,
			const Chem2D_cState &Ur,
			const int &Preconditioning, 
			const double &deltax) {
   return (FluxRoe_x(Ul.W(), Ur.W(),Preconditioning,deltax));
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
Chem2D_cState FluxRoe_n(const Chem2D_pState &Wl,
			const Chem2D_pState &Wr,
			const Vector2D &norm_dir,
			const int &Preconditioning,
			const double &delta_n ) {

 
  double cos_angle, sin_angle, delta_rotated;
  Chem2D_pState Wl_rotated, Wr_rotated;
  Chem2D_cState Flux, Flux_rotated;


    /* Determine the direction cosine's for the frame
       rotation. */

    cos_angle = norm_dir.x; 
    sin_angle = norm_dir.y;

    /* Apply the frame rotation and evaluate left and right
       solution states in the local rotated frame defined
       by the unit normal vector. */
    Wl_rotated.Copy(Wl);
    Wr_rotated.Copy(Wr);
   
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
  
    Flux_rotated = FluxRoe_x(Wl_rotated, Wr_rotated, 
			     Preconditioning,delta_n);

    
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

Chem2D_cState FluxRoe_n(const Chem2D_cState &Ul,
			const Chem2D_cState &Ur,
			const Vector2D &norm_dir, 
			const int &Preconditioning,
			const double &delta_n) {
    return (FluxRoe_n(Ul.W(), Ur.W(), norm_dir,
		      Preconditioning,delta_n));
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
Chem2D_cState Viscous_FluxArithmetic_n(const Chem2D_cState &Ul,
   			               const Chem2D_pState &dWdx_l,
			               const Chem2D_pState &dWdy_l,
				       const Chem2D_cState &Ur,
   			               const Chem2D_pState &dWdx_r,
			               const Chem2D_pState &dWdy_r,
				       const Vector2D &norm_dir) {

  Chem2D_cState Fx, Fy;
  Vector2D i = Vector2D(1,0), j = Vector2D(0,1);

  Fx = HALF*(Ul.Viscous_Flux_x(dWdx_l) + Ur.Viscous_Flux_x(dWdx_r));
  Fy = HALF*(Ul.Viscous_Flux_y(dWdy_l) + Ur.Viscous_Flux_y(dWdy_r));

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
Chem2D_cState Viscous_Flux_n(const Chem2D_pState &W,
			     const Chem2D_pState &dWdx,
			     const Chem2D_pState &dWdy,
			     const int Axisymmetric,
			     const Vector2D X,
			     const Vector2D &norm_dir) {

  Chem2D_cState U;
  double mu, kappa, div_v, radius;
  double Temperature, Rmix, Cp;
  double mu_t, kappa_t, Dm_t;
  Vector2D grad_T;
  Vector2D i = Vector2D(1,0), j = Vector2D(0,1);
  Tensor2D strain_rate;

  //Molecular transport properties
  mu = W.mu();
  kappa = W.kappa();
  Temperature = W.T();
  Rmix = W.Rtot();
  Cp = W.Cp();
 
  //Turbulence model transport properties
  if (W.flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
      W.flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
    mu_t = W.eddy_viscosity();
    kappa_t =  mu_t*Cp/W.Pr_turb();
    Dm_t = W.Dm_turb();
    
  
  } /* endif */

  //Conserved variables from primitive
  U.rho = W.rho;
  U.rhov = W.rhov();
  U.E = W.E();
  U.rhok = W.rho*W.k;
  U.rhoomega = W.rho*W.omega;
  //U.tau = W.tau;
  //U.lambda = W.lambda;

  for(int i=0; i<W.ns; i++){
    U.rhospec[i].c = W.rho*W.spec[i].c;
  } 
	
  //Strain rate (+dilatation)
  div_v = dWdx.v.x + dWdy.v.y;
  if (Axisymmetric == 2) {
    if(X.x < MICRO) {
      radius = MICRO;
    } else {
      radius = X.x;
    }
    div_v += W.v.x/radius;
    
  } else if (Axisymmetric == 1) {
    if(X.y < MICRO) {
      radius = MICRO;
    } else { 
       radius = X.y;
    }
    div_v += W.v.y/radius;
    
  } /* endif */
  
  strain_rate.xx = dWdx.v.x-div_v/THREE;
  strain_rate.xy = HALF*(dWdx.v.y + dWdy.v.x);
  strain_rate.yy = dWdy.v.y-div_v/THREE;
  
  if (Axisymmetric == 0) {
    strain_rate.zz = -(strain_rate.xx + strain_rate.yy); 
  } else if (Axisymmetric == 2) {
      strain_rate.zz = W.v.x/radius-div_v/THREE;
  } else if (Axisymmetric == 1) {
      strain_rate.zz = W.v.y/radius-div_v/THREE;
  } /* endif */
  
  //Temperature gradient
  //dT/dx = 1/rho*R *( dP/dx - P/rho * drho/dx)
  grad_T.x = (ONE/(W.rho*Rmix))*(dWdx.p - (W.p/W.rho)*dWdx.rho);
  grad_T.y = (ONE/(W.rho*Rmix))*(dWdy.p - (W.p/W.rho)*dWdy.rho);
 
  // Molecular (Laminar) diffusion of species
  // for each of the "n" species
  for( int k=0; k<U.ns; k++){
    /***************** Diffusion coefficients **********************/
    // using global Schmidt number relation Sc = mu/rho*Ds
    U.rhospec[k].diffusion_coef = mu/U.Schmidt[k];
    /***************** mass fraction gradients *********************/
    U.rhospec[k].gradc.x = U.rho * dWdx.spec[k].c;
    U.rhospec[k].gradc.y = U.rho * dWdy.spec[k].c;
  }
  
  //Molecular (laminar) stress tensor
  U.tau = TWO*mu*strain_rate;
  
  //Turbulent (Reynolds) Stresses
  if (W.flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
      W.flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
    U.lambda = (TWO*mu_t)*strain_rate;
    U.lambda.xx -= (TWO/THREE)*U.rhok; 
    U.lambda.yy -= (TWO/THREE)*U.rhok; 
    U.lambda.zz -= (TWO/THREE)*U.rhok; 
  } /* endif */
  //Molecular (laminar) heat flux
  //Thermal conduction, q = - kappa * grad(T)
  U.qflux = - kappa*grad_T;
  //Thermal diffusion, q -= rho * sum ( hs * Ds *gradcs)
  U.qflux -= U.rho*U.thermal_diffusion(Temperature);  
  
  //Turbulent heat flux
  //Thermal conduction, q = - kappa * grad(T)
  if (W.flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
      W.flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
    U.theta = - kappa_t*grad_T;
    //Thermal Diffusion, q -= rho * sum ( hs * Ds *gradcs)   
    for (int k=0; k<U.ns; k++) {
      U.theta -= Dm_t*U.rhospec[k].gradc*
	(U.specdata[k].Enthalpy(Temperature)+U.specdata[k].Heatofform());
    }
  } /* endif */
  //Xinfeng Note: the lambda = Reynolds_Stress(... ... ) may be computed/assigned outside
  //To see which way is more general ...
  // maybe modifying the viscous calculations some sort...
  //  lambda = Reynolds_Stress(dWdx,dWdy,Axisymmetric, X);
  //tau = Laminar_Stress(dWdx,dWdy,Axisymmetric, X);					       
  //Return Viscous Fluxes

  return (U.Viscous_Flux_x(dWdx)*dot(i,norm_dir) 
	  + U.Viscous_Flux_y(dWdy)*dot(j,norm_dir));

}

/**********************************************************************
 * Routine: WallShearStress                                           *
 *                                                                    *
 * This routine computes and returns the shear stress at a wall.      *
 *                                                                    *
 **********************************************************************/
double WallShearStress(const Chem2D_pState &W1,
		       const Vector2D &X1,
		       const Vector2D &X2,
		       const Vector2D &X3,
		       const Vector2D &norm_dir) {

  double l21, l32, l13, A, dWdn;
  Vector2D n21, n32, n13;
  Chem2D_pState W2, W3, W_face, dWdx, dWdy;

  // Initialze W2 and W3.
  W2.Vacuum(); W2.rho = W1.rho; W2.p = W1.p;
  W3.Vacuum(); W3.rho = W1.rho; W3.p = W1.p;
  for(int i=0; i<W1.ns; i++){
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
