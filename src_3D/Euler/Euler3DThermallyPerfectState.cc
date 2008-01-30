/*! \file  Euler3DThermallyPerfectState.cc
\brief  Definition of member functions for 3D Euler solution state classes
        for a thermally perfect non-reactive or combusting mixture
*/

/* Include Euler3DThermallyPerfectState header file. */

#ifndef _EULER3D_THERMALLYPERFECT_STATE_INCLUDED
#include "Euler3DThermallyPerfectState.h"
#endif // EULER3D_THERMALLYPERFECT_STATE_INCLUDED   

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState -- Create storage and assign various static values. *
 ***************************************************************************************/
int Euler3D_ThermallyPerfect_pState::ns = 1; 
int Euler3D_ThermallyPerfect_pState::num_vars = NUM_EULER3D_VAR_SANS_SPECIES; 
double* Euler3D_ThermallyPerfect_pState::_temp_values=NULL;
double* Euler3D_ThermallyPerfect_pState::_diff_coeff=NULL;
NASARP1311data* Euler3D_ThermallyPerfect_pState::specdata=NULL;
Reaction_set Euler3D_ThermallyPerfect_pState::React;
double Euler3D_ThermallyPerfect_pState::low_temp_range = 200.0;
double Euler3D_ThermallyPerfect_pState::high_temp_range = 300.0;
int Euler3D_ThermallyPerfect_pState::debug_level = 1; 
double* Euler3D_ThermallyPerfect_pState::Schmidt=NULL;

/***************************************************************************************
 * Euler3D_ThermallyPerfect_cState -- Create storage and assign various static values. *
 ***************************************************************************************/
int Euler3D_ThermallyPerfect_cState::ns = 1; 
int Euler3D_ThermallyPerfect_cState::num_vars = NUM_EULER3D_VAR_SANS_SPECIES;
double* Euler3D_ThermallyPerfect_cState::_temp_values=NULL;
double* Euler3D_ThermallyPerfect_cState::_diff_coeff=NULL;
NASARP1311data* Euler3D_ThermallyPerfect_cState::specdata=NULL;   
double Euler3D_ThermallyPerfect_cState::low_temp_range = 200.0;
double Euler3D_ThermallyPerfect_cState::high_temp_range = 300.0;
double* Euler3D_ThermallyPerfect_cState::Schmidt=NULL;
int Euler3D_ThermallyPerfect_cState::debug_level = 1; 

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState member functions                                    *
 ***************************************************************************************/

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::set_species_data -- Assign species data.           *
 ***************************************************************************************/
void Euler3D_ThermallyPerfect_pState::set_species_data(const int &n,
                                                       const string *S,
                                                       const char *PATH,
                                                       const int &debug, 
                                                       const double &Mr, 
                                                       const double* Sc,
                                                       const int &trans_data) { 

   // Set the new number of species. 
   ns = n;
   num_vars = ns + NUM_EULER3D_VAR_SANS_SPECIES;
   
   // Create static memory for NASA data.
   Deallocate_static();
   Allocate_static();

   // Create memory for species data (static or dynamic as requested).
   Deallocate();
   Allocate();

   // Read in and assign appropriate NASA data for each species to be used.
   for (int i = 0; i < ns; i++) {
     //overwrite default data  
     specdata[i].Getdata(S[i], PATH, trans_data);  
     Schmidt[i] = Sc[i];
   } /* endfor */

   // Set initial values for the species.
   for (int i=0; i < ns; i++) {
      spec[i].c = ONE/ns; 
   } /* endfor */

   // Set data temperature ranges for mixture calculations.
   Temp_low_range(); 
   Temp_high_range(); 
   
   // Set debug information level.
   debug_level = debug;   

}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::set_species_values -- Set initial values for       *
 *                                                        species mass fractions.      *
 ***************************************************************************************/
void Euler3D_ThermallyPerfect_pState::set_initial_values(void) {
   Deallocate();
   Allocate();
   for(int i = 0; i < ns; i++) {
      spec[i].c = ONE/ns ; 
   } /* endfor */
}

void Euler3D_ThermallyPerfect_pState::set_initial_values(const double &value) {
   Deallocate();
   Allocate();
   for (int i = 0; i < ns; i++) {
      spec[i].c = value ; 
   } /* endfor */
}

void Euler3D_ThermallyPerfect_pState::set_initial_values(const double *cfrac) {
  Deallocate();
  Allocate();
  for (int i = 0; i < ns; i++) {
    spec[i].c = cfrac[i];
  } /* endfor */
}

void Euler3D_ThermallyPerfect_pState::set_initial_values(const Species *mfrac) {
   Deallocate();
   Allocate();
   for (int i = 0; i < ns; i++) {
      spec[i] = mfrac[i];
   } /* endfor */
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::Copy -- Makes a copy of solution state vector.     *
 ***************************************************************************************/
void Euler3D_ThermallyPerfect_pState::Copy(const Euler3D_ThermallyPerfect_pState &W) {
  rho = W.rho;
  v = W.v; 
  p = W.p;  
  for (int i=0; i<ns; i++) {
    spec[i] = W.spec[i];
  } /* endfor */
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::Temp_low_range -- Set the lower bound for valid T. *
 ***************************************************************************************/
void Euler3D_ThermallyPerfect_pState::Temp_low_range(void) {  
   // Get max of the min temperature of the lowest region
   double temp = specdata[0].Low_range();
   for (int i = 0; i < ns; i++) {
      temp = max(specdata[i].Low_range(),temp);
   } /* endfor */
   low_temp_range = temp;  
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::Temp_high_range -- Set the upper bound for valid T.*
 ***************************************************************************************/
void Euler3D_ThermallyPerfect_pState::Temp_high_range(void) {
   // Get min of the max temperature of the highest region
   double temp = specdata[0].High_range();
   for (int i = 0; i < ns; i++) {
      temp = min(specdata[i].High_range(),temp);
   } /* endfor */
   high_temp_range = temp;  
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::negative_speccheck -- Check species mass fractions.*
 ***************************************************************************************/
bool Euler3D_ThermallyPerfect_pState::negative_speccheck(void) {
   double sum = ZERO;
   double LOCAL_TOL = MICRO; 
   //-------- Negative Check ------------//
   for (int i=0; i < ns-1; i++){
      if (spec[i].c > ONE) { //check for > 1.0 
         spec[i].c = ONE;
      } else if (spec[i].c < ZERO) { //check for -ve
         if(spec[i].c > -LOCAL_TOL) {  //check for small -ve and set to ZERO 
            spec[i].c = ZERO;
         } else {
            spec[i].c = ZERO;
         } /* endif */
      } /* endif */
      sum += spec[i].c;
   } /* endfor */
   spec[ns-1].c = max(ONE- sum, ZERO);
   sum += spec[ns-1].c;
   for(int i = 0; i < ns; i++) {
      spec[i].c = spec[i].c*(ONE/sum);
   } /* endfor */
   return true;
}

/*****************************************************************************************
 * Euler3D_ThermallyPerfect_pState::Realizable_Solution_Check -- Check physical validity *
 *                                                               validity of solution    *
 *                                                               state.                  *
 *****************************************************************************************/
bool Euler3D_ThermallyPerfect_pState::Realizable_Solution_Check(void) {
   if (rho <= ZERO || !negative_speccheck() || p <= ZERO) {
      cout << "\n " << CFFC_Name() 
           << " ERROR: Primitive solution state has a negative density, pressure,"
           << " and/or species mass fraction.\n";
      return false;
   } else {
      return true;
   } /* endif */
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::Mass -- Return mixture molecular mass (kg/mol).    *
 ***************************************************************************************/
double Euler3D_ThermallyPerfect_pState::Mass(void) const{
   // = 1 / sum(mass fration/mol_mass)
   double sum = ZERO;
   for (int i=0; i < ns; i++) {
      sum += spec[i].c/specdata[i].Mol_mass();
   } /* endfor */
   return ONE/sum;
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::Rtot -- Return mixture gas constant  J/(kg*K).     *
 ***************************************************************************************/
double Euler3D_ThermallyPerfect_pState::Rtot(void) {
   // = sum ( mass fraction * species gas constant)
   double sum = ZERO;
   for (int i = 0; i < ns; i++) {
      sum += spec[i].c * specdata[i].Rs();
   } /* endfor */
   return sum;
}

double Euler3D_ThermallyPerfect_pState::Rtot(void) const {
   // = sum ( mass fraction * species gas constant)
   double sum = 0.0;
   for (int i = 0; i < ns; i++) {
      sum += spec[i].c * specdata[i].Rs();
   } /* endfor */
   return sum;
}

/*****************************************************************************************
 * Euler3D_ThermallyPerfect_pState::HeatofFormation -- Return mixture heat of formation. *
 *****************************************************************************************/
double Euler3D_ThermallyPerfect_pState::HeatofFormation(void) {
   double sum = ZERO;
   for (int i = 0; i < ns; i++) {
      sum += spec[i].c*specdata[i].Heatofform();
   } /* endfor */
   return sum;
}

double Euler3D_ThermallyPerfect_pState::HeatofFormation(void) const {
   double sum = ZERO;
   for (int i = 0; i < ns; i++) {
      sum += spec[i].c*specdata[i].Heatofform();
   } /* endfor */
   return sum;
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_pState::Cp -- Return mixture heat capacity (const pressure) J/(kg*K). *
 **************************************************************************************************/
double Euler3D_ThermallyPerfect_pState::Cp(void) const {
   // = sum ( mass fraction * species Cp) 
   double Temp = T();
   double sum = ZERO;
   for (int i = 0; i < ns; i++) {
      sum += spec[i].c*specdata[i].HeatCapacity_p(Temp);
   } /* endfor */
   return sum;
}

double Euler3D_ThermallyPerfect_pState::Cp(const double& TEMP) const {
   // = sum ( mass fraction * species Cp) 
   double sum = 0.0;
   for (int i = 0; i < ns; i++) {
      sum += spec[i].c*specdata[i].HeatCapacity_p(TEMP);
   } /* endfor */
   return sum;
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_pState::Cv -- Return mixture heat capacity (const volume) J/(kg*K).   *
 **************************************************************************************************/
double Euler3D_ThermallyPerfect_pState::Cv(void) const {
   // = sum ( mass fraction * species Cv)  
   double Temp = T();
   double sum = ZERO;
   for (int i = 0; i < ns; i++) {
      sum += spec[i].c*specdata[i].HeatCapacity_v(Temp);
   } /* endfor */
   return sum;
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_pState::g -- Return mixture specific heat ratio.                      *
 **************************************************************************************************/
double Euler3D_ThermallyPerfect_pState::g(void) const{
   // = Cp / Cv  
   return Cp()/Cv();
}

/**************************************************
  polytropic heat ratio mixture gamma J/(kg*K)
  assuming T=200K as the temperature.
***************************************************/
double Euler3D_ThermallyPerfect_pState::gamma_guess(void) const{  
  double sum1 = ZERO;     double sum2 = ZERO;
  double gamma_s = ZERO;  double Temp = 200.0;
  for (int i = 0; i < ns; i++) {
     sum1 += spec[i].c*specdata[i].Rs();
     gamma_s = ONE/(ONE - specdata[i].Rs()/ specdata[i].HeatCapacity_p(Temp));
     sum2 += ((spec[i].c*specdata[i].Rs()) / (gamma_s - ONE)); 
  } /* endfor */
  return (ONE + sum1/sum2);
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_pState::e -- Return mixture absolute internal energy.                 *
 **************************************************************************************************/
// etotal = sensible & chemical
double Euler3D_ThermallyPerfect_pState::e(void) const{
   // = sum (mass fraction * species e) 
   double sum = ZERO;
   double Temp = T();
   for (int i = 0; i < ns; i++) { //(Enthalpy(Temp) - (R/mol_mass)*Temp)
      sum += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform() - specdata[i].Rs()*Temp);
   } /* endfor */
   return sum;
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_pState::es -- Return sensible internal energy.                        *
 **************************************************************************************************/
// internal energy with no heat of formation included (sensible)
double Euler3D_ThermallyPerfect_pState::es(void) const{
   // = sum (mass fraction * species e) 
   double sum = ZERO;
   double Temp = T();
   for (int i=0; i < ns; i++) { //(Enthalpy(Temp) - (R/mol_mass)*Temp)
      sum += spec[i].c*(specdata[i].Enthalpy(Temp) - specdata[i].Rs()*Temp);
   } /* endfor */
   return sum;
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_pState::eref -- Return reference internal energy.                     *
 **************************************************************************************************/
// reference internal energy e + heat of formation - offset
double Euler3D_ThermallyPerfect_pState::eref(void)const{
   // = sum (mass fraction * species e) 
   double sum = ZERO;
   double Temp = T();
   for (int i = 0; i < ns; i++) { //(Enthalpy(Temp) - (R/mol_mass)*Temp)
      sum += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform() - specdata[i].DeltaHref() - 
                        specdata[i].Rs()*Temp);
   } /* endfor */
   return sum;
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_pState::h -- Return mixture absolute internal enthalpy.               *
 **************************************************************************************************/
// Specific absolute enthalpy = mass fraction * (hsensible + heatofform)
double Euler3D_ThermallyPerfect_pState::h(void) const{
   // = sum (mass fraction * species h) 
   double sum = ZERO;  
   double Temp = T();
   for (int i = 0; i < ns; i++) {
      sum += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform());
   } /* endfor */
   return sum;
}

double Euler3D_ThermallyPerfect_pState::h(const double &Temp) const{
   // = sum (mass fraction * species h) 
   double sum = ZERO;  
   for (int i = 0; i < ns; i++) {
      sum += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform());
   } /* endfor */
   return sum;
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_pState::hs -- Return sensible internal enthalpy.                      *
 **************************************************************************************************/
double Euler3D_ThermallyPerfect_pState::hs(void) const{
   // = sum (mass fraction * species h) 
   double sum = 0.0;  
   double Temp = T();
   for (int i = 0; i < ns; i++) {
      sum += spec[i].c*(specdata[i].Enthalpy(Temp));
   } /* endfor */
   return sum;
}

double Euler3D_ThermallyPerfect_pState::hs(const double &Temp) const{
   // = sum (mass fraction * species h) 
   double sum = 0.0;  
   for (int i = 0; i < ns; i++) {
      sum += spec[i].c*(specdata[i].Enthalpy(Temp));
   } /* endfor */
   return sum;
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_pState::href -- Return reference internal enthalpy.                   *
 **************************************************************************************************/
// mass fraction * (hsensible + heat of formation - offset)
double Euler3D_ThermallyPerfect_pState::href(void)const{
   // = sum (mass fraction * species h) 
   double sum = 0.0;  
   double Temp = T();
   for (int i = 0; i < ns; i++) { 
      sum += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform() - specdata[i].DeltaHref());
   } /* endfor */
   return sum;
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_pState::hprime -- Return dhdT = Cp.                                   *
 **************************************************************************************************/
/**************************************************
   Derivative of specific enthalpy dh/dT
   actually is just Cp as Cp = (dh/dT)_p
***************************************************/
double Euler3D_ThermallyPerfect_pState::hprime(void) const{
   double sum = 0.0; 
   double Temp = T();
   for (int i = 0; i < ns; i++) {
      sum += spec[i].c*specdata[i].Enthalpy_prime(Temp);
   } /* endfor */
   return (sum);
}

double Euler3D_ThermallyPerfect_pState::hprime(double &Temp) const{
   double sum = 0.0;  
   for(int i = 0; i < ns; i++) {
      sum += spec[i].c*specdata[i].Enthalpy_prime(Temp);
   } /* endfor */
   return (sum);
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_pState::E -- Return total energy of the mixture.                      *
 **************************************************************************************************/
double Euler3D_ThermallyPerfect_pState::E(void) const{
   // E = rho*(e + velocity^2 )  
   return (rho*(e() + HALF*v.sqr())); 
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_pState::H -- Return total enthalpy of the mixture.                    *
 **************************************************************************************************/
double Euler3D_ThermallyPerfect_pState::H(void) const{
  // H = rho*(h + velocity^2)
  double TotalEnthalpy = ZERO;
  TotalEnthalpy = (rho*(h() + HALF*v.sqr()));
  return (TotalEnthalpy);
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_pState::Hs -- Return total sensible enthalpy of the mixture.          *
 **************************************************************************************************/
double Euler3D_ThermallyPerfect_pState::Hs(void) const{
  // Hs = rho*(hs + velocity^2)
  double H_S =  rho*(hs() + HALF*v.sqr() );
  return (H_S);
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::rhov -- Return mixture momentum vector.            *
 ***************************************************************************************/
Vector3D Euler3D_ThermallyPerfect_pState::rhov(void) const{
   return rho*v; 
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::T -- Return mixture temperature.                   *
 ***************************************************************************************/
double Euler3D_ThermallyPerfect_pState::T(void) const{
  return p/(rho*Rtot());
}

// Temperature derived from specific sensible enthalpy
double Euler3D_ThermallyPerfect_pState::T(double &h_s) const{
   double T = ZERO;
   //--------- Initial Guess ------------------------------//
   //using a polytropic gas assumption with gamma@200;
   double Tguess = (h_s/Rtot())*(ONE/(ONE/(gamma_guess() - ONE) +ONE));
   //--------- Newtons method to get T --------------------//
   int numit =0;
   double Tmin = low_temp_range;
   double Tmax = high_temp_range;
   //check for start value
   if (Tguess > Tmin && Tguess < Tmax) {
      T=Tguess;
   } else { 
      T=Tmin;
   } /* endif */

   double fa = hs(Tmin) - h_s;
   double fn = hs(T) - h_s; 
   double dfn = hprime(T);
   
   while (fabs(Tmax-Tmin) > CONV_TOLERANCE && 
          fabs(fn) > CONV_TOLERANCE && 
          numit <20 && 
          T >= low_temp_range) {  
      if (T >= Tmin && T <= Tmax) {
         T = T - fn/dfn;
         if (T >= Tmax) T = HALF*(Tmax - Tmin);	
         //Bisection
      } else {
         T = HALF*(Tmax - Tmin);
      } 
      fn = hs(T) - h_s;  
      dfn = hprime(T); 
      //change bisection range
      if (fa*fn <=ZERO) {
         Tmax = T;
      } else {
         Tmin = T;
         fa = fn;
      } /* endif */
      numit++;
   } /* endwhile */
   if (numit>=19 || T <= low_temp_range) {
      T = max(Tguess,low_temp_range); 
      if (debug_level) { 
         cout << "\nTemperature didn't converge in Euler3D_ThermallyPerfect_pState::T(void)";
         cout << " with polytopic Tguess " << Tguess <<", or lower than Tmin "
              << low_temp_range <<" using " << T;
      } /* endif */
   } /* endif */
   return T;
}
    
/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::a -- Return mixture sound speed.                   *
 ***************************************************************************************/
double Euler3D_ThermallyPerfect_pState::a(void){
   double a_temp;
   double RTOT= Rtot();
   a_temp  =  (g()*(Rtot()*T()));
   return sqrt(a_temp);
}

double Euler3D_ThermallyPerfect_pState::a(void) const{
   double a_temp;
   double RTOT= Rtot();
   a_temp = (g()*(Rtot()*T()));
   return sqrt(a_temp);
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::M -- Return mixture Mach number.                   *
 ***************************************************************************************/
double Euler3D_ThermallyPerfect_pState::M(void) const{
   return (v.abs()/a());
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::SpecCon -- Return ith species mass concentration.  *
 ***************************************************************************************/
double Euler3D_ThermallyPerfect_pState::SpecCon(int i) const {
  //returned in kg/m^3 / kg/mol => mol/m^3
  return (rho)*spec[i].c/(specdata[i].Mol_mass());
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::Gibbs -- Return Gibbs free energy of mixture.      *
 ***************************************************************************************/
/******* GIBBS Free Energy ********************************
  Eqn. (10.84) (10.87) Anderson
  Gs = Hs - TS
  Gs(ps=1) = Gs - R_UNIVERSAL*T*ln(ps)  //for data not at 1atm
  ps = cs(M/Ms)*p
***********************************************************/
double Euler3D_ThermallyPerfect_pState::Gibbs(int species) const{
   double Temp = T(); 
   return specdata[species].GibbsFree(Temp);
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::diedip -- Returns de/dp.                           *
 ***************************************************************************************/
double Euler3D_ThermallyPerfect_pState::diedip(void) const{
   double Rlocal = Rtot();
   return (hprime() - Rlocal)/(rho*Rlocal);
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::diedirho -- Returns de/drho.                       *
 ***************************************************************************************/
double Euler3D_ThermallyPerfect_pState::diedirho(void) const{
   double Rlocal = Rtot();
   return -p*(hprime() - Rlocal)/(rho*rho*Rlocal);
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::U -- Return conserved solution state vector.       *
 ***************************************************************************************/
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::U(void) {
   Euler3D_ThermallyPerfect_cState Temp;
   Temp.rho = rho;
   Temp.rhov = rhov();
   Temp.E = E();
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i] = rho*spec[i];
   } /* endfor */
   return Temp;
}

Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::U(void) const {
   Euler3D_ThermallyPerfect_cState Temp;
   Temp.rho = rho;
   Temp.rhov = rhov();
   Temp.E = E();
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i] = rho*spec[i];
   } /* endfor */
   return Temp;
}

Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::
U(const Euler3D_ThermallyPerfect_pState &W) {
   Euler3D_ThermallyPerfect_cState Temp;
   Temp.rho = W.rho;
   Temp.rhov = W.rhov();
   Temp.E = W.E();
   for (int i = 0; i < W.ns; i++) {
     Temp.rhospec[i] = W.rho*W.spec[i];
   } /* endfor */
   return Temp;
}

/***********************************************************************
 ***************** INVISCID FLUX VECTORS *******************************
 ***********************************************************************/

/************************************************************************
 * Euler3D_ThermallyPerfect_pState::F -- Inviscid flux (x-direction).   *
 ************************************************************************/
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::F(void) {
  Euler3D_ThermallyPerfect_cState Temp;
  Temp.rho = rho*v.x;
  Temp.rhov.x = rho*sqr(v.x) + p;
  Temp.rhov.y = rho*v.x*v.y;
  Temp.rhov.z = rho*v.x*v.z;
  Temp.E = v.x*H();
  for (int i = 0; i < ns; i++) {
     Temp.rhospec[i].c = rho*v.x*spec[i].c;
  } /* endfor */
  return (Temp);
}

Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::F(void) const {
  Euler3D_ThermallyPerfect_cState Temp;
  Temp.rho = rho*v.x;
  Temp.rhov.x = rho*sqr(v.x) + p;
  Temp.rhov.y = rho*v.x*v.y;
  Temp.rhov.z = rho*v.x*v.z;
  Temp.E = v.x*H();
  for (int i = 0; i < ns; i++) {
     Temp.rhospec[i].c = rho*v.x*spec[i].c;
  } /*endfor */
  return (Temp);
}

/************************************************************************
 * Euler3D_ThermallyPerfect_pState::Fx -- Inviscid flux (x-direction).  *
 ************************************************************************/
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::Fx(void) {
  Euler3D_ThermallyPerfect_cState Temp;
  Temp.rho = rho*v.x;
  Temp.rhov.x = rho*sqr(v.x) + p;
  Temp.rhov.y = rho*v.x*v.y;
  Temp.rhov.z = rho*v.x*v.z;
  Temp.E = v.x*H();
  for (int i = 0; i < ns; i++) {
     Temp.rhospec[i].c = rho*v.x*spec[i].c;
  } /* endfor */
  return (Temp);
}

Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::Fx(void) const {
  Euler3D_ThermallyPerfect_cState Temp;
  Temp.rho = rho*v.x;
  Temp.rhov.x = rho*sqr(v.x) + p;
  Temp.rhov.y = rho*v.x*v.y;
  Temp.rhov.z = rho*v.x*v.z;
  Temp.E = v.x*H();
  for (int i = 0; i < ns; i++) {
     Temp.rhospec[i].c = rho*v.x*spec[i].c;
  } /*endfor */
  return (Temp);
}

/************************************************************************
 * Euler3D_ThermallyPerfect_pState::Fy -- Inviscid flux (y-direction).  *
 ************************************************************************/
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::Fy(void) {
  Euler3D_ThermallyPerfect_cState Temp;
  Temp.rho = rho*v.y;
  Temp.rhov.x = rho*v.x*v.y;
  Temp.rhov.y = rho*sqr(v.y) + p; 
  Temp.rhov.z = rho*v.y*v.z;
  Temp.E = v.y*H();
  for (int i = 0; i < ns; i++) {
     Temp.rhospec[i].c = rho*v.y*spec[i].c;
  } /* endfor */
  return (Temp);
}

Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::Fy(void) const {
  Euler3D_ThermallyPerfect_cState Temp;
  Temp.rho = rho*v.y;
  Temp.rhov.x = rho*v.x*v.y;
  Temp.rhov.y = rho*sqr(v.y) + p; 
  Temp.rhov.z = rho*v.y*v.z;
  Temp.E = v.y*H();
  for (int i = 0; i < ns; i++) {
     Temp.rhospec[i].c = rho*v.y*spec[i].c;
  } /*endfor */
  return (Temp);
}

/************************************************************************
 * Euler3D_ThermallyPerfect_pState::Fz -- Inviscid flux (z-direction).  *
 ************************************************************************/
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::Fz(void) {
  Euler3D_ThermallyPerfect_cState Temp;
  Temp.rho = rho*v.z;
  Temp.rhov.x = rho*v.x*v.z;
  Temp.rhov.y = rho*v.y*v.z;
  Temp.rhov.z = rho*sqr(v.z) + p;
  Temp.E = v.z*H();
  for (int i = 0; i < ns; i++) {
     Temp.rhospec[i].c = rho*v.z*spec[i].c;
  } /* endfor */
  return (Temp);
}

Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::Fz(void) const {
  Euler3D_ThermallyPerfect_cState Temp;
  Temp.rho = rho*v.z;
  Temp.rhov.x = rho*v.x*v.z;
  Temp.rhov.y = rho*v.y*v.z;
  Temp.rhov.z = rho*sqr(v.z) + p;
  Temp.E = v.z*H();
  for (int i = 0; i < ns; i++) {
     Temp.rhospec[i].c = rho*v.z*spec[i].c;
  } /*endfor */
  return (Temp);
}

/*******************************************************************
 ***************** INVISCID FLUX JACOBIANS *************************
 *******************************************************************/

/***********************************************************************************
 * Euler3D_ThermallyPerfect_pState::dFxdU -- Invisicd flux Jacobian (x-direction). *     
 ***********************************************************************************/
void Euler3D_ThermallyPerfect_pState::dFxdU(DenseMatrix &dFdU){

  // SHOULD DOUBLE CHECK THIS !!!
  double Temp = T();
  double Rt = Rtot();
  double C_p = Cp();
  double ht = h();
  double denominator = (C_p/Rt - ONE);
  double phi = ZERO; 

  for(int i=0; i<ns-1; i++){ 
    phi += spec[i].c*(specdata[i].Enthalpy(Temp) + 
                      specdata[i].Heatofform() - 
                      C_p*Temp*specdata[i].Rs()/Rt -
                      (specdata[ns-1].Enthalpy(Temp)+specdata[ns-1].Heatofform() -
                      C_p*Temp*specdata[ns-1].Rs()/Rt));   
  }
  dFdU(0,1) += ONE;

  dFdU(1,0) += ( (C_p/Rt)*( - v.x*v.x) + HALF*(THREE*v.x*v.x + v.y+v.y) - 
               ht + C_p*Temp + phi )/denominator;
  dFdU(1,1) += v.x*(TWO*C_p/Rt-THREE)/denominator; 
  dFdU(1,2) -= v.y/denominator;
  dFdU(1,3) -= v.z/denominator;  
  dFdU(1,4) += ONE/denominator;

  dFdU(2,0) -= v.x*v.y;
  dFdU(2,1) += v.y;
  dFdU(2,2) += v.x;

  dFdU(3,0) -= v.x*v.y;
  dFdU(3,1) += v.z;
  dFdU(3,3) += v.x;

  dFdU(4,0) += v.x*( v.x*v.x + v.y+v.y +v.z*v.z + C_p*Temp - 
               (C_p/Rt)*( HALF*(v.x*v.x + v.y+v.y + v.z*v.z) + ht) + phi)/denominator;

  dFdU(4,1) += ht + HALF*(v.x*v.x + v.y*v.y) - v.x*v.x/denominator;
  dFdU(4,2) -= v.x*v.y/denominator;
  dFdU(4,3) -= v.x*v.z/denominator;
  dFdU(4,4) += v.x*C_p/denominator/Rt;

  //Species
  int nv_san_spec = NumVarSansSpecies();

  for(int i = 0; i<ns-1; i++){ 
    dFdU(1,nv_san_spec+i) -= (specdata[i].Enthalpy(Temp) + specdata[i].Heatofform() - 
                              C_p*Temp*specdata[i].Rs()/Rt -
			      (specdata[ns-1].Enthalpy(Temp)+specdata[ns-1].Heatofform() - 
                              C_p*Temp*specdata[ns-1].Rs()/Rt))/denominator; 
    dFdU(3,nv_san_spec+i) = v.x*dFdU(1,nv_san_spec+i);    
    dFdU(nv_san_spec+i, 0) -= spec[i].c*v.x ;
    dFdU(nv_san_spec+i, 1) += spec[i].c ;
    dFdU(nv_san_spec+i,nv_san_spec+i) += v.x ;        
  }

}

/*******************************************************************
 ***************** EIGENVALUES *************************************
 *******************************************************************/

/******************************************************************************
 * Euler3D_ThermallyPerfect_pState::lambda -- Eigenvalue(s) (x-direction).    *
 ******************************************************************************/
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::lambda(void) {
  double cc = a();
  Euler3D_ThermallyPerfect_pState Temp;
  Temp.rho = v.x - cc;
  Temp.v.x = v.x;
  Temp.v.y = v.x;
  Temp.v.z = v.x;
  Temp.p = v.x + cc;
  for (int i = 0; i < ns; i++) {
     Temp.spec[i].c = v.x;
  } /* endfor */
  return (Temp);
}

Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::lambda(void) const {
  double cc = a();
  Euler3D_ThermallyPerfect_pState Temp;
  Temp.rho = v.x - cc;
  Temp.v.x = v.x;
  Temp.v.y = v.x;
  Temp.v.z = v.x;
  Temp.p = v.x + cc;
  for (int i = 0; i < ns; i++) {
     Temp.spec[i].c = v.x;
  } /* endfor */
  return (Temp);
}

/******************************************************************************
 * Euler3D_ThermallyPerfect_pState::lambda_x -- Eigenvalue(s) (x-direction).  *
 ******************************************************************************/
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::lambda_x(void) {
  double cc = a();
  Euler3D_ThermallyPerfect_pState Temp;
  Temp.rho = v.x - cc;
  Temp.v.x = v.x;
  Temp.v.y = v.x;
  Temp.v.z = v.x;
  Temp.p = v.x + cc;
  for (int i = 0; i < ns; i++) {
     Temp.spec[i].c = v.x;
  } /* endfor */
  return (Temp);
}

Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::lambda_x(void) const {
  double cc = a();
  Euler3D_ThermallyPerfect_pState Temp;
  Temp.rho = v.x - cc;
  Temp.v.x = v.x;
  Temp.v.y = v.x;
  Temp.v.z = v.x;
  Temp.p = v.x + cc;
  for (int i = 0; i < ns; i++) {
     Temp.spec[i].c = v.x;
  } /* endfor */
  return (Temp);
}

/********************************************************************
 ***************** EIGENVECTORS *************************************
 ********************************************************************/

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::rc -- Conserved right eigenvector (x-direction).   *
 ***************************************************************************************/
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::rc(const int &index) {
   double cc;
   switch(index){  
   case 1: 
      cc = a();
      return (Euler3D_ThermallyPerfect_cState(ONE, 
                                              v.x-cc, v.y, v.z, 
                                              H()/rho-v.x*cc, 
                                              spec));
   case 2:
      return (Euler3D_ThermallyPerfect_cState(ONE, 
                                              v.x, v.y, v.z, 
                                              H()/rho-Cp()*T(), 
                                              spec));
   case 3:
      return (Euler3D_ThermallyPerfect_cState(ZERO, 
                                              ZERO, rho, ZERO, 
                                              rho*v.y, 
                                              ZERO));
   case 4: 
      return (Euler3D_ThermallyPerfect_cState(ZERO, 
                                              ZERO, ZERO, rho, 
                                              rho*v.z, 
                                              ZERO));
   case 5: 
      cc = a();
      return (Euler3D_ThermallyPerfect_cState(ONE, 
                                              v.x+cc, v.y, v.z, 
                                              H()/rho+v.x*cc, 
                                              spec));
   default: 
      Euler3D_ThermallyPerfect_cState rr(ZERO);
      double RTOT = Rtot();
      double TEMP = p/(rho*RTOT);      
      rr.E = rho*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP) + 
              specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Heatofform() - 
              Cp(TEMP)*TEMP*specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs()/RTOT);
      double PHI = specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP)-
                   specdata[num_vars- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP)-
                   (Cp(TEMP) -RTOT)*TEMP*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs() - 
                   specdata[num_vars - (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs())/RTOT;
      rr.rhospec[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].c = rho*PHI; 
      return rr;
   };
}

Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::rc(const int &index) const {
   double cc;
   switch(index){  
   case 1: 
      cc = a();
      return (Euler3D_ThermallyPerfect_cState(ONE, 
                                              v.x-cc, v.y, v.z, 
                                              H()/rho-v.x*cc, 
                                              spec));
   case 2:
      return (Euler3D_ThermallyPerfect_cState(ONE, 
                                              v.x, v.y, v.z, 
                                              H()/rho-Cp()*T(), 
                                              spec));
   case 3:
      return (Euler3D_ThermallyPerfect_cState(ZERO, 
                                              ZERO, rho, ZERO, 
                                              rho*v.y, 
                                              ZERO));
   case 4: 
      return (Euler3D_ThermallyPerfect_cState(ZERO, 
                                              ZERO, ZERO, rho, 
                                              rho*v.z, 
                                              ZERO));
   case 5: 
      cc = a();
      return (Euler3D_ThermallyPerfect_cState(ONE, 
                                              v.x+cc, v.y, v.z, 
                                              H()/rho+v.x*cc, 
                                              spec));
   default: 
      Euler3D_ThermallyPerfect_cState rr(ZERO);
      double RTOT = Rtot();
      double TEMP = p/(rho*RTOT);      
      rr.E = rho*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP) + 
              specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Heatofform() - 
              Cp(TEMP)*TEMP*specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs()/RTOT);
      double PHI = specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP)-
                   specdata[num_vars- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP)-
                   (Cp(TEMP) -RTOT)*TEMP*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs() - 
                   specdata[num_vars - (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs())/RTOT;
      rr.rhospec[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].c = rho*PHI; 
      return rr;
   };
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::rc_x -- Conserved right eigenvector (x-direction). *
 ***************************************************************************************/
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::rc_x(const int &index) {
   double cc;
   switch(index){  
   case 1: 
      cc = a();
      return (Euler3D_ThermallyPerfect_cState(ONE, 
                                              v.x-cc, v.y, v.z, 
                                              H()/rho-v.x*cc, 
                                              spec));
   case 2:
      return (Euler3D_ThermallyPerfect_cState(ONE, 
                                              v.x, v.y, v.z, 
                                              H()/rho-Cp()*T(), 
                                              spec));
   case 3:
      return (Euler3D_ThermallyPerfect_cState(ZERO, 
                                              ZERO, rho, ZERO, 
                                              rho*v.y, 
                                              ZERO));
   case 4: 
      return (Euler3D_ThermallyPerfect_cState(ZERO, 
                                              ZERO, ZERO, rho, 
                                              rho*v.z, 
                                              ZERO));
   case 5: 
      cc = a();
      return (Euler3D_ThermallyPerfect_cState(ONE, 
                                              v.x+cc, v.y, v.z, 
                                              H()/rho+v.x*cc, 
                                              spec));
   default: 
      Euler3D_ThermallyPerfect_cState rr(ZERO);
      double RTOT = Rtot();
      double TEMP = p/(rho*RTOT);      
      rr.E = rho*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP) + 
              specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Heatofform() - 
              Cp(TEMP)*TEMP*specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs()/RTOT);
      double PHI = specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP)-
                   specdata[num_vars- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP)-
                   (Cp(TEMP) -RTOT)*TEMP*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs() - 
                   specdata[num_vars - (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs())/RTOT;
      rr.rhospec[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].c = rho*PHI; 
      return rr;
   };
}

Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::rc_x(const int &index) const {
   double cc;
   switch(index){  
   case 1: 
      cc = a();
      return (Euler3D_ThermallyPerfect_cState(ONE, 
                                              v.x-cc, v.y, v.z, 
                                              H()/rho-v.x*cc, 
                                              spec));
   case 2:
      return (Euler3D_ThermallyPerfect_cState(ONE, 
                                              v.x, v.y, v.z, 
                                              H()/rho-Cp()*T(), 
                                              spec));
   case 3:
      return (Euler3D_ThermallyPerfect_cState(ZERO, 
                                              ZERO, rho, ZERO, 
                                              rho*v.y, 
                                              ZERO));
   case 4: 
      return (Euler3D_ThermallyPerfect_cState(ZERO, 
                                              ZERO, ZERO, rho, 
                                              rho*v.z, 
                                              ZERO));
   case 5: 
      cc = a();
      return (Euler3D_ThermallyPerfect_cState(ONE, 
                                              v.x+cc, v.y, v.z, 
                                              H()/rho+v.x*cc, 
                                              spec));
   default: 
      Euler3D_ThermallyPerfect_cState rr(ZERO);
      double RTOT = Rtot();
      double TEMP = p/(rho*RTOT);      
      rr.E = rho*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP) + 
              specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Heatofform() - 
              Cp(TEMP)*TEMP*specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs()/RTOT);
      double PHI = specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP)-
                   specdata[num_vars- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP)-
                   (Cp(TEMP) -RTOT)*TEMP*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs() - 
                   specdata[num_vars - (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs())/RTOT;
      rr.rhospec[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].c = rho*PHI; 
      return rr;
   };
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::lp -- Primitive left eigenvector (x-direction).    *
 ***************************************************************************************/
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::lp(const int &index) {
   double cc;
   switch(index){  
   case 1 :
      cc = a();
      return (Euler3D_ThermallyPerfect_pState(ZERO, 
                                              -HALF*rho/cc, ZERO, ZERO, 
                                              HALF/(cc*cc), 
                                              ZERO));
   case 2:
      cc = a();
      return (Euler3D_ThermallyPerfect_pState(ONE, 
                                              ZERO, ZERO, ZERO, 
                                              -ONE/(cc*cc), 
                                              ZERO));
   case 3:
      return (Euler3D_ThermallyPerfect_pState(ZERO, 
                                              ZERO, ONE, ZERO, 
                                              ZERO, 
                                              ZERO));
   case 4 :
      return (Euler3D_ThermallyPerfect_pState(ZERO, 
                                              ZERO, ZERO, ONE, 
                                              ZERO, 
                                              ZERO));
   case 5 :
      cc = a();
      return (Euler3D_ThermallyPerfect_pState(ZERO, 
                                              HALF*rho/cc, ZERO, ZERO, 
                                              HALF/(cc*cc), 
                                              ZERO));
   default:
      Euler3D_ThermallyPerfect_pState ll(ZERO);
      ll.spec[index-(NUM_EULER3D_VAR_SANS_SPECIES+1)].c = ONE;
      return ll;
   };
}

// Primitive Left Eigenvector -- (x-direction)
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::lp(const int &index) const {
   double cc;
   switch(index){  
   case 1 :
      cc = a();
      return (Euler3D_ThermallyPerfect_pState(ZERO, 
                                              -HALF*rho/cc, ZERO, ZERO, 
                                              HALF/(cc*cc), 
                                              ZERO));
   case 2:
      cc = a();
      return (Euler3D_ThermallyPerfect_pState(ONE, 
                                              ZERO, ZERO, ZERO, 
                                              -ONE/(cc*cc), 
                                              ZERO));
   case 3:
      return (Euler3D_ThermallyPerfect_pState(ZERO, 
                                              ZERO, ONE, ZERO, 
                                              ZERO, 
                                              ZERO));
   case 4 :
      return (Euler3D_ThermallyPerfect_pState(ZERO, 
                                              ZERO, ZERO, ONE, 
                                              ZERO, 
                                              ZERO));
   case 5 :
      cc = a();
      return (Euler3D_ThermallyPerfect_pState(ZERO, 
                                              HALF*rho/cc, ZERO, ZERO, 
                                              HALF/(cc*cc), 
                                              ZERO));
   default:
      Euler3D_ThermallyPerfect_pState ll(ZERO);
      ll.spec[index-(NUM_EULER3D_VAR_SANS_SPECIES+1)].c = ONE;
      return ll;
   };
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::lp_x -- Primitive left eigenvector (x-direction).  *
 ***************************************************************************************/
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::lp_x(const int &index) {
   double cc;
   switch(index){  
   case 1 :
      cc = a();
      return (Euler3D_ThermallyPerfect_pState(ZERO, 
                                              -HALF*rho/cc, ZERO, ZERO, 
                                              HALF/(cc*cc), 
                                              ZERO));
   case 2:
      cc = a();
      return (Euler3D_ThermallyPerfect_pState(ONE, 
                                              ZERO, ZERO, ZERO, 
                                              -ONE/(cc*cc), 
                                              ZERO));
   case 3:
      return (Euler3D_ThermallyPerfect_pState(ZERO, 
                                              ZERO, ONE, ZERO, 
                                              ZERO, 
                                              ZERO));
   case 4 :
      return (Euler3D_ThermallyPerfect_pState(ZERO, 
                                              ZERO, ZERO, ONE, 
                                              ZERO, 
                                              ZERO));
   case 5 :
      cc = a();
      return (Euler3D_ThermallyPerfect_pState(ZERO, 
                                              HALF*rho/cc, ZERO, ZERO, 
                                              HALF/(cc*cc), 
                                              ZERO));
   default:
      Euler3D_ThermallyPerfect_pState ll(ZERO);
      ll.spec[index-(NUM_EULER3D_VAR_SANS_SPECIES+1)].c = ONE;
      return ll;
   };
}

// Primitive Left Eigenvector -- (x-direction)
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::lp_x(const int &index) const {
   double cc;
   switch(index){  
   case 1 :
      cc = a();
      return (Euler3D_ThermallyPerfect_pState(ZERO, 
                                              -HALF*rho/cc, ZERO, ZERO, 
                                              HALF/(cc*cc), 
                                              ZERO));
   case 2:
      cc = a();
      return (Euler3D_ThermallyPerfect_pState(ONE, 
                                              ZERO, ZERO, ZERO, 
                                              -ONE/(cc*cc), 
                                              ZERO));
   case 3:
      return (Euler3D_ThermallyPerfect_pState(ZERO, 
                                              ZERO, ONE, ZERO, 
                                              ZERO, 
                                              ZERO));
   case 4 :
      return (Euler3D_ThermallyPerfect_pState(ZERO, 
                                              ZERO, ZERO, ONE, 
                                              ZERO, 
                                              ZERO));
   case 5 :
      cc = a();
      return (Euler3D_ThermallyPerfect_pState(ZERO, 
                                              HALF*rho/cc, ZERO, ZERO, 
                                              HALF/(cc*cc), 
                                              ZERO));
   default:
      Euler3D_ThermallyPerfect_pState ll(ZERO);
      ll.spec[index-(NUM_EULER3D_VAR_SANS_SPECIES+1)].c = ONE;
      return ll;
   };
}


/************************************************************************
 *************** NUMERICAL FLUX FUNCTIONS *******************************
 ************************************************************************/

/******************************************************************************************
 * Euler3D_ThermallyPerfect_pState::RoeAverage -- Roe-averaged primitive solution state.  *
 ******************************************************************************************/
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::
RoeAverage(const Euler3D_ThermallyPerfect_pState &Wl,
           const Euler3D_ThermallyPerfect_pState &Wr) {

   double Hl, Hr, srhol, srhor;
   double Ha, ha;
   Euler3D_ThermallyPerfect_pState Temp;

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
   Temp.v.z = (srhol*Wl.v.z+srhor*Wr.v.z)/(srhol+srhor);
   for (int i = 0; i < Wl.ns; i++) {
      Temp.spec[i].c = (srhol*Wl.spec[i].c + srhor*Wr.spec[i].c)/(srhol+srhor);
   } /* endif */

   Ha = (srhol*Hl+srhor*Hr)/(srhol+srhor);
   ha = Ha - HALF*(sqr(Temp.v.x)+sqr(Temp.v.y)+sqr(Temp.v.z));
   
   double TEMP = Temp.T(ha);
   Temp.p = Temp.rho*TEMP*Temp.Rtot();

   /* Return the Roe-averged state. */

   return (Temp);     
      
}

/******************************************************************************************
 * Euler3D_ThermallyPerfect_pState::FluxHLLE_x -- HLLE flux function, x-direction flux.   *
 ******************************************************************************************/
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::
FluxHLLE_x(const Euler3D_ThermallyPerfect_pState &Wl,
           const Euler3D_ThermallyPerfect_pState &Wr) {
   
   double wavespeed_l, wavespeed_r;
   Euler3D_ThermallyPerfect_pState Wa, lambdas_l, lambdas_r, lambdas_a;
   Euler3D_ThermallyPerfect_cState Flux, dUrl;
   
   /* Evaluate the Roe-average primitive solution state. */
   
   Wa = RoeAverage(Wl, Wr);

   /* Evaluate the jumps in the conserved solution states. */

   dUrl = Wr.U() - Wl.U();

   /* Evaluate the left, right, and average state eigenvalues. */
   
   lambdas_l = Wl.lambda_x();
   lambdas_r = Wr.lambda_x();
   lambdas_a = Wa.lambda_x();
    
   /* Determine the intermediate state flux. */
   
   wavespeed_l = min(lambdas_l[1],
                     lambdas_a[1]);
   wavespeed_r = max(lambdas_r[5],
                     lambdas_a[5]);

   wavespeed_l = min(wavespeed_l, ZERO);
   wavespeed_r = max(wavespeed_r, ZERO);
   
   if (wavespeed_l >= ZERO) {
      Flux = Wl.F();
   } else if (wavespeed_r <= ZERO) {
      Flux = Wr.F();
   } else {
      Flux = ( (wavespeed_r*Wl.F()-wavespeed_l*Wr.F())
              +(wavespeed_l*wavespeed_r)*dUrl)/
             (wavespeed_r-wavespeed_l);
   } /* endif */

   /* Return solution flux. */
   
   return (Flux);
   
}

Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::
FluxHLLE_x(const Euler3D_ThermallyPerfect_cState &Ul,
           const Euler3D_ThermallyPerfect_cState &Ur) {

   return (FluxHLLE_x(Ul.W(), Ur.W()));

}

/******************************************************************************************
 * Euler3D_ThermallyPerfect_pState::FluxHLLE_n -- HLLE flux function, n-direction flux.   *
 ******************************************************************************************/
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::
FluxHLLE_n(const Euler3D_ThermallyPerfect_pState &Wl,
           const Euler3D_ThermallyPerfect_pState &Wr,
           const Vector3D &norm_dir) {

   // Determine the left and right solution states in the rotate frame.
   Euler3D_ThermallyPerfect_pState Wl_rot(Wl.Rotate(norm_dir));
   Euler3D_ThermallyPerfect_pState Wr_rot(Wr.Rotate(norm_dir));
    
   // Evaluate the intermediate state solution flux in the rotated frame.
   Euler3D_ThermallyPerfect_cState Flux_rot = FluxHLLE_x(Wl_rot, Wr_rot);
 
   // Return numerical flux in un-rotated frame.
   return (Flux_rot.RotateBack(norm_dir));

}

Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::
FluxHLLE_n(const Euler3D_ThermallyPerfect_cState &Ul,
           const Euler3D_ThermallyPerfect_cState &Ur,
           const Vector3D &norm_dir) {

   return (FluxHLLE_n(Ul.W(), Ur.W(), norm_dir));

}

/******************************************************************************************
 * Euler3D_ThermallyPerfect_pState::FluxRoe_x -- Roe flux function, x-direction flux.     *
 ******************************************************************************************/
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::
FluxRoe_x(const  Euler3D_ThermallyPerfect_pState &Wl,  
          const  Euler3D_ThermallyPerfect_pState &Wr) {
   
   Euler3D_ThermallyPerfect_pState Wa, dWrl, wavespeeds, 
                                   lambdas_l, lambdas_r, lambdas_a;
   Euler3D_ThermallyPerfect_cState Flux;
   
   /* Evaluate the Roe-average primitive solution state. */    

   Wa = RoeAverage(Wl, Wr);
   
   /* Evaluate the jumps in the primitive solution states. */

   dWrl = Wr-Wl;
   
   /* Evaluate the left, right, and average state eigenvalues. */

   lambdas_l = Wl.lambda_x();
   lambdas_r = Wr.lambda_x();
   lambdas_a = Wa.lambda_x();
       
   /* Determine the intermediate state flux. */

   if (Wa.v.x >= ZERO) {
      Flux = Wl.F();   
      wavespeeds = Wa.lambda_minus(lambdas_a,
                                   lambdas_l,
                                   lambdas_r);
       
      for (int i=1 ; i <= Wl.num_vars; i++) {
         if (wavespeeds[i] < ZERO) {
            Flux += wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
         } /* endif */
      } /* endfor */

   } else {
      Flux = Wr.F();
      wavespeeds = Wa.lambda_plus(lambdas_a,
                                  lambdas_l,
                                  lambdas_r);
      
      for (int i=1; i <= Wl.num_vars; i++) {
         if (wavespeeds[i] > ZERO) {
            Flux -= wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i); 
         } /* endif */
      } /* endfor */
   } /* endif */
    
   /* Return solution flux. */    

   return (Flux);    
   
}
   
/******************************************************************************************
 * Euler3D_ThermallyPerfect_pState::FluxRoe_n -- Roe flux function, n-direction flux.     *
 ******************************************************************************************/
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::
FluxRoe_n(const Euler3D_ThermallyPerfect_pState &Wl,
          const Euler3D_ThermallyPerfect_pState &Wr,
          const Vector3D &norm_dir) {

   // Determine the left and right solution states in the rotate frame.
   Euler3D_ThermallyPerfect_pState Wl_rot(Wl.Rotate(norm_dir));
   Euler3D_ThermallyPerfect_pState Wr_rot(Wr.Rotate(norm_dir));
    
   // Evaluate the intermediate state solution flux in the rotated frame.
   Euler3D_ThermallyPerfect_cState Flux_rot = FluxRoe_x(Wl_rot, Wr_rot);
 
   // Return numerical flux in un-rotated frame.
   return (Flux_rot.RotateBack(norm_dir));

}

/******************************************************************************************
 * Euler3D_ThermallyPerfect_pState::lambda_minus -- Negative wave speeds determined using *
 *                                                  Harten entropy fix.                   *
 ******************************************************************************************/
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::
lambda_minus(const Euler3D_ThermallyPerfect_pState &lambdas_a,
             const Euler3D_ThermallyPerfect_pState &lambdas_l,
             const Euler3D_ThermallyPerfect_pState &lambdas_r) {
  
   Euler3D_ThermallyPerfect_pState W_temp;
   
   W_temp.rho = HartenFixNeg(lambdas_a[1], lambdas_l[1], lambdas_r[1]);
   W_temp.v.x = HALF*(lambdas_a[2]-fabs(lambdas_a[2]));
   W_temp.v.y = HALF*(lambdas_a[3]-fabs(lambdas_a[3]));
   W_temp.v.z = HALF*(lambdas_a[4]-fabs(lambdas_a[4]));
   W_temp.p = HartenFixNeg(lambdas_a[5], lambdas_l[5], lambdas_r[5]);
   
   for (int i = NUM_EULER3D_VAR_SANS_SPECIES+1; i <= W_temp.num_vars; i++) {
      W_temp.spec[i-(NUM_EULER3D_VAR_SANS_SPECIES+1)].c = 
          HALF*(lambdas_a[i]-fabs(lambdas_a[i]));
   } /* endfor */
   
   return (W_temp);

}

/******************************************************************************************
 * Euler3D_ThermallyPerfect_pState::lambda_plus -- Positive wave speeds determined using  *
 *                                                 Harten entropy fix.                    *
 ******************************************************************************************/
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::
lambda_plus(const Euler3D_ThermallyPerfect_pState &lambdas_a,
            const Euler3D_ThermallyPerfect_pState &lambdas_l,
            const Euler3D_ThermallyPerfect_pState &lambdas_r) {
   
   Euler3D_ThermallyPerfect_pState W_temp;

   W_temp.rho = HartenFixPos(lambdas_a[1], lambdas_l[1], lambdas_r[1]);
   W_temp.v.x = HALF*(lambdas_a[2]+fabs(lambdas_a[2]));
   W_temp.v.y = HALF*(lambdas_a[3]+fabs(lambdas_a[3]));
   W_temp.v.z = HALF*(lambdas_a[4]+fabs(lambdas_a[4]));
   W_temp.p = HartenFixPos(lambdas_a[5], lambdas_l[5], lambdas_r[5]);
  
   for(int i = NUM_EULER3D_VAR_SANS_SPECIES+1; i <= W_temp.num_vars; i++) {
      W_temp.spec[i-(NUM_EULER3D_VAR_SANS_SPECIES+1)].c = 
          HALF*(lambdas_a[i]+fabs(lambdas_a[i]));
   } /* endif */

   return (W_temp);

}

/******************************************************************************************
 * Euler3D_ThermallyPerfect_pState::HLLE_wavespeeds -- Returns the lambda plus and lambda *
 *                    minus wave speeds for rotated Riemann problem aligned with norm_dir *    
 *                    given unroated solution states Wl and Wr.                           *
 ******************************************************************************************/
Vector2D Euler3D_ThermallyPerfect_pState::HLLE_wavespeeds(const Euler3D_ThermallyPerfect_pState &Wl,
							  const Euler3D_ThermallyPerfect_pState &Wr,
							  const Vector3D &norm_dir) {

    Vector2D wavespeed;
    Euler3D_ThermallyPerfect_pState Wa_n, lambdas_l, lambdas_r, lambdas_a;  //Lots of TEMPS
    Euler3D_ThermallyPerfect_pState Wl_rotated(Wl.Rotate(norm_dir));
    Euler3D_ThermallyPerfect_pState Wr_rotated(Wr.Rotate(norm_dir));

    /* Evaluate the Roe-average primitive solution state. */                           
    Wa_n = Wa_n.RoeAverage(Wl_rotated, Wr_rotated);
    
    /* Evaluate the left, right, and average state eigenvalues. */
    lambdas_l = Wl_rotated.lambda_x();
    lambdas_r = Wr_rotated.lambda_x();
    lambdas_a = Wa_n.lambda_x();

    /* Determine the intermediate state flux. */
    wavespeed.x = min(lambdas_l[1],
                      lambdas_a[1]);
    wavespeed.y = max(lambdas_r[5],
                      lambdas_a[5]);
 
    wavespeed.x = min(wavespeed.x, ZERO); //lambda minus
    wavespeed.y = max(wavespeed.y, ZERO); //lambda plus 

    return (wavespeed);

}

/******************************************************************************************
 * Euler3D_ThermallyPerfect_pState::Rotate -- Returns a rotated primitive state aligned   *
 *                                            with a local x-axis in the norm_dir.        *
 ******************************************************************************************/
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::Rotate(const Vector3D &norm_dir) const {

   double cos_psi, sin_psi, cos_phi, sin_phi, cos_theta, sin_theta;
   cos_phi = ONE;
   sin_phi = ZERO;
   if (fabs(fabs(norm_dir.x)-ONE) < TOLER) {
      cos_psi = norm_dir.x/fabs(norm_dir.x);
      sin_psi = ZERO;
      cos_theta = ONE;
      sin_theta = ZERO;
   } else {
      cos_psi = norm_dir.x;
      sin_psi = sqrt(norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z);
      cos_theta = norm_dir.y/sqrt(norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z);
      sin_theta = norm_dir.z/sqrt(norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z);
   } /* endif */

   return Euler3D_ThermallyPerfect_pState(rho,
                                          (cos_psi*cos_phi-cos_theta*sin_phi*sin_psi)*v.x + 
                                          (cos_psi*sin_phi+cos_theta*cos_phi*sin_psi)*v.y + 
                                          (sin_psi*sin_theta)*v.z,
                                          (-sin_psi*cos_phi-cos_theta*sin_phi*cos_psi)*v.x + 
                                          (-sin_psi*sin_phi+cos_theta*cos_phi*cos_psi)*v.y + 
                                          (cos_psi*sin_theta)*v.z,
                                          (sin_theta*sin_phi)*v.x + 
                                          (-sin_theta*cos_phi)*v.y + 
                                          (cos_theta)*v.z,
                                          p,
                                          spec);

}

/******************************************************************************************
 * Euler3D_ThermallyPerfect_pState::RotateBack -- Returns an un-rotated primitive state   *
 *                                            re-alinged from the x-axis of the global    *
 *                                            problem.                                    *
 ******************************************************************************************/
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::RotateBack(const Vector3D &norm_dir) const {

   double cos_psi, sin_psi, cos_phi, sin_phi, cos_theta, sin_theta;
   cos_phi = ONE;
   sin_phi = ZERO;
   if (fabs(fabs(norm_dir.x)-ONE) < TOLER) {
      cos_psi = norm_dir.x/fabs(norm_dir.x);
      sin_psi = ZERO;
      cos_theta = ONE;
      sin_theta = ZERO;
   } else {
      cos_psi = norm_dir.x;
      sin_psi = sqrt(norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z);
      cos_theta = norm_dir.y/sqrt(norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z);
      sin_theta = norm_dir.z/sqrt(norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z);
   } /* endif */

   return Euler3D_ThermallyPerfect_pState(rho,
                                          (cos_psi*cos_phi-cos_theta*sin_phi*sin_psi)*v.x + 
                                          (-sin_psi*cos_phi-cos_theta*sin_phi*cos_psi)*v.y + 
                                          (sin_theta*sin_phi)*v.z,
                                          (cos_psi*sin_phi+cos_theta*cos_phi*sin_psi)*v.x + 
                                          (-sin_psi*sin_phi+cos_theta*cos_phi*cos_psi)*v.y + 
                                          (-sin_theta*cos_phi)*v.z,
                                          (sin_theta*sin_psi)*v.x + 
                                          (sin_theta*cos_psi)*v.y + 
                                          (cos_theta)*v.z,
                                          p,
                                          spec);

}

/************************************************************************
 ******************** BOUNDARY CONDITIONS *******************************
 ************************************************************************/

/******************************************************************************************
 * Euler3D_ThermallyPerfect_pState::Reflect -- Return reflected solution state.           *
 ******************************************************************************************/
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::
Reflect(const Euler3D_ThermallyPerfect_pState &W,
        const Vector3D &norm_dir) {

   Vector3D ur_norm, ur_tang, vr_tot;
   Euler3D_ThermallyPerfect_pState Temp; Temp.Copy(W);
   
   ur_norm = dot(W.v, norm_dir)*norm_dir;
   ur_tang = W.v - ur_norm;
   
   ur_norm = -ur_norm;
   vr_tot = ur_norm + ur_tang;
   
   Temp.v = vr_tot;
 
   return (Temp);
       
}

/******************************************************************************************
 * Euler3D_ThermallyPerfect_pState::MovingWall -- Return moving wall boundary state.      *
 ******************************************************************************************/
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::
MovingWall(const Euler3D_ThermallyPerfect_pState &Win,
           const Euler3D_ThermallyPerfect_pState &Wout,
           const Vector3D &norm_dir, 
           const Vector3D &wall_velocity,
           const Vector3D &pressure_gradient,
           const int &TEMPERATURE_BC_FLAG) {

   Euler3D_ThermallyPerfect_pState Temp; Temp.Copy(Win);
   
   if (wall_velocity == Vector3D_ZERO) {
      Temp.v = -Win.v;

   } else {
      double  Wall_velocity_tang ;
      Vector3D ur_norm, ur_tang, vr_tot, uw_tang;
      
      ur_norm = dot(Win.v, norm_dir)*norm_dir;
      ur_tang = Win.v - ur_norm;
      
      uw_tang = wall_velocity - dot(norm_dir,wall_velocity)*norm_dir;
      
      ur_norm = -ur_norm;
      ur_tang = 2.0*uw_tang - ur_tang;
      vr_tot = ur_norm + ur_tang;
      
      Temp.v = vr_tot;
   } /* endif */
   
   //  Fixed Wall Temperature or constant extrapolation for adiabatic.
   if (TEMPERATURE_BC_FLAG == FIXED_TEMPERATURE_WALL) {
      if (pressure_gradient != Vector3D_ZERO) {
         Temp.rho = Wout.p/(Temp.Rtot()*Wout.T());
      } else {
         Temp.rho = Win.p/(Temp.Rtot()*Wout.T());
      } /* endif */
   } /* endif */
   
   return (Temp);

}

/******************************************************************************************
 * Euler3D_ThermallyPerfect_pState::NoSlip -- Return no-slip wall boundary state.         *
 ******************************************************************************************/
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::
NoSlip(const Euler3D_ThermallyPerfect_pState &Win,
       const Euler3D_ThermallyPerfect_pState &Wout,
       const Vector3D &norm_dir,
       const Vector3D &pressure_gradient,
       const int &TEMPERATURE_BC_FLAG) {
   
   return (MovingWall(Win, 
                      Wout, 
                      norm_dir,
                      Vector3D_ZERO,
                      pressure_gradient,
                      TEMPERATURE_BC_FLAG));
   
}

/************************************************************************
 *************************** SOURCE TERMS *******************************
 ************************************************************************/

/**********************************************************************************************
 * Euler3D_ThermallyPerfect_pState::Schemistry -- Finite-rate chemical reaction source terms. *
 **********************************************************************************************/
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::
Schemistry(int &REACT_SET_FLAG) const {

   Euler3D_ThermallyPerfect_cState U_temp;     
   U_temp.Vacuum();

   // Adds concentration rate of change for species 1->N
   if (REACT_SET_FLAG != NO_REACTIONS) {
      React.omega(U_temp, *this);  
   } /* endif */

   return U_temp;

}

/********************************************************************************************
 * Euler3D_ThermallyPerfect_pState::dSchemistrydU -- Jacobian of chemical reaction rate     *
 *                                                   source vector.                         *
 ********************************************************************************************/
void Euler3D_ThermallyPerfect_pState::dSchemistrydU(DenseMatrix &dSwdU) const {

   React.dSwdU<Euler3D_ThermallyPerfect_pState,
               Euler3D_ThermallyPerfect_cState>(dSwdU, *this, false);
}


// Max Diagonal of Jacobian for CFL
double Euler3D_ThermallyPerfect_pState::dSchemistrydU_max_diagonal(void) const {

   double max_diagonal = ONE;
   DenseMatrix dSwdU(num_vars-1, num_vars-1);
   dSwdU.zero();

   React.dSwdU<Euler3D_ThermallyPerfect_pState,Euler3D_ThermallyPerfect_cState>(dSwdU, *this, true);
   for (int i = 0; i < num_vars-1; i++) {
      max_diagonal = max(max_diagonal,fabs(dSwdU(i,i)));
   } /* endfor */

   return max_diagonal;
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_cState::set_species_data -- Assign species data.           *
 ***************************************************************************************/
void Euler3D_ThermallyPerfect_cState::set_species_data(const int &n, 
                                                       const string *S, 
                                                       const char *PATH,
                                                       const int &debug, 
                                                       const double &Mr, 
                                                       const double* Sc,
                                                       const int &trans_data) { 

   // Set the new number of species. 
   ns = n;
   num_vars = ns + NUM_EULER3D_VAR_SANS_SPECIES;

   // Create static memory for NASA data.
   Deallocate_static();
   Allocate_static();

   // Create memory for species data.
   Deallocate();
   Allocate();

   // Read in and assign appropriate NASA data for each species to be used.
   for (int i = 0; i < ns; i++) {
      //overwrite default data  
      specdata[i].Getdata(S[i], PATH, trans_data);
      Schmidt[i] = Sc[i];  
   } /* endfor */

   // Set initial values for the species.
   for (int i = 0; i < ns; i++) {
      rhospec[i].c = rho/ns; 
   } /* endfor */

   // Set data temperature ranges for mixture calculations.
   Temp_low_range();
   Temp_high_range();

   // Set debug information level.
   debug_level = debug;
   
}      

/***************************************************************************************
 * Euler3D_ThermallyPerfect_cState::set_species_values -- Set initial values for       *
 *                                                        species mass fractions.      *
 ***************************************************************************************/
void Euler3D_ThermallyPerfect_cState::set_initial_values(void) {
  Deallocate();
  Allocate();
  for (int i = 0; i < ns; i++) {
    rhospec[i].c = rho/ns; 
  } /* endfor */
}

void  Euler3D_ThermallyPerfect_cState::set_initial_values(const double &value) {
   Deallocate();
   Allocate();
   for (int i = 0; i < ns; i++) {
      rhospec[i].c = value; 
   } /* endfor */
}

void Euler3D_ThermallyPerfect_cState::set_initial_values(const double *rhomfrac) {
   Deallocate();
   Allocate();
   for (int i = 0; i < ns; i++) {
      rhospec[i].c = rhomfrac[i];
   } /* endfor */
}

void Euler3D_ThermallyPerfect_cState::set_initial_values(const Species *rhomfrac) {
   Deallocate();
   Allocate();
   for (int i = 0; i < ns; i++) {
      rhospec[i] = rhomfrac[i];
   } /* endfor */
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_cState::Copy -- Makes a copy of solution state vector.     *
 ***************************************************************************************/
void Euler3D_ThermallyPerfect_cState::Copy(const Euler3D_ThermallyPerfect_cState &U){
   rho = U.rho;
   rhov = U.rhov; 
   E = U.E; 
   for (int i = 0; i < ns; i++) { 
      rhospec[i] = U.rhospec[i];
   } /* endfor */
}
/***************************************************************************************
 * Euler3D_ThermallyPerfect_cState::Temp_low_range -- Set the lower bound for valid T. *
 ***************************************************************************************/
void Euler3D_ThermallyPerfect_cState::Temp_low_range(void) {
  double temp = specdata[0].Low_range();
  for (int i = 0; i < ns; i++) {
    temp = max(specdata[i].Low_range(),temp);
  } /* endfor */
  low_temp_range = temp;  
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_cState::Temp_high_range -- Set the upper bound for valid T.*
 ***************************************************************************************/
void Euler3D_ThermallyPerfect_cState::Temp_high_range(void) {
  double temp = specdata[0].High_range();
  for (int i = 0; i < ns; i++) {
    temp = min(specdata[i].High_range(),temp);
  } /* endfor */
  high_temp_range = temp;  
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_cState::negative_speccheck -- Check species mass fractions.*
 ***************************************************************************************/
bool Euler3D_ThermallyPerfect_cState::negative_speccheck(void) {
   double sum = ZERO;
   double temp = ZERO;
   double LOCAL_TOL = MICRO; 
   //-------- Negative Check ------------//
   for (int i = 0; i < ns-1; i++) {
      temp = rhospec[i].c/rho;
      if (temp>ONE) { //check for > 1.0
         rhospec[i].c =rho;
         temp = ONE;
      } else if (temp < ZERO) {  //check for -ve
         if (temp > -LOCAL_TOL) {  //check for small -ve and set to ZERO 
            rhospec[i].c = ZERO;
            temp = ZERO;
         } else {
            rhospec[i].c = ZERO;
            temp = ZERO;
         } /* endif */
      } /* endif */
      sum += temp;
   } /* endfor */
   temp = max(ONE-sum, ZERO);
   sum += temp;
   rhospec[ns-1].c = rho*temp;
   for (int i = 0; i < ns; i++) {
      rhospec[i].c = rhospec[i].c*(ONE/sum);
   } /* endfor */
   return true;
}

/*****************************************************************************************
 * Euler3D_ThermallyPerfect_cState::Realizable_Solution_Check -- Check physical validity *
 *                                                               validity of solution    *
 *                                                               state.                  *
 *****************************************************************************************/
bool Euler3D_ThermallyPerfect_cState::Realizable_Solution_Check(void) {
   if (rho <= ZERO || !negative_speccheck() || es() <= ZERO) {
      cout << "\n " << CFFC_Name() 
           << " ERROR: Conservative solution state has a negative density, energy,"
           << " and/or species mass fraction.\n";
      return false;
   } else {
      return true;
   } /* endif */
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_cState::Rtot -- Return mixture gas constant  J/(kg*K).     *
 ***************************************************************************************/
double Euler3D_ThermallyPerfect_cState::Rtot(void) const {
  // = sum ( mass fraction * species gas constant)
  double sum = ZERO;
  for ( int i = 0; i < ns; i++) {
    sum += rhospec[i].c * specdata[i].Rs();
  } /* endfor */
  return (sum/rho);
}

/*****************************************************************************************
 * Euler3D_ThermallyPerfect_cState::HeatofFormation -- Return mixture heat of formation. *
 *****************************************************************************************/
double Euler3D_ThermallyPerfect_cState::HeatofFormation(void) const {
   double sum = ZERO;
   for (int i = 0; i < ns; i++) {
      sum += rhospec[i].c*specdata[i].Heatofform();
   } /* endfor */
   return (sum/rho);
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_cState::Cp -- Return mixture heat capacity (const pressure) J/(kg*K). *
 **************************************************************************************************/
double Euler3D_ThermallyPerfect_cState::Cp(void) const {
  // = sum ( mass fraction * species Cp) 
  double Temp = T();
  double sum = ZERO;
  for (int i = 0; i < ns; i++) {
    sum += rhospec[i].c*specdata[i].HeatCapacity_p(Temp);
  } /* endif */
  return (sum/rho);
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_cState::Cv -- Return mixture heat capacity (const volume) J/(kg*K).   *
 **************************************************************************************************/
double Euler3D_ThermallyPerfect_cState::Cv(void) const {
  // = sum ( mass fraction * species Cv)  
  double Temp = T();
  double sum = ZERO;
  for (int i = 0; i < ns; i++) {
    sum += rhospec[i].c*specdata[i].HeatCapacity_v(Temp);
  } /* endfor */
  return (sum/rho);
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_cState::g -- Return mixture specific heat ratio.                      *
 **************************************************************************************************/
double Euler3D_ThermallyPerfect_cState::g(void) const {
  // = Cp / Cv  
  return Cp()/Cv();
}

/**************************************************
  polytropic heat ratio mixture gamma J/(kg*K)
  assuming T=200K as the temperature.
***************************************************/
double Euler3D_ThermallyPerfect_cState::gamma_guess(void) const {
  double sum1 = ZERO;
  double sum2 = ZERO;
  double gamma_s = ZERO;  
  double Temp = 200.0;
  for (int i = 0; i < ns; i++) {
    sum1 += (rhospec[i].c/rho)*specdata[i].Rs();
    gamma_s = ONE/(ONE - specdata[i].Rs()/ specdata[i].HeatCapacity_p(Temp));
    sum2 += (((rhospec[i].c/rho)*specdata[i].Rs()) / (gamma_s - ONE)); 
  } /* endfor */
  return (ONE + sum1/sum2);
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_cState::e -- Return mixture absolute internal energy.                 *
 **************************************************************************************************/
double Euler3D_ThermallyPerfect_cState::e(void) const {
  return ((E - HALF*rhov.sqr()/rho)/rho);
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_cState::es -- Return sensible internal energy.                        *
 **************************************************************************************************/
double Euler3D_ThermallyPerfect_cState::es(void) const {
  return ((E - HALF*rhov.sqr()/rho)/rho-HeatofFormation());
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_cState::h -- Return mixture absolute internal enthalpy.               *
 **************************************************************************************************/
double Euler3D_ThermallyPerfect_cState::h(void) const {
  return (e()+p()/rho);
}

double Euler3D_ThermallyPerfect_cState::h(const double &Temp) const {
  // = sum (mass fraction * species h) 
  double sum = ZERO;  
  for (int i = 0; i < ns; i++) {
    sum += rhospec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform());
  } /* endfor */
  return (sum/rho);
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_cState::hs -- Return sensible internal enthalpy.                      *
 **************************************************************************************************/
double Euler3D_ThermallyPerfect_cState::hs(void) const {
  return (es()+p()/rho);
}

double Euler3D_ThermallyPerfect_cState::hs(const double &Temp) const {
  // = sum (mass fraction * species h) 
  double sum = ZERO;  
  for (int i = 0; i < ns; i++) {
    sum += rhospec[i].c*(specdata[i].Enthalpy(Temp));
  } /* endfor */
  return (sum/rho);
}

/**************************************************************************************************
 * Euler3D_ThermallyPerfect_cState::hprime -- Return dhdT = Cp.                                   *
 **************************************************************************************************/
/**************************************************
   Derivative of specific enthalpy dh/dT
   actually is just Cp as Cp = (dh/dT)_p
***************************************************/
double Euler3D_ThermallyPerfect_cState::hprime(const double &Temp) const {
  double sum = ZERO;  
  for (int i = 0; i < ns; i++) {
    sum += rhospec[i].c*specdata[i].Enthalpy_prime(Temp);
  } /* endfor */
  return (sum/rho);
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_cState::T -- Return mixture temperature.                   *
 ***************************************************************************************/
/********************* Temperature ***********************************
  This gets complicated when trying to obtain the primitive variables 
  from the conserved as E is a function of Temperature taken from 
  a nonlinear equation (ie. a polynomial).   Thus is can't be rearranged
  to find T and p, so a iterative Newtons method has to used. Simple 
  enough, but can be extremely expensive in terms of computation time.
  Oh well ce la vie...

  Also note the "E" is actually rho*E is. E = (rho *(e + HALF*v^2))
  in the conserved state.

  If more that 20 iterations are taken than the look jumps out 
  and gives a warning but will continue.  On average it almost always 
  converges in less than 5 iterations with "tolerance" which is defined
  in the header.
**********************************************************************/
double Euler3D_ThermallyPerfect_cState::T(void) const {
  double T = ZERO;
  double RTOT = Rtot();  
  //--------- Initial Guess ------------------------------//
  //using a polytropic gas assumption with gamma@200;
  double Tguess = (gamma_guess() - ONE)*(E - HALF*rhov.sqr()/rho)/(rho*RTOT);
  //--------- global Newtons method to get T ---------------//
  double A = (E - HALF*rhov*rhov/rho )/rho;
 
  int numit =0;
  double Tmin = low_temp_range;
  double Tmax = high_temp_range;

  //check for start value
  if (Tguess > Tmin && Tguess < Tmax) {
    T=Tguess;
  } else {
    T=Tmin;
  } /* endif */
  
  double fa = h(Tmin) - Tmin*RTOT - A;
  double fn = h(T) - T*RTOT - A;
  double dfn = hprime(T) - RTOT;
  while(fabs(Tmax-Tmin) > CONV_TOLERANCE && 
        fabs(fn) > CONV_TOLERANCE && 
        numit<20 && 
        T >= low_temp_range) {    
    // Newton 
    if (T >= Tmin && T <= Tmax) {
      T = T - fn/dfn;
      if (T >= Tmax) T = HALF*(Tmax - Tmin);	
      //Bisection
    } else {
      T = HALF*(Tmax - Tmin);
    } /* endif */

    //evaluate function and derivative
    fn = h(T) - T*RTOT - A;
    dfn = hprime(T) - RTOT;  

    //change bisection range
    if (fa*fn <=ZERO) {
      Tmax = T;
    } else {
      Tmin = T;
      fa = fn;
    } /* endif */
    numit++;
  } /* endwhile */

  if (numit>=19 || T <= low_temp_range){
    T = max(Tguess,low_temp_range); 	
    if (debug_level) { 
      cout<<"\nTemperature didn't converge in Euler3D_ThermallyPerfect_cState::T(void)";
      cout<<" with polytopic Tguess "<<Tguess<<", or lower than Tmin "<<low_temp_range<<" using "<<T;
    } /* endif */
  }

  return T;

} 

/***************************************************************************************
 * Euler3D_ThermallyPerfect_cState::v -- Return mixture velocity.                      *
 ***************************************************************************************/
Vector3D Euler3D_ThermallyPerfect_cState::v(void) const {
  return (rhov/rho);
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_cState::p -- Return mixture pressure.                      *
 ***************************************************************************************/
double Euler3D_ThermallyPerfect_cState::p(void) const {
  return (rho*Rtot()*T()); 
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_cState::a -- Return mixture sound speed.                   *
 ***************************************************************************************/
double Euler3D_ThermallyPerfect_cState::a(void) const {
  double a_temp;
  double RTOT= Rtot();
  double Temp= T();
  a_temp = (g()*(Rtot()*T()));
  return sqrt(a_temp);
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_cState::M -- Return mixture Mach number.                   *
 ***************************************************************************************/
double Euler3D_ThermallyPerfect_cState::M(void) const {
   return (rhov.abs()/(rho*a()));
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_cState::sum_species -- Sum species mass fractions.         *
 ***************************************************************************************/
double Euler3D_ThermallyPerfect_cState::sum_species(void) const {
   double sum = ZERO;
   for(int i = 0; i < ns-1; i++) {
      sum += rhospec[i].c;
   } /* endfor */
   return sum/rho;
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState::U -- Return primitive solution state vector.       *
 ***************************************************************************************/
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_cState::W(void) {
   Euler3D_ThermallyPerfect_pState Temp;
   Temp.rho = rho;
   Temp.v = v();  
   Temp.p = p();
   for (int i = 0; i < ns; i++) {
      Temp.spec[i] = rhospec[i]/rho;
   } /* endfor */
   return Temp;
}

Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_cState::W(void) const {
   Euler3D_ThermallyPerfect_pState Temp;
   Temp.rho = rho;
   Temp.v = v();  
   Temp.p = p();
   for (int i = 0; i < ns; i++) {
      Temp.spec[i] = rhospec[i]/rho;
   } /* endfor */
   return Temp;
}

Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_cState::
W(const Euler3D_ThermallyPerfect_cState &U) const {
  Euler3D_ThermallyPerfect_pState Temp;
  Temp.rho = U.rho;
  Temp.v = U.v();  
  Temp.p = U.p();
  for (int i = 0; i < U.ns; i++) {
    Temp.spec[i] = U.rhospec[i]/U.rho;
  } /* endfor */
  return Temp;
}

/******************************************************************************************
 * Euler3D_ThermallyPerfect_cState::Rotate -- Returns a rotated primitive state aligned   *
 *                                            with a local x-axis in the norm_dir.        *
 ******************************************************************************************/
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_cState::Rotate(const Vector3D &norm_dir) const {

   double cos_psi, sin_psi, cos_phi, sin_phi, cos_theta, sin_theta;
   cos_phi = ONE;
   sin_phi = ZERO;
   if (fabs(fabs(norm_dir.x)-ONE) < TOLER) {
      cos_psi = norm_dir.x/fabs(norm_dir.x);
      sin_psi = ZERO;
      cos_theta = ONE;
      sin_theta = ZERO;
   } else {
      cos_psi = norm_dir.x;
      sin_psi = sqrt(norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z);
      cos_theta = norm_dir.y/sqrt(norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z);
      sin_theta = norm_dir.z/sqrt(norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z);
   } /* endif */

   return Euler3D_ThermallyPerfect_cState(rho,
                                          (cos_psi*cos_phi-cos_theta*sin_phi*sin_psi)*rhov.x + 
                                          (cos_psi*sin_phi+cos_theta*cos_phi*sin_psi)*rhov.y + 
                                          (sin_psi*sin_theta)*rhov.z,
                                          (-sin_psi*cos_phi-cos_theta*sin_phi*cos_psi)*rhov.x + 
                                          (-sin_psi*sin_phi+cos_theta*cos_phi*cos_psi)*rhov.y + 
                                          (cos_psi*sin_theta)*rhov.z,
                                          (sin_theta*sin_phi)*rhov.x + 
                                          (-sin_theta*cos_phi)*rhov.y + 
                                          (cos_theta)*rhov.z,
                                          E,
                                          rhospec);

}

/******************************************************************************************
 * Euler3D_ThermallyPerfect_cState::RotateBack -- Returns an un-rotated primitive state   *
 *                                            re-alinged from the x-axis of the global    *
 *                                            problem.                                    *
 ******************************************************************************************/
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_cState::RotateBack(const Vector3D &norm_dir) const {

   double cos_psi, sin_psi, cos_phi, sin_phi, cos_theta, sin_theta;
   cos_phi = ONE;
   sin_phi = ZERO;
   if (fabs(fabs(norm_dir.x)-ONE) < TOLER) {
      cos_psi = norm_dir.x/fabs(norm_dir.x);
      sin_psi = ZERO;
      cos_theta = ONE;
      sin_theta = ZERO;
   } else {
      cos_psi = norm_dir.x;
      sin_psi = sqrt(norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z);
      cos_theta = norm_dir.y/sqrt(norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z);
      sin_theta = norm_dir.z/sqrt(norm_dir.y*norm_dir.y + norm_dir.z*norm_dir.z);
   } /* endif */

   return Euler3D_ThermallyPerfect_cState(rho,
                                          (cos_psi*cos_phi-cos_theta*sin_phi*sin_psi)*rhov.x + 
                                          (-sin_psi*cos_phi-cos_theta*sin_phi*cos_psi)*rhov.y + 
                                          (sin_theta*sin_phi)*rhov.z,
                                          (cos_psi*sin_phi+cos_theta*cos_phi*sin_psi)*rhov.x + 
                                          (-sin_psi*sin_phi+cos_theta*cos_phi*cos_psi)*rhov.y + 
                                          (-sin_theta*cos_phi)*rhov.z,
                                          (sin_theta*sin_psi)*rhov.x + 
                                          (sin_theta*cos_psi)*rhov.y + 
                                          (cos_theta)*rhov.z,
                                          E,
                                          rhospec);

}
