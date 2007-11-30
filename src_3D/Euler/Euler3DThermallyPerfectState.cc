/********************** Euler3DThermallyPerfectState.cc **************
  This file defines the various member functions of the 
  Euler3D thermally perfect gaseous mixture class.

   assosicated files:
           Euler3DThermallyPerfectState.h    

*********************************************************************/

// Include required CFFC header files

#ifndef _EULER3D_THERMALLYPERFECT_STATE_INCLUDED
#include "Euler3DThermallyPerfectState.h"
#endif // EULER3D_THERMALLYPERFECT_STATE_INCLUDED   

/***********************************************************/
/**Static Data Member Intialziation**/
int Euler3D_ThermallyPerfect_pState::ns = 1; 
int Euler3D_ThermallyPerfect_pState::num_vars = NUM_EULER3D_VAR_SANS_SPECIES; 
NASARP1311data* Euler3D_ThermallyPerfect_pState::specdata=NULL;
Reaction_set Euler3D_ThermallyPerfect_pState::React;
double Euler3D_ThermallyPerfect_pState::low_temp_range = 200.0;
double Euler3D_ThermallyPerfect_pState::high_temp_range = 300.0;
int Euler3D_ThermallyPerfect_pState::debug_level = 1; 
double* Euler3D_ThermallyPerfect_pState::Schmidt=NULL;

int Euler3D_ThermallyPerfect_cState::ns = 1; 
int Euler3D_ThermallyPerfect_cState::num_vars = NUM_EULER3D_VAR_SANS_SPECIES;
NASARP1311data* Euler3D_ThermallyPerfect_cState::specdata=NULL;   
double Euler3D_ThermallyPerfect_cState::low_temp_range = 200.0;
double Euler3D_ThermallyPerfect_cState::high_temp_range = 300.0;
double* Euler3D_ThermallyPerfect_cState::Schmidt=NULL;
int Euler3D_ThermallyPerfect_cState::debug_level = 1; 

/************* set_species_data ***************************/
//set Global data for Species (STATIC, ie. so only call once! Euler3DInput.h)
/* Get max of the min temperature of the lowest region
   and min of the max temperature of the highest region**/
void Euler3D_ThermallyPerfect_pState::Temp_low_range(void){  
   double temp = specdata[0].Low_range();
   for(int i=0; i<ns; i++){
      temp = max(specdata[i].Low_range(),temp);
   }
   low_temp_range = temp;  
}

void Euler3D_ThermallyPerfect_pState::Temp_high_range(void){
   double temp = specdata[0].High_range();
   for(int i=0; i<ns; i++){
      temp = min(specdata[i].High_range(),temp);
   }
   high_temp_range = temp;  
}

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
   specdata = new NASARP1311data[ns]; 
   Schmidt = new double[ns];

   // Create memory for species data.
   Deallocate();
   spec = new Species[ns];

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

/*****************************************************************
 *****************************************************************
 ** _pState::Sw -- Chemical Reaction Rate Source Terms.     **
 **                                                             **
 ** Using the Reaction class to get the source terms for the    ** 
 ** specific "Reaction_set".                                    ** 
 *****************************************************************
 *****************************************************************/
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::Sw(
   int &REACT_SET_FLAG, const int flow_type) const {
   Euler3D_ThermallyPerfect_cState NEW;     
   NEW.Vacuum();
   bool test = negative_speccheck();
   //Adds concentration rate of change for species 1->N
   if( REACT_SET_FLAG != NO_REACTIONS){
      React.omega(NEW,*this, flow_type, ZERO, ZERO);  
   }
   return NEW;
}

/************* Chemical Source Term Jacobian ****************************/
void Euler3D_ThermallyPerfect_pState::dSwdU(DenseMatrix &dSwdU) const {
   React.dSwdU<Euler3D_ThermallyPerfect_pState,Euler3D_ThermallyPerfect_cState>(dSwdU, *this, false);
}

/************* Max Diagonal of Jacobian for CFL **************************/
double  Euler3D_ThermallyPerfect_pState::dSwdU_max_diagonal(void) const {
   double max_diagonal = ONE;
   DenseMatrix dSwdU(num_vars-1, num_vars-1);
   dSwdU.zero();
   bool test = negative_speccheck();
   React.dSwdU<Euler3D_ThermallyPerfect_pState,Euler3D_ThermallyPerfect_cState>(dSwdU, *this, true);
   for(int i=0; i <  num_vars-1; i++){
      max_diagonal = max(max_diagonal,fabs(dSwdU(i,i)));
   }
   return max_diagonal;
}

/*  Get max and min temperature ranges for data */
void Euler3D_ThermallyPerfect_cState::Temp_low_range(void){  
  double temp = specdata[0].Low_range();
  for(int i=0; i<ns; i++){
    temp = max(specdata[i].Low_range(),temp);
  } 
  low_temp_range = temp;  
}

void Euler3D_ThermallyPerfect_cState::Temp_high_range(void){
  double temp = specdata[0].High_range();
  for(int i=0; i<ns; i++){
    temp = min(specdata[i].High_range(),temp);
  } 
  high_temp_range = temp;  
}

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
   specdata = new NASARP1311data[ns];
   Schmidt = new double[ns];

   // Create memory for species data.
   Deallocate();
   rhospec = new Species[ns];

   // Read in and assign appropriate NASA data for each species to be used.
   for (int i = 0; i < ns; i++){
      //overwrite default data  
      specdata[i].Getdata(S[i], PATH, trans_data);
      Schmidt[i] = Sc[i];  
   } /* endfor */

   // Set initial values for the species.
   for (int i = 0; i < ns; i++){
      rhospec[i].c = rho/ns; 
   } /* endfor */

   // Set data temperature ranges for mixture calculations.
   Temp_low_range();
   Temp_high_range();

   // Set debug information level.
   debug_level = debug;
   
}      

/**************************************************
  mixture molecular mass (kg/mol)
***************************************************/
double Euler3D_ThermallyPerfect_pState::Mass() const{
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
double Euler3D_ThermallyPerfect_pState::Rtot(){
   // = sum ( mass fraction * species gas constant)
   double sum = 0.0;
   for(int i=0; i<ns; i++){
      sum += spec[i].c * specdata[i].Rs();
   }
   return sum;
}

double Euler3D_ThermallyPerfect_pState::Rtot() const{
   // = sum ( mass fraction * species gas constant)
   double sum = 0.0;
   for(int i=0; i<ns; i++){
      sum += spec[i].c * specdata[i].Rs();
   }
   return sum;
}

/****************************************************
  Thermal Conductivity - Mason & Saxena (1958)  W/(m*K)
****************************************************/
double Euler3D_ThermallyPerfect_pState::kappa(void) const{
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
  
  //  or Coffee and Heimerl (1981)
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
  mixture Heat Capacity (const pressure) J/(kg*K)
***************************************************/
double Euler3D_ThermallyPerfect_pState::Cp(void) const{
   // = sum ( mass fraction * species Cp) 
   double Temp = T();
   double sum = 0.0;
   for(int i=0; i<ns; i++){
      sum += spec[i].c*specdata[i].HeatCapacity_p(Temp);
   }
   return sum;
}

double Euler3D_ThermallyPerfect_pState::Cp(const double& TEMP) const{
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
double Euler3D_ThermallyPerfect_pState::Cv(void) const{
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
double Euler3D_ThermallyPerfect_pState::g(void) const{
   // = Cp / Cv  
   return Cp()/Cv();
}
/**************************************************
  Specific Internal Energy
***************************************************/
// etotal = sensible & chemical
double Euler3D_ThermallyPerfect_pState::e(void) const{
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
double Euler3D_ThermallyPerfect_pState::eref(void)const{
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
double Euler3D_ThermallyPerfect_pState::es(void) const{
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
double Euler3D_ThermallyPerfect_pState::h(void) const{
   // = sum (mass fraction * species h) 
   double sum = 0.0;  
   double Temp = T();
   for(int i=0; i<ns; i++){
      sum += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform());
      
   }
   return sum;
}

double Euler3D_ThermallyPerfect_pState::h(const double &Temp) const{
   // = sum (mass fraction * species h) 
   double sum = 0.0;  
   for(int i=0; i<ns; i++){
      sum += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform());
   }
   return sum;
}

double Euler3D_ThermallyPerfect_pState::href(void)const{
   // = sum (mass fraction * species h) 
   double sum = 0.0;  
   double Temp = T();
   for(int i=0; i<ns; i++){ 
      sum += spec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform()
                        - specdata[i].DeltaHref() );
   }
   return sum;
}

double Euler3D_ThermallyPerfect_pState::hs(void) const{
   // = sum (mass fraction * species h) 
   double sum = 0.0;  
   double Temp = T();
   for(int i=0; i<ns; i++){
            
      sum += spec[i].c*(specdata[i].Enthalpy(Temp));
   }
   return sum;
}

double Euler3D_ThermallyPerfect_pState::hs(const double &Temp) const{
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
double Euler3D_ThermallyPerfect_pState::hprime() const{
   double sum = 0.0; 
   double Temp = T();
   for(int i=0; i<ns; i++){
      sum += spec[i].c*specdata[i].Enthalpy_prime(Temp);
   }
   return (sum);
}

double Euler3D_ThermallyPerfect_pState::hprime(double &Temp) const{
   double sum = 0.0;  
   for(int i=0; i<ns; i++){
      sum += spec[i].c*specdata[i].Enthalpy_prime(Temp);
   }
   return (sum);
}

/**************************************************
  Total Energy
***************************************************/
double Euler3D_ThermallyPerfect_pState::E(void) const{
   // E = rho*(e + velocity^2 )  
   return (rho*(e() + HALF*v.sqr())); 
}

/**************************************************
  Total Enthalpy
***************************************************/
double Euler3D_ThermallyPerfect_pState::H(void) const{
  // H = h + velocity^2 
  double TotalEnthalpy = ZERO;
  TotalEnthalpy = (rho*(h() + HALF*v.sqr() ));
  return (TotalEnthalpy);
}

double Euler3D_ThermallyPerfect_pState::Hs(void) const{
  // Hs = hs + velocity^2
  double H_S =  rho*(hs() + HALF*v.sqr() );
  return (H_S);
}

/**************************************************
  polytropic heat ratio mixture gamma J/(kg*K)
  assuming T=200K as the temperature.
***************************************************/
double Euler3D_ThermallyPerfect_pState::gamma_guess(void) const{  
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
   Temperature derived from specific sensible enthalpy
***************************************************/
double Euler3D_ThermallyPerfect_pState::T(double &h_s) const{
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
         cout<<"\nTemperature didn't converge in Euler3D_ThermallyPerfect_cState::T(void)";
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
double Euler3D_ThermallyPerfect_pState::a(void){
   double sum;
   double RTOT= Rtot();
   
//  one way to obtain the speed of sound
//   sum = (p/rho)*(RTOT/( hprime() - RTOT) + ONE);
   
// also could use:  
    sum  =  (g()*(Rtot()*T()));
   return sqrt(sum);
}

double Euler3D_ThermallyPerfect_pState::a(void) const{
   double sum;
   double RTOT= Rtot();
   
//  one way to obtain the speed of sound
//   sum = (p/rho)*(RTOT/( hprime() - RTOT) + ONE);
// also could use 
    sum =  (g()*(Rtot()*T()));
   return sqrt(sum);
}

/*******************************************************************
 ***************** FLUXES ******************************************
 *******************************************************************/
/********************************************************
 * Euler3D_ThermallyPerfect_pState::F -- Inviscid flux (x-direction).   *
 ********************************************************/
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::F(void) {
  Euler3D_ThermallyPerfect_cState Temp;
  
  Temp.rho = rho*v.x;
  Temp.rhov.x = rho*sqr(v.x) + p ;
  Temp.rhov.y = rho*v.x*v.y;
  Temp.rhov.z = rho*v.x*v.z;
  Temp.E = v.x*H();
  
  //multispecies transport
  for(int i=0; i<ns; i++){
     Temp.rhospec[i].c = rho*v.x*spec[i].c;
  }
  return (Temp);
}


Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::F(void) const{
  Euler3D_ThermallyPerfect_cState Temp;
  
  Temp.rho = rho*v.x;
  Temp.rhov.x = rho*sqr(v.x) + p ;
  Temp.rhov.y = rho*v.x*v.y;
  Temp.rhov.z = rho*v.x*v.z;
  Temp.E = v.x*H();
  
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
 * Euler3D_ThermallyPerfect_pState::lambda -- Eigenvalue(s) (x-direction).    *
 ************************************************************/

Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::lambda_x(void) {
  double c = a();
  Euler3D_ThermallyPerfect_pState Temp;
  Temp.rho = v.x - c;
  Temp.v.x = v.x;
  Temp.v.y = v.x;
  Temp.v.z = v.x;
  Temp.p = v.x + c;
  for(int i=0; i<ns; i++){
     Temp.spec[i].c = v.x;
  }
 

  
  return (Temp);
}
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::lambda_x(void) const {
  double c = a();
  Euler3D_ThermallyPerfect_pState Temp;
  Temp.rho = v.x - c;
  Temp.v.x = v.x;
  Temp.v.y = v.x;
  Temp.v.z = v.x;
  Temp.p = v.x + c;
  for(int i=0; i<ns; i++){
     Temp.spec[i].c = v.x;
  }
 

  
  return (Temp);
}


Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::rc_x(
   const int &index) const {
  
   switch(index){  
   case 1: 
      return (Euler3D_ThermallyPerfect_cState(
                 ONE, v.x-a(), v.y, v.z, H()/rho-v.x*a(), spec));
   case 2:
      return (Euler3D_ThermallyPerfect_cState(
                 ONE, v.x, v.y, v.z, H()/rho-Cp()*T(), spec)); 
   case 3:
      return (Euler3D_ThermallyPerfect_cState(
                 ZERO, ZERO, rho, ZERO, rho*v.y, ZERO));
   case 4: 
      return (Euler3D_ThermallyPerfect_cState(
                 ZERO, ZERO, ZERO, rho, rho*v.z, ZERO));
   case 5: 
      return (Euler3D_ThermallyPerfect_cState(
                 ONE, v.x+a(), v.y, v.z, H()/rho+v.x*a(), spec));
      
   default: 
      Euler3D_ThermallyPerfect_cState NEW(ZERO);
      double RTOT = Rtot();
      double TEMP = p/(rho*RTOT);      
      NEW.E = rho*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP) + 
              specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Heatofform() - 
              Cp(TEMP)*TEMP*specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs()/RTOT);
      double PHI =specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP)-
                  specdata[num_vars- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP)-
                  (Cp(TEMP) -RTOT)*TEMP*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs() - 
                  specdata[num_vars - (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs())/RTOT;
            
      NEW.rhospec[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].c = rho*PHI; 
            
      return NEW;

   };
   

}

// Primitive Left Eigenvector -- (x-direction)
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::lp_x(const int &index) {
 
   switch(index){  
   case 1 :
      return (Euler3D_ThermallyPerfect_pState(
                 ZERO, -HALF*rho/a(), ZERO, ZERO, HALF/(a()*a()),  ZERO));
   case 2:
      
      return (Euler3D_ThermallyPerfect_pState(
                 ONE, ZERO, ZERO, ZERO, -ONE/(a()*a()),  ZERO));
   case 3:
      return  (Euler3D_ThermallyPerfect_pState(ZERO, ZERO, ONE, ZERO, ZERO, ZERO));
   case 4 :
      return  (Euler3D_ThermallyPerfect_pState(ZERO, ZERO, ZERO, ONE, ZERO, ZERO));
   case 5 :
      
      return (Euler3D_ThermallyPerfect_pState(
                 ZERO, HALF*rho/a(), ZERO, ZERO, HALF/(a()*a()), ZERO));
   default:
      Euler3D_ThermallyPerfect_pState NEW(ZERO);
      NEW.spec[index-(NUM_EULER3D_VAR_SANS_SPECIES+1)].c = ONE;
      return NEW;
      
   };
    
}

Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::rc_x(
   const int &index) {
   
   switch(index){  
   case 1: 
      return (Euler3D_ThermallyPerfect_cState(
                 ONE, v.x-a(), v.y, v.z, H()/rho-v.x*a(), spec));
   case 2:
      return (Euler3D_ThermallyPerfect_cState(
                 ONE, v.x, v.y, v.z, H()/rho-Cp()*T(), spec)); 
   case 3:
      return (Euler3D_ThermallyPerfect_cState(
                 ZERO, ZERO, rho, ZERO, rho*v.y, ZERO));
   case 4: 
      return (Euler3D_ThermallyPerfect_cState(
                 ZERO, ZERO, ZERO, rho, rho*v.z, ZERO));
   case 5: 
      return (Euler3D_ThermallyPerfect_cState(
                 ONE, v.x+a(), v.y, v.z, H()/rho+v.x*a(), spec));
      
   default: 
      Euler3D_ThermallyPerfect_cState NEW(ZERO);
      double RTOT = Rtot();
      double TEMP = p/(rho*RTOT);      
      NEW.E = rho*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP) + 
              specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Heatofform() - 
              Cp(TEMP)*TEMP*specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs()/RTOT);
      double PHI = specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP)-
                   specdata[num_vars- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP)-
                   (Cp(TEMP) -RTOT)*TEMP*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs() - 
                   specdata[num_vars - (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs())/RTOT;
            
      NEW.rhospec[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].c = rho*PHI; 
            
      return NEW;

   };
   
 
}

// Primitive Left Eigenvector -- (x-direction)
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::lp_x(const int &index) const {

   switch(index){  
   case 1 :
      return (Euler3D_ThermallyPerfect_pState(
                 ZERO, -HALF*rho/a(), ZERO, ZERO, HALF/(a()*a()),  ZERO));
   case 2:
      
      return (Euler3D_ThermallyPerfect_pState(
                 ONE, ZERO, ZERO, ZERO, -ONE/(a()*a()),  ZERO));
   case 3:
      return  (Euler3D_ThermallyPerfect_pState(ZERO, ZERO, ONE, ZERO, ZERO, ZERO));
   case 4 :
      return  (Euler3D_ThermallyPerfect_pState(ZERO, ZERO, ZERO, ONE, ZERO, ZERO));
   case 5 :
      
      return (Euler3D_ThermallyPerfect_pState(
                 ZERO, HALF*rho/a(), ZERO, ZERO, HALF/(a()*a()), ZERO));
   default:
      Euler3D_ThermallyPerfect_pState NEW(ZERO);
      NEW.spec[index-(NUM_EULER3D_VAR_SANS_SPECIES+1)].c = ONE;
      return NEW;
      
   };
    
  
}


/**************************************************************************
 * Euler3D_ThermallyPerfect_pState -- Binary arithmetic operators.        *
 **************************************************************************/


Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::operator +(const Euler3D_ThermallyPerfect_pState &W) const { 
   if(ns == W.ns){ //check that species are equal 
      Euler3D_ThermallyPerfect_pState Temp(W.rho,W.v,W.p);
      Temp.Copy(*this);
      Temp += W;
      return Temp;
   }else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
  
}

//------------------ Subtraction ------------------------//
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::operator -(const Euler3D_ThermallyPerfect_pState &W) const{
   if(ns == W.ns){ //check that species are equal
      Euler3D_ThermallyPerfect_pState Temp(W.rho,W.v,W.p);
      Temp.Copy(*this);
      Temp -= W;
      return Temp;
   }else {
         cerr<<"\n Mismatch in number of species \n";
         exit(1);
      }

   }

//---------------- Scalar Multiplication ------------------//
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::operator *(
   const double &a) const{

   Euler3D_ThermallyPerfect_pState Temp(rho,v,p);
   Temp.Copy(*this);
   Temp.rho = rho*a;  Temp.v = v*a; Temp.p = p*a;
  for( int i=0; i<ns; i++){
    Temp.spec[i] = spec[i]*a;
  } 
  
  return(Temp);
}

Euler3D_ThermallyPerfect_pState operator *(const double &a, 
                                           const Euler3D_ThermallyPerfect_pState &W){

   Euler3D_ThermallyPerfect_pState Temp;
  //Temp.Copy(W);
   Temp.rho = W.rho*a;  Temp.v = W.v*a; Temp.p = W.p*a;
   
   for( int i=0; i<W.ns; i++){
      Temp.spec[i] = W.spec[i]*a;
   } 

  return(Temp);

}

//--------------- Scalar Division ------------------------//
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::operator /(const double &a) const {

   Euler3D_ThermallyPerfect_pState Temp(rho,v,p);
   Temp.Copy(*this);
   Temp.rho = rho/a; Temp.v = v/a; Temp.p = p/a; 
  
  for(int i=0; i<ns; i++){
     Temp.spec[i] = spec[i]/a; 
  } 

  return(Temp);

}

//----------------- Inner Product ------------------------//
double Euler3D_ThermallyPerfect_pState::operator *(const Euler3D_ThermallyPerfect_pState &W) const{

   double sum=0.0;
   if(ns == W.ns){ //check that species are equal
      for(int i=0; i<ns; i++){
         sum += spec[i]*W.spec[i];
      }  
      return (rho*W.rho + v*W.v + p*W.p + sum);
   } else {
      cerr<<"\n Mismatch in number of species \n";
      exit(1);
   }

}

//----------- solution state product operator ------------//
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::operator ^(const Euler3D_ThermallyPerfect_pState &W) const{

     if(ns == W.ns){ //check that species are equal
        Euler3D_ThermallyPerfect_pState Temp(rho,v,p);
        Temp.Copy(*this);
        Temp.rho = rho*W.rho;
        Temp.v.x = v.x*W.v.x;
        Temp.v.y = v.y*W.v.y;
        Temp.v.z = v.z*W.v.z;
        Temp.p = p*W.p;
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
Euler3D_ThermallyPerfect_pState& Euler3D_ThermallyPerfect_pState::operator =(const Euler3D_ThermallyPerfect_pState &W){

   //self assignment protection
   if( this != &W){   
      //copy assignment
      rho = W.rho;
      v = W.v; 
      p = W.p; 

    if ( ns == W.ns){
      for(int i=0; i<ns; i++){
	spec[i] = W.spec[i];
      }   

    }   
    else {
      cerr<<"\n Mismatch in number of species \n ";
    }  
  }
  return (*this);

}

/*************************************************************************
 * Euler3D_ThermallyPerfect_pState -- Shortcut arithmetic operators.     *
 *************************************************************************/
Euler3D_ThermallyPerfect_pState& Euler3D_ThermallyPerfect_pState::operator +=(const Euler3D_ThermallyPerfect_pState &W){

   rho += W.rho;
   v += W.v; 
   p += W.p; 
   
   for( int i=0; i<ns; i++){
      spec[i] += W.spec[i];
   } 
   
   return (*this);

}

Euler3D_ThermallyPerfect_pState& Euler3D_ThermallyPerfect_pState::operator -=(const Euler3D_ThermallyPerfect_pState &W) {

   rho -= W.rho;
   v -= W.v;
   p -= W.p;
   
   for(int i=0; i<ns; i++){
      spec[i] -= W.spec[i];
   }  
   return (*this); 

}

/**************************************************************************
 * Euler3D_ThermallyPerfect_pState -- Unary arithmetic operators.         *
 **************************************************************************/
Euler3D_ThermallyPerfect_pState operator +(const Euler3D_ThermallyPerfect_pState &W) {  
   return (Euler3D_ThermallyPerfect_pState(W.rho,W.v,W.p,W.spec));
}

Euler3D_ThermallyPerfect_pState operator -(const Euler3D_ThermallyPerfect_pState &W) {
   Species *spt= new Species[W.ns];
   for(int i=0; i<W.ns; i++){
      spt[i] = -W.spec[i]; 
   }  
   
   Euler3D_ThermallyPerfect_pState Temp(-W.rho,-W.v, -W.p, spt);
   
   delete[] spt;
   return(Temp);
}

/********************************************************
 * Euler3D_ThermallyPerfect_pState -- Relational operators.               *
 ********************************************************/
int operator ==(const Euler3D_ThermallyPerfect_pState &W1, 
                const Euler3D_ThermallyPerfect_pState &W2) {
   if(W1.ns == W2.ns){ //check that species are equal
      bool Temp;
      for(int i=0; i<W1.ns; i++){
         if( W1.spec[i] == W2.spec[i] ){
            Temp = true;
         } else {
            Temp = false;
            break;
      }  
         return (W1.rho == W2.rho && W1.v == W2.v && W1.p == W2.p
                 && Temp == true);
    }
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
}

int operator !=(const Euler3D_ThermallyPerfect_pState &W1, 
                const Euler3D_ThermallyPerfect_pState &W2) {
   if(W1.ns == W2.ns){ //check that species are equal
      bool Temp = true;
      for(int i=0; i<W1.ns; i++){
         if( W1.spec[i] != W2.spec[i] ){
            Temp = false;
            break;
         } 
         return (W1.rho != W2.rho || W1.v != W2.v || W1.p != W2.p 
                 || Temp != true);
    }
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
}

/* Euler3D_ThermallyPerfect_pState -- Input-output operators.*/
ostream &operator << (ostream &out_file, 
                      const Euler3D_ThermallyPerfect_pState &W) {
   out_file.precision(10);
   out_file.setf(ios::scientific);
   out_file << " " << W.rho  << " " << W.v.x << " " << W.v.y 
            << " " << W.v.z << " " << W.p;
   for( int i=0; i<W.ns; i++){
      out_file<<" "<<W.spec[i];
   }
   out_file.unsetf(ios::scientific);
   return (out_file);
}

istream &operator >> (istream &in_file, Euler3D_ThermallyPerfect_pState &W) {
   in_file.setf(ios::skipws);
   in_file >> W.rho >> W.v.x >> W.v.y >> W.v.z >>W.p;
   //W.set_initial_values();
   for( int i=0; i<W.ns; i++){
      in_file>>W.spec[i];
   }
   in_file.unsetf(ios::skipws);
   return (in_file);
}


//  mixture gas constant  J/(kg*K)
double Euler3D_ThermallyPerfect_cState::Rtot() const{
  // = sum ( mass fraction * species gas constant)
  double sum = 0.0;
  for(int i=0; i<ns; i++){
    sum += rhospec[i].c * specdata[i].Rs();
  }
 
  return (sum/rho);
}

//  mixture Heat Capacity (const pressure) J/(kg*K)
double Euler3D_ThermallyPerfect_cState::Cp(void) const{
  // = sum ( mass fraction * species Cp) 
  double Temp = T();
  double sum = 0.0;
  for(int i=0; i<ns; i++){
    sum += rhospec[i].c*specdata[i].HeatCapacity_p(Temp);
  }
  return (sum/rho);
}

//  mixture Heat Capacity (const volume) J/(kg*K)
double Euler3D_ThermallyPerfect_cState::Cv(void) const{
  // = sum ( mass fraction * species Cv)  
  double Temp = T();
  double sum = 0.0;
  for(int i=0; i<ns; i++){
    sum += rhospec[i].c*specdata[i].HeatCapacity_v(Temp);
  }
  return (sum/rho);
}

//  mixture Heat Ratio gamma J/(kg*K)
double Euler3D_ThermallyPerfect_cState::g(void) const{
  // = Cp / Cv  
  return Cp()/Cv();
}

//  Specific Internal Energy
double Euler3D_ThermallyPerfect_cState::e(void) const{
  // = sum (mass fraction * species e) 
  double sum = 0.0;
  double Temp = T();
  for(int i=0; i<ns; i++){ //(Enthalpy(Temp) - (R/mol_mass)*Temp)
     sum += rhospec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform() 
                          - specdata[i].Rs()*Temp);
          
  }
  return (sum/rho);
}

double Euler3D_ThermallyPerfect_cState::es(void) const{
  // = sum (mass fraction * species e) 
  double sum = 0.0;
  double Temp = T();
  for(int i=0; i<ns; i++){ //(Enthalpy(Temp) - (R/mol_mass)*Temp)
    sum += rhospec[i].c*(specdata[i].Enthalpy(Temp) - specdata[i].Rs()*Temp);
  }
  return (sum/rho);
}

// Specific Absolute enthalpy
double Euler3D_ThermallyPerfect_cState::h(const double &Temp) const{
  // = sum (mass fraction * species h) 
 double sum = 0.0;  
 for(int i=0; i<ns; i++){
   sum += rhospec[i].c*(specdata[i].Enthalpy(Temp) + specdata[i].Heatofform());
 }
 return (sum/rho);
}

double Euler3D_ThermallyPerfect_cState::hs(const double &Temp) const{
  // = sum (mass fraction * species h) 
 double sum = 0.0;  
 for(int i=0; i<ns; i++){
   sum += rhospec[i].c*(specdata[i].Enthalpy(Temp));
 }
 return (sum/rho);
}

//   Derivative of specific enthalpy dh/dT
//   actually is just Cp as Cp = (dh/dT)_p
double Euler3D_ThermallyPerfect_cState::hprime(const double &Temp) const{
 double sum = 0.0;  
 for(int i=0; i<ns; i++){
   sum += rhospec[i].c*specdata[i].Enthalpy_prime(Temp);
 }
 return (sum/rho);
}

/************* Mixture Heats of Formation *********/
double Euler3D_ThermallyPerfect_cState::heatofform(void) const{ 
  double sum = 0.0;
  for(int i=0; i<ns; i++){ 
    sum += rhospec[i].c*specdata[i].Heatofform();
  }
  return (sum);
}

//  polytropic heat ratio mixture gamma J/(kg*K)
//  assuming T=273K as the temperature.
double Euler3D_ThermallyPerfect_cState::gamma_guess(void) const{
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

  Also note the "E" is actually rho*E is. E = (rho *(e + HALF*v^2+k))
  in the conserved state.

  If more that 20 iterations are taken than the look jumps out 
  and gives a warning but will continue.  On average it almost always 
  converges in less than 5 iterations with "tolerance" which is defined
  in the header.
**********************************************************************/
double Euler3D_ThermallyPerfect_cState::T(void) const{
  double T = ZERO;
  double RTOT = Rtot();  
  //--------- Initial Guess ------------------------------//
  //using a polytropic gas assumption with gamma@200;
  double Tguess = (gamma_guess() - ONE)*(E - HALF*rhov.sqr()/rho)/(rho*RTOT);
  //--------- global newtons method to get T ---------------//
  double A = (E - HALF*rhov*rhov/rho )/rho;
 
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
      cout<<"\nTemperature didn't converge in Euler3D_ThermallyPerfect_cState::T(void)";
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
double Euler3D_ThermallyPerfect_cState::a(void) const{
  double sum;
  double RTOT= Rtot();
  double Temp= T();
 
  //  one way to obtain the speed of sound 
  //sum = RTOT*Temp*(RTOT/( hprime(Temp) - RTOT) + ONE);
  sum  =  (g()*(Rtot()*T()));
 
  return sqrt(sum);
}

//----------------- Addition -----------------------------//
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_cState::operator +(
   const Euler3D_ThermallyPerfect_cState &U) const{ 
   
   //f(ns == U.ns){ //check that species are equal   
      Euler3D_ThermallyPerfect_cState Temp(U.rho,U.rhov,U.E);
      Temp.Copy(*this);
      Temp += U;
      return Temp;
  //  } else {
//       cerr<<"\n Mismatch in number of species \n";
//       exit(1);
//    }
}

//------------------ Subtraction ------------------------//
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_cState::operator -(
   const Euler3D_ThermallyPerfect_cState &U) const{
   // if(ns == U.ns){ //check that species are equal
   Euler3D_ThermallyPerfect_cState Temp(U.rho,U.rhov,U.E);
   Temp.Copy(*this);
   Temp -= U;
   return Temp;
// } else {
//     cerr<<"\n Mismatch in number of species \n";
//     exit(1);
//   }
}


//---------------- Scalar Multiplication ------------------//
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_cState::operator *(
   const double &a) const{
   Euler3D_ThermallyPerfect_cState Temp(rho,rhov,E);
   Temp.Copy(*this);
   Temp.rho = rho*a;  Temp.rhov = rhov*a; Temp.E = E*a;
   for( int i=0; i<ns; i++){
      Temp.rhospec[i] = rhospec[i]*a;
   } 
   
   return(Temp);
}

Euler3D_ThermallyPerfect_cState operator *(
   const double &a, 
   const Euler3D_ThermallyPerfect_cState &U){
   Euler3D_ThermallyPerfect_cState Temp;
   Temp.rho = U.rho*a;  Temp.rhov = U.rhov*a; Temp.E = U.E*a;
   for( int i=0; i<U.ns; i++){
      Temp.rhospec[i] = U.rhospec[i]*a;
   } 
   return(Temp);
}
//--------------- Scalar Division ------------------------//
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_cState::operator /(
   const double &a) const {
   Euler3D_ThermallyPerfect_cState Temp(rho,rhov,E);
   Temp.Copy(*this);
   Temp.rho = rho/a; Temp.rhov = rhov/a; Temp.E = E/a;

  for(int i=0; i<ns; i++){
     Temp.rhospec[i] = rhospec[i]/a; 
  } 

  return(Temp);
}

//----------------- Inner Product ------------------------//
double Euler3D_ThermallyPerfect_cState::operator *(
   const Euler3D_ThermallyPerfect_cState &U) const{
   double sum=0.0;
   //if(ns == U.ns){ //check that species are equal
      for(int i=0; i<ns; i++){
         sum += rhospec[i]*U.rhospec[i];
      }  
      return (rho*U.rho + rhov*U.rhov + E*U.E + sum);
   // } else {
//       cerr<<"\n Mismatch in number of species \n";
//       exit(1);
//    }
}

//----------- solution state product operator ------------//
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_cState::operator ^(
   const Euler3D_ThermallyPerfect_cState &U) const {
   // if(ns == U.ns){ //check that species are equal
      Euler3D_ThermallyPerfect_cState Temp(rho,rhov,E);
      Temp.Copy(*this);
      Temp.rho = rho*U.rho;
      Temp.rhov.x = rhov.x*U.rhov.x;
      Temp.rhov.y = rhov.y*U.rhov.y;
      Temp.rhov.z = rhov.z*U.rhov.z;
      Temp.E = E*U.E;

      for(int i=0; i<ns; i++){
         Temp.rhospec[i] = rhospec[i]*U.rhospec[i];
      }  
      return(Temp);
  //  } else {
//       cerr<<"\n Mismatch in number of species \n";
//       exit(1);
//    }
}




//----------------- Assignment ----------------------------//
Euler3D_ThermallyPerfect_cState& Euler3D_ThermallyPerfect_cState::operator =(
   const Euler3D_ThermallyPerfect_cState &U){
   //self assignment protection
   if( this != &U){   
      //copy assignment
      rho = U.rho;
      rhov = U.rhov; 
      E = U.E; 
      
      if ( ns == U.ns){
         for(int i=0; i<ns; i++){
            rhospec[i] = U.rhospec[i];
         }
         
      } else {
         cerr<<"\n Mismatch in number of species \n ";
      }  
   }
   return (*this);
}


/********************************************************
 * Euler3D_ThermallyPerfect_cState -- Shortcut arithmetic operators.     *
 ********************************************************/
Euler3D_ThermallyPerfect_cState& Euler3D_ThermallyPerfect_cState::operator +=(
   const Euler3D_ThermallyPerfect_cState &U){
   rho += U.rho;
   rhov += U.rhov; 
   E += U.E;
   
   for( int i=0; i<ns; i++){
      rhospec[i] += U.rhospec[i];
   } 

   return (*this);
}

Euler3D_ThermallyPerfect_cState& Euler3D_ThermallyPerfect_cState::operator -=(
   const Euler3D_ThermallyPerfect_cState &U) {
  
   rho -= U.rho;
   rhov -= U.rhov;
   E -= U.E;
   
  for(int i=0; i<ns; i++){
    rhospec[i] -= U.rhospec[i];
  }  

  return (*this); 
}

/********************************************************
 * Euler3D_ThermallyPerfect_cState -- Unary arithmetic operators.        *
 ********************************************************/
Euler3D_ThermallyPerfect_cState operator -(
   const Euler3D_ThermallyPerfect_cState &U) {
   Species *spt= new Species[U.ns];
   for(int i=0; i<U.ns; i++){
      spt[i] = -U.rhospec[i]; 
   }  
   
   Euler3D_ThermallyPerfect_cState Temp(-U.rho,-U.rhov,-U.E, spt);
  
  delete[] spt;
  return(Temp);
}

/********************************************************
 * Euler3D_ThermallyPerfect_cState -- Relational operators.              *
 ********************************************************/
int operator ==(const Euler3D_ThermallyPerfect_cState &U1, 
                const Euler3D_ThermallyPerfect_cState &U2) {
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
                 &&Temp == true);
      }
   } else {
      cerr<<"\n Mismatch in number of species \n";
      exit(1);
   }
}

int operator !=(const Euler3D_ThermallyPerfect_cState &U1, 
                const Euler3D_ThermallyPerfect_cState &U2) {
   if(U1.ns == U2.ns){ //check that species are equal
      bool Temp = true;
      for(int i=0; i<U1.ns; i++){
         if( U1.rhospec[i] != U2.rhospec[i] ){
            Temp = false;
            break;
         } 
         return (U1.rho != U2.rho || U1.rhov != U2.rhov || U1.E != U2.E 
                 || Temp != true);
      }
   } else {
      cerr<<"\n Mismatch in number of species \n";
      exit(1);
   }
}

/********************************************************
 * Euler3D_ThermallyPerfect_cState -- Input-output operators.            *
 ********************************************************/
ostream &operator << (ostream &out_file, const Euler3D_ThermallyPerfect_cState &U) {
  //out_file.precision(20);
  out_file.setf(ios::scientific);
  out_file << " " << U.rho  << " " << U.rhov.x << " " << U.rhov.y << " " 
           << U.rhov.z << " "<< U.E;
  for( int i=0; i<U.ns; i++){
     out_file<<" "<<U.rhospec[i];
  } 
  out_file.unsetf(ios::scientific);
  return (out_file);
}

istream &operator >> (istream &in_file, Euler3D_ThermallyPerfect_cState &U) {
   in_file.setf(ios::skipws);
   in_file >> U.rho >> U.rhov.x >> U.rhov.y >> U.rhov.z >> U.E;
  //U.set_initial_values();
  for( int i=0; i<U.ns; i++){
    in_file>>U.rhospec[i]; 
  } 

  in_file.unsetf(ios::skipws);
  return (in_file);
}


/********************************************************
 * Routine: RoeAverage (Roe Averages)                   *
 *                                                      *
 * This function returns the Roe-averaged primitive     *
 * solution state given left and right primitive        *
 * solution variables.  See Roe (1981).                 *
 *                                                      *
 ********************************************************/
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::RoeAverage(
   const Euler3D_ThermallyPerfect_pState &Wl,
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

   
   for(int i=0; i<Wl.ns; i++){
      Temp.spec[i].c = (srhol*Wl.spec[i].c + srhor*Wr.spec[i].c)/(srhol+srhor);
   }


   //Temp.p = (srhol*Wl.p + srhor*Wr.p)/(srhol+srhor);
   Ha = (srhol*Hl+srhor*Hr)/(srhol+srhor);
   ha = Ha - HALF*(sqr(Temp.v.x)+sqr(Temp.v.y)+sqr(Temp.v.z));

   
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
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::FluxHLLE_x(
   const Euler3D_ThermallyPerfect_pState &Wl,
   const Euler3D_ThermallyPerfect_pState &Wr) {
   
   double wavespeed_l, wavespeed_r;

    
    /* solnvec in  Wa (lambdas_l, lambdas_r, lambdas_a) 
       is allocated using new  */ 
   Euler3D_ThermallyPerfect_pState Wa, lambdas_l, lambdas_r, lambdas_a;
   Euler3D_ThermallyPerfect_cState Flux, dUrl;
   
   int num_vars = Wl.num_vars;
      
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
   wavespeed_r = max(lambdas_r[num_vars-lambdas_r.ns],
                     lambdas_a[num_vars-lambdas_a.ns]);


  wavespeed_l = min(wavespeed_l, ZERO);
  wavespeed_r = max(wavespeed_r, ZERO);
   
   if (wavespeed_l >= ZERO) {
      Flux = Wl.F();
   } else if (wavespeed_r <= ZERO) {
      Flux = Wr.F();
   } else {
      Flux =   ( (wavespeed_r*Wl.F()-wavespeed_l*Wr.F())
                 +(wavespeed_l*wavespeed_r)*dUrl)/
         (wavespeed_r-wavespeed_l);
   } /* endif */
   

   /* Return solution flux. */
   
   return (Flux);
   
}

Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::FluxHLLE_x(
   const Euler3D_ThermallyPerfect_cState &Ul,
   const Euler3D_ThermallyPerfect_cState &Ur) {
   return (FluxHLLE_x(Ul.W(), Ur.W()));
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
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::FluxHLLE_n(
   const Euler3D_ThermallyPerfect_pState &Wl,
   const Euler3D_ThermallyPerfect_pState &Wr,
   const Vector3D &norm_dir) {

 
 double sin_beta, cos_beta, sin_alpha, cos_alpha;
   
  //  //solnvec in  Wl_rotated (Wr_rotated) is allocated using new 
//    Euler3D_ThermallyPerfect_pState Wl_rotated, Wr_rotated;
//    Euler3D_ThermallyPerfect_cState Flux, Flux_rotated;
   
//    /* Determine the direction cosine's for the frame
//       rotation. */
//    sin_beta = norm_dir.z;
//    cos_beta = sqrt(sqr(norm_dir.x)+sqr(norm_dir.y));
//    if(cos_beta == ZERO){
//       sin_alpha = ZERO;
//       cos_alpha = ZERO;
//    }else{
//       sin_alpha = norm_dir.y/cos_beta;
//       cos_alpha = norm_dir.x/cos_beta;
//    }
   
//    /* Apply the frame rotation and evaluate left and right
//       solution states in the local rotated frame defined
//       by the unit normal vector. */
//    Wl_rotated.Copy(Wl);
//    Wr_rotated.Copy(Wr);


//    Wl_rotated.rho = Wl.rho;
//    Wl_rotated.v.x = (Wl.v.x*cos_alpha +
//                      Wl.v.y*sin_alpha)*cos_beta+ Wl.v.z*sin_beta;
//    Wl_rotated.v.y =  (Wl.v.x*cos_alpha +
//                       Wl.v.y*sin_alpha)*sin_beta+ Wl.v.z*cos_beta;
//    Wl_rotated.v.z = ZERO;
//    Wl_rotated.p = Wl.p;
    
//    Wr_rotated.rho = Wr.rho;
//    Wr_rotated.v.x = (Wr.v.x*cos_alpha +
//                      Wr.v.y*sin_alpha)*cos_beta+ Wr.v.z*sin_beta;
//    Wr_rotated.v.y =  (Wr.v.x*cos_alpha +
//                       Wr.v.y*sin_alpha)*sin_beta+ Wr.v.z*cos_beta;
    
//    Wr_rotated.v.z = ZERO;
//    Wr_rotated.p = Wr.p;

//     /* Evaluate the intermediate state solution 
//        flux in the rotated frame. */

//    Flux_rotated = FluxHLLE_x(Wl_rotated, Wr_rotated);
  
//    /* Rotate back to the original Cartesian reference
//       frame and return the solution flux. */
  
//    Flux.Copy(Flux_rotated);

   
//    Flux.rho = Flux_rotated.rho;
//    Flux.rhov.x = (Flux_rotated.rhov.x*cos_beta + Flux_rotated.rhov.y*sin_beta)*cos_alpha;
//    Flux.rhov.y = (Flux_rotated.rhov.x*cos_beta + Flux_rotated.rhov.y*sin_beta)*sin_alpha;
//    Flux.rhov.z = (Flux_rotated.rhov.x*sin_beta + Flux_rotated.rhov.y* cos_beta);
   
//    Flux.E = Flux_rotated.E;
  
//    Flux.zero_non_sol();
 

//    return (Flux);
  
 double Wl_ur_norm, Wl_ur_tang;
 double Wr_ur_norm, Wr_ur_tang ;
 double Wr_ur_tang_z;
 
 Vector3D Flux_rotated_x, Flux_rotated_tang_y, Flux_rotated_tang_z ;
 Vector3D Wl_ur_tang_vector, Wr_ur_tang_vector;
 Vector3D Wl_ur_tang_unit_vector, Wr_ur_tang_unit_vector;
 Vector3D Wr_ur_tang_z_vector, Wr_ur_tang_z_unit_vector;
 
 //solnvec in  Wl_rotated (Wr_rotated) is allocated using new 
 Euler3D_ThermallyPerfect_pState Wl_rotated, Wr_rotated;
 Euler3D_ThermallyPerfect_cState Flux, Flux_rotated;
 
 
 /* Apply the frame rotation and evaluate left and right
    solution states in the local rotated frame defined
    by the unit normal vector. */
 Wl_rotated.Copy(Wl);
 Wr_rotated.Copy(Wr);
 
 // Left state velocity in rotated frame
 Wl_ur_norm = dot(Wl.v, norm_dir);
 Wl_ur_tang = abs(Wl.v - Wl_ur_norm*norm_dir);
 Wl_ur_tang_vector = (Wl.v - Wl_ur_norm*norm_dir);
 if(Wl_ur_tang != ZERO){
    Wl_ur_tang_unit_vector =  Wl_ur_tang_vector/Wl_ur_tang;
 }else{
    Wl_ur_tang_unit_vector= Vector3D_ZERO;
 }
 
 Wl_rotated.rho = Wl.rho;
 Wl_rotated.v.x = Wl_ur_norm ;
 Wl_rotated.v.y = Wl_ur_tang;
 Wl_rotated.v.z = ZERO;
 Wl_rotated.p = Wl.p;
 // Right state velocity in rotated frame
 Wr_ur_norm = dot(Wr.v, norm_dir);
 Wr_ur_tang_vector = Wr.v - Wr_ur_norm*norm_dir;
 Wr_ur_tang = abs(Wr.v - Wr_ur_norm*norm_dir);
 if( Wr_ur_tang != ZERO){
    Wr_ur_tang_unit_vector =  Wr_ur_tang_vector/Wr_ur_tang ;
 }else{
    Wr_ur_tang_unit_vector= Vector3D_ZERO;  
 }
 
 Wr_rotated.rho = Wr.rho;
 Wr_rotated.v.x = Wr_ur_norm;
 Wr_rotated.v.y = dot( Wr_ur_tang_vector, Wl_ur_tang_unit_vector);
 Wr_rotated.v.z = abs( Wr_ur_tang_vector -Wr_rotated.v.y* Wl_ur_tang_unit_vector);
 Wr_rotated.p = Wr.p;


 Wr_ur_tang_z = abs(Wr_ur_tang_vector-Wr_rotated.v.y* Wl_ur_tang_unit_vector);
 Wr_ur_tang_z_vector = Wr_ur_tang_vector -Wr_rotated.v.y* Wl_ur_tang_unit_vector;
 if(Wr_ur_tang_z !=ZERO){
    Wr_ur_tang_z_unit_vector = Wr_ur_tang_z_vector /Wr_ur_tang_z ;
 }else{
    Wr_ur_tang_z_unit_vector = Vector3D_ZERO;
    
 }
 

    /* Evaluate the intermediate state solution 
       flux in the rotated frame. */

 Flux_rotated = FluxHLLE_x(Wl_rotated, Wr_rotated);
  
   /* Rotate back to the original Cartesian reference
      frame and return the solution flux. */
  
 Flux.Copy(Flux_rotated);


 
 Flux_rotated_x = Flux.rhov.x*norm_dir;
 Flux_rotated_tang_y = Flux.rhov.y* Wl_ur_tang_unit_vector ;
 Flux_rotated_tang_z = Flux.rhov.z* Wr_ur_tang_z_unit_vector;
 
      
 Flux.rhov =  Flux_rotated_x + Flux_rotated_tang_y+ Flux_rotated_tang_z;
 
 Flux.zero_non_sol();
 

 return (Flux);

}

Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::FluxHLLE_n(
   const Euler3D_ThermallyPerfect_cState &Ul,
   const Euler3D_ThermallyPerfect_cState &Ur,
   const Vector3D &norm_dir) {
   return (FluxHLLE_n(Ul.W(), Ur.W(), norm_dir));
}

/********************************************************
 * Routine: HartenFixPos (Harten Entropy Fix)           *
 *                                                      *
 * This function returns the positive parts of the      *
 * corrected elemental wave speeds or eigenvalues       *
 * according to the entropy fix of Harten (1983).       *
 *                                                       *
 ********************************************************/
Euler3D_ThermallyPerfect_pState HartenFixPos(
   const Euler3D_ThermallyPerfect_pState &lambdas_a,
   const Euler3D_ThermallyPerfect_pState &lambdas_l,
   const Euler3D_ThermallyPerfect_pState &lambdas_r) {
   
   Euler3D_ThermallyPerfect_pState NEW;
   NEW.rho = HartenFixPos(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
   NEW.v.x = HALF*(lambdas_a[2]+fabs(lambdas_a[2]));
   NEW.v.y = HALF*(lambdas_a[3]+fabs(lambdas_a[3]));
   NEW.v.z = HALF*(lambdas_a[4]+fabs(lambdas_a[4]));

   NEW.p = HartenFixPos(lambdas_a[5],lambdas_l[5],lambdas_r[5]);
  
   for( int i=NUM_EULER3D_VAR_SANS_SPECIES+1; i<=NEW.num_vars; i++){
      NEW.spec[i-(NUM_EULER3D_VAR_SANS_SPECIES+1)].c = 
         HALF*(lambdas_a[i]+fabs(lambdas_a[i]));
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
Euler3D_ThermallyPerfect_pState HartenFixNeg(
   const Euler3D_ThermallyPerfect_pState &lambdas_a,
   const Euler3D_ThermallyPerfect_pState &lambdas_l,
   const Euler3D_ThermallyPerfect_pState &lambdas_r) {
  
   Euler3D_ThermallyPerfect_pState NEW;
   
   NEW.rho = HartenFixNeg(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
   NEW.v.x = HALF*(lambdas_a[2]-fabs(lambdas_a[2]));
   NEW.v.y = HALF*(lambdas_a[3]-fabs(lambdas_a[3]));
   NEW.v.z = HALF*(lambdas_a[4]-fabs(lambdas_a[4]));
   NEW.p = HartenFixNeg(lambdas_a[5],lambdas_l[5],lambdas_r[5]);
   
   for( int i=NUM_EULER3D_VAR_SANS_SPECIES+1; i<=NEW.num_vars; i++){
      NEW.spec[i-(NUM_EULER3D_VAR_SANS_SPECIES+1)].c = 
         HALF*(lambdas_a[i]-fabs(lambdas_a[i]));
   }
   
   return (NEW);
}

// Flux Roe -- based on Harten fix 
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::FluxRoe_x(
   const  Euler3D_ThermallyPerfect_pState &Wl,  
   const  Euler3D_ThermallyPerfect_pState &Wr){
   
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
      
      wavespeeds = HartenFixNeg(lambdas_a,
                                lambdas_l,
                                lambdas_r);
      
      
      
       
      for (int i=1 ; i <= Wl.num_vars; i++) {
         if (wavespeeds[i] < ZERO) {
            Flux += wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
            
            
         }
      } 
   } else {
      
      Flux = Wr.F();
      
      wavespeeds = HartenFixPos(lambdas_a,
                                lambdas_l,
                                lambdas_r);
      
      for (int i=1; i <= Wl.num_vars; i++) {
         if (wavespeeds[i] > ZERO) {
            Flux -= wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i); 
            
            
         }
      } 
   } 
    
   /* Return solution flux. */    
   return (Flux);    
    
   
}
   
Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::FluxRoe_n(
   const Euler3D_ThermallyPerfect_pState &Wl,
   const Euler3D_ThermallyPerfect_pState &Wr,
   const Vector3D &norm_dir){
 
 double Wl_ur_norm, Wl_ur_tang;
 double Wr_ur_norm, Wr_ur_tang ;
 double Wr_ur_tang_z;
 
 Vector3D Flux_rotated_x, Flux_rotated_tang_y, Flux_rotated_tang_z ;
 Vector3D Wl_ur_tang_vector, Wr_ur_tang_vector;
 Vector3D Wl_ur_tang_unit_vector, Wr_ur_tang_unit_vector;
 Vector3D Wr_ur_tang_z_vector, Wr_ur_tang_z_unit_vector;
 
 //solnvec in  Wl_rotated (Wr_rotated) is allocated using new 
 Euler3D_ThermallyPerfect_pState Wl_rotated, Wr_rotated;
 Euler3D_ThermallyPerfect_cState Flux, Flux_rotated;
 
 
 /* Apply the frame rotation and evaluate left and right
    solution states in the local rotated frame defined
    by the unit normal vector. */
 Wl_rotated.Copy(Wl);
 Wr_rotated.Copy(Wr);
 
 // Left state velocity in rotated frame
 Wl_ur_norm = dot(Wl.v, norm_dir);
 Wl_ur_tang = abs(Wl.v - Wl_ur_norm*norm_dir);
 Wl_ur_tang_vector = (Wl.v - Wl_ur_norm*norm_dir);
 if(Wl_ur_tang != ZERO){
    Wl_ur_tang_unit_vector =  Wl_ur_tang_vector/Wl_ur_tang;
 }else{
    Wl_ur_tang_unit_vector= Vector3D_ZERO;
 }
 
 Wl_rotated.rho = Wl.rho;
 Wl_rotated.v.x = Wl_ur_norm ;
 Wl_rotated.v.y = Wl_ur_tang;
 Wl_rotated.v.z = ZERO;
 Wl_rotated.p = Wl.p;
 // Right state velocity in rotated frame
 Wr_ur_norm = dot(Wr.v, norm_dir);
 Wr_ur_tang_vector = Wr.v - Wr_ur_norm*norm_dir;
 Wr_ur_tang = abs(Wr.v - Wr_ur_norm*norm_dir);
 if( Wr_ur_tang != ZERO){
    Wr_ur_tang_unit_vector =  Wr_ur_tang_vector/Wr_ur_tang ;
 }else{
    Wr_ur_tang_unit_vector= Vector3D_ZERO;  
 }
 
 Wr_rotated.rho = Wr.rho;
 Wr_rotated.v.x = Wr_ur_norm;
 Wr_rotated.v.y = dot( Wr_ur_tang_vector, Wl_ur_tang_unit_vector);
 Wr_rotated.v.z = abs( Wr_ur_tang_vector -Wr_rotated.v.y* Wl_ur_tang_unit_vector);
 Wr_rotated.p = Wr.p;


 Wr_ur_tang_z = abs(Wr_ur_tang_vector-Wr_rotated.v.y* Wl_ur_tang_unit_vector);
 Wr_ur_tang_z_vector = Wr_ur_tang_vector -Wr_rotated.v.y* Wl_ur_tang_unit_vector;
 if(Wr_ur_tang_z !=ZERO){
    Wr_ur_tang_z_unit_vector = Wr_ur_tang_z_vector /Wr_ur_tang_z ;
 }else{
    Wr_ur_tang_z_unit_vector = Vector3D_ZERO;
    
 }
 

    /* Evaluate the intermediate state solution 
       flux in the rotated frame. */

 Flux_rotated = FluxRoe_x(Wl_rotated, Wr_rotated);
  
   /* Rotate back to the original Cartesian reference
      frame and return the solution flux. */
  
 Flux.Copy(Flux_rotated);


 
 Flux_rotated_x = Flux.rhov.x*norm_dir;
 Flux_rotated_tang_y = Flux.rhov.y* Wl_ur_tang_unit_vector ;
 Flux_rotated_tang_z = Flux.rhov.z* Wr_ur_tang_z_unit_vector;
 
      
 Flux.rhov =  Flux_rotated_x + Flux_rotated_tang_y+ Flux_rotated_tang_z;
 
 Flux.zero_non_sol();
 
 return (Flux);

}

/********************************************************
 * Routine: Reflect                                     *
 *                                                      *
 * This function returns the reflected solution state   *
 * in a given direction given the primitive solution    *
 * variables and the unit normal vector in the          *
 * direction of interest.                               *
 *                                                      *
 ********************************************************/
Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::Reflect(
   const Euler3D_ThermallyPerfect_pState &W,
   const Vector3D &norm_dir) {

   Vector3D ur_norm, ur_tang, vr_tot;
   
   Euler3D_ThermallyPerfect_pState Temp;
   Temp.Copy(W);
   
   ur_norm = dot(W.v, norm_dir)*norm_dir;
   ur_tang = W.v - ur_norm;
   
   ur_norm = -ur_norm;
   vr_tot = ur_norm + ur_tang;
   
   Temp.v = vr_tot;
 
    return (Temp);
       
}

Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::Moving_Wall(
   const Euler3D_ThermallyPerfect_pState &Win,
   const Euler3D_ThermallyPerfect_pState &Wout,
   const Vector3D &norm_dir, 
   const Vector3D &wall_velocity,
   const Vector3D &pressure_gradient,
   const int &TEMPERATURE_BC_FLAG) {

   Euler3D_ThermallyPerfect_pState Temp;
   Temp.Copy(Win);
   
   if(wall_velocity == Vector3D_ZERO){
      Temp.v = -Win.v;
      
   }else{
      
      double  Wall_velocity_tang ;
      Vector3D ur_norm, ur_tang, vr_tot, uw_tang;
      
      
      ur_norm = dot(Win.v, norm_dir)*norm_dir;
      ur_tang = Win.v - ur_norm;
      
      
      uw_tang = wall_velocity - dot(norm_dir,wall_velocity)*norm_dir;
      
      ur_norm = -ur_norm;
      ur_tang = 2.0*uw_tang - ur_tang;
      vr_tot = ur_norm + ur_tang;
      
      Temp.v = vr_tot;
   
      
   }
   
   //  Fixed Wall Temperature or constant extrapolation for Adiabatic
   if(TEMPERATURE_BC_FLAG == FIXED_TEMPERATURE_WALL){
      if (pressure_gradient != Vector3D_ZERO){
         Temp.rho = Wout.p/(Temp.Rtot()*Wout.T());
      }else{
         
         Temp.rho = Win.p/(Temp.Rtot()*Wout.T());
      }
      
   }
   
   return (Temp);

}

Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::No_Slip(
   const Euler3D_ThermallyPerfect_pState &Win,
   const Euler3D_ThermallyPerfect_pState &Wout,
   const Vector3D &norm_dir,
   const Vector3D &pressure_gradient,
   const int &TEMPERATURE_BC_FLAG) {
   
   return (Moving_Wall(Win, Wout, norm_dir,Vector3D_ZERO,pressure_gradient,TEMPERATURE_BC_FLAG));
   
}
