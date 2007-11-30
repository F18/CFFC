/********************** FANS3DThermallyPerfectState.cc **************
  This file defines the various member functions of the 
  FANS3D thermally perfect gaseous mixture class.

   assosicated files:
           FANS3DThermallyPerfectState.h    

*********************************************************************/

// Include required CFFC header files

#ifndef _FANS3D_THERMALLYPERFECT_STATE_INCLUDED
#include "FANS3DThermallyPerfectState.h"
#endif // _FANS3D_THERMALLYPERFECT_STATE_INCLUDED   

/****************************************************
  Speed of sound using 
  a^2 = dip/dirho + p/rho^2( die/dip)^-1
  from eigenvalue analysis using e =f(p,rho)
****************************************************/
double FANS3D_ThermallyPerfect_KOmega_pState::a(void){
   double sum;
   sum = sqr(Euler3D_ThermallyPerfect_pState::a());
   sum += 2.0/3.0*k*Euler3D_ThermallyPerfect_pState::g();
   return sqrt(sum);
}

double FANS3D_ThermallyPerfect_KOmega_pState::a(void) const{
   double sum;
   sum = sqr(Euler3D_ThermallyPerfect_pState::a());
   sum += 2.0/3.0*k*Euler3D_ThermallyPerfect_pState::g();
   return sqrt(sum);
}

void FANS3D_ThermallyPerfect_KOmega_pState::set_species_data(const int &n,
                                                             const string *S,
                                                             const char *PATH,
                                                             const int &debug, 
                                                             const double &Mr, 
                                                             const double* Sc,
                                                             const int &trans_data){ 

   //Deallocate_static();
   Deallocate(); //Clean up memory before changing ns
   ns = n;
   
   num_vars = ns + NUM_EULER3D_VAR_SANS_SPECIES+ NUM_FANS3D_VAR_EXTRA;
   
   //read in NASA data for each species to be used 
   specdata = new NASARP1311data[ns]; 
   Schmidt = new double[ns];
   
   for(int i=0; i<ns; i++){
      specdata[i].Getdata(S[i], PATH, trans_data);  
      Schmidt[i] = Sc[i];
   } 

   //set data temperature ranges for mixture
   Temp_low_range(); 
   Temp_high_range(); 
   
   //Set Debug Information level
   debug_level = debug;
   
   //setup initial array for mass fractions

   spec = new Species[ns];
   for(int i=0; i<ns; i++){
      spec[i].c = ONE/ns ; 
   }
}   

void FANS3D_ThermallyPerfect_KOmega_cState::set_species_data(const int &n, 
                                                             const string *S, 
                                                             const char *PATH,
                                                             const int &debug, 
                                                             const double &Mr, 
                                                             const double* Sc,
                                                             const int &trans_data){ 
 
   //Deallocate_static();
   Deallocate(); //Clean up memory before changing ns
 
   ns =n; 
   num_vars = ns + NUM_EULER3D_VAR_SANS_SPECIES+ NUM_FANS3D_VAR_EXTRA;

   //read in NASA data for each species to be used
   specdata = new NASARP1311data[ns];
   Schmidt = new double[ns];

   for(int i=0; i<ns; i++){
      //overwrite default data  
     specdata[i].Getdata(S[i], PATH, trans_data);
      Schmidt[i] = Sc[i];  
   }  

   //set data temperature ranges for mixture
   Temp_low_range();
   Temp_high_range();

   //Set Debug Information level
   debug_level = debug;

   //setup initial array for mass fractions
   rhospec = new Species[ns];
   for(int i=0; i<ns; i++){
      rhospec[i].c = rho/ns; 
   }
   
}      

/********************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::F -- 
Inviscid flux (x-direction).   *
*********************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::F(void){
   
   // call the inviscid flux calculation in Euler3D_ThermallyPerfect_pState class  
   Euler3D_ThermallyPerfect_cState Temp;
   Temp = Euler3D_ThermallyPerfect_pState::F();
   
   Temp.rhov.x += 2.0/3.0*rho*k;
   
   Temp.E += v.x*(2.0/3.0*rho*k);
   
   double Temp_rhok = rho*v.x*k;
   double Temp_rhoomega = rho*v.x*omega;
   
   return (FANS3D_ThermallyPerfect_KOmega_cState(Temp.rho, 
                                                 Temp.rhov, 
                                                 Temp.E, 
                                                 Temp_rhok, 
                                                 Temp_rhoomega,
                                                 Temp.rhospec));
}

FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::F(void) const{
   
   // call the inviscid flux calculation in Euler3D_ThermallyPerfect_pState class  
   Euler3D_ThermallyPerfect_cState Temp;
   Temp = Euler3D_ThermallyPerfect_pState::F();
  
   Temp.rhov.x += 2.0/3.0*rho*k;
   
   Temp.E += v.x*(2.0/3.0*rho*k);
   
   double Temp_rhok = rho*v.x*k;
   double Temp_rhoomega = rho*v.x*omega;
 
   return (FANS3D_ThermallyPerfect_KOmega_cState(
              Temp.rho, Temp.rhov, Temp.E, Temp_rhok, Temp_rhoomega,Temp.rhospec ));
}

/********************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::Fv    -- 
 Viscous flux (x-direction).   
 ********************************************************/

FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::Fv(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx,
                                                                                const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                                                                                const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) const{
   
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.Vacuum();
   
   double Temperature = T();
   double mu_t = eddy_viscosity();
   
   Tensor3D molecular_stress, Reynolds_stress;
   Vector3D heat_flux;
   
   molecular_stress = tau(dWdx, dWdy, dWdz);
   Reynolds_stress = lambda(dWdx, dWdy, dWdz);
   
   /*----------------- Thermal Diffusion ---------------------------*/
   for (int index = 0; index<ns; ++index)
     spec[index].diffusion_coef = (mu()/rho)/Schmidt[index];
   
   heat_flux = qflux(dWdx, dWdy, dWdz) + qflux_t(dWdx, dWdy, dWdz);
   
   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux -=  rho*thermal_diffusion( );
   
   Temp.rhov.x = molecular_stress.xx + Reynolds_stress.xx + 2.0/3.0*rho*k;
   Temp.rhov.y = molecular_stress.xy + Reynolds_stress.xy;
   Temp.rhov.z = molecular_stress.xz + Reynolds_stress.xz;
   Temp.E = v.x*(molecular_stress.xx + Reynolds_stress.xx + 2.0/3.0*rho*k) +
      v.y*(molecular_stress.xy + Reynolds_stress.xy)  
      + v.z*(molecular_stress.xz + Reynolds_stress.xz) - heat_flux.x 
      +(mu() + mu_t*k_omega_model.sigma_star)*dWdx.k;
   Temp.rhok = (mu()+mu_t*k_omega_model.sigma_star)*dWdx.k;
   Temp.rhoomega = (mu()+mu_t*k_omega_model.sigma)*dWdx.omega;
   
   // species transport 
   for (int index = 0; index<ns; ++index){
      Temp.rhospec[index] = rho*(spec[index].diffusion_coef+Dm_t()) *dWdx.spec[index].c;
   }
   return (Temp);
}


FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::Gv(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx,
                                                                                const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                                                                                const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) const{

   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.Vacuum();
   
   double Temperature = T();
   double mu_t = eddy_viscosity();
   
   Tensor3D molecular_stress, Reynolds_stress;
   Vector3D heat_flux;
   
   molecular_stress = tau(dWdx, dWdy, dWdz);
   Reynolds_stress = lambda(dWdx, dWdy, dWdz);
   
   /*----------------- Thermal Diffusion ---------------------------*/
   for (int index = 0; index<ns; ++index)
     spec[index].diffusion_coef = (mu()/rho)/Schmidt[index];
   
   heat_flux = qflux(dWdx, dWdy, dWdz) + qflux_t(dWdx, dWdy, dWdz);
   
   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux -=  rho*thermal_diffusion( );
   
   Temp.rhov.x = molecular_stress.xy + Reynolds_stress.xy;
   Temp.rhov.y = molecular_stress.yy + Reynolds_stress.yy + 2.0/3.0*rho*k ;
   Temp.rhov.z = molecular_stress.yz + Reynolds_stress.yz;
   Temp.E = v.x*(molecular_stress.xy + Reynolds_stress.xy) +
      v.y*(molecular_stress.yy + Reynolds_stress.yy + 2.0/3.0*rho*k)  
      + v.z*(molecular_stress.yz + Reynolds_stress.yz) - heat_flux.y 
      +(mu() + mu_t*k_omega_model.sigma_star)*dWdy.k;
   Temp.rhok = (mu()+mu_t*k_omega_model.sigma_star)*dWdy.k;
   Temp.rhoomega = (mu()+mu_t*k_omega_model.sigma)*dWdy.omega;
   
   // species transport 
   for (int index = 0; index<ns; ++index){
      Temp.rhospec[index] = rho*(spec[index].diffusion_coef+Dm_t()) *dWdy.spec[index].c;
   }
   return (Temp);
}

FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::Hv(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx,
                                                                                const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                                                                                const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) const{

   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.Vacuum();
   
   double Temperature = T();
   double mu_t = eddy_viscosity();
   
   Tensor3D molecular_stress, Reynolds_stress;
   Vector3D heat_flux;
   
   molecular_stress.zero();
   Reynolds_stress.zero();
   heat_flux.zero();
   molecular_stress = tau(dWdx, dWdy, dWdz);
   Reynolds_stress = lambda(dWdx, dWdy, dWdz);
   
   /*----------------- Thermal Diffusion ---------------------------*/
   for (int index = 0; index <ns; ++index)
     spec[index].diffusion_coef = (mu()/rho)/Schmidt[index];
   
   heat_flux = qflux(dWdx, dWdy, dWdz) + qflux_t(dWdx, dWdy, dWdz);
   
   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux -=  rho*thermal_diffusion();
   
   Temp.rhov.x = molecular_stress.xz + Reynolds_stress.xz;
   Temp.rhov.y = molecular_stress.yz + Reynolds_stress.yz;
   Temp.rhov.z = molecular_stress.zz + Reynolds_stress.zz + 2.0/3.0*rho*k;
   Temp.E = v.x*(molecular_stress.xz + Reynolds_stress.xz) +
      v.y*(molecular_stress.yz + Reynolds_stress.yz)  
      + v.z*(molecular_stress.zz + Reynolds_stress.zz + 2.0/3.0*rho*k) - heat_flux.z 
      +(mu() + mu_t*k_omega_model.sigma_star)*dWdz.k;
   Temp.rhok = (mu()+mu_t*k_omega_model.sigma_star)*dWdz.k;
   Temp.rhoomega = (mu()+mu_t*k_omega_model.sigma)*dWdz.omega;
   
   // species transport 
   for (int index = 0; index<ns; ++index){
      Temp.rhospec[index] = rho*(spec[index].diffusion_coef+Dm_t()) *dWdz.spec[index].c;
   }
   return (Temp);
}


Tensor3D FANS3D_ThermallyPerfect_KOmega_pState::rotation_tensor(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                                                                const FANS3D_ThermallyPerfect_KOmega_pState &dWdy, 
                                                                const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) const{
   
   Tensor3D rotation;
   rotation.zero();
   
   rotation.xy = -1.0/2.0*(dWdx.v.y - dWdy.v.x);
   rotation.xz = 1.0/2.0*(dWdz.v.x - dWdx.v.z);
   rotation.yz = -1.0/2.0*(dWdy.v.z - dWdz.v.y);
     
   return (rotation);
}

Tensor3D FANS3D_ThermallyPerfect_KOmega_pState::strain_rate(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                                                            const FANS3D_ThermallyPerfect_KOmega_pState &dWdy, 
                                                            const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) const{
 
   Tensor3D strain_rate;
   strain_rate.zero();
   
   strain_rate.xx = dWdx.v.x;
   strain_rate.xy = 1.0/2.0*(dWdy.v.x + dWdx.v.y);
   strain_rate.yy = dWdy.v.y;
   strain_rate.yz = 1.0/2.0*(dWdz.v.y + dWdy.v.z);
   strain_rate.zz = dWdz.v.z;
   strain_rate.xz = 1.0/2.0*(dWdx.v.z + dWdz.v.x);
   
   return (strain_rate);
}

double FANS3D_ThermallyPerfect_KOmega_cState::a(void) const{
   double sum;
   
   sum = sqr(Euler3D_ThermallyPerfect_cState::a());
   
   sum += 2.0/3.0*rhok/rho*Euler3D_ThermallyPerfect_cState::g();
   
   return sqrt(sum);
}
  
/***************** EIGENVALUES *************************************
 *******************************************************************/

/************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState::lambda -- Eigenvalue(s) (x-direction).    *
 ************************************************************************************/
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::lambda_x(void){
  
   double c = a();
   Euler3D_ThermallyPerfect_pState Temp;
   
   Temp = Euler3D_ThermallyPerfect_pState::lambda_x();
 
   return (FANS3D_ThermallyPerfect_KOmega_pState(v.x-c, Temp.v, v.x+c, v.x, v.x,Temp.spec));
   
}

FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::lambda_x(void) const {
  
   double c = a();
   Euler3D_ThermallyPerfect_pState Temp;
   
   Temp = Euler3D_ThermallyPerfect_pState::lambda_x();
   
   return (FANS3D_ThermallyPerfect_KOmega_pState(v.x-c, Temp.v, v.x+c, v.x, v.x, Temp.spec));
   
}

/********************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState -- Binary arithmetic operators.        *
 ********************************************************************************/
//------------------ Addition ------------------------//
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::operator +(const FANS3D_ThermallyPerfect_KOmega_pState &W) const{ 
   
   FANS3D_ThermallyPerfect_KOmega_pState Temp;
   Temp.Copy(*this);
   Temp += W;
   return Temp;
}

//------------------ Subtraction ------------------------//
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::operator -(const FANS3D_ThermallyPerfect_KOmega_pState &W) const{
   FANS3D_ThermallyPerfect_KOmega_pState Temp;
   Temp.Copy(*this);
   Temp -= W;
   return Temp;
}

//---------------- Scalar Multiplication ------------------//
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::operator *(const double &a) const{
   FANS3D_ThermallyPerfect_KOmega_pState Temp;
   Temp.Copy(*this);
   Temp.rho = rho*a;  Temp.v = v*a; Temp.p = p*a;
   Temp.k = k*a; Temp.omega = omega*a;
   for( int i=0; i<ns; i++){
      Temp.spec[i] = spec[i]*a;
   } 
   return(Temp);
}

FANS3D_ThermallyPerfect_KOmega_pState operator *(const double &a, const FANS3D_ThermallyPerfect_KOmega_pState &W){
   FANS3D_ThermallyPerfect_KOmega_pState Temp;
   //Temp.Copy(W);
   Temp.rho = W.rho*a;  Temp.v = W.v*a; Temp.p = W.p*a;
   Temp.k = W.k*a; Temp.omega = W.omega*a;
   for( int i=0; i<W.ns; i++){
      Temp.spec[i] = W.spec[i]*a;
   } 
   return(Temp);
}

//--------------- Scalar Division ------------------------//
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::operator /(const double &a) const {
   FANS3D_ThermallyPerfect_KOmega_pState Temp(rho,v,p, k, omega);
   Temp.Copy(*this);
   Temp.rho = rho/a; Temp.v = v/a; Temp.p = p/a; 
   Temp.k = k/a; Temp.omega = omega/a;
   for(int i=0; i<ns; i++){
      Temp.spec[i] = spec[i]/a; 
   } 
   return(Temp);
}

//----------------- Inner Product ------------------------//
double FANS3D_ThermallyPerfect_KOmega_pState::operator *(const FANS3D_ThermallyPerfect_KOmega_pState &W) const{
   double sum=0.0;
   for(int i=0; i<ns; i++){
      sum += spec[i]*W.spec[i];
   }  
   return (rho*W.rho + v*W.v + p*W.p + k*W.k + omega*W.omega + sum);
}

//----------- solution state product operator ------------//
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::operator ^(const FANS3D_ThermallyPerfect_KOmega_pState &W) const{
   FANS3D_ThermallyPerfect_KOmega_pState Temp(rho,v,p, k, omega);
   Temp.Copy(*this);
   Temp.rho = rho*W.rho;
   Temp.v.x = v.x*W.v.x;
   Temp.v.y = v.y*W.v.y;
   Temp.v.z = v.z*W.v.z;
   Temp.p = p*W.p;
   Temp.k = k*W.k;
   Temp.omega = omega*W.omega;
   for(int i=0; i<ns; i++){
      Temp.spec[i] = spec[i]*W.spec[i];
   }  
   return(Temp);
}

//----------------- Assignment ----------------------------//
FANS3D_ThermallyPerfect_KOmega_pState& FANS3D_ThermallyPerfect_KOmega_pState::operator =(const FANS3D_ThermallyPerfect_KOmega_pState &W){
   //self assignment protection
   if (this != &W){   
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
      }   
      else {
         cerr<<"\n Mismatch in number of species \n ";
      }  
   }
   return (*this);
}

/*******************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState -- Shortcut arithmetic operators.     *
 *******************************************************************************/
FANS3D_ThermallyPerfect_KOmega_pState& FANS3D_ThermallyPerfect_KOmega_pState::operator +=(const FANS3D_ThermallyPerfect_KOmega_pState &W){
   rho += W.rho;
   v += W.v; 
   p += W.p; 
   k += W.k;
   omega += W.omega;
   for( int i=0; i<ns; i++){
      spec[i] += W.spec[i];
   } 
   return (*this);
}

FANS3D_ThermallyPerfect_KOmega_pState& FANS3D_ThermallyPerfect_KOmega_pState::operator -=(const FANS3D_ThermallyPerfect_KOmega_pState &W) {
   rho -= W.rho;
   v -= W.v;
   p -= W.p;
   k -= W.k;
   omega -= W.omega;
   for(int i=0; i<ns; i++){
      spec[i] -= W.spec[i];
   }  
   return (*this); 
}

/********************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState -- Relational operators.               *
 ********************************************************************************/
int operator ==(const FANS3D_ThermallyPerfect_KOmega_pState &W1, 
                const FANS3D_ThermallyPerfect_KOmega_pState &W2) {
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
                 && W1.k == W2.k && W1.omega == W2.omega
                 && Temp == true);
      }
   } else {
      cerr<<"\n Mismatch in number of species \n";
      exit(1);
   }
}

int operator !=(const FANS3D_ThermallyPerfect_KOmega_pState &W1, 
                const FANS3D_ThermallyPerfect_KOmega_pState &W2) {
   if(W1.ns == W2.ns){ //check that species are equal
      bool Temp = true;
      for(int i=0; i<W1.ns; i++){
         if( W1.spec[i] != W2.spec[i] ){
            Temp = false;
            break;
         } 
         return (W1.rho != W2.rho || W1.v != W2.v || W1.p != W2.p 
                 || W1.k != W2.k || W1.omega != W2.omega
                 || Temp != true);
      }
   } else {
      cerr<<"\n Mismatch in number of species \n";
      exit(1);
   }
}

/* FANS3D_ThermallyPerfect_KOmega_pState -- Input-output operators.*/
ostream &operator << (ostream &out_file, 
                      const FANS3D_ThermallyPerfect_KOmega_pState &W) {
   out_file.precision(10);
   out_file.setf(ios::scientific);
   out_file << " " << W.rho  << " " << W.v.x << " " << W.v.y 
            << " " << W.v.z << " " << W.p << " "<<W.k<< " "<<W.omega;
   for( int i=0; i<W.ns; i++){
      out_file<<" "<<W.spec[i];
   }
   out_file.unsetf(ios::scientific);
   return (out_file);
}

istream &operator >> (istream &in_file, FANS3D_ThermallyPerfect_KOmega_pState &W) {
   in_file.setf(ios::skipws);
   in_file >> W.rho >> W.v.x >> W.v.y >> W.v.z >>W.p>>W.k>>W.omega;
   //W.set_initial_values();
   for( int i=0; i<W.ns; i++){
      in_file>>W.spec[i];
   }
   in_file.unsetf(ios::skipws);
   return (in_file);
}


//----------------- Addition -----------------------------//
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_cState::operator +(const FANS3D_ThermallyPerfect_KOmega_cState &U) const{ 
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.Copy(*this);
   Temp += U;
   return Temp;
}

//------------------ Subtraction ------------------------//
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_cState::operator -(const FANS3D_ThermallyPerfect_KOmega_cState &U) const{
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.Copy(*this);
   Temp -= U;
   return Temp;
}

//---------------- Scalar Multiplication ------------------//
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_cState::operator *(const double &a) const{
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.Copy(*this);
   Temp.rho = rho*a;  Temp.rhov = rhov*a; Temp.E = E*a;
   Temp.rhok = rhok*a; Temp.rhoomega = rhoomega*a;
   for( int i=0; i<ns; i++){
      Temp.rhospec[i] = rhospec[i]*a;
   } 
   return(Temp);
}

FANS3D_ThermallyPerfect_KOmega_cState operator *(const double &a, const FANS3D_ThermallyPerfect_KOmega_cState &U){
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.rho = U.rho*a;  Temp.rhov = U.rhov*a; Temp.E = U.E*a;
   Temp.rhok = U.rhok*a; Temp.rhoomega = U.rhoomega*a;
   for( int i=0; i<U.ns; i++){
      Temp.rhospec[i] = U.rhospec[i]*a;
   } 
   return(Temp);
}
//--------------- Scalar Division ------------------------//
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_cState::operator /(const double &a) const {
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.Copy(*this);
   Temp.rho = rho/a; Temp.rhov = rhov/a; Temp.E = E/a;
   Temp.rhok = rhok/a; Temp.rhoomega = rhoomega/a;
   for(int i=0; i<ns; i++){
      Temp.rhospec[i] = rhospec[i]/a; 
   } 
   return(Temp);
}

//----------------- Inner Product ------------------------//
double FANS3D_ThermallyPerfect_KOmega_cState::operator *(const FANS3D_ThermallyPerfect_KOmega_cState &U) const{
   double sum=0.0;
   for(int i=0; i<ns; i++){
      sum += rhospec[i]*U.rhospec[i];
   }  
   return (rho*U.rho + rhov*U.rhov + E*U.E + rhok*U.rhok + rhoomega*U.rhoomega + sum);
}

//----------- solution state product operator ------------//
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_cState::operator ^(const FANS3D_ThermallyPerfect_KOmega_cState &U) const {
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.Copy(*this);
   Temp.rho = rho*U.rho;
   Temp.rhov.x = rhov.x*U.rhov.x;
   Temp.rhov.y = rhov.y*U.rhov.y;
   Temp.rhov.z = rhov.z*U.rhov.z;
   Temp.E = E*U.E;
   Temp.rhok = rhok*U.rhok;
   Temp.rhoomega = rhoomega*U.rhoomega;
   for(int i=0; i<ns; i++){
      Temp.rhospec[i] = rhospec[i]*U.rhospec[i];
   }  
   return(Temp);
}

//----------------- Assignment ----------------------------//
FANS3D_ThermallyPerfect_KOmega_cState& FANS3D_ThermallyPerfect_KOmega_cState::operator =(
   const FANS3D_ThermallyPerfect_KOmega_cState &U){
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
      } else {
         cerr<<"\n Mismatch in number of species \n ";
      }  
   }
   return (*this);
}


/*******************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState -- Shortcut arithmetic operators.     *
 *******************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState& FANS3D_ThermallyPerfect_KOmega_cState::operator +=(const FANS3D_ThermallyPerfect_KOmega_cState &U){
   rho += U.rho;
   rhov += U.rhov; 
   E += U.E;
   rhok += U.rhok;
   rhoomega += U.rhoomega;
   for( int i=0; i<ns; i++){
      rhospec[i] += U.rhospec[i];
   } 
   return (*this);
}

FANS3D_ThermallyPerfect_KOmega_cState& FANS3D_ThermallyPerfect_KOmega_cState::operator -=(const FANS3D_ThermallyPerfect_KOmega_cState &U) {
   rho -= U.rho;
   rhov -= U.rhov;
   E -= U.E;
   rhok -= U.rhok;
   rhoomega -= U.rhoomega;
   for(int i=0; i<ns; i++){
     rhospec[i] -= U.rhospec[i];
   }  
   return (*this); 
}

/*******************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState -- Unary arithmetic operators.        *
 *******************************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState operator -(const FANS3D_ThermallyPerfect_KOmega_cState &U) {
   Species *spt= new Species[U.ns];
   for(int i=0; i<U.ns; i++){
      spt[i] = -U.rhospec[i]; 
   }  
   FANS3D_ThermallyPerfect_KOmega_cState Temp(
      -U.rho,-U.rhov,-U.E, -U.rhok, -U.rhoomega, spt);
   delete[] spt;
   return(Temp);
}

/********************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState -- Relational operators.              *
 ********************************************************/
int operator ==(const FANS3D_ThermallyPerfect_KOmega_cState &U1, 
                const FANS3D_ThermallyPerfect_KOmega_cState &U2) {
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
                 && U1.rhok == U2.rhok && U1.rhoomega == U2.rhoomega
                 &&Temp == true);
      }
   } else {
      cerr<<"\n Mismatch in number of species \n";
      exit(1);
   }
}

int operator !=(const FANS3D_ThermallyPerfect_KOmega_cState &U1, 
                const FANS3D_ThermallyPerfect_KOmega_cState &U2) {
   if(U1.ns == U2.ns){ //check that species are equal
      bool Temp = true;
      for(int i=0; i<U1.ns; i++){
         if( U1.rhospec[i] != U2.rhospec[i] ){
            Temp = false;
            break;
         } 
         return (U1.rho != U2.rho || U1.rhov != U2.rhov || U1.E != U2.E 
                 || U1.rhok != U2.rhok || U1.rhoomega !=U2.rhoomega
                 || Temp != true);
      }
   } else {
      cerr<<"\n Mismatch in number of species \n";
      exit(1);
   }
}

/*******************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState -- Input-output operators.            *
 *******************************************************************************/
ostream &operator << (ostream &out_file, const FANS3D_ThermallyPerfect_KOmega_cState &U){
  //out_file.precision(20);
  out_file.setf(ios::scientific);
  out_file << " " << U.rho  << " " << U.rhov.x << " " << U.rhov.y << " " 
           << U.rhov.z << " "<< U.E<< " "<<U.rhok<< " "<<U.rhoomega;
  for( int i=0; i<U.ns; i++){
     out_file<<" "<<U.rhospec[i];
  } 
  out_file.unsetf(ios::scientific);
  return (out_file);
}

istream &operator >> (istream &in_file, FANS3D_ThermallyPerfect_KOmega_cState &U){
   in_file.setf(ios::skipws);
   in_file >> U.rho >> U.rhov.x >> U.rhov.y >> U.rhov.z >> U.E >> U.rhok >> U.rhoomega;
  //U.set_initial_values();
   for( int i=0; i<U.ns; i++){
     in_file>>U.rhospec[i]; 
   } 
   in_file.unsetf(ios::skipws);
   return (in_file);
}

FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::rc_x(const int &index) const {
   switch(index){  
   case 1:
      return (FANS3D_ThermallyPerfect_KOmega_cState(
              ONE, v.x-a(), v.y, v.z, (H()+2.0/3.0*rho*k)/rho-v.x*a(), k, omega,spec));
   case 2:
      return (FANS3D_ThermallyPerfect_KOmega_cState(
              ONE, v.x, v.y, v.z, (H() + 2.0/3.0*rho*k)/rho-Cp()*T(), k, omega, spec)); 
   case 3:
      return (FANS3D_ThermallyPerfect_KOmega_cState(
              ZERO, ZERO, rho, ZERO, rho*v.y, ZERO, ZERO, ZERO));
   case 4:
      return (FANS3D_ThermallyPerfect_KOmega_cState(
              ZERO, ZERO, ZERO, rho, rho*v.z, ZERO, ZERO, ZERO));
   case 5:
    
      return (FANS3D_ThermallyPerfect_KOmega_cState(
              ONE, v.x+a(), v.y, v.z, (H()+2.0/3.0*rho*k)/rho+v.x*a(), k, omega, spec));
   case 6:
    
      return (FANS3D_ThermallyPerfect_KOmega_cState(
              ZERO, ZERO, ZERO, ZERO, rho*(ONE-TWO/THREE*(ONE/(g()-ONE))), rho, ZERO, ZERO));
   case 7:
      return (FANS3D_ThermallyPerfect_KOmega_cState(
              ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, rho, ZERO));
   default :
      FANS3D_ThermallyPerfect_KOmega_cState NEW(ZERO);
      double RTOT = Rtot();
      double TEMP = p/(rho*RTOT);      
      NEW.E = rho*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Enthalpy(TEMP) + 
                   specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Heatofform() 
                   - Cp(TEMP)*TEMP*specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Rs()/RTOT);
      double PHI =specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Enthalpy(TEMP)-
         specdata[num_vars-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Enthalpy(TEMP)-
         (Cp(TEMP) -RTOT)*TEMP*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Rs() 
                                - specdata[num_vars - (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Rs())/RTOT;
      NEW.rhospec[index - (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].c = rho*PHI; 
      return NEW;
   };
   
}
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::rc_x(const int &index) {
   switch(index){  
   case 1:
      return (FANS3D_ThermallyPerfect_KOmega_cState(
              ONE, v.x-a(), v.y, v.z, (H()+2.0/3.0*rho*k)/rho-v.x*a(), k, omega,spec));
   case 2:
      return (FANS3D_ThermallyPerfect_KOmega_cState(
                 ONE, v.x, v.y, v.z, (H() + 2.0/3.0*rho*k)/rho-Cp()*T(), k, omega, spec)); 
   case 3:
      return (FANS3D_ThermallyPerfect_KOmega_cState(
                 ZERO, ZERO, rho, ZERO, rho*v.y, ZERO, ZERO, ZERO));
   case 4:
      return (FANS3D_ThermallyPerfect_KOmega_cState(
              ZERO, ZERO, ZERO, rho, rho*v.z, ZERO, ZERO, ZERO));
   case 5:
    
      return (FANS3D_ThermallyPerfect_KOmega_cState(
              ONE, v.x+a(), v.y, v.z, (H()+2.0/3.0*rho*k)/rho+v.x*a(), k, omega, spec));
   case 6:
    
      return (FANS3D_ThermallyPerfect_KOmega_cState(
              ZERO, ZERO, ZERO, ZERO, rho*(ONE-TWO/THREE*(ONE/(g()-ONE))), rho, ZERO, ZERO));
   case 7:
    
      return (FANS3D_ThermallyPerfect_KOmega_cState(
              ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, rho, ZERO));
   default :
      FANS3D_ThermallyPerfect_KOmega_cState NEW(ZERO);
      double RTOT = Rtot();
      double TEMP = p/(rho*RTOT);      
      NEW.E = rho*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Enthalpy(TEMP) + specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Heatofform() 
                   - Cp(TEMP)*TEMP*specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Rs()/RTOT);
      double PHI =specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Enthalpy(TEMP)-specdata[num_vars- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Enthalpy(TEMP)-
         (Cp(TEMP) -RTOT)*TEMP*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Rs() - specdata[num_vars - (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].Rs())/RTOT;
      NEW.rhospec[index - (NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].c = rho*PHI; 
      return NEW;
   };
   
}

// Primitive Left Eigenvector -- (x-direction)
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::lp_x(const int &index) const {
   switch(index){  
   case 1:
      return (FANS3D_ThermallyPerfect_KOmega_pState(
               k/(THREE*a()*a()), -HALF*rho/a(), ZERO, ZERO, HALF/(a()*a()), rho/(THREE*a()*a()), ZERO, ZERO));
   case 2:
      return (FANS3D_ThermallyPerfect_KOmega_pState(
              ONE - TWO*k/(THREE*a()*a()), ZERO, ZERO, ZERO, -ONE/(a()*a()), -TWO*rho/(THREE*a()*a()), ZERO, ZERO));
   case 3 :
      return  (FANS3D_ThermallyPerfect_KOmega_pState(ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO, ZERO));
   case 4 :
      return  (FANS3D_ThermallyPerfect_KOmega_pState(ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO));
   case 5 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(
              k/(THREE*a()*a()), HALF*rho/a(), ZERO, ZERO, HALF/(a()*a()), rho/(THREE*a()*a()), ZERO, ZERO));
   case 6 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(
              ZERO, ZERO, ZERO, ZERO, ZERO, ONE,ZERO, ZERO));
   case 7 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(
              ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO));
   default :
      FANS3D_ThermallyPerfect_KOmega_pState NEW(ZERO);
      NEW.spec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].c = ONE;
      return NEW;
   };
   
}

FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::lp_x(const int &index) {
   switch(index){  
   case 1:
      return (FANS3D_ThermallyPerfect_KOmega_pState(
                k/(THREE*a()*a()), -HALF*rho/a(), ZERO, ZERO, HALF/(a()*a()), rho/(THREE*a()*a()), ZERO, ZERO));
   case 2:
      return (FANS3D_ThermallyPerfect_KOmega_pState(
                 ONE - TWO*k/(THREE*a()*a()), ZERO, ZERO, ZERO, -ONE/(a()*a()), -TWO*rho/(THREE*a()*a()), ZERO, ZERO));
   case 3 :
      return  (FANS3D_ThermallyPerfect_KOmega_pState(ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO, ZERO));
   case 4 :
      return  (FANS3D_ThermallyPerfect_KOmega_pState(ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO));
   case 5 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(
                 k/(THREE*a()*a()), HALF*rho/a(), ZERO, ZERO, HALF/(a()*a()), rho/(THREE*a()*a()), ZERO, ZERO));
   case 6 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(
                 ZERO, ZERO, ZERO, ZERO, ZERO, ONE,ZERO, ZERO));
   case 7 :
      return (FANS3D_ThermallyPerfect_KOmega_pState(
                 ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO));
   default :
      FANS3D_ThermallyPerfect_KOmega_pState NEW(ZERO);
      NEW.spec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].c = ONE;
      return NEW;
   };
}

double FANS3D_ThermallyPerfect_KOmega_cState::T(void) const{
   double T = ZERO;
   double RTOT = Rtot();  
   //--------- Initial Guess ------------------------------//
   //using a polytropic gas assumption with gamma@200;
   double Tguess = (gamma_guess() - ONE)*(E - HALF*rhov.sqr()/rho-rhok)/(rho*RTOT);
   //--------- global newtons method to get T ---------------//
   double A = (E - HALF*rhov*rhov/rho -rhok)/rho;
 
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
         cout<<"\nTemperature didn't converge in FANS3D_ThermallyPerfect_KOmega_cState::T(void)";
         cout<<" with polytopic Tguess "<<Tguess<<", or lower than Tmin "<<low_temp_range<<" using "<<T;
      }   
   }
   return T;
} 

double FANS3D_ThermallyPerfect_KOmega_pState::E(void) const {   
   double E_euler = Euler3D_ThermallyPerfect_pState::E();
   return (E_euler + rho*k);
}

double FANS3D_ThermallyPerfect_KOmega_pState::H(void) const{
   double H_euler = Euler3D_ThermallyPerfect_pState::H();
   return (H_euler + rho*k);
}

double FANS3D_ThermallyPerfect_KOmega_pState::Hs(void) const{
   // Hs = hs + velocity^2 +k;
   double H_S =  rho*(Euler3D_ThermallyPerfect_pState::hs() + HALF*v.sqr() +k );
   return (H_S);
}


double FANS3D_ThermallyPerfect_KOmega_pState::eddy_viscosity(void) const{
   return (rho*k/max(TOLER,omega));
}

double FANS3D_ThermallyPerfect_KOmega_pState::Pr_t(void) const{
   return (0.9); 
}

double FANS3D_ThermallyPerfect_KOmega_pState::Sc_t(void) const{
   return (1.0);
}

double FANS3D_ThermallyPerfect_KOmega_pState::kappa_t(void) const{
  return (eddy_viscosity()*Cp()/Pr_t()); 
}

double FANS3D_ThermallyPerfect_KOmega_pState::Dm_t(void) const{
   return (eddy_viscosity()/(rho*Sc_t()));
}     

double FANS3D_ThermallyPerfect_KOmega_cState::eddy_viscosity(void) const{
   return (rho*rhok/max(TOLER,rhoomega));
}

double FANS3D_ThermallyPerfect_KOmega_cState::Pr_t(void) const{
   return (0.9);
}

double FANS3D_ThermallyPerfect_KOmega_cState::Sc_t(void) const{
   return (1.0);
}

double FANS3D_ThermallyPerfect_KOmega_cState::Dm_t(void) const{
   return (eddy_viscosity()/(rho*Sc_t()));
}

double FANS3D_ThermallyPerfect_KOmega_cState::p() const{
   return (rho*Rtot()*T());
}

double FANS3D_ThermallyPerfect_KOmega_cState::k() const{
   return (rhok/rho);
}

double FANS3D_ThermallyPerfect_KOmega_cState::omega() const{
   return (rhoomega/rho);
}

/********************************************************
 * Routine: RoeAverage (Roe Averages)                   *
 *                                                      *
 * This function returns the Roe-averaged primitive     *
 * solution state given left and right primitive        *
 * solution variables.  See Roe (1981).                 *
 *                                                      *
 ********************************************************/
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::RoeAverage(const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
                                                                                        const FANS3D_ThermallyPerfect_KOmega_pState &Wr) {

   double Hl, Hr, srhol, srhor;
   double Ha, ha;
   FANS3D_ThermallyPerfect_KOmega_pState Temp;

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
   Temp.k = (srhol*Wl.k+srhor*Wr.k)/(srhol+srhor);
   Temp.omega = (srhol*Wl.omega+srhor*Wr.omega)/(srhol+srhor);

   for(int i=0; i<Wl.ns; i++){
      Temp.spec[i].c = (srhol*Wl.spec[i].c + srhor*Wr.spec[i].c)/(srhol+srhor);
   }

   //Temp.p = (srhol*Wl.p + srhor*Wr.p)/(srhol+srhor);
   Ha = (srhol*Hl+srhor*Hr)/(srhol+srhor);
   ha = Ha - HALF*(sqr(Temp.v.x)+sqr(Temp.v.y)+sqr(Temp.v.z)) - Temp.k;
   
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
FANS3D_ThermallyPerfect_KOmega_cState  FANS3D_ThermallyPerfect_KOmega_pState::FluxHLLE_x(const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
                                                                                         const FANS3D_ThermallyPerfect_KOmega_pState &Wr) {
   
   double wavespeed_l, wavespeed_r;
 
   FANS3D_ThermallyPerfect_KOmega_pState Wa, lambdas_l, lambdas_r, lambdas_a;
   FANS3D_ThermallyPerfect_KOmega_cState Flux, dUrl;
   
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
   wavespeed_r = max(lambdas_r[num_vars-NUM_FANS3D_VAR_EXTRA-lambdas_r.ns],
                     lambdas_a[num_vars-NUM_FANS3D_VAR_EXTRA-lambdas_a.ns]);
   
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

FANS3D_ThermallyPerfect_KOmega_cState  FANS3D_ThermallyPerfect_KOmega_pState::FluxHLLE_x(const FANS3D_ThermallyPerfect_KOmega_cState &Ul,
                                                                                         const FANS3D_ThermallyPerfect_KOmega_cState &Ur) {
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
FANS3D_ThermallyPerfect_KOmega_cState  FANS3D_ThermallyPerfect_KOmega_pState::FluxHLLE_n(const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
                                                                                         const FANS3D_ThermallyPerfect_KOmega_pState &Wr,
                                                                                         const Vector3D &norm_dir) {

   double Wl_ur_norm, Wl_ur_tang;
   double Wr_ur_norm, Wr_ur_tang ;
   double Wr_ur_tang_z;
   
   Vector3D Flux_rotated_x, Flux_rotated_tang_y, Flux_rotated_tang_z ;
   Vector3D Wl_ur_tang_vector, Wr_ur_tang_vector;
   Vector3D Wl_ur_tang_unit_vector, Wr_ur_tang_unit_vector;
   Vector3D Wr_ur_tang_z_vector, Wr_ur_tang_z_unit_vector;
   
   FANS3D_ThermallyPerfect_KOmega_pState Wl_rotated, Wr_rotated;
   FANS3D_ThermallyPerfect_KOmega_cState Flux, Flux_rotated;
   
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

FANS3D_ThermallyPerfect_KOmega_cState  FANS3D_ThermallyPerfect_KOmega_pState::FluxHLLE_n(const FANS3D_ThermallyPerfect_KOmega_cState &Ul,
	      	                                                                         const FANS3D_ThermallyPerfect_KOmega_cState &Ur,
                                                                                         const Vector3D &norm_dir) {
    return (FluxHLLE_n(Ul.W(), Ur.W(), norm_dir));
}

/********************************************************
 * Routine: HartenFixPos (Harten Entropy Fix)           *
 *                                                      *
 * This function returns the positive parts of the      *
 * corrected elemental wave speeds or eigenvalues       *
 * according to the entropy fix of Harten (1983).       *
 *                                                      *
 ********************************************************/
FANS3D_ThermallyPerfect_KOmega_pState HartenFixPos(const FANS3D_ThermallyPerfect_KOmega_pState &lambdas_a,
                                                   const FANS3D_ThermallyPerfect_KOmega_pState &lambdas_l,
                                                   const FANS3D_ThermallyPerfect_KOmega_pState &lambdas_r) {
   
   FANS3D_ThermallyPerfect_KOmega_pState NEW;
   NEW.rho = HartenFixPos(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
   NEW.v.x = HALF*(lambdas_a[2]+fabs(lambdas_a[2]));
   NEW.v.y = HALF*(lambdas_a[3]+fabs(lambdas_a[3]));
   NEW.v.z = HALF*(lambdas_a[4]+fabs(lambdas_a[4]));
   NEW.p = HartenFixPos(lambdas_a[5],lambdas_l[5],lambdas_r[5]);
   NEW.k = HALF*(lambdas_a[6] + fabs(lambdas_a[6]));
   NEW.omega = HALF*(lambdas_a[7] +fabs(lambdas_a[7]));
  
   for( int i=(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA +1); i<=NEW.num_vars; i++){
      NEW.spec[i-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].c = 
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
FANS3D_ThermallyPerfect_KOmega_pState HartenFixNeg(const FANS3D_ThermallyPerfect_KOmega_pState &lambdas_a,
                                                   const FANS3D_ThermallyPerfect_KOmega_pState &lambdas_l,
                                                   const FANS3D_ThermallyPerfect_KOmega_pState &lambdas_r) {
  
   FANS3D_ThermallyPerfect_KOmega_pState NEW;
   NEW.rho = HartenFixNeg(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
   NEW.v.x = HALF*(lambdas_a[2]-fabs(lambdas_a[2]));
   NEW.v.y = HALF*(lambdas_a[3]-fabs(lambdas_a[3]));
   NEW.v.z = HALF*(lambdas_a[4]-fabs(lambdas_a[4]));
   NEW.p = HartenFixNeg(lambdas_a[5],lambdas_l[5],lambdas_r[5]);
   NEW.k = HALF*(lambdas_a[6] -fabs(lambdas_a[6]));
   NEW.omega = HALF*(lambdas_a[7] - fabs(lambdas_a[7]));
   
   for( int i = (NUM_EULER3D_VAR_SANS_SPECIES  + NUM_FANS3D_VAR_EXTRA+1); i<=NEW.num_vars; i++){
      NEW.spec[i-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA+1)].c = 
         HALF*(lambdas_a[i]-fabs(lambdas_a[i]));
   }
   
   return (NEW);
}

// Flux Roe -- based on Harten fix 
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::FluxRoe_x(const  FANS3D_ThermallyPerfect_KOmega_pState &Wl,  
                                                                                        const  FANS3D_ThermallyPerfect_KOmega_pState &Wr){
   
   FANS3D_ThermallyPerfect_KOmega_pState Wa, dWrl, wavespeeds, 
                                         lambdas_l, lambdas_r, lambdas_a;
   FANS3D_ThermallyPerfect_KOmega_cState Flux;

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
   
FANS3D_ThermallyPerfect_KOmega_cState  FANS3D_ThermallyPerfect_KOmega_pState::FluxRoe_n(const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
                                                                                        const FANS3D_ThermallyPerfect_KOmega_pState &Wr,
                                                                                        const Vector3D &norm_dir){
   
   double Wl_ur_norm, Wl_ur_tang;
   double Wr_ur_norm, Wr_ur_tang ;
   double Wr_ur_tang_z;
   
   Vector3D Flux_rotated_x, Flux_rotated_tang_y, Flux_rotated_tang_z ;
   Vector3D Wl_ur_tang_vector, Wr_ur_tang_vector;
   Vector3D Wl_ur_tang_unit_vector, Wr_ur_tang_unit_vector;
   Vector3D Wr_ur_tang_z_vector, Wr_ur_tang_z_unit_vector;
   
   FANS3D_ThermallyPerfect_KOmega_pState Wl_rotated, Wr_rotated;
   FANS3D_ThermallyPerfect_KOmega_cState Flux, Flux_rotated;
   
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
FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::Reflect(const FANS3D_ThermallyPerfect_KOmega_pState &W,
                                                                                     const Vector3D &norm_dir) {
   
   Vector3D ur_norm, ur_tang, vr_tot;
   
   FANS3D_ThermallyPerfect_KOmega_pState Temp;
   Temp = W;
   
   ur_norm = dot(W.v, norm_dir)*norm_dir;
   ur_tang = W.v - ur_norm;
   
   ur_norm = -ur_norm;
   vr_tot = ur_norm + ur_tang;
   
   Temp.v = vr_tot;
   
   return (Temp);
       
}


FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::Moving_Wall(const FANS3D_ThermallyPerfect_KOmega_pState &Win,
                                                                                         const FANS3D_ThermallyPerfect_KOmega_pState &Wout,
                                                                                         const Vector3D &norm_dir, 
                                                                                         const Vector3D &wall_velocity,
                                                                                         const Vector3D &pressure_gradient,
                                                                                         const int &TEMPERATURE_BC_FLAG) {
   
   FANS3D_ThermallyPerfect_KOmega_pState Temp;
   Temp = Win;
   
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
  
  Temp.k = ZERO; 
  /* Fixed Wall Temperature or constant extrapolation for Adiabatic */
  if(TEMPERATURE_BC_FLAG == FIXED_TEMPERATURE_WALL){
     if (pressure_gradient != Vector3D_ZERO){
        Temp.rho = Wout.p/(Temp.Rtot()*Wout.T());
     }else{
        Temp.rho = Temp.p/(Temp.Rtot()*Wout.T());
     }
  }
  return (Temp);
}

FANS3D_ThermallyPerfect_KOmega_pState  FANS3D_ThermallyPerfect_KOmega_pState::No_Slip(
   const FANS3D_ThermallyPerfect_KOmega_pState &Win,
   const FANS3D_ThermallyPerfect_KOmega_pState &Wout,
   const Vector3D &norm_dir,
   const Vector3D &pressure_gradient,
   const int &TEMPERATURE_BC_FLAG) {
   
   return (Moving_Wall(Win, Wout, norm_dir, Vector3D_ZERO, pressure_gradient, TEMPERATURE_BC_FLAG));
   
}

// molecular stress tensor
Tensor3D FANS3D_ThermallyPerfect_KOmega_pState::lambda(
   const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
   const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
   const FANS3D_ThermallyPerfect_KOmega_pState &dWdz){
   
   Tensor3D Reynolds_stress;
   double mu_t = eddy_viscosity();
   
   Reynolds_stress.xx = 1.0/3.0 *mu_t*(4.0*dWdx.v.x - 2.0*dWdy.v.y - 2.0*dWdz.v.z) - 2.0/3.0*rho*k;
   Reynolds_stress.xy = mu_t*(dWdx.v.y + dWdy.v.x);
   Reynolds_stress.xz = mu_t*(dWdx.v.z + dWdz.v.x);
   Reynolds_stress.yy = 1.0/3.0 *mu_t*(4.0*dWdy.v.y - 2.0*dWdx.v.x - 2.0*dWdz.v.z) - 2.0/3.0*rho*k;
   Reynolds_stress.yz = mu_t*(dWdy.v.z + dWdz.v.y);
   Reynolds_stress.zz = 1.0/3.0 *mu_t*(4.0*dWdz.v.z - 2.0*dWdx.v.x - 2.0*dWdy.v.y) - 2.0/3.0*rho*k;
   
   return (Reynolds_stress);
}

Tensor3D FANS3D_ThermallyPerfect_KOmega_pState::lambda(
   const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
   const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
   const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) const {

   Tensor3D Reynolds_stress;
   double mu_t = eddy_viscosity();
   
   Reynolds_stress.xx = 1.0/3.0 *mu_t*(4.0*dWdx.v.x - 2.0*dWdy.v.y - 2.0*dWdz.v.z) - 2.0/3.0*rho*k;
   Reynolds_stress.xy = mu_t*(dWdx.v.y + dWdy.v.x);
   Reynolds_stress.xz = mu_t*(dWdx.v.z + dWdz.v.x); 
   Reynolds_stress.yy = 1.0/3.0 *mu_t*(4.0*dWdy.v.y - 2.0*dWdx.v.x - 2.0*dWdz.v.z) - 2.0/3.0*rho*k;
   Reynolds_stress.yz = mu_t*(dWdy.v.z + dWdz.v.y);
   Reynolds_stress.zz = 1.0/3.0 *mu_t*(4.0*dWdz.v.z - 2.0*dWdx.v.x - 2.0*dWdy.v.y) - 2.0/3.0*rho*k;
 
   return (Reynolds_stress);
}

// heat flux: Fourier's law of heat conduction
Vector3D FANS3D_ThermallyPerfect_KOmega_pState::qflux_t(
    const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
    const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
    const FANS3D_ThermallyPerfect_KOmega_pState &dWdz){
  
    double Rmix = Rtot();
    
    Vector3D heat_flux, gradient_T;
    
    /* Temperature gradients from using the chain rule 
       with the ideal gas law (P=rho*R*T) 
       dT/dx = 1/(rho*R) *( dP/dx - P/rho * drho/dx) */
     
    gradient_T.x = (1.0/(rho*Rmix)) * (dWdx.p -(p/rho)*dWdx.rho);
    gradient_T.y = (1.0/(rho*Rmix)) * (dWdy.p -(p/rho)*dWdy.rho);
    gradient_T.z = (1.0/(rho*Rmix)) * (dWdz.p -(p/rho)*dWdz.rho);
    
    heat_flux = -kappa_t()*gradient_T;
        
    return (heat_flux);
}

Vector3D FANS3D_ThermallyPerfect_KOmega_pState::qflux_t(
    const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
    const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
    const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) const{
  
    double Rmix = Rtot();
    
    Vector3D heat_flux, gradient_T;
    
    /* Temperature gradients from using the chain rule 
       with the ideal gas law (P=rho*R*T) 
       dT/dx = 1/(rho*R) *( dP/dx - P/rho * drho/dx) */
     
    gradient_T.x = (1.0/(rho*Rmix)) * (dWdx.p -(p/rho)*dWdx.rho);
    gradient_T.y = (1.0/(rho*Rmix)) * (dWdy.p -(p/rho)*dWdy.rho);
    gradient_T.z = (1.0/(rho*Rmix)) * (dWdz.p -(p/rho)*dWdz.rho);
    
    heat_flux = -kappa_t()*gradient_T;
    
    return (heat_flux);
}

// viscous flux
FANS3D_ThermallyPerfect_KOmega_cState  FANS3D_ThermallyPerfect_KOmega_pState::FluxViscous_n(
   const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
   const FANS3D_ThermallyPerfect_KOmega_pState &Wr,
   const FANS3D_ThermallyPerfect_KOmega_pState &Wc,
   const FANS3D_ThermallyPerfect_KOmega_pState &Wc_Neigbor,
   const FANS3D_ThermallyPerfect_KOmega_pState &dWdx,
   const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
   const FANS3D_ThermallyPerfect_KOmega_pState &dWdz,
   const FANS3D_ThermallyPerfect_KOmega_pState &dWdx_Neigbor,
   const FANS3D_ThermallyPerfect_KOmega_pState &dWdy_Neigbor,
   const FANS3D_ThermallyPerfect_KOmega_pState &dWdz_Neigbor,
   const Vector3D &norm, const Vector3D &ts, const double &deltad,
   const double &Volume, const double &Volume_Neigbor){
   
   // construct the gradients on the cell interface (surface) 
   // based on Hybrid Average Gradient-Diamond-Path Approach
   // Grad Phi = (Phi_c_neigbor - Phi_c)/ds *norm/(norm . ts)
   //            + [bar (Grad Phi) - bar (Grad Phi) . ts norm/(norm.ts)]

   // weighted factor based on volume
   double alpha = Volume/(Volume + Volume_Neigbor);
   
   FANS3D_ThermallyPerfect_KOmega_pState dWdx_Weighted, 
      dWdy_Weighted, dWdz_Weighted, dWdx_face, 
      dWdy_face, dWdz_face, Grad_middle_term;

   FANS3D_ThermallyPerfect_KOmega_pState W_face;
   
   dWdx_Weighted = alpha*dWdx + (1.0 - alpha)*dWdx_Neigbor;
   dWdy_Weighted = alpha*dWdy + (1.0 - alpha)*dWdy_Neigbor;
   dWdz_Weighted = alpha*dWdz + (1.0 - alpha)*dWdz_Neigbor;

   dWdx_face  =dWdx_Weighted ;
   dWdy_face  =dWdy_Weighted ;
   dWdz_face  =dWdz_Weighted ; 
   
   // a weighted term  
   //Grad_middle_term = dWdx_Weighted*ts.x + dWdy_Weighted*ts.y +
   //                   dWdz_Weighted*ts.z;
      
 // gradients of primitive variables on the face
 //  dWdx_face = (Wc_Neigbor - Wc)/deltad *norm.x/dot(norm, ts) + (
 //     dWdx_Weighted -  Grad_middle_term*norm.x/dot(norm, ts));
 //  dWdy_face = (Wc_Neigbor - Wc)/deltad *norm.y/dot(norm, ts) + (
 //     dWdy_Weighted -  Grad_middle_term*norm.y/dot(norm, ts));
 //  dWdz_face = (Wc_Neigbor - Wc)/deltad *norm.z/dot(norm, ts) + (
 //     dWdz_Weighted -  Grad_middle_term*norm.z/dot(norm, ts));
   
   // W_face = HALF*(Wl + Wr);
   W_face = HALF*(Wc + Wc_Neigbor);

   return (W_face.Fv(dWdx_face, dWdy_face, dWdz_face)*norm.x +
           W_face.Gv(dWdx_face, dWdy_face, dWdz_face)*norm.y +
           W_face.Hv(dWdx_face, dWdy_face, dWdz_face)*norm.z);
}

/* Turbulence model source term */  
FANS3D_ThermallyPerfect_KOmega_cState  FANS3D_ThermallyPerfect_KOmega_pState::Src_t(
   const FANS3D_ThermallyPerfect_KOmega_pState &Wc,
   const FANS3D_ThermallyPerfect_KOmega_pState &dWdx,
   const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
   const FANS3D_ThermallyPerfect_KOmega_pState &dWdz) {
   
   double mu_t, production;
   FANS3D_ThermallyPerfect_KOmega_cState Temp; Temp.Vacuum();
   Tensor3D Reynolds_stress;
   Tensor3D rotation_t, strain_rate_t;
   
   //Turbulence model eddy viscosity
   mu_t = Wc.eddy_viscosity();
   Reynolds_stress = Wc.lambda(dWdx, dWdy, dWdz);
  
   production = Reynolds_stress.xx*dWdx.v.x + 
      Reynolds_stress.xy*(dWdy.v.x + dWdx.v.y) +
      Reynolds_stress.yy*dWdy.v.y +
      Reynolds_stress.xz*(dWdz.v.x + dWdx.v.z) +
      Reynolds_stress.yz*(dWdz.v.y + dWdy.v.z) +
      Reynolds_stress.zz*dWdz.v.z;
   
   rotation_t = Wc.rotation_tensor(dWdx, dWdy, dWdz);
   strain_rate_t = Wc.strain_rate(dWdx, dWdy, dWdz);
   
   Temp.rhok = production - Wc.k_omega_model.beta_star*Wc.k_omega_model.f_betastar(
               dWdx.k, dWdy.k, dWdz.k, dWdx.omega, dWdy.omega, dWdz.omega, Wc.omega)*
               Wc.rho*Wc.k*Wc.omega;
   Temp.rhoomega = Wc.k_omega_model.alpha*(Wc.omega/max(Wc.k, TOLER))*production - 
                   Wc.k_omega_model.beta*Wc.k_omega_model.f_beta(rotation_t, strain_rate_t, Wc.omega)*
                   Wc.rho*Wc.omega*Wc.omega;
   //(1989)
   //Temp.rhok = production - Wc.k_omega_model.f_beta_star_const*Wc.k_omega_model.beta_star*Wc.rho*Wc.k*Wc.omega;
   //Temp.rhoomega = Wc.k_omega_model.alpha*(Wc.omega/max(Wc.k, TOLER))*production -
   //                Wc.k_omega_model.f_beta_const*Wc.k_omega_model.beta*Wc.rho*Wc.omega*Wc.omega;
   
   return (Temp);
} 



/*****************************************************************
 ** _pState::Sw -- Chemical Reaction Rate Source Terms.     **
 **                                                             **
 ** Using the Reaction class to get the source terms for the    ** 
 ** specific "Reaction_set".                                    ** 
 *****************************************************************
 *****************************************************************/
FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_pState::Sw(
   int &REACT_SET_FLAG, const int flow_type) const {
   FANS3D_ThermallyPerfect_KOmega_cState NEW;     
   NEW.Vacuum();
   bool test = negative_speccheck();
   //Adds concentration rate of change for species 1->N
   if( REACT_SET_FLAG != NO_REACTIONS){
      React.omega(NEW,*this, flow_type, k_omega_model.beta_star, omega);  
   }
   return NEW;
}

int FANS3D_ThermallyPerfect_KOmega_cState::Unphysical_Properties_Check(
   FANS3D_ThermallyPerfect_KOmega_cState &Uw,
   FANS3D_ThermallyPerfect_KOmega_cState &Ue,
   FANS3D_ThermallyPerfect_KOmega_cState &Ut,
   FANS3D_ThermallyPerfect_KOmega_cState &Ub, 
   FANS3D_ThermallyPerfect_KOmega_cState &Un,
   FANS3D_ThermallyPerfect_KOmega_cState &Us, 
   const int n) {
   
   int index;
   if (rhok <= ZERO ){
      //check all neigbours' rhok values and take the average of 
      //those positive value
      double average_k=ZERO; 
      index=0;
      if(Uw.rhok>ZERO){
         index +=1;
         average_k+= Uw.rhok;
      }
      if(Ue.rhok>ZERO){
         index +=1;
         average_k+=Ue.rhok;
      }
      if(Ut.rhok>ZERO){
         index +=1;
         average_k+=Ut.rhok;
      }
      if(Ub.rhok>ZERO){
         index +=1;
         average_k+=Ub.rhok;
      }
      if(Un.rhok>ZERO){
         index +=1;
         average_k+=Un.rhok;
      }
      if(Us.rhok>ZERO){
         index +=1;
         average_k+=Us.rhok;
      }
      if(index==0){ 
         rhok = rho*TOLER;
      }else{
         rhok = max(rho*TOLER, average_k/double(index));
         
      }//end check rhok
      
   }
      
   if (rhoomega<= rho*MILLI){
      double average_omega=ZERO; 
      index = 0;
      if(Uw.rhoomega>ZERO){
         index +=1;
         average_omega += Uw.rhoomega;
      }
      if(Ue.rhoomega>ZERO){
         index +=1;
         average_omega+=Ue.rhoomega;
      }
      if(Ut.rhoomega>ZERO){
         index +=1;
         average_omega+=Ut.rhoomega;
      }
      if(Ub.rhoomega>ZERO){
         index +=1;
         average_omega+=Ub.rhoomega;
      }
      if(Un.rhoomega>ZERO){
         index +=1;
         average_omega+=Un.rhoomega;
      }
      if(Us.rhoomega>ZERO){
         index +=1;
         average_omega+=Us.rhoomega;
      }
      average_omega = average_omega/double(index); 
      rhoomega = average_omega;
      
   } // end check rhoomega
   
//    if (rho <= ZERO || !negative_speccheck(n) || es() <= ZERO ||
//        rhok <= ZERO || rhoomega <= ZERO) {    
//       cout << "\n " << CFFC_Name() 
//            << " FANS3D ERROR: Negative Density || Energy || mass fractions || "
//            << " Turbulent kinetic energy || Dissipation rate: \n";
      
//       return false;
//    }else{
      
//       return true;
      
//    }
} 
 
