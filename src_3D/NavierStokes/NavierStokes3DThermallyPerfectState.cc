/********************** NavierStokes3DThermallyPerfectState.cc ******
  This file defines the various member functions of the 
  NavierStokes thermally perfect gaseous mixture class.

   assosicated files:
           NavierStokes3DThermallyPerfectState.h    

*********************************************************************/

// Include required CFFC header files.

#ifndef _NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED
#include "NavierStokes3DThermallyPerfectState.h"
#endif // NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED   

//set Global data for Species (STATIC, ie. so only call once! NavierStokes3DInput.h)
/* Get max of the min temperature of the lowest region
   and min of the max temperature of the highest region**/

void NavierStokes3D_ThermallyPerfect_pState::set_species_data(
                 NavierStokes3D_ThermallyPerfect_pState &Wo,
                 const int &n,
                 const string *S,
                 const char *PATH,
                 const int &debug, 
                 const double &Mr, 
                 const double* Sc){ 

   //Deallocate_static();
   Wo.Deallocate(); //Clean up memory before changing ns
   
   ns = n;
   NUM_VAR_3D = ns + NUM_EULER3D_VAR_SANS_SPECIES;
   
   //read in NASA data for each species to be used 
   specdata = new NASARP1311data[ns]; 
   Schmidt = new double[ns];
   // cout<<"\n Get by thse two "<<endl;
   
   for(int i=0; i<ns; i++){
     specdata[i].Getdata(S[i], PATH);  
      Schmidt[i] = Sc[i];
   } 

   //set data temperature ranges for mixture
   Wo.Temp_low_range(); 
   Wo.Temp_high_range(); 
   
   //Set Debug Information level
   debug_level = debug;
   
   // set initial values for the species
   Wo.spec = new Species[ns];
   for(int i=0; i<ns; i++){
      Wo.spec[i].c = ONE/ns ; 
   }
   
}   

void NavierStokes3D_ThermallyPerfect_cState::set_species_data(
                      NavierStokes3D_ThermallyPerfect_cState &Uo, 
                      const int &n, 
                      const string *S, 
                      const char *PATH,
                      const int &debug, 
                      const double &Mr, 
                      const double* Sc){ 

   // Deallocate_static();
   Uo.Deallocate(); //Clean up memory before changing ns
  
   ns =n; 
   
   NUM_VAR_3D = ns + NUM_EULER3D_VAR_SANS_SPECIES;

   //read in NASA data for each species to be used
   specdata = new NASARP1311data[ns];
   Schmidt = new double[ns];

   for(int i=0; i<ns; i++){
      //overwrite default data  
     specdata[i].Getdata(S[i], PATH);
      Schmidt[i] = Sc[i];  
   }  

   //set data temperature ranges for mixture
   Uo.Temp_low_range();
   Uo.Temp_high_range();

   //Set Debug Information level
   debug_level = debug;
   
   //setup initial array for mass fractions
   Uo.rhospec = new Species[ns];
   for(int i=0; i<ns; i++){
      Uo.rhospec[i].c = Uo.rho/ns; 
   }
   
}      

/********************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::F -- 
Inviscid flux (x-direction).   *
*********************************************************************/
NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_pState::F(void) {
   
   // call the inviscid flux calculation in NavierStokes3D_ThermallyPerfect_pState class  
  Euler3D_ThermallyPerfect_cState Temp;
   Temp = Euler3D_ThermallyPerfect_pState::F(); 

   return (NavierStokes3D_ThermallyPerfect_cState(Temp.rho, Temp.rhov, Temp.E, Temp.rhospec ));
 
}
NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_pState::F(void) const {
   
   // call the inviscid flux calculation in NavierStokes3D_ThermallyPerfect_pState class  
   Euler3D_ThermallyPerfect_cState Temp;
   Temp = Euler3D_ThermallyPerfect_pState::F();
   return (NavierStokes3D_ThermallyPerfect_cState(Temp.rho, Temp.rhov, Temp.E, Temp.rhospec ));
 
}

/********************************************************
 * NavierStokes3D_ThermallyPerfect_pState::Fv         -- 
Viscous flux (x-direction).   *
/********************************************************/
NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_pState::Fv(
   const NavierStokes3D_ThermallyPerfect_pState &dWdx,
   const NavierStokes3D_ThermallyPerfect_pState &dWdy,
   const NavierStokes3D_ThermallyPerfect_pState &dWdz){
   
   NavierStokes3D_ThermallyPerfect_cState Temp;
   Temp.Vacuum();
   
   double Temperature = T();
   
   Tensor3D molecular_stress;
   Vector3D heat_flux;
   
   molecular_stress = tau(dWdx, dWdy, dWdz);
   heat_flux = qflux(dWdx, dWdy, dWdz);
   
   /*----------------- Thermal Diffusion ---------------------------*/
   for (int k = 0; k<ns; ++k)
      spec[k].diffusion_coef = mu()/Schmidt[k];

   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux -=  rho*thermal_diffusion( );
   
   Temp.rhov.x = molecular_stress.xx;
   Temp.rhov.y = molecular_stress.xy;
   Temp.rhov.z = molecular_stress.xz;
   Temp.E = v.x*molecular_stress.xx + v.y*molecular_stress.xy 
      + v.z* molecular_stress.xz - heat_flux.x;
   
   // species transport 
   for (int k = 0; k<ns; ++k){
      Temp.rhospec[k] = rho*spec[k].diffusion_coef*dWdx.spec[k].c;
   }
   
   return (Temp);
   
}


NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_pState::Gv(
   const NavierStokes3D_ThermallyPerfect_pState &dWdx,
   const NavierStokes3D_ThermallyPerfect_pState &dWdy,
   const NavierStokes3D_ThermallyPerfect_pState &dWdz){
   

 NavierStokes3D_ThermallyPerfect_cState Temp;
   Temp.Vacuum();
   
   double Temperature = T();
   
   Tensor3D molecular_stress;
   Vector3D heat_flux;
   
   molecular_stress = tau(dWdx, dWdy, dWdz);
   heat_flux = qflux(dWdx, dWdy, dWdz);
   
   /*----------------- Thermal Diffusion ---------------------------*/
   for (int k = 0; k<ns; ++k)
      spec[k].diffusion_coef = mu()/Schmidt[k];

   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux -=  rho*thermal_diffusion();
   
   Temp.rhov.x = molecular_stress.xy;
   Temp.rhov.y = molecular_stress.yy;
   Temp.rhov.z = molecular_stress.yz;
   
   Temp.E = v.x*molecular_stress.xy+ v.y*molecular_stress.yy
      + v.z* molecular_stress.yz- heat_flux.y;
   
   
   // species transport 
   for (int k = 0; k<ns; ++k){
      Temp.rhospec[k] = rho*spec[k].diffusion_coef*dWdy.spec[k].c;
   }
   
   return (Temp);
   
   
}

NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_pState::Hv(
   const NavierStokes3D_ThermallyPerfect_pState &dWdx,
   const NavierStokes3D_ThermallyPerfect_pState &dWdy,
   const NavierStokes3D_ThermallyPerfect_pState &dWdz) {
   
   NavierStokes3D_ThermallyPerfect_cState Temp;
   Temp.Vacuum();
   
   double Temperature = T();
   
   Tensor3D molecular_stress;
   Vector3D heat_flux;
   
   molecular_stress = tau(dWdx, dWdy, dWdz);
   heat_flux = qflux(dWdx, dWdy, dWdz);
   
   /*----------------- Thermal Diffusion ---------------------------*/
   for (int k = 0; k<ns; ++k)
      spec[k].diffusion_coef = mu()/Schmidt[k];

   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux -=  rho*thermal_diffusion();
   
   Temp.rhov.x = molecular_stress.xz;
   Temp.rhov.y = molecular_stress.yz;
   Temp.rhov.z = molecular_stress.zz;
   Temp.E = v.x*molecular_stress.xz+ v.y*molecular_stress.yz
      + v.z* molecular_stress.zz - heat_flux.z;
   
   // species transport 
   for (int k = 0; k<ns; ++k){
      Temp.rhospec[k] = rho*spec[k].diffusion_coef*dWdz.spec[k].c;
   }
      return (Temp);
   
}


NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_pState::Fv(
   const NavierStokes3D_ThermallyPerfect_pState &dWdx,
   const NavierStokes3D_ThermallyPerfect_pState &dWdy,
   const NavierStokes3D_ThermallyPerfect_pState &dWdz) const{
   
   NavierStokes3D_ThermallyPerfect_cState Temp;
   Temp.Vacuum();
   
   double Temperature = T();
   
   Tensor3D molecular_stress;
   Vector3D heat_flux;
   
   molecular_stress = tau(dWdx, dWdy, dWdz);
   heat_flux = qflux(dWdx, dWdy, dWdz);
   
   /*----------------- Thermal Diffusion ---------------------------*/
   for (int k = 0; k<ns; ++k)
      spec[k].diffusion_coef = mu()/Schmidt[k];

   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux -=  rho*thermal_diffusion( );
   
   Temp.rhov.x = molecular_stress.xx;
   Temp.rhov.y = molecular_stress.xy;
   Temp.rhov.z = molecular_stress.xz;
   Temp.E = v.x*molecular_stress.xx + v.y*molecular_stress.xy 
      + v.z* molecular_stress.xz - heat_flux.x;
   
   // species transport 
   for (int k = 0; k<ns; ++k){
      Temp.rhospec[k] = rho*spec[k].diffusion_coef*dWdx.spec[k].c;
   }
   
   return (Temp);
   
}


NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_pState::Gv(
   const NavierStokes3D_ThermallyPerfect_pState &dWdx,
   const NavierStokes3D_ThermallyPerfect_pState &dWdy,
   const NavierStokes3D_ThermallyPerfect_pState &dWdz) const{
   

 NavierStokes3D_ThermallyPerfect_cState Temp;
   Temp.Vacuum();
   
   double Temperature = T();
   
   Tensor3D molecular_stress;
   Vector3D heat_flux;
   
   molecular_stress = tau(dWdx, dWdy, dWdz);
   heat_flux = qflux(dWdx, dWdy, dWdz);
   
   /*----------------- Thermal Diffusion ---------------------------*/
   for (int k = 0; k<ns; ++k)
      spec[k].diffusion_coef = mu()/Schmidt[k];

   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux -=  rho*thermal_diffusion();
   
   Temp.rhov.x = molecular_stress.xy;
   Temp.rhov.y = molecular_stress.yy;
   Temp.rhov.z = molecular_stress.yz;
   
   Temp.E = v.x*molecular_stress.xy+ v.y*molecular_stress.yy
      + v.z* molecular_stress.yz- heat_flux.y;
   
   
   // species transport 
   for (int k = 0; k<ns; ++k){
      Temp.rhospec[k] = rho*spec[k].diffusion_coef*dWdy.spec[k].c;
   }
   
   return (Temp);
   
   
}

NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_pState::Hv(
   const NavierStokes3D_ThermallyPerfect_pState &dWdx,
   const NavierStokes3D_ThermallyPerfect_pState &dWdy,
   const NavierStokes3D_ThermallyPerfect_pState &dWdz) const{
   
   NavierStokes3D_ThermallyPerfect_cState Temp;
   Temp.Vacuum();
   
   double Temperature = T();
   
   Tensor3D molecular_stress;
   Vector3D heat_flux;
   
   molecular_stress = tau(dWdx, dWdy, dWdz);
   heat_flux = qflux(dWdx, dWdy, dWdz);
   
   /*----------------- Thermal Diffusion ---------------------------*/
   for (int k = 0; k<ns; ++k)
      spec[k].diffusion_coef = mu()/Schmidt[k];

   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux -=  rho*thermal_diffusion();
   
   Temp.rhov.x = molecular_stress.xz;
   Temp.rhov.y = molecular_stress.yz;
   Temp.rhov.z = molecular_stress.zz;
   Temp.E = v.x*molecular_stress.xz+ v.y*molecular_stress.yz
      + v.z* molecular_stress.zz - heat_flux.z;
   
   // species transport 
   for (int k = 0; k<ns; ++k){
      Temp.rhospec[k] = rho*spec[k].diffusion_coef*dWdz.spec[k].c;
   }
   
   return (Temp);
   
}

/***************** EIGENVALUES *************************************
 *******************************************************************/

/************************************************************
 * NavierStokes3D_ThermallyPerfect_pState::lambda -- Eigenvalue(s) (x-direction).    *
 ************************************************************/
NavierStokes3D_ThermallyPerfect_pState NavierStokes3D_ThermallyPerfect_pState::lambda_x(void){
  
   double c = a();
   Euler3D_ThermallyPerfect_pState Temp;
   
   Temp = Euler3D_ThermallyPerfect_pState::lambda_x();
 
   
   return (NavierStokes3D_ThermallyPerfect_pState(Temp.rho, Temp.v, Temp.p, Temp.spec));
   
}

NavierStokes3D_ThermallyPerfect_pState NavierStokes3D_ThermallyPerfect_pState::lambda_x(void) const {
  
   double c = a();
   Euler3D_ThermallyPerfect_pState Temp;
   
   Temp = Euler3D_ThermallyPerfect_pState::lambda_x();
   
   return (NavierStokes3D_ThermallyPerfect_pState(Temp.rho, Temp.v, Temp.p, Temp.spec));
   
}

/********************************************************
 * NavierStokes3D_ThermallyPerfect_pState -- Binary arithmetic operators.        *
 ********************************************************/
NavierStokes3D_ThermallyPerfect_pState NavierStokes3D_ThermallyPerfect_pState::operator +(
   const NavierStokes3D_ThermallyPerfect_pState &W) const{ 
   
   NavierStokes3D_ThermallyPerfect_pState Temp(W.rho,W.v,W.p);
   Temp.Copy(*this);
   Temp += W;
   return Temp;
   
}

//------------------ Subtraction ------------------------//
NavierStokes3D_ThermallyPerfect_pState NavierStokes3D_ThermallyPerfect_pState::operator -(
   const NavierStokes3D_ThermallyPerfect_pState &W) const{
   
   NavierStokes3D_ThermallyPerfect_pState Temp(W.rho,W.v,W.p);
   Temp.Copy(*this);
   Temp -= W;
   
   return Temp;
   
}

//---------------- Scalar Multiplication ------------------//
NavierStokes3D_ThermallyPerfect_pState NavierStokes3D_ThermallyPerfect_pState::operator *(
   const double &a) const{
   NavierStokes3D_ThermallyPerfect_pState Temp(rho,v,p);
   Temp.Copy(*this);
   Temp.rho = rho*a;  Temp.v = v*a; Temp.p = p*a;
   for( int i=0; i<ns; i++){
      Temp.spec[i] = spec[i]*a;
   } 
   
   return(Temp);
}

NavierStokes3D_ThermallyPerfect_pState operator *(
   const double &a, 
   const NavierStokes3D_ThermallyPerfect_pState &W){
   NavierStokes3D_ThermallyPerfect_pState Temp;
   //Temp.Copy(W);
   Temp.rho = W.rho*a;  Temp.v = W.v*a; Temp.p = W.p*a;
   
   for( int i=0; i<W.ns; i++){
      Temp.spec[i] = W.spec[i]*a;
   } 

   return(Temp);
}

//--------------- Scalar Division ------------------------//
NavierStokes3D_ThermallyPerfect_pState NavierStokes3D_ThermallyPerfect_pState::operator /(
   const double &a) const {
   NavierStokes3D_ThermallyPerfect_pState Temp(rho,v,p);
   Temp.Copy(*this);
   Temp.rho = rho/a; Temp.v = v/a; Temp.p = p/a; 
   
   for(int i=0; i<ns; i++){
      Temp.spec[i] = spec[i]/a; 
   } 

   return(Temp);
}

//----------------- Inner Product ------------------------//
double NavierStokes3D_ThermallyPerfect_pState::operator *(
   const NavierStokes3D_ThermallyPerfect_pState &W) const{
   double sum=0.0;
   
   for(int i=0; i<ns; i++){
      sum += spec[i]*W.spec[i];
   }  
   return (rho*W.rho + v*W.v + p*W.p + sum);
   
}

//----------- solution state product operator ------------//
NavierStokes3D_ThermallyPerfect_pState NavierStokes3D_ThermallyPerfect_pState::operator ^(
   const NavierStokes3D_ThermallyPerfect_pState &W) const{
   
   NavierStokes3D_ThermallyPerfect_pState Temp(rho,v,p);
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
   
}

//----------------- Assignment ----------------------------//
NavierStokes3D_ThermallyPerfect_pState& NavierStokes3D_ThermallyPerfect_pState::operator =(
   const NavierStokes3D_ThermallyPerfect_pState &W){
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
/********************************************************
 * NavierStokes3D_ThermallyPerfect_pState -- Shortcut arithmetic operators.     *
 ********************************************************/
NavierStokes3D_ThermallyPerfect_pState& NavierStokes3D_ThermallyPerfect_pState::operator +=(
   const NavierStokes3D_ThermallyPerfect_pState &W){
   rho += W.rho;
   v += W.v; 
   p += W.p; 
   
   for( int i=0; i<ns; i++){
      spec[i] += W.spec[i];
   } 
   
   return (*this);
}

NavierStokes3D_ThermallyPerfect_pState& NavierStokes3D_ThermallyPerfect_pState::operator -=(
   const NavierStokes3D_ThermallyPerfect_pState &W) {
   rho -= W.rho;
   v -= W.v;
   p -= W.p;
   
   for(int i=0; i<ns; i++){
      spec[i] -= W.spec[i];
   }  
   
   return (*this); 
}

/********************************************************
 * NavierStokes3D_ThermallyPerfect_pState -- Relational operators.               *
 ********************************************************/
int operator ==(const NavierStokes3D_ThermallyPerfect_pState &W1, 
                const NavierStokes3D_ThermallyPerfect_pState &W2) {
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

int operator !=(const NavierStokes3D_ThermallyPerfect_pState &W1, 
                const NavierStokes3D_ThermallyPerfect_pState &W2) {
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

/* NavierStokes3D_ThermallyPerfect_pState -- Input-output operators.*/
ostream &operator << (ostream &out_file, 
                      const NavierStokes3D_ThermallyPerfect_pState &W) {
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

istream &operator >> (istream &in_file, NavierStokes3D_ThermallyPerfect_pState &W) {
   in_file.setf(ios::skipws);
   in_file >> W.rho >> W.v.x >> W.v.y >> W.v.z >>W.p;
   //W.set_initial_values();
   for( int i=0; i<W.ns; i++){
      in_file>>W.spec[i];
   }
   in_file.unsetf(ios::skipws);
   return (in_file);
}


//----------------- Addition -----------------------------//
NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_cState::operator +(
   const NavierStokes3D_ThermallyPerfect_cState &U) const{ 
   
   NavierStokes3D_ThermallyPerfect_cState Temp(U.rho,U.rhov,U.E);
   Temp.Copy(*this);
   Temp += U;
   return Temp;
   
}

//------------------ Subtraction ------------------------//
NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_cState::operator -(
   const NavierStokes3D_ThermallyPerfect_cState &U) const{
   NavierStokes3D_ThermallyPerfect_cState Temp(U.rho,U.rhov,U.E);
   Temp.Copy(*this);
   Temp -= U;
   return Temp;

}

//---------------- Scalar Multiplication ------------------//
NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_cState::operator *(
   const double &a) const{
   NavierStokes3D_ThermallyPerfect_cState Temp(rho,rhov,E);
   Temp.Copy(*this);
   Temp.rho = rho*a;  Temp.rhov = rhov*a; Temp.E = E*a;
   for( int i=0; i<ns; i++){
      Temp.rhospec[i] = rhospec[i]*a;
   } 
   
   return(Temp);
}

NavierStokes3D_ThermallyPerfect_cState operator *(
   const double &a, 
   const NavierStokes3D_ThermallyPerfect_cState &U){
   NavierStokes3D_ThermallyPerfect_cState Temp;
   Temp.rho = U.rho*a;  Temp.rhov = U.rhov*a; Temp.E = U.E*a;
   for( int i=0; i<U.ns; i++){
      Temp.rhospec[i] = U.rhospec[i]*a;
   } 
   return(Temp);
}
//--------------- Scalar Division ------------------------//
NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_cState::operator /(
   const double &a) const {
   NavierStokes3D_ThermallyPerfect_cState Temp(rho,rhov,E);
   Temp.Copy(*this);
   Temp.rho = rho/a; Temp.rhov = rhov/a; Temp.E = E/a;
   
   for(int i=0; i<ns; i++){
      Temp.rhospec[i] = rhospec[i]/a; 
   } 
   
   return(Temp);
}

//----------------- Inner Product ------------------------//
double NavierStokes3D_ThermallyPerfect_cState::operator *(
   const NavierStokes3D_ThermallyPerfect_cState &U) const{
   double sum=0.0;
   
   for(int i=0; i<ns; i++){
      sum += rhospec[i]*U.rhospec[i];
   }  
   return (rho*U.rho + rhov*U.rhov + E*U.E + sum);
   
}

//----------- solution state product operator ------------//
NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_cState::operator ^(
   const NavierStokes3D_ThermallyPerfect_cState &U) const {
   
   NavierStokes3D_ThermallyPerfect_cState Temp(rho,rhov,E);
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
   
}

//----------------- Assignment ----------------------------//
NavierStokes3D_ThermallyPerfect_cState& NavierStokes3D_ThermallyPerfect_cState::operator =(
   const NavierStokes3D_ThermallyPerfect_cState &U){
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
 * NavierStokes3D_ThermallyPerfect_cState -- Shortcut arithmetic operators.     *
 ********************************************************/
NavierStokes3D_ThermallyPerfect_cState& NavierStokes3D_ThermallyPerfect_cState::operator +=(
   const NavierStokes3D_ThermallyPerfect_cState &U){
   rho += U.rho;
   rhov += U.rhov; 
   E += U.E;
   
   for( int i=0; i<ns; i++){
      rhospec[i] += U.rhospec[i];
   } 

   return (*this);
}

NavierStokes3D_ThermallyPerfect_cState& NavierStokes3D_ThermallyPerfect_cState::operator -=(
   const NavierStokes3D_ThermallyPerfect_cState &U) {
  
   rho -= U.rho;
   rhov -= U.rhov;
   E -= U.E;
   
  for(int i=0; i<ns; i++){
    rhospec[i] -= U.rhospec[i];
  }  

  return (*this); 
}

/********************************************************
 * NavierStokes3D_ThermallyPerfect_cState -- Unary arithmetic operators.        *
 ********************************************************/
NavierStokes3D_ThermallyPerfect_cState operator -(
   const NavierStokes3D_ThermallyPerfect_cState &U) {
   Species *spt= new Species[U.ns];
   for(int i=0; i<U.ns; i++){
      spt[i] = -U.rhospec[i]; 
   }  
   
   NavierStokes3D_ThermallyPerfect_cState Temp(-U.rho,-U.rhov,-U.E, spt);
  
  delete[] spt;
  return(Temp);
}

/********************************************************
 * NavierStokes3D_ThermallyPerfect_cState -- Relational operators.              *
 ********************************************************/
int operator ==(const NavierStokes3D_ThermallyPerfect_cState &U1, 
                const NavierStokes3D_ThermallyPerfect_cState &U2) {
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

int operator !=(const NavierStokes3D_ThermallyPerfect_cState &U1, 
                const NavierStokes3D_ThermallyPerfect_cState &U2) {
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
 * NavierStokes3D_ThermallyPerfect_cState -- Input-output operators.            *
 ********************************************************/
ostream &operator << (ostream &out_file, const NavierStokes3D_ThermallyPerfect_cState &U) {
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

istream &operator >> (istream &in_file, NavierStokes3D_ThermallyPerfect_cState &U) {
   in_file.setf(ios::skipws);
   in_file >> U.rho >> U.rhov.x >> U.rhov.y >> U.rhov.z >> U.E;
  //U.set_initial_values();
  for( int i=0; i<U.ns; i++){
    in_file>>U.rhospec[i]; 
  } 

  in_file.unsetf(ios::skipws);
  return (in_file);
}

NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_pState::rc_x(
   const int &index) const {
   switch(index){  
   case 1: 
      return (NavierStokes3D_ThermallyPerfect_cState(
                 ONE, v.x-a(), v.y, v.z, H()/rho-v.x*a(), spec));
   case 2:
      return (NavierStokes3D_ThermallyPerfect_cState(
                 ONE, v.x, v.y, v.z, H()/rho-Cp()*T(), spec)); 
   case 3:
      return (NavierStokes3D_ThermallyPerfect_cState(
                 ZERO, ZERO, rho, ZERO, rho*v.y, ZERO));
   case 4: 
      return (NavierStokes3D_ThermallyPerfect_cState(
                 ZERO, ZERO, ZERO, rho, rho*v.z, ZERO));
   case 5: 
      return (NavierStokes3D_ThermallyPerfect_cState(
                 ONE, v.x+a(), v.y, v.z, H()/rho+v.x*a(), spec));
      
   default: 
      NavierStokes3D_ThermallyPerfect_cState NEW(ZERO);
      double RTOT = Rtot();
      double TEMP = p/(rho*RTOT);      
      NEW.E = rho*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP) + specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Heatofform() 
                   - Cp(TEMP)*TEMP*specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs()/RTOT);
      double PHI =specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP)-specdata[NUM_VAR_3D- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP)-
         (Cp(TEMP) -RTOT)*TEMP*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs() - specdata[NUM_VAR_3D - (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs())/RTOT;
            
      NEW.rhospec[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].c = rho*PHI; 
            
            return NEW;

   };
   
 
}

NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_pState::rc_x(
   const int &index) {
   switch(index){  
   case 1: 
      return (NavierStokes3D_ThermallyPerfect_cState(
                 ONE, v.x-a(), v.y, v.z, H()/rho-v.x*a(), spec));
   case 2:
      return (NavierStokes3D_ThermallyPerfect_cState(
                 ONE, v.x, v.y, v.z, H()/rho-Cp()*T(), spec)); 
   case 3:
      return (NavierStokes3D_ThermallyPerfect_cState(
                 ZERO, ZERO, rho, ZERO, rho*v.y, ZERO));
   case 4: 
      return (NavierStokes3D_ThermallyPerfect_cState(
                 ZERO, ZERO, ZERO, rho, rho*v.z, ZERO));
   case 5: 
      
      return (NavierStokes3D_ThermallyPerfect_cState(
                 ONE, v.x+a(), v.y, v.z, H()/rho+v.x*a(), spec));
      
   default: 
      NavierStokes3D_ThermallyPerfect_cState NEW(ZERO);
      double RTOT = Rtot();
      double TEMP = p/(rho*RTOT);      
      NEW.E = rho*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP) + specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Heatofform() 
                   - Cp(TEMP)*TEMP*specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs()/RTOT);
      double PHI =specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP)-specdata[NUM_VAR_3D- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP)-
         (Cp(TEMP) -RTOT)*TEMP*(specdata[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs() - specdata[NUM_VAR_3D - (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs())/RTOT;
      
      NEW.rhospec[index- (NUM_EULER3D_VAR_SANS_SPECIES+1)].c = rho*PHI; 
      
      return NEW;
      
   };
   
 
}

// Primitive Left Eigenvector -- (x-direction)
NavierStokes3D_ThermallyPerfect_pState NavierStokes3D_ThermallyPerfect_pState::lp_x(
   const int &index) const {
   switch(index){  
   case 1 :
      return (NavierStokes3D_ThermallyPerfect_pState(
                 ZERO, -HALF*rho/a(), ZERO, ZERO, HALF/(a()*a()),  ZERO));
   case 2:
      
      return (NavierStokes3D_ThermallyPerfect_pState(
                 ONE, ZERO, ZERO, ZERO, -ONE/(a()*a()),  ZERO));
   case 3:
      return  (NavierStokes3D_ThermallyPerfect_pState(ZERO, ZERO, ONE, ZERO, ZERO, ZERO));
   case 4 :
      return  (NavierStokes3D_ThermallyPerfect_pState(ZERO, ZERO, ZERO, ONE, ZERO, ZERO));
   case 5 :
      
      return (NavierStokes3D_ThermallyPerfect_pState(
                 ZERO, HALF*rho/a(), ZERO, ZERO, HALF/(a()*a()), ZERO));
   default:
      NavierStokes3D_ThermallyPerfect_pState NEW(ZERO);
      NEW.spec[index-(NUM_EULER3D_VAR_SANS_SPECIES+1)].c = ONE;
      return NEW;
      
   };
    
    
}

NavierStokes3D_ThermallyPerfect_pState NavierStokes3D_ThermallyPerfect_pState::lp_x(
   const int &index) {
   switch(index){
   case 1 :
      return (NavierStokes3D_ThermallyPerfect_pState(
                 ZERO, -HALF*rho/a(), ZERO, ZERO, HALF/(a()*a()),  ZERO));
   case 2:

      return (NavierStokes3D_ThermallyPerfect_pState(
                 ONE, ZERO, ZERO, ZERO, -ONE/(a()*a()),  ZERO));
   case 3:
      return  (NavierStokes3D_ThermallyPerfect_pState(ZERO, ZERO, ONE, ZERO, ZERO, ZERO));
   case 4 :
      return  (NavierStokes3D_ThermallyPerfect_pState(ZERO, ZERO, ZERO, ONE, ZERO, ZERO));
   case 5 :

      return (NavierStokes3D_ThermallyPerfect_pState(
                 ZERO, HALF*rho/a(), ZERO, ZERO, HALF/(a()*a()), ZERO));
   default:
      NavierStokes3D_ThermallyPerfect_pState NEW(ZERO);
      NEW.spec[index-(NUM_EULER3D_VAR_SANS_SPECIES+1)].c = ONE;
      return NEW;
      
   };
    
    
}

// NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_pState::rc_x(
//    const int &index) {
//    assert( index >= 1 && index <= NUM_VAR_3D);
//    if(index == 1){
//       double c = a(); 
//       return (NavierStokes3D_ThermallyPerfect_cState(
//                  ONE, v.x-c, v.y, v.z, H()/rho-v.x*c, spec));
//    } else if(index == 2) {
//       return (NavierStokes3D_ThermallyPerfect_cState(
//                  ONE, v.x, v.y, v.z, H()/rho-Cp()*T(), spec)); 
//    } else if(index == 3) {
//       return (NavierStokes3D_ThermallyPerfect_cState(
//                  ZERO, ZERO, rho, ZERO, rho*v.y, ZERO));
//    } else if(index == 4) {
//       return (NavierStokes3D_ThermallyPerfect_cState(
//                  ZERO, ZERO, ZERO, rho, rho*v.z, ZERO));
//    }else if(index == 5) {
//       double c = a(); 
//       return (NavierStokes3D_ThermallyPerfect_cState(
//                  ONE, v.x+c, v.y, v.z, H()/rho+v.x*c, spec));
      
//    } else{ 
//       for(int i = NUM_EULER3D_VAR_SANS_SPECIES+1; i<=NUM_VAR_3D; i++){
//          if(index == i){
//             NavierStokes3D_ThermallyPerfect_cState NEW(ZERO);

//             double RTOT = Rtot();
//             double TEMP = p/(rho*RTOT);      
//             NEW.E = rho*(specdata[i- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP) + specdata[i- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Heatofform() 
//                          - Cp(TEMP)*TEMP*specdata[i- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs()/RTOT);
//             double PHI =specdata[i- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP)-specdata[NUM_VAR_3D- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Enthalpy(TEMP)-
//                (Cp(TEMP) -RTOT)*TEMP*(specdata[i- (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs() - specdata[NUM_VAR_3D - (NUM_EULER3D_VAR_SANS_SPECIES+1)].Rs())/RTOT;
            
//             NEW.rhospec[i- (NUM_EULER3D_VAR_SANS_SPECIES+1)].c = rho*PHI; 
            
//             return NEW;
//             break;
//          }
//       }
//    }
// }

// // Primitive Left Eigenvector -- (x-direction)
// NavierStokes3D_ThermallyPerfect_pState NavierStokes3D_ThermallyPerfect_pState::lp_x(
//    const int &index) const {
//    assert( index >= 1 && index <= NUM_VAR_3D );
//    if(index == 1){
//       double c = a(); 
//       return (NavierStokes3D_ThermallyPerfect_pState(
//                  ZERO, -HALF*rho/c, ZERO, ZERO, HALF/(c*c),  ZERO));
//    } else if(index == 2) {
//       double c = a(); 
//       return (NavierStokes3D_ThermallyPerfect_pState(
//                  ONE, ZERO, ZERO, ZERO, -ONE/(c*c),  ZERO));
//    } else if(index == 3) {
//       return  (NavierStokes3D_ThermallyPerfect_pState(ZERO, ZERO, ONE, ZERO, ZERO, ZERO));
//    } else if(index == 4) {
//       return  (NavierStokes3D_ThermallyPerfect_pState(ZERO, ZERO, ZERO, ONE, ZERO, ZERO));
//    } else if(index == 5) {  
//       double c = a(); 
//       return (NavierStokes3D_ThermallyPerfect_pState(
//                  ZERO, HALF*rho/c, ZERO, ZERO, HALF/(c*c), ZERO));
//    } else{ 
//       for(int i=NUM_EULER3D_VAR_SANS_SPECIES+1; i<=NUM_VAR_3D; i++){
// 	if(index == i){
// 	  NavierStokes3D_ThermallyPerfect_pState NEW(ZERO);
// 	  NEW.spec[i-(NUM_EULER3D_VAR_SANS_SPECIES+1)].c = ONE;
// 	  return NEW;
// 	  break;
// 	}
//       }
//     } 

// }



/********************************************************
 * Routine: RoeAverage (Roe Averages)                   *
 *                                                      *
 * This function returns the Roe-averaged primitive     *
 * solution state given left and right primitive        *
 * solution variables.  See Roe (1981).                 *
 *                                                      *
 ********************************************************/
NavierStokes3D_ThermallyPerfect_pState NavierStokes3D_ThermallyPerfect_pState::RoeAverage(
   const NavierStokes3D_ThermallyPerfect_pState &Wl,
   const NavierStokes3D_ThermallyPerfect_pState &Wr) {


   double Hl, Hr, srhol, srhor;
   double Ha, ha;
   NavierStokes3D_ThermallyPerfect_pState Temp;

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
NavierStokes3D_ThermallyPerfect_cState NavierStokes3D_ThermallyPerfect_pState::FluxHLLE_x(
   const NavierStokes3D_ThermallyPerfect_pState &Wl,
   const NavierStokes3D_ThermallyPerfect_pState &Wr) {
   
   double wavespeed_l, wavespeed_r;

 
    /* solnvec in  Wa (lambdas_l, lambdas_r, lambdas_a) 
       is allocated using new  */ 
   NavierStokes3D_ThermallyPerfect_pState Wa, lambdas_l, lambdas_r, lambdas_a;
   NavierStokes3D_ThermallyPerfect_cState Flux, dUrl;
   
   int NUM_VAR_3D = Wl.NUM_VAR_3D;
      
 
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
   wavespeed_r = max(lambdas_r[NUM_VAR_3D-lambdas_r.ns],
                     lambdas_a[NUM_VAR_3D-lambdas_a.ns]);

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

NavierStokes3D_ThermallyPerfect_cState  NavierStokes3D_ThermallyPerfect_pState::FluxHLLE_x(
   const NavierStokes3D_ThermallyPerfect_cState &Ul,
   const NavierStokes3D_ThermallyPerfect_cState &Ur) {
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
NavierStokes3D_ThermallyPerfect_cState  NavierStokes3D_ThermallyPerfect_pState::FluxHLLE_n(
   const NavierStokes3D_ThermallyPerfect_pState &Wl,
   const NavierStokes3D_ThermallyPerfect_pState &Wr,
   const Vector3D &norm_dir) {

   double Wl_ur_norm, Wl_ur_tang;
   double Wr_ur_norm, Wr_ur_tang ;
   double Wr_ur_tang_z;
   
   Vector3D Flux_rotated_x, Flux_rotated_tang_y, Flux_rotated_tang_z ;
   Vector3D Wl_ur_tang_vector, Wr_ur_tang_vector;
   Vector3D Wl_ur_tang_unit_vector, Wr_ur_tang_unit_vector;
   Vector3D Wr_ur_tang_z_vector, Wr_ur_tang_z_unit_vector;
   
   //solnvec in  Wl_rotated (Wr_rotated) is allocated using new 
   NavierStokes3D_ThermallyPerfect_pState Wl_rotated, Wr_rotated;
   NavierStokes3D_ThermallyPerfect_cState Flux, Flux_rotated;
   
   
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


NavierStokes3D_ThermallyPerfect_cState  NavierStokes3D_ThermallyPerfect_pState::FluxHLLE_n(const NavierStokes3D_ThermallyPerfect_cState &Ul,
	      	          const NavierStokes3D_ThermallyPerfect_cState &Ur,
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
NavierStokes3D_ThermallyPerfect_pState  HartenFixPos(
   const NavierStokes3D_ThermallyPerfect_pState &lambdas_a,
   const NavierStokes3D_ThermallyPerfect_pState &lambdas_l,
   const NavierStokes3D_ThermallyPerfect_pState &lambdas_r) {
   
   NavierStokes3D_ThermallyPerfect_pState NEW;
   NEW.rho = HartenFixPos(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
   NEW.v.x = HALF*(lambdas_a[2]+fabs(lambdas_a[2]));
   NEW.v.y = HALF*(lambdas_a[3]+fabs(lambdas_a[3]));
   NEW.v.z = HALF*(lambdas_a[4]+fabs(lambdas_a[4]));

   NEW.p = HartenFixPos(lambdas_a[5],lambdas_l[5],lambdas_r[5]);
  
   for( int i=(NUM_EULER3D_VAR_SANS_SPECIES +1); i<=NEW.NUM_VAR_3D; i++){
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
NavierStokes3D_ThermallyPerfect_pState  HartenFixNeg(
   const NavierStokes3D_ThermallyPerfect_pState &lambdas_a,
   const NavierStokes3D_ThermallyPerfect_pState &lambdas_l,
   const NavierStokes3D_ThermallyPerfect_pState &lambdas_r) {
  
   NavierStokes3D_ThermallyPerfect_pState NEW;
   
   NEW.rho = HartenFixNeg(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
   NEW.v.x = HALF*(lambdas_a[2]-fabs(lambdas_a[2]));
   NEW.v.y = HALF*(lambdas_a[3]-fabs(lambdas_a[3]));
   NEW.v.z = HALF*(lambdas_a[4]-fabs(lambdas_a[4]));
   NEW.p = HartenFixNeg(lambdas_a[5],lambdas_l[5],lambdas_r[5]);
   
   for( int i = (NUM_EULER3D_VAR_SANS_SPECIES +1); i<=NEW.NUM_VAR_3D; i++){
      NEW.spec[i-(NUM_EULER3D_VAR_SANS_SPECIES+1)].c = 
         HALF*(lambdas_a[i]-fabs(lambdas_a[i]));
   }
   
   return (NEW);
}


// Flux Roe -- based on Harten fix 
NavierStokes3D_ThermallyPerfect_cState  NavierStokes3D_ThermallyPerfect_pState::FluxRoe_x(
   const  NavierStokes3D_ThermallyPerfect_pState &Wl,  
   const  NavierStokes3D_ThermallyPerfect_pState &Wr){
   
   NavierStokes3D_ThermallyPerfect_pState Wa, dWrl, wavespeeds, 
      lambdas_l, lambdas_r, lambdas_a;
   NavierStokes3D_ThermallyPerfect_cState Flux;

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
      
          
      for (int i=1 ; i <= Wl.NUM_VAR_3D; i++) {
         if (wavespeeds[i] < ZERO) {
            Flux += wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
            
            
         }
      }

      
   } else {
      
   
  
      Flux = Wr.F();
      
      wavespeeds = HartenFixPos(lambdas_a,
                                lambdas_l,
                                lambdas_r);
      
      
      for (int i=1; i <= Wl.NUM_VAR_3D; i++) {
         if (wavespeeds[i] > ZERO) {
            
            Flux -= wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i); 
            
         }
      }

      
   } 
    
   /* Return solution flux. */    
   return (Flux);    
    
   
}
   
NavierStokes3D_ThermallyPerfect_cState  NavierStokes3D_ThermallyPerfect_pState::FluxRoe_n(
   const NavierStokes3D_ThermallyPerfect_pState &Wl,
   const NavierStokes3D_ThermallyPerfect_pState &Wr,
   const Vector3D &norm_dir){
   
   double Wl_ur_norm, Wl_ur_tang;
   double Wr_ur_norm, Wr_ur_tang ;
   double Wr_ur_tang_z;
   
   Vector3D Flux_rotated_x, Flux_rotated_tang_y, Flux_rotated_tang_z ;
   Vector3D Wl_ur_tang_vector, Wr_ur_tang_vector;
   Vector3D Wl_ur_tang_unit_vector, Wr_ur_tang_unit_vector;
   Vector3D Wr_ur_tang_z_vector, Wr_ur_tang_z_unit_vector;
   
   //solnvec in  Wl_rotated (Wr_rotated) is allocated using new 
   NavierStokes3D_ThermallyPerfect_pState Wl_rotated, Wr_rotated;
   NavierStokes3D_ThermallyPerfect_cState Flux, Flux_rotated;
   
   
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
NavierStokes3D_ThermallyPerfect_pState  NavierStokes3D_ThermallyPerfect_pState::Reflect(
   const NavierStokes3D_ThermallyPerfect_pState &W,
   const Vector3D &norm_dir) {
   

   Vector3D ur_norm, ur_tang, vr_tot;
   
   NavierStokes3D_ThermallyPerfect_pState Temp;
   Temp = W;
   
   ur_norm = dot(W.v, norm_dir)*norm_dir;
   ur_tang = W.v - ur_norm;
   
   ur_norm = -ur_norm;
   vr_tot = ur_norm + ur_tang;
   
   Temp.v = vr_tot;
   
   return (Temp);
       
}


NavierStokes3D_ThermallyPerfect_pState  NavierStokes3D_ThermallyPerfect_pState::Moving_Wall(
   const NavierStokes3D_ThermallyPerfect_pState &Win,
   const NavierStokes3D_ThermallyPerfect_pState &Wout,
   const Vector3D &norm_dir, 
   const Vector3D &wall_velocity,
   const Vector3D &pressure_gradient,
   const int &TEMPERATURE_BC_FLAG) {
   
   NavierStokes3D_ThermallyPerfect_pState Temp;
   Temp = Win;
   
   if(wall_velocity ==  Vector3D_ZERO){
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
   
   
   // Fixed Wall Temperature or constant extrapolation for Adiabatic
   if(TEMPERATURE_BC_FLAG == FIXED_TEMPERATURE_WALL){
      if (pressure_gradient != Vector3D_ZERO){
         Temp.rho = Wout.p/(Temp.Rtot()*Wout.T());
      }else{
         Temp.rho = Win.p/(Temp.Rtot()*Wout.T());
      }
     
   }
    
   return (Temp);
   
}

NavierStokes3D_ThermallyPerfect_pState  NavierStokes3D_ThermallyPerfect_pState::No_Slip(
   const NavierStokes3D_ThermallyPerfect_pState &Win,
   const NavierStokes3D_ThermallyPerfect_pState &Wout,
   const Vector3D &norm_dir, 
   const Vector3D &pressure_gradient,
   const int &TEMPERATURE_BC_FLAG) {
   
   return (Moving_Wall(Win, Wout, norm_dir, Vector3D_ZERO,pressure_gradient,TEMPERATURE_BC_FLAG));
   
}


// molecular stress tensor
Tensor3D NavierStokes3D_ThermallyPerfect_pState::tau(
   const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
   const NavierStokes3D_ThermallyPerfect_pState &dWdy,
   const NavierStokes3D_ThermallyPerfect_pState &dWdz){

   Tensor3D molecular_stress;
   
   molecular_stress.xx = 1.0/3.0 *mu()*(4.0*dWdx.v.x - 2.0*dWdy.v.y - 2.0*dWdz.v.z);
   molecular_stress.xy = mu()*(dWdx.v.y + dWdy.v.x);
   molecular_stress.xz = mu()*(dWdx.v.z + dWdz.v.x);
   molecular_stress.yy = 1.0/3.0 *mu()*(4.0*dWdy.v.y - 2.0*dWdx.v.x - 2.0*dWdz.v.z);
   molecular_stress.yz = mu()*(dWdy.v.z + dWdz.v.y);
   molecular_stress.zz = 1.0/3.0 *mu()*(4.0*dWdz.v.z - 2.0*dWdx.v.x - 2.0*dWdy.v.y);
    
   return (molecular_stress);
   
}

Tensor3D NavierStokes3D_ThermallyPerfect_pState::tau(
   const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
   const NavierStokes3D_ThermallyPerfect_pState &dWdy,
   const NavierStokes3D_ThermallyPerfect_pState &dWdz) const {

   Tensor3D molecular_stress;
   
   molecular_stress.xx = 1.0/3.0 *mu()*(4.0*dWdx.v.x - 2.0*dWdy.v.y - 2.0*dWdz.v.z);
   molecular_stress.xy = mu()*(dWdx.v.y + dWdy.v.x);
   molecular_stress.xz = mu()*(dWdx.v.z + dWdz.v.x);
   molecular_stress.yy = 1.0/3.0 *mu()*(4.0*dWdy.v.y - 2.0*dWdx.v.x - 2.0*dWdz.v.z);
   molecular_stress.yz = mu()*(dWdy.v.z + dWdz.v.y);
   molecular_stress.zz = 1.0/3.0 *mu()*(4.0*dWdz.v.z - 2.0*dWdx.v.x - 2.0*dWdy.v.y);

   return (molecular_stress);
   
}

// heat flux: Fourier's law of heat conduction
 Vector3D NavierStokes3D_ThermallyPerfect_pState::qflux(
    const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
    const NavierStokes3D_ThermallyPerfect_pState &dWdy,
    const NavierStokes3D_ThermallyPerfect_pState &dWdz){
  
    double Rmix = Rtot();
    
    Vector3D heat_flux, gradient_T;
    
  /* Temperature gradients from using the chain rule 
     with the ideal gas law (P=rho*R*T) 
     dT/dx = 1/(rho*R) *( dP/dx - P/rho * drho/dx) */
     
    gradient_T.x = (1.0/(rho*Rmix)) * (dWdx.p -(p/rho)*dWdx.rho);
    gradient_T.y = (1.0/(rho*Rmix)) * (dWdy.p -(p/rho)*dWdy.rho);
    gradient_T.z = (1.0/(rho*Rmix)) * (dWdz.p -(p/rho)*dWdz.rho);
    
    heat_flux = -kappa()*gradient_T;
        
    return (heat_flux);
    
   
 }
 Vector3D NavierStokes3D_ThermallyPerfect_pState::qflux(
    const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
    const NavierStokes3D_ThermallyPerfect_pState &dWdy,
    const NavierStokes3D_ThermallyPerfect_pState &dWdz) const{
  
    double Rmix = Rtot();
    
    Vector3D heat_flux, gradient_T;
    
  /* Temperature gradients from using the chain rule 
     with the ideal gas law (P=rho*R*T) 
     dT/dx = 1/(rho*R) *( dP/dx - P/rho * drho/dx) */
     
    gradient_T.x = (1/(rho*Rmix)) * (dWdx.p -(p/rho)*dWdx.rho);
    gradient_T.y = (1/(rho*Rmix)) * (dWdy.p -(p/rho)*dWdy.rho);
    gradient_T.z = (1/(rho*Rmix)) * (dWdz.p -(p/rho)*dWdz.rho);
    
    heat_flux = -kappa()*gradient_T;
    
    return (heat_flux);
       
 }

// viscous flux
NavierStokes3D_ThermallyPerfect_cState  NavierStokes3D_ThermallyPerfect_pState::FluxViscous_n(
   const NavierStokes3D_ThermallyPerfect_pState &Wl,
   const NavierStokes3D_ThermallyPerfect_pState &Wr,
   const NavierStokes3D_ThermallyPerfect_pState &Wc,
   const NavierStokes3D_ThermallyPerfect_pState &Wc_Neigbor,
   const NavierStokes3D_ThermallyPerfect_pState &dWdx,
   const NavierStokes3D_ThermallyPerfect_pState &dWdy,
   const NavierStokes3D_ThermallyPerfect_pState &dWdz,
   const NavierStokes3D_ThermallyPerfect_pState &dWdx_Neigbor,
   const NavierStokes3D_ThermallyPerfect_pState &dWdy_Neigbor,
   const NavierStokes3D_ThermallyPerfect_pState &dWdz_Neigbor,
   const Vector3D &norm, const Vector3D &ts, const double &deltad,
   const double &Volume, const double &Volume_Neigbor){
   
   // construct the gradients on the cell interface (surface) 
   // based on Hybrid Average Gradient-Diamond-Path Approach
   // Grad Phi = (Phi_c_neigbor - Phi_c)/ds *norm/(norm . ts)
   //            + [bar (Grad Phi) - bar (Grad Phi) . ts norm/(norm.ts)]

   // weighted factor based on volume
   double alpha = Volume/(Volume + Volume_Neigbor);
  
   
   NavierStokes3D_ThermallyPerfect_pState dWdx_Weighted, 
      dWdy_Weighted, dWdz_Weighted, dWdx_face, 
      dWdy_face, dWdz_face, Grad_middle_term;

   NavierStokes3D_ThermallyPerfect_pState W_face;
   
   dWdx_Weighted = alpha*dWdx + (1.0 - alpha)*dWdx_Neigbor;
   dWdy_Weighted = alpha*dWdy + (1.0 - alpha)*dWdy_Neigbor;
   dWdz_Weighted = alpha*dWdz + (1.0 - alpha)*dWdz_Neigbor;
   
   
   // a weighted term  
   Grad_middle_term = dWdx_Weighted*ts.x + dWdy_Weighted*ts.y +
      dWdz_Weighted*ts.z;
      
   // gradients of primitive variables on the face
   dWdx_face = (Wc_Neigbor - Wc)/deltad *norm.x/dot(norm, ts) + (
      dWdx_Weighted -  Grad_middle_term*norm.x/dot(norm, ts));
   dWdy_face = (Wc_Neigbor - Wc)/deltad *norm.y/dot(norm, ts) + (
      dWdy_Weighted -  Grad_middle_term*norm.y/dot(norm, ts));
   dWdz_face = (Wc_Neigbor - Wc)/deltad *norm.z/dot(norm, ts) + (
      dWdz_Weighted -  Grad_middle_term*norm.z/dot(norm, ts));
   
   W_face = HALF*(Wl + Wr);
   
   
   return (W_face.Fv(dWdx_face, dWdy_face, dWdz_face)*norm.x +
           W_face.Gv(dWdx_face, dWdy_face, dWdz_face)*norm.y +
           W_face.Hv(dWdx_face, dWdy_face, dWdz_face)*norm.z);
   


   
   
   
}//endofgradients

  
