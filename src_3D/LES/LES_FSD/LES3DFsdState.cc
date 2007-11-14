/********************** LES3DFsdState.cc **********************
  Constructors for LES3Dstate which handles all the physical
  variables and mixture rules associated with multi-species
  chemically reacting flows.

   assosicated files:
           LES3DFsdState.h    

***************************************************************/

#ifndef _LES3DFSD_STATE_INCLUDED
#include "LES3DFsdState.h"
#endif // LES3DFSD_STATE_INCLUDED   

/************* set_species_data ***************************/

double LES3DFsd_pState::Mref = 0.1;
double LES3DFsd_cState::Mref = 0.1;

/****************************************************
  Speed of sound using 
  a^2 = dip/dirho + p/rho^2( die/dip)^-1
  from eigenvalue analysis using e =f(p,rho)
****************************************************/

double LES3DFsd_pState::a(void){
   double sum;
   sum = sqr(Euler3D_ThermallyPerfect_pState::a());
   sum += 2.0/3.0*k*Euler3D_ThermallyPerfect_pState::g();
   return sqrt(sum);
}

double LES3DFsd_pState::a(void) const{
   double sum;
   sum = sqr(Euler3D_ThermallyPerfect_pState::a());
   sum += 2.0/3.0*k*Euler3D_ThermallyPerfect_pState::g();
   return sqrt(sum);
}

void LES3DFsd_pState::set_species_data(const int &n,
                                       const string *S,
                                       const char *PATH,
                                       const int &debug, 
                                       const double &Mr, 
                                       const double* Sc,
                                       const int &trans_data) { 
   ns = n;
   num_vars = NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3D_VAR_EXTRA;
   //read in NASA data for each species to be used 
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
    Mref = Mr;
   // Set data temperature ranges for mixture calculations.
   Temp_low_range(); 
   Temp_high_range(); 
   // Set debug information level.
   debug_level = debug;
}   

void LES3DFsd_cState::set_species_data(const int &n, 
                                       const string *S, 
                                       const char *PATH,
                                       const int &debug, 
                                       const double &Mr, 
                                       const double* Sc,
                                       const int &trans_data) { 
   ns = n;
   num_vars = NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3D_VAR_EXTRA;
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
   Mref = Mr;
   // Set data temperature ranges for mixture calculations.
   Temp_low_range();
   Temp_high_range();
   // Set debug information level.
   debug_level = debug;
}      

/**********************************************************
     LES3DFsd_pState::F -- Inviscid flux (x-direction).   
***********************************************************/
LES3DFsd_cState LES3DFsd_pState::F(void) {
   // call the inviscid flux calculation in Euler3D_ThermallyPerfect_pState class  
   Euler3D_ThermallyPerfect_cState Temp;
   Temp = Euler3D_ThermallyPerfect_pState::F();
   Temp.rhov.x += 2.0/3.0*rho*k;
   Temp.E += v.x*(2.0/3.0*rho*k);
   double Temp_rhoC = rho*v.x*C;
   double Temp_rhoFsd = rho*v.x*Fsd;
   double Temp_rhok = rho*v.x*k;
   return (LES3DFsd_cState(Temp.rho, Temp.rhov, Temp.E, Temp_rhoC, Temp_rhoFsd, Temp_rhok ));//,Temp.rhospec ));
}

LES3DFsd_cState LES3DFsd_pState::F(void) const {
   // call the inviscid flux calculation in Euler3D_ThermallyPerfect_pState class  
   Euler3D_ThermallyPerfect_cState Temp;
   Temp = Euler3D_ThermallyPerfect_pState::F();
   Temp.rhov.x += 2.0/3.0*rho*k;
   Temp.E += v.x*(2.0/3.0*rho*k);
   double Temp_rhoC = rho*v.x*C;
   double Temp_rhoFsd = rho*v.x*Fsd;
   double Temp_rhok = rho*v.x*k;
   return (LES3DFsd_cState(Temp.rho, Temp.rhov, Temp.E, Temp_rhoC, Temp_rhoFsd, Temp_rhok ));//,Temp.rhospec ));
}

/********************************************************
   LES3DFsd_pState::Fv -- Viscous flux (x-direction).   
 ********************************************************/
LES3DFsd_cState LES3DFsd_pState::Fv(const LES3DFsd_pState &dWdx,
                                    const LES3DFsd_pState &dWdy,
                                    const LES3DFsd_pState &dWdz,
                                    const int &Flow_Type,
                                    const double &Volume) const{
   LES3DFsd_cState Temp;
   Temp.Vacuum();
   double Temperature = T();
   double mu_t = eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,Volume);
   Tensor3D molecular_stress, subfilter_stress;
   Vector3D heat_flux;
   molecular_stress = tau(dWdx,dWdy,dWdz);
   subfilter_stress = lambda(dWdx,dWdy,dWdz,Flow_Type,Volume);
   /*----------------- Thermal Diffusion ---------------------------*/
   heat_flux.x = qflux(dWdx,dWdy,dWdz).x+qflux_t(dWdx,dWdy,dWdz,Flow_Type,Volume).x;
   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux.x -= rho*thermal_diffusion().x;
   // SFS thermal diffusion
   heat_flux.x -= rho*thermal_diffusion_t(dWdx,dWdy,dWdx,Flow_Type,Volume).x;

   Temp.rhov.x = molecular_stress.xx + subfilter_stress.xx + 2.0/3.0*rho*k;
   Temp.rhov.y = molecular_stress.xy + subfilter_stress.xy;
   Temp.rhov.z = molecular_stress.xz + subfilter_stress.xz;
   Temp.E = v.x*(molecular_stress.xx + subfilter_stress.xx + 2.0/3.0*rho*k) +
            v.y*(molecular_stress.xy + subfilter_stress.xy) + 
            v.z*(molecular_stress.xz + subfilter_stress.xz) - heat_flux.x; 
   Temp.rhoC = mu_t*dWdx.C/Sc_t();
   Temp.rhoFsd = mu_t*dWdx.Fsd/Sc_t();
   Temp.rhok = (mu()+mu_t/0.25/*Pr_t()*/)*dWdx.k;
   return (Temp);
}

LES3DFsd_cState LES3DFsd_pState::Gv(const LES3DFsd_pState &dWdx,
                                    const LES3DFsd_pState &dWdy,
                                    const LES3DFsd_pState &dWdz,
                                    const int &Flow_Type,
                                    const double &Volume) const{
   LES3DFsd_cState Temp;
   Temp.Vacuum();
   double Temperature = T();
   double mu_t = eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,Volume);
   Tensor3D molecular_stress, subfilter_stress;
   Vector3D heat_flux;
   molecular_stress = tau(dWdx,dWdy,dWdz);
   subfilter_stress = lambda(dWdx,dWdy,dWdz,Flow_Type,Volume);
   /*----------------- Thermal Diffusion ---------------------------*/
   heat_flux.y = qflux(dWdx,dWdy,dWdz).y+qflux_t(dWdx,dWdy,dWdz,Flow_Type,Volume).y;
   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux.y -= rho*thermal_diffusion().y;
   // SFS thermal diffusion
   heat_flux.y -= rho*thermal_diffusion_t(dWdx,dWdy,dWdz,Flow_Type,Volume).y;

   Temp.rhov.x = molecular_stress.xy + subfilter_stress.xy;
   Temp.rhov.y = molecular_stress.yy + subfilter_stress.yy + 2.0/3.0*rho*k ;
   Temp.rhov.z = molecular_stress.yz + subfilter_stress.yz;
   Temp.E = v.x*(molecular_stress.xy + subfilter_stress.xy) +
            v.y*(molecular_stress.yy + subfilter_stress.yy + 2.0/3.0*rho*k) + 
            v.z*(molecular_stress.yz + subfilter_stress.yz) - heat_flux.y; 
   Temp.rhoC = mu_t*dWdy.C/Sc_t();
   Temp.rhoFsd = mu_t*dWdy.Fsd/Sc_t();
   Temp.rhok = (mu()+mu_t/0.25/*Pr_t()*/)*dWdy.k;
   return (Temp);
}

LES3DFsd_cState LES3DFsd_pState::Hv(const LES3DFsd_pState &dWdx,
                                    const LES3DFsd_pState &dWdy,
                                    const LES3DFsd_pState &dWdz,
                                    const int &Flow_Type,
                                    const double &Volume) const{
   LES3DFsd_cState Temp;
   Temp.Vacuum();
   double Temperature = T();
   double mu_t = eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,Volume);
   Tensor3D molecular_stress, subfilter_stress;
   Vector3D heat_flux;
   molecular_stress = tau(dWdx,dWdy,dWdz);
   subfilter_stress = lambda(dWdx,dWdy,dWdz,Flow_Type,Volume);
   /*----------------- Thermal Diffusion ---------------------------*/
   heat_flux.z = qflux(dWdx,dWdy,dWdz).z+qflux_t(dWdx,dWdy,dWdz,Flow_Type,Volume).z;
   // q -= rho * sum ( hs * Ds *gradcs)   
   heat_flux.z -= rho*thermal_diffusion().z;
   // SFS thermal diffusion
   heat_flux.z -= rho*thermal_diffusion_t(dWdx,dWdy,dWdz,Flow_Type,Volume).z;

   Temp.rhov.x = molecular_stress.xz + subfilter_stress.xz;
   Temp.rhov.y = molecular_stress.yz + subfilter_stress.yz;
   Temp.rhov.z = molecular_stress.zz + subfilter_stress.zz + 2.0/3.0*rho*k;
   Temp.E = v.x*(molecular_stress.xz + subfilter_stress.xz) +
            v.y*(molecular_stress.yz + subfilter_stress.yz) + 
            v.z*(molecular_stress.zz + subfilter_stress.zz + 2.0/3.0*rho*k) - heat_flux.z ;
   Temp.rhoC = mu_t*dWdz.C/Sc_t();
   Temp.rhoFsd = mu_t*dWdz.Fsd/Sc_t();
   Temp.rhok = (mu()+mu_t/0.25/*Pr_t()*/)*dWdz.k;
   return (Temp);
}

Tensor3D LES3DFsd_pState::rotation_tensor(const LES3DFsd_pState &dWdx, 
                                          const LES3DFsd_pState &dWdy, 
                                          const LES3DFsd_pState &dWdz) const{
   Tensor3D rotation;
   rotation.zero();
   rotation.xy = -1.0/2.0*(dWdx.v.y - dWdy.v.x);
   rotation.xz = 1.0/2.0*(dWdz.v.x - dWdx.v.z);
   rotation.yz = -1.0/2.0*(dWdy.v.z - dWdz.v.y);
   return (rotation);
}

double LES3DFsd_cState::a(void) const{
   double sum;
   sum = sqr(Euler3D_ThermallyPerfect_cState::a());
   sum += 2.0/3.0*rhok/rho*Euler3D_ThermallyPerfect_cState::g();
   return sqrt(sum);
}

/***************************************************************
     LES3DFsd_pState::lambda -- Eigenvalue(s) (x-direction).    
 ***************************************************************/
LES3DFsd_pState LES3DFsd_pState::lambda_x(void){
   double c = a();
   return (LES3DFsd_pState(v.x-c, v.x, v.x, v.x, v.x+c, v.x, v.x, v.x));
}

LES3DFsd_pState LES3DFsd_pState::lambda_x(void) const {
   double c = a();
   return (LES3DFsd_pState(v.x-c, v.x, v.x, v.x, v.x+c, v.x, v.x, v.x));
}

/*********************************************************************
              LES3DFsd_pState -- Binary arithmetic operators.        
 *********************************************************************/
LES3DFsd_pState LES3DFsd_pState::operator +(const LES3DFsd_pState &W) const{ 
   LES3DFsd_pState Temp;
   Temp.Copy(*this);
   Temp += W;
   return Temp;
}

//------------------ Subtraction ------------------------//
LES3DFsd_pState LES3DFsd_pState::operator -(const LES3DFsd_pState &W) const{
   LES3DFsd_pState Temp;
   Temp.Copy(*this);
   Temp -= W;
   return Temp;
}

//---------------- Scalar Multiplication ------------------//
LES3DFsd_pState LES3DFsd_pState::operator *(const double &a) const{
   LES3DFsd_pState Temp;
   Temp.Copy(*this);
   Temp.rho = rho*a;  Temp.v = v*a; Temp.p = p*a;
   Temp.C = C*a; Temp.Fsd = Fsd*a; Temp.k = k*a; 
   for( int i=0; i<ns; i++){
      Temp.spec[i] = spec[i]*a;
   } 
   return(Temp);
}

LES3DFsd_pState operator *(const double &a, 
                           const LES3DFsd_pState &W){
   LES3DFsd_pState Temp;
   //Temp.Copy(W);
   Temp.rho = W.rho*a;  Temp.v = W.v*a; Temp.p = W.p*a;
   Temp.C = W.C*a; Temp.Fsd = W.Fsd*a; Temp.k = W.k*a; 
   for( int i=0; i<W.ns; i++){
      Temp.spec[i] = W.spec[i]*a;
   } 
   return(Temp);
}

//--------------- Scalar Division ------------------------//
LES3DFsd_pState LES3DFsd_pState::operator /(const double &a) const {
  LES3DFsd_pState Temp(rho,v,p,C,Fsd,k);
   Temp.Copy(*this);
   Temp.rho = rho/a; Temp.v = v/a; Temp.p = p/a; 
   Temp.C = C/a; Temp.Fsd = Fsd/a; Temp.k = k/a; 
   for(int i=0; i<ns; i++){
      Temp.spec[i] = spec[i]/a; 
   } 
   return(Temp);
}

//----------------- Inner Product ------------------------//
double LES3DFsd_pState::operator *(const LES3DFsd_pState &W) const{
   double sum=0.0;
   for(int i=0; i<ns; i++){
      sum += spec[i]*W.spec[i];
   }  
   return (rho*W.rho + v*W.v + p*W.p + C*W.C + Fsd*W.Fsd + k*W.k + sum);
}

//----------- solution state product operator ------------//
LES3DFsd_pState LES3DFsd_pState::operator ^(const LES3DFsd_pState &W) const{
  LES3DFsd_pState Temp(rho,v,p,C,Fsd,k);
   Temp.Copy(*this);
   Temp.rho = rho*W.rho;
   Temp.v.x = v.x*W.v.x;
   Temp.v.y = v.y*W.v.y;
   Temp.v.z = v.z*W.v.z;
   Temp.p = p*W.p;
   Temp.C = C*W.C;
   Temp.Fsd = Fsd*W.Fsd;
   Temp.k = k*W.k;
   for(int i=0; i<ns; i++){
      Temp.spec[i] = spec[i]*W.spec[i];
   }  
   return(Temp);
}

//----------------- Assignment ----------------------------//
LES3DFsd_pState& LES3DFsd_pState::operator =(const LES3DFsd_pState &W){
   //self assignment protection
   if( this != &W){   
      //copy assignment
      rho = W.rho;
      v = W.v; 
      p = W.p; 
      C = W.C;
      Fsd = W.Fsd;
      k = W.k;
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

/**********************************************************************
            LES3DFsd_pState -- Shortcut arithmetic operators.           
 **********************************************************************/
LES3DFsd_pState& LES3DFsd_pState::operator +=(const LES3DFsd_pState &W){
   rho += W.rho;
   v += W.v; 
   p += W.p; 
   C += W.C;
   Fsd += W.Fsd;
   k += W.k;
   for( int i=0; i<ns; i++){
      spec[i] += W.spec[i];
   } 
   return (*this);
}

LES3DFsd_pState& LES3DFsd_pState::operator -=(const LES3DFsd_pState &W) {
   rho -= W.rho;
   v -= W.v;
   p -= W.p;
   C -= W.C;
   Fsd -= W.Fsd;
   k -= W.k;
   for(int i=0; i<ns; i++){
      spec[i] -= W.spec[i];
   }  
   return (*this); 
}

/****************************************************************
              LES3DFsd_pState -- Relational operators.          
 ****************************************************************/
int operator ==(const LES3DFsd_pState &W1, 
                const LES3DFsd_pState &W2) {
   if(W1.ns == W2.ns){ //check that species are equal
      bool Temp;
      for(int i=0; i<W1.ns; i++){
         if( W1.spec[i] == W2.spec[i] ){
            Temp = true;
         } else {
            Temp = false;
            break;
         }  
         return (W1.rho == W2.rho && W1.v == W2.v && W1.p == W2.p &&
                 W1.C == W2.C && W1.Fsd == W2.Fsd && W1.k == W2.k &&
                 Temp == true);
      }
   } else {
      cerr<<"\n Mismatch in number of species \n";
      exit(1);
   }
}

int operator !=(const LES3DFsd_pState &W1, 
                const LES3DFsd_pState &W2) {
   if(W1.ns == W2.ns){ //check that species are equal
      bool Temp = true;
      for(int i=0; i<W1.ns; i++){
         if( W1.spec[i] != W2.spec[i] ){
            Temp = false;
            break;
         } 
         return (W1.rho != W2.rho || W1.v != W2.v || W1.p != W2.p ||
                 W1.C != W2.C || W1.Fsd != W2.Fsd || W1.k != W2.k ||  
                 Temp != true);
      }
   } else {
      cerr<<"\n Mismatch in number of species \n";
      exit(1);
   }
}

/* LES3DFsd_pState -- Input-output operators.*/
ostream &operator << (ostream &out_file, 
                      const LES3DFsd_pState &W) {
   out_file.precision(10);
   out_file.setf(ios::scientific);
   out_file <<" "<<W.rho<<" "<<W.v.x<<" "<<W.v.y<<" "<< W.v.z
            <<" "<<W.p<<" "<<W.C<<" "<<W.Fsd<<" "<<W.k;
   for(int i=0; i<W.ns; i++){
      out_file<<" "<<W.spec[i];
   }
   out_file.unsetf(ios::scientific);
   return (out_file);
}

istream &operator >> (istream &in_file, LES3DFsd_pState &W) {
   in_file.setf(ios::skipws);
   in_file >> W.rho >> W.v.x >> W.v.y >> W.v.z >> W.p >> W.C >> W.Fsd >> W.k;
   //W.set_initial_values();
   for( int i=0; i<W.ns; i++){
      in_file>>W.spec[i];
   }
   in_file.unsetf(ios::skipws);
   return (in_file);
}


//----------------- Addition ----------------------------//
LES3DFsd_cState LES3DFsd_cState::operator +(const LES3DFsd_cState &U) const{ 
   LES3DFsd_cState Temp;
   Temp.Copy(*this);
   Temp += U;
   return Temp;
}

//------------------ Subtraction ------------------------//
LES3DFsd_cState LES3DFsd_cState::operator -(const LES3DFsd_cState &U) const{
   LES3DFsd_cState Temp;
   Temp.Copy(*this);
   Temp -= U;
   return Temp;
}

//---------------- Scalar Multiplication ------------------//
LES3DFsd_cState LES3DFsd_cState::operator *(const double &a) const{
   LES3DFsd_cState Temp;
   Temp.Copy(*this);
   Temp.rho = rho*a; Temp.rhov = rhov*a; Temp.E = E*a;
   Temp.rhoC = rhoC*a; Temp.rhoFsd = rhoFsd*a; Temp.rhok = rhok*a;
   for( int i=0; i<ns; i++){
      Temp.rhospec[i] = rhospec[i]*a;
   } 
   return(Temp);
}

LES3DFsd_cState operator *(const double &a, const LES3DFsd_cState &U){
   LES3DFsd_cState Temp;
   Temp.rho = U.rho*a;  Temp.rhov = U.rhov*a; Temp.E = U.E*a;
   Temp.rhoC = U.rhoC*a; Temp.rhoFsd = U.rhoFsd*a; Temp.rhok = U.rhok*a;
   for( int i=0; i<U.ns; i++){
      Temp.rhospec[i] = U.rhospec[i]*a;
   } 
   return(Temp);
}
//--------------- Scalar Division ------------------------//
LES3DFsd_cState LES3DFsd_cState::operator /(const double &a) const {
   LES3DFsd_cState Temp;
   Temp.Copy(*this);
   Temp.rho = rho/a; Temp.rhov = rhov/a; Temp.E = E/a;
   Temp.rhoC = rhoC/a; Temp.rhoFsd = rhoFsd/a; Temp.rhok = rhok/a;
   for(int i=0; i<ns; i++){
      Temp.rhospec[i] = rhospec[i]/a; 
   } 
   return(Temp);
}

//----------------- Inner Product ------------------------//
double LES3DFsd_cState::operator *(const LES3DFsd_cState &U) const{
   double sum=0.0;
   for(int i=0; i<ns; i++){
      sum += rhospec[i]*U.rhospec[i];
   }  
   return (rho*U.rho+rhov*U.rhov+E*U.E+rhoC*U.rhoC+rhoFsd*U.rhoFsd+rhok*U.rhok+sum);
}

//----------- solution state product operator ------------//
LES3DFsd_cState LES3DFsd_cState::operator ^(const LES3DFsd_cState &U) const {
   LES3DFsd_cState Temp;
   Temp.Copy(*this);
   Temp.rho = rho*U.rho;
   Temp.rhov.x = rhov.x*U.rhov.x;
   Temp.rhov.y = rhov.y*U.rhov.y;
   Temp.rhov.z = rhov.z*U.rhov.z;
   Temp.E = E*U.E;
   Temp.rhoC = rhoC*U.rhoC;
   Temp.rhoFsd = rhoFsd*U.rhoFsd;
   Temp.rhok = rhok*U.rhok;
   for(int i=0; i<ns; i++){
      Temp.rhospec[i] = rhospec[i]*U.rhospec[i];
   }  
   return(Temp);
}

//----------------- Assignment ----------------------------//
LES3DFsd_cState& LES3DFsd_cState::operator =(const LES3DFsd_cState &U){
   //self assignment protection
   if( this != &U){   
      //copy assignment
      rho = U.rho;
      rhov = U.rhov; 
      E = U.E; 
      rhoC = U.rhoC;
      rhoFsd = U.rhoFsd;
      rhok = U.rhok;
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

/*******************************************************************
         LES3DFsd_cState -- Shortcut arithmetic operators.          
 *******************************************************************/
LES3DFsd_cState& LES3DFsd_cState::operator +=(const LES3DFsd_cState &U){
   rho += U.rho;
   rhov += U.rhov; 
   E += U.E;
   rhoC += U.rhoC;
   rhoFsd += U.rhoFsd;
   rhok += U.rhok;
   for( int i=0; i<ns; i++){
      rhospec[i] += U.rhospec[i];
   } 
   return (*this);
}

LES3DFsd_cState& LES3DFsd_cState::operator -=(const LES3DFsd_cState &U) {
   rho -= U.rho;
   rhov -= U.rhov;
   E -= U.E;
   rhoC -= U.rhoC;
   rhoFsd -= U.rhoFsd;
   rhok -= U.rhok;
  for(int i=0; i<ns; i++){
    rhospec[i] -= U.rhospec[i];
  }  
  return (*this); 
}

/********************************************************************
           LES3DFsd_cState -- Unary arithmetic operators.            
 ********************************************************************/
LES3DFsd_cState operator -(const LES3DFsd_cState &U) {
   Species *spt= new Species[U.ns];
   for(int i=0; i<U.ns; i++){
      spt[i] = -U.rhospec[i]; 
   }  
   LES3DFsd_cState Temp(-U.rho,-U.rhov,-U.E,-U.rhoC,-U.rhoFsd,-U.rhok,spt);
  delete[] spt;
  return(Temp);
}

/****************************************************************
            LES3DFsd_cState -- Relational operators.            
 ****************************************************************/
int operator ==(const LES3DFsd_cState &U1, 
                const LES3DFsd_cState &U2) {
   if(U1.ns == U2.ns){ //check that species are equal
      bool Temp;
      for(int i=0; i<U1.ns; i++){
         if( U1.rhospec[i] == U2.rhospec[i] ){
            Temp = true;
         } else {
            Temp = false;
            break;
         }  
         return (U1.rho == U2.rho && U1.rhov == U2.rhov && U1.E == U2.E &&
                 U1.rhoC == U2.rhoC && U1.rhoFsd == U2.rhoFsd && U1.rhok == U2.rhok && 
                 Temp == true);
      }
   } else {
      cerr<<"\n Mismatch in number of species \n";
      exit(1);
   }
}

int operator !=(const LES3DFsd_cState &U1, 
                const LES3DFsd_cState &U2) {
   if(U1.ns == U2.ns){ //check that species are equal
      bool Temp = true;
      for(int i=0; i<U1.ns; i++){
         if( U1.rhospec[i] != U2.rhospec[i] ){
            Temp = false;
            break;
         } 
         return (U1.rho != U2.rho || U1.rhov != U2.rhov || U1.E != U2.E ||
                 U1.rhoC != U2.rhoC || U1.rhoFsd !=U2.rhoFsd || U1.rhok != U2.rhok || 
                 Temp != true);
      }
   } else {
      cerr<<"\n Mismatch in number of species \n";
      exit(1);
   }
}

/******************************************************************
              LES3DFsd_cState -- Input-output operators.         
 ******************************************************************/
ostream &operator << (ostream &out_file, const LES3DFsd_cState &U) {
  out_file.setf(ios::scientific);
  out_file <<" "<<U.rho<<" "<<U.rhov.x<<" "<<U.rhov.y<<" "<<U.rhov.z
           <<" "<<U.E<<" "<<U.rhoC<<" "<<U.rhoFsd<<" "<<U.rhok;
  for( int i=0; i<U.ns; i++){
     out_file<<" "<<U.rhospec[i];
  } 
  out_file.unsetf(ios::scientific);
  return (out_file);
}

istream &operator >> (istream &in_file, LES3DFsd_cState &U) {
   in_file.setf(ios::skipws);
   in_file >>U.rho>>U.rhov.x>>U.rhov.y>>U.rhov.z>>U.E>>U.rhoC>>U.rhoFsd>>U.rhok;
  for( int i=0; i<U.ns; i++){
    in_file>>U.rhospec[i]; 
  } 
  in_file.unsetf(ios::skipws);
  return (in_file);
}

LES3DFsd_cState LES3DFsd_pState::rc_x(const int &index) const {
   double eta_fsd = Progvar_Species_Grad();  
   switch(index){  
   case 1:
     return (LES3DFsd_cState(ONE, v.x-a(), v.y, v.z, H()/rho-v.x*a(), C, Fsd, k));
   case 2:
     return (LES3DFsd_cState(ONE, v.x, v.y, v.z, H()/rho-a()*a()/(g()-1.0), C, Fsd, k)); 
   case 3:
     return (LES3DFsd_cState(ZERO, ZERO, rho, ZERO, rho*v.y, ZERO, ZERO, ZERO));
   case 4:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, rho, rho*v.z, ZERO, ZERO, ZERO));
   case 5:
     return (LES3DFsd_cState(ONE, v.x+a(), v.y, v.z, H()/rho+v.x*a(), C, Fsd, k));
   case 6:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, rho*eta_fsd, rho, ZERO, ZERO));
   case 7:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, rho, ZERO));
   case 8:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, 5.0*rho/3.0, ZERO, ZERO, rho));
    };
}

LES3DFsd_cState LES3DFsd_pState::rc_x(const int &index) {
   double eta_fsd = Progvar_Species_Grad();  
   switch(index){  
   case 1:
     return (LES3DFsd_cState(ONE, v.x-a(), v.y, v.z, H()/rho-v.x*a(), C, Fsd, k));
   case 2:
     return (LES3DFsd_cState(ONE, v.x, v.y, v.z, H()/rho-a()*a()/(g()-1.0), C, Fsd, k)); 
   case 3:
     return (LES3DFsd_cState(ZERO, ZERO, rho, ZERO, rho*v.y, ZERO, ZERO, ZERO));
   case 4:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, rho, rho*v.z, ZERO, ZERO, ZERO));
   case 5:
     return (LES3DFsd_cState(ONE, v.x+a(), v.y, v.z, H()/rho+v.x*a(), C, Fsd, k));
   case 6:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, rho*eta_fsd, rho, ZERO, ZERO));
   case 7:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, rho, ZERO));
   case 8:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, 5.0*rho/3.0, ZERO, ZERO, rho));
   };
}

// Primitive Left Eigenvector -- (x-direction)
LES3DFsd_pState LES3DFsd_pState::lp_x(const int &index) const {
   double eta_fsd = Progvar_Species_Grad();  
   switch(index){  
   case 1:
      return (LES3DFsd_pState(ZERO, -HALF*rho/a(), ZERO, ZERO, HALF/(a()*a()), ZERO, ZERO, ZERO));
   case 2:
      return (LES3DFsd_pState(ONE, ZERO, ZERO, ZERO, -ONE/(a()*a()), ZERO, ZERO, ZERO));
   case 3 :
      return (LES3DFsd_pState(ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO, ZERO));
   case 4 :
      return (LES3DFsd_pState(ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO));
   case 5 :
      return (LES3DFsd_pState(ZERO, HALF*rho/a(), ZERO, ZERO, HALF/(a()*a()), ZERO, ZERO, ZERO));
   case 6 :
     return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO, ZERO));
   case 7 :
     return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO));
   case 8 :
     return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE));
   };
}

LES3DFsd_pState LES3DFsd_pState::lp_x(const int &index) {
   double eta_fsd = Progvar_Species_Grad();  
   switch(index){  
   case 1:
      return (LES3DFsd_pState(ZERO, -HALF*rho/a(), ZERO, ZERO, HALF/(a()*a()), ZERO, ZERO, ZERO));
   case 2:
      return (LES3DFsd_pState(ONE, ZERO, ZERO, ZERO, -ONE/(a()*a()), ZERO, ZERO, ZERO));
   case 3 :
      return (LES3DFsd_pState(ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO, ZERO));
   case 4 :
      return (LES3DFsd_pState(ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO));
   case 5 :
       return (LES3DFsd_pState(ZERO, HALF*rho/a(), ZERO, ZERO, HALF/(a()*a()), ZERO, ZERO, ZERO));
   case 6 :
     return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO, ZERO));
   case 7 :
     return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO));
   case 8 :
     return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE));
   };
}

// as defined by E.Turkel (1999)
double LES3DFsd_pState::Mr2(const double &deltax, const double &lengthx, const double &dTime) const {
  double c = a();
  double MR2 = min(max((v.sqr()/(c*c)),Mref*Mref),ONE);
  // Need deltax which is based on cell spacing 
    MR2 = pow(max(sqrt(MR2*c*c), mu()/(rho*deltax)),2.0)/(c*c);
    double MR_uns = (lengthx/(PI*dTime))/c;  // ZERO;
// double MR_vis;
//     if (flow_type_flag == FLOWTYPE_INVISCID) {
//       MR_vis = ZERO;
//     } else {
//       MR_vis = (mu()/(rho*deltax))/c;
//     }  
   //  MR2 = max(MR_inv*MR_inv, MR_vis*MR_vis);
    MR2 = max(MR2, MR_uns*MR_uns);
    MR2 = min(max(MR2, Mref*Mref), ONE);
    //cout << "L = " << lengthx <<" dTime = " << dTime << endl;
    //cout << "PI = " << PI <<" MR_uns = "<< MR_uns <<" MR2 = " << MR2 << endl;
  return (MR2);
}
/*******************************************************************
    Low Mach Number Preconditioner Based on Weiss & Smith (1995)   
 *******************************************************************/
// For CFL Calculation
double LES3DFsd_pState::u_plus_aprecon(const double &u, double &deltax, const double &lengthx, const double &dTime) const {
  double Temp = T();
  double c = a();
  double UR2 = Mr2(deltax, lengthx, dTime)*c*c;
  double alpha = HALF*( ONE - (ONE/(Rtot()*Temp) - ONE/(Cp()*Temp))*UR2);
  // uprime + cprime
  return ( u*(ONE - alpha) + sqrt(alpha*alpha*u*u + UR2) );
}

// For eigenvalues and eigenvectors return preconditioned velocity and soundspeed ie. u' and c'
void LES3DFsd_pState::u_a_precon(const double &UR2, double &uprimed, double &cprimed) const{
//   double Temp = T();
//   double alpha = HALF*( ONE - (ONE/(Rtot()*Temp) - ONE/(Cp()*Temp))*UR2);
  double alpha = HALF*(ONE - ONE/(a()*a())*UR2);
  uprimed = v.x*(ONE - alpha);
  cprimed = sqrt(alpha*alpha*v.x*v.x + UR2); 
}

LES3DFsd_cState LES3DFsd_pState::rc_x_precon(const int &index, const double &MR2) const {
    double uprimed,cprimed;
    double c = a();
    double eta_fsd = Progvar_Species_Grad();
    u_a_precon(MR2*a()*a(),uprimed,cprimed);
   switch(index){  
   case 1:
     return (LES3DFsd_cState(ONE, (uprimed-cprimed)/MR2, v.y, v.z, h()+(v.sqr()/MR2)/TWO+FIVE*k/THREE-(v.x*cprimed)/MR2, C, Fsd, k));
   case 2:
     return (LES3DFsd_cState(ONE, v.x, v.y, v.z, (H()/rho-c*c/(g()-ONE)), C, Fsd, k)); 
   case 3:
     return (LES3DFsd_cState(ZERO, ZERO, rho, ZERO, rho*v.y, ZERO, ZERO, ZERO));
   case 4:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, rho, rho*v.z, ZERO, ZERO, ZERO));
   case 5:
     return (LES3DFsd_cState(ONE, (uprimed+cprimed)/MR2, v.y, v.z, h()+(v.sqr()/MR2)/TWO+FIVE*k/THREE+(v.x*cprimed)/MR2, C, Fsd, k));
   case 6:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, rho*eta_fsd, rho, ZERO, ZERO));
   case 7:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, rho, ZERO));
   case 8:
     return (LES3DFsd_cState(ZERO, ZERO, ZERO, ZERO, FIVE*rho/THREE, ZERO, ZERO, rho));
   };
}

LES3DFsd_pState LES3DFsd_pState::lp_x_precon(const int &index,const double &MR2) const {
    double uprimed,cprimed;
    double c = a();
    u_a_precon(MR2*a()*a(),uprimed,cprimed);
   switch(index){  
   case 1:
      return (LES3DFsd_pState(ZERO, -HALF*rho*MR2/cprimed, ZERO, ZERO, (-uprimed+cprimed+v.x)/(TWO*cprimed*c*c), ZERO, ZERO, ZERO));
   case 2:
      return (LES3DFsd_pState(ONE, ZERO, ZERO, ZERO, -ONE/(c*c), ZERO, ZERO, ZERO));
   case 3 :
      return (LES3DFsd_pState(ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO, ZERO));
   case 4 :
      return (LES3DFsd_pState(ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO));
   case 5 :
      return (LES3DFsd_pState(ZERO, HALF*rho*MR2/cprimed, ZERO, ZERO, (uprimed+cprimed-v.x)/(TWO*cprimed*c*c), ZERO, ZERO, ZERO));
   case 6 :
      return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO, ZERO));
   case 7 :
      return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO));
   case 8 :
      return (LES3DFsd_pState(ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE));
   };
}

void LES3DFsd_pState::Low_Mach_Number_Preconditioner(DenseMatrix &P, const double &deltax, 
                                                  const double &lengthx,const double &dTime) const{  
  double Temp = T();
  double Rmix = Rtot();
  double enthalpy = h();
  double CP = Cp();
  double c = a();
  double pt = pmodified();
  double theta = (ONE/(Mr2(deltax,lengthx,dTime)*c*c) + (g()-ONE)/(c*c));// + ONE/(CP*Temp));
  double eta_fsd = Progvar_Species_Grad();
  double phi = C*eta_fsd;
  double alpha = theta*pt/rho;
  double alpham1 = alpha - ONE;
  double Omega = (Rmix - CP)*pt/(rho*Rmix);
  double beta = enthalpy - CP*pt/(rho*Rmix) - phi;
  double V = HALF*v.sqr();
  P.zero();

  P(0,0) = (alpha*(beta-V)+pt/rho+V+Rmix*Temp-enthalpy+phi)/Omega;
  P(0,1) = v.x*alpham1/Omega;
  P(0,2) = v.y*alpham1/Omega;
  P(0,3) = v.z*alpham1/Omega;
  P(0,4) = -alpham1/Omega;
  P(0,5) = eta_fsd*alpham1/Omega;
  P(0,7) = 5.0*alpham1/(3.0*Omega);
  P(1,0) = v.x*(beta-V)*alpham1/Omega;
  P(1,1) = v.x*v.x*alpham1/Omega+1.0;
  P(1,2) = v.x*v.y*alpham1/Omega;
  P(1,3) = v.x*v.z*alpham1/Omega;
  P(1,4) = -v.x*alpham1/Omega;
  P(1,5) = v.x*eta_fsd*alpham1/Omega;
  P(1,7) = 5.0*v.x*alpham1/(3.0*Omega);
  P(2,0) = v.y*(beta-V)*alpham1/Omega;
  P(2,1) = v.x*v.y*alpham1/Omega;
  P(2,2) = v.y*v.y*alpham1/Omega+1.0;
  P(2,3) = v.y*v.z*alpham1/Omega;
  P(2,4) = -v.y*alpham1/Omega;
  P(2,5) = v.y*eta_fsd*alpham1/Omega;
  P(2,7) = 5.0*v.y*alpham1/(3.0*Omega);
  P(3,0) = v.z*(beta-V)*alpham1/Omega;
  P(3,1) = v.x*v.z*alpham1/Omega;
  P(3,2) = v.y*v.z*alpham1/Omega+1.0;
  P(3,3) = v.z*v.z*alpham1/Omega;
  P(3,4) = -v.z*alpham1/Omega;
  P(3,5) = v.z*eta_fsd*alpham1/Omega;
  P(3,7) = 5.0*v.z*alpham1/(3.0*Omega);
  P(4,0) = (enthalpy+V+5.0*k/3.0)*(beta-V)*alpham1/Omega;
  P(4,1) = v.x*(enthalpy+V+5.0*k/3.0)*alpham1/Omega;
  P(4,2) = v.y*(enthalpy+V+5.0*k/3.0)*alpham1/Omega;
  P(4,3) = v.z*(enthalpy+V+5.0*k/3.0)*alpham1/Omega;
  P(4,4) = -(alpha*(enthalpy+V+5.0*k/3.0)-V-5.0*k/3.0-Rmix*Temp-beta-phi)/Omega;
  P(4,5) = eta_fsd*(enthalpy+V+5.0*k/3.0)*alpham1/Omega;
  P(4,7) = 5.0*(enthalpy+V+5.0*k/3.0)*alpham1/(3.0*Omega);
  P(5,0) = C*(beta-V)*alpham1/Omega;
  P(5,1) = C*v.x*alpham1/Omega;
  P(5,2) = C*v.y*alpham1/Omega;
  P(5,3) = C*v.z*alpham1/Omega;
  P(5,4) = -C*alpham1/Omega;
  P(5,5) = C*eta_fsd*alpham1/Omega+1.0;
  P(5,7) = 5.0*C*alpham1/(3.0*Omega);
  P(6,0) = Fsd*(beta-V)*alpham1/Omega;
  P(6,1) = Fsd*v.x*alpham1/Omega;
  P(6,2) = Fsd*v.y*alpham1/Omega;
  P(6,3) = Fsd*v.z*alpham1/Omega;
  P(6,4) = -Fsd*alpham1/Omega;
  P(6,5) = Fsd*eta_fsd*alpham1/Omega;
  P(6,6) = ONE;
  P(6,7) = 5.0*Fsd*alpham1/(3.0*Omega);
  P(7,0) = k*(beta-V)*alpham1/Omega;
  P(7,1) = k*v.x*alpham1/Omega;
  P(7,2) = k*v.y*alpham1/Omega;
  P(7,3) = k*v.z*alpham1/Omega;
  P(7,4) = -k*alpham1/Omega;
  P(7,5) = k*eta_fsd*alpham1/Omega;
  P(7,7) = 5.0*k*alpham1/(3.0*Omega);
}

void LES3DFsd_pState::Low_Mach_Number_Preconditioner_Inverse(DenseMatrix &Pinv, const double &deltax, 
                                                             const double &lengthx,const double &dTime) const{  
  double Temp = T();
  double Rmix = Rtot();
  double enthalpy = h();
  double CP = Cp();
  double c = a();
  double pt = pmodified();
  double theta = (ONE/(Mr2(deltax,lengthx,dTime)*c*c) + (g()-ONE)/(c*c));// + ONE/(CP*Temp));  

  double eta_fsd = Progvar_Species_Grad();
  double phi = C*eta_fsd;
  double AA = pt*(rho*Rmix-theta*pt*CP);
  double BB = Rmix*rho*(theta*pt-rho);
  double EE = HALF*v.sqr() - enthalpy + phi;
  double CC = EE + CP*pt/(rho*Rmix); 
  double DD = HALF*v.sqr() + enthalpy;
  Pinv.zero();    

  Pinv(0,0) = rho*Rmix*(theta*pt*EE-rho*CC+pt)/AA;
  Pinv(0,1) = -v.x*BB/AA;
  Pinv(0,2) = -v.y*BB/AA;
  Pinv(0,3) = -v.z*BB/AA;
  Pinv(0,4) = BB/AA;
  Pinv(0,5) = -eta_fsd*BB/AA;
  Pinv(0,7) = -5.0*BB/(3.0*AA);
  Pinv(1,0) = v.x*CC*BB/AA;
  Pinv(1,1) = rho*Rmix/AA*(pt+rho*v.x*v.x-theta*pt*(v.x*v.x+CP*Temp));
  Pinv(1,2) = -v.x*v.y*BB/AA;
  Pinv(1,3) = -v.x*v.z*BB/AA;
  Pinv(1,4) = v.x*BB/AA;    
  Pinv(1,5) = -v.x*eta_fsd*BB/AA;
  Pinv(1,7) = -5.0*v.x*BB/(3.0*AA);
  Pinv(2,0) = v.y*CC*BB/AA;
  Pinv(2,1) = -v.x*v.y*BB/AA;
  Pinv(2,2) = rho*Rmix/AA*(pt+v.y*v.y*rho-theta*pt*(v.y*v.y+CP*Temp));
  Pinv(2,3) = -v.z*v.y*BB/AA;
  Pinv(2,4) = v.y*BB/AA;  
  Pinv(2,5) = -v.y*eta_fsd*BB/AA;
  Pinv(2,7) = -5.0*v.y*BB/(3.0*AA);
  Pinv(3,0) = v.z*CC*BB/AA;
  Pinv(3,1) = -v.x*v.z*BB/AA;
  Pinv(3,2) = -v.y*v.z*BB/AA;
  Pinv(3,3) = rho*Rmix/AA*(pt+v.z*v.z*rho-theta*pt*(v.z*v.z+CP*Temp));
  Pinv(3,4) = v.z*BB/AA;  
  Pinv(3,5) = -v.z*eta_fsd*BB/AA;
  Pinv(3,7) = -5.0*v.z*BB/(3.0*AA);
  Pinv(4,0) = DD*CC*BB/AA;
  Pinv(4,1) = -v.x*DD*BB/AA;
  Pinv(4,2) = -v.y*DD*BB/AA;
  Pinv(4,3) = -v.z*DD*BB/AA;
  Pinv(4,4) = rho*Rmix/AA*(theta*pt*(DD-CP*Temp)-rho*DD+pt);
  Pinv(4,5) = -DD*eta_fsd*BB/AA;
  Pinv(4,7) = -5.0*DD*BB/(3.0*AA);
  Pinv(5,0) = C*CC*BB/AA;
  Pinv(5,1) = -C*v.x*BB/AA;
  Pinv(5,2) = -C*v.y*BB/AA;
  Pinv(5,3) = -C*v.z*BB/AA;
  Pinv(5,4) = C*BB/AA;
  Pinv(5,5) = 1.0 - C*eta_fsd*BB/AA;
  Pinv(5,7) = -5.0*C*BB/(3.0*AA);
  Pinv(6,0) = Fsd*CC*BB/AA;
  Pinv(6,1) = -Fsd*v.x*BB/AA;
  Pinv(6,2) = -Fsd*v.y*BB/AA;
  Pinv(6,3) = -Fsd*v.z*BB/AA;
  Pinv(6,4) = Fsd*BB/AA;
  Pinv(6,5) = -Fsd*eta_fsd*BB/AA;
  Pinv(6,6) = ONE;
  Pinv(6,7) = -5.0*Fsd*BB/(3.0*AA);
  Pinv(7,0) = k*CC*BB/AA;
  Pinv(7,1) = -k*v.x*BB/AA;
  Pinv(7,2) = -k*v.y*BB/AA;
  Pinv(7,3) = -k*v.z*BB/AA;
  Pinv(7,4) = k*BB/AA;
  Pinv(7,5) = -k*eta_fsd*BB/AA;
  Pinv(7,7) = ONE-5.0*k/(3.0*AA);
}

void LES3DFsd_pState::dWdU(DenseMatrix &dWdU){
  dWdU.zero();
  double alpha = Rtot()/(Cp()-Rtot());
  double eta_fsd = Progvar_Species_Grad();      

  dWdU(0,0) = ONE;
  dWdU(1,0) = -v.x/rho;
  dWdU(1,1) = ONE/rho;
  dWdU(2,0) = -v.y/rho;
  dWdU(2,2) = -ONE/rho;
  dWdU(3,0) = -v.z/rho;
  dWdU(3,3) = -ONE/rho;
  dWdU(4,0) = alpha*(v.sqr()/TWO-h()+Cp()*T()+C*eta_fsd);
  dWdU(4,1) = -v.x*alpha;
  dWdU(4,2) = -v.y*alpha;
  dWdU(4,3) = -v.z*alpha;
  dWdU(4,4) = alpha;
  dWdU(4,5) = -eta_fsd*alpha;
  dWdU(4,7) = -FIVE*alpha/THREE;
  dWdU(5,0) = -C/rho;
  dWdU(5,5) = ONE/rho;
  dWdU(6,0) = -Fsd/rho;
  dWdU(6,6) = ONE/rho;
  dWdU(7,0) = -k/rho;
  dWdU(7,7) = ONE/rho;
}

void LES3DFsd_pState::SemiImplicitSourceJacobi(const LES3DFsd_pState &dWdx, 
                                               const LES3DFsd_pState &dWdy,
                                               const LES3DFsd_pState &dWdz,
                                               const double &d_dWdx_dW, 
                                               const double &d_dWdy_dW,
                                               const double &d_dWdz_dW,
                                               DenseMatrix &dStdW,
                                               const int &Flow_Type,
                                               const double &Volume)const{

     dStdW.zero();
     double filter = filter_width(Volume);
     double Laminar_Flame_Speed = 0.3837;
     double Reactants_Density = 1.13;
     double k_fsd = SFS_Kinetic_Energy_Fsd(dWdx,dWdy,dWdz,Flow_Type,Volume); 
     double kappa_fsd = Efficiency_Function_Fsd(dWdx,dWdy,dWdz,Flow_Type,Volume);   
     double tau_fsd = HeatRelease_Parameter();
     double Cv = 0.086, Cs = 0.018;
     double beta_fsd = 1.0;

     double t1,t4,t5,t6,t7,t8,t9;
     double t10,t11,t12,t14,t15,t16,t17,t18,t19;
     double t22,t23,t24,t27,t28;
     double t30,t31,t32,t33,t34,t36,t37,t38,t39;
     double t41,t42,t43,t44,t45,t47;
     double t50,t51,t52,t53,t54,t55,t59;
     double t60,t61,t65,t66,t67;
     double t73,t77;
     double t81,t82,t83,t86,t87,t89;
     double t90,t95,t97;
     double t100,t106,t108;
     double t111,t117,t119;
     double t122,t128,t129;
     double t131,t132,t135,t137,t139;
     double t140,t142,t143,t145,t146,t148,t149;
     double t151,t152,t157,t158;
     double t165,t167;
     double t184,t188;
     double t190,t193,t195,t198;
     double t200,t207;
     double t246,t248;
     double t254,t258;
     double t262;
     double t279;
     double t291,t295,t296;
     double t301,t302,t304,t308;
     double t317;
     double t321,t326;
     double t334,t336;
     double t347;
     double t352;
     double t364;
     double t376,t379;
     double t390,t394,t398; 
     double t403;

     t1 = Reactants_Density*Laminar_Flame_Speed;
     t4 = dWdx.C;//cx(c);
      t5 = t4*t4;
      t6 = dWdy.C;//cy(c);
      t7 = t6*t6;
      t8 = dWdz.C;//cz(c);
      t9 = t8*t8;
      t10 = t5+t7+t9;
      t11 = 1/t10;
      t12 = t5*t11;
      t14 = t7*t11;
      t15 = t14/3.0;
      t16 = t9*t11;
      t17 = t16/3.0;
      t18 = 2.0/3.0-2.0/3.0*t12+t15+t17;
      t19 = dWdx.v.x;//Ux(U);
      t22 = t12/3.0;
      t23 = 2.0/3.0-2.0/3.0*t14+t22+t17;
      t24 = dWdy.v.y;//Vy(V);
      t27 = 2.0/3.0-2.0/3.0*t16+t22+t15;
      t28 = dWdz.v.z;//Gz(G);
      t30 = t4*t11;
      if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY){
      t31 = dWdy.v.x;//Uy(U);
      t32 = dWdx.v.y;//Vx(V);
      }else{
      t31 = dWdx.v.y;//Vx(V);
      t32 = dWdy.v.x;//Uy(U);
      }
      t33 = t31+t32;
      t34 = t6*t33;
      if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY){
      t36 = dWdz.v.x;//Uz(U);
      t37 = dWdx.v.z;//Gx(G);
      }else{
      t36 = dWdx.v.z;//Gx(G);
      t37 = dWdz.v.x;//Uz(U);
      }
      t38 = t36+t37;
      t39 = t8*t38;
      t41 = t6*t11;
      if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY){
      t42 = dWdz.v.y;//Vz(V);
      t43 = dWdy.v.z;//Gy(G);
      }else{
      t42 = dWdy.v.z;//Gy(G);
      t43 = dWdz.v.y;//Vz(V);
      }
      t44 = t42+t43;
      t45 = t8*t44;
      t47 = t18*t19+t23*t24+t27*t28-t30*t34-t30*t39-t41*t45;
      t50 = 1.0+tau_fsd*C;//c;
      t51 = sqrt(t10);
      t52 = 1/t51;
      t53 = t4*t52;
      t54 = dWdx.Fsd;//Fsdx(Fsd);
      t55 = d_dWdx_dW;//diff(rhox(rho),rho);
      t59 = t6*t52;
      t60 = dWdy.Fsd;//Fsdy(Fsd);
      t61 = d_dWdy_dW;//diff(rhoy(rho),rho);
      t65 = t8*t52;
      t66 = dWdz.Fsd;//Fsdz(Fsd);
      t67 = d_dWdz_dW;//diff(rhoz(rho),rho);
      t73 = tau_fsd*Fsd;
      t77 = -t5*t52-t7*t52-t9*t52;
      t81 = sqrt(k_fsd);
      t82 = kappa_fsd*t81;
      t83 = 1/filter;
      t86 = beta_fsd*Laminar_Flame_Speed;
      t87 = Fsd*Fsd;
      t89 = 1.0-C;
      t90 = 1/t89;
      t95 = d_dWdx_dW;//diff(Ux(U),U);
      t97 = d_dWdy_dW;//diff(Uy(U),U);
      t100 = d_dWdz_dW;//diff(Uz(U),U);
      t106 = d_dWdy_dW;//diff(Vy(V),V);
      t108 = d_dWdx_dW;//diff(Vx(V),V);
      t111 = d_dWdz_dW;//diff(Vz(V),V);
      t117 = d_dWdz_dW;//diff(Gz(G),G);
      t119 = d_dWdx_dW;//diff(Gx(G),G);
      t122 = d_dWdy_dW;//diff(Gy(G),G);
      t128 = d_dWdx_dW;//diff(cx(c),c);
      t129 = t30*t128;
      t131 = t10*t10;
      t132 = 1/t131;
      t135 = d_dWdy_dW;//diff(cy(c),c);
      t137 = d_dWdz_dW;//diff(cz(c),c);
      t139 = t4*t128+t6*t135+t8*t137;
      t140 = 2.0*t5*t132*t139;
      t142 = t41*t135;
      t143 = 2.0/3.0*t142;
      t145 = 2.0*t7*t132*t139;
      t146 = t145/3.0;
      t148 = t8*t11*t137;
      t149 = 2.0/3.0*t148;
      t151 = 2.0*t9*t132*t139;
      t152 = t151/3.0;
      t157 = 2.0/3.0*t129;
      t158 = t140/3.0;
      t165 = t128*t11;
      t167 = t4*t132;
      t184 = (-4.0/3.0*t129+2.0/3.0*t140+t143-t146+t149-t152)*t19+(-4.0/3.0*
t142+2.0/3.0*t145+t157-t158+t149-t152)*t24+(-4.0/3.0*t148+2.0/3.0*t151+t157-
t158+t143-t146)*t28-t165*t34+2.0*t167*t34*t139-t30*t135*t33-t165*t39+2.0*t167*
t39*t139-t30*t137*t38-t135*t11*t45+2.0*t6*t132*t45*t139-t41*t137*t44;
      t188 = dWdx.rho;//rhox(rho);
      t190 = rho*t54+Fsd*t188;
      t193 = dWdy.rho;//rhoy(rho);
      t195 = rho*t60+Fsd*t193;
      t198 = dWdz.rho;//rhoz(rho);
      t200 = rho*t66+Fsd*t198;
      t207 = 1/t51/t10;
      t246 = rho*rho;
      t248 = t89*t89;
      t254 = d_dWdx_dW;//diff(Fsdx(Fsd),Fsd);
      t258 = d_dWdy_dW;//diff(Fsdy(Fsd),Fsd);
      t262 = d_dWdz_dW;//diff(Fsdz(Fsd),Fsd);
      if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K){
      t279 = 1/t81;
      t291 = pow(k_fsd,0.15);
      t295 = Cv*t81;
      t296 = filter*t95;
      t301 = 0.2*t24;
      t302 = 0.2*t28;
      t304 = filter*(0.4*t19-t301-t302);
      t308 = 0.6666666667*rho*k_fsd;
      t317 = filter*t33;
      t321 = filter*t38;
      t326 = filter*t106;
      t334 = 0.2E1*t19;
      t336 = filter*(0.4*t24-t334-t302);
      t347 = filter*t44;
      t352 = filter*t117;
      t364 = filter*(0.4*t28-t334-t301);
      t376 = Cv*t279;
      t379 = 0.6666666667*rho;
      t390 = t33*t33;
      t394 = t38*t38;
      t398 = t44*t44;
      t403 = pow(k_fsd,0.5);
      }
      dStdW(5,0) = t1*Fsd;
      dStdW(5,6) = t1*rho;
      dStdW(6,0) = t47*Fsd-Laminar_Flame_Speed*(t50*(-t53*(t54+Fsd*t55)-t59*(t60+Fsd
*t61)-t65*(t66+Fsd*t67))+t73*t77)+t82*Fsd*t83-2.0*t86*t87*rho*t90;
      dStdW(6,1) = (t18*t95-t30*t6*t97-t30*t8*t100)*Fsd*rho;
      dStdW(6,2) = (t23*t106-t30*t6*t108-t41*t8*t111)*Fsd*rho;
      dStdW(6,3) = (t27*t117-t30*t8*t119-t41*t8*t122)*Fsd*rho;
      dStdW(6,5) = t184*Fsd*rho-Laminar_Flame_Speed*(tau_fsd*(-t53*t190-t59*t195-t65
*t200)+t50*(-t128*t52*t190+t4*t207*t190*t139-t135*t52*t195+t6*t207*t195*t139-
t137*t52*t200+t8*t207*t200*t139)+t73*rho*(-2.0*t53*t128+t5*t207*t139-2.0*t59*
t135+t7*t207*t139-2.0*t65*t137+t9*t207*t139))-t86*t87*t246/t248;
      dStdW(6,6) = t47*rho-Laminar_Flame_Speed*(t50*(-t53*(rho*t254+t188)-t59*(rho*
t258+t193)-t65*(rho*t262+t198))+tau_fsd*rho*t77)+t82*rho*t83-2.0*t86*Fsd*t246*
t90;
      if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K){
      dStdW(6,7) = kappa_fsd*t279*Fsd*rho*t83/2.0;
      dStdW(7,0) = -0.6666666667*k_fsd*t19-0.6666666667*k_fsd*t24-0.6666666667*k_fsd*t28-
Cs*t291*t83;
      dStdW(7,1) = 0.2666666667E1*t295*t296*t19+(0.6666666667*t295*t304-t308)*
t95-0.1333333333E1*t295*t296*t24-0.1333333333E1*t295*t296*t28+0.4E1*t295*t317*
t97+0.4E1*t295*t321*t100;
      dStdW(7,2) = -0.1333333333E1*t295*t326*t19+0.2666666667E1*t295*t326*t24+
(0.6666666667*t295*t336-t308)*t106-0.1333333333E1*t295*t326*t28+0.4E1*t295*t317
*t108+0.4E1*t295*t347*t111;
      dStdW(7,3) = -0.1333333333E1*t295*t352*t19-0.1333333333E1*t295*t352*t24+
0.2666666667E1*t295*t352*t28+(0.6666666667*t295*t364-t308)*t117+0.4E1*t295*t321
*t119+0.4E1*t295*t347*t122;
      dStdW(7,7) = (0.3333333334*t376*t304-t379)*t19+(0.3333333334*t376*t336-
t379)*t24+(0.3333333334*t376*t364-t379)*t28+0.1E1*t376*filter*t390+0.1E1*t376*
filter*t394+0.1E1*t376*filter*t398-0.15E1*Cs*rho*t403*t83;
      }
}

double LES3DFsd_cState::T(void) const{

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
         cout<<"\nTemperature didn't converge in LES3DFsd_cState::T(void)";
         cout<<" with polytopic Tguess "<<Tguess<<", or lower than Tmin "<<low_temp_range<<" using "<<T;
      }   
   }
   return T;
} 

double LES3DFsd_pState::E(void) const {   
   return ( rho*(e()+HALF*v.sqr()+k) );   
}

double LES3DFsd_pState::H(void) const{
   return ( rho*(h()+HALF*v.sqr()+FIVE*k/THREE) );   
}

double LES3DFsd_pState::Hs(void) const{
   return ( rho*(hs()+HALF*v.sqr()+FIVE*k/THREE) );   
}

/************* Turbulence modified pressure **************/
double LES3DFsd_pState::pmodified(void) const { 
   return ( p + 2.0*rho*k/3.0 );
}

/*********** Characteristic filter width ************/ 
double LES3DFsd_pState::filter_width(const double &Volume) const{
  return pow(Volume,1.0/3.0); 
}

double LES3DFsd_pState::abs_strain_rate(const LES3DFsd_pState &dWdx, 
                                        const LES3DFsd_pState &dWdy, 
                                        const LES3DFsd_pState &dWdz) const{
   Tensor3D strain_rate; 
   strain_rate.zero();
   strain_rate.xx = dWdx.v.x;
   strain_rate.yy = dWdy.v.y;
   strain_rate.zz = dWdz.v.z;
   strain_rate.xy = 0.5*(dWdy.v.x + dWdx.v.y);
   strain_rate.yz = 0.5*(dWdz.v.y + dWdy.v.z);
   strain_rate.xz = 0.5*(dWdx.v.z + dWdz.v.x);
  // S[i,j]*S[i,j]
  double SS = sqr(strain_rate.xx) + sqr(strain_rate.yy) + sqr(strain_rate.zz)
            + 2.0*(sqr(strain_rate.xy) + sqr(strain_rate.yz) + sqr(strain_rate.xz));
  // sqrt(2*S*S)
  return sqrt(2.0*SS);
}

double LES3DFsd_pState::eddy_viscosity(const LES3DFsd_pState &dWdx,
 				       const LES3DFsd_pState &dWdy,
				       const LES3DFsd_pState &dWdz,
                                       const int &Flow_Type, 
                                       const double &Volume) const{
    double filter = filter_width(Volume);
  if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY) {
    double Cs = 0.18;
    return(rho*sqr(Cs*filter)*abs_strain_rate(dWdx,dWdy,dWdz));
  }else if(Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K) {
    double Cv = 0.086;
    return(rho*Cv*sqrt(k)*filter);
  }
}

/******************************************************
  Calculating the thermal diffusion component of 
  the heat flux vector (qflux)
   sum( hs * Ds * grad cs)
 *******************************************************/
Vector3D LES3DFsd_pState::thermal_diffusion_t(const LES3DFsd_pState &dWdx,
   	                                      const LES3DFsd_pState &dWdy,
				              const LES3DFsd_pState &dWdz,
                                              const int &Flow_Type, 
                                              const double &Volume) const{
   Vector3D sum;
   sum.zero();
   double Temp = T();
   double Dm_t = eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,Volume)/(rho*Sc_t());
   for(int i=0; i<ns; i++){
     sum +=(/*specdata[i].Enthalpy(Temp) +*/ specdata[i].Heatofform())*Dm_t*spec[i].gradc;
   }
   return sum;
}


double LES3DFsd_pState::Pr_t(void) const{
   return (0.9); 
}

double LES3DFsd_pState::Sc_t(void) const{
   return (1.0);
}

double LES3DFsd_cState::p() const{
   return (rho*Rtot()*T());
}

double LES3DFsd_cState::k() const{
   return (rhok/rho);
}

double LES3DFsd_cState::C() const{
   return (rhoC/rho);
}

double LES3DFsd_cState::Fsd() const{
   return (rhoFsd/rho);
}

/********************************************************
 * Routine: RoeAverage (Roe Averages)                   *
 *                                                      *
 * This function returns the Roe-averaged primitive     *
 * solution state given left and right primitive        *
 * solution variables.  See Roe (1981).                 *
 *                                                      *
 ********************************************************/
LES3DFsd_pState LES3DFsd_pState::RoeAverage(const LES3DFsd_pState &Wl,
                                            const LES3DFsd_pState &Wr) {

   double Hl, Hr, srhol, srhor;
   double Ha, ha;
   LES3DFsd_pState Temp;

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
   Temp.C = (srhol*Wl.C+srhor*Wr.C)/(srhol+srhor);
   Temp.Fsd = (srhol*Wl.Fsd+srhor*Wr.Fsd)/(srhol+srhor);
   Temp.k = (srhol*Wl.k+srhor*Wr.k)/(srhol+srhor);

   for(int i=0; i<Wl.ns; i++){
      Temp.spec[i].c = (srhol*Wl.spec[i].c + srhor*Wr.spec[i].c)/(srhol+srhor);
   }

   //Temp.p = (srhol*Wl.p + srhor*Wr.p)/(srhol+srhor);
   Ha = (srhol*Hl+srhor*Hr)/(srhol+srhor);
   ha = Ha - HALF*(sqr(Temp.v.x)+sqr(Temp.v.y)+sqr(Temp.v.z)) - FIVE*Temp.k/THREE;
   
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
LES3DFsd_cState  LES3DFsd_pState::FluxHLLE_x(const LES3DFsd_pState &Wl,
                                             const LES3DFsd_pState &Wr) {
   
   double wavespeed_l, wavespeed_r;
 
    /* solnvec in  Wa (lambdas_l, lambdas_r, lambdas_a) 
       is allocated using new  */ 
   LES3DFsd_pState Wa, lambdas_l, lambdas_r, lambdas_a;
   LES3DFsd_cState Flux, dUrl;
   int NUM_VAR_3D = Wl.num_vars;
   
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
   wavespeed_r = max(lambdas_r[NUM_VAR_3D-NUM_LES3D_VAR_EXTRA-lambdas_r.ns],
                     lambdas_a[NUM_VAR_3D-NUM_LES3D_VAR_EXTRA-lambdas_a.ns]);
   
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

LES3DFsd_cState  LES3DFsd_pState::FluxHLLE_x(const LES3DFsd_cState &Ul,
                                             const LES3DFsd_cState &Ur) {
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
LES3DFsd_cState  LES3DFsd_pState::FluxHLLE_n(const LES3DFsd_pState &Wl,
                                             const LES3DFsd_pState &Wr,
                                             const Vector3D &norm_dir) {

   double Wl_ur_norm, Wl_ur_tang;
   double Wr_ur_norm, Wr_ur_tang ;
   double Wr_ur_tang_z;
   
   Vector3D Flux_rotated_x, Flux_rotated_tang_y, Flux_rotated_tang_z ;
   Vector3D Wl_ur_tang_vector, Wr_ur_tang_vector;
   Vector3D Wl_ur_tang_unit_vector, Wr_ur_tang_unit_vector;
   Vector3D Wr_ur_tang_z_vector, Wr_ur_tang_z_unit_vector;
   
   //solnvec in  Wl_rotated (Wr_rotated) is allocated using new 
   LES3DFsd_pState Wl_rotated, Wr_rotated;
   LES3DFsd_cState Flux, Flux_rotated;
   
   
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

LES3DFsd_cState  LES3DFsd_pState::FluxHLLE_n(const LES3DFsd_cState &Ul,
                                             const LES3DFsd_cState &Ur,
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
LES3DFsd_pState HartenFixPos(const LES3DFsd_pState &lambdas_a,
                             const LES3DFsd_pState &lambdas_l,
                             const LES3DFsd_pState &lambdas_r) {
   
   LES3DFsd_pState NEW;
   NEW.rho = HartenFixPos(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
   NEW.v.x = HALF*(lambdas_a[2]+fabs(lambdas_a[2]));
   NEW.v.y = HALF*(lambdas_a[3]+fabs(lambdas_a[3]));
   NEW.v.z = HALF*(lambdas_a[4]+fabs(lambdas_a[4]));
   NEW.p = HartenFixPos(lambdas_a[5],lambdas_l[5],lambdas_r[5]);
   NEW.C = HALF*(lambdas_a[6] + fabs(lambdas_a[6]));
   NEW.Fsd = HALF*(lambdas_a[7] +fabs(lambdas_a[7]));
   NEW.k = HALF*(lambdas_a[8] + fabs(lambdas_a[8]));

   for( int i=(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3D_VAR_EXTRA +1); i<=NEW.num_vars; i++){
      NEW.spec[i-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3D_VAR_EXTRA+1)].c = 
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
LES3DFsd_pState HartenFixNeg(const LES3DFsd_pState &lambdas_a,
                             const LES3DFsd_pState &lambdas_l,
                             const LES3DFsd_pState &lambdas_r) {

   LES3DFsd_pState NEW;
   NEW.rho = HartenFixNeg(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
   NEW.v.x = HALF*(lambdas_a[2]-fabs(lambdas_a[2]));
   NEW.v.y = HALF*(lambdas_a[3]-fabs(lambdas_a[3]));
   NEW.v.z = HALF*(lambdas_a[4]-fabs(lambdas_a[4]));
   NEW.p = HartenFixNeg(lambdas_a[5],lambdas_l[5],lambdas_r[5]);
   NEW.C = HALF*(lambdas_a[6] -fabs(lambdas_a[6]));
   NEW.Fsd = HALF*(lambdas_a[7] - fabs(lambdas_a[7]));
   NEW.k = HALF*(lambdas_a[8] -fabs(lambdas_a[8]));
   
   for( int i = (NUM_EULER3D_VAR_SANS_SPECIES  + NUM_LES3D_VAR_EXTRA+1); i<=NEW.num_vars; i++){
      NEW.spec[i-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3D_VAR_EXTRA+1)].c = 
         HALF*(lambdas_a[i]-fabs(lambdas_a[i]));
   }

   return (NEW);
}


// Flux Roe -- based on Harten fix 
LES3DFsd_cState  LES3DFsd_pState::FluxRoe_x(const  LES3DFsd_pState &Wl,  
                                            const  LES3DFsd_pState &Wr){
   
   LES3DFsd_pState Wa, dWrl, wavespeeds, 
      lambdas_l, lambdas_r, lambdas_a;
   LES3DFsd_cState Flux;

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
      for (int i=1 ; i < Wl.num_vars; i++) {
         if (wavespeeds[i] < ZERO) {
            Flux += wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i);
          }
      }
    } else {
      Flux = Wr.F();
      wavespeeds = HartenFixPos(lambdas_a,
                                lambdas_l,
                                lambdas_r);
      for (int i=1; i < Wl.num_vars; i++) {
         if (wavespeeds[i] > ZERO) {
            Flux -= wavespeeds[i]*(Wa.lp_x(i)*dWrl)*Wa.rc_x(i); 
  //          cout<<"\n Wa.lp_x(i) = "<<Wa.lp_x(i)<<" wavespeeds[i]= "<<wavespeeds[i]<<endl;
         }
      }
   } 
    
   /* Return solution flux. */    
   return (Flux);    
    
   
}
   
LES3DFsd_cState  LES3DFsd_pState::FluxRoe_n(const LES3DFsd_pState &Wl,
                                            const LES3DFsd_pState &Wr,
                                            const Vector3D &norm_dir){
   
   double Wl_ur_norm, Wl_ur_tang;
   double Wr_ur_norm, Wr_ur_tang ;
   double Wr_ur_tang_z;
   
   Vector3D Flux_rotated_x, Flux_rotated_tang_y, Flux_rotated_tang_z ;
   Vector3D Wl_ur_tang_vector, Wr_ur_tang_vector;
   Vector3D Wl_ur_tang_unit_vector, Wr_ur_tang_unit_vector;
   Vector3D Wr_ur_tang_z_vector, Wr_ur_tang_z_unit_vector;
   
   //solnvec in  Wl_rotated (Wr_rotated) is allocated using new 
   LES3DFsd_pState Wl_rotated, Wr_rotated;
   LES3DFsd_cState Flux, Flux_rotated;
   
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
LES3DFsd_cState LES3DFsd_pState::FluxAUSMplus_up_x(const LES3DFsd_pState &Wl,
                                                   const LES3DFsd_pState &Wr) {
 
  LES3DFsd_cState Flux;
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
 
  ahalf = HALF*(Wl.a() + Wr.a());
  rhohalf = HALF*(Wl.rho + Wr.rho); 
   
 
  // Determine the left and right state Mach numbers based on the
  // intermediate state sound speed:
  Ml = Wl.v.x/ahalf;
  Mr = Wr.v.x/ahalf;
 
  // Determine the reference Mach number, scaling function and coefficient
  double M2_bar, M2_ref, fa;
  M2_bar = (Wl.v.x*Wl.v.x + Wr.v.x*Wr.v.x)/(TWO*ahalf*ahalf);
  M2_ref = min(ONE, max(M2_bar, Wl.Mref*Wl.Mref));
  if (M2_ref > ONE || M2_ref < 0.0) cout << "\nM2_ref out of range in AUSM+-up.";
  //fa = sqrt(M2_ref)*(TWO - sqrt(M2_ref));
  fa = sqrt(sqr(ONE - M2_ref)*M2_bar + FOUR*M2_ref)/(ONE + M2_ref);
  if (fa > ONE || fa <= ZERO) cout << "\nfa out of range in AUSM+-up.";
  alpha = (3.0/16.0)*(-4.0 + 5.0*fa*fa);
  if (alpha < (-3.0/4.0)  ||  alpha > (3.0/16.0)) cout << "\nalpha out of range in AUSM+-up.";
 
 
  // Determine the left state split Mach number:
  if (fabs(Ml) >= ONE) {
    Mplus = 0.5*(Ml+fabs(Ml));
    pplus = 0.5*(Ml+fabs(Ml))/Ml;
  } else {
    Mplus = 0.25*sqr(Ml+1.0) * (1.0 - 16.0*beta*(-0.25*sqr(Ml-1.0)));
    pplus = 0.25*sqr(Ml+1.0) * ((2.0 - Ml) - 16.0*alpha*Ml*(-0.25*sqr(Ml-1.0)));
  }
  
  // Determine the right state split Mach number:
  if (fabs(Mr) >= ONE) {
    Mminus = 0.5*(Mr-fabs(Mr));
    pminus = 0.5*(Mr-fabs(Mr))/Mr;        
  } else {
    Mminus = -0.25*sqr(Mr-1.0) * (1.0 + 16.0*beta*0.25*sqr(Mr+1.0));
    pminus = -0.25*sqr(Mr-1.0) * ((-2.0 - Mr) + 16.0*alpha*Mr*0.25*sqr(Mr+1.0));
  } 

  // Determine the left state split Mach number:
//   if (fabs(Ml) >= ONE) {
//     Mplus = Mplus_1(Ml);
//     pplus = Mplus_1(Ml)/Ml;
//   } else {
//     Mplus = Mplus_2(Ml) * (1.0 - 16.0*beta*Mminus_2(Ml));
//     pplus = Mplus_2(Ml) * ((2.0 - Ml) - 16.0*alpha*Ml*Mminus_2(Ml));
//   }
  
  // Determine the right state split Mach number:
//   if (fabs(Mr) >= ONE) {
//     Mminus = Mminus_1(Mr);
//     pminus = Mminus_1(Mr)/Mr;        
//   } else {
//     Mminus = Mminus_2(Mr) * (1.0 + 16.0*beta*Mplus_2(Mr));
//     pminus = Mminus_2(Mr) * ((-2.0 - Mr) + 16.0*alpha*Mr*Mplus_2(Mr));
//   } 
 
   // Determine the intermediate state Mach number, pressure and mass flux:
  Mhalf = Mplus + Mminus
    - (Kp/fa)*max((ONE - sigma*M2_bar), ZERO)*(Wr.p - Wl.p)/(rhohalf*ahalf*ahalf);
 
  phalf = pplus*Wl.p + pminus*Wr.p
    - Ku*pplus*pminus*TWO*rhohalf*(fa*ahalf)*(Wr.v.x - Wl.v.x);
 
  mass_flux_half = (Mhalf > ZERO) ? ahalf*Mhalf*Wl.rho : ahalf*Mhalf*Wr.rho; 
  
  // Determine the intermediate state convective solution flux:
  if (mass_flux_half > ZERO) {
    Flux.rho = ONE;
    Flux.rhov.x = Wl.v.x; 
    Flux.rhov.y = Wl.v.y; 
    Flux.E = Wl.H()/Wl.rho;
    Flux.rhoC = Wl.C;
    Flux.rhoFsd = Wl.Fsd;
    Flux.rhok = Wl.k;

//     if(Wl.nscal > 0){
//       for(int i=0; i<Wl.nscal; ++i){
// 	Flux.rhoscalar[i] = Wl.scalar[i];
//       }
//     }
//     for(int i=0; i<Wl.ns; ++i){
//       Flux.rhospec[i].c = Wl.spec[i].c;
//     }
  } else {
    Flux.rho = ONE;
    Flux.rhov.x = Wr.v.x; 
    Flux.rhov.y = Wr.v.y; 
    Flux.E = Wr.H()/Wr.rho;
    Flux.rhoC = Wr.C;
    Flux.rhoFsd = Wr.Fsd;
    Flux.rhok = Wr.k;

//     if(Wr.nscal > 0){
//       for(int i=0; i<Wr.nscal; ++i){
// 	Flux.rhoscalar[i] = Wr.scalar[i];
//       }
//     }
//     for(int i=0; i<Wr.ns; ++i){
//       Flux.rhospec[i].c = Wr.spec[i].c;
//     }
  } //end if

  Flux = mass_flux_half*Flux;
 
  // Add the pressure contribution to the intermediate state solution flux:
  Flux[2] += phalf;
  
  // Return solution flux.
  return Flux;

}

// LES3DFsd_cState FluxAUSMplus_up_x(const LES3DFsd_cState &Ul,
//                                   const LES3DFsd_cState &Ur) {
//   return FluxAUSMplus_up_x(Ul.W(),Ur.W());
// }
 
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
LES3DFsd_cState LES3DFsd_pState::FluxAUSMplus_up_n(const LES3DFsd_pState &Wl,
                                                   const LES3DFsd_pState &Wr,
                                                   const Vector3D &norm_dir) {
 

   double Wl_ur_norm, Wl_ur_tang;
   double Wr_ur_norm, Wr_ur_tang ;
   double Wr_ur_tang_z;
   
   Vector3D Flux_rotated_x, Flux_rotated_tang_y, Flux_rotated_tang_z ;
   Vector3D Wl_ur_tang_vector, Wr_ur_tang_vector;
   Vector3D Wl_ur_tang_unit_vector, Wr_ur_tang_unit_vector;
   Vector3D Wr_ur_tang_z_vector, Wr_ur_tang_z_unit_vector;
   
   //solnvec in  Wl_rotated (Wr_rotated) is allocated using new 
   LES3DFsd_pState Wl_rotated, Wr_rotated;
   LES3DFsd_cState Flux, Flux_rotated;
   
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

//   double cos_angle, sin_angle;
//   LES3DFsd_pState Wl_rotated, Wr_rotated;
//   LES3DFsd_cState Flux, Flux_rotated;
 
//   /* Determine the direction cosine's for the frame
//      rotation. */
 
//   cos_angle = norm_dir.x; 
//   sin_angle = norm_dir.y;

//   /* Apply the frame rotation and evaluate left and right
//      solution states in the local rotated frame defined
//      by the unit normal vector. */
//   Wl_rotated.Copy(Wl);
//   Wr_rotated.Copy(Wr);
   
//   Wl_rotated.v.x = Wl.v.x*cos_angle + Wl.v.y*sin_angle;
//   Wl_rotated.v.y = - Wl.v.x*sin_angle + Wl.v.y*cos_angle;
 
//   Wr_rotated.v.x = Wr.v.x*cos_angle + Wr.v.y*sin_angle;
//   Wr_rotated.v.y = - Wr.v.x*sin_angle + Wr.v.y*cos_angle;
 
   
//   // Evaluate the intermediate state solution flux in the rotated frame.
//   Flux_rotated = FluxAUSMplus_up(Wl_rotated,Wr_rotated);
  
 
//   // Rotate back to the original Cartesian reference frame and return
//   // the solution flux.
//   Flux.Copy(Flux_rotated);
 
//   Flux.rhov.x = Flux_rotated.rhov.x*cos_angle - Flux_rotated.rhov.y*sin_angle;
//   Flux.rhov.y = Flux_rotated.rhov.x*sin_angle + Flux_rotated.rhov.y*cos_angle;
 
//   Flux.zero_non_sol();

  // Return the solution flux.  
  return Flux;
 
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
LES3DFsd_pState LES3DFsd_pState::Reflect(const LES3DFsd_pState &W,
                                         const Vector3D &norm_dir) {
   Vector3D ur_norm, ur_tang, vr_tot;
   
   LES3DFsd_pState Temp;
   Temp = W;
   
   ur_norm = dot(W.v, norm_dir)*norm_dir;
   ur_tang = W.v - ur_norm;
   
   ur_norm = -ur_norm;
   vr_tot = ur_norm + ur_tang;
   
   Temp.v = vr_tot;
   
   return (Temp);
}

LES3DFsd_pState LES3DFsd_pState::MovingWall(const LES3DFsd_pState &Win,
                                            const LES3DFsd_pState &Wout,
                                            const Vector3D &norm_dir, 
                                            const Vector3D &wall_velocity,
                                            const Vector3D &pressure_gradient,
                                            const int &TEMPERATURE_BC_FLAG) {
   LES3DFsd_pState Temp;
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

LES3DFsd_pState  LES3DFsd_pState::NoSlip(const LES3DFsd_pState &Win,
                                         const LES3DFsd_pState &Wout,
                                         const Vector3D &norm_dir,
                                         const Vector3D &pressure_gradient,
                                         const int &TEMPERATURE_BC_FLAG) {
   return (MovingWall(Win, Wout, norm_dir, Vector3D_ZERO, pressure_gradient, TEMPERATURE_BC_FLAG));
}

// molecular stress tensor
Tensor3D LES3DFsd_pState::lambda(const LES3DFsd_pState &dWdx, 
                                 const LES3DFsd_pState &dWdy,
                                 const LES3DFsd_pState &dWdz,
                                 const int &Flow_Type,
                                 const double &Volume){

   
   Tensor3D subfilter_stress;
   double mu_t = eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,filter_width(Volume));
   
   subfilter_stress.xx = 1.0/3.0*mu_t*(4.0*dWdx.v.x - 2.0*dWdy.v.y - 2.0*dWdz.v.z) - 2.0/3.0*rho*k;
   subfilter_stress.xy = mu_t*(dWdx.v.y + dWdy.v.x);
   subfilter_stress.xz = mu_t*(dWdx.v.z + dWdz.v.x);
   subfilter_stress.yy = 1.0/3.0*mu_t*(4.0*dWdy.v.y - 2.0*dWdx.v.x - 2.0*dWdz.v.z) - 2.0/3.0*rho*k;
   subfilter_stress.yz = mu_t*(dWdy.v.z + dWdz.v.y);
   subfilter_stress.zz = 1.0/3.0*mu_t*(4.0*dWdz.v.z - 2.0*dWdx.v.x - 2.0*dWdy.v.y) - 2.0/3.0*rho*k;
   
   return (subfilter_stress);
 
}

Tensor3D LES3DFsd_pState::lambda(const LES3DFsd_pState &dWdx, 
                                 const LES3DFsd_pState &dWdy,
                                 const LES3DFsd_pState &dWdz,
                                 const int &Flow_Type,
                                 const double &Volume) const {

   Tensor3D subfilter_stress;
   double mu_t = eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,filter_width(Volume));
   
   subfilter_stress.xx = 1.0/3.0*mu_t*(4.0*dWdx.v.x - 2.0*dWdy.v.y - 2.0*dWdz.v.z) - 2.0/3.0*rho*k;
   subfilter_stress.xy = mu_t*(dWdx.v.y + dWdy.v.x);
   subfilter_stress.xz = mu_t*(dWdx.v.z + dWdz.v.x); 
   subfilter_stress.yy = 1.0/3.0*mu_t*(4.0*dWdy.v.y - 2.0*dWdx.v.x - 2.0*dWdz.v.z) - 2.0/3.0*rho*k;
   subfilter_stress.yz = mu_t*(dWdy.v.z + dWdz.v.y);
   subfilter_stress.zz = 1.0/3.0*mu_t*(4.0*dWdz.v.z - 2.0*dWdx.v.x - 2.0*dWdy.v.y) - 2.0/3.0*rho*k;
 
   return (subfilter_stress);
   
}

// heat flux: Fourier's law of heat conduction

 Vector3D LES3DFsd_pState::qflux_t(const LES3DFsd_pState &dWdx, 
                                   const LES3DFsd_pState &dWdy,
                                   const LES3DFsd_pState &dWdz,
                                   const int &Flow_Type,
                                   const double &Volume){
    double Rmix = Rtot();
    double kappa_t=eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,filter_width(Volume))*Cp()/Pr_t();
    Vector3D heat_flux, gradient_T;
  /* Temperature gradients from using the chain rule with the ideal gas law (P=rho*R*T) 
     dT/dx = 1/(rho*R) *( dP/dx - P/rho * drho/dx) */
    gradient_T.x = (1.0/(rho*Rmix)) * (dWdx.p -(p/rho)*dWdx.rho);
    gradient_T.y = (1.0/(rho*Rmix)) * (dWdy.p -(p/rho)*dWdy.rho);
    gradient_T.z = (1.0/(rho*Rmix)) * (dWdz.p -(p/rho)*dWdz.rho);
    
    heat_flux = -kappa_t*gradient_T;
    return (heat_flux);
 }

 Vector3D LES3DFsd_pState::qflux_t(const LES3DFsd_pState &dWdx, 
                                   const LES3DFsd_pState &dWdy,
                                   const LES3DFsd_pState &dWdz,
                                   const int &Flow_Type,
                                   const double &Volume) const{
  
    double Rmix = Rtot();
    double kappa_t=eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,filter_width(Volume))*Cp()/Pr_t();
    Vector3D heat_flux, gradient_T;
  /* Temperature gradients from using the chain rule with the ideal gas law (P=rho*R*T) 
     dT/dx = 1/(rho*R) *( dP/dx - P/rho * drho/dx) */
    gradient_T.x = (1.0/(rho*Rmix)) * (dWdx.p -(p/rho)*dWdx.rho);
    gradient_T.y = (1.0/(rho*Rmix)) * (dWdy.p -(p/rho)*dWdy.rho);
    gradient_T.z = (1.0/(rho*Rmix)) * (dWdz.p -(p/rho)*dWdz.rho);

    heat_flux = -kappa_t*gradient_T;

    return (heat_flux);
 }

// viscous flux
LES3DFsd_cState  LES3DFsd_pState::FluxViscous_n(const LES3DFsd_pState &Wl,
                                                const LES3DFsd_pState &Wr,
                                                const LES3DFsd_pState &Wc,
                                                const LES3DFsd_pState &Wc_Neigbor,
                                                const LES3DFsd_pState &dWdx,
                                                const LES3DFsd_pState &dWdy,
                                                const LES3DFsd_pState &dWdz,
                                                const LES3DFsd_pState &dWdx_Neigbor,
                                                const LES3DFsd_pState &dWdy_Neigbor,
                                                const LES3DFsd_pState &dWdz_Neigbor,
                                                const Vector3D &norm, const Vector3D &ts, 
                                                const double &deltad, const double &Volume, 
                                                const double &Volume_Neigbor, const int &Flow_Type){
   
   // construct the gradients on the cell interface (surface) 
   // based on Hybrid Average Gradient-Diamond-Path Approach
   // Grad Phi = (Phi_c_neigbor - Phi_c)/ds *norm/(norm . ts)
   //            + [bar (Grad Phi) - bar (Grad Phi) . ts norm/(norm.ts)]
  
   // weighted factor based on volume
   double alpha = Volume/(Volume + Volume_Neigbor);
   
   LES3DFsd_pState dWdx_Weighted, dWdy_Weighted, dWdz_Weighted, 
                   dWdx_face, dWdy_face, dWdz_face, Grad_middle_term;

   LES3DFsd_pState W_face;
   
   dWdx_Weighted = alpha*dWdx + (1.0 - alpha)*dWdx_Neigbor;
   dWdy_Weighted = alpha*dWdy + (1.0 - alpha)*dWdy_Neigbor;
   dWdz_Weighted = alpha*dWdz + (1.0 - alpha)*dWdz_Neigbor;

   dWdx_face = dWdx_Weighted ;
   dWdy_face = dWdy_Weighted ;
   dWdz_face = dWdz_Weighted ; 
   
  // W_face = HALF*(Wl + Wr);
   W_face = HALF*(Wc + Wc_Neigbor);

   return (W_face.Fv(dWdx_face, dWdy_face, dWdz_face, Flow_Type, Volume)*norm.x +
	   W_face.Gv(dWdx_face, dWdy_face, dWdz_face, Flow_Type, Volume)*norm.y +
	   W_face.Hv(dWdx_face, dWdy_face, dWdz_face, Flow_Type, Volume)*norm.z);
   
}

double LES3DFsd_pState::HeatRelease_Parameter(void) const {
  double Adiabatic_Temperature = 2218.0;
  return (Adiabatic_Temperature/298.0-1.0);
}

double LES3DFsd_pState::SFS_Kinetic_Energy_Fsd(const LES3DFsd_pState &dWdx,
                                               const LES3DFsd_pState &dWdy,
                                               const LES3DFsd_pState &dWdz,
                                               const int &Flow_Type,
                                               const double &Volume) const {
  if ( Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY ) {
    double CI = 0.005;
    return (0.0);//CI*sqr(filter_width(Volume)*abs_strain_rate(dWdx,dWdy,dWdz)));
  } else if ( Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K ) {
    return (k);
  }
}

double LES3DFsd_pState::Efficiency_Function_Fsd(const LES3DFsd_pState &dWdx,
                                                const LES3DFsd_pState &dWdy,
                                                const LES3DFsd_pState &dWdz,
                                                const int &Flow_Type,
                                                const double &Volume) const {
  double filter, kappa_fsd, k_fsd;
  double Laminar_Flame_Speed = 0.3837;
  double Laminar_Flame_Thickness = 4.4E-04;
  k_fsd = SFS_Kinetic_Energy_Fsd(dWdx,dWdy,dWdz,Flow_Type,Volume);
  filter = filter_width(Volume);
  kappa_fsd = 0.75*exp(-1.2/pow(sqrt(k_fsd)/Laminar_Flame_Speed,0.3))*pow(filter/Laminar_Flame_Thickness,2.0/3.0);
  return(kappa_fsd);
}

double LES3DFsd_pState::Progvar_Species_Grad(void) const {
  double Temp, stoich, ratio, f_ub, eta_fsd;
  double Equivalence_Ratio = 1.0;
  Temp = p/(rho*Rtot());
  if ( React.reactset_flag == CH4_1STEP ||
       React.reactset_flag == C3H8_1STEP ) {
  if ( React.reactset_flag == CH4_1STEP ){
    stoich = 2.0*specdata[1].Mol_mass()/specdata[0].Mol_mass();
    ratio = specdata[2].Mol_mass()/(specdata[2].Mol_mass()+2.0*specdata[3].Mol_mass());
    f_ub = specdata[0].Mol_mass()/(specdata[0].Mol_mass()+2.0*specdata[1].Mol_mass()+7.52*specdata[4].Mol_mass());
  }else if ( React.reactset_flag == C3H8_1STEP ){
    stoich = 5.0*specdata[1].Mol_mass()/specdata[0].Mol_mass();
    ratio = 3.0*specdata[2].Mol_mass()/(3.0*specdata[2].Mol_mass()+4.0*specdata[3].Mol_mass());
    f_ub = specdata[0].Mol_mass()/(specdata[0].Mol_mass()+5.0*specdata[1].Mol_mass()+18.8*specdata[4].Mol_mass());
  }
    eta_fsd = (specdata[0].Enthalpy(Temp)+specdata[0].Heatofform()-Cp(Temp)*Temp*specdata[0].Rs()/Rtot())*(-f_ub)
             +(specdata[1].Enthalpy(Temp)+specdata[1].Heatofform()-Cp(Temp)*Temp*specdata[1].Rs()/Rtot())*(-stoich*f_ub/Equivalence_Ratio)
	     +(specdata[2].Enthalpy(Temp)+specdata[2].Heatofform()-Cp(Temp)*Temp*specdata[2].Rs()/Rtot())*((1.0+stoich/Equivalence_Ratio)*f_ub*ratio)
	     +(specdata[3].Enthalpy(Temp)+specdata[3].Heatofform()-Cp(Temp)*Temp*specdata[3].Rs()/Rtot())*((1.0+stoich/Equivalence_Ratio)*f_ub*(1.0-ratio));
  }else if ( React.reactset_flag == H2O2_1STEP ){
    stoich = specdata[1].Mol_mass()/2.0/specdata[0].Mol_mass();
    f_ub = 2.0*specdata[0].Mol_mass()/(2.0*specdata[0].Mol_mass()+specdata[1].Mol_mass()+3.76*specdata[3].Mol_mass());
    eta_fsd = (specdata[0].Enthalpy(Temp)+specdata[0].Heatofform()-Cp(Temp)*Temp*specdata[0].Rs()/Rtot())*(-f_ub)
             +(specdata[1].Enthalpy(Temp)+specdata[1].Heatofform()-Cp(Temp)*Temp*specdata[1].Rs()/Rtot())*(-stoich*f_ub/Equivalence_Ratio)
             +(specdata[2].Enthalpy(Temp)+specdata[2].Heatofform()-Cp(Temp)*Temp*specdata[2].Rs()/Rtot())*((1.0+stoich/Equivalence_Ratio)*f_ub);
  }
  return (eta_fsd);
}

double LES3DFsd_pState::Reaction_Rate_Fsd(const LES3DFsd_pState &dWdx,
                                          const LES3DFsd_pState &dWdy,
                                          const LES3DFsd_pState &dWdz) const {
     double tau_fsd = HeatRelease_Parameter();
     double Laminar_Flame_Speed = 0.3837;
     double Reactants_Density = 1.13;
     if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ){//&& dWdy.C != ZERO && dWdz.C != ZERO ) {
     //       cout<<"reac=    "<<Reactants_Density*Laminar_Flame_Speed*Fsd*rho<<endl; 
       return ( Reactants_Density*Laminar_Flame_Speed*Fsd*rho );//-tau_fsd*Laminar_Flame_Speed*(rho*(1-2*C)*(dWdx.C+dWdy.C+dWdz.C)+C*(1-C)*(dWdx.rho+dWdy.rho+dWdz.rho)) );
    }else{
     return ( 0.0 );
    }
}

double LES3DFsd_pState::M_x(const LES3DFsd_pState &dWdx,
                            const LES3DFsd_pState &dWdy,
                            const LES3DFsd_pState &dWdz) const {
    double Mx, tau_fsd;
    tau_fsd = HeatRelease_Parameter();
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO ) {//&& Fsd != ZERO ) {
      Mx = dWdx.C/sqrt(sqr(dWdx.C)+sqr(dWdy.C)+sqr(dWdz.C));
//  cout<<"Mx=  "<<Mx<<"  c= "<<C<<"   dWdx.c=  "<<dWdx.C<<"   dWdy.c=  "<<dWdy.C<<"   dWdz.c=   "<<dWdz.C<<endl;
//  Mx = -((1.0+tau_fsd)*(1-exp(-0.2*1.28))/sqr(1.0+tau_fsd*C)+exp(-0.2*1.28))*dWdx.C/Fsd/rho;
//  Mx = -(1.0+tau_fsd)*dWdx.C/sqr(1.0+tau_fsd*C)/Fsd/rho;
//  Mx = -dWdx.C/Fsd/rho;
//  if ( Mx < -1.0 ) { Mx = -1.0; }
    return Mx;
     }else {
     return (0.0);
    }
}

double LES3DFsd_pState::M_y(const LES3DFsd_pState &dWdx,
                            const LES3DFsd_pState &dWdy,
                            const LES3DFsd_pState &dWdz) const {
    double My, tau_fsd;
    tau_fsd = HeatRelease_Parameter();
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ){//&& dWdy.C != ZERO && dWdz.C != ZERO && Fsd != ZERO ) {
    My = -dWdy.C/sqrt(sqr(dWdx.C)+sqr(dWdy.C)+sqr(dWdz.C));
//  My = -((1.0+tau_fsd)*(1-exp(-0.2*1.28))/sqr(1.0+tau_fsd*C)+exp(-0.2*1.28))*dWdy.C/Fsd/rho;
//  My = -(1.0+tau_fsd)*dWdy.C/sqr(1.0+tau_fsd*C)/Fsd/rho;
//  My = -dWdy.C/Fsd/rho;
//  if ( My < -1.0 ) { My = -1.0; }
    return My;
    }else {
    return (0.0);
   }
}

double LES3DFsd_pState::M_z(const LES3DFsd_pState &dWdx,
                            const LES3DFsd_pState &dWdy,
                            const LES3DFsd_pState &dWdz) const {

    double Mz, tau_fsd;
    tau_fsd = HeatRelease_Parameter();
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ){//&& dWdy.C != ZERO && dWdz.C != ZERO && Fsd != ZERO ) {
    Mz = -dWdz.C/sqrt(sqr(dWdx.C)+sqr(dWdy.C)+sqr(dWdz.C));
//  Mx = -((1.0+tau_fsd)*(1-exp(-0.2*1.28))/sqr(1.0+tau_fsd*C)+exp(-0.2*1.28))*dWdx.C/Fsd/rho;
//  Mx = -(1.0+tau_fsd)*dWdx.C/sqr(1.0+tau_fsd*C)/Fsd/rho;
//  Mx = -dWdx.C/Fsd/rho;
//  if ( Mx < -1.0 ) { Mx = -1.0; }
    return Mz;
    }else {
    return (0.0);
   }
}

double LES3DFsd_pState::Resolved_Strain(const LES3DFsd_pState &dWdx,
                                        const LES3DFsd_pState &dWdy,
                                        const LES3DFsd_pState &dWdz) const {
  double Mx, My, Mz, n_xx, n_yy, n_zz, n_xy, n_xz, n_yz, alpha_fsd;
  double resolved_strain_xx, resolved_strain_yy, resolved_strain_zz, resolved_strain_xy, resolved_strain_xz, resolved_strain_yz;

    Mx = M_x(dWdx,dWdy,dWdz);
    My = M_y(dWdx,dWdy,dWdz);
    Mz = M_z(dWdx,dWdy,dWdz);
    alpha_fsd = ONE - sqr(Mx) - sqr(My) - sqr(Mz);
    n_xx = sqr(Mx)+ONE/THREE*alpha_fsd;
    n_yy = sqr(My)+ONE/THREE*alpha_fsd;
    n_zz = sqr(Mz)+ONE/THREE*alpha_fsd;
    n_xy = Mx*My;
    n_xz = Mx*Mz;
    n_yz = My*Mz;

    if ( C < 0.999 && C > 0.001 &&  dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO ) {
    resolved_strain_xx = (ONE - n_xx)*dWdx.v.x*Fsd*rho;
    resolved_strain_yy = (ONE - n_yy)*dWdy.v.y*Fsd*rho;
    resolved_strain_zz = (ONE - n_zz)*dWdz.v.z*Fsd*rho;
    resolved_strain_xy = -n_xy*(dWdx.v.y + dWdy.v.x)*Fsd*rho;
    resolved_strain_xz = -n_xz*(dWdx.v.z + dWdz.v.x)*Fsd*rho;
    resolved_strain_yz = -n_yz*(dWdz.v.y + dWdy.v.y)*Fsd*rho;
    return (resolved_strain_xx + resolved_strain_yy + resolved_strain_zz + resolved_strain_xy + resolved_strain_xz + resolved_strain_yz);
    }else{
    return (0.0);
   }
}
double LES3DFsd_pState::Resolved_Propagation_Curvature(const LES3DFsd_pState &dWdx,
                                                       const LES3DFsd_pState &dWdy,
                                                       const LES3DFsd_pState &dWdz) const {
  double tau_fsd, Mx, My, Mz, resolved_propagation_curvature_x, resolved_propagation_curvature_y, resolved_propagation_curvature_z;
    tau_fsd = HeatRelease_Parameter();
    double Laminar_Flame_Speed = 0.3837;
    Mx = M_x(dWdx,dWdy,dWdz);
    My = M_y(dWdx,dWdy,dWdz);
    Mz = M_z(dWdx,dWdy,dWdz);
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO ) {
    resolved_propagation_curvature_x = -Laminar_Flame_Speed*(ONE+tau_fsd*C)*Mx*(rho*dWdx.Fsd+Fsd*dWdx.rho)-Laminar_Flame_Speed*tau_fsd*Fsd*rho*Mx*dWdx.C;
    resolved_propagation_curvature_y = -Laminar_Flame_Speed*(ONE+tau_fsd*C)*My*(rho*dWdy.Fsd+Fsd*dWdy.rho)-Laminar_Flame_Speed*tau_fsd*Fsd*rho*My*dWdy.C;
    resolved_propagation_curvature_z = -Laminar_Flame_Speed*(ONE+tau_fsd*C)*Mz*(rho*dWdz.Fsd+Fsd*dWdz.rho)-Laminar_Flame_Speed*tau_fsd*Fsd*rho*Mz*dWdz.C;
    return ( resolved_propagation_curvature_x + resolved_propagation_curvature_y + resolved_propagation_curvature_z );
    }else{
    return (0.0);
   }
}

double LES3DFsd_pState::SFS_Strain(const LES3DFsd_pState &dWdx,
                                   const LES3DFsd_pState &dWdy,
                                   const LES3DFsd_pState &dWdz,
                                   const int &Flow_Type,
                                   const double &Volume) const {
  double filter, kappa_fsd, k_fsd;
    filter = filter_width(Volume);
    k_fsd = SFS_Kinetic_Energy_Fsd(dWdx,dWdy,dWdz,Flow_Type,Volume); 
    kappa_fsd = Efficiency_Function_Fsd(dWdx,dWdy,dWdz,Flow_Type,Volume);   
   if ( C < 0.999 && C > 0.001 ) {
      return ( kappa_fsd*sqrt(k_fsd)*Fsd*rho/filter );
   }else{
    return (0.0);
   }
}

double LES3DFsd_pState::SFS_Curvature(const LES3DFsd_pState &dWdx,
                                      const LES3DFsd_pState &dWdy,
                                      const LES3DFsd_pState &dWdz) const {

    double Mx, My, Mz, alpha_fsd, beta_fsd;
    double Laminar_Flame_Speed = 0.3837;
    beta_fsd = 1.0;
    Mx = M_x(dWdx,dWdy,dWdz);
    My = M_y(dWdx,dWdy,dWdz);
    Mz = M_z(dWdx,dWdy,dWdz);
    alpha_fsd = ONE - sqr(Mx) - sqr(My) - sqr(Mz);
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO && Fsd != ZERO) {
    return ( -beta_fsd*Laminar_Flame_Speed*sqr(Fsd*rho)/(ONE-C) );

//     if ( Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_CHARLETTE ){
//     double tau_fsd = HeatRelease_Parameter();
//     double c_bar = (1.0+tau_fsd)*C/(1.0+tau_fsd*C);
//     return(-beta_fsd*Laminar_Flame_Speed*(Fsd-(1+tau_fsd)*sqrt(sqr(dWdx.C)+sqr(dWdy.C))/sqr(1+tau_fsd*C))*Fsd/c_bar/(1-c_bar));
   }else{
    return (0.0);
   }
}

double LES3DFsd_pState::M_xx(const LES3DFsd_pState &dWdx,
                             const LES3DFsd_pState &dWdy,
                             const LES3DFsd_pState &dWdz,
                             const LES3DFsd_pState &d_dWdx_dx,
                             const LES3DFsd_pState &d_dWdx_dy,
                             const LES3DFsd_pState &d_dWdx_dz) const {
    double Mxx, magnitude_C;
    magnitude_C = sqrt(sqr(dWdx.C)+sqr(dWdy.C)+sqr(dWdz.C));
    //    cout<<"d_dwdx_dx.c=   "<<d_dWdx_dx.C<<"  dwdx.c=   "<<dWdx.C<<endl;
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO ) { 
      Mxx = d_dWdx_dx.C/magnitude_C-dWdx.C*(dWdx.C*d_dWdx_dx.C+dWdy.C*d_dWdx_dy.C+dWdz.C*d_dWdx_dz.C)/pow(magnitude_C,3.0);
    return ( Mxx );
   }else{
    return (0.0);
   }
}

double LES3DFsd_pState::M_yy(const LES3DFsd_pState &dWdx,
                             const LES3DFsd_pState &dWdy,
                             const LES3DFsd_pState &dWdz,
                             const LES3DFsd_pState &d_dWdy_dy,
                             const LES3DFsd_pState &d_dWdx_dy,
                             const LES3DFsd_pState &d_dWdy_dz) const {
    double Myy, magnitude_C;
    magnitude_C = sqrt(sqr(dWdx.C)+sqr(dWdy.C)+sqr(dWdz.C));
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO ) { 
      Myy = -d_dWdy_dy.C/magnitude_C+dWdy.C*(dWdx.C*d_dWdx_dy.C+dWdy.C*d_dWdy_dy.C+dWdz.C*d_dWdy_dz.C)/pow(magnitude_C,3.0);
    return ( Myy );
   }else{
    return (0.0);
   }
}

double LES3DFsd_pState::M_zz(const LES3DFsd_pState &dWdx,
                             const LES3DFsd_pState &dWdy,
                             const LES3DFsd_pState &dWdz,
                             const LES3DFsd_pState &d_dWdz_dz,
                             const LES3DFsd_pState &d_dWdx_dz,
                             const LES3DFsd_pState &d_dWdy_dz) const {
    double Mzz, magnitude_C;
    magnitude_C = sqrt(sqr(dWdx.C)+sqr(dWdy.C)+sqr(dWdz.C));
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO ) { 
      Mzz = -d_dWdz_dz.C/magnitude_C+dWdz.C*(dWdx.C*d_dWdx_dz.C+dWdy.C*d_dWdy_dz.C+dWdz.C*d_dWdz_dz.C)/pow(magnitude_C,3.0);
    return ( Mzz );
   }else{
    return (0.0);
   }
}

double LES3DFsd_pState::M_xy(const LES3DFsd_pState &dWdx,
                             const LES3DFsd_pState &dWdy,
                             const LES3DFsd_pState &dWdz,
                             const LES3DFsd_pState &d_dWdy_dy,
                             const LES3DFsd_pState &d_dWdx_dy,
                             const LES3DFsd_pState &d_dWdy_dz) const {
    double Mxy, magnitude_C;
    magnitude_C = sqrt(sqr(dWdx.C)+sqr(dWdy.C)+sqr(dWdz.C));
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO ) { 
      Mxy = -d_dWdx_dy.C/magnitude_C+dWdx.C*(dWdx.C*d_dWdx_dy.C+dWdy.C*d_dWdy_dy.C+dWdz.C*d_dWdy_dz.C)/pow(magnitude_C,3.0);
    return ( Mxy );
   }else{
    return (0.0);
   }
}

double LES3DFsd_pState::M_xz(const LES3DFsd_pState &dWdx,
                             const LES3DFsd_pState &dWdy,
                             const LES3DFsd_pState &dWdz,
                             const LES3DFsd_pState &d_dWdz_dz,
                             const LES3DFsd_pState &d_dWdx_dz,
                             const LES3DFsd_pState &d_dWdy_dz) const {
    double Mxz, magnitude_C;
    magnitude_C = sqrt(sqr(dWdx.C)+sqr(dWdy.C)+sqr(dWdz.C));
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO ) { 
      Mxz = -d_dWdx_dz.C/magnitude_C+dWdx.C*(dWdx.C*d_dWdx_dz.C+dWdy.C*d_dWdy_dz.C+dWdz.C*d_dWdz_dz.C)/pow(magnitude_C,3.0);
    return ( Mxz );
   }else{
    return (0.0);
   }
}

double LES3DFsd_pState::M_yz(const LES3DFsd_pState &dWdx,
                             const LES3DFsd_pState &dWdy,
                             const LES3DFsd_pState &dWdz,
                             const LES3DFsd_pState &d_dWdz_dz,
                             const LES3DFsd_pState &d_dWdx_dz,
                             const LES3DFsd_pState &d_dWdy_dz) const {
    double Myz, magnitude_C;
    magnitude_C = sqrt(sqr(dWdx.C)+sqr(dWdy.C)+sqr(dWdz.C));
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO ) { 
      Myz = -d_dWdy_dz.C/magnitude_C+dWdy.C*(dWdx.C*d_dWdx_dz.C+dWdy.C*d_dWdy_dz.C+dWdz.C*d_dWdz_dz.C)/pow(magnitude_C,3.0);
    return ( Myz );
   }else{
    return (0.0);
   }
}

double LES3DFsd_pState::Resolved_Curvature(const LES3DFsd_pState &dWdx,
                                           const LES3DFsd_pState &dWdy,
                                           const LES3DFsd_pState &dWdz,
                                           const LES3DFsd_pState &d_dWdx_dx,
                                           const LES3DFsd_pState &d_dWdy_dy,
                                           const LES3DFsd_pState &d_dWdz_dz,
                                           const LES3DFsd_pState &d_dWdx_dy,
                                           const LES3DFsd_pState &d_dWdx_dz,
                                           const LES3DFsd_pState &d_dWdy_dz) const {
  double tau_fsd, Mxx, Myy, Mzz, resolved_curvature_xx, resolved_curvature_yy, resolved_curvature_zz;
   tau_fsd = HeatRelease_Parameter();
   double Laminar_Flame_Speed = 0.3837;
   Mxx = M_xx(dWdx,dWdy,dWdz,d_dWdx_dx,d_dWdx_dy,d_dWdx_dz);
   Myy = M_yy(dWdx,dWdy,dWdz,d_dWdy_dy,d_dWdx_dy,d_dWdy_dz);
   Mzz = M_zz(dWdx,dWdy,dWdz,d_dWdz_dz,d_dWdx_dz,d_dWdy_dz);

   if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C!= ZERO && dWdz.C != ZERO ) {
   resolved_curvature_xx = Laminar_Flame_Speed*(1.0+tau_fsd*C)*Fsd*rho*Mxx;
   resolved_curvature_yy = Laminar_Flame_Speed*(1.0+tau_fsd*C)*Fsd*rho*Myy;
   resolved_curvature_zz = Laminar_Flame_Speed*(1.0+tau_fsd*C)*Fsd*rho*Mzz;
   return ( resolved_curvature_xx + resolved_curvature_yy + resolved_curvature_zz);
  }else{
   return (0.0);
  }
}

double LES3DFsd_pState::Resolved_Propagation(const LES3DFsd_pState &dWdx,
                                             const LES3DFsd_pState &dWdy,
                                             const LES3DFsd_pState &dWdz,
                                             const LES3DFsd_pState &d_dWdx_dx,
                                             const LES3DFsd_pState &d_dWdy_dy,
                                             const LES3DFsd_pState &d_dWdz_dz,
                                             const LES3DFsd_pState &d_dWdx_dy,
                                             const LES3DFsd_pState &d_dWdx_dz,
                                             const LES3DFsd_pState &d_dWdy_dz) const {

  if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C!= ZERO && dWdz.C != ZERO ) {
     return ( Resolved_Propagation_Curvature(dWdx,dWdy,dWdz)
             -Resolved_Curvature(dWdx,dWdy,dWdz,d_dWdx_dx,d_dWdy_dy,d_dWdz_dz,d_dWdx_dy,d_dWdx_dz,d_dWdy_dz) );
  }else{
   return (0.0);
  }
}

double LES3DFsd_pState::Resolved_Convection_Progvar (const LES3DFsd_pState &dWdx,
                                                     const LES3DFsd_pState &dWdy,
                                                     const LES3DFsd_pState &dWdz) const {

    double resolved_convection_progvar_x, resolved_convection_progvar_y, resolved_convection_progvar_z;
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C!= ZERO && dWdz.C != ZERO ) {
    resolved_convection_progvar_x = -(dWdx.rho*v.x*C+rho*dWdx.v.x*C+rho*v.x*dWdx.C);
    resolved_convection_progvar_y = -(dWdy.rho*v.y*C+rho*dWdy.v.y*C+rho*v.y*dWdy.C);
    resolved_convection_progvar_z = -(dWdz.rho*v.z*C+rho*dWdz.v.z*C+rho*v.z*dWdz.C);
    return( resolved_convection_progvar_x+resolved_convection_progvar_y+resolved_convection_progvar_z );
   }else{
    return (0.0);
   }
}

double LES3DFsd_pState::Resolved_Convection_Fsd (const LES3DFsd_pState &dWdx,
                                                 const LES3DFsd_pState &dWdy,
                                                 const LES3DFsd_pState &dWdz) const {

    double resolved_convection_fsd_x, resolved_convection_fsd_y, resolved_convection_fsd_z;
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C!= ZERO && dWdz.C != ZERO ) {
    resolved_convection_fsd_x = -(dWdx.rho*v.x*Fsd+rho*dWdx.v.x*Fsd+rho*v.x*dWdx.Fsd);
    resolved_convection_fsd_y = -(dWdy.rho*v.y*Fsd+rho*dWdy.v.y*Fsd+rho*v.y*dWdy.Fsd);
    resolved_convection_fsd_z = -(dWdz.rho*v.z*Fsd+rho*dWdz.v.z*Fsd+rho*v.z*dWdz.Fsd);
    return( resolved_convection_fsd_x+resolved_convection_fsd_y+resolved_convection_fsd_z );
   }else{
    return (0.0);
   }
}

double LES3DFsd_pState::NGT_Progvar (const LES3DFsd_pState &dWdx,
                                     const LES3DFsd_pState &dWdy,
                                     const LES3DFsd_pState &dWdz) const {

    double tau_fsd, NGT_progvar_x, NGT_progvar_y, NGT_progvar_z;
    tau_fsd = HeatRelease_Parameter();
    double Laminar_Flame_Speed = 0.3837;

    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C!= ZERO && dWdz.C != ZERO ) {
    NGT_progvar_x = -tau_fsd*Laminar_Flame_Speed*(rho*(1-2*C)*dWdx.C+C*(1-C)*dWdx.rho);
    NGT_progvar_y = -tau_fsd*Laminar_Flame_Speed*(rho*(1-2*C)*dWdy.C+C*(1-C)*dWdy.rho);
    NGT_progvar_z = -tau_fsd*Laminar_Flame_Speed*(rho*(1-2*C)*dWdz.C+C*(1-C)*dWdz.rho);
    return ( NGT_progvar_x+NGT_progvar_y+NGT_progvar_z );
   }else{
    return (0.0);
   }
}

double LES3DFsd_pState::NGT_Fsd (const LES3DFsd_pState &dWdx,
                                 const LES3DFsd_pState &dWdy,
                                 const LES3DFsd_pState &dWdz,
                                 const LES3DFsd_pState &d_dWdx_dx,
                                 const LES3DFsd_pState &d_dWdy_dy,
                                 const LES3DFsd_pState &d_dWdz_dz,
                                 const LES3DFsd_pState &d_dWdx_dy,
                                 const LES3DFsd_pState &d_dWdx_dz,
                                 const LES3DFsd_pState &d_dWdy_dz) const {

    double tau_fsd, Mx, My, Mz, Mxx, Myy, Mzz, NGT_fsd_x, NGT_fsd_y, NGT_fsd_z;
    tau_fsd = HeatRelease_Parameter();
    double Laminar_Flame_Speed = 0.3837;
    Mx = M_x(dWdx,dWdy,dWdz);
    My = M_y(dWdx,dWdy,dWdz);
    Mz = M_z(dWdx,dWdy,dWdz);
    Mxx = M_xx(dWdx,dWdy,dWdz,d_dWdx_dx,d_dWdx_dy,d_dWdx_dz);
    Myy = M_yy(dWdx,dWdy,dWdz,d_dWdy_dy,d_dWdx_dy,d_dWdy_dz);
    Mzz = M_zz(dWdx,dWdy,dWdz,d_dWdz_dz,d_dWdx_dz,d_dWdy_dz);
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C!= ZERO && dWdz.C != ZERO ) {
    NGT_fsd_x = -tau_fsd*Laminar_Flame_Speed*((0.5-C)*(Fsd*Mx*dWdx.rho+rho*Mx*dWdx.Fsd+rho*Fsd*Mxx)-rho*Fsd*Mx*dWdx.C);
    NGT_fsd_y = -tau_fsd*Laminar_Flame_Speed*((0.5-C)*(Fsd*My*dWdy.rho+rho*My*dWdy.Fsd+rho*Fsd*Myy)-rho*Fsd*My*dWdy.C);
    NGT_fsd_z = -tau_fsd*Laminar_Flame_Speed*((0.5-C)*(Fsd*Mz*dWdz.rho+rho*Mz*dWdz.Fsd+rho*Fsd*Mzz)-rho*Fsd*Mz*dWdz.C);
    return ( NGT_fsd_x+NGT_fsd_y+NGT_fsd_z );
   }else{
    return (0.0);
   } 
}

double LES3DFsd_pState::SFS_Diffusion_Progvar (const LES3DFsd_pState &dWdx,
                                               const LES3DFsd_pState &dWdy,
                                               const LES3DFsd_pState &dWdz,
                                               const LES3DFsd_pState &d_dWdx_dx,
                                               const LES3DFsd_pState &d_dWdy_dy,
                                               const LES3DFsd_pState &d_dWdz_dz,
                                               const LES3DFsd_pState &d_dWdx_dy,
                                               const LES3DFsd_pState &d_dWdx_dz,
                                               const LES3DFsd_pState &d_dWdy_dz,
                                               const int &Flow_Type,
                                               const double &Volume) const {
  double grad_eddyviscosity_x, grad_eddyviscosity_y, grad_eddyviscosity_z, sfs_diffusion_progvar_x, sfs_diffusion_progvar_y, sfs_diffusion_progvar_z;
    double Cv = 0.0184, Schmidt_sfs = 1.0;
    grad_eddyviscosity_x = HALF*Cv*filter_width(Volume)*dWdx.k/sqrt(k);
    grad_eddyviscosity_y = HALF*Cv*filter_width(Volume)*dWdy.k/sqrt(k);
    grad_eddyviscosity_z = HALF*Cv*filter_width(Volume)*dWdz.k/sqrt(k);
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C!= ZERO && dWdz.C != ZERO ) {
    sfs_diffusion_progvar_x = (eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,Volume)*dWdx.C*dWdx.rho+rho*dWdx.C*grad_eddyviscosity_x+rho*eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,Volume)*d_dWdx_dx.C)/Schmidt_sfs;
    sfs_diffusion_progvar_y = (eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,Volume)*dWdy.C*dWdy.rho+rho*dWdy.C*grad_eddyviscosity_y+rho*eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,Volume)*d_dWdy_dy.C)/Schmidt_sfs;
    sfs_diffusion_progvar_z = (eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,Volume)*dWdz.C*dWdz.rho+rho*dWdz.C*grad_eddyviscosity_z+rho*eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,Volume)*d_dWdz_dz.C)/Schmidt_sfs;
    return ( sfs_diffusion_progvar_x+sfs_diffusion_progvar_y+sfs_diffusion_progvar_z );
   }else{
    return (0.0);
   }
 }

double LES3DFsd_pState::SFS_Diffusion_Fsd (const LES3DFsd_pState &dWdx,
                                           const LES3DFsd_pState &dWdy,
                                           const LES3DFsd_pState &dWdz,
                                           const LES3DFsd_pState &d_dWdx_dx,
                                           const LES3DFsd_pState &d_dWdy_dy,
                                           const LES3DFsd_pState &d_dWdz_dz,
                                           const LES3DFsd_pState &d_dWdx_dy,
                                           const LES3DFsd_pState &d_dWdx_dz,
                                           const LES3DFsd_pState &d_dWdy_dz,
                                           const int &Flow_Type,
                                           const double &Volume) const {
     double grad_eddyviscosity_x, grad_eddyviscosity_y, grad_eddyviscosity_z, sfs_diffusion_fsd_x, sfs_diffusion_fsd_y, sfs_diffusion_fsd_z;
     double Cv = 0.0184, Schmidt_sfs = 1.0;
    grad_eddyviscosity_x = HALF*Cv*filter_width(Volume)*dWdx.k/sqrt(k);
    grad_eddyviscosity_y = HALF*Cv*filter_width(Volume)*dWdy.k/sqrt(k);
    grad_eddyviscosity_z = HALF*Cv*filter_width(Volume)*dWdz.k/sqrt(k);
    if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C!= ZERO && dWdz.C != ZERO ) {
     sfs_diffusion_fsd_x = (eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,Volume)*dWdx.Fsd*dWdx.rho+rho*dWdx.Fsd*grad_eddyviscosity_x+rho*eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,Volume)*d_dWdx_dx.Fsd)/Schmidt_sfs;
     sfs_diffusion_fsd_y = (eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,Volume)*dWdy.Fsd*dWdy.rho+rho*dWdy.Fsd*grad_eddyviscosity_y+rho*eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,Volume)*d_dWdy_dy.Fsd)/Schmidt_sfs;
     sfs_diffusion_fsd_z = (eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,Volume)*dWdx.Fsd*dWdx.rho+rho*dWdx.Fsd*grad_eddyviscosity_x+rho*eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,Volume)*d_dWdx_dx.Fsd)/Schmidt_sfs;
     return ( sfs_diffusion_fsd_x+sfs_diffusion_fsd_y+sfs_diffusion_fsd_z );
    }else{
     return (0.0);
    }
 }

double LES3DFsd_pState::Heat_Release_Strain (const LES3DFsd_pState &dWdx,
                                             const LES3DFsd_pState &dWdy,
                                             const LES3DFsd_pState &dWdz,
                                             const LES3DFsd_pState &d_dWdx_dx,
                                             const LES3DFsd_pState &d_dWdy_dy,
                                             const LES3DFsd_pState &d_dWdz_dz,
                                             const LES3DFsd_pState &d_dWdx_dy,
                                             const LES3DFsd_pState &d_dWdx_dz,
 	                                     const LES3DFsd_pState &d_dWdy_dz) const {
     double tau_fsd, Mxx, Myy, Mzz, heat_release_strain_xx, heat_release_strain_yy, heat_release_strain_zz;
     tau_fsd = HeatRelease_Parameter();
     double Laminar_Flame_Speed = 0.3837;
     Mxx = M_xx(dWdx,dWdy,dWdz,d_dWdx_dx,d_dWdx_dy,d_dWdx_dz);
     Myy = M_yy(dWdx,dWdy,dWdz,d_dWdy_dy,d_dWdx_dy,d_dWdy_dz);
     Mzz = M_zz(dWdx,dWdy,dWdz,d_dWdz_dz,d_dWdx_dz,d_dWdy_dz);
     if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C!= ZERO && dWdz.C != ZERO ) {
     heat_release_strain_xx = (0.5-C)*tau_fsd*Laminar_Flame_Speed*Fsd*rho*Mxx;
     heat_release_strain_yy = (0.5-C)*tau_fsd*Laminar_Flame_Speed*Fsd*rho*Myy;
     heat_release_strain_zz = (0.5-C)*tau_fsd*Laminar_Flame_Speed*Fsd*rho*Mzz;
     return ( heat_release_strain_xx + heat_release_strain_yy + heat_release_strain_zz );
    }else{
     return (0.0);
    }
}

double LES3DFsd_pState::Net_Rate_Change_Progvar (const LES3DFsd_pState &dWdx,
                                                 const LES3DFsd_pState &dWdy,
                                                 const LES3DFsd_pState &dWdz,
                                                 const LES3DFsd_pState &d_dWdx_dx,
                                                 const LES3DFsd_pState &d_dWdy_dy,
                                                 const LES3DFsd_pState &d_dWdz_dz,
                                                 const LES3DFsd_pState &d_dWdx_dy,
                                                 const LES3DFsd_pState &d_dWdx_dz,
						 const LES3DFsd_pState &d_dWdy_dz,
                                                 const int &Flow_Type,
                                                 const double &Volume) const {
  return(  Resolved_Convection_Progvar(dWdx,dWdy,dWdz)
	   +SFS_Diffusion_Progvar(dWdx,dWdy,dWdz,d_dWdx_dx,d_dWdy_dy,d_dWdz_dz,d_dWdx_dy,d_dWdx_dz,d_dWdy_dz,Flow_Type,Volume)
	   +Reaction_Rate_Fsd(dWdx,dWdy,dWdz) );
}
           

double LES3DFsd_pState::Net_Rate_Change_Fsd (const LES3DFsd_pState &dWdx,
                                             const LES3DFsd_pState &dWdy,
                                             const LES3DFsd_pState &dWdz,
                                             const LES3DFsd_pState &d_dWdx_dx,
                                             const LES3DFsd_pState &d_dWdy_dy,
                                             const LES3DFsd_pState &d_dWdz_dz,
                                             const LES3DFsd_pState &d_dWdx_dy,
                                             const LES3DFsd_pState &d_dWdx_dz,
			  		     const LES3DFsd_pState &d_dWdy_dz,
                                             const int &Flow_Type,
                                             const double &Volume) const {
  return(  Resolved_Convection_Fsd(dWdx,dWdy,dWdz)
          +SFS_Diffusion_Fsd(dWdx,dWdy,dWdz,d_dWdx_dx,d_dWdy_dy,d_dWdz_dz,d_dWdx_dy,d_dWdx_dz,d_dWdy_dz,Flow_Type,Volume)
	  +Resolved_Strain(dWdx,dWdy,dWdz)
	  +Resolved_Propagation_Curvature(dWdx,dWdy,dWdz)
	  +SFS_Strain(dWdx,dWdy,dWdz,Flow_Type,Volume)
	  +SFS_Curvature(dWdx,dWdy,dWdz) );
}
/*
double LES3DFsd_pState::KPP_Tur_Speed(const LES3DFsd_pState &dWdx,
                                    const LES3DFsd_pState &dWdy,
                                    const int &Flow_Type,
                                    const double &Volume) const {
    double Mx = M_x(dWdx,dWdy);
    if ( C < 0.99 && C > 0.01 && dWdx.C != ZERO && dWdy.C != ZERO && dWdz.C != ZERO ) {
    return (-Mx*Laminar_Flame_Speed+2.0*sqrt(SFS_Strain(dWdx,dWdy,Flow_Type,Volume)*eddy_viscosity(dWdx,dWdy,Volume,Flow_Type)));
   }else{
    return(0.0);
   }
}
*/
double LES3DFsd_pState::K_equ_sources(const LES3DFsd_pState &dWdx,
                                      const LES3DFsd_pState &dWdy,
                                      const LES3DFsd_pState &dWdz,
                                      const int &Flow_Type,
                                      const double &Volume) const {

   double production, dissipation, source;
   double mu_t;
   Tensor3D subfilter_stress;
   double Ceps = 0.845;
   //Turbulence model eddy viscosity
   mu_t = eddy_viscosity(dWdx,dWdy,dWdz,Flow_Type,Volume);
   subfilter_stress = lambda(dWdx, dWdy, dWdz, Flow_Type, Volume);
  
     if ( C < 0.999 && C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C!= ZERO && dWdz.C != ZERO ) {
   production = subfilter_stress.xx*dWdx.v.x + 
      subfilter_stress.xy*(dWdy.v.x + dWdx.v.y) +
      subfilter_stress.yy*dWdy.v.y +
      subfilter_stress.xz*(dWdz.v.x + dWdx.v.z) +
      subfilter_stress.yz*(dWdz.v.y + dWdy.v.z) +
      subfilter_stress.zz*dWdz.v.z;

   dissipation = Ceps*rho*pow(k, 1.5)/filter_width(Volume);
   source = production - dissipation;
   return(source);
    }else{
     return (0.0);
    }
}

/* Turbulence model source term */  
LES3DFsd_cState  LES3DFsd_pState::Src_t(const LES3DFsd_pState &Wc,
                                        const LES3DFsd_pState &dWdx,
                                        const LES3DFsd_pState &dWdy,
                                        const LES3DFsd_pState &dWdz,
                                        const int &Flow_Type,
                                        const double &Volume) {
   
  double  resolved_strain, resolved_propagation_curvature, sfs_strain, sfs_curvature;
  LES3DFsd_cState Temp; Temp.Vacuum();

  if ( Wc.C < 0.999 && Wc.C > 0.001 && dWdx.C != ZERO ) {//&& dWdy.C != ZERO && dWdz.C != ZERO ) {

    // Reaction Rate for Progress Variable Equation --- Source term
      Temp.rhoC = Wc.Reaction_Rate_Fsd(dWdx,dWdy,dWdz);
      //      cout<<"rhoc=   "<<Wc.Fsd*Wc.rho<<endl;
    // FSD Equation Source Term
      resolved_strain = Wc.Resolved_Strain(dWdx,dWdy,dWdz);

      resolved_propagation_curvature = Wc.Resolved_Propagation_Curvature(dWdx,dWdy,dWdz);

      sfs_strain = Wc.SFS_Strain(dWdx,dWdy,dWdz,Flow_Type,Volume);

      sfs_curvature = Wc.SFS_Curvature(dWdx,dWdy,dWdz);

      Temp.rhoFsd = resolved_strain + resolved_propagation_curvature + sfs_strain + sfs_curvature ;
      if (Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K) {
    // k-equation Source Term
    Temp.rhok = Wc.K_equ_sources(dWdx,dWdy,dWdz,Flow_Type,Volume);
      }

    }
  return (Temp);


} 

/*************************************************
         Premixed combustion
*************************************************/

/***** Mass fractions for CH4 one step *****/

LES3DFsd_pState LES3DFsd_pState::premixed_mfrac(const LES3DFsd_pState &Wo,const LES3DFsd_pState &W){

       LES3DFsd_pState temp;
       temp.Copy (*this);       

      
  double unburnt_fuel_c, burnt_fuel_c, tol;
  double burnt_oxygen_c, stoich_ratio, phi;
  double c_products, products_ratio;    // c_prod = c_CO2 + c_H2O
  
  tol = MICRO;          // local tolerance
  if ( React.reactset_flag == CH4_1STEP ){
  stoich_ratio = 2.0*temp.specdata[1].Mol_mass()/temp.specdata[0].Mol_mass();  //      4.0;   // stoichiometric O2/CH4 mass ratio
  temp.spec[4].c = Wo.spec[4].c;
  }else if ( React.reactset_flag == C3H8_1STEP ){
  stoich_ratio = 5.0*temp.specdata[1].Mol_mass()/temp.specdata[0].Mol_mass();  //      3.6;   // stoichiometric O2/C3H8 mass ratio
  temp.spec[4].c = Wo.spec[4].c;
  }else if ( React.reactset_flag == H2O2_1STEP ){
  stoich_ratio = temp.specdata[1].Mol_mass()/2.0/temp.specdata[0].Mol_mass();  //      8.0;   // stoichiometric O2/H2 mass ratio
  temp.spec[3].c = Wo.spec[3].c;
  }
  //  cout<<"c=   "<<W.C<<endl;
    if ( temp.C >= ONE ) {
      temp.C = ONE;
    }
    if ( temp.C <= ZERO ) {
      temp.C = ZERO;
    }
    temp.spec[0].c = ( ONE - temp.C )*Wo.spec[0].c;
    //    cout<<"C=   "<<temp.C<<"   CH4=  "<<temp.spec[0].c<<endl;  
  // equivalence ratio
    phi = 1.0;//Equivalence_Ratio;//stoich_ratio/(Wo.spec[1].c / Wo.spec[0].c);

  if(fabs(phi-ONE)<1.0E-2){//  || (phi-ONE)<-1.0E-3){
     phi = ONE;
   }

  // check for negative or small fuel mass fraction
  if(temp.spec[0].c < tol){
    temp.spec[0].c = ZERO; }
  
  // lean mixture(phi < 1) => excessive O2
  if(phi < ONE){  
    unburnt_fuel_c = Wo.spec[0].c; // initial fuel mass fraction
    burnt_fuel_c = ZERO;
    burnt_oxygen_c = (ONE/phi - ONE)*unburnt_fuel_c*stoich_ratio;
    if(temp.spec[0].c == ZERO){
      temp.spec[1].c = burnt_oxygen_c;      
    }else{
      temp.spec[1].c = temp.spec[0].c * stoich_ratio + burnt_oxygen_c;///phi;
    }  

  // rich mixture(phi > 1) => excessive CH4
  }else if(phi > ONE){  
    unburnt_fuel_c = Wo.spec[0].c;  // initial fuel mass fraction
    burnt_oxygen_c = ZERO;
    burnt_fuel_c = (ONE - ONE/phi)*unburnt_fuel_c;///stoich_ratio;
    temp.spec[0].c = temp.C*(burnt_fuel_c - unburnt_fuel_c) + unburnt_fuel_c;
    if(temp.spec[0].c <= burnt_fuel_c){
      temp.spec[1].c = burnt_oxygen_c;
    }else{
      temp.spec[1].c = (temp.spec[0].c - burnt_fuel_c) * stoich_ratio;///phi;
    }

  // stoichiometric mixture(phi = 1)
  }else{ 
    burnt_fuel_c = ZERO;
    burnt_oxygen_c = ZERO;
    temp.spec[1].c = temp.spec[0].c * stoich_ratio; 
   }

  if ( React.reactset_flag == CH4_1STEP ||
       React.reactset_flag == C3H8_1STEP ){
    // mass fractions of products
       c_products = ONE - (temp.spec[0].c + temp.spec[1].c + temp.spec[4].c);
  if ( React.reactset_flag == CH4_1STEP ){
       products_ratio = temp.specdata[2].Mol_mass()/(temp.specdata[2].Mol_mass()+2.0*temp.specdata[3].Mol_mass());
  }else if ( React.reactset_flag == C3H8_1STEP ){
       products_ratio = 3.0*temp.specdata[2].Mol_mass()/(3.0*temp.specdata[2].Mol_mass()+4.0*temp.specdata[3].Mol_mass());
  }
       temp.spec[2].c = products_ratio*c_products; // CO2 mass fraction
       temp.spec[3].c = c_products-temp.spec[2].c;  // H2O mass fraction   
       temp.spec[1].c = ONE-c_products - temp.spec[0].c - temp.spec[4].c ;      // O2 mass fraction
  //  temp.spec[4].c = ONE-temp.spec[0].c-temp.spec[1].c-temp.spec[2].c-temp.spec[3].c;
  }else if ( React.reactset_flag == H2O2_1STEP ){
       temp.spec[2].c = ONE - temp.spec[0].c - temp.spec[1].c - temp.spec[3].c;// H2O mass fraction 
  }

 /*
double mu_f = 1.0;
double mu_O2 = 2.0;
double mu_N2 = 2.0*3.76;
double mu_CO2 = 1.0;
double mu_H2O = 2.0;
double phi = 1.0;
double stoich_ratio, mf_ub, mO2_ub, mN2, yf_ub,  yO2_ub, yN2, mf_b, mO2_b, mCO2_b, mH2O_b,
       yf_b, yO2_b, yCO2_b, yH2O_b;
double sum_ub, sum_b, sum_y_ub, sum_y_b;
	
stoich_ratio = mu_O2*temp.specdata[1].Mol_mass()/(mu_f*temp.specdata[0].Mol_mass());

mf_ub = 1.0;
mO2_ub = stoich_ratio/phi;
mN2 = (mu_N2*temp.specdata[4].Mol_mass())/(mu_f*temp.specdata[0].Mol_mass())/phi;
sum_ub = mf_ub + mO2_ub + mN2;
yf_ub = mf_ub/sum_ub;
yO2_ub = stoich_ratio/phi*yf_ub;
yN2 = mN2/sum_ub;
sum_y_ub = yf_ub + yO2_ub + yN2;
mf_b = max(0.0, mf_ub-mO2_ub/stoich_ratio);
mO2_b = max(0.0, mO2_ub-mf_ub*stoich_ratio);
mCO2_b = (mf_ub-mf_b)*(mu_CO2*temp.specdata[2].Mol_mass())/(mu_f*temp.specdata[0].Mol_mass());
mH2O_b = (mf_ub-mf_b)*(mu_H2O*temp.specdata[3].Mol_mass())/(mu_f*temp.specdata[0].Mol_mass());
sum_b = mf_b + mO2_b + mCO2_b + mH2O_b + mN2;
yf_b = mf_b/sum_b;
yO2_b = mO2_b/sum_b;
yCO2_b = mCO2_b/sum_b;
yH2O_b = mH2O_b/sum_b;
sum_y_b = yf_b + yO2_b + yCO2_b + yH2O_b + yN2;

temp.spec[0].c =  yf_ub + ( yf_b - yf_ub ) * temp.C;
temp.spec[1].c = (stoich_ratio/phi)*temp.spec[0].c;
temp.spec[2].c = (1.0-temp.spec[0].c-temp.spec[1].c-yN2)
	            *(mu_CO2*temp.specdata[2].Mol_mass())/(mu_CO2*temp.specdata[2].Mol_mass()+mu_H2O*temp.specdata[3].Mol_mass());
temp.spec[3].c = (1.0-temp.spec[0].c-temp.spec[1].c-yN2)
	            *(mu_H2O*temp.specdata[3].Mol_mass())/(mu_CO2*temp.specdata[2].Mol_mass()+mu_H2O*temp.specdata[3].Mol_mass());
temp.spec[4].c = yN2; 
 */	    

  double suma;
  suma = 0.0;
  for(int i=0; i<ns; i++){
    suma = suma + temp.spec[i].c;
  }
  for(int i=0; i<ns; i++){
      temp.spec[i].c = temp.spec[i].c*(ONE/suma);
  }

 for (int i=0; i<ns; i++){
   if ( temp.spec[i].c< ZERO ){
     temp.spec[i].c = ZERO;
   }
 }

 return(temp);
}

/*************************************************
         Premixed combustion
*************************************************/


/***** Mass fractions for CH4 one step *****/

LES3DFsd_cState  LES3DFsd_cState::premixed_mfrac(const LES3DFsd_pState &Wo){
 
  LES3DFsd_cState temp;
  temp.Copy(*this);

  
  double unburnt_fuel_c, burnt_fuel_c, tol, sum;
  double burnt_oxygen_c, stoich_ratio, phi;
  double c_products, products_ratio;    // c_prod = c_CO2 + c_H2O

  tol = MICRO;          // local tolerance
  if ( Wo.React.reactset_flag == CH4_1STEP ){
  stoich_ratio = 2.0*temp.specdata[1].Mol_mass()/temp.specdata[0].Mol_mass();  //      4.0;   // stoichiometric O2/CH4 mass ratio
  temp.rhospec[4].c = temp.rho*Wo.spec[4].c;
  }else if ( Wo.React.reactset_flag == C3H8_1STEP ){
  stoich_ratio = 5.0*temp.specdata[1].Mol_mass()/temp.specdata[0].Mol_mass();  //      4.0;   // stoichiometric O2/C3H8 mass ratio
  temp.rhospec[4].c = temp.rho*Wo.spec[4].c;
  }else if ( Wo.React.reactset_flag == H2O2_1STEP ){
  stoich_ratio = temp.specdata[1].Mol_mass()/2.0/temp.specdata[0].Mol_mass();  //      8.0;   // stoichiometric O2/H2 mass ratio
  temp.rhospec[3].c = temp.rho*Wo.spec[3].c;
  }

  if ( temp.rhoC >= temp.rho ) {
      temp.rhoC = temp.rho;
    }
    if ( temp.rhoC <= ZERO ) {
      temp.rhoC = ZERO;
    }
    temp.rhospec[0].c = ( ONE - temp.rhoC/temp.rho )*Wo.spec[0].c*temp.rho;
  // equivalence ratio
    phi = 1.0;//Equivalence_Ratio;//stoich_ratio/(Wo.spec[1].c / Wo.spec[0].c);

  if(fabs(phi-ONE)<1.0E-2){//  || (phi-ONE)<-1.0E-3){
     phi = ONE;
   }

  // check for negative or small fuel mass fraction
  if(temp.rho < ZERO){
    cout << "Negative density: " << temp.rho << endl;
  }
  
  if(temp.rhospec[0].c/temp.rho < ZERO )
    { temp.rhospec[0].c = ZERO; }
  if(temp.rhospec[0].c/temp.rho > Wo.spec[0].c)
    { temp.rhospec[0].c = Wo.spec[0].c * temp.rho; }
  
  // lean mixture(phi < 1) => excessive O2
  if(phi < ONE){  
     unburnt_fuel_c = Wo.spec[0].c; // initial fuel mass fraction
     burnt_fuel_c = ZERO;
     burnt_oxygen_c = (ONE/phi - ONE)*unburnt_fuel_c*stoich_ratio;
  if(temp.rhospec[0].c < tol){
     temp.rhospec[0].c = ZERO;
     temp.rhospec[1].c = temp.rho*burnt_oxygen_c;      
  }else{
     temp.rhospec[1].c = temp.rhospec[0].c * stoich_ratio + temp.rho * burnt_oxygen_c;///phi;
    }  

  // rich mixture(phi > 1) => excessive CH4
  }else if(phi > ONE){  
    unburnt_fuel_c = Wo.spec[0].c;  // initial fuel mass fraction
    burnt_oxygen_c = ZERO;
    burnt_fuel_c = (ONE - ONE/phi)*unburnt_fuel_c;///stoich_ratio;
    temp.rhospec[0].c = temp.rhoC*(burnt_fuel_c - unburnt_fuel_c) + unburnt_fuel_c*temp.rho;
 if(temp.rhospec[0].c <= temp.rho*burnt_fuel_c){
    temp.rhospec[1].c = temp.rho*burnt_oxygen_c;
  }else{
    temp.rhospec[1].c = (temp.rhospec[0].c - burnt_fuel_c*temp.rho)* stoich_ratio;///phi;
    }
  // stoichiometric mixture(phi = 1)
  }else{
    burnt_fuel_c = 0.0;
    burnt_oxygen_c = 0.0;
    temp.rhospec[1].c = temp.rhospec[0].c * stoich_ratio;
 }

  if ( Wo.React.reactset_flag == CH4_1STEP ||
       Wo.React.reactset_flag == C3H8_1STEP ){
  // mass fractions of products
       c_products = ONE - (temp.rhospec[0].c + temp.rhospec[1].c + temp.rhospec[4].c)/temp.rho;
  if ( Wo.React.reactset_flag == CH4_1STEP ){
       products_ratio = temp.specdata[2].Mol_mass()/(temp.specdata[2].Mol_mass()+2.0*temp.specdata[3].Mol_mass());   
  }else if ( Wo.React.reactset_flag == C3H8_1STEP ){
       products_ratio = 3.0*temp.specdata[2].Mol_mass()/(3.0*temp.specdata[2].Mol_mass()+4.0*temp.specdata[3].Mol_mass());
    }   
    temp.rhospec[2].c = products_ratio*c_products*temp.rho; // CO2 mass fraction
    temp.rhospec[3].c = temp.rho*c_products - temp.rhospec[2].c;      // H2O mass fraction
    temp.rhospec[1].c = temp.rho*(ONE-c_products - (temp.rhospec[0].c+temp.rhospec[4].c)/temp.rho);      // O2 mass fraction
    //    temp.rhospec[4].c = temp.rho - temp.rhospec[0].c - temp.rhospec[1].c - temp.rhospec[2].c - temp.rhospec[3].c; 
  }else if ( Wo.React.reactset_flag == CH4_1STEP ){
    temp.rhospec[2].c = c_products*temp.rho; // H2O mass fraction  
  }
  //    }  
    /*       
double mu_f = 1.0;
double mu_O2 = 2.0;
double mu_N2 = 2.0*3.76;
double mu_CO2 = 1.0;
double mu_H2O = 2.0;
double stoich_ratio, mf_ub, mO2_ub, mN2, yf_ub,  yO2_ub, yN2, mf_b, mO2_b, mCO2_b, mH2O_b,
       yf_b, yO2_b, yCO2_b, yH2O_b;
double sum_ub, sum_b, sum_y_ub, sum_y_b;
	
stoich_ratio = mu_O2*temp.specdata[1].Mol_mass()/(mu_f*temp.specdata[0].Mol_mass());
 double phi = 1.0;
//double  phi = stoich_ratio/(Wo.spec[1].c / Wo.spec[0].c);

//  if((phi-ONE)<1.0E-3  || (phi-ONE)<-1.0E-3){
//    phi = ONE;
//  }

  // check for negative or small fuel mass fraction
  if(temp.rho < ZERO){
    cout << "Negative density: " << temp.rho << endl;
  }
  
mf_ub = 1.0;
mO2_ub = stoich_ratio/phi;
mN2 = (mu_N2*temp.specdata[4].Mol_mass())/(mu_f*temp.specdata[0].Mol_mass())/phi;
sum_ub = mf_ub + mO2_ub + mN2;
yf_ub = mf_ub/sum_ub;
yO2_ub = stoich_ratio/phi*yf_ub;
yN2 = mN2/sum_ub;
sum_y_ub = yf_ub + yO2_ub + yN2;
mf_b = max(0.0, mf_ub-mO2_ub/stoich_ratio);
mO2_b = max(0.0, mO2_ub-mf_ub*stoich_ratio);
mCO2_b = (mf_ub-mf_b)*(mu_CO2*temp.specdata[2].Mol_mass())/(mu_f*temp.specdata[0].Mol_mass());
mH2O_b = (mf_ub-mf_b)*(mu_H2O*temp.specdata[3].Mol_mass())/(mu_f*temp.specdata[0].Mol_mass());
sum_b = mf_b + mO2_b + mCO2_b + mH2O_b + mN2;
yf_b = mf_b/sum_b;
yO2_b = mO2_b/sum_b;
yCO2_b = mCO2_b/sum_b;
yH2O_b = mH2O_b/sum_b;
sum_y_b = yf_b + yO2_b + yCO2_b + yH2O_b + yN2;
double ratio1 = (mu_CO2*temp.specdata[2].Mol_mass())/(mu_CO2*temp.specdata[2].Mol_mass()+mu_H2O*temp.specdata[3].Mol_mass());
double ratio2 = (mu_H2O*temp.specdata[3].Mol_mass())/(mu_CO2*temp.specdata[2].Mol_mass()+mu_H2O*temp.specdata[3].Mol_mass());
 temp.rhospec[0].c = (yf_ub + ( yf_b - yf_ub ) * temp.rhoC/temp.rho)*temp.rho; 		    
 temp.rhospec[1].c = (stoich_ratio/phi)*temp.rhospec[0].c; 		    
 temp.rhospec[2].c = (temp.rho-temp.rhospec[0].c-temp.rhospec[1].c-yN2*temp.rho)*ratio1; 		    
 temp.rhospec[3].c = (temp.rho-temp.rhospec[0].c-temp.rhospec[1].c-yN2*temp.rho-temp.rhospec[2].c);//*ratio2; 		    
 // temp.rhospec[4].c = yN2*temp.rho;
 */

  double suma;
  suma = 0.0;
  for(int i=0; i<ns; i++){
    suma = suma + temp.rhospec[i].c/temp.rho;
  }
  for(int i=0; i<ns; i++){
      temp.rhospec[i].c = temp.rhospec[i].c*(ONE/suma);
  }

//    if ( Flow_Type == FLOWTYPE_LAMINAR_C || 
//         Flow_Type == FLOWTYPE_LAMINAR_C_FSD ||
//         Flow_Type == FLOWTYPE_TURBULENT_LES_C ||
//         Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD ||
//         Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K ||
//         Flow_Type == FLOWTYPE_FROZEN_TURBULENT_LES_C_FSD)  {
//         temp.rhoC = min(temp.rho, max(ZERO, temp.rho*(Wo.spec[0].c-temp.rhospec[0].c/temp.rho)/Wo.spec[0].c));

//    }
  
  for (int i=0; i<ns; i++){
    if (temp.rhospec[i].c/temp.rho < ZERO){
      temp.rhospec[i].c = ZERO;
    }
  }

  return temp;

}

