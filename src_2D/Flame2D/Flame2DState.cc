/////////////////////////////////////////////////////////////////////
///
/// \file Flame2DState.cc
/// 
/// \author Marc R.J. Charest
/// 
/// \brief This header file contains the class definition for 
///        describing the state of the chemically reacting, 
///        multispecies ideal gas mixture.
///
/////////////////////////////////////////////////////////////////////
#include "Flame2DState.h"

/**
 * Initialization of static variables.
 */
#ifndef STATIC_NUMBER_OF_SPECIES
int Flame2D_State :: n = 0;
int Flame2D_State :: ns = 0;
#endif
bool Flame2D_State::reacting = false;
double Flame2D_State::Mref = 0.5;
double Flame2D_State::gravity_z = -9.81;
double* Flame2D_State::y = NULL;
double* Flame2D_pState::r = NULL;
double* Flame2D_pState::dihdic = NULL;
double* Flame2D_pState::hs_i = NULL;


/////////////////////////////////////////////////////////////////////
/// Static Setup Functions
/////////////////////////////////////////////////////////////////////

/****************************************************
 * Mixture molecular mass [kg/mol]
 ****************************************************/
void Flame2D_pState::setMixture(const string &mech_name,
				const string &mech_file) {

  // call mixture object setup functin
  Mixture::setMixture(mech_name, mech_file);

  // assign static values
#ifndef STATIC_NUMBER_OF_SPECIES
  ns = Mixture::nSpecies();
  n = ns + NUM_FLAME2D_VAR_SANS_SPECIES;
#endif

  //allocate static memory and load the species data  
  AllocateStatic();

}

/**********************************************************************
 * Flame2D_State::set_gravity -- Set the acceleration due to gravity  *
 *                               in m/s^2.  It acts downwards in the  *
 *                               z-dir (g <= 0)                       *
 **********************************************************************/
void Flame2D_State::set_gravity(const double &g) { // [m/s^2]

  // if gravity is acting upwards (wrong way)
  if (g>0) {
    cerr<<"\n Flame2D_pState::set_gravity() - Gravity acting upwards!!!! \n";
    exit(1);
    
  // gravity acting downwards (-ve), OK
  } else {
    gravity_z = g;
  }
}

/////////////////////////////////////////////////////////////////////
/// Flux Functions
/////////////////////////////////////////////////////////////////////

/****************************************************
 * X-Dir Conserved Flux
 ****************************************************/  
void Flame2D_pState::Fx(Flame2D_State &FluxX) const {
  FluxX.rho() = rho()*vx();
  FluxX.rhovx() = rho()*sqr(vx()) + p();
  FluxX.rhovy() = rho()*vx()*vy();
  FluxX.E() = vx()*H();
  
  //multispecies transport
  for(int i=0; i<ns; i++){
    FluxX.rhoc(i) = rho()*vx()*c(i);
  }
}

void Flame2D_pState::addFx(Flame2D_State &FluxX, const double& mult) const {
  FluxX.rho() += mult*rho()*vx();
  FluxX.rhovx() += mult*rho()*sqr(vx()) + p();
  FluxX.rhovy() += mult*rho()*vx()*vy();
  FluxX.E() += mult*vx()*H();
  
  //multispecies transport
  for(int i=0; i<ns; i++){
    FluxX.rhoc(i) += mult*rho()*vx()*c(i);
  }
}

Flame2D_State Flame2D_pState::Fx(void) const {
  Flame2D_State FluxX;
  FluxX.rho() = rho()*vx();
  FluxX.rhovx() = rho()*sqr(vx()) + p();
  FluxX.rhovy() = rho()*vx()*vy();
  FluxX.E() = vx()*H();
  
  //multispecies transport
  for(int i=0; i<ns; i++){
    FluxX.rhoc(i) = rho()*vx()*c(i);
  }

  return FluxX;
}

/*****************************************************************
 * Viscous fluxes  (laminar flow)                                * 
 * Viscous fluxes  (turbulent flows) are defined in single block * 
 ****************************************************************/
void Flame2D_pState::Viscous_Flux_x(const Flame2D_State &dWdx,
				    const Vector2D &qflux,
				    const Tensor2D &tau,
				    Flame2D_State &Flux,
				    const double& mult) const{
 
  //Flux.rho() += mult * ZERO;
  Flux.rhovx() += mult * tau.xx;
  Flux.rhovy() += mult * tau.xy;
  Flux.E() += mult * ( -qflux.x + vx()*tau.xx + vy()*tau.xy );
  //rho * Diffusion_Coef * grad cn 
  for( int i=0; i<ns; i++){
    Flux.rhoc(i) += mult * (Diffusion_coef(i) * dWdx.c(i)); 
  }

}

void Flame2D_pState::Viscous_Flux_y(const Flame2D_State &dWdy,
				    const Vector2D &qflux,
				    const Tensor2D &tau, 
				    Flame2D_State &Flux,
				    const double& mult) const {

  //Flux.rho() += mult * ZERO;
  Flux.rhovx() += mult * tau.xy; 
  Flux.rhovy() += mult * tau.yy;
  Flux.E() += mult * ( -qflux.y + vx()*tau.xy + vy()*tau.yy );
  //rho * Diffusion_Coef * grad cn 
  for( int i=0; i<ns; i++){
    Flux.rhoc(i) += mult * (Diffusion_coef(i) * dWdy.c(i));     
  }

}


/////////////////////////////////////////////////////////////////////
/// EIGENVALUES
/////////////////////////////////////////////////////////////////////

/************************************************************
 * Flame2D_pState::lambda -- Eigenvalue(s) (x-direction).    *
 ************************************************************/
void Flame2D_pState::lambda_x(Flame2D_State &lambdas) const {
  double aa(a());
  lambdas.rho() = vx() - aa;
  lambdas.vx() = vx();
  lambdas.vy() = vx();
  lambdas.p() = vx() + aa;
  for(int i=0; i<ns; i++) lambdas.c(i) = vx();
}


/************************************************************
 * Low Mach Number Preconditioned Eigenvalue(s)             *
 * Flame2D_pState::lambda_preconditioned_x                   *
 ************************************************************/
void Flame2D_pState::lambda_preconditioned_x(Flame2D_State &lambdas, 
					     const double &MR2) const {

  double uprimed, cprimed;
  double aa(a());
  u_a_precon(MR2*aa*aa,uprimed,cprimed);

  lambdas.rho() = uprimed - cprimed;
  lambdas.vx() = vx();
  lambdas.vy() = vx();
  lambdas.p() = uprimed + cprimed;
  for(int i=0; i<ns; i++){
    lambdas.c(i) = vx();
  }
}


/////////////////////////////////////////////////////////////////////
/// EIGENVECTORS
/////////////////////////////////////////////////////////////////////

// Conserved Right Eigenvector -- (x-direction)
Flame2D_State Flame2D_pState::rc_x(const int &index) const {

  if(index == 1){
    double aa(a()); 
    return (Flame2D_State(ONE, vx()-aa, vy(), H()/rho()-vx()*aa, c()));
  } else if(index == 2) {
    return (Flame2D_State(ONE, vx(), vy(), H()/rho()-Cp()*T(), c())); 
  } else if(index == 3) {
    return (Flame2D_State(ZERO, ZERO, rho(), rho()*vy(), ZERO));
  } else if(index == 4) {
    double aa(a());
    return (Flame2D_State(ONE, vx()+aa, vy(), H()/rho()+vx()*aa, c()));
  } else{ 
    static Flame2D_State tmp;
    tmp.Vacuum();
    int k( (index-1)-NUM_FLAME2D_VAR_SANS_SPECIES );
    tmp.E() = rho()*dihdic[k];  // -> note - dihdic needs to be loaded first outside
    tmp.rhoc(k) = rho();
    return tmp;
  }
}

// Primitive Left Eigenvector -- (x-direction)
Flame2D_State Flame2D_pState::lp_x(const int &index) const {
 
  if(index == 1){
    double aa(a());
    return (Flame2D_State(ZERO, -HALF*rho()/aa, ZERO, HALF/(aa*aa), ZERO));
  } else if(index == 2) {
    double aa(a());
    return (Flame2D_State(ONE, ZERO, ZERO, -ONE/(aa*aa), ZERO));
  } else if(index == 3) {
    return  (Flame2D_State(ZERO, ZERO, ONE, ZERO, ZERO));
  } else if(index == 4) {  
     double aa(a());
     return (Flame2D_State(ZERO, HALF*rho()/aa, ZERO, HALF/(aa*aa), ZERO));
  } else{ 
    static Flame2D_State tmp;
    tmp.Vacuum();
    int k( (index-1)-NUM_FLAME2D_VAR_SANS_SPECIES );
    tmp.rhoc(k) = ONE;
    return tmp;
  } 
  
}

/////////////////////////////////////////////////////////////////////
/// Low Mach Number Preconditioner Based on Weiss & Smith (1995)
/////////////////////////////////////////////////////////////////////

// Conserved Right Eigenvector -- (x-direction)
Flame2D_State Flame2D_pState::rc_x_precon(const int &index, const double &MR2) const {
  
  if(index == 1){
    double aa(a());
    double uprimed,cprimed;
    u_a_precon(MR2*aa*aa,uprimed,cprimed);
    return (Flame2D_State(ONE, 
			  (uprimed-cprimed)/MR2,
			  vy(),			
			  h()+HALF*(vsqr()/MR2) - (vx()*cprimed)/MR2,
			  c()));
    
  } else if(index == 2) {
    return (Flame2D_State(ONE, vx(), vy(), (h()-Cp()*T()) + HALF*vsqr(), c()));

  } else if(index == 3) {
    return (Flame2D_State(ZERO, ZERO, rho(), rho()*vy(),ZERO));
    
  } else if(index == 4) { 
    double aa(a());
    double uprimed,cprimed;
    u_a_precon(MR2*aa*aa,uprimed,cprimed);
    return (Flame2D_State(ONE,
			  (uprimed+cprimed)/MR2,
			  vy(), 
			  h()+HALF*(vsqr()/MR2) + (vx()*cprimed)/MR2,
			  c()));
    
  } else{ 
    static Flame2D_State tmp;
    tmp.Vacuum();
    int k( (index-1)-NUM_FLAME2D_VAR_SANS_SPECIES );
    tmp.E() = rho()*dihdic[k];  // -> note - dihdic needs to be loaded first outside
    tmp.rhoc(k) = rho();
    return tmp;
   }

}

// Primitive Left Eigenvector -- (x-direction)
Flame2D_State Flame2D_pState::lp_x_precon(const int &index, const double &MR2) const {

  if(index == 1){
    double aa(a());
    double uprimed,cprimed;
    u_a_precon(MR2*aa*aa,uprimed,cprimed);
    return (Flame2D_State(ZERO, 
			  -HALF*rho()*MR2/cprimed, 
			  ZERO,
			  (-uprimed+cprimed + vx())/(TWO*cprimed*aa*aa),
			  ZERO));
  } else if(index == 2) {
    double aa(a());
    return (Flame2D_State(ONE, ZERO, ZERO, -ONE/(aa*aa),ZERO));
  } else if(index == 3) {
    return  (Flame2D_State(ZERO, ZERO, ONE, ZERO,ZERO));
  } else if(index == 4) {  
    double aa(a());
    double uprimed,cprimed;
    u_a_precon(MR2*aa*aa,uprimed,cprimed);
    return (Flame2D_State(ZERO, 
			  HALF*rho()*MR2/cprimed, 
			  ZERO,
			  (uprimed+cprimed - vx())/(TWO*cprimed*aa*aa),
			  ZERO));
  } else{ 
    static Flame2D_State tmp;
    tmp.Vacuum();
    tmp.c(index-(NUM_FLAME2D_VAR_SANS_SPECIES+1)) = ONE;
    return tmp;
  } 
}

void Flame2D_pState::Flux_Dissipation_precon(const double &MR2, 
					     const Flame2D_State &dWrl, 
					     const Flame2D_State &wavespeeds, 
					     Flame2D_State &Flux_dissipation) const {

  //declares
  static Flame2D_State lp;
  static Flame2D_State rc;
  double aa(a());
  double uprimed,cprimed;
  double tmp, dv;

  // compute preconditioned velocity and soundspeed
  u_a_precon(MR2*aa*aa,uprimed,cprimed);

  //
  // Loop over all the waves, and compute the flux dissipation
  //
  Flux_dissipation.Vacuum();
  for ( int i = 1 ; i <= n ; i++ ) {

    // Compute the left primitive eigenvector
    if(i == 1){
      lp.rho() = ZERO;
      lp.vx() = -HALF*rho()*MR2/cprimed;
      lp.vy() = ZERO;
      lp.p() = (-uprimed+cprimed + vx())/(TWO*cprimed*aa*aa);
      lp.c(ZERO);
    } else if(i == 2) {
      lp.rho() = ONE;
      lp.vx() = ZERO;
      lp.vy() = ZERO;
      lp.p() = -ONE/(aa*aa);
      lp.c(ZERO);
    } else if(i == 3) {
      lp.rho() = ZERO;
      lp.vx() = ZERO;
      lp.vy() = ONE;
      lp.p() = ZERO;
      lp.c(ZERO);
    } else if(i == 4) {  
      lp.rho() = ZERO;
      lp.vx() = HALF*rho()*MR2/cprimed;
      lp.vy() = ZERO;
      lp.p() = (uprimed+cprimed - vx())/(TWO*cprimed*aa*aa);
      lp.c(ZERO);
    } else{ 
      lp.rho() = ZERO;
      lp.vx() = ZERO;
      lp.vy() = ZERO;
      lp.p() = ZERO;
      lp.c(ZERO);
      lp.c(i-NUM_FLAME2D_VAR_SANS_SPECIES-1) = ONE;
    } 
    
    // Compute the right conserved eigenvector, and store it in-place
    // 
    if(i == 1){
      rc.rho() = ONE;
      rc.rhovx() = (uprimed-cprimed)/MR2;
      rc.rhovy() = vy();
      rc.E() = h()+HALF*(vsqr()/MR2) - (vx()*cprimed)/MR2;
      rc.rhoc(c());
    } else if(i == 2) {
      rc.rho() = ONE;
      rc.rhovx() = vx();
      rc.rhovy() = vy();
      rc.E() = (h()-Cp()*T()) + HALF*vsqr();
      rc.rhoc(c());
    } else if(i == 3) {
      rc.rho() = ZERO;
      rc.rhovx() = ZERO;
      rc.rhovy() = rho();
      rc.E() = rho()*vy();
      rc.rhoc(ZERO);
    } else if(i == 4) {  
      rc.rho() = ONE;
      rc.rhovx() = (uprimed+cprimed)/MR2;
      rc.rhovy() = vy();
      rc.E() = h()+HALF*(vsqr()/MR2) + (vx()*cprimed)/MR2;
      rc.rhoc(c());
    } else{ 
      rc.rho() = ZERO;
      rc.rhovx() = ZERO;
      rc.rhovy() = ZERO;
      rc.E() = rho()*dihdic[i-NUM_FLAME2D_VAR_SANS_SPECIES-1];
      rc.rhoc(ZERO);
      rc.rhoc(i-NUM_FLAME2D_VAR_SANS_SPECIES-1) = rho();
    } 

    // Compute the wavestrengh
    dv = lp * dWrl;

    // compute the flux dissipation
    tmp = HALF*wavespeeds[i]*dv;
    for ( int k = 1 ; k <= n ; k++ ) Flux_dissipation[k] -= tmp*rc[k];

  } // endfor

}

/****************************************************
 * For CFL calculation
 ****************************************************/
double Flame2D_pState::u_plus_aprecon(const double &u, 
				      const int &flow_type_flag, 
				      const double &deltax) const {
  double Temp(T());
  double aa(a());
  double UR2( Mr2(flow_type_flag,deltax)*aa*aa );
  double alpha( HALF*( ONE - (ONE/(Rtot()*Temp) - ONE/(Cp()*Temp))*UR2) );
  
  // uprime + cprime
  return ( u*(ONE - alpha) + sqrt(alpha*alpha*u*u + UR2) );
}

/****************************************************
 * For eignvalues and eigenvectors return preconditioned
 * velocity and soundspeed ie. u' and c'
 ****************************************************/
void Flame2D_pState::u_a_precon(const double &UR2, double &uprimed, double &cprimed) const{
  
  double Temp(T());
  double alpha( HALF*( ONE - (ONE/(Rtot()*Temp) - ONE/(Cp()*Temp))*UR2) );

  uprimed = vx()*(ONE - alpha);
  cprimed = sqrt(alpha*alpha*vx()*vx() + UR2); 

} 

/****************************************************
 * as defined by E.Turkel (1999)
 ****************************************************/
double Flame2D_pState::Mr2(const int    &flow_type_flag, 
                          const double &deltax) const {
  
  double aa(a());
  double MR2 = min(max((vsqr()/(aa*aa)),Mref*Mref),ONE);
  
  // Need deltax which is based on cell spacing 
  if(flow_type_flag){ 
    MR2 = pow(max(sqrt(MR2*aa*aa), mu()/(rho()*deltax)),2.0)/(aa*aa);
  }

  return (MR2);
}

/************************************************************
 * ORIGINAL LAMINAR ONLY
 * (will be replaced by version commented out below 
 *  when corrected)
 ************************************************************/
void Flame2D_pState::Low_Mach_Number_Preconditioner(::DenseMatrix &P,
						    const int &Viscous_flag, 
						    const double &deltax ) const{  
  // Note: make sure you loaded the dihdic array outside
  double Temp( T() );
  double Rmix( Rtot() );
  double enthalpy( h() );
  double CP( Cp() );
  double aa( a() );
  double theta( ONE/(Mr2(Viscous_flag,deltax)*aa*aa) + ONE/(CP*Temp) );
 
  double phi = ZERO;   
  for(int j=0; j<ns-NSm1; j++){   
    phi += c(j)*dihdic[j];
  }		      

  double alpha( theta*p()/rho() );
  double alpham1( alpha - ONE );
  double Omega( (Rmix - CP)*p()/(rho()*Rmix) );
  double beta( enthalpy - CP*p()/(rho()*Rmix) - phi );
  double V( HALF*vsqr() );

  P(0,0) = (alpha*(beta-V)+V+Rmix*Temp-enthalpy+phi)/Omega;
  P(0,1) = vx()*alpham1/Omega;
  P(0,2) = vy()*alpham1/Omega;
  P(0,3) = -alpham1/Omega;
  P(1,0) = vx()*(beta-V)*alpham1/Omega;
  P(1,1) = vx()*vx()*alpham1/Omega+1.0;
  P(1,2) = vx()*vy()*alpham1/Omega;
  P(1,3) = -vx()*alpham1/Omega;
  P(2,0) = vy()*(beta-V)*alpham1/Omega;
  P(2,1) = vx()*vy()*alpham1/Omega;
  P(2,2) = vy()*vy()*alpham1/Omega+1.0;
  P(2,3) = -vy()*alpham1/Omega;
  P(3,0) = (enthalpy+V)*(beta-V)*alpham1/Omega;
  P(3,1) = vx()*(enthalpy+V)*alpham1/Omega;
  P(3,2) = vy()*(enthalpy+V)*alpham1/Omega;
  P(3,3) = -(alpha*(enthalpy+V)-V-Rmix*Temp-beta-phi)/Omega;

  int NUM_VAR = NUM_FLAME2D_VAR_SANS_SPECIES;

  //Multispecies
  for(int j=0; j<ns-NSm1; j++){  
     
    P(0,j+NUM_VAR) = dihdic[j]*alpham1/Omega;
    P(1,j+NUM_VAR) = vx()*dihdic[j]*alpham1/Omega;
    P(2,j+NUM_VAR) = vy()*dihdic[j]*alpham1/Omega;
    P(3,j+NUM_VAR) = dihdic[j]*(V+enthalpy)*alpham1/Omega;	
    for(int i=0; i<ns-NSm1; i++){ 
      if(i==j){ 
	P(i+NUM_VAR,0) = (c(i))*(beta-V)*alpham1/Omega;
	P(i+NUM_VAR,1) = (c(i))*vx()*alpham1/Omega;
	P(i+NUM_VAR,2) = (c(i))*vy()*alpham1/Omega;
	P(i+NUM_VAR,3) = -(c(i))*alpham1/Omega;
	//diagonal
	P(i+NUM_VAR,j+NUM_VAR) = c(i)*dihdic[j]*alpham1/Omega+ONE;
      }
      else {
	P(i+NUM_VAR,j+NUM_VAR) = c(i)*dihdic[j]*alpham1/Omega;
      }
    }       
  }

//   cout<<"\n Pin with Mref \n"<<Mref<<endl<<P;
  //P.zero(); P.identity(); //SET TO Mref=1.0;


}

/************************************************************/
void Flame2D_pState::Low_Mach_Number_Preconditioner_Inverse(::DenseMatrix &Pinv,	
							    const int &Viscous_flag, 
							    const double &deltax ) const{  
  // Note: make sure you loaded the dihdic array outside
  double Temp( T() );
  double Rmix( Rtot() );
  double enthalpy( h() );
  double CP( Cp() );
  double aa( a() );
  double theta( ONE/(Mr2(Viscous_flag,deltax)*aa*aa) + ONE/(CP*Temp) );

  double phi = ZERO;
  for(int j=0; j<ns-NSm1; j++){ 
    phi += c(j)*dihdic[j];
  }

  double AA( p()*(rho()*Rmix-theta*p()*CP) );
  double BB( Rmix*rho()*(theta*p()-rho()) );
  double EE( HALF*vsqr() - enthalpy + phi );
  double CC( EE + CP*Temp );
  double DD( HALF*vsqr() + enthalpy );
  
  Pinv(0,0) = rho()*Rmix/AA*(theta*p()*EE-rho()*CC+p());
  Pinv(0,1) = -vx()*BB/AA;
  Pinv(0,2) = -vy()*BB/AA;
  Pinv(0,3) = BB/AA;
  Pinv(1,0) = vx()*CC*BB/AA;
  Pinv(1,1) = rho()*Rmix/AA*(p()+rho()*vx()*vx()-theta*p()*(vx()*vx()+CP*Temp));
  Pinv(1,2) = -vx()*vy()*BB/AA;
  Pinv(1,3) = vx()*BB/AA;    
  Pinv(2,0) = vy()*CC*BB/AA;
  Pinv(2,1) = -vx()*vy()*BB/AA;
  Pinv(2,2) = rho()*Rmix/AA*(p()+vy()*vy()*rho()-theta*p()*(vy()*vy()+CP*Temp));
  Pinv(2,3) = vy()*BB/AA;  
  Pinv(3,0) = DD*CC*BB/AA;
  Pinv(3,1) = -vx()*DD*BB/AA;
  Pinv(3,2) = -vy()*DD*BB/AA;
  Pinv(3,3) = rho()*Rmix/AA*(theta*p()*(DD-CP*Temp)-rho()*DD+p());

  Pinv(4,4) = ONE; //fixes so it can work without Turbulence !!!
  Pinv(5,5) = ONE;

  int NUM_VAR = NUM_FLAME2D_VAR_SANS_SPECIES;
  
  //Multispecies
  for(int j=0; j<ns-NSm1; j++){   
    Pinv(0,j+NUM_VAR) = -dihdic[j]*BB/AA;
    Pinv(1,j+NUM_VAR) = -vx()*dihdic[j]*BB/AA;
    Pinv(2,j+NUM_VAR) = -vy()*dihdic[j]*BB/AA;
    Pinv(3,j+NUM_VAR) = -dihdic[j]*BB*DD/AA;
    for(int i=0; i<ns-NSm1; i++){  
      if(i==j){
	Pinv(i+NUM_VAR,0) = (c(i))*CC*BB/AA;
	Pinv(i+NUM_VAR,1) = -(c(i))*vx()*BB/AA;
	Pinv(i+NUM_VAR,2) = -(c(i))*vy()*BB/AA;
	Pinv(i+NUM_VAR,3) = (c(i))*BB/AA;
	//diagonal	
	Pinv(i+NUM_VAR,j+NUM_VAR) = 1.0 - c(i)*dihdic[j]*BB/AA ;
      }
      else {
	Pinv(i+NUM_VAR,j+NUM_VAR) = -c(i)*dihdic[j]*BB/AA;
      } 
    }   
  }  

  //Pinv.zero(); Pinv.identity(); //SET TO Mref=1.0;

}


/////////////////////////////////////////////////////////////////////
/// Source Terms
/////////////////////////////////////////////////////////////////////


/****************************************************
 * Axisymmetric Source Terms
 ****************************************************/  
void Flame2D_pState::Sa_inviscid(Flame2D_State &S, 
				 const Vector2D &X, 
				 const int Axisymmetric, 
				 const double& mult) const {
  //x is radial
  if (Axisymmetric == AXISYMMETRIC_X) {
    S.rho() -= mult*rho()*vx()/X.x;
    S.rhovx() -= mult*rho()*vx()*vx()/X.x; 
    S.rhovy() -= mult*rho()*vx()*vy()/X.x;
    S.E() -= mult*vx()*H()/X.x;
    //species contributions
    for(int i=0; i<ns;i++){
      S.rhoc(i) -= mult*rho()*vx()*c(i)/X.x;         //correct for ns-1 ????
    }

  //y is radial
  } else if (Axisymmetric == AXISYMMETRIC_Y) {
    S.rho() -= mult*rho()*vy()/X.y;
    S.rhovx() -= mult*rho()*vx()*vy()/X.y; 
    S.rhovy() -= mult*rho()*vy()*vy()/X.y;
    S.E() -= mult*vy()*H()/X.y;
    //species contributions
    for(int i=0; i<ns; i++){
      S.rhoc(i) -= mult*rho()*vy()*c(i)/X.y;
    }
  }


}
/***************************************************************
 * Axisymmetric flow source terms (Viscous)                    *  
 ***************************************************************/
void Flame2D_pState::Sa_viscous(Flame2D_State &S, 
				const Flame2D_State &dWdx,
				const Flame2D_State &dWdy,
				const Vector2D &X, 
				const int Axisymmetric, 
				const double& mult) const {
  
  // compute laminar stresses and heat flux vector
  static Vector2D qflux;
  static Tensor2D tau;
  Viscous_Quantities(dWdx, dWdy, Axisymmetric, X, qflux, tau);


  // Determine axisymmetric source terms
  if (Axisymmetric == AXISYMMETRIC_X) {
    S.rhovx() += mult * (tau.xx - tau.zz)/X.x;
    S.rhovy() += mult * tau.xy/X.x;
    S.E() += mult * (-qflux.x + vx()*tau.xx + vy()*tau.xy)/X.x;
    for(int i=0; i<ns;i++){
      S.rhoc(i) += mult * rho()*Diffusion_coef(i)*dWdx.c(i)/X.x;
    }
    
  } else if (Axisymmetric == AXISYMMETRIC_Y) {
    S.rhovx() += mult * tau.xy/X.y;
    S.rhovy() += mult * (tau.xx - tau.zz)/X.y;
    S.E() += mult * (-qflux.y + vx()*tau.xy + vy()*tau.yy)/X.y;
    for(int i=0; i<ns;i++){
      S.rhoc(i) += mult * rho()*Diffusion_coef(i)*dWdy.c(i)/X.y;
    }
    
  } /* endif */
  
}

/****************************************************
 * Source terms associated with finite-rate chemistry
 ****************************************************/  
void Flame2D_pState::Sw(Flame2D_State &S, 
			const double& mult) const {
  Mix.getRates( p(), c(), r );
  for(int i=0; i<ns; i++) S.rhoc(i) += mult*r[i];
}

/****************************************************
 * Source terms associated with gravitational forces
 ****************************************************/
void Flame2D_pState::Sg(Flame2D_State &S, 
			const double& mult) const {
  //Gravity only in the z or axial direction
  //  z|   
  //   |___ r
  //   
  //S.rho() += 0.0;
  //S.rhovx() += 0.0;  //rho*g_r
  S.rhovy() += mult*rho()*gravity_z;
  S.E() += mult*rho()*gravity_z*vy();
}


/////////////////////////////////////////////////////////////////////
/// Boundary Conditions
/////////////////////////////////////////////////////////////////////


/********************************************************
 * Routine: Reflect                                     *
 *                                                      *
 * This function returns the reflected solution state   *
 * in a given direction given the primitive solution    *
 * variables and the unit normal vector in the          *
 * direction of interest.                               *
 *                                                      *
 ********************************************************/
void Flame2D_pState::Reflect(const Flame2D_pState &W,
			     const Vector2D &norm_dir) {

    Copy(W);
 
    /* Determine the direction cosine's for the frame
       rotation. */
    double cos_angle(norm_dir.x), sin_angle(norm_dir.y);


    /* Apply the frame rotation and calculate the primitive
       solution state variables in the local rotated frame
       defined by the unit normal vector. */
    double ur(   W.vx()*cos_angle + W.vy()*sin_angle );
    double vr( - W.vx()*sin_angle + W.vy()*cos_angle );

    /* Reflect the normal velocity in the rotated frame. */
    ur = -ur;

    /* Rotate back to the original Cartesian reference frame. */
    vx() = ur*cos_angle - vr*sin_angle;
    vy() = ur*sin_angle + vr*cos_angle;
  }

/********************************************************
 * Routine: FREE SLIP                                   *
 ********************************************************/
void Flame2D_pState::Free_Slip(const Flame2D_pState &Win,
			       const Flame2D_pState &Wout,
			       const Vector2D &norm_dir,
			       const int &TEMPERATURE_BC_FLAG) { 
  
  /***** MOHAMMED et al *******/
  //Flame2D_pState Temp(Wout);
  //Temp.v = Win.v;
  
  /**** DAY & BELL ********/
  Reflect(Win,norm_dir);
 
  //Fixed Temperature
  if(TEMPERATURE_BC_FLAG == FIXED_TEMPERATURE_WALL){
    setTemperature(Wout.T());
  }
}

/********************************************************
 * Routine: NO SLIP                                     *
 ********************************************************/
void Flame2D_pState::No_Slip(const Flame2D_pState &Win,
			     const Flame2D_pState &Wout,
			     const Vector2D &norm_dir,
			     const int &TEMPERATURE_BC_FLAG) {  
  
  Moving_Wall(Win,Wout,norm_dir,ZERO,TEMPERATURE_BC_FLAG);

}

/********************************************************
 * Routine: Moving_Wall
 ********************************************************/
void Flame2D_pState::Moving_Wall(const Flame2D_pState &Win,
				 const Flame2D_pState &Wout,
				 const Vector2D &norm_dir, 
				 const double &wall_velocity,
				 const int &TEMPERATURE_BC_FLAG) {

  Copy(Win);
  
  /* Determine the direction cosine's for the frame
     rotation. */
  double cos_angle( norm_dir.x ), sin_angle( norm_dir.y );

  /* Apply the frame rotation and calculate the primitive
     solution state variables in the local rotated frame
     defined by the unit normal vector. */
  double ur(   Win.vx()*cos_angle +  Win.vy()*sin_angle );
  double vr( - Win.vx()*sin_angle +  Win.vy()*cos_angle );
  
  /* Use the Roeaveraged value to find the correct tangential velocity */
  ur = -ur;
  vr = -TWO*wall_velocity - vr;

  /* Rotate back to the original Cartesin reference frame. */  
  vx() = ur*cos_angle - vr*sin_angle;
  vy() = ur*sin_angle + vr*cos_angle;
   
 //   /* Fixed Wall Temperature */
//    if(TEMPERATURE_BC_FLAG == FIXED_TEMPERATURE_WALL){
//      Temp.rho = Temp.p/(Temp.Rtot()*Wout.T());
//    }
   
   //ADIABATIC_WALL -> constant extrapolation for Adiabatic */

}


/********************************************************
 * Routine: BC_Flame_Inflow                             *
 *        - 1DFlame Inflow conditions                   *
 *                                                      *
 ********************************************************/
void Flame2D_pState::BC_1DFlame_Inflow(const Flame2D_pState &Wi,      //ICu
				       const Flame2D_pState &Wo,      //Ghost
				       const Flame2D_pState &Woutlet, //ICl
				       const Vector2D &norm_dir){ 
  //fixed Wo
  Copy(Wo);
  vy() = ZERO;

  //Constant Extrapolate Wi 
  //Calculate upstream velocity to balance flame ie. mass flow rate 
  // rho_1*u_1 = rho_2*u_2  
  //relaxation
  //Wnew.v.x = Wi.v.x + 0.5*(Woutlet.rho*Woutlet.v.x/Wo.rho - Wi.v.x);

  //no relaxation
  vx() = Woutlet.rho()*Woutlet.vx()/Wo.rho();

  if( (vx() < 0.1) || (vx() > (Wo.vx()+0.5)) ){
    vx() = Wo.vx();
  }

}
/********************************************************
 * Routine: BC_Flame_Outflow                            *
 *        - 1DFlame outflow conditions                  *
 *                                                      *
 ********************************************************/
void Flame2D_pState::BC_1DFlame_Outflow(const Flame2D_pState &Wi,     //ICu
					const Flame2D_pState &Wo,     //Ghost
					const Flame2D_pState &Winlet, //ICl
					const Vector2D &norm_dir){

  //Constant extrapolate ( zero gradient)
  Copy(Wi);
  vy() = ZERO;

  //Calculate Pressure assuming constant mass flow rate
  //and Wo.p == Winput.p (constant pressure initial condition)
  double sum( Wi.rho()*Wi.vx()*(Wi.vx() - Winlet.vx()) );
  if( sum < ZERO){
    setPressure( Wo.p() );
  } else {
    //no relaxation
    setPressure( Wo.p() - sum );

    // relaxation
    //Wnew.p = Wi.p + 0.5*( (Wo.p-sum) - Wi.p);
  }

 
}

/********************************************************
 * Routine: BC_2DFlame_Inflow                           *
 *        - 2DFlame Inflow conditions                   *
 *                                                      *
 ********************************************************/
void Flame2D_pState::BC_2DFlame_Inflow(const Flame2D_pState &Wi,
				       const Flame2D_pState &Wo,
				       const Vector2D &norm_dir){ 

  //fixed rho, v, p, and species
  Copy(Wo);
  vx() = Wi.vx();

//   Flame2D_pState Wnew(Wi);  
//   Wnew.p = Wo.p;         //fix pressure & V velocity
//   Wnew.v.y = Wo.v.y;
 
}

/********************************************************
 * Routine: BC_2DFlame_Outflow                          *
 *        - 2DFlame outflow conditions                  *
 *                                                      *
 ********************************************************/
void Flame2D_pState::BC_2DFlame_Outflow(const Flame2D_pState &Wi, 
					const Flame2D_pState &Wo,
					const Vector2D &norm_dir){
  //2D Coreflame OUTFLOW hold pressure
  Copy(Wi);  
  if(vy() < ZERO){ 
    vy() = ZERO;
  }

  // set the new mixture state
  setPressure( Wo.p() );

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
void Flame2D_pState::BC_Characteristic_Pressure(const Flame2D_pState &Wi,
						const Flame2D_pState &Wo,
						const Vector2D &norm_dir) {

  //USED FOR BUMP FLOW EXIT
  Flame2D_pState Wi_rotated(Wi), Wo_rotated(Wo);
  double ub_rotated, vb_rotated;
  
  /* Determine the direction cosine's for the frame
     rotation. */
  
  double cos_angle( norm_dir.x ), sin_angle( norm_dir.y );
  
  /* Apply the frame rotation and evaluate interior and 
     imposed boundary solution states in the local rotated 
       frame defined by the unit normal vector. */

  Wi_rotated.vx() =   Wi.vx()*cos_angle + Wi.vy()*sin_angle;
  Wi_rotated.vy() = - Wi.vx()*sin_angle + Wi.vy()*cos_angle;
 
  Wo_rotated.vx() =   Wo.vx()*cos_angle + Wo.vy()*sin_angle;
  Wo_rotated.vy() = - Wo.vx()*sin_angle + Wo.vy()*cos_angle;
 
  /* Determine the Mach number at the interior node. */
  double mi( Wi_rotated.vx()/Wi_rotated.a() );
  
  /* Boundary condition for supersonic outflow. */
  if (mi >= ONE) {
    Copy(Wi);
    
    /* Boundary condition for subsonic outflow. 
       Pressure specified. */
  } else if (mi >= ZERO) {
    Copy(Wi_rotated);
    p() = Wo_rotated.p();
    rho() = Wi_rotated.rho()*pow(p()/Wi_rotated.p(), ONE/Wi_rotated.g());
    setGas();

    double ab( a() );
    ub_rotated = Wi_rotated.vx() + TWO*(Wi_rotated.a()-ab)/(Wi_rotated.g() -ONE);
    vb_rotated = Wi_rotated.vy();
    
    /* Boundary condition for subsonic inflow. 
       Pressure specified. */
  } else if (mi >= -ONE) {
    Copy(Wo_rotated);
    
    double ab( a() );
    ub_rotated = Wi_rotated.vx() + TWO*(Wi_rotated.a()-ab)*(Wo_rotated.g() -ONE);
    vb_rotated = Wo_rotated.vy();
     
    /* Boundary condition for supersonic inflow.  */
  } else {
    Copy(Wo);
  } 

  if( fabs(mi) <=ONE){
    //rotate frame back to cartesian
    vx() = ub_rotated*cos_angle - vb_rotated*sin_angle;
    vy() = ub_rotated*sin_angle + vb_rotated*cos_angle;
  }

}

/////////////////////////////////////////////////////////////////////
/// External Flux Function Functions
/////////////////////////////////////////////////////////////////////


/********************************************************
 * Routine: HartenFixPos (Harten Entropy Fix)           *
 *                                                      *
 * This function returns the positive parts of the      *
 * corrected elemental wave speeds or eigenvalues       *
 * according to the entropy fix of Harten (1983).       *
 *                                                       *
 ********************************************************/
void Flame2D_State :: HartenFix_Pos(const Flame2D_State &lambdas_a,
				    const Flame2D_State &lambdas_l,
				    const Flame2D_State &lambdas_r) {
  rho() = HartenFixPos(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
  vx() = HALF*(lambdas_a[2]+fabs(lambdas_a[2]));
  vy() = HALF*(lambdas_a[3]+fabs(lambdas_a[3]));
  p() = HartenFixPos(lambdas_a[4],lambdas_l[4],lambdas_r[4]);
  
  for( int i=0; i<ns; i++){
    c(i) = HALF*(lambdas_a.c(i)+fabs(lambdas_a.c(i)));
  }
}

/********************************************************
 * Routine: HartenFixNeg (Harten Entropy Fix)           *
 *                                                      *
 * This function returns the negative parts of the      *
 * corrected elemental wave speeds or eigenvalues       *
 * according to the entropy fix of Harten (1983).       *
 *                                                      *
 ********************************************************/
void Flame2D_State :: HartenFix_Neg(const Flame2D_State &lambdas_a,
				    const Flame2D_State &lambdas_l,
				    const Flame2D_State &lambdas_r) {
  rho() = HartenFixNeg(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
  vx() = HALF*(lambdas_a[2]-fabs(lambdas_a[2]));
  vy() = HALF*(lambdas_a[3]-fabs(lambdas_a[3]));
  p() = HartenFixNeg(lambdas_a[4],lambdas_l[4],lambdas_r[4]);
  
  for( int i=0; i<ns; i++){
    c(i) = HALF*(lambdas_a.c(i)-fabs(lambdas_a.c(i)));
  }
}

/********************************************************
 * Routine: HartenFixAbs (Harten Entropy Fix)           *
 *                                                      *
 * This function returns the absolute values of the     *
 * corrected elemental wave speeds or eigenvalues       *
 * according to the entropy fix of Harten (1983).       *
 *                                                      *
 ********************************************************/
void Flame2D_State :: HartenFix_Abs(const Flame2D_State &lambdas_a,
				    const Flame2D_State &lambdas_l,
				    const Flame2D_State &lambdas_r) {
  rho() = HartenFixAbs(lambdas_a[1],lambdas_l[1],lambdas_r[1]);
  vx() = fabs(lambdas_a[2]);
  vy() = fabs(lambdas_a[3]);
  p() = HartenFixAbs(lambdas_a[4],lambdas_l[4],lambdas_r[4]);
  
  for( int i=0; i<ns; i++){
    c(i) = fabs(lambdas_a.c(i));
  }
}


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
void Flame2D_pState :: RoeAverage(const Flame2D_pState &Wl,
				  const Flame2D_pState &Wr) {

    /* Determine the left and right state specific enthalpies
       and square roots of the density. */
    double Hl( Wl.H()/Wl.rho() );
    double Hr( Wr.H()/Wr.rho() );
    double srhol( sqrt(Wl.rho()) );
    double srhor( sqrt(Wr.rho()) );

    /* Determine the appropriate Roe averages. */
    rho() = srhol*srhor;
    vx() = (srhol*Wl.vx()+srhor*Wr.vx())/(srhol+srhor);
    vy() = (srhol*Wl.vy()+srhor*Wr.vy())/(srhol+srhor);
    for(int i=0; i<ns; i++){
      c(i) = (srhol*Wl.c(i) + srhor*Wr.c(i))/(srhol+srhor);
    }
 
    double Ha( (srhol*Hl+srhor*Hr)/(srhol+srhor) );
    double ha( Ha - HALF*(sqr(vx())+sqr(vy())) );

    //double TEMP = Temp.T(ha);
    // set pressure
    setEnthalpy(ha);
  
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
void Flame2D_State :: FluxHLLE_x(const Flame2D_pState &Wl,
				 const Flame2D_pState &Wr) {

    double wavespeed_l, wavespeed_r;
    static Flame2D_pState Wa;
    static Flame2D_State lambdas_l, lambdas_r, lambdas_a;

    /* Evaluate the Roe-average primitive solution state. */   
    Wa.RoeAverage(Wl, Wr);

    /* Evaluate the left, right, and average state eigenvalues. */     
    Wl.lambda_x(lambdas_l);
    Wr.lambda_x(lambdas_r);
    Wa.lambda_x(lambdas_a);
    
    /* Determine the intermediate state flux. */
    wavespeed_l = min(lambdas_l[1], lambdas_a[1]);

    //wavespeed_r = max(lambdas_r[NUM_VAR_FLAME2D],
    //                  lambdas_a[NUM_VAR_FLAME2D]);
    wavespeed_r = max(lambdas_r[4],lambdas_a[4]); //MARKTHIS XINFENG
    //   wavespeed_r = max(lambdas_r[NUM_VAR_FLAME2D-lambdas_r.ns],
    //                       lambdas_a[NUM_VAR_FLAME2D-lambdas_a.ns]); 
    wavespeed_l = min(wavespeed_l, ZERO);
    wavespeed_r = max(wavespeed_r, ZERO);
  

    if (wavespeed_l >= ZERO) {
      Wl.Fx(*this);
    } else if (wavespeed_r <= ZERO) {
      Wr.Fx(*this);
    } else { 
      static Flame2D_State dUrl, Fl, Fr;
      dUrl.Delta( Wr.U(), Wl.U() ); // Evaluate the jumps in the conserved solution states.
      Wr.Fx(Fr);
      Wl.Fx(Fl);
      double a(wavespeed_l*wavespeed_r), b(wavespeed_r-wavespeed_l);
      for (int i=0; i<n; i++) {
	x[i] = ( (wavespeed_r*Fl.x[i]-wavespeed_l*Fr.x[i]) + a*dUrl.x[i] ) / b;
      }
    } /* endif */

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
void Flame2D_State :: FluxHLLE_n(const Flame2D_pState &Wl,
				 const Flame2D_pState &Wr,
				 const Vector2D &norm_dir) {

  static Flame2D_pState Wl_rotated, Wr_rotated;
  Wl_rotated.Copy(Wl), Wr_rotated.Copy(Wr);
  
  /* Determine the direction cosine's for the frame
     rotation. */
  
  double cos_angle( norm_dir.x), sin_angle( norm_dir.y );
  
  /* Apply the frame rotation and evaluate left and right
     solution states in the local rotated frame defined
     by the unit normal vector. */
  //    Wl_rotated.Copy(Wl);
  //     Wr_rotated.Copy(Wr);
  
  
  Wl_rotated.setVelocity( Wl.vx()*cos_angle + Wl.vy()*sin_angle,
			  - Wl.vx()*sin_angle + Wl.vy()*cos_angle );
  
  Wr_rotated.setVelocity( Wr.vx()*cos_angle + Wr.vy()*sin_angle,
			  - Wr.vx()*sin_angle + Wr.vy()*cos_angle );
  
  /* Evaluate the intermediate state solution 
     flux in the rotated frame. */
  
  FluxHLLE_x(Wl_rotated, Wr_rotated);
  double ur(rhovx()), vr(rhovy());
  
  /* Rotate back to the original Cartesian reference
     frame and return the solution flux. */
  rhovx() = ur*cos_angle - vr*sin_angle;
  rhovy() = ur*sin_angle + vr*cos_angle;
  
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
void Flame2D_State :: FluxLinde(const Flame2D_pState &Wl,
				const Flame2D_pState &Wr) {

  double wavespeed_l, wavespeed_r;
  static Flame2D_pState Wa;
  static Flame2D_State lambdas_l, lambdas_r, lambdas_a;
  
  /* Evaluate the Roe-average primitive solution state. */   
  Wa.RoeAverage(Wl, Wr);
  
  /* Evaluate the left, right, and average state eigenvalues. */
  Wl.lambda_x(lambdas_l);
  Wr.lambda_x(lambdas_r);
  Wa.lambda_x(lambdas_a);
  
  /* Determine the intermediate state flux. */
  wavespeed_l = min(lambdas_l[1], lambdas_a[1]);
  //   wavespeed_r = max(lambdas_r[NUM_VAR_FLAME2D-lambdas_r.ns],
  //                       lambdas_a[NUM_VAR_FLAME2D-lambdas_a.ns]);
  wavespeed_r = max(lambdas_r[4],lambdas_a[4]);
  
  if (wavespeed_l >= ZERO) {
    Wl.Fx(*this);

  } else if (wavespeed_r <= ZERO) {
    Wr.Fx(*this);

  } else {

    static Flame2D_State Fr, Fl;
    static Flame2D_State dFrl, dUrl;
    const Flame2D_pState& Waa = Wa; // need to force const member functions
    double alpha;
    Wr.Fx(Fr);  Wl.Fx(Fl);

    dUrl.Delta( Wr.U(), Wl.U() );
    dFrl.Delta( Fr, Fl );

    double wavespeed_m( Waa.vx() );
    double rhoa( Waa.rho() );
    double ca( Waa.a() );
    double dU( fabs(dUrl.rho())/rhoa+
	       fabs(dUrl.rhovx())/(rhoa*ca)+
	       fabs(dUrl.rhovy())/(rhoa*ca)+
	       fabs(dUrl.E())/(rhoa*ca*ca) );
    if (dU <= TOLER) {
      alpha = ZERO;
    } else {
      dU = ONE/dU;
      static Flame2D_State dFwave;
      dFwave.Delta( dFrl, wavespeed_m*dUrl );
      alpha = ONE - (fabs(dFwave.rho())/(rhoa*ca)+
		     fabs(dFwave.rhovx())/(rhoa*ca*ca)+
		     fabs(dFwave.rhovy())/(rhoa*ca*ca)+
		     fabs(dFwave.E())/(rhoa*ca*ca*ca))*dU;
      alpha = max(ZERO, alpha);
      
    } /* endif */
    
    // compute the flux
    double a(wavespeed_l*wavespeed_r), 
           b(wavespeed_r-wavespeed_l),
           c(ONE-(ONE-max(wavespeed_m/wavespeed_r, wavespeed_m/wavespeed_l))*alpha);
    for (int i=0; i<n; i++) {
      x[i] =( (wavespeed_r*Fl.x[i]-wavespeed_l*Fr.x[i]) + a*c*dUrl.x[i]) / b;
    }

  } /* endif */
  
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
void Flame2D_State :: FluxLinde_n(const Flame2D_pState &Wl,
				  const Flame2D_pState &Wr,
				  const Vector2D &norm_dir) {

  static Flame2D_pState Wl_rotated, Wr_rotated;
  Wl_rotated.Copy(Wl);  Wr_rotated.Copy(Wr);
  
  /* Determine the direction cosine's for the frame
     rotation. */
  
  double cos_angle( norm_dir.x ); 
  double sin_angle( norm_dir.y );
  
  /* Apply the frame rotation and evaluate left and right
     solution states in the local rotated frame defined
     by the unit normal vector. */
  //     Wl_rotated.Copy(Wl);
  //     Wr_rotated.Copy(Wr);
  
  Wl_rotated.setVelocity( Wl.vx()*cos_angle + Wl.vy()*sin_angle,
			  - Wl.vx()*sin_angle + Wl.vy()*cos_angle );
  
  Wr_rotated.setVelocity( Wr.vx()*cos_angle + Wr.vy()*sin_angle,
			  - Wr.vx()*sin_angle + Wr.vy()*cos_angle );
  
  /* Evaluate the intermediate state solution 
     flux in the rotated frame. */ 
  
  FluxLinde(Wl_rotated, Wr_rotated);
  double ur(rhovx()), vr(rhovy());
  
  /* Rotate back to the original Cartesian reference
     frame and return the solution flux. */
  
  rhovx() = ur*cos_angle - vr*sin_angle;
  rhovy() = ur*sin_angle + vr*cos_angle;

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
void Flame2D_State :: FluxRoe_x(const Flame2D_pState &Wl,
				const Flame2D_pState &Wr,
				const int &Preconditioning,
				const int &flow_type_flag,
				const double &deltax) {
  
  static Flame2D_pState Wa;
  static Flame2D_State dWrl, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
  
  /* Evaluate the Roe-average primitive solution state. */    
  Wa.RoeAverage(Wl, Wr);
  const Flame2D_pState& Waa = Wa;

  /* Evaluate the jumps in the primitive solution states. */
  dWrl.Delta(Wr,Wl);
  
  // load the static dihdic array first before computing 
  // eigenvectors
  Waa.load_dihdic();

  /* Evaluate the left, right, and average state eigenvalues. */
  if(!Preconditioning){
    Wl.lambda_x(lambdas_l);
    Wr.lambda_x(lambdas_r);
    Waa.lambda_x(lambdas_a);

    /* Determine the intermediate state flux. */
    if (Waa.vx() >= ZERO) {
      Wl.Fx(*this);   
      wavespeeds.HartenFix_Neg(lambdas_a,
			       lambdas_l,
			       lambdas_r);
      
      for (int i=0 ; i < n-NSm1; i++) {
	if (wavespeeds[i+1] < ZERO) {
	  *this += wavespeeds[i+1]*((Waa.lp_x(i+1)*dWrl)*Waa.rc_x(i+1));
	}
      } 
    } else {
      Wr.Fx(*this);
      wavespeeds.HartenFix_Pos(lambdas_a,
			       lambdas_l,
			       lambdas_r);
      for (int i=0; i < n-NSm1; i++) {
	if (wavespeeds[i+1] > ZERO) {
	  *this -= wavespeeds[i+1]*((Waa.lp_x(i+1)*dWrl)*Waa.rc_x(i+1));
	}
      } 
    } 
    
    /******* LOW MACH NUMBER PRECONDITIONING ********************/
    /* Evaluate the left, right, and average state eigenvalues. */
  } else if(Preconditioning){

    //calculating Mr^2 and passing to save computation,
    //not conceptually nice but saves from recalculating
    double MR2a( Waa.Mr2(flow_type_flag,deltax) );
    Wl.lambda_preconditioned_x(lambdas_l, Wl.Mr2(flow_type_flag,deltax)); 
    Wr.lambda_preconditioned_x(lambdas_r, Wr.Mr2(flow_type_flag,deltax));
    Waa.lambda_preconditioned_x(lambdas_a, MR2a);
    
    // Evaluate the jumps in the primitive solution states.
    wavespeeds.HartenFix_Abs(lambdas_a,
			     lambdas_l,
			     lambdas_r);
     
    // Evaluate the low-Mach-number local preconditioner for the Roe-averaged state.
    static ::DenseMatrix P(n-NSm1,n-NSm1);
    P.zero();  //RESET MATRIX TO ZERO!!!
    Waa.Low_Mach_Number_Preconditioner(P,flow_type_flag,deltax);
    
    // Determine the intermediate state flux.
    (*this).Vacuum();
    Wl.addFx(*this, HALF);
    Wr.addFx(*this, HALF);
    
    // compute dissipation flux
    static Flame2D_State Flux_dissipation;
    Waa.Flux_Dissipation_precon( MR2a, dWrl, wavespeeds, Flux_dissipation);
    
    // Add preconditioned upwind dissipation flux.
    for ( int i = 0 ; i < n-NSm1 ; i++ ) {
      for ( int j = 0 ; j < n-NSm1 ; j++ ) {
	(*this)[i+1] += P(i,j)*Flux_dissipation[j+1];
      } 
    } 
     
  } else {
    cerr<<"\n Not a valid Preconditioner "<<Preconditioning;
    exit(1);
  }  
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
void Flame2D_State :: FluxRoe_n(const Flame2D_pState &Wl,
				const Flame2D_pState &Wr,
				const Vector2D &norm_dir,
				const int &Preconditioning,
				const int &flow_type_flag,
				const double &delta_n ) {

  static Flame2D_pState Wl_rotated, Wr_rotated;
  Wl_rotated.Copy(Wl);  Wr_rotated.Copy(Wr);

  /* Determine the direction cosine's for the frame rotation. */

  double cos_angle( norm_dir.x ), sin_angle( norm_dir.y );

  /* Apply the frame rotation and evaluate left and right
     solution states in the local rotated frame defined
     by the unit normal vector. */
  
  Wl_rotated.setVelocity( Wl.vx()*cos_angle +  Wl.vy()*sin_angle,
			  - Wl.vx()*sin_angle + Wl.vy()*cos_angle );
  
  Wr_rotated.setVelocity( Wr.vx()*cos_angle + Wr.vy()*sin_angle,
			  - Wr.vx()*sin_angle + Wr.vy()*cos_angle );
  
  /* Evaluate the intermediate state solution 
     flux in the rotated frame. */
  
  FluxRoe_x(Wl_rotated, Wr_rotated, 
	    Preconditioning,flow_type_flag,delta_n);
  double ur(rhovx()), vr(rhovy());
  
  /* Rotate back to the original Cartesian reference
     frame and return the solution flux. */
  
  rhovx() = ur*cos_angle - vr*sin_angle;
  rhovy() = ur*sin_angle + vr*cos_angle;

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
void Flame2D_State::FluxAUSMplus_up(const Flame2D_pState &Wl,
				    const Flame2D_pState &Wr) {

  double beta(0.125), sigma(1.0), Kp(0.25), Ku(0.5)/*0.75*/;
  //double al, ar, atilde_l, atilde_r;

  // Determine the intermediate state sound speed and density:
  // al = sqrt(Wl.H()/Wl.rho)*sqrt(TWO*(Wl.g() - ONE)/(Wl.g() + ONE));
  // ar = sqrt(Wr.H()/Wr.rho)*sqrt(TWO*(Wr.g() - ONE)/(Wr.g() + ONE));
  // atilde_l = sqr(al)/max(al, Wl.v.x);    
  // atilde_r = sqr(ar)/max(ar, -Wr.v.x);
  // ahalf = min(atilde_l, atilde_r);
  double ahalf( HALF*(Wl.a() + Wr.a()) );
  double rhohalf( HALF*(Wl.rho() + Wr.rho()) ); 
  

  // Determine the left and right state Mach numbers based on the
  // intermediate state sound speed:
  double Ml( Wl.vx()/ahalf );
  double Mr( Wr.vx()/ahalf );

  // Determine the reference Mach number, scaling function and coefficient
  double M2_bar( (Wl.vx()*Wl.vx() + Wr.vx()*Wr.vx())/(TWO*ahalf*ahalf) );
  double M2_ref( min(ONE, max(M2_bar, Wl.Mref*Wl.Mref)) );
  if (M2_ref > ONE || M2_ref < 0.0) cout << "\nM2_ref out of range";
  //fa = sqrt(M2_ref)*(TWO - sqrt(M2_ref));
  double fa( sqrt(sqr(ONE - M2_ref)*M2_bar + FOUR*M2_ref)/(ONE + M2_ref) );
  if (fa > ONE || fa <= ZERO) cout << "\nfa out of range";
  double alpha( (3.0/16.0)*(-4.0 + 5.0*fa*fa) );
  if (alpha < (-3.0/4.0)  ||  alpha > (3.0/16.0)) cout << "\nalpha out of range";


  // Determine the left state split Mach number:
  double Mplus, pplus;
  if (fabs(Ml) >= ONE) {
    Mplus = Mplus_1(Ml);
    pplus = Mplus_1(Ml)/Ml;
  } else {
    Mplus = Mplus_2(Ml) * (1.0 - 16.0*beta*Mminus_2(Ml));
    pplus = Mplus_2(Ml) * ((2.0 - Ml) - 16.0*alpha*Ml*Mminus_2(Ml));
  }


  // Determine the right state split Mach number:
  double Mminus, pminus;
  if (fabs(Mr) >= ONE) {
    Mminus = Mminus_1(Mr);
    pminus = Mminus_1(Mr)/Mr;        
  } else {
    Mminus = Mminus_2(Mr) * (1.0 + 16.0*beta*Mplus_2(Mr));
    pminus = Mminus_2(Mr) * ((-2.0 - Mr) + 16.0*alpha*Mr*Mplus_2(Mr));
  } 


  // Determine the intermediate state Mach number, pressure and mass flux:
  double Mhalf( Mplus + Mminus - 
		(Kp/fa)*max((ONE - sigma*M2_bar), ZERO)*(Wr.p() - Wl.p())/(rhohalf*ahalf*ahalf) );

  double phalf( pplus*Wl.p() + pminus*Wr.p() - 
		Ku*pplus*pminus*TWO*rhohalf*(fa*ahalf)*(Wr.vx() - Wl.vx()) );

  double mass_flux_half( (Mhalf > ZERO) ? ahalf*Mhalf*Wl.rho() : ahalf*Mhalf*Wr.rho() );


  // Determine the intermediate state convective solution flux:
  if (mass_flux_half  > ZERO) {
    rho() = ONE;
    rhovx() = Wl.vx(); 
    rhovy() = Wl.vy(); 
    E() = Wl.H()/Wl.rho();
    for(int i=0; i<ns; ++i){
      rhoc(i) = Wl.c(i);
    }
    
  } else {
    rho() = ONE;
    rhovx() = Wr.vx(); 
    rhovy() = Wr.vy(); 
    E() = Wr.H()/Wr.rho();
    for(int i=0; i<ns; ++i){
      rhoc(i) = Wr.c(i);
    }

  } //end if

  // scale the flux
  for(int i=0; i<n; ++i) x[i] *= mass_flux_half;

  // Add the pressure contribution to the intermediate state solution flux:
  rhovx() += phalf;

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
void Flame2D_State::FluxAUSMplus_up_n(const Flame2D_pState &Wl,
				      const Flame2D_pState &Wr,
				      const Vector2D &norm_dir) {


  static Flame2D_pState Wl_rotated, Wr_rotated;
  Wl_rotated.Copy(Wl);  Wr_rotated.Copy(Wr);

  /* Determine the direction cosine's for the frame rotation. */

  double cos_angle( norm_dir.x ), sin_angle( norm_dir.y );

  /* Apply the frame rotation and evaluate left and right
     solution states in the local rotated frame defined
     by the unit normal vector. */
  
  Wl_rotated.setVelocity( Wl.vx()*cos_angle +  Wl.vy()*sin_angle,
			  - Wl.vx()*sin_angle + Wl.vy()*cos_angle );
  
  Wr_rotated.setVelocity( Wr.vx()*cos_angle + Wr.vy()*sin_angle,
			  - Wr.vx()*sin_angle + Wr.vy()*cos_angle );
  
  /* Evaluate the intermediate state solution 
     flux in the rotated frame. */
  
  FluxAUSMplus_up(Wl_rotated, Wr_rotated);
  double ur(rhovx()), vr(rhovy());
  
  /* Rotate back to the original Cartesian reference
     frame and return the solution flux. */
  
  rhovx() = ur*cos_angle - vr*sin_angle;
  rhovy() = ur*sin_angle + vr*cos_angle;

}


/////////////////////////////////////////////////////////////////////
/// Viscous Reconstruction Functions
/////////////////////////////////////////////////////////////////////

/**********************************************************************
 * Routine: ViscousFluxArithmetic_n                                   *
 *                                                                    *
 * This function returns the intermediate state solution viscous flux *
 * calculated by the arithmetic mean of the cell-centred flux terms   *
 * of the neighbouring cells.                                         *
 *                                                                    *
 **********************************************************************/
void Flame2D_State::Viscous_FluxArithmetic_n(const Flame2D_pState &Wl,
					     const Flame2D_State &dWdx_l,
					     const Flame2D_State &dWdy_l,
					     const Vector2D &qflux_l,
					     const Tensor2D &tau_l,
					     const Flame2D_pState &Wr,
					     const Flame2D_State &dWdx_r,
					     const Flame2D_State &dWdy_r,
					     const Vector2D &qflux_r,
					     const Tensor2D &tau_r,
					     const Vector2D &norm_dir,
					     const double &mult) {
  // the cosines
  static const Vector2D i(1,0), j(0,1);
  double nx( dot(i,norm_dir) ), ny( dot(j,norm_dir) );

  // Fx = 0.5 * ( Wl.Viscous_Flux_x + Wr.Viscous_Flux_x)
  // Fy = 0.5 * ( Wl.Viscous_Flux_y + Wr.Viscous_Flux_y)
  // Fn = (Fx*dot(i,norm_dir) + Fy*dot(j,norm_dir))

  // Fx Constribution
  if (fabs(nx)>TOLER) {
    Wl.Viscous_Flux_x(dWdx_l, qflux_l, tau_l, *this, mult*HALF*nx);
    Wr.Viscous_Flux_x(dWdx_r, qflux_r, tau_r, *this, mult*HALF*nx);
  }
  // Fy contribution
  if (fabs(ny)>TOLER) {
    Wl.Viscous_Flux_y(dWdy_l, qflux_l, tau_l, *this, mult*HALF*ny);
    Wr.Viscous_Flux_y(dWdy_r, qflux_r, tau_r, *this, mult*HALF*ny);
  }

}

/**********************************************************************
 * Routine: Viscous_Flux_n                                            *
 *                                                                    *
 * This function returns the intermediate state solution viscous flux *
 * given the primitive variable solution state and the gradients of   *
 * the primitive variables.                                           *
 *                                                                    *
 **********************************************************************/
void Flame2D_State::Viscous_Flux_n(const Flame2D_pState &W,
				   const Flame2D_pState &dWdx,
				   const Flame2D_pState &dWdy,
				   const int Axisymmetric,
				   const Vector2D X,
				   const Vector2D &norm_dir,
				   const double &mult) {

  static Vector2D qflux;
  static Tensor2D tau;

  // the cosines
  static const Vector2D i(1,0), j(0,1);
  double nx( dot(i,norm_dir) ), ny( dot(j,norm_dir) );

  // compute the stress tensor and heat flux vector
  W.Viscous_Quantities( dWdx, dWdy, Axisymmetric, X, qflux, tau );
  
  // Fn = (Fx*dot(i,norm_dir) + Fy*dot(j,norm_dir))

  // Fx Constribution
  if (fabs(nx)>TOLER) {
    W.Viscous_Flux_x(dWdx, qflux, tau, *this, mult*nx);
  }
  // Fy contribution
  if (fabs(ny)>TOLER) {
    W.Viscous_Flux_y(dWdy, qflux, tau, *this, mult*ny);
  }

}

/**********************************************************************
 * Routine: ViscousFluxHybrid_n                                       *
 *                                                                    *
 * This function returns the intermediate state solution viscous flux *
 * calculated by the arithmetic mean of the cell-centred flux terms   *
 * of the neighbouring cells.                                         *
 * Also returns the face gradients in dWd[xy] (passed by reference)   *
 *                                                                    *
 **********************************************************************/
void Flame2D_State::Viscous_FluxHybrid_n(const Flame2D_pState &W,
					Flame2D_State &dWdx, 
					Flame2D_State &dWdy,
					const Vector2D &X,
					const Flame2D_pState &Wl,
					const Flame2D_State &dWdx_l,
					const Flame2D_State &dWdy_l,
					const Vector2D &Xl,
					const Flame2D_pState &Wr,
					const Flame2D_State &dWdx_r,
					const Flame2D_State &dWdy_r,
					const Vector2D &Xr,
					const int &Axisymmetric,
					const Vector2D &norm_dir,
					const double &mult) {

  static Vector2D qflux;
  static Tensor2D tau;

  // the cosines
  static const Vector2D i(1,0), j(0,1);
  double nx( dot(i,norm_dir) ), ny( dot(j,norm_dir) );

  // compute the distances
  Vector2D dX( Xr-Xl );
  double ds( dX.abs() );
  dX /= ds;

  // the dot product
  double dotp( dot(norm_dir,dX) );

  // Compute the Cartesian components of the intermediate state
  // solution viscous flux.
  double dWdx_ave, dWdy_ave, dWds;
  for (int k=0; k<n; k++) {
    dWdx_ave = HALF*(dWdx_l.x[k] + dWdx_r.x[k]);
    dWdy_ave = HALF*(dWdy_l.x[k] + dWdy_r.x[k]);
    dWds = (Wr.x[k]-Wl.x[k])/ds;
    dWdx.x[k] = dWdx_ave + (dWds - dWdx_ave*dX.x)*norm_dir.x/dotp;
    dWdy.x[k] = dWdy_ave + (dWds - dWdy_ave*dX.y)*norm_dir.y/dotp;
  }

  // compute the stress tensor and heat flux vector
  W.Viscous_Quantities( dWdx, dWdy, Axisymmetric, X, qflux, tau );
  
  // Fn = (Fx*dot(i,norm_dir) + Fy*dot(j,norm_dir))
  // Fx Constribution
  if (fabs(nx)>TOLER) {
    W.Viscous_Flux_x(dWdx, qflux, tau, *this, mult*nx);
  }
  // Fy contribution
  if (fabs(ny)>TOLER) {
    W.Viscous_Flux_y(dWdy, qflux, tau, *this, mult*ny);
  }

}


/////////////////////////////////////////////////////////////////////
/// Exact Test Case Solution Functions
/////////////////////////////////////////////////////////////////////

/**********************************************************************
 * Routine: RinglebFlow                                               *
 *                                                                    *
 * This function returns the exact solution to Ringleb's flow for the *
 * location X.                                                        *
 *                                                                    *
 **********************************************************************/
void Flame2D_pState::RinglebFlow(const Vector2D X) {

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
void Flame2D_pState::ViscousChannelFlow(const Vector2D X,
					const double Vwall,
					const double dp) {
  vx()  = ( (HALF/mu())*(-dp/0.2)*
	    (pow(X.y,TWO) - (0.001*0.001/4.0)) + 
	    Vwall*(X.y/0.001 + 0.5) );
  vy()  = ZERO; 
  rho() = DENSITY_STDATM;
  p()   = PRESSURE_STDATM + dp*(ONE - X.x/0.20);
  setGas();
}

/**********************************************************************
 * Routine: FlatPlate                                                 *
 *                                                                    *
 * This function returns the exact solution for the flow over a flat  *
 * (adiabatic) plate (Blasius solution) at a given the position and   *
 * the freestream flow velocity.                                      *
 *                                                                    *
 **********************************************************************/
void Flame2D_pState::FlatPlate(const Vector2D X,
			       double &eta,
			       double &f,
			       double &fp,
			       double &fpp) {
  
  // initialize
  double fo = ZERO; 
  double sign = ONE;
  double dn = 0.000005;
  eta = ZERO;
  f = ZERO; fp = ZERO; fpp = 0.33206;

  // Return upstream conditions before flat plate.
  if (X.x < ZERO) return;

  // Return upstream conditions with zero velocity at the leading edge
  // of the plate.
  if (X.x < TOLER) {
    vx() = 0.0;
    vy() = 0.0;
    return;
  }

  // Determine the dimensionless similarity coordinate, eta:
  eta = X.y*sqrt(vx()/(X.x*mu()/rho()));

  // If eta is greater than 8.4, for the sake of expediency, use linear
  // extrapolation to determine the value of f (fp = ONE and fpp = ZERO)
  // given the tabulated value at 8.4 (note, the analytic solution is 
  // linear in this region).
  if (eta > 8.4) {
    fp = ONE; fpp = ZERO; f = 6.67923 + fp*(eta - 8.4);
    vx() = fp*vx();
    vy() = HALF*(eta*fp-f);
    return;
  }

  // Compute the Blasius solution using a fourth-order Runge-Kutta method.
  double k1, k2, k3, k4;
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
  vx() = fp*vx();
  //W.v.y = HALF*(eta*fp-f);

  //if (eta <= 8.8) cout << eta << " " << f << " " << fp << " " << fpp << endl;
  return;
}


/**********************************************************************
 * Routine: WallShearStress                                           *
 *                                                                    *
 * This routine computes and returns the shear stress at a wall.      *
 *                                                                    *
 **********************************************************************/
double Flame2D_pState::WallShearStress(const Vector2D &X1,
				       const Vector2D &X2,
				       const Vector2D &X3,
				       const Vector2D &norm_dir) const {

  double l21, l32, l13, A, dWdn;
  Vector2D n21, n32, n13;
  Flame2D_State W2, W3, W_face, dWdx, dWdy;

  // Initialze W2 and W3.
  W2.Vacuum(); W2.rho() = rho(); W2.p() = p();
  W3.Vacuum(); W3.rho() = rho(); W3.p() = p();
  for(int i=0; i<ns; i++){
    W2.c(i) = c(i);
    W3.c(i) = c(i);
  }

  // Determine the lengths and normals of te faces and the 
  // areas of the regions of Green-Gauss integration.
  l21 = abs(X2-X1); n21 = Vector2D((X2.y-X1.y),-(X2.x-X1.x))/l21;
  l32 = abs(X3-X2); n32 = Vector2D((X3.y-X2.y),-(X3.x-X2.x))/l32;
  l13 = abs(X1-X3); n13 = Vector2D((X1.y-X3.y),-(X1.x-X3.x))/l13;
  A = HALF*((X2-X1)^(X3-X1));
  // Compute Green-Gauss integration on left triangle.
  W_face = HALF*(W2+*this)*l21;
  dWdx = W_face*n21.x;
  dWdy = W_face*n21.y;
  W_face = HALF*(W3+W2)*l32;
  dWdx += W_face*n32.x;
  dWdy += W_face*n32.y;
  W_face = HALF*(*this+W3)*l13;
  dWdx += W_face*n13.x;
  dWdy += W_face*n13.y;
  dWdx = dWdx/A;
  dWdy = dWdy/A;
  // Determine the normal gradient.
  dWdn = dWdy.vx();//dot(Vector2D(dWdx.v.x,dWdy.v.y),norm_dir);

  // Return the wall shear stress.
  return mu()*dWdn;

}
