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
int Flame2D_State :: ns_eqn = 0;
#endif
bool    Flame2D_State  :: reacting = false;
double  Flame2D_State  :: Mref = 0.5;
double  Flame2D_State  :: gravity_z = -9.81;
double* Flame2D_State  :: y = NULL;
double* Flame2D_pState :: r = NULL;
double* Flame2D_pState :: h_i = NULL;
double* Flame2D_pState :: cp_i = NULL;

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
  ns_eqn = ns - NSm1;
#endif

  // determine if this is a reacting case
  if (Mixture::nReactions()>0) reacting = true;
  else reacting = false;

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

void Flame2D_pState::Fx(Flame2D_State &FluxX, const double& mult) const {
  FluxX.rho() = mult * rho()*vx();
  FluxX.rhovx() = mult * (rho()*sqr(vx()) + p());
  FluxX.rhovy() = mult * rho()*vx()*vy();
  FluxX.E() = mult * vx()*H();
  
  //multispecies transport
  for(int i=0; i<ns; i++){
    FluxX.rhoc(i) = mult * rho()*vx()*c(i);
  }
}

void Flame2D_pState::addFx(Flame2D_State &FluxX, const double& mult) const {
  FluxX.rho() += mult * rho()*vx();
  FluxX.rhovx() += mult * (rho()*sqr(vx()) + p());
  FluxX.rhovy() += mult * rho()*vx()*vy();
  FluxX.E() += mult * vx()*H();
  
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

/********************************************************
 * Chem2D_pState::dFIdU -- Invisicd Flux Jacobian       *
 ********************************************************/
void Flame2D_pState::dFIdU(DenseMatrix &dFdU) {
  
  //cout<<"\n USING DFIDU \n";
  
  // Note: make sure you loaded the dihdic array and phi outside
  updateDihdDic();
  const double phi( Phi() );
  const double Temp( T() );
  const double Rt( Rtot() );
  const double C_p( Cp() );
  const double ht( h() );
  const double A( C_p/Rt );
  const double denominator( A - ONE );
  const double vy2( vy()*vy() );
  const double vx2( vx()*vx() );
  const double vsqrd( vsqr() );
  const double vx_vy( vx()*vy() );

  dFdU(0,1) = ONE;
  dFdU(1,0) = ( A*( - vx2 ) + HALF*(THREE*vx2 + vy2) - 
		ht + C_p*Temp + phi ) / denominator;
  dFdU(1,1) = vx()*(TWO*A-THREE)/denominator; 
  dFdU(1,2) = - vy()/denominator;
  dFdU(1,3) = ONE/denominator;
  dFdU(2,0) = - vx_vy;
  dFdU(2,1) = vy();
  dFdU(2,2) = vx();
  dFdU(3,0) = vx()*( vsqrd + C_p*Temp - 
		     A*( HALF*vsqrd + ht ) + phi ) / denominator;
  dFdU(3,1) = ht + HALF*vsqrd - vx2/denominator;
  dFdU(3,2) = - vx_vy/denominator;
  dFdU(3,3) = vx()*A/denominator;
  //Species
  const int NUM_VAR = NUM_FLAME2D_VAR_SANS_SPECIES;
  for(int i = 0; i<ns_eqn; i++){ 
    dFdU(1,NUM_VAR+i) = - dihdic(i)/denominator; 
    dFdU(3,NUM_VAR+i) = - vx()*dihdic(i)/denominator;    
    dFdU(NUM_VAR+i, 0) = - c(i)*vx() ;
    dFdU(NUM_VAR+i, 1) = c(i);
    dFdU(NUM_VAR+i,NUM_VAR+i) = vx();        
  }
 
}

/********************************************************
 * Chem2D_pState::dFIdU -- Invisicd Flux Jacobian       *
 ********************************************************/
// Finite differnece check of dFxdU
void Flame2D_pState::dFIdU_FD(DenseMatrix &dFdU) const {

  // declares
  const int NUM_VAR = dFdU.get_n();
  static Flame2D_State UU;
  static Flame2D_State A,C;
  static Flame2D_pState B,D;
  const double perturb = 5e-6;
  double a;

  // set UU
  getU(UU);

  // finitie differences
  for(int jcol=0; jcol<NUM_VAR; jcol++){    
    A = UU;  C = UU;
    if( jcol <NUM_FLAME2D_VAR_SANS_SPECIES) {
      A[jcol+1] += perturb*max(ONE,UU[jcol+1]); 
      C[jcol+1] -= perturb*max(ONE,UU[jcol+1]);
    } else {                                       //enforce sum(ci) = 1;
      a = perturb*max(ONE,UU[jcol+1]); 
      A[jcol+1] += a;
      A[NUM_VAR+1] -= a;      
      C[jcol+1] -= a;
      C[NUM_VAR+1] += a;
    }   
    B.setU(A);   D.setU(C);    
    B.Fx(A);     D.Fx(C);
    for(int irow=0; irow<NUM_VAR; irow++){
      dFdU(irow,jcol) = ( A[irow+1] - C[irow+1])/(TWO*perturb*max(ONE, UU[jcol+1]));      
    }
  } 

}


/********************************************************
 * Chem2D_pState::dFIdW -- Invisicd Flux Jacobian       *
 ********************************************************/
void Flame2D_pState::dFIdW(DenseMatrix &dFdW, 
			   const double& mult) {
    
  // Note: make sure you loaded the dihdic array and phi outside
  updateDihdDic();
  double Temp( T() );
  double Rt( Rtot() );
  double C_p( Cp() );
  double ht( h() );
  double A( ht + HALF*vsqr() );
  double rho_vx( rho()*vx() );
  double rho_vy( rho()*vy() );
  double vx_vx( vx()*vx() );
  double vx_vy( vx()*vy() );
  double vy_vy( vy()*vy() );
    
  dFdW(0,0) = mult * vx();
  dFdW(0,1) = mult * rho();
  dFdW(1,0) = mult * vx_vx;
  dFdW(1,1) = mult * TWO*rho_vx; 
  dFdW(1,3) = mult * ONE;
  dFdW(2,0) = mult * vx_vy;
  dFdW(2,1) = mult * rho_vy;
  dFdW(2,2) = mult * rho_vx;
  dFdW(3,0) = mult * ( A*vx() - C_p*Temp*vx() );
  dFdW(3,1) = mult * ( rho()*vx_vx + rho()*A );
  dFdW(3,2) = mult * rho()*vx_vy;
  dFdW(3,3) = mult * vx()*C_p/Rt;
    
  //Species
  const int NUM_VAR = NUM_FLAME2D_VAR_SANS_SPECIES;
  for(int i = 0; i<ns_eqn; i++){ 
    dFdW(3,NUM_VAR+i) = mult * rho_vx*dihdic(i);    
      
    dFdW(NUM_VAR+i, 0) = mult * c(i)*vx();
    dFdW(NUM_VAR+i, 1) = mult * rho()*c(i);
    dFdW(NUM_VAR+i,NUM_VAR+i) = mult * rho_vx;        
  }

}

/********************************************************
 * Chem2D_pState::dFIdW -- Invisicd Flux Jacobian       *
 ********************************************************/
// Finite differnece check of dFxdW
void Flame2D_pState::dFIdW_FD(DenseMatrix &dFdW, const double &mult) const {

  // declares
  const int NUM_VAR = dFdW.get_n();
  static Flame2D_pState A,C;
  static Flame2D_State B,D;
  const double perturb = 5e-6;
  double a;

  // finitie differences
  for(int jcol=0; jcol<NUM_VAR; jcol++){    
    A.Copy(*this);  C.Copy(*this);
    if( jcol <NUM_FLAME2D_VAR_SANS_SPECIES) {
      A[jcol+1] += perturb*max(ONE,(*this)[jcol+1]); 
      C[jcol+1] -= perturb*max(ONE,(*this)[jcol+1]);
    } else {                                       //enforce sum(ci) = 1;
      a = perturb*max(ONE,(*this)[jcol+1]); 
      A[jcol+1] += a;
      A[NUM_VAR+1] -= a;      
      C[jcol+1] -= a;
      C[NUM_VAR+1] += a;
    }   
    A.Fx(B);   C.Fx(D);
    for(int irow=0; irow<NUM_VAR; irow++){
      dFdW(irow,jcol) = mult*( B[irow+1] - D[irow+1])/(TWO*perturb*max(ONE, (*this)[jcol+1]));      
    }
  } 

}

  
/*****************************************************************
 * Viscous fluxes  (laminar flow)                                * 
 * Viscous fluxes  (turbulent flows) are defined in single block * 
 ****************************************************************/
void Flame2D_pState::Viscous_Flux_x(const Flame2D_State &dWdx,
				    const Vector2D &qflux,
				    const Tensor2D &tau,
				    const Vector2D &Vcorr,
				    Flame2D_State &Flux,
				    const double& mult) const{
 
  //Flux.rho() += mult * ZERO;
  Flux.rhovx() += mult * tau.xx;
  Flux.rhovy() += mult * tau.xy;
  Flux.E() += mult * ( -qflux.x + vx()*tau.xx + vy()*tau.xy );
  //rho * Diffusion_Coef * grad cn 
  for( int i=0; i<ns; i++){
    Flux.rhoc(i) += mult * rho() * Diffusion_coef(i) * dWdx.c(i);
    Flux.rhoc(i) -= mult * rho() * c(i) * Vcorr.x; // <- diffusion correction term
  }

}

void Flame2D_pState::Viscous_Flux_y(const Flame2D_State &dWdy,
				    const Vector2D &qflux,
				    const Tensor2D &tau,
				    const Vector2D &Vcorr,
				    Flame2D_State &Flux,
				    const double& mult) const {

  //Flux.rho() += mult * ZERO;
  Flux.rhovx() += mult * tau.xy; 
  Flux.rhovy() += mult * tau.yy;
  Flux.E() += mult * ( -qflux.y + vx()*tau.xy + vy()*tau.yy );
  //rho * Diffusion_Coef * grad cn 
  for( int i=0; i<ns; i++){
    Flux.rhoc(i) += mult * rho() * Diffusion_coef(i) * dWdy.c(i);
    Flux.rhoc(i) -= mult * rho() * c(i) * Vcorr.y; // <- diffusion correction term
  }

}


/********************************************************
 * Chem2D_pState::dFvdWf -- Viscous Flux Jacobian       *
 ********************************************************/
void Flame2D_pState :: dFvdWf_dGvdWf( DenseMatrix &dFvdWf, 
				      DenseMatrix &dGvdWf, 
				      const Flame2D_State &dWdx, 
				      const Flame2D_State &dWdy, 
				      const int &Axisymmetric, 
				      const double &radius ) {

  // get references to gradients
  double *dhdT = cp_i;
  double *h = h_i;
  Cp_and_h( dhdT, h );
  const double *dcdx   = dWdx.c();
  const double *dcdy   = dWdy.c();
  const double &drhodx = dWdx.rho();
  const double &drhody = dWdy.rho();
  const double &dpdx   = dWdx.p();
  const double &dpdy   = dWdy.p();
  const double &dUdx   = dWdx.vx();
  const double &dUdy   = dWdy.vx();
  const double &dVdx   = dWdx.vy();
  const double &dVdy   = dWdy.vy();

   
  //---------------------------------------------------------------
  // Compute Jacobian terms
  //---------------------------------------------------------------

  // get some quantities of interest
  updateTransport();
  const double &rho_ = rho();
  const double &p_ = p();
  const double &vx_ = vx();
  const double &vy_ = vy();
  const double Kappa( kappa() );
  const double Mu( mu() );
  const double Rmix( Rtot() );
  const double rho1R(rho_*Rmix);
  const double rho2R(rho_*rho_*Rmix);
  const double rho3R(rho_*rho_*rho_*Rmix);
  const double mu_vx( Mu*vx_ );
  const double mu_vy( Mu*vy_ );

  /*********************** X - DIRECTION **************************************/

  dFvdWf(1, 5) = FOUR/THREE*Mu;
  dFvdWf(1, 8) = - TWO/THREE*Mu;
  dFvdWf(2, 6) = Mu;
  dFvdWf(2, 7) = Mu;  
   
  double Sum_q(ZERO);
  double Sum_dq(ZERO);
  for(int Num = 0; Num<ns; Num++){
    Sum_q +=  dhdT[Num]*Diffusion_coef(Num)*dcdx[Num]/Rmix;
    Sum_dq -= dhdT[Num]*Diffusion_coef(Num)*dcdx[Num]*p_/rho1R;
  } 

  dFvdWf(3,0) = ( Kappa*( -dpdx/rho2R + TWO*p_*drhodx/rho3R ) + 
		  Sum_dq );
  dFvdWf(3,1) = TWO*Mu*(TWO/THREE*dUdx - dVdy/THREE);
  dFvdWf(3,2) = Mu*(dUdy+dVdx);
  dFvdWf(3,3) = -drhodx/rho2R*Kappa + Sum_q;
  dFvdWf(3,4) = - p_/rho2R*Kappa;
  dFvdWf(3,5) = FOUR/THREE*mu_vx;
  dFvdWf(3,6) = mu_vy;
  dFvdWf(3,7) = mu_vy;
  dFvdWf(3,8) = - TWO/THREE*mu_vx;
  dFvdWf(3,9) = Kappa/rho1R; 

  // Axisymmetric
  if(Axisymmetric == AXISYMMETRIC_Y){
    dFvdWf(1,2) = - TWO/THREE*Mu/radius;
    dFvdWf(3,1) -=  TWO/THREE*mu_vy/radius;
    dFvdWf(3,2) -=  TWO/THREE*mu_vx/radius;
  } else if(Axisymmetric == AXISYMMETRIC_X){    
    dFvdWf(1,1) = - TWO/THREE*Mu/radius;
    dFvdWf(3,1) -=  FOUR/THREE*mu_vx/radius;
  }
 
  //multispecies
  for(int Num = 0; Num<ns_eqn; Num++){
    dFvdWf(3,10+Num) = rho_*Diffusion_coef(Num)*h[Num];  //+  H3???
    dFvdWf(4+Num, 10+Num) = rho_*Diffusion_coef(Num);
  }
   
  /*********************** Y - DIRECTION **************************************/

  dGvdWf(1, 6) = Mu;
  dGvdWf(1, 7) = Mu; 
  dGvdWf(2, 5) = - TWO/THREE*Mu;
  dGvdWf(2, 8) = FOUR/THREE*Mu;

  Sum_q = ZERO;  Sum_dq = ZERO;
  for(int Num = 0; Num<ns; Num++){
    Sum_q +=  dhdT[Num]*Diffusion_coef(Num)*dcdy[Num]/Rmix;
    Sum_dq -= dhdT[Num]*Diffusion_coef(Num)*dcdy[Num]*p_/rho1R;
  } 
   
  dGvdWf(3,0) = ( Kappa*(-dpdy/rho2R + TWO*p_*drhody/rho3R) + 
		  Sum_dq ); 
  dGvdWf(3,1) = Mu*(dUdy + dVdx);
  dGvdWf(3,2) = TWO*Mu*( TWO/THREE*dVdy - dUdx/THREE );
  dGvdWf(3,3) = -drhody/rho2R*Kappa + Sum_q;
  dGvdWf(3,4) = - p_/rho2R*Kappa;
  dGvdWf(3,5) = - TWO/THREE*mu_vy;
  dGvdWf(3,6) = mu_vx;
  dGvdWf(3,7) = mu_vx;
  dGvdWf(3,8) = FOUR/THREE*mu_vy;
  dGvdWf(3,9) = Kappa/rho1R; 

  //Axisymmetric 
  if(Axisymmetric == AXISYMMETRIC_Y){    
    dGvdWf(2,2) = - TWO/THREE*Mu/radius;
    dGvdWf(3,2) -=  FOUR/THREE*mu_vy/radius;
  } else if(Axisymmetric == AXISYMMETRIC_X){
    dGvdWf(2,1) = - TWO/THREE*Mu/radius;
    dGvdWf(3,1) -=  TWO/THREE*mu_vy/radius;
    dGvdWf(3,2) -=  TWO/THREE*mu_vx/radius;
  }

  //multispecies
  for(int Num = 0; Num<ns_eqn; Num++){
    dGvdWf(3,10+Num) = rho_*Diffusion_coef(Num)*h[Num];
    dGvdWf(4+Num, 10+Num) = rho_*Diffusion_coef(Num);
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

/****************************************************
 * Primitive Left Eigenvector -- (x-direction)
 ****************************************************/
void  Flame2D_pState::lp_x(const int &i, Flame2D_State& lp) const {
  double aa(a());
  // Compute the left primitive eigenvector
  if(i == 1){
    lp.rho() = ZERO;
    lp.vx() = -HALF*rho()/aa;
    lp.vy() = ZERO;
    lp.p() =  HALF/(aa*aa);
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
    lp.vx() = HALF*rho()/aa;
    lp.vy() = ZERO;
    lp.p() = HALF/(aa*aa);
    lp.c(ZERO);
  } else{ 
    lp.rho() = ZERO;
    lp.vx() = ZERO;
    lp.vy() = ZERO;
    lp.p() = ZERO;
    lp.c(ZERO);
    lp.c(i-NUM_FLAME2D_VAR_SANS_SPECIES-1) = ONE;
  } 
}

/****************************************************
 * Conserved Right Eigenvector -- (x-direction)
 ****************************************************/
void Flame2D_pState::rc_x(const int &i, 
			  Flame2D_State& rc) {
  double aa(a());
  updateDihdDic();
  // Compute the right conserved eigenvector
  if(i == 1){
    rc.rho() = ONE;
    rc.rhovx() = vx()-aa;
    rc.rhovy() = vy();
    rc.E() = H()/rho()-vx()*aa;
    rc.rhoc(c());
  } else if(i == 2) {
    rc.rho() = ONE;
    rc.rhovx() = vx();
    rc.rhovy() = vy();
    rc.E() = H()/rho()-Cp()*T();
    rc.rhoc(c());
  } else if(i == 3) {
    rc.rho() = ZERO;
    rc.rhovx() = ZERO;
    rc.rhovy() = rho();
    rc.E() = rho()*vy();
    rc.rhoc(ZERO);
  } else if(i == 4) {  
    rc.rho() = ONE;
    rc.rhovx() = vx()+aa;
    rc.rhovy() = vy();
    rc.E() = H()/rho()+vx()*aa;
    rc.rhoc(c());
  } else{ 
    rc.rho() = ZERO;
    rc.rhovx() = ZERO;
    rc.rhovy() = ZERO;
    rc.E() = rho()*dihdic(i-NUM_FLAME2D_VAR_SANS_SPECIES-1);
    rc.rhoc(ZERO);
    rc.rhoc(i-NUM_FLAME2D_VAR_SANS_SPECIES-1) = rho();
  } 
}


/****************************************************
 * Flux Dissipation  = wavespeed*((lp_x*dWrl)*rc_x) -- (x-direction)
 ****************************************************/
void Flame2D_pState::Flux_Dissipation(const int &i,
				      const Flame2D_State &dWrl, 
				      const double &wavespeed, 
				      Flame2D_State &Flux,
				      const double &mult) {

  // NOTE:  dihdic needs to be loaded first outside

  //declares
  static Flame2D_State lp;
  static Flame2D_State rc;
  double tmp, dv;

  // Compute the left primitive eigenvector
  // and the right conserved eigenvector, and store it in-place
  lp_x(i, lp);
  rc_x(i, rc);

  // Compute the wavestrengh
  dv = lp * dWrl;
  
  // compute the flux dissipation
  tmp = wavespeed*dv;
  for ( int k = 1 ; k <= n ; k++ ) Flux[k] += mult*tmp*rc[k];
  
}

/****************************************************
 * Flux Dissipation Jacobian 
 * Jac(irow, jcol) -= HALF*wavespeeds[i]*lp[jcol+1]*rc[irow+1];
 ****************************************************/
void Flame2D_pState::Flux_Dissipation_Jac(DenseMatrix &Jac,
					  const Flame2D_State &wavespeeds, 
					  const double &mult) {

  // NOTE:  dihdic needs to be loaded first outside
  const int NUM_VAR = Jac.get_n();
    
  //declares
  static Flame2D_State lp;
  static Flame2D_State rc;

  //
  // Loop through each wavespeed and each element of Jacobian(i,j)        
  //
  for (int i=1; i <= NUM_VAR; i++) { 
      
    // compute left and right eigenvectors
    lp_x(i, lp);
    rc_x(i, rc);
      
    // compute i,j element
    for(int irow =0; irow< NUM_VAR; irow++)
      for(int jcol =0; jcol< NUM_VAR; jcol++)
	Jac(irow, jcol) -= mult*wavespeeds[i]*lp[jcol+1]*rc[irow+1];
  } 
      
}

/////////////////////////////////////////////////////////////////////
/// Low Mach Number Preconditioner Based on Weiss & Smith (1995)
/////////////////////////////////////////////////////////////////////

/****************************************************
 * Primitive Left Eigenvector -- (x-direction)
 ****************************************************/
void Flame2D_pState::lp_x_precon(const int &i, 
				 const double &MR2, 
				 const double &uprimed, 
				 const double &cprimed, 
				 Flame2D_State& lp) const {
  double aa(a());
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
}

/****************************************************
 * Conserved Right Eigenvector -- (x-direction)
 ****************************************************/
void Flame2D_pState::rc_x_precon(const int &i, 
				 const double &MR2, 
				 const double &uprimed, 
				 const double &cprimed, 
				 Flame2D_State& rc) {
  double aa(a());
  updateDihdDic();
  // Compute the right conserved eigenvector
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
    rc.E() = rho()*dihdic(i-NUM_FLAME2D_VAR_SANS_SPECIES-1);
    rc.rhoc(ZERO);
    rc.rhoc(i-NUM_FLAME2D_VAR_SANS_SPECIES-1) = rho();
  } 
}



/****************************************************
 * Flux Dissipation  = wavespeed*((lp_x*dWrl)*rc_x) -- (x-direction)
 ****************************************************/
void Flame2D_pState::Flux_Dissipation_precon(const double &MR2, 
					     const Flame2D_State &dWrl, 
					     const Flame2D_State &wavespeeds, 
					     Flame2D_State &Flux_dissipation) {

  // NOTE:  dihdic needs to be loaded first outside

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
    // Compute the right conserved eigenvector, and store it in-place
    lp_x_precon(i, MR2, uprimed, cprimed, lp);
    rc_x_precon(i, MR2, uprimed, cprimed, rc);

    // Compute the wavestrengh
    dv = lp * dWrl;

    // compute the flux dissipation
    tmp = HALF*wavespeeds[i]*dv;
    for ( int k = 1 ; k <= n ; k++ ) Flux_dissipation[k] -= tmp*rc[k];

  } // endfor

}

/****************************************************
 * Flux Dissipation Jacobian 
 * Jac(irow, jcol) -= HALF*wavespeeds[i]*lp[jcol+1]*rc[irow+1];
 ****************************************************/
void Flame2D_pState::Flux_Dissipation_Jac_precon(DenseMatrix &Jac,
						 const double &MR2,
						 const Flame2D_State &wavespeeds, 
						 const double &mult) {

  // NOTE:  dihdic needs to be loaded first outside
  const int NUM_VAR = Jac.get_n();
    
  //declares
  static Flame2D_State lp;
  static Flame2D_State rc;
  double aa(a());
  double uprimed,cprimed;

  // compute preconditioned velocity and soundspeed
  u_a_precon(MR2*aa*aa,uprimed,cprimed);

  //
  // Calculate the preconditioned upwind dissipation flux.
  //
  for (int i=1; i <=  NUM_VAR; i++) {
      
    // compute left and right eigenvectors
    lp_x_precon(i, MR2, uprimed, cprimed, lp);
    rc_x_precon(i, MR2, uprimed, cprimed, rc);
      
    // compute i,j element
    for(int irow =0; irow < NUM_VAR; irow++)
      for(int jcol =0; jcol < NUM_VAR; jcol++)	   
	Jac(irow, jcol) -= HALF*wavespeeds[i]*lp[jcol+1]*rc[irow+1];
  }
    
}


/****************************************************
 * For CFL calculation
 ****************************************************/
double Flame2D_pState::u_plus_aprecon(const double &u, 
				      const int &flow_type_flag, 
				      const double &deltax) {
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
double Flame2D_pState::Mr2(const int &flow_type_flag, 
			   const double &deltax) {
  
  double aa(a());
  double MR2 = min(max((vsqr()/(aa*aa)),Mref*Mref),ONE);
  
  // Need deltax which is based on cell spacing 
  if(flow_type_flag){ 
    updateViscosity();
    MR2 = pow(max(sqrt(MR2*aa*aa), mu()/(rho()*deltax)),2.0)/(aa*aa);
  }

  return (MR2);
}


/************************************************************
 * ORIGINAL LAMINAR ONLY
 * (will be replaced by version commented out below 
 *  when corrected)
 ************************************************************/
void Flame2D_pState::Low_Mach_Number_Preconditioner(DenseMatrix &P,
						    const int &Viscous_flag, 
						    const double &deltax ) {  
  // Note: make sure you loaded the dihdic array and phi outside
  updateDihdDic();
  double phi( Phi() );
  double Temp( T() );
  double Rmix( Rtot() );
  double enthalpy( h() );
  double CP( Cp() );
  double aa( a() );
  double theta( ONE/(Mr2(Viscous_flag,deltax)*aa*aa) + ONE/(CP*Temp) );
 
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
  for(int j=0; j<ns_eqn; j++){  
     
    P(0,j+NUM_VAR) = dihdic(j)*alpham1/Omega;
    P(1,j+NUM_VAR) = vx()*dihdic(j)*alpham1/Omega;
    P(2,j+NUM_VAR) = vy()*dihdic(j)*alpham1/Omega;
    P(3,j+NUM_VAR) = dihdic(j)*(V+enthalpy)*alpham1/Omega;	
    for(int i=0; i<ns_eqn; i++){ 
      if(i==j){ 
	P(i+NUM_VAR,0) = (c(i))*(beta-V)*alpham1/Omega;
	P(i+NUM_VAR,1) = (c(i))*vx()*alpham1/Omega;
	P(i+NUM_VAR,2) = (c(i))*vy()*alpham1/Omega;
	P(i+NUM_VAR,3) = -(c(i))*alpham1/Omega;
	//diagonal
	P(i+NUM_VAR,j+NUM_VAR) = c(i)*dihdic(j)*alpham1/Omega+ONE;
      }
      else {
	P(i+NUM_VAR,j+NUM_VAR) = c(i)*dihdic(j)*alpham1/Omega;
      }
    }       
  }

  //   cout<<"\n Pin with Mref \n"<<Mref<<endl<<P;
  //P.zero(); P.identity(); //SET TO Mref=1.0;


}

/************************************************************/
void Flame2D_pState::Low_Mach_Number_Preconditioner_Inverse(DenseMatrix &Pinv,	
							    const int &Viscous_flag, 
							    const double &deltax ) {  
  // Note: make sure you loaded the dihdic array and phi outside
  updateDihdDic();
  double phi( Phi() );
  double Temp( T() );
  double Rmix( Rtot() );
  double enthalpy( h() );
  double CP( Cp() );
  double aa( a() );
  double theta( ONE/(Mr2(Viscous_flag,deltax)*aa*aa) + ONE/(CP*Temp) );

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
  for(int j=0; j<ns_eqn; j++){   
    Pinv(0,j+NUM_VAR) = -dihdic(j)*BB/AA;
    Pinv(1,j+NUM_VAR) = -vx()*dihdic(j)*BB/AA;
    Pinv(2,j+NUM_VAR) = -vy()*dihdic(j)*BB/AA;
    Pinv(3,j+NUM_VAR) = -dihdic(j)*BB*DD/AA;
    for(int i=0; i<ns_eqn; i++){  
      if(i==j){
	Pinv(i+NUM_VAR,0) = (c(i))*CC*BB/AA;
	Pinv(i+NUM_VAR,1) = -(c(i))*vx()*BB/AA;
	Pinv(i+NUM_VAR,2) = -(c(i))*vy()*BB/AA;
	Pinv(i+NUM_VAR,3) = (c(i))*BB/AA;
	//diagonal	
	Pinv(i+NUM_VAR,j+NUM_VAR) = 1.0 - c(i)*dihdic(j)*BB/AA ;
      }
      else {
	Pinv(i+NUM_VAR,j+NUM_VAR) = -c(i)*dihdic(j)*BB/AA;
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


/****************************************************************
 * Axisymmetric Source Term Jacboian (Inviscid)                 * 
 ****************************************************************/
void Flame2D_pState::dSa_idU( DenseMatrix &dSa_IdU, 
			      const Vector2D &X, 
			      const int Axisymmetric ) {

  // Note: make sure you loaded the dihdic array and phi outside
  updateDihdDic();
  double phi( Phi() );
  double enthalpy( h() );
  double CP( Cp() );
  double RTOT( Rtot() );
  double Temp( T() );
  double vsqrd( vsqr() );

  // axisymmetric case 1  -- y radial direction 
  if(Axisymmetric == AXISYMMETRIC_Y){
    cout<<"\n ISSUES IN dSa_idU AXISYMMETRIC_Y \n"; exit(1);

    // axisymmetric case 1  -- x radial direction 
  } else if(Axisymmetric == AXISYMMETRIC_X){ 
    
    dSa_IdU(0,1) -= ONE/X.x;
    dSa_IdU(1,0) += vx()*vx()/X.x;
    dSa_IdU(1,1) -= TWO*vx()/X.x;
    dSa_IdU(2,0) += vx()*vy()/X.x;
    dSa_IdU(2,1) -= vy()/X.x;
    dSa_IdU(2,2) -= vx()/X.x;
    
    double a( RTOT/(CP-RTOT) );
    double b( CP/(CP-RTOT) );
    dSa_IdU(3,0) -= a * ( vx()*(phi + vsqrd) + 
			  vx()*CP*( p()/rho() - enthalpy - HALF*vsqrd )/RTOT ) / X.x;
    dSa_IdU(3,1) -= ( enthalpy + b*HALF*vsqrd - 
		      a*HALF*(THREE*vx()*vx() + vy()*vy()) ) / X.x;
    dSa_IdU(3,2) += a*vx()*vy()/X.x;
    dSa_IdU(3,3) -= b*vx()/X.x;      
    
    //Multispecies terms
    const int NUM_VAR = NUM_FLAME2D_VAR_SANS_SPECIES;
    double d( CP/RTOT - ONE );
    for(int i=0; i<ns_eqn;i++){
      dSa_IdU(3,i+NUM_VAR) += vx()*dihdic(i)/d/X.x; 
      dSa_IdU(NUM_VAR+i,0) += vx()*c(i)/X.x;
      dSa_IdU(NUM_VAR+i,1) -= c(i)/X.x;
      dSa_IdU(NUM_VAR+i,NUM_VAR+i) -= vx()/X.x;
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
				const double& mult) {
  
  // compute laminar stresses and heat flux vector
  static Vector2D qflux;
  static Tensor2D tau;
  static Vector2D Vcorr;
  Viscous_Quantities(dWdx, dWdy, Axisymmetric, X, qflux, tau, Vcorr);


  // Determine axisymmetric source terms
  if (Axisymmetric == AXISYMMETRIC_X) {
    S.rhovx() += mult * (tau.xx - tau.zz)/X.x;
    S.rhovy() += mult * tau.xy/X.x;
    S.E() += mult * (-qflux.x + vx()*tau.xx + vy()*tau.xy)/X.x;
    for(int i=0; i<ns;i++){
      S.rhoc(i) += mult * rho()*Diffusion_coef(i)*dWdx.c(i)/X.x;
      S.rhoc(i) -= mult * rho() * c(i) * Vcorr.x / X.x; // <- diffusion correction term
    }
    
  } else if (Axisymmetric == AXISYMMETRIC_Y) {
    S.rhovx() += mult * tau.xy/X.y;
    S.rhovy() += mult * (tau.xx - tau.zz)/X.y;
    S.E() += mult * (-qflux.y + vx()*tau.xy + vy()*tau.yy)/X.y;
    for(int i=0; i<ns;i++){
      S.rhoc(i) += mult * rho()*Diffusion_coef(i)*dWdy.c(i)/X.y;
      S.rhoc(i) -= mult * rho() * c(i) * Vcorr.y / X.y; // <- diffusion correction term
    }
    
  } /* endif */
  
}


/**************************************************************** 
 * Axisymmetric Source Term Jacboian (Viscous)                  *  
 ****************************************************************/
void Flame2D_pState::dSa_vdW(DenseMatrix &dSa_VdW,
			     const Flame2D_State &dWdx,
			     const Flame2D_State &dWdy,
			     const Vector2D &X, 
			     const int Axisymmetric,
			     const double &d_dWdx_dW, 
			     const double &d_dWdy_dW) {

  updateTransport();
  const int NUM_VAR = NUM_FLAME2D_VAR_SANS_SPECIES;
  double Rmix( Rtot() );
  double Temp( T() );
  double Mu( mu() );
  double Kappa( kappa() );
  double p_( p() );
  double rho_( rho() );

  //-----------------------------------------------------------------
  // axisymmetric case 1  -- x radial direction 
  //-----------------------------------------------------------------
  if(Axisymmetric == AXISYMMETRIC_X){   

    // radius
    double radius( X.x );

    // get species properties
    //Cp(cp_i); // <- individual species specific heats
    //h(h_i); // <- individual species enthalpies
    Cp_and_h( cp_i, h_i );


    dSa_VdW(2,2) += TWO*Mu*(d_dWdx_dW - ONE/radius)/radius;    
    dSa_VdW(2,1) += Mu*d_dWdy_dW/radius;
    dSa_VdW(2,2) += Mu*d_dWdx_dW/radius ;
    
    //energy
    double Sum_q(ZERO), Sum_dq(ZERO), Sum_dhi(ZERO);
    double tmp;
    for(int i = 0; i<ns; i++){
      tmp = cp_i[i]*Diffusion_coef(i)*dWdx.c(i);
      Sum_q += tmp/Rmix;      // - ns-1 ???
      Sum_dq -= tmp*Temp;
      //Sum_dhi += cp_i[i]*Temp*rho_*Diffusion_coef(i)*dWdx.c(i)/(Rmix); // <- For full NS-1 consistency
    } 

    dSa_VdW(3,0) += ( Kappa*( -(dWdx.p() - p_/rho_*dWdx.rho())/rho_ + 
			      (p_*dWdx.rho()/(rho_*rho_) - p_*d_dWdx_dW/rho_) ) / (rho_*Rmix) + 
		      Sum_dq ) / radius;
    dSa_VdW(3,1) += ( Mu*vy()*d_dWdy_dW + 
		      TWO*Mu*(TWO/THREE*dWdx.vx() - dWdx.vy()/THREE - vx()/(THREE*radius)) +
		      TWO*vx()*Mu*(TWO/THREE*d_dWdx_dW-ONE/(THREE*radius)) ) / radius;
    dSa_VdW(3,2) += ( Mu*(dWdy.vx() + dWdx.vy()) + 
		      vy()*Mu*d_dWdx_dW - TWO/THREE*Mu*vx()*d_dWdy_dW ) / radius;
    dSa_VdW(3,3) += ( Kappa*(d_dWdx_dW - dWdx.rho()/rho_)/(rho_*Rmix) + Sum_q )/radius;
   
    //Multispecies
    double rhoD;
    for(int i = 0; i<ns_eqn; i++){ 
      rhoD = rho_*Diffusion_coef(i);
      dSa_VdW(3,NUM_VAR+i) += ( rhoD * h_i[i] * d_dWdx_dW ) / radius;  
      //- specdata[Num].Rs()*Sum_dhi)/radius; // <- For full NS-1 consistency
      dSa_VdW(NUM_VAR+i,NUM_VAR+i) += rhoD*d_dWdx_dW/radius;
    }
    

    //-----------------------------------------------------------------
    // axisymmetric case 1  -- y radial direction 
    //-----------------------------------------------------------------
  } else if(Axisymmetric == AXISYMMETRIC_Y){ 
    cout<<"\n AXISYMMETRIC_Y NOT DONE YET ";
  }
       
}


/****************************************************
 * Source terms associated with finite-rate chemistry
 ****************************************************/  
void Flame2D_pState::Sw(Flame2D_State &S, 
			const double& mult) const {
  Mix.getRates( p(), c(), r );
  for(int i=0; i<ns_eqn; i++) S.rhoc(i) += mult*r[i];
}

/****************************************************
 * Chemistry source term jacobian
 ****************************************************/  
void Flame2D_pState::dSwdU( DenseMatrix &dSdU ) const {
  Mix.dSwdU( dSdU, rho(), p(), c(), NUM_FLAME2D_VAR_SANS_SPECIES, NSm1 );
}


/****************************************************
 * Max diagonal of chemistry source jacobian
 ****************************************************/  
double Flame2D_pState::dSwdU_max_diagonal(void) const {
  return Mix.dSwdU_max_diagonal( rho(), p(), c() );
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

/****************************************************
 * Source terms associated with gravitational forces
 ****************************************************/
void Flame2D_pState::dSgdU(DenseMatrix &dSgdU) const {
  dSgdU(2,0) += gravity_z;
  dSgdU(3,2) += gravity_z;
}

/////////////////////////////////////////////////////////////////////
/// Inviscid Flux Jacobians
/////////////////////////////////////////////////////////////////////

/************************************************************************
 * Chem2D_pState::dWdU -- Primitive/Conserved transformation Jacobian   *
 ************************************************************************/
void Flame2D_pState::dWdU(DenseMatrix &dWdQ) {

  // Note: make sure you loaded the dihdic array and phi outside
  updateDihdDic();
  double phi( Phi() );
  double Temp( T() );
  double Rt( Rtot() );
  double C_p( Cp() );
  double denominator( C_p/Rt - ONE );

  dWdQ(0,0) = ONE;
  dWdQ(1,0) = -vx()/rho();
  dWdQ(1,1) = ONE/rho();
  dWdQ(2,0) = -vy()/rho();
  dWdQ(2,2) = ONE/rho(); 

  dWdQ(3,0) = (HALF*vsqr() - h() + C_p*Temp + phi)/denominator;
  dWdQ(3,1) = -vx()/denominator;
  dWdQ(3,2) = -vy()/denominator;
  dWdQ(3,3) = ONE/denominator;

  int NUM_VAR = NUM_FLAME2D_VAR_SANS_SPECIES;

  //Species
  for(int i=0; i<ns_eqn; i++){   
    dWdQ(3, NUM_VAR+i) = - dihdic(i)/denominator;
    dWdQ(NUM_VAR+i, 0) = - c(i)/rho();
    dWdQ(NUM_VAR+i, NUM_VAR+i) = ONE/rho();
  }

}

/************************************************************************
 * Chem2D_pState::dWdU -- Primitive/Conserved transformation Jacobian   *
 ************************************************************************/
// Finite differnece check of dWdU
// shows error in (3,0) ie dp/rho due to pertubing rho and cState T() calc.
void Flame2D_pState::dWdU_FD(DenseMatrix &dWdQ) const {

  // declares
  const int NUM_VAR = dWdQ.get_n();
  static Flame2D_State UU;
  static Flame2D_State A,C;
  static Flame2D_pState B,D;
  const double perturb = 5e-6;
  double a;

  // set UU
  getU(UU);

  for(int jcol=0; jcol<NUM_VAR; jcol++){    
    A = UU;    C = UU;
    if( jcol <NUM_FLAME2D_VAR_SANS_SPECIES) {
      A[jcol+1] += perturb*max(ONE,A[jcol+1]); 
      C[jcol+1] -= perturb*max(ONE,C[jcol+1]);
    } else {                                       //enforce sum(ci) = 1;
      a =  perturb*max(ONE,A[jcol+1]); 
      A[jcol+1] += a;
      A[NUM_VAR+1] -= a;      
      C[jcol+1] -= a;
      C[NUM_VAR+1] += a;
    }
    B.setU(A);   D.setU(C);    
    for(int irow=0; irow<NUM_VAR; irow++){
      dWdQ(irow,jcol) = ( B[irow+1] - D[irow+1])/(TWO*perturb*max(ONE,UU[jcol+1]));     
    }
  } 

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
    cout << "\nAdjusting Inflow\n";
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
    cout << "\nAdjusting Outflow\n";
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
void Flame2D_pState :: HLLE_wavespeeds(Flame2D_pState &Wl,
				       Flame2D_pState &Wr,
				       const Vector2D &norm_dir,
				       Vector2D &wavespeed) {

  // declares
  static Flame2D_pState Wa;
  static Flame2D_State lambdas_l, lambdas_r, lambdas_a;
  int NUM_VAR_CHEM2D = ( NumVar() );

  // temporary storage
  double ur( ((const Flame2D_pState&)Wr).vx() ), 
    vr( ((const Flame2D_pState&)Wr).vy() ), 
    ul( ((const Flame2D_pState&)Wl).vx() ), 
    vl( ((const Flame2D_pState&)Wl).vy() );

  /* Use rotated values to calculate eignvalues */
  Wl.Rotate(norm_dir);
  Wr.Rotate(norm_dir);

  /* Evaluate the Roe-average primitive solution state. */
  Wa.RoeAverage(Wl, Wr);
    
  /* Evaluate the left, right, and average state eigenvalues. */
  Wl.lambda_x(lambdas_l);
  Wr.lambda_x(lambdas_r);
  Wa.lambda_x(lambdas_a);

  /* Determine the intermediate state flux. */
  wavespeed.x = min(lambdas_l[1],
		    lambdas_a[1]);   //u-a
  wavespeed.y = max(lambdas_r[4],
		    lambdas_a[4]);   //u+a
  wavespeed.x = min(wavespeed.x, ZERO); //lambda minus
  wavespeed.y = max(wavespeed.y, ZERO); //lambda plus 

  /* Rotate Back -> avoid roundoff by setting the exact values */
  Wl.setVelocity( ul, vl );
  Wr.setVelocity( ur, vr );


  //    cout<< "lambda - = "<<wavespeed.x<< "lambda + = "<<wavespeed.y<<endl;

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
    dUrl.DeltaU( Wr, Wl ); // Evaluate the jumps in the conserved solution states.
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
void Flame2D_State :: FluxHLLE_n(Flame2D_pState &Wl,
				 Flame2D_pState &Wr,
				 const Vector2D &norm_dir) {

  /* Apply the frame rotation and evaluate left and right
     solution states in the local rotated frame defined
     by the unit normal vector. */

  double ur( ((const Flame2D_pState&)Wr).vx() ), 
    vr( ((const Flame2D_pState&)Wr).vy() ), 
    ul( ((const Flame2D_pState&)Wl).vx() ), 
    vl( ((const Flame2D_pState&)Wl).vy() );

  /* Use rotated values to calculate eignvalues */
  Wl.Rotate(norm_dir);
  Wr.Rotate(norm_dir);

  /* Evaluate the intermediate state solution 
     flux in the rotated frame. */
  
  FluxHLLE_x(Wl, Wr);
  
  /* Rotate back to the original Cartesian reference
     frame and return the solution flux. */

  RotateBack(norm_dir);

  /* Rotate Back -> avoid roundoff by setting the exact values */

  Wl.setVelocity( ul, vl );
  Wr.setVelocity( ur, vr );

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

    dUrl.DeltaU( Wr, Wl );
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
void Flame2D_State :: FluxLinde_n(Flame2D_pState &Wl,
				  Flame2D_pState &Wr,
				  const Vector2D &norm_dir) {

  /* Apply the frame rotation and evaluate left and right
     solution states in the local rotated frame defined
     by the unit normal vector. */

  double ur( ((const Flame2D_pState&)Wr).vx() ), 
    vr( ((const Flame2D_pState&)Wr).vy() ), 
    ul( ((const Flame2D_pState&)Wl).vx() ), 
    vl( ((const Flame2D_pState&)Wl).vy() );

  /* Use rotated values to calculate eignvalues */
  Wl.Rotate(norm_dir);
  Wr.Rotate(norm_dir);

  /* Evaluate the intermediate state solution 
     flux in the rotated frame. */
  
  FluxLinde(Wl, Wr);
  
  /* Rotate back to the original Cartesian reference
     frame and return the solution flux. */

  RotateBack(norm_dir);

  /* Rotate Back -> avoid roundoff by setting the exact values */

  Wl.setVelocity( ul, vl );
  Wr.setVelocity( ur, vr );

 
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
void Flame2D_State :: FluxRoe_x(Flame2D_pState &Wl,
				Flame2D_pState &Wr,
				const int &Preconditioning,
				const int &flow_type_flag,
				const double &deltax) {
  
  static Flame2D_pState Wa;
  static Flame2D_State dWrl, wavespeeds, lambdas_l, lambdas_r, lambdas_a;
  
  // Evaluate the Roe-average primitive solution state.
  Wa.RoeAverage(Wl, Wr);
  const Flame2D_pState& Waa = Wa;

  // Evaluate the jumps in the primitive solution states.
  dWrl.Delta(Wr,Wl);
  
  //---------------------------------------------------------------
  // No Preconditioning
  //---------------------------------------------------------------
  if(!Preconditioning){
    // Evaluate the left, right, and average state eigenvalues.

    Wl.lambda_x(lambdas_l);
    Wr.lambda_x(lambdas_r);
    Wa.lambda_x(lambdas_a);

    // Determine the intermediate state flux.
    if (Waa.vx() >= ZERO) {
      Wl.Fx(*this);   
      wavespeeds.HartenFix_Neg(lambdas_a,
			       lambdas_l,
			       lambdas_r);
      
      for (int i=0 ; i < NumEqn(); i++) {
	if (wavespeeds[i+1] < ZERO) {
	  //*this += wavespeeds[i+1]*((Waa.lp_x(i+1)*dWrl)*Waa.rc_x(i+1));
	  Wa.Flux_Dissipation(i+1, dWrl, wavespeeds[i+1], *this, 1.0);
	}
      } 
    } else {
      Wr.Fx(*this);
      wavespeeds.HartenFix_Pos(lambdas_a,
			       lambdas_l,
			       lambdas_r);
      for (int i=0; i < NumEqn(); i++) {
	if (wavespeeds[i+1] > ZERO) {
	  //*this -= wavespeeds[i+1]*((Waa.lp_x(i+1)*dWrl)*Waa.rc_x(i+1));
	  Wa.Flux_Dissipation(i+1, dWrl, wavespeeds[i+1], *this, -1.0);
	}
      } 
    } 
    
    //---------------------------------------------------------------
    // LOW MACH NUMBER PRECONDITIONING
    //---------------------------------------------------------------
  } else if(Preconditioning){
    // Evaluate the left, right, and average state eigenvalues.

    //calculating Mr^2 and passing to save computation,
    //not conceptually nice but saves from recalculating
    double MR2a( Wa.Mr2(flow_type_flag,deltax) );
    Wl.lambda_preconditioned_x(lambdas_l, Wl.Mr2(flow_type_flag,deltax)); 
    Wr.lambda_preconditioned_x(lambdas_r, Wr.Mr2(flow_type_flag,deltax));
    Wa.lambda_preconditioned_x(lambdas_a, MR2a);
    
    // Evaluate the jumps in the primitive solution states.
    wavespeeds.HartenFix_Abs(lambdas_a,
			     lambdas_l,
			     lambdas_r);
     
    // Evaluate the low-Mach-number local preconditioner for the Roe-averaged state.
    static DenseMatrix P(NumEqn(),NumEqn(),ZERO);
    //P.zero(); // <- no need, always writing to the same spot
    Wa.Low_Mach_Number_Preconditioner(P, flow_type_flag, deltax);
    
    // Determine the intermediate state flux.
    Wl.Fx(*this, HALF);
    Wr.addFx(*this, HALF);

    // compute dissipation flux
    static Flame2D_State Flux_dissipation;
    Wa.Flux_Dissipation_precon( MR2a, dWrl, wavespeeds, Flux_dissipation);
    
    // Add preconditioned upwind dissipation flux.
    for ( int i = 0 ; i < NumEqn() ; i++ ) {
      for ( int j = 0 ; j < NumEqn() ; j++ ) {
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
void Flame2D_State :: FluxRoe_n(Flame2D_pState &Wl,
				Flame2D_pState &Wr,
				const Vector2D &norm_dir,
				const int &Preconditioning,
				const int &flow_type_flag,
				const double &delta_n ) {

  /* Apply the frame rotation and evaluate left and right
     solution states in the local rotated frame defined
     by the unit normal vector. */

  double ur( ((const Flame2D_pState&)Wr).vx() ), 
    vr( ((const Flame2D_pState&)Wr).vy() ), 
    ul( ((const Flame2D_pState&)Wl).vx() ), 
    vl( ((const Flame2D_pState&)Wl).vy() );

  /* Use rotated values to calculate eignvalues */
  Wl.Rotate(norm_dir);
  Wr.Rotate(norm_dir);

  /* Evaluate the intermediate state solution 
     flux in the rotated frame. */
  
  FluxRoe_x(Wl, Wr, Preconditioning,flow_type_flag,delta_n);
 
  /* Rotate back to the original Cartesian reference
     frame and return the solution flux. */

  RotateBack(norm_dir);

  /* Rotate Back -> avoid roundoff by setting the exact values */

  Wl.setVelocity( ul, vl );
  Wr.setVelocity( ur, vr );

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
void Flame2D_State::FluxAUSMplus_up_n(Flame2D_pState &Wl,
				      Flame2D_pState &Wr,
				      const Vector2D &norm_dir) {


  /* Apply the frame rotation and evaluate left and right
     solution states in the local rotated frame defined
     by the unit normal vector. */

  double ur( ((const Flame2D_pState&)Wr).vx() ), 
    vr( ((const Flame2D_pState&)Wr).vy() ), 
    ul( ((const Flame2D_pState&)Wl).vx() ), 
    vl( ((const Flame2D_pState&)Wl).vy() );

  /* Use rotated values to calculate eignvalues */
  Wl.Rotate(norm_dir);
  Wr.Rotate(norm_dir);

  /* Evaluate the intermediate state solution 
     flux in the rotated frame. */
  
  FluxAUSMplus_up(Wl, Wr);

  /* Rotate back to the original Cartesian reference
     frame and return the solution flux. */

  RotateBack(norm_dir);

  /* Rotate Back -> avoid roundoff by setting the exact values */

  Wl.setVelocity( ul, vl );
  Wr.setVelocity( ur, vr );

}


/////////////////////////////////////////////////////////////////////
/// Viscous Reconstruction Functions
/////////////////////////////////////////////////////////////////////


/**********************************************************************
 * Routine: Viscous_Flux_n                                            *
 *                                                                    *
 * This function returns the intermediate state solution viscous flux *
 * given the primitive variable solution state and the gradients of   *
 * the primitive variables.                                           *
 *                                                                    *
 **********************************************************************/
void Flame2D_State::Viscous_Flux_n(Flame2D_pState &W,
				   const Flame2D_State &dWdx,
				   const Flame2D_State &dWdy,
				   const int Axisymmetric,
				   const Vector2D X,
				   const Vector2D &norm_dir,
				   const double &mult) {

  static Vector2D qflux;
  static Tensor2D tau;
  static Vector2D Vcorr;

  // the cosines
  static const Vector2D i(1,0), j(0,1);
  double nx( dot(i,norm_dir) ), ny( dot(j,norm_dir) );

  // compute the stress tensor and heat flux vector
  W.Viscous_Quantities( dWdx, dWdy, Axisymmetric, X, qflux, tau, Vcorr );
  
  // Fn = (Fx*dot(i,norm_dir) + Fy*dot(j,norm_dir))

  // Fx Constribution
  if (fabs(nx)>TOLER) {
    W.Viscous_Flux_x(dWdx, qflux, tau, Vcorr, *this, mult*nx);
  }
  // Fy contribution
  if (fabs(ny)>TOLER) {
    W.Viscous_Flux_y(dWdy, qflux, tau, Vcorr, *this, mult*ny);
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
void Flame2D_State::Viscous_FluxHybrid_n(Flame2D_pState &W,
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
  static Vector2D Vcorr;

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
  W.Viscous_Quantities( dWdx, dWdy, Axisymmetric, X, qflux, tau, Vcorr );
  
  // Fn = (Fx*dot(i,norm_dir) + Fy*dot(j,norm_dir))
  // Fx Constribution
  if (fabs(nx)>TOLER) {
    W.Viscous_Flux_x(dWdx, qflux, tau, Vcorr, *this, mult*nx);
  }
  // Fy contribution
  if (fabs(ny)>TOLER) {
    W.Viscous_Flux_y(dWdy, qflux, tau, Vcorr, *this, mult*ny);
  }

}

/////////////////////////////////////////////////////////////////////
/// Jump Conditions
/////////////////////////////////////////////////////////////////////

/**********************************************************************
 * Routine: FlameJumpLowMach_x                                        *
 *                                                                    *
 * This routine solves a 1D premixed flame in the x-direction using   *
 * jump condition equations, subject to a low Mach number assumption, *
 * that is,                                                           *
 *                                                                    *
 *   rho1*T1=rho2*T2  (see Poinsot/Veynante, "Numerical Combustion")  *
 *                                                                    *
 * Returns the burnt mixture state vector.                            *
 *                                                                    *
 **********************************************************************/
void Flame2D_pState::FlameJumpLowMach(const Flame2D_pState &Wu) {

  // Product species' mass fractions should already be set

  // declares
  double rho2, u2, T2, p2, Cp2;      // Unknown variables.
  double f[3];                       // LHS vector.
  double delta_vars[3];              // Change in variables.
  double detJ, detf0, detf1, detf2;  // Determinants.

  // iteration parameters
  double norm(MILLION), norm_tolerance(1.0e-6);
  int iter_count, iter_max(20);

  // Initial guess.
  rho2 = DENSITY_STDATM;
  u2   = TWO*Wu.vx();
  T2   = 1000;

  //
  // Apply Newton's method to solve for rho2, u2, T2.
  //
  for( iter_count=0; iter_count<iter_max; iter_count++ ) {

    f[0] = -(rho2*u2 - Wu.rho()*Wu.vx());
    f[1] = -(rho2*T2 - Wu.rho()*Wu.T());
    f[2] = -( (h(T2)+HALF*sqr(u2)) - (Wu.h()+HALF*sqr(Wu.vx())) );

    Cp2 = Cp(T2);

    detJ  = u2*(-u2*rho2) - T2*(rho2*Cp2);
    detf0 = -rho2*(f[0]*u2-f[2]*rho2) + Cp2*(-f[1]*rho2);
    detf1 = u2*(f[1]*Cp2-f[2]*rho2) - T2*(f[0]*Cp2);
    detf2 = u2*(-u2*f[1]) - T2*(f[2]*rho2-f[0]*u2);

    delta_vars[0] = detf0/detJ;
    delta_vars[1] = detf1/detJ;
    delta_vars[2] = detf2/detJ;

    rho2 = rho2 + delta_vars[0];
    u2   = u2 + delta_vars[1];
    T2   = T2 + delta_vars[2];

    norm = sqrt( sqr(delta_vars[0]) + sqr(delta_vars[1]) + sqr(delta_vars[2]) );

    // Newton method converged.
    if (norm <= norm_tolerance) break;

  }

  // Check if converged
  if (norm > norm_tolerance || norm!=norm) {
    cerr << "\nFlame2D_pState::FlameJumpLowMach_x() - Newton's method did not converge.\n";
    exit(-1);
  }


  // assign
  rho() = rho2;
  vx() = u2;
  vy() = Wu.vy();   // Jump in tangential velocity is zero.
  p() = Wu.p() + Wu.rho()*sqr(Wu.vx())*(ONE - vx()/Wu.vx());
  setGas();

}


/////////////////////////////////////////////////////////////////////
/// Exact Test Case Solution Functions
/////////////////////////////////////////////////////////////////////


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
  updateViscosity();
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
void Flame2D_pState::FlatPlate(const Flame2D_pState &Winf,
			       const Vector2D X,
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

  // set bulk flow conditions
  Copy(Winf);

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
  updateViscosity();
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
				       const Vector2D &norm_dir) {

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
  updateViscosity();
  return mu()*dWdn;

}
