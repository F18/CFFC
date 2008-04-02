/**
 * \file  LES3DPolytropicState.cc
 * \brief Definition of member functions for 3D LES solution classes
 *        associated with solution of flows of a polytropic gas using 
 *        a large eddy simulation (LES) technique.
 */

/* Include LES3DPolytropicState header file. */

#ifndef _LES3D_POLYTROPIC_STATE_INCLUDED
#include "LES3DPolytropicState.h"
#endif // LES3D_POLYTROPIC_STATE_INCLUDED   


/*----------------------------------------------------------------------------*
 *                LESS3D_Polytropic_pState member functions                   *
 *----------------------------------------------------------------------------*/

SFS_model_Parameters LES3D_Polytropic_pState::SFS_model = {SFS_MODEL_SMAGORINSKY,
                                                           10000000};
Filter_Parameters LES3D_Polytropic_pState::filter = {2,
                                                     FILTER_TYPE_IMPLICIT};

void LES3D_Polytropic_pState::Set_LES_parameters(SFS_model_Parameters &SFS_model_,Filter_Parameters &filter_){
    SFS_model = SFS_model_;
    filter = filter_;
}



/*
 * Turbulent transport coefficients
 * --------------------------------
 */

/**
 * LES3D_Polytropic_pState::mu_t
 * Return dynamic eddy (turbulent) viscosity using Smagorinsky model
 */
double LES3D_Polytropic_pState::mu_t(const LES3D_Polytropic_pState &dWdx,
                                     const LES3D_Polytropic_pState &dWdy,
                                     const LES3D_Polytropic_pState &dWdz,
                                     const double &Volume) {
    return(rho*nu_t(dWdx,dWdy,dWdz,Volume));
}

/**
 * LES3D_Polytropic_pState::nu_t
 * Return kinematic eddy (turbulent) viscosity using Smagorinsky model
 */
double LES3D_Polytropic_pState::nu_t(const LES3D_Polytropic_pState &dWdx,
                                     const LES3D_Polytropic_pState &dWdy,
                                     const LES3D_Polytropic_pState &dWdz,
                                     const double &Volume) {
    double Delta = filter_width(Volume);
    double Cs = SFS_model.smagorinsky_coefficient;
    return(sqr(Cs*Delta)*abs_strain_rate(dWdx,dWdy,dWdz));
}

/**
 * LES3D_Polytropic_pState::kappa_t
 * Return turbulent thermal conductivity.     
 */
double LES3D_Polytropic_pState::kappa_t(const LES3D_Polytropic_pState &dWdx,
                                        const LES3D_Polytropic_pState &dWdy,
                                        const LES3D_Polytropic_pState &dWdz, 
                                        const double &Volume) {
  return (mu_t(dWdx,dWdy,dWdz,Volume)*Cp/Pr_t()); 
}

/**
 * LES3D_Polytropic_pState::Pr_t
 * Return turbulent Prandtl number.
 */
double LES3D_Polytropic_pState::Pr_t(void) {
   return (0.9); 
}

/**
 * LES3D_Polytropic_pState::filter_width
 * Return LES characteristic filter width
 */
double LES3D_Polytropic_pState::filter_width(const double &Volume) const {
  return (filter.FGR*pow(Volume,1.0/3.0)); 
}








/*
 * Subfilter scale turbulent stress tensor and heat flux vector
 * ------------------------------------------------------------
 */

/**
 * LES3D_Polytropic_pState::tau_t
 * Return subfilter scale (turbulent) stress tensor using Smagorinsky model
 */
Tensor3D LES3D_Polytropic_pState::tau_t(const LES3D_Polytropic_pState &dWdx, 
                                        const LES3D_Polytropic_pState &dWdy,
                                        const LES3D_Polytropic_pState &dWdz,
                                        const double &Volume) {
    Tensor3D SFS_stress;
    double k = SFS_Kinetic_Energy(dWdx,dWdy,dWdz,Volume);
    double mu_t_ = mu_t(dWdx,dWdy,dWdz,Volume);
    SFS_stress.xx = mu_t_*(FOUR*dWdx.v.x - TWO*dWdy.v.y - TWO*dWdz.v.z)/THREE - (TWO/THREE)*rho*k;
    SFS_stress.xy = mu_t_*(dWdx.v.y + dWdy.v.x);
    SFS_stress.xz = mu_t_*(dWdx.v.z + dWdz.v.x);
    SFS_stress.yy = mu_t_*(FOUR*dWdy.v.y - TWO*dWdx.v.x - TWO*dWdz.v.z)/THREE - (TWO/THREE)*rho*k;
    SFS_stress.yz = mu_t_*(dWdy.v.z + dWdz.v.y);
    SFS_stress.zz = mu_t_*(FOUR*dWdz.v.z - TWO*dWdx.v.x - TWO*dWdy.v.y)/THREE - (TWO/THREE)*rho*k;
    return (SFS_stress);
}

Tensor3D LES3D_Polytropic_pState::tau_t(const double &mu_t,
                                        const LES3D_Polytropic_pState &dWdx, 
                                        const LES3D_Polytropic_pState &dWdy,
                                        const LES3D_Polytropic_pState &dWdz,
                                        const double &Volume) {
    Tensor3D SFS_stress;
    double k = SFS_Kinetic_Energy(dWdx,dWdy,dWdz,Volume);
    SFS_stress.xx = mu_t*(FOUR*dWdx.v.x - TWO*dWdy.v.y - TWO*dWdz.v.z)/THREE - (TWO/THREE)*rho*k;
    SFS_stress.xy = mu_t*(dWdx.v.y + dWdy.v.x);
    SFS_stress.xz = mu_t*(dWdx.v.z + dWdz.v.x);
    SFS_stress.yy = mu_t*(FOUR*dWdy.v.y - TWO*dWdx.v.x - TWO*dWdz.v.z)/THREE - (TWO/THREE)*rho*k;
    SFS_stress.yz = mu_t*(dWdy.v.z + dWdz.v.y);
    SFS_stress.zz = mu_t*(FOUR*dWdz.v.z - TWO*dWdx.v.x - TWO*dWdy.v.y)/THREE - (TWO/THREE)*rho*k;
    return (SFS_stress);
}

/**
 * LES3D_Polytropic_pState::tau_t_x
 * Return components LES3D_Polytropic_pState::tau_t in the x-direction.                           
 */
Tensor3D LES3D_Polytropic_pState::tau_t_x(const LES3D_Polytropic_pState &dWdx, 
                                          const LES3D_Polytropic_pState &dWdy,
                                          const LES3D_Polytropic_pState &dWdz,
                                          const double &Volume) {
    Tensor3D SFS_stress;
    double k = SFS_Kinetic_Energy(dWdx,dWdy,dWdz,Volume);
    double mu_t_ = mu_t(dWdx,dWdy,dWdz,Volume);
    SFS_stress.xx = mu_t_*(FOUR*dWdx.v.x - TWO*dWdy.v.y - TWO*dWdz.v.z)/THREE - (TWO/THREE)*rho*k;
    SFS_stress.xy = mu_t_*(dWdx.v.y + dWdy.v.x);
    SFS_stress.xz = mu_t_*(dWdx.v.z + dWdz.v.x);
    return (SFS_stress);
}

Tensor3D LES3D_Polytropic_pState::tau_t_x(const double &mu_t,
                                          const LES3D_Polytropic_pState &dWdx, 
                                          const LES3D_Polytropic_pState &dWdy,
                                          const LES3D_Polytropic_pState &dWdz,
                                          const double &Volume) {
    Tensor3D SFS_stress;
    double k = SFS_Kinetic_Energy(dWdx,dWdy,dWdz,Volume);
    SFS_stress.xx = mu_t*(FOUR*dWdx.v.x - TWO*dWdy.v.y - TWO*dWdz.v.z)/THREE - (TWO/THREE)*rho*k;
    SFS_stress.xy = mu_t*(dWdx.v.y + dWdy.v.x);
    SFS_stress.xz = mu_t*(dWdx.v.z + dWdz.v.x);
    return (SFS_stress);
}

/**
 * LES3D_Polytropic_pState::tau_t_y
 * Return components LES3D_Polytropic_pState::tau_t in the y-direction.                           
 */
Tensor3D LES3D_Polytropic_pState::tau_t_y(const LES3D_Polytropic_pState &dWdx, 
                                          const LES3D_Polytropic_pState &dWdy,
                                          const LES3D_Polytropic_pState &dWdz,
                                          const double &Volume) {
    Tensor3D SFS_stress;
    double k = SFS_Kinetic_Energy(dWdx,dWdy,dWdz,Volume);
    double mu_t_ = mu_t(dWdx,dWdy,dWdz,Volume);
    SFS_stress.xy = mu_t_*(dWdx.v.y + dWdy.v.x);
    SFS_stress.yy = mu_t_*(FOUR*dWdy.v.y - TWO*dWdx.v.x - TWO*dWdz.v.z)/THREE - (TWO/THREE)*rho*k;
    SFS_stress.yz = mu_t_*(dWdy.v.z + dWdz.v.y);
    return (SFS_stress);
}

Tensor3D LES3D_Polytropic_pState::tau_t_y(const double &mu_t,
                                          const LES3D_Polytropic_pState &dWdx, 
                                          const LES3D_Polytropic_pState &dWdy,
                                          const LES3D_Polytropic_pState &dWdz,
                                          const double &Volume) {
    Tensor3D SFS_stress;
    double k = SFS_Kinetic_Energy(dWdx,dWdy,dWdz,Volume);
    SFS_stress.xy = mu_t*(dWdx.v.y + dWdy.v.x);
    SFS_stress.yy = mu_t*(FOUR*dWdy.v.y - TWO*dWdx.v.x - TWO*dWdz.v.z)/THREE - (TWO/THREE)*rho*k;
    SFS_stress.yz = mu_t*(dWdy.v.z + dWdz.v.y);
    return (SFS_stress);
}

/**
 * LES3D_Polytropic_pState::tau_t_z
 * Return components LES3D_Polytropic_pState::tau_t in the z-direction.                           
 */
Tensor3D LES3D_Polytropic_pState::tau_t_z(const LES3D_Polytropic_pState &dWdx, 
                                          const LES3D_Polytropic_pState &dWdy,
                                          const LES3D_Polytropic_pState &dWdz,
                                          const double &Volume) {
    Tensor3D SFS_stress;
    double k = SFS_Kinetic_Energy(dWdx,dWdy,dWdz,Volume);
    double mu_t_ = mu_t(dWdx,dWdy,dWdz,Volume);
    SFS_stress.xz = mu_t_*(dWdx.v.z + dWdz.v.x);
    SFS_stress.yz = mu_t_*(dWdy.v.z + dWdz.v.y);
    SFS_stress.zz = mu_t_*(FOUR*dWdz.v.z - TWO*dWdx.v.x - TWO*dWdy.v.y)/THREE - (TWO/THREE)*rho*k;
    return (SFS_stress);
}

Tensor3D LES3D_Polytropic_pState::tau_t_z(const double &mu_t,
                                          const LES3D_Polytropic_pState &dWdx, 
                                          const LES3D_Polytropic_pState &dWdy,
                                          const LES3D_Polytropic_pState &dWdz,
                                          const double &Volume) {
    Tensor3D SFS_stress;
    double k = SFS_Kinetic_Energy(dWdx,dWdy,dWdz,Volume);
    SFS_stress.xz = mu_t*(dWdx.v.z + dWdz.v.x);
    SFS_stress.yz = mu_t*(dWdy.v.z + dWdz.v.y);
    SFS_stress.zz = mu_t*(FOUR*dWdz.v.z - TWO*dWdx.v.x - TWO*dWdy.v.y)/THREE - (TWO/THREE)*rho*k;
    return (SFS_stress);
}

/**
 * LES3D_Polytropic_pState::q_t
 * Return turbulent heat flux vector.       
 */
Vector3D LES3D_Polytropic_pState::q_t(const LES3D_Polytropic_pState &dWdx, 
                                      const LES3D_Polytropic_pState &dWdy,
                                      const LES3D_Polytropic_pState &dWdz,
                                      const double &Volume) {
    double kappa_t_ = kappa_t(dWdx,dWdy,dWdz,Volume);
    Vector3D heat_flux;
    heat_flux.x = -kappa_t_*(1.0/(rho*R)) * (dWdx.p -(p/rho)*dWdx.rho);
    heat_flux.y = -kappa_t_*(1.0/(rho*R)) * (dWdy.p -(p/rho)*dWdy.rho);
    heat_flux.z = -kappa_t_*(1.0/(rho*R)) * (dWdz.p -(p/rho)*dWdz.rho);
    return (heat_flux);
}

Vector3D LES3D_Polytropic_pState::q_t(const double &kappa_t,
                                      const LES3D_Polytropic_pState &dWdx, 
                                      const LES3D_Polytropic_pState &dWdy,
                                      const LES3D_Polytropic_pState &dWdz,
                                      const double &Volume) {
    Vector3D heat_flux;
    heat_flux.x = -kappa_t*(1.0/(rho*R)) * (dWdx.p -(p/rho)*dWdx.rho);
    heat_flux.y = -kappa_t*(1.0/(rho*R)) * (dWdy.p -(p/rho)*dWdy.rho);
    heat_flux.z = -kappa_t*(1.0/(rho*R)) * (dWdz.p -(p/rho)*dWdz.rho);
    return (heat_flux);
}

/**
 * LES3D_Polytropic_pState::q_t_x
 * Return component of LES3D_Polytropic_pState::q_t in the x-direction.   
 */
Vector3D LES3D_Polytropic_pState::q_t_x(const LES3D_Polytropic_pState &dWdx, 
                                        const LES3D_Polytropic_pState &dWdy,
                                        const LES3D_Polytropic_pState &dWdz,
                                        const double &Volume) {
    double kappa_t_ = kappa_t(dWdx,dWdy,dWdz,Volume);
    Vector3D heat_flux;
    heat_flux.x = -kappa_t_*(1.0/(rho*R)) * (dWdx.p -(p/rho)*dWdx.rho);
    return (heat_flux);
}

Vector3D LES3D_Polytropic_pState::q_t_x(const double &kappa_t,
                                        const LES3D_Polytropic_pState &dWdx, 
                                        const LES3D_Polytropic_pState &dWdy,
                                        const LES3D_Polytropic_pState &dWdz,
                                        const double &Volume) {
    Vector3D heat_flux;
    heat_flux.x = -kappa_t*(1.0/(rho*R)) * (dWdx.p -(p/rho)*dWdx.rho);
    return (heat_flux);
}

/**
 * LES3D_Polytropic_pState::q_t_y
 * Return component of LES3D_Polytropic_pState::q_t in the y-direction.   
 */
Vector3D LES3D_Polytropic_pState::q_t_y(const LES3D_Polytropic_pState &dWdx, 
                                        const LES3D_Polytropic_pState &dWdy,
                                        const LES3D_Polytropic_pState &dWdz,
                                        const double &Volume) {
    double kappa_t_ = kappa_t(dWdx,dWdy,dWdz,Volume);
    Vector3D heat_flux;
    heat_flux.y = -kappa_t_*(1.0/(rho*R)) * (dWdy.p -(p/rho)*dWdy.rho);
    return (heat_flux);

}

Vector3D LES3D_Polytropic_pState::q_t_y(const double &kappa_t,
                                        const LES3D_Polytropic_pState &dWdx, 
                                        const LES3D_Polytropic_pState &dWdy,
                                        const LES3D_Polytropic_pState &dWdz,
                                        const double &Volume) {
    Vector3D heat_flux;
    heat_flux.y = -kappa_t*(1.0/(rho*R)) * (dWdy.p -(p/rho)*dWdy.rho);
    return (heat_flux);
}

/**
 * LES3D_Polytropic_pState::q_t_z
 * Return component of LES3D_Polytropic_pState::q_t in the z-direction.   
 */
Vector3D LES3D_Polytropic_pState::q_t_z(const LES3D_Polytropic_pState &dWdx, 
                                        const LES3D_Polytropic_pState &dWdy,
                                        const LES3D_Polytropic_pState &dWdz,
                                        const double &Volume) {
    double kappa_t_ = kappa_t(dWdx,dWdy,dWdz,Volume);
    Vector3D heat_flux;
    heat_flux.z = -kappa_t_*(1.0/(rho*R)) * (dWdz.p -(p/rho)*dWdz.rho);
    return (heat_flux);
}

Vector3D LES3D_Polytropic_pState::q_t_z(const double &kappa_t,
                                        const LES3D_Polytropic_pState &dWdx, 
                                        const LES3D_Polytropic_pState &dWdy,
                                        const LES3D_Polytropic_pState &dWdz,
                                        const double &Volume) {
    Vector3D heat_flux;
    heat_flux.z = -kappa_t*(1.0/(rho*R)) * (dWdz.p -(p/rho)*dWdz.rho);
    return (heat_flux);
}











/*
 * Viscous flux vectors
 * --------------------
 */

/**
 * LES3D_Polytropic_pState::Fv
 * Viscous flux (x-direction) 
 */
LES3D_Polytropic_cState LES3D_Polytropic_pState::Fv(const LES3D_Polytropic_pState &dWdx,
                                                    const LES3D_Polytropic_pState &dWdy,
                                                    const LES3D_Polytropic_pState &dWdz,
                                                    const double &Volume) {
   Tensor3D tau;
   Vector3D q;
   double mu_ = mu();
   double mu_t_ = mu_t(dWdx,dWdy,dWdz,Volume);
   double kappa_ = kappa();
   double kappa_t_ = kappa_t(dWdx,dWdy,dWdz,Volume);
   double k = SFS_Kinetic_Energy(dWdx,dWdy,dWdz,Volume);
 
   /* Calculate "tau + tau_t" by giving both mu and mu_t.
    * Warning: Doing this adds a term "-2/3*rho*k" too much
    *          on the diagonal element! This must be corrected.
    */
   tau = tau_t_x(mu_ + mu_t_, dWdx, dWdy, dWdz, Volume);
   tau.xx += TWO/THREE*rho*k;
   q = q_t_x(kappa_ + kappa_t_, dWdx, dWdy, dWdz, Volume);
   return (LES3D_Polytropic_cState(ZERO,tau.xx,tau.xy,tau.xz,v.x*tau.xx+v.y*tau.xy+v.z*tau.xz-q.x));
}

/**
 * LES3D_Polytropic_pState::Fvx
 * Viscous flux (x-direction) 
 */
LES3D_Polytropic_cState LES3D_Polytropic_pState::Fvx(const LES3D_Polytropic_pState &dWdx,
                                                     const LES3D_Polytropic_pState &dWdy,
                                                     const LES3D_Polytropic_pState &dWdz,
                                                     const double &Volume) {
    Tensor3D tau;
    Vector3D q;
    double mu_ = mu();
    double mu_t_ = mu_t(dWdx,dWdy,dWdz,Volume);
    double kappa_ = kappa();
    double kappa_t_ = kappa_t(dWdx,dWdy,dWdz,Volume);
    double k = SFS_Kinetic_Energy(dWdx,dWdy,dWdz,Volume);

    /* Calculate "tau + tau_t" by giving both mu and mu_t.
     * Warning: Doing this adds a term "-2/3*rho*k" too much
     *          on the diagonal element! This must be corrected.
     */
    tau = tau_t_x(mu_ + mu_t_, dWdx, dWdy, dWdz, Volume);
    tau.xx += TWO/THREE*rho*k;
    q = q_t_x(kappa_ + kappa_t_, dWdx, dWdy, dWdz, Volume);
    return (LES3D_Polytropic_cState(ZERO,tau.xx,tau.xy,tau.xz,v.x*tau.xx+v.y*tau.xy+v.z*tau.xz-q.x));
}

/**
 * LES3D_Polytropic_pState::Fvy
 * Viscous flux (y-direction) 
 */
LES3D_Polytropic_cState LES3D_Polytropic_pState::Fvy(const LES3D_Polytropic_pState &dWdx,
                                                     const LES3D_Polytropic_pState &dWdy,
                                                     const LES3D_Polytropic_pState &dWdz,
                                                     const double &Volume) {
    Tensor3D tau;
    Vector3D q;
    double mu_ = mu();
    double mu_t_ = mu_t(dWdx,dWdy,dWdz,Volume);
    double kappa_ = kappa();
    double kappa_t_ = kappa_t(dWdx,dWdy,dWdz,Volume);
    double k = SFS_Kinetic_Energy(dWdx,dWdy,dWdz,Volume);
    
    /* Calculate "tau + tau_t" by giving both mu and mu_t.
     * Warning: Doing this adds a term "-2/3*rho*k" too much
     *          on the diagonal element! This must be corrected.
     */
    tau = tau_t_y(mu_ + mu_t_, dWdx, dWdy, dWdz, Volume);
    tau.yy += TWO/THREE*rho*k;
    q = q_t_y(kappa_ + kappa_t_, dWdx, dWdy, dWdz, Volume);
    return (LES3D_Polytropic_cState(ZERO,tau.xy,tau.yy,tau.yz,v.x*tau.xy+v.y*tau.yy+v.z*tau.yz-q.y));
}

/**
 * LES3D_Polytropic_pState::Fvz
 * Viscous flux (z-direction) 
 */
LES3D_Polytropic_cState LES3D_Polytropic_pState::Fvz(const LES3D_Polytropic_pState &dWdx,
                                                     const LES3D_Polytropic_pState &dWdy,
                                                     const LES3D_Polytropic_pState &dWdz,
                                                     const double &Volume) {
    Tensor3D tau;
    Vector3D q;
    double mu_ = mu();
    double mu_t_ = mu_t(dWdx,dWdy,dWdz,Volume);
    double kappa_ = kappa();
    double kappa_t_ = kappa_t(dWdx,dWdy,dWdz,Volume);
    double k = SFS_Kinetic_Energy(dWdx,dWdy,dWdz,Volume);
    
    /* Calculate "tau + tau_t" by giving both mu and mu_t.
     * Warning: Doing this adds a term "-2/3*rho*k" too much
     *          on the diagonal element! This must be corrected.
     */
    tau = tau_t_z(mu_ + mu_t_, dWdx, dWdy, dWdz, Volume);
    tau.zz += TWO/THREE*rho*k;
    q = q_t_z(kappa_ + kappa_t_, dWdx, dWdy, dWdz, Volume);
    return (LES3D_Polytropic_cState(ZERO,tau.xz,tau.yz,tau.zz,v.x*tau.xz+v.y*tau.yz+v.z*tau.zz-q.z));
}


// n-direction viscous flux and heat flux
LES3D_Polytropic_cState  LES3D_Polytropic_pState::FluxViscous_n(const LES3D_Polytropic_pState &Wl,
                                                                const LES3D_Polytropic_pState &Wr,
                                                                const LES3D_Polytropic_pState &Wc,
                                                                const LES3D_Polytropic_pState &Wc_Neighbour,
                                                                const LES3D_Polytropic_pState &dWdx,
                                                                const LES3D_Polytropic_pState &dWdy,
                                                                const LES3D_Polytropic_pState &dWdz,
                                                                const LES3D_Polytropic_pState &dWdx_Neighbour,
                                                                const LES3D_Polytropic_pState &dWdy_Neighbour,
                                                                const LES3D_Polytropic_pState &dWdz_Neighbour,
                                                                const Vector3D &norm, 
                                                                const Vector3D &ts, 
                                                                const double &deltad,
                                                                const double &Volume, 
                                                                const double &Volume_Neighbour){
    
	// construct the gradients on the cell interface (surface) 
	// based on Hybrid Average Gradient-Diamond-Path Approach
	// Grad Phi = (Phi_c_neigbor - Phi_c)/ds *norm/(norm . ts)
	//            + [bar (Grad Phi) - bar (Grad Phi) . ts norm/(norm.ts)]
	
	// weighted factor based on volume
	double alpha = Volume/(Volume + Volume_Neighbour);
	
	LES3D_Polytropic_pState dWdx_Weighted, 
    dWdy_Weighted, dWdz_Weighted, dWdx_face, 
    dWdy_face, dWdz_face, Grad_middle_term;
	
	LES3D_Polytropic_pState W_face;
	
	dWdx_Weighted = alpha*dWdx + (1.0 - alpha)*dWdx_Neighbour;
	dWdy_Weighted = alpha*dWdy + (1.0 - alpha)*dWdy_Neighbour;
	dWdz_Weighted = alpha*dWdz + (1.0 - alpha)*dWdz_Neighbour;
	
    // Evaluate a weighted term for solution gradients  
	Grad_middle_term = dWdx_Weighted*ts.x + dWdy_Weighted*ts.y + dWdz_Weighted*ts.z;
	
    // Evaluate gradients of primitive variables on the face
	dWdx_face = (Wc_Neighbour - Wc)/deltad *norm.x/dot(norm, ts) + (dWdx_Weighted -  Grad_middle_term*norm.x/dot(norm, ts));
	dWdy_face = (Wc_Neighbour - Wc)/deltad *norm.y/dot(norm, ts) + (dWdy_Weighted -  Grad_middle_term*norm.y/dot(norm, ts));
	dWdz_face = (Wc_Neighbour - Wc)/deltad *norm.z/dot(norm, ts) + (dWdz_Weighted -  Grad_middle_term*norm.z/dot(norm, ts));
	
    // Determine face solution state.   
    W_face = HALF*(Wl + Wr);
	
    // Evaluate viscous flux
    if (fabs(norm.y) < TOLER && fabs(norm.z) < TOLER) {
        return (W_face.Fvx(dWdx_face, dWdy_face, dWdz_face,Volume)*norm.x);
        
    } else if (fabs(norm.x) < TOLER && fabs(norm.z) < TOLER) {
        return (W_face.Fvy(dWdx_face, dWdy_face, dWdz_face,Volume)*norm.y);
        
    } else if (fabs(norm.x) < TOLER && fabs(norm.y) < TOLER) {
        return (W_face.Fvz(dWdx_face, dWdy_face, dWdz_face,Volume)*norm.z);
        
    } else {
        return (W_face.Fvx(dWdx_face, dWdy_face, dWdz_face,Volume)*norm.x +
                W_face.Fvy(dWdx_face, dWdy_face, dWdz_face,Volume)*norm.y +
                W_face.Fvz(dWdx_face, dWdy_face, dWdz_face,Volume)*norm.z);
    }
}






/*
 * Turbulence Model Source Terms
 * ----------------------------- 
 */

/**
 * LES3D_Polytropic_pState::Enstrophy
 * Returns the enstrophy
 */
double LES3D_Polytropic_pState::Enstrophy(const LES3D_Polytropic_pState &dWdx, 
                                          const LES3D_Polytropic_pState &dWdy, 
                                          const LES3D_Polytropic_pState &dWdz) const{
   return HALF*sqr(vorticity(dWdx,dWdy,dWdz));
}

/**
 * LES3D_Polytropic_pState::strain_rate
 * Returns strain rate tensor.       
 */
Tensor3D LES3D_Polytropic_pState::strain_rate(const LES3D_Polytropic_pState &dWdx, 
                                              const LES3D_Polytropic_pState &dWdy, 
                                              const LES3D_Polytropic_pState &dWdz) const {
    Tensor3D S;
    S.zero();
    S.xx = dWdx.v.x;
    S.yy = dWdy.v.y;
    S.zz = dWdz.v.z;
    S.xy = HALF*(dWdy.v.x + dWdx.v.y);
    S.xz = HALF*(dWdz.v.x + dWdx.v.z);
    S.yz = HALF*(dWdz.v.y + dWdy.v.z);
    
    return (S);
}

/**
 * LES3D_Polytropic_pState::abs_strain_rate
 * Absolute value of strain rate tensor      *
 */
double LES3D_Polytropic_pState::abs_strain_rate(const LES3D_Polytropic_pState &dWdx, 
                                                const LES3D_Polytropic_pState &dWdy, 
                                                const LES3D_Polytropic_pState &dWdz) const{

   Tensor3D S = strain_rate(dWdx, dWdy, dWdz);
   // S[i,j]*S[i,j]
   double SS = sqr(S.xx) + sqr(S.yy) + sqr(S.zz)
               + TWO*(sqr(S.xy) + sqr(S.yz) + sqr(S.xz));
   // sqrt(2*S*S)
   return sqrt(TWO*SS);

}

/**
 * LES3D_Polytropic_pState::vorticity
 * Returns vorticity vector. 
 */
Vector3D LES3D_Polytropic_pState::vorticity(const LES3D_Polytropic_pState &dWdx, 
                                            const LES3D_Polytropic_pState &dWdy, 
                                            const LES3D_Polytropic_pState &dWdz) const {
    Vector3D vorticity_vector;
    
    vorticity_vector.x = dWdy.v.z - dWdz.v.y;
    vorticity_vector.y = -(dWdx.v.z - dWdz.v.x);
    vorticity_vector.z = dWdx.v.y - dWdy.v.x;
    
    return (vorticity_vector);
}

/**
 * LES3D_Polytropic_pState::Q_criterion
 * A measure to view coherent turbulence structures
 */
double LES3D_Polytropic_pState::Q_criterion(const LES3D_Polytropic_pState &dWdx, 
                                            const LES3D_Polytropic_pState &dWdy, 
                                            const LES3D_Polytropic_pState &dWdz) const{
    
    return 0.25*(sqr(vorticity(dWdx,dWdy,dWdz)) - sqr(abs_strain_rate(dWdx,dWdy,dWdz)));
}


/**
 * LES3D_Polytropic_pState::dissipation
 * Subfilter scale kinetic energy
 */
double LES3D_Polytropic_pState::viscous_dissipation(const LES3D_Polytropic_pState &dWdx,
                                                    const LES3D_Polytropic_pState &dWdy,
                                                    const LES3D_Polytropic_pState &dWdz) {

    return nu()*sqr(abs_strain_rate(dWdx,dWdy,dWdz));
    
}

/**
 * LES3D_Polytropic_pState::SFS_dissipation
 * Subfilter scale kinetic energy
 */
double LES3D_Polytropic_pState::SFS_dissipation(const LES3D_Polytropic_pState &dWdx,
                                                const LES3D_Polytropic_pState &dWdy,
                                                const LES3D_Polytropic_pState &dWdz,
                                                const double &Volume) {
    Tensor3D S = strain_rate(dWdx,dWdy,dWdz);
    Tensor3D tau = tau_t(dWdx,dWdy,dWdz,Volume);
    
    // tau[i,j]*S[i,j]
    double tau_x_S = tau.xx*S.xx + tau.yy*S.yy + tau.zz*S.zz
                     + TWO*(tau.xy*S.xy + tau.yz*S.yz + tau.xz*S.xz);
    return tau_x_S;    
}


/**
 * LES3D_Polytropic_pState::SFS_Kinetic_Energy
 * Subfilter scale kinetic energy
 */
double LES3D_Polytropic_pState::SFS_Kinetic_Energy(const LES3D_Polytropic_pState &dWdx,
                                                   const LES3D_Polytropic_pState &dWdy,
                                                   const LES3D_Polytropic_pState &dWdz,
                                                   const double &Volume) {
 //   double CI = 0.005;
//    return (CI*sqr(filter_width(Volume)*abs_strain_rate(dWdx,dWdy,dWdz)));
    
    double Ck = 0.05;  /* Yoshizawa 1985 */
    double Delta = filter_width(Volume);
    return ( sqr(nu_t(dWdx,dWdy,dWdz,Volume)) / sqr(Ck*Delta) );

}
















/*----------------------------------------------------------------------------*
 *                LESS3D_Polytropic_cState member functions                   *
 *----------------------------------------------------------------------------*/



SFS_model_Parameters LES3D_Polytropic_cState::SFS_model = {SFS_MODEL_SMAGORINSKY,
                                                           0.18};
Filter_Parameters LES3D_Polytropic_cState::filter = {2,
                                                     FILTER_TYPE_IMPLICIT};

void LES3D_Polytropic_cState::Set_LES_parameters(SFS_model_Parameters &SFS_model_, Filter_Parameters &filter_){
    SFS_model = SFS_model_;
    filter = filter_;
}


/*
 * Turbulent transport coefficients
 * --------------------------------
 */

/**
 * LES3D_Polytropic_cState::mu_t
 * Return dynamic eddy (turbulent) viscosity using Smagorinsky model
 */
double LES3D_Polytropic_cState::mu_t(const LES3D_Polytropic_pState &dWdx,
                                     const LES3D_Polytropic_pState &dWdy,
                                     const LES3D_Polytropic_pState &dWdz,
                                     const double &Volume) {
    return(rho*nu_t(dWdx,dWdy,dWdz,Volume));
}

/**
 * LES3D_Polytropic_cState::mu_t
 * Return kinematic eddy (turbulent) viscosity using Smagorinsky model
 */
double LES3D_Polytropic_cState::nu_t(const LES3D_Polytropic_pState &dWdx,
                                     const LES3D_Polytropic_pState &dWdy,
                                     const LES3D_Polytropic_pState &dWdz,
                                     const double &Volume) {
    double Delta = filter_width(Volume);
    double Cs = SFS_model.smagorinsky_coefficient;
    return(sqr(Cs*Delta)*abs_strain_rate(dWdx,dWdy,dWdz));
}

/**
 * LES3D_Polytropic_cState::kappa_t
 * Return turbulent thermal conductivity.     
 */
double LES3D_Polytropic_cState::kappa_t(const LES3D_Polytropic_pState &dWdx,
                                        const LES3D_Polytropic_pState &dWdy,
                                        const LES3D_Polytropic_pState &dWdz,
                                        const double &Volume) {
  return (mu_t(dWdx,dWdy,dWdz,Volume)*Cp/Pr_t()); 
}


/**
 * LES3D_Polytropic_cState::Pr_t
 * Return turbulent Prandtl number.
 */
double LES3D_Polytropic_cState::Pr_t(void) {
   return (0.9); 
}


/**
 * LES3D_Polytropic_cState::filter_width
 * Return LES characteristic filter width
 */
double LES3D_Polytropic_cState::filter_width(const double &Volume) const {
  return (filter.FGR*pow(Volume,1.0/3.0)); 
}








/*
 * Turbulence Model Source Terms
 * ----------------------------- 
 */

/**
 * LES3D_Polytropic_cState::strain_rate
 * Returns strain rate tensor.       
 */
Tensor3D LES3D_Polytropic_cState::strain_rate(const LES3D_Polytropic_pState &dWdx, 
                                              const LES3D_Polytropic_pState &dWdy, 
                                              const LES3D_Polytropic_pState &dWdz) const {
    Tensor3D S;
    S.zero();
    S.xx = dWdx.v.x;
    S.yy = dWdy.v.y;
    S.zz = dWdz.v.z;
    S.xy = HALF*(dWdy.v.x + dWdx.v.y);
    S.xz = HALF*(dWdz.v.x + dWdx.v.z);
    S.yz = HALF*(dWdz.v.y + dWdy.v.z);
    
    return (S);
}

/**
 * LES3D_Polytropic_cState::abs_strain_rate
 * Absolute value of strain rate tensor
 */
double LES3D_Polytropic_cState::abs_strain_rate(const LES3D_Polytropic_pState &dWdx, 
                                                const LES3D_Polytropic_pState &dWdy, 
                                                const LES3D_Polytropic_pState &dWdz) const{

    Tensor3D S = strain_rate(dWdx,dWdy,dWdz);
    // S[i,j]*S[i,j]
    double SS = sqr(S.xx) + sqr(S.yy) + sqr(S.zz)
                + TWO*(sqr(S.xy) + sqr(S.yz) + sqr(S.xz));
    // sqrt(2*S*S)
    return sqrt(TWO*SS);
}

/**
 * LES3D_Polytropic_pState::vorticity
 * Returns vorticity vector. 
 */
Vector3D LES3D_Polytropic_cState::vorticity(const LES3D_Polytropic_pState &dWdx, 
                                            const LES3D_Polytropic_pState &dWdy, 
                                            const LES3D_Polytropic_pState &dWdz) const {
    Vector3D vorticity_vector;
    
    vorticity_vector.x = dWdy.v.z - dWdz.v.y;
    vorticity_vector.y = -(dWdx.v.z - dWdz.v.x);
    vorticity_vector.z = dWdx.v.y - dWdy.v.x;
    
    return (vorticity_vector);
}

/**
 * LES3D_Polytropic_cState::Q_criterion
 * A measure for coherent turbulence structures
 */
double LES3D_Polytropic_cState::Q_criterion(const LES3D_Polytropic_pState &dWdx, 
                                            const LES3D_Polytropic_pState &dWdy, 
                                            const LES3D_Polytropic_pState &dWdz) const{
    
    return 0.25*(sqr(vorticity(dWdx,dWdy,dWdz)) - sqr(abs_strain_rate(dWdx,dWdy,dWdz)));
}


/**
 * LES3D_Polytropic_cState::SFS_Kinetic_Energy
 * Subfilter scale kinetic energy
 */
double LES3D_Polytropic_cState::SFS_Kinetic_Energy(const LES3D_Polytropic_pState &dWdx,
                                                   const LES3D_Polytropic_pState &dWdy,
                                                   const LES3D_Polytropic_pState &dWdz,
                                                   const double &Volume) {
//    double CI = 0.005;
//    return (CI*sqr(filter_width(Volume)*abs_strain_rate(dWdx,dWdy,dWdz)));
    double Ck = 0.05;  /* Yoshizawa 1985 */
    double Delta = filter_width(Volume);
    return ( sqr(nu_t(dWdx,dWdy,dWdz,Volume)) / sqr(Ck*Delta) );
}