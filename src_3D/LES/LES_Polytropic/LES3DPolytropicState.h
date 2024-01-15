/*! \file LES3DPolytropicState.h
 * 	\brief	Header file defining the Favre-filtered Navier-Stokes solution state 
 *              classes associated with solution of turbulent 
 *              flows of a polytropic gas using a large eddy simulation
 *              (LES) technique.
 * .
 */

#ifndef _LES3D_POLYTROPIC_STATE_INCLUDED 
#define _LES3D_POLYTROPIC_STATE_INCLUDED

#ifndef _TURBULENCE_MODELLING_INCLUDED
#include "../../TurbulenceModelling/TurbulenceModelling.h"
#endif // TURBULENCE_MODELLING_INCLUDED  

/* Include header file for base solution classes from which the classes are derived. */

#ifndef _NAVIERSTOKES3D_POLYTROPIC_STATE_INCLUDED
#include "../../NavierStokes/NavierStokes3DPolytropicState.h"
#endif  //NAVIERSTOKES3D_POLYTROPIC_STATE_INCLUDED


/* Define the classes. */

class LES3D_Polytropic_cState;
class LES3D_Polytropic_pState;

struct SFS_model_Parameters {
    int model;
    double smagorinsky_coefficient;
};

struct Filter_Parameters {
    double FGR;
    int type;
};


/*!
 * Class: LES3D_Polytropic_pState
 *
 * \brief Primitive state solution class for 3D Favre-filtered Navier-Stokes
 *        equations associated with solution of compressible turbulent 
 *        flows of a polytropic gas using a large eddy simulation
 *        (LES) technique in conjunction with.
 *
 * Member functions
 *  - rho              -- Return density (kg/m^3)
 *  - v                -- Return flow velocity (m/s)
 *  - p                -- Return pressure (Pa, N/m^2)
 *  - p_t              -- Return turbulence modified pressure (Pa, N/m^2)
 *  - k                -- Return subfilter scale turbulent kinetic energy (m^2/s^2)
 *  - C                -- Return reaction rate progress variable (0 <= C <= 1)
 *  - spec             -- Return array of species mass fraction data
 *  - Mass             -- Return mixture molecular mass (kg/mol)
 *  - Rtot             -- Return mixture gas constant (J/(kg*K))
 *  - HeatofFormation  -- Return heat of formation for the mixture
 *  - Cp               -- Return specific heat at constant pressure for mixture (J/(kg*K))
 *  - Cv               -- Return specific heat at constant volume for mixture (J/(kg*K))
 *  - g                -- Return specific heat ratio for mixture
 *  - e                -- Return mixture absolute internal energy (J/kg)
 *  - h                -- Return mixture absolute specific enthalpy (J/kg)
 *  - E                -- Return total mixture energy (J/kg)
 *  - H                -- Return total mixture enthalpy (J/kg)
 *  - rhov             -- Return momentum of mixture (kg/(m^2*s))
 *  - T                -- Return mixture temperature (K)
 *  - a                -- Return sound speed of mixture including turbulent kinetic energy (m/s)
 *  - M                -- Return Mach number for mixture
 *  - mu               -- Return mixture dynamic viscosity
 *  - kappa            -- Return mixture coefficient of thermal conductivity
 *  - Pr               -- Return mixture Prandtl number
 *  - tau              -- Return (molecular) fluid stress tensor
 *  - tau_x            -- Return components of (molecular) fluid stress tensor in x-direction
 *  - tau_y            -- Return components of (molecular) fluid stress tensor in y-direction
 *  - tau_z            -- Return components of (molecular) fluid stress tensor in z-direction
 *  - q                -- Return (molecular) heat flux vector
 *  - q_x              -- Return component of (molecular) heat flux vector in x-direction
 *  - q_y              -- Return component of (molecular) heat flux vector in y-direction
 *  - q_z              -- Return component of (molecular) heat flux vector in z-direction
 *  - mu_t             -- Return subfilter scale turbulent viscosity
 *  - kappa_t          -- Return subfilter scale turbulent thermal conductivity
 *  - Ds_t             -- Return subfilter scale species turbulent mass diffusion coefficient
 *  - Pr_t             -- Return subfilter scale turbulent Prandtl number
 *  - Sc_t             -- Return subfilter scale turbulent Schmidt number
 *  - Le_t             -- Return subfilter scale turbulent Lewis number
 *  - tau_t            -- Return subfilter scale turbulent stress tensor
 *  - tau_t_x          -- Return components of subfilter scale turbulent stress tensor in x-direction
 *  - tau_t_y          -- Return components of subfilter scale turbulent stress tensor in y-direction
 *  - tau_t_z          -- Return components of subfilter scale turbulent stress tensor in z-direction
 *  - q_t              -- Return subfilter scale turbulent heat flux vector
 *  - q_t_x            -- Return component of subfilter scale turbulent heat flux vector in x-direction
 *  - q_t_y            -- Return component of subfilter scale turbulent heat flux vector in y-direction
 *  - q_t_z            -- Return component of subfilter scale turbulent heat flux vector in z-direction
 *  - U                -- Return conserved solution state
 *  - F, Fx            -- Return x-direction inviscid solution flux
 *  - Fy               -- Return y-direction inviscid solution flux
 *  - Fz               -- Return z-direction inviscid solution flux
 *  - Fv, Fvx          -- Return x-direction viscous solution flux
 *  - Fvy              -- Return y-direction viscous solution flux
 *  - Fvz              -- Return z-direction viscous solution flux
 *  - lambda           -- Return x-direction eigenvalue
 *  - rc               -- Return x-direction conserved right eigenvector
 *  - lp               -- Return x-direction primitive left eigenvector
 *  - lambda_x         -- Return x-direction eigenvalue
 *  - rc_x             -- Return x-direction conserved right eigenvector
 *  - lp_x             -- Return x-direction primitive left eigenvector
 *  - RoeAverage       -- Return Roe-average solution state vector
 *  - FluxHLLE_x       -- Return HLLE numerical solution flux in x-direction
 *  - FluxHLLE_n       -- Return HLLE numerical solution flux in n-direction
 *  - FluxRoe_x        -- Return Roe numerical solution flux in x-direction
 *  - FluxRoe_n        -- Return Roe numerical solution flux in n-direction
 *  - lambda_minus     -- Return negative eigenvalues, applying Harten entropy fix
 *  - lambda_plus      -- Return positive eigenvalues, applying Harten entropy fix
 *  - Reflect          -- Return reflected solution state after application of reflection BC
 *  - MovingWall       -- Return wall solution state after application of moving wall BC
 *  - NoSlip           -- Return wall solution state after application of no-slip BC
 *  - FluxViscous_n    -- Returns viscous flux in n-direction
 *
 * Member operators \n
 *      W -- a primitive solution state \n
 *      c -- a scalar (double)
 *
 *  - W = W;
 *  - c = W[i];
 *  - W = W + W;
 *  - W = W - W;
 *  - c = W * W; (inner product)
 *  - W = c * W;
 *  - W = W * c;
 *  - W = W / c;
 *  - W = W ^ W;
 *  - W = +W;
 *  - W = -W;
 *  - W += W;
 *  - W -= W;
 *  - W *= W;
 *  - W /= W;
 *  - W == W;
 *  - W != W;
 *  - cout << W; (output function)
 *  - cin  >> W; (input function)
 *  \nosubgrouping
 */
class LES3D_Polytropic_pState : public NavierStokes3D_Polytropic_pState {
public:

/** @name Constructors and destructors */
/*        ---------------------------- */
//@{
    LES3D_Polytropic_pState() : NavierStokes3D_Polytropic_pState() { }
    LES3D_Polytropic_pState(const double &rho, 
                            const Vector3D &v, 
                            const double &p) : NavierStokes3D_Polytropic_pState(rho,v,p) { }
    LES3D_Polytropic_pState(const double &rho, 
                            const double &vx, const double &vy, const double &vz, 
                            const double &p) : NavierStokes3D_Polytropic_pState(rho,vx,vy,vz,p) { }
    LES3D_Polytropic_pState(const double &value) : NavierStokes3D_Polytropic_pState(value) { }
    LES3D_Polytropic_pState(const LES3D_Polytropic_pState &W) { Copy(W); }

    LES3D_Polytropic_pState(const NavierStokes3D_Polytropic_pState &W) : NavierStokes3D_Polytropic_pState(W) { } 
    LES3D_Polytropic_pState(const Euler3D_Polytropic_pState &W) : NavierStokes3D_Polytropic_pState(W) { } 
//@}

/** @name Operators */
/*        --------- */
//@{    
    LES3D_Polytropic_pState& operator =(const LES3D_Polytropic_pState &W){
        if(this != &W)
            NavierStokes3D_Polytropic_pState::operator=(W);
        return *this;
    }
//@}
   
/** @name SFS model and filter */
/*        -------------------- */
//@{
    static Filter_Parameters filter;
    static SFS_model_Parameters SFS_model;
    void Set_LES_parameters(SFS_model_Parameters &SFS_model_, Filter_Parameters &filter_);
//@}    
    
    
    static double E_return(LES3D_Polytropic_pState &W) {
        return W.E();
    }
    
/** @name Turbulent transport coefficients */
/*        -------------------------------- */
//@{
    //! Subfilter scale turbulent dynamic viscosity
    double mu_t(const LES3D_Polytropic_pState &dWdx, 
                const LES3D_Polytropic_pState &dWdy, 
                const LES3D_Polytropic_pState &dWdz,
                const double &Volume);
    
    //! Subfilter scale turbulent dynamic viscosity adaptor
    double mu_t(const LES3D_Polytropic_pState &dWdx, 
                const LES3D_Polytropic_pState &dWdy, 
                const LES3D_Polytropic_pState &dWdz,
                const int Flow_Type, 
                const double &Volume) {
        return mu_t(dWdx,dWdy,dWdz,Volume);
    }
    
    //! Subfilter scale turbulent kinematic viscosity
    double nu_t(const LES3D_Polytropic_pState &dWdx, 
                const LES3D_Polytropic_pState &dWdy, 
                const LES3D_Polytropic_pState &dWdz,
                const double &Volume);

    //! Turbulent thermal conductivity
    double kappa_t(const LES3D_Polytropic_pState &dWdx, 
                   const LES3D_Polytropic_pState &dWdy, 
                   const LES3D_Polytropic_pState &dWdz,
                   const double &Volume);      

    //! Turbulent Prandtl number
    double Pr_t(void)  { return 0.9; }

    //! LES filter width
    double filter_width(const double &Volume) const { return (filter.FGR*pow(Volume,1.0/3.0)); }
//@}

/** @name Subfilter scale turbulent stress tensor and heat flux vector */
/*        ------------------------------------------------------------ */
//@{
   //! Returns subfilter scale turbulent stress tensor 
   Tensor3D tau_t(const LES3D_Polytropic_pState &dWdx, 
                  const LES3D_Polytropic_pState &dWdy,
                  const LES3D_Polytropic_pState &dWdz,
                  const double &Volume);

   //! Returns subfilter scale turbulent stress tensor 
   Tensor3D tau_t(const double &mu_t_temp, 
                  const LES3D_Polytropic_pState &dWdx, 
                  const LES3D_Polytropic_pState &dWdy,
                  const LES3D_Polytropic_pState &dWdz,
                  const double &Volume);

   //! Returns components of subfilter scale turbulent stress tensor in x-direction 
   Tensor3D tau_t_x(const LES3D_Polytropic_pState &dWdx, 
                    const LES3D_Polytropic_pState &dWdy,
                    const LES3D_Polytropic_pState &dWdz,
                    const double &Volume);

   //! Returns components of subfilter scale turbulent stress tensor in x-direction
   Tensor3D tau_t_x(const double &mu_t_temp, 
                    const LES3D_Polytropic_pState &dWdx, 
                    const LES3D_Polytropic_pState &dWdy,
                    const LES3D_Polytropic_pState &dWdz,
                    const double &Volume);

   //! Returns components of subfilter scale turbulent stress tensor in y-direction 
   Tensor3D tau_t_y(const LES3D_Polytropic_pState &dWdx, 
                    const LES3D_Polytropic_pState &dWdy,
                    const LES3D_Polytropic_pState &dWdz,
                    const double &Volume);

   //! Returns components of subfilter scale turbulent stress tensor in y-direction
   Tensor3D tau_t_y(const double &mu_t_temp, 
                    const LES3D_Polytropic_pState &dWdx, 
                    const LES3D_Polytropic_pState &dWdy,
                    const LES3D_Polytropic_pState &dWdz,
                    const double &Volume);

   //! Returns components of subfilter scale turbulent stress tensor in z-direction 
   Tensor3D tau_t_z(const LES3D_Polytropic_pState &dWdx, 
                    const LES3D_Polytropic_pState &dWdy,
                    const LES3D_Polytropic_pState &dWdz,
                    const double &Volume);

   //! Returns components of subfilter scale turbulent stress tensor in z-direction
   Tensor3D tau_t_z(const double &mu_t_temp, 
                    const LES3D_Polytropic_pState &dWdx, 
                    const LES3D_Polytropic_pState &dWdy,
                    const LES3D_Polytropic_pState &dWdz,
                    const double &Volume);

   //! Returns subfilter scale turbulent heat flux vector 
   Vector3D q_t(const LES3D_Polytropic_pState &dWdx, 
                const LES3D_Polytropic_pState &dWdy,
                const LES3D_Polytropic_pState &dWdz,
                const double &Volume);

   //! Returns subfilter scale turbulent heat flux vector 
   Vector3D q_t(const double &kappa_t_temp,
                const LES3D_Polytropic_pState &dWdx, 
                const LES3D_Polytropic_pState &dWdy,
                const LES3D_Polytropic_pState &dWdz,
                const double &Volume);

   //! Returns component of subfilter scale turbulent heat flux vector in x-direction
   Vector3D q_t_x(const LES3D_Polytropic_pState &dWdx, 
                  const LES3D_Polytropic_pState &dWdy,
                  const LES3D_Polytropic_pState &dWdz,
                  const double &Volume);

   //! Returns component of subfilter scale turbulent heat flux vector in x-direction
   Vector3D q_t_x(const double &kappa_t_temp,
                  const LES3D_Polytropic_pState &dWdx, 
                  const LES3D_Polytropic_pState &dWdy,
                  const LES3D_Polytropic_pState &dWdz,
                  const double &Volume);

   //! Returns component of subfilter scale turbulent heat flux vector in y-direction
   Vector3D q_t_y(const LES3D_Polytropic_pState &dWdx, 
                  const LES3D_Polytropic_pState &dWdy,
                  const LES3D_Polytropic_pState &dWdz,
                  const double &Volume);

   //! Returns component of subfilter scale turbulent heat flux vector in y-direction
   Vector3D q_t_y(const double &kappa_t_temp,
                  const LES3D_Polytropic_pState &dWdx, 
                  const LES3D_Polytropic_pState &dWdy,
                  const LES3D_Polytropic_pState &dWdz,
                  const double &Volume);

   //! Returns component of subfilter scale turbulent heat flux vector in z-direction
   Vector3D q_t_z(const LES3D_Polytropic_pState &dWdx, 
                  const LES3D_Polytropic_pState &dWdy,
                  const LES3D_Polytropic_pState &dWdz,
                  const double &Volume);

   //! Returns component of subfilter scale turbulent heat flux vector in z-direction
   Vector3D q_t_z(const double &kappa_t_temp,
                  const LES3D_Polytropic_pState &dWdx, 
                  const LES3D_Polytropic_pState &dWdy,
                  const LES3D_Polytropic_pState &dWdz,
                  const double &Volume);
//@}
    

/** @name Viscous flux vectors */
/*        --------------------- */
//@{
   //! x-direction viscous solution flux + SFS flux
   LES3D_Polytropic_cState Fv(const LES3D_Polytropic_pState &dWdx,
                              const LES3D_Polytropic_pState &dWdy,
                              const LES3D_Polytropic_pState &dWdz,
                              const double &Volume);

   //! x-direction viscous solution flux  + SFS flux
   LES3D_Polytropic_cState Fvx(const LES3D_Polytropic_pState &dWdx,
                               const LES3D_Polytropic_pState &dWdy,
                               const LES3D_Polytropic_pState &dWdz,
                               const double &Volume);

   //! y-direction viscous solution flux  + SFS flux
   LES3D_Polytropic_cState Fvy(const LES3D_Polytropic_pState &dWdx,
                               const LES3D_Polytropic_pState &dWdy,
                               const LES3D_Polytropic_pState &dWdz,
                               const double &Volume);

   //! z-direction viscous solution flux + SFS flux
   LES3D_Polytropic_cState Fvz(const LES3D_Polytropic_pState &dWdx,
                               const LES3D_Polytropic_pState &dWdy,
                               const LES3D_Polytropic_pState &dWdz,
                               const double &Volume);
    
    static LES3D_Polytropic_cState FluxViscous_n(const LES3D_Polytropic_pState &Wl,
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
                                          const double &Volume_Neighbour);
    
    static LES3D_Polytropic_cState FluxViscous_HighOrder_n(LES3D_Polytropic_pState &W_face,
                                                           const LES3D_Polytropic_pState &dWdx_face,
                                                           const LES3D_Polytropic_pState &dWdy_face,
                                                           const LES3D_Polytropic_pState &dWdz_face,
                                                           const double &Volume,
                                                           const Vector3D &norm);
                
//@}


/** @name Turbulence Model Source Terms */
/*        ----------------------------- */
//@{
    //! Enstrophy
    double Enstrophy(const LES3D_Polytropic_pState &dWdx, 
                     const LES3D_Polytropic_pState &dWdy, 
                     const LES3D_Polytropic_pState &dWdz) const;
    
    //! Returns strain rate tensor.       
    Tensor3D strain_rate(const LES3D_Polytropic_pState &dWdx, 
                         const LES3D_Polytropic_pState &dWdy, 
                         const LES3D_Polytropic_pState &dWdz) const;
    
    //! Absolute value of strain rate
    double abs_strain_rate(const LES3D_Polytropic_pState &dWdx, 
                           const LES3D_Polytropic_pState &dWdy, 
                           const LES3D_Polytropic_pState &dWdz) const;
    
    //! Returns vorticity vector. 
    Vector3D vorticity(const LES3D_Polytropic_pState &dWdx, 
                       const LES3D_Polytropic_pState &dWdy, 
                       const LES3D_Polytropic_pState &dWdz) const;  
    
    //! Q_criterion for turbulence structures
    double Q_criterion(const LES3D_Polytropic_pState &dWdx, 
                       const LES3D_Polytropic_pState &dWdy, 
                       const LES3D_Polytropic_pState &dWdz) const;
    
    //! Dissipation based on strain rate (model)
    double viscous_dissipation(const LES3D_Polytropic_pState &dWdx,
                               const LES3D_Polytropic_pState &dWdy,
                               const LES3D_Polytropic_pState &dWdz);
    
    //! Subfilter scale dissipation
    double SFS_dissipation(const LES3D_Polytropic_pState &dWdx,
                           const LES3D_Polytropic_pState &dWdy,
                           const LES3D_Polytropic_pState &dWdz,
                           const double &Volume);
    
    //! Subfilter scale kinetic energy
    double SFS_Kinetic_Energy(const LES3D_Polytropic_pState &dWdx,
                              const LES3D_Polytropic_pState &dWdy,
                              const LES3D_Polytropic_pState &dWdz,
                              const double &Volume);
    
    //! Subfilter scale kinetic energy adaptor
    double SFS_Kinetic_Energy(const LES3D_Polytropic_pState &dWdx,
                              const LES3D_Polytropic_pState &dWdy,
                              const LES3D_Polytropic_pState &dWdz,
                              const int Flow_Type,
                              const double &Volume) {
        return SFS_Kinetic_Energy(dWdx,dWdy,dWdz,Volume);
    }

    
 //@}

};
  
/*!
 * Class: LES3D_Polytropic_cState
 *
 * \brief Conserved state solution class for 3D Favre-filtered Navier-Stokes
 *        equations associated with solution of premixed compressible turbulent 
 *        combusting flows of a thermally perfect using a large eddy simulation
 *        (LES) technique in conjunction with a flame surface density (FSD)
 *        subfilter scale model.
 *
 * Member functions
 *  - rho      -- Return mixture density (kg/m^3)
 *  - rhov     -- Return mixture momentum (kg/(m^2-s))   
 *  - E        -- Return mixture total energy (J/kg)
 *  - rhoC     -- Return density of reaction rate progress variable (0 <= C <= 1)
 *  - rhoFsD   -- Return premixed flame surface density (FSD)
 *  - rhok     -- Return turbulent kinetic energy (kg/m-s^2)
 *  - rhospec  -- Return array of species density data
 *
 * Member operators \n
 *      U -- a conserved solution state \n
 *      c -- a scalar (double)
 *
 *  - U = U;
 *  - c = U[i];
 *  - U = U + U;
 *  - U = U - U;
 *  - c = U * U; (inner product)
 *  - U = c * U;
 *  - U = U * c;
 *  - U = U / c;
 *  - U = U ^ U;
 *  - U = +U;
 *  - U = -U;
 *  - U += U;
 *  - U -= U;
 *  - U *= U;
 *  - U /= U;
 *  - U == U;
 *  - U != U;
 *  - cout << U; (output function)
 *  - cin  >> U; (input function)
 *  \nosubgrouping
 */
class LES3D_Polytropic_cState : public NavierStokes3D_Polytropic_cState {
public:
    
/** @name Constructors and destructors */
/*        ---------------------------- */
//@{
    LES3D_Polytropic_cState() : NavierStokes3D_Polytropic_cState() { }
    LES3D_Polytropic_cState(const double &rho, 
                            const Vector3D &rhov, 
                            const double &E) : NavierStokes3D_Polytropic_cState(rho,rhov,E) { }
    LES3D_Polytropic_cState(const double &rho, 
                            const double &rhovx, const double &rhovy, const double &rhovz, 
                            const double &E) : NavierStokes3D_Polytropic_cState(rho,rhovx,rhovy,rhovz,E) { }
    LES3D_Polytropic_cState(const double &value) : NavierStokes3D_Polytropic_cState(value) { }
    LES3D_Polytropic_cState(const LES3D_Polytropic_cState &U) {Copy(U); }
    
    LES3D_Polytropic_cState(const NavierStokes3D_Polytropic_cState &U) : NavierStokes3D_Polytropic_cState(U) { }
    LES3D_Polytropic_cState(const Euler3D_Polytropic_cState &U) : NavierStokes3D_Polytropic_cState(U) { } 
//@}

/** @name Operators */
/*        --------- */
//@{    
    LES3D_Polytropic_cState& operator =(const LES3D_Polytropic_cState &U){
        if(this != &U)
            NavierStokes3D_Polytropic_cState::operator=(U);
        return *this;
    }
//@}
    
/** @name SFS model and filter */
/*        -------------------- */
//@{
    static Filter_Parameters filter;
    static SFS_model_Parameters SFS_model;
    void Set_LES_parameters(SFS_model_Parameters &SFS_model_, Filter_Parameters &filter_);
    //@}     
   
/** @name Turbulent transport coefficients */
/*        -------------------------------- */
//@{
    //! Eddy (turbulent) viscosity dynamic
    double mu_t(const LES3D_Polytropic_pState &dWdx, 
                const LES3D_Polytropic_pState &dWdy, 
                const LES3D_Polytropic_pState &dWdz,
                const double &Volume);
    
    //! Eddy (turbulent) viscosity kinematic
    double nu_t(const LES3D_Polytropic_pState &dWdx, 
                const LES3D_Polytropic_pState &dWdy, 
                const LES3D_Polytropic_pState &dWdz,
                const double &Volume);

    //! Turbulent) thermal conductivity
    double kappa_t(const LES3D_Polytropic_pState &dWdx, 
                   const LES3D_Polytropic_pState &dWdy, 
                   const LES3D_Polytropic_pState &dWdz,
                   const double &Volume);     

    //! Turbulent Prandtl number
    double Pr_t(void) { return 0.9; }

    //! LES filter width
    double filter_width(const double &Volume) const { return (filter.FGR*pow(Volume,1.0/3.0)); }

//@}


/** @name Turbulence Model Source Terms */
/*        ----------------------------- */
//@{
    //! Returns strain rate tensor.       
    Tensor3D strain_rate(const LES3D_Polytropic_pState &dWdx, 
                         const LES3D_Polytropic_pState &dWdy, 
                         const LES3D_Polytropic_pState &dWdz) const;
    
    //! Absolute value of strain rate
    double abs_strain_rate(const LES3D_Polytropic_pState &dWdx, 
                          const LES3D_Polytropic_pState &dWdy, 
                          const LES3D_Polytropic_pState &dWdz) const;
    
    //! Returns vorticity vector. 
    Vector3D vorticity(const LES3D_Polytropic_pState &dWdx, 
                       const LES3D_Polytropic_pState &dWdy, 
                       const LES3D_Polytropic_pState &dWdz) const;  
    
    //! Q_criterion for turbulence structures
    double Q_criterion(const LES3D_Polytropic_pState &dWdx, 
                       const LES3D_Polytropic_pState &dWdy, 
                       const LES3D_Polytropic_pState &dWdz) const;
        
    //! Subfilter scale kinetic energy    
    double SFS_Kinetic_Energy(const LES3D_Polytropic_pState &dWdx,
                              const LES3D_Polytropic_pState &dWdy,
                              const LES3D_Polytropic_pState &dWdz,
                              const double &Volume);
//@} 
};

#endif // _LES3DFSD_STATE_INCLUDED
