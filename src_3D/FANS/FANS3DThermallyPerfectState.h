/*! \file FANS3DThermallyPerfectState.h
 * 	\brief	Header file defining the Favre-averaged Navier-Stokes solution state 
 *              classes associated with solution of compressible turbulent flows 
 *              of a thermally perfect non-reactive or combusting mixture.
 */

#ifndef _FANS3D_THERMALLYPERFECT_STATE_INCLUDED 
#define _FANS3D_THERMALLYPERFECT_STATE_INCLUDED

/* Include header file for base solution classes from which the classes are derived. */

#ifndef _NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED
#include "../NavierStokes/NavierStokes3DThermallyPerfectState.h"
#endif  //NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED

/* Include turbulence modelling header file. */

#ifndef _TURBULENCE_MODELLING_INCLUDED
#include "../TurbulenceModelling/TurbulenceModelling.h"
#endif  //TURBULENCE_MODELLING_INCLUDED

/* Define the classes. */

class FANS3D_ThermallyPerfect_KOmega_pState;
class FANS3D_ThermallyPerfect_KOmega_cState;

#define NUM_FANS3D_VAR_EXTRA 2  //rho, v(3), p

/*!
 * Class: FANS3D_ThermallyPerfect_KOmega_pState
 *
 * \brief Primitive state solution class for 3D Favre-average Navier-Stokes
 *        equations with k-omega turbulence model governing flows of thermally 
 *        perfect non-reactive and combusting mixtures.
 *
 * Member functions
 *  - rho           -- Return density (kg/m^3)
 *  - v             -- Return flow velocity (m/s)
 *  - p             -- Return pressure (Pa, N/m^2)
 *  - k             -- Return turbulent kinetic energy (m^2/s^2)
 *  - omega         -- Return specific dissipation rate for turbulent kinetic energy (1/s)
 *  - spec          -- Return array of species mass fraction data
 *  - Mass          -- Return mixture molecular mass (kg/mol)
 *  - Rtot          -- Return mixture gas constant (J/(kg*K))
 *  - Cp            -- Return specific heat at constant pressure for mixture (J/(kg*K))
 *  - Cv            -- Return specific heat at constant volume for mixture (J/(kg*K))
 *  - g             -- Return specific heat ratio for mixture
 *  - e             -- Return mixture absolute internal energy (J/kg)
 *  - es            -- Return mixture sensible internal energy (J/kg)
 *  - h             -- Return mixture absolute specific enthalpy (J/kg)
 *  - hs            -- Return mixture sensible specific enthalpy (J/kg)
 *  - E             -- Return total mixture energy (J/kg)
 *  - H             -- Return total mixture enthalpy (J/kg)
 *  - rhov          -- Return momentum of mixture (kg/(m^2*s))
 *  - T             -- Return mixture temperature (K)
 *  - a             -- Return sound speed of mixture (m/s)
 *  - a_t           -- Return sound speed of mixture including turbulent kinetic energy (m/s)
 *  - M             -- Return Mach number for mixture
 *  - Gibbs         -- Return species Gibbs free energy
 *  - mu            -- Return mixture dynamic viscosity
 *  - kappa         -- Return mixture coefficient of thermal conductivity
 *  - Ds            -- Return species mass diffusion coefficient
 *  - Pr            -- Return mixture Prandtl number
 *  - Sc            -- Return species Schmidt number
 *  - Le            -- Return species Lewis number
 *  - tau           -- Return (molecular) fluid stress tensor
 *  - tau_x         -- Return components of (molecular) fluid stress tensor in x-direction
 *  - tau_y         -- Return components of (molecular) fluid stress tensor in y-direction
 *  - tau_z         -- Return components of (molecular) fluid stress tensor in z-direction
 *  - q             -- Return (molecular) heat flux vector
 *  - q_x           -- Return component of (molecular) heat flux vector in x-direction
 *  - q_y           -- Return component of (molecular) heat flux vector in y-direction
 *  - q_z           -- Return component of (molecular) heat flux vector in z-direction
 *  - mu_t          -- Return eddy (turbulent) viscosity
 *  - kappa_t       -- Return turbulent thermal conductivity
 *  - Ds_t          -- Return species turbulent mass diffusion coefficient
 *  - Pr_t          -- Return turbulent Prandtl number
 *  - Sc_t          -- Return turbulent Schmidt number
 *  - Le_t          -- Return turbulent Lewis number
 *  - tau_t         -- Return Reynolds (turbulent) stress tensor
 *  - tau_t_x       -- Return components of Reynolds (turbulent) stress tensor in x-direction
 *  - tau_t_y       -- Return components of Reynolds (turbulent) stress tensor in y-direction
 *  - tau_t_z       -- Return components of Reynolds (turbulent) stress tensor in z-direction
 *  - q_t           -- Return turbulent heat flux vector
 *  - q_t_x         -- Return component of turbulent heat flux vector in x-direction
 *  - q_t_y         -- Return component of turbulent heat flux vector in y-direction
 *  - q_t_z         -- Return component of turbulent heat flux vector in z-direction
 *  - U             -- Return conserved solution state
 *  - F, Fx         -- Return x-direction inviscid solution flux
 *  - Fy            -- Return y-direction inviscid solution flux
 *  - Fz            -- Return z-direction inviscid solution flux
 *  - Fv, Fvx       -- Return x-direction viscous solution flux
 *  - Fvy           -- Return y-direction viscous solution flux
 *  - Fvz           -- Return z-direction viscous solution flux
 *  - Schemistry    -- Return source terms associated with finite-rate chemistry
 *  - Sturbulence   -- Return source terms associated with turbulence modelling
 *  - lambda        -- Return x-direction eigenvalue
 *  - rc            -- Return x-direction conserved right eigenvector
 *  - lp            -- Return x-direction primitive left eigenvector
 *  - lambda_x      -- Return x-direction eigenvalue
 *  - rc_x          -- Return x-direction conserved right eigenvector
 *  - lp_x          -- Return x-direction primitive left eigenvector
 *  - RoeAverage    -- Return Roe-average solution state vector
 *  - FluxHLLE_x    -- Return HLLE numerical solution flux in x-direction
 *  - FluxHLLE_n    -- Return HLLE numerical solution flux in n-direction
 *  - FluxRoe_x     -- Return Roe numerical solution flux in x-direction
 *  - FluxRoe_n     -- Return Roe numerical solution flux in n-direction
 *  - lambda_minus -- Return negative eigenvalues, applying Harten entropy fix
 *  - lambda_plus  -- Return positive eigenvalues, applying Harten entropy fix
 *  - Reflect       -- Return reflected solution state after application of reflection BC
 *  - MovingWall    -- Return wall solution state after application of moving wall BC
 *  - NoSlip        -- Return wall solution state after application of no-slip BC
 *  - FluxViscous_n -- Returns viscous flux in n-direction
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
class FANS3D_ThermallyPerfect_KOmega_pState : public NavierStokes3D_ThermallyPerfect_pState {
  public:
   double                      k;             //!< Turbulent kinetic energy (m^2/s^2)
   double                      omega;         //!< Specific dissipation rate for turbulent energy (1/s)
   Turbulence_Model_k_omega    k_omega_model; //!< k-omega turbulence model data

/** @name Constructors and destructors */
/*        ---------------------------- */
//@{
   //! Default creation constructor (assign default values)
   FANS3D_ThermallyPerfect_KOmega_pState(void): NavierStokes3D_ThermallyPerfect_pState(), 
     k(ZERO), omega(MILLION) { }

   //! Constructor from base class (allows return of derived type)
   FANS3D_ThermallyPerfect_KOmega_pState(const NavierStokes3D_ThermallyPerfect_pState &W1) : 
     NavierStokes3D_ThermallyPerfect_pState(W1), k(ZERO), omega(MILLION) { }

   //! Constructor from base class
   FANS3D_ThermallyPerfect_KOmega_pState(const NavierStokes3D_ThermallyPerfect_pState &W1,
                                         const double &ke, 
                                         const double &Omega) : 
     NavierStokes3D_ThermallyPerfect_pState(W1), k(ke), omega(Omega) { }

   //! Assignment constructor
   FANS3D_ThermallyPerfect_KOmega_pState(const double &value) :
     NavierStokes3D_ThermallyPerfect_pState(value), k(value), omega(value) { }
   
   //! Assignment constructor
   FANS3D_ThermallyPerfect_KOmega_pState(const double &d, 
                                         const Vector3D &V,
                                         const double &pre, 
                                         const double &ke, 
                                         const double &Omega) :
     NavierStokes3D_ThermallyPerfect_pState(d, V, pre), k(ke), omega(Omega) { }

   //! Assignment constructor
   FANS3D_ThermallyPerfect_KOmega_pState(const double &d, 
                                         const Vector3D &V,
                                         const double &pre, 
                                         const double &ke, 
                                         const double &Omega, 
                                         const double &frac) :
     NavierStokes3D_ThermallyPerfect_pState(d, V, pre, frac), k(ke), omega(Omega) { }
   
   //! Assignment constructor
   FANS3D_ThermallyPerfect_KOmega_pState(const double &d, 
                                         const Vector3D &V,
                                         const double &pre, 
                                         const double &ke, 
                                         const double &Omega, 
                                         const Species *mfrac) :
     NavierStokes3D_ThermallyPerfect_pState(d, V, pre, mfrac), k(ke), omega(Omega) { }

   //! Assignment constructor
   FANS3D_ThermallyPerfect_KOmega_pState(const double &d, 
                                         const Vector3D &V,
                                         const double &pre, 
                                         const double *mfrac) :
     NavierStokes3D_ThermallyPerfect_pState(d, V, pre, mfrac) {
     k = ONE; omega = MILLION; set_initial_values(mfrac); 
   }

   //! Assignment constructor
   FANS3D_ThermallyPerfect_KOmega_pState(const double &d, 
                                         const Vector3D &V,
                                         const double &pre, 
                                         const double &ke, 
                                         const double &Omega, 
                                         const double *mfrac) :
     NavierStokes3D_ThermallyPerfect_pState(d, V, pre, mfrac), k(ke), omega(Omega) {
     set_initial_values(mfrac); 
   }

   //! Assignment constructor
   FANS3D_ThermallyPerfect_KOmega_pState(const double &d, 
                                         const double &vx,
                                         const double &vy, 
                                         const double &vz,
                                         const double &ke, 
                                         const double &Omega,
                                         const double &pre) :
     NavierStokes3D_ThermallyPerfect_pState(d, vx, vy, vz, pre), k(ke), omega(Omega) { }

   //! Assignment constructor
   FANS3D_ThermallyPerfect_KOmega_pState(const double &d, 
                                         const double &vx,
                                         const double &vy, 
                                         const double &vz,
                                         const double &pre, 
                                         const double &ke, 
                                         const double &Omega,
                                         const double &frac) :
     NavierStokes3D_ThermallyPerfect_pState(d, vx, vy, vz, pre, frac), k(ke), omega(Omega) { }
   
   //! Assignment constructor
   FANS3D_ThermallyPerfect_KOmega_pState(const double &d, 
                                         const double &vx,
                                         const double &vy, 
                                         const double &vz,
                                         const double &pre, 
                                         const double &ke, 
                                         const double &Omega, 
                                         const Species *mfrac) :
     NavierStokes3D_ThermallyPerfect_pState(d, vx, vy, vz, pre, mfrac), k(ke), omega(Omega) { }
      
   //! Copy constructor (this is needed for the operator overload returns)
   FANS3D_ThermallyPerfect_KOmega_pState(const FANS3D_ThermallyPerfect_KOmega_pState &W) {
     species_null(); rho = DENSITY_STDATM; set_initial_values(); Copy(W);
   }

   //! Default destructor
   ~FANS3D_ThermallyPerfect_KOmega_pState() {
      Deallocate();
   }
//@}

/** @name Some useful operators */
/*        --------------------- */
//@{
   //! Copy solution state (cheaper than = operator)
   void Copy(const FANS3D_ThermallyPerfect_KOmega_pState &W) {
      Euler3D_ThermallyPerfect_pState::Copy(W); k = W.k; omega = W.omega;
   }

   //! Assigns a vacuum solution state
   void Vacuum(void) {
      Euler3D_ThermallyPerfect_pState::Vacuum(); k = ZERO; omega = ZERO;
   }

   //! Check for physical validity of the solution vector
   bool Realizable_Solution_Check(void);
//@}

/** @name Thermodynamic and other state functions */
/*        --------------------------------------- */
//@{
   //! Total energy of mixture
   double E(void) const;

   //! Total enthalpy of mixture
   double H(void) const;

   //! Total sensible enthalpy of mixture
   double Hs(void) const;

   //! Mixture sound speed including turbulent kinetic energy
   double a_t(void);                
   //! Mixture sound speed including turbulent kinetic energy
   double a_t(void) const;
//@}
   
/** @name Turbulent transport coefficients */
/*        -------------------------------- */
//@{
   //! Eddy (turbulent) viscosity
   double mu_t(void);      

   //! Turbulent) thermal conductivity
   double kappa_t(void);      

   //! Species turbulent diffusion coefficient
   double Ds_t(const int &i);
   //! Species turbulent diffusion coefficient
   double Ds_t(const int &i,
               const double &mu_t_temp);

   //! Turbulent Prandtl number
   double Pr_t(void);      

   //! Turbulent Schmidt number
   double Sc_t(void);      

   //! Turbulent Lewis number
   double Le_t(void); 
//@}
  
/** @name Turbulent (Reynolds) stress tensor and turbulent heat flux vector */
/*        ----------------------------------------------------------------- */
//@{
   //! Returns (Reynolds) turbulent stress tensor 
   Tensor3D tau_t(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                  const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                  const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);

   //! Returns (Reynolds) turbulent stress tensor 
   Tensor3D tau_t(const double &mu_t_temp, 
                  const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                  const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                  const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);

   //! Returns components of (Reynolds) turbulent stress tensor in x-direction 
   Tensor3D tau_t_x(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                    const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                    const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);

   //! Returns components of (Reynolds) turbulent stress tensor in x-direction
   Tensor3D tau_t_x(const double &mu_t_temp, 
                    const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                    const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                    const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);

   //! Returns components of (Reynolds) turbulent stress tensor in y-direction 
   Tensor3D tau_t_y(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                    const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                    const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);

   //! Returns components of (Reynolds) turbulent stress tensor in y-direction
   Tensor3D tau_t_y(const double &mu_t_temp, 
                    const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                    const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                    const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);

   //! Returns components of (Reynolds) turbulent stress tensor in z-direction 
   Tensor3D tau_t_z(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                    const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                    const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);

   //! Returns components of (Reynolds) turbulent stress tensor in z-direction
   Tensor3D tau_t_z(const double &mu_t_temp, 
                    const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                    const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                    const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);

   //! Returns turbulent heat flux vector 
   Vector3D q_t(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);

   //! Returns turbulent heat flux vector 
   Vector3D q_t(const double &kappa_t_temp,
                const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);

   //! Returns component of turbulent heat flux vector in x-direction
   Vector3D q_t_x(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                  const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                  const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);

   //! Returns component of turbulent heat flux vector in x-direction
   Vector3D q_t_x(const double &kappa_t_temp,
                  const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                  const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                  const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);

   //! Returns component of turbulent heat flux vector in y-direction
   Vector3D q_t_y(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                  const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                  const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);

   //! Returns component of turbulent heat flux vector in y-direction
   Vector3D q_t_y(const double &kappa_t_temp,
                  const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                  const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                  const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);

   //! Returns component of turbulent heat flux vector in z-direction
   Vector3D q_t_z(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                  const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                  const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);

   //! Returns component of turbulent heat flux vector in z-direction
   Vector3D q_t_z(const double &kappa_t_temp,
                  const FANS3D_ThermallyPerfect_KOmega_pState &dWdx, 
                  const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                  const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);
//@}

/** @name Conserved solution state */ 
/*        ------------------------ */
//@{
   //! Returns conserved solution state
   FANS3D_ThermallyPerfect_KOmega_cState U(void);
   //! Returns conserved solution state
   FANS3D_ThermallyPerfect_KOmega_cState U(void) const;
   //! Returns conserved solution state
   FANS3D_ThermallyPerfect_KOmega_cState U(const FANS3D_ThermallyPerfect_KOmega_pState &W);  
//@}
   
/** @name Inviscid Flux Vectors */
/*        --------------------- */
//@{
   //! x-direction inviscid solution flux
   FANS3D_ThermallyPerfect_KOmega_cState F(void) ;
   //! x-direction inviscid solution flux
   FANS3D_ThermallyPerfect_KOmega_cState F(void) const ;

   //! x-direction inviscid solution flux
   FANS3D_ThermallyPerfect_KOmega_cState Fx(void);
   //! x-direction inviscid solution flux
   FANS3D_ThermallyPerfect_KOmega_cState Fx(void) const;

   //! y-direction inviscid solution flux
   FANS3D_ThermallyPerfect_KOmega_cState Fy(void);
   //! y-direction inviscid solution flux
   FANS3D_ThermallyPerfect_KOmega_cState Fy(void) const;

   //! z-direction inviscid solution flux
   FANS3D_ThermallyPerfect_KOmega_cState Fz(void);
   //! z-direction inviscid solution flux
   FANS3D_ThermallyPerfect_KOmega_cState Fz(void) const;
//@}

/** @name Viscous flux vectors */
/*        --------------------- */
//@{
   //! x-direction viscous solution flux
   FANS3D_ThermallyPerfect_KOmega_cState Fv(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx,
                                            const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                                            const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);

   //! x-direction viscous solution flux
   FANS3D_ThermallyPerfect_KOmega_cState Fvx(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx,
                                             const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                                             const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);

   //! y-direction viscous solution flux
   FANS3D_ThermallyPerfect_KOmega_cState Fvy(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx,
                                             const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                                             const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);

   //! z-direction viscous solution flux
   FANS3D_ThermallyPerfect_KOmega_cState Fvz(const FANS3D_ThermallyPerfect_KOmega_pState &dWdx,
                                             const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                                             const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);
//@}

/** @name Eigenvalue(s) and eigenvector(s) (x-direction) */
/*        ---------------------------------------------- */
//@{
   //! x-direction eigenvalues
   FANS3D_ThermallyPerfect_KOmega_pState lambda(void);
   //! x-direction eigenvalues
   FANS3D_ThermallyPerfect_KOmega_pState lambda(void) const;

   //! x-direction eigenvalues
   FANS3D_ThermallyPerfect_KOmega_pState lambda_x(void);
   //! x-direction eigenvalues
   FANS3D_ThermallyPerfect_KOmega_pState lambda_x(void) const;

   //! x-direction conservative eigenvectors
   FANS3D_ThermallyPerfect_KOmega_cState rc(const int &index) ;
   //! x-direction conservative eigenvectors
   FANS3D_ThermallyPerfect_KOmega_cState rc(const int &index) const;

   //! x-direction conservative eigenvectors
   FANS3D_ThermallyPerfect_KOmega_cState rc_x(const int &index) ;
   //! x-direction conservative eigenvectors
   FANS3D_ThermallyPerfect_KOmega_cState rc_x(const int &index) const;

   //! x-direction primitive eignenvectors
   FANS3D_ThermallyPerfect_KOmega_pState lp(const int &index) ;
   //! x-direction primitive eignenvectors
   FANS3D_ThermallyPerfect_KOmega_pState lp(const int &index) const;

   //! x-direction primitive eignenvectors
   FANS3D_ThermallyPerfect_KOmega_pState lp_x(const int &index) ;
   //! x-direction primitive eignenvectors
   FANS3D_ThermallyPerfect_KOmega_pState lp_x(const int &index) const;  
//@}

/** @name Numerical Flux Functions */
/*        ------------------------ */
//@{
   //! Returns Roe-averaged primitive solution state
   static FANS3D_ThermallyPerfect_KOmega_pState RoeAverage(const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
                                                           const FANS3D_ThermallyPerfect_KOmega_pState &Wr);

   //! Returns Harten-Lax-van-Leer-Einfeldt (HLLE) x-direction flux
   static FANS3D_ThermallyPerfect_KOmega_cState FluxHLLE_x(const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
                                                           const FANS3D_ThermallyPerfect_KOmega_pState &Wr);

   //! Returns Harten-Lax-van-Leer-Einfeldt (HLLE) x-direction flux
   static FANS3D_ThermallyPerfect_KOmega_cState FluxHLLE_x(const FANS3D_ThermallyPerfect_KOmega_cState &Ul,
                                                           const FANS3D_ThermallyPerfect_KOmega_cState &Ur);

   //! Returns Harten-Lax-van-Leer-Einfeldt (HLLE) flux in n-direction
   static FANS3D_ThermallyPerfect_KOmega_cState FluxHLLE_n(const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
                                                           const FANS3D_ThermallyPerfect_KOmega_pState &Wr,
                                                           const Vector3D &norm_dir);

   //! Returns Harten-Lax-van-Leer-Einfeldt (HLLE) flux in n-direction
   static FANS3D_ThermallyPerfect_KOmega_cState FluxHLLE_n(const FANS3D_ThermallyPerfect_KOmega_cState &Ul,
                                                           const FANS3D_ThermallyPerfect_KOmega_cState &Ur,
                                                           const Vector3D &norm_dir);

   //! Returns Roe flux x-direction flux
   static FANS3D_ThermallyPerfect_KOmega_cState FluxRoe_x(const FANS3D_ThermallyPerfect_KOmega_pState &Wl, 
                                                          const FANS3D_ThermallyPerfect_KOmega_pState &Wr);

   //! Returns Roe flux in n-direction
   static FANS3D_ThermallyPerfect_KOmega_cState FluxRoe_n(const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
                                                          const FANS3D_ThermallyPerfect_KOmega_pState &Wr,
                                                          const Vector3D &norm_dir);
   
   //! Returns negative waves speeds (eigenvalues) using Harten entropy fix
   FANS3D_ThermallyPerfect_KOmega_pState lambda_minus(const FANS3D_ThermallyPerfect_KOmega_pState &lambda_a,
                                                      const FANS3D_ThermallyPerfect_KOmega_pState &lambda_l,
                                                      const FANS3D_ThermallyPerfect_KOmega_pState &lambda_r);

   //! Returns positive waves speeds (eigenvalues) using Harten entropy fix
   FANS3D_ThermallyPerfect_KOmega_pState lambda_plus(const FANS3D_ThermallyPerfect_KOmega_pState &lambda_a,
                                                     const FANS3D_ThermallyPerfect_KOmega_pState &lambda_l,
                                                     const FANS3D_ThermallyPerfect_KOmega_pState &lambda_r);
//@}

/** @name Numerical Evaluation of Viscous Fluxes */
/*        -------------------------------------- */
//@{
   //! Returns viscous flux in n-direction
   static FANS3D_ThermallyPerfect_KOmega_cState FluxViscous_n(const FANS3D_ThermallyPerfect_KOmega_pState &Wl,
                                                              const FANS3D_ThermallyPerfect_KOmega_pState &Wr,
                                                              const FANS3D_ThermallyPerfect_KOmega_pState &W1,
                                                              const FANS3D_ThermallyPerfect_KOmega_pState &W2,
                                                              const FANS3D_ThermallyPerfect_KOmega_pState &dWdx1,
                                                              const FANS3D_ThermallyPerfect_KOmega_pState &dWdy1,
                                                              const FANS3D_ThermallyPerfect_KOmega_pState &dWdz1,
                                                              const FANS3D_ThermallyPerfect_KOmega_pState &dWdx2,
                                                              const FANS3D_ThermallyPerfect_KOmega_pState &dWdy2,
                                                              const FANS3D_ThermallyPerfect_KOmega_pState &dWdz2,
                                                              const Vector3D &norm, 
                                                              const Vector3D &ts, 
                                                              const double &deltad, 
                                                              const double &Volume, 
                                                              const double &Volume_Neigbor);
//@}

/** @name Boundary Conditions */
/*        ------------------- */
//@{
   //! Return reflected solution state after application of reflection BC
   static FANS3D_ThermallyPerfect_KOmega_pState Reflect(const FANS3D_ThermallyPerfect_KOmega_pState &W,
                                                        const Vector3D &norm_dir);

   //! Return wall solution state after application of moving wall BC
   static FANS3D_ThermallyPerfect_KOmega_pState MovingWall(const FANS3D_ThermallyPerfect_KOmega_pState &Win,
                                                           const FANS3D_ThermallyPerfect_KOmega_pState &Wout,
                                                           const Vector3D &norm_dir,
                                                           const Vector3D &wall_velocity,
                                                           const Vector3D &pressure_gradient,
                                                           const int &TEMPERATURE_BC_FLAG);

   //! Return wall solution state after application of no-slip BC
   static FANS3D_ThermallyPerfect_KOmega_pState NoSlip(const FANS3D_ThermallyPerfect_KOmega_pState &Win, 
                                                       const FANS3D_ThermallyPerfect_KOmega_pState &Wout, 
                                                       const Vector3D &norm_dir,  
                                                       const Vector3D &pressure_gradient,
                                                       const int &TEMPERATURE_BC_FLAG);
//@}
   
/** @name Turbulence Model Source Terms */
/*        ----------------------------- */
//@{
   static FANS3D_ThermallyPerfect_KOmega_cState Sturbulence(FANS3D_ThermallyPerfect_KOmega_pState &Wc,
                                                            const FANS3D_ThermallyPerfect_KOmega_pState &dWdx,
                                                            const FANS3D_ThermallyPerfect_KOmega_pState &dWdy,
                                                            const FANS3D_ThermallyPerfect_KOmega_pState &dWdz);
//@}

/** @name Operators */
/*        --------- */
//@{
   //! Index operator
   double &operator[](int index);
   //! Index operator
   const double &operator[](int index) const;
   
   //! Binary addition operator
   FANS3D_ThermallyPerfect_KOmega_pState operator +(const FANS3D_ThermallyPerfect_KOmega_pState &W) const;

   //! Binary subtraction operator
   FANS3D_ThermallyPerfect_KOmega_pState operator -(const FANS3D_ThermallyPerfect_KOmega_pState &W) const;

   //! Binary multiplication operator
   FANS3D_ThermallyPerfect_KOmega_pState operator *(const double &a) const;

   //! Binary multiplication operator
   friend FANS3D_ThermallyPerfect_KOmega_pState operator *(const double &a, 
                                                           const FANS3D_ThermallyPerfect_KOmega_pState &W);

   //! Binary multiplication operator
   double operator *(const FANS3D_ThermallyPerfect_KOmega_pState &W) const;

   //! Binary division operator
   FANS3D_ThermallyPerfect_KOmega_pState operator /(const double &a) const;

   //! Binary vector product operator
   FANS3D_ThermallyPerfect_KOmega_pState operator ^(const FANS3D_ThermallyPerfect_KOmega_pState &W) const;

   //! Assignment operator
   FANS3D_ThermallyPerfect_KOmega_pState& operator =(const FANS3D_ThermallyPerfect_KOmega_pState &W);

   //! Shortcut addition operator
   FANS3D_ThermallyPerfect_KOmega_pState& operator +=(const FANS3D_ThermallyPerfect_KOmega_pState &W);

   //! Shortcut subtraction operator
   FANS3D_ThermallyPerfect_KOmega_pState& operator -=(const FANS3D_ThermallyPerfect_KOmega_pState &W);

   //! Unary subtraction operators
   friend FANS3D_ThermallyPerfect_KOmega_pState operator -(const FANS3D_ThermallyPerfect_KOmega_pState &W);
  
   //! Equal relational operator
   friend int operator ==(const FANS3D_ThermallyPerfect_KOmega_pState &W1,
                          const FANS3D_ThermallyPerfect_KOmega_pState &W2);

   //! Not equal relational operator
   friend int operator !=(const FANS3D_ThermallyPerfect_KOmega_pState &W1,
                          const FANS3D_ThermallyPerfect_KOmega_pState &W2);

   //! Output stream operator
   friend ostream& operator << (ostream &out_file,
                                const FANS3D_ThermallyPerfect_KOmega_pState &W);

   //! Input stream operator
   friend istream& operator >> (istream &in_file,
                                FANS3D_ThermallyPerfect_KOmega_pState &W);
//@}
};

/*!
 * Class: FANS3D_ThermallyPerfect_KOmega_pState
 *
 * \brief Conserved state solution class for 3D Favre-average Navier-Stokes
 *        equations with k-omega turbulence model governing flows of thermally 
 *        perfect non-reactive and combusting mixtures.
 *
 * Member functions
 *  - rho      -- Return mixture density (kg/m^3)
 *  - rhov     -- Return mixture momentum (kg/(m^2-s))   
 *  - E        -- Return mixture total energy (J/kg)
 *  - rhok     -- Return turbulent kinetic energy (kg/m-s^2)
 *  - rhoomega -- Return specific dissipation rate for turbulent kinetic energy (kg/m^3-s)
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
class FANS3D_ThermallyPerfect_KOmega_cState : public NavierStokes3D_ThermallyPerfect_cState {
  public:
   double                      rhok;          //!< Turbulent kinetic energy (kg/m-s^2)
   double                      rhoomega;      //!< Specific dissipation rate for turbulent energy (kg/m^3-s)
   Turbulence_Model_k_omega    k_omega_model; //!< k-omega turbulence model data
   
/** @name Constructors and destructors */
/*        ---------------------------- */
//@{
   //! Default creation constructor (assign default values)
   FANS3D_ThermallyPerfect_KOmega_cState(void) : NavierStokes3D_ThermallyPerfect_cState(), rhok(ZERO) {
      rhoomega = rho*MILLION;
   }
   
   //! Constructor from base class (allows return of derived type)
   FANS3D_ThermallyPerfect_KOmega_cState(const NavierStokes3D_ThermallyPerfect_cState &U1) : 
     NavierStokes3D_ThermallyPerfect_cState(U1), rhok(ZERO), rhoomega() {
     rhoomega = rho*MILLION;
   }

   //! Constructor from base class
     FANS3D_ThermallyPerfect_KOmega_cState(const NavierStokes3D_ThermallyPerfect_cState &U1,
                                           const double &dk,
                                           const double &domega) : 
     NavierStokes3D_ThermallyPerfect_cState(U1), rhok(dk), rhoomega(domega) { }

   //! Assignment constructor
   FANS3D_ThermallyPerfect_KOmega_cState(const double &value) : 
     NavierStokes3D_ThermallyPerfect_cState(value), rhok(value), rhoomega(value) { }
   
   //! Assignment constructor
   FANS3D_ThermallyPerfect_KOmega_cState(const double &d, 
                                         const Vector3D &V,
                                         const double &En, 
                                         const double &dk,
                                         const double &domega) :
     NavierStokes3D_ThermallyPerfect_cState(d, V, En), rhok(dk), rhoomega(domega) { }

   //! Assignment constructor
   FANS3D_ThermallyPerfect_KOmega_cState(const double &d, 
                                         const Vector3D &V,
                                         const double &En, 
                                         const double &dk,
                                         const double &domega, 
                                         const double &frac) :
     NavierStokes3D_ThermallyPerfect_cState(d, V, En, frac), rhok(dk), rhoomega(domega) { }

   //! Assignment constructor
   FANS3D_ThermallyPerfect_KOmega_cState(const double &d, 
                                         const Vector3D &dV,
                                         const double &En,
                                         const double &dk,
                                         const double &domega,
                                         const Species *rhomfrac) :
     NavierStokes3D_ThermallyPerfect_cState(d, dV, En, rhomfrac), rhok(dk), rhoomega(domega) { }

   //! Assignment constructor
   FANS3D_ThermallyPerfect_KOmega_cState(const double &d, 
                                         const double &vx,
                                         const double &vy, 
                                         const double &vz,
                                         const double &En, 
                                         const double &dk,
                                         const double &domega) :
     NavierStokes3D_ThermallyPerfect_cState(d, vx, vy, vz, En), rhok(dk), rhoomega(domega) { }
   
   //! Assignment constructor
   FANS3D_ThermallyPerfect_KOmega_cState(const double &d, 
                                         const double &vx, 
                                         const double &vy, 
                                         const double &vz, 
                                         const double &En, 
                                         const double &dk,
                                         const double &domega,
                                         const double &rhomfrac) :
     NavierStokes3D_ThermallyPerfect_cState(d, vx, vy, vz, En, rhomfrac), rhok(dk), rhoomega(domega) { }

   //! Assignment constructor
   FANS3D_ThermallyPerfect_KOmega_cState(const double &d, 
                                         const double &vx,
                                         const double &vy, 
                                         const double &vz,
                                         const double &En, 
                                         const double &dk,
                                         const double &domega,
                                         const Species *rhomfrac) :
     NavierStokes3D_ThermallyPerfect_cState(d, vx, vy, vz, En, rhomfrac), rhok(dk), rhoomega(domega) { }
 
   //! Copy constructor (this is needed for the operator overload returns)
   FANS3D_ThermallyPerfect_KOmega_cState(const FANS3D_ThermallyPerfect_KOmega_cState &U) {
     rhospec_null(); rho = DENSITY_STDATM; set_initial_values(); Copy(U);
   }
   
   //! Default destructor
   ~FANS3D_ThermallyPerfect_KOmega_cState() {
      Deallocate();
   }
//@}
                           
/** @name Some useful operators */
/*        --------------------- */
//@{
   //! Copy solution state (cheaper than = operator)
   void Copy(const FANS3D_ThermallyPerfect_KOmega_cState &U) {
      Euler3D_ThermallyPerfect_cState::Copy(U); rhok = U.rhok; rhoomega = U.rhoomega;
   }

   //! Assigns a vacuum solution state
   void Vacuum(void) {
      Euler3D_ThermallyPerfect_cState::Vacuum(); rhok = ZERO; rhoomega = ZERO;
   }

   //! Check for physical validity of the solution vector
   bool Realizable_Solution_Check(void);
//@}

/** @name Thermodynamic and other state functions */
/*        --------------------------------------- */
//@{
   //! Mixture pressure
   double p(void) const;

   //! Mixture temperature
   double T(void) const;

   //! Mixture sound speed including turbulent kinetic energy
   double a_t(void) const;

   //! Turbulent kinetic energy
   double k(void) const;

   //! Specific dissipation rate for turbulent energy
   double omega(void) const;
//@}
   
/** @name Turbulent transport coefficients */
/*        -------------------------------- */
//@{
   //! Eddy (turbulent) viscosity
   double mu_t(void);

   //! Turbulent) thermal conductivity
   double kappa_t(void);     

   //! Species turbulent diffusion coefficient
   double Ds_t(const int &i);
   //! Species turbulent diffusion coefficient
   double Ds_t(const int &i,
               const double &mu_t_temp);

   //! Turbulent Prandtl number
   double Pr_t(void);

   //! Turbulent Schmidt number
   double Sc_t(void);

   //! Turbulent Lewis number
   double Le_t(void);
//@}

/** @name Primitive solution state */ 
/*        ------------------------ */
//@{
   //! Returns primitive solution state
   FANS3D_ThermallyPerfect_KOmega_pState W(void) ; 
   //! Returns primitive solution state
   FANS3D_ThermallyPerfect_KOmega_pState W(void) const;
   //! Returns primitive solution state
   FANS3D_ThermallyPerfect_KOmega_pState W(const FANS3D_ThermallyPerfect_KOmega_cState &U) const;
//@}

/** @name Operators */
/*        --------- */
//@{
   //! Index operator
   double &operator[](int index);
   //! Index operator
   const double &operator[](int index) const;
   
   //! Binary addition operator
   FANS3D_ThermallyPerfect_KOmega_cState operator +(const FANS3D_ThermallyPerfect_KOmega_cState &U) const;

   //! Binary subtraction operator
   FANS3D_ThermallyPerfect_KOmega_cState operator -(const FANS3D_ThermallyPerfect_KOmega_cState &U) const;

   //! Binary multiplication operator
   FANS3D_ThermallyPerfect_KOmega_cState operator *(const double &a) const;

   //! Binary multiplication operator
   friend FANS3D_ThermallyPerfect_KOmega_cState operator *(const double &a, 
                                                           const FANS3D_ThermallyPerfect_KOmega_cState &U);

   //! Binary multiplication operator
   double operator *(const FANS3D_ThermallyPerfect_KOmega_cState &U) const;

   //! Binary division operator
   FANS3D_ThermallyPerfect_KOmega_cState operator /(const double &a) const;

   //! Binary vector product operator
   FANS3D_ThermallyPerfect_KOmega_cState operator ^(const FANS3D_ThermallyPerfect_KOmega_cState &U) const;

   //! Assignment operator
   FANS3D_ThermallyPerfect_KOmega_cState& operator =(const FANS3D_ThermallyPerfect_KOmega_cState &U);

   //! Shortcut addition operator
   FANS3D_ThermallyPerfect_KOmega_cState& operator +=(const FANS3D_ThermallyPerfect_KOmega_cState &U);

   //! Shortcut subtraction operator
   FANS3D_ThermallyPerfect_KOmega_cState& operator -=(const FANS3D_ThermallyPerfect_KOmega_cState &U);

   //! Unary subtraction operators
   friend FANS3D_ThermallyPerfect_KOmega_cState operator -(const FANS3D_ThermallyPerfect_KOmega_cState &U);
  
   //! Equal relational operator
   friend int operator ==(const FANS3D_ThermallyPerfect_KOmega_cState &U1,
                          const FANS3D_ThermallyPerfect_KOmega_cState &U2);

   //! Not equal relational operator
   friend int operator !=(const FANS3D_ThermallyPerfect_KOmega_cState &U1,
                          const FANS3D_ThermallyPerfect_KOmega_cState &U2);

   //! Output stream operator
   friend ostream& operator << (ostream &out_file,
                                const FANS3D_ThermallyPerfect_KOmega_cState &U);

   //! Input stream operator
   friend istream& operator >> (istream &in_file,
                                FANS3D_ThermallyPerfect_KOmega_cState &U);
//@}
};


/***************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState member functions                              *
 ***************************************************************************************/

/***************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState -- Index operators.                           *
 ***************************************************************************************/
inline double& FANS3D_ThermallyPerfect_KOmega_pState::operator[](int index) {
   switch(index){  
   case 1:
      return rho;    
   case 2:
      return v.x;
   case 3:
      return v.y;
   case 4:
      return v.z;
   case 5:
      return p;
   case 6:
      return k;
   case 7:
      return omega;
   default :
      return spec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA +1)].c;
   };
}

inline const double& FANS3D_ThermallyPerfect_KOmega_pState::operator[](int index) const {
   switch(index){  
   case 1:
      return rho;    
   case 2:
      return v.x;
   case 3:
      return v.y;
   case 4:
      return v.z;
   case 5:
      return p;
   case 6:
      return k;
   case 7:
      return omega;
   default :
      return spec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA +1)].c;     
   };
}

/***************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState -- Binary arithmetic operators.   *
 ***************************************************************************/
//------------------ Addition ------------------------//
inline FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::
operator +(const FANS3D_ThermallyPerfect_KOmega_pState &W) const { 
   FANS3D_ThermallyPerfect_KOmega_pState Temp;
   Temp.Copy(*this);
   Temp += W;
   return Temp;
}

//------------------ Subtraction ------------------------//
inline FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::
operator -(const FANS3D_ThermallyPerfect_KOmega_pState &W) const {
   FANS3D_ThermallyPerfect_KOmega_pState Temp;
   Temp.Copy(*this);
   Temp -= W;
   return Temp;
}

//---------------- Scalar Multiplication ------------------//
inline FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::
operator *(const double &a) const {
   FANS3D_ThermallyPerfect_KOmega_pState Temp;
   Temp.Copy(*this);
   Temp.rho = rho*a;  
   Temp.v = v*a; 
   Temp.p = p*a;
   Temp.k = k*a; 
   Temp.omega = omega*a;
   for (int i = 0; i < ns; i++) {
      Temp.spec[i] = spec[i]*a;
   } /* endfor */
   return(Temp);
}

inline FANS3D_ThermallyPerfect_KOmega_pState operator *(const double &a, 
                                                        const FANS3D_ThermallyPerfect_KOmega_pState &W) {
   FANS3D_ThermallyPerfect_KOmega_pState Temp;
   //Temp.Copy(W);
   Temp.rho = W.rho*a;  
   Temp.v = W.v*a; 
   Temp.p = W.p*a;
   Temp.k = W.k*a; 
   Temp.omega = W.omega*a;
   for(int i = 0; i < W.ns; i++) {
      Temp.spec[i] = W.spec[i]*a;
   } /* endfor */
   return(Temp);
}

//--------------- Scalar Division ------------------------//
inline FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::
operator /(const double &a) const {
   FANS3D_ThermallyPerfect_KOmega_pState Temp(rho,v,p, k, omega);
   Temp.Copy(*this);
   Temp.rho = rho/a; 
   Temp.v = v/a; 
   Temp.p = p/a; 
   Temp.k = k/a; 
   Temp.omega = omega/a;
   for (int i = 0; i < ns; i++) {
      Temp.spec[i] = spec[i]/a; 
   } /* endfor */
   return(Temp);
}

//----------------- Inner Product ------------------------//
inline double FANS3D_ThermallyPerfect_KOmega_pState::
operator *(const FANS3D_ThermallyPerfect_KOmega_pState &W) const {
   double sum=0.0;
   for (int i=0; i < ns; i++) {
      sum += spec[i]*W.spec[i];
   } /* endfor */
   return (rho*W.rho + v*W.v + p*W.p + k*W.k + omega*W.omega + sum);
}

//----------- Solution state product operator ------------//
inline FANS3D_ThermallyPerfect_KOmega_pState FANS3D_ThermallyPerfect_KOmega_pState::
operator ^(const FANS3D_ThermallyPerfect_KOmega_pState &W) const {
   FANS3D_ThermallyPerfect_KOmega_pState Temp(rho,v,p, k, omega);
   Temp.Copy(*this);
   Temp.rho = rho*W.rho;
   Temp.v.x = v.x*W.v.x;
   Temp.v.y = v.y*W.v.y;
   Temp.v.z = v.z*W.v.z;
   Temp.p = p*W.p;
   Temp.k = k*W.k;
   Temp.omega = omega*W.omega;
   for (int i = 0; i < ns; i++){
      Temp.spec[i] = spec[i]*W.spec[i];
   } /* endfor */
   return(Temp);
}

//----------------- Assignment ----------------------------//
inline FANS3D_ThermallyPerfect_KOmega_pState& FANS3D_ThermallyPerfect_KOmega_pState::
operator =(const FANS3D_ThermallyPerfect_KOmega_pState &W) {
   //self assignment protection
   if (this != &W){   
      //copy assignment
      rho = W.rho;
      v = W.v; 
      p = W.p; 
      k = W.k;
      omega = W.omega;
      for (int i = 0; i < ns; i++) {
         spec[i] = W.spec[i];
      } /* endfor */ 
   } /* endif */
   return (*this);
}

/***************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState -- Shortcut arithmetic operators. *
 ***************************************************************************/
inline FANS3D_ThermallyPerfect_KOmega_pState& FANS3D_ThermallyPerfect_KOmega_pState::
operator +=(const FANS3D_ThermallyPerfect_KOmega_pState &W) {
   rho += W.rho;
   v += W.v; 
   p += W.p; 
   k += W.k;
   omega += W.omega;
   for (int i = 0; i < ns; i++) {
      spec[i].c += W.spec[i].c;
   } /* endfor */
   return (*this);
}

inline FANS3D_ThermallyPerfect_KOmega_pState& FANS3D_ThermallyPerfect_KOmega_pState::
operator -=(const FANS3D_ThermallyPerfect_KOmega_pState &W) {
   rho -= W.rho;
   v -= W.v;
   p -= W.p;
   k -= W.k;
   omega -= W.omega;
   for (int i = 0; i < ns; i++) {
      spec[i].c -= W.spec[i].c;
   } /* endfor */
   return (*this); 
}

/***************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState -- Unary arithmetic operators.    *
 ***************************************************************************/
inline FANS3D_ThermallyPerfect_KOmega_pState operator -(const FANS3D_ThermallyPerfect_KOmega_pState &W) {
   Species *spt= new Species[W.ns];
   for (int i = 0; i < W.ns; i++) {
      spt[i] = -W.spec[i]; 
   } /* endfor */ 
   FANS3D_ThermallyPerfect_KOmega_pState Temp(-W.rho, -W.v, -W.p, W.k, -W.omega, spt);
   delete[] spt;
   return(Temp);
}

/***************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState -- Relational operators.          *
 ***************************************************************************/
inline int operator ==(const FANS3D_ThermallyPerfect_KOmega_pState &W1, 
                       const FANS3D_ThermallyPerfect_KOmega_pState &W2) {
   bool Temp;
   for (int i = 0; i < W1.ns; i++) {
      if (W1.spec[i] == W2.spec[i]) {
         Temp = true;
      } else {
         Temp = false;
         break;
      } /* endif */
   } /* endfor */
   return (W1.rho == W2.rho && 
           W1.v == W2.v && 
           W1.p == W2.p && 
           W1.k == W2.k && 
           W1.omega == W2.omega && 
           Temp == true);
}

inline int operator !=(const FANS3D_ThermallyPerfect_KOmega_pState &W1, 
                       const FANS3D_ThermallyPerfect_KOmega_pState &W2) {
   bool Temp = true;
   for (int i = 0; i < W1.ns; i++) {
      if (W1.spec[i] != W2.spec[i]) {
         Temp = false;
         break;
      } /* endif */
   } /* endfor */
   return (W1.rho != W2.rho || 
           W1.v != W2.v || 
           W1.p != W2.p || 
           W1.k != W2.k || 
           W1.omega != W2.omega || 
           Temp != true);
}

/***************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_pState -- Input-output operators.        *
 ***************************************************************************/
inline ostream &operator << (ostream &out_file, 
                             const FANS3D_ThermallyPerfect_KOmega_pState &W) {
   out_file.precision(10);
   out_file.setf(ios::scientific);
   out_file << " " << W.rho  << " " << W.v.x << " " << W.v.y 
            << " " << W.v.z << " " << W.p << " "<<W.k<< " "<< W.omega;
   for (int i = 0; i < W.ns; i++) {
      out_file<<" "<<W.spec[i];
   } /* endfor */
   out_file.unsetf(ios::scientific);
   return (out_file);
}

inline istream &operator >> (istream &in_file, 
                             FANS3D_ThermallyPerfect_KOmega_pState &W) {
   in_file.setf(ios::skipws);
   in_file >> W.rho >> W.v.x >> W.v.y >> W.v.z >> W.p>> W.k>> W.omega;
   //W.set_initial_values();
   for (int i = 0; i < W.ns; i++) {
      in_file>>W.spec[i];
   } /* endfor */
   in_file.unsetf(ios::skipws);
   return (in_file);
}

/***************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState member functions                              *
 ***************************************************************************************/

/***************************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState -- Index operators.                           *
 ***************************************************************************************/
//index operators
inline double& FANS3D_ThermallyPerfect_KOmega_cState::operator[](int index) {
   switch(index){  
   case 1:
      return rho;    
   case 2:
      return rhov.x;
   case 3:
      return rhov.y;
   case 4:
      return rhov.z;
   case 5:
      return E;
   case 6:
      return rhok;
   case 7:
      return rhoomega;
   default :
      return rhospec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA +1)].c;     
   };
}

inline const double& FANS3D_ThermallyPerfect_KOmega_cState::operator[](int index) const {
   switch(index){  
   case 1:
      return rho;    
   case 2:
      return rhov.x;
   case 3:
      return rhov.y;
   case 4:
      return rhov.z;
   case 5:
      return E;
   case 6:
      return rhok;
   case 7:
      return rhoomega;
   default :
      return rhospec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_FANS3D_VAR_EXTRA +1)].c; 
   };
}

/***************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState -- Binary arithmetic operators.   *
 ***************************************************************************/
//------------------ Addition ------------------------//
inline FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_cState::
operator +(const FANS3D_ThermallyPerfect_KOmega_cState &U) const { 
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.Copy(*this);
   Temp += U;
   return Temp;
}

//------------------ Subtraction ------------------------//
inline FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_cState::
operator -(const FANS3D_ThermallyPerfect_KOmega_cState &U) const {
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.Copy(*this);
   Temp -= U;
   return Temp;
}

//---------------- Scalar Multiplication ------------------//
inline FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_cState::
operator *(const double &a) const {
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.Copy(*this);
   Temp.rho = rho*a;  
   Temp.rhov = rhov*a; 
   Temp.E = E*a;
   Temp.rhok = rhok*a; 
   Temp.rhoomega = rhoomega*a;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i] = rhospec[i]*a;
   } /* endfor */
   return(Temp);
}

inline FANS3D_ThermallyPerfect_KOmega_cState operator *(const double &a, 
                                                        const FANS3D_ThermallyPerfect_KOmega_cState &U) {
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.rho = U.rho*a;  Temp.rhov = U.rhov*a; Temp.E = U.E*a;
   Temp.rhok = U.rhok*a; Temp.rhoomega = U.rhoomega*a;
   for (int i = 0; i < U.ns; i++) {
      Temp.rhospec[i] = U.rhospec[i]*a;
   } /* endfor */
   return(Temp);
}

//--------------- Scalar Division ------------------------//
inline FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_cState::
operator /(const double &a) const {
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.Copy(*this);
   Temp.rho = rho/a; Temp.rhov = rhov/a; Temp.E = E/a;
   Temp.rhok = rhok/a; Temp.rhoomega = rhoomega/a;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i] = rhospec[i]/a; 
   } /* endfor */
   return(Temp);
}

//----------------- Inner Product ------------------------//
inline double FANS3D_ThermallyPerfect_KOmega_cState::
operator *(const FANS3D_ThermallyPerfect_KOmega_cState &U) const {
   double sum=0.0;
   for (int i = 0; i < ns; i++) {
      sum += rhospec[i]*U.rhospec[i];
   } /* endfor */
   return (rho*U.rho + rhov*U.rhov + E*U.E + rhok*U.rhok + rhoomega*U.rhoomega + sum);
}

//----------- Solution state product operator ------------//
inline FANS3D_ThermallyPerfect_KOmega_cState FANS3D_ThermallyPerfect_KOmega_cState::
operator ^(const FANS3D_ThermallyPerfect_KOmega_cState &U) const {
   FANS3D_ThermallyPerfect_KOmega_cState Temp;
   Temp.Copy(*this);
   Temp.rho = rho*U.rho;
   Temp.rhov.x = rhov.x*U.rhov.x;
   Temp.rhov.y = rhov.y*U.rhov.y;
   Temp.rhov.z = rhov.z*U.rhov.z;
   Temp.E = E*U.E;
   Temp.rhok = rhok*U.rhok;
   Temp.rhoomega = rhoomega*U.rhoomega;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i] = rhospec[i]*U.rhospec[i];
   } /* endfor */
   return(Temp);
}

//----------------- Assignment ----------------------------//
inline FANS3D_ThermallyPerfect_KOmega_cState& FANS3D_ThermallyPerfect_KOmega_cState::
operator =(const FANS3D_ThermallyPerfect_KOmega_cState &U) {
   //self assignment protection
   if (this != &U) {   
      //copy assignment
      rho = U.rho;
      rhov = U.rhov; 
      E = U.E; 
      rhok = U.rhok;
      rhoomega = U.rhoomega;
      for (int i = 0; i < ns; i++) {
         rhospec[i] = U.rhospec[i];
      } /* endfor */
   } /* endif */
   return (*this);
}

/***************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState -- Shortcut arithmetic operators. *
 ***************************************************************************/
inline FANS3D_ThermallyPerfect_KOmega_cState& FANS3D_ThermallyPerfect_KOmega_cState::
operator +=(const FANS3D_ThermallyPerfect_KOmega_cState &U) {
   rho += U.rho;
   rhov += U.rhov; 
   E += U.E;
   rhok += U.rhok;
   rhoomega += U.rhoomega;
   for (int i = 0; i < ns; i++) {
      rhospec[i].c += U.rhospec[i].c;
   } /* endfor */
   return (*this);
}

inline FANS3D_ThermallyPerfect_KOmega_cState& FANS3D_ThermallyPerfect_KOmega_cState::
operator -=(const FANS3D_ThermallyPerfect_KOmega_cState &U) {
   rho -= U.rho;
   rhov -= U.rhov;
   E -= U.E;
   rhok -= U.rhok;
   rhoomega -= U.rhoomega;
   for (int i = 0; i < ns; i++) {
     rhospec[i].c -= U.rhospec[i].c;
   } /* endfor */
   return (*this); 
}

/***************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState -- Unary arithmetic operators.    *
 ***************************************************************************/
inline FANS3D_ThermallyPerfect_KOmega_cState operator -(const FANS3D_ThermallyPerfect_KOmega_cState &U) {
   Species *spt= new Species[U.ns];
   for (int i=0; i < U.ns; i++) {
      spt[i] = -U.rhospec[i]; 
   } /* endfor */
   FANS3D_ThermallyPerfect_KOmega_cState Temp(-U.rho,-U.rhov,-U.E, -U.rhok, -U.rhoomega, spt);
   delete[] spt;
   return(Temp);
}

/***************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState -- Relational operators.          *
 ***************************************************************************/
inline int operator ==(const FANS3D_ThermallyPerfect_KOmega_cState &U1, 
                       const FANS3D_ThermallyPerfect_KOmega_cState &U2) {
   bool Temp;
   for (int i = 0; i < U1.ns; i++) {
      if (U1.rhospec[i] == U2.rhospec[i] ) {
         Temp = true;
      } else {
         Temp = false;
         break;
      } /* endif */
   } /* endfor */
   return (U1.rho == U2.rho && 
           U1.rhov == U2.rhov && 
           U1.E == U2.E && 
           U1.rhok == U2.rhok && 
           U1.rhoomega == U2.rhoomega &&
           Temp == true);
}

inline int operator !=(const FANS3D_ThermallyPerfect_KOmega_cState &U1, 
                       const FANS3D_ThermallyPerfect_KOmega_cState &U2) {
   bool Temp = true;
   for (int i = 0; i < U1.ns; i++) {
      if( U1.rhospec[i] != U2.rhospec[i] ){
         Temp = false;
         break;
      } /* endif */
   } /* endfor */
   return (U1.rho != U2.rho || 
           U1.rhov != U2.rhov || 
           U1.E != U2.E ||
           U1.rhok != U2.rhok || 
           U1.rhoomega !=U2.rhoomega || 
           Temp != true);
}

/***************************************************************************
 * FANS3D_ThermallyPerfect_KOmega_cState -- Input-output operators.        *
 ***************************************************************************/
inline ostream &operator << (ostream &out_file, 
                             const FANS3D_ThermallyPerfect_KOmega_cState &U) {
  //out_file.precision(20);
  out_file.setf(ios::scientific);
  out_file << " " << U.rho  << " " << U.rhov.x << " " << U.rhov.y << " " 
           << U.rhov.z << " " << U.E << " " << U.rhok << " " << U.rhoomega;
  for (int i = 0; i < U.ns; i++) {
     out_file<<" "<<U.rhospec[i];
  } /* endfor */
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, 
                             FANS3D_ThermallyPerfect_KOmega_cState &U) {
   in_file.setf(ios::skipws);
   in_file >> U.rho >> U.rhov.x >> U.rhov.y >> U.rhov.z >> U.E >> U.rhok >> U.rhoomega;
   for (int i = 0; i < U.ns; i++) {
     in_file>>U.rhospec[i]; 
   } /* endfor */
   in_file.unsetf(ios::skipws);
   return (in_file);
}

#endif // _FANS3D_THERMALLYPERFECT_STATE_INCLUDED 
