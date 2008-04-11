/*! \file LES3DThickenedFlameState.h
 * 	\brief	Header file defining the Favre-filtered Navier-Stokes solution state 
 *              classes associated with solution of premixed compressible turbulent 
 *              combusting flows of a thermally perfect using a large eddy simulation
 *              (LES) technique in conjunction with a thickened flame
 *              subfilter scale model.
 * .
 */

#ifndef _LES3DTF_STATE_INCLUDED 
#define _LES3DTF_STATE_INCLUDED

/* Include header file for base solution classes from which the classes are derived. */

#ifndef _NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED
#include "../../NavierStokes/NavierStokes3DThermallyPerfectState.h"
#endif  //NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED

/* Include turbulence modelling header file. */

#ifndef _TURBULENCE_MODELLING_INCLUDED
#include "../../TurbulenceModelling/TurbulenceModelling.h"
#endif  //TURBULENCE_MODELLING_INCLUDED

/* Include power law header file. */

#ifndef _POWER_LAW_INCLUDED
#include "PowerLaw.h"
#endif // _POWER_LAW_INCLUDED

/* Define the classes. */

class LES3DTF_pState;
class LES3DTF_cState;

#define NUM_LES3DTF_VAR_EXTRA 1  // k


/* Define some constants. */
const double   GRAVITY_Z(-9.81);      // acceleration due to gravity [m/s^2]

const double   Cv_CONSTANT(0.086);    // constant used in k-equation (eddy viscosity)
const double   Ceps_CONSTANT(0.845);  // constant used in k-equation (dissipation term)

const double   Cs_CONSTANT(0.18);     // Smagorinsky constant
const double   CI_CONSTANT(0.005);    // constant used in Yoshizawa's model for the SFS k

const double   TWO_THIRDS(2.0/3.0); 
const double   FIVE_THIRDS(5.0/3.0);

/*!
 * Class: LES3DTF_pState
 *
 * \brief Primitive state solution class for 3D Favre-filtered Navier-Stokes
 *        equations associated with solution of premixed compressible turbulent 
 *        combusting flows of a thermally perfect using a large eddy simulation
 *        (LES) technique in conjunction with a thickened flame
 *        subfilter scale model.
 *
 * Member functions
 *  - rho              -- Return density (kg/m^3)
 *  - v                -- Return flow velocity (m/s)
 *  - p                -- Return pressure (Pa, N/m^2)
 *  - p_t              -- Return turbulence modified pressure (Pa, N/m^2)
 *  - k                -- Return subfilter scale turbulent kinetic energy (m^2/s^2)
 *  - spec             -- Return array of species mass fraction data
 *  - Mass             -- Return mixture molecular mass (kg/mol)
 *  - Rtot             -- Return mixture gas constant (J/(kg*K))
 *  - HeatofFormation  -- Return heat of formation for the mixture
 *  - Cp               -- Return specific heat at constant pressure for mixture (J/(kg*K))
 *  - Cv               -- Return specific heat at constant volume for mixture (J/(kg*K))
 *  - g                -- Return specific heat ratio for mixture
 *  - e                -- Return mixture absolute internal energy (J/kg)
 *  - es               -- Return mixture sensible internal energy (J/kg)
 *  - h                -- Return mixture absolute specific enthalpy (J/kg)
 *  - hs               -- Return mixture sensible specific enthalpy (J/kg)
 *  - E                -- Return total mixture energy (J/kg)
 *  - H                -- Return total mixture enthalpy (J/kg)
 *  - rhov             -- Return momentum of mixture (kg/(m^2*s))
 *  - T                -- Return mixture temperature (K)
 *  - a                -- Return sound speed of mixture (m/s)
 *  - a_t              -- Return sound speed of mixture including turbulent kinetic energy (m/s)
 *  - M                -- Return Mach number for mixture
 *  - Gibbs            -- Return species Gibbs free energy
 *  - mu               -- Return mixture dynamic viscosity
 *  - kappa            -- Return mixture coefficient of thermal conductivity
 *  - Ds               -- Return species mass diffusion coefficient
 *  - Pr               -- Return mixture Prandtl number
 *  - Sc               -- Return species Schmidt number
 *  - Le               -- Return species Lewis number
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
 *  - Schemistry       -- Return source terms associated with finite-rate chemistry
 *  - Sturbchem        -- Return source terms associated with modelled turbulence-chemistry interaction
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
class LES3DTF_pState : public NavierStokes3D_ThermallyPerfect_pState {
  public:
   static double           _laminar_flame_speed;  //!< Propagation speed of laminar premixed flame
   static double       _laminar_flame_thickness;  //!< Thickness of laminar premixed flame
   static double                       _TFactor;  //!< Maximum thickening factor
   static double                  _filter_width;  //!< Constant filter width

   double                     k; //!< Subfilter scale turbulent kinetic energy (m^2/s^2)

   PowerLaw               flame; //!< Object containing the thickening and wrinkling factors. 
 
   static double           Mref; //!< Reference Mach number for low-Mach-number precondtioning (normally set to incoming freestream Mach)  


/** @name Constructors and destructors */
/*        ---------------------------- */
//@{
   //! Default creation constructor (assign default values)
   LES3DTF_pState(void): 
     NavierStokes3D_ThermallyPerfect_pState(), k(ZERO) { }

   //! Constructor from base class (allows return of derived type)
   LES3DTF_pState(const NavierStokes3D_ThermallyPerfect_pState &W1) : 
     NavierStokes3D_ThermallyPerfect_pState(W1), k(ZERO) { }

   //! Constructor from base class
   LES3DTF_pState(const NavierStokes3D_ThermallyPerfect_pState &W1,
                  const double &ke) : 
     NavierStokes3D_ThermallyPerfect_pState(W1), k(ke) { }

   //! Assignment constructor
   LES3DTF_pState(const double &value) :
     NavierStokes3D_ThermallyPerfect_pState(value), k(value) { }
   
   //! Assignment constructor
   LES3DTF_pState(const double &d, 
                  const Vector3D &V,
                  const double &pre, 
                  const double &ke) :
     NavierStokes3D_ThermallyPerfect_pState(d, V, pre), k(ke) { }

   //! Assignment constructor
   LES3DTF_pState(const double &d, 
                  const Vector3D &V,
                  const double &pre, 
                  const double &ke, 
                  const double &frac) :
     NavierStokes3D_ThermallyPerfect_pState(d, V, pre, frac), k(ke) { }
   
   //! Assignment constructor
   LES3DTF_pState(const double &d, 
                  const Vector3D &V,
                  const double &pre, 
                  const double &ke, 
                  const Species *mfrac) :
     NavierStokes3D_ThermallyPerfect_pState(d, V, pre, mfrac), k(ke) { }

   //! Assignment constructor
   LES3DTF_pState(const double &d, 
                  const Vector3D &V,
                  const double &pre, 
                  const double *mfrac) :
     NavierStokes3D_ThermallyPerfect_pState(d, V, pre, mfrac), k(ZERO) { set_initial_values(mfrac); }

   //! Assignment constructor
   LES3DTF_pState(const double &d, 
                  const Vector3D &V,
                  const double &pre, 
                  const double &ke, 
                  const double *mfrac) :
     NavierStokes3D_ThermallyPerfect_pState(d, V, pre, mfrac), k(ke) { set_initial_values(mfrac); }

   //! Assignment constructor
   LES3DTF_pState(const double &d, 
                  const double &vx,
                  const double &vy, 
                  const double &vz,
                  const double &ke, 
                  const double &pre) :
     NavierStokes3D_ThermallyPerfect_pState(d, vx, vy, vz, pre), k(ke) { }

   //! Assignment constructor
   LES3DTF_pState(const double &d, 
                  const double &vx,
                  const double &vy, 
                  const double &vz,
                  const double &pre, 
                  const double &ke, 
                  const double &frac) :
     NavierStokes3D_ThermallyPerfect_pState(d, vx, vy, vz, pre, frac), k(ke) { }
   
   //! Assignment constructor
   LES3DTF_pState(const double &d, 
                  const double &vx,
                  const double &vy, 
                  const double &vz,
                  const double &pre, 
                  const double &ke, 
                  const Species *mfrac) :
     NavierStokes3D_ThermallyPerfect_pState(d, vx, vy, vz, pre, mfrac), k(ke) { }

   //! Copy constructor (this is needed for the operator overload returns)
   LES3DTF_pState(const  LES3DTF_pState &W) { 
     species_null();
     set_initial_values(); 
     Copy(W); 
   }

   //! Default destructor
   ~LES3DTF_pState(void) {
      Deallocate();
   }
//@}

/** @name Some useful operators */
/*        --------------------- */
//@{
   //! Assign the species data (needs to be called only once as it is static)
   void set_species_data(const int &n,
                         const string *S,
                         const char *PATH,
                         const int &debug, 
                         const double &Mr, 
                         const double* Sc,
                         const int &trans_data) {
      Euler3D_ThermallyPerfect_pState::set_species_data(n, S, PATH, debug, Mr, Sc, trans_data);
      num_vars = NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA + ns;
   }

   //! Assign modelling data (needs to be called only once as it is static)
   void set_modelling_data(const double &laminar_flame_speed,
			   const double &laminar_flame_thickness,
			   const double &thickening_factor,
			   const double &filter_width) {
     _laminar_flame_speed = laminar_flame_speed;
     _laminar_flame_thickness = laminar_flame_thickness;
     _TFactor = thickening_factor;
     _filter_width = filter_width;    
   }

   //! Copy solution state (cheaper than = operator)
   void Copy(const LES3DTF_pState &W);

   //! Assigns a vacuum solution state
   void Vacuum(void) { 
     rho = ZERO, v.zero(); p = ZERO; k = ZERO;
     for(int i=0; i<ns; ++i)  spec[i].Vacuum();
   }

   //! Check for physical validity of scalars
   void Realizable_Scalar_Check(void) {
     if ( k < MICRO ) { k = ZERO; }
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

   //! Mixture sound speed
   double a(void) const;

   //! Modified pressure
   double p_t(void) const;

   //! Mixture sound speed including turbulent kinetic energy
   double a_t(void);

   //! Mixture sound speed including turbulent kinetic energy 
   double a_t(void) const;

   //! Species mass fractions
   void premixed_mfrac(void);

//@}
  
/** @name Turbulent transport coefficients */
/*        -------------------------------- */
//@{
   //! Subfilter scale turbulent viscosity
   double mu_t(const LES3DTF_pState &dWdx, 
               const LES3DTF_pState &dWdy, 
               const LES3DTF_pState &dWdz,
               const int Flow_Type, 
               const double &Volume);

   //! Turbulent thermal conductivity
   double kappa_t(const LES3DTF_pState &dWdx, 
                  const LES3DTF_pState &dWdy, 
                  const LES3DTF_pState &dWdz,
                  const int Flow_Type, 
                  const double &Volume);

   //! Turbulent thermal conductivity
   double kappa_t(const double &mu_t_temp);

   //! Species turbulent diffusion coefficient
   double Ds_t(const LES3DTF_pState &dWdx, 
               const LES3DTF_pState &dWdy, 
               const LES3DTF_pState &dWdz,
               const int Flow_Type, 
               const double &Volume);

   //! Species turbulent diffusion coefficient
   double Ds_t(const double &mu_t_temp);

   //! Turbulent Prandtl number
   double Pr_t(void);

   //! Turbulent Schmidt number
   double Sc_t(void);      

   //! Turbulent Lewis number
   double Le_t(void); 

   //! LES filter width
   double filter_width(const double &Volume) const;
   double filter_width(void) const;
//@}

/** @name Subfilter scale turbulent stress tensor and heat flux vector */
/*        ------------------------------------------------------------ */
//@{
   //! Returns subfilter scale turbulent stress tensor 
   Tensor3D tau_t(const LES3DTF_pState &dWdx,
                  const LES3DTF_pState &dWdy,
                  const LES3DTF_pState &dWdz,
                  const int Flow_Type,
                  const double &Volume);

   //! Returns subfilter scale turbulent stress tensor 
   Tensor3D tau_t(const double &mu_t_temp, 
                  const LES3DTF_pState &dWdx, 
                  const LES3DTF_pState &dWdy,
                  const LES3DTF_pState &dWdz,
                  const int Flow_Type, 
                  const double &Volume);

   //! Returns components of subfilter scale turbulent stress tensor in x-direction
   Tensor3D tau_t_x(const double &mu_t_temp, 
                    const LES3DTF_pState &dWdx, 
                    const LES3DTF_pState &dWdy,
                    const LES3DTF_pState &dWdz,
                    const int Flow_Type, 
                    const double &Volume);

   //! Returns components of subfilter scale turbulent stress tensor in y-direction
   Tensor3D tau_t_y(const double &mu_t_temp, 
                    const LES3DTF_pState &dWdx, 
                    const LES3DTF_pState &dWdy,
                    const LES3DTF_pState &dWdz,
                    const int Flow_Type, 
                    const double &Volume);

   //! Returns components of subfilter scale turbulent stress tensor in z-direction
   Tensor3D tau_t_z(const double &mu_t_temp, 
                    const LES3DTF_pState &dWdx, 
                    const LES3DTF_pState &dWdy,
                    const LES3DTF_pState &dWdz,
                    const int Flow_Type, 
                    const double &Volume);

   //! Returns subfilter scale turbulent heat flux vector 
   Vector3D q_t(const double &kappa_t_temp,
                const LES3DTF_pState &dWdx, 
                const LES3DTF_pState &dWdy,
                const LES3DTF_pState &dWdz,
                const int Flow_Type, 
                const double &Volume);

   //! Returns component of subfilter scale turbulent heat flux vector in x-direction
   Vector3D q_t_x(const double &kappa_t_temp,
                  const LES3DTF_pState &dWdx, 
                  const LES3DTF_pState &dWdy,
                  const LES3DTF_pState &dWdz,
                  const int Flow_Type, 
                  const double &Volume);

   //! Returns component of subfilter scale turbulent heat flux vector in y-direction
   Vector3D q_t_y(const double &kappa_t_temp,
                  const LES3DTF_pState &dWdx, 
                  const LES3DTF_pState &dWdy,
                  const LES3DTF_pState &dWdz,
                  const int Flow_Type, 
                  const double &Volume);

   //! Returns component of subfilter scale turbulent heat flux vector in z-direction
   Vector3D q_t_z(const double &kappa_t_temp,
                  const LES3DTF_pState &dWdx, 
                  const LES3DTF_pState &dWdy,
                  const LES3DTF_pState &dWdz,
                  const int Flow_Type, 
                  const double &Volume);

   //! Returns subfilter thermal diffusion flux vector (due to species diffusion processes) 
   Vector3D thermal_diffusion_t(const double &mu_t_temp,
				const LES3DTF_pState &dWdx, 
				const LES3DTF_pState &dWdy,
				const LES3DTF_pState &dWdz);

   //! Returns components of subfilter thermal diffusion flux vector in x-direction (due to species diffusion processes) 
   Vector3D thermal_diffusion_t_x(const double &mu_t_temp,
				  const LES3DTF_pState &dWdx, 
				  const LES3DTF_pState &dWdy,
				  const LES3DTF_pState &dWdz);

   //! Returns components of subfilter thermal diffusion flux vector in y-direction (due to species diffusion processes) 
   Vector3D thermal_diffusion_t_y(const double &mu_t_temp,
				  const LES3DTF_pState &dWdx, 
				  const LES3DTF_pState &dWdy,
				  const LES3DTF_pState &dWdz);

   //! Returns components of subfilter thermal diffusion flux vector in z-direction (due to species diffusion processes) 
   Vector3D thermal_diffusion_t_z(const double &mu_t_temp,
				  const LES3DTF_pState &dWdx, 
				  const LES3DTF_pState &dWdy,
				  const LES3DTF_pState &dWdz);

//@}

/** @name Conserved solution state */ 
/*        ------------------------ */
//@{
   //! Returns conserved solution state
   LES3DTF_cState U(void);
   //! Returns conserved solution state
   LES3DTF_cState U(void) const;
   //! Returns conserved solution state
   LES3DTF_cState U(const LES3DTF_pState &W); 
//@}
   
/** @name Inviscid Flux Vectors */
/*        --------------------- */
//@{
   //! x-direction inviscid solution flux
   LES3DTF_cState F(void) const ;

   //! x-direction inviscid solution flux
   LES3DTF_cState Fx(void) const;

   //! y-direction inviscid solution flux
   LES3DTF_cState Fy(void) const;

   //! z-direction inviscid solution flux
   LES3DTF_cState Fz(void) const;
//@}

/** @name Viscous flux vectors */
/*        --------------------- */
//@{
   //! x-direction viscous solution flux
   LES3DTF_cState Fv(const LES3DTF_pState &dWdx,
                     const LES3DTF_pState &dWdy,
                     const LES3DTF_pState &dWdz,
                     const int Flow_Type,
                     const double &Volume);

   //! x-direction viscous solution flux
   LES3DTF_cState Fvx(const LES3DTF_pState &dWdx,
                      const LES3DTF_pState &dWdy,
                      const LES3DTF_pState &dWdz,
                      const int Flow_Type,
                      const double &Volume);

   //! y-direction viscous solution flux
   LES3DTF_cState Fvy(const LES3DTF_pState &dWdx,
                      const LES3DTF_pState &dWdy,
                      const LES3DTF_pState &dWdz,
                      const int Flow_Type,
                      const double &Volume);

   //! z-direction viscous solution flux
   LES3DTF_cState Fvz(const LES3DTF_pState &dWdx,
                      const LES3DTF_pState &dWdy,
                      const LES3DTF_pState &dWdz,
                      const int Flow_Type,
                      const double &Volume);
//@}

/** @name Eigenvalue(s) and eigenvector(s) (x-direction) */
/*        ---------------------------------------------- */
//@{
   //! x-direction eigenvalues
   LES3DTF_pState lambda(void) const;
   
   //! x-direction eigenvalues
   LES3DTF_pState lambda_x(void) const;

   //! x-direction conservative eigenvectors
   LES3DTF_cState rc(const int &index) const;

   //! x-direction conservative eigenvectors
   LES3DTF_cState rc_x(const int &index) const;

   //! x-direction primitive eignenvectors
   LES3DTF_pState lp(const int &index) const;

   //! x-direction primitive eignenvectors
   LES3DTF_pState lp_x(const int &index) const;  
//@}

/** @name Eigenvector(s) and transform matrix for preconditioner */
/*        ------------------------------------------------------ */
//@{
   //! Return the square value of Mach number
   double Mr2(const double &deltax, 
              const double &lengthx, 
              const double &dTime);

   //! Preconditioned velocity
   double u_plus_aprecon(const double &u, 
                         const double &deltax, 
                         const double &lengthx, 
                         const double &dTime);

   //! Preconditioned velocity and sound speed
   void u_a_precon(const double &UR,
                   double &uprimed, 
                   double &cprimed) const;

   //! x-direction conservative eigenvectors
   LES3DTF_cState rc_x_precon(const int &index, 
                               const double &MR2) const; 

   //! x-direction primitive eigenvectors
   LES3DTF_pState lp_x_precon(const int &index, 
                               const double &MR2) const; 

   //! Preconditioner matrix
   void Low_Mach_Number_Preconditioner(DenseMatrix &P, 
                                       const double &deltax, 
                                       const double &lengthx, 
                                       const double &dTime); 

   //! Preconditioner matrix inverse
   void Low_Mach_Number_Preconditioner_Inverse(DenseMatrix &Pinv, 
                                               const double &deltax, 
                                               const double &lengthx, 
                                               const double &dTime); 
//@}

/** @name Numerical Flux Functions */
/*        ------------------------ */
//@{
   //! Returns Roe-averaged primitive solution state
   static LES3DTF_pState RoeAverage(const LES3DTF_pState &Wl,
                                    const LES3DTF_pState &Wr);

   //! Returns Harten-Lax-van-Leer-Einfeldt (HLLE) x-direction flux
   static LES3DTF_cState FluxHLLE_x(const LES3DTF_pState &Wl,
                                    const LES3DTF_pState &Wr);

   //! Returns Harten-Lax-van-Leer-Einfeldt (HLLE) flux in n-direction
   static LES3DTF_cState FluxHLLE_n(const LES3DTF_pState &Wl,
                                    const LES3DTF_pState &Wr,
                                    const Vector3D &norm_dir);

   //! Returns Roe flux x-direction flux
   static LES3DTF_cState FluxRoe_x(const LES3DTF_pState &Wl, 
                                   const LES3DTF_pState &Wr);

   //! Returns Roe flux in n-direction
   static LES3DTF_cState FluxRoe_n(const LES3DTF_pState &Wl,
                                   const LES3DTF_pState &Wr,
                                   const Vector3D &norm_dir);

   //! Returns AUSMplus_up flux x-direction flux
   static LES3DTF_cState FluxAUSMplus_up_x(const LES3DTF_pState &Wl,
                                           const LES3DTF_pState &Wr);

   //! Returns AUSMplus_up flux in n-direction
   static LES3DTF_cState FluxAUSMplus_up_n(const LES3DTF_pState &Wl,
                                           const LES3DTF_pState &Wr,
                                           const Vector3D &norm_dir);

   //! Returns negative waves speeds (eigenvalues) using Harten entropy fix
   LES3DTF_pState lambda_minus(const LES3DTF_pState &lambda_a,
                               const LES3DTF_pState &lambda_l,
                               const LES3DTF_pState &lambda_r);

   //! Returns positive waves speeds (eigenvalues) using Harten entropy fix
   LES3DTF_pState lambda_plus(const LES3DTF_pState &lambda_a,
                              const LES3DTF_pState &lambda_l,
                              const LES3DTF_pState &lambda_r);

   //! HLLE wavespeeds in n-direction given 2 primitive states and a direction
   static Vector2D HLLE_wavespeeds(const LES3DTF_pState &Wl,
                                   const LES3DTF_pState &Wr,
                                   const Vector3D &norm_dir);
  
   //! Returns rotated primitive state aligned with local x-axis in the norm_dir
   LES3DTF_pState Rotate(const Vector3D &norm_dir) const;				  
  
   //! Returns un-rotated primitive state aligned with x-axis of global problem
   LES3DTF_pState RotateBack(const Vector3D &norm_dir) const;
//@}

/** @name Numerical Evaluation of Viscous Fluxes */
/*        -------------------------------------- */
//@{
   //! Returns viscous flux in n-direction
   static LES3DTF_cState FluxViscous_n(const LES3DTF_pState &Wl,
                                       const LES3DTF_pState &Wr,
                                       const LES3DTF_pState &W1,
                                       const LES3DTF_pState &W2,
                                       const LES3DTF_pState &dWdx1,
                                       const LES3DTF_pState &dWdy1,
                                       const LES3DTF_pState &dWdz1,
                                       const LES3DTF_pState &dWdx2,
                                       const LES3DTF_pState &dWdy2,
                                       const LES3DTF_pState &dWdz2,
                                       const Vector3D &norm, 
                                       const Vector3D &ts, 
                                       const double &deltad, 
                                       const double &Volume, 
                                       const double &Volume_Neighbour, 
                                       const int Flow_Type);
//@}

/** @name Boundary Conditions */
/*        ------------------- */
//@{
   //! Return reflected solution state after application of reflection BC
   static LES3DTF_pState Reflect(const LES3DTF_pState &W,
                                 const Vector3D &norm_dir);

   //! Return wall solution state after application of moving wall BC
   static LES3DTF_pState MovingWall(const LES3DTF_pState &Win,
                                    const LES3DTF_pState &Wout,
                                    const Vector3D &norm_dir,				 
                                    const Vector3D &wall_velocity,
                                    const Vector3D &pressure_gradient,
                                    const int &TEMPERATURE_BC_FLAG);

   //! Return wall solution state after application of no-slip BC
   static LES3DTF_pState NoSlip(const LES3DTF_pState &Win, 
                                const LES3DTF_pState &Wout, 
                                const Vector3D &norm_dir,  
                                const Vector3D &pressure_gradient,
                                const int &TEMPERATURE_BC_FLAG);
//@}

/** @name Turbulence Model Source Terms */
/*        ----------------------------- */
//@{
   //! Enstrophy
   double Enstrophy(const LES3DTF_pState &dWdx, 
                    const LES3DTF_pState &dWdy, 
                    const LES3DTF_pState &dWdz) const;

   //! Q-criterion to detect coherent structures of turbulence
   double Q_criterion(const LES3DTF_pState &dWdx, 
		      const LES3DTF_pState &dWdy, 
		      const LES3DTF_pState &dWdz) const;

   //! Absolute value of strain rate
   double abs_strain_rate(const LES3DTF_pState &dWdx, 
                          const LES3DTF_pState &dWdy, 
                          const LES3DTF_pState &dWdz) const;

   //! Subfilter scale kinetic energy
   double SFS_Kinetic_Energy(const LES3DTF_pState &dWdx,
                             const LES3DTF_pState &dWdy,
                             const LES3DTF_pState &dWdz,
                             const int Flow_Type,
                             const double &Volume);

   //! Source terms associated with finite-rate chemistry 
   LES3DTF_cState Sw(const int &REACT_SET_FLAG);
   double dSwdU_max_diagonal(void) const;

   //! Source terms associated with gravitational forces
   LES3DTF_cState Sg(void) const;

   //! Source term for k-equation
   double K_equ_sources(const LES3DTF_pState &dWdx,
                        const LES3DTF_pState &dWdy,
                        const LES3DTF_pState &dWdz,
                        const int Flow_Type,
                        const double &Volume);

   //! Source terms for k-eq model
   static LES3DTF_cState Sturbchem(LES3DTF_pState &Wc,
                                   const LES3DTF_pState &dWdx,
                                   const LES3DTF_pState &dWdy,
                                   const LES3DTF_pState &dWdz,
                                   const int Flow_Type,
                                   const double &Volume);

//@}

/** @name Operators */
/*        --------- */
//@{
   //! Index operator
   double &operator[](int index);
   //! Index operator
   const double &operator[](int index) const;

   //! Binary addition operator
   LES3DTF_pState operator +(const LES3DTF_pState &W) const;

   //! Binary subtraction operator
   LES3DTF_pState operator -(const LES3DTF_pState &W) const;

   //! Binary multiplication operator
   LES3DTF_pState operator *(const double &a) const;

   //! Binary multiplication operator
   friend LES3DTF_pState operator *(const double &a, const LES3DTF_pState &W);

   //! Binary multiplication operator
   double operator *(const LES3DTF_pState &W) const;

   //! Binary division operator
   LES3DTF_pState operator /(const double &a) const;

   //! Binary vector product operator
   LES3DTF_pState operator ^(const LES3DTF_pState &W) const;
   
   //! Assignment operator
   LES3DTF_pState& operator =(const LES3DTF_pState &W);

   //! Shortcut addition operator
   LES3DTF_pState& operator +=(const LES3DTF_pState &W);

   //! Shortcut subtraction operator
   LES3DTF_pState& operator -=(const LES3DTF_pState &W);
  
   //! Equal relational operator
   friend int operator ==(const LES3DTF_pState &W1,
                          const LES3DTF_pState &W2);

   //! Not equal relational operator
   friend int operator !=(const LES3DTF_pState &W1,
                          const LES3DTF_pState &W2);
  
   //! Output stream operator
   friend ostream& operator << (ostream &out_file,
                                const LES3DTF_pState &W);

   //! Input stream operator
   friend istream& operator >> (istream &in_file,
                                LES3DTF_pState &W);
//@}
};
  
/*!
 * Class: LES3DTF_cState
 *
 * \brief Conserved state solution class for 3D Favre-filtered Navier-Stokes
 *        equations associated with solution of premixed compressible turbulent 
 *        combusting flows of a thermally perfect using a large eddy simulation
 *        (LES) technique in conjunction with a thickened flame
 *        subfilter scale model.
 *
 * Member functions
 *  - rho      -- Return mixture density (kg/m^3)
 *  - rhov     -- Return mixture momentum (kg/(m^2-s))   
 *  - E        -- Return mixture total energy (J/kg)
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
class LES3DTF_cState : public NavierStokes3D_ThermallyPerfect_cState {
  public:
   static double           _laminar_flame_speed;  //!< Propagation speed of laminar premixed flame
   static double       _laminar_flame_thickness;  //!< Thickness of laminar premixed flame
   static double                       _TFactor;  //!< Maximum thickening factor
   static double                  _filter_width;  //!< Constant filter width

   double                  rhok; //!< Subfilter scale turbulent kinetic energy (kg/m-s^2)

   PowerLaw               flame; //!< Object containing the thickening and wrinkling factors. 

   static double           Mref; //!< Reference Mach number for low-Mach-number precondtioning (normally set to incoming freestream Mach)


/** @name Constructors and destructors */
/*        ---------------------------- */
//@{
   //! Default creation constructor (assign default values)
   LES3DTF_cState(void) : NavierStokes3D_ThermallyPerfect_cState(), rhok(ZERO) { }
   
   //! Constructor from base class (allows return of derived type)
   LES3DTF_cState(const NavierStokes3D_ThermallyPerfect_cState &U1) : 
     NavierStokes3D_ThermallyPerfect_cState(U1), rhok(ZERO) { }

   //! Constructor from base class
     LES3DTF_cState(const NavierStokes3D_ThermallyPerfect_cState &U1,
                    const double &dk) : 
     NavierStokes3D_ThermallyPerfect_cState(U1), rhok(dk) { }

   //! Assignment constructor
   LES3DTF_cState(const double &value) : 
     NavierStokes3D_ThermallyPerfect_cState(value), rhok(value) { }
   
   //! Assignment constructor
   LES3DTF_cState(const double &d, 
                  const Vector3D &V,
                  const double &En, 
                  const double &dk) :
     NavierStokes3D_ThermallyPerfect_cState(d, V, En), rhok(dk) { }

   //! Assignment constructor
   LES3DTF_cState(const double &d, 
                  const Vector3D &V,
                  const double &En, 
                  const double &dk,
                  const double &frac) :
     NavierStokes3D_ThermallyPerfect_cState(d, V, En, frac), rhok(dk) { }

   //! Assignment constructor
   LES3DTF_cState(const double &d, 
                  const Vector3D &dV,
                  const double &En,
                  const double &dk,
                  const Species *rhomfrac) :
     NavierStokes3D_ThermallyPerfect_cState(d, dV, En, rhomfrac), rhok(dk) { }

   //! Assignment constructor
   LES3DTF_cState(const double &d, 
                  const double &vx,
                  const double &vy, 
                  const double &vz,
                  const double &En, 
                  const double &dk) :
     NavierStokes3D_ThermallyPerfect_cState(d, vx, vy, vz, En), rhok(dk) { }
   
   //! Assignment constructor
   LES3DTF_cState(const double &d, 
                  const double &vx, 
                  const double &vy, 
                  const double &vz, 
                  const double &En, 
                  const double &dk,
                  const double &rhomfrac) :
     NavierStokes3D_ThermallyPerfect_cState(d, vx, vy, vz, En, rhomfrac), rhok(dk) { }

   //! Assignment constructor
   LES3DTF_cState(const double &d, 
                  const double &vx,
                  const double &vy, 
                  const double &vz,
                  const double &En, 
                  const double &dk,
                  const Species *rhomfrac) :
     NavierStokes3D_ThermallyPerfect_cState(d, vx, vy, vz, En, rhomfrac), rhok(dk) { }
   
   //! Copy constructor (this is needed for the operator overload returns)
   LES3DTF_cState(const LES3DTF_cState &U) {
     rhospec_null();
     rho = DENSITY_STDATM;
     set_initial_values(); 
     Copy(U);
   }
                           
   //! Default destructor
   ~LES3DTF_cState(void) {
      Deallocate();
   }
//@}
                           
/** @name Some useful operators */
/*        --------------------- */
//@{
   //! Assign the species data (needs to be called only once as it is static)
   void set_species_data(const int &n,
                         const string *S,
                         const char *PATH,
                         const int &debug, 
                         const double &Mr, 
                         const double* Sc,
                         const int &trans_data) {
      Euler3D_ThermallyPerfect_cState::set_species_data(n, S, PATH, debug, Mr, Sc, trans_data);
      num_vars = NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA + ns;
   }

   //! Assign modelling data (needs to be called only once as it is static)
   void set_modelling_data(const double &laminar_flame_speed,
			   const double &laminar_flame_thickness,
			   const double &thickening_factor,
			   const double &filter_width) {
     _laminar_flame_speed = laminar_flame_speed;
     _laminar_flame_thickness = laminar_flame_thickness;
     _TFactor = thickening_factor;
     _filter_width = filter_width;    
   }

   //! Copy solution state (cheaper than = operator)
   void Copy(const LES3DTF_cState &U);

   //! Assigns a vacuum solution state
   void Vacuum(void) { 
     rho = ZERO; rhov.zero(); E=ZERO; rhok = ZERO;
     for(int i=0; i<ns; ++i) rhospec[i].Vacuum();
   }

   //! Sum N-1 species
   double sum_species(void) {
     double sum(ZERO);
     for(int i=0; i<ns-1; ++i) sum += rhospec[i].c;
     return sum/rho;
   }

   //! Check for physical validity of scalars
   void Realizable_Scalar_Check(void) {
     if ( rhok < MICRO ) { rhok = ZERO; }
   }

   //! Check for physical validity of the solution vector
   bool Realizable_Solution_Check(void);
//@}

/** @name Thermodynamic and other state functions */
/*        --------------------------------------- */
//@{
   //! Mixture absolute (sensible+chemical) internal energy
   double e(void) const;            

   //! Mixture sensible internal energy
   double es(void) const;           

   //! Mixture absolute specific enthalpy  
   double h(void) const; 
   //! Mixture absolute specific enthalpy  
   double h(const double &T) const;  

   //! Mixture sensible specific enthalpy
   double hs(void) const;
   //! Mixture sensible specific enthalpy
   double hs(const double &T) const; 

   //! Mixture temperature
   double T(void) const;

   //! Turbulence modified pressure
   double p_t(void) const;

   //! Mixture sound speed
   double a(void) const;

   //! Mixture sound speed including turbulent kinetic energy
   double a_t(void) const;

   //! Turbulent kinetic energy
   double k(void) const;

//@}

/** @name Turbulent transport coefficients */
/*        -------------------------------- */
//@{
   //! LES filter width
   double filter_width(const double &Volume) const;
   double filter_width(void) const;

   //! Absolute value of strain rate
/*    double abs_strain_rate(const LES3DTF_pState &dWdx,  */
/*                           const LES3DTF_pState &dWdy,  */
/*                           const LES3DTF_pState &dWdz) const; */
//@}

/** @name Primitive solution state */ 
/*        ------------------------ */
//@{
   //! Returns primitive solution state
   LES3DTF_pState W(void);
   //! Returns primitive solution state
   LES3DTF_pState W(void) const;
   //! Returns primitive solution state
   LES3DTF_pState W(const LES3DTF_cState &U) const;
//@}

/** @name Numerical Flux Functions */
/*        ------------------------ */
//@{
   //! Returns rotated conserved state aligned with x-axis in the norm_dir
   LES3DTF_cState Rotate(const Vector3D &norm_dir) const;

   //! Returns un-rotated conserved state aligned with x-axis of the global problem
   LES3DTF_cState RotateBack(const Vector3D &norm_dir) const;
//@}

/** @name Operators */
/*        --------- */
//@{
   //! Index operator
   double &operator[](int index);
   //! Index operator
   const double &operator[](int index) const;

   //! Binary addition operator
   LES3DTF_cState operator +(const LES3DTF_cState &U) const;

   //! Binary subtraction operator
   LES3DTF_cState operator -(const LES3DTF_cState &U) const;

   //! Binary multiplication operator
   LES3DTF_cState operator *(const double &a) const;

   //! Binary multiplication operator
   friend LES3DTF_cState operator *(const double &a, const LES3DTF_cState &U);

   //! Binary multiplication operator
   double operator *(const LES3DTF_cState &U) const;

   //! Binary division operator
   LES3DTF_cState operator /(const double &a) const;

   //! Binary vector product operator
   LES3DTF_cState operator ^(const LES3DTF_cState &U) const;
   
   //! Assignment operator
   LES3DTF_cState& operator =(const LES3DTF_cState &U);

   //! Shortcut addition operator
   LES3DTF_cState& operator +=(const LES3DTF_cState &U);

   //! Shortcut subtraction operator
   LES3DTF_cState& operator -=(const LES3DTF_cState &U);
  
   //! Equal relational operator
   friend int operator ==(const LES3DTF_cState &U1,
                          const LES3DTF_cState &U2);

   //! Not equal relational operator
   friend int operator !=(const LES3DTF_cState &U1,
                          const LES3DTF_cState &U2);
  
   //! Output stream operator
   friend ostream& operator << (ostream &out_file,
                                const LES3DTF_cState &U);

   //! Input stream operator
   friend istream& operator >> (istream &in_file,
                                LES3DTF_cState &U);
//@}
};

/***************************************
 * LES3DTF_pState member functions    *
 ***************************************/

/***************************************
 * LES3DTF_pState -- Index operators. *
 ***************************************/
inline double& LES3DTF_pState::operator[](int index) {
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
   default :
     if (index <= num_vars) {
        return spec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].c;
      } else if (index == num_vars+1) {
	return flame.WF;
      } else {
	return flame.TF;
      }
   };
}

inline const double& LES3DTF_pState::operator[](int index) const {
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
   default :
     if (index <= num_vars) {
        return spec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].c;     
      } else if (index == num_vars+1) {
	return flame.WF;
      } else {
	return flame.TF;
      }
   };
}

/***************************************************************************
 * LES3DTF_pState -- Binary arithmetic operators.                          *
 ***************************************************************************/
//------------------ Addition ------------------------//
inline LES3DTF_pState LES3DTF_pState::operator +(const LES3DTF_pState &W) const { 
   LES3DTF_pState Temp(*this);
   Temp += W;
   return Temp;
}

//------------------ Subtraction ------------------------//
inline LES3DTF_pState LES3DTF_pState::operator -(const LES3DTF_pState &W) const {
   LES3DTF_pState Temp(*this);
   Temp -= W;
   return Temp;
}

//---------------- Scalar Multiplication ------------------//
inline LES3DTF_pState LES3DTF_pState::operator *(const double &a) const {
   LES3DTF_pState Temp(*this);
   Temp.rho = rho*a;  
   Temp.v = v*a; 
   Temp.p = p*a;
   Temp.k = k*a; 
   for (int i = 0; i < ns; ++i) {
      Temp.spec[i] = spec[i]*a;
   } /* endfor */
   Temp.flame = flame*a; 
   return(Temp);
}

inline LES3DTF_pState operator *(const double &a, 
                                 const LES3DTF_pState &W) {
   LES3DTF_pState Temp;
   //Temp.Copy(W);
   Temp.rho = W.rho*a;  
   Temp.v = W.v*a; 
   Temp.p = W.p*a;
   Temp.k = W.k*a; 
   for(int i = 0; i < W.ns; ++i) {
      Temp.spec[i] = W.spec[i]*a;
   } /* endfor */
   Temp.flame = W.flame*a;
   return(Temp);
}

//--------------- Scalar Division ------------------------//
inline LES3DTF_pState LES3DTF_pState::operator /(const double &a) const {
   LES3DTF_pState Temp(*this);
   Temp.rho = rho/a; 
   Temp.v = v/a; 
   Temp.p = p/a; 
   Temp.k = k/a; 
   for (int i = 0; i < ns; ++i) {
      Temp.spec[i] = spec[i]/a; 
   } /* endfor */
   Temp.flame = flame/a;
   return(Temp);
}

//----------------- Inner Product ------------------------//
inline double LES3DTF_pState::operator *(const LES3DTF_pState &W) const {
   double sum(0.0);
   for (int i=0; i < ns; ++i) {
      sum += spec[i]*W.spec[i];
   } /* endfor */
   return (rho*W.rho + v*W.v + p*W.p + k*W.k  + sum);
}

//----------- Solution state product operator ------------//
inline LES3DTF_pState LES3DTF_pState::operator ^(const LES3DTF_pState &W) const {
   LES3DTF_pState Temp(*this);
   Temp.rho = rho*W.rho;
   Temp.v.x = v.x*W.v.x;
   Temp.v.y = v.y*W.v.y;
   Temp.v.z = v.z*W.v.z;
   Temp.p = p*W.p;
   Temp.k = k*W.k;
   for (int i = 0; i < ns; ++i){
      Temp.spec[i] = spec[i]*W.spec[i];
   } /* endfor */
   Temp.flame.WF = flame.WF*W.flame.WF;
   Temp.flame.TF = flame.TF*W.flame.TF;
   return(Temp);
}

//----------------- Assignment ----------------------------//
inline LES3DTF_pState& LES3DTF_pState::operator =(const LES3DTF_pState &W) {
   //self assignment protection
   if (this != &W){   
     //copy assignment
     Copy(W);
   } 
   return (*this);
}

/***************************************************************************
 * LES3DTF_pState -- Shortcut arithmetic operators.                        *
 ***************************************************************************/
inline LES3DTF_pState& LES3DTF_pState::operator +=(const LES3DTF_pState &W) {
   rho += W.rho;
   v += W.v; 
   p += W.p; 
   k += W.k;
   for (int i = 0; i < ns; ++i) {
      spec[i].c += W.spec[i].c;
   } /* endfor */
   flame += W.flame;
   return (*this);
}

inline LES3DTF_pState& LES3DTF_pState::operator -=(const LES3DTF_pState &W) {
   rho -= W.rho;
   v -= W.v;
   p -= W.p;
   k -= W.k;
   for (int i = 0; i < ns; ++i) {
      spec[i].c -= W.spec[i].c;
   } /* endfor */
   flame -= W.flame;
   return (*this); 
}

/***************************************************************************
 * LES3DTF_pState -- Unary arithmetic operators.                           *
 ***************************************************************************/
inline LES3DTF_pState operator -(const LES3DTF_pState &W) {
   Species *spt= new Species[W.ns];
   for (int i = 0; i < W.ns; ++i) {
      spt[i] = -W.spec[i]; 
   } /* endfor */ 
   LES3DTF_pState Temp(-W.rho, -W.v, -W.p, -W.k, spt);
   Temp.flame = -W.flame;
   delete[] spt;
   return(Temp);
}

/***************************************************************************
 * LES3DTF_pState -- Relational operators.                                 *
 ***************************************************************************/
inline int operator ==(const LES3DTF_pState &W1, 
                       const LES3DTF_pState &W2) {
   bool Temp;
   for (int i = 0; i < W1.ns; ++i) {
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
           Temp == true);
}

inline int operator !=(const LES3DTF_pState &W1, 
                       const LES3DTF_pState &W2) {
   bool Temp = true;
   for (int i = 0; i < W1.ns; ++i) {
      if (W1.spec[i] != W2.spec[i]) {
         Temp = false;
         break;
      } /* endif */
   } /* endfor */
   return (W1.rho != W2.rho || 
           W1.v != W2.v || 
           W1.p != W2.p || 
           W1.k != W2.k || 
           Temp != true);
}

/****************************************************
 * LES3DTF_pState -- Input-output operators.        *
 ****************************************************/
inline ostream &operator << (ostream &out_file, 
                             const LES3DTF_pState &W) {
   out_file.precision(10);
   out_file.setf(ios::scientific);
   out_file << " " << W.rho  << " " << W.v.x << " " << W.v.y 
            << " " << W.v.z << " " << W.p << " "<<W.k;
   for (int i = 0; i < W.ns; ++i) {
      out_file<<" "<<W.spec[i];
   } /* endfor */
   out_file << " " << W.flame;
   out_file.unsetf(ios::scientific);
   return (out_file);
}

inline istream &operator >> (istream &in_file, 
                             LES3DTF_pState &W) {
   in_file.setf(ios::skipws);
   in_file >> W.rho >> W.v.x >> W.v.y >> W.v.z >> W.p>> W.k;
   //W.set_initial_values();
   for (int i = 0; i < W.ns; ++i) {
      in_file>>W.spec[i];
   } /* endfor */
   in_file >> W.flame;
   in_file.unsetf(ios::skipws);
   return (in_file);
}

/****************************************************************
 * LES3DTF_cState member functions                              *
 ****************************************************************/

/****************************************************************
 * LES3DTF_cState -- Index operators.                           *
 ****************************************************************/
//index operators
inline double& LES3DTF_cState::operator[](int index) {
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
   default :
     if (index <= num_vars) {
        return rhospec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].c;     
      } else if (index == num_vars+1) {
	return flame.WF;
      } else {
	return flame.TF;
      }
   };
}

inline const double& LES3DTF_cState::operator[](int index) const {
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
   default :
     if (index <= num_vars) {
        return rhospec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3DTF_VAR_EXTRA+1)].c; 
      } else if (index == num_vars+1) {
	return flame.WF;
      } else {
	return flame.TF;
      }
   };
}

/***************************************************************************
 * LES3DTF_cState -- Binary arithmetic operators.                          *
 ***************************************************************************/
//------------------ Addition ------------------------//
inline LES3DTF_cState LES3DTF_cState::operator +(const LES3DTF_cState &U) const { 
   LES3DTF_cState Temp(*this);
   Temp += U;
   return Temp;
}

//------------------ Subtraction ------------------------//
inline LES3DTF_cState LES3DTF_cState::operator -(const LES3DTF_cState &U) const {
   LES3DTF_cState Temp(*this);
   Temp -= U;
   return Temp;
}

//---------------- Scalar Multiplication ------------------//
inline LES3DTF_cState LES3DTF_cState::operator *(const double &a) const {
   LES3DTF_cState Temp(*this);
   Temp.rho = rho*a;  
   Temp.rhov = rhov*a; 
   Temp.E = E*a;
   Temp.rhok = rhok*a; 
   for (int i = 0; i < ns; ++i) {
      Temp.rhospec[i] = rhospec[i]*a;
   } /* endfor */
   Temp.flame = flame*a;
   return(Temp);
}

inline LES3DTF_cState operator *(const double &a, 
                                 const LES3DTF_cState &U) {
   LES3DTF_cState Temp;
   Temp.rho = U.rho*a;  
   Temp.rhov = U.rhov*a; 
   Temp.E = U.E*a;
   Temp.rhok = U.rhok*a; 
   for (int i = 0; i < U.ns; ++i) {
      Temp.rhospec[i] = U.rhospec[i]*a;
   } /* endfor */
   Temp.flame = U.flame*a;
   return(Temp);
}

//--------------- Scalar Division ------------------------//
inline LES3DTF_cState LES3DTF_cState::operator /(const double &a) const {
   LES3DTF_cState Temp(*this);
   Temp.rho = rho/a; 
   Temp.rhov = rhov/a; 
   Temp.E = E/a;
   Temp.rhok = rhok/a; 
   for (int i = 0; i < ns; ++i) {
      Temp.rhospec[i] = rhospec[i]/a; 
   } /* endfor */
   Temp.flame = flame/a;
   return(Temp);
}

//----------------- Inner Product ------------------------//
inline double LES3DTF_cState::operator *(const LES3DTF_cState &U) const {
   double sum(0.0);
   for (int i = 0; i < ns; ++i) {
      sum += rhospec[i]*U.rhospec[i];
   } /* endfor */
   return (rho*U.rho + rhov*U.rhov + E*U.E + rhok*U.rhok + sum);
}

//----------- Solution state product operator ------------//
inline LES3DTF_cState LES3DTF_cState::operator ^(const LES3DTF_cState &U) const {
   LES3DTF_cState Temp(*this);
   Temp.rho = rho*U.rho;
   Temp.rhov.x = rhov.x*U.rhov.x;
   Temp.rhov.y = rhov.y*U.rhov.y;
   Temp.rhov.z = rhov.z*U.rhov.z;
   Temp.E = E*U.E;
   Temp.rhok = rhok*U.rhok;
   for (int i = 0; i < ns; ++i) {
      Temp.rhospec[i] = rhospec[i]*U.rhospec[i];
   } /* endfor */
   Temp.flame.WF = flame.WF*U.flame.WF;
   Temp.flame.TF = flame.TF*U.flame.TF;
   return(Temp);
}

//----------------- Assignment ----------------------------//
inline LES3DTF_cState& LES3DTF_cState::operator =(const LES3DTF_cState &U) {
   //self assignment protection
   if (this != &U) {   
     //copy assignment
     Copy(U);
   } 
   return (*this);
}

/****************************************************
 * LES3DTF_cState -- Shortcut arithmetic operators. *
 ****************************************************/
inline LES3DTF_cState& LES3DTF_cState::operator +=(const LES3DTF_cState &U) {
   rho += U.rho;
   rhov += U.rhov; 
   E += U.E;
   rhok += U.rhok;
   for (int i = 0; i < ns; ++i) {
      rhospec[i].c += U.rhospec[i].c;
   } /* endfor */
   flame += U.flame;
   return (*this);
}

inline LES3DTF_cState& LES3DTF_cState::operator -=(const LES3DTF_cState &U) {
   rho -= U.rho;
   rhov -= U.rhov;
   E -= U.E;
   rhok -= U.rhok;
   for (int i = 0; i < ns; ++i) {
     rhospec[i].c -= U.rhospec[i].c;
   } /* endfor */
   flame -= U.flame;
   return (*this); 
}

/****************************************************
 * LES3DTF_cState -- Unary arithmetic operators.    *
 ****************************************************/
inline LES3DTF_cState operator -(const LES3DTF_cState &U) {
   Species *spt= new Species[U.ns];
   for (int i=0; i < U.ns; ++i) {
      spt[i] = -U.rhospec[i]; 
   } /* endfor */
   LES3DTF_cState Temp(-U.rho,-U.rhov,-U.E, -U.rhok, spt);
   Temp.flame = -U.flame;
   delete[] spt;
   return(Temp);
}

/****************************************************
 * LES3DTF_cState -- Relational operators.          *
 ****************************************************/
inline int operator ==(const LES3DTF_cState &U1, 
                       const LES3DTF_cState &U2) {
   bool Temp;
   for (int i = 0; i < U1.ns; ++i) {
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
           Temp == true);
}

inline int operator !=(const LES3DTF_cState &U1, 
                       const LES3DTF_cState &U2) {
   bool Temp = true;
   for (int i = 0; i < U1.ns; ++i) {
      if( U1.rhospec[i] != U2.rhospec[i] ){
         Temp = false;
         break;
      } /* endif */
   } /* endfor */
   return (U1.rho != U2.rho || 
           U1.rhov != U2.rhov || 
           U1.E != U2.E ||
           U1.rhok != U2.rhok || 
           Temp != true);
}

/****************************************************
 * LES3DTF_cState -- Input-output operators.        *
 ****************************************************/
inline ostream &operator << (ostream &out_file, 
                             const LES3DTF_cState &U) {
  //out_file.precision(20);
  out_file.setf(ios::scientific);
  out_file << " " << U.rho  << " " << U.rhov.x << " " << U.rhov.y << " " 
           << U.rhov.z << " " << U.E << " " << U.rhok;
  for (int i = 0; i < U.ns; ++i) {
     out_file<<" "<<U.rhospec[i];
  } /* endfor */
  out_file << " " << U.flame;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, 
                             LES3DTF_cState &U) {
   in_file.setf(ios::skipws);
   in_file >> U.rho >> U.rhov.x >> U.rhov.y >> U.rhov.z >> U.E >> U.rhok;
   for (int i = 0; i < U.ns; ++i) {
     in_file>>U.rhospec[i]; 
   } /* endfor */
   in_file >> U.flame;
   in_file.unsetf(ios::skipws);
   return (in_file);
}

#endif // _LES3DTF_STATE_INCLUDED
