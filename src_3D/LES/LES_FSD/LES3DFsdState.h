/*! \file LES3DFsdState.h
 * 	\brief	Header file defining the Favre-filtered Navier-Stokes solution state 
 *              classes associated with solution of premixed compressible turbulent 
 *              combusting flows of a thermally perfect using a large eddy simulation
 *              (LES) technique in conjunction with a flame surface density (FSD)
 *              subfilter scale model.
 * .
 */

#ifndef _LES3DFSD_STATE_INCLUDED 
#define _LES3DFSD_STATE_INCLUDED

/* Include header file for base solution classes from which the classes are derived. */

#ifndef _NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED
#include "../../NavierStokes/NavierStokes3DThermallyPerfectState.h"
#endif  //NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED

/* Include turbulence modelling header file. */

#ifndef _TURBULENCE_MODELLING_INCLUDED
#include "../../TurbulenceModelling/TurbulenceModelling.h"
#endif  //TURBULENCE_MODELLING_INCLUDED

/* Define the classes. */

class LES3DFsd_cState;
class LES3DFsd_pState;

#define NUM_LES3D_VAR_EXTRA 3  //k,C,Fsd

/*!
 * Class: LES3DFsd_pState
 *
 * \brief Primitive state solution class for 3D Favre-filtered Navier-Stokes
 *        equations associated with solution of premixed compressible turbulent 
 *        combusting flows of a thermally perfect using a large eddy simulation
 *        (LES) technique in conjunction with a flame surface density (FSD)
 *        subfilter scale model.
 *
 * Member functions
 *  - rho              -- Return density (kg/m^3)
 *  - v                -- Return flow velocity (m/s)
 *  - p                -- Return pressure (Pa, N/m^2)
 *  - p_t               -- Return turbulence modified pressure (Pa, N/m^2)
 *  - k                -- Return subfilter scale turbulent kinetic energy (m^2/s^2)
 *  - C                -- Return reaction rate progress variable (0 <= C <= 1)
 *  - Fsd              -- Return premixed flame surface density
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
 *  - Sturbulence      -- Return source terms associated with turbulence modelling
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
class LES3DFsd_pState : public NavierStokes3D_ThermallyPerfect_pState {
  protected:
   static double        _fuel_equivalence_ratio;  //!< Fuel equivalence ratio
   static double    _unburnt_fuel_mass_fraction;  //!< Mass fraction of unburnt fuel (premixed flame)
   static double             _reactants_density;  //!< Reactants Density
   static double           _laminar_flame_speed;  //!< Propagation speed of laminar premixed flame
   static double       _laminar_flame_thickness;  //!< Thickness of laminar premixed flame
   static double   _adiabatic_flame_temperature;  //!< Adiabitc flame temperature

  public:
   double                     C; //!< Reaction rate progress variable (0 <= C <= 1)
   double                   Fsd; //!< Flame surface density (FSD) of premixed flame
   double                     k; //!< Subfilter scale turbulent kinetic energy (m^2/s^2)

   static double           Mref; //!< Reference Mach number for low-Mach-number precondtioning (normally set to incoming freestream Mach)

/** @name Constructors and destructors */
/*        ---------------------------- */
//@{
   //! Default creation constructor (assign default values)
   LES3DFsd_pState(void) : NavierStokes3D_ThermallyPerfect_pState(), C(ONE), Fsd(MILLION), k(ZERO) {
     premixed_mfrac(); }
   
   //! Constructor from base class (allows return of derived type)
   LES3DFsd_pState(const NavierStokes3D_ThermallyPerfect_pState &W1) : 
     NavierStokes3D_ThermallyPerfect_pState(W1), C(ONE), Fsd(MILLION), k(ZERO) { }

   //! Constructor from base class
   LES3DFsd_pState(const NavierStokes3D_ThermallyPerfect_pState &W1,
                   const double &cc,
                   const double &ff, 
                   const double &kk) : 
     NavierStokes3D_ThermallyPerfect_pState(W1), C(cc), Fsd(ff), k(kk) { 
     premixed_mfrac(); }

   //! Assignment constructor
   LES3DFsd_pState(const double &value):
     NavierStokes3D_ThermallyPerfect_pState(value), C(value), Fsd(value), k(value) {
     premixed_mfrac(); }
   
   //! Assignment constructor
   LES3DFsd_pState(const double &d, 
                   const Vector3D &V,
                   const double &pre, 
                   const double &cc, 
                   const double &ff, 
                   const double &kk):
     NavierStokes3D_ThermallyPerfect_pState(d, V, pre), C(cc), Fsd(ff), k(kk) {
     premixed_mfrac(); }

   //! Assignment constructor
   LES3DFsd_pState(const double &d,
                   const double &vx,
                   const double &vy,
                   const double &vz,
                   const double &pre,
                   const double &cc, 
                   const double &ff,
                   const double &kk):
     NavierStokes3D_ThermallyPerfect_pState(d, vx, vy, vz, pre), C(cc), Fsd(ff), k(kk) { 
     premixed_mfrac(); }

   //! Copy constructor (this is needed for the operator overload returns)
   LES3DFsd_pState(const  LES3DFsd_pState &W) {
      species_null(); rho = DENSITY_STDATM; set_initial_values(); Copy(W);
   }

   //! Default destructor
   ~LES3DFsd_pState(void) {
      Deallocate();
   }
//@}

/** @name Some useful operators */
/*        --------------------- */
//@{
   //! Copy solution state (cheaper than = operator)
   void Copy(const LES3DFsd_pState &W);

   //! Assigns a vacuum solution state
   void Vacuum(void) { 
     rho = ZERO, v.zero(); p = ZERO; C = ZERO; Fsd = ZERO; k = ZERO;
   }

   //! Check for physical validity of the progress variable
   void Realizable_C_Check(void) {
      if (C > ONE) {
        C = ONE;
      } else if (C < ZERO) {
        C = ZERO;
      } /* endif */
   }

   //! Check for physical validity of the progress variable
   void Realizable_Mass_Fraction_Check(void) {
      if (C > ONE) {
        C = ONE;
        premixed_mfrac();
      } else if (C < ZERO) {
        C = ZERO;
        premixed_mfrac();
      } /* endif */
   }

   //! Check for physical validity of scalars
   bool Realizable_Scalar_Check(void) {
     double LOCAL_TOL = MICRO;
     if ( C < LOCAL_TOL ) { C = ZERO; premixed_mfrac(); }
     if ( Fsd < LOCAL_TOL ) { Fsd = ZERO; }
     if ( k < LOCAL_TOL ) { k = ZERO; }
     if ( C > 0.999 || C < 0.001 ) { Fsd = ZERO; }
     return true;
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
   double mu_t(const LES3DFsd_pState &dWdx, 
               const LES3DFsd_pState &dWdy, 
               const LES3DFsd_pState &dWdz,
               const int Flow_Type, 
               const double &Volume);

   //! Turbulent thermal conductivity
   double kappa_t(const LES3DFsd_pState &dWdx, 
                  const LES3DFsd_pState &dWdy, 
                  const LES3DFsd_pState &dWdz,
                  const int Flow_Type, 
                  const double &Volume);      

   //! Species turbulent diffusion coefficient
   double Ds_t(const int i,
               const LES3DFsd_pState &dWdx, 
               const LES3DFsd_pState &dWdy, 
               const LES3DFsd_pState &dWdz,
               const int Flow_Type, 
               const double &Volume);

   //! Species turbulent diffusion coefficient
   double Ds_t(const int i,
               const double &mu_t_temp);

   //! Turbulent Prandtl number
   double Pr_t(void);

   //! Turbulent Schmidt number
   double Sc_t(void);      

   //! Turbulent Lewis number
   double Le_t(void); 

   //! LES filter width
   double filter_width(const double &Volume) const;
//@}

/** @name Subfilter scale turbulent stress tensor and heat flux vector */
/*        ------------------------------------------------------------ */
//@{
   //! Returns subfilter scale turbulent stress tensor 
   Tensor3D tau_t(const LES3DFsd_pState &dWdx, 
                  const LES3DFsd_pState &dWdy,
                  const LES3DFsd_pState &dWdz,
                  const int Flow_Type, 
                  const double &Volume);

   //! Returns subfilter scale turbulent stress tensor 
   Tensor3D tau_t(const double &mu_t_temp, 
                  const LES3DFsd_pState &dWdx, 
                  const LES3DFsd_pState &dWdy,
                  const LES3DFsd_pState &dWdz,
                  const int Flow_Type, 
                  const double &Volume);

   //! Returns components of subfilter scale turbulent stress tensor in x-direction 
   Tensor3D tau_t_x(const LES3DFsd_pState &dWdx, 
                    const LES3DFsd_pState &dWdy,
                    const LES3DFsd_pState &dWdz,
                    const int Flow_Type, 
                    const double &Volume);

   //! Returns components of subfilter scale turbulent stress tensor in x-direction
   Tensor3D tau_t_x(const double &mu_t_temp, 
                    const LES3DFsd_pState &dWdx, 
                    const LES3DFsd_pState &dWdy,
                    const LES3DFsd_pState &dWdz,
                    const int Flow_Type, 
                    const double &Volume);

   //! Returns components of subfilter scale turbulent stress tensor in y-direction 
   Tensor3D tau_t_y(const LES3DFsd_pState &dWdx, 
                    const LES3DFsd_pState &dWdy,
                    const LES3DFsd_pState &dWdz,
                    const int Flow_Type, 
                    const double &Volume);

   //! Returns components of subfilter scale turbulent stress tensor in y-direction
   Tensor3D tau_t_y(const double &mu_t_temp, 
                    const LES3DFsd_pState &dWdx, 
                    const LES3DFsd_pState &dWdy,
                    const LES3DFsd_pState &dWdz,
                    const int Flow_Type, 
                    const double &Volume);

   //! Returns components of subfilter scale turbulent stress tensor in z-direction 
   Tensor3D tau_t_z(const LES3DFsd_pState &dWdx, 
                    const LES3DFsd_pState &dWdy,
                    const LES3DFsd_pState &dWdz,
                    const int Flow_Type, 
                    const double &Volume);

   //! Returns components of subfilter scale turbulent stress tensor in z-direction
   Tensor3D tau_t_z(const double &mu_t_temp, 
                    const LES3DFsd_pState &dWdx, 
                    const LES3DFsd_pState &dWdy,
                    const LES3DFsd_pState &dWdz,
                    const int Flow_Type, 
                    const double &Volume);

   //! Returns subfilter scale turbulent heat flux vector 
   Vector3D q_t(const LES3DFsd_pState &dWdx, 
                const LES3DFsd_pState &dWdy,
                const LES3DFsd_pState &dWdz,
                const int Flow_Type, 
                const double &Volume);

   //! Returns subfilter scale turbulent heat flux vector 
   Vector3D q_t(const double &kappa_t_temp,
                const LES3DFsd_pState &dWdx, 
                const LES3DFsd_pState &dWdy,
                const LES3DFsd_pState &dWdz,
                const int Flow_Type, 
                const double &Volume);

   //! Returns component of subfilter scale turbulent heat flux vector in x-direction
   Vector3D q_t_x(const LES3DFsd_pState &dWdx, 
                  const LES3DFsd_pState &dWdy,
                  const LES3DFsd_pState &dWdz,
                  const int Flow_Type, 
                  const double &Volume);

   //! Returns component of subfilter scale turbulent heat flux vector in x-direction
   Vector3D q_t_x(const double &kappa_t_temp,
                  const LES3DFsd_pState &dWdx, 
                  const LES3DFsd_pState &dWdy,
                  const LES3DFsd_pState &dWdz,
                  const int Flow_Type, 
                  const double &Volume);

   //! Returns component of subfilter scale turbulent heat flux vector in y-direction
   Vector3D q_t_y(const LES3DFsd_pState &dWdx, 
                  const LES3DFsd_pState &dWdy,
                  const LES3DFsd_pState &dWdz,
                  const int Flow_Type, 
                  const double &Volume);

   //! Returns component of subfilter scale turbulent heat flux vector in y-direction
   Vector3D q_t_y(const double &kappa_t_temp,
                  const LES3DFsd_pState &dWdx, 
                  const LES3DFsd_pState &dWdy,
                  const LES3DFsd_pState &dWdz,
                  const int Flow_Type, 
                  const double &Volume);

   //! Returns component of subfilter scale turbulent heat flux vector in z-direction
   Vector3D q_t_z(const LES3DFsd_pState &dWdx, 
                  const LES3DFsd_pState &dWdy,
                  const LES3DFsd_pState &dWdz,
                  const int Flow_Type, 
                  const double &Volume);

   //! Returns component of subfilter scale turbulent heat flux vector in z-direction
   Vector3D q_t_z(const double &kappa_t_temp,
                  const LES3DFsd_pState &dWdx, 
                  const LES3DFsd_pState &dWdy,
                  const LES3DFsd_pState &dWdz,
                  const int Flow_Type, 
                  const double &Volume);
//@}

/** @name Conserved solution state */ 
/*        ------------------------ */
//@{
   //! Returns conserved solution state
   LES3DFsd_cState U(void);
   //! Returns conserved solution state
   LES3DFsd_cState U(void)const;
   //! Returns conserved solution state
   LES3DFsd_cState U(const LES3DFsd_pState &W); 
//@}
   
/** @name Inviscid Flux Vectors */
/*        --------------------- */
//@{
   //! x-direction inviscid solution flux
   LES3DFsd_cState F(void) ;
   //! x-direction inviscid solution flux
   LES3DFsd_cState F(void) const ;

   //! x-direction inviscid solution flux
   LES3DFsd_cState Fx(void);
   //! x-direction inviscid solution flux
   LES3DFsd_cState Fx(void) const;

   //! y-direction inviscid solution flux
   LES3DFsd_cState Fy(void);
   //! y-direction inviscid solution flux
   LES3DFsd_cState Fy(void) const;

   //! z-direction inviscid solution flux
   LES3DFsd_cState Fz(void);
   //! z-direction inviscid solution flux
   LES3DFsd_cState Fz(void) const;
//@}

/** @name Viscous flux vectors */
/*        --------------------- */
//@{
   //! x-direction viscous solution flux
   LES3DFsd_cState Fv(const LES3DFsd_pState &dWdx,
                      const LES3DFsd_pState &dWdy,
                      const LES3DFsd_pState &dWdz,
                      const int Flow_Type,
                      const double &Volume);

   //! x-direction viscous solution flux
   LES3DFsd_cState Fvx(const LES3DFsd_pState &dWdx,
                       const LES3DFsd_pState &dWdy,
                       const LES3DFsd_pState &dWdz,
                       const int Flow_Type,
                       const double &Volume);

   //! y-direction viscous solution flux
   LES3DFsd_cState Fvy(const LES3DFsd_pState &dWdx,
                       const LES3DFsd_pState &dWdy,
                       const LES3DFsd_pState &dWdz,
                       const int Flow_Type,
                       const double &Volume);

   //! z-direction viscous solution flux
   LES3DFsd_cState Fvz(const LES3DFsd_pState &dWdx,
                       const LES3DFsd_pState &dWdy,
                       const LES3DFsd_pState &dWdz,
                       const int Flow_Type,
                       const double &Volume);
//@}

/** @name Eigenvalue(s) and eigenvector(s) (x-direction) */
/*        ---------------------------------------------- */
//@{
   //! x-direction eigenvalues
   LES3DFsd_pState lambda(void);
   //! x-direction eigenvalues
   LES3DFsd_pState lambda(void) const;

   //! x-direction eigenvalues
   LES3DFsd_pState lambda_x(void);
   //! x-direction eigenvalues
   LES3DFsd_pState lambda_x(void) const;

   //! x-direction conservative eigenvectors
   LES3DFsd_cState rc(const int &index) ;
   //! x-direction conservative eigenvectors
   LES3DFsd_cState rc(const int &index) const;

   //! x-direction conservative eigenvectors
   LES3DFsd_cState rc_x(const int &index) ;
   //! x-direction conservative eigenvectors
   LES3DFsd_cState rc_x(const int &index) const;

   //! x-direction primitive eignenvectors
   LES3DFsd_pState lp(const int &index) ;
   //! x-direction primitive eignenvectors
   LES3DFsd_pState lp(const int &index) const;

   //! x-direction primitive eignenvectors
   LES3DFsd_pState lp_x(const int &index) ;
   //! x-direction primitive eignenvectors
   LES3DFsd_pState lp_x(const int &index) const;  
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
   LES3DFsd_cState rc_x_precon(const int &index, 
                               const double &MR2) const; 

   //! x-direction primitive eigenvectors
   LES3DFsd_pState lp_x_precon(const int &index, 
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

/** @name Semi-Impilicit source term Jacobians */
/*        ------------------------------------ */
//@{
   //! Transform Jacobian for primitive to conservative
   void dWdU(DenseMatrix &dWdU);

   //! Return the source term Jacobian for SemiImplicit time marching
   void SemiImplicitSourceJacobi(const LES3DFsd_pState &dWdx, 
                                 const LES3DFsd_pState &dWdy,
                                 const LES3DFsd_pState &dWdz,
                                 const double &d_dWdx_dW, 
                                 const double &d_dWdy_dW,
                                 const double &d_dWdz_dW,
                                 DenseMatrix &dStdW,
                                 const int Flow_Type,
                                 const double &Volume);
//@}

/** @name Numerical Flux Functions */
/*        ------------------------ */
//@{
   //! Returns Roe-averaged primitive solution state
   static LES3DFsd_pState RoeAverage(const LES3DFsd_pState &Wl,
                                     const LES3DFsd_pState &Wr);

   //! Returns Harten-Lax-van-Leer-Einfeldt (HLLE) x-direction flux
   static LES3DFsd_cState FluxHLLE_x(const LES3DFsd_pState &Wl,
                                     const LES3DFsd_pState &Wr);
   //! Returns Harten-Lax-van-Leer-Einfeldt (HLLE) x-direction flux
   static LES3DFsd_cState FluxHLLE_x(const LES3DFsd_cState &Ul,
                                     const LES3DFsd_cState &Ur);

    //! Returns Harten-Lax-van-Leer-Einfeldt (HLLE) flux in n-direction
   static LES3DFsd_cState FluxHLLE_n(const LES3DFsd_pState &Wl,
                                     const LES3DFsd_pState &Wr,
                                     const Vector3D &norm_dir);
    //! Returns Harten-Lax-van-Leer-Einfeldt (HLLE) flux in n-direction
   static LES3DFsd_cState FluxHLLE_n(const LES3DFsd_cState &Ul,
                                     const LES3DFsd_cState &Ur,
                                     const Vector3D &norm_dir);

   //! Returns Roe flux x-direction flux
   static LES3DFsd_cState FluxRoe_x(const LES3DFsd_pState &Wl, 
                                    const LES3DFsd_pState &Wr);
   //! Returns Roe flux x-direction flux
   static LES3DFsd_cState FluxRoe_x(const LES3DFsd_cState &Ul,
				    const LES3DFsd_cState &Ur);

   //! Returns Roe flux in n-direction
   static LES3DFsd_cState FluxRoe_n(const LES3DFsd_pState &Wl,
                                    const LES3DFsd_pState &Wr,
                                    const Vector3D &norm_dir);

   //! Returns AUSMplus_up flux x-direction flux
   static LES3DFsd_cState FluxAUSMplus_up_x(const LES3DFsd_pState &Wl,
                                            const LES3DFsd_pState &Wr);
   //! Returns AUSMplus_up flux x-direction flux
   static LES3DFsd_cState FluxAUSMplus_up_x(const LES3DFsd_cState &Ul,
                                            const LES3DFsd_cState &Ur);

   //! Returns AUSMplus_up flux in n-direction
   static LES3DFsd_cState FluxAUSMplus_up_n(const LES3DFsd_pState &Wl,
                                            const LES3DFsd_pState &Wr,
                                            const Vector3D &norm_dir);

   //! Returns negative waves speeds (eigenvalues) using Harten entropy fix
   LES3DFsd_pState lambda_minus(const LES3DFsd_pState &lambda_a,
                                const LES3DFsd_pState &lambda_l,
                                const LES3DFsd_pState &lambda_r);

   //! Returns positive waves speeds (eigenvalues) using Harten entropy fix
   LES3DFsd_pState lambda_plus(const LES3DFsd_pState &lambda_a,
                               const LES3DFsd_pState &lambda_l,
                               const LES3DFsd_pState &lambda_r);

   //! HLLE wavespeeds in n-direction given 2 primitive states and a direction
   static Vector2D HLLE_wavespeeds(const LES3DFsd_pState &Wl,
                                   const LES3DFsd_pState &Wr,
                                   const Vector3D &norm_dir);
  
   //! Returns rotated primitive state aligned with local x-axis in the norm_dir
   LES3DFsd_pState Rotate(const Vector3D &norm_dir) const;				  
  
   //! Returns un-rotated primitive state aligned with x-axis of global problem
   LES3DFsd_pState RotateBack(const Vector3D &norm_dir) const;
//@}

/** @name Numerical Evaluation of Viscous Fluxes */
/*        -------------------------------------- */
//@{
   //! Returns viscous flux in n-direction
   static LES3DFsd_cState FluxViscous_n(const LES3DFsd_pState &Wl,
                                        const LES3DFsd_pState &Wr,
                                        const LES3DFsd_pState &W1,
                                        const LES3DFsd_pState &W2,
                                        const LES3DFsd_pState &dWdx1,
                                        const LES3DFsd_pState &dWdy1,
                                        const LES3DFsd_pState &dWdz1,
                                        const LES3DFsd_pState &dWdx2,
                                        const LES3DFsd_pState &dWdy2,
                                        const LES3DFsd_pState &dWdz2,
                                        const Vector3D &norm, 
                                        const Vector3D &ts, 
                                        const double &deltad, 
                                        const double &Volume, 
                                        const double &Volume_Neigbor, 
                                        const int Flow_Type);
//@}

/** @name Boundary Conditions */
/*        ------------------- */
//@{
   //! Return reflected solution state after application of reflection BC
   static LES3DFsd_pState Reflect(const LES3DFsd_pState &W,
                                  const Vector3D &norm_dir);

   //! Return wall solution state after application of moving wall BC
   static LES3DFsd_pState MovingWall(const LES3DFsd_pState &Win,
                                     const LES3DFsd_pState &Wout,
                                     const Vector3D &norm_dir,				 
                                     const Vector3D &wall_velocity,
                                     const Vector3D &pressure_gradient,
                                     const int &TEMPERATURE_BC_FLAG);

   //! Return wall solution state after application of no-slip BC
   static LES3DFsd_pState NoSlip(const LES3DFsd_pState &Win, 
                                 const LES3DFsd_pState &Wout, 
                                 const Vector3D &norm_dir,  
                                 const Vector3D &pressure_gradient,
                                 const int &TEMPERATURE_BC_FLAG);
//@}

/** @name Turbulence Model Source Terms */
/*        ----------------------------- */
//@{
   //! Absolute value of strain rate
   double abs_strain_rate(const LES3DFsd_pState &dWdx, 
                          const LES3DFsd_pState &dWdy, 
                          const LES3DFsd_pState &dWdz) const;

   //! Gradient of strain rate
   Vector3D grad_abs_strain_rate(const LES3DFsd_pState &dWdx, 
                                 const LES3DFsd_pState &dWdy, 
                                 const LES3DFsd_pState &dWdz,
                                 const LES3DFsd_pState &d_dWdx_dx,
                                 const LES3DFsd_pState &d_dWdy_dy,
                                 const LES3DFsd_pState &d_dWdz_dz,
                                 const LES3DFsd_pState &d_dWdx_dy,
                                 const LES3DFsd_pState &d_dWdx_dz,
		                 const LES3DFsd_pState &d_dWdy_dz) const;
   //! Heat release parameter
   double HeatRelease_Parameter(void) const;

   //! Subfilter scale kinetic energy
   double SFS_Kinetic_Energy_Fsd(const LES3DFsd_pState &dWdx,
                                 const LES3DFsd_pState &dWdy,
                                 const LES3DFsd_pState &dWdz,
                                 const int Flow_Type,
                                 const double &Volume);

   //! Efficiency function for subfilter scale strain term
   double Efficiency_Function_Fsd(const LES3DFsd_pState &dWdx,
                                  const LES3DFsd_pState &dWdy,
                                  const LES3DFsd_pState &dWdz,
                                  const int Flow_Type,
                                  const double &Volume);
 
   //! Gradients of progress variable to species mass fractions
   double Progvar_Species_Grad(void) const;

   //! Reaction rate
   double Reaction_Rate_Fsd(const LES3DFsd_pState &dWdx,
                            const LES3DFsd_pState &dWdy,
                            const LES3DFsd_pState &dWdz);

   //! x-direction surface averaged normal 
   double M_x(const LES3DFsd_pState &dWdx,
              const LES3DFsd_pState &dWdy,
              const LES3DFsd_pState &dWdz) const;
   //! y-direction surface averaged normal 
   double M_y(const LES3DFsd_pState &dWdx,
              const LES3DFsd_pState &dWdy,
              const LES3DFsd_pState &dWdz) const;
   //! z-direction surface averaged normal 
   double M_z(const LES3DFsd_pState &dWdx,
              const LES3DFsd_pState &dWdy,
              const LES3DFsd_pState &dWdz) const;

   //! Resolved strain term
   double Resolved_Strain(const LES3DFsd_pState &dWdx,
                          const LES3DFsd_pState &dWdy,
                          const LES3DFsd_pState &dWdz) const;

   //! Resolved propagation and curvature term
   double Resolved_Propagation_Curvature(const LES3DFsd_pState &dWdx,
                                         const LES3DFsd_pState &dWdy,
                                         const LES3DFsd_pState &dWdz) const;

   //! Subfilter scale strain term
   double SFS_Strain(const LES3DFsd_pState &dWdx,
                     const LES3DFsd_pState &dWdy,
                     const LES3DFsd_pState &dWdz,
                     const int Flow_Type,
                     const double &Volume);

   //! Subfilter scale curvature term
   double SFS_Curvature(const LES3DFsd_pState &dWdx,
                        const LES3DFsd_pState &dWdy,
                        const LES3DFsd_pState &dWdz) const;

   //! Gradient of surface averged normal
   double M_xx(const LES3DFsd_pState &dWdx,
               const LES3DFsd_pState &dWdy,
               const LES3DFsd_pState &dWdz,
               const LES3DFsd_pState &d_dWdx_dx,
               const LES3DFsd_pState &d_dWdx_dy,
               const LES3DFsd_pState &d_dWdx_dz) const;
   //! Gradient of surface averged normal
   double M_yy(const LES3DFsd_pState &dWdx,
               const LES3DFsd_pState &dWdy,
               const LES3DFsd_pState &dWdz,
               const LES3DFsd_pState &d_dWdy_dy,
               const LES3DFsd_pState &d_dWdx_dy,
               const LES3DFsd_pState &d_dWdy_dz) const;
   //! Gradient of surface averged normal
   double M_zz(const LES3DFsd_pState &dWdx,
               const LES3DFsd_pState &dWdy,
               const LES3DFsd_pState &dWdz,
               const LES3DFsd_pState &d_dWdz_dz,
               const LES3DFsd_pState &d_dWdx_dz,
               const LES3DFsd_pState &d_dWdy_dz) const;
   //! Gradient of surface averged normal
   double M_xy(const LES3DFsd_pState &dWdx,
               const LES3DFsd_pState &dWdy,
               const LES3DFsd_pState &dWdz,
               const LES3DFsd_pState &d_dWdy_dy,
               const LES3DFsd_pState &d_dWdx_dy,
               const LES3DFsd_pState &d_dWdy_dz) const;
   //! Gradient of surface averged normal
   double M_xz(const LES3DFsd_pState &dWdx,
               const LES3DFsd_pState &dWdy,
               const LES3DFsd_pState &dWdz,
               const LES3DFsd_pState &d_dWdz_dz,
               const LES3DFsd_pState &d_dWdx_dz,
               const LES3DFsd_pState &d_dWdy_dz) const;
   //! Gradient of surface averged normal
   double M_yz(const LES3DFsd_pState &dWdx,
               const LES3DFsd_pState &dWdy,
               const LES3DFsd_pState &dWdz,
               const LES3DFsd_pState &d_dWdz_dz,
               const LES3DFsd_pState &d_dWdx_dz,
               const LES3DFsd_pState &d_dWdy_dz) const;

   //! Resolved curvature term
   double Resolved_Curvature(const LES3DFsd_pState &dWdx,
                             const LES3DFsd_pState &dWdy,
                             const LES3DFsd_pState &dWdz,
                             const LES3DFsd_pState &d_dWdx_dx,
                             const LES3DFsd_pState &d_dWdy_dy,
                             const LES3DFsd_pState &d_dWdz_dz,
                             const LES3DFsd_pState &d_dWdx_dy,
                             const LES3DFsd_pState &d_dWdx_dz,
                             const LES3DFsd_pState &d_dWdy_dz) const;

   //! Resolved propagation term
   double Resolved_Propagation(const LES3DFsd_pState &dWdx,
                               const LES3DFsd_pState &dWdy,
                               const LES3DFsd_pState &dWdz,
                               const LES3DFsd_pState &d_dWdx_dx,
                               const LES3DFsd_pState &d_dWdy_dy,
                               const LES3DFsd_pState &d_dWdz_dz,
                               const LES3DFsd_pState &d_dWdx_dy,
                               const LES3DFsd_pState &d_dWdx_dz,
                               const LES3DFsd_pState &d_dWdy_dz) const;

   //! Resolved convection term for progress variable
   double Resolved_Convection_Progvar(const LES3DFsd_pState &dWdx,
                                      const LES3DFsd_pState &dWdy,
                                      const LES3DFsd_pState &dWdz) const;

   //! Resolved convection term for flame surface density
   double Resolved_Convection_Fsd(const LES3DFsd_pState &dWdx,
                                  const LES3DFsd_pState &dWdy,
                                  const LES3DFsd_pState &dWdz) const;

   //! Non-gradient term for progress variable
   double NGT_Progvar(const LES3DFsd_pState &dWdx,
                      const LES3DFsd_pState &dWdy,
                      const LES3DFsd_pState &dWdz) const;

   //! Non-gradient term for flame surface density
   double NGT_Fsd(const LES3DFsd_pState &dWdx,
                  const LES3DFsd_pState &dWdy,
                  const LES3DFsd_pState &dWdz,
                  const LES3DFsd_pState &d_dWdx_dx,
                  const LES3DFsd_pState &d_dWdy_dy,
                  const LES3DFsd_pState &d_dWdz_dz,
                  const LES3DFsd_pState &d_dWdx_dy,
                  const LES3DFsd_pState &d_dWdx_dz,
                  const LES3DFsd_pState &d_dWdy_dz) const;

   //! Subfilter scale term for progress variable
   double SFS_Diffusion_Progvar(const LES3DFsd_pState &dWdx,
                                const LES3DFsd_pState &dWdy,
                                const LES3DFsd_pState &dWdz,
                                const LES3DFsd_pState &d_dWdx_dx,
                                const LES3DFsd_pState &d_dWdy_dy,
                                const LES3DFsd_pState &d_dWdz_dz,
                                const LES3DFsd_pState &d_dWdx_dy,
                                const LES3DFsd_pState &d_dWdx_dz,
                                const LES3DFsd_pState &d_dWdy_dz,
                                const int Flow_Type,
                                const double &Volume);

   //! Subfilter scale term for flame surface density
   double SFS_Diffusion_Fsd(const LES3DFsd_pState &dWdx,
                            const LES3DFsd_pState &dWdy,
                            const LES3DFsd_pState &dWdz,
                            const LES3DFsd_pState &d_dWdx_dx,
                            const LES3DFsd_pState &d_dWdy_dy,
                            const LES3DFsd_pState &d_dWdz_dz,
                            const LES3DFsd_pState &d_dWdx_dy,
                            const LES3DFsd_pState &d_dWdx_dz,
             	            const LES3DFsd_pState &d_dWdy_dz,
                            const int Flow_Type,
                            const double &Volume); 

   //! Heat release strain term
   double Heat_Release_Strain(const LES3DFsd_pState &dWdx,
                              const LES3DFsd_pState &dWdy,
                              const LES3DFsd_pState &dWdz,
                              const LES3DFsd_pState &d_dWdx_dx,
                              const LES3DFsd_pState &d_dWdy_dy,
                              const LES3DFsd_pState &d_dWdz_dz,
                              const LES3DFsd_pState &d_dWdx_dy,
                              const LES3DFsd_pState &d_dWdx_dz,
                              const LES3DFsd_pState &d_dWdy_dz);


   //! Net rate change for progress variable
   double Net_Rate_Change_Progvar(const LES3DFsd_pState &dWdx,
                                  const LES3DFsd_pState &dWdy,
                                  const LES3DFsd_pState &dWdz,
                                  const LES3DFsd_pState &d_dWdx_dx,
                                  const LES3DFsd_pState &d_dWdy_dy,
                                  const LES3DFsd_pState &d_dWdz_dz,
                                  const LES3DFsd_pState &d_dWdx_dy,
                                  const LES3DFsd_pState &d_dWdx_dz,
		                  const LES3DFsd_pState &d_dWdy_dz,
                                  const int Flow_Type,
                                  const double &Volume);

   //! Net rate change for flame surface density
   double Net_Rate_Change_Fsd(const LES3DFsd_pState &dWdx,
                              const LES3DFsd_pState &dWdy,
                              const LES3DFsd_pState &dWdz,
                              const LES3DFsd_pState &d_dWdx_dx,
                              const LES3DFsd_pState &d_dWdy_dy,
                              const LES3DFsd_pState &d_dWdz_dz,
                              const LES3DFsd_pState &d_dWdx_dy,
                              const LES3DFsd_pState &d_dWdx_dz,
		              const LES3DFsd_pState &d_dWdy_dz,
                              const int Flow_Type,
                              const double &Volume);

   //! Soure term for k-equation
   double K_equ_sources(const LES3DFsd_pState &dWdx,
                        const LES3DFsd_pState &dWdy,
                        const LES3DFsd_pState &dWdz,
                        const int Flow_Type,
                        const double &Volume);

   //! Source terms for C-Fsd-k model
   static LES3DFsd_cState Sturbulence(LES3DFsd_pState &Wc,
                                      const LES3DFsd_pState &dWdx,
                                      const LES3DFsd_pState &dWdy,
                                      const LES3DFsd_pState &dWdz,
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
   LES3DFsd_pState operator +(const LES3DFsd_pState &W) const;

   //! Binary subtraction operator
   LES3DFsd_pState operator -(const LES3DFsd_pState &W) const;

   //! Binary multiplication operator
   LES3DFsd_pState operator *(const double &a) const;

   //! Binary multiplication operator
   friend LES3DFsd_pState operator *(const double &a, const LES3DFsd_pState &W);

   //! Binary multiplication operator
   double operator *(const LES3DFsd_pState &W) const;

   //! Binary division operator
   LES3DFsd_pState operator /(const double &a) const;

   //! Binary vector product operator
   LES3DFsd_pState operator ^(const LES3DFsd_pState &W) const;
   
   //! Assignment operator
   LES3DFsd_pState& operator =(const LES3DFsd_pState &W);

   //! Shortcut addition operator
   LES3DFsd_pState& operator +=(const LES3DFsd_pState &W);

   //! Shortcut subtraction operator
   LES3DFsd_pState& operator -=(const LES3DFsd_pState &W);
  
   //! Equal relational operator
   friend int operator ==(const LES3DFsd_pState &W1,
                          const LES3DFsd_pState &W2);

   //! Not equal relational operator
   friend int operator !=(const LES3DFsd_pState &W1,
                          const LES3DFsd_pState &W2);
  
   //! Output stream operator
   friend ostream& operator << (ostream &out_file,
                                const LES3DFsd_pState &W);

   //! Input stream operator
   friend istream& operator >> (istream &in_file,
                                LES3DFsd_pState &W);
//@}
};
  
/*!
 * Class: LES3DFsd_cState
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
class LES3DFsd_cState : public NavierStokes3D_ThermallyPerfect_cState {
  protected:
   static double        _fuel_equivalence_ratio;  //!< Fuel equivalence ratio
   static double    _unburnt_fuel_mass_fraction;  //!< Mass fraction of unburnt fuel (premixed flame)
   static double             _reactants_density;  //!< Reactants Density
   static double           _laminar_flame_speed;  //!< Propagation speed of laminar premixed flame
   static double       _laminar_flame_thickness;  //!< Thickness of laminar premixed flame
   static double   _adiabatic_flame_temperature;  //!< Adiabitc flame temperature

  public:
   double                  rhoC; //!< Density of reaction rate progress variable (0 <= C <= 1)
   double                rhoFsd; //!< Flame surface density (FSD) of premixed flame
   double                  rhok; //!< Subfilter scale turbulent kinetic energy (kg/m-s^2)

   static double           Mref; //!< Reference Mach number for low-Mach-number precondtioning (normally set to incoming freestream Mach)

/** @name Constructors and destructors */
/*        ---------------------------- */
//@{
   //! Default creation constructor (assign default values)
   LES3DFsd_cState(): NavierStokes3D_ThermallyPerfect_cState(), rhoC(rho*ONE), rhoFsd(MILLION), rhok(ZERO) {
     premixed_mfrac(); }
   
   //! Constructor from base class (allows return of derived type)
   LES3DFsd_cState(const NavierStokes3D_ThermallyPerfect_cState &U1) : 
     NavierStokes3D_ThermallyPerfect_cState(U1), rhoC(rho*ONE), rhoFsd(MILLION), rhok(ZERO) { }

   //! Constructor from base class
   LES3DFsd_cState(const NavierStokes3D_ThermallyPerfect_cState &U1,
                   const double &cc,
                   const double &ff, 
                   const double &kk) : 
     NavierStokes3D_ThermallyPerfect_cState(U1), rhoC(cc), rhoFsd(ff), rhok(kk) {
     premixed_mfrac(); }

   //! Assignment constructor
   LES3DFsd_cState(const double &value): 
     NavierStokes3D_ThermallyPerfect_cState(value), rhoC(value), rhoFsd(value), rhok(value) {
     premixed_mfrac(); }
   
   //! Assignment constructor
   LES3DFsd_cState(const double &d, 
                   const Vector3D &V,
                   const double &En, 
                   const double &cc, 
                   const double &ff, 
                   const double &kk):
     NavierStokes3D_ThermallyPerfect_cState(d, V, En), rhoC(cc), rhoFsd(ff), rhok(kk) { 
     premixed_mfrac(); }

   //! Assignment constructor
   LES3DFsd_cState(const double &d, 
                   const double &vx,
                   const double &vy, 
                   const double &vz,
                   const double &En, 
                   const double &cc, 
                   const double &ff, 
                   const double &kk):
     NavierStokes3D_ThermallyPerfect_cState(d, vx, vy, vz, En), rhoC(cc), rhoFsd(ff), rhok(kk) {
     premixed_mfrac(); }
   
   //! Copy constructor (this is needed for the operator overload returns)
   LES3DFsd_cState(const LES3DFsd_cState &U) {
     rhospec_null(); rho = DENSITY_STDATM; set_initial_values(); Copy(U);
   }
                           
   //! Default destructor
   ~LES3DFsd_cState(void) {
      Deallocate();
   }
//@}
                           
/** @name Some useful operators */
/*        --------------------- */
//@{
   void Copy(const LES3DFsd_cState &U);

   void Vacuum(void) { 
     rho = ZERO; rhov.zero(); E=ZERO; rhoC = ZERO; rhoFsd = ZERO; rhok = ZERO;
   }

   //! Check for physical validity of the progress variable
   void Realizable_C_Check(void) {
      if (rhoC/rho > ONE) {
        rhoC = rho*ONE;
      } else if (rhoC/rho < ZERO) {
        rhoC = ZERO;
      } /* endif */
   }

   //! Check for physical validity of the progress variable
   void Realizable_Mass_Fraction_Check(void) {
      if (rhoC/rho > ONE) {
        rhoC = ONE*rho;
        premixed_mfrac();
      } else if (rhoC < ZERO) {
        rhoC = ZERO;
        premixed_mfrac();
      } /* endif */
   }

   //! Check for physical validity of scalars
   bool Realizable_Scalar_Check(void) {
     double LOCAL_TOL = MICRO;
     if ( rhoC < LOCAL_TOL ) { rhoC = ZERO; premixed_mfrac(); }
     if ( rhoFsd < LOCAL_TOL ) { rhoFsd = ZERO; }
     if ( rhok < LOCAL_TOL ) { rhok = ZERO; }
     if ( rhoC/rho > 0.999 || rhoC/rho < 0.001 ) { rhoFsd = ZERO; }
     return true;
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

   //! Mixture sound speed including turbulent kinetic energy
   double a_t(void) const;

   //! Progress variable
   double C(void) const;

   //! Flame surface density
   double Fsd(void) const;

   //! Turbulent kinetic energy
   double k(void) const;

   //! Species mass fractions
   void premixed_mfrac(void);
//@}

/** @name Turbulent transport coefficients */
/*        -------------------------------- */
//@{
   //! Eddy (turbulent) viscosity
   double mu_t(const LES3DFsd_pState &dWdx, 
               const LES3DFsd_pState &dWdy, 
               const LES3DFsd_pState &dWdz,
               const int Flow_Type, 
               const double &Volume);

   //! Turbulent) thermal conductivity
   double kappa_t(const LES3DFsd_pState &dWdx, 
                  const LES3DFsd_pState &dWdy, 
                  const LES3DFsd_pState &dWdz,
                  const int Flow_Type, 
                  const double &Volume);     

   //! Species turbulent diffusion coefficient
   double Ds_t(const int i,
               const LES3DFsd_pState &dWdx, 
               const LES3DFsd_pState &dWdy, 
               const LES3DFsd_pState &dWdz,
               const int Flow_Type, 
               const double &Volume);
   //! Species turbulent diffusion coefficient
   double Ds_t(const int i,
               const double &mu_t_temp);

   //! Turbulent Prandtl number
   double Pr_t(void);

   //! Turbulent Schmidt number
   double Sc_t(void);

   //! Turbulent Lewis number
   double Le_t(void);
 
   //! LES filter width
   double filter_width(const double &Volume) const;

   //! Absolute value of strain rate
   double abs_strain_rate(const LES3DFsd_pState &dWdx, 
                          const LES3DFsd_pState &dWdy, 
                          const LES3DFsd_pState &dWdz) const;
//@}

/** @name Primitive solution state */ 
/*        ------------------------ */
//@{
   //! Returns primitive solution state
   LES3DFsd_pState W(void) ; 
   //! Returns primitive solution state
   LES3DFsd_pState W(void) const;
   //! Returns primitive solution state
   LES3DFsd_pState W(const LES3DFsd_cState &U) const;
//@}

/** @name Numerical Flux Functions */
/*        ------------------------ */
//@{
   //! Returns rotated conserved state aligned with x-axis in the norm_dir
   LES3DFsd_cState Rotate(const Vector3D &norm_dir) const;

   //! Returns un-rotated conserved state aligned with x-axis of the global problem
   LES3DFsd_cState RotateBack(const Vector3D &norm_dir) const;
//@}

/** @name Operators */
/*        --------- */
//@{
   //! Index operator
   double &operator[](int index);
   //! Index operator
   const double &operator[](int index) const;

   //! Binary addition operator
   LES3DFsd_cState operator +(const LES3DFsd_cState &U) const;

   //! Binary subtraction operator
   LES3DFsd_cState operator -(const LES3DFsd_cState &U) const;

   //! Binary multiplication operator
   LES3DFsd_cState operator *(const double &a) const;

   //! Binary multiplication operator
   friend LES3DFsd_cState operator *(const double &a, const LES3DFsd_cState &U);

   //! Binary multiplication operator
   double operator *(const LES3DFsd_cState &U) const;

   //! Binary division operator
   LES3DFsd_cState operator /(const double &a) const;

   //! Binary vector product operator
   LES3DFsd_cState operator ^(const LES3DFsd_cState &U) const;
   
   //! Assignment operator
   LES3DFsd_cState& operator =(const LES3DFsd_cState &U);

   //! Shortcut addition operator
   LES3DFsd_cState& operator +=(const LES3DFsd_cState &U);

   //! Shortcut subtraction operator
   LES3DFsd_cState& operator -=(const LES3DFsd_cState &U);
  
   //! Equal relational operator
   friend int operator ==(const LES3DFsd_cState &U1,
                          const LES3DFsd_cState &U2);

   //! Not equal relational operator
   friend int operator !=(const LES3DFsd_cState &U1,
                          const LES3DFsd_cState &U2);
  
   //! Output stream operator
   friend ostream& operator << (ostream &out_file,
                                const LES3DFsd_cState &U);

   //! Input stream operator
   friend istream& operator >> (istream &in_file,
                                LES3DFsd_cState &U);
//@}
};

/***************************************
 * LES3DFsd_pState member functions    *
 ***************************************/

/***************************************
 * LES3DFsd_pState -- Index operators. *
 ***************************************/
inline double& LES3DFsd_pState::operator[](int index) {
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
      return C;
   case 7:
      return Fsd;
   case 8:
      return k;
   default :
      return spec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3D_VAR_EXTRA +1)].c;
   };
}

inline const double& LES3DFsd_pState::operator[](int index) const {
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
      return C;
   case 7:
      return Fsd;
   case 8:
      return k;
   default :
      return spec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3D_VAR_EXTRA +1)].c;
   };
}

/***************************************************
 * LES3DFsd_pState -- Binary arithmetic operators. *
 ***************************************************/
//------------------ Addition ------------------------//
inline LES3DFsd_pState LES3DFsd_pState::operator +(const LES3DFsd_pState &W) const { 
   LES3DFsd_pState Temp;
   Temp.Copy(*this);
   Temp += W;
   return Temp;
}

//------------------ Subtraction ------------------------//
inline LES3DFsd_pState LES3DFsd_pState::operator -(const LES3DFsd_pState &W) const {
   LES3DFsd_pState Temp;
   Temp.Copy(*this);
   Temp -= W;
   return Temp;
}

//---------------- Scalar Multiplication ------------------//
inline LES3DFsd_pState LES3DFsd_pState::operator *(const double &a) const {
   LES3DFsd_pState Temp;
   Temp.Copy(*this);
   Temp.rho = rho*a;  
   Temp.v = v*a; 
   Temp.p = p*a;
   Temp.C = C*a; 
   Temp.Fsd = Fsd*a; 
   Temp.k = k*a; 
/*    for( int i=0; i<ns; i++){ */
/*       Temp.spec[i] = spec[i]*a; */
/*    } /\* endfor *\/  */
   return (Temp);
}

inline LES3DFsd_pState operator *(const double &a, 
                                  const LES3DFsd_pState &W) {
   LES3DFsd_pState Temp;
   Temp.rho = W.rho*a;  
   Temp.v = W.v*a; 
   Temp.p = W.p*a;
   Temp.C = W.C*a; 
   Temp.Fsd = W.Fsd*a; 
   Temp.k = W.k*a; 
/*    for( int i=0; i<W.ns; i++){ */
/*       Temp.spec[i] = W.spec[i]*a; */
/*    } /\* endfor *\/  */
   return (Temp);
}

//--------------- Scalar Division ------------------------//
inline LES3DFsd_pState LES3DFsd_pState::operator /(const double &a) const {
   LES3DFsd_pState Temp(rho,v,p,C,Fsd,k);
   Temp.Copy(*this);
   Temp.rho = rho/a; 
   Temp.v = v/a; 
   Temp.p = p/a; 
   Temp.C = C/a; 
   Temp.Fsd = Fsd/a; 
   Temp.k = k/a; 
/*    for(int i=0; i<ns; i++){ */
/*       Temp.spec[i] = spec[i]/a;  */
/*    } /\* endfor *\/ */ 
   return(Temp);
}

//----------------- Inner Product ------------------------//
inline double LES3DFsd_pState::operator *(const LES3DFsd_pState &W) const {
   double sum=0.0;
/*    for(int i=0; i<ns; i++){ */
/*       sum += spec[i]*W.spec[i]; */
/*    } */  
   return (rho*W.rho + v*W.v + p*W.p + C*W.C + Fsd*W.Fsd + k*W.k + sum);
}

//----------- Solution state product operator ------------//
inline LES3DFsd_pState LES3DFsd_pState::operator ^(const LES3DFsd_pState &W) const {
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
/*    for(int i=0; i<ns; i++){ */
/*       Temp.spec[i] = spec[i]*W.spec[i]; */
/*    } /\* endfor *\/  */ 
   return(Temp);
}

//----------------- Assignment ----------------------------//
inline LES3DFsd_pState& LES3DFsd_pState::operator =(const LES3DFsd_pState &W) {
   //self assignment protection
   if( this != &W){   
      //copy assignment
      rho = W.rho;
      v = W.v; 
      p = W.p; 
      C = W.C;
      Fsd = W.Fsd;
      k = W.k;
/*       for(int i=0; i<ns; i++){ */
/*          spec[i] = W.spec[i]; */
/*        } /\* endfor *\/ */   
   } /* endif */   
   return (*this);
}

/******************************************************
 * LES3DFsd_pState -- Shortcut arithmetic operators.  *         
 ******************************************************/
inline LES3DFsd_pState& LES3DFsd_pState::operator +=(const LES3DFsd_pState &W) {
   rho += W.rho;
   v += W.v; 
   p += W.p; 
   C += W.C;
   Fsd += W.Fsd;
   k += W.k;
/*    for (int i=0; i<ns; i++) { */
/*       spec[i] += W.spec[i]; */
/*    } /\* endfor *\/ */ 
   return (*this);
}

inline LES3DFsd_pState& LES3DFsd_pState::operator -=(const LES3DFsd_pState &W) {
   rho -= W.rho;
   v -= W.v;
   p -= W.p;
   C -= W.C;
   Fsd -= W.Fsd;
   k -= W.k;
/*    for (int i=0; i<ns; i++) { */
/*       spec[i] -= W.spec[i]; */
/*    } /\* endfor *\/ */  
   return (*this); 
}

/*****************************************************
 * LES3DFsd_pState -- Unary arithmetic operators.    *
 *****************************************************/
inline LES3DFsd_pState operator -(const LES3DFsd_pState &W) {
   return (LES3DFsd_pState(-W.rho, -W.v, -W.p, -W.C, -W.Fsd, -W.k));
}

/********************************************
 * LES3DFsd_pState -- Relational operators. *         
 ********************************************/
inline int operator ==(const LES3DFsd_pState &W1, 
                       const LES3DFsd_pState &W2) {
      return (W1.rho == W2.rho && 
              W1.v == W2.v && 
              W1.p == W2.p &&
              W1.C == W2.C && 
              W1.Fsd == W2.Fsd && 
              W1.k == W2.k);
}

inline int operator !=(const LES3DFsd_pState &W1, 
                       const LES3DFsd_pState &W2) {
   return (W1.rho != W2.rho || 
           W1.v != W2.v || 
           W1.p != W2.p ||
           W1.C != W2.C || 
           W1.Fsd != W2.Fsd || 
           W1.k != W2.k);
}

/*************************************************** 
 *  LES3DFsd_pState -- Input-output operators.     *
 ***************************************************/
inline ostream &operator << (ostream &out_file, 
                             const LES3DFsd_pState &W) {
   out_file.precision(10);
   out_file.setf(ios::scientific);
   out_file <<" "<< W.rho << " " << W.v.x << " " << W.v.y << " " << W.v.z
            <<" "<< W.p << " " << W.C << " " << W.Fsd << " " << W.k;
   for (int i=0; i < W.ns; i++) {
      out_file<<" "<<W.spec[i];
   } /* endfor */
   out_file.unsetf(ios::scientific);
   return (out_file);
}

inline istream &operator >> (istream &in_file, 
                             LES3DFsd_pState &W) {
   in_file.setf(ios::skipws);
   in_file >> W.rho >> W.v.x >> W.v.y >> W.v.z >> W.p >> W.C >> W.Fsd >> W.k;
   for (int i=0; i<W.ns; i++) {
      in_file>>W.spec[i];
   } /* endfor */
   in_file.unsetf(ios::skipws);
   return (in_file);
}

/***************************************
 * LES3DFsd_cState member functions    *
 ***************************************/

/***************************************
 * LES3DFsd_pState -- Index operators. *
 ***************************************/
inline double& LES3DFsd_cState::operator[](int index) {
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
      return rhoC;
   case 7:
      return rhoFsd;
   case 8:
      return rhok;
   default :
      return rhospec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3D_VAR_EXTRA +1)].c;
   };
}

inline const double& LES3DFsd_cState::operator[](int index) const {
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
      return rhoC;
   case 7:
      return rhoFsd;
   case 8:
      return rhok;
   default :
      return rhospec[index-(NUM_EULER3D_VAR_SANS_SPECIES + NUM_LES3D_VAR_EXTRA +1)].c;
   };
}

//----------------- Addition ----------------------------//
inline LES3DFsd_cState LES3DFsd_cState::operator +(const LES3DFsd_cState &U) const { 
   LES3DFsd_cState Temp;
   Temp.Copy(*this);
   Temp += U;
   return Temp;
}

//------------------ Subtraction ------------------------//
inline LES3DFsd_cState LES3DFsd_cState::operator -(const LES3DFsd_cState &U) const {
   LES3DFsd_cState Temp;
   Temp.Copy(*this);
   Temp -= U;
   return Temp;
}

//---------------- Scalar Multiplication ------------------//
inline LES3DFsd_cState LES3DFsd_cState::operator *(const double &a) const {
   LES3DFsd_cState Temp;
   Temp.Copy(*this);
   Temp.rho = rho*a; 
   Temp.rhov = rhov*a; 
   Temp.E = E*a;
   Temp.rhoC = rhoC*a; 
   Temp.rhoFsd = rhoFsd*a; 
   Temp.rhok = rhok*a;
/*    for( int i=0; i<ns; i++){ */
/*       Temp.rhospec[i] = rhospec[i]*a; */
/*    } /\* endfor *\/  */
   return (Temp);
}

inline LES3DFsd_cState operator *(const double &a, const LES3DFsd_cState &U) {
   LES3DFsd_cState Temp;
   Temp.rho = U.rho*a;  
   Temp.rhov = U.rhov*a; 
   Temp.E = U.E*a;
   Temp.rhoC = U.rhoC*a; 
   Temp.rhoFsd = U.rhoFsd*a; 
   Temp.rhok = U.rhok*a;
/*    for( int i=0; i<U.ns; i++){ */
/*       Temp.rhospec[i] = U.rhospec[i]*a; */
/*    } /\* endfor *\/  */
   return (Temp);
}

//--------------- Scalar Division ------------------------//
inline LES3DFsd_cState LES3DFsd_cState::operator /(const double &a) const {
   LES3DFsd_cState Temp;
   Temp.Copy(*this);
   Temp.rho = rho/a; 
   Temp.rhov = rhov/a; 
   Temp.E = E/a;
   Temp.rhoC = rhoC/a; 
   Temp.rhoFsd = rhoFsd/a; 
   Temp.rhok = rhok/a;
/*    for(int i=0; i<ns; i++){ */
/*       Temp.rhospec[i] = rhospec[i]/a;  */
/*    } /\* endfor *\/ */ 
   return(Temp);
}

//----------------- Inner Product ------------------------//
inline double LES3DFsd_cState::operator *(const LES3DFsd_cState &U) const{
   double sum=0.0;
/*    for(int i=0; i<ns; i++){ */
/*       sum += rhospec[i]*U.rhospec[i]; */
/*    } /\* endfor *\/ */  
   return (rho*U.rho+rhov*U.rhov+E*U.E+rhoC*U.rhoC+rhoFsd*U.rhoFsd+rhok*U.rhok+sum);
}

//----------- Solution state product operator ------------//
inline LES3DFsd_cState LES3DFsd_cState::operator ^(const LES3DFsd_cState &U) const {
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
/*    for(int i=0; i<ns; i++){ */
/*       Temp.rhospec[i] = rhospec[i]*U.rhospec[i]; */
/*    } /\* endfor *\/ */  
   return(Temp);
}

//----------------- Assignment ----------------------------//
inline LES3DFsd_cState& LES3DFsd_cState::operator =(const LES3DFsd_cState &U) {
   //self assignment protection
   if( this != &U){   
      //copy assignment
      rho = U.rho;
      rhov = U.rhov; 
      E = U.E; 
      rhoC = U.rhoC;
      rhoFsd = U.rhoFsd;
      rhok = U.rhok;
/*       for(int i=0; i<ns; i++){ */
/*           rhospec[i] = U.rhospec[i]; */
/*       } /\* endfor * */
   } /* endif */
   return (*this);
}

/*******************************************************************
 *        LES3DFsd_cState -- Shortcut arithmetic operators.        *  
 *******************************************************************/
inline LES3DFsd_cState& LES3DFsd_cState::operator +=(const LES3DFsd_cState &U){
   rho += U.rho;
   rhov += U.rhov; 
   E += U.E;
   rhoC += U.rhoC;
   rhoFsd += U.rhoFsd;
   rhok += U.rhok;
/*    for(int i=0; i<ns; i++) { */
/*       rhospec[i] += U.rhospec[i]; */
/*    } /\* endfor *\/ */
   return (*this);
}

inline LES3DFsd_cState& LES3DFsd_cState::operator -=(const LES3DFsd_cState &U) {
   rho -= U.rho;
   rhov -= U.rhov;
   E -= U.E;
   rhoC -= U.rhoC;
   rhoFsd -= U.rhoFsd;
   rhok -= U.rhok;
/*    for (int i=0; i<ns; i++) { */
/*      rhospec[i] -= U.rhospec[i]; */
/*    } /\* endfor *\/  */
   return (*this); 
}

/***************************************************
 * LES3DFsd_cState -- Unary arithmetic operators.  *          
 ***************************************************/
inline LES3DFsd_cState operator -(const LES3DFsd_cState &U) {
   return (LES3DFsd_cState(-U.rho,-U.rhov,-U.E,-U.rhoC,-U.rhoFsd,-U.rhok));
}

/*********************************************
 * LES3DFsd_cState -- Relational operators.  *          
 *********************************************/
inline int operator ==(const LES3DFsd_cState &U1, 
                       const LES3DFsd_cState &U2) {
   return (U1.rho == U2.rho && 
           U1.rhov == U2.rhov && 
           U1.E == U2.E &&
           U1.rhoC == U2.rhoC && 
           U1.rhoFsd == U2.rhoFsd && 
           U1.rhok == U2.rhok);
}

inline int operator !=(const LES3DFsd_cState &U1, 
                       const LES3DFsd_cState &U2) {
   return (U1.rho != U2.rho || 
           U1.rhov != U2.rhov || 
           U1.E != U2.E ||
           U1.rhoC != U2.rhoC || 
           U1.rhoFsd !=U2.rhoFsd || 
           U1.rhok != U2.rhok);
}

/*************************************************
 *  LES3DFsd_cState -- Input-output operators.   *      
 *************************************************/
inline ostream &operator << (ostream &out_file, const LES3DFsd_cState &U) {
  out_file.setf(ios::scientific);
  out_file << " " << U.rho << " " << U.rhov.x << " " << U.rhov.y << " " << U.rhov.z
           << " " << U.E << " " << U.rhoC << " " << U.rhoFsd << " " << U.rhok;
  for (int i=0; i<U.ns; i++) {
     out_file<<" "<<U.rhospec[i];
  } /* endfor */
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, LES3DFsd_cState &U) {
   in_file.setf(ios::skipws);
   in_file >> U.rho >> U.rhov.x >> U.rhov.y >> U.rhov.z >> U.E >> U.rhoC >> U.rhoFsd >> U.rhok;
  for (int i=0; i<U.ns; i++) {
    in_file>>U.rhospec[i]; 
  } /* endfor */
  in_file.unsetf(ios::skipws);
  return (in_file);
}

#endif // _LES3DFSD_STATE_INCLUDED
