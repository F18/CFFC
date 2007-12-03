/*! \file NavierStokes3DThermallyPerfectState.h
 * 	\brief	Header file defining the Navier-Stokes solution state classes 
 *              associated with solution of compressible viscous flows 
 *              of a thermally perfect non-reactive or combusting mixture.
 */

#ifndef _NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED 
#define _NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED

/* Include header file for base solution classes from which the classes are derived. */

#ifndef _EULER3D_THERMALLYPERFECT_STATE_INCLUDED
#include "../Euler/Euler3DThermallyPerfectState.h"
#endif  //EULER3D_THERMALLYPERFECT_STATE_INCLUDED

/* Define the classes. */

class NavierStokes3D_ThermallyPerfect_pState;
class NavierStokes3D_ThermallyPerfect_cState;

/*!
 * Class: NavierStokes3D_ThermallyPerfect_pState
 *
 * \brief Primitive state solution class for 3D Navier-Stokes equations
 *        governing flows of thermally perfect non-reactive and 
 *        combusting mixtures.
 *
 * Member functions
 *  - rho           -- Return density (kg/m^3)
 *  - v             -- Return flow velocity (m/s)
 *  - p             -- Return pressure (Pa, N/m^2)
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
 *  - U             -- Return conserved solution state
 *  - F, Fx         -- Return x-direction inviscid solution flux
 *  - Fy            -- Return y-direction inviscid solution flux
 *  - Fz            -- Return z-direction inviscid solution flux
 *  - Fv, Fvx       -- Return x-direction viscous solution flux
 *  - Fvy           -- Return y-direction viscous solution flux
 *  - Fvz           -- Return z-direction viscous solution flux
 *  - Schemistry    -- Return source terms associated with finite-rate chemistry
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
 *  - lambda_minus  -- Return negative eigenvalues, applying Harten entropy fix
 *  - lambda_plus   -- Return positive eigenvalues, applying Harten entropy fix
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
class NavierStokes3D_ThermallyPerfect_pState : public Euler3D_ThermallyPerfect_pState {
  public:

/** @name Constructors and destructors */
/*        ---------------------------- */
//@{
   //! Default creation constructor (assign default values)
   NavierStokes3D_ThermallyPerfect_pState(void) : Euler3D_ThermallyPerfect_pState() { }

   //! Constructor from base class (allows return of derived type)
   NavierStokes3D_ThermallyPerfect_pState(const Euler3D_ThermallyPerfect_pState &W1): 
     Euler3D_ThermallyPerfect_pState(W1) { }
   
   //! Assignment constructor
   NavierStokes3D_ThermallyPerfect_pState(const double &value) :
     Euler3D_ThermallyPerfect_pState(value) { }
   
   //! Assignment constructor
   NavierStokes3D_ThermallyPerfect_pState(const double &d, 
                                          const Vector3D &V,
                                          const double &pre) :
     Euler3D_ThermallyPerfect_pState(d, V, pre) { }
   
   //! Assignment constructor
   NavierStokes3D_ThermallyPerfect_pState(const double &d, 
                                          const Vector3D &V,
                                          const double &pre, 
                                          const double &frac) :
     Euler3D_ThermallyPerfect_pState(d, V, pre, frac) { }
   
   //! Assignment constructor
   NavierStokes3D_ThermallyPerfect_pState(const double &d, 
                                          const Vector3D &V,
                                          const double &pre, 
                                          const Species *mfrac) :
     Euler3D_ThermallyPerfect_pState(d, V, pre, mfrac) { }
   
   //! Assignment constructor
   NavierStokes3D_ThermallyPerfect_pState(const double &d, 
                                          const Vector3D &V,
                                          const double &pre, 
                                          const double *mfrac):
     Euler3D_ThermallyPerfect_pState(d, V, pre, mfrac) { }

   //! Assignment constructor
   NavierStokes3D_ThermallyPerfect_pState(const double &d, 
                                          const double &vx,
                                          const double &vy, 
                                          const double &vz,
                                          const double &pre) :
     Euler3D_ThermallyPerfect_pState(d, vx, vy, vz, pre) { }

   //! Assignment constructor
   NavierStokes3D_ThermallyPerfect_pState(const double &d, 
                                          const double &vx,
                                          const double &vy, 
                                          const double &vz,
                                          const double &pre, 
                                          const double &frac) :
     Euler3D_ThermallyPerfect_pState(d, vx, vy, vz, pre, frac) { }
   
   //! Assignment constructor
   NavierStokes3D_ThermallyPerfect_pState(const double &d, 
                                          const double &vx,
                                          const double &vy,  
                                          const double &vz,
                                          const double &pre, 
                                          const Species *mfrac):
     Euler3D_ThermallyPerfect_pState(d, vx, vy, vz, pre, mfrac) { }
     
   //! Copy constructor (this is needed for the operator overload returns)
   NavierStokes3D_ThermallyPerfect_pState(const NavierStokes3D_ThermallyPerfect_pState &W) {
     species_null(); rho = DENSITY_STDATM; set_initial_values(); Copy(W);
   }

   //! Default destructor
   ~NavierStokes3D_ThermallyPerfect_pState() {
      Deallocate();
   }
//@}
                
/** @name Transport coefficients */
/*        ---------------------- */
//@{
   //! Mixture dynamic viscosity
   double mu(void);    
   double dmudT(void);

   //! Mixture thermal conductivity
   double kappa(void);
   double dkappadT(void);

   //! Species mass diffusion coefficient
   double Ds(const int & i);
   //! Species mass diffusion coefficient
   double Ds(const int & i,
             const double &mu_temp);

   //! Mixture Prandtl number
   double Pr(void);

   //! Species Schmidt number
   double Sc(const int &);

   //! Species Lewis number
   double Le(const int &);
//@}

/** @name Viscous stress tensor and heat flux vector */
/*        ------------------------------------------ */
//@{
   //! Returns (molecular) fluid stress tensor 
   Tensor3D tau(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                const NavierStokes3D_ThermallyPerfect_pState &dWdz);

   //! Returns (molecular) fluid stress tensor 
   Tensor3D tau(const double &mu_temp,
                const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                const NavierStokes3D_ThermallyPerfect_pState &dWdz);

   //! Returns components of (molecular) fluid stress tensor in x-direction 
   Tensor3D tau_x(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                  const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                  const NavierStokes3D_ThermallyPerfect_pState &dWdz);

   //! Returns components of (molecular) fluid stress tensor in x-direction
   Tensor3D tau_x(const double &mu_temp,
                  const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                  const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                  const NavierStokes3D_ThermallyPerfect_pState &dWdz);

   //! Returns components of (molecular) fluid stress tensor in y-direction 
   Tensor3D tau_y(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                  const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                  const NavierStokes3D_ThermallyPerfect_pState &dWdz);

   //! Returns components of (molecular) fluid stress tensor in y-direction
   Tensor3D tau_y(const double &mu_temp,
                  const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                  const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                  const NavierStokes3D_ThermallyPerfect_pState &dWdz);

   //! Returns components of (molecular) fluid stress tensor in z-direction 
   Tensor3D tau_z(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                  const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                  const NavierStokes3D_ThermallyPerfect_pState &dWdz);

   //! Returns components of (molecular) fluid stress tensor in z-direction
   Tensor3D tau_z(const double &mu_temp,
                  const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                  const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                  const NavierStokes3D_ThermallyPerfect_pState &dWdz);

   //! Returns (molecular) heat flux vector 
   Vector3D q(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
              const NavierStokes3D_ThermallyPerfect_pState &dWdy,
              const NavierStokes3D_ThermallyPerfect_pState &dWdz);
   
   //! Returns (molecular) heat flux vector 
   Vector3D q(const double &kappa_temp,
              const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
              const NavierStokes3D_ThermallyPerfect_pState &dWdy,
              const NavierStokes3D_ThermallyPerfect_pState &dWdz);

   //! Returns component of (molecular) heat flux vector in x-direction
   Vector3D q_x(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                const NavierStokes3D_ThermallyPerfect_pState &dWdz);
   
   //! Returns component of (molecular) heat flux vector in x-direction
   Vector3D q_x(const double &kappa_temp,
                const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                const NavierStokes3D_ThermallyPerfect_pState &dWdz);

   //! Returns component of (molecular) heat flux vector in y-direction
   Vector3D q_y(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                const NavierStokes3D_ThermallyPerfect_pState &dWdz);
   
   //! Returns component of (molecular) heat flux vector in y-direction
   Vector3D q_y(const double &kappa_temp,
                const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                const NavierStokes3D_ThermallyPerfect_pState &dWdz);

   //! Returns component of (molecular) heat flux vector in z-direction
   Vector3D q_z(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                const NavierStokes3D_ThermallyPerfect_pState &dWdz);
   
   //! Returns component of (molecular) heat flux vector in z-direction
   Vector3D q_z(const double &kappa_temp,
                const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                const NavierStokes3D_ThermallyPerfect_pState &dWdz);

   //! Returns thermal diffusion flux vector (due to species diffusion processes) 
   Vector3D thermal_diffusion(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                              const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                              const NavierStokes3D_ThermallyPerfect_pState &dWdz) const;

   //! Returns components of thermal diffusion flux vector in x-direction (due to species diffusion processes) 
   Vector3D thermal_diffusion_x(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                                const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                                const NavierStokes3D_ThermallyPerfect_pState &dWdz) const;

   //! Returns components of thermal diffusion flux vector in y-direction (due to species diffusion processes) 
   Vector3D thermal_diffusion_y(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                                const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                                const NavierStokes3D_ThermallyPerfect_pState &dWdz) const;

   //! Returns components of thermal diffusion flux vector in z-direction (due to species diffusion processes) 
   Vector3D thermal_diffusion_z(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                                const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                                const NavierStokes3D_ThermallyPerfect_pState &dWdz) const;

   //! Returns the strain rate tensor
   Tensor3D strain_rate(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                        const NavierStokes3D_ThermallyPerfect_pState &dWdy, 
                        const NavierStokes3D_ThermallyPerfect_pState &dWdz) const;

   //! Returns the flow vorticity tensor
   Tensor3D vorticity(const NavierStokes3D_ThermallyPerfect_pState &dWdx, 
                      const NavierStokes3D_ThermallyPerfect_pState &dWdy, 
                      const NavierStokes3D_ThermallyPerfect_pState &dWdz) const;
//@}

/** @name Viscous flux vectors */
/*        --------------------- */
//@{
   //! x-direction viscous solution flux
   NavierStokes3D_ThermallyPerfect_cState Fv(const NavierStokes3D_ThermallyPerfect_pState &dWdx,
                                             const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                                             const NavierStokes3D_ThermallyPerfect_pState &dWdz);

   //! x-direction viscous solution flux
   NavierStokes3D_ThermallyPerfect_cState Fvx(const NavierStokes3D_ThermallyPerfect_pState &dWdx,
                                              const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                                              const NavierStokes3D_ThermallyPerfect_pState &dWdz);

   //! y-direction viscous solution flux
   NavierStokes3D_ThermallyPerfect_cState Fvy(const NavierStokes3D_ThermallyPerfect_pState &dWdx,
                                              const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                                              const NavierStokes3D_ThermallyPerfect_pState &dWdz);

   //! z-direction viscous solution flux
   NavierStokes3D_ThermallyPerfect_cState Fvz(const NavierStokes3D_ThermallyPerfect_pState &dWdx,
                                              const NavierStokes3D_ThermallyPerfect_pState &dWdy,
                                              const NavierStokes3D_ThermallyPerfect_pState &dWdz);
//@}

/** @name Numerical Evaluation of Viscous Fluxes */
/*        -------------------------------------- */
//@{
   //! Returns viscous flux in n-direction
   static NavierStokes3D_ThermallyPerfect_cState FluxViscous_n(const NavierStokes3D_ThermallyPerfect_pState &Wl,
                                                               const NavierStokes3D_ThermallyPerfect_pState &Wr,
                                                               const NavierStokes3D_ThermallyPerfect_pState &W1 ,
                                                               const NavierStokes3D_ThermallyPerfect_pState &W2,
                                                               const NavierStokes3D_ThermallyPerfect_pState &dWdx1,
                                                               const NavierStokes3D_ThermallyPerfect_pState &dWdy1,
                                                               const NavierStokes3D_ThermallyPerfect_pState &dWdz1,
                                                               const NavierStokes3D_ThermallyPerfect_pState &dWdx2,
                                                               const NavierStokes3D_ThermallyPerfect_pState &dWdy2,
                                                               const NavierStokes3D_ThermallyPerfect_pState &dWdz2,
                                                               const Vector3D &norm, 
                                                               const Vector3D &ts, 
                                                               const double &deltad, 
                                                               const double &Volume, 
                                                               const double &Volume_Neigbor);
//@}
  
/** @name Operators */
/*        --------- */
//@{
   //! Binary multiplication operator
   friend NavierStokes3D_ThermallyPerfect_pState operator *(const double &a, 
                                                            const NavierStokes3D_ThermallyPerfect_pState &W);
  
   //! Unary subtraction operators
   friend NavierStokes3D_ThermallyPerfect_pState operator -(const NavierStokes3D_ThermallyPerfect_pState &W);

   //! Equal relational operator
   friend int operator ==(const NavierStokes3D_ThermallyPerfect_pState &W1,
                          const NavierStokes3D_ThermallyPerfect_pState &W2);

   //! Not equal relational operator
   friend int operator !=(const NavierStokes3D_ThermallyPerfect_pState &W1,
                          const NavierStokes3D_ThermallyPerfect_pState &W2);

   //! Output stream operator
   friend ostream& operator << (ostream &out_file,
                                const NavierStokes3D_ThermallyPerfect_pState &W);

   //! Input stream operator
   friend istream& operator >> (istream &in_file,
                                NavierStokes3D_ThermallyPerfect_pState &W);
//@}
};

/*!
 * Class: NavierStokes3D_TheramllyPerfect_cState
 *
 * \brief Conserved state solution class for 3D Navier-Stokes equations
 *        governing flows of thermally perfect non-reactive and 
 *        combusting mixtures.
 *
 * Member functions
 *  - rho     -- Return mixture density (kg/m^3)
 *  - rhov    -- Return mixture momentum (kg/(m^2-s))   
 *  - E       -- Return mixture total energy (J/kg)
 *  - rhospec -- Return array of species density data
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
class NavierStokes3D_ThermallyPerfect_cState : public Euler3D_ThermallyPerfect_cState {
  public:

/** @name Constructors and destructors */
/*        ---------------------------- */
//@{
   //! Default creation constructor (assign default values)
   NavierStokes3D_ThermallyPerfect_cState(void): Euler3D_ThermallyPerfect_cState() { }

   //! Constructor from base class (allows return of derived type)
   NavierStokes3D_ThermallyPerfect_cState(const Euler3D_ThermallyPerfect_cState &U1) : 
     Euler3D_ThermallyPerfect_cState(U1) { }
 
   //! Assignment constructor
   NavierStokes3D_ThermallyPerfect_cState(const double &value) : 
     Euler3D_ThermallyPerfect_cState(value) { }
  
   //! Assignment constructor
   NavierStokes3D_ThermallyPerfect_cState(const double &d, 
                                          const Vector3D &V,
                                          const double &En) :
      Euler3D_ThermallyPerfect_cState(d, V, En) { }

   //! Assignment constructor
   NavierStokes3D_ThermallyPerfect_cState(const double &d, 
                                          const Vector3D &V,
                                          const double &En, 
                                          const double &frac) :
      Euler3D_ThermallyPerfect_cState(d, V, En, frac) {}

   //! Assignment constructor
   NavierStokes3D_ThermallyPerfect_cState(const double &d, 
                                          const Vector3D &dV,
                                          const double &En, 
                                          const Species *rhomfrac) :
      Euler3D_ThermallyPerfect_cState(d, dV, En, rhomfrac) { }

   //! Assignment constructor
   NavierStokes3D_ThermallyPerfect_cState(const double &d, 
                                          const double &vx,
                                          const double &vy, 
                                          const double &vz,
                                          const double &En) :
      Euler3D_ThermallyPerfect_cState(d, vx, vy, vz, En) { }
   
   //! Assignment constructor
   NavierStokes3D_ThermallyPerfect_cState(const double &d, 
                                          const double &vx, 
                                          const double &vy, 
                                          const double &vz, 
                                          const double &En,	
                                          const double &rhomfrac) :
      Euler3D_ThermallyPerfect_cState(d, vx, vy, vz, En, rhomfrac) { }

   //! Assignment constructor
   NavierStokes3D_ThermallyPerfect_cState(const double &d, 
                                          const double &vx,
                                          const double &vy, 
                                          const double &vz,
                                          const double &En,
                                          const Species *rhomfrac) :
      Euler3D_ThermallyPerfect_cState(d, vx, vy, vz, En, rhomfrac) { }

   //! Copy constructor (this is needed for the operator overload returns)
   NavierStokes3D_ThermallyPerfect_cState(const NavierStokes3D_ThermallyPerfect_cState &U) { 
     rhospec_null(); rho = DENSITY_STDATM; set_initial_values(); Copy(U); 
   }

   //! Default destructor
   ~NavierStokes3D_ThermallyPerfect_cState() {
      Deallocate();
   }
//@}
                          
/** @name Transport coefficients */
/*        ---------------------- */
//@{
   //! Mixture dynamic viscosity
   double mu(void);    
   double dmudT(void);

   //! Mixture thermal conductivity
   double kappa(void);
   double dkappadT(void);

   //! Species mass diffusion coefficient
   double Ds(const int &i);
   //! Species mass diffusion coefficient
   double Ds(const int &i,
             const double &mu_temp);
//@}

/** @name Operators */
/*        --------- */
//@{
   //! Binary multiplication operator
   friend NavierStokes3D_ThermallyPerfect_cState operator *(const double &a, 
                                                            const NavierStokes3D_ThermallyPerfect_cState &U);
  
   //! Unary subtraction operators
   friend NavierStokes3D_ThermallyPerfect_cState operator -(const NavierStokes3D_ThermallyPerfect_cState &U);

   //! Equal relational operator
   friend int operator ==(const NavierStokes3D_ThermallyPerfect_cState &U1,
                          const NavierStokes3D_ThermallyPerfect_cState &U2);

   //! Not equal relational operator
   friend int operator !=(const NavierStokes3D_ThermallyPerfect_cState &U1,
                          const NavierStokes3D_ThermallyPerfect_cState &U2);

   //! Output stream operator
   friend ostream& operator << (ostream &out_file,
                                const NavierStokes3D_ThermallyPerfect_cState &U);

   //! Input stream operator
   friend istream& operator >> (istream &in_file,
                                NavierStokes3D_ThermallyPerfect_cState &U);
//@}
};

/****************************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState member functions.                             *
 ****************************************************************************************/

/*********************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState -- Binary arithmetic operators.        *
 *********************************************************************************/
//---------------- Scalar Multiplication ------------------//
inline NavierStokes3D_ThermallyPerfect_pState operator *(const double &a, 
                                                         const NavierStokes3D_ThermallyPerfect_pState &W) {
   NavierStokes3D_ThermallyPerfect_pState Temp;
   //Temp.Copy(W);
   Temp.rho = W.rho*a;  Temp.v = W.v*a; Temp.p = W.p*a;
   for (int i = 0; i < W.ns; i++) {
      Temp.spec[i] = W.spec[i]*a;
   } /* endfor */
   return(Temp);
}

/*********************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState -- Unary arithmetic operators.         *
 *********************************************************************************/
inline NavierStokes3D_ThermallyPerfect_pState operator -(const NavierStokes3D_ThermallyPerfect_pState &W) {
   Species *spt= new Species[W.ns];
   for (int i = 0; i < W.ns; i++) {
      spt[i] = -W.spec[i]; 
   } /* endfor */ 
   NavierStokes3D_ThermallyPerfect_pState Temp(-W.rho,-W.v, -W.p, spt);
   delete[] spt;
   return(Temp);
}

/*********************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState -- Relational operators.               *
 *********************************************************************************/
inline int operator ==(const NavierStokes3D_ThermallyPerfect_pState &W1, 
                       const NavierStokes3D_ThermallyPerfect_pState &W2) {
   bool Temp;
   for (int i = 0; i < W1.ns; i++) {
      if (W1.spec[i] == W2.spec[i]) {
         Temp = true;
      } else {
         Temp = false;
         break;
      } /* endif */
   } /* endfor */
   return (W1.rho == W2.rho && W1.v == W2.v && W1.p == W2.p && Temp == true);
}

inline int operator !=(const NavierStokes3D_ThermallyPerfect_pState &W1, 
                       const NavierStokes3D_ThermallyPerfect_pState &W2) {
   bool Temp = true;
   for (int i = 0; i < W1.ns; i++) {
      if (W1.spec[i] != W2.spec[i]) {
         Temp = false;
         break;
      } /* endif */
   } /* endfor */
   return (W1.rho != W2.rho || W1.v != W2.v || W1.p != W2.p || Temp != true);
}

/*********************************************************************************
 * NavierStokes3D_ThermallyPerfect_pState -- Input-output operators.             *
 *********************************************************************************/
inline ostream &operator << (ostream &out_file, 
                             const NavierStokes3D_ThermallyPerfect_pState &W) {
   out_file.precision(10);
   out_file.setf(ios::scientific);
   out_file << " " << W.rho  << " " << W.v.x << " " << W.v.y 
            << " " << W.v.z << " " << W.p;
   for (int i = 0; i < W.ns; i++) {
      out_file<<" "<<W.spec[i];
   } /* endfor */
   out_file.unsetf(ios::scientific);
   return (out_file);
}

inline istream &operator >> (istream &in_file, 
                             NavierStokes3D_ThermallyPerfect_pState &W) {
   in_file.setf(ios::skipws);
   in_file >> W.rho >> W.v.x >> W.v.y >> W.v.z >>W.p;
   //W.set_initial_values();
   for (int i=0; i < W.ns; i++) {
      in_file>>W.spec[i];
   } /* endfor */
   in_file.unsetf(ios::skipws);
   return (in_file);
}

/****************************************************************************************
 * NavierStokes3D_ThermallyPerfect_cState member functions.                             *
 ****************************************************************************************/

/*********************************************************************************
 * NavierStokes3D_ThermallyPerfect_cState -- Binary arithmetic operators.        *
 *********************************************************************************/
//---------------- Scalar Multiplication ------------------//
inline NavierStokes3D_ThermallyPerfect_cState operator *(const double &a, 
                                                         const NavierStokes3D_ThermallyPerfect_cState &U) {
   NavierStokes3D_ThermallyPerfect_cState Temp;
   Temp.rho = U.rho*a;  
   Temp.rhov = U.rhov*a; 
   Temp.E = U.E*a;
   for (int i = 0; i < U.ns; i++) {
      Temp.rhospec[i] = U.rhospec[i]*a;
   } /* endfor */
   return(Temp);
}

/*********************************************************************************
 * NavierStokes3D_ThermallyPerfect_cState -- Unary arithmetic operators.         *
 *********************************************************************************/
inline NavierStokes3D_ThermallyPerfect_cState operator -(const NavierStokes3D_ThermallyPerfect_cState &U) {
   Species *spt= new Species[U.ns];
   for (int i = 0; i < U.ns; i++) {
      spt[i] = -U.rhospec[i]; 
   } /* endfor */
   NavierStokes3D_ThermallyPerfect_cState Temp(-U.rho,-U.rhov,-U.E, spt);
   delete[] spt;
   return(Temp);
}

/*********************************************************************************
 * NavierStokes3D_ThermallyPerfect_cState -- Relational operators.               *
 *********************************************************************************/
inline int operator ==(const NavierStokes3D_ThermallyPerfect_cState &U1, 
                       const NavierStokes3D_ThermallyPerfect_cState &U2) {
   bool Temp;
   for (int i = 0; i < U1.ns; i++) {
      if (U1.rhospec[i] == U2.rhospec[i] ){
         Temp = true;
      } else {
         Temp = false;
         break;
      } /* endif */
   } /* endfor */
   return (U1.rho == U2.rho && U1.rhov == U2.rhov && U1.E == U2.E && Temp == true);
}

inline int operator !=(const NavierStokes3D_ThermallyPerfect_cState &U1, 
                       const NavierStokes3D_ThermallyPerfect_cState &U2) {
   bool Temp = true;
   for (int i = 0; i < U1.ns; i++) {
      if (U1.rhospec[i] != U2.rhospec[i] ){
         Temp = false;
         break;
      } /* endif */
    } /* endfor */
   return (U1.rho != U2.rho || U1.rhov != U2.rhov || U1.E != U2.E || Temp != true);
}

/*********************************************************************************
 * NavierStokes3D_ThermallyPerfect_cState -- Input-output operators.             *
 *********************************************************************************/
inline ostream &operator << (ostream &out_file, 
                             const NavierStokes3D_ThermallyPerfect_cState &U) {
   //out_file.precision(20);
   out_file.setf(ios::scientific);
   out_file << " " << U.rho  << " " << U.rhov.x << " " << U.rhov.y << " " << U.rhov.z << " "<< U.E;
   for (int i = 0; i < U.ns; i++) {
      out_file << " " << U.rhospec[i];
   } /* endfor */
   out_file.unsetf(ios::scientific);
   return (out_file);
}

inline istream &operator >> (istream &in_file, 
                             NavierStokes3D_ThermallyPerfect_cState &U) {
   in_file.setf(ios::skipws);
   in_file >> U.rho >> U.rhov.x >> U.rhov.y >> U.rhov.z >> U.E;
   //U.set_initial_values();
   for (int i = 0; i < U.ns; i++) {
     in_file>>U.rhospec[i]; 
   } /* endfor */
   in_file.unsetf(ios::skipws);
   return (in_file);
}

#endif // _NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED 
