/*! \file NavierStokes3DPolytropicState.h
 * 	\brief	Header file defining the Navier-Stokes solution state classes 
 *              associated with solution of compressible viscous flows 
 *              of a polytropic gas.
 */

#ifndef _NAVIERSTOKES3D_POLYTROPOIC_STATE_INCLUDED 
#define _NAVIERSTOKES3D_POLYTROPIC_STATE_INCLUDED

/* Include header file for base solution classes from which the classes are derived. */

#ifndef _EULER3D_POLYTROPIC_STATE_INCLUDED
#include "../Euler/Euler3DPolytropicState.h"
#endif  //EULER3D_POLYTROPIC_STATE_INCLUDED

/* Define the classes. */

class NavierStokes3D_Polytropic_pState;
class NavierStokes3D_Polytropic_cState;

/*!
 * Class: NavierStokes3D_Polytropic_pState
 *
 * \brief Primitive state solution class for 3D Navier-Stokes equations
 *        governing flows of polytropic gases.
 *
 
 * Member functions
 *  - rho           -- Return density.
 *  - v             -- Return flow velocity.
 *  - p             -- Return pressure.
 *  - g             -- Return specific heat ratio.
 *  - gm1           -- Return g-1
 *  - gm1i          -- Return 1/(g-1).
 *  - R             -- Return gas constant.
 *  - setgas        -- Set gas constants.
 *  - T             -- Return temperature.
 *  - e             -- Return specific internal energy.
 *  - E             -- Return total energy.
 *  - h             -- Return specific enthalpy.
 *  - H             -- Return total enthalpy.
 *  - a             -- Return sound speed.
 *  - a2            -- Return sound speed square.
 *  - M             -- Return Mach number.
 *  - s             -- Return specific entropy.
 *  - rhov          -- Return momentum.
 *  - To            -- Return stagnation temperature.
 *  - po            -- Return stagnation pressure.
 *  - ao            -- Return stagnation sound speed.
 *  - ho            -- Return stagnation enthalpy.
 *  - U             -- Return conserved solution state.
 *  - F, Fx         -- Return x-direction solution flux.
 *  - Fy            -- Return y-direction solution flux.
 *  - Fz            -- Return z-direction solution flux.
 *  - dFxdU         -- Return x-direction jacobian.
 *  - dFydU         -- Return y-direction jacobian.
 *  - dFzdU         -- Return z-direction jacobian.
 *  - dUdW          -- Return solution variable jacobian.
 *  - dWdU          -- Return solution variable jacobian.
 *  - lambda        -- Return x-direction eigenvalue.
 *  - lambda_x      -- Return x-direction eigenvalue.
 *  - lambda_y      -- Return y-direction eigenvalue.
 *  - lambda_z      -- Return z-direction eigenvalue.
 *  - rp, rp_x      -- Return x-direction primitive right eigenvector.
 *  - rp_y          -- Return y-direction primitive right eigenvector.
 *  - rp_z          -- Return z-direction primitive right eigenvector.
 *  - rc, rc_x      -- Return x-direction conserved right eigenvector.
 *  - rc_y          -- Return y-direction conserved right eigenvector.
 *  - rc_z          -- Return z-direction conserved right eigenvector.
 *  - lp, lp_x      -- Return x-direction primitive left eigenvector.
 *  - lp_y          -- Return y-direction primitive left eigenvector.
 *  - lp_z          -- Return z-direction primitive left eigenvector.
 *  - RoeAverage    -- Return Roe average.
 *  - lambda_minus  -- Return negative eigenvalues, applying Harten entropy fix
 *  - lambda_plus   -- Return positive eigenvalues, applying Harten entropy fix
 *  - FluxRoe_x     -- Return Roe's solution of Riemann problem in x-direction.
 *  - FluxRoe_y     -- Return Roe's solution of Riemann problem in y-direction.
 *  - FluxRoe_z     -- Return Roe's solution of Riemann problem in z-direction.
 *  - FluxRoe_n     -- Return Roe's solution of Riemann problem in n-direction.
 *  - FluxHLLE_x    -- Return HLLE solution of Riemann problem in x-direction.
 *  - FluxHLLE_y    -- Return HLLE solution of Riemann problem in y-direction.
 *  - FluxHLLE_z    -- Return HLLE solution of Riemann problem in z-direction.
 *  - FluxHLLE_n    -- Return HLLE solution of Riemann problem in n-direction.
 *  - Reflect       -- Return reflected solution state after application of reflection BC
 *  - MovingWall    -- Return wall solution state after application of moving wall BC
 *  - NoSlip        -- Return wall solution state after application of no-slip BC
 * Extra in NavierStokes:
 *  - Cp            -- Return specific heat at constant pressure (J/(kg*K))
 *  - Cv            -- Return specific heat at constant volume (J/(kg*K))
 *  - mu            -- Return mixture dynamic viscosity
 *  - kappa         -- Return mixture coefficient of thermal conductivity
 *  - tau           -- Return (molecular) fluid stress tensor
 *  - tau_x         -- Return components of (molecular) fluid stress tensor in x-direction
 *  - tau_y         -- Return components of (molecular) fluid stress tensor in y-direction
 *  - tau_z         -- Return components of (molecular) fluid stress tensor in z-direction
 *  - q             -- Return (molecular) heat flux vector
 *  - q_x           -- Return component of (molecular) heat flux vector in x-direction
 *  - q_y           -- Return component of (molecular) heat flux vector in y-direction
 *  - q_z           -- Return component of (molecular) heat flux vector in z-direction 
 *  - Fv, Fvx       -- Return x-direction viscous solution flux
 *  - Fvy           -- Return y-direction viscous solution flux
 *  - Fvz           -- Return z-direction viscous solution flux
 *  - FluxViscous_n -- Returns viscous flux in n-direction
 *
 * Member operators \n
 *  W         -- a primitive solution state \n
 *  c         -- a scalar (double)
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
 *  - W = -W;
 *  - W += W;
 *  - W -= W;
 *  - W *= W;
 *  - W /= W;
 *  - W == W;
 *  - W != W;
 *  - W <= W;
 *  - W >= W;
 *  - W < W;
 *  - W > W;
 *  - cout << W; (output function)
 *  - cin  >> W; (input function)
 *  \nosubgrouping
 */

class NavierStokes3D_Polytropic_pState: public Euler3D_Polytropic_pState {
public :
    
    /*
     * Constructors
     * ------------
     */
	NavierStokes3D_Polytropic_pState() : Euler3D_Polytropic_pState() { }
	NavierStokes3D_Polytropic_pState(const double &rho, 
									 const Vector3D &v, 
									 const double &p) : Euler3D_Polytropic_pState(rho,v,p) { }
	NavierStokes3D_Polytropic_pState(const double &rho, 
									 const double &vx, const double &vy, const double &vz, 
									 const double &p) : Euler3D_Polytropic_pState(rho,vx,vy,vz,p) { }
    NavierStokes3D_Polytropic_pState(const double &value) : Euler3D_Polytropic_pState(value) { }
    
	NavierStokes3D_Polytropic_pState(const Euler3D_Polytropic_pState &W) : Euler3D_Polytropic_pState(W) { }
//    NavierStokes3D_Polytropic_pState(const Euler3D_Polytropic_cState &U) : Euler3D_Polytropic_pState(U) { }
	
    /*
     * Gas specific constants
     * ----------------------
     */
	static double v1;	//!< Viscosity law coefficient (only for Navier-Stokes)
    static double v2;	//!< Viscosity law coefficient (only for Navier-Stokes)
	static double v3;	//!< Viscosity law coefficient (only for Navier-Stokes)
	static double v4;	//!< Viscosity law coefficient (only for Navier-Stokes)
    static double v5;	//!< Viscosity law coefficient (only for Navier-Stokes)
	static double Cp;   //!< Heat capacity with constant pressure
	static double Cv;   //!< Heat capacity with constant volume
    
    /* 
     * Set gas specific constants
     * --------------------------
     */
    void setgas(void);              //!< Set gas to air
	void setgas(char *str_ptr);     //!< Set gas constants
	
	/*
     * Viscous- and Heat Fluxes
     * ------------------------
     */
    
    //! Viscous flux and Heat flux in x-direction
	NavierStokes3D_Polytropic_cState Fvx(const NavierStokes3D_Polytropic_pState &dWdx,
                                         const NavierStokes3D_Polytropic_pState &dWdy,
                                         const NavierStokes3D_Polytropic_pState &dWdz);
    
    NavierStokes3D_Polytropic_cState Fv(const NavierStokes3D_Polytropic_pState &dWdx,
                                        const NavierStokes3D_Polytropic_pState &dWdy,
                                        const NavierStokes3D_Polytropic_pState &dWdz);

    //! Viscous flux and Heat flux in y-direction
	NavierStokes3D_Polytropic_cState Fvy(const NavierStokes3D_Polytropic_pState &dWdx,
                                         const NavierStokes3D_Polytropic_pState &dWdy,
                                         const NavierStokes3D_Polytropic_pState &dWdz);

    //! Viscous flux and Heat flux in z-direction
	NavierStokes3D_Polytropic_cState Fvz(const NavierStokes3D_Polytropic_pState &dWdx,
										const NavierStokes3D_Polytropic_pState &dWdy,
										const NavierStokes3D_Polytropic_pState &dWdz);
    
    //! Viscous flux and Heat flux in n-direction
	static NavierStokes3D_Polytropic_cState FluxViscous_n(const NavierStokes3D_Polytropic_pState &Wl,
												   const NavierStokes3D_Polytropic_pState &Wr,
												   const NavierStokes3D_Polytropic_pState &W1 ,
												   const NavierStokes3D_Polytropic_pState &W2,
												   const NavierStokes3D_Polytropic_pState &dWdx1,
												   const NavierStokes3D_Polytropic_pState &dWdy1,
												   const NavierStokes3D_Polytropic_pState &dWdz1,
												   const NavierStokes3D_Polytropic_pState &dWdx2,
												   const NavierStokes3D_Polytropic_pState &dWdy2,
												   const NavierStokes3D_Polytropic_pState &dWdz2,
												   const Vector3D &norm, 
												   const Vector3D &ts, 
												   const double &deltad, 
												   const double &Volume, 
												   const double &Volume_Neigbour);
	/*
     * Navier-Stokes related functions
     * -------------------------------
     */
    
    double mu(void);        //!< Dynamic Viscosity
    double nu(void);        //!< Kinematic Viscosity
	double kappa(void);     //!< Thermal Conductivity
    
	//! Molecular stress tensor
	Tensor3D tau(const NavierStokes3D_Polytropic_pState &dWdx, 
                 const NavierStokes3D_Polytropic_pState &dWdy,
                 const NavierStokes3D_Polytropic_pState &dWdz);
	
	Tensor3D tau_x(const NavierStokes3D_Polytropic_pState &dWdx, 
                   const NavierStokes3D_Polytropic_pState &dWdy,
                   const NavierStokes3D_Polytropic_pState &dWdz);
	
	Tensor3D tau_y(const NavierStokes3D_Polytropic_pState &dWdx, 
                   const NavierStokes3D_Polytropic_pState &dWdy,
                   const NavierStokes3D_Polytropic_pState &dWdz);

	Tensor3D tau_z(const NavierStokes3D_Polytropic_pState &dWdx, 
                   const NavierStokes3D_Polytropic_pState &dWdy,
                   const NavierStokes3D_Polytropic_pState &dWdz);
	
	//! Heat flux vector 
	Vector3D q(const NavierStokes3D_Polytropic_pState &dWdx, 
               const NavierStokes3D_Polytropic_pState &dWdy,
               const NavierStokes3D_Polytropic_pState &dWdz);
	
	Vector3D q_x(const NavierStokes3D_Polytropic_pState &dWdx, 
                 const NavierStokes3D_Polytropic_pState &dWdy,
                const NavierStokes3D_Polytropic_pState &dWdz);
	
	Vector3D q_y(const NavierStokes3D_Polytropic_pState &dWdx, 
                 const NavierStokes3D_Polytropic_pState &dWdy,
                 const NavierStokes3D_Polytropic_pState &dWdz);
	
	Vector3D q_z(const NavierStokes3D_Polytropic_pState &dWdx, 
                 const NavierStokes3D_Polytropic_pState &dWdy,
                 const NavierStokes3D_Polytropic_pState &dWdz);
};












/* -------------------------------------------------------------------------- *
 *                     NavierStokes3D_Polytropic_cState                       *
 * -------------------------------------------------------------------------- */
/*! \brief This is the conservative state class for Navier-Stokes 3D 
 *			Polytropic flows.
 *
 *	This class inherits and builds on the class Euler3D_Polytropic_cState.
 *  \nosubgrouping
 */
class NavierStokes3D_Polytropic_cState: public Euler3D_Polytropic_cState {
public :
    
    /*
     * Constructors
     * ------------
     */
	NavierStokes3D_Polytropic_cState() : Euler3D_Polytropic_cState() { }
	NavierStokes3D_Polytropic_cState(const double &rho, 
									 const Vector3D &rhov, 
									 const double &E) : Euler3D_Polytropic_cState(rho,rhov,E) { }
	NavierStokes3D_Polytropic_cState(const double &rho, 
									 const double &rhovx, const double &rhovy, const double &rhovz, 
									 const double &E) : Euler3D_Polytropic_cState(rho,rhovx,rhovy,rhovz,E) { }
    NavierStokes3D_Polytropic_cState(const double &value) : Euler3D_Polytropic_cState(value) { }
    
	NavierStokes3D_Polytropic_cState(const Euler3D_Polytropic_cState &U) : Euler3D_Polytropic_cState(U) { }
 //   NavierStokes3D_Polytropic_cState(const Euler3D_Polytropic_pState &W) : Euler3D_Polytropic_cState(W) { }
	
    /*
     * Gas specific constants
     * ----------------------
     */
	static double v1;	//!< Viscosity law coefficient (only for Navier-Stokes)
    static double v2;	//!< Viscosity law coefficient (only for Navier-Stokes)
	static double v3;	//!< Viscosity law coefficient (only for Navier-Stokes)
	static double v4;	//!< Viscosity law coefficient (only for Navier-Stokes)
    static double v5;	//!< Viscosity law coefficient (only for Navier-Stokes)
	static double Cp;   //!< Heat capacity with constant pressure
	static double Cv;   //!< Heat capacity with constant volume
    
    /* 
     * Set gas specific constants
     * --------------------------
     */
    void setgas(void);              //!< Set gas to air
	void setgas(char *str_ptr);     //!< Set gas constants
		
};

#endif // _NAVIERSTOKES3D_POLYTROPIC_STATE_INCLUDED

