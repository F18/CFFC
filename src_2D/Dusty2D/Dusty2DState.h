/**********************************************************************
 * Dusty2DState.h: Header file defining 2D Dusty solution state       *
 *                 classes.                                           *
 **********************************************************************/

#ifndef _DUSTY2D_STATE_INCLUDED
#define _DUSTY2D_STATE_INCLUDED

// Include required C++ libraries.
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

// Include math macro, CFD, 2D vector, and gas constant header files.

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _MATRIX_INCLUDED
#include "../Math/Matrix.h"
#endif // _MATRIX_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _VECTOR2D_INCLUDED
#include "../Math/Vector2D.h"
#endif //_VECTOR2D_INCLUDED

#ifndef _TENSOR2D_INCLUDED
#include "../Math/Tensor2D.h"
#endif //_TENSOR2D_INCLUDED

#ifndef _GAS_CONSTANTS_INCLUDED
#include "../Physics/GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED

// Include the particle state and components class header files.

#ifndef _PARTICLE2D_STATE_INCLUDED
#include "Particle2DState.h"
#endif // _PARTICLE2D_STATE_INCLUDED

#ifndef _PARTICLE2D_COMPONENTS_INCLUDED
#include "Particle2DComponents.h"
#endif // _PARTICLE2D_COMPONENTS_INCLUDED

// Include the electrostatic state class header file.

#ifndef _ELECTROSTATIC2D_STATE_INCLUDED
#include "../Electrostatic2D/Electrostatic2DState.h"
#endif // _ELECTROSTATIC2D_STATE_INCLUDED

// Define parameters for the various drag laws.
#define DRAG_LAW_STOKES                     10
#define DRAG_LAW_OSEEN                      11
#define DRAG_LAW_KLYACHKO_PUTNAM            12
#define DRAG_LAW_TURTON_LEVENSPIEL          13
#define DRAG_LAW_HAIDER_LEVENSPIEL          14

// Define the number of base variables (gas-phase).
#define	NUM_VAR_BASE  6

// Define the classes.

class Dusty2D_cState;

/*!
 * Class: Dusty2D_pState
 *
 * @brief Primitive variable solution state class definition for a 
 * gas-particle flow.
 *
 * Primitive variable solution state class definition for a gas-particle
 * flow.  The particle-phase is modelled using an Eulerian formulation
 * and is asumed to be inert, disperse and dilute particle-phase.  The
 * gas-phase can be inviscid, laminar, or turbulent.  The two-equation 
 * k-omega turbulence model is used to achieve turbulence closure.
 * Electrostatic effects on the particle-phase can also be included.
 *
 * \verbatim
 * Member functions
 *     rho      -- Gas-phase density.
 *     v        -- Gas-phase velocity.
 *     p        -- Gas-phase pressure.
 *     k        -- Gas turbulent kinetic energy.
 *     omega    -- Gas specific dissipation rate.
 *     tau      -- Viscous stress tensor (laminar and Reynolds).
 *     q        -- Heat flux vector (laminar and turbulent).
 *     g        -- Gas specific heat ratio.
 *     gm1      -- Return g-1.
 *     gm1i     -- Return 1/(g-1).
 *     R        -- Gas constant.
 *     v1,v2,v3,v4,v5 -- Viscosity law coefficients.
 *     cp       -- Specific heat at constant pressure.
 *     T        -- Gas-phase temperature.
 *     e        -- Gas-phase specific internal energy.
 *     E        -- Gas-phase total energy.
 *     h        -- Gas specific enthalpy.
 *     H        -- Gas-phase total enthalpy.
 *     a        -- Gas-phase sound speed.
 *     a2       -- Gas-phase sound speed squared.
 *     M        -- Gas-phase Mach number.
 *     s        -- Gas-phase specific entropy.
 *     dv       -- Gas-phase momentum.
 *     To       -- Gas-phase stagnation temperature.
 *     po       -- Gas-phase stagnation pressure.
 *     ao       -- Gas-phase stagnation sound speed.
 *     ho       -- Gas-phase stagnation enthalpy.
 *     mu       -- Gas-phase dynamic viscosity.
 *     nu       -- Gas-phase kinematic viscosity.
 *     kappa    -- Gas-phase thermal conductivity.
 *     Pr       -- Prandtl number.
 *     meanfreepath -- Gas-phase mean free path.
 *
 *     Wp       -- Array of particle primitive states.
 *     NUM_CMP_PART -- Number of particle-phase components.
 *     NUM_VAR_PART -- Number of particle-phase variables.
 *
 *     muT      -- Turbulent eddy dynamic viscosity.
 *     nuT      -- Turbulent eddy kinematic viscosity.
 *     kappaT   -- Turbulent eddy thermal conductivity.
 *     PrT      -- Turbulent Prandtl number.
 *     epsilon  -- Return the turbulent eddy dissipation.
 *     c        -- Turbulence modified sound speed.
 *     cc       -- Turbulence modified sound speed squared.
 *     pmodified -- Turbulence modified pressure.
 *
 *     rhos     -- Propellant density.
 *     n        -- Propellant burning rate constant.
 *     beta     -- Propellant burning rate coefficient.
 *     Tf       -- Propellant flame temperature.
 *     Ts       -- Propellant surface temperature.
 *     alphas   -- Propellant particle concentration.
 *
 *     NUM_VAR_DUSTY2D -- Total number of variables.
 *
 *     tauv     -- Momentum transfer relaxation time scale.
 *     tauT     -- Heat transfer relaxation time scale.
 *     Rep      -- Slip Reynolds number.
 *     fRep     -- Stokes drag correction.
 *     Kn       -- Particle Knudsen number.
 *     Cc       -- Cunningham correction.
 *     Nu       -- Nusselt number.
 *     burningrate -- Burning rate.
 *
 *     U        -- Return the conserved solution state.
 *     dUdW     -- Return the Jacobian of the conserved solution
 *                 variables with respect to the primitive solution
 *                 variables.
 *     dWdU     -- Return the Jacobian of the primitive solution
 *                 variables with respect to the conserved solution
 *                 variables.
 *
 *     F        -- Return x-direction inviscid solution flux.
 *     dFdU     -- Return the Jacobian of the inviscid solution flux
 *                 vector with respect to the conserved solution
 *                 variables.
 *
 *     Gx       -- Return x-direction viscous solution flux.
 *     Gy       -- Return y-direction viscous solution flux.
 *     dGxdU    -- Return the Jacobian of the x-direction viscous
 *                 solution flux vector with respect to the conserved
 *                 solution variables.
 *     dGydU    -- Return the Jacobian of the y-direction viscous
 *                 solution flux vector with respect to the conserved
 *                 solution variables.
 *
 *     lambda_x -- Return x-direction eigenvalue(s).
 *     rp_x     -- Return primitive right eigenvector (x-direction).
 *     rc_x     -- Return conserved right eigenvector (x-direction).
 *     lp_x     -- Return primitive left eigenvector (x-direction).
 *
 *     Si       -- Return inviscid axisymmetric source term vector.
 *     dSidU    -- Return the Jacobian of the inviscid axisymmetric
 *                 source term vector with respect to the conserved
 *                 solution variables.
 *
 *     Sv       -- Return viscous axisymmetric source term vector.
 *     dSvdU    -- Return the Jacobian of the viscous axisymmetric 
 *                 source term vector with respect to the conserved
 *                 solution variables.
 *
 *     Sp       -- Return phase interaction source term vector.
 *     dSpdU    -- Return the Jacobian of the phase interaction
 *                 source term vector with respect to the conserved
 *                 solution variables.
 *
 *     Se       -- Return electric field source term vector.
 *     dSedU    -- Return the Jacobian of the electric field source
 *                 term vector with respect to the conserved solution
 *                 variables.
 *
 * Member operators
 *      W -- a primitive solution state
 *      c -- a scalar (double)
 *
 * W = W;
 * c = W[i];
 * W = W + W;
 * W = W - W;
 * c = W * W; (inner product)
 * W = c * W;
 * W = W * c;
 * W = W / c;
 * W = W ^ W; (a useful product)
 * W = +W;
 * W = -W;
 * W += W;
 * W -= W;
 * W == W;
 * W != W;
 * cout << W; (output function)
 * cin  >> W; (input function)
 * \endverbatim
 */
class Dusty2D_pState {
 private:
 public:
  //@{ @name Gas-phase primitive variables:
  double                       rho; //!< Gas-phase density.
  Vector2D                       v; //!< Gas-phase velocity (2D vector).
  double                         p; //!< Gas-phase pressure.
  double                         k; //!< Gas-phase turbulent kinetic energy.
  double                     omega; //!< Gas-phase specific dissipation rate.
  Tensor2D                     tau; //!< Viscous stess tensor (laminar and turbulent).
  Vector2D                       q; //!< Heat flux vector (laminar and turbulent).
  static double                  g; //!< Specific heat ratio.
  static double                gm1; //!< g-1
  static double               gm1i; //!< 1/(g-1)
  static double                  R; //!< Gas constant.
  static double                 cp; //!< Specific heat at constant pressure.
  static double v1, v2, v3, v4, v5; //!< Viscosity law coefficients.
  static int             flow_type; //!< Flow-type indicator (inviscid, laminar, or k-omega turbulent).
  //@}

  //@{ @name Particle-phase primitive variables:
  Particle2D_pComponents        Wp; //!< Array of particle primitive components
  static int          NUM_CMP_PART; //!< Total number of particle components.
  static int          NUM_VAR_PART; //!< Total number of particle variables.
  static int              drag_law; //!< Drag law definition.
  static double d1, d2, d3, d4, d5; //!< Drag law coefficients.
  static double                 cm; //!< Characteristic particle specific heat.
  static double                 dp; //!< Characteristic particle diameter.
  static double               rhop; //!< Characteristic particle density.
  static double                 mp; //!< Characteristic particle mass.
  static double                 qe; //!< Characteristic particle charge.
  //@}

  //@{ @name Turbulence boundary-layer constants:
  static double            yplus_o; //!< Transition between viscous sublayer and log layer.
  static double                  C; //!< Surface roughness coefficient.
  static double         von_karman; //!< Von Karman constant.
  static double     yplus_sublayer; //!< Sublayer dimensionless wall distance.
  static double yplus_buffer_layer; //!< Buffer layer dimensionless wall distance.
  static double  yplus_outer_layer; //!< Outer layer dimensionless wall distance.
  //@}

  //@{ @name k-omega closure coefficients:
  static double                PrT; //!< Turbulent Prandtl number.
  static double                Cmu;
  static double           beta_k_o; //!< Destruction of turbulent kinetic energy closure coefficient (a.k.a beta_star).
  static double       beta_omega_o; //!< Destruction of specific dissipation closure coefficient (a.k.a beta).
  static double            sigma_k;
  static double        sigma_omega;
  static double              alpha;
  static double                 xi; //!< Compressiblity correction coefficient.
  static double                Mto; //!< Compressiblity correction coefficient.
  //@}

  //@{ @name Propellant variables:
  static double               rhos; //!< Propellant density.
  static double               beta; //!< Propellant burning rate constant.
  static double                  n; //!< Propellant burning rate coefficiens.
  static double                 Tf; //!< Propellant flame temperature.
  static double                 Ts; //!< Propellant surface temperature.
  static double             alphas; //!< Propellant particle concentration.
  //@}

  //@{ @name General purpose static variables:
  static int       NUM_VAR_DUSTY2D; //<! Total number of variables.
  static int               *cindex; //<! Look up table for the particle-phase component indicing.
  static int               *pindex; //<! Look up table for the particle-phase state indicing.
  //@}

  //@{ @name Creation, copy, and assignment constructors.

  //! Creation constructor.
  Dusty2D_pState(void) {
    allocate();
    rho = DENSITY_STDATM; p = PRESSURE_STDATM;
    //v = Vector2D_ZERO; tau.zero(); q.zero();
    //Standard_Atmosphere();
  }

  //! Copy constructor.
  Dusty2D_pState(const Dusty2D_pState &W) {
    allocate();
    rho = W.rho; v.x = W.v.x; v.y = W.v.y; p = W.p; k = W.k; omega = W.omega;
    if (NUM_VAR_PART) Wp = W.Wp;
  }

  //! Copy constructor.
  Dusty2D_pState(const Dusty2D_cState &U);

  //! Assignment constructor.
  Dusty2D_pState(const double &dens,
		 const Vector2D &V,
		 const double &pre) {
    allocate();
    rho = dens; v.x = V.x; v.y = V.y; p = pre; k = ZERO; omega = ZERO;
  }
  
  //! Assignment constructor.
  Dusty2D_pState(const double &dens,
		 const double &vx,
		 const double &vy,
		 const double &pre) {
    allocate();
    rho = dens; v.x = vx; v.y = vy; p = pre; k = ZERO; omega = ZERO;
  }

  //! Assignment constructor.
  Dusty2D_pState(const double &dens,
		 const Vector2D &V,
		 const double &pre,
		 const Particle2D_pComponents &Wcomp) {
    allocate();
    rho = dens; v.x = V.x; v.y = V.y; p = pre; k = ZERO; omega = ZERO;
    Wp = Wcomp;
  }

  //! Assignment constructor.
  Dusty2D_pState(const double &dens,
		 const double &vx,
		 const double &vy,
		 const double &pre,
		 const Particle2D_pComponents &Wcomp) {
    allocate();
    rho = dens; v.x = vx; v.y = vy; p = pre; k = ZERO; omega = ZERO;
    Wp = Wcomp;
  }

  //! Assignment constructor.
  Dusty2D_pState(const double &dens,
		 const Vector2D &V,
		 const double &pre,
		 const double &kk,
		 const double &omga) {
    allocate();
    rho = dens; v.x = V.x; v.y = V.y; p = pre; k = kk; omega = omga;
  }
  
  //! Assignment constructor.
  Dusty2D_pState(const double &dens,
		 const double &vx,
		 const double &vy,
		 const double &pre,
		 const double &kk,
		 const double &omga) {
    allocate();
    rho = dens; v.x = vx; v.y = vy; p = pre; k = kk; omega = omga;
  }

  //! Assignment constructor.
  Dusty2D_pState(const double &dens,
		 const Vector2D &V,
		 const double &pre,
		 const Particle2D_pComponents &Wcomp,
		 const double &kk,
		 const double &omga) {
    allocate();
    rho = dens; v.x = V.x; v.y = V.y; p = pre; k = kk; omega = omga;
    Wp = Wcomp;
  }

  //! Assignment constructor.
  Dusty2D_pState(const double &dens,
		 const double &vx,
		 const double &vy,
		 const double &pre,
		 const Particle2D_pComponents &Wcomp,
		 const double &kk,
		 const double &omga) {
    allocate();
    rho = dens; v.x = vx; v.y = vy; p = pre; k = kk; omega = omga;
    Wp = Wcomp;
  }

  //! Destructor.
  ~Dusty2D_pState(void) { deallocate(); }

  //@}
  
  //! Return the number of variables.
  int NumVar(void) const { return NUM_VAR_DUSTY2D; }

  //@{ @name Allocation and deallocation functions.
  //! Allocation function.
  void allocate(void);
  //! Deallocation function.
  void deallocate(void);
  //@}

  //@{ @name Set static variables.
  void set_static_variables(void);
  void set_static_variables(char *gas_type,
			    const int &FlowType,
			    const double &C_constant,
			    const double &von_karman_constant,
			    const double &yplus_sub,
			    const double &yplus_buffer,
			    const double &yplus_outer,
			    char *propellant_type);
  void set_static_variables(char *gas_type,
			    char *particle_type,
			    const int &num_cmp_part,
			    const int &draglaw,
			    const int &FlowType,
			    const double &C_constant,
			    const double &von_karman_constant,
			    const double &yplus_sub,
			    const double &yplus_buffer,
			    const double &yplus_outer,
			    char *propellant_type);
  void set_gas(char *gas_type);
  void set_particles(char *particle_type, const int &draglaw);
  void set_turbulence(const double &C_constant,
		      const double &von_karman,
		      const double &yplus_sub,
		      const double &yplus_buffer,
		      const double &yplus_outer);
  void set_propellant(char *propellant_type);
  void set_indicing(void);
  //@}

  //@{ @name Useful operators.
  //! Copy operator.
  void Copy(const Dusty2D_pState &W);

  //! Vacuum operator.
  void Vacuum(void);

  //! Standard atmosphere operator.
  void Standard_Atmosphere(void);

  //! Check for unphysical state properties.
  int Unphysical_Properties(void) const;
  //@}

  //@{ @name Gas-phase related functions.

  //! Gas-phase temperature.
  double T(void) const;
  
  //! Gas-phase specific internal energy.
  double e(void) const;
  
  //! Gas-phase total energy.
  double E(void) const;
  
  //! Gas-phase specific enthalpy.
  double h(void) const;
  
  //! Gas-phase total enthalpy.
  double H(void) const;

  //! Gas-phase sound speed.
  double a(void) const;

  //! Gas-phase sound speed squared.
  double a2(void) const;

  //! Gas-phase Mach number.
  double M(void) const;

  //! Gas-phase specific entropy.
  double s(void) const;

  //! Gas-phase momentum.
  Vector2D dv(void) const;

  //! Gas-phase momentum.
  double dv(const Vector2D &n) const;
  
  //! Gas-phase stagnation temperature.
  double To(void) const;

  //! Gas-phase stagnation pressure.
  double po(void) const;

  //! Gas-phase stagnation sound speed.
  double ao(void) const;

  //! Gas-phase stagnation enthalpy.
  double ho(void) const;

  //! Gas-phase dynamic viscosity.
  double mu(void) const;

  //! Gas-phase kinematic viscosity.
  double nu(void) const;

  //! Gas-phase thermal heat conductivity.
  double kappa(void) const;

  //! Prandtl number.
  double Pr(void) const;

  //! Gas-phase mean free path.
  double meanfreepath(void) const;
  //@}

  //@{ @name Phase interaction functions.
  //! Momentum transfer relaxation time.
  double tauv(void) const;

  //! Heat transfer relaxation time.
  double tauT(void) const;

  //! Slip Reynolds number.
  double Rep(const Vector2D &u) const;

  //! Stokes drag correction.
  double fRep(const Vector2D &u) const;

  //! Particle Knudsen number.
  double Kn(void) const;

  //! Cunningham correction.
  double Cc(void) const;

  //! Nusselt number.
  double Nu(const Vector2D &u) const;

  //! Particle loading factor.
  double chi(void) const;

  //! Particle-phase volume fraction.
  double zetap(void) const;

  //! Particle-phase mass fraction.
  double phip(void) const;

  //! Particle-phase multi-velocity component switch.
  void MultiVelocity_Switch(void);
  //@}

  //@{ @name Turbulence related functions.
  //! Return to the turbulent kinetic energy.
  double dk(void) const;

  //! Return to the turbulent specific dissipation.
  double domega(void) const;

  //! Return the turbulent eddy dissipation.
  double epsilon(void) const;

  //! Return the turbulent eddy dissipation.
  double depsilon(void) const;

  //! Return the turbulent length scale.
  double ell(void) const;

  //! Return the tubulent Mach number.
  double Mt(void) const;

  //! Return the tubulent Mach number squared.
  double Mt2(void) const;

  //! Turbulent eddy dynamic viscosity.
  double muT(void) const;

  //! Turbulent eddy kinematic viscosity.
  double nuT(void) const;

  //! Turbulent eddy thermal heat conductivity.
  double kappaT(void) const;

  //! Turbulence modified sound speed.
  double c(void) const;

  //! Turbulence modified sound speed squared.
  double c2(void) const;

  //! Turbulence modified pressure.
  double pmodified(void) const;

  double beta_k(const Dusty2D_pState &dWdx,const Dusty2D_pState &dWdy) const;
  double beta_omega(const Dusty2D_pState &dWdx, const Dusty2D_pState &dWdy) const;
  double f_beta_k(const Dusty2D_pState &dWdx, const Dusty2D_pState &dWdy) const;
  double f_beta_omega(const Dusty2D_pState &dWdx, const Dusty2D_pState &dWdy) const;
  double chi_k(const Dusty2D_pState &dWdx, const Dusty2D_pState &dWdy) const;
  double chi_omega(const Dusty2D_pState &dWdx, const Dusty2D_pState &dWdy) const;
  //@}

  //@{ @name Solid propellant related functions.
  //! Burning rate.
  double burningrate(void) const;
  //@}

  //@{ @name Conserved solution state.
  Dusty2D_cState U(void) const;
  Dusty2D_cState U(const Dusty2D_pState &W) const;
  friend Dusty2D_cState U(const Dusty2D_cState &W);
  //@}

  //! Jacobian of the conserved solution variables with respect to the
  //! primitive solution variables.
  void dUdW(DenseMatrix &dUdW) const;

  //! Jacobian of the primitive solution variables with respect to the
  //! conserved solution variables.
  void dWdU(DenseMatrix &dWdU) const;

  //@{ @name Inviscid solution flux (x-direction) and Jacobian.
  Dusty2D_cState F(void) const;
  Dusty2D_cState F(const Vector2D &V) const;
  void dFdU(DenseMatrix &dFdU) const;
  //@}

  //@{ @name Viscous solution fluxes and Jacobians.
  Dusty2D_cState Gx(const Dusty2D_pState &dWdx) const;
  Dusty2D_cState Gy(const Dusty2D_pState &dWdy) const;
  //Dusty2D_cState dGxdU(???) const;
  //Dusty2D_cState dGydU(???) const;
  //@}

  //! Compute viscous stress tensor and heat flux vector.
  void ComputeViscousTerms(const Dusty2D_pState &dWdx,
			   const Dusty2D_pState &dWdy,
			   const Vector2D &X,
			   const int &Axisymmetric,
			   const int &adiabatic_flag);

  //@{ @name Eigenstructure.
  //! Eigenvalue(s) (x-direction).
  Dusty2D_pState lambda_x(void) const;

  //! Eigenvalue(s) (x-direction).
  Dusty2D_pState lambda_x(const Vector2D &V) const;

  //! Primitive right eigenvector (x-direction).
  Dusty2D_pState rp_x(int index) const;

  //! Conserved right eigenvector (x-direction).
  Dusty2D_cState rc_x(int index) const;

  //! Primitive left eigenvector (x-direction).
  Dusty2D_pState lp_x(int index) const;
  //@}

  //@{ @name Include all source vectors and Jacobians.
  Dusty2D_cState S(const Vector2D &X,
		   const Dusty2D_pState &dWdx,
		   const Dusty2D_pState &dWdy,
		   const int &Axisymmetric) const;
  void dSdU(DenseMatrix &dSdU,
	    const Vector2D &X,
	    const Dusty2D_pState &dWdx,
	    const Dusty2D_pState &dWdy,
	    const int &Axisymmetric) const;
  //@}

  //@{ @name Inviscid axisymmetric flow source vector and Jacobian.
  Dusty2D_cState Si(const Vector2D &X) const;
  void dSidU(DenseMatrix &dSidU, const Vector2D &X) const;
  //@}

  //@{ @name Viscous axisymmetric flow source vector and Jacobian.
  Dusty2D_cState Sv(const Vector2D &X,
		    const Dusty2D_pState &dWdy) const;
  void dSvdU(DenseMatrix &dSvdU, const Dusty2D_pState &dWdy, const Vector2D &X) const;
  //@}

  //@{ @name Turbulent source term vector and Jacobian.
  Dusty2D_cState St(const Vector2D &X,
		    const Dusty2D_pState &dWdx,
		    const Dusty2D_pState &dWdy,
		    const int &Axisymmetric) const;
  void dStdU(DenseMatrix &dStdU,
	     const Vector2D &X,
	     const Dusty2D_pState &dWdx,
	     const Dusty2D_pState &dWdy,
	     const int &Axisymmetric) const;
  //@}

  //@{ @name Phase interaction source term vector and Jacobian.
  Dusty2D_cState Sp(const int &interaction_flag) const;
  void dSpdU(DenseMatrix &dSpdU, const int &interaction_flag) const;
  //@}

  //@{ @name Electrostatic force source term vector and Jacobian.
  Dusty2D_cState Se(const Electrostatic2DState &We) const;
  void dSedU(DenseMatrix &dSpdU, const Electrostatic2DState &We) const;
  //@}

  //@{ @name Index operator.
  double &operator[](int index) {
    //assert(index >= 1 && index <= NUM_VAR_DUSTY2D);
    if (index <= NUM_VAR_BASE) {
      switch(index) {
      case 1 :
	return rho;
      case 2 :
	return v.x;
      case 3 :
	return v.y;
      case 4 :
	return p;
      case 5 :
	return k;
      case 6 :
	return omega;
      };
    } else if (index <= NUM_VAR_BASE + NUM_VAR_PART) {
      return Wp[cindex[index]][pindex[index]];
    }
    // Default return, this is never reached.
    return rho;
  }
  
  const double &operator[](int index) const {
    //assert(index >= 1 && index <= NUM_VAR_DUSTY2D);
    if (index <= NUM_VAR_BASE) {
      switch(index) {
      case 1 :
	return rho;
      case 2 :
	return v.x;
      case 3 :
	return v.y;
      case 4 :
	return p;
      case 5 :
	return k;
      case 6 :
	return omega;
      };
    } else if (index <= NUM_VAR_BASE + NUM_VAR_PART) {
      return Wp[cindex[index]][pindex[index]];
    }
    // Default return, this is never reached.
    return rho;
  }
  //@}

  //@{ @name Binary arithmetic operators.
  Dusty2D_pState operator +(const Dusty2D_pState &W) const;
  Dusty2D_pState operator -(const Dusty2D_pState &W) const;
  double operator *(const Dusty2D_pState &W1) const;
  Dusty2D_pState operator *(const double &a) const;
  friend Dusty2D_pState operator *(const double &a, const Dusty2D_pState &W);
  Dusty2D_pState operator /(const double &a) const;
  Dusty2D_pState operator ^(const Dusty2D_pState &W) const;
  //@}

  //@{ @name Assignment operator.
  Dusty2D_pState& operator =(const Dusty2D_pState &W);
  //@}

  //@{ @name Unary arithmetic operators.
  //Dusty2D_pState operator +(const Dusty2D_pState &W);
  friend Dusty2D_pState operator -(const Dusty2D_pState &W);
  //@}

  //@{ @name Shortcut arithmetic operators.
  Dusty2D_pState &operator +=(const Dusty2D_pState &W);
  Dusty2D_pState &operator -=(const Dusty2D_pState &W);
  Dusty2D_pState &operator *=(const double &a);
  Dusty2D_pState &operator /=(const double &a);
  //@}

  //@{ @name Relational operators.
  friend int operator ==(const Dusty2D_pState &W1, const Dusty2D_pState &W2);
  friend int operator !=(const Dusty2D_pState &W1, const Dusty2D_pState &W2);
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const Dusty2D_pState &W);
  friend istream &operator >> (istream &in_file, Dusty2D_pState &W);
  //@}

  void output_labels(ostream &out_file) {
    out_file << "\"rho\" \\ \n"
	     << "\"vx\" \\ \n"
	     << "\"vy\" \\ \n"
	     << "\"p\" \\ \n"
	     << "\"T\" \\ \n"
	     << "\"M\" \\ \n"
	     << "\"H\" \\ \n"
	     << "\"s\" \\ \n";
    if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      out_file << "\"k\" \\ \n"
	       << "\"omega\" \\ \n"
	       << "\"epsilon\" \\ \n"
	       << "\"ell\" \\ \n"
	       << "\"p_modified\" \\ \n";
    }
    if (NUM_VAR_PART) {
      out_file << "\"sigma\" \\ \n"
	       << "\"ux\" \\ \n"
	       << "\"uy\" \\ \n"
	       << "\"Tp\" \\ \n"
	       << "\"chi\" \\ \n"
	       << "\"zetap\" \\ \n"
	       << "\"phip\" \\ \n";
      out_file << "\"ln(sigmap)\" \\ \n";
      if (NUM_CMP_PART == PARTICLE2D_MULTI_VELOCITY_FORMULATION) {
	for (int npc = 0; npc < NUM_CMP_PART; npc++) {
	  out_file << "\"sigma" << npc+1 << "\" \\ \n"
		   << "\"ux" << npc+1 << "\" \\ \n"
		   << "\"uy" << npc+1 << "\" \\ \n"
		   << "\"Tp" << npc+1 << "\" \\ \n";
	}
      }
    }
  }

  void output_data(ostream &out_file, double dummy1, double dummy2) { //dummies needed for compatibility with NS turbulent
    out_file << " " << rho << " " << v.x << " " << v.y << " " << p
	     << " " << T() << " " << M() << " " << H() << " " << s();
    if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      out_file << " " << k << " " << omega << " " << epsilon()
	       << " " << ell() << " " << pmodified();
    }
    if (NUM_VAR_PART) {
      if (Wp.sigma() > TOLER) {
	out_file << " " << Wp.sigma()
		 << Wp.u()
		 << " " << Wp.Tp()
		 << " " << chi()
		 << " " << zetap()
		 << " " << phip();
      } else {
	out_file << " " << ZERO << " " << ZERO << " " << ZERO << " " << ZERO
		 << " " << ZERO << " " << ZERO << " " << ZERO;
      }
      out_file << " " << log(max(Wp.sigma(),TOLER));
      if (NUM_CMP_PART == PARTICLE2D_MULTI_VELOCITY_FORMULATION) {
	for (int npc = 0; npc < NUM_CMP_PART; npc++) {
	  out_file << Wp[npc];
	}
      }
    }
  }

  int analytically_inverted_relaxation() { //this is needed for embeddedboundaries with gaussian2D
    return 0;
  }
  void relax(double deltat, int stage, const Dusty2D_pState &W) {return;} //this is needed for embeddedboundaries with gaussian2D

  double pressure() const {return p;} //added for compatibility with embeddedboundaries2D

};

/*!
 * Class: Dusty2D_cState
 *
 * @brief Conserved variable solution state class definition for a 
 * gas-particle flow.
 *
 * Conserved variable solution state class definition for a gas-particle
 * flow.  The particle-phase is modelled using an Eulerian formulation
 * and is asumed to be inert, disperse and dilute particle-phase.  The
 * gas-phase can be inviscid, laminar, or turbulent.  The two-equation 
 * k-omega turbulence model is used to achieve turbulence closure.
 * Electrostatic effects on the particle-phase can also be included.
 *
 * \verbatim
 * Member functions
 *     rho      -- Gas-phase density.
 *     dv       -- Gas-phase momentum.
 *     E        -- Gas-phase total energy.
 *     dk       -- Gas-phase total turbulent kinetic energy.
 *     domega   -- Gas-phase total specific dissipation rate.
 *     tau      -- Viscous stress tensor (laminar and Reynolds).
 *     q        -- Heat flux vector (laminar and turbulent).
 *     g        -- Gas-phase specific heat ratio.
 *     gm1      -- Return g-1.
 *     gm1i     -- Return 1/(g-1).
 *     R        -- Gas gas constant.
 *     cp       -- Specific heat at constant pressure.
 *     v1,v2,v3,v4,v5 -- Viscosity law coefficients.
 *     v        -- Gas-phase flow velocity.
 *     p        -- Gas-phase pressure.
 *     T        -- Gas-phase temperature.
 *     e        -- Gas-phase specific internal energy.
 *     h        -- Gas-phase specific enthalpy.
 *     H        -- Gas-phase total enthalpy.
 *     a        -- Gas-phase sound speed.
 *     a2       -- Gas-phase sound speed squared.
 *     M        -- Gas-phase Mach number.
 *     s        -- Gas-phase specific entropy.
 *     To       -- Gas-phase stagnation temperature.
 *     po       -- Gas-phase stagnation pressure.
 *     ao       -- Gas-phase stagnation sound speed.
 *     ho       -- Gas-phase stagnation enthalpy.
 *     mu       -- Gas-phase dynamic viscosity.
 *     nu       -- Gas-phase kinematic viscosity.
 *     kappa    -- Gas-phase thermal conductivity.
 *     Pr       -- Prandtl number.
 *     meanfreepath -- Gas-phase mean free path.
 *
 *     Up       -- Array of particle conservative states.
 *     NUM_CMP_PART -- Number of particle components.
 *     NUM_VAR_PART -- Total number of particle variables.
 *
 *     dk       -- Return or write to the turbulent kinetic energy.
 *     domega   -- Return or write to the turbulent specific
 *                 dissipation.
 *     depsilon -- Return or write to the turbulent eddy dissipation.
 *     muT      -- Turbulent eddy dynamic viscosity.
 *     nuT      -- Turbulent eddy kinematic viscosity.
 *     kappaT   -- Turbulent eddy thermal conductivity.
 *     PrT      -- Turbulent Prandtl number.
 *     c        -- Turbulence modified sound speed.
 *     c2       -- Turbulence modified sound speed squared.
 *     pmodified -- Turbulence modified pressure.
 *
 *     rhos     -- Propellant density.
 *     n        -- Propellant burning rate constant.
 *     beta     -- Propellant burning rate coefficient.
 *     Tf       -- Propellant flame temperature.
 *     Ts       -- Propellant surface temperature.
 *     alphas   -- Propellant particle concentration.
 *
 *     NUM_VAR_DUSTY2D -- Total number of variables.
 *
 *     tauv     -- Momentum transfer relaxation time scale.
 *     tauT     -- Heat transfer relaxation time scale.
 *     Rep      -- Slip Reynolds number.
 *     fRep     -- Stokes drag correction.
 *     Kn       -- Particle Knudsen number.
 *     Cc       -- Cunningham correction.
 *     Nu       -- Nusselt number.
 *     burningrate -- Burning rate.
 *
 *     W        -- Return primitive solution state.
 *     dUdW     -- Return the Jacobian of the conserved solution
 *                 variables with respect to the primitive solution
 *                 variables.
 *     dWdU     -- Return the Jacobian of the primitive solution
 *                 variables with respect to the conserved solution
 *                 variables.
 *
 *     F        -- Return x-direction solution flux.
 *     dFdU     -- Return the Jacobian of the inviscid solution flux
 *                 vector with respect to the conserved solution
 *                 variables.
 *
 *     Gx       -- Return x-direction viscous solution flux.
 *     Gy       -- Return y-direction viscous solution flux.
 *     dGxdU    -- Return the Jacobian of the x-direction viscous
 *                 solution flux vector with respect to the conserved
 *                 solution variables.
 *     dGydU    -- Return the Jacobian of the y-direction viscous
 *                 solution flux vector with respect to the conserved
 *                 solution variables.
 *
 *     lambda_x -- Return x-direction eigenvalue(s).
 *     rp_x     -- Return primitive right eigenvector (x-direction).
 *     rc_x     -- Return conserved right eigenvector (x-direction).
 *     lp_x     -- Return primitive left eigenvector (x-direction).
 *
 * Member operators
 *      U -- a primitive solution state
 *      c -- a scalar (double)
 *
 * U = U;
 * c = U[i];
 * U = U + U;
 * U = U - U;
 * c = U * U; (inner product)
 * U = c * U;
 * U = U * c;
 * U = U / c;
 * U = U ^ U; (a useful product)
 * U = +U;
 * U = -U;
 * U += U;
 * U -= U;
 * U == U;
 * U != U;
 * cout << U; (output function)
 * cin  >> U; (input function)
 * \endverbatim
 */
class Dusty2D_cState {
  private:
  public:
  //@{ @name Gas-phase conservative variables:
  double                       rho; //!< Gas-phase density.
  Vector2D                      dv; //!< Gas-phase momentum.
  double                         E; //!< Gas-phase total energy.
  double                        dk; //!< Gas-phase total turbulent kinetic energy.
  double                    domega; //!< Gas-phase total turbulent specific dissipation rate.
  Tensor2D                     tau; //!< Viscous stess tensor (laminar and turbulent).
  Vector2D                       q; //!< Heat flux vector (laminar and turbulent).
  static double                  g; //!< Specific heat ratio.
  static double                gm1; //!< g-1
  static double               gm1i; //!< 1/(g-1)
  static double                  R; //!< Gas constant.
  static double                 cp; //!< Specific heat at constant pressure.
  static double v1, v2, v3, v4, v5; //!< Viscosity law coefficients.
  static int             flow_type; //!< Flow-type indicator (inviscid, laminar, or k-omega turbulent).
  //@}

  //@{ @name Particle-phase conservative variables:
  Particle2D_cComponents        Up; //!< Array of particle conserved components.
  static int          NUM_CMP_PART; //!< Number of particle components.
  static int          NUM_VAR_PART; //!< Number of particle variables.
  static int              drag_law; //!< Drag law definition.
  static double d1, d2, d3, d4, d5; //!< Drag law coefficients.
  static double                 cm; //!< Characteristic particle specific heat.
  static double                 dp; //!< Characteristic particle diameter.
  static double               rhop; //!< Characteristic particle density.
  static double                 mp; //!< Characteristic particle mass.
  static double                 qe; //!< Characteristic particle charge.
  //@}

  //@{ @name Turbulent boundary-layer constants:
  static double            yplus_o; //!< Transition between viscous sublayer and log layer.
  static double                  C; //!< Surface roughness coefficient.
  static double         von_karman; //!< Von Karman constant.
  static double     yplus_sublayer; //!< Sublayer dimensionless wall distance.
  static double yplus_buffer_layer; //!< Buffer layer dimensionless wall distance.
  static double  yplus_outer_layer; //!< Outer layer dimensionless wall distance.
  //@}

  //@{ @name k-omega closure coefficients:
  static double                PrT; //!< Turbulent Prandtl number.
  static double                Cmu;
  static double           beta_k_o; //!< Destruction of turbulent kinetic energy closure coefficient (a.k.a beta_star).
  static double       beta_omega_o; //!< Destruction of specific dissipation closure coefficient (a.k.a beta).
  static double            sigma_k;
  static double        sigma_omega;
  static double              alpha;
  static double                 xi; //!< Compressiblity correction coefficient.
  static double                Mto; //!< Compressiblity correction coefficient.
  //@}

  //@{ @name Propellant variables:
  static double               rhos; //!< Propellant density.
  static double               beta; //!< Propellant burning rate constant.
  static double                  n; //!< Propellant burning rate coefficient.
  static double                 Tf; //!< Propellant flame temperature.
  static double                 Ts; //!< Propellant surface temperature.
  static double             alphas; //!< Propellant particle concentration.
  //@}

  //@{ @name General purpose static variables:
  static int       NUM_VAR_DUSTY2D; //!< Total number of variables.
  static int               *cindex; //!< Look up table for the particle-phase component indicing.
  static int               *pindex; //!< Look up table for the particle-phase state indicing.
  //@}

  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  Dusty2D_cState(void) {
    allocate();
    Standard_Atmosphere();
  }

  //! Copy constructor.
  Dusty2D_cState(const Dusty2D_cState &U) {
    allocate();
    rho = U.rho; dv.x = U.dv.x; dv.y = U.dv.y; E = U.E; dk = U.dk; domega = U.domega;
    if (NUM_VAR_PART) Up = U.Up;
  }

  //! Copy constructor.
  Dusty2D_cState(const Dusty2D_pState &W);

  //! Assignment constructor.
  Dusty2D_cState(const double &dens,
		 const Vector2D &dV,
		 const double &Etot) {
    allocate();
    rho = dens; dv.x = dV.x; dv.y = dV.y; E = Etot;
  }

  //! Assignment constructor.
  Dusty2D_cState(const double &dens,
		 const double &dvx,
		 const double &dvy,
		 const double &Etot) {
    allocate();
    rho = dens; dv.x = dvx; dv.y = dvy; E = Etot;
  }

  //! Assignment constructor.
  Dusty2D_cState(const double &dens,
		 const Vector2D &dV,
		 const double &Etot,
		 const Particle2D_cComponents &Ucomp) {
    allocate();
    rho = dens; dv.x = dV.x; dv.y = dV.y; E = Etot;
    Up = Ucomp;
  }

  //! Assignment constructor.
  Dusty2D_cState(const double &dens,
		 const double &dvx,
		 const double &dvy,
		 const double &Etot,
		 const Particle2D_cComponents &Ucomp) {
    allocate();
    rho = dens; dv.x = dvx; dv.y = dvy; E = Etot;
    Up = Ucomp;
  }

  //! Assignment constructor.
  Dusty2D_cState(const double &dens,
		 const Vector2D &dV,
		 const double &Etot,
		 const double &dkdk,
		 const double &domga) {
    allocate();
    rho = dens; dv.x = dV.x; dv.y = dV.y; E = Etot; dk = dkdk; domega = domga;
  }

  //! Assignment constructor.
  Dusty2D_cState(const double &dens,
		 const double &dvx,
		 const double &dvy,
		 const double &Etot,
		 const double &dkdk,
		 const double &domga) {
    allocate();
    rho = dens; dv.x = dvx; dv.y = dvy; E = Etot; dk = dkdk; domega = domga;
  }

  //! Assignment constructor.
  Dusty2D_cState(const double &dens,
		 const Vector2D &dV,
		 const double &Etot,
		 const Particle2D_cComponents &Ucomp,
		 const double &dkdk,
		 const double &domga) {
    allocate();
    rho = dens; dv.x = dV.x; dv.y = dV.y; E = Etot; dk = dkdk; domega = domga;
    Up = Ucomp;
  }

  //! Assignment constructor.
  Dusty2D_cState(const double &dens,
		 const double &dvx,
		 const double &dvy,
		 const double &Etot,
		 const Particle2D_cComponents &Ucomp,
		 const double &dkdk,
		 const double &domga) {
    allocate();
    rho = dens; dv.x = dvx; dv.y = dvy; E = Etot; dk = dkdk; domega = domga;
    Up = Ucomp;
  }

  //! Destructor.
  ~Dusty2D_cState(void) { deallocate(); }
  //@}
  
  //! Return the number of variables.
  int NumVar(void) const { return NUM_VAR_DUSTY2D; }

  //@{ @name Allocation and deallocation functions.
  //! Allocation function.
  void allocate(void);
  //! Deallocation functions.
  void deallocate(void);
  //@}

  //@{ @name Set static variables.
  void set_static_variables(void);
  void set_static_variables(char *gas_type,
			    const int &FlowType,
			    const double &C_constant,
			    const double &von_karman_constant,
			    const double &yplus_sub,
			    const double &yplus_buffer,
			    const double &yplus_outer,
			    char *propellant_type);
  void set_static_variables(char *gas_type,
			    char *particle_type,
			    const int &num_cmp_part,
			    const int &draglaw,
			    const int &FlowType,
			    const double &C_constant,
			    const double &von_karman_constant,
			    const double &yplus_sub,
			    const double &yplus_buffer,
			    const double &yplus_outer,
			    char *propellant_type);
  void set_gas(char *gas_type);
  void set_particles(char *particle_type, const int &draglaw);
  void set_turbulence(const double &C_constant,
		      const double &von_karman_constant,
		      const double &yplus_sub,
		      const double &yplus_buffer,
		      const double &yplus_outer);
  void set_propellant(char *propellant_type);
  void set_indicing(void);
  //@}

  //@{ @name Useful operators.
  //! Copy operator.
  void Copy(const Dusty2D_cState &U);

  //! Vacuum operator.
  void Vacuum(void);

  //! Standard atmosphere operator.
  void Standard_Atmosphere(void);

  //! Check for unphysical state properties.
  int Unphysical_Properties(void) const;

  //! Copy variables solved by multigrid only.
  void Copy_Multigrid_State_Variables(const Dusty2D_cState &Ufine);

  //! Zero variables not-solved by multigrid.
  void Zero_Non_Multigrid_State_Variables(void);
  //@}

  //@{ @name Gas-phase related functions.

  //! Gas-phase flow velocity.
  Vector2D v(void) const;

  //! Gas-phase flow velocity.
  double v(const Vector2D &n) const;

  //! Gas-phase pressure.
  double p(void) const;

  //! Gas-phase temperature.
  double T(void) const;

  //! Gas-phase specific internal energy.
  double e(void) const;

  //! Gas-phase specific enthalpy.
  double h(void) const;

  //! Gas-phase total enthalpy.
  double H(void) const;

  //! Gas-phase sound speed.
  double a(void) const;

  //! Gas-phase sound speed squared.
  double a2(void) const;

  //! Gas-phase Mach number.
  double M(void) const;

  //! Gas-phase specific entropy.
  double s(void) const;

  //! Gas-phase stagnation temperature.
  double To(void) const;

  //! Gas-phase stagnation pressure.
  double po(void) const;

  //! Gas-phase stagnation sound speed.
  double ao(void) const;

  //! Gas-phase stagnation enthalpy.
  double ho(void) const;

  //! Gas-phase dynamic viscosity.
  double mu(void) const;

  //! Gas-phase kinematic viscosity.
  double nu(void) const;

  //! Gas-phase thermal heat conductivity.
  double kappa(void) const;

  //! Prandtl number.
  double Pr(void) const;

  //! Gas-phase mean free path.
  double meanfreepath(void) const;
  //@}

  //@{ @name Phase interaction functions.
  //! Momentum transfer relaxation time.
  double tauv(void) const;

  //! Heat transfer relaxation time.
  double tauT(void) const;

  //! Slip Reynolds number.
  double Rep(const Vector2D &u) const;

  //! Stokes drag correction.
  double fRep(const Vector2D &u) const;

  //! Particle Knudsen number.
  double Kn(void) const;

  //! Cunningham correction.
  double Cc(void) const;

  //! Nusselt number.
  double Nu(const Vector2D &u) const;

  //! Particle loading factor.
  double chi(void) const;

  //! Particle-phase volume fraction.
  double zetap(void) const;

  //! Particle-phase mass fraction.
  double phip(void) const;
  //@}

  //@{ @name Turbulence (k-omega) related functions.
  //! Return the turbulent kinetic energy.
  double k(void) const;

  //! Return the turbulent specific dissipation.
  double omega(void) const;

  //! Return the turbulent eddy dissipation.
  double depsilon(void) const;

  //! Return the turbulent eddy dissipation.
  double epsilon(void) const;

  //! Return the turbulent length scale.
  double ell(void) const;

  //! Return the tubulent Mach number.
  double Mt(void) const;

  //! Return the tubulent Mach number squared.
  double Mt2(void) const;

  //! Turbulent eddy dynamic viscosity.
  double muT(void) const;

  //! Turbulent eddy kinematic viscosity.
  double nuT(void) const;

  //! Turbulent eddy thermal heat conductivity.
  double kappaT(void) const;

  //! Turbulence modified sound speed.
  double c(void) const;

  //! Turbulence modified sound speed squared.
  double c2(void) const;

  //! Turbulence modified pressure.
  double pmodified(void) const;

  double beta_k(const Dusty2D_pState &dWdx, const Dusty2D_pState &dWdy) const;
  double beta_omega(const Dusty2D_pState &dWdx, const Dusty2D_pState &dWdy) const;
  double f_beta_k(const Dusty2D_pState &dWdx, const Dusty2D_pState &dWdy) const;
  double f_beta_omega(const Dusty2D_pState &dWdx, const Dusty2D_pState &dWdy) const;
  double chi_k(const Dusty2D_pState &dWdx, const Dusty2D_pState &dWdy) const;
  double chi_omega(const Dusty2D_pState &dWdx, const Dusty2D_pState &dWdy) const;
  //@}

  //@{ @name Solid propellant related functions.
  //! Burning rate.
  double burningrate(void) const;
  //@}

  //@{ @name Primitive solution state.
  Dusty2D_pState W(void) const;
  Dusty2D_pState W(const Dusty2D_cState &U) const;
  friend Dusty2D_pState W(const Dusty2D_cState &U);
  //@}

  //! Jacobian of the conserved solution variables with respect to the
  //! primitive solution variables.
  void dUdW(DenseMatrix &dUdW) const;

  //! Jacobian of the primitive solution variables with respect to the
  //! conserved solution variables.
  void dWdU(DenseMatrix &dWdU) const;

  //@{ @name Inviscid solution flux (x-direction) and Jacobian.
  Dusty2D_cState F(void) const;
  Dusty2D_cState F(const Vector2D &V) const;
  void dFdU(DenseMatrix &dFdU) const;
  //@}

  //@{ @name Viscous solution flux and Jacobians.
  Dusty2D_cState Gx(const Dusty2D_pState &dWdx) const;
  Dusty2D_cState Gy(const Dusty2D_pState &dWdy) const;
  //Dusty2D_cState dGxdU(???) const;
  //Dusty2D_cState dGydU(???) const;
  //@}

  //! Compute viscous stress tensor and heat flux vector.
  void ComputeViscousTerms(const Dusty2D_pState &dWdx,
			   const Dusty2D_pState &dWdy,
			   const Vector2D &X,
			   const int &Axisymmetric,
			   const int &adiabatic_flag);

  //@{ @name Eigenstructure.
  //! Eigenvalue(s) (x-direction).
  Dusty2D_pState lambda_x(void) const;
  Dusty2D_pState lambda_x(const Vector2D &V) const;

  //! Primitive right eigenvector (x-direction).
  Dusty2D_pState rp_x(int index) const;

  //! Conserved right eigenvector (x-direction).
  Dusty2D_cState rc_x(int index) const;

  //! Primitive left eigenvector (x-direction).
  Dusty2D_pState lp_x(int index) const;
  //@}

  //@{ @name Index operator.
  double &operator[](int index) {
    //assert(index >= 1 && index <= NUM_VAR_DUSTY2D);
    if (index <= NUM_VAR_BASE) {
      switch(index) {
      case 1 :
	return rho;
      case 2 :
	return dv.x;
      case 3 :
	return dv.y;
      case 4 :
	return E;
      case 5 :
	return dk;
      case 6 :
	return domega;
      };
    } else if (index <= NUM_VAR_BASE + NUM_VAR_PART) {
      return Up[cindex[index]][pindex[index]];
    }
    // Default return, this is never reached.
    return rho;
  }

  const double &operator[](int index) const {
    //assert(index >= 1 && index <= NUM_VAR_DUSTY2D);
    if (index <= NUM_VAR_BASE) {
      switch(index) {
      case 1 :
	return rho;
      case 2 :
	return dv.x;
      case 3 :
	return dv.y;
      case 4 :
	return E;
      case 5 :
	return dk;
      case 6 :
	return domega;
      };
    } else if (index <= NUM_VAR_BASE + NUM_VAR_PART) {
      return Up[cindex[index]][pindex[index]];
    }
    // Default return, this is never reached.
    return rho;
  }
  //@}

  //@{ @name Binary arithmetic operators.
  Dusty2D_cState operator +(const Dusty2D_cState &U) const;
  Dusty2D_cState operator -(const Dusty2D_cState &U) const;
  double operator *(const Dusty2D_cState &U) const;
  Dusty2D_cState operator *(const double &a) const;
  friend Dusty2D_cState operator *(const double &a, const Dusty2D_cState &U);
  Dusty2D_cState operator /(const double &a) const;
  Dusty2D_cState operator ^(const Dusty2D_cState &U) const;
  //@}

  //@{ @name Assignment operator.
  Dusty2D_cState& operator =(const Dusty2D_cState &U);
  //@}

  //@{ @name Unary arithmetic operators.
  //Dusty2D_cState operator +(const Dusty2D_cState &U);
  friend Dusty2D_cState operator -(const Dusty2D_cState &U);
  //@}

  //@{ @name Shortcut arithmetic operators.
  Dusty2D_cState &operator +=(const Dusty2D_cState &U);
  Dusty2D_cState &operator -=(const Dusty2D_cState &U);
  Dusty2D_cState &operator *=(const double &a);
  Dusty2D_cState &operator /=(const double &a);
  //@}

  //@{ @name Relational operators.
  friend int operator ==(const Dusty2D_cState &U1, const Dusty2D_cState &U2);
  friend int operator !=(const Dusty2D_cState &U1, const Dusty2D_cState &U2);
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const Dusty2D_cState &U);
  friend istream &operator >> (istream &in_file, Dusty2D_cState &U);
  //@}

};

/**********************************************************************
 * Dusty2D_pState::allocate -- Allocate memory for data arrays.       *
 **********************************************************************/
inline void Dusty2D_pState::allocate(void) {
  if (NUM_VAR_PART) Wp.allocate();
}

/**********************************************************************
 *  Dusty2D_pState::deallocate -- Deallocate memory.                  *
 **********************************************************************/
inline void Dusty2D_pState::deallocate(void) {
  if (NUM_VAR_PART) Wp.deallocate();
}

/**********************************************************************
 * Dusty2D_pState::Copy -- Copy operator.                             *
 **********************************************************************/
inline void Dusty2D_pState::Copy(const Dusty2D_pState &W) {
  rho = W.rho; v.x = W.v.x; v.y = W.v.y; p = W.p; k = W.k; omega = W.omega;
  if (NUM_VAR_PART) Wp.Copy(W.Wp);
}

/**********************************************************************
 * Dusty2D_pState::Vacuum -- Vacuum operator.                         *
 **********************************************************************/
inline void Dusty2D_pState::Vacuum(void) {
  rho = ZERO; v.x = ZERO; v.y = ZERO; p = ZERO; k = ZERO; omega = ZERO; tau.zero(); q.zero();
  if (NUM_VAR_PART) Wp.Vacuum();
}

/**********************************************************************
 * Dusty2D_pState::Standard_Atmosphere -- Standard atmosphere operator.*
 **********************************************************************/
inline void Dusty2D_pState::Standard_Atmosphere(void) {
  rho = DENSITY_STDATM; v.x = ZERO; v.y = ZERO; p = PRESSURE_STDATM; k = ZERO; omega = ZERO; tau.zero(); q.zero();
  if (NUM_VAR_PART) Wp.Vacuum();
}

/**********************************************************************
 * Dusty2D_pState::Unphysical_Properties -- Check for unphysical      *
 *                                          state properties.         *
 **********************************************************************/
inline int Dusty2D_pState::Unphysical_Properties(void) const {
  if (rho <= ZERO || p <= ZERO || E() <= ZERO) return 1;
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) if (k < ZERO || omega < ZERO) return 1;
  if (NUM_VAR_PART) return Wp.Unphysical_Properties();
  return 0;
}

/**********************************************************************
 * Dusty2D_pState::set_static_variables -- Set all static variables.  *
 **********************************************************************/
inline void Dusty2D_pState::set_static_variables(void) {
  // Set gas constants.
  set_gas("AIR");
  // Set particle constants.
  NUM_CMP_PART = 0;
  NUM_VAR_PART = 0;
  // Set total number of variables.
  NUM_VAR_DUSTY2D = NUM_VAR_BASE + NUM_VAR_PART;
  // Allocate data arrays.
  allocate();
  // Set indicing look up tables.
  set_indicing();
  // Set the flow type.
  flow_type = FLOWTYPE_INVISCID;
  // Set turbulence constants.
  set_turbulence(ZERO,ZERO,ZERO,ZERO,ZERO);
  // Set propellant type.
  set_propellant("AP_HTPB");
}

inline void Dusty2D_pState::set_static_variables(char *gas_type,
						 const int &FlowType,
						 const double &C_constant,
						 const double &von_karman_constant,
						 const double &yplus_sub,
						 const double &yplus_buffer,
						 const double &yplus_outer,
						 char *propellant_type) {
  // Deallocate memory arrays.
  deallocate();
  // Set gas constants.
  set_gas(gas_type);
  // Set particle constants.
  NUM_CMP_PART = 0;
  NUM_VAR_PART = 0;
  // Set total number of variables.
  NUM_VAR_DUSTY2D = NUM_VAR_BASE + NUM_VAR_PART;
  // Allocate data arrays.
  allocate();
  // Set indicing look up tables.
  set_indicing();
  // Set the flow type.
  flow_type = FlowType;
  // Set turbulence constants.
  set_turbulence(C_constant,von_karman_constant,yplus_sub,yplus_buffer,yplus_outer);
  // Set propellant type.
  set_propellant(propellant_type);
}

inline void Dusty2D_pState::set_static_variables(char *gas_type,
						 char *particle_type,
						 const int &num_cmp_part,
						 const int &draglaw,
						 const int &FlowType,
						 const double &C_constant,
						 const double &von_karman_constant,
						 const double &yplus_sub,
						 const double &yplus_buffer,
						 const double &yplus_outer,
						 char *propellant_type) {
  // Deallocate memory arrays.
  deallocate();
  // Set gas constants.
  set_gas(gas_type);
  // Set particle constants.
  NUM_CMP_PART = num_cmp_part;
  NUM_VAR_PART = NUM_CMP_PART*NUM_VAR_PARTICLE2D;
  if (NUM_VAR_PART) {
    set_particles(particle_type,draglaw);
    Wp.set_particle_components(num_cmp_part);
  }
  // Set total number of variables.
  NUM_VAR_DUSTY2D = NUM_VAR_BASE + NUM_VAR_PART;
  // Allocate data arrays.
  allocate();
  // Set indicing look up tables.
  set_indicing();
  // Set the flow type.
  flow_type = FlowType;
  // Set turbulence constants.
  set_turbulence(C_constant,von_karman_constant,yplus_sub,yplus_buffer,yplus_outer);
  // Set propellant type.
  set_propellant(propellant_type);
}

/**********************************************************************
 * Dusty2D_pState::set_gas -- Set gas-phase static variables.         *
 **********************************************************************/
inline void Dusty2D_pState::set_gas(char *gas_type) {
  if (strcmp(gas_type,"AIR") == 0) {
    g  = GAMMA_AIR;
    R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    v1 = AIR_c1; v2 = AIR_c2; v3 = AIR_c3; v4 = AIR_c4; v5 = AIR_c5;
  } else if (strcmp(gas_type,"H2") == 0) {
    g  = GAMMA_H2;
    R  = R_UNIVERSAL/(MOLE_WT_H2*MILLI);
    v1 = H2_c1; v2 = H2_c2; v3 = H2_c3; v4 = H2_c4; v5 = H2_c5;
  } else if (strcmp(gas_type,"HE") == 0) {
    g  = GAMMA_HE;
    R  = R_UNIVERSAL/(MOLE_WT_HE*MILLI);
    v1 = HE_c1; v2 = HE_c2; v3 = HE_c3; v4 = HE_c4; v5 = HE_c5;
  } else if (strcmp(gas_type,"N2") == 0) {
    g  = GAMMA_N2;
    R  = R_UNIVERSAL/(MOLE_WT_N2*MILLI);
    v1 = N2_c1; v2 = N2_c2; v3 = N2_c3; v4 = N2_c4; v5 = N2_c5;
  } else if (strcmp(gas_type,"O2") == 0) {
    g  = GAMMA_O2;
    R  = R_UNIVERSAL/(MOLE_WT_O2*MILLI);
    v1 = O2_c1; v2 = O2_c2; v3 = O2_c3; v4 = O2_c4; v5 = O2_c5;
  } else if (strcmp(gas_type,"AP_HTPB") == 0) {
    g  = GAMMA_APHTPB;
    R  = R_UNIVERSAL/(MOLE_WT_APHTPB*MILLI);
    v1 = APHTPB_c1; v2 = APHTPB_c2; v3 = APHTPB_c3; v4 = APHTPB_c4; v5 = APHTPB_c5;
  } else {
    g  = GAMMA_AIR;
    R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    v1 = AIR_c1; v2 = AIR_c2; v3 = AIR_c3; AIR_c4;
  }
  gm1  = g - 1.0;
  gm1i = 1.0/gm1;
  cp = g*R*gm1i;
}

/**********************************************************************
 * Dusty2D_pState::set_particles -- Set the particle constants and    *
 *                                  the drag law coefficients.        *
 **********************************************************************/
inline void Dusty2D_pState::set_particles(char *particle_type, const int &draglaw) {
  // Set the particle type.
  if (strcmp(particle_type,"AP_HTPB") == 0) {
    cm   = CM_APHTPB;
    dp   = DP_APHTPB;
    rhop = RHOP_APHTPB;
    mp   = MP_APHTPB;
    qe   = ZERO;
  } else if (strcmp(particle_type,"Glass_Beads") == 0) {
    cm   = CM_GLASS_BEADS;
    dp   = DP_GLASS_BEADS;
    rhop = RHOP_GLASS_BEADS;
    mp   = MP_GLASS_BEADS;
    qe   = QE_GLASS_BEADS;
  } else if (strcmp(particle_type,"SCIEX_TINY") == 0) {
    cm   = CM_SCIEX_TINY;
    dp   = DP_SCIEX_TINY;
    rhop = RHOP_SCIEX_TINY;
    mp   = MP_SCIEX_TINY;
    qe   = QE_SCIEX_TINY;
  } else if (strcmp(particle_type,"SCIEX_SMALL") == 0) {
    cm   = CM_SCIEX_SMALL;
    dp   = DP_SCIEX_SMALL;
    rhop = RHOP_SCIEX_SMALL;
    mp   = MP_SCIEX_SMALL;
    qe   = QE_SCIEX_SMALL;
  } else if (strcmp(particle_type,"SCIEX_MEDIUM") == 0) {
    cm   = CM_SCIEX_MEDIUM;
    dp   = DP_SCIEX_MEDIUM;
    rhop = RHOP_SCIEX_MEDIUM;
    mp   = MP_SCIEX_MEDIUM;
    qe   = QE_SCIEX_MEDIUM;
  } else if (strcmp(particle_type,"SCIEX_LARGE") == 0) {
    cm   = CM_SCIEX_LARGE;
    dp   = DP_SCIEX_LARGE;
    rhop = RHOP_SCIEX_LARGE;
    mp   = MP_SCIEX_LARGE;
    qe   = QE_SCIEX_LARGE;
  } else if (strcmp(particle_type,"Heavy_AP_HTPB") == 0) {
    cm   = CM_HEAVY_APHTPB;
    dp   = DP_HEAVY_APHTPB;
    rhop = RHOP_HEAVY_APHTPB;
    mp   = MP_HEAVY_APHTPB;
    qe   = QE_HEAVY_APHTPB;
  } else {
    cm   = CM_APHTPB;
    dp   = DP_APHTPB;
    rhop = RHOP_APHTPB;
    mp   = MP_APHTPB;
    qe   = ZERO;
  }
  // Set the drag law.
  drag_law = draglaw;
  switch(drag_law) {
  case DRAG_LAW_STOKES :
    d1 = 0.0;
    d2 = 0.0;
    d3 = 0.0;
    d4 = 0.0;
    d5 = 0.0;
    break;
  case DRAG_LAW_OSEEN :
    d1 = 3.0/16.0;
    d2 = 1.0;
    d3 = 0.0;
    d4 = 0.0;
    d5 = 0.0;
    break;
  case DRAG_LAW_KLYACHKO_PUTNAM :
    d1 = 1.0/6.0;
    d2 = 2.0/3.0;
    d3 = 0.0;
    d4 = 0.0;
    d5 = 0.0;
    break;
  case DRAG_LAW_TURTON_LEVENSPIEL :
    d1 = 0.173;
    d2 = 0.657;
    d3 = 0.413;
    d4 = 16300.0;
    d5 = -1.09;
    break;
  case DRAG_LAW_HAIDER_LEVENSPIEL :
    d1 = 0.1806;
    d2 = 0.6459;
    d3 = 0.4251;
    d4 = 6880.95;
    d5 = -1.0;
    break;
  };
}

/**********************************************************************
 * Dusty2D_pState::set_turbulence -- Set the turbulence static        *
 *                                   variables.                       *
 **********************************************************************/
inline void Dusty2D_pState::set_turbulence(const double &C_constant,
					   const double &von_karman_constant,
					   const double &yplus_sub,
					   const double &yplus_buffer,
					   const double &yplus_outer) {
  // k-omega closure coefficients:
  Cmu = 0.09;
  beta_k_o = 0.09;
  beta_omega_o = 0.072;
  sigma_k = 0.50;
  sigma_omega = 0.50;
  alpha = 0.52;
  // Turbulent boundary-layer constants:
  C = C_constant;
  von_karman = von_karman_constant;
  yplus_sublayer = yplus_sub;
  yplus_buffer_layer = yplus_buffer;
  yplus_outer_layer = yplus_outer;
  //yplus_o = Iterative_Sub_To_Log_Layer_Transition_Point(C,von_karman);
  double f, df, E = exp(von_karman*C);
  // Set the initial guess.
  yplus_o = 10.0;
  // Iterate to determine the transition point.
  do {
    f = von_karman*yplus_o - log(E*yplus_o);
    df = von_karman - ONE/max(yplus_o,TOLER*TOLER);
    yplus_o -= f/df;
  } while(fabs(f) >= 0.00000001);
}

/**********************************************************************
 * Dusty2D_pState::set_propellant -- Set propellant static variables. *
 **********************************************************************/
inline void Dusty2D_pState::set_propellant(char *propellant_type) {
  if (strcmp(propellant_type,"AP_HTPB") == 0) {
    rhos   = RHOS_APHTPB;
    n      = N_APHTPB;
    beta   = BETA_APHTPB;
    Tf     = TF_APHTPB;
    Ts     = TS_APHTPB;
    alphas = ALPHAS_APHTPB;
  } else if (strcmp(propellant_type,"QUICK_AP_HTPB") == 0) {
    rhos   = RHOS_QUICK;
    n      = N_QUICK;
    beta   = BETA_QUICK;
    Tf     = TF_QUICK;
    Ts     = TS_QUICK;
    alphas = ALPHAS_QUICK;
  } else if (strcmp(propellant_type,"PLAID_AP_HTPB") == 0) {
    rhos   = RHOS_PLAID;
    n      = N_PLAID;
    beta   = BETA_PLAID;
    Tf     = TF_PLAID;
    Ts     = TS_PLAID;
    alphas = ALPHAS_PLAID;
  } else {
    rhos   = RHOS_APHTPB;
    n      = N_APHTPB;
    beta   = BETA_APHTPB;
    Tf     = TF_APHTPB;
    Ts     = TS_APHTPB;
    alphas = ALPHAS_APHTPB;
  }
}

/**********************************************************************
 * Dusty2D_pState::set_indicing -- Set indicing look up tables.       *
 **********************************************************************/
inline void Dusty2D_pState::set_indicing(void) {
  int nc, np;
  if (cindex != NULL) { delete []cindex; cindex = NULL; }
  if (pindex != NULL) { delete []pindex; pindex = NULL; }
  cindex = new int[NUM_VAR_DUSTY2D+1]; cindex[0] = 0;
  pindex = new int[NUM_VAR_DUSTY2D+1]; pindex[0] = 0;
  for (int index = 1; index <= NUM_VAR_DUSTY2D; index++) {
    if (index <= NUM_VAR_BASE) {
      // Gas-phase variables.
      cindex[index] = index;
      pindex[index] = index;
    } else if (index <= NUM_VAR_BASE+NUM_VAR_PART) {
      // Particle-phase variables.
      Wp.get_indices(index-NUM_VAR_BASE-1,nc,np);
      cindex[index] = nc;
      pindex[index] = np;
    }
  }
}

/**********************************************************************
 * Dusty2D_pState::T -- Gas-phase temperature.                        *
 **********************************************************************/
inline double Dusty2D_pState::T(void) const {
  //assert(rho > ZERO);
  return p/(rho*R);
}

/**********************************************************************
 * Dusty2D_pState::e -- Gas-phase specific internal energy.           *
 **********************************************************************/
inline double Dusty2D_pState::e(void) const {
  //assert(rho > ZERO);
  return p/(gm1*rho);
}

/**********************************************************************
 * Dusty2D_pState::E -- Gas-phase total energy.                       *
 **********************************************************************/
inline double Dusty2D_pState::E(void) const {
  return p*gm1i + HALF*rho*v.sqr() + dk();
}

/**********************************************************************
 * Dusty2D_pState::h -- Gas-phase specific enthalpy.                  *
 **********************************************************************/
inline double Dusty2D_pState::h(void) const {
  //assert(rho > ZERO);
  return g*gm1i*p/rho + HALF*v.sqr() + k;
}

/**********************************************************************
 * Dusty2D_pState::H -- Gas-phase total enthalpy.                     *
 **********************************************************************/
inline double Dusty2D_pState::H(void) const {
  return g*gm1i*p + HALF*rho*v.sqr() + dk();
}

/**********************************************************************
 * Dusty2D_pState::a -- Gas-phase sound speed.                        *
 **********************************************************************/
inline double Dusty2D_pState::a(void) const {
  //assert(rho > ZERO);
  return sqrt(g*p/rho);
}

/**********************************************************************
 * Dusty2D_pState::a2 -- Gas-phase sound speed squared.               *
 **********************************************************************/
inline double Dusty2D_pState::a2(void) const {
  //assert(rho > ZERO);
  return g*p/rho;
}

/**********************************************************************
 * Dusty2D_pState::M -- Gas-phase Mach number.                        *
 **********************************************************************/
inline double Dusty2D_pState::M(void) const {
  //assert(rho > ZERO && p > ZERO);
  return abs(v)/a();
}

/**********************************************************************
 * Dusty2D_pState::s -- Gas-phase specific entropy.                   *
 **********************************************************************/
inline double Dusty2D_pState::s(void) const {
  //assert(rho > ZERO && p > ZERO);
  return R*gm1i*log(p/pow(rho,g));
}

/**********************************************************************
 * Dusty2D_pState::dv -- Gas-phase momentum.                          *
 **********************************************************************/
inline Vector2D Dusty2D_pState::dv(void) const {
  return rho*v;
}

/**********************************************************************
 * Dusty2D_pState::dv -- Gas-phase momentum.                          *
 **********************************************************************/
inline double Dusty2D_pState::dv(const Vector2D &n) const {
  return rho*(v*n);
}

/**********************************************************************
 * Dusty2D_pState::To -- Gas-phase stagnation temperature.            *
 **********************************************************************/
inline double Dusty2D_pState::To(void) const {
  //assert(rho > ZERO && p > ZERO);
  return (p/(rho*R))*(ONE+HALF*gm1*v.sqr()/(g*p/rho));
}

/**********************************************************************
 * Dusty2D_pState::po -- Gas-phase stagnation pressure.               *
 **********************************************************************/
inline double Dusty2D_pState::po(void) const {
  //assert(rho > ZERO && p > ZERO);
  return p*pow(ONE+HALF*gm1*v.sqr()/(g*p/rho),g*gm1i);
}

/**********************************************************************
 * Dusty2D_pState::ao -- Gas-phase stagnation sound speed.            *
 **********************************************************************/
inline double Dusty2D_pState::ao(void) const {
  //assert(rho > ZERO && p > ZERO);
  return sqrt((g*p/rho)*(1.0+HALF*gm1*v.sqr()/(g*p/rho)));
}

/**********************************************************************
 * Dusty2D_pState::ho -- Gas-phase stagnation enthalpy.               *
 **********************************************************************/
inline double Dusty2D_pState::ho(void) const {
  //assert(rho > ZERO && p > ZERO);
  return (g*p/(gm1*rho) + HALF*v.sqr())*(1.0+HALF*gm1*v.sqr()/(g*p/rho));
}

/**********************************************************************
 * Dusty2D_pState::mu -- Gas-phase dynamic viscosity.                 *
 **********************************************************************/
inline double Dusty2D_pState::mu(void) const {
  return mu_gottlieb(v1,v2,v3,v4,v5,T());
}

/**********************************************************************
 * Dusty2D_pState::nu -- Gas-phase kinematic viscosity.               *
 **********************************************************************/
inline double Dusty2D_pState::nu(void) const {
  return mu()/rho;
}

/**********************************************************************
 * Dusty2D_pState::kappa -- Gas-phase thermal heat conductivity.      *
 **********************************************************************/
inline double Dusty2D_pState::kappa(void) const {
  return kappa_gottlieb(v1,v2,v3,v4,v5,T(),g,cp);
}

/**********************************************************************
 * Dusty2D_pState::Pr -- Prandtl number.                              *
 **********************************************************************/
inline double Dusty2D_pState::Pr(void) const {
  //assert(kappa() > ZERO);
  return cp*mu()/kappa();
}

/**********************************************************************
 * Dusty2D_pState::meanfreepath -- Gas-phase mean free path.          *
 **********************************************************************/
inline double Dusty2D_pState::meanfreepath(void) const {
  //assert(rho > ZERO && T() > ZERO);
  return 16.0*mu()/(5.0*rho*sqrt(2.0*PI*R*T()));
}

/**********************************************************************
 * Dusty2D_pState::tauv -- Momentum transfer relaxation time-scale.   *
 **********************************************************************/
inline double Dusty2D_pState::tauv(void) const {
  return mp/(3.0*PI*dp*mu());
}

/**********************************************************************
 * Dusty2D_pState::tauT -- Heat transfer relaxation time-scale.       *
 **********************************************************************/
inline double Dusty2D_pState::tauT(void) const {
  return mp*cp/(2.0*PI*dp*kappa());
}

/**********************************************************************
 * Dusty2D_pState::Rep -- Slip Reynolds number.                       *
 **********************************************************************/
inline double Dusty2D_pState::Rep(const Vector2D &u) const {
  return (rho*dp/mu())*sqrt((v.x-u.x)*(v.x-u.x) + (v.y-u.y)*(v.y-u.y));
}

/**********************************************************************
 * Dusty2D_pState::fRep -- Stokes drag correction.                    *
 **********************************************************************/
inline double Dusty2D_pState::fRep(const Vector2D &u) const {
  double rep = max(Rep(u),TOLER*TOLER);
  return 1.0 + d1*pow(rep,d2) + (rep/24.0)*d3/(1.0 + d4*pow(rep,d5));
  //return 1.0 + pow(Rep(u),2.0/3.0)/SIX;
}

/**********************************************************************
 * Dusty2D_pState::Kn -- Particle Knudsen number.                     *
 **********************************************************************/
inline double Dusty2D_pState::Kn(void) const {
  return meanfreepath()/dp;
}

/**********************************************************************
 * Dusty2D_pState::Cc -- Cunningham correction.                       *
 **********************************************************************/
inline double Dusty2D_pState::Cc(void) const {
  return 1.0 + Kn()*(2.492 + 0.84*exp(-1.74/Kn()));
}

/**********************************************************************
 * Dusty2D_pState::Nu -- Nusselt number.                              *
 **********************************************************************/
inline double Dusty2D_pState::Nu(const Vector2D &u) const {
  //return 2.0 + 0.6*pow(Pr(),1.0/3.0)*sqrt(Rep(u));
  return 1.0 + 0.3*pow(Pr(),1.0/3.0)*sqrt(Rep(u));
}

/**********************************************************************
 * Dusty2D_pState::chi -- Particle-phase loading factor.              *
 **********************************************************************/
inline double Dusty2D_pState::chi(void) const {
  return Wp.sigma()/rho;
}

/**********************************************************************
 * Dusty2D_pState::zetap -- Particle-phase volume fraction.           *
 **********************************************************************/
inline double Dusty2D_pState::zetap(void) const {
  return Wp.sigma()/rhop;
}

/**********************************************************************
 * Dusty2D_pState::phip -- Particle-phase mass fraction.              *
 **********************************************************************/
inline double Dusty2D_pState::phip(void) const {
  return Wp.sigma()/((1.0-zetap())*rho + Wp.sigma());;
}

/**********************************************************************
 * Dusty2D_pState::MultiVelocity_Switch                               *
 **********************************************************************/
// inline void Dusty2D_pState::MultiVelocity_Switch(void) {
//   assert(NUM_CMP_PART == PARTICLE2D_MULTI_VELOCITY_FORMULATION);
//   Particle2D_pComponents Wcomp;
//   // Consider each component of the particle family and switch
//   // compoments if necessary.
//   for (int nc = 0; nc < NUM_CMP_PART; nc++) {
//     if (Wp[nc].sigma >= NANO) {
//       if (fabs(Wp[nc].u.x) < NANO) Wp[nc].u.x = ZERO;
//       if (fabs(Wp[nc].u.y) < NANO) Wp[nc].u.y = ZERO;
//       if (Wp[nc].u.x > ZERO && Wp[nc].u.y >= ZERO) {
// 	Wcomp[0] += Particle2D_pState(Wp[nc].sigma,Wp[nc].sigma*Wp[nc].u.x,Wp[nc].sigma*Wp[nc].u.y,Wp[nc].sigma*Wp[nc].Tp);
//       } else if (Wp[nc].u.x <= ZERO && Wp[nc].u.y > ZERO) {
// 	Wcomp[1] += Particle2D_pState(Wp[nc].sigma,Wp[nc].sigma*Wp[nc].u.x,Wp[nc].sigma*Wp[nc].u.y,Wp[nc].sigma*Wp[nc].Tp);
//       } else if (Wp[nc].u.x < ZERO && Wp[nc].u.y <= ZERO) {
// 	Wcomp[2] += Particle2D_pState(Wp[nc].sigma,Wp[nc].sigma*Wp[nc].u.x,Wp[nc].sigma*Wp[nc].u.y,Wp[nc].sigma*Wp[nc].Tp);
//       } else if (Wp[nc].u.x >= ZERO && Wp[nc].u.y < ZERO) {
// 	Wcomp[3] += Particle2D_pState(Wp[nc].sigma,Wp[nc].sigma*Wp[nc].u.x,Wp[nc].sigma*Wp[nc].u.y,Wp[nc].sigma*Wp[nc].Tp);
//       }
//     }
//   }
//   // Divide each (mass) state variable by the accumulated density.
//   for (int nc = 0; nc < NUM_CMP_PART; nc++) {
//     if (Wcomp[nc].sigma < NANO) {
//       Wcomp[nc].Vacuum();
//     } else {
//       Wcomp[nc].u.x = Wcomp[nc].u.x/Wcomp[nc].sigma;
//       Wcomp[nc].u.y = Wcomp[nc].u.y/Wcomp[nc].sigma;
//       Wcomp[nc].Tp = Wcomp[nc].Tp/Wcomp[nc].sigma;
//     }
//   }
//   // Copy the reset component state.
//   Wp = Wcomp;
// }

/**********************************************************************
 * Dusty2D_pState::dk -- Gas total turbulent kinetic energy.          *
 **********************************************************************/
inline double Dusty2D_pState::dk(void) const {
  return rho*k;
}

/**********************************************************************
 * Dusty2D_pState::domega -- Gas total turbulent specific dissipation.*
 **********************************************************************/
inline double Dusty2D_pState::domega(void) const {
  return rho*omega;
}

/**********************************************************************
 * Dusty2D_pState::epsilon -- Gas specific turbulent eddy dissipation.*
 **********************************************************************/
inline double Dusty2D_pState::epsilon(void) const {
  return beta_k_o*k*omega;
}

/**********************************************************************
 * Dusty2D_pState::depsilon -- Gas total turbulent eddy dissipation.  *
 **********************************************************************/
inline double Dusty2D_pState::depsilon(void) const {
  return rho*epsilon();
}

/**********************************************************************
 * Dusty2D_pState::ell -- Gas turbulent length scale.                 *
 **********************************************************************/
inline double Dusty2D_pState::ell(void) const {
  return sqrt(k)/max(omega,NANO);
}

/**********************************************************************
 * Dusty2D_pState::Mt -- Gas turbulent Mach number.                   *
 **********************************************************************/
inline double Dusty2D_pState::Mt(void) const {
  return sqrt(TWO*k/a2());
}

/**********************************************************************
 * Dusty2D_pState::Mt -- Gas turbulent Mach number.                   *
 **********************************************************************/
inline double Dusty2D_pState::Mt2(void) const {
  return TWO*k/a2();
}

/**********************************************************************
 * Dusty2D_pState::muT -- Gas turbulent eddy dynamic viscosity.       *
 **********************************************************************/
inline double Dusty2D_pState::muT(void) const {
  return rho*nuT();
}

/**********************************************************************
 * Dusty2D_pState::nuT -- Gas turbulent eddy kinematic viscosity.     *
 **********************************************************************/
inline double Dusty2D_pState::nuT(void) const {
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) return k/max(omega,TOLER);
  return ZERO;
}

/**********************************************************************
 * Dusty2D_pState::kappaT -- Gas turbulent eddy thermal heat          *
 *                           conductivity.                            *
 **********************************************************************/
inline double Dusty2D_pState::kappaT(void) const {
  return muT()*cp/PrT;
}

/**********************************************************************
 * Dusty2D_pState::c -- Turbulence modified sound speed.              *
 **********************************************************************/
inline double Dusty2D_pState::c(void) const {
  return sqrt(c2());
}

/**********************************************************************
 * Dusty2D_pState::c2 -- Turbulence modified sound speed squared.     *
 **********************************************************************/
inline double Dusty2D_pState::c2(void) const {
  //assert(rho > ZERO);
  return a2() + (2.0/3.0)*g*k;
}

/**********************************************************************
 * Dusty2D_pState::pmodified -- Turbulence modified pressure.         *
 **********************************************************************/
inline double Dusty2D_pState::pmodified(void) const {
  //assert(rho > ZERO);
  return p + (2.0/3.0)*dk();
}

/**********************************************************************
 * Dusty2D_pState::beta_k -- k-omega auxilary relation.               *
 **********************************************************************/
inline double Dusty2D_pState::beta_k(const Dusty2D_pState &dWdx,
				     const Dusty2D_pState &dWdy) const {
  return beta_k_o*f_beta_k(dWdx,dWdy)*(ONE + xi*(Mt2()-sqr(Mto))*heaviside(Mt()-Mto));
}

/**********************************************************************
 * Dusty2D_pState::beta_omega -- k-omega auxilary relation.           *
 **********************************************************************/
inline double Dusty2D_pState::beta_omega(const Dusty2D_pState &dWdx,
					 const Dusty2D_pState &dWdy) const {
  return beta_omega_o*f_beta_omega(dWdx,dWdy) - beta_k_o*f_beta_k(dWdx,dWdy)*xi*(Mt2()-sqr(Mto))*heaviside(Mt()-Mto);
}

/**********************************************************************
 * Dusty2D_pState::f_beta_k -- k-omega auxilary relation.             *
 **********************************************************************/
inline double Dusty2D_pState::f_beta_k(const Dusty2D_pState &dWdx,
				       const Dusty2D_pState &dWdy) const {
  double chi = chi_k(dWdx,dWdy);
  if (chi <= ZERO) return ONE;
  return (ONE + 680.0*chi*chi)/(ONE + 400.0*chi*chi);
}

/**********************************************************************
 * Dusty2D_pState::f_beta_omega -- k-omega auxilary relation.         *
 **********************************************************************/
inline double Dusty2D_pState::f_beta_omega(const Dusty2D_pState &dWdx,
					   const Dusty2D_pState &dWdy) const {
  double chi = chi_omega(dWdx,dWdy);
  return (ONE + 70.0*chi)/(ONE + 80.0*chi);
}

/**********************************************************************
 * Dusty2D_pState::chi_k -- k-omega auxilary relation.                *
 **********************************************************************/
inline double Dusty2D_pState::chi_k(const Dusty2D_pState &dWdx,
				    const Dusty2D_pState &dWdy) const {
  return (dWdx.k*dWdx.omega + dWdy.k*dWdy.omega)/max(TOLER,cube(omega));
}

/**********************************************************************
 * Dusty2D_pState::chi_omega -- k-omega auxilary relation.            *
 **********************************************************************/
inline double Dusty2D_pState::chi_omega(const Dusty2D_pState &dWdx,
					const Dusty2D_pState &dWdy) const {
  return 0.25*fabs(sqr(dWdx.v.y - dWdy.v.x)*(dWdx.v.x + dWdy.v.y)/max(TOLER,cube(beta_omega_o*omega)));
}

/**********************************************************************
 * Dusty2D_pState::burningrate -- Solid propellent burning rate.      *
 **********************************************************************/
inline double Dusty2D_pState::burningrate(void) const {
  return -beta*pow(p,n);
}

/**********************************************************************
 * Dusty2D_pState -- Binary arithmetic operators.                     *
 **********************************************************************/
inline Dusty2D_pState Dusty2D_pState::operator +(const Dusty2D_pState &W) const {
  Dusty2D_pState Wtemp(*this);
  Wtemp += W;
  return Wtemp;
}

inline Dusty2D_pState Dusty2D_pState::operator -(const Dusty2D_pState &W) const {
  Dusty2D_pState Wtemp(*this);
  Wtemp -= W;
  return Wtemp;
}

// Inner product operator.
inline double Dusty2D_pState::operator *(const Dusty2D_pState &W) const {
  double sum;
  sum = rho*W.rho + v.x*W.v.x + v.y*W.v.y + p*W.p + k*W.k + omega*W.omega;
  if (NUM_VAR_PART) sum += Wp*W.Wp;
  return sum;
}

inline Dusty2D_pState Dusty2D_pState::operator *(const double &a) const {
  Dusty2D_pState Wtemp(*this);
  Wtemp.rho *= a; Wtemp.v.x *= a; Wtemp.v.y *= a; Wtemp.p *= a; Wtemp.k *= a; Wtemp.omega *= a;
  if (NUM_VAR_PART) Wtemp.Wp *= a;
  return Wtemp;
}

inline Dusty2D_pState operator *(const double &a, const Dusty2D_pState &W) {
  Dusty2D_pState Wtemp(W);
  Wtemp.rho *= a; Wtemp.v.x *= a; Wtemp.v.y *= a; Wtemp.p *= a; Wtemp.k *= a; Wtemp.omega *= a;
  if (W.NUM_VAR_PART) Wtemp.Wp *= a;
  return Wtemp;
}

inline Dusty2D_pState Dusty2D_pState::operator /(const double &a) const {
  Dusty2D_pState Wtemp(*this);
  Wtemp.rho /= a; Wtemp.v.x /= a; Wtemp.v.y /= a; Wtemp.p /= a; Wtemp.k /= a; Wtemp.omega /= a;
  if (NUM_VAR_PART) Wtemp.Wp /= a;
  return Wtemp;
}

// A useful solution state product operator.
inline Dusty2D_pState Dusty2D_pState::operator ^(const Dusty2D_pState &W) const {
  Dusty2D_pState Wtemp;
  Wtemp.rho = rho*W.rho; Wtemp.v.x = v.x*W.v.x; Wtemp.v.y = v.y*W.v.y; Wtemp.p = p*W.p;
  Wtemp.k = k*W.k; Wtemp.omega = omega*W.omega;
  if (NUM_VAR_PART) Wtemp.Wp = Wp^W.Wp;  
  return Wtemp;
}

/**********************************************************************
 * Dusty2D_pState -- Assignment operator.                             *
 **********************************************************************/
inline Dusty2D_pState& Dusty2D_pState::operator =(const Dusty2D_pState &W) {
  //if (this != &W) {
    rho = W.rho; v.x = W.v.x; v.y = W.v.y; p = W.p; k = W.k; omega = W.omega;
    if (NUM_VAR_PART) Wp = W.Wp;  
    //}
  return *this;
}

/**********************************************************************
 * Dusty2D_pState -- Unary arithmetic operators.                      *
 **********************************************************************/
//inline Dusty2D_pState operator +(const Dusty2D_pState &W) {
//return W;
//}

inline Dusty2D_pState operator -(const Dusty2D_pState &W) {
  Dusty2D_pState Wtemp;
  Wtemp.rho = -W.rho; Wtemp.v.x = -W.v.x; Wtemp.v.y = -W.v.y; Wtemp.p = -W.p; Wtemp.k = -W.k; Wtemp.omega = -W.omega;
  if (W.NUM_VAR_PART) Wtemp.Wp = -W.Wp;  
  return Wtemp;
}

/**********************************************************************
 * Dusty2D_pState -- Shortcut arithmetic operators.                   *
 **********************************************************************/
inline Dusty2D_pState& Dusty2D_pState::operator +=(const Dusty2D_pState &W) {
  rho += W.rho; v.x += W.v.x; v.y += W.v.y; p += W.p; k += W.k; omega += W.omega;
  if (NUM_VAR_PART) Wp += W.Wp;  
  return *this;
}

inline Dusty2D_pState& Dusty2D_pState::operator -=(const Dusty2D_pState &W) {
  rho -= W.rho; v.x -= W.v.x; v.y -= W.v.y; p -= W.p; k -= W.k; omega -= W.omega;
  if (NUM_VAR_PART) Wp -= W.Wp;
  return *this;
}

inline Dusty2D_pState& Dusty2D_pState::operator *=(const double &a) {
  rho *= a; v.x *= a; v.y *= a; p *= a; k *= a; omega *= a;
  if (NUM_VAR_PART) Wp *= a;
  return *this;
}

inline Dusty2D_pState& Dusty2D_pState::operator /=(const double &a) {
  rho /= a; v.x /= a; v.y /= a; p /= a; k /= a; omega /= a;
  if (NUM_VAR_PART) Wp /= a;
  return *this;
}

/**********************************************************************
 * Dusty2D_pState -- Relational operators.                            *
 **********************************************************************/
inline int operator ==(const Dusty2D_pState &W1, const Dusty2D_pState &W2) {
  int check = 1;
  for (int nv = 1; nv <= W1.NUM_VAR_DUSTY2D; nv++) if (W1[nv] != W2[nv]) check = 0;
  return check;
}

inline int operator !=(const Dusty2D_pState &W1, const Dusty2D_pState &W2) {
  int check = 0;
  for (int nv = 1; nv <= W1.NUM_VAR_DUSTY2D; nv++) if (W1[nv] == W2[nv]) check = 1;
  return check;
}

/**********************************************************************
 * Dusty2D_pState -- Input-output operators.                          *
 **********************************************************************/
inline ostream &operator << (ostream &out_file, const Dusty2D_pState &W) {
  out_file.setf(ios::scientific);
  out_file << " " << W.rho << " " << W.v.x << " " << W.v.y << " " << W.p << " " << W.k << " " << W.omega;
  out_file.unsetf(ios::scientific);
  if (W.NUM_VAR_PART) out_file << W.Wp;
  return out_file;
}

inline istream &operator >> (istream &in_file, Dusty2D_pState &W) {
  in_file.setf(ios::skipws);
  in_file >> W.rho >> W.v.x >> W.v.y >> W.p >> W.k >> W.omega;
  in_file.unsetf(ios::skipws);
  if (W.NUM_VAR_PART) in_file >> W.Wp;
  return in_file;
}

/**********************************************************************
 * Dusty2D_cState::allocate -- Allocate memory for data arrays.       *
 **********************************************************************/
inline void Dusty2D_cState::allocate(void) {
  if (NUM_VAR_PART) Up.allocate();
}

/**********************************************************************
 * Dusty2D_cState::deallocate -- Deallocate memory.                   *
 **********************************************************************/
inline void Dusty2D_cState::deallocate(void) {
  if (NUM_VAR_PART) Up.deallocate();
}

/**********************************************************************
 * Dusty2D_cState::Copy -- Copy operator.                             *
 **********************************************************************/
inline void Dusty2D_cState::Copy(const Dusty2D_cState &U) {
  rho = U.rho; dv.x = U.dv.x; dv.y = U.dv.y; E = U.E; dk = U.dk; domega = U.domega;
  if (NUM_VAR_PART) Up.Copy(U.Up);
}

/**********************************************************************
 * Dusty2D_cState::Vacuum -- Vacuum state.                            *
 **********************************************************************/
inline void Dusty2D_cState::Vacuum(void) {
  rho = ZERO; dv.x = ZERO; dv.y = ZERO; E = ZERO; dk = ZERO; domega = ZERO; tau.zero(); q.zero();
  if (NUM_VAR_PART) Up.Vacuum();
}

/**********************************************************************
 * Dusty2D_cState::Standard_Atmosphere -- Standard atmosphere state.  *
 **********************************************************************/
inline void Dusty2D_cState::Standard_Atmosphere(void) {
  rho = DENSITY_STDATM; dv.x = ZERO; dv.y = ZERO; E = PRESSURE_STDATM/(GAMMA_AIR-1.0);
  dk = ZERO; domega = ZERO; tau.zero(); q.zero();
  if (NUM_VAR_PART) Up.Vacuum();
}

/**********************************************************************
 * Dusty2D_cState::Unphysical_Properties -- Check for unphysical      *
 *                                          state properties.         *
 **********************************************************************/
inline int Dusty2D_cState::Unphysical_Properties(void) const {
  if (rho <= ZERO || E <= ZERO || e() <= ZERO) return 1;
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) if (dk < ZERO || domega < ZERO) return 1;
  if (NUM_VAR_PART) return Up.Unphysical_Properties();
  return 0;
}

/**********************************************************************
 * Dusty2D_cState::Copy_Multigrid_State_Variables --                  *
 *                           Copy variables solved by multigrid only. *
 **********************************************************************/
inline void Dusty2D_cState::Copy_Multigrid_State_Variables(const Dusty2D_cState &Ufine) {
  rho = Ufine.rho; dv.x = Ufine.dv.x; dv.y = Ufine.dv.y; E = Ufine.E;
}

/**********************************************************************
 * Dusty2D_cState::Zero_Non_Multigrid_State_Variables --              *
 *                            Zero variables not-solved by multigrid. *
 **********************************************************************/
inline void Dusty2D_cState::Zero_Non_Multigrid_State_Variables(void) {
  dk = ZERO; domega = ZERO;
}

/**********************************************************************
 * Dusty2D_cState::set_static_variables -- Set all static variables.  *
 **********************************************************************/
inline void Dusty2D_cState::set_static_variables(void) {
  // Set gas constants.
  set_gas("AIR");
  // Set particle constants.
  NUM_CMP_PART = 0;
  NUM_VAR_PART = 0;
  // Set total number of variables.
  NUM_VAR_DUSTY2D = NUM_VAR_BASE + NUM_VAR_PART;
  // Allocate data arrays.
  allocate();
  // Set indicing look up tables.
  set_indicing();
  // Set the flow type.
  flow_type = FLOWTYPE_INVISCID;
  // Set turbulence constants.
  set_turbulence(ZERO,ZERO,ZERO,ZERO,ZERO);
  // Set propellant type.
  set_propellant("AP_HTPB");
}


inline void Dusty2D_cState::set_static_variables(char *gas_type,
						 const int &FlowType,
						 const double &C_constant,
						 const double &von_karman_constant,
						 const double &yplus_sub,
						 const double &yplus_buffer,
						 const double &yplus_outer,
						 char *propellant_type) {
  // Deallocate memory arrays.
  deallocate();
  // Set gas constants.
  set_gas(gas_type);
  // Set particle constants.
  NUM_CMP_PART = 0;
  NUM_VAR_PART = 0;
  // Set total number of variables.
  NUM_VAR_DUSTY2D = NUM_VAR_BASE + NUM_VAR_PART;
  // Allocate data arrays.
  allocate();
  // Set indicing look up tables.
  set_indicing();
  // Set the flow type.
  flow_type = FlowType;
  // Set turbulence constants.
  set_turbulence(C_constant,von_karman_constant,yplus_sub,yplus_buffer,yplus_outer);
  // Set propellant type.
  set_propellant(propellant_type);
}

inline void Dusty2D_cState::set_static_variables(char *gas_type,
						 char *particle_type,
						 const int &num_cmp_part,
						 const int &draglaw,
						 const int &FlowType,
						 const double &C_constant,
						 const double &von_karman_constant,
						 const double &yplus_sub,
						 const double &yplus_buffer,
						 const double &yplus_outer,
						 char *propellant_type) {
  // Deallocate memory arrays.
  deallocate();
  // Set gas constants.
  set_gas(gas_type);
  // Set particle constants.
  NUM_CMP_PART = num_cmp_part;
  NUM_VAR_PART = NUM_CMP_PART*NUM_VAR_PARTICLE2D;
  if (NUM_VAR_PART) {
    set_particles(particle_type,drag_law);
    Up.set_particle_components(num_cmp_part);
  }
  // Set total number of variables.
  NUM_VAR_DUSTY2D = NUM_VAR_BASE + NUM_VAR_PART;
  // Allocate data arrays.
  allocate();
  // Set indicing look up tables.
  set_indicing();
  // Set the flow type.
  flow_type = FlowType;
  // Set turbulence constants.
  set_turbulence(C_constant,von_karman_constant,yplus_sub,yplus_buffer,yplus_outer);
  // Set propellant type.
  set_propellant(propellant_type);
}

/**********************************************************************
 * Dusty2D_cState::set_gas -- Set gas-phase static variables.         *
 **********************************************************************/
inline void Dusty2D_cState::set_gas(char *gas_type) {
  if (strcmp(gas_type,"AIR") == 0) {
    g  = GAMMA_AIR;
    R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    v1 = AIR_c1; v2 = AIR_c2; v3 = AIR_c3; v4 = AIR_c4; v5 = AIR_c5;
  } else if (strcmp(gas_type,"H2") == 0) {
    g  = GAMMA_H2;
    R  = R_UNIVERSAL/(MOLE_WT_H2*MILLI);
    v1 = H2_c1; v2 = H2_c2; v3 = H2_c3; v4 = H2_c4; v5 = H2_c5;
  } else if (strcmp(gas_type,"HE") == 0) {
    g  = GAMMA_HE;
    R  = R_UNIVERSAL/(MOLE_WT_HE*MILLI);
    v1 = HE_c1; v2 = HE_c2; v3 = HE_c3; v4 = HE_c4; v5 = HE_c5;
  } else if (strcmp(gas_type,"N2") == 0) {
    g  = GAMMA_N2;
    R  = R_UNIVERSAL/(MOLE_WT_N2*MILLI);
    v1 = N2_c1; v2 = N2_c2; v3 = N2_c3; v4 = N2_c4; v5 = N2_c5;
  } else if (strcmp(gas_type,"O2") == 0) {
    g  = GAMMA_O2;
    R  = R_UNIVERSAL/(MOLE_WT_O2*MILLI);
    v1 = O2_c1; v2 = O2_c2; v3 = O2_c3; v4 = O2_c4; v5 = O2_c5;
  } else if (strcmp(gas_type,"AP_HTPB") == 0) {
    g  = GAMMA_APHTPB;
    R  = R_UNIVERSAL/(MOLE_WT_APHTPB*MILLI);
    v1 = APHTPB_c1; v2 = APHTPB_c2; v3 = APHTPB_c3; v4 = APHTPB_c4; v5 = APHTPB_c5;
  } else {
    g  = GAMMA_AIR;
    R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    v1 = AIR_c1; v2 = AIR_c2; v3 = AIR_c3; v4 = AIR_c4; v5 = AIR_c5;
  }
  gm1  = g - 1.0;
  gm1i = 1.0/gm1;
  cp = g*R*gm1i;
}

/**********************************************************************
 * Dusty2D_cState::set_particles -- Set the particle constants and    *
 *                                  the drag law coefficients.        *
 **********************************************************************/
inline void Dusty2D_cState::set_particles(char *particle_type,
					  const int &draglaw) {
  // Set the particle type.
  if (strcmp(particle_type,"AP_HTPB") == 0) {
    cm   = CM_APHTPB;
    dp   = DP_APHTPB;
    rhop = RHOP_APHTPB;
    mp   = MP_APHTPB;
    qe   = ZERO;
  } else if (strcmp(particle_type,"Glass_Beads") == 0) {
    cm   = CM_GLASS_BEADS;
    dp   = DP_GLASS_BEADS;
    rhop = RHOP_GLASS_BEADS;
    mp   = MP_GLASS_BEADS;
    qe   = QE_GLASS_BEADS;
  } else if (strcmp(particle_type,"SCIEX_TINY") == 0) {
    cm   = CM_SCIEX_TINY;
    dp   = DP_SCIEX_TINY;
    rhop = RHOP_SCIEX_TINY;
    mp   = MP_SCIEX_TINY;
    qe   = QE_SCIEX_TINY;
  } else if (strcmp(particle_type,"SCIEX_SMALL") == 0) {
    cm   = CM_SCIEX_SMALL;
    dp   = DP_SCIEX_SMALL;
    rhop = RHOP_SCIEX_SMALL;
    mp   = MP_SCIEX_SMALL;
    qe   = QE_SCIEX_SMALL;
  } else if (strcmp(particle_type,"SCIEX_MEDIUM") == 0) {
    cm   = CM_SCIEX_MEDIUM;
    dp   = DP_SCIEX_MEDIUM;
    rhop = RHOP_SCIEX_MEDIUM;
    mp   = MP_SCIEX_MEDIUM;
    qe   = QE_SCIEX_MEDIUM;
  } else if (strcmp(particle_type,"SCIEX_LARGE") == 0) {
    cm   = CM_SCIEX_LARGE;
    dp   = DP_SCIEX_LARGE;
    rhop = RHOP_SCIEX_LARGE;
    mp   = MP_SCIEX_LARGE;
    qe   = QE_SCIEX_LARGE;
  } else if (strcmp(particle_type,"Heavy_AP_HTPB") == 0) {
    cm   = CM_HEAVY_APHTPB;
    dp   = DP_HEAVY_APHTPB;
    rhop = RHOP_HEAVY_APHTPB;
    mp   = MP_HEAVY_APHTPB;
    qe   = QE_HEAVY_APHTPB;
  } else {
    cm   = CM_APHTPB;
    dp   = DP_APHTPB;
    rhop = RHOP_APHTPB;
    mp   = MP_APHTPB;
    qe   = ZERO;
  }
  // Set the drag law.
  drag_law = draglaw;
  switch(drag_law) {
  case DRAG_LAW_STOKES :
    d1 = 0.0;
    d2 = 0.0;
    d3 = 0.0;
    d4 = 0.0;
    d5 = 0.0;
    break;
  case DRAG_LAW_OSEEN :
    d1 = 3.0/16.0;
    d2 = 1.0;
    d3 = 0.0;
    d4 = 0.0;
    d5 = 0.0;
    break;
  case DRAG_LAW_KLYACHKO_PUTNAM :
    d1 = 1.0/6.0;
    d2 = 2.0/3.0;
    d3 = 0.0;
    d4 = 0.0;
    d5 = 0.0;
    break;
  case DRAG_LAW_TURTON_LEVENSPIEL :
    d1 = 0.173;
    d2 = 0.657;
    d3 = 0.413;
    d4 = 16300.0;
    d5 = -1.09;
    break;
  case DRAG_LAW_HAIDER_LEVENSPIEL :
    d1 = 0.1806;
    d2 = 0.6459;
    d3 = 0.4251;
    d4 = 6880.95;
    d5 = -1.0;
    break;
  };
}

/**********************************************************************
 * Dusty2D_cState::set_turbulence -- Set the turbulence static        *
 *                                   variables.                       *
 **********************************************************************/
inline void Dusty2D_cState::set_turbulence(const double &C_constant,
					   const double &von_karman_constant,
					   const double &yplus_sub,
					   const double &yplus_buffer,
					   const double &yplus_outer) {
  // k-omega closure coefficients:
  Cmu = 0.09;
  beta_k_o = 0.09;
  beta_omega_o = 0.072;
  sigma_k = 0.50;
  sigma_omega = 0.50;
  alpha = 0.52;
  // Turbulent boundary-layer constants:
  C = C_constant;
  von_karman = von_karman_constant;
  yplus_sublayer = yplus_sub;
  yplus_buffer_layer = yplus_buffer;
  yplus_outer_layer = yplus_outer;
  //yplus_o = Iterative_Sub_To_Log_Layer_Transition_Point(C,von_karman);
  double f, df, E = exp(von_karman*C);
  // Set the initial guess.
  yplus_o = 10.0;
  // Iterate to determine the transition point.
  do {
    f = von_karman*yplus_o - log(E*yplus_o);
    df = von_karman - ONE/max(yplus_o,TOLER*TOLER);
    yplus_o -= f/df;
  } while(fabs(f) >= 0.00000001);
}

/**********************************************************************
 * Dusty2D_cState::set_propellant -- Set propellant static variables. *
 **********************************************************************/
inline void Dusty2D_cState::set_propellant(char *propellant_type) {
  if (strcmp(propellant_type,"AP_HTPB") == 0) {
    rhos   = RHOS_APHTPB;
    n      = N_APHTPB;
    beta   = BETA_APHTPB;
    Tf     = TF_APHTPB;
    Ts     = TS_APHTPB;
    alphas = ALPHAS_APHTPB;
  } else if (strcmp(propellant_type,"QUICK_AP_HTPB") == 0) {
    rhos   = RHOS_QUICK;
    n      = N_QUICK;
    beta   = BETA_QUICK;
    Tf     = TF_QUICK;
    Ts     = TS_QUICK;
    alphas = ALPHAS_QUICK;
  } else if (strcmp(propellant_type,"PLAID_AP_HTPB") == 0) {
    rhos   = RHOS_PLAID;
    n      = N_PLAID;
    beta   = BETA_PLAID;
    Tf     = TF_PLAID;
    Ts     = TS_PLAID;
    alphas = ALPHAS_PLAID;
  } else {
    rhos   = RHOS_APHTPB;
    n      = N_APHTPB;
    beta   = BETA_APHTPB;
    Tf     = TF_APHTPB;
    Ts     = TS_APHTPB;
    alphas = ALPHAS_APHTPB;
  }
}

/**********************************************************************
 * Dusty2D_cState::set_indicing -- Set indicing look up tables.       *
 **********************************************************************/
inline void Dusty2D_cState::set_indicing(void) {
  int nc, np;
  if (cindex != NULL) { delete []cindex; cindex = NULL; }
  if (pindex != NULL) { delete []pindex; pindex = NULL; }
  cindex = new int[NUM_VAR_DUSTY2D+1]; cindex[0] = 0;
  pindex = new int[NUM_VAR_DUSTY2D+1]; pindex[0] = 0;
  for (int index = 1; index <= NUM_VAR_DUSTY2D; index++) {
    if (index <= NUM_VAR_BASE) {
      // Gas-phase variables.
      cindex[index] = index;
      pindex[index] = index;
    } else if (index <= NUM_VAR_BASE+NUM_VAR_PART) {
      // Particle-phase variables.
      Up.get_indices(index-NUM_VAR_BASE-1,nc,np);
      cindex[index] = nc;
      pindex[index] = np;
    }
  }
}

/**********************************************************************
 * Dusty2D_cState::v -- Gas-phase flow velocity.                      *
 **********************************************************************/
inline Vector2D Dusty2D_cState::v(void) const {
  //assert(rho > ZERO);
  return dv/rho;
}

/**********************************************************************
 * Dusty2D_cState::v -- Gas-phase flow velocity.                      *
**********************************************************************/
inline double Dusty2D_cState::v(const Vector2D &n) const {
  //assert(rho > ZERO);
  return (dv*n)/rho;
}

/**********************************************************************
 * Dusty2D_cState::p -- Gas-phase pressure.                           *
 **********************************************************************/
inline double Dusty2D_cState::p(void) const {
  //assert(rho > ZERO);
  return gm1*(E - HALF*dv.sqr()/rho - dk);
}

/**********************************************************************
 * Dusty2D_cState::T -- Gas-phase temperature.                        *
 **********************************************************************/
inline double Dusty2D_cState::T(void) const {
  //assert(rho > ZERO);
  return p()/(rho*R);
}

/**********************************************************************
 * Dusty2D_cState::e -- Gas-phase specific internal energy.           *
 **********************************************************************/
inline double Dusty2D_cState::e(void) const {
  //assert(rho > ZERO);
  return p()/(gm1*rho);
}

/**********************************************************************
 * Dusty2D_cState::h -- Gas-phase specific enthalpy.                  *
 **********************************************************************/
inline double Dusty2D_cState::h(void) const {
  //assert(rho > ZERO);
  return g*gm1i*p()/rho + HALF*dv.sqr()/sqr(rho) + k();
}

/**********************************************************************
 * Dusty2D_cState::H -- Gas-phase total enthalpy.                     *
 **********************************************************************/
inline double Dusty2D_cState::H(void) const {
  return g*gm1i*p() + HALF*dv.sqr()/rho + dk;
}

/**********************************************************************
 * Dusty2D_cState::a -- Gas-phase sound speed.                        *
 **********************************************************************/
inline double Dusty2D_cState::a(void) const {
  //assert(rho > ZERO);
  return sqrt(g*p()/rho);
}

/**********************************************************************
 * Dusty2D_cState::a2 -- Gas-phase sound speed squared.               *
 **********************************************************************/
inline double Dusty2D_cState::a2(void) const {
  //assert(rho > ZERO);
  return g*p()/rho;
}

/**********************************************************************
 * Dusty2D_cState::M -- Gas-phase Mach number.                        *
 **********************************************************************/
inline double Dusty2D_cState::M(void) const {
  //assert(rho > ZERO);
  return abs(v())/a();
}

/**********************************************************************
 * Dusty2D_cState::s -- Gas-phase specific entropy.                   *
 **********************************************************************/
inline double Dusty2D_cState::s(void) const {
  //assert(rho > ZERO);
  return R*gm1i*log(p()/pow(rho,g));
  //return R*gm1i*log(gm1*(E - HALF*dv.sqr()/rho)/pow(rho,g));
}

/**********************************************************************
 * Dusty2D_cState::To -- Gas-phase stagnation temperature.            *
 **********************************************************************/
inline double Dusty2D_cState::To(void) const {
  //assert(rho > ZERO);
  return (gm1*(E - HALF*dv.sqr()/rho)/(rho*R))*
	 (1.0+HALF*gm1*dv.sqr()/(rho*rho*g*gm1*(E/rho - HALF*dv.sqr()/sqr(rho))));
}

/**********************************************************************
 * Dusty2D_cState::po -- Gas-phase stagnation pressure.               *
 **********************************************************************/
inline double Dusty2D_cState::po(void) const {
  return (gm1*(E - HALF*dv.sqr()/rho))*
	  pow(1.0+HALF*gm1*dv.sqr()/(rho*rho*g*gm1*(E/rho - HALF*dv.sqr()/sqr(rho))),g*gm1i);
}

/**********************************************************************
 * Dusty2D_cState::ao -- Gas-phase stagnation sound speed.            *
 **********************************************************************/
inline double Dusty2D_cState::ao(void) const {
  return sqrt((g*gm1*(E/rho - HALF*dv.sqr()/sqr(rho)))*
	      (1.0+HALF*gm1*dv.sqr()/(rho*rho*g*gm1*(E/rho - HALF*dv.sqr()/sqr(rho)))));
}

/**********************************************************************
 * Dusty2D_cState::ho -- Gas-phase stagnation enthalpy.               *
 **********************************************************************/
inline double Dusty2D_cState::ho(void) const {
  return (g*E/rho - gm1*HALF*dv.sqr()/sqr(rho))*
         (1.0+HALF*gm1*dv.sqr()/(rho*rho*g*gm1*(E/rho - HALF*dv.sqr()/sqr(rho))));
}

/**********************************************************************
 * Dusty2D_cState::mu -- Gas-phase dynamic viscosity.                 *
 **********************************************************************/
inline double Dusty2D_cState::mu(void) const {
  return mu_gottlieb(v1,v2,v3,v4,v5,T());
}

/**********************************************************************
 * Dusty2D_cState::nu -- Gas-phase kinematic viscosity.               *
 **********************************************************************/
inline double Dusty2D_cState::nu(void) const {
  return mu()/rho;
}

/**********************************************************************
 * Dusty2D_cState::kappa -- Gas-phase thermal heat conductivity.      *
 **********************************************************************/
inline double Dusty2D_cState::kappa(void) const {
  return kappa_gottlieb(v1,v2,v3,v4,v5,T(),g,cp);
}

/**********************************************************************
 * Dusty2D_cState::Pr -- Prandtl number.                              *
 **********************************************************************/
inline double Dusty2D_cState::Pr(void) const {
  return cp*mu()/kappa();
}

/**********************************************************************
 * Dusty2D_cState:meanfreepath -- Gas-phase mean free path.           *
 **********************************************************************/
inline double Dusty2D_cState::meanfreepath(void) const {
  //assert(rho > ZERO && T() > ZERO);
  return 16.0*mu()/(5.0*rho*sqrt(2.0*PI*R*T()));
}

/**********************************************************************
 * Dusty2D_cState::tauv -- Momentum transfer relaxation time-scale.   *
 **********************************************************************/
inline double Dusty2D_cState::tauv(void) const {
  return mp/(3.0*PI*dp*mu());
}

/**********************************************************************
 * Dusty2D_cState::tauT -- Heat transfer relaxation time-scale.       *
 **********************************************************************/
inline double Dusty2D_cState::tauT(void) const {
  return mp*cp/(2.0*PI*dp*kappa());
}

/**********************************************************************
 * Dusty2D_cState::Rep -- Slip Reynolds number.                       *
 **********************************************************************/
inline double Dusty2D_cState::Rep(const Vector2D &u) const {
  return (rho*dp/mu())*sqrt((v().x-u.x)*(v().x-u.x) + (v().y-u.y)*(v().y-u.y));
}

/**********************************************************************
 * Dusty2D_cState::fRep -- Stokes drag correction.                    *
 **********************************************************************/
inline double Dusty2D_cState::fRep(const Vector2D &u) const {
  double rep = Rep(u);
  return 1.0 + d1*pow(rep,d2) + (rep/24.0)*d3/(1.0 + d4*pow(rep,d5));
  //return 1.0 + pow(Rep(u),2.0/3.0)/SIX;
}

/**********************************************************************
 * Dusty2D_cState::Kn -- Particle Knudsen number.                     *
 **********************************************************************/
inline double Dusty2D_cState::Kn(void) const {
  return meanfreepath()/dp;
}

/**********************************************************************
 * Dusty2D_cState::Cc -- Cunningham correction.                       *
 **********************************************************************/
inline double Dusty2D_cState::Cc(void) const {
  return 1.0 + Kn()*(2.492 + 0.84*exp(-1.74/Kn()));
}

/**********************************************************************
 * Dusty2D_cState::Nu -- Nusselt number.                              *
 **********************************************************************/
inline double Dusty2D_cState::Nu(const Vector2D &u) const {
  //return 2.0 + 0.6*pow(Pr(),1.0/3.0)*sqrt(Rep(u));
  return 1.0 + 0.3*pow(Pr(),1.0/3.0)*sqrt(Rep(u));
}

/**********************************************************************
 * Dusty2D_cState::chi -- Particle-phase loading factor.              *
 **********************************************************************/
inline double Dusty2D_cState::chi(void) const {
  return Up.sigma()/rho;
}

/**********************************************************************
 * Dusty2D_cState::zetap -- Particle-phase volume fraction.           *
 **********************************************************************/
inline double Dusty2D_cState::zetap(void) const {
  return Up.sigma()/rhop;
}

/**********************************************************************
 * Dusty2D_cState::phip -- Particle-phase mass fraction.              *
 **********************************************************************/
inline double Dusty2D_cState::phip(void) const {
  return Up.sigma()/((1.0-zetap())*rho + Up.sigma());;
}

/**********************************************************************
 * Dusty2D_cState::dk -- Gas specific turbulent kinetic energy.       *
 **********************************************************************/
inline double Dusty2D_cState::k(void) const {
  return dk/rho;
}

/**********************************************************************
 * Dusty2D_cState::depsilon -- Gas total turbulent eddy dissipation.  *
 **********************************************************************/
inline double Dusty2D_cState::depsilon(void) const {
  return beta_k_o*dk*domega/rho;
}

/**********************************************************************
 * Dusty2D_cState::depsilon -- Gas specific turbulent eddy            *
 *                             dissipation.                           *
 **********************************************************************/
inline double Dusty2D_cState::epsilon(void) const {
  return depsilon()/rho;
}

/**********************************************************************
 * Dusty2D_cState::omega -- Gas specific turbulent dissipation rate.  *
 **********************************************************************/
inline double Dusty2D_cState::omega(void) const {
  return domega/rho;
}

/**********************************************************************
 * Dusty2D_cState::ell -- Return the turbulent length scale.          *
 **********************************************************************/
inline double Dusty2D_cState::ell(void) const {
  return sqrt(k())/max(omega(),NANO);
}

/**********************************************************************
 * Dusty2D_cState::Mt -- Return the turbulent Mach number.            *
 **********************************************************************/
inline double Dusty2D_cState::Mt(void) const {
  return sqrt(TWO*k()/a2());
}

/**********************************************************************
 * Dusty2D_cState::Mt2 -- Gas turbulent Mach number squared.          *
 **********************************************************************/
inline double Dusty2D_cState::Mt2(void) const {
  return TWO*k()/a2();
}

/**********************************************************************
 * Dusty2D_cState::muT -- Turbulent eddy dynamic viscosity.           *
 **********************************************************************/
inline double Dusty2D_cState::muT(void) const {
  return rho*nuT();
}

/**********************************************************************
 * Dusty2D_cState::nuT -- Turbulent eddy kinematic viscosity.         *
 **********************************************************************/
inline double Dusty2D_cState::nuT(void) const {
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) return k()/max(omega(),TOLER);
  return ZERO;
}

/**********************************************************************
 * Dusty2D_cState::kappaT -- Turbulent eddy thermal heat conductivity.*
 **********************************************************************/
inline double Dusty2D_cState::kappaT(void) const {
  return muT()*cp/PrT;
}

/**********************************************************************
 * Dusty2D_cState::c -- Turbulence modified sound speed.              *
 **********************************************************************/
inline double Dusty2D_cState::c(void) const {
  return sqrt(c2());
}

/**********************************************************************
 * Dusty2D_cState::c2 -- Turbulence modified sound speed squared.     *
 **********************************************************************/
inline double Dusty2D_cState::c2(void) const {
  return a2() + (2.0/3.0)*g*k();
}

/**********************************************************************
 * Dusty2D_cState::pmodified -- Turbulence modified pressure.         *
 **********************************************************************/
inline double Dusty2D_cState::pmodified(void) const {
  return p() + (2.0/3.0)*dk;
}

/**********************************************************************
 * Dusty2D_cState::beta_k -- k-omega auxilary relation.               *
 **********************************************************************/
inline double Dusty2D_cState::beta_k(const Dusty2D_pState &dWdx,
					    const Dusty2D_pState &dWdy) const {
  return beta_k_o*f_beta_k(dWdx,dWdy)*(ONE + xi*(Mt2()-sqr(Mto))*heaviside(Mt()-Mto));
}

/**********************************************************************
 * Dusty2D_cState::beta_omega -- k-omega auxilary relation.           *
 **********************************************************************/
inline double Dusty2D_cState::beta_omega(const Dusty2D_pState &dWdx,
						const Dusty2D_pState &dWdy) const {
  return beta_omega_o*f_beta_omega(dWdx,dWdy) - beta_k_o*f_beta_k(dWdx,dWdy)*xi*(Mt2()-sqr(Mto))*heaviside(Mt()-Mto);
}

/**********************************************************************
 * Dusty2D_cState::f_beta_k -- k-omega auxilary relation.             *
 **********************************************************************/
inline double Dusty2D_cState::f_beta_k(const Dusty2D_pState &dWdx,
					      const Dusty2D_pState &dWdy) const {
  double chi = chi_k(dWdx,dWdy);
  if (chi <= ZERO) return ONE;
  return (ONE + 680.0*chi*chi)/(ONE + 400.0*chi*chi);
}

/**********************************************************************
 * Dusty2D_cState::f_beta_omega -- k-omega auxilary relation.         *
 **********************************************************************/
inline double Dusty2D_cState::f_beta_omega(const Dusty2D_pState &dWdx,
						  const Dusty2D_pState &dWdy) const {
  double chi = chi_omega(dWdx,dWdy);
  return (ONE + 70.0*chi)/(ONE + 80.0*chi);
}

/**********************************************************************
 * Dusty2D_cState::chi_k -- k-omega auxilary relation.                *
 **********************************************************************/
inline double Dusty2D_cState::chi_k(const Dusty2D_pState &dWdx,
				    const Dusty2D_pState &dWdy) const {
  //return 0.0;
  return (dWdx.k*dWdx.omega + dWdy.k*dWdy.omega)/max(TOLER,cube(omega()));
}

/**********************************************************************
 * Dusty2D_cState::chi_omega -- k-omega auxilary relation.            *
 **********************************************************************/
inline double Dusty2D_cState::chi_omega(const Dusty2D_pState &dWdx,
					const Dusty2D_pState &dWdy) const {
  //return 0.0;
  return 0.25*fabs(sqr(dWdx.v.y - dWdy.v.x)*(dWdx.v.x + dWdy.v.y)/max(TOLER,cube(beta_omega_o*omega())));
}

/**********************************************************************
 * Dusty2D_cState::burningrate -- Solid propellent burning rate.      *
 **********************************************************************/
inline double Dusty2D_cState::burningrate(void) const {
  return -beta*pow(p(),n);
}

/**********************************************************************
 * Dusty2D_cState -- Binary arithmetic operators.                     *
 **********************************************************************/
inline Dusty2D_cState Dusty2D_cState::operator +(const Dusty2D_cState &U) const {
  Dusty2D_cState Utemp(*this);
  Utemp += U;
  return Utemp;
}

inline Dusty2D_cState Dusty2D_cState::operator -(const Dusty2D_cState &U) const {
  Dusty2D_cState Utemp(*this);
  Utemp -= U;
  return Utemp;
}

// Inner product operator.
inline double Dusty2D_cState::operator *(const Dusty2D_cState &U) const {
  double sum;
  sum = rho*U.rho + dv.x*U.dv.x + dv.y*U.dv.y + E*U.E + dk*U.dk + domega*U.domega;
  if (NUM_VAR_PART) sum += Up*U.Up;
  return sum;
}

inline Dusty2D_cState Dusty2D_cState::operator *(const double &a) const {
  Dusty2D_cState Utemp(*this);
  Utemp.rho *= a; Utemp.dv.x *= a; Utemp.dv.y *= a; Utemp.E *= a; Utemp.dk *= a; Utemp.domega *= a;
  if (NUM_VAR_PART) Utemp.Up *= a;
  return Utemp;
}

inline Dusty2D_cState operator *(const double &a, const Dusty2D_cState &U) {
  Dusty2D_cState Utemp(U);
  Utemp.rho *= a; Utemp.dv.x *= a; Utemp.dv.y *= a; Utemp.E *= a; Utemp.dk *= a; Utemp.domega *= a;
  if (U.NUM_VAR_PART) Utemp.Up *= a;
  return Utemp;
}

inline Dusty2D_cState Dusty2D_cState::operator /(const double &a) const {
  Dusty2D_cState Utemp(*this);
  Utemp.rho /= a; Utemp.dv.x /= a; Utemp.dv.y /= a; Utemp.E /= a; Utemp.dk /= a; Utemp.domega /= a;
  if (NUM_VAR_PART) Utemp.Up /= a;
  return Utemp;
}

// A useful solution state product operator.
inline Dusty2D_cState Dusty2D_cState::operator ^(const Dusty2D_cState &U) const {
  Dusty2D_cState Utemp;
  Utemp.rho = rho*U.rho;
  Utemp.dv.x = dv.x*U.dv.x;
  Utemp.dv.y = dv.y*U.dv.y;
  Utemp.E = E*U.E;
  Utemp.dk = dk*U.dk;
  Utemp.domega = domega*U.domega;
  if (NUM_VAR_PART) Utemp.Up = Up^U.Up;
  return Utemp;
}

/**********************************************************************
 * Dusty2D_cState -- Assignment operator.                             *
 **********************************************************************/
inline Dusty2D_cState& Dusty2D_cState::operator =(const Dusty2D_cState &U) {
  //if (this != &U) {
    rho = U.rho; dv.x = U.dv.x; dv.y = U.dv.y; E = U.E; dk = U.dk; domega = U.domega;
    if (NUM_VAR_PART) Up = U.Up;
    //}
  return *this;
}

/**********************************************************************
 * Dusty2D_cState -- Unary arithmetic operators.                      *
 **********************************************************************/
//inline Dusty2D_cState operator +(const Dusty2D_cState &U) {
//return U;
//}

inline Dusty2D_cState operator -(const Dusty2D_cState &U) {
  Dusty2D_cState Utemp;
  Utemp.rho = -U.rho; Utemp.dv.x = -U.dv.x; Utemp.dv.y = -U.dv.y; Utemp.E = -U.E;
  Utemp.dk = -U.dk; Utemp.domega = -U.domega;
  if (U.NUM_VAR_PART) Utemp.Up = -U.Up;
  return Utemp;
}

/**********************************************************************
 * Dusty2D_cState -- Shortcut arithmetic operators.                   *
 **********************************************************************/
inline Dusty2D_cState& Dusty2D_cState::operator +=(const Dusty2D_cState &U) {
  rho += U.rho; dv.x += U.dv.x; dv.y += U.dv.y; E += U.E; dk += U.dk; domega += U.domega;
  if (NUM_VAR_PART) Up += U.Up;
  return *this;
}

inline Dusty2D_cState& Dusty2D_cState::operator -=(const Dusty2D_cState &U) {
  rho -= U.rho; dv.x -= U.dv.x; dv.y -= U.dv.y; E -= U.E; dk -= U.dk; domega -= U.domega;
  if (NUM_VAR_PART) Up -= U.Up;
  return *this;
}

inline Dusty2D_cState& Dusty2D_cState::operator *=(const double &a) {
  rho *= a; dv.x *= a; dv.y *= a; E *= a; dk *= a; domega *= a;
  if (NUM_VAR_PART) Up *= a;
  return *this;
}

inline Dusty2D_cState& Dusty2D_cState::operator /=(const double &a) {
  rho /= a; dv.x /= a; dv.y /= a; E /= a; dk /= a; domega /= a;
  if (NUM_VAR_PART) Up /= a;
  return *this;
}

/**********************************************************************
 * Dusty2D_cState -- Relational operators.                            *
 **********************************************************************/
inline int operator ==(const Dusty2D_cState &U1, const Dusty2D_cState &U2) {
  int check = 1;
  for (int nv = 1; nv <= U1.NUM_VAR_DUSTY2D; nv++) if (U1[nv] != U2[nv]) check = 0;
  return check;
}

inline int operator !=(const Dusty2D_cState &U1, const Dusty2D_cState &U2) {
  int check = 0;
  for (int nv = 1; nv <= U1.NUM_VAR_DUSTY2D; nv++) if (U1[nv] == U2[nv]) check = 1;
  return check;
}

/**********************************************************************
 * Dusty2D_cState -- Input-output operators.                          *
 **********************************************************************/
inline ostream &operator << (ostream &out_file, const Dusty2D_cState &U) {
  out_file.setf(ios::scientific);
  out_file << " " << U.rho << " " << U.dv.x << " " << U.dv.y << " " << U.E << " " << U.dk << " " << U.domega;
  out_file.unsetf(ios::scientific);
  if (U.NUM_VAR_PART) out_file << U.Up;
  return out_file;
}

inline istream &operator >> (istream &in_file, Dusty2D_cState &U) {
  in_file.setf(ios::skipws);
  in_file >> U.rho >> U.dv.x >> U.dv.y >> U.E >> U.dk >> U.domega;
  in_file.unsetf(ios::skipws);
  if (U.NUM_VAR_PART) in_file >> U.Up;
  return in_file;
}

/**********************************************************************
 * Dusty2D_pState::Dusty2D_pState -- Constructor.                     *
 **********************************************************************/
inline Dusty2D_pState::Dusty2D_pState(const Dusty2D_cState &U) {
  rho = U.rho; v.x = U.v().x; v.y = U.v().y; p = U.p(); k = U.k(); omega = U.omega();
  if (NUM_VAR_PART) Wp = Particle2D_pComponents(U.Up,cm);
}

/**********************************************************************
 * Dusty2D_pState::U -- Conserved solution state.                     *
 **********************************************************************/
inline Dusty2D_cState Dusty2D_pState::U(void) const {
  return U(*this);
}

inline Dusty2D_cState Dusty2D_pState::U(const Dusty2D_pState &W) const {
  Dusty2D_cState Utemp;
  Utemp.rho = W.rho; Utemp.dv.x = W.dv().x; Utemp.dv.y = W.dv().y; Utemp.E = W.E(); Utemp.dk = W.dk(); Utemp.domega = W.domega();
  if (NUM_VAR_PART) Utemp.Up = Particle2D_cComponents(W.Wp,cm);
  return Utemp;
}

inline Dusty2D_cState U(const Dusty2D_pState &W) {
  Dusty2D_cState Utemp;
  Utemp.rho = W.rho; Utemp.dv.x = W.dv().x; Utemp.dv.y = W.dv().y; Utemp.E = W.E(); Utemp.dk = W.dk(); Utemp.domega = W.domega();
  if (Utemp.NUM_VAR_PART) Utemp.Up = Particle2D_cComponents(W.Wp,W.cm);
  return Utemp;
}

/**********************************************************************
 * Dusty2D_pState::dUdW -- Jacobian of the conserved solution         *
 *                         variables with respect to the primitive    *
 *                         solution variables.                        *
 **********************************************************************/
inline void Dusty2D_pState::dUdW(DenseMatrix &dUdW) const {
  dUdW(0,0) += ONE;
  dUdW(1,0) += v.x;
  dUdW(1,1) += rho;
  dUdW(2,0) += v.y;
  dUdW(2,2) += rho;
  dUdW(3,0) += HALF*v.sqr();
  dUdW(3,1) += rho*v.x;
  dUdW(3,2) += rho*v.y;
  dUdW(3,3) += gm1i;
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    dUdW(3,0) += k;
    dUdW(3,4) += rho;
    dUdW(4,0) += k;
    dUdW(4,4) += rho;
    dUdW(5,0) += omega;
    dUdW(5,5) += rho;
  }
  // Particle-phase contribution to the transformation Jacobian.
  if (NUM_VAR_PART) Wp.dUdW(dUdW,NUM_VAR_BASE,cm);
}

/**********************************************************************
 * Dusty2D_pState::dWdU -- Jacobian of the primitive solution         *
 *                         variables with respect to the conserved    *
 *                         solution variables.                        *
 **********************************************************************/
inline void Dusty2D_pState::dWdU(DenseMatrix &dWdU) const {
  dWdU(0,0) += ONE;
  dWdU(1,0) -= v.x/rho;
  dWdU(1,1) += ONE/rho;
  dWdU(2,0) -= v.y/rho;
  dWdU(2,2) += ONE/rho;
  dWdU(3,0) += HALF*gm1*v.sqr();
  dWdU(3,1) -= gm1*v.x;
  dWdU(3,2) -= gm1*v.y;
  dWdU(3,3) += gm1;
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    dWdU(3,4) -= gm1;
    dWdU(4,0) -= k/rho;
    dWdU(4,4) += ONE/rho;
    dWdU(5,0) -= omega/rho;
    dWdU(5,5) += ONE/rho;
  }
  // Particle-phase contribution to the transformation Jacobian.
  if (NUM_VAR_PART) Wp.dWdU(dWdU,NUM_VAR_BASE,cm);
}

/**********************************************************************
 * Dusty2D_pState::F -- Solution inviscid flux (x-direction).         *
 **********************************************************************/
inline Dusty2D_cState Dusty2D_pState::F(void) const {
  Dusty2D_cState F;
  F[1] = rho*v.x;
  F[2] = rho*sqr(v.x) + p;
  F[3] = rho*v.x*v.y;
  F[4] = v.x*H();
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    F[2] += (2.0/3.0)*dk();
    F[4] += v.x*(2.0/3.0)*dk();
    F[5] = rho*v.x*k;
    F[6] = rho*v.x*omega;
  }
  return F;
}

inline Dusty2D_cState Dusty2D_pState::F(const Vector2D &V) const {
  Dusty2D_cState F;
  F[1] = rho*(v.x - V.x);
  F[2] = rho*(v.x - V.x)*v.x + p;
  F[3] = rho*(v.x - V.x)*v.y;
  F[4] = (v.x - V.x)*E() + v.x*p;
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    F[2] += (2.0/3.0)*dk();
    F[4] += v.x*(2.0/3.0)*dk();
    F[5] = rho*(v.x - V.x)*k;
    F[6] = rho*(v.x - V.x)*omega;
  }
  return F;
}

/**********************************************************************
 * Dusty2D_pState::dFdU -- Jacobian of the inviscid solution flux     *
 *                         with respect to the conserved solution     *
 *                         variables.                                 *
 **********************************************************************/
inline void Dusty2D_pState::dFdU(DenseMatrix &dFdU) const {
  dFdU(0,1) += 1.0;
  dFdU(1,0) += HALF*gm1*(sqr(v.x) + sqr(v.y)) - sqr(v.x);
  dFdU(1,1) -= v.x*(g-3.0);
  dFdU(1,2) -= v.y*gm1;
  dFdU(1,3) += gm1;
  dFdU(2,0) -= v.x*v.y;
  dFdU(2,1) += v.y;
  dFdU(2,2) += v.x;
  dFdU(3,0) -= (h() - HALF*gm1*(sqr(v.x)+sqr(v.y)))*v.x + (2.0/3.0)*rho*v.x*k;
  dFdU(3,1) +=  h() - gm1*sqr(v.x) + (2.0/3.0)*k;
  dFdU(3,2) -= v.x*v.y*gm1;
  dFdU(3,3) += v.x*g;
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    dFdU(1,4) += 5.0/3.0 - g;
    dFdU(3,4) -= (3.0*g - 5.0)*v.x/3.0;
    dFdU(4,0) -= v.x*k;
    dFdU(4,1) += k;
    dFdU(4,4) += v.x;
    dFdU(5,0) -= v.x*omega;
    dFdU(5,1) += omega;
    dFdU(5,5) += v.x;
  }
  // Particle-phase contribution to the Jacobian of the inviscid solution flux.
  if (NUM_VAR_PART) Wp.dFdU(dFdU,NUM_VAR_BASE,cm);
}

/**********************************************************************
 * Dusty2D_pState::Gx, Gy -- Solution viscous fluxes.                 *
 **********************************************************************/
inline Dusty2D_cState Dusty2D_pState::Gx(const Dusty2D_pState &dWdx) const {
  Dusty2D_cState G;
  G[1] = ZERO;
  G[2] = tau.xx; 
  G[3] = tau.xy;
  G[4] = - q.x + v.x*tau.xx + v.y*tau.xy;
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    G[4] += (mu()+sigma_k*muT())*dWdx.k;
    G[5] = (mu()+sigma_k*muT())*dWdx.k;
    G[6] = (mu()+sigma_omega*muT())*dWdx.omega;
  }
  return G;
}

inline Dusty2D_cState Dusty2D_pState::Gy(const Dusty2D_pState &dWdy) const {
  Dusty2D_cState G;
  G[1] = ZERO;
  G[2] = tau.xy; 
  G[3] = tau.yy;
  G[4] = - q.y + v.x*tau.xy + v.y*tau.yy;
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    G[4] += (mu()+sigma_k*muT())*dWdy.k;
    G[5] = (mu()+sigma_k*muT())*dWdy.k;
    G[6] = (mu()+sigma_omega*muT())*dWdy.omega;
  }
  return G;
}

/**********************************************************************
 * Dusty2D_pState::ComputeViscousTerms -- Compute viscous stress      *
 *                                        tensor and heat flux vector.*
 **********************************************************************/
inline void Dusty2D_pState::ComputeViscousTerms(const Dusty2D_pState &dWdx,
						const Dusty2D_pState &dWdy,
						const Vector2D &X,
						const int &Axisymmetric,
						const int &adiabatic_flag) {
  double div, radius, mumu, kap;
  mumu = mu() + muT(); kap = kappa() + kappaT();
  if (Axisymmetric) radius = max(X.y,TOLER);
  // Divergence of the velocity field.
  div = dWdx.v.x + dWdy.v.y;
  if (Axisymmetric) div += v.y/radius;
  // Stress tensor.
  tau.xx = 2.0*mumu*(dWdx.v.x - div/3.0);
  tau.xy = mumu*(dWdy.v.x + dWdx.v.y);
  tau.yy = 2.0*mumu*(dWdy.v.y - div/3.0);
  if (Axisymmetric) tau.zz = 2.0*mumu*(v.y/radius - div/3.0);
  else tau.zz = ZERO;
  // Heat flux components.
  //if (adiabatic_flag) {
  //q = Vector2D_ZERO;
  //} else {
  q.x = -kap*(dWdx.p - (p/rho)*dWdx.rho)/(rho*R);
  q.y = -kap*(dWdy.p - (p/rho)*dWdy.rho)/(rho*R);
  //}
}

/**********************************************************************
 * Dusty2D_pState::lambda_x -- Eigenvalue(s) (x-direction).           *
 **********************************************************************/
inline Dusty2D_pState Dusty2D_pState::lambda_x(void) const {
  return Dusty2D_pState(v.x-c(),v.x,v.x,v.x+c(),v.x,v.x);
}

inline Dusty2D_pState Dusty2D_pState::lambda_x(const Vector2D &V) const {
  return Dusty2D_pState(v.x-V.x-c(),v.x-V.x,v.x-V.x,v.x-V.x+c(),v.x-V.x,v.x-V.x);
}

/**********************************************************************
 * Dusty2D_pState::rp_x -- Primitive right eigenvector (x-direction). *
 **********************************************************************/
inline Dusty2D_pState Dusty2D_pState::rp_x(int index) const {
  //assert(index >= 1 && index <= NUM_VAR_DUSTY2D);
  switch(index) {
  case 1 :
    return Dusty2D_pState(ONE,-c()/rho,ZERO,c2()-(2.0/3.0)*k,ZERO,ZERO);
    break;
  case 2 :
    return Dusty2D_pState(ONE,ZERO,ZERO,-(2.0/3.0)*k,ZERO,ZERO);
    break;
  case 3 :
    return Dusty2D_pState(ZERO,ZERO,ONE,ZERO,ZERO,ZERO);
    break;
  case 4 :
    return Dusty2D_pState(ONE,c()/rho,ZERO,c2()-(2.0/3.0)*k,ZERO,ZERO);
    break;
  case 5 :
    return Dusty2D_pState(ZERO,ZERO,ZERO,-(2.0/3.0)*rho,ONE,ZERO);
    break;
  case 6 :
    return Dusty2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE);
    break;
  };
  return Dusty2D_pState(ONE,-c()/rho,ZERO,c2()-(2.0/3.0)*k,ZERO,ZERO);
}

/**********************************************************************
 * Dusty2D_pState::rc_x -- Conserved right eigenvector (x-direction). *
 **********************************************************************/
inline Dusty2D_cState Dusty2D_pState::rc_x(int index) const {
  //assert(index >= 1 && index <= NUM_VAR_DUSTY2D);
  switch(index) {
  case 1 :
    return Dusty2D_cState(ONE,v.x-c(),v.y,h()-c()*v.x+(2.0/3.0)*k,k,omega);
    break;
  case 2 :
    return Dusty2D_cState(ONE,v.x,v.y,HALF*v.sqr()+gm1i*(g-5.0/3.0)*k,k,omega);
    break;
  case 3 :
    return Dusty2D_cState(ZERO,ZERO,rho,rho*v.y,ZERO,ZERO);
    break;
  case 4 :
    return Dusty2D_cState(ONE,v.x+c(),v.y,h()+c()*v.x+(2.0/3.0)*k,k,omega);
    break;
  case 5 :
    return Dusty2D_cState(ZERO,ZERO,ZERO,rho*gm1i*(g-5.0/3.0),rho,ZERO);
    break;
  case 6 :
    return Dusty2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,rho);
    break;
  };
  return Dusty2D_cState(ONE,v.x-c(),v.y,h()-c()*v.x+(2.0/3.0)*k,k,omega);
}

/**********************************************************************
 * Dusty2D_pState::lp_x -- Primitive left eigenvector (x-direction).  *
 **********************************************************************/
inline Dusty2D_pState Dusty2D_pState::lp_x(int index) const {
  //assert(index >= 1 && index <= NUM_VAR_DUSTY2D);
  switch(index) {
  case 1 :
    return Dusty2D_pState(k/(3.0*c2()),-HALF*rho/c(),ZERO,HALF/c2(),rho/(3.0*c2()),ZERO);
    break;
  case 2 :
    return Dusty2D_pState(ONE-(2.0/3.0)*k/c2(),ZERO,ZERO,-ONE/c2(),-(2.0/3.0)*rho/c2(),ZERO);
    break;
  case 3 :
    return Dusty2D_pState(ZERO,ZERO,ONE,ZERO,ZERO,ZERO);
    break;
  case 4 :
    return Dusty2D_pState(k/(3.0*c2()),HALF*rho/c(),ZERO,HALF/c2(),rho/(3.0*c2()),ZERO);
    break;
  case 5 :
    return Dusty2D_pState(ZERO,ZERO,ZERO,ZERO,ONE,ZERO);
    break;
  case 6 :
    return Dusty2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE);
    break;
  };
  return Dusty2D_pState(k/(3.0*c2()),-HALF*rho/c(),ZERO,HALF/c2(),rho/(3.0*c2()),ZERO);
}

/**********************************************************************
 * Dusty2D_pState::S -- Include all source term vectors and Jacobians.*
 **********************************************************************/
inline Dusty2D_cState Dusty2D_pState::S(const Vector2D &X,
					const Dusty2D_pState &dWdx,
					const Dusty2D_pState &dWdy,
					const int &Axisymmetric) const {
  Dusty2D_cState Sall; Sall.Vacuum();
  // Include the axisymmetric source terms if required.
  if (Axisymmetric) {
    Sall = Si(X);
    if (flow_type) Sall += Sv(X,dWdy);
  }
  // Include the phase-interaction source term if required.
  if (NUM_VAR_PART) Sall += Sp(ON);
  // Include the turbulence model source term if required.
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) Sall += St(X,dWdx,dWdy,Axisymmetric);
  // Return the total source term vector.
  return Sall; 
}

inline void Dusty2D_pState::dSdU(DenseMatrix &dSdU,
				 const Vector2D &X,
				 const Dusty2D_pState &dWdx,
				 const Dusty2D_pState &dWdy,
				 const int &Axisymmetric) const {
  // Include the axisymmetric source Jacobians.
  if (Axisymmetric) {
    dSidU(dSdU,X);
    if (flow_type) dSvdU(dSdU,dWdy,X);
  }
  // Include the phase-interaction source Jacobians.
  dSpdU(dSdU,ON);
  // Include the turbulence model source term if required.
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) dStdU(dSdU,X,dWdx,dWdy,Axisymmetric);
}

/**********************************************************************
 * Dusty2D_pState::Si -- Inviscid axisymmetric source terms and       *
 *                       Jacobian.                                    *
 **********************************************************************/
inline Dusty2D_cState Dusty2D_pState::Si(const Vector2D &X) const {
  Dusty2D_cState S;
  S[1] = -rho*v.y/X.y;
  S[2] = -rho*v.x*v.y/X.y;
  S[3] = -rho*sqr(v.y)/X.y;
  S[4] = -v.y*(H()+(2.0/3.0)*dk())/X.y;
  S[5] = -v.y*dk()/X.y;
  S[6] = -v.y*domega()/X.y;
  if (NUM_VAR_PART) S.Up = Wp.Sa(X,cm);
  return S;
}

inline void Dusty2D_pState::dSidU(DenseMatrix &dSidU, const Vector2D &X) const {
  dSidU(0,2) -= 1.0/X.y;
  dSidU(1,0) += v.x*v.y/X.y;
  dSidU(1,1) -= v.y/X.y;
  dSidU(1,2) -= v.x/X.y;
  dSidU(2,0) += v.y*v.y/X.y;
  dSidU(2,2) -= 2.0*v.y/X.y;
  dSidU(3,0) += v.y*(h()+(2.0/3.0)*k)/X.y;
  dSidU(3,1) += gm1*v.x*v.y/X.y;
  dSidU(3,2) -= (h() - gm1*v.y*v.y + (2.0/3.0)*k)/X.y;
  dSidU(3,3) -= g*v.y/X.y;
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    dSidU(3,4) += (3.0*g-5.0)*v.y/(3.0*X.y);
    dSidU(4,0) += k*v.y/X.y;
    dSidU(4,2) -= k/X.y;
    dSidU(4,4) -= v.y/X.y;
    dSidU(5,0) += omega*v.y/X.y;
    dSidU(5,2) -= omega/X.y;
    dSidU(5,5) -= v.y/X.y;
  }
  // Particle-phase contribution to the inviscid axisymmetric source term.
  if (NUM_VAR_PART) Wp.dSadU(dSidU,X,NUM_VAR_BASE,cm);
}

/**********************************************************************
 * Dusty2D_pState::Sv -- Viscous axisymmetric flow source term vector *
 *                       and Jacobian.                                *
 **********************************************************************/
inline Dusty2D_cState Dusty2D_pState::Sv(const Vector2D &X,
					 const Dusty2D_pState &dWdy) const {
  Dusty2D_cState S;
  S[1] = ZERO;
  S[2] = tau.xy/X.y;
  S[3] = (tau.yy - tau.zz)/X.y;
  S[4] = (-q.y + v.x*tau.xy + v.y*tau.yy+(mu()+sigma_k*muT())*dWdy.k)/X.y;
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    S[5] = (mu()+sigma_k*muT())*dWdy.k/X.y;
    S[6] = (mu()+sigma_omega*muT())*dWdy.omega/X.y;
  }
  return S;
}

inline void Dusty2D_pState::dSvdU(DenseMatrix &dSvdU, const Dusty2D_pState &dWdy, const Vector2D &X) const {
  dSvdU(0,0) += ZERO;
  dSvdU(1,1) -= ZERO;
}

/**********************************************************************
 * Dusty2D_pState::Sp -- Phase interaction source term and Jacobian.  *
 **********************************************************************/
inline Dusty2D_cState Dusty2D_pState::Sp(const int &interaction_flag) const {
  Dusty2D_cState S; S[1] = ZERO; S[4] = ZERO;
  double tau_v, tau_T;
  if (NUM_VAR_PART) {
    for (int nc = 0; nc < NUM_CMP_PART; nc++) {
      if (Wp[nc].sigma > ZERO) {
	tau_v = tauv()*Cc()/fRep(Wp[nc].u);
	tau_T = tauT()/Nu(Wp[nc].u);
	S.Up[nc][2] = (Wp[nc].sigma/tau_v)*(v.x - Wp[nc].u.x);
	S.Up[nc][3] = (Wp[nc].sigma/tau_v)*(v.y - Wp[nc].u.y);
 	S.Up[nc][4] = (Wp[nc].sigma*cp/tau_T)*(T() - Wp[nc].Tp);
// 	if (interaction_flag == 1) {
 	  S[2] -= S.Up[nc][2];
 	  S[3] -= S.Up[nc][3];
 	  S[4] -= S.Up[nc][4];
// 	}
      }
    }
  }
  return S;
}

inline void Dusty2D_pState::dSpdU(DenseMatrix &dSpdU, const int &interaction_flag) const {
  double tau_v, tau_T;
  int ii;
  if (NUM_VAR_PART) {
    for (int nc = 0; nc < NUM_CMP_PART; nc++) {
      ii = NUM_VAR_BASE + nc*NUM_VAR_PARTICLE2D;
      if (Wp[nc].sigma > ZERO) {
	tau_v = tauv()*Cc()/fRep(Wp[nc].u);
	tau_T = tauT()/Nu(Wp[nc].u);
	if (interaction_flag == 1) {
	  dSpdU(1   ,0   ) += Wp[nc].sigma*v.x/(tau_v*rho);
	  dSpdU(1   ,1   ) -= Wp[nc].sigma/(tau_v*rho);
	  dSpdU(1   ,ii  ) -= v.x/tau_v;
	  dSpdU(1   ,ii+1) += 1.0/tau_v;
	  dSpdU(2   ,0   ) += Wp[nc].sigma*v.y/(tau_v*rho);
	  dSpdU(2   ,2   ) -= Wp[nc].sigma/(tau_v*rho);
	  dSpdU(2   ,ii  ) -= v.y/tau_v;
	  dSpdU(2   ,ii+2) += 1.0/tau_v;
	  dSpdU(3   ,0   ) += Wp[nc].sigma*(Wp[nc].u.x*v.x + Wp[nc].u.y*v.y)/(tau_v*rho);
	  dSpdU(3   ,1   ) -= Wp[nc].sigma*Wp[nc].u.x/(tau_v*rho);
	  dSpdU(3   ,2   ) -= Wp[nc].sigma*Wp[nc].u.y/(tau_v*rho);
	  dSpdU(3   ,ii  ) -= (1.0/tau_v - HALF*(cp/cm)/tau_T)*(sqr(Wp[nc].u.x)+sqr(Wp[nc].u.y)) + cp*T()/tau_T;
	  dSpdU(3   ,ii+1) += (1.0/tau_v - (cp/cm)/tau_T)*Wp[nc].u.x - (v.x-Wp[nc].u.x)/tau_v;
	  dSpdU(3   ,ii+2) += (1.0/tau_v - (cp/cm)/tau_T)*Wp[nc].u.y - (v.y-Wp[nc].u.y)/tau_v;
	  dSpdU(3   ,ii+3) += (cp/cm)/tau_T;
	}
	dSpdU(ii+1,0   ) -= Wp[nc].sigma*v.x/(tau_v*rho);
	dSpdU(ii+1,1   ) += Wp[nc].sigma/(tau_v*rho);
	dSpdU(ii+1,ii  ) += v.x/tau_v;
	dSpdU(ii+1,ii+1) -= 1.0/tau_v;
	dSpdU(ii+2,0   ) -= Wp[nc].sigma*v.y/(tau_v*rho);
	dSpdU(ii+2,2   ) += Wp[nc].sigma/(tau_v*rho);
	dSpdU(ii+2,ii  ) += v.y/tau_v;
	dSpdU(ii+2,ii+2) -= 1.0/tau_v;
	dSpdU(ii+3,0   ) -= Wp[nc].sigma*(Wp[nc].u.x*v.x + Wp[nc].u.y*v.y)/(tau_v*rho);
	dSpdU(ii+3,1   ) += Wp[nc].sigma*Wp[nc].u.x/(tau_v*rho);
	dSpdU(ii+3,2   ) += Wp[nc].sigma*Wp[nc].u.y/(tau_v*rho);
	dSpdU(ii+3,ii  ) += (1.0/tau_v - HALF*(cp/cm)/tau_T)*(sqr(Wp[nc].u.x)+sqr(Wp[nc].u.y)) + cp*T()/tau_T;
	dSpdU(ii+3,ii+1) -= (1.0/tau_v - (cp/cm)/tau_T)*Wp[nc].u.x - (v.x-Wp[nc].u.x)/tau_v;
	dSpdU(ii+3,ii+2) -= (1.0/tau_v - (cp/cm)/tau_T)*Wp[nc].u.y - (v.y-Wp[nc].u.y)/tau_v;
	dSpdU(ii+3,ii+3) -= (cp/cm)/tau_T;
      }
    }
  }
}

/**********************************************************************
 * Dusty2D_pState::St -- Turbulent source term vector and Jacobian.   *
 **********************************************************************/
inline Dusty2D_cState Dusty2D_pState::St(const Vector2D &X,
					 const Dusty2D_pState &dWdx,
					 const Dusty2D_pState &dWdy,
					 const int &Axisymmetric) const {
  Dusty2D_cState S; S.rho = ZERO; S.E = ZERO; //S.Vacuum();
  double production, mut;
  Tensor2D lambda;
  mut = muT()/(mu() + muT());
  lambda.xx = mut*tau.xx - (2.0/3.0)*dk();
  lambda.xy = mut*tau.xy;
  lambda.yy = mut*tau.yy - (2.0/3.0)*dk();
  production = lambda.xx*dWdx.v.x + lambda.xy*(dWdy.v.x + dWdx.v.y) + lambda.yy*dWdy.v.y;
  if (Axisymmetric) {
    lambda.zz = mut*tau.zz - (2.0/3.0)*dk();
    production += lambda.zz*v.y/max(X.y,TOLER);
  }
  S.dk = production-beta_k(dWdx,dWdy)*dk()*omega;
  S.domega = alpha*(omega/max(k,TOLER))*production-beta_omega(dWdx,dWdy)*rho*sqr(omega);
  return S;
}

inline void Dusty2D_pState::dStdU(DenseMatrix &dStdU,
				  const Vector2D &X,
				  const Dusty2D_pState &dWdx,
				  const Dusty2D_pState &dWdy,
				  const int &Axisymmetric) const {
  dStdU(0,0) += ZERO;
  dStdU(1,1) -= ZERO;
}

/**********************************************************************
 * Dusty2D_pState::Se -- Electric-field source term and Jacobian.     *
 **********************************************************************/
inline Dusty2D_cState Dusty2D_pState::Se(const Electrostatic2DState &We) const {
  Dusty2D_cState S; S.rho = ZERO; S.E = ZERO; //S.Vacuum();
  S.Up = Wp.Se(We.E,cm);
  return S;
}

inline void Dusty2D_pState::dSedU(DenseMatrix &dSedU, const Electrostatic2DState &We) const {
  Wp.dSedU(dSedU,NUM_VAR_BASE,We.E,cm);
}

/**********************************************************************
 * Dusty2D_cState::Dusty2D_cState -- Constructor.                     *
 **********************************************************************/
inline Dusty2D_cState::Dusty2D_cState(const Dusty2D_pState &W) {
  rho = W.rho; dv = W.dv(); E = W.E(); dk = W.dk(); domega = W.domega();
  if (NUM_VAR_PART) Up = Particle2D_cComponents(W.Wp,W.cm);
}

/**********************************************************************
 * Dusty2D_cState::W -- Primitive solution state.                     *
 **********************************************************************/
inline Dusty2D_pState Dusty2D_cState::W(void) const {
  return W(*this);
}

inline Dusty2D_pState Dusty2D_cState::W(const Dusty2D_cState &U) const {
  Dusty2D_pState Wtemp;
  Wtemp.rho = U.rho; Wtemp.v.x = U.v().x; Wtemp.v.y = U.v().y; Wtemp.p = U.p(); Wtemp.k = U.k(); Wtemp.omega = U.omega();
  if (NUM_VAR_PART) Wtemp.Wp = Particle2D_pComponents(U.Up,U.cm);
  return Wtemp;
}

inline Dusty2D_pState W(const Dusty2D_cState &U) {
  Dusty2D_pState Wtemp;
  Wtemp.rho = U.rho; Wtemp.v.x = U.v().x; Wtemp.v.y = U.v().y; Wtemp.p = U.p(); Wtemp.k = U.k(); Wtemp.omega = U.omega();
  if (Wtemp.NUM_VAR_PART) Wtemp.Wp = Particle2D_pComponents(U.Up,U.cm);
  return Wtemp;
}

/**********************************************************************
 * Dusty2D_cState::dUdW -- Jacobian of the conserved solution         *
 *                         variables with respect to the primitive    *
 *                         solution variables.                        *
 **********************************************************************/
inline void Dusty2D_cState::dUdW(DenseMatrix &dUdW) const {
  dUdW(0,0) += 1.0;
  dUdW(1,0) += dv.x/rho;
  dUdW(1,1) += rho;
  dUdW(2,0) += dv.y/rho;
  dUdW(2,2) += rho;
  dUdW(3,0) += HALF*dv.sqr()/(rho*rho);
  dUdW(3,1) += dv.x;
  dUdW(3,2) += dv.y;
  dUdW(3,3) += gm1i;
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    dUdW(3,0) += k();
    dUdW(3,4) += rho;
    dUdW(4,0) += k();
    dUdW(4,4) += rho;
    dUdW(5,0) += omega();
    dUdW(5,5) += rho;
  }
  // Particle-phase contribution to the transformation Jacobian.
  if (NUM_VAR_PART) Up.dUdW(dUdW,NUM_VAR_BASE,cm);
}

/**********************************************************************
 * Dusty2D_cState::dWdU -- Jacobian of the primitive solution         *
 *                         variables with respect to the conserved    *
 *                         solution variables.                        *
 **********************************************************************/
inline void Dusty2D_cState::dWdU(DenseMatrix &dWdU) const {
  dWdU(0,0) += 1.0;
  dWdU(1,0) -= dv.x/(rho*rho);
  dWdU(1,1) += 1.0/rho;
  dWdU(2,0) -= dv.y/(rho*rho);
  dWdU(2,2) += 1.0/rho;
  dWdU(3,0) += HALF*gm1*dv.sqr()/(rho*rho);
  dWdU(3,1) -= gm1*dv.x/rho;
  dWdU(3,2) -= gm1*dv.y/rho;
  dWdU(3,3) += gm1;
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    dWdU(3,4) -= gm1;
    dWdU(4,0) -= k()/rho;
    dWdU(4,4) += ONE/rho;
    dWdU(5,0) -= omega()/rho;
    dWdU(5,5) += ONE/rho;
  }
  // Particle-phase contribution to the transformation Jacobian.
  if (NUM_VAR_PART) Up.dWdU(dWdU,NUM_VAR_BASE,cm);
}

/**********************************************************************
 * Dusty2D_cState::F -- Solution inviscid flux (x-direction).         *
 **********************************************************************/
inline Dusty2D_cState Dusty2D_cState::F(void) const {
  Dusty2D_cState F;
  F[1] = dv.x;
  F[2] = sqr(dv.x)/rho + p();
  F[3] = dv.x*dv.y/rho;
  F[4] = dv.x*H()/rho;
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    F[2] += (2.0/3.0)*dk;
    F[4] += v().x*(2.0/3.0)*dk;
    F[5] = v().x*dk;
    F[6] = v().x*domega;
  }
  return F;
}

inline Dusty2D_cState Dusty2D_cState::F(const Vector2D &V) const {
  Dusty2D_cState F;
  double vx = v().x;
  F.rho  = rho*(vx - V.x);
  F.dv.x = (vx - V.x)*dv.x + p();
  F.dv.y = (vx - V.x)*dv.y;
  F.E    = (vx - V.x)*E + vx*p();
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    F[2] += (2.0/3.0)*dk;
    F[4] += vx*(2.0/3.0)*dk;
    F[5] = (vx - V.x)*dk;
    F[6] = (vx - V.x)*domega;
  }
  return F;
}

/**********************************************************************
 * Dusty2D_cState::dFdU -- Jacobian of the inviscid solution flux     *
 *                         with respect to the conserved solution     *
 *                         variables.                                 *
 **********************************************************************/
inline void Dusty2D_cState::dFdU(DenseMatrix &dFdU) const {
  dFdU(0,1) += 1.0;
  dFdU(1,0) += HALF*gm1*(sqr(v().x) + sqr(v().y)) - sqr(v().x);
  dFdU(1,1) -= v().x*(g-3.0);
  dFdU(1,2) -= v().y*gm1;
  dFdU(1,3) += gm1;
  dFdU(2,0) -= v().x*v().y;
  dFdU(2,1) += v().y;
  dFdU(2,2) += v().x;
  dFdU(3,0) -= (h() - HALF*gm1*(sqr(v().x)+sqr(v().y)))*v().x + (2.0/3.0)*dv.x*k();
  dFdU(3,1) +=  h() - gm1*sqr(v().x) + (2.0/3.0)*k();
  dFdU(3,2) -= v().x*v().y*gm1;
  dFdU(3,3) += v().x*g;
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    dFdU(1,4) += 5.0/3.0 - g;
    dFdU(3,4) -= (g - 5.0/3.0)*v().x;
    dFdU(4,0) -= v().x*k();
    dFdU(4,1) += k();
    dFdU(4,4) += v().x;
    dFdU(5,0) -= v().x*omega();
    dFdU(5,1) += omega();
    dFdU(5,5) += v().x;
  }
  // Particle-phase contribution to the inviscid solution flux Jacobian.
  if (NUM_VAR_PART) Up.dFdU(dFdU,NUM_VAR_BASE);
}

/**********************************************************************
 * Dusty2D_cState::Gx, Gy -- Solution viscous fluxes.                 *
 **********************************************************************/
inline Dusty2D_cState Dusty2D_cState::Gx(const Dusty2D_pState &dWdx) const {
  Dusty2D_cState G;
  G[1] = ZERO;
  G[2] = tau.xx; 
  G[3] = tau.xy;
  G[4] = - q.x + v().x*tau.xx + v().y*tau.xy;
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    G[4] += (mu()+sigma_k*muT())*dWdx.k;
    G[5] = (mu()+sigma_k*muT())*dWdx.k;
    G[6] = (mu()+sigma_omega*muT())*dWdx.omega;
  }
  return G;
}

inline Dusty2D_cState Dusty2D_cState::Gy(const Dusty2D_pState &dWdy) const {
  Dusty2D_cState G;
  G[1] = ZERO;
  G[2] = tau.xy; 
  G[3] = tau.yy;
  G[4] = - q.y + v().x*tau.xy + v().y*tau.yy;
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    G[4] += (mu()+sigma_k*muT())*dWdy.k;
    G[5] = (mu()+sigma_k*muT())*dWdy.k;
    G[6] = (mu()+sigma_omega*muT())*dWdy.omega;
  }
  return G;
}

/**********************************************************************
 * Dusty2D_cState::ComputeViscousTerms -- Compute viscous stress      *
 *                                        tensor and heat flux vector.*
 **********************************************************************/
inline void Dusty2D_cState::ComputeViscousTerms(const Dusty2D_pState &dWdx,
						const Dusty2D_pState &dWdy,
						const Vector2D &X,
						const int &Axisymmetric,
						const int &adiabatic_flag) {
  double div, radius, mumu, kap;
  mumu = mu() + muT(); kap = kappa() + kappaT();
  if (Axisymmetric) radius = max(X.y,TOLER);
  // Divergence of the velocity field.
  div = dWdx.v.x + dWdy.v.y;
  if (Axisymmetric) div += v().y/radius;
  // Stress tensor.
  tau.xx = 2.0*mumu*(dWdx.v.x - div/3.0);
  tau.xy = mumu*(dWdy.v.x + dWdx.v.y);
  tau.yy = 2.0*mumu*(dWdy.v.y - div/3.0);
  if (Axisymmetric) tau.zz = 2.0*mumu*(v().y/radius - div/3.0);
  else tau.zz = ZERO;
  // Heat flux components.
  //if (adiabatic_flag) {
  //q = Vector2D_ZERO;
  //} else {
  q.x = -kap*(dWdx.p - (p()/rho)*dWdx.rho)/(rho*R);
  q.y = -kap*(dWdy.p - (p()/rho)*dWdy.rho)/(rho*R);
  //}
}

/**********************************************************************
 * Dusty2D_cState::lambda_x -- Eigenvalue(s) (x-direction).           *
 **********************************************************************/
inline Dusty2D_pState Dusty2D_cState::lambda_x(void) const {
  double vx = v().x, cc = c();
  return Dusty2D_pState(vx-cc,vx,vx,vx-cc,vx,vx);
}

inline Dusty2D_pState Dusty2D_cState::lambda_x(const Vector2D &V) const {
  double vx = v().x, cc = c();
  return Dusty2D_pState(vx-V.x-cc,vx-V.x,vx-V.x,vx-V.x+cc,vx-V.x,vx-V.x);
}

/**********************************************************************
 * Dusty2D_pState::rp_x -- Primitive right eigenvector (x-direction). *
 **********************************************************************/
inline Dusty2D_pState Dusty2D_cState::rp_x(int index) const {
  //assert(index >= 1 && index <= NUM_VAR_DUSTY2D);
  switch(index) {
  case 1 :
    return Dusty2D_pState(ONE,-c()/rho,ZERO,c2()-(2.0/3.0)*k(),ZERO,ZERO);
    break;
  case 2 :
    return Dusty2D_pState(ONE,ZERO,ZERO,-(2.0/3.0)*k(),ZERO,ZERO);
    break;
  case 3 :
    return Dusty2D_pState(ZERO,ZERO,ONE,ZERO,ZERO,ZERO);
    break;
  case 4 :
    return Dusty2D_pState(ONE,c()/rho,ZERO,c2()-(2.0/3.0)*k(),ZERO,ZERO);
    break;
  case 5 :
    return Dusty2D_pState(ZERO,ZERO,ZERO,-(2.0/3.0)*rho,ONE,ZERO);
    break;
  case 6 :
    return Dusty2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE);
    break;
  };
  return Dusty2D_pState(ONE,-c()/rho,ZERO,c2()-(2.0/3.0)*k(),ZERO,ZERO);
}

/**********************************************************************
 * Dusty2D_cState::rc_x -- Conserved right eigenvector (x-direction). *
 **********************************************************************/
inline Dusty2D_cState Dusty2D_cState::rc_x(int index) const {
  //assert(index >= 1 && index <= NUM_VAR_DUSTY2D);
  switch(index) {
  case 1 :
    return Dusty2D_cState(ONE,v().x-c(),v().y,h()-c()*v().x+(2.0/3.0)*k(),k(),omega());
    break;
  case 2 :
    return Dusty2D_cState(ONE,v().x,v().y,HALF*v().sqr()+gm1i*(g-5.0/3.0)*k(),k(),omega());
    break;
  case 3 :
    return Dusty2D_cState(ZERO,ZERO,rho,dv.y,ZERO,ZERO);
    break;
  case 4 :
    return Dusty2D_cState(ONE,v().x+c(),v().y,h()+c()*v().x+(2.0/3.0)*k(),k(),omega());
    break;
  case 5 :
    return Dusty2D_cState(ZERO,ZERO,ZERO,rho*gm1i*(g-5.0/3.0),rho,ZERO);
    break;
  case 6 :
    return Dusty2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,rho);
    break;
  };
  return Dusty2D_cState(ONE,v().x-c(),v().y,h()-c()*v().x+(2.0/3.0)*k(),k(),omega());
}

/**********************************************************************
 * Dusty2D_cState::lp_x -- Primitive left eigenvector (x-direction).  *
 **********************************************************************/
inline Dusty2D_pState Dusty2D_cState::lp_x(int index) const {
  //assert(index >= 1 && index <= NUM_VAR_DUSTY2D);
  switch(index) {
  case 1 :
    return Dusty2D_pState(k()/(3.0*c2()),-HALF*rho/c(),ZERO,HALF/c2(),rho/(3.0*c2()),ZERO);
    break;
  case 2 :
    return Dusty2D_pState(ONE-(2.0/3.0)*k()/c2(),ZERO,ZERO,-ONE/c2(),-(2.0/3.0)*rho/c2(),ZERO);
    break;
  case 3 :
    return Dusty2D_pState(ZERO,ZERO,ONE,ZERO,ZERO,ZERO);
    break;
  case 4 :
    return Dusty2D_pState(k()/(3.0*c2()),HALF*rho/c(),ZERO,HALF/c2(),rho/(3.0*c2()),ZERO);
    break;
  case 5 :
    return Dusty2D_pState(ZERO,ZERO,ZERO,ZERO,ONE,ZERO);
    break;
  case 6 :
    return Dusty2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE);
    break;
  };
  return Dusty2D_pState(k()/(3.0*c2()),-HALF*rho/c(),ZERO,HALF/c2(),rho/(3.0*c2()),ZERO);
}

/**********************************************************************
 * Dusty2DState -- External subroutines.                              *
 **********************************************************************/

extern Dusty2D_pState Riemann(const Dusty2D_pState &Wl,
	      	              const Dusty2D_pState &Wr);

extern Dusty2D_pState RoeAverage(const Dusty2D_pState &Wl,
	      	                 const Dusty2D_pState &Wr);

extern Dusty2D_pState Rotate(const Dusty2D_pState &W,
	      	             const Vector2D &norm_dir);

extern Dusty2D_cState Rotate(const Dusty2D_cState &U,
	      	             const Vector2D &norm_dir);

extern DenseMatrix RotationMatrix(const Dusty2D_pState &Wdum,
				  const Vector2D &norm_dir);

extern DenseMatrix iRotationMatrix(const Dusty2D_pState &Wdum,
				   const Vector2D &norm_dir);

extern Dusty2D_pState Translate(const Dusty2D_pState &W,
				const Vector2D &V);

extern Dusty2D_pState Reflect(const Dusty2D_pState &W,
	      	              const Vector2D &norm_dir);

extern Dusty2D_pState Absorb(const Dusty2D_pState &W,
  			     const Vector2D &norm_dir);

extern Dusty2D_pState WallViscousHeatFlux(const Dusty2D_pState &W,
					  const Vector2D &norm_dir);

extern Dusty2D_pState WallViscousIsothermal(const Dusty2D_pState &W,
					    const Vector2D &norm_dir,
					    const double &Twall);

extern Dusty2D_pState MovingWallHeatFlux(const Dusty2D_pState &W,
					 const Vector2D &norm_dir,
					 const double &Vwall);

extern Dusty2D_pState MovingWallIsothermal(const Dusty2D_pState &W,
					   const Vector2D &norm_dir,
					   const double &Vwall,
					   const double &Twall);

extern Dusty2D_pState BurningSurface(const Dusty2D_pState &W,
				     const Vector2D &norm_dir);

extern Dusty2D_pState RinglebFlow(const Dusty2D_pState &Wdum,
				  const Vector2D &X);

extern Dusty2D_pState RinglebFlow(const Dusty2D_pState &Wdum,
				  const Vector2D &X,
				  double &q, double &k);

extern Dusty2D_pState RinglebFlowAverageState(const Dusty2D_pState &Wdum,
					      const Vector2D &Y1,
					      const Vector2D &Y2,
					      const Vector2D &Y3,
					      const Vector2D &Y4);

extern Dusty2D_pState ViscousChannelFlow(const Dusty2D_pState &Wdum,
					 const Vector2D X,
					 const Vector2D Vwall,
					 const double dp,
					 const double length,
					 const double height);

extern Dusty2D_pState ViscousChannelFlowVelocity(const Dusty2D_pState &Wdum,
						 const Vector2D X,
						 const Vector2D Vwall,
						 const double dp,
						 const double length,
						 const double height);

extern void ViscousChannelFlowDerivatives(const Dusty2D_pState &W,
					  Dusty2D_pState &dWdx,
					  Dusty2D_pState &dWdy,
					  const Vector2D X,
					  const Vector2D Vwall,
					  const double dp,
					  const double length,
					  const double height);

extern Dusty2D_pState ViscousPipeFlow(const Dusty2D_pState &Wdum,
				      const Vector2D X,
				      const double dp,
				      const double length,
				      const double radius);

extern Dusty2D_pState TurbulentPipeFlow(const Dusty2D_pState &Wo,
					const Vector2D X,
					const double dp,
					const double length,
					const double radius);

extern Dusty2D_pState FlatPlate(const Dusty2D_pState &Winf,
				const Vector2D &X,
				const double &plate_length,
				double &eta,
				double &f,
				double &fp,
				double &fpp);

extern Dusty2D_pState DrivenCavityFlow(const Dusty2D_pState &Wo,
				       const double &l,
				       const double &Re);

extern Dusty2D_pState BackwardFacingStep(const Dusty2D_pState &Wo,
					 const Vector2D &X,
					 const double &h,
					 const double &ho,
					 const double &Re,
					 const double &M);

extern Dusty2D_pState BC_Characteristic(const Dusty2D_pState &Wi,
                                        const Dusty2D_pState &Wo,
	      	                        const Vector2D &norm_dir);

extern Dusty2D_pState BC_Characteristic_Pressure(const Dusty2D_pState &Wi,
                                                 const Dusty2D_pState &Wo,
	      	                                 const Vector2D &norm_dir);

extern Dusty2D_pState BC_Characteristic_Mach_Number(const Dusty2D_pState &Wi,
                                                    const Dusty2D_pState &Wo,
	      	                                    const Vector2D &norm_dir);

extern Dusty2D_pState WaveSpeedPos(const Dusty2D_pState &lambda_a,
                                   const Dusty2D_pState &lambda_l,
                                   const Dusty2D_pState &lambda_r);

extern Dusty2D_pState WaveSpeedNeg(const Dusty2D_pState &lambda_a,
                                   const Dusty2D_pState &lambda_l,
                                   const Dusty2D_pState &lambda_r);

extern Dusty2D_pState WaveSpeedAbs(const Dusty2D_pState &lambda_a,
                                   const Dusty2D_pState &lambda_l,
                                   const Dusty2D_pState &lambda_r);

extern Dusty2D_pState HartenFixPos(const Dusty2D_pState &lambda_a,
                                   const Dusty2D_pState &lambda_l,
                                   const Dusty2D_pState &lambda_r);

extern Dusty2D_pState HartenFixNeg(const Dusty2D_pState &lambda_a,
                                   const Dusty2D_pState &lambda_l,
                                   const Dusty2D_pState &lambda_r);

extern Dusty2D_pState HartenFixAbs(const Dusty2D_pState &lambda_a,
                                   const Dusty2D_pState &lambda_l,
                                   const Dusty2D_pState &lambda_r);

extern Dusty2D_cState FluxGodunov_n(const Dusty2D_pState &Wl,
	      	                    const Dusty2D_pState &Wr,
                                    const Vector2D &norm_dir);

extern Dusty2D_cState FluxGodunov_n(const Dusty2D_cState &Ul,
	      	                    const Dusty2D_cState &Ur,
                                    const Vector2D &norm_dir);

extern Dusty2D_cState FluxGodunov_MB_n(const Dusty2D_pState &Wl,
				       const Dusty2D_pState &Wr,
				       const Vector2D &V,
				       const Vector2D &norm_dir);

extern Dusty2D_cState FluxRoe(const Dusty2D_pState &Wl,
	      	              const Dusty2D_pState &Wr);

extern Dusty2D_cState FluxRoe(const Dusty2D_cState &Ul,
	      	              const Dusty2D_cState &Ur);

extern Dusty2D_cState FluxRoe_n(const Dusty2D_pState &Wl,
	      	                const Dusty2D_pState &Wr,
                                const Vector2D &norm_dir);

extern Dusty2D_cState FluxRoe_n(const Dusty2D_cState &Ul,
	      	                const Dusty2D_cState &Ur,
                                const Vector2D &norm_dir);

extern Dusty2D_cState FluxRoe_MB(const Dusty2D_pState &Wl,
				 const Dusty2D_pState &Wr,
				 const Vector2D &V);

extern Dusty2D_cState FluxRoe_MB(const Dusty2D_cState &Ul,
				 const Dusty2D_cState &Ur,
				 const Vector2D &V);

extern Dusty2D_cState FluxRoe_MB_n(const Dusty2D_pState &Wl,
				   const Dusty2D_pState &Wr,
				   const Vector2D &V,
				   const Vector2D &norm_dir);

extern Dusty2D_cState FluxRoe_MB_n(const Dusty2D_cState &Ul,
				   const Dusty2D_cState &Ur,
				   const Vector2D &V,
				   const Vector2D &norm_dir);

extern Dusty2D_cState FluxRusanov(const Dusty2D_pState &Wl,
	      	                  const Dusty2D_pState &Wr);

extern Dusty2D_cState FluxRusanov(const Dusty2D_cState &Ul,
	      	                  const Dusty2D_cState &Ur);

extern Dusty2D_cState FluxRusanov_n(const Dusty2D_pState &Wl,
	      	                    const Dusty2D_pState &Wr,
                                    const Vector2D &norm_dir);

extern Dusty2D_cState FluxRusanov_n(const Dusty2D_cState &Ul,
	      	                    const Dusty2D_cState &Ur,
                                    const Vector2D &norm_dir);

extern Dusty2D_cState FluxHLLE(const Dusty2D_pState &Wl,
	      	               const Dusty2D_pState &Wr);

extern Dusty2D_cState FluxHLLE(const Dusty2D_cState &Ul,
	      	               const Dusty2D_cState &Ur);

extern Dusty2D_cState FluxHLLE_n(const Dusty2D_pState &Wl,
	      	                 const Dusty2D_pState &Wr,
                                 const Vector2D &norm_dir);

extern Dusty2D_cState FluxHLLE_n(const Dusty2D_cState &Ul,
	      	                 const Dusty2D_cState &Ur,
                                 const Vector2D &norm_dir);

extern Dusty2D_cState FluxHLLE_MB(const Dusty2D_pState &Wl,
				  const Dusty2D_pState &Wr,
				  const Vector2D &V);

extern Dusty2D_cState FluxHLLE_MB(const Dusty2D_cState &Ul,
				  const Dusty2D_cState &Ur,
				  const Vector2D &V);

extern Dusty2D_cState FluxHLLE_MB_n(const Dusty2D_pState &Wl,
				    const Dusty2D_pState &Wr,
				    const Vector2D &V,
				    const Vector2D &norm_dir);

extern Dusty2D_cState FluxHLLE_MB_n(const Dusty2D_cState &Ul,
				    const Dusty2D_cState &Ur,
				    const Vector2D &V,
				    const Vector2D &norm_dir);

extern Vector2D HLLE_wavespeeds(const Dusty2D_pState &Wl,
                                const Dusty2D_pState &Wr,
                                const Vector2D &norm_dir);

extern Dusty2D_cState FluxHLLL(const Dusty2D_pState &Wl,
			       const Dusty2D_pState &Wr);

extern Dusty2D_cState FluxHLLL(const Dusty2D_cState &Ul,
			       const Dusty2D_cState &Ur);

extern Dusty2D_cState FluxHLLL_n(const Dusty2D_pState &Wl,
				 const Dusty2D_pState &Wr,
				 const Vector2D &norm_dir);

extern Dusty2D_cState FluxHLLL_n(const Dusty2D_cState &Ul,
				 const Dusty2D_cState &Ur,
				 const Vector2D &norm_dir);

extern Dusty2D_cState FluxHLLC(const Dusty2D_pState &Wl,
	      	               const Dusty2D_pState &Wr);

extern Dusty2D_cState FluxHLLC(const Dusty2D_cState &Ul,
	      	               const Dusty2D_cState &Ur);

extern Dusty2D_cState FluxHLLC_n(const Dusty2D_pState &Wl,
	      	                 const Dusty2D_pState &Wr,
                                 const Vector2D &norm_dir);

extern Dusty2D_cState FluxHLLC_n(const Dusty2D_cState &Ul,
	      	                 const Dusty2D_cState &Ur,
                                 const Vector2D &norm_dir);

extern Dusty2D_cState FluxVanLeer(const Dusty2D_pState &Wl,
				  const Dusty2D_pState &Wr);

extern Dusty2D_cState FluxVanLeer(const Dusty2D_cState &Wl,
				  const Dusty2D_cState &Wr);

extern Dusty2D_cState FluxVanLeer_n(const Dusty2D_pState &Wl,
				    const Dusty2D_pState &Wr,
				    const Vector2D &norm_dir);

extern Dusty2D_cState FluxVanLeer_n(const Dusty2D_cState &Wl,
				    const Dusty2D_cState &Wr,
				    const Vector2D &norm_dir);

extern Dusty2D_cState FluxVanLeer_MB(const Dusty2D_pState &Wl,
				     const Dusty2D_pState &Wr,
				     const Vector2D &V);

extern Dusty2D_cState FluxVanLeer_MB_n(const Dusty2D_pState &Wl,
				       const Dusty2D_pState &Wr,
				       const Vector2D &V,
				       const Vector2D &norm_dir);

extern Dusty2D_cState FluxAUSM(const Dusty2D_pState &Wl,
			       const Dusty2D_pState &Wr);

extern Dusty2D_cState FluxAUSM(const Dusty2D_cState &Wl,
			       const Dusty2D_cState &Wr);

extern Dusty2D_cState FluxAUSM_n(const Dusty2D_pState &Wl,
				 const Dusty2D_pState &Wr,
				 const Vector2D &norm_dir);

extern Dusty2D_cState FluxAUSM_n(const Dusty2D_cState &Wl,
				 const Dusty2D_cState &Wr,
				 const Vector2D &norm_dir);

extern Dusty2D_cState FluxAUSMplus(const Dusty2D_pState &Wl,
				   const Dusty2D_pState &Wr);

extern Dusty2D_cState FluxAUSMplus(const Dusty2D_cState &Wl,
				   const Dusty2D_cState &Wr);

extern Dusty2D_cState FluxAUSMplus_n(const Dusty2D_pState &Wl,
				     const Dusty2D_pState &Wr,
				     const Vector2D &norm_dir);

extern Dusty2D_cState FluxAUSMplus_n(const Dusty2D_cState &Wl,
				     const Dusty2D_cState &Wr,
				     const Vector2D &norm_dir);

extern Dusty2D_cState FluxSaurel_n(const Dusty2D_pState &Wl,
				   const Dusty2D_pState &Wr,
				   const Vector2D &norm_dir);

extern Dusty2D_cState FluxSaurel_MB_n(const Dusty2D_pState &Wl,
				      const Dusty2D_pState &Wr,
				      const Vector2D &V,
				      const Vector2D &norm_dir);

extern Dusty2D_cState FluxMultiVelocity_n(const Dusty2D_pState &Wl,
					  const Dusty2D_pState &Wr,
					  const Vector2D &norm_dir);

extern Dusty2D_cState FluxMultiVelocity_MB_n(const Dusty2D_pState &Wl,
					     const Dusty2D_pState &Wr,
					     const Vector2D &V,
					     const Vector2D &norm_dir);

extern Dusty2D_cState ViscousFlux_n(const Vector2D &X,
				    Dusty2D_pState &W,
				    const Dusty2D_pState &dWdx,
				    const Dusty2D_pState &dWdy,
				    const Vector2D &norm_dir,
				    const int &Axisymmetric,
				    const int &adiabatic_flag);

extern Dusty2D_cState ViscousFluxDiamondPath_n(const Vector2D &X,
					       const Vector2D &Xl, const Dusty2D_pState &Wl,
					       const Vector2D &Xd, const Dusty2D_pState &Wd,
					       const Vector2D &Xr, const Dusty2D_pState &Wr,
					       const Vector2D &Xu, const Dusty2D_pState &Wu,
					       const Vector2D &norm_dir,
					       const int &Axisymmetric,
					       const int &stencil_flag);

extern Dusty2D_cState ViscousFluxHybrid_n(const Vector2D &X,
					  Dusty2D_pState &W,
					  const Vector2D &X1,
					  const Dusty2D_pState &W1,
					  const Dusty2D_pState &dW1dx,
					  const Dusty2D_pState &dW1dy,
					  const Vector2D &X2,
					  const Dusty2D_pState &W2,
					  const Dusty2D_pState &dW2dx,
					  const Dusty2D_pState &dW2dy,
					  const Vector2D &norm_dir,
					  const int &Axisymmetric,
					  const int &adiabatic_flag);

extern double ShearStress(const Dusty2D_pState &W,
			  const Dusty2D_pState &dWdx,
			  const Dusty2D_pState &dWdy,
			  const Vector2D &nhat);

extern double WallShearStress(const Dusty2D_pState &W1,
			      const Vector2D &X1,
			      const Vector2D &X2,
			      const Vector2D &X3,
			      const Vector2D &nhat);

#endif // _DUSTY2D_STATE_INCLUDED
