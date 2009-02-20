/**********************************************************************
 * NavierStokes2DState.h: Header file defining 2D NavierStokes        *
 *                        solution state classes.                     *
 **********************************************************************/

#ifndef _NAVIERSTOKES2D_STATE_INCLUDED
#define _NAVIERSTOKES2D_STATE_INCLUDED

// Include required C++ libraries.
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

/* Using std namespace functions */
using namespace std;

/* Include CFFC header files */
#include "../Math/Math.h"  /* Include math macro header files. */
#include "../CFD/CFD.h"    /* Include CFD header files. */
#include "../Math/Matrix.h"  /* Include matrix header files. */
#include "../Math/Vector2D.h" /* Include vector 2D header files. */
#include "../Math/Tensor2D.h" /* Include tensor 2D header files. */
#include "../Physics/GasConstants.h" /* Include gas constant header files. */
#include "../Physics/SolidConstants.h" /* Include solid constant header files. */
#include "../Utilities/Utilities.h"

// Define the number of variables.
#define	NUM_VAR_NAVIERSTOKES2D  8

// Define the classes.

class NavierStokes2D_cState;

/*!
 * Class: NavierStokes2D_pState
 *
 * @brief Primitive variable solution state class definition for a 
 *        laminar or turbulent compressible gas-flow.
 *
 * Primitive variable solution state class definition for a laminar or
 * turbulent (k-omega) gas-flow.
 *
 * \verbatim
 * Member functions
 *     rho      -- Gas density.
 *     v        -- Gas velocity.
 *     p        -- Gas pressure.
 *     k        -- Gas turbulent kinetic energy.
 *     omega    -- Gas specific dissipation rate.
 *     ke       -- Gas internal energy variance.
 *     ee       -- Gas energy variance dissipation rate.
 *     tau      -- Viscous stress tensor (laminar and Reynolds).
 *     q        -- Heat flux vector (laminar and turbulent).
 *     g        -- Gas specific heat ratio.
 *     gm1      -- Return g-1.
 *     gm1i     -- Return 1/(g-1).
 *     R        -- Gas constant.
 *     v1,v2,v3,v4,v5 -- Viscosity law coefficients.
 *     cp       -- Specific heat at constant pressure.
 *     cv       -- Specific heat at constant volume.
 *     T        -- Gas temperature.
 *     e        -- Gas specific internal energy.
 *     E        -- Gas total energy.
 *     h        -- Gas specific enthalpy.
 *     H        -- Gas total enthalpy.
 *     a        -- Gas sound speed.
 *     a2       -- Gas sound speed squared.
 *     M        -- Gas Mach number.
 *     s        -- Gas specific entropy.
 *     dv       -- Gas momentum.
 *     To       -- Gas stagnation temperature.
 *     po       -- Gas stagnation pressure.
 *     ao       -- Gas stagnation sound speed.
 *     ho       -- Gas stagnation enthalpy.
 *     mu       -- Gas dynamic viscosity.
 *     nu       -- Gas kinematic viscosity.
 *     kappa    -- Gas thermal conductivity.
 *     Pr       -- Prandtl number.
 *     meanfreepath -- Gas mean free path.
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
 *     sigmav   -- Surface roughness parameter.
 *     lw       -- Length scale of the turbulence due to mass injection.
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
 *     dFdW     -- Return the Jacobian of the inviscid solution flux
 *                 vector with respect to the primitive solution
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
 *     St       -- Return turbulent source term vector.
 *     dStdU    -- Return the Jacobian of the turbulent source term
 *                 vector with respect to the conserved solution 
 *                 variables.
 *
 *     dSvpdU   -- Return the Jacobian of the source terms of the two
 *                 additional equations for the variable turbulent Prandtl 
 *                 number with respect to the conserved solution variables.
 *
 *    productionK -- Return production term for the Ke PDE.
 *    diff      -- difference between Major production and destruction terms.
 *    deriv2    -- derivative of internal energy, squared
 *    productionE1 -- Return production term for the ee PDE.
 *    productionE1_1 -- Return the first part of production term E1 for the ee PDE.
 *    productionE1_2 -- Return the second part of production term E1 for the ee PDE.
 *    productionE2 -- Return production term for the ee PDE.
 *    D1         -- First destruction term coefficient in ee PDE.
 *    D2         -- Second destruction term coefficient in ee PDE.
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
class NavierStokes2D_pState {
 private:
 public:
  //@{ @name Primitive variables and associated constants:
  double                       rho; //!< Gas density.
  Vector2D                       v; //!< Gas velocity (2D vector).
  double                         p; //!< Gas pressure.
  double                         k; //!< Gas turbulent kinetic energy.
  double                     omega; //!< Gas specific dissipation rate.
  double                        ke; //!< Gas internal energy variance.
  double                        ee; //!< Gas energy variance dissipation rate.
  Tensor2D                     tau; //!< Viscous stess tensor (laminar and turbulent).
  Vector2D                       q; //!< Heat flux vector (laminar and turbulent).
  static double                  g; //!< Specific heat ratio.
  static double                gm1; //!< g-1
  static double               gm1i; //!< 1/(g-1)
  static double                  R; //!< Gas constant.
  static double                 cp; //!< Specific heat at constant pressure.
  static double                 cv; //!< Specific heat at constant volume.
  static double v1, v2, v3, v4, v5; //!< Viscosity law coefficients.
  static int             flow_type; //!< Flow-type indicator (inviscid, laminar, or k-omega turbulent).
  //@}

  //@{ @name Turbulence boundary-layer constants:
  static int      Transition_Model; //!< Transition model indicator (Off, Wilcox, Menter).
  static int Compressibility_Effect;//!< Compressibility Correction (Off, Sarkar, Wilcox, Zeman).
  static int      Variable_Prandtl; //! Variable Prandtl number indicatot (Off, On)
  static double            yplus_o; //!< Transition between viscous sublayer and log layer.
  static double                  C; //!< Surface roughness coefficient.
  static double         von_karman; //!< Von Karman constant.
  static double     yplus_sublayer; //!< Sublayer dimensionless wall distance.
  static double yplus_buffer_layer; //!< Buffer layer dimensionless wall distance.
  static double  yplus_outer_layer; //!< Outer layer dimensionless wall distance.
  //@}

  //@{ @name k-omega closure coefficients:
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
  static double             sigmav; //!< Surface roughness parameter.
  static double                 lw; //!< Length scale of the turbulence due to mass injection.
  //@}

  //@{ @Wilcox Transition model variables:
  static double             Rbeta;
  static double                Rk;
  static double            Romega;     
  static double           alpha_o; 
  static double        sigma_star;
  static double      sigma_Wilcox; 
  static double       beta_Wilcox;
  //@}

  //@{Variable Prandtl number coefficients
  static double               Cd1;
  static double               Cd2;
  static double               Cd3;
  static double               Cd4;
  static double               Cd5;
  static double         sigma_k_e;
  static double        sigma_ep_e;
  static double            A_plus;
  static double         C1_lambda;
  static double          C_lambda;

  //@}

  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  NavierStokes2D_pState(void) {
    Standard_Atmosphere();
  }

  //! Copy constructor.
  NavierStokes2D_pState(const NavierStokes2D_pState &W) {
    rho = W.rho; v.x = W.v.x; v.y = W.v.y; p = W.p; k = W.k; omega = W.omega; ke = W.ke; ee = W.ee;
  }

  //! Copy constructor.
  NavierStokes2D_pState(const NavierStokes2D_cState &U);

  //! Assignment constructor.
  NavierStokes2D_pState(const double &dens,
			const Vector2D &V,
			const double &pre) {
    rho = dens; v.x = V.x; v.y = V.y; p = pre; k = ZERO; omega = ZERO; ke = ZERO; ee = ZERO;
  }
  
  //! Assignment constructor.
  NavierStokes2D_pState(const double &dens,
			const double &vx,
			const double &vy,
			const double &pre) {
    rho = dens; v.x = vx; v.y = vy; p = pre; k = ZERO; omega = ZERO; ke = ZERO; ee = ZERO;
  }

  //! Assignment constructor.
  NavierStokes2D_pState(const double &dens,
			const Vector2D &V,
			const double &pre,
			const double &kk,
			const double &omga) {
    rho = dens; v.x = V.x; v.y = V.y; p = pre; k = kk; omega = omga; ke = ZERO; ee = ZERO;
  }

  //! Assignment constructor.
  NavierStokes2D_pState(const double &dens,
			const double &vx,
			const double &vy,
			const double &pre,
			const double &kk,
			const double &omga) {
    rho = dens; v.x = vx; v.y = vy; p = pre; k = kk; omega = omga; ke = ZERO; ee = ZERO;
  }

  //! Assignment constructor.
  NavierStokes2D_pState(const double &dens,
			const Vector2D &V,
			const double &pre,
			const double &kk,
			const double &omga,
			const double &kee,
			const double &epsi_e) {
    rho = dens; v.x = V.x; v.y = V.y; p = pre; k = kk; omega = omga; ke = kee; ee = epsi_e;
  }

  //! Assignment constructor.
  NavierStokes2D_pState(const double &dens,
			const double &vx,
			const double &vy,
			const double &pre,
			const double &kk,
			const double &omga,
			const double &kee,
			const double &epsi_e) {
    rho = dens; v.x = vx; v.y = vy; p = pre; k = kk; omega = omga; ke = kee; ee = epsi_e;
  }

  //! Value Constructor
  explicit NavierStokes2D_pState(const double &Val);

  //! Destructor.
  ~NavierStokes2D_pState(void) { }
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
			    char *propellant_type,
			    const int &i_Compressibility_Effect,
			    const int &i_Transition_Model,
			    const int &i_Variable_Prandtl,
			    const double &sigv,
			    const double &lengthw,
			    const double &C_lamb,
			    const double &Cd_1,
			    const double &Cd_4,
			    const double &Cd_5);
  void set_gas(char *gas_type);
  void set_turbulence(const double &C_constant,
		      const double &von_karman,
		      const double &yplus_sub,
		      const double &yplus_buffer,
		      const double &yplus_outer,
                      const int &i_Compressibility_Effect,
		      const int &i_Transition_Model,
		      const int &i_Variable_Prandtl,
		      const double &C_lamb,
		      const double &Cd_1,
		      const double &Cd_4,
		      const double &Cd_5);
  void set_propellant(char *propellant_type,
		      const double &sigv,
		      const double &lengthw);
  //@}

  //@{ @name Useful operators.
  //! Return the number of variables.
  static int NumVar(void) { return NUM_VAR_NAVIERSTOKES2D; }

  //! Number of active solution state variables (i.e. it depends on the flow type)
  static int NumVarActive(void);

  //! Copy operator.
  void Copy(const NavierStokes2D_pState &W);

  //! Vacuum operator.
  void Vacuum(void);

  //! Standard atmosphere operator.
  void Standard_Atmosphere(void);

  //! One operator. Set the solution to ONE.
  void One(void);

  //! Check for unphysical state properties.
  int Unphysical_Properties(void) const;
  //@}

  //@{ @name State functions.
  //! Gas temperature.
  double T(void) const;
  
  //! Gas specific internal energy.
  double e(void) const;
  
  //! Gas total energy.
  double E(void) const;
  
  //! Gas specific enthalpy.
  double h(void) const;
  
  //! Gas total enthalpy.
  double H(void) const;

  //! Gas sound speed.
  double a(void) const;

  //! Gas sound speed squared.
  double a2(void) const;

  //! Gas Mach number.
  double M(void) const;

  //! Gas specific entropy.
  double s(void) const;

  //! Gas momentum.
  Vector2D dv(void) const;

  //! Gas momentum.
  double dv(const Vector2D &n) const;
  
  //! Gas stagnation temperature.
  double To(void) const;

  //! Gas stagnation pressure.
  double po(void) const;

  //! Gas stagnation sound speed.
  double ao(void) const;

  //! Gas stagnation enthalpy.
  double ho(void) const;

  //! Gas dynamic viscosity.
  double mu(void) const;

  //! Gas kinematic viscosity.
  double nu(void) const;

  //! Gas thermal heat conductivity.
  double kappa(void) const;

  //! Prandtl number.
  double Pr(void) const;

  //! Thermal diffusivity.
  double Alpha(void) const;

  //! Gas mean free path.
  double meanfreepath(void) const;
  //@}

  //@{ @name Turbulence related functions.
  //! Return the total turbulent kinetic energy.
  double dk(void) const;

  //! Return the total turbulent specific dissipation.
  double domega(void) const;

  //! Return the total internal energy variance.
  double dke(void) const;

  //! Return the internal energy variance rate.
  double dee(void) const;

  //! Return the turbulent eddy dissipation.
  double epsilon(void) const;

  //! Return the turbulent eddy dissipation.
  double depsilon(void) const;

  //! Return the turbulent length scale.
  double ell(void) const;

  //! Return the turbulent Mach number.
  double Mt(void) const;

  //! Return the turbulent Mach number squared.
  double Mt2(void) const;

  //! Return the turbulent Reynolds number
  double ReT(void) const;

  //! Return the alpha0 star coefficient(Wilcox transition model)
  double alpha_o_star(void) const;

  //! Return the beta star coefficient(Wilcox transition model)
  double beta_star(void) const;

  //! Return the alpha coefficient(Wilcox transition model)
  double alpha_Wilcox(void) const;

  //! Return the alpha star coefficient(Wilcox transition model)
  double alpha_star_Wilcox(void) const;

  //! Turbulent eddy dynamic viscosity.
  double muT(void) const;

  //! Turbulent Prandtl number.
  double PrT(const double &ywall,const double &yplus) const;
  
  //! Turbulent thermal diffusivity
  double alphaT(const double &ywall,const double &yplus) const;

  //! Turbulent eddy kinematic viscosity.
  double nuT(void) const;

  //! Turbulent eddy thermal heat conductivity.
  double kappaT(const double &ywall,const double &yplus) const;

  //! Turbulence modified sound speed.
  double c(void) const;

  //! Turbulence modified sound speed squared.
  double c2(void) const;

  //! Turbulence modified pressure.
  double pmodified(void) const;

  //! Kolmogorov velocity scale.  
  double U_e(void) const;

  //! Near Wall damping function (variable Prandtl number model).
  double f_lambda(const double &ywall)const;

  //! Near Wall damping function (variable Prandtl number model).
  double f_mu(const double &ywall)const;
  
  //! Near Wall damping function (variable Prandtl number model).
  double xi_et(const double &ywall, const NavierStokes2D_pState &dWdx, const NavierStokes2D_pState &dWdy)const;
  
  //! Beta_star value for the k equation.
  double beta_k(const NavierStokes2D_pState &dWdx,const NavierStokes2D_pState &dWdy) const;
  
  //! Beta value for the omega equation.
  double beta_omega(const NavierStokes2D_pState &dWdx, const NavierStokes2D_pState &dWdy) const;
 
  //! F_beta_star function for the k equation.
  double f_beta_k(const NavierStokes2D_pState &dWdx, const NavierStokes2D_pState &dWdy) const;

  //! F_beta function for the omega equation.
  double f_beta_omega(const NavierStokes2D_pState &dWdx, const NavierStokes2D_pState &dWdy) const;

  //! Sarkar dilatation dissipation correction function.
  double F_Sarkar(void) const;

  //! Zeman dilatation dissipation correction function.
  double F_Zeman(void) const;

  //! Wilcox dilatation dissipation correction function.
  double F_Wilcox(void) const;
  
  //! Chi_k in the f_beta_star function in k equation.
  double chi_k(const NavierStokes2D_pState &dWdx, const NavierStokes2D_pState &dWdy) const;

  //! Chi_Omega in the f_beta function in the omega equation.
  double chi_omega(const NavierStokes2D_pState &dWdx, const NavierStokes2D_pState &dWdy) const;
  //@}

  //@{ @name Solid propellant related functions.
  //! Burning rate.
  double burningrate(void) const;
  //@}

  //@{ @name Conserved solution state.
  NavierStokes2D_cState U(void) const;
  NavierStokes2D_cState U(const NavierStokes2D_pState &W) const;
  friend NavierStokes2D_cState U(const NavierStokes2D_pState &W);
  //@}

  //! Jacobian of the conserved solution variables with respect to the
  //! primitive solution variables.
  void dUdW(DenseMatrix &dUdW) const;

  //! Jacobian of the primitive solution variables with respect to the
  //! conserved solution variables.
  void dWdU(DenseMatrix &dWdU) const;

  //@{ @name Inviscid solution flux (x-direction) and Jacobian.
  NavierStokes2D_cState F(void) const;
  NavierStokes2D_cState F(const Vector2D &V) const;
  void dFdU(DenseMatrix &dFdU) const;
  void dFdW(DenseMatrix &dFdW) const;
  //! @brief Inviscid solution flux in the normal direction
  NavierStokes2D_cState Fn(const Vector2D & normal_dir) const;
  //! @brief Calculate flux in the provided normal direction for a given solution state
  friend NavierStokes2D_cState Fn(const NavierStokes2D_pState &W, const Vector2D & normal_dir);
  //@}

  //@{ @name Viscous solution fluxes and Jacobians.
  NavierStokes2D_cState Gx(const NavierStokes2D_pState &dWdx,const double &ywall,const double &yplus) const;
  NavierStokes2D_cState Gy(const NavierStokes2D_pState &dWdy,const double &ywall,const double &yplus) const;
  //@}

  //! Compute viscous stress tensor and heat flux vector.
  void ComputeViscousTerms(const NavierStokes2D_pState &dWdx,
			   const NavierStokes2D_pState &dWdy,
			   const Vector2D &X,
			   const int &Axisymmetric,
			   const int &adiabatic_flag,
			   const double &ywall,
			   const double &yplus);

  //@{ @name Eigenstructure.
  //! Eigenvalue(s) (x-direction).
  NavierStokes2D_pState lambda_x(void) const;

  //! Eigenvalue(s) (x-direction).
  NavierStokes2D_pState lambda_x(const Vector2D &V) const;

  //! Primitive right eigenvector (x-direction).
  NavierStokes2D_pState rp_x(int index) const;

  //! Conserved right eigenvector (x-direction).
  NavierStokes2D_cState rc_x(int index) const;

  //! Primitive left eigenvector (x-direction).
  NavierStokes2D_pState lp_x(int index) const;
  //@}

  //@{ @name Include all source vectors and Jacobians.
  NavierStokes2D_cState S(const Vector2D &X,
                          const NavierStokes2D_pState &W,
			  const NavierStokes2D_pState &dWdx,
			  const NavierStokes2D_pState &dWdy,
			  const int &Axisymmetric,
			  const double &ywall,
		          const double &yplus) const;

  void dSdU(DenseMatrix &dSdU,
	    const Vector2D &X,
	    const NavierStokes2D_pState &dWdx,
	    const NavierStokes2D_pState &dWdy,
	    const int &Axisymmetric) const;
  //@}

  //@{ @name Inviscid axisymmetric flow source vector and Jacobian.
  NavierStokes2D_cState Si(const Vector2D &X) const;
  void dSidU(DenseMatrix &dSidU, const Vector2D &X) const;
  //@}

  //@{ @name Viscous axisymmetric flow source vector and Jacobian.
  NavierStokes2D_cState Sv(const Vector2D &X,
			   const NavierStokes2D_pState &dWdy,
			   const double &ywall,
		           const double &yplus) const;
  void dSvdU(DenseMatrix &dSvdU, const Vector2D &X, const NavierStokes2D_pState &dWdy) const;
  //@}

  //@{ @name Turbulent source term vector and Jacobian.
  NavierStokes2D_cState St(const Vector2D &X,
                           const NavierStokes2D_pState &W,
			   const NavierStokes2D_pState &dWdx,
			   const NavierStokes2D_pState &dWdy,
			   const int &Axisymmetric,
			   const double &ywall,
			   const double &yplus) const;

  void dStdU(DenseMatrix &dStdU,
	     const Vector2D &X,
	     const NavierStokes2D_pState &dWdx,
	     const NavierStokes2D_pState &dWdy,
	     const int &Axisymmetric) const;

  void dSvpdU(DenseMatrix &dSvpdU,
	      const Vector2D &X,
	      const NavierStokes2D_pState &dWdx,
	      const NavierStokes2D_pState &dWdy,
	      const double &d_dWdx_dW, 
	      const double &d_dWdy_dW,
	      const int &Axisymmetric,
	      const double &ywall,
	      const double &yplus) const;
  //@}

  double deriv2(const NavierStokes2D_pState &dWdx,
		const NavierStokes2D_pState &dWdy) const;
  
  double diff(const NavierStokes2D_pState &dWdx,
	      const NavierStokes2D_pState &dWdy,
	      const double &ywall,
	      const double &yplus) const;
  
  double productionE1(const NavierStokes2D_pState &dWdx,
		     const NavierStokes2D_pState &dWdy,
		     const double &ywall,
		     const double &yplus) const;

  double productionE1_1(const NavierStokes2D_pState &dWdx,
			const NavierStokes2D_pState &dWdy,
			const double &ywall,
			const double &yplus) const;
  
  double productionE1_2(const NavierStokes2D_pState &dWdx,
			const NavierStokes2D_pState &dWdy,
			const double &ywall,
			const double &yplus) const;
  
  double productionE2(const double &production_cap) const;
  double D1(void) const;
  double D2(void) const;

  double productionK(const NavierStokes2D_pState &dWdx,
		     const NavierStokes2D_pState &dWdy,
		     const double &ywall,
		     const double &yplus) const;
  
  //@{ @name Index operator.
  double &operator[](int index) {
    //assert(index >= 1 && index <= NUM_VAR_NAVIERSTOKES2D);
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
    case 7 :
      return ke;
    case 8 :
      return ee;
    };
    // Default return, this is never reached.
    return rho;
  }
  
  const double &operator[](int index) const {
    //assert(index >= 1 && index <= NUM_VAR_NAVIERSTOKES2D);
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
    case 7 :
      return ke;
    case 8 :
      return ee;
    };
    // Default return, this is never reached.
    return rho;
  }
  //@}

  //@{ @name Binary arithmetic operators.
  NavierStokes2D_pState operator +(const NavierStokes2D_pState &W) const;
  NavierStokes2D_pState operator -(const NavierStokes2D_pState &W) const;
  double operator *(const NavierStokes2D_pState &W) const;
  NavierStokes2D_pState operator *(const double &a) const;
  friend NavierStokes2D_pState operator *(const double &a, const NavierStokes2D_pState &W){ return W*a; }
  NavierStokes2D_pState operator /(const double &a) const;
  NavierStokes2D_pState operator /(const NavierStokes2D_pState &W);
  NavierStokes2D_pState operator ^(const NavierStokes2D_pState &W) const;
  friend NavierStokes2D_pState max(const NavierStokes2D_pState &W1, const NavierStokes2D_pState &W2 );
  friend NavierStokes2D_pState min(const NavierStokes2D_pState &W1, const NavierStokes2D_pState &W2 );
  //@}

  //@{ @name Assignment operator.
  NavierStokes2D_pState& operator =(const NavierStokes2D_pState &W);
  //@}

  //@{ @name Unary arithmetic operators.
  friend NavierStokes2D_pState operator +(const NavierStokes2D_pState &W);
  friend NavierStokes2D_pState operator -(const NavierStokes2D_pState &W);
  friend NavierStokes2D_pState fabs(const NavierStokes2D_pState &W);
  friend NavierStokes2D_pState sqr(const NavierStokes2D_pState &W);
  //@}

  //@{ @name Shortcut arithmetic operators.
  NavierStokes2D_pState &operator +=(const NavierStokes2D_pState &W);
  NavierStokes2D_pState &operator -=(const NavierStokes2D_pState &W);
  NavierStokes2D_pState &operator /=(const NavierStokes2D_pState &W);
  NavierStokes2D_pState &operator *=(const NavierStokes2D_pState &W);
  NavierStokes2D_pState &operator *=(const double &a);
  NavierStokes2D_pState &operator /=(const double &a);
  //@}

  //@{ @name Relational operators.
  friend int operator ==(const NavierStokes2D_pState &W1, const NavierStokes2D_pState &W2);
  friend int operator !=(const NavierStokes2D_pState &W1, const NavierStokes2D_pState &W2);
  friend bool operator >=(const NavierStokes2D_pState& W1, const NavierStokes2D_pState& W2);
  friend bool operator <=(const NavierStokes2D_pState& W1, const NavierStokes2D_pState& W2);
  friend bool operator <(const NavierStokes2D_pState &W1, const NavierStokes2D_pState &W2);
  friend bool operator >(const NavierStokes2D_pState &W1, const NavierStokes2D_pState &W2);
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const NavierStokes2D_pState &W);
  friend istream &operator >> (istream &in_file, NavierStokes2D_pState &W);
  //@}

  //@{ @name Output functions.
  void output_labels(ostream &out_file);

  void output_data(ostream &out_file, const double &ywall,const double &yplus);
  //@}

  int analytically_inverted_relaxation() { //this is needed for embeddedboundaries with gaussian2D
    return 0;
  }
  void relax(double deltat, int stage, const NavierStokes2D_pState &W) {
    //this is needed for embeddedboundaries with gaussian2D
    return;
  }

  double pressure() const {return p;} //added for compatibility with embeddedboundaries2D

};

/*!
 * Class: NavierStokes2D_cState
 *
 * @brief Conserved variable solution state class definition for a 
 *        laminar or turbulent compressible gas-flow.
 *
 * Conserved variable solution state class definition for a laminar or
 * turbulent (k-omega) gas-flow.
 *
 * \verbatim
 * Member functions
 *     rho      -- Gas density.
 *     dv       -- Gas momentum.
 *     E        -- Gas total energy.
 *     dk       -- Gas total turbulent kinetic energy.
 *     domega   -- Gas total specific dissipation rate.
 *     dke      -- Gas total internal energy variance.
 *     dee      -- Gas total internal energy variance dissipation rate. 
 *     tau      -- Viscous stress tensor (laminar and Reynolds).
 *     q        -- Heat flux vector (laminar and turbulent).
 *     g        -- Gas specific heat ratio.
 *     gm1      -- Return g-1.
 *     gm1i     -- Return 1/(g-1).
 *     R        -- Gas gas constant.
 *     cp       -- Specific heat at constant pressure.
 *     cv       -- Specific heat at constant volume.
 *     v1,v2,v3,v4,v5 -- Viscosity law coefficients.
 *     v        -- Gas flow velocity.
 *     p        -- Gas pressure.
 *     T        -- Gas temperature.
 *     e        -- Gas specific internal energy.
 *     h        -- Gas specific enthalpy.
 *     H        -- Gas total enthalpy.
 *     a        -- Gas sound speed.
 *     a2       -- Gas sound speed squared.
 *     M        -- Gas Mach number.
 *     s        -- Gas specific entropy.
 *     To       -- Gas stagnation temperature.
 *     po       -- Gas stagnation pressure.
 *     ao       -- Gas stagnation sound speed.
 *     ho       -- Gas stagnation enthalpy.
 *     mu       -- Gas dynamic viscosity.
 *     nu       -- Gas kinematic viscosity.
 *     kappa    -- Gas thermal conductivity.
 *     Pr       -- Prandtl number.
 *     meanfreepath -- Gas mean free path.
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
 *     sigmav   -- Surface roughness parameter.
 *     lw       -- Length scale of the turbulence due to mass injection.
 *
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
class NavierStokes2D_cState {
  private:
  public:
  //@{ @name Gas conservative variables:
  double                       rho; //!< Gas density.
  Vector2D                      dv; //!< Gas momentum.
  double                         E; //!< Gas total energy.
  double                        dk; //!< Gas total turbulent kinetic energy.
  double                    domega; //!< Gas total turbulent specific dissipation rate.
  double                   dke; //!< Gas total internal energy variance.
  double                   dee; //!< Gas total internal energy variance dissipation rate.
  Tensor2D                     tau; //!< Viscous stess tensor (laminar and turbulent).
  Vector2D                       q; //!< Heat flux vector (laminar and turbulent).
  static double                  g; //!< Specific heat ratio.
  static double                gm1; //!< g-1
  static double               gm1i; //!< 1/(g-1)
  static double                  R; //!< Gas constant.
  static double                 cp; //!< Specific heat at constant pressure.
  static double                 cv; //!< Specific heat at constant volume.
  static double v1, v2, v3, v4, v5; //!< Viscosity law coefficients.
  static int             flow_type; //!< Flow-type indicator (inviscid, laminar, or k-omega turbulent).
  //@}

  //@{ @name Turbulent boundary-layer constants:
  static int      Transition_Model; //!< Transition model indicator (Off, Wilcox, Menter).
  static int Compressibility_Effect;//!< Compressibility Correction (Off, Wilcox, Zeman, Sarkar).
  static int      Variable_Prandtl; //! Variable Prandtl number indicatot (Off, On).
  static double            yplus_o; //!< Transition between viscous sublayer and log layer.
  static double                  C; //!< Surface roughness coefficient.
  static double         von_karman; //!< Von Karman constant.
  static double     yplus_sublayer; //!< Sublayer dimensionless wall distance.
  static double yplus_buffer_layer; //!< Buffer layer dimensionless wall distance.
  static double  yplus_outer_layer; //!< Outer layer dimensionless wall distance.
  //@}

  //@{ @name k-omega closure coefficients:
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
  static double             sigmav; //!< Surface roughness parameter.
  static double                 lw; //!< Length scale of the turbulence due to mass injection.
  //@}

  //@{ @Wilcox Transition model variables:
  static double             Rbeta;
  static double                Rk;
  static double            Romega;
  static double           alpha_o; 
  static double        sigma_star;
  static double      sigma_Wilcox; 
  static double       beta_Wilcox;     
  //@}

  //@{Variable Prandtl number coefficients
  static double               Cd1;
  static double               Cd2;
  static double               Cd3;
  static double               Cd4;
  static double               Cd5;
  static double         sigma_k_e;
  static double        sigma_ep_e;
  static double            A_plus;
  static double         C1_lambda;
  static double          C_lambda;
  //@}

  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  NavierStokes2D_cState(void) {
    Standard_Atmosphere();
  }

  //! Copy constructor.
  NavierStokes2D_cState(const NavierStokes2D_cState &U) {
    rho = U.rho; dv.x = U.dv.x; dv.y = U.dv.y; E = U.E; dk = U.dk; domega = U.domega; dke = U.dke; dee = U.dee;
  }

  //! Copy constructor.
  NavierStokes2D_cState(const NavierStokes2D_pState &W);

  //! Value Constructor
  explicit NavierStokes2D_cState(const double &Val);

  //! Assignment constructor.
  NavierStokes2D_cState(const double &dens,
			const Vector2D &dV,
			const double &Etot) {
    rho = dens; dv.x = dV.x; dv.y = dV.y; E = Etot; dk = ZERO; domega = ZERO; dke = ZERO; dee = ZERO;
  }

  //! Assignment constructor.
  NavierStokes2D_cState(const double &dens,
			const double &dvx,
			const double &dvy,
			const double &Etot) {
    rho = dens; dv.x = dvx; dv.y = dvy; E = Etot; dk = ZERO; domega = ZERO; dke = ZERO; dee = ZERO;
  }

  //! Assignment constructor.
  NavierStokes2D_cState(const double &dens,
			const Vector2D &dV,
			const double &Etot,
			const double &dkdk,
			const double &domga) {
    rho = dens; dv.x = dV.x; dv.y = dV.y; E = Etot; dk = dkdk; domega = domga; dke = ZERO; dee = ZERO;
  }

  //! Assignment constructor.
  NavierStokes2D_cState(const double &dens,
			const double &dvx,
			const double &dvy,
			const double &Etot,
			const double &dkdk,
			const double &domga) {
    rho = dens; dv.x = dvx; dv.y = dvy; E = Etot; dk = dkdk; domega = domga; dke = ZERO; dee = ZERO;
  }

  //! Assignment constructor.
  NavierStokes2D_cState(const double &dens,
			const Vector2D &dV,
			const double &Etot,
			const double &dkdk,
			const double &domga,
			const double &dkedke,
			const double &deedee) {
    rho = dens; dv.x = dV.x; dv.y = dV.y; E = Etot; dk = dkdk; domega = domga; dke = dkedke; dee = deedee;
  }

  //! Assignment constructor.
  NavierStokes2D_cState(const double &dens,
			const double &dvx,
			const double &dvy,
			const double &Etot,
			const double &dkdk,
			const double &domga,
			const double &dkedke,
			const double &deedee) {
    rho = dens; dv.x = dvx; dv.y = dvy; E = Etot; dk = dkdk; domega = domga; dke = dkedke; dee = deedee;
  }

  //! Destructor.
  ~NavierStokes2D_cState(void) { }
  //@}

  //! Return the number of variables.
  static int NumVar(void) { return NUM_VAR_NAVIERSTOKES2D; }

  //! Number of active solution state variables (i.e. it depends on the flow type)
  static int NumVarActive(void);

  //@{ @name Set static variables.
  void set_static_variables(void);
  void set_static_variables(char *gas_type,
			    const int &FlowType,
			    const double &C_constant,
			    const double &von_karman_constant,
			    const double &yplus_sub,
			    const double &yplus_buffer,
			    const double &yplus_outer,
			    char *propellant_type,
			    const int &i_Compressibility_Effect,
			    const int &i_Transition_Model,
			    const int &i_Variable_Prandtl,
			    const double &sigv,
			    const double &lengthw,
			    const double &C_lamb,
			    const double &Cd_1,
			    const double &Cd_4,
			    const double &Cd_5);
  void set_gas(char *gas_type);
  void set_turbulence(const double &C_constant,
		      const double &von_karman_constant,
		      const double &yplus_sub,
		      const double &yplus_buffer,
		      const double &yplus_outer,
                      const int &i_Compressibility_Effect,
		      const int &i_Transition_Model,
		      const int &i_Variable_Prandtl,
		      const double &C_lamb,
		      const double &Cd_1,
		      const double &Cd_4,
		      const double &Cd_5);
  void set_propellant(char *propellant_type,
		      const double &sigv,
		      const double &lengthw);
  //@}

  //@{ @name Useful operators.
  //! Copy operator.
  void Copy(const NavierStokes2D_cState &U);

  //! Vacuum operator.
  void Vacuum(void);

  //! Standard atmosphere operator.
  void Standard_Atmosphere(void);

  //! One operator. Set the solution to ONE.
  void One(void);

  //! Check for unphysical state properties.
  int Unphysical_Properties(void) const;
  //@}

  //@{ @name Functions required for multigrid.
  //! Copy variables solved by multigrid only.
  void Copy_Multigrid_State_Variables(const NavierStokes2D_cState &Ufine);

  //! Zero variables not-solved by multigrid.
  void Zero_Non_Multigrid_State_Variables(void);
  //@}

  //@{ @name Gas related functions.
  //! Gas flow velocity.
  Vector2D v(void) const;

  //! Gas flow velocity.
  double v(const Vector2D &n) const;

  //! Gas pressure.
  double p(void) const;

  //! Gas temperature.
  double T(void) const;

  //! Gas specific internal energy.
  double e(void) const;

  //! Gas specific enthalpy.
  double h(void) const;

  //! Gas total enthalpy.
  double H(void) const;

  //! Gas sound speed.
  double a(void) const;

  //! Gas sound speed squared.
  double a2(void) const;

  //! Gas Mach number.
  double M(void) const;

  //! Gas specific entropy.
  double s(void) const;

  //! Gas stagnation temperature.
  double To(void) const;

  //! Gas stagnation pressure.
  double po(void) const;

  //! Gas stagnation sound speed.
  double ao(void) const;

  //! Gas stagnation enthalpy.
  double ho(void) const;

  //! Gas dynamic viscosity.
  double mu(void) const;

  //! Gas kinematic viscosity.
  double nu(void) const;

  //! Gas thermal heat conductivity.
  double kappa(void) const;

  //! Prandtl number.
  double Pr(void) const;

  //! Thermal diffusivity.
  double Alpha(void) const;

  //! Gas mean free path.
  double meanfreepath(void) const;
  //@}

  //@{ @name Turbulence (k-omega) related functions.
  //! Return the turbulent kinetic energy.
  double k(void) const;

  //! Return the turbulent specific dissipation.
  double omega(void) const;

  //! Return the internal energy variance rate.
  double ke(void) const;

  //! Return the internal energy variance dissipation rate.
  double ee(void) const;

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

  //! Return the turbulent Reynolds number
  double ReT(void) const;

  //! Return the alpha0 star coefficient(Wilcox transition model)
  double alpha_o_star(void) const;

  //! Return the beta star coefficient(Wilcox transition model)
  double beta_star(void) const;

  //! Return the alpha coefficient(Wilcox transition model)
  double alpha_Wilcox(void) const;

  //! Return the alpha star coefficient(Wilcox transition model)
  double alpha_star_Wilcox(void) const;

  //! Turbulent eddy dynamic viscosity.
  double muT(void) const;

  //! Turbulent Prandtl number.
  double PrT(const NavierStokes2D_pState &W,const double &ywall,const double &yplus) const;

  //! Turbulent thermal diffusivity
  double alphaT(const NavierStokes2D_pState &W,const double &ywall,const double &yplus) const;

  //! Turbulent eddy kinematic viscosity.
  double nuT(void) const;

  //! Turbulent eddy thermal heat conductivity.
  double kappaT(const NavierStokes2D_pState &W, const double &ywall,const double &yplus) const;

  //! Turbulence modified sound speed.
  double c(void) const;

  //! Turbulence modified sound speed squared.
  double c2(void) const;

  //! Turbulence modified pressure.
  double pmodified(void) const;

  //! Kolmogorov velocity scale.  
  double U_e(void) const;

  //! Near Wall damping function (variable Prandtl number model).
  double f_lambda(const double &ywall)const;

  //! Near Wall damping function (variable Prandtl number model).
  double f_mu(const double &ywall)const;

  //! Near Wall damping function (variable Prandtl number model).
  double xi_et(const double &ywall, const NavierStokes2D_pState &dWdx, const NavierStokes2D_pState &dWdy)const;

  //! Beta_star value for the k equation.
  double beta_k(const NavierStokes2D_pState &dWdx, const NavierStokes2D_pState &dWdy) const;

  //! Beta value for the omega equation.
  double beta_omega(const NavierStokes2D_pState &dWdx, const NavierStokes2D_pState &dWdy) const;
  
  //! F_beta_star function for the k equation.
  double f_beta_k(const NavierStokes2D_pState &dWdx, const NavierStokes2D_pState &dWdy) const;

  //! F_beta function for the omega equation.
  double f_beta_omega(const NavierStokes2D_pState &dWdx, const NavierStokes2D_pState &dWdy) const;

  //! Sarkar dilatation dissipation correction function.
  double F_Sarkar(void) const;

  //! Zeman dilatation dissipation correction function.
  double F_Zeman(void) const;

  //! Wilcox dilatation dissipation correction function.
  double F_Wilcox(void) const;
  
  //! Chi_k in the f_beta_star function in k equation.
  double chi_k(const NavierStokes2D_pState &dWdx, const NavierStokes2D_pState &dWdy) const;

  //! Chi_Omega in the f_beta function in the omega equation.
  double chi_omega(const NavierStokes2D_pState &dWdx, const NavierStokes2D_pState &dWdy) const;
  //@}

  //@{ @name Solid propellant related functions.
  //! Burning rate.
  double burningrate(void) const;
  //@}

  //@{ @name Primitive solution state.
  NavierStokes2D_pState W(void) const;
  NavierStokes2D_pState W(const NavierStokes2D_cState &U) const;
  friend NavierStokes2D_pState W(const NavierStokes2D_cState &U);
  //@}

  //! Jacobian of the conserved solution variables with respect to the
  //! primitive solution variables.
  void dUdW(DenseMatrix &dUdW) const;

  //! Jacobian of the primitive solution variables with respect to the
  //! conserved solution variables.
  void dWdU(DenseMatrix &dWdU) const;

  //@{ @name Inviscid solution flux (x-direction) and Jacobian.
  NavierStokes2D_cState F(void) const;
  NavierStokes2D_cState F(const Vector2D &V) const;
  void dFdU(DenseMatrix &dFdU) const;
  //@}

  //@{ @name Viscous solution flux and Jacobians.
  NavierStokes2D_cState Gx(const NavierStokes2D_pState &dWdx,const NavierStokes2D_pState &W,
			   const double &ywall,const double &yplus) const;
  NavierStokes2D_cState Gy(const NavierStokes2D_pState &dWdy,const NavierStokes2D_pState &W,
			   const double &ywall,const double &yplus) const;
  //@}

  //! Compute viscous stress tensor and heat flux vector.
  void ComputeViscousTerms(const NavierStokes2D_pState &dWdx,
			   const NavierStokes2D_pState &dWdy,
			   const NavierStokes2D_pState &W,
			   const Vector2D &X,
			   const int &Axisymmetric,
			   const int &adiabatic_flag,
			   const double &ywall,
			   const double &yplus);

  //@{ @name Eigenstructure.
  //! Eigenvalue(s) (x-direction).
  NavierStokes2D_pState lambda_x(void) const;

  //! Eigenvalue(s) (x-direction).
  NavierStokes2D_pState lambda_x(const Vector2D &V) const;

  //! Primitive right eigenvector (x-direction).
  NavierStokes2D_pState rp_x(int index) const;

  //! Conserved right eigenvector (x-direction).
  NavierStokes2D_cState rc_x(int index) const;

  //! Primitive left eigenvector (x-direction).
  NavierStokes2D_pState lp_x(int index) const;
  //@}

  //@{ @name Index operator.
  double &operator[](int index) {
    //assert(index >= 1 && index <= NUM_VAR_NAVIERSTOKES2D);
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
    case 7 :
      return dke;
    case 8 :
      return dee;
    };
    // Default return, this is never reached.
    return rho;
  }

  const double &operator[](int index) const {
    //assert(index >= 1 && index <= NUM_VAR_NAVIERSTOKES2D);
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
    case 7 :
      return dke;
    case 8 :
      return dee;
    };
    // Default return, this is never reached.
    return rho;
  }
  //@}

  //@{ @name Binary arithmetic operators.
  NavierStokes2D_cState operator +(const NavierStokes2D_cState &U) const;
  NavierStokes2D_cState operator -(const NavierStokes2D_cState &U) const;
  double operator *(const NavierStokes2D_cState &U) const;
  NavierStokes2D_cState operator *(const double &a) const;
  friend NavierStokes2D_cState operator *(const double &a, const NavierStokes2D_cState &U){ return U*a; }
  NavierStokes2D_cState operator /(const double &a) const;
  NavierStokes2D_cState operator ^(const NavierStokes2D_cState &U) const;
  //@}

  //@{ @name Assignment operator.
  NavierStokes2D_cState& operator =(const NavierStokes2D_cState &U);
  //@}

  //@{ @name Unary arithmetic operators.
  friend NavierStokes2D_cState operator +(const NavierStokes2D_cState &U);
  friend NavierStokes2D_cState operator -(const NavierStokes2D_cState &U);
  friend NavierStokes2D_cState fabs(const NavierStokes2D_cState &U);
  friend NavierStokes2D_cState sqr(const NavierStokes2D_cState &U);
  //@}

  //@{ @name Shortcut arithmetic operators.
  NavierStokes2D_cState &operator +=(const NavierStokes2D_cState &U);
  NavierStokes2D_cState &operator -=(const NavierStokes2D_cState &U);
  NavierStokes2D_cState &operator *=(const double &a);
  NavierStokes2D_cState &operator /=(const double &a);
  //@}

  //@{ @name Relational operators.
  friend int operator ==(const NavierStokes2D_cState &U1, const NavierStokes2D_cState &U2);
  friend int operator !=(const NavierStokes2D_cState &U1, const NavierStokes2D_cState &U2);
  friend bool operator >=(const NavierStokes2D_cState& U1, const NavierStokes2D_cState& U2);
  friend bool operator <=(const NavierStokes2D_cState& U1, const NavierStokes2D_cState& U2);
  friend bool operator <(const NavierStokes2D_cState &U1, const NavierStokes2D_cState &U2);
  friend bool operator >(const NavierStokes2D_cState &U1, const NavierStokes2D_cState &U2);  
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const NavierStokes2D_cState &U);
  friend istream &operator >> (istream &in_file, NavierStokes2D_cState &U);
  //@}

};


/*****************************************************************************************************
 *****************************************************************************************************
 *                                                                                                   *
 *                 NavierStokes2D_pState CLASS MEMBER FUNCTIONS IMPLEMENTATION                       *
 *                                                                                                   *
 *****************************************************************************************************
 ****************************************************************************************************/

/**********************************************************************
 * NavierStokes2D_pState::NumVarActive --                             *
 * Returns number of state variables based on the flow type.          *
 **********************************************************************/
inline int NavierStokes2D_pState::NumVarActive(void){
  switch(flow_type){
  case FLOWTYPE_INVISCID:
    return 4;			// Four primitive variables (rho, v.x, v.y, p)
    break;

  case FLOWTYPE_LAMINAR:
    return 4;                   // Four primitive variables (rho, v.x, v.y, p)
    break;

  case FLOWTYPE_TURBULENT_RANS_K_OMEGA:
    if (Variable_Prandtl){
      return NUM_VAR_NAVIERSTOKES2D; // Eight primitive variables (rho, v.x, v.y, p, k, omega, ke, ee )
    } else {
      return 6;			// Six primitive variables (rho, v.x, v.y, p, k, omega)
    }
    break;

  default:
    return NUM_VAR_NAVIERSTOKES2D;
  }
}

/********************************************
 * NavierStokes2D_pState Value Constructor. *
 *******************************************/
inline NavierStokes2D_pState::NavierStokes2D_pState(const double &Val):
  rho(Val), v(Val), p(Val), k(Val), omega(Val), ke(Val), ee(Val) {
  // Nothing
}

/**********************************************************************
 * NavierStokes2D_pState::Copy -- Copy operator.                      *
 **********************************************************************/
inline void NavierStokes2D_pState::Copy(const NavierStokes2D_pState &W) {
  rho = W.rho; v.x = W.v.x; v.y = W.v.y; p = W.p;
  k = W.k; omega = W.omega; ke = W.ke; ee = W.ee;
}

/**********************************************************************
 * NavierStokes2D_pState::Vacuum -- Vacuum operator.                  *
 **********************************************************************/
inline void NavierStokes2D_pState::Vacuum(void) {
  rho = ZERO; v.x = ZERO; v.y = ZERO; p = ZERO;
  k = ZERO; omega = ZERO; ke = ZERO; ee = ZERO;
  tau.zero(); q.zero();
}

/**********************************************************************
 * NavierStokes2D_pState::Standard_Atmosphere -- Standard atmosphere  *
 *                                               operator.            *
 **********************************************************************/
inline void NavierStokes2D_pState::Standard_Atmosphere(void) {
  rho = DENSITY_STDATM; v.x = ZERO; v.y = ZERO; p = PRESSURE_STDATM;
  k = ZERO; omega = ZERO; ke = ZERO; ee = ZERO;
  tau.zero(); q.zero();
}

/****************************************************************
 * NavierStokes2D_pState::One -- One operator.                  *
 ***************************************************************/
inline void NavierStokes2D_pState::One(void) {
  rho = ONE; v.x = ONE; v.y = ONE; p = ONE;
  k = ONE; omega = ONE; ke = ONE; ee = ONE;
}

/**********************************************************************
 * NavierStokes2D_pState::Unphysical_Properties -- Check for          *
 *                                                 unphysical state   *
 *                                                 properties.        *
 **********************************************************************/
inline int NavierStokes2D_pState::Unphysical_Properties(void) const {
  if (rho <= ZERO || p <= ZERO || E() <= ZERO) return 1;
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) if (k < ZERO || omega < ZERO) return 1;
  if (Variable_Prandtl == ON) if (ke < ZERO || ee < ZERO) return 1;
  return 0;
}

/**********************************************************************
 * NavierStokes2D_pState::set_static_variables -- Set all static      *
 *                                                variables.          *
 **********************************************************************/
inline void NavierStokes2D_pState::set_static_variables(void) {
  // Set gas constants.
  set_gas("AIR");
  // Set the flow type.
  flow_type = FLOWTYPE_LAMINAR;
  // Set turbulence constants.
  set_turbulence(ZERO,ZERO,ZERO,ZERO,ZERO,0,0,0,ZERO,ZERO,ZERO,ZERO);
  // Set propellant type.
  set_propellant("AP_HTPB",ZERO,ZERO);
}

inline void NavierStokes2D_pState::set_static_variables(char *gas_type,
							const int &FlowType,
							const double &C_constant,
							const double &von_karman_constant,
							const double &yplus_sub,
							const double &yplus_buffer,
							const double &yplus_outer,
							char *propellant_type,
							const int &i_Compressibility_Effect,
                                                        const int &i_Transition_Model,
							const int &i_Variable_Prandtl,
							const double &sigv,
							const double &lengthw,
							const double &C_lamb,
							const double &Cd_1,
							const double &Cd_4,
							const double &Cd_5) {
  // Set gas constants.
  set_gas(gas_type);
  // Set the flow type.
  flow_type = FlowType;
  // Set turbulence constants.
  set_turbulence(C_constant,von_karman_constant,
		 yplus_sub,yplus_buffer,yplus_outer,
		 i_Compressibility_Effect,
		 i_Transition_Model,
		 i_Variable_Prandtl,
		 C_lamb,Cd_1,
		 Cd_4,Cd_5);
  // Set propellant type.
  set_propellant(propellant_type,
		 sigv,
		 lengthw);
}

/**********************************************************************
 * NavierStokes2D_pState::set_gas -- Set gas static variables.        *
 **********************************************************************/
inline void NavierStokes2D_pState::set_gas(char *gas_type) {
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
  } else if (strcmp(gas_type,"AIR-viscous") == 0) {
    g  = GAMMA_AIR;
    R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    v1 = 0.0000056; v2 = AIR_c2; v3 = AIR_c3; v4 = AIR_c4; v5 = AIR_c5;
  } else if (strcmp(gas_type,"AIR-constant") == 0) {
    // This AIR has the viscosity and thermal conductivity constant and equal to the ones at standard temperature (i.e. 288.15 K)
    g  = GAMMA_AIR;
    R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    v1 = 0.0; v2 = 0.0; v3 = 0.0; v4 = 0.0; v5 = 1.83255359613465e-5;
  } else {
    g  = GAMMA_AIR;
    R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    v1 = AIR_c1; v2 = AIR_c2; v3 = AIR_c3; AIR_c4;
  }
  gm1 = g-ONE;
  gm1i = ONE/gm1;
  cp = g*R*gm1i;
  cv = R*gm1i;
}

/**********************************************************************
 * NavierStokes2D_pState::set_turbulence -- Set the turbulence static *
 *                                          variables.                *
 **********************************************************************/
inline void NavierStokes2D_pState::set_turbulence(const double &C_constant,
						  const double &von_karman_constant,
						  const double &yplus_sub,
						  const double &yplus_buffer,
						  const double &yplus_outer,
						  const int &i_Compressibility_Effect,
						  const int &i_Transition_Model,
						  const int &i_Variable_Prandtl,
						  const double &C_lamb,
						  const double &Cd_1,
						  const double &Cd_4,
						  const double &Cd_5) {

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
  Compressibility_Effect = i_Compressibility_Effect;
  Transition_Model = i_Transition_Model;
  Variable_Prandtl = i_Variable_Prandtl;
  C_lambda = C_lamb;
  Cd1 = Cd_1;
  Cd4 = Cd_4;
  Cd5 = Cd_5;
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
 * NavierStokes2D_pState::set_propellant -- Set propellant static     *
 *                                          variables.                *
 **********************************************************************/
inline void NavierStokes2D_pState::set_propellant(char *propellant_type,
						  const double &sigv,
						  const double &lengthw) {
  if (strcmp(propellant_type,"AP_HTPB") == 0) {
    rhos   = RHOS_APHTPB;
    n      = N_APHTPB;
    beta   = BETA_APHTPB;
    Tf     = TF_APHTPB;
    Ts     = TS_APHTPB;
  } else if (strcmp(propellant_type,"QUICK_AP_HTPB") == 0) {
    rhos   = RHOS_QUICK;
    n      = N_QUICK;
    beta   = BETA_QUICK;
    Tf     = TF_QUICK;
    Ts     = TS_QUICK;
  } else if (strcmp(propellant_type,"PLAID_AP_HTPB") == 0) {
    rhos   = RHOS_PLAID;
    n      = N_PLAID;
    beta   = BETA_PLAID;
    Tf     = TF_PLAID;
    Ts     = TS_PLAID;
  } else {
    rhos   = RHOS_APHTPB;
    n      = N_APHTPB;
    beta   = BETA_APHTPB;
    Tf     = TF_APHTPB;
    Ts     = TS_APHTPB;
  }
  sigmav = sigv;
  lw = lengthw;
}

/**********************************************************************
 * NavierStokes2D_pState::T -- Gas temperature.                       *
 **********************************************************************/
inline double NavierStokes2D_pState::T(void) const {
  //assert(rho > ZERO);
  return p/(rho*R);
}

/**********************************************************************
 * NavierStokes2D_pState::e -- Gas specific internal energy.          *
 **********************************************************************/
inline double NavierStokes2D_pState::e(void) const {
  //assert(rho > ZERO);
  return p/(gm1*rho);
}

/**********************************************************************
 * NavierStokes2D_pState::E -- Gas total energy.                      *
 **********************************************************************/
inline double NavierStokes2D_pState::E(void) const {
  return p*gm1i + HALF*rho*v.sqr() + dk();
}

/**********************************************************************
 * NavierStokes2D_pState::h -- Gas specific enthalpy.                 *
 **********************************************************************/
inline double NavierStokes2D_pState::h(void) const {
  //assert(rho > ZERO);
  return g*gm1i*p/rho + HALF*v.sqr() + k;
}

/**********************************************************************
 * NavierStokes2D_pState::H -- Gas total enthalpy.                    *
 **********************************************************************/
inline double NavierStokes2D_pState::H(void) const {
  return g*gm1i*p + HALF*rho*v.sqr() + dk();
}

/**********************************************************************
 * NavierStokes2D_pState::a -- Gas sound speed.                       *
 **********************************************************************/
inline double NavierStokes2D_pState::a(void) const {
  //assert(rho > ZERO);
  return sqrt(g*p/rho);
}

/**********************************************************************
 * NavierStokes2D_pState::a2 -- Gas sound speed squared.              *
 **********************************************************************/
inline double NavierStokes2D_pState::a2(void) const {
  //assert(rho > ZERO);
  return g*p/rho;
}

/**********************************************************************
 * NavierStokes2D_pState::M -- Gas Mach number.                       *
 **********************************************************************/
inline double NavierStokes2D_pState::M(void) const {
  //assert(rho > ZERO && p > ZERO);
  return abs(v)/a();
}

/**********************************************************************
 * NavierStokes2D_pState::s -- Gas specific entropy.                  *
 **********************************************************************/
inline double NavierStokes2D_pState::s(void) const {
  //assert(rho > ZERO && p > ZERO);
  return R*gm1i*log(p/pow(rho,g));
}

/**********************************************************************
 * NavierStokes2D_pState::dv -- Gas momentum.                         *
 **********************************************************************/
inline Vector2D NavierStokes2D_pState::dv(void) const {
  return rho*v;
}

/**********************************************************************
 * NavierStokes2D_pState::dv -- Gas momentum.                         *
 **********************************************************************/
inline double NavierStokes2D_pState::dv(const Vector2D &n) const {
  return rho*(v*n);
}

/**********************************************************************
 * NavierStokes2D_pState::To -- Gas stagnation temperature.           *
 **********************************************************************/
inline double NavierStokes2D_pState::To(void) const {
  //assert(rho > ZERO && p > ZERO);
  return (p/(rho*R))*(ONE+HALF*gm1*v.sqr()/(g*p/rho));
}

/**********************************************************************
 * NavierStokes2D_pState::po -- Gas stagnation pressure.              *
 **********************************************************************/
inline double NavierStokes2D_pState::po(void) const {
  //assert(rho > ZERO && p > ZERO);
  return p*pow(ONE+HALF*gm1*v.sqr()/(g*p/rho),g*gm1i);
}

/**********************************************************************
 * NavierStokes2D_pState::ao -- Gas stagnation sound speed.           *
 **********************************************************************/
inline double NavierStokes2D_pState::ao(void) const {
  //assert(rho > ZERO && p > ZERO);
  return sqrt((g*p/rho)*(ONE+HALF*gm1*v.sqr()/(g*p/rho)));
}

/**********************************************************************
 * NavierStokes2D_pState::ho -- Gas stagnation enthalpy.              *
 **********************************************************************/
inline double NavierStokes2D_pState::ho(void) const {
  //assert(rho > ZERO && p > ZERO);
  return (g*p/(gm1*rho) + HALF*v.sqr())*(ONE+HALF*gm1*v.sqr()/(g*p/rho));
}

/**********************************************************************
 * NavierStokes2D_pState::mu -- Gas dynamic viscosity.                *
 **********************************************************************/
inline double NavierStokes2D_pState::mu(void) const {
  return mu_gottlieb(v1,v2,v3,v4,v5,T());
}

/**********************************************************************
 * NavierStokes2D_pState::nu -- Gas kinematic viscosity.              *
 **********************************************************************/
inline double NavierStokes2D_pState::nu(void) const {
  return mu()/rho;
}

/**********************************************************************
 * NavierStokes2D_pState::kappa -- Gas thermal heat conductivity.     *
 **********************************************************************/
inline double NavierStokes2D_pState::kappa(void) const {
  return kappa_gottlieb(v1,v2,v3,v4,v5,T(),g,cp);
}

/**********************************************************************
 * NavierStokes2D_pState::Pr -- Prandtl number.                       *
 **********************************************************************/
inline double NavierStokes2D_pState::Pr(void) const {
  //assert(kappa() > ZERO);
  return cp*mu()/kappa();
}

/**********************************************************************
 * NavierStokes2D_pState::Alpha -- Thermal diffusivity.               *
 **********************************************************************/
inline double NavierStokes2D_pState::Alpha(void) const {
  return nu()/Pr();
}

/**********************************************************************
 * NavierStokes2D_pState::meanfreepath -- Gas mean free path.         *
 **********************************************************************/
inline double NavierStokes2D_pState::meanfreepath(void) const {
  //assert(rho > ZERO && T() > ZERO);
  return 16.0*mu()/(5.0*rho*sqrt(2.0*PI*R*T()));
}

/**********************************************************************
 * NavierStokes2D_pState::dk -- Total turbulent kinetic energy.       *
 **********************************************************************/
inline double NavierStokes2D_pState::dk(void) const {
  return rho*k;
}

/**********************************************************************
 * NavierStokes2D_pState::domega -- Total turbulent specific          *
 *                                  dissipation.                      *
 **********************************************************************/
inline double NavierStokes2D_pState::domega(void) const {
  return rho*omega;
}

/**********************************************************************
 * NavierStokes2D_pState::dke -- Total internal energy variance.      *
 **********************************************************************/
inline double NavierStokes2D_pState::dke(void) const {
  return rho*ke;
}

/**********************************************************************
 * NavierStokes2D_pState::dee -- Total internal energy variance       *
 *                               dissipation rate.                    *
 **********************************************************************/
inline double NavierStokes2D_pState::dee(void) const {
  return rho*ee;
}

/**********************************************************************
 * NavierStokes2D_pState::epsilon -- Gas specific turbulent eddy      *
 *                                   dissipation.                     *
 **********************************************************************/
inline double NavierStokes2D_pState::epsilon(void) const {
  return beta_k_o*k*omega;
}

/**********************************************************************
 * NavierStokes2D_pState::depsilon -- Total turbulent eddy            *
 *                                    dissipation.                    *
 **********************************************************************/
inline double NavierStokes2D_pState::depsilon(void) const {
  return rho*epsilon();
}

/**********************************************************************
 * NavierStokes2D_pState::ell -- Turbulent length scale.              *
 **********************************************************************/
inline double NavierStokes2D_pState::ell(void) const {
  return sqrt(k)/max(omega,NANO);
}

/**********************************************************************
 * NavierStokes2D_pState::Mt -- Turbulent Mach number.                *
 **********************************************************************/
inline double NavierStokes2D_pState::Mt(void) const {
  return sqrt(TWO*k/a2());
}

/**********************************************************************
 * NavierStokes2D_pState::Mt2 -- Turbulent Mach number squared.       *
 **********************************************************************/
inline double NavierStokes2D_pState::Mt2(void) const {
  return TWO*k/a2();
}

/**********************************************************************
 * NavierStokes2D_pState::ReT -- Gas turbulent Reynolds number.       *
 **********************************************************************/
inline double NavierStokes2D_pState::ReT(void) const {
  return k/max(omega,TOLER)/nu();
}

/**********************************************************************
 * NavierStokes2D_pState::alpha_o_star -- Wilcox's Transition model.  *
 **********************************************************************/
inline double NavierStokes2D_pState::alpha_o_star(void) const {
  return beta_Wilcox / THREE; 
}

/**********************************************************************
 * NavierStokes2D_pState::beta_star --  Wilcox's Transition model.    *
 **********************************************************************/
inline double NavierStokes2D_pState::beta_star(void) const {
  return NINE/100.0*((FIVE/18.0)+pow(ReT()/Rbeta,4.0))/(ONE+pow(ReT()/Rbeta,4.0)) ;
}

/**********************************************************************
 * NavierStokes2D_pState::alpha --  Wilcox's Transition model.        *
 **********************************************************************/
inline double NavierStokes2D_pState::alpha_Wilcox(void) const {
  return FIVE/NINE * (alpha_o + (ReT()/Romega))/(ONE + (ReT()/Romega))/alpha_star_Wilcox();
}

/**********************************************************************
 * NavierStokes2D_pState::alpha_star --  Wilcox's Transition model.   *
 **********************************************************************/
inline double NavierStokes2D_pState::alpha_star_Wilcox(void) const{
  return (alpha_o_star() + (ReT())/(Rk))/(ONE + (ReT())/(Rk));
}

/**********************************************************************
 * NavierStokes2D_pState::muT -- Turbulent eddy dynamic viscosity.    *
 **********************************************************************/
inline double NavierStokes2D_pState::muT(void) const {
  return rho*nuT();
}

/**********************************************************************
 * NavierStokes2D_pState::PrT -- Turbulent Prandtl number.            *
 **********************************************************************/
inline double NavierStokes2D_pState::PrT(const double &ywall,const double &yplus) const {
  if (Variable_Prandtl == OFF){
    double Vpr = 0.9;
    return Vpr;
  }else {    
    /************************************************************************************
     * The original formula given in the Calhoon, Brinckman paper (AIAA-2006-1452)      *
     * uses Cmu and Fmu() for calculating the PrT, where Fmu is a damping function used *
     * by the K-epsilon model. We replace the product Cmu*Fmu() by the equivalent in    * 
     * terms of muT,k and epsilon (eqn 11) in the paper.                                *
     ************************************************************************************/
    double CmuFmu = muT()*epsilon()/rho/sqr(max(k,TOLER*TOLER));
    double Vpr = CmuFmu/C_lambda/f_lambda(ywall)*sqrt(k*ee/max(epsilon(),TOLER*TOLER)/max(ke,TOLER*TOLER));  
    return Vpr;
  } 
}

/**********************************************************************
 * NavierStokes2D_pState::alphaT -- Turbulent thermal diffusivity.    *
 **********************************************************************/
inline double NavierStokes2D_pState::alphaT(const double &ywall,const double &yplus) const {
 return nuT()/max(PrT(ywall,yplus),TOLER);
}

/**********************************************************************
 * NavierStokes2D_pState::U_e -- Kolmogorov velocity scale.           *
 **********************************************************************/
inline double NavierStokes2D_pState::U_e(void) const {
  return sqrt(sqrt(nu()*epsilon()));
}

/**********************************************************************
 * NavierStokes2D_pState::f_lambda -- Wall damping function(Variable  * 
 *                                     Prandtl number model).         *
 **********************************************************************/
inline double NavierStokes2D_pState::f_lambda(const double &ywall) const {
  double Re_t = sqr(k)/nu()/max(epsilon(),TOLER);
  double f_et = exp(-sqr(Re_t/80.0));
  double y_star = U_e()*ywall/nu();
  return f_et*C1_lambda/sqrt(sqrt(max(ReT(),TOLER)/beta_k_o)) + sqr(ONE - exp(-y_star/A_plus));
}

/**********************************************************************
 * NavierStokes2D_pState::xi_et -- Wall damping function(Variable     * 
 *                                     Prandtl number model).         *
 **********************************************************************/
inline double NavierStokes2D_pState::xi_et(const double &ywall,
					   const NavierStokes2D_pState &dWdx,
					   const NavierStokes2D_pState &dWdy) const {
  double deriv = (ONE/(g-1)/pow(rho,TWO))*(rho*dWdx.p-p*dWdx.rho);
  double Re_t = sqr(k)/nu()/max(epsilon(),TOLER);
  double f_et = exp(-sqr(Re_t/80.0));
  double Pe_star = -sqrt(TWO/THREE*k*ke)*deriv;
  double ee_star = ee - Alpha()*ke/max(sqr(ywall),TOLER);
  double e_cap = epsilon() - TWO*nu()/FOUR/max(k,TOLER)*sqr(dWdy.k);
  double ee_cap = ee - Alpha()/FOUR/max(ke,TOLER)*sqr(dWdy.ke);
  return f_et*rho*((Cd4-FOUR)*ee_cap*ee/max(ke,TOLER)+Cd5*e_cap*ee/max(k,TOLER)-sqr(ee_star)/max(ke,TOLER)+(TWO-Cd1-Cd2*Pr())*ee*Pe_star/max(ke,TOLER));
}

/**********************************************************************
 * NavierStokes2D_pState::f_mu -- Wall damping function(Variable      * 
 *                                  Prandtl number model).            *
 **********************************************************************/
inline double NavierStokes2D_pState::f_mu(const double &ywall) const {
  double Rek;
  Rek = rho*sqrt(k)*ywall/mu();
  return (ONE + 4.0/pow(max(ReT(),TOLER)/beta_k_o,THREE/FOUR))*tanh(Rek/125.0);
}

/**********************************************************************
 * NavierStokes2D_pState::nuT -- Turbulent eddy kinematic viscosity.  *
 **********************************************************************/
inline double NavierStokes2D_pState::nuT(void) const {
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA){
    if (Transition_Model == TRANSITION_WILCOX) return alpha_o_star()*k/max(omega,TOLER);
    return k/max(omega,TOLER);
  }  
  return ZERO;
}

/**********************************************************************
 * NavierStokes2D_pState::kappaT -- Turbulent eddy thermal heat       *
 *                                  conductivity.                     *
 **********************************************************************/
inline double NavierStokes2D_pState::kappaT(const double &ywall,const double &yplus) const {
  return muT()*cp/max(PrT(ywall,yplus),TOLER);
}

/**********************************************************************
 * NavierStokes2D_pState::c -- Turbulence modified sound speed.       *
 **********************************************************************/
inline double NavierStokes2D_pState::c(void) const {
  return sqrt(c2());
}

/**********************************************************************
 * NavierStokes2D_pState::c2 -- Turbulence modified sound speed       *
 *                              squared.                              *
 **********************************************************************/
inline double NavierStokes2D_pState::c2(void) const {
  //assert(rho > ZERO);
  return a2() + (2.0/3.0)*g*k;
}

/**********************************************************************
 * NavierStokes2D_pState::pmodified -- Turbulence modified pressure.  *
 **********************************************************************/
inline double NavierStokes2D_pState::pmodified(void) const {
  //assert(rho > ZERO);
  return p + (2.0/3.0)*dk();
}

/**********************************************************************
 * NavierStokes2D_pState::beta_k -- k-omega auxilary relation.        *
 **********************************************************************/
inline double NavierStokes2D_pState::beta_k(const NavierStokes2D_pState &dWdx,
					    const NavierStokes2D_pState &dWdy) const {
  double B, xi;
  if (Transition_Model == TRANSITION_WILCOX){
    B =  beta_star();
  }else B = beta_k_o*f_beta_k(dWdx,dWdy);
  if (Compressibility_Effect == OFF){
    return B;
  }else if (Compressibility_Effect == COMPRESSIBILITY_CORRECTION_SARKAR){
    xi = ONE;
    return B*( ONE + xi*F_Sarkar() );
  }else if (Compressibility_Effect == COMPRESSIBILITY_CORRECTION_ZEMAN){
    xi = THREE/FOUR;
    return B*( ONE + xi*F_Zeman() );
  }else if (Compressibility_Effect == COMPRESSIBILITY_CORRECTION_WILCOX){
    xi = THREE/TWO;
    return B*( ONE + xi*F_Wilcox() );
  }
}

/**********************************************************************
 * NavierStokes2D_pState::beta_omega -- k-omega auxilary relation.    *
 **********************************************************************/
inline double NavierStokes2D_pState::beta_omega(const NavierStokes2D_pState &dWdx,
						const NavierStokes2D_pState &dWdy) const {
  double B1,B2;
  if (Transition_Model == TRANSITION_WILCOX){
    B1 =  beta_Wilcox;
    B2 = beta_star();
  } else {
    B1 = beta_omega_o*f_beta_omega(dWdx,dWdy);
    B2 = beta_k_o*f_beta_k(dWdx,dWdy);
  }
  if (Compressibility_Effect == OFF){
    return B1;
  }else if (Compressibility_Effect == COMPRESSIBILITY_CORRECTION_SARKAR){
    double xi = ONE;
    return B1 - B2*xi*F_Sarkar();
  }else if (Compressibility_Effect == COMPRESSIBILITY_CORRECTION_ZEMAN){
    double xi = THREE/FOUR;
    return B1 - B2*xi*F_Zeman();
  }else if (Compressibility_Effect == COMPRESSIBILITY_CORRECTION_WILCOX){
    double xi = THREE/TWO;
    return B1 - B2*xi*F_Wilcox();
  }
}

/**********************************************************************
 * NavierStokes2D_pState::f_beta_k -- k-omega auxilary relation.      *
 **********************************************************************/
inline double NavierStokes2D_pState::f_beta_k(const NavierStokes2D_pState &dWdx,
					      const NavierStokes2D_pState &dWdy) const {
  double chi = chi_k(dWdx,dWdy);
  if (chi <= ZERO) return ONE;
  return (ONE + 680.0*chi*chi)/(ONE + 400.0*chi*chi);
}

/**********************************************************************
 * NavierStokes2D_pState::f_beta_omega -- k-omega auxilary relation.  *
 **********************************************************************/
inline double NavierStokes2D_pState::f_beta_omega(const NavierStokes2D_pState &dWdx,
						  const NavierStokes2D_pState &dWdy) const {
  double chi = chi_omega(dWdx,dWdy);
  return (ONE + 70.0*chi)/(ONE + 80.0*chi);
}

/***************************************************************************
 * NavierStokes2D_pState::F_Sarkar -- k-omega auxilary relation.           *
 ***************************************************************************/
inline double NavierStokes2D_pState::F_Sarkar(void) const {
  return Mt2() ; 
}

/**************************************************************************
 * NavierStokes2D_pState::F_Zeman -- k-omega auxilary relation.           *
 **************************************************************************/
inline double NavierStokes2D_pState::F_Zeman(void) const {
  double Mt0 = 0.10*sqrt(TWO/(g+ONE)) ;
  double Omega = 0.60; 
  return (ONE-exp(-HALF*(g+1)*sqr((sqrt(Mt2()) - Mt0)/Omega)))*heaviside(Mt()-Mt0);
}

/***************************************************************************
 * NavierStokes2D_pState::F_Wilcox -- k-omega auxilary relation.           *
 ***************************************************************************/
inline double NavierStokes2D_pState::F_Wilcox(void) const {
  double Mt0 = ONE/FOUR; 
  return (Mt2()-sqr(Mt0))*heaviside(Mt()-Mt0);
}

/**********************************************************************
 * NavierStokes2D_pState::chi_k -- k-omega auxilary relation.         *
 **********************************************************************/
inline double NavierStokes2D_pState::chi_k(const NavierStokes2D_pState &dWdx,
					   const NavierStokes2D_pState &dWdy) const {
  return (dWdx.k*dWdx.omega + dWdy.k*dWdy.omega)/max(TOLER,cube(omega));
}

/**********************************************************************
 * NavierStokes2D_pState::chi_omega -- k-omega auxilary relation.     *
 **********************************************************************/
inline double NavierStokes2D_pState::chi_omega(const NavierStokes2D_pState &dWdx,
					       const NavierStokes2D_pState &dWdy) const {
  return 0.25*fabs(sqr(dWdx.v.y - dWdy.v.x)*(dWdx.v.x + dWdy.v.y)/max(TOLER,cube(beta_omega_o*omega)));
}

/**********************************************************************
 * NavierStokes2D_pState::burningrate -- Solid propellent burning rate.*
 **********************************************************************/
inline double NavierStokes2D_pState::burningrate(void) const {
  return -beta*pow(p,n);
}

/**********************************************************************
 * NavierStokes2D_pState -- Binary arithmetic operators.              *
 **********************************************************************/
inline NavierStokes2D_pState NavierStokes2D_pState::operator +(const NavierStokes2D_pState &W) const {
  return NavierStokes2D_pState(rho+W.rho,v.x+W.v.x,v.y+W.v.y,p+W.p,k+W.k,omega+W.omega,ke+W.ke,ee+W.ee);
}

inline NavierStokes2D_pState NavierStokes2D_pState::operator -(const NavierStokes2D_pState &W) const {
  return NavierStokes2D_pState(rho-W.rho,v.x-W.v.x,v.y-W.v.y,p-W.p,k-W.k,omega-W.omega,ke-W.ke,ee-W.ee);
}

// Inner product operator.
inline double NavierStokes2D_pState::operator *(const NavierStokes2D_pState &W) const {
  return rho*W.rho + v.x*W.v.x + v.y*W.v.y + p*W.p + k*W.k + omega*W.omega + ke*W.ke + ee*W.ee;
}

inline NavierStokes2D_pState NavierStokes2D_pState::operator *(const double &a) const {
  return NavierStokes2D_pState(rho*a,v.x*a,v.y*a,p*a,k*a,omega*a,ke*a,ee*a);
}

inline NavierStokes2D_pState NavierStokes2D_pState::operator /(const double &a) const {
  return NavierStokes2D_pState(rho/a,v.x/a,v.y/a,p/a,k/a,omega/a,ke/a,ee/a);
}

inline NavierStokes2D_pState NavierStokes2D_pState::operator /(const NavierStokes2D_pState &W) {
  return NavierStokes2D_pState(rho/W.rho,v.x/W.v.x,v.y/W.v.y,p/W.p,k/W.k,omega/W.omega,ke/W.ke,ee/W.ee);
}

// A useful solution state product operator.
inline NavierStokes2D_pState NavierStokes2D_pState::operator ^(const NavierStokes2D_pState &W) const {
  return NavierStokes2D_pState(rho*W.rho,v.x*W.v.x,v.y*W.v.y,p*W.p,k*W.k,omega*W.omega,ke*W.ke,ee*W.ee);
}

/*!
 * Compute maximum between 2 states. 
 * Return the state of maximum values.
 */
inline NavierStokes2D_pState max(const NavierStokes2D_pState &W1, const NavierStokes2D_pState &W2 ){
  return NavierStokes2D_pState(max(W1.rho,W2.rho),max(W1.v.x,W2.v.x),max(W1.v.y,W2.v.y),max(W1.p,W2.p),
			       max(W1.k,W2.k), max(W1.omega,W2.omega), max(W1.ke,W2.ke), max(W1.ee,W2.ee));
}

/*!
 * Compute minimum between 2 states. 
 * Return the state of minimum values.
 */
inline NavierStokes2D_pState min(const NavierStokes2D_pState &W1, const NavierStokes2D_pState &W2 ){
  return NavierStokes2D_pState(min(W1.rho,W2.rho),min(W1.v.x,W2.v.x),min(W1.v.y,W2.v.y),min(W1.p,W2.p),
			       min(W1.k,W2.k), min(W1.omega,W2.omega), min(W1.ke,W2.ke), min(W1.ee,W2.ee));
}


/**********************************************************************
 * NavierStokes2D_pState -- Assignment operator.                      *
 **********************************************************************/
inline NavierStokes2D_pState& NavierStokes2D_pState::operator =(const NavierStokes2D_pState &W) {
  // handle self-assignment
  if(this == &W) return *this;

  rho = W.rho; v.x = W.v.x; v.y = W.v.y; p = W.p; k = W.k; omega = W.omega; ke = W.ke; ee = W.ee;
  return *this;
}

/**********************************************************************
 * NavierStokes2D_pState -- Unary arithmetic operators.               *
 **********************************************************************/
inline NavierStokes2D_pState operator +(const NavierStokes2D_pState &W) {
  return W;
}

inline NavierStokes2D_pState operator -(const NavierStokes2D_pState &W) {
  return NavierStokes2D_pState(-W.rho,-W.v.x,-W.v.y,-W.p,-W.k,-W.omega,-W.ke,-W.ee);
}

inline NavierStokes2D_pState fabs(const NavierStokes2D_pState &W){
  return NavierStokes2D_pState(fabs(W.rho), fabs(W.v.x), fabs(W.v.y), fabs(W.p),
			       fabs(W.k), fabs(W.omega), fabs(W.ke), fabs(W.ee));
}

inline NavierStokes2D_pState sqr(const NavierStokes2D_pState &W){
  return NavierStokes2D_pState(sqr(W.rho), sqr(W.v.x), sqr(W.v.y), sqr(W.p),
			       sqr(W.k), sqr(W.omega), sqr(W.ke), sqr(W.ee));
}

/**********************************************************************
 * NavierStokes2D_pState -- Shortcut arithmetic operators.            *
 **********************************************************************/
inline NavierStokes2D_pState& NavierStokes2D_pState::operator +=(const NavierStokes2D_pState &W) {
  rho += W.rho; v.x += W.v.x; v.y += W.v.y; p += W.p; k += W.k; omega += W.omega; ke += W.ke; ee += W.ee;
  return *this;
}

inline NavierStokes2D_pState& NavierStokes2D_pState::operator -=(const NavierStokes2D_pState &W) {
  rho -= W.rho; v.x -= W.v.x; v.y -= W.v.y; p -= W.p; k -= W.k; omega -= W.omega; ke -= W.ke; ee -= W.ee;
  return *this;
}

inline NavierStokes2D_pState& NavierStokes2D_pState::operator /=(const NavierStokes2D_pState &W) {
  rho /= W.rho; v.x /= W.v.x; v.y /= W.v.y; p /= W.p; k /= W.k; omega /= W.omega; ke /= W.ke; ee /= W.ee;
  return *this;
}

inline NavierStokes2D_pState& NavierStokes2D_pState::operator *=(const NavierStokes2D_pState &W) {
  rho *= W.rho; v.x *= W.v.x; v.y *= W.v.y; p *= W.p; k *= W.k; omega *= W.omega; ke *= W.ke; ee *= W.ee;
  return *this;
}

inline NavierStokes2D_pState& NavierStokes2D_pState::operator *=(const double &a) {
  rho *= a; v.x *= a; v.y *= a; p *= a; k *= a; omega *= a; ke *= a; ee *= a;
  return *this;
}

inline NavierStokes2D_pState& NavierStokes2D_pState::operator /=(const double &a) {
  rho /= a; v.x /= a; v.y /= a; p /= a; k /= a; omega /= a; ke /= a; ee /= a;
  return *this;
}

/**********************************************************************
 * NavierStokes2D_pState -- Relational operators.                     *
 **********************************************************************/
inline int operator ==(const NavierStokes2D_pState &W1, const NavierStokes2D_pState &W2) {
  return (W1.rho == W2.rho && W1.v == W2.v         && W1.p == W2.p   && 
	  W1.k == W2.k     && W1.omega == W2.omega && W1.ke == W2.ke && W1.ee == W2.ee);
}

inline int operator !=(const NavierStokes2D_pState &W1, const NavierStokes2D_pState &W2) {
  return (W1.rho != W2.rho || W1.v != W2.v         || W1.p != W2.p || 
	  W1.k != W2.k     || W1.omega != W2.omega || W1.ke != W2.ke || W1.ee != W2.ee);
}

inline bool operator <=(const NavierStokes2D_pState &W1, const NavierStokes2D_pState &W2) {
  return (W1.rho<=W2.rho && W1.v.x<=W2.v.x && W1.v.y<=W2.v.y && W1.p<=W2.p &&
	  W1.k<=W2.k && W1.omega<=W2.omega && W1.ke<=W2.ke && W1.ee<=W2.ee );
}

inline bool operator >=(const NavierStokes2D_pState &W1, const NavierStokes2D_pState &W2) {
  return (W1.rho>=W2.rho && W1.v.x>=W2.v.x && W1.v.y>=W2.v.y && W1.p>=W2.p &&
	  W1.k>=W2.k && W1.omega>=W2.omega && W1.ke>=W2.ke && W1.ee>=W2.ee );
}

inline bool operator <(const NavierStokes2D_pState &W1, const NavierStokes2D_pState &W2) {
  return (W1.rho<W2.rho && W1.v.x<W2.v.x && W1.v.y<W2.v.y && W1.p<W2.p &&
	  W1.k<W2.k && W1.omega<W2.omega && W1.ke<W2.ke && W1.ee<W2.ee );
}

inline bool operator >(const NavierStokes2D_pState &W1, const NavierStokes2D_pState &W2) {
  return (W1.rho>W2.rho && W1.v.x>W2.v.x && W1.v.y>W2.v.y && W1.p>W2.p &&
	  W1.k>W2.k && W1.omega>W2.omega && W1.ke>W2.ke && W1.ee>W2.ee );
}

/**********************************************************************
 * NavierStokes2D_pState -- Input-output operators.                   *
 **********************************************************************/
inline ostream &operator << (ostream &out_file, const NavierStokes2D_pState &W) {
  out_file.setf(ios::scientific);
  out_file << " " << W.rho << " " << W.v.x << " " << W.v.y << " " << W.p 
	   << " " << W.k << " " << W.omega<< " " << W.ke << " " << W.ee;
  out_file.unsetf(ios::scientific);
  return out_file;
}

inline istream &operator >> (istream &in_file, NavierStokes2D_pState &W) {
  in_file.setf(ios::skipws);
  in_file >> W.rho >> W.v.x >> W.v.y >> W.p >> W.k >> W.omega >> W.ke >> W.ee;
  in_file.unsetf(ios::skipws);
  return in_file;
}

inline void NavierStokes2D_pState::output_labels(ostream &out_file) {
  out_file << "\"rho\" \\ \n"
	   << "\"u\" \\ \n"
	   << "\"v\" \\ \n"
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
    if (Variable_Prandtl == ON) {
      out_file << "\"ke\" \\ \n"
	       << "\"ee\" \\ \n"
	       << "\"PrT\" \\ \n";
    }
  }
}

inline void NavierStokes2D_pState::output_data(ostream &out_file,
					       const double &ywall,
					       const double &yplus) {
  out_file << " " << rho << " " << v.x << " " << v.y << " " << p
	   << " " << T() << " " << M() << " " << H() << " " << s();
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    out_file << " " << k << " " << omega << " " << epsilon()
	     << " " << ell() << " " << pmodified();
    if (Variable_Prandtl == ON) {
      out_file << " " << ke << " " << ee<< " "<<PrT(ywall,yplus);
    }
  }
}


/*****************************************************************************************************
 *****************************************************************************************************
 *                                                                                                   *
 *                 NavierStokes2D_cState CLASS MEMBER FUNCTIONS IMPLEMENTATION                       *
 *                                                                                                   *
 *****************************************************************************************************
 ****************************************************************************************************/

/**********************************************************************
 * NavierStokes2D_cState::NumVarActive --                             *
 * Returns number of state variables based on the flow type.          *
 **********************************************************************/
inline int NavierStokes2D_cState::NumVarActive(void){
  switch(flow_type){
  case FLOWTYPE_INVISCID:
    return 4;			// Four conserved variables (rho, dv.x, dv.y, E)
    break;

  case FLOWTYPE_LAMINAR:
    return 4;                   // Four conserved variables (rho, dv.x, dv.y, E)
    break;

  case FLOWTYPE_TURBULENT_RANS_K_OMEGA:
    if (Variable_Prandtl){
      return NUM_VAR_NAVIERSTOKES2D;  // Eight conserved quantities (rho, dv.x, dv.y, E, dk, domega, dke, dee)
    } else {
      return 6;			// Six conserved quantities (rho, dv.x, dv.y, E, dk, domega)
    }
    break;

  default:
    return NUM_VAR_NAVIERSTOKES2D;
  }
}

/********************************************
 * NavierStokes2D_cState Value Constructor. *
 *******************************************/
inline NavierStokes2D_cState::NavierStokes2D_cState(const double &Val):
  rho(Val), dv(Val), E(Val), dk(Val), domega(Val), dke(Val), dee(Val) {
  // Nothing
}

/**********************************************************************
 * NavierStokes2D_cState::Copy -- Copy operator.                      *
 **********************************************************************/
inline void NavierStokes2D_cState::Copy(const NavierStokes2D_cState &U) {
  rho = U.rho; dv.x = U.dv.x; dv.y = U.dv.y; E = U.E; 
  dk = U.dk; domega = U.domega; dke = U.dke; dee = U.dee;
}

/**********************************************************************
 * NavierStokes2D_cState::Vacuum -- Vacuum state.                     *
 **********************************************************************/
inline void NavierStokes2D_cState::Vacuum(void) {
  rho = ZERO; dv.x = ZERO; dv.y = ZERO; E = ZERO;
  dk = ZERO; domega = ZERO; dke = ZERO; dee = ZERO;
  tau.zero(); q.zero();
}

/**********************************************************************
 * NavierStokes2D_cState::Standard_Atmosphere -- Standard atmosphere  *
 *                                               state.               *
 **********************************************************************/
inline void NavierStokes2D_cState::Standard_Atmosphere(void) {
  rho = DENSITY_STDATM; dv.x = ZERO; dv.y = ZERO;
  E = PRESSURE_STDATM/(GAMMA_AIR-ONE); tau.zero(); q.zero();
  dk = ZERO; domega = ZERO; dke = ZERO; dee = ZERO;
}

/****************************************************************
 * NavierStokes2D_cState::One -- One state.                     *
 ***************************************************************/
inline void NavierStokes2D_cState::One(void) {
  rho = ONE; dv.x = ONE; dv.y = ONE; E = ONE;
  dk = ONE; domega = ONE; dke = ONE; dee = ONE;
}

/**********************************************************************
 * NavierStokes2D_cState::Unphysical_Properties -- Check for          *
 *                                                 unphysical state   *
 *                                                 properties.        *
 **********************************************************************/
inline int NavierStokes2D_cState::Unphysical_Properties(void) const {
  if (rho <= ZERO || E <= ZERO || e() <= ZERO) return 1;
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) if (dk < ZERO || domega < ZERO) return 1;
  if (Variable_Prandtl == ON) if (dke < ZERO || dee < ZERO) return 1;
  return 0;
}

/**********************************************************************
 * NavierStokes2D_cState::Copy_Multigrid_State_Variables --           *
 *                           Copy variables solved by multigrid only. *
 **********************************************************************/
inline void NavierStokes2D_cState::Copy_Multigrid_State_Variables(const NavierStokes2D_cState &Ufine) {
  rho = Ufine.rho; dv.x = Ufine.dv.x; dv.y = Ufine.dv.y; E = Ufine.E;
}

/**********************************************************************
 * NavierStokes2D_cState::Zero_Non_Multigrid_State_Variables --       *
 *                            Zero variables not-solved by multigrid. *
 **********************************************************************/
inline void NavierStokes2D_cState::Zero_Non_Multigrid_State_Variables(void) {
  dk = ZERO; domega = ZERO; dke = ZERO; dee = ZERO;
}

/**********************************************************************
 * NavierStokes2D_cState::set_static_variables -- Set all static      *
 *                                                variables.          *
 **********************************************************************/
inline void NavierStokes2D_cState::set_static_variables(void) {
  // Set gas constants.
  set_gas("AIR");
  // Set the flow type.
  flow_type = FLOWTYPE_LAMINAR;
  // Set turbulence constants.
  set_turbulence(ZERO,ZERO,ZERO,ZERO,ZERO,0,0,0,ZERO,ZERO,ZERO,ZERO);
  // Set propellant type.
  set_propellant("AP_HTPB",ZERO,ZERO);
}

inline void NavierStokes2D_cState::set_static_variables(char *gas_type,
							const int &FlowType,
							const double &C_constant,
							const double &von_karman_constant,
							const double &yplus_sub,
							const double &yplus_buffer,
							const double &yplus_outer,
							char *propellant_type,
                                                        const int &i_Compressibility_Effect,
                                                        const int &i_Transition_Model,  
                                                        const int &i_Variable_Prandtl,
							const double &sigv,
							const double &lengthw,
							const double &C_lamb,
							const double &Cd_1,
							const double &Cd_4,
							const double &Cd_5) {
  // Set gas constants.
  set_gas(gas_type);
  // Set the flow type.
  flow_type = FlowType;
  // Set turbulence constants.
  set_turbulence(C_constant,von_karman_constant,
		 yplus_sub,yplus_buffer,yplus_outer,
		 i_Compressibility_Effect,
		 i_Transition_Model,
		 i_Variable_Prandtl,
		 C_lamb,Cd_1,
		 Cd_4,Cd_5);
  // Set propellant type.
  set_propellant(propellant_type,
		 sigv,
		 lengthw);
}

/**********************************************************************
 * NavierStokes2D_cState::set_gas -- Set gas static variables.         *
 **********************************************************************/
inline void NavierStokes2D_cState::set_gas(char *gas_type) {
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
  } else if (strcmp(gas_type,"AIR-viscous") == 0) {
    g  = GAMMA_AIR;
    R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    v1 = 0.0000056; v2 = AIR_c2; v3 = AIR_c3; v4 = AIR_c4; v5 = AIR_c5;
  } else {
    g  = GAMMA_AIR;
    R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    v1 = AIR_c1; v2 = AIR_c2; v3 = AIR_c3; v4 = AIR_c4; v5 = AIR_c5;
  }
  gm1 = g-ONE;
  gm1i = ONE/gm1;
  cp = g*R*gm1i;
  cv = R*gm1i;
}

/**********************************************************************
 * NavierStokes2D_cState::set_turbulence -- Set the turbulence static *
 *                                          variables.                *
 **********************************************************************/
inline void NavierStokes2D_cState::set_turbulence(const double &C_constant,
						  const double &von_karman_constant,
						  const double &yplus_sub,
						  const double &yplus_buffer,
						  const double &yplus_outer,
						  const int &i_Compressibility_Effect,
						  const int &i_Transition_Model,
						  const int &i_Variable_Prandtl,
						  const double &C_lamb,
						  const double &Cd_1,
						  const double &Cd_4,
						  const double &Cd_5) {
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
  Compressibility_Effect = i_Compressibility_Effect;  
  Transition_Model = i_Transition_Model;
  Variable_Prandtl = i_Variable_Prandtl;
  C_lambda = C_lamb;
  Cd1 = Cd_1;
  Cd4 = Cd_4;
  Cd5 = Cd_5;
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
 * NavierStokes2D_cState::set_propellant -- Set propellant static     *
 *                                          variables.                *
 **********************************************************************/
inline void NavierStokes2D_cState::set_propellant(char *propellant_type,
						  const double &sigv,
						  const double &lengthw) {
  if (strcmp(propellant_type,"AP_HTPB") == 0) {
    rhos   = RHOS_APHTPB;
    n      = N_APHTPB;
    beta   = BETA_APHTPB;
    Tf     = TF_APHTPB;
    Ts     = TS_APHTPB;
  } else if (strcmp(propellant_type,"QUICK_AP_HTPB") == 0) {
    rhos   = RHOS_QUICK;
    n      = N_QUICK;
    beta   = BETA_QUICK;
    Tf     = TF_QUICK;
    Ts     = TS_QUICK;
  } else if (strcmp(propellant_type,"PLAID_AP_HTPB") == 0) {
    rhos   = RHOS_PLAID;
    n      = N_PLAID;
    beta   = BETA_PLAID;
    Tf     = TF_PLAID;
    Ts     = TS_PLAID;
  } else {
    rhos   = RHOS_APHTPB;
    n      = N_APHTPB;
    beta   = BETA_APHTPB;
    Tf     = TF_APHTPB;
    Ts     = TS_APHTPB;
  }
  sigmav = sigv;
  lw = lengthw;
}

/**********************************************************************
 * NavierStokes2D_cState::v -- Gas flow velocity.                     *
 **********************************************************************/
inline Vector2D NavierStokes2D_cState::v(void) const {
  //assert(rho > ZERO);
  return dv/rho;
}

/**********************************************************************
 * NavierStokes2D_cState::v -- Gas flow velocity.                     *
**********************************************************************/
inline double NavierStokes2D_cState::v(const Vector2D &n) const {
  //assert(rho > ZERO);
  return (dv*n)/rho;
}

/**********************************************************************
 * NavierStokes2D_cState::p -- Gas pressure.                          *
 **********************************************************************/
inline double NavierStokes2D_cState::p(void) const {
  //assert(rho > ZERO);
  return gm1*(E - HALF*dv.sqr()/rho - dk);
}

/**********************************************************************
 * NavierStokes2D_cState::T -- Gas temperature.                       *
 **********************************************************************/
inline double NavierStokes2D_cState::T(void) const {
  //assert(rho > ZERO);
  return p()/(rho*R);
}

/**********************************************************************
 * NavierStokes2D_cState::e -- Gas specific internal energy.          *
 **********************************************************************/
inline double NavierStokes2D_cState::e(void) const {
  //assert(rho > ZERO);
  return p()/(gm1*rho);
}

/**********************************************************************
 * NavierStokes2D_cState::h -- Gas specific enthalpy.                 *
 **********************************************************************/
inline double NavierStokes2D_cState::h(void) const {
  //assert(rho > ZERO);
  return g*gm1i*p()/rho + HALF*dv.sqr()/sqr(rho) + k();
}

/**********************************************************************
 * NavierStokes2D_cState::H -- Gas total enthalpy.                    *
 **********************************************************************/
inline double NavierStokes2D_cState::H(void) const {
  return g*gm1i*p() + HALF*dv.sqr()/rho + dk;
}

/**********************************************************************
 * NavierStokes2D_cState::a -- Gas sound speed.                       *
 **********************************************************************/
inline double NavierStokes2D_cState::a(void) const {
  //assert(rho > ZERO);
  return sqrt(g*p()/rho);
}

/**********************************************************************
 * NavierStokes2D_cState::a2 -- Gas sound speed squared.              *
 **********************************************************************/
inline double NavierStokes2D_cState::a2(void) const {
  //assert(rho > ZERO);
  return g*p()/rho;
}

/**********************************************************************
 * NavierStokes2D_cState::M -- Gas Mach number.                       *
 **********************************************************************/
inline double NavierStokes2D_cState::M(void) const {
  //assert(rho > ZERO);
  return abs(v())/a();
}

/**********************************************************************
 * NavierStokes2D_cState::s -- Gas specific entropy.                  *
 **********************************************************************/
inline double NavierStokes2D_cState::s(void) const {
  //assert(rho > ZERO);
  return R*gm1i*log(p()/pow(rho,g));
  //return R*gm1i*log(gm1*(E - HALF*dv.sqr()/rho)/pow(rho,g));
}

/**********************************************************************
 * NavierStokes2D_cState::To -- Gas stagnation temperature.           *
 **********************************************************************/
inline double NavierStokes2D_cState::To(void) const {
  //assert(rho > ZERO);
  return (gm1*(E - HALF*dv.sqr()/rho)/(rho*R))*
	 (ONE+HALF*gm1*dv.sqr()/(rho*rho*g*gm1*(E/rho - HALF*dv.sqr()/sqr(rho))));
}

/**********************************************************************
 * NavierStokes2D_cState::po -- Gas stagnation pressure.              *
 **********************************************************************/
inline double NavierStokes2D_cState::po(void) const {
  return (gm1*(E - HALF*dv.sqr()/rho))*
	  pow(ONE+HALF*gm1*dv.sqr()/(rho*rho*g*gm1*(E/rho - HALF*dv.sqr()/sqr(rho))),g*gm1i);
}

/**********************************************************************
 * NavierStokes2D_cState::ao -- Gas stagnation sound speed.           *
 **********************************************************************/
inline double NavierStokes2D_cState::ao(void) const {
  return sqrt((g*gm1*(E/rho - HALF*dv.sqr()/sqr(rho)))*
	      (ONE+HALF*gm1*dv.sqr()/(rho*rho*g*gm1*(E/rho - HALF*dv.sqr()/sqr(rho)))));
}

/**********************************************************************
 * NavierStokes2D_cState::ho -- Gas stagnation enthalpy.              *
 **********************************************************************/
inline double NavierStokes2D_cState::ho(void) const {
  return (g*E/rho - gm1*HALF*dv.sqr()/sqr(rho))*
         (ONE+HALF*gm1*dv.sqr()/(rho*rho*g*gm1*(E/rho - HALF*dv.sqr()/sqr(rho))));
}

/**********************************************************************
 * NavierStokes2D_cState::mu -- Gas dynamic viscosity.                *
 **********************************************************************/
inline double NavierStokes2D_cState::mu(void) const {
  return mu_gottlieb(v1,v2,v3,v4,v5,T());
}

/**********************************************************************
 * NavierStokes2D_cState::nu -- Gas kinematic viscosity.              *
 **********************************************************************/
inline double NavierStokes2D_cState::nu(void) const {
  return mu()/rho;
}

/**********************************************************************
 * NavierStokes2D_cState::kappa -- Gas thermal heat conductivity.     *
 **********************************************************************/
inline double NavierStokes2D_cState::kappa(void) const {
  return kappa_gottlieb(v1,v2,v3,v4,v5,T(),g,cp);
}

/**********************************************************************
 * NavierStokes2D_cState::Pr -- Prandtl number.                       *
 **********************************************************************/
inline double NavierStokes2D_cState::Pr(void) const {
  return cp*mu()/kappa();
}

/**********************************************************************
 * NavierStokes2D_cState::Alpha -- Thermal diffusivity.               *
 **********************************************************************/
inline double NavierStokes2D_cState::Alpha(void) const {
  return nu()/Pr();
}

/**********************************************************************
 * NavierStokes2D_cState:meanfreepath -- Gas mean free path.          *
 **********************************************************************/
inline double NavierStokes2D_cState::meanfreepath(void) const {
  //assert(rho > ZERO && T() > ZERO);
  return 16.0*mu()/(5.0*rho*sqrt(2.0*PI*R*T()));
}

/**********************************************************************
 * NavierStokes2D_cState::k -- Specific turbulent kinetic energy.     *
 **********************************************************************/
inline double NavierStokes2D_cState::k(void) const {
  return dk/rho;
}

/**********************************************************************
 * NavierStokes2D_cState::depsilon -- Total turbulent eddy            *
 *                                    dissipation.                    *
 **********************************************************************/
inline double NavierStokes2D_cState::depsilon(void) const {
  return beta_k_o*dk*domega/rho;
}

/**********************************************************************
 * NavierStokes2D_cState::epsilon -- Specific turbulent eddy          *
 *                                   dissipation.                     *
 **********************************************************************/
inline double NavierStokes2D_cState::epsilon(void) const {
  return depsilon()/rho;
}

/**********************************************************************
 * NavierStokes2D_cState::omega -- Specific turbulent dissipation     *
 *                                 rate.                              *
 **********************************************************************/
inline double NavierStokes2D_cState::omega(void) const {
  return domega/rho;
}

/**********************************************************************
 * NavierStokes2D_cState::ke -- Internal energy variance rate         *
 **********************************************************************/
inline double NavierStokes2D_cState::ke(void) const {
  return dke/rho;
}

/**********************************************************************
 * NavierStokes2D_cState::ee -- Gas internal energy variance          *
 *                              dissipation rate.                     *   
 **********************************************************************/
inline double NavierStokes2D_cState::ee(void) const {
  return dee/rho;
}

/**********************************************************************
 * NavierStokes2D_cState::ell -- Return the turbulent length scale.   *
 **********************************************************************/
inline double NavierStokes2D_cState::ell(void) const {
  return sqrt(k())/max(omega(),NANO);
}

/**********************************************************************
 * NavierStokes2D_cState::Mt -- Return the turbulent Mach number.     *
 **********************************************************************/
inline double NavierStokes2D_cState::Mt(void) const {
  return sqrt(TWO*k()/a2());
}

/**********************************************************************
 * NavierStokes2D_cState::Mt2 -- Turbulent Mach number squared.       *
 **********************************************************************/
inline double NavierStokes2D_cState::Mt2(void) const {
  return TWO*k()/a2();
}

/**********************************************************************
 * NavierStokes2D_cState::ReT -- Gas turbulent Reynolds number.        *
 **********************************************************************/
inline double NavierStokes2D_cState::ReT(void) const {
  return k()/max(omega(),TOLER)/nu();
}

/**********************************************************************
 * NavierStokes2D_cState::alpha_o_star -- Wilcox's Transition model.  *
 **********************************************************************/
inline double NavierStokes2D_cState:: alpha_o_star(void) const {
  return beta_Wilcox / THREE; 
}

/**********************************************************************
 * NavierStokes2D_cState::beta_star --  Wilcox's Transition model.    *
 **********************************************************************/
inline double NavierStokes2D_cState::beta_star(void) const {
  return NINE/100.0*((FIVE/18.0)+pow(ReT()/Rbeta,4.0))/(ONE+pow(ReT()/Rbeta,4.0)) ;
}

/**********************************************************************
 * NavierStokes2D_cState::alpha --  Wilcox's Transition model.        *
 **********************************************************************/
inline double NavierStokes2D_cState::alpha_Wilcox(void) const {
  return FIVE/NINE * (alpha_o + (ReT()/Romega))/(ONE + (ReT()/Romega))/alpha_star_Wilcox();
}

/**********************************************************************
 * NavierStokes2D_cState::alpha_star --  Wilcox's Transition model.   *
 **********************************************************************/
inline double NavierStokes2D_cState:: alpha_star_Wilcox(void) const{
  return (alpha_o_star() + (ReT())/(Rk))/(ONE + (ReT())/(Rk));

}

/**********************************************************************
 * NavierStokes2D_cState::muT -- Turbulent eddy dynamic viscosity.    *
 **********************************************************************/
inline double NavierStokes2D_cState::muT(void) const {
  return rho*nuT();
}

/**********************************************************************
 * NavierStokes2D_cState::PrT -- Gas turbulent Prandtl number.        *
 **********************************************************************/
inline double NavierStokes2D_cState::PrT(const NavierStokes2D_pState &W,const double &ywall,const double &yplus) const {
  if (Variable_Prandtl == OFF){
    return 0.9;
  }else {
    /************************************************************************************
     * The original formula given in the Calhoon, Brinckman paper (AIAA-2006-1452)      *
     * uses Cmu and Fmu() for calculating the PrT, where Fmu is a damping function used *
     * by the K-epsilon model. We replace the product Cmu*Fmu() by the equivalent in    * 
     * terms of muT,k and epsilon (eqn 11) in the paper.                                *
     ************************************************************************************/
    double CmuFmu = muT()*epsilon()/rho/sqr(k());
    double Vpr;
    //  if(yplus<yplus_sublayer){
    //	Vpr = 0.2;}else{
    Vpr = CmuFmu/C_lambda/max(f_lambda(ywall),TOLER)*sqrt(k()*ee()/max(epsilon(),TOLER)/max(ke(),TOLER));
    // }
    return Vpr;
  }
}

/**********************************************************************
 * NavierStokes2D_cState::alphaT -- Turbulent thermal diffusivity.    *
 **********************************************************************/
inline double NavierStokes2D_cState::alphaT(const NavierStokes2D_pState &W,const double &ywall,const double &yplus) const {
  return nuT()/max(PrT(W,ywall,yplus),TOLER);
}

/**********************************************************************
 * NavierStokes2D_pState::U_e -- Kolmogorov velocity scale.           *
 **********************************************************************/
inline double NavierStokes2D_cState::U_e(void) const {
  return sqrt(sqrt(nu()*epsilon()));
}

/**********************************************************************
 * NavierStokes2D_cState::f_lambda -- Wall damping function(Variable  * 
 *                                     Prandtl number model).         *
 **********************************************************************/
inline double NavierStokes2D_cState::f_lambda(const double &ywall) const {
  double Re_t = sqr(k())/nu()/epsilon();
  double f_et = exp(-sqr(Re_t/80.0));
  double y_star = U_e()*ywall/nu();
  return f_et*C1_lambda/sqrt(sqrt(max(ReT(),TOLER)/beta_k_o))+sqr(ONE-exp(-y_star/A_plus));
}

/**********************************************************************
 * NavierStokes2D_cState::xi_et -- Wall damping function(Variable     * 
 *                                     Prandtl number model).         *
 **********************************************************************/
inline double NavierStokes2D_cState::xi_et(const double &ywall,
					   const NavierStokes2D_pState &dWdx, 
					   const NavierStokes2D_pState &dWdy) const {
  double deriv = (ONE/R/pow(rho,TWO))*(rho*dWdx.p-p()*dWdx.rho);
  double Re_t = sqr(k())/nu()/max(epsilon(),TOLER);
  double f_et = exp(-sqr(Re_t/80.0));
  double Pe_star = -sqrt(TWO/THREE*k()*ke())*deriv;
  double ee_star = ee() - Alpha()*ke()/sqr(ywall);
  double e_cap = epsilon() - TWO*nu()/FOUR/k()*sqr(dWdy.k);
  double ee_cap = ee() - Alpha()/FOUR/max(ke(),TOLER)*sqr(dWdy.ke);
  return f_et*rho*((Cd4-4.0)*ee_cap*ee()/max(ke(),TOLER) + Cd5*e_cap*ee()/max(k(),TOLER) - sqr(ee_star)/max(ke(),TOLER) + (TWO-Cd1-Cd2*Pr())*ee()*Pe_star/max(ke(),TOLER));
}

/**********************************************************************
 * NavierStokes2D_pState::f_mu -- Wall damping function(Variable      * 
 *                                     Prandtl number model).         *
 **********************************************************************/
inline double NavierStokes2D_cState::f_mu(const double &ywall) const {
  double Rek = rho*sqrt(k())*ywall/mu();
  return (ONE + FOUR/pow(ReT()/beta_k_o,THREE/FOUR))*tanh(Rek/125.0);
}

/**********************************************************************
 * NavierStokes2D_cState::nuT -- Turbulent eddy kinematic viscosity.  *
 **********************************************************************/
inline double NavierStokes2D_cState::nuT(void) const {
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA){
    if (Transition_Model == TRANSITION_WILCOX) return alpha_o_star()*k()/max(omega(),TOLER);
    return k()/max(omega(),TOLER);
  }
  return ZERO;
}

/**********************************************************************
 * NavierStokes2D_cState::kappaT -- Turbulent eddy thermal heat       *
 *                                  conductivity.                     *
 **********************************************************************/
inline double NavierStokes2D_cState::kappaT(const NavierStokes2D_pState &W,const double &ywall,const double &yplus) const {
  return muT()*cp/max(PrT(W,ywall,yplus),TOLER);
}

/**********************************************************************
 * NavierStokes2D_cState::c -- Turbulence modified sound speed.       *
 **********************************************************************/
inline double NavierStokes2D_cState::c(void) const {
  return sqrt(c2());
}

/**********************************************************************
 * NavierStokes2D_cState::c2 -- Turbulence modified sound speed       *
 *                              squared.                              *
 **********************************************************************/
inline double NavierStokes2D_cState::c2(void) const {
  return a2() + (2.0/3.0)*g*k();
}

/**********************************************************************
 * NavierStokes2D_cState::pmodified -- Turbulence modified pressure.  *
 **********************************************************************/
inline double NavierStokes2D_cState::pmodified(void) const {
  return p() + (2.0/3.0)*dk;
}

/**********************************************************************
 * NavierStokes2D_cState::beta_k -- k-omega auxilary relation.        *
 **********************************************************************/
inline double NavierStokes2D_cState::beta_k(const NavierStokes2D_pState &dWdx,
					    const NavierStokes2D_pState &dWdy) const {
  double B, xi;
  if (Transition_Model == TRANSITION_WILCOX) B = beta_star();
  else B = beta_k_o*f_beta_k(dWdx,dWdy);

  if (Compressibility_Effect == OFF){
    return B;
  }else if (Compressibility_Effect == COMPRESSIBILITY_CORRECTION_SARKAR){
    xi = ONE;
    return B * ( ONE + xi*F_Sarkar() );
  }else if (Compressibility_Effect == COMPRESSIBILITY_CORRECTION_ZEMAN){
    xi = THREE/FOUR;
    return B * ( ONE + xi*F_Zeman() );
  }else if (Compressibility_Effect == COMPRESSIBILITY_CORRECTION_WILCOX){
    xi = THREE/TWO;
    return B * ( ONE + xi*F_Wilcox() );
  }
}

/**********************************************************************
 * NavierStokes2D_cState::beta_omega -- k-omega auxilary relation.    *
 **********************************************************************/
inline double NavierStokes2D_cState::beta_omega(const NavierStokes2D_pState &dWdx,
						const NavierStokes2D_pState &dWdy) const {
  double B1, B2; 
  if (Transition_Model == TRANSITION_WILCOX){
    B1 = beta_Wilcox;
    B2 = beta_star();
  } else {
    B1 = beta_omega_o*f_beta_omega(dWdx,dWdy);
    B2 = beta_k_o*f_beta_k(dWdx,dWdy);
  }
  
  if (Compressibility_Effect == OFF){
    return B1;
  } else if (Compressibility_Effect == COMPRESSIBILITY_CORRECTION_SARKAR){
    double xi = ONE;
    return B1 - B2*xi*F_Sarkar();
  }else if (Compressibility_Effect == COMPRESSIBILITY_CORRECTION_ZEMAN){
    double xi = THREE/FOUR;
    return B1 - B2*xi*F_Zeman();
  }else if (Compressibility_Effect == COMPRESSIBILITY_CORRECTION_WILCOX){
    double xi = THREE/TWO;
    return B1 - B2*xi*F_Wilcox();
  }
}

/**********************************************************************
 * NavierStokes2D_cState::f_beta_k -- k-omega auxilary relation.      *
 **********************************************************************/
inline double NavierStokes2D_cState::f_beta_k(const NavierStokes2D_pState &dWdx,
					      const NavierStokes2D_pState &dWdy) const {
  double chi = chi_k(dWdx,dWdy);
  if (chi <= ZERO) return ONE;
  return (ONE + 680.0*chi*chi)/(ONE + 400.0*chi*chi);
}

/**********************************************************************
 * NavierStokes2D_cState::f_beta_omega -- k-omega auxilary relation.  *
 **********************************************************************/
inline double NavierStokes2D_cState::f_beta_omega(const NavierStokes2D_pState &dWdx,
						  const NavierStokes2D_pState &dWdy) const {
  double chi = chi_omega(dWdx,dWdy);
  return (ONE + 70.0*chi)/(ONE + 80.0*chi);
}

/***************************************************************************
 * NavierStokes2D_cState::F_Sarkar -- k-omega auxilary relation.           *
 ***************************************************************************/
inline double NavierStokes2D_cState::F_Sarkar(void) const {
  return Mt2() ; 
}

/**************************************************************************
 * NavierStokes2D_cState::F_Zeman -- k-omega auxilary relation.           *
 **************************************************************************/
inline double NavierStokes2D_cState::F_Zeman(void) const {
  double Mt0 = 0.10*sqrt(TWO/(g+1)) ;
  double Omega = 0.60; 
  return (ONE-exp(-HALF*(g+1)*sqr((sqrt(Mt2()) - Mt0)/Omega)))*heaviside(Mt()-Mt0);
}

/***************************************************************************
 * NavierStokes2D_cState::F_Wilcox -- k-omega auxilary relation.           *
 ***************************************************************************/
inline double NavierStokes2D_cState::F_Wilcox(void) const {
  double Mt0 = ONE/FOUR;
  return (Mt2() - sqr(Mt0))*heaviside(Mt()-Mt0);
}

/**********************************************************************
 * NavierStokes2D_cState::chi_k -- k-omega auxilary relation.         *
 **********************************************************************/
inline double NavierStokes2D_cState::chi_k(const NavierStokes2D_pState &dWdx,
					   const NavierStokes2D_pState &dWdy) const {
  //return 0.0;
  return (dWdx.k*dWdx.omega + dWdy.k*dWdy.omega)/max(TOLER,cube(omega()));
}

/**********************************************************************
 * NavierStokes2D_cState::chi_omega -- k-omega auxilary relation.     *
 **********************************************************************/
inline double NavierStokes2D_cState::chi_omega(const NavierStokes2D_pState &dWdx,
					       const NavierStokes2D_pState &dWdy) const {
  //return 0.0;
  return 0.25*fabs(sqr(dWdx.v.y - dWdy.v.x)*(dWdx.v.x + dWdy.v.y)/max(TOLER,cube(beta_omega_o*omega())));
}

/***********************************************************************
 * NavierStokes2D_cState::burningrate -- Solid propellent burning rate.*
 ***********************************************************************/
inline double NavierStokes2D_cState::burningrate(void) const {
  return -beta*pow(p(),n);
}

/**********************************************************************
 * NavierStokes2D_cState -- Binary arithmetic operators.              *
 **********************************************************************/
inline NavierStokes2D_cState NavierStokes2D_cState::operator +(const NavierStokes2D_cState &U) const {
  return NavierStokes2D_cState(rho+U.rho,dv.x+U.dv.x,dv.y+U.dv.y,E+U.E,dk+U.dk,domega+U.domega,dke+U.dke,dee+U.dee);
}

inline NavierStokes2D_cState NavierStokes2D_cState::operator -(const NavierStokes2D_cState &U) const {
  return NavierStokes2D_cState(rho-U.rho,dv.x-U.dv.x,dv.y-U.dv.y,E-U.E,dk-U.dk,domega-U.domega,dke-U.dke,dee-U.dee);
}

// Inner product operator.
inline double NavierStokes2D_cState::operator *(const NavierStokes2D_cState &U) const {
  return rho*U.rho + dv.x*U.dv.x + dv.y*U.dv.y + E*U.E + dk*U.dk + domega*U.domega + dke*U.dke + dee*U.dee;
}

inline NavierStokes2D_cState NavierStokes2D_cState::operator *(const double &a) const {
  return NavierStokes2D_cState(rho*a,dv.x*a,dv.y*a,E*a,dk*a,domega*a,dke*a,dee*a);
}

inline NavierStokes2D_cState NavierStokes2D_cState::operator /(const double &a) const {
  return NavierStokes2D_cState(rho/a,dv.x/a,dv.y/a,E/a,dk/a,domega/a,dke/a,dee/a);
}

// A useful solution state product operator.
inline NavierStokes2D_cState NavierStokes2D_cState::operator ^(const NavierStokes2D_cState &U) const {
  return NavierStokes2D_cState(rho*U.rho,dv.x*U.dv.x,dv.y*U.dv.y,E*U.E,dk*U.dk,domega*U.domega,dke*U.dke,dee*U.dee);
}

/**********************************************************************
 * NavierStokes2D_cState -- Assignment operator.                      *
 **********************************************************************/
inline NavierStokes2D_cState& NavierStokes2D_cState::operator =(const NavierStokes2D_cState &U) {
  // handle self-assignment
  if (this == &U) return *this;

  rho = U.rho; dv.x = U.dv.x; dv.y = U.dv.y; E = U.E;
  dk = U.dk; domega = U.domega; dke = U.dke; dee = U.dee;

  return *this;
}

/**********************************************************************
 * NavierStokes2D_cState -- Unary arithmetic operators.               *
 **********************************************************************/
inline NavierStokes2D_cState operator +(const NavierStokes2D_cState &U) {
  return U;
}

inline NavierStokes2D_cState operator -(const NavierStokes2D_cState &U) {
  return NavierStokes2D_cState(-U.rho,-U.dv.x,-U.dv.y,-U.E,-U.dk,-U.domega,-U.dke,-U.dee);
}

inline NavierStokes2D_cState fabs(const NavierStokes2D_cState &U){
  return NavierStokes2D_cState(fabs(U.rho),fabs(U.dv.x),fabs(U.dv.y),fabs(U.E),
			       fabs(U.dk),fabs(U.domega),fabs(U.dke),fabs(U.dee));
}

inline NavierStokes2D_cState sqr(const NavierStokes2D_cState &U){
  return NavierStokes2D_cState(sqr(U.rho),sqr(U.dv.x),sqr(U.dv.y),sqr(U.E),
			       sqr(U.dk),sqr(U.domega),sqr(U.dke),sqr(U.dee));
}


/**********************************************************************
 * NavierStokes2D_cState -- Shortcut arithmetic operators.            *
 **********************************************************************/
inline NavierStokes2D_cState& NavierStokes2D_cState::operator +=(const NavierStokes2D_cState &U) {
  rho += U.rho; dv.x += U.dv.x; dv.y += U.dv.y; E += U.E; dk += U.dk; domega += U.domega; dke += U.dke; dee += U.dee;
  return *this;
}

inline NavierStokes2D_cState& NavierStokes2D_cState::operator -=(const NavierStokes2D_cState &U) {
  rho -= U.rho; dv.x -= U.dv.x; dv.y -= U.dv.y; E -= U.E; dk -= U.dk; domega -= U.domega; dke -= U.dke; dee -= U.dee;
  return *this;
}

inline NavierStokes2D_cState& NavierStokes2D_cState::operator *=(const double &a) {
  rho *= a; dv.x *= a; dv.y *= a; E *= a; dk *= a; domega *= a; dke *= a; dee *= a;
  return *this;
}

inline NavierStokes2D_cState& NavierStokes2D_cState::operator /=(const double &a) {
  rho /= a; dv.x /= a; dv.y /= a; E /= a; dk /= a; domega /= a; dke /= a; dee /= a;
  return *this;
}

/**********************************************************************
 * NavierStokes2D_cState -- Relational operators.                     *
 **********************************************************************/
inline int operator ==(const NavierStokes2D_cState &U1, const NavierStokes2D_cState &U2) {
  return (U1.rho == U2.rho && U1.dv == U2.dv && U1.E == U2.E && U1.dk == U2.dk && 
	  U1.domega == U2.domega && U1.dke == U2.dke && U1.dee == U2.dee);
}

inline int operator !=(const NavierStokes2D_cState &U1, const NavierStokes2D_cState &U2) {
  return (U1.rho != U2.rho || U1.dv != U2.dv || U1.E != U2.E || U1.dk != U2.dk || 
	  U1.domega != U2.domega || U1.dke != U1.dke || U1.dee != U2.dee);
}

inline bool operator >=(const NavierStokes2D_cState& U1, const NavierStokes2D_cState& U2){
  return (U1.rho >= U2.rho && U1.dv.x >= U2.dv.x && U1.dv.y >= U2.dv.y && U1.E >= U2.E && 
	  U1.dk >= U2.dk && U1.domega >= U2.domega && U1.dke >= U1.dke && U1.dee >= U2.dee);
}

inline bool operator <=(const NavierStokes2D_cState& U1, const NavierStokes2D_cState& U2){
  return (U1.rho <= U2.rho && U1.dv.x <= U2.dv.x && U1.dv.y <= U2.dv.y && U1.E <= U2.E && 
	  U1.dk <= U2.dk && U1.domega <= U2.domega && U1.dke <= U1.dke && U1.dee <= U2.dee);
}

inline bool operator <(const NavierStokes2D_cState &U1, const NavierStokes2D_cState &U2){
  return (U1.rho < U2.rho && U1.dv.x < U2.dv.x && U1.dv.y < U2.dv.y && U1.E < U2.E &&
	  U1.dk < U2.dk && U1.domega < U2.domega && U1.dke < U1.dke && U1.dee < U2.dee);  
}

inline bool operator >(const NavierStokes2D_cState &U1, const NavierStokes2D_cState &U2){
  return (U1.rho > U2.rho && U1.dv.x > U2.dv.x && U1.dv.y > U2.dv.y && U1.E > U2.E &&
	  U1.dk > U2.dk && U1.domega > U2.domega && U1.dke > U1.dke && U1.dee > U2.dee);
}


/**********************************************************************
 * NavierStokes2D_cState -- Input-output operators.                   *
 **********************************************************************/
inline ostream &operator << (ostream &out_file, const NavierStokes2D_cState &U) {
  out_file.setf(ios::scientific);
  out_file << " " << U.rho << " " << U.dv.x << " " << U.dv.y << " " << U.E 
	   << " " << U.dk << " " << U.domega << " " << U.dke << " " << U.dee;
  out_file.unsetf(ios::scientific);
  return out_file;
}

inline istream &operator >> (istream &in_file, NavierStokes2D_cState &U) {
  in_file.setf(ios::skipws);
  in_file >> U.rho >> U.dv.x >> U.dv.y >> U.E >> U.dk >> U.domega >> U.dke >> U.dee;
  in_file.unsetf(ios::skipws);
  return in_file;
}

/**********************************************************************
 * NavierStokes2D_pState::NavierStokes2D_pState -- Constructor.       *
 **********************************************************************/
inline NavierStokes2D_pState::NavierStokes2D_pState(const NavierStokes2D_cState &U) {
  rho = U.rho; v.x = U.v().x; v.y = U.v().y; p = U.p(); k = U.k(); omega = U.omega(); ke = U.ke(); ee = U.ee();
}

/**********************************************************************
 * NavierStokes2D_pState::U -- Conserved solution state.              *
 **********************************************************************/
inline NavierStokes2D_cState NavierStokes2D_pState::U(void) const {
  return NavierStokes2D_cState(rho,dv(),E(),dk(),domega(),dke(),dee());
}

inline NavierStokes2D_cState NavierStokes2D_pState::U(const NavierStokes2D_pState &W) const {
  return W.U();
}

inline NavierStokes2D_cState U(const NavierStokes2D_pState &W) {
  return W.U();
}

/**********************************************************************
 * NavierStokes2D_pState::dUdW -- Jacobian of the conserved solution  *
 *                                variables with respect to the       *
 *                                primitive solution variables .      *
 **********************************************************************/
inline void NavierStokes2D_pState::dUdW(DenseMatrix &dUdW) const {
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
    if (Variable_Prandtl == ON) {
      dUdW(6,0) += ke;
      dUdW(6,6) += rho;
      dUdW(7,0) += ee;
      dUdW(7,7) += rho;
    }
  }
}

/**********************************************************************
 * NavierStokes2D_pState::dWdU -- Jacobian of the primitive solution  *
 *                                variables with respect to the       *
 *                                conserved solution variables.       *
 **********************************************************************/
inline void NavierStokes2D_pState::dWdU(DenseMatrix &dWdU) const {
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
    if(Variable_Prandtl == ON) {
      dWdU(6,0) -= ke/rho;
      dWdU(6,6) += ONE/rho;
      dWdU(7,0) -= ee/rho;
      dWdU(7,7) += ONE/rho;
    }
  }
}

/**********************************************************************
 * NavierStokes2D_pState::F -- Solution inviscid flux (x-direction).  *
 **********************************************************************/
inline NavierStokes2D_cState NavierStokes2D_pState::F(void) const {
  if (Variable_Prandtl == ON){
  return NavierStokes2D_cState(rho*v.x,
			       rho*sqr(v.x)+p+(2.0/3.0)*dk(),
			       rho*v.x*v.y,
			       v.x*H()+v.x*(2.0/3.0)*dk(),
			       rho*v.x*k,
			       rho*v.x*omega,
			       rho*v.x*ke,
			       rho*v.x*ee);
  }else{
    return NavierStokes2D_cState(rho*v.x,
				 rho*sqr(v.x)+p+(2.0/3.0)*dk(),
				 rho*v.x*v.y,
				 v.x*H()+v.x*(2.0/3.0)*dk(),
				 rho*v.x*k,
				 rho*v.x*omega,
				 ZERO,
				 ZERO);
  }
}
  
inline NavierStokes2D_cState NavierStokes2D_pState::F(const Vector2D &V) const {
  if (Variable_Prandtl == ON){
  return NavierStokes2D_cState(rho*(v.x-V.x),
			       rho*(v.x-V.x)*v.x+p+(2.0/3.0)*dk(),
			       rho*(v.x-V.x)*v.y,
			       (v.x-V.x)*E()+ v.x*(p+(2.0/3.0)*dk()),
			       rho*(v.x-V.x)*k,
			       rho*(v.x-V.x)*omega,
			       rho*(v.x-V.x)*ke,
			       rho*(v.x-V.x)*ee);
  }else{
  return NavierStokes2D_cState(rho*(v.x-V.x),
			       rho*(v.x-V.x)*v.x+p+(2.0/3.0)*dk(),
			       rho*(v.x-V.x)*v.y,
			       (v.x-V.x)*E()+ v.x*(p+(2.0/3.0)*dk()),
			       rho*(v.x-V.x)*k,
			       rho*(v.x-V.x)*omega,
			       ZERO, ZERO);
  }
}

/*
 * Calculate the inviscid flux in the normal direction 
 * based on the current solution state.
 *
 * \param normal_dir vector defining the normal direction
 */
inline NavierStokes2D_cState NavierStokes2D_pState::Fn(const Vector2D &normal_dir) const {
  double V_dot_n(v.x*normal_dir.x + v.y*normal_dir.y);  
  double k_Term(0.66666666666666666667*dk()); // Turbulent kinetic energy term, (0.66666666666666666667 = 2/3)

  if (Variable_Prandtl == ON){
    return NavierStokes2D_cState(rho*V_dot_n,
				 rho*v.x*V_dot_n + (p + k_Term)*normal_dir.x,
				 rho*v.y*V_dot_n + (p + k_Term)*normal_dir.y,
				 V_dot_n*H() + V_dot_n*k_Term,
				 rho*V_dot_n*k,
				 rho*V_dot_n*omega,
				 rho*V_dot_n*ke,
				 rho*V_dot_n*ee);
  }else{
    return NavierStokes2D_cState(rho*V_dot_n,
				 rho*v.x*V_dot_n + (p + k_Term)*normal_dir.x,
				 rho*v.y*V_dot_n + (p + k_Term)*normal_dir.y,
				 V_dot_n*H() + V_dot_n*k_Term,
				 rho*V_dot_n*k,
				 rho*V_dot_n*omega,
				 ZERO,
				 ZERO);
  }
}

/*
 * Calculate the inviscid flux in the normal direction 
 * based on the given solution state.
 */
inline NavierStokes2D_cState Fn(const NavierStokes2D_pState &W, const Vector2D & normal_dir){ 
  return W.Fn(normal_dir); 
}

/**********************************************************************
 * NavierStokes2D_pState::dFdU -- Jacobian of the inviscid solution   *
 *                                flux with respect to the conserved  *
 *                                solution variables.                 *
 **********************************************************************/
inline void NavierStokes2D_pState::dFdU(DenseMatrix &dFdU) const {
  dFdU(0,1) += ONE;
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
    if(Variable_Prandtl == ON) {    
      dFdU(6,0) -= v.x*ke;
      dFdU(6,1) += ke;
      dFdU(6,6) += v.x;
      dFdU(7,0) -= v.x*ee;
      dFdU(7,1) += ee;
      dFdU(7,7) += v.x;
    }
  }
}

/**********************************************************************
 * NavierStokes2D_pState::dFdW -- Jacobian of the inviscid solution       *
 *                            flux with respect to the primitive      *
 *                            solution variables.                     *
 **********************************************************************/
inline void NavierStokes2D_pState::dFdW(DenseMatrix &dFdW) const {
  dFdW(0,0) += v.x;
  dFdW(0,1) += rho;
  dFdW(1,0) += v.x*v.x;
  dFdW(1,1) += TWO*rho*v.x; 
  dFdW(1,3) += ONE;
  dFdW(2,0) += v.x*v.y;
  dFdW(2,1) += rho*v.y;
  dFdW(2,2) += rho*v.x;
  dFdW(3,0) += HALF*(v.x*v.x+v.y*v.y)*v.x;
  dFdW(3,1) += rho*v.x*v.x+rho*(g*p/rho/gm1+HALF*(v.x*v.x+v.y*v.y));
  dFdW(3,2) += rho*v.x*v.y;
  dFdW(3,3) += v.x*g/gm1;
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    dFdW(1,0) += TWO/THREE*k;
    dFdW(1,4) += TWO/THREE*rho;
    dFdW(3,0) += FIVE/THREE*k*v.x;
    dFdW(3,1) += FIVE/THREE*rho*k; 
    dFdW(3,4) += FIVE/THREE*rho*v.x;
    dFdW(4,0) += v.x*k;
    dFdW(4,1) += rho*k;  
    dFdW(5,0) += v.x*omega;
    dFdW(5,1) += rho*omega;
  }
  dFdW(4,4) += rho*v.x;
  dFdW(5,5) += rho*v.x;
}

/**********************************************************************
 * NavierStokes2D_pState::Gx, Gy -- Solution viscous fluxes.          *
 **********************************************************************/
inline NavierStokes2D_cState NavierStokes2D_pState::Gx(const NavierStokes2D_pState &dWdx,
						       const double &ywall,const double &yplus) const {
  if (Variable_Prandtl == ON){
    return NavierStokes2D_cState(ZERO,
				 tau.xx,
				 tau.xy,
				 -q.x+v.x*tau.xx+v.y*tau.xy+(mu()+sigma_k*muT())*dWdx.k,
				 (mu()+sigma_k*muT())*dWdx.k,
				 (mu()+sigma_omega*muT())*dWdx.omega,
				 rho*(Alpha()+alphaT(ywall,yplus)/sigma_k_e)*dWdx.ke,
				 rho*(Alpha()+alphaT(ywall,yplus)/sigma_ep_e)*dWdx.ee);
  }else{
    return NavierStokes2D_cState(ZERO,
				 tau.xx,
				 tau.xy,
				 -q.x+v.x*tau.xx+v.y*tau.xy+(mu()+sigma_k*muT())*dWdx.k,
				 (mu()+sigma_k*muT())*dWdx.k,
				 (mu()+sigma_omega*muT())*dWdx.omega,
				 ZERO, 
				 ZERO);
  }
}

inline NavierStokes2D_cState NavierStokes2D_pState::Gy(const NavierStokes2D_pState &dWdy, 
						       const double &ywall,const double &yplus) const {
  if (Variable_Prandtl == ON){
    return NavierStokes2D_cState(ZERO,
				 tau.xy,
				 tau.yy,
				 -q.y+v.x*tau.xy+v.y*tau.yy+(mu()+sigma_k*muT())*dWdy.k,
				 (mu()+sigma_k*muT())*dWdy.k,
				 (mu()+sigma_omega*muT())*dWdy.omega,
				 rho*(Alpha()+alphaT(ywall,yplus)/sigma_k_e)*dWdy.ke,
				 rho*(Alpha()+alphaT(ywall,yplus)/sigma_ep_e)*dWdy.ee);
  }else{
    return NavierStokes2D_cState(ZERO,
				 tau.xy,
				 tau.yy,
				 -q.y+v.x*tau.xy+v.y*tau.yy+(mu()+sigma_k*muT())*dWdy.k,
				 (mu()+sigma_k*muT())*dWdy.k,
				 (mu()+sigma_omega*muT())*dWdy.omega,
				 ZERO, 
				 ZERO);}
}

/**********************************************************************
 * NavierStokes2D_pState::ComputeViscousTerms -- Compute viscous      *
 *                                               stress tensor and    *
 *                                               heat flux vector.    *
 **********************************************************************/
inline void NavierStokes2D_pState::ComputeViscousTerms(const NavierStokes2D_pState &dWdx,
						       const NavierStokes2D_pState &dWdy,
						       const Vector2D &X,
						       const int &Axisymmetric,
						       const int &adiabatic_flag,
						       const double &ywall,
						       const double &yplus) {

  double div, radius, mumu, kap;

  // Total (i.e. molecular + turbulent) dynamic viscosity and thermal conductivity
  mumu = mu() + muT();
  kap = kappa() + kappaT(ywall,yplus);


  if (Axisymmetric){
    radius = max(X.y,TOLER);
  }

  // Divergence of the velocity field.
  div = dWdx.v.x + dWdy.v.y;
  if (Axisymmetric){
    div += v.y/radius;
  }

  // Stress tensor.
  div /= 3.0;	// divide divergence by three
  tau.xx = 2.0*mumu*(dWdx.v.x - div);
  tau.xy = mumu*(dWdy.v.x + dWdx.v.y);
  tau.yy = 2.0*mumu*(dWdy.v.y - div);
  if (Axisymmetric){
    tau.zz = 2.0*mumu*(v.y/radius - div);
  } else {
    tau.zz = ZERO;
  }

  // Heat flux
  q.x = -kap*(dWdx.p - (p/rho)*dWdx.rho)/(rho*R);
  q.y = -kap*(dWdy.p - (p/rho)*dWdy.rho)/(rho*R);

}

/**********************************************************************
 * NavierStokes2D_pState::lambda_x -- Eigenvalue(s) (x-direction).    *
 **********************************************************************/
inline NavierStokes2D_pState NavierStokes2D_pState::lambda_x(void) const {
  return NavierStokes2D_pState(v.x-c(),v.x,v.x,v.x+c(),v.x,v.x,v.x,v.x);
}

inline NavierStokes2D_pState NavierStokes2D_pState::lambda_x(const Vector2D &V) const {
  return NavierStokes2D_pState(v.x-V.x-c(),v.x-V.x,v.x-V.x,v.x-V.x+c(),v.x-V.x,v.x-V.x,v.x -V.x,v.x-V.x);
}

/**********************************************************************
 * NavierStokes2D_pState::rp_x -- Primitive right eigenvector         *
 *                                (x-direction).                      *
 **********************************************************************/
inline NavierStokes2D_pState NavierStokes2D_pState::rp_x(int index) const {
  //assert(index >= 1 && index <= NUM_VAR_NAVIERSTOKES2D);
  switch(index) {
  case 1 :
    return NavierStokes2D_pState(ONE,-c()/rho,ZERO,c2()-(2.0/3.0)*k,ZERO,ZERO,ZERO,ZERO);
    break;
  case 2 :
    return NavierStokes2D_pState(ONE,ZERO,ZERO,-(2.0/3.0)*k,ZERO,ZERO,ZERO,ZERO);
    break;
  case 3 :
    return NavierStokes2D_pState(ZERO,ZERO,ONE,ZERO,ZERO,ZERO,ZERO,ZERO);
    break;
  case 4 :
    return NavierStokes2D_pState(ONE,c()/rho,ZERO,c2()-(2.0/3.0)*k,ZERO,ZERO,ZERO,ZERO);
    break;
  case 5 :
    return NavierStokes2D_pState(ZERO,ZERO,ZERO,-(2.0/3.0)*rho,ONE,ZERO,ZERO,ZERO);
    break;
  case 6 :
    return NavierStokes2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO,ZERO);
    break;
  case 7 :
    return NavierStokes2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO);
    break;
  case 8 :
    return NavierStokes2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE);
    break;
  };
  return NavierStokes2D_pState(ONE,-c()/rho,ZERO,c2()-(2.0/3.0)*k,ZERO,ZERO,ZERO,ZERO);
}

/**********************************************************************
 * NavierStokes2D_pState::rc_x -- Conserved right eigenvector         *
 *                                (x-direction).                      *
 **********************************************************************/
inline NavierStokes2D_cState NavierStokes2D_pState::rc_x(int index) const {
  //assert(index >= 1 && index <= NUM_VAR_NAVIERSTOKES2D);
  switch(index) {
  case 1 :
    return NavierStokes2D_cState(ONE,v.x-c(),v.y,h()-c()*v.x+(2.0/3.0)*k,k,omega,ke,ee);
    break;
  case 2 :
    return NavierStokes2D_cState(ONE,v.x,v.y,HALF*v.sqr()+gm1i*(g-5.0/3.0)*k,k,omega,ke,ee);
    break;
  case 3 :
    return NavierStokes2D_cState(ZERO,ZERO,rho,rho*v.y,ZERO,ZERO,ZERO,ZERO);
    break;
  case 4 :
    return NavierStokes2D_cState(ONE,v.x+c(),v.y,h()+c()*v.x+(2.0/3.0)*k,k,omega,ke,ee);
    break;
  case 5 :
    return NavierStokes2D_cState(ZERO,ZERO,ZERO,rho*gm1i*(g-5.0/3.0),rho,ZERO,ZERO,ZERO);
    break;
  case 6 :
    return NavierStokes2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,rho,ZERO,ZERO);
    break;
  case 7 :
    return NavierStokes2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,rho,ZERO);
    break;
  case 8 :
    return NavierStokes2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,rho);
    break;
  };
  return NavierStokes2D_cState(ONE,v.x-c(),v.y,h()-c()*v.x+(2.0/3.0)*k,k,omega,ke,ee);
}

/**********************************************************************
 * NavierStokes2D_pState::lp_x -- Primitive left eigenvector          *
 *                                (x-direction).                      *
 **********************************************************************/
inline NavierStokes2D_pState NavierStokes2D_pState::lp_x(int index) const {
  //assert(index >= 1 && index <= NUM_VAR_NAVIERSTOKES2D);
  switch(index) {
  case 1 :
    return NavierStokes2D_pState(k/(3.0*c2()),-HALF*rho/c(),ZERO,HALF/c2(),rho/(3.0*c2()),ZERO,ZERO,ZERO);
    break;
  case 2 :
    return NavierStokes2D_pState(ONE-(2.0/3.0)*k/c2(),ZERO,ZERO,-ONE/c2(),-(2.0/3.0)*rho/c2(),ZERO,ZERO,ZERO);
    break;
  case 3 :
    return NavierStokes2D_pState(ZERO,ZERO,ONE,ZERO,ZERO,ZERO,ZERO,ZERO);
    break;
  case 4 :
    return NavierStokes2D_pState(k/(3.0*c2()),HALF*rho/c(),ZERO,HALF/c2(),rho/(3.0*c2()),ZERO,ZERO,ZERO);
    break;
  case 5 :
    return NavierStokes2D_pState(ZERO,ZERO,ZERO,ZERO,ONE,ZERO,ZERO,ZERO);
    break;
  case 6 :
    return NavierStokes2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO,ZERO);
    break;
  case 7 :
    return NavierStokes2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO);
    break;
  case 8 :
    return NavierStokes2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE);
    break;
  };
  return NavierStokes2D_pState(k/(3.0*c2()),-HALF*rho/c(),ZERO,HALF/c2(),rho/(3.0*c2()),ZERO,ZERO,ZERO);
}

/**********************************************************************
 * NavierStokes2D_pState::S -- Include all source term vectors and    *
 *                             Jacobians.                             *
 **********************************************************************/
inline NavierStokes2D_cState NavierStokes2D_pState::S(const Vector2D &X,
                                                      const NavierStokes2D_pState &W,  
						      const NavierStokes2D_pState &dWdx,
						      const NavierStokes2D_pState &dWdy,
						      const int &Axisymmetric,
						      const double &ywall,
						      const double &yplus) const {
  NavierStokes2D_cState Sall; Sall.Vacuum();
  // Include the axisymmetric source terms if required.
  if (Axisymmetric) {
    Sall = Si(X);
    if (flow_type) Sall += Sv(X,dWdy,ywall,yplus);
  }

  // Include the turbulence model source term if required.
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) Sall += St(X,W,dWdx,dWdy,Axisymmetric,ywall,yplus);
  // Return the total source term vector.
  return Sall; 
}

inline void NavierStokes2D_pState::dSdU(DenseMatrix &dSdU,
					const Vector2D &X,
					const NavierStokes2D_pState &dWdx,
					const NavierStokes2D_pState &dWdy,
					const int &Axisymmetric) const {
  // Include the axisymmetric source Jacobians.
  if (Axisymmetric) {
    dSidU(dSdU,X);
    if (flow_type) dSvdU(dSdU,X,dWdy);
  }
  // Include the turbulence model source term if required.
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) dStdU(dSdU,X,dWdx,dWdy,Axisymmetric);
}

/**********************************************************************
 * NavierStokes2D_pState::Si -- Inviscid axisymmetric source terms    *
 *                              and Jacobian.                         * 
 *                              The axisymmetric convention:          * 
 *                              V.y = radial velocity                 * 
 *                              V.x = axial velocity                  * 
 *                              X.y = Radial direction                *
 **********************************************************************/
inline NavierStokes2D_cState NavierStokes2D_pState::Si(const Vector2D &X) const {
  return NavierStokes2D_cState(-rho*v.y/X.y,
			       -rho*v.x*v.y/X.y,
			       -rho*sqr(v.y)/X.y,
			       -v.y*(H()+(2.0/3.0)*dk())/X.y,
			       -v.y*dk()/X.y,
			       -v.y*domega()/X.y,
			       -v.y*dke()/X.y,
			       -v.y*dee()/X.y);
}

inline void NavierStokes2D_pState::dSidU(DenseMatrix &dSidU, const Vector2D &X) const {
  dSidU(0,2) -= ONE/X.y;
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
    if (Variable_Prandtl == ON) {
      dSidU(6,0) += ke*v.y/X.y;
      dSidU(6,2) -= ke/X.y;
      dSidU(6,6) -= v.y/X.y;
      dSidU(7,0) += ee*v.y/X.y;
      dSidU(7,2) -= ee/X.y;
      dSidU(7,7) -= v.y/X.y; 
    }
  }
}

/**********************************************************************
 * NavierStokes2D_pState::Sv -- Viscous axisymmetric flow source term *
 *                              vector and Jacobian.                  *
 **********************************************************************/
inline NavierStokes2D_cState NavierStokes2D_pState::Sv(const Vector2D &X,
						       const NavierStokes2D_pState &dWdy,
						       const double &ywall,
						       const double &yplus) const {
  if (Variable_Prandtl == ON){
  return NavierStokes2D_cState(ZERO,
			       tau.xy/X.y,
			       (tau.yy-tau.zz)/X.y,
			       (-q.y+v.x*tau.xy +v.y*tau.yy+(mu()+sigma_k*muT())*dWdy.k)/X.y,
			       (mu()+sigma_k*muT())*dWdy.k/X.y,
			       (mu()+sigma_omega*muT())*dWdy.omega/X.y,
			       rho*(Alpha()+alphaT(ywall,yplus)/sigma_k_e)*dWdy.ke/X.y,
			       rho*(Alpha()+alphaT(ywall,yplus)/sigma_ep_e)*dWdy.ee/X.y);
  }else{
    return NavierStokes2D_cState(ZERO,
				 tau.xy/X.y,
				 (tau.yy-tau.zz)/X.y,
				 (-q.y+v.x*tau.xy +v.y*tau.yy+(mu()+sigma_k*muT())*dWdy.k)/X.y,
				 (mu()+sigma_k*muT())*dWdy.k/X.y,
				 (mu()+sigma_omega*muT())*dWdy.omega/X.y,
				 ZERO, ZERO);
  }
}

inline void NavierStokes2D_pState::dSvdU(DenseMatrix &dSvdU,
					 const Vector2D &X,
					 const NavierStokes2D_pState &dWdy) const {
  dSvdU(0,0) += ZERO;
  dSvdU(1,1) -= ZERO;
}

/**********************************************************************
 * NavierStokes2D_pState::St -- Turbulent source term vector and      *
 *                              Jacobian.                             *
 **********************************************************************/
inline NavierStokes2D_cState NavierStokes2D_pState::St(const Vector2D &X,
						       const NavierStokes2D_pState &W,
						       const NavierStokes2D_pState &dWdx,
						       const NavierStokes2D_pState &dWdy,
						       const int &Axisymmetric,
						       const double &ywall,
						       const double &yplus) const {
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
  if (Variable_Prandtl == ON){
    double production_cap, deriv2;
    double alp1=2.5, alp2=2.0, lam=0.2;
    double Mt_cap = max(Mt()-lam,ZERO);
    //deriv2 is the square of the derivative of the specific internal energy
    deriv2 = sqr(gm1i/sqr(rho))*(sqr(rho*dWdx.p-p*dWdx.rho) + sqr(rho*dWdy.p-p*dWdy.rho));
    production_cap = production; //-alp1*sqr(Mt_cap)*production-alp2*sqr(Mt_cap)*rho*epsilon();
    return NavierStokes2D_cState(ZERO,ZERO,ZERO,ZERO,
				 production-beta_k(dWdx,dWdy)*dk()*omega,
				 alpha*(omega/max(k,TOLER))*production-beta_omega(dWdx,dWdy)*rho*sqr(omega),
				 productionK(dWdx,dWdy,ywall,yplus) - TWO*rho*ee,
				 productionE1(dWdx,dWdy,ywall,yplus) + productionE2(production_cap) - (Cd4*ee/max(ke,TOLER) + Cd5*epsilon()/max(k,TOLER))*rho*ee );//+ xi_et(ywall,dWdx,dWdy));
  }else{
    return NavierStokes2D_cState(ZERO,ZERO,ZERO,ZERO,
				 production-beta_k(dWdx,dWdy)*dk()*omega,
				 alpha*(omega/max(k,TOLER))*production-beta_omega(dWdx,dWdy)*rho*sqr(omega),
				 ZERO, ZERO);
  }
}

inline void NavierStokes2D_pState::dStdU(DenseMatrix &dStdU,
					 const Vector2D &X,
					 const NavierStokes2D_pState &dWdx,
					 const NavierStokes2D_pState &dWdy,
					 const int &Axisymmetric) const {
  dStdU(0,0) += ZERO;
  dStdU(1,1) -= ZERO;
}

inline void NavierStokes2D_pState::dSvpdU(DenseMatrix &dSvpdU,
					  const Vector2D &X,
					  const NavierStokes2D_pState &dWdx,
					  const NavierStokes2D_pState &dWdy,
					  const double &d_dWdx_dW, 
					  const double &d_dWdy_dW,
					  const int &Axisymmetric,
					  const double &ywall,
					  const double &yplus) const {

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

  /***************************************************************
   * In the production terms,                                    *   
   *                                                             *   
   *      / de \2  / 1    1  \2 /    d p     d rho \2            *
   *      | -- | = | -  * -  |  |rho --- - p ----- |             *
   *      \ dx /   \rho  g-1 /  \    d x     d x   /             *
   *                                                             *
   *                                                             *
   *                                                             * 
   * In the following expressions:                               *
   *                                                             *
   *          /    d p     d rho \2                              *
   *    fx =  |rho --- - p ----- |                               *
   *          \    d x     d x   /                               *
   *                                                             *
   *          d fx              d fy                             *
   *  DfxDp = ----- , DfxDrho = -----                            * 
   *          d p               d rho                            *
   *                                                             *
   * and similiar notation for all the Y derivatives.            *
   *                                                             *
   * A, B,'coeff' and 'coeff1' are just to simplify expressions  * 
   * as some terms appear in every term.                         *
   ***************************************************************/

  double A = alphaT(ywall,yplus)*Cd1*sqr(gm1i);   
  double B = TWO*alphaT(ywall,yplus)*sqr(gm1i);
  double fx = sqr(rho*dWdx.p-p*dWdx.rho); 
  double fy = sqr(rho*dWdy.p-p*dWdy.rho);
  double DfxDp = TWO*sqrt(fx)*(rho*d_dWdx_dW - dWdx.rho);
  double DfyDp = TWO*sqrt(fy)*(rho*d_dWdy_dW - dWdy.rho);
  double DfxDrho = TWO*sqrt(fx)*(dWdx.p - p*d_dWdx_dW);
  double DfyDrho = TWO*sqrt(fy)*(dWdy.p - p*d_dWdy_dW);
  double coeff = A*ee/pow(rho,3.0)/max(ke,TOLER);
  double coeff1 = B/pow(rho,3.0);

  dSvpdU(6,0) += coeff1/rho*(-THREE*(fx+fy)+rho*((DfxDrho+DfyDrho)+gm1*v.sqr()*(DfxDp+DfyDp)/TWO));
  dSvpdU(6,1) -= coeff1*v.x*gm1*(DfxDp+DfyDp);
  dSvpdU(6,2) -= coeff1*v.y*gm1*(DfxDp+DfyDp);
  dSvpdU(6,3) += coeff1*gm1*(DfxDp+DfyDp);
  dSvpdU(6,4) -= coeff1*gm1*(DfxDp+DfyDp);
  dSvpdU(6,5) += ZERO;
  dSvpdU(6,6) += ZERO;
  dSvpdU(6,7) -= TWO;

  dSvpdU(7,0) += coeff/rho*(-THREE*(fx+fy)+rho*((DfxDrho+DfyDrho)+gm1*v.sqr()*(DfxDp+DfyDp)))+Cd5*beta_k_o*omega*ee;
  dSvpdU(7,1) -= coeff*v.x*gm1*(DfxDp+DfyDp);
  dSvpdU(7,2) -= coeff*v.y*gm1*(DfxDp+DfyDp);
  dSvpdU(7,3) += coeff*gm1*(DfxDp+DfyDp);
  dSvpdU(7,4) -= (coeff*gm1*(DfxDp+DfyDp) + Cd3*production*ee/rho/max(sqr(k),TOLER));
  dSvpdU(7,5) -= Cd5*beta_k_o*ee;
  dSvpdU(7,6) -= coeff*(fx+fy)/rho/max(ke,TOLER) + Cd4*sqr(ee/max(ke,TOLER));
  dSvpdU(7,7) += A*(fx + fy)/sqr(sqr(rho))/max(ke,TOLER) - TWO*Cd4*ee/max(ke,TOLER)-Cd5*beta_k_o*omega + Cd3*production/rho/max(k,TOLER); 
}


/**********************************************************************
 * NavierStokes2D_pState::deriv2 -- derivative square.                *
 **********************************************************************/

inline double NavierStokes2D_pState::deriv2(const NavierStokes2D_pState &dWdx,
					    const NavierStokes2D_pState &dWdy) const {
  return sqr(gm1i/sqr(rho))*(sqr(rho*dWdx.p-p*dWdx.rho) + sqr(rho*dWdy.p-p*dWdy.rho));  
}

/**********************************************************************
 * NavierStokes2D_pState::diff -- difference between P and D terms.   *
 **********************************************************************/

inline double NavierStokes2D_pState::diff(const NavierStokes2D_pState &dWdx,
					  const NavierStokes2D_pState &dWdy,
					  const double &ywall,
					  const double &yplus) const {
  return productionE1(dWdx,dWdy,ywall, yplus) - D1() - D2(); 
}

/**********************************************************************
 * NavierStokes2D_pState::D1 -- ee destruction term.                 *
 **********************************************************************/

inline double NavierStokes2D_pState::D1(void) const {
  return Cd4*ee/max(ke,TOLER)*rho*ee;
}

/**********************************************************************
 * NavierStokes2D_pState::D2 -- ee destruction term.                  *
 **********************************************************************/

inline double NavierStokes2D_pState::D2(void) const {
  return Cd5*epsilon()/max(k,TOLER)*rho*ee;
}

/**********************************************************************
 * NavierStokes2D_pState::productionK -- Ke source term.              *
 **********************************************************************/

inline double NavierStokes2D_pState::productionK(const NavierStokes2D_pState &dWdx,
						 const NavierStokes2D_pState &dWdy,
						 const double &ywall,
						 const double &yplus) const {
  double deriv2 = sqr(gm1i/sqr(rho))*(sqr(rho*dWdx.p-p*dWdx.rho) + sqr(rho*dWdy.p-p*dWdy.rho));
  return TWO*rho*alphaT(ywall,yplus)*deriv2;
}

/**********************************************************************
 * NavierStokes2D_pState::productionE1 -- ee source term 1.           *
 **********************************************************************/

inline double NavierStokes2D_pState::productionE1(const NavierStokes2D_pState &dWdx,
						  const NavierStokes2D_pState &dWdy,
						  const double &ywall,
						  const double &yplus) const {
  double deriv2 = sqr(gm1i/sqr(rho))*(sqr(rho*dWdx.p-p*dWdx.rho) + sqr(rho*dWdy.p-p*dWdy.rho));
  return rho*alphaT(ywall,yplus)*(Cd1*ee/max(ke,TOLER)+Cd2*epsilon()/max(k,TOLER))*deriv2 ;
}

/**********************************************************************
 * NavierStokes2D_pState::productionE1_1 -- ee source term 1_1.       *
 **********************************************************************/

inline double NavierStokes2D_pState::productionE1_1(const NavierStokes2D_pState &dWdx,
						    const NavierStokes2D_pState &dWdy,
						    const double &ywall,
						    const double &yplus) const {
  double deriv2 = sqr(gm1i/sqr(rho))*(sqr(rho*dWdx.p-p*dWdx.rho) + sqr(rho*dWdy.p-p*dWdy.rho));  
  return rho*alphaT(ywall,yplus)*(Cd1*ee/max(ke,TOLER))*deriv2 ;  
}

/**********************************************************************
 * NavierStokes2D_pState::productionE1_2 -- ee source term 1_2.       *
 **********************************************************************/

inline double NavierStokes2D_pState::productionE1_2(const NavierStokes2D_pState &dWdx,
						    const NavierStokes2D_pState &dWdy,
						    const double &ywall,
						    const double &yplus) const {
  double deriv2 = sqr(gm1i/sqr(rho))*(sqr(rho*dWdx.p-p*dWdx.rho) + sqr(rho*dWdy.p-p*dWdy.rho));  
  return rho*alphaT(ywall,yplus)*(Cd2*epsilon()/max(k,TOLER))*deriv2 ;  
}

/**********************************************************************
 * NavierStokes2D_pState::productionE2 -- ee source term 2.           *
 **********************************************************************/

inline double NavierStokes2D_pState::productionE2(const double &production_cap) const {
  //return 0.0;
  return Cd3*production_cap*ee/max(k,TOLER);
}

/**********************************************************************
 * NavierStokes2D_cState::NavierStokes2D_cState -- Constructor.       *
 **********************************************************************/
inline NavierStokes2D_cState::NavierStokes2D_cState(const NavierStokes2D_pState &W) {
  rho = W.rho; dv = W.dv(); E = W.E(); dk = W.dk(); domega = W.domega(); dke = W.dke();dee = W.dee();
}

/**********************************************************************
 * NavierStokes2D_cState::W -- Primitive solution state.              *
 **********************************************************************/
inline NavierStokes2D_pState NavierStokes2D_cState::W(void) const {
  return NavierStokes2D_pState(rho,v().x,v().y,p(),k(),omega(),ke(),ee());
}

inline NavierStokes2D_pState NavierStokes2D_cState::W(const NavierStokes2D_cState &U) const {
  return U.W();
}

inline NavierStokes2D_pState W(const NavierStokes2D_cState &U) {
  return U.W();
}

/**********************************************************************
 * NavierStokes2D_cState::dUdW -- Jacobian of the conserved solution  *
 *                                variables with respect to the       *
 *                                primitive solution variables.       *
 **********************************************************************/
inline void NavierStokes2D_cState::dUdW(DenseMatrix &dUdW) const {
  dUdW(0,0) += ONE;
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
    if (Variable_Prandtl == ON) {
      dUdW(6,0) += ke();
      dUdW(6,6) += rho;
      dUdW(7,0) += ee();
      dUdW(7,7) += rho;
    }
  }
}

/**********************************************************************
 * NavierStokes2D_cState::dWdU -- Jacobian of the primitive solution  *
 *                         variables with respect to the conserved    *
 *                         solution variables.                        *
 **********************************************************************/
inline void NavierStokes2D_cState::dWdU(DenseMatrix &dWdU) const {
  dWdU(0,0) += ONE;
  dWdU(1,0) -= dv.x/(rho*rho);
  dWdU(1,1) += ONE/rho;
  dWdU(2,0) -= dv.y/(rho*rho);
  dWdU(2,2) += ONE/rho;
  dWdU(3,0) += HALF*gm1*dv.sqr()/(rho*rho);
  dWdU(3,1) -= gm1*dv.x/rho;
  dWdU(3,2) -= gm1*dv.y/rho;
  dWdU(3,3) += gm1;
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
    dWdU(3,4) -= gm1;
    dWdU(4,0) -= k()/(rho*rho);
    dWdU(4,4) += ONE/rho;
    dWdU(5,0) -= omega()/(rho*rho);
    dWdU(5,5) += ONE/rho;
    if (Variable_Prandtl == ON) {
      dWdU(6,0) -= ke()/(rho*rho);
      dWdU(6,6) += ONE/rho;
      dWdU(7,0) -= ee()/(rho*rho);
      dWdU(7,7) += ONE/rho;
    }
  }
}

/**********************************************************************
 * NavierStokes2D_cState::F -- Solution inviscid flux (x-direction).  *
 **********************************************************************/
inline NavierStokes2D_cState NavierStokes2D_cState::F(void) const {
    if (Variable_Prandtl == ON) {
      return NavierStokes2D_cState(dv.x,
				   sqr(dv.x)/rho+p()+(2.0/3.0)*dk,
				   dv.x*dv.y/rho,
				   dv.x*H()/rho+v().x*(2.0/3.0)*dk,
				   v().x*dk,
				   v().x*domega,
				   v().x*dke,
				   v().x*dee);
    }else{
      return NavierStokes2D_cState(dv.x,
				   sqr(dv.x)/rho+p()+(2.0/3.0)*dk,
				   dv.x*dv.y/rho,
				   dv.x*H()/rho+v().x*(2.0/3.0)*dk,
				   v().x*dk,
				   v().x*domega,
				   ZERO,ZERO);
    }
}

inline NavierStokes2D_cState NavierStokes2D_cState::F(const Vector2D &V) const {
  double vx = v().x;
  if (Variable_Prandtl == ON) {
    return NavierStokes2D_cState(rho*(vx-V.x),
				 (vx-V.x)*dv.x+p()+(2.0/3.0)*dk,
				 (vx-V.x)*dv.y,
				 (vx-V.x)*E+ vx*(p()+(2.0/3.0)*dk),
				 (vx-V.x)*dk,
				 (vx-V.x)*domega,
				 (vx-V.x)*dke,
				 (vx-V.x)*dee);
  }else{
    return NavierStokes2D_cState(rho*(vx-V.x),
				 (vx-V.x)*dv.x+p()+(2.0/3.0)*dk,
				 (vx-V.x)*dv.y,
				 (vx-V.x)*E+ vx*(p()+(2.0/3.0)*dk),
				 (vx-V.x)*dk,
				 (vx-V.x)*domega,
				 ZERO, ZERO);
  }
}

/**********************************************************************
 * NavierStokes2D_cState::dFdU -- Jacobian of the inviscid solution   *
 *                                flux with respect to the conserved  *
 *                                solution variables.                 *
 **********************************************************************/
inline void NavierStokes2D_cState::dFdU(DenseMatrix &dFdU) const {
  dFdU(0,1) += ONE;
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
    if (Variable_Prandtl == ON) {
      dFdU(6,0) -= v().x*ke();
      dFdU(6,1) += ke();
      dFdU(6,6) += v().x;
      dFdU(7,0) -= v().x*ee();
      dFdU(7,1) += ee();
      dFdU(7,7) += v().x;
    }
  }
}

/**********************************************************************
 * NavierStokes2D_cState::Gx, Gy -- Solution viscous fluxes.          *
 **********************************************************************/
inline NavierStokes2D_cState NavierStokes2D_cState::Gx(const NavierStokes2D_pState &dWdx,
						       const NavierStokes2D_pState &W,
						       const double &ywall,const double &yplus) const{
  if (Variable_Prandtl == ON ){
    return NavierStokes2D_cState(ZERO,
				 tau.xx,
				 tau.xy,
				 -q.x+v().x*tau.xx+v().y*tau.xy+(mu()+sigma_k*muT())*dWdx.k,
				 (mu()+sigma_k*muT())*dWdx.k,
				 (mu()+sigma_omega*muT())*dWdx.omega,
				 rho*(Alpha()+alphaT(W,ywall,yplus)/sigma_k_e)*dWdx.ke,
				 rho*(Alpha()+alphaT(W,ywall,yplus)/sigma_ep_e)*dWdx.ee);
  }else{
    return NavierStokes2D_cState(ZERO,
				 tau.xx,
				 tau.xy,
				 -q.x+v().x*tau.xx+v().y*tau.xy+(mu()+sigma_k*muT())*dWdx.k,
				 (mu()+sigma_k*muT())*dWdx.k,
				 (mu()+sigma_omega*muT())*dWdx.omega,
				 ZERO, ZERO);
  }
}

inline NavierStokes2D_cState NavierStokes2D_cState::Gy(const NavierStokes2D_pState &dWdy,
						       const NavierStokes2D_pState &W,
						       const double &ywall,const double &yplus)const{
  if (Variable_Prandtl == ON){
    return NavierStokes2D_cState(ZERO,
				 tau.xy,
				 tau.yy,
				 -q.y+v().x*tau.xy+v().y*tau.yy+(mu()+sigma_k*muT())*dWdy.k,
				 (mu()+sigma_k*muT())*dWdy.k,
				 (mu()+sigma_omega*muT())*dWdy.omega,
				 rho*(Alpha()+alphaT(W,ywall,yplus)/sigma_k_e)*dWdy.ke,
				 rho*(Alpha()+alphaT(W,ywall,yplus)/sigma_ep_e)*dWdy.ee);
  }else{    
    return NavierStokes2D_cState(ZERO,
				 tau.xy,
				 tau.yy,
				 -q.y+v().x*tau.xy+v().y*tau.yy+(mu()+sigma_k*muT())*dWdy.k,
				 (mu()+sigma_k*muT())*dWdy.k,
				 (mu()+sigma_omega*muT())*dWdy.omega, 
				 ZERO, ZERO);
  }
}


/**********************************************************************
 * NavierStokes2D_cState::ComputeViscousTerms -- Compute viscous      *
 *                                               stress tensor and    *
 *                                               heat flux vector.    *
 **********************************************************************/
inline void NavierStokes2D_cState::ComputeViscousTerms(const NavierStokes2D_pState &dWdx,
						       const NavierStokes2D_pState &dWdy,
						       const NavierStokes2D_pState &W,
						       const Vector2D &X,
						       const int &Axisymmetric,
						       const int &adiabatic_flag,
						       const double &ywall,
						       const double &yplus) {

  double div, radius, mumu, kap;

  // Total (i.e. molecular + turbulent) dynamic viscosity and thermal conductivity
  mumu = mu() + muT();
  kap = kappa() + kappaT(W,ywall,yplus);

  if (Axisymmetric){
    radius = max(X.y,TOLER);
  }

  // Divergence of the velocity field.
  div = dWdx.v.x + dWdy.v.y;
  if (Axisymmetric){
    div += v().y/radius;
  }

  // Stress tensor.
  div /= 3.0;	// divide divergence by three
  tau.xx = 2.0*mumu*(dWdx.v.x - div);
  tau.xy = mumu*(dWdy.v.x + dWdx.v.y);
  tau.yy = 2.0*mumu*(dWdy.v.y - div);
  if (Axisymmetric){
    tau.zz = 2.0*mumu*(v().y/radius - div);
  } else {
    tau.zz = ZERO;
  }

  // Heat flux
  q.x = -kap*(dWdx.p - (p()/rho)*dWdx.rho)/(rho*R);
  q.y = -kap*(dWdy.p - (p()/rho)*dWdy.rho)/(rho*R);
}

/**********************************************************************
 * NavierStokes2D_cState::lambda_x -- Eigenvalue(s) (x-direction).    *
 **********************************************************************/
inline NavierStokes2D_pState NavierStokes2D_cState::lambda_x(void) const {
  double vx = v().x, cc = c();
  return NavierStokes2D_pState(vx-cc,vx,vx,vx+cc,vx,vx,vx,vx);
}

inline NavierStokes2D_pState NavierStokes2D_cState::lambda_x(const Vector2D &V) const {
  double vx = v().x, cc = c();
  return NavierStokes2D_pState(vx-V.x-cc,vx-V.x,vx-V.x,vx-V.x+cc,vx-V.x,vx-V.x,vx-V.x,vx-V.x);
}

/**********************************************************************
 * NavierStokes2D_pState::rp_x -- Primitive right eigenvector         *
 *                                (x-direction).                      *
 **********************************************************************/
inline NavierStokes2D_pState NavierStokes2D_cState::rp_x(int index) const {
  //assert(index >= 1 && index <= NUM_VAR_NAVIERSTOKES2D);
  switch(index) {
  case 1 :
    return NavierStokes2D_pState(ONE,-c()/rho,ZERO,c2()-(2.0/3.0)*k(),ZERO,ZERO,ZERO,ZERO);
    break;
  case 2 :
    return NavierStokes2D_pState(ONE,ZERO,ZERO,-(2.0/3.0)*k(),ZERO,ZERO,ZERO,ZERO);
    break;
  case 3 :
    return NavierStokes2D_pState(ZERO,ZERO,ONE,ZERO,ZERO,ZERO,ZERO,ZERO);
    break;
  case 4 :
    return NavierStokes2D_pState(ONE,c()/rho,ZERO,c2()-(2.0/3.0)*k(),ZERO,ZERO,ZERO,ZERO);
    break;
  case 5 :
    return NavierStokes2D_pState(ZERO,ZERO,ZERO,-(2.0/3.0)*rho,ONE,ZERO,ZERO,ZERO);
    break;
  case 6 :
    return NavierStokes2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO,ZERO);
    break;
  case 7 :
    return NavierStokes2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO);
    break;
  case 8 :
    return NavierStokes2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE);
    break;
  };
  return NavierStokes2D_pState(ONE,-c()/rho,ZERO,c2()-(2.0/3.0)*k(),ZERO,ZERO,ZERO,ZERO);
}

/**********************************************************************
 * NavierStokes2D_cState::rc_x -- Conserved right eigenvector         *
 *                                (x-direction).                      *
 **********************************************************************/
inline NavierStokes2D_cState NavierStokes2D_cState::rc_x(int index) const {
  //assert(index >= 1 && index <= NUM_VAR_NAVIERSTOKES2D);
  switch(index) {
  case 1 :
    return NavierStokes2D_cState(ONE,v().x-c(),v().y,h()-c()*v().x+(2.0/3.0)*k(),k(),omega(),ke(),ee());
    break;
  case 2 :
    return NavierStokes2D_cState(ONE,v().x,v().y,HALF*v().sqr()+gm1i*(g-5.0/3.0)*k(),k(),omega(),ke(),ee());
    break;
  case 3 :
    return NavierStokes2D_cState(ZERO,ZERO,rho,dv.y,ZERO,ZERO,ZERO,ZERO);
    break;
  case 4 :
    return NavierStokes2D_cState(ONE,v().x+c(),v().y,h()+c()*v().x+(2.0/3.0)*k(),k(),omega(),ke(),ee());
    break;
  case 5 :
    return NavierStokes2D_cState(ZERO,ZERO,ZERO,rho*gm1i*(g-5.0/3.0),rho,ZERO,ZERO,ZERO);
    break;
  case 6 :
    return NavierStokes2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,rho,ZERO,ZERO);
    break;
  case 7 :
    return NavierStokes2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,rho,ZERO);
    break;
  case 8 :
    return NavierStokes2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,rho);
    break;
  };
  return NavierStokes2D_cState(ONE,v().x-c(),v().y,h()-c()*v().x+(2.0/3.0)*k(),k(),omega(),ke(),ee());
}

/**********************************************************************
 * NavierStokes2D_cState::lp_x -- Primitive left eigenvector          *
 *                                (x-direction).                      *
 **********************************************************************/
inline NavierStokes2D_pState NavierStokes2D_cState::lp_x(int index) const {
  //assert(index >= 1 && index <= NUM_VAR_NAVIERSTOKES2D);
  switch(index) {
  case 1 :
    return NavierStokes2D_pState(k()/(3.0*c2()),-HALF*rho/c(),ZERO,HALF/c2(),rho/(3.0*c2()),ZERO,ZERO,ZERO);
    break;
  case 2 :
    return NavierStokes2D_pState(ONE-(2.0/3.0)*k()/c2(),ZERO,ZERO,-ONE/c2(),-(2.0/3.0)*rho/c2(),ZERO,ZERO,ZERO);
    break;
  case 3 :
    return NavierStokes2D_pState(ZERO,ZERO,ONE,ZERO,ZERO,ZERO,ZERO,ZERO);
    break;
  case 4 :
    return NavierStokes2D_pState(k()/(3.0*c2()),HALF*rho/c(),ZERO,HALF/c2(),rho/(3.0*c2()),ZERO,ZERO,ZERO);
    break;
  case 5 :
    return NavierStokes2D_pState(ZERO,ZERO,ZERO,ZERO,ONE,ZERO,ZERO,ZERO);
    break;
  case 6 :
    return NavierStokes2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO,ZERO);
    break;
  case 7 :
    return NavierStokes2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE,ZERO);
    break;
  case 8 :
    return NavierStokes2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,ONE);
    break;
  };
  return NavierStokes2D_pState(k()/(3.0*c2()),-HALF*rho/c(),ZERO,HALF/c2(),rho/(3.0*c2()),ZERO,ZERO,ZERO);
}

/**********************************************************************
 * NavierStokes2DState -- External subroutines.                       *
 **********************************************************************/

extern NavierStokes2D_pState Riemann(const NavierStokes2D_pState &Wl,
				     const NavierStokes2D_pState &Wr);

extern NavierStokes2D_pState Riemann_Wrs(const NavierStokes2D_pState &Wl,
					 const NavierStokes2D_pState &Wr,
					 NavierStokes2D_pState &Wls,
					 NavierStokes2D_pState &Wrs);

extern NavierStokes2D_pState RoeAverage(const NavierStokes2D_pState &Wl,
					const NavierStokes2D_pState &Wr);

extern NavierStokes2D_pState Rotate(const NavierStokes2D_pState &W,
				    const Vector2D &norm_dir);

extern NavierStokes2D_cState Rotate(const NavierStokes2D_cState &U,
				    const Vector2D &norm_dir);

extern NavierStokes2D_pState Translate(const NavierStokes2D_pState &W,
				       const Vector2D &V);

extern NavierStokes2D_pState Reflect(const NavierStokes2D_pState &W,
				     const Vector2D &norm_dir);

extern NavierStokes2D_pState Mirror(const NavierStokes2D_pState &W,
				    const Vector2D &norm_dir);

extern NavierStokes2D_pState WallViscousHeatFlux(const NavierStokes2D_pState &W,
						 const Vector2D &norm_dir);

extern NavierStokes2D_pState WallViscousIsothermal(const NavierStokes2D_pState &W,
						   const Vector2D &norm_dir,
						   const double &Twall);

extern NavierStokes2D_pState MovingWallHeatFlux(const NavierStokes2D_pState &W,
						const Vector2D &norm_dir,
						const double &Vwall);

extern NavierStokes2D_pState MovingWallIsothermal(const NavierStokes2D_pState &W,
						  const Vector2D &norm_dir,
						  const double &Vwall,
						  const double &Twall);

extern NavierStokes2D_pState Reflect(const NavierStokes2D_pState &W,
				     const Vector2D &norm_dir,
				     const Vector2D &V);

extern NavierStokes2D_pState WallViscousHeatFlux(const NavierStokes2D_pState &W,
						 const Vector2D &norm_dir,
						 const Vector2D &V);

extern NavierStokes2D_pState WallViscousIsothermal(const NavierStokes2D_pState &W,
						   const Vector2D &norm_dir,
						   const Vector2D &V,
						   const double &Twall);

extern NavierStokes2D_pState BurningSurface(const NavierStokes2D_pState &W,
					    const Vector2D &norm_dir);

extern NavierStokes2D_pState MassInjection2(const NavierStokes2D_pState &W,
					    const Vector2D &X,
					    const Vector2D &norm_dir,
					    const double &Twall);

extern NavierStokes2D_pState MassInjection(const NavierStokes2D_pState &W,
					   const Vector2D &norm_dir);

extern NavierStokes2D_pState RinglebFlow(const NavierStokes2D_pState &Wdum,
					 const Vector2D &X);

extern NavierStokes2D_pState RinglebFlow(const NavierStokes2D_pState &Wdum,
					 const Vector2D &X,
					 double &q, double &k);

extern NavierStokes2D_pState RinglebFlowAverageState(const NavierStokes2D_pState &Wdum,
						     const Vector2D &Y1,
						     const Vector2D &Y2,
						     const Vector2D &Y3,
						     const Vector2D &Y4);

extern NavierStokes2D_pState ViscousChannelFlow(const NavierStokes2D_pState &Wdum,
						const Vector2D X,
						const Vector2D Vwall,
						const double dp,
						const double length,
						const double height);

extern NavierStokes2D_pState ViscousChannelFlowVelocity(const NavierStokes2D_pState &Wdum,
							const Vector2D X,
							const Vector2D Vwall,
							const double dp,
							const double length,
							const double height);

extern NavierStokes2D_pState ViscousPipeFlow(const NavierStokes2D_pState &Wdum,
					     const Vector2D X,
					     const double dp,
					     const double length,
					     const double radius);

extern NavierStokes2D_pState TurbulentPipeFlow(const NavierStokes2D_pState &Wo,
					       const Vector2D X,
					       const double dp,
					       const double length,
					       const double radius,
					       const double ReN);

extern NavierStokes2D_pState FlatPlate(const NavierStokes2D_pState &Winf,
				       const Vector2D &X,
				       const double &plate_length,
				       double &eta,
				       double &f,
				       double &fp,
				       double &fpp);

extern NavierStokes2D_pState DrivenCavityFlow(const NavierStokes2D_pState &Wo,
					      const double &l,
					      const double &Re);

extern NavierStokes2D_pState BackwardFacingStep(const NavierStokes2D_pState &Wo,
						const Vector2D &X,
						const double &h,
						const double &ho,
						const double &Re,
						const double &M);

extern NavierStokes2D_pState BC_Characteristic(const NavierStokes2D_pState &Wi,
					       const NavierStokes2D_pState &Wo,
					       const Vector2D &norm_dir);

extern NavierStokes2D_pState BC_Characteristic_Pressure(const NavierStokes2D_pState &Wi,
							const NavierStokes2D_pState &Wo,
							const Vector2D &norm_dir);

extern NavierStokes2D_pState BC_Characteristic_Mach_Number(const NavierStokes2D_pState &Wi,
							   const NavierStokes2D_pState &Wo,
							   const Vector2D &norm_dir);

Vector2D HLLE_wavespeeds(const NavierStokes2D_pState &Wl,
    const NavierStokes2D_pState &Wr,
    const Vector2D &norm_dir);

extern NavierStokes2D_pState WaveSpeedPos(const NavierStokes2D_pState &lambda_a,
					  const NavierStokes2D_pState &lambda_l,
					  const NavierStokes2D_pState &lambda_r);

extern NavierStokes2D_pState WaveSpeedNeg(const NavierStokes2D_pState &lambda_a,
					  const NavierStokes2D_pState &lambda_l,
					  const NavierStokes2D_pState &lambda_r);

extern NavierStokes2D_pState WaveSpeedAbs(const NavierStokes2D_pState &lambda_a,
					  const NavierStokes2D_pState &lambda_l,
					  const NavierStokes2D_pState &lambda_r);

extern NavierStokes2D_pState HartenFixPos(const NavierStokes2D_pState &lambda_a,
					  const NavierStokes2D_pState &lambda_l,
					  const NavierStokes2D_pState &lambda_r);

extern NavierStokes2D_pState HartenFixNeg(const NavierStokes2D_pState &lambda_a,
					  const NavierStokes2D_pState &lambda_l,
					  const NavierStokes2D_pState &lambda_r);

extern NavierStokes2D_pState HartenFixAbs(const NavierStokes2D_pState &lambda_a,
					  const NavierStokes2D_pState &lambda_l,
					  const NavierStokes2D_pState &lambda_r);

extern NavierStokes2D_cState FluxGodunov_n(const NavierStokes2D_pState &Wl,
					   const NavierStokes2D_pState &Wr,
					   const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxGodunov_n(const NavierStokes2D_cState &Ul,
					   const NavierStokes2D_cState &Ur,
					   const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxGodunov_Wrs_n(const NavierStokes2D_pState &Wl,
					       const NavierStokes2D_pState &Wr,
					       const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxGodunov_Wrs_n(const NavierStokes2D_cState &Ul,
					       const NavierStokes2D_cState &Ur,
					       const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxGodunov_MB_n(const NavierStokes2D_pState &Wl,
					      const NavierStokes2D_pState &Wr,
					      const Vector2D &V,
					      const Vector2D &norm_dir);

extern NavierStokes2D_pState StateGodunov_n(const NavierStokes2D_pState &Wl,
					    const NavierStokes2D_pState &Wr,
					    const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxRoe(const NavierStokes2D_pState &Wl,
				     const NavierStokes2D_pState &Wr);

extern NavierStokes2D_cState FluxRoe(const NavierStokes2D_cState &Ul,
				     const NavierStokes2D_cState &Ur);

extern NavierStokes2D_cState FluxRoe_n(const NavierStokes2D_pState &Wl,
				       const NavierStokes2D_pState &Wr,
				       const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxRoe_n(const NavierStokes2D_cState &Ul,
				       const NavierStokes2D_cState &Ur,
				       const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxRoe_MB(const NavierStokes2D_pState &Wl,
					const NavierStokes2D_pState &Wr,
					const Vector2D &V);

extern NavierStokes2D_cState FluxRoe_MB(const NavierStokes2D_cState &Ul,
					const NavierStokes2D_cState &Ur,
					const Vector2D &V);

extern NavierStokes2D_cState FluxRoe_MB_n(const NavierStokes2D_pState &Wl,
					  const NavierStokes2D_pState &Wr,
					  const Vector2D &V,
					  const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxRoe_MB_n(const NavierStokes2D_cState &Ul,
					  const NavierStokes2D_cState &Ur,
					  const Vector2D &V,
					  const Vector2D &norm_dir);

extern NavierStokes2D_cState StateRoe(const NavierStokes2D_pState &Wl,
				      const NavierStokes2D_pState &Wr);

extern NavierStokes2D_cState StateRoe_n(const NavierStokes2D_pState &Wl,
					const NavierStokes2D_pState &Wr,
					const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxRusanov(const NavierStokes2D_pState &Wl,
					 const NavierStokes2D_pState &Wr);

extern NavierStokes2D_cState FluxRusanov(const NavierStokes2D_cState &Ul,
					 const NavierStokes2D_cState &Ur);

extern NavierStokes2D_cState FluxRusanov_n(const NavierStokes2D_pState &Wl,
					   const NavierStokes2D_pState &Wr,
					   const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxRusanov_n(const NavierStokes2D_cState &Ul,
					   const NavierStokes2D_cState &Ur,
					   const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxHLLE(const NavierStokes2D_pState &Wl,
				      const NavierStokes2D_pState &Wr);

extern NavierStokes2D_cState FluxHLLE(const NavierStokes2D_cState &Ul,
				      const NavierStokes2D_cState &Ur);

extern NavierStokes2D_cState FluxHLLE_n(const NavierStokes2D_pState &Wl,
					const NavierStokes2D_pState &Wr,
					const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxHLLE_n(const NavierStokes2D_cState &Ul,
					const NavierStokes2D_cState &Ur,
					const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxHLLE_MB(const NavierStokes2D_pState &Wl,
					 const NavierStokes2D_pState &Wr,
					 const Vector2D &V);

extern NavierStokes2D_cState FluxHLLE_MB(const NavierStokes2D_cState &Ul,
					 const NavierStokes2D_cState &Ur,
					 const Vector2D &V);

extern NavierStokes2D_cState FluxHLLE_MB_n(const NavierStokes2D_pState &Wl,
					   const NavierStokes2D_pState &Wr,
					   const Vector2D &V,
					   const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxHLLE_MB_n(const NavierStokes2D_cState &Ul,
					   const NavierStokes2D_cState &Ur,
					   const Vector2D &V,
					   const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxHLLL(const NavierStokes2D_pState &Wl,
				      const NavierStokes2D_pState &Wr);

extern NavierStokes2D_cState FluxHLLL(const NavierStokes2D_cState &Ul,
				      const NavierStokes2D_cState &Ur);

extern NavierStokes2D_cState FluxHLLL_n(const NavierStokes2D_pState &Wl,
					const NavierStokes2D_pState &Wr,
					const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxHLLL_n(const NavierStokes2D_cState &Ul,
					const NavierStokes2D_cState &Ur,
					const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxHLLC(const NavierStokes2D_pState &Wl,
				      const NavierStokes2D_pState &Wr);

extern NavierStokes2D_cState FluxHLLC(const NavierStokes2D_cState &Ul,
				      const NavierStokes2D_cState &Ur);

extern NavierStokes2D_cState FluxHLLC_n(const NavierStokes2D_pState &Wl,
					const NavierStokes2D_pState &Wr,
					const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxHLLC_n(const NavierStokes2D_cState &Ul,
					const NavierStokes2D_cState &Ur,
					const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxVanLeer(const NavierStokes2D_pState &Wl,
					 const NavierStokes2D_pState &Wr);

extern NavierStokes2D_cState FluxVanLeer(const NavierStokes2D_cState &Wl,
					 const NavierStokes2D_cState &Wr);

extern NavierStokes2D_cState FluxVanLeer_n(const NavierStokes2D_pState &Wl,
					   const NavierStokes2D_pState &Wr,
					   const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxVanLeer_n(const NavierStokes2D_cState &Wl,
					   const NavierStokes2D_cState &Wr,
					   const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxVanLeer_MB(const NavierStokes2D_pState &Wl,
					    const NavierStokes2D_pState &Wr,
					    const Vector2D &V);

extern NavierStokes2D_cState FluxVanLeer_MB_n(const NavierStokes2D_pState &Wl,
					      const NavierStokes2D_pState &Wr,
					      const Vector2D &V,
					      const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxAUSM(const NavierStokes2D_pState &Wl,
				      const NavierStokes2D_pState &Wr);

extern NavierStokes2D_cState FluxAUSM(const NavierStokes2D_cState &Wl,
				      const NavierStokes2D_cState &Wr);

extern NavierStokes2D_cState FluxAUSM_n(const NavierStokes2D_pState &Wl,
					const NavierStokes2D_pState &Wr,
					const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxAUSM_n(const NavierStokes2D_cState &Wl,
					const NavierStokes2D_cState &Wr,
					const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxAUSMplus(const NavierStokes2D_pState &Wl,
					  const NavierStokes2D_pState &Wr);

extern NavierStokes2D_cState FluxAUSMplus(const NavierStokes2D_cState &Wl,
					  const NavierStokes2D_cState &Wr);

extern NavierStokes2D_cState FluxAUSMplus_n(const NavierStokes2D_pState &Wl,
					    const NavierStokes2D_pState &Wr,
					    const Vector2D &norm_dir);

extern NavierStokes2D_cState FluxAUSMplus_n(const NavierStokes2D_cState &Wl,
					    const NavierStokes2D_cState &Wr,
					    const Vector2D &norm_dir);

extern NavierStokes2D_cState ViscousFlux_n(const Vector2D &X,
					   NavierStokes2D_pState &W,
					   const NavierStokes2D_pState &dWdx,
					   const NavierStokes2D_pState &dWdy,
					   const Vector2D &norm_dir,
					   const int &Axisymmetric,
					   const int &adiabatic_flag,
					   const double &ywall,
				           const double &yplus);

extern NavierStokes2D_cState ViscousFluxDiamondPath_n(const Vector2D &X,
						      const Vector2D &Xl, const NavierStokes2D_pState &Wl,
						      const Vector2D &Xd, const NavierStokes2D_pState &Wd,
						      const Vector2D &Xr, const NavierStokes2D_pState &Wr,
						      const Vector2D &Xu, const NavierStokes2D_pState &Wu,
						      const Vector2D &norm_dir,
						      const int &Axisymmetric,
						      const int &stencil_flag,
						      const double &ywall,
						      const double &yplus);

extern NavierStokes2D_cState ViscousFluxHybrid_n(const Vector2D &X,
						 NavierStokes2D_pState &W,
						 const Vector2D &X1,
						 const NavierStokes2D_pState &W1,
						 const NavierStokes2D_pState &dW1dx,
						 const NavierStokes2D_pState &dW1dy,
						 const Vector2D &X2,
						 const NavierStokes2D_pState &W2,
						 const NavierStokes2D_pState &dW2dx,
						 const NavierStokes2D_pState &dW2dy,
						 const Vector2D &norm_dir,
						 const int &Axisymmetric,
						 const double &ywall,
			                         const double &yplus);

extern NavierStokes2D_cState ViscousFluxDiamondPath_n(const Vector2D &X,
						      const Vector2D &Xl, const NavierStokes2D_pState &Wl,
						      const Vector2D &Xd, const NavierStokes2D_pState &Wd,
						      const Vector2D &Xr, const NavierStokes2D_pState &Wr,
						      const Vector2D &Xu, const NavierStokes2D_pState &Wu,
						      const Vector2D &norm_dir,
						      const int &Axisymmetric,
						      const int &stencil_flag,
						      const double &ywall,
						      const double &yplus,
						      NavierStokes2D_pState &dWdx, NavierStokes2D_pState &dWdy);

extern NavierStokes2D_cState ViscousFluxHybrid_n(const Vector2D &X,
						 NavierStokes2D_pState &W,
						 const Vector2D &X1,
						 const NavierStokes2D_pState &W1,
						 const NavierStokes2D_pState &dW1dx,
						 const NavierStokes2D_pState &dW1dy,
						 const Vector2D &X2,
						 const NavierStokes2D_pState &W2,
						 const NavierStokes2D_pState &dW2dx,
						 const NavierStokes2D_pState &dW2dy,
						 const Vector2D &norm_dir,
						 const int &Axisymmetric,
						 const double &ywall,
			                         const double &yplus,
						 NavierStokes2D_pState &dWdx, NavierStokes2D_pState &dWdy);

extern double ShearStress(const NavierStokes2D_pState &W,
			  const NavierStokes2D_pState &dWdx,
			  const NavierStokes2D_pState &dWdy,
			  const Vector2D &nhat);

extern double WallShearStress(const NavierStokes2D_pState &W1,
			      const Vector2D &X1,
			      const Vector2D &X2,
			      const Vector2D &X3,
			      const Vector2D &nhat);

extern double WallShearStress2(const Vector2D &X,
			       const Vector2D &X1,
			       const NavierStokes2D_pState &W1,
			       const NavierStokes2D_pState &dW1dx,
			       const NavierStokes2D_pState &dW1dy,
			       const Vector2D &nhat);

#endif // _NAVIERSTOKES2D_STATE_INCLUDED
