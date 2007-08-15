/**********************************************************************
 * HighTemp2DState.h: Header file defining 2D HighTemperature         *
 *                        solution state classes.                     *
 **********************************************************************/

#ifndef _HIGHTEMP2D_STATE_INCLUDED
#define _HIGHTEMP2D_STATE_INCLUDED

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

#include "../Math/Math.h"
#include "../Math/Matrix.h"
#include "../EquationOfState/EquationOfState.h"
#include "../CFD/CFD.h"
#include "../Math/Vector2D.h"
#include "../Math/Tensor2D.h"
#include "../Physics/GasConstants.h"
#include "../Physics/SolidConstants.h"

#define	NUM_VAR_HIGHTEMP2D  6

/* Define some functions for AUSMplusUP flux calc. */
// M+1
inline double Mplus_1(double M)
  { return 0.5*(M + fabs(M)); }

// M-1
inline double Mminus_1(double M)
  { return 0.5*(M - fabs(M)); }

// M+2
inline double Mplus_2(double M)
  { return 0.25*sqr(M + 1.0); }

// M-2
inline double Mminus_2(double M)
  { return -0.25*sqr(M - 1.0); }


// Define the classes.

class HighTemp2D_cState;

/*!
 * Class: HighTemp2D_pState
 *
 * @brief Primitive variable solution state class definition for a 
 *        laminar or turbulent compressible gas-flow.
 *
 * Primitive variable solution state class definition for a laminar or
 * turbulent (k-omega) gas-flow.
 *
 * \verbatim
 * Member functions
 *     air_ideal -- Ideal Equation of State variable
 *     air_tgas -- High Temperature Equation of State variable
 *     rho      -- Gas density.
 *     v        -- Gas velocity.
 *     p        -- Gas pressure.
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
 *     cv       -- Specific heat at constant volume.
 *     T        -- Gas temperature.
 *     e        -- Gas specific internal energy.
 *     E        -- Gas total energy.
 *     h        -- Gas specific enthalpy.
 *     H        -- Gas total enthalpy.
 *     a        -- Gas sound speed.
 *     a2       -- Gas sound speed squared.
 *     M        -- Gas Mach number.
 *     Mref     -- Reference Mach number
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
 *     meanfreepath -- Gas mean free path. - REMOVED
 *     gamma    -- High temperature specific heat ratio
 * 
 *     dpde     -- Partial derivative of pressure wrt internal energy     
 *     dpdrho   -- Partial derivative of pressure wrt density
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
 *
 *     U        -- Return the conserved solution state.
 *     dUdW     -- Return the Jacobian of the conserved solution
 *                 variables with respect to the primitive solution
 *                 variables.
 *     dWdU     -- Return the Jacobian of the primitive solution
 *                 variables with respect to the conserved solution
 *                 variables.
 *
 *     Fx        -- Return x-direction inviscid solution flux.
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
class HighTemp2D_pState {
 private:
 public:
  //@{ @name Primitive variables and associated constants:
  double                       rho; //!< Gas density.
  Vector2D                       v; //!< Gas velocity (2D vector).
  double                         p; //!< Gas pressure.
  double                         k; //!< Gas turbulent kinetic energy.
  double                     omega; //!< Gas specific dissipation rate.
  Tensor2D                     tau; //!< Viscous stess tensor (laminar and turbulent).
  Vector2D                       q; //!< Heat flux vector (laminar and turbulent).
  static double                  g; //!< Specific heat ratio.
  static double                gm1; //!< g-1
  static double               gm1i; //!< 1/(g-1)
  static double                  R; //!< Gas constant.
  static double                 cp; //!< Specific heat at constant pressure.
  static double                 cv; //!< Specific heat at constant volume.
  static double     v1,v2,v3,v4,v5; //!< Viscosity law coefficients.
  static double               Mref; //!< Reference Mach Number
  static int             flow_type; //!< Flow-type indicator (inviscid, laminar, or k-omega turbulent).
  static int              eos_type; //!< Equation of State-type indicator (ideal, tgas(high temperature)).
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
  static double           beta_k_o;
  static double       beta_omega_o;
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
  //@}

  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  HighTemp2D_pState(void) {
    Standard_Atmosphere();
  }

  //! Copy constructor.
  HighTemp2D_pState(const HighTemp2D_pState &W) {
    rho = W.rho; v.x = W.v.x; v.y = W.v.y; p = W.p; k = W.k; omega = W.omega;
  }

  //! Copy constructor.
  HighTemp2D_pState(const HighTemp2D_cState &U);

  //! Assignment constructor.
  HighTemp2D_pState(const double &dens,
		    const Vector2D &V,
		    const double &pre) {
    rho = dens; v.x = V.x; v.y = V.y; p = pre; k = ZERO; omega = ZERO;
  }
  
  //! Assignment constructor.
  HighTemp2D_pState(const double &dens,
			const double &vx,
			const double &vy,
			const double &pre) {
    rho = dens; v.x = vx; v.y = vy; p = pre; k = ZERO; omega = ZERO;
  }

  //! Assignment constructor.
  HighTemp2D_pState(const double &dens,
			const Vector2D &V,
			const double &pre,
			const double &kk,
			const double &omga) {
    rho = dens; v.x = V.x; v.y = V.y; p = pre; k = kk; omega = omga;
  }

  //! Assignment constructor.
  HighTemp2D_pState(const double &dens,
		    const double &vx,
		    const double &vy,
		    const double &pre,
		    const double &kk,
		    const double &omga) {
    rho = dens; v.x = vx; v.y = vy; p = pre; k = kk; omega = omga;
  }

  //! Destructor.
  ~HighTemp2D_pState(void) { }
  //@}

  //! Return the number of variables.
  int NumVar(void) { return NUM_VAR_HIGHTEMP2D; }

  //@{ @name Set static variables.
  static void set_static_variables(void);
  static void set_static_variables(char *gas_type,
			    const int &EOSType,
			    const int &FlowType,
			    const double &C_constant,
			    const double &von_karman_constant,
			    const double &yplus_sub,
			    const double &yplus_buffer,
			    const double &yplus_outer,
			    char *propellant_type);
  static void set_gas(char *gas_type);
  static void set_turbulence(const double &C_constant,
		      const double &von_karman,
		      const double &yplus_sub,
		      const double &yplus_buffer,
		      const double &yplus_outer);
  static void set_propellant(char *propellant_type);
  //@}

  //@{ @name Useful operators.
  //! Copy operator.
  void Copy(const HighTemp2D_pState &W);

  //! Vacuum operator.
  void Vacuum(void);

  //! Standard atmosphere operator.
  void Standard_Atmosphere(void);

  //! Check for unphysical state properties.
  int Unphysical_Properties(void) const;
  //@}

  //@{ @name Gas related functions.
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

  //! Specific Heat ratio for high-temperature
  // double gamma(void) const;
	//  gamma() is not used so I removed it. Alistair Wood. Wed Aug 01 2007.

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

  //! Gas mean free path.
  //double meanfreepath(void) const;
  //@}

  //! High Temperature dp/de cell-centred calculation.
  double dpde(void) const;

  //! High Temperature dp/drho cell-centred calculation.
  double dpdrho(void) const;

	double dhdrho(void) const;
	double dhdp(void) const;

  //! High Temperature dp/drho cell-centred calculation.
  double dTdrho(void) const;

  //! High Temperature dp/drho cell-centred calculation.
  double dTdp(void) const;

  //! High Temperature dp/drho cell-centred calculation.
  double ddTdrho(void) const;

  //! High Temperature dp/drho cell-centred calculation.
  double ddTdp(void) const;

  //! High Temperature dp/drho cell-centred calculation.
  double ddTdpdrho(void) const;

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

  double beta_k(const HighTemp2D_pState &dWdx,const HighTemp2D_pState &dWdy) const;
  double beta_omega(const HighTemp2D_pState &dWdx, const HighTemp2D_pState &dWdy) const;
  double f_beta_k(const HighTemp2D_pState &dWdx, const HighTemp2D_pState &dWdy) const;
  double f_beta_omega(const HighTemp2D_pState &dWdx, const HighTemp2D_pState &dWdy) const;
  double chi_k(const HighTemp2D_pState &dWdx, const HighTemp2D_pState &dWdy) const;
  double chi_omega(const HighTemp2D_pState &dWdx, const HighTemp2D_pState &dWdy) const;
  //@}

  //@{ @name Solid propellant related functions.
  //! Burning rate.
  double burningrate(void) const;
  //@}

  //@{ @name Conserved solution state.
  HighTemp2D_cState U(void) const;
  HighTemp2D_cState U(const HighTemp2D_pState &W) const;
  friend HighTemp2D_cState U(const HighTemp2D_cState &W);
  //@}

  //! Jacobian of the conserved solution variables with respect to the
  //! primitive solution variables.
  void dUdW(DenseMatrix &dUdW) const;

  //! Jacobian of the primitive solution variables with respect to the
  //! conserved solution variables.
  void dWdU(DenseMatrix &dWdU) const;

  //@{ @name Inviscid solution flux (x-direction) and Jacobian.
  //  HighTemp2D_cState F(void) const;
  // HighTemp2D_cState F(const Vector2D &V) const;
  HighTemp2D_cState Fx(void) const;
  HighTemp2D_cState Fx(const Vector2D &V) const;
  //HighTemp2D_cState Fx(double dpde, double dpdrho) const;
  void dFdU(DenseMatrix &dFdU) const;
  void dFdW(DenseMatrix &dFdW) const;
  //@}

  //@{ @name Viscous solution fluxes and Jacobians.
  HighTemp2D_cState Gx(const HighTemp2D_pState &dWdx) const;
  HighTemp2D_cState Gy(const HighTemp2D_pState &dWdy) const;
  //HighTemp2D_cState dGxdU(???) const;
  //HighTemp2D_cState dGydU(???) const;
  //@}

  //! Compute viscous stress tensor and heat flux vector.
  void ComputeViscousTerms(const HighTemp2D_pState &dWdx,
			   const HighTemp2D_pState &dWdy,
			   const Vector2D &X,
			   int Axisymmetric);

  //@{ @name Eigenstructure.
  //! Eigenvalue(s) (x-direction).
  HighTemp2D_pState lambda_x(void) const;

  //! Eigenvalue(s) (x-direction).
  HighTemp2D_pState lambda_x(double aAvg) const;

  //! Eigenvalue(s) (x-direction).
  HighTemp2D_pState lambda_x(const Vector2D &V) const;

  //! Primitive right eigenvector (x-direction).
  HighTemp2D_pState rp_x(int index) const;

  //! Conserved right eigenvector (x-direction).
  HighTemp2D_cState rc_x(int index) const;
  HighTemp2D_cState rc_x(int index, double dpde, double dpdrho,double cAvg) const;

  //! Primitive left eigenvector (x-direction).
  HighTemp2D_pState lp_x(int index) const;
  HighTemp2D_pState lp_x(int index, double cAvg) const;
  //@}

  //@{ @name Include all source vectors and Jacobians.
  HighTemp2D_cState S(const Vector2D &X,
			  const HighTemp2D_pState &dWdx,
			  const HighTemp2D_pState &dWdy,
			  const int &Axisymmetric) const;
  void dSdU(DenseMatrix &dSdU,
	    const Vector2D &X,
	    const HighTemp2D_pState &dWdx,
	    const HighTemp2D_pState &dWdy,
	    const int &Axisymmetric) const;
  //@}

  //@{ @name Inviscid axisymmetric flow source vector and Jacobian.
  HighTemp2D_cState Si(const Vector2D &X) const;
  void dSidU(DenseMatrix &dSidU, const Vector2D &X) const;
  //@}

  //@{ @name Viscous axisymmetric flow source vector and Jacobian.
  HighTemp2D_cState Sv(const Vector2D &X,
			   const HighTemp2D_pState &dWdy) const;
  void dSvdU(DenseMatrix &dSvdU, const Vector2D &X, const HighTemp2D_pState &dWdy) const;
  //@}

  //@{ @name Turbulent source term vector and Jacobian.
  HighTemp2D_cState St(const Vector2D &X,
			   const HighTemp2D_pState &dWdx,
			   const HighTemp2D_pState &dWdy,
			   const int &Axisymmetric) const;
  void dStdU(DenseMatrix &dStdU,
	     const Vector2D &X,
	     const HighTemp2D_pState &dWdx,
	     const HighTemp2D_pState &dWdy,
	     const int &Axisymmetric) const;
  //@}

  //@{ @name Index operator.
  double &operator[](int index) {
    //assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
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
    // Default return, this is never reached.
    return rho;
  }
  
  const double &operator[](int index) const {
    //assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
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
    // Default return, this is never reached.
    return rho;
  }
  //@}

  //@{ @name Binary arithmetic operators.
  HighTemp2D_pState operator +(const HighTemp2D_pState &W) const;
  HighTemp2D_pState operator -(const HighTemp2D_pState &W) const;
  double operator *(const HighTemp2D_pState &W) const;
  HighTemp2D_pState operator *(const double &a) const;
  friend HighTemp2D_pState operator *(const double &a, const HighTemp2D_pState &W);
  HighTemp2D_pState operator /(const double &a) const;
  HighTemp2D_pState operator ^(const HighTemp2D_pState &W) const;
  //@}

  //@{ @name Assignment operator.
  HighTemp2D_pState& operator =(const HighTemp2D_pState &W);
  //@}

  //@{ @name Unary arithmetic operators.
  //HighTemp2D_pState operator +(const HighTemp2D_pState &W);
  friend HighTemp2D_pState operator -(const HighTemp2D_pState &W);
  //@}

  //@{ @name Shortcut arithmetic operators.
  HighTemp2D_pState &operator +=(const HighTemp2D_pState &W);
  HighTemp2D_pState &operator -=(const HighTemp2D_pState &W);
  HighTemp2D_pState &operator *=(const double &a);
  HighTemp2D_pState &operator /=(const double &a);
  //@}

  //@{ @name Relational operators.
  friend int operator ==(const HighTemp2D_pState &W1, const HighTemp2D_pState &W2);
  friend int operator !=(const HighTemp2D_pState &W1, const HighTemp2D_pState &W2);
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const HighTemp2D_pState &W);
  friend istream &operator >> (istream &in_file, HighTemp2D_pState &W);
  //@}

	/**********************************************************************
	 * Routine: Rotate                                                    *
	 *                                                                    *
	 * This function returns the solution in the local rotated frame.     *
	 *                                                                    *
	 **********************************************************************/
	void Rotate(const HighTemp2D_pState &W, const Vector2D &norm_dir);

  void output_labels(ostream &out_file) {
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
    }
  }

  void output_data(ostream &out_file) {
    out_file << " " << rho << " " << v.x << " " << v.y << " " << p
	     << " " << T() << " " << M() << " " << H() << " " << s();
    if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      out_file << " " << k << " " << omega << " " << epsilon()
	       << " " << ell() << " " << pmodified();
    }
  }

};

/*!
 * Class: HighTemp2D_cState
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
 *     meanfreepath -- Gas mean free path. - REMOVED
 *     gamma    -- High temperature specific heat ratio
 *
 *     dpde     -- Partial derivative of pressure wrt internal energy     
 *     dpdrho   -- Partial derivative of pressure wrt density
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
class HighTemp2D_cState {
  private:
  public:
  //@{ @name Gas conservative variables:
  double                   rho; //!< Gas density.
  Vector2D                  dv; //!< Gas momentum.
  double                     E; //!< Gas total energy.
  double                    dk; //!< Gas total turbulent kinetic energy.
  double                domega; //!< Gas total turbulent specific dissipation rate.
  Tensor2D                 tau; //!< Viscous stess tensor (laminar and turbulent).
  Vector2D                   q; //!< Heat flux vector (laminar and turbulent).
  static double              g; //!< Specific heat ratio.
  static double            gm1; //!< g-1
  static double           gm1i; //!< 1/(g-1)
  static double              R; //!< Gas constant.
  static double             cp; //!< Specific heat at constant pressure.
  static double             cv; //!< Specific heat at constant volume.
  static double v1,v2,v3,v4,v5; //!< Viscosity law coefficients.
  static int         flow_type; //!< Flow-type indicator (inviscid, laminar, or k-omega turbulent).
  static int             eos_type; //!< Eqyation of State-type indicator (ideal, tgas(high temperature)).
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
  static double           beta_k_o;
  static double       beta_omega_o;
  static double            sigma_k;
  static double        sigma_omega;
  static double              alpha;
  static double                 xi; //!< Compressiblity correction coefficient.
  static double                Mto; //!< Compressiblity correction coefficient.
  //@}

  //@{ @name Propellant variables:
  static double           rhos; //!< Propellant density.
  static double           beta; //!< Propellant burning rate constant.
  static double              n; //!< Propellant burning rate coefficient.
  static double             Tf; //!< Propellant flame temperature.
  static double             Ts; //!< Propellant surface temperature.
  //@}

  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  HighTemp2D_cState(void) {
    Standard_Atmosphere();
  }

  //! Copy constructor.
  HighTemp2D_cState(const HighTemp2D_cState &U) {
    rho = U.rho; dv.x = U.dv.x; dv.y = U.dv.y; E = U.E; dk = U.dk; domega = U.domega;
  }

  //! Copy constructor.
  HighTemp2D_cState(const HighTemp2D_pState &W);

  //! Assignment constructor.
  HighTemp2D_cState(const double &dens,
			const Vector2D &dV,
			const double &Etot) {
    rho = dens; dv.x = dV.x; dv.y = dV.y; E = Etot; dk = ZERO; domega = ZERO;
  }

  //! Assignment constructor.
  HighTemp2D_cState(const double &dens,
			const double &dvx,
			const double &dvy,
			const double &Etot) {
    rho = dens; dv.x = dvx; dv.y = dvy; E = Etot; dk = ZERO; domega = ZERO;
  }

  //! Assignment constructor.
  HighTemp2D_cState(const double &dens,
			const Vector2D &dV,
			const double &Etot,
			const double &dkdk,
			const double &domga) {
    rho = dens; dv.x = dV.x; dv.y = dV.y; E = Etot; dk = dkdk; domega = domga;
  }

  //! Assignment constructor.
  HighTemp2D_cState(const double &dens,
			const double &dvx,
			const double &dvy,
			const double &Etot,
			const double &dkdk,
			const double &domga) {
    rho = dens; dv.x = dvx; dv.y = dvy; E = Etot; dk = dkdk; domega = domga;
  }

  //! Destructor.
  ~HighTemp2D_cState(void) { }
  //@}
  
  //! Return the number of variables.
  int NumVar(void) { return NUM_VAR_HIGHTEMP2D; }

  //@{ @name Set static variables.
  static void set_static_variables(void);
  static void set_static_variables(char *gas_type,
			    const int &EOSType,
			    const int &FlowType,
			    const double &C_constant,
			    const double &von_karman_constant,
			    const double &yplus_sub,
			    const double &yplus_buffer,
			    const double &yplus_outer,
			    char *propellant_type);
  static void set_gas(char *gas_type);
  static void set_turbulence(const double &C_constant,
		      const double &von_karman_constant,
		      const double &yplus_sub,
		      const double &yplus_buffer,
		      const double &yplus_outer);
  static void set_propellant(char *propellant_type);
  //@}

  //@{ @name Useful operators.
  //! Copy operator.
  void Copy(const HighTemp2D_cState &U);

  //! Vacuum operator.
  void Vacuum(void);

  //! Standard atmosphere operator.
  void Standard_Atmosphere(void);

  //! Check for unphysical state properties.
  int Unphysical_Properties(void) const;

  //! Copy variables solved by multigrid only.
  void Copy_Multigrid_State_Variables(const HighTemp2D_cState &Ufine);

  //! Zero variables not-solved by multigrid.
  void Zero_Non_Multigrid_State_Variables(void);
  //@}

  //@{ @name Gas related functions.
  //! Gas flow velocity.
  Vector2D v(void) const;

  //! Specific Heat ratio for high-temperature
  //double gamma(void) const;
	// gamma() is not used so I removed it. Alistair Wood. Wed Aug 01 2007.

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

  //! Gas mean free path.
  //double meanfreepath(void) const;
  //@}

  //! High Temperature dp/de cell-centred calculation.
  double dpde(void) const;

  //! High Temperature dp/drho cell-centred calculation.
  double dpdrho(void) const;

  //! High Temperature dp/drho cell-centred calculation.
  double dTdrho(void) const;

  //! High Temperature dp/drho cell-centred calculation.
  double dTdp(void) const;

  //! High Temperature dp/drho cell-centred calculation.
  double ddTdrho(void) const;

  //! High Temperature dp/drho cell-centred calculation.
  double ddTdp(void) const;

  //! High Temperature dp/drho cell-centred calculation.
  double ddTdpdrho(void) const;

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

  double beta_k(const HighTemp2D_pState &dWdx, const HighTemp2D_pState &dWdy) const;
  double beta_omega(const HighTemp2D_pState &dWdx, const HighTemp2D_pState &dWdy) const;
  double f_beta_k(const HighTemp2D_pState &dWdx, const HighTemp2D_pState &dWdy) const;
  double f_beta_omega(const HighTemp2D_pState &dWdx, const HighTemp2D_pState &dWdy) const;
  double chi_k(const HighTemp2D_pState &dWdx, const HighTemp2D_pState &dWdy) const;
  double chi_omega(const HighTemp2D_pState &dWdx, const HighTemp2D_pState &dWdy) const;
  //@}

  //@{ @name Solid propellant related functions.
  //! Burning rate.
  double burningrate(void) const;
  //@}

  //@{ @name Primitive solution state.
  HighTemp2D_pState W(void) const;
  HighTemp2D_pState W(const HighTemp2D_cState &U) const;
  friend HighTemp2D_pState W(const HighTemp2D_cState &U);
  //@}

  //! Jacobian of the conserved solution variables with respect to the
  //! primitive solution variables.
  void dUdW(DenseMatrix &dUdW) const;

  //! Jacobian of the primitive solution variables with respect to the
  //! conserved solution variables.
  void dWdU(DenseMatrix &dWdU) const;

  //@{ @name Inviscid solution flux (x-direction) and Jacobian.
  HighTemp2D_cState Fx(void) const;
  HighTemp2D_cState Fx(const Vector2D &V) const;
  //HighTemp2D_cState Fx(double dpde, double dpdrho) const; 

	// dFdU is possibly not working right now.
	//  -- Alistair Wood Apr 15 2007 
  //void xdFdU(DenseMatrix &dFdU) const;

  //@}

  //@{ @name Viscous solution flux and Jacobians.
  HighTemp2D_cState Gx(const HighTemp2D_pState &dWdx) const;
  HighTemp2D_cState Gy(const HighTemp2D_pState &dWdy) const;
  //HighTemp2D_cState dGxdU(???) const;
  //HighTemp2D_cState dGydU(???) const;
  //@}

  //! Compute viscous stress tensor and heat flux vector.
  void ComputeViscousTerms(const HighTemp2D_pState &dWdx,
			   const HighTemp2D_pState &dWdy,
			   const Vector2D &X,
			   int Axisymmetric);

  //@{ @name Eigenstructure.
  //! Eigenvalue(s) (x-direction).
  HighTemp2D_pState lambda_x(void) const;
 
  //! Eigenvalue(s) (x-direction).
  HighTemp2D_pState lambda_x(double aAvg) const;

  //! Eigenvalue(s) (x-direction).
  HighTemp2D_pState lambda_x(const Vector2D &V) const;

  //! Primitive right eigenvector (x-direction).
  HighTemp2D_pState rp_x(int index) const;

  //! Conserved right eigenvector (x-direction).
  HighTemp2D_cState rc_x(int index) const;
  HighTemp2D_cState rc_x(int index, double dpde, double dpdrho, double cAvg) const;

  //! Primitive left eigenvector (x-direction).
  HighTemp2D_pState lp_x(int index) const;
  HighTemp2D_pState lp_x(int index, double cAvg) const;
  //@}

  //@{ @name Index operator.
  double &operator[](int index) {
    //assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
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
    // Default return, this is never reached.
    return rho;
  }

  const double &operator[](int index) const {
    //assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
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
    // Default return, this is never reached.
    return rho;
  }
  //@}

  //@{ @name Binary arithmetic operators.
  HighTemp2D_cState operator +(const HighTemp2D_cState &U) const;
  HighTemp2D_cState operator -(const HighTemp2D_cState &U) const;
  double operator *(const HighTemp2D_cState &U) const;
  HighTemp2D_cState operator *(const double &a) const;
  friend HighTemp2D_cState operator *(const double &a, const HighTemp2D_cState &U);
  HighTemp2D_cState operator /(const double &a) const;
  HighTemp2D_cState operator ^(const HighTemp2D_cState &U) const;
  //@}

  //@{ @name Assignment operator.
  HighTemp2D_cState& operator =(const HighTemp2D_cState &U);
  //@}

  //@{ @name Unary arithmetic operators.
  //HighTemp2D_cState operator +(const HighTemp2D_cState &U);
  friend HighTemp2D_cState operator -(const HighTemp2D_cState &U);
  //@}

  //@{ @name Shortcut arithmetic operators.
  HighTemp2D_cState &operator +=(const HighTemp2D_cState &U);
  HighTemp2D_cState &operator -=(const HighTemp2D_cState &U);
  HighTemp2D_cState &operator *=(const double &a);
  HighTemp2D_cState &operator /=(const double &a);
  //@}

  //@{ @name Relational operators.
  friend int operator ==(const HighTemp2D_cState &U1, const HighTemp2D_cState &U2);
  friend int operator !=(const HighTemp2D_cState &U1, const HighTemp2D_cState &U2);
  //@}

  //@{ @name Input-output operators.
  friend ostream &operator << (ostream &out_file, const HighTemp2D_cState &U);
  friend istream &operator >> (istream &in_file, HighTemp2D_cState &U);
  //@}

	/**********************************************************************
	 * Routine: Rotate                                                    *
	 *                                                                    *
	 * This function returns the solution in the local rotated frame.     *
	 *                                                                    *
	 **********************************************************************/
	void Rotate(const HighTemp2D_cState &U, const Vector2D &norm_dir);

};

/**********************************************************************
 * HighTemp2D_pState::Copy -- Copy operator.                          *
 **********************************************************************/
inline void HighTemp2D_pState::Copy(const HighTemp2D_pState &W) {
  rho = W.rho; v.x = W.v.x; v.y = W.v.y; p = W.p; k = W.k; omega = W.omega;
}

/**********************************************************************
 * HighTemp2D_pState::Vacuum -- Vacuum operator.                  *
 **********************************************************************/
inline void HighTemp2D_pState::Vacuum(void) {
  rho = ZERO; v.x = ZERO; v.y = ZERO; p = ZERO; k = ZERO; omega = ZERO; tau.zero(); q.zero();
}

/**********************************************************************
 * HighTemp2D_pState::Standard_Atmosphere -- Standard atmosphere  *
 *                                               operator.            *
 **********************************************************************/
inline void HighTemp2D_pState::Standard_Atmosphere(void) {
  rho = DENSITY_STDATM; v.x = ZERO; v.y = ZERO; p = PRESSURE_STDATM; k = ZERO; omega = ZERO; tau.zero(); q.zero();
}

/**********************************************************************
 * HighTemp2D_pState::Unphysical_Properties -- Check for          *
 *                                                 unphysical state   *
 *                                                 properties.        *
 **********************************************************************/
inline int HighTemp2D_pState::Unphysical_Properties(void) const {
  if (rho <= ZERO || p <= ZERO || E() <= ZERO){
    cout<<" rho = "<<rho<<" p = "<<p<<" E() is "<<E()<<endl;
    return 1;
  }
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) if (k < ZERO || omega < ZERO) return 1;
  return 0; 
}

/**********************************************************************
 * HighTemp2D_pState::set_static_variables -- Set all static      *
 *                                                variables.          *
 **********************************************************************/
inline void HighTemp2D_pState::set_static_variables(void) {
  // Set gas constants.
  set_gas("HTAIR");
  // Set the Equation of State type.
  eos_type = EOS_TGAS;
  // Set the flow type.
  flow_type = FLOWTYPE_LAMINAR;
  // Set turbulence constants.
  set_turbulence(ZERO,ZERO,ZERO,ZERO,ZERO);
  // Set propellant type.
  set_propellant("AP_HTPB");
}

inline void HighTemp2D_pState::set_static_variables(char *gas_type,
						    const int &EOSType,
						    const int &FlowType,
						    const double &C_constant,
						    const double &von_karman_constant,
						    const double &yplus_sub,
						    const double &yplus_buffer,
						    const double &yplus_outer,
						    char *propellant_type) {
  // Set gas constants.
  set_gas(gas_type);
  // Set the equation of state type.
  eos_type = EOSType;
  // Set the flow type.
  flow_type = FlowType;
  // Set turbulence constants.
  set_turbulence(C_constant,von_karman_constant,yplus_sub,yplus_buffer,yplus_outer);
  // Set propellant type.
  set_propellant(propellant_type);
}

/**********************************************************************
 * HighTemp2D_pState::set_gas -- Set gas static variables.        *
 **********************************************************************/
inline void HighTemp2D_pState::set_gas(char *gas_type) {
  if (strcmp(gas_type,"AIR") == 0) {
    g  = GAMMA_AIR;
    R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    v1 = AIR_c1; v2 = AIR_c2; v3 = AIR_c3; v4 = AIR_c4; v5 = AIR_c5;
  } else if (strcmp(gas_type,"HTAIR") == 0) {
    //g  = GAMMA_AIR;
    R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    //v1 = AIR_c1; v2 = AIR_c2; v3 = AIR_c3; v4 = AIR_c4; v5 = AIR_c5;
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
    v1 = AIR_c1; v2 = AIR_c2; v3 = AIR_c3; AIR_c4;
  }

  if (strcmp(gas_type,"HTAIR") != 0) {
    gm1 = g-ONE;
    gm1i = ONE/gm1;
    cp = g*R*gm1i;
    cv = R*gm1i;
  } 
}

/**********************************************************************
 * HighTemp2D_pState::set_turbulence -- Set the turbulence static *
 *                                          variables.                *
 **********************************************************************/
inline void HighTemp2D_pState::set_turbulence(const double &C_constant,
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
 * HighTemp2D_pState::set_propellant -- Set propellant static     *
 *                                          variables.                *
 **********************************************************************/
inline void HighTemp2D_pState::set_propellant(char *propellant_type) {
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
}

//  /**********************************************************************
//   * HighTemp2D_pState::gamma -- Specific Heat Ratio for HTair.         *
//   **********************************************************************/
//  inline double HighTemp2D_pState::gamma(void) const {
//    double temp, g2;
//    temp = Tgas_temp(p,rho);
//    g2 = Tgas_gamma(p,temp);
//    //cout<<"in g function, p state: gamma = "<<g2<<endl; 
//    return g2;
//  }
//
//  gamma() is not used so I removed it. Alistair Wood. Wed Aug 01 2007.

/**********************************************************************
 * HighTemp2D_pState::T -- Gas temperature.                       *
 **********************************************************************/
inline double HighTemp2D_pState::T(void) const {
//    assert(rho > ZERO);
//    double t1 = p/(rho*R);
//    double t2 = Tgas_temp(p,rho);
//  	cout << "eos_type == ";
//  	switch (eos_type) {
//  		case EOS_TGAS:  cout << "EOS_TGAS";  break;
//  		case EOS_IDEAL: cout << "EOS_IDEAL"; break;
//  	}
//  	cout << " using hightemp T = "<<t2<<" and ideal T = "<<t1<<endl; 
 
  switch(eos_type){
  case EOS_TGAS:
    return Tgas_temp(p,rho);
    break;
  case EOS_IDEAL:
    return p/(rho*R);
    break; 
  };
  return p/(rho*R);
}

/**********************************************************************
 * HighTemp2D_pState::e -- Gas specific internal energy.          *
 **********************************************************************/
inline double HighTemp2D_pState::e(void) const {
  //assert(rho > ZERO);

 switch(eos_type){
  case EOS_TGAS:
// Tgas_h is a straight look-up.
// Tgas_e is a solve. Please avoid Tgas_e.
// specific e is [energy/mass]. 
// pressure is [energy/volume].
// p/rho is [energy/mass].
    return Tgas_h(p, rho) - p/rho;
    break;
  case EOS_IDEAL:
    return p/(gm1*rho);
    break;
  };
  return p/(gm1*rho);
}

/**********************************************************************
 * HighTemp2D_pState::E -- Gas total energy.                      *
 **********************************************************************/
inline double HighTemp2D_pState::E(void) const {

  switch(eos_type){
   case EOS_TGAS:
    return rho*e() + HALF*rho*v.sqr() + dk();
    break;
   case EOS_IDEAL:
    return p*gm1i + HALF*rho*v.sqr() + dk();
    break;
  };
  return p*gm1i + HALF*rho*v.sqr() + dk();
}

/**********************************************************************
 * HighTemp2D_pState::h -- Gas specific enthalpy.                 *
 **********************************************************************/
inline double HighTemp2D_pState::h(void) const {
  //assert(rho > ZERO);
  //double h1, h2;
  //h1 = Tgas_h(p,rho) + HALF*v.sqr() + k; 
  //h2 = g*gm1i*p/rho + HALF*v.sqr() + k;  
  
  switch(eos_type){
   case EOS_TGAS:
     //cout<<" in EOS_TGAS, h1 = "<<h1<<" and h2 = "<<h2<<endl;
    return Tgas_h(p,rho) + HALF*v.sqr() + k; 
    break;
   case EOS_IDEAL:
    return g*gm1i*p/rho + HALF*v.sqr() + k;
    break;
  };
  return g*gm1i*p/rho + HALF*v.sqr() + k;
}

/**********************************************************************
 * HighTemp2D_pState::H -- Gas total enthalpy.                    *
 **********************************************************************/
inline double HighTemp2D_pState::H(void) const {
	return h() * rho;
}

/**********************************************************************
 * HighTemp2D_pState::a -- Gas sound speed.                       *
 **********************************************************************/
inline double HighTemp2D_pState::a(void) const {

 switch(eos_type){
   case EOS_TGAS:
    return Tgas_a_from_e_rho(e(), rho);
    break;
   case EOS_IDEAL:
    return sqrt(g*p/rho);
    break;
  };
  return sqrt(g*p/rho);
}

/**********************************************************************
 * HighTemp2D_pState::a2 -- Gas sound speed squared.              *
 **********************************************************************/
inline double HighTemp2D_pState::a2(void) const {
  //assert(rho > ZERO);
  
  switch(eos_type){
   case EOS_TGAS:
		{
			double ax = a(); // in case someone implemented sqr() as a preprocessor define.
			return sqr(ax); 
		 }
    break;
   case EOS_IDEAL:
    return g*p/rho;
    break;
  };
  return g*p/rho;
}

/**********************************************************************
 * HighTemp2D_pState::M -- Gas Mach number.                       *
 **********************************************************************/
inline double HighTemp2D_pState::M(void) const {
  //assert(rho > ZERO && p > ZERO);
  return abs(v)/a();
}

/**********************************************************************
 * HighTemp2D_pState::s -- Gas specific entropy.                  *
 **********************************************************************/
inline double HighTemp2D_pState::s(void) const {
  //assert(rho > ZERO && p > ZERO);

  //double s1, s2;
  // s1 =R*gm1i*log(p/pow(rho,g));
  // s2 = Tgas_s(e(),rho);    

  switch(eos_type){
   case EOS_TGAS:
     // cout<<" in EOS_TGAS, s1 = "<<s1<<" and s2 = "<<s2<<endl;
    return Tgas_s(e(),rho);
    break;
   case EOS_IDEAL:
    return R*gm1i*log(p/pow(rho,g));
    break;
  };
  return R*gm1i*log(p/pow(rho,g));
}

/**********************************************************************
 * HighTemp2D_pState::dv -- Gas momentum.                         *
 **********************************************************************/
inline Vector2D HighTemp2D_pState::dv(void) const {
  return rho*v;
}

/**********************************************************************
 * HighTemp2D_pState::dv -- Gas momentum.                         *
 **********************************************************************/
inline double HighTemp2D_pState::dv(const Vector2D &n) const {
  return rho*(v*n);
}

/**********************************************************************
 * HighTemp2D_pState::To -- Gas stagnation temperature.           *
 **********************************************************************/
inline double HighTemp2D_pState::To(void) const {
  //assert(rho > ZERO && p > ZERO);
  //cout<<"To is called in computation"<<endl;
  return (p/(rho*R))*(ONE+HALF*gm1*v.sqr()/(g*p/rho));
}

/**********************************************************************
 * HighTemp2D_pState::po -- Gas stagnation pressure.              *
 **********************************************************************/
inline double HighTemp2D_pState::po(void) const {
  //assert(rho > ZERO && p > ZERO);
  //cout<<"po is called in computation"<<endl;
  return p*pow(ONE+HALF*gm1*v.sqr()/(g*p/rho),g*gm1i);
}

/**********************************************************************
 * HighTemp2D_pState::ao -- Gas stagnation sound speed.           *
 **********************************************************************/
inline double HighTemp2D_pState::ao(void) const {
  //assert(rho > ZERO && p > ZERO);
  //cout<<"ao is called in computation"<<endl;
  return sqrt((g*p/rho)*(ONE+HALF*gm1*v.sqr()/(g*p/rho)));
}

/**********************************************************************
 * HighTemp2D_pState::ho -- Gas stagnation enthalpy.              *
 **********************************************************************/
inline double HighTemp2D_pState::ho(void) const {
  //assert(rho > ZERO && p > ZERO);
  //cout<<"ho is called in computation"<<endl;
  return (g*p/(gm1*rho) + HALF*v.sqr())*(ONE+HALF*gm1*v.sqr()/(g*p/rho));
}

/**********************************************************************
 * HighTemp2D_pState::mu -- Gas dynamic viscosity.                *
 **********************************************************************/
inline double HighTemp2D_pState::mu(void) const {
  // double temp = Tgas_temp(p,rho);
  // double mu1, mu2;
  // mu1 = Tgas_mu(temp,rho);
  //mu2 = mu_gottlieb(v1,v2,v3,v4,v5,T());

 switch(eos_type){
   case EOS_TGAS:
     // cout<<" in EOS_TGAS, mu1 = "<<mu1<<" and mu2 = "<<mu2<<endl;
    return Tgas_mu(Tgas_temp(p,rho), rho);
    break;
   case EOS_IDEAL:
    return mu_gottlieb(v1,v2,v3,v4,v5,T());
    break;
  };
  return mu_gottlieb(v1,v2,v3,v4,v5,T());
}

/**********************************************************************
 * HighTemp2D_pState::nu -- Gas kinematic viscosity.              *
 **********************************************************************/
inline double HighTemp2D_pState::nu(void) const {
  return mu()/rho;
}

/**********************************************************************
 * HighTemp2D_pState::kappa -- Gas thermal heat conductivity.     *
 **********************************************************************/
inline double HighTemp2D_pState::kappa(void) const {
  
  /*
  double k1, k2;
  k1 =  Tgas_kappa(e(),rho);
  k2 = kappa_gottlieb(v1,v2,v3,v4,v5,T(),g,cp);
  cout<<" in EOS_TGAS, kappa pState, k1 = "<<k1<<" and k2 = "<<k2<<endl;
  */

 switch(eos_type){
   case EOS_TGAS:
     //cout<<" in EOS_TGAS, k1 = "<<k1<<" and k2 = "<<k2<<" and k3 ="<<k3<<endl;
    return Tgas_kappa(e(),rho);
    break;
   case EOS_IDEAL:
    return kappa_gottlieb(v1,v2,v3,v4,v5,T(),g,cp);
    break;
  };
  return kappa_gottlieb(v1,v2,v3,v4,v5,T(),g,cp);
}

/**********************************************************************
 * HighTemp2D_pState::Pr -- Prandtl number.                       *
 **********************************************************************/
inline double HighTemp2D_pState::Pr(void) const {
  //assert(kappa() > ZERO);
 
  switch(eos_type){
   case EOS_TGAS:
    return Tgas_Pr(Tgas_temp(p,rho), rho);
    break;
   case EOS_IDEAL:
    return cp*mu()/kappa();
    break;
  };
  return cp*mu()/kappa();
}

/**********************************************************************
 * HighTemp2D_pState::meanfreepath -- Gas mean free path.         *
 **********************************************************************/
/*
inline double HighTemp2D_pState::meanfreepath(void) const {
  //assert(rho > ZERO && T() > ZERO);
  return 16.0*mu()/(5.0*rho*sqrt(2.0*PI*R*T()));
}
*/
/**********************************************************************
 * HighTemp2D_pState::dk -- Gas total turbulent kinetic energy.   *
 **********************************************************************/
inline double HighTemp2D_pState::dk(void) const {
  return rho*k;
}

/**********************************************************************
 * HighTemp2D_pState::domega -- Gas total turbulent specific      *
 *                                  dissipation.                      *
 **********************************************************************/
inline double HighTemp2D_pState::domega(void) const {
  return rho*omega;
}

/**********************************************************************
 * HighTemp2D_pState::epsilon -- Gas specific turbulent eddy      *
 *                                   dissipation.                     *
 **********************************************************************/
inline double HighTemp2D_pState::epsilon(void) const {
  return beta_k_o*k*omega;
}

/**********************************************************************
 * HighTemp2D_pState::depsilon -- Gas total turbulent eddy        *
 *                                    dissipation.                    *
 **********************************************************************/
inline double HighTemp2D_pState::depsilon(void) const {
  return rho*epsilon();
}

/**********************************************************************
 * HighTemp2D_pState::ell -- Gas turbulent length scale.          *
 **********************************************************************/
inline double HighTemp2D_pState::ell(void) const {
  return sqrt(k)/max(omega,NANO);
}

/**********************************************************************
 * HighTemp2D_pState::Mt -- Gas turbulent Mach number.            *
 **********************************************************************/
inline double HighTemp2D_pState::Mt(void) const {
  return sqrt(TWO*k/a2());
}

/**********************************************************************
 * HighTemp2D_pState::Mt2 -- Turbulent Mach number squared.       *
 **********************************************************************/
inline double HighTemp2D_pState::Mt2(void) const {
  return TWO*k/a2();
}

/**********************************************************************
 * HighTemp2D_pState::muT -- Gas turbulent eddy dynamic viscosity.*
 **********************************************************************/
inline double HighTemp2D_pState::muT(void) const {
  return rho*nuT();
}

/**********************************************************************
 * HighTemp2D_pState::nuT -- Gas turbulent eddy kinematic viscosity.*
 **********************************************************************/
inline double HighTemp2D_pState::nuT(void) const {
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) return k/max(omega,TOLER);
  return ZERO;
}

/**********************************************************************
 * HighTemp2D_pState::kappaT -- Gas turbulent eddy thermal heat   *
 *                                  conductivity.                     *
 **********************************************************************/
inline double HighTemp2D_pState::kappaT(void) const {
  if (flow_type != FLOWTYPE_TURBULENT_RANS_K_OMEGA) return ZERO;

  switch(eos_type){
   case EOS_TGAS:
     //cout<<"cp ideal is = "<<cp<<" and cpHT is = "<<cpHT<<endl;
		 {
		double temp = Tgas_temp(p,rho);
		double cpHT = Tgas_cp(p,temp); 
		// There is an incorrect assumption of a calorically perfect gas 
		// in here somewhere.
    return muT()*cpHT/PrT;  
		 }
    break;
   case EOS_IDEAL:
    return muT()*cp/PrT;
    break;
  };
  return muT()*cp/PrT;
}

/**********************************************************************
 * HighTemp2D_pState::dpde -- dp/de, cell-centred                 *
 **********************************************************************/
inline double HighTemp2D_pState::dpde(void) const {
  
  double elocal, e1, e2; // dppde;
  double ex = e();
  e1 = HTONEPLUST*ex;
  e2 = HTONEMINT*ex;
   
  return (Tgas_p(e1,rho) - Tgas_p(e2,rho))/(2.0*HTTOL*ex);
}

/**********************************************************************
 * HighTemp2D_pState::dpdrho -- dp/drho, cell-centred                 *
 **********************************************************************/
inline double HighTemp2D_pState::dpdrho(void) const {
  
  double rho1, rho2;
  rho1 = HTONEPLUST*rho;
  rho2 = HTONEMINT*rho;
	double ex = e();
   
  return (Tgas_p(ex,rho1) - Tgas_p(ex,rho2))/(2.0*HTTOL*rho);
}

inline double HighTemp2D_pState::dhdrho(void) const {
  return (Tgas_h(p, HTONEPLUST*rho) - Tgas_h(p, HTONEMINT*rho))/(2.0*HTTOL*rho);
}

inline double HighTemp2D_pState::dhdp(void) const { 
  return (Tgas_h(HTONEPLUST*p, rho) - Tgas_h(HTONEMINT*p, rho))/(2.0*HTTOL*p);
}

/**********************************************************************
 * HighTemp2D_pState::dTdp -- dT/dp, cell-centred                 *
 **********************************************************************/
inline double HighTemp2D_pState::dTdp(void) const {
  
  double p1, p2; 
  p1 = HTONEPLUST*p;
  p2 = HTONEMINT*p;  
   
  return (Tgas_temp(p1,rho) - Tgas_temp(p2,rho))/(2.0*HTTOL*p);
}

/**********************************************************************
 * HighTemp2D_pState::dTdrho -- dT/drho, cell-centred                 *
 **********************************************************************/
inline double HighTemp2D_pState::dTdrho(void) const {
  
  double rho1, rho2; 
  rho1 = HTONEPLUST*rho;
  rho2 = HTONEMINT*rho;
   
  return (Tgas_temp(p,rho1) - Tgas_temp(p,rho2))/(2.0*HTTOL*rho);
}
 
/**********************************************************************
 * HighTemp2D_pState::ddTdp -- d^2T/dp^2, cell-centred                 *
 **********************************************************************/
inline double HighTemp2D_pState::ddTdp(void) const {
  
  double p1, p2; 
  p1 = HTONEPLUST*p;
  p2 = HTONEMINT*p;
   
  return (Tgas_temp(p1,rho) -2.0*Tgas_temp(p,rho)+ Tgas_temp(p2,rho))/(pow(HTTOL*p,2.0));
}

/**********************************************************************
 * HighTemp2D_pState::ddTdrho -- d^2T/drho^2, cell-centred                 *
 **********************************************************************/
inline double HighTemp2D_pState::ddTdrho(void) const {
  
  double rho1, rho2; 
  rho1 = HTONEPLUST*rho;
  rho2 = HTONEMINT*rho;
   
  return (Tgas_temp(p,rho1) -2.0*Tgas_temp(p,rho)+ Tgas_temp(p,rho2))/(pow(HTTOL*rho,2.0));
}

/**********************************************************************
 * HighTemp2D_pState::ddTdrho -- d^2T/dpdrho, cell-centred                 *
 **********************************************************************/
inline double HighTemp2D_pState::ddTdpdrho(void) const {
  
  double rho1, rho2, p1, p2;
  double dTdrhoP1, dTdrhoP2;
  rho1 = HTONEPLUST*rho;
  rho2 = HTONEMINT*rho;
  p1 = HTONEPLUST*p;
  p2 = HTONEMINT*p;
   
  dTdrhoP1 = (Tgas_temp(p1,rho1) - Tgas_temp(p1,rho2))/(2.0*HTTOL*rho);
  dTdrhoP2 = (Tgas_temp(p2,rho1) - Tgas_temp(p2,rho2))/(2.0*HTTOL*rho);
 
  return (dTdrhoP1 - dTdrhoP2)/(2.0*HTTOL*p);
}

/**********************************************************************
 * HighTemp2D_pState::c -- Turbulence modified sound speed.       *
 **********************************************************************/
inline double HighTemp2D_pState::c(void) const {
  return sqrt(c2());
}

/**********************************************************************
 * HighTemp2D_pState::c2 -- Turbulence modified sound speed       *
 *                              squared.                              *
 **********************************************************************/
inline double HighTemp2D_pState::c2(void) const {
  //assert(rho > ZERO);

  //double c1, c4;
  //c1 = a2() + (2.0/3.0)*g*k; 
  //c3 = a2() + (2.0/3.0)*ght*k;
  //c4 = a2() + (2.0/3.0)*k + (2.0*k*dpde())/(3.0*rho); 

  //double ght;
  //ght = gamma();

  //double dppdrho;
  //dppdrho = dpdrho();
  //cout<<"dpdrho = "<<dppdrho<<endl;

  switch(eos_type){
   case EOS_TGAS:
     //cout<<" in EOS_TGAS, c4 = "<<c4<<" and c1 ideal = "<<c1<<endl;
     //return a2() + (2.0/3.0)*ght*k;
     return a2() + (2.0/3.0)*k + (2.0*k*dpde())/(3.0*rho);
    break;
   case EOS_IDEAL:
    return a2() + (2.0/3.0)*g*k;
    break;
  };
  return a2() + (2.0/3.0)*g*k;
}

/**********************************************************************
 * HighTemp2D_pState::pmodified -- Turbulence modified pressure.  *
 **********************************************************************/
inline double HighTemp2D_pState::pmodified(void) const {
  //assert(rho > ZERO);
  return p + (2.0/3.0)*dk();
}

/**********************************************************************
 * HighTemp2D_pState::beta_k -- k-omega auxilary relation.        *
 **********************************************************************/
inline double HighTemp2D_pState::beta_k(const HighTemp2D_pState &dWdx,
					    const HighTemp2D_pState &dWdy) const {
   return beta_k_o*f_beta_k(dWdx,dWdy)*(ONE + xi*(Mt2()-sqr(Mto))*heaviside(Mt()-Mto));
  // return beta_k_o*f_beta_k(dWdx,dWdy);
}

/**********************************************************************
 * HighTemp2D_pState::beta_omega -- k-omega auxilary relation.    *
 **********************************************************************/
inline double HighTemp2D_pState::beta_omega(const HighTemp2D_pState &dWdx,
						const HighTemp2D_pState &dWdy) const {
  return beta_omega_o*f_beta_omega(dWdx,dWdy) - beta_k_o*f_beta_k(dWdx,dWdy)*xi*(Mt2()-sqr(Mto))*heaviside(Mt()-Mto);
  // return beta_omega_o*f_beta_omega(dWdx,dWdy);
}

/**********************************************************************
 * HighTemp2D_pState::f_beta_k -- k-omega auxilary relation.      *
 **********************************************************************/
inline double HighTemp2D_pState::f_beta_k(const HighTemp2D_pState &dWdx,
					      const HighTemp2D_pState &dWdy) const {
  double chi = chi_k(dWdx,dWdy);
  if (chi <= ZERO) return ONE;
  return (ONE + 680.0*chi*chi)/(ONE + 400.0*chi*chi);
}

/**********************************************************************
 * HighTemp2D_pState::f_beta_omega -- k-omega auxilary relation.  *
 **********************************************************************/
inline double HighTemp2D_pState::f_beta_omega(const HighTemp2D_pState &dWdx,
						  const HighTemp2D_pState &dWdy) const {
  double chi = chi_omega(dWdx,dWdy);
  return (ONE + 70.0*chi)/(ONE + 80.0*chi);
}

/**********************************************************************
 * HighTemp2D_pState::chi_k -- k-omega auxilary relation.         *
 **********************************************************************/
inline double HighTemp2D_pState::chi_k(const HighTemp2D_pState &dWdx,
					   const HighTemp2D_pState &dWdy) const {
  //return 0.0;
  return (dWdx.k*dWdx.omega + dWdy.k*dWdy.omega)/max(TOLER,cube(omega));
}

/**********************************************************************
 * HighTemp2D_pState::chi_omega -- k-omega auxilary relation.     *
 **********************************************************************/
inline double HighTemp2D_pState::chi_omega(const HighTemp2D_pState &dWdx,
					       const HighTemp2D_pState &dWdy) const {
  //return 0.0;
  return 0.25*fabs(sqr(dWdx.v.y - dWdy.v.x)*(dWdx.v.x + dWdy.v.y)/max(TOLER,cube(beta_omega_o*omega)));
}

/**********************************************************************
 * HighTemp2D_pState::burningrate -- Solid propellent burning rate.*
 **********************************************************************/
inline double HighTemp2D_pState::burningrate(void) const {
  return -beta*pow(p,n);
}

/**********************************************************************
 * HighTemp2D_pState -- Binary arithmetic operators.              *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_pState::operator +(const HighTemp2D_pState &W) const {
  return HighTemp2D_pState(rho+W.rho,v.x+W.v.x,v.y+W.v.y,p+W.p,k+W.k,omega+W.omega);
}

inline HighTemp2D_pState HighTemp2D_pState::operator -(const HighTemp2D_pState &W) const {
  return HighTemp2D_pState(rho-W.rho,v.x-W.v.x,v.y-W.v.y,p-W.p,k-W.k,omega-W.omega);
}

// Inner product operator.
inline double HighTemp2D_pState::operator *(const HighTemp2D_pState &W) const {
  return rho*W.rho + v.x*W.v.x + v.y*W.v.y + p*W.p + k*W.k + omega*W.omega;
}

inline HighTemp2D_pState HighTemp2D_pState::operator *(const double &a) const {
  return HighTemp2D_pState(rho*a,v.x*a,v.y*a,p*a,k*a,omega*a);
}

inline HighTemp2D_pState operator *(const double &a, const HighTemp2D_pState &W) {
  return HighTemp2D_pState(W.rho*a,W.v.x*a,W.v.y*a,W.p*a,W.k*a,W.omega*a);
}

inline HighTemp2D_pState HighTemp2D_pState::operator /(const double &a) const {
  return HighTemp2D_pState(rho/a,v.x/a,v.y/a,p/a,k/a,omega/a);
}

// A useful solution state product operator.
inline HighTemp2D_pState HighTemp2D_pState::operator ^(const HighTemp2D_pState &W) const {
  return HighTemp2D_pState(rho*W.rho,v.x*W.v.x,v.y*W.v.y,p*W.p,k*W.k,omega*W.omega);
}

/**********************************************************************
 * HighTemp2D_pState -- Assignment operator.                      *
 **********************************************************************/
inline HighTemp2D_pState& HighTemp2D_pState::operator =(const HighTemp2D_pState &W) {
  rho = W.rho; v.x = W.v.x; v.y = W.v.y; p = W.p; k = W.k; omega = W.omega;
  return *this;
}

/**********************************************************************
 * HighTemp2D_pState -- Unary arithmetic operators.               *
 **********************************************************************/
//inline HighTemp2D_pState operator +(const HighTemp2D_pState &W) {
//return W;
//}

inline HighTemp2D_pState operator -(const HighTemp2D_pState &W) {
  return HighTemp2D_pState(-W.rho,-W.v.x,-W.v.y,-W.p,-W.k,-W.omega);
}

/**********************************************************************
 * HighTemp2D_pState -- Shortcut arithmetic operators.            *
 **********************************************************************/
inline HighTemp2D_pState& HighTemp2D_pState::operator +=(const HighTemp2D_pState &W) {
  rho += W.rho; v.x += W.v.x; v.y += W.v.y; p += W.p; k += W.k; omega += W.omega;
  return *this;
}

inline HighTemp2D_pState& HighTemp2D_pState::operator -=(const HighTemp2D_pState &W) {
  rho -= W.rho; v.x -= W.v.x; v.y -= W.v.y; p -= W.p; k -= W.k; omega -= W.omega;
  return *this;
}

inline HighTemp2D_pState& HighTemp2D_pState::operator *=(const double &a) {
  rho *= a; v.x *= a; v.y *= a; p *= a; k *= a; omega *= a;
  return *this;
}

inline HighTemp2D_pState& HighTemp2D_pState::operator /=(const double &a) {
  rho /= a; v.x /= a; v.y /= a; p /= a; k /= a; omega /= a;
  return *this;
}

/**********************************************************************
 * HighTemp2D_pState -- Relational operators.                     *
 **********************************************************************/
inline int operator ==(const HighTemp2D_pState &W1, const HighTemp2D_pState &W2) {
  return (W1.rho == W2.rho && W1.v == W2.v && W1.p == W2.p && W1.k == W2.k && W1.omega == W2.omega);
}

inline int operator !=(const HighTemp2D_pState &W1, const HighTemp2D_pState &W2) {
  return (W1.rho != W2.rho || W1.v != W2.v || W1.p != W2.p || W1.k != W2.k || W1.omega != W2.omega);
}

/**********************************************************************
 * HighTemp2D_pState -- Input-output operators.                   *
 **********************************************************************/
inline ostream &operator << (ostream &out_file, const HighTemp2D_pState &W) {
  out_file.setf(ios::scientific);
  out_file << " " << W.rho << " " << W.v.x << " " << W.v.y << " " << W.p << " " << W.k << " " << W.omega;
  out_file.unsetf(ios::scientific);
  return out_file;
}

inline istream &operator >> (istream &in_file, HighTemp2D_pState &W) {
  in_file.setf(ios::skipws);
  in_file >> W.rho >> W.v.x >> W.v.y >> W.p >> W.k >> W.omega;
  in_file.unsetf(ios::skipws);
  return in_file;
}

/**********************************************************************
 * HighTemp2D_cState::Copy -- Copy operator.                      *
 **********************************************************************/
inline void HighTemp2D_cState::Copy(const HighTemp2D_cState &U) {
  rho = U.rho; dv.x = U.dv.x; dv.y = U.dv.y; E = U.E; dk = U.dk; domega = U.domega;
}

/**********************************************************************
 * HighTemp2D_cState::Vacuum -- Vacuum state.                     *
 **********************************************************************/
inline void HighTemp2D_cState::Vacuum(void) {
  rho = ZERO; dv.x = ZERO; dv.y = ZERO; E = ZERO; dk = ZERO; domega = ZERO; tau.zero(); q.zero();
}

/**********************************************************************
 * HighTemp2D_cState::Standard_Atmosphere -- Standard atmosphere  *
 *                                               state.               *
 **********************************************************************/
inline void HighTemp2D_cState::Standard_Atmosphere(void) {
  rho = DENSITY_STDATM; dv.x = ZERO; dv.y = ZERO;
  E = PRESSURE_STDATM/(GAMMA_AIR-ONE); tau.zero(); q.zero();
  dk = ZERO; domega = ZERO;
}

/**********************************************************************
 * HighTemp2D_cState::Unphysical_Properties -- Check for          *
 *                                                 unphysical state   *
 *                                                 properties.        *
 **********************************************************************/
inline int HighTemp2D_cState::Unphysical_Properties(void) const {
  if (rho <= ZERO || E <= ZERO || e() <= ZERO){
    cout<<" rho = "<<rho<<" E = "<<E<<" e() is "<<e()<<endl;
    return 1;
  }
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) if (dk < ZERO || domega < ZERO) return 1;
  return 0;
}

/**********************************************************************
 * HighTemp2D_cState::Copy_Multigrid_State_Variables --           *
 *                           Copy variables solved by multigrid only. *
 **********************************************************************/
inline void HighTemp2D_cState::Copy_Multigrid_State_Variables(const HighTemp2D_cState &Ufine) {
  rho = Ufine.rho; dv.x = Ufine.dv.x; dv.y = Ufine.dv.y; E = Ufine.E;
}

/**********************************************************************
 * HighTemp2D_cState::Zero_Non_Multigrid_State_Variables --       *
 *                            Zero variables not-solved by multigrid. *
 **********************************************************************/
inline void HighTemp2D_cState::Zero_Non_Multigrid_State_Variables(void) {
  dk = ZERO; domega = ZERO;
}

/**********************************************************************
 * HighTemp2D_cState::set_static_variables -- Set all static      *
 *                                                variables.          *
 **********************************************************************/
inline void HighTemp2D_cState::set_static_variables(void) {
  // Set gas constants.
  set_gas("HTAIR");
  // Set the EOS Type
  eos_type = EOS_TGAS;
  // Set the flow type.
  flow_type = FLOWTYPE_LAMINAR;
  // Set turbulence constants.
  set_turbulence(ZERO,ZERO,ZERO,ZERO,ZERO);
  // Set propellant type.
  set_propellant("AP_HTPB");
}

inline void HighTemp2D_cState::set_static_variables(char *gas_type,
						        const int &EOSType,
							const int &FlowType,
							const double &C_constant,
							const double &von_karman_constant,
							const double &yplus_sub,
							const double &yplus_buffer,
							const double &yplus_outer,
							char *propellant_type) {
  // Set gas constants.
  set_gas(gas_type);
  // Set Equation of State Type
  eos_type = EOSType;
  // Set the flow type.
  flow_type = FlowType;
  // Set turbulence constants.
  set_turbulence(C_constant,von_karman_constant,yplus_sub,yplus_buffer,yplus_outer);
  // Set propellant type.
  set_propellant(propellant_type);
}

/**********************************************************************
 * HighTemp2D_cState::set_gas -- Set gas static variables.         *
 **********************************************************************/
inline void HighTemp2D_cState::set_gas(char *gas_type) {
  if (strcmp(gas_type,"AIR") == 0) {
    g  = GAMMA_AIR;
    R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    v1 = AIR_c1; v2 = AIR_c2; v3 = AIR_c3; v4 = AIR_c4; v5 = AIR_c5;
  }else if (strcmp(gas_type,"HTAIR") == 0) {
    //g  = GAMMA_AIR;
    //R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    //v1 = AIR_c1; v2 = AIR_c2; v3 = AIR_c3; v4 = AIR_c4; v5 = AIR_c5;
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
 
 if (strcmp(gas_type,"HTAIR") != 0) {
    gm1 = g-ONE;
    gm1i = ONE/gm1;
    cp = g*R*gm1i;
    cv = R*gm1i;
    //cout<<"gas_type is not HTAIR"<<endl;
 } //else cout<<"HTAIR gas, no cp, cv, gm1, gm1i defined!"<<endl;
}

/**********************************************************************
 * HighTemp2D_cState::set_turbulence -- Set the turbulence static *
 *                                          variables.                *
 **********************************************************************/
inline void HighTemp2D_cState::set_turbulence(const double &C_constant,
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
 * HighTemp2D_cState::set_propellant -- Set propellant static     *
 *                                          variables.                *
 **********************************************************************/
inline void HighTemp2D_cState::set_propellant(char *propellant_type) {
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
}

/**********************************************************************
 * HighTemp2D_cState::v -- Gas flow velocity.                     *
 **********************************************************************/
inline Vector2D HighTemp2D_cState::v(void) const {
  //assert(rho > ZERO);
  return dv/rho;
}

/**********************************************************************
 * HighTemp2D_cState::v -- Gas flow velocity.                     *
**********************************************************************/
inline double HighTemp2D_cState::v(const Vector2D &n) const {
  //assert(rho > ZERO);
  return (dv*n)/rho;
}

/**********************************************************************
 * HighTemp2D_cState::p -- Gas pressure.                          *
 **********************************************************************/
inline double HighTemp2D_cState::p(void) const {
  //assert(rho > ZERO);

  switch(eos_type){
   case EOS_TGAS:
		 {
    return Tgas_p(e(), rho);
		 }
    break;
   case EOS_IDEAL:
    return gm1*(E - HALF*dv.sqr()/rho - dk);
    break;
  }
  return gm1*(E - HALF*dv.sqr()/rho - dk);
}

//  /**********************************************************************
//   * HighTemp2D_cState::g -- Specific Heat Ratio for high temp air.     *
//   **********************************************************************/
//  inline double HighTemp2D_cState::gamma(void) const{
//    double temp, g2;
//    temp = Tgas_temp(p(),rho);
//    g2 = Tgas_gamma(p(),temp);
//    //cout<<"in g function, c state class: gamma = "<<g2<<endl;
//    return g2; 
//  }
//
//  gamma() is not used so I removed it. Alistair Wood. Wed Aug 01 2007.

/**********************************************************************
 * HighTemp2D_cState::T -- Gas temperature.                       *
 **********************************************************************/
inline double HighTemp2D_cState::T(void) const {
  //assert(rho > ZERO);

 //double t1, t2;
  //t1 =  Tgas_temp(p(),rho);
  //t2 = p()/(rho*R);

 switch(eos_type){
  case EOS_TGAS:
    //cout<<"EOS_TGAS true, using t2 ="<<t2<<" and t1 = "<<t1<<endl; 
    return Tgas_temp(p(),rho);
    break;
  case EOS_IDEAL:
    return p()/(rho*R);
    break; 
  };
  return p()/(rho*R);
}

/**********************************************************************
 * HighTemp2D_cState::e -- Gas specific internal energy.          *
 **********************************************************************/
inline double HighTemp2D_cState::e(void) const {
  //assert(rho > ZERO);

	return (E - HALF*dv.sqr()/rho - dk)/rho;

	// Previously we were jumping hoops (calling three(?) equation-of-
	// state functions) to determine internal energy. But energy is
	// a conserved variable so I rewrote this function to simply
	// return total energy minus kinetic energy. 
	//   -- Alistair Wood Tue Jul 10 2007 
  
//   double temp;
//    temp = Tgas_temp(p(),rho);
//    // double e1, e2;
//    // e1 =  Tgas_e(p(),temp);
//    //e2 = gm1i*p()/rho;
//  
//   switch(eos_type){
//    case EOS_TGAS:
//      // cout<<" in EOS_TGAS, e1 = "<<e1<<endl;
//      return Tgas_e(p(),temp);
//      break;
//    case EOS_IDEAL:
//      return gm1i*p()/rho;
//      break;
//    };
//    return gm1i*p()/rho;
}

/**********************************************************************
 * HighTemp2D_cState::h -- Gas specific total enthalpy.                 *
 **********************************************************************/
inline double HighTemp2D_cState::h(void) const {
  //assert(rho > ZERO);

	return H() / rho;

	// See the comment for HighTemp2D_cState::H().
	//   -- Alistair Wood Tue Jul 10 2007 

//    //double h1, h2;
//    //h1 = Tgas_h(p(),rho)+ HALF*dv.sqr()/sqr(rho) + k();
//    //h2 = g*gm1i*p()/rho + HALF*dv.sqr()/sqr(rho) + k(); 
//   
//   switch(eos_type){
//    case EOS_TGAS:
//      //cout<<" in EOS_TGAS, h1 = "<<h1<<" and h2 = "<<h2<<endl;
//      return Tgas_h(p(),rho)+ HALF*dv.sqr()/sqr(rho) + k();
//      break;
//    case EOS_IDEAL:
//      return g*gm1i*p()/rho + HALF*dv.sqr()/sqr(rho) + k();
//      break;
//    };
//    return g*gm1i*p()/rho + HALF*dv.sqr()/sqr(rho) + k();
}

/**********************************************************************
 * HighTemp2D_cState::H -- Gas total enthalpy.                    *
 **********************************************************************/
inline double HighTemp2D_cState::H(void) const {

	return E + p();

	// Energy (per volume) is a conserved variable and so is given
	// by E. Enthalpy is a definition which is not a function of
	// the gas type: H = E + p [J / m^3]. Below is what we used to do.
	//   -- Alistair Tue Jul 10 2007 

//    double H1, H2;
//    // H1 = rho*(Tgas_h(p(),rho))+ HALF*dv.sqr()/rho + dk;  
//    // H2 = g*gm1i*p() + HALF*dv.sqr()/rho + dk;
//  
//    switch(eos_type){
//    case EOS_TGAS:
//      //   cout<<" in EOS_TGAS, H1 = "<<H1<<" and H2 = "<<H2<<endl;
//      return rho*(Tgas_h(p(),rho))+ HALF*dv.sqr()/rho + dk; 
//      break;
//    case EOS_IDEAL:
//      return g*gm1i*p() + HALF*dv.sqr()/rho + dk;
//      break;
//    };
//    return g*gm1i*p() + HALF*dv.sqr()/rho + dk;
}

/**********************************************************************
 * HighTemp2D_cState::a -- Gas sound speed.                       *
 **********************************************************************/
inline double HighTemp2D_cState::a(void) const {
  //assert(rho > ZERO);
  
 switch(eos_type){
   case EOS_TGAS:  return Tgas_a_from_e_rho(e(), rho); break;
   case EOS_IDEAL: return sqrt(g*p()/rho);             break;
  };
  return sqrt(g*p()/rho);
}

/**********************************************************************
 * HighTemp2D_cState::a2 -- Gas sound speed squared.              *
 **********************************************************************/
inline double HighTemp2D_cState::a2(void) const {
  //assert(rho > ZERO);
  
  switch(eos_type){
   case EOS_TGAS:
		{
			double ax = a();
			return sqr(ax);
		 }
    break;
   case EOS_IDEAL:
    return g*p()/rho;
    break;
  };
  return g*p()/rho;
}

/**********************************************************************
 * HighTemp2D_cState::M -- Gas Mach number.                       *
 **********************************************************************/
inline double HighTemp2D_cState::M(void) const {
  //assert(rho > ZERO);
  return abs(v())/a();
}

/**********************************************************************
 * HighTemp2D_cState::s -- Gas specific entropy.                  *
 **********************************************************************/
inline double HighTemp2D_cState::s(void) const {
  //assert(rho > ZERO);

  //double s1, s2;
  // s1 =R*gm1i*log(p()/pow(rho,g));
  //s2 = Tgas_s(e(),rho);    

  switch(eos_type){
   case EOS_TGAS:
     // cout<<" in EOS_TGAS, s1 = "<<s1<<" and s2 = "<<s2<<endl;
    return Tgas_s(e(),rho);
    break;
   case EOS_IDEAL:
    return R*gm1i*log(p()/pow(rho,g));
    break;
  };

  return R*gm1i*log(p()/pow(rho,g));
  //return R*gm1i*log(gm1*(E - HALF*dv.sqr()/rho)/pow(rho,g));
}

/**********************************************************************
 * HighTemp2D_cState::To -- Gas stagnation temperature.           *
 **********************************************************************/
inline double HighTemp2D_cState::To(void) const {
  //assert(rho > ZERO);
  return (gm1*(E - HALF*dv.sqr()/rho)/(rho*R))*
	 (ONE+HALF*gm1*dv.sqr()/(rho*rho*g*gm1*(E/rho - HALF*dv.sqr()/sqr(rho))));
}

/**********************************************************************
 * HighTemp2D_cState::po -- Gas stagnation pressure.              *
 **********************************************************************/
inline double HighTemp2D_cState::po(void) const {
  return (gm1*(E - HALF*dv.sqr()/rho))*
	  pow(ONE+HALF*gm1*dv.sqr()/(rho*rho*g*gm1*(E/rho - HALF*dv.sqr()/sqr(rho))),g*gm1i);
}

/**********************************************************************
 * HighTemp2D_cState::ao -- Gas stagnation sound speed.           *
 **********************************************************************/
inline double HighTemp2D_cState::ao(void) const {
  return sqrt((g*gm1*(E/rho - HALF*dv.sqr()/sqr(rho)))*
	      (ONE+HALF*gm1*dv.sqr()/(rho*rho*g*gm1*(E/rho - HALF*dv.sqr()/sqr(rho)))));
}

/**********************************************************************
 * HighTemp2D_cState::ho -- Gas stagnation enthalpy.              *
 **********************************************************************/
inline double HighTemp2D_cState::ho(void) const {
  return (g*E/rho - gm1*HALF*dv.sqr()/sqr(rho))*
         (ONE+HALF*gm1*dv.sqr()/(rho*rho*g*gm1*(E/rho - HALF*dv.sqr()/sqr(rho))));
}

/**********************************************************************
 * HighTemp2D_cState::mu -- Gas dynamic viscosity.                *
 **********************************************************************/
inline double HighTemp2D_cState::mu(void) const {

  //double temp = Tgas_temp(p(),rho);
  //double mu1, mu2;
  //mu1 = Tgas_mu(temp,rho);
  //mu2 = mu_gottlieb(v1,v2,v3,v4,v5,T());

  switch(eos_type){
   case EOS_TGAS:
     //cout<<" in EOS_TGAS, mu1 = "<<mu1<<" and mu2 = "<<mu2<<endl;
    return Tgas_mu(Tgas_temp(p(),rho), rho);
    break;
   case EOS_IDEAL:
    return mu_gottlieb(v1,v2,v3,v4,v5,T());
    break;
  };


  return mu_gottlieb(v1,v2,v3,v4,v5,T());
}

/**********************************************************************
 * HighTemp2D_cState::nu -- Gas kinematic viscosity.              *
 **********************************************************************/
inline double HighTemp2D_cState::nu(void) const {
  return mu()/rho;
}

/**********************************************************************
 * HighTemp2D_cState::kappa -- Gas thermal heat conductivity.     *
 **********************************************************************/
inline double HighTemp2D_cState::kappa(void) const {
  /* 
  double k1, k2;
  k1 =  Tgas_kappa(e(),rho);
  k2 = kappa_gottlieb(v1,v2,v3,v4,v5,T(),g,cp);
  */
  switch(eos_type){
   case EOS_TGAS:
     //cout<<" in EOS_TGAS, k1 = "<<k1<<" and k2 = "<<k2<<endl;
    return Tgas_kappa(e(),rho);
    break;
   case EOS_IDEAL:
    return kappa_gottlieb(v1,v2,v3,v4,v5,T(),g,cp);
    break;
  };
  return kappa_gottlieb(v1,v2,v3,v4,v5,T(),g,cp);
}

/**********************************************************************
 * HighTemp2D_cState::Pr -- Prandtl number.                       *
 **********************************************************************/
inline double HighTemp2D_cState::Pr(void) const {
 
  // double temp = Tgas_temp(p(),rho);
  /* double Pr1, Pr2;
  Pr1 = Tgas_Pr(temp,rho);
  Pr2 = cp*mu()/kappa(); 
  */
  switch(eos_type){
   case EOS_TGAS:
     //cout<<" in EOS_TGAS, Pr1 = "<<Pr1<<" and Pr2 = "<<Pr2<<endl;
    return Tgas_Pr(Tgas_temp(p(),rho), rho);
    break; 
   case EOS_IDEAL:
    return cp*mu()/kappa();
    break;
  };
  return cp*mu()/kappa();
}

/**********************************************************************
 * HighTemp2D_cState:meanfreepath -- Gas mean free path.          *
 **********************************************************************/
/*
inline double HighTemp2D_cState::meanfreepath(void) const {
  //assert(rho > ZERO && T() > ZERO);
  return 16.0*mu()/(5.0*rho*sqrt(2.0*PI*R*T()));
}
*/

/**********************************************************************
 * HighTemp2D_cState::dpde -- dp/de, cell-centred                 *
 **********************************************************************/
inline double HighTemp2D_cState::dpde(void) const {
  
  double e1, e2; // dppde;
  e1 = HTONEPLUST*e();
  e2 = HTONEMINT*e();
   
  return (Tgas_p(e1,rho) - Tgas_p(e2,rho))/(2.0*HTTOL*e());
}

/**********************************************************************
 * HighTemp2D_cState::dpdrho -- dp/drho, cell-centred                 *
 **********************************************************************/
inline double HighTemp2D_cState::dpdrho(void) const {
  
  double rho1, rho2, dppdrho;
  rho1 = HTONEPLUST*rho;
  rho2 = HTONEMINT*rho;
   
  return (Tgas_p(e(),rho1) - Tgas_p(e(),rho2))/(2.0*HTTOL*rho);
}

/**********************************************************************
 * HighTemp2D_pState::dTdp -- dT/dp, cell-centred                 *
 **********************************************************************/
inline double HighTemp2D_cState::dTdp(void) const {
  
	double px = p();
  double p1 = HTONEPLUST*px;
  double p2 = HTONEMINT*px;
   
  return (Tgas_temp(p1,rho) - Tgas_temp(p2,rho))/(2.0*HTTOL*px);
}

/**********************************************************************
 * HighTemp2D_pState::dTdrho -- dT/drho, cell-centred                 *
 **********************************************************************/
inline double HighTemp2D_cState::dTdrho(void) const {
  
  double rho1, rho2; 
  rho1 = HTONEPLUST*rho;
  rho2 = HTONEMINT*rho;
	double px = p();
   
  return (Tgas_temp(px,rho1) - Tgas_temp(px,rho2))/(2.0*HTTOL*rho);
}
 
/**********************************************************************
 * HighTemp2D_pState::ddTdp -- d^2T/dp^2, cell-centred                 *
 **********************************************************************/
inline double HighTemp2D_cState::ddTdp(void) const {
  
  double p1, p2; 
	double px = p();
  p1 = HTONEPLUST*px;
  p2 = HTONEMINT*px;
   
  return (Tgas_temp(p1,rho) -2.0*Tgas_temp(px,rho)+ Tgas_temp(p2,rho))/sqr(HTTOL*px);
}

/**********************************************************************
 * HighTemp2D_pState::ddTdrho -- d^2T/drho^2, cell-centred                 *
 **********************************************************************/
inline double HighTemp2D_cState::ddTdrho(void) const {
  
  double rho1, rho2; 
	double px = p();
  rho1 = HTONEPLUST*rho;
  rho2 = HTONEMINT*rho;
   
  return (Tgas_temp(px,rho1) -2.0*Tgas_temp(px,rho)+ Tgas_temp(px,rho2))/sqr(HTTOL*rho);
}

/**********************************************************************
 * HighTemp2D_pState::ddTdrho -- d^2T/dpdrho, cell-centred                 *
 **********************************************************************/
inline double HighTemp2D_cState::ddTdpdrho(void) const {
  
  double rho1, rho2, p1, p2;
  double dTdrhoP1, dTdrhoP2;
	double px = p();
  rho1 = HTONEPLUST*rho;
  rho2 = HTONEMINT*rho;
  p1 = HTONEPLUST*px;
  p2 = HTONEMINT*px;
   
  dTdrhoP1 = (Tgas_temp(p1,rho1) - Tgas_temp(p1,rho2))/(2.0*HTTOL*rho);
  dTdrhoP2 = (Tgas_temp(p2,rho1) - Tgas_temp(p2,rho2))/(2.0*HTTOL*rho);
 
    return (dTdrhoP1 - dTdrhoP2)/(2.0*HTTOL*px);
}

/**********************************************************************
 * HighTemp2D_cState::dk -- Gas specific turbulent kinetic energy.*
 **********************************************************************/
inline double HighTemp2D_cState::k(void) const {
  return dk/rho;
}

/**********************************************************************
 * HighTemp2D_cState::depsilon -- Gas total turbulent eddy        *
 *                                    dissipation.                    *
 **********************************************************************/
inline double HighTemp2D_cState::depsilon(void) const {
  return beta_k_o*dk*domega/rho;
}

/**********************************************************************
 * HighTemp2D_cState::depsilon -- Gas specific turbulent eddy     *
 *                                    dissipation.                    *
 **********************************************************************/
inline double HighTemp2D_cState::epsilon(void) const {
  return depsilon()/rho;
}

/**********************************************************************
 * HighTemp2D_cState::omega -- Gas specific turbulent dissipation *
 *                                 rate.                              *
 **********************************************************************/
inline double HighTemp2D_cState::omega(void) const {
  return domega/rho;
}

/**********************************************************************
 * HighTemp2D_cState::ell -- Return the turbulent length scale.   *
 **********************************************************************/
inline double HighTemp2D_cState::ell(void) const {
  return sqrt(k())/max(omega(),NANO);
}

/**********************************************************************
 * HighTemp2D_cState::Mt -- Return the turbulent Mach number.     *
 **********************************************************************/
inline double HighTemp2D_cState::Mt(void) const {
  return sqrt(TWO*k()/a2());
}

/**********************************************************************
 * HighTemp2D_cState::Mt2 -- Turbulent Mach number squared.       *
 **********************************************************************/
inline double HighTemp2D_cState::Mt2(void) const {
  return TWO*k()/a2();
}

/**********************************************************************
 * HighTemp2D_cState::muT -- Turbulent eddy dynamic viscosity.    *
 **********************************************************************/
inline double HighTemp2D_cState::muT(void) const {
  return rho*nuT();
}

/**********************************************************************
 * HighTemp2D_cState::nuT -- Turbulent eddy kinematic viscosity.  *
 **********************************************************************/
inline double HighTemp2D_cState::nuT(void) const {
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) return k()/max(omega(),TOLER);
  return ZERO;
}

/**********************************************************************
 * HighTemp2D_cState::kappaT -- Turbulent eddy thermal heat       *
 *                                  conductivity.                     *
 **********************************************************************/
inline double HighTemp2D_cState::kappaT(void) const {
  if (flow_type != FLOWTYPE_TURBULENT_RANS_K_OMEGA) return ZERO;

  switch(eos_type){
   case EOS_TGAS:
		 {
			 double px = p();
		double temp = Tgas_temp(px, rho);
		double cpHT = Tgas_cp(px, temp); 
     //cout<<"cp ideal is = "<<cp<<" and cpHT is = "<<cpHT<<endl;
    return muT()*cpHT/PrT;
		 }
    break;
   case EOS_IDEAL:
    return muT()*cp/PrT;
    break;
  };
  return muT()*cp/PrT;
}

/**********************************************************************
 * HighTemp2D_cState::c -- Turbulence modified sound speed.       *
 **********************************************************************/
inline double HighTemp2D_cState::c(void) const {
  return sqrt(c2());
}

/**********************************************************************
 * HighTemp2D_cState::c2 -- Turbulence modified sound speed       *
 *                              squared.                              *
 **********************************************************************/
inline double HighTemp2D_cState::c2(void) const {
  //assert(rho > ZERO);

  //double c1, c4;
  //c1 = a2() + (2.0/3.0)*g*k(); 
  //c3 = a2() + (2.0/3.0)*ght*k();
  //c4 = a2() + (2.0/3.0)*k() + (2.0*k()*dpde())/(3.0*rho); 

  //double ght;
  //ght = gamma();

  switch(eos_type){
   case EOS_TGAS:
     //cout<<" in EOS_TGAS cState, c4 = "<<c4<<" and c1 ideal = "<<c1<<endl;
     //return a2() + (2.0/3.0)*ght*k();
     return a2() + (2.0/3.0)*k() + (2.0*k()*dpde())/(3.0*rho);
    break;
   case EOS_IDEAL:
    return a2() + (2.0/3.0)*g*k();
    break;
  };
  return a2() + (2.0/3.0)*g*k();
}

/**********************************************************************
 * HighTemp2D_cState::pmodified -- Turbulence modified pressure.  *
 **********************************************************************/
inline double HighTemp2D_cState::pmodified(void) const {
  return p() + (2.0/3.0)*dk;
}

/**********************************************************************
 * HighTemp2D_cState::beta_k -- k-omega auxilary relation.        *
 **********************************************************************/
inline double HighTemp2D_cState::beta_k(const HighTemp2D_pState &dWdx,
					    const HighTemp2D_pState &dWdy) const {
  return beta_k_o*f_beta_k(dWdx,dWdy)*(ONE + xi*(Mt2()-sqr(Mto))*heaviside(Mt()-Mto));
  // return beta_k_o*f_beta_k(dWdx,dWdy);
}

/**********************************************************************
 * HighTemp2D_cState::beta_omega -- k-omega auxilary relation.    *
 **********************************************************************/
inline double HighTemp2D_cState::beta_omega(const HighTemp2D_pState &dWdx,
						const HighTemp2D_pState &dWdy) const {
  return beta_omega_o*f_beta_omega(dWdx,dWdy) - beta_k_o*f_beta_k(dWdx,dWdy)*xi*(Mt2()-sqr(Mto))*heaviside(Mt()-Mto);
  // return beta_omega_o*f_beta_omega(dWdx,dWdy);
}

/**********************************************************************
 * HighTemp2D_cState::f_beta_k -- k-omega auxilary relation.      *
 **********************************************************************/
inline double HighTemp2D_cState::f_beta_k(const HighTemp2D_pState &dWdx,
					      const HighTemp2D_pState &dWdy) const {
  double chi = chi_k(dWdx,dWdy);
  if (chi <= ZERO) return ONE;
  return (ONE + 680.0*chi*chi)/(ONE + 400.0*chi*chi);
}

/**********************************************************************
 * HighTemp2D_cState::f_beta_omega -- k-omega auxilary relation.  *
 **********************************************************************/
inline double HighTemp2D_cState::f_beta_omega(const HighTemp2D_pState &dWdx,
						  const HighTemp2D_pState &dWdy) const {
  double chi = chi_omega(dWdx,dWdy);
  return (ONE + 70.0*chi)/(ONE + 80.0*chi);
}

/**********************************************************************
 * HighTemp2D_cState::chi_k -- k-omega auxilary relation.         *
 **********************************************************************/
inline double HighTemp2D_cState::chi_k(const HighTemp2D_pState &dWdx,
					   const HighTemp2D_pState &dWdy) const {
  //return 0.0;
  return (dWdx.k*dWdx.omega + dWdy.k*dWdy.omega)/max(TOLER,cube(omega()));
}

/**********************************************************************
 * HighTemp2D_cState::chi_omega -- k-omega auxilary relation.     *
 **********************************************************************/
inline double HighTemp2D_cState::chi_omega(const HighTemp2D_pState &dWdx,
					       const HighTemp2D_pState &dWdy) const {
  //return 0.0;
  return 0.25*fabs(sqr(dWdx.v.y - dWdy.v.x)*(dWdx.v.x + dWdy.v.y)/max(TOLER,cube(beta_omega_o*omega())));
}

/**********************************************************************
 * HighTemp2D_cState::burningrate -- Solid propellent burning rate.*
 **********************************************************************/
inline double HighTemp2D_cState::burningrate(void) const {
  return -beta*pow(p(),n);
}

/**********************************************************************
 * HighTemp2D_cState -- Binary arithmetic operators.              *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_cState::operator +(const HighTemp2D_cState &U) const {
  return HighTemp2D_cState(rho+U.rho,dv.x+U.dv.x,dv.y+U.dv.y,E+U.E,dk+U.dk,domega+U.domega);
}

inline HighTemp2D_cState HighTemp2D_cState::operator -(const HighTemp2D_cState &U) const {
  return HighTemp2D_cState(rho-U.rho,dv.x-U.dv.x,dv.y-U.dv.y,E-U.E,dk-U.dk,domega-U.domega);
}

// Inner product operator.
inline double HighTemp2D_cState::operator *(const HighTemp2D_cState &U) const {
  return rho*U.rho + dv.x*U.dv.x + dv.y*U.dv.y + E*U.E + dk*U.dk + domega*U.domega;
}

inline HighTemp2D_cState HighTemp2D_cState::operator *(const double &a) const {
  return HighTemp2D_cState(rho*a,dv.x*a,dv.y*a,E*a,dk*a,domega*a);
}

inline HighTemp2D_cState operator *(const double &a, const HighTemp2D_cState &U) {
  return HighTemp2D_cState(U.rho*a,U.dv.x*a,U.dv.y*a,U.E*a,U.dk*a,U.domega*a);
}

inline HighTemp2D_cState HighTemp2D_cState::operator /(const double &a) const {
  return HighTemp2D_cState(rho/a,dv.x/a,dv.y/a,E/a,dk/a,domega/a);
}

// A useful solution state product operator.
inline HighTemp2D_cState HighTemp2D_cState::operator ^(const HighTemp2D_cState &U) const {
  return HighTemp2D_cState(rho*U.rho,dv.x*U.dv.x,dv.y*U.dv.y,E*U.E,dk*U.dk,domega*U.domega);
}

/**********************************************************************
 * HighTemp2D_cState -- Assignment operator.                      *
 **********************************************************************/
inline HighTemp2D_cState& HighTemp2D_cState::operator =(const HighTemp2D_cState &U) {
  //if (this != &U) {
  rho = U.rho; dv.x = U.dv.x; dv.y = U.dv.y; E = U.E; dk = U.dk; domega = U.domega;
  //}
  return *this;
}

/**********************************************************************
 * HighTemp2D_cState -- Unary arithmetic operators.               *
 **********************************************************************/
//inline HighTemp2D_cState operator +(const HighTemp2D_cState &U) {
//return U;
//}

inline HighTemp2D_cState operator -(const HighTemp2D_cState &U) {
  return HighTemp2D_cState(-U.rho,-U.dv.x,-U.dv.y,-U.E,-U.dk,-U.domega);
}

/**********************************************************************
 * HighTemp2D_cState -- Shortcut arithmetic operators.            *
 **********************************************************************/
inline HighTemp2D_cState& HighTemp2D_cState::operator +=(const HighTemp2D_cState &U) {
  rho += U.rho; dv.x += U.dv.x; dv.y += U.dv.y; E += U.E; dk += U.dk; domega += U.domega;
  return *this;
}

inline HighTemp2D_cState& HighTemp2D_cState::operator -=(const HighTemp2D_cState &U) {
  rho -= U.rho; dv.x -= U.dv.x; dv.y -= U.dv.y; E -= U.E; dk -= U.dk; domega -= U.domega;
  return *this;
}

inline HighTemp2D_cState& HighTemp2D_cState::operator *=(const double &a) {
  rho *= a; dv.x *= a; dv.y *= a; E *= a; dk *= a; domega *= a;
  return *this;
}

inline HighTemp2D_cState& HighTemp2D_cState::operator /=(const double &a) {
  rho /= a; dv.x /= a; dv.y /= a; E /= a; dk /= a; domega /= a;
  return *this;
}

/**********************************************************************
 * HighTemp2D_cState -- Relational operators.                     *
 **********************************************************************/
inline int operator ==(const HighTemp2D_cState &U1, const HighTemp2D_cState &U2) {
  return (U1.rho == U2.rho && U1.dv == U2.dv && U1.E == U2.E && U1.dk == U2.dk && U1.domega == U2.domega);
}

inline int operator !=(const HighTemp2D_cState &U1, const HighTemp2D_cState &U2) {
  return (U1.rho != U2.rho || U1.dv != U2.dv || U1.E != U2.E || U1.dk != U2.dk || U1.domega != U2.domega);
}

/**********************************************************************
 * HighTemp2D_cState -- Input-output operators.                   *
 **********************************************************************/
inline ostream &operator << (ostream &out_file, const HighTemp2D_cState &U) {
  out_file.setf(ios::scientific);
  out_file << " " << U.rho << " " << U.dv.x << " " << U.dv.y << " " << U.E << " " << U.dk << " " << U.domega;
  out_file.unsetf(ios::scientific);
  return out_file;
}

inline istream &operator >> (istream &in_file, HighTemp2D_cState &U) {
  in_file.setf(ios::skipws);
  in_file >> U.rho >> U.dv.x >> U.dv.y >> U.E >> U.dk >> U.domega;
  in_file.unsetf(ios::skipws);
  return in_file;
}

/**********************************************************************
 * Routine: Rotate                                                    *
 *                                                                    *
 * This function returns the solution in the local rotated frame.     *
 *                                                                    *
 **********************************************************************/
inline void HighTemp2D_cState::Rotate(const HighTemp2D_cState &U, const Vector2D &norm_dir) {
	rho = U.rho;
	dv.x =  U.dv.x*norm_dir.x+U.dv.y*norm_dir.y;
	dv.y = -U.dv.x*norm_dir.y+U.dv.y*norm_dir.x;
	E = U.E;
	dk = U.dk;
	domega = U.domega;
}

/**********************************************************************
 * HighTemp2D_pState::HighTemp2D_pState -- Constructor.       *
 **********************************************************************/
inline HighTemp2D_pState::HighTemp2D_pState(const HighTemp2D_cState &U) {
  rho = U.rho; v.x = U.v().x; v.y = U.v().y; p = U.p(); k = U.k(); omega = U.omega();
}

/**********************************************************************
 * HighTemp2D_pState::U -- Conserved solution state.              *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_pState::U(void) const {
  return U(*this);
}

inline HighTemp2D_cState HighTemp2D_pState::U(const HighTemp2D_pState &W) const {
  return HighTemp2D_cState(W.rho,W.dv(),W.E(),W.dk(),W.domega());
}

inline HighTemp2D_cState U(const HighTemp2D_pState &W) {
  return HighTemp2D_cState(W.rho,W.dv(),W.E(),W.dk(),W.domega());
}

/**********************************************************************
 * HighTemp2D_pState::dUdW -- Jacobian of the conserved solution  *
 *                                variables with respect to the       *
 *                                primitive solution variables .      *
 **********************************************************************/
inline void HighTemp2D_pState::dUdW(DenseMatrix &dUdW) const {

  //  cout<<"in dUdW, assigning variables"<<endl;
  switch (eos_type) {
  case EOS_IDEAL :
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
    break;
  case EOS_TGAS : { 
    double dpdex = dpde();
    double dpdrhox = dpdrho();
    double ex = e();
    dUdW(0,0) += ONE;
    dUdW(1,0) += v.x;
    dUdW(1,1) += rho;
    dUdW(2,0) += v.y;
    dUdW(2,2) += rho;
    dUdW(3,0) += ex + HALF*v.sqr() - rho*(dpdrhox/dpdex);
    dUdW(3,1) += rho*v.x;
    dUdW(3,2) += rho*v.y;
    dUdW(3,3) += rho/dpdex;
    if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      dUdW(3,0) += k;
      dUdW(3,4) += rho;
      dUdW(4,0) += k;
      dUdW(4,4) += rho;
      dUdW(5,0) += omega;
      dUdW(5,5) += rho;
    }
	}  break;
  }
}

/**********************************************************************
 * HighTemp2D_pState::dWdU -- Jacobian of the primitive solution  *
 *                                variables with respect to the       *
 *                                conserved solution variables.       *
 **********************************************************************/
inline void HighTemp2D_pState::dWdU(DenseMatrix &dWdU) const {

  //cout<<"in dWdU, assigning variables"<<endl;
  switch (eos_type) {
   case EOS_IDEAL :
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
     break;
   case EOS_TGAS : {
   double dpdex = dpde();
   double dpdrhox = dpdrho();
   double ex = e();
    dWdU(0,0) += ONE;
    dWdU(1,0) -= v.x/rho;
    dWdU(1,1) += ONE/rho;
    dWdU(2,0) -= v.y/rho;
    dWdU(2,2) += ONE/rho;
    dWdU(3,0) += dpdex/rho*HALF*(v.sqr() - 2.0*ex) + dpdrhox;
    dWdU(3,1) -= dpdex*v.x;
    dWdU(3,2) -= dpdex*v.y;
    dWdU(3,3) += dpdex/rho;
    if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      dWdU(3,4) -= dpdex/rho;
      dWdU(4,0) -= k/rho;
      dWdU(4,4) += ONE/rho;
      dWdU(5,0) -= omega/rho;
      dWdU(5,5) += ONE/rho;
    }
    } break;
  }
}

/**********************************************************************
 * HighTemp2D_pState::Fx -- Solution inviscid flux (x-direction).  *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_pState::Fx(void) const {
  return HighTemp2D_cState(rho*v.x,rho*sqr(v.x)+p+(2.0/3.0)*dk(),rho*v.x*v.y,v.x*H()+v.x*(2.0/3.0)*dk(),rho*v.x*k,rho*v.x*omega);
}

inline HighTemp2D_cState HighTemp2D_pState::Fx(const Vector2D &V) const {
  return HighTemp2D_cState(rho*(v.x-V.x),
			       rho*(v.x-V.x)*v.x+p+(2.0/3.0)*dk(),
			       rho*(v.x-V.x)*v.y,
			       (v.x-V.x)*E()+ v.x*(p+(2.0/3.0)*dk()),
			       rho*(v.x-V.x)*k,
			       rho*(v.x-V.x)*omega);
}
//NOT NEEDED!!!
/**********************************************************************
 * HighTemp2D_pState::Fx -- Solution inviscid flux (x-direction) HT only!! *
 **********************************************************************/
//inline HighTemp2D_cState HighTemp2D_pState::Fx(double dpde, double dpdrho) const {
//  return HighTemp2D_cState(rho*v.x,rho*sqr(v.x)+p+(2.0/3.0)*dk(),rho*v.x*v.y,v.x*H()+v.x*(2.0/3.0)*dk(),rho*v.x*k,rho*v.x*omega);
//}

/**********************************************************************
 * HighTemp2D_pState::dFdU -- Jacobian of the inviscid solution   *
 *                                flux with respect to the conserved  *
 *                                solution variables.                 *
 **********************************************************************/
inline void HighTemp2D_pState::dFdU(DenseMatrix &dFdU) const {

  switch (eos_type) {

   case EOS_IDEAL :
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
     }
   break;

   case EOS_TGAS : {
		double dhdpx = dhdp();
		double dhdrhox = dhdrho();
		double hx = Tgas_h(p, rho);
		
		dFdU(0,1) += 1.0;
		dFdU(1,0) += -(2.0*dhdpx*rho-3.0)/(dhdpx*rho-1.0)*v.x*v.x/2.0+1/(dhdpx*rho-1.0)*v.y*v.y/2.0-(hx+rho*dhdrhox)/(dhdpx*rho-1.0);
		dFdU(1,1) += (2.0*dhdpx*rho-3.0)/(dhdpx*rho-1.0)*v.x;
		dFdU(1,2) += -v.y/(dhdpx*rho-1.0);
		dFdU(1,3) += 1/(dhdpx*rho-1.0);
		dFdU(2,0) += -v.x*v.y;
		dFdU(2,1) += v.y;
		dFdU(2,2) += v.x;
		dFdU(3,0) += -(dhdpx*rho-2.0)/(dhdpx*rho-1.0)*v.x*v.x*v.x/2.0+(-(dhdpx*rho-2.0)/(dhdpx*rho-1.0)*v.y*v.y/2.0-rho*(dhdrhox+dhdpx*hx)/(dhdpx*rho-1.0))*v.x;
		dFdU(3,1) += (dhdpx*rho-3.0)/(dhdpx*rho-1.0)*v.x*v.x/2.0+hx+v.y*v.y/2.0;
		dFdU(3,2) += -v.y/(dhdpx*rho-1.0)*v.x;
		dFdU(3,3) += rho*v.x*dhdpx/(dhdpx*rho-1.0);

// Dagmara's take one:
//       dFdU(0,1) += ONE;
//       dFdU(1,0) += dpde()*(sqr(v.x) + sqr(v.y) - 2.0*e())/(2.0*rho) + dpdrho() - sqr(v.x);     
//       dFdU(1,1) -= 2.0*v.x - v.x*dpde()/rho;
//       dFdU(1,2) -= v.y*dpde()/rho;
//       dFdU(1,3) += dpde()/rho;
//       dFdU(2,0) -= v.x*v.y;
//       dFdU(2,1) += v.y;
//       dFdU(2,2) += v.x;
//       dFdU(3,0) -= v.x*(h() - dpdrho() - (dpde()*HALF/rho)*(sqr(v.x)+sqr(v.y)-2.0*e())+ (2.0/3.0)*k);
//       dFdU(3,1) +=  h() - dpde()*sqr(v.x)/rho + (2.0/3.0)*k;
//       dFdU(3,2) -= v.x*v.y*dpde()/rho;
//       dFdU(3,3) += v.x*dpde()/rho + ONE;
     if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
     // May not work -- Alistair Wood
       dFdU(1,4) += 2.0/3.0 - dpde()/rho; 
       dFdU(3,4) -= v.x*(2.0/3.0 - dpde()/rho);
       dFdU(4,0) -= v.x*k;
       dFdU(4,1) += k;
       dFdU(4,4) += v.x;
       dFdU(5,0) -= v.x*omega;
       dFdU(5,1) += omega;
       dFdU(5,5) += v.x;
     }
  } break;
  }
  
}

#if defined(ALI_CHECK_HIGHTEMP) && defined(COMPILING_DRDU)
bool ali_dump_diffs_global = true;
int ali_dump_diffs_cpu = 0;
#endif

 /**********************************************************************
 * HighTemp2D_pState::dFdW -- Jacobian of the inviscid solution   *
 *                                flux with respect to the primitive  *
 *                                solution variables.                 *
 **********************************************************************/
inline void HighTemp2D_pState::dFdW(DenseMatrix &dFdW) const
{
	//  The mass and momentum inviscid fluxes written in terms of the
	//  primitive variables do not require a gas equation of state. As
	//  such, the first three rows of dFdW do not depend on the
	//  equation of state.
	//  
	//  This is in contrast to dFdU where an equation of state is needed
	//  to relate, for example, pressure in the momentum flux to energy.

	dFdW(0,0) = v.x;
	dFdW(0,1) = rho;
	dFdW(1,0) = v.x*v.x;
	dFdW(1,1) = TWO*rho*v.x; 
	dFdW(1,3) = ONE;
	dFdW(2,0) = v.x*v.y;
	dFdW(2,1) = rho*v.y;
	dFdW(2,2) = rho*v.x;
	switch (eos_type) {
		case EOS_IDEAL :
			dFdW(3,0) = HALF*(v.x*v.x+v.y*v.y)*v.x;
			dFdW(3,1) = rho*v.x*v.x+ rho*(g*p/rho/gm1 + HALF*(v.x*v.x+v.y*v.y));
			dFdW(3,2) = rho*v.x*v.y;
			dFdW(3,3) = v.x*g/gm1;
			break;
 		case EOS_TGAS : {

			//  The function h() returns the total enthalpy (i.e.
			//  including kinetic energy). Here I want internal
			//  enthalpy. So I call Tgas_h directly. It feels like
			//  there is a better way to do this. You think on it.
			//   -- Alistair Wood Mon Apr 09 2007 

			double dhdpx = dhdp();
			double dhdrhox = dhdrho();
			double hx = Tgas_h(p, rho);
			dFdW(3,0) = v.x*v.x*v.x/2.0+(hx+v.y*v.y/2.0+rho*dhdrhox)*v.x;
			dFdW(3,1) = 3.0/2.0*rho*v.x*v.x+rho*v.y*v.y/2.0+rho*hx;
			dFdW(3,2) = rho*v.x*v.y;
			dFdW(3,3) = rho*v.x*dhdpx;
			} break;
	}
	if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
		dFdW(1,0) += TWO/THREE*k;
		dFdW(1,4) = TWO/THREE*rho;
		dFdW(3,0) += FIVE/THREE*k*v.x;
		dFdW(3,1) += FIVE/THREE*rho*k; 
		dFdW(3,4) += FIVE/THREE*rho*v.x;
		dFdW(4,0) = v.x*k;
		dFdW(4,1) = rho*k;  
		dFdW(5,0) = v.x*omega;
		dFdW(5,1) = rho*omega;
	}  
	dFdW(4,4) = rho*v.x;
	dFdW(5,5) = rho*v.x;  

	// This function ends here. The rest is debugging for my sanity.
	//   -- Alistair Wood Mon Apr 09 2007 

//  #if defined(ALI_CHECK_HIGHTEMP) && defined(COMPILING_DRDU)
//  	if (ali_dump_diffs_global) {
//  
//  		char fname[100];
//  		sprintf(fname, "inviscid_dump_cpu%.3d.txt", ali_dump_diffs_cpu);
//  		ofstream fout;
//  		fout.open(fname, ios::app);
//  
//  		double ival[4] = {0.0}, hval[4] = {0.0};
//  
//  		ival[0] = HALF*(v.x*v.x+v.y*v.y)*v.x;
//  		ival[1] = rho*v.x*v.x+ rho*(g*p/rho/gm1 + HALF*(v.x*v.x+v.y*v.y));
//  		ival[2] = rho*v.x*v.y;
//  		ival[3] = v.x*g/gm1;
//  
//  		hval[0] = v.x*v.x*v.x/2.0+(Tgas_h(p, rho)+v.y*v.y/2.0+rho*dhdrho())*v.x;
//  		hval[1] = 3.0/2.0*rho*v.x*v.x+rho*v.y*v.y/2.0+rho*Tgas_h(p, rho);
//  		hval[2] = rho*v.x*v.y;
//  		hval[3] = rho*v.x*dhdp();
//  
//  #define QQTOL 0.000001
//  		fout << setprecision(14) << scientific;
//  		fout << "p = " << p << endl;
//  		fout << "rho = " << rho << endl;
//  		fout << "g = " << g << endl;
//  		fout << "gm1 = " << gm1 << endl;
//  		fout << "g p / rho / gm1 = " << g * p / rho / gm1 << endl;
//  		fout << "Tgas_h = " << Tgas_h(p, rho) << endl;
//  
//  		// this should be small:
//  		fout << "((g p / rho / gm1) - Tgas_h) / Tgas_h = " << (g*p/rho/gm1 - Tgas_h(p, rho))/Tgas_h(p, rho) << endl << endl;
//  
//  		fout << "- g p / rho / rho / gm1 = " << - g * p / rho / rho / gm1 << endl;
//  		fout << "dhdrho() = " << dhdrho() << endl;
//  
//  		// this should be small too:
//  		fout << "((- g p / rho / rho / gm1) - dhdrho()) / dhdrho() = " << ((- g * p / rho / rho / gm1) - dhdrho()) / dhdrho() << endl << endl;
//  
//  		fout << "rho * dhdrho() = " << rho * dhdrho() << endl;
//  
//  		// in the ideal case this is identically zero. 
//  		fout << "Tgas_h + rho * dhdrho() = " << Tgas_h(p, rho) + rho * dhdrho() << "  ideally zero ..." << endl;
//  
//  		// just for my sanity:
//  		fout << "\nh / (rho - delta^2) = " << Tgas_h(p, rho) / (rho - sqr(HTTOL)) << "  for my sanity\n\n";
//  
//  		// this is what skews entry (3, 0) so much. 
//  		fout << "u (Tgas_h + rho * dhdrho())  = " << v.x*(Tgas_h(p, rho) + rho * dhdrho()) << "  what skews (3,0)" << endl;
//  
//  		fout << "u = " << v.x << endl;
//  		fout << "v = " << v.y << endl;
//  		fout << "\n";
//  		fout << "entry    ideal       hightemp    diff normalized by ideal value\n";
//  		for (int qi = 0; qi < 4; qi++) {
//  			if (qi == 2) { continue; }
//  			fout << setw(4) << qi;
//  			fout << fixed << setprecision(3) << setw(13) << ival[qi];
//  			fout << fixed << setprecision(3) << setw(13) << hval[qi];
//  			fout << scientific << setprecision(2) << setw(12) << ( (hval[qi] - ival[qi]) / max(fabs(ival[qi]), QQTOL) );
//  			fout << endl;
//  		}
//  		fout << "\n--\n\n";
//  		fout.close();
//  	} // if (ali_dump_diffs_global)
//  
//  #endif // ALI_CHECK_HIGHTEMP 

}

/**********************************************************************
 * HighTemp2D_pState::Gx, Gy -- Solution viscous fluxes.          *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_pState::Gx(const HighTemp2D_pState &dWdx) const {
  return HighTemp2D_cState(ZERO,tau.xx,tau.xy,-q.x+v.x*tau.xx+v.y*tau.xy+(mu()+sigma_k*muT())*k,(mu()+sigma_k*muT())*k,(mu()+sigma_omega*muT())*omega);
}

inline HighTemp2D_cState HighTemp2D_pState::Gy(const HighTemp2D_pState &dWdy) const {
  return HighTemp2D_cState(ZERO,tau.xy,tau.yy,-q.y+v.x*tau.xy+v.y*tau.yy+(mu()+sigma_k*muT())*k,(mu()+sigma_k*muT())*k,(mu()+sigma_omega*muT())*omega);
}

/**********************************************************************
 * HighTemp2D_pState::ComputeViscousTerms -- Compute viscous      *
 *                                               stress tensor and    *
 *                                               heat flux vector.    *
 **********************************************************************/
inline void HighTemp2D_pState::ComputeViscousTerms(const HighTemp2D_pState &dWdx,
						       const HighTemp2D_pState &dWdy,
						       const Vector2D &X,
						       int Axisymmetric)
{
  double radius;
  if (Axisymmetric) radius = max(X.y,TOLER);
  // Divergence of the velocity field.
  double div = dWdx.v.x + dWdy.v.y;
  if (Axisymmetric) div += v.y/radius;

  // Stress tensor.
  double mumu = mu() + muT(); 
  tau.xx = 2.0*mumu*(dWdx.v.x - div/3.0);
  tau.xy = mumu*(dWdy.v.x + dWdx.v.y);
  tau.yy = 2.0*mumu*(dWdy.v.y - div/3.0);
  if (Axisymmetric) tau.zz = 2.0*mumu*(v.y/radius - div/3.0);
  else tau.zz = ZERO;

	//  Heat flux components:
	//
	//  Previously, an "adiabatic_flag" was passed to this function; if
	//  it was set then we did:
	//  
	//    q = Vector2D_ZERO;
	//  
	//  But Jai Sachdev realised that we only want the normal (to the
	//  wall) component to be zero. To do this we had two choices:
	//    1. keep the adiabatic flag and do something like: rotate, set
	//       the x-component to zero, rotate back, or
	//    2. remove the adiabatic flag and trust that the application of
	//       the boudary conditions would do the right thing.
	//  The second option was selected.
	//      -- Alistair Wood Thu Mar 15 2007 
	double kap = kappa() + kappaT();
	switch (eos_type) {
		case EOS_TGAS: {
			double dTdp_local = dTdp();
			double dTdrho_local = dTdrho();
			q.x = -kap*(dWdx.p*dTdp_local + dWdx.rho*dTdrho_local);
			q.y = -kap*(dWdy.p*dTdp_local + dWdy.rho*dTdrho_local);
			} break;
		case EOS_IDEAL:
			q.x = -kap*(dWdx.p - (p/rho)*dWdx.rho)/(rho*R);
			q.y = -kap*(dWdy.p - (p/rho)*dWdy.rho)/(rho*R);
			break;
	}
}

/**********************************************************************
 * HighTemp2D_pState::lambda_x -- Eigenvalue(s) (x-direction).    *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_pState::lambda_x(void) const {
  return HighTemp2D_pState(v.x-c(),v.x,v.x,v.x+c(),v.x,v.x);
}

inline HighTemp2D_pState HighTemp2D_pState::lambda_x(const Vector2D &V) const {
  return HighTemp2D_pState(v.x-V.x-c(),v.x-V.x,v.x-V.x,v.x-V.x+c(),v.x-V.x,v.x-V.x);
}

/**********************************************************************
 * HighTemp2D_pState::lambda_x -- Eigenvalue(s) (x-direction).    *
 * For High-Temperature Air, average State only                   *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_pState::lambda_x(double aAvg) const {
  return HighTemp2D_pState(v.x-aAvg,v.x,v.x,v.x+aAvg,v.x,v.x);
}

/**********************************************************************
 * HighTemp2D_pState::rp_x -- Primitive right eigenvector         *
 *                                (x-direction).                      *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_pState::rp_x(int index) const {
  //assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
  switch(index) {
  case 1 :
    return HighTemp2D_pState(ONE,-c()/rho,ZERO,c2()-(2.0/3.0)*k,ZERO,ZERO);
    break;
  case 2 :
    return HighTemp2D_pState(ONE,ZERO,ZERO,-(2.0/3.0)*k,ZERO,ZERO);
    break;
  case 3 :
    return HighTemp2D_pState(ZERO,ZERO,ONE,ZERO,ZERO,ZERO);
    break;
  case 4 :
    return HighTemp2D_pState(ONE,c()/rho,ZERO,c2()-(2.0/3.0)*k,ZERO,ZERO);
    break;
  case 5 :
    return HighTemp2D_pState(ZERO,ZERO,ZERO,-(2.0/3.0)*rho,ONE,ZERO);
    break;
  case 6 :
    return HighTemp2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE);
    break;
  };
  return HighTemp2D_pState(ONE,-c()/rho,ZERO,c2()-(2.0/3.0)*k,ZERO,ZERO);
}

/**********************************************************************
 * HighTemp2D_pState::rc_x -- Conserved right eigenvector         *
 *                                (x-direction).                      *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_pState::rc_x(int index) const {
  //assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
  switch(index) {
  case 1 :
    return HighTemp2D_cState(ONE,v.x-c(),v.y,h()-c()*v.x+(2.0/3.0)*k,k,omega);
    break;
  case 2 :
    //switch(eos_type){
    //case EOS_TGAS:
    if (eos_type == EOS_TGAS)
      return HighTemp2D_cState(ONE,v.x,v.y,h()+2.0*k/3.0-rho*c()*c()/dpde(),k,omega); 
    else
      return HighTemp2D_cState(ONE,v.x,v.y,HALF*v.sqr()+gm1i*(g-5.0/3.0)*k,k,omega);
    break;
  case 3 :
    return HighTemp2D_cState(ZERO,ZERO,rho,rho*v.y,ZERO,ZERO);
    break;
  case 4 :
    return HighTemp2D_cState(ONE,v.x+c(),v.y,h()+c()*v.x+(2.0/3.0)*k,k,omega);
    break;
  case 5 :
    if (eos_type == EOS_TGAS)
      return HighTemp2D_cState(ZERO,ZERO,ZERO,rho*(-2.0*rho/(3.0*dpde())+ONE),rho,ZERO); 
    else 
      return HighTemp2D_cState(ZERO,ZERO,ZERO,rho*gm1i*(g-5.0/3.0),rho,ZERO);
    break;
  case 6 :
    return HighTemp2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,rho);
    break;
  };
  return HighTemp2D_cState(ONE,v.x-c(),v.y,h()-c()*v.x+(2.0/3.0)*k,k,omega);
}

/**********************************************************************
 * HighTemp2D_pState::rc_x -- Conserved right eigenvector         *
 *                                (x-direction).                      *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_pState::rc_x(int index,double dpde,double dpdrho,double cAvg) const {
  //assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
  switch(index) {
  case 1 :
    //return HighTemp2D_cState(ONE,v.x-cAvg,v.y,HALF*v.sqr()+k-cAvg*v.x+e()-rho*(2.0*k+3.0*cAvg*cAvg-3.0*dpdrho)/(3.0*dpde),k,omega);
    return HighTemp2D_cState(ONE,v.x-cAvg,v.y,h()-cAvg*v.x+(2.0/3.0)*k,k,omega);
    break;
  case 2 :
    //return HighTemp2D_cState(ONE,v.x,ZERO,HALF*(v.x*v.x-v.y*v.y)+e()-rho*dpdrho/dpde,ZERO,ZERO);
    return HighTemp2D_cState(ONE,v.x,v.y,h()+2.0*k/3.0-rho*cAvg*cAvg/dpde,k,omega);
    break;
  case 3 :
    return HighTemp2D_cState(ZERO,ZERO,rho,rho*v.y,ZERO,ZERO);
    break;
  case 4 :
    //return HighTemp2D_cState(ONE,v.x+cAvg,v.y,HALF*v.sqr()+k+cAvg*v.x+e()-rho*(2.0*k+3.0*cAvg*cAvg-3.0*dpdrho)/(3.0*dpde),k,omega);
    return HighTemp2D_cState(ONE,v.x+cAvg,v.y,h()+cAvg*v.x+(2.0/3.0)*k,k,omega);
    break;
  case 5 :
    return HighTemp2D_cState(ZERO,ZERO,ZERO,rho*(-2.0*rho/(3.0*dpde)+ONE),rho,ZERO);
    break;
  case 6 :
    return HighTemp2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,rho);
    break;
  };
  return HighTemp2D_cState(ZERO,ZERO,rho,rho*v.y,ZERO,ZERO);
}

/**********************************************************************
 * HighTemp2D_pState::lp_x -- Primitive left eigenvector          *
 *                                (x-direction).                      *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_pState::lp_x(int index) const {
  //assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
  switch(index) {
  case 1 :
    return HighTemp2D_pState(k/(3.0*c2()),-HALF*rho/c(),ZERO,HALF/c2(),rho/(3.0*c2()),ZERO);
    break;
  case 2 :
    return HighTemp2D_pState(ONE-(2.0/3.0)*k/c2(),ZERO,ZERO,-ONE/c2(),-(2.0/3.0)*rho/c2(),ZERO);
    break;
  case 3 :
    return HighTemp2D_pState(ZERO,ZERO,ONE,ZERO,ZERO,ZERO);
    break;
  case 4 :
    return HighTemp2D_pState(k/(3.0*c2()),HALF*rho/c(),ZERO,HALF/c2(),rho/(3.0*c2()),ZERO);
    break;
  case 5 :
    return HighTemp2D_pState(ZERO,ZERO,ZERO,ZERO,ONE,ZERO);
    break;
  case 6 :
    return HighTemp2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE);
    break;
  };
  return HighTemp2D_pState(k/(3.0*c2()),-HALF*rho/c(),ZERO,HALF/c2(),rho/(3.0*c2()),ZERO);
}

/**********************************************************************
 * HighTemp2D_pState::lp_x -- Primitive left eigenvector          *
 *  ONLY 4 GLAISTER FLUX      (x-direction).                      *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_pState::lp_x(int index, double cAvg) const {
  //assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
  
  switch(index) {
  case 1 :
    return HighTemp2D_pState(k/(3.0*cAvg*cAvg),-HALF*rho/cAvg,ZERO,HALF/(cAvg*cAvg),rho/(3.0*cAvg*cAvg),ZERO);
    break;
  case 2 :
    return HighTemp2D_pState(ONE-(2.0/3.0)*k/(cAvg*cAvg),ZERO,ZERO,-ONE/(cAvg*cAvg),-(2.0/3.0)*rho/(cAvg*cAvg),ZERO);
    break;
  case 3 :
    return HighTemp2D_pState(ZERO,ZERO,ONE,ZERO,ZERO,ZERO);
    break;
  case 4 :
    return HighTemp2D_pState(k/(3.0*cAvg*cAvg),HALF*rho/cAvg,ZERO,HALF/(cAvg*cAvg),rho/(3.0*cAvg*cAvg),ZERO);
    break;
  case 5 :
    return HighTemp2D_pState(ZERO,ZERO,ZERO,ZERO,ONE,ZERO);
    break;
  case 6 :
    return HighTemp2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE);
    break;
  };
  return HighTemp2D_pState(k/(3.0*cAvg*cAvg),-HALF*rho/cAvg,ZERO,HALF/(cAvg*cAvg),rho/(3.0*cAvg*cAvg),ZERO);
  
}

/**********************************************************************
 * HighTemp2D_pState::S -- Include all source term vectors and    *
 *                             Jacobians.                             *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_pState::S(const Vector2D &X,
						      const HighTemp2D_pState &dWdx,
						      const HighTemp2D_pState &dWdy,
						      const int &Axisymmetric) const {
  HighTemp2D_cState Sall; Sall.Vacuum();
  // Include the axisymmetric source terms if required.
  if (Axisymmetric) {
    Sall = Si(X);
    if (flow_type) Sall += Sv(X,dWdy);
  }
  // Include the turbulence model source term if required.
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) Sall += St(X,dWdx,dWdy,Axisymmetric);
  // Return the total source term vector.
  return Sall; 
}

inline void HighTemp2D_pState::dSdU(DenseMatrix &dSdU,
					const Vector2D &X,
					const HighTemp2D_pState &dWdx,
					const HighTemp2D_pState &dWdy,
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
 * HighTemp2D_pState::Si -- Inviscid axisymmetric source terms    *
 *                              and Jacobian.                         *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_pState::Si(const Vector2D &X) const {
  return HighTemp2D_cState(-rho*v.y/X.y,-rho*v.x*v.y/X.y,-rho*sqr(v.y)/X.y,-v.y*(H()+(2.0/3.0)*dk())/X.y,-v.y*dk()/X.y,-v.y*domega()/X.y);
}

inline void HighTemp2D_pState::dSidU(DenseMatrix &dSidU, const Vector2D &X) const {
  //cout<<"Si function uses gm1!! in terms for dSidU"<<endl;
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
  }
}

/**********************************************************************
 * HighTemp2D_pState::Sv -- Viscous axisymmetric flow source term *
 *                              vector and Jacobian.                  *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_pState::Sv(const Vector2D &X,
						       const HighTemp2D_pState &dWdy) const {
  return HighTemp2D_cState(ZERO,
			       tau.xy/X.y,
			       (tau.yy-tau.zz)/X.y,
			       (-q.y+v.x*tau.xy+v.y*tau.yy+(mu()+sigma_k*muT())*dWdy.k)/X.y,
			       (mu()+sigma_k*muT())*dWdy.k/X.y,
			       (mu()+sigma_omega*muT())*dWdy.omega/X.y);
}

inline void HighTemp2D_pState::dSvdU(DenseMatrix &dSvdU,
					 const Vector2D &X,
					 const HighTemp2D_pState &dWdy) const {
  dSvdU(0,0) += ZERO;
  dSvdU(1,1) -= ZERO;
}

/**********************************************************************
 * HighTemp2D_pState::St -- Turbulent source term vector and      *
 *                              Jacobian.                             *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_pState::St(const Vector2D &X,
						       const HighTemp2D_pState &dWdx,
						       const HighTemp2D_pState &dWdy,
						       const int &Axisymmetric) const {
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
  return HighTemp2D_cState(ZERO,ZERO,ZERO,ZERO,
			       production-beta_k(dWdx,dWdy)*dk()*omega,
			       alpha*(omega/max(k,TOLER))*production-beta_omega(dWdx,dWdy)*rho*sqr(omega));
}

inline void HighTemp2D_pState::dStdU(DenseMatrix &dStdU,
					 const Vector2D &X,
					 const HighTemp2D_pState &dWdx,
					 const HighTemp2D_pState &dWdy,
					 const int &Axisymmetric) const {
  dStdU(0,0) += ZERO;
  dStdU(1,1) -= ZERO;
}

/**********************************************************************
 * Routine: Rotate                                                    *
 *                                                                    *
 * This function returns the solution in the local rotated frame.     *
 *                                                                    *
 **********************************************************************/
inline void HighTemp2D_pState::Rotate(const HighTemp2D_pState &W, const Vector2D &norm_dir) {
	rho = W.rho;
	v.x =  W.v.x*norm_dir.x+W.v.y*norm_dir.y;
	v.y = -W.v.x*norm_dir.y+W.v.y*norm_dir.x;
	p = W.p;
	k = W.k;
	omega = W.omega;
}

// Should only be called from the preconditioner. AW. Tue Aug 07 2007.
inline HighTemp2D_pState Rotate(const HighTemp2D_pState &W,
		const Vector2D &norm_dir) 
{
	HighTemp2D_pState W_rotated;
	W_rotated.Rotate(W, norm_dir);
	return W_rotated;
}

/**********************************************************************
 * HighTemp2D_cState::HighTemp2D_cState -- Constructor.       *
 **********************************************************************/
inline HighTemp2D_cState::HighTemp2D_cState(const HighTemp2D_pState &W) {
  rho = W.rho; dv = W.dv(); E = W.E(); dk = W.dk(); domega = W.domega();
}

/**********************************************************************
 * HighTemp2D_cState::W -- Primitive solution state.              *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_cState::W(void) const {
  return W(*this);
}

inline HighTemp2D_pState HighTemp2D_cState::W(const HighTemp2D_cState &U) const {
  return HighTemp2D_pState(U.rho,U.v().x,U.v().y,U.p(),U.k(),U.omega());
}

inline HighTemp2D_pState W(const HighTemp2D_cState &U) {
  return HighTemp2D_pState(U.rho,U.v().x,U.v().y,U.p(),U.k(),U.omega());
}

/**********************************************************************
 * HighTemp2D_cState::dUdW -- Jacobian of the conserved solution  *
 *                                variables with respect to the       *
 *                                primitive solution variables.       *
 **********************************************************************/
inline void HighTemp2D_cState::dUdW(DenseMatrix &dUdW) const {
  //cout<<"in dUdW, assigning variables"<<endl;
  switch (eos_type) {
  case EOS_IDEAL :
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
    }
    break;
  case EOS_TGAS : {
  double dpdex = dpde();
  double dpdrhox = dpdrho();
    dUdW(0,0) += ONE;
    dUdW(1,0) += v().x;
    dUdW(1,1) += rho;
    dUdW(2,0) += v().y;
    dUdW(2,2) += rho;
    dUdW(3,0) += e() + HALF*v().sqr() - rho*(dpdrhox/dpdex);
    dUdW(3,1) += rho*v().x;
    dUdW(3,2) += rho*v().y;
    dUdW(3,3) += rho/dpdex;
    if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      dUdW(3,0) += k();
      dUdW(3,4) += rho;
      dUdW(4,0) += k();
      dUdW(4,4) += rho;
      dUdW(5,0) += omega();
      dUdW(5,5) += rho;
    }
    } break;
  }
}

/**********************************************************************
 * HighTemp2D_cState::dWdU -- Jacobian of the primitive solution  *
 *                         variables with respect to the conserved    *
 *                         solution variables.                        *
 **********************************************************************/
inline void HighTemp2D_cState::dWdU(DenseMatrix &dWdU) const {
    //cout<<"in dWdU, assigning variables"<<endl;
  switch (eos_type) {
  case EOS_IDEAL :
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
      dWdU(4,0) -= k()/rho;
      dWdU(4,4) += ONE/rho;
      dWdU(5,0) -= omega()/rho;
      dWdU(5,5) += ONE/rho;
    }
    break;
  case EOS_TGAS : {
  double dpdex = dpde();
  double dpdrhox = dpdrho();
    dWdU(0,0) += ONE;
    dWdU(1,0) -= dv.x/(rho*rho);
    dWdU(1,1) += ONE/rho;
    dWdU(2,0) -= dv.y/(rho*rho);
    dWdU(2,2) += ONE/rho;
    dWdU(3,0) += dpdex/rho*HALF*(dv.sqr()/rho - 2.0*e()) + dpdrhox;
    dWdU(3,1) -= dpdex*dv.x/rho;
    dWdU(3,2) -= dpdex*dv.y/rho;
    dWdU(3,3) += dpdex/rho;
    if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
      dWdU(3,4) -= dpdex/rho;
      dWdU(4,0) -= k()/rho;
      dWdU(4,4) += ONE/rho;
      dWdU(5,0) -= omega()/rho;
      dWdU(5,5) += ONE/rho;
    }
    } break;
  }
}

/**********************************************************************
 * HighTemp2D_cState::F -- Solution inviscid flux (x-direction).  *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_cState::Fx(void) const {
  return HighTemp2D_cState(dv.x,sqr(dv.x)/rho+p()+(2.0/3.0)*dk,dv.x*dv.y/rho,dv.x*H()/rho+v().x*(2.0/3.0)*dk,v().x*dk,v().x*domega);
}

inline HighTemp2D_cState HighTemp2D_cState::Fx(const Vector2D &V) const {
  double vx = v().x;
  return HighTemp2D_cState(rho*(vx-V.x),(vx-V.x)*dv.x+p()+(2.0/3.0)*dk,
			       (vx-V.x)*dv.y,(vx-V.x)*E+ vx*(p()+(2.0/3.0)*dk),
			       (vx-V.x)*dk,(vx-V.x)*domega);
}

//  Because I was experiencing problems with the preconditioner in
//  the Newton-Krylov-Schwarz code, I rewrote the high-temperature
//  part of the primitive dFdU function. The problems with the
//  preconditioner went away. I was not paying close enough
//  attention to conclusively blame the first version.
//  
//  In this rewrite I used enthalpy instead of internal energy.
//  This meant that I called functions such as dhdp instead of
//  functions such as dpde. For the primitive state this is good
//  since a curve fit of enthalpy as a function of pressure and
//  density is available (and this curve fit does not use a "damped
//  Newton's method iteration scheme" unlike the curve fit which
//  returns internal energy as a function of pressure and
//  temperature).
//  
//  Given that I had problems with the primitive dFdU function I
//  would suspect that this conservative dFdU function could be
//  broken as well. However, I would not like to simply copy my
//  rewrite from the primitive function since this would be
//  inefficient. For example, to calculate dhdp we would first need
//  to calculate pressure but pressure is not a member of the
//  conservative solution vector. However, dpde would be efficient
//  since internal energy is a member of the conservative solution
//  vector (at least it should be - why does HighTemp2D_cState::e()
//  call curve fits? Isn't energy right there?). Also, there is a
//  curve fit (which does not use Newton's method) that returns
//  pressure given internal energy and density (and density is also
//  a member of the conservative solution vector).
//  
//  But no one actually calls the conservative dFdU function.
//  So why not just eliminate it?
//
//  Alistair Wood
//  Sunday, April 15, 2007
//
//  /**********************************************************************
//   * HighTemp2D_cState::dFdU -- Jacobian of the inviscid solution   *
//   *                                flux with respect to the conserved  *
//   *                                solution variables.                 *
//   **********************************************************************/
//  inline void HighTemp2D_cState::xdFdU(DenseMatrix &dFdU) const {
//    //cout<<"in cState dFdU gm1 used"<<endl;
//    //cout<<"in dFdU conservative, assigning variables"<<endl;
//    switch(eos_type) {
//     case EOS_IDEAL :
//       dFdU(0,1) += ONE;
//       dFdU(1,0) += HALF*gm1*(sqr(v().x) + sqr(v().y)) - sqr(v().x);
//       dFdU(1,1) -= v().x*(g-3.0);
//       dFdU(1,2) -= v().y*gm1;
//       dFdU(1,3) += gm1;
//       dFdU(2,0) -= v().x*v().y;
//       dFdU(2,1) += v().y;
//       dFdU(2,2) += v().x;
//       dFdU(3,0) -= (h() - HALF*gm1*(sqr(v().x)+sqr(v().y)))*v().x + (2.0/3.0)*dv.x*k();
//       dFdU(3,1) +=  h() - gm1*sqr(v().x) + (2.0/3.0)*k();
//       dFdU(3,2) -= v().x*v().y*gm1;
//       dFdU(3,3) += v().x*g;
//       if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
//         dFdU(1,4) += 5.0/3.0 - g;
//         dFdU(3,4) -= (g - 5.0/3.0)*v().x;
//         dFdU(4,0) -= v().x*k();
//         dFdU(4,1) += k();
//         dFdU(4,4) += v().x;
//         dFdU(5,0) -= v().x*omega();
//         dFdU(5,1) += omega();
//         dFdU(5,5) += v().x;
//       }
//     break;
//     case EOS_TGAS :
//       dFdU(0,1) += ONE;
//       dFdU(1,0) += dpde()*(sqr(v().x) + sqr(v().y) - 2.0*e())/(2.0*rho) + dpdrho() - sqr(v().x);     
//       dFdU(1,1) -= 2.0*v().x - v().x*dpde()/rho;
//       dFdU(1,2) -= v().y*dpde()/rho;
//       dFdU(1,3) += dpde()/rho;
//       dFdU(2,0) -= v().x*v().y;
//       dFdU(2,1) += v().y;
//       dFdU(2,2) += v().x;
//       dFdU(3,0) -= v().x*(h() - dpdrho() - (dpde()*HALF/rho)*(sqr(v().x)+sqr(v().y)-2.0*e())+ (2.0/3.0)*k());
//       dFdU(3,1) +=  h() - dpde()*sqr(v().x)/rho + (2.0/3.0)*k();
//       dFdU(3,2) -= v().x*v().y*dpde()/rho;
//       dFdU(3,3) += v().x*dpde()/rho + ONE;
//       if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) {
//         dFdU(1,4) += 2.0/3.0 - dpde()/rho;
//         dFdU(3,4) -= v().x*(2.0/3.0 - dpde()/rho);
//         dFdU(4,0) -= v().x*k();
//         dFdU(4,1) += k();
//         dFdU(4,4) += v().x;
//         dFdU(5,0) -= v().x*omega();
//         dFdU(5,1) += omega();
//         dFdU(5,5) += v().x;
//       }
//    }
//  }

/**********************************************************************
 * HighTemp2D_cState::Gx, Gy -- Solution viscous fluxes.          *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_cState::Gx(const HighTemp2D_pState &dWdx) const {
  return HighTemp2D_cState(ZERO,tau.xx,tau.xy,-q.x+v().x*tau.xx+v().y*tau.xy+(mu()+sigma_k*muT())*dWdx.k,(mu()+sigma_k*muT())*dWdx.k,(mu()+sigma_omega*muT())*dWdx.omega);
}

inline HighTemp2D_cState HighTemp2D_cState::Gy(const HighTemp2D_pState &dWdy) const {
  return HighTemp2D_cState(ZERO,tau.xy,tau.yy,-q.y+v().x*tau.xy+v().y*tau.yy+(mu()+sigma_k*muT())*dWdy.k,(mu()+sigma_k*muT())*dWdy.k,(mu()+sigma_omega*muT())*dWdy.omega);
}

/**********************************************************************
 * HighTemp2D_cState::ComputeViscousTerms -- Compute viscous      *
 *                                               stress tensor and    *
 *                                               heat flux vector.    *
 **********************************************************************/
inline void HighTemp2D_cState::ComputeViscousTerms(const HighTemp2D_pState &dWdx,
						       const HighTemp2D_pState &dWdy,
						       const Vector2D &X,
						       int Axisymmetric) 
{
  double radius;
  if (Axisymmetric) radius = max(X.y,TOLER);
  // Divergence of the velocity field.
  double div = dWdx.v.x + dWdy.v.y;
  if (Axisymmetric) div += v().y/radius;

  // Stress tensor.
  double mumu = mu() + muT();
  tau.xx = 2.0*mumu*(dWdx.v.x - div/3.0);
  tau.xy = mumu*(dWdy.v.x + dWdx.v.y);
  tau.yy = 2.0*mumu*(dWdy.v.y - div/3.0);
  if (Axisymmetric) tau.zz = 2.0*mumu*(v().y/radius - div/3.0);
  else tau.zz = ZERO;

	//  Heat flux components:
	//
	//  Previously, an "adiabatic_flag" was passed to this function; if
	//  it was set then we did:
	//  
	//    q = Vector2D_ZERO;
	//  
	//  But Jai Sachdev realised that we only want the normal (to the
	//  wall) component to be zero. To do this we had two choices:
	//    1. keep the adiabatic flag and do something like: rotate, set
	//       the x-component to zero, rotate back, or
	//    2. remove the adiabatic flag and trust that the application of
	//       the boudary conditions would do the right thing.
	//  The second option was selected.
	//      -- Alistair Wood Thu Mar 15 2007 
	double kap = kappa() + kappaT();
	switch (eos_type) {
		case EOS_TGAS: {
			double dTdp_local = dTdp();
			double dTdrho_local = dTdrho();
			q.x = -kap*(dWdx.p*dTdp_local + dWdx.rho*dTdrho_local);
			q.y = -kap*(dWdy.p*dTdp_local + dWdy.rho*dTdrho_local);
			} break;
		case EOS_IDEAL:
			q.x = -kap*(dWdx.p - (p()/rho)*dWdx.rho)/(rho*R);
			q.y = -kap*(dWdy.p - (p()/rho)*dWdy.rho)/(rho*R);
			break;
	}
}

/**********************************************************************
 * HighTemp2D_cState::lambda_x -- Eigenvalue(s) (x-direction).    *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_cState::lambda_x(void) const {
  double vx = v().x, cc = c();
  return HighTemp2D_pState(vx-cc,vx,vx,vx+cc,vx,vx);
}

inline HighTemp2D_pState HighTemp2D_cState::lambda_x(const Vector2D &V) const {
  double vx = v().x, cc = c();
  return HighTemp2D_pState(vx-V.x-cc,vx-V.x,vx-V.x,vx-V.x+cc,vx-V.x,vx-V.x);
}

/**********************************************************************
 * HighTemp2D_pState::lambda_x -- Eigenvalue(s) (x-direction).    *
 * For High-Temperature Air, average State only                   *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_cState::lambda_x(double aAvg) const {
  double vx = v().x;
  return HighTemp2D_pState(vx-aAvg,vx,vx,vx+aAvg,vx,vx);
}

/**********************************************************************
 * HighTemp2D_pState::rp_x -- Primitive right eigenvector         *
 *                                (x-direction).                      *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_cState::rp_x(int index) const {
  //assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
  switch(index) {
  case 1 :
    return HighTemp2D_pState(ONE,-c()/rho,ZERO,c2()-(2.0/3.0)*k(),ZERO,ZERO);
    break;
  case 2 :
    return HighTemp2D_pState(ONE,ZERO,ZERO,-(2.0/3.0)*k(),ZERO,ZERO);
    break;
  case 3 :
    return HighTemp2D_pState(ZERO,ZERO,ONE,ZERO,ZERO,ZERO);
    break;
  case 4 :
    return HighTemp2D_pState(ONE,c()/rho,ZERO,c2()-(2.0/3.0)*k(),ZERO,ZERO);
    break;
  case 5 :
    return HighTemp2D_pState(ZERO,ZERO,ZERO,-(2.0/3.0)*rho,ONE,ZERO);
    break;
  case 6 :
    return HighTemp2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE);
    break;
  };
  return HighTemp2D_pState(ONE,-c()/rho,ZERO,c2()-(2.0/3.0)*k(),ZERO,ZERO);
}

/**********************************************************************
 * HighTemp2D_cState::rc_x -- Conserved right eigenvector         *
 *                                (x-direction).                      *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_cState::rc_x(int index) const {
  //assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
  switch(index) {
  case 1 :
    return HighTemp2D_cState(ONE,v().x-c(),v().y,h()-c()*v().x+(2.0/3.0)*k(),k(),omega());
    break;
  case 2 : 
    if (eos_type == EOS_TGAS)
      return HighTemp2D_cState(ONE,v().x,v().y,h()+2.0*k()/3.0-rho*c()*c()/dpde(),k(),omega()); 
    else
      return HighTemp2D_cState(ONE,v().x,v().y,HALF*v().sqr()+gm1i*(g-5.0/3.0)*k(),k(),omega());
    break;
  case 3 :
    //made change from v.y to rho*v.y from before 
    return HighTemp2D_cState(ZERO,ZERO,rho,dv.y,ZERO,ZERO);
    break;
  case 4 :
    return HighTemp2D_cState(ONE,v().x+c(),v().y,h()+c()*v().x+(2.0/3.0)*k(),k(),omega());
    break;
  case 5 :
    if (eos_type == EOS_TGAS)
      return HighTemp2D_cState(ZERO,ZERO,ZERO,rho*(-2.0*rho/(3.0*dpde())+ONE),rho,ZERO); 
    else 
      return HighTemp2D_cState(ZERO,ZERO,ZERO,rho*gm1i*(g-5.0/3.0),rho,ZERO);
    break;
  case 6 :
    return HighTemp2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,rho);
    break;
  };
  return HighTemp2D_cState(ONE,v().x-c(),v().y,h()-c()*v().x+(2.0/3.0)*k(),k(),omega());
}

/**********************************************************************
 * HighTemp2D_cState::rc_x -- Conserved right eigenvector         *
 *           GLAISTER FLUX ONLY!     (x-direction).                      *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_cState::rc_x(int index, double dpde, double dpdrho, double cAvg) const {
  //assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
 switch(index) {
  case 1 :
    return HighTemp2D_cState(ONE,v().x-cAvg,v().y,h()-cAvg*v().x+(2.0/3.0)*k(),k(),omega());
   break;
  case 2 :
    //return HighTemp2D_cState(ONE,v().x,ZERO,HALF*(v().x*v().x-v().y*v().y)+e()-rho*dpdrho/dpde,ZERO,ZERO);
    return HighTemp2D_cState(ONE,v().x,v().y,h()+2.0*k()/3.0-rho*cAvg*cAvg/dpde,k(),omega());
    break;
  case 3 :
    return HighTemp2D_cState(ZERO,ZERO,rho,rho*v().y,ZERO,ZERO);
    break;
  case 4 :
    return HighTemp2D_cState(ONE,v().x+cAvg,v().y,h()+cAvg*v().x+(2.0/3.0)*k(),k(),omega());
    break;
  case 5 :
    return HighTemp2D_cState(ZERO,ZERO,ZERO,rho*(-2.0*rho/(3.0*dpde)+ONE),rho,ZERO);
    break;
  case 6 :
    return HighTemp2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,rho);
    break;
  };
  return HighTemp2D_cState(ONE,v().x-c(),v().y,h()-c()*v().x+(2.0/3.0)*k(),k(),omega());  
}

/**********************************************************************
 * HighTemp2D_cState::lp_x -- Primitive left eigenvector          *
 *                                (x-direction).                      *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_cState::lp_x(int index) const {
  //assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
  switch(index) {
  case 1 :
    return HighTemp2D_pState(k()/(3.0*c2()),-HALF*rho/c(),ZERO,HALF/c2(),rho/(3.0*c2()),ZERO);
    break;
  case 2 :
    return HighTemp2D_pState(ONE-(2.0/3.0)*k()/c2(),ZERO,ZERO,-ONE/c2(),-(2.0/3.0)*rho/c2(),ZERO);
    break;
  case 3 :
    return HighTemp2D_pState(ZERO,ZERO,ONE,ZERO,ZERO,ZERO);
    break;
  case 4 :
    return HighTemp2D_pState(k()/(3.0*c2()),HALF*rho/c(),ZERO,HALF/c2(),rho/(3.0*c2()),ZERO);
    break;
  case 5 :
    return HighTemp2D_pState(ZERO,ZERO,ZERO,ZERO,ONE,ZERO);
    break;
  case 6 :
    return HighTemp2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE);
    break;
  };
  return HighTemp2D_pState(k()/(3.0*c2()),-HALF*rho/c(),ZERO,HALF/c2(),rho/(3.0*c2()),ZERO);
}

/**********************************************************************
 * HighTemp2D_pState::lp_x -- Primitive left eigenvector          *
 *  ONLY 4 GLAISTER FLUX      (x-direction).                      *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_cState::lp_x(int index, double cAvg) const {
  //assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
  // cout<<"in lp_x conservative, using vectors"<<endl;
  
 switch(index) {
  case 1 :
    return HighTemp2D_pState(k()/(3.0*cAvg*cAvg),-HALF*rho/cAvg,ZERO,HALF/(cAvg*cAvg),rho/(3.0*cAvg*cAvg),ZERO);
    break;
  case 2 :
    return HighTemp2D_pState(ONE-(2.0/3.0)*k()/(cAvg*cAvg),ZERO,ZERO,-ONE/(cAvg*cAvg),-(2.0/3.0)*rho/(cAvg*cAvg),ZERO);
    break;
  case 3 :
    return HighTemp2D_pState(ZERO,ZERO,ONE,ZERO,ZERO,ZERO);
    break;
  case 4 :
    return HighTemp2D_pState(k()/(3.0*cAvg*cAvg),HALF*rho/cAvg,ZERO,HALF/(cAvg*cAvg),rho/(3.0*cAvg*cAvg),ZERO);
    break;
  case 5 :
    return HighTemp2D_pState(ZERO,ZERO,ZERO,ZERO,ONE,ZERO);
    break;
  case 6 :
    return HighTemp2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE);
    break;
  };
  return HighTemp2D_pState(k()/(3.0*cAvg*cAvg),-HALF*rho/cAvg,ZERO,HALF/(cAvg*cAvg),rho/(3.0*cAvg*cAvg),ZERO);
}

/********************************************************
 * Useful 2D HighTemp state constants.                     *
 ********************************************************/
const HighTemp2D_pState HighTemp2D_W_STDATM(DENSITY_STDATM,
				      Vector2D_ZERO, PRESSURE_STDATM);

/**********************************************************************
 * HighTemp2DState -- External subroutines.                       *
 **********************************************************************/

extern HighTemp2D_pState Riemann(const HighTemp2D_pState &Wl,
				     const HighTemp2D_pState &Wr);

extern HighTemp2D_pState RoeAverage(const HighTemp2D_pState &Wl,
					const HighTemp2D_pState &Wr);

extern HighTemp2D_pState Translate(const HighTemp2D_pState &W,
				       const Vector2D &V);

extern HighTemp2D_pState Reflect(const HighTemp2D_pState &W,
				     const Vector2D &norm_dir);

extern HighTemp2D_pState Mirror(const HighTemp2D_pState &W,
				    const Vector2D &norm_dir);

extern HighTemp2D_pState WallViscousHeatFlux(const HighTemp2D_pState &W,
						 const Vector2D &norm_dir);

extern HighTemp2D_pState WallViscousIsothermal(const HighTemp2D_pState &W,
						   const Vector2D &norm_dir,
						   const double &Twall);

extern HighTemp2D_pState MovingWallHeatFlux(const HighTemp2D_pState &W,
						const Vector2D &norm_dir,
						const double &Vwall);

extern HighTemp2D_pState MovingWallIsothermal(const HighTemp2D_pState &W,
						  const Vector2D &norm_dir,
						  const double &Vwall,
						  const double &Twall);

extern HighTemp2D_pState BurningSurface(const HighTemp2D_pState &W,
				    const Vector2D &norm_dir);

extern HighTemp2D_pState RinglebFlow(const HighTemp2D_pState &Wdum,
					 const Vector2D &X);

extern HighTemp2D_pState RinglebFlow(const HighTemp2D_pState &Wdum,
					 const Vector2D &X,
					 double &q, double &k);

extern HighTemp2D_pState RinglebFlowAverageState(const HighTemp2D_pState &Wdum,
						     const Vector2D &Y1,
						     const Vector2D &Y2,
						     const Vector2D &Y3,
						     const Vector2D &Y4);

extern HighTemp2D_pState ViscousChannelFlow(const HighTemp2D_pState &Wdum,
						const Vector2D X,
						const Vector2D Vwall,
						const double dp,
						const double length,
						const double height);

extern HighTemp2D_pState ViscousChannelFlowVelocity(const HighTemp2D_pState &Wdum,
							const Vector2D X,
							const Vector2D Vwall,
							const double dp,
							const double length,
							const double height);

extern HighTemp2D_pState ViscousPipeFlow(const HighTemp2D_pState &Wdum,
					     const Vector2D X,
					     const double dp,
					     const double length,
					     const double radius);

extern HighTemp2D_pState TurbulentPipeFlow(const HighTemp2D_pState &Wo,
					       const Vector2D X,
					       const double dp,
					       const double length,
					       const double radius);

extern HighTemp2D_pState FlatPlate(const HighTemp2D_pState &Winf,
				       const Vector2D &X,
				       const double &plate_length,
				       double &eta,
				       double &f,
				       double &fp,
				       double &fpp);

extern HighTemp2D_pState DrivenCavityFlow(const HighTemp2D_pState &Wo,
					      const double &l,
					      const double &Re);

extern HighTemp2D_pState BackwardFacingStep(const HighTemp2D_pState &Wo,
						const Vector2D &X,
						const double &h,
						const double &ho,
						const double &Re,
						const double &M);

extern HighTemp2D_pState BC_Characteristic(const HighTemp2D_pState &Wi,
					       const HighTemp2D_pState &Wo,
					       const Vector2D &norm_dir);

extern HighTemp2D_pState BC_Characteristic_Pressure(const HighTemp2D_pState &Wi,
							const HighTemp2D_pState &Wo,
							const Vector2D &norm_dir);

extern HighTemp2D_pState BC_Characteristic_Mach_Number(const HighTemp2D_pState &Wi,
							   const HighTemp2D_pState &Wo,
							   const Vector2D &norm_dir);

extern HighTemp2D_pState WaveSpeedPos(const HighTemp2D_pState &lambda_a,
					  const HighTemp2D_pState &lambda_l,
					  const HighTemp2D_pState &lambda_r);

extern HighTemp2D_pState WaveSpeedNeg(const HighTemp2D_pState &lambda_a,
					  const HighTemp2D_pState &lambda_l,
					  const HighTemp2D_pState &lambda_r);

extern HighTemp2D_pState WaveSpeedAbs(const HighTemp2D_pState &lambda_a,
					  const HighTemp2D_pState &lambda_l,
					  const HighTemp2D_pState &lambda_r);

extern HighTemp2D_pState HartenFixPos(const HighTemp2D_pState &lambda_a,
					  const HighTemp2D_pState &lambda_l,
					  const HighTemp2D_pState &lambda_r);

extern HighTemp2D_pState HartenFixNeg(const HighTemp2D_pState &lambda_a,
					  const HighTemp2D_pState &lambda_l,
					  const HighTemp2D_pState &lambda_r);

extern HighTemp2D_pState HartenFixAbs(const HighTemp2D_pState &lambda_a,
					  const HighTemp2D_pState &lambda_l,
					  const HighTemp2D_pState &lambda_r);

extern HighTemp2D_cState FluxGodunov_n(const HighTemp2D_pState &Wl,
					   const HighTemp2D_pState &Wr,
					   const Vector2D &norm_dir);

extern HighTemp2D_cState FluxGodunov_n(const HighTemp2D_cState &Ul,
					   const HighTemp2D_cState &Ur,
					   const Vector2D &norm_dir);

extern HighTemp2D_cState FluxGodunov_MB_n(const HighTemp2D_pState &Wl,
					      const HighTemp2D_pState &Wr,
					      const Vector2D &V,
					      const Vector2D &norm_dir);

extern HighTemp2D_pState StateGodunov_n(const HighTemp2D_pState &Wl,
					    const HighTemp2D_pState &Wr,
					    const Vector2D &norm_dir);

extern HighTemp2D_cState FluxRoe(const HighTemp2D_pState &Wl,
				     const HighTemp2D_pState &Wr);

extern HighTemp2D_cState FluxRoe(const HighTemp2D_cState &Ul,
				     const HighTemp2D_cState &Ur);

extern HighTemp2D_cState FluxRoe_n(const HighTemp2D_pState &Wl,
				       const HighTemp2D_pState &Wr,
				       const Vector2D &norm_dir);

extern HighTemp2D_cState FluxRoe_n(const HighTemp2D_cState &Ul,
				       const HighTemp2D_cState &Ur,
				       const Vector2D &norm_dir);

extern HighTemp2D_cState FluxGlaister(const HighTemp2D_pState &Wl,
				     const HighTemp2D_pState &Wr);

extern HighTemp2D_cState FluxGlaister(const HighTemp2D_cState &Ul,
				     const HighTemp2D_cState &Ur);

extern HighTemp2D_cState FluxGlaister_n(const HighTemp2D_pState &Wl,
				       const HighTemp2D_pState &Wr,
				       const Vector2D &norm_dir);

extern HighTemp2D_cState FluxGlaister_n(const HighTemp2D_cState &Ul,
				       const HighTemp2D_cState &Ur,
				       const Vector2D &norm_dir);

extern HighTemp2D_cState FluxAUSMplusUP(const HighTemp2D_pState &Wl,
				      const HighTemp2D_pState &Wr);

extern HighTemp2D_cState FluxAUSMplusUP(const HighTemp2D_cState &Ul,
				      const HighTemp2D_cState &Ur);

extern HighTemp2D_cState FluxAUSMplusUP_n(const HighTemp2D_pState &Wl,
					const HighTemp2D_pState &Wr,
					const Vector2D &norm_dir);

extern HighTemp2D_cState FluxAUSMplusUP_n(const HighTemp2D_cState &Ul,
					const HighTemp2D_cState &Ur,
					const Vector2D &norm_dir);


extern HighTemp2D_cState FluxRoe_MB(const HighTemp2D_pState &Wl,
					const HighTemp2D_pState &Wr,
					const Vector2D &V);

extern HighTemp2D_cState FluxRoe_MB(const HighTemp2D_cState &Ul,
					const HighTemp2D_cState &Ur,
					const Vector2D &V);

extern HighTemp2D_cState FluxRoe_MB_n(const HighTemp2D_pState &Wl,
					  const HighTemp2D_pState &Wr,
					  const Vector2D &V,
					  const Vector2D &norm_dir);

extern HighTemp2D_cState FluxRoe_MB_n(const HighTemp2D_cState &Ul,
					  const HighTemp2D_cState &Ur,
					  const Vector2D &V,
					  const Vector2D &norm_dir);

extern HighTemp2D_cState FluxRusanov(const HighTemp2D_pState &Wl,
					 const HighTemp2D_pState &Wr);

extern HighTemp2D_cState FluxRusanov(const HighTemp2D_cState &Ul,
					 const HighTemp2D_cState &Ur);

extern HighTemp2D_cState FluxRusanov_n(const HighTemp2D_pState &Wl,
					   const HighTemp2D_pState &Wr,
					   const Vector2D &norm_dir);

extern HighTemp2D_cState FluxRusanov_n(const HighTemp2D_cState &Ul,
					   const HighTemp2D_cState &Ur,
					   const Vector2D &norm_dir);

extern HighTemp2D_cState FluxHLLE(const HighTemp2D_pState &Wl,
				      const HighTemp2D_pState &Wr);

extern HighTemp2D_cState FluxHLLE(const HighTemp2D_cState &Ul,
				      const HighTemp2D_cState &Ur);

extern HighTemp2D_cState FluxHLLE_n(const HighTemp2D_pState &Wl,
					const HighTemp2D_pState &Wr,
					const Vector2D &norm_dir);

extern HighTemp2D_cState FluxHLLE_n(const HighTemp2D_cState &Ul,
					const HighTemp2D_cState &Ur,
					const Vector2D &norm_dir);

extern HighTemp2D_cState FluxGHLLE(const HighTemp2D_pState &Wl,
				      const HighTemp2D_pState &Wr);

extern HighTemp2D_cState FluxGHLLE(const HighTemp2D_cState &Ul,
				      const HighTemp2D_cState &Ur);

extern HighTemp2D_cState FluxGHLLE_n(const HighTemp2D_pState &Wl,
					const HighTemp2D_pState &Wr,
					const Vector2D &norm_dir);

extern HighTemp2D_cState FluxGHLLE_n(const HighTemp2D_cState &Ul,
					const HighTemp2D_cState &Ur,
					const Vector2D &norm_dir);

Vector2D HLLE_wavespeeds(
		const HighTemp2D_pState &Wl,
		const HighTemp2D_pState &Wr,
		const Vector2D &norm_dir);

Vector2D GHLLE_wavespeeds(
		const HighTemp2D_pState &Wl,
		const HighTemp2D_pState &Wr,
		const Vector2D &norm_dir);

extern HighTemp2D_cState FluxHLLL(const HighTemp2D_pState &Wl,
				      const HighTemp2D_pState &Wr);

extern HighTemp2D_cState FluxHLLL(const HighTemp2D_cState &Ul,
				      const HighTemp2D_cState &Ur);

extern HighTemp2D_cState FluxHLLL_n(const HighTemp2D_pState &Wl,
					const HighTemp2D_pState &Wr,
					const Vector2D &norm_dir);

extern HighTemp2D_cState FluxHLLL_n(const HighTemp2D_cState &Ul,
					const HighTemp2D_cState &Ur,
					const Vector2D &norm_dir);

extern HighTemp2D_cState FluxGHLLL(const HighTemp2D_pState &Wl,
				      const HighTemp2D_pState &Wr);

extern HighTemp2D_cState FluxGHLLL(const HighTemp2D_cState &Ul,
				      const HighTemp2D_cState &Ur);

extern HighTemp2D_cState FluxGHLLL_n(const HighTemp2D_pState &Wl,
					const HighTemp2D_pState &Wr,
					const Vector2D &norm_dir);

extern HighTemp2D_cState FluxGHLLL_n(const HighTemp2D_cState &Ul,
					const HighTemp2D_cState &Ur,
					const Vector2D &norm_dir);

extern HighTemp2D_cState FluxHLLC(const HighTemp2D_pState &Wl,
				      const HighTemp2D_pState &Wr);

extern HighTemp2D_cState FluxHLLC(const HighTemp2D_cState &Ul,
				      const HighTemp2D_cState &Ur);

extern HighTemp2D_cState FluxHLLC_n(const HighTemp2D_pState &Wl,
					const HighTemp2D_pState &Wr,
					const Vector2D &norm_dir);

extern HighTemp2D_cState FluxHLLC_n(const HighTemp2D_cState &Ul,
					const HighTemp2D_cState &Ur,
					const Vector2D &norm_dir);

extern HighTemp2D_cState FluxVanLeer(const HighTemp2D_pState &Wl,
					 const HighTemp2D_pState &Wr);

extern HighTemp2D_cState FluxVanLeer(const HighTemp2D_cState &Wl,
					 const HighTemp2D_cState &Wr);

extern HighTemp2D_cState FluxVanLeer_n(const HighTemp2D_pState &Wl,
					   const HighTemp2D_pState &Wr,
					   const Vector2D &norm_dir);

extern HighTemp2D_cState FluxVanLeer_n(const HighTemp2D_cState &Wl,
					   const HighTemp2D_cState &Wr,
					   const Vector2D &norm_dir);

extern HighTemp2D_cState FluxVanLeer_MB(const HighTemp2D_pState &Wl,
					    const HighTemp2D_pState &Wr,
					    const Vector2D &V);

extern HighTemp2D_cState FluxVanLeer_MB_n(const HighTemp2D_pState &Wl,
					      const HighTemp2D_pState &Wr,
					      const Vector2D &V,
					      const Vector2D &norm_dir);

extern HighTemp2D_cState FluxAUSM(const HighTemp2D_pState &Wl,
				      const HighTemp2D_pState &Wr);

extern HighTemp2D_cState FluxAUSM(const HighTemp2D_cState &Wl,
				      const HighTemp2D_cState &Wr);

extern HighTemp2D_cState FluxAUSM_n(const HighTemp2D_pState &Wl,
					const HighTemp2D_pState &Wr,
					const Vector2D &norm_dir);

extern HighTemp2D_cState FluxAUSM_n(const HighTemp2D_cState &Wl,
					const HighTemp2D_cState &Wr,
					const Vector2D &norm_dir);

extern HighTemp2D_cState FluxAUSMplus(const HighTemp2D_pState &Wl,
					  const HighTemp2D_pState &Wr);

extern HighTemp2D_cState FluxAUSMplus(const HighTemp2D_cState &Wl,
					  const HighTemp2D_cState &Wr);

extern HighTemp2D_cState FluxAUSMplus_n(const HighTemp2D_pState &Wl,
					    const HighTemp2D_pState &Wr,
					    const Vector2D &norm_dir);

extern HighTemp2D_cState FluxAUSMplus_n(const HighTemp2D_cState &Wl,
					    const HighTemp2D_cState &Wr,
					    const Vector2D &norm_dir);

extern HighTemp2D_cState ViscousFlux_n(const Vector2D &X,
					   HighTemp2D_pState &W,
					   const HighTemp2D_pState &dWdx,
					   const HighTemp2D_pState &dWdy,
					   const Vector2D &norm_dir,
					   int Axisymmetric);

extern double ShearStress(const HighTemp2D_pState &W,
			  const HighTemp2D_pState &dWdx,
			  const HighTemp2D_pState &dWdy,
			  const Vector2D &nhat);

extern double WallShearStress(const HighTemp2D_pState &W1,
			      const Vector2D &X1,
			      const Vector2D &X2,
			      const Vector2D &X3,
			      const Vector2D &nhat);

extern double WallShearStress2(const Vector2D &X,
			       const Vector2D &X1,
			       const HighTemp2D_pState &W1,
			       const HighTemp2D_pState &dW1dx,
			       const HighTemp2D_pState &dW1dy,
			       const Vector2D &nhat);

#endif // _HIGHTEMP2D_STATE_INCLUDED
