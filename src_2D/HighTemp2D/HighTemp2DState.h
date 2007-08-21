/**********************************************************************
 * HighTemp2DState.h: Header file defining 2D HighTemperature         *
 *                    solution state classes.                         *
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

// For high temperature stuff. HTTOL is variable so 
// that we can vary it without recompiling.
// HTTOL is declared here and defined (and initialized) 
// in HighTemp2DState.cc.
extern double HTTOL;

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
 *     k        -- Turbulent kinetic energy.
 *     omega    -- Specific dissipation rate.
 *     tau      -- Viscous stress tensor (laminar and Reynolds).
 *     q        -- Heat flux vector (laminar and turbulent).
 *     g        -- Specific heat ratio.
 *     gm1      -- Return g-1.
 *     gm1i     -- Return 1/(g-1).
 *     R        -- Gas constant.
 *     v1,v2,v3,v4,v5 -- Viscosity law coefficients.
 *     cp       -- Specific heat at constant pressure.
 *     cv       -- Specific heat at constant volume.
 *     T        -- Temperature.
 *     e        -- Specific internal energy.
 *     E        -- Total energy.
 *     h        -- Specific enthalpy.
 *     H        -- Total enthalpy.
 *     a        -- Sound speed.
 *     a2       -- Sound speed squared.
 *     M        -- Mach number.
 *     Mref     -- Reference Mach number
 *     s        -- Specific entropy.
 *     dv       -- Gas momentum.
 *     To       -- Stagnation temperature.
 *     po       -- Stagnation pressure.
 *     ao       -- Stagnation sound speed.
 *     ho       -- Stagnation enthalpy.
 *     mu       -- Dynamic viscosity.
 *     nu       -- Kinematic viscosity.
 *     kappa    -- Thermal conductivity.
 *     Pr       -- Prandtl number.
 *     gamma    -- Specific heat ratio
 *     dpdrho   -- Partial derivative of pressure wrt density.
 *     dpde     -- Partial derivative of pressure wrt internal energy.
 *     dhdrho   -- Partial derivative of internal enthalpy wrt density
 *     dhdp     -- Partial derivative of internal enthalpy wrt pressure.
 *     dTdrho   -- Partial derivative of temperature wrt density.
 *     dTdp     -- Partial derivative of temperature wrt pressure.
 *     muT      -- Turbulent eddy dynamic viscosity.
 *     nuT      -- Turbulent eddy kinematic viscosity.
 *     kappaT   -- Turbulent eddy thermal conductivity.
 *     PrT      -- Turbulent Prandtl number.
 *     epsilon  -- Return the turbulent eddy dissipation.
 *     c        -- Turbulence modified sound speed.
 *     c2       -- Turbulence modified sound speed squared.
 *     pmodified -- Turbulence modified pressure.
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
 *     RoeAverage -- Return Roe average solution state.
 *     RoeAverage_SoundSpeed - Returns Roe average sound speed.
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
  double                         k; //!< Turbulent kinetic energy.
  double                     omega; //!< Specific dissipation rate.
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
			           const double &yplus_outer);
  static void set_gas(char *gas_type);
  static void set_turbulence(const double &C_constant,
		             const double &von_karman,
		             const double &yplus_sub,
		             const double &yplus_buffer,
		             const double &yplus_outer);
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

  //@{ @name Related thermodynamic functions.
  //! Temperature.
  double T(void) const;
  
  //! Specific internal energy.
  double e(void) const;
  
  //! Total energy.
  double E(void) const;
  
  //! Specific enthalpy.
  double h(void) const;
  
  //! Total enthalpy.
  double H(void) const;

  //! Sound speed.
  double a(void) const;

  //! Sound speed squared.
  double a2(void) const;

  //! Mach number.
  double M(void) const;

  //! Specific entropy.
  double s(void) const;

  //! Gas momentum.
  Vector2D dv(void) const;

  //! Gas momentum.
  double dv(const Vector2D &n) const;

  //! Specific heat ratio.
  double gamma(void) const;

  //! Stagnation temperature.
  double To(void) const;

  //! Stagnation pressure.
  double po(void) const;

  //! Stagnation sound speed.
  double ao(void) const;

  //! Stagnation enthalpy.
  double ho(void) const;

  //! Dynamic viscosity.
  double mu(void) const;

  //! Kinematic viscosity.
  double nu(void) const;

  //! Thermal heat conductivity.
  double kappa(void) const;

  //! Prandtl number.
  double Pr(void) const;

  //! dp/drho for high-temperature equation of state.
  double dpdrho(void) const;

  //! dp/de for high-temperature EOS.
  double dpde(void) const;

  //! dh/drho for high-temperature equation of state.
  double dhdrho(void) const;

  //! dh/dp for high-temperature equation of state.
  double dhdp(void) const;

  //! dT/drho for high-temperature equation of state.
  double dTdrho(void) const;

  //! dT/dp for high-temperature equation of state.
  double dTdp(void) const;

  //! for high-temperature equation of state.
  double ddTdrho(void) const;
  double ddTdp(void) const;
  double ddTdpdrho(void) const;

  //@{ @name Turbulence related functions.
  //! Return the turbulent kinetic energy.
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
  //HighTemp2D_cState F(void) const;
  //HighTemp2D_cState F(const Vector2D &V) const;
  HighTemp2D_cState Fx(void) const;
  HighTemp2D_cState Fx(const Vector2D &V) const;
  void dFdU(DenseMatrix &dFdU) const;
  void dFdW(DenseMatrix &dFdW) const;
  //@}

  //@{ @name Viscous solution fluxes and Jacobians.
  HighTemp2D_cState Gx(const HighTemp2D_pState &dWdx) const;
  HighTemp2D_cState Gy(const HighTemp2D_pState &dWdy) const;
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
  void dSvdU(DenseMatrix &dSvdU, 
             const Vector2D &X, 
             const HighTemp2D_pState &dWdy) const;
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

  //@{ @name Roe-average state.
  void RoeAverage(const HighTemp2D_pState &Wl,
                  const HighTemp2D_pState &Wr);
  void RoeAverage(HighTemp2D_pState &Wa,
                  const HighTemp2D_pState &Wl,
                  const HighTemp2D_pState &Wr) const;
  void RoeAverage_SoundSpeed(double &c,
                             double &dpdrho,
                             double &dpde,
                             const HighTemp2D_pState &Wa,
                             const HighTemp2D_pState &Wl,
                             const HighTemp2D_pState &Wr) const;
  //@}

  //@{ @name Index operator.
  double &operator[](int index) {
    assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
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
  }
  
  const double &operator[](int index) const {
    assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
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

  //@{ @name Rotates the solution in the local rotated frame.
  void Rotate(const HighTemp2D_pState &W, const Vector2D &norm_dir);
  //@}

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
 *     E        -- Total energy.
 *     dk       -- Total turbulent kinetic energy.
 *     domega   -- Total specific dissipation rate.
 *     tau      -- Viscous stress tensor (laminar and Reynolds).
 *     q        -- Heat flux vector (laminar and turbulent).
 *     g        -- Specific heat ratio.
 *     gm1      -- Return g-1.
 *     gm1i     -- Return 1/(g-1).
 *     R        -- Gas gas constant.
 *     cp       -- Specific heat at constant pressure.
 *     cv       -- Specific heat at constant volume.
 *     v1,v2,v3,v4,v5 -- Viscosity law coefficients.
 *     v        -- Flow velocity.
 *     p        -- Pressure.
 *     T        -- Temperature.
 *     e        -- Specific internal energy.
 *     h        -- Specific enthalpy.
 *     H        -- Total enthalpy.
 *     a        -- Sound speed.
 *     a2       -- Sound speed squared.
 *     M        -- Mach number.
 *     s        -- Specific entropy.
 *     To       -- Stagnation temperature.
 *     po       -- Stagnation pressure.
 *     ao       -- Stagnation sound speed.
 *     ho       -- Stagnation enthalpy.
 *     mu       -- Dynamic viscosity.
 *     nu       -- Kinematic viscosity.
 *     kappa    -- Thermal conductivity.
 *     Pr       -- Prandtl number.
 *     gamma    -- Specific heat ratio.
 *     dpdrho   -- Partial derivative of pressure wrt density.
 *     dpde     -- Partial derivative of pressure wrt internal energy.
 *     dhdrho   -- Partial derivative of internal enthalpy wrt density
 *     dhdp     -- Partial derivative of internal enthalpy wrt pressure.
 *     dTdrho   -- Partial derivative of temperature wrt density.
 *     dTdp     -- Partial derivative of temperature wrt pressure.
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
 *     W        -- Return primitive solution state.
 *     dUdW     -- Return the Jacobian of the conserved solution
 *                 variables with respect to the primitive solution
 *                 variables.
 *     dWdU     -- Return the Jacobian of the primitive solution
 *                 variables with respect to the conserved solution
 *                 variables.
 *     F        -- Return x-direction solution flux.
 *     Gx       -- Return x-direction viscous solution flux.
 *     Gy       -- Return y-direction viscous solution flux.
 *     dGxdU    -- Return the Jacobian of the x-direction viscous
 *                 solution flux vector with respect to the conserved
 *                 solution variables.
 *     dGydU    -- Return the Jacobian of the y-direction viscous
 *                 solution flux vector with respect to the conserved
 *                 solution variables.
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
  double                     E; //!< Total energy.
  double                    dk; //!< Total turbulent kinetic energy.
  double                domega; //!< Total turbulent specific dissipation rate.
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
			           const double &yplus_outer);
  static void set_gas(char *gas_type);
  static void set_turbulence(const double &C_constant,
		             const double &von_karman_constant,
		             const double &yplus_sub,
		             const double &yplus_buffer,
		             const double &yplus_outer);
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

  //@{ @name Gas thermodynamic functions.
  //! Flow velocity.
  Vector2D v(void) const;

  //! Flow speed.
  double v(const Vector2D &n) const;

  //! Pressure.
  double p(void) const;

  //! Temperature.
  double T(void) const;

  //! Specific internal energy.
  double e(void) const;

  //! Specific enthalpy.
  double h(void) const;

  //! Total enthalpy.
  double H(void) const;

  //! Sound speed.
  double a(void) const;

  //! Sound speed squared.
  double a2(void) const;

  //! Mach number.
  double M(void) const;

  //! Specific entropy.
  double s(void) const;

  //! Stagnation temperature.
  double To(void) const;

  //! Stagnation pressure.
  double po(void) const;

  //! Stagnation sound speed.
  double ao(void) const;

  //! Stagnation enthalpy.
  double ho(void) const;

  //! Dynamic viscosity.
  double mu(void) const;

  //! Kinematic viscosity.
  double nu(void) const;

  //! Specific heat ratio.
  double gamma(void) const;

  //! Thermal heat conductivity.
  double kappa(void) const;

  //! Prandtl number.
  double Pr(void) const;

  //! dp/drho for high-temperature equation of state.
  double dpdrho(void) const;

  //! dp/de for high-temperature equation of state.
  double dpde(void) const;

  //! dh/drho for high-temperature equation of state.
  double dhdrho(void) const;

  //! dh/dp for high-temperature equation of state.
  double dhdp(void) const;

  //! dT/drho for high-temperature equation of state.
  double dTdrho(void) const;

  //! dT/dp for high-temperature equation of state.
  double dTdp(void) const;

  //! for high-temperature equation of state.
  double ddTdrho(void) const;
  double ddTdp(void) const;
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

  // dFdU is possibly not working right now.
  //  -- Alistair Wood Apr 15 2007 
  //void dFdU(DenseMatrix &dFdU) const;

  //@}

  //@{ @name Viscous solution flux and Jacobians.
  HighTemp2D_cState Gx(const HighTemp2D_pState &dWdx) const;
  HighTemp2D_cState Gy(const HighTemp2D_pState &dWdy) const;
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
    assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
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
  }

  const double &operator[](int index) const {
    assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
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


  //@{ @name Rotates the solution in the local rotated frame.
  void Rotate(const HighTemp2D_cState &U, const Vector2D &norm_dir);
  //@}

};

/**********************************************************************
 * HighTemp2D_pState::Copy -- Copy operator.                          *
 **********************************************************************/
inline void HighTemp2D_pState::Copy(const HighTemp2D_pState &W) {
  rho = W.rho; v.x = W.v.x; v.y = W.v.y; p = W.p; k = W.k; omega = W.omega;
}

/**********************************************************************
 * HighTemp2D_pState::Vacuum -- Vacuum operator.                      *
 **********************************************************************/
inline void HighTemp2D_pState::Vacuum(void) {
  rho = ZERO; v.x = ZERO; v.y = ZERO; p = ZERO; k = ZERO; omega = ZERO; tau.zero(); q.zero();
}

/**********************************************************************
 * HighTemp2D_pState::Standard_Atmosphere -- Standard atmosphere      *
 *                                               operator.            *
 **********************************************************************/
inline void HighTemp2D_pState::Standard_Atmosphere(void) {
  rho = DENSITY_STDATM; v.x = ZERO; v.y = ZERO; p = PRESSURE_STDATM; k = ZERO; omega = ZERO; tau.zero(); q.zero();
}

/**********************************************************************
 * HighTemp2D_pState::Unphysical_Properties -- Check for              *
 *                                             unphysical state       *
 *                                             properties.            *
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
 * HighTemp2D_pState::set_static_variables -- Set all static          *
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
}

inline void HighTemp2D_pState::set_static_variables(char *gas_type,
						    const int &EOSType,
						    const int &FlowType,
						    const double &C_constant,
						    const double &von_karman_constant,
						    const double &yplus_sub,
						    const double &yplus_buffer,
						    const double &yplus_outer) {
  // Set gas constants.
  set_gas(gas_type);
  // Set the equation of state type.
  eos_type = EOSType;
  // Set the flow type.
  flow_type = FlowType;
  // Set turbulence constants.
  set_turbulence(C_constant,von_karman_constant,yplus_sub,yplus_buffer,yplus_outer);
}

/**********************************************************************
 * HighTemp2D_pState::set_gas -- Set gas static variables.            *
 **********************************************************************/
inline void HighTemp2D_pState::set_gas(char *gas_type) {
  if (strcmp(gas_type,"AIR") == 0) {
    g  = GAMMA_AIR;
    R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    v1 = AIR_c1; v2 = AIR_c2; v3 = AIR_c3; v4 = AIR_c4; v5 = AIR_c5;
  } else if (strcmp(gas_type,"HTAIR") == 0) {
    R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
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
 * HighTemp2D_pState::set_turbulence -- Set the turbulence static     *
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
 * HighTemp2D_pState::T -- Temperature.                               *
 **********************************************************************/
inline double HighTemp2D_pState::T(void) const {
  switch(eos_type) {
    case EOS_TGAS:
      return Tgas_temp(p, rho);
    case EOS_IDEAL:
    default:
      return p/(rho*R);
  };
}

/**********************************************************************
 * HighTemp2D_pState::e -- Specific internal energy.                  *
 **********************************************************************/
inline double HighTemp2D_pState::e(void) const {
  switch(eos_type) {
    case EOS_TGAS:
      return Tgas_h(p, rho) - p/rho;
    case EOS_IDEAL:
    default:
      return p/(gm1*rho);
  };
}

/**********************************************************************
 * HighTemp2D_pState::E -- Total energy.                              *
 **********************************************************************/
inline double HighTemp2D_pState::E(void) const {
  switch(eos_type) {
    case EOS_TGAS:
      return rho*e() + HALF*rho*v.sqr() + dk();
    case EOS_IDEAL:
    default:
      return p*gm1i + HALF*rho*v.sqr() + dk();
  };
}

/**********************************************************************
 * HighTemp2D_pState::h -- Specific enthalpy.                         *
 **********************************************************************/
inline double HighTemp2D_pState::h(void) const {
  switch(eos_type) {
    case EOS_TGAS:
      return Tgas_h(p,rho) + HALF*v.sqr() + k; 
    case EOS_IDEAL:
    default: 
      return g*gm1i*p/rho + HALF*v.sqr() + k;
  };
}

/**********************************************************************
 * HighTemp2D_pState::H -- Total enthalpy.                            *
 **********************************************************************/
inline double HighTemp2D_pState::H(void) const {
  return h() * rho;
}

/**********************************************************************
 * HighTemp2D_pState::a -- Sound speed.                               *
 **********************************************************************/
inline double HighTemp2D_pState::a(void) const {
  switch(eos_type) {
    case EOS_TGAS:
      //return sqrt(p*dpde()/sqr(rho)+dpdrho());
      return Tgas_a_from_e_rho(e(), rho);
    case EOS_IDEAL:
    default:
      return sqrt(g*p/rho);
  };
}

/**********************************************************************
 * HighTemp2D_pState::a2 -- Sound speed squared.                      *
 **********************************************************************/
inline double HighTemp2D_pState::a2(void) const {
  switch(eos_type) {
    case EOS_TGAS:
      //return p*dpde()/sqr(rho)+dpdrho();
      return sqr(Tgas_a_from_e_rho(e(), rho));
    case EOS_IDEAL:
    default:
      return g*p/rho;
  };
}

/**********************************************************************
 * HighTemp2D_pState::M -- Mach number.                               *
 **********************************************************************/
inline double HighTemp2D_pState::M(void) const {
  return abs(v)/a();
}

/**********************************************************************
 * HighTemp2D_pState::s -- Specific entropy.                          *
 **********************************************************************/
inline double HighTemp2D_pState::s(void) const {
  switch(eos_type) {
    case EOS_TGAS:
      return Tgas_s(e(),rho);
    case EOS_IDEAL:
    default:
      return R*gm1i*log(p/pow(rho,g));
  };
}

/**********************************************************************
 * HighTemp2D_pState::dv -- Gas momentum.                             *
 **********************************************************************/
inline Vector2D HighTemp2D_pState::dv(void) const {
  return rho*v;
}

/**********************************************************************
 * HighTemp2D_pState::dv -- Gas momentum.                             *
 **********************************************************************/
inline double HighTemp2D_pState::dv(const Vector2D &n) const {
  return rho*(v*n);
}

/**********************************************************************
 * HighTemp2D_pState::To -- Stagnation temperature.                   *
 **********************************************************************/
inline double HighTemp2D_pState::To(void) const {
  return (p/(rho*R))*(ONE+HALF*gm1*v.sqr()/(g*p/rho));
}

/**********************************************************************
 * HighTemp2D_pState::po -- Stagnation pressure.                      *
 **********************************************************************/
inline double HighTemp2D_pState::po(void) const {
  return p*pow(ONE+HALF*gm1*v.sqr()/(g*p/rho),g*gm1i);
}

/**********************************************************************
 * HighTemp2D_pState::ao -- Stagnation sound speed.                   *
 **********************************************************************/
inline double HighTemp2D_pState::ao(void) const {
  return sqrt((g*p/rho)*(ONE+HALF*gm1*v.sqr()/(g*p/rho)));
}

/**********************************************************************
 * HighTemp2D_pState::ho -- Stagnation enthalpy.                      *
 **********************************************************************/
inline double HighTemp2D_pState::ho(void) const {
  return (g*p/(gm1*rho) + HALF*v.sqr())*(ONE+HALF*gm1*v.sqr()/(g*p/rho));
}

/**********************************************************************
 * HighTemp2D_pState::mu -- Dynamic viscosity.                        *
 **********************************************************************/
inline double HighTemp2D_pState::mu(void) const {
  switch(eos_type) {
    case EOS_TGAS:
      return Tgas_mu(Tgas_temp(p,rho), rho);
    case EOS_IDEAL:
    default:
      return mu_gottlieb(v1,v2,v3,v4,v5,T()); 
  };
}

/**********************************************************************
 * HighTemp2D_pState::nu -- Kinematic viscosity.                      *
 **********************************************************************/
inline double HighTemp2D_pState::nu(void) const {
  return mu()/rho;
}

/**********************************************************************
 * HighTemp2D_pState::gamma -- Specific heat ratio.                   *
 **********************************************************************/
inline double HighTemp2D_pState::gamma(void) const {
  switch(eos_type) {
    case EOS_TGAS: {
      double temp, g2;
      temp = Tgas_temp(p,rho);
      g2 = Tgas_gamma(p,temp);
      return g2; }
    case EOS_IDEAL:
    default:
      return g;
  };
}

/**********************************************************************
 * HighTemp2D_pState::kappa -- Thermal heat conductivity.             *
 **********************************************************************/
inline double HighTemp2D_pState::kappa(void) const {
  switch(eos_type) {
    case EOS_TGAS:
      return Tgas_kappa(e(),rho);
    case EOS_IDEAL:
    default:
      return kappa_gottlieb(v1,v2,v3,v4,v5,T(),g,cp);
  };
}

/**********************************************************************
 * HighTemp2D_pState::Pr -- Prandtl number.                           *
 **********************************************************************/
inline double HighTemp2D_pState::Pr(void) const {
  switch(eos_type){
    case EOS_TGAS:
      return Tgas_Pr(Tgas_temp(p,rho), rho);
    case EOS_IDEAL:
    default:
      return cp*mu()/kappa();
  };
}

/**********************************************************************
 * HighTemp2D_pState::dk -- Total turbulent kinetic energy.           *
 **********************************************************************/
inline double HighTemp2D_pState::dk(void) const {
  return rho*k;
}

/**********************************************************************
 * HighTemp2D_pState::domega -- Total turbulent specific dissipation. *
 **********************************************************************/
inline double HighTemp2D_pState::domega(void) const {
  return rho*omega;
}

/**********************************************************************
 * HighTemp2D_pState::epsilon -- Specific turbulent eddy dissipation. *
 **********************************************************************/
inline double HighTemp2D_pState::epsilon(void) const {
  return beta_k_o*k*omega;
}

/**********************************************************************
 * HighTemp2D_pState::depsilon -- Total turbulent eddy dissipation.   *
 **********************************************************************/
inline double HighTemp2D_pState::depsilon(void) const {
  return rho*epsilon();
}

/**********************************************************************
 * HighTemp2D_pState::ell -- Turbulent length scale.                  *
 **********************************************************************/
inline double HighTemp2D_pState::ell(void) const {
  return sqrt(k)/max(omega, NANO);
}

/**********************************************************************
 * HighTemp2D_pState::Mt -- Turbulent Mach number.                    *
 **********************************************************************/
inline double HighTemp2D_pState::Mt(void) const {
  return sqrt(TWO*k/a2());
}

/**********************************************************************
 * HighTemp2D_pState::Mt2 -- Turbulent Mach number squared.           *
 **********************************************************************/
inline double HighTemp2D_pState::Mt2(void) const {
  return TWO*k/a2();
}

/**********************************************************************
 * HighTemp2D_pState::muT -- Turbulent eddy dynamic viscosity.        *
 **********************************************************************/
inline double HighTemp2D_pState::muT(void) const {
  return rho*nuT();
}

/**********************************************************************
 * HighTemp2D_pState::nuT -- Turbulent eddy kinematic viscosity.      *
 **********************************************************************/
inline double HighTemp2D_pState::nuT(void) const {
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) return k/max(omega,TOLER);
  return ZERO;
}

/**********************************************************************
 * HighTemp2D_pState::kappaT -- Turbulent eddy thermal heat           *
 *                              conductivity.                         *
 **********************************************************************/
inline double HighTemp2D_pState::kappaT(void) const {
  if (flow_type != FLOWTYPE_TURBULENT_RANS_K_OMEGA) return ZERO;
  switch(eos_type) {
    case EOS_TGAS:
      return muT()*Tgas_cp(p,Tgas_temp(p,rho))/PrT;  
    case EOS_IDEAL:
    default:
      return muT()*cp/PrT;
  };
}

/**********************************************************************
 * HighTemp2D_pState::dpdrho -- dp/drho                               *
 **********************************************************************/
inline double HighTemp2D_pState::dpdrho(void) const {
  switch(eos_type) {
    case EOS_TGAS: {
      double ex; ex = e();
      return (Tgas_p(ex, rho*(ONE+HTTOL)) - 
              Tgas_p(ex, rho*(ONE-HTTOL)))/(TWO*HTTOL*rho); }
    case EOS_IDEAL:
    default:
      return p/rho;
  };
}

/**********************************************************************
 * HighTemp2D_pState::dpde -- dp/de                                   *
 **********************************************************************/
inline double HighTemp2D_pState::dpde(void) const {
  switch(eos_type) {
    case EOS_TGAS: {
      double ex; ex = e();
      return (Tgas_p(ex*(ONE+HTTOL), rho) - 
              Tgas_p(ex*(ONE-HTTOL), rho))/(TWO*HTTOL*ex); }
    case EOS_IDEAL:
    default:
      return rho*gm1;
  };
}

/**********************************************************************
 * HighTemp2D_pState::dhdrho -- dh/drho                               *
 **********************************************************************/
inline double HighTemp2D_pState::dhdrho(void) const {
  switch(eos_type) {
    case EOS_TGAS:
      return (Tgas_h(p, rho*(ONE+HTTOL)) -  
              Tgas_h(p, rho*(ONE-HTTOL)))/(TWO*HTTOL*rho);
    case EOS_IDEAL:
    default:
      return -g*gm1i*p/sqr(rho);
  };
}

/**********************************************************************
 * HighTemp2D_pState::dhdp -- dh/dp                                   *
 **********************************************************************/
inline double HighTemp2D_pState::dhdp(void) const { 
  switch(eos_type) {
    case EOS_TGAS:
      return (Tgas_h(p*(ONE+HTTOL), rho) - 
              Tgas_h(p*(ONE-HTTOL), rho))/(TWO*HTTOL*p);
    case EOS_IDEAL:
    default:
      return g*gm1i/rho;
  };
}

/**********************************************************************
 * HighTemp2D_pState::dTdrho -- dT/drho                               *
 **********************************************************************/
inline double HighTemp2D_pState::dTdrho(void) const {
  switch(eos_type) {
    case EOS_TGAS:
      return (Tgas_temp(p, rho*(ONE+HTTOL)) - 
              Tgas_temp(p, rho*(ONE-HTTOL)))/(TWO*HTTOL*rho);
    case EOS_IDEAL:
    default:
      return -p/(sqr(rho)*R);
  };
}

/**********************************************************************
 * HighTemp2D_pState::dTdp -- dT/dp                                   *
 **********************************************************************/
inline double HighTemp2D_pState::dTdp(void) const {
  switch(eos_type) {
    case EOS_TGAS:
      return (Tgas_temp(p*(ONE+HTTOL), rho) -
              Tgas_temp(p*(ONE-HTTOL), rho))/(TWO*HTTOL*p);
    case EOS_IDEAL:
    default:
      return ONE/(rho*R);
  };
}

/**********************************************************************
 * HighTemp2D_pState::ddTdrho -- d^2T/drho^2                          *
 **********************************************************************/
inline double HighTemp2D_pState::ddTdrho(void) const {
  switch(eos_type) {
    case EOS_TGAS:
      return (Tgas_temp(p, rho*(ONE+HTTOL)) -
              TWO*Tgas_temp(p,rho) + 
              Tgas_temp(p, rho*(ONE-HTTOL)))/sqr(HTTOL*rho);
    case EOS_IDEAL:
    default:
      return TWO*p/(cube(rho)*R);
  };
}

/**********************************************************************
 * HighTemp2D_pState::ddTdp -- d^2T/dp^2                              *
 **********************************************************************/
inline double HighTemp2D_pState::ddTdp(void) const {
  switch(eos_type) {
    case EOS_TGAS:
      return (Tgas_temp(p*(ONE+HTTOL), rho) -
              TWO*Tgas_temp(p, rho) +
              Tgas_temp(p*(ONE-HTTOL), rho))/sqr(HTTOL*p);
    case EOS_IDEAL:
    default:
      return ZERO;
  };
}

/**********************************************************************
 * HighTemp2D_pState::ddTdrho -- d^2T/dpdrho                          *
 **********************************************************************/
inline double HighTemp2D_pState::ddTdpdrho(void) const {
  switch(eos_type) {
    case EOS_TGAS: {
      double rho1, rho2, p1, p2;
      double dTdp1, dTdp2;
      rho1 = rho*(ONE-HTTOL);
      rho2 = rho*(ONE+HTTOL);
      p1 = p*(ONE-HTTOL);
      p2 = p*(ONE+HTTOL);
      dTdp1 = (Tgas_temp(p2, rho1) - Tgas_temp(p1, rho1))/(TWO*HTTOL*p);
      dTdp2 = (Tgas_temp(p2, rho2) - Tgas_temp(p1, rho2))/(TWO*HTTOL*p);
      return (dTdp2 - dTdp1)/(TWO*HTTOL*rho); }
    case EOS_IDEAL:
    default:
      return -ONE/(sqr(rho)*R);
  };
}

/**********************************************************************
 * HighTemp2D_pState::c -- Turbulence modified sound speed.           *
 **********************************************************************/
inline double HighTemp2D_pState::c(void) const {
  return sqrt(c2());
}

/**********************************************************************
 * HighTemp2D_pState::c2 -- Turbulence modified sound speed squared.  *
 **********************************************************************/
inline double HighTemp2D_pState::c2(void) const {
  switch(eos_type){
    case EOS_TGAS:
      return a2() + (TWO/THREE)*k + (TWO*k*dpde())/(THREE*rho);
    case EOS_IDEAL:
    default:
      return a2() + (TWO/THREE)*g*k;
  };
}

/**********************************************************************
 * HighTemp2D_pState::pmodified -- Turbulence modified pressure.      *
 **********************************************************************/
inline double HighTemp2D_pState::pmodified(void) const {
  return p + (2.0/3.0)*dk();
}

/**********************************************************************
 * HighTemp2D_pState::beta_k -- k-omega auxilary relation.            *
 **********************************************************************/
inline double HighTemp2D_pState::beta_k(const HighTemp2D_pState &dWdx,
                                        const HighTemp2D_pState &dWdy) const {
   return beta_k_o*f_beta_k(dWdx,dWdy)*
          (ONE + xi*(Mt2()-sqr(Mto))*heaviside(Mt()-Mto));
}

/**********************************************************************
 * HighTemp2D_pState::beta_omega -- k-omega auxilary relation.        *
 **********************************************************************/
inline double HighTemp2D_pState::beta_omega(const HighTemp2D_pState &dWdx,
                                            const HighTemp2D_pState &dWdy) const {
  return beta_omega_o*f_beta_omega(dWdx,dWdy) -
         beta_k_o*f_beta_k(dWdx,dWdy)*xi*(Mt2()-sqr(Mto))*heaviside(Mt()-Mto);
}

/**********************************************************************
 * HighTemp2D_pState::f_beta_k -- k-omega auxilary relation.          *
 **********************************************************************/
inline double HighTemp2D_pState::f_beta_k(const HighTemp2D_pState &dWdx,
                                          const HighTemp2D_pState &dWdy) const {
  double chi = chi_k(dWdx,dWdy);
  if (chi <= ZERO) return ONE;
  return (ONE + 680.0*chi*chi)/(ONE + 400.0*chi*chi);
}

/**********************************************************************
 * HighTemp2D_pState::f_beta_omega -- k-omega auxilary relation.      *
 **********************************************************************/
inline double HighTemp2D_pState::f_beta_omega(const HighTemp2D_pState &dWdx,
                                              const HighTemp2D_pState &dWdy) const {
  double chi = chi_omega(dWdx,dWdy);
  return (ONE + 70.0*chi)/(ONE + 80.0*chi);
}

/**********************************************************************
 * HighTemp2D_pState::chi_k -- k-omega auxilary relation.             *
 **********************************************************************/
inline double HighTemp2D_pState::chi_k(const HighTemp2D_pState &dWdx,
                                       const HighTemp2D_pState &dWdy) const {
  return (dWdx.k*dWdx.omega + dWdy.k*dWdy.omega)/max(TOLER,cube(omega));
}

/**********************************************************************
 * HighTemp2D_pState::chi_omega -- k-omega auxilary relation.         *
 **********************************************************************/
inline double HighTemp2D_pState::chi_omega(const HighTemp2D_pState &dWdx,
                                           const HighTemp2D_pState &dWdy) const {
  return 0.25*fabs(sqr(dWdx.v.y - dWdy.v.x)*(dWdx.v.x + dWdy.v.y)/
         max(TOLER,cube(beta_omega_o*omega)));
}

/**********************************************************************
 * HighTemp2D_pState -- Binary arithmetic operators.                  *
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
 * HighTemp2D_pState -- Assignment operator.                          *
 **********************************************************************/
inline HighTemp2D_pState& HighTemp2D_pState::operator =(const HighTemp2D_pState &W) {
  rho = W.rho; v.x = W.v.x; v.y = W.v.y; p = W.p; k = W.k; omega = W.omega;
  return *this;
}

/**********************************************************************
 * HighTemp2D_pState -- Unary arithmetic operators.                   *
 **********************************************************************/
//inline HighTemp2D_pState operator +(const HighTemp2D_pState &W) {
//return W;
//}

inline HighTemp2D_pState operator -(const HighTemp2D_pState &W) {
  return HighTemp2D_pState(-W.rho,-W.v.x,-W.v.y,-W.p,-W.k,-W.omega);
}

/**********************************************************************
 * HighTemp2D_pState -- Shortcut arithmetic operators.                *
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
 * HighTemp2D_pState -- Relational operators.                         *
 **********************************************************************/
inline int operator ==(const HighTemp2D_pState &W1, const HighTemp2D_pState &W2) {
  return (W1.rho == W2.rho && W1.v == W2.v && W1.p == W2.p && 
          W1.k == W2.k && W1.omega == W2.omega);
}

inline int operator !=(const HighTemp2D_pState &W1, const HighTemp2D_pState &W2) {
  return (W1.rho != W2.rho || W1.v != W2.v || W1.p != W2.p || 
          W1.k != W2.k || W1.omega != W2.omega);
}

/**********************************************************************
 * HighTemp2D_pState -- Input-output operators.                       *
 **********************************************************************/
inline ostream &operator << (ostream &out_file, const HighTemp2D_pState &W) {
  out_file.setf(ios::scientific);
  out_file << " " << W.rho << " " << W.v.x << " " << W.v.y << " " << W.p 
           << " " << W.k << " " << W.omega;
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
 * HighTemp2D_cState::Copy -- Copy operator.                          *
 **********************************************************************/
inline void HighTemp2D_cState::Copy(const HighTemp2D_cState &U) {
  rho = U.rho; dv.x = U.dv.x; dv.y = U.dv.y; E = U.E; 
  dk = U.dk; domega = U.domega;
}

/**********************************************************************
 * HighTemp2D_cState::Vacuum -- Vacuum state.                         *
 **********************************************************************/
inline void HighTemp2D_cState::Vacuum(void) {
  rho = ZERO; dv.x = ZERO; dv.y = ZERO; E = ZERO; dk = ZERO; 
  domega = ZERO; tau.zero(); q.zero();
}

/**********************************************************************
 * HighTemp2D_cState::Standard_Atmosphere -- Standard atmosphere      *
 *                                               state.               *
 **********************************************************************/
inline void HighTemp2D_cState::Standard_Atmosphere(void) {
  rho = DENSITY_STDATM; dv.x = ZERO; dv.y = ZERO;
  E = PRESSURE_STDATM/(GAMMA_AIR-ONE); tau.zero(); q.zero();
  dk = ZERO; domega = ZERO;
}

/**********************************************************************
 * HighTemp2D_cState::Unphysical_Properties -- Check for              *
 *                                             unphysical state       *
 *                                             properties.            *
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
 * HighTemp2D_cState::Copy_Multigrid_State_Variables --               *
 *                           Copy variables solved by multigrid only. *
 **********************************************************************/
inline void HighTemp2D_cState::Copy_Multigrid_State_Variables(const HighTemp2D_cState &Ufine) {
  rho = Ufine.rho; dv.x = Ufine.dv.x; dv.y = Ufine.dv.y; E = Ufine.E;
}

/**********************************************************************
 * HighTemp2D_cState::Zero_Non_Multigrid_State_Variables --           *
 *                            Zero variables not-solved by multigrid. *
 **********************************************************************/
inline void HighTemp2D_cState::Zero_Non_Multigrid_State_Variables(void) {
  dk = ZERO; domega = ZERO;
}

/**********************************************************************
 * HighTemp2D_cState::set_static_variables -- Set all static          *
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
}

inline void HighTemp2D_cState::set_static_variables(char *gas_type,
						    const int &EOSType,
                                                    const int &FlowType,
                                                    const double &C_constant,
                                                    const double &von_karman_constant,
                                                    const double &yplus_sub,
                                                    const double &yplus_buffer,
                                                    const double &yplus_outer) {
  // Set gas constants.
  set_gas(gas_type);
  // Set Equation of State Type
  eos_type = EOSType;
  // Set the flow type.
  flow_type = FlowType;
  // Set turbulence constants.
  set_turbulence(C_constant,von_karman_constant,yplus_sub,yplus_buffer,yplus_outer);
}

/**********************************************************************
 * HighTemp2D_cState::set_gas -- Set gas static variables.            *
 **********************************************************************/
inline void HighTemp2D_cState::set_gas(char *gas_type) {
  if (strcmp(gas_type,"AIR") == 0) {
    g  = GAMMA_AIR;
    R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
    v1 = AIR_c1; v2 = AIR_c2; v3 = AIR_c3; v4 = AIR_c4; v5 = AIR_c5;
  }else if (strcmp(gas_type,"HTAIR") == 0) {
    R  = R_UNIVERSAL/(MOLE_WT_AIR*MILLI);
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
 }
}

/**********************************************************************
 * HighTemp2D_cState::set_turbulence -- Set the turbulence static     *
 *                                      variables.                    *
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
 * HighTemp2D_cState::v -- Gas flow velocity.                         *
 **********************************************************************/
inline Vector2D HighTemp2D_cState::v(void) const {
  return dv/rho;
}

/*********************************************************************
 * HighTemp2D_cState::v -- Gas flow velocity.                        *
**********************************************************************/
inline double HighTemp2D_cState::v(const Vector2D &n) const {
  return (dv*n)/rho;
}

/*********************************************************************
 * HighTemp2D_cState::p -- Pressure.                                 *
 *********************************************************************/
inline double HighTemp2D_cState::p(void) const {
  switch(eos_type){
    case EOS_TGAS:
      return Tgas_p(e(), rho);
    case EOS_IDEAL:
    default:
      return gm1*(E - HALF*dv.sqr()/rho - dk);
  }
}

/**********************************************************************
 * HighTemp2D_cState::T -- Temperature.                               *
 **********************************************************************/
inline double HighTemp2D_cState::T(void) const {
  switch(eos_type){
    case EOS_TGAS:
      return Tgas_temp(p(),rho);
    case EOS_IDEAL:
    default: 
      return p()/(rho*R);
  };
}

/**********************************************************************
 * HighTemp2D_cState::e -- Specific internal energy.                  *
 **********************************************************************/
inline double HighTemp2D_cState::e(void) const {
  return (E - HALF*dv.sqr()/rho - dk)/rho;
}

/**********************************************************************
 * HighTemp2D_cState::h -- Specific total enthalpy.                   *
 **********************************************************************/
inline double HighTemp2D_cState::h(void) const {
  return H() / rho;
}

/**********************************************************************
 * HighTemp2D_cState::H -- Total enthalpy.                            *
 **********************************************************************/
inline double HighTemp2D_cState::H(void) const {
  return E + p();
}

/**********************************************************************
 * HighTemp2D_cState::a -- Sound speed.                               *
 **********************************************************************/
inline double HighTemp2D_cState::a(void) const {
  switch(eos_type){
    case EOS_TGAS:  
      //return sqrt(p()*dpde()/sqr(rho)+dpdrho());
      return Tgas_a_from_e_rho(e(), rho); 
    case EOS_IDEAL: 
    default:
      return sqrt(g*p()/rho);
  };
}

/**********************************************************************
 * HighTemp2D_cState::a2 -- Sound speed squared.                      *
 **********************************************************************/
inline double HighTemp2D_cState::a2(void) const {
  switch(eos_type){
    case EOS_TGAS:
      //return p()*dpde()/sqr(rho)+dpdrho();
      return sqr(Tgas_a_from_e_rho(e(), rho));
    case EOS_IDEAL:
    default:
      return g*p()/rho;
  };
}

/**********************************************************************
 * HighTemp2D_cState::M -- Mach number.                               *
 **********************************************************************/
inline double HighTemp2D_cState::M(void) const {
  return abs(v())/a();
}

/**********************************************************************
 * HighTemp2D_cState::s -- Specific entropy.                          *
 **********************************************************************/
inline double HighTemp2D_cState::s(void) const {
  switch(eos_type){
    case EOS_TGAS:
      return Tgas_s(e(),rho);
    case EOS_IDEAL:
    default:
      return R*gm1i*log(p()/pow(rho,g));
  };
}

/**********************************************************************
 * HighTemp2D_cState::To -- Stagnation temperature.                   *
 **********************************************************************/
inline double HighTemp2D_cState::To(void) const {
  return (gm1*(E - HALF*dv.sqr()/rho)/(rho*R))*
	 (ONE+HALF*gm1*dv.sqr()/(rho*rho*g*gm1*(E/rho - HALF*dv.sqr()/sqr(rho))));
}

/**********************************************************************
 * HighTemp2D_cState::po -- Stagnation pressure.                      *
 **********************************************************************/
inline double HighTemp2D_cState::po(void) const {
  return (gm1*(E - HALF*dv.sqr()/rho))*
	  pow(ONE+HALF*gm1*dv.sqr()/(rho*rho*g*gm1*(E/rho - HALF*dv.sqr()/sqr(rho))),g*gm1i);
}

/**********************************************************************
 * HighTemp2D_cState::ao -- Stagnation sound speed.                   *
 **********************************************************************/
inline double HighTemp2D_cState::ao(void) const {
  return sqrt((g*gm1*(E/rho - HALF*dv.sqr()/sqr(rho)))*
	      (ONE+HALF*gm1*dv.sqr()/(rho*rho*g*gm1*(E/rho - HALF*dv.sqr()/sqr(rho)))));
}

/**********************************************************************
 * HighTemp2D_cState::ho -- Stagnation enthalpy.                      *
 **********************************************************************/
inline double HighTemp2D_cState::ho(void) const {
  return (g*E/rho - gm1*HALF*dv.sqr()/sqr(rho))*
         (ONE+HALF*gm1*dv.sqr()/(rho*rho*g*gm1*(E/rho - HALF*dv.sqr()/sqr(rho))));
}

/**********************************************************************
 * HighTemp2D_cState::mu -- Gas dynamic viscosity.                    *
 **********************************************************************/
inline double HighTemp2D_cState::mu(void) const {
  switch(eos_type) {
    case EOS_TGAS:
      return Tgas_mu(Tgas_temp(p(),rho), rho);
    case EOS_IDEAL:
    default:
      return mu_gottlieb(v1,v2,v3,v4,v5,T());
  };
}

/**********************************************************************
 * HighTemp2D_cState::nu -- Kinematic viscosity.                      *
 **********************************************************************/
inline double HighTemp2D_cState::nu(void) const {
  return mu()/rho;
}

/**********************************************************************
 * HighTemp2D_cState::gamma -- Specific heat ratio.                   *
 **********************************************************************/
inline double HighTemp2D_cState::gamma(void) const{
  switch(eos_type) {
    case EOS_TGAS: {
      double px, temp, g2; px = p();
      temp = Tgas_temp(px, rho);
      g2 = Tgas_gamma(px, temp);
      return g2; }
    case EOS_IDEAL:
    default:
      return g;
  };
}

/**********************************************************************
 * HighTemp2D_cState::kappa -- Thermal heat conductivity.             *
 **********************************************************************/
inline double HighTemp2D_cState::kappa(void) const {
  switch(eos_type){
    case EOS_TGAS:
      return Tgas_kappa(e(),rho);
    case EOS_IDEAL:
    default:
      return kappa_gottlieb(v1,v2,v3,v4,v5,T(),g,cp);
  };
}

/**********************************************************************
 * HighTemp2D_cState::Pr -- Prandtl number.                           *
 **********************************************************************/
inline double HighTemp2D_cState::Pr(void) const {
  switch(eos_type){
    case EOS_TGAS:
      return Tgas_Pr(Tgas_temp(p(),rho), rho);
    case EOS_IDEAL:
    default:
      return cp*mu()/kappa();
  };
}

/**********************************************************************
 * HighTemp2D_cState::dpdrho -- dp/drho                               *
 **********************************************************************/
inline double HighTemp2D_cState::dpdrho(void) const {
  switch(eos_type) {
    case EOS_TGAS: {
      double ee; ee = e();
      return (Tgas_p(ee, rho*(ONE+HTTOL)) - 
              Tgas_p(ee, rho*(ONE-HTTOL)))/(TWO*HTTOL*rho); }
    case EOS_IDEAL:
    default:
      return p()/rho;
  };
}

/**********************************************************************
 * HighTemp2D_cState::dpde -- dp/de                                   *
 **********************************************************************/
inline double HighTemp2D_cState::dpde(void) const {
  switch(eos_type) {
    case EOS_TGAS: {
      double ee; ee = e();
      return (Tgas_p(ee*(ONE+HTTOL), rho) - 
              Tgas_p(ee*(ONE-HTTOL), rho))/(TWO*HTTOL*ee); }
    case EOS_IDEAL:
    default:
      return rho*gm1;
  };
}

/**********************************************************************
 * HighTemp2D_cState::dhdrho -- dh/drho                               *
 **********************************************************************/
inline double HighTemp2D_cState::dhdrho(void) const {
  switch(eos_type) {
    case EOS_TGAS: {
      double px; px = p();
      return (Tgas_h(px, rho*(ONE+HTTOL)) -  
              Tgas_h(px, rho*(ONE-HTTOL)))/(TWO*HTTOL*rho); }
    case EOS_IDEAL:
    default:
      return -g*gm1i*p()/sqr(rho);
  };

}

/**********************************************************************
 * HighTemp2D_cState::dhdp -- dh/dp                                   *
 **********************************************************************/
inline double HighTemp2D_cState::dhdp(void) const { 
  switch(eos_type) {
    case EOS_TGAS: {
      double px; px = p();
      return (Tgas_h(px*(ONE+HTTOL), rho) - 
              Tgas_h(px*(ONE-HTTOL), rho))/(TWO*HTTOL*px); }
    case EOS_IDEAL:
    default:
      return g*gm1i/rho;
  };
}

/**********************************************************************
 * HighTemp2D_pState::dTdrho -- dT/drho                               *
 **********************************************************************/
inline double HighTemp2D_cState::dTdrho(void) const {
  switch(eos_type) {
    case EOS_TGAS: {
      double px; px = p();
      return (Tgas_temp(px, rho*(ONE+HTTOL)) - 
              Tgas_temp(px, rho*(ONE-HTTOL)))/(TWO*HTTOL*rho); }
    case EOS_IDEAL:
    default:
      return -p()/(sqr(rho)*R);
  };
}

/**********************************************************************
 * HighTemp2D_pState::dTdp -- dT/dp                                   *
 **********************************************************************/
inline double HighTemp2D_cState::dTdp(void) const {
  switch(eos_type) {
    case EOS_TGAS: {
      double px; px = p();
      return (Tgas_temp(px*(ONE+HTTOL), rho) - 
              Tgas_temp(px*(ONE-HTTOL), rho))/(TWO*HTTOL*px); }
    case EOS_IDEAL:
    default:
      return ONE/(rho*R);
  };
}

/**********************************************************************
 * HighTemp2D_pState::ddTdrho -- d^2T/drho^2                          *
 **********************************************************************/
inline double HighTemp2D_cState::ddTdrho(void) const {
  switch(eos_type) {
    case EOS_TGAS: {
      double px; px = p();
      return (Tgas_temp(px, rho*(ONE+HTTOL)) -
              TWO*Tgas_temp(px, rho) + 
              Tgas_temp(px, rho*(ONE-HTTOL)))/sqr(HTTOL*rho); }
    case EOS_IDEAL:
    default:
      return TWO*p()/(cube(rho)*R);
  };
}

/**********************************************************************
 * HighTemp2D_pState::ddTdp -- d^2T/dp^2                              *
 **********************************************************************/
inline double HighTemp2D_cState::ddTdp(void) const {
  switch(eos_type) {
    case EOS_TGAS: {
      double px; px = p();
      return (Tgas_temp(px*(ONE+HTTOL), rho) -
              TWO*Tgas_temp(px, rho) + 
              Tgas_temp(px*(ONE-HTTOL), rho))/sqr(HTTOL*px); }
    case EOS_IDEAL:
    default:
      return ZERO;
  };
}

/**********************************************************************
 * HighTemp2D_pState::ddTdrho -- d^2T/dpdrho                          *
 **********************************************************************/
inline double HighTemp2D_cState::ddTdpdrho(void) const {
  switch(eos_type) {
    case EOS_TGAS: {
      double rho1, rho2, p1, p2;
      double dTdp1, dTdp2; 
      double px;
      px = p();
      rho1 = rho*(ONE-HTTOL);
      rho2 = rho*(ONE+HTTOL);
      p1 = px*(ONE-HTTOL);
      p2 = px*(ONE+HTTOL);
      dTdp1 = (Tgas_temp(p2, rho1) - Tgas_temp(p1, rho1))/(TWO*HTTOL*px);
      dTdp2 = (Tgas_temp(p2, rho2) - Tgas_temp(p1, rho2))/(TWO*HTTOL*px);
      return (dTdp2 - dTdp1)/(TWO*HTTOL*rho); }
    case EOS_IDEAL:
    default:
      return -ONE/(sqr(rho)*R);
  };
}

/**********************************************************************
 * HighTemp2D_cState::dk -- Specific turbulent kinetic energy.        *
 **********************************************************************/
inline double HighTemp2D_cState::k(void) const {
  return dk/rho;
}

/**********************************************************************
 * HighTemp2D_cState::depsilon -- Total turbulent eddy dissipation.   *
 **********************************************************************/
inline double HighTemp2D_cState::depsilon(void) const {
  return beta_k_o*dk*domega/rho;
}

/**********************************************************************
 * HighTemp2D_cState::depsilon -- Specific turbulent eddy dissipation.*
 **********************************************************************/
inline double HighTemp2D_cState::epsilon(void) const {
  return depsilon()/rho;
}

/**********************************************************************
 * HighTemp2D_cState::omega -- Specific turbulent dissipation rate.   *
 **********************************************************************/
inline double HighTemp2D_cState::omega(void) const {
  return domega/rho;
}

/**********************************************************************
 * HighTemp2D_cState::ell -- Turbulent length scale.                  *
 **********************************************************************/
inline double HighTemp2D_cState::ell(void) const {
  return sqrt(k())/max(omega(), NANO);
}

/**********************************************************************
 * HighTemp2D_cState::Mt -- Turbulent Mach number.                    *
 **********************************************************************/
inline double HighTemp2D_cState::Mt(void) const {
  return sqrt(TWO*k()/a2());
}

/**********************************************************************
 * HighTemp2D_cState::Mt2 -- Turbulent Mach number squared.           *
 **********************************************************************/
inline double HighTemp2D_cState::Mt2(void) const {
  return TWO*k()/a2();
}

/**********************************************************************
 * HighTemp2D_cState::muT -- Turbulent eddy dynamic viscosity.        *
 **********************************************************************/
inline double HighTemp2D_cState::muT(void) const {
  return rho*nuT();
}

/**********************************************************************
 * HighTemp2D_cState::nuT -- Turbulent eddy kinematic viscosity.      *
 **********************************************************************/
inline double HighTemp2D_cState::nuT(void) const {
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) return k()/max(omega(),TOLER);
  return ZERO;
}

/**********************************************************************
 * HighTemp2D_cState::kappaT -- Turbulent eddy thermal heat           *
 *                              conductivity.                         *
 **********************************************************************/
inline double HighTemp2D_cState::kappaT(void) const {
  if (flow_type != FLOWTYPE_TURBULENT_RANS_K_OMEGA) return ZERO;
  switch(eos_type) {
    case EOS_TGAS: {
      double px; px = p();
      return muT()*Tgas_cp(px,Tgas_temp(px,rho))/PrT; }
    case EOS_IDEAL:
    default:
      return muT()*cp/PrT;
  };
}

/**********************************************************************
 * HighTemp2D_cState::c -- Turbulence modified sound speed.           *
 **********************************************************************/
inline double HighTemp2D_cState::c(void) const {
  return sqrt(c2());
}

/**********************************************************************
 * HighTemp2D_cState::c2 -- Turbulence modified sound speed squared.  *
 **********************************************************************/
inline double HighTemp2D_cState::c2(void) const {
  switch(eos_type){
    case EOS_TGAS:
      return a2() + (TWO/THREE)*k() + (TWO*k()*dpde())/(THREE*rho);
    case EOS_IDEAL:
    default:
      return a2() + (TWO/THREE)*g*k();
  };
}

/**********************************************************************
 * HighTemp2D_cState::pmodified -- Turbulence modified pressure.      *
 **********************************************************************/
inline double HighTemp2D_cState::pmodified(void) const {
  return p() + (2.0/3.0)*dk;
}

/**********************************************************************
 * HighTemp2D_cState::beta_k -- k-omega auxilary relation.            *
 **********************************************************************/
inline double HighTemp2D_cState::beta_k(const HighTemp2D_pState &dWdx,
		                        const HighTemp2D_pState &dWdy) const {
  return beta_k_o*f_beta_k(dWdx,dWdy)*(ONE + xi*(Mt2()-sqr(Mto))*heaviside(Mt()-Mto));
}

/**********************************************************************
 * HighTemp2D_cState::beta_omega -- k-omega auxilary relation.        *
 **********************************************************************/
inline double HighTemp2D_cState::beta_omega(const HighTemp2D_pState &dWdx,
				            const HighTemp2D_pState &dWdy) const {
  return beta_omega_o*f_beta_omega(dWdx,dWdy) -
         beta_k_o*f_beta_k(dWdx,dWdy)*xi*(Mt2()-sqr(Mto))*heaviside(Mt()-Mto);
}

/**********************************************************************
 * HighTemp2D_cState::f_beta_k -- k-omega auxilary relation.          *
 **********************************************************************/
inline double HighTemp2D_cState::f_beta_k(const HighTemp2D_pState &dWdx,
					  const HighTemp2D_pState &dWdy) const {
  double chi = chi_k(dWdx,dWdy);
  if (chi <= ZERO) return ONE;
  return (ONE + 680.0*chi*chi)/(ONE + 400.0*chi*chi);
}

/**********************************************************************
 * HighTemp2D_cState::f_beta_omega -- k-omega auxilary relation.      *
 **********************************************************************/
inline double HighTemp2D_cState::f_beta_omega(const HighTemp2D_pState &dWdx,
                                              const HighTemp2D_pState &dWdy) const {
  double chi = chi_omega(dWdx,dWdy);
  return (ONE + 70.0*chi)/(ONE + 80.0*chi);
}

/**********************************************************************
 * HighTemp2D_cState::chi_k -- k-omega auxilary relation.             *
 **********************************************************************/
inline double HighTemp2D_cState::chi_k(const HighTemp2D_pState &dWdx,
				       const HighTemp2D_pState &dWdy) const {
  return (dWdx.k*dWdx.omega + dWdy.k*dWdy.omega)/max(TOLER,cube(omega()));
}

/**********************************************************************
 * HighTemp2D_cState::chi_omega -- k-omega auxilary relation.         *
 **********************************************************************/
inline double HighTemp2D_cState::chi_omega(const HighTemp2D_pState &dWdx,
					       const HighTemp2D_pState &dWdy) const {
  return 0.25*fabs(sqr(dWdx.v.y - dWdy.v.x)*(dWdx.v.x + dWdy.v.y)/
         max(TOLER,cube(beta_omega_o*omega())));
}

/**********************************************************************
 * HighTemp2D_cState -- Binary arithmetic operators.                  *
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
 * HighTemp2D_cState -- Assignment operator.                          *
 **********************************************************************/
inline HighTemp2D_cState& HighTemp2D_cState::operator =(const HighTemp2D_cState &U) {
  rho = U.rho; dv.x = U.dv.x; dv.y = U.dv.y; E = U.E; dk = U.dk; domega = U.domega;
  return *this;
}

/**********************************************************************
 * HighTemp2D_cState -- Unary arithmetic operators.                   *
 **********************************************************************/
//inline HighTemp2D_cState operator +(const HighTemp2D_cState &U) {
//return U;
//}

inline HighTemp2D_cState operator -(const HighTemp2D_cState &U) {
  return HighTemp2D_cState(-U.rho,-U.dv.x,-U.dv.y,-U.E,-U.dk,-U.domega);
}

/**********************************************************************
 * HighTemp2D_cState -- Shortcut arithmetic operators.                *
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
 * HighTemp2D_cState -- Relational operators.                         *
 **********************************************************************/
inline int operator ==(const HighTemp2D_cState &U1, const HighTemp2D_cState &U2) {
  return (U1.rho == U2.rho && U1.dv == U2.dv && U1.E == U2.E && 
          U1.dk == U2.dk && U1.domega == U2.domega);
}

inline int operator !=(const HighTemp2D_cState &U1, const HighTemp2D_cState &U2) {
  return (U1.rho != U2.rho || U1.dv != U2.dv || U1.E != U2.E || 
          U1.dk != U2.dk || U1.domega != U2.domega);
}

/**********************************************************************
 * HighTemp2D_cState -- Input-output operators.                       *
 **********************************************************************/
inline ostream &operator << (ostream &out_file, const HighTemp2D_cState &U) {
  out_file.setf(ios::scientific);
  out_file << " " << U.rho << " " << U.dv.x << " " << U.dv.y << " " << U.E << " " 
           << U.dk << " " << U.domega;
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
 * HighTemp2D_cState::Rotate -- Rotates solution state.               *
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
 * HighTemp2D_pState::HighTemp2D_pState -- Constructor.               *
 **********************************************************************/
inline HighTemp2D_pState::HighTemp2D_pState(const HighTemp2D_cState &U) {
  rho = U.rho; v.x = U.v().x; v.y = U.v().y; p = U.p(); k = U.k(); omega = U.omega();
}

/**********************************************************************
 * HighTemp2D_pState::U -- Conserved solution state.                  *
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
 * HighTemp2D_pState::dUdW -- Jacobian of the conserved solution      *
 *                                variables with respect to the       *
 *                                primitive solution variables .      *
 **********************************************************************/
inline void HighTemp2D_pState::dUdW(DenseMatrix &dUdW) const {
  switch (eos_type) {
    case EOS_TGAS: {
      double dpdex; dpdex = dpde();
      double dpdrhox; dpdrhox = dpdrho();
      double ex; ex = e();
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
      break; }
    case EOS_IDEAL:
    default: {
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
      break; }
  }
}

/**********************************************************************
 * HighTemp2D_pState::dWdU -- Jacobian of the primitive solution      *
 *                                variables with respect to the       *
 *                                conserved solution variables.       *
 **********************************************************************/
inline void HighTemp2D_pState::dWdU(DenseMatrix &dWdU) const {
  switch (eos_type) {
    case EOS_TGAS: {
      double dpdex; dpdex = dpde();
      double dpdrhox; dpdrhox = dpdrho();
      double ex; ex = e();
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
      break; }
    case EOS_IDEAL:
    default: {
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
      break; }
  }
}

/**********************************************************************
 * HighTemp2D_pState::Fx -- Solution inviscid flux (x-direction).     *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_pState::Fx(void) const {
  return HighTemp2D_cState(rho*v.x,rho*sqr(v.x)+p+(2.0/3.0)*dk(),
                           rho*v.x*v.y,v.x*H()+v.x*(2.0/3.0)*dk(),
                           rho*v.x*k,rho*v.x*omega);
}

inline HighTemp2D_cState HighTemp2D_pState::Fx(const Vector2D &V) const {
  return HighTemp2D_cState(rho*(v.x-V.x),
			   rho*(v.x-V.x)*v.x+p+(2.0/3.0)*dk(),
			   rho*(v.x-V.x)*v.y,
			   (v.x-V.x)*E()+ v.x*(p+(2.0/3.0)*dk()),
			   rho*(v.x-V.x)*k,
			   rho*(v.x-V.x)*omega);
}

/**********************************************************************
 * HighTemp2D_pState::dFdU -- Jacobian of the inviscid solution       *
 *                            flux with respect to the conserved      *
 *                            solution variables.                     *
 **********************************************************************/
inline void HighTemp2D_pState::dFdU(DenseMatrix &dFdU) const {
  switch (eos_type) {
    case EOS_TGAS: {
      double dhdpx; dhdpx = dhdp();
      double dhdrhox; dhdrhox = dhdrho();
      double hx; hx = Tgas_h(p, rho);
      dFdU(0,1) += 1.0;
      dFdU(1,0) += -(2.0*dhdpx*rho-3.0)/
                   (dhdpx*rho-1.0)*v.x*v.x/2.0+1/(dhdpx*rho-1.0)*v.y*v.y/2.0-
                   (hx+rho*dhdrhox)/(dhdpx*rho-1.0);
      dFdU(1,1) += (2.0*dhdpx*rho-3.0)/(dhdpx*rho-1.0)*v.x;
      dFdU(1,2) += -v.y/(dhdpx*rho-1.0);
      dFdU(1,3) += 1/(dhdpx*rho-1.0);
      dFdU(2,0) += -v.x*v.y;
      dFdU(2,1) += v.y;
      dFdU(2,2) += v.x;
      dFdU(3,0) += -(dhdpx*rho-2.0)/(dhdpx*rho-1.0)*v.x*v.x*v.x/2.0+
                  (-(dhdpx*rho-2.0)/(dhdpx*rho-1.0)*v.y*v.y/2.0-
                  rho*(dhdrhox+dhdpx*hx)/(dhdpx*rho-1.0))*v.x;
      dFdU(3,1) += (dhdpx*rho-3.0)/(dhdpx*rho-1.0)*v.x*v.x/2.0+hx+v.y*v.y/2.0;
      dFdU(3,2) += -v.y/(dhdpx*rho-1.0)*v.x;
      dFdU(3,3) += rho*v.x*dhdpx/(dhdpx*rho-1.0);
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
      break; }
    case EOS_IDEAL:
    default: {
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
      break; }
  }
}

#if defined(ALI_CHECK_HIGHTEMP) && defined(COMPILING_DRDU)
bool ali_dump_diffs_global = true;
int ali_dump_diffs_cpu = 0;
#endif

/**********************************************************************
 * HighTemp2D_pState::dFdW -- Jacobian of the inviscid solution       *
 *                            flux with respect to the primitive      *
 *                            solution variables.                     *
 **********************************************************************/
inline void HighTemp2D_pState::dFdW(DenseMatrix &dFdW) const {
  dFdW(0,0) += v.x;
  dFdW(0,1) += rho;
  dFdW(1,0) += v.x*v.x;
  dFdW(1,1) += TWO*rho*v.x; 
  dFdW(1,3) += ONE;
  dFdW(2,0) += v.x*v.y;
  dFdW(2,1) += rho*v.y;
  dFdW(2,2) += rho*v.x;
  switch (eos_type) {
    case EOS_TGAS: {
      double dhdpx; dhdpx = dhdp();
      double dhdrhox; dhdrhox = dhdrho();
      double hx; hx = Tgas_h(p, rho);
      dFdW(3,0) += v.x*v.x*v.x/2.0+(hx+v.y*v.y/2.0+rho*dhdrhox)*v.x;
      dFdW(3,1) += 3.0/2.0*rho*v.x*v.x+rho*v.y*v.y/2.0+rho*hx;
      dFdW(3,2) += rho*v.x*v.y;
      dFdW(3,3) += rho*v.x*dhdpx;
      break; }
    case EOS_IDEAL:
    default: {
      dFdW(3,0) += HALF*(v.x*v.x+v.y*v.y)*v.x;
      dFdW(3,1) += rho*v.x*v.x+ rho*(g*p/rho/gm1 + HALF*(v.x*v.x+v.y*v.y));
      dFdW(3,2) += rho*v.x*v.y;
      dFdW(3,3) += v.x*g/gm1;
      break; }
  }
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
 * HighTemp2D_pState::Gx, Gy -- Solution viscous fluxes.              *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_pState::Gx(const HighTemp2D_pState &dWdx) const {
  return HighTemp2D_cState(ZERO,tau.xx,tau.xy,-q.x+v.x*tau.xx+v.y*tau.xy+(mu()+sigma_k*muT())*k,
                           (mu()+sigma_k*muT())*k,(mu()+sigma_omega*muT())*omega);
}

inline HighTemp2D_cState HighTemp2D_pState::Gy(const HighTemp2D_pState &dWdy) const {
  return HighTemp2D_cState(ZERO,tau.xy,tau.yy,-q.y+v.x*tau.xy+v.y*tau.yy+(mu()+sigma_k*muT())*k,
                           (mu()+sigma_k*muT())*k,(mu()+sigma_omega*muT())*omega);
}

/**********************************************************************
 * HighTemp2D_pState::ComputeViscousTerms -- Compute viscous          *
 *                                           stress tensor and        *
 *                                           heat flux vector.        *
 **********************************************************************/
inline void HighTemp2D_pState::ComputeViscousTerms(const HighTemp2D_pState &dWdx,
						   const HighTemp2D_pState &dWdy,
						   const Vector2D &X,
						   int Axisymmetric) {
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

  double kap = kappa() + kappaT();
  switch (eos_type) {
    case EOS_TGAS: {
      double dTdp_local; dTdp_local = dTdp();
      double dTdrho_local; dTdrho_local = dTdrho();
      q.x = -kap*(dWdx.p*dTdp_local + dWdx.rho*dTdrho_local);
      q.y = -kap*(dWdy.p*dTdp_local + dWdy.rho*dTdrho_local);
      break; }
    case EOS_IDEAL:
    default: {
      q.x = -kap*(dWdx.p - (p/rho)*dWdx.rho)/(rho*R);
      q.y = -kap*(dWdy.p - (p/rho)*dWdy.rho)/(rho*R);
      break; }
  }
}

/**********************************************************************
 * HighTemp2D_pState::lambda_x -- Eigenvalue(s) (x-direction).        *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_pState::lambda_x(void) const {
  return HighTemp2D_pState(v.x-c(),v.x,v.x,v.x+c(),v.x,v.x);
}

inline HighTemp2D_pState HighTemp2D_pState::lambda_x(const Vector2D &V) const {
  return HighTemp2D_pState(v.x-V.x-c(),v.x-V.x,v.x-V.x,v.x-V.x+c(),v.x-V.x,v.x-V.x);
}

/**********************************************************************
 * HighTemp2D_pState::lambda_x -- Eigenvalue(s) (x-direction).        *
 * For High-Temperature Air, average State only                       *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_pState::lambda_x(double aAvg) const {
  return HighTemp2D_pState(v.x-aAvg,v.x,v.x,v.x+aAvg,v.x,v.x);
}

/**********************************************************************
 * HighTemp2D_pState::rp_x -- Primitive right eigenvector             *
 *                            (x-direction).                          *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_pState::rp_x(int index) const {
  assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
  switch(index) {
  case 1 :
    return HighTemp2D_pState(ONE,-c()/rho,ZERO,c2()-(2.0/3.0)*k,ZERO,ZERO);
  case 2 :
    return HighTemp2D_pState(ONE,ZERO,ZERO,-(2.0/3.0)*k,ZERO,ZERO);
  case 3 :
    return HighTemp2D_pState(ZERO,ZERO,ONE,ZERO,ZERO,ZERO);
  case 4 :
    return HighTemp2D_pState(ONE,c()/rho,ZERO,c2()-(2.0/3.0)*k,ZERO,ZERO);
  case 5 :
    return HighTemp2D_pState(ZERO,ZERO,ZERO,-(2.0/3.0)*rho,ONE,ZERO);
  case 6 :
  default:
    return HighTemp2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE);
  };
}

/**********************************************************************
 * HighTemp2D_pState::rc_x -- Conserved right eigenvector             *
 *                            (x-direction).                          *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_pState::rc_x(int index) const {
  assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
  switch(index) {
  case 1 :
    return HighTemp2D_cState(ONE,v.x-c(),v.y,h()-c()*v.x+(2.0/3.0)*k,k,omega);
  case 2 : {
    switch(eos_type) {
      case EOS_TGAS:
        return HighTemp2D_cState(ONE,v.x,v.y,h()+2.0*k/3.0-rho*c()*c()/dpde(),k,omega);
      case EOS_IDEAL:
        return HighTemp2D_cState(ONE,v.x,v.y,HALF*v.sqr()+gm1i*(g-5.0/3.0)*k,k,omega);
      default:
        return HighTemp2D_cState(ONE,v.x,v.y,HALF*v.sqr()+gm1i*(g-5.0/3.0)*k,k,omega);
    } }
  case 3 :
    return HighTemp2D_cState(ZERO,ZERO,rho,rho*v.y,ZERO,ZERO);
  case 4 :
    return HighTemp2D_cState(ONE,v.x+c(),v.y,h()+c()*v.x+(2.0/3.0)*k,k,omega);
  case 5 : {
    switch(eos_type){
      case EOS_TGAS:
        return HighTemp2D_cState(ZERO,ZERO,ZERO,rho*(-2.0*rho/(3.0*dpde())+ONE),rho,ZERO); 
      case EOS_IDEAL:
        return HighTemp2D_cState(ZERO,ZERO,ZERO,rho*gm1i*(g-5.0/3.0),rho,ZERO);
      default:
        return HighTemp2D_cState(ZERO,ZERO,ZERO,rho*gm1i*(g-5.0/3.0),rho,ZERO);
    } }
  case 6 :
  default:
    return HighTemp2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,rho);
  };
}

/**********************************************************************
 * HighTemp2D_pState::rc_x -- Conserved right eigenvector             *
 *                            (x-direction).                          *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_pState::rc_x(int index,double dpde,double dpdrho,double cAvg) const {
  assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
  switch(index) {
  case 1 :
    return HighTemp2D_cState(ONE,v.x-cAvg,v.y,h()-cAvg*v.x+(2.0/3.0)*k,k,omega);
  case 2 :
    return HighTemp2D_cState(ONE,v.x,v.y,h()+2.0*k/3.0-rho*cAvg*cAvg/dpde,k,omega);
  case 3 :
    return HighTemp2D_cState(ZERO,ZERO,rho,rho*v.y,ZERO,ZERO);
  case 4 :
    return HighTemp2D_cState(ONE,v.x+cAvg,v.y,h()+cAvg*v.x+(2.0/3.0)*k,k,omega);
  case 5 :
    return HighTemp2D_cState(ZERO,ZERO,ZERO,rho*(-2.0*rho/(3.0*dpde)+ONE),rho,ZERO);
  case 6 :
  default:
    return HighTemp2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,rho);
  };
}

/**********************************************************************
 * HighTemp2D_pState::lp_x -- Primitive left eigenvector              *
 *                            (x-direction).                          *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_pState::lp_x(int index) const {
  assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
  switch(index) {
  case 1 :
    return HighTemp2D_pState(k/(3.0*c2()),-HALF*rho/c(),ZERO,HALF/c2(),rho/(3.0*c2()),ZERO);
  case 2 :
    return HighTemp2D_pState(ONE-(2.0/3.0)*k/c2(),ZERO,ZERO,-ONE/c2(),-(2.0/3.0)*rho/c2(),ZERO);
  case 3 :
    return HighTemp2D_pState(ZERO,ZERO,ONE,ZERO,ZERO,ZERO);
  case 4 :
    return HighTemp2D_pState(k/(3.0*c2()),HALF*rho/c(),ZERO,HALF/c2(),rho/(3.0*c2()),ZERO);
  case 5 :
    return HighTemp2D_pState(ZERO,ZERO,ZERO,ZERO,ONE,ZERO);
  case 6 :
  default:
    return HighTemp2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE);
  };
}

/**********************************************************************
 * HighTemp2D_pState::lp_x -- Primitive left eigenvector              *
 *  ONLY 4 GLAISTER FLUX      (x-direction).                          *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_pState::lp_x(int index, double cAvg) const {
  assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
  switch(index) {
  case 1 :
    return HighTemp2D_pState(k/(3.0*cAvg*cAvg),-HALF*rho/cAvg,ZERO,HALF/(cAvg*cAvg),
                            rho/(3.0*cAvg*cAvg),ZERO);
  case 2 :
    return HighTemp2D_pState(ONE-(2.0/3.0)*k/(cAvg*cAvg),ZERO,ZERO,
                             -ONE/(cAvg*cAvg),-(2.0/3.0)*rho/(cAvg*cAvg),ZERO);
  case 3 :
    return HighTemp2D_pState(ZERO,ZERO,ONE,ZERO,ZERO,ZERO);
  case 4 :
    return HighTemp2D_pState(k/(3.0*cAvg*cAvg),HALF*rho/cAvg,ZERO,HALF/(cAvg*cAvg),
                             rho/(3.0*cAvg*cAvg),ZERO);
  case 5 :
    return HighTemp2D_pState(ZERO,ZERO,ZERO,ZERO,ONE,ZERO);
  case 6 :
  default:
    return HighTemp2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE);
  };
}

/**********************************************************************
 * HighTemp2D_pState::S -- Include all source term vectors and        *
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
 * HighTemp2D_pState::Si -- Inviscid axisymmetric source terms        *
 *                          and Jacobian.                             *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_pState::Si(const Vector2D &X) const {
  return HighTemp2D_cState(-rho*v.y/X.y,-rho*v.x*v.y/X.y,-rho*sqr(v.y)/X.y,
                           -v.y*(H()+(2.0/3.0)*dk())/X.y,-v.y*dk()/X.y,-v.y*domega()/X.y);
}

inline void HighTemp2D_pState::dSidU(DenseMatrix &dSidU, const Vector2D &X) const {
  // Appears to be incomplete!  CPTG.
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
 * HighTemp2D_pState::Sv -- Viscous axisymmetric flow source term     *
 *                          vector and Jacobian.                      *
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
 * HighTemp2D_pState::St -- Turbulent source term vector and          *
 *                          Jacobian.                                 *
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
			   alpha*(omega/max(k,TOLER))*
                           production-beta_omega(dWdx,dWdy)*rho*sqr(omega));
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
 * HighTemp2D_pState::RoeAverage -- Determine Roe-average state.      *
 **********************************************************************/
inline void HighTemp2D_pState::RoeAverage(const HighTemp2D_pState &Wl,
                                          const HighTemp2D_pState &Wr) {
  switch (eos_type) {
      case EOS_TGAS: {
      double hl, hr, el, er, srhol, srhor, ha, ea;
      el = Wl.e();
      er = Wr.e();
      hl = el + Wl.p/Wl.rho;
      hr = er + Wr.p/Wr.rho;
      srhol = sqrt(Wl.rho);
      srhor = sqrt(Wr.rho);
      rho = srhol*srhor;
      v.x = (srhol*Wl.v.x + srhor*Wr.v.x)/(srhol+srhor);
      v.y = (srhol*Wl.v.y + srhor*Wr.v.y)/(srhol+srhor);
      k = (srhol*Wl.k + srhor*Wr.k)/(srhol+srhor);
      omega = (srhol*Wl.omega + srhor*Wr.omega)/(srhol+srhor);
      ha = (srhol*hl + srhor*hr)/(srhol+srhor);
      ea = (srhol*el + srhor*er)/(srhol+srhor);
      p = rho*(ha - ea);
      //p = Tgas_p(ea, rho);
      break; }
    case EOS_IDEAL:
    default: {
      double hl, hr, srhol, srhor, aa2, ha;
      hl = Wl.h(); 
      hr = Wr.h();
      srhol = sqrt(Wl.rho); 
      srhor = sqrt(Wr.rho);
      rho = srhol*srhor;
      v.x = (srhol*Wl.v.x + srhor*Wr.v.x)/(srhol+srhor);
      v.y = (srhol*Wl.v.y + srhor*Wr.v.y)/(srhol+srhor);
      k = (srhol*Wl.k + srhor*Wr.k)/(srhol+srhor);
      omega = (srhol*Wl.omega + srhor*Wr.omega)/(srhol+srhor);
      ha  = (srhol*hl + srhor*hr)/(srhol+srhor);
      aa2 = gm1*(ha - HALF*(sqr(v.x) + sqr(v.y)) - k);
      p = rho*aa2/g;   
      break; }
  }
}

inline void HighTemp2D_pState::RoeAverage(HighTemp2D_pState &Wa,
                                          const HighTemp2D_pState &Wl,
                                          const HighTemp2D_pState &Wr) const {
  Wa.RoeAverage(Wl, Wr);
}


/**********************************************************************
 * HighTemp2D_pState::RoeAverage_SoundSpeed -- Determines the         *
 *      Roe-average sound speed.  Also returns dpdrho & dpde.         *
 **********************************************************************/
inline void HighTemp2D_pState::RoeAverage_SoundSpeed(double &c_a,
                                                     double &dpdrho_a,
                                                     double &dpde_a,
                                                     const HighTemp2D_pState &Wa,
                                                     const HighTemp2D_pState &Wl,
                                                     const HighTemp2D_pState &Wr) const {
  switch (eos_type) {
    case EOS_TGAS: {
      int no_density_change, no_energy_change;
      double el, er, ea, 
             p2l, p1l, p2r, p1r, perdl, peldr,
  	     de, drho;
      el = Wl.e();
      er = Wr.e();
      ea = Wa.e();
      drho = (Wr.rho-Wl.rho);
      de = (er-el);
      if (abs(drho)<(TWO*HTTOL*Wa.rho))
        no_density_change = 1;
      else no_density_change = 0;
      if (abs(de)<(TWO*HTTOL*ea))
        no_energy_change = 1;
      else no_energy_change = 0;
      if (!no_density_change && !no_energy_change) {
         perdl = Tgas_p(er, Wl.rho);
         peldr = Tgas_p(el, Wr.rho);
         dpdrho_a = HALF*(Wr.p + peldr - perdl - Wl.p)/drho;
         dpde_a = HALF*(Wr.p + perdl - peldr - Wl.p)/de;
      } else if (no_density_change &&  no_energy_change) {
         p2l = Tgas_p(ea, Wa.rho*(ONE+HTTOL));
         p1l = Tgas_p(ea, Wa.rho*(ONE-HTTOL)); 
         dpdrho_a = (p2l-p1l)/(TWO*HTTOL*Wa.rho);
         p2l = Tgas_p(ea*(ONE+HTTOL), Wa.rho);
         p1l = Tgas_p(ea*(ONE-HTTOL), Wa.rho);    
         dpde_a = (p2l-p1l)/(TWO*HTTOL*ea);
      } else if (no_density_change) {
         perdl = Tgas_p(er, Wl.rho);
         peldr = Tgas_p(el, Wr.rho);
         dpde_a = HALF*(Wr.p + perdl - peldr - Wl.p)/de;
         p2l = Tgas_p(el, Wa.rho*(ONE+HTTOL));
         p1l = Tgas_p(el, Wa.rho*(ONE-HTTOL)); 
         p2r = Tgas_p(er, Wa.rho*(ONE+HTTOL));
         p1r = Tgas_p(er, Wa.rho*(ONE-HTTOL));
         dpdrho_a = HALF*((p2l-p1l)/(TWO*HTTOL*Wa.rho)+
		          (p2r-p1r)/(TWO*HTTOL*Wa.rho)); 
      } else {
         perdl = Tgas_p(er, Wl.rho);
         peldr = Tgas_p(el, Wr.rho);
         dpdrho_a = HALF*(Wr.p + peldr - perdl - Wl.p)/drho;
         p2l = Tgas_p(ea*(ONE+HTTOL), Wl.rho);
         p1l = Tgas_p(ea*(ONE-HTTOL), Wl.rho); 
         p2r = Tgas_p(ea*(ONE+HTTOL), Wr.rho);
         p1r = Tgas_p(ea*(ONE-HTTOL), Wr.rho);
         dpde_a = HALF*((p2l-p1l)/(TWO*HTTOL*ea)+
		        (p2r-p1r)/(TWO*HTTOL*ea));
      } /* endif */
      c_a = sqrt(Wa.p*dpde_a/(Wa.rho*Wa.rho)+dpdrho_a+
                 TWO*Wa.k/THREE + (TWO*Wa.k*dpde_a)/(THREE*Wa.rho));
      break; }
    case EOS_IDEAL:
    default: {
      c_a = Wa.c();
      dpdrho_a = Wa.dpdrho();
      dpde_a = Wa.dpde();
      break; }
  }
}

/**********************************************************************
 * HighTemp2D_pState::Rotate -- Rotates solution state.               *
 **********************************************************************/
inline void HighTemp2D_pState::Rotate(const HighTemp2D_pState &W, 
                                      const Vector2D &norm_dir) {
  rho = W.rho;
  v.x =  W.v.x*norm_dir.x+W.v.y*norm_dir.y;
  v.y = -W.v.x*norm_dir.y+W.v.y*norm_dir.x;
  p = W.p;
  k = W.k;
  omega = W.omega;
}

/**********************************************************************
 * HighTemp2D_cState::HighTemp2D_cState -- Constructor.               *
 **********************************************************************/
inline HighTemp2D_cState::HighTemp2D_cState(const HighTemp2D_pState &W) {
  rho = W.rho; dv = W.dv(); E = W.E(); dk = W.dk(); domega = W.domega();
}

/**********************************************************************
 * HighTemp2D_cState::W -- Primitive solution state.                  *
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
 * HighTemp2D_cState::dUdW -- Jacobian of the conserved solution      *
 *                            variables with respect to the           *
 *                            primitive solution variables.           *
 **********************************************************************/
inline void HighTemp2D_cState::dUdW(DenseMatrix &dUdW) const {
  switch (eos_type) {
    case EOS_TGAS: {
      double dpdex; dpdex = dpde();
      double dpdrhox; dpdrhox = dpdrho();
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
      break; }
    case EOS_IDEAL:
    default: {
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
      break; }
  }
}

/**********************************************************************
 * HighTemp2D_cState::dWdU -- Jacobian of the primitive solution      *
 *                            variables with respect to the conserved *
 *                            solution variables.                     *
 **********************************************************************/
inline void HighTemp2D_cState::dWdU(DenseMatrix &dWdU) const {
  switch (eos_type) {
    case EOS_TGAS: {
      double dpdex; dpdex = dpde();
      double dpdrhox; dpdrhox = dpdrho();
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
      break; }
    case EOS_IDEAL:
    default: {
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
      break; }
  }
}

/**********************************************************************
 * HighTemp2D_cState::F -- Solution inviscid flux (x-direction).  *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_cState::Fx(void) const {
  return HighTemp2D_cState(dv.x,sqr(dv.x)/rho+p()+(2.0/3.0)*dk,dv.x*dv.y/rho,
                           dv.x*H()/rho+v().x*(2.0/3.0)*dk,v().x*dk,v().x*domega);
}

inline HighTemp2D_cState HighTemp2D_cState::Fx(const Vector2D &V) const {
  double vx = v().x;
  return HighTemp2D_cState(rho*(vx-V.x),(vx-V.x)*dv.x+p()+(2.0/3.0)*dk,
			       (vx-V.x)*dv.y,(vx-V.x)*E+ vx*(p()+(2.0/3.0)*dk),
			       (vx-V.x)*dk,(vx-V.x)*domega);
}

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
 * HighTemp2D_cState::Gx, Gy -- Solution viscous fluxes.              *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_cState::Gx(const HighTemp2D_pState &dWdx) const {
  return HighTemp2D_cState(ZERO,tau.xx,tau.xy,
                           -q.x+v().x*tau.xx+v().y*tau.xy+(mu()+sigma_k*muT())*dWdx.k,
                           (mu()+sigma_k*muT())*dWdx.k,
                           (mu()+sigma_omega*muT())*dWdx.omega);
}

inline HighTemp2D_cState HighTemp2D_cState::Gy(const HighTemp2D_pState &dWdy) const {
  return HighTemp2D_cState(ZERO,tau.xy,tau.yy,
                           -q.y+v().x*tau.xy+v().y*tau.yy+(mu()+sigma_k*muT())*dWdy.k,
                           (mu()+sigma_k*muT())*dWdy.k,
                           (mu()+sigma_omega*muT())*dWdy.omega);
}

/**********************************************************************
 * HighTemp2D_cState::ComputeViscousTerms -- Compute viscous          *
 *                                           stress tensor and        *
 *                                           heat flux vector.        *
 **********************************************************************/
inline void HighTemp2D_cState::ComputeViscousTerms(const HighTemp2D_pState &dWdx,
						   const HighTemp2D_pState &dWdy,
						   const Vector2D &X,
						   int Axisymmetric) {
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

  double kap = kappa() + kappaT();
  switch (eos_type) {
    case EOS_TGAS: {
      double dTdp_local; dTdp_local = dTdp();
      double dTdrho_local; dTdrho_local = dTdrho();
      q.x = -kap*(dWdx.p*dTdp_local + dWdx.rho*dTdrho_local);
      q.y = -kap*(dWdy.p*dTdp_local + dWdy.rho*dTdrho_local);
      break; }
    case EOS_IDEAL:
    default: {
      q.x = -kap*(dWdx.p - (p()/rho)*dWdx.rho)/(rho*R);
      q.y = -kap*(dWdy.p - (p()/rho)*dWdy.rho)/(rho*R);
      break; }
  }
}

/**********************************************************************
 * HighTemp2D_cState::lambda_x -- Eigenvalue(s) (x-direction).        *
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
 * HighTemp2D_pState::lambda_x -- Eigenvalue(s) (x-direction).        *
 * For High-Temperature Air, average State only                       *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_cState::lambda_x(double aAvg) const {
  double vx = v().x;
  return HighTemp2D_pState(vx-aAvg,vx,vx,vx+aAvg,vx,vx);
}

/**********************************************************************
 * HighTemp2D_pState::rp_x -- Primitive right eigenvector             *
 *                            (x-direction).                          *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_cState::rp_x(int index) const {
  assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
  switch(index) {
  case 1 :
    return HighTemp2D_pState(ONE,-c()/rho,ZERO,c2()-(2.0/3.0)*k(),ZERO,ZERO);
  case 2 :
    return HighTemp2D_pState(ONE,ZERO,ZERO,-(2.0/3.0)*k(),ZERO,ZERO);
  case 3 :
    return HighTemp2D_pState(ZERO,ZERO,ONE,ZERO,ZERO,ZERO);
  case 4 :
    return HighTemp2D_pState(ONE,c()/rho,ZERO,c2()-(2.0/3.0)*k(),ZERO,ZERO);
  case 5 :
    return HighTemp2D_pState(ZERO,ZERO,ZERO,-(2.0/3.0)*rho,ONE,ZERO);
  case 6 :
  default:
    return HighTemp2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE);
  };
}

/**********************************************************************
 * HighTemp2D_cState::rc_x -- Conserved right eigenvector             *
 *                            (x-direction).                          *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_cState::rc_x(int index) const {
  assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
  switch(index) {
    case 1 :
      return HighTemp2D_cState(ONE,v().x-c(),v().y,h()-c()*v().x+(2.0/3.0)*k(),k(),omega());
    case 2 : {
      switch(eos_type) {
        case EOS_TGAS:
          return HighTemp2D_cState(ONE,v().x,v().y,h()+2.0*k()/3.0-rho*c()*c()/dpde(),
                                   k(),omega()); 
        case EOS_IDEAL:
        default:
          return HighTemp2D_cState(ONE,v().x,v().y,HALF*v().sqr()+gm1i*(g-5.0/3.0)*k(),
                                   k(),omega());
      } }
    case 3 :
      return HighTemp2D_cState(ZERO,ZERO,rho,dv.y,ZERO,ZERO);
    case 4 :
      return HighTemp2D_cState(ONE,v().x+c(),v().y,h()+c()*v().x+(2.0/3.0)*k(),k(),omega());
    case 5 : {
      switch(eos_type) {
        case EOS_TGAS:
          return HighTemp2D_cState(ZERO,ZERO,ZERO,rho*(-2.0*rho/(3.0*dpde())+ONE),rho,ZERO); 
        case EOS_IDEAL:
        default:
          return HighTemp2D_cState(ZERO,ZERO,ZERO,rho*gm1i*(g-5.0/3.0),rho,ZERO);
      } }
    case 6 :
    default:
      return HighTemp2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,rho);
  };
}

/**********************************************************************
 * HighTemp2D_cState::rc_x -- Conserved right eigenvector             *
 *           GLAISTER FLUX ONLY!     (x-direction).                   *
 **********************************************************************/
inline HighTemp2D_cState HighTemp2D_cState::rc_x(int index, double dpde, double dpdrho, double cAvg) const {
  assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
  switch(index) {
    case 1 :
      return HighTemp2D_cState(ONE,v().x-cAvg,v().y,h()-cAvg*v().x+(2.0/3.0)*k(),k(),omega());
    case 2 :
      return HighTemp2D_cState(ONE,v().x,v().y,h()+2.0*k()/3.0-rho*cAvg*cAvg/dpde,k(),omega());
    case 3 :
      return HighTemp2D_cState(ZERO,ZERO,rho,rho*v().y,ZERO,ZERO);
    case 4 :
      return HighTemp2D_cState(ONE,v().x+cAvg,v().y,h()+cAvg*v().x+(2.0/3.0)*k(),k(),omega());
    case 5 :
      return HighTemp2D_cState(ZERO,ZERO,ZERO,rho*(-2.0*rho/(3.0*dpde)+ONE),rho,ZERO);
    case 6 :
    default:
      return HighTemp2D_cState(ZERO,ZERO,ZERO,ZERO,ZERO,rho);
  };
}

/**********************************************************************
 * HighTemp2D_cState::lp_x -- Primitive left eigenvector              *
 *                            (x-direction).                          *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_cState::lp_x(int index) const {
  assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
  switch(index) {
    case 1 :
      return HighTemp2D_pState(k()/(3.0*c2()),-HALF*rho/c(),ZERO,HALF/c2(),
                               rho/(3.0*c2()),ZERO);
    case 2 :
      return HighTemp2D_pState(ONE-(2.0/3.0)*k()/c2(),ZERO,ZERO,-ONE/c2(),
                               -(2.0/3.0)*rho/c2(),ZERO);
    case 3 :
      return HighTemp2D_pState(ZERO,ZERO,ONE,ZERO,ZERO,ZERO);
    case 4 :
      return HighTemp2D_pState(k()/(3.0*c2()),HALF*rho/c(),ZERO,HALF/c2(),
                               rho/(3.0*c2()),ZERO);
    case 5 :
      return HighTemp2D_pState(ZERO,ZERO,ZERO,ZERO,ONE,ZERO);
    case 6 :
    default:
      return HighTemp2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE);
  };
}

/**********************************************************************
 * HighTemp2D_pState::lp_x -- Primitive left eigenvector              *
 *  ONLY 4 GLAISTER FLUX      (x-direction).                          *
 **********************************************************************/
inline HighTemp2D_pState HighTemp2D_cState::lp_x(int index, double cAvg) const {
  assert(index >= 1 && index <= NUM_VAR_HIGHTEMP2D);
  switch(index) {
    case 1 :
      return HighTemp2D_pState(k()/(3.0*cAvg*cAvg),-HALF*rho/cAvg,ZERO,HALF/(cAvg*cAvg),
                               rho/(3.0*cAvg*cAvg),ZERO);
    case 2 :
      return HighTemp2D_pState(ONE-(2.0/3.0)*k()/(cAvg*cAvg),ZERO,ZERO,-ONE/(cAvg*cAvg),                                          -(2.0/3.0)*rho/(cAvg*cAvg),ZERO);
    case 3 :
      return HighTemp2D_pState(ZERO,ZERO,ONE,ZERO,ZERO,ZERO);
    case 4 :
      return HighTemp2D_pState(k()/(3.0*cAvg*cAvg),HALF*rho/cAvg,ZERO,HALF/(cAvg*cAvg),
                               rho/(3.0*cAvg*cAvg),ZERO);
    case 5 :
      return HighTemp2D_pState(ZERO,ZERO,ZERO,ZERO,ONE,ZERO);
    case 6 :
    default:
      return HighTemp2D_pState(ZERO,ZERO,ZERO,ZERO,ZERO,ONE);
  };
}

/**********************************************************************
 * Useful 2D HighTemp state constants.                                *
 **********************************************************************/
const HighTemp2D_pState HighTemp2D_W_STDATM(DENSITY_STDATM,
                                            Vector2D_ZERO, PRESSURE_STDATM);

/**********************************************************************
 * HighTemp2DState -- External subroutines.                           *
 **********************************************************************/

extern HighTemp2D_pState Riemann(const HighTemp2D_pState &Wl,
                                 const HighTemp2D_pState &Wr);

extern HighTemp2D_pState RoeAverage(const HighTemp2D_pState &Wl,
			            const HighTemp2D_pState &Wr);

extern HighTemp2D_pState Translate(const HighTemp2D_pState &W,
			           const Vector2D &V);

extern HighTemp2D_pState Reflect(const HighTemp2D_pState &W,
			         const Vector2D &norm_dir);

extern HighTemp2D_pState Rotate(const HighTemp2D_pState &W, 
                                const Vector2D &norm_dir);

extern HighTemp2D_pState Mirror(const HighTemp2D_pState &W,
			        const Vector2D &norm_dir);

extern HighTemp2D_pState WallViscousHeatFlux(const HighTemp2D_pState &W,
				             const Vector2D &norm_dir);

extern HighTemp2D_pState WallViscousIsothermal(const HighTemp2D_pState &W,
					       const Vector2D &norm_dir,
					       const double &Twall);

extern HighTemp2D_pState BurningSurface(const HighTemp2D_pState &W,
                                        const Vector2D &norm_dir);

extern HighTemp2D_pState MovingWallHeatFlux(const HighTemp2D_pState &W,
					    const Vector2D &norm_dir,
					    const double &Vwall);

extern HighTemp2D_pState MovingWallIsothermal(const HighTemp2D_pState &W,
					      const Vector2D &norm_dir,
				              const double &Vwall,
				              const double &Twall);

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

extern Vector2D HLLE_wavespeeds(const HighTemp2D_pState &Wl,
		                const HighTemp2D_pState &Wr,
		                const Vector2D &norm_dir);

extern Vector2D GHLLE_wavespeeds(const HighTemp2D_pState &Wl,
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
