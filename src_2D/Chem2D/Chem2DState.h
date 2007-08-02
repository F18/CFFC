/****************** Chem2DState.h **************************************
  This class defines the state variables and constructors for the 
  2D Axisymmertric Navier-Stokes with Multiple Species solution.

  TODO:  
         - be careful with set_species_data() and set_initial_values(),
           due to dynamic memory 
    
         - Avoid cState.T() as much as possible as it is $$$ to calculate, pass
           the temperature from pState if possible as it is cheap (p/rhoR)
             
              
   NOTES: - The Static "specdata" has to be deleted outside the class
            in order to avoid seg faults. In Chem2D it is done by calling
            the deallocate function from the Input destructor as it is kept
            until the end of the program.
                   
***********************************************************************/

#ifndef _CHEM2D_STATE_INCLUDED 
#define _CHEM2D_STATE_INCLUDED

// define early so can be used in other classes
// yeah I know its a fudge, but it works...
class Chem2D_cState;
class Chem2D_pState;

// Required C++ libraries
#include <iostream>

using namespace std;

// Required CFDkit+caboodle header files

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _MATRIX_INCLUDED
#include "../Math/Matrix.h"
#endif // _MATRIX_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _TENSOR2D_INCLUDED
#include "../Math/Tensor2D.h"
#endif //_TENSOR2D_INCLUDED

#ifndef _VECTOR2D_INCLUDED
#include "../Math/Vector2D.h"
#endif //_VECTOR2D_INCLUDED

#ifndef _GAS_CONSTANTS_INCLUDED
#include "../Physics/GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED

// CHEM2D Specific headers
#ifndef _SPECIES_INCLUDED
#include "Species.h"
#endif //_SPECIES_INCLUDED

#ifndef _NASARP1311_DATA_INCLUDED
#include "NASARP1311data.h"
#endif

#ifndef _REACTIONS_INCLUDED
#include "Reactions.h"
#endif 

//Temperature convergence tolerance in
//Chem2D_cState::T(void)
// these should be moved to CFD.h or Math.h
#define CONV_TOLERANCE 1e-8
#define gravity_z -9.81  //m/s^2 acceleration due to gravity

//number of fixed variables in the Chem2D class
#define NUM_CHEM2D_VAR_SANS_SPECIES 6  //rho, v(2), p, k and omega

/*!
 * Class: Chem2D_pState
 *
 * @brief  Primitive variable solution state class definition for an
 *         inviscid, laminar, or turbulent chemically reacting gas-flow.
 *
 * Primitive variable solution state class definition for an inviscid, 
 * laminar, or turbulent chemically reacting gas-flow.
 *
 * \verbatim
 * Density:   rho  kg/m^3
 * Velocity:  v    m/s
 * Pressure:  p    Pa (N/m^2)
 * 
 * Molecular Mass:          M    kg/mol
 * Species Gas Constant:    Rs   J/(kg*K)
 * 
 * Temperature(Kelvin) Dependent data: f(T) 
 *  
 * Heat Capacity (const Pressure):  Cp  J/(kg*K)
 * Heat Capacity (const Volume):    Cv  J/(kg*K)
 * Specific Heat Ratio:             g
 * Specific Enthalpy:               h   J/kg
 * Specific Internal Energy:        e   J/kg
 * Total Enthalpy:                  H   J/kg 
 * Total Internal Energy:           E   J/kg
 *
 * Viscosity:                       mu  kg/(m*s) N*s/m^2  
 * Thermal Conductivity:            k   N/(s*K)  W.(m*K)
 * \endverbatim
 */
class Chem2D_pState {
private: 
protected:
public:
  //@{ @name Primitive variables and associated constants:
  double                      rho; //!< Density.
  Vector2D                      v; //!< Flow velocity (2D)  
  double                        p; //!< Pressure.
  Species                   *spec; //!< Species class c[ns]
  double                        k; //!< Turbulent kinetic energy.
  double                    omega; //!< Turbulent specific dissipation rate.
  Tensor2D                    tau; //!< Shear Stress Tensor
  Vector2D                  qflux; //!< Heat Flux Vector  
  Tensor2D                 lambda; //!< Reynolds Stress Tensor
  Vector2D                  theta; //!< Turbulent Heat Flux Vector  
  static int                   ns; //!< number of species
  static NASARP1311data *specdata; //!< Global Species Data
  static double          *Schmidt; //!< Schmidt Number for each species
  static Reaction_set       React; //!< Global Reaction Data
  static double    low_temp_range; //!< Low temp data range
  static double   high_temp_range; //!< High temp data range
  static int       NUM_VAR_CHEM2D; //!< Number of Chem2d variables (6+ns)
  static int            flow_type; //!< Flow-type indicator (inviscid, laminar, or k-omega turbulent).
  static double              Mref; //!< Mref for Precondtioning (normally set to incoming freestream Mach)
  static int          debug_level; //!< Debug level flag (0=none, 1,2,..)
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
  static double              alpha;
  static double              sigma;
  static double         sigma_star;
  static double               beta;
  static double             f_beta;
  static double          beta_star;
  static double        f_beta_star;
  static double          Coeff_edm;
  //@}
 
  //@{ @name Creation, copy, and assignment constructors.
  Chem2D_pState(){rho = DENSITY_STDATM; v.zero(); p = PRESSURE_STDATM; 
                  k = ZERO; omega = ZERO; tau.zero(); qflux.zero(); 
                  lambda.zero(); theta.zero(); spec = NULL;  set_initial_values(); }

  Chem2D_pState(const double &value)
               {rho =value; v.x=value, v.y=value; p = value; 
                k = ZERO; omega = ZERO; tau.zero(); qflux.zero(); lambda.zero(); theta.zero(); 
                spec = NULL;  set_initial_values(value); }

  Chem2D_pState(const double &d, const Vector2D &V, const double &pre)
                {rho = d; v = V; p = pre;  k = ZERO; omega = ZERO; 
                 tau.zero(); qflux.zero(); lambda.zero(); theta.zero();
		 spec = NULL;  set_initial_values(); }

  Chem2D_pState(const double &d, const double &vx, const double &vy, const double &pre)
                {rho = d; v.x = vx; v.y = vy; p = pre; k = ZERO; omega = ZERO; 
                 tau.zero(); qflux.zero(); lambda.zero(); theta.zero();
		 spec = NULL;  set_initial_values(); }

  Chem2D_pState(const double &d, const double &vx, const double &vy, const double &pre, const double &value)
                {rho = d; v.x=vx; v.y=vy; p = pre; k = ZERO; omega = ZERO; 
                 tau.zero(); qflux.zero(); lambda.zero(); theta.zero();
		 spec = NULL;  set_initial_values(value); }

  Chem2D_pState(const double &d, const double &vx, const double &vy, const double &pre, 
                const double &kk, const double &oomega,	const double &value)
                {rho = d; v.x=vx; v.y=vy; p = pre; k = kk; omega = oomega; 
                 tau.zero(); qflux.zero(); lambda.zero(); theta.zero();
		 spec = NULL;  set_initial_values(value); }

  Chem2D_pState(const double &d, const double &vx, const double &vy, const double &pre, 
                const double &kk, const double &oomega) 
                {rho = d; v.x=vx; v.y=vy; p = pre; k = kk; omega = oomega; tau.zero(); qflux.zero();
		lambda.zero(); theta.zero(); spec = NULL;  set_initial_values(); }

  Chem2D_pState(const double &d, const Vector2D &V, const double &pre, const double &kk, const double & oomega)
                {rho = d; v = V; p = pre; k = kk; omega = oomega; tau.zero(); qflux.zero();
		lambda.zero(); theta.zero(); spec = NULL;  set_initial_values(); } 

  Chem2D_pState(const double &d, const double &vx, const double &vy, const double &pre, 
                const double &kk, const double & oomega, Species *mfrac) 
                {rho = d; v.x=vx; v.y=vy; p = pre; k = kk; omega = oomega; 
                 tau.zero(); qflux.zero(); lambda.zero(); theta.zero(); 
		 spec = NULL;  set_initial_values(mfrac); }

  Chem2D_pState(const double &d, const Vector2D &V, const double &pre,
                const double &kk, const double & oomega, Species *mfrac) 
                {rho = d; v = V; p = pre; k = kk; omega = oomega; tau.zero(); qflux.zero();
		 lambda.zero(); theta.zero(); spec = NULL; set_initial_values(mfrac); }

  //this is needed for the operator overload returns!!!!
  Chem2D_pState(const Chem2D_pState &W)
                {spec = NULL; rho = DENSITY_STDATM; set_initial_values(); Copy(W);}

  //! Default destructor.
  ~Chem2D_pState(){ Deallocate(); }

  //! Deallocation for static specdata.
  void Deallocate_static(void){ if(specdata != NULL) delete[] specdata; 
                                          specdata = NULL; 
				if(Schmidt != NULL) delete[] Schmidt; 
				           Schmidt = NULL;
                              }

  void Deallocate(void){ if(spec != NULL) delete[] spec; spec = NULL; }
  //@}

  //@{ @name Set static variables.
  //! Read in ns species data, call only once as its static
  void set_species_data(const int &, const string *, const char *,
			const int &, const double&, const double *);

  //! Set initial data values predominately used internally !!!
  void set_initial_values();
  void set_initial_values(const double &value);
  void set_initial_values(double *cfrac);
  void set_initial_values(Species *mfrac);

  //! Set flow type static variable.
  void set_flow_type(const int &Flow_Type) { flow_type = Flow_Type; }

  //! Set turbulence static variables.
  void set_turbulence_variables(const double &C_constant,
				const double &von_karman_constant,
				const double &yplus_sub,
				const double &yplus_buffer,
				const double &yplus_oute);
  //@}

  //@{ @name Useful operators.
  //! Copy operator.
  void Copy(const Chem2D_pState &W);

  //! Vacuum operator.
  void Vacuum(){ rho=ZERO; v.x=ZERO; v.y=ZERO; p=ZERO; k=ZERO; omega = ZERO; 
    for(int i=0; i<ns; i++){
      spec[i].Vacuum();
    }
    tau.zero();  qflux.zero();
    lambda.zero(); theta.zero(); 

  }

  //! Zero non-solution quantities.
  void zero_non_sol(){
    for(int i = 0; i < ns; i++){
      spec[i].gradc.zero();
      spec[i].diffusion_coef=ZERO;
    }
    tau.zero();  qflux.zero();
    lambda.zero(); theta.zero();
  }

  //! Set the lowest value of the temperature range.
  void Temp_low_range();

  //! Set the highest value of the temperature range.
  void Temp_high_range();

  //! Check for unphysical state properties.
  int Unphysical_Properties(void) const;
  //@}

  //@{ @name Mixing Rules.

  // The following constructors return "total" physical
  // parameters based on mixture rules for each.

  //! Mixture molecular mass.
  double Mass(void) const;
  //! Mixture gas constant.
  double Rtot(void);
  double Rtot(void) const;
  //! Mixture heat capacity (Pressure constant)
  double Cp(void) const;
  double Cp(const double &TEMP) const;
  //! Mixture heat capacity (Volume constant)
  double Cv(void) const;
  //! Mixture heat capacity ratio
  double g(void) const;
  //! Mixture absolute (sensible+chemical) internal energy
  double e(void) const;
  double eref(void) const;
  //! Mixture sensible internal energy.
  double es(void) const;
  //! Mixture specific enthalpy.
  double h(void) const;
  double h(const double &T) const;
  double href(void) const;
  double hs(void) const;
  double hs(const double &T) const;
  //! Mixture total internal energy.
  double E(void) const;
  //! Mixture total enthalpy.
  double H(void) const;
  //! Mixture total enthalpy.
  double Hs(void) const;
  //! Mixture viscosity.
  double mu(void) const;
  //! Mixture thermal conductivity.
  double kappa(void) const;
  double hprime(void) const;  
  double hprime(double &Temp) const;

  // Mixture diffusion coefficient 
  //double Diffusion_coef() const;
  //@}

  //@{ @name State functions.
  //! Momentum.
  Vector2D rhov(void) const;
  //! Temperature.
  double T(void) const;
  //! Determine temperature given the sensible enthalpy.
  double T(double &h_s) const;
  //! Initial gamma for iterative schemes.
  double gamma_guess(void) const;
  //! Speed of sound.
  double a(void);
  double a(void) const;
  //! Check for -ve mass fractions and sets small -ve c's to ZERO.
  bool negative_speccheck(void) const;
  //! Check for unphysical state properties.
  bool Unphysical_Properties_Check(Chem2D_cState &U, const int n) const;
  //! Species i concentration (rho*c/mol_mass).
  double SpecCon(int i) const;
  //! Gibbs Free Energy (H-TS) for species.
  double Gibbs(int species) const;
  //@}

  //@{ @name Turbulence related functions.
  double eddy_viscosity(void) const;      
  double Pr_turb(void) const;      
  double Sc_turb(void) const;      
  double Kappa_turb(void) const;      
  double Dm_turb(void) const;
  double omega_sublayer_BC(const double &y) const;
  //@}

  //@{ @name Dimensionless Parameters.
  double Schmidt_No(const int &) const;
  double Prandtl() const;
  double Lewis(const int &) const;
  //@}

  //@{ @name Temperature Derivatives.
  double diedip() const;
  double diedirho() const;
  double dmudT(void) const;
  //@}

  //@{ @name Strain rate, laminar stress, and Reynolds stress tensors.
  Tensor2D Strain_Rate(const Chem2D_pState &dWdx,
		       const Chem2D_pState &dWdy,
		       const int Axisymmetric,
		       const Vector2D X);

  Tensor2D Laminar_Stress(const Chem2D_pState &dWdx,
			  const Chem2D_pState &dWdy,
			  const int Axisymmetric,
			  const Vector2D X);

  Tensor2D Reynolds_Stress(const Chem2D_pState &dWdx,
			   const Chem2D_pState &dWdy,
			   const int Axisymmetric,
			   const Vector2D X);

  Tensor2D Strain_Rate(const Chem2D_pState &W,
		       const Chem2D_pState &dWdx,
		       const Chem2D_pState &dWdy,
		       const int Axisymmetric,
		       const Vector2D X);

  Tensor2D Laminar_Stress(const Chem2D_pState &W,
			  const Chem2D_pState &dWdx,
			  const Chem2D_pState &dWdy,
			  const int Axisymmetric,
			  const Vector2D X);

  Tensor2D Reynolds_Stress(const Chem2D_pState &W,
			   const Chem2D_pState &dWdx,
			   const Chem2D_pState &dWdy,
			   const int Axisymmetric,
			   const Vector2D X);
  //@}

  //@{ @name Heat Flux vector thermal Diffusion.
  Vector2D thermal_diffusion(void) const;
  //@}

  //@{ @name Conserved solution state.
  Chem2D_cState U(void) const;
  Chem2D_cState U(const Chem2D_pState &W) const;
  friend Chem2D_cState U(const Chem2D_pState &W);
  //@}

  //@{ @name Fluxes.
  //! Inviscid Solution flux (x-direction).
  Chem2D_cState Fx(void) const;
  
  //! For rotated frame.
  void dFIdU(DenseMatrix &dFdU);
  void dFIdU(DenseMatrix &dFdU) const;
  friend void dFIdU(DenseMatrix &dFdU, const Chem2D_pState &W);
  
  //! Compute the dWdU.
  void dWdU(DenseMatrix &dWdQ);
  //@}

  //@{ @name Eigenstructure.
  //! Eigenvalue(s) (x-direction).
  Chem2D_pState lambda_x(void) const;

  //! Precondtioned Eigenvalues.
  Chem2D_pState lambda_preconditioned_x(const double &MR2) const;

  //! Conserved right eigenvectors (x-direction).
  Chem2D_cState rc_x(const int &index) const;
  //! Primitive left eigenvectors (x-direction).
  Chem2D_pState lp_x(const int &index) const;
   
  //! Conserved right preconditioned eigenvectors (x-direction).
  Chem2D_cState rc_x_precon(const int &index,const double &MR2) const;
  //! Primitive left preconditioned eigenvectors (x-direction).
  Chem2D_pState lp_x_precon(const int &index,const double &MR2) const;
  //@}
  
  //@{ @name Preconditioner.
  double u_plus_aprecon(const double &u,
			const double &deltax) const;
  void u_a_precon(const double &UR,double &uprimed, double &cprimed) const;
  double Mr2(const double &deltax) const;
  void Low_Mach_Number_Preconditioner(DenseMatrix &P,
				      const double &deltax) const; 
  void Low_Mach_Number_Preconditioner_Inverse(DenseMatrix &Pinv,
					      const double &deltax) const; 
  //@}

  //@{ @name Axisymmetric Source Terms.
  //! Inviscid axisymmetric coordinate system source term.
  Chem2D_cState Sa_inviscid(const Vector2D &X,
                            const int Axisymmetric) const;

  //! Vicous axisymmetric coordinate system source term.
  Chem2D_cState Sa_viscous(const Chem2D_pState &dWdx, 
                           const Chem2D_pState &dWdy, 
                           const Vector2D &X,
                           const int Axisymmetric);
  //! Inviscid axisymmetric coordinate system Jacobian.
  void dSa_idU(DenseMatrix &dSa_idU,
	       const Vector2D &X,
	       const int Axisymmetric) const;
  //! Vicous axisymmetric coordinate system Jacobian.
  void dSa_vdU(DenseMatrix &dSa_VdU,
	       DenseMatrix &dWdQ,
	       const Chem2D_pState &dWdx,
	       const Chem2D_pState &dWdy,
	       const Vector2D &X, 
	       const int Axisymmetric, 
	       const double d_dWdx_dW,
	       const double d_dWdy_dW) const;
  //@}

  //@{ @name Source term and Jacobian associated with finite-rate chemistry.
  Chem2D_cState Sw(int &REACT_SET_FLAG, const int &Flow_Type) const;
  void dSwdU(DenseMatrix &dSwdU) const;
  double dSwdU_max_diagonal(const int &Preconditioned,
			    const double &delta_n) const;

  //@{ @name Source term and Jacobian associated with gravitational forces.
  Chem2D_cState Sg(void) const;
  void dSgdU(DenseMatrix &dSgdU) const;
  //@}

  //@{ @name Source terms associated with turbulence model.
  Chem2D_cState S_turbulence_model(const Chem2D_pState &dWdx, 
                                   const Chem2D_pState &dWdy, 
                                   const Vector2D &X,
                                   const int Axisymmetric);
  //@}

  /**************** Operators Overloading ********************/

  //@{ @name Index operator.
  double &operator[](int index);
  const double &operator[](int index) const;
  //@}

  //@{ @name Binary arithmetic operators.
  Chem2D_pState operator +(const Chem2D_pState &W) const;
  Chem2D_pState operator -(const Chem2D_pState &W) const;
  Chem2D_pState operator *(const double &a) const;
  friend Chem2D_pState operator *(const double &a, const Chem2D_pState &W);
  Chem2D_pState operator /(const double &a) const;
  double operator *(const Chem2D_pState &W) const;
  Chem2D_pState operator ^(const Chem2D_pState &W) const;
  //@}

  //@{ @name Assignment Operator.
  Chem2D_pState& operator =(const Chem2D_pState &W); 
  //@}

  //@{ @name Shortcut arithmetic operators.
  Chem2D_pState& operator +=(const Chem2D_pState &W);
  Chem2D_pState& operator -=(const Chem2D_pState &W);
  Chem2D_pState& operator *=(const double &a);
  Chem2D_pState& operator /=(const double &a);
  //@}

  //@{ @name Unary arithmetic operators.
  //Chem2D_pState operator +(const Chem2D_pState &W);
  friend Chem2D_pState operator -(const Chem2D_pState &W);
  //@}

  //@{ @name Relational operators.
  friend int operator ==(const Chem2D_pState &W1, const Chem2D_pState &W2);
  friend int operator !=(const Chem2D_pState &W1, const Chem2D_pState &W2);
  //@}

  //@{ @name Input-output operators.
  friend ostream& operator << (ostream &out_file, const Chem2D_pState &W);
  friend istream& operator >> (istream &in_file,  Chem2D_pState &W);
  //@}

};


/*!
 * Class: Chem2D_cState
 *
 * @brief Conserved variable solution state class definition for an
 *        inviscid, laminar, or turbulent chemically reacting gas-flow.
 *
 * Conserved variable solution state class definition for an inviscid, 
 * laminar, or turbulent chemically reacting gas-flow.
 */
class Chem2D_cState {
private: 
protected:
public:
  //@{ @name Cerved variables and associated constants:
  double                      rho; //!< Density.
  Vector2D                   rhov; //!< Momentum (2D)
  double                        E; //!< Total Energy (rho *(e + HALF*v^2))
  Species                *rhospec; //!< Species class using (rho*c[ns])
  double                     rhok; //!< Total turbulent kinetic energy.
  double                 rhoomega; //!< Total turbulent specific dissipation rate.
  Tensor2D                    tau; //!< Shear Stress Tensor
  Vector2D                  qflux; //!< Heat Flux Vector  
  Tensor2D                 lambda; //!< Reynolds Stress Tensor
  Vector2D                  theta; //!< Turbulent Heat Flux Vector  
  static int                   ns; //!< number of species
  static NASARP1311data *specdata; //!< Global species data 
  static double          *Schmidt; //!< Schmidt Number for each species
  static double    low_temp_range; //!< Low temp data range
  static double   high_temp_range; //!< High temp data range
  static int       NUM_VAR_CHEM2D; //!< Number of Chem2d variables (4+ns)
  static int            flow_type; //!< Flow-type indicator (inviscid, laminar, or k-omega turbulent).
  static double              Mref; //!< Mref for Precondtioning (normally set to incoming freestream Mach)
  static int          debug_level; //!< Debug level flag (0=none, 1,2,..) 
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
  static double alpha;
  static double sigma;
  static double sigma_star;
  static double beta;
  static double f_beta;
  static double beta_star;
  static double f_beta_star;
  static double Coeff_edm;
  //@}

  //@{ @name Creation, copy, and assignment constructors.
  Chem2D_cState(){rho = DENSITY_STDATM; rhov.zero(); E = PRESSURE_STDATM/(rho*(0.4)); 
                  rhok = ZERO; rhoomega = ZERO; tau.zero(); qflux.zero(); lambda.zero(); 
                  theta.zero(); rhospec = NULL; set_initial_values(); }   

  Chem2D_cState(const double &d, const double &vx, const double &vy, const double &En,
		const double &dk, const double &domega )
                {rho = d; rhov.x=vx; rhov.y=vy; E = En; rhok = dk; rhoomega = domega;  
		tau.zero(); qflux.zero(); lambda.zero(); theta.zero(); rhospec = NULL;  set_initial_values(); }

  Chem2D_cState(const double &value){ rho = value; rhov.x=value; rhov.y=value;
                                      E = value; rhok = ZERO; rhoomega = ZERO; tau.zero(); qflux.zero();
                                     lambda.zero(); theta.zero();  rhospec = NULL; set_initial_values(value); }   

  Chem2D_cState(const double &d, const double &vx, const double &vy, const double &En)
                {rho = d; rhov.x=vx; rhov.y=vy; E = En; rhok = ZERO; rhoomega = ZERO; tau.zero(); qflux.zero();
		lambda.zero(); theta.zero(); rhospec = NULL;  set_initial_values(); }

  Chem2D_cState(const double &d, const Vector2D &V, const double &En, const double &dk,const double &domega)
                {rho = d; rhov = V; E = En; rhok = dk; rhoomega = domega; tau.zero(); qflux.zero();
		lambda.zero(); theta.zero(); rhospec = NULL;  set_initial_values(); }

  Chem2D_cState(const double &d, const double &vx, const double &vy, const double &En,	const double &dk,
		const double &domega, Species *rhomfrac) 
                {rho = d; rhov.x=vx; rhov.y=vy; E = En; rhok = dk; rhoomega = domega;  tau.zero(); qflux.zero();
		lambda.zero(); theta.zero(); rhospec = NULL;  set_initial_values(rhomfrac); }

  Chem2D_cState(const double &d, const Vector2D &V, const double &En)
                {rho = d; rhov = V; E = En; rhok = ZERO; rhoomega = ZERO; tau.zero(); qflux.zero();
		 lambda.zero(); theta.zero(); rhospec = NULL;  set_initial_values(); } 

  Chem2D_cState(const double &d, const double &vx, const double &vy, const double &En, const double &value) 
                {rho = d; rhov.x=vx; rhov.y=vy; E = En; rhok = ZERO; rhoomega = ZERO; tau.zero(); qflux.zero();
                 lambda.zero(); theta.zero(); rhospec = NULL;  set_initial_values(value); }

  Chem2D_cState(const double &d, const double &vx, const double &vy, const double &En, 
                const double &kk, const double &oomega,	const double &value) 
                {rho = d; rhov.x=vx; rhov.y=vy; E = En; rhok = kk; rhoomega = oomega; tau.zero(); qflux.zero();
                 lambda.zero(); theta.zero(); rhospec = NULL;  set_initial_values(value); }

  Chem2D_cState(const double &d, const double &vx, const double &vy, const double &En, Species *rhomfrac) 
                {rho = d; rhov.x=vx; rhov.y=vy; E = En; rhok = ZERO; rhoomega = ZERO; tau.zero(); qflux.zero();
		lambda.zero(); theta.zero(); rhospec = NULL;  set_initial_values(rhomfrac); }

  Chem2D_cState(const double &d, const Vector2D &V, const double &En,  const double &dk, 
		const double &domega, Species *rhomfrac) 
                {rho = d; rhov = V; E = En; rhok = dk; rhoomega = domega;  tau.zero(); qflux.zero();
		lambda.zero(); theta.zero();
		rhospec = NULL;  set_initial_values(rhomfrac); }

  //this is needed for the operator overload returns!!!!
  Chem2D_cState(const Chem2D_cState &U)
               { rhospec = NULL; rho = DENSITY_STDATM; set_initial_values();  Copy(U); }

  //! Default destructor.
  ~Chem2D_cState(){ Deallocate(); }

  //! Deallocation for static specdata.
  void Deallocate_static(void){ if(specdata != NULL) delete[] specdata; 
                                  specdata = NULL; 
				if(Schmidt != NULL) delete[] Schmidt; 
				  Schmidt = NULL; 
                              }
  void Deallocate(void){ if(rhospec != NULL) delete[] rhospec; rhospec = NULL;  }
  //@}

  //@{ @name Set static variables.
  //! Read in ns species data, call only once as its static
  void set_species_data(const int &,const string *,const char *,
			const int &, const double&,const double *);
  
  //! Set initial data values predominately used internally 
  void set_initial_values();
  void set_initial_values(const double &value);
  void set_initial_values(double *rhomfrac);
  void set_initial_values(Species *rhomfrac);

  //! Set flow type static variable.
  void set_flow_type(const int &Flow_Type) { flow_type = Flow_Type; }

  //! Set turbulence static variables.
  void set_turbulence_variables(const double &C_constant,
				const double &von_karman_constant,
				const double &yplus_sub,
				const double &yplus_buffer,
				const double &yplus_outer);
  //@}

  //@{ @name Useful operators.
  //! Copy operator.
  void Copy(const Chem2D_cState &U);

  //! Vacuum operator.
  void Vacuum(){
    rho=ZERO; rhov.x=ZERO; rhov.y=ZERO; E=ZERO; rhok = ZERO; rhoomega = ZERO; 
    for(int i=0; i<ns; i++){
      rhospec[i].Vacuum();
    }
    tau.zero();  qflux.zero(); lambda.zero(); theta.zero(); 
  }  

  //! Zero non-solution quantities.
  void zero_non_sol(){
    for(int i=0; i<ns; i++){
      rhospec[i].gradc.zero();
      rhospec[i].diffusion_coef=ZERO;
    }
    tau.zero();  qflux.zero();
    lambda.zero(); theta.zero(); 
 
  }  

  //! Set the lowest value of the temperature range.
  void Temp_low_range();

  //! Set the highest value of the temperature range.
  void Temp_high_range();

  //! Check for unphysical state properties.
  int Unphysical_Properties(void) const;
  //@}
 
  //@{ @name Functions required for multigrid.
  //! Copy variables solved by multigrid only.
  void Copy_Multigrid_State_Variables(const Chem2D_cState &Ufine);

  //! Zero variables not-solved by multigrid.
  void Zero_Non_Multigrid_State_Variables(void);
  //@}

  //@{ @name Mixing Rules.

  // The following constructors return "total" physical
  // parameters based on mixture rules for each.

  // Mixture molecular mass
  //double Mass(void);
  double Rtot(void) const; 
  // Mixture heat capacity (Pressure constant)
  //double Cp(void) const;
  // Mixture heat capacity (Volume constant)
  //double Cv(void) const;
  // Mixture heat capacity ratio
  //double g(void) const;
  double gamma_guess(void) const;   //mixture specifc heat ratio
  //! Mixture specific internal energy
  double e(void) const;
  //! Mixture specific sensible internal energy.
  double es(void) const;
  //! Mixture specific enthalpy
  double h(const double &T) const;
  //! Mixture specific enthalpy
  double hs(const double &T) const;
  double hprime(const double &T) const; 
  //! Mixture heat of formation.
  double heatofform(void) const;
  // Mixture total internal energy
  //double E(void) const;
  // Mixture total enthalpy
  //double H(void);
  //! Mixture viscosity.
  double mu(void) const;
  //! Mixture thermal conductivity.
  double kappa(void) const;
  //@}

  //@{ @name State functions.
  //! Velocity.
  Vector2D v(void) const;
  //! Pressure.
  double p(void) const;
  //! Turbulent kinetic energy.
  double k(void) const;
  //! Turbulent specific dissipation.
  double omega(void) const;
  //! Temperature
  double T(void) const;
  //! Speed of sound.
  double a(void) const;
  //! Check for -ve mass fractions and sets small -ve c's to ZERO.
  bool negative_speccheck(const int &step) const;
  //! Check for unphysical state properties.
  bool Unphysical_Properties_Check(const int n) const;   
  double sum_species(void) const;
  //@}

  //@{ @name Temperature Derivatives.
  double dmudT(void) const;
  //@}

  //@{ @name Heat Flux vector thermal Diffusion.
  Vector2D thermal_diffusion(const double &Temp) const;
  //@}

  //@{ @name Turbulence related functions.
  double eddy_viscosity(void) const;      
  double Pr_turb(void) const;      
  double Sc_turb(void) const;      
  double Kappa_turb(void) const;     
  double Dm_turb(void) const;  //!< Molecular diffusivity.
  double omega_sublayer_BC(const double &y) const;
  //@}

  //@{ @name Primitive solution state.
  Chem2D_pState W(void) const;
  Chem2D_pState W(const Chem2D_cState &U) const;
  friend Chem2D_pState W(const Chem2D_cState &U);
  //@}

  //@{ @name Viscous Flux (laminar+turbulent)
  Chem2D_cState Viscous_Flux_x(const Chem2D_pState &dWdx) const;
  Chem2D_cState Viscous_Flux_y(const Chem2D_pState &dWdy) const;
  //@}

  //@{ @name Low Mach number preconditioner.
  void Low_Mach_Number_Preconditioner(DenseMatrix &P,
				      const double &deltax) const; 
  void Low_Mach_Number_Preconditioner_Inverse(DenseMatrix &Pinv,
					      const double &deltax) const; 
  //@}

  //@{ @name Index operators.
  double &operator[](int index);
  const double &operator[](int index) const;
  //@}

  //@{ @name Binary arithmetic operators.
  Chem2D_cState operator +(const Chem2D_cState &U) const;
  Chem2D_cState operator -(const Chem2D_cState &U) const;
  Chem2D_cState operator *(const double &a) const;
  friend Chem2D_cState operator *(const double &a, const Chem2D_cState &U);
  Chem2D_cState operator /(const double &a) const;
  double operator *(const Chem2D_cState &U) const;
  Chem2D_cState operator ^(const Chem2D_cState &U) const;
  //@}

  //@{ @name Assignment Operator.
  Chem2D_cState& operator =(const Chem2D_cState &U); 
  //@}

  //@{ @name Shortcut arithmetic operators.
  Chem2D_cState& operator +=(const Chem2D_cState &U);
  Chem2D_cState& operator -=(const Chem2D_cState &U);
  Chem2D_cState& operator *=(const double &a);
  Chem2D_cState& operator /=(const double &a);
  //@}

  //@{ @name Unary arithmetic operators.
  //Chem2D_cState operator +(const Chem2D_cState &U);
  friend Chem2D_cState operator -(const Chem2D_cState &U);
  //@}

  //@{ @name Relational operators.
  friend int operator ==(const Chem2D_cState &U1, const Chem2D_cState &U2);
  friend int operator !=(const Chem2D_cState &U1, const Chem2D_cState &U2);
  //@}

  //@{ @name Input-output operators.
  friend ostream& operator << (ostream &out_file, const Chem2D_cState &U);
  friend istream& operator >> (istream &in_file,  Chem2D_cState &U);
  //@}

};

/**************************************************************************
********************* CHEM2D_PSTATE CONSTRUCTORS **************************
***************************************************************************/

/**************************************************************************
   Set up mass fractions memory and initial values, ie. Species Class  
**************************************************************************/
inline void Chem2D_pState::set_initial_values(){
  Deallocate();
  spec = new Species[ns];
  for(int i=0; i<ns; i++){
    spec[i].c = ONE/ns ; 
  }
}

inline void  Chem2D_pState::set_initial_values(const double &value){
  Deallocate();
  spec = new Species[ns];
  for(int i=0; i<ns; i++){
    spec[i].c = value ; 
  }
}

//user specified
inline void Chem2D_pState::set_initial_values(double *cfrac){
  Deallocate();
  spec = new Species[ns];
  for(int i=0; i<ns; i++){
    spec[i].c = cfrac[i];
  }
}

//another set using species class
inline void Chem2D_pState::set_initial_values(Species *mfrac){
  Deallocate();
  spec = new Species[ns];
  for(int i=0; i<ns; i++){
    spec[i] = mfrac[i];
  }
}

/**********************************************************************
 * Chem2D_pState::set_turbulence_variables -- Set the turbulence      *
 *                                            static variables.       *
 **********************************************************************/
inline void Chem2D_pState::set_turbulence_variables(const double &C_constant,
						    const double &von_karman_constant,
						    const double &yplus_sub,
						    const double &yplus_buffer,
						    const double &yplus_outer) {
  C = C_constant;
  von_karman = von_karman_constant;
  yplus_sublayer = yplus_sub;
  yplus_buffer_layer = yplus_buffer;
  yplus_outer_layer = yplus_outer;
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

/*****************  Momentum *******************************/
inline Vector2D Chem2D_pState::rhov(void) const{
  return rho*v; 
} 

/********************** Prandtl ****************************/
inline double Chem2D_pState::Prandtl(void) const{
  //Pr = Cp*mu/k
  return Cp()*mu()/kappa();
}

/********************** Schmidt ****************************/
inline double Chem2D_pState::Schmidt_No(const int &i) const{
  if(spec[i].diffusion_coef != ZERO){
    return mu()/(rho*spec[i].diffusion_coef);
  } else {
    return Schmidt[i];
  }
  
}

/********************** Lewis *****************************/
inline double Chem2D_pState::Lewis(const int &i) const{
  if(spec[i].diffusion_coef != ZERO){
    return kappa()/(rho*Cp()*spec[i].diffusion_coef);
  }
  return ZERO;
}

/******* Mixture Diffusion Coefficient ********************/
// inline double Chem2D_pState::Diffusion_coef(void) const{
//   double sum=ZERO;
//   for(int i=0; i<ns; i++){
//     sum += spec[i].c * spec[i].diffusion_coef;
//   }
//   return sum;
// }

/************* Temperature ********************************/
inline double Chem2D_pState::T(void) const{
  return p/(rho*Rtot());
}
//strain rate tensor 
inline Tensor2D Chem2D_pState::Strain_Rate(const Chem2D_pState &dWdx,
					   const Chem2D_pState &dWdy,
					   const int Axisymmetric,
					   const Vector2D X){

  Tensor2D strain_rate;

  double r, div_v;
 
  /***************** Strain rate (+dilatation) **********************/	
  div_v = dWdx.v.x + dWdy.v.y;
  if (Axisymmetric == 2) {
    r = X.x; 
    div_v += v.x/r;
  } else if (Axisymmetric == 1) {
    r = X.y;
    div_v += v.y/r;
  } /* endif */

  strain_rate.xx = dWdx.v.x-div_v/THREE;
  strain_rate.xy = HALF*(dWdx.v.y + dWdy.v.x);
  strain_rate.yy = dWdy.v.y-div_v/THREE;
  
  if (Axisymmetric == 0) {
    strain_rate.zz = -(strain_rate.xx + strain_rate.yy); 
  } else if (Axisymmetric == 2) {
    strain_rate.zz = v.x/r-div_v/THREE;
  } else if (Axisymmetric == 1) {
    strain_rate.zz = v.y/r-div_v/THREE;
  } /* endif */
  
  return strain_rate;
  
}
//laminar (molecular) fluid stress
inline Tensor2D Chem2D_pState::Laminar_Stress(const Chem2D_pState &dWdx,
					      const Chem2D_pState &dWdy,
					      const int Axisymmetric,
					      const Vector2D X){
 /*  Tensor2D laminar_stress; */
 
/*   laminar_stress = mu()* Strain_Rate(dWdx,dWdy,Axisymmetric,X); */
  
/*   return laminar_stress; */

  tau = mu()* Strain_Rate(dWdx,dWdy,Axisymmetric,X);
  
  return tau;
  
}
// Reynolds stress
inline Tensor2D Chem2D_pState::Reynolds_Stress(const Chem2D_pState &dWdx,
					       const Chem2D_pState &dWdy,
					       const int Axisymmetric,
					       const Vector2D X){
/*   Tensor2D reynolds_stress; */
/*   double mu_t; */
/*   mu_t = eddy_viscosity(); */
  
/*   reynolds_stress = mu_t*Strain_Rate(dWdx,dWdy,Axisymmetric, X); */
  
/*   return reynolds_stress; */


 
  double mu_t;
  mu_t = eddy_viscosity();
  
  lambda = mu_t*Strain_Rate(dWdx,dWdy,Axisymmetric, X);
  
  return lambda;



}
//The passing parameter const Chem2D_pState &W is necessary when using the 
//solution parameter on the face (for computing viscous fluxes) 
inline Tensor2D Chem2D_pState::Strain_Rate(const Chem2D_pState &W,
					   const Chem2D_pState &dWdx,
					   const Chem2D_pState &dWdy,
					   const int Axisymmetric,
					   const Vector2D X){

  Tensor2D strain_rate;

  double  r, div_v;
  
   /***************** Strain rate (+dilatation) **********************/	
  div_v = dWdx.v.x + dWdy.v.y;
  if (Axisymmetric == 2) {
    r = X.x; 
    div_v += W.v.x/r;
  } else if (Axisymmetric == 1) {
    r = X.y;
    div_v += W.v.y/r;
  } /* endif */

  strain_rate.xx = dWdx.v.x-div_v/THREE;
  strain_rate.xy = HALF*(dWdx.v.y + dWdy.v.x);
  strain_rate.yy = dWdy.v.y-div_v/THREE;
  
  if (Axisymmetric == 0) {
    strain_rate.zz = -(strain_rate.xx + strain_rate.yy); 
  } else if (Axisymmetric == 2) {
    strain_rate.zz = W.v.x/r-div_v/THREE;
  } else if (Axisymmetric == 1) {
    strain_rate.zz = W.v.y/r-div_v/THREE;
  } /* endif */
  
  return strain_rate;
  
}
//laminar (molecular) fluid stress
inline Tensor2D Chem2D_pState::Laminar_Stress(const Chem2D_pState &W,
					      const Chem2D_pState &dWdx,
					      const Chem2D_pState &dWdy,
					      const int Axisymmetric,
					      const Vector2D X){
  Tensor2D laminar_stress;
  double mu;
  mu = W.mu();
 
  laminar_stress = mu* Strain_Rate(W,dWdx,dWdy,Axisymmetric,X);
  
  return laminar_stress;
  
}
// Reynolds stress
inline Tensor2D Chem2D_pState::Reynolds_Stress(const Chem2D_pState &W,
					       const Chem2D_pState &dWdx,
					       const Chem2D_pState &dWdy,
					       const int Axisymmetric,
					       const Vector2D X){
  Tensor2D reynolds_stress;
  double mu_t;
  mu_t = W.eddy_viscosity();

  reynolds_stress = mu_t*Strain_Rate(W,dWdx,dWdy,Axisymmetric, X);

  return reynolds_stress;
}
//Check for unphysical properties
/**********************************************************/
/* If unphysical properties and using global timestepping */ 
/* stop simulation                                        */
/**********************************************************/ 

inline bool Chem2D_pState::Unphysical_Properties_Check(Chem2D_cState &U, const int n) const{
  if ((U.flow_type == FLOWTYPE_INVISCID ||
       U.flow_type == FLOWTYPE_LAMINAR) &&
      (U.rho <= ZERO ||!U.negative_speccheck(n) ||U.es() <= ZERO)) {
    cout << "\n " << CFDkit_Name() 
	 << " Chem2D ERROR: Negative Density || Energy || mass fractions: \n"
	 << U << "\n";
    return false;
  }
  if ((U.flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON ||
       U.flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) &&
      (U.rho <= ZERO || !U.negative_speccheck(n) ||U.es() <= ZERO ||
       U.rhok < ZERO || U.rhoomega < ZERO)) {
    cout << "\n " << CFDkit_Name() 
	 << " Chem2D ERROR: Negative Density || Energy || mass fractions || Turbulent kinetic energy || : \n"
	 << " Dissipation rate per unit turbulent kinetic energy || : \n"
	 <<U << "\n";
    return false;
  }else{
    return true;
  }  
} 
/**************************************************************
  Check for -ve mass fractions and set small -ve values
  to ZERO. Then check that no mass is lost and that 
  the mass fractions still sum = 1

  Return "true" if passes
  and "false" if failed

  Should add "strict" and "anything goes" flags
  currently only set to give warnings, but still 
  continues.

***************************************************************/
inline bool Chem2D_pState::negative_speccheck(void) const{
  double LOCAL_TOL = 0.0001; 
  //-------- Negative Check ------------//
  for(int i=0; i<ns; i++){
    if(spec[i].c < ZERO){  //check for -ve
      if(spec[i].c > -LOCAL_TOL){  //check for small -ve and set to ZERO 
	spec[i].c = ZERO;
      } else {
	spec[i].c = ZERO;
	if(debug_level){ 	
	  cout <<"\n -ve mass fraction in "<<specdata[i].Speciesname()<<" "<<
	       spec[i].c<<" greater than allowed tolerance of "<<-LOCAL_TOL; 
	}
      }
    } 
  } 
  return(1);
}

/***** Species Concentrations ******************************/
inline double Chem2D_pState::SpecCon(int i) const{
  //returned in kg/m^3 / kg/mol => mol/m^3
  return (rho)*spec[i].c/(specdata[i].Mol_mass());
}

/******* GIBBS Free Energy ********************************
  Eqn. (10.84) (10.87) Anderson
  Gs = Hs - TS
  Gs(ps=1) = Gs - R_UNIVERSAL*T*ln(ps)  //for data not at 1atm
  ps = cs(M/Ms)*p
***********************************************************/
inline double Chem2D_pState::Gibbs(int species) const{
  double Temp = T(); 
//   double ps = spec[species].c*(Mass()/specdata[species].Mol_mass())*p;
//   if(ps != ZERO){  //check ps as log(0) -> INF
//     return specdata[species].Enthalpy_mol(Temp)  
//       - Temp*specdata[species].Entropy_mol(Temp) 
//       - R_UNIVERSAL*Temp*log(ps);
//   } 
//   else{
    return specdata[species].Enthalpy_mol(Temp)  
      - Temp*specdata[species].Entropy_mol(Temp);
    //  }
}

/************** Temperature Derivatives *******************/
inline double Chem2D_pState::diedip() const{
  double Rlocal = Rtot();
  return (hprime() - Rlocal)/(rho*Rlocal);
}

inline double Chem2D_pState::diedirho() const{
  double Rlocal = Rtot();
  return -p*(hprime() - Rlocal)/(rho*rho*Rlocal);
}

/**************** Copy *************************************/
inline void Chem2D_pState::Copy(const Chem2D_pState &W){
  rho = W.rho;
  v = W.v; 
  p = W.p;  
  k = W.k;
  omega = W.omega;
  for( int i=0; i<ns; i++){
    spec[i] = W.spec[i];
  }
  tau = W.tau;
  qflux = W.qflux;
  lambda = W.lambda;
  theta = W.theta;
}
//----------------- Index Operator ------------------------/
//should probably inline these guys for speed
inline double& Chem2D_pState::operator[](int index) {
  
    assert( index >= 1 && index <= NUM_VAR_CHEM2D );
    if(index == 1){
      return (rho);
    } else if(index ==  2) {
      return (v.x);
    } else if(index ==  3) {
      return (v.y);
    } else if(index ==  4) {
      return (p);
    } else if(index ==   5) {
      return (k);
    } else if(index ==   6) {
      return (omega);
    } else {
      return spec[index-7].c;
    }

}

inline const double& Chem2D_pState::operator[](int index) const {
  
    assert( index >= 1 && index <= NUM_VAR_CHEM2D );
    if(index == 1){
      return (rho);
    } else if(index ==  2) {
      return (v.x);
    } else if(index ==  3) {
      return (v.y);
    } else if(index ==  4) {
      return (p);
    }else if(index ==   5 ) {
      return (k);
    }else if(index ==   6) {
      return (omega);
    } else{
      return spec[index-7].c;
    }
 
}
/**************************************************************
  Get max of the min temperature of the lowest region
  and min of the max temperature of the highest region
***************************************************************/
inline void Chem2D_pState::Temp_low_range(void){  
  double temp = specdata[0].Low_range();
  for(int i=0; i<ns; i++){
    temp = max(specdata[i].Low_range(),temp);
  }
  low_temp_range = temp;  
}

inline void Chem2D_pState::Temp_high_range(void){
  double temp = specdata[0].High_range();
  for(int i=0; i<ns; i++){
    temp = min(specdata[i].High_range(),temp);
  }
  high_temp_range = temp;  
}



/********************************************************
 * Chem2D_pState::U -- Conserved solution state.        *
 ********************************************************/
inline Chem2D_cState Chem2D_pState::U(void) const {
  return U(*this);
}

inline Chem2D_cState Chem2D_pState::U(const Chem2D_pState &W) const{
  if(ns == W.ns){ //check that species are equal   
    Chem2D_cState Temp;
    Temp.rho = W.rho;
    Temp.rhov = W.rhov();
    Temp.E = W.E();
    for(int i=0; i<W.ns; i++){
      Temp.rhospec[i] = W.rho*W.spec[i];
      Temp.rhospec[i].gradc = W.rho*W.spec[i].gradc;
      Temp.rhospec[i].diffusion_coef = W.rho*W.spec[i].diffusion_coef;
    } 
    Temp.rhok = W.rho*W.k;
    Temp.rhoomega = W.rho*W.omega;
    Temp.tau = W.tau;
    Temp.qflux = W.qflux; 
    Temp.lambda = W.lambda;
    Temp.theta = W.theta; 
    return Temp;
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  } 
}

inline Chem2D_cState U(const Chem2D_pState &W) {
  Chem2D_cState Temp;
  Temp.rho = W.rho;
  Temp.rhov = W.rhov();
  Temp.E = W.E();
  for(int i=0; i<W.ns; i++){
    Temp.rhospec[i] = W.rho*W.spec[i];
    Temp.rhospec[i].gradc = W.rho*W.spec[i].gradc;
    Temp.rhospec[i].diffusion_coef = W.rho*W.spec[i].diffusion_coef;
  }  
  Temp.rhok = W.rho*W.k;
  Temp.rhoomega = W.rho*W.omega;
  Temp.tau = W.tau;
  Temp.qflux = W.qflux; 
  Temp.lambda = W.lambda;
  Temp.theta = W.theta; 
  return Temp;
}

inline void dFIdU(DenseMatrix &dFdU, const Chem2D_pState &W) {

  double Rsum = W.Rtot();
  double Temp = W.T();
  double term1 = W.Cp() - Rsum;
  double term2 = W.Cp();
  double g = W.g();
  double gm1 = g-ONE;
  double gm1i = ONE/gm1;

  dFdU(0,1) += ONE;
  dFdU(1,0) += HALF*(sqr(W.v.x)*(g-THREE)+sqr(W.v.y)*(g-1));
  dFdU(1,1) -= W.v.x*(g-THREE);
  dFdU(1,2) -= W.v.y*gm1;
  dFdU(1,3) += gm1;
  dFdU(2,0) -= W.v.x*W.v.y;
  dFdU(2,1) += W.v.y;
  dFdU(2,2) += W.v.x;

  double phi = ZERO;   
  for(int i=0; i<W.ns-1; i++) {
    phi += W.spec[i].c*(W.specdata[i].Enthalpy(Temp) + W.specdata[i].Heatofform());   
  }
  dFdU(3,0) += W.v.x*((HALF*g -1)*(W.v.x*W.v.x +W.v.y+W.v.y) - phi);
  dFdU(3,1) += (ONE -g)*W.v.x + W.H()/W.rho;
  dFdU(3,2) -= W.v.x*W.v.y*gm1;
  dFdU(3,3) += g*W.v.x; 
 
  for(int i = 0; i<(W.ns-1); i++) {
    dFdU(1,NUM_CHEM2D_VAR_SANS_SPECIES+i) -=(g-ONE)*(W.specdata[i].Enthalpy(Temp)+W.specdata[i].Heatofform())- g*W.specdata[i].Rs()*Temp;
    dFdU(3,NUM_CHEM2D_VAR_SANS_SPECIES+i) = W.v.x*dFdU(1,NUM_CHEM2D_VAR_SANS_SPECIES+i);
    dFdU(NUM_CHEM2D_VAR_SANS_SPECIES+i, 0) -= W.spec[i].c*W.v.x ;
    dFdU(NUM_CHEM2D_VAR_SANS_SPECIES+i, 1) += W.spec[i].c ;
    dFdU(NUM_CHEM2D_VAR_SANS_SPECIES+i,NUM_CHEM2D_VAR_SANS_SPECIES+i) += W.v.x ;    
  }

  if (W.flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
      W.flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON){ 
    dFdU(3,0) -= W.v.x*W.k;
    dFdU(3,4) += (ONE-g)*W.v.x; 
    dFdU(4,0) = -W.k*W.v.x;
    dFdU(4,1) = W.k;
    dFdU(4,4) = W.v.x;
    dFdU(5,0) = -W.omega*W.v.x;
    dFdU(5,1) = W.omega;
    dFdU(5,5) = W.v.x;
  }

}

inline void Chem2D_pState::dWdU(DenseMatrix &dWdQ) {
  

  dWdQ(0,0) = ONE;
  dWdQ(1,0) = -v.x/rho;
  dWdQ(1,1) = ONE/rho;
  dWdQ(2,0) = -v.y/rho;
  dWdQ(2,2) = ONE/rho;

  double Temp = T();
  double Rt = Rtot();
  double C_p = Cp();
  double C_v1 = Rt-C_p;

  
  dWdQ(3,0) = -ONE/TWO *sqr(v.x+v.y)/(rho*C_v1);
  dWdQ(3,1) = v.x/(rho*C_v1);
  dWdQ(3,2) = v.y/(rho*C_v1);
  dWdQ(3,3) = -ONE/(rho*C_v1);

  int NUM_VAR = NUM_CHEM2D_VAR_SANS_SPECIES;
  for(int i=0; i<(ns-1);i++){
    
    dWdQ(3, NUM_VAR+i) = (specdata[i].Enthalpy(Temp)+specdata[i].Heatofform()-specdata[i].Rs()*Temp)/(rho*C_v1);
    dWdQ(NUM_VAR+i, 0) = -spec[i].c/rho;
    dWdQ(NUM_VAR+i, NUM_VAR+i) = ONE/rho;
  }
  
  if (flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
      flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
    
    dWdQ(3,4) = ONE/(rho*C_v1);
    
    dWdQ(4,0) = -k/rho;
    dWdQ(5,0) = -omega/rho;
    dWdQ(4,4) =  ONE/rho;
    dWdQ(5,5)=  ONE/rho;
  
  }

}

/**************************************************************************
********************* CHEM2D_CSTATE CONSTRUCTORS **************************
***************************************************************************/

/**************************************************************
  Set up mass fractions memory and initial values, ie. Species Class  
****************************************************************/
inline void Chem2D_cState::set_initial_values(){
  Deallocate();
  rhospec = new Species[ns];
  for(int i=0; i<ns; i++){
    rhospec[i].c = rho/ns; 
  }
}

inline void  Chem2D_cState::set_initial_values(const double &value){
  Deallocate();
  rhospec = new Species[ns];
  for(int i=0; i<ns; i++){
    rhospec[i].c = value; 
  }
}

//user specified
inline void Chem2D_cState::set_initial_values(double *rhomfrac){
  Deallocate();
  rhospec = new Species[ns];
  for(int i=0; i<ns; i++){
    rhospec[i].c = rhomfrac[i];
  }
}

//another set using species class
inline void Chem2D_cState::set_initial_values(Species *rhomfrac){
  Deallocate();
  rhospec = new Species[ns];
  for(int i=0; i<ns; i++){
    rhospec[i] = rhomfrac[i];
  }
}

/**********************************************************************
 * Chem2D_cState::set_turbulence_variables -- Set the turbulence      *
 *                                            static variables.       *
 **********************************************************************/
inline void Chem2D_cState::set_turbulence_variables(const double &C_constant,
						    const double &von_karman_constant,
						    const double &yplus_sub,
						    const double &yplus_buffer,
						    const double &yplus_outer) {
  C = C_constant;
  von_karman = von_karman_constant;
  yplus_sublayer = yplus_sub;
  yplus_buffer_layer = yplus_buffer;
  yplus_outer_layer = yplus_outer;
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

/**************** Copy *************************************/
inline void Chem2D_cState::Copy(const Chem2D_cState &U){
  rho = U.rho;
  rhov = U.rhov; 
  E = U.E; 
  rhok = U.rhok;
  rhoomega = U.rhoomega;

  for( int i=0; i<ns; i++){ 
    rhospec[i] = U.rhospec[i];
  } 
  tau = U.tau;
  qflux = U.qflux; 
  lambda = U.lambda;
  theta = U.theta; 
}

/**********************************************************************
 * Chem2D_cState::Unphysical_Properties -- Check for unphysical state *
 *                                          properties.               *
 **********************************************************************/
inline int Chem2D_cState::Unphysical_Properties(void) const {
  if (rho <= ZERO || !negative_speccheck(ns) || es() <= ZERO || rhok < ZERO || rhoomega < ZERO) return 1;
  return 0;
}

/**********************************************************************
 * Chem2D_cState::Copy_Multigrid_State_Variables --                   *
 *                           Copy variables solved by multigrid only. *
 **********************************************************************/
inline void Chem2D_cState::Copy_Multigrid_State_Variables(const Chem2D_cState &Ufine) {
  Copy(Ufine);
  Zero_Non_Multigrid_State_Variables();
}

/**********************************************************************
 * Chem2D_cState::Zero_Non_Multigrid_State_Variables --               *
 *                            Zero variables not-solved by multigrid. *
 **********************************************************************/
inline void Chem2D_cState::Zero_Non_Multigrid_State_Variables(void) {
  rhok = ZERO; rhoomega = ZERO;
}

/**************** Velocity *********************************/
inline Vector2D Chem2D_cState::v() const{
  return (rhov/rho);
}

/**************** Pressure *********************************/
inline double Chem2D_cState::p() const{
  return (rho*Rtot()*T());
}
/**************** kinetic energy *********************************/
inline double Chem2D_cState::k() const{
  return (rhok/rho);
}

/**************** omega *********************************/
inline double Chem2D_cState::omega() const{
  return (rhoomega/rho);
}

//----------------- Index Operator ---------------------//
inline double& Chem2D_cState::operator[](int index) {
  assert( index >= 1 && index <= NUM_VAR_CHEM2D );

    if(index == 1){
      return (rho);
    } else if(index ==  2) {
      return (rhov.x);
    } else if(index ==  3) {
      return (rhov.y);
    } else if(index ==  4) {
      return (E);
    } else if(index ==  5 ) {
      return (rhok);
    }else if(index ==   6) {
      return (rhoomega);
    } else {
      for(int i=7; i<=(NUM_VAR_CHEM2D); i++){
	if(index ==i){
	  return rhospec[i-7].c;
	  break;
	}
      }
    }
 
}

inline const double& Chem2D_cState::operator[](int index) const{
  assert( index >= 1 && index <= NUM_VAR_CHEM2D );
 
    if(index == 1){
      return (rho);
    } else if(index ==  2) {
      return (rhov.x);
    } else if(index ==  3) {
      return (rhov.y);
    } else if(index ==  4) {
      return (E);
    } else if(index==   5){
      return (rhok);
    }else if(index==    6){
      return (rhoomega);
    }else{
      for(int i=7; i<=(NUM_VAR_CHEM2D); i++){
	if(index ==i){
	  return rhospec[i-7].c;
	  break;
	}
      }
    }
 
}

/**************************************************************
  Check for -ve mass fractions and set small -ve values
  to ZERO. Then check that no mass is lost and that 
  the mass fractions still sum = 1

  Return "true" if passes
  and "false" if failed

  Should add "strict" and "anything goes" flags
  currently only set to give warnings, but still 
  continues.

***************************************************************/
inline bool Chem2D_cState::negative_speccheck(const int &step) const{
  double sum = 0.0;
  double temp = 0.0;

  double LOCAL_TOL = 0.0001; 
  //-------- Negative Check ------------//
  for(int i=0; i<ns; i++){
    temp = rhospec[i].c/rho;
    if(temp < ZERO){  //check for -ve
      if(temp > -LOCAL_TOL){  //check for small -ve and set to ZERO 
	rhospec[i].c = ZERO;
      } else {
	if(debug_level){ 	
	  cout<<"\n -ve mass fraction in "<<specdata[i].Speciesname()<<" "<<
	    temp<<" greater than allowed tolerance of "<<-LOCAL_TOL; 
	}
	if( step < 10){
	  return false;
	} else {
	  rhospec[i].c = ZERO;
	  cout<<"\n Still -ve mass fractions, but continuing anyway ";
	  return true;
	}
      }
    } 
    sum += temp;
  } 
  
  //--------- Sum Check --------------//
  if( sum < ONE - LOCAL_TOL || sum > ONE + LOCAL_TOL){
    cout.precision(10);
    if(debug_level){ 
      cout<<"\n WARNING IN MASS FRACTION SUM != 1  ->"<<sum; 
    }
    if(step <10){
      return false;
    } else {
      cout<<"\n Still mass fractions sum !=1, but continuing anyway ";
      return true;
    }
  } else {
    return true;
  }
}

inline bool Chem2D_cState::Unphysical_Properties_Check(const int n) const{
  if ((flow_type == FLOWTYPE_INVISCID ||
       flow_type == FLOWTYPE_LAMINAR) &&
      (rho <= ZERO || !negative_speccheck(n) ||es() <= ZERO)) {
    cout << "\n " << CFDkit_Name() 
	 << " Chem2D ERROR: Negative Density || Energy || mass fractions: \n";
    return false;
  }
  if ((flow_type == FLOWTYPE_TURBULENT_RANS_K_EPSILON ||
       flow_type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) &&
      (rho <= ZERO || !negative_speccheck(n) ||es() <= ZERO ||
       rhok < ZERO ||rhoomega < ZERO)) {
    cout << "\n " << CFDkit_Name() 
	 << " Chem2D ERROR: Negative Density || Energy || mass fractions || Turbulent kinetic energy || : \n"
	 << " Dissipation rate per unit turbulent kinetic energy || : \n";
    return false;
  } else {
    return true ;
  }
  
} 

/**************************************************************
  Sum N-1 species.
***************************************************************/
inline double Chem2D_cState::sum_species(void) const{
  double sum=ZERO;
  for(int i=0; i<ns-1; i++){
    sum += rhospec[i].c;
  }
  return sum/rho;
}

/**************************************************************
  Get max and min temperature ranges for data
***************************************************************/
inline void Chem2D_cState::Temp_low_range(void){  
  double temp = specdata[0].Low_range();
  for(int i=0; i<ns; i++){
    temp = max(specdata[i].Low_range(),temp);
  } 
  low_temp_range = temp;  
}

inline void Chem2D_cState::Temp_high_range(void){
  double temp = specdata[0].High_range();
  for(int i=0; i<ns; i++){
    temp = min(specdata[i].High_range(),temp);
  } 
  high_temp_range = temp;  
}

/***************************************************************
 * Chem2D_cState::Sa -- Axisymmetric flow source terms.  *
 ***************************************************************/
// inline Chem2D_cState Chem2D_cState::Sa(const Vector2D &X) const{
//   Chem2D_cState Temp;
//   return Temp;
// }


/********************************************************
 * Chem2D_cState::W -- Primitive solution state.       *
 ********************************************************/
inline Chem2D_pState Chem2D_cState::W(void) const {
  return W(*this);
}

inline Chem2D_pState Chem2D_cState::W(const Chem2D_cState &U) const{
  if(ns == U.ns){ //check that species are equal   
    Chem2D_pState Temp;
    Temp.rho = U.rho;
    Temp.v = U.v();  
    Temp.p = U.p();
    Temp.k = U.k();
    Temp.omega = U.omega();
    for(int i=0; i<U.ns; i++){
      Temp.spec[i] = U.rhospec[i]/U.rho;
      Temp.spec[i].gradc = U.rhospec[i].gradc/U.rho;
      Temp.spec[i].diffusion_coef = U.rhospec[i].diffusion_coef/U.rho;
    }
    
    Temp.tau = U.tau;
    Temp.qflux = U.qflux; 
    Temp.lambda = U.lambda;
    Temp.theta = U.theta; 
   
    return Temp;
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  } 
}

inline Chem2D_pState W(const Chem2D_cState &U) {
  Chem2D_pState Temp;
  Temp.rho = U.rho;
  Temp.v = U.v();
  Temp.p = U.p();
  Temp.k = U.k();
  Temp.omega = U.omega();
  for(int i=0; i<U.ns; i++){
    Temp.spec[i] = U.rhospec[i]/U.rho;
    Temp.spec[i].gradc = U.rhospec[i].gradc/U.rho;
    Temp.spec[i].diffusion_coef = U.rhospec[i].diffusion_coef/U.rho;
  }
  Temp.tau = U.tau;
  Temp.qflux = U.qflux;
  Temp.lambda = U.lambda;
  Temp.theta = U.theta;
  
 
  return Temp;
}


/*********************************************************************************
 *********************************************************************************
 *********************************************************************************
 *********************** EXTERNAL FUNCTIONS **************************************
 *********************************************************************************
 *********************************************************************************
 *********************************************************************************/

/********************************************************
 * External Boundary Conditions Functions               *
 ********************************************************/
extern Chem2D_pState Reflect(const Chem2D_pState &W,
	      	              const Vector2D &norm_dir);

extern Chem2D_pState Free_Slip(const Chem2D_pState &Win,
			     const Chem2D_pState &Wout,
			     const Vector2D &norm_dir,
			     const int &TEMPERATURE_BC_FLAG);

extern Chem2D_pState No_Slip(const Chem2D_pState &Win,
			     const Chem2D_pState &Wout,
			     const Vector2D &norm_dir,
			     const int &TEMPERATURE_BC_FLAG);

extern Chem2D_pState Moving_Wall(const Chem2D_pState &Win,
				 const Chem2D_pState &Wout,
				 const Vector2D &norm_dir,				 
				 const double &wall_velocity,
				 const int &TEMPERATURE_BC_FLAG);

extern Chem2D_pState BC_Characteristic_Pressure(const Chem2D_pState &Wi,
						const Chem2D_pState &Wo,
						const Vector2D &norm_dir);

extern Chem2D_pState BC_Flame_Inflow(const Chem2D_pState &Wi,
				     const Chem2D_pState &Wo, 
				     const Chem2D_pState &Woutlet,
				     const Vector2D &norm_dir);

extern Chem2D_pState BC_Flame_Outflow(const Chem2D_pState &Wi,
				      const Chem2D_pState &Wo,
				      const Chem2D_pState &Winlet,
				      const Vector2D &norm_dir);

/*******************************************************
 * Exact Test Case Solution Functions                  *
 *******************************************************/
extern Chem2D_pState RinglebFlow(const Chem2D_pState &Wdum,
				  const Vector2D X);

extern Chem2D_pState ViscousChannelFlow(const Chem2D_pState &Wdum,
					const Vector2D X,
					const double Vwall,
					const double dp);

extern Chem2D_pState FlatPlate(const Chem2D_pState &Winf,
			       const Vector2D X,
			       double &eta,
			       double &f,
			       double &fp,
			       double &fpp);

/*******************************************************
 * External Flux Function Functions                    *
 *******************************************************/
extern Chem2D_pState RoeAverage(const Chem2D_pState &Wl,
	      	                 const Chem2D_pState &Wr);

// HLLE
extern Chem2D_cState FluxHLLE_x(const Chem2D_pState &Wl,
				const Chem2D_pState &Wr,
				const int &Preconditioning);

extern Chem2D_cState FluxHLLE_x(const Chem2D_cState &Ul,
				const Chem2D_cState &Ur,
				const int &Preconditioning);
  
extern Chem2D_cState FluxHLLE_n(const Chem2D_pState &Wl,
				const Chem2D_pState &Wr,
				const Vector2D &norm_dir,
				const int &Preconditioning);

extern Chem2D_cState FluxHLLE_n(const Chem2D_cState &Ul,
				const Chem2D_cState &Ur,
				const Vector2D &norm_dir,
				const int &Preconditioning);

// Linde
extern Chem2D_cState FluxLinde(const Chem2D_pState &Wl,
			       const Chem2D_pState &Wr);

extern Chem2D_cState FluxLinde(const Chem2D_cState &Ul,
			       const Chem2D_cState &Ur);


extern Chem2D_cState FluxLinde_n(const Chem2D_pState &Wl,
				 const Chem2D_pState &Wr,
				 const Vector2D &norm_dir);

extern Chem2D_cState FluxLinde_n(const Chem2D_cState &Ul,
				 const Chem2D_cState &Ur,
				 const Vector2D &norm_dir);

// Roe
extern Chem2D_pState WaveSpeedPos(const Chem2D_pState &lambda_a,
				  const Chem2D_pState &lambda_l,
				  const Chem2D_pState &lambda_r);
			

extern Chem2D_pState WaveSpeedNeg(const Chem2D_pState &lambda_a,
				  const Chem2D_pState &lambda_l,
				  const Chem2D_pState &lambda_r);

extern Chem2D_pState WaveSpeedAbs(const Chem2D_pState &lambda_a,
				  const Chem2D_pState &lambda_l,
				  const Chem2D_pState &lambda_r);
 
extern Chem2D_pState HartenFixPos(const Chem2D_pState &lambda_a,
				  const Chem2D_pState &lambda_l,
				  const Chem2D_pState &lambda_r );
				

extern Chem2D_pState HartenFixNeg(const Chem2D_pState &lambda_a,
				  const Chem2D_pState &lambda_l,
				  const Chem2D_pState &lambda_r);
				  

extern Chem2D_pState HartenFixAbs(const Chem2D_pState &lambdas_a,
				  const Chem2D_pState &lambdas_l,
				  const Chem2D_pState &lambdas_r );
				

extern Chem2D_cState FluxRoe_x(const Chem2D_pState &Wl,
			       const Chem2D_pState &Wr,
			       const int &Preconditioning, 
			       const double &deltax);

extern Chem2D_cState FluxRoe_x(const Chem2D_cState &Ul,
			       const Chem2D_cState &Ur,
			       const int &Preconditioning, 
			       const double &deltax);

extern Chem2D_cState FluxRoe_n(const Chem2D_pState &Wl,
			       const Chem2D_pState &Wr,
			       const Vector2D &norm_dir,
			       const int &Preconditioning, 
			       const double &deltax);

extern Chem2D_cState FluxRoe_n(const Chem2D_cState &Ul,
			       const Chem2D_cState &Ur,
			       const Vector2D &norm_dir,
			       const int &Preconditioning, 
			       const double &deltax);

/* Viscous Solution flux (laminar+turbulent) */
extern Chem2D_cState Viscous_FluxArithmetic_n(const Chem2D_cState &Ul,
   			                      const Chem2D_pState &dWdx_l,
			                      const Chem2D_pState &dWdy_l,
				              const Chem2D_cState &Ur,
   			                      const Chem2D_pState &dWdx_r,
			                      const Chem2D_pState &dWdy_r,
				              const Vector2D &norm_dir);

extern Chem2D_cState Viscous_Flux_n(const Chem2D_pState &W,
				    const Chem2D_pState &dWdx,
				    const Chem2D_pState &dWdy,
				    const int Axisymmetric,
				    const Vector2D X,			
				    const Vector2D &norm_dir);

extern double WallShearStress(const Chem2D_pState &W1,
			      const Vector2D &X1,
			      const Vector2D &X2,
			      const Vector2D &X3,
			      const Vector2D &norm_dir);

/*********************************************************/
extern Chem2D_pState Rotate(const Chem2D_pState &W,
	      	             const Vector2D &norm_dir);

extern DenseMatrix Rotation_Matrix2(Vector2D nface, int Size, int A_matrix); 

extern Vector2D HLLE_wavespeeds(const Chem2D_pState &Wl,
                                const Chem2D_pState &Wr,
                                const Vector2D &norm_dir); 

#endif //end _CHEM2D_STATE_INCLUDED 
