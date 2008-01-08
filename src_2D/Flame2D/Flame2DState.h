/////////////////////////////////////////////////////////////////////
///
/// \file Flame2DState.h
/// 
/// \author Marc R.J. Charest
/// 
/// \brief This header file contains the class definition for 
///        describing the state of the chemically reacting, 
///        multispecies ideal gas mixture.
///
/////////////////////////////////////////////////////////////////////
#ifndef _FLAME2D_STATE_INCLUDED 
#define _FLAME2D_STATE_INCLUDED

// Required C++ libraries
#include <iostream>
#include <cstdlib> 
using namespace std;

#include "../Math/Math.h"
#include "../Math/Matrix.h"
#include "../CFD/CFD.h"
#include "../Math/Tensor2D.h"
#include "../Math/Vector2D.h"
#include "../Physics/GasConstants.h"

// FLAME2D Specific headers
#include "Mixture.h"
#include "Soot2DState.h"

/////////////////////////////////////////////////////////////////////
/// Defines
/////////////////////////////////////////////////////////////////////

// Used in negative_speccheck for species round off (was MICRO)
#define SPEC_TOLERANCE  1e-8

// number of fixed variables in the Flame2D class
#define NUM_FLAME2D_VAR_SANS_SPECIES 4  //rho, v(2), p

// if static species && static soot var, the total number of variables
#ifdef STATIC_NUMBER_OF_SPECIES
#ifdef STATIC_NUM_SOOT2D_VAR
#define STATIC_NUM_FLAME2D_VAR NUM_FLAME2D_VAR_SANS_SPECIES + STATIC_NUMBER_OF_SPECIES + STATIC_NUM_SOOT2D_VAR
#endif
#endif

/////////////////////////////////////////////////////////////////////
/// Class Definitions
/////////////////////////////////////////////////////////////////////
// define early so can be used in other classes
class Flame2D_State;
class Flame2D_pState;

/**
 * Class: Flame2D_State
 *
 * @brief  Base solution state class definition for an
 *         inviscid, laminar, or turbulent chemically reacting gas-flow.
 *
 */
class Flame2D_State {

public:
  
  /******************* Constructors/Destructors *********************/
  Flame2D_State() { Nullify(); Allocate();  }
  Flame2D_State(const double &d, const double &vvx, const double &vvy,
		const double &pre, const double &val) { 
    Nullify(); Allocate(); 
    rho()=d; p()=pre; vx()=vvx; vy()=vvy; c(val); sc(0.0);
  }
  Flame2D_State(const double &d, const double &vvx, const double &vvy,
		const double &pre, const double* val) {
    Nullify(); Allocate(); 
    rho()=d; p()=pre; vx()=vvx; vy()=vvy; c(val); sc(0.0);
  }
  Flame2D_State(const Flame2D_State &X) { Nullify(); Allocate(); Copy(X);  }
  ~Flame2D_State() { Deallocate(); }
  void Copy( const Flame2D_State &X ) { for (int i=0; i<n; i++) x[i] = X.x[i]; }

  
  /*********************** Allocate / Deallocate ********************/
  void Allocate();
  void Deallocate();

  /*********************** Initializers *****************************/
  void Nullify(); 
  void Vacuum() { Zero(); }
  void Zero() { for (int i=0; i<n; i++) x[i] = 0.0; }

  /************************* Accessors ******************************/
  //@{ @name Primitive variables:
  //!< Density.
  double rho(void) const { return x[iRho]; }
  double& rho(void) { return x[iRho]; }

  //!< Flow velocity (2D)
  double vx(void) const { return x[iVx]; }
  double& vx(void) { return x[iVx]; }
  double vy(void) const { return x[iVy]; }
  double& vy(void) { return x[iVy]; }
  void v(const Vector2D &V) { x[iVx]=V.x; x[iVy]=V.y; };
  double vsqr(void) const { return x[iVx]*x[iVx] + x[iVy]*x[iVy]; };
  double vabs(void) const { return sqrt(x[iVx]*x[iVx] + x[iVy]*x[iVy]); };
  Vector2D v(void) const { return Vector2D(x[iVx], x[iVy]); };

  //!< Pressure.
  double p(void) const { return x[iPress]; }
  double& p(void) { return x[iPress]; }

  //!< Species mass fractions
  const double* c(void) const { return &x[iSpec]; }
  double c(const int&i) const { return x[iSpec+i]; }
  double* c(void) { return &x[iSpec]; }
  void c(const double *mfrac) { for (int i=0; i<ns; i++) x[iSpec+i]=mfrac[i]; }
  void c(const double &val)  { for (int i=0; i<ns; i++) x[iSpec+i]=val; }
  double& c(const int&i)  { return x[iSpec+i]; }

  //!< Soot scalars
  const double* sc(void) const { return &x[iSoot]; }
  double sc(const int&i) const { return x[iSoot+i]; }
  double* sc(void) { return &x[iSoot]; }
  void sc(const double *s) { for (int i=0; i<nsc; i++) x[iSoot+i]=s[i]; }
  void sc(const double &val)  { for (int i=0; i<nsc; i++) x[iSoot+i]=val; }
  double& sc(const int&i)  { return x[iSoot+i]; }
  //@}

  //@{ @name Conserved variables:
  //!< Momentum (2D)
  double rhovx(void) const { return x[iVx]; }
  double& rhovx(void) { return x[iVx]; }
  double rhovy(void) const { return x[iVy]; }
  double& rhovy(void) { return x[iVy]; }
  void rhov(const Vector2D &rhoV) { x[iVx]=rhoV.x; x[iVy]=rhoV.y; };
  double rhovsqr(void) const { return (x[iVx]*x[iVx] + x[iVy]*x[iVy]); };
  Vector2D rhov(void) const { return Vector2D(x[iVx], x[iVy]); };

  //!< Total Energy (rho *(e + HALF*v^2))
  double E(void) const { return x[iPress]; }
  double& E(void) { return x[iPress]; }

  //!< Species mass fractions times density
  const double* rhoc(void) const { return &x[iSpec]; }
  double rhoc(const int&i) const { return x[iSpec+i]; }
  double* rhoc(void) { return &x[iSpec]; }
  void rhoc(const double *mfrac) { for (int i=0; i<ns; i++) x[iSpec+i]=mfrac[i]; }
  void rhoc(const double &val)  { for (int i=0; i<ns; i++) x[iSpec+i]=val; }
  double& rhoc(const int&i)  { return x[iSpec+i]; }

  //!< Soot scalars times density
  const double* rhosc(void) const { return &x[iSoot]; }
  double rhosc(const int&i) const { return x[iSoot+i]; }
  double* rhosc(void) { return &x[iSoot]; }
  void rhosc(const double *s) { for (int i=0; i<nsc; i++) x[iSoot+i]=s[i]; }
  void rhosc(const double &val)  { for (int i=0; i<nsc; i++) x[iSoot+i]=val; }
  double& rhosc(const int&i)  { return x[iSoot+i]; }
  //@}


  /******************** Static Functions ****************************/
  //! return the number of variables - number of species
  static int NumVarSansSpecies(void) { return NUM_FLAME2D_VAR_SANS_SPECIES; }

  //! return the number of variables
  static int NumVar(void) { return n; }

  //! return the number of species
  static int NumSpecies(void) { return ns; }

  //! is the mixture reacting
  static bool isReacting(void) { return reacting; }

  //! is the mixture sooting
  static bool isSooting(void) { return soot_flag; }

  //! set static variables
  static void setNonReacting(void) { reacting = false; };
  static void set_gravity(const double &g);
  static void set_Mref(const double &Ma) { Mref = Ma; };

  //! print var list
  static void outputTecplotVarList(ostream &out,
				   const string &who, 
				   const string &prefix = "");

  /************************ Helper Functions ************************/
  void add( const Flame2D_State &U, const double &mult=1.0);
  void set( const Flame2D_State &U, const double &mult=1.0);
  void Average( const Flame2D_State &U1, const Flame2D_State &U2 );
  void Delta(const Flame2D_State& Ur, const Flame2D_State& Ul);
  void DeltaU(const Flame2D_pState& Wr, const Flame2D_pState& Wl);
  void Rotate(const Vector2D &norm_dir);
  void RotateBack(const Vector2D &norm_dir);
  void AdjustDensity(void);
  void ForceSpecSumToZero(void);
  void ForceSpecMin(void);

  /************************* Fluxes *********************************/
  //! Harten entropy fixes
  void HartenFix_Pos(const Flame2D_State &lambdas_a, 
		     const Flame2D_State &lambdas_l,
		     const Flame2D_State &lambdas_r);
  void HartenFix_Neg(const Flame2D_State &lambdas_a, 
		     const Flame2D_State &lambdas_l,
		     const Flame2D_State &lambdas_r);
  void HartenFix_Abs(const Flame2D_State &lambdas_a, 
		     const Flame2D_State &lambdas_l,
		     const Flame2D_State &lambdas_r);
  //! HLLE flux functions
  void FluxHLLE_x(const Flame2D_pState &Wl, const Flame2D_pState &Wr);
  void FluxHLLE_n(Flame2D_pState &Wl, 
		  Flame2D_pState &Wr,
		  const Vector2D &norm_dir);
  //! Linde flux functions
  void FluxLinde(const Flame2D_pState &Wl, const Flame2D_pState &Wr);
  void FluxLinde_n(Flame2D_pState &Wl, 
		   Flame2D_pState &Wr,
		   const Vector2D &norm_dir);
  //! Roe Flux Functions
  void FluxRoe_x(Flame2D_pState &Wl,
		 Flame2D_pState &Wr,
		 const int &Preconditioning,
		 const int &flow_type_flag,
		 const double &deltax);
  void FluxRoe_n(Flame2D_pState &Wl,
		 Flame2D_pState &Wr,
		 const Vector2D &norm_dir,
		 const int &Preconditioning,
		 const int &flow_type_flag,
		 const double &delta_n );
  //! AUSM flux functions
  void FluxAUSMplus_up(const Flame2D_pState &Wl,
		       const Flame2D_pState &Wr);
  void FluxAUSMplus_up_n(Flame2D_pState &Wl,
			 Flame2D_pState &Wr,
			 const Vector2D &norm_dir);
  //! Viscous Flux Functions
  void Viscous_Flux_n(Flame2D_pState &W,
		      const Flame2D_State &dWdx,
		      const Flame2D_State &dWdy,
		      const int Axisymmetric,
		      const Vector2D X,
		      const Vector2D &norm_dir, 
		      const double &mult=1.0);
  void Viscous_FluxHybrid_n(Flame2D_pState &W,
			    Flame2D_State &dWdx, 
			    Flame2D_State &dWdy,
			    const Vector2D &X,
			    const Flame2D_pState &Wl,
			    const Flame2D_State &dWdx_l,
			    const Flame2D_State &dWdy_l,
			    const Vector2D &Xl,
			    const Flame2D_pState &Wr,
			    const Flame2D_State &dWdx_r,
			    const Flame2D_State &dWdy_r,
			    const Vector2D &Xr,
			    const int &Axisymmetric,
			    const Vector2D &norm_dir,
			    const double &mult=1.0);

  /*********************** Checking *********************************/
  bool isPhysical(const int &harshness);
  bool speciesOK(const int &harshness);

  /********************* Operators Overloading **********************/
  // Index operator
  double &operator[](int index) { return x[index-1]; };
  const double &operator[](int index) const { return x[index-1]; };

  // Binary arithmetic operators.
  Flame2D_State operator +(const Flame2D_State &X) const;
  Flame2D_State operator -(const Flame2D_State &X) const;
  Flame2D_State operator *(const double &a) const;
  friend Flame2D_State operator *(const double &a, const Flame2D_State &X);
  Flame2D_State operator /(const double &a) const;
  double operator *(const Flame2D_State &X) const;
  Flame2D_State operator ^(const Flame2D_State &X) const;

  // Assignment Operator.
  Flame2D_State& operator =(const Flame2D_State &X); 

  // Shortcut arithmetic operators.
  Flame2D_State& operator +=(const Flame2D_State &X);
  Flame2D_State& operator -=(const Flame2D_State &X);
 
  // Unary arithmetic operators.
  //Flame2D_pState operator +(const Flame2D_pState &W);
  friend Flame2D_State operator -(const Flame2D_State &X);

  // Input-output operators.
  friend ostream& operator << (ostream &out_file, const Flame2D_State &X);
  friend istream& operator >> (istream &in_file,  Flame2D_State &X);


  /********************* Private Objects ****************************/
protected:

#ifdef STATIC_NUM_FLAME2D_VAR
  double                  x[STATIC_NUM_FLAME2D_VAR];
#else
  double                 *x; //!< solution state array
#endif
  static int              n; //!< Total number of vars
  static int             ns; //!< Number of species
  static int            nsc; //!< Number of soot scalars
  static int           ngas; //!< Number of variables associated with the gas phase
  static bool      reacting; //!< boolean indicating whether gas is reacting
  static double        Mref; //!< Mref for Precondtioning (normally set to incoming freestream Mach)
  static double   gravity_z; //!< m/s^2 acceleration due to gravity  
  static int      soot_flag; //!< flag indicating soot model
  static double*          y; //!< temporary storage for mass fractions

  //! solution state array indices
  static const int iRho   = 0; //!< density
  static const int iVx    = 1; //!< vx / rhovx 
  static const int iVy    = 2; //!< vy / rhovy
  static const int iPress = 3; //!< Pressure / energy
  static const int iSpec  = 4; //!< species array index
  static       int iSoot;      //!< soot scalar array index
    
};



/*!
 * Class: Flame2D_pState
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
class Flame2D_pState : public Flame2D_State {


  /********************* Private Objects ****************************/
private:
  Mixture Mix;         //!< Mixture Object
  static double* r;    //!< Temporary storage for reaction rates
  static double* h_i;  //!< Temporary storage for h(i)
  static double* cp_i; //!< Temporary storage for cp(i)
  static Soot2D_State  soot; //!< soot model object

  /********************** Constructors/Destructors ******************/
public:

  Flame2D_pState() {
    rho()=DENSITY_STDATM; p()=PRESSURE_STDATM;
    vx()=0.0; vy()=0.0; c(1.0/ns);
  }

  Flame2D_pState(const double &d, const Vector2D &V, 
		 const double &pre) 
  { rho()=d; p()=pre; v(V); c(1.0/ns); sc(0.0); setGas(); }

  Flame2D_pState(const double &d, const double &vvx, 
		 const double &vvy, const double &pre)
  { rho()=d; p()=pre; vx()=vvx; vy()=vvy; c(1.0/ns); sc(0.0); setGas(); }
  
  Flame2D_pState(const double &d, const double &vvx, 
		 const double &vvy, const double &pre, 
		 const double *mfrac) 
  { rho()=d; p()=pre; vx()=vvx; vy()=vvy; c(mfrac); sc(0.0); setGas(); }
  
  Flame2D_pState(const double &d, const Vector2D &V, 
		 const double &pre, const double *mfrac) 
  { rho()=d; p()=pre; v(V); c(mfrac); sc(0.0); setGas(); }

  Flame2D_pState(const Flame2D_pState &W) { Copy(W); }

  Flame2D_pState(const Flame2D_State &W) { 
    Copy(W);
  }

  ~Flame2D_pState() { Deallocate(); }

  void Copy( const Flame2D_pState &W ) { 
    for (int i=0; i<n; i++) x[i] = W.x[i];
    Mix = W.Mix;
  }
  void Copy( const Flame2D_State &W ) { 
    for(int i=0; i<n; i++) x[i] = W[i+1];
    setGas();
  }

  /******************* Static Functions *****************************/
public:

  //! initial mixture setup function
  static void setMixture(const string &mech_name,
			 const string &mech_file,
			 const int soot_model=0);
  //! set constant schmidt number
  static void setConstantSchmidt(const double* Sc) { Mixture::setConstantSchmidt(Sc); };
  //! Static memory allocator/deallocator  
  static void AllocateStatic(void);
  static void DeallocateStatic(void);
  //! The species names
  static string speciesName(const int&i) { return Mixture::speciesName(i); };
  static void speciesNames(string *names) { Mixture::getSpeciesNames(names); };
  //! The species index
  static int speciesIndex(const string&name) { return Mixture::speciesIndex(name); };
  //! The mechanism name
  static string mechName(void) { return Mixture::mechName(); };


  /********************* Accessors **********************************/
public:

  //!< Density.
  double rho(void) const { return x[iRho]; }

  //!< Flow velocity (2D)
  double vx(void) const { return x[iVx]; }
  double vy(void) const { return x[iVy]; }
  double vsqr(void) const { return x[iVx]*x[iVx] + x[iVy]*x[iVy]; };
  double vabs(void) const { return sqrt(x[iVx]*x[iVx] + x[iVy]*x[iVy]); };
  Vector2D v(void) const { return Vector2D(x[iVx], x[iVy]); };

  //!< Pressure.
  double p(void) const { return x[iPress]; }

  //!< Species mass fractions
  const double* c(void) const { return &x[iSpec]; }
  double c(const int&i) const { return x[iSpec+i]; }

  //!< Soot scalars
  const double* sc(void) const { return &x[iSoot]; }
  double sc(const int&i) const { return x[iSoot+i]; }

  //!< Momentum (2D)
  double rhovx(void) const { return rho()*vx(); };
  double rhovy(void) const { return rho()*vy(); };
  double rhovsqr(void) const { return rhovx()*rhovx()+rhovy()*rhovy(); }
  Vector2D rhov(void) const { return Vector2D(rhovx(), rhovy()); };

  //!< Species mass fractions times density
  double rhoc(const int&i) const { return rho()*c(i); };

  //!< Soot scalars times density
  double rhosc(const int&i) const { return rho()*sc(i); }

private:

  //!< Density.
  double& rho(void) { return x[iRho]; }

  //!< Flow velocity (2D)
  double& vx(void) { return x[iVx]; }
  double& vy(void) { return x[iVy]; }
  void v(const Vector2D &V) { x[iVx]=V.x; x[iVy]=V.y; };

  //!< Pressure.
  double& p(void) { return x[iPress]; }

  //!< Species mass fractions
  double* c(void) { return &x[iSpec]; }
  void c(const double *mfrac) { for (int i=0; i<ns; i++) x[iSpec+i]=mfrac[i]; }
  void c(const double &val)  { for (int i=0; i<ns; i++) x[iSpec+i]=val; }
  double& c(const int&i)  { return x[iSpec+i]; }

  //!< Soot scalars
  double* sc(void) { return &x[iSoot]; }
  void sc(const double *s) { for (int i=0; i<nsc; i++) x[iSoot+i]=s[i]; }
  void sc(const double &val)  { for (int i=0; i<nsc; i++) x[iSoot+i]=val; }
  double& sc(const int&i)  { return x[iSoot+i]; }

  //@{ @name Conserved variables:
  //!< Note: override these so they don't work:

  //!< Momentum (2D)
  double& rhovx(void) { assert(0); }
  double& rhovy(void) { assert(0); }

  //!< Total Energy (rho *(e + HALF*v^2))
  double& E(void) { assert(0); }

  //!< Species mass fractions times density
  const double* rhoc(void) const { assert(0); }
  double* rhoc(void) { assert(0); }
  void rhoc(const double *mfrac) { assert(0); }
  void rhoc(const double &val)  { assert(0); }
  double& rhoc(const int&i)  { assert(0); }

  //!< Soot scalars times density
  const double* rhosc(void) const { assert(0); }
  double* rhosc(void) { assert(0); }
  void rhosc(const double *s) { assert(0); }
  void rhosc(const double &val)  { assert(0); }
  double& rhosc(const int&i)  { assert(0); }

  //@}

  /******************** State Setup Functions ***********************/
public:

  //!< Set gas state throught temperature, pressure, mass fracs
  void setState_TPY(const double &Temp, const double &Press, const double *y);
  //!< Set velocity
  void setVelocity(const double &vvx, const double &vvy) { vx() = vvx; vy() = vvy; };
  void setVelocityX(const double &vvx) { vx() = vvx; };
  void setVelocityY(const double &vvy) { vy() = vvy; };
  //!< Set density through temperature
  void setTemperature(const double &Temp);
  //!< Set pressure through various means
  void setPressure(const double &Press);
  void setEnergy(const double &en);
  void setEnthalpy(const double &h);
  //!< Equilibrium composition
  void equilibrate_HP(void) { Mix.equilibrate_HP(rho(), p(), c()); };
  void equilibrate_TP(void) { Mix.equilibrate_TP(rho(), p(), c()); };

private:

  //!< Set mixture object
  void setGas(void) { Mix.setState_DPY(rho(), p(), c()); };

  /*************** Primitive / Conserved Transformation *************/
public:

  //!< conserved to primitive transforation
  void setU(const Flame2D_State &U);
  void setW(const Flame2D_pState &W){ if( this != &W)  Copy(W); }
  void setW(const Flame2D_State &W){ Copy(W); }
  void getU(Flame2D_State &U) const;
  Flame2D_State U(void) const;

  /******************** Misc Quantities of Interest *****************/
  //!< Total Energy (rho *(e + HALF*v^2))
  double E(void) const { return (rho()*(e() + 0.5*vsqr())); };
  //!< Total Enthalpy ( H = h + velocity^2 )
  double H(void) const { return (rho()*(h() + 0.5*vsqr())); };
  double Hs(void) const { return (rho()*(hs() + 0.5*vsqr())); };
  //!< Speed of sound 
  double a(void) const { return sqrt(g()*Rtot()*T()); }

  /****************** Mixture Object Wrappers ***********************/
  //!< Temperature [K]
  double T(void) const { return Mix.temperature(); };
  //!< Molar mass [kg/kmole]
  double MW(void) const { return Mix.molarMass(); };
  static void MW(double*MWs) { Mixture::getMolarMasses(MWs); };
  //!< Gas constant [ J/(kg K) ]
  double Rtot(void) const { return Mix.gasConstant(); };
  //!< Internal Energy (et = es + echem) [ J/kg ]
  double e(void) const { return Mix.internalEnergy(); };
  double es(void) const { return Mix.internalEnergySens(); };
  //!< enthalpy  (ht = hs + hf) [ J/kg ]
  double h(void) const { return Mix.enthalpy(); };
  double h(const double &Temp) const { return Mix.enthalpy(Temp, p(), c()); };
  double hs(void) const { return Mix.enthalpySens(); };
  void h(double*hi) const { Mix.getEnthalpy(p(), c(), hi); };
  //!< heat capacity [ J / (kg K) ]
  double hprime(void) const { return Mix.heatCapacity_p(); };
  double Cp(const double &Temp) const { return Mix.heatCapacity_p(Temp, p(), c()); };
  double Cp(void) const { return Mix.heatCapacity_p(); };
  double Cv(void) const { return Mix.heatCapacity_v(); };
  void Cp(double*cpi) const { Mix.getHeatCapacity_p(p(), c(), cpi); };
  //!< species heat capacities [J/(kg K)] and enthalpies (ht = hs + hf) [ J/kg ]
  void Cp_and_h( double*cp, double*h ) const {
    Mix.get_cp_and_h( p(), c(), cp, h );
  }
  //! ideal gas ratio
  double g(void) const { return Mix.heatRatio(); };
  //!< dynamic viscosity [ kg m^-1 s^-1]
  double mu(void) const { return Mix.viscosity(); };
  //! Mixture Thermal conducitivity [W/(m K)]
  double kappa(void) const { return Mix.thermalCond(); };
  //! Mixture-averaged diffusion coefficients [m^2/s].
  double Diffusion_coef(const int &i) const { return Mix.speciesDiffCoef(i); };
  //!< Dimensionaless Quantities
  double Sc(const int &i) const { return Mix.schmidt(rho(), i); };
  double Pr(void) const { return Mix.prandtl(); };
  double Le(const int &i) const { return Mix.lewis(rho(), i); };
  double Re(const double &l) { updateViscosity(); return rho()*vabs()*l/mu(); };
  //!< Derivatives
  double diedip() const { return (hprime() - Rtot())/(rho()*Rtot()); };
  double diedirho() const { return -p()*(hprime() - Rtot())/(rho()*rho()*Rtot()); };
  double dihdic(const int &i) const { return Mix.DihdDiy(i); };
  double Phi(void) const { return Mix.Phi(); };
  //! update viscosity
  void updateViscosity(void) { Mix.updateViscosity(rho(), c()); };
  //! update related transport properties -> kappa and D_i
  void updateTransport(void) { Mix.updateTransport(rho(), c()); };
  //! update derivative of species enthalpies wrt mass fracs
  double updateDihdDic(void) { Mix.updateDihdDic( rho(), c() ); };

  /************ Strain rate tensor, laminar stress tensor ***********/
  void Strain_Rate(const Flame2D_State &dWdx,
		   const Flame2D_State &dWdy,
		   const int Axisymmetric,
		   const Vector2D &X,
		   Tensor2D &strain_rate) const;
  void Laminar_Stress(const Flame2D_State &dWdx,
		      const Flame2D_State &dWdy,
		      const int Axisymmetric,
		      const Vector2D &X,
		      Tensor2D &laminar_stress);
  void Viscous_Quantities(const Flame2D_State &dWdx,
			  const Flame2D_State &dWdy,
			  const int Axisymmetric,
			  const Vector2D &X,
			  Vector2D &qflux,
			  Tensor2D &tau,
			  Vector2D &Vcorr);
  double WallShearStress(const Vector2D &X1,
			 const Vector2D &X2,
			 const Vector2D &X3,
			 const Vector2D &norm_dir);
  Vector2D DiffusionVel( const double &dcdx,
			 const double &dcdy,
			 const int &k ) const;

  /*********************** Helper Functions *************************/
  void Reconstruct( const Flame2D_pState &Wc, 
		    const Flame2D_State &phi, 
		    const Flame2D_State &dWdx, 
		    const Flame2D_State &dWdy,
		    const Vector2D &dX, 
		    const double &mult=1.0);
  void Average( const Flame2D_State &W1, const Flame2D_State &W2 );
  void Average( const Flame2D_pState &W1, const Flame2D_pState &W2 );

  /********************* Eigenvalues / Eigenvectors *****************/
  //! No preconditioning
  void lambda_x(Flame2D_State &lambdas) const;
  void rc_x(const int &index, Flame2D_State& rc);
  void lp_x(const int &index, Flame2D_State& lp) const;
  void Flux_Dissipation(const int &i,
			const Flame2D_State &dWrl, 
			const double &wavespeed, 
			Flame2D_State &Flux,
			const double &mult=1.0);
  void Flux_Dissipation_Jac(DenseMatrix &Jac,
			    const Flame2D_State &wavespeeds, 
			    const double &mult=1.0);
  //! Weiss-Smith Preconditioning
  void lambda_preconditioned_x(Flame2D_State &lambdas, const double &MR2) const;
  void rc_x_precon(const int &index, 
		   const double &MR2, 
		   const double &uprimed, 
		   const double &cprimed, 
		   Flame2D_State& rc);
  void lp_x_precon(const int &index, 
		   const double &MR2, 
		   const double &uprimed, 
		   const double &cprimed, 
		   Flame2D_State& lp) const;
  void Flux_Dissipation_precon(const double &MR2, 
			       const Flame2D_State &dWrl, 
			       const Flame2D_State &wavespeeds, 
			       Flame2D_State &Flux_dissipation);
  void Flux_Dissipation_Jac_precon(DenseMatrix &Jac,
				   const double &MR2,
				   const Flame2D_State &wavespeeds, 
				   const double &mult=1.0);

  /*********************** Fluxes ***********************************/
  //! Inviscid flux functions
  void Fx(Flame2D_State &FluxX, const double& mult=1.0) const;
  Flame2D_State Fx(void) const;
  void addFx(Flame2D_State &FluxX, const double& mult=1.0) const;
  //! Roe Average state
  void RoeAverage(const Flame2D_pState &Wl, const Flame2D_pState &Wr);
  //! HLLE wavespeeds
  static void HLLE_wavespeeds(Flame2D_pState &Wl,
			      Flame2D_pState &Wr,
			      const Vector2D &norm_dir,
			      Vector2D &wavespeeds);
  //! Viscous flux
  void Viscous_Flux_x(const Flame2D_State &dWdx,
		      const Vector2D &qflux,
		      const Tensor2D &tau,
		      const Vector2D &Vcorr,
		      Flame2D_State &Flux, 
		      const double& mult=1.0) const;
  void Viscous_Flux_y(const Flame2D_State &dWdy,
		      const Vector2D &qflux,
		      const Tensor2D &tau,
		      const Vector2D &Vcorr,
		      Flame2D_State &Flux, 
		      const double& mult=1.0) const;

  /************************ Flux Jacobians **************************/
  void dWdU(DenseMatrix &dWdQ);
  void dWdU_FD(DenseMatrix &dWdQ) const;
  void dFIdU(DenseMatrix &dFdU);
  void dFIdU_FD(DenseMatrix &dFdU) const;
  void dFIdW(DenseMatrix &dFdW, const double& mult=1.0);
  void dFIdW_FD(DenseMatrix &dFdW, const double& mult=1.0) const;
  void dFvdWf_dGvdWf( DenseMatrix &dFvdWf, 
		      DenseMatrix &dGvdWf, 
		      const Flame2D_State &dWdx, 
		      const Flame2D_State &dWdy, 
		      const int &Axisymmetric, 
		      const double &radius );

  /********************** Source Terms ******************************/
  //! Inviscid component of axisymetric source term
  void Sa_inviscid(Flame2D_State &S, const Vector2D &X, 
		   const int Axisymmetric, const double& mult=1.0) const;
  void dSa_idU( DenseMatrix &dSa_IdU, 
		const Vector2D &X, 
		const int Axisymmetric );
  //! Viscous component of axisymetric source term
  void Sa_viscous(Flame2D_State &S, 
		  const Flame2D_State &dWdx,
		  const Flame2D_State &dWdy,
		  const Vector2D &X, 
		  const int Axisymmetric, 
		  const double& mult=1.0);
  void dSa_vdW( DenseMatrix &dSa_VdW,
		const Flame2D_State &dWdx,
		const Flame2D_State &dWdy,
		const Vector2D &X, 
		const int Axisymmetric,
		const double &d_dWdx_dW, 
		const double &d_dWdy_dW);
  
  //! Source terms associated with finite-rate chemistry
  void Sw(Flame2D_State &S, const double& mult=1.0) const;
  void dSwdU(DenseMatrix &dSdU) const;
  double dSwdU_max_diagonal(void) const;

  //! Source terms associated with gravitational forces
  void Sg(Flame2D_State &S, const double& mult=1.0) const;
  void dSgdU(DenseMatrix &dSgdU) const;

  //! Source terms associated with soot
  void Ssoot(Flame2D_State &S, const double& mult=1.0) const;

  /********************* Preconditioning ****************************/
  double u_plus_aprecon(const double &u,const int &flow_type_flag,
			const double &deltax);
  void u_a_precon(const double &UR,double &uprimed, double &cprimed) const;
  double Mr2(const int &flow_type_flag, const double &deltax);
  void Low_Mach_Number_Preconditioner(DenseMatrix &P,
				      const int &Viscous_flag, 
				      const double &deltax );
  void Low_Mach_Number_Preconditioner_Inverse(DenseMatrix &Pinv,	
					      const int &Viscous_flag, 
					      const double &deltax );
  
  /******************** Boundary Conditions *************************/
  void Reflect( const Flame2D_pState &W,
		const Vector2D &norm_dir );
  void Free_Slip( const Flame2D_pState &Win,
		  const Flame2D_pState &Wout,
		  const Vector2D &norm_dir,
		  const int &TEMPERATURE_BC_FLAG );
  void No_Slip(const Flame2D_pState &Win,
	       const Flame2D_pState &Wout,
	       const Vector2D &norm_dir,
	       const int &TEMPERATURE_BC_FLAG);
  void Moving_Wall(const Flame2D_pState &Win,
		   const Flame2D_pState &Wout,
		   const Vector2D &norm_dir, 
		   const double &wall_velocity,
		   const int &TEMPERATURE_BC_FLAG);
  void BC_1DFlame_Inflow(const Flame2D_pState &Wi,
			 const Flame2D_pState &Wo,
			 const Flame2D_pState &Woutlet,
			 const Vector2D &norm_dir);
  void BC_1DFlame_Outflow(const Flame2D_pState &Wi,
			  const Flame2D_pState &Wo,
			  const Flame2D_pState &Winlet,
			  const Vector2D &norm_dir);
  void BC_2DFlame_Inflow(const Flame2D_pState &Wi,
			 const Flame2D_pState &Wo,
			 const Vector2D &norm_dir);
  void BC_2DFlame_Outflow(const Flame2D_pState &Wi, 
			  const Flame2D_pState &Wo,
			  const Vector2D &norm_dir);
  void BC_Characteristic_Pressure(const Flame2D_pState &Wi,
				  const Flame2D_pState &Wo,
				  const Vector2D &norm_dir);

  /****************** Jump Conditions *******************************/
  void FlameJumpLowMach(const Flame2D_pState &Wu);

  /********************** Exact Solutions ***************************/
  void ViscousChannelFlow(const Vector2D X,
			  const double Vwall,
			  const double dp);
  void FlatPlate(const Flame2D_pState &Winf,
		 const Vector2D X,
		 double &eta,
		 double &f,
		 double &fp,
		 double &fpp);

  /************************* Checking *******************************/
  //! compute sum of mass fractions
  double SpecSum(void) const;
  //! compute sum of reaction rates
  double OmegaSum(void) const;
  //! compute sum of diffusion velocities
  Vector2D DiffSum(const double* dcdx, 
		   const double* dcdy) const;

private:
  //! overload so can't be used
  bool isPhysical(const int &harshness)  {assert(0);};
  bool speciesOK(const int &harshness) {assert(0);};
    
  /******************** Operators Overloading ***********************/
public:
  // Assignment Operator.
  Flame2D_pState& operator =(const Flame2D_pState &W); 
  Flame2D_pState& operator =(const Flame2D_State &W); 

  // Input-output operators.
  friend ostream& operator << (ostream &out_file, const Flame2D_pState &X);
  friend istream& operator >> (istream &in_file,  Flame2D_pState &X);

};



/////////////////////////////////////////////////////////////////////
/// Allocators / Deallocators
/////////////////////////////////////////////////////////////////////

/****************************************************
 * Allocator/deallocator for object memory.
 ****************************************************/
inline void Flame2D_State :: Allocate() {
#ifndef STATIC_NUM_FLAME2D_VAR
  Deallocate();
  if (n>0) { x = new double[n]; }
#endif
};

inline void Flame2D_State :: Deallocate() { 
#ifndef STATIC_NUM_FLAME2D_VAR
  if (x!=NULL) { delete[] x; x = NULL; } 
#endif
};

inline void Flame2D_State :: Nullify() { 
#ifndef STATIC_NUM_FLAME2D_VAR
  x = NULL;
#endif
};


/****************************************************
 * Allocator/deallocator for static memory.
 ****************************************************/
inline void Flame2D_pState :: AllocateStatic() {
  if (ns>0) { 
    r = new double[ns];
    h_i = new double[ns];
    cp_i = new double[ns];
    y = new double[ns];
  }
};

inline void Flame2D_pState :: DeallocateStatic() { 
  if (r!=NULL) { delete[] r; r = NULL; } 
  if (h_i!=NULL) { delete[] h_i; h_i = NULL; } 
  if (cp_i!=NULL) { delete[] cp_i; cp_i = NULL; } 
  if (y!=NULL) { delete[] y; y = NULL; } 
  Mixture::DeallocateStatic();
};


/////////////////////////////////////////////////////////////////////
/// Operator Overloads
/////////////////////////////////////////////////////////////////////

/********************************************************
 * Flame2D_State -- Binary arithmetic operators.        *
 ********************************************************/
//----------------- Addition -----------------------------//
inline Flame2D_State Flame2D_State::operator +(const Flame2D_State &X) const{    
  Flame2D_State Temp(*this);
  Temp += X;
  return Temp;
}

//------------------ Subtraction ------------------------//
inline Flame2D_State Flame2D_State::operator -(const Flame2D_State &X) const{
  Flame2D_State Temp(*this);
  Temp -= X;
  return Temp;
}

//---------------- Scalar Multiplication ------------------//
inline Flame2D_State Flame2D_State::operator *(const double &a) const{
  Flame2D_State Temp(*this);
  for( int i=0; i<n; i++) Temp.x[i] = x[i]*a; 
  return(Temp);
}

inline Flame2D_State operator *(const double &a, const Flame2D_State &X){
  Flame2D_State Temp;
  for( int i=0; i<X.n; i++) Temp.x[i] = X.x[i]*a;
  return(Temp);
}

//--------------- Scalar Division ------------------------//
inline Flame2D_State Flame2D_State::operator /(const double &a) const {
  Flame2D_State Temp(*this);
  for(int i=0; i<n; i++) Temp.x[i] = x[i]/a;
  return(Temp);
}

//----------------- Inner Product ------------------------//
inline double Flame2D_State::operator *(const Flame2D_State &X) const{
  double sum(0.0);
  for(int i=0; i<n; i++) sum += x[i]*X.x[i];
  return (sum);
}

//----------- solution state product operator ------------//
inline Flame2D_State Flame2D_State::operator ^( const Flame2D_State &X) const {
  Flame2D_State Temp(*this);
  for(int i=0; i<n; i++) Temp.x[i] = x[i]*X.x[i];
  return(Temp);
}

//----------------- Assignment ----------------------------//
inline Flame2D_State& Flame2D_State::operator =(const Flame2D_State &X){
  //self assignment protection
  if( this != &X){ 
    Copy(X);
  }
  return (*this);
}
inline Flame2D_pState& Flame2D_pState::operator =(const Flame2D_pState &W){
  if( this != &W) Copy(W);
  return (*this);
}
inline Flame2D_pState& Flame2D_pState::operator =(const Flame2D_State &W){
  Copy(W);
  return (*this);
}

/********************************************************
 * Flame2D_pState -- Shortcut arithmetic operators.     *
 ********************************************************/
inline Flame2D_State& Flame2D_State::operator +=(const Flame2D_State &X){
  for( int i=0; i<n; i++)  x[i] += X.x[i];
  return (*this);
}

inline Flame2D_State& Flame2D_State::operator -=(const Flame2D_State &X) {
  for(int i=0; i<n; i++) x[i] -= X.x[i];
  return (*this); 
}

/********************************************************
 * Flame2D_pState -- Unary arithmetic operators.         *
 ********************************************************/
inline Flame2D_State operator -(const Flame2D_State &X) {
  Flame2D_State Temp(X);
  for(int i=0; i<X.n; i++)  Temp.x[i] = -X.x[i]; 
  return(Temp);
}



/********************************************************
 * Flame2D_pState -- Input-output operators.            *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const Flame2D_State &X) {
  out_file.precision(10);
  out_file.setf(ios::scientific);
  for( int i=0; i<X.n; i++) out_file<<" "<<X.x[i];
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, Flame2D_State &X) {
  in_file.setf(ios::skipws);
  for( int i=0; i<X.n; i++) in_file>>X.x[i];
  in_file.unsetf(ios::skipws);
  return (in_file);
}

inline ostream &operator << (ostream &out_file, const Flame2D_pState &W) {
  out_file.precision(10);
  out_file.setf(ios::scientific);
  for( int i=0; i<W.NumVar(); i++) out_file<<" "<<W[i+1];
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, Flame2D_pState &W) {
  in_file.setf(ios::skipws);
  for( int i=0; i<W.NumVar(); i++) in_file>>W[i+1];
  in_file.unsetf(ios::skipws);
  W.setGas();
  return (in_file);
}


/////////////////////////////////////////////////////////////////////
/// Helper Functions
/////////////////////////////////////////////////////////////////////

/****************************************************
 * Flame2D_State::add() -- Add state to current one
 ****************************************************/
inline void Flame2D_State::add( const Flame2D_State &U, 
				const double &mult) {
  for (int i=0; i<n; i++) x[i] += mult*U.x[i];
}

/****************************************************
 * Flame2D_State::set() -- Set current state
 ****************************************************/
inline void Flame2D_State::set( const Flame2D_State &U, 
				const double &mult) {
  for (int i=0; i<n; i++) x[i] = mult*U.x[i];
}

/****************************************************
 * Flame2D_State::Average() -- Set average of two states
 ****************************************************/
inline void Flame2D_State::Average( const Flame2D_State &U1, 
				    const Flame2D_State &U2 ) {
  for (int i=0; i<n; i++) x[i] = 0.5*(U1.x[i]+U2.x[i]);
}

/****************************************************
 * Flame2D_State::Delta() -- Set difference of two states
 ****************************************************/
inline void Flame2D_State::Delta(const Flame2D_State& Ur, 
				 const Flame2D_State& Ul) {
  for (int i=0; i<n; i++) x[i] = Ur.x[i] - Ul.x[i];
}


/********************************************************
 * Flame2D_State -- DeltaU: Change in conserved solution*
 ********************************************************/
inline void Flame2D_State::DeltaU(const Flame2D_pState& Wr, 
				  const Flame2D_pState& Wl) {
  rho() = Wr.rho() - Wl.rho();
  rhovx() = Wr.rhovx() - Wl.rhovx();
  rhovy() = Wr.rhovy() - Wl.rhovy();
  E() = Wr.E() - Wl.E();
  for (int i=0; i<ns; i++) rhoc(i) = Wr.rhoc(i) - Wl.rhoc(i);
  for (int i=0; i<nsc; i++) rhosc(i) = Wr.rhosc(i) - Wl.rhosc(i);
}

/********************************************************
 * Flame2D_State -- Rotate: Rotate the solution state   *
 ********************************************************/
inline void Flame2D_State::Rotate( const Vector2D &norm_dir ) {
  double u(vx()), v(vy());
  vx() =    u*norm_dir.x + v*norm_dir.y;
  vy() =  - u*norm_dir.y + v*norm_dir.x;
}

inline void Flame2D_State::RotateBack( const Vector2D &norm_dir ) {
  double ur(vx()), vr(vy());
  vx() = ur*norm_dir.x - vr*norm_dir.y;
  vy() = ur*norm_dir.y + vr*norm_dir.x;
}

/****************************************************
 * Flame2D_pState::Reconstruct()
 *
 * First order piecewise linear solution reconstruction
 ****************************************************/
inline void Flame2D_pState::Reconstruct( const Flame2D_pState &Wc, 
					 const Flame2D_State &phi, 
					 const Flame2D_State &dWdx, 
					 const Flame2D_State &dWdy,
					 const Vector2D &dX, 
					 const double &mult) {
  for (int i=0; i<n; i++)
    x[i] = Wc.x[i] + mult*(phi[i+1]*dWdx[i+1]*dX.x + phi[i+1]*dWdy[i+1]*dX.y);
  setGas();
}

/****************************************************
 * Flame2D_pState::Average()
 *
 * Set the state as the average of two states
 ****************************************************/
inline void Flame2D_pState::Average( const Flame2D_State &W1, const Flame2D_State &W2 ) {
  for (int i=0; i<n; i++) x[i] = 0.5*(W1[i+1]+W2[i+1]);
  setGas();
}

inline void Flame2D_pState::Average( const Flame2D_pState &W1, const Flame2D_pState &W2 ) {
  for (int i=0; i<n; i++) x[i] = 0.5*(W1.x[i]+W2.x[i]);
  setGas();
}

/****************************************************
 * Adjust density to be consistent
 ****************************************************/
inline void Flame2D_State::AdjustDensity(void) { 
  double sum(0.0); 
  for (int i=0; i<ns; i++) sum += rhoc(i);
  rho() = sum;
}

/****************************************************
 * Divide error among species, forcing sum to zero.
 ****************************************************/
inline void Flame2D_State::ForceSpecSumToZero(void) { 
  double sum( 0.0 );
  for (int i=0; i<ns; i++) sum += c(i);
  sum /= ns;
  for (int i=0; i<ns; i++) c(i) -= sum;
}


/****************************************************
 * Set species values to minimum.
 ****************************************************/
inline void Flame2D_State::ForceSpecMin(void) { 
  double cmin( MILLION );
  for (int i=0; i<ns; i++) cmin = min( cmin, c(i) );
  for (int i=0; i<ns; i++) c(i) = cmin;
}


/////////////////////////////////////////////////////////////////////
/// Primitive State Setup Functions
/////////////////////////////////////////////////////////////////////

/****************************************************
 * Flame2D_pState::setState_TPY()
 *
 * Set the state through temperature, pressure and 
 * mass fractions
 ****************************************************/
inline void Flame2D_pState::setState_TPY(const double &Temp, 
					 const double &Press, 
					 const double *y) {    
  p() = Press;
  c( y );
  Mix.setState_TPY(Temp, Press, y); 
  rho() = p() / (Mix.gasConstant()*Temp);
};

/****************************************************
 * Flame2D_pState::setTemperature()
 *
 * Change the state temperature -> changes density
 * Temperature in [K]
 ****************************************************/
inline void Flame2D_pState::setTemperature(const double &Temp) {
  Mix.setState_TPY(Temp, p(), c());
  rho() = p() / (Mix.gasConstant()*Temp);
};

/****************************************************
 * Flame2D_pState::setPressure()
 *
 * Change the state pressure [ Pressure in Pa ]
 ****************************************************/
inline void Flame2D_pState::setPressure(const double &Press) {    
  p() = Press;
  setGas(); 
};

/****************************************************
 * Flame2D_pState::setEnergy()
 *
 * Set the state energy -> changes pressure
 * Total Energy (es + ef) in [J/kg]
 ****************************************************/
inline void Flame2D_pState::setEnergy(const double &en) { 
  Mix.setState_DEY(rho(), en, c()); 
  p() = rho()*Mix.gasConstant()*Mix.temperature();
};

/****************************************************
 * Flame2D_pState::setEnthalpy()
 * 
 * Set the state enthalpy -> changes pressure
 * Total enthalpy (hs + hf) in [J/kg]
 ****************************************************/
inline void Flame2D_pState::setEnthalpy(const double &h) { 
  Mix.setState_DHY(rho(), h, c()); 
  p() = rho()*Mix.gasConstant()*Mix.temperature();
};


/////////////////////////////////////////////////////////////////////
/// Primitive / Conserved Transformation
/////////////////////////////////////////////////////////////////////
/****************************************************
 * Flame2D_pState::setU()
 *
 * Various functions to set the primitive state from 
 * the conserved state.
 ****************************************************/
inline void Flame2D_pState::setU(const Flame2D_State &U) {
  rho() = U.rho();
  vx() = U.rhovx()/U.rho();
  vy() = U.rhovy()/U.rho();
  for(int i=0; i<ns; i++) c(i) = U.rhoc(i)/U.rho();
  for(int i=0; i<nsc; i++) sc(i) = U.rhosc(i)/U.rho();
  double e( U.E()/U.rho() - 0.5*vsqr() );
  setEnergy(e);
}

/****************************************************
 * Flame2D_pState::getU()
 *
 * Various functions to get the conserved state from 
 * the primitive state.
 ****************************************************/
inline void Flame2D_pState::getU(Flame2D_State &U) const {
  U.rho() = rho();
  U.rhovx() = rhovx();
  U.rhovy() = rhovy();
  U.E() = E();
  for(int i=0; i<ns; i++) U.rhoc(i) = rhoc(i);
  for(int i=0; i<nsc; i++) U.rhosc(i) = rhosc(i);
}

inline Flame2D_State Flame2D_pState::U(void) const {
  static Flame2D_State U;
  getU(U);
  return U;
}


/////////////////////////////////////////////////////////////////////
/// Laminar Quantities
/////////////////////////////////////////////////////////////////////

/********************************************************
 * Flame2D_pState -- Strain rate tensor                 *
 ********************************************************/
inline void Flame2D_pState::Strain_Rate(const Flame2D_State &dWdx,
					const Flame2D_State &dWdy,
					const int Axisymmetric,
					const Vector2D &X,
					Tensor2D &strain_rate) const {
  double radius;

  // Strain rate (+dilatation)	
  double div_v( dWdx.vx() + dWdy.vy() );
  if (Axisymmetric == AXISYMMETRIC_X) {
    radius = (X.x < MICRO) ? MICRO : X.x;    //fabs(X.x) ??
    div_v += vx()/radius;
  } else if (Axisymmetric == AXISYMMETRIC_Y) {    
    radius = (X.y < MICRO) ? MICRO : X.y;
    div_v += vy()/radius;
  } 

  strain_rate.xx = dWdx.vx()-div_v/THREE;
  strain_rate.xy = HALF*(dWdx.vy() + dWdy.vx());
  strain_rate.yy = dWdy.vy()-div_v/THREE;

  if (Axisymmetric == PLANAR) {
    strain_rate.zz = -(strain_rate.xx + strain_rate.yy); 
  } else if (Axisymmetric == AXISYMMETRIC_X) {
    strain_rate.zz = vx()/radius-div_v/THREE;
  } else if (Axisymmetric == AXISYMMETRIC_Y) {
    strain_rate.zz = vy()/radius-div_v/THREE;
  } 
}

/********************************************************
 * Flame2D_pState -- Laminar (molecular) fluid stress   *
 ********************************************************/
inline void Flame2D_pState::Laminar_Stress(const Flame2D_State &dWdx,
					   const Flame2D_State &dWdy,
					   const int Axisymmetric,
					   const Vector2D &X,
					   Tensor2D &laminar_stress) {
  
  Strain_Rate(dWdx, dWdy, Axisymmetric, X, laminar_stress);

  updateViscosity();
  double f( TWO*mu() ); 
  laminar_stress.xx *= f;
  laminar_stress.xy *= f;
  laminar_stress.yy *= f;
  laminar_stress.zz *= f;
}



/******************************************************
 Compute viscous heat flux and stresses
*******************************************************/
inline void Flame2D_pState::Viscous_Quantities(const Flame2D_State &dWdx,
					       const Flame2D_State &dWdy,
					       const int Axisymmetric,
					       const Vector2D &X,
					       Vector2D &qflux,
					       Tensor2D &tau,
					       Vector2D &Vcorr) {

  // declares
  static Vector2D grad_T;
  updateTransport();

  /**************** Temperature gradient ***************************/
  /* Temperature gradients from using the chain rule 
     with the ideal gas law (P=rho*R*T) 
     dT/dx = 1/rho*R *( dP/dx - P/rho * drho/dx) */
  grad_T.x = (dWdx.p() - (p()/rho())*dWdx.rho());
  grad_T.y = (dWdy.p() - (p()/rho())*dWdy.rho());
  grad_T *= ONE/(rho()*Rtot());
  
  /***************** Molecular (Laminar) Stresses ******************/
  Laminar_Stress(dWdx, dWdy, Axisymmetric, X, tau);
  
  /************* Molecular (Laminar) Heat flux Vector **************/
  /****************** Thermal Conduction ***************************
     q = - kappa * grad(T)                                         */
  qflux = grad_T;
  qflux *= - kappa();
  /****************** Thermal Diffusion ****************************/
  // q -= rho * sum ( hi * Di *gradci)  
  h(h_i); // <- individual species enthalpies
  double tmp;
  for(int i=0; i<ns; i++){ 
    tmp = rho() * h_i[i] * Diffusion_coef(i);
    qflux.x -= tmp * dWdx.c(i);
    qflux.y -= tmp * dWdy.c(i);
  }
  /****************** Diffusion Correction Velocity ****************/
  //  Vc = \sum_{i=0}^N { D_k * dY_i/dx_i }
  Vcorr.zero();
  for(int i=0; i<ns; i++){
    tmp = Diffusion_coef(i);
    Vcorr.x += tmp*dWdx.c(i);
    Vcorr.y += tmp*dWdy.c(i);
  }

}

/******************************************************
 * Compute diffusion velocities for individual species
 ******************************************************/
inline Vector2D Flame2D_pState::DiffusionVel( const double &dcdx,
					      const double &dcdy,
					      const int &k ) const {
  Vector2D Vk( - Diffusion_coef(k) * dcdx, 
	       - Diffusion_coef(k) * dcdy );
  return Vk;
}



/////////////////////////////////////////////////////////////////////
/// Checks
/////////////////////////////////////////////////////////////////////

/**************************************************************
 * Check the sum of the species mass factions, sum(Y_k) = 1
 **************************************************************/
inline double Flame2D_pState::SpecSum(void) const {
  double sum( 0.0 );
  for ( int i=0; i<ns; i++) sum += c(i);
  return sum;
}

/**************************************************************
 * Check the sum of the species reaction rates, sum(omega_k) = 0
 **************************************************************/
inline double Flame2D_pState::OmegaSum(void) const { 
  double sum( 0.0 );
  Mix.getRates( p(), c(), r );
  for ( int i=0; i<ns; i++) sum += r[i];
  return sum;
}

/**************************************************************
 * Check the sum of the species diffusion velocities, sum(Y_k*V_k) = 0
 **************************************************************/
inline Vector2D Flame2D_pState::DiffSum(const double* dcdx, 
					const double* dcdy) const { 
  Vector2D Vsum(0.0, 0.0);
  for ( int i=0; i<ns; i++) { 
    Vsum += DiffusionVel( dcdx[i], dcdy[i], i );
  }
  return Vsum;
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

  This is meant for a conserved state only

***************************************************************/
inline bool Flame2D_State::speciesOK(const int &harshness) {

  // declares
  double yi;
  double sum(ZERO);

  //
  // Loop over the species
  //
  for(int i=0; i<ns; i++){


    yi = rhoc(i)/rho();

    //=================================================================
    // check for -ve
    //=================================================================
    if(yi < ZERO){

      //---------------------------------------------------------------
      // a. -> check for small -ve and set to ZERO 
      //---------------------------------------------------------------
      if(yi > -SPEC_TOLERANCE){
	rhoc(i) = ZERO;
	yi = ZERO;

      //---------------------------------------------------------------
      // b. ->  report error depending upon harshness 
      //---------------------------------------------------------------
      } else {
	
	// display error
#ifdef _DEBUG
	cout<<"\ncState -ve mass fraction in "<<speciesName(i)<<" "<<temp
	    <<" greater than allowed tolerance of "<<-SPEC_TOLERANCE; 
#endif

	// if harshness less than ten, return in error
	if( harshness < 10){
	  return false;

	// else, just force to zero
	} else { 
	  rhoc(i) = ZERO;
	  yi = ZERO;

	  // display error
#ifdef _DEBUG
	  cout<<"\ncState rhospec["<<i<<"] = " << rhoc(i) 
	      << " -ve mass fraction larger than tolerance,"
	      <<" but setting to zero and continuing anyway. ";
#endif
	} // endif - harshness

      //---------------------------------------------------------------
      } // endif - yi<-TOL
      //---------------------------------------------------------------

      
    //=================================================================
    } // endif - y<0
    //=================================================================

    // add the contribution to the sum
    sum += rhoc(i);

  } // enfor - species


  // Distribute error over all the species
  // sum /= rho();
  // for(int i=0; i<ns; i++) rhoc(i) /= sum;

  // Give error to density
  // if ( fabs(sum-rho_)>NANO ) rho() = sum;
  rho() = sum;

  // SUCCESS!
  return true;
}



/****************************************************
 * Check if the properties seam physical
 * This is meant for a conserved state only
 ****************************************************/
inline bool Flame2D_State::isPhysical(const int &harshness) {
  
  // check for nan's, inf's etc.... debugging !!!
#ifdef _DEBUG
  for( int i=0; i<n; i++){ 
    if( x[i] != x[i]){ cout<<"\n nan's in solution, variable "<<i<<endl; exit(1); return false; }
  }
#endif

  // Get sensible internal energy
  // E = rho*(e + HALF*v.sqr()); 
  //   double e_sens(ZERO);
  //   if (rho()>ZERO) {
  //     for (int i=0; i<ns; i++) y[i] = rhoc(i)/rho();
  //     e_sens = (E()/rho() - HALF*rhovsqr()/rho()/rho()) -  Mixture::heatFormation(y);
  //   }

  // check properties
  if (rho() <= ZERO || !speciesOK(harshness) /*|| e_sens <= ZERO*/) {
    cout << "\n " << CFFC_Name() 
	 << " Flame2D ERROR: Negative Density || Energy || mass fractions: \n" << *this <<endl;
    return false;
  }
  return true;

}


/////////////////////////////////////////////////////////////////////
/// Static Functions
/////////////////////////////////////////////////////////////////////

/****************************************************
 * Setup mixture data
 ****************************************************/
inline void Flame2D_pState::setMixture(const string &mech_name,
				       const string &mech_file,
				       const int soot_model) {

  // deallocate first
  DeallocateStatic();

  // call mixture object setup functin
  Mixture::setMixture(mech_name, mech_file);

  // determine the total number of varibles associated with the gas phase
  ns = Mixture::nSpecies();
  ngas = NUM_FLAME2D_VAR_SANS_SPECIES + ns;

  //
  // if there is soot, set the relevant parameters
  //
  if ( soot_flag = soot_model ) { 

    // get the molar masses, species names, and set the soot object
    double *MWs   = new double[ns];
    string *names = new string[ns];
    MW( MWs );
    speciesNames( names );
    soot.setGasPhase( MWs, ns, names);
    soot.setModelParams( soot_flag );

    // set the relevant Flame2D object params
    iSoot = ngas;
    nsc = soot.NumVar();

  //
  // No soot
  //
  } else {
    iSoot = -1;
    nsc = 0;
  }

  // the total number state variables
  n = ngas + nsc;

  // determine if this is a reacting case
  if (Mixture::nReactions()>0) reacting = true;
  else reacting = false;

  //allocate static memory and load the species data  
  AllocateStatic();

}

/**********************************************************************
 * Flame2D_State::set_gravity -- Set the acceleration due to gravity  *
 *                               in m/s^2.  It acts downwards in the  *
 *                               z-dir (g <= 0)                       *
 **********************************************************************/
inline void Flame2D_State::set_gravity(const double &g) { // [m/s^2]

  // if gravity is acting upwards (wrong way)
  if (g>0) {
    cerr<<"\n Flame2D_pState::set_gravity() - Gravity acting upwards!!!! \n";
    exit(1);
    
  // gravity acting downwards (-ve), OK
  } else {
    gravity_z = g;
  }
}

/****************************************************
 * Print tecplot variable list
 ****************************************************/
inline void Flame2D_State::outputTecplotVarList(ostream &out,
						const string &who, 
						const string &prefix) {

  // Print out U - conserved var list
  if (who == "U") {
    out << "\"" << prefix << "rho\" \\ \n"
	<< "\"" << prefix << "rhou\" \\ \n"
	<< "\"" << prefix << "rhov\" \\ \n"
	<< "\"" << prefix << "rhoE\" \\ \n";
    for(int i =0; i<ns; i++) 
      out <<"\"" << prefix << "rhoc"<<Flame2D_pState::speciesName(i)<<"\" \\ \n";
    for(int i =0; i<nsc; i++) 
      out <<"\"" << prefix << "rhosc" <<i<<"\" \\ \n";

  // or print out W - primitive var list
  } else {
    out << "\"" << prefix << "rho\" \\ \n"
	<< "\"" << prefix << "u\" \\ \n"
	<< "\"" << prefix << "v\" \\ \n"
	<< "\"" << prefix << "p\" \\ \n";
    for(int i =0; i<ns; i++) 
      out <<"\"" << prefix << "c"<<Flame2D_pState::speciesName(i)<<"\" \\ \n";
    for(int i =0; i<nsc; i++) 
      out <<"\"" << prefix << "sc"<<i<<"\" \\ \n";
  }

}


#endif //end _FLAME2D_STATE_INCLUDED 
