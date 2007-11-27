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


//Temperature convergence tolerance in
//Flame2D_cState::T(void)
// these should be moved to CFD.h or Math.h
#define CONV_TOLERANCE  1e-8   //Tolerance used for temperature convergence
#define SPEC_TOLERANCE  1e-8   //Used in negative_speccheck for species round off (was MICRO)
#define TLOWERBOUNDS   200.0   //Uncoment this fix a lower bounds in T calculation

//number of fixed variables in the Flame2D class
#define NUM_FLAME2D_VAR_SANS_SPECIES 4  //rho, v(2), p
enum INDICES { RHO_ = 0,
	       VX_ = 1,
	       VY_ = 2,
	       PRESS_ = 3,
	       SPEC_ = 4 };

// if static species, the total number of variables
#ifdef STATIC_NUMBER_OF_SPECIES
#define STATIC_NUM_FLAME2D_VAR   NUM_FLAME2D_VAR_SANS_SPECIES+STATIC_NUMBER_OF_SPECIES
#endif

//DEBUGGING FLAG FOR FIGUREING OUT proper NS-1 setup.
#define _NS_MINUS_ONE

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

  /**
   * Public Members
   */
public:
  
  /************** Constructors/Destructors ***********/
  Flame2D_State() { Nullify(); Allocate();  }
  Flame2D_State(const double &val) { Nullify(); Allocate(); for (int i=0; i<n; i++) x[i] = val; }
  Flame2D_State(const double &d, const double &vvx, const double &vvy,
		const double &pre, const double &val) { 
    Nullify(); Allocate(); 
    x[RHO_] = d;
    x[VX_] = vvx;
    x[VY_] = vvy;
    x[PRESS_] = pre;
    for (int i=0; i<ns; i++) x[SPEC_+i] = val; 
  }
  Flame2D_State(const double &d, const double &vvx, const double &vvy,
		const double &pre, const double* val) {
    Nullify(); Allocate(); 
    x[RHO_] = d;
    x[VX_] = vvx;
    x[VY_] = vvy;
    x[PRESS_] = pre;
    for (int i=0; i<ns; i++) x[SPEC_+1] = val[i]; 
  }
  Flame2D_State(const Flame2D_State &X) { Nullify(); Allocate(); Copy(X);  }
  ~Flame2D_State() { Deallocate(); }
  void Copy( const Flame2D_State &X ) { for (int i=0; i<n; i++) x[i] = X.x[i]; }

  
  /**************** Allocate / Deallocate ********************/
  void Allocate();
  void Deallocate();
  void Nullify(); 

  /**************** Initializers *****************************/
  void Vacuum() { Zero(); }
  void Zero() { for (int i=0; i<n; i++) x[i] = 0.0; }

  /*************** Accessors ***************************/
  //@{ @name Primitive variables:
  //!< Density.
  double rho(void) const { return x[RHO_]; }
  double& rho(void) { return x[RHO_]; }

  //!< Flow velocity (2D)
  double vx(void) const { return x[VX_]; }
  double& vx(void) { return x[VX_]; }
  double vy(void) const { return x[VY_]; }
  double& vy(void) { return x[VY_]; }
  void v(const Vector2D &V) { x[VX_]=V.x; x[VY_]=V.y; };
  double vsqr(void) const { return x[VX_]*x[VX_] + x[VY_]*x[VY_]; };
  double vabs(void) const { return sqrt(x[VX_]*x[VX_] + x[VY_]*x[VY_]); };
  Vector2D v(void) const { return Vector2D(x[VX_], x[VY_]); };

  //!< Pressure.
  double p(void) const { return x[PRESS_]; }
  double& p(void) { return x[PRESS_]; }

  //!< Species mass fractions
  const double* c(void) const { return &x[SPEC_]; }
  double c(const int&i) const { return x[SPEC_+i]; }
  void c(const double *mfrac) { for (int i=0; i<ns; i++) x[SPEC_+i]=mfrac[i]; }
  void c(const double &val)  { for (int i=0; i<ns; i++) x[SPEC_+i]=val; }
  double& c(const int&i)  { return x[SPEC_+i]; }
  //@}

  //@{ @name Conserved variables:
  //!< Momentum (2D)
  double rhovx(void) const { return x[VX_]; }
  double& rhovx(void) { return x[VX_]; }
  double rhovy(void) const { return x[VY_]; }
  double& rhovy(void) { return x[VY_]; }
  void rhov(const Vector2D &rhoV) { x[VX_]=rhoV.x; x[VY_]=rhoV.y; };
  double rhovsqr(void) const { return (x[VX_]*x[VX_] + x[VY_]*x[VY_]); };
  Vector2D rhov(void) const { return Vector2D(x[VX_], x[VY_]); };

  //!< Total Energy (rho *(e + HALF*v^2))
  double E(void) const { return x[PRESS_]; }
  double& E(void) { return x[PRESS_]; }

  //!< Species mass fractions times density
  double rhoc(const int&i) const { return x[SPEC_+i]; }
  void rhoc(const double *mfrac) { for (int i=0; i<ns; i++) x[SPEC_+i]=mfrac[i]; }
  void rhoc(const double &val)  { for (int i=0; i<ns; i++) x[SPEC_+i]=val; }
  double& rhoc(const int&i)  { return x[SPEC_+i]; }
  //@}


  /*************** Static Functions ********************/
  // return the number of variables - number of species
  static int NumVarSansSpecies() { return NUM_FLAME2D_VAR_SANS_SPECIES; }

   // return the number of variables
  static int NumVar() { return n; }

   // return the number of species
  static int NumSpecies() { return ns; }

  //! set static variables
  static void set_gravity(const double &g);
  static void set_Mref(const double &Ma) { Mref = Ma; };
    
  /******************* Fluxes ***************************/
  void Delta(const Flame2D_State& Ur, const Flame2D_State& Ul) {
    for (int i=0; i<n; i++) x[i] = Ur.x[i] - Ul.x[i];
  };
  void HartenFix_Pos(const Flame2D_State &lambdas_a, const Flame2D_State &lambdas_l,
		     const Flame2D_State &lambdas_r);
  void HartenFix_Neg(const Flame2D_State &lambdas_a, const Flame2D_State &lambdas_l,
		     const Flame2D_State &lambdas_r);
  void HartenFix_Abs(const Flame2D_State &lambdas_a, const Flame2D_State &lambdas_l,
		     const Flame2D_State &lambdas_r);
  void HLLEAverage( const double& wavespeed_r, const double& wavespeed_l, 
		    const Flame2D_State& Fr, const Flame2D_State& Fl, 
		    const Flame2D_State& dUrl );
  void FluxHLLE_x(const Flame2D_pState &Wl, const Flame2D_pState &Wr,
		  const int &Preconditioning);
  void FluxHLLE_n(const Flame2D_pState &Wl, const Flame2D_pState &Wr,
		  const Vector2D &norm_dir, const int &Preconditioning);
  void LindeAverage( const double& wavespeed_r, const double& wavespeed_l, 
		     const double& wavespeed_m, const double& alpha,
		     const Flame2D_State& Fr, const Flame2D_State& Fl, 
		     const Flame2D_State& dUrl );
  void FluxLinde(const Flame2D_pState &Wl, const Flame2D_pState &Wr);
  void FluxLinde_n(const Flame2D_pState &Wl, const Flame2D_pState &Wr,
		   const Vector2D &norm_dir);
  void FluxRoe_x(const Flame2D_pState &Wl,
		 const Flame2D_pState &Wr,
		 const int &Preconditioning,
		 const int &flow_type_flag,
		 const double &deltax);
  void FluxRoe_n(const Flame2D_pState &Wl,
		 const Flame2D_pState &Wr,
		 const Vector2D &norm_dir,
		 const int &Preconditioning,
		 const int &flow_type_flag,
		 const double &delta_n );

  /**************** Operators Overloading ********************/
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


  /**
   * Private Objects
   */
protected:

#ifdef STATIC_NUMBER_OF_SPECIES
  static const int  n = STATIC_NUM_FLAME2D_VAR;
  static const int ns = STATIC_NUMBER_OF_SPECIES;
  double            x[n];
#else
  static int        n; //!< Total number of vars
  static int       ns; //!< Number of species
  double           *x;
#endif
#ifdef _NS_MINUS_ONE
  static const int  NSm1 = 1;
#else
  static const int  NSm1 = 0;
#endif
  static double        Mref; //!< Mref for Precondtioning (normally set to incoming freestream Mach)
  static double   gravity_z; //!< m/s^2 acceleration due to gravity  

};

/////////////////////////////////////////////////////////////////////
/// Member Functinos
/////////////////////////////////////////////////////////////////////


/****************************************************
 * Allocator/deallocator for object memory.
 ****************************************************/
inline void Flame2D_State :: Allocate() {
#ifndef STATIC_NUMBER_OF_SPECIES
  Deallocate();
  if (n>0) { x = new double[n]; }
#endif
  };

inline void Flame2D_State :: Deallocate() { 
#ifndef STATIC_NUMBER_OF_SPECIES
  if (x!=NULL) { delete[] x; x = NULL; } 
#endif
};

inline void Flame2D_State :: Nullify() { 
#ifndef STATIC_NUMBER_OF_SPECIES
  x = NULL;
#endif
};

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


  /************** Private Objects **********************/
private:
  Mixture Mix;       //!< Mixture Object
  static double* r;  //!< Temporary storage for reaction rates
  static double* dihdic;  //!< Temporary storage for dihdic

public:
  /************** Constructors/Destructors ***********/
  Flame2D_pState() {
    rho()=DENSITY_STDATM; p()=PRESSURE_STDATM;
    vx()=0.0; vy()=0.0; c(1.0/ns);
  }

  Flame2D_pState(const double &d, const Vector2D &V, 
		 const double &pre)
  { rho()=d; p()=pre; v(V); c(1.0/ns); setGas(); }

  Flame2D_pState(const double &d, const double &vvx, 
		const double &vvy, const double &pre)
  { rho()=d; p()=pre; vx()=vvx; vy()=vvy; c(1.0/ns); setGas(); }
  
  Flame2D_pState(const double &d, const double &vvx, 
		const double &vvy, const double &pre, 
		const double *mfrac)
  { rho()=d; p()=pre; vx()=vvx; vy()=vvy; c(mfrac); setGas(); }
  
  Flame2D_pState(const double &d, const Vector2D &V, 
		const double &pre, const double *mfrac)
  { rho()=d; p()=pre; v(V); c(mfrac); setGas(); }

  Flame2D_pState(const Flame2D_pState &W) { Copy(W); }

  Flame2D_pState(const Flame2D_State &W) { Copy(W);  setGas(); }

  ~Flame2D_pState() { }

  void Copy( const Flame2D_pState &W ) { 
    for (int i=0; i<n; i++) x[i] = W.x[i];
    Mix = W.Mix;
  }

  /*************** VACUUM OPERATOR *********************/

  /*************** Static Functions ********************/
  //! initial mixture setup function
  static void setMixture(const string &mech_name,
			 const string &mech_file);

  //! Static memory allocator/deallocator  
  static void AllocateStatic(void);
  static void DeallocateStatic(void);
  //! The species names
  static string speciesName(const int&i) { return Mixture::speciesName(i); };

  /*************** Accessors ***************************/
  //!< Density.
  double rho(void) const { return x[RHO_]; }

  //!< Flow velocity (2D)
  double vx(void) const { return x[VX_]; }
  double vy(void) const { return x[VY_]; }
  double vsqr(void) const { return x[VX_]*x[VX_] + x[VY_]*x[VY_]; };
  double vabs(void) const { return sqrt(x[VX_]*x[VX_] + x[VY_]*x[VY_]); };
  Vector2D v(void) const { return Vector2D(x[VX_], x[VY_]); };

  //!< Pressure.
  double p(void) const { return x[PRESS_]; }

  //!< Species mass fractions
  const double* c(void) const { return &x[SPEC_]; }
  double c(const int&i) const { return x[SPEC_+i]; }

  //!< Momentum (2D)
  double rhovx(void) const { return rho()*vx(); };
  double rhovy(void) const { return rho()*vy(); };
  double rhovsqr(void) const { return rhovx()*rhovx()+rhovy()*rhovy(); }
  Vector2D rhov(void) const { return Vector2D(rhovx(), rhovy()); };

  //!< Species mass fractions times density
  double rhoc(const int&i) const { return rho()*c(i); };

  /*************** Assignment Objects ******************/
private:
  //!< Density.
  double& rho(void) { return x[RHO_]; }
  double& vx(void) { return x[VX_]; }

  //!< Flow velocity (2D)
  double& vy(void) { return x[VY_]; }
  void v(const Vector2D &V) { x[VX_]=V.x; x[VY_]=V.y; };

  //!< Pressure.
  double& p(void) { return x[PRESS_]; }

  //!< Species mass fractions
  void c(const double *mfrac) { for (int i=0; i<ns; i++) x[SPEC_+i]=mfrac[i]; }
  void c(const double &val)  { for (int i=0; i<ns; i++) x[SPEC_+i]=val; }
  double& c(const int&i)  { return x[SPEC_+i]; }

  //@{ @name Conserved variables (override these so they don't work):
  //!< Momentum (2D)
  double& rhovx(void) { assert(1); }
  double& rhovy(void) { assert(1); }

  //!< Total Energy (rho *(e + HALF*v^2))
  double& E(void) { assert(1); }

  //!< Species mass fractions times density
  void rhoc(const double *mfrac) { assert(1); }
  void rhoc(const double &val)  { assert(1); }
  double& rhoc(const int&i)  { assert(1); }
  //@}

  /*************** Assignment Objects ******************/
  //!< Set gas state 
public:
  void setState_TPY(const double &Temp, const double &Press, const double *y) {    
    p() = Press;
    c( y );
    Mix.setState_TPY(Temp, Press, y); 
    rho() = p() / (Mix.gasConstant()*Temp);
  };
  void setVelocity(const double &vvx, const double &vvy) {    
    vx() = vvx;
    vy() = vvy;
  };
  double load_dihdic(void) const {
    Mix.getDihdDc( p(), c(), dihdic );
#ifdef _NS_MINUS_ONE
    for (int i=0; i<ns; i++) dihdic[i] -= dihdic[ns-1];
#endif
  };

private:
  void setGas(void) { Mix.setState_DPY(rho(), p(), c()); };
  void setGas_T(const double&Temp) { 
    Mix.setState_TPY(Temp, p(), c());
    rho() = p() / (Mix.gasConstant()*Temp);
  };
  void setGas_E(const double &en) { 
    Mix.setState_DEY(rho(), en, c()); 
    p() = rho()*Mix.gasConstant()*Mix.temperature();
  };
  void setGas_H(const double &h) { 
    Mix.setState_DHY(rho(), h, c()); 
    p() = rho()*Mix.gasConstant()*Mix.temperature();
  };

public:
 /********* Primitive / Conserved Transformation ******/
  //!< conserved to primitive transforation
  void setU(const Flame2D_State &U) {
    rho() = U.rho();
    vx() = U.rhovx()/U.rho();
    vy() = U.rhovy()/U.rho();
    for(int i=0; i<ns; i++) c(i) = U.rhoc(i)/U.rho();
    double e( U.E()/U.rho() - 0.5*U.rhovsqr()/(U.rho()*U.rho()) );
    setGas_E(e);
  }
  void setW(const Flame2D_pState &W){ if( this != &W)  Copy(W); }
  void setW(const Flame2D_State &W){
    rho() = W.rho();
    vx() = W.vx();
    vy() = W.vy();
    p() = W.p();
    for(int i=0; i<ns; i++) c(i) = W.c(i);
    setGas();
  }
  void getU(Flame2D_State &U) const {
    U.rho() = rho();
    U.rhovx() = rhovx();
    U.rhovy() = rhovy();
    U.E() = E();
    for(int i=0; i<ns; i++) U.rhoc(i) = rhoc(i);
  }
  Flame2D_State U(void) const {
    Flame2D_State U;
    getU(U);
    return U;
  }


  /************ Misc Quantities of Interest *************/
  //!< Total Energy (rho *(e + HALF*v^2))
  double E(void) const { return (rho()*(e() + 0.5*vsqr())); };

  //!< Total Enthalpy ( H = h + velocity^2 )
  double H(void) const { return (rho()*(h() + 0.5*vsqr())); };
  double Hs(void) const { return (rho()*(hs() + 0.5*vsqr())); };

  //!< Speed of sound 
  double a(void) const { return sqrt(g()*Rtot()*T()); }

  /************ Mixture Object Wrappers *****************/
  //!< Temperature
  double T(void) const { return Mix.temperature(); };
  //!< Gas constant
  double Rtot(void) const { return Mix.gasConstant(); };
  //!< Internal Energy
  double e(void) const { return Mix.internalEnergy(); };
  double es(void) const { return Mix.internalEnergySens(); };
  //!< enthalpy
  double h(void) const { return Mix.enthalpy(); };
  double hs(void) const { return Mix.enthalpySens(); };
  double hprime(void) const { return Mix.enthalpyPrime(); };
  //!< heat capacity
  double Cp(void) const { return Mix.heatCapacity_p(); };
  double Cv(void) const { return Mix.heatCapacity_v(); };
  double g(void) const { return Mix.heatRatio(); };
  //!< transport
  double mu(void) const { return Mix.gasConstant(); };
  double kappa(void) const { return Mix.gasConstant(); };
  double Diffusion_coef(const int &i) const { return Mix.gasConstant(); };
  //!< Dimensionaless
  double Sc(const int &i) const { return Mix.schmidt(rho(), i); };
  double Pr(void) const { return Mix.prandtl(); };
  double Le(const int &i) const { return Mix.lewis(rho(), i); };
  double Re(const double &l) const { return rho()*vabs()*l/mu(); };
  //!< Derivatives
  double diedip() const { return (hprime() - Rtot())/(rho()*Rtot()); };
  double diedirho() const { return -p()*(hprime() - Rtot())/(rho()*rho()*Rtot()); };

  /************* Eigenvalues / Eigenvectors *************/
  void lambda_x(Flame2D_State &lambdas) const;
  Flame2D_State rc_x(const int &index) const;
  Flame2D_State lp_x(const int &index) const;
  void lambda_preconditioned_x(Flame2D_State &lambdas, const double &MR2) const;
  Flame2D_State rc_x_precon(const int &index, const double &MR2) const;
  Flame2D_State lp_x_precon(const int &index, const double &MR2) const;

  /******************* Fluxes ***************************/
  void Fx(Flame2D_State &FluxX) const;
  Flame2D_State Fx(void) const;
  void addFx(Flame2D_State &FluxX, const double& mult) const;
  void RoeAverage(const Flame2D_pState &Wl, const Flame2D_pState &Wr);

  /********** Axisymmetric Source Terms *****************/
  void Sa_inviscid(Flame2D_State &S, const Vector2D &X, 
		   const int Axisymmetric, const double& mult=1.0) const;

  /* Source terms associated with finite-rate chemistry */
  void Sw(Flame2D_State &S, const double& mult=1.0) const;

  /** Source terms associated with gravitational forces */
  void Sg(Flame2D_State &S, const double& mult=1.0) const;

  /*************** Preconditioning **********************/
  double u_plus_aprecon(const double &u,const int &flow_type_flag,
			const double &deltax) const;
  void u_a_precon(const double &UR,double &uprimed, double &cprimed) const;
  double Mr2(const int &flow_type_flag, const double &deltax) const;
  void Low_Mach_Number_Preconditioner(::DenseMatrix &P,
				      const int &Viscous_flag, 
				      const double &deltax ) const;
  void Low_Mach_Number_Preconditioner_Inverse(::DenseMatrix &Pinv,	
					      const int &Viscous_flag, 
					      const double &deltax ) const;
  
  /************** Boundary Conditions *******************/
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
    
  /************* Operators Overloading ******************/
  // Assignment Operator.
  Flame2D_pState& operator =(const Flame2D_pState &W); 
 
  // Input-output operators.
  friend ostream& operator << (ostream &out_file, const Flame2D_pState &X);
  friend istream& operator >> (istream &in_file,  Flame2D_pState &X);



};



/////////////////////////////////////////////////////////////////////
/// LOCAL MEMBER FUNCTIONS
/////////////////////////////////////////////////////////////////////

/****************************************************
 * Allocator/deallocator for static memory.
 ****************************************************/
inline void Flame2D_pState :: AllocateStatic() {
  if (ns>0) { 
    r = new double[ns];
    dihdic = new double[ns];
  }
};

inline void Flame2D_pState :: DeallocateStatic() { 
  if (r!=NULL) { delete[] r; r = NULL; } 
  if (dihdic!=NULL) { delete[] dihdic; dihdic = NULL; } 
};


/********************************************************
 * Flame2D_State -- Binary arithmetic operators.        *
 ********************************************************/
//----------------- Assignment ----------------------------//
inline Flame2D_pState& Flame2D_pState::operator =(const Flame2D_pState &W){
  if( this != &W) Copy(W);
  return (*this);
}

/********************************************************
 * Flame2D_pState -- Input-output operators.            *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const Flame2D_pState &W) {
  out_file.precision(10);
  out_file.setf(ios::scientific);
  for( int i=0; i<W.NumVar(); i++) out_file<<" "<<W[i+1];
  out_file << endl << W.Mix;
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


/*********************************************************************************
 *********************************************************************************
 *********************************************************************************
 *********************** EXTERNAL FUNCTIONS **************************************
 *********************************************************************************
 *********************************************************************************
 *********************************************************************************/



#endif //end _FLAME2D_STATE_INCLUDED 
