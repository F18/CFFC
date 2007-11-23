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

// define early so can be used in other classes
class Flame2D_cState;
class Flame2D_pState;

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

// If you define this variable, the number of species will be
// predetermined for faster calculations.., however it is not as general 
#define STATIC_NUMBER_OF_SPECIES 36 //2 AIR, 6 2STEP_CH4
#ifdef STATIC_NUMBER_OF_SPECIES
#define STATIC_NUM_FLAME2D_VAR   NUM_FLAME2D_VAR_SANS_SPECIES+STATIC_NUMBER_OF_SPECIES
#endif

//DEBUGGING FLAG FOR FIGUREING OUT proper NS-1 setup.
#define _NS_MINUS_ONE

/////////////////////////////////////////////////////////////////////
/// Class Definitions
/////////////////////////////////////////////////////////////////////
class Flame2D_State;
class Flame2D_pState;
class Flame2D_cState;


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
  Flame2D_State() { Allocate();  }
  Flame2D_State(const double &val) { Allocate(); for (int i=0; i<n; i++) x[i] = val; }
  Flame2D_State(const Flame2D_State &X) { Allocate(); Copy(X);  }
  ~Flame2D_State() { Deallocate(); }
  void Copy( const Flame2D_State &X ) { for (int i=0; i<n; i++) x[i] = X.x[i]; }


  /**************** Allocate / Deallocate ********************/
  void  Allocate();
  void  Deallocate();

  /**************** Initializers *****************************/
  void Vacuum() { Zero(); }
  void Zero() { for (int i=0; i<n; i++) x[i] = 0.0; }

  /*************** Static Functions ********************/
  // return the number of variables - number of species
  static int NumVarSansSpecies() { return NUM_FLAME2D_VAR_SANS_SPECIES; }

   // return the number of variables
  static int NumVar() { return n; }

   // return the number of species
  static int NumSpecies() { return ns; }

  //! set acceleration due to gravity
  static void set_gravity(const double &g);

  /**************** Operators Overloading ********************/
  // Index operator
  double &operator[](int index);
  const double &operator[](int index) const;

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
  static int       ns; //!< Number of 
  double           *x;
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

/********************************************************
 * Flame2D_State -- Index operators.                    *
 ********************************************************/
inline double& Flame2D_State::operator[](int index) {  
  return x[index-1];
}

inline const double& Flame2D_State::operator[](int index) const {  
  return x[index-1];
}

/********************************************************
 * Flame2D_State -- Binary arithmetic operators.        *
 ********************************************************/
//----------------- Addition -----------------------------//
Flame2D_State Flame2D_State::operator +(const Flame2D_State &X) const{    
  Flame2D_State Temp(*this);
  Temp += X;
  return Temp;
}

//------------------ Subtraction ------------------------//
Flame2D_State Flame2D_State::operator -(const Flame2D_State &X) const{
    Flame2D_State Temp(*this);
    Temp -= X;
    return Temp;
}

//---------------- Scalar Multiplication ------------------//
Flame2D_State Flame2D_State::operator *(const double &a) const{
  Flame2D_State Temp(*this);
  for( int i=0; i<n; i++) Temp.x[i] = x[i]*a; 
  return(Temp);
}

Flame2D_State operator *(const double &a, const Flame2D_State &X){
  Flame2D_State Temp;
  for( int i=0; i<X.n; i++) Temp.x[i] = X.x[i]*a;
  return(Temp);
}

//--------------- Scalar Division ------------------------//
Flame2D_State Flame2D_State::operator /(const double &a) const {
  Flame2D_State Temp(*this);
  for(int i=0; i<n; i++) Temp.x[i] = x[i]/a;
  return(Temp);
}

//----------------- Inner Product ------------------------//
double Flame2D_State::operator *(const Flame2D_State &X) const{
  double sum(0.0);
  for(int i=0; i<n; i++) sum += x[i]*X.x[i];
  return (sum);
}

//----------- solution state product operator ------------//
Flame2D_State Flame2D_State::operator ^( const Flame2D_State &X) const {
    Flame2D_State Temp(*this);
    for(int i=0; i<n; i++) Temp.x[i] = x[i]*X.x[i];
    return(Temp);
}

//----------------- Assignment ----------------------------//
Flame2D_State& Flame2D_State::operator =(const Flame2D_State &X){
  //self assignment protection
  if( this != &X){ 
    Copy(X);
  }
  return (*this);
}

/********************************************************
 * Flame2D_pState -- Shortcut arithmetic operators.     *
 ********************************************************/
Flame2D_State& Flame2D_State::operator +=(const Flame2D_State &X){
  for( int i=0; i<n; i++)  x[i] += X.x[i];
  return (*this);
}

Flame2D_State& Flame2D_State::operator -=(const Flame2D_State &X) {
  for(int i=0; i<n; i++) x[i] -= X.x[i];
  return (*this); 
}

/********************************************************
 * Flame2D_pState -- Unary arithmetic operators.         *
 ********************************************************/
Flame2D_State operator -(const Flame2D_State &X) {
  Flame2D_State Temp(X);
  for(int i=0; i<X.n; i++)  Temp.x[i] = -X.x[i]; 
  return(Temp);
}



/********************************************************
 * Flame2D_pState -- Input-output operators.            *
 ********************************************************/
ostream &operator << (ostream &out_file, const Flame2D_State &X) {
  out_file.precision(10);
  out_file.setf(ios::scientific);
  for( int i=0; i<X.n; i++) out_file<<" "<<X.x[i];
  out_file.unsetf(ios::scientific);
  return (out_file);
}

istream &operator >> (istream &in_file, Flame2D_State &X) {
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


  /**
   * Public Objects
   */
public:


  /**
   * Private Objects
   */
private:
  double T;   //!< Temperature
  Mixture M;  //!< Mixture Object

  /**
   * Public Members
   */
public:

  /************** Constructors/Destructors ***********/
  Flame2D_pState() : T(0.0) {
    rho(DENSITY_STDATM); p(PRESSURE_STDATM);
    vx(0.0); vy(0.0); c(1/ns); T = T();
  }

  Flame2D_pState(const double &d, const Vector2D &V, 
		 const double &pre)
  { rho(d); p(pre); v(V); c(1/ns); T = T(); }

  Flame2D_pState(const double &d, const double &vvx, 
		const double &vvy, const double &pre)
  { rho(d); p(pre); vx(vvx); vy(vvy); c(1/ns); T = T(); }
  
  Flame2D_pState(const double &d, const double &vvx, 
		const double &vvy, const double &pre, 
		const double *mfrac) : T(0.0)
  { rho(d); p(pre); vx(vvx); vy(vvy); c(mfrac); T = T(); }
  
  Flame2D_pState(const double &d, const Vector2D &V, 
		const double &pre, const double *mfrac) : T(0.0)
  { rho(d); p(pre); v(V); c(mfrac); T = T(); }

  Flame2D_pState(const Flame2D_pState &W) { Copy(W); }
  Flame2D_pState(const Flame2D_cState &U) { W(U); }
  Flame2D_pState(const Flame2D_State &X) { Copy(X); }
  ~Flame2D_pState() { }

  /*************** Accessors ***************************/
  //@{ @name Primitive variables:
  //!< Density.
  double rho(void) const { return x[RHO_]; }
  void rho(const double &d) { x[RHO_]=d; }

  //!< Flow velocity (2D)
  double vx(void) const { return x[VX_]; }
  void vx(const double &vvx) { x[VX_]=vvx; }
  double vy(void) const { return x[VY_]; }
  void vy(const double &vvy) { x[VY_]=vvy; }
  void v(const Vector2D &V) { x[VX_]=V.x; x[VY_]=V.y; };
  double vsqr(void) const { return x[VX_]*x[VX_] + x[VY_]*x[VY_]; };

  //!< Pressure.
  double p(void) const { return x[PRESS_]; }
  void p(const double &pr) { x[PRESS_]=pr; }

  //!< Species mass fractions
  double c(const int&i) const { return x[SPEC_+i]; }
  void c(const double *mfrac) { for (int i=0; i<ns; i++) x[SPEC_+i]=mfrac[i]; }
  void c(const double &val)  { for (int i=0; i<ns; i++) x[SPEC_+i]=val; }
  void c(const int&i, const double &val)  { x[SPEC_+i]=val; }
  //@}

  /********* Primitive / Conserved Transformation ******/
  //!< conserved to primitive transforation
  void W(const Flame2D_cState &U);

  //!< Momentum (2D)
  double rhovx(void) const { return rho()*vx(); };
  double rhovy(void) const { return rho()*vy(); };

  //!< Species mass fractions times density
  double rhoc(const int&i) const { return rho()*c(i); };

  //!< Total Energy (rho *(e + HALF*v^2))
  double E(void) const { return (rho()*(e() + 0.5*vsqr())); };

  //!< Total Enthalpy ( H = h + velocity^2 )
  double H(void) const { return (rho()*(h() + 0.5*vsqr())); };
  double Hs(void) const { return (rho()*(hs() + 0.5*vsqr())); };

  //!< Speed of sound 
  double a(void) const { return sqrt(g()*Rtot()*T); }


  /************ Mixture Object Wrappers *****************/

  //!< Gas constant
  double Rtot(void) const { return M.gasConstant(); };
  //!< Internal Energy
  double e(void) const { return M.internalEnergy(T); };
  double es(void) const { return M.internalEnergySens(T); };
  //!< enthalpy
  double h(void) const { return M.enthalpy(); };
  double hs(void) const { return M.enthalpySens(); };
  //!< heat capacity
  double Cp(void) const { return M.heatCapacity_p(); };
  double Cv(void) const { return M.heatCapacity_v(); };
  double g(void) const { return M.heatRatio(); };
  //!< transport
  double mu(void) const { return M.gasConstant(); };
  double kappa(void) const { return M.gasConstant(); };
  double Diffusion_coef(const int &i) const { return M.gasConstant(); };
  //!< Dimensionaless
  double Sc(const int &i) const { return M.schmidt(rho(), i); };
  double Pr(void) const { return M.prandtl(); };
  double Le(const int &i) const { return M.lewis(rho(), i); };



  /*************** VACUUM OPERATOR *********************/


  /*************** Static Functions ********************/


};


/*!
 * Class: Flame2D_cState
 *
 * @brief  Conserved variable solution state class definition for an
 *         inviscid, laminar, or turbulent chemically reacting gas-flow.
 *
 * Conserved variable solution state class definition for an inviscid, 
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
class Flame2D_cState : public Flame2D_State {

  /**
   * Public Objects
   */
public:


  /**
   * Private Objects
   */
private:


  /**
   * Public Members
   */
public:

  /************** Constructors/Destructors ***********/
  Flame2D_cState() {
    rho(DENSITY_STDATM); E(PRESSURE_STDATM/(rho()*0.4));
    rhovx(0.0); rhovy(0.0); rhoc(0.0);
  }

  Flame2D_cState(const double &value)
  { rho(value); E(value); rhovx(value); rhovy(value); rhoc(0.0); }
  
  Flame2D_cState(const double &d, const Vector2D &rhoV, 
		 const double &en) 
  { rho(d); E(en); rhov(rhoV); rhoc(0.0); }
  
  Flame2D_cState(const double &d, const double &rhovvx, 
		 const double &rhovvy, const double &en) 
  { rho(d); E(en); rhovx(rhovvx); rhovy(rhovvy); rhoc(0.0); }

  Flame2D_cState(const double &d, const Vector2D &rhoV, 
		 const double &en, const double &value) 
  { rho(d); E(en); rhov(rhoV); rhoc(value); }

  Flame2D_cState(const double &d, const double &rhovvx, 
		 const double &rhovvy, const double &en, 
		 const double &value)
  { rho(d); E(en); rhovx(rhovvx); rhovy(rhovvy); rhoc(value); }
  
  Flame2D_cState(const double &d, const double &rhovvx, 
		 const double &rhovvy, const double &en, 
		 const double *mfrac) 
  { rho(d); E(en); rhovx(rhovvx); rhovy(rhovvy); rhoc(mfrac); }
  
  Flame2D_cState(const double &d, const Vector2D &rhoV, 
		 const double &en, const double *mfrac) 
  { rho(d); E(en); rhov(rhoV); rhoc(mfrac); }

  Flame2D_cState(const Flame2D_cState &U) { Copy(U); }
  Flame2D_cState(const Flame2D_pState &W) { U(W); }
  Flame2D_cState(const Flame2D_State &X) { Copy(X); }
  ~Flame2D_cState() { }

  /*************** Accessors ***************************/
  //@{ @name Conserved variables:
  //!< Density.
  const double& rho(void) const { return x[RHO_]; }
  void rho(const double &d) { x[RHO_]=d; }

  //!< Momentum (2D)
  double rhovx(void) const { return x[VX_]; }
  void rhovx(const double &rhovvx) { x[VX_]=rhovvx; }
  double rhovy(void) const { return x[VY_]; }
  void rhovy(const double &rhovvy) { x[VY_]=rhovvy; }
  void rhov(const Vector2D &rhoV) { x[VX_]=rhoV.x; x[VY_]=rhoV.y; };
  double vsqr(void) const { return (x[VX_]*x[VX_] + x[VY_]*x[VY_])/(x[RHO_]*x[RHO_]); };

  //!< Total Energy (rho *(e + HALF*v^2))
  double E(void) const { return x[PRESS_]; }
  void E(const double &en) { x[PRESS_]=en; }

  //!< Species mass fractions times density
  double rhoc(const int&i) const { return x[SPEC_+i]; }
  void rhoc(const double *mfrac) { for (int i=0; i<ns; i++) x[SPEC_+i]=mfrac[i]; }
  void rhoc(const double &val)  { for (int i=0; i<ns; i++) x[SPEC_+i]=val; }
  void rhoc(const int&i, const double &val)  { x[SPEC_+i]=val; }
  //@}

  /********* Primitive / Conserved Transformation ******/
  //!< primitive to conserved transforation
  void U(const Flame2D_pState &W);

  //!< Flow velocity (2D)
  double vx(void) const { return x[VX_]/x[RHO_]; }
  double vy(void) const { return x[VY_]/x[RHO_]; }

  //!< Pressure.
  double p(void) const { return (rho()*Rtot()*T()); };

  //!< Species mass fractions
  double c(const int&i) const { return x[SPEC_+i]/x[RHO_]; }

  //!< Speed of sound 
  double a(void) const { return sqrt(g()*Rtot()*T()); };

  //!< Temperature 
  double T(void) const;

  /************ Mixture Object Wrappers *****************/

  //!< Gas constant
  double Rtot(void) const;
  //!< Internal Energy
  double e(void) const;
  double es(void) const;
  //!< enthalpy
  double h(void) const;
  double hs(void) const;
  //!< heat capacity
  double Cp(void) const;
  double Cv(void) const;
  double g(void) const;
  //!< transport
  double mu(void) const;
  double kappa(void) const;
  double Diffusion_coef(const int &i) const;



  /*************** VACUUM OPERATOR *********************/


  /*************** Static Functions ********************/


 };

/////////////////////////////////////////////////////////////////////
/// LOCAL MEMBER FUNCTIONS
/////////////////////////////////////////////////////////////////////

/****************************************************
 * Primitive / conserved transformations.
 ****************************************************/

inline void Flame2D_cState::U(const Flame2D_pState &W) {
  rho(W.rho());
  rhovx(W.rhovx());
  rhovy(W.rhovy());
  E(W.E());
  for(int i=0; i<W.NumSpecies(); i++) rhoc(i, W.rhoc(i));
}

inline void Flame2D_pState::W(const Flame2D_cState &U) {
  rho(U.rho());
  vx(U.vx());
  p(U.p());
  for(int i=0; i<U.NumSpecies(); i++) c(i, U.c(i));
}



#endif //end _FLAME2D_STATE_INCLUDED 
