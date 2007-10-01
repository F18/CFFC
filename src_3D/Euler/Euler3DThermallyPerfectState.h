/* Euler3DThermallyPerfectState.h:  Header file defining the Euler solution state
                                    classes associated with the solution of 
                                    a thermally perfect non-reactive or combusting
                                    mixture. */

#ifndef _EULER3D_THERMALLYPERFECT_STATE_INCLUDED 
#define _EULER3D_THERMALLYPERFECT_STATE_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

/* Include required CFFC header files. */

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _MATRIX_INCLUDED
#include "../Math/Matrix.h"
#endif // _MATRIX_INCLUDED

#ifndef _TENSOR3D_INCLUDED
#include "../Math/Tensor3D.h"
#endif //_TENSOR3D_INCLUDED

#ifndef _VECTOR3D_INCLUDED
#include "../Math/Vector3D.h"
#endif //_VECTOR3D_INCLUDED

#ifndef _GAS_CONSTANTS_INCLUDED
#include "../Physics/GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED

#ifndef _SPECIES_INCLUDED
#include "../Physics/Species.h"
#endif //_SPECIES_INCLUDED

#ifndef _NASARP1311_DATA_INCLUDED
#include "../Physics/NASAData/NASARP1311data.h"
#endif // _NASARP1311_DATA_INCLUDED

#ifndef _REACTIONS_INCLUDED
#include "../Reactions/Reactions.h"
#endif // _REACTIONS_INCLUDED

/* Define the class. */

class Euler3D_ThermallyPerfect_cState;
class Euler3D_ThermallyPerfect_pState;

//number of fixed variables in the Euler3D_ThermallyPerfect class
#define NUM_EULER3D_VAR_SANS_SPECIES 5  //rho, v(3), p

#define CONV_TOLERANCE 1e-8

/*!
 * Class: Euler3D_ThermallyPerfect_pState
 *
 * @brief Euler 3D primitive solution state class.
 *
 * This class defines the primitive solution state for the
 * Euler equations governing three-dimensional flows of
 * thermally perfect non-reactive and combusting mixtures.
 *
 * \verbatim
 * Member functions
 *      d       -- Return density
 *      v       -- Return flow velocity
 *      p       -- Return pressure
 *      g       -- Return specific heat ratio
 *      gm1     -- Return g-1
 *      gm1i    -- Return 1/(g-1)
 *      R       -- Return gas constant
 *      T       -- Return temperature
 *      e       -- Return total energy
 *      E       -- Return interal energy
 *      h       -- Return specific enthalpy
 *      H       -- Return total enthalpy
 *      a       -- Return sound speed
 *      a2      -- Return sound speed square
 *      M       -- Return Mach number
 *      s       -- Return specific entropy
 *      dv      -- Return momentum
 *      To      -- Return stagnation temperature
 *      po      -- Return stagnation pressure
 *      ao      -- Return stagnation sound speed
 *      ho      -- Return stagnation enthalpy
 *      U       -- Return conserved solution state
 *      F       -- Return x-direction solution flux
 *      G       -- Return y-direction solution flux
 *      H       -- Return z-direction solution flux
 *      Fn      -- Return n-direction solution flux
 *      lambda  -- Return eigenvalue
 *      rp      -- Return primitive right eigenvector
 *      rc      -- Return conserved right eigenvector
 *      lp      -- Return primitive left eigenvector
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
 * W = +W;
 * W = -W;
 * W += W;
 * W -= W;
 * W == W;
 * W != W;
 * cout << W; (output function)
 * cin  >> W; (input function)
 *
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
 * \endverbatim
 */
class Euler3D_ThermallyPerfect_pState {
  public:
   double         rho;  
   Vector3D         v;   
   double           p;   
   Species      *spec;   
   
   static int num_vars; 
   static int ns;               
   static NASARP1311data *specdata;  
   static double *Schmidt;
   static Reaction_set React;       
   static double low_temp_range;    
   static double high_temp_range;   
   static int debug_level;        
   
  //@{ @name Constructors and desctructors:
  //! Default constructor (assign default values).
   Euler3D_ThermallyPerfect_pState(): 
    rho(DENSITY_STDATM), p(PRESSURE_STDATM), spec(NULL) {
      v.zero(); set_initial_values();
   }
   
   Euler3D_ThermallyPerfect_pState(const double &value): 
    rho(value), p(value), spec(NULL) {
      v.x = value;  v.y = value; v.z = value; 
      set_initial_values(value);
   }
         
   Euler3D_ThermallyPerfect_pState(const double &d, 
                                   const Vector3D &V,
                                   const double &pre):
    rho(d), v(V), p(pre), spec(NULL) {
      set_initial_values();
   }

   Euler3D_ThermallyPerfect_pState(const double &d, 
                                   const Vector3D &V, 
                                   const double &pre, 
                                   const double &frac):
    rho(d), v(V), p(pre), spec(NULL) {
      set_initial_values(frac);
   }
             
   Euler3D_ThermallyPerfect_pState(const double &d, 
                                   const double &vx, 
                                   const double &vy, 
                                   const double &vz, 
                                   const double &pre, 
                                   const double &frac):
    rho(d), p(pre), spec(NULL) {
      v.x = vx; v.y = vy; v.z = vz; 
      set_initial_values(frac); 
   }
                 
   Euler3D_ThermallyPerfect_pState(const double &d, 
                                   const double &vx,
                                   const double &vy, 
                                   const double &vz,
                                   const double &pre):
    rho(d), p(pre), spec(NULL) {
      v.x = vx;  v.y = vy; v.z = vz;
      set_initial_values(); 
   }
   
   Euler3D_ThermallyPerfect_pState(const double &d, 
                                   const double &vx, 
                                   const double &vy,  
                                   const double &vz,
                                   const double &pre, 
                                   Species *mfrac):
    rho(d), p (pre), spec(NULL) {
      v.x=vx; v.y=vy; v.z = vz; 
      set_initial_values(mfrac); 
   }
  
   Euler3D_ThermallyPerfect_pState(const double &d, 
                                   const Vector3D &V, 
                                   const double &pre, 
                                   Species *mfrac):
    rho(d), v(V), p(pre), spec(NULL) {
      set_initial_values(mfrac); 
   }
  
   Euler3D_ThermallyPerfect_pState(const double &d, 
                                   const Vector3D &V,
	 	                   const double &pre, 
                                   double *mfrac):
    rho(d), v(V), p(pre), spec(NULL) {
      set_initial_values(mfrac); 
   }

   // This is needed for the operator overload returns!!!!
   Euler3D_ThermallyPerfect_pState(const Euler3D_ThermallyPerfect_pState &W) {
      spec = NULL; rho = DENSITY_STDATM; 
      set_initial_values(); Copy(W);
   }
     
   //! Default destructor
   ~Euler3D_ThermallyPerfect_pState() {
      Deallocate(); 
   }
   //@}

   //@{ @name Other Member functions:
   void Deallocate(void){
      if(spec != NULL) delete[] spec;
      spec = NULL;
   }

   void Deallocate_static(void){
      if (specdata != NULL) delete[] specdata;
      specdata = NULL;
      if (Schmidt != NULL) delete[] Schmidt; 
      Schmidt = NULL;
   }

   // Sets the species data, needs to be called only once as it is static
   void set_species_data(const int &n,
                         const string *S,
                         const char *PATH,
                         const int &debug, 
                         const double &Mr, 
                         const double* Sc,
                         const int &trans_data);
     
   //Set initial data values
   void set_initial_values();
   void set_initial_values(const double &value);
   void set_initial_values(double *cfrac);
   void set_initial_values(Species *mfrac);

   //Copy construtor, cheaper than = operator
   void Copy(const Euler3D_ThermallyPerfect_pState &W);

   // return the number of variables - number of species
   int NumVarSansSpecies() const { return num_vars - ns; }

   // Vacuum operator
   void Vacuum(){ 
     rho=ZERO; v.x=ZERO; v.y=ZERO; v.z=ZERO;  p = ZERO;
       
     for(int i=0; i<ns; i++){
        spec[i].Vacuum();
     }
   }
  
   void zero_non_sol(){
     for(int i=0; i<ns; i++){
           spec[i].gradc.zero();
           spec[i].diffusion_coef=ZERO;
        }
   }  

   // Set Data Temperature Ranges 
   void Temp_low_range();     
   void Temp_high_range(); 

   /***************** Mixing Rules ************************
   The physical parameters based on mixture rules for each.
   ********************************************************/
   //mixture molecular mass
   double Mass(void) const;   
   double Rtot(void);        
   double Rtot(void) const;   
   double Cp(void) const;    
   double Cp(const double& TEMP) const;
   double Cv(void) const;     
   double g(void) const;    
   //mixture absolute (sensible+chemical) internal energy
   double e(void) const;     
   double eref(void) const;
   //mixture sensible internal energy  
   double es(void) const; 
   //mixture specific enthalpy   
   double h(void) const;      
   double h(const double &T) const;
   double href(void) const;
   double hs(void) const;
   double hs(const double &T) const;
   double E(void) const;      
   double H(void) const;     
   double Hs(void) const;    
   double mu(void) const;    
   double kappa(void) const; 
   double hprime(void) const;  
   double hprime(double &Temp) const;
   Vector3D rhov(void) const;      
   double T(void) const;         
   //Determine temperature knowing sensible enthalpy 
   double T(double &h_s) const;   
   double gamma_guess(void) const; 
   double a(void);                
   double a(void) const;

   bool negative_speccheck(void) const;
   //Species i concentration (rho*c/mol_mass)
   double SpecCon(int i) const;     
   //Gibbs Free Energy (H-TS) for species
   double Gibbs(int species) const;

   double Schmidt_No(const int &) const;
   double Prandtl() const;
   double Lewis(const int &) const;

   //Temperature Derivatives 
   double diedip() const;
   double diedirho() const;
   double dmudT(void) const;
   double dkappadT(void) const;

   // Conserved solution state. /
   Euler3D_ThermallyPerfect_cState U(void); 
   Euler3D_ThermallyPerfect_cState U(void)const; 
   Euler3D_ThermallyPerfect_cState U(const Euler3D_ThermallyPerfect_pState &W); 
 
   //Fluxes /
   Euler3D_ThermallyPerfect_cState F(void);
   Euler3D_ThermallyPerfect_cState F(void) const;
   Euler3D_ThermallyPerfect_cState F(const Euler3D_ThermallyPerfect_pState &W);

   /* Eigenvalue(s), Eigenvectors (x-direction). */ 
   Euler3D_ThermallyPerfect_pState lambda_x(void);
   Euler3D_ThermallyPerfect_cState rc_x(const int &index); 
   Euler3D_ThermallyPerfect_pState lp_x(const int &index) ; 
   Euler3D_ThermallyPerfect_pState lambda_x(void) const;
   Euler3D_ThermallyPerfect_cState rc_x(const int &index) const; 
   Euler3D_ThermallyPerfect_pState lp_x(const int &index) const; 
  
   /* Index operator */
   double &operator[](int index);
   const double &operator[](int index) const;

   /* Binary arithmetic operators. */
   Euler3D_ThermallyPerfect_pState operator +(const Euler3D_ThermallyPerfect_pState &W) const;
   Euler3D_ThermallyPerfect_pState operator -(const Euler3D_ThermallyPerfect_pState &W) const;
   Euler3D_ThermallyPerfect_pState operator *(const double &a) const;
   friend Euler3D_ThermallyPerfect_pState operator *(const double &a, const Euler3D_ThermallyPerfect_pState &W);
   Euler3D_ThermallyPerfect_pState operator /(const double &a) const;
   double operator *(const Euler3D_ThermallyPerfect_pState &W) const;
   Euler3D_ThermallyPerfect_pState operator ^(const Euler3D_ThermallyPerfect_pState &W) const;
  
   /* Assignment Operator. */ 
   Euler3D_ThermallyPerfect_pState& operator =(const Euler3D_ThermallyPerfect_pState &W); 
  
   /* Shortcut arithmetic operators. */
   Euler3D_ThermallyPerfect_pState& operator +=(const Euler3D_ThermallyPerfect_pState &W);
   Euler3D_ThermallyPerfect_pState& operator -=(const Euler3D_ThermallyPerfect_pState &W);
  
   // Flux functions
   static Euler3D_ThermallyPerfect_pState RoeAverage(const Euler3D_ThermallyPerfect_pState &Wl,
                                                     const Euler3D_ThermallyPerfect_pState &Wr);

   static Euler3D_ThermallyPerfect_cState FluxHLLE_x(const Euler3D_ThermallyPerfect_pState &Wl,
                                                     const Euler3D_ThermallyPerfect_pState &Wr);
   static Euler3D_ThermallyPerfect_cState FluxHLLE_x(const Euler3D_ThermallyPerfect_cState &Ul,
                                                     const Euler3D_ThermallyPerfect_cState &Ur);
   static  Euler3D_ThermallyPerfect_cState FluxHLLE_n(const Euler3D_ThermallyPerfect_pState &Wl,
                                                      const Euler3D_ThermallyPerfect_pState &Wr,
                                                      const Vector3D &norm_dir);
   static Euler3D_ThermallyPerfect_cState FluxHLLE_n(const Euler3D_ThermallyPerfect_cState &Ul,
                                                     const Euler3D_ThermallyPerfect_cState &Ur,
                                                     const Vector3D &norm_dir);

   static Euler3D_ThermallyPerfect_cState FluxRoe_x(const Euler3D_ThermallyPerfect_pState &Wl, 
                                                    const Euler3D_ThermallyPerfect_pState &Wr);
   static Euler3D_ThermallyPerfect_cState FluxRoe_n(const Euler3D_ThermallyPerfect_pState &Wl,
                                                    const Euler3D_ThermallyPerfect_pState &Wr,
                                                    const Vector3D &norm_dir);

   friend Euler3D_ThermallyPerfect_pState HartenFixNeg(const Euler3D_ThermallyPerfect_pState  &lambda_a,
                                                       const Euler3D_ThermallyPerfect_pState  &lambda_l,
                                                       const Euler3D_ThermallyPerfect_pState  &lambda_r);
   friend Euler3D_ThermallyPerfect_pState HartenFixPos(const Euler3D_ThermallyPerfect_pState  &lambda_a,
                                                       const Euler3D_ThermallyPerfect_pState  &lambda_l,
                                                       const Euler3D_ThermallyPerfect_pState  &lambda_r);
				  
   // Boundary Conditions
   static Euler3D_ThermallyPerfect_pState Reflect(const Euler3D_ThermallyPerfect_pState &W, 
                                                  const Vector3D &norm_dir);

  static Euler3D_ThermallyPerfect_pState Moving_Wall(const Euler3D_ThermallyPerfect_pState &Win,
                                                     const Euler3D_ThermallyPerfect_pState &Wout,
                                                     const Vector3D &norm_dir,				 
                                                     const Vector3D &wall_velocity,
                                                     const Vector3D &pressure_gradient,
                                                     const int &TEMPERATURE_BC_FLAG);

   static Euler3D_ThermallyPerfect_pState No_Slip(const Euler3D_ThermallyPerfect_pState &Win, 
                                                  const Euler3D_ThermallyPerfect_pState &Wout,
                                                  const Vector3D &norm_dir,  
                                                  const Vector3D &pressure_gradient,
                                                  const int &TEMPERATURE_BC_FLAG);

   /****** Source terms associated with finite-rate chemistry ******/
   Euler3D_ThermallyPerfect_cState Sw(int &REACT_SET_FLAG) const;
   void dSwdU(DenseMatrix &dSwdU) const; //Jacobian
   double dSwdU_max_diagonal(void) const;
   //@}
  
   /* Unary arithmetic operators. */
   //friend Euler3D_ThermallyPerfect_pState operator -(
   //   const Euler3D_ThermallyPerfect_pState &W);
  
   /* Relational operators. */
   friend int operator ==(const Euler3D_ThermallyPerfect_pState &W1, 
                          const Euler3D_ThermallyPerfect_pState &W2);
   friend int operator !=(const Euler3D_ThermallyPerfect_pState &W1, 
                          const Euler3D_ThermallyPerfect_pState &W2);

   //@{ @name Input-output operators:
   friend ostream& operator << (ostream &out_file, 
                                const Euler3D_ThermallyPerfect_pState &W);
   friend istream& operator >> (istream &in_file,  
                                Euler3D_ThermallyPerfect_pState &W);
   //@}
 
};

/*!
 * Class: Euler3D_ThermallyPerfect_cState
 *
 * @brief Euler 3D conservative solution state class.
 *
 * This class defines the primitive solution state for the
 * Euler equations governing three-dimensional flows of
 * thermally perfect non-reactive and combusting mixtures.
 *
 */
class Euler3D_ThermallyPerfect_cState {
  public:
   double         rho;  
   Vector3D      rhov;   
   double           E;  
   Species   *rhospec;   

   static int num_vars;  
   static int      ns;           
   static NASARP1311data *specdata;  
   static double low_temp_range;    
   static double high_temp_range;   
   static int debug_level;          
   static double *Schmidt;          
 
   //constructors
   Euler3D_ThermallyPerfect_cState(): 
    rho(DENSITY_STDATM), E(PRESSURE_STDATM/(rho*(0.4))), rhospec(NULL) {
      rhov.zero(); set_initial_values(); 
   }   

   Euler3D_ThermallyPerfect_cState(const double &value): 
    rho(value), E(value), rhospec(NULL) {
      rhov.x = value; rhov.y = value; rhov.z = value;
      set_initial_values(value); 
   }   
  
   Euler3D_ThermallyPerfect_cState(const double &d, 
                                   const double &vx, 
                                   const double &vy, 
                                   const double &vz,
                                   const double &En ):
    rho(d), E(En), rhospec(NULL) {
      rhov.x=vx; rhov.y=vy; rhov.z=vz; 
      set_initial_values(); 
   }
  
   Euler3D_ThermallyPerfect_cState(const double &d, 
                                   const double &vx, 
                                   const double &vy, 
                                   const double &vz, 
                                   const double &En,	
                                   Species *rhomfrac): 
    rho(d), E(En), rhospec(NULL) { 
      rhov.x=vx; rhov.y=vy; rhov.z = vz;
      set_initial_values(rhomfrac); 
   }

   Euler3D_ThermallyPerfect_cState(const double &d, 
                                   const double &vx, 
                                   const double &vy, 
                                   const double &vz, 
                                   const double &En,	
                                   const double &rhomfrac): 
    rho(d), E(En), rhospec(NULL) { 
      rhov.x=vx; rhov.y=vy; rhov.z = vz;
      set_initial_values(rhomfrac); 
   }
  
   Euler3D_ThermallyPerfect_cState(const double &d, 
                                   const Vector3D &V, 
                                   const double &En):
    rho(d), rhov(V), E(En), rhospec (NULL) {
      set_initial_values();
   } 

   Euler3D_ThermallyPerfect_cState(const double &d, 
                                   const Vector3D &V, 
                                   const double &En, 
                                   const double &frac):
    rho(d), rhov(V), E(En), rhospec (NULL) {
      set_initial_values(frac);
   } 

   Euler3D_ThermallyPerfect_cState(const double &d, 
                                   const Vector3D &dV, 
                                   const double &En, 
                                   Species *rhomfrac):
    rho(d), rhov(dV), E(En), rhospec(NULL){
      set_initial_values(rhomfrac);
   }

   // This is needed for the operator overload returns!!!!
   Euler3D_ThermallyPerfect_cState(const Euler3D_ThermallyPerfect_cState &U) { 
      rhospec = NULL; rho = DENSITY_STDATM; 
      set_initial_values();
      Copy(U); 
   }
  
   // Sets the species data, needs to be called only once as it is static
   void set_species_data(const int &n,
                         const string *S,
                         const char *PATH,
                         const int &debug, 
                         const double &Mr, 
                         const double* Sc,
                         const int &trans_data);

   //Set initial data values  
   void set_initial_values();
   void set_initial_values(const double &value);
   void set_initial_values(double *rhomfrac);
   void set_initial_values(Species *rhomfrac);
     
   // Copy operator
   void Copy(const Euler3D_ThermallyPerfect_cState &U);

   /*************** VACUUM OPERATOR *********************/
   void Vacuum(){
      rho=ZERO; rhov.x=ZERO; rhov.y=ZERO; rhov.z=ZERO; E=ZERO;
      for(int i=0; i<ns; i++){
         rhospec[i].Vacuum();
      }
   }  

   void zero_non_sol(){
     for(int i=0; i<ns; i++){
       rhospec[i].gradc.zero();
       rhospec[i].diffusion_coef=ZERO;
     }
   }

   /* Set Data Temperature Ranges */
   void Temp_low_range();     
   void Temp_high_range(); 

   // Mixing Rules 
   double Rtot(void) const; 
   double Cp(void) const;  
   double Cv(void) const;  
   double g(void) const;   
   double gamma_guess(void) const;  
   double e(void) const;            
   double es(void) const;           
   double h(const double &T) const;  
   double hs(const double &T) const; 
   double hprime(const double &T) const; 
   double heatofform(void) const;      
   double mu(void) const;             
   double kappa(void) const;         
   Vector3D v(void) const;   
   double p(void) const;    
   double T(void) const;      
   double a(void) const;     
   bool negative_speccheck(  const Euler3D_ThermallyPerfect_cState &Uo, const int &flag) const; 
   double sum_species(void) const;
 
   // Temperature Derivatives 
   double dmudT(void) const;

   Euler3D_ThermallyPerfect_pState W(void) ;
   Euler3D_ThermallyPerfect_pState W(void) const;
   Euler3D_ThermallyPerfect_pState W(const Euler3D_ThermallyPerfect_cState &U) const;
   friend Euler3D_ThermallyPerfect_pState W(const Euler3D_ThermallyPerfect_cState &U);

   /* Index operators */
   double &operator[](int index);
   const double &operator[](int index) const;
 
   /* Binary arithmetic operators. */
   Euler3D_ThermallyPerfect_cState operator +(const Euler3D_ThermallyPerfect_cState &U) const;
   Euler3D_ThermallyPerfect_cState operator -(const Euler3D_ThermallyPerfect_cState &U) const;
   Euler3D_ThermallyPerfect_cState operator *(const double &a) const;
   friend Euler3D_ThermallyPerfect_cState operator *(const double &a, const Euler3D_ThermallyPerfect_cState &U);
   Euler3D_ThermallyPerfect_cState operator /(const double &a) const;
   double operator *(const Euler3D_ThermallyPerfect_cState &U) const;
   Euler3D_ThermallyPerfect_cState operator ^(const Euler3D_ThermallyPerfect_cState &U) const;

   /* Assignment Operator. */ 
   Euler3D_ThermallyPerfect_cState& operator =(const Euler3D_ThermallyPerfect_cState &U); 

   /* Shortcut arithmetic operators. */
   Euler3D_ThermallyPerfect_cState& operator +=(const Euler3D_ThermallyPerfect_cState &U);
   Euler3D_ThermallyPerfect_cState& operator -=(const Euler3D_ThermallyPerfect_cState &U);
      
   /* Unary arithmetic operators. */
   friend Euler3D_ThermallyPerfect_cState operator -(const Euler3D_ThermallyPerfect_cState &U);
  
   /* Relational operators. */
   friend int operator ==(const Euler3D_ThermallyPerfect_cState &U1, 
                          const Euler3D_ThermallyPerfect_cState &U2);
   friend int operator !=(const Euler3D_ThermallyPerfect_cState &U1, 
                          const Euler3D_ThermallyPerfect_cState &U2);
  
   /* Input-output operators. */
   friend ostream& operator << (ostream &out_file, 
                                const Euler3D_ThermallyPerfect_cState &U);
   friend istream& operator >> (istream &in_file,  
                                Euler3D_ThermallyPerfect_cState &U);

   /* Destructors */
   void Deallocate_static(void) {
     if (specdata != NULL) delete[] specdata;
     specdata = NULL;
     if(Schmidt != NULL)  delete[] Schmidt;
     Schmidt = NULL;
   }

   void Deallocate(void) {
     if (rhospec != NULL) delete[] rhospec;
     rhospec = NULL;
   }
  
   ~Euler3D_ThermallyPerfect_cState(){ 
     Deallocate();
   }
  
};

/**************************************************************************
*Euler3D_ThermallyPerfect_pState member functions  ************************
***************************************************************************/

//   Set up mass fractions memory and initial values, ie. Species Class  
inline void Euler3D_ThermallyPerfect_pState::set_initial_values(){
   Deallocate();
   spec = new Species[ns];
   for(int i=0; i<ns; i++){
      spec[i].c = ONE/ns ; 
   }
}

inline void  Euler3D_ThermallyPerfect_pState::set_initial_values(const double &value){
   Deallocate();
   spec = new Species[ns];
   for(int i=0; i<ns; i++){
      spec[i].c = value ; 
   }
}

//user specified
inline void Euler3D_ThermallyPerfect_pState::set_initial_values(double *cfrac){
   Deallocate();
  spec = new Species[ns];
  for(int i=0; i<ns; i++){
    spec[i].c = cfrac[i];
  }
}

//another set using species class
inline void Euler3D_ThermallyPerfect_pState::set_initial_values(Species *mfrac){
   Deallocate();
   spec = new Species[ns];
   for(int i=0; i<ns; i++){
      spec[i] = mfrac[i];
   }
}

//momentum
inline Vector3D Euler3D_ThermallyPerfect_pState::rhov(void) const{
   return rho*v; 
}
 
inline double Euler3D_ThermallyPerfect_pState::mu(void) const{
 double sum =0.0;
  double Temp = T();
  double  *vis = new double[ns];

  for(int i=0; i<ns; i++){
    double phi = 0.0;
    for (int j=0; j<ns; j++){
      if(i == 0){
        vis[j] = specdata[j].Viscosity(Temp);
      }
      phi += (spec[j].c / specdata[j].Mol_mass())*
        pow(ONE + sqrt(vis[i]/vis[j])*
            pow(specdata[j].Mol_mass()/specdata[i].Mol_mass(),0.25),2.0)/
        sqrt(EIGHT*(ONE +specdata[i].Mol_mass()/specdata[j].Mol_mass()));
    }
    sum += (spec[i].c * vis[i]) /
      (specdata[i].Mol_mass() * phi);
  }

  delete[] vis;

  return sum;
}

//Prandtl number
inline double Euler3D_ThermallyPerfect_pState::Prandtl(void) const{
   //Pr = Cp*mu/k
   return Cp()*mu()/kappa();
}

//Schmidt no.
inline double Euler3D_ThermallyPerfect_pState::Schmidt_No(const int &i) const{
   if(spec[i].diffusion_coef != ZERO){
      return mu()/(rho*spec[i].diffusion_coef);
   } else {
      return Schmidt[i];
   }
   
}

//Lewis number
inline double Euler3D_ThermallyPerfect_pState::Lewis(const int &i) const{
   if(spec[i].diffusion_coef != ZERO){
      return kappa()/(rho*Cp()*spec[i].diffusion_coef);
   }
   return ZERO;
}

//Temperature
inline double Euler3D_ThermallyPerfect_pState::T(void) const{

  return p/(rho*Rtot());

}
//Species concentrations
inline double Euler3D_ThermallyPerfect_pState::SpecCon(int i) const{
  //returned in kg/m^3 / kg/mol => mol/m^3
  return (rho)*spec[i].c/(specdata[i].Mol_mass());
}

/******* GIBBS Free Energy ********************************
  Eqn. (10.84) (10.87) Anderson
  Gs = Hs - TS
  Gs(ps=1) = Gs - R_UNIVERSAL*T*ln(ps)  //for data not at 1atm
  ps = cs(M/Ms)*p
***********************************************************/
inline double Euler3D_ThermallyPerfect_pState::Gibbs(int species) const{
   double Temp = T(); 
   return specdata[species].GibbsFree(Temp);
}

// Temperature Derivatives 
inline double Euler3D_ThermallyPerfect_pState::diedip() const{
   double Rlocal = Rtot();
   return (hprime() - Rlocal)/(rho*Rlocal);
}

inline double Euler3D_ThermallyPerfect_pState::diedirho() const{
   double Rlocal = Rtot();
   return -p*(hprime() - Rlocal)/(rho*rho*Rlocal);
}

//Copy
inline void Euler3D_ThermallyPerfect_pState::Copy(
   const Euler3D_ThermallyPerfect_pState &W){
  rho = W.rho;
  v = W.v; 
  p = W.p;  

  for( int i=0; i<ns; i++){
    spec[i] = W.spec[i];
  }
}

//----------------- Index Operator ------------------------/
inline double& Euler3D_ThermallyPerfect_pState::operator[](int index) {
   switch(index){  
   case 1:
      return rho;    
   case 2:
      return v.x;
   case 3:
      return v.y;
   case 4:
      return v.z;
   case 5:
      return p;
   default :
      return spec[index-NUM_EULER3D_VAR_SANS_SPECIES-1].c;
   };
}

inline const double& Euler3D_ThermallyPerfect_pState::operator[](int index) const {
   switch(index){  
   case 1:
      return rho;    
   case 2:
      return v.x;
   case 3:
      return v.y;
   case 4:
      return v.z;
   case 5:
      return p;
   default :
      return spec[index-NUM_EULER3D_VAR_SANS_SPECIES-1].c;
   };
}

inline Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::U(void){
   Euler3D_ThermallyPerfect_cState Temp;
   Temp.rho = rho;
   Temp.rhov = rhov();
   Temp.E = E();
   for(int i=0; i<ns; i++){
      Temp.rhospec[i] = rho*spec[i];
      Temp.rhospec[i].gradc = rho*spec[i].gradc;
      Temp.rhospec[i].diffusion_coef = rho*spec[i].diffusion_coef;
   }
 
   return Temp;
}


inline Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::U(void)const{
  
   Euler3D_ThermallyPerfect_cState Temp;
   Temp.rho = rho;
   Temp.rhov = rhov();
   Temp.E = E();
   for(int i=0; i<ns; i++){
      Temp.rhospec[i] = rho*spec[i];
      Temp.rhospec[i].gradc = rho*spec[i].gradc;
      Temp.rhospec[i].diffusion_coef = rho*spec[i].diffusion_coef;
   }
 
   return Temp;
}

inline Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_pState::U(const Euler3D_ThermallyPerfect_pState &W){
  if(ns == W.ns){ //check that species are equal
    Euler3D_ThermallyPerfect_cState Temp;
    Temp.rho = W.rho;
    Temp.rhov = W.rhov();
    Temp.E = W.E();
    for(int i=0; i<W.ns; i++){
      Temp.rhospec[i] = W.rho*W.spec[i];
      Temp.rhospec[i].gradc = W.rho*W.spec[i].gradc;
      Temp.rhospec[i].diffusion_coef = W.rho*W.spec[i].diffusion_coef;
    }
 
    return Temp;
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  }
}

inline Euler3D_ThermallyPerfect_cState U(const Euler3D_ThermallyPerfect_pState &W) {
   Euler3D_ThermallyPerfect_cState Temp;
   Temp.rho = W.rho;
   Temp.rhov = W.rhov();
   Temp.E = W.E();
   for(int i=0; i<W.ns; i++){
      Temp.rhospec[i] = W.rho*W.spec[i];
      Temp.rhospec[i].gradc = W.rho*W.spec[i].gradc;
      Temp.rhospec[i].diffusion_coef = W.rho*W.spec[i].diffusion_coef;
   }
   return Temp;
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
inline bool Euler3D_ThermallyPerfect_pState::negative_speccheck(void) const{
   double sum = ZERO;
   
   double LOCAL_TOL = MICRO; 
   //-------- Negative Check ------------//
   for(int i=0; i<ns-1; i++){
      if(spec[i].c > ONE){ //check for > 1.0 
         spec[i].c = ONE;
      } else if(spec[i].c < ZERO){ //check for -ve
         if(spec[i].c > -LOCAL_TOL){  //check for small -ve and set to ZERO 
            spec[i].c = ZERO;
         } else {
            spec[i].c = ZERO;
         }
      }
      sum += spec[i].c;
   }
   spec[ns-1].c = max(ONE- sum, ZERO);
   sum += spec[ns-1].c;
   
   for(int i=0; i<ns; i++){
      spec[i].c = spec[i].c*(ONE/sum);
   }
   
   return true;
   
}

/**************************************************************************
*Euler3D_ThermallyPerfect_cState member functions  ************************
***************************************************************************/
//   Set up mass fractions memory and initial values, ie. Species Class  

inline void Euler3D_ThermallyPerfect_cState::set_initial_values(){

   Deallocate();
  rhospec = new Species[ns];
  for(int i=0; i<ns; i++){
    rhospec[i].c = rho/ns; 
  }
}

inline void  Euler3D_ThermallyPerfect_cState::set_initial_values(
   const double &value){

   Deallocate();
   rhospec = new Species[ns];
   for(int i=0; i<ns; i++){
      rhospec[i].c = value; 
   }

}

//user specified
inline void Euler3D_ThermallyPerfect_cState::set_initial_values(double *rhomfrac){
   Deallocate();
   rhospec = new Species[ns];
   for(int i=0; i<ns; i++){
      rhospec[i].c = rhomfrac[i];
   }
   
}

//another set using species class
inline void Euler3D_ThermallyPerfect_cState::set_initial_values(Species *rhomfrac){
   
   Deallocate();
   rhospec = new Species[ns];
   for(int i=0; i<ns; i++){
      rhospec[i] = rhomfrac[i];
   }
}

/* Copy */
inline void Euler3D_ThermallyPerfect_cState::Copy(
   const Euler3D_ThermallyPerfect_cState &U){
   rho = U.rho;
   rhov = U.rhov; 
   E = U.E; 
   
   for( int i=0; i<ns; i++){ 
      rhospec[i] = U.rhospec[i];
   } 
 
}
//velocity
inline Vector3D Euler3D_ThermallyPerfect_cState::v() const{
  return (rhov/rho);
}
//pressure
inline double Euler3D_ThermallyPerfect_cState::p() const{
  return (rho*Rtot()*T());
}
//index operators
inline double& Euler3D_ThermallyPerfect_cState::operator[](int index) {

  switch(index){  
   case 1:
      return rho;    
   case 2:
      return rhov.x;
   case 3:
      return rhov.y;
   case 4:
      return rhov.z;
   case 5:
      return E;
   default :
      return rhospec[index-NUM_EULER3D_VAR_SANS_SPECIES-1].c;
   };
    
}

inline const double& Euler3D_ThermallyPerfect_cState::operator[](int index) const{
  
   switch(index){  
   case 1:
      return rho;    
   case 2:
      return rhov.x;
   case 3:
      return rhov.y;
   case 4:
      return rhov.z;
   case 5:
      return E;
   default :
      return rhospec[index-NUM_EULER3D_VAR_SANS_SPECIES-1].c;
   };
 
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
inline bool Euler3D_ThermallyPerfect_cState::negative_speccheck(
   const Euler3D_ThermallyPerfect_cState &Uo, const int &flag) const{
   
   double sum = ZERO;
   double temp = ZERO;
   double LOCAL_TOL = MICRO; 
   //-------- Negative Check ------------//
   for(int i=0; i<ns-1; i++){
      temp = rhospec[i].c/rho;
      if(temp>ONE){ //check for > 1.0
         rhospec[i].c =rho;
         temp = ONE;
         
      }else if (temp < ZERO){  //check for -ve
         if(temp > -LOCAL_TOL){  //check for small -ve and set to ZERO 
            rhospec[i].c = ZERO;
            temp = ZERO;
         } else {
            rhospec[i].c = ZERO;
            temp = ZERO;
         }
      } 
      sum += temp;
   } 
   temp = max(ONE- sum, ZERO);
   sum += temp;
   rhospec[ns-1].c = rho*temp;
   
   for(int i=0; i<ns; i++){
      rhospec[i].c = rhospec[i].c*(ONE/sum);
   }
   
   return true;

}


inline double Euler3D_ThermallyPerfect_cState::mu(void) const{
   double sum = ZERO;

   double Temp = T();
   
  for(int i=0; i<ns; i++){
    double phi = 0.0;
    for (int j=0; j<ns; j++){
      phi += ((rhospec[j].c/rho) / specdata[j].Mol_mass())*
        pow(1.0 + sqrt(specdata[i].Viscosity(Temp)/specdata[j].Viscosity(Temp))*
            pow(specdata[j].Mol_mass()/specdata[i].Mol_mass(),0.25),2.0)/
       sqrt(8.0*(1.0 +specdata[i].Mol_mass()/specdata[j].Mol_mass()));
    } 
    sum += ((rhospec[i].c/rho)* specdata[i].Viscosity(Temp) ) /
      (specdata[i].Mol_mass() * phi);
  }  

  return sum;
}

/*  Sum N-1 species.*/
inline double Euler3D_ThermallyPerfect_cState::sum_species(void) const{
   double sum=ZERO;
   for(int i=0; i<ns-1; i++){
      sum += rhospec[i].c;
   }
   return sum/rho;
}

inline Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_cState::W(void){
  
   Euler3D_ThermallyPerfect_pState Temp;
   Temp.rho = rho;
   Temp.v = v();  
   Temp.p = p();
   
   for(int i=0; i<ns; i++){
      Temp.spec[i] = rhospec[i]/rho;
      Temp.spec[i].gradc = rhospec[i].gradc/rho;
      Temp.spec[i].diffusion_coef = rhospec[i].diffusion_coef/rho;
   }
   
   return Temp;
}

inline Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_cState::W(void)const{
  
   Euler3D_ThermallyPerfect_pState Temp;
   Temp.rho = rho;
   Temp.v = v();  
   Temp.p = p();
   
   for(int i=0; i<ns; i++){
      Temp.spec[i] = rhospec[i]/rho;
      Temp.spec[i].gradc = rhospec[i].gradc/rho;
      Temp.spec[i].diffusion_coef = rhospec[i].diffusion_coef/rho;
   }
   
   return Temp;
}

inline Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_cState::W(
   const Euler3D_ThermallyPerfect_cState &U) const{
  if(ns == U.ns){ //check that species are equal   
    Euler3D_ThermallyPerfect_pState Temp;
    Temp.rho = U.rho;
    Temp.v = U.v();  
    Temp.p = U.p();
  
    for(int i=0; i<U.ns; i++){
      Temp.spec[i] = U.rhospec[i]/U.rho;
      Temp.spec[i].gradc = U.rhospec[i].gradc/U.rho;
      Temp.spec[i].diffusion_coef = U.rhospec[i].diffusion_coef/U.rho;
    }
  
    return Temp;
  } else {
    cerr<<"\n Mismatch in number of species \n";
    exit(1);
  } 
}

inline Euler3D_ThermallyPerfect_pState W(const Euler3D_ThermallyPerfect_cState &U) {
  Euler3D_ThermallyPerfect_pState Temp;
  Temp.rho = U.rho;
  Temp.v = U.v();
  Temp.p = U.p();

  for(int i=0; i<U.ns; i++){
    Temp.spec[i] = U.rhospec[i]/U.rho;
    Temp.spec[i].gradc = U.rhospec[i].gradc/U.rho;
    Temp.spec[i].diffusion_coef = U.rhospec[i].diffusion_coef/U.rho;
  }
 
  return Temp;
}

#endif // _EULER3D_THERMALLYPERFECT_STATE_INCLUDED
